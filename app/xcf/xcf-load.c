/* GIMP - The GNU Image Manipulation Program
 * Copyright (C) 1995 Spencer Kimball and Peter Mattis
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "config.h"

#include <string.h>
#include <zlib.h>

#include <cairo.h>
#include <gegl.h>
#include <gegl-plugin.h>
#include <gdk-pixbuf/gdk-pixbuf.h>

#include "libgimpbase/gimpbase.h"
#include "libgimpcolor/gimpcolor.h"
#include "libgimpconfig/gimpconfig.h"

#include "core/core-types.h"

#include "config/gimpcoreconfig.h"

#include "gegl/gimp-babl.h"
#include "gegl/gimp-gegl-tile-compat.h"

#include "core/gimp.h"
#include "core/gimpcontainer.h"
#include "core/gimpdrawable-filters.h"
#include "core/gimpdrawable-private.h" /* eek */
#include "core/gimpdrawablefilter.h"
#include "core/gimpfilterstack.h"
#include "core/gimpgrid.h"
#include "core/gimpgrouplayer.h"
#include "core/gimpimage.h"
#include "core/gimpimage-color-profile.h"
#include "core/gimpimage-colormap.h"
#include "core/gimpimage-grid.h"
#include "core/gimpimage-guides.h"
#include "core/gimpimage-metadata.h"
#include "core/gimpimage-private.h"
#include "core/gimpimage-sample-points.h"
#include "core/gimpimage-symmetry.h"
#include "core/gimpimage-undo.h"
#include "core/gimpitemlist.h"
#include "core/gimpitemstack.h"
#include "core/gimplayer-floating-selection.h"
#include "core/gimplayer-new.h"
#include "core/gimplayermask.h"
#include "core/gimpparasitelist.h"
#include "core/gimpprogress.h"
#include "core/gimpselection.h"
#include "core/gimpsymmetry.h"
#include "core/gimptemplate.h"
#include "core/gimpunit.h"

#include "operations/layer-modes/gimp-layer-modes.h"

#include "path/gimpanchor.h"
#include "path/gimpstroke.h"
#include "path/gimpbezierstroke.h"
#include "path/gimppath.h"
#include "path/gimppath-compat.h"

#include "text/gimptextlayer.h"
#include "text/gimptextlayer-xcf.h"

#include "xcf-private.h"
#include "xcf-load.h"
#include "xcf-read.h"
#include "xcf-seek.h"
#include "xcf-utils.h"

#include "gimp-log.h"
#include "gimp-intl.h"


#define MAX_XCF_PARASITE_DATA_LEN (256L * 1024 * 1024)

/* #define GIMP_XCF_PATH_DEBUG */

/* Filters can not be created until a layer is attached
 * to an image, so we use this struct to store relevant
 * information until then */
typedef struct
{
  GeglNode              *operation;
  gchar                 *name;
  gchar                 *icon_name;
  gchar                 *operation_name;
  gchar                 *op_version;
  GimpChannel           *mask;
  gboolean               is_visible;
  gdouble                opacity;
  GimpLayerMode          paint_mode;
  GimpLayerColorSpace    blend_space;
  GimpLayerColorSpace    composite_space;
  GimpLayerCompositeMode composite_mode;
  gboolean               clip;
  GimpFilterRegion       region;

  gboolean               unsupported_operation;
} FilterData;

static void            xcf_load_add_masks     (GimpImage     *image);
static void            xcf_load_add_effects   (XcfInfo       *info,
                                               GimpImage     *image);
static gboolean        xcf_load_image_props   (XcfInfo       *info,
                                               GimpImage     *image);
static gboolean        xcf_load_layer_props   (XcfInfo       *info,
                                               GimpImage     *image,
                                               GimpLayer    **layer,
                                               GList        **item_path,
                                               gboolean      *apply_mask,
                                               gboolean      *edit_mask,
                                               gboolean      *show_mask,
                                               guint32       *text_layer_flags,
                                               guint32       *group_layer_flags);
static gboolean        xcf_check_layer_props  (XcfInfo       *info,
                                               GList        **item_path,
                                               gboolean      *is_group_layer,
                                               gboolean      *is_text_layer);
static gboolean        xcf_load_channel_props (XcfInfo       *info,
                                               GimpImage     *image,
                                               GimpChannel  **channel,
                                               gboolean       is_mask);
static gboolean        xcf_load_effect_props  (XcfInfo       *info,
                                               FilterData    *filter);
static gboolean        xcf_load_path_props    (XcfInfo       *info,
                                               GimpImage     *image,
                                               GimpPath     **paths);
static gboolean        xcf_load_prop          (XcfInfo       *info,
                                               PropType      *prop_type,
                                               guint32       *prop_size);
static GimpLayer     * xcf_load_layer         (XcfInfo       *info,
                                               GimpImage     *image,
                                               GList        **item_path,
                                               gint          *n_broken_effects);
static GimpChannel   * xcf_load_channel       (XcfInfo       *info,
                                               GimpImage     *image);
static FilterData    * xcf_load_effect        (XcfInfo       *info,
                                               GimpImage     *image,
                                               GimpDrawable  *drawable);
static void            xcf_load_free_effect   (FilterData    *data);
static void            xcf_load_free_effects  (GList         *effects);
static GimpPath      * xcf_load_path          (XcfInfo       *info,
                                               GimpImage     *image);
static GimpLayerMask * xcf_load_layer_mask    (XcfInfo       *info,
                                               GimpImage     *image);
static gboolean        xcf_load_buffer        (XcfInfo       *info,
                                               GeglBuffer    *buffer);
static gboolean        xcf_load_level         (XcfInfo       *info,
                                               GeglBuffer    *buffer);
static gboolean        xcf_load_tile          (XcfInfo       *info,
                                               GeglBuffer    *buffer,
                                               GeglRectangle *tile_rect,
                                               const Babl    *format);
static gboolean        xcf_load_tile_rle      (XcfInfo       *info,
                                               GeglBuffer    *buffer,
                                               GeglRectangle *tile_rect,
                                               const Babl    *format,
                                               gint           data_length);
static gboolean        xcf_load_tile_zlib     (XcfInfo       *info,
                                               GeglBuffer    *buffer,
                                               GeglRectangle *tile_rect,
                                               const Babl    *format,
                                               gint           data_length);
static GimpParasite  * xcf_load_parasite      (XcfInfo       *info);
static gboolean        xcf_load_old_paths     (XcfInfo       *info,
                                               GimpImage     *image);
static gboolean        xcf_load_old_path      (XcfInfo       *info,
                                               GimpImage     *image);
static gboolean        xcf_load_old_vectors   (XcfInfo       *info,
                                               GimpImage     *image);
static gboolean        xcf_load_old_vector    (XcfInfo       *info,
                                               GimpImage     *image);

static gboolean        xcf_skip_unknown_prop  (XcfInfo       *info,
                                               gsize          size);

static gboolean        xcf_item_path_is_parent (GList        *path,
                                                GList        *parent_path);
static void            xcf_fix_item_path       (GimpLayer    *layer,
                                                GList       **path,
                                                GList        *broken_paths);

#define xcf_progress_update(info) G_STMT_START  \
  {                                             \
    if (info->progress)                         \
      gimp_progress_pulse (info->progress);     \
  } G_STMT_END


GimpImage *
xcf_load_image (Gimp     *gimp,
                XcfInfo  *info,
                GError  **error)
{
  GimpImage          *image = NULL;
  const GimpParasite *parasite;
  gboolean            has_metadata = FALSE;
  goffset             saved_pos;
  goffset             offset;
  gint                width;
  gint                height;
  gint                image_type;
  GimpPrecision       precision = GIMP_PRECISION_U8_NON_LINEAR;
  gint                num_successful_elements = 0;
  gint                n_broken_layers         = 0;
  gint                n_broken_channels       = 0;
  gint                n_broken_paths          = 0;
  gint                n_broken_effects        = 0;
  GList              *broken_paths            = NULL;
  GList              *group_layers            = NULL;
  GList              *syms;
  GList              *iter;

  /* read in the image width, height and type */
  xcf_read_int32 (info, (guint32 *) &width, 1);
  xcf_read_int32 (info, (guint32 *) &height, 1);
  xcf_read_int32 (info, (guint32 *) &image_type, 1);
  if (image_type < GIMP_RGB || image_type > GIMP_INDEXED)
    goto hard_error;

  /* Be lenient with corrupt image dimensions.
   * Hopefully layer dimensions will be valid. */
  if (width <= 0 || height <= 0 ||
      width > GIMP_MAX_IMAGE_SIZE || height > GIMP_MAX_IMAGE_SIZE)
    {
      GIMP_LOG (XCF, "Invalid image size %d x %d, setting to 1x1.", width, height);
      width  = 1;
      height = 1;
    }

  if (info->file_version >= 4)
    {
      gint p;

      xcf_read_int32 (info, (guint32 *) &p, 1);

      if (info->file_version == 4)
        {
          switch (p)
            {
            case 0: precision = GIMP_PRECISION_U8_NON_LINEAR;  break;
            case 1: precision = GIMP_PRECISION_U16_NON_LINEAR; break;
            case 2: precision = GIMP_PRECISION_U32_LINEAR;     break;
            case 3: precision = GIMP_PRECISION_HALF_LINEAR;    break;
            case 4: precision = GIMP_PRECISION_FLOAT_LINEAR;   break;
            default:
              goto hard_error;
            }
        }
      else if (info->file_version == 5 ||
               info->file_version == 6)
        {
          switch (p)
            {
            case 100: precision = GIMP_PRECISION_U8_LINEAR;        break;
            case 150: precision = GIMP_PRECISION_U8_NON_LINEAR;    break;
            case 200: precision = GIMP_PRECISION_U16_LINEAR;       break;
            case 250: precision = GIMP_PRECISION_U16_NON_LINEAR;   break;
            case 300: precision = GIMP_PRECISION_U32_LINEAR;       break;
            case 350: precision = GIMP_PRECISION_U32_NON_LINEAR;   break;
            case 400: precision = GIMP_PRECISION_HALF_LINEAR;      break;
            case 450: precision = GIMP_PRECISION_HALF_NON_LINEAR;  break;
            case 500: precision = GIMP_PRECISION_FLOAT_LINEAR;     break;
            case 550: precision = GIMP_PRECISION_FLOAT_NON_LINEAR; break;
            default:
              goto hard_error;
            }
        }
      else
        {
          precision = p;
        }
    }

  GIMP_LOG (XCF, "version=%d, width=%d, height=%d, image_type=%d, precision=%d",
            info->file_version, width, height, image_type, precision);

  if (! gimp_babl_is_valid (image_type, precision))
    {
      gimp_message_literal (gimp, G_OBJECT (info->progress),
                            GIMP_MESSAGE_ERROR,
                            _("Invalid image mode and precision combination."));
      goto hard_error;
    }

  image = gimp_create_image (gimp, width, height, image_type, precision,
                             FALSE);

  gimp_image_undo_disable (image);

  xcf_progress_update (info);

  /* read the image properties */
  if (! xcf_load_image_props (info, image))
    goto hard_error;

  GIMP_LOG (XCF, "image props loaded");

  /* Order matters for item sets. */
  info->layer_sets = g_list_reverse (info->layer_sets);
  info->channel_sets = g_list_reverse (info->channel_sets);

  /* check for simulation intent parasite */
  parasite = gimp_image_parasite_find (GIMP_IMAGE (image),
                                       "image-simulation-intent");
  if (parasite)
    {
      guint32           parasite_size;
      const guint8     *intent;
      GimpImagePrivate *private = GIMP_IMAGE_GET_PRIVATE (image);

      intent = (const guint8 *) gimp_parasite_get_data (parasite, &parasite_size);

      if (parasite_size == 1)
        {
          if (*intent != GIMP_COLOR_RENDERING_INTENT_PERCEPTUAL            &&
              *intent != GIMP_COLOR_RENDERING_INTENT_RELATIVE_COLORIMETRIC &&
              *intent != GIMP_COLOR_RENDERING_INTENT_SATURATION            &&
              *intent != GIMP_COLOR_RENDERING_INTENT_ABSOLUTE_COLORIMETRIC)
            {
              gimp_message (info->gimp, G_OBJECT (info->progress),
                            GIMP_MESSAGE_ERROR,
                            "Unknown simulation rendering intent: %d",
                            *intent);
            }
          else
            {
              gimp_image_set_simulation_intent (image,
                                                (GimpColorRenderingIntent) *intent);
            }
        }
      else
        {
          gimp_message (info->gimp, G_OBJECT (info->progress),
                        GIMP_MESSAGE_ERROR,
                        "Invalid simulation intent data");
        }

      gimp_parasite_list_remove (private->parasites,
                                 gimp_parasite_get_name (parasite));
    }


/* check for simulation bpc parasite */
  parasite = gimp_image_parasite_find (GIMP_IMAGE (image),
                                       "image-simulation-bpc");
  if (parasite)
    {
      guint32           parasite_size;
      const guint8     *bpc;
      gboolean          status  = FALSE;
      GimpImagePrivate *private = GIMP_IMAGE_GET_PRIVATE (image);

      bpc = (const guint8 *) gimp_parasite_get_data (parasite, &parasite_size);

      if (parasite_size == 1)
        {
          if (*bpc)
            status = TRUE;

          gimp_image_set_simulation_bpc (image, status);
        }
      else
        {
          gimp_message (info->gimp, G_OBJECT (info->progress),
                        GIMP_MESSAGE_ERROR,
                        "Invalid simulation bpc data");
        }

      gimp_parasite_list_remove (private->parasites,
                                 gimp_parasite_get_name (parasite));
    }

  /* check for a GimpGrid parasite */
  parasite = gimp_image_parasite_find (GIMP_IMAGE (image),
                                       gimp_grid_parasite_name ());
  if (parasite)
    {
      GimpGrid *grid = gimp_grid_from_parasite (parasite);

      if (grid)
        {
          GimpImagePrivate *private = GIMP_IMAGE_GET_PRIVATE (image);

          gimp_parasite_list_remove (private->parasites,
                                     gimp_parasite_get_name (parasite));

          gimp_image_set_grid (GIMP_IMAGE (image), grid, FALSE);
          g_object_unref (grid);
        }
    }

  /* check for a metadata parasite */
  parasite = gimp_image_parasite_find (GIMP_IMAGE (image),
                                       "gimp-image-metadata");
  if (parasite)
    {
      GimpImagePrivate *private  = GIMP_IMAGE_GET_PRIVATE (image);
      GimpMetadata     *metadata = NULL;
      gchar            *meta_string;
      guint32           parasite_data_size;

      meta_string = (gchar *) gimp_parasite_get_data (parasite, &parasite_data_size);
      if (meta_string)
        {
          meta_string = g_strndup (meta_string, parasite_data_size);
          metadata = gimp_metadata_deserialize (meta_string);
          g_free (meta_string);
        }

      if (metadata)
        {
          has_metadata = TRUE;

          gimp_image_set_metadata (image, metadata, FALSE);
          g_object_unref (metadata);
        }

      gimp_parasite_list_remove (private->parasites,
                                 gimp_parasite_get_name (parasite));
    }

  /* check for symmetry parasites */
  syms = gimp_image_symmetry_list ();
  for (iter = syms; iter; iter = g_list_next (iter))
    {
      GType  type = (GType) iter->data;
      gchar *parasite_name = gimp_symmetry_parasite_name (type);

      parasite = gimp_image_parasite_find (image,
                                           parasite_name);
      g_free (parasite_name);
      if (parasite)
        {
          GimpSymmetry *sym = gimp_symmetry_from_parasite (parasite,
                                                           image,
                                                           type);

          if (sym)
            {
              GimpImagePrivate *private = GIMP_IMAGE_GET_PRIVATE (image);

              gimp_parasite_list_remove (private->parasites,
                                         gimp_parasite_get_name (parasite));

              gimp_image_symmetry_add (image, sym);

              g_signal_emit_by_name (sym, "active-changed", NULL);
              if (sym->active)
                gimp_image_set_active_symmetry (image, type);
            }
        }
    }
  g_list_free (syms);

  /* migrate the old "exif-data" parasite */
  parasite = gimp_image_parasite_find (GIMP_IMAGE (image),
                                       "exif-data");
  if (parasite)
    {
      GimpImagePrivate *private = GIMP_IMAGE_GET_PRIVATE (image);

      if (has_metadata)
        {
          g_printerr ("xcf-load: inconsistent metadata discovered: XCF file "
                      "has both 'gimp-image-metadata' and 'exif-data' "
                      "parasites, dropping old 'exif-data'\n");
        }
      else
        {
          GimpMetadata *metadata = gimp_image_get_metadata (image);
          GError       *my_error = NULL;
          const guchar *parasite_data;
          guint32       parasite_data_size;

          parasite_data = (const guchar *) gimp_parasite_get_data (parasite, &parasite_data_size);

          if (metadata)
            g_object_ref (metadata);
          else
            metadata = gimp_metadata_new ();

          if (! gimp_metadata_set_from_exif (metadata,
                                             parasite_data, parasite_data_size,
                                             &my_error))
            {
              gimp_message (gimp, G_OBJECT (info->progress),
                            GIMP_MESSAGE_WARNING,
                            _("Corrupt 'exif-data' parasite discovered.\n"
                              "Exif data could not be migrated: %s"),
                            my_error->message);
              g_clear_error (&my_error);
            }
          else
            {
              gimp_image_set_metadata (image, metadata, FALSE);
            }

          g_object_unref (metadata);
        }

      gimp_parasite_list_remove (private->parasites,
                                 gimp_parasite_get_name (parasite));
    }

  /* migrate the old "gimp-metadata" parasite */
  parasite = gimp_image_parasite_find (GIMP_IMAGE (image),
                                       "gimp-metadata");
  if (parasite)
    {
      GimpImagePrivate *private    = GIMP_IMAGE_GET_PRIVATE (image);
      const gchar      *xmp_data;
      guint32           xmp_length;

      xmp_data = (gchar *) gimp_parasite_get_data (parasite, &xmp_length);

      if (has_metadata)
        {
          g_printerr ("xcf-load: inconsistent metadata discovered: XCF file "
                      "has both 'gimp-image-metadata' and 'gimp-metadata' "
                      "parasites, dropping old 'gimp-metadata'\n");
        }
      else if (xmp_length < 14 ||
               strncmp (xmp_data, "GIMP_XMP_1", 10) != 0)
        {
          gimp_message (gimp, G_OBJECT (info->progress),
                        GIMP_MESSAGE_WARNING,
                        _("Corrupt 'gimp-metadata' parasite discovered.\n"
                          "XMP data could not be migrated."));
        }
      else
        {
          GimpMetadata *metadata = gimp_image_get_metadata (image);
          GError       *my_error = NULL;

          if (metadata)
            g_object_ref (metadata);
          else
            metadata = gimp_metadata_new ();

          if (! gimp_metadata_set_from_xmp (metadata,
                                            (const guint8 *) xmp_data + 10,
                                            xmp_length - 10,
                                            &my_error))
            {
              /* XMP metadata from 2.8.x or earlier can be really messed up.
               * Let's make the message more user friendly so they will
               * understand that we can't do anything about it.
               * See issue #987. */
              gimp_message (gimp, G_OBJECT (info->progress),
                            GIMP_MESSAGE_WARNING,
                            _("Corrupt XMP metadata saved by an older version of "
                              "GIMP could not be converted and will be ignored.\n"
                              "If you don't know what XMP is, you most likely don't "
                              "need it. Reported error: %s."),
                            my_error->message);
              g_clear_error (&my_error);
            }
          else
            {
              gimp_image_set_metadata (image, metadata, FALSE);
            }

          g_object_unref (metadata);
        }

      gimp_parasite_list_remove (private->parasites,
                                 gimp_parasite_get_name (parasite));
    }

  /* check for a gimp-xcf-compatibility-mode parasite */
  parasite = gimp_image_parasite_find (GIMP_IMAGE (image),
                                       "gimp-xcf-compatibility-mode");
  if (parasite)
    {
      GimpImagePrivate *private = GIMP_IMAGE_GET_PRIVATE (image);

      /* just ditch it, it's unused but shouldn't be re-saved */
      gimp_parasite_list_remove (private->parasites,
                                 gimp_parasite_get_name (parasite));
    }

  xcf_progress_update (info);

  while (TRUE)
    {
      GimpLayer *layer;
      GList     *item_path = NULL;

      /* read in the offset of the next layer */
      if (xcf_read_offset (info, &offset, 1) < info->bytes_per_offset)
        {
          GIMP_LOG (XCF, "Failed to read layer offset"
                    " at offset: %" G_GOFFSET_FORMAT, info->cp);
          break;
        }

      /* if the offset is 0 then we are at the end
       *  of the layer list.
       */
      if (offset == 0)
        break;

      /* save the current position as it is where the
       *  next layer offset is stored.
       */
      saved_pos = info->cp;

      if (offset < saved_pos)
        {
          GIMP_LOG (XCF, "Invalid layer offset: %" G_GOFFSET_FORMAT
                    " at offset: %" G_GOFFSET_FORMAT, offset, saved_pos);
          goto error;
        }

      /* seek to the layer offset */
      if (! xcf_seek_pos (info, offset, NULL))
        goto error;

      /* read in the layer */
      layer = xcf_load_layer (info, image, &item_path, &n_broken_effects);
      if (! layer)
        {
          n_broken_layers++;

          if (! xcf_seek_pos (info, saved_pos, NULL))
            {
              if (item_path)
                g_list_free (item_path);

              goto error;
            }

          /* Don't just stop at the first broken layer. Load as much as
           * possible.
           */
          if (! item_path)
            {
              GimpContainer *layers = gimp_image_get_layers (image);

              item_path = g_list_prepend (NULL,
                                          GUINT_TO_POINTER (gimp_container_get_n_children (layers)));

              broken_paths = g_list_prepend (broken_paths, item_path);
            }

          continue;
        }

      if (broken_paths && item_path)
        {
          /* Item paths may be a problem when layers are missing. */
          xcf_fix_item_path (layer, &item_path, broken_paths);
        }

      num_successful_elements++;

      xcf_progress_update (info);

      /* suspend layer-group size updates */
      if (GIMP_IS_GROUP_LAYER (layer))
        {
          GimpGroupLayer *group = GIMP_GROUP_LAYER (layer);

          group_layers = g_list_prepend (group_layers, group);

          gimp_group_layer_suspend_resize (group, FALSE);
        }

      /* add the layer to the image if its not the floating selection */
      if (layer != info->floating_sel)
        {
          GimpContainer *layers = gimp_image_get_layers (image);
          GimpContainer *container;
          GimpLayer     *parent;

          if (item_path)
            {
              if (info->floating_sel)
                {
                  /* there is a floating selection, but it will get
                   * added after all layers are loaded, so toplevel
                   * layer indices are off-by-one. Adjust item paths
                   * accordingly:
                   */
                  gint toplevel_index;

                  toplevel_index = GPOINTER_TO_UINT (item_path->data);

                  toplevel_index--;

                  item_path->data = GUINT_TO_POINTER (toplevel_index);
                }

              parent = GIMP_LAYER
                (gimp_item_stack_get_parent_by_path (GIMP_ITEM_STACK (layers),
                                                     item_path,
                                                     NULL));

              container = gimp_viewable_get_children (GIMP_VIEWABLE (parent));

              g_list_free (item_path);
            }
          else
            {
              parent    = NULL;
              container = layers;
            }

          gimp_image_add_layer (image, layer,
                                parent,
                                gimp_container_get_n_children (container),
                                FALSE);
        }

      /* restore the saved position so we'll be ready to
       *  read the next offset.
       */
      if (! xcf_seek_pos (info, saved_pos, NULL))
        goto error;
    }

  /* resume layer-group size updates, in reverse order */
  for (iter = group_layers; iter; iter = g_list_next (iter))
    {
      GimpGroupLayer *group = iter->data;

      gimp_group_layer_resume_resize (group, FALSE);
    }
  g_clear_pointer (&group_layers, g_list_free);

  if (broken_paths)
    {
      g_list_free_full (broken_paths, (GDestroyNotify) g_list_free);
      broken_paths = NULL;
    }

  while (TRUE)
    {
      GimpChannel *channel;

      /* read in the offset of the next channel */
      if (xcf_read_offset (info, &offset, 1) < info->bytes_per_offset)
        {
          GIMP_LOG (XCF, "Failed to read channel offset"
                    " at offset: %" G_GOFFSET_FORMAT, info->cp);
          break;
        }

      /* if the offset is 0 then we are at the end
       *  of the channel list.
       */
      if (offset == 0)
        break;

      /* save the current position as it is where the
       *  next channel offset is stored.
       */
      saved_pos = info->cp;

      if (offset < saved_pos)
        {
          GIMP_LOG (XCF, "Invalid channel offset: %" G_GOFFSET_FORMAT
                    " at offset: % "G_GOFFSET_FORMAT, offset, saved_pos);
          goto error;
        }

      /* seek to the channel offset */
      if (! xcf_seek_pos (info, offset, NULL))
        goto error;

      /* read in the channel */
      channel = xcf_load_channel (info, image);
      if (!channel)
        {
          n_broken_channels++;
          GIMP_LOG (XCF, "Failed to load channel.");

          if (! xcf_seek_pos (info, saved_pos, NULL))
            goto error;

          continue;
        }

      num_successful_elements++;

      xcf_progress_update (info);

      /* add the channel to the image if its not the selection */
      if (channel != gimp_image_get_mask (image))
        gimp_image_add_channel (image, channel,
                                NULL, /* FIXME tree */
                                gimp_container_get_n_children (gimp_image_get_channels (image)),
                                FALSE);

      /* restore the saved position so we'll be ready to
       *  read the next offset.
       */
      if (! xcf_seek_pos (info, saved_pos, NULL))
        goto error;
    }

  if (n_broken_layers == 0 && n_broken_channels == 0)
    {
      xcf_load_add_masks (image);
      xcf_load_add_effects (info, image);
    }

  if (info->floating_sel && info->floating_sel_drawable)
    {
      /* we didn't fix the loaded floating selection's format before
       * because we didn't know if it needed the layer space
       */
      if (GIMP_IS_LAYER (info->floating_sel_drawable) &&
          gimp_drawable_is_gray (GIMP_DRAWABLE (info->floating_sel)))
        gimp_layer_fix_format_space (info->floating_sel, TRUE, FALSE);

      floating_sel_attach (info->floating_sel, info->floating_sel_drawable);
    }

  if (info->file_version >= 18)
    {
      while (TRUE)
        {
          GimpPath *path;

          /* read in the offset of the next path */
          if (xcf_read_offset (info, &offset, 1) < info->bytes_per_offset)
            {
              GIMP_LOG (XCF, "Failed to read path offset"
                        " at offset: %" G_GOFFSET_FORMAT, info->cp);
              break;
            }

          /* if the offset is 0 then we are at the end
           *  of the path list.
           */
          if (offset == 0)
            break;

          /* save the current position as it is where the
           *  next channel offset is stored.
           */
          saved_pos = info->cp;

          if (offset < saved_pos)
            {
              GIMP_LOG (XCF, "Invalid path offset: %" G_GOFFSET_FORMAT
                        " at offset: % "G_GOFFSET_FORMAT, offset, saved_pos);
              goto error;
            }

          /* seek to the path offset */
          if (! xcf_seek_pos (info, offset, NULL))
            goto error;

          /* read in the path */
          path = xcf_load_path (info, image);
          if (! path)
            {
              n_broken_paths++;
              GIMP_LOG (XCF, "Failed to load path.");

              if (! xcf_seek_pos (info, saved_pos, NULL))
                goto error;

              continue;
            }

          num_successful_elements++;

          xcf_progress_update (info);

          gimp_image_add_path (image, path,
                               NULL, /* can't be a tree */
                               gimp_container_get_n_children (gimp_image_get_paths (image)),
                               FALSE);

          /* restore the saved position so we'll be ready to
           *  read the next offset.
           */
          if (! xcf_seek_pos (info, saved_pos, NULL))
            goto error;
        }
    }

  if (info->selected_layers)
    {
      gimp_image_set_selected_layers (image, info->selected_layers);
      g_clear_pointer (&info->selected_layers, g_list_free);
    }

  if (info->selected_channels)
    gimp_image_set_selected_channels (image, info->selected_channels);

  if (info->selected_paths)
    gimp_image_set_selected_paths (image, info->selected_paths);

  /* We don't have linked items concept anymore. We transform formerly
   * linked items into stored sets of named items instead.
   */
  if (info->linked_layers)
    {
      GimpItemList *set;

      set = gimp_item_list_named_new (image, GIMP_TYPE_LAYER,
                                      _("Linked Layers"),
                                      info->linked_layers);
      gimp_image_store_item_set (image, set);
      g_clear_pointer (&info->linked_layers, g_list_free);
    }
  if (info->linked_channels)
    {
      GimpItemList *set;

      set = gimp_item_list_named_new (image, GIMP_TYPE_CHANNEL,
                                      _("Linked Channels"),
                                      info->linked_channels);
      gimp_image_store_item_set (image, set);
      g_clear_pointer (&info->linked_channels, g_list_free);
    }
  if (info->linked_paths)
    {
      /* It is kind of ugly but paths are really implemented as
       * exception in our XCF spec and building over it seems like a
       * mistake. Since I'm seriously not sure this would be much of an
       * issue, I'll let it as it for now.
       * Note that it's still possible to multi-select paths. It's only
       * not possible to store these selections.
       *
       * Only warn for more than 1 linked path. Less is kind of
       * pointless and doesn't deserve worrying people for no reason.
       */
      if (g_list_length (info->linked_paths) > 1)
        g_printerr ("xcf: some paths were linked. "
                    "GIMP does not support linked paths since version 3.0.\n");

#if 0
      GimpItemList *set;

      set = gimp_item_list_named_new (image, GIMP_TYPE_PATH,
                                      _("Linked Paths"),
                                      info->linked_paths);
      gimp_image_store_item_set (image, set);
#endif
      g_clear_pointer (&info->linked_paths, g_list_free);
    }

  for (iter = g_list_last (info->layer_sets); iter; iter = iter->prev)
    {
      if (iter->data)
        gimp_image_store_item_set (image, iter->data);
    }
  g_list_free (info->layer_sets);

  for (iter = g_list_last (info->channel_sets); iter; iter = iter->prev)
    {
      if (iter->data)
        gimp_image_store_item_set (image, iter->data);
    }
  g_list_free (info->channel_sets);

  if (info->file)
    gimp_image_set_file (image, info->file);

  if (info->tattoo_state > 0)
    gimp_image_set_tattoo_state (image, info->tattoo_state);

  if (n_broken_layers > 0 || n_broken_channels > 0 || n_broken_paths > 0 || n_broken_effects)
    goto error;

  gimp_image_undo_enable (image);

  return image;

 error:
  if (num_successful_elements == 0)
    goto hard_error;

  g_clear_pointer (&group_layers, g_list_free);

  if (broken_paths)
    {
      g_list_free_full (broken_paths, (GDestroyNotify) g_list_free);
      broken_paths = NULL;
    }

  gimp_message_literal (gimp, G_OBJECT (info->progress), GIMP_MESSAGE_WARNING,
                        _("This XCF file is corrupt!  I have loaded as much "
                          "of it as I can, but it is incomplete."));

  xcf_load_add_masks (image);
  xcf_load_add_effects (info, image);

  gimp_image_undo_enable (image);

  return image;

 hard_error:
  g_clear_pointer (&group_layers, g_list_free);

  if (broken_paths)
    {
      g_list_free_full (broken_paths, (GDestroyNotify) g_list_free);
      broken_paths = NULL;
    }

  g_set_error_literal (error, G_FILE_ERROR, G_FILE_ERROR_FAILED,
                       _("This XCF file is corrupt!  I could not even "
                         "salvage any partial image data from it."));

  g_clear_object (&image);

  return NULL;
}

static void
xcf_load_add_masks (GimpImage *image)
{
  GList *layers;
  GList *list;

  layers = gimp_image_get_layer_list (image);

  for (list = layers; list; list = g_list_next (list))
    {
      GimpLayer     *layer = list->data;
      GimpLayerMask *mask;

      mask = g_object_get_data (G_OBJECT (layer), "gimp-layer-mask");

      if (mask)
        {
          gboolean apply_mask;
          gboolean edit_mask;
          gboolean show_mask;

          apply_mask = GPOINTER_TO_INT (g_object_get_data (G_OBJECT (layer),
                                                           "gimp-layer-mask-apply"));
          edit_mask = GPOINTER_TO_INT (g_object_get_data (G_OBJECT (layer),
                                                          "gimp-layer-mask-edit"));
          show_mask = GPOINTER_TO_INT (g_object_get_data (G_OBJECT (layer),
                                                          "gimp-layer-mask-show"));

          gimp_layer_add_mask (layer, mask, FALSE, NULL);

          gimp_layer_set_apply_mask (layer, apply_mask, FALSE);
          gimp_layer_set_edit_mask  (layer, edit_mask);
          gimp_layer_set_show_mask  (layer, show_mask, FALSE);

          g_object_set_data (G_OBJECT (layer), "gimp-layer-mask",       NULL);
          g_object_set_data (G_OBJECT (layer), "gimp-layer-mask-apply", NULL);
          g_object_set_data (G_OBJECT (layer), "gimp-layer-mask-edit",  NULL);
          g_object_set_data (G_OBJECT (layer), "gimp-layer-mask-show",  NULL);
        }
    }

  g_list_free (layers);
}

static void
xcf_load_add_effects (XcfInfo   *info,
                      GimpImage *image)
{
  GList *layers;
  GList *list;

  layers = gimp_image_get_layer_list (image);

  for (list = layers; list; list = g_list_next (list))
    {
      GimpLayer *layer = list->data;
      GList     *effects_nodes;

      effects_nodes = g_object_get_data (G_OBJECT (layer), "gimp-layer-effects");

      if (effects_nodes)
        {
          GList *iter;

          for (iter = effects_nodes; iter; iter = iter->next)
            {
              FilterData *data = iter->data;

              if (! data->unsupported_operation)
                {
                  GimpDrawableFilter *filter = NULL;

                  if (! data->icon_name)
                    data->icon_name = g_strdup ("gimp-gegl");

                  filter = gimp_drawable_filter_new (GIMP_DRAWABLE (layer),
                                                     data->name, data->operation,
                                                     data->icon_name);

                  gimp_drawable_filter_set_opacity (filter, data->opacity);
                  gimp_drawable_filter_set_mode (filter, data->paint_mode,
                                                 data->blend_space,
                                                 data->composite_space,
                                                 data->composite_mode);
                  gimp_drawable_filter_set_clip (filter, data->clip);
                  gimp_drawable_filter_set_region (filter, data->region);

                  g_object_set (filter,
                                "mask", data->mask,
                                NULL);

                  gimp_drawable_filter_apply_with_mask (filter, data->mask,
                                                        NULL);

                  gimp_drawable_filter_commit (filter, TRUE, NULL, FALSE);

                  gimp_drawable_filter_layer_mask_freeze (filter);

                  gimp_filter_set_active (GIMP_FILTER (filter), data->is_visible);

                  g_object_unref (filter);
                }
            }

          g_object_set_data (G_OBJECT (layer), "gimp-layer-effects", NULL);
        }
    }

  g_list_free (layers);
}

static gboolean
xcf_load_image_props (XcfInfo   *info,
                      GimpImage *image)
{
  PropType prop_type;
  guint32  prop_size;

  while (TRUE)
    {
      if (! xcf_load_prop (info, &prop_type, &prop_size))
        return FALSE;

      switch (prop_type)
        {
        case PROP_END:
          return TRUE;

        case PROP_COLORMAP:
          {
            guint32 n_colors;
            guchar  cmap[GIMP_IMAGE_COLORMAP_SIZE];

            xcf_read_int32 (info, &n_colors, 1);

            if (n_colors > (GIMP_IMAGE_COLORMAP_SIZE / 3))
              {
                gimp_message (info->gimp, G_OBJECT (info->progress),
                              GIMP_MESSAGE_ERROR,
                              "Maximum colormap size (%d) exceeded",
                              GIMP_IMAGE_COLORMAP_SIZE);
                return FALSE;
              }

            if (info->file_version == 0)
              {
                gint i;

                gimp_message_literal (info->gimp, G_OBJECT (info->progress),
                                      GIMP_MESSAGE_WARNING,
                                      _("XCF warning: version 0 of XCF file format\n"
                                        "did not save indexed colormaps correctly.\n"
                                        "Substituting grayscale map."));

                if (! xcf_seek_pos (info, info->cp + n_colors, NULL))
                  return FALSE;

                for (i = 0; i < n_colors; i++)
                  {
                    cmap[i * 3 + 0] = i;
                    cmap[i * 3 + 1] = i;
                    cmap[i * 3 + 2] = i;
                  }
              }
            else
              {
                xcf_read_int8 (info, cmap, n_colors * 3);
              }

            /* only set color map if image is indexed, this is just
             * sanity checking to make sure gimp doesn't end up with
             * an image state that is impossible.
             */
            if (gimp_image_get_base_type (image) == GIMP_INDEXED)
              _gimp_image_set_colormap (image, cmap, n_colors, FALSE);

            GIMP_LOG (XCF, "prop colormap n_colors=%d", n_colors);
          }
          break;

        case PROP_COMPRESSION:
          {
            guint8 compression;

            xcf_read_int8 (info, (guint8 *) &compression, 1);

            if ((compression != COMPRESS_NONE) &&
                (compression != COMPRESS_RLE) &&
                (compression != COMPRESS_ZLIB) &&
                (compression != COMPRESS_FRACTAL))
              {
                gimp_message (info->gimp, G_OBJECT (info->progress),
                              GIMP_MESSAGE_ERROR,
                              "Unknown compression type: %d",
                              (gint) compression);
                return FALSE;
              }

            info->compression = compression;

            gimp_image_set_xcf_compression (image,
                                            compression >= COMPRESS_ZLIB);

            GIMP_LOG (XCF, "prop compression=%d", compression);
          }
          break;

        case PROP_GUIDES:
          {
            GimpImagePrivate *private = GIMP_IMAGE_GET_PRIVATE (image);
            gint32            position;
            gint8             orientation;
            gint              i, nguides;

            nguides = prop_size / (4 + 1);
            for (i = 0; i < nguides; i++)
              {
                xcf_read_int32 (info, (guint32 *) &position,    1);
                xcf_read_int8  (info, (guint8 *)  &orientation, 1);

                /* Some very old XCF had -1 guides which have been
                 * skipped since 2003 (commit 909a28ced2).
                 * Then XCF up to version 14 only had positive guide
                 * positions.
                 * Since XCF 15 (GIMP 3.0), off-canvas guides became a
                 * thing.
                 */
                if (info->file_version < 15 && position < 0)
                  continue;

                GIMP_LOG (XCF, "prop guide orientation=%d position=%d",
                          orientation, position);

                switch (orientation)
                  {
                  case XCF_ORIENTATION_HORIZONTAL:
                    if (info->file_version < 15 && position > gimp_image_get_height (image))
                      gimp_message (info->gimp, G_OBJECT (info->progress),
                                    GIMP_MESSAGE_WARNING,
                                    "Ignoring off-canvas horizontal guide (position %d) in XCF %d file",
                                    position, info->file_version);
                    else
                      gimp_image_add_hguide (image, position, FALSE);
                    break;

                  case XCF_ORIENTATION_VERTICAL:
                    if (info->file_version < 15 && position > gimp_image_get_width (image))
                      gimp_message (info->gimp, G_OBJECT (info->progress),
                                    GIMP_MESSAGE_WARNING,
                                    "Ignoring off-canvas vertical guide (position %d) in XCF %d file",
                                    position, info->file_version);
                    else
                      gimp_image_add_vguide (image, position, FALSE);
                    break;

                  default:
                    gimp_message_literal (info->gimp, G_OBJECT (info->progress),
                                          GIMP_MESSAGE_WARNING,
                                          "Guide orientation out of range in XCF file");
                    continue;
                  }
              }

            /*  this is silly as the order of guides doesn't really matter,
             *  but it restores the list to its original order, which
             *  cannot be wrong  --Mitch
             */
            private->guides = g_list_reverse (private->guides);
          }
          break;

        case PROP_SAMPLE_POINTS:
          {
            gint n_sample_points, i;

            n_sample_points = prop_size / (5 * 4);
            for (i = 0; i < n_sample_points; i++)
              {
                GimpSamplePoint   *sample_point;
                gint32             x, y;
                GimpColorPickMode  pick_mode;
                guint32            padding[2] = { 0, };

                xcf_read_int32 (info, (guint32 *) &x,         1);
                xcf_read_int32 (info, (guint32 *) &y,         1);
                xcf_read_int32 (info, (guint32 *) &pick_mode, 1);
                xcf_read_int32 (info, (guint32 *) padding,    2);

                GIMP_LOG (XCF, "prop sample point x=%d y=%d mode=%d",
                          x, y, pick_mode);

                if (pick_mode > GIMP_COLOR_PICK_MODE_LAST)
                  pick_mode = GIMP_COLOR_PICK_MODE_PIXEL;

                sample_point = gimp_image_add_sample_point_at_pos (image,
                                                                   x, y, FALSE);
                gimp_image_set_sample_point_pick_mode (image, sample_point,
                                                       pick_mode, FALSE);
              }
          }
          break;

        case PROP_OLD_SAMPLE_POINTS:
          {
            gint32 x, y;
            gint   i, n_sample_points;

            /* if there are already sample points, we loaded the new
             * prop before
             */
            if (gimp_image_get_sample_points (image))
              {
                if (! xcf_skip_unknown_prop (info, prop_size))
                  return FALSE;

                break;
              }

            n_sample_points = prop_size / (4 + 4);
            for (i = 0; i < n_sample_points; i++)
              {
                xcf_read_int32 (info, (guint32 *) &x, 1);
                xcf_read_int32 (info, (guint32 *) &y, 1);

                GIMP_LOG (XCF, "prop old sample point x=%d y=%d", x, y);

                gimp_image_add_sample_point_at_pos (image, x, y, FALSE);
              }
          }
          break;

        case PROP_RESOLUTION:
          {
            gfloat xres, yres;

            xcf_read_float (info, &xres, 1);
            xcf_read_float (info, &yres, 1);

            GIMP_LOG (XCF, "prop resolution x=%f y=%f", xres, yres);

            if (xres < GIMP_MIN_RESOLUTION || xres > GIMP_MAX_RESOLUTION ||
                yres < GIMP_MIN_RESOLUTION || yres > GIMP_MAX_RESOLUTION)
              {
                GimpTemplate *template = image->gimp->config->default_image;

                gimp_message_literal (info->gimp, G_OBJECT (info->progress),
                                      GIMP_MESSAGE_WARNING,
                                      "Warning, resolution out of range in XCF file");
                xres = gimp_template_get_resolution_x (template);
                yres = gimp_template_get_resolution_y (template);
              }

            gimp_image_set_resolution (image, xres, yres);
          }
          break;

        case PROP_TATTOO:
          {
            xcf_read_int32 (info, &info->tattoo_state, 1);

            GIMP_LOG (XCF, "prop tattoo state=%d", info->tattoo_state);
          }
          break;

        case PROP_PARASITES:
          {
            goffset base = info->cp;

            while (info->cp - base < prop_size)
              {
                GimpParasite *p     = xcf_load_parasite (info);
                GError       *error = NULL;

                if (! p)
                  {
                    gimp_message (info->gimp, G_OBJECT (info->progress),
                                  GIMP_MESSAGE_WARNING,
                                  "Invalid image parasite found. "
                                  "Possibly corrupt XCF file.");

                    xcf_seek_pos (info, base + prop_size, NULL);
                    continue;
                  }

                if (! gimp_image_parasite_validate (image, p, &error))
                  {
                    gimp_message (info->gimp, G_OBJECT (info->progress),
                                  GIMP_MESSAGE_WARNING,
                                  "Warning, invalid image parasite in XCF file: %s",
                                  error->message);
                    g_clear_error (&error);
                  }
                else
                  {
                    gimp_image_parasite_attach (image, p, FALSE);
                  }

                gimp_parasite_free (p);
              }

            if (info->cp - base != prop_size)
              gimp_message_literal (info->gimp, G_OBJECT (info->progress),
                                    GIMP_MESSAGE_WARNING,
                                    "Error while loading an image's parasites");
          }
          break;

        case PROP_UNIT:
          {
            guint32 unit_index;

            xcf_read_int32 (info, &unit_index, 1);

            GIMP_LOG (XCF, "prop unit=%d", unit_index);

            if (unit_index <= GIMP_UNIT_PIXEL || unit_index >= GIMP_UNIT_END)
              {
                gimp_message_literal (info->gimp, G_OBJECT (info->progress),
                                      GIMP_MESSAGE_WARNING,
                                      "Warning, unit out of range in XCF file, "
                                      "falling back to inches");
                unit_index = GIMP_UNIT_INCH;
              }

            gimp_image_set_unit (image, gimp_unit_get_by_id (unit_index));
          }
          break;

        case PROP_PATHS:
          {
            goffset base = info->cp;

            if (info->file_version >= 18)
              gimp_message (info->gimp, G_OBJECT (info->progress),
                            GIMP_MESSAGE_WARNING,
                            "XCF %d file should not contain PROP_PATHS image properties",
                            info->file_version);

            if (! xcf_load_old_paths (info, image))
              xcf_seek_pos (info, base + prop_size, NULL);
          }
          break;

        case PROP_USER_UNIT:
          {
            gchar     *unit_strings[5] = { 0 };
            float      factor;
            guint32    digits;
            GimpUnit  *unit;
            GList     *iter;
            gint       n_fields = 3;
            gint       i;

            xcf_read_float  (info, &factor,      1);
            xcf_read_int32  (info, &digits,      1);

            /* Depending on XCF version, read more or less strings. */
            if (info->file_version < 21)
              n_fields = 5;
            xcf_read_string (info, unit_strings, n_fields);

            for (i = 0; i < n_fields; i++)
              if (unit_strings[i] == NULL)
                unit_strings[i] = g_strdup ("");

            for (iter = info->gimp->user_units; iter; iter = iter->next)
              {
                unit = iter->data;
                /* if the factor and the name match some unit in unitrc,
                 * use the unitrc unit
                 */
                if (ABS (gimp_unit_get_factor (unit) - factor) < 1e-5 &&
                    (strcmp (unit_strings[0], gimp_unit_get_name (unit)) == 0 ||
                     (info->file_version < 21 &&
                      strcmp (unit_strings[4], gimp_unit_get_name (unit)) == 0)))
                  {
                    break;
                  }
              }

            if (iter == NULL)
              /* No match. Create a temporary unit set with deletion
               * flag.
               */
              unit = _gimp_unit_new (info->gimp,
                                     unit_strings[4] && strlen (unit_strings[4]) > 0 ? unit_strings[4] : unit_strings[0],
                                     (gdouble) factor,
                                     digits,
                                     unit_strings[1],
                                     unit_strings[2]);

            gimp_image_set_unit (image, unit);

            for (i = 0; i < n_fields; i++)
              g_free (unit_strings[i]);
          }
         break;

        case PROP_VECTORS:
          {
            goffset base = info->cp;

            if (info->file_version >= 18)
              gimp_message (info->gimp, G_OBJECT (info->progress),
                            GIMP_MESSAGE_WARNING,
                            "XCF %d file should not contain PROP_VECTORS image properties",
                            info->file_version);

            if (xcf_load_old_vectors (info, image))
              {
                if (base + prop_size != info->cp)
                  {
                    g_printerr ("Mismatch in PROP_VECTORS size: "
                                "skipping %" G_GOFFSET_FORMAT " bytes.\n",
                                base + prop_size - info->cp);
                    xcf_seek_pos (info, base + prop_size, NULL);
                  }
              }
            else
              {
                /* skip silently since we don't understand the format and
                 * xcf_load_old_vectors already explained what was wrong
                 */
                xcf_seek_pos (info, base + prop_size, NULL);
              }
          }
          break;

        case PROP_ITEM_SET:
          {
            GimpItemList *set       = NULL;
            gchar        *label;
            GType         item_type = 0;
            guint32       itype;
            guint32       method;

            xcf_read_int32  (info, &itype, 1);
            xcf_read_int32  (info, &method, 1);
            xcf_read_string (info, &label, 1);

            if (itype == 0)
              item_type = GIMP_TYPE_LAYER;
            else
              item_type = GIMP_TYPE_CHANNEL;

            if (itype > 1)
              {
                g_printerr ("xcf: unsupported item set '%s' type: %d (skipping)\n",
                            label ? label : "unnamed", itype);
                /* Only case where we break because we wouldn't even
                 * know where to categorize the item set anyway. */
                break;
              }
            else if (label == NULL)
              {
                g_printerr ("xcf: item set without a name or pattern (skipping)\n");
              }
            else if (method != G_MAXUINT32 && method > GIMP_SELECT_GLOB_PATTERN)
              {
                g_printerr ("xcf: unsupported item set '%s' selection method attribute: 0x%x (skipping)\n",
                            label, method);
              }
            else
              {
                if (method == G_MAXUINT32)
                  {
                    /* Don't use gimp_item_list_named_new() because it
                     * doesn't allow NULL items (it would try to get the
                     * selected items instead).
                     */
                    set = g_object_new (GIMP_TYPE_ITEM_LIST,
                                        "image",      image,
                                        "name",       label,
                                        "is-pattern", FALSE,
                                        "item-type",  item_type,
                                        "items",      NULL,
                                        NULL);
                  }
                else
                  {
                    set = gimp_item_list_pattern_new (image, item_type,
                                                      method, label);
                  }
              }

            /* Note: we are still adding invalid item sets as NULL on
             * purpose, in order not to break order-base association
             * between PROP_ITEM_SET and PROP_ITEM_SET_ITEM.
             */
            if (item_type == GIMP_TYPE_LAYER)
              info->layer_sets = g_list_prepend (info->layer_sets, set);
            else
              info->channel_sets = g_list_prepend (info->channel_sets, set);
          }
          break;

        default:
#ifdef GIMP_UNSTABLE
          g_printerr ("unexpected/unknown image property: %d (skipping)\n",
                      prop_type);
#endif
          if (! xcf_skip_unknown_prop (info, prop_size))
            return FALSE;
          break;
        }
    }

  return FALSE;
}

static gboolean
xcf_load_layer_props (XcfInfo    *info,
                      GimpImage  *image,
                      GimpLayer **layer,
                      GList     **item_path,
                      gboolean   *apply_mask,
                      gboolean   *edit_mask,
                      gboolean   *show_mask,
                      guint32    *text_layer_flags,
                      guint32    *group_layer_flags)
{
  PropType prop_type;
  guint32  prop_size;

  while (TRUE)
    {
      if (! xcf_load_prop (info, &prop_type, &prop_size))
        return FALSE;

      switch (prop_type)
        {
        case PROP_END:
          return TRUE;

        case PROP_ACTIVE_LAYER:
          {
            if (g_list_index (info->selected_layers, *layer) < 0)
              info->selected_layers = g_list_prepend (info->selected_layers, *layer);
            else
              gimp_message_literal (info->gimp, G_OBJECT (info->progress),
                                    GIMP_MESSAGE_WARNING,
                                    "Invalid duplicate selected layer detected");
          }
          break;

        case PROP_FLOATING_SELECTION:
          info->floating_sel = *layer;
          xcf_read_offset (info, &info->floating_sel_offset, 1);
          break;

        case PROP_OPACITY:
          {
            guint32 opacity;

            xcf_read_int32 (info, &opacity, 1);

            gimp_layer_set_opacity (*layer, (gdouble) opacity / 255.0, FALSE);
          }
          break;

        case PROP_FLOAT_OPACITY:
          {
            gfloat opacity;

            xcf_read_float (info, &opacity, 1);

            gimp_layer_set_opacity (*layer, opacity, FALSE);
          }
          break;

        case PROP_VISIBLE:
          {
            gboolean visible;

            xcf_read_int32 (info, (guint32 *) &visible, 1);

            gimp_item_set_visible (GIMP_ITEM (*layer), visible, FALSE);
          }
          break;

        case PROP_LINKED:
          {
            gboolean linked;

            xcf_read_int32 (info, (guint32 *) &linked, 1);

            if (linked)
              info->linked_layers = g_list_prepend (info->linked_layers, *layer);
          }
          break;

        case PROP_COLOR_TAG:
          {
            GimpColorTag color_tag;

            xcf_read_int32 (info, (guint32 *) &color_tag, 1);

            gimp_item_set_color_tag (GIMP_ITEM (*layer), color_tag, FALSE);
          }
          break;

        case PROP_LOCK_CONTENT:
          {
            gboolean lock_content;

            xcf_read_int32 (info, (guint32 *) &lock_content, 1);

            if (gimp_item_can_lock_content (GIMP_ITEM (*layer)))
              gimp_item_set_lock_content (GIMP_ITEM (*layer),
                                          lock_content, FALSE);
          }
          break;

        case PROP_LOCK_ALPHA:
          {
            gboolean lock_alpha;

            xcf_read_int32 (info, (guint32 *) &lock_alpha, 1);

            if (gimp_layer_can_lock_alpha (*layer))
              gimp_layer_set_lock_alpha (*layer, lock_alpha, FALSE);
          }
          break;

        case PROP_LOCK_POSITION:
          {
            gboolean lock_position;

            xcf_read_int32 (info, (guint32 *) &lock_position, 1);

            if (gimp_item_can_lock_position (GIMP_ITEM (*layer)))
              gimp_item_set_lock_position (GIMP_ITEM (*layer),
                                           lock_position, FALSE);
          }
          break;

        case PROP_LOCK_VISIBILITY:
          {
            gboolean lock_visibility;

            xcf_read_int32 (info, (guint32 *) &lock_visibility, 1);

            if (gimp_item_can_lock_visibility (GIMP_ITEM (*layer)))
              gimp_item_set_lock_visibility (GIMP_ITEM (*layer),
                                             lock_visibility, FALSE);
          }
          break;

        case PROP_APPLY_MASK:
          xcf_read_int32 (info, (guint32 *) apply_mask, 1);
          break;

        case PROP_EDIT_MASK:
          xcf_read_int32 (info, (guint32 *) edit_mask, 1);
          break;

        case PROP_SHOW_MASK:
          xcf_read_int32 (info, (guint32 *) show_mask, 1);
          break;

        case PROP_OFFSETS:
          {
            gint32 offset_x;
            gint32 offset_y;

            xcf_read_int32 (info, (guint32 *) &offset_x, 1);
            xcf_read_int32 (info, (guint32 *) &offset_y, 1);

            if (offset_x < -GIMP_MAX_IMAGE_SIZE ||
                offset_x > GIMP_MAX_IMAGE_SIZE)
              {
                g_printerr ("unexpected item offset_x (%d) in XCF, "
                            "setting to 0\n", offset_x);
                offset_x = 0;
              }

            if (offset_y < -GIMP_MAX_IMAGE_SIZE ||
                offset_y > GIMP_MAX_IMAGE_SIZE)
              {
                g_printerr ("unexpected item offset_y (%d) in XCF, "
                            "setting to 0\n", offset_y);
                offset_y = 0;
              }

            gimp_item_set_offset (GIMP_ITEM (*layer), offset_x, offset_y);
          }
          break;

        case PROP_MODE:
          {
            GimpLayerMode mode;

            xcf_read_int32 (info, (guint32 *) &mode, 1);

            if (mode == GIMP_LAYER_MODE_OVERLAY_LEGACY)
              mode = GIMP_LAYER_MODE_SOFTLIGHT_LEGACY;

            gimp_layer_set_mode (*layer, mode, FALSE);
          }
          break;

        case PROP_BLEND_SPACE:
          {
            gint32 blend_space;

            xcf_read_int32 (info, (guint32 *) &blend_space, 1);

            /* if blend_space < 0 it was originally AUTO, and its negative is
             * the actual value AUTO used to map to at the time the file was
             * saved.  if AUTO still maps to the same value, keep using AUTO
             * for the property; otherwise, use the concrete value.
             */
            if (blend_space < 0)
              {
                GimpLayerMode mode = gimp_layer_get_mode (*layer);

                blend_space = -blend_space;

                if (blend_space == gimp_layer_mode_get_blend_space (mode))
                  blend_space = GIMP_LAYER_COLOR_SPACE_AUTO;
                else
                  GIMP_LOG (XCF, "BLEND_SPACE: AUTO => %d", blend_space);
              }

            gimp_layer_set_blend_space (*layer, blend_space, FALSE);
          }
          break;

        case PROP_COMPOSITE_SPACE:
          {
            gint32 composite_space;

            xcf_read_int32 (info, (guint32 *) &composite_space, 1);

            /* if composite_space < 0 it was originally AUTO, and its negative
             * is the actual value AUTO used to map to at the time the file was
             * saved.  if AUTO still maps to the same value, keep using AUTO
             * for the property; otherwise, use the concrete value.
             */
            if (composite_space < 0)
              {
                GimpLayerMode mode = gimp_layer_get_mode (*layer);

                composite_space = -composite_space;

                if (composite_space == gimp_layer_mode_get_composite_space (mode))
                  composite_space = GIMP_LAYER_COLOR_SPACE_AUTO;
                else
                  GIMP_LOG (XCF, "COMPOSITE_SPACE: AUTO => %d", composite_space);
              }

            gimp_layer_set_composite_space (*layer, composite_space, FALSE);
          }
          break;

        case PROP_COMPOSITE_MODE:
          {
            gint32 composite_mode;

            xcf_read_int32 (info, (guint32 *) &composite_mode, 1);

            /* if composite_mode < 0 it was originally AUTO, and its negative
             * is the actual value AUTO used to map to at the time the file was
             * saved.  if AUTO still maps to the same value, keep using AUTO
             * for the property; otherwise, use the concrete value.
             */
            if (composite_mode < 0)
              {
                GimpLayerMode mode = gimp_layer_get_mode (*layer);

                composite_mode = -composite_mode;

                if (composite_mode == gimp_layer_mode_get_composite_mode (mode))
                  composite_mode = GIMP_LAYER_COMPOSITE_AUTO;
                else
                  GIMP_LOG (XCF, "COMPOSITE_MODE: AUTO => %d", composite_mode);
              }

            gimp_layer_set_composite_mode (*layer, composite_mode, FALSE);
          }
          break;

        case PROP_TATTOO:
          {
            GimpTattoo tattoo;

            xcf_read_int32 (info, (guint32 *) &tattoo, 1);

            gimp_item_set_tattoo (GIMP_ITEM (*layer), tattoo);
          }
          break;

        case PROP_PARASITES:
          {
            goffset base = info->cp;

            while (info->cp - base < prop_size)
              {
                GimpParasite *p     = xcf_load_parasite (info);
                GError       *error = NULL;

                if (! p)
                  return FALSE;

                if (! gimp_item_parasite_validate (GIMP_ITEM (*layer), p, &error))
                  {
                    gimp_message (info->gimp, G_OBJECT (info->progress),
                                  GIMP_MESSAGE_WARNING,
                                  "Warning, invalid layer parasite in XCF file: %s",
                                  error->message);
                    g_clear_error (&error);
                  }
                else
                  {
                    gimp_item_parasite_attach (GIMP_ITEM (*layer), p, FALSE);
                  }

                gimp_parasite_free (p);
              }

            if (info->cp - base != prop_size)
              gimp_message_literal (info->gimp, G_OBJECT (info->progress),
                                    GIMP_MESSAGE_WARNING,
                                    "Error while loading a layer's parasites");
          }
          break;

        case PROP_TEXT_LAYER_FLAGS:
          xcf_read_int32 (info, text_layer_flags, 1);
          break;

        case PROP_GROUP_ITEM:
          {
            GimpLayer *group;
            gboolean   is_selected_layer;

            /* We're going to delete *layer, Don't leave its pointers
             * in @info.  After that, we'll restore them back with the
             * new pointer. See bug #767873.
             */
            is_selected_layer = (g_list_find (info->selected_layers, *layer ) != NULL);
            if (is_selected_layer)
              info->selected_layers = g_list_remove (info->selected_layers, *layer);

            if (*layer == info->floating_sel)
              info->floating_sel = NULL;

            group = gimp_group_layer_new (image);

            gimp_object_set_name (GIMP_OBJECT (group),
                                  gimp_object_get_name (*layer));

            g_object_ref_sink (*layer);
            g_object_unref (*layer);
            *layer = group;

            if (is_selected_layer)
              info->selected_layers = g_list_prepend (info->selected_layers, *layer);

            /* Don't restore info->floating_sel because group layers
             * can't be floating selections
             */
          }
          break;

        case PROP_ITEM_PATH:
          {
            goffset  base = info->cp;
            GList   *path = NULL;

            while (info->cp - base < prop_size)
              {
                guint32 index;

                if (xcf_read_int32 (info, &index, 1) != 4)
                  {
                    g_list_free (path);
                    return FALSE;
                  }

                path = g_list_append (path, GUINT_TO_POINTER (index));
              }

            *item_path = path;
          }
          break;

        case PROP_GROUP_ITEM_FLAGS:
          xcf_read_int32 (info, group_layer_flags, 1);
          break;

        case PROP_ITEM_SET_ITEM:
            {
              GimpItemList *set;
              guint32       n;

              xcf_read_int32 (info, &n, 1);
              set = g_list_nth_data (info->layer_sets, n);
              if (set == NULL)
                g_printerr ("xcf: layer '%s' cannot be added to unknown layer set at index %d (skipping)\n",
                            gimp_object_get_name (*layer), n);
              else if (! g_type_is_a (G_TYPE_FROM_INSTANCE (*layer),
                                      gimp_item_list_get_item_type (set)))
                g_printerr ("xcf: layer '%s' cannot be added to item set '%s' with item type %s (skipping)\n",
                            gimp_object_get_name (*layer), gimp_object_get_name (set),
                            g_type_name (gimp_item_list_get_item_type (set)));
              else if (gimp_item_list_is_pattern (set, NULL))
                g_printerr ("xcf: layer '%s' cannot be added to pattern item set '%s' (skipping)\n",
                            gimp_object_get_name (*layer), gimp_object_get_name (set));
              else
                gimp_item_list_add (set, GIMP_ITEM (*layer));
            }
          break;

        default:
#ifdef GIMP_UNSTABLE
          g_printerr ("unexpected/unknown layer property: %d (skipping)\n",
                      prop_type);
#endif
          if (! xcf_skip_unknown_prop (info, prop_size))
            return FALSE;
          break;
        }
    }

  return FALSE;
}

static gboolean
xcf_check_layer_props (XcfInfo    *info,
                       GList     **item_path,
                       gboolean   *is_group_layer,
                       gboolean   *is_text_layer)
{
  PropType prop_type;
  guint32  prop_size;

  g_return_val_if_fail (*is_group_layer == FALSE, FALSE);
  g_return_val_if_fail (*is_text_layer  == FALSE, FALSE);

  while (TRUE)
    {
      if (! xcf_load_prop (info, &prop_type, &prop_size))
        return FALSE;

      switch (prop_type)
        {
        case PROP_END:
          return TRUE;

        case PROP_TEXT_LAYER_FLAGS:
          *is_text_layer = TRUE;

          if (! xcf_skip_unknown_prop (info, prop_size))
            return FALSE;
          break;

        case PROP_GROUP_ITEM:
        case PROP_GROUP_ITEM_FLAGS:
          *is_group_layer = TRUE;

          if (! xcf_skip_unknown_prop (info, prop_size))
            return FALSE;
          break;

        case PROP_ITEM_PATH:
          {
            goffset  base = info->cp;
            GList   *path = NULL;

            while (info->cp - base < prop_size)
              {
                guint32 index;

                if (xcf_read_int32 (info, &index, 1) != 4)
                  {
                    g_list_free (path);
                    return FALSE;
                  }

                path = g_list_append (path, GUINT_TO_POINTER (index));
              }

            *item_path = path;
          }
          break;

        case PROP_ACTIVE_LAYER:
        case PROP_FLOATING_SELECTION:
        case PROP_OPACITY:
        case PROP_FLOAT_OPACITY:
        case PROP_VISIBLE:
        case PROP_LINKED:
        case PROP_COLOR_TAG:
        case PROP_LOCK_CONTENT:
        case PROP_LOCK_ALPHA:
        case PROP_LOCK_POSITION:
        case PROP_LOCK_VISIBILITY:
        case PROP_APPLY_MASK:
        case PROP_EDIT_MASK:
        case PROP_SHOW_MASK:
        case PROP_OFFSETS:
        case PROP_MODE:
        case PROP_BLEND_SPACE:
        case PROP_COMPOSITE_SPACE:
        case PROP_COMPOSITE_MODE:
        case PROP_TATTOO:
        case PROP_PARASITES:
        case PROP_ITEM_SET_ITEM:
          if (! xcf_skip_unknown_prop (info, prop_size))
            return FALSE;
          /* Just ignore for now. */
          break;

        default:
#ifdef GIMP_UNSTABLE
          g_printerr ("unexpected/unknown layer property: %d (skipping)\n",
                      prop_type);
#endif
          if (! xcf_skip_unknown_prop (info, prop_size))
            return FALSE;
          break;
        }
    }

  return FALSE;
}

static gboolean
xcf_load_channel_props (XcfInfo      *info,
                        GimpImage    *image,
                        GimpChannel **channel,
                        gboolean      is_mask)
{
  PropType prop_type;
  guint32  prop_size;

  while (TRUE)
    {
      if (! xcf_load_prop (info, &prop_type, &prop_size))
        return FALSE;

      switch (prop_type)
        {
        case PROP_END:
          return TRUE;

        case PROP_ACTIVE_CHANNEL:
          info->selected_channels = g_list_prepend (info->selected_channels, *channel);
          break;

        case PROP_SELECTION:
          {
            GimpChannel *mask;
            GList       *iter;

            if (is_mask)
              {
                /* PROP_SELECTION is not valid for masks, and we have to avoid
                 * overwriting the channel.
                 */
                continue;
              }

            if (*channel == gimp_image_get_mask (image))
              {
                /* PROP_SELECTION was already seen once for this
                 * channel. Let's silently ignore the second identical
                 * property to avoid a double free.
                 */
                continue;
              }
            else if (gimp_image_get_mask (image) != NULL &&
                     ! gimp_channel_is_empty (gimp_image_get_mask (image)))
              {
                /* This would happen when PROP_SELECTION was already set
                 * on a previous channel. This is a minor case of data
                 * loss (we don't know which selection was the right one
                 * and we drop the non-first ones), and also means it's
                 * a broken XCF, though it's not a major bug either. So
                 * let's go with a stderr print.
                 */
                g_printerr ("PROP_SELECTION property was set on 2 channels (skipping)\n");
                continue;
              }

            /* We're going to delete *channel, Don't leave its pointer
             * in @info. See bug #767873.
             */
            for (iter = info->selected_channels; iter; iter = iter->next)
              if (*channel == iter->data)
                {
                  info->selected_channels = g_list_delete_link (info->selected_channels, iter);
                  break;
                }

            mask =
              gimp_selection_new (image,
                                  gimp_item_get_width  (GIMP_ITEM (*channel)),
                                  gimp_item_get_height (GIMP_ITEM (*channel)));
            gimp_image_take_mask (image, mask);

            gimp_drawable_steal_buffer (GIMP_DRAWABLE (mask),
                                        GIMP_DRAWABLE (*channel));
            g_object_unref (*channel);
            *channel = mask;

            /* Don't restore info->selected_channels because the
             * selection can't be the active channel
             */
          }
          break;

        case PROP_OPACITY:
          {
            guint32 opacity;

            xcf_read_int32 (info, &opacity, 1);

            gimp_channel_set_opacity (*channel, opacity / 255.0, FALSE);
          }
          break;

        case PROP_FLOAT_OPACITY:
          {
            gfloat opacity;

            xcf_read_float (info, &opacity, 1);

            gimp_channel_set_opacity (*channel, opacity, FALSE);
          }
          break;

        case PROP_VISIBLE:
          {
            gboolean visible;

            xcf_read_int32 (info, (guint32 *) &visible, 1);

            gimp_item_set_visible (GIMP_ITEM (*channel), visible, FALSE);
          }
          break;

        case PROP_COLOR_TAG:
          {
            GimpColorTag color_tag;

            xcf_read_int32 (info, (guint32 *) &color_tag, 1);

            gimp_item_set_color_tag (GIMP_ITEM (*channel), color_tag, FALSE);
          }
          break;

        case PROP_LINKED:
          {
            gboolean linked;

            xcf_read_int32 (info, (guint32 *) &linked, 1);

            if (linked)
              info->linked_channels = g_list_prepend (info->linked_channels, *channel);
          }
          break;

        case PROP_LOCK_CONTENT:
          {
            gboolean lock_content;

            xcf_read_int32 (info, (guint32 *) &lock_content, 1);

            if (gimp_item_can_lock_content (GIMP_ITEM (*channel)))
              gimp_item_set_lock_content (GIMP_ITEM (*channel),
                                          lock_content, FALSE);
          }
          break;

        case PROP_LOCK_POSITION:
          {
            gboolean lock_position;

            xcf_read_int32 (info, (guint32 *) &lock_position, 1);

            if (gimp_item_can_lock_position (GIMP_ITEM (*channel)))
              gimp_item_set_lock_position (GIMP_ITEM (*channel),
                                           lock_position, FALSE);
          }
          break;

        case PROP_LOCK_VISIBILITY:
          {
            gboolean lock_visibility;

            xcf_read_int32 (info, (guint32 *) &lock_visibility, 1);

            if (gimp_item_can_lock_visibility (GIMP_ITEM (*channel)))
              gimp_item_set_lock_visibility (GIMP_ITEM (*channel),
                                             lock_visibility, FALSE);
          }
          break;

        case PROP_SHOW_MASKED:
          {
            gboolean show_masked;

            xcf_read_int32 (info, (guint32 *) &show_masked, 1);

            gimp_channel_set_show_masked (*channel, show_masked);
          }
          break;

        case PROP_COLOR:
          {
            guchar  col[3];
            gdouble opacity;

            /* Load existing opacity */
            opacity = gimp_channel_get_opacity (*channel);

            xcf_read_int8 (info, (guint8 *) col, 3);

            gegl_color_set_pixel ((*channel)->color, babl_format ("R'G'B' u8"), col);

            gimp_channel_set_opacity (*channel, opacity, FALSE);
          }
          break;

        case PROP_FLOAT_COLOR:
          {
            gfloat  col[3];
            gdouble opacity;

            /* Load existing opacity */
            opacity = gimp_channel_get_opacity (*channel);

            xcf_read_float (info, col, 3);

            /* TODO: is the channel color in sRGB or in the image's color space? */
            gegl_color_set_pixel ((*channel)->color, babl_format ("R'G'B' float"), col);

            gimp_channel_set_opacity (*channel, opacity, FALSE);
          }
          break;

        case PROP_TATTOO:
          {
            GimpTattoo tattoo;

            xcf_read_int32 (info, (guint32 *) &tattoo, 1);

            gimp_item_set_tattoo (GIMP_ITEM (*channel), tattoo);
          }
          break;

        case PROP_PARASITES:
          {
            goffset base = info->cp;

            while ((info->cp - base) < prop_size)
              {
                GimpParasite *p     = xcf_load_parasite (info);
                GError       *error = NULL;

                if (! p)
                  return FALSE;

                if (! gimp_item_parasite_validate (GIMP_ITEM (*channel), p,
                                                    &error))
                  {
                    gimp_message (info->gimp, G_OBJECT (info->progress),
                                  GIMP_MESSAGE_WARNING,
                                  "Warning, invalid channel parasite in XCF file: %s",
                                  error->message);
                    g_clear_error (&error);
                  }
                else
                  {
                    gimp_item_parasite_attach (GIMP_ITEM (*channel), p, FALSE);
                  }

                gimp_parasite_free (p);
              }

            if (info->cp - base != prop_size)
              gimp_message_literal (info->gimp, G_OBJECT (info->progress),
                                    GIMP_MESSAGE_WARNING,
                                    "Error while loading a channel's parasites");
          }
          break;

        case PROP_ITEM_SET_ITEM:
            {
              GimpItemList *set;
              guint32       n;

              xcf_read_int32 (info, &n, 1);
              set = g_list_nth_data (info->channel_sets, n);
              if (set == NULL)
                g_printerr ("xcf: unknown channel set: %d (skipping)\n", n);
              else if (! g_type_is_a (G_TYPE_FROM_INSTANCE (*channel),
                                      gimp_item_list_get_item_type (set)))
                g_printerr ("xcf: channel '%s' cannot be added to item set '%s' with item type %s (skipping)\n",
                            gimp_object_get_name (*channel), gimp_object_get_name (set),
                            g_type_name (gimp_item_list_get_item_type (set)));
              else
                gimp_item_list_add (set, GIMP_ITEM (*channel));
            }
          break;

        default:
#ifdef GIMP_UNSTABLE
          g_printerr ("unexpected/unknown channel property: %d (skipping)\n",
                      prop_type);
#endif
          if (! xcf_skip_unknown_prop (info, prop_size))
            return FALSE;
          break;
        }
    }

  return FALSE;
}

static gboolean
xcf_load_effect_props (XcfInfo      *info,
                       FilterData   *filter)
{
  PropType prop_type;
  guint32  prop_size;

  while (TRUE)
    {
      goffset next_prop;

      if (! xcf_load_prop (info, &prop_type, &prop_size))
        return FALSE;

      next_prop = info->cp + prop_size;
      switch (prop_type)
        {
        case PROP_END:
          return TRUE;

        case PROP_VISIBLE:
          {
            gboolean visible;

            xcf_read_int32 (info, (guint32 *) &visible, 1);
            filter->is_visible = visible;
          }
          break;

        case PROP_FLOAT_OPACITY:
          {
            gfloat opacity;

            xcf_read_float (info, &opacity, 1);

            filter->opacity = (gdouble) opacity;
          }
          break;

        case PROP_MODE:
          {
            GimpLayerMode mode;

            xcf_read_int32 (info, (guint32 *) &mode, 1);

            if (mode == GIMP_LAYER_MODE_OVERLAY_LEGACY)
              mode = GIMP_LAYER_MODE_SOFTLIGHT_LEGACY;

            filter->paint_mode = mode;
          }
          break;

        case PROP_BLEND_SPACE:
          {
            gint32 blend_space;

            xcf_read_int32 (info, (guint32 *) &blend_space, 1);

            /* TODO: Revisit when blend space can be set
             * on filter effects */
            blend_space = GIMP_LAYER_COLOR_SPACE_AUTO;

            filter->blend_space = blend_space;
          }
          break;

        case PROP_COMPOSITE_SPACE:
          {
            gint32 composite_space;

            xcf_read_int32 (info, (guint32 *) &composite_space, 1);

            /* TODO: Revisit when composite space can be set
             * on filter effects */
            composite_space = GIMP_LAYER_COLOR_SPACE_AUTO;

            filter->composite_space = composite_space;
          }
          break;

        case PROP_COMPOSITE_MODE:
          {
            gint32 composite_mode;

            xcf_read_int32 (info, (guint32 *) &composite_mode, 1);

            /* TODO: Revisit when composite mode can be set
             * on filter effects */
            composite_mode = GIMP_LAYER_COMPOSITE_AUTO;

            filter->composite_mode = composite_mode;
          }
          break;

        case PROP_FILTER_REGION:
          {
            GimpFilterRegion region;

            xcf_read_int32 (info, (guint32 *) &region, 1);

            filter->region = region;
          }
          break;

        case PROP_FILTER_ARGUMENT:
          {
            GParamSpec *pspec;
            gchar      *filter_prop_name;
            guint32     filter_type       = FILTER_PROP_UNKNOWN;
            GValue      filter_prop_value = G_VALUE_INIT;
            gboolean    valid_prop_value  = TRUE;

            xcf_read_string (info, &filter_prop_name, 1);
            xcf_read_int32 (info, (guint32 *) &filter_type, 1);

            /* Check if valid property first */
            if (! (pspec = gegl_operation_find_property (filter->operation_name,
                                                         filter_prop_name)))
              {
                gimp_message (info->gimp, G_OBJECT (info->progress),
                              GIMP_MESSAGE_WARNING,
                              _("XCF Warning: version %s of filter \"%s\" does not "
                                "have the %s property. The property was ignored and "
                                "the filter may not render properly.\n"
                                "This should not happen. "
                                "You should report the issue to the filter's developers."),
                              filter->op_version ? filter->op_version : "0:0",
                              filter->operation_name, filter_prop_name);
                valid_prop_value = FALSE;
                goto set_or_seek_node_property;
              }

            switch (filter_type)
              {
                case FILTER_PROP_INT:
                case FILTER_PROP_UINT:
                case FILTER_PROP_ENUM:
                  {
                    guint32 value;

                    xcf_read_int32 (info, (guint32 *) &value, 1);
                    g_value_init (&filter_prop_value, pspec->value_type);
                    if (filter_type == FILTER_PROP_INT)
                      g_value_set_int (&filter_prop_value, value);
                    else if (filter_type == FILTER_PROP_UINT)
                      g_value_set_uint (&filter_prop_value, value);
                    else
                      g_value_set_enum (&filter_prop_value, value);
                  }
                  break;

                case FILTER_PROP_BOOL:
                  {
                    gboolean value;

                    xcf_read_int32 (info, (guint32 *) &value, 1);
                    g_value_init (&filter_prop_value, G_TYPE_BOOLEAN);
                    g_value_set_boolean (&filter_prop_value, value);
                  }
                  break;

                case FILTER_PROP_FLOAT:
                  {
                    gfloat value;

                    xcf_read_float (info, &value, 1);
                    g_value_init (&filter_prop_value, G_TYPE_FLOAT);
                    g_value_set_float (&filter_prop_value, value);
                  }
                  break;

                case FILTER_PROP_STRING:
                  {
                    gchar *value;

                    xcf_read_string (info, &value, 1);
                    g_value_init (&filter_prop_value, G_TYPE_STRING);
                    g_value_set_string (&filter_prop_value, value);

                    g_free (value);
                  }
                  break;

                case FILTER_PROP_CONFIG:
                  {
                    GObject *config;
                    gchar   *serialized;
                    GError  *error = NULL;

                    if (! g_type_is_a (pspec->value_type, GIMP_TYPE_CONFIG))
                      {
                        gimp_message (info->gimp, G_OBJECT (info->progress),
                                      GIMP_MESSAGE_WARNING,
                                      "XCF Warning: property '%s' of filter '%s' is a %s, which does not implement GimpConfig interface.\n",
                                      filter_prop_name, filter->operation_name, g_type_name (pspec->value_type));
                        valid_prop_value = FALSE;
                        break;
                      }

                    xcf_read_string (info, &serialized, 1);
                    if (! (serialized && strlen (serialized)))
                      {
                        gimp_message (info->gimp, G_OBJECT (info->progress),
                                      GIMP_MESSAGE_WARNING,
                                      "XCF Warning: failure to deserialize config object for property '%s' of filter '%s'.\n"
                                      "Serialized config could not be read.",
                                      filter_prop_name, filter->operation_name);
                        valid_prop_value = FALSE;
                        break;
                      }

                    g_value_init (&filter_prop_value, pspec->value_type);
                    config = g_object_new (pspec->value_type, NULL);
                    if (gimp_config_deserialize_string (GIMP_CONFIG (config), serialized, -1, NULL, &error))
                      {
                        g_value_set_object (&filter_prop_value, config);
                      }
                    else
                      {
                        gimp_message (info->gimp, G_OBJECT (info->progress),
                                      GIMP_MESSAGE_WARNING,
                                      "XCF Warning: failure to deserialize config object for property '%s' of filter '%s': %s\n"
                                      "Serialized config was: %s",
                                      filter_prop_name, filter->operation_name,
                                      error->message, serialized);
                        valid_prop_value = FALSE;
                      }

                    g_object_unref (config);
                    g_free (serialized);
                  }
                  break;

                case FILTER_PROP_COLOR:
                  {
                    GeglColor  *color = gegl_color_new (NULL);
                    const Babl *format;
                    gchar      *encoding;
                    guint8     *data  = NULL;
                    gint        data_length;
                    gint        profile_data_length;

                    g_value_init (&filter_prop_value, GEGL_TYPE_COLOR);

                    xcf_read_string (info, &encoding, 1);
                    if (! babl_format_exists (encoding))
                      {
                        gimp_message (info->gimp, G_OBJECT (info->progress),
                                      GIMP_MESSAGE_WARNING,
                                      "XCF Warning: format \"%s\" for "
                                      "property '%s' of filter '%s' is "
                                      "invalid. The color was discarded.",
                                      encoding, filter->operation_name,
                                      filter_prop_name);

                        g_free (encoding);
                        valid_prop_value = FALSE;
                        break;
                      }

                    format = babl_format (encoding);
                    g_free (encoding);

                    xcf_read_int32 (info, (guint32 *) &data_length, 1);
                    if (data_length != babl_format_get_bytes_per_pixel (format))
                      {
                        gimp_message (info->gimp, G_OBJECT (info->progress),
                                      GIMP_MESSAGE_WARNING,
                                      "XCF Warning: format \"%s\" for "
                                      "property '%s' of filter '%s' expected "
                                      "%d bpp, but color was serialized as %d "
                                      "bpp. The color was discarded.",
                                      babl_get_name (format), filter->operation_name,
                                      filter_prop_name,
                                      babl_format_get_bytes_per_pixel (format),
                                      data_length);

                        valid_prop_value = FALSE;
                        break;
                      }

                    data = g_new (guint8, data_length);
                    xcf_read_int8 (info, data, data_length);

                    xcf_read_int32 (info, (guint32 *) &profile_data_length, 1);
                    if (profile_data_length > 0)
                      {
                        const Babl       *space = NULL;
                        GimpColorProfile *profile;
                        guint8           *profile_data;
                        GError           *error = NULL;

                        profile_data = g_new (guint8, profile_data_length);
                        xcf_read_int8 (info, profile_data, profile_data_length);
                        profile = gimp_color_profile_new_from_icc_profile (profile_data,
                                                                           profile_data_length,
                                                                           &error);

                        if (profile)
                          {
                            space = gimp_color_profile_get_space (profile,
                                                                  GIMP_COLOR_RENDERING_INTENT_RELATIVE_COLORIMETRIC,
                                                                  &error);

                            if (! space)
                              {
                                gimp_message (info->gimp, G_OBJECT (info->progress),
                                              GIMP_MESSAGE_WARNING,
                                              "XCF Warning: failed to create "
                                              "Babl space for serialized color "
                                              "from profile");

                                g_free (data);
                                valid_prop_value = FALSE;
                                break;
                              }
                            g_object_unref (profile);
                          }
                        else
                          {
                            gimp_message (info->gimp, G_OBJECT (info->progress),
                                          GIMP_MESSAGE_WARNING,
                                          "XCF Warning: Invalid profile for "
                                          "serialized color.");
                          }

                        format = babl_format_with_space (babl_format_get_encoding (format), space);
                        g_free (profile_data);
                      }

                      gegl_color_set_pixel (color, format, data);
                      g_value_set_object (&filter_prop_value, color);

                      g_free (data);
                      g_object_unref (color);
                  }
                  break;

                default:
                  gimp_message (info->gimp, G_OBJECT (info->progress),
                                GIMP_MESSAGE_WARNING,
                                "XCF Warning: property '%s' of filter '%s' holds unsupported type %s.\n",
                                filter_prop_name, filter->operation_name, g_type_name (pspec->value_type));
                  valid_prop_value = FALSE;
                  break;
              }

set_or_seek_node_property:
            if (valid_prop_value)
              gegl_node_set_property (filter->operation, filter_prop_name,
                                      &filter_prop_value);
            else
              xcf_seek_pos (info, next_prop, NULL);

            g_value_unset (&filter_prop_value);
            g_free (filter_prop_name);
          }
          break;

        case PROP_FILTER_CLIP:
          {
            gboolean clip;

            xcf_read_int32 (info, (guint32 *) &clip, 1);
            filter->clip = clip;
          }
          break;

        default:
#ifdef GIMP_UNSTABLE
          g_printerr ("unexpected/unknown effect property: %d (skipping)\n",
                      prop_type);
#endif
          if (! xcf_skip_unknown_prop (info, prop_size))
            return FALSE;
        }
    }

  return FALSE;
}

static gboolean
xcf_load_path_props (XcfInfo    *info,
                     GimpImage  *image,
                     GimpPath  **path)
{
  PropType prop_type;
  guint32  prop_size;

  while (TRUE)
    {
      if (! xcf_load_prop (info, &prop_type, &prop_size))
        return FALSE;

      switch (prop_type)
        {
        case PROP_END:
          return TRUE;

        case PROP_SELECTED_PATH:
          info->selected_paths = g_list_prepend (info->selected_paths, *path);
          break;

        case PROP_VISIBLE:
          {
            gboolean visible;

            xcf_read_int32 (info, (guint32 *) &visible, 1);

            gimp_item_set_visible (GIMP_ITEM (*path), visible, FALSE);
          }
          break;

        case PROP_COLOR_TAG:
          {
            GimpColorTag color_tag;

            xcf_read_int32 (info, (guint32 *) &color_tag, 1);

            gimp_item_set_color_tag (GIMP_ITEM (*path), color_tag, FALSE);
          }
          break;

        case PROP_LOCK_CONTENT:
          {
            gboolean lock_content;

            xcf_read_int32 (info, (guint32 *) &lock_content, 1);

            if (gimp_item_can_lock_content (GIMP_ITEM (*path)))
              gimp_item_set_lock_content (GIMP_ITEM (*path),
                                          lock_content, FALSE);
          }
          break;

        case PROP_LOCK_POSITION:
          {
            gboolean lock_position;

            xcf_read_int32 (info, (guint32 *) &lock_position, 1);

            if (gimp_item_can_lock_position (GIMP_ITEM (*path)))
              gimp_item_set_lock_position (GIMP_ITEM (*path),
                                           lock_position, FALSE);
          }
          break;

        case PROP_LOCK_VISIBILITY:
          {
            gboolean lock_visibility;

            xcf_read_int32 (info, (guint32 *) &lock_visibility, 1);

            if (gimp_item_can_lock_visibility (GIMP_ITEM (*path)))
              gimp_item_set_lock_visibility (GIMP_ITEM (*path),
                                             lock_visibility, FALSE);
          }
          break;

        case PROP_TATTOO:
          {
            GimpTattoo tattoo;

            xcf_read_int32 (info, (guint32 *) &tattoo, 1);

            gimp_item_set_tattoo (GIMP_ITEM (*path), tattoo);
          }
          break;

        case PROP_PARASITES:
          {
            goffset base = info->cp;

            while ((info->cp - base) < prop_size)
              {
                GimpParasite *p     = xcf_load_parasite (info);
                GError       *error = NULL;

                if (! p)
                  return FALSE;

                if (! gimp_item_parasite_validate (GIMP_ITEM (*path), p,
                                                    &error))
                  {
                    gimp_message (info->gimp, G_OBJECT (info->progress),
                                  GIMP_MESSAGE_WARNING,
                                  "Warning, invalid path parasite in XCF file: %s",
                                  error->message);
                    g_clear_error (&error);
                  }
                else
                  {
                    gimp_item_parasite_attach (GIMP_ITEM (*path), p, FALSE);
                  }

                gimp_parasite_free (p);
              }

            if (info->cp - base != prop_size)
              gimp_message_literal (info->gimp, G_OBJECT (info->progress),
                                    GIMP_MESSAGE_WARNING,
                                    "Error while loading a path's parasites");
          }
          break;

#if 0
        case PROP_ITEM_SET_ITEM:
            {
              GimpItemList *set;
              guint32       n;

              xcf_read_int32 (info, &n, 1);
              set = g_list_nth_data (info->path_sets, n);
              if (set == NULL)
                g_printerr ("xcf: unknown path set: %d (skipping)\n", n);
              else if (! g_type_is_a (G_TYPE_FROM_INSTANCE (*path),
                                      gimp_item_list_get_item_type (set)))
                g_printerr ("xcf: path '%s' cannot be added to item set '%s' with item type %s (skipping)\n",
                            gimp_object_get_name (*path), gimp_object_get_name (set),
                            g_type_name (gimp_item_list_get_item_type (set)));
              else
                gimp_item_list_add (set, GIMP_ITEM (*path));
            }
          break;
#endif

        default:
#ifdef GIMP_UNSTABLE
          g_printerr ("unexpected/unknown path property: %d (skipping)\n",
                      prop_type);
#endif
          if (! xcf_skip_unknown_prop (info, prop_size))
            return FALSE;
          break;
        }
    }

  return FALSE;
}

static gboolean
xcf_load_prop (XcfInfo  *info,
               PropType *prop_type,
               guint32  *prop_size)
{
  if (G_UNLIKELY (xcf_read_int32 (info, (guint32 *) prop_type, 1) != 4))
    return FALSE;

  if (G_UNLIKELY (xcf_read_int32 (info, (guint32 *) prop_size, 1) != 4))
    return FALSE;

  GIMP_LOG (XCF, "prop type=%d size=%u", *prop_type, *prop_size);

  return TRUE;
}

static GimpLayer *
xcf_load_layer (XcfInfo    *info,
                GimpImage  *image,
                GList     **item_path,
                gint       *n_broken_effects)
{
  GimpLayer         *layer;
  GimpLayerMask     *layer_mask;
  goffset            hierarchy_offset;
  goffset            effects_offset = 0;
  goffset            layer_mask_offset;
  gboolean           apply_mask = TRUE;
  gboolean           edit_mask  = FALSE;
  gboolean           show_mask  = FALSE;
  GList             *selected;
  GList             *linked;
  GList             *filter_data_list = NULL;
  gboolean           floating;
  guint32            group_layer_flags = 0;
  guint32            text_layer_flags = 0;
  gint               filter_count = 0;
  gint               width;
  gint               height;
  gint               type;
  GimpImageBaseType  base_type;
  gboolean           has_alpha;
  const Babl        *format;
  gboolean           is_fs_drawable;
  gchar             *name;
  goffset            cur_offset;

  /* check and see if this is the drawable the floating selection
   *  is attached to. if it is then we'll do the attachment in our caller.
   */
  is_fs_drawable = (info->cp == info->floating_sel_offset);

  /* read in the layer width, height, type and name */
  xcf_read_int32  (info, (guint32 *) &width,  1);
  xcf_read_int32  (info, (guint32 *) &height, 1);
  xcf_read_int32  (info, (guint32 *) &type,   1);
  xcf_read_string (info,             &name,   1);

  GIMP_LOG (XCF, "width=%d, height=%d, type=%d, name='%s'",
            width, height, type, name);

  switch (type)
    {
    case GIMP_RGB_IMAGE:
      base_type = GIMP_RGB;
      has_alpha = FALSE;
      break;

    case GIMP_RGBA_IMAGE:
      base_type = GIMP_RGB;
      has_alpha = TRUE;
      break;

    case GIMP_GRAY_IMAGE:
      base_type = GIMP_GRAY;
      has_alpha = FALSE;
      break;

    case GIMP_GRAYA_IMAGE:
      base_type = GIMP_GRAY;
      has_alpha = TRUE;
      break;

    case GIMP_INDEXED_IMAGE:
      base_type = GIMP_INDEXED;
      has_alpha = FALSE;
      break;

    case GIMP_INDEXEDA_IMAGE:
      base_type = GIMP_INDEXED;
      has_alpha = TRUE;
      break;

    default:
      g_free (name);
      return NULL;
    }

  if (width <= 0 || height <= 0 ||
      width > GIMP_MAX_IMAGE_SIZE || height > GIMP_MAX_IMAGE_SIZE)
    {
      gboolean is_group_layer = FALSE;
      gboolean is_text_layer  = FALSE;
      goffset  saved_pos;

      saved_pos = info->cp;
      /* Load item path and check if this is a group or text layer. */
      xcf_check_layer_props (info, item_path, &is_group_layer, &is_text_layer);
      if ((is_text_layer || is_group_layer) &&
          xcf_seek_pos (info, saved_pos, NULL))
        {
          /* Something is wrong, but leave a chance to the layer because
           * anyway group and text layer depends on their contents.
           */
          width = height = 1;
          g_clear_pointer (item_path, g_list_free);
        }
      else
        {
          g_free (name);
          return NULL;
        }
    }

  if (base_type == GIMP_GRAY)
    {
      /* do not use gimp_image_get_layer_format() because it might
       * be the floating selection of a channel or mask
       */
      format = gimp_image_get_format (image, base_type,
                                      gimp_image_get_precision (image),
                                      has_alpha,
                                      NULL /* we will fix the space later */);
    }
  else
    {
      format = gimp_image_get_layer_format (image, has_alpha);
    }

  /* create a new layer */
  layer = gimp_layer_new (image, width, height,
                          format, name,
                          GIMP_OPACITY_OPAQUE, GIMP_LAYER_MODE_NORMAL);
  g_free (name);
  if (! layer)
    return NULL;

  /* read in the layer properties */
  if (! xcf_load_layer_props (info, image, &layer, item_path,
                              &apply_mask, &edit_mask, &show_mask,
                              &text_layer_flags, &group_layer_flags))
    goto error;

  GIMP_LOG (XCF, "layer props loaded");

  xcf_progress_update (info);

  /* call the evil text layer hack that might change our layer pointer */
  selected = g_list_find (info->selected_layers, layer);
  linked   = g_list_find (info->linked_layers, layer);
  floating = (info->floating_sel == layer);

  if (gimp_text_layer_xcf_load_hack (&layer))
    {
      gimp_text_layer_set_xcf_flags (GIMP_TEXT_LAYER (layer),
                                     text_layer_flags);

      if (selected)
        {
          info->selected_layers = g_list_delete_link (info->selected_layers, selected);
          info->selected_layers = g_list_prepend (info->selected_layers, layer);
        }
      if (linked)
        {
          info->linked_layers = g_list_delete_link (info->linked_layers, linked);
          info->linked_layers = g_list_prepend (info->linked_layers, layer);
        }
      if (floating)
        info->floating_sel = layer;
    }

  /* if this is not the floating selection, we can fix the layer's
   * space already now, the function will do nothing if we already
   * created the layer with the right format
   */
  if (! floating && base_type == GIMP_GRAY)
    gimp_layer_fix_format_space (layer, FALSE, FALSE);

  /* read the hierarchy and layer mask offsets */
  cur_offset = info->cp;
  xcf_read_offset (info, &hierarchy_offset,  1);
  xcf_read_offset (info, &layer_mask_offset, 1);
  if (info->file_version >= 20)
    xcf_read_offset (info, &effects_offset,  1);

  /* read in the hierarchy (ignore it for group layers, both as an
   * optimization and because the hierarchy's extents don't match
   * the group layer's tiles)
   */
  if (! gimp_viewable_get_children (GIMP_VIEWABLE (layer)))
    {
      if (hierarchy_offset < cur_offset)
        {
          GIMP_LOG (XCF, "Invalid layer hierarchy offset!");
          goto error;
        }
      if (! xcf_seek_pos (info, hierarchy_offset, NULL))
        goto error;

      GIMP_LOG (XCF, "loading buffer");

      if (! xcf_load_buffer (info,
                             gimp_drawable_get_buffer (GIMP_DRAWABLE (layer))))
        goto error;

      GIMP_LOG (XCF, "buffer loaded");

      xcf_progress_update (info);
    }
  else
    {
      gboolean expanded = group_layer_flags & XCF_GROUP_ITEM_EXPANDED;

      gimp_viewable_set_expanded (GIMP_VIEWABLE (layer), expanded);
    }

  cur_offset += info->bytes_per_offset;

  /* read in the layer mask */
  if (layer_mask_offset != 0)
    {
      if (layer_mask_offset < cur_offset)
        {
          GIMP_LOG (XCF, "Invalid layer mask offset!");
          goto error;
        }
      if (! xcf_seek_pos (info, layer_mask_offset, NULL))
        goto error;

      layer_mask = xcf_load_layer_mask (info, image);
      if (! layer_mask)
        goto error;

      cur_offset += info->bytes_per_offset;

      xcf_progress_update (info);

      /* don't add the layer mask yet, that won't work for group
       * layers which update their size automatically; instead
       * attach it so it can be added when all layers are loaded
       */
      g_object_set_data_full (G_OBJECT (layer), "gimp-layer-mask",
                              g_object_ref_sink (layer_mask),
                              (GDestroyNotify) g_object_unref);
      g_object_set_data (G_OBJECT (layer), "gimp-layer-mask-apply",
                         GINT_TO_POINTER (apply_mask));
      g_object_set_data (G_OBJECT (layer), "gimp-layer-mask-edit",
                         GINT_TO_POINTER (edit_mask));
      g_object_set_data (G_OBJECT (layer), "gimp-layer-mask-show",
                         GINT_TO_POINTER (show_mask));
    }
  else
    {
      /* If no layer mask, move ahead to layer effects */
      cur_offset += info->bytes_per_offset;
    }

  /* read in any layer effects and effect masks */
  while (effects_offset != 0)
    {
      FilterData *filter_data;

      if (effects_offset < cur_offset)
        {
          GIMP_LOG (XCF, "Invalid effect offset: %" G_GOFFSET_FORMAT
                    " at offset: %" G_GOFFSET_FORMAT, effects_offset, cur_offset);
          goto error;
        }

      /* seek to the effect offset */
      if (! xcf_seek_pos (info, effects_offset, NULL))
        goto error;

      filter_data = xcf_load_effect (info, image, GIMP_DRAWABLE (layer));
      if (! filter_data)
        {
          (*n_broken_effects)++;
        }
      else
        {
          filter_data_list = g_list_prepend (filter_data_list, filter_data);
          filter_count++;
        }
      xcf_progress_update (info);

      /* restore the saved position so we'll be ready to
       * read the next offset.
       */
      cur_offset += info->bytes_per_offset;
      if (! xcf_seek_pos (info, cur_offset, NULL))
        goto error;

      /* read in the offset of the next effect */
      if (xcf_read_offset (info, &effects_offset, 1) < info->bytes_per_offset)
        {
          GIMP_LOG (XCF, "Failed to read effects offset"
                    " at offset: %" G_GOFFSET_FORMAT, info->cp);
          break;
        }
    }

  if (filter_count > 0)
    g_object_set_data_full (G_OBJECT (layer), "gimp-layer-effects", filter_data_list,
                            (GDestroyNotify) xcf_load_free_effects);

  /* attach the floating selection... */
  if (is_fs_drawable)
    info->floating_sel_drawable = GIMP_DRAWABLE (layer);

  return layer;

 error:
  info->selected_layers = g_list_remove (info->selected_layers, layer);

  if (info->floating_sel == layer)
    info->floating_sel = NULL;

  if (info->floating_sel_drawable == GIMP_DRAWABLE (layer))
    info->floating_sel_drawable = NULL;

  g_object_unref (layer);

  return NULL;
}

static GimpChannel *
xcf_load_channel (XcfInfo   *info,
                  GimpImage *image)
{
  GimpChannel *channel;
  goffset      hierarchy_offset;
  gint         width;
  gint         height;
  gboolean     is_fs_drawable;
  gchar       *name;
  goffset      cur_offset;

  /* check and see if this is the drawable the floating selection
   *  is attached to. if it is then we'll do the attachment in our caller.
   */
  is_fs_drawable = (info->cp == info->floating_sel_offset);

  /* read in the layer width, height and name */
  xcf_read_int32 (info, (guint32 *) &width,  1);
  xcf_read_int32 (info, (guint32 *) &height, 1);
  if (width <= 0 || height <= 0 ||
      width > GIMP_MAX_IMAGE_SIZE || height > GIMP_MAX_IMAGE_SIZE)
    {
      GIMP_LOG (XCF, "Invalid channel size %d x %d.", width, height);
      return NULL;
    }

  xcf_read_string (info, &name, 1);
  GIMP_LOG (XCF, "Channel width=%d, height=%d, name='%s'",
            width, height, name);

  /* create a new channel */
  channel = gimp_channel_new (image, width, height, name, NULL);
  g_free (name);
  if (!channel)
    return NULL;

  /* read in the channel properties */
  if (! xcf_load_channel_props (info, image, &channel, FALSE))
    goto error;

  xcf_progress_update (info);

  /* read the hierarchy offset */
  cur_offset = info->cp;
  xcf_read_offset (info, &hierarchy_offset, 1);

  if (hierarchy_offset < cur_offset)
    {
      GIMP_LOG (XCF, "Invalid hierarchy offset!");
      goto error;
    }

  /* read in the hierarchy */
  if (! xcf_seek_pos (info, hierarchy_offset, NULL))
    goto error;

  if (! xcf_load_buffer (info,
                         gimp_drawable_get_buffer (GIMP_DRAWABLE (channel))))
    goto error;

  xcf_progress_update (info);

  if (is_fs_drawable)
    info->floating_sel_drawable = GIMP_DRAWABLE (channel);

  return channel;

 error:
  /* don't unref the selection of a partially loaded XCF */
  if (channel != gimp_image_get_mask (image))
    {
      GList *iter;

      for (iter = info->selected_channels; iter; iter = iter->next)
        if (channel == iter->data)
          {
            info->selected_channels = g_list_delete_link (info->selected_channels, iter);
            break;
          }

      if (info->floating_sel_drawable == GIMP_DRAWABLE (channel))
        info->floating_sel_drawable = NULL;

      g_object_unref (channel);
    }

  return NULL;
}

/* Comes from gegl/operation/gegl-operations.h which is not public. */
GType gegl_operation_gtype_from_name (const gchar *name);

static FilterData *
xcf_load_effect (XcfInfo      *info,
                 GimpImage    *image,
                 GimpDrawable *drawable)
{
  FilterData  *filter;
  GimpChannel *effect_mask = NULL;
  goffset      mask_offset = 0;
  gchar       *string;
  GType        op_type;

  filter = g_new0 (FilterData, 1);

  /* Effect name */
  xcf_read_string (info, &string, 1);
  filter->name = string;

  /* Effect icon */
  xcf_read_string (info, &string, 1);
  filter->icon_name = string;

  /* Effect operation */
  xcf_read_string (info, &string, 1);
  filter->operation_name = string;

  op_type = gegl_operation_gtype_from_name (filter->operation_name);
  if (g_type_is_a (op_type, GEGL_TYPE_OPERATION_SINK))
    {
      /* Forbid filters which directly write into files. These should
       * not be creatable through GIMP UI, but just in case someone
       * builds one such XCF file (maliciously or by mistake/through a
       * bug), let's prevent this filter to overwrite random files on
       * load.
       */
      filter->unsupported_operation = TRUE;

      gimp_message (info->gimp, G_OBJECT (info->progress),
                    GIMP_MESSAGE_WARNING,
                    /* TODO: localize after string freeze. */
                    "XCF Warning: the \"%s\" (%s) filter is "
                    "unsafe. It was discarded.",
                    filter->name, filter->operation_name);

      return filter;
    }
  else if (g_strcmp0 (filter->operation_name, "gegl:gegl") == 0 &&
           g_getenv ("GIMP_ALLOW_GEGL_GRAPH_LAYER_EFFECT") == NULL)
    {
      filter->unsupported_operation = TRUE;

      gimp_message (info->gimp, G_OBJECT (info->progress),
                    GIMP_MESSAGE_WARNING,
                    /* TODO: localize after string freeze. */
                    "XCF Warning: the \"%s\" (%s) filter is unsafe. It was discarded.\n"
                    "For development purpose, set environment variable GIMP_ALLOW_GEGL_GRAPH_LAYER_EFFECT.",
                    filter->name, filter->operation_name);

      return filter;
    }

  if (info->file_version >= 22)
    {
      xcf_read_string (info, &string, 1);
      filter->op_version = string;
    }

  if (! gegl_has_operation (filter->operation_name) ||
      ! g_strcmp0 (filter->operation_name, "gegl:nop"))
    {
      filter->unsupported_operation = TRUE;

      if (! g_strcmp0 (filter->operation_name, "gegl:nop"))
        gimp_message (info->gimp, G_OBJECT (info->progress),
                      GIMP_MESSAGE_WARNING,
                      "XCF Warning: A filter was saved as a "
                      "gegl:nop. This should not happen. Please "
                      "report this to the developers.");
      else
        gimp_message (info->gimp, G_OBJECT (info->progress),
                      GIMP_MESSAGE_WARNING,
                      "XCF Warning: the \"%s\" (%s) filter is "
                      "not installed. It was discarded.",
                      filter->name, filter->operation_name);

      return filter;
    }

  if (filter->op_version &&
      g_strcmp0 (gegl_operation_get_op_version (filter->operation_name),
                 filter->op_version) != 0)
    {
      filter->unsupported_operation = TRUE;

      gimp_message (info->gimp, G_OBJECT (info->progress),
                    GIMP_MESSAGE_WARNING,
                    _("XCF Warning: version '%s' of filter '%s' is unsupported. "
                      "The filter was discarded.\n"
                      "It either means that you are using an old version of the filter or "
                      "that it was updated without proper version management. "
                      "In the latter case, you should report the issue to the filter's developers."),
                    filter->op_version, filter->operation_name);

      return filter;
    }

  /* We'll set valid properties individually */
  filter->operation = gegl_node_new ();
  gegl_node_set (filter->operation,
                 "operation", filter->operation_name,
                 NULL);

  /* read in the effect properties */
  if (! xcf_load_effect_props (info, &(*filter)))
    goto error;

  /* seek to the effect mask offset */
  xcf_read_offset (info, &mask_offset, 1);

  if (mask_offset != 0)
    {
      if (! xcf_seek_pos (info, mask_offset, NULL))
        goto error;

      effect_mask = xcf_load_channel (info, image);
      if (! effect_mask)
        goto error;
    }

  filter->mask = effect_mask;

  xcf_progress_update (info);

  return filter;

  error:

  xcf_load_free_effect (filter);

  return NULL;
}

static void
xcf_load_free_effect (FilterData *data)
{
  g_free (data->name);
  g_free (data->icon_name);
  g_free (data->operation_name);
  g_free (data->op_version);
  g_clear_object (&data->operation);
  g_clear_object (&data->mask);
  g_free (data);
}

static void
xcf_load_free_effects (GList *effects)
{
  g_list_free_full (effects, (GDestroyNotify) xcf_load_free_effect);
}

/* The new path structure since XCF 18. */
static GimpPath *
xcf_load_path (XcfInfo   *info,
               GimpImage *image)
{
  GimpPath *path = NULL;
  gchar    *name;
  guint32   version;
  guint32   plength;
  guint32   num_strokes;
  goffset   base;
  gint      i;

  /* read in the path name. */
  xcf_read_string (info, &name,   1);

  GIMP_LOG (XCF, "Path name='%s'", name);

  /* create a new path */
  path = gimp_path_new (image, name);
  g_free (name);
  if (! path)
    return NULL;

  /* Read the path's payload size. */
  xcf_read_int32 (info, (guint32 *) &plength, 1);
  base = info->cp;

  /* read in the path properties */
  if (! xcf_load_path_props (info, image, &path))
    goto error;

  GIMP_LOG (XCF, "path props loaded");

  xcf_progress_update (info);

  xcf_read_int32 (info, (guint32 *) &version, 1);

  if (version != 1)
    {
      gimp_message (info->gimp, G_OBJECT (info->progress),
                    GIMP_MESSAGE_WARNING,
                    "Unknown path version: %d (skipping)", version);
      goto error;
    }

  /* Read the number of strokes. */
  xcf_read_int32  (info, &num_strokes, 1);

  for (i = 0; i < num_strokes; i++)
    {
      guint32      stroke_type_id;
      guint32      closed;
      guint32      num_axes;
      guint32      num_control_points;
      guint32      type;
      gfloat       coords[13] = GIMP_COORDS_DEFAULT_VALUES;
      GimpStroke  *stroke;
      gint         j;

      GimpValueArray *control_points;
      GValue          value  = G_VALUE_INIT;
      GimpAnchor      anchor = { { 0, } };
      GType           stroke_type;

      g_value_init (&value, GIMP_TYPE_ANCHOR);

      xcf_read_int32 (info, &stroke_type_id,     1);
      xcf_read_int32 (info, &closed,             1);
      xcf_read_int32 (info, &num_axes,           1);
      xcf_read_int32 (info, &num_control_points, 1);

#ifdef GIMP_XCF_PATH_DEBUG
      g_printerr ("stroke_type: %d, closed: %d, num_axes %d, len %d\n",
                  stroke_type_id, closed, num_axes, num_control_points);
#endif

      switch (stroke_type_id)
        {
        case XCF_STROKETYPE_BEZIER_STROKE:
          stroke_type = GIMP_TYPE_BEZIER_STROKE;
          break;

        default:
          g_printerr ("skipping unknown stroke type\n");
          xcf_seek_pos (info,
                        info->cp + 4 * num_axes * num_control_points,
                        NULL);
          continue;
        }

      if (num_axes < 2 || num_axes > 6)
        {
          g_printerr ("bad number of axes in stroke description\n");
          goto error;
        }

      control_points = gimp_value_array_new (num_control_points);

      anchor.selected = FALSE;

      for (j = 0; j < num_control_points; j++)
        {
          xcf_read_int32 (info, &type,  1);
          xcf_read_float (info, coords, num_axes);

          anchor.type              = type;
          anchor.position.x        = coords[0];
          anchor.position.y        = coords[1];
          anchor.position.pressure = coords[2];
          anchor.position.xtilt    = coords[3];
          anchor.position.ytilt    = coords[4];
          anchor.position.wheel    = coords[5];

          g_value_set_boxed (&value, &anchor);
          gimp_value_array_append (control_points, &value);

#ifdef GIMP_XCF_PATH_DEBUG
          g_printerr ("Anchor: %d, (%f, %f, %f, %f, %f, %f)\n", type,
                      coords[0], coords[1], coords[2], coords[3],
                      coords[4], coords[5]);
#endif
        }

      g_value_unset (&value);

      stroke = g_object_new (stroke_type,
                             "closed",         closed,
                             "control-points", control_points,
                             NULL);

      gimp_path_stroke_add (path, stroke);

      g_object_unref (stroke);
      gimp_value_array_unref (control_points);
    }


  if (plength != info->cp - base)
    {
      gimp_message (info->gimp, G_OBJECT (info->progress),
                    GIMP_MESSAGE_WARNING,
                    "Path payload size does not match stored size (skipping)");
      goto error;
    }

  return path;

error:

  xcf_seek_pos (info, base + plength, NULL);
  g_clear_object (&path);

  return NULL;
}

static GimpLayerMask *
xcf_load_layer_mask (XcfInfo   *info,
                     GimpImage *image)
{
  GimpLayerMask *layer_mask;
  GimpChannel   *channel;
  GList         *iter;
  goffset        hierarchy_offset;
  gint           width;
  gint           height;
  gboolean       is_fs_drawable;
  gchar         *name;
  GeglColor     *color = gegl_color_new ("black");
  goffset        cur_offset;

  /* check and see if this is the drawable the floating selection
   *  is attached to. if it is then we'll do the attachment in our caller.
   */
  is_fs_drawable = (info->cp == info->floating_sel_offset);

  /* read in the layer width, height and name */
  xcf_read_int32 (info, (guint32 *) &width,  1);
  xcf_read_int32 (info, (guint32 *) &height, 1);
  if (width <= 0 || height <= 0 ||
      width > GIMP_MAX_IMAGE_SIZE || height > GIMP_MAX_IMAGE_SIZE)
    {
      GIMP_LOG (XCF, "Invalid layer mask size %d x %d.", width, height);
      return NULL;
    }

  xcf_read_string (info, &name, 1);
  GIMP_LOG (XCF, "Layer mask width=%d, height=%d, name='%s'",
            width, height, name);

  /* create a new layer mask */
  layer_mask = gimp_layer_mask_new (image, width, height, name, color);
  g_object_unref (color);
  g_free (name);
  if (! layer_mask)
    return NULL;

  /* read in the layer_mask properties */
  channel = GIMP_CHANNEL (layer_mask);
  if (! xcf_load_channel_props (info, image, &channel, TRUE))
    goto error;

  xcf_progress_update (info);

  /* read the hierarchy offset */
  cur_offset = info->cp;
  xcf_read_offset (info, &hierarchy_offset, 1);

  if (hierarchy_offset < cur_offset)
    {
      GIMP_LOG (XCF, "Invalid hierarchy offset!");
      goto error;
    }

  /* read in the hierarchy */
  if (! xcf_seek_pos (info, hierarchy_offset, NULL))
    goto error;

  if (! xcf_load_buffer (info,
                         gimp_drawable_get_buffer (GIMP_DRAWABLE (layer_mask))))
    goto error;

  xcf_progress_update (info);

  /* attach the floating selection... */
  if (is_fs_drawable)
    info->floating_sel_drawable = GIMP_DRAWABLE (layer_mask);

  return layer_mask;

 error:
  for (iter = info->selected_channels; iter; iter = iter->next)
    if (layer_mask == iter->data)
      {
        info->selected_channels = g_list_delete_link (info->selected_channels, iter);
        break;
      }

  if (info->floating_sel_drawable == GIMP_DRAWABLE (layer_mask))
    info->floating_sel_drawable = NULL;

  g_object_unref (layer_mask);

  return NULL;
}

static gboolean
xcf_load_buffer (XcfInfo    *info,
                 GeglBuffer *buffer)
{
  const Babl *format;
  goffset     offset;
  gint        width;
  gint        height;
  gint        bpp;
  goffset     cur_offset;

  format = gegl_buffer_get_format (buffer);

  xcf_read_int32 (info, (guint32 *) &width,  1);
  xcf_read_int32 (info, (guint32 *) &height, 1);
  xcf_read_int32 (info, (guint32 *) &bpp,    1);

  /* make sure the values in the file correspond to the values
   *  calculated when the GeglBuffer was created.
   */
  if (width  != gegl_buffer_get_width (buffer)  ||
      height != gegl_buffer_get_height (buffer) ||
      bpp    != babl_format_get_bytes_per_pixel (format))
    return FALSE;

  cur_offset = info->cp;
  xcf_read_offset (info, &offset, 1); /* top level */

  if (offset < cur_offset)
    {
      GIMP_LOG (XCF, "Invalid buffer offset!");
      return FALSE;
    }

  /* seek to the level offset */
  if (! xcf_seek_pos (info, offset, NULL))
    return FALSE;

  /* read in the level */
  if (! xcf_load_level (info, buffer))
    return FALSE;

  /* discard levels below first.
   */

  return TRUE;
}


static gboolean
xcf_load_level (XcfInfo    *info,
                GeglBuffer *buffer)
{
  const Babl *format;
  gint        bpp;
  goffset     saved_pos;
  goffset     offset;
  goffset     offset2;
  goffset     max_data_length;
  gint        n_tile_rows;
  gint        n_tile_cols;
  guint       ntiles;
  gint        width;
  gint        height;
  gint        i;
  gint        fail;

  format = gegl_buffer_get_format (buffer);
  bpp    = babl_format_get_bytes_per_pixel (format);

  xcf_read_int32 (info, (guint32 *) &width,  1);
  xcf_read_int32 (info, (guint32 *) &height, 1);

  if (width  != gegl_buffer_get_width (buffer) ||
      height != gegl_buffer_get_height (buffer))
    return FALSE;

  /* maximal allowable size of on-disk tile data.  make it somewhat bigger than
   * the uncompressed tile size, to allow for the possibility of negative
   * compression.
   */
  max_data_length = XCF_TILE_WIDTH * XCF_TILE_HEIGHT * bpp *
                    XCF_TILE_MAX_DATA_LENGTH_FACTOR /* = 1.5, currently */;

  /* read in the first tile offset.
   *  if it is '0', then this tile level is empty
   *  and we can simply return.
   */
  if (xcf_read_offset (info, &offset, 1) < info->bytes_per_offset)
    {
      GIMP_LOG (XCF, "Failed to read tile offset"
                " at offset: %" G_GOFFSET_FORMAT, info->cp);
      return FALSE;
    }
  if (offset == 0)
    return TRUE;

  n_tile_rows = gimp_gegl_buffer_get_n_tile_rows (buffer, XCF_TILE_HEIGHT);
  n_tile_cols = gimp_gegl_buffer_get_n_tile_cols (buffer, XCF_TILE_WIDTH);

  ntiles = n_tile_rows * n_tile_cols;
  for (i = 0; i < ntiles; i++)
    {
      GeglRectangle rect;

      fail = FALSE;

      if (offset == 0)
        {
          gimp_message_literal (info->gimp, G_OBJECT (info->progress),
                                GIMP_MESSAGE_ERROR,
                                "not enough tiles found in level");
          return FALSE;
        }

      /* save the current position as it is where the
       *  next tile offset is stored.
       */
      saved_pos = info->cp;

      /* read in the offset of the next tile so we can calculate the amount
       * of data needed for this tile
       */
      if (xcf_read_offset (info, &offset2, 1) < info->bytes_per_offset)
        {
          GIMP_LOG (XCF, "Failed to read tile offset"
                    " at offset: %" G_GOFFSET_FORMAT, info->cp);
          return FALSE;
        }

      /* if the offset is 0 then we need to read in the maximum possible
       * allowing for negative compression
       */
      if (offset2 == 0)
        offset2 = offset + max_data_length;

      /* seek to the tile offset */
      if (! xcf_seek_pos (info, offset, NULL))
        return FALSE;

      if (offset2 < offset || offset2 - offset > max_data_length)
        {
          gimp_message (info->gimp, G_OBJECT (info->progress),
                        GIMP_MESSAGE_ERROR,
                        "invalid tile data length: %" G_GOFFSET_FORMAT,
                        offset2 - offset);
          return FALSE;
        }

      /* get buffer rectangle to write to */
      gimp_gegl_buffer_get_tile_rect (buffer,
                                      XCF_TILE_WIDTH, XCF_TILE_HEIGHT,
                                      i, &rect);

      GIMP_LOG (XCF, "loading tile %d/%d", i + 1, ntiles);

      /* read in the tile */
      switch (info->compression)
        {
        case COMPRESS_NONE:
          if (! xcf_load_tile (info, buffer, &rect, format))
            fail = TRUE;
          break;
        case COMPRESS_RLE:
          if (! xcf_load_tile_rle (info, buffer, &rect, format,
                                   offset2 - offset))
            fail = TRUE;
          break;
        case COMPRESS_ZLIB:
          if (! xcf_load_tile_zlib (info, buffer, &rect, format,
                                    offset2 - offset))
            fail = TRUE;
          break;
        case COMPRESS_FRACTAL:
          g_printerr ("xcf: fractal compression unimplemented. "
                      "Possibly corrupt XCF file.");
          fail = TRUE;
          break;
        default:
          g_printerr ("xcf: unknown compression. "
                      "Possibly corrupt XCF file.");
          fail = TRUE;
          break;
        }

      if (fail)
        return FALSE;

      GIMP_LOG (XCF, "loaded tile %d/%d", i + 1, ntiles);

      /* restore the saved position so we'll be ready to
       *  read the next offset.
       */
      if (!xcf_seek_pos (info, saved_pos, NULL))
        return FALSE;

      /* read in the offset of the next tile */
      if (xcf_read_offset (info, &offset, 1) < info->bytes_per_offset)
        {
          GIMP_LOG (XCF, "Failed to read tile offset"
                    " at offset: %" G_GOFFSET_FORMAT, info->cp);
          return FALSE;
        }
    }

  if (offset != 0)
    {
      gimp_message (info->gimp, G_OBJECT (info->progress), GIMP_MESSAGE_ERROR,
                    "encountered garbage after reading level: %" G_GOFFSET_FORMAT,
                    offset);
      return FALSE;
    }

  return TRUE;
}

static gboolean
xcf_load_tile (XcfInfo       *info,
               GeglBuffer    *buffer,
               GeglRectangle *tile_rect,
               const Babl    *format)
{
  gint    bpp       = babl_format_get_bytes_per_pixel (format);
  gint    tile_size = bpp * tile_rect->width * tile_rect->height;
  guchar *tile_data = g_alloca (tile_size);

  if (info->file_version <= 11)
    {
      xcf_read_int8 (info, tile_data, tile_size);
    }
  else
    {
      gint n_components = babl_format_get_n_components (format);

      xcf_read_component (info, bpp / n_components, tile_data,
                          tile_size / bpp * n_components);
    }

  if (! xcf_data_is_zero (tile_data, tile_size))
    {
      gegl_buffer_set (buffer, tile_rect, 0, format, tile_data,
                       GEGL_AUTO_ROWSTRIDE);
    }

  return TRUE;
}

static gboolean
xcf_load_tile_rle (XcfInfo       *info,
                   GeglBuffer    *buffer,
                   GeglRectangle *tile_rect,
                   const Babl    *format,
                   gint           data_length)
{
  gint    bpp       = babl_format_get_bytes_per_pixel (format);
  gint    tile_size = bpp * tile_rect->width * tile_rect->height;
  guchar *tile_data = g_alloca (tile_size);
  guchar  nonzero   = FALSE;
  gsize   bytes_read;
  gint    i;
  guchar *xcfdata;
  guchar *xcfodata;
  guchar *xcfdatalimit;

  /* Workaround for bug #357809: avoid crashing on g_malloc() and skip
   * this tile (return TRUE without storing data) as if it did not
   * contain any data.  It is better than returning FALSE, which would
   * skip the whole hierarchy while there may still be some valid
   * tiles in the file.
   */
  if (data_length <= 0)
    return TRUE;

  xcfdata = xcfodata = g_alloca (data_length);

  /* we have to read directly instead of xcf_read_* because we may be
   * reading past the end of the file here
   */
  g_input_stream_read_all (info->input, xcfdata, data_length,
                           &bytes_read, NULL, NULL);
  info->cp += bytes_read;

  if (bytes_read == 0)
    return TRUE;

  xcfdatalimit = &xcfodata[bytes_read - 1];

  for (i = 0; i < bpp; i++)
    {
      guchar *data  = tile_data + i;
      gint    size  = tile_rect->width * tile_rect->height;
      gint    count = 0;
      guchar  val;
      gint    length;
      gint    j;

      while (size > 0)
        {
          if (xcfdata > xcfdatalimit)
            {
              goto bogus_rle;
            }

          val = *xcfdata++;

          length = val;
          if (length >= 128)
            {
              length = 255 - (length - 1);
              if (length == 128)
                {
                  if (xcfdata >= xcfdatalimit)
                    {
                      goto bogus_rle;
                    }

                  length = (*xcfdata << 8) + xcfdata[1];
                  xcfdata += 2;
                }

              count += length;
              size -= length;

              if (size < 0)
                {
                  goto bogus_rle;
                }

              if (&xcfdata[length-1] > xcfdatalimit)
                {
                  goto bogus_rle;
                }

              while (length-- > 0)
                {
                  *data = *xcfdata++;
                  nonzero |= *data;
                  data += bpp;
                }
            }
          else
            {
              length += 1;
              if (length == 128)
                {
                  if (xcfdata >= xcfdatalimit)
                    {
                      goto bogus_rle;
                    }

                  length = (*xcfdata << 8) + xcfdata[1];
                  xcfdata += 2;
                }

              count += length;
              size -= length;

              if (size < 0)
                {
                  goto bogus_rle;
                }

              if (xcfdata > xcfdatalimit)
                {
                  goto bogus_rle;
                }

              val = *xcfdata++;
              nonzero |= val;

              for (j = 0; j < length; j++)
                {
                  *data = val;
                  data += bpp;
                }
            }
        }
    }

  if (nonzero)
    {
      if (info->file_version >= 12)
        {
          gint n_components = babl_format_get_n_components (format);

          xcf_read_from_be (bpp / n_components, tile_data,
                            tile_size / bpp * n_components);
        }

      gegl_buffer_set (buffer, tile_rect, 0, format, tile_data,
                       GEGL_AUTO_ROWSTRIDE);
    }

  return TRUE;

 bogus_rle:
  return FALSE;
}

static gboolean
xcf_load_tile_zlib (XcfInfo       *info,
                    GeglBuffer    *buffer,
                    GeglRectangle *tile_rect,
                    const Babl    *format,
                    gint           data_length)
{
  z_stream  strm;
  int       action;
  int       status;
  gint      bpp       = babl_format_get_bytes_per_pixel (format);
  gint      tile_size = bpp * tile_rect->width * tile_rect->height;
  guchar   *tile_data = g_alloca (tile_size);
  gsize     bytes_read;
  guchar   *xcfdata;

  /* Workaround for bug #357809: avoid crashing on g_malloc() and skip
   * this tile (return TRUE without storing data) as if it did not
   * contain any data.  It is better than returning FALSE, which would
   * skip the whole hierarchy while there may still be some valid
   * tiles in the file.
   */
  if (data_length <= 0)
    return TRUE;

  xcfdata = g_alloca (data_length);

  /* we have to read directly instead of xcf_read_* because we may be
   * reading past the end of the file here
   */
  g_input_stream_read_all (info->input, xcfdata, data_length,
                           &bytes_read, NULL, NULL);
  info->cp += bytes_read;

  if (bytes_read == 0)
    return TRUE;

  strm.next_out  = tile_data;
  strm.avail_out = tile_size;

  strm.zalloc    = Z_NULL;
  strm.zfree     = Z_NULL;
  strm.opaque    = Z_NULL;
  strm.next_in   = xcfdata;
  strm.avail_in  = bytes_read;

  /* Initialize the stream decompression. */
  status = inflateInit (&strm);
  if (status != Z_OK)
    return FALSE;

  action = Z_NO_FLUSH;

  while (status == Z_OK)
    {
      if (strm.avail_in == 0)
        {
          action = Z_FINISH;
        }

      status = inflate (&strm, action);

      if (status == Z_STREAM_END)
        {
          /* All the data was successfully decoded. */
          break;
        }
      else if (status == Z_BUF_ERROR)
        {
          g_printerr ("xcf: decompressed tile bigger than the expected size.");
          inflateEnd (&strm);
          return FALSE;
        }
      else if (status != Z_OK)
        {
          g_printerr ("xcf: tile decompression failed: %s", zError (status));
          inflateEnd (&strm);
          return FALSE;
        }
    }

  if (! xcf_data_is_zero (tile_data, tile_size))
    {
      if (info->file_version >= 12)
        {
          gint n_components = babl_format_get_n_components (format);

          xcf_read_from_be (bpp / n_components, tile_data,
                            tile_size / bpp * n_components);
        }

      gegl_buffer_set (buffer, tile_rect, 0, format, tile_data,
                       GEGL_AUTO_ROWSTRIDE);
    }

  inflateEnd (&strm);

  return TRUE;
}

static GimpParasite *
xcf_load_parasite (XcfInfo *info)
{
  GimpParasite *parasite = NULL;
  gchar        *name;
  guint32       flags;
  guint32       size, size_read;
  gpointer      data;

  xcf_read_string (info, &name,  1);
  xcf_read_int32  (info, &flags, 1);
  xcf_read_int32  (info, &size,  1);

  GIMP_LOG (XCF, "Parasite name: %s, flags: %d, size: %d", name, flags, size);

  if (size > MAX_XCF_PARASITE_DATA_LEN)
    {
      g_printerr ("Maximum parasite data length (%ld bytes) exceeded. "
                  "Possibly corrupt XCF file.", MAX_XCF_PARASITE_DATA_LEN);
      g_free (name);
      return NULL;
    }

  if (!name)
    {
      g_printerr ("Parasite has no name! Possibly corrupt XCF file.\n");
      return NULL;
    }

  data = g_new (gchar, size);
  size_read = xcf_read_int8 (info, data, size);

  if (size_read != size)
    {
      g_printerr ("Incorrect parasite data size: read %u bytes instead of %u. "
                  "Possibly corrupt XCF file.\n",
                  size_read, size);
    }
  else
    {
      parasite = gimp_parasite_new (name, flags, size, data);
    }

  g_free (name);
  g_free (data);

  return parasite;
}

/* Old paths are the PROP_PATHS property, even older than PROP_VECTORS. */
static gboolean
xcf_load_old_paths (XcfInfo   *info,
                    GimpImage *image)
{
  guint32   num_paths;
  guint32   last_selected_row;
  GimpPath *active_path;

  xcf_read_int32 (info, &last_selected_row, 1);
  xcf_read_int32 (info, &num_paths,         1);

  GIMP_LOG (XCF, "Number of old paths: %u", num_paths);

  while (num_paths-- > 0)
    if (! xcf_load_old_path (info, image))
      return FALSE;

  active_path =
    GIMP_PATH (gimp_container_get_child_by_index (gimp_image_get_paths (image),
                                                     last_selected_row));

  if (active_path)
    {
      GList *list = g_list_prepend (NULL, active_path);
      gimp_image_set_selected_paths (image, list);
      g_list_free (list);
    }

  return TRUE;
}

static gboolean
xcf_load_old_path (XcfInfo   *info,
                   GimpImage *image)
{
  gchar                  *name;
  guint32                 locked;
  guint8                  state;
  guint32                 closed;
  guint32                 num_points;
  guint32                 version; /* changed from num_paths */
  GimpTattoo              tattoo = 0;
  GimpPath               *path;
  GimpPathCompatPoint    *points;
  gint                    i;

  xcf_read_string (info, &name,       1);
  xcf_read_int32  (info, &locked,     1);
  xcf_read_int8   (info, &state,      1);
  xcf_read_int32  (info, &closed,     1);
  xcf_read_int32  (info, &num_points, 1);
  xcf_read_int32  (info, &version,    1);

  if (version == 2)
    {
      guint32 dummy;

      /* Had extra type field and points are stored as doubles */
      xcf_read_int32 (info, (guint32 *) &dummy, 1);
    }
  else if (version == 3)
    {
      guint32 dummy;

      /* Has extra tattoo field */
      xcf_read_int32 (info, (guint32 *) &dummy,  1);
      xcf_read_int32 (info, (guint32 *) &tattoo, 1);
    }
  else if (version != 1)
    {
      g_printerr ("Unknown path type (version: %u). Possibly corrupt XCF file.\n", version);

      g_free (name);
      return FALSE;
    }

  /* skip empty compatibility paths */
  if (num_points == 0)
    {
      g_free (name);
      return FALSE;
    }

  points = g_new0 (GimpPathCompatPoint, num_points);

  for (i = 0; i < num_points; i++)
    {
      if (version == 1)
        {
          gint32 x;
          gint32 y;

          xcf_read_int32 (info, &points[i].type, 1);
          xcf_read_int32 (info, (guint32 *) &x,  1);
          xcf_read_int32 (info, (guint32 *) &y,  1);

          points[i].x = x;
          points[i].y = y;
        }
      else
        {
          gfloat x;
          gfloat y;

          xcf_read_int32 (info, &points[i].type, 1);
          xcf_read_float (info, &x,              1);
          xcf_read_float (info, &y,              1);

          points[i].x = x;
          points[i].y = y;
        }
    }

  path = gimp_path_compat_new (image, name, points, num_points, closed);

  g_free (name);
  g_free (points);

  if (locked)
    info->linked_paths = g_list_prepend (info->linked_paths, path);

  if (tattoo)
    gimp_item_set_tattoo (GIMP_ITEM (path), tattoo);

  gimp_image_add_path (image, path,
                       NULL, /* can't be a tree */
                       gimp_container_get_n_children (gimp_image_get_paths (image)),
                       FALSE);

  return TRUE;
}

/* Old vectors are the PROP_VECTORS property up to all GIMP 2.10 versions. */
static gboolean
xcf_load_old_vectors (XcfInfo   *info,
                      GimpImage *image)
{
  guint32   version;
  guint32   active_index;
  guint32   num_paths;
  GimpPath *active_path;

#ifdef GIMP_XCF_PATH_DEBUG
  g_printerr ("xcf_load_old_vectors\n");
#endif

  xcf_read_int32 (info, &version, 1);

  if (version != 1)
    {
      gimp_message (info->gimp, G_OBJECT (info->progress),
                    GIMP_MESSAGE_WARNING,
                    "Unknown path version: %d (skipping)", version);
      return FALSE;
    }

  xcf_read_int32 (info, &active_index, 1);
  xcf_read_int32 (info, &num_paths,    1);

#ifdef GIMP_XCF_PATH_DEBUG
  g_printerr ("%d paths (active: %d)\n", num_paths, active_index);
#endif

  while (num_paths-- > 0)
    if (! xcf_load_old_vector (info, image))
      return FALSE;

  /* FIXME tree */
  active_path =
    GIMP_PATH (gimp_container_get_child_by_index (gimp_image_get_paths (image),
                                                  active_index));

  if (active_path)
    {
      GList *list = g_list_prepend (NULL, active_path);
      gimp_image_set_selected_paths (image, list);
      g_list_free (list);
    }

#ifdef GIMP_XCF_PATH_DEBUG
  g_printerr ("xcf_load_old_vectors: loaded %d bytes\n", info->cp - base);
#endif
  return TRUE;
}

static gboolean
xcf_load_old_vector (XcfInfo   *info,
                     GimpImage *image)
{
  gchar      *name;
  GimpTattoo  tattoo = 0;
  guint32     visible;
  guint32     linked;
  guint32     num_parasites;
  guint32     num_strokes;
  GimpPath   *path;
  gint        i;

#ifdef GIMP_XCF_PATH_DEBUG
  g_printerr ("xcf_load_old_vector\n");
#endif

  xcf_read_string (info, &name,          1);
  xcf_read_int32  (info, &tattoo,        1);
  xcf_read_int32  (info, &visible,       1);
  xcf_read_int32  (info, &linked,        1);
  xcf_read_int32  (info, &num_parasites, 1);
  xcf_read_int32  (info, &num_strokes,   1);

#ifdef GIMP_XCF_PATH_DEBUG
  g_printerr ("name: %s, tattoo: %d, visible: %d, linked: %d, "
              "num_parasites %d, num_strokes %d\n",
              name, tattoo, visible, linked, num_parasites, num_strokes);
#endif

  path = gimp_path_new (image, name);
  g_free (name);

  gimp_item_set_visible (GIMP_ITEM (path), visible, FALSE);
  if (linked)
    info->linked_paths = g_list_prepend (info->linked_paths, path);

  if (tattoo)
    gimp_item_set_tattoo (GIMP_ITEM (path), tattoo);

  for (i = 0; i < num_parasites; i++)
    {
      GimpParasite *parasite = xcf_load_parasite (info);
      GError       *error    = NULL;

      if (! parasite)
        return FALSE;

      if (! gimp_item_parasite_validate (GIMP_ITEM (path), parasite, &error))
        {
          gimp_message (info->gimp, G_OBJECT (info->progress),
                        GIMP_MESSAGE_WARNING,
                        "Warning, invalid vectors parasite in XCF file: %s",
                        error->message);
          g_clear_error (&error);
        }
      else
        {
          gimp_item_parasite_attach (GIMP_ITEM (path), parasite, FALSE);
        }

      gimp_parasite_free (parasite);
    }

  for (i = 0; i < num_strokes; i++)
    {
      guint32      stroke_type_id;
      guint32      closed;
      guint32      num_axes;
      guint32      num_control_points;
      guint32      type;
      gfloat       coords[13] = GIMP_COORDS_DEFAULT_VALUES;
      GimpStroke  *stroke;
      gint         j;

      GimpValueArray *control_points;
      GValue          value  = G_VALUE_INIT;
      GimpAnchor      anchor = { { 0, } };
      GType           stroke_type;

      g_value_init (&value, GIMP_TYPE_ANCHOR);

      xcf_read_int32 (info, &stroke_type_id,     1);
      xcf_read_int32 (info, &closed,             1);
      xcf_read_int32 (info, &num_axes,           1);
      xcf_read_int32 (info, &num_control_points, 1);

#ifdef GIMP_XCF_PATH_DEBUG
      g_printerr ("stroke_type: %d, closed: %d, num_axes %d, len %d\n",
                  stroke_type_id, closed, num_axes, num_control_points);
#endif

      switch (stroke_type_id)
        {
        case XCF_STROKETYPE_BEZIER_STROKE:
          stroke_type = GIMP_TYPE_BEZIER_STROKE;
          break;

        default:
          g_printerr ("skipping unknown stroke type\n");
          xcf_seek_pos (info,
                        info->cp + 4 * num_axes * num_control_points,
                        NULL);
          continue;
        }

      if (num_axes < 2 || num_axes > 6)
        {
          g_printerr ("bad number of axes in stroke description\n");
          return FALSE;
        }

      control_points = gimp_value_array_new (num_control_points);

      anchor.selected = FALSE;

      for (j = 0; j < num_control_points; j++)
        {
          xcf_read_int32 (info, &type,  1);
          xcf_read_float (info, coords, num_axes);

          anchor.type              = type;
          anchor.position.x        = coords[0];
          anchor.position.y        = coords[1];
          anchor.position.pressure = coords[2];
          anchor.position.xtilt    = coords[3];
          anchor.position.ytilt    = coords[4];
          anchor.position.wheel    = coords[5];

          g_value_set_boxed (&value, &anchor);
          gimp_value_array_append (control_points, &value);

#ifdef GIMP_XCF_PATH_DEBUG
          g_printerr ("Anchor: %d, (%f, %f, %f, %f, %f, %f)\n", type,
                      coords[0], coords[1], coords[2], coords[3],
                      coords[4], coords[5]);
#endif
        }

      g_value_unset (&value);

      stroke = g_object_new (stroke_type,
                             "closed",         closed,
                             "control-points", control_points,
                             NULL);

      gimp_path_stroke_add (path, stroke);

      g_object_unref (stroke);
      gimp_value_array_unref (control_points);
    }

  gimp_image_add_path (image, path,
                       NULL, /* FIXME tree */
                       gimp_container_get_n_children (gimp_image_get_paths (image)),
                       FALSE);

  return TRUE;
}

static gboolean
xcf_skip_unknown_prop (XcfInfo *info,
                       gsize   size)
{
  guint8 buf[16];
  guint  amount;

  while (size > 0)
    {
      if (g_input_stream_is_closed (info->input))
        return FALSE;

      amount = MIN (16, size);
      amount = xcf_read_int8 (info, buf, amount);
      if (amount == 0)
        return FALSE;

      size -= amount;
    }

  return TRUE;
}

static gboolean
xcf_item_path_is_parent (GList *path,
                         GList *parent_path)
{
  GList *iter        = path;
  GList *parent_iter = parent_path;

  if (g_list_length (parent_path) >= g_list_length (path))
    return FALSE;

  while (iter && parent_iter)
    {
      if (iter->data != parent_iter->data)
        return FALSE;

      iter = iter->next;
      parent_iter = parent_iter->next;
    }

  return TRUE;
}

static void
xcf_fix_item_path (GimpLayer  *layer,
                   GList     **path,
                   GList      *broken_paths)
{
  GList *iter;

  for (iter = broken_paths; iter; iter = iter->next)
    {
      if (xcf_item_path_is_parent (*path, iter->data))
        {
          /* Not much to do when the absent path is a parent. */
          g_printerr ("%s: layer '%s' moved to layer tree root because of missing parent.",
                      G_STRFUNC, gimp_object_get_name (layer));
          g_clear_pointer (path, g_list_free);
          return;
        }
    }

  /* Check if a parent of path, or path itself is on the same
   * tree level as any broken path; and if so, and if the broken path is
   * in a lower position in the item group, decrement it.
   */
  for (iter = broken_paths; iter; iter = iter->next)
    {
      GList *broken_path = iter->data;
      GList *iter1       = *path;
      GList *iter2       = broken_path;

      if (g_list_length (broken_path) > g_list_length (*path))
        continue;

      while (iter1 && iter2)
        {
          if (iter2->next && iter1->data != iter2->data)
            /* Paths diverged before reaching iter2 leaf. */
            break;

          if (iter2->next)
            {
              iter1 = iter1->next;
              iter2 = iter2->next;
              continue;
            }

          if (GPOINTER_TO_UINT (iter2->data) < GPOINTER_TO_UINT (iter1->data))
            iter1->data = GUINT_TO_POINTER (GPOINTER_TO_UINT (iter1->data) - 1);

          break;
        }
    }
}
