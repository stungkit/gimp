/*   GIMP - The GNU Image Manipulation Program
 *   Copyright (C) 1995 Spencer Kimball and Peter Mattis
 *
 *   Amiga IFF plug-in
 *
 *   Copyright (C) 2023 Alex S.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#include "config.h"

#include <string.h>
#include <errno.h>

#include <glib/gstdio.h>

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>

#include <libilbm/ilbm.h>
#include <libilbm/ilbmimage.h>
#include <libilbm/interleave.h>
#include <libilbm/byterun.h>

#include "libgimp/stdplugins-intl.h"


#define LOAD_PROC      "file-iff-load"
#define PLUG_IN_BINARY "file-iff"
#define PLUG_IN_ROLE   "gimp-file-iff"

typedef struct _Iff      Iff;
typedef struct _IffClass IffClass;

struct _Iff
{
  GimpPlugIn      parent_instance;
};

struct _IffClass
{
  GimpPlugInClass parent_class;
};


#define IFF_TYPE  (iff_get_type ())
#define IFF (obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), IFF_TYPE, Iff))

GType                   iff_get_type         (void) G_GNUC_CONST;


static GList          * iff_query_procedures (GimpPlugIn           *plug_in);
static GimpProcedure  * iff_create_procedure (GimpPlugIn           *plug_in,
                                              const gchar          *name);

static GimpValueArray * iff_load             (GimpProcedure        *procedure,
                                              GimpRunMode           run_mode,
                                              GFile                *file,
                                              const GimpValueArray *args,
                                              gpointer              run_data);

static GimpImage      * load_image           (GFile                *file,
                                              GObject              *config,
                                              GimpRunMode           run_mode,
                                              GError              **error);

static void             deleave_indexed_row  (IFF_UByte            *bitplanes,
                                              guchar               *pixel_row,
                                              gint                  width,
                                              gint                  nPlanes);

static void             pbm_row              (IFF_UByte            *bitplanes,
                                              guchar               *pixel_row,
                                              gint                  width);


G_DEFINE_TYPE (Iff, iff, GIMP_TYPE_PLUG_IN)

GIMP_MAIN (IFF_TYPE)
DEFINE_STD_SET_I18N


static void
iff_class_init (IffClass *klass)
{
  GimpPlugInClass *plug_in_class = GIMP_PLUG_IN_CLASS (klass);

  plug_in_class->query_procedures = iff_query_procedures;
  plug_in_class->create_procedure = iff_create_procedure;
  plug_in_class->set_i18n         = STD_SET_I18N;
}

static void
iff_init (Iff *iff)
{
}

static GList *
iff_query_procedures (GimpPlugIn *plug_in)
{
  GList *list = NULL;

  list = g_list_append (list, g_strdup (LOAD_PROC));

  return list;
}

static GimpProcedure *
iff_create_procedure (GimpPlugIn  *plug_in,
                      const gchar *name)
{
  GimpProcedure *procedure = NULL;

  if (! strcmp (name, LOAD_PROC))
    {
      procedure = gimp_load_procedure_new (plug_in, name,
                                           GIMP_PDB_PROC_TYPE_PLUGIN,
                                           iff_load, NULL, NULL);

      gimp_procedure_set_menu_label (procedure,
                                     N_("Amiga IFF"));

      gimp_procedure_set_documentation (procedure,
                                        _("Load file in the IFF file format"),
                                        _("Load file in the IFF file format"),
                                        name);
      gimp_procedure_set_attribution (procedure,
                                      "Alex S.",
                                      "Alex S.",
                                      "2023");

      gimp_file_procedure_set_mime_types (GIMP_FILE_PROCEDURE (procedure),
                                          "image/iff");
      gimp_file_procedure_set_extensions (GIMP_FILE_PROCEDURE (procedure),
                                          "iff,ilbm,acbm");
      gimp_file_procedure_set_magics (GIMP_FILE_PROCEDURE (procedure),
                                      "0,string,FORM");
    }

  return procedure;
}

static GimpValueArray *
iff_load (GimpProcedure        *procedure,
          GimpRunMode           run_mode,
          GFile                *file,
          const GimpValueArray *args,
          gpointer              run_data)
{
  GimpProcedureConfig *config;
  GimpValueArray      *return_vals;
  GimpImage           *image;
  GError              *error = NULL;

  gegl_init (NULL, NULL);

  config = gimp_procedure_create_config (procedure);
  gimp_procedure_config_begin_run (config, NULL, run_mode, args);

  image = load_image (file, G_OBJECT (config), run_mode, &error);

  if (! image)
    return gimp_procedure_new_return_values (procedure,
                                             GIMP_PDB_EXECUTION_ERROR,
                                             error);

  gimp_procedure_config_end_run (config, GIMP_PDB_SUCCESS);
  g_object_unref (config);

  return_vals = gimp_procedure_new_return_values (procedure,
                                                  GIMP_PDB_SUCCESS,
                                                  NULL);

  GIMP_VALUES_SET_IMAGE (return_vals, 1, image);

  return return_vals;
}

static GimpImage *
load_image (GFile        *file,
            GObject      *config,
            GimpRunMode   run_mode,
            GError      **error)
{
  GimpImage   *image = NULL;
  GimpLayer   *layer;
  GeglBuffer  *buffer;
  FILE        *fp;
  guint        imagesLength;
  IFF_Chunk   *chunk;
  ILBM_Image **iff_image;

  fp = g_fopen (g_file_peek_path (file), "rb");

  if (! fp)
    {
      g_set_error (error, G_FILE_ERROR, g_file_error_from_errno (errno),
                   _("Could not open '%s' for reading: %s"),
                   gimp_file_get_utf8_name (file), g_strerror (errno));
      return NULL;
    }

  fclose (fp);

  chunk = ILBM_read (g_file_peek_path (file));
  iff_image = ILBM_extractImages (chunk, &imagesLength);

  for (gint i = 0; i < imagesLength; i++)
    {
      ILBM_Image        *true_image   = iff_image[i];
      /* Struct representing bitmap header properties */
      ILBM_BitMapHeader *bitMapHeader = true_image->bitMapHeader;
      /* Struct containing the color palette */
      ILBM_ColorMap     *colorMap     = true_image->colorMap;
      IFF_UByte         *bitplanes;
      guchar             gimp_cmap[768];  /* Max index is (2^nplanes) - 1 */
      gint               width;
      gint               height;
      gint               nPlanes;
      gint               row_length;
      gint               y_height     = 0;

      if (! true_image || ! bitMapHeader)
        {
          g_message (_("Invalid or missing ILBM image"));
          return image;
        }

      /* Convert ACBM files to ILBM format */
      if (ILBM_imageIsACBM (true_image))
        ILBM_convertACBMToILBM (true_image);

      width      = bitMapHeader->w;
      height     = bitMapHeader->h;
      nPlanes    = bitMapHeader->nPlanes;
      row_length = (width + 15) / 16;

      /* Load palette if it exists */
      if (colorMap)
        {
          for (gint j = 0; j < colorMap->colorRegisterLength; j++)
            {
              gimp_cmap[j * 3]     = colorMap->colorRegister[j].red;
              gimp_cmap[j * 3 + 1] = colorMap->colorRegister[j].green;
              gimp_cmap[j * 3 + 2] = colorMap->colorRegister[j].blue;
            }
        }

      /* Ignored if already uncompressed */
      ILBM_unpackByteRun (true_image);

      image = gimp_image_new (width, height, GIMP_INDEXED);

      layer = gimp_layer_new (image, _("Background"), width, height,
                              GIMP_INDEXED_IMAGE, 100,
                              gimp_image_get_default_new_layer_mode (image));
      gimp_image_insert_layer (image, layer, NULL, 0);

      buffer = gimp_drawable_get_buffer (GIMP_DRAWABLE (layer));

      bitplanes = true_image->body->chunkData;
      /* Loading rows */
      for (gint j = 0; j < height; j++)
        {
          guchar *pixel_row;

          pixel_row = g_malloc (width * sizeof (guchar));

          /* PBM uses one byte per pixel index */
          if (ILBM_imageIsPBM (true_image))
            pbm_row (bitplanes, pixel_row, width);
          else
            deleave_indexed_row (bitplanes, pixel_row, width, nPlanes);

          bitplanes += (row_length * 2 * nPlanes);

          gegl_buffer_set (buffer, GEGL_RECTANGLE (0, y_height, width, 1), 0,
                           NULL, pixel_row, GEGL_AUTO_ROWSTRIDE);

          y_height++;
          g_free (pixel_row);
        }

      gimp_image_set_colormap (image, gimp_cmap,
                               colorMap->colorRegisterLength);

      g_object_unref (buffer);
    }
  ILBM_freeImages (iff_image, imagesLength);

  return image;
}

static void
deleave_indexed_row (IFF_UByte            *bitplanes,
                     guchar               *pixel_row,
                     gint                  width,
                     gint                  nPlanes)
{
  guchar index[width];
  gint   row_length = ((width + 15) / 16) * 2;

  /* Initialize index array */
  for (gint i = 0; i < width; i++)
    index[i] = 0;

  /* Deleave rows */
  for (gint i = 0; i < row_length; i++)
    {
      for (gint j = 0; j < 8; j++)
        {
          guint8 bitmask = (1 << (8 - j)) - (1 << (7 - j));

          for (gint k = 0; k < nPlanes; k++)
            {
               guint8 update = (1 << (k + 1)) - (1 << (k));

               if (bitplanes[i + (row_length * k)] & bitmask)
                 index[j + (i * 8)] += update;
            }
        }
    }

  /* Associate palette with pixels */
  for (gint i = 0; i < width; i++)
    {
      pixel_row[i] = index[i];
    }
}

static void
pbm_row (IFF_UByte *bitplanes,
         guchar    *pixel_row,
         gint       width)
{
  for (gint i = 0; i < width; i++)
    pixel_row[i] = bitplanes[i];
}
