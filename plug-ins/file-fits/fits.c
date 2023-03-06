/* GIMP - The GNU Image Manipulation Program
 * Copyright (C) 1995 Spencer Kimball and Peter Mattis
 * FITS file plugin
 * reading and writing code Copyright (C) 1997 Peter Kirchgessner
 * e-mail: peter@kirchgessner.net, WWW: http://www.kirchgessner.net
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
 *
 */

/* Event history:
 * V 1.00, PK, 05-May-97: Creation
 * V 1.01, PK, 19-May-97: Problem with compilation on Irix fixed
 * V 1.02, PK, 08-Jun-97: Bug with saving gray images fixed
 * V 1.03, PK, 05-Oct-97: Parse rc-file
 * V 1.04, PK, 12-Oct-97: No progress bars for non-interactive mode
 * V 1.05, nn, 20-Dec-97: Initialize image_ID in run()
 * V 1.06, PK, 21-Nov-99: Internationalization
 *                        Fix bug with gimp_export_image()
 *                        (moved it from load to save)
 * V 1.07, PK, 16-Aug-06: Fix problems with internationalization
 *                        (writing 255,0 instead of 255.0)
 *                        Fix problem with not filling up properly last record
 */

#include "config.h"

#include <string.h>
#include <errno.h>

#include <glib/gstdio.h>

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>

#include "fits-io.h"
#include <fitsio.h>

#include "libgimp/stdplugins-intl.h"


#define LOAD_PROC      "file-fits-load"
#define SAVE_PROC      "file-fits-save"
#define PLUG_IN_BINARY "file-fits"
#define PLUG_IN_ROLE   "gimp-file-fits"


typedef struct _Fits      Fits;
typedef struct _FitsClass FitsClass;

struct _Fits
{
  GimpPlugIn      parent_instance;
};

struct _FitsClass
{
  GimpPlugInClass parent_class;
};


#define FITS_TYPE  (fits_get_type ())
#define FITS (obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), FITS_TYPE, Fits))

GType                   fits_get_type         (void) G_GNUC_CONST;

typedef struct
{
  gint   naxis;
  glong  naxisn[3];
  gint   bitpix;
  gint   bpp;
  gint   datatype;
} FitsHduData;

static GList          * fits_query_procedures (GimpPlugIn           *plug_in);
static GimpProcedure  * fits_create_procedure (GimpPlugIn           *plug_in,
                                               const gchar          *name);

static GimpValueArray * fits_load             (GimpProcedure        *procedure,
                                               GimpRunMode           run_mode,
                                               GFile                *file,
                                               const GimpValueArray *args,
                                               gpointer              run_data);
static GimpValueArray * fits_save             (GimpProcedure        *procedure,
                                               GimpRunMode           run_mode,
                                               GimpImage            *image,
                                               gint                  n_drawables,
                                               GimpDrawable        **drawables,
                                               GFile                *file,
                                               const GimpValueArray *args,
                                               gpointer              run_data);

static GimpImage      * load_image            (GFile              *file,
                                               GObject            *config,
                                               GimpRunMode         run_mode,
                                               GError            **error);
static gint             save_image            (GFile              *file,
                                               GimpImage          *image,
                                               GimpDrawable       *drawable,
                                               GError            **error);

static FitsHduList    * create_fits_header    (FitsFile           *ofp,
                                               guint               width,
                                               guint               height,
                                               guint               channels,
                                               guint               bitpix);

static gint             save_fits             (FitsFile           *ofp,
                                               GimpImage          *image,
                                               GimpDrawable       *drawable);

static GimpImage      * create_new_image      (GFile              *file,
                                               guint               pagenum,
                                               guint               width,
                                               guint               height,
                                               GimpImageBaseType   itype,
                                               GimpImageType       dtype,
                                               GimpPrecision       iprecision,
                                               GimpLayer         **layer,
                                               GeglBuffer        **buffer);

static gboolean         load_dialog           (GimpProcedure      *procedure,
                                               GObject            *config);
static void             show_fits_errors      (gint                status);


G_DEFINE_TYPE (Fits, fits, GIMP_TYPE_PLUG_IN)

GIMP_MAIN (FITS_TYPE)
DEFINE_STD_SET_I18N


static void
fits_class_init (FitsClass *klass)
{
  GimpPlugInClass *plug_in_class = GIMP_PLUG_IN_CLASS (klass);

  plug_in_class->query_procedures = fits_query_procedures;
  plug_in_class->create_procedure = fits_create_procedure;
  plug_in_class->set_i18n         = STD_SET_I18N;
}

static void
fits_init (Fits *fits)
{
}

static GList *
fits_query_procedures (GimpPlugIn *plug_in)
{
  GList *list = NULL;

  list = g_list_append (list, g_strdup (LOAD_PROC));
  list = g_list_append (list, g_strdup (SAVE_PROC));

  return list;
}

static GimpProcedure *
fits_create_procedure (GimpPlugIn  *plug_in,
                       const gchar *name)
{
  GimpProcedure *procedure = NULL;

  if (! strcmp (name, LOAD_PROC))
    {
      procedure = gimp_load_procedure_new (plug_in, name,
                                           GIMP_PDB_PROC_TYPE_PLUGIN,
                                           fits_load, NULL, NULL);

      gimp_procedure_set_menu_label (procedure,
                                     N_("Flexible Image Transport System"));

      gimp_procedure_set_documentation (procedure,
                                        _("Load file of the FITS file format"),
                                        _("Load file of the FITS file format "
                                          "(Flexible Image Transport System)"),
                                        name);
      gimp_procedure_set_attribution (procedure,
                                      "Peter Kirchgessner",
                                      "Peter Kirchgessner (peter@kirchgessner.net)",
                                      "1997");

      gimp_file_procedure_set_mime_types (GIMP_FILE_PROCEDURE (procedure),
                                          "image/x-fits");
      gimp_file_procedure_set_extensions (GIMP_FILE_PROCEDURE (procedure),
                                          "fit,fits");
      gimp_file_procedure_set_magics (GIMP_FILE_PROCEDURE (procedure),
                                      "0,string,SIMPLE");

      GIMP_PROC_AUX_ARG_INT (procedure, "replace",
                             _("Replacement for undefined pixels"),
                             _("Replacement for undefined pixels"),
                             0, 255, 0,
                             G_PARAM_READWRITE);

      GIMP_PROC_AUX_ARG_INT (procedure, "use-data-min-max",
                             _("Pixel value scaling"),
                             _("Use DATAMIN/DATAMAX-scaling if possible"),
                             0, 1, 0,
                             G_PARAM_READWRITE);
    }
  else if (! strcmp (name, SAVE_PROC))
    {
      procedure = gimp_save_procedure_new (plug_in, name,
                                           GIMP_PDB_PROC_TYPE_PLUGIN,
                                           fits_save, NULL, NULL);

      gimp_procedure_set_image_types (procedure, "RGB, GRAY, INDEXED");

      gimp_procedure_set_menu_label (procedure,
                                     N_("Flexible Image Transport System"));

      gimp_procedure_set_documentation (procedure,
                                        _("Export file in the FITS file format"),
                                        _("FITS exporting handles all image "
                                          "types except those with alpha channels."),
                                        name);
      gimp_procedure_set_attribution (procedure,
                                      "Peter Kirchgessner",
                                      "Peter Kirchgessner (peter@kirchgessner.net)",
                                      "1997");

      gimp_file_procedure_set_mime_types (GIMP_FILE_PROCEDURE (procedure),
                                          "image/x-fits");
      gimp_file_procedure_set_extensions (GIMP_FILE_PROCEDURE (procedure),
                                          "fit,fits");
    }

  return procedure;
}

static GimpValueArray *
fits_load (GimpProcedure        *procedure,
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

  if (run_mode == GIMP_RUN_INTERACTIVE)
    {
      if (! load_dialog (procedure, G_OBJECT (config)))
        return gimp_procedure_new_return_values (procedure,
                                                 GIMP_PDB_CANCEL,
                                                 NULL);
    }

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

static GimpValueArray *
fits_save (GimpProcedure        *procedure,
           GimpRunMode           run_mode,
           GimpImage            *image,
           gint                  n_drawables,
           GimpDrawable        **drawables,
           GFile                *file,
           const GimpValueArray *args,
           gpointer              run_data)
{
  GimpPDBStatusType  status = GIMP_PDB_SUCCESS;
  GimpExportReturn   export = GIMP_EXPORT_CANCEL;
  GError            *error = NULL;

  gegl_init (NULL, NULL);

  switch (run_mode)
    {
    case GIMP_RUN_INTERACTIVE:
    case GIMP_RUN_WITH_LAST_VALS:
      gimp_ui_init (PLUG_IN_BINARY);

      export = gimp_export_image (&image, &n_drawables, &drawables, "FITS",
                                  GIMP_EXPORT_CAN_HANDLE_RGB  |
                                  GIMP_EXPORT_CAN_HANDLE_GRAY |
                                  GIMP_EXPORT_CAN_HANDLE_INDEXED);

      if (export == GIMP_EXPORT_CANCEL)
        return gimp_procedure_new_return_values (procedure,
                                                 GIMP_PDB_CANCEL,
                                                 NULL);
      break;

    default:
      break;
    }

  if (n_drawables != 1)
    {
      g_set_error (&error, G_FILE_ERROR, 0,
                   _("FITS format does not support multiple layers."));

      return gimp_procedure_new_return_values (procedure,
                                               GIMP_PDB_CALLING_ERROR,
                                               error);
    }

  if (! save_image (file, image, drawables[0], &error))
    {
      status = GIMP_PDB_EXECUTION_ERROR;
    }

  if (export == GIMP_EXPORT_EXPORT)
    {
      gimp_image_delete (image);
      g_free (drawables);
    }

  return gimp_procedure_new_return_values (procedure, status, error);
}

static GimpImage *
load_image (GFile        *file,
            GObject      *config,
            GimpRunMode   run_mode,
            GError      **error)
{
  GimpImage         *image       = NULL;
  GimpLayer         *layer;
  GeglBuffer        *buffer;
  FILE              *fp;
  fitsfile          *ifp;
  FitsHduData        hdu;
  gint               n_pics;
  gint               count       = 1;
  gint               width;
  gint               height;
  gint               row_length;
  int                status      = 0;
  glong              fpixel[3]   = {1, 1, 1};
  GimpImageBaseType  itype       = GIMP_GRAY;
  GimpImageType      dtype       = GIMP_GRAYA_IMAGE;
  GimpPrecision      iprecision  = GIMP_PRECISION_U16_NON_LINEAR;
  const Babl        *type        = NULL;
  const Babl        *format      = NULL;
  gdouble           *pixels;
  gdouble            datamin     = 1.0E30f;
  gdouble            datamax     = -1.0E30f;
  gint               channels    = 1;
  gint               replace;
  gdouble            replace_val = 0;
  gboolean           use_datamin;

  g_object_get (config,
                "replace",          &replace,
                "use-data-min-max", &use_datamin,
                NULL);

  fp = g_fopen (g_file_peek_path (file), "rb");

  if (! fp)
    {
      g_set_error (error, G_FILE_ERROR, g_file_error_from_errno (errno),
                   _("Could not open '%s' for reading: %s"),
                   gimp_file_get_utf8_name (file), g_strerror (errno));
      return NULL;
    }

  fclose (fp);

  fits_open_diskfile (&ifp, g_file_peek_path (file), READONLY, &status);
  if (status)
    show_fits_errors (status);

  if (! ifp)
    {
      g_set_error (error, G_FILE_ERROR, G_FILE_ERROR_FAILED,
                   "%s", _("Error during open of FITS file"));
      return NULL;
    }

  /* Get first item */
  fits_get_num_hdus (ifp, &n_pics, &status);

  if (status)
    show_fits_errors (status);

  if (n_pics <= 0)
    {
      g_set_error (error, G_FILE_ERROR, G_FILE_ERROR_FAILED,
                   "%s", _("FITS file keeps no displayable images"));
      fits_close_file (ifp, &status);
      return NULL;
    }

  while (count <= n_pics)
    {
      hdu.naxis = 0;

      fits_movabs_hdu (ifp, count, NULL, &status);

      if (status)
        show_fits_errors (status);

      fits_get_img_param (ifp, 3, &hdu.bitpix, &hdu.naxis, hdu.naxisn,
                          &status);

      width  = hdu.naxisn[0];
      height = hdu.naxisn[1];

      if (status)
        show_fits_errors (status);

      /* Skip if invalid dimensions; possibly header data */
      if (hdu.naxis < 2)
        {
          count++;
          continue;
        }

      type = babl_type ("double");
      switch (hdu.bitpix)
      {
        case 8:
          iprecision = GIMP_PRECISION_U8_LINEAR;
          if (replace)
            replace_val = 255;
          break;
        case 16:
          iprecision = GIMP_PRECISION_U16_NON_LINEAR;
          if (replace)
            replace_val = G_MAXSHORT;
          break;
        case 32:
          iprecision = GIMP_PRECISION_U32_LINEAR;
          if (replace)
            replace_val = G_MAXINT;
          break;
        case -32:
          iprecision = GIMP_PRECISION_FLOAT_LINEAR;
          if (replace)
            replace_val = G_MAXFLOAT;
          break;
        case -64:
          iprecision = GIMP_PRECISION_DOUBLE_LINEAR;
          if (replace)
            replace_val = G_MAXDOUBLE;
          break;
      }

      if (hdu.naxis == 2)
        {
          itype = GIMP_GRAY;
          dtype = GIMP_GRAYA_IMAGE;
          format = babl_format_new (babl_model ("Y'"),
                                    type,
                                    babl_component ("Y'"),
                                    NULL);
        }
      else if (hdu.naxisn[2])
        {
          /* Original RGB format */
          if (hdu.naxisn[0] == 3)
            {
              width  = hdu.naxisn[1];
              height = hdu.naxisn[2];
            }
          channels = 3;

          itype  = GIMP_RGB;
          dtype  = GIMP_RGB_IMAGE;
          format = babl_format_new (babl_model ("R'G'B'"),
                                    type,
                                    babl_component ("R'"),
                                    babl_component ("G'"),
                                    babl_component ("B'"),
                                    NULL);
        }

      /* If RGB FITS image, we need to increase the size by the number of channels */
      pixels = (gdouble *) malloc (width * sizeof (gdouble) * channels);

      if (! image)
        {
          image = create_new_image (file, count, width, height,
                                    itype, dtype, iprecision,
                                    &layer, &buffer);
        }
      else
        {
          layer = gimp_layer_new (image, _("FITS HDU"), width, height,
                                  dtype, 100,
                                  gimp_image_get_default_new_layer_mode (image));
          gimp_image_insert_layer (image, layer, NULL, 0);
          buffer = gimp_drawable_get_buffer (GIMP_DRAWABLE (layer));
        }

      row_length = width * channels;

      /* Calculate min/max pixel value for normalizing */
      for (fpixel[1] = height; fpixel[1] >= 1; fpixel[1]--)
        {
          if (fits_read_pix (ifp, TDOUBLE, fpixel, row_length, NULL,
                             pixels, NULL, &status))
            break;

          for (gint ii = 0; ii < row_length; ii++)
            {
              if (pixels[ii] < datamin)
                datamin = pixels[ii];

              if (pixels[ii] > datamax)
                datamax = pixels[ii];
            }
        }

      if (status)
        show_fits_errors (status);

      /* Read pixel values in */
      for (fpixel[1] = height; fpixel[1] >= 1; fpixel[1]--)
        {
          gdouble *temp =
            (gdouble *) malloc (width * sizeof (gdouble) * channels);

          if (fits_read_pix (ifp, TDOUBLE, fpixel, row_length, &replace_val,
                             pixels, NULL, &status))
            break;

          if (datamin < datamax)
            {
              for (gint ii = 0; ii < row_length; ii++)
                pixels[ii] = (pixels[ii] - datamin) / (datamax - datamin);
            }

          if (hdu.naxisn[2] && hdu.naxisn[2] == 3) /* Packed RGB format */
            {
              /* Cover planes to RGB format */
              for (gint ii = 0; ii < (row_length / 3); ii++)
                {
                  temp[(ii * 3)]     = pixels[ii];
                  temp[(ii * 3) + 1] = pixels[ii + (row_length / 3)];
                  temp[(ii * 3) + 2] = pixels[ii + ((row_length / 3) * 2)];
                }
            }
          else
            {
              temp = pixels;
            }

          gegl_buffer_set (buffer,
                           GEGL_RECTANGLE (0, height - fpixel[1],
                                           width, 1), 0,
                           format, temp, GEGL_AUTO_ROWSTRIDE);
        }
      if (status)
        show_fits_errors (status);

      g_object_unref (buffer);
      if (pixels)
        g_free (pixels);

      count++;
    }

  /* As there might be different sized layers,
   * we need to resize the canvas afterwards */
  gimp_image_resize_to_layers (image);

  fits_close_file (ifp, &status);

  if (! image)
    {
      g_set_error (error, G_FILE_ERROR, G_FILE_ERROR_FAILED,
                   "%s", _("FITS file keeps no displayable images"));
      fits_close_file (ifp, &status);
      return NULL;
    }

  return image;
}

static gint
save_image (GFile         *file,
            GimpImage     *image,
            GimpDrawable  *drawable,
            GError       **error)
{
  FitsFile      *ofp;
  GimpImageType  drawable_type;
  gint           retval;

  drawable_type = gimp_drawable_type (drawable);

  /*  Make sure we're not exporting an image with an alpha channel  */
  if (gimp_drawable_has_alpha (drawable))
    {
      g_set_error (error, G_FILE_ERROR, G_FILE_ERROR_FAILED,
                   "%s",
                   _("FITS export cannot handle images with alpha channels"));
      return FALSE;
    }

  switch (drawable_type)
    {
    case GIMP_INDEXED_IMAGE: case GIMP_INDEXEDA_IMAGE:
    case GIMP_GRAY_IMAGE:    case GIMP_GRAYA_IMAGE:
    case GIMP_RGB_IMAGE:     case GIMP_RGBA_IMAGE:
      break;
    default:
      g_message (_("Cannot operate on unknown image types."));
      return (FALSE);
      break;
    }

  gimp_progress_init_printf (_("Exporting '%s'"),
                             gimp_file_get_utf8_name (file));

  /* Open the output file. */
  ofp = fits_open (g_file_peek_path (file), "w");

  if (! ofp)
    {
      g_set_error (error, G_FILE_ERROR, g_file_error_from_errno (errno),
                   _("Could not open '%s' for writing: %s"),
                   gimp_file_get_utf8_name (file), g_strerror (errno));
      return (FALSE);
    }

  retval = save_fits (ofp, image, drawable);

  fits_close (ofp);

  return retval;
}

/* Create an image. Sets layer_ID, drawable and rgn. Returns image_ID */
static GimpImage *
create_new_image (GFile              *file,
                  guint               pagenum,
                  guint               width,
                  guint               height,
                  GimpImageBaseType   itype,
                  GimpImageType       dtype,
                  GimpPrecision       iprecision,
                  GimpLayer         **layer,
                  GeglBuffer        **buffer)
{
  GimpImage *image;

  image = gimp_image_new_with_precision (width, height, itype, iprecision);

  gimp_image_undo_disable (image);
  *layer = gimp_layer_new (image, _("Background"), width, height,
                           dtype, 100,
                           gimp_image_get_default_new_layer_mode (image));
  gimp_image_insert_layer (image, *layer, NULL, 0);

  *buffer = gimp_drawable_get_buffer (GIMP_DRAWABLE (*layer));

  return image;
}

static FitsHduList *
create_fits_header (FitsFile *ofp,
                    guint     width,
                    guint     height,
                    guint     channels,
                    guint     bitpix)
{
  FitsHduList *hdulist;
  gint         print_ctype3 = 0; /* The CTYPE3-card may not be FITS-conforming */

  static const char *ctype3_card[] =
  {
    NULL, NULL, NULL,  /* bpp = 0: no additional card */
    "COMMENT Image type within GIMP: GIMP_GRAY_IMAGE",
    NULL,
    NULL,
    "COMMENT Image type within GIMP: GIMP_GRAYA_IMAGE (gray with alpha channel)",
    "COMMENT Sequence for NAXIS3   : GRAY, ALPHA",
    "CTYPE3  = 'GRAYA   '           / GRAY IMAGE WITH ALPHA CHANNEL",
    "COMMENT Image type within GIMP: GIMP_RGB_IMAGE",
    "COMMENT Sequence for NAXIS3   : RED, GREEN, BLUE",
    "CTYPE3  = 'RGB     '           / RGB IMAGE",
    "COMMENT Image type within GIMP: GIMP_RGBA_IMAGE (rgb with alpha channel)",
    "COMMENT Sequence for NAXIS3   : RED, GREEN, BLUE, ALPHA",
    "CTYPE3  = 'RGBA    '           / RGB IMAGE WITH ALPHA CHANNEL"
  };

  hdulist = fits_add_hdu (ofp);
  if (hdulist == NULL)
    return NULL;

  hdulist->used.simple  = 1;
  hdulist->bitpix       = bitpix;
  hdulist->naxis        = (channels == 1) ? 2 : 3;
  hdulist->naxisn[0]    = width;
  hdulist->naxisn[1]    = height;
  hdulist->naxisn[2]    = channels;
  hdulist->used.datamin = 1;
  hdulist->datamin      = 0.0;
  hdulist->used.datamax = 1;
  hdulist->used.bzero   = 1;
  hdulist->bzero        = 0.0;
  hdulist->used.bscale  = 1;
  hdulist->bscale       = 1.0;

  switch (bitpix)
    {
    case 8:
      hdulist->datamax = 255;
      break;
    case 16:
      hdulist->datamax = 65535;
      break;
    case 32:
      hdulist->datamax = 4294967295.0; /* .0 to silence gcc */
      break;
    case -32:
      hdulist->datamax = 1.0;
      break;
    case -64:
      hdulist->datamax = 1.0;
      break;
    default:
      return NULL;
    }

  fits_add_card (hdulist, "");
  fits_add_card (hdulist,
                 "HISTORY THIS FITS FILE WAS GENERATED BY GIMP USING FITSRW");
  fits_add_card (hdulist, "");
  fits_add_card (hdulist,
                 "COMMENT FitsRW is (C) Peter Kirchgessner (peter@kirchgessner.net), but available");
  fits_add_card (hdulist,
                 "COMMENT under the GNU general public licence.");
  fits_add_card (hdulist,
                 "COMMENT For sources see http://www.kirchgessner.net");
  fits_add_card (hdulist, "");
  fits_add_card (hdulist, ctype3_card[channels * 3]);

  if (ctype3_card[channels * 3 + 1] != NULL)
    fits_add_card (hdulist, ctype3_card[channels * 3 + 1]);

  if (print_ctype3 && (ctype3_card[channels * 3 + 2] != NULL))
    fits_add_card (hdulist, ctype3_card[channels * 3 + 2]);

  fits_add_card (hdulist, "");

  return hdulist;
}


/* Save direct colors (GRAY, GRAYA, RGB, RGBA) */
static gint
save_fits (FitsFile     *ofp,
           GimpImage    *image,
           GimpDrawable *drawable)
{
  gint           height, width, i, j, channel, channelnum;
  gint           tile_height, bpp, bpsl, bitpix, bpc;
  long           nbytes;
  guchar        *data, *src;
  GeglBuffer    *buffer;
  const Babl    *format, *type;
  FitsHduList   *hdu;

  buffer = gimp_drawable_get_buffer (drawable);

  width  = gegl_buffer_get_width  (buffer);
  height = gegl_buffer_get_height (buffer);

  format = gegl_buffer_get_format (buffer);
  type   = babl_format_get_type (format, 0);

  if (type == babl_type ("u8"))
    {
      bitpix = 8;
    }
  else if (type == babl_type ("u16"))
    {
      bitpix = 16;
    }
  else if (type == babl_type ("u32"))
    {
      bitpix = 32;
    }
  else if (type == babl_type ("half"))
    {
      bitpix = -32;
      type = babl_type ("float");
    }
  else if (type == babl_type ("float"))
    {
      bitpix = -32;
    }
  else if (type == babl_type ("double"))
    {
      bitpix = -64;
    }
  else
    {
      return FALSE;
    }

  switch (gimp_drawable_type (drawable))
    {
    case GIMP_GRAY_IMAGE:
      format = babl_format_new (babl_model ("Y'"),
                                type,
                                babl_component ("Y'"),
                                NULL);
      break;

    case GIMP_GRAYA_IMAGE:
      format = babl_format_new (babl_model ("Y'A"),
                                type,
                                babl_component ("Y'"),
                                babl_component ("A"),
                                NULL);
      break;

    case GIMP_RGB_IMAGE:
    case GIMP_INDEXED_IMAGE:
      format = babl_format_new (babl_model ("R'G'B'"),
                                type,
                                babl_component ("R'"),
                                babl_component ("G'"),
                                babl_component ("B'"),
                                NULL);
      break;

    case GIMP_RGBA_IMAGE:
    case GIMP_INDEXEDA_IMAGE:
      format = babl_format_new (babl_model ("R'G'B'A"),
                                type,
                                babl_component ("R'"),
                                babl_component ("G'"),
                                babl_component ("B'"),
                                babl_component ("A"),
                                NULL);
      break;
    }

  channelnum = babl_format_get_n_components (format);
  bpp        = babl_format_get_bytes_per_pixel (format);

  bpc  = bpp / channelnum; /* Bytes per channel */
  bpsl = width * bpp;      /* Bytes per scanline */

  tile_height = gimp_tile_height ();

  /* allocate a buffer for retrieving information from the pixel region  */
  src = data = (guchar *) g_malloc (width * height * bpp);

  hdu = create_fits_header (ofp, width, height, channelnum, bitpix);
  if (hdu == NULL)
    return FALSE;

  if (fits_write_header (ofp, hdu) < 0)
    return FALSE;

  nbytes = 0;
  for (channel = 0; channel < channelnum; channel++)
    {
      for (i = 0; i < height; i++)
        {
          if ((i % tile_height) == 0)
            {
              gint scan_lines;

              scan_lines = (i + tile_height-1 < height) ?
                            tile_height : (height - i);

              gegl_buffer_get (buffer,
                               GEGL_RECTANGLE (0, height - i - scan_lines,
                                               width, scan_lines), 1.0,
                               format, data,
                               GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);

              src = data + bpsl * (scan_lines - 1) + channel * bpc;
            }

          if (channelnum == 1 && bitpix == 8)  /* One channel and 8 bit? Write the scanline */
            {
              fwrite (src, bpc, width, ofp->fp);
              src += bpsl;
            }
          else           /* Multiple channels or high bit depth */
            {
            /* Write out bytes for current channel */
            /* FIXME: Don't assume a little endian arch */
            switch (bitpix)
              {
              case 8:
                for (j = 0; j < width; j++)
                  {
                    putc (*src, ofp->fp);
                    src += bpp;
                  }
                break;
              case 16:
                for (j = 0; j < width; j++)
                  {
                    *((guint16*)src) += 32768;
                    putc (*(src + 1), ofp->fp);
                    putc (*(src + 0), ofp->fp);
                    src += bpp;
                  }
                break;
              case 32:
                for (j = 0; j < width; j++)
                  {
                    *((guint32*)src) += 2147483648.0; /* .0 to silence gcc */
                    putc (*(src + 3), ofp->fp);
                    putc (*(src + 2), ofp->fp);
                    putc (*(src + 1), ofp->fp);
                    putc (*(src + 0), ofp->fp);
                    src += bpp;
                  }
                break;
              case -32:
                for (j = 0; j < width; j++)
                  {
                    putc (*(src + 3), ofp->fp);
                    putc (*(src + 2), ofp->fp);
                    putc (*(src + 1), ofp->fp);
                    putc (*(src + 0), ofp->fp);
                    src += bpp;
                  }
                break;
              case -64:
                for (j = 0; j < width; j++)
                  {
                    putc (*(src + 7), ofp->fp);
                    putc (*(src + 6), ofp->fp);
                    putc (*(src + 5), ofp->fp);
                    putc (*(src + 4), ofp->fp);
                    putc (*(src + 3), ofp->fp);
                    putc (*(src + 2), ofp->fp);
                    putc (*(src + 1), ofp->fp);
                    putc (*(src + 0), ofp->fp);
                    src += bpp;
                  }
                break;
              default:
                return FALSE;
              }
            }

          nbytes += width * bpc;
          src -= 2 * bpsl;

          if ((i % 20) == 0)
            gimp_progress_update ((gdouble) (i + channel * height) /
                                  (gdouble) (height * channelnum));
        }
    }

  nbytes = nbytes % FITS_RECORD_SIZE;
  if (nbytes)
    {
      while (nbytes++ < FITS_RECORD_SIZE)
        putc (0, ofp->fp);
    }

  g_free (data);

  g_object_unref (buffer);

  gimp_progress_update (1.0);

  if (ferror (ofp->fp))
    {
      g_message (_("Write error occurred"));
      return FALSE;
    }

  return TRUE;
}


/*  Load interface functions  */

static gboolean
load_dialog (GimpProcedure *procedure,
             GObject       *config)
{
  GtkWidget    *dialog;
  GtkListStore *store;
  GtkWidget    *frame;
  gboolean      run;

  gimp_ui_init (PLUG_IN_BINARY);

  dialog = gimp_procedure_dialog_new (procedure,
                                      GIMP_PROCEDURE_CONFIG (config),
                                      _("Open FITS File"));

  gimp_window_set_transient (GTK_WINDOW (dialog));

  store = gimp_int_store_new (_("Black"), 0,
                              _("White"), 255,
                              NULL);
  frame = gimp_procedure_dialog_get_int_radio (GIMP_PROCEDURE_DIALOG (dialog),
                                               "replace",
                                               GIMP_INT_STORE (store));
  gtk_widget_set_margin_bottom (frame, 12);

  store = gimp_int_store_new (_("Automatic"),          0,
                              _("By DATAMIN/DATAMAX"), 1,
                              NULL);
  frame = gimp_procedure_dialog_get_int_radio (GIMP_PROCEDURE_DIALOG (dialog),
                                               "use-data-min-max",
                                               GIMP_INT_STORE (store));
  gtk_widget_set_margin_bottom (frame, 12);

  gimp_procedure_dialog_fill (GIMP_PROCEDURE_DIALOG (dialog),
                              NULL);

  gtk_widget_show (dialog);

  run = gimp_procedure_dialog_run (GIMP_PROCEDURE_DIALOG (dialog));

  gtk_widget_destroy (dialog);

  return run;
}

static void
show_fits_errors (gint status)
{
  gchar status_str[FLEN_STATUS];

  /* Write out error messages of FITS-Library */
  fits_get_errstatus (status, status_str);
  g_message ("FITS: %s\n", status_str);
}
