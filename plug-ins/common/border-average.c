/* borderaverage 0.01 - image processing plug-in for GIMP.
 *
 * Copyright (C) 1998 Philipp Klaus (webmaster@access.ch)
 *
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

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>

#include "libgimp/stdplugins-intl.h"


#define PLUG_IN_PROC   "plug-in-borderaverage"
#define PLUG_IN_BINARY "border-average"
#define PLUG_IN_ROLE   "gimp-border-average"


typedef struct _BorderAverage      BorderAverage;
typedef struct _BorderAverageClass BorderAverageClass;

struct _BorderAverage
{
  GimpPlugIn      parent_instance;
};

struct _BorderAverageClass
{
  GimpPlugInClass parent_class;
};


#define BORDER_AVERAGE_TYPE  (border_average_get_type ())
#define BORDER_AVERAGE (obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), BORDER_AVERAGE_TYPE, BorderAverage))


GType                   border_average_get_type               (void) G_GNUC_CONST;

static GList          * border_average_query_procedures       (GimpPlugIn           *plug_in);
static GimpProcedure  * border_average_create_procedure       (GimpPlugIn           *plug_in,
                                                               const gchar          *name);

static GimpValueArray * border_average_run                    (GimpProcedure        *procedure,
                                                               GimpRunMode           run_mode,
                                                               GimpImage            *image,
                                                               gint                  n_drawables,
                                                               GimpDrawable        **drawables,
                                                               const GimpValueArray *args,
                                                               gpointer              run_data);


static void      borderaverage        (GObject      *config,
                                       GeglBuffer   *buffer,
                                       GimpDrawable *drawable,
                                       GimpRGB      *result);

static gboolean  borderaverage_dialog (GimpProcedure *procedure,
                                       GObject       *config,
                                       GimpImage     *image,
                                       GimpDrawable  *drawable);

static void      add_new_color        (const guchar *buffer,
                                       gint         *cube,
                                       gint          bucket_expo);


G_DEFINE_TYPE (BorderAverage, border_average, GIMP_TYPE_PLUG_IN)

GIMP_MAIN (BORDER_AVERAGE_TYPE)
DEFINE_STD_SET_I18N


static void
border_average_class_init (BorderAverageClass *klass)
{
  GimpPlugInClass *plug_in_class = GIMP_PLUG_IN_CLASS (klass);

  plug_in_class->query_procedures = border_average_query_procedures;
  plug_in_class->create_procedure = border_average_create_procedure;
  plug_in_class->set_i18n         = STD_SET_I18N;
}

static void
border_average_init (BorderAverage *film)
{
}

static GList *
border_average_query_procedures (GimpPlugIn *plug_in)
{
  return g_list_append (NULL, g_strdup (PLUG_IN_PROC));
}

static GimpProcedure *
border_average_create_procedure (GimpPlugIn  *plug_in,
                               const gchar *name)
{
  GimpProcedure *procedure = NULL;

  if (! strcmp (name, PLUG_IN_PROC))
    {
      procedure = gimp_image_procedure_new (plug_in, name,
                                            GIMP_PDB_PROC_TYPE_PLUGIN,
                                            border_average_run, NULL, NULL);

      gimp_procedure_set_image_types (procedure, "RGB*");
      gimp_procedure_set_sensitivity_mask (procedure,
                                           GIMP_PROCEDURE_SENSITIVE_DRAWABLE);

      gimp_procedure_set_menu_label (procedure, _("_Border Average..."));
      gimp_procedure_add_menu_path (procedure, "<Image>/Colors/Info");

      gimp_procedure_set_documentation (procedure,
                                        _("Set foreground to the average color of the image border"),
                                        "",
                                        name);
      gimp_procedure_set_attribution (procedure,
                                      "Philipp Klaus",
                                      "Internet Access AG",
                                      "1998");

      GIMP_PROC_ARG_INT (procedure, "thickness",
                         _("_Thickness"),
                         _("Border size to take in count"),
                         0, G_MAXINT, 3,
                         G_PARAM_READWRITE);
      GIMP_PROC_AUX_ARG_UNIT (procedure, "thickness-unit",
                              _("Thickness unit of measure"),
                              _("Border size unit of measure"),
                              TRUE, TRUE, GIMP_UNIT_PIXEL,
                              GIMP_PARAM_READWRITE);
      GIMP_PROC_ARG_INT (procedure, "bucket-exponent",
                         _("Bucket Si_ze"),
                         _("Bits for bucket size (default=4: 16 Levels)"),
                         0, G_MAXINT, 4,
                         G_PARAM_READWRITE);

      GIMP_PROC_VAL_RGB (procedure, "borderaverage",
                         _("The average color of the specified border."),
                         _("The average color of the specified border."),
                         TRUE, NULL,
                         G_PARAM_READWRITE);
    }

  return procedure;
}


static GimpValueArray *
border_average_run (GimpProcedure        *procedure,
                    GimpRunMode           run_mode,
                    GimpImage            *image,
                    gint                  n_drawables,
                    GimpDrawable        **drawables,
                    const GimpValueArray *args,
                    gpointer              run_data)
{
  GimpProcedureConfig *config;
  GimpDrawable        *drawable;
  GimpValueArray      *return_vals = NULL;
  GimpPDBStatusType    status = GIMP_PDB_SUCCESS;
  GimpRGB              result_color = { 0.0, };
  GeglBuffer          *buffer;

  gegl_init (NULL, NULL);

  if (n_drawables != 1)
    {
      GError *error = NULL;

      g_set_error (&error, GIMP_PLUG_IN_ERROR, 0,
                   _("Procedure '%s' only works with one drawable."),
                   PLUG_IN_PROC);

      return gimp_procedure_new_return_values (procedure,
                                               GIMP_PDB_CALLING_ERROR,
                                               error);
    }
  else
    {
      drawable = drawables[0];
    }

  config = gimp_procedure_create_config (procedure);
  gimp_procedure_config_begin_run (config, NULL, run_mode, args);

  buffer = gimp_drawable_get_buffer (drawable);

  switch (run_mode)
    {
    case GIMP_RUN_INTERACTIVE:
      if (! borderaverage_dialog (procedure, G_OBJECT (config),
                                  image, drawable))
        status = GIMP_PDB_EXECUTION_ERROR;
      break;

    case GIMP_RUN_NONINTERACTIVE:
      if (gimp_value_array_length (args) != 2)
        status = GIMP_PDB_CALLING_ERROR;
      break;

    default:
      break;
    }

  if (status == GIMP_PDB_SUCCESS)
    {
      /*  Make sure that the drawable is RGB color  */
      if (gimp_drawable_is_rgb (drawable))
        {
          gimp_progress_init ( _("Border Average"));
          borderaverage (G_OBJECT (config), buffer, drawable, &result_color);

          if (run_mode != GIMP_RUN_NONINTERACTIVE)
            gimp_context_set_foreground (&result_color);
        }
      else
        {
          status = GIMP_PDB_EXECUTION_ERROR;
        }
    }

  g_object_unref (buffer);

  return_vals = gimp_procedure_new_return_values (procedure, status, NULL);

  gimp_procedure_config_end_run (config, GIMP_PDB_SUCCESS);
  g_object_unref (config);

  if (status == GIMP_PDB_SUCCESS)
    GIMP_VALUES_SET_RGB (return_vals, 1, &result_color);

  return return_vals;
}


static void
borderaverage (GObject      *config,
               GeglBuffer   *buffer,
               GimpDrawable *drawable,
               GimpRGB      *result)
{
  gint            x, y, width, height;
  gint            max;
  guchar          r, g, b;
  gint            bucket_num, bucket_expo, bucket_rexpo;
  gint           *cube;
  gint            i, j, k;
  GeglRectangle   border[4];
  gint            borderaverage_thickness;
  gint            borderaverage_bucket_exponent;

  g_object_get (config,
                "thickness", &borderaverage_thickness,
                "bucket-exponent", &borderaverage_bucket_exponent,
                NULL);

  if (! gimp_drawable_mask_intersect (drawable, &x, &y, &width, &height))
    {
      gimp_rgba_set_uchar (result, 0, 0, 0, 255);
      return;
    }

  /* allocate and clear the cube before */
  bucket_expo = borderaverage_bucket_exponent;
  bucket_rexpo = 8 - bucket_expo;
  cube = g_new (gint, 1 << (bucket_rexpo * 3));
  bucket_num = 1 << bucket_rexpo;

  for (i = 0; i < bucket_num; i++)
    {
      for (j = 0; j < bucket_num; j++)
        {
          for (k = 0; k < bucket_num; k++)
            {
              cube[(i << (bucket_rexpo << 1)) + (j << bucket_rexpo) + k] = 0;
            }
        }
    }

  /* Top */
  border[0].x =       x;
  border[0].y =       y;
  border[0].width =   width;
  border[0].height =  borderaverage_thickness;

  /* Bottom */
  border[1].x =       x;
  border[1].y =       y + height - borderaverage_thickness;
  border[1].width =   width;
  border[1].height =  borderaverage_thickness;

  /* Left */
  border[2].x =       x;
  border[2].y =       y + borderaverage_thickness;
  border[2].width =   borderaverage_thickness;
  border[2].height =  height - 2 * borderaverage_thickness;

  /* Right */
  border[3].x =       x + width - borderaverage_thickness;
  border[3].y =       y + borderaverage_thickness;
  border[3].width =   borderaverage_thickness;
  border[3].height =  height - 2 * borderaverage_thickness;

  /* Fill the cube */
  for (i = 0; i < 4; i++)
    {
      if (border[i].width > 0 && border[i].height > 0)
        {
          GeglBufferIterator *gi;

          gi = gegl_buffer_iterator_new (buffer, &border[i], 0, babl_format ("R'G'B' u8"),
                                         GEGL_ACCESS_READWRITE, GEGL_ABYSS_NONE, 1);

          while (gegl_buffer_iterator_next (gi))
            {
              guint   k;
              guchar *data;

              data = (guchar*) gi->items[0].data;

              for (k = 0; k < gi->length; k++)
                {
                  add_new_color (data + k * 3,
                                 cube,
                                 bucket_expo);
                }
            }
        }
    }

  max = 0; r = 0; g = 0; b = 0;

  /* get max of cube */
  for (i = 0; i < bucket_num; i++)
    {
      for (j = 0; j < bucket_num; j++)
        {
          for (k = 0; k < bucket_num; k++)
            {
              if (cube[(i << (bucket_rexpo << 1)) +
                      (j << bucket_rexpo) + k] > max)
                {
                  max = cube[(i << (bucket_rexpo << 1)) +
                             (j << bucket_rexpo) + k];
                  r = (i<<bucket_expo) + (1<<(bucket_expo - 1));
                  g = (j<<bucket_expo) + (1<<(bucket_expo - 1));
                  b = (k<<bucket_expo) + (1<<(bucket_expo - 1));
                }
            }
        }
    }

  /* return the color */
  gimp_rgba_set_uchar (result, r, g, b, 255);

  g_free (cube);
}

static void
add_new_color (const guchar *buffer,
               gint         *cube,
               gint          bucket_expo)
{
  guchar r, g, b;
  gint   bucket_rexpo;

  bucket_rexpo = 8 - bucket_expo;
  r = buffer[0] >> bucket_expo;
  g = buffer[1] >> bucket_expo;
  b = buffer[2] >> bucket_expo;
  cube[(r << (bucket_rexpo << 1)) + (g << bucket_rexpo) + b]++;
}

static gboolean
borderaverage_dialog (GimpProcedure *procedure,
                      GObject       *config,
                      GimpImage     *image,
                      GimpDrawable  *drawable)
{
  GtkWidget    *dialog;
  GtkListStore *store;
  GtkWidget    *size_entry;
  gboolean      run;
  gdouble       xres, yres;
  GeglBuffer   *buffer = NULL;

  gimp_ui_init (PLUG_IN_BINARY);

  dialog = gimp_procedure_dialog_new (procedure,
                                      GIMP_PROCEDURE_CONFIG (config),
                                      _("Border Average"));

  gimp_dialog_set_alternative_button_order (GTK_DIALOG (dialog),
                                           GTK_RESPONSE_OK,
                                           GTK_RESPONSE_CANCEL,
                                           -1);

  gimp_window_set_transient (GTK_WINDOW (dialog));

  gimp_procedure_dialog_get_label (GIMP_PROCEDURE_DIALOG (dialog),
                                   "border-size-label",
                                   _("Border Size"));

  /*  Get the image resolution  */
  gimp_image_get_resolution (image, &xres, &yres);
  size_entry = gimp_procedure_dialog_get_size_entry (GIMP_PROCEDURE_DIALOG (dialog),
                                                     "thickness", TRUE,
                                                     "thickness-unit", "%a",
                                                     GIMP_SIZE_ENTRY_UPDATE_SIZE,
                                                     xres);

  /*  set the size (in pixels) that will be treated as 0% and 100%  */
  buffer = gimp_drawable_get_buffer (drawable);
  if (buffer)
    gimp_size_entry_set_size (GIMP_SIZE_ENTRY (size_entry), 0, 0.0,
                              MIN (gegl_buffer_get_width (buffer),
                                   gegl_buffer_get_height (buffer)));
  gimp_size_entry_set_refval_boundaries (GIMP_SIZE_ENTRY (size_entry), 0,
                                         1.0, 256.0);

  gimp_procedure_dialog_fill_frame (GIMP_PROCEDURE_DIALOG (dialog),
                                    "border-size-frame", "border-size-label",
                                    FALSE, "thickness");

  gimp_procedure_dialog_get_label (GIMP_PROCEDURE_DIALOG (dialog),
                                   "bucket-size-label",
                                   _("Number of Colors"));

  store = gimp_int_store_new ("1",   0, "2",   1, "4",   2, "8",   3,
                              "16",  4, "32",  5, "64",  6, "128", 7,
                              "256", 8, NULL);
  gimp_procedure_dialog_get_int_combo (GIMP_PROCEDURE_DIALOG (dialog),
                                       "bucket-exponent",
                                       GIMP_INT_STORE (store));

  gimp_procedure_dialog_fill_frame (GIMP_PROCEDURE_DIALOG (dialog),
                                    "bucket-size-frame", "bucket-size-label",
                                    FALSE, "bucket-exponent");

  gimp_procedure_dialog_fill (GIMP_PROCEDURE_DIALOG (dialog),
                              "border-size-frame",
                              "bucket-size-frame",
                              NULL);

  gtk_widget_show (dialog);

  run = gimp_procedure_dialog_run (GIMP_PROCEDURE_DIALOG (dialog));

  gtk_widget_destroy (dialog);

  return run;
}
