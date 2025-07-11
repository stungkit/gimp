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

#include <gegl.h>
#include <gtk/gtk.h>

#include "libgimpbase/gimpbase.h"
#include "libgimpwidgets/gimpwidgets.h"

#include "actions-types.h"

#include "core/gimp.h"
#include "core/gimp-edit.h"
#include "core/gimpbuffer.h"
#include "core/gimpcontainer.h"
#include "core/gimpdrawable.h"
#include "core/gimpdrawable-edit.h"
#include "core/gimpfilloptions.h"
#include "core/gimplayer.h"
#include "core/gimplayer-new.h"
#include "core/gimplayermask.h"
#include "core/gimplist.h"
#include "core/gimpimage.h"
#include "core/gimpimage-undo.h"

#include "path/gimppath-import.h"

#include "widgets/gimpclipboard.h"
#include "widgets/gimphelp-ids.h"
#include "widgets/gimpdialogfactory.h"
#include "widgets/gimpmessagebox.h"
#include "widgets/gimpmessagedialog.h"
#include "widgets/gimpwidgets-utils.h"
#include "widgets/gimpwindowstrategy.h"

#include "display/gimpdisplay.h"
#include "display/gimpdisplayshell.h"
#include "display/gimpdisplayshell-transform.h"

#include "tools/gimptools-utils.h"
#include "tools/tool_manager.h"

#include "actions.h"
#include "edit-commands.h"

#include "gimp-intl.h"


/*  local function prototypes  */

static gboolean   check_drawable_alpha               (GimpDrawable  *drawable,
                                                      gpointer       data);
static void       edit_paste                         (GimpDisplay   *display,
                                                      GimpPasteType  paste_type,
                                                      gboolean       merged,
                                                      gboolean       try_svg);
static void       cut_named_buffer_callback          (GtkWidget     *widget,
                                                      const gchar   *name,
                                                      gpointer       data);
static void       copy_named_buffer_callback         (GtkWidget     *widget,
                                                      const gchar   *name,
                                                      gpointer       data);
static void       copy_named_visible_buffer_callback (GtkWidget     *widget,
                                                      const gchar   *name,
                                                      gpointer       data);


/*  public functions  */

void
edit_undo_cmd_callback (GimpAction *action,
                        GVariant   *value,
                        gpointer    data)
{
  GimpImage   *image;
  GimpDisplay *display;
  return_if_no_image (image, data);
  return_if_no_display (display, data);

  if (tool_manager_undo_active (image->gimp, display) ||
      gimp_image_undo (image))
    {
      gimp_image_flush (image);
    }
}

void
edit_redo_cmd_callback (GimpAction *action,
                        GVariant   *value,
                        gpointer    data)
{
  GimpImage   *image;
  GimpDisplay *display;
  return_if_no_image (image, data);
  return_if_no_display (display, data);

  if (tool_manager_redo_active (image->gimp, display) ||
      gimp_image_redo (image))
    {
      gimp_image_flush (image);
    }
}

void
edit_strong_undo_cmd_callback (GimpAction *action,
                               GVariant   *value,
                               gpointer    data)
{
  GimpImage *image;
  return_if_no_image (image, data);

  if (gimp_image_strong_undo (image))
    gimp_image_flush (image);
}

void
edit_strong_redo_cmd_callback (GimpAction *action,
                               GVariant   *value,
                               gpointer    data)
{
  GimpImage *image;
  return_if_no_image (image, data);

  if (gimp_image_strong_redo (image))
    gimp_image_flush (image);
}

void
edit_undo_clear_cmd_callback (GimpAction *action,
                              GVariant   *value,
                              gpointer    data)
{
  GimpImage     *image;
  GimpUndoStack *undo_stack;
  GimpUndoStack *redo_stack;
  GtkWidget     *widget;
  GtkWidget     *dialog;
  gchar         *size;
  gint64         memsize;
  gint64         guisize;
  return_if_no_image (image, data);
  return_if_no_widget (widget, data);

  dialog = gimp_message_dialog_new (_("Clear Undo History"),
                                    GIMP_ICON_DIALOG_WARNING,
                                    widget,
                                    GTK_DIALOG_MODAL |
                                    GTK_DIALOG_DESTROY_WITH_PARENT,
                                    gimp_standard_help_func,
                                    GIMP_HELP_EDIT_UNDO_CLEAR,

                                    _("_Cancel"), GTK_RESPONSE_CANCEL,
                                    _("Cl_ear"),  GTK_RESPONSE_OK,

                                    NULL);

  gimp_dialog_set_alternative_button_order (GTK_DIALOG (dialog),
                                           GTK_RESPONSE_OK,
                                           GTK_RESPONSE_CANCEL,
                                           -1);

  g_signal_connect_object (gtk_widget_get_toplevel (widget), "unmap",
                           G_CALLBACK (gtk_widget_destroy),
                           dialog, G_CONNECT_SWAPPED);

  g_signal_connect_object (image, "disconnect",
                           G_CALLBACK (gtk_widget_destroy),
                           dialog, G_CONNECT_SWAPPED);

  gimp_message_box_set_primary_text (GIMP_MESSAGE_DIALOG (dialog)->box,
                                     _("Really clear image's undo history?"));

  undo_stack = gimp_image_get_undo_stack (image);
  redo_stack = gimp_image_get_redo_stack (image);

  memsize =  gimp_object_get_memsize (GIMP_OBJECT (undo_stack), &guisize);
  memsize += guisize;
  memsize += gimp_object_get_memsize (GIMP_OBJECT (redo_stack), &guisize);
  memsize += guisize;

  size = g_format_size (memsize);

  gimp_message_box_set_text (GIMP_MESSAGE_DIALOG (dialog)->box,
                             _("Clearing the undo history of this "
                               "image will gain %s of memory."), size);
  g_free (size);

  if (gimp_dialog_run (GIMP_DIALOG (dialog)) == GTK_RESPONSE_OK)
    {
      gimp_image_undo_disable (image);
      gimp_image_undo_enable (image);
      gimp_image_flush (image);
    }

  gtk_widget_destroy (dialog);
}

void
edit_cut_cmd_callback (GimpAction *action,
                       GVariant   *value,
                       gpointer    data)
{
  GimpImage    *image;
  GList        *drawables;
  GList        *iter;
  GimpObject   *cut;
  GError       *error = NULL;
  return_if_no_drawables (image, drawables, data);

  for (iter = drawables; iter; iter = iter->next)
    if (! check_drawable_alpha (iter->data, data))
      {
        g_list_free (drawables);
        return;
      }

  cut = gimp_edit_cut (image, drawables, action_data_get_context (data), &error);

  if (cut)
    {
      GimpDisplay *display = action_data_get_display (data);

      if (display)
        {
          gchar *msg;

          if (GIMP_IS_IMAGE (cut))
            msg = g_strdup_printf (ngettext ("Cut layer to the clipboard.",
                                             "Cut %d layers to the clipboard.",
                                             g_list_length (drawables)),
                                   g_list_length (drawables));
          else
            msg = g_strdup (_("Cut pixels to the clipboard."));

          gimp_message_literal (image->gimp,
                                G_OBJECT (display), GIMP_MESSAGE_INFO,
                                msg);
          g_free (msg);
        }

      gimp_image_flush (image);
    }
  else
    {
      gimp_message_literal (image->gimp,
                            G_OBJECT (action_data_get_display (data)),
                            GIMP_MESSAGE_WARNING,
                            error->message);
      g_clear_error (&error);
    }
  g_list_free (drawables);
}

void
edit_copy_cmd_callback (GimpAction *action,
                        GVariant   *value,
                        gpointer    data)
{
  GimpImage    *image;
  GList        *drawables;
  GimpObject   *copy;
  GError       *error = NULL;
  return_if_no_drawables (image, drawables, data);

  copy = gimp_edit_copy (image, drawables, action_data_get_context (data),
                         FALSE, &error);

  if (copy)
    {
      GimpDisplay *display = action_data_get_display (data);

      if (display)
        gimp_message_literal (image->gimp,
                              G_OBJECT (display), GIMP_MESSAGE_INFO,
                              GIMP_IS_IMAGE (copy) ?
                              _("Copied layer to the clipboard.") :
                              _("Copied pixels to the clipboard."));

      gimp_image_flush (image);
    }
  else
    {
      gimp_message_literal (image->gimp,
                            G_OBJECT (action_data_get_display (data)),
                            GIMP_MESSAGE_WARNING,
                            error->message);
      g_clear_error (&error);
    }

  g_list_free (drawables);
}

void
edit_copy_visible_cmd_callback (GimpAction *action,
                                GVariant   *value,
                                gpointer    data)
{
  GimpImage *image;
  GError    *error = NULL;
  return_if_no_image (image, data);

  if (gimp_edit_copy_visible (image, action_data_get_context (data), &error))
    {
      GimpDisplay *display = action_data_get_display (data);

      if (display)
        gimp_message_literal (image->gimp,
                              G_OBJECT (display), GIMP_MESSAGE_INFO,
                              _("Copied pixels to the clipboard."));

      gimp_image_flush (image);
    }
  else
    {
      gimp_message_literal (image->gimp,
                            G_OBJECT (action_data_get_display (data)),
                            GIMP_MESSAGE_WARNING,
                            error->message);
      g_clear_error (&error);
    }
}

void
edit_paste_cmd_callback (GimpAction *action,
                         GVariant   *value,
                         gpointer    data)
{
  GimpImage     *image;
  GimpDisplay   *display    = action_data_get_display (data);
  GimpPasteType  paste_type = (GimpPasteType) g_variant_get_int32 (value);
  GimpPasteType  converted_type;
  GList         *drawables;
  gboolean       merged     = FALSE;

  if (paste_type == GIMP_PASTE_TYPE_NEW_LAYER_OR_FLOATING)
    {
      if (! display || ! gimp_display_get_image (display))
        {
          edit_paste_as_new_image_cmd_callback (action, value, data);
          return;
        }
    }

  return_if_no_image (image, data);

  if (! display)
    return;

  switch (paste_type)
    {
    case GIMP_PASTE_TYPE_FLOATING:
    case GIMP_PASTE_TYPE_FLOATING_IN_PLACE:
    case GIMP_PASTE_TYPE_FLOATING_INTO:
    case GIMP_PASTE_TYPE_FLOATING_INTO_IN_PLACE:
      edit_paste (display, paste_type, merged, TRUE);
      break;

    case GIMP_PASTE_TYPE_NEW_LAYER:
    case GIMP_PASTE_TYPE_NEW_LAYER_IN_PLACE:
      edit_paste (display, paste_type, merged, FALSE);
      break;

    case GIMP_PASTE_TYPE_NEW_MERGED_LAYER_OR_FLOATING:
    case GIMP_PASTE_TYPE_NEW_MERGED_LAYER_OR_FLOATING_IN_PLACE:
      merged = TRUE;
    case GIMP_PASTE_TYPE_NEW_LAYER_OR_FLOATING:
    case GIMP_PASTE_TYPE_NEW_LAYER_OR_FLOATING_IN_PLACE:
      drawables = gimp_image_get_selected_drawables (image);

      if (drawables &&
         (g_list_length (drawables) == 1) &&
          GIMP_IS_LAYER_MASK (drawables->data))
        {
          converted_type = (paste_type == GIMP_PASTE_TYPE_NEW_LAYER_OR_FLOATING ||
                            paste_type == GIMP_PASTE_TYPE_NEW_MERGED_LAYER_OR_FLOATING) ?
                            GIMP_PASTE_TYPE_FLOATING :
                            GIMP_PASTE_TYPE_FLOATING_IN_PLACE;

          edit_paste (display, converted_type, merged, TRUE);
        }
      else
        {
          converted_type = (paste_type == GIMP_PASTE_TYPE_NEW_LAYER_OR_FLOATING ||
                            paste_type == GIMP_PASTE_TYPE_NEW_MERGED_LAYER_OR_FLOATING) ?
                            GIMP_PASTE_TYPE_NEW_LAYER :
                            GIMP_PASTE_TYPE_NEW_LAYER_IN_PLACE;

          edit_paste (display, converted_type, merged, TRUE);
        }
      g_list_free (drawables);

      break;
    }
}

void
edit_paste_as_new_image_cmd_callback (GimpAction *action,
                                      GVariant   *value,
                                      gpointer    data)
{
  Gimp       *gimp;
  GtkWidget  *widget;
  GimpObject *paste;
  GimpImage  *image = NULL;
  return_if_no_gimp (gimp, data);
  return_if_no_widget (widget, data);

  paste = gimp_clipboard_get_object (gimp);

  if (paste)
    {
      image = gimp_edit_paste_as_new_image (gimp, paste,
                                            action_data_get_context (data));
      g_object_unref (paste);
    }

  if (image)
    {
      gimp_create_display (gimp, image, gimp_unit_pixel (), 1.0,
                           G_OBJECT (gimp_widget_get_monitor (widget)));
      g_object_unref (image);
    }
  else
    {
      gimp_message_literal (gimp, NULL, GIMP_MESSAGE_WARNING,
                            _("There is no image data in the clipboard "
                              "to paste."));
    }
}

void
edit_named_cut_cmd_callback (GimpAction *action,
                             GVariant   *value,
                             gpointer    data)
{
  GimpImage *image;
  GtkWidget *widget;
  GtkWidget *dialog;
  return_if_no_image (image, data);
  return_if_no_widget (widget, data);

  dialog = gimp_query_string_box (_("Cut Named"), widget,
                                  gimp_standard_help_func,
                                  GIMP_HELP_BUFFER_CUT,
                                  _("Enter a name for this buffer"),
                                  NULL,
                                  G_OBJECT (image), "disconnect",
                                  cut_named_buffer_callback,
                                  image, NULL);
  gtk_widget_show (dialog);
}

void
edit_named_copy_cmd_callback (GimpAction *action,
                              GVariant   *value,
                              gpointer    data)
{
  GimpImage *image;
  GtkWidget *widget;
  GtkWidget *dialog;
  return_if_no_image (image, data);
  return_if_no_widget (widget, data);

  dialog = gimp_query_string_box (_("Copy Named"), widget,
                                  gimp_standard_help_func,
                                  GIMP_HELP_BUFFER_COPY,
                                  _("Enter a name for this buffer"),
                                  NULL,
                                  G_OBJECT (image), "disconnect",
                                  copy_named_buffer_callback,
                                  image, NULL);
  gtk_widget_show (dialog);
}

void
edit_named_copy_visible_cmd_callback (GimpAction *action,
                                      GVariant   *value,
                                      gpointer    data)
{
  GimpImage *image;
  GtkWidget *widget;
  GtkWidget *dialog;
  return_if_no_image (image, data);
  return_if_no_widget (widget, data);

  dialog = gimp_query_string_box (_("Copy Visible Named"), widget,
                                  gimp_standard_help_func,
                                  GIMP_HELP_BUFFER_COPY,
                                  _("Enter a name for this buffer"),
                                  NULL,
                                  G_OBJECT (image), "disconnect",
                                  copy_named_visible_buffer_callback,
                                  image, NULL);
  gtk_widget_show (dialog);
}

void
edit_named_paste_cmd_callback (GimpAction *action,
                               GVariant   *value,
                               gpointer    data)
{
  Gimp      *gimp;
  GtkWidget *widget;
  return_if_no_gimp (gimp, data);
  return_if_no_widget (widget, data);

  gimp_window_strategy_show_dockable_dialog (GIMP_WINDOW_STRATEGY (gimp_get_window_strategy (gimp)),
                                             gimp,
                                             gimp_dialog_factory_get_singleton (),
                                             gimp_widget_get_monitor (widget),
                                             "gimp-buffer-list|gimp-buffer-grid");
}

void
edit_clear_cmd_callback (GimpAction *action,
                         GVariant   *value,
                         gpointer    data)
{
  GimpImage *image;
  GList     *drawables;
  GList     *iter;

  return_if_no_drawables (image, drawables, data);

  for (iter = drawables; iter; iter = iter->next)
    /* Return if any has a locked alpha. */
    if (! check_drawable_alpha (iter->data, data))
      {
        g_list_free (drawables);
        return;
      }

  gimp_image_undo_group_start (image, GIMP_UNDO_GROUP_PAINT,
                               _("Clear"));

  for (iter = drawables; iter; iter = iter->next)
    if (! gimp_viewable_get_children (GIMP_VIEWABLE (iter->data)) &&
        ! gimp_item_is_content_locked (GIMP_ITEM (iter->data), NULL))
      gimp_drawable_edit_clear (iter->data, action_data_get_context (data));

  gimp_image_undo_group_end (image);
  gimp_image_flush (image);
  g_list_free (drawables);
}

void
edit_fill_cmd_callback (GimpAction *action,
                        GVariant   *value,
                        gpointer    data)
{
  GimpImage       *image;
  GList           *drawables;
  GList           *iter;
  GimpFillType     fill_type;
  GimpFillOptions *options;
  GError          *error = NULL;
  return_if_no_drawables (image, drawables, data);

  fill_type = (GimpFillType) g_variant_get_int32 (value);

  options = gimp_fill_options_new (action_data_get_gimp (data), NULL, FALSE);

  if (gimp_fill_options_set_by_fill_type (options,
                                          action_data_get_context (data),
                                          fill_type, &error))
    {
      gimp_image_undo_group_start (image, GIMP_UNDO_GROUP_PAINT,
                                   gimp_fill_options_get_undo_desc (options));

      for (iter = drawables; iter; iter = iter->next)
        gimp_drawable_edit_fill (iter->data, options, NULL);

      gimp_image_undo_group_end (image);
      gimp_image_flush (image);
    }
  else
    {
      gimp_message_literal (image->gimp, NULL, GIMP_MESSAGE_WARNING,
                            error->message);
      g_clear_error (&error);
    }

  g_list_free (drawables);
  g_object_unref (options);
}


/*  private functions  */

static gboolean
check_drawable_alpha (GimpDrawable *drawable,
                      gpointer      data)
{
  GimpLayer *locked_layer = NULL;

  if (gimp_drawable_has_alpha (drawable) &&
      GIMP_IS_LAYER (drawable)           &&
      gimp_layer_is_alpha_locked (GIMP_LAYER (drawable), &locked_layer))
    {
      Gimp        *gimp    = action_data_get_gimp    (data);
      GimpDisplay *display = action_data_get_display (data);

      if (gimp && display)
        {
          gimp_message_literal (
            gimp, G_OBJECT (display), GIMP_MESSAGE_WARNING,
            _("A selected layer's alpha channel is locked."));

          gimp_tools_blink_lock_box (gimp, GIMP_ITEM (locked_layer));
        }

      return FALSE;
    }

  return TRUE;
}

static void
edit_paste (GimpDisplay   *display,
            GimpPasteType  paste_type,
            gboolean       merged,
            gboolean       try_svg)
{
  GimpImage  *image = gimp_display_get_image (display);
  GimpObject *paste;

  g_return_if_fail (paste_type != GIMP_PASTE_TYPE_NEW_LAYER_OR_FLOATING          &&
                    paste_type != GIMP_PASTE_TYPE_NEW_LAYER_OR_FLOATING_IN_PLACE &&
                    paste_type != GIMP_PASTE_TYPE_NEW_MERGED_LAYER_OR_FLOATING   &&
                    paste_type != GIMP_PASTE_TYPE_NEW_MERGED_LAYER_OR_FLOATING_IN_PLACE);

  if (try_svg)
    {
      gchar *svg;
      gsize  svg_size;

      svg = gimp_clipboard_get_svg (display->gimp, &svg_size);

      if (svg)
        {
          if (gimp_path_import_buffer (image, svg, svg_size,
                                       TRUE, FALSE,
                                       GIMP_IMAGE_ACTIVE_PARENT, -1,
                                       NULL, NULL))
            {
              gimp_image_flush (image);
            }

          g_free (svg);

          return;
        }
    }

  paste = gimp_clipboard_get_object (display->gimp);

  if (paste)
    {
      GimpDisplayShell *shell     = gimp_display_get_shell (display);
      GList            *drawables = gimp_image_get_selected_drawables (image);
      GList            *pasted_layers;
      gint              x, y;
      gint              width, height;

      if (g_list_length (drawables) != 1 ||
          (paste_type != GIMP_PASTE_TYPE_NEW_LAYER &&
           paste_type != GIMP_PASTE_TYPE_NEW_LAYER_IN_PLACE))
        {
          if (g_list_length (drawables) != 1)
            {
              gimp_message_literal (display->gimp, G_OBJECT (display),
                                    GIMP_MESSAGE_INFO,
                                    _("Pasted as new layer because the "
                                      "target is not a single layer or layer mask."));
            }
          else if (gimp_viewable_get_children (GIMP_VIEWABLE (drawables->data)))
            {
              gimp_message_literal (display->gimp, G_OBJECT (display),
                                    GIMP_MESSAGE_INFO,
                                    _("Pasted as new layer because the "
                                      "target is a layer group."));
            }
          else if (gimp_item_is_content_locked (GIMP_ITEM (drawables->data), NULL))
            {
              gimp_message_literal (display->gimp, G_OBJECT (display),
                                    GIMP_MESSAGE_INFO,
                                    _("Pasted as new layer because the "
                                      "target's pixels are locked."));
            }

          /* the actual paste-type conversion happens in gimp_edit_paste() */
        }

      gimp_display_shell_untransform_viewport (
        shell,
        ! gimp_display_shell_get_infinite_canvas (shell),
        &x, &y, &width, &height);

      if ((pasted_layers = gimp_edit_paste (image, drawables, paste, paste_type,
                                            gimp_get_user_context (display->gimp),
                                            merged, x, y, width, height)))
        {
          gimp_image_set_selected_layers (image, pasted_layers);

          g_list_free (pasted_layers);
          gimp_image_flush (image);
        }

      g_list_free (drawables);
      g_object_unref (paste);
    }
  else
    {
      gimp_message_literal (display->gimp, G_OBJECT (display),
                            GIMP_MESSAGE_WARNING,
                            _("There is no image data in the clipboard "
                              "to paste."));
    }
}

static void
cut_named_buffer_callback (GtkWidget   *widget,
                           const gchar *name,
                           gpointer     data)
{
  GimpImage *image     = GIMP_IMAGE (data);
  GList     *drawables = gimp_image_get_selected_drawables (image);
  GError    *error     = NULL;

  if (! drawables)
    {
      gimp_message_literal (image->gimp, NULL, GIMP_MESSAGE_WARNING,
                            _("There are no selected layers or channels to cut from."));
      return;
    }

  if (! (name && strlen (name)))
    name = _("(Unnamed Buffer)");

  if (gimp_edit_named_cut (image, name, drawables,
                           gimp_get_user_context (image->gimp), &error))
    {
      gimp_image_flush (image);
    }
  else
    {
      gimp_message_literal (image->gimp, NULL, GIMP_MESSAGE_WARNING,
                            error->message);
      g_clear_error (&error);
    }
  g_list_free (drawables);
}

static void
copy_named_buffer_callback (GtkWidget   *widget,
                            const gchar *name,
                            gpointer     data)
{
  GimpImage *image     = GIMP_IMAGE (data);
  GList     *drawables = gimp_image_get_selected_drawables (image);
  GError    *error     = NULL;

  if (! drawables)
    {
      gimp_message_literal (image->gimp, NULL, GIMP_MESSAGE_WARNING,
                            _("There are no selected layers or channels to copy from."));
      return;
    }

  if (! (name && strlen (name)))
    name = _("(Unnamed Buffer)");

  if (gimp_edit_named_copy (image, name, drawables,
                            gimp_get_user_context (image->gimp), &error))
    {
      gimp_image_flush (image);
    }
  else
    {
      gimp_message_literal (image->gimp, NULL, GIMP_MESSAGE_WARNING,
                            error->message);
      g_clear_error (&error);
    }
  g_list_free (drawables);
}

static void
copy_named_visible_buffer_callback (GtkWidget   *widget,
                                    const gchar *name,
                                    gpointer     data)
{
  GimpImage *image = GIMP_IMAGE (data);
  GError    *error = NULL;

  if (! (name && strlen (name)))
    name = _("(Unnamed Buffer)");

  if (gimp_edit_named_copy_visible (image, name,
                                    gimp_get_user_context (image->gimp),
                                    &error))
    {
      gimp_image_flush (image);
    }
  else
    {
      gimp_message_literal (image->gimp, NULL, GIMP_MESSAGE_WARNING,
                            error->message);
      g_clear_error (&error);
    }
}
