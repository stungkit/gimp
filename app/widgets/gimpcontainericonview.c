/* GIMP - The GNU Image Manipulation Program
 * Copyright (C) 1995 Spencer Kimball and Peter Mattis
 *
 * gimpcontainericonview.c
 * Copyright (C) 2010 Michael Natterer <mitch@gimp.org>
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

#include "libgimpwidgets/gimpwidgets.h"

#include "widgets-types.h"

#include "core/gimp.h"
#include "core/gimpcontainer.h"
#include "core/gimpcontext.h"
#include "core/gimpviewable.h"

#include "gimpcellrendererviewable.h"
#include "gimpcontainertreestore.h"
#include "gimpcontainericonview.h"
#include "gimpcontainerview.h"
#include "gimpdnd.h"
#include "gimpviewrenderer.h"
#include "gimpwidgets-utils.h"


struct _GimpContainerIconViewPrivate
{
  GimpViewRenderer *dnd_renderer;

  gulong            color_scheme_handler_id;
};


static void          gimp_container_icon_view_view_iface_init   (GimpContainerViewInterface  *iface);

static void          gimp_container_icon_view_constructed       (GObject                     *object);
static void          gimp_container_icon_view_finalize          (GObject                     *object);

static void          gimp_container_icon_view_unmap             (GtkWidget                   *widget);
static gboolean      gimp_container_icon_view_popup_menu        (GtkWidget                   *widget);

static void          gimp_container_icon_view_set_container     (GimpContainerView           *view,
                                                                 GimpContainer               *container);
static void          gimp_container_icon_view_set_context       (GimpContainerView           *view,
                                                                 GimpContext                 *context);
static void          gimp_container_icon_view_set_selection_mode(GimpContainerView           *view,
                                                                 GtkSelectionMode             mode);
static void          gimp_container_icon_view_set_view_size     (GimpContainerView           *view);
static gboolean      gimp_container_icon_view_set_selected      (GimpContainerView           *view,
                                                                 GList                       *items);
static gint          gimp_container_icon_view_get_selected        (GimpContainerView      *view,
                                                                   GList                 **items);

static gpointer      gimp_container_icon_view_insert_item       (GimpContainerView           *view,
                                                                 GimpViewable                *viewable,
                                                                 gpointer                     parent_insert_data,
                                                                 gint                         index);
static void          gimp_container_icon_view_remove_item       (GimpContainerView           *view,
                                                                 GimpViewable                *viewable,
                                                                 gpointer                     insert_data);
static void          gimp_container_icon_view_reorder_item      (GimpContainerView           *view,
                                                                 GimpViewable                *viewable,
                                                                 gint                         new_index,
                                                                 gpointer                     insert_data);
static void          gimp_container_icon_view_rename_item       (GimpContainerView           *view,
                                                                 GimpViewable                *viewable,
                                                                 gpointer                     insert_data);
static void          gimp_container_icon_view_clear_items       (GimpContainerView           *view);

static void          gimp_container_icon_view_invalidate        (GimpContainerIconView       *view);

static void          gimp_container_icon_view_selection_changed (GtkIconView                 *view,
                                                                 GimpContainerIconView       *icon_view);
static void          gimp_container_icon_view_item_activated    (GtkIconView                 *view,
                                                                 GtkTreePath                 *path,
                                                                 GimpContainerIconView       *icon_view);
static gboolean      gimp_container_icon_view_button_press      (GtkWidget                   *widget,
                                                                 GdkEventButton              *bevent,
                                                                 GimpContainerIconView       *icon_view);
static gboolean      gimp_container_icon_view_tooltip           (GtkWidget                   *widget,
                                                                 gint                         x,
                                                                 gint                         y,
                                                                 gboolean                     keyboard_tip,
                                                                 GtkTooltip                  *tooltip,
                                                                 GimpContainerIconView       *icon_view);

static GimpViewable * gimp_container_icon_view_drag_viewable    (GtkWidget    *widget,
                                                                 GimpContext **context,
                                                                 gpointer      data);
static GdkPixbuf   * gimp_container_icon_view_drag_pixbuf       (GtkWidget *widget,
                                                                 gpointer   data);
static gboolean      gimp_container_icon_view_get_1_selected    (GimpContainerIconView  *icon_view,
                                                                 GtkTreeIter            *iter);

static void          gimp_container_icon_view_trigger_redraw    (GimpContainerIconView *view);


G_DEFINE_TYPE_WITH_CODE (GimpContainerIconView, gimp_container_icon_view,
                         GIMP_TYPE_CONTAINER_BOX,
                         G_ADD_PRIVATE (GimpContainerIconView)
                         G_IMPLEMENT_INTERFACE (GIMP_TYPE_CONTAINER_VIEW,
                                                gimp_container_icon_view_view_iface_init))

#define parent_class gimp_container_icon_view_parent_class

static GimpContainerViewInterface *parent_view_iface = NULL;


static void
gimp_container_icon_view_class_init (GimpContainerIconViewClass *klass)
{
  GObjectClass   *object_class = G_OBJECT_CLASS (klass);
  GtkWidgetClass *widget_class = GTK_WIDGET_CLASS (klass);

  object_class->constructed = gimp_container_icon_view_constructed;
  object_class->finalize    = gimp_container_icon_view_finalize;

  widget_class->unmap       = gimp_container_icon_view_unmap;
  widget_class->popup_menu  = gimp_container_icon_view_popup_menu;

  gtk_widget_class_set_css_name (widget_class, "GimpContainerIconView");
}

static void
gimp_container_icon_view_view_iface_init (GimpContainerViewInterface *iface)
{
  parent_view_iface = g_type_interface_peek_parent (iface);

  if (! parent_view_iface)
    parent_view_iface = g_type_default_interface_peek (GIMP_TYPE_CONTAINER_VIEW);

  iface->set_container      = gimp_container_icon_view_set_container;
  iface->set_context        = gimp_container_icon_view_set_context;
  iface->set_selection_mode = gimp_container_icon_view_set_selection_mode;
  iface->set_view_size      = gimp_container_icon_view_set_view_size;
  iface->set_selected       = gimp_container_icon_view_set_selected;
  iface->get_selected       = gimp_container_icon_view_get_selected;

  iface->insert_item        = gimp_container_icon_view_insert_item;
  iface->remove_item        = gimp_container_icon_view_remove_item;
  iface->reorder_item       = gimp_container_icon_view_reorder_item;
  iface->rename_item        = gimp_container_icon_view_rename_item;
  iface->clear_items        = gimp_container_icon_view_clear_items;

  iface->insert_data_free = (GDestroyNotify) gtk_tree_iter_free;
}

static void
gimp_container_icon_view_init (GimpContainerIconView *icon_view)
{
  GimpContainerBox *box = GIMP_CONTAINER_BOX (icon_view);

  icon_view->priv = gimp_container_icon_view_get_instance_private (icon_view);

  icon_view->priv->color_scheme_handler_id = 0;

  gimp_container_tree_store_columns_init (icon_view->model_columns,
                                          &icon_view->n_model_columns);

  gtk_scrolled_window_set_shadow_type (GTK_SCROLLED_WINDOW (box->scrolled_win),
                                       GTK_SHADOW_IN);
  gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (box->scrolled_win),
                                  GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
}

static void
gimp_container_icon_view_constructed (GObject *object)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (object);
  GimpContainerView     *view      = GIMP_CONTAINER_VIEW (object);
  GimpContainerBox      *box       = GIMP_CONTAINER_BOX (object);

  G_OBJECT_CLASS (parent_class)->constructed (object);

  icon_view->model = gimp_container_tree_store_new (view,
                                                    icon_view->n_model_columns,
                                                    icon_view->model_columns);

  icon_view->view = g_object_new (GTK_TYPE_ICON_VIEW,
                                  "model",          icon_view->model,
                                  "row-spacing",    0,
                                  "column-spacing", 0,
                                  "margin",         0,
                                  "item-padding",   1,
                                  "has-tooltip",    TRUE,
                                  NULL);
  g_object_unref (icon_view->model);

  gtk_container_add (GTK_CONTAINER (box->scrolled_win),
                     GTK_WIDGET (icon_view->view));
  gtk_widget_show (GTK_WIDGET (icon_view->view));

  gimp_container_view_set_dnd_widget (view, GTK_WIDGET (icon_view->view));

  icon_view->renderer_cell = gimp_cell_renderer_viewable_new ();
  gtk_cell_layout_pack_start (GTK_CELL_LAYOUT (icon_view->view),
                              icon_view->renderer_cell,
                              FALSE);
  gtk_cell_layout_set_attributes (GTK_CELL_LAYOUT (icon_view->view),
                                  icon_view->renderer_cell,
                                  "renderer", GIMP_CONTAINER_TREE_STORE_COLUMN_RENDERER,
                                  NULL);

  gimp_container_tree_store_add_renderer_cell (GIMP_CONTAINER_TREE_STORE (icon_view->model),
                                               icon_view->renderer_cell, -1);

  g_signal_connect (icon_view->view, "selection-changed",
                    G_CALLBACK (gimp_container_icon_view_selection_changed),
                    icon_view);
  g_signal_connect (icon_view->view, "item-activated",
                    G_CALLBACK (gimp_container_icon_view_item_activated),
                    icon_view);
  g_signal_connect (icon_view->view, "query-tooltip",
                    G_CALLBACK (gimp_container_icon_view_tooltip),
                    icon_view);
}

static void
gimp_container_icon_view_finalize (GObject *object)
{
  //GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (object);

  G_OBJECT_CLASS (parent_class)->finalize (object);
}

static void
gimp_container_icon_view_unmap (GtkWidget *widget)
{
  //GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (widget);

  GTK_WIDGET_CLASS (parent_class)->unmap (widget);
}

static gboolean
gimp_container_icon_view_popup_menu (GtkWidget *widget)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (widget);
  GtkTreeIter            iter;
  GtkTreePath           *path;
  GdkRectangle           rect;

  if (! gimp_container_icon_view_get_1_selected (icon_view, &iter))
    return FALSE;

  path = gtk_tree_model_get_path (icon_view->model, &iter);
  gtk_icon_view_get_cell_rect (icon_view->view, path, NULL, &rect);
  gtk_tree_path_free (path);

  return gimp_editor_popup_menu_at_rect (GIMP_EDITOR (widget),
                                         gtk_widget_get_window (GTK_WIDGET (icon_view->view)),
                                         &rect, GDK_GRAVITY_CENTER, GDK_GRAVITY_NORTH_WEST,
                                         NULL);
}

GtkWidget *
gimp_container_icon_view_new (GimpContainer *container,
                              GimpContext   *context,
                              gint           view_size,
                              gint           view_border_width)
{
  GimpContainerIconView *icon_view;
  GimpContainerView     *view;

  g_return_val_if_fail (container == NULL || GIMP_IS_CONTAINER (container),
                        NULL);
  g_return_val_if_fail (context == NULL || GIMP_IS_CONTEXT (context), NULL);
  g_return_val_if_fail (view_size > 0 &&
                        view_size <= GIMP_VIEWABLE_MAX_PREVIEW_SIZE, NULL);
  g_return_val_if_fail (view_border_width >= 0 &&
                        view_border_width <= GIMP_VIEW_MAX_BORDER_WIDTH,
                        NULL);

  icon_view = g_object_new (GIMP_TYPE_CONTAINER_ICON_VIEW, NULL);

  view = GIMP_CONTAINER_VIEW (icon_view);

  gimp_container_view_set_view_size (view, view_size, 0 /* ignore border */);

  if (container)
    gimp_container_view_set_container (view, container);

  if (context)
    gimp_container_view_set_context (view, context);

  return GTK_WIDGET (icon_view);
}


/*  GimpContainerView methods  */

static void
gimp_container_icon_view_set_container (GimpContainerView *view,
                                        GimpContainer     *container)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (view);
  GimpContainer         *old_container;

  old_container = gimp_container_view_get_container (view);

  if (old_container)
    {
      if (! container)
        {
          if (gimp_dnd_viewable_source_remove (GTK_WIDGET (icon_view->view),
                                               gimp_container_get_child_type (old_container)))
            {
              if (GIMP_VIEWABLE_CLASS (g_type_class_peek (gimp_container_get_child_type (old_container)))->get_size)
                gimp_dnd_pixbuf_source_remove (GTK_WIDGET (icon_view->view));

              gtk_drag_source_unset (GTK_WIDGET (icon_view->view));
            }

          g_signal_handlers_disconnect_by_func (icon_view->view,
                                                gimp_container_icon_view_button_press,
                                                icon_view);
        }
    }
  else if (container)
    {
      if (gimp_dnd_drag_source_set_by_type (GTK_WIDGET (icon_view->view),
                                            GDK_BUTTON1_MASK | GDK_BUTTON2_MASK,
                                            gimp_container_get_child_type (container),
                                            GDK_ACTION_COPY))
        {
          gimp_dnd_viewable_source_add (GTK_WIDGET (icon_view->view),
                                        gimp_container_get_child_type (container),
                                        gimp_container_icon_view_drag_viewable,
                                        icon_view);

          if (GIMP_VIEWABLE_CLASS (g_type_class_peek (gimp_container_get_child_type (container)))->get_size)
            gimp_dnd_pixbuf_source_add (GTK_WIDGET (icon_view->view),
                                        gimp_container_icon_view_drag_pixbuf,
                                        icon_view);
        }

      g_signal_connect_object (icon_view->view, "button-press-event",
                               G_CALLBACK (gimp_container_icon_view_button_press),
                               icon_view,
                               0);
    }

  parent_view_iface->set_container (view, container);
}

static void
gimp_container_icon_view_set_context (GimpContainerView *view,
                                      GimpContext       *context)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (view);

  parent_view_iface->set_context (view, context);

  if (icon_view->model)
    gimp_container_tree_store_set_context (GIMP_CONTAINER_TREE_STORE (icon_view->model),
                                           context);

  if (context != NULL)
    {
      if (icon_view->priv->color_scheme_handler_id == 0)
        icon_view->priv->color_scheme_handler_id =
          g_signal_connect_object (context->gimp->config,
                                   "notify::theme-color-scheme",
                                   G_CALLBACK (gimp_container_icon_view_trigger_redraw),
                                   icon_view, G_CONNECT_SWAPPED);
    }
}

static void
gimp_container_icon_view_set_selection_mode (GimpContainerView *view,
                                             GtkSelectionMode   mode)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (view);

  gtk_icon_view_set_selection_mode (icon_view->view, mode);

  parent_view_iface->set_selection_mode (view, mode);
}

static void
gimp_container_icon_view_set_view_size (GimpContainerView *view)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (view);

  if (icon_view->model)
    gimp_container_tree_store_set_view_size (GIMP_CONTAINER_TREE_STORE (icon_view->model));

  if (icon_view->view)
    {
      gtk_icon_view_set_columns (icon_view->view, -1);
      gtk_icon_view_set_item_width (icon_view->view, -1);

      /* ugly workaround to force the icon view to invalidate all its
       * cached icon sizes
       */
      gtk_icon_view_set_item_orientation (icon_view->view,
                                          GTK_ORIENTATION_VERTICAL);
      gtk_icon_view_set_item_orientation (icon_view->view,
                                          GTK_ORIENTATION_HORIZONTAL);
    }
}

static gboolean
gimp_container_icon_view_set_selected (GimpContainerView *view,
                                       GList             *items)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (view);
  GList                 *list;

  if (items)
    {
      GList *paths = NULL;

      for (list = items; list; list = g_list_next (list))
        {
          GtkTreeIter *iter;
          GtkTreePath *path;

          iter = _gimp_container_view_lookup (GIMP_CONTAINER_VIEW (view),
                                              list->data);
          if (! iter)
            /* This may happen when the GimpContainerIconView only
             * shows a subpart of the whole icons. We don't select
             * what is not shown.
             */
            continue;

          path = gtk_tree_model_get_path (icon_view->model, iter);

          paths = g_list_prepend (paths, path);
        }

      paths = g_list_reverse (paths);

      g_signal_handlers_block_by_func (icon_view->view,
                                       gimp_container_icon_view_selection_changed,
                                       icon_view);

      gtk_icon_view_unselect_all (icon_view->view);

      for (list = paths; list; list = g_list_next (list))
        {
          gtk_icon_view_select_path (icon_view->view, list->data);
        }

      if (paths)
        {
          gtk_icon_view_set_cursor (icon_view->view, paths->data, NULL, FALSE);
          gtk_icon_view_scroll_to_path (icon_view->view, paths->data,
                                        FALSE, 0.0, 0.0);
        }

      g_signal_handlers_unblock_by_func (icon_view->view,
                                         gimp_container_icon_view_selection_changed,
                                         icon_view);

      g_list_free_full (paths, (GDestroyNotify) gtk_tree_path_free);

      _gimp_container_view_selection_changed (GIMP_CONTAINER_VIEW (icon_view));
    }
  else
    {
      gtk_icon_view_unselect_all (icon_view->view);
    }

  return TRUE;
}

static gint
gimp_container_icon_view_get_selected (GimpContainerView  *view,
                                       GList             **items)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (view);
  GList                 *paths;
  gint                   n_selected;

  paths = gtk_icon_view_get_selected_items (icon_view->view);
  n_selected = g_list_length (paths);

  if (items)
    {
      GList *list;

      for (list = paths; list; list = g_list_next (list))
        {
          GtkTreeIter       iter;
          GimpViewRenderer *renderer;

          gtk_tree_model_get_iter (GTK_TREE_MODEL (icon_view->model), &iter,
                                   (GtkTreePath *) list->data);

          gtk_tree_model_get (icon_view->model, &iter,
                              GIMP_CONTAINER_TREE_STORE_COLUMN_RENDERER, &renderer,
                              -1);

          if (renderer->viewable)
            *items = g_list_prepend (*items, renderer->viewable);

          g_object_unref (renderer);
        }

      *items = g_list_reverse (*items);
    }

  g_list_free_full (paths, (GDestroyNotify) gtk_tree_path_free);

  return n_selected;
}

static gpointer
gimp_container_icon_view_insert_item (GimpContainerView *view,
                                      GimpViewable      *viewable,
                                      gpointer           parent_insert_data,
                                      gint               index)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (view);
  GtkTreeIter           *iter;

  iter = gimp_container_tree_store_insert_item (GIMP_CONTAINER_TREE_STORE (icon_view->model),
                                                viewable,
                                                parent_insert_data,
                                                index);

  if (parent_insert_data)
    {
#if 0
      GtkTreePath *path = gtk_tree_model_get_path (icon_view->model, iter);

      gtk_icon_view_expand_to_path (icon_view->view, path);

      gtk_tree_path_free (path);
#endif
    }

  return iter;
}

static void
gimp_container_icon_view_remove_item (GimpContainerView *view,
                                      GimpViewable      *viewable,
                                      gpointer           insert_data)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (view);

  gimp_container_tree_store_remove_item (GIMP_CONTAINER_TREE_STORE (icon_view->model),
                                         viewable,
                                         insert_data);
}

static void
gimp_container_icon_view_reorder_item (GimpContainerView *view,
                                       GimpViewable      *viewable,
                                       gint               new_index,
                                       gpointer           insert_data)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (view);
  GtkTreeIter           *iter      = (GtkTreeIter *) insert_data;
  gboolean               selected  = FALSE;

  if (iter)
    {
      GtkTreeIter selected_iter;

      selected = gimp_container_icon_view_get_1_selected (icon_view,
                                                          &selected_iter);

      if (selected)
        {
          GimpViewRenderer *renderer;

          gtk_tree_model_get (icon_view->model, &selected_iter,
                              GIMP_CONTAINER_TREE_STORE_COLUMN_RENDERER, &renderer,
                              -1);

          if (renderer->viewable != viewable)
            selected = FALSE;

          g_object_unref (renderer);
        }
    }

  gimp_container_tree_store_reorder_item (GIMP_CONTAINER_TREE_STORE (icon_view->model),
                                          viewable,
                                          new_index,
                                          iter);

  if (selected)
    gimp_container_view_set_1_selected (view, viewable);
}

static void
gimp_container_icon_view_rename_item (GimpContainerView *view,
                                      GimpViewable      *viewable,
                                      gpointer           insert_data)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (view);
  GtkTreeIter           *iter      = (GtkTreeIter *) insert_data;

  gimp_container_tree_store_rename_item (GIMP_CONTAINER_TREE_STORE (icon_view->model),
                                         viewable,
                                         iter);
}

static void
gimp_container_icon_view_clear_items (GimpContainerView *view)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (view);

  gimp_container_tree_store_clear_items (GIMP_CONTAINER_TREE_STORE (icon_view->model));

  parent_view_iface->clear_items (view);
}

static void
gimp_container_icon_view_invalidate (GimpContainerIconView *view)
{
  GtkTreeIter iter;
  gboolean    iter_valid;

  for (iter_valid = gtk_tree_model_get_iter_first (view->model, &iter);
       iter_valid;
       iter_valid = gtk_tree_model_iter_next (view->model, &iter))
    {
      GimpViewRenderer *renderer;

      gtk_tree_model_get (view->model, &iter,
                          GIMP_CONTAINER_TREE_STORE_COLUMN_RENDERER, &renderer,
                          -1);

      gimp_view_renderer_invalidate (renderer);

      g_object_unref (renderer);
    }
}


/*  callbacks  */

static void
gimp_container_icon_view_selection_changed (GtkIconView           *gtk_icon_view,
                                            GimpContainerIconView *icon_view)
{
  _gimp_container_view_selection_changed (GIMP_CONTAINER_VIEW (icon_view));
}

static void
gimp_container_icon_view_item_activated (GtkIconView           *view,
                                         GtkTreePath           *path,
                                         GimpContainerIconView *icon_view)
{
  GtkTreeIter       iter;
  GimpViewRenderer *renderer;

  gtk_tree_model_get_iter (icon_view->model, &iter, path);

  gtk_tree_model_get (icon_view->model, &iter,
                      GIMP_CONTAINER_TREE_STORE_COLUMN_RENDERER, &renderer,
                      -1);

  _gimp_container_view_item_activated (GIMP_CONTAINER_VIEW (icon_view),
                                       renderer->viewable);

  g_object_unref (renderer);
}

static gboolean
gimp_container_icon_view_button_press (GtkWidget             *widget,
                                       GdkEventButton        *bevent,
                                       GimpContainerIconView *icon_view)
{
  GimpContainerView *container_view = GIMP_CONTAINER_VIEW (icon_view);
  GtkTreePath       *path;
  gboolean           handled        = GDK_EVENT_PROPAGATE;

  icon_view->priv->dnd_renderer = NULL;

  path = gtk_icon_view_get_path_at_pos (GTK_ICON_VIEW (widget),
                                        bevent->x, bevent->y);

  if (path)
    {
      GimpViewRenderer *renderer;
      GtkTreeIter       iter;

      gtk_tree_model_get_iter (icon_view->model, &iter, path);

      gtk_tree_model_get (icon_view->model, &iter,
                          GIMP_CONTAINER_TREE_STORE_COLUMN_RENDERER, &renderer,
                          -1);

      icon_view->priv->dnd_renderer = renderer;

      if (gdk_event_triggers_context_menu ((GdkEvent *) bevent))
        {
          /* If the clicked item is not selected, it becomes the new
           * selection. Otherwise, we use the current selection. This
           * allows to not break multiple selection when right-clicking.
           */
          if (! gimp_container_view_is_item_selected (container_view,
                                                      renderer->viewable))
            {
              gimp_container_view_set_1_selected (container_view,
                                                  renderer->viewable);
            }

          /* Show the context menu. */
          handled = gimp_editor_popup_menu_at_pointer (GIMP_EDITOR (icon_view),
                                                       (GdkEvent *) bevent);
        }
      /* Else LMB down or similar.  Propagate. */

      g_object_unref (renderer);

      gtk_tree_path_free (path);
    }
  else
    {
      /* Button down outside any item. */
      if (gdk_event_triggers_context_menu ((GdkEvent *) bevent))
        {
          (void) gimp_editor_popup_menu_at_pointer (GIMP_EDITOR (icon_view),
                                                    (GdkEvent *) bevent);
          /* Discard result.  Does not actually popup menu except for managed views, e.g. dockable. */
        }

      /* Else LMB down or similar.
       * Say we handled event even though we did nothing.
       * Otherwise, selection is set to None.
       */
      handled = GDK_EVENT_STOP;
    }

  return handled;
}

static gboolean
gimp_container_icon_view_tooltip (GtkWidget             *widget,
                                  gint                   x,
                                  gint                   y,
                                  gboolean               keyboard_tip,
                                  GtkTooltip            *tooltip,
                                  GimpContainerIconView *icon_view)
{
  GimpViewRenderer *renderer;
  GtkTreeIter       iter;
  GtkTreePath      *path;
  gboolean          show_tip = FALSE;

  if (! gtk_icon_view_get_tooltip_context (GTK_ICON_VIEW (widget), &x, &y,
                                           keyboard_tip,
                                           NULL, &path, &iter))
    return FALSE;

  gtk_tree_model_get (icon_view->model, &iter,
                      GIMP_CONTAINER_TREE_STORE_COLUMN_RENDERER, &renderer,
                      -1);

  if (renderer)
    {
      gchar *desc;
      gchar *tip;

      desc = gimp_viewable_get_description (renderer->viewable, &tip);

      if (tip)
        {
          gtk_tooltip_set_text (tooltip, tip);
          gtk_icon_view_set_tooltip_cell (GTK_ICON_VIEW (widget), tooltip, path,
                                          icon_view->renderer_cell);

          show_tip = TRUE;

          g_free (tip);
        }

      g_free (desc);
      g_object_unref (renderer);
    }

  gtk_tree_path_free (path);

  return show_tip;
}

static GimpViewable *
gimp_container_icon_view_drag_viewable (GtkWidget    *widget,
                                        GimpContext **context,
                                        gpointer      data)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (data);

  if (context)
    *context = gimp_container_view_get_context (GIMP_CONTAINER_VIEW (data));

  if (icon_view->priv->dnd_renderer)
    return icon_view->priv->dnd_renderer->viewable;

  return NULL;
}

static GdkPixbuf *
gimp_container_icon_view_drag_pixbuf (GtkWidget *widget,
                                      gpointer   data)
{
  GimpContainerIconView *icon_view = GIMP_CONTAINER_ICON_VIEW (data);
  GimpViewRenderer      *renderer  = icon_view->priv->dnd_renderer;
  gint                   width;
  gint                   height;

  if (renderer && gimp_viewable_get_size (renderer->viewable, &width, &height))
    {
      GeglColor *color      = NULL;
      GeglColor *background = NULL;
      GdkPixbuf *pixbuf;

      if (renderer->context)
        {
          gboolean follow_theme = FALSE;

          g_object_get (renderer->context->gimp->config,
                        "viewables-follow-theme", &follow_theme,
                        NULL);
          if (follow_theme)
            {
              GtkStyleContext *style;
              GdkRGBA         *fg_color = NULL;
              GdkRGBA         *bg_color = NULL;

              style = gtk_widget_get_style_context (widget);
              gtk_style_context_get (style, gtk_style_context_get_state (style),
                                     GTK_STYLE_PROPERTY_COLOR,            &fg_color,
                                     GTK_STYLE_PROPERTY_BACKGROUND_COLOR, &bg_color,
                                     NULL);
              if (fg_color && bg_color)
                {
                  color = gegl_color_new (NULL);
                  gegl_color_set_rgba_with_space (color,
                                                  fg_color->red, fg_color->green, fg_color->blue, 1.0,
                                                  NULL);

                  background = gegl_color_new (NULL);
                  gegl_color_set_rgba_with_space (background,
                                                  bg_color->red, bg_color->green, bg_color->blue, 1.0,
                                                  NULL);
                }
              g_clear_pointer (&fg_color, gdk_rgba_free);
              g_clear_pointer (&bg_color, gdk_rgba_free);
            }
        }

      pixbuf = gimp_viewable_get_new_pixbuf (renderer->viewable,
                                             renderer->context,
                                             width, height, color, background);

      g_clear_object (&color);
      g_clear_object (&background);

      return pixbuf;
    }

  return NULL;
}

static gboolean
gimp_container_icon_view_get_1_selected (GimpContainerIconView *icon_view,
                                         GtkTreeIter           *iter)
{
  GList    *selected_items;
  gboolean  retval;

  selected_items = gtk_icon_view_get_selected_items (icon_view->view);

  if (g_list_length (selected_items) == 1)
    {
      gtk_tree_model_get_iter (GTK_TREE_MODEL (icon_view->model), iter,
                               selected_items->data);

      retval = TRUE;
    }
  else
    {
      retval = FALSE;
    }

  g_list_free_full (selected_items, (GDestroyNotify) gtk_tree_path_free);

  return retval;
}

static void
gimp_container_icon_view_trigger_redraw (GimpContainerIconView *view)
{
  gimp_container_icon_view_invalidate (view);
  gtk_widget_queue_draw (GTK_WIDGET (view));
}
