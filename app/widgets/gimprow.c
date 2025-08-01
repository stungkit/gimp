/* GIMP - The GNU Image Manipulation Program
 * Copyright (C) 1995 Spencer Kimball and Peter Mattis
 *
 * gimprow.c
 * Copyright (C) 2001-2006 Michael Natterer <mitch@gimp.org>
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

#include <gegl.h>
#include <gtk/gtk.h>

#include "widgets-types.h"

#include "core/gimpviewable.h"

#include "gimprow.h"

#include "gimp-intl.h"


enum
{
  EDIT_NAME,
  LAST_SIGNAL
};

enum
{
  PROP_0,
  PROP_VIEWABLE,
  N_PROPS
};

static GParamSpec *obj_props[N_PROPS] = { NULL, };


typedef struct _GimpRowPrivate GimpRowPrivate;

struct _GimpRowPrivate
{
  GimpViewable *viewable;

  GtkWidget    *icon;

  GtkWidget    *label;
  GtkWidget    *popover;
  GtkGesture   *label_long_press;
};

#define GET_PRIVATE(obj) \
  ((GimpRowPrivate *) gimp_row_get_instance_private ((GimpRow *) obj))


static void       gimp_row_constructed           (GObject          *object);
static void       gimp_row_dispose               (GObject          *object);
static void       gimp_row_set_property          (GObject          *object,
                                                  guint             property_id,
                                                  const GValue     *value,
                                                  GParamSpec       *pspec);
static void       gimp_row_get_property          (GObject          *object,
                                                  guint             property_id,
                                                  GValue           *value,
                                                  GParamSpec       *pspec);

static void       gimp_row_real_edit_name        (GimpRow          *row);
static gboolean   gimp_row_real_name_edited      (GimpRow          *row,
                                                  const gchar      *new_name);

static void       gimp_row_viewable_icon_changed (GimpViewable     *viewable,
                                                  const GParamSpec *pspec,
                                                  GimpRow          *row);
static void       gimp_row_viewable_name_changed (GimpViewable     *viewable,
                                                  GimpRow          *row);

static void       gimp_row_label_long_pressed    (GtkGesture       *gesture,
                                                  gdouble           x,
                                                  gdouble           y,
                                                  GimpRow          *row);
static void       gimp_row_rename_entry_activate (GtkEntry         *entry,
                                                  GimpRow          *row);


G_DEFINE_TYPE_WITH_PRIVATE (GimpRow, gimp_row, GTK_TYPE_LIST_BOX_ROW)

#define parent_class gimp_row_parent_class

static guint row_signals[LAST_SIGNAL] = { 0 };


static void
gimp_row_class_init (GimpRowClass *klass)
{
  GObjectClass  *object_class = G_OBJECT_CLASS (klass);
  GtkBindingSet *binding_set;

  object_class->constructed  = gimp_row_constructed;
  object_class->dispose      = gimp_row_dispose;
  object_class->set_property = gimp_row_set_property;
  object_class->get_property = gimp_row_get_property;

  klass->edit_name           = gimp_row_real_edit_name;
  klass->name_edited         = gimp_row_real_name_edited;

  row_signals[EDIT_NAME] =
    g_signal_new ("edit-name",
                  G_TYPE_FROM_CLASS (klass),
                  G_SIGNAL_RUN_LAST | G_SIGNAL_ACTION,
                  G_STRUCT_OFFSET (GimpRowClass, edit_name),
                  NULL, NULL, NULL,
                  G_TYPE_NONE, 0);

  binding_set = gtk_binding_set_by_class (klass);

  gtk_binding_entry_add_signal (binding_set, GDK_KEY_F2, 0,
                                "edit-name", 0);

  obj_props[PROP_VIEWABLE] =
    g_param_spec_object ("viewable",
                         NULL, NULL,
                         GIMP_TYPE_VIEWABLE,
                         G_PARAM_READWRITE |
                         G_PARAM_CONSTRUCT);

  g_object_class_install_properties (object_class, N_PROPS, obj_props);
}

static void
gimp_row_init (GimpRow *row)
{
}

static void
gimp_row_constructed (GObject *object)
{
  GimpRowPrivate *priv = GET_PRIVATE (object);
  GtkWidget      *box;
  GtkWidget      *ebox;
  const gchar    *icon_name   = NULL;
  const gchar    *description = NULL;

  G_OBJECT_CLASS (parent_class)->constructed (object);

  if (priv->viewable)
    {
      icon_name   = gimp_viewable_get_icon_name (priv->viewable);
      description = gimp_viewable_get_description (priv->viewable, NULL);
    }

  box = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 6);
  gtk_container_add (GTK_CONTAINER (object), box);
  gtk_widget_set_visible (box, TRUE);

  priv->icon = gtk_image_new_from_icon_name (icon_name, GTK_ICON_SIZE_BUTTON);
  gtk_box_pack_start (GTK_BOX (box), priv->icon, FALSE, FALSE, 0);
  gtk_widget_set_visible (priv->icon, TRUE);

  ebox = gtk_event_box_new ();
  gtk_box_pack_start (GTK_BOX (box), ebox, TRUE, TRUE, 0);
  gtk_widget_show (ebox);

  priv->label = gtk_label_new (description);
  gtk_label_set_xalign (GTK_LABEL (priv->label), 0.0);
  gtk_container_add (GTK_CONTAINER (ebox), priv->label);
  gtk_widget_set_visible (priv->label, TRUE);

  priv->label_long_press = gtk_gesture_long_press_new (ebox);
  gtk_event_controller_set_propagation_phase (GTK_EVENT_CONTROLLER (priv->label_long_press),
                                              GTK_PHASE_TARGET);

  g_signal_connect (priv->label_long_press, "pressed",
                    G_CALLBACK (gimp_row_label_long_pressed),
                    object);
}

static void
gimp_row_dispose (GObject *object)
{
  GimpRowPrivate *priv = GET_PRIVATE (object);

  gimp_row_set_viewable (GIMP_ROW (object), NULL);

  g_clear_pointer (&priv->popover, gtk_widget_destroy);

  G_OBJECT_CLASS (parent_class)->dispose (object);
}

static void
gimp_row_set_property (GObject      *object,
                       guint         property_id,
                       const GValue *value,
                       GParamSpec   *pspec)
{
  switch (property_id)
    {
    case PROP_VIEWABLE:
      gimp_row_set_viewable (GIMP_ROW (object), g_value_get_object (value));
      break;

    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
      break;
    }
}

static void
gimp_row_get_property (GObject    *object,
                       guint       property_id,
                       GValue     *value,
                       GParamSpec *pspec)
{
  switch (property_id)
    {
    case PROP_VIEWABLE:
      g_value_set_object (value, gimp_row_get_viewable (GIMP_ROW (object)));
      break;

    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
      break;
    }
}

static void
gimp_row_real_edit_name (GimpRow *row)
{
  GimpRowPrivate *priv = GET_PRIVATE (row);
  GtkWidget      *box;
  GtkWidget      *label;
  GtkWidget      *entry;
  const gchar    *default_name;
  gchar          *title;

  if (! priv->viewable ||
      ! gimp_viewable_is_name_editable (priv->viewable))
    {
      gtk_widget_error_bell (GTK_WIDGET (row));
      return;
    }

  priv->popover = gtk_popover_new (priv->label);
  gtk_popover_set_modal (GTK_POPOVER (priv->popover), TRUE);

  g_object_add_weak_pointer (G_OBJECT (priv->popover),
                             (gpointer) &priv->popover);

  g_signal_connect (priv->popover, "closed",
                    G_CALLBACK (gtk_widget_destroy),
                    NULL);

  box = gtk_box_new (GTK_ORIENTATION_VERTICAL, 6);
  gtk_container_add (GTK_CONTAINER (priv->popover), box);
  gtk_widget_show (box);

  default_name = GIMP_VIEWABLE_GET_CLASS (priv->viewable)->default_name;

  title = g_strdup_printf (_("Rename %s"), default_name);
  label = gtk_label_new (title);
  g_free (title);
  gtk_box_pack_start (GTK_BOX (box), label, FALSE, FALSE, 0);
  gtk_widget_show (label);

  entry = gtk_entry_new ();
  gtk_entry_set_text (GTK_ENTRY (entry), gimp_object_get_name (priv->viewable));
  gtk_box_pack_start (GTK_BOX (box), entry, FALSE, FALSE, 0);
  gtk_widget_show (entry);

  g_signal_connect (entry, "activate",
                    G_CALLBACK (gimp_row_rename_entry_activate),
                    row);

  gtk_popover_popup (GTK_POPOVER (priv->popover));
}

static gboolean
gimp_row_real_name_edited (GimpRow     *row,
                           const gchar *new_name)
{
  GimpRowPrivate *priv = GET_PRIVATE (row);

  gimp_object_set_name (GIMP_OBJECT (priv->viewable), new_name);

  return TRUE;
}


/*  public functions  */

GtkWidget *
gimp_row_new (GimpViewable *viewable)
{
  g_return_val_if_fail (GIMP_IS_VIEWABLE (viewable), NULL);

  return g_object_new (GIMP_TYPE_ROW,
                       "viewable", viewable,
                       NULL);
}

void
gimp_row_set_viewable (GimpRow      *row,
                       GimpViewable *viewable)
{
  GimpRowPrivate *priv = GET_PRIVATE (row);

  g_return_if_fail (GIMP_IS_ROW (row));
  g_return_if_fail (viewable == NULL || GIMP_IS_VIEWABLE (viewable));

  if (priv->viewable)
    {
      g_signal_handlers_disconnect_by_func (priv->viewable,
                                            gimp_row_viewable_icon_changed,
                                            row);
      g_signal_handlers_disconnect_by_func (priv->viewable,
                                            gimp_row_viewable_name_changed,
                                            row);
    }

  g_set_object (&priv->viewable, viewable);

  if (priv->viewable)
    {
      g_signal_connect (priv->viewable,
                        GIMP_VIEWABLE_GET_CLASS (viewable)->name_changed_signal,
                        G_CALLBACK (gimp_row_viewable_name_changed),
                        row);
      gimp_row_viewable_name_changed (priv->viewable, row);

      g_signal_connect (priv->viewable, "notify::icon-name",
                        G_CALLBACK (gimp_row_viewable_icon_changed),
                        row);
      gimp_row_viewable_icon_changed (priv->viewable, NULL, row);
    }

  g_object_notify_by_pspec (G_OBJECT (row), obj_props[PROP_VIEWABLE]);
}

GimpViewable *
gimp_row_get_viewable (GimpRow *row)
{
  g_return_val_if_fail (GIMP_IS_ROW (row), NULL);

  return GET_PRIVATE (row)->viewable;
}

GtkWidget *
gimp_row_create (gpointer item,
                 gpointer user_data)
{
  return gimp_row_new (item);
}


/*  private functions  */

static void
gimp_row_viewable_icon_changed (GimpViewable     *viewable,
                                const GParamSpec *pspec,
                                GimpRow          *row)
{
  GimpRowPrivate *priv = GET_PRIVATE (row);

  if (priv->icon)
    gtk_image_set_from_icon_name (GTK_IMAGE (priv->icon),
                                  gimp_viewable_get_icon_name (viewable),
                                  GTK_ICON_SIZE_BUTTON);
}

static void
gimp_row_viewable_name_changed (GimpViewable *viewable,
                                GimpRow      *row)
{
  GimpRowPrivate *priv = GET_PRIVATE (row);

  if (priv->label)
    gtk_label_set_text (GTK_LABEL (priv->label),
                        gimp_viewable_get_description (viewable, NULL));
}

static void
gimp_row_label_long_pressed (GtkGesture *gesture,
                             gdouble     x,
                             gdouble     y,
                             GimpRow    *row)
{
  g_signal_emit (row, row_signals[EDIT_NAME], 0);
}

static void
gimp_row_rename_entry_activate (GtkEntry *entry,
                                GimpRow  *row)
{
  GimpRowPrivate *priv = GET_PRIVATE (row);
  const gchar    *old_name;
  const gchar    *new_name;

  old_name = gimp_object_get_name (priv->viewable);
  new_name = gtk_entry_get_text (entry);

  if (g_strcmp0 (old_name, new_name))
    {
      if (! GIMP_ROW_GET_CLASS (row)->name_edited (row, new_name))
        {
          gtk_widget_error_bell (GTK_WIDGET (row));
        }
    }

  gtk_popover_popdown (GTK_POPOVER (priv->popover));
}
