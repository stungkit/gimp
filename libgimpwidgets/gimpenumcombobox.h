/* LIBGIMP - The GIMP Library
 * Copyright (C) 1995-1997 Peter Mattis and Spencer Kimball
 *
 * gimpenumcombobox.h
 * Copyright (C) 2004  Sven Neumann <sven@gimp.org>
 *
 * This library is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library.  If not, see
 * <https://www.gnu.org/licenses/>.
 */

#pragma once

#if !defined (__GIMP_WIDGETS_H_INSIDE__) && !defined (GIMP_WIDGETS_COMPILATION)
#error "Only <libgimpwidgets/gimpwidgets.h> can be included directly."
#endif

#include <libgimpwidgets/gimpintcombobox.h>

G_BEGIN_DECLS

#define GIMP_TYPE_ENUM_COMBO_BOX (gimp_enum_combo_box_get_type ())
G_DECLARE_DERIVABLE_TYPE (GimpEnumComboBox,
                          gimp_enum_combo_box,
                          GIMP, ENUM_COMBO_BOX,
                          GimpIntComboBox)


struct _GimpEnumComboBoxClass
{
  GimpIntComboBoxClass  parent_class;

  /* Padding for future expansion */
  void (* _gimp_reserved0) (void);
  void (* _gimp_reserved1) (void);
  void (* _gimp_reserved2) (void);
  void (* _gimp_reserved3) (void);
  void (* _gimp_reserved4) (void);
  void (* _gimp_reserved5) (void);
  void (* _gimp_reserved6) (void);
  void (* _gimp_reserved7) (void);
  void (* _gimp_reserved8) (void);
  void (* _gimp_reserved9) (void);
};


GtkWidget * gimp_enum_combo_box_new              (GType             enum_type);
GtkWidget * gimp_enum_combo_box_new_with_model   (GimpEnumStore    *enum_store);

void        gimp_enum_combo_box_set_icon_prefix  (GimpEnumComboBox *combo_box,
                                                  const gchar      *icon_prefix);

G_END_DECLS
