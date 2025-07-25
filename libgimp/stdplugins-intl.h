/* LIBGIMP - The GIMP Library
 * Copyright (C) 1995-1997 Peter Mattis and Spencer Kimball
 *
 * stdplugins-intl.h
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

#ifndef GETTEXT_PACKAGE
#error "config.h must be included prior to stdplugins-intl.h"
#endif

#include <glib/gi18n.h>


#define DEFINE_STD_SET_I18N                                          \
static gboolean                                                      \
set_i18n (GimpPlugIn   *plug_in,                                     \
          const gchar  *procedure_name,                              \
          gchar       **gettext_domain,                              \
          gchar       **catalog_dir)                                 \
{                                                                    \
  *gettext_domain = g_strdup (GETTEXT_PACKAGE"-std-plug-ins");       \
  return TRUE;                                                       \
};

#define STD_SET_I18N set_i18n
