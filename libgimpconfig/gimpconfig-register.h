/* LIBGIMP - The GIMP Library
 * Copyright (C) 1995-2003 Peter Mattis and Spencer Kimball
 *
 * gimpconfig-register.c
 * Copyright (C) 2008-2019 Michael Natterer <mitch@gimp.org>
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

#if !defined (__GIMP_CONFIG_H_INSIDE__) && !defined (GIMP_CONFIG_COMPILATION)
#error "Only <libgimpconfig/gimpconfig.h> can be included directly."
#endif

G_BEGIN_DECLS

/* For information look into the C source or the html documentation */


GType   gimp_config_type_register (GType         parent_type,
                                   const gchar  *type_name,
                                   GParamSpec  **pspecs,
                                   gint          n_pspecs);

G_END_DECLS
