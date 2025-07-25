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

#pragma once

#include "gimpcloneoptions.h"


#define GIMP_TYPE_PERSPECTIVE_CLONE_OPTIONS            (gimp_perspective_clone_options_get_type ())
#define GIMP_PERSPECTIVE_CLONE_OPTIONS(obj)            (G_TYPE_CHECK_INSTANCE_CAST ((obj), GIMP_TYPE_PERSPECTIVE_CLONE_OPTIONS, GimpPerspectiveCloneOptions))
#define GIMP_PERSPECTIVE_CLONE_OPTIONS_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST ((klass), GIMP_TYPE_PERSPECTIVE_CLONE_OPTIONS, GimpPerspectiveCloneOptionsClass))
#define GIMP_IS_PERSPECTIVE_CLONE_OPTIONS(obj)         (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GIMP_TYPE_PERSPECTIVE_CLONE_OPTIONS))
#define GIMP_IS_PERSPECTIVE_CLONE_OPTIONS_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), GIMP_TYPE_PERSPECTIVE_CLONE_OPTIONS))
#define GIMP_PERSPECTIVE_CLONE_OPTIONS_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj), GIMP_TYPE_PERSPECTIVE_CLONE_OPTIONS, GimpPerspectiveCloneOptionsClass))


typedef struct _GimpPerspectiveCloneOptionsClass GimpPerspectiveCloneOptionsClass;

struct _GimpPerspectiveCloneOptions
{
  GimpCloneOptions          paint_instance;

  GimpPerspectiveCloneMode  clone_mode;
};

struct _GimpPerspectiveCloneOptionsClass
{
  GimpCloneOptionsClass  parent_class;
};


GType   gimp_perspective_clone_options_get_type (void) G_GNUC_CONST;
