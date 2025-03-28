/* LIBGIMP - The GIMP Library
 * Copyright (C) 1995-2003 Peter Mattis and Spencer Kimball
 *
 * gimpdrawablefilter_pdb.c
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

/* NOTE: This file is auto-generated by pdbgen.pl */

#include "config.h"

#include "stamp-pdbgen.h"

#include "gimp.h"


/**
 * SECTION: gimpdrawablefilter
 * @title: gimpdrawablefilter
 * @short_description: Operations on drawable filters.
 *
 * Operations on drawable filters: creation, editing.
 **/


/**
 * gimp_drawable_filter_id_is_valid:
 * @filter_id: The filter ID to check.
 *
 * Returns %TRUE if the drawable filter ID is valid.
 *
 * This procedure checks if the given drawable filter ID is valid and
 * refers to an existing filter.
 *
 * Returns: Whether the filter ID is valid.
 *
 * Since: 3.0
 **/
gboolean
gimp_drawable_filter_id_is_valid (gint filter_id)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  gboolean valid = FALSE;

  args = gimp_value_array_new_from_types (NULL,
                                          G_TYPE_INT, filter_id,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-id-is-valid",
                                               args);
  gimp_value_array_unref (args);

  if (GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS)
    valid = GIMP_VALUES_GET_BOOLEAN (return_vals, 1);

  gimp_value_array_unref (return_vals);

  return valid;
}

/**
 * gimp_drawable_filter_new:
 * @drawable: The drawable.
 * @operation_name: The GEGL operation's name.
 * @name: (nullable): The effect name.
 *
 * Create a new drawable filter.
 *
 * This procedure creates a new filter for the specified operation on
 * @drawable.
 * The new effect still needs to be either added or merged to @drawable
 * later. Add the effect non-destructively with
 * [method@Gimp.Drawable.append_filter].
 * Currently only layers can have non-destructive effects. The effects
 * must be merged for all other types of drawable.
 *
 * Returns: (transfer none): The newly created filter.
 *
 * Since: 3.0
 **/
GimpDrawableFilter *
gimp_drawable_filter_new (GimpDrawable *drawable,
                          const gchar  *operation_name,
                          const gchar  *name)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  GimpDrawableFilter *filter = NULL;

  args = gimp_value_array_new_from_types (NULL,
                                          GIMP_TYPE_DRAWABLE, drawable,
                                          G_TYPE_STRING, operation_name,
                                          G_TYPE_STRING, name,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-new",
                                               args);
  gimp_value_array_unref (args);

  if (GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS)
    filter = GIMP_VALUES_GET_DRAWABLE_FILTER (return_vals, 1);

  gimp_value_array_unref (return_vals);

  return filter;
}

/**
 * gimp_drawable_filter_get_name:
 * @filter: The filter whose name you want.
 *
 * Get a drawable filter's name.
 *
 * This procedure returns the specified filter's name.
 * Since it is not possible to set a drawable filter's name yet, this
 * will be the operation's name. Eventually this filter's name will be
 * a free form field so do not rely on this information for any
 * processing.
 *
 * Returns: (transfer full): The filter's name.
 *          The returned value must be freed with g_free().
 *
 * Since: 3.0
 **/
gchar *
gimp_drawable_filter_get_name (GimpDrawableFilter *filter)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  gchar *name = NULL;

  args = gimp_value_array_new_from_types (NULL,
                                          GIMP_TYPE_DRAWABLE_FILTER, filter,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-get-name",
                                               args);
  gimp_value_array_unref (args);

  if (GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS)
    name = GIMP_VALUES_DUP_STRING (return_vals, 1);

  gimp_value_array_unref (return_vals);

  return name;
}

/**
 * gimp_drawable_filter_get_operation_name:
 * @filter: The filter whose operation name you want.
 *
 * Get a drawable filter's operation name.
 *
 * This procedure returns the specified filter's operation name.
 *
 * Returns: (transfer full): The filter's operation name.
 *          The returned value must be freed with g_free().
 *
 * Since: 3.0
 **/
gchar *
gimp_drawable_filter_get_operation_name (GimpDrawableFilter *filter)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  gchar *name = NULL;

  args = gimp_value_array_new_from_types (NULL,
                                          GIMP_TYPE_DRAWABLE_FILTER, filter,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-get-operation-name",
                                               args);
  gimp_value_array_unref (args);

  if (GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS)
    name = GIMP_VALUES_DUP_STRING (return_vals, 1);

  gimp_value_array_unref (return_vals);

  return name;
}

/**
 * gimp_drawable_filter_get_visible:
 * @filter: The filter.
 *
 * Get the visibility of the specified filter.
 *
 * This procedure returns the specified filter's visibility.
 *
 * Returns: The filter visibility.
 *
 * Since: 3.0
 **/
gboolean
gimp_drawable_filter_get_visible (GimpDrawableFilter *filter)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  gboolean visible = FALSE;

  args = gimp_value_array_new_from_types (NULL,
                                          GIMP_TYPE_DRAWABLE_FILTER, filter,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-get-visible",
                                               args);
  gimp_value_array_unref (args);

  if (GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS)
    visible = GIMP_VALUES_GET_BOOLEAN (return_vals, 1);

  gimp_value_array_unref (return_vals);

  return visible;
}

/**
 * gimp_drawable_filter_set_visible:
 * @filter: The filter.
 * @visible: The new filter visibility.
 *
 * Set the visibility of the specified filter.
 *
 * This procedure sets the specified filter's visibility.
 * The drawable won't be immediately rendered. Use
 * [method@Gimp.Drawable.update] to trigger an update.
 *
 * Returns: TRUE on success.
 *
 * Since: 3.0
 **/
gboolean
gimp_drawable_filter_set_visible (GimpDrawableFilter *filter,
                                  gboolean            visible)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  gboolean success = TRUE;

  args = gimp_value_array_new_from_types (NULL,
                                          GIMP_TYPE_DRAWABLE_FILTER, filter,
                                          G_TYPE_BOOLEAN, visible,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-set-visible",
                                               args);
  gimp_value_array_unref (args);

  success = GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS;

  gimp_value_array_unref (return_vals);

  return success;
}

/**
 * gimp_drawable_filter_get_opacity:
 * @filter: The filter.
 *
 * Get the opacity of the specified filter.
 *
 * This procedure returns the specified filter's opacity.
 *
 * Returns: The filter's opacity.
 *
 * Since: 3.0
 **/
gdouble
gimp_drawable_filter_get_opacity (GimpDrawableFilter *filter)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  gdouble opacity = 0.0;

  args = gimp_value_array_new_from_types (NULL,
                                          GIMP_TYPE_DRAWABLE_FILTER, filter,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-get-opacity",
                                               args);
  gimp_value_array_unref (args);

  if (GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS)
    opacity = GIMP_VALUES_GET_DOUBLE (return_vals, 1);

  gimp_value_array_unref (return_vals);

  return opacity;
}

/**
 * gimp_drawable_filter_get_blend_mode:
 * @filter: The filter.
 *
 * Get the blending mode of the specified filter.
 *
 * This procedure returns the specified filter's mode.
 *
 * Returns: The effect blending mode.
 *
 * Since: 3.0
 **/
GimpLayerMode
gimp_drawable_filter_get_blend_mode (GimpDrawableFilter *filter)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  GimpLayerMode mode = 0;

  args = gimp_value_array_new_from_types (NULL,
                                          GIMP_TYPE_DRAWABLE_FILTER, filter,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-get-blend-mode",
                                               args);
  gimp_value_array_unref (args);

  if (GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS)
    mode = GIMP_VALUES_GET_ENUM (return_vals, 1);

  gimp_value_array_unref (return_vals);

  return mode;
}

/**
 * _gimp_drawable_filter_update:
 * @filter: The filter.
 * @propnames: (array zero-terminated=1): Array of property names.
 * @propvalues: Array of values, one per property in propnames.
 * @opacity: The filter's opacity.
 * @blend_mode: The effect blending mode.
 * @blend_space: The effect blending space.
 * @composite_mode: The layer composite mode.
 * @composite_space: The effect composite space.
 * @auxinputnames: (array zero-terminated=1): Array of aux input pads.
 * @auxinputs: (element-type GimpDrawable) (array zero-terminated=1): Array of drawables, one per auxinputnames.
 *
 * Update the settings of the specified filter.
 *
 * This procedure updates the settings of the specified filter all at
 * once.
 * In particular, update will be frozen and will happen only once for
 * all changed settings.
 * This PDB function is internal, meant to be private and its arguments
 * will likely change as filters evolve. It should not be used.
 *
 * Returns: TRUE on success.
 *
 * Since: 3.0
 **/
gboolean
_gimp_drawable_filter_update (GimpDrawableFilter      *filter,
                              const gchar            **propnames,
                              const GimpValueArray    *propvalues,
                              gdouble                  opacity,
                              GimpLayerMode            blend_mode,
                              GimpLayerColorSpace      blend_space,
                              GimpLayerCompositeMode   composite_mode,
                              GimpLayerColorSpace      composite_space,
                              const gchar            **auxinputnames,
                              const GimpDrawable     **auxinputs)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  gboolean success = TRUE;

  args = gimp_value_array_new_from_types (NULL,
                                          GIMP_TYPE_DRAWABLE_FILTER, filter,
                                          G_TYPE_STRV, propnames,
                                          GIMP_TYPE_VALUE_ARRAY, propvalues,
                                          G_TYPE_DOUBLE, opacity,
                                          GIMP_TYPE_LAYER_MODE, blend_mode,
                                          GIMP_TYPE_LAYER_COLOR_SPACE, blend_space,
                                          GIMP_TYPE_LAYER_COMPOSITE_MODE, composite_mode,
                                          GIMP_TYPE_LAYER_COLOR_SPACE, composite_space,
                                          G_TYPE_STRV, auxinputnames,
                                          GIMP_TYPE_CORE_OBJECT_ARRAY, auxinputs,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-update",
                                               args);
  gimp_value_array_unref (args);

  success = GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS;

  gimp_value_array_unref (return_vals);

  return success;
}

/**
 * _gimp_drawable_filter_get_number_arguments:
 * @filter: The filter.
 *
 * Queries for the number of arguments on the specified filter.
 *
 * This procedure returns the number of arguments on the specified
 * filter.
 * For specific information on each input argument, use
 * gimp_drawable_filter_get_argument().
 *
 * Returns: The number of input arguments.
 *
 * Since: 3.0
 **/
gint
_gimp_drawable_filter_get_number_arguments (GimpDrawableFilter *filter)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  gint num_args = 0;

  args = gimp_value_array_new_from_types (NULL,
                                          GIMP_TYPE_DRAWABLE_FILTER, filter,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-get-number-arguments",
                                               args);
  gimp_value_array_unref (args);

  if (GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS)
    num_args = GIMP_VALUES_GET_INT (return_vals, 1);

  gimp_value_array_unref (return_vals);

  return num_args;
}

/**
 * _gimp_drawable_filter_get_pspec:
 * @filter: The filter.
 * @arg_num: The argument number.
 *
 * Queries for information on the specified filter's argument.
 *
 * This procedure returns the #GParamSpec of filter's argument.
 *
 * Returns: (transfer full): The GParamSpec of the argument.
 *          The returned value must be freed with g_param_spec_unref().
 *
 * Since: 3.0
 **/
GParamSpec *
_gimp_drawable_filter_get_pspec (GimpDrawableFilter *filter,
                                 gint                arg_num)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  GParamSpec *param_spec = NULL;

  args = gimp_value_array_new_from_types (NULL,
                                          GIMP_TYPE_DRAWABLE_FILTER, filter,
                                          G_TYPE_INT, arg_num,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-get-pspec",
                                               args);
  gimp_value_array_unref (args);

  if (GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS)
    param_spec = GIMP_VALUES_DUP_PARAM (return_vals, 1);

  gimp_value_array_unref (return_vals);

  return param_spec;
}

/**
 * _gimp_drawable_filter_get_arguments:
 * @filter: The filter.
 * @values: (out): The values of the arguments in same order.
 *
 * Returns the currently set filter arguments.
 *
 * This procedure returns the filter's arguments.
 *
 * Returns: (array zero-terminated=1) (transfer full):
 *          The names of the arguments.
 *          The returned value must be freed with g_strfreev().
 *
 * Since: 3.0
 **/
gchar **
_gimp_drawable_filter_get_arguments (GimpDrawableFilter  *filter,
                                     GimpValueArray     **values)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  gchar **argnames = NULL;

  args = gimp_value_array_new_from_types (NULL,
                                          GIMP_TYPE_DRAWABLE_FILTER, filter,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-get-arguments",
                                               args);
  gimp_value_array_unref (args);

  if (GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS)
    {
      argnames = GIMP_VALUES_DUP_STRV (return_vals, 1);
      *values = g_value_dup_boxed (gimp_value_array_index (return_vals, 2));
    }

  gimp_value_array_unref (return_vals);

  return argnames;
}

/**
 * gimp_drawable_filter_delete:
 * @filter: The filter to delete.
 *
 * Delete a drawable filter.
 *
 * This procedure deletes the specified filter. This must not be done
 * if the drawable whose this filter was applied to was already deleted
 * or if the drawable was already removed from the image.
 * Do not use anymore the @filter object after having deleted it.
 *
 * Returns: TRUE on success.
 *
 * Since: 3.0
 **/
gboolean
gimp_drawable_filter_delete (GimpDrawableFilter *filter)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  gboolean success = TRUE;

  args = gimp_value_array_new_from_types (NULL,
                                          GIMP_TYPE_DRAWABLE_FILTER, filter,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-drawable-filter-delete",
                                               args);
  gimp_value_array_unref (args);

  success = GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS;

  gimp_value_array_unref (return_vals);

  return success;
}
