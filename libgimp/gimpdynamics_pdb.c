/* LIBGIMP - The GIMP Library
 * Copyright (C) 1995-2003 Peter Mattis and Spencer Kimball
 *
 * gimpdynamics_pdb.c
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
 * SECTION: gimpdynamics
 * @title: gimpdynamics
 * @short_description: Operations related to paint dynamics.
 *
 * Operations related to paint dynamics.
 **/


/**
 * gimp_dynamics_refresh:
 *
 * Refresh current paint dynamics. This function always succeeds.
 *
 * This procedure retrieves all paint dynamics currently in the user's
 * paint dynamics path and updates the paint dynamics dialogs
 * accordingly.
 *
 * Returns: TRUE on success.
 *
 * Since: 2.8
 **/
gboolean
gimp_dynamics_refresh (void)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  gboolean success = TRUE;

  args = gimp_value_array_new_from_types (NULL,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-dynamics-refresh",
                                               args);
  gimp_value_array_unref (args);

  success = GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS;

  gimp_value_array_unref (return_vals);

  return success;
}

/**
 * gimp_dynamics_get_name_list:
 * @filter: (nullable): An optional regular expression used to filter the list.
 *
 * Retrieve the list of loaded paint dynamics.
 *
 * This procedure returns a list of the paint dynamics that are
 * currently available.
 *
 * Returns: (array zero-terminated=1) (transfer full):
 *          The list of paint dynamics names.
 *          The returned value must be freed with g_strfreev().
 *
 * Since: 2.8
 **/
gchar **
gimp_dynamics_get_name_list (const gchar *filter)
{
  GimpValueArray *args;
  GimpValueArray *return_vals;
  gchar **dynamics_list = NULL;

  args = gimp_value_array_new_from_types (NULL,
                                          G_TYPE_STRING, filter,
                                          G_TYPE_NONE);

  return_vals = _gimp_pdb_run_procedure_array (gimp_get_pdb (),
                                               "gimp-dynamics-get-name-list",
                                               args);
  gimp_value_array_unref (args);

  if (GIMP_VALUES_GET_ENUM (return_vals, 0) == GIMP_PDB_SUCCESS)
    dynamics_list = GIMP_VALUES_DUP_STRV (return_vals, 1);

  gimp_value_array_unref (return_vals);

  return dynamics_list;
}
