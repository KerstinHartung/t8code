/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_gtest_geometry_handling.cxx
 * This file contains tests for the geometry handling of t8code.
 */

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_handler.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_zero.hxx>

/* In this file we collect tests for t8code's cmesh geometry module.
 * These tests are
 *  - geometry.geometry_linear:  Check that linear geometry has correct name and dimension.
 *  - geometry.geometry_zero:    Check that zero geometry has correct name and dimension.
 *  - test_geometry.cmesh_geometry: 
 *  - test_geometry.cmesh_geometry_unique: Check that we can access the geometry via the tree id if
 *                                   we only use one geometry and did not specify tree ids for it.
 *                                   In this case t8code should automatically associate this geometry to all trees.
 *  - test_geometry.geom_handler_register: Tests the geometry_handler register and find interface.
 */
/* TODO: 
  * - Add a test for the jacobian, as soon as its implemented in parameterized test geometry.cmesh_geometry_linear.
  */

template <typename T>
class geometry_handling: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {
    dim = GetParam ();
    geom = new T (dim);
  }
  int dim;
  T *geom;
};

using MyTypes = ::testing::Types<t8_geometry_linear, t8_geometry_linear_axis_aligned, t8_geometry_occ,
                                 t8_geometry_analytic, t8_geometry_zero>;
TYPED_TEST_SUITE (geometry_handling, MyTypes);

/* Check that the geometries for dimensions 0,1,2,3
 * has the correct name and dimension. */
TYPED_TEST (geometry_handling, geometry_name_and_handling)
{
  char name[BUFSIZ];
  snprintf (name, BUFSIZ, "t8_geom_linear_%i", dim);
  ASSERT_EQ (strcmp (linear_geom.t8_geom_get_name (), name), 0)
    << "Linear geometry of dim " << dim << "has wrong name. Expected " << name << " got "
    << linear_geom.t8_geom_get_name ();
  ASSERT_EQ (dim, linear_geom.t8_geom_get_dimension ())
    << "Linear geometry of dim " << dim << "has wrong dimension: " << linear_geom.t8_geom_get_dimension () << ".";
}

TEST (test_geometry, cmesh_geometry)
{
  t8_cmesh_t cmesh;

  t8_debugf ("Testing cmesh tree geometry set/get.\n");

  /* Build a simple 2 tree cmesh and set geometries for the trees. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  /* Register the linear geometry and zero geometry to this cmesh. */
  auto linear_geom = t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 2);
  ;
  auto zero_geom = t8_cmesh_register_geometry<t8_geometry_zero> (cmesh, 2);
  /* Set the id geometry for the trees. */
  t8_cmesh_set_tree_geometry (cmesh, 0, linear_geom);
  t8_cmesh_set_tree_geometry (cmesh, 1, zero_geom);
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  /* Check that we can get the geometry back over the tree id. */
  const t8_geometry *found_geom = t8_cmesh_get_tree_geometry (cmesh, 0);
  ASSERT_EQ (found_geom->t8_geom_get_hash (), linear_geom->t8_geom_get_hash ())
    << "Could not find linear tree geometry at tree 0.";

  found_geom = t8_cmesh_get_tree_geometry (cmesh, 1);
  ASSERT_EQ (found_geom->t8_geom_get_hash (), zero_geom->t8_geom_get_hash ())
    << "Could not find linear tree geometry at tree 1.";
  /* clean-up */
  t8_cmesh_destroy (&cmesh);
}

TEST (test_geometry, cmesh_geometry_unique)
{
  t8_cmesh_t cmesh;

  t8_debugf ("Testing cmesh tree geometry get with unique geometry.\n");

  /* Build a simple 2 tree cmesh and set geometries for the trees. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  /* Register the linear_geometry to this cmesh. */
  auto provided_geom = t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 2);
  ;
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  /* Check that we can get the geometry back over the tree id.
   * This must now work even though we did not register the geometry for 
   * this tree. Since we only have one geometry. */
  auto found_geom = t8_cmesh_get_tree_geometry (cmesh, 0);
  ASSERT_TRUE (found_geom != nullptr) << "Could not find any geometry.";
  ASSERT_EQ (found_geom->t8_geom_get_hash (), provided_geom->t8_geom_get_hash ())
    << "Could not find cmesh tree geometry.";

  /* clean-up */
  t8_cmesh_destroy (&cmesh);
}

TEST (test_geometry, geom_handler_register)
{
  t8_geometry_handler geom_handler;
  const t8_geometry *found_geom;

  t8_debugf ("Testing geometry handler register.\n");

  /* For each dimension build the zero geometry and register it.
   * We then commit the handler and check that we can find the geometries. */
  for (int idim = 0; idim <= 3; ++idim) {
    /* Register the geometry. */
    geom_handler.register_geometry<t8_geometry_zero> (idim);
  }

  /* Check find geometry. */
  for (int idim = 0; idim < 3; ++idim) {
    t8_geometry_zero zero_geom (idim);
    std::string name;

    /* Get the name of this geometry. */
    name = zero_geom.t8_geom_get_name ();

    t8_debugf ("Name of geometry: %s.\n", name);

    /* Find the geometry by name. */
    found_geom = geom_handler.get_geometry (name);
    ASSERT_TRUE (found_geom != NULL) << "No geometry found.";
    ASSERT_EQ (found_geom->t8_geom_get_name (), name) << "Could not find identity geometry.";
  }
  /* Try to find a different geometry via the name. Must return nullptr. */
  found_geom = geom_handler.get_geometry ("random_name34823412414");
  ASSERT_TRUE (found_geom == nullptr) << "Found a geometry that should not exist.";

  /* Try to find a different geometry via the hash. Must return nullptr. */
  found_geom = geom_handler.get_geometry (std::hash<std::string> {}("random_name34823412414"));
  ASSERT_TRUE (found_geom == nullptr) << "Found a geometry that should not exist.";
}
