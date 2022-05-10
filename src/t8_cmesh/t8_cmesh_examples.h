/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

/**
 * \file t8_cmesh_examples.h
 * We provide example coarse meshes in this file 
 */

#ifndef T8_CMESH_EXAMPLES
#define T8_CMESH_EXAMPLES
#include <t8_cmesh.h>

T8_EXTERN_C_BEGIN ();

/** Constructs a cmesh from a given p4est_connectivity structure.
 * \param[in]       conn       The p4est connectivity.
 * \param[in]       comm       mpi communicator to be used with the new cmesh.
 * \param[in]       do_partition Flag whether the cmesh should be partitioned or not.
 * \return          A t8_cmesh structure that holds the same connectivity information
 *                  as \a conn.
 */
t8_cmesh_t          t8_cmesh_new_from_p4est (p4est_connectivity_t * conn,
                                             sc_MPI_Comm comm,
                                             int do_partition);

/** Constructs a cmesh from a given p8est_connectivity structure.
 * \param[in]       conn       The p8est connectivity.
 * \param[in]       comm       mpi communicator to be used with the new cmesh.
 * \param[in]       do_dup     Flag whether the communicator shall be duplicated or not.
 * \param[in]       do_partition Flag whether the cmesh should be partitioned or not.
 * \return          A t8_cmesh structure that holds the same connectivity information
 *                  as \a conn.
 */
t8_cmesh_t          t8_cmesh_new_from_p8est (p8est_connectivity_t * conn,
                                             sc_MPI_Comm comm,
                                             int do_partition);

/* TODO: it could possibly be a problem that we do not set the dimension of
 * the cmesh. This could i.e. be difficult when we combine an empty cmesh with
 * a non-empty one. */
/** Construct a cmesh that has no trees. We do not know a special use case,
 * this function is merely for debugging and to show the possibility.
 * \param [in]      comm       mpi communicator to be used with the new cmesh.
 * \param [in]      do_partition Flag whether the cmesh should be partitioned or not.
 * \param [in]      dimension  An empty cmesh requires a dimension nevertheless 0 <= \a dimension <= 4.
 * \return                     A committed t8_cmesh structure that has no trees.
 */
t8_cmesh_t          t8_cmesh_new_empty (sc_MPI_Comm comm, int do_partition,
                                        int dimension);

/** Constructs a cmesh that consists only of one tree of a given element class.
 * \param [in]      eclass     The element class.
 * \param [in]      comm       mpi communicator to be used with the new cmesh.
 * \param [in]      do_dup     Flag whether the communicator shall be duplicated or not.
 * \return          A committed t8_cmesh structure with one tree of class \a eclass.
 */
t8_cmesh_t          t8_cmesh_new_from_class (t8_eclass_t eclass,
                                             sc_MPI_Comm comm);

/** Construct a hypercube forest from one primitive tree class.
 * \param [in] eclass       This element class determines the dimension and
 *                          the number of trees needed to construct a cube.
 * \param [in] comm         The mpi communicator to be used.
 * \param [in] do_bcast     If this flag is nonzero the cmesh is only constructed
 *                          on processor 0 and then broadcasted to the other
 *                          processors in \a comm.
 *                          TODO: this parameter will be moved to internal.
 * \param [in] do_partition Create a partitioned cmesh.
 * \param [in] periodic     If true, the coarse mesh will be periodic in each direction.
 *                          Not possible with \a eclass pyramid.
 * TODO: Add periodic flags for each dimension.
 */
t8_cmesh_t          t8_cmesh_new_hypercube (t8_eclass_t eclass,
                                            sc_MPI_Comm comm,
                                            int do_bcast, int do_partition,
                                            int periodic);

/** Hybercube with 6 Tets, 6 Prism, 4 Hex. */
/* TODO: Document */
t8_cmesh_t          t8_cmesh_new_hypercube_hybrid (sc_MPI_Comm comm,
                                                   int do_partition,
                                                   int periodic);

/** Construct a unit interval/square/cube coarse mesh that is periodic in each direction.
 * Element class?
 * Hypercube?
 * TODO: redundant, remove.
 * \param [in] comm         The mpi communicator to use.
 * \param [in] dim          The dimension of the forest, 1, 2 or 3.
 * \return                  A valid cmesh, as if _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_periodic (sc_MPI_Comm comm, int dim);

/** Construct a unit square of two triangles that is periodic in x and y.
 * \param [in] comm         The mpi communicator to use.
 * \return                  A valid cmesh, as if _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_periodic_tri (sc_MPI_Comm comm);

/** Construct a unit square of two quads and four triangles that is periodic in x and y.
 * \param [in] comm         The mpi communicator to use.
 * \return                  A valid cmesh, as if _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_periodic_hybrid (sc_MPI_Comm comm);

/** Construct a unit interval coarse mesh that consists of 3 trees and is
 * periodic.
 * \param [in] comm         The mpi communicator to use.
 * \return                  A valid cmesh, as is _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_periodic_line_more_trees (sc_MPI_Comm comm);

/** Construct a mesh consisting of a given number of same type trees.
 * \param [in] eclass       This element class determines the dimension and
 *                          the type trees used.
 * \param [in] num_trees    The number of trees to use.
 * \param [in] comm         The MPI_Communicator used to commit the cmesh.
 * \return                  A valid cmesh, as if _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_bigmesh (t8_eclass_t eclass, int num_trees,
                                          sc_MPI_Comm comm);

/** Construct a forest of three connected askew lines
  * \param [in] comm         The mpi communicator to use.
  * \return                  A valid cmesh, as if _init and _commit had been called.
  */
t8_cmesh_t          t8_cmesh_new_line_zigzag (sc_MPI_Comm comm);

/** Construct a forest of num_of_prisms connected prism, all with one edge in 0,
  * except for num_of_prisms = 2, then the return is the hypercube mesh
  * \param [in] comm        The mpi communicator to use.
  * \param [in] num_of_prisms The number of prisms to be used.
  * \return                 A valid cmesh, as if _init and _commit had been called.
  */
t8_cmesh_t          t8_cmesh_new_prism_cake (sc_MPI_Comm comm,
                                             int num_of_prisms);

/** Construct a single deformed prism
  * \param [in] comm        The mpi communicator to use.
  * \return                 A valid cmesh; as if _init and _commit had been called.*/
t8_cmesh_t          t8_cmesh_new_prism_deformed (sc_MPI_Comm comm);

/** Construct a single deformed pyramid
 * \param [in] comm       The mpi communicator to use.
 * \return                 A valid cmesh; as if _init and _commit had been called.*/
t8_cmesh_t          t8_cmesh_new_pyramid_deformed (sc_MPI_Comm comm);

/** Construct a forest of six connected noncannoical oriented prisms
  * \param [in] comm        The mpi communicator to use.
  * \return                 A valid cmesh, as if _init and _commit had been called.
  */
t8_cmesh_t          t8_cmesh_new_prism_cake_funny_oriented (sc_MPI_Comm comm);

/** Construct a forest of six connected noncannoical oriented prisms
  * \param [in] comm        The mpi communicator to use.
  * \return                 A valid cmesh, as if _init and _commit had been called.
  */
t8_cmesh_t          t8_cmesh_new_prism_geometry (sc_MPI_Comm comm);

/** Create a partitoned cmesh of quads whose local trees are given by an
 * num_x by num_y brick connectivity from p4est
 * or a num_x by num_y by num_z brick connectivity from p8est.
 * num_x and num_y and num_z can be different for different MPI ranks.
 * \param [in] num_x       The number of trees in x direction for this rank. Must be >= 0.
 * \param [in] num_y       The number of trees in y direction for this rank. Must be >= 0.
 * \param [in] num_y       The number of trees in z direction for this rank. Must be >= 0.
 *                         If nonzero, the cmesh is 3 dimensional.
 * \param [in] x_periodic  If nonzero, the local brick connectivity is periodic in x direction.
 * \param [in] y_periodic  If nonzero, the local brick connectivity is periodic in y direction.
 * \param [in] y_periodic  If nonzero and \a num_z > 0, the local brick connectivity is periodic in z direction.
 * \param [in] comm        The MPI communicator used to commit the cmesh.
 * \return                 A committed and partitioned cmesh. The process local trees
 *                         form a \a num_x by \a num_y (by \a num_z) brick.
 * It is possible for num_x or num_y to be set to zero. In this case the local part
 * of the cmesh will be empty.
 * If num_z is set to zero, the cmesh is 2 dimensional.
 */
t8_cmesh_t          t8_cmesh_new_disjoint_bricks (t8_gloidx_t num_x,
                                                  t8_gloidx_t num_y,
                                                  t8_gloidx_t num_z,
                                                  int x_periodic,
                                                  int y_periodic,
                                                  int z_periodic,
                                                  sc_MPI_Comm comm);

/** Construct a tetrahedral cmesh that has all possible face to face
 * connections and orientations.
 * This cmesh is used for testing and debugging.
 * \param [in] comm        The MPI communicator used to commit the cmesh.
 * \return                 A committed and replicated cmesh of 24 tetrahedron trees
 *                         in which each (face -> face, orientation) face connection
 *                         is set at least once.
 *                         Note that most faces in this cmesh are boundary faces.
 */
t8_cmesh_t          t8_cmesh_new_tet_orientation_test (sc_MPI_Comm comm);

/** Construct a hybrid cmesh with 2 tets, 2 prism, 1 hex.
 * This cmesh is used for testing and debugging.
 * \param [in] comm        The MPI communicator used to commit the cmesh.
 * \return                 A committed and replicated hybrid cmesh of 5 trees.
 */
t8_cmesh_t          t8_cmesh_new_hybrid_gate (sc_MPI_Comm comm);

/** Construct a hybrid cmesh with 2 tets, 2 prism, 1 hex and all are deformed.
 * This cmesh is used for testing and debugging.
 * \param [in] comm        The MPI communicator used to commit the cmesh.
 * \return                 A committed and replicated hybrid cmesh of 5 trees.
 */
t8_cmesh_t          t8_cmesh_new_hybrid_gate_deformed (sc_MPI_Comm comm);

/** Construct a full hybrig cmesh, with 1 hex, 1 pyra, 1 prism and 1 tet
 * This cmesh is used for testing and debugging.
 * \param [in] comm        The MPI communicator used to commit the cmesh.
 * \return                 A committed and replicated hybrid cmesh of 4 trees.
 */
t8_cmesh_t          t8_cmesh_new_full_hybrid (sc_MPI_Comm comm);

/** Construct a mesh out of num_of_pyra many pyramids. They form a circle, face 0 is
 * connected with face 1 of the next pyramid.
 * \param [in] comm         The MPI communicator used to commit the cmesh
 * \param [in] num_of_pyra  The number of pyramids to construct. Should be larger than 2
 * \return                  A cmesh with num_of_pyra many pyramids
*/
t8_cmesh_t          t8_cmesh_new_pyramid_cake (sc_MPI_Comm comm,
                                               int num_of_pyra);

/** Construct a bigger mesh, consisting of many cubes made by pyramids
 * \param [in] comm         The MPI communicator used to commit the cmesh
 * \param [in] num_cubes    The number of cubes of pyramids
 * return                   A cmesh with \a num_cubes many hypercubes*/
t8_cmesh_t          t8_cmesh_new_long_brick_pyramid (sc_MPI_Comm comm,
                                                     int num_cubes);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_EXAMPLES */
