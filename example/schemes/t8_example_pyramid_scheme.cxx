/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

#include <t8_eclass.h>
#include <t8_element_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_bits.c>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid.h>

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_dpyramid_t       root;
  t8_dpyramid_root (&root);

  t8_linearidx_t      leafs_on_level =
    t8_dpyramid_num_descendants_at_leveldiff (&root, 12);
  t8_debugf ("leafs_on_level: %lu\n", leafs_on_level);
  t8_dpyramid_init_linear_id (&root, 12, leafs_on_level / 2);
  t8_dpyramid_debug_print (&root);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
