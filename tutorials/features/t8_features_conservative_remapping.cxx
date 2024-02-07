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
                                                                                                                                
/* See also: https://github.com/holke/t8code/wiki/Step-2---Creating-a-uniform-forest                                            
 *                                                                                                                              
 * After we learned how to create a cmesh in step1, we will                                                                     
 * now build our first partitioned forest, get its local and global                                                             
 * element count, and output it into .vtu files.                                                                                
 *                                                                                                                              
 * When we create a forest from a coarse mesh, the forest will always be                                                        
 * uniform (every element has the same refinement level) and can then be adapted                                                
 * later (see the following steps).                                                                                             
 * Together with the cmesh, we also need a refinement scheme. This scheme tells the                                             
 * forest how elements of each shape (t8_eclass_t) are refined, what their neighbor                                             
 * are etc.                                                                                                                     
 * The default scheme in t8_schemes/t8_default/t8_default_cxx.hxx provides an implementation for                                
 * all element shapes that t8code supports (with pyramids currently under construction).                                        
 *                                                                                                                              
 * How you can experiment here:                                                                                                 
 *  - Use Paraview to visualize the output files.                                                                               
 *  - Execute this program with different numbers of processes.                                                                 
 *  - Change the initial refinement level.                                                                                      
 *  - Use a different cmesh (See step1).                                                                                        
 *  - Look into t8_forest.h and try to get different information about the                                                      
 *    forest (for example the number of local trees).                                                                           
 */                                                                                                                             
                                                                        
#include <t8.h>                                     /* General t8code header, always include this. */
#include <t8_cmesh.h>                               /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>             /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>            /* forest definition and general interface. */
#include <t8_forest/t8_forest_io.h>                 /* forest io interface. */
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* default refinement scheme. */
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h> /* linear geometry of the cmesh */
#include <t8_forest/t8_forest_conservative_remapping.cxx>

/* Builds cmesh of 2 prisms that build up a unit cube. 
 * See step1 for a detailed description.
 * \param [in] comm   MPI Communicator to use.
 * \return            The coarse mesh.
 */
static t8_cmesh_t
t8_step2_build_prismcube_coarse_mesh (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;

  /* Build a coarse mesh of 2 prism trees that form a cube. */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_PRISM, comm, 0, 0, 0);
  t8_global_productionf (" [step2] Constructed coarse mesh with 2 prism trees.\n");

  return cmesh;
}


/* Definition of a two dimensional mesh with linear geometry and periodic boundaries.
 * The mesh consists of four triangles and two quads.
 *
 *  This is how the cmesh looks like. The numbers are the tree numbers:
 *
 *   +---+---+
 *   |   |5 /|
 *   | 3 | / |
 *   |   |/ 4|
 *   +---+---+
 *   |1 /|   |
 *   | / | 2 |
 *   |/0 |   |
 *   +---+---+
 *  */
t8_cmesh_t
t8_cmesh_new_periodic_hybrid_2d (sc_MPI_Comm comm)
{

  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[60] = { 0,   0,   0,                                        /* tree 0, triangle */
                          0.5, 0,   0, 0.5, 0.5, 0, 0,   0,   0,              /* tree 1, triangle */
                          0.5, 0.5, 0, 0,   0.5, 0, 0.5, 0,   0,              /* tree 2, quad */
                          1,   0,   0, 0.5, 0.5, 0, 1,   0.5, 0, 0,   0.5, 0, /* tree 3, quad */
                          0.5, 0.5, 0, 0,   1,   0, 0.5, 1,   0, 0.5, 0.5, 0, /* tree 4, triangle */
                          1,   0.5, 0, 1,   1,   0, 0.5, 0.5, 0,              /* tree 5, triangle */
                          1,   1,   0, 0.5, 1,   0 };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_geometry_c *linear_geom = t8_geometry_linear_new (2);
  t8_cmesh_register_geometry (cmesh, linear_geom); /* Use linear geometry */

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 5, T8_ECLASS_TRIANGLE);

  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 4);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 30, 4);
  t8_cmesh_set_tree_vertices (cmesh, 4, vertices + 42, 3);
  t8_cmesh_set_tree_vertices (cmesh, 5, vertices + 51, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
  t8_cmesh_set_join (cmesh, 0, 2, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 0, 3, 2, 3, 0);

  t8_cmesh_set_join (cmesh, 1, 3, 0, 2, 1);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 1, 0);

  t8_cmesh_set_join (cmesh, 2, 4, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 5, 2, 0, 1);

  t8_cmesh_set_join (cmesh, 3, 5, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 0, 0, 0);

  t8_cmesh_set_join (cmesh, 4, 5, 1, 2, 0);

  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

/* Build a uniform forest on a cmesh 
 * using the default refinement scheme.
 * \param [in] comm   MPI Communicator to use.
 * \param [in] cmesh  The coarse mesh to use.
 * \param [in] level  The initial uniform refinement level.
 * \return            A uniform forest with the given refinement level that is
 *                    partitioned across the processes in \a comm.
 */
static t8_forest_t
t8_step2_build_uniform_forest (sc_MPI_Comm comm, t8_cmesh_t cmesh, int level)
{
  t8_forest_t forest;
  t8_scheme_cxx_t *scheme;

  /* Create the refinement scheme. */
  scheme = t8_scheme_new_default_cxx ();
  /* Creat the uniform forest. */
  forest = t8_forest_new_uniform (cmesh, scheme, level, 0, comm);

  return forest;
}

/* Write vtk (or more accurately vtu) files of the forest.
 * \param [in] forest   A forest.
 * \param [in] prefix   A string that is used as a prefix of the output files.
 * 
 * This will create the file prefix.pvtu
 * and additionally one file prefix_MPIRANK.vtu per MPI rank.
 */
static void
t8_step2_write_forest_vtk (t8_forest_t forest, const char *prefix)
{
  t8_forest_write_vtk (forest, prefix);
}

/* Destroy a forest. This will free all allocated memory.
 * \param [in] forest    A forest.
 * NOTE: This will also free the memory of the scheme and the cmesh, since
 *       the forest took ownership of them.
 *       If we do not want this behaviour, but want to reuse for example the cmesh,
 *       we need to call t8_cmesh_ref (cmesh) before passing it to t8_forest_new_uniform.
 */
static void
t8_step2_destroy_forest (t8_forest_t forest)
{
  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  int mpiret;
  sc_MPI_Comm comm;
  t8_cmesh_t cmesh_old, cmesh_new;
  t8_forest_t forest_old, forest_new;
  /* The prefix for our output files. */
  const char *prefix = "t8_step2_uniform_forest";
  /* The uniform refinement level of the forest. */
  const int level_old = 3, level_new = 2;
  t8_locidx_t local_num_elements;
  t8_gloidx_t global_num_elements;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* Print a message on the root process. */
  t8_global_productionf (" [step2] \n");
  t8_global_productionf (" [step2] Hello, this is the step2 example of t8code.\n");
  t8_global_productionf (" [step2] In this example we build our first uniform forest and output it to vtu files.\n");
  t8_global_productionf (" [step2] \n");

  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;
  /* Create the cmesh from step1 */
  cmesh_old = t8_step2_build_prismcube_coarse_mesh (comm);
  cmesh_new = t8_cmesh_new_periodic_hybrid_2d (comm);

  /* Build the uniform forest, it is automatically partitioned among the processes. */
  forest_old = t8_step2_build_uniform_forest (comm, cmesh_old, level_old);
  forest_new = t8_step2_build_uniform_forest (comm, cmesh_new, level_new);

  /* Perform planar conservative remapping. */
  t8_forest_conservative_remapping_planar(forest_old, forest_new);

  /* Get the local number of elements. */
  local_num_elements = t8_forest_get_local_num_elements (forest_old);
  /* Get the global number of elements. */
  global_num_elements = t8_forest_get_global_num_elements (forest_old);

  /* Print information on the forest. */
  t8_global_productionf (" [step2] Created uniform forest.\n");
  t8_global_productionf (" [step2] Refinement level:\t\t\t%i\n", level_old);
  t8_global_productionf (" [step2] Local number of elements:\t\t%i\n", local_num_elements);
  t8_global_productionf (" [step2] Global number of elements:\t%li\n", global_num_elements);

  /* Write forest to vtu files. */
  t8_step2_write_forest_vtk (forest_old, prefix);
  t8_global_productionf (" [step2] Wrote forest to vtu files:\t%s*\n", prefix);

  /* Destroy the forest. */
  t8_step2_destroy_forest (forest_old);
  t8_global_productionf (" [step2] Destroyed forest.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
