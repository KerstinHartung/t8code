#include <t8_eclass.h>
#include <t8_cmesh.h>
#include "t8_cmesh/t8_cmesh_trees.h"
#include <t8_cmesh_readmshfile.h>

static int
t8_test_supported_msh_file (t8_cmesh_t cmesh)
{
  int                 retval;
  int                 read_node = 1;
  int                 check_neigh_elem = 1;
  int                 check_element_type;
  double             *vertices;
  t8_gloidx_t         num_gtree;
  t8_locidx_t         ltree_id;
  t8_locidx_t         ltree_it;
  t8_locidx_t         lnum_trees;
  t8_eclass_t         tree_class;

  /* Description of the properties of the example msh-files. */
  const int           number_elements = 4;
  const int           dimension = 2;
  const t8_eclass_t   elem_type = T8_ECLASS_TRIANGLE;
  /* *INDENT-OFF* */
  int vertex[6][2] = {
                       {0, 0},
                       {2, 0},
                       {4, 0},
                       {1, 2},
                       {3, 2},
                       {2, 4} };

  int elements[4][3] = {
                        {0, 1, 3},
                        {1, 4, 3},
                        {1, 2, 4},
                        {3, 4, 5} };

  int face_neigh_elem[4][3] = {
                                {1, -1,-1},
                                {3, 0, 2},
                                {-1, 1, -1},
                                {-1, -1, 1} };
  /* *INDENT-ON* */

  if (cmesh == NULL) {
    /* If the cmesh is NULL. */
    return 0;
  }
  else {
    /* Checks if the cmesh was comitted. */
    retval = t8_cmesh_is_committed (cmesh);
    SC_CHECK_ABORT (retval == 1, "Cmesh commit failed.");
    /* Checks for face consistency. */
    retval = t8_cmesh_trees_is_face_consistend (cmesh, cmesh->trees);
    SC_CHECK_ABORT (retval == 1, "Cmesh face consistency failed.");
    /* Checks if the number of elements was read correctly. */
    num_gtree = t8_cmesh_get_num_trees (cmesh);
    SC_CHECK_ABORT (num_gtree == number_elements,
                    "Number of elements in msh-file was read incorrectly.");
    /* Number of local trees. */
    lnum_trees = t8_cmesh_get_num_local_trees (cmesh);
    /* Iterate through the local elements and check if they were read properly. */
    /* All trees should be local to the master rank. */
    for (ltree_it = 0; ltree_it < lnum_trees; ltree_it++) {
      tree_class = t8_cmesh_get_tree_class (cmesh, ltree_it);
      check_element_type = t8_eclass_compare (tree_class, elem_type);
      SC_CHECK_ABORT (check_element_type == 0,
                      "Element type in msh-file was read incorrectly.");
      /* Get pointer to the vertices of the tree. */
      vertices = t8_cmesh_get_tree_vertices (cmesh, ltree_it);
      /* Checking the msh-files elements and nodes. */
      for (int i = 0; i < 3; i++) {
        /* Checks if x and y coordinate of the nodes are not read correctly. */
        if (!((vertex[elements[ltree_it][i]][0] == (int) vertices[3 * i])
              && (vertex[elements[ltree_it][i]][1] ==
                  (int) vertices[(3 * i) + 1]))) {
          read_node = 0;
          SC_CHECK_ABORT (read_node == 1, "Node was read incorrectly.");
          /* Return error code, if the nodes are not read correctly. */
          return -1;
        }
        /* Checks whether the face neighbor elements are not read correctly. */
        ltree_id =
          t8_cmesh_get_face_neighbor (cmesh, ltree_it, i, NULL, NULL);
        if (!(ltree_id == face_neigh_elem[ltree_it][i])) {
          check_neigh_elem = 0;
          SC_CHECK_ABORT (check_neigh_elem == 1,
                          "The face neigbhor element in the example test file was not read correctly.");
          /* Return error code, if the face neighbor elements are not read correctly. */
          return -1;
        }
      }
    }
    /* If the checks were performed correctly. */
    return 1;
  }
}

int
main (int argc, char **argv)
{

  int                 mpiret, retval;
  t8_cmesh_t          cmesh;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

/* The mesh is unpartitioned and all trees are local to the master rank, no other rank has any trees. */

  t8_global_productionf ("Testing: Reading of msh-files.\n");

  /* Testing supported msh-file version 2 (as example). */
  cmesh =
    t8_cmesh_from_msh_file ("test/test_msh_file_vers2_ascii", 1,
                            sc_MPI_COMM_WORLD, 2, 0);
  retval = t8_test_supported_msh_file (cmesh);
  if (retval != 1) {
    t8_global_errorf ("Reading of supported msh-file version failed.\n");
    goto end_test;
  }

#if 1
  /* Testing unsupported version of msh-files (bin-format). */
  cmesh =
    t8_cmesh_from_msh_file ("test/test_msh_file_vers2_bin", 1,
                            sc_MPI_COMM_WORLD, 2, 0);
  retval = t8_test_supported_msh_file (cmesh);
  if (retval != 0) {
    t8_global_errorf ("Rejecting of unsupported msh-files failed.\n");
    goto end_test;
  }

  /* Testing unsupported version of msh-files (version 4). */
  cmesh =
    t8_cmesh_from_msh_file ("test/test_msh_file_vers4_ascii", 1,
                            sc_MPI_COMM_WORLD, 2, 0);
  retval = t8_test_supported_msh_file (cmesh);
  if (retval != 0) {
    t8_global_errorf ("Rejecting of unsupported msh-files failed.\n");
    goto end_test;
  }
#endif
  t8_global_productionf ("Done testing reading of msh-files\n");

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
/* If the test was succesful. */
  return 0;

/* In case anything failed, the test is going to be closed. */
end_test:
  t8_debugf ("The test was not succesful.\n");
  if (t8_cmesh_is_committed (cmesh)) {
    t8_cmesh_destroy (&cmesh);
  }
  return 1;
}
