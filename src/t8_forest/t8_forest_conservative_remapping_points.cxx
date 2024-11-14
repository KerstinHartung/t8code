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

/** file t8_forest_vtk.h
 */

/* TODO: Document this file */
// qKH and cKH to indicate my questions and comments

#ifndef T8_FOREST_CONSERVATIVE_REMAPPING
#define T8_FOREST_CONSERVATIVE_REMAPPING

#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_element_cxx.hxx>
#include <vector>

/*
 Approach: 
 - generate points within old elements
 - loop through new elements
   - and search for old elements inside (ensure that not counting double)
   - and calculate weights of old elements and sum to get new value
*/

T8_EXTERN_C_BEGIN ();

void
t8_next_element(int index_in, int &index_out, int direction, t8_element_shape_t element_shape);

/* Search query, a set of points. Stores information of new cell and intersecting old cells */
// 
typedef struct
{
  std::vector<t8_locidx_t> intersection_cell_indices;
  std::vector<t8_locidx_t> intersection_cell_treeid;
  double volume_cell;
  double volume_intersection_cells;
  double value_cell;

} t8_cell_new_t;

/* For old elements store:
   - list of points that exist within (and if they have been accounted for yet)
   - volume of cell
   - volume that is accounted for in new cells already (should be equal to volume_cell in the end)
*/
typedef struct
{
 std::vector<double*> coord_point;
 std::vector<int*> point_counted;
 t8_locidx_t cell_treeid;
 double volume_cell;
 double volume_cell_counted;
} t8_cell_old_t;

/* Additional user data that we process during search.
 * For each element we count the number of particles that it contains
 * and we count the total number of elements that we constructed during search. */
typedef struct
{
  double volume_cell;

  //t8_cell_corners_t *corners_per_cell;   /* List of each cell including the corners. */
  //Liste (Anzahl Zellen Zielgitter) von dynamischen structuren fuer die Zell-IDs des Ursprungsgitters

} t8_search_user_data_t;

/*
 This routine should generate points equally inside an element
 Missing:
  - approach to generate points efficiently depending on shape of element 
*/
static void
t8_forest_get_element_points_inside (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                 std::vector<double*> &out_coords)
{
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
  element_shape = scheme->t8_element_shape (element);
  int ipoint;
  for( int isize = 0; isize < t8_element_corner_ref_coords[element_shape][0]; isize++ )
  {
    for (int jsize = 0; jsize < t8_element_corner_ref_coords[element_shape][1]; jsize++)
    {
      for (int ksize = 0; ksize < t8_element_corner_ref_coords[element_shape][2]; ksize++)
      {
        // fill vector
        const double *ref_coords = t8_element_corner_ref_coords[element_shape][idim];
        // TODO: add check if this new point is inside the respective element - otherwise continue to next iteration
        ipoint++;
        out_coords.push_back(0); 
        out_coords.at(ipoint) = new double[3];
        // add corner elements, function requires pre-existing memory
        t8_forest_element_from_ref_coords (forest, ltreeid, element, ref_coords, 1, out_coords.at(ipoint), NULL );
      }
    }
  }
}

/*
 * The search callback.
 * It will be called once per element and generally decides whether or not
 * to continue the search with the children of the element.
 * Since we will continue as long as there are points left to consider,
 * we always return 1 here.
 * The search will then only stop when no queries are active (thus, no points
 * could be in this element) or the element is a leaf element.
 */
static int
t8_search_points_element_callback (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, const int is_leaf,
                             t8_element_array_t *leaf_elements, t8_locidx_t tree_leaf_index, void *query,
                             size_t query_index)
{
  T8_ASSERT (query == NULL);

  /* Get a pointer to our user data. */
  t8_search_user_data_t *user_data = (t8_search_user_data_t *) t8_forest_get_user_data (forest);
  // todoKH
  // enable assertion of user data again
  //T8_ASSERT (user_data != NULL);
  return 1;
}

/* The query callback. This will be called for each element once per active query
 * (= particle that may be inside this element).
 * The return value determines whether or not the query remains active for the children
 * of this element.
 * In our example this is the case if the particle is inside the element.
 * The additional input parameter 'is_leaf' will be true if the given element is a leaf
 * element (= an actual element in our forest, not a 'virtual' element in the hierarchy).
 * If the element is a leaf and the particle is contained in it, then we will increase
 * a counter for this element by one.
 * These counters are provided in an sc_array as user data of the input forest.
 */
// Note: forest will be forest_new
static int
t8_search_points_query_callback (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                   const int is_leaf, t8_element_array_t *leaf_elements, t8_locidx_t tree_leaf_index,
                                   void *query, size_t query_index)
{
  int point_is_inside_element;
  t8_cell_new_t *build_new_cell = (t8_cell_new_t *) query;
 // T8_ASSERT (user_data != NULL);
  T8_ASSERT (query != NULL);
  /* Numerical tolerance for the is_inside_element check. */
  const double tolerance = 1e-8;
  t8_search_user_data_t *user_data = (t8_search_user_data_t *) t8_forest_get_user_data (forest);
  /* Ensure user_data is present. */
  T8_ASSERT (user_data != NULL);

  // search for old points in new element, i.e. forest parameters below refers to forest_new
  // loop through old element points
  int points_counted = 0;
  for( size_t ipoint = 0; ipoint < points_of_cell->coord_point.size(); ipoint++ )
  {
    point_is_inside_element = t8_forest_element_point_inside (forest, ltreeid, element, 
        points_of_cell->coord_point.at(ipoint), tolerance);
    if (point_is_inside_element && points_of_cell->point_counted.at(ipoint)==0) {
      if (is_leaf) {
        t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
        //add index of the cell to list of new element
        build_new_cell->intersection_cell_indices.push_back( element_index );
        build_new_cell->intersection_cell_treeid.push_back( ltreeid );
        // store that this point was counted now
        points_of_cell->point_counted.at(ipoint) = 1;
        point_of_cell->point_cou
        // stop loop now to avoid adding the same element several times
      }
      // TODO what to do here?? where should "return 1" be, now that points are very much independent
      /* The particles is inside the element (finding one corner is sufficient). This query should remain active.
       * If this element is not a leaf the search will continue with its children. */
      return 1;
    }
  }
  if (!point_is_inside_element){
    /* The particle is not inside the element. Deactivate this query.
     * If no active queries are left, the search will stop for this element and its children. */
    return 0;
  }
}

/* For old forest - build points that span elements */
static sc_array *
t8_build_points( t8_forest_t forest )
{
  int ielem;
  double test[3];
  t8_locidx_t itree;
  sc_array *points = sc_array_new_count (sizeof (t8_cell_old_t), t8_forest_get_global_num_elements( forest ));
  std::cout << t8_forest_get_num_local_trees (forest) <<"  "<<t8_forest_get_global_num_elements( forest )<<std::endl;
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
    //cKH: Returns a pointer to an array element indexed by a plain int. (note on sc_array_index_int)
      t8_cell_old_t *points
        = (t8_cell_old_t *) sc_array_index_int (points, ielem);
      t8_locidx_t looptree=itree;
      auto element = t8_forest_get_element ( forest, ielem_tree, &looptree);
      //Generate points inside of an element
      t8_forest_get_element_points_inside(forest, looptree, element, points->coord_point);
      points->volume_cell = t8_forest_element_volume( forest, looptree, element );
    }
  }
  return points;
}

void t8_free_points( t8_forest_t forest, sc_array *points)
{
  int ielem;
  t8_locidx_t itree;
  std::cout << t8_forest_get_num_local_trees (forest) <<"  "<<t8_forest_get_global_num_elements( forest )<<std::endl;
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
      t8_cell_corners_t *corner
        = (t8_cell_corners_t *) sc_array_index_int (corners, ielem);
      t8_locidx_t looptree=itree;
      auto element = t8_forest_get_element ( forest, ielem_tree, &looptree);
      const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, looptree);
      const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
      t8_element_shape_t element_shape = scheme->t8_element_shape (element);
      int num_corner = t8_eclass_num_vertices[element_shape];
      for( int icorner = 0; icorner < num_corner; icorner++ ){
         delete corner->coordinates.at(icorner);
      }
    }
  }
}

void
t8_forest_conservative_remapping_planar( t8_forest_t forest_old, t8_forest_t forest_new )
{
  int ielem;
  t8_element_t *ielem_t;
  t8_locidx_t itree;
  //0) find corners of elements of new forest in old forest elements
  sc_array *corners = t8_build_corners(forest_new);
  printf("corners built\n");
  ielem=0;
  t8_cell_corners_t *corner= (t8_cell_corners_t *) sc_array_index_int (corners, ielem);
  std::cout << points->coordinates.at(2)[1]<<"\n";
  //1) find elements of old forest which contain corners of new elements (via corners)
  //   the search already loops through all elements of new forest
  t8_forest_search(forest_new, t8_search_points_element_callback, t8_search_points_query_callback, points);
  //2) iterate throught all new cells
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest_new); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest_new, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    // use later on again!! for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
    for (t8_locidx_t ielem_tree = 0; ielem_tree < 2; ielem_tree++, ielem++) {
      t8_cell_corners_t* corner
        = (t8_cell_corners_t *) sc_array_index_int (corners, ielem);
      // ( t8_forest_t forest1, t8_forest_t forest2, sc_array *corners, t8_element_t *elem2 )
      // cKH: I changed the arguments here - instead of elem1 (elem of the old forest) the
      //      corners information is passed (or one element of corners)
      //3) determine intersection (based on indices of old elements containing new corners))
      // qKH: what is difference ielem and ielem_t
      ielem_t = (t8_element_t *)  sc_array_index (corners, ielem);
      t8_locidx_t looptree=itree;
      auto element = t8_forest_get_element ( forest_new, ielem_tree, &looptree);
      // qKH: one cell_intersection function for all types, or varying?
      std::cout<<"start intersection\n";
      cell_intersection(forest_old, forest_new, corner, element, itree);
      //4) Compare volume
      // add volume as output argument
      //point_cloud_to_volume(points, ...);
    }
  }

  // cKH
  // the next few lines can probably be integrated in the above nested loop
  //loop over trees and number of elements of new forest
      // sum volume and check if volume_new_remap = volume_new
      // i.e. if volume_cell and volume_intersection_cells are the same

  // continue only for new elements which don't fulfill condition, i.e. volume_new=volume_old_sum
  // need to ensure that elements are not counted twice!!!
  // check also which other steps might be covered by tolerance
        //5) Schnittpunkte von Kanten berechnen
        //6) Zellen des alten Forests suchen, in dem die Kanten enthalten sind
            // find old elements which contain point cloud points of new element
            // and which were not already stored with flag 1
            // store those indices with flag 2

        //7) Schnitt der Zellen bilden
            // for all points in point cloud:
            //   search for neighbour of old cell (or go through all cells)
            //   i) that contains point and where one corner of new element contained
            //   ii) search for intersection of edges (if i) no fulfilled)

        //8) Volumen vergleichen

          //9) Eckpunkte der alten Zellen in der neuen Zelle suchen
          //10) Zellen suchen, die diese Eckpunkte als Eckpunkte haben
             // find for remaining cells - revert search and check if old completely part of
             //   new element
             // ensure that cells not already captured as type 1 or 2

          //11) Schnitt der Zellen bilden

          //12) Volumen vergleichen -> nicht gleich? -> Fehler

    //13) Berechnen des Wertes der neuen Zelle

   //14) free memory
  t8_free_corners( forest_new, corners);

}

T8_EXTERN_C_END ();

/* Volumen berechnen */
//

#endif /* !T8_FOREST_CONSERVATIVE_REMAPPING */
