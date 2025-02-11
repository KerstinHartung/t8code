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
#include <t8_element.hxx>
#include <vector>

/*
 Approach:
 - generate points within old elements
   - applying the batch-based search it seems easiest to store all points in one long list, including their
     metadata (i.e. element of origin in old forest
   - for now decided against just creating general grid of data, to allow even volume distribution within elements (although
     this is not really the case yet, only for hexahedra))
 - loop through new elements
   - and search for old elements inside (ensure that not counting double)
   - and calculate weights of old elements and sum to get new value
*/

/* TODO
   - for now just one value of user data - advance so that several can be processed at once?
     or is it better to do this only fully after the search, based on stored volume contributions (i.e. weights)
     and indices of old elements that make contributions? - probably more efficient to do it like that
     (although this means more calculations are required)

*/

T8_EXTERN_C_BEGIN ();

/* Search query, a set of points. Stores information of old cells
   - list of points that exist within (and if they have been accounted for yet)
   - volume of cell
   - volume that is accounted for in new cells already (should be equal to volume_cell in the end)
*/
typedef struct
{
 std::vector<double*> coord_point;
 std::vector<int> point_counted;
 t8_locidx_t cell_treeid;
 t8_locidx_t cell_index;
 double volume_cell;
 double volume_cell_part;
} t8_cell_old_t;

// Additional user data that we process during search.
// Everything related to the new mesh, i.e. its volume and volume re-created from old elements
// Potential also the value calculated for the new field, although this potentially is done separate to the search
// The list of old elements associated
typedef struct
{
  double volume_cell;
  double volume_cell_remapped;
  double value;
  std::vector<t8_locidx_t> element_index;
  std::vector<t8_locidx_t> element_treeid;
} t8_search_user_data_t;

/*
 This routine should generate points equally inside an element (here from old forest)
 Missing:
  - approach to generate points efficiently depending on shape of element
*/
static void
t8_forest_get_element_points_inside (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                 int number_points, int ipoint,
                                 std::vector<double*> &out_coords)
{
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
  t8_element_shape_t element_shape = scheme->t8_element_shape (element);
  int num_corner = t8_eclass_num_vertices[element_shape];
  int *point_is_inside_element;
  const double tolerance = 1e-8;
  out_coords.push_back(0);
  ipoint = 0;
  double step = 1/number_points;
  out_coords.at(ipoint) = new double[3];
  for( int isize = 0; isize < t8_element_corner_ref_coords[element_shape][num_corner][0]*number_points; isize++ )
  {
    for (int jsize = 0; jsize < t8_element_corner_ref_coords[element_shape][num_corner][1]*number_points; jsize++)
    {
      for (int ksize = 0; ksize < t8_element_corner_ref_coords[element_shape][num_corner][2]*number_points; ksize++)
      {
        // fill vector
        //const double *ref_coords = t8_element_corner_ref_coords[element_shape][T8_ECLASS_MAX_CORNERS][idim];
        const double ref_coords[3] = {step*(0.5+isize),(0.5+jsize)*step,(0.5+ksize)*step};
        // add corner elements, function requires pre-existing memory
        t8_forest_element_from_ref_coords (forest, ltreeid, element, ref_coords, 1, out_coords.at(ipoint));
        t8_forest_element_points_inside (forest, ltreeid, element,
                out_coords.at(ipoint), 1, point_is_inside_element, tolerance);
        if (point_is_inside_element){
          ipoint++;
          out_coords.at(ipoint) = new double[3];
        }
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
t8_search_points_element_callback (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                   const int is_leaf, const t8_element_array_t *leaf_elements,
                                   const t8_locidx_t tree_leaf_index)
{
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
// Note: this was updated to t8code v3, but only to compile, not yet checked for reasonability
static void
t8_search_points_query_callback (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                   const int is_leaf, const t8_element_array_t *leaf_elements,
                                   const t8_locidx_t tree_leaf_index, sc_array_t *queries, sc_array_t *query_indices,
                                   int *query_matches, const size_t num_active_queries)
                                   //void *query, size_t query_index)
{
  //t8_cell_old_t *points_of_cell = (t8_cell_old_t *) queries;
  //T8_ASSERT (query != NULL);
  /* Build an array of all point-coords, used for t8_forest_element_point_batch_inside */
  double *coords = T8_ALLOC (double, 3 * num_active_queries);
  for (size_t point_iter = 0; point_iter < num_active_queries; point_iter++) {
    /* Get the query at the current query-index (point_iter in this case). */
    // pointer_iter=point_id only in the beginning, will change once number of queries reduces
    const size_t point_id = *(size_t *) sc_array_index_int (query_indices, point_iter);
    /* Cast the query into a particle*/
    t8_cell_old_t *point
      = (t8_cell_old_t *) sc_array_index ((sc_array_t *) queries, point_id);
    /* extract the coordinates of the particle struct */
    coords[3 * point_iter] = point->coord_point.at(0)[0];
    coords[3 * point_iter + 1] = point->coord_point.at(0)[1];
    coords[3 * point_iter + 2] = point->coord_point.at(0)[2];
  }
  t8_search_user_data_t *build_new_cell = (t8_search_user_data_t *) t8_forest_get_user_data (forest);
  T8_ASSERT (build_new_cell != NULL);
  /* Numerical tolerance for the is_inside_element check. */
  const double tolerance = 1e-8;

  // search for old points in new element, i.e. forest parameters below refers to forest_new
  // loop through old element points
//  int no_points_inside = true;
//  for( size_t ipoint = 0; ipoint < points_of_cell->coord_point.size(); ipoint++ )
//  {
//    t8_forest_element_points_inside (forest, ltreeid, element,
//        points_of_cell->coord_point.at(ipoint), 1, num_active_queries, tolerance);
//    // point of old element is inside new and not already counted
//    if (point_is_inside_element && points_of_cell->point_counted.at(ipoint)==0) {
//      if (is_leaf) {
//        //add index of the old element to list of new element, to calculate new value later on
//        build_new_cell->element_index.push_back( points_of_cell->cell_index );
//        build_new_cell->element_treeid.push_back( points_of_cell->cell_treeid );
//        //add volume, as correctness test
//        build_new_cell->volume_cell_remapped += points_of_cell->volume_cell_part;
//        // store that this point was counted now
//        points_of_cell->point_counted.at(ipoint) = 1;
//      }
//      no_points_inside = false;
//      // TODO what to do here?? where should "return 1" be, now that points are very much independent
//      /* The particles is inside the element (finding one corner is sufficient). This query should remain active.
//       * If this element is not a leaf the search will continue with its children. */
//  //    return 1;
//    }
//  }
  //if (no_points_inside || sum(points_of_cell->point_counted) == points_of_cell->coord_point.size()){
    /* No points of this element are inside the new element or all were already matches. Deactivate this query.
     * If no active queries are left, the search will stop for this element and its children. */
 //   return 0;
 // }
}

/* For old forest - build points that span elements */
static sc_array *
t8_build_points( t8_forest_t forest, int nr_points )
{
  int ielem;
  int ipoint;
  double test[3];
  t8_locidx_t itree;
  sc_array *points = sc_array_new_count (sizeof (t8_cell_old_t), t8_forest_get_global_num_elements( forest ));
  std::cout << t8_forest_get_num_local_trees (forest) <<"  "<<t8_forest_get_global_num_elements( forest )<<std::endl;
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
    //cKH: Returns a pointer to an array element indexed by a plain int. (note on sc_array_index_int)
      t8_cell_old_t *point
        = (t8_cell_old_t *) sc_array_index_int (points, ielem);
      t8_locidx_t looptree=itree;
      auto element = t8_forest_get_element ( forest, ielem_tree, &looptree);
      //Generate points inside of an element
      t8_forest_get_element_points_inside(forest, looptree, element, nr_points, ipoint, point->coord_point);
      point->cell_treeid = itree;
      point->cell_index = ielem_tree;
      point->volume_cell = t8_forest_element_volume( forest, looptree, element );
      point->volume_cell_part = point->volume_cell/ipoint;
      point->point_counted.resize(ipoint);
    }
  }
  return points;
}

void t8_free_points( t8_forest_t forest, sc_array *points, int nr_points)
{
  t8_locidx_t itree;
  int ielem;
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest, itree);
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
      t8_cell_old_t *point
        = (t8_cell_old_t *) sc_array_index_int (points, ielem);
      for( int ipoint = 0; ipoint < nr_points; ipoint++ ){
         //todo: not all points actually allocated, need additional check
         delete point->coord_point.at(ipoint);
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
  int nr_points;
  //0) construct poinst in old elements
  sc_array *points = t8_build_points(forest_old, nr_points);
  printf("points built\n");
  // debugging
  ielem=0;
  t8_cell_old_t *point= (t8_cell_old_t *) sc_array_index_int (points, ielem);
  std::cout << point->coord_point.at(2)[1]<<"\n";
  //1) find elements of new forest which contain points in old elements
  //   the search already loops through all elements of new forest
  t8_forest_search(forest_new, t8_search_points_element_callback, t8_search_points_query_callback, points);
 ////2) iterate throught all new cells
 //for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest_new); itree++) {
 //  const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest_new, itree);
 //  /* Inner loop: Iteration over the elements of the local tree */
 //  // use later on again!! for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
 //  for (t8_locidx_t ielem_tree = 0; ielem_tree < 2; ielem_tree++, ielem++) {
 //    t8_cell_corners_t* corner
 //      = (t8_cell_corners_t *) sc_array_index_int (corners, ielem);
 //    // ( t8_forest_t forest1, t8_forest_t forest2, sc_array *corners, t8_element_t *elem2 )
 //    // cKH: I changed the arguments here - instead of elem1 (elem of the old forest) the
 //    //      corners information is passed (or one element of corners)
 //    //3) determine intersection (based on indices of old elements containing new corners))
 //    // qKH: what is difference ielem and ielem_t
 //    ielem_t = (t8_element_t *)  sc_array_index (corners, ielem);
 //    t8_locidx_t looptree=itree;
 //    auto element = t8_forest_get_element ( forest_new, ielem_tree, &looptree);
 //    // qKH: one cell_intersection function for all types, or varying?
 //    std::cout<<"start intersection\n";
 //    cell_intersection(forest_old, forest_new, corner, element, itree);
 //    //4) Compare volume
 //    // add volume as output argument
 //    //point_cloud_to_volume(points, ...);
 //  }
 //}
  t8_free_points( forest_new, points, nr_points);
}

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_CONSERVATIVE_REMAPPING */
