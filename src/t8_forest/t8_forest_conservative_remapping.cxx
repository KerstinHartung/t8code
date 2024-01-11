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
#include <vector>

T8_EXTERN_C_BEGIN ();

/* Search query, a set of points. Stores information of new cell and intersecting old cells */
typedef struct
{ 
  std::vector<double[3]> coordinates;/* The corners of our cell. */
  int num_corners;
  std::vector<t8_locidx_t> *intersection_cell_indices;
  std::vector<t8_locidx_t> *intersection_cell_treeid;
  double volume_cell;
  double volume_intersection_cells;

} t8_cell_corners_t;

/* Additional user data that we process during search.
 * For each element we count the number of particles that it contains
 * and we count the total number of elements that we constructed during search. */
typedef struct
{
  t8_cell_corners_t *corners_per_cell;   /* List of each cell including the corners. */
  //Liste (Anzahl Zellen Zielgitter) von dynamischen structuren fuer die Zell-IDs des Ursprungsgitters
 
} t8_search_user_data_t;

/*
 * The search callback.
 * It will be called once per element and generally decides whether or not
 * to continue the search with the children of the element.
 * Since we will continue as long as there are corners left to consider,
 * we always return 1 here.
 * The search will then only stop when no queries are active (thus, no corners
 * could be in this element) or the element is a leaf element.
 */
static int
t8_search_callback (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, const int is_leaf,
                             t8_element_array_t *leaf_elements, t8_locidx_t tree_leaf_index, void *query,
                             size_t query_index)
{
  T8_ASSERT (query == NULL);

  /* Get a pointer to our user data. */
  t8_search_user_data_t *user_data = (t8_search_user_data_t *) t8_forest_get_user_data (forest);
  T8_ASSERT (user_data != NULL);
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
static int
t8_tutorial_search_query_callback (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                   const int is_leaf, t8_element_array_t *leaf_elements, t8_locidx_t tree_leaf_index,
                                   void *query, size_t query_index)
{
  int corner_is_inside_element;
  t8_cell_corners_t *corners_of_cell = (t8_cell_corners_t *) query;
  //T8_ASSERT (user_data != NULL);
  T8_ASSERT (query != NULL);
  /* Numerical tolerance for the is_inside_element check. */
  const double tolerance = 1e-8;

  for( int icorner = 0; icorner < corners_of_cell->num_corners ; icorner++ )
  {
    corner_is_inside_element = t8_forest_element_point_inside (forest, ltreeid, element, corners_of_cell->coordinates[icorner], tolerance);
    if (corner_is_inside_element) {
      if (is_leaf) {
        t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
        //add index of the cell to 
        // qKH: what is missing in the above comment?
        corners_of_cell->intersection_cell_indices->push_back( element_index );
        corners_of_cell->intersection_cell_treeid->push_back( ltreeid );
      }
      /* The particles is inside the element. This query should remain active.
       * If this element is not a leaf the search will continue with its children. */
      return 1;
    }
    /* The particle is not inside the element. Deactivate this query.
     * If no active queries are left, the search will stop for this element and its children. */
    return 0;
  }
}

// qKH
/*
 * forest1 describes the original grid
 * forest2 describes the new grid, on which data should be remapped
 */
// qKH: I added "corners" as argument, or is this available another way?
void //rueckgabewert anpassen oder als pointer uebergeben - Punktwolke
cell_intersection( t8_forest_t forest1, t8_forest_t forest2, sc_array *corners, t8_element_t *elem2 )
{
  // 
  //t8_forest_element_coordinate (t8_forest_t forest, t8_locidx_t ltree_id, const t8_element_t *element, int corner_number,     
  //                              double *coordinates) 
  // following usage in t8_forest_element_diam
  t8_eclass_t tree_class;                                                                                                       
  t8_eclass_scheme_c *ts;
 
  double coordinates[3]; 
  int num_corners;

  // cKH
  /* Approach
   * - start at first corner of forest1 (old) - obtained from search for points
   * - store this coordinate as first member of point cloud
   * - check in direction of next corner (pick order randomly): crossing edge of forest2?
   *    - check by calculating t8_vertex_point_inside between all neighboring corners of forest2
   *        if next (or do we only have "one of the points"?) point of corners inside => no
   *    - if yes: the crossing point is next member of point cloud
             - continue along the edge of elem2 to next corner of elem2
             - check if corner inside of elem1 (t8_forest_element_point_inside)
                 - if yes: this is the next member of the point cloud
                 - if no: go in direction of other boardering corner of elem2, add this to point cloud
   *    - if no: the next corner of elem1 is the next member of point cloud
   * - continue either on edge of elem1 or elem2 until next corner and then repeat previous step
   * - stop when arrived at first corner again
   * - pass point cloud list 
   */

  /* cKH: some old stuff
  tree_class = t8_forest_get_tree_class (forest2, ltreeid);
  ts = t8_forest_get_eclass_scheme (forest2, tree_class);
  T8_ASSERT (ts->t8_element_is_valid (elem2));
  num_corners = ts->t8_element_num_corners (elem2);

  for (i = 0; i < num_corners; i++) {
    //(t8_forest_t forest,t8_locidx_t ltree_id,const t8_element_t *element,int corner_number,double *coordinates)
    t8_forest_element_coordinate (forest2, ltreeid, elem2, i, coordinates);
  }
  */

  // cKH: alternative approach assuming struct "corners" containing information for one new cell: which 
  //      old cells are intersected
  /*
   - get coordinates of first element of corners_of_cell->intersection_cell_indices and add to point cloud
   - then follow the steps listed above, always calculate coordinates for new elements 
   - fill volume_cell (based on t8_forest_element_volume)
  */

// cKH: potentially useful functions to determine next member of point cloud
// t8_forest_element_line_length: distance between two corners
/* t8_triangle_point_inside: const double p_0[3], const double v[3], const double w[3], const double point[3],
 *                          const double tolerance)
 * \param[in] v           The vector from p_0 to p_1 (second vertex in the triangle)
 * \param[in] w           The vector from p_0 to p_2 (third vertex in the triangle)
*/


//ltree_id in t8_cell_corners_t -> uebergeben
// qKH: but doesn't "corners" already store the ltree_id value after the callback (for the old forest)?

//Eckpunkte des Schnittes zurueckgeben

}

/*
 * We are looking at the intersection between two convex cells. 
 * Thus, the result is again a convex shape. We can build a Tet-Mesh out of these
 * points. Therefore, we chose one point and search for the three points, that are closest
 * to this point. If these are not planar, they are the first tet.
 * From each surface we now search for the point, that is closest to this. The triangle and 
 * the closest point again are a tet. We have to test, if the tet intersects another
 * cell. Then it is no tet of the mesh and the tested triangle is an outer surface and has no
 * neighbor.
 */
void 
point_cloud_to_volume( std::vector<double[3]> points )
{
  // triangulation https://www.kiv.zcu.cz/site/documents/verejne/vyzkum/publikace/technicke-zpravy/2002/tr-2002-02.pdf
  // qKH: Use pre-existing library like qhull for this? (needs to be available for commerical use)
}

// qKH
/*
 - Is t8_build_corners complete? I guess I am not sure what its purpose is.
 - A "corners" array is created but never filled? 
 - "corner" is filled with some indices but not returned from the function? What am I missing?
 - are num_corners and coordinates filled here (or should be)?
 - answer: yes, return is missing here
 - Chiara: set num_corners to "corners" struct here
 */
static sc_array *
t8_build_corners( t8_forest_t forest )
{
  int itree, ielem;
  sc_array *corners = sc_array_new_count (sizeof (t8_cell_corners_t), t8_forest_get_global_num_elements( forest ));
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
      t8_cell_corners_t *corner
        = (t8_cell_corners_t *) sc_array_index_int (corners, ielem);
    }
  }
}

void
t8_forest_conservative_remapping_planar( t8_forest_t forest_old, t8_forest_t forest_new )
{
  int itree, ielem;
  sc_array *corners = t8_build_corners(forest_new);
  //1) Iterieren ueber Zellen des neuen Forests
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest_new); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest_new, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
      //2) Zellen des alten Forests suchen, in dem die Eckpunkte der jeweiligen Zelle enthalten sind
      t8_forest_search(forest_old, t8_search_callback, t8_tutorial_search_query_callback, corners );
      // ( t8_forest_t forest1, t8_forest_t forest2, sc_array *corners, t8_element_t *elem2 )
      // qKH: does it make sense to have this in one loop or is it possible to access a previous
      //      ielem in the corners structure? antyhing against merging these associated calls
      //      which cycle over new forest's elements?
      // cKH: I changed the arguments here - instead of elem1 (elem of the old forest) the 
      //      corners information is passed
      cell_intersection(forest_old, forest_new, corners, ielem_tree)
      // add volume as output argument
      point_cloud_to_volume(points, ...)
    }
  }
  // cKH
  // the next few lines can probably be integrated in the above nested loop
  // volume_new_remap = 0 (can this be stored in volume_intersection_cells?)
  //loop over trees and number of elements of new forest
      //volume_new = t8_forest_element_volume (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element)
      //for elements of corners (i.e. old elements that intersect with new forest)
          // calculate intersection of corners and forest_new elements
          // calculate volume 
      // sum volume and check if volume_new_remap = volume_new 
      // i.e. if volume_cell and volume_intersection_cells are the same?


      //3) Schnitt der Zellen bilden

      //4) Volumen vergleichen

        //5) Schnittpunkte von Kanten berechnen

        //6) Zellen des alten Forests suchen, in dem die Schnittpunkte enthalten sind

        //7) Schnitt der Zellen bilden

        //8) Volumen vergleichen

          //9) Eckpunkte der alten Zellen in der neuen Zelle suchen

          //10) Zellen suchen, die diese Eckpunkte als Eckpunkte haben

          //11) Schnitt der Zellen bilden

          //12) Volumen vergleichen -> nicht gleich? -> Fehler
  
    //13) Berechnen des Wertes der neuen Zelle  
}

T8_EXTERN_C_END ();

/* Volumen berechnen */
//

#endif /* !T8_FOREST_CONSERVATIVE_REMAPPING */
