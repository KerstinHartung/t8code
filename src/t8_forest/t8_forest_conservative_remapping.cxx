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
#include <t8_element_cxx.hxx>
#include <vector>

T8_EXTERN_C_BEGIN ();

/* Search query, a set of points. Stores information of new cell and intersecting old cells */
typedef struct
{
  std::vector<double*> coordinates;/* The corners of our cell. */
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

static void
t8_forest_get_element_nodes (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                 std::vector<double *> out_coords)
{
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
  const t8_element_shape_t element_shape = scheme->t8_element_shape (element);
  for( size_t icorner = 0; icorner < element_shape; icorner++ )
  {
    const double *ref_coords = t8_element_corner_ref_coords[element_shape][icorner];
    t8_forest_element_from_ref_coords (forest, ltreeid, element, ref_coords, 1, out_coords[icorner], NULL );
  }
}

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

  for( int icorner = 0; icorner < corners_of_cell->coordinates.size(); icorner++ )
  {
    corner_is_inside_element = t8_forest_element_point_inside (forest, ltreeid, element, corners_of_cell->coordinates[icorner], tolerance);
    if (corner_is_inside_element) {
      if (is_leaf) {
        t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
        //add index of the cell to
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

void //rueckgabewert anpassen oder als pointer uebergeben - Punktwolke
cell_intersection( t8_forest_t forest_old, t8_forest_t forest_new, t8_cell_corners_t *corner,
        t8_element_t *elem_new, t8_locidx_t treeid_new )
{
  // information about old forest
  double coordinates[3];
  // todo update to variable number of corners
  int order_corner_old[] = {0,1,3,2};
  int order_corner_new[] = {0,1,3,2};
  const double tolerance = 1e-8;
  std::vector<double*> point_cloud;

  // qKH
  /*
    - why is the number of corners equal to coordinates.size()?
   */


  // cKH
  /* Approach
   * - start at first corner of forest_old
   * - iterate over corners and check if inside of new element
   *   - if yes: add to point_cloud and continue with next forest_old element's corner
   *   - if  no (and already found one point inside): look for cross section point between new and old element
   *      - iterate over corners of new element untill cross section found
   *      - add cross section point to point_cloud
   *      - check in which direction next point of new forest element of old element
   *      - check for further corners of new element in old element and add to point cloud
   *      - once no longer inside of old element: revert cross section point check (i.e.
   *        cycle through old element corners until intersection with current new element corners found)
   *   - if  no: continue
   * - stop when arrived at first corner again
   * - pass point cloud list
   * - fill volume_cell (based on t8_forest_element_volume)
   * - add check for parallel orientation
  */


  //number of cells of old forest within new forest element
  forest_old_nr_elem = corners->intersection_cell_indices.size();
  //number of corners of new forest element
  ncorner         = corners->coordinates.size();

  //todo: add loop over forest_old_nr_elem (independent outside)
  int_ind = 0  // index to iterate through old elements intersection with one new element

  //Calculate the corner elements of element of the current old tree
  t8_forest_get_element_nodes(forest_old, corner->intersection_cell_treeid[int_ind],
        corner->intersection_cell_indices[int_ind], coordinates );
  // find which corners/elements of old forest lie within the new element
  old_inside = false;
  // if start not from 0 but from variable => can call this again later as function with updated
  // start index variable
  icorner_old = 0
  for (int icorn_old=icorner_old;icorn_old<coordinates.size();icorn_old++)
  {
     icornnext = order_corner_old[icorn_old];
     corner_is_inside_element =
       t8_forest_element_point_inside (forest_new, treeid_new, elem_new, coordinates[icornnext], tolerance);
     if (old_inside .and. !corner_is_inside_elemnt){
       old_inside=false;
       // previous one was inside, next one isn't -> look for intersection old and new mesh element
       finding_intersec_fromold=1;
       direction_new=1;
       icornnow=order_corner_old[icorn_old-1];
       while(loop_new){
         // find intersection between two neighboring points in old and new grid
         intersect_point = t8_vec_segxseg(coordinates[icornnow], coordinates[icornnext],
            corner->coordinates[order_corner_new[icorn_new]], corner->coordinates[order_corner_new[icorn_new+1]],
            tolerance);
         // check if two edges intersect in exactly one point
         if (sizeof(intersect_point)/sizeof(double)==1){
           // add intersection to list of points
           point_cloud.push(intersect_point);
           // check in which direction next new corner is in old element
           // is it even possible that the direction will be changed????
           // qKH: is treeid and cell_index from corner used correctlz here?
           corner_new1_is_inside_element =
            t8_forest_element_point_inside (forest_old, corner->intersection_cell_treeid[int_ind], 
                corner->intersection_cell_indices[int_ind],
                corner->coordinates[order_corner_new[icorn_new]], tolerance);
           corner_new2_is_inside_element =
            t8_forest_element_point_inside (forest_old, corner->intersection_cell_treeid[int_ind],
                corner->intersection_cell_indices[int_ind],
                corner->coordinates[order_corner_new[icorn_new+1]], tolerance);
           // in one direction should be inside of old element -> add that corner to point cloud
           // this "inside point" determines the order in which to search for next "inside points"
           if (corner_new1_is_inside_element){
             point_cloud.push(corner->coordinates[order_corner_new[icorn_new]]);
             direction_new = -1;
             // qKH: can this even occur? is the the order not always staying constant between elements?
           }else if(corner_new2_is_inside_element){
             point_cloud.push(corner->coordinates[order_corner_new[icorn_new+1]]);
             icorn_new++;
             direction_new = 1;
           }else if (corner_new1_is_inside_element .and. corner_new2_is_inside_element){
             // something went wrong earlier!!
           }
           corner_new_is_inside_element = 1;
           // finding more corners of new forest inside of old element (if present)
           while (corner_new_is_inside_element){
             icorn_new++;
             if (direction_new==-1){
               icorn_new-=2;
             }
             // is there another new corner in old element?
             corner_new_is_inside_element = 
               t8_forest_element_point_inside (forest_old, corner->intersection_cell_treeid[int_ind],
                  corner->intersection_cell_indices[int_ind],
                  corner->coordinates[order_corner_new[icorn_new]], tolerance);
             if(corner_new_is_inside_element){
               point_cloud.push(corner->coordinates[order_corner_new[icorn_new]])
             }else{
               corner_new_is_inside_element = 0;
               // resetting indices if no more "inside points" found
               icorn_new-=1;
               if (direction_new==-1){
                 icorn_new+=2;
               }
             }
           } // while (corner_new_is_inside_element)
           finding_intersec_fromold=0;
         } // if (sizeof(intersect_point)/sizeof(double)==1) 
         if (finding_intersec_fromold){
           icorn_new++;
           // todo check that icorn_new doesn't get too large
         }
         if (!finding_intersec_fromold){
           icornnow  = order_corner_old[icorn_old];
           icorn_old++;
           icornnext = order_corner_old[icorn_old];
           // todo: add check that icorn_old doesn't get too large
         }
       } // while(loop_new)
     } // if (old_inside .and. !corner_is_inside_elemnt){
     // found corner of old forest in new element
     if (corner_is_inside_element){
       point_cloud.push(coordinates[icornnext]);
       old_inside=true;
     }
  }
  // Add check to see if just all old element corners within new element? 

//ltree_id in t8_cell_corners_t -> uebergeben
// qKH: but doesn't "corners" already store the ltree_id value after the callback (for the old forest)?

//Eckpunkte des Schnittes zurueckgeben

}

void
t8_next_element(index_in, index_out, direction)
{
  // direction +1 or -1
  // assuming element shape for now
  int order_corners_new[] = {0,1,3,2};
  size_t nr = sizeof(order_corners_new)/sizeof(int);
  int index;
  if (direction==1){
    for (index=0, index<nr, index++){
      if (index==index_in){
        // if last element in round
        if ((index+1)==nr){
          index_out = order_corners_new[0];
        }else{
          index_out = order_corners_new[index_in+1];
        }
      }
    }
  }else if(direction=-1){
    for (index=nr, index>0, index-=1){
      if (index==index_in){
        // if first element in round
        if ((index-1)==0){
          index_out = order_corners_new[nr];
        }else{
          index_out = order_corners_new[index_in-1];
        }
      }
    }
  }
}

static std::vector<t8_locidx_t>
t8_vec_segxseg(const double vec_a[3]; const double vec_b[3]; const double vec_c[3]; const double vec_d[3]; double tol)
{
  /*
   From Zoltan Csati for Matlab,  21/11/2018
   Determine intersection of segment 1 (from vec_a to vec_b) and segment 2 (from vec_c to vec_d) with a tolerance tol.
   Return the result P=
   [ ]: the segments don't intersect
   [vec_e]: the segments intersect in one point
   [vec_e, vec_f]: if the segments overlap

   From https://github.com/w8r/orourke-compc/blob/master/segseg/segseg.c 
   ---------------------------------------------------------------------
   SegSegInt: Finds the point of intersection p between two closed
   segments ab and cd.  Returns p and a char with the following meaning:
   'e': The segments collinearly overlap, sharing a point.
   'v': An endpoint (vertex) of one segment is on the other segment,
        but 'e' doesn't hold.
   '1': The segments intersect properly (i.e., they share a point and
        neither 'v' nor 'e' holds).
   '0': The segments do not intersect (i.e., they share no points).
   Note that two collinear segments that share just one point, an endpoint
   of each, returns 'e' rather than 'v' as one might expect.
   ---------------------------------------------------------------------

   Main algorithm is from
   O'Rourke, J. Computational Geometry in C, 2nd Edition, 1998.

   https://github.com/w8r/orourke-compc/blob/master/segseg/segseg.c
   and  Zoltan Csati for Matlab,  21/11/2018
   as basis
  */
   double  s, t;       /* The two parameters of the parametric eqns. */
   double num, denom;  /* Numerator and denoninator of equations. */
   std::vector<t8_locidx_t> *P;
   //char code = '?';    /* Return char characterizing intersection. */
   int X, Y;

   X=1, Y=1;

   denom = vec_a[X] * (double)( vec_d[Y] - vec_c[Y] ) +
           vec_b[X] * (double)( vec_c[Y] - vec_d[Y] ) +
           vec_d[X] * (double)( vec_b[Y] - vec_a[Y] ) +
           vec_c[X] * (double)( vec_a[Y] - vec_b[Y] );
   /* If denom is zero, then segments are parallel: handle separately. */
   if (denom < tol)  return  ParallelInt(vec_a, vec_b, vec_c, vec_d, P);
   
   num =    vec_a[X] * (double)( vec_d[Y] - vec_c[Y] ) +
            vec_c[X] * (double)( vec_a[Y] - vec_d[Y] ) +
            vec_d[X] * (double)( vec_c[Y] - vec_a[Y] );
   // update num == denom also with tolerance?
   //if ( (num < tol) || (num == denom) ) code = 'v';
   s = num / denom;
   printf("num=%lf, denom=%lf, s=%lf\n", num, denom, s);

   num = -( vec_a[X] * (double)( vec_c[Y] - vec_b[Y] ) +
            vec_b[X] * (double)( vec_a[Y] - vec_c[Y] ) +
            vec_c[X] * (double)( vec_b[Y] - vec_a[Y] ) );
   //if ( (num < tol) || (num == denom) ) code = 'v';
   t = num / denom;
   printf("num=%lf, denom=%lf, t=%lf\n", num, denom, t);

   if      ( (-tol < s) && (s < (1.0+tol)) &&
             (-tol < t) && (t < (1.0+tol)) ){
   //  code = '1';
   //else if ( (0.0 > s) || (s > 1.0) ||
   //          (0.0 > t) || (t > 1.0) )
   //  code = '0';

     P[X] = vec_a[X] + s * ( vec_b[X] - vec_a[X] );
     P[Y] = vec_a[Y] + s * ( vec_b[Y] - vec_a[Y] );
   }

   return P;
}

/*
 *  if ( (num == 0.0) || (num == denom) ) code = 'v';
   s = num / denom;
   printf("num=%lf, denom=%lf, s=%lf\n", num, denom, s);We are looking at the intersection between two convex cells.
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

static sc_array *
t8_build_corners( t8_forest_t forest )
{
  int ielem;
  t8_locidx_t itree;
  sc_array *corners = sc_array_new_count (sizeof (t8_cell_corners_t), t8_forest_get_global_num_elements( forest ));
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
    //cKH: Returns a pointer to an array element indexed by a plain int. (note on sc_array_index_int)
      t8_cell_corners_t *corner
        = (t8_cell_corners_t *) sc_array_index_int (corners, ielem);
      auto element = t8_forest_get_element ( forest, ielem_tree, &itree);

      //Calculate the corner elements of an element
      t8_forest_get_element_nodes(forest, itree, element, corner->coordinates );
      corner->volume_cell = t8_forest_element_volume( forest, itree, element );
    }
  }
  return corners;
}

void
t8_forest_conservative_remapping_planar( t8_forest_t forest_old, t8_forest_t forest_new )
{
  int ielem;
  t8_locidx_t itree;
  sc_array *corners = t8_build_corners(forest_new);
  //0) Fill structure corners with relevant information (i.e. cells of old forest which contain corners
  //   of new forest)
  //   Search already loops through all elements of new forest
  t8_forest_search(forest_old, t8_search_callback, t8_tutorial_search_query_callback, corners);
  //1) Iterieren ueber Zellen des neuen Forests
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest_new); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest_new, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
      t8_cell_corners_t *corner
        = (t8_cell_corners_t *) sc_array_index_int (corners, ielem);
      // ( t8_forest_t forest1, t8_forest_t forest2, sc_array *corners, t8_element_t *elem2 )
      // cKH: I changed the arguments here - instead of elem1 (elem of the old forest) the
      //      corners information is passed (or one element of corners)
      //3) Schnitt der Zellen bilden
      cell_intersection(forest_old, forest_new, corner, ielem, itree)
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
