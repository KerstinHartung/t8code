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

T8_EXTERN_C_BEGIN ();

void                                                                                                                            
t8_next_element(int index_in, int index_out, int direction, t8_element_shape_t element_shape); 

/* Search query, a set of points. Stores information of new cell and intersecting old cells */
typedef struct
{
  /* assumption: element shape is globally constant */
  t8_element_shape_t element_shape_new;
  t8_element_shape_t element_shape_old;
  std::vector<double*> coordinates;/* The corners of our cell. */
  std::vector<t8_locidx_t> *intersection_cell_indices;
  std::vector<t8_locidx_t> *intersection_cell_treeid;
  double volume_cell;
  double volume_intersection_cells;

} t8_cell_corners_t;

const double
  t8_element_corner_order_2D[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS]
  = { { /* T8_ECLASS_VERTEX */
        { 0 } },
      { /* T8_ECLASS_LINE */
        { 0 }, { 1 } },
      { /* T8_ECLASS_QUAD */
        { 0 }, { 1 }, { 3 }, { 2 } },
      { /* T8_ECLASS_TRIANGLE */
        { 0 }, { 1 }, { 2 } } };
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
                                 std::vector<double *> out_coords, t8_element_shape_t element_shape)
{
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
  element_shape = scheme->t8_element_shape (element);
  int num_corner = t8_eclass_num_vertices[element_shape];
  out_coords.reserve(num_corner);
  /*
   Thread 1 "t8_features_con" received signal SIGSEGV, Segmentation fault.
   0x00005555555d3105 in t8_geom_triangular_interpolation (coefficients=0x7fffffffe0d0, corner_values=0x5555557c0810, corner_value_dim=3, interpolation_dim=2, 
   evaluated_function=0x0) at ../src/t8_geometry/t8_geometry_helpers.c:77
      evaluated_function[i_dim] = temp[i_dim];
  */
  //double local[3];
  /*
     [libsc 0] Abort: Assertion 'i >= 0 && (size_t) i < array->elem_count'
     [libsc 0] Abort: ../sc/src/sc_containers.h:497
   */
  //static std::vector<double *> out_coords = [](){ std::vector<double*> out_coords; out_coords.reserve(3); return out_coords; }();
  for( size_t icorner = 0; icorner < num_corner; icorner++ )
  {
    const double *ref_coords = t8_element_corner_ref_coords[element_shape][icorner];
    t8_forest_element_from_ref_coords (forest, ltreeid, element, ref_coords, 1, out_coords[icorner], NULL );
    //t8_forest_element_from_ref_coords (forest, ltreeid, element, ref_coords, 1, local, NULL);
    //out_coords.push_back(local);
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
  //for(std::vector<double*>::size_type icorner = 0; icorner < corners_of_cell->coordinates.size(); icorner++ )
  {
    corner_is_inside_element = t8_forest_element_point_inside (forest, ltreeid, element, corners_of_cell->coordinates[icorner], tolerance);
    if (corner_is_inside_element) {
      if (is_leaf) {
        t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
        //add index of the cell to
        //Kh: use . instead of -> for push_back
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

//static std::vector<t8_locidx_t>*
//std::vector<double>
//double
void
t8_vec_segxseg(const double vec_a[3], const double vec_b[3], const double vec_c[3], const double vec_d[3], double tol, double P[3])
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
   //std::vector<double> P[2];
   //double P[2];
   //char code = '?';    /* Return char characterizing intersection. */
   int X, Y;

   X=1, Y=1;

   denom = vec_a[X] * (double)( vec_d[Y] - vec_c[Y] ) +
           vec_b[X] * (double)( vec_c[Y] - vec_d[Y] ) +
           vec_d[X] * (double)( vec_b[Y] - vec_a[Y] ) +
           vec_c[X] * (double)( vec_a[Y] - vec_b[Y] );
   /* If denom is zero, then segments are parallel: handle separately. */
   //if (denom < tol)  return  ParallelInt(vec_a, vec_b, vec_c, vec_d, P);

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

   //return P;
}

void //rueckgabewert anpassen oder als pointer uebergeben - Punktwolke
cell_intersection( t8_forest_t forest_old, t8_forest_t forest_new, t8_cell_corners_t *corner,
        const t8_element_t *elem_new, t8_locidx_t treeid_new )
{
  // information about old forest
  std::vector<double*> coordinates;
  // todo update to variable number of corners
  int order_corner_old[] = {0,1,3,2};
  int order_corner_new[] = {0,1,3,2};
  const double tolerance = 1e-8;
  std::vector<double*> point_cloud;
  bool old_inside, corner_is_inside_element, corner_new_is_inside_element;
  bool corner_new1_is_inside_element, corner_new2_is_inside_element;
  bool finding_intersec_fromold, loop_new;
  //t8_locidx_t int_ind;
  int int_ind;
  int forest_old_nr_elem, icorner_old, icorn_old_cur;
  int icorn_old_prev,icorn_n0, direction_new,icorn_old, icorn_new;
  int icorn_new_cur, icorn_new_next;
  double intersect_point[3];

  // qKH
  /*
    - why is the number of corners equal to coordinates.size()?
   */


  /* todo
   - add outer loop over all old elements with at least one corner inside
   - before I implement the next point: check with Chiara that I understand
     tolerance of "finding point inside" correctly, i.e. finding points on the
     face of new element
   - add case check (and thus return char for test): if points are "inside"
     but in reality sharing a segment then can either have after an "e" intersect
     a) next one outside or b) next one inside, in case a) abort since not a true
    "inside" element
   - use internal data on ordering of corners, currently hardcoded for 4 corners
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
  // qKH: why not .size()? here?
  forest_old_nr_elem = corner->intersection_cell_indices->size();

  //todo: add loop over forest_old_nr_elem (independent outside)
  int_ind = 0;  // index to iterate through old elements intersection with one new element


  //corner->intersection_cell_treeid[int_ind];

  //Calculate the corner elements of element of the current old tree
  // add to function t8_build_corners_elem
  t8_locidx_t ltree = corner->intersection_cell_treeid->at(int_ind);
  t8_locidx_t element_index = corner->intersection_cell_indices->at(int_ind);
  t8_tree_t tree = t8_forest_get_tree(forest_old, ltree);
  t8_element_t* element = t8_forest_get_tree_element(tree, element_index);
  t8_forest_get_element_nodes(forest_old, ltree, element, coordinates, corner->element_shape_old);
  //t8_forest_get_element_nodes(forest_old, corner->intersection_cell_treeid[int_ind],
  //       corner->intersection_cell_indices[int_ind], coordinates );
  //t8_forest_get_element_nodes(forest_old, (t8_locidx_t*) corner->intersection_cell_treeid[int_ind],
  //      (t8_element_t*) corner->intersection_cell_indices[int_ind], coordinates );
  // find which corners/elements of old forest lie within the new element
  //qKH: do we need this, i.e. might this happend several times? otherwise can just
  //     check length of point_inside/=0
  old_inside = 1;
  // if start not from 0 but from variable => can call this again later as function with updated
  // start index variable
  icorner_old = 0;
  loop_new = 1;
  for (icorn_old=icorner_old;icorn_old< (int) t8_eclass_num_vertices[corner->element_shape_old];icorn_old++){
     icorn_old_cur = t8_element_corner_order_2D[corner->element_shape_old][icorn_old];
     corner_is_inside_element =
       t8_forest_element_point_inside (forest_new, treeid_new, elem_new, coordinates[icorn_old_cur], tolerance);
     if (old_inside && !corner_is_inside_element){
       old_inside=false;
       // previous one was inside, next one isn't -> look for intersection old and new mesh element
       // first: check from old->new
       finding_intersec_fromold=1;
       // direction check is maybe not required
       direction_new=1;
       t8_next_element(icorn_old, icorn_old_prev, -1, corner->element_shape_old);
       icorn_new = 0;
       while(loop_new){
         /* look for edges of new element, starting from edge 0/face f2 */
         icorn_new_cur = t8_element_corner_order_2D[corner->element_shape_new][icorn_new];
         t8_next_element(icorn_new, icorn_new_next, direction_new, corner->element_shape_new);
         // find intersection between two neighboring points in old and new grid
         t8_vec_segxseg(coordinates[icorn_old_prev], coordinates[icorn_old_cur],
            corner->coordinates[icorn_new_cur], corner->coordinates[icorn_new_next],
            tolerance, intersect_point);
         // check if two edges intersect in exactly one point
         if (sizeof(intersect_point)/sizeof(double)==1){
           // add intersection to list of points
           point_cloud.push_back(intersect_point);
           // check in which direction next new corner is in old element
           // is it even possible that the direction will be changed????
           /* I think the "next" vertex should always be inside the old element */
           // qKH: is treeid and cell_index from corner used correctly here?
           t8_locidx_t element_index = corner->intersection_cell_indices->at(int_ind);
           t8_tree_t tree = t8_forest_get_tree(forest_old, corner->intersection_cell_treeid->at(int_ind));
           t8_element_t* element = t8_forest_get_tree_element(tree, element_index);
           corner_new1_is_inside_element =
            t8_forest_element_point_inside (forest_old, (t8_locidx_t) corner->intersection_cell_treeid->at(int_ind),
                element,
                corner->coordinates[icorn_new_cur], tolerance);
           corner_new2_is_inside_element =
            t8_forest_element_point_inside (forest_old, (t8_locidx_t) corner->intersection_cell_treeid->at(int_ind),
                element,
                corner->coordinates[icorn_new_next], tolerance);
           // in one direction should be inside of old element -> add that corner to point cloud
           // this "inside point" determines the order in which to search for next "inside points"
           if (corner_new1_is_inside_element){
             point_cloud.push_back(corner->coordinates[icorn_new_cur]);
             direction_new = -1;
             // qKH: can this even occur? is the the order not always staying constant between elements?
           }else if(corner_new2_is_inside_element){
             point_cloud.push_back(corner->coordinates[icorn_new_next]);
             icorn_new++;
             direction_new = 1;
           }else if (corner_new1_is_inside_element && corner_new2_is_inside_element){
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
               t8_forest_element_point_inside (forest_old, corner->intersection_cell_treeid->at(int_ind),
                  element,
                  corner->coordinates[icorn_new_cur], tolerance);
             if(corner_new_is_inside_element){
               //qkh:  error: no matching function for call to ‘std::vector<double [3]>::push_back(double*&)
               //point_cloud->push_back(corner->coordinates[order_corner_new[icorn_new]]);
             }else{
               corner_new_is_inside_element = 0;
               // resetting indices if no more "inside points" found
               icorn_new-=1;
               if (direction_new==-1){
                 icorn_new+=2;
               }
             }
           } // while (corner_new_is_inside_element)
           // now are inside new, so check new->old
           finding_intersec_fromold=0;
         } // if (sizeof(intersect_point)/sizeof(double)==1)
         if (finding_intersec_fromold){
           icorn_new++;
           if ((std::vector<double*>::size_type) icorn_new == corner->coordinates.size()) loop_new = 0;
           // todo check that icorn_new doesn't get too large
         }
         if (!finding_intersec_fromold){
           icorn_old_prev  = icorn_old_cur;
           icorn_old++;
           if ((std::vector<double*>::size_type) icorn_old == coordinates.size()) loop_new = 0;
           // todo: add check that icorn_old doesn't get too large
         }
       } // while(loop_new)
     } // if (old_inside && !corner_is_inside_element){
     // found corner of old forest in new element
     if (corner_is_inside_element){
       point_cloud.push_back(coordinates[icorn_old_cur]);
       old_inside=true;
     }
  }
  // Add check to see if just all old element corners within new element?

//ltree_id in t8_cell_corners_t -> uebergeben
// qKH: but doesn't "corners" already store the ltree_id value after the callback (for the old forest)?

//Eckpunkte des Schnittes zurueckgeben

}



void
t8_next_element(int index_in, int index_out, int direction, t8_element_shape_t element_shape)
{
  // direction +1 or -1
  size_t nr = t8_eclass_num_vertices[element_shape];
  int index;
  if (direction==1){
    for (index=0; index<(int) nr; index++){
      if (index==index_in){
        // if last element in round
        if ((index+1)==(int) nr){
          index_out = t8_element_corner_order_2D[element_shape][0];
        }else{
          index_out = t8_element_corner_order_2D[element_shape][index_in+1];
        }
      }
    }
  }else if(direction==-1){
    for (index=nr; index>0; index-=1){
      if (index==index_in){
        // if first element in round
        if ((index-1)==0){
          index_out = t8_element_corner_order_2D[element_shape][nr];
        }else{
          index_out = t8_element_corner_order_2D[element_shape][index_in-1];
        }
      }
    }
  }
}

//void
//t8_next_element(int index_in, int index_out, int direction)
//{
//  // direction +1 or -1
//  // assuming element shape for now
//  int order_corners_new[] = {0,1,3,2};
//  size_t nr = sizeof(order_corners_new)/sizeof(int);
//  int index;
//  if (direction==1){
//    for (index=0; index<(int) nr; index++){
//      if (index==index_in){
//        // if last element in round
//        if ((index+1)==(int) nr){
//          index_out = order_corners_new[0];
//        }else{
//          index_out = order_corners_new[index_in+1];
//        }
//      }
//    }
//  }else if(direction==-1){
//    for (index=nr; index>0; index-=1){
//      if (index==index_in){
//        // if first element in round
//        if ((index-1)==0){
//          index_out = order_corners_new[nr];
//        }else{
//          index_out = order_corners_new[index_in-1];
//        }
//      }
//    }
//  }
//}

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
      t8_forest_get_element_nodes(forest, itree, element, corner->coordinates, corner->element_shape_new );
      corner->volume_cell = t8_forest_element_volume( forest, itree, element );
    }
  }
  return corners;
}

void
t8_forest_conservative_remapping_planar( t8_forest_t forest_old, t8_forest_t forest_new )
{
  int ielem;
  t8_element_t *ielem_t;
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
      t8_cell_corners_t* corner
        = (t8_cell_corners_t *) sc_array_index_int (corners, ielem);
      // ( t8_forest_t forest1, t8_forest_t forest2, sc_array *corners, t8_element_t *elem2 )
      // cKH: I changed the arguments here - instead of elem1 (elem of the old forest) the
      //      corners information is passed (or one element of corners)
      //3) Schnitt der Zellen bilden
      //ielem_t = (t8_element_t *) ielem;
      ielem_t = (t8_element_t *)  sc_array_index (corners, ielem);
      cell_intersection(forest_old, forest_new, corner, ielem_t, itree);
      // add volume as output argument
      //point_cloud_to_volume(points, ...);
    }
  }

  // cKH
  // the next few lines can probably be integrated in the above nested loop
  // volume_new_remap = 0 (can this be stored in volume_intersection_cells?)
  //loop over trees and number of elements of new forest     //volume_new = t8_forest_element_volume (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element)
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
