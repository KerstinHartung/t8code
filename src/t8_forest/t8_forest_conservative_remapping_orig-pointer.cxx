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
TODO
- finish cell_intersection
- add volume calculation from point cloud (maybe via qhull)
- add another search algorithm t8_search_edges_element_callback
  - looking for edges via finding elements that contain one of intersection points
    and where not already considered (i.e. don't contain corner of new element
- volume check in between)
- add another search algorithm t8_search_contained_element_callback
  - looking for new elements which completely contain old ones
- maybe as a last step - if still not complete volume - look for edges which
  just intersect but don't contain any edge points in either direction
  (because these are hardest to find)
*/

T8_EXTERN_C_BEGIN ();

void
t8_next_element(int index_in, int &index_out, int direction, t8_element_shape_t element_shape);

/* Search query, a set of points. Stores information of new cell and intersecting old cells */
typedef struct
{
  /* assumption: element shape is globally constant */
  t8_element_shape_t element_shape_new;
  t8_element_shape_t element_shape_old;
  std::vector<double*> coordinates;
  std::vector<t8_locidx_t> intersection_cell_indices;
  std::vector<t8_locidx_t> intersection_cell_treeid;
  std::vector<int> intersection_type;
  /*
    intersection_type can be
    - 1: corners contained
    - 2: edges contained
    - 3: new element completely inside of old element
  */
  double volume_cell;
  double volume_intersection_cells;

} t8_cell_corners_t;

const double
  t8_element_corner_order_2D[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS]
  = { { /* T8_ECLASS_VERTEX */
         0  } ,
      { /* T8_ECLASS_LINE */
         0, 1  },
      { /* T8_ECLASS_QUAD */
         0, 1, 3, 2  },
      { /* T8_ECLASS_TRIANGLE */
         0, 1, 2  },
      { /* T8_ECLASS_HEX */
         0, 1, 3, 2 } };

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
                                 std::vector<double*> &out_coords, t8_element_shape_t &element_shape)
{
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
  element_shape = scheme->t8_element_shape (element);
  int num_corner = t8_eclass_num_vertices[element_shape];
  //out_coords.reserve(num_corner);
  //double test[3]; - pass as argument to routine to be able to access content outside
  for( int icorner = 0; icorner < num_corner; icorner++ )
  {
    // fill vector
    out_coords.push_back(0);
    const double *ref_coords = t8_element_corner_ref_coords[element_shape][icorner];
    //double local[3] = {0,0,0};
    //out_coords.at(icorner) = local;
    //todoKH: need to release memory after usage if this approach is kept
    out_coords.at(icorner) = new double[3];
    t8_forest_element_from_ref_coords (forest, ltreeid, element, ref_coords, 1, out_coords.at(icorner), NULL );
    //t8_forest_element_from_ref_coords (forest, ltreeid, element, ref_coords, 1, test, NULL);
    //out_coords.at(icorner)[0] = test[0];
  }
  //  std::cout << "size of coordinates "<<out_coords.size()<<"\n";
  //  std::cout << *out_coords.at(0)<<"\n";
  //std::cout << out_coords.at(0)[0]<< " "<< out_coords.at(1)[1]<< " "<< out_coords.at(2)[2]<<'\n';
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
t8_search_corners_element_callback (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, const int is_leaf,
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
static int
t8_search_corners_query_callback (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                   const int is_leaf, t8_element_array_t *leaf_elements, t8_locidx_t tree_leaf_index,
                                   void *query, size_t query_index)
{
  int corner_is_inside_element;
  t8_cell_corners_t *corners_of_cell = (t8_cell_corners_t *) query;
 // T8_ASSERT (user_data != NULL);
  T8_ASSERT (query != NULL);
  /* Numerical tolerance for the is_inside_element check. */
  const double tolerance = 1e-8;

  // search for new corners in old element, i.e.  forest below is forest_old
  for( size_t icorner = 0; icorner < corners_of_cell->coordinates.size(); icorner++ )
  {
    corner_is_inside_element = t8_forest_element_point_inside (forest, ltreeid, element, corners_of_cell->coordinates.at(icorner), tolerance);
    if (corner_is_inside_element) {
      if (is_leaf) {
        t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
        //add index of the cell to list of new element
        //Kh: use . instead of -> for push_back
        corners_of_cell->intersection_cell_indices.push_back( element_index );
        corners_of_cell->intersection_cell_treeid.push_back( ltreeid );
        corners_of_cell->intersection_type.push_back( 1 );
        // qkh: stop here automatically since return 1 (i.e. leave function) next
        //      leaf elements are not searched again, right?
        //      need to avoid adding the same element several times
        // stop loop now to avoid adding the same element several times
      }
      /* The particles is inside the element (finding one corner is sufficient). This query should remain active.
       * If this element is not a leaf the search will continue with its children. */
      return 1;
    }
  }
  // qkh: okay to move this out of loop? otherwise only first corner is tested, right?
  if (!corner_is_inside_element){
    /* The particle is not inside the element. Deactivate this query.
     * If no active queries are left, the search will stop for this element and its children. */
    return 0;
  }
}

//static std::vector<t8_locidx_t>*
//std::vector<double>
//double
void
t8_vec_segxseg(const double vec_a[3], const double vec_b[3], const double vec_c[3], const double vec_d[3], double tol, double P[3], char code)
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
   code = '?';    /* Return char characterizing intersection. */
   int X, Y, Z;

   X=0, Y=1, P= new double(3);

   switch (X+Y+Z){
     case 1:
       Z = 2;
       break;
     case 2:
       Z = 1;
       break;
     case 3:
       Z = 0;
       break;
     default:
       abort;
   }

   std::cout<<"start segxseg\n";
   /*
   printf("veca=%lf, %lf, %lf\n",vec_a[0], vec_a[1], vec_a[2]);
   printf("vecb=%lf, %lf, %lf\n",vec_b[0], vec_b[1], vec_b[2]);
   printf("vecc=%lf, %lf, %lf\n",vec_c[0], vec_c[1], vec_c[2]);
   printf("vecd=%lf, %lf, %lf\n",vec_d[0], vec_d[1], vec_d[2]);
   */
   denom = vec_a[X] * (double)( vec_d[Y] - vec_c[Y] ) +
           vec_b[X] * (double)( vec_c[Y] - vec_d[Y] ) +
           vec_d[X] * (double)( vec_b[Y] - vec_a[Y] ) +
           vec_c[X] * (double)( vec_a[Y] - vec_b[Y] );
   /* If denom is zero, then segments are parallel: handle separately. */
   // denom == 0
   if (std::abs(denom) < tol) code = 'p'; //  return  ParallelInt(vec_a, vec_b, vec_c, vec_d, P);

   num =    vec_a[X] * (double)( vec_d[Y] - vec_c[Y] ) +
            vec_c[X] * (double)( vec_a[Y] - vec_d[Y] ) +
            vec_d[X] * (double)( vec_c[Y] - vec_a[Y] );
   // denom == 0 || num == denom
   if ( (std::abs(num) < tol) || (std::abs(num - denom) < tol) ) code = 'v';
   s = num / denom;
   //printf("num=%lf, denom=%lf, s=%lf\n", num, denom, s);

   num = -( vec_a[X] * (double)( vec_c[Y] - vec_b[Y] ) +
            vec_b[X] * (double)( vec_a[Y] - vec_c[Y] ) +
            vec_c[X] * (double)( vec_b[Y] - vec_a[Y] ) );
   //  denom == 0 || num == denom
   if ( (std::abs(num) < tol) || (std::abs(num - denom) < tol) ) code = 'v';
   t = num / denom;
   //printf("num=%lf, denom=%lf, t=%lf\n", num, denom, t);

   // including tolerance - try to capture points by code = '1'
   if        ( (-tol < s) && (s < (1.0+tol)) &&
             (-tol < t) && (t < (1.0+tol)) ){
     code = '1';
   }else if ( ((-tol > s) || (s > (1.0+tol)) ||
             (-tol > t) || (t > (1.0+tol))) &&
             std::isfinite(s) && std::isfinite(t)){
     code = '0';
   }

   P[X] = vec_a[X] + s * ( vec_b[X] - vec_a[X] );
   P[Y] = vec_a[Y] + s * ( vec_b[Y] - vec_a[Y] );
   P[Z] = vec_a[Z];
   // ensure that Z-axis is the same for all points?

   if (code == '1') printf("P=%lf, %lf, %lf\n", P[X], P[Y], P[Z]);
   //printf("code=%c\n", code);
   //return P;
}

// determine intersection based on indices of old forest elements overlapping currently
// selected new element
void //rueckgabewert anpassen oder als pointer uebergeben - Punktwolke
cell_intersection( t8_forest_t forest_old, t8_forest_t forest_new, t8_cell_corners_t *corner,
        const t8_element_t *elem_new, t8_locidx_t treeid_new )
{
  // information about old forest
  std::vector<double*> coordinates;
  const double tolerance = 1e-8;
  std::vector<double*> point_cloud;
  bool new_inside, corner_is_inside_element, corner_new_is_inside_element;
  bool finding_intersec_fromold, search_cross;
  //t8_locidx_t int_ind;
  int int_ind, is_3D;
  int forest_old_nr_elem, icorner_new, icorn_new_cur;
  int icorn_new_prev,icorn_new, icorn_old;
  int icorn_old_cur, icorn_new_next, icorn_old_next;
  double intersect_point[3];
  char intersect_state;
  double search_z;
  int z_shift;

  /* todo
   - update routine: first started with old corners - but probably makes more sense to start with new corner
     as this is anyways what we already found
   - add outer loop over all old elements with at least one corner inside
   - add case check (and thus return char for test): if points are "inside"
     but in reality only sharing a segment then can either have after an "e" intersect
     a) next one outside or b) next one inside, in case a) abort since not a true
    "inside" element
   - add check for parallel orientation
  */

  // cKH
  /*
   * - start with first old element (eventually have loop over old elements here)
   * - start at first corner of forest_new
   * - iterate over corners and check if inside of old element
   *   - if yes: add to point_cloud and continue with next forest_new element's corner
   *   - if no (and already found one point inside): look for cross section point between new and old element
   *      - iterate over corners of old element untill cross section found
   *      - add cross section point to point_cloud
   *      - (check in which direction next point of old forest element inside of new  element) - decided always counterclockwise
   *      - check for further corners of old element in new element and add to point cloud
   *      - once no longer inside of new element: revert cross section point check (i.e.
   *        cycle through new element corners until intersection with current old element corners found)
   *   - if no: continue
   * - stop when arrived at first corner again
   * - pass point cloud list
  */

  //number of elements of old forest within new forest element
  forest_old_nr_elem = corner->intersection_cell_indices.size();

  //todo: add loop over forest_old_nr_elem (independent outside)
  int_ind = 0;  // index to iterate through old elements intersecting with one new element
  //todo: need to set up point_cloud with index storing which element
  //      of old forest is curently handled, so that find associated volume and value
  //      also store volume of this element there
  //      when calculating volume can add all compontens together and build array
  //      of new elements and their user_data

  //Calculate the corner elements of the current old tree
  // add to function t8_build_corners_elem(corner, int_ind, forest_old)
  t8_locidx_t ltree = corner->intersection_cell_treeid.at(int_ind);
  t8_locidx_t element_index = corner->intersection_cell_indices.at(int_ind);
  t8_tree_t tree = t8_forest_get_tree(forest_old, ltree);
  t8_element_t* element = t8_forest_get_tree_element(tree, element_index);
  t8_forest_get_element_nodes(forest_old, ltree, element, coordinates, corner->element_shape_old);
  // find which corners/elements of old forest lie within the new element
  //qKH: do we need this, i.e. might this happen several times? otherwise can just
  //     check length of point_inside/=0
  new_inside = true;
  search_cross = true;
  // if start not from 0 but from variable => can call this again later as function with updated
  // start index variable
  icorner_new = 0;
  // need cycling through to icorner_new again?
  // yes, probably - because start search only truly from corner which is inside, that might be the
  // the first element but then old_inside still TRUE
  // check if 2D or 3D
  is_3D = 0;
  // re-initialize search_z
  search_z=-1;
  z_shift=0;
  //if (t8_get_eclass_scheme[corner->element_shape_old]){
  //
 // }

  // use z_shift also in main part!!
  for (icorn_new=icorner_new;icorn_new< (int) t8_eclass_num_vertices[corner->element_shape_new];icorn_new++){
    icorn_new_cur = t8_element_corner_order_2D[corner->element_shape_new][icorn_new];
    std::cout<<"in loop "<<icorn_new<<" of "<<t8_eclass_num_vertices[corner->element_shape_new]<<"\n";
    // store index alongside search? then don't need to execute this exact same call twice?
    corner_is_inside_element =
      t8_forest_element_point_inside (forest_old, ltree, element, corner->coordinates.at(icorn_new_cur), tolerance);
    if (new_inside && !corner_is_inside_element){
      new_inside = false;
      std::cout<<"switch search order\n";
      // previous one was inside old element, next one isn't -> look for intersection old and new mesh element
      // first: check new->old
      finding_intersec_fromold = false;
      t8_next_element(icorn_new, icorn_new_prev, -1, corner->element_shape_new);
      // start searching for intersection at same vertical level
      icorn_old = z_shift;
      /* look for edges of old element, starting from edge 0/face f2 */
      icorn_old_cur = t8_element_corner_order_2D[corner->element_shape_old][icorn_old];
      t8_next_element(icorn_old, icorn_old_next, +1, corner->element_shape_old);
      while(search_cross){
        // find intersection between two neighboring points in old and new grid
        t8_vec_segxseg(corner->coordinates.at(icorn_new_prev), corner->coordinates.at(icorn_new_cur),
          coordinates.at(icorn_old_cur), coordinates.at(icorn_old_next),
          tolerance, intersect_point, intersect_state);
        // check if two edges intersect in exactly one point
        if (intersect_state=='1'){
          // add intersection to list of points
          point_cloud.push_back(intersect_point);
          /* I think the "next" vertex should always be inside the old element - otherwise cycle through*/
          // qKH: is treeid and cell_index from corner used correctly here?
          if (finding_intersec_fromold){
            /* next few lines in seperate function based on int_ind?*/
            corner_new_is_inside_element =
              t8_forest_element_point_inside (forest_old, ltree, element, corner->coordinates.at(icorn_new_next), tolerance);
          }else{ //finding from new to old
            corner_new_is_inside_element =
              t8_forest_element_point_inside (forest_new, treeid_new, elem_new, coordinates[icorn_old_cur], tolerance);
          }
          // have an intersection and new element is inside old
          if(corner_new_is_inside_element){
            // add element to point cloud, advancing index in next step
            if (finding_intersec_fromold){
              point_cloud.push_back(corner->coordinates.at(icorn_new_next));
            }else{
              point_cloud.push_back(coordinates[icorn_old_cur]);
              search_cross = false;
              // leave while loop!!
            }
          } // if(corner_new_is_inside_element)
          std::cout<<"after corner_new_inside check\n";
          //}else{
            // otherwise this means that old element is much smaller and lies between two vertices
            /* find two crossing points, no corners of new part of cross section
             for this will switch search order, i.e. look where new crosses old
             staying inside of while(search_cross)
            */
          //  finding_intersec_fromold=!finding_intersec_fromold;
          //}
          // finding more corners of old forest element inside of new element (if present)
          // return to main intersection loop again afterwards
          // don't do for "from old" - then go to outside of search_cross again
          if (finding_intersec_fromold){
            while (corner_new_is_inside_element){
              icorn_old+=2; // added icorn_new_next before
              icorn_old_cur = t8_element_corner_order_2D[corner->element_shape_new][icorn_old];
              corner_new_is_inside_element =
                t8_forest_element_point_inside (forest_new,treeid_new, elem_new, coordinates[icorn_old_cur], tolerance);
              if(corner_new_is_inside_element){
                point_cloud.push_back(coordinates[icorn_old_cur]);
              }else{
                // resetting indices to last point inside if no more "inside points" found
                icorn_old-=1;
              }
            } // while (corner_new_is_inside_element)
          } // if (finding_intersec_fromold)
          // now are inside old, so check old->new
          if ((finding_intersec_fromold == 0) and (!search_cross)){
            // abort - something went wrong
          }
          /* get here directly if old much smaller than new and no
             new corners inside old element
             Assumption: if there is an intersection next element should be inside
                         and never previous one
           */
          finding_intersec_fromold = 0;
        } // if (sizeof(intersect_point)/sizeof(double)==1)
        if (finding_intersec_fromold){
          icorn_old++;
          // todo check that icorn_new doesn't get too large
          // check: corrct max index?
          if ((std::vector<double*>::size_type) icorn_old == (t8_eclass_num_vertices[corner->element_shape_old]-1)){
            // do we need one last round?
            search_cross = false;
          }
          icorn_old_cur = t8_element_corner_order_2D[corner->element_shape_new][icorn_old];
          t8_next_element(icorn_old, icorn_old_next, +1, corner->element_shape_old);
        }
        if (!finding_intersec_fromold){
          icorn_new_prev  = icorn_new_cur;
          icorn_new++;
          icorn_new_cur = t8_element_corner_order_2D[corner->element_shape_new][icorn_new];
          if ((std::vector<double*>::size_type) icorn_new == (t8_eclass_num_vertices[corner->element_shape_new]-1)){
            // I think we shouldn't get here
            search_cross = false;
          }
        }
      } // while(search_cross)
    } // if (new_inside && !corner_is_inside_element){
    // found corner of new forest in old element
    if (corner_is_inside_element){
      std::cout<<"corner is inside element\n";
      // if not first corner check that on same level
      if ((is_3D) && (search_z!=-1)){
        search_z = corner->coordinates.at(icorn_new_cur)[2];
        if (icorn_new_cur>3) z_shift = 4;
      }
      if (search_z == *corner->coordinates.at(icorn_new_cur)){
        point_cloud.push_back(corner->coordinates.at(icorn_new_cur));
        new_inside=1;
      }
    }
  }
  // Add check to see if just all old element corners within new element?

//ltree_id in t8_cell_corners_t -> uebergeben
// qKH: but doesn't "corners" already store the ltree_id value after the callback (for the old forest)?

//Eckpunkte des Schnittes zurueckgeben
}



void
t8_next_element(int index_in, int &index_out, int direction, t8_element_shape_t element_shape)
{
  // direction +1 or -1
  size_t nr = t8_eclass_num_vertices[element_shape];
  int index;

  //std::cout<<"calc next element for index_in:"<<index_in << " and direction "<< direction<<" and number of vertices "<<nr <<" \n";
  index_out = 0;
  if (element_shape==T8_ECLASS_HEX){
    if (index_in>3){
      index_in = index_in - 4;
      index_out = index_out + 4;
      nr = 4;
    }else{
      nr = 4;
    }
  }
  if (direction==1){
    for (index=0; index<(int) nr; index++){
      if (index_in==t8_element_corner_order_2D[element_shape][index]){
        // if last element in round
        if ((index+1)==nr){
          index_out = index_out + t8_element_corner_order_2D[element_shape][0];
          break;
        }else{
          index_out = index_out + t8_element_corner_order_2D[element_shape][index+1];
          break;
        }
      }
    }
  }else if(direction==-1){
    for (index=nr; index>0; index-=1){
      if (index_in==t8_element_corner_order_2D[element_shape][index]){
        // if first element in round
        if ((index)==0){
          index_out = index_out + t8_element_corner_order_2D[element_shape][nr-1];
          break;
        }else{
          index_out = index_out + t8_element_corner_order_2D[element_shape][index-1];
          break;
        }
      }
    }
  }
  //std::cout<<"index_in and out: "<<index_in<<" "<<index_out<<" \n";
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
  double test[3];
  t8_locidx_t itree;
  sc_array *corners = sc_array_new_count (sizeof (t8_cell_corners_t), t8_forest_get_global_num_elements( forest ));
  std::cout << t8_forest_get_num_local_trees (forest) <<"  "<<t8_forest_get_global_num_elements( forest )<<std::endl;
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
    //cKH: Returns a pointer to an array element indexed by a plain int. (note on sc_array_index_int)
      t8_cell_corners_t *corner
        = (t8_cell_corners_t *) sc_array_index_int (corners, ielem);
      t8_locidx_t looptree=itree;
      auto element = t8_forest_get_element ( forest, ielem_tree, &looptree);
      //Calculate the corner elements of an element
      t8_forest_get_element_nodes(forest, looptree, element, corner->coordinates, corner->element_shape_new);
      //t8_forest_get_element_nodes(forest, looptree, element, corner->coordinates, corner->element_shape_new,test );
      //if (ielem<4){
      //  printf("output get_element_nodes \n ");
        //std::cout << test[1]<<"\n";
        //std::cout << corner->coordinates.size()<<"\n";
      //  std::cout << ielem<< "  " << corner->coordinates.at(2)[1] <<"\n";
        //std::cout << ielem<< "  " <<"\n";
      //}
      corner->volume_cell = t8_forest_element_volume( forest, looptree, element );
    }
  }
  //t8_cell_corners_t *corner
  //      = (t8_cell_corners_t *) sc_array_index_int (corners, 0);
  //std::cout << corner->coordinates.at(0)[1]<<"\n";
  return corners;
}

void t8_free_corners( t8_forest_t forest, sc_array *corners)
{
  int ielem;
  t8_locidx_t itree;
  std::cout << t8_forest_get_num_local_trees (forest) <<"  "<<t8_forest_get_global_num_elements( forest )<<std::endl;
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
    //cKH: Returns a pointer to an array element indexed by a plain int. (note on sc_array_index_int)
      t8_cell_corners_t *corner
        = (t8_cell_corners_t *) sc_array_index_int (corners, ielem);
      t8_locidx_t looptree=itree;
      auto element = t8_forest_get_element ( forest, ielem_tree, &looptree);
   // t8_forest_get_element_nodes(forest, looptree, element, corner->coordinates, corner->element_shape_new);
   //t8_forest_get_element_nodes (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
   //                              std::vector<double*> &out_coords, t8_element_shape_t &element_shape)
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
  std::cout << corner->coordinates.at(2)[1]<<"\n";
  //1) find elements of old forest which contain corners of new elements (via corners)
  //   the search already loops through all elements of new forest
  t8_forest_search(forest_old, t8_search_corners_element_callback, t8_search_corners_query_callback, corners);
  //2) iterate throught all new cells
  for (itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees (forest_new); itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest_new, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
      t8_cell_corners_t* corner
        = (t8_cell_corners_t *) sc_array_index_int (corners, ielem);
      // ( t8_forest_t forest1, t8_forest_t forest2, sc_array *corners, t8_element_t *elem2 )
      // cKH: I changed the arguments here - instead of elem1 (elem of the old forest) the
      //      corners information is passed (or one element of corners)
      //3) determine intersection (based on indices of old elements containing new corners))
      // qKH: what is difference ielem and ielem_t
      ielem_t = (t8_element_t *)  sc_array_index (corners, ielem);
      if (ielem<4){
        std::cout << "in cons_remap"<<"\n";
        std::cout << corner->coordinates.at(1)[1]<<"\n";
      }
      t8_locidx_t looptree=itree;
      auto element = t8_forest_get_element ( forest_new, ielem_tree, &looptree);
      // qKH: one cell_intersection function for all types, or varying?
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
