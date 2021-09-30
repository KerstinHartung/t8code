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

/** \file t8_default_common_cxx.hxx
 * We provide some functions that are useful across element classes.
 */

#ifndef T8_DEFAULT_COMMON_CXX_HXX
#define T8_DEFAULT_COMMON_CXX_HXX

#include <t8_element_cxx.hxx>

/* Macro to check whether a pointer (VAR) to a base class, comes from an
 * implementation of a child class (TYPE). */
#define T8_COMMON_IS_TYPE(VAR, TYPE) \
  ((dynamic_cast<TYPE> (VAR)) != NULL)

class               t8_default_scheme_common_c:public t8_eclass_scheme_c
{
public:
  /** Destructor for all default schemes */
  virtual ~ t8_default_scheme_common_c ();

  /** Compute the number of corners of a given element. */
  virtual int         t8_element_num_corners (const t8_element_t * elem);

  /** Allocate space for a bunch of elements. */
  virtual void        t8_element_new (int length, t8_element_t ** elem);

  /** Deallocate space for a bunch of elements. */
  virtual void        t8_element_destroy (int length, t8_element_t ** elem);

  /** Return the shape of an element */
  virtual t8_element_shape_t t8_element_shape (const t8_element_t * elem);

  /** Count how many leaf descendants of a given uniform level an element would produce.
   * \param [in] t     The element to be checked.
   * \param [in] level A refinement level.
   * \return Suppose \a t is uniformly refined up to level \a level. The return value
   * is the resulting number of elements (of the given level).
   * Each default element (except pyramids) refines into 2^{dim * (level - level(t))}
   * children.
   */
  virtual t8_gloidx_t t8_element_count_leafs (const t8_element_t * t,
                                              int level);

  /** Compute the maximum number of siblings of an element or any descendants of it. That is the number of 
   * Children of its parent.
   * \param [in] elem The element.
   * \return          The maximum number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   */
  virtual int         t8_element_max_num_siblings (const t8_element_t *
                                                   elem) const;

  /** Compute the number of siblings of an element. That is the number of 
   * Children of its parent.
   * \param [in] elem The element.
   * \return          The number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   */
  virtual int         t8_element_num_siblings (const t8_element_t *
                                               elem) const;

  /** Count how many leaf descendants of a given uniform level the root element will produce.
   * \param [in] level A refinement level.
   * \return The value of \ref t8_element_count_leafs if the input element
   *      is the root (level 0) element.
   */
  virtual t8_gloidx_t t8_element_count_leafs_from_root (int level);

  /** The common implementation of the general function for the default scheme
   * has no effect. This function literally does nothing.
   * The tri, tet and prism scheme override this implementation with a function that
   * stores the type of the element in \a outdata.
   *  \param [in] elem A valid element
   *  \param [in] indata Is ignored. Can be NULL.
   *  \param [out] outdata Is ignored. Can be NULL.
   * \note Calling this function has no effect. See the specialized implementations in
   * t8_default_tri_cxx.hxx, t8_default_tet_cxx.hxx and t8_default_prism_cxx.hxx.
   */
  virtual void        t8_element_general_function (const t8_element_t * elem,
                                                   const void *indata,
                                                   void *outdata);

  /** Construct a subelement */
  virtual void        t8_element_to_subelement (const t8_element_t * elem,
                                                int type, t8_element_t * c[]);

  /* TODO: comment */
  virtual int         t8_element_test_if_subelement (const t8_element * elem);

  /** Determine the number of subelements, used to remove hanging nodes from a element of a given type */
  virtual int         t8_element_get_number_of_subelements (int
                                                            subelement_type,
                                                            const
                                                            t8_element_t *
                                                            elem);

  /* TODO: comment */
  virtual int         t8_element_get_subelement_type (const
                                                      t8_element * elem);

  /** TODO: comment */
  virtual int         t8_element_get_subelement_id (const
                                                    t8_element * elem);

  /* TODO: comment */
  virtual void        t8_element_get_element_data (const t8_element * elem,
                                                   int anchor_node[],
                                                   int level[],
                                                   int subelement_data[]);

  /** TODO: comment */
  virtual int         t8_element_find_neighbor_in_transition_cell (const t8_element_t * elem, 
                                                                   const t8_element_t *neigh, 
                                                                   int elem_face);

  /** TODO: comment */
  virtual int         t8_element_adjust_subelement_neighbor_index (const
                                                                   t8_element_t
                                                                   * elem,
                                                                   const
                                                                   t8_element_t
                                                                   * neigh,
                                                                   int
                                                                   elem_index,
                                                                   int
                                                                   elem_face);

  /** This function will determine the location of a specific subelement in the parent element.
   *  Since different subelement types are possible, it is a priori not known where for example the
   *  subelement with id 3 is located. 
   *  \param [in] elem A valid subelement
   *  \param [out] An array, whose entries are face_number, split and sub_face_id
   *                      face_number: the face number, the given subelement is adjacent to (value between 0 and 3)
   *                      split: whether there is a hanging node at the face, the subelement is adjacent to 
   *                             (value 0 if there is not hanging node and 1 if there is one)
   *                      sub_face_id: if there is a hanging node at the face, it is important to know if the given 
   *                                   subelement is the first or the second subelement at this face
   *                                   (value 0 if it is the first and 1 if it is the second)
   *  The information in the location can be used to automatically determine the verticies of any subelement.
   *  Since this function is only used to determine the vertices of subelements, it can be declared as a private/protected function.
   */
  virtual void        t8_element_get_location_of_subelement (const
                                                             t8_element_t *
                                                             elem,
                                                             int location[]);
};

#endif /* !T8_DEFAULT_COMMON_CXX_HXX */
