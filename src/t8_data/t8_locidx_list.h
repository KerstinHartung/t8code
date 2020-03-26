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

/** \file t8_locidx_list.h
 * We define a linked list of t8_locidx_t as a wrapper around sc_list_t.
 */

#ifndef T8_LOCIDX_LIST_H
#define T8_LOCIDX_LIST_H

#include <t8.h>
#include <sc_shmem.h>

typedef struct t8_shmem_array *t8_shmem_array_t;

T8_EXTERN_C_BEGIN ();

/** Defines a linked list of t8_locidx_t entries.
 * The entries are allocated and stored in a mempool.
 */
typedef struct
{
  sc_list_t           list;
                  /**< The linked list */
  sc_mempool_t        indices;
                        /**< The mempool that stores the indices */
} t8_locidx_list_t;

/** Allocate a new t8_locidx_list object and return a poitner to it.
 * \return A pointer to a newly allocated t8_locidx_list_t
 */
t8_locidx_list_t   *t8_locidx_list_new ();

/** Initialize an allocated t8_locidx_list
 * \param [in, out] Pointer to an allocated t8_locidx_list_t
 * \note If this list is already filled the all its elements are free'd
 */
void                t8_locidx_list_init (t8_locidx_list_t * list);

/** Query whether a t8_locidx_list_t is initialized.
 * \param [in] list A pointer to a list.
 * \return True (non-zero) if and only if the list is properly initialized.
 * \note Calling this with \a list = NULL is valid and will return 0.
 */
int                 t8_locidx_list_is_initialized (const t8_locidx_list_t *
                                                   list);

/** Remove the last element in a list and return its entry.
 * \param [in,out] list An initliazed list with at least one entry.
 *                      On return its last entry will have been removed.
 * \return  The last entry of \a list.
 * \note It is illegal to call this function if \a list does not have any elements.
 */
t8_locidx_t         t8_locidx_list_pop (t8_locidx_list_t * list);

/** Returns the number of entries in a list.
 * \param [in] list An initialized list.
 * \return          The number of entries in \a list.
 */
size_t              t8_locidx_list_count (const t8_locidx_list_t * list);

/** Append a new entry to the end of a list.
 * \param [in,out] list An initialized list.
 * \param [in]     entry The new entry.
 */
void                t8_locidx_list_append (t8_locidx_list_t * list,
                                           const t8_locidx_t entry);

/** Free all elements of a list.
 * \param [in, out] list An initialized list.
 * After calling this function \a list will still be initialized but
 * all its elements will have been free'd and its count will be zero.
 */
void                t8_locidx_list_reset (t8_locidx_list_t * list);

T8_EXTERN_C_END ();

#endif /* !T8_LOCIDX_LIST_H */
