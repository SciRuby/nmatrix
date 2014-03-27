/////////////////////////////////////////////////////////////////////
// = NMatrix
//
// A linear algebra library for scientific computation in Ruby.
// NMatrix is part of SciRuby.
//
// NMatrix was originally inspired by and derived from NArray, by
// Masahiro Tanaka: http://narray.rubyforge.org
//
// == Copyright Information
//
// SciRuby is Copyright (c) 2010 - 2014, Ruby Science Foundation
// NMatrix is Copyright (c) 2012 - 2014, John Woods and the Ruby Science Foundation
//
// Please see LICENSE.txt for additional copyright notices.
//
// == Contributing
//
// By contributing source code to SciRuby, you agree to be bound by
// our Contributor Agreement:
//
// * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
//
// == list.c
//
// List-of-lists n-dimensional matrix storage. Uses singly-linked
// lists.
/*
 * Standard Includes
 */

#include <ruby.h>
#include <algorithm> // std::min
#include <iostream>
#include <vector>
#include <list>

/*
 * Project Includes
 */

#include "../../types.h"

#include "../../data/data.h"

#include "../dense/dense.h"
#include "../common.h"
#include "list.h"

#include "../../math/math.h"
#include "../../util/sl_list.h"

/*
 * Macros
 */

/*
 * Global Variables
 */


extern "C" {
static void slice_set_single(LIST_STORAGE* dest, LIST* l, void* val, size_t* coords, size_t* lengths, size_t n);
static void __nm_list_storage_unregister_temp_value_list(std::list<VALUE*>& temp_vals);
static void __nm_list_storage_unregister_temp_list_list(std::list<LIST*>& temp_vals, size_t recursions);
}

namespace nm { namespace list_storage {

/*
 * Forward Declarations
 */

class RecurseData {
public:
  // Note that providing init_obj argument does not override init.
  RecurseData(const LIST_STORAGE* s, VALUE init_obj__ = Qnil) : ref(s), actual(s), shape_(s->shape), offsets(s->dim, 0), init_(s->default_val), init_obj_(init_obj__) {
    while (actual->src != actual) {
      for (size_t i = 0; i < s->dim; ++i) // update offsets as we recurse
        offsets[i] += actual->offset[i];
      actual = reinterpret_cast<LIST_STORAGE*>(actual->src);
    }
    nm_list_storage_register(actual);
    nm_list_storage_register(ref);
    actual_shape_ = actual->shape;

    if (init_obj_ == Qnil) {
      init_obj_ = s->dtype == nm::RUBYOBJ ? *reinterpret_cast<VALUE*>(s->default_val) : rubyobj_from_cval(s->default_val, s->dtype).rval;
    }
    nm_register_value(init_obj_);
  }

  ~RecurseData() {
    nm_unregister_value(init_obj_);
    nm_list_storage_unregister(ref);
    nm_list_storage_unregister(actual);
  }

  dtype_t dtype() const { return ref->dtype; }


  size_t dim() const { return ref->dim; }

  size_t ref_shape(size_t rec) const {
    return shape_[ref->dim - rec - 1];
  }

  size_t* copy_alloc_shape() const {
    size_t* new_shape = NM_ALLOC_N(size_t, ref->dim);
    memcpy(new_shape, shape_, sizeof(size_t)*ref->dim);
    return new_shape;
  }

  size_t actual_shape(size_t rec) const {
    return actual_shape_[actual->dim - rec - 1];
  }

  size_t offset(size_t rec) const {
    return offsets[ref->dim - rec - 1];
  }

  void* init() const {
    return init_;
  }

  VALUE init_obj() const { return init_obj_; }

  LIST* top_level_list() const {
    return reinterpret_cast<LIST*>(actual->rows);
  }

  const LIST_STORAGE* ref;
  const LIST_STORAGE* actual;

  size_t* shape_; // of ref
  size_t* actual_shape_;
protected:
  std::vector<size_t> offsets; // relative to actual
  void* init_;
  VALUE init_obj_;

};


template <typename LDType, typename RDType>
static LIST_STORAGE* cast_copy(const LIST_STORAGE* rhs, nm::dtype_t new_dtype);

template <typename LDType, typename RDType>
static bool eqeq_r(RecurseData& left, RecurseData& right, const LIST* l, const LIST* r, size_t rec);

template <typename SDType, typename TDType>
static bool eqeq_empty_r(RecurseData& s, const LIST* l, size_t rec, const TDType* t_init);

/*
 * Recursive helper for map_merged_stored_r which handles the case where one list is empty and the other is not.
 */
static void map_empty_stored_r(RecurseData& result, RecurseData& s, LIST* x, const LIST* l, size_t rec, bool rev, const VALUE& t_init) {
  if (s.dtype() == nm::RUBYOBJ) {
    nm_list_storage_register_list(l, rec);
  }
  if (result.dtype() == nm::RUBYOBJ) {
    nm_list_storage_register_list(x, rec);
  }

  NODE *curr  = l->first,
       *xcurr = NULL;

  // For reference matrices, make sure we start in the correct place.
  size_t offset   = s.offset(rec);
  size_t x_shape  = s.ref_shape(rec);

  while (curr && curr->key < offset) {  curr = curr->next;  }
  if (curr && curr->key - offset >= x_shape) curr = NULL;

  if (rec) {
    std::list<LIST*> temp_vals;
    while (curr) {
      LIST* val = nm::list::create();
      map_empty_stored_r(result, s, val, reinterpret_cast<const LIST*>(curr->val), rec-1, rev, t_init);

      if (!val->first) nm::list::del(val, 0);
      else {
        nm_list_storage_register_list(val, rec-1);
	temp_vals.push_front(val);
        nm::list::insert_helper(x, xcurr, curr->key - offset, val);
      } 
      curr = curr->next;
      if (curr && curr->key - offset >= x_shape) curr = NULL;
    }
    __nm_list_storage_unregister_temp_list_list(temp_vals, rec-1);
  } else {
    std::list<VALUE*> temp_vals;
    while (curr) {
      VALUE val, s_val = rubyobj_from_cval(curr->val, s.dtype()).rval;
      if (rev) val = rb_yield_values(2, t_init, s_val);
      else     val = rb_yield_values(2, s_val, t_init);

      nm_register_value(val);

      if (rb_funcall(val, rb_intern("!="), 1, result.init_obj()) == Qtrue) {
        xcurr = nm::list::insert_helper(x, xcurr, curr->key - offset, val);
        temp_vals.push_front(reinterpret_cast<VALUE*>(xcurr->val));
        nm_register_value(*reinterpret_cast<VALUE*>(xcurr->val));
      }
      nm_unregister_value(val);

      curr = curr->next;
      if (curr && curr->key - offset >= x_shape) curr = NULL;
    }
    __nm_list_storage_unregister_temp_value_list(temp_vals);
  }

  if (s.dtype() == nm::RUBYOBJ){
    nm_list_storage_unregister_list(l, rec);
  }
  if (result.dtype() == nm::RUBYOBJ) {
    nm_list_storage_unregister_list(x, rec);
  }

}


/*
 * Recursive helper function for nm_list_map_stored
 */
static void map_stored_r(RecurseData& result, RecurseData& left, LIST* x, const LIST* l, size_t rec) {
  if (left.dtype() == nm::RUBYOBJ) {
    nm_list_storage_register_list(l, rec);
  }
  if (result.dtype() == nm::RUBYOBJ) {
    nm_list_storage_register_list(x, rec);
  }
  NODE *lcurr = l->first,
       *xcurr = x->first;

  // For reference matrices, make sure we start in the correct place.
  while (lcurr && lcurr->key < left.offset(rec))  {  lcurr = lcurr->next;  }

  if (lcurr && lcurr->key - left.offset(rec) >= result.ref_shape(rec))  lcurr = NULL;

  if (rec) {
    std::list<LIST*> temp_vals;
    while (lcurr) {
      size_t key;
      LIST*  val = nm::list::create();
      map_stored_r(result, left, val, reinterpret_cast<const LIST*>(lcurr->val), rec-1);
      key        = lcurr->key - left.offset(rec);
      lcurr      = lcurr->next;

      if (!val->first) nm::list::del(val, 0); // empty list -- don't insert
      else {
        nm_list_storage_register_list(val, rec-1);
        temp_vals.push_front(val);
        xcurr = nm::list::insert_helper(x, xcurr, key, val);
      }
      if (lcurr && lcurr->key - left.offset(rec) >= result.ref_shape(rec)) lcurr = NULL;
    }
    __nm_list_storage_unregister_temp_list_list(temp_vals, rec-1);
  } else {
    std::list<VALUE*> temp_vals;
    while (lcurr) {
      size_t key;
      VALUE  val;

      val   = rb_yield_values(1, rubyobj_from_cval(lcurr->val, left.dtype()).rval);
      key   = lcurr->key - left.offset(rec);
      lcurr = lcurr->next;

      if (!rb_equal(val, result.init_obj())) {
        xcurr = nm::list::insert_helper(x, xcurr, key, val);
        temp_vals.push_front(reinterpret_cast<VALUE*>(xcurr->val));
        nm_register_value(*reinterpret_cast<VALUE*>(xcurr->val));
      }

      if (lcurr && lcurr->key - left.offset(rec) >= result.ref_shape(rec)) lcurr = NULL;
    }
    __nm_list_storage_unregister_temp_value_list(temp_vals);
  }

  if (left.dtype() == nm::RUBYOBJ) {
    nm_list_storage_unregister_list(l, rec);
  }
  if (result.dtype() == nm::RUBYOBJ) {
    nm_list_storage_unregister_list(x, rec);
  }
}



/*
 * Recursive helper function for nm_list_map_merged_stored
 */
static void map_merged_stored_r(RecurseData& result, RecurseData& left, RecurseData& right, LIST* x, const LIST* l, const LIST* r, size_t rec) {
  if (left.dtype() == nm::RUBYOBJ) {
    nm_list_storage_register_list(l, rec);
  }
  if (right.dtype() == nm::RUBYOBJ) {
    nm_list_storage_register_list(r, rec);
  }
  if (result.dtype() == nm::RUBYOBJ) {
    nm_list_storage_register_list(x, rec);
  }


  NODE *lcurr = l->first,
       *rcurr = r->first,
       *xcurr = x->first;

  // For reference matrices, make sure we start in the correct place.
  while (lcurr && lcurr->key < left.offset(rec))  {  lcurr = lcurr->next;  }
  while (rcurr && rcurr->key < right.offset(rec)) {  rcurr = rcurr->next;  }

  if (rcurr && rcurr->key - right.offset(rec) >= result.ref_shape(rec)) rcurr = NULL;
  if (lcurr && lcurr->key - left.offset(rec) >= result.ref_shape(rec))  lcurr = NULL;

  if (rec) {
    std::list<LIST*> temp_vals;
    while (lcurr || rcurr) {
      size_t key;
      LIST*  val = nm::list::create();

      if (!rcurr || (lcurr && (lcurr->key - left.offset(rec) < rcurr->key - right.offset(rec)))) {
        map_empty_stored_r(result, left, val, reinterpret_cast<const LIST*>(lcurr->val), rec-1, false, right.init_obj());
        key   = lcurr->key - left.offset(rec);
        lcurr = lcurr->next;
      } else if (!lcurr || (rcurr && (rcurr->key - right.offset(rec) < lcurr->key - left.offset(rec)))) {
        map_empty_stored_r(result, right, val, reinterpret_cast<const LIST*>(rcurr->val), rec-1, true, left.init_obj());
        key   = rcurr->key - right.offset(rec);
        rcurr = rcurr->next;
      } else { // == and both present
        map_merged_stored_r(result, left, right, val, reinterpret_cast<const LIST*>(lcurr->val), reinterpret_cast<const LIST*>(rcurr->val), rec-1);
        key   = lcurr->key - left.offset(rec);
        lcurr = lcurr->next;
        rcurr = rcurr->next;
      }


      if (!val->first) nm::list::del(val, 0); // empty list -- don't insert
      else {
        nm_list_storage_register_list(val, rec-1);
        temp_vals.push_front(val);
        xcurr = nm::list::insert_helper(x, xcurr, key, val);
      }
      if (rcurr && rcurr->key - right.offset(rec) >= result.ref_shape(rec)) rcurr = NULL;
      if (lcurr && lcurr->key - left.offset(rec) >= result.ref_shape(rec)) lcurr = NULL;
    }
    __nm_list_storage_unregister_temp_list_list(temp_vals, rec-1);
  } else {
    std::list<VALUE*> temp_vals;
    while (lcurr || rcurr) {
      size_t key;
      VALUE  val;

      if (!rcurr || (lcurr && (lcurr->key - left.offset(rec) < rcurr->key - right.offset(rec)))) {
        val   = rb_yield_values(2, rubyobj_from_cval(lcurr->val, left.dtype()).rval, right.init_obj());
        key   = lcurr->key - left.offset(rec);
        lcurr = lcurr->next;
      } else if (!lcurr || (rcurr && (rcurr->key - right.offset(rec) < lcurr->key - left.offset(rec)))) {
	      val   = rb_yield_values(2, left.init_obj(), rubyobj_from_cval(rcurr->val, right.dtype()).rval);
        key   = rcurr->key - right.offset(rec);
        rcurr = rcurr->next;
      } else { // == and both present
        val   = rb_yield_values(2, rubyobj_from_cval(lcurr->val, left.dtype()).rval, rubyobj_from_cval(rcurr->val, right.dtype()).rval);
        key   = lcurr->key - left.offset(rec);
        lcurr = lcurr->next;
        rcurr = rcurr->next;
      }

      nm_register_value(val);

      if (rb_funcall(val, rb_intern("!="), 1, result.init_obj()) == Qtrue) {
        xcurr = nm::list::insert_helper(x, xcurr, key, val);
        temp_vals.push_front(reinterpret_cast<VALUE*>(xcurr->val));
        nm_register_value(*reinterpret_cast<VALUE*>(xcurr->val));
      }

      nm_unregister_value(val);

      if (rcurr && rcurr->key - right.offset(rec) >= result.ref_shape(rec)) rcurr = NULL;
      if (lcurr && lcurr->key - left.offset(rec) >= result.ref_shape(rec)) lcurr = NULL;
    }
    __nm_list_storage_unregister_temp_value_list(temp_vals);
  }

  if (left.dtype() == nm::RUBYOBJ) {
    nm_list_storage_unregister_list(l, rec);
  }
  if (right.dtype() == nm::RUBYOBJ) {
    nm_list_storage_unregister_list(r, rec);
  }
  if (result.dtype() == nm::RUBYOBJ) {
    nm_list_storage_unregister_list(x, rec);
  }
}


/*
 * Recursive function, sets multiple values in a matrix from multiple source values. Also handles removal; returns true
 * if the recursion results in an empty list at that level (which signals that the current parent should be removed).
 */
template <typename D>
static bool slice_set(LIST_STORAGE* dest, LIST* l, size_t* coords, size_t* lengths, size_t n, D* v, size_t v_size, size_t& v_offset) {
  using nm::list::node_is_within_slice;
  using nm::list::remove_by_node;
  using nm::list::find_preceding_from_list;
  using nm::list::insert_first_list;
  using nm::list::insert_first_node;
  using nm::list::insert_after;
  size_t* offsets = dest->offset;

  nm_list_storage_register(dest);
  if (dest->dtype == nm::RUBYOBJ) {
    nm_register_values(reinterpret_cast<VALUE*>(v), v_size);
    nm_list_storage_register_list(l, dest->dim - n - 1);
  }

  // drill down into the structure
  NODE* prev = find_preceding_from_list(l, coords[n] + offsets[n]);
  NODE* node = NULL;
  if (prev) node = prev->next && node_is_within_slice(prev->next, coords[n] + offsets[n], lengths[n]) ? prev->next : NULL;
  else      node = node_is_within_slice(l->first, coords[n] + offsets[n], lengths[n]) ? l->first : NULL;

  if (dest->dim - n > 1) {
    size_t i    = 0;
    size_t key  = i + offsets[n] + coords[n];

    // Make sure we have an element to work with
    if (!node) {
      if (!prev) {
        node = insert_first_list(l, key, nm::list::create());
      } else {
        node = insert_after(prev, key, nm::list::create());
      }
    }

    // At this point, it's guaranteed that there is a list here matching key.
    std::list<LIST*> temp_lists;
    while (node) {
      // Recurse down into the list. If it returns true, it's empty, so we need to delete it.
      bool remove_parent = slice_set(dest, reinterpret_cast<LIST*>(node->val), coords, lengths, n+1, v, v_size, v_offset);
      if (dest->dtype == nm::RUBYOBJ) {
        temp_lists.push_front(reinterpret_cast<LIST*>(node->val));
        nm_list_storage_register_list(reinterpret_cast<LIST*>(node->val), dest->dim - n - 2);
      }
      if (remove_parent) {
        NM_FREE(remove_by_node(l, prev, node));
        if (prev) node = prev->next ? prev->next : NULL;
        else      node = l->first   ? l->first   : NULL;
      } else {  // move forward
        prev = node;
        node = node_is_within_slice(prev->next, key-i, lengths[n]) ? prev->next : NULL;
      }

      ++i; ++key;

      if (i >= lengths[n]) break;

      // Now do we need to insert another node here? Or is there already one?
      if (!node) {
        if (!prev) {
          node = insert_first_list(l, key, nm::list::create());
        } else {
          node = insert_after(prev, key, nm::list::create());
        }
      }
    }
    __nm_list_storage_unregister_temp_list_list(temp_lists, dest->dim - n - 2);

  } else {

    size_t i    = 0;
    size_t key  = i + offsets[n] + coords[n];
    std::list<VALUE*> temp_vals;
    while (i < lengths[n]) {
      // Make sure we have an element to work with
      if (v_offset >= v_size) v_offset %= v_size;

      if (node) {
        if (node->key == key) {
          if (v[v_offset] == *reinterpret_cast<D*>(dest->default_val)) { // remove zero value

            NM_FREE(remove_by_node(l, (prev ? prev : l->first), node));

            if (prev) node = prev->next ? prev->next : NULL;
            else      node = l->first   ? l->first   : NULL;

          } else { // edit directly
            *reinterpret_cast<D*>(node->val) = v[v_offset];
            prev = node;
            node = node->next ? node->next : NULL;
          }
        } else if (node->key > key) {
          D* nv = NM_ALLOC(D); *nv = v[v_offset++];
          if (dest->dtype == nm::RUBYOBJ) {
            nm_register_value(*reinterpret_cast<VALUE*>(nv));
            temp_vals.push_front(reinterpret_cast<VALUE*>(nv));
          }

          if (prev) node = insert_after(prev, key, nv);
          else      node = insert_first_node(l, key, nv, sizeof(D));

          prev = node;
          node = prev->next ? prev->next : NULL;
        }
      } else { // no node -- insert a new one
        D* nv = NM_ALLOC(D); *nv = v[v_offset++];
        if (dest->dtype == nm::RUBYOBJ) {
          nm_register_value(*reinterpret_cast<VALUE*>(nv));
          temp_vals.push_front(reinterpret_cast<VALUE*>(nv));
        }
        if (prev) node = insert_after(prev, key, nv);
        else      node = insert_first_node(l, key, nv, sizeof(D));

        prev = node;
        node = prev->next ? prev->next : NULL;
      }

      ++i; ++key;
    }
    __nm_list_storage_unregister_temp_value_list(temp_vals);
  }

  if (dest->dtype == nm::RUBYOBJ) {
    nm_unregister_values(reinterpret_cast<VALUE*>(v), v_size);
    nm_list_storage_unregister_list(l, dest->dim - n - 1);
  }
  nm_list_storage_unregister(dest);

  return (l->first) ? false : true;
}


template <typename D>
void set(VALUE left, SLICE* slice, VALUE right) {
  NM_CONSERVATIVE(nm_register_value(left));
  NM_CONSERVATIVE(nm_register_value(right));
  LIST_STORAGE* s = NM_STORAGE_LIST(left);
  
  std::pair<NMATRIX*,bool> nm_and_free =
    interpret_arg_as_dense_nmatrix(right, NM_DTYPE(left));

  // Map the data onto D* v.
  D*     v;
  size_t v_size = 1;

  if (nm_and_free.first) {
    DENSE_STORAGE* t = reinterpret_cast<DENSE_STORAGE*>(nm_and_free.first->storage);
    v                = reinterpret_cast<D*>(t->elements);
    v_size           = nm_storage_count_max_elements(t);

  } else if (TYPE(right) == T_ARRAY) {
    nm_register_nmatrix(nm_and_free.first);
    v_size = RARRAY_LEN(right);
    v      = NM_ALLOC_N(D, v_size);
    if (NM_DTYPE(left) == nm::RUBYOBJ)
        nm_register_values(reinterpret_cast<VALUE*>(v), v_size);

    for (size_t m = 0; m < v_size; ++m) {
      rubyval_to_cval(rb_ary_entry(right, m), s->dtype, &(v[m]));
    }
    if (NM_DTYPE(left) == nm::RUBYOBJ)
        nm_unregister_values(reinterpret_cast<VALUE*>(v), v_size);

  } else {
    nm_register_nmatrix(nm_and_free.first);
    v = reinterpret_cast<D*>(rubyobj_to_cval(right, NM_DTYPE(left)));
  }

  if (v_size == 1 && *v == *reinterpret_cast<D*>(s->default_val)) {
    if (*reinterpret_cast<D*>(nm_list_storage_get(s, slice)) != *reinterpret_cast<D*>(s->default_val)) {
      nm::list::remove_recursive(s->rows, slice->coords, s->offset, slice->lengths, 0, s->dim);
    }
  } else if (slice->single) {
    slice_set_single(s, s->rows, reinterpret_cast<void*>(v), slice->coords, slice->lengths, 0);
  } else {
    size_t v_offset = 0;
    slice_set<D>(s, s->rows, slice->coords, slice->lengths, 0, v, v_size, v_offset);
  }


  // Only free v if it was allocated in this function.
  if (nm_and_free.first) {
    if (nm_and_free.second) {
      nm_delete(nm_and_free.first);
    }
  } else {
    NM_FREE(v);
    nm_unregister_nmatrix(nm_and_free.first);
  }
  NM_CONSERVATIVE(nm_unregister_value(left));
  NM_CONSERVATIVE(nm_unregister_value(right));
}

/*
 * Used only to set a default initial value.
 */
template <typename D>
void init_default(LIST_STORAGE* s) {
  s->default_val = NM_ALLOC(D);
  *reinterpret_cast<D*>(s->default_val) = 0;
}


}} // end of namespace list_storage

extern "C" {

/*
 * Functions
 */


////////////////
// Lifecycle //
///////////////


/*
 * Creates a list-of-lists(-of-lists-of-lists-etc) storage framework for a
 * matrix.
 *
 * Note: The pointers you pass in for shape and init_val become property of our
 * new storage. You don't need to free them, and you shouldn't re-use them.
 */
LIST_STORAGE* nm_list_storage_create(nm::dtype_t dtype, size_t* shape, size_t dim, void* init_val) {
  LIST_STORAGE* s = NM_ALLOC( LIST_STORAGE );

  s->dim   = dim;
  s->shape = shape;
  s->dtype = dtype;

  s->offset = NM_ALLOC_N(size_t, s->dim);
  memset(s->offset, 0, s->dim * sizeof(size_t));

  s->rows  = nm::list::create();
  if (init_val)
    s->default_val = init_val;
  else {
    DTYPE_TEMPLATE_TABLE(nm::list_storage::init_default, void, LIST_STORAGE*)
    ttable[dtype](s);
  }
  s->count = 1;
  s->src = s;

  return s;
}

/*
 * Destructor for list storage.
 */
void nm_list_storage_delete(STORAGE* s) {
  if (s) {
    LIST_STORAGE* storage = (LIST_STORAGE*)s;
    if (storage->count-- == 1) {
      nm::list::del( storage->rows, storage->dim - 1 );

      NM_FREE(storage->shape);
      NM_FREE(storage->offset);
      NM_FREE(storage->default_val);
      NM_FREE(s);
    }
  }
}

/*
 * Destructor for a list storage reference slice.
 */
void nm_list_storage_delete_ref(STORAGE* s) {
  if (s) {
    LIST_STORAGE* storage = (LIST_STORAGE*)s;

    nm_list_storage_delete( reinterpret_cast<STORAGE*>(storage->src ) );
    NM_FREE(storage->shape);
    NM_FREE(storage->offset);
    NM_FREE(s);
  }
}

/*
 * GC mark function for list storage.
 */
void nm_list_storage_mark(STORAGE* storage_base) {
  LIST_STORAGE* storage = (LIST_STORAGE*)storage_base;

  if (storage && storage->dtype == nm::RUBYOBJ) {
    rb_gc_mark(*((VALUE*)(storage->default_val)));
    nm::list::mark(storage->rows, storage->dim - 1);
  }
}

static void __nm_list_storage_unregister_temp_value_list(std::list<VALUE*>& temp_vals) {
  for (std::list<VALUE*>::iterator it = temp_vals.begin(); it != temp_vals.end(); ++it) {
    nm_unregister_value(**it);
  }
}

static void __nm_list_storage_unregister_temp_list_list(std::list<LIST*>& temp_vals, size_t recursions) {
  for (std::list<LIST*>::iterator it = temp_vals.begin(); it != temp_vals.end(); ++it) {
    nm_list_storage_unregister_list(*it, recursions);
  }
}

void nm_list_storage_register_node(const NODE* curr) {
  nm_register_value(*reinterpret_cast<VALUE*>(curr->val));      
}

void nm_list_storage_unregister_node(const NODE* curr) {
  nm_unregister_value(*reinterpret_cast<VALUE*>(curr->val));      
}

/**
 * Gets rid of all instances of a given node in the registration list.
 * Sometimes a node will get deleted and replaced deep in a recursion, but
 * further up it will still get registered.  This leads to a potential read
 * after free during the GC marking.  This function completely clears out a
 * node so that this won't happen.
 */
void nm_list_storage_completely_unregister_node(const NODE* curr) {
  nm_completely_unregister_value(*reinterpret_cast<VALUE*>(curr->val));
}

void nm_list_storage_register_list(const LIST* list, size_t recursions) {
  NODE* next;
  if (!list) return;
  NODE* curr = list->first;

  while (curr != NULL) {
    next = curr->next;
    if (recursions == 0) {
      nm_list_storage_register_node(curr);
    } else {
      nm_list_storage_register_list(reinterpret_cast<LIST*>(curr->val), recursions - 1);
    }
    curr = next;
  }
}

void nm_list_storage_unregister_list(const LIST* list, size_t recursions) {
  NODE* next;
  if (!list) return;
  NODE* curr = list->first;

  while (curr != NULL) {
    next = curr->next;
    if (recursions == 0) {
      nm_list_storage_unregister_node(curr);
    } else {
      nm_list_storage_unregister_list(reinterpret_cast<LIST*>(curr->val), recursions - 1);
    }
    curr = next;
  }
}

void nm_list_storage_register(const STORAGE* s) {
  const LIST_STORAGE* storage = reinterpret_cast<const LIST_STORAGE*>(s);
  if (storage && storage->dtype == nm::RUBYOBJ) {
    nm_register_value(*reinterpret_cast<VALUE*>(storage->default_val));
    nm_list_storage_register_list(storage->rows, storage->dim - 1);
  }
}

void nm_list_storage_unregister(const STORAGE* s) {
  const LIST_STORAGE* storage = reinterpret_cast<const LIST_STORAGE*>(s);
  if (storage && storage->dtype == nm::RUBYOBJ) {
    nm_unregister_value(*reinterpret_cast<VALUE*>(storage->default_val));
    nm_list_storage_unregister_list(storage->rows, storage->dim - 1);
  }
}

///////////////
// Accessors //
///////////////

/*
 * Documentation goes here.
 */
static NODE* list_storage_get_single_node(LIST_STORAGE* s, SLICE* slice) {
  size_t r;
  LIST*  l = s->rows;
  NODE*  n;

  for (r = 0; r < s->dim; r++) {
    n = nm::list::find(l, s->offset[r] + slice->coords[r]);
    if (n)  l = reinterpret_cast<LIST*>(n->val);
    else return NULL;
  }

  return n;
}


/*
 * Recursive helper function for each_with_indices, based on nm_list_storage_count_elements_r.
 * Handles empty/non-existent sublists.
 */
static void each_empty_with_indices_r(nm::list_storage::RecurseData& s, size_t rec, VALUE& stack) {
  VALUE empty  = s.dtype() == nm::RUBYOBJ ? *reinterpret_cast<VALUE*>(s.init()) : s.init_obj();
  NM_CONSERVATIVE(nm_register_value(stack));

  if (rec) {
    for (long index = 0; index < s.ref_shape(rec); ++index) {
      // Don't do an unshift/shift here -- we'll let that be handled in the lowest-level iteration (recursions == 0)
      rb_ary_push(stack, LONG2NUM(index));
      each_empty_with_indices_r(s, rec-1, stack);
      rb_ary_pop(stack);
    }
  } else {
    rb_ary_unshift(stack, empty);
    for (long index = 0; index < s.ref_shape(rec); ++index) {
      rb_ary_push(stack, LONG2NUM(index));
      rb_yield_splat(stack);
      rb_ary_pop(stack);
    }
    rb_ary_shift(stack);
  }
  NM_CONSERVATIVE(nm_unregister_value(stack));
}

/*
 * Recursive helper function for each_with_indices, based on nm_list_storage_count_elements_r.
 */
static void each_with_indices_r(nm::list_storage::RecurseData& s, const LIST* l, size_t rec, VALUE& stack) {
  if (s.dtype() == nm::RUBYOBJ)
    nm_list_storage_register_list(l, rec);
  NM_CONSERVATIVE(nm_register_value(stack));
  NODE*  curr  = l->first;

  size_t offset = s.offset(rec);
  size_t shape  = s.ref_shape(rec);

  while (curr && curr->key < offset) curr = curr->next;
  if (curr && curr->key - offset >= shape) curr = NULL;


  if (rec) {
    for (long index = 0; index < shape; ++index) { // index in reference
      rb_ary_push(stack, LONG2NUM(index));
      if (!curr || index < curr->key - offset) {
        each_empty_with_indices_r(s, rec-1, stack);
      } else { // index == curr->key - offset
        each_with_indices_r(s, reinterpret_cast<const LIST*>(curr->val), rec-1, stack);
        curr = curr->next;
      }
      rb_ary_pop(stack);
    }
  } else {
    for (long index = 0; index < shape; ++index) {

      rb_ary_push(stack, LONG2NUM(index));

      if (!curr || index < curr->key - offset) {
        rb_ary_unshift(stack, s.dtype() == nm::RUBYOBJ ? *reinterpret_cast<VALUE*>(s.init()) : s.init_obj());

      } else { // index == curr->key - offset
        rb_ary_unshift(stack, s.dtype() == nm::RUBYOBJ ? *reinterpret_cast<VALUE*>(curr->val) : rubyobj_from_cval(curr->val, s.dtype()).rval);

        curr = curr->next;
      }
      rb_yield_splat(stack);

      rb_ary_shift(stack);
      rb_ary_pop(stack);
    }
  }
  NM_CONSERVATIVE(nm_unregister_value(stack));
  if (s.dtype() == nm::RUBYOBJ)
    nm_list_storage_unregister_list(l, rec);
}


/*
 * Recursive helper function for each_stored_with_indices, based on nm_list_storage_count_elements_r.
 */
static void each_stored_with_indices_r(nm::list_storage::RecurseData& s, const LIST* l, size_t rec, VALUE& stack) {
  if (s.dtype() == nm::RUBYOBJ)
    nm_list_storage_register_list(l, rec);
  NM_CONSERVATIVE(nm_register_value(stack));
  
  NODE* curr = l->first;

  size_t offset = s.offset(rec);
  size_t shape  = s.ref_shape(rec);

  while (curr && curr->key < offset) { curr = curr->next; }
  if (curr && curr->key - offset >= shape) curr = NULL;

  if (rec) {
    while (curr) {

      rb_ary_push(stack, LONG2NUM(static_cast<long>(curr->key - offset)));
      each_stored_with_indices_r(s, reinterpret_cast<const LIST*>(curr->val), rec-1, stack);
      rb_ary_pop(stack);

      curr = curr->next;
      if (curr && curr->key - offset >= shape) curr = NULL;
    }
  } else {
    while (curr) {
      rb_ary_push(stack, LONG2NUM(static_cast<long>(curr->key - offset))); // add index to end

      // add value to beginning
      rb_ary_unshift(stack, s.dtype() == nm::RUBYOBJ ? *reinterpret_cast<VALUE*>(curr->val) : rubyobj_from_cval(curr->val, s.dtype()).rval);
      // yield to the whole stack (value, i, j, k, ...)
      rb_yield_splat(stack);

      // remove the value
      rb_ary_shift(stack);

      // remove the index from the end
      rb_ary_pop(stack);

      curr = curr->next;
      if (curr && curr->key - offset >= shape) curr = NULL;
    }
  }
  NM_CONSERVATIVE(nm_unregister_value(stack));
  if (s.dtype() == nm::RUBYOBJ)
    nm_list_storage_unregister_list(l, rec);
}



/*
 * Each/each-stored iterator, brings along the indices.
 */
VALUE nm_list_each_with_indices(VALUE nmatrix, bool stored) {

  NM_CONSERVATIVE(nm_register_value(nmatrix));

  // If we don't have a block, return an enumerator.
  RETURN_SIZED_ENUMERATOR_PRE
  NM_CONSERVATIVE(nm_unregister_value(nmatrix));
  RETURN_SIZED_ENUMERATOR(nmatrix, 0, 0, 0);

  nm::list_storage::RecurseData sdata(NM_STORAGE_LIST(nmatrix));

  VALUE stack = rb_ary_new();

  if (stored) each_stored_with_indices_r(sdata, sdata.top_level_list(), sdata.dim() - 1, stack);
  else        each_with_indices_r(sdata, sdata.top_level_list(), sdata.dim() - 1, stack);

  NM_CONSERVATIVE(nm_unregister_value(nmatrix));
  return nmatrix;
}


/*
 * map merged stored iterator. Always returns a matrix containing RubyObjects which probably needs to be casted.
 */
VALUE nm_list_map_stored(VALUE left, VALUE init) {
  NM_CONSERVATIVE(nm_register_value(left));
  NM_CONSERVATIVE(nm_register_value(init));

  bool scalar = false;

  LIST_STORAGE *s   = NM_STORAGE_LIST(left);

  // For each matrix, if it's a reference, we want to deal directly with the original (with appropriate offsetting)
  nm::list_storage::RecurseData sdata(s);

  void* scalar_init = NULL;

  //if (!rb_block_given_p()) {
  //  rb_raise(rb_eNotImpError, "RETURN_SIZED_ENUMERATOR probably won't work for a map_merged since no merged object is created");
  //}
  // If we don't have a block, return an enumerator.
  RETURN_SIZED_ENUMERATOR_PRE
  NM_CONSERVATIVE(nm_unregister_value(left));
  NM_CONSERVATIVE(nm_unregister_value(init));
  RETURN_SIZED_ENUMERATOR(left, 0, 0, 0); // FIXME: Test this. Probably won't work. Enable above code instead.

  // Figure out default value if none provided by the user
  if (init == Qnil) {
    nm_unregister_value(init);
    init = rb_yield_values(1, sdata.init_obj());
    nm_register_value(init);
  }
	// Allocate a new shape array for the resulting matrix.
  void* init_val = NM_ALLOC(VALUE);
  memcpy(init_val, &init, sizeof(VALUE));
  nm_register_value(*reinterpret_cast<VALUE*>(init_val));

  NMATRIX* result = nm_create(nm::LIST_STORE, nm_list_storage_create(nm::RUBYOBJ, sdata.copy_alloc_shape(), s->dim, init_val));
  LIST_STORAGE* r = reinterpret_cast<LIST_STORAGE*>(result->storage);
  nm::list_storage::RecurseData rdata(r, init);
  nm_register_nmatrix(result);
  map_stored_r(rdata, sdata, rdata.top_level_list(), sdata.top_level_list(), sdata.dim() - 1);

  VALUE to_return = Data_Wrap_Struct(CLASS_OF(left), nm_mark, nm_delete, result);

  nm_unregister_nmatrix(result);
  nm_unregister_value(*reinterpret_cast<VALUE*>(init_val));
  NM_CONSERVATIVE(nm_unregister_value(init));
  NM_CONSERVATIVE(nm_unregister_value(left));

  return to_return;
}


/*
 * map merged stored iterator. Always returns a matrix containing RubyObjects which probably needs to be casted.
 */
VALUE nm_list_map_merged_stored(VALUE left, VALUE right, VALUE init) {
  NM_CONSERVATIVE(nm_register_value(left));
  NM_CONSERVATIVE(nm_register_value(right));
  NM_CONSERVATIVE(nm_register_value(init));

  bool scalar = false;

  LIST_STORAGE *s   = NM_STORAGE_LIST(left),
               *t;

  // For each matrix, if it's a reference, we want to deal directly with the original (with appropriate offsetting)
  nm::list_storage::RecurseData sdata(s);

  void* scalar_init = NULL;

  // right might be a scalar, in which case this is a scalar operation.
  if (TYPE(right) != T_DATA || (RDATA(right)->dfree != (RUBY_DATA_FUNC)nm_delete && RDATA(right)->dfree != (RUBY_DATA_FUNC)nm_delete_ref)) {
    nm::dtype_t r_dtype = Upcast[NM_DTYPE(left)][nm_dtype_min(right)];
    scalar_init         = rubyobj_to_cval(right, r_dtype); // make a copy of right

    t                   = reinterpret_cast<LIST_STORAGE*>(nm_list_storage_create(r_dtype, sdata.copy_alloc_shape(), s->dim, scalar_init));
    scalar              = true;
  } else {
    t                   = NM_STORAGE_LIST(right); // element-wise, not scalar.
  }

  //if (!rb_block_given_p()) {
  //  rb_raise(rb_eNotImpError, "RETURN_SIZED_ENUMERATOR probably won't work for a map_merged since no merged object is created");
  //}
  // If we don't have a block, return an enumerator.
  RETURN_SIZED_ENUMERATOR_PRE
  NM_CONSERVATIVE(nm_unregister_value(left));
  NM_CONSERVATIVE(nm_unregister_value(right));
  NM_CONSERVATIVE(nm_unregister_value(init));
  RETURN_SIZED_ENUMERATOR(left, 0, 0, 0); // FIXME: Test this. Probably won't work. Enable above code instead.

  // Figure out default value if none provided by the user
  nm::list_storage::RecurseData& tdata = *(new nm::list_storage::RecurseData(t)); //FIXME: this is a hack to make sure that we can run the destructor before nm_list_storage_delete(t) below.
  if (init == Qnil) {
    nm_unregister_value(init);
    init = rb_yield_values(2, sdata.init_obj(), tdata.init_obj());
    nm_register_value(init);
  }

  // Allocate a new shape array for the resulting matrix.
  void* init_val = NM_ALLOC(VALUE);
  memcpy(init_val, &init, sizeof(VALUE));
  nm_register_value(*reinterpret_cast<VALUE*>(init_val));

  NMATRIX* result = nm_create(nm::LIST_STORE, nm_list_storage_create(nm::RUBYOBJ, sdata.copy_alloc_shape(), s->dim, init_val));
  LIST_STORAGE* r = reinterpret_cast<LIST_STORAGE*>(result->storage);
  nm::list_storage::RecurseData rdata(r, init);
  map_merged_stored_r(rdata, sdata, tdata, rdata.top_level_list(), sdata.top_level_list(), tdata.top_level_list(), sdata.dim() - 1);

  delete &tdata;
  // If we are working with a scalar operation
  if (scalar) nm_list_storage_delete(t);

  VALUE to_return = Data_Wrap_Struct(CLASS_OF(left), nm_mark, nm_delete, result);

  nm_unregister_value(*reinterpret_cast<VALUE*>(init_val));

  NM_CONSERVATIVE(nm_unregister_value(init));
  NM_CONSERVATIVE(nm_unregister_value(right));
  NM_CONSERVATIVE(nm_unregister_value(left));

  return to_return;
}


/*
 * Copy a slice of a list matrix into a regular list matrix.
 */
static LIST* slice_copy(const LIST_STORAGE* src, LIST* src_rows, size_t* coords, size_t* lengths, size_t n) {
  nm_list_storage_register(src);
  void *val = NULL;
  int key;
  
  LIST* dst_rows = nm::list::create();
  NODE* src_node = src_rows->first;
  std::list<VALUE*> temp_vals;
  std::list<LIST*> temp_lists;
  while (src_node) {
    key = src_node->key - (src->offset[n] + coords[n]);
    
    if (key >= 0 && (size_t)key < lengths[n]) {
      if (src->dim - n > 1) {
        val = slice_copy( src,
                          reinterpret_cast<LIST*>(src_node->val),
                          coords,
                          lengths,
                          n + 1    );
        if (val) {
          if (src->dtype == nm::RUBYOBJ) {
            nm_list_storage_register_list(reinterpret_cast<LIST*>(val), src->dim - n - 2);
            temp_lists.push_front(reinterpret_cast<LIST*>(val));
          }
          nm::list::insert_copy(dst_rows, false, key, val, sizeof(LIST));
        }
      } else { // matches src->dim - n > 1
        if (src->dtype == nm::RUBYOBJ) {
          nm_register_value(*reinterpret_cast<VALUE*>(src_node->val));
          temp_vals.push_front(reinterpret_cast<VALUE*>(src_node->val));
        }
        nm::list::insert_copy(dst_rows, false, key, src_node->val, DTYPE_SIZES[src->dtype]);
      }
    }
    src_node = src_node->next;
 }
  if (src->dtype == nm::RUBYOBJ) {
    __nm_list_storage_unregister_temp_list_list(temp_lists, src->dim - n - 2);
    __nm_list_storage_unregister_temp_value_list(temp_vals);
  }
  nm_list_storage_unregister(src);
  return dst_rows;
}

/*
 * Documentation goes here.
 */
void* nm_list_storage_get(const STORAGE* storage, SLICE* slice) {
  LIST_STORAGE* s = (LIST_STORAGE*)storage;
  LIST_STORAGE* ns = NULL;

  nm_list_storage_register(s);

  if (slice->single) {
    NODE* n = list_storage_get_single_node(s, slice);
    nm_list_storage_unregister(s);
    return (n ? n->val : s->default_val);

  } else {
    void *init_val = NM_ALLOC_N(char, DTYPE_SIZES[s->dtype]);
    memcpy(init_val, s->default_val, DTYPE_SIZES[s->dtype]);
    if (s->dtype == nm::RUBYOBJ)
      nm_register_value(*reinterpret_cast<VALUE*>(init_val));

    size_t *shape = NM_ALLOC_N(size_t, s->dim);
    memcpy(shape, slice->lengths, sizeof(size_t) * s->dim);

    ns = nm_list_storage_create(s->dtype, shape, s->dim, init_val);
  
    ns->rows = slice_copy(s, s->rows, slice->coords, slice->lengths, 0);

    if (s->dtype == nm::RUBYOBJ)
      nm_unregister_value(*reinterpret_cast<VALUE*>(init_val));
    nm_list_storage_unregister(s);

    return ns;
  }
}

/*
 * Get the contents of some set of coordinates. Note: Does not make a copy!
 * Don't free!
 */
void* nm_list_storage_ref(const STORAGE* storage, SLICE* slice) {
  LIST_STORAGE* s = (LIST_STORAGE*)storage;
  LIST_STORAGE* ns = NULL;
  nm_list_storage_register(s);

  //TODO: It needs a refactoring.
  if (slice->single) {
    NODE* n = list_storage_get_single_node(s, slice);
    nm_list_storage_unregister(s);
    return (n ? n->val : s->default_val);
  } 
  else {
    ns              = NM_ALLOC( LIST_STORAGE );
    
    ns->dim         = s->dim;
    ns->dtype       = s->dtype;
    ns->offset      = NM_ALLOC_N(size_t, ns->dim);
    ns->shape       = NM_ALLOC_N(size_t, ns->dim);

    for (size_t i = 0; i < ns->dim; ++i) {
      ns->offset[i] = slice->coords[i] + s->offset[i];
      ns->shape[i]  = slice->lengths[i];
    }

    ns->rows        = s->rows;
    ns->default_val = s->default_val;
    
    s->src->count++;
    ns->src         = s->src;
    nm_list_storage_unregister(s);
    return ns;
  }
}


/*
 * Recursive function, sets multiple values in a matrix from a single source value.
 */
static void slice_set_single(LIST_STORAGE* dest, LIST* l, void* val, size_t* coords, size_t* lengths, size_t n) {
  nm_list_storage_register(dest);
  if (dest->dtype == nm::RUBYOBJ) {
    nm_register_value(*reinterpret_cast<VALUE*>(val));
    nm_list_storage_register_list(l, dest->dim - n - 1);
  }

  // drill down into the structure
  NODE* node = NULL;
  if (dest->dim - n > 1) {
    std::list<LIST*> temp_nodes; 
    for (size_t i = 0; i < lengths[n]; ++i) {

      size_t key = i + dest->offset[n] + coords[n];

      if (!node) {
        node = nm::list::insert(l, false, key, nm::list::create()); // try to insert list
      } else if (!node->next || (node->next && node->next->key > key)) {
        node = nm::list::insert_after(node, key, nm::list::create());
      } else {
        node = node->next; // correct rank already exists.
      }

      if (dest->dtype == nm::RUBYOBJ) {
        temp_nodes.push_front(reinterpret_cast<LIST*>(node->val));
        nm_list_storage_register_list(reinterpret_cast<LIST*>(node->val), dest->dim - n - 2);
      }

      // cast it to a list and recurse
      slice_set_single(dest, reinterpret_cast<LIST*>(node->val), val, coords, lengths, n + 1);
    }
    __nm_list_storage_unregister_temp_list_list(temp_nodes, dest->dim - n - 2);
  } else {
    std::list<VALUE*> temp_vals;
    for (size_t i = 0; i < lengths[n]; ++i) {

      size_t key = i + dest->offset[n] + coords[n];

      if (!node)  {
        node = nm::list::insert_copy(l, true, key, val, DTYPE_SIZES[dest->dtype]);
      } else {
        node = nm::list::replace_insert_after(node, key, val, true, DTYPE_SIZES[dest->dtype]);
      }
      if (dest->dtype == nm::RUBYOBJ) {
        temp_vals.push_front(reinterpret_cast<VALUE*>(node->val));
        nm_register_value(*reinterpret_cast<VALUE*>(node->val));
      }
    }
    __nm_list_storage_unregister_temp_value_list(temp_vals);
  }

  nm_list_storage_unregister(dest);
  if (dest->dtype == nm::RUBYOBJ) {
    nm_unregister_value(*reinterpret_cast<VALUE*>(val));
    nm_list_storage_unregister_list(l, dest->dim - n - 1);
  }
}



/*
 * Set a value or values in a list matrix.
 */
void nm_list_storage_set(VALUE left, SLICE* slice, VALUE right) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::list_storage::set, void, VALUE, SLICE*, VALUE)
  ttable[NM_DTYPE(left)](left, slice, right);
}


/*
 * Insert an entry directly in a row (not using copy! don't free after).
 *
 * Returns a pointer to the insertion location.
 *
 * TODO: Allow this function to accept an entire row and not just one value -- for slicing
 */
NODE* nm_list_storage_insert(STORAGE* storage, SLICE* slice, void* val) {
  LIST_STORAGE* s = (LIST_STORAGE*)storage;
  nm_list_storage_register(s);
  if (s->dtype == nm::RUBYOBJ)
    nm_register_value(*reinterpret_cast<VALUE*>(val));
  // Pretend dims = 2
  // Then coords is going to be size 2
  // So we need to find out if some key already exists
  size_t r;
  NODE*  n;
  LIST*  l = s->rows;

  // drill down into the structure
  for (r = 0; r < s->dim -1; ++r) {
    n = nm::list::insert(l, false, s->offset[r] + slice->coords[s->dim - r], nm::list::create());
    l = reinterpret_cast<LIST*>(n->val);
  }

  nm_list_storage_unregister(s);
  if (s->dtype == nm::RUBYOBJ)
    nm_unregister_value(*reinterpret_cast<VALUE*>(val));

  return nm::list::insert(l, true, s->offset[r] + slice->coords[r], val);
}

/*
 * Remove an item or slice from list storage.
 */
void nm_list_storage_remove(STORAGE* storage, SLICE* slice) {
  LIST_STORAGE* s = (LIST_STORAGE*)storage;

  // This returns a boolean, which will indicate whether s->rows is empty.
  // We can safely ignore it, since we never want to delete s->rows until
  // it's time to destroy the LIST_STORAGE object.
  nm::list::remove_recursive(s->rows, slice->coords, s->offset, slice->lengths, 0, s->dim);
}

///////////
// Tests //
///////////

/*
 * Comparison of contents for list storage.
 */
bool nm_list_storage_eqeq(const STORAGE* left, const STORAGE* right) {
	NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::list_storage::eqeq_r, bool, nm::list_storage::RecurseData& left, nm::list_storage::RecurseData& right, const LIST* l, const LIST* r, size_t rec)

  nm::list_storage::RecurseData ldata(reinterpret_cast<const LIST_STORAGE*>(left)),
                                rdata(reinterpret_cast<const LIST_STORAGE*>(right));

	return ttable[left->dtype][right->dtype](ldata, rdata, ldata.top_level_list(), rdata.top_level_list(), ldata.dim()-1);
}

//////////
// Math //
//////////


/*
 * List storage matrix multiplication.
 */
STORAGE* nm_list_storage_matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector) {
  free(resulting_shape);
  rb_raise(rb_eNotImpError, "multiplication not implemented for list-of-list matrices");
  return NULL;
  //DTYPE_TEMPLATE_TABLE(dense_storage::matrix_multiply, NMATRIX*, STORAGE_PAIR, size_t*, bool);

  //return ttable[reinterpret_cast<DENSE_STORAGE*>(casted_storage.left)->dtype](casted_storage, resulting_shape, vector);
}


/*
 * List storage to Hash conversion. Uses Hashes with default values, so you can continue to pretend
 * it's a sparse matrix.
 */
VALUE nm_list_storage_to_hash(const LIST_STORAGE* s, const nm::dtype_t dtype) {
  nm_list_storage_register(s);
  // Get the default value for the list storage.
  VALUE default_value = rubyobj_from_cval(s->default_val, dtype).rval;
  nm_list_storage_unregister(s);
  // Recursively copy each dimension of the matrix into a nested hash.
  return nm_list_copy_to_hash(s->rows, dtype, s->dim - 1, default_value);
}

/////////////
// Utility //
/////////////

/*
 * Recursively count the non-zero elements in a list storage object.
 */
size_t nm_list_storage_count_elements_r(const LIST* l, size_t recursions) {
  size_t count = 0;
  NODE* curr = l->first;
  
  if (recursions) {
    while (curr) {
      count += nm_list_storage_count_elements_r(reinterpret_cast<const LIST*>(curr->val), recursions - 1);
      curr   = curr->next;
    }
    
  } else {
    while (curr) {
      ++count;
      curr = curr->next;
    }
  }
  
  return count;
}

/*
 * Count non-diagonal non-zero elements.
 */
size_t nm_list_storage_count_nd_elements(const LIST_STORAGE* s) {
  NODE *i_curr, *j_curr;
  size_t count = 0;
  
  if (s->dim != 2) {
  	rb_raise(rb_eNotImpError, "non-diagonal element counting only defined for dim = 2");
  }

  for (i_curr = s->rows->first; i_curr; i_curr = i_curr->next) {
    int i = i_curr->key - s->offset[0];
    if (i < 0 || i >= (int)s->shape[0]) continue;

    for (j_curr = ((LIST*)(i_curr->val))->first; j_curr; j_curr = j_curr->next) {
      int j = j_curr->key - s->offset[1];
      if (j < 0 || j >= (int)s->shape[1]) continue;

      if (i != j)  	++count;
    }
  }
  
  return count;
}

/////////////////////////
// Copying and Casting //
/////////////////////////
//
/*
 * List storage copy constructor C access.
 */

LIST_STORAGE* nm_list_storage_copy(const LIST_STORAGE* rhs) {
  nm_list_storage_register(rhs);
  size_t *shape = NM_ALLOC_N(size_t, rhs->dim);
  memcpy(shape, rhs->shape, sizeof(size_t) * rhs->dim);
  
  void *init_val = NM_ALLOC_N(char, DTYPE_SIZES[rhs->dtype]);
  memcpy(init_val, rhs->default_val, DTYPE_SIZES[rhs->dtype]);

  LIST_STORAGE* lhs = nm_list_storage_create(rhs->dtype, shape, rhs->dim, init_val);
  nm_list_storage_register(lhs);

  lhs->rows = slice_copy(rhs, rhs->rows, lhs->offset, lhs->shape, 0);

  nm_list_storage_unregister(rhs);
  nm_list_storage_unregister(lhs);
  return lhs;
}

/*
 * List storage copy constructor C access with casting.
 */
STORAGE* nm_list_storage_cast_copy(const STORAGE* rhs, nm::dtype_t new_dtype, void* dummy) {
  NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::list_storage::cast_copy, LIST_STORAGE*, const LIST_STORAGE* rhs, nm::dtype_t new_dtype);

  return (STORAGE*)ttable[new_dtype][rhs->dtype]((LIST_STORAGE*)rhs, new_dtype);
}


/*
 * List storage copy constructor for transposing.
 */
STORAGE* nm_list_storage_copy_transposed(const STORAGE* rhs_base) {
  rb_raise(rb_eNotImpError, "list storage transpose not yet implemented");
  return NULL;
}


} // end of extern "C" block


/////////////////////////
// Templated Functions //
/////////////////////////


namespace nm {
namespace list_storage {


/*
 * List storage copy constructor for changing dtypes.
 */
template <typename LDType, typename RDType>
static LIST_STORAGE* cast_copy(const LIST_STORAGE* rhs, dtype_t new_dtype) {
  nm_list_storage_register(rhs);
  // allocate and copy shape
  size_t* shape = NM_ALLOC_N(size_t, rhs->dim);
  memcpy(shape, rhs->shape, rhs->dim * sizeof(size_t));

  // copy default value
  LDType* default_val = NM_ALLOC_N(LDType, 1);
  *default_val = *reinterpret_cast<RDType*>(rhs->default_val);

  LIST_STORAGE* lhs = nm_list_storage_create(new_dtype, shape, rhs->dim, default_val);
  //lhs->rows         = nm::list::create();

  nm_list_storage_register(lhs);
  // TODO: Needs optimization. When matrix is reference it is copped twice.
  if (rhs->src == rhs) 
    nm::list::cast_copy_contents<LDType, RDType>(lhs->rows, rhs->rows, rhs->dim - 1);
  else {
    LIST_STORAGE *tmp = nm_list_storage_copy(rhs);
    nm_list_storage_register(tmp);
    nm::list::cast_copy_contents<LDType, RDType>(lhs->rows, tmp->rows, rhs->dim - 1);
    nm_list_storage_unregister(tmp);
    nm_list_storage_delete(tmp);
  }
  nm_list_storage_unregister(lhs);
  nm_list_storage_unregister(rhs);
  return lhs;
}


/*
 * Recursive helper function for eqeq. Note that we use SDType and TDType instead of L and R because this function
 * is a re-labeling. That is, it can be called in order L,R or order R,L; and we don't want to get confused. So we
 * use S and T to denote first and second passed in.
 */
template <typename SDType, typename TDType>
static bool eqeq_empty_r(RecurseData& s, const LIST* l, size_t rec, const TDType* t_init) {
  NODE* curr  = l->first;

  // For reference matrices, make sure we start in the correct place.
  while (curr && curr->key < s.offset(rec)) {  curr = curr->next;  }
  if (curr && curr->key - s.offset(rec) >= s.ref_shape(rec)) curr = NULL;

  if (rec) {
    while (curr) {
      if (!eqeq_empty_r<SDType,TDType>(s, reinterpret_cast<const LIST*>(curr->val), rec-1, t_init)) return false;
      curr = curr->next;

      if (curr && curr->key - s.offset(rec) >= s.ref_shape(rec)) curr = NULL;
    }
  } else {
    while (curr) {
      if (*reinterpret_cast<SDType*>(curr->val) != *t_init) return false;
      curr = curr->next;

      if (curr && curr->key - s.offset(rec) >= s.ref_shape(rec)) curr = NULL;
    }
  }
  return true;
}



/*
 * Do these two list matrices of the same dtype have exactly the same contents (accounting for default_vals)?
 *
 * This function is recursive.
 */
template <typename LDType, typename RDType>
static bool eqeq_r(RecurseData& left, RecurseData& right, const LIST* l, const LIST* r, size_t rec) {
  NODE *lcurr = l->first,
       *rcurr = r->first;

  // For reference matrices, make sure we start in the correct place.
  while (lcurr && lcurr->key < left.offset(rec)) {  lcurr = lcurr->next;  }
  while (rcurr && rcurr->key < right.offset(rec)) {  rcurr = rcurr->next;  }
  if (rcurr && rcurr->key - right.offset(rec) >= left.ref_shape(rec)) rcurr = NULL;
  if (lcurr && lcurr->key - left.offset(rec) >= left.ref_shape(rec)) lcurr = NULL;

  bool compared = false;

  if (rec) {

    while (lcurr || rcurr) {

      if (!rcurr || (lcurr && (lcurr->key - left.offset(rec) < rcurr->key - right.offset(rec)))) {
        if (!eqeq_empty_r<LDType,RDType>(left, reinterpret_cast<const LIST*>(lcurr->val), rec-1, reinterpret_cast<const RDType*>(right.init()))) return false;
        lcurr   = lcurr->next;
      } else if (!lcurr || (rcurr && (rcurr->key - right.offset(rec) < lcurr->key - left.offset(rec)))) {
        if (!eqeq_empty_r<RDType,LDType>(right, reinterpret_cast<const LIST*>(rcurr->val), rec-1, reinterpret_cast<const LDType*>(left.init()))) return false;
        rcurr   = rcurr->next;
      } else { // keys are == and both present
        if (!eqeq_r<LDType,RDType>(left, right, reinterpret_cast<const LIST*>(lcurr->val), reinterpret_cast<const LIST*>(rcurr->val), rec-1)) return false;
        lcurr   = lcurr->next;
        rcurr   = rcurr->next;
      }
      if (rcurr && rcurr->key - right.offset(rec) >= right.ref_shape(rec)) rcurr = NULL;
      if (lcurr && lcurr->key - left.offset(rec)  >= left.ref_shape(rec)) lcurr = NULL;
      compared = true;
    }
  } else {
    while (lcurr || rcurr) {

      if (rcurr && rcurr->key - right.offset(rec) >= left.ref_shape(rec)) rcurr = NULL;
      if (lcurr && lcurr->key - left.offset(rec) >= left.ref_shape(rec)) lcurr = NULL;

      if (!rcurr || (lcurr && (lcurr->key - left.offset(rec) < rcurr->key - right.offset(rec)))) {
        if (*reinterpret_cast<LDType*>(lcurr->val) != *reinterpret_cast<const RDType*>(right.init())) return false;
        lcurr         = lcurr->next;
      } else if (!lcurr || (rcurr && (rcurr->key - right.offset(rec) < lcurr->key - left.offset(rec)))) {
        if (*reinterpret_cast<RDType*>(rcurr->val) != *reinterpret_cast<const LDType*>(left.init())) return false;
        rcurr         = rcurr->next;
      } else { // keys == and both left and right nodes present
        if (*reinterpret_cast<LDType*>(lcurr->val) != *reinterpret_cast<RDType*>(rcurr->val)) return false;
        lcurr         = lcurr->next;
        rcurr         = rcurr->next;
      }
      if (rcurr && rcurr->key - right.offset(rec) >= right.ref_shape(rec)) rcurr = NULL;
      if (lcurr && lcurr->key - left.offset(rec)  >= left.ref_shape(rec)) lcurr = NULL;
      compared = true;
    }
  }

  // Final condition: both containers are empty, and have different default values.
  if (!compared && !lcurr && !rcurr) return *reinterpret_cast<const LDType*>(left.init()) == *reinterpret_cast<const RDType*>(right.init());
  return true;
}


}} // end of namespace nm::list_storage

extern "C" {
  /*
   * call-seq:
   *     __list_to_hash__ -> Hash
   *
   * Create a Ruby Hash from a list NMatrix.
   *
   * This is an internal C function which handles list stype only.
   */
  VALUE nm_to_hash(VALUE self) {
    return nm_list_storage_to_hash(NM_STORAGE_LIST(self), NM_DTYPE(self));
  }

  /*
   * call-seq:
   *     __list_default_value__ -> ...
   *
   * Get the default_value property from a list matrix.
   */
  VALUE nm_list_default_value(VALUE self) {
    NM_CONSERVATIVE(nm_register_value(self));
    VALUE to_return = (NM_DTYPE(self) == nm::RUBYOBJ) ? *reinterpret_cast<VALUE*>(NM_DEFAULT_VAL(self)) : rubyobj_from_cval(NM_DEFAULT_VAL(self), NM_DTYPE(self)).rval;
    NM_CONSERVATIVE(nm_unregister_value(self));
    return to_return;
  }
} // end of extern "C" block
