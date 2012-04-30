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
// SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
// NMatrix is Copyright (c) 2012, Ruby Science Foundation
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

#ifndef LIST_C
# define LIST_C

#include <ruby.h>

#include "nmatrix.h"

extern VALUE nm_eStorageTypeError;
extern nm_eqeq_t ElemEqEq;


/* Calculate the max number of elements in the list storage structure, based on shape and rank */
inline size_t count_storage_max_elements(const STORAGE* s) {
  return count_dense_storage_elements((DENSE_STORAGE*)s);
}

static size_t count_list_storage_elements_r(const LIST* l, size_t recursions) {
  size_t count = 0;
  NODE* curr = l->first;
  if (recursions) {
    while (curr) {
      count += count_list_storage_elements_r(curr->val, recursions-1);
      curr = curr->next;
    }
  } else {
    while (curr) {
      ++count;
      curr = curr->next;
    }
  }
  return count;
}


// Count non-zero elements. See also count_list_storage_nd_elements.
size_t count_list_storage_elements(const LIST_STORAGE* s) {
  return count_list_storage_elements_r(s->rows, s->rank-1);
}


// Count non-diagonal non-zero elements
size_t count_list_storage_nd_elements(const LIST_STORAGE* s) {
  NODE *i_curr, *j_curr;
  size_t count = 0;
  if (s->rank != 2) rb_raise(rb_eNotImpError, "non-diagonal element counting only defined for rank = 2");

  for (i_curr = s->rows->first; i_curr; i_curr = i_curr->next) {
    for (j_curr = ((LIST*)(i_curr->val))->first; j_curr; j_curr = j_curr->next) {
      if (i_curr->key != j_curr->key) ++count;
    }
  }
  return count;
}


/* Finds the node that should go before whatever key we request, whether or not that key is present */
static NODE* list_find_preceding_from(NODE* prev, size_t key) {
  NODE* curr = prev->next;

  if (!curr || key <= curr->key) return prev;
  return list_find_preceding_from(curr, key);
}


/* Finds a node or the one immediately preceding it if it doesn't exist */
static NODE* list_find_nearest_from(NODE* prev, size_t key) {
  NODE* f;

  if (prev && prev->key == key) return prev;

  f = list_find_preceding_from(prev, key);

  if (!f->next) return f;
  else if (key == f->next->key) return f->next;
  else return prev;
}


/* Finds the node or, if not present, the node that it should follow.
 * NULL indicates no preceding node. */
//static NODE* list_find_nearest(LIST* list, size_t key) {
//  return list_find_nearest_from(list->first, key);
//}


/* Find some element in the list and return the node ptr for that key. */
static NODE* list_find(LIST* list, size_t key) {
  NODE* f;
  if (!list->first) return NULL; // empty list -- does not exist

  // see if we can find it.
  f = list_find_nearest_from(list->first, key);
  if (!f || f->key == key) return f;
  return NULL;
}



/* Get the contents of some set of coordinates. Note: Does not make a copy! Don't free! */
void* list_storage_get(LIST_STORAGE* s, SLICE* slice) {
  //LIST_STORAGE* s = (LIST_STORAGE*)(t);
  size_t r;
  NODE*  n;
  LIST*  l = s->rows;

  for (r = s->rank; r > 1; --r) {
    n = list_find(l, slice->coords[s->rank - r]);
    if (n)  l = n->val;
    else return s->default_val;
  }

  n = list_find(l, slice->coords[s->rank - r]);
  if (n) return n->val;
  else   return s->default_val;
}


/* Returns the value pointer (not the node) for some key. Note that it doesn't free the memory
 * for the value stored in the node -- that pointer gets returned! Only the node is destroyed.
 */
static void* list_remove(LIST* list, size_t key) {
  NODE *f, *rm;
  void* val;

  if (!list->first || list->first->key > key) return NULL; // empty list or def. not present

  if (list->first->key == key) {
    val = list->first->val;
    rm  = list->first;
    list->first = rm->next;
    free(rm);
    return val;
  }

  f = list_find_preceding_from(list->first, key);
  if (!f || !f->next) return NULL; // not found, end of list

  if (f->next->key == key) {
    // remove the node
    rm      = f->next;
    f->next = rm->next;

    // get the value and free the memory for the node
    val = rm->val;
    free(rm);
    return val;
  }

  return NULL; // not found, middle of list
}


/// TODO: Speed up removal.
void* list_storage_remove(LIST_STORAGE* s, SLICE* slice) {
  int r;
  NODE  *n = NULL;
  LIST*  l = s->rows;
  void*  rm = NULL;

  // keep track of where we are in the traversals
  NODE** stack = ALLOCA_N( NODE*, s->rank - 1 );

  for (r = (int)(s->rank); r > 1; --r) {
    n = list_find(l, slice->coords[s->rank - r]); // does this row exist in the matrix?

    if (!n) { // not found
      free(stack);
      return NULL;
    } else { // found
      stack[s->rank - r]    = n;
      l                     = n->val;
    }
  }

  rm = list_remove(l, slice->coords[s->rank - r]);

  // if we removed something, we may now need to remove parent lists
  if (rm) {
    for (r = (int)(s->rank) - 2; r >= 0; --r) { // walk back down the stack
      if (((LIST*)(stack[r]->val))->first == NULL)
        free(list_remove(stack[r]->val, slice->coords[r]));
      else
        break; // no need to continue unless we just deleted one.
    }
  }

  return rm;
}


/* Creates an empty linked list */
static LIST* create_list() {
  LIST* list;
  //if (!(list = malloc(sizeof(LIST)))) return NULL;
  list = ALLOC( LIST );

  //fprintf(stderr, "    create_list LIST: %p\n", list);

  list->first = NULL;
  return list;
}


static NODE* list_insert_after(NODE* node, size_t key, void* val) {
  NODE* ins;

  //if (!(ins = malloc(sizeof(NODE)))) return NULL;
  ins = ALLOC(NODE);

  // insert 'ins' between 'node' and 'node->next'
  ins->next  = node->next;
  node->next = ins;

  // initialize our new node
  ins->key   = key;
  ins->val   = val;

  return ins;
}



/* Given a list and a key/value-ptr pair, create a node (and return that node).
 * If NULL is returned, it means insertion failed.
 * If the key already exists in the list, replace tells it to delete the old value
 * and put in your new one. !replace means delete the new value.
 */
static NODE* list_insert(LIST* list, bool replace, size_t key, void* val) {
  NODE *ins;

  if (list->first == NULL) {                        // List is empty
    //if (!(ins = malloc(sizeof(NODE)))) return NULL;
    ins = ALLOC(NODE);
    ins->next             = NULL;
    ins->val              = val;
    ins->key              = key;
    list->first           = ins;
    return ins;

  } else if (key < list->first->key) {              // Goes at the beginning of the list
    //if (!(ins = malloc(sizeof(NODE)))) return NULL;
    ins = ALLOC(NODE);
    ins->next             = list->first;
    ins->val              = val;
    ins->key              = key;
    list->first           = ins;
    return ins;
  }

  // Goes somewhere else in the list.
  ins = list_find_nearest_from(list->first, key);

  if (ins->key == key) {
    // key already exists
    if (replace) {
      free(ins->val);
      ins->val = val;
    } else free(val);
    return ins;

  } else return list_insert_after(ins, key, val);

}



// TODO: Allow this function to accept an entire row and not just one value -- for slicing
void* list_storage_insert(LIST_STORAGE* s, SLICE* slice, void* val) {
  // Pretend ranks = 2
  // Then coords is going to be size 2
  // So we need to find out if some key already exists
  size_t r;
  NODE*  n;
  LIST*  l = s->rows;

  // drill down into the structure
  for (r = s->rank; r > 1; --r) {
    n = list_insert(l, false, slice->coords[s->rank - r], create_list());
    l = n->val;
  }

  n = list_insert(l, true, slice->coords[s->rank - r], val);
  return n->val;
}

// Creates a list-of-lists(-of-lists-of-lists-etc) storage framework for a matrix.
//
// Note: The pointers you pass in for shape and init_val become property of our new
// storage. You don't need to free them, and you shouldn't re-use them.
LIST_STORAGE* create_list_storage(int8_t dtype, size_t* shape, size_t rank, void* init_val) {
  LIST_STORAGE* s;

  s = ALLOC( LIST_STORAGE );

  s->rank  = rank;
  s->shape = shape;
  s->dtype = dtype;

  s->rows  = create_list();

  s->default_val = init_val;

  return s;
}


static void cast_copy_list_contents(LIST* lhs, LIST* rhs, int8_t lhs_dtype, int8_t rhs_dtype, size_t recursions) {
  NODE *lcurr, *rcurr;

  if (rhs->first) {
    // copy head node
    rcurr = rhs->first;
    lcurr = lhs->first = ALLOC( NODE );

    while (rcurr) {
      lcurr->key = rcurr->key;

      if (recursions == 0) { // contents is some kind of value
        lcurr->val = ALLOC_N(char, nm_sizeof[lhs_dtype]);
        //fprintf(stderr, "    create_val: %p\n", lcurr->val);

        if (lhs_dtype == rhs_dtype) memcpy(lcurr->val, rcurr->val, nm_sizeof[lhs_dtype]);
        else                        SetFuncs[lhs_dtype][rhs_dtype](1, lcurr->val, 0, rcurr->val, 0);

      } else { // contents is a list
        lcurr->val = ALLOC( LIST );
        //fprintf(stderr, "    create_list: %p\n", lcurr->val);

        cast_copy_list_contents(lcurr->val, rcurr->val, lhs_dtype, rhs_dtype, recursions-1);
      }
      if (rcurr->next) lcurr->next = ALLOC( NODE );
      else             lcurr->next = NULL;

      lcurr = lcurr->next;
      rcurr = rcurr->next;
    }
  } else {
    lhs->first = NULL;
  }
}


/* Deletes the linked list and all of its contents. If you want to delete a list inside of a list,
 * set recursions to 1. For lists inside of lists inside of the list, set it to 2; and so on.
 * Setting it to 0 is for no recursions.
 */
static void delete_list(LIST* list, size_t recursions) {
  NODE* next;
  NODE* curr = list->first;

  while (curr != NULL) {
    next = curr->next;

    if (recursions == 0) {
      //fprintf(stderr, "    free_val: %p\n", curr->val);
      free(curr->val);
    } else {
      //fprintf(stderr, "    free_list: %p\n", list);
      delete_list(curr->val, recursions-1);
    }

    free(curr);
    curr = next;
  }
  //fprintf(stderr, "    free_list: %p\n", list);
  free(list);
}


// Copy dense into lists recursively
//
// TODO: This works, but could probably be cleaner (do we really need to pass coords around?)
static bool cast_copy_list_contents_dense(LIST* lhs, const char* rhs, void* zero, int8_t l_dtype, int8_t r_dtype, size_t* pos, size_t* coords, const size_t* shape, size_t rank, size_t recursions) {
  NODE *prev;
  LIST *sub_list;
  bool added = false, added_list = false;
  void* insert_value;

  for (coords[rank-1-recursions] = 0; coords[rank-1-recursions] < shape[rank-1-recursions]; ++coords[rank-1-recursions], ++(*pos)) {
    //fprintf(stderr, "(%u)\t<%u, %u>: ", recursions, coords[0], coords[1]);

    if (recursions == 0) {  // create nodes
      if (!ElemEqEq[r_dtype][0]((char*)rhs + (*pos)*nm_sizeof[r_dtype], zero, 1, nm_sizeof[r_dtype])) { // is not zero
        //fprintf(stderr, "inserting value\n");

        // Create a copy of our value that we will insert in the list
        insert_value = ALLOC_N(char, nm_sizeof[l_dtype]);
        cast_copy_value_single(insert_value, rhs + (*pos)*nm_sizeof[r_dtype], l_dtype, r_dtype);

        if (!lhs->first) prev = list_insert(lhs, false, coords[rank-1-recursions], insert_value);
        else prev = list_insert_after(prev, coords[rank-1-recursions], insert_value);
        added = true;
      } //else fprintf(stderr, "zero\n");
      // no need to do anything if the element is zero
    } else { // create lists
      //fprintf(stderr, "inserting list\n");
      // create a list as if there's something in the row in question, and then delete it if nothing turns out to be there
      sub_list = create_list();

      added_list = cast_copy_list_contents_dense(sub_list, rhs, zero, l_dtype, r_dtype, pos, coords, shape, rank, recursions-1);

      if (!added_list)    { delete_list(sub_list, recursions-1); fprintf(stderr, "deleting list\n"); }// nothing added
      else if (!lhs->first) prev = list_insert(lhs, false, coords[rank-1-recursions], sub_list);
      else                  prev = list_insert_after(prev, coords[rank-1-recursions], sub_list);

      // added = (added || added_list);

    }
  }

  coords[rank-1-recursions] = 0;
  --(*pos);

  return added;
}



LIST_STORAGE* copy_list_storage(LIST_STORAGE* rhs) {
  LIST_STORAGE* lhs;
  size_t* shape;
  void* default_val = ALLOC_N(char, nm_sizeof[rhs->dtype]);

  //fprintf(stderr, "copy_list_storage\n");

  // allocate and copy shape
  shape = ALLOC_N(size_t, rhs->rank);
  memcpy(shape, rhs->shape, rhs->rank * sizeof(size_t));
  memcpy(default_val, rhs->default_val, nm_sizeof[rhs->dtype]);

  lhs = create_list_storage(rhs->dtype, shape, rhs->rank, default_val);

  if (lhs) {
    lhs->rows = create_list();
    cast_copy_list_contents(lhs->rows, rhs->rows, rhs->dtype, rhs->dtype, rhs->rank - 1);
  } else free(shape);

  return lhs;
}


LIST_STORAGE* cast_copy_list_storage(LIST_STORAGE* rhs, int8_t new_dtype) {
  LIST_STORAGE* lhs;
  size_t* shape;
  void* default_val = ALLOC_N(char, nm_sizeof[rhs->dtype]);

  //fprintf(stderr, "copy_list_storage\n");

  // allocate and copy shape
  shape = ALLOC_N(size_t, rhs->rank);
  memcpy(shape, rhs->shape, rhs->rank * sizeof(size_t));

  // copy default value
  if (new_dtype == rhs->dtype)    memcpy(default_val, rhs->default_val, nm_sizeof[rhs->dtype]);
  else                            SetFuncs[new_dtype][rhs->dtype](1, default_val, 0, rhs->default_val, 0);

  lhs = create_list_storage(new_dtype, shape, rhs->rank, default_val);

  lhs->rows = create_list();
  cast_copy_list_contents(lhs->rows, rhs->rows, new_dtype, rhs->dtype, rhs->rank - 1);

  return lhs;
}



LIST_STORAGE* scast_copy_list_dense(const DENSE_STORAGE* rhs, int8_t l_dtype) {
  LIST_STORAGE* lhs;
  size_t pos = 0;
  void* l_default_val = ALLOC_N(char, nm_sizeof[l_dtype]);
  void* r_default_val = ALLOCA_N(char, nm_sizeof[rhs->dtype]); // clean up when finished with this function

  // allocate and copy shape and coords
  size_t *shape = ALLOC_N(size_t, rhs->rank), *coords = ALLOC_N(size_t, rhs->rank);
  memcpy(shape, rhs->shape, rhs->rank * sizeof(size_t));
  memset(coords, 0, rhs->rank * sizeof(size_t));

  // set list default_val to 0
  if (l_dtype == NM_ROBJ) *(VALUE*)l_default_val = INT2FIX(0);
  else                    memset(l_default_val, 0, nm_sizeof[l_dtype]);

  // need test default value for comparing to elements in dense matrix
  if (rhs->dtype == l_dtype)      r_default_val = l_default_val;
  else if (rhs->dtype == NM_ROBJ) *(VALUE*)r_default_val = INT2FIX(0);
  else                            memset(r_default_val, 0, nm_sizeof[rhs->dtype]);

  lhs = create_list_storage(l_dtype, shape, rhs->rank, l_default_val);

  lhs->rows = create_list();
  cast_copy_list_contents_dense(lhs->rows, rhs->elements, r_default_val, l_dtype, rhs->dtype, &pos, coords, rhs->shape, rhs->rank, rhs->rank - 1);

  return lhs;
}


LIST_STORAGE* scast_copy_list_yale(const YALE_STORAGE* rhs, int8_t l_dtype) {
  LIST_STORAGE* lhs;
  NODE *last_added, *last_row_added = NULL;
  LIST* curr_row;
  y_size_t ija, ija_next, i, jj;
  bool add_diag;
  void* default_val = ALLOC_N(char, nm_sizeof[l_dtype]);
  void* R_ZERO = (char*)(rhs->a) + rhs->shape[0]*nm_sizeof[rhs->dtype];
  void* insert_val;

  // allocate and copy shape
  size_t *shape = ALLOC_N(size_t, rhs->rank);
  shape[0] = rhs->shape[0]; shape[1] = rhs->shape[1];

  // copy default value from the zero location in the Yale matrix
  SetFuncs[l_dtype][rhs->dtype](1, default_val, 0, R_ZERO, 0);

  lhs = create_list_storage(l_dtype, shape, rhs->rank, default_val);

  if (rhs->rank != 2)
    rb_raise(nm_eStorageTypeError, "can only convert matrices of rank 2 from yale");

  // Walk through rows and columns as if RHS were a dense matrix
  for (i = 0; i < rhs->shape[0]; ++i) {

    // Get boundaries of beginning and end of row
    YaleGetIJA(ija, rhs, i);
    YaleGetIJA(ija_next, rhs, i+1);

    // Are we going to need to add a diagonal for this row?
    if (ElemEqEq[rhs->dtype][0]((char*)(rhs->a) + i*nm_sizeof[rhs->dtype], R_ZERO, 1, nm_sizeof[rhs->dtype])) add_diag = false; // zero
    else add_diag = true; // nonzero diagonal


    if (ija < ija_next || add_diag) {

      curr_row = create_list();
      last_added = NULL;

      while (ija < ija_next) {
        YaleGetIJA(jj, rhs, ija); // what column number is this?

        // Is there a nonzero diagonal item between the previously added item and the current one?
        if (jj > i && add_diag) {
          // Allocate and copy insertion value
          insert_val = ALLOC_N(char, nm_sizeof[l_dtype]);
          SetFuncs[l_dtype][rhs->dtype](1, insert_val, 0, (char*)(rhs->a) + i*nm_sizeof[rhs->dtype], 0);

          // insert the item in the list at the appropriate location
          if (last_added) last_added = list_insert_after(last_added, i, insert_val);
          else            last_added = list_insert(curr_row, false, i, insert_val);

          add_diag = false; // don't add again!
        }

        // now allocate and add the current item
        insert_val = ALLOC_N(char, nm_sizeof[l_dtype]);
        SetFuncs[l_dtype][rhs->dtype](1, insert_val, 0, (char*)(rhs->a) + ija*nm_sizeof[rhs->dtype], 0);

        if (last_added) last_added = list_insert_after(last_added, jj, insert_val);
        else            last_added = list_insert(curr_row, false, jj, insert_val);

        ++ija; // move to next entry in Yale matrix
      }

      if (add_diag) { // still haven't added the diagonal.
        insert_val = ALLOC_N(char, nm_sizeof[l_dtype]);
        SetFuncs[l_dtype][rhs->dtype](1, insert_val, 0, (char*)(rhs->a) + i*nm_sizeof[rhs->dtype], 0);

        // insert the item in the list at the appropriate location
        if (last_added) last_added = list_insert_after(last_added, i, insert_val);
        else            last_added = list_insert(curr_row, false, i, insert_val);
      }

      // Now add the list at the appropriate location
      if (last_row_added) last_row_added = list_insert_after(last_row_added, i, curr_row);
      else                last_row_added = list_insert(lhs->rows, false, i, curr_row);
    }

  } // end of walk through rows

  return lhs;
}



// Do all values in a list == some value?
static bool list_eqeq_value(const LIST* l, const void* v, int8_t dtype, size_t recursions, size_t* checked) {
  NODE *next, *curr = l->first;

  while (curr) {
    next = curr->next;

    if (recursions == 0) {
      ++(*checked);
      if (!ElemEqEq[dtype][0](curr->val, v, 1, nm_sizeof[dtype])) return false;
    } else if (!list_eqeq_value(curr->val, v, dtype, recursions-1, checked))
      return false;

    curr = next;
  }
  return true;
}


// Are all values in the two lists equal? If one is missing a value, but the other isn't, does the value in the list match
// the default value?
static bool list_eqeq_list(const LIST* left, const LIST* right, const void* left_val, const void* right_val, int8_t dtype, size_t recursions, size_t* checked) {
  NODE *lnext, *lcurr = left->first, *rnext, *rcurr = right->first;

  //fprintf(stderr, "list_eqeq_list: recursions=%d\n", recursions);

  if (lcurr) lnext = lcurr->next;
  if (rcurr) rnext = rcurr->next;

  while (lcurr && rcurr) {

    if (lcurr->key == rcurr->key) {   // MATCHING KEYS
      if (recursions == 0) {
        ++(*checked);
        if (!ElemEqEq[dtype][0](lcurr->val, rcurr->val, 1, nm_sizeof[dtype])) return false;
      } else if (!list_eqeq_list(lcurr->val, rcurr->val, left_val, right_val, dtype, recursions-1, checked))
        return false;

      // increment both iterators
      rcurr = rnext;
      if (rcurr) rnext = rcurr->next;
      lcurr = lnext;
      if (lcurr) lnext = lcurr->next;

    } else if (lcurr->key < rcurr->key) { // NON-MATCHING KEYS

      if (recursions == 0) {
        // compare left entry to right default value
        ++(*checked);
        if (!ElemEqEq[dtype][0](lcurr->val, right_val, 1, nm_sizeof[dtype])) return false;
      } else if (!list_eqeq_value(lcurr->val, right_val, dtype, recursions-1, checked))
        return false;

      // increment left iterator
      lcurr = lnext;
      if (lcurr) lnext = lcurr->next;

    } else { // if (rcurr->key < lcurr->key)

      if (recursions == 0) {
        // compare right entry to left default value
        ++(*checked);
        if (!ElemEqEq[dtype][0](rcurr->val, left_val, 1, nm_sizeof[dtype])) return false;
      } else if (!list_eqeq_value(rcurr->val, left_val, dtype, recursions-1, checked))
        return false;

      // increment right iterator
      rcurr = rnext;
      if (rcurr) rnext = rcurr->next;
    }

  }

  // One final check, in case we get to the end of one list but not the other one.
  if (lcurr) { // nothing left in right-hand list
    if (!ElemEqEq[dtype][0](lcurr->val, right_val, 1, nm_sizeof[dtype])) return false;
  } else if (rcurr) { // nothing left in left-hand list
    if (!ElemEqEq[dtype][0](rcurr->val, left_val, 1, nm_sizeof[dtype])) return false;
  }

  // Nothing different between the two lists -- but make sure after this return that you compare the default values themselves,
  // if we haven't visited every value in the two matrices.
  return true;
}


// Do these two dense matrices of the same dtype have exactly the same contents?
bool list_storage_eqeq(const LIST_STORAGE* left, const LIST_STORAGE* right) {

  // in certain cases, we need to keep track of the number of elements checked.
  size_t num_checked  = 0,
         max_elements = count_storage_max_elements((STORAGE*)left);

  if (!left->rows->first) {
    // fprintf(stderr, "!left->rows true\n");
    // Easy: both lists empty -- just compare default values
    if (!right->rows->first) return ElemEqEq[left->dtype][0](left->default_val, right->default_val, 1, nm_sizeof[left->dtype]);

    // Left empty, right not empty. Do all values in right == left->default_val?
    if (!list_eqeq_value(right->rows, left->default_val, left->dtype, left->rank-1, &num_checked)) return false;

    // If the matrix isn't full, we also need to compare default values.
    if (num_checked < max_elements) return ElemEqEq[left->dtype][0](left->default_val, right->default_val, 1, nm_sizeof[left->dtype]);

  } else if (!right->rows->first) {
    // fprintf(stderr, "!right->rows true\n");
    // Right empty, left not empty. Do all values in left == right->default_val?
    if (!list_eqeq_value(left->rows, right->default_val, left->dtype, left->rank-1, &num_checked)) return false;

    // If the matrix isn't full, we also need to compare default values.
    if (num_checked < max_elements) return ElemEqEq[left->dtype][0](left->default_val, right->default_val, 1, nm_sizeof[left->dtype]);

  } else {
    // fprintf(stderr, "both matrices have entries\n");
    // Hardest case. Compare lists node by node. Let's make it simpler by requiring that both have the same default value
    if (!list_eqeq_list(left->rows, right->rows, left->default_val, right->default_val, left->dtype, left->rank-1, &num_checked)) return false;
    if (num_checked < max_elements) return ElemEqEq[left->dtype][0](left->default_val, right->default_val, 1, nm_sizeof[left->dtype]);
  }

  return true;
}


void delete_list_storage(LIST_STORAGE* s) {
  if (s) {
    //fprintf(stderr, "* Deleting list storage rows at %p\n", s->rows);
    delete_list( s->rows, s->rank - 1 );

    //fprintf(stderr, "  Deleting list storage shape at %p\n", s->shape);
    free(s->shape);
    //fprintf(stderr, "  Deleting list storage default_val at %p\n", s->default_val);
    free(s->default_val);
    //fprintf(stderr, "  Deleting list storage at %p\n", s);
    free(s);
  }
}


static void mark_list(LIST* list, size_t recursions) {
  NODE* next;
  NODE* curr = list->first;

  while (curr != NULL) {
    next = curr->next;
    if (recursions == 0)  rb_gc_mark(*((VALUE*)(curr->val)));
    else                  mark_list(curr->val, recursions-1);
    curr = next;
  }
}


void mark_list_storage(void* m) {
  LIST_STORAGE* storage;

  if (m) {
    storage = (LIST_STORAGE*)(((NMATRIX*)m)->storage);
    if (storage && storage->dtype == NM_ROBJ) {
      rb_gc_mark(*((VALUE*)(storage->default_val)));
      mark_list(storage->rows, storage->rank - 1);
    }
  }
}


#endif
