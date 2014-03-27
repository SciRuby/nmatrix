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
// == sl_list.cpp
//
// Singly-linked list implementation

/*
 * Standard Includes
 */

#include <ruby.h>

/*
 * Project Includes
 */

#include "types.h"

#include "data/data.h"

#include "sl_list.h"

#include "storage/list/list.h"

namespace nm { namespace list {

/*
 * Macros
 */

#ifndef RHASH_SET_IFNONE
#define RHASH_SET_IFNONE(h, v) (RHASH(h)->ifnone = (v))
#endif

/*
 * Global Variables
 */
 

/*
 * Forward Declarations
 */

/*
 * Functions
 */

////////////////
// Lifecycle //
///////////////

/*
 * Creates an empty linked list.
 */
LIST* create(void) {
  LIST* list = NM_ALLOC( LIST );
  list->first = NULL;
  return list;
}

/*
 * Deletes the linked list and all of its contents. If you want to delete a
 * list inside of a list, set recursions to 1. For lists inside of lists inside
 *  of the list, set it to 2; and so on. Setting it to 0 is for no recursions.
 */
void del(LIST* list, size_t recursions) {
  NODE* next;
  NODE* curr = list->first;

  while (curr != NULL) {
    next = curr->next;

    if (recursions == 0) {
      //fprintf(stderr, "    free_val: %p\n", curr->val);
      nm_list_storage_completely_unregister_node(curr);
      NM_FREE(curr->val);
      
    } else {
      //fprintf(stderr, "    free_list: %p\n", list);
      del((LIST*)curr->val, recursions - 1);
    }

    NM_FREE(curr);
    curr = next;
  }
  //fprintf(stderr, "    free_list: %p\n", list);
  NM_FREE(list);
}

/*
 * Documentation goes here.
 */
void mark(LIST* list, size_t recursions) {
  NODE* next;
  NODE* curr = list->first;

  while (curr != NULL) {
    next = curr->next;
    
    if (recursions == 0) {
    	rb_gc_mark(*((VALUE*)(curr->val)));
    	
    } else {
    	mark((LIST*)curr->val, recursions - 1);
    }
    
    curr = next;
  }
}

///////////////
// Accessors //
///////////////


/*
 * Given a list, insert key/val as the first entry in the list. Does not do any
 * checks, just inserts.
 */
NODE* insert_first_node(LIST* list, size_t key, void* val, size_t val_size) {
  NODE* ins   = NM_ALLOC(NODE);
  ins->next   = list->first;

  void* val_copy = NM_ALLOC_N(char, val_size);
  memcpy(val_copy, val, val_size);

  ins->val    = reinterpret_cast<void*>(val_copy);
  ins->key    = key;
  list->first = ins;

  return ins;
}

NODE* insert_first_list(LIST* list, size_t key, LIST* l) {
  NODE* ins   = NM_ALLOC(NODE);
  ins->next   = list->first;

  ins->val    = reinterpret_cast<void*>(l);
  ins->key    = key;
  list->first = ins;

  return ins;
}


/* 
 * Given a list and a key/value-ptr pair, create a node (and return that node).
 * If NULL is returned, it means insertion failed.
 * If the key already exists in the list, replace tells it to delete the old
 * value and put in your new one. !replace means delete the new value.
 */
NODE* insert(LIST* list, bool replace, size_t key, void* val) {
  NODE *ins;

  if (list->first == NULL) {
  	// List is empty
  	
    //if (!(ins = malloc(sizeof(NODE)))) return NULL;
    ins = NM_ALLOC(NODE);
    ins->next             = NULL;
    ins->val              = val;
    ins->key              = key;
    list->first           = ins;
    
    return ins;

  } else if (key < list->first->key) {
  	// Goes at the beginning of the list
  	
    //if (!(ins = malloc(sizeof(NODE)))) return NULL;
    ins = NM_ALLOC(NODE);
    ins->next             = list->first;
    ins->val              = val;
    ins->key              = key;
    list->first           = ins;
    
    return ins;
  }

  // Goes somewhere else in the list.
  ins = find_nearest_from(list->first, key);

  if (ins->key == key) {
    // key already exists
    if (replace) {
      nm_list_storage_completely_unregister_node(ins);
      NM_FREE(ins->val);
      ins->val = val;
    } else {
      NM_FREE(val);
    }
    
    return ins;

  } else {
  	return insert_after(ins, key, val);
  }
}



/*
 * Documentation goes here.
 */
NODE* insert_after(NODE* node, size_t key, void* val) {
  //if (!(ins = malloc(sizeof(NODE)))) return NULL;
  NODE* ins = NM_ALLOC(NODE);

  // insert 'ins' between 'node' and 'node->next'
  ins->next  = node->next;
  node->next = ins;

  // initialize our new node
  ins->key   = key;
  ins->val   = val;

  return ins;
}


/*
 * Insert a new node immediately after +node+, or replace the existing one if its key is a match.
 */
NODE* replace_insert_after(NODE* node, size_t key, void* val, bool copy, size_t copy_size) {
  if (node->next && node->next->key == key) {

    // Should we copy into the current one or free and insert?
    if (copy) memcpy(node->next->val, val, copy_size);
    else {
      NM_FREE(node->next->val);
      node->next->val = val;
    }

    return node->next;

  } else { // no next node, or if there is one, it's greater than the current key

    if (copy) {
      void* val_copy = NM_ALLOC_N(char, copy_size);
      memcpy(val_copy, val, copy_size);
      return insert_after(node, key, val_copy);
    } else {
      return insert_after(node, key, val);
    }

  }
}



/*
 * Functions analogously to list::insert but this inserts a copy of the value instead of the original.
 */
NODE* insert_copy(LIST *list, bool replace, size_t key, void *val, size_t size) {
  void *copy_val = NM_ALLOC_N(char, size);
  memcpy(copy_val, val, size);

  return insert(list, replace, key, copy_val);
}


/*
 * Returns the value pointer for some key. Doesn't free the memory for that value. Doesn't require a find operation,
 * assumes finding has already been done. If rm is the first item in the list, prev should be NULL.
 */
void* remove_by_node(LIST* list, NODE* prev, NODE* rm) {
  if (!prev)  list->first = rm->next;
  else        prev->next  = rm->next;

  void* val   = rm->val;
  NM_FREE(rm);

  return val;
}


/*
 * Returns the value pointer (not the node) for some key. Note that it doesn't
 * free the memory for the value stored in the node -- that pointer gets
 * returned! Only the node is destroyed.
 */
void* remove_by_key(LIST* list, size_t key) {
  NODE *f, *rm;
  void* val;

  if (!list->first || list->first->key > key) { // empty list or def. not present
  	return NULL;
  }

  if (list->first->key == key) {
    val = list->first->val;
    rm  = list->first;
    
    list->first = rm->next;
    NM_FREE(rm);
    
    return val;
  }

  f = find_preceding_from_node(list->first, key);
  if (!f || !f->next) { // not found, end of list
  	return NULL;
  }

  if (f->next->key == key) {
    // remove the node
    rm      = f->next;
    f->next = rm->next;

    // get the value and free the memory for the node
    val = rm->val;
    NM_FREE(rm);

    return val;
  }

  return NULL; // not found, middle of list
}


bool node_is_within_slice(NODE* n, size_t coord, size_t len) {
  if (!n) return false;
  if (n->key >= coord && n->key < coord + len) return true;
  else return false;
}


/*
 * Recursive removal of lists that may contain sub-lists. Stores the value ultimately removed in rm.
 */
bool remove_recursive(LIST* list, const size_t* coords, const size_t* offsets, const size_t* lengths, size_t r, const size_t& dim) {
//  std::cerr << "remove_recursive: " << r << std::endl;
  // find the current coordinates in the list
  NODE* prev    = find_preceding_from_list(list, coords[r] + offsets[r]);
  NODE* n;
  if (prev) n  = prev->next && node_is_within_slice(prev->next, coords[r] + offsets[r], lengths[r]) ? prev->next : NULL;
  else      n  = node_is_within_slice(list->first, coords[r] + offsets[r], lengths[r]) ? list->first : NULL;

  if (r < dim-1) { // nodes here are lists

    while (n) {
      // from that sub-list, call remove_recursive.
      bool remove_parent = remove_recursive(reinterpret_cast<LIST*>(n->val), coords, offsets, lengths, r+1, dim);

      if (remove_parent) { // now empty -- so remove the sub-list
//        std::cerr << r << ": removing parent list at " << n->key << std::endl;
        NM_FREE(remove_by_node(list, prev, n));

        if (prev) n  = prev->next && node_is_within_slice(prev->next, coords[r] + offsets[r], lengths[r]) ? prev->next : NULL;
        else      n  = node_is_within_slice(list->first, coords[r] + offsets[r], lengths[r]) ? list->first : NULL;
      } else {
        // Move forward to next node (list at n still exists)
        prev         = n;
        n            = prev->next && node_is_within_slice(prev->next, coords[r] + offsets[r], lengths[r]) ? prev->next : NULL;
      }

      // Iterate to next one.
      if (prev) n  = prev->next && node_is_within_slice(prev->next, coords[r] + offsets[r], lengths[r]) ? prev->next : NULL;
      else      n  = node_is_within_slice(list->first, coords[r] + offsets[r], lengths[r]) ? list->first : NULL;
    }

  } else { // nodes here are not lists, but actual values

    while (n) {
//      std::cerr << r << ": removing node at " << n->key << std::endl;
      NM_FREE(remove_by_node(list, prev, n));

      if (prev) n  = prev->next && node_is_within_slice(prev->next, coords[r] + offsets[r], lengths[r]) ? prev->next : NULL;
      else      n  = node_is_within_slice(list->first, coords[r] + offsets[r], lengths[r]) ? list->first : NULL;
    }
  }

  if (!list->first) return true; // if current list is now empty, signal its removal

  return false;
}

///////////
// Tests //
///////////


/////////////
// Utility //
/////////////

/*
 * Find some element in the list and return the node ptr for that key.
 */
NODE* find(LIST* list, size_t key) {
  NODE* f;
  if (!list->first) {
  	// empty list -- does not exist
  	return NULL;
  }

  // see if we can find it.
  f = find_nearest_from(list->first, key);
  
  if (!f || f->key == key) {
  	return f;
  }
  
  return NULL;
}



/*
 * Find some element in the list and return the node ptr for that key.
 */
NODE* find_with_preceding(LIST* list, size_t key, NODE*& prev) {
  if (!prev) prev = list->first;
  if (!prev) return NULL; // empty list, does not exist

  if (prev->key == key) {
    NODE* n = prev;
    prev    = NULL;
    return n;
  }

  while (prev->next && prev->next->key < key) {
    prev = prev->next;
  }

  return prev->next;
}




/*
 * Finds the node that should go before whatever key we request, whether or not
 * that key is present.
 */
NODE* find_preceding_from_node(NODE* prev, size_t key) {
  NODE* curr = prev->next;

  if (!curr || key <= curr->key) {
  	return prev;
  	
  } else {
  	return find_preceding_from_node(curr, key);
  }
}


/*
 * Returns NULL if the key being sought is first in the list or *should* be first in the list but is absent. Otherwise
 * returns the previous node to where that key is or should be.
 */
NODE* find_preceding_from_list(LIST* l, size_t key) {
  NODE* n = l->first;
  if (!n || n->key >= key)  return NULL;
  else                      return find_preceding_from_node(n, key);
}

/*
 * Finds the node or, if not present, the node that it should follow. NULL
 * indicates no preceding node.
 */
NODE* find_nearest(LIST* list, size_t key) {
  return find_nearest_from(list->first, key);
}

/*
 * Finds a node or the one immediately preceding it if it doesn't exist.
 */
NODE* find_nearest_from(NODE* prev, size_t key) {
  NODE* f;

  if (prev && prev->key == key) {
  	return prev;
  }

  f = find_preceding_from_node(prev, key);

  if (!f->next) { // key exceeds final node; return final node.
  	return f;
  	
  } else if (key == f->next->key) { // node already present; return location
  	return f->next;

  } else {
  	return f;
  }
}

/////////////////////////
// Copying and Casting //
/////////////////////////


/*
 * Copy the contents of a list.
 */
template <typename LDType, typename RDType>
void cast_copy_contents(LIST* lhs, const LIST* rhs, size_t recursions) {
  NODE *lcurr, *rcurr;

  if (rhs->first) {
    // copy head node
    rcurr = rhs->first;
    lcurr = lhs->first = NM_ALLOC( NODE );

    while (rcurr) {
      lcurr->key = rcurr->key;

      if (recursions == 0) {
      	// contents is some kind of value

        lcurr->val = NM_ALLOC( LDType );

        *reinterpret_cast<LDType*>(lcurr->val) = *reinterpret_cast<RDType*>( rcurr->val );

      } else {
      	// contents is a list

        lcurr->val = NM_ALLOC( LIST );

        cast_copy_contents<LDType, RDType>(
          reinterpret_cast<LIST*>(lcurr->val),
          reinterpret_cast<LIST*>(rcurr->val),
          recursions-1
        );
      }

      if (rcurr->next) {
      	lcurr->next = NM_ALLOC( NODE );

      } else {
      	lcurr->next = NULL;
      }

      lcurr = lcurr->next;
      rcurr = rcurr->next;
    }

  } else {
    lhs->first = NULL;
  }
}

}} // end of namespace nm::list

extern "C" {

  /*
   * C access for copying the contents of a list.
   */
  void nm_list_cast_copy_contents(LIST* lhs, const LIST* rhs, nm::dtype_t lhs_dtype, nm::dtype_t rhs_dtype, size_t recursions) {
    LR_DTYPE_TEMPLATE_TABLE(nm::list::cast_copy_contents, void, LIST*, const LIST*, size_t);

    ttable[lhs_dtype][rhs_dtype](lhs, rhs, recursions);
  }

  /*
   * Sets up a hash with an appropriate default values. That means that if recursions == 0, the default value is default_value,
   * but if recursions == 1, the default value is going to be a hash with default value of default_value, and if recursions == 2,
   * the default value is going to be a hash with default value of hash with default value of default_value, and so on.
   * In other words, it's recursive.
   */
  static VALUE empty_list_to_hash(const nm::dtype_t dtype, size_t recursions, VALUE default_value) {
    VALUE h = rb_hash_new();
    if (recursions) {
      RHASH_SET_IFNONE(h, empty_list_to_hash(dtype, recursions-1, default_value));
    } else {
      RHASH_SET_IFNONE(h, default_value);
    }
    return h;
  }


  /*
   * Copy a list to a Ruby Hash
   */
  VALUE nm_list_copy_to_hash(const LIST* l, const nm::dtype_t dtype, size_t recursions, VALUE default_value) {

    // Create a hash with default values appropriately specified for a sparse matrix.
    VALUE h = empty_list_to_hash(dtype, recursions, default_value);

    if (l->first) {
      NODE* curr = l->first;

      while (curr) {

        size_t key = curr->key;

        if (recursions == 0) { // content is some kind of value
          rb_hash_aset(h, INT2FIX(key), rubyobj_from_cval(curr->val, dtype).rval);
        } else { // content is a list
          rb_hash_aset(h, INT2FIX(key), nm_list_copy_to_hash(reinterpret_cast<const LIST*>(curr->val), dtype, recursions-1, default_value));
        }

        curr = curr->next;

      }

    }

    return h;
  }


} // end of extern "C" block

