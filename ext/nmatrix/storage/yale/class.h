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
// == class.h
//
// Object-oriented interface for Yale.
//

#ifndef YALE_CLASS_H
# define YALE_CLASS_H

#include "../dense/dense.h"
#include "math/transpose.h"
#include "yale.h"

namespace nm {


/*
 * This class is basically an intermediary for YALE_STORAGE objects which enables us to treat it like a C++ object. It
 * keeps the src pointer as its s, along with other relevant slice information.
 *
 * It's useful for creating iterators and such. It isn't responsible for allocating or freeing its YALE_STORAGE* pointers.
 */

template <typename D>
class YaleStorage {
public:
  YaleStorage(const YALE_STORAGE* storage)
   : s(reinterpret_cast<YALE_STORAGE*>(storage->src)),
     slice(storage != storage->src),
     slice_shape(storage->shape),
     slice_offset(storage->offset)
  {
    nm_yale_storage_register(storage->src);
  }

  YaleStorage(const STORAGE* storage)
   : s(reinterpret_cast<YALE_STORAGE*>(storage->src)),
     slice(storage != storage->src),
     slice_shape(storage->shape),
     slice_offset(storage->offset)
  {
    nm_yale_storage_register(reinterpret_cast<STORAGE*>(storage->src));
  }

  ~YaleStorage() {
    nm_yale_storage_unregister(s);
  }

  /* Allows us to do YaleStorage<uint8>::dtype() to get an nm::dtype_t */
  static nm::dtype_t dtype() {
    return nm::ctype_to_dtype_enum<D>::value_type;
  }


  bool is_ref() const { return slice; }

  inline D* default_obj_ptr() { return &(a(s->shape[0])); }
  inline D& default_obj() { return a(s->shape[0]); }
  inline const D& default_obj() const { return a(s->shape[0]); }
  inline const D& const_default_obj() const { return a(s->shape[0]); }


  /*
   * Return a Ruby VALUE representation of default_obj()
   */
  VALUE const_default_value() const {
    return nm::yale_storage::nm_rb_dereference(a(s->shape[0]));
  }

  inline size_t* ija_p()       const       { return reinterpret_cast<size_t*>(s->ija); }
  inline const size_t& ija(size_t p) const { return ija_p()[p]; }
  inline size_t& ija(size_t p)             { return ija_p()[p]; }
  inline D* a_p()         const       { return reinterpret_cast<D*>(s->a); }
  inline const D& a(size_t p) const   { return a_p()[p]; }
  inline D& a(size_t p)               { return a_p()[p]; }

  bool real_row_empty(size_t i) const { return ija(i+1) - ija(i) == 0 ? true : false; }

  inline size_t* shape_p()        const { return slice_shape;      }
  inline size_t  shape(uint8_t d) const { return slice_shape[d];   }
  inline size_t* real_shape_p() const { return s->shape;           }
  inline size_t  real_shape(uint8_t d) const { return s->shape[d]; }
  inline size_t* offset_p()     const { return slice_offset;       }
  inline size_t  offset(uint8_t d) const { return slice_offset[d]; }
  inline size_t  capacity() const { return s->capacity;            }
  inline size_t  size() const { return ija(real_shape(0));         }


  /*
   * Returns true if the value at apos is the default value.
   * Mainly used for determining if the diagonal contains zeros.
   */
  bool is_pos_default_value(size_t apos) const {
    return (a(apos) == const_default_obj());
  }

  /*
   * Given a size-2 array of size_t, representing the shape, determine
   * the maximum size of YaleStorage arrays.
   */
  static size_t max_size(const size_t* shape) {
    size_t result = shape[0] * shape[1] + 1;
    if (shape[0] > shape[1])
      result += shape[0] - shape[1];
    return result;
  }

  /*
   * Minimum size of Yale Storage arrays given some shape.
   */
  static size_t min_size(const size_t* shape) {
    return shape[0]*2 + 1;
  }

  /*
   * This is the guaranteed maximum size of the IJA/A arrays of the matrix given its shape.
   */
  inline size_t real_max_size() const {
    return YaleStorage<D>::max_size(real_shape_p());
  }

  // Binary search between left and right in IJA for column ID real_j. Returns left if not found.
  size_t real_find_pos(size_t left, size_t right, size_t real_j, bool& found) const {
    if (left > right) {
      found = false;
      return left;
    }

    size_t mid   = (left + right) / 2;
    size_t mid_j = ija(mid);

    if (mid_j == real_j) {
      found = true;
      return mid;
    } else if (mid_j > real_j)  return real_find_pos(left, mid - 1, real_j, found);
    else                        return real_find_pos(mid + 1, right, real_j, found);
  }

  // Binary search between left and right in IJA for column ID real_j. Essentially finds where the slice should begin,
  // with no guarantee that there's anything in there.
  size_t real_find_left_boundary_pos(size_t left, size_t right, size_t real_j) const {
    if (left > right) return right;
    if (ija(left) >= real_j) return left;

    size_t mid   = (left + right) / 2;
    size_t mid_j = ija(mid);

    if (mid_j == real_j)      return mid;
    else if (mid_j > real_j)  return real_find_left_boundary_pos(left, mid, real_j);
    else                      return real_find_left_boundary_pos(mid + 1, right, real_j);
  }

  // Binary search between left and right in IJA for column ID real_j. Essentially finds where the slice should begin,
  // with no guarantee that there's anything in there.
  size_t real_find_right_boundary_pos(size_t left, size_t right, size_t real_j) const {
    if (left > right) return right;
    if (ija(right) <= real_j) return right;

    size_t mid   = (left + right) / 2;
    size_t mid_j = ija(mid);

    if (mid_j == real_j)      return mid;
    else if (mid_j > real_j)  return real_find_right_boundary_pos(left, mid, real_j);
    else                      return real_find_right_boundary_pos(mid + 1, right, real_j);
  }


  // Binary search for coordinates i,j in the slice. If not found, return -1.
  std::pair<size_t,bool> find_pos(const std::pair<size_t,size_t>& ij) const {
    size_t left   = ija(ij.first + offset(0));
    size_t right  = ija(ij.first + offset(0) + 1) - 1;

    std::pair<size_t, bool> result;
    result.first = real_find_pos(left, right, ij.second + offset(1), result.second);
    return result;
  }

  // Binary search for coordinates i,j in the slice, and return the first position >= j in row i.
  size_t find_pos_for_insertion(size_t i, size_t j) const {
    size_t left   = ija(i + offset(0));
    size_t right  = ija(i + offset(0) + 1) - 1;

    // Check that the right search point is valid. rflbp will check to make sure the left is valid relative to left.
    if (right > ija(real_shape(0))) {
      right = ija(real_shape(0))-1;
    }
    size_t result = real_find_left_boundary_pos(left, right, j + offset(1));
    return result;
  }

  typedef yale_storage::basic_iterator_T<D,D,YaleStorage<D> >              basic_iterator;
  typedef yale_storage::basic_iterator_T<D,const D,const YaleStorage<D> >  const_basic_iterator;

  typedef yale_storage::stored_diagonal_iterator_T<D,D,YaleStorage<D> >              stored_diagonal_iterator;
  typedef yale_storage::stored_diagonal_iterator_T<D,const D,const YaleStorage<D> >  const_stored_diagonal_iterator;

  typedef yale_storage::iterator_T<D,D,YaleStorage<D> >                iterator;
  typedef yale_storage::iterator_T<D,const D,const YaleStorage<D> >    const_iterator;


  friend class yale_storage::row_iterator_T<D,D,YaleStorage<D> >;
  typedef yale_storage::row_iterator_T<D,D,YaleStorage<D> >             row_iterator;
  typedef yale_storage::row_iterator_T<D,const D,const YaleStorage<D> > const_row_iterator;

  typedef yale_storage::row_stored_iterator_T<D,D,YaleStorage<D>,row_iterator>    row_stored_iterator;
  typedef yale_storage::row_stored_nd_iterator_T<D,D,YaleStorage<D>,row_iterator> row_stored_nd_iterator;
  typedef yale_storage::row_stored_iterator_T<D,const D,const YaleStorage<D>,const_row_iterator>       const_row_stored_iterator;
  typedef yale_storage::row_stored_nd_iterator_T<D,const D,const YaleStorage<D>,const_row_iterator>    const_row_stored_nd_iterator;
  typedef std::pair<row_iterator,row_stored_nd_iterator>                                               row_nd_iter_pair;

  // Variety of iterator begin and end functions.
  iterator begin(size_t row = 0)                      {      return iterator(*this, row);                 }
  iterator row_end(size_t row)                        {      return begin(row+1);                         }
  iterator end()                                      {      return iterator(*this, shape(0));            }
  const_iterator cbegin(size_t row = 0) const         {      return const_iterator(*this, row);           }
  const_iterator crow_end(size_t row) const           {      return cbegin(row+1);                        }
  const_iterator cend() const                         {      return const_iterator(*this, shape(0));      }

  stored_diagonal_iterator sdbegin(size_t d = 0)      {      return stored_diagonal_iterator(*this, d);   }
  stored_diagonal_iterator sdend()                    {
    return stored_diagonal_iterator(*this, std::min( shape(0) + offset(0), shape(1) + offset(1) ) - std::max(offset(0), offset(1)) );
  }
  const_stored_diagonal_iterator csdbegin(size_t d = 0) const { return const_stored_diagonal_iterator(*this, d); }
  const_stored_diagonal_iterator csdend() const        {
    return const_stored_diagonal_iterator(*this, std::min( shape(0) + offset(0), shape(1) + offset(1) ) - std::max(offset(0), offset(1)) );
  }
  row_iterator ribegin(size_t row = 0)                {      return row_iterator(*this, row);             }
  row_iterator riend()                                {      return row_iterator(*this, shape(0));        }
  const_row_iterator cribegin(size_t row = 0) const   {      return const_row_iterator(*this, row);       }
  const_row_iterator criend() const                   {      return const_row_iterator(*this, shape(0));  }


  /*
   * Get a count of the ndnz in the slice as if it were its own matrix.
   */
  size_t count_copy_ndnz() const {
    if (!slice) return s->ndnz; // easy way -- not a slice.
    size_t count = 0;

    // Visit all stored entries.
    for (const_row_iterator it = cribegin(); it != criend(); ++it){
      for (auto jt = it.begin(); jt != it.end(); ++jt) {
        if (it.i() != jt.j() && *jt != const_default_obj()) ++count;
      }
    }

    return count;
  }

  /*
   * Returns the iterator for i,j or snd_end() if not found.
   */
/*  stored_nondiagonal_iterator find(const std::pair<size_t,size_t>& ij) {
    std::pair<size_t,bool> find_pos_result = find_pos(ij);
    if (!find_pos_result.second) return sndend();
    else return stored_nondiagonal_iterator(*this, ij.first, find_pos_result.first);
  } */

  /*
   * Returns a stored_nondiagonal_iterator pointing to the location where some coords i,j should go, or returns their
   * location if present.
   */
  /*std::pair<row_iterator, row_stored_nd_iterator> lower_bound(const std::pair<size_t,size_t>& ij)  {
    row_iterator it            = ribegin(ij.first);
    row_stored_nd_iterator jt  = it.lower_bound(ij.second);
    return std::make_pair(it,jt);
  } */

  class multi_row_insertion_plan {
  public:
    std::vector<size_t>   pos;
    std::vector<int>      change;
    int                   total_change; // the net change occurring
    size_t                num_changes;  // the total number of rows that need to change size
    multi_row_insertion_plan(size_t rows_in_slice) : pos(rows_in_slice), change(rows_in_slice), total_change(0), num_changes(0) { }

    void add(size_t i, const std::pair<int,size_t>& change_and_pos) {
      pos[i]        = change_and_pos.second;
      change[i]     = change_and_pos.first;
      total_change += change_and_pos.first;
      if (change_and_pos.first != 0) num_changes++;
    }
  };


  /*
   * Find all the information we need in order to modify multiple rows.
   */
  multi_row_insertion_plan insertion_plan(row_iterator i, size_t j, size_t* lengths, D* const v, size_t v_size) const {
    multi_row_insertion_plan p(lengths[0]);

    // v_offset is our offset in the array v. If the user wants to change two elements in each of three rows,
    // but passes an array of size 3, we need to know that the second insertion plan must start at position
    // 2 instead of 0; and then the third must start at 1.
    size_t v_offset = 0;
    for (size_t m = 0; m < lengths[0]; ++m, ++i) {
      p.add(m, i.single_row_insertion_plan(j, lengths[1], v, v_size, v_offset));
    }

    return p;
  }



  /*
   * Insert entries in multiple rows. Slice-setting.
   */
  void insert(row_iterator i, size_t j, size_t* lengths, D* const v, size_t v_size) {
    // Expensive pre-processing step: find all the information we need in order to do insertions.
    multi_row_insertion_plan p = insertion_plan(i, j, lengths, v, v_size);

    // There are more efficient ways to do this, but this is the low hanging fruit version of the algorithm.
    // Here's the full problem: http://stackoverflow.com/questions/18753375/algorithm-for-merging-short-lists-into-a-long-vector
    // --JW

    bool resize = false;
    size_t sz = size();
    if (p.num_changes > 1) resize = true; // TODO: There are surely better ways to do this, but I've gone for the low-hanging fruit
    else if (sz + p.total_change > capacity() || sz + p.total_change <= capacity() / nm::yale_storage::GROWTH_CONSTANT) resize = true;

    if (resize) {
      update_resize_move_insert(i.i() + offset(0), j + offset(1), lengths, v, v_size, p);
    } else {

      // Make the necessary modifications, which hopefully can be done in-place.
      size_t v_offset = 0;
      //int accum       = 0;
      for (size_t ii = 0; ii < lengths[0]; ++ii, ++i) {
        i.insert(row_stored_nd_iterator(i, p.pos[ii]), j, lengths[1], v, v_size, v_offset);
      }
    }
  }


  /*
   * Most Ruby-centric insert function. Accepts coordinate information in slice,
   * and value information of various types in +right+. This function must evaluate
   * +right+ and determine what other functions to call in order to properly handle
   * it.
   */
  void insert(SLICE* slice, VALUE right) {

    NM_CONSERVATIVE(nm_register_value(right));

    std::pair<NMATRIX*,bool> nm_and_free =
      interpret_arg_as_dense_nmatrix(right, dtype());
    // Map the data onto D* v

    D*     v;
    size_t v_size = 1;

    if (nm_and_free.first) {
      DENSE_STORAGE* s = reinterpret_cast<DENSE_STORAGE*>(nm_and_free.first->storage);
      v       = reinterpret_cast<D*>(s->elements);
      v_size  = nm_storage_count_max_elements(s);

    } else if (TYPE(right) == T_ARRAY) {
      v_size = RARRAY_LEN(right);
      v      = NM_ALLOC_N(D, v_size);
      if (dtype() == nm::RUBYOBJ) {
       nm_register_values(reinterpret_cast<VALUE*>(v), v_size);
      }
      for (size_t m = 0; m < v_size; ++m) {
        rubyval_to_cval(rb_ary_entry(right, m), s->dtype, &(v[m]));
      }
      if (dtype() == nm::RUBYOBJ) {
       nm_unregister_values(reinterpret_cast<VALUE*>(v), v_size);
      }

    } else {
      v = reinterpret_cast<D*>(rubyobj_to_cval(right, dtype()));
    }

    row_iterator i = ribegin(slice->coords[0]);

    if (slice->single || (slice->lengths[0] == 1 && slice->lengths[1] == 1)) { // single entry
      i.insert(slice->coords[1], *v);
    } else if (slice->lengths[0] == 1) { // single row, multiple entries
      i.insert(slice->coords[1], slice->lengths[1], v, v_size);
    } else { // multiple rows, unknown number of entries
      insert(i, slice->coords[1], slice->lengths, v, v_size);
    }

    // Only free v if it was allocated in this function.
    if (nm_and_free.first) {
      if (nm_and_free.second) {
        nm_delete(nm_and_free.first);
      }
    } else NM_FREE(v);

    NM_CONSERVATIVE(nm_unregister_value(right));
  }


  /*
   * Remove an entry from an already found non-diagonal position.
   */
  row_iterator erase(row_iterator it, const row_stored_nd_iterator& position) {
    it.erase(position);
    return it;
  }


  /*
   * Remove an entry from the matrix at the already-located position. If diagonal, just sets to default; otherwise,
   * actually removes the entry.
   */
  row_iterator erase(row_iterator it, const row_stored_iterator& jt) {
    it.erase((const row_stored_nd_iterator&)jt);
    return it;
  }


  row_iterator insert(row_iterator it, row_stored_iterator position, size_t j, const D& val) {
    it.insert(position, j, val);
    return it;
  }


  /*
   * Insert an element in column j, using position's p() as the location to insert the new column. i and j will be the
   * coordinates. This also does a replace if column j is already present.
   *
   * Returns true if a new entry was added and false if an entry was replaced.
   *
   * Pre-conditions:
   *   - position.p() must be between ija(real_i) and ija(real_i+1), inclusive, where real_i = i + offset(0)
   *   - real_i and real_j must not be equal
   */
  row_iterator insert(row_iterator it, row_stored_nd_iterator position, size_t j, const D& val) {
    it.insert(position, j, val);
    return it;
  }


  /*
   * Insert n elements v in columns j, using position as a guide. i gives the starting row. If at any time a value in j
   * decreases,
   */
  /*bool insert(stored_iterator position, size_t n, size_t i, size_t* j, DType* v) {

  } */

  /*
   * A pseudo-insert operation, since the diagonal portion of the A array is constant size.
   */
  stored_diagonal_iterator insert(stored_diagonal_iterator position, const D& val) {
    *position = val;
    return position;
  }


/*  iterator insert(iterator position, size_t j, const D& val) {
    if (position.real_i() == position.real_j()) {
      s->a(position.real_i()) = val;
      return position;
    } else {
      row_iterator it = ribegin(position.i());
      row_stored_nd_iterator position = it.ndbegin(j);
      return insert(it, position, j, val);
    }
  }*/




  /*
   * Returns a pointer to the location of some entry in the matrix.
   *
   * This is needed for backwards compatibility. We don't really want anyone
   * to modify the contents of that pointer, because it might be the ZERO location.
   *
   * TODO: Change all storage_get functions to return a VALUE once we've put list and
   * dense in OO mode. ???
   */
  inline D* get_single_p(SLICE* slice) {
    size_t real_i = offset(0) + slice->coords[0],
           real_j = offset(1) + slice->coords[1];

    if (real_i == real_j)
      return &(a(real_i));

    if (ija(real_i) == ija(real_i+1))
      return default_obj_ptr(); // zero pointer

    // binary search for a column's location
    std::pair<size_t,bool> p = find_pos(std::make_pair(slice->coords[0], slice->coords[1]));
    if (p.second)
      return &(a(p.first));
                       // not found: return default
    return default_obj_ptr(); // zero pointer
  }


  /*
   * Allocate a reference pointing to s. Note that even if +this+ is a reference,
   * we can create a reference within it.
   *
   * Note: Make sure you NM_FREE() the result of this call. You can't just cast it
   * directly into a YaleStorage<D> class.
   */
  YALE_STORAGE* alloc_ref(SLICE* slice) {
    YALE_STORAGE* ns  = NM_ALLOC( YALE_STORAGE );

    ns->dim           = s->dim;
    ns->offset        = NM_ALLOC_N(size_t, ns->dim);
    ns->shape         = NM_ALLOC_N(size_t, ns->dim);

    for (size_t d = 0; d < ns->dim; ++d) {
      ns->offset[d]   = slice->coords[d]  + offset(d);
      ns->shape[d]    = slice->lengths[d];
    }

    ns->dtype         = s->dtype;
    ns->a             = a_p();
    ns->ija           = ija_p();

    ns->src           = s;
    s->count++;

    ns->ndnz          = 0;
    ns->capacity      = 0;

    return ns;
  }


  /*
   * Allocates and initializes the basic struct (but not IJA or A vectors).
   */
  static YALE_STORAGE* alloc(size_t* shape, size_t dim = 2) {
    YALE_STORAGE* s = NM_ALLOC( YALE_STORAGE );

    s->ndnz         = 0;
    s->dtype        = dtype();
    s->shape        = shape;
    s->offset       = NM_ALLOC_N(size_t, dim);
    for (size_t d = 0; d < dim; ++d)
      s->offset[d]  = 0;
    s->dim          = dim;
    s->src          = reinterpret_cast<STORAGE*>(s);
    s->count        = 1;

    return s;
  }


  /*
   * Create basic storage of same dtype as YaleStorage<D>. Allocates it,
   * reserves necessary space, but doesn't fill structure at all.
   */
  static YALE_STORAGE* create(size_t* shape, size_t reserve) {

    YALE_STORAGE* s = alloc( shape, 2 );
    size_t max_sz   = YaleStorage<D>::max_size(shape),
           min_sz   = YaleStorage<D>::min_size(shape);

    if (reserve < min_sz) {
      s->capacity = min_sz;
    } else if (reserve > max_sz) {
      s->capacity = max_sz;
    } else {
      s->capacity = reserve;
    }

    s->ija = NM_ALLOC_N( size_t, s->capacity );
    s->a   = NM_ALLOC_N( D,      s->capacity );

    return s;
  }


  /*
   * Clear out the D portion of the A vector (clearing the diagonal and setting
   * the zero value).
   */
  static void clear_diagonal_and_zero(YALE_STORAGE& s, D* init_val = NULL) {
    D* a  = reinterpret_cast<D*>(s.a);

    // Clear out the diagonal + one extra entry
    if (init_val) {
      for (size_t i = 0; i <= s.shape[0]; ++i)
        a[i] = *init_val;
    } else {
      for (size_t i = 0; i <= s.shape[0]; ++i)
        a[i] = 0;
    }
  }


  /*
   * Empty the matrix by initializing the IJA vector and setting the diagonal to 0.
   *
   * Called when most YALE_STORAGE objects are created.
   *
   * Can't go inside of class YaleStorage because YaleStorage creation requires that
   * IJA already be initialized.
   */
  static void init(YALE_STORAGE& s, D* init_val) {
    size_t IA_INIT = s.shape[0] + 1;
    for (size_t m = 0; m < IA_INIT; ++m) {
      s.ija[m] = IA_INIT;
    }

    clear_diagonal_and_zero(s, init_val);
  }


  /*
   * Make a very basic allocation. No structure or copy or anything. It'll be shaped like this
   * matrix.
   *
   * TODO: Combine this with ::create()'s ::alloc(). These are redundant.
   */
   template <typename E>
   YALE_STORAGE* alloc_basic_copy(size_t new_capacity, size_t new_ndnz) const {
     nm::dtype_t new_dtype = nm::ctype_to_dtype_enum<E>::value_type;
     YALE_STORAGE* lhs     = NM_ALLOC( YALE_STORAGE );
     lhs->dim              = s->dim;
     lhs->shape            = NM_ALLOC_N( size_t, lhs->dim );

     lhs->shape[0]         = shape(0);
     lhs->shape[1]         = shape(1);

     lhs->offset           = NM_ALLOC_N( size_t, lhs->dim );

     lhs->offset[0]        = 0;
     lhs->offset[1]        = 0;

     lhs->capacity         = new_capacity;
     lhs->dtype            = new_dtype;
     lhs->ndnz             = new_ndnz;
     lhs->ija              = NM_ALLOC_N( size_t, new_capacity );
     lhs->a                = NM_ALLOC_N( E,      new_capacity );
     lhs->src              = lhs;
     lhs->count            = 1;

     return lhs;
   }


  /*
   * Make a full matrix structure copy (entries remain uninitialized). Remember to NM_FREE()!
   */
  template <typename E>
  YALE_STORAGE* alloc_struct_copy(size_t new_capacity) const {
    YALE_STORAGE* lhs     = alloc_basic_copy<E>(new_capacity, count_copy_ndnz());
    // Now copy the IJA contents
    if (slice) {
      rb_raise(rb_eNotImpError, "cannot copy struct due to different offsets");
    } else {
      for (size_t m = 0; m < size(); ++m) {
        lhs->ija[m] = ija(m); // copy indices
      }
    }
    return lhs;
  }


  /*
   * Copy this slice (or the full matrix if it isn't a slice) into a new matrix which is already allocated, ns.
   */
  template <typename E, bool Yield=false>
  void copy(YALE_STORAGE& ns) const {
    //nm::dtype_t new_dtype = nm::ctype_to_dtype_enum<E>::value_type;
    // get the default value for initialization (we'll re-use val for other copies after this)
    E val = static_cast<E>(const_default_obj());

    // initialize the matrix structure and clear the diagonal so we don't have to
    // keep track of unwritten entries.
    YaleStorage<E>::init(ns, &val);

    E* ns_a    = reinterpret_cast<E*>(ns.a);
    size_t sz  = shape(0) + 1; // current used size of ns
    nm_yale_storage_register(&ns);

    // FIXME: If diagonals line up, it's probably faster to do this with stored diagonal and stored non-diagonal iterators
    for (const_row_iterator it = cribegin(); it != criend(); ++it) {
      for (auto jt = it.begin(); !jt.end(); ++jt) {
        if (it.i() == jt.j()) {
          if (Yield)  ns_a[it.i()] = rb_yield(~jt);
          else        ns_a[it.i()] = static_cast<E>(*jt);
        } else if (*jt != const_default_obj()) {
          if (Yield)  ns_a[sz]     = rb_yield(~jt);
          else        ns_a[sz]     = static_cast<E>(*jt);
          ns.ija[sz]    = jt.j();
          ++sz;
        }
      }
      ns.ija[it.i()+1]  = sz;
    }
    nm_yale_storage_unregister(&ns);

    //ns.ija[shape(0)] = sz;                // indicate end of last row
    ns.ndnz          = sz - shape(0) - 1; // update ndnz count
  }


  /*
   * Allocate a casted copy of this matrix/reference. Remember to NM_FREE() the result!
   *
   * If Yield is true, E must be nm::RubyObject, and it will call an rb_yield upon the stored value.
   */
  template <typename E, bool Yield = false>
  YALE_STORAGE* alloc_copy() const {
    //nm::dtype_t new_dtype = nm::ctype_to_dtype_enum<E>::value_type;

    YALE_STORAGE* lhs;
    if (slice) {
      size_t* xshape    = NM_ALLOC_N(size_t, 2);
      xshape[0]         = shape(0);
      xshape[1]         = shape(1);
      size_t ndnz       = count_copy_ndnz();
      size_t reserve    = shape(0) + ndnz + 1;

//      std::cerr << "reserve = " << reserve << std::endl;

      lhs               = YaleStorage<E>::create(xshape, reserve);

      // FIXME: This should probably be a throw which gets caught outside of the object.
      if (lhs->capacity < reserve)
        rb_raise(nm_eStorageTypeError, "conversion failed; capacity of %lu requested, max allowable is %lu", reserve, lhs->capacity);

      // Fill lhs with what's in our current matrix.
      copy<E, Yield>(*lhs);
    } else {
      // Copy the structure and setup the IJA structure.
      lhs               = alloc_struct_copy<E>(s->capacity);

      E* la = reinterpret_cast<E*>(lhs->a);

      nm_yale_storage_register(lhs);
      for (size_t m = 0; m < size(); ++m) {
        if (Yield) {
	  la[m] = rb_yield(nm::yale_storage::nm_rb_dereference(a(m)));
	}
        else       la[m] = static_cast<E>(a(m));
      }
      nm_yale_storage_unregister(lhs);

    }

    return lhs;
  }

  /*
   * Allocate a transposed copy of the matrix
   */
  /*
   * Allocate a casted copy of this matrix/reference. Remember to NM_FREE() the result!
   *
   * If Yield is true, E must be nm::RubyObject, and it will call an rb_yield upon the stored value.
   */
  template <typename E, bool Yield = false>
  YALE_STORAGE* alloc_copy_transposed() const {

    if (slice) {
      rb_raise(rb_eNotImpError, "please make a copy before transposing");
    } else {
      // Copy the structure and setup the IJA structure.
      size_t* xshape    = NM_ALLOC_N(size_t, 2);
      xshape[0]         = shape(1);
      xshape[1]         = shape(0);

      // Take a stab at the number of non-diagonal stored entries we'll have.
      size_t reserve    = size() - xshape[1] + xshape[0];
      YALE_STORAGE* lhs = YaleStorage<E>::create(xshape, reserve);
      E r_init          = static_cast<E>(const_default_obj());
      YaleStorage<E>::init(*lhs, &r_init);

      nm::yale_storage::transpose_yale<D,E,true,true>(shape(0), shape(1), ija_p(), ija_p(), a_p(), const_default_obj(),
                                                      lhs->ija, lhs->ija, reinterpret_cast<E*>(lhs->a), r_init);
      return lhs;
    }

    return NULL;
  }


  /*
   * Comparison between two matrices. Does not check size and such -- assumption is that they are the same shape.
   */
  template <typename E>
  bool operator==(const YaleStorage<E>& rhs) const {
    for (size_t i = 0; i < shape(0); ++i) {
      typename YaleStorage<D>::const_row_iterator li = cribegin(i);
      typename YaleStorage<E>::const_row_iterator ri = rhs.cribegin(i);

      size_t j = 0; // keep track of j so we can compare different defaults

      auto lj = li.begin();
      auto rj = ri.begin();
      while (!lj.end() || !rj.end()) {
        if (lj < rj) {
          if (*lj != rhs.const_default_obj()) return false;
          ++lj;
        } else if (rj < lj) {
          if (const_default_obj() != *rj)     return false;
          ++rj;
        } else { // rj == lj
          if (*lj != *rj) return false;
          ++lj;
          ++rj;
        }
        ++j;
      }

      // if we skip an entry (because it's an ndnz in BOTH matrices), we need to compare defaults.
      // (We know we skipped if lj and rj hit end before j does.)
      if (j < shape(1) && const_default_obj() != rhs.const_default_obj()) return false;

      ++li;
      ++ri;
    }

    return true;
  }

  /*
   * Necessary for element-wise operations. The return dtype will be nm::RUBYOBJ.
   */
  template <typename E>
  VALUE map_merged_stored(VALUE klass, nm::YaleStorage<E>& t, VALUE r_init) const {
    nm_register_value(r_init);
    VALUE s_init    = const_default_value(),
          t_init    = t.const_default_value();
    nm_register_value(s_init);
    nm_register_value(t_init);
    
    // Make a reasonable approximation of the resulting capacity
    size_t s_ndnz   = count_copy_ndnz(),
           t_ndnz   = t.count_copy_ndnz();
    size_t reserve  = shape(0) + std::max(s_ndnz, t_ndnz) + 1;

    size_t* xshape  = NM_ALLOC_N(size_t, 2);
    xshape[0]       = shape(0);
    xshape[1]       = shape(1);

    YALE_STORAGE* rs= YaleStorage<nm::RubyObject>::create(xshape, reserve);

    if (r_init == Qnil) {
      nm_unregister_value(r_init);
      r_init       = rb_yield_values(2, s_init, t_init);
      nm_register_value(r_init);
    }

    nm::RubyObject r_init_obj(r_init);

    // Prepare the matrix structure
    YaleStorage<nm::RubyObject>::init(*rs, &r_init_obj);
    NMATRIX* m     = nm_create(nm::YALE_STORE, reinterpret_cast<STORAGE*>(rs));
    nm_register_nmatrix(m);
    VALUE result   = Data_Wrap_Struct(klass, nm_mark, nm_delete, m);
    nm_unregister_nmatrix(m);
    nm_register_value(result);
    nm_unregister_value(r_init);

    RETURN_SIZED_ENUMERATOR_PRE
    nm_unregister_value(result);
    nm_unregister_value(t_init);
    nm_unregister_value(s_init);
    // No obvious, efficient way to pass a length function as the fourth argument here:
    RETURN_SIZED_ENUMERATOR(result, 0, 0, 0);

    // Create an object for us to iterate over.
    YaleStorage<nm::RubyObject> r(rs);

    // Walk down our new matrix, inserting values as we go.
    for (size_t i = 0; i < xshape[0]; ++i) {
      YaleStorage<nm::RubyObject>::row_iterator   ri = r.ribegin(i);
      typename YaleStorage<D>::const_row_iterator si = cribegin(i);
      typename YaleStorage<E>::const_row_iterator ti = t.cribegin(i);

      auto sj = si.begin();
      auto tj = ti.begin();
      auto rj = ri.ndbegin();

      while (sj != si.end() || tj != ti.end()) {
        VALUE  v;
        size_t j;

        if (sj < tj) {
          v = rb_yield_values(2, ~sj, t_init);
          j = sj.j();
          ++sj;
        } else if (tj < sj) {
          v = rb_yield_values(2, s_init, ~tj);
          j = tj.j();
          ++tj;
        } else {
          v = rb_yield_values(2, ~sj, ~tj);
          j = sj.j();
          ++sj;
          ++tj;
        }

        // FIXME: This can be sped up by inserting all at the same time
        // since it's a new matrix. But that function isn't quite ready
        // yet.
        if (j == i) r.a(i) = v;
        else        rj     = ri.insert(rj, j, v);
        //RB_P(rb_funcall(result, rb_intern("yale_ija"), 0));
      }
    }
    nm_unregister_value(result);
    nm_unregister_value(t_init);
    nm_unregister_value(s_init);

    return result;
  }

protected:
  /*
   * Update row sizes starting with row i
   */
  void update_real_row_sizes_from(size_t real_i, int change) {
    ++real_i;
    for (; real_i <= real_shape(0); ++real_i) {
      ija(real_i) += change;
    }
  }


  /*
   * Like move_right, but also involving a resize. This updates row sizes as well. This version also takes a plan for
   * multiple rows, and tries to do them all in one copy. It's used for multi-row slice-setting.
   *
   * This also differs from update_resize_move in that it resizes to the exact requested size instead of reserving space.
   */
  void update_resize_move_insert(size_t real_i, size_t real_j, size_t* lengths, D* const v, size_t v_size, multi_row_insertion_plan p) {
    size_t sz      = size(); // current size of the storage vectors
    size_t new_cap = sz + p.total_change;

    if (new_cap > real_max_size()) {
      NM_FREE(v);
      rb_raise(rb_eStandardError, "resize caused by insertion of size %d (on top of current size %lu) would have caused yale matrix size to exceed its maximum (%lu)", p.total_change, sz, real_max_size());
    }

    if (s->dtype == nm::RUBYOBJ) {
      nm_register_values(reinterpret_cast<VALUE*>(v), v_size);
    }

    size_t* new_ija     = NM_ALLOC_N( size_t,new_cap );
    D* new_a            = NM_ALLOC_N( D,     new_cap );

    // Copy unchanged row pointers first.
    size_t m = 0;
    for (; m <= real_i; ++m) {
      new_ija[m]        = ija(m);
      new_a[m]          = a(m);
    }

    // Now copy unchanged locations in IJA and A.
    size_t q = real_shape(0)+1; // q is the copy-to position.
    size_t r = real_shape(0)+1; // r is the copy-from position.
    for (; r < p.pos[0]; ++r, ++q) {
      new_ija[q]        = ija(r);
      new_a[q]          = a(r);
    }

    // For each pos and change in the slice, copy the information prior to the insertion point. Then insert the necessary
    // information.
    size_t v_offset = 0;
    int accum = 0; // keep track of the total change as we go so we can update row information.
    for (size_t i = 0; i < lengths[0]; ++i, ++m) {
      for (; r < p.pos[i]; ++r, ++q) {
        new_ija[q]      = ija(r);
        new_a[q]        = a(r);
      }

      // Insert slice data for a single row.
      for (size_t j = 0; j < lengths[1]; ++j, ++v_offset) {
        if (v_offset >= v_size) v_offset %= v_size;

        if (j + real_j == i + real_i) { // modify diagonal
          new_a[real_i + i] = v[v_offset];
        } else if (v[v_offset] != const_default_obj()) {
          new_ija[q]        = j + real_j;
          new_a[q]          = v[v_offset];
          ++q; // move on to next q location
        }

        if (r < ija(real_shape(0)) && ija(r) == j + real_j) ++r; // move r forward if the column matches.
      }

      // Update the row pointer for the current row.
      accum                += p.change[i];
      new_ija[m]            = ija(m) + accum;
      new_a[m]              = a(m); // copy diagonal for this row
    }

    // Now copy everything subsequent to the last insertion point.
    for (; r < size(); ++r, ++q) {
      new_ija[q]            = ija(r);
      new_a[q]              = a(r);
    }

    // Update the remaining row pointers and copy remaining diagonals
    for (; m <= real_shape(0); ++m) {
      new_ija[m]            = ija(m) + accum;
      new_a[m]              = a(m);
    }

    s->capacity = new_cap;

    NM_FREE(s->ija);
    NM_FREE(s->a);

    if (s->dtype == nm::RUBYOBJ) {
      nm_unregister_values(reinterpret_cast<VALUE*>(v), v_size);
    }   

    s->ija      = new_ija;
    s->a        = reinterpret_cast<void*>(new_a);
  }




  /*
   * Like move_right, but also involving a resize. This updates row sizes as well.
   */
  void update_resize_move(row_stored_nd_iterator position, size_t real_i, int n) {
    size_t sz      = size(); // current size of the storage vectors
    size_t new_cap = n > 0 ? capacity() * nm::yale_storage::GROWTH_CONSTANT
                           : capacity() / nm::yale_storage::GROWTH_CONSTANT;
    size_t max_cap = real_max_size();

    if (new_cap > max_cap) {
      new_cap = max_cap;
      if (sz + n > max_cap)
        rb_raise(rb_eStandardError, "resize caused by insertion/deletion of size %d (on top of current size %lu) would have caused yale matrix size to exceed its maximum (%lu)", n, sz, real_max_size());
    }

    if (new_cap < sz + n) new_cap = sz + n;

    size_t* new_ija     = NM_ALLOC_N( size_t,new_cap );
    D* new_a            = NM_ALLOC_N( D,     new_cap );

    // Copy unchanged row pointers first.
    for (size_t m = 0; m <= real_i; ++m) {
      new_ija[m]        = ija(m);
      new_a[m]          = a(m);
    }

    // Now update row pointers following the changed row as we copy the additional values.
    for (size_t m = real_i + 1; m <= real_shape(0); ++m) {
      new_ija[m]        = ija(m) + n;
      new_a[m]          = a(m);
    }

    // Copy all remaining prior to insertion/removal site
    for (size_t m = real_shape(0) + 1; m < position.p(); ++m) {
      new_ija[m]        = ija(m);
      new_a[m]          = a(m);
    }

    // Copy all subsequent to insertion/removal site
    size_t m = position.p();
    if (n < 0) m -= n;

    for (; m < sz; ++m) {
      new_ija[m+n]      = ija(m);
      new_a[m+n]        = a(m);
    }

    if (s->dtype == nm::RUBYOBJ) {
      nm_yale_storage_register_a(new_a, new_cap);
    }

    s->capacity = new_cap;

    NM_FREE(s->ija);
    NM_FREE(s->a);

    if (s->dtype == nm::RUBYOBJ) {
      nm_yale_storage_unregister_a(new_a, new_cap);
    }

    s->ija      = new_ija;
    s->a        = reinterpret_cast<void*>(new_a);
  }


  /*
   * Move elements in the IJA and A arrays by n (to the right).
   * Does not update row sizes.
   */
  void move_right(row_stored_nd_iterator position, size_t n) {
    size_t sz = size();
    for (size_t m = 0; m < sz - position.p(); ++m) {
      ija(sz+n-1-m) = ija(sz-1-m);
      a(sz+n-1-m)   = a(sz-1-m);
    }
  }

  /*
   * Move elements in the IJA and A arrays by n (to the left). Here position gives
   * the location to move to, and they should come from n to the right.
   */
  void move_left(row_stored_nd_iterator position, size_t n) {
    size_t sz = size();
    for (size_t m = position.p() + n; m < sz; ++m) {   // work backwards
      ija(m-n)      = ija(m);
      a(m-n)        = a(m);
    }
  }

  YALE_STORAGE* s;
  bool          slice;
  size_t*       slice_shape;
  size_t*       slice_offset;
};

} // end of nm namespace

#endif // YALE_CLASS_H
