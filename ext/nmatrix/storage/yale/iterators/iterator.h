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
// == iterator.h
//
// Iterate over yale as if dense
//

#ifndef YALE_ITERATORS_ITERATOR_H
# define YALE_ITERATORS_ITERATOR_H

#include <type_traits>
#include <typeinfo>

namespace nm { namespace yale_storage {

/*
 * Iterator for traversing matrix class as if it were dense (visits each entry in order).
 */
template <typename D,
          typename RefType,
          typename YaleRef = typename std::conditional<
            std::is_const<RefType>::value,
            const nm::YaleStorage<D>,
            nm::YaleStorage<D>
          >::type>
class iterator_T : public basic_iterator_T<D,RefType,YaleRef> {
  using basic_iterator_T<D,RefType,YaleRef>::i_;
  using basic_iterator_T<D,RefType,YaleRef>::p_;
  using basic_iterator_T<D,RefType,YaleRef>::y;
  using basic_iterator_T<D,RefType,YaleRef>::offset;
  using basic_iterator_T<D,RefType,YaleRef>::shape;
  using basic_iterator_T<D,RefType,YaleRef>::ija;

protected:
  size_t j_; // These are relative to the slice.

public:
  // Create an iterator. May select the row since this is O(1).
  iterator_T(YaleRef& obj, size_t ii = 0)
  : basic_iterator_T<D,RefType,YaleRef>(obj, ii, obj.ija(ii + obj.offset(0))), j_(0)
  {
    // advance to the beginning of the row
    if (obj.offset(1) > 0)
      p_ = y.find_pos_for_insertion(i_,j_);
  }

  // Prefix ++
  iterator_T<D,RefType,YaleRef>& operator++() {
    size_t prev_j = j_++;
    if (j_ >= shape(1)) {
      j_ = 0;
      ++i_;

      // Do a binary search to find the beginning of the slice
      p_ = offset(0) > 0 ? y.find_pos_for_insertion(i_,j_) : ija(i_);
    } else {
      // If the last j was actually stored in this row of the matrix, need to advance p.

      if (!y.real_row_empty(i_ + offset(0)) && ija(p_) <= prev_j + offset(1)) ++p_;  // this test is the same as real_ndnz_exists
    }

    return *this;
  }

  iterator_T<D,RefType,YaleRef> operator++(int dummy) const {
    iterator_T<D,RefType,YaleRef> iter(*this);
    return ++iter;
  }

  virtual bool operator!=(const iterator_T<D,RefType,YaleRef>& rhs) const {
    return this->dense_location() != rhs.dense_location();
  }

  virtual bool operator==(const iterator_T<D,RefType,YaleRef>& rhs) const {
    return this->dense_location() == rhs.dense_location();
  }

  bool operator<(const iterator_T<D,RefType,YaleRef>& rhs) const {
    return this->dense_location() < rhs.dense_location();
  }

  bool operator>(const iterator_T<D,RefType,YaleRef>& rhs) const {
    return this->dense_location() > rhs.dense_location();
  }

  virtual bool diag() const { return i_ + offset(0) == j_ + offset(1); }

  // De-reference
  RefType& operator*() {
    if (diag())                                                                return y.a( i_ + offset(0) );
    else if (p_ >= ija(i_+offset(0)+1))                                        return y.const_default_obj();
    else if (!y.real_row_empty(i_ + offset(0)) && ija(p_) == j_ + offset(1))   return y.a( p_ );
    else                                                                       return y.const_default_obj();
  }

  RefType& operator*() const {
    if (diag())                                                                return y.a( i_ + offset(0) );
    else if (p_ >= ija(i_+offset(0)+1))                                        return y.const_default_obj();
    else if (!y.real_row_empty(i_ + offset(0)) && ija(p_) == j_ + offset(1))   return y.a( p_ );
    else                                                                       return y.const_default_obj();
  }

  virtual size_t j() const { return j_; }
};


} } // end of namespace nm::yale_storage

#endif // YALE_ITERATORS_ITERATOR_H