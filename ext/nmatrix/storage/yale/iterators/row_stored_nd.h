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
// == row_stored_nd.h
//
// Yale storage row-by-row nondiagonal-storage iterator
//

#ifndef YALE_ITERATORS_ROW_STORED_ND_H
# define YALE_ITERATORS_ROW_STORED_ND_H

#include <type_traits>
#include <typeinfo>
#include <stdexcept>

namespace nm { namespace yale_storage {

/*
 * Constants
 */
const float GROWTH_CONSTANT = 1.5;


/*
 * Forward declarations
 */
template <typename D, typename RefType, typename YaleRef> class row_iterator_T;

/*
 * Iterator for visiting each stored element in a row, including diagonals.
 */
template <typename D,
          typename RefType,
          typename YaleRef = typename std::conditional<
            std::is_const<RefType>::value,
            const nm::YaleStorage<D>,
            nm::YaleStorage<D>
          >::type,
          typename RowRef = typename std::conditional<
            std::is_const<RefType>::value,
            const row_iterator_T<D,RefType,YaleRef>,
            row_iterator_T<D,RefType,YaleRef>
          >::type>
class row_stored_nd_iterator_T {
protected:
  RowRef& r;
  size_t p_;

public:

  row_stored_nd_iterator_T(RowRef& row, size_t pp)
  : r(row),
    p_(pp)        // do we start at the diagonal?
  {
  }

  // DO NOT IMPLEMENT THESE FUNCTIONS. They prevent C++ virtual slicing
  //template <typename T> row_stored_nd_iterator_T(T const& rhs);
  //template <typename T> row_stored_nd_iterator_T<D,RefType,YaleRef,RowRef> const& operator=(T const& rhs);

  // Next two functions are to ensure we can still cast between nd iterators.
  row_stored_nd_iterator_T(row_stored_nd_iterator_T<D,RefType,YaleRef,RowRef> const& rhs)
  : r(rhs.r), p_(rhs.p_)
  { }

  row_stored_nd_iterator_T<D,RefType,YaleRef,RowRef> const& operator=(row_stored_nd_iterator_T<D,RefType,YaleRef,RowRef> const& rhs) {
    if (&r != &(rhs.r))
      throw std::logic_error("can't assign iterator from another row iterator");
    p_ = rhs.p_;
    return *this;
  }

  virtual size_t p() const { return p_; }

  virtual bool end() const {
    return p_ > r.p_last;
  }

  row_stored_nd_iterator_T<D,RefType,YaleRef,RowRef>& operator++() {
    if (end()) throw std::out_of_range("cannot increment row stored iterator past end of stored row");
    ++p_;

    return *this;
  }

  row_stored_nd_iterator_T<D,RefType,YaleRef,RowRef> operator++(int dummy) const {
    row_stored_nd_iterator_T<D,RefType,YaleRef,RowRef> r(*this);
    return ++r;
  }

  virtual size_t j() const {
    if (end()) throw std::out_of_range("cannot dereference (get j()) for an end pointer");
    return r.ija(p_) - r.offset(1);
  }

  // Need to declare all row_stored_nd_iterator_T friends of each other.
  template <typename E, typename ERefType, typename EYaleRef, typename ERowRef> friend class row_stored_nd_iterator_T;


  virtual bool operator==(const row_stored_nd_iterator_T<D,RefType>& rhs) const {
    if (r.i() != rhs.r.i())     return false;
    if (end())                  return rhs.end();
    else if (rhs.end())         return false;
    return j() == rhs.j();
  }

  // There is something wrong with this function.
  virtual bool operator!=(const row_stored_nd_iterator_T<D,RefType>& rhs) const {
    if (r.i() != rhs.r.i()) return true;
    if (end())              return !rhs.end();
    else if (rhs.end())     return true;
    return j() != rhs.j();
  }

  template <typename E, typename ERefType = typename std::conditional<std::is_const<RefType>::value, const E, E>::type>
  bool operator<(const row_stored_nd_iterator_T<E,ERefType>& rhs) const {
    if (r < rhs.r)      return true;
    if (r > rhs.r)      return false;

    // r == rhs.r
    if (end())        return false;
    if (rhs.end())    return true;
    return j() < rhs.j();
  }

  // De-reference the iterator
  RefType& operator*()       {
    return r.a(p_);
  }

  RefType& operator*() const {
    return r.a(p_);
  }

  // Ruby VALUE de-reference
  VALUE operator~() const {
    return nm_rb_dereference<D>(**this);
  }

  inline virtual VALUE rb_j() const { return LONG2NUM(j()); }

};



} } // end of namespace nm::yale_storage

#endif // YALE_ITERATORS_ROW_STORED_ND_H
