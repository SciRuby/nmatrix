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
// == row_stored.h
//
// Iterator for traversing a single stored row of a matrix (needed
// for row.h). FIXME: This is not as efficient as it could be; it uses
// two binary searches to find the beginning and end of each slice.
// The end search shouldn't be necessary, but I couldn't make it
// work without it, and eventually decided my dissertation should
// be a priority.
//

#ifndef YALE_ITERATORS_ROW_STORED_H
# define YALE_ITERATORS_ROW_STORED_H

#include <stdexcept>

namespace nm { namespace yale_storage {


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
class row_stored_iterator_T : public row_stored_nd_iterator_T<D,RefType,YaleRef,RowRef> {
protected:
  using row_stored_nd_iterator_T<D,RefType,YaleRef,RowRef>::r;
  using row_stored_nd_iterator_T<D,RefType,YaleRef,RowRef>::p_;
  bool d_visited, d;

public:

  // end_ is necessary for the logic when a row is empty other than the diagonal. If we just
  // relied on pp == last_p+1, it'd look like these empty rows were actually end() iterators.
  // So we have to actually mark end_ by telling it to ignore that diagonal visitation.
  row_stored_iterator_T(RowRef& row, size_t pp, bool end_ = false)
  : row_stored_nd_iterator_T<D,RefType,YaleRef,RowRef>(row, pp),
    d_visited(!row.has_diag()), // if the row has no diagonal, just marked it as visited.
    d(r.is_diag_first() && !end_)        // do we start at the diagonal?
  {
  }

  /* Diagonal constructor. Puts us on the diagonal (unless end is true) */
  /*row_stored_iterator_T(RowRef& row, bool end_, size_t j)
  : row_stored_nd_iterator_T<D,RefType,YaleRef,RowRef>(row.ndfind(j)),
    d_visited(false),
    d(!end_ && j + row.offset(1) == row.real_i())
  { }*/

  virtual bool diag() const {
    return d;
  }

  virtual bool end() const {
    return !d && p_ > r.p_last;
  }

  row_stored_iterator_T<D,RefType,YaleRef,RowRef>& operator++() {
    if (end()) throw std::out_of_range("cannot increment row stored iterator past end of stored row");
    if (d) {
      d_visited = true;
      d         = false;
    } else {
      ++p_;
      // Are we at a diagonal?
      // If we hit the end or reach a point where j > diag_j, and still
      // haven't visited the diagonal, we should do so before continuing.
      if (!d_visited && (end() || j() > r.diag_j())) {
        d = true;
      }
    }

    return *this;
  }

  row_stored_iterator_T<D,RefType,YaleRef,RowRef> operator++(int dummy) const {
    row_stored_iterator_T<D,RefType,YaleRef,RowRef> r(*this);
    return ++r;
  }

  size_t j() const {
    if (end()) throw std::out_of_range("cannot dereference an end pointer");
    return (d ? r.p_diag() : r.ija(p_)) - r.offset(1);
  }

  // Need to declare all row_stored_iterator_T friends of each other.
  template <typename E, typename ERefType, typename EYaleRef, typename ERowRef> friend class row_stored_iterator_T;

  // De-reference the iterator
  RefType& operator*()       {
    return d ? r.a(r.p_diag()) : r.a(p_);
  }

  RefType& operator*() const {
    return d ? r.a(r.p_diag()) : r.a(p_);
  }

  // Ruby VALUE de-reference
  VALUE operator~() const {
    return nm_rb_dereference<D>(**this);
  }

};

}} // end of namespace nm::yale_storage

#endif // YALE_ITERATORS_ROW_STORED_H