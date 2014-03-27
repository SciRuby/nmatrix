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
// == stored_diagonal_iterator.h
//
// Yale storage diagonal-storage iterator
//

#ifndef YALE_ITERATORS_STORED_DIAGONAL_H
# define YALE_ITERATORS_STORED_DIAGONAL_H

#include <type_traits>
#include <typeinfo>

namespace nm { namespace yale_storage {

/*
 * Iterate across the stored diagonal.
 */
template <typename D,
          typename RefType,
          typename YaleRef = typename std::conditional<
            std::is_const<RefType>::value,
            const nm::YaleStorage<D>,
            nm::YaleStorage<D>
          >::type>
class stored_diagonal_iterator_T : public basic_iterator_T<D,RefType,YaleRef> {
  using basic_iterator_T<D,RefType,YaleRef>::p_;
  using basic_iterator_T<D,RefType,YaleRef>::y;
  using basic_iterator_T<D,RefType,YaleRef>::offset;
  using basic_iterator_T<D,RefType,YaleRef>::shape;
public:
  stored_diagonal_iterator_T(YaleRef& obj, size_t d = 0)
  : basic_iterator_T<D,RefType,YaleRef>(obj,                // y
                   std::max(obj.offset(0), obj.offset(1)) + d - obj.offset(0), // i_
                   std::max(obj.offset(0), obj.offset(1)) + d) // p_
  {
//      std::cerr << "sdbegin: d=" << d << ", p_=" << p_ << ", i()=" << i() << ", j()=" << j() << std::endl;
    // p_ can range from max(y.offset(0), y.offset(1)) to min(y.real_shape(0), y.real_shape(1))
  }


  size_t d() const {
    return p_ - std::max(offset(0), offset(1));
  }

  stored_diagonal_iterator_T<D,RefType,YaleRef>& operator++() {
    if (i() < shape(0)) ++p_;
    return *this;
  }

  stored_diagonal_iterator_T<D,RefType,YaleRef> operator++(int dummy) const {
    stored_diagonal_iterator_T<D,RefType,YaleRef> iter(*this);
    return ++iter;
  }

  // Indicates if we're at the end of the iteration.
  bool end() const {
    return p_ >= std::min( shape(0) + offset(0), shape(1) + offset(1) );
  }

  // i() and j() are how we know if we're past-the-end. i will be shape(0) and j will be 0.
  size_t i() const {
    return p_ - offset(0);
  }

  size_t j() const {
    return p_ - offset(1);
  }


  template <typename E, typename ERefType = typename std::conditional<std::is_const<RefType>::value, const E, E>::type>
  bool operator!=(const stored_diagonal_iterator_T<E,ERefType>& rhs) const { return d() != rhs.d(); }

  template <typename E, typename ERefType = typename std::conditional<std::is_const<RefType>::value, const E, E>::type>
  bool operator==(const stored_diagonal_iterator_T<E,ERefType>& rhs) const { return !(*this != rhs); }

  template <typename E, typename ERefType = typename std::conditional<std::is_const<RefType>::value, const E, E>::type>
  bool operator<(const stored_diagonal_iterator_T<E,ERefType>& rhs) const {  return d() < rhs.d(); }

  template <typename E, typename ERefType = typename std::conditional<std::is_const<RefType>::value, const E, E>::type>
  bool operator<=(const stored_diagonal_iterator_T<E,ERefType>& rhs) const {
    return d() <= rhs.d();
  }

  template <typename E, typename ERefType = typename std::conditional<std::is_const<RefType>::value, const E, E>::type>
  bool operator>(const stored_diagonal_iterator_T<E,ERefType>& rhs) const {
    return d() > rhs.d();
  }

  template <typename E, typename ERefType = typename std::conditional<std::is_const<RefType>::value, const E, E>::type>
  bool operator>=(const stored_diagonal_iterator_T<E,ERefType>& rhs) const {
    return d() >= rhs.d();
  }

  RefType& operator*() { return y.a(p_); }
  RefType& operator*() const { return y.a(p_); }

};

} } // end of namespace nm::yale_storage

#endif // YALE_ITERATORS_STORED_DIAGONAL_H