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
// == base.h
//
// Yale storage pure virtual basic_iterator class.
//

#ifndef YALE_ITERATORS_BASE_H
# define YALE_ITERATORS_BASE_H

#include <type_traits>
#include <typeinfo>
#include <stdexcept>

namespace nm {

template <typename D> class YaleStorage;

namespace yale_storage {

template <typename D>
VALUE nm_rb_dereference(D const& v) {
  return nm::RubyObject(v).rval;
}

template <>
VALUE nm_rb_dereference<nm::RubyObject>(nm::RubyObject const& v) {
  return v.rval;
}

/*
 * Iterator base class (pure virtual).
 */
template <typename D,
          typename RefType,
          typename YaleRef = typename std::conditional<
            std::is_const<RefType>::value,
            const nm::YaleStorage<D>,
            nm::YaleStorage<D>
          >::type>
class basic_iterator_T {

protected:
  YaleRef& y;
  size_t i_;
  size_t p_;

public:
  size_t offset(size_t d) const { return y.offset(d); }
  size_t shape(size_t d) const { return y.shape(d); }
  size_t real_shape(size_t d) const { return y.real_shape(d); }

  size_t dense_location() const {
    return i()*shape(1) + j();
  }

  template <typename T = typename std::conditional<std::is_const<RefType>::value, const size_t, size_t>::type>
  T& ija(size_t pp) const { return y.ija(pp); }

  template <typename T = typename std::conditional<std::is_const<RefType>::value, const size_t, size_t>::type>
  T& ija(size_t pp) { return y.ija(pp); }

  virtual bool diag() const {
    return p_ < std::min(y.real_shape(0), y.real_shape(1));
  }
  virtual bool done_with_diag() const {
    return p_ == std::min(y.real_shape(0), y.real_shape(1));
  }
  virtual bool nondiag() const {
    return p_ > std::min(y.real_shape(0), y.real_shape(1));
  }

  basic_iterator_T(YaleRef& obj, size_t ii = 0, size_t pp = 0) : y(obj), i_(ii), p_(pp) { }

  basic_iterator_T<D,RefType,YaleRef>& operator=(const basic_iterator_T<D,RefType,YaleRef>& rhs) {
    if (&y != &(rhs.y)) throw std::logic_error("can only be used on iterators with the same matrix");
    i_ = rhs.i_;
    p_ = rhs.p_;
    return *this;
  }

  virtual inline size_t i() const { return i_; }
  virtual size_t j() const = 0;

  virtual inline VALUE rb_i() const { return LONG2NUM(i()); }
  virtual inline VALUE rb_j() const { return LONG2NUM(j()); }

  virtual size_t real_i() const { return offset(0) + i(); }
  virtual size_t real_j() const { return offset(1) + j(); }
  virtual size_t p() const { return p_; }
  virtual bool real_ndnz_exists() const { return !y.real_row_empty(real_i()) && ija(p_) == real_j(); }

  virtual RefType& operator*() = 0;
  virtual RefType& operator*() const = 0;


  // Ruby VALUE de-reference
  inline VALUE operator~() const {
    return nm_rb_dereference<D>(**this);
  //virtual VALUE operator~() const {
  //  if (typeid(D) == typeid(RubyObject)) return (**this); // FIXME: return rval instead, faster;
  //  else return RubyObject(*(*this)).rval;
  }

  virtual bool operator==(const std::pair<size_t,size_t>& ij) {
    if (p() >= ija(real_shape(0))) return false;
    else return i() == ij.first && j() == ij.second;
  }

  virtual bool operator==(const basic_iterator_T<D,RefType,YaleRef>& rhs) const {
    return i() == rhs.i() && j() == rhs.j();
  }
  virtual bool operator!=(const basic_iterator_T<D,RefType,YaleRef>& rhs) const {
    return i() != rhs.i() || j() != rhs.j();
  }
};


} } // end of namespace nm::yale_storage

#endif // YALE_ITERATORS_BASE_H