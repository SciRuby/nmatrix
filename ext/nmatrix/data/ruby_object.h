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
// == ruby_object.h
//
// Functions and classes for dealing with Ruby objects.

#ifndef RUBY_OBJECT_H
#define RUBY_OBJECT_H

/*
 * Standard Includes
 */

#include <ruby.h>
#include <iostream>
#include <type_traits>

/*
 * Project Includes
 */

#include "ruby_constants.h"

/*
 * Macros
 */
#define NM_RUBYVAL_IS_NUMERIC(val)                (FIXNUM_P(val) or (TYPE(val) == T_FLOAT) or (TYPE(val) == T_COMPLEX) or (TYPE(val) == T_RATIONAL))
#define NMATRIX_CHECK_TYPE(val) \
	if (TYPE(val) != T_DATA || (RDATA(val)->dfree != (RUBY_DATA_FUNC)nm_delete && RDATA(val)->dfree != (RUBY_DATA_FUNC)nm_delete_ref)) \
		rb_raise(rb_eTypeError, "Expected NMatrix on left-hand side of operation.");

/*
 * Classes and Functions
 */

namespace nm {
template<typename T, typename U>
struct made_from_same_template : std::false_type {}; 
 
template<template<typename> class Templ, typename Arg1, typename Arg2>
struct made_from_same_template<Templ<Arg1>, Templ<Arg2>> : std::true_type {};

class RubyObject {
	public:
	VALUE rval;
	
	/*
	 * Value constructor.
	 */
	inline RubyObject(VALUE ref = Qnil) : rval(ref) {}
	
	/*
	 * Complex number constructor.
	 */
	template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
	inline RubyObject(const Complex<FloatType>& other) : rval(rb_complex_new(rb_float_new(other.r), rb_float_new(other.i))) {}
	
	/*
	 * Rational number constructor.
	 */
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline RubyObject(const Rational<IntType>& other) : rval(rb_rational_new(INT2FIX(other.n), INT2FIX(other.d))) {}
	
	/*
	 * Integer constructor.
	 *
	 * Does not work as a template.
	 */
	inline RubyObject(uint8_t other)  : rval(INT2FIX(other)) {}
	inline RubyObject(int8_t other)   : rval(INT2FIX(other)) {}
	inline RubyObject(int16_t other)  : rval(INT2FIX(other)) {}
	inline RubyObject(uint16_t other) : rval(INT2FIX(other)) {}
	inline RubyObject(int32_t other)  : rval(INT2FIX(other)) {}
	// there is no uint32_t here because that's a Ruby VALUE type, and we need the compiler to treat that as a VALUE.
	inline RubyObject(int64_t other)  : rval(INT2FIX(other)) {}
//	inline RubyObject(uint64_t other) : rval(INT2FIX(other)) {}


	/*
	 * Float constructor.
	 *
	 * Does not work as a template.
	 */
	inline RubyObject(float other)   : rval(rb_float_new(other)) {}
	inline RubyObject(double other)  : rval(rb_float_new(other)) {}

  /*
   * Operators for converting RubyObjects to other C types.
   */

#define RETURN_OBJ2NUM(mac)   if (this->rval == Qtrue) return 1; else if (this->rval == Qfalse) return 0; else return mac(this->rval);

  inline operator int8_t()  const { RETURN_OBJ2NUM(NUM2INT)         }
  inline operator uint8_t() const { RETURN_OBJ2NUM(NUM2UINT)        }
  inline operator int16_t() const { RETURN_OBJ2NUM(NUM2INT)         }
  inline operator uint16_t() const { RETURN_OBJ2NUM(NUM2UINT)       }
  inline operator int32_t() const { RETURN_OBJ2NUM(NUM2LONG)        }
  inline operator VALUE() const { return rval; }
  //inline operator uint32_t() const { return NUM2ULONG(this->rval);      }
  inline operator int64_t() const { RETURN_OBJ2NUM(NUM2LONG)        }
  inline operator uint64_t() const { RETURN_OBJ2NUM(NUM2ULONG)      }
  inline operator double()   const { RETURN_OBJ2NUM(NUM2DBL)        }
  inline operator float()  const { RETURN_OBJ2NUM(NUM2DBL)          }

  inline operator Rational32() const { return this->to<Rational32>();   }
  inline operator Rational64() const { return this->to<Rational64>();   }
  inline operator Rational128() const { return this->to<Rational128>(); }

  inline operator Complex64() const { return this->to<Complex64>(); }
  inline operator Complex128() const { return this->to<Complex128>(); }
  /*
	 * Copy constructors.
	 */
	inline RubyObject(const RubyObject& other) : rval(other.rval) {}

  /*
   * Inverse operator.
   */
	inline RubyObject inverse() const {
	  rb_raise(rb_eNotImpError, "RubyObject#inverse needs to be implemented");
	}

	/*
	 * Absolute value.
	 */
	inline RubyObject abs() const {
	  return RubyObject(rb_funcall(this->rval, rb_intern("abs"), 0));
	}

	/*
	 * Binary operator definitions.
	 */
	
	inline RubyObject operator+(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->rval, nm_rb_add, 1, other.rval));
	}

	inline RubyObject& operator+=(const RubyObject& other) {
    this->rval = rb_funcall(this->rval, nm_rb_add, 1, other.rval);
    return *this;
	}

	inline RubyObject operator-(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->rval, nm_rb_sub, 1, other.rval));
	}

	inline RubyObject& operator-=(const RubyObject& other) {
    this->rval = rb_funcall(this->rval, nm_rb_sub, 1, other.rval);
    return *this;
	}
	
	inline RubyObject operator*(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->rval, nm_rb_mul, 1, other.rval));
	}

	inline RubyObject& operator*=(const RubyObject& other) {
    this->rval = rb_funcall(this->rval, nm_rb_mul, 1, other.rval);
    return *this;
	}
	
	inline RubyObject operator/(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->rval, nm_rb_div, 1, other.rval));
	}

	inline RubyObject& operator/=(const RubyObject& other) {
    this->rval = rb_funcall(this->rval, nm_rb_div, 1, other.rval);
    return *this;
	}
	
	inline RubyObject operator%(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->rval, nm_rb_percent, 1, other.rval));
	}
	
	inline bool operator>(const RubyObject& other) const {
		return rb_funcall(this->rval, nm_rb_gt, 1, other.rval) == Qtrue;
	}
	
	inline bool operator<(const RubyObject& other) const {
		return rb_funcall(this->rval, nm_rb_lt, 1, other.rval) == Qtrue;
	}

	template <typename OtherType>
	inline bool operator<(const OtherType& other) const {
		return *this < RubyObject(other);
	}
	
	inline bool operator==(const RubyObject& other) const {
		return rb_funcall(this->rval, nm_rb_eql, 1, other.rval) == Qtrue;
	}

	template <typename OtherType>
	inline bool operator==(const OtherType& other) const {
		return *this == RubyObject(other);
	}
	
	inline bool operator!=(const RubyObject& other) const {
		return rb_funcall(this->rval, nm_rb_neql, 1, other.rval) == Qtrue;
	}

	template <typename OtherType>
	inline bool operator!=(const OtherType& other) const {
		return *this != RubyObject(other);
	}
	
	inline bool operator>=(const RubyObject& other) const {
		return rb_funcall(this->rval, nm_rb_gte, 1, other.rval) == Qtrue;
	}

	template <typename OtherType>
	inline bool operator>=(const OtherType& other) const {
		return *this >= RubyObject(other);
	}
	
	inline bool operator<=(const RubyObject& other) const {
		return rb_funcall(this->rval, nm_rb_lte, 1, other.rval) == Qtrue;
	}

	template <typename OtherType>
	inline bool operator<=(const OtherType& other) const {
		return *this <= RubyObject(other);
	}

	////////////////////////////
	// RUBY-NATIVE OPERATIONS //
	////////////////////////////
/*
	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator==(const NativeType& other) const {
		return *this == RubyObject(other);
	}

  template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator!=(const NativeType& other) const {
		return *this != RubyObject(other);
	}
*/
	//////////////////////////////
	// RUBY-RATIONAL OPERATIONS //
	//////////////////////////////

	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator==(const Rational<IntType>& other) const {
		return *this == RubyObject(other);
	}

  template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator!=(const Rational<IntType>& other) const {
		return *this != RubyObject(other);
	}

	//////////////////////////////
	// RUBY-COMPLEX OPERATIONS //
	//////////////////////////////

	template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
	inline bool operator==(const Complex<FloatType>& other) const {
		return *this == RubyObject(other);
	}

  template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
	inline bool operator!=(const Complex<FloatType>& other) const {
		return *this != RubyObject(other);
	}

	/*
	 * Convert a Ruby object to an integer.
	 */
	template <typename IntType>
	inline typename std::enable_if<std::is_integral<IntType>::value, IntType>::type to(void) {
		return NUM2INT(this->rval);
	}
	
	/*
	 * Convert a Ruby object to a floating point number.
	 */
	template <typename FloatType>
	inline typename std::enable_if<std::is_floating_point<FloatType>::value, FloatType>::type to(void) {
		return NUM2DBL(this->rval);
	}
	
	/*
	 * Convert a Ruby object to a complex number.
	 */
	template <typename ComplexType>
	inline typename std::enable_if<made_from_same_template<ComplexType, Complex64>::value, ComplexType>::type to(void) const {
		if (FIXNUM_P(this->rval) or TYPE(this->rval) == T_FLOAT or TYPE(this->rval) == T_RATIONAL) {
			return ComplexType(NUM2DBL(this->rval));
			
		} else if (TYPE(this->rval) == T_COMPLEX) {
			return ComplexType(NUM2DBL(rb_funcall(this->rval, nm_rb_real, 0)), NUM2DBL(rb_funcall(this->rval, nm_rb_imag, 0)));
			
		} else {
			rb_raise(rb_eTypeError, "Invalid conversion to Complex type.");
		}
	}
	
	/*
	 * Convert a Ruby object to a rational number.
	 */
	template <typename RationalType>
	inline typename std::enable_if<made_from_same_template<RationalType, Rational32>::value, RationalType>::type to(void) const {
		if (FIXNUM_P(this->rval) or TYPE(this->rval) == T_FLOAT or TYPE(this->rval) == T_COMPLEX) {
			return RationalType(NUM2INT(this->rval));
			
		} else if (TYPE(this->rval) == T_RATIONAL) {
			return RationalType(NUM2INT(rb_funcall(this->rval, nm_rb_numer, 0)), NUM2INT(rb_funcall(this->rval, nm_rb_denom, 0)));
			
		} else {
			rb_raise(rb_eTypeError, "Invalid conversion to Rational type.");
		}
	}

};

// Negative operator
inline RubyObject operator-(const RubyObject& rhs) {
  return RubyObject(rb_funcall(rhs.rval, nm_rb_negate, 0));
}


////////////////////////////
// NATIVE-RUBY OPERATIONS //
////////////////////////////

template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline RubyObject operator/(const NativeType left, const RubyObject& right) {
  return RubyObject(left) / right;
}

template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator==(const NativeType left, const RubyObject& right) {
  return RubyObject(left) == right;
}

template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator!=(const NativeType left, const RubyObject& right) {
  return RubyObject(left) != right;
}

template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator<=(const NativeType left, const RubyObject& right) {
  return RubyObject(left) <= right;
}

template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator>=(const NativeType left, const RubyObject& right) {
  return RubyObject(left) >= right;
}

template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator<(const NativeType left, const RubyObject& right) {
  return RubyObject(left) < right;
}

template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator>(const NativeType left, const RubyObject& right) {
  return RubyObject(left) > right;
}


/////////////////////////////
// COMPLEX-RUBY OPERATIONS //
/////////////////////////////

template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline bool operator==(const Complex<FloatType>& left, const RubyObject& right) {
	return RubyObject(left) == right;
}

template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline bool operator!=(const Complex<FloatType>& left, const RubyObject& right) {
	return RubyObject(left) != right;
}

template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline bool operator<=(const Complex<FloatType>& left, const RubyObject& right) {
	return RubyObject(left) <= right;
}

template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline bool operator>=(const Complex<FloatType>& left, const RubyObject& right) {
	return RubyObject(left) >= right;
}

template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline bool operator<(const Complex<FloatType>& left, const RubyObject& right) {
	return RubyObject(left) < right;
}

template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline bool operator>(const Complex<FloatType>& left, const RubyObject& right) {
	return RubyObject(left) > right;
}



//////////////////////////////
// RATIONAL-RUBY OPERATIONS //
//////////////////////////////

template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator==(const Rational<IntType>& left, const RubyObject& right) {
	return RubyObject(left) == right;
}

template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator!=(const Rational<IntType>& left, const RubyObject& right) {
	return RubyObject(left) != right;
}

template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator>=(const Rational<IntType>& left, const RubyObject& right) {
	return RubyObject(left) >= right;
}

template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator<=(const Rational<IntType>& left, const RubyObject& right) {
	return RubyObject(left) <= right;
}

template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator<(const Rational<IntType>& left, const RubyObject& right) {
	return RubyObject(left) < right;
}

template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator>(const Rational<IntType>& left, const RubyObject& right) {
	return RubyObject(left) > right;
}

inline std::ostream& operator<<(std::ostream& out, const RubyObject& rhs) {
  out << "RUBYOBJECT" << std::flush; // FIXME: Try calling inspect or something on the Ruby object if we really need to debug it.
  return out;
}

} // end of namespace nm

namespace std {
  inline nm::RubyObject abs(const nm::RubyObject& obj) {
    return obj.abs();
  }


  inline nm::RubyObject sqrt(const nm::RubyObject& obj) {
    VALUE cMath = rb_const_get(rb_cObject, rb_intern("Math"));
    return nm::RubyObject(rb_funcall(cMath, rb_intern("sqrt"), 1, obj.rval));
  }
}

#endif // RUBY_OBJECT_H
