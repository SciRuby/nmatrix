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
// == rational.h
//
// Functions and classes for dealing with rational numbers.

#ifndef RATIONAL_H
#define RATIONAL_H

/*
 * Standard Includes
 */

#include <type_traits>
#include <ruby.h>
#include <iostream>

/*
 * Project Includes
 */

#include "types.h"

#include "util/util.h"

/*
 * Macros
 */

/*
 * Types
 */

namespace nm {

template <typename Type> class Rational;

typedef Rational<int16_t>	Rational32;
typedef Rational<int32_t>	Rational64;
typedef Rational<int64_t>	Rational128;

/*
 * Data
 */

/*
 * Classes and Functions
 */

template <typename Type>
class Rational {
	public:
	// The numerator and denominator of the rational number.
	Type n;
	Type d;
	
	/*
	 * Default constructor.
	 */
	inline Rational(Type num = 0, Type den = 1) : n(num), d(den) {}
	
	/*
	 * Copy constructors.
	 */
	template <typename OtherType>
	inline Rational(const Rational<OtherType>& other) : n(other.n), d(other.d) {}

	template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
	inline Rational(const Complex<FloatType>& other) : n(0), d(1) {
	  rb_raise(rb_eNotImpError, "cannot convert from complex to rational");
	}

	Rational(const RubyObject& other);

  /*
   * Rational inverse function -- creates a copy, but inverted.
   */
  inline Rational<Type> inverse() const {
    return Rational<Type>(this->d, this->n);
  }

	/*
	 * Binary operator definitions for varous types.
	 */

	//////////////////////////////////
	// Rational-Rational Operations //
	//////////////////////////////////
	
	template <typename OtherType>
	inline Rational<Type> operator+(const Rational<OtherType>& other) const {
		Rational<Type> result((this->n * other.d) + (other.n * this->d), this->d * other.d);
		
		long simplify = gcf<Type>(result.n, result.d);

		result.n /= simplify;
		result.d /= simplify;
		
		return result;
	}

	template <typename OtherType>
	inline Rational<Type>& operator+=(const Rational<OtherType>& other) {
    this->n = (this->n * other.d) + (other.n * this->d);
    this->d = this->d * other.d;

    long simplify = gcf<Type>(this->n, this->d);

    this->n /= simplify;
    this->d /= simplify;

    return *this;
	}
	
	template <typename OtherType>
	inline Rational<Type> operator-(const Rational<OtherType>& other) const {
		Rational<Type> result((this->n * other.d) - (other.n * this->d), this->d * other.d);
		
		long simplify = gcf<Type>(result.n, result.d);
		
		result.n /= simplify;
		result.d /= simplify;
		
		return result;
	}

	template <typename OtherType>
	inline Rational<Type>& operator-=(const Rational<OtherType>& other) {
    this->n = (this->n * other.d) - (other.n * this->d);
    this->d = this->d * other.d;

    long simplify = gcf<Type>(this->n, this->d);

    this->n /= simplify;
    this->d /= simplify;

    return *this;
	}
	
	template <typename OtherType>
	inline Rational<Type> operator*(const Rational<OtherType>& other) const {
		int g1 = gcf<Type>(this->n, other.d);
		int g2 = gcf<Type>(this->d, other.n);
		
		return Rational<Type>((this->n / g1) * (other.n / g2), (this->d / g2) * (other.d / g1));
	}


	template <typename OtherType>
	inline Rational<Type>& operator*=(const Rational<OtherType>& other) {
    int g1 = gcf<Type>(this->n, other.d);
    int g2 = gcf<Type>(this->d, other.n);

    this->n = (this->n / g1) * (other.n / g2);
    this->d = (this->d / g2) * (other.d / g1);

    return *this;
	}

	
	template <typename OtherType>
	inline Rational<Type> operator/(const Rational<OtherType>& other) const {
		return *this * Rational<OtherType>(other.d, other.n);
	}

	template <typename OtherType>
	inline Rational<Type> operator/=(const Rational<OtherType>& other) {
	  *this *= Rational<OtherType>(other.d, other.n);
	  return *this;
	}
	
	template <typename OtherType>
	inline Rational<Type> operator%(const Rational<OtherType>& other) const {
		long floor_div      = (this->n * other.n) / (this->d * other.d);
		Rational<Type> prod	= other * Rational<long>(floor_div, 1);
		
		return Rational<long>(this->n, other.n) - prod;
	}
	
	template <typename OtherType>
	inline bool operator<(const Rational<OtherType>& other) const {
		return (this->n * other.d) < (other.n * this->d);
	}
	
	template <typename OtherType>
	inline bool operator>(const Rational<OtherType>& other) const {
		return (this->n * other.d) > (other.n * this->d);
	}
	
	template <typename OtherType>
	inline bool operator==(const Rational<OtherType>& other) const {
		return (this->n == other.n) && (this->d == other.d);
	}
	
	template <typename OtherType>
	inline bool operator!=(const Rational<OtherType>& other) const {
		return !(*this == other);
	}
	
	template <typename OtherType>
	inline bool operator<=(const Rational<OtherType>& other) const {
		return (*this < other) || (*this == other);
	}
	
	template <typename OtherType>
	inline bool operator>=(const Rational<OtherType>& other) const {
		return (*this > other) || (*this == other);
	}
	
	template <typename OtherType>
	inline operator Rational<OtherType> () const {
		return Rational<OtherType>((OtherType)this->n, (OtherType)this->d);
	}
	
	////////////////////////////////
	// Rational-Native Operations //
	////////////////////////////////
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline Rational<Type> operator+(const IntType& other) const {
		return *this + Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline Rational<Type> operator-(const IntType& other) const {
		return *this - Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline Rational<Type> operator*(const IntType& other) const {
		return *this * Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline Rational<Type> operator/(const IntType& other) const {
		return *this / Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline Rational<Type> operator%(const IntType& other) const {
		return *this % Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator<(const IntType& other) const {
		return *this < Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator>(const IntType& other) const {
		return *this > Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator==(const IntType& other) const {
		return *this == Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator!=(const IntType& other) const {
		return *this != Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator<=(const IntType& other) const {
		return *this <= Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator>=(const IntType& other) const {
		return *this >= Rational<Type>(other);
	}
	
	template <typename NumType, typename = typename std::enable_if<std::is_arithmetic<NumType>::value>::type>
	inline operator NumType () const {
		return (NumType)this->n / (NumType)this->d;
	}
	
	/*
	 * Special casting operator for Complex numbers.
	 */
	template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
	inline operator Rational<FloatType> () const {
		return Rational<FloatType>(((FloatType)this->n) / ((FloatType)this->d));
	}
};

// Negative operator
template <typename Type, typename = typename std::enable_if<std::is_integral<Type>::value>::type>
inline Rational<Type> operator-(const Rational<Type>& rhs) {
  return Rational<Type>(-rhs.n, rhs.d);
}

////////////////////////////////
// Native-Rational Operations //
////////////////////////////////

/*
 * Integer Math
 */

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline Rational<RationalType> operator+(const IntType& left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) + right;
}

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline Rational<RationalType> operator-(const IntType& left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) - right;
}

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline Rational<RationalType> operator*(const IntType& left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) * right;
}

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline Rational<RationalType> operator/(const IntType& left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) / right;
}

/*
 * Floating Point Math
 */

template <typename FloatType, typename RationalType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline FloatType operator+(const FloatType& left, const Rational<RationalType>& right) {
	return left + (FloatType)right;
}

template <typename FloatType, typename RationalType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline FloatType operator-(const FloatType& left, const Rational<RationalType>& right) {
	return left - (FloatType)right;
}

template <typename FloatType, typename RationalType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline FloatType operator*(const FloatType& left, const Rational<RationalType>& right) {
	return left * (FloatType)right;
}

template <typename FloatType, typename RationalType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline FloatType operator/(const FloatType& left, const Rational<RationalType>& right) {
	return left / (FloatType)right;
}

/*
 * Comparisons
 */

template <typename NativeType, typename RationalType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator<(const NativeType left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) < right;
	return (left * right.d) < right.n;
}

template <typename NativeType, typename RationalType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator>(const NativeType left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) > right;
	return (left * right.d) > right.n;
}

template <typename IntType, typename RationalType>
inline bool operator==(const typename std::enable_if<std::is_integral<IntType>::value, IntType>::type left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) == right;
	return (left * right.d) == right.n;
}

template <typename FloatType, typename RationalType>
inline bool operator==(const typename std::enable_if<std::is_floating_point<FloatType>::value, FloatType>::type left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) == right;
	return FP_EQUAL(left, ((FloatType)right));
}

template <typename NativeType, typename RationalType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator!=(const NativeType left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) != right;
	return !(left == right);
}

template <typename NativeType, typename RationalType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator<=(const NativeType left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) <= right;
	return (left < right) or (left == right);
}

template <typename NativeType, typename RationalType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator>=(const NativeType left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) >= right;
	return (left > right) or (left == right); 
}

template <typename Type>
inline std::ostream& operator<<(std::ostream& out, const Rational<Type>& rhs) {
  out << rhs.n << "/" << rhs.d << std::flush;
  return out;
}

} // end of namespace nm

namespace std {
  template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
  nm::Rational<IntType> abs(const nm::Rational<IntType>& value) {
    if (value.n >= 0) return value;
    return nm::Rational<IntType>(-value.n, value.d);
  }

  template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
  nm::Rational<IntType> sqrt(const nm::Rational<IntType>& value) {
    nm::Rational<IntType> result(std::sqrt(value.n), std::sqrt(value.d));
    if (value * value == result)      return result;
    else                              rb_raise(rb_eArgError, "square root of the given rational is not rational");
  }
}

#endif // RATIONAL_H
