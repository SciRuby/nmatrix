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
// == complex.h
//
// Functions and classes for dealing with complex numbers.

#ifndef COMPLEX_H
#define COMPLEX_H

/*
 * Standard Includes
 */

#include <type_traits>
#include <iostream>
#include <cmath>

/*
 * Project Includes
 */

#include "types.h"

/*
 * Macros
 */

/*
 * Types
 */
namespace nm {

class RubyObject;
template <typename IntType> class Rational;
template <typename Type> class Complex;

typedef Complex<float32_t> Complex64;
typedef Complex<float64_t> Complex128;

/*
 * Data
 */
 
/*
 * Classes and Functions
 */

template <typename Type>
class Complex {
	public:
	// The real and immaginary parts of the complex number.
	Type r;
	Type i;
	
	/*
	 * Default constructor.
	 */
	inline Complex(Type real = 0, Type imaginary = 0) : r(real), i(imaginary) {}
	
	/*
	 * Copy constructors.
	 */
	template <typename ComplexType>
	inline Complex(const Complex<ComplexType>& other) : r(other.r), i(other.i) {}

	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline Complex(const Rational<IntType>& other) : r(Type(other.n) / Type(other.d)), i(0) {}

  Complex(const RubyObject& other);

  /*
   * Complex conjugate function -- creates a copy, but inverted.
   */
  inline Complex<Type> conjugate() const {
    return Complex<Type>(this->r, -(this->i));
  }

  /*
   * Complex inverse function -- creates a copy, but inverted.
   *
   * FIXME: Check that this doesn't duplicate functionality of NativeType / Complex<Type>
   */
  inline Complex<Type> inverse() const {
    Complex<Type> conj = conjugate();
    Type denom = this->r * this->r + this->i * this->i;
    return Complex<Type>(conj.r / denom, conj.i / denom);
  }



	/*
	 * Binary operator definitions for various types.
	 */
	
	////////////////////////////////
	// Complex-Complex Operations //
	////////////////////////////////
	
	template <typename OtherType>
	inline Complex<Type> operator+(const Complex<OtherType>& other) const {
		return Complex<Type>(this->r + other.r, this->i + other.i);
	}

  template <typename OtherType>
  inline Complex<Type>& operator+=(const Complex<OtherType>& other) {
    this->r += other.r;
    this->i += other.i;
    return *this;
  }

	template <typename OtherType>
	inline Complex<Type> operator-(const Complex<OtherType>& other) const {
		return Complex<Type>(this->r - other.r, this->i - other.i);
	}
	
	template <typename OtherType>
	inline Complex<Type> operator*(const Complex<OtherType>& other) const {
		return Complex<Type>(this->r * other.r - this->i * other.i, this->r * other.i + this->i * other.r);
	}

  template <typename OtherType>
  inline Complex<Type>& operator*=(const Complex<OtherType>& other) {
    this->r = this->r * other.r - this->i * other.i;
    this->i = this->r * other.i + this->i * other.r;
    return *this;
  }

	template <typename OtherType>
	inline Complex<Type> operator/(const Complex<OtherType>& other) const {
		Type new_r, new_i;
		Type denom = other.i * other.i + other.r * other.r;
		
		new_r = (this->r * other.r + this->i * other.i) / denom;
		new_i = (this->i * other.r - this->r * other.i) / denom;
		
		return Complex<Type>(new_r, new_i);
	}
	
	template <typename OtherType>
	inline bool operator<(const Complex<OtherType>& other) const {
		return (this->r < other.r) || ((this->r <= other.r) && (this->i < other.i));
	}
	
	template <typename OtherType>
	inline bool operator>(const Complex<OtherType>& other) const {
		return (this->r > other.r) || ((this->r >= other.r) && (this->i > other.i));
	}
	
	template <typename OtherType>
	inline bool operator==(const Complex<OtherType>& other) const {
		return FP_EQUAL(this->r, other.r) && FP_EQUAL(this->i, other.i);
	}
	
	template <typename OtherType>
	inline bool operator!=(const Complex<OtherType>& other) const {
		return !(*this == other);
	}
	
	template <typename OtherType>
	inline bool operator<=(const Complex<OtherType>& other) const {
		return (*this < other) || (*this == other);
	}
	
	template <typename OtherType>
	inline bool operator>=(const Complex<OtherType>& other) const {
		return (*this > other) || (*this == other);
	}
	
	template <typename OtherType>
	inline operator Complex<OtherType> () const {
		return Complex<OtherType>((OtherType)this->r, (OtherType)this->i);
	}

	/////////////////////////////////
	// Complex-Rational Operations //
	/////////////////////////////////
	
	template <typename RationalType>
	inline Complex<Type> operator+(const Rational<RationalType>& other) const {
		return *this + Complex<Type>(other);
	}
	
	template <typename RationalType>
	inline Complex<Type> operator-(const Rational<RationalType>& other) const {
		return *this - Complex<Type>(other);
	}
	
	template <typename RationalType>
	inline Complex<Type> operator*(const Rational<RationalType>& other) const {
		return *this * Complex<Type>(other);
	}
	
	template <typename RationalType>
	inline Complex<Type> operator/(const Rational<RationalType>& other) const {
		return *this / Complex<Type>(other);
	}
	
	template <typename RationalType, typename = typename std::enable_if<std::is_integral<RationalType>::value>::type>
	inline bool operator!=(const Rational<RationalType>& other) const {
	  return *this != Complex<Type>(other);
	}

  template <typename RationalType, typename = typename std::enable_if<std::is_integral<RationalType>::value>::type>
	inline bool operator==(const Rational<RationalType>& other) const {
	  return *this == Complex<Type>(other);
	}

	///////////////////////////////
	// Complex-Native Operations //
	///////////////////////////////
	
	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline Complex<Type> operator+(const NativeType& other) const {
		return *this + Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline Complex<Type> operator-(const NativeType& other) const {
		return *this - Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline Complex<Type> operator*(const NativeType& other) const {
		return *this * Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline Complex<Type> operator/(const NativeType& other) const {
		return *this / Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator<(const NativeType& other) const {
		return *this < Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator>(const NativeType& other) const {
		return *this > Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator==(const NativeType& other) const {
		return *this == Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator!=(const NativeType& other) const {
		return *this != Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator<=(const NativeType& other) const {
		return *this <= Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator>=(const NativeType& other) const {
		return *this >= Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline operator NativeType () const {
		return (NativeType)this->r;
	}
};


/////////////////////////////////
// Rational-Complex Operations //
/////////////////////////////////


template <typename IntType, typename ComplexType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator==(const Rational<IntType>& left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) == right;
}

template <typename IntType, typename ComplexType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator!=(const Rational<IntType>& left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) != right;
}


///////////////////////////////
// Native-Complex Operations //
///////////////////////////////

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline Complex<ComplexType> operator+(const NativeType& left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) + right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline Complex<ComplexType> operator-(const NativeType& left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) - right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline Complex<ComplexType> operator*(const NativeType& left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) * right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline Complex<ComplexType> operator/(const NativeType& left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) / right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator<(const NativeType left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) < right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator>(const NativeType left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) > right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator==(const NativeType left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) == right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator!=(const NativeType left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) != right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator<=(const NativeType left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) <= right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator>=(const NativeType left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) >= right;
}

template <typename Type>
inline std::ostream& operator<<(std::ostream& out, const Complex<Type>& rhs) {
  out << "(" << rhs.r << "," << rhs.i << "i)" << std::flush;
  return out;
}

// Negative operator
template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline Complex<IntType> operator-(const Complex<IntType>& rhs) {
  return Complex<IntType>(-rhs.r, -rhs.i);
}

} // end of namespace nm

namespace std {
  template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
  nm::Complex<FloatType> piecewise_abs(const nm::Complex<FloatType>& value) {
    return nm::Complex<FloatType>(value.r < 0 ? -value.r : value.r,
                                  value.i < 0 ? -value.i : value.i);
  }

  template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
  nm::Complex<FloatType> real_abs(const nm::Complex<FloatType>& value) {
    return nm::Complex<FloatType>(value.r < 0 ? -value.r : value.r,
                                  value.i);
  }

  template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
  nm::Complex<FloatType> imag_abs(const nm::Complex<FloatType>& value) {
    return nm::Complex<FloatType>(value.r,
                                  value.i < 0 ? -value.i : value.i);
  }

  template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
  double abs(const nm::Complex<FloatType>& value) {
    return std::sqrt(double(value.r)*double(value.r) + double(value.i)*double(value.i));
  }
}

#endif // COMPLEX_H
