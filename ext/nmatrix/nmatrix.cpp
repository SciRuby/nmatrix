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
// == nmatrix.cpp
//
// Main C++ source file for NMatrix. Contains Init_nmatrix and most Ruby
// instance and class methods for NMatrix. Also responsible for calling Init
// methods on related modules.

/*
 * Standard Includes
 */

extern "C" {
#if defined HAVE_CBLAS_H
  #include <cblas.h>
#elif defined HAVE_ATLAS_CBLAS_H
  #include <atlas/cblas.h>
#endif

#if defined HAVE_CLAPACK_H
  #include <clapack.h>
#elif defined HAVE_ATLAS_CLAPACK_H
  #include <atlas/clapack.h>
#endif
}

#include <ruby.h>
#include <algorithm> // std::min
#include <fstream>

/*
 * Project Includes
 */
#include "nmatrix_config.h"

#include "types.h"
#include "data/data.h"
#include "math/math.h"
#include "util/io.h"
#include "storage/storage.h"
#include "storage/list/list.h"
#include "storage/yale/yale.h"

#include "nmatrix.h"

#include "ruby_constants.h"

/*
 * Ruby internals
 */


/*
 * Macros
 */


/*
 * Global Variables
 */

namespace nm {

  /*
   * This function is pulled out separately so it can be called for hermitian matrix writing, which also uses it.
   */
  template <typename DType>
  size_t write_padded_dense_elements_upper(std::ofstream& f, DENSE_STORAGE* storage, symm_t symm) {
    // Write upper triangular portion. Assume 2D square matrix.
    DType* elements = reinterpret_cast<DType*>(storage->elements);
    size_t length = storage->shape[0];

    size_t bytes_written = 0;

    for (size_t i = 0; i < length; ++i) { // which row are we on?

      f.write( reinterpret_cast<const char*>( &(elements[ i*(length + 1) ]) ),
               (length - i) * sizeof(DType) );

      bytes_written += (length - i) * sizeof(DType);
    }
    return bytes_written;
  }

  /*
   * We need to specialize for Hermitian matrices. The next six functions accomplish that specialization, basically
   * by ensuring that non-complex matrices cannot read or write hermitians (which would cause big problems).
   */
  template <typename DType>
  size_t write_padded_dense_elements_herm(std::ofstream& f, DENSE_STORAGE* storage, symm_t symm) {
    rb_raise(rb_eArgError, "cannot write a non-complex matrix as hermitian");
  }

  template <>
  size_t write_padded_dense_elements_herm<Complex64>(std::ofstream& f, DENSE_STORAGE* storage, symm_t symm) {
    return write_padded_dense_elements_upper<Complex64>(f, storage, symm);
  }

  template <>
  size_t write_padded_dense_elements_herm<Complex128>(std::ofstream& f, DENSE_STORAGE* storage, symm_t symm) {
    return write_padded_dense_elements_upper<Complex128>(f, storage, symm);
  }

  template <typename DType>
  void read_padded_dense_elements_herm(DType* elements, size_t length) {
    rb_raise(rb_eArgError, "cannot read a non-complex matrix as hermitian");
  }

  template <>
  void read_padded_dense_elements_herm(Complex64* elements, size_t length) {
    for (size_t i = 0; i < length; ++i) {
      for (size_t j = i+1; j < length; ++j) {
        elements[j * length + i] = elements[i * length + j].conjugate();
      }
    }
  }

  template <>
  void read_padded_dense_elements_herm(Complex128* elements, size_t length) {
    for (size_t i = 0; i < length; ++i) {
      for (size_t j = i+1; j < length; ++j) {
        elements[j * length + i] = elements[i * length + j].conjugate();
      }
    }
  }

  /*
   * Read the elements of a dense storage matrix from a binary file, padded to 64-bits.
   *
   * storage should already be allocated. No initialization necessary.
   */
  template <typename DType>
  void read_padded_dense_elements(std::ifstream& f, DENSE_STORAGE* storage, nm::symm_t symm) {
    size_t bytes_read = 0;

    if (symm == nm::NONSYMM) {
      // Easy. Simply read the whole elements array.
      size_t length = nm_storage_count_max_elements(reinterpret_cast<STORAGE*>(storage));
      f.read(reinterpret_cast<char*>(storage->elements), length * sizeof(DType) );

      bytes_read += length * sizeof(DType);
    } else if (symm == LOWER) {

      // Read lower triangular portion and initialize remainder to 0
      DType* elements = reinterpret_cast<DType*>(storage->elements);
      size_t length = storage->shape[0];

      for (size_t i = 0; i < length; ++i) { // which row?

        f.read( reinterpret_cast<char*>(&(elements[i * length])), (i + 1) * sizeof(DType) );

        // need to zero-fill the rest of the row.
        for (size_t j = i+1; j < length; ++j)
          elements[i * length + j] = 0;

        bytes_read += (i + 1) * sizeof(DType);
      }
    } else {

      DType* elements = reinterpret_cast<DType*>(storage->elements);
      size_t length = storage->shape[0];

      for (size_t i = 0; i < length; ++i) { // which row?
        f.read( reinterpret_cast<char*>(&(elements[i * (length + 1)])), (length - i) * sizeof(DType) );

        bytes_read += (length - i) * sizeof(DType);
      }

      if (symm == SYMM) {
        for (size_t i = 0; i < length; ++i) {
          for (size_t j = i+1; j < length; ++j) {
            elements[j * length + i] = elements[i * length + j];
          }
        }
      } else if (symm == SKEW) {
        for (size_t i = 0; i < length; ++i) {
          for (size_t j = i+1; j < length; ++j) {
            elements[j * length + i] = -elements[i * length + j];
          }
        }
      } else if (symm == HERM) {
        read_padded_dense_elements_herm<DType>(elements, length);

      } else if (symm == UPPER) { // zero-fill the rest of the rows
        for (size_t i = 0; i < length; ++i) {
          for(size_t j = i+1; j < length; ++j) {
            elements[j * length + i] = 0;
          }
        }
      }

    }

    // Ignore any padding.
    if (bytes_read % 8) f.ignore(bytes_read % 8);
  }

  template <typename DType>
  void write_padded_yale_elements(std::ofstream& f, YALE_STORAGE* storage, size_t length, nm::symm_t symm) {
    if (symm != nm::NONSYMM) rb_raise(rb_eNotImpError, "Yale matrices can only be read/written in full form");

    // Keep track of bytes written for each of A and IJA so we know how much padding to use.
    size_t bytes_written = length * sizeof(DType);

    // Write A array
    f.write(reinterpret_cast<const char*>(storage->a), bytes_written);

    // Padding
    int64_t zero = 0;
    f.write(reinterpret_cast<const char*>(&zero), bytes_written % 8);

    bytes_written = length * sizeof(IType);
    f.write(reinterpret_cast<const char*>(storage->ija), bytes_written);

    // More padding
    f.write(reinterpret_cast<const char*>(&zero), bytes_written % 8);
  }


  template <typename DType>
  void read_padded_yale_elements(std::ifstream& f, YALE_STORAGE* storage, size_t length, nm::symm_t symm) {
    if (symm != NONSYMM) rb_raise(rb_eNotImpError, "Yale matrices can only be read/written in full form");

    size_t bytes_read = length * sizeof(DType);
    f.read(reinterpret_cast<char*>(storage->a), bytes_read);

    int64_t padding = 0;
    f.read(reinterpret_cast<char*>(&padding), bytes_read % 8);

    bytes_read = length * sizeof(IType);
    f.read(reinterpret_cast<char*>(storage->ija), bytes_read);

    f.read(reinterpret_cast<char*>(&padding), bytes_read % 8);
  }

  /*
   * Write the elements of a dense storage matrix to a binary file, padded to 64-bits.
   */
  template <typename DType>
  void write_padded_dense_elements(std::ofstream& f, DENSE_STORAGE* storage, nm::symm_t symm) {
    size_t bytes_written = 0;

    if (symm == nm::NONSYMM) {
      // Simply write the whole elements array.
      size_t length = nm_storage_count_max_elements(storage);
      f.write(reinterpret_cast<const char*>(storage->elements), length * sizeof(DType));

      bytes_written += length * sizeof(DType);

    } else if (symm == nm::LOWER) {

      // Write lower triangular portion. Assume 2D square matrix.
      DType* elements = reinterpret_cast<DType*>(storage->elements);
      size_t length = storage->shape[0];
      for (size_t i = 0; i < length; ++i) { // which row?

        f.write( reinterpret_cast<const char*>( &(elements[i * length]) ),
                 (i + 1) * sizeof(DType) );

        bytes_written += (i + 1) * sizeof(DType);
      }
    } else if (symm == nm::HERM) {
      bytes_written += write_padded_dense_elements_herm<DType>(f, storage, symm);
    } else { // HERM, UPPER, SYMM, SKEW
      bytes_written += write_padded_dense_elements_upper<DType>(f, storage, symm);
    }

    // Padding
    int64_t zero = 0;
    f.write(reinterpret_cast<const char*>(&zero), bytes_written % 8);
  }

} // end of namespace nm

extern "C" {
  #include "ruby_nmatrix.c"
} // end of extern "C"
