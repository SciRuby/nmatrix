//need to define a few things before including the real lapacke.h

#include "data/data.h" //needed because this is where our complex types are defined

//tell LAPACKE to use our complex types
#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_float nm::Complex64
#define lapack_complex_double nm::Complex128

//define name-mangling scheme for FORTRAN functions
//ADD_ means that the symbol dgemm_ is associated with the fortran
//function DGEMM
#define ADD_

//now we can include the real lapacke.h
#include "lapacke.h"
