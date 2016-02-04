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
// == cblas_templaces_lapacke.h
//
// Define template functions for calling CBLAS functions in the
// nm::math::lapacke namespace.
//

#ifndef CBLAS_TEMPLATES_LAPACK_H
#define CBLAS_TEMPLATES_LAPACK_H

//includes so we have access to internal implementations
#include "math/rotg.h"
#include "math/rot.h"
#include "math/asum.h"
#include "math/nrm2.h"
#include "math/imax.h"
#include "math/scal.h"
#include "math/gemv.h"
#include "math/gemm.h"
#include "math/trsm.h"

namespace nm { namespace math { namespace lapacke {
 
//Add cblas templates in the correct namespace
#include "math/cblas_templates_core.h"

}}}

#endif
