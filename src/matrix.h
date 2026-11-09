// src/matrix.h
//
// Helpers shared by the C++ sources.

#ifndef FDATEST_MATRIX_H
#define FDATEST_MATRIX_H

#include <cpp11.hpp>

// An nrow x ncol double matrix. Rf_allocMatrix() computes the length as
// R_xlen_t, while cpp11's matrix constructor multiplies nrow * ncol in int,
// which overflows beyond INT_MAX elements.
inline cpp11::sexp new_double_matrix(int nrow, int ncol) {
  return cpp11::safe[Rf_allocMatrix](REALSXP, nrow, ncol);
}

#endif
