#ifndef __CHECKS_H
#define __CHECKS_H

#include "inmost.h"

/// Check that this will give M-matrix for i-th unknown.
/// @param i Position of the diagonal unknown.
/// @param flux Computed flux expression with derivatives.
/// @param name Name of the flux to be printed out.
/// @param print Print out errors.
bool check_flux_properties(INMOST_DATA_ENUM_TYPE i, const INMOST::variable & flux, const std::string & name, bool print = true);


///Check monotonicity of the matrix.
/// @param A Input matrix to be checked.
/// @param print Print out errors.
bool check_matrix_properties(const INMOST::Sparse::Matrix & A, bool print = true);


#endif
