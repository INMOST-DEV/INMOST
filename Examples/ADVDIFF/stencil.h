#ifndef __STENCIL_H
#define __STENCIL_H

#include "inmost.h"

enum StencilType
{
	CONVECTION,
	DIFFUSION
};

///Structure that helps to input information for stencil computation.
struct Stencil
{
	INMOST::Storage::real             v[3];     //< direction
	INMOST::Storage::real_array       coefs;    //< stencil coefficients
	INMOST::Storage::reference_array  elems;    //< stencil elements
	INMOST::Storage::real *           rhs;      //< right hand side
	INMOST::Storage::real             sign;     //< triplet sign
	StencilType                       type;     //< type of requested stencil
	bool                              computed; //< indicates whether the stencil was found
	bool                              nonnegative; //< indicates that all the coefficients are positive
};

/// Find triplet in current cell, currently uses nearest neighbours.
///TODO: allow for consideration of multiple directions of gradient.
/// @param cK Cell in which to seek the gradient.
/// @param compute A set of stencils to be computed.
/// @param tag_BC Tag for boundary conditions.
/// @param tag_K Tag for diffusion tensor.
/// @param tag_iT Tag for interpolation tensor (must be private to process in OpenMP).
/// @param tag_iC Tag for interpolation correction (must be private to process in OpenMP).
/// @param tag_U Tag for velocity value.
/// @param boundary_marker Marker for boundary faces.
/// @param bridge Type of elements to be used as bridge to cells in layers.
/// @param regularization Regularizer for degenerate diffusion.
/// @param max_layers Maximum number of layers to be considered.
bool find_stencils(INMOST::Cell cK, 
	               std::vector<Stencil> & compute, 
				   INMOST::Tag tag_BC, 
				   INMOST::Tag tag_K, 
				   INMOST::Tag tag_iT, 
				   INMOST::Tag tag_iC,
				   //INMOST::Tag tag_U,
				   INMOST::MarkerType boundary_marker,
				   INMOST::ElementType bridge,
				   INMOST::Storage::real regularization,
				   int max_layers);

#endif //__STENCIL_H