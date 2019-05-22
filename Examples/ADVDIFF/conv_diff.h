#ifndef __CONV_DIFF_H
#define __CONV_DIFF_H

#include "inmost.h"
#include "limited_average.h"

///Comment me
class ConvectionDiffusion
{
	INMOST::Mesh * m;
	//tags for convection
	INMOST::Tag tag_CONV_CB;	// Coefficients for upwind corrector at back cell
	INMOST::Tag tag_CONV_EB;	// Elements for upwind corrector at back cell
	INMOST::Tag tag_CONV_RB;	// Right hand side for upwind corrector at back cell
	INMOST::Tag tag_CONV_CF;	// Coefficients for upwind corrector at front cell
	INMOST::Tag tag_CONV_EF;	// Elements for upwind corrector at front cell
	INMOST::Tag tag_CONV_RF;	// Right hand side for upwind corrector at front cell
	INMOST::Tag tag_CONV_CU;	// Upwind coefficient
	INMOST::Tag tag_CONV_RU;	// Right hand side for upwind
	INMOST::Tag tag_CONV_EU;	// Upstream cell
	INMOST::Tag tag_CONV_F;		// A flag indicating is there
								// a full nonlinear part  (internal     = 2),
								// a one-sided correction (on outlet BC = 1),
								// no correction          (on inlet BC  = 0)
	INMOST::Tag tag_CONV_VU;	// Two vectors for upwind correctors at back and front cells
	//tags for diffusion
	INMOST::Tag tag_DIFF_CB;	// Coefficients for transversal correction at back cell
	INMOST::Tag tag_DIFF_EB;	// Elements for transversal correction at back cell
	INMOST::Tag tag_DIFF_RB;	// Right hand side for two transversal correction at back cell
	INMOST::Tag tag_DIFF_CF;	// Coefficients for transversal correction at front cell
	INMOST::Tag tag_DIFF_EF;	// Elements for transversal correction at front cell
	INMOST::Tag tag_DIFF_RF;	// Right hand side for two transversal correction at back cell
	INMOST::Tag tag_DIFF_CT;	// Two-point transmissibility
	INMOST::Tag tag_DIFF_RT;	// Two-point right hand side
	INMOST::Tag tag_DIFF_F;		// A flag indicating is there
								// a full nonlinear part   (internal = 3),
								// no nonlinear part       (internal = 2),
								// a one-sided correction  (on BC    = 1),
								// no one-sided correction (on BC    = 0)
	INMOST::Tag tag_DIFF_VT;	// Vector for transversal correction to two-point part


	INMOST::Tag tag_U;  // Normal velocity vector on faces
	INMOST::Tag tag_K;  // Diffusion tensor
	INMOST::Tag tag_BC; // Boundary conditions
	

	INMOST::MarkerType build_flux; // defines wheather to build flux approximation on interface in parallel
	INMOST::MarkerType boundary_face; //defines boundary faces in parallel

	bool initialization_failure; //there is a failure during intialization

	bool perform_correction_diffusion; //add nonlinear correction to two-point flux
	bool perform_correction_convection; //add nonlinear correction to single point upstream 
public:
	ConvectionDiffusion(INMOST::Mesh * _m, INMOST::Tag _tag_U, INMOST::Tag _tag_K, INMOST::Tag _tag_BC, INMOST::MarkerType _boundary_face, bool correct_convection, bool correct_diffusion);
	bool Failure() const;
	~ConvectionDiffusion();
	void DiffusionFluxes(INMOST::abstract_variable & param, INMOST::Face fKL, INMOST::variable & fluxT, INMOST::variable & corrK, INMOST::variable & corrL, INMOST::variable & corrK_cK, INMOST::variable & corrL_cL) const;
	void AdvectionFluxes(INMOST::abstract_variable & param, INMOST::Face fKL, INMOST::variable & fluxT, INMOST::variable & corrK, INMOST::variable & corrL, INMOST::variable & corrK_cK, INMOST::variable & corrL_cL) const;
	void Averaging(SchemeType scheme_type, INMOST_DATA_REAL_TYPE regularization, const INMOST::variable & fluxT, const INMOST::variable & corrK, const INMOST::variable & corrL, const INMOST::variable & corrK_cK, const INMOST::variable & corrL_cL, INMOST::variable & fluxK, INMOST::variable & fluxL, bool boundary) const;
	
	bool BuildFlux(INMOST::Face f) const;
};

#endif //__CONV_DIFF_H