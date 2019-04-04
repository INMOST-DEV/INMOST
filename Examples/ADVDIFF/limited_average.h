#ifndef __LIMITED_AVERAGE_H
#define __LIMITED_AVERAGE_H

#include "inmost.h"

/// Specification of scheme type.
enum SchemeType
{
	MPFA, //< Linear scheme with algebraic mean.
	NTPFA, //< With account for derivatives in convex combination. Fast but can violate positivity.
	NTPFA_PICARD, //< Do not account for derivatives in convex combination. Always positive with good right hand side.
	NMPFA_VAN_LEER, //< Satisfy discrete maximum principle on convergence, non-Lipschitz-continuous before regularization. With account for derivatives in convex combination.
	NMPFA_VAN_LEER_PICARD, //< Satisfy discrete maximum principle on convergence, non-Lipschitz-continuous before regularization. No account for derivatives in convex combination.
	NMPFA_VAN_LEER_DUALFLUX, //< Satisfy discrete maximum principle on each iteration, non-Lipschitz-continuous before regularization. No account for derivatives in convex combination.
	NMPFA_VAN_ALBADA, //< Should satisfy TVD, continuous, single flux definition. With account for derivatives in convex combination.
	NMPFA_VAN_ALBADA_PICARD, //< Should satisfy TVD, continuous, single flux definition. No account for derivatives in convex combination.
	NMPFA_VAN_ALBADA_DUALFLUX, //< Should satisfy TVD, continuous, dual flux definition. No account for derivatives in convex combination.
	NMPFA_VAN_ALBADA_CORRECTED, //< Should satisfy discrete maximum principle, non-Lipschitz-continuous before regularization. With account for derivatives in convex combination.
	NMPFA_VAN_ALBADA_CORRECTED_PICARD, //< Should satisfy discrete maximum principle, non-Lipschitz-continuous before regularization. No account for derivatives in convex combination.
	NMPFA_VAN_ALBADA_CORRECTED_DUALFLUX, //< Should satisfy discrete maximum principle, non-Lipschitz-continuous before regularization. No account for derivatives in convex combination.
	NMPFA_QUADRATIC, //< Should satisfy TVD
	NMPFA_QUADRATIC_CORRECTED //< Should satisfy discrete maximum principle
};


/// Limited averaging function, computes and returns the fluxes inside.
/// @param fluxK Flux for back cell.
/// @param fluxL Flux for front cell.
/// @param outK Flux for the back cell.
/// @param outL FLux for the front cell. May be different from outK.
/// @param regularization Regularization parameter in computation of expression
/// @param scheme_type Method for the flux computation, see SchemeType.
void LimitedAverageNonlinearMultiPoint(const INMOST::variable & fluxK, 
									   const INMOST::variable & fluxL,
									   INMOST::variable & outK, 
									   INMOST::variable & outL,
									   INMOST::Storage::real regularization,
									   SchemeType scheme_type);

/// Function that performs non-linear weighted convex combination of two fluxes in
/// a style of non-linear two-point flux approximation.
/// Non-linear two-point flux approximation has positivity-preserving property.
/// Depending on the template type it will account or not for derivatives in
/// convex combination, which may result in faster convergence but loss of 
/// positivity-preserving property. 
/// @param fluxK Flux for back cell.
/// @param fluxL Flux for front cell.
/// @param fluxK_cK Flux for back cell only with back cell.
/// @param fluxK_cL Flux for back cell only with front cell.
/// @param fluxL_cK Flux for front cell only with back cell.
/// @param fluxL_cL Flux for front cell only with front cell.
/// @param outK Flux for the back cell.
/// @param outL FLux for the front cell. May be different from outK.
/// @param regularization Regularization parameter in computation of expression.
void LimitedAverageNonlinearTwoPoint(const INMOST::variable & fluxK,
									 const INMOST::variable & fluxL,
									 const INMOST::variable & fluxK_cK,
									 const INMOST::variable & fluxL_cL,
									 INMOST::variable & outK, 
									 INMOST::variable & outL,
									 INMOST::Storage::real regularization,
									 SchemeType scheme_type);


#endif