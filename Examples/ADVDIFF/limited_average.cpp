#include "limited_average.h"

using namespace INMOST;

//shortcuts
typedef Storage::real       real;
typedef Storage::integer    integer;
typedef Storage::enumerator enumerator;

/// Function that performs non-linear weighted convex combination of two fluxes in
/// a style of non-linear two-point flux approximation.
/// Non-linear two-point flux approximation has positivity-preserving property.
/// Depending on the template type it will account or not for derivatives in
/// convex combination, which may result in faster convergence but loss of 
/// positivity-preserving property. 
/// @param fluxK Flux for back cell.
/// @param fluxL Flux for front cell.
/// @param fluxK_cK Flux for back cell only part with back cell.
/// @param fluxK_cL Flux for back cell only part with front cell.
/// @param fluxL_cK Flux for front cell only part with back cell.
/// @param fluxL_cL Flux for front cell only part with front cell.
/// @param outK Flux for the back cell.
/// @param outL FLux for the front cell. May be different from outK.
/// @param regularization Regularization parameter in computation of expression.
template<typename VarType> 
void LimitedAverageTwoPoint(const INMOST::variable & fluxK,
	                        const INMOST::variable & fluxL,
							const INMOST::variable & fluxK_cK,
							const INMOST::variable & fluxL_cL,
							INMOST::variable & outK, 
							INMOST::variable & outL,
							real regularization)
{
	VarType muK, muL, div; //convex combination that elemenates additional entries
	assign(muL, -(fluxK - fluxK_cK));
	assign(muK, +(fluxL - fluxL_cL));
	// diffusion and advection fluxes come with a different sign,
	// this trick should regularize both. Regularization is needed
	// in order to avoid tottaly zero flux with zero derivatives.
	if( muK < 0.0 || muL < 0.0 )
	{
		muK -= regularization;
		muL -= regularization;
	}
	else
	{
		muK += regularization;
		muL += regularization;
	}
	div = muK+muL;
	muK /= div;
	muL /= div;
	// two-point transmissibilities for each cell
	outK = outL = (fluxK_cK*muK) + (fluxL_cL*muL);
}

/// Function performs non-linear weighted convex combination for two fluxes with Van-Albada limited averages.
/// Van-Albada limited averages function is supposed to enforce TVD property. 
/// It does not lead to always positive approximation which is required to satisfy discrete maximum principle on each iteration.
/// This funciton returns single flux definition, thus it is conservative on each iteration.
/// Depending on the template type it will account or not for derivatives in convex combination.
/// @param inK Flux at one side of the face.
/// @param inL Flux at another side of the face.
/// @param outK Limited conservative flux for one side of the face.
/// @param outL Limited conservative flux for another side of the face.
/// @param regularization Regularization parameter in computation of expression.
template<typename VarType> 
void LimitedAverageVanAlbada(const variable & inK, const variable & inL, variable & outK, variable & outL, real regularization)
{
	VarType vK, vL; //values on each side
	VarType muK, muL; //multipliers for each side
	VarType denomenator;
	assign(vK,inK);
	assign(vL,inL);
	muK = vL*vL + regularization;
	muL = vK*vK + regularization;
	//compute nonlinear coefficients
	denomenator = muK + muL;
	//multipliers for each side
	muK /= denomenator;
	muL /= denomenator;
	//output
	outK = outL = muK*inK + muL*inL;
}


/// Function performs non-linear weighted convex combination for two fluxes with coefficients in power of 4.
/// This limited averages function is supposed to enforce TVD property. 
/// It does not lead to always positive approximation which is required to satisfy discrete maximum principle on each iteration.
/// This funciton returns single flux definition, thus it is conservative on each iteration.
/// Depending on the template type it will account or not for derivatives in convex combination.
/// @param inK Flux at one side of the face.
/// @param inL Flux at another side of the face.
/// @param outK Limited conservative flux for one side of the face.
/// @param outL Limited conservative flux for another side of the face.
/// @param regularization Regularization parameter in computation of expression.
template<typename VarType> 
void LimitedAverageQuadratic(const variable & inK, const variable & inL, variable & outK, variable & outL, real regularization)
{
	VarType vK, vL; //values on each side
	VarType muK, muL; //multipliers for each side
	VarType denomenator;
	assign(vK,inK);
	assign(vL,inL);
	muK = vL*vL*vL*vL + regularization;
	muL = vK*vK*vK*vK + regularization;
	//compute nonlinear coefficients
	denomenator = muK + muL;
	//multipliers for each side
	muK /= denomenator;
	muL /= denomenator;
	//output
	outK = outL = muK*inK + muL*inL;
}

/// Function performs non-linear weighted convex combination for two fluxes with Van-Albada limited averages.
/// Van-Albada limited averages function is supposed to enforce TVD property. 
/// It does not lead to always positive approximation which is required to satisfy discrete maximum principle on each iteration.
/// This funciton returns dual flux definition, thus it is not conservative on each iteration.
/// No large difference from LimitedAverageVanAlbada if we account for all derivative.
/// @param inK Flux at one side of the face.
/// @param inL Flux at another side of the face.
/// @param outK Limited conservative flux for one side of the face.
/// @param outL Limited conservative flux for another side of the face.
/// @param regularization Regularization parameter in computation of expression.
void LimitedAverageVanAlbadaDual(const variable & inK, const variable & inL, variable & outK, variable & outL, real regularization)
{
	real vK, vL; //values on each side
	real cK, cL; //multipliers for each side
	real denomenator;
	assign(vK,inK);
	assign(vL,inL);
	//compute nonlinear coefficients
	denomenator = vK*vK + vL*vL + 2*regularization;
	//multipliers for each side
	cK = (vL*(vK + vL) + regularization*2)/denomenator;
	cL = (vK*(vK + vL) + regularization*2)/denomenator;
	//output
	outK = cK*inK;
	outL = cL*inL;
}

/// Function performs non-linear weighted convex combination for two fluxes with Van-Albada limited averages with correction.
/// Corrected Van-Albada limited averages function is supposed to enforce discrete maximum principle property on convergence. 
/// This funciton returns single flux definition, thus it is conservative on each iteration.
/// Depending on the template type it will account or not for derivatives in convex combination.
/// @param inK Flux at one side of the face.
/// @param inL Flux at another side of the face.
/// @param outK Limited conservative flux for one side of the face.
/// @param outL Limited conservative flux for another side of the face.
/// @param regularization Regularization parameter in computation of expression.
template<typename VarType> 
void LimitedAverageVanAlbadaCorrected(const variable & inK, const variable & inL, variable & outK, variable & outL, real regularization)
{
	VarType vK, vL; //values on each side
	VarType muK, muL; //multipliers for each side
	VarType denomenator;
	VarType sign; //sign of the product of fluxes
	VarType corr; //correction
	assign(vK,inK);
	assign(vL,inL);
	sign = soft_sign(vK*vL,regularization);
	//sign = 2.0/(1+exp(-vK*vL))-1.0;
	corr = 0.5*(1.0-sign);
	muK = vL*vL - corr*vK*vL + regularization;
	muL = vK*vK - corr*vK*vL + regularization;
	//compute nonlinear coefficients
	denomenator = muK + muL;
	//multipliers for each side
	muK /= denomenator;
	muL /= denomenator;
	//output
	outK = outL = muK*inK + muL*inL;
}


/// Function performs non-linear weighted convex combination for two fluxes with coefficients in power of 4and correction.
/// Corrected quadratic limited averages function is supposed to enforce discrete maximum principle property on convergence. 
/// This funciton returns single flux definition, thus it is conservative on each iteration.
/// Depending on the template type it will account or not for derivatives in convex combination.
/// @param inK Flux at one side of the face.
/// @param inL Flux at another side of the face.
/// @param outK Limited conservative flux for one side of the face.
/// @param outL Limited conservative flux for another side of the face.
/// @param regularization Regularization parameter in computation of expression.
template<typename VarType> 
void LimitedAverageQuadraticCorrected(const variable & inK, const variable & inL, variable & outK, variable & outL, real regularization)
{
	VarType vK, vL; //values on each side
	VarType muK, muL; //multipliers for each side
	VarType denomenator;
	VarType sign; //sign of the product of fluxes
	VarType corr; //correction
	assign(vK,inK);
	assign(vL,inL);
	sign = soft_sign(vK*vL,regularization);
	//sign = 2.0/(1+exp(-vK*vL))-1.0;
	corr = 0.5*(1.0-sign);
	muK = vL*vL*vL*vL - corr*vL*vK*vK*vK + regularization;
	muL = vK*vK*vK*vK - corr*vK*vL*vL*vL + regularization;
	//compute nonlinear coefficients
	denomenator = muK + muL;
	//multipliers for each side
	muK /= denomenator;
	muL /= denomenator;
	//output
	outK = outL = muK*inK + muL*inL;
}

/// Function performs non-linear weighted convex combination for two fluxes with Van-Albada limited averages 
/// with correction and dual flux definition.
/// Not conservative on each iteration but should satisfy discrete maximum principle on each iteration (up to regularization).
/// @param inL Flux at another side of the face.
/// @param outK Limited conservative flux for one side of the face.
/// @param outL Limited conservative flux for another side of the face.
/// @param regularization Regularization parameter in computation of expression.
void LimitedAverageVanAlbadaCorrectedDual(const variable & inK, const variable & inL, variable & outK, variable & outL, real regularization)
{
	real vK, vL; //values on each side
	real cK, cL; //multipliers for each side
	real sign, corr, denomenator;
	vK = get_value(inK);
	vL = get_value(inL);
	//positivity corrector
	sign = soft_sign(vK*vL,regularization);
	corr = 0.5*(1.0 - sign);
	//compute nonlinear coefficients
	denomenator = vK*vK + vL*vL - 2*corr*vK*vL + 2*regularization;
	//multipliers for each side
	cK = (vL*(vK + vL)*(1-corr) + regularization*2)/denomenator;
	cL = (vK*(vK + vL)*(1-corr) + regularization*2)/denomenator;
	//output
	outK = cK*inK;
	outL = cL*inL;
}


/// Function performs non-linear weighted convex combination for two fluxes with Van-Leer limited averages.
/// Van-Leer limited averages function is supposed to enforce discrete maximum principle property. 
/// This funciton returns single flux definition, thus it is conservative on each iteration.
/// Depending on the template type it will account or not for derivatives in
/// convex combination.
/// This method will satsify discrete maximum principle only on convergence, but not on each
/// nonlinear iteration. To satisfy discrete maximum principle on each iteration see dual flux 
/// Van-Leer limited averages function.
/// @param inK Flux at one side of the face.
/// @param inL Flux at another side of the face.
/// @param outK Limited conservative flux for one side of the face.
/// @param outL Limited conservative flux for another side of the face.
/// @param regularization Regularization parameter in computation of expression.
template<typename VarType> 
void LimitedAverageVanLeer(const variable & inK, const variable & inL, variable & outK, variable & outL, real regularization)
{
	VarType muK, muL; //coefficients of convex combination
	VarType denomenator;
	assign(muK,inL);
	assign(muL,inK);
	muK = soft_fabs(muK,regularization);
	muL = soft_fabs(muL,regularization);
	denomenator = muK + muL;
	//multipliers for each side
	muK /= denomenator;
	muL /= denomenator;
	//output
	outK = outL = muK*inK + muL*inL;
}

/// Function performs non-linear weighted convex combination for two fluxes with Van-Leer limited averages
/// with dual flux definition. It is not conservative but satisfy discrete maximum principle on each iteration.
/// @param inK Flux at one side of the face.
/// @param inL Flux at another side of the face.
/// @param outK Limited conservative flux for one side of the face.
/// @param outL Limited conservative flux for another side of the face.
/// @param regularization Regularization parameter in computation of expression.
void LimitedAverageVanLeerDual(const variable & inK, const variable & inL, variable & outK, variable & outL, real regularization)
{
	//currently no account for variations
	real vK, vL; //values on each side
	real absK, absL; //absolute values of input
	real muK, muL; //coefficients of convex combination
	real cK, cL; //multipliers for each side
	vK = get_value(inK);
	vL = get_value(inL);
	if( vK*vL > 0.0 )
	{
		//calculate absolute values for nonlinear weighting
		absK = soft_fabs(vK,regularization);
		absL = soft_fabs(vL,regularization);
		//compute nonlinear coefficients
		muK = absL/(absK+absL);
		muL = absK/(absK+absL);
		//multipliers for each side
		cK = 2*muK;
		cL = 2*muL;
		//output
		outK = cK*inK;
		outL = cL*inL;
	}
	else outK = outL = 0.0;
}

/// Multi-point flux approximation method.
/// Computes wighted combination based on provided coefficients cK and cL.
/// @param cK Some coefficient for the left flux that will form convex combination. Can be just 1.
/// @param cL Some coefficient for the right flux that will form convex combination. Can be just 1.
/// @param inK Flux at one side of the face.
/// @param inL Flux at another side of the face.
/// @param outK Limited conservative flux for one side of the face.
/// @param outL Limited conservative flux for another side of the face.
/// @param regularization Regularization parameter in computation of expression.
void LimitedAverageMultiPoint(real cK, real cL, const variable & inK, const variable & inL, variable & outK, variable & outL, real regularization)
{
	real muK, muL;
	//take abs just in case
	cK = fabs(cK)+regularization;
	cL = fabs(cL)+regularization;
	muK = cL/(cK+cL);
	muL = cK/(cK+cL);
	outK = outL = muK*inK+muL*inL;
}


void LimitedAverageNonlinearMultiPoint(const INMOST::variable & corrK, 
									   const INMOST::variable & corrL, 
									   INMOST::variable & outK, 
									   INMOST::variable & outL,
									   INMOST::Storage::real regularization,
									   SchemeType scheme_type)
{
	if( scheme_type == NMPFA_VAN_ALBADA )
		LimitedAverageVanAlbada<variable>(corrK,corrL,outK,outL,regularization);
	else if( scheme_type == NMPFA_VAN_ALBADA_PICARD )
		LimitedAverageVanAlbada<real>(corrK,corrL,outK,outL,regularization);
	else if( scheme_type == NMPFA_VAN_ALBADA_DUALFLUX )
		LimitedAverageVanAlbadaDual(corrK,corrL,outK,outL,regularization);
	else if( scheme_type == NMPFA_VAN_ALBADA_CORRECTED )
		LimitedAverageVanAlbadaCorrected<variable>(corrK,corrL,outK,outL,regularization);
	else if( scheme_type == NMPFA_VAN_ALBADA_CORRECTED_PICARD )
		LimitedAverageVanAlbadaCorrected<real>(corrK,corrL,outK,outL,regularization);
	else if( scheme_type == NMPFA_VAN_ALBADA_CORRECTED_DUALFLUX )
		LimitedAverageVanAlbadaCorrectedDual(corrK,corrL,outK,outL,regularization);
	else if( scheme_type == NMPFA_VAN_LEER )
		LimitedAverageVanLeer<variable>(corrK,corrL,outK,outL,regularization);
	else if( scheme_type == NMPFA_VAN_LEER_PICARD )
		LimitedAverageVanLeer<real>(corrK,corrL,outK,outL,regularization);
	else if( scheme_type == NMPFA_VAN_LEER_DUALFLUX )
		LimitedAverageVanLeerDual(corrK,corrL,outK,outL,regularization);
	else if( scheme_type == NMPFA_QUADRATIC )
		LimitedAverageQuadratic<variable>(corrK,corrL,outK,outL,regularization);
	else if( scheme_type == NMPFA_QUADRATIC_CORRECTED )
		LimitedAverageQuadraticCorrected<variable>(corrK,corrL,outK,outL,regularization);
	else if( scheme_type == MPFA )
		LimitedAverageMultiPoint(1,1,corrK,corrL,outK,outL,regularization);
}

void LimitedAverageNonlinearTwoPoint(const INMOST::variable & fluxK,
									 const INMOST::variable & fluxL,
									 const INMOST::variable & fluxK_cK,
									 const INMOST::variable & fluxL_cL,
									 INMOST::variable & outK, 
									 INMOST::variable & outL,
									 INMOST::Storage::real regularization,
									 SchemeType scheme_type)
{
	if( scheme_type == NTPFA_PICARD )
		LimitedAverageTwoPoint<real>(fluxK,fluxL,fluxK_cK,fluxL_cL,outK,outL,regularization);
	else if( scheme_type == NTPFA )
		LimitedAverageTwoPoint<variable>(fluxK,fluxL,fluxK_cK,fluxL_cL,outK,outL,regularization);
}





