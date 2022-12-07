#ifndef INMOST_OPERATOR_INCLUDED
#define INMOST_OPERATOR_INCLUDED

#include "inmost_variable.h"

#if defined(USE_AUTODIFF) && defined(USE_MESH) && defined(USE_SOLVER)

namespace INMOST
{
	class Model;
	/// This class is responsible to unite access to various point-wise
	/// implementations of discrete operators, such as grad, curl.
	/// \todo
	/// 1. Different types of operators:
	///    time-stepping,
	///    local point-wise (curl,grad on element),
	///    global integrators (div,curl on domain),
	///    interpolators,
	///    inter-mesh interpolators.
	///    Each has it's own functions. Implementation should be
	///    flexible enough to prevent limitation.
	/// 2. Ultimately operators should stack together:
	///    for staggered incompressible navier-stokes:
	/// Time(nU) + Projection(Divergence(ConvectionDiffusion(nU,\mu,Reconstruction(nU)))) - Grad(P) = f
	/// Divergence(nU) = 0
	///
	class AbstractOperator
	{
	public:
		/// Destroy all the data of the operator.
		virtual ~AbstractOperator() {};
		/// Initialize all the data necessery to evalute the operator.
		virtual bool Initialize(Model & m) = 0;
		/// Let operator prepare data on the mesh before evaluation.
		virtual bool PrepareIterations() = 0;
		/// Setup coupling with unknowns of otheer models
		virtual bool SetupCoupling(Model& P) { return true; }
		/// Check, whether we need to compute operator on this element.
		virtual bool isValid(const Storage & e) const = 0;
		/// Provides input domain of the operator. (TODO)
		virtual std::pair<ElementType,MarkerType> GetUnknownDomain() const = 0;
		/// Provides output domain of the operator. (TODO)
		virtual std::pair<ElementType, MarkerType> GetOperatorDomain() const = 0;
		/// Compute expression of the opertor (TODO)
		virtual vMatrix Evaluate(const Storage& e) const = 0;
	};
}
		
#endif //USE_AUTODIFF && USE_MESH

#endif //INMOST_OPERATOR_INCLUDED
