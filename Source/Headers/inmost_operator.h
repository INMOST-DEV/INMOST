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
		/// Operator may introduce it's unknowns on this step
		virtual bool PrepareEntries(Model & m) = 0;
		/// Initialize all the data necessery to evalute the operator.
		virtual bool Initialize(Model & m) = 0;
		/// Let operator prepare data on the mesh before evaluation.
		virtual bool Prepare() = 0;
		/// This allows operator to update introduced variables.
		/// May be useful for mimetic finite differences operators.
		virtual bool UpdateSolution(const Sparse::Vector & V, double alpha) = 0;
		/// This allows operator to update time-dependent variables.
		virtual bool UpdateTimeStep() = 0;
		/// Check, whether we need to compute operator on this element.
		virtual bool isValid(const Storage & e) const = 0;
		/// Provide data necessery for initialization of the operator.
		virtual bool SetTag(std::string name, Tag t) = 0;
		/// Provide data necessery for initialization of the operator.
		virtual bool SetMarker(std::string name, MarkerType marker) = 0;
	};
	
	class AbstractSpaceOperator : virtual public AbstractOperator
	{
	public:
		virtual ~AbstractSpaceOperator() {}
		/// Rebuild the discretization of the operator due to update of boundary conditions.
		//virtual bool UpdateBoundaryCondition() = 0;
	};
	
	struct EquationPosition
	{
		const AbstractEntry & entry;
		INMOST_DATA_ENUM_TYPE pos;
		EquationPosition(const AbstractEntry & entry, INMOST_DATA_ENUM_TYPE pos) : entry(entry), pos(pos) {}
		INMOST_DATA_ENUM_TYPE GetPosition(const Storage & e) const {return entry.Index(e,pos);}
	};
	
	class AbstractScalarOperator : virtual public AbstractSpaceOperator
	{
	public:
		virtual ~AbstractScalarOperator() {}
		/// The main function provided by operator.
		/// @param R Residual that is filled in.
		/// @param e Element on which operator is evaluated.
		/// @param param Expression to be evaluated by operator.
		/// @param eqpos Equation position in residual that operator affects.
		/// @return Success or failure of operator application on provied expression.
		virtual bool Evaluate(Residual & R, const Storage & e, const abstract_variable & param, const EquationPosition & eqpos) const = 0;
	};
	
	class AbstractVectorOperator : virtual public AbstractSpaceOperator
	{
	public:
		virtual ~AbstractVectorOperator() {}
		/// The main function provided by operator.
		/// @param expr Operator may have multiple expression types, each is accessable by number.
		/// @param param Expression to be evaluated by operator.
		/// @return Result of operator application on provied expression in the form of 3x1 matrix.
		virtual void Evaluate(Residual & R, const Storage & e, const stored_variable_expression param[3], const EquationPosition eqpos[3]) const = 0;
	};
	
	class AbstractTensorOperator : virtual public AbstractSpaceOperator
	{
	public:
		virtual ~AbstractTensorOperator() {}
		/// The main function provided by operator.
		/// @param expr Operator may have multiple expression types, each is accessable by number.
		/// @param param Expression to be evaluated by operator.
		/// @return Result of operator application on provied expression in the form of 3x3 matrix.
		virtual bool Evaluate(Residual & R, const Storage & e, const stored_variable_expression param[9], const EquationPosition eqpos[9]) const = 0;
	};
	
	class AbstractVoigtOperator : virtual public AbstractSpaceOperator
	{
	public:
		virtual ~AbstractVoigtOperator() {}
		/// The main function provided by operator.
		/// @param expr Operator may have multiple expression types, each is accessable by number.
		/// @param param Expression to be evaluated by operator.
		/// @return Result of operator application on provied expression in the form of 6x1 matrix.
		virtual bool Evaluate(Residual & R, const Storage & e, const stored_variable_expression param[6], const EquationPosition eqpos[6]) const = 0;
	};
}
		
#endif //USE_AUTODIFF && USE_MESH

#endif //INMOST_OPERATOR_INCLUDED
