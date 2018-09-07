#ifndef INMOST_NONLINEAR_INCLUDED
#define INMOST_NONLINEAR_INCLUDED

#include "inmost_common.h"
#include "inmost_autodiff.h"

#if defined(USE_NONLINEAR)
namespace INMOST
{
	
	typedef INMOST_DATA_BULK_TYPE RequestedAction;
	
	static const RequestedAction COMPUTE_FUNCTION = 0x01;
	static const RequestedAction COMPUTE_JACOBIAN = 0x02;
	static const RequestedAction COMPUTE_HESSIAN  = 0x04;
	static const RequestedAction FINISHED         = 0x08;
	
	class NonlinearSolver
	{
		//Automatizator & aut;
	public:
		//NonlinearSolver(Automatizator & aut) : aut(aut) {}
		NonlinearSolver() {}
		NonlinearSolver(const NonlinearSolver & b) /*: aut(b.aut)*/ {}
		NonlinearSolver & operator =(NonlinearSolver const & b) {/*aut = b.aut;*/ return *this;}
		~NonlinearSolver() {}
		
		RequestedAction GetAction() const;
		INMOST_DATA_REAL_TYPE GetResidual() const;
		INMOST_DATA_ENUM_TYPE GetIterations() const;
		std::string GetReason() const;
	};
	
};
#endif

#endif //INMOST_NONLINEAR_INCLUDED
