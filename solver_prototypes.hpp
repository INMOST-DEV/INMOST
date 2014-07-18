#pragma once
#ifndef __SOLVER_PROTOTYPES__
#define __SOLVER_PROTOTYPES__

//TODO:
// allow user to provide new RHS or new initial solution

#include "inmost_solver.h"

using namespace INMOST;

class Method
{
public:
	virtual INMOST_DATA_REAL_TYPE & RealParameter(std::string name) = 0;
	virtual INMOST_DATA_ENUM_TYPE & EnumParameter(std::string name) = 0;
	virtual bool Initialize() = 0;
	virtual bool isInitialized() = 0;
	virtual bool Finalize() = 0;
	virtual bool isFinalized() = 0;
	virtual bool Solve(Solver::Vector & input, Solver::Vector & output) = 0;
	virtual bool ReplaceMAT(Solver::Matrix & A) = 0; //provide matrix
	virtual bool ReplaceRHS(Solver::Vector & b) = 0; //apply modification such as rescaling or reordering to the new right hand side
	virtual bool ReplaceSOL(Solver::Vector & x) = 0; //apply modification such as rescaling or reordering to the new solution
	virtual void Copy(const Method * other) = 0;
	virtual Method * Duplicate() {return NULL;}
	virtual ~Method() {}
};

class IterativeMethod : public Method
{
public:
	virtual INMOST_DATA_ENUM_TYPE GetIterations() = 0;
	virtual INMOST_DATA_REAL_TYPE GetResidual() = 0;
	virtual std::string GetReason() = 0;
	virtual ~IterativeMethod() {}
};

#endif
