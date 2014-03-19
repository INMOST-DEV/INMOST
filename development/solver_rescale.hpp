#pragma once
#ifndef __SOLVER_RESCALE__
#define __SOLVER_RESCALE__

//TODO:
//test implementation
//implement RSS scaling from Documents/Read/solver/scaling/2007KaporinRJNAMM.pdf

#include "inmost_solver.h"
#include "solver_prototypes.hpp"

using namespace INMOST;

class RAS_rescale : public Method //Row-Alternating Scaling for 2-norm
{
	Solver::Matrix * Alink;
	Solver::Vector * xlink;
	Solver::Vector * blink;
	Solver::OrderInfo * info;
	Solver::Vector DL,DR;
	INMOST_DATA_ENUM_TYPE iters;
	INMOST_DATA_REAL_TYPE subst, eps;
	bool init;
public:

	RAS_rescale(Solver::Matrix & A, Solver::Vector & x, Solver::Vector & b,  Solver::OrderInfo & Minfo)
		:Alink(&A),xlink(&x),blink(&b),info(&Minfo)
	{
		iters = 5;
		eps = 1.0e-54;
		subst = 1.0;
	}
	void Copy(const Method * other)
	{
		const RAS_rescale * b = dynamic_cast<const RAS_rescale *>(other);
		assert(b != NULL);
		Alink = b->Alink;
		info = b->info;
		DL = b->DL;
		DR = b->DR;
		iters = b->iters;
		eps = b->eps;
		subst = b->subst;
	}
	RAS_rescale(const RAS_rescale & other) 
	{
		Copy(&other);
	}
	RAS_rescale & operator =(RAS_rescale const & other)
	{
		Copy(&other);
		return *this;
	}
	INMOST_DATA_ENUM_TYPE & EnumParameter(std::string name)
	{
		if( name == "iters" )
			return iters;
		throw -1;
	}
	INMOST_DATA_REAL_TYPE & RealParameter(std::string name)
	{
		if( name == "substitute" )
			return subst;
		else if( name == "epsilon" )
			return eps;
		throw -1;
	}
	bool Initialize()
	{
		if (isInitialized()) Finalize();
		INMOST_DATA_ENUM_TYPE mobeg,moend, vlocbeg, vlocend, vbeg, vend,k, end, iter,r;
		info->GetOverlapRegion(info->GetRank(),mobeg,moend);
		info->GetLocalRegion(info->GetRank(),vlocbeg,vlocend);
		info->GetVectorRegion(vbeg,vend);
		DL.SetInterval(mobeg,moend);
		info->PrepareVector(DR);
		for(k = mobeg; k < moend; k++)
		{
			for(Solver::Row::iterator q = (*Alink)[k].Begin(); q != (*Alink)[k].End(); ++q)
				DL[k] += q->second*q->second;
			if( DL[k] < eps ) DL[k] = 1.0/subst; else DL[k] = 1.0/DL[k];
		}
		for(iter = 0; iter < iters; iter++)
		{
			std::fill(DR.Begin(),DR.End(),0.0);
			for(k = vlocbeg; k < vlocend; k++)
				for(Solver::Row::iterator q = (*Alink)[k].Begin(); q != (*Alink)[k].End(); ++q) DR[q->first] += DL[k]*q->second*q->second;
			info->Accumulate(DR);
			info->Update(DR);
			for(k = vbeg; k < vend; k++)
				if( DR[k] < eps ) DR[k] = 1.0/subst; else DR[k] = 1.0/DR[k];
			std::fill(DL.Begin(),DL.End(),0.0);
			for(k = mobeg; k < moend; k++)
				for(Solver::Row::iterator q = (*Alink)[k].Begin(); q != (*Alink)[k].End(); ++q) DL[k] += DR[q->first]*q->second*q->second;
			for(k = mobeg; k < moend; k++)
				if( DL[k] < eps ) DL[k] = 1.0/subst; else DL[k] = 1.0/DL[k];
		}
		for(k = mobeg; k < moend; k++) DL[k] = sqrt(DL[k]);
		for(k =  vbeg; k <  vend; k++) DR[k] = sqrt(DR[k]);
		//Rescale matrix
		for(k = mobeg; k < moend; k++)
		{
			for(Solver::Row::iterator q = (*Alink)[k].Begin(); q != A[k].End(); ++q)
				q->second = DL[k] * q->second * DR[q->first];
		}
		for(k = mobeg; k < moend; k++) (*blink)[k] *= DL[k];
		info->Update(*blink);
		for(k = vbeg; k < vend; k++) (*xlink)[k] /= DR[k];
		init = true;
		return true;
	}
	bool isInitialized() { return init; }
	bool ReplaceMAT(Solver::Matrix & A) { if (isInitialized()) Finalize();  Alink = &A; return true; }
	bool ReplaceRHS(Solver::Vector & b)
	{
		INMOST_DATA_ENUM_TYPE mobeg,moend,k;
		info->GetOverlapRegion(info->GetRank(),mobeg,moend);
		//restore old rhs
		for(k = mobeg; k < moend; k++) (*blink)[k] /= DL[k];
		//select new rhs
		blink = &b;
		//rescale new rhs
		for(k = mobeg; k < moend; k++) (*blink)[k] *= DL[k];
		info->Update(*blink);
		return true;
	}
	bool ReplaceSOL(Solver::Vector & x)
	{
		INMOST_DATA_ENUM_TYPE vbeg,vend,k;
		info->GetVectorRegion(vbeg,vend);
		//restore old sol
		for(k = vbeg; k < vend; k++) (*xlink)[k] *= DR[k];
		//select new sol
		xlink = &x;
		//rescale new sol
		for(k = vbeg; k < vend; k++) (*xlink)[k] /= DR[k];
		return true;
	}
	bool Solve(Solver::Vector & x, Solver::Vector & out)
	{
		// A x = b
		// DL A DR DR^{-1} x = DL b
		//rescale vector x by DR
		INMOST_DATA_ENUM_TYPE vbeg,vend;
		info->GetVectorRegion(vbeg,vend);
		for(INMOST_DATA_ENUM_TYPE k = vbeg; k != vend; k++) out[k] = x[k] * DR[k];
		return true;
	}
	bool Finalize()
	{
		if (!isFinalized())
		{
			INMOST_DATA_ENUM_TYPE mobeg, moend, k, end, r, vbeg, vend;
			info->GetOverlapRegion(info->GetRank(), mobeg, moend);
			info->GetVectorRegion(vbeg, vend);
			for (k = mobeg; k < moend; k++)
			{
				Solver::Row & Ak = (*Alink)[k];
				end = Ak.Size();
				for (r = 0; r < end; r++)
					Ak.GetValue(r) = Ak.GetValue(r) / DL[k] / DR[Ak.GetIndex(r)];
			}
			for (k = mobeg; k < moend; k++) (*blink)[k] /= DL[k];
			info->Update(*blink);
			for (k = vbeg; k < vend; k++) (*xlink)[k] *= DR[k];
			init = false;
		}
		return true;
	}
	bool isFinalized() { return !init; }
	~RAS_rescale()
	{
		if (!isFinalized()) Finalize();
	}
	Method * Duplicate() { return new RAS_rescale(*this); }
};

#endif