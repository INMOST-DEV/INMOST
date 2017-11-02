
#ifndef __SOLVER_MLMTILUC2__
#define __SOLVER_MLMTILUC2__
#include <iomanip>

#include "inmost_solver.h"
#include "../solver_prototypes.hpp"

class MLMTILUC_preconditioner : public Method
{
	typedef struct Interval_t
	{
		INMOST_DATA_ENUM_TYPE first, last;
		INMOST_DATA_ENUM_TYPE Size() { return last - first; }
	} Interval;
	typedef dynarray<INMOST_DATA_ENUM_TYPE,256> levels_t;
	
	std::vector<Sparse::Row::entry> LU_Entries, B_Entries;
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LU_Diag;
	interval<INMOST_DATA_ENUM_TYPE, Interval> U_Address, L_Address, B_Address;
	std::vector<interval<INMOST_DATA_ENUM_TYPE, Interval> *> F_Address, E_Address;
	std::vector<Sparse::Row::entry> E_Entries, F_Entries;
	levels_t level_size; //remember size of each level
	std::vector<Interval> level_interval;
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> temp; // temporal place for solve phase
	//reordering information
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > ddP,ddQ;
	INMOST_DATA_REAL_TYPE condestL, condestU;	
	INMOST_DATA_ENUM_TYPE estimator;
	INMOST_DATA_REAL_TYPE iluc2_tau;
	INMOST_DATA_REAL_TYPE pivot_cond;
	INMOST_DATA_REAL_TYPE pivot_diag;
	INMOST_DATA_REAL_TYPE tau, eps;
	INMOST_DATA_ENUM_TYPE sciters;
	Sparse::Matrix * Alink;
	Solver::OrderInfo * info;
	bool init;
	void DumpMatrix(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
					std::vector<Sparse::Row::entry> & Entries,
					INMOST_DATA_ENUM_TYPE wmbeg, INMOST_DATA_ENUM_TYPE wmend,
					std::string file_name);
	void CheckOrder(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
					std::vector<Sparse::Row::entry> & Entries,
					INMOST_DATA_ENUM_TYPE rbeg, INMOST_DATA_ENUM_TYPE rend);
	void inversePQ(INMOST_DATA_ENUM_TYPE wbeg,
				   INMOST_DATA_ENUM_TYPE wend, 
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ);

	void applyPQ(INMOST_DATA_ENUM_TYPE wbeg,
				 INMOST_DATA_ENUM_TYPE wend,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ);
	void ReorderEF(INMOST_DATA_ENUM_TYPE wbeg,
				   INMOST_DATA_ENUM_TYPE wend,
				   interval<INMOST_DATA_ENUM_TYPE, bool> & donePQ,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ);
public:
	INMOST_DATA_ENUM_TYPE & EnumParameter(std::string name);
	INMOST_DATA_REAL_TYPE & RealParameter(std::string name);
	void Copy(const Method * other);
	MLMTILUC_preconditioner(const MLMTILUC_preconditioner & other);
	MLMTILUC_preconditioner & operator =(MLMTILUC_preconditioner const & other);
	MLMTILUC_preconditioner(Solver::OrderInfo & info);
	bool isInitialized();
	bool isFinalized();
	bool Initialize();
	bool Finalize();
	void ApplyB(double alpha, Sparse::Vector & x, double beta, Sparse::Vector & y);
	int Descend(int level, Sparse::Vector & inout);
	int Ascend(int level, Sparse::Vector & inout);
	bool Solve(Sparse::Vector & input, Sparse::Vector & output);
	bool ReplaceMAT(Sparse::Matrix & A);
	bool ReplaceSOL(Sparse::Vector & x);
	bool ReplaceRHS(Sparse::Vector & b);
	Method * Duplicate();
	~MLMTILUC_preconditioner();
};

#endif //__SOLVER_MLMTILUC2__
