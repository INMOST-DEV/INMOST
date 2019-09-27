
#ifndef __SOLVER_MTILUC2__
#define __SOLVER_MTILUC2__
#include <iomanip>

#include "inmost_solver.h"
#include "../solver_prototypes.hpp"

class MTILUC_preconditioner : public Method
{
	typedef struct Interval_t
	{
		INMOST_DATA_ENUM_TYPE first, last;
		INMOST_DATA_ENUM_TYPE Size() { return last - first; }
	} Interval;
	typedef std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> coord;
	typedef std::pair<INMOST_DATA_REAL_TYPE, coord > wgt_coord;
	typedef std::vector< wgt_coord > wgt_coords;
	typedef dynarray<INMOST_DATA_ENUM_TYPE,256> levels_t;
	
	//Sparse::Vector div;
	array<Sparse::Row::entry> LU_Entries;
	std::vector<Sparse::Row::entry> B_Entries;
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LU_Diag;
	interval<INMOST_DATA_ENUM_TYPE, Interval> U_Address, L_Address, B_Address;
	
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> temp; // temporal place for solve phase
	//reordering information
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > ddP,ddQ;
	INMOST_DATA_REAL_TYPE condestL, condestU;	
	INMOST_DATA_ENUM_TYPE estimator;
	INMOST_DATA_REAL_TYPE iluc2_tau;
	INMOST_DATA_REAL_TYPE tau, eps;
	INMOST_DATA_ENUM_TYPE sciters,verbosity;
	Sparse::Matrix * Alink;
	Solver::OrderInfo * info;
	bool init;
	void DumpMatrix(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
									array<Sparse::Row::entry> & Entries,
									INMOST_DATA_ENUM_TYPE wmbeg, INMOST_DATA_ENUM_TYPE wmend,
									std::string file_name);
	void CheckOrder(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
					 array<Sparse::Row::entry> & Entries, 
					 INMOST_DATA_ENUM_TYPE rbeg, INMOST_DATA_ENUM_TYPE rend);
	void SwapEntries(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
					 array<Sparse::Row::entry> & Entries, 
					 INMOST_DATA_ENUM_TYPE rbeg, INMOST_DATA_ENUM_TYPE rend, 
					 INMOST_DATA_ENUM_TYPE k, INMOST_DATA_ENUM_TYPE j);
	void SwapLine(interval<INMOST_DATA_ENUM_TYPE, Interval> & Line, INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j);

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
public:
	INMOST_DATA_ENUM_TYPE & EnumParameter(std::string name);
	INMOST_DATA_REAL_TYPE & RealParameter(std::string name);
	void Copy(const Method * other);
	MTILUC_preconditioner(const MTILUC_preconditioner & other);
	MTILUC_preconditioner & operator =(MTILUC_preconditioner const & other);
	MTILUC_preconditioner(Solver::OrderInfo & info);
	bool isInitialized();
	bool isFinalized();
	bool Initialize();
	bool Finalize();
	void ApplyB(double alpha, Sparse::Vector & x, double beta, Sparse::Vector & y);
	int Descend(Sparse::Vector & input, Sparse::Vector & output);
	int Ascend(Sparse::Vector & input, Sparse::Vector & output);
	bool Solve(Sparse::Vector & input, Sparse::Vector & output);
	bool ReplaceMAT(Sparse::Matrix & A);
	bool ReplaceSOL(Sparse::Vector & x);
	bool ReplaceRHS(Sparse::Vector & b);
	Method * Duplicate();
	~MTILUC_preconditioner();
};

#endif //__SOLVER_MTILUC2__
