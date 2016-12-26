
#ifndef __SOLVER_DDPQILUC2__
#define __SOLVER_DDPQILUC2__
#include <iomanip>

#include "inmost_solver.h"
#include "../solver_prototypes.hpp"

class ILUC_preconditioner : public Method
{
	typedef struct Interval_t
	{
		INMOST_DATA_ENUM_TYPE first, last;
		INMOST_DATA_ENUM_TYPE Size() { return last - first; }
	} Interval;
	typedef std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> coord;
	typedef std::pair<INMOST_DATA_REAL_TYPE, coord > wgt_coord;
	typedef std::vector< wgt_coord > wgt_coords;
	typedef struct row_col_t
	{
		Sparse::Row row, col;
		INMOST_DATA_REAL_TYPE diag;
	} row_col;
	typedef dynarray<INMOST_DATA_ENUM_TYPE,256> levels_t;
	//result of multilevel preconditioner
	//        |LDU  F |
	//  A  =  |       |
	//        | E   C |
	// LU decomposition is stored in skyline format with reversed direction
	//
	// For the next step of multilevel factorization we have to compute
	// S = C - (E U) D (L F)
	//
	// U is stored by rows
	// L is stored by columns
	// F is stored by rows
	// E is stored by rows
	// C is stored by rows
	//
	// During LF calculation F is traversed transposed (by columns) and
	// each column is solved with L. 
	// LF is obtained by columns.
	//
	// During EU calculating E is traversed by rows and
	// each row is solved with U. Then re
	//
	// For the faster multiplication of 
	std::vector<Sparse::Row::entry> LU_Entries, B_Entries;
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LU_Diag;
	interval<INMOST_DATA_ENUM_TYPE, Interval> U_Address, L_Address, B_Address;
	std::vector<Sparse::Row::entry> E_Entries, F_Entries;
	std::vector<interval<INMOST_DATA_ENUM_TYPE, Interval> *> E_Address;
	interval<INMOST_DATA_ENUM_TYPE, Interval> F_Address;
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> temp; // temporal place for solve phase
	levels_t level_size; //remember size of each level
	std::vector<Interval> level_interval;
	//reordering information
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > ddP,ddQ;
	INMOST_DATA_ENUM_TYPE reorder_nnz, ddpq_tau_adapt, estimator;
	INMOST_DATA_REAL_TYPE ddpq_tau, iluc2_tau;
	INMOST_DATA_REAL_TYPE tau, eps;
	INMOST_DATA_ENUM_TYPE sciters;
	Sparse::Matrix * Alink;
	Solver::OrderInfo * info;
	bool init;
	void DumpMatrix(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
									std::vector<Sparse::Row::entry> & Entries,
									INMOST_DATA_ENUM_TYPE wmbeg, INMOST_DATA_ENUM_TYPE wmend,
									std::string file_name);
	void SwapEntries(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
					 std::vector<Sparse::Row::entry> & Entries, 
					 INMOST_DATA_ENUM_TYPE rbeg, INMOST_DATA_ENUM_TYPE rend, 
					 INMOST_DATA_ENUM_TYPE k, INMOST_DATA_ENUM_TYPE j);
	void SwapLine(interval<INMOST_DATA_ENUM_TYPE, Interval> & Line, INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j);
	void SwapE(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j);
	
	void ReorderEF(INMOST_DATA_ENUM_TYPE mobeg, 
					INMOST_DATA_ENUM_TYPE cbeg,
					INMOST_DATA_ENUM_TYPE cend,
					INMOST_DATA_ENUM_TYPE wend, 
					interval<INMOST_DATA_ENUM_TYPE, bool> & donePQ, 
					interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
					interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ);
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
	ILUC_preconditioner(const ILUC_preconditioner & other);
	ILUC_preconditioner & operator =(ILUC_preconditioner const & other);
	ILUC_preconditioner(Solver::OrderInfo & info);
	bool isInitialized();
	bool isFinalized();
	bool Initialize();
	bool Finalize();
	void Multiply(int level, Sparse::Vector & input, Sparse::Vector & output);
	void ApplyLU(int level, Sparse::Vector & input,Sparse::Vector & output);
	void ApplyB(int level, double alpha, Sparse::Vector & x, double beta, Sparse::Vector & y);
	int Descend(int level, Sparse::Vector & input, Sparse::Vector & output);
	int Ascend(int level, Sparse::Vector & input, Sparse::Vector & output);
	bool Solve(Sparse::Vector & input, Sparse::Vector & output);
	bool ReplaceMAT(Sparse::Matrix & A);
	bool ReplaceSOL(Sparse::Vector & x);
	bool ReplaceRHS(Sparse::Vector & b);
	Method * Duplicate();
	~ILUC_preconditioner();
};

#endif //__SOLVER_DDPQILUC2__
