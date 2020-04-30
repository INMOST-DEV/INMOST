
#ifndef __SOLVER_MLMTILUC2__
#define __SOLVER_MLMTILUC2__
#include <iomanip>

#include "inmost_solver.h"
#include "../solver_prototypes.hpp"

/*
class Range
{
	INMOST_DATA_ENUM_TYPE first, last;
public:
	Range(const Range & b):first(b.first), last(b.last) {}
	Range(INMOST_DATA_ENUM_TYPE first, INMOST_DATA_ENUM_TYPE last) : first(first), last(last) {}
	Range & operator =(Range const & b) {first = b.first; last = b.last; return *this;}
	~Range() {}
	INMOST_DATA_ENUM_TYPE Size() { return last - first; }
	INMOST_DATA_ENUM_TYPE First() {return first;}
	INMOST_DATA_ENUM_TYPE Last() {return last;}
};

class CSR_Row
{
	Range & range;
	std::vector<Sparse::Row::entry> & entries;
public:
}

class CSR_Matrix
{
	interval<INMOST_DATA_ENUM_TYPE, Range> Address;
	std::vector<Sparse::Row::entry> Entries;
public:
};
*/

class MLMTILUC_preconditioner : public Method
{
	typedef struct Block_t
	{
		INMOST_DATA_ENUM_TYPE row_start, row_end;
		INMOST_DATA_ENUM_TYPE col_start, col_end;
		bool separator;
		Block_t(INMOST_DATA_ENUM_TYPE rows, 
				INMOST_DATA_ENUM_TYPE rowe,
				INMOST_DATA_ENUM_TYPE cols,
				INMOST_DATA_ENUM_TYPE cole,
				bool separator)
		: row_start(rows), row_end(rowe),
		  col_start(cols), col_end(cole),
		  separator(separator) {}
		Block_t(INMOST_DATA_ENUM_TYPE rows, 
				INMOST_DATA_ENUM_TYPE rowe,
				INMOST_DATA_ENUM_TYPE cols,
				INMOST_DATA_ENUM_TYPE cole)
		: row_start(rows), row_end(rowe),
		  col_start(cols), col_end(cole),
		  separator(false) {}
		INMOST_DATA_ENUM_TYPE RowSize() const {return row_end - row_start;}
		INMOST_DATA_ENUM_TYPE ColSize() const {return col_end - col_start;}
		// Block_t(const Block_t & b)
		// : row_start(b.row_start), row_end(b.row_end),
		//   col_start(b.col_start), col_end(b.col_end) {}
		// Block_t & operator =(Block_t const & b)
		// {
		// 	row_start = b.row_start;
		// 	row_end = b.row_end;
		// 	col_start = b.col_start;
		// 	col_end = b.col_end;
		// 	return *this;
		// }
	} Block;
	typedef struct Interval_t
	{
		INMOST_DATA_ENUM_TYPE first, last;
		INMOST_DATA_ENUM_TYPE Size() const { return last - first; }
	} Interval;
	typedef dynarray<INMOST_DATA_ENUM_TYPE,256> levels_t;
	std::vector<Sparse::Row::entry> L_Entries, U_Entries, B_Entries;
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
	INMOST_DATA_ENUM_TYPE sciters, verbosity;
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
	INMOST_DATA_REAL_TYPE AddListOrdered(INMOST_DATA_ENUM_TYPE cbeg,
										 Interval & Address,
										 std::vector<Sparse::Row::entry> & Entries,
										 INMOST_DATA_REAL_TYPE coef,
										 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
										 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
										 INMOST_DATA_REAL_TYPE droptol);
	INMOST_DATA_REAL_TYPE AddListUnordered(INMOST_DATA_ENUM_TYPE & Sbeg,
										   Interval & Address,
										   std::vector<Sparse::Row::entry> & Entries,
										   INMOST_DATA_REAL_TYPE coef,
										   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
										   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
										   INMOST_DATA_REAL_TYPE droptol);
	void OrderList(INMOST_DATA_ENUM_TYPE & Sbeg, 
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
				   std::vector<INMOST_DATA_ENUM_TYPE> & indices);
	void ScaleList(INMOST_DATA_REAL_TYPE coef,
				   INMOST_DATA_ENUM_TYPE Sbeg, 
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues);
	INMOST_DATA_REAL_TYPE Estimator1(INMOST_DATA_ENUM_TYPE k,
									 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
									 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
									 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Est,
									 INMOST_DATA_REAL_TYPE & mu_update);
	INMOST_DATA_REAL_TYPE Estimator2(INMOST_DATA_ENUM_TYPE k,
									 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
									 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
									 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Est,
									 INMOST_DATA_REAL_TYPE & mu_update);
	void EstimatorUpdate(INMOST_DATA_ENUM_TYPE k,
						 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
						 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
						 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Est,
						 INMOST_DATA_REAL_TYPE & mu_update);
	void DiagonalUpdate(INMOST_DATA_ENUM_TYPE k,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Diag,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndecesL,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValuesL,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeceU,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValuesU);
	void ClearList(INMOST_DATA_ENUM_TYPE beg,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues);
	void PrepareGraph(INMOST_DATA_ENUM_TYPE wbeg,
					  INMOST_DATA_ENUM_TYPE wend,
					  const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
					  const std::vector<Sparse::Row::entry> & Entries,
					  interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G);
	// tG should be preallocated by number of columns, that may be wider then wbeg:wend
	void PrepareGraphTranspose(INMOST_DATA_ENUM_TYPE wbeg,
							   INMOST_DATA_ENUM_TYPE wend,
							   const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G,
							   interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG);
	// pG should be preallocated by wbeg:wend, tind computed with PrepareGraphTranspose
	void PrepareGraphProduct(INMOST_DATA_ENUM_TYPE wbeg,
							 INMOST_DATA_ENUM_TYPE wend,
							 const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G,
							 const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG,
							 interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & pG);
	// moves wbeg:wend subset of rows from tG_in to tG_out and filters entries outside wbeg:wend
	// tG_in is not preserved to save memory
	void FilterGraphProduct(const Block & b,
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				 		  	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
							interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & pG_in,
							interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & pG_out);
	void FilterGraphTranspose(const Block & b,
							  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		  	  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ,
							  interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG_in,
							  interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG_out);
	void FilterGraph(const Block & b,
					 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				 	 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
					 interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G_in,
					 interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G_out);
	// finds permutation that separates matrix into blocks by running GreedyDissection recursively
	void NestedDissection(INMOST_DATA_ENUM_TYPE wbeg, 
					 	  INMOST_DATA_ENUM_TYPE wend,
						  const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						  const std::vector<Sparse::Row::entry> & Entries,
						  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
						  std::vector<Block> & blocks,
						  INMOST_DATA_ENUM_TYPE max_size);
	void KwayDissection(INMOST_DATA_ENUM_TYPE wbeg, 
					 	INMOST_DATA_ENUM_TYPE wend,
						const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						const std::vector<Sparse::Row::entry> & Entries,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
						std::vector<Block> & blocks,
						int parts);
	// finds permutation that separates matrix into blocks
	void GreedyDissection(const Block & b,
						  const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G,
						  const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG,
						  const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & pG,
						  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
						  std::vector<Block> & blocks, int kway_parts);
						//   int wgt_sep, int wgt_blk, int kway_parts);
	//compute interval of column indices
	void ColumnInterval(INMOST_DATA_ENUM_TYPE wbeg, 
						 INMOST_DATA_ENUM_TYPE wend,
						 const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						 const std::vector<Sparse::Row::entry> & Entries,
						 INMOST_DATA_ENUM_TYPE & cbeg,
						 INMOST_DATA_ENUM_TYPE & cend);
	// column permutations
	// has to call functions to permute scaling and E,F blocks
	void ReorderMatrixQ(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						std::vector<Sparse::Row::entry> & Entries,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ);

	// rows and columns permutations
	// has to call functions to permute scaling and E,F blocks
	void ReorderMatrixPQ(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						 std::vector<Sparse::Row::entry> & Entries,
						 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
						 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP);
	
	void PrepareTranspose(INMOST_DATA_ENUM_TYPE cbeg,
						  INMOST_DATA_ENUM_TYPE cend,
						  interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						  std::vector<Sparse::Row::entry> & Entries,
						  interval<INMOST_DATA_ENUM_TYPE, std::vector< std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> > > & Indices );
	void DumpGraph(std::string name, interval<INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G);

	int Thread();
	int Threads();
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
