
#ifndef __SOLVER_MLMTILUC2__
#define __SOLVER_MLMTILUC2__
#include <iomanip>

#include "inmost_solver.h"
#include "../solver_prototypes.hpp"


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
		INMOST_DATA_ENUM_TYPE first, last, thr;
		INMOST_DATA_ENUM_TYPE Size() const { return last - first; }
	} Interval;
	typedef std::vector<INMOST_DATA_ENUM_TYPE> levels_t;
	std::vector< std::vector<Sparse::Row::entry> > L_Entries, U_Entries;
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LU_Diag;
	interval<INMOST_DATA_ENUM_TYPE, Interval> U_Address, L_Address;
	std::vector<interval<INMOST_DATA_ENUM_TYPE, Interval> *> F_Address, E_Address;
	std::vector< std::vector<Sparse::Row::entry> > E_Entries;
	std::vector< std::vector<Sparse::Row::entry> > F_Entries;
	levels_t level_size; //remember size of each level
	std::vector< std::vector<Block> > level_blocks;
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
	__INLINE bool check_zero(INMOST_DATA_REAL_TYPE u) { return tau + u == tau; }
	//__INLINE static bool check_zero(INMOST_DATA_REAL_TYPE u) { return 1 + u == 1; }
	//__INLINE static bool check_zero(INMOST_DATA_REAL_TYPE u) { return u == 0; }
	//__INLINE static bool check_zero(INMOST_DATA_REAL_TYPE u) { return fabs(u) < 1.0e-13; }
	void CheckBlock(const Block& b,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
		const std::vector< std::vector<Sparse::Row::entry> >& Entries,
		const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ,
		INMOST_DATA_ENUM_TYPE sepbeg,
		INMOST_DATA_ENUM_TYPE sepend,
		std::string file, int line);
	void CheckBlock(const Block& b,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
		const std::vector< std::vector<Sparse::Row::entry> >& Entries,
		INMOST_DATA_ENUM_TYPE sepbeg,
		INMOST_DATA_ENUM_TYPE sepend,
		std::string file, int line);
	void DumpMatrix(const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
					const std::vector< std::vector<Sparse::Row::entry> > & Entries,
					INMOST_DATA_ENUM_TYPE wmbeg, INMOST_DATA_ENUM_TYPE wmend,
					std::string file_name);
	INMOST_DATA_ENUM_TYPE ComputeNonzeroes(const Block& b, 
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
		const std::vector< std::vector<Sparse::Row::entry> >& Entries);
	INMOST_DATA_ENUM_TYPE ComputeNonzeroes(const Block& b,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
		const std::vector< std::vector<Sparse::Row::entry> >& Entries,
		const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ);
	void DumpMatrixBlock( const Block & b, 
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
		const std::vector< std::vector<Sparse::Row::entry> >& Entries,
		const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ,
		std::string file_name);
	void DumpMatrixBlock(const Block& b, 
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
		const std::vector< std::vector<Sparse::Row::entry> >& Entries,
		std::string file_name);
	void CheckOrder(const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
					const std::vector< std::vector<Sparse::Row::entry> >& Entries,
					INMOST_DATA_ENUM_TYPE rbeg, INMOST_DATA_ENUM_TYPE rend);
	void inversePQ(INMOST_DATA_ENUM_TYPE wbeg,
				   INMOST_DATA_ENUM_TYPE wend, 
				   const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				   const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ);

	void applyPQ(INMOST_DATA_ENUM_TYPE wbeg,
				 INMOST_DATA_ENUM_TYPE wend,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
				 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ);
	void ReorderEF(INMOST_DATA_ENUM_TYPE wbeg,
				   INMOST_DATA_ENUM_TYPE wend,
				   interval<INMOST_DATA_ENUM_TYPE, bool> & donePQ,
				   const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				   const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ);
	INMOST_DATA_REAL_TYPE AddListOrdered(INMOST_DATA_ENUM_TYPE cbeg,
										 const Interval & Address,
										 const std::vector<Sparse::Row::entry> & Entries,
										 INMOST_DATA_REAL_TYPE coef,
										 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
										 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
										 INMOST_DATA_REAL_TYPE droptol);
	INMOST_DATA_REAL_TYPE AddListUnordered(INMOST_DATA_ENUM_TYPE & Sbeg,
										   const Interval & Address,
										   const std::vector<Sparse::Row::entry> & Entries,
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
									 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
									 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
									 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Est,
									 INMOST_DATA_REAL_TYPE & mu_update);
	INMOST_DATA_REAL_TYPE Estimator2(INMOST_DATA_ENUM_TYPE k,
									 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
									 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
									 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Est,
									 INMOST_DATA_REAL_TYPE & mu_update);
	void EstimatorUpdate(INMOST_DATA_ENUM_TYPE k,
						 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
						 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
						 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Est,
						 INMOST_DATA_REAL_TYPE & mu_update);
	void DiagonalUpdate(INMOST_DATA_ENUM_TYPE k,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Diag,
						const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndecesL,
						const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValuesL,
						const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeceU,
						const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValuesU);
	void ClearList(INMOST_DATA_ENUM_TYPE beg,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues);
	void PrepareGraph(INMOST_DATA_ENUM_TYPE wbeg,
					  INMOST_DATA_ENUM_TYPE wend,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
		const std::vector< std::vector<Sparse::Row::entry> >& Entries,
		interval<INMOST_DATA_ENUM_TYPE, Interval>& G_Address,
		std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& G_Entries);
	// tG should be preallocated by number of columns, that may be wider then wbeg:wend
	void PrepareGraphTranspose(INMOST_DATA_ENUM_TYPE wbeg,
							   INMOST_DATA_ENUM_TYPE wend,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& G_Address,
		const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& G_Entries,
		interval<INMOST_DATA_ENUM_TYPE, Interval>& tG_Address,
		std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& tG_Entries);
	// pG should be preallocated by wbeg:wend, tind computed with PrepareGraphTranspose
	void PrepareGraphProduct(INMOST_DATA_ENUM_TYPE wbeg,
							 INMOST_DATA_ENUM_TYPE wend,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& G_Address,
		const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& G_Entries,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& tG_Address,
		const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& tG_Entries,
		interval<INMOST_DATA_ENUM_TYPE, Interval>& pG_Address,
		std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& pG_Entries);
	// moves wbeg:wend subset of rows from tG_in to tG_out and filters entries outside wbeg:wend
	// tG_in is not preserved to save memory
	void FilterGraphProduct(const Block & b,
							const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& invP,
							const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localP,
							const interval<INMOST_DATA_ENUM_TYPE, Interval>& in_pG_Address,
							const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > in_pG_Entries,
							interval<INMOST_DATA_ENUM_TYPE, Interval>& out_pG_Address,
							std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > out_pG_Entries);
	void FilterGraphTranspose(const Block & b,
							const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localP,
							const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& invQ,
							const interval<INMOST_DATA_ENUM_TYPE, Interval>& in_tG_Address,
							const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > in_tG_Entries,
							interval<INMOST_DATA_ENUM_TYPE, Interval>& out_tG_Address,
							std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > out_tG_Entries);
	void FilterGraph(const Block & b,
		const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& invP,
		const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& in_G_Address,
		const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > in_G_Entries,
		interval<INMOST_DATA_ENUM_TYPE, Interval>& out_G_Address,
		std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > out_G_Entries);
	// finds permutation that separates matrix into blocks by running GreedyDissection recursively
	void NestedDissection(INMOST_DATA_ENUM_TYPE wbeg, 
					 	  INMOST_DATA_ENUM_TYPE wend,
						  const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						  const std::vector< std::vector<Sparse::Row::entry> > & Entries,
						  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
						  std::vector<Block> & blocks,
						  INMOST_DATA_ENUM_TYPE max_size);
	void KwayDissection(INMOST_DATA_ENUM_TYPE wbeg, 
					 	INMOST_DATA_ENUM_TYPE wend,
						const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						const std::vector< std::vector<Sparse::Row::entry> > & Entries,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
						std::vector<Block> & blocks,
						int parts);
	void KwaySymmetricDissection(INMOST_DATA_ENUM_TYPE wbeg,
		INMOST_DATA_ENUM_TYPE wend,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
		const std::vector< std::vector<Sparse::Row::entry> >& Entries,
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localP,
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ,
		std::vector<Block>& blocks,
		int parts);
	// finds permutation that separates matrix into blocks
	void GreedyDissection(const Block & b,
						  const interval< INMOST_DATA_ENUM_TYPE, Interval >& G_Address,
						  const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& G_Entries,
						  const interval< INMOST_DATA_ENUM_TYPE, Interval >& tG_Address,
						  const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& tG_Entries,
						  const interval< INMOST_DATA_ENUM_TYPE, Interval >& pG_Address,
						  const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& pG_Entries,
						  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
						  std::vector<Block> & blocks, int kway_parts);
						//   int wgt_sep, int wgt_blk, int kway_parts);
	//compute interval of column indices
	void ColumnInterval(INMOST_DATA_ENUM_TYPE wbeg, 
						 INMOST_DATA_ENUM_TYPE wend,
						 const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						 const std::vector< std::vector<Sparse::Row::entry> > & Entries,
						 INMOST_DATA_ENUM_TYPE & cbeg,
						 INMOST_DATA_ENUM_TYPE & cend);
	// column permutations
	// has to call functions to permute scaling and E,F blocks
	void ReorderMatrixQ(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						std::vector< std::vector<Sparse::Row::entry> > & Entries,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ);

	// rows and columns permutations
	// has to call functions to permute scaling and E,F blocks
	void ReorderMatrixPQ(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						 std::vector< std::vector<Sparse::Row::entry> > & Entries,
						 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
						 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP);
	
	void PrepareTranspose(INMOST_DATA_ENUM_TYPE cbeg,
						  INMOST_DATA_ENUM_TYPE cend,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
		const std::vector< std::vector<Sparse::Row::entry> >& Entries,
		interval<INMOST_DATA_ENUM_TYPE, Interval> G_Address,
		std::vector< std::vector< std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> > > G_Entries);
						  
						  
	void MaximalTransversal(INMOST_DATA_ENUM_TYPE rbeg,
							INMOST_DATA_ENUM_TYPE rend,
							INMOST_DATA_ENUM_TYPE cbeg,
							INMOST_DATA_ENUM_TYPE cend,
							const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
							const std::vector< std::vector<Sparse::Row::entry> > & Entries,
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & DL,
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & DR);
	void FillColumnGaps(INMOST_DATA_ENUM_TYPE cbeg,
						INMOST_DATA_ENUM_TYPE cend,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ);
	
	void SymmetricGraph(INMOST_DATA_ENUM_TYPE wbeg,
						INMOST_DATA_ENUM_TYPE wend,
						const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						const std::vector< std::vector<Sparse::Row::entry> > & Entries,
						const interval<INMOST_DATA_ENUM_TYPE, bool>& Pivot,
						std::vector<INMOST_DATA_ENUM_TYPE> & xadj,
						std::vector< INMOST_DATA_ENUM_TYPE > & adjncy);
	void GraphWeights(	INMOST_DATA_ENUM_TYPE wbeg,
						INMOST_DATA_ENUM_TYPE wend,
						const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
						const std::vector< std::vector<Sparse::Row::entry> >& Entries,
						const interval<INMOST_DATA_ENUM_TYPE, bool>& Pivot,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DL,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DR,
						std::vector<INMOST_DATA_REAL_TYPE>& wadj);
							
	void ReorderColumns(INMOST_DATA_ENUM_TYPE wbeg,
						INMOST_DATA_ENUM_TYPE wend,
						interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						std::vector< std::vector<Sparse::Row::entry> > & Entries,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & DR,
						Sparse::Vector & DR0);

	void ReorderSystem(	INMOST_DATA_ENUM_TYPE wbeg,
						INMOST_DATA_ENUM_TYPE wend,
						interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
						std::vector< std::vector<Sparse::Row::entry> >& Entries,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localP,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DL,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DR,
						Sparse::Vector& DL0,
						Sparse::Vector& DR0);

	void SinkhornScaling(	INMOST_DATA_ENUM_TYPE wbeg,
							INMOST_DATA_ENUM_TYPE wend,
							interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
							std::vector<Sparse::Row::entry>& Entries,
							interval<INMOST_DATA_ENUM_TYPE, bool>& Pivot,
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DL,
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DR,
							INMOST_DATA_REAL_TYPE p);
	void IdominanceScaling(	INMOST_DATA_ENUM_TYPE wbeg,
							INMOST_DATA_ENUM_TYPE wend,
							const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
							std::vector< std::vector<Sparse::Row::entry> >& Entries,
							const interval<INMOST_DATA_ENUM_TYPE, bool>& Pivot,
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DL,
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DR);
	
	void PivotScaling(	INMOST_DATA_ENUM_TYPE wbeg,
						INMOST_DATA_ENUM_TYPE wend,
						const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
						std::vector< std::vector<Sparse::Row::entry> >& Entries,
						const interval<INMOST_DATA_ENUM_TYPE, bool>& Pivot,
						const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DL,
						const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DR);

	void EFScaling(	INMOST_DATA_ENUM_TYPE wbeg,
					INMOST_DATA_ENUM_TYPE wend,
					INMOST_DATA_ENUM_TYPE mobeg,
					INMOST_DATA_ENUM_TYPE moend,
					const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DL,
					const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DR);
						
						
	void WRCMOrdering(	INMOST_DATA_ENUM_TYPE wbeg,
						INMOST_DATA_ENUM_TYPE wend,
						std::vector<INMOST_DATA_ENUM_TYPE> & xadj, 
						std::vector<INMOST_DATA_ENUM_TYPE> & adjncy,
						std::vector<INMOST_DATA_REAL_TYPE> & wadj,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ);

	void RCMOrdering(	INMOST_DATA_ENUM_TYPE wbeg,
						INMOST_DATA_ENUM_TYPE wend,
						std::vector<INMOST_DATA_ENUM_TYPE> & xadj, 
						std::vector<INMOST_DATA_ENUM_TYPE> & adjncy,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ);				
	void MetisOrdering(	INMOST_DATA_ENUM_TYPE wbeg,
						INMOST_DATA_ENUM_TYPE wend,
						std::vector<INMOST_DATA_ENUM_TYPE> & xadj, 
						std::vector<INMOST_DATA_ENUM_TYPE> & adjncy,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ);

	void Factorize(INMOST_DATA_ENUM_TYPE cbeg,
					INMOST_DATA_ENUM_TYPE cend,
					interval<INMOST_DATA_ENUM_TYPE, Interval>& A_Address,
					std::vector< std::vector<Sparse::Row::entry> > & A_Entries,
					interval<INMOST_DATA_ENUM_TYPE, Interval>& L2_Address,
					std::vector< std::vector<Sparse::Row::entry> >& L2_Entries,
					interval<INMOST_DATA_ENUM_TYPE, Interval>& U2_Address,
					std::vector< std::vector<Sparse::Row::entry> >& U2_Entries,
					interval<INMOST_DATA_ENUM_TYPE, bool>& Pivot,
					bool block_pivot,
					INMOST_DATA_REAL_TYPE& NuLout,
					INMOST_DATA_REAL_TYPE& NuUout,
					double& testimator);

	void CheckColumnGaps(const Block& b, const interval<INMOST_DATA_ENUM_TYPE, Interval>& A_Address, const std::vector< std::vector<Sparse::Row::entry> >& A_Entries);
	
	void DumpGraph(std::string name, const interval<INMOST_DATA_ENUM_TYPE, Interval>& G_Address,
		const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > G_Entries);

	int Thread();
	int Threads();

	// Sparse::Vector div;
	INMOST_DATA_ENUM_TYPE progress_all, progress_cur;
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
	int Descend(INMOST_DATA_ENUM_TYPE level, Sparse::Vector & inout);
	int Ascend(INMOST_DATA_ENUM_TYPE level, Sparse::Vector & inout);
	bool Solve(Sparse::Vector & input, Sparse::Vector & output);
	bool ReplaceMAT(Sparse::Matrix & A);
	bool ReplaceSOL(Sparse::Vector & x);
	bool ReplaceRHS(Sparse::Vector & b);
	Method * Duplicate();
	~MLMTILUC_preconditioner();
};

#endif //__SOLVER_MLMTILUC2__
