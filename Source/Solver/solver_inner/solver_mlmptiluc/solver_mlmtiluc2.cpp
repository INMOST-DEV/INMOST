#define _CRT_SECURE_NO_WARNINGS
#include "inmost_solver.h"
#if defined(USE_SOLVER)
#include "solver_mlmtiluc2.hpp"
#include <sstream>
#include <deque>
#include <iomanip>

#include "../../../Misc/utils.h"
//#define USE_OMP

using namespace INMOST;

#ifndef DEFAULT_TAU
#define DEFAULT_TAU 0.01
#endif

#if defined(USE_OMP)
#define USE_OMP_FACT
#endif //USE_OMP

//~ #undef USE_OMP


//~ #define REORDER_RCM
#define REORDER_WRCM
//~ #define REORDER_BRCM
#if defined(USE_SOLVER_METIS)
//~ #define REORDER_METIS_ND
#endif
#if defined(USE_SOLVER_MONDRIAAN)
//#define REORDER_MONDRIAAN
#endif
//#define REORDER_ZOLTAN_HUND


static bool run_mpt = true;
static bool rescale_b = true;
static bool allow_pivot = true;
const INMOST_DATA_ENUM_TYPE UNDEF = ENUMUNDEF, EOL = ENUMUNDEF - 1;

//~ #define PREMATURE_DROPPING

//~ #define EQUALIZE_1NORM
//~ #define EQUALIZE_2NORM
#define EQUALIZE_IDOMINANCE

#define PIVOT_THRESHOLD
#define PIVOT_THRESHOLD_VALUE 1.0e-9
//#define DIAGONAL_PERTURBATION
#define DIAGONAL_PERTURBATION_REL 1.0e-7
#define DIAGONAL_PERTURBATION_ABS 1.0e-10
//#define ILUC2
//~ #define ILUC2_SCHUR

//#define PIVOT_COND_DEFAULT 0.1/tau
#define PIVOT_COND_DEFAULT 1.0e+2
#define PIVOT_DIAG_DEFAULT 1.0e+5
//~ #define SCHUR_DROPPING_LF
//~ #define SCHUR_DROPPING_EU
#define SCHUR_DROPPING_S
#define SCHUR_DROPPING_LF_PREMATURE
#define SCHUR_DROPPING_EU_PREMATURE
// #define SCHUR_DROPPING_S_PREMATURE
#define CONDITION_PIVOT
// #define NNZ_PIVOT

#define NNZ_GROWTH_PARAM 1.05

#if defined(REORDER_METIS_ND)
#define METIS_EXPORT
#include "metis.h"
#endif
#if defined(REORDER_ZOLTAN_HUND)
#include <zoltan.h>
#endif
#if defined(REORDER_MONDRIAAN)
#include <Mondriaan.h>
#endif

	template<typename T>
	std::ostream & fmt(std::ostream & in, const T & value, int width, int precision)
	{
		std::ios save(NULL);
		save.copyfmt(in);
		in << std::fixed << std::setw(width) << std::setprecision(precision) << value;
		std::cout.copyfmt(save);
		return in;
	}


	void MLMTILUC_preconditioner::ReorderEF(INMOST_DATA_ENUM_TYPE wbeg,
											INMOST_DATA_ENUM_TYPE wend,
											interval<INMOST_DATA_ENUM_TYPE, bool> & donePQ,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ)
	{
		INMOST_DATA_ENUM_TYPE i, k, l;
		if( !E_Address.empty() )
		{
			for(k = wbeg; k < wend; ++k) donePQ[k] = false;
			i = wbeg;
			while (i < wend)
			{
				if (donePQ[i]) i++;
				else
				{
					if (localP[i] != i)
					{
						k = i;
						do
						{
							for(l = 0; l < E_Address.size(); ++l)
							{
								Interval t = E_Address[l]->at(i);
								E_Address[l]->at(i) = E_Address[l]->at(localP[k]);
								E_Address[l]->at(localP[k]) = t;
							}
							donePQ[localP[k]] = true;
							k = localP[k];
						} while (k != i);
					}
					donePQ[i] = true;
				}
			}
		}
		if( !F_Address.empty() )
		{
			for(k = wbeg; k < wend; ++k) donePQ[k] = false;
			i = wbeg;
			while (i < wend)
			{
				if (donePQ[i]) i++;
				else
				{
					if (localQ[i] != i)
					{
						k = i;
						do
						{
							for(l = 0; l < F_Address.size(); ++l)
							{
								Interval t = F_Address[l]->at(i);
								F_Address[l]->at(i) = F_Address[l]->at(localQ[k]);
								F_Address[l]->at(localQ[k]) = t;
							}
							donePQ[localQ[k]] = true;
							k = localQ[k];
						} while (k != i);
					}
					donePQ[i] = true;
				}
			}
		}
	}


	void MLMTILUC_preconditioner::DumpMatrix(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
											 std::vector<Sparse::Row::entry> & Entries,
											 INMOST_DATA_ENUM_TYPE wmbeg, INMOST_DATA_ENUM_TYPE wmend,
											 std::string file_name)
	{
		INMOST_DATA_REAL_TYPE norm1 = 0, norm2 = 0, max = -1.0e54, min = 1.0e54, minabs = 1.0e54;
		INMOST_DATA_REAL_TYPE vrowmax, diag, mindiag = 1.0e54, maxdiag = -1.0e54, maxabsdiag = -1.0e54, minabsdiag = 1.0e54;
		INMOST_DATA_ENUM_TYPE nnz = 0, dominant_rows = 0, dominant_cols = 0, irowmax = 0, nodiag = 0;
		INMOST_DATA_ENUM_TYPE addrbeg = Address.get_interval_beg();
		INMOST_DATA_ENUM_TYPE addrend = Address.get_interval_end();
		interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> vcolmax(wmbeg,wmend,0);
		interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> icolmax(wmbeg,wmend,ENUMUNDEF);
		for (INMOST_DATA_ENUM_TYPE k = wmbeg; k < wmend; ++k) 
		{
			if( k < addrbeg || k >= addrend) continue;
			nnz += Address[k].Size();
			vrowmax = 0;

			bool diag_found = false;
			diag = 0;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it)
			{
				norm1 += fabs(Entries[it].second);
				norm2 += Entries[it].second*Entries[it].second;
				if( fabs(Entries[it].second) > vrowmax ) 
				{
					vrowmax = fabs(Entries[it].second);
					irowmax = Entries[it].first;
				}

				if( fabs(Entries[it].second) > vcolmax[Entries[it].first] ) 
				{
					vcolmax[Entries[it].first] = fabs(Entries[it].second);
					icolmax[Entries[it].first] = k;
				}
				if( Entries[it].second > max ) max = Entries[it].second;
				if( Entries[it].second < min ) min = Entries[it].second;
				if( fabs(Entries[it].second) < fabs(minabs) ) minabs = Entries[it].second;

				if( Entries[it].first == k )
				{
					diag_found = true;
					diag = Entries[it].second;
				}
			}

			if( diag_found )
			{
				if( mindiag > diag ) mindiag = diag;
				if( maxdiag < diag ) maxdiag = diag;
				if( minabsdiag > fabs(diag) ) minabsdiag = fabs(diag);
				if( maxabsdiag < fabs(diag) ) maxabsdiag = fabs(diag);
			}
			else nodiag++;

			if( irowmax == k ) ++dominant_rows;
		}

		for (INMOST_DATA_ENUM_TYPE k = wmbeg; k < wmend; ++k) if( icolmax[k] == k ) ++dominant_cols;

		std::cout << "Writing matrix to " << file_name.c_str() << std::endl;
		std::fstream fout(file_name.c_str(),std::ios::out);
		fout << "%%MatrixMarket matrix coordinate real general" << std::endl;
		fout << "% maximum " << max << std::endl;
		fout << "% minimum " << min << std::endl;
		fout << "% absolute minimum " << minabs << std::endl;
		fout << "% A 1-norm  " << norm1 << std::endl;
		fout << "% A 2-norm  " << sqrt(norm2) << std::endl;
		fout << "% mean 1-norm  " << norm1/(wmend-wmbeg) << std::endl;
		fout << "% mean 2-norm  " << sqrt(norm2/(wmend-wmbeg)) << std::endl;
		fout << "% dominant rows  " << dominant_rows << std::endl;
		fout << "% dominant cols  " << dominant_cols << std::endl;
		fout << "% maximal diagonal value " << maxdiag << std::endl;
		fout << "% minimal diagonal value " << mindiag << std::endl;
		fout << "% absolute maximal diagonal value " << maxabsdiag << std::endl;
		fout << "% absolute minimal diagonal value " << minabsdiag << std::endl;
		fout << "% no diagonal value  " << nodiag << std::endl;
		fout << "% true matrix indices interval " << wmbeg << ":" << wmend << std::endl;
		fout << std::scientific;
		
		//fout.close(); return;
		
		fout << wmend-wmbeg << " " << wmend-wmbeg << " " << nnz << std::endl;;
		for (INMOST_DATA_ENUM_TYPE k = wmbeg; k < wmend; ++k)
		{
			if( k < addrbeg || k >= addrend) continue;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it)
				fout << (k-wmbeg+1) << " " << (Entries[it].first-wmbeg+1) << " " << Entries[it].second << std::endl;
		}
		fout.close();
	}
	void MLMTILUC_preconditioner::CheckOrder(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
											 std::vector<Sparse::Row::entry> & Entries,
											 INMOST_DATA_ENUM_TYPE rbeg, INMOST_DATA_ENUM_TYPE rend)
	{
		INMOST_DATA_ENUM_TYPE i,r;
		for (i = rbeg; i < rend; i++) if( Address[i].first < Address[i].last )
		{
			//check ordered on entry
			bool ordered = true;
			for (r = Address[i].first; r < Address[i].last-1; r++)
			{
				if( !(Entries[r].first < Entries[r+1].first) )
				{
					ordered = false;
					break;
				}
			}
			if( !ordered )
			{
				std::cout << "Row " << i << " not ordered: " << std::endl;
				std::cout << "Interval: " << Address[i].first << ":" << Address[i].last << std::endl;
				for (r = Address[i].first; r < Address[i].last; r++)
				{
					std::cout << "(" << Entries[r].first << "," << Entries[r].second << ") ";
				}
				std::cout << std::endl;
				throw -1;
			}
			bool nan = false;
			for (r = Address[i].first; r < Address[i].last; r++)
			{
				if( Entries[r].second != Entries[r].second )
				{
					nan = true;
					break;
				}
			}
			if( nan )
			{
				std::cout << "Row " << i << " contains nan: " << std::endl;
				std::cout << "Interval: " << Address[i].first << ":" << Address[i].last << std::endl;
				for (r = Address[i].first; r < Address[i].last; r++)
				{
					std::cout << "(" << Entries[r].first << "," << Entries[r].second << ") ";
				}
				std::cout << std::endl;
				throw -1;
			}
		}
	}
	
	void MLMTILUC_preconditioner::inversePQ(INMOST_DATA_ENUM_TYPE wbeg,
											INMOST_DATA_ENUM_TYPE wend,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ)
	{
		//inverse reordering
		// in invPQ numbers indicate where to get current column
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) invP[localP[k]] = k;
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) invQ[localQ[k]] = k;
	}
	void MLMTILUC_preconditioner::applyPQ(INMOST_DATA_ENUM_TYPE wbeg,
										  INMOST_DATA_ENUM_TYPE wend,
										  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
										  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
										  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
										  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ)
	{
		INMOST_DATA_ENUM_TYPE k;
		// compute reordering in global P,Q, we need it to compute reordering in vector during solve phase
		for (k = wbeg; k < wend; ++k)
		{
			localP[k] = ddP[invP[k]];
			localQ[k] = ddQ[invQ[k]];
		}
		// update reordering in global P,Q
		for (k = wbeg; k < wend; ++k)
		{
			ddP[k] = localP[k];
			ddQ[k] = localQ[k];
		}
	}
	INMOST_DATA_REAL_TYPE MLMTILUC_preconditioner::AddListOrdered(INMOST_DATA_ENUM_TYPE cbeg, 
																  Interval & Address,
																  std::vector<Sparse::Row::entry> & Entries,
																  INMOST_DATA_REAL_TYPE coef,
																  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
																  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
																  INMOST_DATA_REAL_TYPE droptol)
	{
		INMOST_DATA_ENUM_TYPE curr = cbeg, next;
		INMOST_DATA_REAL_TYPE drop_sum = 0;
		for (INMOST_DATA_ENUM_TYPE it = Address.first; it < Address.last; ++it)
		{
			//~ assert(!Pivot[i]); // pivoted rows should be empty
			//~ assert(j >= k); // all indices are supposed to be to the right
			INMOST_DATA_ENUM_TYPE j = Entries[it].first;
			INMOST_DATA_REAL_TYPE u = coef*Entries[it].second;
			//eliminate values
			if (LineIndeces[j] != UNDEF) //there is an entry in the list
				LineValues[j] += u;
			else if( fabs(u) > droptol )
			{
				next = curr;
				while (next < j)
				{
					curr = next;
					next = LineIndeces[curr];
				}
				assert(curr < j);
				assert(j < next);
				LineIndeces[curr] = j;
				LineIndeces[j] = next;
				LineValues[j] = u; 
			}
			else drop_sum += u;
		}
		return drop_sum;
	}
	INMOST_DATA_REAL_TYPE MLMTILUC_preconditioner::AddListUnordered(INMOST_DATA_ENUM_TYPE & Sbeg,
																	Interval & Address,
																	std::vector<Sparse::Row::entry> & Entries,
																	INMOST_DATA_REAL_TYPE coef,
																	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
																	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
																	INMOST_DATA_REAL_TYPE droptol)
	{
		INMOST_DATA_REAL_TYPE drop_sum = 0;
		for (INMOST_DATA_ENUM_TYPE it = Address.first; it < Address.last; ++it)
		{
			//~ assert(!Pivot[i]); // pivoted rows should be empty
			//~ assert(j >= k); // all indices are supposed to be to the right
			INMOST_DATA_ENUM_TYPE j = Entries[it].first;
			INMOST_DATA_REAL_TYPE u = coef*Entries[it].second;
			//eliminate values
			if (LineIndeces[j] != UNDEF) //there is an entry in the list
				LineValues[j] += u;
			else if( fabs(u) > droptol )
			{
				LineValues[j] = u;
				LineIndeces[j] = Sbeg;
				Sbeg = j;
			}
			else drop_sum += u;
		}
		return drop_sum;
	}
	void MLMTILUC_preconditioner::OrderList(INMOST_DATA_ENUM_TYPE & Sbeg,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
											std::vector<INMOST_DATA_ENUM_TYPE> & indices)
	{
		indices.clear();
		INMOST_DATA_ENUM_TYPE i = Sbeg;
		while(i != EOL)
		{
			indices.push_back(i);
			i = LineIndeces[i];
		}
		if( !indices.empty() )
		{
			std::sort(indices.begin(),indices.end());
			Sbeg = indices[0];
			for(size_t qt = 1; qt < indices.size(); ++qt)
				LineIndeces[indices[qt-1]] = indices[qt];
			LineIndeces[indices.back()] = EOL;
		}
	}
	
	void MLMTILUC_preconditioner::ScaleList(INMOST_DATA_REAL_TYPE coef,
											INMOST_DATA_ENUM_TYPE beg,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues)
	{										
		INMOST_DATA_ENUM_TYPE i = beg;
		while (i != EOL)
		{
			LineValues[i] *= coef;
			i = LineIndeces[i];
		}
	}
	INMOST_DATA_REAL_TYPE MLMTILUC_preconditioner::Estimator1(INMOST_DATA_ENUM_TYPE k,
															  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
															  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
															  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Est,
															  INMOST_DATA_REAL_TYPE & mu_update)
	{
		INMOST_DATA_REAL_TYPE mup = 1.0 - Est[k];
		INMOST_DATA_REAL_TYPE mum = -1.0 - Est[k];
		INMOST_DATA_REAL_TYPE smup = 0, smum = 0;
		INMOST_DATA_ENUM_TYPE i = LineIndeces[k];
		while (i != EOL)
		{
			smup += fabs(Est[i] + LineValues[i] * mup);
			smum += fabs(Est[i] + LineValues[i] * mum);
			i = LineIndeces[i];
		}
		if (smup > smum) mu_update = mup; else mu_update = mum;
		return std::max(fabs(mup),fabs(mum));
	}
	
	
	INMOST_DATA_REAL_TYPE MLMTILUC_preconditioner::Estimator2(INMOST_DATA_ENUM_TYPE k,
															  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
															  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
															  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Est,
															  INMOST_DATA_REAL_TYPE & mu_update)
	{
		INMOST_DATA_REAL_TYPE mup = 1.0 - Est[k];
		INMOST_DATA_REAL_TYPE mum = -1.0 - Est[k];
		INMOST_DATA_ENUM_TYPE np = 0, nm = 0;
		//start from the element next after diagonal position
		INMOST_DATA_ENUM_TYPE i = LineIndeces[k];
		INMOST_DATA_REAL_TYPE v, vp, vm;
		while (i != EOL)
		{
			v = Est[i];
			vp = fabs(v + LineValues[i] * mup);
			vm = fabs(v + LineValues[i] * mum);
			v = fabs(v);
			if (vp > std::max<INMOST_DATA_REAL_TYPE>(2 * v, 0.5)) np++;
			if (vm > std::max<INMOST_DATA_REAL_TYPE>(2 * v, 0.5)) nm++;
			if (std::max<INMOST_DATA_REAL_TYPE>(2 * vp, 0.5) < v) np--;
			if (std::max<INMOST_DATA_REAL_TYPE>(2 * vm, 0.5) < v) nm--;
			i = LineIndeces[i];
		}
		if (np > nm) mu_update = mup; else mu_update = mum;
		return std::max(fabs(mup), fabs(mum));
	}
	
	void MLMTILUC_preconditioner::EstimatorUpdate(INMOST_DATA_ENUM_TYPE k,
												  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
												  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
												  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Est,
												  INMOST_DATA_REAL_TYPE & mu_update)
	{
		INMOST_DATA_ENUM_TYPE i = LineIndeces[k];
		while (i != EOL)
		{
			Est[i] += LineValues[i] * mu_update;
			i = LineIndeces[i];
		}
	}
	
	void MLMTILUC_preconditioner::DiagonalUpdate(INMOST_DATA_ENUM_TYPE k,
												 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Diag,
												 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndecesL,
												 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValuesL,
												 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndecesU,
												 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValuesU)
	{
		INMOST_DATA_ENUM_TYPE Li = LineIndecesL[k];
		INMOST_DATA_ENUM_TYPE Ui = LineIndecesU[k];
		while (Li != EOL && Ui != EOL)
		{
			if( Ui > Li ) Li = LineIndecesL[Li];
			else if( Li > Ui ) Ui = LineIndecesU[Ui];
			else
			{
				assert(Ui > k && Li > k && Ui == Li);
				//~ std::cout << "Diag[" << Ui << "] " << Diag[Ui] << " - " << LineValuesL[Li]*LineValuesU[Ui]*Diag[k] << " = ";
				Diag[Ui] -= LineValuesL[Li]*LineValuesU[Ui]*Diag[k];
				//~ std::cout << Diag[Ui] << " from " << Li << "," << Ui << std::endl;
				Li = LineIndecesL[Li];
				Ui = LineIndecesU[Ui];
			}
		}
	}
	
	void MLMTILUC_preconditioner::ClearList(INMOST_DATA_ENUM_TYPE beg,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues)
	{
		INMOST_DATA_ENUM_TYPE i = beg, j;
		while (i != EOL)
		{
			j = LineIndeces[i];
			LineValues[i] = 0.0; //clean values after use
			LineIndeces[i] = UNDEF; //clean indeces after use
			i = j;
		}
	}
	
	void MLMTILUC_preconditioner::PrepareTranspose(INMOST_DATA_ENUM_TYPE cbeg,
												   INMOST_DATA_ENUM_TYPE cend,
												   interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
												   std::vector<Sparse::Row::entry> & Entries,
												   interval<INMOST_DATA_ENUM_TYPE, std::vector< std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> > > & Indices )
	{
		for(INMOST_DATA_ENUM_TYPE i = cbeg; i < cend; ++i) //row-index
		{
			for(INMOST_DATA_ENUM_TYPE j = Address[i].first; j < Address[i].last; ++j) //entry index
				Indices[Entries[j].first].push_back( std::make_pair(i,j) ); // Entries[j].first is column index, record row-index, entry index
		}
	}

	void MLMTILUC_preconditioner::PrepareGraph(INMOST_DATA_ENUM_TYPE wbeg,
					  						   INMOST_DATA_ENUM_TYPE wend,
					  						   const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
					  						   const std::vector<Sparse::Row::entry> & Entries,
					  						   interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G)
	{
		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) 
		{
			G[k].reserve(Address[k].Size());
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it) 
				G[k].push_back(Entries[it].first);
		}
	}

	void MLMTILUC_preconditioner::PrepareGraphTranspose(INMOST_DATA_ENUM_TYPE wbeg,
							   							INMOST_DATA_ENUM_TYPE wend,
							   							const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G,
							   							interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG)
	{
		//preallocate tG???
		interval< INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > nnz(tG.get_interval_beg(),tG.get_interval_end(),0);
		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) 
			for (INMOST_DATA_ENUM_TYPE it = 0; it < G[k].size(); ++it) 
				nnz[G[k][it]]++;
		for(INMOST_DATA_ENUM_TYPE k = tG.get_interval_beg(); k < tG.get_interval_end(); ++k) 
			tG[k].reserve(nnz[k]);
		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) 
		{
			for (INMOST_DATA_ENUM_TYPE it = 0; it < G[k].size(); ++it) 
				tG[G[k][it]].push_back(k);
		}
	}

	void MLMTILUC_preconditioner::PrepareGraphProduct(INMOST_DATA_ENUM_TYPE wbeg,
													  INMOST_DATA_ENUM_TYPE wend,
							 						  const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G,
							 						  const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG,
							 						  interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & pG)
	{
#if defined(USE_OMP_FACT)
#pragma omp parallel
#endif
		{
			interval< INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > List(wbeg,wend,ENUMUNDEF);
#if defined(USE_OMP_FACT)
#pragma omp for
#endif
			for(INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
			{
				INMOST_DATA_ENUM_TYPE Beg = EOL, nnz = 0;
				// go over connection of k-th row
				std::swap(List[k] = k,Beg);
				for(INMOST_DATA_ENUM_TYPE it = 0; it < G[k].size(); ++it)
				{
					INMOST_DATA_ENUM_TYPE kt = G[k][it];
					for(INMOST_DATA_ENUM_TYPE jt = 0; jt < tG[kt].size(); ++jt)
					{
						INMOST_DATA_ENUM_TYPE lt = tG[kt][jt];
						//insert to linked list
						if( List[lt] == ENUMUNDEF )
							std::swap(List[lt] = lt,Beg), nnz++;
					}
				}
				// wadj_sep[k-wbeg] = nnz;
				pG[k].reserve(nnz);
				INMOST_DATA_ENUM_TYPE it = Beg, jt;
				while( it != static_cast<INMOST_DATA_ENUM_TYPE>(k) )
				{
					// if( it != k )
						pG[k].push_back(it);
					jt = List[it]; //save next entry
					List[it] = ENUMUNDEF; //clear list
					it = jt; //restore next
				}
				List[k] = ENUMUNDEF;
			}
		}
	}

	int MLMTILUC_preconditioner::Thread()
	{
#if defined(USE_OMP)
		return omp_get_thread_num();
#else
		return 0;
#endif
	}

	int MLMTILUC_preconditioner::Threads()
	{
#if defined(USE_OMP)
		return omp_get_max_threads();
#else
		return 1;
#endif
	}

	void MLMTILUC_preconditioner::ColumnInterval(INMOST_DATA_ENUM_TYPE wbeg, 
						 						 INMOST_DATA_ENUM_TYPE wend,
						 						 const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						 						 const std::vector<Sparse::Row::entry> & Entries,
						 						 INMOST_DATA_ENUM_TYPE & cbeg,
						 						 INMOST_DATA_ENUM_TYPE & cend)
	{
		cbeg = ENUMUNDEF, cend = 0;
		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) 
		{
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it) 
			{
				INMOST_DATA_ENUM_TYPE jt = Entries[it].first;
				if( jt > cend ) cend = jt;
				if( jt < cbeg ) cbeg = jt;
			}
		}
		cend++;
	}

	void MLMTILUC_preconditioner::FilterGraph(const Block & b,
											  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				 		  					  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
											  interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G_in,
											  interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G_out)
	{
		INMOST_DATA_ENUM_TYPE wbeg = b.row_start, wend = b.row_end;
		INMOST_DATA_ENUM_TYPE cbeg = b.col_start, cend = b.col_end;
		G_out.set_interval_beg(wbeg);
		G_out.set_interval_end(wend);
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for(INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
		{
			// std::swap(G_out[k],G_in[invP[k]]); //invP is where to get the row
			G_out[k] = G_in[invP[k]];
			INMOST_DATA_ENUM_TYPE it = 0, jt = 0;
			while(it < G_out[k].size())
			{
				G_out[k][jt] = localQ[G_out[k][it++]]; //localQ is the new column position
				if( cbeg <= G_out[k][jt] && G_out[k][jt] < cend ) jt++; //column is inside block
			}
			G_out[k].resize(jt);
		}
	}

	void MLMTILUC_preconditioner::FilterGraphTranspose(const Block & b,
											  		  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		  					  		  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ,
											  		  interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG_in,
											  		  interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG_out)
	{
		INMOST_DATA_ENUM_TYPE wbeg = b.row_start, wend = b.row_end;
		INMOST_DATA_ENUM_TYPE cbeg = b.col_start, cend = b.col_end;
		tG_out.set_interval_beg(cbeg);
		tG_out.set_interval_end(cend);
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for(INMOST_DATA_INTEGER_TYPE k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
		{
			// std::swap(tG_out[k],tG_in[invQ[k]]); //invQ is where to get the column
			tG_out[k] = tG_in[invQ[k]];
			INMOST_DATA_ENUM_TYPE it = 0, jt = 0;
			while(it < tG_out[k].size())
			{
				tG_out[k][jt] = localP[tG_out[k][it++]]; //localP is the new row position
				if( wbeg <= tG_out[k][jt] && tG_out[k][jt] < wend ) jt++; //row is inside block
			}
			tG_out[k].resize(jt);
		}
	}

	void MLMTILUC_preconditioner::FilterGraphProduct(const Block & b,
													 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				 		  							 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
													 interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & pG_in,
													 interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & pG_out)
	{
		INMOST_DATA_ENUM_TYPE wbeg = b.row_start, wend = b.row_end;
		pG_out.set_interval_beg(wbeg);
		pG_out.set_interval_end(wend);
#if defined(USE_OMP)
#pragma omp parallel for
#endif	
		for(INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
		{
			// std::swap(pG_out[k],pG_in[invP[k]]); //invP is where to get the row
			pG_out[k] = pG_in[invP[k]];
			INMOST_DATA_ENUM_TYPE it = 0, jt = 0;
			while(it < pG_out[k].size())
			{
				pG_out[k][jt] = localP[pG_out[k][it++]]; //localP is the new row position
				if( wbeg <= pG_out[k][jt] && pG_out[k][jt] < wend ) jt++; //row is inside block
			}
			pG_out[k].resize(jt);
		}
	}

	void MLMTILUC_preconditioner::DumpGraph(std::string name, interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G)
	{
		INMOST_DATA_ENUM_TYPE wbeg, wend, cbeg, cend, nnz = 0, side;
		wbeg = G.get_interval_beg();
		wend = G.get_interval_end();
		cbeg = ENUMUNDEF;
		cend = 0;
		std::cout << __FUNCTION__ << " " << name << std::endl;
		std::cout << "wbeg " << wbeg << " wend " << wend << std::endl;
		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
		{
			for(INMOST_DATA_ENUM_TYPE it = 0; it < G[k].size(); ++it)
			{
				INMOST_DATA_ENUM_TYPE jt = G[k][it];
				if( cbeg > jt ) cbeg = jt;
				if( cend < jt ) cend = jt;
				nnz++;
			}
		}
		std::cout << "cbeg " << cbeg << " cend " << cend << std::endl;
		side = std::max(wend-wbeg,cend-cbeg);
		std::ofstream file(name.c_str());
		file << "%%MatrixMarket matrix coordinate real general\n";
		file << "%graph only, values are units\n";
		file << "%writing as square graph, some rows/cols may be missing\n";
		file << side << " " << side << " " << nnz << std::endl;
		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
		{
			for(INMOST_DATA_ENUM_TYPE it = 0; it < G[k].size(); ++it)
			{
				INMOST_DATA_ENUM_TYPE jt = G[k][it];
				file << k-wbeg+1 << " " << jt-cbeg+1 << " 1\n";
			}
		}
		file.close();
	}

	void MLMTILUC_preconditioner::NestedDissection(INMOST_DATA_ENUM_TYPE wbeg,
												   INMOST_DATA_ENUM_TYPE wend,
												   const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						 						   const std::vector<Sparse::Row::entry> & Entries,
						 						   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		 						   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
												   std::vector<Block> & blocks,
												   INMOST_DATA_ENUM_TYPE max_size)
	{
		const int kway_parts = 2;
		double /*timer = Timer(),*/ total_time = Timer();
		INMOST_DATA_ENUM_TYPE cbeg, cend, sep, blks;
		ColumnInterval(wbeg,wend,Address,Entries,cbeg,cend);

		interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > G(wbeg,wend), tG(cbeg,cend), pG(wbeg,wend), G2, tG2, pG2;
		interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > * rG = &G, * rtG = &tG, * rpG = &pG;
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>  P(wbeg,wend), Q(cbeg,cend), invP(wbeg,wend), invQ(cbeg,cend);
		std::vector<Block> blocks_new;

		// std::cout << __FUNCTION__ << " row " << wbeg << ":" << wend << " col " << cbeg << ":" << cend << std::endl;
		
		//~ timer = Timer();

		PrepareGraph(wbeg,wend,Address,Entries,G);
		PrepareGraphTranspose(wbeg,wend,G,tG);

		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) localP[k] = invP[k] = k;
		for(INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k) localQ[k] = invQ[k] = k;
		
		// std::cout << "prepare tG time " << Timer() - timer << std::endl;
		// if( wgt_sep )
		{
			//~ timer = Timer();
			PrepareGraphProduct(wbeg,wend,G,tG,pG);
			// std::cout << "prepare pG time " << Timer() - timer << std::endl;
		}

		// DumpGraph("G.mtx",G);
		// DumpGraph("tG.mtx",tG);
		// DumpGraph("pG.mtx",pG);
		blocks.push_back(Block(wbeg,wend,cbeg,cend));

		//find largest non-separator block
		INMOST_DATA_ENUM_TYPE cur = 0;//, cnt = 0;
		do
		{
			// std::cout << "separate block " << cur << " rows " << blocks[cur].RowSize() << " cols " << blocks[cur].ColSize() << std::endl;
			GreedyDissection(blocks[cur],*rG,*rtG,*rpG,P,Q,blocks_new,kway_parts);//,1,0,0);
			// compute reordering in global P,Q, we need it to compute reordering in vector during solve phase
			for (INMOST_DATA_ENUM_TYPE k = blocks[cur].row_start; k < blocks[cur].row_end; ++k) localP[invP[k]] = P[k];
			for (INMOST_DATA_ENUM_TYPE k = blocks[cur].col_start; k < blocks[cur].col_end; ++k) localQ[invQ[k]] = Q[k];
			//replace blocks
			blocks.insert(blocks.erase(blocks.begin()+cur),blocks_new.begin(),blocks_new.end());
			blocks_new.clear();
			
			for(INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k) if( !blocks[k].separator )
			{
				if( cur == ENUMUNDEF || blocks[cur].RowSize() < blocks[k].RowSize() )
					cur = k;
			}
			// std::cout << "next selected block " << cur << " rows " << blocks[cur].RowSize() << " cols " << blocks[cur].ColSize() << std::endl;
			// break;
			//prepare reduced graphs
			//TODO: do not use entire initial graph!
			//inverse permutation

			for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) invP[localP[k]] = k;
			for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k) invQ[localQ[k]] = k;

			// if( transpose )
			// {
			// 	FilterGraph(blocks[cur],invP,localQ,tG,G2);
			// 	FilterGraphTranspose(blocks[cur],localP,invQ,G,tG2);
			// 	FilterGraphProduct(blocks[cur],invP,localP,pG,pG2);
			// }
			// else
			{
				FilterGraph(blocks[cur],invP,localQ,G,G2);
				FilterGraphTranspose(blocks[cur],localP,invQ,tG,tG2);
				FilterGraphProduct(blocks[cur],invP,localP,pG,pG2);
			}


			rG = &G2;
			rtG = &tG2;
			rpG = &pG2;
			// DumpGraph("rG.mtx",G2);
			// DumpGraph("rtG.mtx",tG2);
			// DumpGraph("rpG.mtx",pG2);
			//break;
			// scanf("%*c");
			// cnt++;
		}
		while( blocks[cur].RowSize() > max_size );//&& cnt < 2);

		blks = sep = 0;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k)
		{
			if( blocks[k].separator ) sep += blocks[k].RowSize();
			else blks++;
			std::cout << (blocks[k].separator?"separator":"block") << "[" << k << "] rows " << blocks[k].row_start << ":" << blocks[k].row_end << "(" << blocks[k].RowSize() << ") cols " << blocks[k].col_start << ":" << blocks[k].col_end << "(" << blocks[k].ColSize() << ")" << std::endl;
		}
		std::cout << "total separator " << sep << "/" << wend-wbeg << " blocks " << blks << std::endl;


		std::cout << __FUNCTION__ << " time " << Timer() - total_time << std::endl;

		std::ofstream file("blocks.txt");
		for(INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k)
			file << blocks[k].row_start << " " << blocks[k].row_end << " " << blocks[k].col_start << " " << blocks[k].col_end << std::endl;
		file.close();
	}

	void MLMTILUC_preconditioner::KwayDissection(INMOST_DATA_ENUM_TYPE wbeg,
												 INMOST_DATA_ENUM_TYPE wend,
												 const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						 						 const std::vector<Sparse::Row::entry> & Entries,
						 						 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		 						 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
												 std::vector<Block> & blocks, int parts)
	{
		double /*timer = Timer(),*/ total_time = Timer();
		INMOST_DATA_ENUM_TYPE cbeg, cend, sep, blks;
		ColumnInterval(wbeg,wend,Address,Entries,cbeg,cend);

		interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > G(wbeg,wend), tG(cbeg,cend), pG(wbeg,wend);
		
		// std::cout << __FUNCTION__ << " row " << wbeg << ":" << wend << " col " << cbeg << ":" << cend << std::endl;
		
		//~ timer = Timer();

		PrepareGraph(wbeg,wend,Address,Entries,G);
		PrepareGraphTranspose(wbeg,wend,G,tG);
		// std::cout << "prepare tG time " << Timer() - timer << std::endl;
		
		//~ timer = Timer();
		PrepareGraphProduct(wbeg,wend,G,tG,pG);
		// std::cout << "prepare pG time " << Timer() - timer << std::endl;
		
		GreedyDissection(Block(wbeg,wend,cbeg,cend),G,tG,pG,localP,localQ,blocks,parts);
			
		blks = sep = 0;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k)
		{
			if( blocks[k].separator ) sep += blocks[k].RowSize();
			else blks++;
			std::cout << (blocks[k].separator?"separator":"block") << "[" << k << "] rows " << blocks[k].row_start << ":" << blocks[k].row_end << "(" << blocks[k].RowSize() << ") cols " << blocks[k].col_start << ":" << blocks[k].col_end << "(" << blocks[k].ColSize() << ")" << std::endl;
			
		}
		std::cout << "total separator " << sep << " blocks " << blks << std::endl;


		std::cout << __FUNCTION__ << " time " << Timer() - total_time << std::endl;

		std::ofstream file("blocks.txt");
		for(INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k)
			file << blocks[k].row_start << " " << blocks[k].row_end << " " << blocks[k].col_start << " " << blocks[k].col_end << std::endl;
		file.close();
	}

	void MLMTILUC_preconditioner::GreedyDissection(const Block & b,
												   const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G,
						  						   const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG,
						  						   const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & pG,
						 						   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		 						   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
												   std::vector<Block> & blocks, int kway_parts)
	{
		// double timer = Timer(), total_time = Timer();
		
		// const bool kway = false;
		bool kway = (kway_parts > 1);
		// const int kway_parts = 4;
		const int upd_sep = 1, upd_blk = 1;
		const int wgt_sep = 1, wgt_blk = 0;

		//int nthreads = Threads();

		// std::cout << __FUNCTION__ <<  " wgt sep " << wgt_sep << " blk " << wgt_blk << " kway " << kway << std::endl;
		
		INMOST_DATA_ENUM_TYPE wbeg = b.row_start, wend = b.row_end, cbeg = b.col_start, cend = b.col_end;
		INMOST_DATA_ENUM_TYPE index_col = cbeg, index_row = wbeg, cur = ENUMUNDEF;//, nxt;
		INMOST_DATA_ENUM_TYPE Beg, Beg_upd, sep_size = 0, row_size = 0, col_size = 0, tot_sep_size = 0;//, nnz;
		INMOST_DATA_ENUM_TYPE index_row_beg = index_row, index_col_beg = index_col;
		INMOST_DATA_REAL_TYPE min_ratio = 1.0e+20, ratio, block_mean, balance_ratio, optsep_ratio;
		INMOST_DATA_ENUM_TYPE min_index = ENUMUNDEF;
		int had_col, had_sep, had_wgt_blk, had_wgt_sep;

		std::vector<int> wadj_blk(wend-wbeg,0), wadj_sep(wend-wbeg,0);
		BinaryHeapCustom<int, std::less<int> > q0(wend-wbeg);
		interval< INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > List(wbeg,wend,ENUMUNDEF);
		interval< INMOST_DATA_ENUM_TYPE, std::pair<int,INMOST_DATA_ENUM_TYPE> > List_upd(wbeg,wend,std::make_pair(0,ENUMUNDEF));
		std::vector<INMOST_DATA_ENUM_TYPE> visits;
		std::vector< std::pair<INMOST_DATA_ENUM_TYPE, bool> > upd;

		if( !kway ) visits.reserve(wend-wbeg);

		// std::cout << "alloc time " << Timer() - timer << std::endl;
		// timer = Timer();

		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) 
		{
			localP[k] = ENUMUNDEF;
			if( wgt_blk ) wadj_blk[k-wbeg] = (int)G[k].size(); //row connection wgt
			if( wgt_sep ) wadj_sep[k-wbeg] = (int)pG[k].size();
		}

		for(INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k) localQ[k] = ENUMUNDEF;


		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) 
			q0.PushHeap(k-wbeg,wgt_sep*wadj_sep[k-wbeg]+wgt_blk*wadj_blk[k-wbeg]); 

		// std::cout << "prepare heap time " << Timer() - timer << std::endl;
		// timer = Timer();
		// exit(-1);
		
		
		// cur = 0;
		// for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) 
		// 	if( wadj_sep[k-wbeg] > wadj_sep[cur-wbeg] ) cur = k;
		// srand(time(NULL));
		// cur = rand()%(wend-wbeg)+wbeg;
		// std::cout << "start at " << cur << std::endl;
		// goto start;
		// std::ofstream rec("rec.txt");
		Beg = EOL;
		Beg_upd = EOL;
		{
			while( !q0.Empty() ) 
			{
				//cur has minimal separator growth size
				cur = q0.PopHeap() + wbeg; //row index
			// start:
				if( localP[cur] != ENUMUNDEF ) continue; //skip separator of other blocks
				had_col = index_col;
				had_sep = sep_size;
				had_wgt_blk = wadj_blk[cur-wbeg]; 
				had_wgt_sep = wadj_sep[cur-wbeg];
				// rec << "(" << cur << "," << wadj_blk[cur-wbeg] << "," << wadj_sep[cur-wbeg] << ")";
				// rec << " row " << index_row << " col " << index_col << " sep " << sep_size;
				localP[cur] = index_row++;
				row_size++;
				if( List[cur] != ENUMUNDEF ) //uncount me from separator
					sep_size--;
				else //include me to linked list to avoid count as separator
				{
					std::swap(List[cur] = cur,Beg);
					//extract me from rows that depend on me
					if( wgt_sep ) upd.push_back(std::make_pair(cur,false));
				}
				// go over connection of enumerated row
				int col_skip = 0;
				int sep_skip = 0;
				for(INMOST_DATA_ENUM_TYPE it = 0; it < G[cur].size(); ++it)
				{
					INMOST_DATA_ENUM_TYPE kt = G[cur][it]; //column index
					if( localQ[kt] == ENUMUNDEF ) 
					{
						localQ[kt] = index_col++;
						col_size++;
						for(INMOST_DATA_ENUM_TYPE jt = 0; jt < tG[kt].size(); ++jt)
						{
							INMOST_DATA_ENUM_TYPE lt = tG[kt][jt]; //row index
							wadj_blk[lt-wbeg]-=upd_blk;
							// if( wgt_blk && q0.Contains(lt-wbeg) )
							// 	q0.ChangeKey(lt-wbeg,wgt_sep*wadj_sep[lt-wbeg]+wgt_blk*wadj_blk[lt-wbeg]);
							if( q0.Contains(lt-wbeg) )
							{
								List_upd[lt].first += wgt_blk*upd_blk;
								if( List_upd[lt].second == ENUMUNDEF )
									std::swap(List_upd[lt].second = lt, Beg_upd);
							}
							// add dependencies to linked-list of separator
							if( List[lt] == ENUMUNDEF ) //new row into separator
							{
								std::swap(List[lt] = lt, Beg);
								sep_size++;
								//extract me from rows that depend on me
								if( wgt_sep ) upd.push_back(std::make_pair(lt,true));
							}
							else sep_skip++;
						}
					}
					else col_skip++;
				}
				// collect updates
				for(INMOST_DATA_ENUM_TYPE it = 0; it < upd.size(); ++it)
				{
					INMOST_DATA_ENUM_TYPE lt = upd[it].first; //row index
					for(INMOST_DATA_ENUM_TYPE ut = 0; ut < pG[lt].size(); ++ut)
					{
						INMOST_DATA_ENUM_TYPE qt = pG[lt][ut]; //row index
						if( localP[qt] == ENUMUNDEF )
						{
							wadj_sep[qt-wbeg]-=upd_sep;
							if( q0.Contains(qt-wbeg) )
							{
								List_upd[qt].first += wgt_sep*upd_sep;
								if( List_upd[qt].second == ENUMUNDEF )
									std::swap(List_upd[qt].second = qt, Beg_upd);
							}
						}
					}
					if( upd[it].second )
					{
						wadj_sep[lt-wbeg]-=upd_sep; //taking me will reduce separator size
						if( q0.Contains(lt-wbeg) )
						{
							List_upd[lt].first += wgt_sep*upd_sep;
							if( List_upd[lt].second == ENUMUNDEF )
								std::swap(List_upd[lt].second = lt, Beg_upd);
						}
					}
				}
				// update weights
				{
					INMOST_DATA_ENUM_TYPE it = Beg_upd,jt;
					while(it != EOL)
					{
						q0.ChangeKey(it-wbeg, q0.GetKey(it-wbeg)-List_upd[it].first);
						jt = List_upd[it].second;
						List_upd[it].first = 0;
						List_upd[it].second = ENUMUNDEF;
						it = jt;
					}
					Beg_upd = EOL;
				}
				upd.clear();
				if( wgt_blk && upd_blk == 1 && had_wgt_blk != index_col - had_col )
					std::cout << "Wrong prediction on column growth: " << had_wgt_blk << " real growth "  << index_col - had_col  << std::endl;

				if( wgt_sep && upd_sep == 1 && had_wgt_sep != ((int)sep_size - had_sep) )
					std::cout << "Wrong prediction on separtor growth: " << had_wgt_sep << " real growth "  << ((int)sep_size - had_sep)  << std::endl;

				//disconnect each row that connects to me to know block growth size
				if( !kway )
				{
					visits.push_back(cur);
					// block_mean = sqrt(((INMOST_DATA_REAL_TYPE)row_size)*((INMOST_DATA_REAL_TYPE)wend-wbeg-row_size-sep_size));
					block_mean = ((INMOST_DATA_REAL_TYPE)std::min(row_size,wend-wbeg-row_size-sep_size));
					// balance_ratio = pow(0.5*((INMOST_DATA_REAL_TYPE)wend-wbeg-sep_size) / block_mean,5);
					balance_ratio = 1;
					optsep_ratio = ((INMOST_DATA_REAL_TYPE)sep_size) / block_mean;
					ratio = optsep_ratio * balance_ratio;
					if( ratio < min_ratio )
					{
						min_ratio = ratio;
						min_index = index_row;
					}
				}
				else
				{
					bool take_next = false;
					// if( !q0.Empty() )
					// {
					// 	nxt = q0.PeekHeap()+wbeg;
					// 	take_next = (wgt_sep && upd_sep == 1 && nxt != ENUMUNDEF && wadj_sep[nxt-wbeg] <= 0);
					// }
					if( (row_size > (wend-wbeg-tot_sep_size-sep_size)/kway_parts) && !take_next )
					// if( col_size > (cend-cbeg)/kway_parts && !take_next )
					// if( row_size > 5*sep_size && !take_next )
					{
						
						INMOST_DATA_ENUM_TYPE it = Beg;//, jt;
						while( it != EOL )
						{
							if( localP[it] == ENUMUNDEF ) 
							{
								localP[it] = ENUMUNDEF-1; //enumerate separator
								tot_sep_size++;
							}
							it = List[it];
						}
						// std::cout << "make new block! total separator " << tot_sep_size << " gained " << sep_size << " rows " << index_row-index_row_beg << " cols " << index_col-index_col_beg << std::endl;
						blocks.push_back(Block(index_row_beg,index_row,index_col_beg,index_col));
						index_row_beg = index_row;
						index_col_beg = index_col;
						col_size = row_size = sep_size = 0;
					}
				}
			}
		}
		// rec.close();
		if( !kway )
		{
			// std::cout << "simulation time " << Timer() - timer << std::endl;
			// timer = Timer();

			//clear linked-list
			{
				//interval< INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > & List = Lists[0];
				INMOST_DATA_ENUM_TYPE it = Beg, jt;
				while( it != EOL )
				{
					jt = List[it]; //save next entry
					List[it] = ENUMUNDEF; //clear list
					it = jt; //restore next
				}
				Beg = EOL;
			}
		
			for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) localP[k] = ENUMUNDEF;
			for(INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k) localQ[k] = ENUMUNDEF;

			index_row_beg = index_row = wbeg;
			index_col_beg = index_col = cbeg;
			//fill first block
			for(INMOST_DATA_ENUM_TYPE k = 0; k < min_index-wbeg; ++k)
			{
				//interval< INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > & List = Lists[0];
				cur = visits[k];
				localP[cur] = index_row++;
				for(INMOST_DATA_ENUM_TYPE it = 0; it < G[cur].size(); ++it)
				{
					INMOST_DATA_ENUM_TYPE kt = G[cur][it];
					if( localQ[kt] == ENUMUNDEF ) 
						localQ[kt] = index_col++;
					for(INMOST_DATA_ENUM_TYPE jt = 0; jt < tG[kt].size(); ++jt)
					{
						INMOST_DATA_ENUM_TYPE lt = tG[kt][jt];
						if( List[lt] == ENUMUNDEF ) //new row into separator
							std::swap(List[lt] = lt,Beg);
					}
				}
			}
			//mark separator
			tot_sep_size = 0;
			{
				//interval< INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > & List = Lists[0];
				INMOST_DATA_ENUM_TYPE it = Beg;//, jt;
				while( it != EOL )
				{
					if( localP[it] == ENUMUNDEF ) 
					{
						localP[it] = ENUMUNDEF-1; //enumerate separator
						tot_sep_size++;
					}
					it = List[it];
				}
			}
			if( index_row_beg != index_row )
			{
				// std::cout << "make new block! total separator " << tot_sep_size  << " rows " << index_row-index_row_beg << " cols " << index_col-index_col_beg << std::endl;
				blocks.push_back(Block(index_row_beg,index_row,index_col_beg,index_col));
				index_row_beg = index_row;
				index_col_beg = index_col;
			}
			//fill second block preserving order and skipping separator
			for(INMOST_DATA_ENUM_TYPE k = min_index-wbeg; k < wend-wbeg; ++k)
			{
				cur = visits[k];
				if( localP[cur] != ENUMUNDEF ) continue;
				localP[cur] = index_row++;
				for(INMOST_DATA_ENUM_TYPE it = 0; it < G[cur].size(); ++it)
				{
					INMOST_DATA_ENUM_TYPE kt = G[cur][it];
					if( localQ[kt] == ENUMUNDEF ) 
						localQ[kt] = index_col++;
				}
			}
		}
		//add last block
		if( index_row_beg != index_row )
		{
			// std::cout << "make new block! total separator " << tot_sep_size << " rows " << index_row-index_row_beg << " cols " << index_col-index_col_beg << std::endl;
			blocks.push_back(Block(index_row_beg,index_row,index_col_beg,index_col));
			index_row_beg = index_row;
			index_col_beg = index_col;
		}

		// std::cout << "blocks: " << blocks.size() << " total separator " << tot_sep_size << "/" << wend-wbeg << std::endl;

		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) if( localP[k] == ENUMUNDEF-1 ) localP[k] = index_row++;
		for(INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k) if( localQ[k] == ENUMUNDEF   ) localQ[k] = index_col++;

		if( index_row_beg != index_row )
		{
			// std::cout << "make separator block! total separator " << tot_sep_size << " rows " << index_row-index_row_beg << " extra cols " << index_col-index_col_beg << std::endl;
			blocks.push_back(Block(index_row_beg,index_row,cbeg,index_col,true));
		}

		//check blocks
#if 0
		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
		{
			INMOST_DATA_ENUM_TYPE r = localP[k]; //new k-th row number
			INMOST_DATA_ENUM_TYPE b = ENUMUNDEF; //block number
			for(INMOST_DATA_ENUM_TYPE it = 0; it < blocks.size(); ++it)
				if( blocks[it].row_start <= r && r < blocks[it].row_end )
				{
					b = it;
					break;
				}
			if( b == ENUMUNDEF )
			{
				std::cout << "cannot find block number for row " << k << " new position " << localP[k] << std::endl;
				continue;
			}
			for(INMOST_DATA_ENUM_TYPE it = 0; it < G[k].size(); ++it)
			{
				INMOST_DATA_ENUM_TYPE c = localQ[G[k][it]]; //new column position
				if( !(blocks[b].col_start <= c && c < blocks[b].col_end) )
				{
					std::cout << "position (" << k << "," << G[k][it] << ")";
					std::cout << " new (" << r << "," << c << ") outside block ";
					std::cout << "(" << blocks[b].row_start << ":" << blocks[b].row_end << ",";
					std::cout << blocks[b].col_start << ":" << blocks[b].col_end << ")" << std::endl;
				}
			}
		}
#endif
		// std::cout << "construction time " << Timer() - timer << std::endl;
		// std::cout << __FUNCTION__ << " time " << Timer() - total_time << std::endl;

		// exit(-1);
		
	}

	
	INMOST_DATA_ENUM_TYPE & MLMTILUC_preconditioner::EnumParameter(std::string name)
	{
		if (name == "scale_iters") return sciters;
		else if( name == "estimator" ) return estimator;
		else if( name == "verbosity" ) return verbosity;
		throw - 1;
	}
	INMOST_DATA_REAL_TYPE & MLMTILUC_preconditioner::RealParameter(std::string name)
	{
		if (name == "tau") return tau;
		else if( name == "tau2" ) return iluc2_tau;
		else if( name == "pivot_cond" ) return pivot_cond;
		else if( name == "pivot_diag" ) return pivot_diag;
		else if( name == "condition_number_L" ) return condestL;
		else if( name == "condition_number_U" ) return condestU;
		throw - 1;
	}
	void MLMTILUC_preconditioner::Copy(const Method * other)
	{
		const MLMTILUC_preconditioner * b = dynamic_cast<const MLMTILUC_preconditioner *>(other);
		assert(b != NULL);
		tau = b->tau;
		iluc2_tau = b->iluc2_tau;
		pivot_cond = b->pivot_cond;
		pivot_diag = b->pivot_diag;
		Alink = b->Alink;
		info = b->info;
		sciters = b->sciters;
		eps = b->eps;
	}
	MLMTILUC_preconditioner::MLMTILUC_preconditioner(const MLMTILUC_preconditioner & other) :Method(other)
	{
		Copy(&other);
	}
	MLMTILUC_preconditioner & MLMTILUC_preconditioner::operator =(MLMTILUC_preconditioner const & other)
	{
		Copy(&other);
		return *this;
	}
	MLMTILUC_preconditioner::MLMTILUC_preconditioner(Solver::OrderInfo & info)
		:tau(DEFAULT_TAU),info(&info)
	{
		condestL = condestU = 1;
		verbosity = 0;
		Alink = NULL;
		init = false;
		sciters = 8;
		eps = 1e-54;
		estimator = 1;
		tau = 1.0e-3;
		iluc2_tau = tau*tau;
		pivot_cond = PIVOT_COND_DEFAULT;
		pivot_diag = PIVOT_DIAG_DEFAULT;
	}
	bool MLMTILUC_preconditioner::isInitialized() { return init; }
	bool MLMTILUC_preconditioner::isFinalized() { return !init; }
	bool MLMTILUC_preconditioner::Initialize()
	{
		if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
		int nthreads = Threads();
        const INMOST_DATA_REAL_TYPE subst = 1.0; (void)subst;
		const INMOST_DATA_REAL_TYPE tol_modif = PIVOT_THRESHOLD_VALUE;
		
		bool block_pivot = false;
		
		INMOST_DATA_ENUM_TYPE wbeg, wend; //working interval
		INMOST_DATA_ENUM_TYPE mobeg, moend; // total interval
		INMOST_DATA_ENUM_TYPE vbeg, vend; // vector interval
		
		INMOST_DATA_INTEGER_TYPE k;
		INMOST_DATA_ENUM_TYPE i, j, Li, Ui, curr, next;
		INMOST_DATA_REAL_TYPE l,u, max_diag = 0, min_diag = 0;
		INMOST_DATA_ENUM_TYPE nzA, nzLU = 0, nzA0;
		Sparse::Vector DL, DR;
		Sparse::Vector DL0,DR0;
		info->GetOverlapRegion(info->GetRank(), mobeg, moend);
		info->GetVectorRegion(vbeg,vend);
		
		

		
		//prepare reordering vectors
		ddP.set_interval_beg(mobeg);
		ddP.set_interval_end(moend);
		ddQ.set_interval_beg(mobeg);
		ddQ.set_interval_end(moend);

		

		//prepare rescaling vectors
		DL.SetInterval(mobeg, moend);
		DR.SetInterval(mobeg, moend);
		DL0.SetInterval(mobeg, moend);
		DR0.SetInterval(mobeg, moend);
		for(k = mobeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(moend); ++k) 
			DL[k] = DR[k] = DL0[k] = DR0[k] = 1.0;
		
		

		for (k = mobeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(moend); k++) 
			ddP[k] = ddQ[k] = k;
		// supplementary data for ddPQ reordering
		interval<INMOST_DATA_ENUM_TYPE, bool> donePQ(mobeg, moend, false);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> invP(mobeg, moend), invQ(mobeg, moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> localP(mobeg, moend, ENUMUNDEF), localQ(mobeg, moend, ENUMUNDEF);
		
		interval<INMOST_DATA_ENUM_TYPE, bool> Pivot(mobeg,moend,false);

		if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

		//supplimentary data structures for ILUC
		INMOST_DATA_ENUM_TYPE L_Beg, U_Beg, Sbeg;
		U_Address.set_interval_beg(mobeg);
		U_Address.set_interval_end(moend);
		L_Address.set_interval_beg(mobeg);
		L_Address.set_interval_end(moend);
		LU_Diag.set_interval_beg(mobeg);
		LU_Diag.set_interval_end(moend);
		std::vector<Sparse::Row::entry> A_Entries;
		
		
		//std::vector<INMOST_DATA_ENUM_TYPE> sort_indeces;
		interval<INMOST_DATA_ENUM_TYPE, Interval> A_Address(mobeg, moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ulist(mobeg, moend), Llist(mobeg, moend), Blist(mobeg, moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Bbeg(mobeg,moend,EOL);
#if defined(NNZ_PIVOT)
		INMOST_DATA_ENUM_TYPE nnzL = 0, nnzU = 0, nnzL_max = 0, nnzU_max = 0;
#endif
		INMOST_DATA_ENUM_TYPE nnzE = 0, nnzF = 0, nnzA = 0;
		INMOST_DATA_REAL_TYPE NuU = 1, NuL = 1, NuD = 1, NuU_max = 1.0, NuL_max = 1.0;
		INMOST_DATA_REAL_TYPE NuU_acc = 1, NuL_acc = 1, NuD_acc = 1;
		
		//supplimentary data structures for returning values of dropped elements
		//~ INMOST_DATA_REAL_TYPE DropLk, DropUk, SumLk, SumUk;
		//~ interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> DropU(mobeg,moend,0.0), DropL(mobeg,moend,0.0);
//~ #if defined(ILUC2)
		//~ interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> DropU2(mobeg,moend,0.0), DropL2(mobeg,moend,0.0);
//~ #endif
		//data structure for linked list

		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> U(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> V(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		//interval<INMOST_DATA_ENUM_TYPE, std::vector< std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> > > TransposeU(mobeg,moend), TransposeL(mobeg,moend);
		
		
		
        double tfactor = 0.0, trescale = 0.0, treorder = 0.0, ttransversal = 0.0, treassamble = 0.0, ttotal, tt, testimator = 0.0, tschur = 0.0, tlocal;
#if defined(REORDER_METIS_ND)
		double tmetisgraph = 0, tmetisnd = 0;
#endif
#if defined(REORDER_RCM) || defined(REORDER_WRCM) || defined(REORDER_BRCM)
		double trcmgraph = 0, trcmorder = 0;
#endif
        double tlfactor, tlrescale, tlreorder, tlreassamble, tlschur;
		ttotal = Timer();

		//(*Alink).Save("M.mtx");
		
		//calculate number of nonzeros
		nzA = 0;
		for (k = mobeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(moend); ++k)
		{
			for (Sparse::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
				if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0) nzA++;
		}
		nzA0 = nzA;

		//sort_indeces.reserve(256);
		A_Entries.resize(nzA);
		//LU_Entries.reserve(nzA*4);

		j = 0;
		for (k = mobeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(moend); ++k)
		{
			A_Address[k].first = j;
			for (Sparse::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
			{
				if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0 )//1.0e-30 )// != 0.0)
					A_Entries[j++] = Sparse::Row::make_entry(r->first, r->second);
			}
			A_Address[k].last = j;
			//assert(A_Address[k].Size() != 0); //singular matrix
		}

		if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
		
		INMOST_DATA_REAL_TYPE Unorm, Lnorm;
		//~ INMOST_DATA_REAL_TYPE Unum, Lnum;

		INMOST_DATA_ENUM_TYPE nzLU2 = 0, nzLU2tot = 0, ndrops = 0;
#if defined(ILUC2)
		
		INMOST_DATA_REAL_TYPE tau2 = iluc2_tau;
		std::vector<Sparse::Row::entry> L2_Entries, U2_Entries;
		interval<INMOST_DATA_ENUM_TYPE, Interval> L2_Address(mobeg, moend+1), U2_Address(mobeg, moend+1);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2list(mobeg, moend, UNDEF), L2list(mobeg, moend, UNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2beg(mobeg, moend, EOL), L2beg(mobeg, moend, EOL);
		//LU2_Entries.reserve(nzA*4);
#else
		INMOST_DATA_REAL_TYPE tau2 = tau;
#endif

		for(k = mobeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(moend); ++k)
		{
			U_Address[k].first = U_Address[k].last = ENUMUNDEF;
			L_Address[k].first = L_Address[k].last = ENUMUNDEF;
#if defined(ILUC2)
			U2_Address[k].first = U2_Address[k].last = ENUMUNDEF;
			L2_Address[k].first = L2_Address[k].last = ENUMUNDEF;
#endif
		}

		if( verbosity > 1 ) 
			std::cout << "nonzeros in A " << nzA << " pivot cond " << pivot_cond << " diag " << pivot_diag << " tau " << tau << " tau2 " << tau2 << std::endl;

		
		int swaps = 0, totswaps = 0;
		//set up working interval
		wbeg = mobeg;
		wend = moend;
		while (wbeg < wend) //iterate into levels until all matrix is factored
		{
			//ddPQ reordering on current Schur complement
			INMOST_DATA_ENUM_TYPE cbeg = wbeg, cend = wend; //next size of factored B block

			if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
			
		//DumpMatrix(A_Address,A_Entries,cbeg,cend,"mat"+to_string(level_size.size())+".mtx");
///////////////////////////////////////////////////////////////////////////////////
///   MAXIMUM TRANSVERSE REORDERING                                             ///
///////////////////////////////////////////////////////////////////////////////////
			if( run_mpt )
			{
				if( verbosity > 1 )
					std::cout << "Reordering with MPT\n";

				ttransversal = Timer();
				INMOST_DATA_ENUM_TYPE ColumnBegin;
				INMOST_DATA_ENUM_TYPE pop_heap_pos;
				INMOST_DATA_REAL_TYPE ShortestPath, AugmentPath;
				INMOST_DATA_ENUM_TYPE PathEnd, Trace, IPermPrev;
				//Sparse::Vector & U = DL;
				//Sparse::Vector & V = DR;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Cmax(wbeg,wend,0.0);
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & Perm = localQ;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & IPerm = localP;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & ColumnList = Llist;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & Parent = Blist;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> AugmentPosition(wbeg,wend,ENUMUNDEF);
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> ColumnPosition(wbeg,wend,ENUMUNDEF);
				std::vector<INMOST_DATA_REAL_TYPE> C_Entries(A_Entries.size());
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Dist(wbeg, wend + 1, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				
				BinaryHeap Heap(&Dist[wbeg],wend-wbeg);
				// arrays U,V,Dist are cleared at the end of schur complement calculation
				//std::fill(U.begin() + wbeg - mobeg, U.begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				//std::fill(V.begin() + wbeg - mobeg, V.begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				//std::fill(Cmax.begin() + wbeg - mobeg, Cmax.begin() + wend - mobeg, 0.0);
				//std::fill(Dist.begin() + wbeg - mobeg, Dist.begin() + wend + 1 - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					Perm[k] = ENUMUNDEF;
					IPerm[k] = ENUMUNDEF;
					Parent[k] = ENUMUNDEF;
					ColumnList[k] = ENUMUNDEF;
				}

				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
				
				// Initial LOG transformation to dual problem and initial extreme match 
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
					{
						i = A_Entries[it].first;
						u = C_Entries[it] = fabs(A_Entries[it].second);
						if( u > Cmax[i] ) Cmax[i] = u;
						//C_Entries.push_back(Sparse::Row::make_entry(i,u));
					}
				}

				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
					{
						i = A_Entries[it].first;
						if( Cmax[i] == 0 || C_Entries[it] == 0 )
							C_Entries[it] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
						else
						{
							C_Entries[it] = log(Cmax[i])-log(C_Entries[it]);
							if( C_Entries[it] < U[i] ) U[i] = C_Entries[it];
						}
					}
				}
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
					{
						u = C_Entries[it] - U[A_Entries[it].first];
						if( u < V[k] ) V[k] = u;
					}
				}
				/// Update cost and match
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
					{
						u = fabs(C_Entries[it] - V[k] - U[A_Entries[it].first]);
						if( u < 1.0e-30 && Perm[A_Entries[it].first] == ENUMUNDEF && IPerm[k] == ENUMUNDEF )
						{
							 Perm[A_Entries[it].first] = k;
							 IPerm[k] = A_Entries[it].first;
							 ColumnPosition[k] = it;
						}
					}
				}
				/// 1-step augmentation
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					if( IPerm[k] == ENUMUNDEF ) //unmatched row
					{
						for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last && IPerm[k] == ENUMUNDEF; ++it)
						{
							u = fabs(C_Entries[it] - V[k] - U[A_Entries[it].first]);
							if( u <= 1.0e-30 )
							{
								Li = Perm[A_Entries[it].first];
								assert(Li != ENUMUNDEF);
								// Search other row in C for 0
								for (INMOST_DATA_ENUM_TYPE Lit = A_Address[Li].first; Lit < A_Address[Li].last; ++Lit)
								{
									u = fabs(C_Entries[Lit]- V[Li] - U[A_Entries[Lit].first]);
									if( u <= 1.0e-30 && Perm[A_Entries[Lit].first] == ENUMUNDEF )
									{
										Perm[A_Entries[it].first] = k;
										IPerm[k] = A_Entries[it].first;
										ColumnPosition[k] = it;
										Perm[A_Entries[Lit].first] = Li;
										IPerm[Li] = A_Entries[Lit].first;
										ColumnPosition[Li] = Lit;
										break;
									}
								}
							}
						}
					}
				}
				/// Weighted bipartite matching
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					if( IPerm[k] != ENUMUNDEF )
						continue;
					Li = k;
					ColumnBegin = EOL;
					Parent[Li] = ENUMUNDEF;
					PathEnd = ENUMUNDEF;
					Trace = k;
					ShortestPath = 0;
					AugmentPath = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
					while(true)
					{
						for (INMOST_DATA_ENUM_TYPE Lit = A_Address[Li].first; Lit < A_Address[Li].last; ++Lit)
						{
							Ui = A_Entries[Lit].first;
							//if( ColumnList[Ui] == k ) continue;
							if( ColumnList[Ui] != ENUMUNDEF ) continue;
							l = fabs(ShortestPath + C_Entries[Lit] - V[Li] - U[Ui]);
							//if( l < 0.0 ) printf("row %d col %d negative l %g Augment %lf Shortest %lf C %lf V %lf U %lf\n",k,Ui,l,AugmentPath,ShortestPath,C_Entries[Lit],V[Li],U[Ui]);
							if( l < 0.0 && l > -1.0e-8 ) l = 0;
							if( l < 0.0 ) continue;
							if( l < AugmentPath )
							{
								if( Perm[Ui] == ENUMUNDEF )
								{
									PathEnd = Ui;
									Trace = Li;
									AugmentPath = l;
									AugmentPosition[Ui] = Lit;
								}
								else if( l < Dist[Ui] )
								{
									Parent[Perm[Ui]] = Li;
									AugmentPosition[Ui] = Lit;
									if( Heap.Contains(Ui-wbeg) )
										Heap.DecreaseKey(Ui-wbeg,l);
									else
										Heap.PushHeap(Ui-wbeg,l);
								}
							}
						}

						pop_heap_pos = Heap.PopHeap();
						if( pop_heap_pos == ENUMUNDEF ) break;
					
						Ui = pop_heap_pos+wbeg;
						ShortestPath = Dist[Ui];

						if( AugmentPath <= ShortestPath ) 
						{
							Dist[Ui] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
							//Heap.increaseKey(Ui,Dist[Ui]);
							break;
						}


						ColumnList[Ui] = ColumnBegin;
						ColumnBegin = Ui;

						Li = Perm[Ui];
						
					}
					if( PathEnd != ENUMUNDEF )
					{
						Ui = ColumnBegin;
						while(Ui != EOL)
						{
							U[Ui] += Dist[Ui] - AugmentPath;
							if( Perm[Ui] != ENUMUNDEF ) V[Perm[Ui]] = C_Entries[ColumnPosition[Perm[Ui]]] - U[Ui];
							Dist[Ui] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
							Li = ColumnList[Ui];
							ColumnList[Ui] = ENUMUNDEF;
							Ui = Li;
						}

						Ui = PathEnd;
						while(Trace != ENUMUNDEF)
						{
							IPermPrev = IPerm[Trace];
							Perm[Ui] = Trace;
							IPerm[Trace] = Ui;

							ColumnPosition[Trace] = AugmentPosition[Ui];
							V[Trace] = C_Entries[ColumnPosition[Trace]] - U[Ui];

							Ui = IPermPrev;
							Trace = Parent[Trace];

						}
						Heap.Clear();
					}
				}
#if defined(USE_OMP_FACT)
#pragma omp parallel for private(k,l,u,i,j)
#endif
				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
				{
					if( V[k] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() ) l = 1;
					else l = exp(V[k]);
					if( U[k] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() || Cmax[k] == 0 ) u = 1;
					else u = exp(U[k])/Cmax[k];
					DL[k] = l;
					DR[k] = u;
					

					bool flip_sign = false;
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt < A_Address[k].last; ++jt)
					{
						i = A_Entries[jt].first;
						j = Perm[A_Entries[jt].first];
						
						if( U[i] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() || Cmax[i] == 0 ) u = 1;
						else u = exp(U[i])/Cmax[i];
						
						if( fabs(l*A_Entries[jt].second*u) > 1 + 1.0e-7 )
						{
							std::cout << "element on row " << k << " col " << A_Entries[jt].first << " value " << A_Entries[jt].second << " u " << u << " l " << l << " U " << U[i] << " V " << V[k] << " Cmax " << Cmax[i] << " scaled " << l*A_Entries[jt].second*u << std::endl;
							exit(-1);
						}
						if( static_cast<INMOST_DATA_INTEGER_TYPE>(j) == k )
						{
							
							if( l*A_Entries[jt].second*u < 0.0 ) flip_sign = true;
						}
					}

					if( flip_sign ) DL[k] *= -1;
				}
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k) localP[k] = ENUMUNDEF;
                
                { //check that there are no gaps in Perm
                    for(k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
                    {
                        if( Perm[k] != ENUMUNDEF )
						{
							assert(localP[Perm[k]] == ENUMUNDEF);
                            localP[Perm[k]] = 0;
						}
                    }
					std::vector<INMOST_DATA_ENUM_TYPE> gaps;
                    for(k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
                        if( localP[k] == ENUMUNDEF )
                            gaps.push_back(k);
					
					//~ std::cout << "@ gaps: " << gaps.size() << std::endl;
                    
                    for(k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
                        if( Perm[k] == ENUMUNDEF )
                        {
                            Perm[k] = gaps.back();
                            gaps.pop_back();
                        }                        
                    for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k) localP[k] = ENUMUNDEF;
                }
				
				ttransversal = Timer() - ttransversal;

				treorder += ttransversal;
				if( verbosity > 1 )
					std::cout << "Time " << ttransversal << "\n";

			}
			else
			{
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					localQ[k] = k;
					DL[k] = DR[k] = 1;
				}
			}
			///  END MAXIMUM TRANSVERSE REORDERING
			tt = Timer();
#if defined(REORDER_METIS_ND)
			idx_t nvtxs = wend-wbeg;
			std::vector<idx_t> xadj(nvtxs+1), adjncy, perm(nvtxs),iperm(nvtxs);
			//adjncy.reserve(nzA*2);
			if( verbosity > 1 )
				std::cout << "Reordering with METIS\n";

			
			for (k = cbeg; k < cend; ++k) if( localQ[k] == ENUMUNDEF ) std::cout << __FILE__ << ":" << __LINE__ << " No column permutation for row " << k << ". Matrix is structurally singular\n";

			for (k = cbeg; k < cend; ++k) localP[k] = k;
			for (k = cbeg; k < cend; ++k) V[localQ[k]] = DR[k]; //if you expirience a bug here then the matrix is structurally singular
			for (k = cbeg; k < cend; ++k) DR[k] = V[k];
			
			for (k = cbeg; k < cend; ++k)
			{
				if( !(localQ[k] >= cbeg && localQ[k] < cend) )
					std::cout << "Bad permutation: " << localQ[k] << " interval [" << cbeg << ":" << cend << "]" << std::endl;
				V[localQ[k]] = DR0[k]; //if you expirience a bug here then the matrix is structurally singular
			}
			for (k = cbeg; k < cend; ++k) DR0[k] = V[k];

			
			
			for(k = wbeg; k < wend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt < A_Address[k].last; ++jt)
				{
					A_Entries[jt].first = localQ[A_Entries[jt].first];
					//A_Entries[jt].second *= DL[k]*DR[A_Entries[jt].first];
					
				}
				std::sort(A_Entries.begin()+A_Address[k].first,A_Entries.begin()+A_Address[k].last);
			}
			ReorderEF(wbeg, wend, donePQ, localP, localQ);
			inversePQ(wbeg,wend,localP,localQ, invP,invQ);
			applyPQ(wbeg, wend, localP, localQ, invP, invQ);

			//DumpMatrix(A_Address,A_Entries,wbeg,wend,"a.mtx");

///////////////////////////////////////////////////////////////////////////////////
//       setup indices for transposed traversal of B block                       //
///////////////////////////////////////////////////////////////////////////////////
			
			tmetisgraph = Timer();

			for (k = wend; k > wbeg; --k)
			{
				//vwgt[k-1] = 1;
				if (A_Address[k-1].Size() > 0)
				{
					i = A_Entries[A_Address[k-1].first].first;
					Blist[k-1] = Bbeg[i];
					Bbeg[i] = k-1;
				}
				Ulist[k-1] = A_Address[k-1].first;
			}
			

			xadj[0] = 0;
			for(i = wbeg; i < wend; ++i)
			{
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[i].first; jt < A_Address[i].last; ++jt)
				{
					if( A_Entries[jt].first != i )
						adjncy.push_back(static_cast<idx_t>(A_Entries[jt].first-wbeg));
				}
				
				Li = Bbeg[i];
				while (Li != EOL)
				{
					if( Li != i )
						adjncy.push_back(static_cast<idx_t>(Li-wbeg));
					Li = Blist[Li];
				}
				
				Li = Bbeg[i];
				while (Li != EOL)
				{
					Ui = Blist[Li];
					Ulist[Li]++;
					if (A_Address[Li].last - Ulist[Li] > 0)
					{
						k = A_Entries[Ulist[Li]].first;
						Blist[Li] = Bbeg[k];
						Bbeg[k] = Li;
					}

					Li = Ui;
				}
				
				std::sort(adjncy.begin()+xadj[i-wbeg],adjncy.end());
				adjncy.resize(std::unique(adjncy.begin()+xadj[i-wbeg],adjncy.end())-adjncy.begin());
				
				xadj[i-wbeg+1] = static_cast<idx_t>(adjncy.size());
			}

			//std::fill(Bbeg.begin(),Bbeg.end(),EOL);
			for(i = mobeg; i < moend; ++i) Bbeg[i] = EOL;

			tmetisgraph = Timer() - tmetisgraph;
			/*
			{
				FILE * f = fopen("metis.mtx","w");
				fprintf(f,"%%MatrixMarket matrix coordinate real general\n");
				fprintf(f,"%d %d %d\n",nvtxs,nvtxs,adjncy.size());
				for(i = wbeg; i < wend; ++i)
				{
					for(j = xadj[i]; j < xadj[i+1]; ++j)
						fprintf(f,"%d %d 1\n",i+1,adjncy[j]+1);
				}
				fclose(f);
			}
			*/
			//A_Address[0].first = 0;
			//for(k = wbeg+1; k < wend; ++k)
			//	A_Address[k].first = A_Address[k-1].last;
			if( !adjncy.empty() )
			{
				idx_t options[METIS_NOPTIONS];
				METIS_SetDefaultOptions(options);
				options[METIS_OPTION_NUMBERING] = 0;
				
				//options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO | METIS_DBG_TIME | METIS_DBG_COARSEN | METIS_DBG_CONNINFO | METIS_DBG_CONTIGINFO | METIS_DBG_IPART | METIS_DBG_MEMORY | METIS_DBG_MOVEINFO | METIS_DBG_REFINE | METIS_DBG_SEPINFO;
				//options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM;
				//options[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP2SIDED;
				//options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE;
				//options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
				//options[METIS_OPTION_NO2HOP] = 0;
				//options[METIS_OPTION_MINCONN] = 0;
				//options[METIS_OPTION_NITER] = 4;
				//METIS_NodeNDP(nvtxs,&xadj[0],&adjncy[0],NULL,4,options,&perm[0],&iperm[0],&sizes[0]);
				tmetisnd = Timer();
				METIS_NodeND(&nvtxs,&xadj[0],&adjncy[0],NULL,options,&perm[0],&iperm[0]);
				tmetisnd = Timer()-tmetisnd;
				
				for(k = wbeg; k < wend; ++k)
				{
					//localP[k] = iperm[k-wbeg]+wbeg;
					//localQ[k] = iperm[k-wbeg]+wbeg;
					localP[perm[k-wbeg]+wbeg] = k;
					localQ[perm[k-wbeg]+wbeg] = k;
				}
			}
			else
			{
				for(k = wbeg; k < wend; ++k)
				{
					localP[k] = k;
					localQ[k] = k;
				}
			}
			cend = wend;
			i = wend;
#elif defined(REORDER_RCM)
			{
				if( verbosity > 1 ) 
					std::cout << "Reordering with RCM\n";
				//create a symmetric graph of the matrix A + A^T
				std::vector<INMOST_DATA_ENUM_TYPE> xadj(wend-wbeg+1), adjncy;
				//adjncy.reserve(nzA*2);

				for (k = cbeg; k < cend; ++k) if( localQ[k] == ENUMUNDEF ) std::cout << __FILE__ << ":" << __LINE__ << " No column permutation for row " << k << ". Matrix is structurally singular\n";

				for (k = cbeg; k < cend; ++k) localP[k] = k;
				for (k = cbeg; k < cend; ++k) V[localQ[k]] = DR[k]; //if you expirience a bug here then the matrix is structurally singular
				for (k = cbeg; k < cend; ++k) DR[k] = V[k];
				
				for (k = cbeg; k < cend; ++k) V[localQ[k]] = DR0[k]; //if you expirience a bug here then the matrix is structurally singular
				for (k = cbeg; k < cend; ++k) DR0[k] = V[k];

				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt < A_Address[k].last; ++jt)
					{
						A_Entries[jt].first = localQ[A_Entries[jt].first];
						//A_Entries[jt].second *= DL[k]*DR[A_Entries[jt].first];
					}
					std::sort(A_Entries.begin()+A_Address[k].first,A_Entries.begin()+A_Address[k].last);
				}
				
				//DumpMatrix(A_Address,A_Entries,cbeg,cend,"A.mtx");
				
				ReorderEF(wbeg, wend, donePQ, localP, localQ);
				inversePQ(wbeg,wend,localP,localQ, invP,invQ);
				applyPQ(wbeg, wend, localP, localQ, invP, invQ);

				trcmgraph = Timer();
				
				for(k = wbeg; k < wend; ++k)
				{
					Ulist[k] = EOL;
					Bbeg[k] = EOL;
					Blist[k] = EOL;
				}
				//std::fill(Ulist.begin()+wbeg-mobeg, Ulist.begin()+wend-mobeg,EOL);
				//std::fill(Bbeg.begin()+wbeg-mobeg, Bbeg.begin()+wend-mobeg,EOL);
				//std::fill(Blist.begin()+wbeg-mobeg, Blist.begin()+wend-mobeg,EOL);
				for (k = wend; k > wbeg; --k)
				{
					//vwgt[k-1] = 1;
					if (A_Address[k-1].Size() > 0)
					{
						i = A_Entries[A_Address[k-1].first].first;
						Blist[k-1] = Bbeg[i];
						Bbeg[i] = k-1;
					}
					Ulist[k-1] = A_Address[k-1].first;
				}
				xadj[0] = 0;
				for(i = wbeg; i < wend; ++i)
				{
					/*
					double Arowmax = 0, Acolmax = 0;
					
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[i].first; jt < A_Address[i].last; ++jt)
						if( A_Entries[jt].first != i )
							Arowmax = std::max(Arowmax,fabs(DL[i]*A_Entries[jt].second*DR[A_Entries[jt].first]));
						
					Li = Bbeg[i];
					while (Li != EOL)
					{
						if( Li != i )
							Acolmax = std::max(Acolmax,DL[A_Entries[A_Address[Li].first].first]*A_Entries[A_Address[Li].first].second*DR[i]);
						Li = Blist[Li];
					}
					*/
					
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[i].first; jt < A_Address[i].last; ++jt)
					{
						if( A_Entries[jt].first != i  )//&& fabs(DL[i]*A_Entries[jt].second*DR[A_Entries[jt].first]) > 0.1*std::min(Acolmax,Arowmax) )
							adjncy.push_back(static_cast<INMOST_DATA_ENUM_TYPE>(A_Entries[jt].first));
					}

					Li = Bbeg[i];
					while (Li != EOL)
					{
						if( Li != i )//&& fabs(DL[A_Entries[A_Address[Li].first].first]*A_Entries[A_Address[Li].first].second*DR[i]) > 0.1*std::min(Acolmax,Arowmax) )
							adjncy.push_back(static_cast<INMOST_DATA_ENUM_TYPE>(Li));
						Li = Blist[Li];
					}

					Li = Bbeg[i];
					while (Li != EOL)
					{
						Ui = Blist[Li];
						Ulist[Li]++;
						if (A_Address[Li].Size() && A_Address[Li].last - Ulist[Li] > 0)
						{
							k = A_Entries[Ulist[Li]].first;
							Blist[Li] = Bbeg[k];
							Bbeg[k] = Li;
						}

						Li = Ui;
					}

					std::sort(adjncy.begin()+xadj[i-wbeg],adjncy.end());
					adjncy.resize(std::unique(adjncy.begin()+xadj[i-wbeg],adjncy.end())-adjncy.begin());

					xadj[i-wbeg+1] = static_cast<INMOST_DATA_ENUM_TYPE>(adjncy.size());
				}

				//std::fill(Bbeg.begin()+wbeg - mobeg, Bbeg.begin()+wend - mobeg,EOL);
				for(k = wbeg; k < wend; ++k) Bbeg[k] = EOL;

				trcmgraph = Timer()-trcmgraph;

				trcmorder = Timer();
				//std::fill(Ulist.begin() + wbeg - mobeg, Ulist.begin() + wend - mobeg, ENUMUNDEF);
				for(k = wbeg; k < wend; ++k) Ulist[k] = ENUMUNDEF;
				//find node with the lowest order
                INMOST_DATA_ENUM_TYPE index = wbeg;
				INMOST_DATA_ENUM_TYPE cur = ENUMUNDEF;
				std::deque<INMOST_DATA_ENUM_TYPE> q;
				std::vector<INMOST_DATA_ENUM_TYPE> conns;
                do
				{
					cur = ENUMUNDEF;
					for(k = wbeg; k < wend && cur == ENUMUNDEF; ++k)
					{
						if( Ulist[k] == ENUMUNDEF )
							cur = k;
					}
					assert(cur != ENUMUNDEF);
					for(k = cur+1; k < wend; ++k) if( Ulist[k] == ENUMUNDEF )
					{
						if( RCM_Comparator(wbeg,xadj)(k,cur) )
							cur = k;
					}
					q.push_back(cur);
					Ulist[cur] = index++;
					while(!q.empty())
					{
						cur = q.front();
						q.pop_front();
						for (INMOST_DATA_ENUM_TYPE it = xadj[cur-wbeg]; it < xadj[cur-wbeg+1]; ++it)
							if( Ulist[adjncy[it]] == ENUMUNDEF )
								conns.push_back(adjncy[it]);
						std::sort(conns.begin(),conns.end(),RCM_Comparator(wbeg,xadj));
						for (k = 0; k < static_cast<INMOST_DATA_ENUM_TYPE>(conns.size()); ++k)
						{
							Ulist[conns[k]] = index++;
							q.push_back(conns[k]);
						}
						conns.clear();
					}
				}
				while( index < wend );

				for(k = wbeg; k < wend; ++k)
					Ulist[k] = wend-(Ulist[k]-wbeg)-1;

				for(k = wbeg; k < wend; ++k)
				{
					localP[k] = Ulist[k];
					localQ[k] = Ulist[k];
				}
				cend = wend;
				i = wend;


				trcmorder = Timer() - trcmorder;
			}
#elif defined(REORDER_WRCM)
			{
				if( verbosity > 1 ) 
					std::cout << "Reordering with WRCM\n";
				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
				//create a symmetric graph of the matrix A + A^T
				std::vector<INMOST_DATA_ENUM_TYPE> xadj(wend-wbeg+1);
				std::vector<INMOST_DATA_REAL_TYPE> wadj(wend-wbeg);
				//~ std::vector< std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> > adjncy;
				std::vector< INMOST_DATA_ENUM_TYPE > adjncy;
				//adjncy.reserve(nzA*2);
				
				trcmgraph = Timer();

				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) if( localQ[k] == ENUMUNDEF ) std::cout << __FILE__ << ":" << __LINE__ << " No column permutation for row " << k << ". Matrix is structurally singular\n";

				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) localP[k] = k;
				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) V[localQ[k]] = DR[k]; //if you expirience a bug here then the matrix is structurally singular
				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) DR[k] = V[k];
				
				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) V[localQ[k]] = DR0[k]; //if you expirience a bug here then the matrix is structurally singular
				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) DR0[k] = V[k];
				

#if defined(USE_OMP_FACT)
#pragma omp parallel for private(k)
#endif
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt < A_Address[k].last; ++jt)
						A_Entries[jt].first = localQ[A_Entries[jt].first];
					std::sort(A_Entries.begin()+A_Address[k].first,A_Entries.begin()+A_Address[k].last);
				}
				
				//DumpMatrix(A_Address,A_Entries,cbeg,cend,"A.mtx");
				
				ReorderEF(wbeg, wend, donePQ, localP, localQ);
				inversePQ(wbeg,wend,localP,localQ, invP,invQ);
				applyPQ(wbeg, wend, localP, localQ, invP, invQ);

				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
				
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					Ulist[k] = EOL;
					Bbeg[k] = EOL;
					Blist[k] = EOL;
				}
				for (k = wend; k > static_cast<INMOST_DATA_INTEGER_TYPE>(wbeg); --k)
				{
					if (A_Address[k-1].Size() > 0)
					{
						i = A_Entries[A_Address[k-1].first].first;
						Blist[k-1] = Bbeg[i];
						Bbeg[i] = k-1;
					}
					Ulist[k-1] = A_Address[k-1].first;
				}
				xadj[0] = 0;
				for(i = wbeg; i < wend; ++i)
				{
					wadj[i-wbeg] = 0;
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[i].first; jt < A_Address[i].last; ++jt)
					{
						if( A_Entries[jt].first != i  )//&& fabs(DL[i]*A_Entries[jt].second*DR[A_Entries[jt].first]) > 0.1*std::min(Acolmax,Arowmax) )
						{
							adjncy.push_back(static_cast<INMOST_DATA_ENUM_TYPE>(A_Entries[jt].first));
							wadj[i-wbeg] += fabs(DL[i]*A_Entries[jt].second*DR[A_Entries[jt].first]);
						}
					}

					Li = Bbeg[i];
					while (Li != EOL)
					{
						if( Li != i )//&& fabs(DL[A_Entries[A_Address[Li].first].first]*A_Entries[A_Address[Li].first].second*DR[i]) > 0.1*std::min(Acolmax,Arowmax) )
						{
							adjncy.push_back(static_cast<INMOST_DATA_ENUM_TYPE>(Li));
							wadj[i-wbeg] += fabs(DL[A_Entries[A_Address[Li].first].first]*A_Entries[A_Address[Li].first].second*DR[i]);
						}
						Li = Blist[Li];
					}

					Li = Bbeg[i];
					while (Li != EOL)
					{
						Ui = Blist[Li];
						Ulist[Li]++;
						if (A_Address[Li].Size() && A_Address[Li].last - Ulist[Li] > 0)
						{
							k = A_Entries[Ulist[Li]].first;
							Blist[Li] = Bbeg[k];
							Bbeg[k] = Li;
						}

						Li = Ui;
					}

					std::sort(adjncy.begin()+xadj[i-wbeg],adjncy.end());
					adjncy.resize(std::unique(adjncy.begin()+xadj[i-wbeg],adjncy.end())-adjncy.begin());

					xadj[i-wbeg+1] = static_cast<INMOST_DATA_ENUM_TYPE>(adjncy.size());
				}
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k) Bbeg[k] = EOL;

				trcmgraph = Timer()-trcmgraph;

				trcmorder = Timer();
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k) Ulist[k] = ENUMUNDEF;
				//find node with the lowest order
                INMOST_DATA_ENUM_TYPE index = wbeg;
				INMOST_DATA_ENUM_TYPE cur = ENUMUNDEF;
				std::deque<INMOST_DATA_ENUM_TYPE> q;
				//~ std::vector< std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> > conns;
				std::vector< INMOST_DATA_ENUM_TYPE > conns;
                do
				{
					cur = ENUMUNDEF;
					for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend) && cur == ENUMUNDEF; ++k)
					{
						if( Ulist[k] == ENUMUNDEF )
							cur = k;
					}
					if( cur == ENUMUNDEF ) break;
					assert(cur != ENUMUNDEF);
					for(k = cur+1; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k) if( Ulist[k] == ENUMUNDEF )
					{
						if( WRCM_Comparator(wbeg,wadj)(k,cur) )
							cur = k;
					}
					q.push_back(cur);
					Ulist[cur] = index++;
					bool finish_block = false;
					while(!q.empty())
					{
						cur = q.front();
						q.pop_front();
						if( finish_block )
						{
							for (INMOST_DATA_ENUM_TYPE it = xadj[cur-wbeg]; it < xadj[cur-wbeg+1]; ++it)
								if( Ulist[adjncy[it]] == ENUMUNDEF )
									Ulist[adjncy[it]] = ENUMUNDEF-1;
						}
						else
						{
							for (INMOST_DATA_ENUM_TYPE it = xadj[cur-wbeg]; it < xadj[cur-wbeg+1]; ++it)
								if( Ulist[adjncy[it]] == ENUMUNDEF )
									conns.push_back(adjncy[it]);
							std::sort(conns.begin(),conns.end(),WRCM_Comparator(wbeg,wadj));
							for (k = 0; k < static_cast<INMOST_DATA_INTEGER_TYPE>(conns.size()); ++k)
							{
								Ulist[conns[k]] = index++;
								q.push_back(conns[k]);
								//~ if( (index-wbeg)%512==0 ) finish_block = true;
							}
							conns.clear();
						}
					}
					
				}
				while( index < wend );
				
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
					if( Ulist[k] == ENUMUNDEF-1 ) Ulist[k] = index++;

				//reverse
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
					Ulist[k] = wend-(Ulist[k]-wbeg)-1;

				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					localP[k] = Ulist[k];//wend - Ulist[k] + 1;
					localQ[k] = Ulist[k];//wend - Ulist[k] + 1;
					//localP[Ulist[k]] = k;
					//localQ[Ulist[k]] = k;
				}
				cend = wend;
				i = wend;
				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

				trcmorder = Timer() - trcmorder;
			}
#elif defined(REORDER_BRCM)
			{
				if( verbosity > 1 ) 
					std::cout << "Reordering with BRCM\n";
				
				trcmgraph = Timer();

				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) if( localQ[k] == ENUMUNDEF ) std::cout << __FILE__ << ":" << __LINE__ << " No column permutation for row " << k << ". Matrix is structurally singular\n";

				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) localP[k] = k;
				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) V[localQ[k]] = DR[k]; //if you expirience a bug here then the matrix is structurally singular
				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) DR[k] = V[k];
				
				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) V[localQ[k]] = DR0[k]; //if you expirience a bug here then the matrix is structurally singular
				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) DR0[k] = V[k];

				for (k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt < A_Address[k].last; ++jt)
					{
						A_Entries[jt].first = localQ[A_Entries[jt].first];
						//A_Entries[jt].second *= DL[k]*DR[A_Entries[jt].first];
					}
				}
				
				//DumpMatrix(A_Address,A_Entries,cbeg,cend,"A.mtx");
				
				ReorderEF(wbeg, wend, donePQ, localP, localQ);
				inversePQ(wbeg,wend,localP,localQ, invP,invQ);
				applyPQ(wbeg, wend, localP, localQ, invP, invQ);

				
				trcmgraph = Timer()-trcmgraph;

				trcmorder = Timer();

				std::vector<Block> blocks;

				 NestedDissection(wbeg,wend,A_Address,A_Entries,localP,localQ,blocks,(wend-wbeg)/64);
				//KwayDissection(wbeg,wend,A_Address,A_Entries,localP,localQ,blocks,64);

				cend = wend;
				i = wend;
				

				trcmorder = Timer() - trcmorder;
			}
#elif defined(REORDER_BRCM)
			{
				if( verbosity > 1 ) 
					std::cout << "Reordering with BRCM\n";
				
				trcmgraph = Timer();

				for (k = cbeg; k < cend; ++k) if( localQ[k] == ENUMUNDEF ) std::cout << __FILE__ << ":" << __LINE__ << " No column permutation for row " << k << ". Matrix is structurally singular\n";

				for (k = cbeg; k < cend; ++k) localP[k] = k;
				for (k = cbeg; k < cend; ++k) V[localQ[k]] = DR[k]; //if you expirience a bug here then the matrix is structurally singular
				for (k = cbeg; k < cend; ++k) DR[k] = V[k];
				
				for (k = cbeg; k < cend; ++k) V[localQ[k]] = DR0[k]; //if you expirience a bug here then the matrix is structurally singular
				for (k = cbeg; k < cend; ++k) DR0[k] = V[k];

				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt < A_Address[k].last; ++jt)
					{
						A_Entries[jt].first = localQ[A_Entries[jt].first];
						//A_Entries[jt].second *= DL[k]*DR[A_Entries[jt].first];
					}
				}
				
				//DumpMatrix(A_Address,A_Entries,cbeg,cend,"A.mtx");
				
				ReorderEF(wbeg, wend, donePQ, localP, localQ);
				inversePQ(wbeg,wend,localP,localQ, invP,invQ);
				applyPQ(wbeg, wend, localP, localQ, invP, invQ);

				
				trcmgraph = Timer()-trcmgraph;

				trcmorder = Timer();
				std::vector<bool> mrk(wend-wbeg,false);
				std::vector<INMOST_DATA_REAL_TYPE> wadj_blk(wend-wbeg,0.0), wadj_sep(wend-wbeg,0.0);
				std::vector<INMOST_DATA_REAL_TYPE> keys0(wend-wbeg,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > tind(wbeg,wend), tind_sep(wbeg,wend);
				double prep_time = Timer();
				for(k = wbeg; k < wend; ++k) 
				{
					Llist[k] = ENUMUNDEF;
					Ulist[k] = ENUMUNDEF;
					Blist[k] = ENUMUNDEF;
					wadj_blk[k-wbeg] = A_Address[k].Size(); //row connection wgt
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it) //if( A_Entries[it].first != k )
						tind[A_Entries[it].first].push_back(k);
				}
				for(k = wbeg; k < wend; ++k) 
				{
					INMOST_DATA_ENUM_TYPE Beg = EOL;
					// go over connection of k-th row
					for(INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
					{
						INMOST_DATA_ENUM_TYPE kt = A_Entries[it].first;
						// assume diagonal in place
						if( Blist[kt] == ENUMUNDEF )
						{
							Blist[kt] = Beg;
							Beg = kt;
						}
						// each row that has same connection
						for(INMOST_DATA_ENUM_TYPE jt = 0; jt < tind[kt].size(); ++jt)
						{
							INMOST_DATA_ENUM_TYPE lt = tind[kt][jt];
							//insert to linked list
							if( Blist[lt] == ENUMUNDEF )
							{
								Blist[lt] = Beg;
								Beg = lt;
							}
						}
					}					
					INMOST_DATA_ENUM_TYPE it = Beg, jt;
					while( it != EOL )
					{
						if( it != k )
						{
							tind_sep[k].push_back(it);
							wadj_sep[it-wbeg]++;
						}
						jt = Blist[it]; //save next entry
						Blist[it] = ENUMUNDEF; //clear list
						it = jt; //restore next
					}
					
				}
				std::cout << "prepare time " << Timer() - prep_time << std::endl;
				
				INMOST_DATA_ENUM_TYPE index_col = wbeg, index_row = wbeg;
				INMOST_DATA_ENUM_TYPE cur = ENUMUNDEF;
				
				BinaryHeapCustom< INMOST_DATA_REAL_TYPE, std::less<INMOST_DATA_REAL_TYPE> > q0(wend-wbeg);
				//~ BinaryHeap q0(&keys0[0],wend-wbeg);
				//std::vector< std::pair<INMOST_DATA_REAL_TYPE,INMOST_DATA_ENUM_TYPE> > conns;
				INMOST_DATA_REAL_TYPE upd_sep = 1, upd_blk = 1;
				INMOST_DATA_REAL_TYPE wgt_sep = 1, wgt_blk = -1;
				std::cout << "wgt sep " << wgt_sep << " blk " << wgt_blk << std::endl;
				// cur = wbeg;
				for(k = wbeg; k < wend; ++k) 
				{
					// if( wadj_sep[k-wbeg] > wadj_sep[cur-wbeg] ) cur = k;
					q0.PushHeap(k-wbeg,wgt_sep*wadj_sep[k-wbeg]+wgt_blk*wadj_blk[k-wbeg]); 
				}

					
				//~ std::cout << "q0 size " << q0.Size() << std::endl;
				
				//Ulist // - row index
				//Llist // - col index
				INMOST_DATA_ENUM_TYPE Beg = EOL, sep_size = 0, row_size = 0, col_size = 0, blk_num = 0, tot_sep_size = 0;
				INMOST_DATA_ENUM_TYPE index_row_beg = index_row, index_col_beg = index_col;
				INMOST_DATA_REAL_TYPE min_ratio = 1.0e+20, ratio;
				INMOST_DATA_ENUM_TYPE min_index = ENUMUNDEF;
				std::vector< std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> > block_row, block_col;
				// goto start;
				while( !q0.Empty() ) 
				{
					//cur has minimal separator growth size
					cur = q0.PopHeap() + wbeg;
				// start:
					int had_col = index_col, had_sep = sep_size;
					if( Ulist[cur] != ENUMUNDEF ) continue; //skip separator of other blocks
					// std::cout << "(" << cur << "," << wadj_blk[cur-wbeg] << "," << wadj_sep[cur-wbeg] << ")";
					// std::cout << " row " << index_row << " col " << index_col << " sep " << sep_size << " update ";
					Ulist[cur] = index_row++;
					row_size++;
					if( Blist[cur] != ENUMUNDEF ) //uncount me from separator
					{
						sep_size--;
						// std::cout << "sep ";
					}
					else //include me to linked list to avoid count as separator
					{
						// std::cout << "free ";
						Blist[cur] = Beg;
						Beg = cur;
						//extract me from rows that depend on me
						for(INMOST_DATA_ENUM_TYPE it = 0; it < tind_sep[cur].size(); ++it)
						{
							INMOST_DATA_ENUM_TYPE kt = tind_sep[cur][it];
							//extract me from rows that depend on me
							wadj_sep[kt-wbeg]-=upd_sep;
							if( wgt_sep && q0.Contains(kt-wbeg) ) //TODO: keep separate linked-list to update keys
								q0.ChangeKey(kt-wbeg, wgt_sep*wadj_sep[kt-wbeg]+wgt_blk*wadj_blk[kt-wbeg]);
						}
					}
					// go over connection of enumerated row
					int col_skip = 0;
					for(INMOST_DATA_ENUM_TYPE it = A_Address[cur].first; it < A_Address[cur].last; ++it)
					{
						INMOST_DATA_ENUM_TYPE kt = A_Entries[it].first;
						if( Llist[kt] == ENUMUNDEF ) 
						{
							Llist[kt] = index_col++;
							col_size++;
							for(INMOST_DATA_ENUM_TYPE jt = 0; jt < tind[kt].size(); ++jt)
							{
								wadj_blk[tind[kt][jt]-wbeg]-=upd_blk;
								if( q0.Contains(kt-wbeg) && wgt_blk )
									q0.ChangeKey(kt-wbeg,wgt_sep*wadj_sep[kt-wbeg]+wgt_blk*wadj_blk[kt-wbeg]);
							}
						}
						else col_skip++;
					}
					// add dependencies to linked-list of separator
					int sep_skip = 0;
					for(INMOST_DATA_ENUM_TYPE it = 0; it < tind_sep[cur].size(); ++it)
					{
						INMOST_DATA_ENUM_TYPE kt = tind_sep[cur][it];
						//insert to linked-list for separator
						if( Blist[kt] == ENUMUNDEF ) //new row into separator
						{
							Blist[kt] = Beg;
							Beg = kt;
							//extract me from rows that depend on me
							for(INMOST_DATA_ENUM_TYPE jt = 0; jt < tind_sep[kt].size(); ++jt)
							{
								INMOST_DATA_ENUM_TYPE qt = tind_sep[kt][jt];
								wadj_sep[qt-wbeg]-=upd_sep;
								//~ if( wadj_sep[qt-wbeg] < 0 && Ulist[qt] == ENUMUNDEF ) std::cout << std::endl << " ! terminal separator " << qt << " wgt " << wadj_sep[qt-wbeg] <<  std::endl;
								if( q0.Contains(qt-wbeg) && wgt_sep ) //TODO: keep separate linked-list to update keys
									q0.ChangeKey(qt-wbeg, wgt_sep*wadj_sep[qt-wbeg]+wgt_blk*wadj_blk[qt-wbeg]);
							}
							wadj_sep[kt-wbeg]-=upd_sep; //taking me will reduce separator size
							//~ if( wadj_sep[kt-wbeg] < 0 && Ulist[kt] == ENUMUNDEF ) std::cout << std::endl << " ! terminal separator " << kt << " wgt " << wadj_sep[kt-wbeg] << std::endl;
							if( q0.Contains(kt-wbeg) && wgt_sep ) //TODO: keep separate linked-list to update keys
								q0.ChangeKey(kt-wbeg, wgt_sep*wadj_sep[kt-wbeg]+wgt_blk*wadj_blk[kt-wbeg]);
							sep_size++;
						}
						else 
						{
							//~ std::cout << "skip " << kt << std::endl;
							sep_skip++;
						}
					}
					//disconnect each row that connects to me to know block growth size
					//~ for(INMOST_DATA_ENUM_TYPE jt = 0; jt < tind[cur].size(); ++jt)
						//~ wadj_blk[tind[cur][jt]-wbeg]--;
					//~ std::cout << " rows+sep " << row_size + sep_size << " cols " << col_size << " cols-rows " << col_size-row_size << " sep " << sep_size << std::endl;
					// std::cout << " rows " << row_size << " cols " << col_size << " sep " << sep_size << std::endl;
					INMOST_DATA_ENUM_TYPE nxt = q0.PeekHeap(), N = wend-wbeg, b1 = row_size, s = sep_size, b2 = N-b1-s;
					//~ if( col_size > row_size + sep_size ) //form new block
					//~ if( (row_size > part_size && !(nxt != ENUMUNDEF && wadj_sep[nxt-wbeg] <= 0)) || (col_size > part_size && !(nxt != ENUMUNDEF && wadj_blk[nxt-wbeg] <= 0)) )
					// if( (row_size > 40*sep_size) && !(upd_sep == 1 && nxt != ENUMUNDEF && wadj_sep[nxt-wbeg] <= 0) )
					if( col_size > (wend-wbeg)/8 && !(nxt != ENUMUNDEF && wadj_sep[nxt-wbeg] <= 0) )
					// ratio = s / sqrt(b1*b2);
					// if( ratio < min_ratio )
					// {
					// 	std::cout << "new minimum ratio b1 " << b1 << " b2 " << b2 << " s " << s << " ratio " << ratio << std::endl;
					// 	min_ratio = ratio;
					// 	min_index = index_row;
					// }
					// if( false ) if( (row_size > (wend-wbeg-tot_sep_size-sep_size)/2) && !(upd_sep == 1 && nxt != ENUMUNDEF && wadj_sep[nxt-wbeg] <= 0) )
					{
						
						INMOST_DATA_ENUM_TYPE it = Beg, jt;
						while( it != EOL )
						{
							if( Ulist[it] == ENUMUNDEF ) 
							{
								Ulist[it] = ENUMUNDEF-1; //enumerate separator
								tot_sep_size++;
							}
							it = Blist[it];
							//clear linked list
							//~ jt = Blist[it];
							//~ Blist[it] = ENUMUNDEF;
							//~ it = jt;
						}
						//~ Beg = EOL;
						std::cout << "make new block! total separator " << tot_sep_size << " gained " << sep_size << " rows " << index_row-index_row_beg << " cols " << index_col-index_col_beg << std::endl;
						blk_num++;
						block_row.push_back( std::make_pair(index_row_beg,index_row) );
						block_col.push_back( std::make_pair(index_col_beg,index_col) );
						index_row_beg = index_row;
						index_col_beg = index_col;

						col_size = row_size = sep_size = 0;
					}
				}
				//add last block
				if( index_row_beg != index_row || index_col_beg != index_col )
				{
					std::cout << "make new block! total separator " << tot_sep_size << " rows " << index_row-index_row_beg << " cols " << index_col-index_col_beg << std::endl;
					blk_num++;
					block_row.push_back( std::make_pair(index_row_beg,index_row) );
					block_col.push_back( std::make_pair(index_col_beg,index_col) );
					index_row_beg = index_row;
					index_col_beg = index_col;
				}
				
				std::cout << "blocks: " << blk_num << " total separator " << tot_sep_size << "/" << wend-wbeg << std::endl;
				
				//~ std::set<INMOST_DATA_ENUM_TYPE> block_cols;
				/*
				std::vector<INMOST_DATA_REAL_TYPE> keys(wend-wbeg,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				BinaryHeap q(&keys[0], wend-wbeg);
				std::vector< std::pair< INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > > blocks;
				INMOST_DATA_ENUM_TYPE block_start = index;
                do
				{
					cur = ENUMUNDEF;
					while( cur == ENUMUNDEF && !q0.Empty() ) 
					{
						cur = q0.PopHeap()+wbeg;
						if( Ulist[cur] != ENUMUNDEF ) cur = ENUMUNDEF;
					}
					
					if( cur == ENUMUNDEF ) break;
					
					q.PushHeap(cur-wbeg,wadj[cur-wbeg]);
					//~ std::cout << "pick (" << cur << "," << wadj[cur-wbeg] << ") heap size " << q.Size() << std::endl;
					for(INMOST_DATA_ENUM_TYPE jt = 0; jt < tind[cur].size(); ++jt)
					{
						wadj[tind[cur][jt]-wbeg]--;
						if( q0.Contains(tind[cur][jt]-wbeg) ) q0.ChangeKey(tind[cur][jt]-wbeg,wadj[tind[cur][jt]-wbeg]);
						if( q.Contains(tind[cur][jt]-wbeg) ) q.DecreaseKey(tind[cur][jt]-wbeg,wadj[tind[cur][jt]-wbeg]);
					}
					
					while(!q.Empty())
					{
						cur = q.PopHeap()+wbeg;
						INMOST_DATA_ENUM_TYPE qs = q.Size();
						INMOST_DATA_REAL_TYPE wc = wadj[cur-wbeg];
						Ulist[cur] = index++;
						for (INMOST_DATA_ENUM_TYPE it = 0; it < tind[cur].size(); ++it)
						{
							INMOST_DATA_ENUM_TYPE kt = tind[cur][it];
							if( Ulist[kt] == ENUMUNDEF )
							{
								if( !q.Contains(kt-wbeg) )
								{ 
									q.PushHeap(kt-wbeg,wadj[kt-wbeg]);
									for(INMOST_DATA_ENUM_TYPE jt = 0; jt < tind[kt].size(); ++jt)
									{
										wadj[tind[kt][jt]-wbeg]--;
										if( q0.Contains(tind[kt][jt]-wbeg) ) q0.ChangeKey(tind[kt][jt]-wbeg,wadj[tind[kt][jt]-wbeg]);
										if( q.Contains(tind[kt][jt]-wbeg) ) q.DecreaseKey(tind[kt][jt]-wbeg,wadj[tind[kt][jt]-wbeg]);
									}
								}
							}
						}
						//~ std::cout << "add (" << cur << "," << wc << ") block_size " << index-block_start << " heap size " << q.Size() << " inc " << q.Size() - qs << std::endl;
						//~ if( index-block_start + q.Size() > (wend-wbeg)/2 )
						if( false )//index-block_start > 10*q.Size() )
						{
							//~ std::cout << "finalize block, size " << index-block_start << " heap size " << q.Size() << std::endl;
							while( !q.Empty() ) 
							{ 
								cur = q.PopHeap()+wbeg;
								if( wadj[cur-wbeg] == 0 )
									Ulist[cur] = index++;
								else
									Ulist[cur] = ENUMUNDEF-1; 
							}
						}
					}
					if( !blocks.empty() && (index-block_start < (blocks.back().second-blocks.back().first)/5) ) //attach to last block
						blocks.back().second = index;
					else if( block_start != index ) //exhausted block
						blocks.push_back(std::make_pair(block_start,index));
					block_start = index;
				}
				while( index < wend );
				
				
				
				for(k = wbeg; k < wend; ++k)
					if( Ulist[k] == ENUMUNDEF-1 ) Ulist[k] = index++;

				blocks.push_back(std::make_pair(block_start,index));
				
				
				std::cout << "total_blocks " << blocks.size()-1 << " separator size " << index-block_start << "/" << wend-wbeg << std::endl;
				for(k = 0; k < blocks.size(); ++k)
					std::cout << "block[" << k << "]: " <<  blocks[k].first << ":" << blocks[k].second << " (" << blocks[k].second-blocks[k].first << ")" << std::endl;
				
				*/
				
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					if( Ulist[k] == ENUMUNDEF-1 ) Ulist[k] = index_row++;
					//~ if( Llist[k] == ENUMUNDEF-1 ) Llist[k] = index_col++;
					//~ if( Ulist[k] == ENUMUNDEF ) std::cout << "skipped row " << k << std::endl;
					if( Llist[k] == ENUMUNDEF ) Llist[k] = index_col++; //std::cout << "skipped col " << k << std::endl;
				}
				//add separator block
				std::cout << "make separator block! total separator " << tot_sep_size << " rows " << index_row-index_row_beg << " extra cols " << index_col-index_col_beg << std::endl;
				blk_num++;
				block_row.push_back( std::make_pair(index_row_beg,index_row) );
				block_col.push_back( std::make_pair(wbeg,index_col) );
				std::ofstream file("blocks.txt");
				for(k = 0; k < blk_num; ++k)
				{
					std::cout << "block[" << k << "] rows " << block_row[k].first << ":" << block_row[k].second << "(" << block_row[k].second - block_row[k].first << ") cols " << block_col[k].first << ":" << block_col[k].second << "(" << block_col[k].second - block_col[k].first << ")" << std::endl;
					file << block_row[k].first << " " << block_row[k].second << " " << block_col[k].first << " " <<  block_col[k].second << std::endl;
				}
				file.close();
				//reverse
				//~ for(k = wbeg; k < wend; ++k)
					//~ Ulist[k] = wend-(Ulist[k]-wbeg)-1;

				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					localP[k] = Ulist[k];//wend - Ulist[k] + 1;
					//~ localP[k] = k;
					localQ[k] = Llist[k];//wend - Ulist[k] + 1;
					//~ localQ[k] = k;
					//localP[Ulist[k]] = k;
					//localQ[Ulist[k]] = k;
				}
				cend = wend;
				i = wend;


				trcmorder = Timer() - trcmorder;
			}
#else
			cend = wend;
			i = wbeg;
#endif
			tt = Timer() - tt;
			treorder += tt;
			if( verbosity > 1 )
			{
				std::cout << "Time " << tt << "\n";
				std::cout << "Reorder\n";
			}
			tt = Timer();
			if (cbeg == cend && cbeg != wend)
			{
				//std::cout << __FILE__ << ":" << __LINE__ << " singular matrix, factored " << mobeg << ".." << cend << " out of " << mobeg << ".." << moend << std::endl;
				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
					LU_Diag[k] = tol_modif;
				break;
			}
			
			//finish reordering
			if (i < wend)
			{
				j = i;
				for (k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					if (localP[k] == ENUMUNDEF) localP[k] = i++;
					if (localQ[k] == ENUMUNDEF) localQ[k] = j++;
				}
			}

			for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) 
			{
				U[localP[k]] = DL[k];
				V[localQ[k]] = DR[k];
			}
			for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) 
			{
				DL[k] = U[k];
				DR[k] = V[k];
			}
			
			for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) 
			{
				if( !(localQ[k] >= cbeg && localQ[k] < cend) ||
				    !(localP[k] >= cbeg && localP[k] < cend))
					std::cout << "Bad permutations P: " << localP[k] << " Q: " << localQ[k] << " interval [" << cbeg << ":" << cend << "]" << std::endl;
				U[localP[k]] = DL0[k];
				V[localQ[k]] = DR0[k];
			}
			for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) 
			{
				DL0[k] = U[k];
				DR0[k] = V[k];
			}

			
			tlreorder = Timer() - tt;
			treorder += tlreorder;
			if( verbosity > 1 )
				std::cout << "Time " << tlreorder << "\n";

			if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
			/// REASSAMBLE
			if( verbosity > 1 ) 
				std::cout << "Reassemble\n";
			//tt = Timer();
			double tt1, tt2,ttt;
			//in localPQ numbers indicate where to put current row/column
			//reorder E,F blocks by swaps
			tt1 = Timer();
			//inverse ordering
			ReorderEF(wbeg,wend, donePQ, localP, localQ);
			inversePQ(wbeg,wend,localP,localQ, invP,invQ);
			tt1 = Timer() - tt1;

			//std::cout << "reorder: " << tt1 << std::endl;
			tt2 = Timer();
			nzA = 0;
			{
				//std::vector<Sparse::Row::entry> B_Entries;
				interval<INMOST_DATA_ENUM_TYPE, Interval> B_Address(A_Address.get_interval_beg(),A_Address.get_interval_end());
				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
				{
					B_Address[k].first = A_Address[invP[k]].first;// static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[invP[k]].first; jt < A_Address[invP[k]].last; ++jt)
					{
						i = localQ[A_Entries[jt].first];
						//u = A_Entries[jt].second;
						A_Entries[jt].first = i;
						//B_Entries.push_back(Sparse::Row::make_entry(i, u));
					}
					B_Address[k].last = A_Address[invP[k]].last;//static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());
					nzA += A_Address[k].Size();
				}
				//B_Entries.push_back(Sparse::Row::make_entry(-1, 0.0));
				A_Address.swap(B_Address);
				//A_Entries.swap(B_Entries);
			}
			if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

#if defined(USE_OMP_FACT)
#pragma omp parallel for private(k)
#endif
			for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
				std::sort(A_Entries.begin() + A_Address[k].first, A_Entries.begin() + A_Address[k].last);
			
			if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

			double sparsity = nzA/(double)(cend-cbeg)/(double)(cend-cbeg);
			
			if( verbosity > 1 )
				std::cout << "nzA " << nzA << " size " << cend-cbeg << " sparsity " << sparsity << "\n";
			
			if ( (sparsity > 0.85) || (cend-cbeg < 32) )
			{
				block_pivot = true;
			}
			
			
			tt2 = Timer() - tt2;

			ttt = Timer();
			applyPQ(wbeg, wend, localP, localQ, invP, invQ);
			tt1 += Timer() - ttt;
			
			//reset localPQ
			
			tlreassamble = Timer() - tt;
			treassamble += tlreassamble;
			if( verbosity > 1 )
				std::cout << "Time " << tlreassamble << "\n";
			/// RESCALING

			if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
	
			//Rescale current B block by Row-Column alternating scaling
			tt = Timer();
			if( rescale_b )
			{
				if( verbosity > 1 )
					std::cout << " rescaling block B, iters " << sciters << std::endl;


#if defined(EQUALIZE_IDOMINANCE)
				std::vector<INMOST_DATA_REAL_TYPE> C_Entries(A_Entries.size());
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> temp(cbeg,cend,0);
#endif
#if defined(USE_OMP_FACT)
#pragma omp parallel
#endif
				{

				
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
					{
						for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
							A_Entries[r].second *= (DL[k] * DR[A_Entries[r].first]);
						
						U[k] = DL[k];
						V[k] = DR[k];
					}
#if defined(EQUALIZE_1NORM)
///////////////////////////////////////////////////////////////////////////////////
///         ROW-COLUMN ALTERNATING SCALING FOR 1-NORM BALANCING                 ///
///////////////////////////////////////////////////////////////////////////////////
					//std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, 0.0);
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) DL[k] = 0.0;
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
					{
						for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
							DL[k] += fabs(A_Entries[r].second);
					}
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
					for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
					{
						//std::fill(DR.Begin()+cbeg-mobeg, DR.Begin()+cend-mobeg, 0.0);
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
						for(k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) DR[k] = 0.0;
					
#if defined(USE_OMP_FACT)
#pragma omp for private(k,u)
#endif
						for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
						{
							for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
							{
								u = DL[k] * fabs(A_Entries[r].second);
#if defined(USE_OMP_FACT)
#pragma omp atomic
#endif
								DR[A_Entries[r].first] += u;
							}
						}
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
						for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
						//std::fill(DL.Begin()+cbeg-mobeg, DL.Begin()+cend-mobeg, 0.0);
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
						for(k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) DL[k] = 0.0;
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
						for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
						{
							for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
								DL[k] += DR[A_Entries[r].first] * fabs(A_Entries[r].second);
						}
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
						for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
					}
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
					{
						for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
							A_Entries[r].second *= DL[k] * DR[A_Entries[r].first];
					}
#elif defined(EQUALIZE_2NORM)
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for(k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) DL[k] = 0.0;
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
					{
						for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
							DL[k] += A_Entries[r].second*A_Entries[r].second;
					}
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
					for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
					{
						//std::fill(DR.Begin()+cbeg-mobeg, DR.Begin()+cend-mobeg, 0.0);
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
						for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) DR[k] = 0.0;
#if defined(USE_OMP_FACT)
#pragma omp for private(k,u)
#endif
						for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
						{
							for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
							{
								u = DL[k] * A_Entries[r].second*A_Entries[r].second;
#if defined(USE_OMP_FACT)
#pragma omp atomic
#endif
								DR[A_Entries[r].first] += u;
							}
						}
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
						for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
						for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) DL[k] = 0.0;
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
						for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
						{
							for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
								DL[k] += DR[A_Entries[r].first] * A_Entries[r].second*A_Entries[r].second;
						}
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
						for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
					}
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
					{
						for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
							A_Entries[r].second *= DL[k] * DR[A_Entries[r].first];
					}
#elif defined(EQUALIZE_IDOMINANCE)
					/// THIS VERSION OF RESCALING INCREASES DIAGONAL DOMINANCE
#if defined(USE_OMP_FACT)
#pragma omp for private(k,u)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
					{
						for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
						{
							u = fabs(A_Entries[r].second);
							if( u > 1.0e-12 )
								C_Entries[r] = -log(u);
							else
								C_Entries[r] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
						}
					}

				
					

					for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
					{
						//std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
						//std::fill(DR.Begin() + cbeg - mobeg, DR.Begin() + cend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
						for(k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) DL[k] = DR[k] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
					
#if defined(USE_OMP_FACT)
//#pragma omp for private(k,i,u)
#pragma omp single
#endif
						for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++) //row number
						{
							for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
							{
								i = A_Entries[r].first; //column number
								if( static_cast<INMOST_DATA_INTEGER_TYPE>(i) != k ) //out of diagonal
								{
									u = C_Entries[r] + temp[k] - temp[i];
									//if( isnan(u) || u != u ) std::cout << __FILE__ << ":" << __LINE__ << " u is " << u << std::endl;
									DL[k] = std::min(DL[k],u);// update Y1
									DR[i] = std::min(DR[i],u);// update Y2
								}
							}
						}

#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
						for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
						{
							if( DR[k] != std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() &&
								DL[k] != std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() )
								temp[k] += (DR[k]-DL[k])*0.5;
						}

					}

#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
					{
						DL[k] = exp(-temp[k]);
						DR[k] = exp(temp[k]);
						if( DL[k] != DL[k] ) DL[k] = 1;
						if( DR[k] != DR[k] ) DR[k] = 1;
						//if( isnan(DL[k]) || DL[k] != DL[k] || fabs(DL[k]) < 1.0e-12 ) std::cout << __FILE__ << ":" << __LINE__ << " DL[" << k << "] is " << DL[k] << std::endl;
						//if( isnan(DR[k]) || DR[k] != DR[k] || fabs(DR[k]) < 1.0e-12 ) std::cout << __FILE__ << ":" << __LINE__ << " DR[" << k << "] is " << DR[k] << std::endl;
					}
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
					{
						for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
							A_Entries[r].second *= (DL[k] * DR[A_Entries[r].first]);
					}
#endif
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
					{
						DL[k] *= U[k];
						DR[k] *= V[k];
					}
				//DumpMatrix(B_Address,B_Entries,cbeg,cend,"mat_equilibration.mtx");
				/// RESCALING DONE
				//stack scaling
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
					for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
					{
						DL0[k] *= DL[k];
						DR0[k] *= DR[k];
					}
				
					//rescale EF
					if( !level_size.empty() )
					{
						INMOST_DATA_ENUM_TYPE first = mobeg, last;
						INMOST_DATA_ENUM_TYPE kbeg, kend;
						for(size_t level = 0; level < level_size.size(); ++level)
						{
							last = first + level_size[level];
							kbeg = std::max(cbeg,last);
							kend = std::min(cend,moend);
							//std::cout << "Rescale EF level " << level << " ";
							//std::cout << "interval [" << first << "," << last << "] ";
							//std::cout << "matrix [" << mobeg << "," << moend << "] ";
							//std::cout << "rescaling [" << cbeg << "," << cend << "] ";
							//std::cout << "selected [" << kbeg << ":" << kend << "]" << std::endl;
#if defined(USE_OMP_FACT)
#pragma omp for private(k)
#endif
							for (k = kbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(kend); ++k) 
							{
								if( F_Address[level]->at(k).Size() )
								{
									for (INMOST_DATA_ENUM_TYPE r = F_Address[level]->at(k).first; r < F_Address[level]->at(k).last; ++r)
										F_Entries[r].second *= DR[k];
								}
								if( E_Address[level]->at(k).Size() )
								{
									for (INMOST_DATA_ENUM_TYPE r = E_Address[level]->at(k).first; r < E_Address[level]->at(k).last; ++r)
										E_Entries[r].second *= DL[k];
								}
							}
							first = last;
						}
					}
				}
				//End rescale B block
				tlrescale = Timer() - tt;
				trescale += tlrescale;
				if( verbosity > 1 )
					std::cout << "Time " << tlrescale << "\n";

			}

			if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
			/*
			{
				std::stringstream s;
				s << "B" << level_size.size() << ".mtx";
				DumpMatrix(A_Address,A_Entries,cbeg,cend,s.str());
				std::cout << "Press any key" << std::endl;
				scanf("%*c");
			}
			*/
			/// FACTORIZATION BEGIN
			{
				//Setup column addressing for B,F, in descending order to keep list ordered

				//~ Clear TransposeB
				//~ PrepareTranspose(cbeg,cend,A_Address,A_Entries,TransposeA);
				for (k = cend; k > static_cast<INMOST_DATA_INTEGER_TYPE>(cbeg); --k)
				{
					if (A_Address[k - 1].Size() > 0)
					{
						i = A_Entries[A_Address[k - 1].first].first;
						if (static_cast<INMOST_DATA_INTEGER_TYPE>(i) < k - 1)
						{
							Blist[k - 1] = Bbeg[i];
							Bbeg[i] = k - 1;
						}
					}
				}

				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

				//supplimentary data structures for condition estimates of L^{-1}, U^{-1}
				INMOST_DATA_REAL_TYPE mu_update_L1, mu_update_L2, mu_update_U1, mu_update_U2;
				INMOST_DATA_REAL_TYPE NuU_old = 1, NuL_old = 1, NuD_old = 1;

				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LineValuesU(cbeg, cend, 0.0), LineValuesL(cbeg, cend, 0.0);
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> LineIndecesU(cbeg, cend + 1, UNDEF), LineIndecesL(cbeg, cend + 1, UNDEF);
				std::vector<INMOST_DATA_ENUM_TYPE> indicesU, indicesL;

				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL1(cbeg, cend, 0.0), EstU1(cbeg, cend, 0.0);
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL2(cbeg, cend, 0.0), EstU2(cbeg, cend, 0.0);
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Bstart(cbeg, cend);

				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ubeg(cbeg, cend, EOL), Lbeg(cbeg, cend, EOL);

				L_Beg = static_cast<INMOST_DATA_ENUM_TYPE>(L_Entries.size());
				U_Beg = static_cast<INMOST_DATA_ENUM_TYPE>(U_Entries.size());
#if defined(ILUC2)
				nzLU2tot += nzLU2;
				nzLU2 = 0;
				L2_Entries.clear();
				U2_Entries.clear();
#endif

				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

				//setup diagonal values and memorize indices
				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
				{
					LU_Diag[k] = 0.0;
					Bstart[k] = A_Address[k].first; //remember B block starts permutation for pivoting
					for (i = A_Address[k].first; i < A_Address[k].last; i++)
					{
						if (static_cast<INMOST_DATA_INTEGER_TYPE>(A_Entries[i].first) == k)
						{
							LU_Diag[k] = A_Entries[i].second;
							break;
						}
					}
				}
				//Initialize data for condition estimator
				tt = Timer();
				if( estimator )
				{
					NuU_acc *= NuU;
					NuL_acc *= NuL;
					NuD_acc *= NuD;
					NuU = NuL = 1.0;
					NuU_max = NuL_max = 1.0;
				}

				max_diag = 0;
				min_diag = 1.0e300;
				NuD = 1.0;
				if( verbosity > 1 )
					std::cout << " starting factorization " << std::endl;
					

				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

				for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
				{
#if defined(USE_OMP_FACT)
#pragma omp parallel sections num_threads(2)
#endif
					{
						// Compute U-part
#if defined(USE_OMP_FACT)
#pragma omp section
#endif
						{
							// Prepare linked list for U-part
							INMOST_DATA_ENUM_TYPE Sbeg, i;
							INMOST_DATA_REAL_TYPE coef;
							//uncompress k-th row
							// add diagonal value first, there shouldn't be values on left from diagonal
							//assert(B_Entries[A_Address[k].first].first == k);
							LineIndecesU[cbeg] = k;
							LineIndecesU[k] = EOL;
							LineValuesU[k] = 0.0;
							if (A_Address[k].Size())
							{
								if (static_cast<INMOST_DATA_INTEGER_TYPE>(A_Entries[A_Address[k].first].first) == k)
									LineValuesU[k] = A_Entries[A_Address[k].first].second;
								else
									LineValuesU[k] = 0.0;
								i = k;
								for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first + (static_cast<INMOST_DATA_INTEGER_TYPE>(A_Entries[A_Address[k].first].first) == k ? 1 : 0); it < A_Address[k].last; ++it)
								{
									LineValuesU[A_Entries[it].first] = A_Entries[it].second;
									i = LineIndecesU[i] = A_Entries[it].first;
								}
								LineIndecesU[i] = EOL;
							}
							Sbeg = k;
							// U-part elimination with L
							i = Lbeg[k];
							while (i != EOL)
							{
								assert(i != UNDEF);
								assert(static_cast<INMOST_DATA_INTEGER_TYPE>(L_Entries[L_Address[i].first].first) == k);
								coef = -L_Entries[L_Address[i].first].second * LU_Diag[i];
								AddListUnordered(Sbeg, U_Address[i], U_Entries, coef, LineIndecesU, LineValuesU, 0.0);
#if defined(ILUC2)
								AddListUnordered(Sbeg, U2_Address[i], U2_Entries, coef, LineIndecesU, LineValuesU, 0.0);
#endif
								i = Llist[i];
							}
#if defined(ILUC2)
							// U-part elimination with second-order L 
							i = L2beg[k];
							while (i != EOL)
							{
								assert(i != UNDEF);
								assert(static_cast<INMOST_DATA_INTEGER_TYPE>(L2_Entries[L2_Address[i].first].first) == k);
								coef = -L2_Entries[L2_Address[i].first].second * LU_Diag[i];
								AddListUnordered(Sbeg, U_Address[i], U_Entries, coef, LineIndecesU, LineValuesU, 0.0);
#if 0
								AddListUnordered(Sbeg, U2_Address[i], U2_Entries, coef, LineIndecesU, LineValuesU, 0.0);
#endif
								i = L2list[i];
							}
#endif
							// Order contents of linked list
							OrderList(Sbeg, LineIndecesU, indicesU);
							// Rescale with diagonal
#if defined(PIVOT_THRESHOLD)
							if (fabs(LineValuesU[k]) < tau2) LineValuesU[k] = (LineValuesU[k] < 0.0 ? -1 : 1) * tau2;
#endif
							ScaleList(1.0 / LineValuesU[k], LineIndecesU[k], LineIndecesU, LineValuesU);
							// Condition estimator for U part
							if (estimator)
							{
								tlocal = Timer();
								NuU_old = NuU;
								NuU = Estimator1(k, LineIndecesU, LineValuesU, EstU1, mu_update_U1);
								NuU = std::max(NuU, Estimator2(k, LineIndecesU, LineValuesU, EstU2, mu_update_U2));
								testimator += Timer() - tlocal;
							} //if( estimator )
						}
						// Compute L-part
#if defined(USE_OMP_FACT)
#pragma omp section
#endif
						{
							// prepearing linked list for L-part
							INMOST_DATA_ENUM_TYPE Sbeg, i, j;
							INMOST_DATA_REAL_TYPE coef;
							//uncompress column
							//insert diagonal value first
							LineIndecesL[cbeg] = k;
							LineIndecesL[k] = EOL;
							LineValuesL[k] = 0.0;
							if (A_Address[k].Size())
							{
								if (static_cast<INMOST_DATA_INTEGER_TYPE>(A_Entries[A_Address[k].first].first) == k)
									LineValuesL[k] = A_Entries[A_Address[k].first].second;
								else
									LineValuesL[k] = 0.0;
								//start from diagonal
								j = k;
								i = Bbeg[k];
								while (i != EOL)
								{
									assert(static_cast<INMOST_DATA_INTEGER_TYPE>(B_Entries[B_Address[i].first].first) == k);
									LineValuesL[i] = A_Entries[A_Address[i].first].second;
									j = LineIndecesL[j] = i;
									i = Blist[i];
								}
								LineIndecesL[j] = EOL;
							}
							Sbeg = cbeg;
							// L-part elimination with U
							i = Ubeg[k];
							while (i != EOL)
							{
								assert(i != UNDEF);
								assert(static_cast<INMOST_DATA_INTEGER_TYPE>(U_Entries[U_Address[i].first].first) == k);
								coef = -U_Entries[U_Address[i].first].second * LU_Diag[i];
								AddListUnordered(Sbeg, L_Address[i], L_Entries, coef, LineIndecesL, LineValuesL, 0.0);
#if defined(ILUC2)
								AddListUnordered(Sbeg, L2_Address[i], L2_Entries, coef, LineIndecesL, LineValuesL, 0.0);
#endif
								i = Ulist[i];
							}
#if defined(ILUC2)
							// L-part elimination with second-order U
							i = U2beg[k];
							while (i != EOL)
							{
								assert(i != UNDEF);
								assert(static_cast<INMOST_DATA_INTEGER_TYPE>(U2_Entries[U2_Address[i].first].first) == k);
								coef = -U2_Entries[U2_Address[i].first].second * LU_Diag[i];
								AddListUnordered(Sbeg, L_Address[i], L_Entries, coef, LineIndecesL, LineValuesL, 0.0);
#if 0
								AddListUnordered(Sbeg, L2_Address[i], L2_Entries, coef, LineIndecesL, LineValuesL, 0.0);
#endif
								i = U2list[i];
							}
#endif
							// Order contents of linked list
							OrderList(Sbeg, LineIndecesL, indicesL);
							// Rescale with diagonal
#if defined(PIVOT_THRESHOLD)
							if (fabs(LineValuesL[k]) < tau2) LineValuesL[k] = (LineValuesL[k] < 0.0 ? -1 : 1) * tau2;
#endif
							ScaleList(1.0 / LineValuesL[k], LineIndecesL[k], LineIndecesL, LineValuesL);
							// Estimate condition for L^{-1}
							if (estimator)
							{
								tlocal = Timer();
								NuL_old = NuL;
								NuL = Estimator1(k, LineIndecesL, LineValuesL, EstL1, mu_update_L1);
								NuL = std::max(NuL, Estimator2(k, LineIndecesL, LineValuesL, EstL2, mu_update_L2));
								testimator += Timer() - tlocal;
							}
						}
					}
					LU_Diag[k] = (LineValuesU[k]+LineValuesL[k])*0.5;
					// Condition estimator for diagonal D
					NuD_old = NuD;
					NuD = std::max<INMOST_DATA_REAL_TYPE>(fabs(LU_Diag[k]),max_diag) / std::min<INMOST_DATA_REAL_TYPE>(fabs(LU_Diag[k]),min_diag);
					//discard line if condition number is too large
#if defined(CONDITION_PIVOT)
					if(allow_pivot && !block_pivot && (NuU > pivot_cond || NuL > pivot_cond || NuD > pivot_diag) )
					{
						//restore condition number
						EstU1[k] = EstL1[k] = 0;
						EstU2[k] = EstL2[k] = 0;
						//Cancel this row
						Pivot[k] = true;
						swaps++;
					}
#endif
#if defined(NNZ_PIVOT)
					if( !Pivot[k] )
					{
						nnzL = 0, nnzU = 0;
						Li = LineIndecesL[k];
						while (Li != EOL)
						{
							l = fabs(LineValuesL[Li]);
							if (l*NuL > tau ) nnzL++;
							Li = LineIndecesL[Li];
						}
						Ui = LineIndecesU[k];
						while (Ui != EOL)
						{
							u = fabs(LineValuesU[Ui]);
							if (u*NuU > tau ) nnzU++;
							Ui = LineIndecesU[Ui];
						}
						if( allow_pivot && !block_pivot && (k != cbeg && nnzL > NNZ_GROWTH_PARAM*nnzL_max && nnzU > NNZ_GROWTH_PARAM*nnzU_max) )
						{
							//restore condition number
							EstU1[k] = EstL1[k] = 0;
							EstU2[k] = EstL2[k] = 0;
							//Cancel this row
							Pivot[k] = true;
							swaps++;
						}
						if( !Pivot[k] )
						{
							nnzL_max = std::max(nnzL,nnzL_max);
							nnzU_max = std::max(nnzU,nnzU_max);
						}
					}
#endif //NNZ_PIVOT
					if( Pivot[k] )
					{
						NuL = NuL_old;
						NuU = NuU_old;
						NuD = NuD_old;
						U_Address[k].first = U_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(U_Entries.size());
						L_Address[k].first = L_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(L_Entries.size());
#if defined(ILU2)
						U2_Address[k].first = U2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(U2_Entries.size());
						L2_Address[k].first = L2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(L2_Entries.size());
#endif
					}
					else
					{
						NuL_max = std::max(NuL,NuL_max);
						NuU_max = std::max(NuU,NuU_max);
						max_diag = std::max<INMOST_DATA_REAL_TYPE>(fabs(LU_Diag[k]),max_diag);
						min_diag = std::min<INMOST_DATA_REAL_TYPE>(fabs(LU_Diag[k]),min_diag);
						// Update values on the whole diagonal with L and U
						//DiagonalUpdate(k, LU_Diag, LineIndecesL, LineValuesL, LineIndecesU, LineValuesU);
#if defined(USE_OMP_FACT)
#pragma omp parallel sections num_threads(2)
#endif
						{
							// reconstruct U-part from linked list
#if defined(USE_OMP_FACT)
#pragma omp section
#endif
							{
								Unorm = 0;
								Ui = LineIndecesU[k];
								while (Ui != EOL)
								{
									Unorm += LineValuesU[Ui] * LineValuesU[Ui];
									Ui = LineIndecesU[Ui];
								}
								Unorm = sqrt(Unorm);
								U_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(U_Entries.size());
#if defined(ILUC2)
								U2_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(U2_Entries.size());
#endif
								Ui = LineIndecesU[k];
								while (Ui != EOL)
								{
									u = fabs(LineValuesU[Ui]);
									if (u * NuU > tau * Unorm) // apply dropping rule
										U_Entries.push_back(Sparse::Row::make_entry(Ui, LineValuesU[Ui]));
#if defined(ILUC2)
									else if (u * NuU > tau2 * Unorm)
										U2_Entries.push_back(Sparse::Row::make_entry(Ui, LineValuesU[Ui]));
#endif
									else ndrops++;
									Ui = LineIndecesU[Ui];
								}
								U_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(U_Entries.size());
#if defined(ILUC2)
								U2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(U2_Entries.size());
#endif
								// Update estimator vectors
								if (estimator)
								{
									tlocal = Timer();
									//U-estimator
									EstimatorUpdate(k, LineIndecesU, LineValuesU, EstU1, mu_update_U1);
									EstimatorUpdate(k, LineIndecesU, LineValuesU, EstU2, mu_update_U2);
									testimator += Timer() - tlocal;
								}
								// Cleanup structures
								ClearList(cbeg, LineIndecesU, LineValuesU);
								// Insert indexes for transposed U-part traversal
								//~ PrepareTranspose(k,k+1,U_Address,U_Entries,TransposeU);
								if (U_Address[k].Size() != 0)
								{
									INMOST_DATA_ENUM_TYPE  i = U_Entries[U_Address[k].first].first;
									assert(static_cast<INMOST_DATA_INTEGER_TYPE>(i) > k);
									Ulist[k] = Ubeg[i];
									Ubeg[i] = k;
								}
								// Shift indexes for transposed U-part traversal for next iteration
								{
									INMOST_DATA_ENUM_TYPE i = Ubeg[k];
									while (i != EOL)
									{
										U_Address[i].first++;
										INMOST_DATA_ENUM_TYPE Li = Ulist[i], Ui;
										if (U_Address[i].Size() > 0)
										{
											Ui = U_Entries[U_Address[i].first].first;
											Ulist[i] = Ubeg[Ui];
											Ubeg[Ui] = i;
										}
										i = Li;
									}
								}
#if defined(ILUC2)
								// Insert indexes for transposed second-order U-part traversal
								//~ PrepareTranspose(k,k+1,U2_Address,U2_Entries,TransposeU2);
								if (U2_Address[k].Size() != 0)
								{
									INMOST_DATA_ENUM_TYPE i = U2_Entries[U2_Address[k].first].first;
									assert(static_cast<INMOST_DATA_INTEGER_TYPE>(i) > k);
									U2list[k] = U2beg[i];
									U2beg[i] = k;
								}
								// Shift indexes for transposed second-order U-part traversal
								{
									INMOST_DATA_ENUM_TYPE i = U2beg[k];
									while (i != EOL)
									{
										U2_Address[i].first++;
										INMOST_DATA_ENUM_TYPE Li = U2list[i], Ui;
										if (U2_Address[i].Size() > 0)
										{
											Ui = U2_Entries[U2_Address[i].first].first;
											U2list[i] = U2beg[Ui];
											U2beg[Ui] = i;
										}
										i = Li;
									}
								}
#endif
							}
							// reconstruct L-part from linked list
#if defined(USE_OMP_FACT)
#pragma omp section
#endif
							{
								Lnorm = 0;
								Li = LineIndecesL[k];
								while (Li != EOL)
								{
									Lnorm += LineValuesL[Li] * LineValuesL[Li];
									Li = LineIndecesL[Li];
								}
								Lnorm = sqrt(Lnorm);
								//insert column to L part
								L_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(L_Entries.size());
#if defined(ILUC2)
								L2_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(L2_Entries.size());
#endif
								Li = LineIndecesL[k];
								while (Li != EOL)
								{
									u = fabs(LineValuesL[Li]);
									if (u * NuL > tau * Lnorm)
										L_Entries.push_back(Sparse::Row::make_entry(Li, LineValuesL[Li]));
#if defined(ILUC2)
									else if (u * NuL > tau2 * Lnorm)
										L2_Entries.push_back(Sparse::Row::make_entry(Li, LineValuesL[Li]));
#endif
									else ndrops++;
									Li = LineIndecesL[Li];
								}
								L_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(L_Entries.size());
#if defined(ILUC2)
								L2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(L2_Entries.size());
#endif
								// Update estimator vectors
								if (estimator)
								{
									tlocal = Timer();
									//L-estimator
									EstimatorUpdate(k, LineIndecesL, LineValuesL, EstL1, mu_update_L1);
									EstimatorUpdate(k, LineIndecesL, LineValuesL, EstL2, mu_update_L2);
									testimator += Timer() - tlocal;
								}
								// Cleanup structures
								ClearList(cbeg, LineIndecesL, LineValuesL);
								// Insert indexes for transposed L-part traversal
								//~ PrepareTranspose(k,k+1,L_Address,L_Entries,TransposeU);
								if (L_Address[k].Size() > 0)
								{
									INMOST_DATA_ENUM_TYPE i = L_Entries[L_Address[k].first].first;
									assert(static_cast<INMOST_DATA_INTEGER_TYPE>(i) > k);
									Llist[k] = Lbeg[i];
									Lbeg[i] = k;
								}
								// Shift indexes for transposed L-part traversal for next iteration
								{
									INMOST_DATA_ENUM_TYPE  i = Lbeg[k];
									while (i != EOL)
									{
										L_Address[i].first++;
										INMOST_DATA_ENUM_TYPE Li = Llist[i], Ui;
										if (L_Address[i].Size() > 0)
										{
											Ui = L_Entries[L_Address[i].first].first;
											Llist[i] = Lbeg[Ui];
											Lbeg[Ui] = i;
										}
										i = Li;
									}
								}
#if defined(ILUC2)
								// Insert indexes for transposed second-order L-part traversal
								//~ PrepareTranspose(k,k+1,L2_Address,L2_Entries,TransposeL2);
								if (L2_Address[k].Size() != 0)
								{
									INMOST_DATA_ENUM_TYPE i = L2_Entries[L2_Address[k].first].first;
									assert(static_cast<INMOST_DATA_INTEGER_TYPE>(i) > k);
									L2list[k] = L2beg[i];
									L2beg[i] = k;
								}
								// Shift indexes for transposed second-order L-part traversal
								{
									INMOST_DATA_ENUM_TYPE i = L2beg[k];
									while (i != EOL)
									{
										L2_Address[i].first++;
										INMOST_DATA_ENUM_TYPE Li = L2list[i], Ui;
										if (L2_Address[i].Size())
										{
											Ui = L2_Entries[L2_Address[i].first].first;
											L2list[i] = L2beg[Ui];
											L2beg[Ui] = i;
										}
										i = Li;
									}
								}
#endif
							}
						}
					}
					// Shift indexes for transposed B matrix traversal
					Li = Bbeg[k];
					while (Li != EOL)
					{
						A_Address[Li].first++;
						Ui = Blist[Li];
						if (A_Address[Li].Size() > 0)
						{
							i = A_Entries[A_Address[Li].first].first;
							if (i < Li)
							{
								if (Bbeg[i] > Li)
								{
									Blist[Li] = Bbeg[i];
									Bbeg[i] = Li;
								}
								else
								{
									curr = next = Bbeg[i];
									while (next < Li)
									{
										curr = next;
										next = Blist[curr];
									}
									assert(curr < Li);
									assert(Li < next);
									Blist[curr] = Li;
									Blist[Li] = next;
								}
							}
						}
						Li = Ui;
					}
					//number of nonzeros
					nzLU += U_Address[k].Size() + L_Address[k].Size() + 1;
#if defined(ILUC2)
					nzLU2 += U2_Address[k].Size() + L2_Address[k].Size();
#endif
					if( verbosity )
					{
						if( k % 100 == 0 )
						{
							std::ios save(NULL);
							save.copyfmt(std::cout);
							if (verbosity == 1)
								std::cout << cend << "/" << moend << " factor " << std::setw(6) << std::fixed << std::setprecision(2) << 100.0f*(k - cbeg) / (float)(cend - cbeg) << "\r" << std::flush;
							else 
							{
								std::cout << std::setw(6) << std::fixed << std::setprecision(2) << 100.0f*(k - cbeg) / (float)(cend - cbeg);
								std::cout << " nnz LU " << std::setw(8) << nzLU;
								std::cout << " LU2 " << std::setw(8) << nzLU2;
								std::cout << " cond L " << std::setprecision(6) << std::setw(10) << NuL;
								std::cout << " D " << std::setw(10) << NuD;
								std::cout << " U " << std::setw(10) << NuU;
								std::cout << " swaps " << std::setw(10) << swaps;
								std::cout << " drops " << std::setw(10) << ndrops;
								std::cout << "\r" << std::flush;
							}
							std::cout.copyfmt(save);
						}
					}
					// iteration done
				}
				if( verbosity > 1 )
				{
					std::cout << "size " << cend-cbeg;
					std::cout << " total nonzeros in A " << nzA << " (sparsity " << nzA/(double)(cend-cbeg)/(double)(cend-cbeg)*100.0 << "%) in LU " << nzLU << " (fill " << nzLU/(double)nzA0*100.0 << "%)";
					std::cout << " in LU2 " << nzLU2;
					std::cout << " conditions L " << NuL_max << " D " << NuD << " U " << NuU_max << " pivot swaps " << swaps << " drops " << ndrops << "            " << std::endl;
				}
				//restore indexes
				A_Address[cbeg].first = Bstart[cbeg];
				U_Address[cbeg].first = U_Beg;
				L_Address[cbeg].first = L_Beg;
#if defined(ILUC2)
				U2_Address[cbeg].first = 0;
				L2_Address[cbeg].first = 0;
#endif
				for (k = cbeg + 1; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
				{
					A_Address[k].first = Bstart[k];
					U_Address[k].first = U_Address[k - 1].last;
					L_Address[k].first = L_Address[k - 1].last;
#if defined(ILUC2)
					U2_Address[k].first = U2_Address[k - 1].last;
					L2_Address[k].first = L2_Address[k - 1].last;
#endif
				}
			}

			
			/// FACTORIZATION COMPLETE
			tlfactor = Timer() - tt;
			tfactor += tlfactor;
			if( verbosity > 1 )
				std::cout << "Time " << tlfactor << "\n";

			if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
			
			//reset the scaling vectors
			for (k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) DL[k] = DR[k] = 1.0;
			//  Check that we have to enter the next level due to pivoting
			if( swaps )
			{
				tlocal = Timer();
				cend = wend-swaps;

				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

				std::vector<Sparse::Row::entry> C_Entries;
				interval<INMOST_DATA_ENUM_TYPE, Interval> C_Address(cend, wend);




				tt = Timer();
				if( verbosity > 1 )
				{
					std::cout << "Total swaps: " << swaps << " interval: " << cend << " " << wend << std::endl;
					std::cout << "Reassemble pivots\n";
				}

				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

				level_size.push_back(cend - wbeg);
				i = wbeg;
				//enumerate entries that we keep first
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k) if( !Pivot[k] )
				{
					localP[k] = i;
					localQ[k] = i;
					++i;
				}
				//enumerate entries that we dropped at the end
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k) if( Pivot[k] )
				{
					localP[k] = i;
					localQ[k] = i;
					++i;
				}
				
				for (k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					U[localP[k]] = DL0[k];
					V[localQ[k]] = DR0[k];
				}
				for (k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					DL0[k] = U[k];
					DR0[k] = V[k];
				}
				
				//inverse the ordering
				ReorderEF(wbeg, wend, donePQ, localP, localQ);
				inversePQ(wbeg,wend,localP,localQ, invP,invQ);
				// Reassamble matrix
				//here F is stored by columns and E by rows
				E_Address.push_back(new interval<INMOST_DATA_ENUM_TYPE, Interval>(cend, wend));
				F_Address.push_back(new interval<INMOST_DATA_ENUM_TYPE, Interval>(cend, wend));
				//find out the space required for each column of F
				//and setup the first and last pointers
				for(k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
					F_Address.back()->at(k).first = F_Address.back()->at(k).last = static_cast<INMOST_DATA_ENUM_TYPE>(F_Entries.size());
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					j = invP[k];
					if( k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend) )
					{
						for (INMOST_DATA_ENUM_TYPE jt = A_Address[j].first; jt < A_Address[j].last; ++jt)
						{
							i = localQ[A_Entries[jt].first];
							if( i >= cend ) //put into F
								F_Address.back()->at(i).last++;
						}
					}
				}
				//Establish the addresses
				INMOST_DATA_ENUM_TYPE offset = 0;
				for(k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					F_Address.back()->at(k).first += offset;
					F_Address.back()->at(k).last  += offset;
					offset += F_Address.back()->at(k).Size();
				}
				F_Entries.resize(F_Address.back()->at(wend-1).last);
				//Reset addresses to advance current fill position
				for(k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
					F_Address.back()->at(k).last = F_Address.back()->at(k).first;
				//reassamble lower-right corner of B into A
				nnzF = nnzE = 0;
				//put lower-left into E and upper-right into F
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					j = invP[k];
					if( k >= static_cast<INMOST_DATA_INTEGER_TYPE>(cend) )
					{
						E_Address.back()->at(k).first = static_cast<INMOST_DATA_ENUM_TYPE>(E_Entries.size());
						C_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(C_Entries.size());
						for (INMOST_DATA_ENUM_TYPE jt = A_Address[j].first; jt < A_Address[j].last; ++jt)
						{
							i = localQ[A_Entries[jt].first];
							u = A_Entries[jt].second;
							if( i >= cend ) //put into A
								C_Entries.push_back(Sparse::Row::make_entry(i, u));
							else //put into E
							{
								nnzE++;
								E_Entries.push_back(Sparse::Row::make_entry(i, u));
							}
						}
						C_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(C_Entries.size());
						E_Address.back()->at(k).last = static_cast<INMOST_DATA_ENUM_TYPE>(E_Entries.size());
						std::sort(C_Entries.begin() + C_Address[k].first, C_Entries.end());
						std::sort(E_Entries.begin() + E_Address.back()->at(k).first, E_Entries.end());
					}
					else
					{
						for (INMOST_DATA_ENUM_TYPE jt = A_Address[j].first; jt < A_Address[j].last; ++jt)
						{
							i = localQ[A_Entries[jt].first];
							u = A_Entries[jt].second;
							if( i >= cend ) //put into F
							{
								nnzF++;
								F_Entries[F_Address.back()->at(i).last++] = Sparse::Row::make_entry(k, u);
							}
						}
					}
				}

				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

				A_Address.clear();
				A_Entries.clear();
				A_Entries.swap(std::vector<Sparse::Row::entry>()); //release memory

				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
				//  Reassamble L and U
				//first change positions
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					j = invP[k];
					//if( j < cend )
					{
						INMOST_DATA_ENUM_TYPE cnt;
						// L-part
						cnt = 0;
						for (INMOST_DATA_ENUM_TYPE jt = L_Address[j].first; jt < L_Address[j].last; ++jt)
						{
							i = localQ[L_Entries[jt].first];
							if( i < cend )
								L_Entries[jt].first = i;
							else
							{
								L_Entries[jt].first = ENUMUNDEF;
								L_Entries[jt].second = 0.0;
								cnt++;
							}
						}
						//this will move all undefined to the end
						std::sort(L_Entries.begin()+L_Address[j].first,L_Entries.begin()+L_Address[j].last);
						//this will remove references to undefined
						L_Address[j].last -= cnt;
						// U-part
						cnt = 0;
						for (INMOST_DATA_ENUM_TYPE jt = U_Address[j].first; jt < U_Address[j].last; ++jt)
						{
							Ui = U_Entries[jt].first;
							i = localQ[Ui];
							if( i < cend )
								U_Entries[jt].first = i;
							else
							{
								U_Entries[jt].first = ENUMUNDEF;
								U_Entries[jt].second = 0.0;
								cnt++;
							}
						}
						//this will move all undefined to the end
						std::sort(U_Entries.begin()+U_Address[j].first,U_Entries.begin()+U_Address[j].last);
						//this will remove references to undefined
						U_Address[j].last -= cnt;
					}
				}
#if defined(ILU2)
				//now the same for second-order LU
				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					j = invP[k];
					//if( j < cend )
					{
						INMOST_DATA_ENUM_TYPE cnt;
						//L-part
						cnt = 0;
						for (INMOST_DATA_ENUM_TYPE jt = L2_Address[j].first; jt < L2_Address[j].last; ++jt)
						{
							i = localQ[L2_Entries[jt].first];
							if( i < cend )
								L2_Entries[jt].first = i;
							else
							{
								L2_Entries[jt].first = ENUMUNDEF;
								L2_Entries[jt].second = 0.0;
								cnt++;
							}
						}
						//this will move all undefined to the end
						std::sort(L2_Entries.begin()+L2_Address[j].first,L2_Entries.begin()+L2_Address[j].last);
						//this will remove references to undefined
						L2_Address[j].last -= cnt;
						//U-part
						cnt = 0;
						for (INMOST_DATA_ENUM_TYPE jt = U2_Address[j].first; jt < U2_Address[j].last; ++jt)
						{
							i = localQ[U2_Entries[jt].first];
							if( i < cend )
								U2_Entries[jt].first = i;
							else
							{
								U2_Entries[jt].first = ENUMUNDEF;
								U2_Entries[jt].second = 0.0;
								cnt++;
							}
						}
						//this will move all undefined to the end
						std::sort(U2_Entries.begin()+U2_Address[j].first,U2_Entries.begin()+U2_Address[j].last);
						//this will remove references to undefined
						U2_Address[j].last -= cnt;
					}
				}
#endif
				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
				//first establish intervals in arrays
				{
					interval<INMOST_DATA_ENUM_TYPE, Interval> tmpL(wbeg,cend), tmpU(wbeg,cend);
#if defined(ILU2)
					interval<INMOST_DATA_ENUM_TYPE, Interval> tmpL2(wbeg, cend), tmpU2(wbeg, cend)
#endif
					interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> tmpD(wbeg,cend);
					for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
					{
						j = invP[k];
						tmpL[k] = L_Address[j];
						tmpU[k] = U_Address[j];
#if defined(ILU2)
						tmpL2[k] = L2_Address[j];
						tmpU2[k] = U2_Address[j];
#endif
						tmpD[k] = LU_Diag[j];
					}
					for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
					{
						L_Address[k] = tmpL[k];
						U_Address[k] = tmpU[k];
#if defined(ILU2)
						L2_Address[k] = tmpL2[k];
						U2_Address[k] = tmpU2[k];
#endif
						LU_Diag[k] = tmpD[k];
					}
				}
				
				applyPQ(wbeg, wend, localP, localQ, invP, invQ);
				tt = Timer()-tt;
				if( verbosity > 1 )
					std::cout << "Time " << tt << "\n";

				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

				int ndrops_lf = 0, ndrops_eu = 0, ndrops_s = 0;
///////////////////////////////////////////////////////////////////////////////////
//  Compute Schur                                                                //
///////////////////////////////////////////////////////////////////////////////////
				//result of multilevel preconditioner
				//        |LDU  F |
				//  A  =  |       |
				//        | E   A |
				//
				// For the next step of multilevel factorization we have to compute
				// S = C - (E U^{-1}) D^{-1} (L^{-1} F)
				//
				//first precompute LF block
				{
					if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

					std::vector< interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> > LineValuesS(nthreads, interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>(mobeg, moend, 0.0));
					std::vector< interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> > LineIndecesS(nthreads, interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>(mobeg, moend + 1, UNDEF));
					std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > indicesS(nthreads);

					std::vector<Sparse::Row::entry>  LFt_Entries, EU_Entries, S_Entries;
					interval<INMOST_DATA_ENUM_TYPE, Interval> LFt_Address(wbeg, cend), EU_Address(cend, wend), S_Address(cend, wend);

					{
						std::vector<Sparse::Row::entry> LF_Entries;
						interval<INMOST_DATA_ENUM_TYPE, Interval> LF_Address(cend, wend);

						tt = Timer();
						if (verbosity > 1)
							std::cout << "Assemble LF\n";

						if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

#if defined(USE_OMP_FACT)
#pragma omp parallel for private(k, Li, u, j, curr, next) reduction(+:ndrops_lf)
#endif
						for (k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
						{
							int nthr = Thread();
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& LineValues = LineValuesS[nthr];
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& LineIndeces = LineIndecesS[nthr];
							INMOST_DATA_REAL_TYPE Fmax = 0;//Fnorm = 0, Fnum = 0;
							Li = cbeg;
							for (INMOST_DATA_ENUM_TYPE jt = F_Address.back()->at(k).first; jt < F_Address.back()->at(k).last; ++jt)
							{
								u = F_Entries[jt].second;
								j = F_Entries[jt].first;
								LineValues[j] = u;
								Li = LineIndeces[Li] = j + 1;
								//Fnorm += u*u;
								//Fnum  ++;
								Fmax = std::max<INMOST_DATA_REAL_TYPE>(Fmax, fabs(u));
							}
							LineIndeces[Li] = EOL;
							//~ if( Fnum ) Fnorm = sqrt(Fnorm/Fnum);
							//~ LFdrop = 0;
							// perform solve with L part
							//compute ~F_i = L^{-1} F_i
							Li = LineIndeces[cbeg];
							//Sbeg = cbeg;
							while (Li != EOL)
							{
								curr = Li;
								if (1 + LineValues[Li - 1] != 1)
								{
									double dtol = tau2 * Fmax / NuL_max;
									for (INMOST_DATA_ENUM_TYPE ru = L_Address[Li - 1].first; ru < L_Address[Li - 1].last; ++ru)
									{
										u = LineValues[Li - 1] * L_Entries[ru].second;
										j = L_Entries[ru].first;
										if (LineIndeces[j + 1] != UNDEF) // There is an entry in the list
											LineValues[j] -= u;
#if defined(SCHUR_DROPPING_LF_PREMATURE)
										else if (fabs(u) > dtol * LU_Diag[j])
#else
										else if (1 + u != 1)
#endif
										{

											next = curr;
											while (next < j + 1)
											{
												curr = next;
												next = LineIndeces[curr];
											}
											assert(curr < j + 1);
											assert(j + 1 < next);
											LineIndeces[curr] = j + 1;
											LineIndeces[j + 1] = next;
											LineValues[j] = -u;

											//LineValues[j] = -u;
											//LineIndeces[j + 1] = Sbeg;
											//Sbeg = j+1;
										}
										//else LFdrop += u/LU_Diag[j];
									}
#if defined(ILUC2) && defined(ILUC2_SCHUR)
									// perform solve with second-order L part
									curr = Li;
									for (INMOST_DATA_ENUM_TYPE ru = L2_Address[Li - 1].first; ru < L2_Address[Li - 1].last; ++ru)
									{
										u = LineValues[Li - 1] * L2_Entries[ru].second;
										j = L2_Entries[ru].first;
										if (LineIndeces[j + 1] != UNDEF) // There is an entry in the list
											LineValues[j] -= u;
#if defined(SCHUR_DROPPING_LF_PREMATURE)
										else if (fabs(u) > dtol * LU_Diag[j])
#else
										else if (1 + u != 1)
#endif
										{

											next = curr;
											while (next < j + 1)
											{
												curr = next;
												next = LineIndeces[curr];
											}
											assert(curr < j + 1);
											assert(j + 1 < next);
											LineIndeces[curr] = j + 1;
											LineIndeces[j + 1] = next;
											LineValues[j] = -u;

											//LineValues[j] = -u;
											//LineIndeces[j + 1] = Sbeg;
											//Sbeg = j+1;
										}
										//else LFdrop += u/LU_Diag[j];
									}
#endif
								}
								// solve iteration done
								Li = LineIndeces[Li];
							}
							// Rescale by diagonal
							Li = LineIndeces[cbeg];
							while (Li != EOL)
							{
								LineValues[Li - 1] /= LU_Diag[Li - 1]; //Li - 1 or k?
								Li = LineIndeces[Li];
							}
							// Compute norm for dropping
#if defined(SCHUR_DROPPING_LF)
							INMOST_DATA_REAL_TYPE LFnorm, LFnum, LFmax, LFmin, LFtau;
							LFnorm = LFnum = 0;
							Li = LineIndeces[cbeg];
							LFmax = 0;
							LFmin = 1.0e+54;
							while (Li != EOL)
							{
								u = fabs(LineValues[Li - 1]);
								LFnorm += u * u;
								if (u > LFmax) LFmax = u;
								if (u < LFmin) LFmin = u;
								LFnum++;
								Li = LineIndeces[Li];
							}
							if (LFnum) LFnorm = sqrt(LFnorm / LFnum);
							//LFnorm = std::min(1.0,LFnorm);
							//LFtau = LFmin + (LFmax - LFmin)*std::min(tau*tau,tau2) / NuL_max;
							//LFtau = LFmin + std::min(tau*tau,tau2)*LFnorm;// / NuL_max;
							LFtau = LFmin + (LFmax - LFmin) * pow(LFnorm / LFmax, 2);
							//LFtau = LFmin + (LFmax - LFmin)*pow(LFnorm/LFmax,4);
#endif
						// Assemble column into matrix
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
							{
								Li = LineIndeces[cbeg];
								LF_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LF_Entries.size());
								while (Li != EOL)
								{
									u = LineValues[Li - 1];
									assert(Li - 1 >= wbeg && Li - 1 < cend);
#if defined(SCHUR_DROPPING_LF)
									if (fabs(u) >= LFtau)
#else
									if (1 + u != 1)
#endif
										LF_Entries.push_back(Sparse::Row::make_entry(Li - 1, u));
									else
									{
										//LFdrop += u;
										ndrops_lf++;
									}
									Li = LineIndeces[Li];
								}
								LF_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LF_Entries.size());
							}
							// if( LFdrop )
							// {
							// 	LFdrop /= (double) LF_Address[k].Size();
							// 	for(j = LF_Address[k].first; j < LF_Address[k].last; ++j)
							// 		LF_Entries[j].second += LFdrop;
							// }
							//assert(std::is_sorted(LF_Entries.begin()+LF_Address[k].first,LF_Entries.end()));
							// clean linked list
							Li = LineIndeces[cbeg];
							while (Li != EOL)
							{
								j = LineIndeces[Li];
								LineIndeces[Li] = UNDEF;
								LineValues[Li - 1] = 0.0;
								Li = j;
							}
							LineIndeces[cbeg] = UNDEF;
							if (verbosity)
							{
								if (k % 100 == 0)
								{
									std::ios save(NULL);
									save.copyfmt(std::cout);
									std::cout << "LF " << std::setw(6) << std::fixed << std::setprecision(2) << 100.f * (k - cend + 1) / (1.f * (wend - cend));
									std::cout << " nnz " << std::setw(10) << LF_Entries.size() << " drops " << std::setw(10) << ndrops_lf;
									std::cout << "\r" << std::flush;
									std::cout.copyfmt(save);
								}
							}
							// iteration done!
						}
						tt = Timer() - tt;
						if (verbosity > 1)
						{
							std::cout << "F nnz " << nnzF << " LF nnz " << LF_Entries.size() << " drops " << ndrops_lf << "         \t\t\n";
							std::cout << "Time " << tt << "\n";
						}

						if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

						// prepearing LF block for transposed traversal
						{
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Flist(cend, wend);
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Fbeg(wbeg, cend);
							tt = Timer();
							if (verbosity > 1)
								std::cout << "Transpose LF\n";
							//std::fill(Fbeg.begin() + wbeg - mobeg, Fbeg.begin() + cend - mobeg, EOL);
							for (k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) Fbeg[k] = EOL;
							for (k = wend; k > static_cast<INMOST_DATA_INTEGER_TYPE>(cend); --k)
							{
								if (LF_Address[k - 1].Size() > 0)
								{
									i = LF_Entries[LF_Address[k - 1].first].first;
									Flist[k - 1] = Fbeg[i];
									Fbeg[i] = k - 1;
								}
								else Flist[k - 1] = EOL;
							}
							//transpone LF matrix
							LFt_Entries.resize(LF_Entries.size());

							if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

							j = 0;
							for (k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); k++)
							{
								// assemble line of transposed LF
								LFt_Address[k].first = j;
								Li = Fbeg[k];
								//Fnorms[k] = 0;
								while (Li != EOL)
								{
									LFt_Entries[j].first = Li;
									LFt_Entries[j].second = LF_Entries[LF_Address[Li].first].second;
									j++;
									Li = Flist[Li];
								}
								LFt_Address[k].last = j;
								// shift indices for transposed traversal
								Li = Fbeg[k];
								while (Li != EOL)
								{
									Ui = Flist[Li];
									LF_Address[Li].first++;
									if (LF_Address[Li].Size() > 0)
									{
										i = LF_Entries[LF_Address[Li].first].first;
										Flist[Li] = Fbeg[i];
										Fbeg[i] = Li;
									}
									Li = Ui;
								}
								// iteration done
							}
						}
						tt = Timer() - tt;
						if (verbosity > 1)
							std::cout << "Time " << tt << "\n";
					}

					if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

					//DumpMatrix(LFt_Address,LFt_Entries,wbeg,wend,"LFt.mtx");
					// EU and Schur
					tt = Timer();
					if (verbosity > 1)
						std::cout << "Assemble EU\n";
#if defined(USE_OMP_FACT)
#pragma omp parallel for private(k, Li, u, i, j, curr, next) reduction(+:ndrops_eu)
#endif
					for (k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
					{
						int nthr = Thread();
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& LineValues = LineValuesS[nthr];
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& LineIndeces = LineIndecesS[nthr];
						double Emax = 0;//Enorm = 0, Enum = 0;
						// no values at E block - unpack C block into Schur 
						Li = cbeg;
						for (j = E_Address.back()->at(k).first; j < E_Address.back()->at(k).last; ++j) //iterate over values of k-th row of E
						{
							u = E_Entries[j].second;
							i = E_Entries[j].first;
							LineValues[i] = u;
							Li = LineIndeces[Li] = i + 1;
							//~ Enorm += u*u;
							//~ Enum  ++;
							Emax = std::max(Emax, fabs(u));
						}
						LineIndeces[Li] = EOL;
						//~ if( Enum ) Enorm = sqrt(Enorm/Enum);
						//~ EUdrop = 0;
						// perform solve with U part
						// compute ~E_i = E_i U^{-1}
						Li = LineIndeces[cbeg];
						while (Li != EOL)
						{
							curr = Li;
							if (1 + LineValues[Li - 1] != 1)
							{
								double dtol = tau2 * Emax / NuU_max;
								for (INMOST_DATA_ENUM_TYPE ru = U_Address[Li - 1].first; ru < U_Address[Li - 1].last; ++ru)
								{
									u = LineValues[Li - 1] * U_Entries[ru].second;
									j = U_Entries[ru].first;
									if (LineIndeces[j + 1] != UNDEF) // There is an entry in the list
										LineValues[j] -= u;
#if defined(SCHUR_DROPPING_EU_PREMATURE)
									else if (fabs(u) > dtol * LU_Diag[j])
#else
									else if (1 + u != 1)
#endif
									{
										next = curr;
										while (next < j + 1)
										{
											curr = next;
											next = LineIndeces[curr];
										}
										assert(curr < j + 1);
										assert(j + 1 < next);
										LineIndeces[curr] = j + 1;
										LineIndeces[j + 1] = next;
										LineValues[j] = -u;
									}
									//else EUdrop += u/LU_Diag[j];
								}
#if defined(ILUC2) && defined(ILUC2_SCHUR)
								// perform solve with second-order U part
								curr = Li;
								for (INMOST_DATA_ENUM_TYPE ru = U2_Address[Li - 1].first; ru < U2_Address[Li - 1].last; ++ru)
								{
									u = LineValues[Li - 1] * U2_Entries[ru].second;
									j = U2_Entries[ru].first;
									if (LineIndeces[j + 1] != UNDEF) // There is an entry in the list
										LineValues[j] -= u;
#if defined(SCHUR_DROPPING_EU_PREMATURE)
									else if (fabs(u) > dtol * LU_Diag[j])
#else
									else if (1 + u != 1)
#endif
									{
										next = curr;
										while (next < j + 1)
										{
											curr = next;
											next = LineIndeces[curr];
										}
										assert(curr < j + 1);
										assert(j + 1 < next);
										LineIndeces[curr] = j + 1;
										LineIndeces[j + 1] = next;
										LineValues[j] = -u;
									}
									//else EUdrop += u/LU_Diag[j];
								}
#endif
							}
							Li = LineIndeces[Li];
							// solve iteration done
						}
						// rescale vector by diagonal *E = ~E_i D^{-1} = E_i U^{-1} D^{-1}
						Li = LineIndeces[cbeg];
						while (Li != EOL)
						{
							LineValues[Li - 1] /= LU_Diag[Li - 1];
							Li = LineIndeces[Li];
						}
						// compute norm
#if defined(SCHUR_DROPPING_EU)
						INMOST_DATA_REAL_TYPE EUnorm, EUnum, EUmax, EUmin, EUtau;
						EUnorm = EUnum = 0;
						EUmax = 0;
						EUmin = 1.0e+54;
						Li = LineIndeces[cbeg];
						while (Li != EOL)
						{
							u = fabs(LineValues[Li - 1]);
							EUnorm += u * u;
							if (u > EUmax) EUmax = u;
							if (u < EUmin) EUmin = u;
							EUnum++;
							Li = LineIndeces[Li];
						}
						if (EUnum) EUnorm = sqrt(EUnorm / EUnum);
						//~ EUtau = EUmin + std::min(tau*tau,tau2)*EUnorm;// / NuU_max;
						//EUnorm = std::min(1.0,EUnorm);
						EUtau = EUmin + (EUmax - EUmin) * pow(EUnorm / EUmax, 2);
						//EUtau = EUmin + (EUmax - EUmin)*pow(EUnorm/EUmax,4);
#endif
						// drop values that do not satisfy tolerances from linked list of line of EU block
						// Assemble column into matrix
						Li = LineIndeces[cbeg];
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
						{
							EU_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(EU_Entries.size());
							while (Li != EOL)
							{
								u = LineValues[Li - 1];
								assert(Li - 1 >= wbeg && Li - 1 < cend);
#if defined(SCHUR_DROPPING_EU)
								if (fabs(u) >= EUtau)
#else
								if (1 + u != u)
#endif
									EU_Entries.push_back(Sparse::Row::make_entry(Li - 1, u));
								else
								{
									//EUdrop += u;
									ndrops_eu++;
								}
								Li = LineIndeces[Li];
							}
							EU_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(EU_Entries.size());
						}
						// if( EUdrop )
						// {
						// 	EUdrop /= (double)EU_Address[k].Size();
						// 	for(INMOST_DATA_ENUM_TYPE r = EU_Address[k].first; r < EU_Address[k].last; ++r)
						// 		EU_Entries[r].second += EUdrop;
						// }
	///////////////////////////////////////////////////////////////////////////////////
	//         clean up linked list                                                  //
	///////////////////////////////////////////////////////////////////////////////////
						Li = LineIndeces[cbeg];
						while (Li != EOL)
						{
							j = LineIndeces[Li];
							LineIndeces[Li] = UNDEF;
							LineValues[Li - 1] = 0.0;
							Li = j;
						}
						LineIndeces[cbeg] = UNDEF;
						if (verbosity)
						{
							if (k % 100 == 0)
							{
								std::ios save(NULL);
								save.copyfmt(std::cout);
								std::cout << "EU " << std::setw(6) << std::fixed << std::setprecision(2) << 100.f * (k - cend + 1) / (1.f * (wend - cend));
								std::cout << " nnz " << std::setw(10) << EU_Entries.size() << " drops " << std::setw(10) << ndrops_eu;
								std::cout << "\r" << std::flush;
								std::cout.copyfmt(save);
							}
						}
					}
					tt = Timer() - tt;
					if (verbosity > 1)
					{
						std::cout << "E nnz " << nnzE << " EU nnz " << EU_Entries.size() << " drops " << ndrops_eu << "          \t\t\n";
						std::cout << "Time " << tt << "\n";
					}

					if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

#if defined(SCHUR_DROPPING_S) || defined(SCHUR_DROPPING_S_PREMATURE)
					{
						///////////////////////////////////////////////////////////////////////////////////
						//         Compute column-norms of schur                                         //
						///////////////////////////////////////////////////////////////////////////////////
						tt = Timer();

						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Scolmax(mobeg, moend, 0.0), Srowmax(mobeg, moend, 0.0);
						{
							std::vector< interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> > Scolmax_local(nthreads, interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>(mobeg, moend, 0.0));

							if (verbosity > 1)
								std::cout << "Compute Schur column norm\n";

#if defined(USE_OMP_FACT)
#pragma omp parallel private(k) 
#endif
							{
								int nthr = Thread();
								interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& LineValues = LineValuesS[nthr];
								interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& LineIndeces = LineIndecesS[nthr];
								interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& Scolmaxloc = Scolmax_local[nthr];
#if defined(USE_OMP_FACT)
#pragma omp for
#endif
								for (k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
								{
									INMOST_DATA_REAL_TYPE v, coef;
									INMOST_DATA_ENUM_TYPE Sbeg = EOL, i, j;
									for (INMOST_DATA_ENUM_TYPE r = C_Address[k].first; r < C_Address[k].last; ++r)
									{
										LineValues[C_Entries[r].first] = C_Entries[r].second;
										LineIndeces[C_Entries[r].first] = Sbeg;
										Sbeg = C_Entries[r].first;
									}
									for (INMOST_DATA_ENUM_TYPE r = EU_Address[k].first; r < EU_Address[k].last; ++r)
									{
										coef = -EU_Entries[r].second * LU_Diag[EU_Entries[r].first];
										AddListUnordered(Sbeg, LFt_Address[EU_Entries[r].first], LFt_Entries, coef, LineIndeces, LineValues, 0);
									}
									i = Sbeg;
									while (i != EOL)
									{
										j = i;
										i = LineIndeces[i];
										v = fabs(LineValues[j]);
										Srowmax[k] = std::max(Srowmax[k], v);
										Scolmaxloc[j] = std::max(Scolmaxloc[j], v);
										LineValues[j] = 0;
										LineIndeces[j] = UNDEF;
									}
									if (verbosity)
									{
										if (k % 100 == 0)
										{
											std::ios save(NULL);
											save.copyfmt(std::cout);
											std::cout << "Schur column norm " << std::setw(6) << std::fixed << std::setprecision(2) << 100.f * (k - cend + 1) / (1.f * (wend - cend));
											std::cout << "\t\t\r" << std::flush;
											std::cout.copyfmt(save);
										}
									}
								}
#if defined(USE_OMP_FACT)
#pragma omp for
#endif
								for (k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
								{
									for (unsigned j = 0; j < Scolmax_local.size(); ++j)
										Scolmax[k] = std::max(Scolmax[k], Scolmax_local[j][k]);
								}
							}
						}

						tt = Timer() - tt;
						if (verbosity > 1)
							std::cout << "\nTime " << tt << "\n";

						if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
#endif //SCHUR_DROPPING_S
						// Construction of Schur complement 
						tt = Timer();
						if (verbosity > 1)
							std::cout << "Construct Schur\n";

						if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

						nnzA = 0;
						for (k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k) nnzA += C_Address[k].Size();
#if defined(USE_OMP_FACT)
#pragma omp parallel for private(k,Sbeg, Li, Ui, i,j, u,l/*, Smin,Smax,Snorm,Snum,Stau,Sdrop*/) reduction(+:ndrops_s)
#endif
						for (k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
						{
							int nthr = Thread();
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& LineValues = LineValuesS[nthr];
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& LineIndeces = LineIndecesS[nthr];
							std::vector<INMOST_DATA_ENUM_TYPE>& indices = indicesS[nthr];
							///////////////////////////////////////////////////////////////////////////////////
							//         unpack line of A matrix to temporary linked list                      //
							///////////////////////////////////////////////////////////////////////////////////
							Sbeg = EOL;
							for (INMOST_DATA_ENUM_TYPE r = C_Address[k].first; r < C_Address[k].last; ++r)
							{
								LineValues[C_Entries[r].first] = C_Entries[r].second;
								LineIndeces[C_Entries[r].first] = Sbeg;
								Sbeg = C_Entries[r].first;
							}
///////////////////////////////////////////////////////////////////////////////////
//        multiply EU and LF blocks and add to temporary linked list             //
///////////////////////////////////////////////////////////////////////////////////
							for (INMOST_DATA_ENUM_TYPE r = EU_Address[k].first; r < EU_Address[k].last; ++r)
							{
								//iterate over corresponding row of LF, add multiplication to list
								l = -EU_Entries[r].second * LU_Diag[EU_Entries[r].first];
								AddListUnordered(Sbeg, LFt_Address[EU_Entries[r].first], LFt_Entries, l, LineIndeces, LineValues, 0);
							}
							OrderList(Sbeg, LineIndeces, indices);
							Ui = Sbeg;
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
							{
								S_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(S_Entries.size());
								while (Ui != EOL)
								{
									u = LineValues[Ui];
#if defined(SCHUR_DROPPING_S)
									if (fabs(u) >= std::min(Srowmax[k], Scolmax[Ui]) * tau2)
#else
									if (1 + u != 1)
#endif
										S_Entries.push_back(Sparse::Row::make_entry(Ui, u));
									else ndrops_s++;
									Ui = LineIndeces[Ui];
								}
								S_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(S_Entries.size());
							}
							ClearList(Sbeg, LineIndeces, LineValues);
							if (verbosity)
							{
								if (k % 100 == 0)
								{
									std::ios save(NULL);
									save.copyfmt(std::cout);
									std::cout << "Schur " << std::setw(6) << std::fixed << std::setprecision(2) << 100.f * (k - cend + 1) / (1.f * (wend - cend));
									std::cout << " nnz " << std::setw(10) << S_Entries.size() << " drop S " << ndrops_s;
									std::cout << "\t\t\r" << std::flush;
									std::cout.copyfmt(save);
								}
							}
						}
						///////////////////////////////////////////////////////////////////////////////////
						//         Schur complement row done                                             //
						///////////////////////////////////////////////////////////////////////////////////
					}

					tt = Timer() - tt;
					if (verbosity > 1)
					{
						std::cout << "Schur nnz " << S_Entries.size() << " drop S " << ndrops_s << "          \t\t\n";
						std::cout << "Time " << tt << "\n";
					}

					if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

					A_Entries.swap(S_Entries);
					A_Address.swap(S_Address);
				}
				// Cleanup arrays for the next iteration
				tt = Timer();
				if( verbosity > 1 )
					std::cout << "Cleanup\n";

				if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

				for(k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					localP[k] = localQ[k] = ENUMUNDEF;
					Bbeg[k] = EOL;
#if defined(ILU2)
					L2beg[k] = U2beg[k] = EOL;
#endif
					U[k] = V[k] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
					Pivot[k] = false;
				}
				for(k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					U_Address[k].first = U_Address[k].last = ENUMUNDEF;
					L_Address[k].first = L_Address[k].last = ENUMUNDEF;
#if defined(ILU2)
					U2_Address[k].first = U2_Address[k].last = ENUMUNDEF;
					L2_Address[k].first = L2_Address[k].last = ENUMUNDEF;
#endif
				}
#if defined(ILU2)
				U2_Entries.clear();
				L2_Entries.clear();
#endif
				if( /*wend-cend < 16 ||*/ A_Entries.empty() )
					block_pivot = true;
				//rescale_b = false;
				//run_mpt = false;
				tt = Timer()-tt;
				tlschur = Timer()-tlocal;
				tschur += tlschur;
				if( verbosity > 1 )
				{
					std::cout << "Time " << tt << "\n";
					std::cout << "Total Schur " << tlschur << "\n";
				}

				wbeg = cend; //there is more work to do
			}
			else
			{
				level_size.push_back(wend - wbeg);
				wbeg = wend;
			}
			totswaps += swaps;
			swaps = 0;
		}
		
		tt = Timer();
		
		if( rescale_b )
		{
			if( verbosity > 1 )
				std::cout << std::endl << " rescaling block B back " << std::endl;

			for (k = mobeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(moend); ++k)
			{
				LU_Diag[k] /= (DL0[k] * DR0[k]);
				for (INMOST_DATA_ENUM_TYPE r = U_Address[k].first; r < U_Address[k].last; ++r) U_Entries[r].second *= (DR0[k] / DR0[U_Entries[r].first]);
				for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r < L_Address[k].last; ++r) L_Entries[r].second *= (DL0[k] / DL0[L_Entries[r].first]);
			}
			
			//rescale EF
			if( !level_size.empty() )
			{
				INMOST_DATA_ENUM_TYPE first = mobeg, last;
				INMOST_DATA_ENUM_TYPE kbeg, kend, r;
				for(size_t level = 0; level < level_size.size(); ++level)
				{
					last = first + level_size[level];
					kbeg = last;//std::max(cbeg,last);
					kend = moend;//std::min(cend,moend);
					//std::cout << "Rescale EF back level " << level << " ";
					//std::cout << "size " << level_size[level] << " ";
					//std::cout << "interval [" << first << "," << last << "] ";
					//std::cout << "selected [" << kbeg << ":" << kend << "]" << std::endl;
					for (k = kbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(kend); ++k) 
					{
						if( F_Address[level]->at(k).Size() )
						{
							for (r = F_Address[level]->at(k).first; r < F_Address[level]->at(k).last; ++r)
								F_Entries[r].second /= DR0[k]*DL0[F_Entries[r].first];
						}
						if( E_Address[level]->at(k).Size() )
						{
							for (r = E_Address[level]->at(k).first; r < E_Address[level]->at(k).last; ++r)
								E_Entries[r].second /= DL0[k]*DR0[E_Entries[r].first];
						}
					}
					first = last;
				}
			}
		}
		/// RESCALING DONE
		tlrescale = Timer() - tt;
		trescale += tlrescale;
		
		level_interval.resize(level_size.size());
		if( !level_size.empty() )
		{
			level_interval[0].first = mobeg;
			level_interval[0].last = mobeg + level_size[0];
			for (k = 1; k < static_cast<INMOST_DATA_INTEGER_TYPE>(level_size.size()); k++)
			{
				level_interval[k].first = level_interval[k - 1].last;
				level_interval[k].last = level_interval[k].first + level_size[k];
			}
		}
		/// MULTILEVEL FACTORIZATION COMPLETE
		ttotal = Timer() - ttotal;
		if( verbosity > 1 )
		{
			std::cout << "total      " << ttotal       << "\n";
			std::cout << "reorder    " << treorder     << " ("; fmt(std::cout,100.0*treorder/ttotal    ,6,2) << "%)\n";
			std::cout << "   mpt     " << ttransversal << " ("; fmt(std::cout,100.0*ttransversal/ttotal,6,2) << "%)\n";
#if defined(REORDER_METIS_ND)
			std::cout << "   graph   " << tmetisgraph  << " ("; fmt(std::cout,100.0*tmetisgraph/ttotal ,6,2) << "%)\n";
			std::cout << "   nd      " << tmetisnd     << " ("; fmt(std::cout,100.0*tmetisnd/ttotal    ,6,2) << "%)\n";
#endif
#if defined(REORDER_RCM)
			std::cout << "   graph   " << trcmgraph    << " ("; fmt(std::cout,100.0*trcmgraph/ttotal   ,6,2) << "%)\n";
			std::cout << "   rcm     " << trcmorder    << " ("; fmt(std::cout,100.0*trcmorder/ttotal   ,6,2) << "%)\n";
#endif
			std::cout << "reassamble " << treassamble  << " ("; fmt(std::cout,100.0*treassamble/ttotal ,6,2) << "%)\n";
			std::cout << "rescale    " << trescale     << " ("; fmt(std::cout,100.0*trescale/ttotal    ,6,2) << "%)\n";
			std::cout << "factor     " << tfactor      << " ("; fmt(std::cout,100.0*tfactor/ttotal     ,6,2) << "%)\n";
			std::cout << "   schur   " << tschur       << " ("; fmt(std::cout,100.0*tschur/ttotal       ,6,2) << "%)\n";
			std::cout << "   cond    " << testimator   << " ("; fmt(std::cout,100.0*testimator/ttotal  ,6,2) << "%)\n";
			std::cout << "nnz A " << nzA << " LU " << nzLU << " LU2 " << nzLU2tot << " swaps " << totswaps << "levels " << level_size.size() << "\n";
		}
		if (verbosity > 1) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
		init = true;
		/*
		{
			std::fstream fout("Diag.mtx",std::ios::out);
			fout << "%%MatrixMarket matrix coordinate real general" << std::endl;
			fout << std::scientific;
			fout << moend-mobeg << " " << moend-mobeg << " " << moend-mobeg << std::endl;
			for(k = mobeg; k < moend; ++k)
				fout << k+1 << " " << k+1 << " " << LU_Diag[k] << std::endl;
			fout.close();
		}
		*/
		// info->PrepareVector(div);
		// std::fill(div.Begin(),div.End(),0);
		// for(k = mobeg; k < moend; k++) div[k] = 1.0;
		// info->Accumulate(div);
		// for(k = mobeg; k < moend; k++) div[k] = 1.0/div[k];

		condestL = NuL;
		condestU = NuU;

		//prepare temporal array
		temp.set_interval_beg(vbeg);
		temp.set_interval_end(vend);
		for (k = vbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(vend); ++k) temp[k] = 0.0;
		return true;
	}
	bool MLMTILUC_preconditioner::Finalize()
	{
		init = false;
		L_Address.clear();
		U_Address.clear();
		L_Entries.clear();
		U_Entries.clear();
		temp.clear();
		ddP.clear();
		ddQ.clear();
		LU_Diag.clear();
		for(int k = 0; k < (int)E_Address.size(); ++k)
			delete E_Address[k];
		for(int k = 0; k < (int)F_Address.size(); ++k)
			delete F_Address[k];
		E_Address.clear();
		F_Address.clear();
		level_size.clear();
		level_interval.clear();
		E_Entries.clear();
		F_Entries.clear();
		return true;
	}
	
	int MLMTILUC_preconditioner::Descend(INMOST_DATA_ENUM_TYPE level, Sparse::Vector & inout)
	{
		INMOST_DATA_ENUM_TYPE k, m, cbeg, cend;
		//current B block interval
		cbeg = level_interval[level].first;
		cend = level_interval[level].last;
		//check that it makes sense to make restriction step
		if( cend < level_interval.back().last )
		{
			// step 1) obtain ~f
			//apply B^{-1} on f
			for (k = cbeg; k < cend; ++k) temp[k] = inout[k];
			//Solve with L first
			for (k = cbeg; k < cend; ++k) if( L_Address[k].Size() )//iterator over columns of L
			{
				for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r < L_Address[k].last; ++r)
					temp[L_Entries[r].first] -= temp[k] * L_Entries[r].second;
			}
			//Solve with diagonal
			for (k = cbeg; k < cend; ++k) temp[k] /= LU_Diag[k];
			//Solve with U
			for (k = cend; k > cbeg; --k) if( U_Address[k - 1].Size() ) //iterator over rows of U
			{
				for (INMOST_DATA_ENUM_TYPE r = U_Address[k - 1].first; r < U_Address[k - 1].last; ++r)
					temp[k - 1] -= temp[U_Entries[r].first] * U_Entries[r].second;
			}
			//now can calculate ~g
			//multiply ~f by row of E
			for (k = cend; k < level_interval.back().last; ++k) if( E_Address[level]->at(k).Size() )
			{
				for (m = E_Address[level]->at(k).first; m < E_Address[level]->at(k).last; ++m)
				{
					assert(E_Entries[m].first < cend);
					inout[k] -= temp[E_Entries[m].first] * E_Entries[m].second;
				}
			}
		}
		return level + 1;
	}
	int MLMTILUC_preconditioner::Ascend(INMOST_DATA_ENUM_TYPE level, Sparse::Vector & inout)
	{
		level = level - 1;
		INMOST_DATA_ENUM_TYPE cbeg,cend, k, m;
		cbeg = level_interval[level].first;
		cend = level_interval[level].last;
		//calculate Fy, iterate over ~g, write to unused input vector
		//for (k = cbeg; k < cend; ++k)
		//{
		//	for (m = F_Address[k].first; m < F_Address[k].last; m++)
		//		output[k] -= output[F_Entries[m].first] * F_Entries[m].second;
		//}
		//check that it makes sense to make prolongation step
		for (k = cend; k < level_interval.back().last; ++k) if( F_Address[level]->at(k).Size() )
		{
				for (m = F_Address[level]->at(k).first; m < F_Address[level]->at(k).last; ++m)
				{
					assert(F_Entries[m].first < cend);
					inout[F_Entries[m].first] -= inout[k] * F_Entries[m].second;
				}
		}
		//perform solve over calculated vector
		//Solve with L first
		for (k = cbeg; k < cend; ++k) if( L_Address[k].Size() )//iterator over columns of L
		{
			for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r < L_Address[k].last; ++r)
				inout[L_Entries[r].first] -= inout[k] * L_Entries[r].second; //r->first alwayse > k
		}
		//Solve with diagonal
		for (k = cbeg; k < cend; ++k) inout[k] /= LU_Diag[k];
		//Solve with U
		for (k = cend; k > cbeg; --k) if( U_Address[k - 1].Size() )//iterator over rows of U
		{
			for (INMOST_DATA_ENUM_TYPE r = U_Address[k - 1].first; r < U_Address[k - 1].last; ++r)
				inout[k - 1] -= inout[U_Entries[r].first] * U_Entries[r].second; // r->first always > k
		}
		return level;
	}
	bool MLMTILUC_preconditioner::Solve(Sparse::Vector & input, Sparse::Vector & output)
	{
		assert(&input != &output);
		//
#if defined(USE_OMP)
#pragma omp single
#endif
		{
			INMOST_DATA_ENUM_TYPE k, mobeg, moend, vbeg, vend;
			info->GetOverlapRegion(info->GetRank(), mobeg, moend);
			info->GetVectorRegion(vbeg, vend);
		
			
			for (k = mobeg; k < moend; k++) temp[k] = input[k];
			
			//the original matrix A was separated by multilevel algorithm to the following form
			//     | B  F |
			// A = |      |
			//     | E  C |
			// In order to apply solve phase on this matrix, we consider the matrix in following sense:
			//     | I         0 | | B   0 | | I    B^{-1} F |
			// A = |             | |       | |               | = L D U
			//     | E B^{-1}  I | | 0   S | | 0           I |
			// where S = C - E B^{-1} F
			// consider the system
			//  | B  F | | u |   | f |
			//  |      | |   | = |   |
			//  | E  C | | y |   | g |
			// then solution is obtained by steps:
			// 1) ~f = B^{-1} f
			// 2) ~g = g - E ~f
			// 3) y = S^{-1} ~g
			// 4) u = ~f - B^{-1} F y
			
			for (k = vbeg; k < mobeg; k++) output[k] = 0;
			for (k = mobeg; k < moend; ++k) output[k] = temp[ddP[k]];//*DL0[k];//*DL0[k];
			for (k = moend; k < vend; k++) output[k] = 0;
			
			INMOST_DATA_ENUM_TYPE level = 0;
			while (level < level_size.size()) level = Descend(level, output);
			while (level > 0) level = Ascend(level, output);
			
			for (k = mobeg; k < moend; ++k) temp[ddQ[k]] = output[k];//*DR0[k];//*DR0[k];
			for (k = mobeg; k < moend; ++k) output[k] = temp[k];
			
			

			//Restrict additive schwartz (maybe do it outside?)
			//May assamble partition of unity instead of restriction before accumulation
			//assembly should be done instead of initialization
			// for (k = vbeg; k < mobeg; ++k) output[k] = 0;
			// for (k = moend; k < vend; ++k) output[k] = 0;
			for (k = vbeg; k < mobeg; ++k) output[k] = 0;
			for (k = moend; k < vend; ++k) output[k] = 0;
			// for (k = mobeg; k < moend; ++k) output[k] *= div[k];
		}
		info->Accumulate(output);
		return true;
	}
	bool MLMTILUC_preconditioner::ReplaceMAT(Sparse::Matrix & A){ if (isInitialized()) Finalize(); Alink = &A; return true; }
	bool MLMTILUC_preconditioner::ReplaceSOL(Sparse::Vector & x) {(void)x; return true; }
	bool MLMTILUC_preconditioner::ReplaceRHS(Sparse::Vector & b) {(void)b; return true; }
	Method * MLMTILUC_preconditioner::Duplicate() { return new MLMTILUC_preconditioner(*this); }
	MLMTILUC_preconditioner::~MLMTILUC_preconditioner()
	{
		if (!isFinalized()) Finalize();
	}
#endif

