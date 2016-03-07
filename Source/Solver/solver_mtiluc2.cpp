#define _CRT_SECURE_NO_WARNINGS
#include <iomanip>
#include "inmost_solver.h"
#include "solver_mtiluc2.hpp"
#include <sstream>
//#define REPORT_ILU
//#undef REPORT_ILU
//#define REPORT_ILU_PROGRESS
//#define REPORT_ILU_END
//#define REPORT_ILU_SUMMARY
//#undef REPORT_ILU_PROGRESS

//#define USE_OMP

using namespace INMOST;

#ifndef DEFAULT_TAU
#define DEFAULT_TAU 0.01
#endif


//#define REORDER_RCM
//#define REORDER_NNZ
#if defined(USE_SOLVER_METIS)
#define REORDER_METIS_ND
#endif
#if defined(USE_SOLVER_MONDRIAAN)
//#define REORDER_MONDRIAAN
#endif
//#define REORDER_DDPQ
//#define REORDER_ZOLTAN_HUND

#define ESTIMATOR
#define ESTIMATOR_REFINE
#define RESCALE_B

#define PREMATURE_DROPPING

#define EQUALIZE_1NORM
//#define EQUALIZE_2NORM
//#define EQUALIZE_IDOMINANCE

#define PIVOT_THRESHOLD
#define PIVOT_THRESHOLD_VALUE 1.0e-9
//#define DIAGONAL_PERTURBATION
#define DIAGONAL_PERTURBATION_REL 1.0e-10
#define DIAGONAL_PERTURBATION_ABS 1.0e-12
//#define DIAGONAL_PIVOT //probably there is some bug
#define DIAGONAL_PIVOT_TAU 0.01
//#define DIAGONAL_PIVOT_COND 1.0e+6
#define ILUC2
#define TRACK_DIAGONAL

#if defined(DIAGONAL_PIVOT) && defined(DIAGONAL_PIVOT_TAU) && !defined(TRACK_DIAGONAL)
#define TRACK_DIAGONAL
#endif


#if defined(REORDER_METIS_ND)
#include "metis.h"
#endif
#if defined(REORDER_ZOLTAN_HUND)
#include <zoltan.h>
#endif
#if defined(REORDER_MONDRIAAN)
#include <Mondriaan.h>
#endif

template<class T> struct make_integer;
template<> struct make_integer<float> {typedef int type;};
template<> struct make_integer<double> {typedef long long type;};

__INLINE static bool compare(INMOST_DATA_REAL_TYPE * a, INMOST_DATA_REAL_TYPE * b)
{
	return (*reinterpret_cast< make_integer<INMOST_DATA_REAL_TYPE>::type * >(a)) <=
		(*reinterpret_cast< make_integer<INMOST_DATA_REAL_TYPE>::type * >(b));
}


class BinaryHeap
{
	
	INMOST_DATA_REAL_TYPE * Base;
	std::vector<INMOST_DATA_REAL_TYPE *> Array;
	std::vector<INMOST_DATA_ENUM_TYPE> Position;
public:
	void Clear()
	{
		while(!Array.empty())
		{
			Position[Array.back()-Base] = ENUMUNDEF;
			Array.pop_back();
		}
	}
	INMOST_DATA_REAL_TYPE * Get(INMOST_DATA_ENUM_TYPE pos) {return Array[pos];}
	INMOST_DATA_ENUM_TYPE GetSize() {return static_cast<INMOST_DATA_ENUM_TYPE>(Array.size());}
	INMOST_DATA_ENUM_TYPE GetPosition(INMOST_DATA_ENUM_TYPE pos)
	{
		return Position[pos];
	}
	INMOST_DATA_ENUM_TYPE DecreaseKey(INMOST_DATA_ENUM_TYPE pos)
	{
		INMOST_DATA_ENUM_TYPE i = Position[pos];
		++i;
		while(i > 1)
		{
			//if((*Array[i-1]) <= (*Array[i/2-1]))
			if( compare(Array[i-1],Array[i/2-1]) )
			{
				Position[(Array[i/2-1]-Base)] = i-1; 
				Position[(Array[i-1]-Base)] = i/2-1; 
				std::swap(Array[i/2-1],Array[i-1]);
			}
			else break;
			i = i/2;
		}
		return i;
	}
	INMOST_DATA_ENUM_TYPE PushHeap(INMOST_DATA_REAL_TYPE * key)
	{
		INMOST_DATA_ENUM_TYPE i = GetSize();
		Array.push_back(key);
		Position[(key-Base)] = i;
		++i;
		while(i > 1)
		{
			//if((*Array[i-1]) <= (*Array[i/2-1]) )
			if( compare(Array[i-1],Array[i/2-1]) )
			{
				Position[(Array[i-1]-Base)] = i/2-1; 
				Position[(Array[i/2-1]-Base)] = i-1; 
				std::swap(Array[i-1],Array[i/2-1]);
			}
			else break;
			i = i/2;
		}
		return i;
	}

	void BalanceHeap(INMOST_DATA_ENUM_TYPE i)
	{
		INMOST_DATA_ENUM_TYPE Index;
		++i;
		while(i <= Array.size()/2)
		{
			if( 2*i+1 > Array.size() )
				Index = 2*i;
			//else if( (*Array[2*i-1]) <= (*Array[2*i+1-1]) )
			else if( compare(Array[2*i-1],Array[2*i+1-1]) )
				Index = 2*i;
			else
				Index = 2*i+1;
			//if(!((*Array[i-1]) <= (*Array[Index-1])))
			if(!compare(Array[i-1],Array[Index-1]))
			{ 
				Position[(Array[i-1]-Base)] = Index-1;
				Position[(Array[Index-1]-Base)] = i-1; 
				std::swap(Array[i-1],Array[Index-1]);
			}
			else break;
			i = Index;
		}
	}
	INMOST_DATA_ENUM_TYPE PopHeap()
	{
		INMOST_DATA_ENUM_TYPE Ret = ENUMUNDEF;
		if(Array.empty()) return Ret;
		Ret = static_cast<INMOST_DATA_ENUM_TYPE>(Array[0]-Base);
		Array[0] = Array.back();
		Array.pop_back();
		Position[Array[0]-Base] = 0;
		Position[Ret] = ENUMUNDEF;
		BalanceHeap(0);	
		return Ret;
	}
	BinaryHeap(INMOST_DATA_REAL_TYPE * Base, INMOST_DATA_ENUM_TYPE Size) 
		: Base(Base)
	{
		Position.resize(Size,ENUMUNDEF);
		Array.reserve(4096);
	}
	~BinaryHeap()
	{
	}
};



	void MTILUC_preconditioner::DumpMatrix(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
									std::vector<Sparse::Row::entry> & Entries,
									INMOST_DATA_ENUM_TYPE wmbeg, INMOST_DATA_ENUM_TYPE wmend,
									std::string file_name)
	{
		INMOST_DATA_REAL_TYPE norm1 = 0, norm2 = 0, max = -1.0e54, min = 1.0e54, minabs = 1.0e54, vrowmax, diag, mindiag = 1.0e54, maxdiag = -1.0e54, maxabsdiag = -1.0e54, minabsdiag = 1.0e54;
		INMOST_DATA_ENUM_TYPE nnz = 0, dominant_rows = 0, dominant_cols = 0, irowmax = 0, nodiag = 0;
		interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> vcolmax(wmbeg,wmend,0);
		interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> icolmax(wmbeg,wmend,ENUMUNDEF);
		for (INMOST_DATA_ENUM_TYPE k = wmbeg; k < wmend; ++k) 
		{
			nnz += Address[k].Size();
			vrowmax = 0;

			bool diag_found = false;
			diag = 0;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it != Address[k].last; ++it)
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
		
		fout << wmend-wmbeg << " " << wmend-wmbeg << " " << nnz << std::endl;;
		for (INMOST_DATA_ENUM_TYPE k = wmbeg; k < wmend; ++k)
		{
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it != Address[k].last; ++it)
				fout << (k-wmbeg+1) << " " << (Entries[it].first-wmbeg+1) << " " << Entries[it].second << std::endl;
		}
		fout.close();
	}
	void MTILUC_preconditioner::SwapEntries(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
					 std::vector<Sparse::Row::entry> & Entries, 
					 INMOST_DATA_ENUM_TYPE rbeg, INMOST_DATA_ENUM_TYPE rend, 
					 INMOST_DATA_ENUM_TYPE k, INMOST_DATA_ENUM_TYPE j)
	{
		INMOST_DATA_ENUM_TYPE i;
		INMOST_DATA_REAL_TYPE u;
		for (i = rbeg; i < rend; i++)
		{
			INMOST_DATA_ENUM_TYPE pk = Address[i].last, pj = Address[i].last, r;
			for (r = Address[i].first; r < Address[i].last; r++)
			{
				if (Entries[r].first >= k) { pk = r; break; }
			}
			for (r = pk; r < Address[i].last; r++)
			{
				if (Entries[r].first >= j) { pj = r; break; }
			}
			if (pk < Address[i].last)
			{
				if (pj < Address[i].last && Entries[pj].first == j)
				{
					if (Entries[pk].first == k)
					{
						u = Entries[pj].second;
						Entries[pj].second = Entries[pk].second;
						Entries[pk].second = u;
					}
					else
					{
						u = Entries[pj].second;
						for (r = pj; r > pk; r--) Entries[r] = Entries[r - 1];
						Entries[pk].first = k;
						Entries[pk].second = u;
					}
				}
				else
				{
					if (Entries[pk].first == k)
					{
						u = Entries[pk].second;
						for (r = pk; r < pj - 1; r++) Entries[r] = Entries[r + 1];
						Entries[pj - 1].first = j;
						Entries[pj - 1].second = u;
					}
				}
			}
		}
	}
	void MTILUC_preconditioner::SwapLine(interval<INMOST_DATA_ENUM_TYPE, Interval> & Line, INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j)
	{
		INMOST_DATA_ENUM_TYPE temp;
		temp = Line[i].first;
		Line[i].first = Line[j].first;
		Line[j].first = temp;
		temp = Line[i].last;
		Line[i].last = Line[j].last;
		Line[j].last = temp;
	}
	
	void MTILUC_preconditioner::inversePQ(INMOST_DATA_ENUM_TYPE wbeg,
				   INMOST_DATA_ENUM_TYPE wend, 
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ)
	{
		//inverse reordering
		// in invPQ numbers indicate where to get current column
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k != wend; ++k) invP[localP[k]] = k;
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k != wend; ++k) invQ[localQ[k]] = k;
	}
	void MTILUC_preconditioner::applyPQ(INMOST_DATA_ENUM_TYPE wbeg,
				 INMOST_DATA_ENUM_TYPE wend,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ)
	{
		INMOST_DATA_ENUM_TYPE k;
		// compute reordering in global P,Q, we need it to compute reordering in vector during solve phase
		for (k = wbeg; k != wend; ++k)
		{
			localP[k] = ddP[invP[k]];
			localQ[k] = ddQ[invQ[k]];
		}
		// update reordering in global P,Q
		for (k = wbeg; k != wend; ++k)
		{
			ddP[k] = localP[k];
			ddQ[k] = localQ[k];
		}
	}
	INMOST_DATA_ENUM_TYPE & MTILUC_preconditioner::EnumParameter(std::string name)
	{
		if (name == "scale_iters") return sciters;
		else if( name == "estimator" ) return estimator;
		throw - 1;
	}
	INMOST_DATA_REAL_TYPE & MTILUC_preconditioner::RealParameter(std::string name)
	{
		if (name == "tau") return tau;
		else if( name == "tau2" ) return iluc2_tau;
    else if( name == "condition_number_L" ) return condestL;
    else if( name == "condition_number_U" ) return condestU;
		throw - 1;
	}
	void MTILUC_preconditioner::Copy(const Method * other)
	{
		const MTILUC_preconditioner * b = dynamic_cast<const MTILUC_preconditioner *>(other);
		assert(b != NULL);
		tau = b->tau;
		Alink = b->Alink;
		info = b->info;
		sciters = b->sciters;
		eps = b->eps;
	}
	MTILUC_preconditioner::MTILUC_preconditioner(const MTILUC_preconditioner & other) :Method(other)
	{
		Copy(&other);
	}
	MTILUC_preconditioner & MTILUC_preconditioner::operator =(MTILUC_preconditioner const & other)
	{
		Copy(&other);
		return *this;
	}
	MTILUC_preconditioner::MTILUC_preconditioner(Solver::OrderInfo & info)
		:tau(DEFAULT_TAU),info(&info)
	{
		Alink = NULL;
		init = false;
		sciters = 8;
		eps = 1e-54;
#if defined(ESTIMATOR)
		estimator = 1;
#else
		estimator = 0;
#endif
		tau = 1.0e-3;
		iluc2_tau = tau*tau;
	}
	bool MTILUC_preconditioner::isInitialized() { return init; }
	bool MTILUC_preconditioner::isFinalized() { return !init; }
	bool MTILUC_preconditioner::Initialize()
	{

		const INMOST_DATA_REAL_TYPE subst = 1.0;
		const INMOST_DATA_REAL_TYPE tol_modif = PIVOT_THRESHOLD_VALUE;
		const INMOST_DATA_ENUM_TYPE UNDEF = ENUMUNDEF, EOL = ENUMUNDEF - 1;

		
		INMOST_DATA_ENUM_TYPE wbeg, wend; //working interval
		INMOST_DATA_ENUM_TYPE mobeg, moend; // total interval
		INMOST_DATA_ENUM_TYPE vbeg, vend; // vector interval
		
		INMOST_DATA_ENUM_TYPE k, i, j, Li, Ui, curr, next;
		INMOST_DATA_REAL_TYPE l,u,udiag, abs_udiag, max_diag, min_diag, mean_diag;
		INMOST_DATA_ENUM_TYPE nzA, nzLU = 0;
		Sparse::Vector DL, DR;
		info->GetOverlapRegion(info->GetRank(), mobeg, moend);
		info->GetVectorRegion(vbeg,vend);
		
		//prepare temporal array
		temp.set_interval_beg(vbeg);
		temp.set_interval_end(vend);
		//for(interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE>::iterator rit = temp.begin();
		//	rit != temp.end(); ++rit) *rit = 0.0;
		std::fill(temp.begin(), temp.end(), 0.0);

		
		//prepare reordering vectors
		ddP.set_interval_beg(mobeg);
		ddP.set_interval_end(moend);
		ddQ.set_interval_beg(mobeg);
		ddQ.set_interval_end(moend);

		

		//prepare rescaling vectors
		DL.SetInterval(mobeg, moend);
		DR.SetInterval(mobeg, moend);
		for(Sparse::Vector::iterator ri = DL.Begin(); ri != DL.End(); ++ri) *ri = 0.0;

		

		for (k = mobeg; k != moend; k++) ddP[k] = ddQ[k] = k;
		// supplementary data for ddPQ reordering
		interval<INMOST_DATA_ENUM_TYPE, bool> donePQ(mobeg, moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> invP(mobeg, moend), invQ(mobeg,moend), localP(mobeg, moend,ENUMUNDEF), localQ(mobeg,moend,ENUMUNDEF);		

		//supplimentary data structures for ILUC
		INMOST_DATA_ENUM_TYPE LU_Beg;
		U_Address.set_interval_beg(mobeg);
		U_Address.set_interval_end(moend);
		L_Address.set_interval_beg(mobeg);
		L_Address.set_interval_end(moend);
		B_Address.set_interval_beg(mobeg);
		B_Address.set_interval_end(moend);
		LU_Diag.set_interval_beg(mobeg);
		LU_Diag.set_interval_end(moend);
		std::vector<Sparse::Row::entry> A_Entries;
		//std::vector<INMOST_DATA_ENUM_TYPE> sort_indeces;
		interval<INMOST_DATA_ENUM_TYPE, Interval> A_Address(mobeg, moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ulist(mobeg, moend), Llist(mobeg, moend), Blist(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ubeg(mobeg, moend,EOL), Lbeg(mobeg, moend,EOL), Bbeg(mobeg,moend,EOL);

		INMOST_DATA_REAL_TYPE NuU = 1, NuL = 1, NuD, NuU_max = 1.0, NuL_max = 1.0;
#if defined(ESTIMATOR)
		//supplimentary data structures for condition estimates of L^{-1}, U^{-1}
		INMOST_DATA_REAL_TYPE mup, mum, smup, smum, NuL1 = 1, NuL2 = 1, NuU1 = 1, NuU2 = 1;
		INMOST_DATA_REAL_TYPE NuU1_old = 1, NuL1_old = 1, NuU2_old = 1, NuL2_old = 1;
		INMOST_DATA_REAL_TYPE NuU1_new = 1, NuL1_new = 1, vp, vm, v;
#if defined(ESTIMATOR_REFINE)
		INMOST_DATA_ENUM_TYPE np, nm;
		INMOST_DATA_REAL_TYPE NuU2_new, NuL2_new;
#endif
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL1(mobeg, moend,0.0), EstU1(mobeg, moend,0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL2(mobeg, moend,0.0), EstU2(mobeg, moend,0.0);
#endif
#if defined(ESTIMATOR) && defined(DIAGONAL_PIVOT_COND)
		INMOST_DATA_REAL_TYPE NuU_tmp, NuL_tmp;
#endif
		//supplimentary data structures for returning values of dropped elements
		//INMOST_DATA_REAL_TYPE DropLk, DropUk;
		//interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> DropU(mobeg,moend,0.0), DropL(mobeg,moend,0.0);
		//data structure for linked list

		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> U(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> V(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Dist(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());		
		
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LineValuesU(mobeg, moend,0.0), LineValuesL(mobeg,moend,0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> LineIndecesU(mobeg, moend+1,UNDEF), LineIndecesL(mobeg,moend+1,UNDEF);
		double tfactor = 0.0, trescale = 0.0, treorder = 0.0, treassamble = 0.0, ttotal, tt;
#if defined(REORDER_METIS_ND)
		double tmetisgraph = 0, tmetisnd = 0;
#endif
        double tlfactor, tlrescale, tlreorder, tlreassamble;
		ttotal = Timer();

		
		//calculate number of nonzeros
		nzA = 0;
		for (k = mobeg; k < moend; ++k)
		{
			for (Sparse::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
				if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0) nzA++;
		}

		//sort_indeces.reserve(256);
		A_Entries.resize(nzA);
		B_Entries.reserve(nzA);
		LU_Entries.reserve(nzA*4);

#if defined(REORDER_NNZ)
		wgt_coords sort_wgts(2*(moend - mobeg));
#endif

#if defined(REORDER_ZOLTAN_HUND)
#endif

#if defined(REORDER_MONDRIAAN)
		{
			sparsematrix sp;
			int numparts = 8;
			int * cols = new int[nzA];
			int * rows = new int[Alink->Size()+1];
			double * values = new double[nzA];
			int cnt = 0;
			rows[0] = 0;
			for (k = mobeg; k < moend; ++k)
			{
				for (Sparse::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
				{
					if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0)
					{
						cols[cnt] = r->first;
						values[cnt] = r->second;
						++cnt;
					}
				}
				rows[k-mobeg+1] = cnt;
			}
			CRSSparseMatrixInit(&sp,Alink->Size(),Alink->Size(),nzA,cols,rows,values,0);
			
			PstartInit(&sp,numparts);
			opts options;
			
			SetDefaultOptions(&options);
			//SetOption(&options,"SplitMethod","KLFM");
			SetOption(&options,"Permute","BBD");
			SetOption(&options,"SplitStrategy","onedimrow");
			DistributeMatrixMondriaan(&sp,numparts,0.1,&options,NULL);
			sp.

			j = 0;
			for (k = mobeg; k < moend; ++k)
			{
				A_Address[k].first = j;
				for (Sparse::Row::iterator r = (*Alink)[sp.row_perm[k-mobeg]].Begin(); r != (*Alink)[sp.row_perm[k-mobeg]].End(); ++r)
				{
					if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0)
						A_Entries[j++] = Sparse::Row::make_entry(sp.col_perm_inv[r->first-mobeg], r->second);
				}
				A_Address[k].last = j;
				assert(A_Address[k].Size() != 0); //singular matrix
			}
			DumpMatrix(A_Address,A_Entries,mobeg,moend,"mondriaan.mtx");
			/*
			sp.MMTypeCode[0] = 'M';
			FILE * f = fopen("mondriaan.mtx","w");
			MMWriteSparseMatrix(&sp,f,"mondriaan",&options);
			fclose(f);
			*/
			MMDeleteSparseMatrix(&sp);
			delete [] cols;
			delete [] rows;
			delete [] values;
		}
#else
		
		j = 0;
		for (k = mobeg; k < moend; ++k)
		{
			A_Address[k].first = j;
			for (Sparse::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
			{
				if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0)
					A_Entries[j++] = Sparse::Row::make_entry(r->first, r->second);
			}
			A_Address[k].last = j;
			assert(A_Address[k].Size() != 0); //singular matrix
		}
#endif

		std::vector<INMOST_DATA_REAL_TYPE> C_Entries(A_Entries.size());

#if defined(ILUC2)
		INMOST_DATA_ENUM_TYPE nzLU2 = 0;
		INMOST_DATA_REAL_TYPE tau2 = iluc2_tau;
		std::vector<Sparse::Row::entry> LU2_Entries;
		interval<INMOST_DATA_ENUM_TYPE, Interval> L2_Address(mobeg, moend+1), U2_Address(mobeg, moend+1);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2list(mobeg, moend, UNDEF), L2list(mobeg, moend, UNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2beg(mobeg, moend, EOL), L2beg(mobeg, moend, EOL);
		LU2_Entries.reserve(nzA*4);
#else
		INMOST_DATA_REAL_TYPE tau2 = tau;
#endif

#if defined(REPORT_ILU)
		std::cout << "nonzeros in A " << nzA << std::endl;
#endif
		//set up working interval
		wbeg = mobeg;
		wend = moend;
		while (wbeg < wend) //iterate into levels until all matrix is factored
		{
			//std::fill(localP.begin() + (wbeg - mobeg), localP.begin() + (wend - mobeg), ENUMUNDEF);
			//std::fill(localQ.begin() + (wbeg - mobeg), localQ.begin() + (wend - mobeg), ENUMUNDEF);
			//ddPQ reordering on current Schur complement
			INMOST_DATA_ENUM_TYPE cbeg = wbeg, cend = wend; //next size of factored B block
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// MAXIMUM TRANSVERSE REORDERING /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			{
				tt = Timer();
				INMOST_DATA_ENUM_TYPE ColumnBegin;
				INMOST_DATA_ENUM_TYPE pop_heap_pos;
				INMOST_DATA_REAL_TYPE ShortestPath, AugmentPath;
				INMOST_DATA_ENUM_TYPE PathEnd, Trace, IPermPrev;
				//Sparse::Vector & U = DL;
				//Sparse::Vector & V = DR;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Cmax(wbeg,wend,0.0);
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & Perm = localQ;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & IPerm = localP;
				//interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & UpdateStack = Ulist;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & ColumnList = Llist;//(wbeg,wend,ENUMUNDEF);
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & Parent = Blist;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> AugmentPosition(wbeg,wend,ENUMUNDEF);
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> ColumnPosition(wbeg,wend,ENUMUNDEF);
				
				BinaryHeap Heap(&Dist[wbeg],wend-wbeg);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Arrays initialization   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				
				//std::fill(U.Begin() + wbeg - mobeg, U.Begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				//std::fill(V.Begin() + wbeg - mobeg, V.Begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				//std::fill(Cmax.begin() + wbeg - mobeg, Cmax.begin() + wend - mobeg, 0.0);
				//std::fill(Dist.begin() + wbeg - mobeg, Dist.begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				std::fill(Perm.begin() + wbeg - mobeg, Perm.begin() + wend - mobeg, ENUMUNDEF);
				std::fill(IPerm.begin() + wbeg - mobeg, IPerm.begin() + wend - mobeg, ENUMUNDEF);
				std::fill(Parent.begin() + wbeg - mobeg, Parent.begin() + wend - mobeg, ENUMUNDEF);
				std::fill(ColumnList.begin() + wbeg - mobeg, ColumnList.begin() + wend - mobeg, ENUMUNDEF);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Initial LOG transformation to dual problem and initial extreme match /////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//double T = Timer();
				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
					{
						i = A_Entries[it].first;
						u = C_Entries[it] = fabs(A_Entries[it].second);
						if( u > Cmax[i] ) Cmax[i] = u;
						//C_Entries.push_back(Sparse::Row::make_entry(i,u));
					}
				}

				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
					{
						i = A_Entries[it].first;
						if( Cmax[i] == 0.0 )
							C_Entries[it] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
						else
						{
							C_Entries[it] = log(Cmax[i])-log(C_Entries[it]);
							//C_Entries[it] = log10(Cmax[i])-log10(C_Entries[it]);
							//C_Entries[it] = fabs(log10(Cmax[i]/C_Entries[it]));
							//C_Entries[it] = fabs(log(Cmax[i]/C_Entries[it]));
							if( C_Entries[it] < U[i] ) U[i] = C_Entries[it];
						}
					}
				}
				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
					{
						u = C_Entries[it] - U[A_Entries[it].first];
						if( u < V[k] ) V[k] = u;
					}
				}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Update cost and match ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// 1-step augmentation   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				for(k = wbeg; k < wend; ++k)
				{
					if( IPerm[k] == ENUMUNDEF ) //unmatched row
					{
						for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last && IPerm[k] == ENUMUNDEF; ++it)
						{
							u = fabs(C_Entries[it] - V[k] - U[A_Entries[it].first]);
							if( u <= 1.0e-30 )
							{
								Li = Perm[A_Entries[it].first];
								assert(Li != ENUMUNDEF);
								// Search other row in C for 0
								for (INMOST_DATA_ENUM_TYPE Lit = A_Address[Li].first; Lit != A_Address[Li].last; ++Lit)
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

				
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Weighted bipartite matching ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for(k = wbeg; k < wend; ++k)
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
						for (INMOST_DATA_ENUM_TYPE Lit = A_Address[Li].first; Lit != A_Address[Li].last; ++Lit)
						{
							Ui = A_Entries[Lit].first;
							//if( ColumnList[Ui] == k ) continue;
							if( ColumnList[Ui] != ENUMUNDEF ) continue;
							l = ShortestPath + C_Entries[Lit] - V[Li] - U[Ui];
							//if( l < 0.0 ) printf("row %d col %d negative l %g Shortest %lf C %lf V %lf U %lf\n",k,Ui,l,C_Entries[Lit],V[Li],U[Ui]);
							if( l < 0.0 && fabs(l) < 1.0e-12 ) l = 0;
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
									Dist[Ui] = l;
									Parent[Perm[Ui]] = Li;
									AugmentPosition[Ui] = Lit;
									if( Heap.GetPosition(Ui-wbeg) != ENUMUNDEF )
										Heap.DecreaseKey(Ui-wbeg);
									else 
										Heap.PushHeap(&Dist[Ui]);
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
						for(Ui = 0; Ui < Heap.GetSize(); ++Ui) *Heap.Get(Ui) = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
						Heap.Clear();
					}
				}
				//printf("Maximum product transversal %lf\n",Timer()-T);
				//fclose(rec);

				//T = Timer();
				for (k = cbeg; k < cend; ++k)
				{
					//B_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());

					bool flip_sign = false;
					
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt != A_Address[k].last; ++jt)
					{
						i = A_Entries[jt].first;
						j = Perm[A_Entries[jt].first];
						
						l = exp(V[k]);
						u = exp(U[i])/Cmax[i];
						DL[k] = l;
						DR[i] = u;
						//l = l*u;
						//l = l*u;//exp(U[i]+V[k])/Cmax[i];
						//l = pow(10,U[i]+V[k])/Cmax[i];
						//l = exp(U[i]+V[k])/Cmax[i];
						//A_Entries[jt].first = j;
						//A_Entries[jt].second *= l;

						
						
						//B_Entries.push_back(Sparse::Row::make_entry(j, l*A_Entries[jt].second));

						if( j == k && (l*A_Entries[jt].second*u) < 0.0 ) flip_sign = true;
					}

					//B_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());

					if( flip_sign )
					{
						DL[k] *= -1;
						//for(Li = B_Address[k].first; Li < B_Address[k].last; ++Li)
						//	B_Entries[Li].second *= -1;
					}
						

					//std::sort(B_Entries.begin() + B_Address[k].first, B_Entries.end());
					//std::sort(A_Entries.begin() + A_Address[k].first, A_Entries.begin() + A_Address[k].last);
				}
				//printf("Reorder matrix %lf\n",Timer()-T);

				//DumpMatrix(A_Address,A_Entries,wbeg,wend,"A.mtx");
				//DumpMatrix(B_Address,B_Entries,wbeg,wend,"MC64.mtx");

				std::fill(localP.begin() + (wbeg - mobeg), localP.begin() + (wend - mobeg), ENUMUNDEF);
				
				//exit(-1);

				tt = Timer() - tt;

				treorder += tt;
			}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// END MAXIMUM TRANSVERSE REORDERING /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined(REORDER_DDPQ)

#if defined(REPORT_ILU)
			std::cout << level_size.size() << " calculate weights" << std::endl;
#endif
			std::fill(DL.Begin() + wbeg - mobeg, DL.Begin() + wend - mobeg, 0.0);
			std::fill(DR.Begin() + wbeg - mobeg, DR.Begin() + wend - mobeg, 0.0);
			std::fill(Ulist.begin() + wbeg - mobeg, Ulist.begin() + wend - mobeg, 0);
			std::fill(Llist.begin() + wbeg - mobeg, Llist.begin() + wend - mobeg, 0);
			for (k = wbeg; k < wend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
				{
					i = A_Entries[it].first;
					u = fabs(A_Entries[it].second);
					//row-col norms
					DL[k] += u; //row norm
					DR[i] += u; //col norm
					//nonzeros
					Ulist[k]++; //row nnz
					Llist[i]++; //col nnz
				}
			}

			INMOST_DATA_REAL_TYPE ddmaxall = 0, ddmaxcur, ddcur, ddmaxmean = 0.0, ddmaxmin = 1.0e+20;
			for (k = wbeg; k < wend; ++k)
			{
				ddmaxcur = 0.0;
				Blist[k] = ENUMUNDEF;
				for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
				{
					i = A_Entries[it].first;
					u = fabs(A_Entries[it].second);
					ddcur = u / (DL[k] + DR[i] - u);
					if (ddmaxcur < ddcur)
					{
						ddmaxcur = ddcur;
						Blist[k] = i;
					}
				}
				temp[k] = ddmaxcur;
				if (ddmaxall < ddmaxcur) ddmaxall = ddmaxcur;
				if (ddmaxmin > ddmaxcur) ddmaxmin = ddmaxcur;
				ddmaxmean += ddmaxcur*ddmaxcur;
			}

			ddmaxmean = sqrt(ddmaxmean / (1.0*(wend - wbeg)));

			
			//first select those elements that are largest over row and column simultaneously
			sort_wgts.clear();
			for (k = wbeg; k < wend; ++k)
			{
				if (Blist[k] != ENUMUNDEF)
				{
					//if ((temp[k]-ddmaxmin)/(ddmaxall-ddmaxmin) > ddpq_tau)
					//	sort_wgts.push_back(wgt_coord((Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
					//sort_wgts.push_back(wgt_coord(temp[k], coord(k, Blist[k])));
					//sort_wgts.push_back(wgt_coord(pow(temp[k],0.25)/(Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
					//sort_wgts.push_back(wgt_coord(1.0/(Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
					sort_wgts.push_back(wgt_coord(temp[k]/(Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
				}
			}
			//if(reorder_nnz) 
				std::sort(sort_wgts.rbegin(), sort_wgts.rend());
#if defined(REPORT_ILU)
			std::cout << level_size.size() << " fill reordering" << std::endl;
#endif
			i = wbeg;
			for (wgt_coords::iterator it = sort_wgts.begin(); it != sort_wgts.end(); ++it)
			{
				if (localP[it->second.first] == ENUMUNDEF && localQ[it->second.second] == ENUMUNDEF)
				{
					localP[it->second.first] = i;
					localQ[it->second.second] = i;
					i++;
				}

				if( it->first / sort_wgts.begin()->first < ddpq_tau ) break;
				if( it-sort_wgts.begin() > 0.01*(moend-mobeg) ) break; //5% of all values
			}
			cend = i;
#elif defined(REORDER_METIS_ND)
			tt = Timer();
			idx_t nvtxs = wend-wbeg;
			std::vector<idx_t> xadj(nvtxs+1), adjncy, perm(nvtxs),iperm(nvtxs);
			adjncy.reserve(nzA*2);
			
			for (k = cbeg; k < cend; ++k) if( localQ[k] == ENUMUNDEF ) printf("%s:%d No column permutation for row %d. Matrix is structurally singular\n",__FILE__,__LINE__,k);

			for (k = cbeg; k < cend; ++k) localP[k] = k;
			for (k = cbeg; k < cend; ++k) V[localQ[k]] = DR[k]; //if you expirience a bug here then the matrix is structurally singular
			for (k = cbeg; k < cend; ++k) DR[k] = V[k];

			
			
			for(k = wbeg; k < wend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt != A_Address[k].last; ++jt)
				{
					A_Entries[jt].first = localQ[A_Entries[jt].first];
					//A_Entries[jt].second *= DL[k]*DR[A_Entries[jt].first];
					
				}
				std::sort(A_Entries.begin()+A_Address[k].first,A_Entries.begin()+A_Address[k].last);
			}

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
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[i].first; jt != A_Address[i].last; ++jt)
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

			std::fill(Bbeg.begin(),Bbeg.end(),EOL);

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
			//printf("before metis\n");
			//METIS_NodeNDP(nvtxs,&xadj[0],&adjncy[0],NULL,4,options,&perm[0],&iperm[0],&sizes[0]);
			tmetisnd = Timer();
			METIS_NodeND(&nvtxs,&xadj[0],&adjncy[0],NULL,options,&perm[0],&iperm[0]);
			tmetisnd = Timer()-tmetisnd;
			//printf("after metis\n");
			for(k = wbeg; k < wend; ++k)
			{
				//localP[k] = iperm[k-wbeg]+wbeg;
				//localQ[k] = iperm[k-wbeg]+wbeg;
				localP[perm[k-wbeg]+wbeg] = k;
				localQ[perm[k-wbeg]+wbeg] = k;
			}
			cend = wend;
			i = wend;

			treorder += Timer() - tt;
#elif defined(REORDER_NNZ)
			//for (k = cbeg; k < cend; ++k) V[k] = DR[localQ[k]];
			//for (k = cbeg; k < cend; ++k) DR[k] = V[k];

			//DumpMatrix(A_Address,A_Entries,cbeg,cend,"mat_original.mtx");
			std::fill(U.begin() + wbeg - mobeg, U.begin() + wend - mobeg, 0.0);
			std::fill(V.begin() + wbeg - mobeg, V.begin() + wend - mobeg, 0.0);
			

			std::fill(Ulist.begin() + wbeg - mobeg, Ulist.begin() + wend - mobeg, 0);
			std::fill(Llist.begin() + wbeg - mobeg, Llist.begin() + wend - mobeg, 0);

			for (k = wbeg; k < wend; ++k) Blist[localQ[k]] = k;

			for (k = wbeg; k < wend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
				{
					i = A_Entries[it].first;
					u = fabs(A_Entries[it].second);
					//row-col norms
					U[k] += u; //row norm
					V[i] += u; //col norm
					//nonzeros
					Ulist[k]++; //row nnz
					Llist[i]++; //col nnz

					if( Blist[i] == k ) temp[k] = u;
				}
			}

			sort_wgts.clear();
			for (k = wbeg; k < wend; ++k)
                sort_wgts.push_back(wgt_coord(Ulist[k]+Llist[Blist[k]]-1, coord(k, Blist[k])));
				//sort_wgts.push_back(wgt_coord((Ulist[k]+Llist[Blist[k]]-1)*(U[k]+V[Blist[k]]), coord(k, Blist[k])));
				//sort_wgts.push_back(wgt_coord((Ulist[k]+Llist[Blist[k]]-1)*(U[k]+V[Blist[k]]-temp[k])/temp[k], coord(k, Blist[k])));
				//sort_wgts.push_back(wgt_coord(1.0/temp[k], coord(k, Blist[k])));
				//sort_wgts.push_back(wgt_coord((U[k]+V[Blist[k]])*(U[k]+V[Blist[k]])/(temp[k]), coord(k, Blist[k])));
			std::sort(sort_wgts.begin(), sort_wgts.end());

			i = wbeg;
			for (wgt_coords::iterator it = sort_wgts.begin(); it != sort_wgts.end(); ++it)
			{
				localP[it->second.first] = i;
				localQ[it->second.second] = i;
				++i;
			}
			cend = i;
#elif defined(REORDER_RCM)
            //compute order of each entry
            std::fill(Ulist.begin() + wbeg - mobeg, Ulist.begin() + wend - mobeg, 0);
            std::fill(Llist.begin() + wbeg - mobeg, Llist.begin() + wend - mobeg, 0);
            for (k = wbeg; k < wend; ++k)
            {
                for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
                {
                    i = A_Entries[it].first;
                    //nonzeros
                    Ulist[k]++; //row nnz
                    Llist[i]++; //col nnz
                }
            }
            //find node with the lowest order
            INMOST_DATA_ENUM_TYPE start = wbeg;
            INMOST_DATA_ENUM_TYPE index = wbeg;
            for(k = wbeg; k < wend; ++k)
                if( Ulist[k] + Llist[Blist[k]] < Ulist[start] + Llist[Blist[start]] )
                    start = k;
            localP[start] = localQ[Blist[start]] = index++;
            
#else
			cend = wend;
			i = wbeg;
#endif
			if (cbeg == cend && cbeg != wend)
			{
				//std::cout << __FILE__ << ":" << __LINE__ << " singular matrix, factored " << mobeg << ".." << cend << " out of " << mobeg << ".." << moend << std::endl;
				for (k = cbeg; k < wend; ++k)
					LU_Diag[k] = tol_modif;
				break;
			}
#if defined(REORDER_DDPQ)
#if defined(REPORT_ILU)
			//std::cout << level_size.size() << " new level " << cbeg << ".." << cend << " out of " << wbeg << ".." << wend << std::endl;
			std::cout << std::setw(8) << level_size.size() << " total " << std::setw(8) << wend - wbeg << " after tau filtering " << std::setw(8) << sort_wgts.size() << " selected " << std::setw(8) << cend - wbeg;
			std::cout << std::scientific << " max " << std::setw(8) << ddmaxall << " mean " << std::setw(8) << ddmaxmean << " min " << std::setw(8) << ddmaxmin << " ratio " << ddmaxall/ddmaxmean;
			std::cout << std::endl;
#endif
#endif
			tt = Timer();
			//finish reordering
			if (i < wend)
			{
				j = i;
				for (k = wbeg; k < wend; ++k)
				{
					if (localP[k] == ENUMUNDEF) localP[k] = i++;
					if (localQ[k] == ENUMUNDEF) localQ[k] = j++;
				}
			}

			for (k = cbeg; k < cend; ++k) 
			{
				U[localP[k]] = DL[k];
				V[localQ[k]] = DR[k];

				//U[k] = DL[localP[k]];
				//V[k] = DR[localQ[k]];
			}
			for (k = cbeg; k < cend; ++k) 
			{
				DL[k] = U[k];
				DR[k] = V[k];
			}

			tlreorder = Timer() - tt;
			treorder += tlreorder;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// REASSAMBLE             ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			tt = Timer();
			double tt1, tt2,ttt;
			//in localPQ numbers indicate where to put current row/column
			//reorder E,F blocks by swaps
			tt1 = Timer();
			//inverse ordering
			inversePQ(wbeg,wend,localP,localQ, invP,invQ);
			tt1 = Timer() - tt1;

			//std::cout << "reorder: " << tt1 << std::endl;
			tt2 = Timer();
			nzA = 0;
			for (k = cbeg; k < cend; ++k)
			{
				B_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[invP[k]].first; jt != A_Address[invP[k]].last; ++jt)
				{
					i = localQ[A_Entries[jt].first];
					u = A_Entries[jt].second;
					B_Entries.push_back(Sparse::Row::make_entry(i, u));
				}
				B_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());
				std::sort(B_Entries.begin() + B_Address[k].first, B_Entries.end());
				nzA += A_Address[k].Size();
			}

			//DumpMatrix(B_Address,B_Entries,cbeg,cend,"mat_reorder.mtx");
			//Setup column addressing for B,F, in descending order to keep list ordered
			for (k = cend; k > cbeg; --k)
			{
				if (B_Address[k - 1].Size() > 0)
				{
					i = B_Entries[B_Address[k - 1].first].first;
					if (i < k-1)
					{
						Blist[k-1] = Bbeg[i];
						Bbeg[i] = k-1;
					}
				}
			}

			tt2 = Timer() - tt2;

			ttt = Timer();
			applyPQ(wbeg, wend, localP, localQ, invP, invQ);
			tt1 += Timer() - ttt;

			tlreassamble = Timer() - tt;
			treassamble += tlreassamble;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// RESCALING                 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Rescale current B block by Row-Column alternating scaling
			tt = Timer();
#if defined(RESCALE_B)
#if defined(REPORT_ILU)
			std::cout << " rescaling block B " << std::endl;
#endif
			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					B_Entries[r].second *= DL[k] * DR[B_Entries[r].first];

				U[k] = DL[k];
				V[k] = DR[k];
			}
			//DumpMatrix(B_Address,B_Entries,cbeg,cend,"MC64_MPTILUC2.mtx");
			/*
			{
				std::stringstream name;
				name << "mc64_" << info->GetRank() << ".mtx";
				DumpMatrix(B_Address,B_Entries,cbeg,cend,name.str());
			}
			*/
			/*
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// COMPUTE GIRSCHGORIN RADIUS ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			INMOST_DATA_REAL_TYPE radii = 0;
			for (k = cbeg; k < cend; k++)
			{
				INMOST_DATA_REAL_TYPE local_radii = 0;
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					if( B_Entries[r].first != k ) local_radii += fabs(B_Entries[r].second);
				if( radii < local_radii ) radii = local_radii;
			}

			printf("Gershgorin's radius after mc64: %e\n",radii);
			*/


			//DumpMatrix(B_Address,B_Entries,cbeg,cend,"mat_mc64.mtx");
#if defined(EQUALIZE_1NORM)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// ROW-COLUMN ALTERNATING SCALING FOR 1-NORM BALANCING////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, 0.0);
			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
					DL[k] += fabs(B_Entries[r].second);//*B_Entries[r].second;
			}
			for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
			{
				std::fill(DR.Begin()+cbeg-mobeg, DR.Begin()+cend-mobeg, 0.0);
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
						DR[B_Entries[r].first] += DL[k] * fabs(B_Entries[r].second);//*B_Entries[r].second;
				}
				for (k = cbeg; k < cend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
				std::fill(DL.Begin()+cbeg-mobeg, DL.Begin()+cend-mobeg, 0.0);
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
						DL[k] += DR[B_Entries[r].first] * fabs(B_Entries[r].second);//*B_Entries[r].second;
				}
				for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			}
			//for (k = cbeg; k < cend; k++) DL[k] = sqrt(DL[k]);
			//for (k = cbeg; k < cend; k++) DR[k] = sqrt(DR[k]);
			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					B_Entries[r].second *= DL[k] * DR[B_Entries[r].first];
			}
#elif defined(EQUALIZE_2NORM)
			std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, 0.0);
			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
					DL[k] += B_Entries[r].second*B_Entries[r].second;
			}
			for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
			{
				std::fill(DR.Begin()+cbeg-mobeg, DR.Begin()+cend-mobeg, 0.0);
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
						DR[B_Entries[r].first] += DL[k] * B_Entries[r].second*B_Entries[r].second;
				}
				for (k = cbeg; k < cend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
				std::fill(DL.Begin()+cbeg-mobeg, DL.Begin()+cend-mobeg, 0.0);
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
						DL[k] += DR[B_Entries[r].first] * B_Entries[r].second*B_Entries[r].second;
				}
				for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			}
			//for (k = cbeg; k < cend; k++) DL[k] = sqrt(DL[k]);
			//for (k = cbeg; k < cend; k++) DR[k] = sqrt(DR[k]);
			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					B_Entries[r].second *= DL[k] * DR[B_Entries[r].first];
			}
#elif defined(EQUALIZE_IDOMINANCE)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////THIS VERSION OF RESCALING INCREASES DIAGONAL DOMINANCE///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			for(k = cbeg; k < cend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
				{
					u = fabs(B_Entries[r].second);
					if( u > 0.0 )
						C_Entries[r] = -log(u);
					else 
						C_Entries[r] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
				}
			}

			
			std::fill(temp.begin() + cbeg - mobeg, temp.begin() + cend - mobeg, 0.0);

			for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
			{
				std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				std::fill(DR.Begin() + cbeg - mobeg, DR.Begin() + cend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());	
				for (k = cbeg; k < cend; k++) //row number
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					{
						i = B_Entries[r].first; //column number
						if( i != k ) //out of diagonal
						{
							u = C_Entries[r] + temp[k] - temp[i];
							DL[k] = std::min(DL[k],u);// update Y1
							DR[i] = std::min(DR[i],u);// update Y2
						}
					}
				}

				for (k = cbeg; k < cend; k++)
				{
					if( DR[k] != std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() && 
						  DL[k] != std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() )
						temp[k] += (DR[k]-DL[k])*0.5;
				}

			}

			for (k = cbeg; k < cend; k++)
			{
				DL[k] = exp(-temp[k]);
				DR[k] = exp(temp[k]);
			}

			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					B_Entries[r].second *= DL[k] * DR[B_Entries[r].first];
			}
#endif
			
			/*
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// COMPUTE GIRSCHGORIN RADIUS ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			radii = 0;
			for (k = cbeg; k < cend; k++)
			{
				INMOST_DATA_REAL_TYPE local_radii = 0;
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					if( B_Entries[r].first != k ) local_radii += fabs(B_Entries[r].second);
				if( radii < local_radii ) radii = local_radii;
			}

			printf("Gershgorin's radius after equilibration: %e\n",radii);
			*/

#if defined(EQUALIZE_1NORM) || defined(EQUALIZE_2NORM) || defined(EQUALIZE_IDOMINANCE)
			for (k = cbeg; k < cend; k++)
			{
				DL[k] *= U[k];
				DR[k] *= V[k];
			}
#else
			for (k = cbeg; k < cend; k++)
			{
				DL[k] = U[k];
				DR[k] = V[k];
			}
#endif

			//DumpMatrix(B_Address,B_Entries,cbeg,cend,"mat_equilibration.mtx");
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// RESCALING DONE ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//End rescale B block
			trescale += Timer() - tt;
#endif
			

			
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// FACTORIZATION BEGIN ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			{
        //int discarded = 0;
#if defined(REPORT_ILU)
        int swaps = 0;
#endif
#if defined(ILUC2)
				nzLU2 = 0;
				LU2_Entries.clear();
#endif
				mean_diag = 0.0;
///////////////////////////////////////////////////////////////////////////////////
//       setup diagonal values and memorize indices                              //
///////////////////////////////////////////////////////////////////////////////////
				//initialize diagonal values
				for (k = cbeg; k < cend; k++)
				{
          LU_Diag[k] = 0.0;
					localQ[k] = B_Address[k].first; //remember B block starts permutation for pivoting
					for (i = B_Address[k].first; i < B_Address[k].last; i++)
					{
						if (B_Entries[i].first == k)
						{
#if defined(TRACK_DIAGONAL)
							LU_Diag[k] = B_Entries[i].second;
#endif
#if defined(DIAGONAL_PIVOT)
							mean_diag += fabs(LU_Diag[k]);
#endif
							break;
						}
					}
				}
				//initialize first row and column
///////////////////////////////////////////////////////////////////////////////////
//       get diagonal value                                                      //
///////////////////////////////////////////////////////////////////////////////////
				assert(B_Entries[B_Address[cbeg].first].first == cbeg);
				if (B_Entries[B_Address[cbeg].first].first == cbeg)
					LU_Diag[cbeg] = B_Entries[B_Address[cbeg].first].second;
				else
				{
					std::cout << __LINE__ << " No diagonal value! " << cbeg << " " << B_Entries[B_Address[cbeg].first].first << std::endl;
					LU_Diag[cbeg] = 0.0;
				}
///////////////////////////////////////////////////////////////////////////////////
//       update for diagonal pivoting                                            //
///////////////////////////////////////////////////////////////////////////////////
#if defined(DIAGONAL_PIVOT)
				mean_diag -= fabs(LU_Diag[cbeg]);
#endif
///////////////////////////////////////////////////////////////////////////////////
//       modify diagonal value according to rules                                //
///////////////////////////////////////////////////////////////////////////////////
#if defined(PIVOT_THRESHOLD)
				if (fabs(LU_Diag[cbeg]) < tol_modif)
				{
					LU_Diag[cbeg] = LU_Diag[cbeg] < 0.0 ? -tol_modif : tol_modif;
				}
#endif
#if defined(DIAGONAL_PERTURBATION)
				{
					LU_Diag[cbeg] = LU_Diag[cbeg] * (1.0 + DIAGONAL_PERTURBATION_REL) + (LU_Diag[cbeg] < 0.0 ? -1.0 : 1.0)*DIAGONAL_PERTURBATION_ABS;
				}
#endif
///////////////////////////////////////////////////////////////////////////////////
//      form first line of U and second-order U                                  //
///////////////////////////////////////////////////////////////////////////////////
				LU_Beg = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
				U_Address[cbeg].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
        U2_Address[cbeg].first = 0;
#endif
				for (INMOST_DATA_ENUM_TYPE r = B_Address[cbeg].first + (B_Entries[B_Address[cbeg].first].first == cbeg ? 1 : 0); r != B_Address[cbeg].last; ++r) 
        {
          u = B_Entries[r].second / LU_Diag[cbeg];
          if( fabs(u) > tau )
					  LU_Entries.push_back(Sparse::Row::make_entry(B_Entries[r].first, u));
#if defined(ILUC2)
          else if( fabs(u) > tau2 )
            LU2_Entries.push_back(Sparse::Row::make_entry(B_Entries[r].first, u));
#endif
        }
        U_Address[cbeg].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());				
#if defined(ILUC2)
        U2_Address[cbeg].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());				
#endif
///////////////////////////////////////////////////////////////////////////////////
//       form first line of L and second-order L                                 //
///////////////////////////////////////////////////////////////////////////////////
        L_Address[cbeg].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
        L2_Address[cbeg].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());				
#endif
				Li = Bbeg[cbeg];
				while (Li != EOL)
				{
          l = B_Entries[B_Address[Li].first].second / LU_Diag[cbeg];
          if( fabs(l) > tau )
					  LU_Entries.push_back(Sparse::Row::make_entry(Li, l));
#if defined(ILUC2)
          else if( fabs(l) > tau2 )
            LU2_Entries.push_back(Sparse::Row::make_entry(Li, l));
#endif
					Li = Blist[Li];
				}
				L_Address[cbeg].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
        L2_Address[cbeg].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());				
#endif
				//update diagonal by factors
///////////////////////////////////////////////////////////////////////////////////
//       update whole diagonal with calculated L and U lines                     //
///////////////////////////////////////////////////////////////////////////////////
#if defined(TRACK_DIAGONAL)
				i = U_Address[cbeg].first;
				j = L_Address[cbeg].first;
				while (i != U_Address[cbeg].last && j != L_Address[cbeg].last)
				{
					Ui = LU_Entries[i].first;
					Li = LU_Entries[j].first;
					if (Ui > Li) j++;
					else if (Ui < Li) i++;
					else
					{
						assert(Ui > cbeg);
#if defined(DIAGONAL_PIVOT)
						mean_diag -= fabs(LU_Diag[Ui]);
#endif
						LU_Diag[Ui] -= LU_Entries[i].second * LU_Entries[j].second * LU_Diag[cbeg];
#if defined(DIAGONAL_PIVOT)
						mean_diag += fabs(LU_Diag[Ui]);
#endif
						i++;
						j++;
					}
				}
#endif
				//update column indexes
///////////////////////////////////////////////////////////////////////////////////
//       setup indices for transposed traversal of B block                       //
///////////////////////////////////////////////////////////////////////////////////
				Li = Bbeg[cbeg];
				while (Li != EOL)
				{
					B_Address[Li].first++;
					Ui = Blist[Li];
					if (B_Address[Li].Size() > 0)
					{
						i = B_Entries[B_Address[Li].first].first;
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
///////////////////////////////////////////////////////////////////////////////////
//       setup indices for transposed traversal of U                             //
///////////////////////////////////////////////////////////////////////////////////
				//form partial list from this row
				if (U_Address[cbeg].Size() > 0)
				{
					i = LU_Entries[U_Address[cbeg].first].first;
					Ubeg[i] = cbeg;
					Ulist[cbeg] = EOL;
				}
///////////////////////////////////////////////////////////////////////////////////
//       setup indices for transposed traversal of L                             //
///////////////////////////////////////////////////////////////////////////////////
				if (L_Address[cbeg].Size() > 0)
				{
					i = LU_Entries[L_Address[cbeg].first].first;
					Lbeg[i] = cbeg;
					Llist[cbeg] = EOL;
				}
				//Initialize data for condition estimator
				tt = Timer();
///////////////////////////////////////////////////////////////////////////////////
//              initialize condition estimator                                   //
///////////////////////////////////////////////////////////////////////////////////
#if defined(ESTIMATOR)
        if( estimator )
        {
				  NuU = NuL = 1.0;
				  if( estimator )
				  {
					  EstU1[cbeg] = EstL1[cbeg] = 0.0;
#if defined(ESTIMATOR_REFINE)
						EstU2[cbeg] = EstL2[cbeg] = 0.0;
#endif
					  for (INMOST_DATA_ENUM_TYPE it = U_Address[cbeg].first; it != U_Address[cbeg].last; ++it) 
						{
							EstU1[LU_Entries[it].first] = LU_Entries[it].second;
#if defined(ESTIMATOR_REFINE)
							EstU2[LU_Entries[it].first] = LU_Entries[it].second;
#endif
						}
					  for (INMOST_DATA_ENUM_TYPE it = L_Address[cbeg].first; it != L_Address[cbeg].last; ++it) 
						{
							EstL1[LU_Entries[it].first] = LU_Entries[it].second;
#if defined(ESTIMATOR_REFINE)
							EstL2[LU_Entries[it].first] = LU_Entries[it].second;
#endif
						}
				  }
        }
#endif
				nzLU += L_Address[cbeg].Size() + U_Address[cbeg].Size() + 1;
				max_diag = min_diag = fabs(LU_Diag[cbeg]);
				NuD = 1.0;
#if defined(REPORT_ILU)
				std::cout << " starting factorization " << std::endl;
#endif
				for (k = cbeg + 1; k < cend; ++k)
				{
///////////////////////////////////////////////////////////////////////////////////
//              starting factorization step                                      //
///////////////////////////////////////////////////////////////////////////////////
					//observe that there is a better diagonal pivot
#if defined(DIAGONAL_PIVOT)
					/*
					double mean_diag2 = 0.0;
					for (i = k; i < cend; i++)
					{
						mean_diag2 += fabs(LU_Diag[i]);
					}
					std::cout << k << "/" << cend << " " << mean_diag << " " << mean_diag2 << std::endl;
					*/
///////////////////////////////////////////////////////////////////////////////////
//              diagonal pivoting algorithm                                      //
///////////////////////////////////////////////////////////////////////////////////
#if defined(DIAGONAL_PIVOT_COND)
          int no_swap_algorithm = 1;
#endif
#if defined(DIAGONAL_PIVOT_TAU)
					if ( k < cend-1 &&  fabs(LU_Diag[k])*(cend-k) < DIAGONAL_PIVOT_TAU*mean_diag)
#else
					if( false )
#endif
					{
#if defined(DIAGONAL_PIVOT_COND)
            if( false )
            {
swap_algorithm:
              no_swap_algorithm--;
            }
#endif
						j = k+1;
						for (i = k + 1; i < cend; i++) if (fabs(LU_Diag[j]) < fabs(LU_Diag[i])) j = i;

            if( j > cend-1 )
            {
              std::cout << "Oops! asking to exchange beyond matrix! k = " << k << " j = " << j << " cend = " << cend << std::endl;
            }
						/*
            if( fabs(LU_Diag[j]) < fabs(LU_Diag[k]) )
            {
              std::cout << "Oops! I'm the best pivot! [" << j << "] " << LU_Diag[j] << " [" << k << "] " << LU_Diag[k] << " mean " << mean_diag/static_cast<INMOST_DATA_REAL_TYPE>(cend-k) << "\t\t\t" << std::endl;
              
            }
            else*/
						{
							//TODO!!!
              ++swaps;
#if defined(REPORT_ILU)
							//std::cout << "Detected, that there is a much better pivot, i'm " << k << " " << LU_Diag[k] << " other " << j << " " << LU_Diag[j] << std::endl;
							//std::cout << "Condition numbers: L " << NuL << " D " << NuD << " U " << NuU << std::endl;
#endif
							//scanf("%*c");
							//std::cout << "But repivoting algorithm not implemented" << std::endl;
							//This algorithm may be quite costly, but the effect is miraculous
							//First correct E
							//SwapEntries(*E_Address.back(), E_Entries, cend, wend, k, j);
							//Then correct F

							
							//SwapLine(F_Address,k, j);
							//Swap lines of B
							SwapLine(B_Address, k, j);
							//Correct B the same way E is corrected
							SwapEntries(B_Address, B_Entries, cbeg, cend, k, j);
							//now reintialize entrie column linked list
							for (i = k; i < cend; ++i) Bbeg[i] = EOL;
							for (i = cend; i > k; --i)
							{
								Li = B_Entries[B_Address[i - 1].first].first;
								if (Li < i - 1)
								{
									Blist[i - 1] = Bbeg[Li];
									Bbeg[Li] = i - 1;
								}
							}
							//now update L
							SwapEntries(L_Address, LU_Entries, cbeg, cend, k, j);
							//reinitialize linked list
							//now update U
							SwapEntries(U_Address, LU_Entries, cbeg, cend, k, j);
							//reinitialize linked list
#if defined(ILUC2)
							//now update L2
							SwapEntries(U2_Address, LU2_Entries, cbeg, cend, k, j);
							//reinitialize linked list
							//now update U2
							SwapEntries(L2_Address, LU2_Entries, cbeg, cend, k, j);
							//reinitialize linked list
#endif
							for (i = k; i < cend; i++)
							{
								Ubeg[i] = EOL;
								Lbeg[i] = EOL;
#if defined(ILUC2)
								U2beg[i] = EOL;
								L2beg[i] = EOL;
#endif
							}
							for (i = cbeg; i < k; i++)
							{
								if (U_Address[i].Size() > 0)
								{
									Li = LU_Entries[U_Address[i].first].first;
									Ulist[i] = Ubeg[Li];
									Ubeg[Li] = i;
								}
								if (L_Address[i].Size() > 0)
								{
									Li = LU_Entries[L_Address[i].first].first;
									Llist[i] = Lbeg[Li];
									Lbeg[Li] = i;
								}
#if defined(ILUC2)
								if (U2_Address[i].Size() > 0)
								{
									Li = LU2_Entries[U2_Address[i].first].first;
									U2list[i] = U2beg[Li];
									U2beg[Li] = i;
								}
								if (L2_Address[i].Size() > 0)
								{
									Li = LU2_Entries[L2_Address[i].first].first;
									L2list[i] = L2beg[Li];
									L2beg[Li] = i;
								}
#endif
							}
							//swap diagonal values
							u = LU_Diag[k];
							LU_Diag[k] = LU_Diag[j];
							LU_Diag[j] = u;
							//swap values in estimator
#if defined(ESTIMATOR)
							if( estimator )
							{
								u = EstL1[k];
								EstL1[k] = EstL1[j];
								EstL1[j] = u;
								u = EstU1[k];
								EstU1[k] = EstU1[j];
								EstU1[j] = u;
#if defined(ESTIMATOR_REFINE)
								u = EstU2[k];
								EstU2[k] = EstU2[j];
								EstU2[j] = u;
								u = EstL2[k];
								EstL2[k] = EstL2[j];
								EstL2[j] = u;
#endif
							}
#endif
							//swap rescaling vectors
							u = DL[k];
							DL[k] = DL[j];
							DL[j] = u;

							u = DR[k];
							DR[k] = DR[j];
							DR[j] = u;
							//update indeces in ddPQ
							Ui = ddP[k];
							ddP[k] = ddP[j];
							ddP[j] = Ui;
							Ui = ddQ[k];
							ddQ[k] = ddQ[j];
							ddQ[j] = Ui;

							//update interval starts for B block to correctly restore intervals
							Ui = localQ[k];
							localQ[k] = localQ[j];
							localQ[j] = Ui;

							Ui = localP[k];
							localP[k] = localP[j];
							localP[j] = Ui;
						}
						/*
						else
						{
							std::cout << "Detected that there are no good pivots anymore! interval " << cbeg << "..." << cend << " current element " << k << std::endl;
							std::cout << "But level rejection algorithm is not implemented" << std::endl;
							//probably should reject all unfactored parts for the next level
						}
						*/
					}
					
#endif
///////////////////////////////////////////////////////////////////////////////////
//            Prepare diagonal value                                             //
///////////////////////////////////////////////////////////////////////////////////
#if defined(TRACK_DIAGONAL)
					udiag = LU_Diag[k]; // LU_Diag is calculated with less dropping
///////////////////////////////////////////////////////////////////////////////////
// Prevent small pivot, hopefully it would work with correct right hand side     //
///////////////////////////////////////////////////////////////////////////////////
#if defined(PIVOT_THRESHOLD)
					if (fabs(udiag) < tol_modif)
					{
						udiag = udiag < 0.0 ? -tol_modif : tol_modif;
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
// Perturb diagonal value, this will shift eigenvalues                           //
///////////////////////////////////////////////////////////////////////////////////
#if defined(DIAGONAL_PERTURBATION)
					udiag = udiag * (1.0 + DIAGONAL_PERTURBATION_REL) + (udiag < 0.0 ? -1.0 : 1.0)*DIAGONAL_PERTURBATION_ABS;
#endif
					//abs_udiag = mean_diag/static_cast<INMOST_DATA_REAL_TYPE>(cend-k);//fabs(LU_Diag[k]);
					abs_udiag = fabs(LU_Diag[k]);
#else
					abs_udiag = 1.0;
#endif
///////////////////////////////////////////////////////////////////////////////////
//            Prepare linked list for U-part                                     //
///////////////////////////////////////////////////////////////////////////////////
					//DropLk = DropUk = 0.0;
					//uncompress k-th row
					// add diagonal value first, there shouldn't be values on left from diagonal
					assert(B_Entries[B_Address[k].first].first == k);
					LineIndecesU[cbeg] = k;
					if (B_Entries[B_Address[k].first].first == k)
						LineValuesU[k] = B_Entries[B_Address[k].first].second;
					else
					{
						std::cout << __LINE__ << " No diagonal value! " << k << " " << B_Entries[B_Address[k].first].first << std::endl;
						LineValuesU[k] = 0.0;
					}

					Ui = k;
					for (INMOST_DATA_ENUM_TYPE it = B_Address[k].first + (B_Entries[B_Address[k].first].first == k ? 1 : 0); it != B_Address[k].last; ++it)
					{
						LineValuesU[B_Entries[it].first] = B_Entries[it].second;
						Ui = LineIndecesU[Ui] = B_Entries[it].first;
					}
					LineIndecesU[Ui] = EOL;
///////////////////////////////////////////////////////////////////////////////////
//                    U-part elimination with L                                  //
///////////////////////////////////////////////////////////////////////////////////
					// perform the elimination from row
					i = Lbeg[k];
					while (i != EOL)
					{
						assert(i != UNDEF);
						assert(LU_Entries[L_Address[i].first].first == k);
						l = LU_Entries[L_Address[i].first].second*LU_Diag[i];
            //if( fabs(l)*NuU /*NuL_old*/ > tau2 *abs_udiag )//global check
#if defined(PREMATURE_DROPPING)
						if( fabs(l) > tau2*tau2 * abs_udiag )
#endif
            {
						  curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = U_Address[i].first; it < U_Address[i].last; ++it)
						  {
							  j = LU_Entries[it].first;
							  u = l*LU_Entries[it].second;
							  //eliminate values
							  if (LineIndecesU[j] != UNDEF) //there is an entry in the list
								  LineValuesU[j] -= u;
							  else //if (fabs(u)*NuU /*NuL_old*/ > tau2 *abs_udiag )//add new entry
							  {
								  next = curr;
								  while (next < j)
								  {
									  curr = next;
									  next = LineIndecesU[curr];
								  }
								  assert(curr < j);
								  assert(j < next);
								  LineIndecesU[curr] = j;
								  LineIndecesU[j] = next;
								  LineValuesU[j] = -u;
							  }
						  }
            }
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//                 second order U-part elimination with L                        //
///////////////////////////////////////////////////////////////////////////////////
						//if( fabs(l)*NuU /*NuL_old*/ > tau*abs_udiag )//global check
#if defined(PREMATURE_DROPPING)
						if( fabs(l) > tau2*abs_udiag )
#endif
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = U2_Address[i].first; it < U2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								u = l*LU2_Entries[it].second;
								//eliminate values
								if (LineIndecesU[j] != UNDEF) //there is an entry in the list
									LineValuesU[j] -= u;
								else //if (fabs(u)*NuU /*NuL_old*/ > tau2*abs_udiag)//add new entry
								{
									next = curr;
									while (next < j)
									{
										curr = next;
										next = LineIndecesU[curr];
									}
									assert(curr < j);
									assert(j < next);
									LineIndecesU[curr] = j;
									LineIndecesU[j] = next;
									LineValuesU[j] = -u;
								}
							}
						}
#endif
						i = Llist[i];
					}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//                 U-part elimination with second-order L                        //
///////////////////////////////////////////////////////////////////////////////////
					i = L2beg[k];
					while (i != EOL)
					{
						assert(i != UNDEF);
						assert(LU2_Entries[L2_Address[i].first].first == k);
						l = LU2_Entries[L2_Address[i].first].second*LU_Diag[i];
#if defined(PREMATURE_DROPPING)
						if( fabs(l) > tau2*abs_udiag )
#endif
            {
						  curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = U_Address[i].first; it < U_Address[i].last; ++it)
						  {
							  j = LU_Entries[it].first;
							  u = l*LU_Entries[it].second;
							  //eliminate values
							  if (LineIndecesU[j] != UNDEF) //there is an entry in the list
								  LineValuesU[j] -= u;
							  else //if (fabs(u)*NuU /*NuL_old*/ > tau2*abs_udiag)//add new entry
							  {
								  next = curr;
								  while (next < j)
								  {
									  curr = next;
									  next = LineIndecesU[curr];
								  }
								  assert(curr < j);
								  assert(j < next);
								  LineIndecesU[curr] = j;
								  LineIndecesU[j] = next;
								  LineValuesU[j] = -u;
							  }
						  }
///////////////////////////////////////////////////////////////////////////////////
//                 second order U-part elimination with second-order L           //
///////////////////////////////////////////////////////////////////////////////////
#if 0
              curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = U2_Address[i].first; it < U2_Address[i].last; ++it)
						  {
                j = LU2_Entries[it].first;
							  u = l*LU2_Entries[it].second;
							  //eliminate values
							  if (LineIndecesU[j] != UNDEF) //there is an entry in the list
								  LineValuesU[j] -= u;
                else //if (fabs(u)*NuU /*NuL_old*/ > tau2*abs_udiag)//add new entry
							  {
								  next = curr;
								  while (next < j)
								  {
									  curr = next;
									  next = LineIndecesU[curr];
								  }
								  assert(curr < j);
								  assert(j < next);
								  LineIndecesU[curr] = j;
								  LineIndecesU[j] = next;
								  LineValuesU[j] = -u;
							  }
						  }
#endif
						}
						i = L2list[i];
					}
#endif
					//std::cout << std::endl;
					//define the diagonal value
//this check will fail due to global check of tolerances, if global check will be removed then this will never fail
//#if defined(DIAGONAL_PIVOT) && !defined(DIAGONAL_PERTURBATION)
//					if (fabs(LU_Diag[k] - LineValues[k]) > 1e-6) std::cout << __LINE__ << " Diagonal value went wrong, good: " << LineValues[k] << " have " << LU_Diag[k] << " line " << k << std::endl;
//#endif
///////////////////////////////////////////////////////////////////////////////////
//                  Retrieve diagonal value                                      //
///////////////////////////////////////////////////////////////////////////////////
#if !defined(TRACK_DIAGONAL)
					udiag = LineValuesU[k];
///////////////////////////////////////////////////////////////////////////////////
// Prevent small pivot, hopefully it would work with correct right hand side     //
///////////////////////////////////////////////////////////////////////////////////
#if defined(PIVOT_THRESHOLD)
					if (fabs(udiag) < tol_modif)
					{
						//std::cout << "udiag too small " << 1.0 /udiag << std::endl;
						udiag = udiag < 0.0 ? -tol_modif : tol_modif;
					}
#endif
#if defined(DIAGONAL_PERTURBATION)
					udiag = udiag * (1.0 + DIAGONAL_PERTURBATION_REL) + (udiag < 0.0 ? -1.0 : 1.0)*DIAGONAL_PERTURBATION_ABS;
#endif
#endif
///////////////////////////////////////////////////////////////////////////////////
//                    Condition estimator for diagonal D                         //
///////////////////////////////////////////////////////////////////////////////////
          if (fabs(udiag) > max_diag) max_diag = fabs(udiag);
					if (fabs(udiag) < min_diag) min_diag = fabs(udiag);
					NuD = max_diag / min_diag;
///////////////////////////////////////////////////////////////////////////////////
//                    Rescale with diagonal                                     //
///////////////////////////////////////////////////////////////////////////////////
          //rescale elements next after diagonal position
					Ui = LineIndecesU[k];
					while (Ui != EOL)
					{
						LineValuesU[Ui] /= udiag;
						Ui = LineIndecesU[Ui];
					}
///////////////////////////////////////////////////////////////////////////////////
//                 prepearing linked list for L-part                             //
///////////////////////////////////////////////////////////////////////////////////
					//uncompress column
					//insert diagonal value first
					LineIndecesL[cbeg] = k;
					if (B_Entries[B_Address[k].first].first == k)
						LineValuesL[k] = B_Entries[B_Address[k].first].second;
					else
					{
						std::cout << __LINE__ << " No diagonal value! " << k << std::endl;
						LineValuesL[k] = 0.0;
					}
					//start from diagonal
					//sort_indeces.clear();
					Ui = k;
					Li = Bbeg[k];
					while (Li != EOL)
					{
						assert(B_Entries[B_Address[Li].first].first == k);
						//sort_indeces.push_back(Li);
						LineValuesL[Li] = B_Entries[B_Address[Li].first].second;
						Ui = LineIndecesL[Ui] = Li;
						Li = Blist[Li];
					}
					/*
					std::sort(sort_indeces.begin(), sort_indeces.end());
					for (Li = 0; Li < sort_indeces.size(); Li++)
					Ui = LineIndeces[Ui] = sort_indeces[Li];
					*/
					LineIndecesL[Ui] = EOL;
///////////////////////////////////////////////////////////////////////////////////
//                 L-part elimination with U                                     //
///////////////////////////////////////////////////////////////////////////////////
					// perform the elimination from column
					i = Ubeg[k];
					while (i != EOL)
					{
						assert(i != UNDEF);
						assert(LU_Entries[U_Address[i].first].first == k);
						u = LU_Entries[U_Address[i].first].second*LU_Diag[i];
            //if( fabs(u)*NuL /*NuU_old*/ > tau2*abs_udiag )//global check
#if defined(PREMATURE_DROPPING)
						if( fabs(u) > tau2*tau2*abs_udiag )
#endif
            {
						  curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = L_Address[i].first; it < L_Address[i].last; ++it)
						  {
							  j = LU_Entries[it].first;
							  l = u*LU_Entries[it].second;
							  //eliminate values
							  if (LineIndecesL[j] != UNDEF) //there is an entry in the list
								  LineValuesL[j] -= l;
							  else //if (fabs(l)*NuL /*NuU_old*/ > tau2*abs_udiag)//add new entry
							  {
								  next = curr;
								  while (next < j)
								  {
									  curr = next;
									  next = LineIndecesL[curr];
								  }
								  assert(curr < j);
								  assert(j < next);
								  LineIndecesL[curr] = j;
								  LineIndecesL[j] = next;
								  LineValuesL[j] = -l;
							  }
						  }
            }
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//                 second-order L-part elimination with U                        //
///////////////////////////////////////////////////////////////////////////////////
						//if( fabs(u)*NuL /*NuU_old*/ > tau*abs_udiag )//global check
#if defined(PREMATURE_DROPPING)
						if( fabs(u) > tau2*abs_udiag )
#endif
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = L2_Address[i].first; it < L2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								l = u*LU2_Entries[it].second;
								if (LineIndecesL[j] != UNDEF) //there is an entry in the list
									LineValuesL[j] -= l;
								else //if (fabs(l)*NuL /*NuU_old*/ > tau2*abs_udiag)//add new entry
								{
									next = curr;
									while (next < j)
									{
										curr = next;
										next = LineIndecesL[curr];
									}
									assert(curr < j);
									assert(j < next);
									LineIndecesL[curr] = j;
									LineIndecesL[j] = next;
									LineValuesL[j] = -l;
								}
							}
						}
#endif
						i = Ulist[i];
					}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//                 L-part elimination with second-order U                        //
///////////////////////////////////////////////////////////////////////////////////
					i = U2beg[k];
					while (i != EOL)
					{
						assert(i != UNDEF);
						assert(LU2_Entries[U2_Address[i].first].first == k);
						u = LU2_Entries[U2_Address[i].first].second*LU_Diag[i];
            //if( fabs(u)*NuL /*NuU_old*/ > tau2*abs_udiag )//global check
#if defined(PREMATURE_DROPPING)
						if( fabs(u) > tau2*abs_udiag )
#endif
            {
						  curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = L_Address[i].first; it != L_Address[i].last; ++it)
						  {
							  j = LU_Entries[it].first;
							  l = u*LU_Entries[it].second;
							  //eliminate values
							  if (LineIndecesL[j] != UNDEF) //there is an entry in the list
								  LineValuesL[j] -= l;
							  else //if (fabs(l)*NuL /*NuU_old*/ > tau2*abs_udiag)//add new entry
							  {
								  next = curr;
								  while (next < j)
								  {
									  curr = next;
									  next = LineIndecesL[curr];
								  }
								  assert(curr < j);
								  assert(j < next);
								  LineIndecesL[curr] = j;
								  LineIndecesL[j] = next;
								  LineValuesL[j] = -l;
                }
						  }
///////////////////////////////////////////////////////////////////////////////////
//           second-order L-part elimination with second-order U                 //
///////////////////////////////////////////////////////////////////////////////////
#if 0
              curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = L2_Address[i].first; it != L2_Address[i].last; ++it)
						  {
							  j = LU2_Entries[it].first;
							  l = u*LU2_Entries[it].second;
							  //eliminate values
							  if (LineIndecesL[j] != UNDEF) //there is an entry in the list
								  LineValuesL[j] -= l;
                else //if (fabs(l)*NuL /*NuU_old*/ > tau2 *abs_udiag)//add new entry
							  {
								  next = curr;
								  while (next < j)
								  {
									  curr = next;
									  next = LineIndecesL[curr];
								  }
								  assert(curr < j);
								  assert(j < next);
								  LineIndecesL[curr] = j;
								  LineIndecesL[j] = next;
								  LineValuesL[j] = -l;
                }
						  }
#endif
						}
						i = U2list[i];
					}
#endif
					//check that diagonal value is the same(must be!)
					//assert(fabs(LineValues[k] / udiag - 1.0) < 1.0e-10);
//this check will fail due to global check of tolerances, if global check will be removed then this will never fail
//#if defined(DIAGONAL_PIVOT) && !defined(DIAGONAL_PERTURBATION)
//					if (fabs(LineValues[k] - udiag) > 1e-9) std::cout << __LINE__ << " Diagonal value went wrong, good: " << udiag << " have " << LineValues[k] << " line " << k << std::endl;
//#endif
///////////////////////////////////////////////////////////////////////////////////
//                    Rescale with diagonal                                      //
///////////////////////////////////////////////////////////////////////////////////
					//rescale line by diagonal
					Li = LineIndecesL[k];
					while (Li != EOL)
					{
						LineValuesL[Li] /= udiag;
						Li = LineIndecesL[Li];
					}
///////////////////////////////////////////////////////////////////////////////////
//                    Condition estimator for U part                             //
///////////////////////////////////////////////////////////////////////////////////
					//estimate condition for U^{-1}
#if defined(ESTIMATOR)
					if( estimator )
					{
						NuU1_old = NuU1;
						mup = 1.0 - EstU1[k];
						mum = -1.0 - EstU1[k];
						smup = smum = 0;
						Ui = LineIndecesU[k];
						while (Ui != EOL)
						{
							v = EstU1[Ui];
							vp = fabs(v + LineValuesU[Ui] * mup);
							vm = fabs(v + LineValuesU[Ui] * mum);
							smup += vp;
							smum += vm;
							Ui = LineIndecesU[Ui];
						}
						if (smup > smum) NuU1 = mup; else NuU1 = mum;
						NuU1_new = std::max(fabs(mup),fabs(mum));
#if defined(ESTIMATOR_REFINE)
						NuU2_old = NuU2;
						mup = 1.0 - EstU2[k];
						mum = -1.0 - EstU2[k];
						np = nm = 0;
						//start from the element next after diagonal position
						Ui = LineIndecesU[k];
						while (Ui != EOL)
						{
							v = EstU2[Ui];
							vp = fabs(v + LineValuesU[Ui] * mup);
							vm = fabs(v + LineValuesU[Ui] * mum);
							v = fabs(v);
							if (vp > std::max(2 * v, 0.5)) np++;
							if (std::max(2 * vp, 0.5) < v) np--;
							if (vm > std::max(2 * v, 0.5)) nm++;
							if (std::max(2 * vm, 0.5) < v) nm--;
							Ui = LineIndecesU[Ui];
						}
						if (np > nm) NuU2 = mup; else NuU2 = mum;
						NuU2_new = std::max(fabs(mup), fabs(mum));
#endif
					} //if( estimator )
#endif //ESTIMATOR
///////////////////////////////////////////////////////////////////////////////////
//                 condition estimation for L part                               //
///////////////////////////////////////////////////////////////////////////////////
					//estimate condition for L^{-1}
#if defined(ESTIMATOR)
					if( estimator )
					{
						NuL1_old = NuL1;
						mup = 1.0 - EstL1[k];
						mum = -1.0 - EstL1[k];
						smup = smum = 0;
						Li = LineIndecesL[k];
						while (Li != EOL)
						{
							v = EstL1[Li];
							vp = fabs(v + LineValuesL[Li] * mup);
							vm = fabs(v + LineValuesL[Li] * mum);
							smup += vp;
							smum += vm;
							Li = LineIndecesL[Li];
						}
						if (smup > smum) NuL1 = mup; else NuL1 = mum;
						NuL1_new = std::max(fabs(mup),fabs(mum));
#if defined(ESTIMATOR_REFINE)
						NuL2_old = NuL2;
						mup = 1.0 - EstL2[k];
						mum = -1.0 - EstL2[k];
						np = nm = 0;
						//start from the element next after diagonal position
						Li = LineIndecesL[k];
						while (Li != EOL)
						{
							v = EstL2[Li];
							vp = fabs(v + LineValuesL[Li] * mup);
							vm = fabs(v + LineValuesL[Li] * mum);
							v = fabs(v);
							if (vp > std::max(2 * v, 0.5)) np++;
							if (std::max(2 * vp, 0.5) < v) np--;
							if (vm > std::max(2 * v, 0.5)) nm++;
							if (std::max(2 * vm, 0.5) < v) nm--;
							Li = LineIndecesL[Li];
						}
						if (np > nm) NuL2 = mup; else NuL2 = mum;
						NuL2_new = std::max(fabs(mup), fabs(mum));
#endif
					}
#endif //ESTIMATOR
///////////////////////////////////////////////////////////////////////////////////
//         discarding current iteration based on condition numbers               //
///////////////////////////////////////////////////////////////////////////////////
#if defined(ESTIMATOR)
          //discard line if condition number is too large
#if defined(DIAGONAL_PIVOT) && defined(DIAGONAL_PIVOT_COND)
#if defined(ESTIMATOR_REFINE)
					NuU_tmp = std::max(NuU1_new,NuU2_new);
					NuL_tmp = std::max(NuL1_new,NuL2_new);
#else
					NuU_tmp = NuU1_new;
					NuL_tmp = NuL1_new;
#endif
          if( k != cend-1 && no_swap_algorithm && (NuU_tmp > DIAGONAL_PIVOT_COND || NuL_tmp > DIAGONAL_PIVOT_COND) )
          {
#if defined(REPORT_ILU)
            std::cout << "Requested pivoting based on condition estimator (L" << NuL1_new << " " << NuL2_new << " U " << NuU1_new << " " << NuU2_new << ")! row " << k << "/" << cend << std::endl;
#endif //REPORT_ILU
            //restore condition number
						NuL1 = NuU1_old;
						NuU1 = NuU1_old;
#if defined(ESTIMATOR_REFINE)
						NuL2 = NuU2_old;
            NuU2 = NuU2_old;
#endif
            //clean up values
            Li = cbeg;
						while (Li != EOL)
						{
							i = LineIndecesL[Li];
							LineValuesL[Li] = 0.0; //clean values after use
							LineIndecesL[Li] = UNDEF; //clean indeces after use
							Li = i;
						}
						Ui = cbeg;
						while (Ui != EOL)
						{
							i = LineIndecesU[Ui];
							LineValuesU[Ui] = 0.0; // clean values after use
							LineIndecesU[Ui] = UNDEF; //clean indeces after use
							Ui = i;
						}
            goto swap_algorithm;
          }
          else
#endif //DIAGONAL_PIVOT / DIAGONAL_PIVOT_COND
          {
#if defined(ESTIMATOR_REFINE)
						NuL = std::max(NuL1_new,NuL2_new);
						NuU = std::max(NuU1_new,NuU2_new);
#else
						NuL = NuL1_new;
						NuU = NuU1_new;
#endif
						NuL_max = std::max(NuL,NuL_max);
						NuU_max = std::max(NuU,NuU_max);
          }
#endif //ESTIMATOR
///////////////////////////////////////////////////////////////////////////////////
//                 reconstruct U-part from linked list                           //
///////////////////////////////////////////////////////////////////////////////////
          //insert line to U part
					U_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
          U2_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
#endif
					Ui = LineIndecesU[k];
					while (Ui != EOL)
					{
						u = fabs(LineValuesU[Ui]);
						if (u*NuU /*NuL_old*/ > tau) // apply dropping rule
							LU_Entries.push_back(Sparse::Row::make_entry(Ui, LineValuesU[Ui]));
#if defined(ILUC2)
						else if (u*NuU /*NuL_old*/ > tau2)
							LU2_Entries.push_back(Sparse::Row::make_entry(Ui, LineValuesU[Ui]));
#endif
						Ui = LineIndecesU[Ui];
					}
					U_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
					U2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
#endif
///////////////////////////////////////////////////////////////////////////////////
//                 reconstruct L-part from linked list                           //
///////////////////////////////////////////////////////////////////////////////////
					//insert column to L part
					Li = LineIndecesL[k];
					L_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
          L2_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
#endif
					while (Li != EOL)
					{
						u = fabs(LineValuesL[Li]);
						if (u*NuL /*NuU_old*/ > tau) //apply dropping 
							LU_Entries.push_back(Sparse::Row::make_entry(Li, LineValuesL[Li]));
#if defined(ILUC2)
						else if (u*NuL /*NuU_old*/ > tau2)
							LU2_Entries.push_back(Sparse::Row::make_entry(Li, LineValuesL[Li]));
#endif
						Li = LineIndecesL[Li];
					}
					L_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
          L2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
#endif
///////////////////////////////////////////////////////////////////////////////////
//                 Update estimator vectors                                      //
///////////////////////////////////////////////////////////////////////////////////
#if defined(ESTIMATOR)
					if( estimator )
					{
						{ //U-estimator
							Ui = LineIndecesU[k];
							while (Ui != EOL)
							{
								EstU1[Ui] += LineValuesU[Ui] * NuU1;
								Ui = LineIndecesU[Ui];
							}
							//choose new condition number
							NuU1 = NuU1_new;//std::max(NuU_new,NuU_old);
							//update estimator vector
#if defined(ESTIMATOR_REFINE)
							Ui = LineIndecesU[k];
							while (Ui != EOL)
							{
								EstU2[Ui] += LineValuesU[Ui] * NuU2;
								Ui = LineIndecesU[Ui];
							}
							//choose new condition number
							NuU2 = NuU2_new;//std::max(NuU_new,NuU_old);
#endif
						}
						{ //L-estimator
							Li = LineIndecesL[k];
						  while (Li != EOL)
						  {
							  EstL1[Li] += LineValuesL[Li] * NuL1;
							  Li = LineIndecesL[Li];
						  }
							NuL1 = NuL1_new;
              //update estimator vector
#if defined(ESTIMATOR_REFINE)
              Li = LineIndecesL[k];
						  while (Li != EOL)
						  {
							  EstL2[Li] += LineValuesL[Li] * NuL2;
							  Li = LineIndecesL[Li];
						  }
              //choose new condition number
              NuL2 = NuL2_new;//std::max(NuL_new,NuL_old);
#endif
						}
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
//                 Cleaning data structures                                      //
///////////////////////////////////////////////////////////////////////////////////
					Li = cbeg;
					while (Li != EOL)
					{
						i = LineIndecesL[Li];
						LineValuesL[Li] = 0.0; //clean values after use
						LineIndecesL[Li] = UNDEF; //clean indeces after use
						Li = i;
					}
					Ui = cbeg;
					while (Ui != EOL)
					{
						i = LineIndecesU[Ui];
						LineValuesU[Ui] = 0.0; // clean values after use
						LineIndecesU[Ui] = UNDEF; //clean indeces after use
						Ui = i;
          }
///////////////////////////////////////////////////////////////////////////////////
//     Update diagonal                                                           //
///////////////////////////////////////////////////////////////////////////////////
#if defined(DIAGONAL_PIVOT)
          mean_diag -= fabs(LU_Diag[k]); //remove current value
#endif
					LU_Diag[k] = udiag;
///////////////////////////////////////////////////////////////////////////////////
//                 Insert indexes for transposed U-part traversal                //
///////////////////////////////////////////////////////////////////////////////////
          //Update linked list for factors
					if (U_Address[k].Size() != 0)
					{
						i = LU_Entries[U_Address[k].first].first;
						assert(i > k);
						Ulist[k] = Ubeg[i];
						Ubeg[i] = k;
					}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//         Insert indexes for transposed second-order U-part traversal           //
///////////////////////////////////////////////////////////////////////////////////
					if (U2_Address[k].Size() != 0)
					{
						i = LU2_Entries[U2_Address[k].first].first;
						assert(i > k);
						U2list[k] = U2beg[i];
						U2beg[i] = k;
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
//                 Insert indexes for transposed L-part traversal                //
///////////////////////////////////////////////////////////////////////////////////
          //update linked list for factors
					if (L_Address[k].Size() > 0)
					{
						i = LU_Entries[L_Address[k].first].first;
						assert(i > k);
						Llist[k] = Lbeg[i];
						Lbeg[i] = k;
					}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//         Insert indexes for transposed second-order L-part traversal           //
///////////////////////////////////////////////////////////////////////////////////
					if (L2_Address[k].Size() != 0)
					{
						i = LU2_Entries[L2_Address[k].first].first;
						assert(i > k);
						L2list[k] = L2beg[i];
						L2beg[i] = k;
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
//     Shift indexes for transposed U-part traversal for next iteration          //
///////////////////////////////////////////////////////////////////////////////////
					//Update indeces in linked lists
					i = Ubeg[k];
					while (i != EOL)
					{
						U_Address[i].first++;
						Li = Ulist[i];
						if (U_Address[i].Size() > 0)
						{
							Ui = LU_Entries[U_Address[i].first].first;
							Ulist[i] = Ubeg[Ui];
							Ubeg[Ui] = i;
						}
						i = Li;
					}
///////////////////////////////////////////////////////////////////////////////////
//     Shift indexes for transposed L-part traversal for next iteration          //
///////////////////////////////////////////////////////////////////////////////////
					i = Lbeg[k];
					while (i != EOL)
					{
						L_Address[i].first++;
						Li = Llist[i];
						if (L_Address[i].Size() > 0)
						{
							Ui = LU_Entries[L_Address[i].first].first;
							Llist[i] = Lbeg[Ui];
							Lbeg[Ui] = i;
						}
						i = Li;
					}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//         Shift indexes for transposed second-order U-part traversal            //
///////////////////////////////////////////////////////////////////////////////////
					i = U2beg[k];
					while (i != EOL)
					{
						U2_Address[i].first++;
						Li = U2list[i];
						if (U2_Address[i].Size() > 0)
						{
							Ui = LU2_Entries[U2_Address[i].first].first;
							U2list[i] = U2beg[Ui];
							U2beg[Ui] = i;
						}
						i = Li;
					}
///////////////////////////////////////////////////////////////////////////////////
//         Shift indexes for transposed second-order L-part traversal            //
///////////////////////////////////////////////////////////////////////////////////
					i = L2beg[k];
					while (i != EOL)
					{
						L2_Address[i].first++;
						Li = L2list[i];
						if (L2_Address[i].Size())
						{
							Ui = LU2_Entries[L2_Address[i].first].first;
							L2list[i] = L2beg[Ui];
							L2beg[Ui] = i;
						}
						i = Li;
					}
					nzLU2 += U2_Address[k].Size() + L2_Address[k].Size();
#endif
///////////////////////////////////////////////////////////////////////////////////
//         Shift indexes for transposed B matrix traversal                       //
///////////////////////////////////////////////////////////////////////////////////
					//update column
					Li = Bbeg[k];
					while (Li != EOL)
					{
						B_Address[Li].first++;
						Ui = Blist[Li];
						if (B_Address[Li].Size() > 0)
						{
							i = B_Entries[B_Address[Li].first].first;
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
///////////////////////////////////////////////////////////////////////////////////
//         Update values on the whole diagonal with L and U                      //
///////////////////////////////////////////////////////////////////////////////////
#if defined(TRACK_DIAGONAL)
					//update diagonal by optained factors
					i = U_Address[k].first;
					j = L_Address[k].first;
					while (i != U_Address[k].last && j != L_Address[k].last)
					{
						Ui = LU_Entries[i].first;
						Li = LU_Entries[j].first;
						if (Ui > Li) j++;
						else if (Ui < Li) i++;
						else
						{
							assert(Ui > k);
#if defined(DIAGONAL_PIVOT)
							mean_diag -= fabs(LU_Diag[Ui]);
#endif
							LU_Diag[Ui] -= LU_Entries[i].second * LU_Entries[j].second * LU_Diag[k];
#if defined(DIAGONAL_PIVOT)
							mean_diag += fabs(LU_Diag[Ui]);
#endif
							i++;
							j++;
						}
					}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//         Update values on the whole diagonal with L and second-order U         //
///////////////////////////////////////////////////////////////////////////////////
					i = U_Address[k].first;
					j = L2_Address[k].first;
					while (i != U_Address[k].last && j != L2_Address[k].last)
					{
						Ui = LU_Entries[i].first;
						Li = LU2_Entries[j].first;
						if (Ui > Li) j++;
						else if (Ui < Li) i++;
						else
						{
							assert(Ui > k);
#if defined(DIAGONAL_PIVOT)
							mean_diag -= fabs(LU_Diag[Ui]);
#endif
							LU_Diag[Ui] -= LU_Entries[i].second * LU2_Entries[j].second * LU_Diag[k];
#if defined(DIAGONAL_PIVOT)
							mean_diag += fabs(LU_Diag[Ui]);
#endif
							i++;
							j++;
						}
					}
///////////////////////////////////////////////////////////////////////////////////
//         Update values on the whole diagonal with second-order L and U         //
///////////////////////////////////////////////////////////////////////////////////
					i = U2_Address[k].first;
					j = L_Address[k].first;
					while (i != U2_Address[k].last && j != L_Address[k].last)
					{
						Ui = LU2_Entries[i].first;
						Li = LU_Entries[j].first;
						if (Ui > Li) j++;
						else if (Ui < Li) i++;
						else
						{
							assert(Ui > k);
#if defined(DIAGONAL_PIVOT)
							mean_diag -= fabs(LU_Diag[Ui]);
#endif
							LU_Diag[Ui] -= LU2_Entries[i].second * LU_Entries[j].second * LU_Diag[k];
#if defined(DIAGONAL_PIVOT)
							mean_diag += fabs(LU_Diag[Ui]);
#endif
							i++;
							j++;
						}
					}
///////////////////////////////////////////////////////////////////////////////////
// Update values on the whole diagonal with second-order L and second-order U    //
///////////////////////////////////////////////////////////////////////////////////
#if 0
					i = U2_Address[k].first;
					j = L2_Address[k].first;
					while (i != U2_Address[k].last && j != L2_Address[k].last)
					{
						Ui = LU2_Entries[i].first;
						Li = LU2_Entries[j].first;
						if (Ui > Li) j++;
						else if (Ui < Li) i++;
						else
						{
							assert(Ui > k);
#if defined(DIAGONAL_PIVOT)
							mean_diag -= fabs(LU_Diag[Ui]);
#endif
							LU_Diag[Ui] -= LU2_Entries[i].second * LU2_Entries[j].second * LU_Diag[k];
#if defined(DIAGONAL_PIVOT)
							mean_diag += fabs(LU_Diag[Ui]);
#endif
							i++;
							j++;
						}
					}
#endif
#endif //ILUC2
#endif //TRACK_DIAGONAL
//#endif //DIAGONAL_PIVOT

					//number of nonzeros
					nzLU += U_Address[k].Size() + L_Address[k].Size() + 1;
					
#if defined(REPORT_ILU)
					if (k % 2000 == 0)
					{
						
						std::cout << std::fixed << std::setprecision(2) << std::setw(6) << 100.0f*(k - cbeg) / (float)(cend - cbeg) << "%";
						std::cout << " nnz LU " << std::setprecision(10) << std::setw(10) << nzLU;
#if defined(ILUC2)
						std::cout << " LU2 " << std::setw(10) << nzLU2;
#endif
						std::cout << std::scientific << std::setprecision(5) << " condition L " << NuL << " D " << NuD << " U " << NuU;
            std::cout << " pivot swaps " << swaps;
						std::cout << "\r"; //carrige return
						std::cout.flush();
            //std::cout << std::endl;
						std::cout << std::setprecision(0);
						//std::cout.setf(0,std::ios::floatfield);
						
						//printf("%6.2f%% nnz LU %8d condition L %10f D %10f U %10f\r", 100.0f*(k - cbeg) / (float)(cend - cbeg), nzLU, NuL, NuD, NuU);
						//fflush(stdout);
					}
#elif defined(REPORT_ILU_PROGRESS)
					if (k % 500 == 0)
					{
						printf("%lu %d/%d factor %6.2f%%\t\t\r",static_cast<unsigned long>(level_size.size()), cend,moend, 100.0f*(k - cbeg) / (float)(cend - cbeg));
						//printf("%6.2f%% nnz LU %8d condition L %10f D %10f U %10f\r", 100.0f*(k - cbeg) / (float)(cend - cbeg), nzLU, NuL, NuD, NuU);
						fflush(stdout);
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
//                       iteration done                                          //
///////////////////////////////////////////////////////////////////////////////////
				}
#if defined(REPORT_ILU_END)
        std::cout << "size " << moend-mobeg;
				std::cout << " total nonzeros in A " << nzA << " in LU " << nzLU;
#if defined(ILUC2)
				std::cout << " in LU2 " << nzLU2;
#endif
				std::cout << " conditions L " << NuL_max << " D " << NuD << " U " << NuU_max << " pivot swaps " << swaps << std::endl;
#endif
			}
			tlfactor = Timer() - tt;
			tfactor += tlfactor;
			//restore indexes
			B_Address[cbeg].first = localQ[cbeg];
			U_Address[cbeg].first = LU_Beg;
			L_Address[cbeg].first = U_Address[cbeg].last;
#if defined(ILUC2)
      U2_Address[cbeg].first = 0;
      L2_Address[cbeg].first = U2_Address[cbeg].last;
#endif
			for (k = cbeg+1; k < cend; ++k)
			{
				B_Address[k].first = localQ[k];
				U_Address[k].first = L_Address[k - 1].last;
				L_Address[k].first = U_Address[k].last;
#if defined(ILUC2)
				U2_Address[k].first = L2_Address[k - 1].last;
				L2_Address[k].first = U2_Address[k].last;
#endif
			}
			
			
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// FACTORIZATION COMPLETE ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//After we have factored the rescaled matrix we must rescale obtained factors
			tt = Timer();

#if defined(RESCALE_B)
#if defined(REPORT_ILU)
			std::cout << " rescaling block B back " << std::endl;
#endif
			for (k = cbeg; k < cend; ++k)
			{
				LU_Diag[k] /= DL[k] * DR[k];
				for (INMOST_DATA_ENUM_TYPE r = U_Address[k].first; r < U_Address[k].last; ++r) LU_Entries[r].second *= DR[k] / DR[LU_Entries[r].first];
				for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r < L_Address[k].last; ++r) LU_Entries[r].second *= DL[k] / DL[LU_Entries[r].first];
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) B_Entries[r].second /= DL[k] * DR[B_Entries[r].first];
#if defined(ILUC2)
				for (INMOST_DATA_ENUM_TYPE r = U2_Address[k].first; r < U2_Address[k].last; ++r) LU2_Entries[r].second *= DR[k] / DR[LU2_Entries[r].first];
				for (INMOST_DATA_ENUM_TYPE r = L2_Address[k].first; r < L2_Address[k].last; ++r) LU2_Entries[r].second *= DL[k] / DL[LU2_Entries[r].first];
#endif
			}
#endif
			tlrescale = Timer() - tt;
			trescale += tlrescale;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// RESCALING DONE   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			/*
//#if defined(REPORT_ILU)
			std::cout << "============================================================================" << std::endl;
			//std::cout << "factor: " << tlrescale + tlfactor << " trescale: " << tlrescale << " tfactor " << tlfactor << std::endl;
			//std::cout << "schur:  " << tlreorder + tlreassamble + tlschur << " treorder: " << tlreorder << " treassamble " << tlreassamble << " tschur: " << tlschur << std::endl;
			std::cout << " time ratio " << ratio << std::endl;
			std::cout << " time ratio2 " << ratio2 << std::endl;
			std::cout << " ddpq_tau " << ddpq_tau << std::endl;
			std::cout <<  cbeg << "..." << cend << " of " << mobeg << "..." << moend << " level " << level_size.size() << std::endl;
			std::cout << "============================================================================" << std::endl;
//#endif
			*/
			wend = wbeg;
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////  FACTORIZATION COMPLETE ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		ttotal = Timer() - ttotal;
#if defined(REPORT_ILU_SUMMARY)
		printf("total      %f\n",ttotal);
		printf("reorder    %f (%6.2f%%)\n", treorder, 100.0*treorder/ttotal);
#if defined(REORDER_METIS_ND)
		printf("metis      graph %f nd %f\n", tmetisgraph, tmetisnd);
#endif
		printf("reassamble %f (%6.2f%%)\n", treassamble, 100.0*treassamble / ttotal);
		printf("rescale    %f (%6.2f%%)\n", trescale, 100.0*trescale / ttotal);
		printf("factor     %f (%6.2f%%)\n", tfactor, 100.0*tfactor / ttotal);
#endif
		init = true;
		/*
		level_interval.resize(level_size.size());
		level_interval[0].first = mobeg;
		level_interval[0].last = mobeg + level_size[0];
		for (k = 1; k < level_size.size(); k++)
		{
			level_interval[k].first = level_interval[k - 1].last;
			level_interval[k].last = level_interval[k].first + level_size[k];
		}
		*/
    condestL = NuL;
    condestU = NuU;
		return true;
	}
	bool MTILUC_preconditioner::Finalize()
	{
		init = false;
		L_Address.clear();
		U_Address.clear();
		LU_Entries.clear();
		B_Entries.clear();
		B_Address.clear();
		temp.clear();
		ddP.clear();
		ddQ.clear();
		LU_Diag.clear();
		return true;
	}
	void MTILUC_preconditioner::ApplyB(double alpha, Sparse::Vector & x, double beta, Sparse::Vector & y) // y = alpha A x + beta y
	{
		INMOST_DATA_ENUM_TYPE k, m, cbeg, cend;
		INMOST_DATA_REAL_TYPE temp;
		info->GetOverlapRegion(info->GetRank(), cbeg, cend);
		for(k = cbeg; k < cend; ++k)
		{
			temp = 0.0;
			for(m = B_Address[k].first; m != B_Address[k].last; ++m)
			{
				temp += B_Entries[m].second*x[B_Entries[m].first];
			}
			y[k] = beta*y[k] + alpha * temp;
		}
	}
	bool MTILUC_preconditioner::Solve(Sparse::Vector & input, Sparse::Vector & output)
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
		
			for (k = vbeg; k < mobeg; k++) temp[k] = 0;
			for (k = mobeg; k < moend; ++k) temp[k] = input[ddP[k]];
			for (k = moend; k < vend; k++) temp[k] = 0;



			//Solve with L first
			for (k = mobeg; k < moend; ++k) //iterator over columns of L
			{
				for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r != L_Address[k].last; ++r)
					temp[LU_Entries[r].first] -= temp[k] * LU_Entries[r].second;
			}
			//Solve with diagonal
			for (k = mobeg; k < moend; ++k) temp[k] /= LU_Diag[k];
			//Solve with U
			for (k = moend; k > mobeg; --k) //iterator over rows of U
			{
				for (INMOST_DATA_ENUM_TYPE r = U_Address[k - 1].first; r != U_Address[k - 1].last; ++r)
					temp[k - 1] -= temp[LU_Entries[r].first] * LU_Entries[r].second;
			}

		
			for (k = mobeg; k < moend; ++k) output[ddQ[k]] = temp[k];

			//Restrict additive schwartz (maybe do it outside?)
			//May assamble partition of unity instead of restriction before accumulation
			//assembly should be done instead of initialization
			for (k = vbeg; k < mobeg; ++k) output[k] = 0;
			for (k = moend; k < vend; ++k) output[k] = 0;
		}
		info->Accumulate(output);
		return true;
	}
	bool MTILUC_preconditioner::ReplaceMAT(Sparse::Matrix & A){ if (isInitialized()) Finalize(); Alink = &A; return true; }
	bool MTILUC_preconditioner::ReplaceSOL(Sparse::Vector & x) {(void)x; return true; }
	bool MTILUC_preconditioner::ReplaceRHS(Sparse::Vector & b) {(void)b; return true; }
	Method * MTILUC_preconditioner::Duplicate() { return new MTILUC_preconditioner(*this); }
	MTILUC_preconditioner::~MTILUC_preconditioner()
	{
		if (!isFinalized()) Finalize();
	}

