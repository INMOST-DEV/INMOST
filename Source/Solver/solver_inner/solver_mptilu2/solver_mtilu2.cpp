#include "inmost_solver.h"
#if defined(USE_SOLVER)
#include "solver_mtilu2.hpp"

#define DEFAULT_TAU 0.005
#define DEFAULT_TAU2 0.00001
//#define LFILL //control, that factorization is not less then fill for ilu2

//select one of the two rescaling techniques
#define RESCALE_EQUALIZE_1NORM //equalize 1-norms of each row and each column to 1
//#define RESCALE_EQUALIZE_2NORM
//#define RESCALE_MAXIMUM_TRANSVERSAL //use rescaling produced by maximum product transversal algorithm that bounds all values between -1 and 1

#define REORDER_MAXIMUM_TRANSVERSAL

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
		Position[Array[0] - Base] = 0;
		Array.pop_back();
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


void MTILU2_preconditioner::DumpMatrix(interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & Address, 
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
			nnz += Address[k+1] - Address[k];
			vrowmax = 0;

			bool diag_found = false;
			diag = 0;
			for (INMOST_DATA_ENUM_TYPE it = Address[k]; it != Address[k+1]; ++it)
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
			for (INMOST_DATA_ENUM_TYPE it = Address[k]; it != Address[k+1]; ++it)
				fout << (k-wmbeg+1) << " " << (Entries[it].first-wmbeg+1) << " " << Entries[it].second << std::endl;
		}
		fout.close();
	}


  INMOST_DATA_REAL_TYPE & MTILU2_preconditioner::RealParameter(std::string name)
	{
		if( name == "tau" ) return tau;
		else if( name == "tau2" ) return tau2;
		throw -1;
	}
	INMOST_DATA_ENUM_TYPE & MTILU2_preconditioner::EnumParameter(std::string name)
	{
		if (name == "fill") return Lfill;
		else if (name == "scale_iters") return sciters;
		throw -1;
	}
	MTILU2_preconditioner::MTILU2_preconditioner(Solver::OrderInfo & info)
		:info(&info),tau(DEFAULT_TAU), tau2(DEFAULT_TAU2)
	{
		Alink = NULL;
		init = false;
		sciters = 12;
		Lfill = 1;
	}




  bool MTILU2_preconditioner::Initialize()
	{
		if (isInitialized()) Finalize();
		assert(Alink != NULL);
		double treorder = 0, trescale = 0, tfactor = 0, ttotal = 0;
		ttotal = Timer();

		Sparse::Matrix & A = *Alink;

		INMOST_DATA_ENUM_TYPE mobeg, moend, vlocbeg, vlocend, vbeg, vend, k, r, j;
		INMOST_DATA_REAL_TYPE leabs, flin, ldiag, udiag, mva;
		INMOST_DATA_ENUM_TYPE curr, foll;
		INMOST_DATA_ENUM_TYPE  ind, jn;
		std::vector<INMOST_DATA_REAL_TYPE> rv;
		std::vector<INMOST_DATA_ENUM_TYPE> ri;
		INMOST_DATA_INTEGER_TYPE prev, ipred;
		const INMOST_DATA_REAL_TYPE tol_modif = 1e-12, eps = 1.0e-54, subst = 1.0;
		const INMOST_DATA_ENUM_TYPE UNDEF = ENUMUNDEF, EOL = ENUMUNDEF - 1;
		//Calculate scaling vectors for matrix (from genebs)
		info->GetOverlapRegion(info->GetRank(), mobeg, moend);
		info->GetLocalRegion(info->GetRank(), vlocbeg, vlocend);
		info->GetVectorRegion(vbeg, vend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> ir(mobeg, moend + 1);
		interval<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_REAL_TYPE> RowValues(vbeg, vend);
		nnz = 0;
		for (Sparse::Matrix::iterator it = (*Alink).Begin(); it != (*Alink).End(); ++it) nnz += it->Size();
#if defined(LFILL)
		interval<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_ENUM_TYPE> RowFill(vbeg, vend);
		//std::fill(RowFill.begin(),RowFill.end(),ENUMUNDEF);
		std::vector<INMOST_DATA_ENUM_TYPE> lfill;
		lfill.reserve(nnz * 4);
#endif
		interval<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_ENUM_TYPE> RowIndeces(vbeg - 1, vend,UNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> B_Address(mobeg,moend+1);
		std::vector<Sparse::Row::entry> B_Entries(nnz);
		Perm.set_interval_beg(mobeg);
		Perm.set_interval_end(moend);
		Sparse::Vector DL("",mobeg,moend), DR("",mobeg,moend);
		std::fill(DL.Begin(),DL.End(),1.0);
		std::fill(DR.Begin(),DR.End(),1.0);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// MAXIMUM TRANSVERSE REORDERING /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined(REORDER_MAXIMUM_TRANSVERSAL)
		{
			treorder = Timer();
			INMOST_DATA_ENUM_TYPE ColumnBegin;
			INMOST_DATA_ENUM_TYPE pop_heap_pos;
			INMOST_DATA_REAL_TYPE ShortestPath, AugmentPath;
			INMOST_DATA_ENUM_TYPE PathEnd, Trace, IPermPrev;
			//Sparse::Vector & U = DL;
			//Sparse::Vector & V = DR;
			std::fill(Perm.begin(),Perm.end(), ENUMUNDEF);
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> U(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> V(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Dist(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Cmax(mobeg,moend,0.0);
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> IPerm(mobeg,moend,ENUMUNDEF);
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> UpdateStack(mobeg,moend);
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> ColumnList(mobeg,moend,ENUMUNDEF);
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Parent(mobeg,moend,ENUMUNDEF);
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> AugmentPosition(mobeg,moend,ENUMUNDEF);
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> ColumnPosition(mobeg,moend,ENUMUNDEF);
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & C_Address = B_Address;
			std::vector<Sparse::Row::entry> C_Entries(nnz);
			BinaryHeap Heap(&Dist[mobeg],moend-mobeg);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Arrays initialization   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				
			//std::fill(U.Begin() + wbeg - mobeg, U.Begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
			//std::fill(V.Begin() + wbeg - mobeg, V.Begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
			//std::fill(Cmax.begin() + wbeg - mobeg, Cmax.begin() + wend - mobeg, 0.0);
			//std::fill(Dist.begin() + wbeg - mobeg, Dist.begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
			
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Initial LOG transformation to dual problem and initial extreme match /////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			//double T = Timer();
			INMOST_DATA_ENUM_TYPE cnt = 0;
			C_Address[mobeg] = 0;
			for(k = mobeg; k < moend; ++k)
			{
				
				for (Sparse::Row::iterator it = A[k].Begin(); it != A[k].End(); ++it)
				{
          if( it->first >= mobeg && it->first < moend ) //filter any extension of the matrix
          {
					  INMOST_DATA_ENUM_TYPE i = C_Entries[cnt].first = it->first;
					  INMOST_DATA_REAL_TYPE u = C_Entries[cnt].second = fabs(it->second);
					  if( u > Cmax[i] ) Cmax[i] = u;
					  ++cnt;
          }
				}
				C_Address[k+1] = cnt;
			}

			for(k = mobeg; k < moend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE it = C_Address[k]; it < C_Address[k+1]; ++it)
				{
					INMOST_DATA_ENUM_TYPE i = C_Entries[it].first;
					if( Cmax[i] == 0.0 || C_Entries[it].second == 0.0 )
						C_Entries[it].second = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
					else
					{
						C_Entries[it].second = log(Cmax[i])-log(C_Entries[it].second);
						//C_Entries[it] = log10(Cmax[i])-log10(C_Entries[it]);
						//C_Entries[it] = fabs(log10(Cmax[i]/C_Entries[it]));
						//C_Entries[it] = fabs(log(Cmax[i]/C_Entries[it]));
						if( C_Entries[it].second < U[i] ) U[i] = C_Entries[it].second;
					}
				}
			}
			for(k = mobeg; k < moend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE it = C_Address[k]; it < C_Address[k+1]; ++it)
				{
					INMOST_DATA_REAL_TYPE u = C_Entries[it].second - U[C_Entries[it].first];
					if( u < V[k] ) V[k] = u;
				}
			}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Update cost and match ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for(k = mobeg; k < moend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE it = C_Address[k]; it < C_Address[k+1]; ++it)
				{
					INMOST_DATA_REAL_TYPE u = fabs(C_Entries[it].second - V[k] - U[C_Entries[it].first]);
					if( u < 1.0e-30 && Perm[C_Entries[it].first] == ENUMUNDEF && IPerm[k] == ENUMUNDEF )
					{
							Perm[C_Entries[it].first] = k;
							IPerm[k] = C_Entries[it].first;
							ColumnPosition[k] = it;
					}
				}
			}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// 1-step augmentation   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			for(k = mobeg; k < moend; ++k)
			{
				if( IPerm[k] == ENUMUNDEF ) //unmatched row
				{
					for (INMOST_DATA_ENUM_TYPE it = C_Address[k]; it < C_Address[k+1] && IPerm[k] == ENUMUNDEF; ++it)
					{
						INMOST_DATA_REAL_TYPE u = fabs(C_Entries[it].second - V[k] - U[C_Entries[it].first]);
						if( u <= 1.0e-30 )
						{
							INMOST_DATA_ENUM_TYPE Li = Perm[C_Entries[it].first];
							assert(Li != ENUMUNDEF);
							// Search other row in C for 0
							for (INMOST_DATA_ENUM_TYPE Lit = C_Address[Li]; Lit < C_Address[Li+1]; ++Lit)
							{
								u = fabs(C_Entries[Lit].second - V[Li] - U[C_Entries[Lit].first]);
								if( u <= 1.0e-30 && Perm[C_Entries[Lit].first] == ENUMUNDEF )
								{
									Perm[C_Entries[it].first] = k;
									IPerm[k] = C_Entries[it].first;
									ColumnPosition[k] = it;
									Perm[C_Entries[Lit].first] = Li;
									IPerm[Li] = C_Entries[Lit].first;
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
			for(k = mobeg; k < moend; ++k)
			{
				if( IPerm[k] != ENUMUNDEF )
					continue;
				INMOST_DATA_ENUM_TYPE Li = k;
				ColumnBegin = EOL;
				Parent[Li] = ENUMUNDEF;
				PathEnd = ENUMUNDEF;
				Trace = k;
				ShortestPath = 0;
				AugmentPath = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
				while(true)
				{
					for (INMOST_DATA_ENUM_TYPE Lit = C_Address[Li]; Lit < C_Address[Li+1]; ++Lit)
					{
						INMOST_DATA_ENUM_TYPE Ui = C_Entries[Lit].first;
						//if( ColumnList[Ui] == k ) continue;
						if( ColumnList[Ui] != ENUMUNDEF ) continue;
						INMOST_DATA_REAL_TYPE l = ShortestPath + C_Entries[Lit].second - V[Li] - U[Ui];
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
								if( Heap.GetPosition(Ui-mobeg) != ENUMUNDEF )
									Heap.DecreaseKey(Ui-mobeg);
								else 
									Heap.PushHeap(&Dist[Ui]);
							}
						}
					}

					pop_heap_pos = Heap.PopHeap();
					if( pop_heap_pos == ENUMUNDEF ) break;
					
					INMOST_DATA_ENUM_TYPE Ui = pop_heap_pos+mobeg;
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
					INMOST_DATA_ENUM_TYPE Ui = ColumnBegin;
					while(Ui != EOL)
					{
						U[Ui] += Dist[Ui] - AugmentPath;
						if( Perm[Ui] != ENUMUNDEF ) V[Perm[Ui]] = C_Entries[ColumnPosition[Perm[Ui]]].second - U[Ui];
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
						V[Trace] = C_Entries[ColumnPosition[Trace]].second - U[Ui];

						Ui = IPermPrev;
						Trace = Parent[Trace];

					}
					for(Ui = 0; Ui < Heap.GetSize(); ++Ui) *Heap.Get(Ui) = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
					Heap.Clear();
				}
			}

			//for(k = mobeg; k != moend; ++k)
			//if( Perm[k] != k ) std::cout << "Column " << k << " go to " << Perm[k] << " inverse " << IPerm[k] << std::endl;
            
            { //fill gaps in perm
                std::fill(IPerm.begin(), IPerm.end(), ENUMUNDEF);
                for(k = mobeg; k < moend; ++k)
                {
                    if( Perm[k] != ENUMUNDEF )
                        IPerm[Perm[k]] = 0;
                }
                std::vector<INMOST_DATA_ENUM_TYPE> gaps;
                for(k = mobeg; k < moend; ++k)
                    if( IPerm[k] == ENUMUNDEF )
                        gaps.push_back(k);
                
                for(k = mobeg; k < moend; ++k)
                    if( Perm[k] == ENUMUNDEF )
                    {
                        Perm[k] = gaps.back();
                        gaps.pop_back();
                    }
            }
            
			

			cnt = 0;
#if defined(RESCALE_MAXIMUM_TRANSVERSAL)
            trescale = Timer();
            for (k = mobeg; k < moend; ++k)
            {
                bool flip_sign = false;
                for(INMOST_DATA_ENUM_TYPE qt = 0; qt < A[k].Size(); ++qt)
                {
                    INMOST_DATA_ENUM_TYPE i = A[k].GetIndex(qt), j = Perm[i];
					INMOST_DATA_REAL_TYPE l,u;
					if( V[k] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() ) l = 1;
					else l = exp(V[k]);
					if( U[i] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() || Cmax[i] == 0 ) u = 1;
					else u = exp(U[i])/Cmax[i];
                    DL[k] = l;
                    DR[i] = u;
                    B_Entries[cnt++] = 	Sparse::Row::make_entry(j, l*u*A[k].GetValue(qt));
                    if( j == k && (B_Entries[cnt-1].second) < 0.0 ) flip_sign = true;
                }
                if( flip_sign )
                {
                    DL[k] *= -1;
                    for(INMOST_DATA_ENUM_TYPE Li = B_Address[k]; Li < B_Address[k+1]; ++Li) B_Entries[Li].second *= -1;
                }
                //std::sort(B_Entries.begin()+C_Address[k],B_Entries.begin()+C_Address[k+1]);
            }
            for (k = mobeg; k < moend; ++k) U[Perm[k]] = DR[k];
            for (k = mobeg; k < moend; ++k) DR[k] = U[k];
            
            trescale = Timer() - trescale;
#else
            for (k = mobeg; k < moend; ++k)
            {
                for(INMOST_DATA_ENUM_TYPE qt = 0; qt < A[k].Size(); ++qt)
                {
                    if( A[k].GetIndex(qt) >= mobeg && A[k].GetIndex(qt) < moend )
                        B_Entries[cnt++] = 	Sparse::Row::make_entry(Perm[A[k].GetIndex(qt)], A[k].GetValue(qt));
                }
                //std::sort(B_Entries.begin()+C_Address[k],B_Entries.begin()+C_Address[k+1]);
            }
#endif

			//for (k = mobeg; k < moend; ++k) IPerm[Perm[k]] = k; //inversePQ
			//for (k = mobeg; k < moend; ++k) ColumnPosition[k] = k;
			//for (k = mobeg; k < moend; ++k) Perm[k] = ColumnPosition[IPerm[k]]; //applyPQ
			//DumpMatrix(B_Address,B_Entries,mobeg,moend,"MC64.mtx");
			treorder = Timer() - treorder;
		}
#else
		{
            treorder = Timer();
            INMOST_DATA_ENUM_TYPE cnt = 0;
            B_Address[mobeg] = 0;
            for(k = mobeg; k < moend; ++k)
            {
                Perm[k] = k;
                for (Sparse::Row::iterator it = A[k].Begin(); it != A[k].End(); ++it)
                {
                    if( it->first >= mobeg && it->first < moend )
                        B_Entries[cnt++] = Sparse::Row::make_entry(it->first,it->second);
                }
                B_Address[k+1] = cnt;
            }
            //DumpMatrix(B_Address,B_Entries,mobeg,moend,"NOMC64.mtx");
            treorder = Timer() - treorder;
		}
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// END MAXIMUM TRANSVERSE REORDERING /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		

		luv.reserve(nnz * 4);
		lui.reserve(nnz * 4);
		rv.reserve(nnz * 16);
		ri.reserve(nnz * 16);
		ilu.set_interval_beg(mobeg);
		ilu.set_interval_end(moend + 1);
		iu.set_interval_beg(mobeg);
		iu.set_interval_end(moend);
		ilu[mobeg] = 0;
		ir[mobeg] = 0;
#if defined(REPORT_ILU)
		std::cout << "Matrix overlap    " << mobeg << ".." << moend << std::endl;
		std::cout << "Local vector part " << vlocbeg << ".." << vlocend << std::endl;
		std::cout << "Entire vector     " << vbeg << ".." << vend << std::endl;
#endif

#if defined(RESCALE_EQUALIZE_1NORM)
		{
			trescale = Timer();
			std::fill(DL.Begin(),DL.End(),0.0);
			for (k = mobeg; k < moend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE rit = B_Address[k]; rit < B_Address[k+1]; ++rit) 
					DL[k] += fabs(B_Entries[rit].second);
				if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			}
			for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
			{
				std::fill(DR.Begin(),DR.End(),0.0);
				for (k = mobeg; k < moend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE rit = B_Address[k]; rit < B_Address[k+1]; ++rit) 
						DR[B_Entries[rit].first] += DL[k] * fabs(B_Entries[rit].second);
				}
				for (k = mobeg; k < moend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
				std::fill(DL.Begin(),DL.End(),0.0);
				for (k = mobeg; k < moend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE rit = B_Address[k]; rit < B_Address[k+1]; ++rit) 
						DL[k] += DR[B_Entries[rit].first] * fabs(B_Entries[rit].second);
				}
				for (k = mobeg; k < moend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			}

			for (k = mobeg; k < moend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE rit = B_Address[k]; rit < B_Address[k+1]; ++rit)
					B_Entries[rit].second = DL[k] * B_Entries[rit].second * DR[B_Entries[rit].first];
			}
			trescale = Timer() - trescale;
		}

		//DumpMatrix(B_Address,B_Entries,mobeg,moend,"rescale.mtx");
#elif defined(RESCALE_EQUALIZE_2NORM)
    {
			trescale = Timer();
			std::fill(DL.Begin(),DL.End(),0.0);
			for (k = mobeg; k < moend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE rit = B_Address[k]; rit < B_Address[k+1]; ++rit) 
					DL[k] += (B_Entries[rit].second)*(B_Entries[rit].second);
				if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			}
			for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
			{
				std::fill(DR.Begin(),DR.End(),0.0);
				for (k = mobeg; k < moend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE rit = B_Address[k]; rit < B_Address[k+1]; ++rit) 
						DR[B_Entries[rit].first] += DL[k] * (B_Entries[rit].second)*(B_Entries[rit].second);
				}
				for (k = mobeg; k < moend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
				std::fill(DL.Begin(),DL.End(),0.0);
				for (k = mobeg; k < moend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE rit = B_Address[k]; rit < B_Address[k+1]; ++rit) 
						DL[k] += DR[B_Entries[rit].first] * (B_Entries[rit].second)*(B_Entries[rit].second);
				}
				for (k = mobeg; k < moend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			}

			for (k = mobeg; k < moend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE rit = B_Address[k]; rit < B_Address[k+1]; ++rit)
					B_Entries[rit].second = DL[k] * (B_Entries[rit].second) * DR[B_Entries[rit].first];
			}
			trescale = Timer() - trescale;
		}
#endif


		tfactor = Timer();
		std::vector<INMOST_DATA_ENUM_TYPE> sort_indeces;
		//INMOST_DATA_ENUM_TYPE nza = 0, nzl = 0, nzu = 0, nzu2 = 0;
		//for(k = mobeg; k != moend; k++) nza += A[k].Size();
		for (k = mobeg; k != moend; k++)
		{
#if defined(REPORT_ILU_PROGRESS)
			if (k % 1000 == 0)
			{
				//std::cout << "precond: " << (double)(k-mobeg)/(double)(moend-mobeg)*100 << "\r";
				//printf("%6.2f nza %12d nzl %12d nzu %12d nzu2 %12d\r", (double)(k-mobeg)/(double)(moend-mobeg)*100,nza,nzl,nzu,nzu2);
				printf("precond: %6.2f\r", (double)(k - mobeg) / (double)(moend - mobeg) * 100);
				fflush(stdout);
			}
#endif
			//Uncompress row
			//row_uncompr
			sort_indeces.clear();
			for (r = B_Address[k]; r < B_Address[k+1]; r++) 
				if (fabs(B_Entries[r].second) > eps)
				{
					RowValues[B_Entries[r].first] = B_Entries[r].second;
#if defined(LFILL)
					RowFill[B_Entries[r].first] = 0;
#endif
					sort_indeces.push_back(B_Entries[r].first);
				}
			std::sort(sort_indeces.begin(), sort_indeces.end());
			prev = static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg)-1;
			ipred = static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg)-1;
			for (r = 0; r < sort_indeces.size(); r++)
			{
				ind = sort_indeces[r];
				RowIndeces[prev] = ind;
				prev = static_cast<INMOST_DATA_INTEGER_TYPE>(ind);
				if (ind <= k) ipred = ind;
			}
			RowIndeces[prev] = EOL;

			if (ipred != static_cast<INMOST_DATA_INTEGER_TYPE>(k))
			{
				RowValues[k] = 0.0;
#if defined(LFILL)
				RowFill[k] = 0;
#endif
				ind = RowIndeces[ipred];
				RowIndeces[ipred] = k;
				RowIndeces[k] = ind;
			}
#if defined(DIAGONAL_PERTURBATION)
			RowValues[k] = RowValues[k]*(1.0+DIAGONAL_PERTURBATION_REL) + (RowValues[k] < 0.0? -1.0 : 1.0)*DIAGONAL_PERTURBATION_REL;
#endif
			//Eliminate lower part
			//elim_lpart
			j = RowIndeces[static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg)-1];
			while (j < k) //until diagonal entry
			{
				assert(lui[iu[j]] == j);
				RowValues[j] *= luv[iu[j]]; //scale by diagonal
				leabs = fabs(RowValues[j]);
#if defined(LFILL)
				if (leabs > tau2*tau2)// introduce a non-zero, if threshold permits
#else
				if (leabs > tau2)// introduce a non-zero, if threshold permits
#endif
				{
					curr = j;
					for (r = iu[j] + 1; r < ilu[j + 1]; r++)
					{
						ind = lui[r];
						if (RowIndeces[ind] != UNDEF) //update without pondering on thresholds
						{
							RowValues[ind] -= RowValues[j] * luv[r];
#if defined(LFILL)
							RowFill[ind] = std::min(lfill[r]+1,RowFill[ind]);
#endif
						}
						else 
						{
							flin = -RowValues[j] * luv[r];
							//insert new value
							foll = curr;
							while (foll < ind)
							{
								curr = foll;
								foll = RowIndeces[curr];
							}
							assert(curr < ind);
							assert(ind < foll);
							RowIndeces[curr] = ind;
							RowIndeces[ind] = foll;
							RowValues[ind] = flin;
#if defined(LFILL)
							RowFill[ind] = lfill[r] + 1;
#endif
						}
						curr = ind;
					}

					if (leabs > tau)
					{
						curr = j;
						for (r = ir[j]; r < ir[j + 1]; r++)
						{
							//ind = U2j.GetIndex(r);
							ind = ri[r];
							if (RowIndeces[ind] != UNDEF) //update without pondering on thresholds
								RowValues[ind] -= RowValues[j] * rv[r];
							else // introduce a non-zero if threshold permits
							{
								flin = -RowValues[j] * rv[r];
								//insert new value
								foll = curr;
								while (foll < ind)
								{
									curr = foll;
									foll = RowIndeces[curr];
								}
								assert(curr < ind);
								assert(ind < foll);
								RowIndeces[curr] = ind;
								RowIndeces[ind] = foll;
								RowValues[ind] = flin;
#if defined(LFILL)
								RowFill[ind] = ENUMUNDEF;
#endif
							}
							curr = ind;
						}
					}
				}
				j = RowIndeces[j];
			}
			// Compress row
			//row_compr
			j = RowIndeces[static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg)-1];
			//find minimum value in row
			ldiag = 0;
			while (j != EOL)
			{
				INMOST_DATA_REAL_TYPE temp = fabs(RowValues[j]);
				ldiag = std::max(ldiag, temp);
				j = RowIndeces[j];
			}
			if (ldiag < tau2)
			{
				ldiag = 1.0 / tau2;
				//std::cout << "ldiag too small " << ldiag << std::endl;
			}
			else
				ldiag = 1.0 / ldiag;

			//if (ldiag > 1000) std::cout << "ldiag is big " << k << " " << ldiag << std::endl;
			//divide all entries on right from the diagonal
			j = k;
			while (j != EOL)
			{
				RowValues[j] *= ldiag;
				j = RowIndeces[j];
			}
			j = RowIndeces[static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg)-1];
			while (j < k)
			{
				mva = fabs(RowValues[j]);
				if (mva > tau2*tau2 )
				{
					if (mva > tau
#if defined(LFILL)
						|| RowFill[j] <= Lfill
#endif
						)
					{
						//L[k][j] = RowValues[j];
						lui.push_back(j); //lui indicates column index of L matrix
						luv.push_back(RowValues[j]); //luv indicates corresponding value
#if defined(LFILL)
						lfill.push_back(RowFill[j]);
#endif
						//nzl++;
					}
				}
				jn = RowIndeces[j];
				RowIndeces[j] = UNDEF;
				j = jn;
			}
			//add last diagonal entry to L matrix
			lui.push_back(k);
			luv.push_back(ldiag);
#if defined(LFILL)
			lfill.push_back(0);
#endif
			//nzl++;

			iu[k] = static_cast<INMOST_DATA_ENUM_TYPE>(luv.size()); //iu points to the first entry of current line of U matrix
			// END of L-part
			if (fabs(RowValues[j]) > tol_modif)
				udiag = 1.0 / RowValues[j];
			else
			{
				//std::cout << "udiag too small " << RowValues[j] << std::endl;
				udiag = (RowValues[j] < 0.0 ? -1.0 : 1.0) / tol_modif;
			}

			//if (fabs(udiag) > 1000) std::cout << "udiag is big " << k << " " << udiag << std::endl;

			jn = RowIndeces[j];
			RowIndeces[j] = UNDEF;
			j = jn;
			//start of U matrix entries
			//add diagonal value for U matrix
			lui.push_back(k);
			luv.push_back(udiag);
#if defined(LFILL)
			lfill.push_back(RowFill[k]);
#endif
			//nzu++;
			while (j != EOL)
			{
				mva = fabs(RowValues[j]);
				if (mva > tau2*tau2)
				{
					if (mva > tau
#if defined(LFILL)
						|| RowFill[j] <= Lfill
#endif
						)
					{
						//add values to U matrix
						lui.push_back(j);
						luv.push_back(RowValues[j]);
#if defined(LFILL)
						lfill.push_back(RowFill[j]);
#endif
						//nzu++;
					}
					else if (mva > tau2)
					{
						//add values to U2 matrix
						ri.push_back(j);
						rv.push_back(RowValues[j]);
						//nzu2++;
					}
				}
				jn = RowIndeces[j];
				RowIndeces[j] = UNDEF;
				j = jn;
			}
			ilu[k + 1] = static_cast<INMOST_DATA_ENUM_TYPE>(luv.size()); //next first entry for L
			ir[k + 1] = static_cast<INMOST_DATA_ENUM_TYPE>(rv.size()); //next first entry for U2
			//END U-part
		}
		//Rescale LU
		for (k = mobeg; k < moend; k++)
		{
			for (r = iu[k] - 1; r > ilu[k]; r--)
				luv[r - 1] /= DL[k];
			luv[iu[k] - 1] *= DL[k]; // L diagonal entry
		}
		for (k = mobeg; k < moend; k++)
		{
			for (r = iu[k] + 1; r < ilu[k + 1]; r++)
				luv[r] /= DR[lui[r]];
			luv[iu[k]] *= DR[k]; // U diagonal entry
		}

		//allocate space for reordering
		temporary.set_interval_beg(vbeg);
		temporary.set_interval_end(vend);

		tfactor = Timer() - tfactor;
		ttotal = Timer() - ttotal;
#if defined(REPORT_ILU_SUMMARY)
		INMOST_DATA_ENUM_TYPE nzu,nzl;
		nzu = 0;
		nzl = 0;
		for(INMOST_DATA_ENUM_TYPE k = mobeg; k < moend; k++)
		{
			nzl += iu[k] - ilu[k];
			nzu += ilu[k+1] - iu[k] - 1;
		}
		std::cout << "      nonzeros in A = " << nnz << std::endl;
		std::cout << "      nonzeros in L = " << nzl - (moend-mobeg) << std::endl;
		std::cout << "      nonzeros in U = " << nzu << std::endl;
		std::cout << "     nonzeros in LU = " << ilu[moend] - 1 << std::endl;
		std::cout << "     nonzeros in U2 = " << ir[moend] - 1 << std::endl;
		//std::cout << __FUNCTION__ << " done" << std::endl;
		printf("total      %f\n",ttotal);
		printf("reorder    %f (%6.2f%%)\n", treorder, 100.0*treorder/ttotal);
#if defined(REORDER_METIS_ND)
		printf("metis      graph %f nd %f\n", tmetisgraph, tmetisnd);
#endif
		printf("rescale    %f (%6.2f%%)\n", trescale, 100.0*trescale / ttotal);
		printf("factor     %f (%6.2f%%)\n", tfactor, 100.0*tfactor / ttotal);

#endif

		/*
		//partition of unity for unrestricted additive schwartz
		info.PrepareVector(div);
		std::fill(div.Begin(),div.End(),0);
		for(k = mobeg; k < moend; k++) div[k] = 1.0;
		info.Accumulate(div);
		for(k = mobeg; k < moend; k++) div[k] = 1.0/div[k];
		*/
		init = true;
		return true;
	}

  bool MTILU2_preconditioner::isInitialized(){ return init; }


  bool MTILU2_preconditioner::Finalize()
	{
		if (!isFinalized())
		{
			luv.clear();
			lui.clear();
			init = false;
		}
		return true;
	}
	bool MTILU2_preconditioner::isFinalized() { return !init; }
	MTILU2_preconditioner::~MTILU2_preconditioner()
	{
		if (!isFinalized()) Finalize();
	}
	void MTILU2_preconditioner::Copy(const Method * other)
	{
		const MTILU2_preconditioner * b = dynamic_cast<const MTILU2_preconditioner *>(other);
		assert(b != NULL);
		info = b->info;
		Alink = b->Alink;
		nnz = b->nnz;
		luv = b->luv;
		lui = b->lui;
		iu = b->iu;
		ilu = b->ilu;
		Perm = b->Perm;
	}


  MTILU2_preconditioner::MTILU2_preconditioner(const MTILU2_preconditioner & other)
		:Method(other)
	{
		Copy(&other);
	}
	MTILU2_preconditioner & MTILU2_preconditioner::operator =(MTILU2_preconditioner const & other)
	{
		Copy(&other);
		return *this;
	}
	bool MTILU2_preconditioner::Solve(Sparse::Vector & input, Sparse::Vector & output)
	{
		assert(isInitialized());
#if defined(USE_OMP)
#pragma omp single
#endif
		{
			INMOST_DATA_ENUM_TYPE mobeg, moend, r, k, vbeg,vend; //, end;
			info->GetOverlapRegion(info->GetRank(),mobeg,moend);
			info->GetVectorRegion(vbeg,vend);
      for(k = vbeg; k < mobeg; k++) temporary[k] = 0;
			for(k = mobeg; k < moend; k++) temporary[k] = input[k];
      for(k = moend; k < vend; k++) temporary[k] = 0;

      //for(k = vbeg; k < vend; k++) temporary[k] = input[k];
      
      //for(k = vbeg; k < vend; k++) temporary[k] = input[Perm[k]];
      //for(k = vbeg; k < vend; k++) temporary[Perm[k]] = input[k];
			for(k = mobeg; k < moend; k++) //iterate over L part
			{
				for(r = iu[k]-1; r > ilu[k]; r--) 
					temporary[k] -= luv[r-1]*temporary[lui[r-1]];
				temporary[k] *= luv[iu[k]-1];
			}
			for(k = moend; k > mobeg; k--) //iterate over U part
			{
				for(r = iu[k-1]+1; r < ilu[k]; r++)
					temporary[k-1] -= luv[r]*temporary[lui[r]];
				temporary[k-1] *= luv[iu[k-1]];
			}
			
      //for (k = mobeg; k < moend; ++k) output[k] = temporary[k];
			for (k = mobeg; k < moend; ++k) if( Perm[k] != ENUMUNDEF ) output[k] = temporary[Perm[k]];
			//for (k = mobeg; k < moend; ++k) output[Perm[k]] = temporary[k];
      for(k = vbeg; k < mobeg; k++) output[k] = 0;
			for(k = moend; k < vend; k++) output[k] = 0;
		}
		info->Accumulate(output);
		return true;
	}
#endif
