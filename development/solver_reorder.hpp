#pragma once
#ifndef __SOLVER_REORDER__
#define __SOLVER_REORDER__

//TODO:
// finish implementation, test
// use parallel merge-sort instead of final std::sort

#include "inmost_solver.h"
#include "solver_prototypes.hpp"

using namespace INMOST;

class ddPQ_reorder : public Method // Two-side unsymmetric reordering
{
	MPI_Datatype wgt_pair_type;
	typedef std::pair<INMOST_DATA_REAL_TYPE,INMOST_DATA_ENUM_TYPE> wgt_pair;
	
	MPI_Datatype col_type;
	typedef struct col_weight_struct
	{
		INMOST_DATA_REAL_TYPE max, sum;
		INMOST_DATA_ENUM_TYPE nnz, pos;
	} col_weight;
	
	/*
	MPI_Op col_op;
	static void col_op_func(void * pa, void * pb, int * len, MPI_Datatype * datatype)
	{
		col_weight * a = (col_weight *)pa;
		col_weight * b = (col_weight *)pb;
		for(int i = 0; i < *len; i++)
		{
			b[i].sum += a[i].sum;
			b[i].nnz += a[i].nnz;
			if( a[i].max > b[i].max ) b[i].max = a[i].max;
		}
	}
	*/
	Solver::Matrix * Alink;
	Solver::Vector * xlink;
	Solver::Vector * blink;
	Solver::OrderInfo * info;
	interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> p,q;
	INMOST_DATA_REAL_TYPE min;
public:
	INMOST_DATA_REAL_TYPE & RealParameter(std::string name) 
	{
		if( name == "minimum" )
			return min;
		throw -1;
	}
	INMOST_DATA_ENUM_TYPE & EnumParameter(std::string name) {throw -1;}
	ddPQ_reorder(Solver::Matrix & A, Solver::Vector & x, Solver::Vector & b,  Solver::OrderInfo & Minfo)
		:Alink(&A),xlink(&x),blink(&b),info(&Minfo)
	{
		min = -1e54;
		
		{
			MPI_Datatype type[5] = 
			{ 
				INMOST_MPI_DATA_REAL_TYPE, 
				INMOST_MPI_DATA_REAL_TYPE, 
				INMOST_MPI_DATA_ENUM_TYPE,
				INMOST_MPI_DATA_ENUM_TYPE, 
				MPI_UB
			};
			int blocklen[5] = { 1, 1, 1, 1, 1};
			MPI_Aint disp[5];
			disp[0] = offsetof(col_weight,max);
			disp[1] = offsetof(col_weight,sum);
			disp[2] = offsetof(col_weight,nnz);
			disp[3] = offsetof(col_weight,pos);
			disp[4] = sizeof(col_weight);
			MPI_Type_create_struct(5, blocklen, disp, type, &col_type);
			MPI_Type_commit(&col_type);
		}
		
		{
			MPI_Datatype type[3] = 
			{
				INMOST_MPI_DATA_REAL_TYPE, 
				INMOST_MPI_DATA_ENUM_TYPE, 
				MPI_UB
			};
			int blocklen[3] = { 1, 1, 1};
			MPI_Aint disp[5];
			disp[0] = offsetof(wgt_pair,first);
			disp[1] = offsetof(wgt_pair,second);
			disp[2] = sizeof(wgt_pair);
			MPI_Type_create_struct(3, blocklen, disp, type, &wgt_pair_type);
			MPI_Type_commit(&wgt_pair_type);
		}

		//MPI_Op_create((MPI_User_function *)col_op_func,1,&col_op);
	}
	void Copy(const Method * other)
	{
		const ddPQ_reorder * b = dynamic_cast<const ddPQ_reorder *>(other);
		assert(b != NULL);
		Alink = b->Alink;
		xlink = b->xlink;
		blink = b->blink;
		info = b->info;
		p = b->p;
		q = b->q;
	}
	Method * Duplicate() { return new ddPQ_reorder(*this); }
	ddPQ_reorder(const ddPQ_reorder & other)
	{
		Copy(&other);
	}
	ddPQ_reorder & operator = (ddPQ_reorder const & other)
	{
		Copy(&other);
		return *this;
	}
	~ddPQ_reorder()
	{
		//MPI_Op_free(&col_op);
		MPI_Type_free(&col_type);
		MPI_Type_free(&wgt_pair_type);
	}
	bool Initialize()
	{
		if (isInitialized()) Finalize();
		//permute matrix, x and right-hand-side
		int ierr;
		INMOST_DATA_ENUM_TYPE mbeg, mend, mzero,msize, k,i,j,l, vbeg, vend;
		std::vector<INMOST_DATA_ENUM_TYPE> sizes(info->GetSize());
		std::vector<int> recvcount(info->GetSize()), displs(info->GetSize());
		info->GetLocalRegion(0,mzero,mend);
		info->GetLocalRegion(info->GetSize()-1,mbeg,msize);
		info->GetLocalRegion(info->GetRank(),mbeg,mend);
		interval<INMOST_DATA_ENUM_TYPE,wgt_pair> col_wgts(mzero,msize), row_wgts(mzero,msize);
		{
			INMOST_DATA_REAL_TYPE maxval = min;
			for(INMOST_DATA_ENUM_TYPE k = 0; k < info->GetSize(); k++)
			{
				info->GetLocalRegion(k,vbeg,vend);
				displs[k] = vbeg;
				recvcount[k] = vend-vbeg;
			}
			info->GetLocalRegion(info->GetRank(),mbeg,mend);
			info->GetVectorRegion(vbeg,vend);
			interval<INMOST_DATA_ENUM_TYPE,wgt_pair> loc_wgt(mbeg,mend);
			interval<INMOST_DATA_ENUM_TYPE,col_weight> loc_col_val(vbeg,vend);
			col_weight init;
			init.max = min;
			init.sum = 0;
			init.nnz = 0;
			init.pos = ENUMUNDEF;
			std::vector<col_weight> send_storage(info->GetSendExchangeSize()),recv_storage(info->GetRecvExchangeSize());
			std::fill(loc_col_val.begin(),loc_col_val.end(),init);
			
			for(k = mbeg; k < mend; k++)
			{
				INMOST_DATA_ENUM_TYPE maxpos = ENUMUNDEF;
				INMOST_DATA_REAL_TYPE rowmax = min, rowsum = 0;
				for(Solver::Row::iterator q = (*Alink)[k].Begin(); q != (*Alink)[k].End(); ++q)
				{
					if( q->second > maxval )
					{
						INMOST_DATA_ENUM_TYPE j = q->first;
						INMOST_DATA_REAL_TYPE vabs = fabs(q->second);
						rowsum += vabs;
						if( rowmax < vabs ) 
						{
							rowmax = vabs;
							maxpos = j;
						}
						
						loc_col_val[j].nnz++;
						loc_col_val[j].sum += vabs;
						if( loc_col_val[j].max < vabs ) 
						{
							loc_col_val[j].max = vabs;
							loc_col_val[j].pos = k;
						}
						
					}
				}
				loc_wgt[k].first = rowmax/rowsum/static_cast<INMOST_DATA_REAL_TYPE>((*Alink)[k].Size());
				loc_wgt[k].second = maxpos;
			}
			//can use merge-sort algorithm here
			MPI_Allgatherv(&loc_wgt[mbeg],mend-mbeg,wgt_pair_type,&row_wgts[mzero],&recvcount[0],&displs[0],wgt_pair_type,info->GetComm());
			
			j = 1, l = 0;
			for(i = 0; i < info->GetSendExchangeArray()[0]; i++)
			{
				MPI_Irecv(&send_storage[l],info->GetSendExchangeArray()[j+1],col_type,info->GetSendExchangeArray()[j],info->GetSendExchangeArray()[j]*info->GetSize()+info->GetRank(),info->GetComm(),&info->GetSendRequests()[i]);
				l += info->GetSendExchangeArray()[j+1];
				j += info->GetSendExchangeArray()[j+1] + 2;
			}
			j = 1, l = 0;
			for(i = 0; i < info->GetRecvExchangeArray()[0]; i++)
			{
				for(k = 0; k < info->GetRecvExchangeArray()[j+1]; k++)
					recv_storage[l+k] = loc_col_val[info->GetRecvExchangeArray()[k+j+2]];
				MPI_Isend(&recv_storage[l],info->GetRecvExchangeArray()[j+1],col_type,info->GetRecvExchangeArray()[j],info->GetRank()*info->GetSize()+info->GetRecvExchangeArray()[j],info->GetComm(),&info->GetRecvRequests()[i]);
				l += info->GetRecvExchangeArray()[j+1];
				j += info->GetRecvExchangeArray()[j+1] + 2;
			}
			if( info->GetSendExchangeArray()[0] > 0 )
			{
				MPI_Waitall(info->GetSendRequestsSize(),info->GetSendRequests(),MPI_STATUSES_IGNORE);
				j = 1, l = 0;
				for(i = 0; i < info->GetSendExchangeArray()[0]; i++)
				{
					for(k = 0; k < info->GetSendExchangeArray()[j+1]; k++)
					{
						loc_col_val[info->GetSendExchangeArray()[k+j+2]].sum += send_storage[l+k].sum;
						loc_col_val[info->GetSendExchangeArray()[k+j+2]].nnz += send_storage[l+k].nnz;
						if( loc_col_val[info->GetSendExchangeArray()[k+j+2]].max < send_storage[l+k].max )
						{
							loc_col_val[info->GetSendExchangeArray()[k+j+2]].max = send_storage[l+k].max;
							loc_col_val[info->GetSendExchangeArray()[k+j+2]].pos = send_storage[l+k].pos;
						}
					}
					l += info->GetSendExchangeArray()[j+1];
					j += info->GetSendExchangeArray()[j+1] + 2;
				}
			}
			if( info->GetRecvExchangeArray()[0] > 0 ) 
				MPI_Waitall(info->GetRecvRequestsSize(),info->GetRecvRequests(),MPI_STATUSES_IGNORE);
			for(k = mbeg; k < mend; k++)
			{
				loc_wgt[k].first = loc_col_val[k].max/loc_col_val[k].sum/static_cast<INMOST_DATA_REAL_TYPE>(loc_col_val[k].nnz);
				loc_wgt[k].second = loc_col_val[k].pos;
			}
			
			//can use merge-sort algorithm here
			MPI_Allgatherv(&loc_wgt[mbeg],mend-mbeg,wgt_pair_type,&col_wgts[mzero],&recvcount[0],&displs[0],wgt_pair_type,info->GetComm());
		}
		{
			std::vector< std::pair<INMOST_DATA_REAL_TYPE, std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE>> > sort_wgts(2*(msize-mzero));
			for(k = mzero; k < msize; k++)
			{
				sort_wgts[2*(k-mzero)].second.first = k;
				sort_wgts[2*(k-mzero)].second.second = row_wgts[k].second;
				sort_wgts[2*(k-mzero)].first = row_wgts[k].first;
				sort_wgts[2*(k-mzero)+1].second.first = col_wgts[k].second;
				sort_wgts[2*(k-mzero)+1].second.second = k;
				sort_wgts[2*(k-mzero)+1].first = col_wgts[k].first;
			}
			//may be replaced by merge-sort
			std::sort(sort_wgts.begin(),sort_wgts.end());
			p.set_interval_beg(mzero);
			p.set_interval_end(msize);
			q.set_interval_beg(mzero);
			q.set_interval_end(msize);
			std::fill(p.begin(),p.end(),ENUMUNDEF);
			std::fill(q.begin(),q.end(),ENUMUNDEF);
			l = mzero;
			for(k = 0; k < static_cast<INMOST_DATA_ENUM_TYPE>(sort_wgts.size()); )
			{
				for(i = 0; i < info->GetSize(); i++) if( recvcount[i] > 0 )
				{
					if( p[sort_wgts[k].second.first] == ENUMUNDEF && q[sort_wgts[k].second.second] == ENUMUNDEF )
					{
						p[sort_wgts[k].second.first] = displs[i];
						q[sort_wgts[k].second.second] = displs[i];
						displs[i]++;
						recvcount[i]--;
						k++;
						l++;
					}
				}
			}
			if( l < msize )
			{
				std::vector<INMOST_DATA_ENUM_TYPE> indicesp,indicesq;
				for(i = 0; i < info->GetSize(); i++)
				{
					while( recvcount[i] > 0 ) 
					{
						indicesp.push_back(displs[i]);
						displs[i]++;
						recvcount[i]--;
					}
				}
				indicesq = indicesp;
				for(k = mzero; k < msize; k++)
				{
					if( p[k] == ENUMUNDEF )
					{
						p[k] = indicesp.back();
						indicesp.pop_back();
					}
					if( q[k] == ENUMUNDEF )
					{
						q[k] = indicesq.back();
						indicesq.pop_back();
					}
				}
			}
		}
		//to whome do I send
		// who is sending to me

		//permute matrix
		//permute x according to Q
		//permute rhs according to P
	}
	bool ReplaceMAT(Solver::Matrix & A) { if (isInitialized()) Finalize();  Alink = &A; return true; }
	bool ReplaceSOL(Solver::Vector & x)
	{
		//permute sol according to Q
	}
	bool ReplaceRHS(Solver::Vector & b)
	{
		//permute rhs according to P
	}
	bool Finalize()
	{
		//permute everything back
	}
	bool isInitialized() { return true; }
	bool isFinalized() { return true; }
	bool Solve(Solver::Vector & x, Solver::Vector & out)
	{
		//compute back permutation from x to out
	}
};

#endif