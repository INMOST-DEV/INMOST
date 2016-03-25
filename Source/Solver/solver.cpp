#define _CRT_SECURE_NO_WARNINGS
#include "inmost_solver.h"
#if defined(USE_SOLVER)
#include "solver_petsc.h"
#include "solver_superlu.h"
#include "solver_ani.h"
#include "solver_ilu2.hpp"
#include "solver_ddpqiluc2.hpp"
#if defined(HAVE_SOLVER_MPTILUC2)
#include "solver_mtiluc2.hpp"
#endif
#if defined(HAVE_SOLVER_MPTILU2)
#include "solver_mtilu2.hpp"
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
#include "solver_fcbiilu2.h"
#endif
#if defined(HAVE_SOLVER_K3BIILU2)
#include "solver_k3biilu2.h"
#endif
#include "solver_bcgsl.hpp"
#include <fstream>
#include <sstream>
#include <stddef.h>

#if defined(USE_SOLVER_TRILINOS)
#if defined(USE_MPI)
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "AztecOO.h"
#include "Ifpack.h"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "BelosEpetraOperator.h"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif
//#define USE_OMP

//#define KSOLVER BCGSL_solver
#define KSOLVER BCGS_solver
//#define ACCELERATED_CONDEST
//#define PRINT_CONDEST

namespace INMOST
{
	static std::string petsc_database_file = "";
	static std::string trilinos_aztec_database_file = "";
	static std::string trilinos_ifpack_database_file = "";
	static std::string trilinos_belos_database_file = "";
	static std::string trilinos_ml_database_file = "";
	static std::string ani_database_file = "";
	static std::string fcbiilu2_database_file = "";
	static std::string k3biilu2_database_file = "";
	static std::string superlu_database_file = "";



#define GUARD_MPI(x) {ierr = x; if( ierr != MPI_SUCCESS ) {char str[4096]; int len; MPI_Error_string(ierr,str,&len); std::cout << #x << " not successfull: " << str << std::endl; MPI_Abort(comm,-1000);}}
#define HASH_TABLE_SIZE 2048

	bool Solver::is_initialized = false;
	bool Solver::is_finalized = false;

  std::string Solver::TypeName(Type t)
  {
    switch(t)
    {
    case INNER_ILU2: return "INNER_ILU2";
    case INNER_DDPQILUC: return "INNER_DDPQILUC";
    case INNER_MPTILUC: return "INNER_MPTILUC";
    case INNER_MPTILU2: return "INNER_MPTILU2";
    case Trilinos_Aztec: return "Trilinos_Aztec";
    case Trilinos_Belos: return "Trilinos_Belos";
    case Trilinos_ML: return "Trilinos_ML";
    case Trilinos_Ifpack: return "Trilinos_Ifpack";
    case PETSc: return "PETSc";
    case ANI: return "ANI";
    case FCBIILU2: return "FCBIILU2";
    case K3BIILU2: return "K3BIILU2";
	case SUPERLU: return "SUPERLU";
    }
    return "Unknown";
  }

	int comparator(const void * pa, const void *pb)
	{
		INMOST_DATA_ENUM_TYPE * a = (INMOST_DATA_ENUM_TYPE *)pa, * b = (INMOST_DATA_ENUM_TYPE *)pb;
		return a[0] - b[0];
	}

	INMOST_DATA_ENUM_TYPE binary_search_pairs(INMOST_DATA_ENUM_TYPE * link, INMOST_DATA_ENUM_TYPE size, INMOST_DATA_ENUM_TYPE find)
	{
		INMOST_DATA_ENUM_TYPE rcur = size >> 1, lcur = 0, mid, chk;
		while( rcur >= lcur )
		{
			mid = lcur + ((rcur-lcur) >> 1);
			chk = mid << 1;
			if( link[chk] < find ) lcur = mid+1;
			else if( link[chk] > find ) rcur = mid-1;
			else return chk;
		}
		return size;
	}
	void Solver::OrderInfo::Integrate(INMOST_DATA_REAL_TYPE * inout, INMOST_DATA_ENUM_TYPE num) const
	{
#if defined(USE_MPI)
		if( GetSize() == 1 ) return;
#if defined(USE_OMP)
#pragma omp single
#endif
		{
			int ierr = 0;
			dynarray<INMOST_DATA_REAL_TYPE,1024> temp(num);
			memcpy(temp.data(),inout,sizeof(INMOST_DATA_REAL_TYPE)*num);
			GUARD_MPI(MPI_Allreduce(temp.data(),inout,num,INMOST_MPI_DATA_REAL_TYPE,MPI_SUM,comm));
		}
#else
		(void) inout;
		(void) num;
#endif
	}

	
	void Solver::OrderInfo::PrepareMatrix(Sparse::Matrix & m, INMOST_DATA_ENUM_TYPE overlap)
	{
		have_matrix = true;
		m.isParallel() = true;
		INMOST_DATA_ENUM_TYPE  two[2];
		INMOST_DATA_ENUM_TYPE mbeg,mend;
		int initial_rank;
#if defined(USE_MPI)
		int ierr = 0;
		if( comm != INMOST_MPI_COMM_WORLD )
		{
			MPI_Comm_free(&comm);
			comm = INMOST_MPI_COMM_WORLD;
		}
		if( m.GetCommunicator() == INMOST_MPI_COMM_WORLD )
			comm = INMOST_MPI_COMM_WORLD;
		else MPI_Comm_dup(m.GetCommunicator(), &comm);
		MPI_Comm_rank(comm,&rank);
		MPI_Comm_size(comm,&size);
#else
		(void) overlap;
		rank = 0;
		size = 1;
#endif
		initial_rank = rank;
		//std::vector<MPI_Request> requests;
		global_overlap.resize(size*2);
		global_to_proc.resize(size+1);
		m.GetInterval(mbeg,mend);
		global_to_proc[0] = 0;
		initial_matrix_begin = mbeg;
		initial_matrix_end = mend;
		two[0] = mbeg;
		two[1] = mend;
#if defined(USE_MPI)
		GUARD_MPI(MPI_Allgather(two,2,INMOST_MPI_DATA_ENUM_TYPE,&global_overlap[0],2,INMOST_MPI_DATA_ENUM_TYPE,comm));
#else
		local_vector_begin = initial_matrix_begin = local_matrix_begin = global_overlap[0] = mbeg;
		local_vector_end   = initial_matrix_end   = local_matrix_end   = global_overlap[1] = mend;
		global_to_proc[1] = mend;
#endif
#if defined(USE_MPI)
		//reorder processors if needed
		{
			//starts of local indexes should appear in asscending order
			bool reorder = false;
			for(int k = 0; k < size-1; k++)
				if( global_overlap[2*k] > global_overlap[2*(k+1)] )
				{
					reorder = true;
					break;
				}
			if( reorder )
			{
				storage_type temp(size*2);
				//assemble array that includes rank
				for(int k = 0; k < size; ++k)
				{
					temp[2*k+0] = global_overlap[2*k];
					temp[2*k+1] = k;
				}
				//sort array
				qsort(&temp[0],size,sizeof(INMOST_DATA_ENUM_TYPE)*2,comparator);
				//create new group
				MPI_Group oldg, newg;
				MPI_Comm newcomm;
				std::vector<int> ranks(size);
				for(int k = 0; k < size; ++k)
					ranks[k] = temp[2*k+1];
				GUARD_MPI(MPI_Comm_group(comm,&oldg));
				GUARD_MPI(MPI_Group_incl(oldg,size,&ranks[0],&newg));
				GUARD_MPI(MPI_Comm_create(comm,newg,&newcomm));
				if( comm != INMOST_MPI_COMM_WORLD )
				{
					GUARD_MPI(MPI_Comm_free(&comm));
				}
				comm = newcomm;
				//compute new rank
				MPI_Comm_rank(comm,&rank);
				//sort array pairs, so we don't need to exchange them again
				qsort(&global_overlap[0],size,sizeof(INMOST_DATA_ENUM_TYPE)*2,comparator);
			}
			//now check that there are no overlaps of local indexes
			//every mend must be equal to every mbeg
			reorder = false;
			for(int k = 0; k < size-1; k++)
				if( global_overlap[2*k+1] != global_overlap[2*(k+1)] )
				{
					//check that end is strictly greater then begin
					if( global_overlap[2*k+1] < global_overlap[2*(k+1)] )
					{
						if( initial_rank == 0 )
						{
							std::cout << __FILE__ << ":" << __LINE__ << " Matrix index intervals are not complete:";
							std::cout << " processor " << k+0 << " interval " << global_overlap[2*(k+0)] << ":" << global_overlap[2*(k+0)+1];
							std::cout << " processor " << k+1 << " interval " << global_overlap[2*(k+1)] << ":" << global_overlap[2*(k+1)+1];
							std::cout << std::endl;
							MPI_Abort(comm,-1000);
						}
					}
					reorder = true;
				}
			if( reorder )
			{
				storage_type old_overlap(global_overlap);
				//move local bounds to get non-overlapping regions
				for(int k = 0; k < size-1; k++)
					while( global_overlap[2*k+1] > global_overlap[2*(k+1)] )
					{
						//move bounds to equalize sizes
						if( global_overlap[2*k+1] - global_overlap[2*k] < global_overlap[2*(k+1)+1] - global_overlap[2*(k+1)] )
							global_overlap[2*k+1]--; //move right bound of the current processor to left
						else
							global_overlap[2*(k+1)]++; //move left bound of the next processor to right
					}
				
				//TODO: if we need to merge overlapping parts of the matrices - do it here
			}
			local_matrix_begin = global_overlap[2*rank+0];
			local_matrix_end   = global_overlap[2*rank+1];
			for(int k = 0; k < size; k++)
				global_to_proc[k+1] = global_overlap[2*k+1];
		}
		MPI_Status stat;
		INMOST_DATA_ENUM_TYPE ext_pos = local_matrix_end;
		//may replace std::map here
		small_hash<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE,HASH_TABLE_SIZE> global_to_local;
		std::vector< std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> > current_global_to_local;
		std::vector< Sparse::Row::entry > send_row_data, recv_row_data;
		std::vector< INMOST_DATA_ENUM_TYPE > send_row_sizes, recv_row_sizes;
		std::vector<INMOST_DATA_ENUM_TYPE> incoming(4*size);
		std::vector<MPI_Request> requests;
		INMOST_DATA_ENUM_TYPE total_send = 0, total_recv = 0;
		INMOST_DATA_ENUM_TYPE local_start = local_matrix_begin, local_end = local_matrix_end;
		for(INMOST_DATA_ENUM_TYPE it = 0; it < overlap+1; it++)
		{
			total_send = 0, total_recv = 0;
			current_global_to_local.clear();
			for(INMOST_DATA_ENUM_TYPE k = local_start; k < local_end; ++k)
			{
				Sparse::Row & r = m[k];
				INMOST_DATA_ENUM_TYPE jend = r.Size(), ind;
				for(INMOST_DATA_ENUM_TYPE j = 0; j < jend; ++j)
				{
					ind = r.GetIndex(j);
					if( ind < local_matrix_begin || ind >= local_matrix_end) 
					{
						INMOST_DATA_ENUM_TYPE & recv_pos = global_to_local[ind];
						if( recv_pos == 0 ) //this number was not assigned yet
						{
							recv_pos = ext_pos++;
							if( it < overlap ) current_global_to_local.push_back(std::make_pair(ind,recv_pos));
						}
					}
				}
			}
			if( it == overlap ) 
				current_global_to_local = global_to_local.serialize();
			std::sort(current_global_to_local.begin(),current_global_to_local.end());
			//if( !current_global_to_local.empty() )
			{
				//check all the indexes that comes from other processors
				//for every processor we need arrays:
				// processor -> (array of index positions where to receive))
				// processor -> (array of index positions from where to send)
				memset(&incoming[0],0,sizeof(INMOST_DATA_ENUM_TYPE)*size*2);
				vector_exchange_recv.clear();
				vector_exchange_recv.push_back(0);
				if( !current_global_to_local.empty() )
				{
					INMOST_DATA_ENUM_TYPE proc_beg = GetProcessor(current_global_to_local.begin()->first), proc_end = GetProcessor(current_global_to_local.rbegin()->first)+1;
					INMOST_DATA_ENUM_TYPE current_ind = 0;
					for(INMOST_DATA_ENUM_TYPE proc = proc_beg; proc < proc_end; proc++)
					{
						bool first = true;
						INMOST_DATA_ENUM_TYPE numind = static_cast<INMOST_DATA_ENUM_TYPE>(vector_exchange_recv.size() + 1);
						while( current_ind < current_global_to_local.size() && current_global_to_local[current_ind].first < global_to_proc[proc+1] )
						{
							INMOST_DATA_ENUM_TYPE k = current_global_to_local[current_ind].first;
							if( first )
							{
								vector_exchange_recv.push_back(proc);
								vector_exchange_recv.push_back(1);
								vector_exchange_recv.push_back(k);
								first = false;
							}
							else
							{
								vector_exchange_recv[numind]++;
								vector_exchange_recv.push_back(k);
							}
							current_ind++;
						}
						if( !first ) 
						{
							incoming[proc]++;
							incoming[proc+size] += vector_exchange_recv[numind];
							vector_exchange_recv[0]++;
						}
					}
				}
				
				GUARD_MPI(MPI_Allreduce(&incoming[0],&incoming[2*size],size*2,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,comm));
				//std::cout << GetRank() << " MPI_Allreduce " << __FILE__ << ":" << __LINE__ << " incoming " << incoming[size*2+rank] << " size " << incoming[size*3+rank] << std::endl;
				//prepare array that helps exchanging vector values
				requests.resize(2*vector_exchange_recv[0] + incoming[size*2+rank]);
				INMOST_DATA_ENUM_TYPE j = 1;
				for(INMOST_DATA_ENUM_TYPE k = 0; k < vector_exchange_recv[0]; k++) //send rows that i want to receive
				{
					total_recv += vector_exchange_recv[j+1];
					GUARD_MPI(MPI_Isend(&vector_exchange_recv[j+1],1,INMOST_MPI_DATA_ENUM_TYPE,vector_exchange_recv[j],size+vector_exchange_recv[j],comm,&requests[k])); //send number of rows
					GUARD_MPI(MPI_Isend(&vector_exchange_recv[j+2],vector_exchange_recv[j+1],INMOST_MPI_DATA_ENUM_TYPE,vector_exchange_recv[j],2*size+vector_exchange_recv[j],comm,&requests[k+vector_exchange_recv[0]])); //send row positions
					j += vector_exchange_recv[j+1] + 2;
				}

				recv_row_sizes.resize(incoming[size*3+rank]);
				vector_exchange_send.resize(1+incoming[size*2+rank]*2+incoming[size*3+rank]);
				vector_exchange_send[0] = 0;
				j = 1;
				for(INMOST_DATA_ENUM_TYPE k = 0; k < incoming[size*2+rank]; k++) //receive rows that others want from me
				{
					INMOST_DATA_ENUM_TYPE msgsize;
					GUARD_MPI(MPI_Recv(&msgsize,1,INMOST_MPI_DATA_ENUM_TYPE,MPI_ANY_SOURCE,size+rank,comm,&stat)); //recv number of rows
					vector_exchange_send[j++] = stat.MPI_SOURCE;
					vector_exchange_send[j++] = msgsize;
					//std::cout << GetRank() << " MPI_Irecv size " << msgsize << " rank " << stat.MPI_SOURCE << " tag " << 2*size+rank << __FILE__ << ":" << __LINE__ << std::endl;
					GUARD_MPI(MPI_Irecv(&vector_exchange_send[j],msgsize,INMOST_MPI_DATA_ENUM_TYPE,stat.MPI_SOURCE,2*size+rank,comm,&requests[2*vector_exchange_recv[0]+k])); //recv rows
					j += msgsize;
					total_send += msgsize;
					vector_exchange_send[0]++;
				}
				assert(total_send == incoming[size*3+rank]);
				assert(vector_exchange_send[0] == incoming[size*2+rank]);
				if( 2*vector_exchange_recv[0] + incoming[size*2+rank] > 0 )
					GUARD_MPI(MPI_Waitall(2*vector_exchange_recv[0] + incoming[size*2+rank],&requests[0],MPI_STATUSES_IGNORE));
			}
			/*
			else
			{
				vector_exchange_recv.resize(1,0);
				vector_exchange_send.resize(1,0);
			}
			*/
			if( it == overlap )
			{
				//std::cout << rank << " reorder " << std::endl;
				//now we need to reorder off-diagonal parts of the matrix
				for(INMOST_DATA_ENUM_TYPE k = local_matrix_begin; k < local_end; ++k)
					for(Sparse::Row::iterator jt = m[k].Begin(); jt != m[k].End(); ++jt)
						if( global_to_local.is_present(jt->first) ) 
							jt->first = global_to_local[jt->first];
						else
						{
							assert(jt->first >= local_matrix_begin);
							assert(jt->first < local_matrix_end);
						}
				local_vector_begin = local_matrix_begin;
				local_vector_end = ext_pos;
				{ 
					// change indexes for recv array
					INMOST_DATA_ENUM_TYPE i,j = 1,k;
					//for(k = 0; k < GetRank(); k++) MPI_Barrier(comm);
					//std::cout << "rank " << GetRank() << std::endl;
					//std::cout << "recv:" << std::endl;
					for(i = 0; i < vector_exchange_recv[0]; i++)
					{
						//std::cout << "proc " << vector_exchange_recv[j] << " size " << vector_exchange_recv[j+1] << std::endl;
						j++; //skip processor number
						for(k = 0; k < vector_exchange_recv[j]; ++k)
						{
							assert(global_to_local.is_present(vector_exchange_recv[j+k+1]));
							vector_exchange_recv[j+k+1] = global_to_local[vector_exchange_recv[j+k+1]];
							assert(vector_exchange_recv[j+k+1] >= local_matrix_end);
						}
						j+=vector_exchange_recv[j]+1; //add vector length + size position
					}
					//check that indexes in send array are in local matrix bounds
					//std::cout << "send:" << std::endl;
#ifndef NDEBUG
					j = 1;
					for(i = 0; i < vector_exchange_send[0]; i++)
					{
						//std::cout << "proc " << vector_exchange_send[j] << " size " << vector_exchange_send[j+1] << std::endl;
						j++; //skip processor number
						for(k = 0; k < vector_exchange_send[j]; ++k)
						{
							assert(vector_exchange_send[j+k+1] >= local_matrix_begin);
							assert(vector_exchange_send[j+k+1] < local_matrix_end);
						}
						j+=vector_exchange_send[j]+1; //add vector length + size position
					}
#endif
					//for(k = GetRank(); k < GetSize(); k++) MPI_Barrier(comm);
				}
				//prepare array local->global
				extended_indexes.resize(local_vector_end-local_matrix_end);
				for(std::vector< std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> >::iterator jt = current_global_to_local.begin(); jt != current_global_to_local.end(); ++jt)
					extended_indexes[jt->second-local_matrix_end] = jt->first;

				send_storage.resize(total_send);
				recv_storage.resize(total_recv);
				send_requests.resize(vector_exchange_send[0]);
				recv_requests.resize(vector_exchange_recv[0]);
				
			}
			else
			{
				send_row_sizes.resize(total_send);
				recv_row_sizes.resize(total_recv);

				INMOST_DATA_ENUM_TYPE j = 1, q = 0, f = 0, total_rows_send = 0, total_rows_recv = 0;
				for(INMOST_DATA_ENUM_TYPE k = 0; k < vector_exchange_recv[0]; k++) //recv sizes of rows
				{
					GUARD_MPI(MPI_Irecv(&recv_row_sizes[q],vector_exchange_recv[j+1],INMOST_MPI_DATA_ENUM_TYPE,vector_exchange_recv[j],3*size+vector_exchange_recv[j],comm,&requests[k]));
					q += vector_exchange_recv[j+1];
					j += vector_exchange_recv[j+1]+2;
				}
				j = 1;
				q = 0;

				for(INMOST_DATA_ENUM_TYPE k = 0; k < vector_exchange_send[0]; k++) //send sizes of rows
				{
					for(INMOST_DATA_ENUM_TYPE r = 0; r < vector_exchange_send[j+1]; r++)
					{
						send_row_sizes[q+r] = m[vector_exchange_send[j+2+r]].Size();
						total_rows_send += m[vector_exchange_send[j+2+r]].Size();
					}
					GUARD_MPI(MPI_Isend(&send_row_sizes[q],vector_exchange_send[j+1],INMOST_MPI_DATA_ENUM_TYPE,vector_exchange_send[j],3*size+rank,comm,&requests[vector_exchange_recv[0]+k])); //recv rows
					//remember processor numbers here
					q += vector_exchange_send[j+1];
					j += vector_exchange_send[j+1]+2;
				}
				send_row_data.clear();
				send_row_data.reserve(total_rows_send);

				
				j = 1;
				for(INMOST_DATA_ENUM_TYPE k = 0; k < vector_exchange_send[0]; k++) //accumulate data in array
				{
					for(INMOST_DATA_ENUM_TYPE r = 0; r < vector_exchange_send[j+1]; r++)
						send_row_data.insert(send_row_data.end(),m[vector_exchange_send[j+2+r]].Begin(),m[vector_exchange_send[j+2+r]].End());
					j += vector_exchange_send[j+1]+2;
				}

				
				//replace by mpi_waitsome
				if( vector_exchange_recv[0]+vector_exchange_send[0] > 0 )
					GUARD_MPI(MPI_Waitall(vector_exchange_recv[0]+vector_exchange_send[0],&requests[0],MPI_STATUSES_IGNORE));

				
				j = 1;
				q = 0;
				for(INMOST_DATA_ENUM_TYPE k = 0; k < vector_exchange_recv[0]; k++) //compute total size of data to receive
				{
					for(INMOST_DATA_ENUM_TYPE r = 0; r < vector_exchange_recv[j+1]; r++)
						total_rows_recv += recv_row_sizes[q+r];
					q += vector_exchange_recv[j+1];
					j += vector_exchange_recv[j+1]+2;
				}
				recv_row_data.resize(total_rows_recv);
				j = 1;
				q = 0;
				f = 0;
				for(INMOST_DATA_ENUM_TYPE k = 0; k < vector_exchange_recv[0]; k++) //receive row data
				{
					INMOST_DATA_ENUM_TYPE local_size = 0;
					for(INMOST_DATA_ENUM_TYPE r = 0; r < vector_exchange_recv[j+1]; r++)
						local_size += recv_row_sizes[q+r];
					GUARD_MPI(MPI_Irecv(&recv_row_data[f],local_size,Sparse::GetRowEntryType(),vector_exchange_recv[j],4*size+vector_exchange_recv[j],comm,&requests[k]));
					q += vector_exchange_recv[j+1];
					j += vector_exchange_recv[j+1]+2;
					f += local_size;
				}

	
				j = 1;
				q = 0;
				f = 0;
				for(INMOST_DATA_ENUM_TYPE k = 0; k < vector_exchange_send[0]; k++) //receive row data
				{
					INMOST_DATA_ENUM_TYPE local_size = 0;
					for(INMOST_DATA_ENUM_TYPE r = 0; r < vector_exchange_send[j+1]; r++)
						local_size += send_row_sizes[q+r];
					GUARD_MPI(MPI_Isend(&send_row_data[f],local_size,Sparse::GetRowEntryType(),vector_exchange_send[j],4*size+rank,comm,&requests[k+vector_exchange_recv[0]]));
					q += vector_exchange_send[j+1];
					j += vector_exchange_send[j+1]+2;
					f += local_size;
				}

	
				local_start = local_end;
				m.SetInterval(local_matrix_begin,ext_pos);
				local_end = ext_pos;
				if( vector_exchange_recv[0]+vector_exchange_send[0] > 0 )
					GUARD_MPI(MPI_Waitall(vector_exchange_recv[0]+vector_exchange_send[0],&requests[0],MPI_STATUSES_IGNORE));
				j = 1;
				q = 0;
				f = 0;
				for(INMOST_DATA_ENUM_TYPE k = 0; k < vector_exchange_recv[0]; k++) //extend matrix
				{
					for(INMOST_DATA_ENUM_TYPE r = 0; r < vector_exchange_recv[j+1]; r++)
					{
						m[global_to_local[vector_exchange_recv[j+2+r]]] = Sparse::Row(&recv_row_data[f],&recv_row_data[f]+recv_row_sizes[q+r]);
						f += recv_row_sizes[q+r];
					}
					q += vector_exchange_recv[j+1];
					j += vector_exchange_recv[j+1]+2;
				}
			}
			//std::cout << it << "/" << overlap << " done" << std::endl;
		}
		two[0] = local_matrix_begin;
		two[1] = local_end;
		GUARD_MPI(MPI_Allgather(two,2,INMOST_MPI_DATA_ENUM_TYPE,&global_overlap[0],2,INMOST_MPI_DATA_ENUM_TYPE,comm));
		//std::cout << __FUNCTION__ << " done" << std::endl;
#endif
	}

	void Solver::OrderInfo::RestoreMatrix(Sparse::Matrix & m)
	{
		//restore matrix size
		m.SetInterval(initial_matrix_begin,initial_matrix_end);
		//restore indexes
		for(Sparse::Matrix::iterator it = m.Begin(); it != m.End(); ++it)
			for(Sparse::Row::iterator jt = it->Begin(); jt != it->End(); ++jt)
				if( jt->first >= initial_matrix_end )
					jt->first = extended_indexes[jt->first-initial_matrix_end];
		m.isParallel() = false;
		have_matrix = false;
#if defined(USE_MPI)
		if( comm != INMOST_MPI_COMM_WORLD )
		{
			MPI_Comm_free(&comm);
			comm = INMOST_MPI_COMM_WORLD;
		}
#endif
		//std::cout << __FUNCTION__ << std::endl;
	}

	Solver::OrderInfo::~OrderInfo()
	{
#if defined(USE_MPI)
		if( comm != INMOST_MPI_COMM_WORLD )
			MPI_Comm_free(&comm);
#endif
	}
	
	void Solver::OrderInfo::Clear()
	{
		global_to_proc.clear();
		global_overlap.clear();
		vector_exchange_recv.clear();
		vector_exchange_send.clear();
		send_storage.clear();
		recv_storage.clear();
		send_requests.clear();
		recv_requests.clear();
		extended_indexes.clear();
		local_vector_begin = local_vector_end = 0;
		initial_matrix_begin = initial_matrix_end = 0;
		local_matrix_begin = local_matrix_end = 0;
		have_matrix = false;
	}
	
	void Solver::OrderInfo::PrepareVector(Sparse::Vector & v) const
	{
		if( !have_matrix ) throw PrepareMatrixFirst;
		v.SetInterval(local_vector_begin,local_vector_end);
		v.isParallel() = true;
	}
	
	void Solver::OrderInfo::RestoreVector(Sparse::Vector & v) const
	{
		assert(have_matrix);
		if( v.isParallel() )
		{
			v.SetInterval(initial_matrix_begin,initial_matrix_end);
			v.isParallel() = false;
		}
	}
	
	Solver::OrderInfo::OrderInfo() : 
	global_to_proc(), 
	global_overlap(),
	vector_exchange_recv(), 
	vector_exchange_send(), 
	send_storage(), 
	recv_storage(), 
	send_requests(), 
	recv_requests(),
	extended_indexes()
	{
		comm = INMOST_MPI_COMM_WORLD;
		rank = 0;
		size = 1;
		initial_matrix_begin = 0;
		initial_matrix_end = 0;
		local_matrix_begin = 0;
		local_matrix_end = 0;
		local_vector_begin = 0;
		local_vector_end = 0;
		have_matrix = false;
	}
	
	Solver::OrderInfo::OrderInfo(const OrderInfo & other) 
		:global_to_proc(other.global_to_proc), global_overlap(other.global_overlap),
		vector_exchange_recv(other.vector_exchange_recv), vector_exchange_send(other.vector_exchange_send), 
		extended_indexes(other.extended_indexes)
	{
#if defined(USE_MPI)
		if( other.comm == INMOST_MPI_COMM_WORLD )
			comm = INMOST_MPI_COMM_WORLD;
		else MPI_Comm_dup(other.comm,&comm);
#else
		comm = other.comm;
#endif
		rank = other.rank;
		size = other.size;
		initial_matrix_begin = other.initial_matrix_begin;
		initial_matrix_end = other.initial_matrix_end;
		local_vector_begin = other.local_vector_begin;
		local_vector_end = other.local_vector_end;
		local_matrix_begin = other.local_matrix_begin;
		local_matrix_end = other.local_matrix_end;
		have_matrix = other.have_matrix; 
		send_storage.resize(other.send_storage.size());
		recv_storage.resize(other.recv_storage.size());
		send_requests.resize(other.send_requests.size());
		recv_requests.resize(other.recv_requests.size());
	}
	
	Solver::OrderInfo & Solver::OrderInfo::operator =(OrderInfo const & other) 
	{
#if defined(USE_MPI)
		if( other.comm == INMOST_MPI_COMM_WORLD )
			comm = INMOST_MPI_COMM_WORLD;
		else MPI_Comm_dup(other.comm,&comm);
#else
		comm = other.comm;
#endif
		global_to_proc = other.global_to_proc; 
		global_overlap = other.global_overlap; 
		vector_exchange_recv = other.vector_exchange_recv;
		vector_exchange_send = other.vector_exchange_send; 
		extended_indexes = other.extended_indexes;
		rank = other.rank;
		size = other.size;
		initial_matrix_begin = other.initial_matrix_begin;
		initial_matrix_end = other.initial_matrix_end;
		local_vector_begin = other.local_vector_begin;
		local_vector_end = other.local_vector_end;
		local_matrix_begin = other.local_matrix_begin;
		local_matrix_end = other.local_matrix_end;
		have_matrix = other.have_matrix; 
		send_storage.resize(other.send_storage.size());
		recv_storage.resize(other.recv_storage.size());
		send_requests.resize(other.send_requests.size());
		recv_requests.resize(other.recv_requests.size());
		return *this;
	}


	INMOST_DATA_ENUM_TYPE Solver::OrderInfo::GetProcessor(INMOST_DATA_ENUM_TYPE gind) const
	{
		assert(have_matrix);
		storage_type::const_iterator find = std::lower_bound(global_to_proc.begin(),global_to_proc.end(),gind);
		assert(find != global_to_proc.end());
		if( (*find) == gind && find+1 != global_to_proc.end()) return static_cast<INMOST_DATA_ENUM_TYPE>(find - global_to_proc.begin());
		else return static_cast<INMOST_DATA_ENUM_TYPE>(find - global_to_proc.begin())-1;
	}
	void Solver::OrderInfo::GetOverlapRegion(INMOST_DATA_ENUM_TYPE proc, INMOST_DATA_ENUM_TYPE & mbeg, INMOST_DATA_ENUM_TYPE & mend) const
	{
		assert(have_matrix); 
		mbeg = global_overlap[proc*2+0];
		mend = global_overlap[proc*2+1];
	}
	void Solver::OrderInfo::GetLocalRegion(INMOST_DATA_ENUM_TYPE proc, INMOST_DATA_ENUM_TYPE & mbeg, INMOST_DATA_ENUM_TYPE & mend) const
	{
		assert(have_matrix);
		mbeg = global_to_proc[proc+0];
		mend = global_to_proc[proc+1];
	}
	
	
	void Solver::OrderInfo::Update(Sparse::Vector & x)
	{
		//std::cout << __FUNCTION__ << " start" << std::endl;
#if defined(USE_MPI)
		if( GetSize() == 1 ) return;
#if defined(USE_OMP)
#pragma omp single
#endif
		{
			//use MPI_Put/MPI_Get to update vector
			assert(x.isParallel()); //the vector was prepared
			INMOST_DATA_ENUM_TYPE i, j = 1, k, l = 0;
			int ierr;
			for(i = 0; i < vector_exchange_recv[0]; i++)
			{
				//std::cout << GetRank() << " MPI_Irecv size " << vector_exchange_recv[j+1] << " dest " << vector_exchange_recv[j] << " tag " << vector_exchange_recv[j]*size+rank << std::endl;
				GUARD_MPI(MPI_Irecv(&recv_storage[l],vector_exchange_recv[j+1],INMOST_MPI_DATA_REAL_TYPE,vector_exchange_recv[j],vector_exchange_recv[j]*size+rank,comm,&recv_requests[i]));
				l += vector_exchange_recv[j+1];
				j += vector_exchange_recv[j+1] + 2;
			}
			j = 1, l = 0;
			for(i = 0; i < vector_exchange_send[0]; i++)
			{
				//std::cout << GetRank() << " MPI_Isend size " << vector_exchange_send[j+1] << " dest " << vector_exchange_send[j] << " tag " << rank*size+vector_exchange_send[j] << std::endl;
				for(k = 0; k < vector_exchange_send[j+1]; k++)
					send_storage[l+k] = x[vector_exchange_send[k+j+2]];
				GUARD_MPI(MPI_Isend(&send_storage[l],vector_exchange_send[j+1],INMOST_MPI_DATA_REAL_TYPE,vector_exchange_send[j],rank*size+vector_exchange_send[j],comm,&send_requests[i]));
				l += vector_exchange_send[j+1];
				j += vector_exchange_send[j+1] + 2;
			}
			if( vector_exchange_recv[0] > 0 )
			{
				GUARD_MPI(MPI_Waitall(static_cast<int>(recv_requests.size()),&recv_requests[0],MPI_STATUSES_IGNORE));
				j = 1, l = 0;
				for(i = 0; i < vector_exchange_recv[0]; i++)
				{
					for(k = 0; k < vector_exchange_recv[j+1]; k++)
						x[vector_exchange_recv[k+j+2]] = recv_storage[l+k];
					l += vector_exchange_recv[j+1];
					j += vector_exchange_recv[j+1] + 2;
				}
			}
			if( vector_exchange_send[0] > 0 ) 
			{
				GUARD_MPI(MPI_Waitall(static_cast<int>(send_requests.size()),&send_requests[0],MPI_STATUSES_IGNORE));
			}
		}
#else
		(void) x;
#endif
		//std::cout << __FUNCTION__ << " end" << std::endl;
	}
	void Solver::OrderInfo::Accumulate(Sparse::Vector & x)
	{
		//std::cout << __FUNCTION__ << " start" << std::endl;
#if defined(USE_MPI)
		if( GetSize() == 1 ) return;
#if defined(USE_OMP)
#pragma omp single
#endif
		{
			//use MPI_Put/MPI_Get to update vector
			assert(x.isParallel()); //the vector was prepared
			INMOST_DATA_ENUM_TYPE i, j = 1, k, l = 0;
			int ierr;
			for(i = 0; i < vector_exchange_send[0]; i++)
			{
				//std::cout << GetRank() << " MPI_Irecv size " << vector_exchange_send[j+1] << " dest " << vector_exchange_send[j] << " tag " << vector_exchange_send[j]*size+rank << std::endl;
				GUARD_MPI(MPI_Irecv(&send_storage[l],vector_exchange_send[j+1],INMOST_MPI_DATA_REAL_TYPE,vector_exchange_send[j],vector_exchange_send[j]*size+rank,comm,&send_requests[i]));
				l += vector_exchange_send[j+1];
				j += vector_exchange_send[j+1] + 2;
			}
			j = 1, l = 0;
			for(i = 0; i < vector_exchange_recv[0]; i++)
			{
				for(k = 0; k < vector_exchange_recv[j+1]; k++)
					recv_storage[l+k] = x[vector_exchange_recv[k+j+2]];
				//std::cout << GetRank() << " MPI_Isend size " << vector_exchange_recv[j+1] << " dest " << vector_exchange_recv[j] << " tag " << rank*size+vector_exchange_recv[j] << std::endl;
				GUARD_MPI(MPI_Isend(&recv_storage[l],vector_exchange_recv[j+1],INMOST_MPI_DATA_REAL_TYPE,vector_exchange_recv[j],rank*size+vector_exchange_recv[j],comm,&recv_requests[i]));
				l += vector_exchange_recv[j+1];
				j += vector_exchange_recv[j+1] + 2;
			}
			if( vector_exchange_send[0] > 0 )
			{
				//std::cout << GetRank() << " Waitall send " << send_requests.size() << std::endl;
				GUARD_MPI(MPI_Waitall(static_cast<int>(send_requests.size()),&send_requests[0],MPI_STATUSES_IGNORE));
				j = 1, l = 0;
				for(i = 0; i < vector_exchange_send[0]; i++)
				{
					for(k = 0; k < vector_exchange_send[j+1]; k++)
						x[vector_exchange_send[k+j+2]] += send_storage[l+k];
					l += vector_exchange_send[j+1];
					j += vector_exchange_send[j+1] + 2;
				}
			}
			if( vector_exchange_recv[0] > 0 ) 
			{
				//std::cout << GetRank() << " Waitall recv " << recv_requests.size() << std::endl;
				GUARD_MPI(MPI_Waitall(static_cast<int>(recv_requests.size()),&recv_requests[0],MPI_STATUSES_IGNORE));
			}
		}
#else
		(void) x;
#endif
		//std::cout << __FUNCTION__ << " end" << std::endl;
	}
	/*
	void Solver::OrderInfo::ScalarProd(Vector const & left, Vector const & right, INMOST_DATA_ENUM_TYPE index_begin, INMOST_DATA_ENUM_TYPE index_end, INMOST_DATA_REAL_TYPE & sum) const
	{
		INMOST_DATA_INTEGER_TYPE ibeg = index_begin, iend = index_end;
#if defined(USE_OMP)
#pragma omp for reduction(+:sum)
#endif
		for(INMOST_DATA_INTEGER_TYPE i = ibeg; i < iend; ++i)
		{
			sum += left[i]*right[i];
		}
		Integrate(&sum,1);
	}
	*/
	
	
	
	



	

	//======================================================================
	//SOLVER SOURCE CODE
	//======================================================================
	void Solver::Initialize(int * argc, char *** argv, const char * database)
	{
		(void)database;
		(void)argc;
		(void)argv;
		if( database != NULL )
		{
			FILE * f = fopen(database,"r");
			if( f )
			{
				//std::fstream file(database,std::ios::in);
				char str[4096];
				//while( !file.eof() && file.getline(str,4096) )
				while( !feof(f) && fgets(str,4096,f) )
				{
					int k = 0, l;
					for(k = 0; k < (int)strlen(str); ++k)
					{
						if( str[k] == ':' ) break;
					}
					if( k == strlen(str) ) continue; //invalid line
					for(l = 0; l < k; ++l) str[l] = tolower(str[l]);
					l = (int)strlen(str)-1; // Right-trim string
					while(l > 0 && isspace(str[l]) ) --l;
					str[l+1] = 0;
					l = k+1;
					while(l < (int)strlen(str) && isspace(str[l]) ) ++l;
					if( l == strlen(str) ) continue; //skip empty entry
					if( !strncmp(str,"petsc",k) )
						petsc_database_file = std::string(str+l);
					else if( !strncmp(str,"trilinos_ifpack",k) )
						trilinos_ifpack_database_file = std::string(str+l);
					else if( !strncmp(str,"trilinos_aztec",k) )
						trilinos_aztec_database_file = std::string(str+l);
					else if( !strncmp(str,"trilinos_ml",k) )
						trilinos_ml_database_file = std::string(str+l);
					else if( !strncmp(str,"trilinos_belos",k) )
						trilinos_belos_database_file = std::string(str+l);
					else if( !strncmp(str,"ani",k) )
						ani_database_file = std::string(str+l);
					else if( !strncmp(str,"fcbiilu2",k) )
						fcbiilu2_database_file = std::string(str+l);
					else if( !strncmp(str,"k3biilu2",k) )
						k3biilu2_database_file = std::string(str+l);
					else if( !strncmp(str,"superlu",k) )
						superlu_database_file = std::string(str+l);
				}
				//file.close();
				fclose(f);
			}
		}
		//std::cout << "PETSc \"" << petsc_database_file << "\"" << std::endl;
		//std::cout << "Trilinos_Ifpack \"" << trilinos_ifpack_database_file << "\"" << std::endl;
#if defined(USE_SOLVER_PETSC)
		SolverInitializePetsc(argc,argv,petsc_database_file.c_str());
#endif
#if defined(USE_SOLVER_SUPERLU)
		SolverInitializeSuperLU(argc,argv,superlu_database_file.c_str());
#endif
#if defined(USE_SOLVER_ANI)
		SolverInitializeAni(argc,argv,ani_database_file.c_str());
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
		SolverInitializeFcbiilu2(argc,argv,fcbiilu2_database_file.c_str());
#endif
#if defined(HAVE_SOLVER_K3BIILU2)
		SolverInitializeK3biilu2(argc,argv,k3biilu2_database_file.c_str());
#endif
#if defined(USE_MPI)
		{
			int flag = 0;
			int ierr = 0;
			MPI_Initialized(&flag);
			if( !flag )
			{
				ierr = MPI_Init(argc,argv);
				if( ierr != MPI_SUCCESS )
				{
					std::cout << __FILE__ << ":" << __LINE__ << "problem in MPI_Init" << std::endl;
				}
			}
		}
#endif
#if defined(USE_SOLVER_SUPERLU)
		
#endif
		Sparse::CreateRowEntryType();
		is_initialized = true;
	}
	
	void Solver::Finalize()
	{
    Sparse::DestroyRowEntryType();
#if defined(USE_SOLVER_PETSC)
		SolverFinalizePetsc();
#endif
#if defined(USE_SOLVER_ANI)
		SolverFinalizeAni();
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
		SolverFinalizeFcbiilu2();
#endif
#if defined(HAVE_SOLVER_K3BIILU2)
		SolverFinalizeK3biilu2();
#endif
#if defined(USE_MPI)
		{
			int flag = 0;
			MPI_Finalized(&flag);
			if( !flag ) 
				MPI_Finalize();
		}
#endif
    is_initialized = false;
	  is_finalized = true;
	}
	void Solver::SetParameterEnum(std::string name, INMOST_DATA_ENUM_TYPE val)
	{
		if( name == "maximum_iterations" )
			maximum_iterations = val;
		else if( name == "rescale_iterations" )
			preconditioner_rescale_iterations = val;
		else if( name == "condition_estimation" )
			preconditioner_condition_estimation = val;
		else if( name == "adapt_ddpq_tolerance" )
			preconditioner_adapt_ddpq_tolerance = val;
		else if( name == "schwartz_overlap" )
			additive_schwartz_overlap = val;
		else if( name == "gmres_substeps" )
			solver_gmres_substeps = val;
		else if( name == "reorder_nonzeros" )
			preconditioner_reorder_nonzero = val;
		/* Block below leads to confusion - parameters would be reset by preset values during setup
		else if( _pack == INNER_ILU2 || _pack == INNER_MLILUC )
		{
			try
			{
				//This leads to confusion - 
				IterativeMethod * method = (IterativeMethod *)solver_data;
				if (name[0] == ':')
					 method->EnumParameter(name.substr(1, name.size() - 1)) = val;
				else if (name == "overlap") additive_schwartz_overlap = val;
				else method->EnumParameter(name) = val;
			} 
			catch(...)
			{
				std::cout << "Parameter " << name << " of intergral type is unknown" << std::endl;
			}
		}
		*/
		else std::cout << "Parameter " << name << " of integral type is unknown" << std::endl;
	}
	void Solver::SetParameterReal(std::string name, INMOST_DATA_REAL_TYPE val)
	{
		if( name == "absolute_tolerance" )
			absolute_tolerance = val;
		else if( name == "relative_tolerance" )
			relative_tolerance = val;
		else if( name == "divergence_tolerance" )
			divergence_tolerance = val;
		else if( name == "drop_tolerance" )
			preconditioner_drop_tolerance = val;
		else if( name == "reuse_tolerance" )
			preconditioner_reuse_tolerance = val;
		else if( name == "ddpq_tolerance" )
			preconditioner_ddpq_tolerance = val;
		else if( name == "fill_level" )
			preconditioner_fill_level = val;
		/* Block below leads to confusion - parameters would be reset by preset values during setup
		else if( _pack == INNER_ILU2 || _pack == INNER_MLILUC )
		{
			try
			{
				IterativeMethod * method = (IterativeMethod *)solver_data;
				if (name[0] == ':')	method->RealParameter(name.substr(1, name.size() - 1)) = val;
				else method->RealParameter(name) = val;
			} 
			catch(...)
			{
				std::cout << "Parameter " << name << " of real type is unknown" << std::endl;
			}
		}
		*/
		else std::cout << "Parameter " << name << " of real type is unknown" << std::endl;
	}
	Solver::Solver(Type pack, std::string _name, INMOST_MPI_Comm _comm)
	{
		additive_schwartz_overlap = DEFAULT_ADDITIVE_SCHWARTZ_OVERLAP;
		maximum_iterations = DEFAULT_MAXIMUM_ITERATIONS;
		absolute_tolerance = DEFAULT_ABSOLUTE_TOLERANCE;
		relative_tolerance = DEFAULT_RELATIVE_TOLERANCE;
		divergence_tolerance = DEFAULT_DIVERGENCE_TOLERANCE;
		preconditioner_ddpq_tolerance = DEFAULT_PRECONDITIONER_DDPQ_TOLERANCE;
		preconditioner_drop_tolerance = DEFAULT_PRECONDITIONER_DROP_TOLERANCE;
		preconditioner_reuse_tolerance = DEFAULT_PRECONDITIONER_REUSE_TOLERANCE;
		preconditioner_reorder_nonzero = DEFAULT_PRECONDITIONER_REORDER_NONZEROS;
		preconditioner_rescale_iterations = DEFAULT_PRECONDITIONER_RESCALE_ITERS;
		preconditioner_fill_level = DEFAULT_PRECONDITIONER_FILL_LEVEL;
		preconditioner_adapt_ddpq_tolerance = DEFAULT_PRECONDITIONER_ADAPT_DDPQ_TOLERANCE;
		preconditioner_condition_estimation = DEFAULT_PRECONDITIONER_CONDITION_ESTIMATION;
		solver_gmres_substeps = DEFAULT_SOLVER_GMRES_SUBSTEPS;
		comm = _comm;
		_pack = pack;
		name = _name;
		solver_data = NULL;
		local_size = global_size = 0;
		last_it = 0;
		last_resid = 0;
		matrix_data = rhs_data = solution_data = precond_data =  NULL;
#if defined(USE_SOLVER_PETSC)
		if( _pack == PETSc )
		{
			SolverInitDataPetsc(&solver_data,_comm,name.c_str());
		}
#endif
#if defined(USE_SOLVER_ANI)
		if( _pack == ANI )
		{
			SolverInitDataAni(&solver_data,_comm,name.c_str());
		}
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
		if( _pack == FCBIILU2 )
		{
			SolverInitDataFcbiilu2(&solver_data,_comm,name.c_str());
		}
#endif
#if defined(HAVE_SOLVER_K3BIILU2)
		if( _pack == K3BIILU2 )
		{
			SolverInitDataK3biilu2(&solver_data,_comm,name.c_str());
		}
#endif
#if defined(USE_SOLVER_SUPERLU)
		if( _pack == SUPERLU )
		{
			SolverInitDataSuperLU(&solver_data);
		}
#endif
#if defined(USE_SOLVER_TRILINOS)
		if( _pack == Trilinos_Aztec || 
			_pack == Trilinos_Belos ||
			_pack == Trilinos_ML ||
			_pack == Trilinos_Ifpack)
		{
			std::pair<std::string,Epetra_LinearProblem> * problem = 
				new std::pair<std::string,Epetra_LinearProblem>;
			problem->first = _name;
			solver_data = static_cast<void *>(problem);
		}
#endif
		if( _pack == INNER_ILU2 || _pack == INNER_DDPQILUC || _pack == INNER_MPTILUC || _pack == INNER_MPTILU2)
		{
			Method * prec;
			if (_pack == INNER_ILU2)
			{
#if defined(__SOLVER_ILU2__)
				prec = new ILU2_preconditioner(info);
#else
				std::cout << "Sorry, ILU2 preconditioner is not included in this release" << std::endl;
				prec = NULL;
#endif
			}
			else if( _pack == INNER_DDPQILUC )
			{
#if defined(__SOLVER_DDPQILUC2__)
				prec = new ILUC_preconditioner(info);
#else
				std::cout << "Sorry, multilevel condition estimation crout-ilu preconditioner with reordering for diagonal dominance is not included in this release" << std::endl;
				prec = NULL;
#endif
			}
			else if( _pack == INNER_MPTILUC )
			{
#if defined(__SOLVER_MTILUC2__)
				prec = new MTILUC_preconditioner(info);
#else
				std::cout << "Sorry, maximum product transverse condition estimation crout-ilu preconditioner is not included in this release" << std::endl;
				prec = NULL;
#endif
			}
			else if( _pack == INNER_MPTILU2 )
			{
#if defined(__SOLVER_MTILU2__)
				prec = new MTILU2_preconditioner(info);
#else
				std::cout << "Sorry, maximum product transverse ilu2 preconditioner is not included in this release" << std::endl;
				prec = NULL;
#endif
			}
            else prec = NULL;
			solver_data = new KSOLVER(prec, info);
		}
	}
	Solver::Solver(const Solver & other)
	{
		comm = other.comm;
#if defined(USE_SOLVER_PETSC)
		if( _pack == PETSc )
		{
			SolverCopyDataPetsc(&solver_data,other.solver_data,comm);
			if( other.matrix_data != NULL ) 
			{
				MatrixCopyDataPetsc(&matrix_data,other.matrix_data);
				SolverSetMatrixPetsc(solver_data,matrix_data,false,false);
			}
			if( other.rhs_data != NULL ) 
				VectorCopyDataPetsc(&rhs_data,other.rhs_data);
			if( other.solution_data != NULL ) 
				VectorCopyDataPetsc(&solution_data,other.solution_data);
		}
#endif
#if defined(USE_SOLVER_ANI)
		if( _pack == ANI )
		{
			SolverCopyDataAni(&solver_data,other.solver_data,comm);
			if( other.matrix_data != NULL ) 
			{
				MatrixCopyDataAni(&matrix_data,other.matrix_data);
				SolverSetMatrixAni(solver_data,matrix_data,false,false);
			}
			if( other.rhs_data != NULL ) 
				VectorCopyDataAni(&rhs_data,other.rhs_data);
			if( other.solution_data != NULL ) 
				VectorCopyDataAni(&solution_data,other.solution_data);
		}
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
		if( _pack == FCBIILU2 )
		{
			SolverCopyDataFcbiilu2(&solver_data,other.solver_data,comm);
			if( other.matrix_data != NULL ) 
			{
				MatrixCopyDataFcbiilu2(&matrix_data,other.matrix_data);
				SolverSetMatrixFcbiilu2(solver_data,matrix_data,false,false);
			}
			if( other.rhs_data != NULL ) 
				VectorCopyDataFcbiilu2(&rhs_data,other.rhs_data);
			if( other.solution_data != NULL ) 
				VectorCopyDataFcbiilu2(&solution_data,other.solution_data);
		}
#endif
#if defined(HAVE_SOLVER_K3BIILU2)
		if( _pack == K3BIILU2 )
		{
			SolverCopyDataK3biilu2(&solver_data,other.solver_data,comm);
			if( other.matrix_data != NULL ) 
			{
				MatrixCopyDataK3biilu2(&matrix_data,other.matrix_data);
				SolverSetMatrixK3biilu2(solver_data,matrix_data,false,false);
			}
			if( other.rhs_data != NULL ) 
				VectorCopyDataK3biilu2(&rhs_data,other.rhs_data);
			if( other.solution_data != NULL ) 
				VectorCopyDataK3biilu2(&solution_data,other.solution_data);
		}
#endif
#if defined(USE_SOLVER_TRILINOS)
		if( _pack == Trilinos_Aztec || 
			_pack == Trilinos_Belos ||
			_pack == Trilinos_ML ||
			_pack == Trilinos_Ifpack)
		{
			throw - 1; //You should not really want to copy solver's information
		}
#endif
		if (_pack == INNER_ILU2 || _pack == INNER_DDPQILUC || _pack == INNER_MPTILUC || _pack == INNER_MPTILU2)
		{
			throw - 1; //You should not really want to copy solver's information
		}
	}
	void Solver::Clear()
	{
		local_size = global_size = 0;
		info.Clear();
#if defined(USE_SOLVER_PETSC)
		if( _pack == PETSc )
		{
			if( matrix_data != NULL ) 
				MatrixDestroyDataPetsc(&matrix_data);
			if( rhs_data != NULL )
				VectorDestroyDataPetsc(&rhs_data);
			if( solution_data != NULL )
				VectorDestroyDataPetsc(&solution_data);
			SolverDestroyDataPetsc(&solver_data);
		}
#endif
#if defined(USE_SOLVER_ANI)
		if( _pack == ANI )
		{
			if( matrix_data != NULL ) 
				MatrixDestroyDataAni(&matrix_data);
			if( rhs_data != NULL )
				VectorDestroyDataAni(&rhs_data);
			if( solution_data != NULL )
				VectorDestroyDataAni(&solution_data);
			SolverDestroyDataAni(&solver_data);
		}
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
		if( _pack == FCBIILU2 )
		{
			if( matrix_data != NULL ) 
				MatrixDestroyDataFcbiilu2(&matrix_data);
			if( rhs_data != NULL )
				VectorDestroyDataFcbiilu2(&rhs_data);
			if( solution_data != NULL )
				VectorDestroyDataFcbiilu2(&solution_data);
			SolverDestroyDataFcbiilu2(&solver_data);
		}
#endif
#if defined(HAVE_SOLVER_K3BIILU2)
		if( _pack == K3BIILU2 )
		{
			if( matrix_data != NULL ) 
				MatrixDestroyDataK3biilu2(&matrix_data);
			if( rhs_data != NULL )
				VectorDestroyDataK3biilu2(&rhs_data);
			if( solution_data != NULL )
				VectorDestroyDataK3biilu2(&solution_data);
			SolverDestroyDataK3biilu2(&solver_data);
		}
#endif
#if defined(USE_SOLVER_TRILINOS)
		if( _pack == Trilinos_Aztec || 
			_pack == Trilinos_Belos ||
			_pack == Trilinos_ML ||
			_pack == Trilinos_Ifpack)
		{
			if( matrix_data != NULL )
			{
				delete static_cast<Epetra_CrsMatrix *>(matrix_data);
				matrix_data = NULL;
			}
			if( solver_data != NULL )
			{
				delete static_cast<std::pair<std::string,Epetra_LinearProblem> *>(solver_data);
				solver_data = NULL;
			}
		}
#endif
#if defined(USE_SOLVER_SUPERLU)
		if(_pack == SUPERLU )
		{
			SolverDestroyDataSuperLU(&solver_data);
			matrix_data = NULL;
		}
#endif
		if (_pack == INNER_ILU2 || _pack == INNER_DDPQILUC || _pack == INNER_MPTILUC || _pack == INNER_MPTILU2)
		{
			if( matrix_data != NULL )
			{
				delete (Sparse::Matrix * )matrix_data;
				matrix_data = NULL;
			}
			if( solver_data != NULL ) 
			{
				delete (Method *)solver_data;
				solver_data = NULL;
			}
		}
	}
	Solver::~Solver()
	{
		Clear();
	}
	Solver & Solver::operator =(Solver const & other)
	{
		if( this != &other )
		{
			comm = other.comm;
#if defined(USE_SOLVER_PETSC)
			if( _pack == PETSc )
			{
				SolverAssignDataPetsc(solver_data,other.solver_data);
				if( other.matrix_data != NULL ) 
				{
					if( matrix_data != NULL )
						MatrixAssignDataPetsc(matrix_data,other.matrix_data);
					else
						MatrixCopyDataPetsc(&matrix_data,other.matrix_data); 
					SolverSetMatrixPetsc(solver_data,matrix_data,false,false);
				}
				if( other.rhs_data != NULL ) 
				{
					if( rhs_data != NULL )
						VectorAssignDataPetsc(rhs_data,other.rhs_data);
					else
						VectorCopyDataPetsc(&rhs_data,other.rhs_data);
				}
				if( other.solution_data != NULL ) 
				{
					if( solution_data != NULL )
						VectorAssignDataPetsc(solution_data,other.solution_data);
					else
						VectorCopyDataPetsc(&solution_data,other.solution_data);
				}
			}
#endif
#if defined(USE_SOLVER_ANI)
			if( _pack == ANI )
			{
				SolverAssignDataAni(solver_data,other.solver_data);
				if( other.matrix_data != NULL ) 
				{
					if( matrix_data != NULL )
						MatrixAssignDataAni(matrix_data,other.matrix_data);
					else
						MatrixCopyDataAni(&matrix_data,other.matrix_data); 
					SolverSetMatrixAni(solver_data,matrix_data,false,false);
				}
				if( other.rhs_data != NULL ) 
				{
					if( rhs_data != NULL )
						VectorAssignDataAni(rhs_data,other.rhs_data);
					else
						VectorCopyDataAni(&rhs_data,other.rhs_data);
				}
				if( other.solution_data != NULL ) 
				{
					if( solution_data != NULL )
						VectorAssignDataAni(solution_data,other.solution_data);
					else
						VectorCopyDataAni(&solution_data,other.solution_data);
				}
			}
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
			if( _pack == FCBIILU2 )
			{
				SolverAssignDataFcbiilu2(solver_data,other.solver_data);
				if( other.matrix_data != NULL ) 
				{
					if( matrix_data != NULL )
						MatrixAssignDataFcbiilu2(matrix_data,other.matrix_data);
					else
						MatrixCopyDataFcbiilu2(&matrix_data,other.matrix_data); 
					SolverSetMatrixFcbiilu2(solver_data,matrix_data,false,false);
				}
				if( other.rhs_data != NULL ) 
				{
					if( rhs_data != NULL )
						VectorAssignDataFcbiilu2(rhs_data,other.rhs_data);
					else
						VectorCopyDataFcbiilu2(&rhs_data,other.rhs_data);
				}
				if( other.solution_data != NULL ) 
				{
					if( solution_data != NULL )
						VectorAssignDataFcbiilu2(solution_data,other.solution_data);
					else
						VectorCopyDataFcbiilu2(&solution_data,other.solution_data);
				}
			}
#endif
#if defined(HAVE_SOLVER_K3BIILU2)
			if( _pack == K3BIILU2 )
			{
				SolverAssignDataK3biilu2(solver_data,other.solver_data);
				if( other.matrix_data != NULL ) 
				{
					if( matrix_data != NULL )
						MatrixAssignDataK3biilu2(matrix_data,other.matrix_data);
					else
						MatrixCopyDataK3biilu2(&matrix_data,other.matrix_data); 
					SolverSetMatrixK3biilu2(solver_data,matrix_data,false,false);
				}
				if( other.rhs_data != NULL ) 
				{
					if( rhs_data != NULL )
						VectorAssignDataK3biilu2(rhs_data,other.rhs_data);
					else
						VectorCopyDataK3biilu2(&rhs_data,other.rhs_data);
				}
				if( other.solution_data != NULL ) 
				{
					if( solution_data != NULL )
						VectorAssignDataK3biilu2(solution_data,other.solution_data);
					else
						VectorCopyDataK3biilu2(&solution_data,other.solution_data);
				}
			}
#endif
#if defined(USE_SOLVER_TRILINOS)
			if( _pack == Trilinos_Aztec || 
				_pack == Trilinos_Belos ||
				_pack == Trilinos_ML ||
				_pack == Trilinos_Ifpack)
			{
				throw - 1; //You should not really want to copy solver's information
			}
#endif
			if (_pack == INNER_ILU2 || _pack == INNER_DDPQILUC || _pack == INNER_MPTILUC || _pack == INNER_MPTILU2)
			{
				throw - 1; //You should not really want to copy solver's information
			}
		}
		return *this;
	}
	
	void Solver::SetMatrix(Sparse::Matrix & A, bool ModifiedPattern, bool OldPreconditioner)
	{
		(void) OldPreconditioner;
		bool ok = false;
#if defined(USE_SOLVER_PETSC)
		if( _pack == PETSc )
		{
			bool modified_pattern = ModifiedPattern;
			//~ if( A.comm != comm ) throw DifferentCommunicatorInSolver;
			if( matrix_data == NULL ) 
			{
				MatrixInitDataPetsc(&matrix_data,A.GetCommunicator(),A.GetName().c_str());
				modified_pattern = true;
				//printf("matrix %p\n",matrix_data);
			}
			INMOST_DATA_ENUM_TYPE mbeg,mend;
			A.GetInterval(mbeg,mend);
			if( modified_pattern )
			{				
				local_size = A.Size();
#if defined(USE_MPI)
				MPI_Allreduce(&local_size,&global_size,1,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,comm);
#else
				global_size = local_size;
#endif
				int max = 0;
				{
					int * diag_nonzeroes = new int[local_size];
					int * off_diag_nonzeroes = new int[local_size];
					unsigned k = 0;
					
					for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); ++it)
					{
						diag_nonzeroes[k] = off_diag_nonzeroes[k] = 0;
						for(INMOST_DATA_ENUM_TYPE i = 0; i < it->Size(); i++)
						{
							INMOST_DATA_ENUM_TYPE index = it->GetIndex(i);
							if( index < mbeg || index > mend-1 ) off_diag_nonzeroes[k]++;
							else diag_nonzeroes[k]++;
						}
						if( diag_nonzeroes[k] + off_diag_nonzeroes[k] > max ) max = diag_nonzeroes[k] + off_diag_nonzeroes[k];
						k++;
					}
					MatrixPreallocatePetsc(matrix_data,local_size,global_size,diag_nonzeroes,off_diag_nonzeroes);
					delete [] diag_nonzeroes;
					delete [] off_diag_nonzeroes;
				}
				if( max > 0 )
				{
					int * col_positions = new int[max];
					double * col_values = new double[max];
					unsigned k = 0, m;
					for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); ++it)
					{
						m = 0;
						for(INMOST_DATA_ENUM_TYPE i = 0; i < it->Size(); i++)
						{
							col_positions[m] = it->GetIndex(i);
							col_values[m] = it->GetValue(i);
							m++;
						}
						MatrixFillPetsc(matrix_data,mbeg+k,m,col_positions,col_values);
						k++;
					}
					delete [] col_positions;
					delete [] col_values;
				}
			}
			else
			{
				unsigned max = 0;
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); ++it)
					if( it->Size() > max ) max = it->Size();
				int * col_positions = new int[max];
				double * col_values = new double[max];
				unsigned k = 0, m;
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); ++it)
				{
					m = 0;
					for(INMOST_DATA_ENUM_TYPE i = 0; i < it->Size(); i++)
					{
						col_positions[m] = it->GetIndex(i);
						col_values[m] = it->GetValue(i);
						m++;
					}
					MatrixFillPetsc(matrix_data,mbeg+k,m,col_positions,col_values);
					k++;
				}
				delete [] col_positions;
				delete [] col_values;
			}
			MatrixFinalizePetsc(matrix_data);
			if( petsc_database_file == "" ) //set parameters if no file is set
			{
				SolverSetDropTolerancePetsc(solver_data,preconditioner_drop_tolerance);
				SolverSetFillLevelPetsc(solver_data,preconditioner_fill_level);
				SolverSetOverlapPetsc(solver_data,additive_schwartz_overlap);
			}
			SolverSetMatrixPetsc(solver_data,matrix_data,modified_pattern,OldPreconditioner);
			ok = true;
		}
#endif
#if defined(USE_SOLVER_ANI)
		if( _pack == ANI )
		{
			bool modified_pattern = ModifiedPattern;
			if( matrix_data == NULL )
			{ 
				MatrixInitDataAni(&matrix_data,A.GetCommunicator(),A.GetName().c_str());
				modified_pattern = true;
			}
			if( modified_pattern )
			{
				global_size = local_size = A.Size();
				
				MatrixDestroyDataAni(&matrix_data);
				MatrixInitDataAni(&matrix_data,A.GetCommunicator(),A.GetName().c_str());
				INMOST_DATA_ENUM_TYPE nnz = 0, k = 0, q = 1, shift = 1;
				int n = A.Size();
				int * ia = (int *)malloc(sizeof(int)*(n+1));
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) nnz += it->Size();
				int * ja = (int *)malloc(sizeof(int)*nnz);
				double * values = (double *) malloc(sizeof(double)*nnz);
				ia[0] = shift;
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++)
				{
					for(Sparse::Row::iterator jt = it->Begin(); jt != it->End(); jt++)
					{
						ja[k] = jt->first + 1;
						values[k] = jt->second;
						k++;
					}
					shift += it->Size();
					ia[q++] = shift;
				}
				MatrixFillAni(matrix_data,n,ia,ja,values);
			}
			else
			{
				INMOST_DATA_ENUM_TYPE nnz = 0, k = 0;
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) nnz += it->Size();
				double * values = (double *)malloc(sizeof(double)*nnz);
				k = 0;
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++)
					for(Sparse::Row::iterator jt = it->Begin(); jt != it->End(); jt++)
						values[k++] = jt->second;
				MatrixFillValuesAni(matrix_data,values);
			}
			MatrixFinalizeAni(matrix_data);
			SolverSetMatrixAni(solver_data,matrix_data,modified_pattern,OldPreconditioner);
			ok = true;
		}
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
		if( _pack == FCBIILU2 )
		{
			bool modified_pattern = ModifiedPattern;
			//~ if( A.comm != comm ) throw DifferentCommunicatorInSolver;
			if( matrix_data == NULL )
			{ 
				MatrixInitDataFcbiilu2(&matrix_data,A.GetCommunicator(),A.GetName().c_str());
				modified_pattern = true;
			}
			if( modified_pattern )
			{
				global_size = local_size = A.Size();
				
				MatrixDestroyDataFcbiilu2(&matrix_data);
				MatrixInitDataFcbiilu2(&matrix_data,A.GetCommunicator(),A.GetName().c_str());
				INMOST_DATA_ENUM_TYPE nnz = 0, k = 0, q = 1, shift = 1;
				int nproc, myid;
#if defined(USE_MPI)
				MPI_Comm_size(A.GetCommunicator(),&nproc);
				MPI_Comm_rank(A.GetCommunicator(),&myid);
#else
				nproc = 1;
				myid = 0;
#endif
				int * ibl = (int *)malloc(sizeof(int)*(nproc+1));
				ibl[0] = 0;
				int n = A.Size();
#if defined(USE_MPI)
				INMOST_DATA_ENUM_TYPE mbeg,mend;
				A.GetInterval(mbeg,mend);
				n = mend - mbeg;
				int block_end = mend;
				MPI_Allgather(&block_end,1,MPI_INT,&ibl[1],1,MPI_INT,A.GetCommunicator());
#else
				ibl[1] = n;
#endif
				int * ia = (int *)malloc(sizeof(int)*(n+1));
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) nnz += it->Size();
				int * ja = (int *)malloc(sizeof(int)*nnz);
				double * values = (double *) malloc(sizeof(double)*nnz);
				ia[0] = shift;
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++)
				{
					for(Row::iterator jt = it->Begin(); jt != it->End(); jt++)
					{
						ja[k] = jt->first + 1;
						values[k] = jt->second;
		                                //std::cout<<"# q="<<q<<" k="<<k<<" ja="<<ja[k]<<" a="<<values[k]<<std::endl;//db!
						k++;
					}
					shift += it->Size();
					ia[q++] = shift;
				}
				MatrixFillFcbiilu2(matrix_data,n,nproc,ibl,ia,ja,values);
			}
			else
			{
				INMOST_DATA_ENUM_TYPE nnz = 0, k = 0;
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) nnz += it->Size();
				double * values = (double *)malloc(sizeof(double)*nnz);
				k = 0;
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++)
					for(Sparse::Row::iterator jt = it->Begin(); jt != it->End(); jt++)
						values[k++] = jt->second;
				MatrixFillValuesFcbiilu2(matrix_data,values);
			}
			MatrixFinalizeFcbiilu2(matrix_data);
			SolverSetMatrixFcbiilu2(solver_data,matrix_data,modified_pattern,OldPreconditioner);
			ok = true;
		}
#endif
#if defined(HAVE_SOLVER_K3BIILU2)
		if( _pack == K3BIILU2 )
		{
			bool modified_pattern = ModifiedPattern;
			//~ if( A.comm != comm ) throw DifferentCommunicatorInSolver;
			if( matrix_data == NULL )
			{ 
				MatrixInitDataK3biilu2(&matrix_data,A.GetCommunicator(),A.GetName().c_str());
				modified_pattern = true;
			}
			if( modified_pattern )
			{
				global_size = local_size = A.Size();
				
				MatrixDestroyDataK3biilu2(&matrix_data);
				MatrixInitDataK3biilu2(&matrix_data,A.GetCommunicator(),A.GetName().c_str());
				INMOST_DATA_ENUM_TYPE nnz = 0, k = 0, q = 1, shift = 0;
				int nproc, myid;
#if defined(USE_MPI)
				MPI_Comm_size(A.GetCommunicator(),&nproc);
				MPI_Comm_rank(A.GetCommunicator(),&myid);
#else
				nproc = 1;
				myid = 0;
#endif
				int * ibl = (int *)malloc(sizeof(int)*(nproc+1));
				ibl[0] = 0;
				int n = A.Size();
#if defined(USE_MPI)
				INMOST_DATA_ENUM_TYPE mbeg,mend;
				A.GetInterval(mbeg,mend);
				n = mend - mbeg;
				int block_end = mend;
				MPI_Allgather(&block_end,1,MPI_INT,&ibl[1],1,MPI_INT,A.GetCommunicator());
#else
				ibl[1] = n;
#endif
				int * ia = (int *)malloc(sizeof(int)*(n+1));
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) nnz += it->Size();
				int * ja = (int *)malloc(sizeof(int)*nnz);
				double * values = (double *) malloc(sizeof(double)*nnz);
				ia[0] = shift;
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++)
				{
					for(Sparse::Row::iterator jt = it->Begin(); jt != it->End(); jt++)
					{
						ja[k] = jt->first + 0;
						values[k] = jt->second;
		                                //std::cout<<"# q="<<q<<" k="<<k<<" ja="<<ja[k]<<" a="<<values[k]<<std::endl;//db!
						k++;
					}
					shift += it->Size();
					ia[q++] = shift;
				}
				MatrixFillK3biilu2(matrix_data,n,nproc,ibl,ia,ja,values);
			}
			else
			{
				INMOST_DATA_ENUM_TYPE nnz = 0, k = 0;
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) nnz += it->Size();
				double * values = (double *)malloc(sizeof(double)*nnz);
				k = 0;
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++)
					for(Sparse::Row::iterator jt = it->Begin(); jt != it->End(); jt++)
						values[k++] = jt->second;
				MatrixFillValuesK3biilu2(matrix_data,values);
			}
			MatrixFinalizeK3biilu2(matrix_data);
			SolverSetMatrixK3biilu2(solver_data,matrix_data,modified_pattern,OldPreconditioner);
			ok = true;
		}
#endif
#if defined(USE_SOLVER_SUPERLU)
		if( _pack == SUPERLU ) //just serial SuperLU
		{
			//check that the run is serial!
			int * ia, * ja, nnz = 0;
			double * a;
			int mbeg = A.GetFirstIndex();
			int mend = A.GetLastIndex();
			int size = 0;
			int * remap = new int[mend-mbeg];
			for(int k = 0; k < mend-mbeg; ++k) 
			{
				Sparse::Row & r = A[k+mbeg];
				if( r.Size() )
				{
					double nrm = 0;
					for(int l = 0; l < r.Size(); ++l)
						nrm += r.GetValue(l)*r.GetValue(l);
					if( nrm )
					{
						std::sort(r.Begin(),r.End());
						nnz += r.Size();
						remap[k] = size;
						size++;
					}
					else remap[k] = -1;
				}
				else remap[k] = -1;
			}
			ia = (int *)malloc(sizeof(int)*(size+1));
			ja = (int *)malloc(sizeof(int)*nnz);
			a = (double *)malloc(sizeof(double)*nnz);
			int q = 0, f = 0;
			ia[0] = 0;
			for(int k = 0; k < mend-mbeg; ++k) if( remap[k] != -1 )
			{
				Sparse::Row & r = A[k+mbeg];
				for(int l = 0; l < r.Size(); ++l)
				{
					if( remap[r.GetIndex(l)-mbeg] != -1 )
					{
						ja[q] = remap[r.GetIndex(l)-mbeg];
						a[q] = r.GetValue(l);
						++q;
					}
					else //if( fabs(a[q]) > 1.0e-9 )
						std::cout << "Matrix has connections to singular rows" << std::endl;
				}
				ia[f+1] = q;
				f++;
			}
			MatrixFillSuperLU(solver_data,size,nnz,ia,ja,a,remap);
			//arrays are freed inside SuperLU
			//delete [] ia;
			//delete [] ja;
			//delete [] a;
			matrix_data = solver_data; //to avoid exception in Solve()
			ok = true;
		}
#endif
#if defined(USE_SOLVER_TRILINOS)
		if( _pack == Trilinos_Aztec || 
			_pack == Trilinos_Belos ||
			_pack == Trilinos_ML ||
			_pack == Trilinos_Ifpack)
		{
#if defined(USE_MPI)
			Epetra_MpiComm Comm(A.GetCommunicator());
#else
			Epetra_SerialComm Comm;
#endif
			//std::cout << Comm.MyPID() << " " << __FUNCTION__ << ":" << __LINE__ << std::endl;
			//std::cout << Comm.MyPID() << " " << "Check modified pattern" << std::endl;
			bool modified_pattern = ModifiedPattern;
			INMOST_DATA_ENUM_TYPE mbeg,mend;
			A.GetInterval(mbeg,mend);
			//std::cout << Comm.MyPID() << " " << "Get interval " << mbeg << ":" << mend << std::endl;
			bool refill = true;
			if( matrix_data == NULL || modified_pattern )
			{
				//std::cout << Comm.MyPID() << " " << "Have to reallocate matrix" << std::endl;
				if( matrix_data != NULL )
				{
					delete static_cast<Epetra_CrsMatrix *>(matrix_data);
					matrix_data = NULL;
				}
				local_size = A.Size();
#if defined(USE_MPI)
				MPI_Allreduce(&local_size,&global_size,1,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,comm);
#else
				global_size = local_size;
#endif
				//std::cout << Comm.MyPID() << " " << "local size " << local_size << " global size " << global_size << std::endl;
				int k;
				int imbeg = static_cast<int>(mbeg);
				int imend = static_cast<int>(mend);
				int * temporary = new int[imend-imbeg];
				for(k = imbeg; k < imend; ++k) temporary[k-imbeg] = k;
				//std::cout << Comm.MyPID() << " " << "Creating Epetra_Map" << std::endl;
				Epetra_Map Map(static_cast<int>(global_size),static_cast<int>(local_size),temporary,0,Comm);
				k = 0;
				for(k = imbeg; k < imend; ++k) temporary[k-imbeg] = A[k].Size();
				//std::cout << Comm.MyPID() << " " << "Creating Epetra_CrsMatrix" << std::endl;
				Epetra_CrsMatrix * Matrix = new Epetra_CrsMatrix(Copy,Map,temporary,true);
				matrix_data = static_cast<void *>(Matrix);
				delete [] temporary;
				refill = false;
			}
			{
				assert(matrix_data);
				Epetra_CrsMatrix * Matrix = static_cast<Epetra_CrsMatrix *>(matrix_data);
				unsigned max = 0;
				for(Sparse::Matrix::iterator it = A.Begin(); it != A.End(); ++it)
					if( it->Size() > max ) max = it->Size();
				int * col_positions = new int[max];
				double * col_values = new double[max];
				for(INMOST_DATA_ENUM_TYPE k = mbeg; k < mend; ++k)
				{
					INMOST_DATA_ENUM_TYPE m = 0;
					for(INMOST_DATA_ENUM_TYPE i = 0; i < A[k].Size(); i++)
					{
						col_positions[m] = A[k].GetIndex(i);
						col_values[m] = A[k].GetValue(i);
						m++;
					}
					int ret;
					if( refill )
						ret = Matrix->ReplaceGlobalValues(k,m,col_values,col_positions);
					else
						ret = Matrix->InsertGlobalValues(k,m,col_values,col_positions);
					if( ret != 0 )
						std::cout << "Problem reported by Trilinos! return code " << ret << std::endl;
				}
				delete [] col_positions;
				delete [] col_values;
				Matrix->FillComplete();
				//std::stringstream str;
				//str << "epetra" << Comm.MyPID() << ".mat";
				//std::fstream file(str.str(),std::ios::out);
				//file << *Matrix;
				//file.close();
				assert(solver_data);
				Epetra_LinearProblem * problem = &static_cast<std::pair<std::string,Epetra_LinearProblem> *>(solver_data)->second;
				problem->SetOperator(Matrix);
			}
			ok = true;
		}
#endif
		if (_pack == INNER_ILU2 || _pack == INNER_DDPQILUC || _pack == INNER_MPTILUC || _pack == INNER_MPTILU2)
		{
			
			Sparse::Matrix * mat = new Sparse::Matrix(A);
			info.PrepareMatrix(*mat, additive_schwartz_overlap);
			IterativeMethod * sol = (IterativeMethod *)solver_data;
			sol->ReplaceMAT(*mat);
			if( matrix_data != NULL ) delete (Sparse::Matrix *)matrix_data;
			matrix_data = (void *)mat;

			sol->RealParameter(":tau") = preconditioner_drop_tolerance;
			sol->RealParameter(":tau2") = preconditioner_reuse_tolerance;
			sol->EnumParameter(":scale_iters") = preconditioner_rescale_iterations;

			if( _pack == INNER_DDPQILUC )
			{
				sol->RealParameter(":ddpq_tau") = preconditioner_ddpq_tolerance;
				sol->EnumParameter(":reorder_nnz") = preconditioner_reorder_nonzero;
				sol->EnumParameter(":estimator") = preconditioner_condition_estimation;
				sol->EnumParameter(":ddpq_tau_adapt") = preconditioner_adapt_ddpq_tolerance;
			}
			else if( _pack == INNER_MPTILUC )
			{
				sol->EnumParameter(":estimator") = preconditioner_condition_estimation;
			}
			else if( _pack == INNER_ILU2 ) sol->EnumParameter(":fill") = static_cast<INMOST_DATA_ENUM_TYPE>(preconditioner_fill_level);

			if( sizeof(KSOLVER) == sizeof(BCGSL_solver) )
				sol->EnumParameter("levels") = solver_gmres_substeps;

			if (!sol->isInitialized()) 
			{
				sol->Initialize();
			}
			ok = true;

		}
		if(!ok) throw NotImplemented;
	}

#if defined(ACCELERATED_CONDEST)
  INMOST_DATA_REAL_TYPE Solver::Condest(INMOST_DATA_REAL_TYPE tol, INMOST_DATA_ENUM_TYPE maxiter)
  {
    if( matrix_data == NULL ) throw MatrixNotSetInSolver;
    if( _pack == INNER_ILU2 || _pack == INNER_DDPQILUC || _pack == INNER_MPTILU2 || _pack == INNER_MPTILUC )
    {
      Sparse::Matrix & A = *(Sparse::Matrix *)matrix_data;
      //IterativeMethod & sol = *(IterativeMethod *)solver_data;
      INMOST_DATA_ENUM_TYPE lbeg, lend, l, iter;
      INMOST_DATA_REAL_TYPE norm, sum[2], norm_prev, lambda_min, lambda_max;
      bool diverged_max = false, diverged_min = false;
      info.GetLocalRegion(info.GetRank(),lbeg,lend);
      Sparse::Vector v, Av;
      info.PrepareVector(v);
      info.PrepareVector(Av);
      //Set v to random
      norm = 0;
      for(l = lbeg; l < lend; ++l)
      {
        v[l] = rand()/static_cast<INMOST_DATA_REAL_TYPE>(RAND_MAX);
        norm += v[l]*v[l];
      }
      info.Integrate(&norm,1);
      norm_prev = 0.0;
      norm = sqrt(norm);
      for(l = lbeg; l < lend; ++l) v[l] /= norm;
      //Compute maximum eigenvalue
      iter = 0;
      while( fabs((norm-norm_prev)/norm) > tol && iter < maxiter)
      {
        info.Update(v);
#if defined(USE_OMP)
#pragma omp parallel
#endif
        A.MatVec(1.0,v,0.0,Av);
        norm_prev = norm;
        sum[0] = sum[1] = 0.0;
        for(l = lbeg; l < lend; ++l) 
        {
          sum[0] += v[l]*Av[l];
          sum[1] += v[l]*v[l];
        }
        info.Integrate(sum,2);
        norm = fabs(sum[0])/sum[1];
        for(l = lbeg; l < lend; ++l)  v[l]  = Av[l]/norm;
#if defined(PRINT_CONDEST)
        std::cout << "iteration " << iter << " norm " << norm << std::endl;
#endif
        iter++;
      }
#if defined(PRINT_CONDEST)
      std::cout << "lambda_max " << norm << std::endl;
#endif
      if( iter == maxiter ) 
      {
        diverged_max = true;
        std::cout << "Max not converged" << std::endl;
      }
      lambda_max = norm;
      //Set v to random
      norm = 0;
      for(l = lbeg; l < lend; ++l)
      {
        v[l] = rand()/static_cast<INMOST_DATA_REAL_TYPE>(RAND_MAX);
        norm += v[l]*v[l];
      }
      info.Integrate(&norm,1);
      norm_prev = 0.0;
      norm = sqrt(norm);
      for(l = lbeg; l < lend; ++l) v[l] /= norm;
      //Compute minimal eigenvalue
      iter = 0;
      while( fabs((norm-norm_prev)/norm) > tol && iter < maxiter)
      {
        info.Update(v);
        Solve(v,Av);
        norm_prev = norm;
        sum[0] = sum[1] = 0;
        for(l = lbeg; l < lend; ++l) 
        {
          sum[0] += v[l]*Av[l];
          sum[1] += v[l]*v[l];
        }
        info.Integrate(sum,2);
        norm = fabs(sum[0])/sum[1];
        for(l = lbeg; l < lend; ++l) v[l] = Av[l]/norm;
#if defined(PRINT_CONDEST)
        std::cout << "iteration " << iter << " norm " << norm << "\t\t" << std::endl;
#endif
        iter++;
      }
#if defined(PRINT_CONDEST)
      std::cout << "lambda_min " << 1.0/norm << std::endl;
#endif
      if( iter == maxiter ) 
      {
        diverged_min = true;
        std::cout << "Min not converged" << std::endl;
      }
      lambda_min = 1.0/norm;
      if( diverged_max || diverged_min )
        return 1.0e+100;
#if defined(PRINT_CONDEST)
      std::cout << "Condest: " << lambda_max /lambda_min << std::endl;
#endif
      return lambda_max/lambda_min;
    }
    else 
    {
      throw -1; //other packages are not supported yet
    }
  }
#else
  INMOST_DATA_REAL_TYPE Solver::Condest(INMOST_DATA_REAL_TYPE tol, INMOST_DATA_ENUM_TYPE maxiter)
  {
    if( matrix_data == NULL ) throw MatrixNotSetInSolver;
    if( _pack == INNER_ILU2 || _pack == INNER_DDPQILUC || _pack == INNER_MPTILU2 || _pack == INNER_MPTILUC )
    {
      Sparse::Matrix & A = *(Sparse::Matrix *)matrix_data;
      //IterativeMethod & sol = *(IterativeMethod *)solver_data;
      INMOST_DATA_ENUM_TYPE lbeg, lend, l, iter;
      INMOST_DATA_REAL_TYPE norm, norm_prev, lambda_min, lambda_max;
      bool diverged_max = false, diverged_min = false;
      info.GetLocalRegion(info.GetRank(),lbeg,lend);
      Sparse::Vector v, Av;
      info.PrepareVector(v);
      info.PrepareVector(Av);
      //Set v to random
      norm = 0;
      for(l = lbeg; l < lend; ++l)
      {
        v[l] = rand()/static_cast<INMOST_DATA_REAL_TYPE>(RAND_MAX);
        norm += v[l]*v[l];
      }
      info.Integrate(&norm,1);
      norm_prev = 0.0;
      norm = sqrt(norm);
      for(l = lbeg; l < lend; ++l) v[l] /= norm;
      //Compute maximum eigenvalue
      iter = 0;
      while( fabs(norm-norm_prev)/norm > tol && iter < maxiter)
      {
        info.Update(v);
#if defined(USE_OMP)
#pragma omp parallel
#endif
        A.MatVec(1.0,v,0.0,Av);
        v.Swap(Av);
        norm_prev = norm;
        norm = 0.0;
        for(l = lbeg; l < lend; ++l) norm += v[l]*v[l];
        info.Integrate(&norm,1);
        norm = sqrt(norm);
        for(l = lbeg; l < lend; ++l) v[l] /= norm;
#if defined(PRINT_CONDEST)
        std::cout << "iteration " << iter << " norm " << norm << std::endl;
#endif
        iter++;
      }
#if defined(PRINT_CONDEST)
      std::cout << "lambda_max " << norm << std::endl;
#endif
      if( iter == maxiter ) 
      {
        norm = std::max(norm,norm_prev);
        //diverged_max = true;
        //std::cout << "Max not converged" << std::endl;
      }
      lambda_max = norm;
      //Set v to random
      norm = 0;
      for(l = lbeg; l < lend; ++l)
      {
        v[l] = rand()/static_cast<INMOST_DATA_REAL_TYPE>(RAND_MAX);
        norm += v[l]*v[l];
      }
      info.Integrate(&norm,1);
      norm_prev = 0.0;
      norm = sqrt(norm);
      for(l = lbeg; l < lend; ++l) v[l] /= norm;
      //Compute minimal eigenvalue
      iter = 0;
      while( fabs(norm-norm_prev)/norm > tol && iter < maxiter)
      {
        info.Update(v);
        Solve(v,Av);
        v.Swap(Av);
        norm_prev = norm;
        norm = 0.0;
        for(l = lbeg; l < lend; ++l) norm += v[l]*v[l];
        info.Integrate(&norm,1);
        norm = sqrt(norm);
        for(l = lbeg; l < lend; ++l) v[l] /= norm;
#if defined(PRINT_CONDEST)
        std::cout << "iteration " << iter << " norm " << norm << "\t\t" << std::endl;
#endif
        iter++;
      }
#if defined(PRINT_CONDEST)
      std::cout << "lambda_min " << 1.0/norm << std::endl;
#endif
      if( iter == maxiter ) 
      {
        norm = std::max(norm,norm_prev);
        //diverged_min = true;
        //std::cout << "Min not converged" << std::endl;
      }
      lambda_min = 1.0/norm;
      if( diverged_max || diverged_min )
        return 1.0e+100;
#if defined(PRINT_CONDEST)
      std::cout << "Condest: " << lambda_max /lambda_min << std::endl;
#endif
      return lambda_max/lambda_min;
    }
    else 
    {
      throw -1; //other packages are not supported yet
    }
  }
#endif
	bool Solver::Solve(Sparse::Vector & RHS, Sparse::Vector & SOL)
	{
		//check the data
		if( matrix_data == NULL ) throw MatrixNotSetInSolver;
		if( RHS.GetCommunicator() != comm || SOL.GetCommunicator() != comm ) throw DifferentCommunicatorInSolver;
		INMOST_DATA_ENUM_TYPE vbeg,vend;
		RHS.GetInterval(vbeg,vend);
		if( RHS.Size() != SOL.Size() )
		{
			if( SOL.Size() == 0 )
			{
				SOL.SetInterval(vbeg,vend);
				for(Sparse::Vector::iterator ri = SOL.Begin(); ri != SOL.End(); ++ri) *ri = 0.0;
			}
			else throw InconsistentSizesInSolver;
		}
		//run the solver
#if defined(USE_SOLVER_PETSC)
		if( _pack == PETSc )
		{
			if( rhs_data == NULL )
				VectorInitDataPetsc(&rhs_data,RHS.GetCommunicator(),RHS.GetName().c_str());
			VectorPreallocatePetsc(rhs_data,local_size,global_size);
			
			if( solution_data == NULL )
				VectorInitDataPetsc(&solution_data,SOL.GetCommunicator(),SOL.GetName().c_str());
			VectorPreallocatePetsc(solution_data,local_size,global_size);
		
			
			int * positions = new int[local_size];
			double * values = new double[local_size];
			{
				unsigned k = 0;
				for(Sparse::Vector::iterator it = RHS.Begin(); it != RHS.End(); ++it)
				{
					positions[k] = vbeg+k;
					values[k] = *it;
					k++;
				}
				VectorFillPetsc(rhs_data,local_size,positions,values);
				VectorFinalizePetsc(rhs_data);
				
				k = 0;
				for(Sparse::Vector::iterator it = SOL.Begin(); it != SOL.End(); ++it)
				{
					values[k] = *it;
					k++;
				}
				VectorFillPetsc(solution_data,local_size,positions,values);
				VectorFinalizePetsc(solution_data);
			}
			if( petsc_database_file == "" ) //set parameters if no file is set
			{
				SolverSetTolerancesPetsc(solver_data,relative_tolerance,absolute_tolerance,divergence_tolerance,maximum_iterations);
			}
			bool result = SolverSolvePetsc(solver_data,rhs_data,solution_data);
			if( result )
			{
				VectorLoadPetsc(solution_data,local_size,positions,values);
				unsigned k = 0;
				for(Sparse::Vector::iterator it = SOL.Begin(); it != SOL.End(); ++it)
				{
					*it = values[k];
					k++;
				}
			}
			delete [] positions;
			delete [] values;
			last_resid = SolverResidualNormPetsc(solver_data);
			last_it = SolverIterationNumberPetsc(solver_data);
			return_reason = std::string(SolverConvergedReasonPetsc(solver_data));
			return result;
		}
#endif
#if defined(USE_SOLVER_ANI)
		if( _pack == ANI )
		{
			if( rhs_data == NULL )
				VectorInitDataAni(&rhs_data,RHS.GetCommunicator(),RHS.GetName().c_str());
			VectorPreallocateAni(rhs_data,local_size);
			
			if( solution_data == NULL )
				VectorInitDataAni(&solution_data,SOL.GetCommunicator(),SOL.GetName().c_str());
			VectorPreallocateAni(solution_data,local_size);
			{
				VectorFillAni(rhs_data,&RHS[vbeg]);
				VectorFinalizeAni(rhs_data);
				
				VectorFillAni(solution_data,&SOL[vbeg]);
				VectorFinalizeAni(solution_data);
			}
			bool result = SolverSolveAni(solver_data,rhs_data,solution_data);
			if( result ) VectorLoadAni(solution_data,&SOL[vbeg]);
			last_resid = SolverResidualNormAni(solver_data);
			last_it = SolverIterationNumberAni(solver_data);
			return_reason = "Unspecified for ANI";
			return result;
		}
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
		if( _pack == FCBIILU2 )
		{
			if( rhs_data == NULL )
				VectorInitDataFcbiilu2(&rhs_data,RHS.GetCommunicator(),RHS.GetName().c_str());
			VectorPreallocateFcbiilu2(rhs_data,local_size);
			
			if( solution_data == NULL )
				VectorInitDataFcbiilu2(&solution_data,SOL.GetCommunicator(),SOL.GetName().c_str());
			VectorPreallocateFcbiilu2(solution_data,local_size);
			{
				VectorFillFcbiilu2(rhs_data,&RHS[vbeg]);
				VectorFinalizeFcbiilu2(rhs_data);
				
				VectorFillFcbiilu2(solution_data,&SOL[vbeg]);
				VectorFinalizeFcbiilu2(solution_data);
			}
			bool result = SolverSolveFcbiilu2(solver_data,rhs_data,solution_data);
			if( result ) VectorLoadFcbiilu2(solution_data,&SOL[vbeg]);
			last_resid = SolverResidualNormFcbiilu2(solver_data);
			last_it = SolverIterationNumberFcbiilu2(solver_data);
			//return_reason = std::string(SolverConvergedReasonFcbiilu2(solver_data));
			return_reason = "Unspecified for FCBIILU2";
			return result;
		}
#endif
#if defined(HAVE_SOLVER_K3BIILU2)
		if( _pack == K3BIILU2 )
		{
			if( rhs_data == NULL )
				VectorInitDataK3biilu2(&rhs_data,RHS.GetCommunicator(),RHS.GetName().c_str());
			VectorPreallocateK3biilu2(rhs_data,local_size);
			
			if( solution_data == NULL )
				VectorInitDataK3biilu2(&solution_data,SOL.GetCommunicator(),SOL.GetName().c_str());
			VectorPreallocateK3biilu2(solution_data,local_size);
			{
				VectorFillK3biilu2(rhs_data,&RHS[vbeg]);
				VectorFinalizeK3biilu2(rhs_data);
				
				VectorFillK3biilu2(solution_data,&SOL[vbeg]);
				VectorFinalizeK3biilu2(solution_data);
			}
			bool result = SolverSolveK3biilu2(solver_data,rhs_data,solution_data);
			if( result ) VectorLoadK3biilu2(solution_data,&SOL[vbeg]);
			last_resid = SolverResidualNormK3biilu2(solver_data);
			last_it = SolverIterationNumberK3biilu2(solver_data);
			//return_reason = std::string(SolverConvergedReasonK3biilu2(solver_data));
			return_reason = "Unspecified for K3BIILU2";
			return result;
		}
#endif
#if defined(USE_SOLVER_TRILINOS)
		if(_pack == Trilinos_Aztec || 
		   _pack == Trilinos_ML ||
		   _pack == Trilinos_Ifpack)
		{
			assert(matrix_data);
			std::string name = static_cast<std::pair<std::string,Epetra_LinearProblem> *>(solver_data)->first;
			Epetra_LinearProblem * problem = &static_cast<std::pair<std::string,Epetra_LinearProblem> *>(solver_data)->second;
			Epetra_CrsMatrix * Matrix = static_cast<Epetra_CrsMatrix *>(matrix_data);
			Epetra_Vector VectorRHS(View,Matrix->Map(),&*RHS.Begin());
			Epetra_Vector VectorSOL(View,Matrix->Map(),&*SOL.Begin());
			problem->SetRHS(&VectorRHS);
			problem->SetLHS(&VectorSOL);
			std::string specific_database_file = "";
			if( _pack == Trilinos_Aztec )
				specific_database_file = trilinos_aztec_database_file;
			else if( _pack == Trilinos_ML )
				specific_database_file = trilinos_ml_database_file;
			else if( _pack == Trilinos_Ifpack )
				specific_database_file = trilinos_ifpack_database_file;

			bool have_params = specific_database_file != "";
			const Teuchos::RCP<Teuchos::ParameterList> top_level_params = Teuchos::createParameterList();
			Teuchos::ParameterList local_list;
			if( have_params ) 
			{
				Teuchos::updateParametersFromXmlFileAndBroadcast(specific_database_file,top_level_params.ptr(),Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(RHS.GetCommunicator())));
				//top_level_params->print(std::cout,0,true,true);
				if( !top_level_params->isSublist(name) )
					have_params = false;
				else
				{
					local_list = top_level_params->sublist(name);
					//std::cout << "Got local list for " << name << std::endl;
				}
			}
			

			AztecOO AztecSolver(*problem);
			
			
			if( have_params && local_list.isSublist("AztecOO"))
			{
				Teuchos::ParameterList AztecOOParams = local_list.sublist("AztecOO");
				if( AztecOOParams.isParameter("Max Iterations") )
				{
					maximum_iterations = AztecOOParams.get<int>("Max Iterations");
					//std::cout << "Got maximum iterations " << maximum_iterations << std::endl;
				}
				if( AztecOOParams.isParameter("Tolerance") )
				{
					relative_tolerance = AztecOOParams.get<double>("Tolerance");
					//std::cout << "Got tolerance " << relative_tolerance << std::endl;
				}
				if( AztecOOParams.isSublist("AztecOO Settings") )
				{
					AztecSolver.SetParameters(AztecOOParams.sublist("AztecOO Settings"));
					//std::cout << "Submit Aztec settings" << std::endl;
				}
			}
			else
			{
				AztecSolver.SetAztecOption(AZ_diagnostics,AZ_none);
				AztecSolver.SetAztecOption(AZ_output,AZ_none);
				AztecSolver.SetAztecOption(AZ_solver,AZ_bicgstab);
				AztecSolver.SetAztecOption(AZ_overlap,additive_schwartz_overlap);
			}
			
			void * precinfo = NULL;
			if( _pack == Trilinos_Aztec )
			{
				if( !have_params )
				{
					AztecSolver.SetAztecParam(AZ_drop,preconditioner_drop_tolerance);
					AztecSolver.SetAztecParam(AZ_ilut_fill,preconditioner_fill_level);
					//AztecSolver.SetAztecOption(AZ_solver,AZ_tfqmr);
					//AztecSolver.SetAztecOption(AZ_precond,AZ_Neumann);
					//AztecSolver.SetAztecOption(AZ_poly_ord,3);
				}
				else
				{
					//should be set automatically
				}
			}
			else if( _pack == Trilinos_ML )
			{
				Teuchos::ParameterList List;
				if( have_params && local_list.isSublist("ML") && local_list.sublist("ML").isSublist("ML Settings") )
				{
					List = local_list.sublist("ML").sublist("ML Settings");
				}
				else
				{
					ML_Epetra::SetDefaults("SA",List);
					List.set("max levels",6);
					List.set("increasing or decreasing","decreasing");
					//List.set("aggreagation: type", "MIS");
					//List.set("coarse: type", "Amesos-KLU");
				}
				ML_Epetra::MultiLevelPreconditioner * Prec = new ML_Epetra::MultiLevelPreconditioner(*Matrix,List,true);
				AztecSolver.SetPrecOperator(Prec);
				precinfo = Prec;
			}
			else if( _pack == Trilinos_Ifpack )
			{
				Ifpack * Factory = new Ifpack();
				Ifpack_Preconditioner * Prec;
				std::string PrecType = "ILU";
				if( have_params && local_list.isSublist("Ifpack") )
				{
					Teuchos::ParameterList ifpacklist = local_list.sublist("Ifpack");
					if( ifpacklist.isParameter("Prec Type") )
					{
						PrecType = ifpacklist.get<std::string>("Prec Type");
						//std::cout << "Got preconditioner type " << PrecType << std::endl;
					}
					if( ifpacklist.isParameter("Overlap") )
					{
						additive_schwartz_overlap = ifpacklist.get<int>("Overlap");
						//std::cout << "Got overlap " << additive_schwartz_overlap << std::endl;
					}
				}
				Prec = Factory->Create(PrecType,Matrix,additive_schwartz_overlap);
				Teuchos::ParameterList List;
				if( have_params && local_list.isSublist("Ifpack") && local_list.sublist("Ifpack").isSublist("Ifpack Settings") )
				{
					List = local_list.sublist("Ifpack").sublist("Ifpack Settings");
					//std::cout << "Submit settings to ifpack" << std::endl;
				}
				else
					List.set("fact: level-of-fill",static_cast<int>(preconditioner_fill_level));
				Prec->SetParameters(List);
				Prec->Initialize();
				Prec->Compute();
				AztecSolver.SetPrecOperator(Prec);
				precinfo = Factory;
			}
			AztecSolver.Iterate(maximum_iterations,relative_tolerance);
			if( _pack == Trilinos_ML )
				delete static_cast<ML_Epetra::MultiLevelPreconditioner *>(precinfo);
			else if( _pack == Trilinos_Ifpack )
				delete static_cast<Ifpack *>(precinfo);
			const double * stats = AztecSolver.GetAztecStatus();
			bool ret = true;
			switch(static_cast<int>(stats[AZ_why]))
			{
			case AZ_normal:
				return_reason = "User requested convergence criteria is satisfied.";
				break;
			case AZ_param:
				return_reason = "User requested option is not availible.";
				ret = false;
				break;
			case AZ_breakdown:
				return_reason = "Numerical breakdown occurred.";
				ret = false;
				break;
			case AZ_loss:
				return_reason = "Numerical loss precision occurred.";
        ret = false;
				break;
			case AZ_ill_cond:
				return_reason = "The Hessenberg matrix within GMRES is illconditioned. "
								"This could be caused by a number "
								"of reasons. For example, the preconditioning "
								"matrix could be nearly singular due to an unstable "
								"factorization (note: pivoting is not implemented "
								"in any of the incomplete factorizations). "
								"Ill-conditioned Hessenberg matrices could also "
								"arise from a singular application matrix. In this "
								"case, GMRES tries to compute a least-squares "
								"solution.";
				ret = false;
				break;
			case AZ_maxits:
				return_reason = "Maximum iterations taken without convergence.";
				ret = false;
				break;
			default:
				{
					std::stringstream str;
					str << "reason code " << static_cast<int>(stats[AZ_why]) << " was not specified in manual by the time, reffer to Trilinos manual.";
					return_reason = str.str();
				}
				break;
			}
			last_it = AztecSolver.NumIters();
			last_resid = AztecSolver.TrueResidual();
			return ret;
		}
		if(_pack == Trilinos_Belos )
		{
			assert(matrix_data);
			std::string name = static_cast<std::pair<std::string,Epetra_LinearProblem> *>(solver_data)->first;
			Epetra_LinearProblem * problem = &static_cast<std::pair<std::string,Epetra_LinearProblem> *>(solver_data)->second;
			Epetra_CrsMatrix * Matrix = static_cast<Epetra_CrsMatrix *>(matrix_data);
			Epetra_Vector VectorRHS(View,Matrix->Map(),&*RHS.Begin());
			Epetra_Vector VectorSOL(View,Matrix->Map(),&*SOL.Begin());
			problem->SetRHS(&VectorRHS);
			problem->SetLHS(&VectorSOL);

			bool have_params = trilinos_belos_database_file != "";
			const Teuchos::RCP<Teuchos::ParameterList> top_level_params = Teuchos::createParameterList();
			if( have_params ) Teuchos::updateParametersFromXmlFileAndBroadcast(trilinos_belos_database_file,top_level_params.ptr(),Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(RHS.GetCommunicator())));
			
			Teuchos::RCP<Teuchos::ParameterList> List = Teuchos::rcp(new Teuchos::ParameterList);
			
			if( have_params && top_level_params->isSublist(name) && top_level_params->sublist(name).isSublist("Belos") )
				*List = top_level_params->sublist(name).sublist("Belos");
			else
			{
				List->set("Num Blocks",100);
				List->set("Block Size",1);
				List->set("Maximum Iterations",static_cast<int>(maximum_iterations));
				List->set("Maximum Restarts",20);
				List->set("Convergence Tolerance",static_cast<double>(relative_tolerance));
			}
			//int verbosity = Belos::Warnings + Belos::Errors + Belos::StatusTestDetails + Belos::Debug + Belos::FinalSummary + Belos::TimingDetails;
			//List->set("Verbosity",verbosity);
			Teuchos::RCP<Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> > Problem = 
				Teuchos::rcp(new Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>);
			Problem->setLHS(Teuchos::rcp_implicit_cast<Epetra_MultiVector,Epetra_Vector>(Teuchos::rcp(&VectorSOL,false)));
			Problem->setRHS(Teuchos::rcp_implicit_cast<Epetra_MultiVector,Epetra_Vector>(Teuchos::rcp(&VectorRHS,false)));
			Problem->setOperator(Teuchos::rcp_implicit_cast<Epetra_Operator,Epetra_CrsMatrix>(Teuchos::rcp(Matrix,false)));

			/*
			//This don't work
				Teuchos::ParameterList MLList;
				ML_Epetra::SetDefaults("SA",MLList);
				MLList.set("max levels",6);
				MLList.set("increasing or decreasing","decreasing");
				//List.set("aggreagation: type", "MIS");
				//List.set("coarse: type", "Amesos-KLU");
				ML_Epetra::MultiLevelPreconditioner * Prec = new ML_Epetra::MultiLevelPreconditioner(*Matrix,MLList,true);
				Prec->ComputePreconditioner();
				Problem->setLeftPrec(Teuchos::rcp_implicit_cast<Epetra_Operator,ML_Epetra::MultiLevelPreconditioner>(Teuchos::rcp(Prec,false)));
			*/

			Problem->setProblem();
			
			Teuchos::RCP< Belos::SolverManager<double,Epetra_MultiVector,Epetra_Operator> > BelosSolver =
			//	Teuchos::rcp(new Belos::BlockCgSolMgr<double,Epetra_MultiVector,Epetra_Operator >);
				Teuchos::rcp(new Belos::BlockGmresSolMgr<double,Epetra_MultiVector,Epetra_Operator >);
			//	Teuchos::rcp(new Belos::PseudoBlockGmresSolMgr<double,Epetra_MultiVector,Epetra_Operator >);
			BelosSolver->setParameters(List);
			BelosSolver->setProblem(Problem);
			Belos::ReturnType ret = BelosSolver->solve();
			last_it = BelosSolver->getNumIters();
			last_resid = BelosSolver->achievedTol();
			
			//delete Prec;

			if( ret == Belos::Converged )
			{
				return_reason = "Converged";
				return true;
			}
			else 
			{
				std::stringstream reason_stream;
				reason_stream << "Diverged ";
				if( BelosSolver->isLOADetected() )
					reason_stream << "loss of accuracy detected ";
				if( BelosSolver->getNumIters() >= 1000 )
					reason_stream << "maximum iteration number reached ";
				if( BelosSolver->achievedTol() > 1.0e+6 )
					reason_stream << "divergence tolerance reached ";
				return_reason = reason_stream.str();
				return false;
			}
			
		}
#endif
#if defined(USE_SOLVER_SUPERLU)
		if(_pack == SUPERLU )
		{
			int size = MatrixSizeSuperLU(solver_data);
			double * inout = new double[size];
			int * remap = MatrixRemapArraySuperLU(solver_data);
			int mbeg = RHS.GetFirstIndex(), mend = RHS.GetLastIndex();
			for(int k = 0; k < mend - mbeg; ++k) if( remap[k] != -1 ) inout[remap[k]] = RHS[k+mbeg];
			bool ret = SolverSolveSuperLU(solver_data,inout);
			for(int k = 0; k < mend - mbeg; ++k) if( remap[k] != -1 ) SOL[k+mbeg] = inout[remap[k]];
			delete [] inout;
			last_it = 1;
			last_resid = 0;
			return_reason = SolverConvergedReasonSuperLU(solver_data);
			return ret;
		}
#endif
		if(_pack == INNER_ILU2 || _pack == INNER_DDPQILUC || _pack == INNER_MPTILUC || _pack == INNER_MPTILU2)
		{
			IterativeMethod * sol = static_cast<IterativeMethod *>(solver_data);
			
			sol->EnumParameter("maxits") = maximum_iterations;
			sol->RealParameter("rtol") = relative_tolerance;
			sol->RealParameter("atol") = absolute_tolerance;
			sol->RealParameter("divtol") = divergence_tolerance;
			
			bool ret = sol->Solve(RHS,SOL);
			last_it = sol->GetIterations();
			last_resid = sol->GetResidual();
			return_reason = sol->GetReason();
			return ret;
		}
		throw NotImplemented;
	}

	std::string Solver::GetReason()
	{
		return return_reason;
	}
	
	INMOST_DATA_ENUM_TYPE Solver::Iterations()
	{
		return last_it;
	}
	INMOST_DATA_REAL_TYPE Solver::Residual()
	{
		return last_resid;
	}

  INMOST_DATA_REAL_TYPE Solver::GetConditionNumberL()
  {
    if(_pack == INNER_MPTILUC )
    {
      return static_cast<IterativeMethod *>(solver_data)->RealParameter(":condition_number_L");
    }
    else return 0.0;
  }
  INMOST_DATA_REAL_TYPE Solver::GetConditionNumberU()
  {
    if(_pack == INNER_MPTILUC )
    {
      return static_cast<IterativeMethod *>(solver_data)->RealParameter(":condition_number_U");
    }
    else return 0.0;
  }

}
#endif

