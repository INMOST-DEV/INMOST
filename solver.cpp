#include "inmost_solver.h"
#if defined(USE_SOLVER)
#include "solver_petsc.h"
#include "solver_ani.h"
#include "solver_ilu2.hpp"
#include "solver_ddpqiluc2.hpp"
#include "solver_bcgsl.hpp"
#include <fstream>
#include <sstream>
#include <bitset>
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

#define KSOLVER BCGSL_solver
//#define KSOLVER BCGS_solver


namespace INMOST
{
	static std::string petsc_database_file = "";
	static std::string trilinos_aztec_database_file = "";
	static std::string trilinos_ifpack_database_file = "";
	static std::string trilinos_belos_database_file = "";
	static std::string trilinos_ml_database_file = "";
	static std::string ani_database_file = "";


#if defined(USE_MPI)
	INMOST_MPI_Type Solver::RowEntryType = MPI_DATATYPE_NULL;
#else
	INMOST_MPI_Type Solver::RowEntryType = 0;
#endif

#define GUARD_MPI(x) {ierr = x; if( ierr != MPI_SUCCESS ) {std::cout << #x << " not successfull " << std::endl; MPI_Abort(comm,-1000);}}
#define HASH_TABLE_SIZE 2048

	bool Solver::is_initialized = false;
	bool Solver::is_finalized = false;

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
		int ierr = 0;
		dynarray<INMOST_DATA_REAL_TYPE,64> temp(num);
		memcpy(temp.data(),inout,sizeof(INMOST_DATA_REAL_TYPE)*num);
		GUARD_MPI(MPI_Allreduce(temp.data(),inout,num,INMOST_MPI_DATA_REAL_TYPE,MPI_SUM,comm));
#else
		(void) inout;
		(void) num;
#endif
	}

	void Solver::OrderInfo::PrepareMatrix(Matrix & m, INMOST_DATA_ENUM_TYPE overlap)
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
		std::vector< Row::entry > send_row_data, recv_row_data;
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
				Solver::Row & r = m[k];
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
					for(Solver::Row::iterator jt = m[k].Begin(); jt != m[k].End(); ++jt)
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
					GUARD_MPI(MPI_Irecv(&recv_row_data[f],local_size,Solver::GetRowEntryType(),vector_exchange_recv[j],4*size+vector_exchange_recv[j],comm,&requests[k]));
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
					GUARD_MPI(MPI_Isend(&send_row_data[f],local_size,Solver::GetRowEntryType(),vector_exchange_send[j],4*size+rank,comm,&requests[k+vector_exchange_recv[0]]));
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
						m[global_to_local[vector_exchange_recv[j+2+r]]] = Solver::Row(&recv_row_data[f],&recv_row_data[f]+recv_row_sizes[q+r]);
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

	void Solver::OrderInfo::RestoreMatrix(Solver::Matrix & m)
	{
		//restore matrix size
		m.SetInterval(initial_matrix_begin,initial_matrix_end);
		//restore indexes
		for(Solver::Matrix::iterator it = m.Begin(); it != m.End(); ++it)
			for(Solver::Row::iterator jt = it->Begin(); jt != it->End(); ++jt)
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
	
	void Solver::OrderInfo::PrepareVector(Vector & v) const
	{
		if( !have_matrix ) throw PrepareMatrixFirst;
		v.SetInterval(local_vector_begin,local_vector_end);
		v.isParallel() = true;
	}
	
	void Solver::OrderInfo::RestoreVector(Vector & v) const
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
	
	Solver::Vector::Vector(std::string _name, INMOST_DATA_ENUM_TYPE start, INMOST_DATA_ENUM_TYPE end, INMOST_MPI_Comm _comm) :data(start,end) 
	{
		comm = _comm;
		name = _name;
		is_parallel = false;
	}
	
	Solver::Vector::Vector(const Vector & other) : data(other.data) 
	{
		comm = other.comm;
		name = other.name;
		is_parallel = other.is_parallel;
	}
	
	Solver::Vector & Solver::Vector::operator =(Vector const & other) 
	{
		comm = other.comm;
		data = other.data; 
		name = other.name; 
		is_parallel = other.is_parallel;
		return *this;
	}
	
	Solver::Vector::~Vector() 
	{
	}
	
	void Solver::OrderInfo::Update(Vector & x)
	{
		//std::cout << __FUNCTION__ << " start" << std::endl;
#if defined(USE_MPI)
		if( GetSize() == 1 ) return;
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
#else
		(void) x;
#endif
		//std::cout << __FUNCTION__ << " end" << std::endl;
	}
	void Solver::OrderInfo::Accumulate(Vector & x)
	{
		//std::cout << __FUNCTION__ << " start" << std::endl;
#if defined(USE_MPI)
		if( GetSize() == 1 ) return;
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
#else
		(void) x;
#endif
		//std::cout << __FUNCTION__ << " end" << std::endl;
	}

	INMOST_DATA_REAL_TYPE Solver::OrderInfo::ScalarProd(Vector const & left, Vector const & right, INMOST_DATA_ENUM_TYPE index_begin, INMOST_DATA_ENUM_TYPE index_end) const
	{
		INMOST_DATA_REAL_TYPE ret = 0;
		for(INMOST_DATA_ENUM_TYPE i = index_begin; i < index_end; ++i)
			ret += left[i]*right[i];
		Integrate(&ret,1);
		return ret;
	}

	
	INMOST_DATA_REAL_TYPE   Solver::Row::RowVec(Solver::Vector & x) const
	{
		INMOST_DATA_REAL_TYPE ret = 0;
		INMOST_DATA_ENUM_TYPE end = Size(),i;
		for(i = 0; i < end; i++) ret = ret + x[GetIndex(i)]*GetValue(i);
		return ret;
	}
	
	void Solver::Matrix::MatVec(INMOST_DATA_REAL_TYPE alpha, Solver::Vector & x, INMOST_DATA_REAL_TYPE beta, Solver::Vector & out) const //y = alpha*A*x + beta * y
	{
		INMOST_DATA_ENUM_TYPE mbeg, mend;
		INMOST_DATA_INTEGER_TYPE ind, imbeg, imend;
		if( out.Empty() )
		{
			INMOST_DATA_ENUM_TYPE vbeg,vend;
			GetInterval(vbeg,vend);
			out.SetInterval(vbeg,vend);
		}
		//CHECK SOMEHOW FOR DEBUG THAT PROVIDED VECTORS ARE OK
		//~ assert(GetFirstIndex() == out.GetFirstIndex());
		//~ assert(Size() == out.Size());
		GetInterval(mbeg,mend);
		imbeg = mbeg;
		imend = mend;
#if defined(USE_OMP)
#pragma omp for private(ind)
#endif
		for(ind = imbeg; ind < imend; ++ind) //iterate rows of matrix
			out[ind] = beta * out[ind] + alpha * (*this)[ind].RowVec(x);
		// outer procedure should update Sout vector, if needed
	}
	
	Solver::Matrix::Matrix(std::string _name, INMOST_DATA_ENUM_TYPE start, INMOST_DATA_ENUM_TYPE end, INMOST_MPI_Comm _comm) 
	:data(start,end) 
	{
		is_parallel = false;
		comm = _comm;
		SetInterval(start,end);
		name = _name;
	}
	
	Solver::Matrix::Matrix(const Matrix & other) :data(other.data) 
	{
		comm = other.comm;
		name = other.name;
	}
	
	Solver::Matrix & Solver::Matrix::operator =(Matrix const & other) 
	{
		comm = other.comm;
		data = other.data; 
		name = other.name; return *this;
	}
	
	Solver::Matrix::~Matrix() 
	{
	}
	
	
	void      Solver::Matrix::MoveRows(INMOST_DATA_ENUM_TYPE from, INMOST_DATA_ENUM_TYPE to, INMOST_DATA_ENUM_TYPE size)
	{
		INMOST_DATA_ENUM_TYPE i = to + size, j = from + size;
		if( size > 0 && to != from )
			while( j != from ) data[--j].MoveRow(data[--i]);
	}

	void     Solver::Matrix::Load(std::string file, INMOST_DATA_ENUM_TYPE mbeg, INMOST_DATA_ENUM_TYPE mend)
	{
		char str[16384];
		std::ifstream input(file.c_str());
		if( input.fail() ) throw -1;
		int state = 0, k;
		INMOST_DATA_ENUM_TYPE mat_size, max_lines, row, col, mat_block;
		INMOST_DATA_REAL_TYPE val;
		int size = 1, rank = 0;
#if defined(USE_MPI)
		if( mend == ENUMUNDEF && mbeg == ENUMUNDEF )
		{
			MPI_Comm_rank(GetCommunicator(),&rank);
			MPI_Comm_size(GetCommunicator(),&size);
		}
#endif
		int line = 0;
		while( !input.getline(str,16384).eof() )
		{
			line++;
			k = 0; while( isspace(str[k]) ) k++;
			if( str[k] == '%' || str[k] == '\0' ) continue;
			std::istringstream istr(str+k);
			switch(state)
			{
			case 0: 
				istr >> mat_size >> mat_size >> max_lines; state = 1; 
				mat_block = mat_size/size;
				if( mbeg == ENUMUNDEF ) mbeg = rank*mat_block;
				if( mend == ENUMUNDEF ) 
				{
					if( rank == size-1 ) mend = mat_size;
					else mend = mbeg+mat_block;
				}
				SetInterval(mbeg,mend);
				//~ std::cout << rank << " my interval " << mbeg << ":" << mend << std::endl;
			break;
			case 1: 
				istr >> row >> col >> val;
				row--; col--;
				if( row >= mbeg && row < mend ) data[row][col] = val; 
			break;
			}
		}
		int nonzero = 0;
		for(iterator it = Begin(); it != End(); ++it) nonzero += it->Size();
		//~ std::cout << rank << " total nonzero " << max_lines << " my nonzero " << nonzero << std::endl;
		input.close();
	}


	void     Solver::Vector::Load(std::string file, INMOST_DATA_ENUM_TYPE mbeg, INMOST_DATA_ENUM_TYPE mend)
	{
		char str[16384];
		std::ifstream input(file.c_str());
		if( input.fail() ) throw -1;
		int state = 0, k;
		INMOST_DATA_ENUM_TYPE vec_size, vec_block, ind = 0;
		INMOST_DATA_REAL_TYPE val;
		int size = 1, rank = 0;
#if defined(USE_MPI)
		if( mend == ENUMUNDEF && mbeg == ENUMUNDEF )
		{
			MPI_Comm_rank(GetCommunicator(),&rank);
			MPI_Comm_size(GetCommunicator(),&size);
		}
#endif
		while( !input.getline(str,16384).eof() )
		{
			k = 0; while( isspace(str[k]) ) k++;
			if( str[k] == '%' || str[k] == '\0' ) continue;
			std::istringstream istr(str+k);
			switch(state)
			{
			case 0: 
				istr >> vec_size; state = 1; 
				vec_block = vec_size/size;
				if( mbeg == ENUMUNDEF ) mbeg = rank*vec_block;
				if( mend == ENUMUNDEF ) 
				{
					if( rank == size-1 ) mend = vec_size;
					else mend = mbeg+vec_block;
				}
				SetInterval(mbeg,mend);
			break;
			case 1: 
				istr >> val;
				if( ind >= mbeg && ind < mend ) data[ind] = val;
				ind++;
			break;
			}
		}
		input.close();
	}

	void     Solver::Vector::Save(std::string file)
	{
		INMOST_DATA_ENUM_TYPE vecsize = Size();
		
#if defined(USE_MPI)
		int rank = 0, size = 1;
		{
			MPI_Comm_rank(GetCommunicator(),&rank);
			MPI_Comm_size(GetCommunicator(),&size);
			INMOST_DATA_ENUM_TYPE temp = vecsize;
			MPI_Allreduce(&temp,&vecsize,1,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,GetCommunicator());
		}
#endif
		std::stringstream rhs(std::ios::in | std::ios::out);
		rhs << std::scientific;
		rhs.precision(15);
		for(iterator it = Begin(); it != End(); ++it) rhs << *it << std::endl;
#if defined(USE_MPI) && defined(USE_MPI_FILE) // Use mpi files
		{ 
			int ierr;
			MPI_File fh;
			MPI_Status stat;
			ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(file.c_str()),MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&fh);
			if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			if( rank == 0 )
			{
				std::stringstream header;
				//header << "% vector " << name << std::endl;
				//header << "% is written by INMOST" << std::endl;
				//header << "% by MPI_File_* api" << std::endl;
				header << vecsize << std::endl;
				std::string header_data(header.str());
				ierr = MPI_File_write_shared(fh,&header_data[0],static_cast<int>(header_data.size()),MPI_CHAR,&stat);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			}
			std::string local_data(rhs.str());
			ierr = MPI_File_write_ordered(fh,&local_data[0],static_cast<int>(local_data.size()),MPI_CHAR,&stat);
			if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			ierr = MPI_File_close(&fh);
			if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
		}
#elif defined(USE_MPI) //USE_MPI alternative
		std::string senddata = rhs.str(), recvdata;
		int sendsize = static_cast<int>(senddata.size());
		std::vector<int> recvsize(size), displ(size);
		MPI_Gather(&sendsize,1,MPI_INT,&recvsize[0],1,MPI_INT,0,GetCommunicator());
		if( rank == 0 )
		{
			int totsize = recvsize[0];
			
			displ[0] = 0;
			for(int i = 1; i < size; i++) 
			{
				totsize += recvsize[i];
				displ[i] = displ[i-1]+recvsize[i-1];
			}
			recvdata.resize(totsize);
		}
		else recvdata.resize(1); //protect from dereferencing null
		MPI_Gatherv(&senddata[0],sendsize,MPI_CHAR,&recvdata[0],&recvsize[0],&displ[0],MPI_CHAR,0,GetCommunicator());
		if( rank == 0 )
		{
			std::fstream output(file.c_str(),std::ios::out);
			output << vecsize << std::endl;
			output << recvdata;
		}
#else
		std::fstream output(file.c_str(),std::ios::out);
		//output << "% vector " << name << std::endl;
		//output << "% is written by INMOST" << std::endl;
		//output << "% by sequential write" << std::endl;
		output << vecsize << std::endl;
		output << rhs.rdbuf();
#endif
	}

	void     Solver::Matrix::Save(std::string file)
	{
		INMOST_DATA_ENUM_TYPE matsize = Size(), nonzero = 0, row = GetFirstIndex()+1;
		
		for(iterator it = Begin(); it != End(); ++it) nonzero += it->Size();
#if defined(USE_MPI)
		int rank = 0, size = 1;
		{
			MPI_Comm_rank(GetCommunicator(),&rank);
			MPI_Comm_size(GetCommunicator(),&size);
			INMOST_DATA_ENUM_TYPE temp_two[2] = {matsize,nonzero}, two[2];
			MPI_Allreduce(temp_two,two,2,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,GetCommunicator());
			matsize = two[0];
			nonzero = two[1];
		}
#endif
		std::stringstream mtx(std::ios::in | std::ios::out);
		mtx << std::scientific;
		mtx.precision(15);
		for(iterator it = Begin(); it != End(); ++it)
		{
			for(Row::iterator jt = it->Begin(); jt != it->End(); ++jt)
				mtx << row << " " << jt->first+1 << " " << jt->second << std::endl;
			++row;
		}
#if defined(USE_MPI) && defined(USE_MPI_FILE) // USE_MPI2?
		{
			int ierr;
			MPI_File fh;
			MPI_Status stat;
			ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(file.c_str()),MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&fh);
			if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			if( rank == 0 )
			{
				std::stringstream header;
				header << "%%MatrixMarket matrix coordinate real general" << std::endl;
				header << "% matrix " << name << std::endl;
				header << "% is written by INMOST" << std::endl;
				header << "% by MPI_File_* api" << std::endl;
				header << matsize << " " << matsize << " " << nonzero << std::endl;
				std::string header_data(header.str());
				ierr = MPI_File_write_shared(fh,&header_data[0],static_cast<int>(header_data.size()),MPI_CHAR,&stat);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			}
			std::string local_data(mtx.str());
			ierr = MPI_File_write_ordered(fh,&local_data[0],static_cast<int>(local_data.size()),MPI_CHAR,&stat);
			if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			ierr = MPI_File_close(&fh);
			if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
		}
#elif defined(USE_MPI)//USE_MPI alternative
		std::string senddata = mtx.str(), recvdata;
		int sendsize = static_cast<int>(senddata.size());
		std::vector<int> recvsize(size), displ(size);
		MPI_Gather(&sendsize,1,MPI_INT,&recvsize[0],1,MPI_INT,0,GetCommunicator());
		if( rank == 0 )
		{
			int totsize = recvsize[0];
			
			displ[0] = 0;
			for(int i = 1; i < size; i++) 
			{
				totsize += recvsize[i];
				displ[i] = displ[i-1]+recvsize[i-1];
			}
			recvdata.resize(totsize);
		}
		else recvdata.resize(1); //protect from dereferencing null
		MPI_Gatherv(&senddata[0],sendsize,MPI_CHAR,&recvdata[0],&recvsize[0],&displ[0],MPI_CHAR,0,GetCommunicator());
		if( rank == 0 )
		{
			std::fstream output(file.c_str(),std::ios::out);
			output << "%%MatrixMarket matrix coordinate real general" << std::endl;
			output << "% matrix " << name << std::endl;
			output << "% is written by INMOST" << std::endl;
			output << "% by MPI_Gather* api and sequential write" << std::endl;
			output << matsize << " " << matsize << " " << nonzero << std::endl;
			output << recvdata;
		}
#else
		std::fstream output(file.c_str(),std::ios::out);
		output << "%%MatrixMarket matrix coordinate real general" << std::endl;
		output << "% matrix " << name << std::endl;
		output << "% is written by INMOST" << std::endl;
		output << "% by sequential write " << std::endl;
		output << matsize << " " << matsize << " " << nonzero << std::endl;
		output << mtx.rdbuf();
#endif
	}

	////////////////////////////////////////////////////////////////////////
	//SOLVER SOURCE CODE
	////////////////////////////////////////////////////////////////////////
	void Solver::Initialize(int * argc, char *** argv, const char * database)
	{
		(void)database;
		(void)argc;
		(void)argv;
		if( database != NULL )
		{
			std::fstream file(database,std::ios::in);
			char str[4096];
			while( !file.eof() && file.getline(str,4096) )
			{
				int k = 0, l;
				for(k = 0; k < (int)strlen(str); ++k)
				{
					if( str[k] == ':' ) break;
				}
				if( k == strlen(str) ) continue; //invalid line
				l = k+1;
				while(l < (int)strlen(str) && isspace(str[l]) ) ++l;
				if( l == strlen(str) ) continue; //skip empty entry
				if( !strncmp(str,"PETSc",k) )
					petsc_database_file = std::string(str+l);
				else if( !strncmp(str,"Trilinos_Ifpack",k) )
					trilinos_ifpack_database_file = std::string(str+l);
				else if( !strncmp(str,"Trilinos_Aztec",k) )
					trilinos_aztec_database_file = std::string(str+l);
				else if( !strncmp(str,"Trilinos_ML",k) )
					trilinos_ml_database_file = std::string(str+l);
				else if( !strncmp(str,"Trilinos_Belos",k) )
					trilinos_belos_database_file = std::string(str+l);
				else if( !strncmp(str,"ANI",k) )
					ani_database_file = std::string(str+l);
			}
		}
		//std::cout << "PETSc \"" << petsc_database_file << "\"" << std::endl;
		//std::cout << "Trilinos_Ifpack \"" << trilinos_ifpack_database_file << "\"" << std::endl;
#if defined(USE_SOLVER_PETSC)
		SolverInitializePetsc(argc,argv,petsc_database_file.c_str());
#endif
#if defined(USE_SOLVER_ANI)
		SolverInitializeAni(argc,argv,ani_database_file.c_str());
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
			MPI_Datatype type[3] = { INMOST_MPI_DATA_ENUM_TYPE, INMOST_MPI_DATA_REAL_TYPE, MPI_UB};
			int blocklen[3] = { 1, 1, 1 };
			MPI_Aint disp[3];
			disp[0] = offsetof(Solver::Row::entry,first);
			disp[1] = offsetof(Solver::Row::entry,second);
			disp[2] = sizeof(Solver::Row::entry);
			ierr = MPI_Type_create_struct(3, blocklen, disp, type, &RowEntryType);
			if( ierr != MPI_SUCCESS )
			{
				std::cout << __FILE__ << ":" << __LINE__ << "problem in MPI_Type_create_struct" << std::endl;
			}
			ierr = MPI_Type_commit(&RowEntryType);
			if( ierr != MPI_SUCCESS )
			{
				std::cout << __FILE__ << ":" << __LINE__ << "problem in MPI_Type_commit" << std::endl;
			}
			is_initialized = true;
		}
#endif
	}
	
	void Solver::Finalize()
	{
#if defined(USE_MPI)
		MPI_Type_free(&RowEntryType);
#endif
#if defined(USE_SOLVER_PETSC)
		SolverFinalizePetsc();
#endif
#if defined(USE_SOLVER_ANI)
		SolverFinalizeAni();
#endif
#if defined(USE_MPI)
		{
			int flag = 0;
			MPI_Finalized(&flag);
			if( !flag ) 
				MPI_Finalize();
			is_initialized = false;
			is_finalized = true;
		}
#endif
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
		else if( name == "divergance_tolerance" )
			divergance_tolerance = val;
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
		divergance_tolerance = DEFAULT_DIVERGENCE_TOLERANCE;
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
		if( _pack == INNER_ILU2 || _pack == INNER_MLILUC )
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
			else
			{
#if defined(__SOLVER_DDPQILUC2__)
				prec = new ILUC_preconditioner(info);
#else
				std::cout << "Sorry, multilevel condition estimation crout-ilu preconditioner with reordering for diagonal dominance is not included in this release" << std::endl;
				prec = NULL;
#endif
			}
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
#if defined(USE_SOLVER_TRILINOS)
		if( _pack == Trilinos_Aztec || 
			_pack == Trilinos_Belos ||
			_pack == Trilinos_ML ||
			_pack == Trilinos_Ifpack)
		{
			throw - 1; //You should not really want to copy solver's information
		}
#endif
		if (_pack == INNER_ILU2 || _pack == INNER_MLILUC)
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
		if (_pack == INNER_ILU2 || _pack == INNER_MLILUC)
		{
			if( matrix_data != NULL )
			{
				delete (Solver::Matrix * )matrix_data;
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
#if defined(USE_SOLVER_TRILINOS)
			if( _pack == Trilinos_Aztec || 
				_pack == Trilinos_Belos ||
				_pack == Trilinos_ML ||
				_pack == Trilinos_Ifpack)
			{
				throw - 1; //You should not really want to copy solver's information
			}
#endif
			if (_pack == INNER_ILU2 || _pack == INNER_MLILUC)
			{
				throw - 1; //You should not really want to copy solver's information
			}
		}
		return *this;
	}
	
	void Solver::SetMatrix(Matrix & A, bool OldPreconditioner)
	{
		(void) OldPreconditioner;
		bool ok = false;
#if defined(USE_SOLVER_PETSC)
		if( _pack == PETSc )
		{
			bool modified_pattern = false;
			for(Matrix::iterator it = A.Begin(); it != A.End() && !modified_pattern; ++it)
				modified_pattern |= it->modified_pattern;
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
					
					for(Matrix::iterator it = A.Begin(); it != A.End(); ++it)
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
					for(Matrix::iterator it = A.Begin(); it != A.End(); ++it)
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
				for(Matrix::iterator it = A.Begin(); it != A.End(); ++it)
					if( it->Size() > max ) max = it->Size();
				int * col_positions = new int[max];
				double * col_values = new double[max];
				unsigned k = 0, m;
				for(Matrix::iterator it = A.Begin(); it != A.End(); ++it)
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
			bool modified_pattern = false;
			for(Matrix::iterator it = A.Begin(); it != A.End() && !modified_pattern; ++it)
				modified_pattern |= it->modified_pattern;
			//~ if( A.comm != comm ) throw DifferentCommunicatorInSolver;
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
				for(Matrix::iterator it = A.Begin(); it != A.End(); it++) nnz += it->Size();
				int * ja = (int *)malloc(sizeof(int)*nnz);
				double * values = (double *) malloc(sizeof(double)*nnz);
				ia[0] = shift;
				for(Matrix::iterator it = A.Begin(); it != A.End(); it++)
				{
					for(Row::iterator jt = it->Begin(); jt != it->End(); jt++)
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
				for(Matrix::iterator it = A.Begin(); it != A.End(); it++) nnz += it->Size();
				double * values = (double *)malloc(sizeof(double)*nnz);
				k = 0;
				for(Matrix::iterator it = A.Begin(); it != A.End(); it++)
					for(Row::iterator jt = it->Begin(); jt != it->End(); jt++)
						values[k++] = jt->second;
				MatrixFillValuesAni(matrix_data,values);
			}
			MatrixFinalizeAni(matrix_data);
			SolverSetMatrixAni(solver_data,matrix_data,modified_pattern,OldPreconditioner);
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
			bool modified_pattern = false;
			for(Matrix::iterator it = A.Begin(); it != A.End() && !modified_pattern; ++it)
				modified_pattern |= it->modified_pattern;
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
				for(Matrix::iterator it = A.Begin(); it != A.End(); ++it)
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
		if (_pack == INNER_ILU2 || _pack == INNER_MLILUC)
		{
			
			Solver::Matrix * mat = new Solver::Matrix(A);
			info.PrepareMatrix(*mat, additive_schwartz_overlap);
			IterativeMethod * sol = (IterativeMethod *)solver_data;
			sol->ReplaceMAT(*mat);
			if( matrix_data != NULL ) delete (Solver::Matrix *)matrix_data;
			matrix_data = (void *)mat;

			sol->RealParameter(":tau") = preconditioner_drop_tolerance;
			sol->RealParameter(":tau2") = preconditioner_reuse_tolerance;
			sol->EnumParameter(":scale_iters") = preconditioner_rescale_iterations;

			if( _pack == INNER_MLILUC )
			{
				sol->RealParameter(":ddpq_tau") = preconditioner_ddpq_tolerance;
				sol->EnumParameter(":reorder_nnz") = preconditioner_reorder_nonzero;
				sol->EnumParameter(":estimator") = preconditioner_condition_estimation;
				sol->EnumParameter(":ddpq_tau_adapt") = preconditioner_adapt_ddpq_tolerance;
			}
			else sol->EnumParameter(":fill") = static_cast<INMOST_DATA_ENUM_TYPE>(preconditioner_fill_level);

			if( sizeof(KSOLVER) == sizeof(BCGSL_solver) )
				sol->EnumParameter("levels") = solver_gmres_substeps;

			if (!sol->isInitialized()) 
			{
				sol->Initialize();
			}
			ok = true;

		}
		for(Matrix::iterator it = A.Begin(); it != A.End(); it++) it->modified_pattern = false;
		if(!ok) throw NotImplemented;
	}
	
	bool Solver::Solve(Vector & RHS, Vector & SOL)
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
				for(Solver::Vector::iterator ri = SOL.Begin(); ri != SOL.End(); ++ri) *ri = 0.0;
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
				for(Vector::iterator it = RHS.Begin(); it != RHS.End(); ++it)
				{
					positions[k] = vbeg+k;
					values[k] = *it;
					k++;
				}
				VectorFillPetsc(rhs_data,local_size,positions,values);
				VectorFinalizePetsc(rhs_data);
				
				k = 0;
				for(Vector::iterator it = SOL.Begin(); it != SOL.End(); ++it)
				{
					values[k] = *it;
					k++;
				}
				VectorFillPetsc(solution_data,local_size,positions,values);
				VectorFinalizePetsc(solution_data);
			}
			if( petsc_database_file == "" ) //set parameters if no file is set
			{
				SolverSetTolerancesPetsc(solver_data,relative_tolerance,absolute_tolerance,divergance_tolerance,maximum_iterations);
			}
			bool result = SolverSolvePetsc(solver_data,rhs_data,solution_data);
			if( result )
			{
				VectorLoadPetsc(solution_data,local_size,positions,values);
				unsigned k = 0;
				for(Vector::iterator it = SOL.Begin(); it != SOL.End(); ++it)
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
				return_reason = "The Hessenberg matrix within GMRES is illconditioned."
								"This could be caused by a number"
								"of reasons. For example, the preconditioning"
								"matrix could be nearly singular due to an unstable"
								"factorization (note: pivoting is not implemented"
								"in any of the incomplete factorizations)."
								"Ill-conditioned Hessenberg matrices could also"
								"arise from a singular application matrix. In this"
								"case, GMRES tries to compute a least-squares"
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
		if(_pack == INNER_ILU2 || _pack == INNER_MLILUC)
		{
			IterativeMethod * sol = static_cast<IterativeMethod *>(solver_data);
			
			sol->EnumParameter("maxits") = maximum_iterations;
			sol->RealParameter("rtol") = relative_tolerance;
			sol->RealParameter("atol") = absolute_tolerance;
			sol->RealParameter("divtol") = divergance_tolerance;
			
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
}
#endif

