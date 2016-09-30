#include <inmost_sparse.h>
#include <inmost_solver.h>
#include "SolverFactory.h"

#include "solver_inner/solver_ilu2/SolverILU2.h"

#if defined(USE_SOLVER_PETSC)
    #include "solver_petsc/SolverPETSc.h"
#endif

namespace INMOST {

    int *Solver::argc = NULL;
    char ***Solver::argv = NULL;
    const char *Solver::database = NULL;
    bool Solver::is_initialized = false;
    bool Solver::is_finalized = false;

    Solver::Solver(std::string solverName, std::string prefix, INMOST_MPI_Comm _comm) {
        this->solver = SolverFactory::getSolver(solverName);
        this->prefix = prefix;

        solver->SetCommunicator(_comm);
        std::string solverDatabasePath = Solver::parseDatabase(solverName);
        solver->Initialize(argc, argv, solverDatabasePath.c_str(), prefix);
    }

    Solver::Solver(const Solver &other) {
        this->solver = SolverFactory::copySolver(other.solver);
        this->prefix = other.prefix;

        solver->SetCommunicator(other.solver->GetCommunicator());
        std::string solverDatabasePath = Solver::parseDatabase(solver->SolverName());
        solver->Initialize(argc, argv, solverDatabasePath.c_str(), prefix);
    }

    Solver& Solver::operator=(const Solver& other) {
        if( this != &other ) {
            this->solver->SetCommunicator(other.solver->GetCommunicator());
            this->prefix = other.prefix;
            this->solver->Assign(other.solver);
        }
        return *this;
    }

    std::string Solver::SolverName() const {
        return solver->SolverName();
    }

    std::string Solver::SolverPrefix() const {
        return prefix;
    }

    void Solver::Initialize(int *argc, char ***argv, const char *database) {
        Solver::argc = argc;
        Solver::argv = argv;
        Solver::database = database;
        Solver::is_initialized = true;
        Solver::is_finalized = false;
#if defined(USE_MPI)
        {
            int flag = 0;
            int ierr = 0;
            MPI_Initialized(&flag);
            if (!flag) {
                ierr = MPI_Init(argc,argv);
                if( ierr != MPI_SUCCESS ) {
                    std::cout << __FILE__ << ":" << __LINE__ << "problem in MPI_Init" << std::endl;
                }
            }
        }
#endif
        //Register all available solvers
        SolverFactory::registerSolver<SolverILU2>("inner_ilu2");
#if defined(USE_SOLVER_PETSC)
        SolverFactory::registerSolver<SolverPETSc>("petsc");
#endif
        Sparse::CreateRowEntryType();
    }

    void Solver::Finalize() {
        Sparse::DestroyRowEntryType();
#if defined(USE_MPI)
    {
        int flag = 0;
        MPI_Finalized(&flag);
        if( !flag ) {
            MPI_Finalize();
        }
    }
#endif
        Solver::is_finalized = true;
        Solver::is_initialized = false;
    }

    bool Solver::isInitialized() {
        return is_initialized;
    }

    bool Solver::isFinalized() {
        return is_finalized;
    }

    void Solver::SetMatrix(Sparse::Matrix & A, bool ModifiedPattern, bool OldPreconditioner) {
        solver->SetMatrix(A, ModifiedPattern, OldPreconditioner);
    }

    bool Solver::Solve(INMOST::Sparse::Vector & RHS, INMOST::Sparse::Vector & SOL) {
        if( !solver->isMatrixSet()) throw MatrixNotSetInSolver;
        if( RHS.GetCommunicator() != solver->GetCommunicator() || SOL.GetCommunicator() != solver->GetCommunicator()) throw DifferentCommunicatorInSolver;
        INMOST_DATA_ENUM_TYPE vbeg,vend;
        RHS.GetInterval(vbeg,vend);
        if (RHS.Size() != SOL.Size()) {
            if (SOL.Size() == 0) {
                SOL.SetInterval(vbeg,vend);
                for(Sparse::Vector::iterator ri = SOL.Begin(); ri != SOL.End(); ++ri) *ri = 0.0;
            } else throw InconsistentSizesInSolver;
        }
        return solver->Solve(RHS, SOL);
    }

        bool Solver::Clear() {
        return solver->Clear();
    }

    INMOST_DATA_REAL_TYPE Solver::GetPropertyReal(std::string property) const {
        return solver->GetPropertyReal(property);
    }

    INMOST_DATA_ENUM_TYPE Solver::GetPropertyEnum(std::string property) const {
        return solver->GetPropertyEnum(property);
    }

    void Solver::SetPropertyReal(std::string property, INMOST_DATA_REAL_TYPE value) {
        solver->SetPropertyReal(property, value);
    }

    void Solver::SetPropertyEnum(std::string property, INMOST_DATA_ENUM_TYPE value) {
        solver->SetPropertyEnum(property, value);
    }

    const INMOST_DATA_ENUM_TYPE Solver::Iterations() const {
        return solver->Iterations();
    }

    const INMOST_DATA_REAL_TYPE Solver::Residual() const {
        return solver->Residual();
    }

    const std::string Solver::ReturnReason() const {
        return solver->ReturnReason();
    }

    Solver::~Solver() {
        solver->Finalize();
        delete solver;
    }

    std::string Solver::parseDatabase(std::string solverName) {
        const char *name = solverName.c_str();
        if( database != NULL ) {
            FILE * f = fopen(database, "r");
            if (f != NULL) {
                char str[4096];
                while( !feof(f) && fgets(str,4096,f)) {
                    int k = 0, l;
                    for(k = 0; k < (int)strlen(str); ++k) {
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
                    if( !strncmp(str, name, k) ) {
                        return std::string(str+l);
                    }
                }
                fclose(f);
            }
        }
        return std::string("");
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
    

}
