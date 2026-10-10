#include "inmost.h"
#include <cstdio>
#include <iostream>
#include <set>

using namespace INMOST;

// Inspect real message matching through MPI's profiling interface, without
// exposing private mesh state. A small bound also exercises exhaustion/reuse.
static bool capture = false;
static int tag_bound = 32767, failures = 0, rank = 0;
static int attribute_queries = 0;
static std::set<int> message_tags;

extern "C" int MPI_Comm_get_attr(MPI_Comm comm, int key, void *value, int *flag)
{
    if (key == MPI_TAG_UB)
    {
        if (capture) ++attribute_queries;
        *static_cast<int **>(value) = &tag_bound;
        *flag = 1;
        return MPI_SUCCESS;
    }
    return PMPI_Comm_get_attr(comm, key, value, flag);
}

extern "C" int MPI_Attr_get(MPI_Comm comm, int key, void *value, int *flag)
{
    return MPI_Comm_get_attr(comm, key, value, flag);
}

extern "C" int MPI_Isend(const void *buf, int count, MPI_Datatype type,
                          int dest, int tag, MPI_Comm comm, MPI_Request *req)
{
    if (capture) message_tags.insert(tag);
    return PMPI_Isend(buf, count, type, dest, tag, comm, req);
}

extern "C" int MPI_Irecv(void *buf, int count, MPI_Datatype type,
                          int source, int tag, MPI_Comm comm, MPI_Request *req)
{
    if (capture) message_tags.insert(tag);
    return PMPI_Irecv(buf, count, type, source, tag, comm, req);
}

static void check(bool condition, const char *message)
{
    if (!condition)
    {
        std::cerr << "Rank " << rank << ": " << message << std::endl;
        ++failures;
    }
}

static void start_capture()
{
    message_tags.clear();
    attribute_queries = 0;
    capture = true;
}

static void fill(Mesh &mesh, Tag tag, int stamp)
{
    for (Mesh::iteratorCell cell = mesh.BeginCell(); cell != mesh.EndCell(); ++cell)
    {
        bool ghost = cell->GetStatus() == Element::Ghost;
        if (tag.GetSize() == ENUMUNDEF)
        {
            Storage::integer_array values = cell->IntegerArray(tag);
            values.resize(ghost ? 0 : 1 + cell->GlobalID() % 3);
            for (size_t k = 0; k < values.size(); ++k)
                values[k] = stamp + cell->GlobalID() + static_cast<int>(k);
        }
        else cell->Integer(tag) = ghost ? -1 : stamp + cell->GlobalID();
    }
}

static void verify(Mesh &mesh, Tag tag, int stamp)
{
    for (Mesh::iteratorCell cell = mesh.BeginCell(); cell != mesh.EndCell(); ++cell)
    {
        if (tag.GetSize() == ENUMUNDEF)
        {
            Storage::integer_array values = cell->IntegerArray(tag);
            check(values.size() == 1 + cell->GlobalID() % 3, "Variable buffer size mismatch");
            for (size_t k = 0; k < values.size(); ++k)
                check(values[k] == stamp + cell->GlobalID() + static_cast<int>(k),
                      "Variable exchange received the wrong message");
        }
        else
        {
            if (cell->Integer(tag) != stamp + cell->GlobalID() && failures < 8)
                std::cerr << "Rank " << rank << " mesh " << mesh.GetMeshName()
                          << " stamp " << stamp << " global ID " << cell->GlobalID()
                          << " actual " << cell->Integer(tag) << std::endl;
            check(cell->Integer(tag) == stamp + cell->GlobalID(),
                  "Fixed exchange received the wrong message");
        }
    }
}

static void setup(Mesh &mesh, const char *file, int size)
{
    mesh.SetCommunicator(MPI_COMM_WORLD);
    if (!rank) mesh.Load(file);
    mesh.AssignGlobalID(CELL);
    Tag destination = mesh.RedistributeTag();
    for (Mesh::iteratorCell cell = mesh.BeginCell(); cell != mesh.EndCell(); ++cell)
        cell->Integer(destination) = cell->GlobalID() % size;
    mesh.Redistribute();
    mesh.ExchangeGhost(1, FACE);
    mesh.CreateTag("tag_test_a", DATA_INTEGER, CELL, NONE, 1);
}

static int data_tag(Mesh &mesh)
{
    mesh.RecomputeParallelStorage(CELL);
    Tag tag = mesh.GetTag("tag_test_a");
    fill(mesh, tag, 100);
    start_capture();
    mesh.ExchangeData(tag, CELL);
    capture = false;
    verify(mesh, tag, 100);
    check(message_tags.size() == 1, "Fixed exchange did not use one data tag");
    check(attribute_queries == 0, "Exchange queried MPI_TAG_UB");
    return message_tags.empty() ? -1 : *message_tags.begin();
}

static void solver_exchange(Mesh &mesh, int size)
{
#if defined(USE_SOLVER)
    Tag tag = mesh.GetTag("tag_test_a");
    fill(mesh, tag, 600);
    Mesh::exchange_data pending;
    start_capture();
    mesh.ExchangeDataBegin(tag, CELL, 0, pending);
    Sparse::Matrix matrix("tag_test_overlap", 2*rank, 2*rank+2);
    Sparse::Vector vector("tag_test_overlap", 2*rank, 2*rank+2);
    for (int i = 2*rank; i < 2*rank+2; ++i)
    {
        matrix[i][i] = 4;
        if (i) matrix[i][i-1] = -1;
        if (i+1 < 2*size) matrix[i][i+1] = -1;
        vector[i] = 1;
    }
    Solver::OrderInfo info;
    info.PrepareMatrix(matrix, 1);
    info.PrepareVector(vector);
    info.Update(vector);
    for (Sparse::Vector::iterator i = vector.Begin(); i != vector.End(); ++i)
        check(*i == 1, "Solver Update received the wrong values");
    int contributions = static_cast<int>(vector.Size()), expected = 0;
    MPI_Allreduce(&contributions, &expected, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    info.Accumulate(vector);
    info.RestoreVector(vector);
    int sum = 0, total = 0;
    for (Sparse::Vector::iterator i = vector.Begin(); i != vector.End(); ++i)
        sum += static_cast<int>(*i);
    MPI_Allreduce(&sum, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    check(total == expected, "Solver Accumulate lost contributions");
    mesh.ExchangeDataEnd(tag, CELL, 0, pending);
    capture = false;
    check(message_tags.count(MPIExchangeTag::SolverUpdate) &&
          message_tags.count(MPIExchangeTag::SolverAccumulate),
          "Solver exchanges did not use their fixed tags");
    verify(mesh, tag, 600);
#else
    (void) mesh;
    (void) size;
#endif
}

static void roundtrip(Mesh &mesh, const std::string &prefix)
{
    const std::string file = prefix + ".pmf";
    // Loading may renumber protected global IDs; preserve the original ID
    // as ordinary user data to check values independently of that numbering.
    Tag original_id = mesh.CreateTag("tag_test_original_id", DATA_INTEGER, CELL, NONE, 1);
    for (Mesh::iteratorCell cell = mesh.BeginCell(); cell != mesh.EndCell(); ++cell)
        cell->Integer(original_id) = cell->GlobalID();
    mesh.SetParallelFileStrategy(0); // Exercise point-to-point PMF I/O.
    start_capture();
    mesh.Save(file);
    capture = false;
    check(message_tags.size() == 1 && message_tags.count(MPIExchangeTag::Pmf),
          "PMF save did not use its fixed tag");
    {
        Mesh loaded;
        loaded.SetCommunicator(MPI_COMM_WORLD);
        loaded.SetParallelFileStrategy(0);
        start_capture();
        loaded.Load(file);
        capture = false;
        check(message_tags.count(MPIExchangeTag::Pmf), "PMF load did not use its fixed tag");
        check(loaded.NumberOfCells() == mesh.NumberOfCells(), "PMF roundtrip lost cells");
        Tag restored_id = loaded.GetTag("tag_test_original_id");
        Tag restored_data = loaded.GetTag("tag_test_a");
        for (Mesh::iteratorCell cell = loaded.BeginCell(); cell != loaded.EndCell(); ++cell)
            check(cell->Integer(restored_data) == 100 + cell->Integer(restored_id),
                  "PMF roundtrip corrupted cell values");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (!rank) std::remove(file.c_str());

#if defined(USE_SOLVER)
    Sparse::Matrix matrix("tag_test_matrix", 2*rank, 2*rank+2);
    Sparse::Vector vector("tag_test_vector", 2*rank, 2*rank+2);
    for (int i = 2*rank; i < 2*rank+2; ++i)
    {
        matrix[i][i] = 4;
        vector[i] = i + 1;
    }
    const std::string matfile = prefix + ".mtx", vecfile = prefix + ".rhs";
    matrix.Save(matfile);
    vector.Save(vecfile);
    Sparse::Matrix restored_matrix;
    Sparse::Vector restored_vector;
    restored_matrix.Load(matfile, 2*rank, 2*rank+2);
    restored_vector.Load(vecfile, 2*rank, 2*rank+2);
    for (int i = 2*rank; i < 2*rank+2; ++i)
    {
        check(restored_matrix[i][i] == 4, "Sparse matrix roundtrip corrupted values");
        check(restored_vector[i] == i+1, "Sparse vector roundtrip corrupted values");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (!rank)
    {
        std::remove(matfile.c_str());
        std::remove(vecfile.c_str());
    }
#endif
}

int main(int argc, char **argv)
{
    Mesh::Initialize(&argc, &argv);
#if defined(USE_SOLVER)
    Solver::Initialize(&argc, &argv, NULL);
#endif
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (argc < 3 || size < 2) MPI_Abort(MPI_COMM_WORLD, 1);
    {
        Mesh first;
        setup(first, argv[1], size);
        Tag a = first.CreateTag("tag_test_a", DATA_INTEGER, CELL, NONE, 1);
        Tag b = first.CreateTag("tag_test_b", DATA_INTEGER, CELL, NONE, 1);
        Tag variable = first.CreateTag("tag_test_variable", DATA_INTEGER, CELL, NONE);
        const int first_tag = data_tag(first);

        fill(first, a, 100); fill(first, b, 200); fill(first, variable, 300);
        Mesh::exchange_data sa, sb, sv;
        start_capture();
        first.ExchangeDataBegin(a, CELL, 0, sa);
        first.ExchangeDataBegin(variable, CELL, 0, sv);
        first.ExchangeDataBegin(b, CELL, 0, sb);
        first.ExchangeDataEnd(b, CELL, 0, sb);
        first.ExchangeDataEnd(variable, CELL, 0, sv);
        first.ExchangeDataEnd(a, CELL, 0, sa);
        capture = false;
        check(attribute_queries == 0, "Asynchronous exchange queried MPI_TAG_UB");
#if !defined(USE_MPI_P2P) || !defined(PREFFER_MPI_P2P)
        check(message_tags.size() == 2 && message_tags.count(first_tag-1) &&
              message_tags.count(first_tag), "Sizes and data did not use separate tags");
#endif
        verify(first, a, 100); verify(first, b, 200); verify(first, variable, 300);

        Mesh *second = new Mesh;
        setup(*second, argv[1], size);
        const int second_tag = data_tag(*second);
        check(second_tag != first_tag, "Different meshes share a tag");
        fill(first, a, 400);
        Tag other_a = second->GetTag("tag_test_a");
        fill(*second, other_a, 500);
        // Mesh order intentionally differs between ranks; their tags isolate them.
        if (rank % 2)
        {
            first.ExchangeDataBegin(a, CELL, 0, sa);
            second->ExchangeDataBegin(other_a, CELL, 0, sb);
        }
        else
        {
            second->ExchangeDataBegin(other_a, CELL, 0, sb);
            first.ExchangeDataBegin(a, CELL, 0, sa);
        }
        second->ExchangeDataEnd(other_a, CELL, 0, sb);
        first.ExchangeDataEnd(a, CELL, 0, sa);
        verify(first, a, 400); verify(*second, other_a, 500);
        delete second;

        for (int i = 0; i < 8; ++i)
        {
            Mesh replacement;
            setup(replacement, argv[1], size);
            check(data_tag(replacement) == second_tag, "Highest freed tag pair was not reused");
            replacement.SetCommunicator(MPI_COMM_WORLD);
            check(data_tag(replacement) == second_tag, "Identical communicator changed tags");
        }
        {
            Mesh empty, assigned, serial;
            empty.SetCommunicator(MPI_COMM_WORLD);
            Mesh copied(empty);
            assigned = empty;
            assigned = serial;
            Mesh replacement;
            setup(replacement, argv[1], size);
            check(data_tag(replacement) == second_tag+4,
                  "Copy/assignment leaked a tag pair or reused a live pair");
        }
        {
            // All ranks must respect the smallest bound, including its last tag.
            tag_bound = second_tag + 2*rank;
            Mesh last;
            setup(last, argv[1], size);
            check(data_tag(last) == second_tag, "Last available tag pair was rejected");
            Mesh exhausted;
            bool rejected = false;
            try { exhausted.SetCommunicator(MPI_COMM_WORLD); }
            catch (ErrorType error) { rejected = error == NoSpaceForMpiTag; }
            check(rejected, "Exhaustion wrapped tags instead of reporting an error");
            tag_bound = 32767;
        }
        {
            Mesh replacement;
            setup(replacement, argv[1], size);
            check(data_tag(replacement) == second_tag, "Exhaustion leaked a tag pair");
        }
        {
            // Local registries can differ: the world allocation must avoid
            // the pair still occupied by a mesh on a single-rank communicator.
            Mesh *local = NULL;
            if (!rank)
            {
                local = new Mesh;
                local->SetCommunicator(MPI_COMM_SELF);
            }
            {
                Mesh shared;
                setup(shared, argv[1], size);
                check(data_tag(shared) == second_tag+2,
                      "Collective allocation ignored another rank's live IDs");
            }
            delete local;
        }
        solver_exchange(first, size);
        data_tag(first);
        roundtrip(first, argv[2]);
    }
    int total = 0;
    MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (!rank && !total) std::cout << "MPI tag isolation, ordering, reuse and I/O OK" << std::endl;
#if defined(USE_SOLVER)
    Solver::Finalize();
#endif
    Mesh::Finalize();
    return total ? 1 : 0;
}
