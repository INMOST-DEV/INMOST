#include "inmost.h"
#include <iostream>
#include <map>

// Use MPI's profiling interface to detect resource leaks even on implementations
// that silently tolerate them at MPI_Finalize.
static std::map<MPI_Win, void *> live_windows;
static int live_allocations = 0;
static int created_windows = 0;
static int failures = 0;

static void check(bool condition, const char *message)
{
    if (!condition)
    {
        std::cerr << message << std::endl;
        ++failures;
    }
}

extern "C" int MPI_Win_create(void *base, MPI_Aint size, int disp_unit,
                               MPI_Info info, MPI_Comm comm, MPI_Win *win)
{
    int err = PMPI_Win_create(base, size, disp_unit, info, comm, win);
    if (err == MPI_SUCCESS)
    {
        live_windows[*win] = base;
        ++created_windows;
    }
    return err;
}

extern "C" int MPI_Win_free(MPI_Win *win)
{
    MPI_Win old = *win;
    int err = PMPI_Win_free(win);
    if (err == MPI_SUCCESS) live_windows.erase(old);
    return err;
}

extern "C" int MPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr)
{
    int err = PMPI_Alloc_mem(size, info, baseptr);
    if (err == MPI_SUCCESS) ++live_allocations;
    return err;
}

extern "C" int MPI_Free_mem(void *base)
{
    for (std::map<MPI_Win, void *>::const_iterator it = live_windows.begin();
         it != live_windows.end(); ++it)
        check(it->second != base, "Window memory freed before MPI_Win_free");
    int err = PMPI_Free_mem(base);
    if (err == MPI_SUCCESS) --live_allocations;
    return err;
}

int main(int argc, char **argv)
{
    INMOST::Mesh::Initialize(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm duplicate, subgroup;
    MPI_Comm_dup(MPI_COMM_WORLD, &duplicate);
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &subgroup);
    {
        INMOST::Mesh mesh;
        mesh.SetCommunicator(MPI_COMM_WORLD);
        mesh.SetParallelFileStrategy(0);
        int initial_windows = created_windows;
        for (int i = 0; i < 8; ++i) mesh.SetCommunicator(MPI_COMM_WORLD);
#if defined(USE_MPI_P2P)
        check(created_windows == initial_windows, "Identical communicator recreated window");
        check(mesh.GetParallelFileStrategy() == 0, "Identical communicator reset settings");
#endif
        mesh.SetCommunicator(duplicate); // Congruent, but not identical.
        mesh.SetCommunicator(subgroup); // A different group and window size.
        check(mesh.GetProcessorsNumber() == (size + (rank % 2 == 0 ? 1 : 0)) / 2,
              "Invalid subgroup size");
        mesh.SetCommunicator(MPI_COMM_WORLD);
        INMOST::Mesh copy(mesh);
        copy.SetCommunicator(MPI_COMM_WORLD);
        INMOST::Mesh assigned;
        assigned = mesh; // Serial -> parallel.
        assigned = copy; // Parallel -> parallel on the same communicator.
        assigned.SetCommunicator(duplicate);
        assigned = mesh; // Parallel -> parallel on a different communicator.
        INMOST::Mesh serial;
        assigned = serial; // Parallel -> serial releases its window.
        assigned.SetCommunicator(MPI_COMM_WORLD);
    }
    check(live_windows.empty(), "Mesh leaked MPI windows");
    check(live_allocations == 0, "Mesh leaked MPI allocations");
    MPI_Comm_free(&subgroup);
    MPI_Comm_free(&duplicate);
    int global_failures = 0;
    MPI_Allreduce(&failures, &global_failures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    INMOST::Mesh::Finalize();
    if (rank == 0 && global_failures == 0)
        std::cout << "MPI window lifecycle and finalization OK" << std::endl;
    return global_failures ? 1 : 0;
}
