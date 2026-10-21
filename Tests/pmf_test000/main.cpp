#include "inmost.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
using namespace INMOST;
static int rank;
static void check(bool ok, const char *message)
{
    if (!ok) { std::cerr << "rank " << rank << ": " << message << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); }
}
void generate(INMOST::Mesh &mesh, int n) {
    using namespace INMOST;
    mesh.SetCommunicator(INMOST_MPI_COMM_WORLD);
    ElementArray<Node> nodes(&mesh);
    for (int z = 0; z <= n; ++z)
        for (int y = 0; y <= n; ++y)
            for (int x = 0; x <= n; ++x) {
                Storage::real xyz[] = {Storage::real(rank) * n + x, Storage::real(y), Storage::real(z)};
                nodes.push_back(mesh.CreateNode(xyz));
            }
    Storage::integer faces[] = {0, 2, 3, 1, 4, 5, 7, 6, 0, 1, 5, 4, 2, 6, 7, 3, 0, 4, 6, 2, 1, 3, 7, 5};
    Storage::integer sizes[] = {4, 4, 4, 4, 4, 4};
    for (int z = 0; z < n; ++z)
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x) {
                ElementArray<Node> cell(&mesh);
                for (int dz = 0; dz < 2; ++dz)
                    for (int dy = 0; dy < 2; ++dy)
                        for (int dx = 0; dx < 2; ++dx)
                            cell.push_back(nodes[((z + dz) * (n + 1) + y + dy) * (n + 1) + x + dx]);
                mesh.CreateCell(cell, faces, sizes, 6, cell);
            }
    mesh.ResolveShared();
    mesh.ExchangeGhost(1, FACE);
    Tag tag = mesh.CreateTag("temperature", DATA_REAL, NODE | CELL, NONE, 1);
    for (auto it = mesh.Begin(NODE | CELL); it != mesh.End(); ++it)
        it->Real(tag) = Storage::real(rank) + 0.125;
}

static std::vector<double> values(Mesh &mesh)
{
    std::vector<double> result;
    Tag tag = mesh.GetTag("temperature");
    for (Mesh::iteratorElement it = mesh.BeginElement(NODE | CELL); it != mesh.EndElement(); ++it)
        result.push_back(it->Real(tag));
    std::sort(result.begin(), result.end());
    return result;
}
static std::vector<char> bytes(const std::string &path)
{
    std::ifstream in(path.c_str(), std::ios::binary);
    check(in.good(), "cannot open saved PMF");
    return std::vector<char>((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}
int main(int argc, char **argv)
{
    Mesh::Initialize(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    check(argc == 2, "expected output prefix");
    const std::string file = std::string(argv[1]) + ".pmf";
    const std::string reference = std::string(argv[1]) + "-reference.pmf";
    {
        // Reuse the same path for large and small meshes: stale tails must go.
        for (int iteration = 0; iteration < 6; ++iteration)
        {
            Mesh source;
            generate(source, iteration % 2 ? 1 : 3);
            source.SetParallelFileStrategy(0);
            source.Save(reference);
            MPI_Barrier(MPI_COMM_WORLD);
            source.SetParallelFileStrategy(1);
            source.Save(file);
            MPI_Barrier(MPI_COMM_WORLD);
            if (!rank) check(bytes(file) == bytes(reference), "MPI-IO output differs from strategy 0");
            for (int strategy = 0; strategy <= 1; ++strategy)
            {
                Mesh restored;
                restored.SetCommunicator(MPI_COMM_WORLD);
                restored.SetParallelFileStrategy(strategy);
                restored.Load(file);
                check(restored.NumberOfNodes() == source.NumberOfNodes(), "lost nodes");
                check(restored.NumberOfCells() == source.NumberOfCells(), "lost cells");
                check(values(restored) == values(source), "changed tag values");
            }
        }
        // Only the last rank receives a short payload; all ranks must throw
        // before mesh reconstruction or collective mesh cleanup can diverge.
        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, const_cast<char *>(file.c_str()), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        MPI_Offset size;
        MPI_File_get_size(fh, &size);
        MPI_File_set_size(fh, size - 1);
        MPI_File_close(&fh);
        MPI_Barrier(MPI_COMM_WORLD);
        Mesh restored;
        restored.SetCommunicator(MPI_COMM_WORLD);
        restored.SetParallelFileStrategy(1);
        bool rejected = false;
        try { restored.Load(file); }
        catch (ErrorType error) { rejected = error == BadFile; }
        check(rejected, "truncated PMF was not rejected on every rank");
        MPI_Barrier(MPI_COMM_WORLD);
        if (!rank) { std::remove(file.c_str()); std::remove(reference.c_str()); }
    }
    Mesh::Finalize();
    return 0;
}
