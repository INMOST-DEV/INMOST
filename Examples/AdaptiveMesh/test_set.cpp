#include "amesh.h"
#include <fstream>
#include <iomanip>
#include <sstream>


using namespace std;
using namespace INMOST;

void add_elem_to_set(AdaptiveMesh& m, ElementSet& set, Element& elem)
{
    set.PutElement(elem);
    m.parent_set[elem] = set.GetHandle();
}

int main(int argc, char ** argv)
{
    Mesh::Initialize(&argc,&argv);

	Mesh m;
    AdaptiveMesh am(m);
    m.Load("grid_for_set.pvtk");

    int size = m.GetProcessorsNumber();
    int rank = m.GetProcessorRank();
    m.SetCommunicator(MPI_COMM_WORLD);

    m.ResolveShared();
    am.UpdateStatus();
    m.Save("test_sets_begin.pvtk");

    m.ResolveModification();
    am.UpdateStatus();
    m.Save("test_sets_resolve.pvtk");
    
    ElementSet base;
    ElementSet C1;
    ElementSet C2;

    // Create sets.
    base = m.CreateSetUnique("BASE").first;
    C1   = m.CreateSetUnique("C1").first;
    C2   = m.CreateSetUnique("C2").first;

    // Set relations for sets
    base.AddChild(C1);
    base.AddChild(C2);

    // Syncronize status of sets
    m.ResolveSets();

    // Fill sets
    for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) 
    {
        int gl = it->GlobalID();
        if (rank == 0 && gl == 0 ) add_elem_to_set(am,base,it->self());
        if (rank == 1 && gl == 16) add_elem_to_set(am,base,it->self());

        if (gl >= 1 && gl <= 8             ) add_elem_to_set(am,C1,it->self());
        if (gl >= 9 && gl <= 17 && gl != 16) add_elem_to_set(am,C2,it->self());
    }

    // PrintSets
    am.PrintSet();

    // Print Mesh cells to beforeSync_rank file
    {
        stringstream ss;
        ss << "beforeSync_" << rank;
        ofstream ofs(ss.str().c_str());
        am.PrintMesh(ofs,1,0,0,0);
    }
    
    // Synchronize sets. 
    am.SynchronizeSet(C1); // C1 = C1 (on 0 rank) union C1 (on 1 rank)
    am.SynchronizeSet(C2); 

    // Print Mesh cells to afterSync_rank file
    {
        stringstream ss;
        ss << "afterSync_" << rank;
        ofstream ofs(ss.str().c_str());
        am.PrintMesh(ofs,1,0,0,0);
    }
    
    // PrintSets
    am.PrintSet();

    am.UpdateStatus();
    m.Save("after_set_sync.pvtk");

	Mesh::Finalize();
    return 0;
}

