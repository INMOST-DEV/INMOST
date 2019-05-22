#include <utility>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <map>
#include "inmost.h"
#include "utils.h"
#include <fstream>
#include <string>
#include "inmost_ani_fem.h"

//#define USE_MPI
//#define USE_PARTITIONER
using namespace INMOST;
#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif



int main(int argc, char *argv[]) {
    int i, j, k;
    unsigned long int memory_test;
    int processRank = 0, processorsCount = 1;

#if defined(USE_MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(INMOST_MPI_COMM_WORLD, &processRank);  // Get the rank of the current process
    MPI_Comm_size(INMOST_MPI_COMM_WORLD, &processorsCount); // Get the total number of processors used
#endif
////////////////////////////////////////////////////////get parameters from arguments///////////////////////////////////
    std::string parametersFileName = "";
    std::string solverName = "inner_ilu2";
    std::string MeshName = "";
    int refineNumber = 0;

    bool parametersFound = false;
    bool typeFound = false;
    bool meshFound = false;
    bool refine = false;

    //Parse argv parameters
    if (argc == 1) goto helpMessage;
    for (i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            helpMessage:
            if (processRank == 0) {
                std::cout << "Help message: " << std::endl;
                std::cout << "Command line options: " << std::endl;
                std::cout << "-d, --database <Solver parameters file name>" << std::endl;
                std::cout << "-t, --type <Solver type name>" << std::endl;
                std::cout << "-am, --animesh <Mesh name>" << std::endl;
                std::cout << "--refine <number of refines>" << std::endl;
                std::cout << "  Available solvers:" << std::endl;
                Solver::Initialize(NULL, NULL, NULL);
                std::vector<std::string> availableSolvers = Solver::getAvailableSolvers();
                for (solvers_names_iterator_t it = availableSolvers.begin(); it != availableSolvers.end(); it++) {
                    std::cout << "      " << *it << std::endl;
                }
                Solver::Finalize();
            }
            return 0;
        }
        //Parameters file name found with -d or --database options
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--database") == 0) {
            if (processRank == 0) {
                std::cout << "Solver parameters file found: " << argv[i + 1] << std::endl;
            }
            parametersFound = true;
            parametersFileName = std::string(argv[i + 1]);
            i++;
            continue;
        }
        //Solver type found with -t ot --type options
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--type") == 0) {
            if (processRank == 0) {
                std::cout << "Solver type index found: " << argv[i + 1] << std::endl;
            }
            typeFound = true;
            solverName = std::string(argv[i + 1]);
            i++;
            continue;
        }

        if (strcmp(argv[i], "-am") == 0 || strcmp(argv[i], "--animesh") == 0) {
            if (processRank == 0) {
                std::cout << "mesh name found: " << argv[i + 1] << std::endl;
            }
            meshFound = true;
            MeshName = std::string(argv[i + 1]);
            i++;
            continue;
        }
        if ( strcmp(argv[i], "--refine") == 0) {
            if (processRank == 0) {
                std::cout << "mesh will be refined " << argv[i + 1] <<" times" << std::endl;
            }
            refine = true;
            refineNumber = atoi(argv[i + 1]);
            i++;
            continue;
        }
    }
    if (!typeFound)
        if (processRank == 0)
            std::cout
                    << "Solver type not found in command line, you can specify solver type with -t or --type option, using INNER_ILU2 solver by default."
                    << std::endl;

    if (!parametersFound)
        if (processRank == 0)
            std::cout << "Parameters not found, you can specify parameter file name with -d or --database option, "<< std::endl;

    if (!meshFound) {
    if (processRank == 0)
        std::cout
                << "Mesh name not found , you can specify mesh name  with -am or --animesh option."
                << std::endl;
        exit(-1);
    }
/////////////////////////////////////////////////////////////////////////////end get parameters/////////////////////////////////

    Mesh::Initialize(&argc, &argv);
    Partitioner::Initialize(&argc, &argv);
    Solver::Initialize(&argc, &argv, parametersFound ? parametersFileName.c_str() : NULL);

    if (!Solver::isSolverAvailable(solverName)) {
        if (processRank == 0) {
            std::cout << "Solver " << solverName << " is not available" << std::endl;
        }
        Mesh::Finalize();
        Partitioner::Finalize();
        Solver::Finalize();
        exit(-1);
    }
    Solver s = Solver(solverName, "test");

    if (processRank == 0) {
        std::cout << "Solving with " << solverName << std::endl;
    }
    double t0, t_init, t_prepare, t_assmble, t_solve;
    BARRIER
    t0 = Timer();

    Mesh *mesh_init;
    mesh_init = new Mesh();
#if defined(USE_MPI)
    mesh_init->SetCommunicator(INMOST_MPI_COMM_WORLD);
#endif
    if(MeshName.find(".ani") != std::string::npos) {
        if (mesh_init->GetProcessorRank() == 0) {
            ReadInmostFromAniFile(mesh_init, MeshName);
            std::cout << " loading success of mesh " << MeshName << std::endl;
        }
    }
    else{
        if(MeshName.find(".pmf") != std::string::npos){
            mesh_init->Load(MeshName);
        }
        else{
            if (processRank == 0)
                std::cout<<"Wrong mesh format"<<std::endl;
            Mesh::Finalize();
            Partitioner::Finalize();
            Solver::Finalize();
            exit(-1);
        }
    }

#ifdef USE_PARTITIONER
    Partitioner *p = new Partitioner(mesh_init);
//     p->SetMethod(Partitioner::Inner_RCM,Partitioner::Partition);
#ifdef USE_PARTITIONER_PARMETIS
    p->SetMethod(Partitioner::Parmetis, Partitioner::Partition);
#elif USE_PARTITIONER_ZOLTAN
    p->SetMethod(Partitioner::Zoltan, Partitioner::Partition);
#else
    p->SetMethod(Partitioner::INNER_KMEANS,Partitioner::Partition);
#endif
    p->Evaluate();
    delete p;
    mesh_init->Redistribute();
    mesh_init->ReorderEmpty(CELL | FACE | EDGE | NODE);
    mesh_init->ExchangeGhost(1, NODE);
#endif
   // mesh_init->Save("Save_Stokes_init.pmf");

    if (processRank == 0) {
        std::cout << " repartitioning success" << std::endl;
    }

    Mesh *mesh_main;
    RepartitionStatistics(mesh_init, std::cout);
    BARRIER
    if (processRank == 0) {
        std::cout << " " << std::endl;
    }
    if(refine){
        mesh_main = new Mesh();
#if defined(USE_MPI)
        mesh_main->SetCommunicator(INMOST_MPI_COMM_WORLD);
#endif
        if (processRank == 0) {
            std::cout<<" refine mesh on "<<mesh_init->GetProcessorsNumber()<< " processors "<<refineNumber<<" times"<<std::endl;
        }
        RefineLocalMesh(mesh_init, mesh_main, refineNumber);
        delete mesh_init;

        mesh_main->ResolveShared();
        mesh_main->ExchangeGhost(1, NODE); // Construct Ghost cells in 1 layers connected via nodes
       // RepartitionStatistics(mesh_main, std::cout);

    }
    else
    {
        mesh_main = mesh_init;
    }

   // return 0;
    if (processRank == 0) {
        std::cout << "Before creating discretization" << std::endl;
    }
   // PrintCurrentMemory(std::cout,mesh_main);

    Tag Label_tag = mesh_main->GetTag(LabelTag);
    M_Assert(Label_tag.isValid(), "We have no valid Label tag");

    ////mark BC for points
    //   Tag BC_Edge_tag_init = mesh_init->CreateTag(EdgeTagBC,DATA_INTEGER,EDGE,EDGE,1);
    //mark nodes due to faces
    //see Ani3D MultiPackage Stokes example
    for (Mesh::iteratorEdge it = mesh_main->BeginEdge(); it != mesh_main->EndEdge(); it++) {
        if (it->GetStatus() != Element::Ghost) {
            it->Integer(Label_tag) = No_BC_mark;
        }
    }
    //mark nodes due to faces
    //see Ani3D MultiPackage Stokes example
    for (Mesh::iteratorNode it = mesh_main->BeginNode(); it != mesh_main->EndNode(); it++) {
        if (it->GetStatus() != Element::Ghost) {
            it->Integer(Label_tag) = No_BC_mark;
        }
    }


    for (Mesh::iteratorFace it = mesh_main->BeginFace(); it != mesh_main->EndFace(); it++)
        if ((it->Integer(Label_tag) != No_BC_mark)) {
            ElementArray<Node> nodes = it->getNodes();

            if (it->Integer(Label_tag) != 5) {
                for (i = 0; i < 3; i++)
                    if (nodes[i].GetStatus() != Element::Ghost)
                    {
                    if (nodes[i].Integer(Label_tag) == No_BC_mark) {
                        nodes[i].Integer(Label_tag) = it->Integer(Label_tag);
                    } else {
                        if (nodes[i].Integer(Label_tag) != it->Integer(Label_tag))
                            nodes[i].Integer(Label_tag) = 9;
                    }
                }
            }
        }

    mesh_main->ExchangeData(Label_tag,NODE,0);
    mesh_main->ExchangeData(Label_tag,EDGE,0);

    //mesh_main->Save("Save_Stokes_mesh.pmf");
    //return 0;
    BARRIER
    t_init = Timer();


   // mesh_main->Save("test_new.pmf");
   // PrintCurrentMemory(std::cout,mesh_main->GetProcessorRank());

    Ani_discretization *discr =  new Ani_discretization(mesh_main);
//    std::cout<<"Processor rank "<< mesh_main->GetProcessorRank()<<" marking success"<<std::endl;
    discr->PrepareProblemAni(4,3,0,0,ANITYPE,REVERSE);
    Sparse::Matrix mat("A"); // Declare the matrix of the linear system to be solved
    Sparse::Vector b("rhs"); // Declare the right-hand side vector
    Sparse::Vector x("sol"); // Declare the solution vector
    x.SetInterval(discr->getMinInd(),discr->getMaxInd()+1);
    double *DDATA = new double[5];
    double visc_coeff = 1.0;
    DDATA[0] = visc_coeff;

    int *IDATA = new int[5];
    BARRIER
    t_prepare = Timer();
  //  PrintCurrentMemory(std::cout,mesh_main->GetProcessorRank());

    discr->Assemble( Label_tag, mat, b, DDATA, 5, IDATA, 5);
    if (processRank == 0) {
        std::cout<<"assembling success"<<std::endl;
    }
    BARRIER
    t_assmble = Timer();
   // if(mesh_main->GetProcessorsNumber() > 1)
        SaveMatrix(mat,discr->MaxInd,INMOST_MPI_COMM_WORLD);
        b.Save("b.rhs");

//    PrintCurrentMemory(std::cout,mesh_main);
    if (processRank == 0) {
        std::cout << "set matrix" << std::endl;
    }

    s.SetMatrix(mat);

    //BARRIER
    if (processRank == 0) {
        std::cout << "solving" << std::endl;
    }
    if(!s.Solve(b,x)){
        std::cout<<"solution failed, residual  "<<s.Residual()<<" iterations "<<s.Iterations()<<" name "<<std::endl;
        return -1;
    }
    BARRIER
    t_solve=Timer();
    //std::cout<<"solution ended, residual  "<<s.Residual()<<" iterations "<<s.Iterations()<<std::endl;
    if (processRank == 0) {
        std::cout<<"solution finished on "<< mesh_main->GetProcessorsNumber() <<"processors, residual  "<<s.Residual()<<" iterations "<<s.Iterations()<<std::endl;

        std::cout << "init " << t_init - t0 << " preparing " << t_prepare - t_init << " assembling "
                  << t_assmble - t_prepare << " solving " << t_solve - t_assmble << std::endl;
    }
    Tag P_sol = mesh_main->CreateTag("P_Sol",DATA_REAL,NODE,NODE,1);
    for (Mesh::iteratorNode it = mesh_main->BeginNode(); it != mesh_main->EndNode(); it++) {
        if (it->GetStatus() != Element::Ghost)
        {
            it->Real(P_sol) = x[it->IntegerArray(discr->NumTags[0])[3]];
        }
    }
    mesh_main->ExchangeData(P_sol,NODE,0);
    std::vector<Tag> tags; tags.push_back(P_sol);


    WriteTags("test_init.vtk",mesh_main,tags);


    //return 0;




  //  if (processRank == 0) {
 //       std::cout<<"solution finished on "<< mesh_main->GetProcessorsNumber() <<"processors, residual  "<<s.Residual()<<" iterations "<<s.Iterations()<<std::endl;
//
      //  std::cout << "init " << t_init - t0 << " preparing " << t_prepare - t_init << " assembling "
 //                 << t_assmble - t_prepare << " solving " << t_solve - t_assmble << std::endl;
  //  }
 //   memory_test = getPeakRSS();
 //   std::cout<<"Processor rank "<< mesh_init->GetProcessorRank()<<" peak memory usage "<<memory_test
 //            << " B, or "<<memory_test/(1024L) <<" KB, or "<<memory_test/(1024L*1024L)<<" MB, or "<<memory_test/(1024L*1024L* 1024L) <<" GB "<<std::endl;
  //  PrintPeakMemory(std::cout,mesh_main);

    Solver::Finalize();
    Partitioner::Finalize();
    Mesh::Finalize();
    return 0;
}

