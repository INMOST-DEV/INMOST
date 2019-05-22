#include <utility>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include "inmost.h"
#include "utils.h"
#include "inmost_ani_fem.h"
//#define USE_MPI
//#define USE_PARTITIONER
using namespace INMOST;
#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif





int main(int argc, char **argv) {
    int i, j, k;
    int processRank = 0, processorsCount = 1;

#if defined(USE_MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);  // Get the rank of the current process
    MPI_Comm_size(MPI_COMM_WORLD, &processorsCount); // Get the total number of processors used
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
                std::cout << "-am, --animesh <AniMesh name>" << std::endl;
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
        if (strcmp(argv[i], "--refine") == 0) {
            if (processRank == 0) {
                std::cout << "mesh will be refined " << argv[i + 1] << " times" << std::endl;
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
            std::cout << "Parameters not found, you can specify parameter file name with -d or --database option, "
                      << std::endl;

    if (!meshFound) {
        if (processRank == 0)
            std::cout
                    << "Mesh name not found , you can specify mesh name  with -am or --animesh option."
                    << std::endl;
        exit(-1);
    }
/////////////////////////////////////////////////////////////////////////////end get parameters/////////////////////////////////
    ////////////////////////////////////// prepare  INMOST //////////////////////////
    Mesh::Initialize(&argc, &argv);
    Partitioner::Initialize(&argc, &argv);
    Solver::Initialize(&argc, &argv, parametersFound ? parametersFileName.c_str() : NULL);

    //std::cout<<"parameter filename  "<<parametersFileName<<"  "<< (parametersFound ? parametersFileName.c_str() : NULL)<<std::endl;
    if (!Solver::isSolverAvailable(solverName)) {
        if (processRank == 0) {
            std::cout << "Solver " << solverName << " is not available" << std::endl;
        }
        Mesh::Finalize();
        Partitioner::Finalize();
        Solver::Finalize();
        exit(-1);
    }


    if (processRank == 0) {
        std::cout << "Solving with " << solverName << std::endl;
    }
    double t0, t_init, t_refine, t_prepare, t_begin_it, t_prepare_it, t_assmble, t_solve, t_exch;

    double sum_t_prepare, sum_t_assmble, sum_t_solve, sum_t_exch;
    BARRIER
    t0 = Timer();
    ////////////////////////////////////// load and repartition mesh //////////////////////////
    Mesh *mesh_init;
    mesh_init = new Mesh();
#if defined(USE_MPI)
    mesh_init->SetCommunicator(INMOST_MPI_COMM_WORLD);
#endif
    if (mesh_init->GetProcessorRank() == 0) {
        ReadInmostFromAniFile(mesh_init, MeshName);
        std::cout << " loading success of mesh " << MeshName << std::endl;
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

    if (processRank == 0) {
        std::cout << " repartitioning success" << std::endl;
    }
    BARRIER
    t_init = Timer();
    Mesh *mesh_main;

////////////////////////////////////////////////////refine mesh , if it is needed/////////////////////////////////////////////
    if (refine) {
        mesh_main = new Mesh();
#if defined(USE_MPI)
        mesh_main->SetCommunicator(INMOST_MPI_COMM_WORLD);
#endif
        if (processRank == 0) {
            std::cout << " refine mesh on " << mesh_init->GetProcessorsNumber() << " processors " << refineNumber
                      << " times" << std::endl;
        }

        RefineLocalMesh(mesh_init, mesh_main, refineNumber);

        delete mesh_init;
        mesh_main->ResolveShared();

        mesh_main->ExchangeGhost(1, NODE); // Construct Ghost cells in 1 layers connected via nodes
    } else {
        mesh_main = mesh_init;
    }
    BARRIER
    t_refine = Timer();
    /////////////////////////////////////////set Boundary conditions, an in Adi3D example////////////////////
    Tag Label_tag = mesh_main->GetTag(LabelTag);
    M_Assert(Label_tag.isValid(), "We have no valid Label tag");
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

    for (Mesh::iteratorFace it = mesh_main->BeginFace(); it != mesh_main->EndFace(); it++) {
        if (it->GetStatus() != Element::Ghost) {
            ElementArray<Node> nodes = it->getNodes();
            for(i=0;i<nodes.size();i++){
                nodes[i].Integer(Label_tag) = std::max(nodes[i].Integer(Label_tag),it->Integer(Label_tag));
            }
        }
    }

    mesh_main->ExchangeData(Label_tag, NODE, 0);
    mesh_main->ExchangeData(Label_tag, EDGE, 0);
    //////////////////////////////////create discretization/////////////////////////////

    Ani_discretization *discr = new Ani_discretization(mesh_main);
    if (mesh_main->GetProcessorRank() == 0)
        std::cout << "before create Discr" << std::endl;

    discr->PrepareProblemAni(1, 1, 0, 0, ANITYPE, REVERSE);
    if (processRank == 0) {
        std::cout << " preparation success" << std::endl;
    }
    int sol_ind = 0;
    int prev_sol_ind = 1;
    //c time integration step
    double DeltaT = 0.015;
    int Itime = 0;
    double Time = 0;
    double FinalTime = 0.090    ;
    double Velocity[3] = {3, 2, 1};
    double D_coeff = 1;
    Tag Sol_tag = mesh_main->CreateTag("U", DATA_REAL, NODE|EDGE, NODE, 1);
    double *DDATA = new double[4];

    DDATA[0] = Velocity[0];//    ! v_x
    DDATA[1] = Velocity[1];//    ! v_y
    DDATA[2] = Velocity[2];//   ! v_y
    DDATA[3] = D_coeff;// ! D
    int *IDATA = new int[1];
    IDATA[0] = 4;


    Solver solv = Solver(solverName, "test");

    Sparse::Matrix mat("A"); // Declare the matrix of the linear system to be solved
    Sparse::Vector b("rhs"); // Declare the right-hand side vector
    Sparse::Vector x("sol"); // Declare the solution vector
    x.SetInterval(discr->getMinInd(), discr->getMaxInd() + 1);
    discr->Assemble(Label_tag, mat, b, DDATA, IDATA[0], IDATA, 1);

    if (processRank == 0) {
        std::cout << " assembling success" << std::endl;
    }
    BARRIER
    t_assmble = Timer();

    mat.Save("A.mtx");
    b.Save("b.rhs");
    solv.SetMatrix(mat);

    //    BARRIER
    if (processRank == 0) {
        std::cout << "solving" << std::endl;
    }
    if (!solv.Solve(b, x)) {
        std::cout << "solution failed, residual  " << solv.Residual() << " iterations " << solv.Iterations()
                  << " name "
                  << std::endl;
        return 0;
    }
    if (processRank == 0) {
        std::cout << "solving success" << std::endl;
    }
    for (Mesh::iteratorNode it = mesh_main->BeginNode(); it != mesh_main->EndNode(); it++)
        if (it->GetStatus() != Element::Ghost) {

            {
                it->Real(Sol_tag) = x[it->IntegerArray(discr->NumTags[0])[0]];
            }
        }

    for (Mesh::iteratorEdge it = mesh_main->BeginEdge(); it != mesh_main->EndEdge(); it++)
        if (it->GetStatus() != Element::Ghost) {

            {
                it->Real(Sol_tag) = x[it->IntegerArray(discr->NumTags[1])[0]];
            }
        }

    mesh_main->ExchangeData(Sol_tag, NODE|EDGE, 0);


    // std::vector<Tag> tags;
    // tags.push_back(Sol_tags[sol_ind]);
    std::vector<Tag> tags2;
    tags2.push_back(Sol_tag);
    //tags2.push_back(Label_tag);
    WriteTags("test_writer", mesh_main, tags2);

    //mesh_main->Save("save_sol.pvtk");
    //return 0;

    return 0;



}

