#include <utility>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include "inmost.h"
#include "../inmost_ani_lib/utils.h"
#include "../inmost_ani_lib/inmost_ani_fem.h"
#define USE_MPI
#define USE_PARTITIONER
using namespace INMOST;
#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif
extern "C" {
 int dbc_(double &x, double &y, double &z, int &label, double *DDATA,int *IDATA,int *iSYS, double* eBC);
};




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
    p->SetMethod(Partitioner::Inner_RCM,Partitioner::Partition);
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
    for (Mesh::iteratorNode it = mesh_main->BeginNode(); it != mesh_main->EndNode(); it++) {
        if (it->GetStatus() != Element::Ghost) {
            if (it->Coords()[0] == 0) it->Integer(Label_tag) = 1;
            if ((it->Coords()[0] == 0) && (it->Coords()[1] >= 0.25) && (it->Coords()[1] <= 0.75) &&
                (it->Coords()[2] >= 0.25) && (it->Coords()[2] <= 0.75))
                it->Integer(Label_tag) = 2;
        }
    }
    mesh_main->ExchangeData(Label_tag, NODE, 0);
    mesh_main->ExchangeData(Label_tag, EDGE, 0);
    //////////////////////////////////create discretization/////////////////////////////

    Ani_discretization *discr = new Ani_discretization(mesh_main);
    if (mesh_main->GetProcessorRank() == 0)
        std::cout << "before create Discr" << std::endl;

    discr->PrepareProblemAni(1, 0, 0, 0, ANITYPE, REVERSE);
    if (processRank == 0) {
        std::cout << " preparation success" << std::endl;
    }
    int sol_ind = 0;
    int prev_sol_ind = 1;
    //c time integration step
    double DeltaT = 0.015;
    int Itime = 0;
    double Time = 0;
    double FinalTime = 0.45;
    double Velocity[3] = {1, 0, 0};
    double D_coeff = 1e-4;

    Tag Sol_tag = mesh_main->CreateTag("U", DATA_REAL, NODE, NODE, 2);
    Tag Sol_tags[2];// = mesh_main->CreateTag("U", DATA_REAL, NODE, NODE, 2);
    Sol_tags[0] = mesh_main->CreateTag("U0", DATA_REAL, NODE, NODE, 1);
    Sol_tags[1] = mesh_main->CreateTag("U1", DATA_REAL, NODE, NODE, 1);


    double *DDATA = new double[5 + mesh_main->NumberOfCells() + mesh_main->NumberOfNodes()];
    DDATA[0] = D_coeff;// ! D
    DDATA[1] = Velocity[0];//    ! v_x
    DDATA[2] = Velocity[1];//    ! v_y
    DDATA[3] = Velocity[2];//   ! v_y
    DDATA[4] = DeltaT;

    for (Mesh::iteratorNode it = mesh_main->BeginNode(); it != mesh_main->EndNode(); it++)
        if (it->GetStatus() != Element::Ghost) {
            //std::cout<<it->GlobalID()<<" "<<sol_ind<<" "<<prev_sol_ind<< std::endl;
            if (it->Integer(Label_tag) > 0) {
                double dd[2], ec[2];
                int id[2];
                double x0 = it->Coords()[0];
                double y = it->Coords()[1];
                double z = it->Coords()[2];
                int lb = it->Integer(Label_tag);
                int iSYS[21];
                dbc_(x0, y, z, lb, dd, id, iSYS, ec);
                it->Real(Sol_tags[sol_ind]) = ec[0];
                it->Real(Sol_tags[prev_sol_ind]) = ec[0];
            } else {
                it->Real(Sol_tags[sol_ind]) = 0.0;
                it->Real(Sol_tags[prev_sol_ind]) = 0.0;
            }
        }
    mesh_main->ExchangeData(Sol_tag, NODE, 0);
    mesh_main->ExchangeData(Sol_tags[0], NODE, 0);
    mesh_main->ExchangeData(Sol_tags[1], NODE, 0);
    //std::cout<<"end set IC"<<std::endl;
    int *IDATA = new int[5];
    IDATA[0] = 5 + mesh_main->NumberOfCells() + mesh_main->NumberOfNodes();

    sum_t_prepare = 0.0;
    sum_t_assmble = 0.0;
    sum_t_solve = 0.0;
    sum_t_exch = 0.0;
    BARRIER
    t_prepare = Timer();

    do {
        BARRIER
        t_begin_it = Timer();
        Solver solv = Solver(solverName, "test");

        Sparse::Matrix mat("A"); // Declare the matrix of the linear system to be solved
        Sparse::Vector b("rhs"); // Declare the right-hand side vector
        Sparse::Vector x("sol"); // Declare the solution vector
        x.SetInterval(discr->getMinInd(), discr->getMaxInd() + 1);


        Itime++;

        Time += DeltaT;
        for (Mesh::iteratorCell it = mesh_main->BeginCell(); it != mesh_main->EndCell(); it++) {
            ElementArray<Node> nodes = it->getNodes();
            DDATA[5 + it->DataLocalID()] = 0.0;
            for (i = 0; i < 4; i++) {
                DDATA[5 + it->DataLocalID()] +=
                        nodes[i].Real(Sol_tags[sol_ind]) / 2.0 - nodes[i].Real(Sol_tags[prev_sol_ind]) / 8.0;
            }
        }

        for (Mesh::iteratorNode it = mesh_main->BeginNode(); it != mesh_main->EndNode(); it++) {
            int ind = it->IntegerArray(discr->NumTags[0])[0];
            DDATA[5 + mesh_main->NumberOfCells() + it->DataLocalID()] =
                    -(it->Real(Sol_tags[prev_sol_ind]) - 4 * it->Real(Sol_tags[sol_ind])) / 3.0;
        }
        BARRIER
        t_prepare_it = Timer();

        discr->Assemble(Label_tag, mat, b, DDATA, IDATA[0], IDATA, 5);

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
        BARRIER
        t_solve = Timer();
        if (processRank == 0) {
            std::cout << processorsCount << " processors for Solver " << solverName;
            std::cout << " with " << solv.Iterations() << " iterations to " << solv.Residual() << " norm" << std::endl;
        }

        for (Mesh::iteratorNode it = mesh_main->BeginNode(); it != mesh_main->EndNode(); it++)
            if (it->GetStatus() != Element::Ghost) {
                if (it->Integer(Label_tag) > 0) {
                    double dd[2], ec[2];
                    int id[2];
                    double x0 = it->Coords()[0];
                    double y = it->Coords()[1];
                    double z = it->Coords()[2];
                    int lb = it->Integer(Label_tag);
                    int iSYS[21];
                    dbc_(x0, y, z, lb, dd, id, iSYS, ec);
                    it->Real(Sol_tags[prev_sol_ind]) = ec[0];
                } else {
                    it->Real(Sol_tags[prev_sol_ind]) = x[it->IntegerArray(discr->NumTags[0])[0]];
                }
            }
        mesh_main->ExchangeData(Sol_tags[prev_sol_ind], NODE, 0);

        i = sol_ind;
        sol_ind = prev_sol_ind;
        prev_sol_ind = i;

        solv.Clear();

        BARRIER
        t_exch = Timer();
        if (processRank == 0) {
            std::cout << "iteration " << Itime << " t_prep "
                      << t_prepare_it - t_begin_it << " t_assmble " << t_assmble - t_prepare_it << " t_solve "
                      << t_solve - t_assmble << " T_exchange " << t_exch - t_solve << std::endl;
        }

        sum_t_prepare += t_prepare_it - t_begin_it;
        sum_t_assmble += t_assmble - t_prepare_it;
        sum_t_solve += t_solve - t_assmble;
        sum_t_exch += t_exch - t_solve;
    } while (Time < FinalTime);


    if (processRank == 0) {
        std::cout << std::endl;
        std::cout << "End of program, Time init " << t_init - t0 << " t_refine " << t_refine - t_init << " t_prepare "
                  << t_prepare - t_refine << std::endl;
        std::cout << "iterations number " << Itime << " t_prep "
                  << sum_t_prepare << " t_assmble " << sum_t_assmble << " t_solve "
                  << sum_t_solve << " T_exchange " << sum_t_exch << std::endl;
    }

    std::vector<Tag> tags2;
    tags2.push_back(Sol_tags[sol_ind]);
    if (mesh_main->GetProcessorsNumber() == 1) {
        WriteTags_vtk("test_writer.vtk", mesh_main, tags2);

    } else {
        WriteTags_pvtk("test_pvtk_writer.pvtk", mesh_main, tags2);
    }

    //mesh_main->Save("save_sol.pvtk");

    return 0;



}

