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
#include "../inmost_ani_lib/inmost_ani_fem.h"
#define USE_MPI
#define USE_PARTITIONER
using namespace INMOST;
#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif


//#define DEBUG_FORTRAN_NUMERATION
extern "C" {
    void drotate_(double *D, int *i, double *theta);
    void drecover_(double *vrt,int * lbE,double *  DDATA,int * IDATA,double* Sol,double * SolCell,double *SolVCell);


//      Subroutine HMFEMrecoverUP(nv, vrt, nt, tet, material, DDATA,IDATA,
//    &                          SOL, SOLatTets, VECatTets, iW, MaxWi)

//    void hmfemrecoverup_(int *nv ,double *vrt,int *nt, int *tet,int *material,double *  DDATA,int * IDATA,
//                         double *SOL, double *SOLatTets, double *VECatTets, int *iW, int *MaxWi);
/*
void refine_(int *number_of_refines ,int *nv ,int *nvmax ,int *nt ,int *ntmax ,int *nb ,int *nbmax,
             double *vrt, int *tet, int *bnd, int *material, int *labelF,
             double *MapMtr, int* Ref2MapMtr, int *iW, int *MaxWi,int* flag_write);

void loadmani_(
        int * nvmax,int * nbmax,int * ntmax,
        int * nv,int * nb,int * nt,
        double * vrt,int * bnd,int * tet,int * labelF,int * material,
        int *nvfix,int * nbfix,int * ntfix,int * ivfix,int * ibfix,int * itfix,
        int *iW,int * iW1, const char * name);
*/
}


int main(int argc, char *argv[]) {
    int i, j, k;
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
    Solver s = Solver(solverName, "Solver");

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
    if (mesh_init->GetProcessorRank() == 0) {
        //  LoadAFT(MeshName,mesh_init);
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

    Mesh *mesh_main;
    RepartitionStatistics(mesh_init, std::cout);

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
        RepartitionStatistics(mesh_main, std::cout);
    } else {
        mesh_main = mesh_init;
    }

    // return 0;
    Tag Label_tag = mesh_main->GetTag(LabelTag);
    M_Assert(Label_tag.isValid(), "We have no valid Label tag");
    Tag Fortran_num_tag = mesh_main->GetTag(FortranNumTag);
    M_Assert(Fortran_num_tag.isValid(), "We have no valid Fortran_num_tag tag");

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


//  return 0;

    for (Mesh::iteratorNode it = mesh_main->BeginNode(); it != mesh_main->EndNode(); it++) {
        if (it->GetStatus() != Element::Ghost) {
            it->Integer(Label_tag) = 0;
        }
    }

    mesh_main->ExchangeData(Label_tag, NODE, 0);
    mesh_main->ExchangeData(Label_tag, EDGE, 0);

    //return 0;
    BARRIER
    t_init = Timer();
    Ani_discretization *discr = new Ani_discretization(mesh_main);
//    std::cout<<"Processor rank "<< mesh_main->GetProcessorRank()<<" marking success"<<std::endl;
    discr->PrepareProblemAni(0, 0, 1, 0, ANITYPE, STRAIGHT);
    Sparse::Matrix mat("A"); // Declare the matrix of the linear system to be solved
    Sparse::Vector b("rhs"); // Declare the right-hand side vector
    Sparse::Vector x("sol"); // Declare the solution vector

    x.SetInterval(discr->getMinInd(), discr->getMaxInd() + 1);

    double *DDATA = new double[9];
    int num_ax;
    double theta;
    for (i = 0; i < 9; i++) DDATA[i] = 0.0;
    DDATA[0] = 100.0;
    DDATA[4] = 10.0;
    DDATA[8] = 1.0;

    theta = M_PI / 3.0;
    num_ax = 1;
    drotate_(DDATA, &num_ax, &theta);

    theta = M_PI / 4.0;
    num_ax = 2;
    drotate_(DDATA, &num_ax, &theta);

    theta = M_PI / 6.0;
    num_ax = 3;
    drotate_(DDATA, &num_ax, &theta);

    std::cout << "test rotate " << std::endl;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            std::cout << DDATA[i * 3 + j] << " ";
        std::cout << std::endl;
    }
    int *IDATA = new int[1];
    IDATA[0] = 9;
    BARRIER
    t_prepare = Timer();
    discr->Assemble(Label_tag, mat, b, DDATA, 9, IDATA, 5);
    if (processRank == 0) {
        std::cout << " assembling success" << std::endl;
    }
    BARRIER
    t_assmble = Timer();
    //  if(mesh_main->GetProcessorsNumber() > 1)
    //    SaveMatrix(mat,discr->MaxInd,INMOST_MPI_COMM_WORLD);

    mat.Save("mat_hybr_cpp.mtx");
  //  return 0;
    //  if (processRank == 0) {
    //std::cout<<"try to solve main task"<<std::endl;}
    //BARRIER
    std::cout << "Processor rank " << mesh_main->GetProcessorRank() << " minInd " << discr->MinInd << " maxInd "
              << discr->MaxInd << " MatrSize " << discr->MatrSize << std::endl;
    if (processRank == 0) {
        std::cout << "set matrix" << std::endl;
    }
    s.SetMatrix(mat);

    //BARRIER
    if (processRank == 0) {
        std::cout << "solving" << std::endl;
    }
    if (!s.Solve(b, x)) {
        std::cout << "solution failed, residual  " << s.Residual() << " iterations " << s.Iterations() << " name "
                  << std::endl;
        return -1;
    }
    BARRIER
    t_solve = Timer();
    //std::cout<<"solution ended, residual  "<<s.Residual()<<" iterations "<<s.Iterations()<<std::endl;
    if (processRank == 0) {
        std::cout << "init " << t_init - t0 << " preparing " << t_prepare - t_init << " assembling "
                  << t_assmble - t_prepare << " solving " << t_solve - t_assmble << std::endl;
    }

    x.Save("test_sol.rhs");
    Tag SolutionFace = mesh_main->CreateTag("SolFace", DATA_REAL, FACE, FACE, 1);

    for (Mesh::iteratorFace it = mesh_main->BeginFace(); it != mesh_main->EndFace(); it++)
        if (it->GetStatus() != Element::Ghost) {
//            std::cout<<it->IntegerArray(discr->NumTags[2])[0]<<" "<<it->Integer(discr->NumTags[2])<<std::endl;
            it->Real(SolutionFace) = x[it->IntegerArray(discr->NumTags[2])[0]];
        }
    Tag Fortran_tag = mesh_main->GetTag(FortranNumTag);
    std::vector<Tag> tags;
    tags.push_back(SolutionFace);
    tags.push_back(Fortran_tag);

    WriteTags("test_init.vtk", mesh_main, tags);


    ///////recover solution/////
    Solver::Finalize();
    Partitioner::Finalize();
    Mesh::Finalize();
    return 0;

}

