#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "inmost.h"
using namespace INMOST;

#ifndef M_PI
#define M_PI 3.141592653589
#endif

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

//shortcuts
typedef Storage::bulk bulk;
typedef Storage::real real;
typedef Storage::integer integer;
typedef Storage::enumerator enumerator;
typedef Storage::real_array real_array;
typedef Storage::var_array var_array;



//#define OPTIMIZATION

int main(int argc,char ** argv)
{
    Solver::Initialize(&argc,&argv,""); // Initialize the solver and MPI activity
#if defined(USE_PARTITIONER)
    Partitioner::Initialize(&argc,&argv); // Initialize the partitioner activity
#endif
    if( argc > 1 )
    {
        /*
        std::cout << "Row: " << sizeof(Sparse::Row) << std::endl;
        std::cout << "HessianRow: " << sizeof(Sparse::HessianRow) << std::endl;
        std::cout << "Matrix: " << sizeof(Sparse::Matrix) << std::endl;
        std::cout << "HessianMatrix: " << sizeof(Sparse::HessianMatrix) << std::endl;
        std::cout << "variable: " << sizeof(multivar_expression) << std::endl;
        Sparse::Row J;
        Sparse::HessianRow H;
        unknown a = unknown(0.5,1);
        unknown b = unknown(0.1,2);
        (a*a*b*(a+b)).GetHessian(1.0,J,1.0,H);
        std::cout << "J: "; J.Print();
        std::cout << "H: "; H.Print();
        J.Clear();
        (a*a*b*(a+b)).GetJacobian(1.0,J);
        std::cout << "J: "; J.Print();
        return 0;
         */
        double ttt; // Variable used to measure timing
        bool repartition = false; // Is it required to redistribute the mesh?
        Mesh * m = new Mesh(); // Create an empty mesh
        { // Load the mesh
            ttt = Timer();
            m->SetCommunicator(INMOST_MPI_COMM_WORLD); // Set the MPI communicator for the mesh
            if( m->GetProcessorRank() == 0 ) std::cout << "Processors: " << m->GetProcessorsNumber() << std::endl;
            if( m->isParallelFileFormat(argv[1]) ) //The format is
            {
                m->Load(argv[1]); // Load mesh from the parallel file format
                repartition = true; // Ask to repartition the mesh
            }
            else if( m->GetProcessorRank() == 0 ) m->Load(argv[1]); // Load mesh from the serial file format
            BARRIER
            if( m->GetProcessorRank() == 0 ) std::cout << "Load the mesh: " << Timer()-ttt << std::endl;
        }
        
        
#if defined(USE_PARTITIONER)
        if (m->GetProcessorsNumber() > 1 )//&& !repartition) // Currently only non-distributed meshes are supported by Inner_RCM partitioner
        {
            { // Compute mesh partitioning
                ttt = Timer();
                Partitioner p(m); //Create Partitioning object
                p.SetMethod(Partitioner::Inner_RCM,repartition ? Partitioner::Repartition : Partitioner::Partition); // Specify the partitioner
                p.Evaluate(); // Compute the partitioner and store new processor ID in the mesh
                BARRIER
                if( m->GetProcessorRank() == 0 ) std::cout << "Evaluate: " << Timer()-ttt << std::endl;
            }
            
            { //Distribute the mesh
                ttt = Timer();
                m->Redistribute(); // Redistribute the mesh data
                m->ReorderEmpty(CELL|FACE|EDGE|NODE); // Clean the data after reordring
                BARRIER
                if( m->GetProcessorRank() == 0 ) std::cout << "Redistribute: " << Timer()-ttt << std::endl;
            }
        }
#endif
        
        { // prepare geometrical data on the mesh
            ttt = Timer();
            Mesh::GeomParam table;
            table[CENTROID]    = CELL | FACE; //Compute averaged center of mass
            table[NORMAL]      = FACE;        //Compute normals
            table[ORIENTATION] = FACE;        //Check and fix normal orientation
            table[MEASURE]     = CELL | FACE; //Compute volumes and areas
            //table[BARYCENTER]  = CELL | FACE; //Compute volumetric center of mass
            m->PrepareGeometricData(table); //Ask to precompute the data
            BARRIER
            if( m->GetProcessorRank() == 0 ) std::cout << "Prepare geometric data: " << Timer()-ttt << std::endl;
        }
        
        // data tags for
        Tag tag_P;  // Pressure
        Tag tag_K;  // Diffusion tensor
        Tag tag_F;  // Forcing term
        Tag tag_BC; // Boundary conditions
        Tag tag_W;  // Gradient matrix acting on harmonic points on faces and returning gradient on faces
        Tag tag_DMP; // Indicates weather local W matrix satisfy DMP condition
        
        if( m->GetProcessorsNumber() > 1 ) //skip for one processor job
        { // Exchange ghost cells
            ttt = Timer();
            m->ExchangeGhost(1,FACE); // Produce layer of ghost cells
            BARRIER
            if( m->GetProcessorRank() == 0 ) std::cout << "Exchange ghost: " << Timer()-ttt << std::endl;
        }
        
        { //initialize data
            if( m->HaveTag("PERM") ) // is diffusion tensor already defined on the mesh? (PERM from permeability)
                tag_K = m->GetTag("PERM"); // get the diffusion tensor
            
            if( !tag_K.isValid() || !tag_K.isDefined(CELL) ) // diffusion tensor was not initialized or was not defined on cells.
            {
                tag_K = m->CreateTag("PERM",DATA_REAL,CELL,NONE,6); // create a new tag for symmetric diffusion tensor K
                for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) ) // loop over mesh cells
                {
                    Cell cell = m->CellByLocalID(q);
                    real_array K = cell->RealArray(tag_K);
                    // assign a symmetric positive definite tensor K
                    K[0] = 1.0; //XX
                    K[1] = 0.0; //XY
                    K[2] = 0.0; //XZ
                    K[3] = 1.0; //YY
                    K[4] = 0.0; //YZ
                    K[5] = 1.0; //ZZ
                }
                
                m->ExchangeData(tag_K,CELL,0); //Exchange diffusion tensor
            }
            
            if( m->HaveTag("PRESSURE") ) //Is there a pressure on the mesh?
                tag_P = m->GetTag("PRESSURE"); //Get the pressure
            
            if( !tag_P.isValid() || !tag_P.isDefined(CELL) ) // Pressure was not initialized or was not defined on nodes
            {
                srand(1); // Randomization
                tag_P = m->CreateTag("PRESSURE",DATA_REAL,CELL|FACE,NONE,1); // Create a new tag for the pressure
                for(Mesh::iteratorElement e = m->BeginElement(CELL|FACE); e != m->EndElement(); ++e) //Loop over mesh cells
                    e->Real(tag_P) = 0;//(rand()*1.0)/(RAND_MAX*1.0); // Prescribe random value in [0,1]
            }
            
            if( !tag_P.isDefined(FACE) )
            {
                tag_P = m->CreateTag("PRESSURE",DATA_REAL,FACE,NONE,1);
                for(Mesh::iteratorElement e = m->BeginElement(FACE); e != m->EndElement(); ++e) //Loop over mesh cells
                    e->Real(tag_P) = 0;//(rand()*1.0)/(RAND_MAX*1.0); // Prescribe random value in [0,1]
            }
            
            
            
            if( m->HaveTag("BOUNDARY_CONDITION") ) //Is there boundary condition on the mesh?
            {
                tag_BC = m->GetTag("BOUNDARY_CONDITION");
                //initialize unknowns at boundary
            }
            m->ExchangeData(tag_P,CELL|FACE,0); //Synchronize initial solution with boundary unknowns
            tag_W = m->CreateTag("nKGRAD",DATA_REAL,CELL,NONE);
            tag_DMP = m->CreateTag("isDMP",DATA_BULK,CELL,NONE);
            ttt = Timer();
            //Assemble gradient matrix W on cells
            int total = 0, dmp = 0;
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:total) reduction(+:dmp)
#endif
            for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
            {
                Cell cell = m->CellByLocalID(q);
                real xP[3]; //center of the cell
                real yF[3]; //center of the face
                real nF[3]; //normal to the face
                real aF; //area of the face
                real vP = cell->Volume(); //volume of the cell
                cell->Centroid(xP);
                ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
                int NF = (int)faces.size(); //number of faces;
                rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),cell->RealArrayDF(tag_K).size()); //get permeability for the cell
                //rMatrix U,S,V;
                //K0.SVD(U,S,V);
                //for(int k = 0; k < 3; ++k) S(k,k) = sqrt(S(k,k));
                //rMatrix K = U*S*V;
                rMatrix nKGRAD(NF,NF), NK(NF,3), R(NF,3), D(NF,NF), U(NF,NF), Areas(NF,1); //big gradient matrix, co-normals, directions
                for(int k = 0; k < NF; ++k) //loop over faces
                {
                    aF = faces[k].Area();
                    faces[k].Centroid(yF);
                    faces[k].OrientedUnitNormal(cell->self(),nF);
                    // assemble matrix of directions
                    R(k,0) = (yF[0]-xP[0])*aF;
                    R(k,1) = (yF[1]-xP[1])*aF;
                    R(k,2) = (yF[2]-xP[2])*aF;
                    // assemble matrix of co-normals
                    rMatrix nK = rMatrix::FromVector(nF,3).Transpose()*K;
                    NK(k,0) = nK(0,0);
                    NK(k,1) = nK(0,1);
                    NK(k,2) = nK(0,2);
                } //end of loop over faces
                nKGRAD = NK*(NK.Transpose()*R).Invert(true).first*NK.Transpose(); //stability part
				nKGRAD += (rMatrix::Unit(NF) - R*(R.Transpose()*R).Invert(true).first*R.Transpose())*(2.0/(static_cast<real>(NF)*vP)*(NK*K.Invert(true).first*NK.Transpose()).Trace());


                /*
				 rMatrix SU,SS,SV;
                 std::cout << "W" << std::endl;
                 nKGRAD.Print();
                 nKGRAD.SVD(SU,SS,SV);
                 std::cout << "U" << std::endl;
                 SU.Print();
                 std::cout << "S" << std::endl;
                 SS.Print();
                 std::cout << "V" << std::endl;
                 SV.Print();
                 std::cout << "Check " << (nKGRAD - SU*SS*SV.Transpose()).FrobeniusNorm() << std::endl;
                 */
                /*
                int rank = 0; //size of matrix U
                
                { //Retrive orthogonal to R matrix D
                    //Symmetric orthogonal matrix
                    rMatrix DUD = (rMatrix::Unit(NF) - R*(R.Transpose()*R).Invert(true).first*R.Transpose());
                    //perfrom singular value decomposition
                    //S should be unity matrix with rank NF-3
                    rMatrix DUD_U,DUD_S,DUD_V;
                    DUD.SVD(DUD_U,DUD_S,DUD_V);
                    //compute the rank
                    for(int q = 0; q < NF; ++q)
                        if( DUD_S(q,q) > 1.0e-2 )
                            rank++;
                    rank = NF-3;
                    if( rank != NF-3)
                    {
                        std::cout << "rank: " << rank << " expected " << NF-3 << std::endl;
                        DUD_S.Print();
                    }
                    //chop matrix to the full rank
                    DUD_S.RemoveSubset(rank,NF,rank,NF);
                    DUD_V.RemoveColumns(rank,NF);
                    //assign the matrix
                    D = DUD_V;
                    U = DUD_S;
                }
                //std::cout << "D" << std::endl;
                //D.Print();
                //std::cout << "U" << std::endl;
                //U.Print();
                //std::cout << "DtR" << std::endl;
                //(D.Transpose()*R).Print();
                
                U *=(2.0/(static_cast<real>(NF)*vP)*(NK*K.Invert(true).first*NK.Transpose()).Trace());
#if defined(OPTIMIZATION)
                { //Make W a Z-matrix
                    vMatrix vL(rank,rank), vD(rank,rank), vW(NF,NF);
                    int unk = 0;
                    // U = L*D*L^T
                    vD.Zero();
                    //diagonal D
                    for(int i = 0; i < rank; ++i)
                    {
                        vD(i,i) = unknown(U(i,i),unk);
                        unk++;
                    }
                    //off-diagonal
                    vL.Zero();
                    for(int i = 0; i < rank; ++i)
                        vL(i,i) = 1.0;
                    for(int i = 1; i < rank; ++i)
                        for(int j = 0; j < i; ++j)
                        {
                            vL(i,j) = unknown(0.0,unk);
                            unk++;
                        }
                    //std::cout << "unknowns: " << unk << std::endl;
                    variable phi,s;
                    //std::cout << "vD" << std::endl;
                    //vD.Print();
                    //std::cout << "vL" << std::endl;
                    //vL.Print();
                    
                    int iter = 0;
                    do
                    { //Optimize U matrix
                        vW = nKGRAD + D*vL*vD*vL.Transpose()*D.Transpose();
                        //construct minimization functional phi(W)
                        phi = 0.0;
                        for(int i = 0; i < NF; ++i)
                        {
                            //phi += 1.0 / (vW(i,i)*vW(i,i));
                            s = vW(i,i)*faces[i].Area();
                            for(int j = 0; j < NF; ++j) if( i != j )
                            {
                                s += vW(i,j)*faces[j].Area();
                                phi += (vW(i,j)+fabs(vW(i,j)))*(vW(i,j)+fabs(vW(i,j)));
                            }
                            phi += (s - fabs(s))*(s - fabs(s));
                        }
                        Sparse::Row & der = phi.GetRow(); //row of derivatives
                        //std::sort(der.Begin(),der.End());
                        //for(int i = 0; i < der.Size(); ++i)
                        //    std::cout << "(" << der.GetIndex(i) << "," << der.GetValue(i) << ") ";
                        //std::cout<<std::endl;
                        int q = 0;
                        real a = 0.00005;
                        real minvD = 1.0e20;
                        //diagonal
                        for(int i = 0; i < rank; ++i)
                        {
                            real d = a*der[q++];
                            //if( vD(i,i)-d > 0.0 )
                            vD(i,i) -= d;
                            if( vD(i,i) < minvD ) minvD = get_value(vD(i,i));
                        }
                        std::cout << "[" << iter << "] phi: " << get_value(phi) << " minD " << minvD << std::endl;
                        //off-diagonal
                        for(int i = 1; i < rank; ++i)
                            for(int j = 0; j < i; ++j)
                            {
                                real d = a*der[q++];
                                vL(i,j) -= d;
                            }
                        iter++;
                        //std::cout << "vD" << std::endl;
                        //vD.Print();
                        //std::cout << "vL" << std::endl;
                        //vL.Print();
                    } while(iter < 100 && phi > 1.0e-3);
                    {
                        vMatrix vU = vL*vD*vL.Transpose();
                        for(int i = 0; i < rank; ++i)
                            for(int j = 0; j < rank; ++j)
                                U(i,j) = get_value(vU(i,j));
                    }
                    //std::cout << "U: " << std::endl;
                    //U.Print();
                    
                    
                }
#endif
                //std::cout << "UDtR" << std::endl;
                //(U*D.Transpose()*R).Print();
                
                nKGRAD += D*U*D.Transpose();
                */
                //std::cout << "W: " << std::endl;
                //nKGRAD.Print();
                bulk & isDMP = cell->Bulk(tag_DMP);
                isDMP = 1;
                for(int k = 0; k < NF; ++k)
                {
                    real row_sum = 0;
                    if( nKGRAD(k,k) < 0.0 ) isDMP = 0;
                    for(int j = 0; j < NF; ++j)
                        row_sum += nKGRAD(k,j);
                    if( row_sum < 0.0 ) isDMP = 0;
                    for(int j = k+1; j < NF; ++j)
                        if( nKGRAD(k,j) > 0.0 )
                            isDMP = 0;
                }
                ++total;
                if( isDMP ) ++dmp;
                real_array W = cell->RealArrayDV(tag_W); //access data structure for gradient matrix in mesh
                W.resize(NF*NF); //resize the structure
                std::copy(nKGRAD.data(),nKGRAD.data()+NF*NF,W.data()); //write down the gradient matrix
            } //end of loop over cells
            std::cout << "Construct W matrix: " << Timer() - ttt << std::endl;
            std::cout << "Satisfy DMP: " << dmp << " out of " << total << std::endl;
            
            if( m->HaveTag("FORCE") ) //Is there force on the mesh?
            {
                tag_F = m->GetTag("FORCE"); //initial force
                assert(tag_F.isDefined(CELL)); //assuming it was defined on cells
            } // end of force
        } //end of initialize data
        
        std::cout << "Initialization done" << std::endl;
        
        
        integer nit = 0;
        ttt = Timer();
        
        { //Main loop for problem solution
            Automatizator aut(m); // declare class to help manage unknowns
            Automatizator::MakeCurrent(&aut);
            dynamic_variable P(aut,aut.RegisterDynamicTag(tag_P,CELL|FACE)); //register pressure as primary unknown
            aut.EnumerateDynamicTags(); //enumerate all primary variables
            
            std::cout << "Enumeration done" << std::endl;
            
            
            Residual R("",aut.GetFirstIndex(),aut.GetLastIndex());
            Sparse::LockService Locks(aut.GetFirstIndex(),aut.GetLastIndex());
            Sparse::AnnotationService Text(aut.GetFirstIndex(),aut.GetLastIndex());
            Sparse::Vector Update  ("",aut.GetFirstIndex(),aut.GetLastIndex()); //vector for update
            {//Annotate matrix
                for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
                {
                    Cell cell = m->CellByLocalID(q);
                    if( cell.GetStatus() != Element::Ghost )
                        Text.SetAnnotation(P.Index(cell),"Cell-centered pressure value");
                }
                for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
                {
                    Face face = m->FaceByLocalID(q);
                    if( face.GetStatus() != Element::Ghost )
                    {
                        if( tag_BC.isValid() && face.HaveData(tag_BC) )
                            Text.SetAnnotation(P.Index(face),"Pressure guided by boundary condition");
                        else
                            Text.SetAnnotation(P.Index(face),"Interface pressure");
                    }
                }
            }
            
            std::cout << "Matrix was annotated" << std::endl;
            
            do
            {
                R.Clear(); //clean up the residual
                double tttt = Timer();
                int total = 0, dmp = 0;
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:total) reduction(+:dmp)
#endif
                for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) ) //loop over cells
                {
                    Cell cell = m->CellByLocalID(q);
                    ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
                    int NF = (int)faces.size();
                    rMatrix nKGRAD(cell->RealArrayDV(tag_W).data(),NF,NF); //Matrix for gradient
                    bulk & isDMP = cell->Bulk(tag_DMP);
                    vMatrix pF(NF,1); //vector of pressure differences on faces
                    vMatrix FLUX(NF,1); //computed flux on faces
                    for(int k = 0; k < NF; ++k)
                        pF(k,0) = (P(faces[k]) - P(cell))*faces[k].Area();
                    FLUX = nKGRAD*pF; //fluxes on faces
                    //Change matrix by nonlinear correction for DMP
                    ++total;
                    
                    //if( !isDMP )
                    if( false )
                    {
                        const real var = 0;//0.25;
                        vMatrix vnKGRAD = nKGRAD;
                        vMatrix card(NF,1);
                        vMatrix dpF(NF,1);
                        variable div = 0;
                        //denominators of corrections
                        //this is sum of absolute values of differences
                        for(int k = 0; k < NF; ++k)
                        {
                            div += FLUX(k,0)*faces[k].Area();;
                            for(int j = 0; j < NF; ++j)
                            {
                                dpF(k,0) += soft_fabs(pF(j,0) - pF(k,0),1.0e-9);
                                card(k,0) += (soft_fabs(pF(j,0) - pF(k,0),1.0e-9)) / (soft_fabs(pF(j,0) - pF(k,0),1.0e-9)+1.0e-9);
                            }
                        }
                        //actual corrections
                        variable beta;
                        for(int k = 0; k < NF; ++k)
                        {
                            for(int j = k+1; j < NF; ++j) //if( nKGRAD(k,j) > 0.0 )
                            {
                                //beta = soft_fabs(FLUX(k,0),1.0e-9)/(dpF(k,0)) + soft_fabs(FLUX(j,0),1.0e-9)/(dpF(j,0));
                                beta = soft_max(soft_fabs(FLUX(k,0)-FLUX(j,0),1.0e-9)/card(k,0),soft_fabs(FLUX(j,0)-FLUX(k,0),1.0e-9)/card(j,0),1.0e-9)/soft_fabs(pF(k,0)-pF(j,0),1.0e-9);
                                //beta = soft_max(soft_fabs(FLUX(k,0),1.0e-9)/card(k,0),soft_fabs(FLUX(j,0),1.0e-9)/card(j,0),1.0e-9)/soft_fabs(pF(k,0)-pF(j,0),1.0e-9);
                                //beta = soft_max(soft_fabs(div,1.0e-9)/card(k,0),soft_fabs(div,1.0e-9)/card(j,0),1.0e-9)/soft_fabs(pF(k,0)-pF(j,0),1.0e-9);
                                beta = variation(beta,var);
                                vnKGRAD(k,j) -= (beta);
                                vnKGRAD(j,k) -= (beta);
                                vnKGRAD(k,k) += (beta);
                                vnKGRAD(j,j) += (beta);
                            }
                        }
                        FLUX = vnKGRAD*pF;
                        bool haveDMP = true;
                        for(int k = 0; k < NF; ++k)
                        {
                            real row_sum = 0;
                            if( vnKGRAD(k,k) < 0.0 ) haveDMP = false;
                            for(int j = 0; j < NF; ++j)
                                row_sum += get_value(vnKGRAD(k,j));
                            if( row_sum < 0.0 ) haveDMP = false;
                            for(int j = k+1; j < NF; ++j)
                                if( vnKGRAD(k,j) > 0.0 )
                                    haveDMP = false;
                        }
                        if( !haveDMP )
                        {
                            //std::cout << "Failed to correct, matrix:" << std::endl;
                            //nKGRAD.Print();
                        }
                        else ++dmp;
                    }
                    else ++dmp;
                    
                    if( cell.GetStatus() != Element::Ghost )
                    {
                        for(int k = 0; k < NF; ++k) //loop over faces of current cell
                            R[P.Index(cell)] += FLUX(k,0)*faces[k].Area();
                    }
                    for(int k = 0; k < NF; ++k) //loop over faces of current cell
                    {
                        if( faces[k].GetStatus() == Element::Ghost ) continue;
                        int index = P.Index(faces[k]);
                        Locks.Lock(index);
                        if( tag_BC.isValid() && faces[k].HaveData(tag_BC) )
                        {
                            real_array BC = faces[k].RealArray(tag_BC);
                            R[index] -= BC[0]*P(faces[k]) + BC[1]*FLUX(k,0) - BC[2];
                        }
                        else
                            R[index] -= FLUX(k,0);
                        Locks.UnLock(index);
                    }
                } //end of loop over cells
                
                std::cout << "Satisfy DMP: " << dmp << " out of " << total << std::endl;
                
                
                if( tag_F.isValid() )
                {
#if defined(USE_OMP)
#pragma omp parallel for
#endif
                    for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
                    {
                        Cell cell = m->CellByLocalID(q);
                        if( cell.GetStatus() == Element::Ghost ) continue;
                        if( cell->HaveData(tag_F) ) R[P.Index(cell)] += cell->Real(tag_F)*cell->Volume();
                    }
                }
                
                std::cout << "assembled in " << Timer() - tttt << "\t\t\t" << std::endl;
                
                
                R.Rescale();
                //R.GetJacobian().Save("jacobian.mtx",&Text);
                //R.GetResidual().Save("residual.mtx");
                
                
                std::cout << "Nonlinear residual: " << R.Norm() << "\t\t" << std::endl;
                
                if( R.Norm() < 1.0e-4 ) break;
                
                //Solver S(Solver::INNER_MPTILUC);
				Solver S(Solver::SUPERLU);
				S.SetParameterReal("relative_tolerance", 1.0e-14);
                S.SetParameterReal("absolute_tolerance", 1.0e-12);
                S.SetParameterReal("drop_tolerance", 1.0e-2);
                S.SetParameterReal("reuse_tolerance", 1.0e-4);
                S.SetMatrix(R.GetJacobian());
                //std::fill(Update.Begin(),Update.End(),0.0);
                if( S.Solve(R.GetResidual(),Update) )
                {
#if defined(USE_OMP)
#pragma omp parallel for
#endif
                    for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
                    {
                        Cell cell = m->CellByLocalID(q);
                        cell->Real(tag_P) -= Update[P.Index(cell)];
                    }
#if defined(USE_OMP)
#pragma omp parallel for
#endif
                    for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
                    {
                        Face face = m->FaceByLocalID(q);
                        face->Real(tag_P) -= Update[P.Index(face)];
                    }
                    m->ExchangeData(tag_P, CELL|FACE, 0);
                    {
                        std::stringstream str;
                        str << "iter" << nit;
                        if( m->GetProcessorsNumber() == 1 )
                            str << ".vtk";
                        else
                            str << ".pvtk";
                        m->Save(str.str());
                    }
                }
                else
                {
                    std::cout << "Unable to solve: " << S.GetReason() << std::endl;
                    break;
                }
                ++nit;
            } while( R.Norm() > 1.0e-4 && nit < 10); //check the residual norm
        }
        std::cout << "Solved problem in " << Timer() - ttt << " seconds with " << nit << " iterations " << std::endl;
        
        if( m->HaveTag("REFERENCE_SOLUTION") )
        {
            Tag tag_E = m->CreateTag("ERRROR",DATA_REAL,CELL,NONE,1);
            Tag tag_R = m->GetTag("REFERENCE_SOLUTION");
            real C, L2, volume;
            C = L2 = volume = 0.0;
            for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
            {
                Cell cell = m->CellByLocalID(q);
                real err = cell->Real(tag_P) - cell->Real(tag_R);
                real vol = cell->Volume();
                if( C < fabs(err) ) C = fabs(err);
                L2 += err*err*vol;
                volume += vol;
                cell->Real(tag_E) = err;
            }
            L2 = sqrt(L2/volume);
            std::cout << "Error on cells, C-norm " << C << " L2-norm " << L2 << std::endl;
            C = L2 = volume = 0.0;
            if( tag_R.isDefined(FACE) )
            {
                tag_E = m->CreateTag("ERRROR",DATA_REAL,FACE,NONE,1);
                for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
                {
                    Face face = m->FaceByLocalID(q);
                    real err = face->Real(tag_P) - face->Real(tag_R);
                    real vol = (face->BackCell()->Volume() + (face->FrontCell().isValid() ? face->FrontCell()->Volume() : 0))*0.5;
                    if( C < fabs(err) ) C = fabs(err);
                    L2 += err*err*vol;
                    volume += vol;
                    face->Real(tag_E) = err;
                }
                L2 = sqrt(L2/volume);
                std::cout << "Error on faces, C-norm " << C << " L2-norm " << L2 << std::endl;
            }
            else std::cout << "Reference solution was not defined on faces" << std::endl;
        }
        
        if( m->GetProcessorsNumber() == 1 )
            m->Save("out.vtk");
        else
            m->Save("out.pvtk");
        
        delete m; //clean up the mesh
    }
    else
    {
        std::cout << argv[0] << " mesh_file" << std::endl;
    }
    
#if defined(USE_PARTITIONER)
    Partitioner::Finalize(); // Finalize the partitioner activity
#endif
    Solver::Finalize(); // Finalize solver and close MPI activity
    return 0;
}
