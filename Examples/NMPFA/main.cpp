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
        Tag tag_W;  // Cell-wise gradient matrix
        
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
            
            if( m->HaveTag("BOUNDARY_CONDITION") ) //Is there boundary condition on the mesh?
            {
                tag_BC = m->GetTag("BOUNDARY_CONDITION");
                //initialize unknowns at boundary
            }
            m->ExchangeData(tag_P,CELL|FACE,0); //Synchronize initial solution with boundary unknowns
            tag_W = m->CreateTag("nKGRAD",DATA_REAL,CELL,NONE);
            ttt = Timer();
            //Assemble gradient matrix W on cells
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				//private memory
				rMatrix y(3,1), //harmonic point
					    yK(3,1), //projection of center of K onto common face
						yL(3,1), //projection of center of L onto common face
						xK(3,1), //center of cell K
						xL(3,1), //center of cell L
						xKL(3,1), //center of faces common to K and L
						nKL(3,1), //normal pointing out of K to L
						lKs(3,1), //co-normal direction in current cell
						lLs(3,1), //co-normal direction in adjacent cell
						KK(3,3), //permeability of current cell
						KL(3,3); //permeability of adjacent cell
				real aKL; //area of the face
				ElementArray<Face> faces; //array of faces of each cell
				int NF; //number of faces
				Cell cK, cL;
#if defined(USE_OMP)
#pragma omp for
#endif
				for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
				{
					cK = m->CellByLocalID(q);
					cK->Centroid(xK.data()); //retrive cell centroid
					faces = cK->getFaces(); //obtain faces of the cell
					NF = (int)faces.size(); //number of faces;
					KK = rMatrix::FromTensor(cK->RealArrayDF(tag_K).data(),cK->RealArrayDF(tag_K).size()); //get permeability for the cell
					rMatrix W(NF,NF), //big cell-wise gradient matrix,
						    NK(NF,3), //co-normals,
							R(NF,3),  //directions to unknowns,
							D(NF,NF), //orthogonal part of the gradient,
							U(NF,NF); //optimized matrix of orthogonal part.
					for(int k = 0; k < NF; ++k) //loop over faces
					{
						aKL = faces[k].Area(); //get area of the face
						faces[k].Centroid(xKL.data()); //get center of the face
						faces[k].OrientedUnitNormal(cK->self(),nKL.data()); //get unit normal pointing out of K into L
						if( !faces[k].Boundary() )
						{
							KL = rMatrix::FromTensor(cL->RealArrayDF(tag_K).data(),cL->RealArrayDF(tag_K).size()); //get permeability for the cell
						}
					} //end of loop over faces
					W = NK*(NK.Transpose()*R).Invert(true).first*NK.Transpose(); //consistency part
					//Retrive orthogonal to R matrix D
					{ 
						//Symmetric orthogonal matrix
						rMatrix DUD = (rMatrix::Unit(NF) - R*(R.Transpose()*R).Invert(true).first*R.Transpose());
						//perfrom singular value decomposition
						//S should be unity matrix with rank NF-3
						rMatrix DUD_U,DUD_S,DUD_V;
						DUD.SVD(DUD_U,DUD_S,DUD_V);
						//compute the rank
						int rank = NF-3;
						//chop matrix to the full rank
						DUD_V.RemoveColumns(rank,NF);
						//assign the matrix
						D = DUD_V;
					}
					//Optimize U matrix
					{
					}
					//Add stability part
					W += U*D.Transpose();
					real_array store_W = cK->RealArrayDV(tag_W); //access data structure for gradient matrix in mesh
					store_W.resize(NF*NF); //resize the structure
					std::copy(W.data(),W.data()+NF*NF,store_W.data()); //write down the gradient matrix
				} //end of loop over cells
			}
            std::cout << "Construct W matrix: " << Timer() - ttt << std::endl;
            
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
#if defined(USE_OMP)
#pragma omp parallel for
#endif
                for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) ) //loop over cells
                {
                    Cell cell = m->CellByLocalID(q);
                    ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
                    int NF = (int)faces.size();
                    rMatrix nKGRAD(cell->RealArrayDV(tag_W).data(),NF,NF); //Matrix for gradient
                    vMatrix pF(NF,1); //vector of pressure differences on faces
                    vMatrix FLUX(NF,1); //computed flux on faces
                    for(int k = 0; k < NF; ++k)
                        pF(k,0) = (P(faces[k]) - P(cell))*faces[k].Area();
                    FLUX = nKGRAD*pF; //fluxes on faces
                    //Change matrix by nonlinear correction for DMP

                    
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
