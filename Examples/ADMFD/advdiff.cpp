#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


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

bool print_niter = false; //save file on nonlinear iterations


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
                p.SetMethod(Partitioner::INNER_KMEANS,repartition ? Partitioner::Repartition : Partitioner::Partition); // Specify the partitioner
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
		Tag tag_U;  // Normal velocity vector on faces

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

			if( m->HaveTag("VELOCITY") )
			{
				tag_U = m->GetTag("VELOCITY");
				assert(tag_U.isDefined(FACE)); //should be normal velocity on face
				assert(tag_U.GetSize() == 1); //only one component
				//check velocity field is div-free
				if( check_div )
				{
					for(int k = 0; k < m->CellLastLocalID(); ++k) if( m->isValidCell(k) )
					{
						real	U, //face normal velocity
						A, //face area
						sgn; //sign of the normal
						Cell c = m->CellByLocalID(k);
						ElementArray<Face> faces = c.getFaces();
						real divU = 0.0; //divergence of the velocity
						for(int q = 0; q < (int)faces.size(); ++q)
						{
							sgn = (faces[q].FaceOrientedOutside(c) ? 1 : -1); //retrive sign from orientation
							U = faces[q].Real(tag_U); //retrive normal velocity
							A = faces[q].Area(); //retrive area
							divU += U*A*sgn; //compute divergence
						}
						//divU /= c->Volume();
						if( fabs(divU) > 1.0e-8 )
							std::cout << "Velocity at cell " << c->LocalID() << " is not divergence-free: " << divU << std::endl;
					}
				}
			}

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
				rMatrix NK, R, Areas, L;
				rMatrix x(1,3), xf(1,3), n(1,3);
				double area; //area of the face
				double volume; //volume of the cell
				double vel;
#if defined(USE_OMP)
#pragma omp for
#endif
				 for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
				 {
					 Mesh * mesh = m;
					 Cell cell = m->CellByLocalID(q);
					 ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
					 int NF = (int)faces.size(); //number of faces;
					 rMatrix W(NF,NF);
					 volume = cell->Volume(); //volume of the cell
					 cell->Centroid(x.data());
					 //get permeability for the cell
					 rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),
													 cell->RealArrayDF(tag_K).size());
					 NK.Resize(NF,3); //co-normals
					 R.Resize(NF,3); //directions
					 L.Resize(NF,NF);
					 Areas.Resize(NF,NF); //areas
					 Areas.Zero();
					 for(int k = 0; k < NF; ++k) //loop over faces
					 {
						 area = faces[k].Area();
						 faces[k].Centroid(xf.data());
						 faces[k].OrientedUnitNormal(cell->self(),n.data());
						 // assemble matrix of directions
						 R(k,k+1,0,3) = (xf-x)*area;
						 // assemble matrix of co-normals
						 NK(k,k+1,0,3) = area*n*K;
						 L(k,k) = area*n.DotProduct(n*K)/n.DotProduct(xf-x);
						 
						 if( tag_U.isValid() )
						 {
							 vel = faces[k]->Real(tag_U)*(faces[k].FaceOrientedOutside(cell)?1:-1);
							 if( vel > 0 ) NK(k,k+1,0,3) -= area*vel*(xf-x);
						 }
					 } //end of loop over faces
					 W = NK*(NK.Transpose()*R).PseudoInvert(1.0e-12)*NK.Transpose(); //stability part
					 W+= L - (L*R)*((L*R).Transpose()*R).PseudoInvert(1.0e-9)*(L*R).Transpose();
					 //access data structure for gradient matrix in mesh
					 real_array store_W = cell->RealArrayDV(tag_W);
					 //resize the structure
					 store_W.resize(NF*NF);
					 //write down the gradient matrix
					 std::copy(W.data(),W.data()+NF*NF,store_W.data());
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
            Automatizator aut; // declare class to help manage unknowns
            Automatizator::MakeCurrent(&aut);
            dynamic_variable P(aut,aut.RegisterTag(tag_P,CELL|FACE)); //register pressure as primary unknown
            aut.EnumerateEntries(); //enumerate all primary variables
            std::cout << "Enumeration done, size " << aut.GetLastIndex() - aut.GetFirstIndex() << std::endl;

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
#pragma omp parallel
#endif
				{
					vMatrix pF; //vector of pressure differences on faces
					vMatrix FLUX; //computed flux on faces
#if defined(USE_OMP)
#pragma omp for
#endif
					for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) ) //loop over cells
					{
						Cell cell = m->CellByLocalID(q);
						ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
						int NF = (int)faces.size();
						
						raMatrix W = raMatrixMake(cell->RealArrayDV(tag_W).data(),NF,NF); //Matrix for gradient
						
						pF.Resize(NF,1);
						FLUX.Resize(NF,1);
						
						for(int k = 0; k < NF; ++k)
							pF(k,0) = (P(faces[k]) - P(cell));
						FLUX = W*pF; //fluxes on faces
						if( cell.GetStatus() != Element::Ghost )
						{
							for(int k = 0; k < NF; ++k) //loop over faces of current cell
								R[P.Index(cell)] += FLUX(k,0);//faces[k].Area();
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
#pragma omp for
#endif
						 for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
						 {
							 Cell cell = m->CellByLocalID(q);
							 if( cell.GetStatus() == Element::Ghost ) continue;
							 if( cell->HaveData(tag_F) ) R[P.Index(cell)] += cell->Real(tag_F)*cell->Volume();
						 }
					 }
				}
				std::cout << "assembled in " << Timer() - tttt << "\t\t\t" << std::endl;

                std::cout << "Nonlinear residual: " << R.Norm() << "\t\t" << std::endl;

                if( R.Norm() < 1.0e-4 ) break;

				//Solver S(Solver::INNER_ILU2);
                Solver S(Solver::INNER_MPTILUC);
				//Solver S("superlu");
                S.SetParameter("relative_tolerance", "1.0e-14");
                S.SetParameter("absolute_tolerance", "1.0e-12");
                S.SetParameter("drop_tolerance", "1.0e-2");
                S.SetParameter("reuse_tolerance", "1.0e-3");

                S.SetMatrix(R.GetJacobian());
				
                if( S.Solve(R.GetResidual(),Update) )
                {
#if defined(USE_OMP)
#pragma omp parallel for
#endif
                    for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
                    {
                        Cell cell = m->CellByLocalID(q);
						if( cell->GetStatus() == Element::Ghost ) continue;
                        cell->Real(tag_P) -= Update[P.Index(cell)];
                    }
#if defined(USE_OMP)
#pragma omp parallel for
#endif
                    for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
                    {
                        Face face = m->FaceByLocalID(q);
						if( face->GetStatus() == Element::Ghost ) continue;
                        face->Real(tag_P) -= Update[P.Index(face)];
                    }
                    m->ExchangeData(tag_P, CELL|FACE, 0);
					
					if( print_niter )
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
                    std::cout << "Unable to solve: " << S.ReturnReason() << std::endl;
                    break;
                }
                ++nit;
            } while( R.Norm() > 1.0e-4 && nit < 10); //check the residual norm
        }
        std::cout << "Solved problem in " << Timer() - ttt << " seconds with " << nit << " iterations " << std::endl;

        if( m->HaveTag("REFERENCE_SOLUTION") )
        {
            Tag tag_E = m->CreateTag("ERROR",DATA_REAL,CELL,NONE,1);
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
			C = m->AggregateMax(C);
			L2 = m->Integrate(L2);
			volume = m->Integrate(volume);
            L2 = sqrt(L2/volume);
            if( m->GetProcessorRank() == 0 ) std::cout << "Error on cells, C-norm " << C << " L2-norm " << L2 << std::endl;
            C = L2 = volume = 0.0;
            if( tag_R.isDefined(FACE) )
            {
                tag_E = m->CreateTag("ERROR",DATA_REAL,FACE,NONE,1);
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
				C = m->AggregateMax(C);
				L2 = m->Integrate(L2);
				volume = m->Integrate(volume);
                L2 = sqrt(L2/volume);
                if( m->GetProcessorRank() == 0 ) std::cout << "Error on faces, C-norm " << C << " L2-norm " << L2 << std::endl;
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
