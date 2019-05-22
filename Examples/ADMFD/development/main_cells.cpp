//this is not finished
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
typedef Storage::real real;
typedef Storage::integer integer;
typedef Storage::enumerator enumerator;
typedef Storage::real_array real_array;
typedef Storage::var_array var_array;


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
			if( m->GetProcessorRank() == 0 )
				std::cout << "Processors: " << m->GetProcessorsNumber() << std::endl;
			if( m->isParallelFileFormat(argv[1]) ) //The format is
			{
				m->Load(argv[1]); // Load mesh from the parallel file format
				repartition = true; // Ask to repartition the mesh
			}
			else if( m->GetProcessorRank() == 0 )
				m->Load(argv[1]); // Load mesh from the serial file format
			BARRIER
			if( m->GetProcessorRank() == 0 ) std::cout << "Load the mesh: " << Timer()-ttt << std::endl;
		}
		
		
#if defined(USE_PARTITIONER)
		if (m->GetProcessorsNumber() > 1 && !repartition) // Currently only non-distributed meshes are supported by Inner_RCM partitioner
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
		Tag tag_H;  // Harmonic points
		Tag tag_W;  // Gradient matrix
		Tag tag_D;  // Entries for scaling matrix D
		Tag tag_T;  // Transmissibility
		Tag tag_Y;  // Transverse component of tensor
		
		if( m->GetProcessorsNumber() > 1 ) //skip for one processor job
		{ // Exchange ghost cells
			ttt = Timer();
			m->ExchangeGhost(1,FACE); // Produce layer of ghost cells
			m->ExchangeData(tag_H,FACE,0); //Synchronize harmonic points
			BARRIER
			if( m->GetProcessorRank() == 0 ) std::cout << "Exchange ghost: " << Timer()-ttt << std::endl;
		}
		
		{ //initialize data
			if( m->HaveTag("PERM") ) // is diffusion tensor already defined on the mesh? (PERM from permeability)
				tag_K = m->GetTag("PERM"); // get the diffusion tensor
			
			if( !tag_K.isValid() || !tag_K.isDefined(CELL) ) // diffusion tensor was not initialized or was not defined on cells.
			{
				tag_K = m->CreateTag("PERM",DATA_REAL,CELL,NONE,6); // create a new tag for symmetric diffusion tensor K
				for( Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell ) // loop over mesh cells
				{
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
				tag_P = m->CreateTag("PRESSURE",DATA_REAL,CELL|FACE,FACE,1); // Create a new tag for the pressure
				for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell) //Loop over mesh cells
					cell->Real(tag_P) = (rand()*1.0)/(RAND_MAX*1.0); // Prescribe random value in [0,1]
			}
			
			
			if( m->HaveTag("BOUNDARY_CONDITION") ) //Is there boundary condition on the mesh?
			{
				tag_BC = m->GetTag("BOUNDARY_CONDITION");
				
			}
			
			//initialize unknowns at boundary
			for(Mesh::iteratorFace face = m->BeginFace(); face != m->EndFace(); ++face) if( face->Boundary() ) //if( face->HaveData(tag_BC) )
				face->Real(tag_P) = 0.0; //create zero entry
			
			m->ExchangeData(tag_P,CELL|FACE,0); //Synchronize initial solution with boundary unknowns
			
			
			//run in a loop to identify boundary pressures so that they enter as primary variables
			
			//3 entries are coordinates and last entry is the coefficient for FrontCell
			tag_H = m->CreateTag("HARMONIC_POINT_",DATA_REAL,FACE,NONE,4);
			tag_Y = m->CreateTag("GAMMA",DATA_REAL,FACE,NONE,3);
			tag_T = m->CreateTag("TRANSMISSIBILITY",DATA_REAL,FACE,NONE,1);
			// Assemble coefficients for harmonic averaging matrix H
			for(Mesh::iteratorFace fKL = m->BeginFace(); fKL != m->EndFace(); ++fKL) //go over faces
			{
				real_array H = fKL->RealArrayDF(tag_H); //access data structure
				real & T = fKL->RealDF(tag_T);
				real_array Y = fKL->RealArrayDF(tag_Y);
				Cell cK = fKL->BackCell(); //get cell K from the back of the normal
				Cell cL = fKL->FrontCell(); //get cell L to the front of the normal
				
				real dK, dL, lK, lL, D;
				real coefK, coefL, coefQ, coefDiv,aF;
				
				rMatrix y(3,1); //harmonic point
				rMatrix yK(3,1); //projection of the center of cell K onto surface
				rMatrix yL(3,1); //projection of the center of cell L onto surface
				rMatrix xK(3,1); //center of cell K
				rMatrix xL(3,1); //center of cell L
				rMatrix xKL(3,1); //point on face
				rMatrix nKL(3,1); //normal to face
				rMatrix Ks(3,1); //reminder of the co-normal vector Kn in cell K when projection to normal is substracted, i.e. (Kn - n n.Kn)
				rMatrix Ls(3,1); //reminder of the co-normal vector Kn in cell L when projection to normal is substracted, i.e. (Kn - n n.Kn)
				rMatrix gamma(3,1); //transverse part of co-normal not accounted with two-point transmissbility
				rMatrix KK = rMatrix::FromTensor(cK->RealArray(tag_K).data(),cK->RealArray(tag_K).size()).Transpose(); //diffusion tensor in cell K
				
				cK->Centroid(xK.data()); //retrive center of cell K
				fKL->Centroid(xKL.data()); //retrive center of face
				fKL->OrientedUnitNormal(cK,nKL.data()); //retrive unit normal to face
				aF = fKL->Area();
				
				
				D = nKL.DotProduct(xKL); // Compute constant in equation for plane (nx,ny,nz,D)
				dK = D-nKL.DotProduct(xK); //compute distance to the center of cell K

				if( dK < 0 ) std::cout << "incorrect sign for dK " << dK << std::endl;
				
				Ks = KK*nKL; //get conormal in cell K
				lK = nKL.DotProduct(Ks); //find projection of conormal onto normal in cell K
				yK = xK + nKL*dK; //compute projection of the center of cell K onto face

				
				if( !cL.isValid() ) //this is boundary face
				{
					
					T = lK / dK * aF;
					
					//here Ks is full co-normal
					gamma = Ks - (xKL - xK) * lK / dK;
					
					y = xKL;
					
					coefL = 0;
				}
				else
				{
					rMatrix KL = rMatrix::FromTensor(cL->RealArray(tag_K).data(),cL->RealArray(tag_K).size()).Transpose(); //diffusion tensor in cell L
				
					cL->Centroid(xL.data()); //retrive center of cell L
				
					Ls = KL*nKL; //get conormal in cell L
					lL = nKL.DotProduct(Ls); //find projection of conormal onto normal in cell L
					
					Ks -= nKL*lK; //obtain reminder in cell K
					Ls -= nKL*lL; //obtain reminder in cell L
				
					dL = D-nKL.DotProduct(xL); //compute distance to the center of cell L
					
					if( dL > 0 ) std::cout << "incorrect sign for dK " << dK << std::endl;
				
					if( dK*dL > 0 ) //check consistency of geometry
						std::cout << "Cell centers are on the same side from the face" << std::endl;
					
					dL = -dL;
				
					if( dK < 1.0e-9 || dL < 1.0e-9 ) std::cout << "Cell center is located on the face, dK " << dK << " dL " << dL << std::endl;
				
					yL = xL - nKL*dL; //compute projection of the center of cell L onto face
				
					//compute coefficients for harmonic point
					coefDiv = (lL*dK+lK*dL);
					coefL = lL*dK/coefDiv;
					coefK = lK*dL/coefDiv;
					coefQ = dK*dL/coefDiv;

					assert(coefDiv >= 0.0);
					assert(coefL >= 0.0);
					assert(coefK >= 0.0);
					assert(coefQ >= 0.0);
				
					T = lL*lK/coefDiv*aF;
				
					gamma = (yK-yL)*T + coefL*Ks + coefK*Ls;
				
					y = yK*coefK + yL*coefL + (Ks - Ls)*coefQ; //calculate position of harmonic point
				}
				
				std::copy(gamma.data(),gamma.data()+3,Y.data()); //store co-normal reminder
				std::copy(y.data(),y.data()+3,H.data()); //store position of harmonic point
				H[3] = coefL; //store multiplier for pressure in FrontCell, the other is coefK = 1 - coefL
			} //end of loop over faces
			
			tag_W = m->CreateTag("GRAD",DATA_REAL,CELL,NONE);
			//tag_D = m->CreateTag("DIAG",DATA_VARIABLE,CELL,NONE);
			//Assemble all-positive gradient matrix W on cells
			for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell)
			{
				//real_array D = cell->VariableArrayDV(tag_D); //resize scaling matrix D for the future use
				real_array W = cell->RealArrayDV(tag_W); //access data structure for gradient matrix in mesh
				real xP[3]; //center of the cell
				real nF[3]; //normal to the face
				real vP = cell->Volume();
				cell->Centroid(xP);
				ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
				int NF = (int)faces.size(); //number of faces;
				rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),cell->RealArrayDF(tag_K).size()); //get permeability for the cell
				rMatrix GRAD(NF,NF), NK(NF,3), R(NF,3); //big gradient matrix, co-normals, directions
				rMatrix H(NF,NF),T(NF,NF);
				GRAD.Zero();
				H.Zero();
				T.Zero();
				for(int k = 0; k < NF; ++k) //loop over faces
				{
					real aF = faces[k].Area();
					real sign = faces[k].BackCell() == cell->self() ? 1.0 : -1.0;
					real_array yF = faces[k]->RealArrayDF(tag_H); //point on face is harmonic point
					real_array gammaF = faces[k]->RealArrayDF(tag_Y);
					//faces[k]->OrientedUnitNormal(cell->self(),nF);
					//rMatrix nK = rMatrix::FromVector(nF,3).Transpose()*K;
					// assemble matrix of directions
					R(k,0) = (yF[0]-xP[0])*aF;
					R(k,1) = (yF[1]-xP[1])*aF;
					R(k,2) = (yF[2]-xP[2])*aF;
					// assemble matrix of co-normals
					NK(k,0) = gammaF[0]*sign;
					NK(k,1) = gammaF[1]*sign;
					NK(k,2) = gammaF[2]*sign;
					//NK(k,0) = nK(0,0);
					//NK(k,1) = nK(0,1);
					//NK(k,2) = nK(0,2);
					H(k,k) = aF * (0.5*(sign+1) - sign*yF[3]);
					T(k,k) = faces[k]->RealDF(tag_T);
				} //end of loop over faces
				GRAD = NK*(NK.Transpose()*R).Invert().first*NK.Transpose(); //stability part
				GRAD += (rMatrix::Unit(NF) - R*(R.Transpose()*R).Invert().first*R.Transpose())*(2.0/(vP*NF)*GRAD.Trace());
				GRAD = H*GRAD*H + T;

				W.resize(NF*NF); //resize the structure
				std::copy(GRAD.data(),GRAD.data()+NF*NF,W.data()); //write down the gradient matrix
				//std::cout << "Grad matrix:" << std::endl;
				//GRAD.Print();
			} //end of loop over cells
		} //end of initialize data
		
		
		
		
		
		{ //Main loop for problem solution
			Automatizator aut(m); // declare class to help manage unknowns
			Automatizator::MakeCurrent(&aut);
			dynamic_variable P(aut,aut.RegisterDynamicTag(tag_P,CELL|FACE)); //register pressure as primary unknown
			aut.EnumerateDynamicTags(); //enumerate all primary variables
			Residual R("",aut.GetFirstIndex(),aut.GetLastIndex());
			
			
			Sparse::Vector Update("",aut.GetFirstIndex(),aut.GetLastIndex()); //vector for update
			
			int nit = 0;
			
			do
			{
				R.Clear(); //clean up the residual
				
				//First we need to evaluate the gradient at each cell for scaling matrix D
				for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell) //loop over cells
				{
					ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
					int NF = (int)faces.size();
					rMatrix GRAD(cell->RealArrayDV(tag_W).data(),NF,NF); //Matrix for gradient
					vMatrix pF(NF,1); //vector of pressures on faces
					vMatrix FLUX(NF,1); //computed fluxes
					for(int k = 0; k < NF; ++k)
					{
						if( faces[k]->Boundary() )
							pF(k,0) = P(faces[k]) - P(cell->self());
						else
							pF(k,0) = P(faces[k]->FrontCell()) - P(cell->self());
					}
					FLUX = GRAD*pF; //fluxes on faces
					
					for(int k = 0; k < NF; ++k)
					{
						R[P.Index(cell->self())] += FLUX(k,0);
						if( faces[k]->Boundary() )
						{
							if( faces[k]->HaveData(tag_BC) )
							{
								real_array bc = faces[k]->RealArray(tag_BC);
								R[P.Index(faces[k])] -= bc[0]*P(faces[k]) + bc[1]*FLUX(k,0) - bc[2];
							}
							else R[P.Index(faces[k])] -= FLUX(k,0);
						}
						else
							R[P.Index(faces[k]->FrontCell())] -= FLUX(k,0);
					}
				} //end of loop over cells
				

				
				if( tag_F.isValid() )
				{
					for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell)
						R[P.Index(cell->self())] += cell->RealDF(tag_F)*cell->Volume();
				}
				
				R.Rescale();
				
				
				
				std::cout << "Nonlinear residual: " << R.Norm() << std::endl;
				
				if( R.Norm() < 1.0e-4 ) break;
				
				Solver S(Solver::INNER_MPTILUC);
				S.SetParameterReal("drop_tolerance",1.0e-2);
				S.SetParameterReal("reuse_tolerance",1.0e-4);
				S.SetParameterEnum("gmres_substeps",4);
				
				R.GetJacobian().Save("jacobian.mtx");
				R.GetResidual().Save("residual.mtx");
				
				S.SetMatrix(R.GetJacobian());
				if( S.Solve(R.GetResidual(),Update) )
				{
					for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell)
						cell->Real(tag_P) -= Update[P.Index(cell->self())];
					for(Mesh::iteratorFace face = m->BeginFace(); face != m->EndFace(); ++face) if( face->HaveData(tag_P) )
						face->Real(tag_P) -= Update[P.Index(face->self())];
				}
				else
				{
					std::cout << "Unable to solve: " << S.GetReason() << std::endl;
					break;
				}
				
			} while( R.Norm() > 1.0e-4 && nit < 10 ); //check the residual norm
			
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
			}
			
			m->Save("out.pmf");
			
			if( m->GetProcessorsNumber() == 1 )
				m->Save("out.vtk");
			else
				m->Save("out.pvtk");

		}
		
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
