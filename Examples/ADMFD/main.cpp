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

const bool hybrid = false;
const double s = 1;

template<typename M> void PrintMatrix(const M & mat)
{
	for(int i = 0; i < mat.Rows(); ++i)
	{
		for(int j = 0; j < mat.Cols(); ++j)
			std::cout << get_value(mat(i,j)) << " ";
		std::cout << std::endl;
	}
}

// A = B
void CopyMatrix(Sparse::Matrix & A, const Sparse::Matrix & B)
{
	A.SetInterval(B.GetFirstIndex(),B.GetLastIndex());
	for(int k = B.GetFirstIndex(); k < B.GetLastIndex(); ++k)
	{
		A[k].Clear();
		for(int l = 0; l < B[k].Size(); ++l)
			A[k].Push(B[k].GetIndex(l),B[k].GetValue(l));
	}
}

// A = B
void CopyVector(Sparse::Vector & A, const Sparse::Vector & B)
{
	A.SetInterval(B.GetFirstIndex(),B.GetLastIndex());
	for(int k = B.GetFirstIndex(); k < B.GetLastIndex(); ++k)
		A[k] = B[k];
}


//A = B+C
void AddMatrix(Sparse::Matrix & A, const Sparse::Matrix & B, const Sparse::Matrix & C)
{
	assert(B.GetFirstIndex() == C.GetFirstIndex());
	assert(B.GetLastIndex() == C.GetLastIndex());
	A.SetInterval(B.GetFirstIndex(),B.GetLastIndex());
	for(int k = B.GetFirstIndex(); k < B.GetLastIndex(); ++k)
	{
		A[k].Clear();
		for(int l = 0; l < B[k].Size(); ++l)
			A[k][B[k].GetIndex(l)] += B[k].GetValue(l);
		for(int l = 0; l < C[k].Size(); ++l)
			A[k][C[k].GetIndex(l)] += C[k].GetValue(l);
	}
}

void MultMatrix(Sparse::Matrix & A, double c)
{
	for(int k = A.GetFirstIndex(); k < A.GetLastIndex(); ++k)
	{
		for(int l = 0; l < A[k].Size(); ++l)
			A[k].GetValue(l) *= c;
	}
}

void MultVector(Sparse::Vector & A, double c)
{
	for(int k = A.GetFirstIndex(); k < A.GetLastIndex(); ++k)
	{
		A[k] *= c;
	}
}


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
		Tag tag_W;  // Local approximation matrix
		Tag tag_L;  // store two-point half fluxes
		Tag tag_Q;  // store flux limiting values
		Tag tag_H;  // harmonic points
		
		
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
			/*
			if( m->HaveTag("REFERENCE_SOLUTION") )
			{
				Tag tag_Pr = m->GetTag("REFERENCE_SOLUTION");
				for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE); it != m->EndElement(); ++it)
					it->Real(tag_P) = it->Real(tag_Pr);
			}
			 */
			
			if( hybrid && !tag_P.isDefined(FACE) )
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
			m->ExchangeData(tag_P,CELL|(hybrid?FACE:NONE),0); //Synchronize initial solution with boundary unknowns
			tag_W = m->CreateTag("W",DATA_REAL,CELL,NONE);
			tag_L = m->CreateTag("L",DATA_REAL,CELL,NONE);
			if( hybrid ) tag_H = m->CreateTag("H",DATA_REAL,FACE,NONE,3);
			tag_Q = m->CreateTag("Q",DATA_VARIABLE,FACE,NONE,2);
			
			ttt = Timer();
			//Assemble gradient matrix W on cells
			int total = 0, dmp = 0;
#if defined(USE_OMP)
#pragma omp parallel reduction(+:total) reduction(+:dmp)
#endif
			{
				rMatrix xc(3,1), //center of the cell
						xn(3,1), //center of neighbour cell
						ys(3,1), //harmonic point on face
						nf(3,1), //normal to the face
						Knf(3,1), //co-normal vector
						r(3,1);  //vector from cell center to face center
				rMatrix W, //approximation matrix
						K(3,3), Kn(3,3), //permeability tensor
						N, //matrix of normals
						R, //matrix of vectors
						V,L; //half-flux transmissibility
				rMatrix x1(3,1), x2(3,1), K1(3,3), K2(3,3);
				real l1,l2,d1,d2;
				real l,T, //half-flux two-point transmissibility
					 d, //distance from cell to face
					 af, //area of the face
					mult_T,mult_xi, mult_r, mult_p;
				int NF; //number of faces
				//find out harmonic points
				if( hybrid )
				{
#if defined(USE_OMP)
#pragma omp for
#endif
					for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
					{
						Face face = m->FaceByLocalID(q);
						real_array H = face->RealArray(tag_H);
						Cell c1 = face->BackCell();
						Cell c2 = face->FrontCell();
						face->Centroid(ys.data());
						if( c2.isValid() )
						{
							face->UnitNormal(nf.data());
							c1.Centroid(x1.data());
							K1 = rMatrix::FromTensor(c1->RealArrayDF(tag_K).data(),c1->RealArrayDF(tag_K).size());
							d1 = nf.DotProduct(ys-x1);
							l1 = nf.DotProduct(K1*nf);
							c2.Centroid(x2.data());
							K2 = rMatrix::FromTensor(c2->RealArrayDF(tag_K).data(),c2->RealArrayDF(tag_K).size());
							d2 = nf.DotProduct(x2-ys);
							l2 = nf.DotProduct(K2*nf);
							r = (l2*d1*x2+l1*d2*x1+d1*d2*(K1-K2)*nf)/(l1*d2+l2*d1);
							H[0] = r(0,0);
							H[1] = r(1,0);
							H[2] = r(2,0);
						}
						else
						{
							H[0] = ys(0,0);
							H[1] = ys(1,0);
							H[2] = ys(2,0);
						}
					}
				}
				//find out approximation matrix
#if defined(USE_OMP)
#pragma omp for
#endif
				for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
				{
					Cell cell = m->CellByLocalID(q);
					real_array store_W = cell->RealArray(tag_W);
					real_array store_L = cell->RealArray(tag_L);
					ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
					NF = (int)faces.size(); //number of faces;
					cell->Centroid(xc.data());
					K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),cell->RealArrayDF(tag_K).size());
					W.Resize(NF,NF);
					N.Resize(NF,3);
					R.Resize(NF,3);
					V.Resize(NF,3);
					L.Resize(NF,NF);
					store_L.resize(NF);
					for(int k = 0; k < NF; ++k) //loop over faces
					{
						mult_T = mult_xi = mult_r = 1;
						mult_p = 0;
						faces[k].OrientedUnitNormal(cell,nf.data());
						if( hybrid )
						{
							ys = rMatrix::FromVector(faces[k].RealArray(tag_H).data(),3);
							
							r = ys-xc;
						}
						else
						{
							Cell n = cell.Neighbour(faces[k]);
							faces[k].Centroid(ys.data());
							if( n.isValid() )
							{
								Kn = rMatrix::FromTensor(n->RealArrayDF(tag_K).data(),n->RealArrayDF(tag_K).size());
								n.Centroid(xn.data());
								d = nf.DotProduct(xn-ys);
								l = nf.DotProduct(Kn*nf);
								r = xn-xc + d/l*(K-Kn)*nf;
								// d = n.r = d1+d2 + d2/l2*(l1-l2) = d1 + d2 + d2*l1/l2 - d2 = (d1*l2 + d2*l1)/l2
								// l1/d = l1*l2/(d1*l2+d2*l1)
							}
							else
							{
								real BC[3] = {0,1,0};
								if( tag_BC.isValid() && faces[k].HaveData(tag_BC) )
								{
									BC[0] = faces[k].RealArray(tag_BC)[0];
									BC[1] = faces[k].RealArray(tag_BC)[1];
									BC[2] = faces[k].RealArray(tag_BC)[2];
								}
								d = nf.DotProduct(ys-xc);
								l = nf.DotProduct(K*nf);
								T = l/d;
								
								mult_r = (BC[0] + s*BC[1]*T);
								mult_T = 1.0/mult_r;
								mult_xi = BC[0]*mult_T;
								mult_p = BC[1];
								
								r = ys - xc;
								
								
							}
						}
						
						Knf = K*nf;
						l = fabs(nf.DotProduct(Knf));
						d = fabs(nf.DotProduct(r));
						af = faces[k].Area();
						
						T = l/d;
						
						// assemble matrix of directions
						R(k,0) = r(0,0)*mult_r + mult_p*(Knf(0,0) - s*T*r(0,0));
						R(k,1) = r(1,0)*mult_r + mult_p*(Knf(1,0) - s*T*r(1,0));
						R(k,2) = r(2,0)*mult_r + mult_p*(Knf(2,0) - s*T*r(2,0));
						// assemble matrix of normals
						N(k,0) = nf(0,0)*af;
						N(k,1) = nf(1,0)*af;
						N(k,2) = nf(2,0)*af;
						
						
						V(k,0) = (Knf(0,0) - s*T*r(0,0))*af*mult_xi;
						V(k,1) = (Knf(1,0) - s*T*r(1,0))*af*mult_xi;
						V(k,2) = (Knf(2,0) - s*T*r(2,0))*af*mult_xi;
						//L(k,k) = l/d*af;
						store_L[k] = mult_T*s*T*af;
						//L(k,k) = mult_T*T*af;
					} //end of loop over faces
					
					W = V*(N.Transpose()*R).Invert(true).first*N.Transpose();
					if( W.CheckNans() )
					{
						std::cout << "Nans on " << cell->LocalID() << std::endl;
						std::cout << "W" << std::endl;
						PrintMatrix(W);
						std::cout << "V" << std::endl;
						PrintMatrix(V);
						std::cout << "R" << std::endl;
						PrintMatrix(R);
					}
					//W = N*K*((N*K).Transpose()*R).Invert(true).first*(N*K).Transpose() - s*L*R*((L*R).Transpose()*R).Invert(true).first*(L*R).Transpose();
					//W = (N*K-s*L*R)*(N.Transpose()*R).Invert(true).first*N.Transpose();
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
			Automatizator aut; // declare class to help manage unknowns
			Automatizator::MakeCurrent(&aut);
			dynamic_variable P(aut,aut.RegisterTag(tag_P,CELL|(hybrid?FACE:NONE))); //register pressure as primary unknown
			aut.EnumerateEntries(); //enumerate all primary variables
			
			std::cout << "Enumeration done, size " << aut.GetLastIndex() - aut.GetFirstIndex() << std::endl;
			
			
			Residual Resid("",aut.GetFirstIndex(),aut.GetLastIndex());
			Sparse::Matrix A("",aut.GetFirstIndex(),aut.GetLastIndex());
			Sparse::Vector x("",aut.GetFirstIndex(),aut.GetLastIndex());
			Sparse::Matrix J_prev("",aut.GetFirstIndex(),aut.GetLastIndex());
			Sparse::Vector x_prev("",aut.GetFirstIndex(),aut.GetLastIndex());
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
				if(hybrid)for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
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
				Resid.Clear(); //clean up the residual
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
					rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),cell->RealArrayDF(tag_K).size());
					rMatrix L = rMatrix::FromDiagonal(cell->RealArrayDV(tag_L).data(),NF);
					rMatrix W(cell->RealArray(tag_W).data(),NF,NF);
					//rMatrix R(cell->RealArrayDV(tag_R).data(),NF,3); //Matrix for directions
					//rMatrix N(cell->RealArrayDV(tag_N).data(),NF,3); //Matrix for directions
					//rMatrix G(cell->RealArrayDV(tag_G).data(),3,NF); //Matrix for gradient
					//rMatrix V(cell->RealArrayDV(tag_V).data(),NF,3); //Matrix for transversal directions
					vMatrix DP(NF,1); //vector of pressure differences on faces
					vMatrix FLUX(NF,1);
					
					if( hybrid )
					{
						for(int k = 0; k < NF; ++k)
							DP(k,0) = (P[faces[k]] - P[cell]);
					}
					else
					{
						for(int k = 0; k < NF; ++k)
						{
							Cell n = cell.Neighbour(faces[k]);
							if( n.isValid() )
								DP(k,0) = (P[n] - P[cell]);
							else if( faces[k].HaveData(tag_BC) )
							{
								real_array BC = faces[k].RealArray(tag_BC);
								DP(k,0) = BC[2] - BC[0]*P[cell];
							}
							else DP(k,0) = 0;
						}
					}
					
					FLUX = W*DP; //fluxes on faces
					
					for(int k = 0; k < NF; ++k) //loop over faces of current cell
					{
						var_array Q = faces[k]->VariableArray(tag_Q);
						int ind = faces[k].BackCell() == cell ? 0 : 1;
						Q[ind] = FLUX(k,0);
					}
				}
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:total) reduction(+:dmp)
#endif
				for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) ) //loop over cells
				{
					Cell cell = m->CellByLocalID(q);
					ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
					int NF = (int)faces.size();
					//rMatrix G(cell->RealArrayDV(tag_G).data(),3,NF); //Matrix for gradient
					//rMatrix V(cell->RealArrayDV(tag_V).data(),NF,3); //Matrix for transversal directions
					rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),cell->RealArrayDF(tag_K).size());
					rMatrix W(cell->RealArray(tag_W).data(),NF,NF);
					
					vMatrix J(NF,NF);
					//vMatrix M(NF,NF);
					//rMatrix W(NF,NF);
					//rMatrix R(cell->RealArrayDV(tag_R).data(),NF,3); //Matrix for directions
					//rMatrix N(cell->RealArrayDV(tag_N).data(),NF,3); //Matrix for directions
					rMatrix L = rMatrix::FromDiagonal(cell->RealArrayDV(tag_L).data(),NF);
					vMatrix DP(NF,1); //vector of pressure differences on faces
					vMatrix FLUX(NF,1); //computed flux on faces
					//vMatrix iM1(3,3), U(3,3), S(3,3), V(3,3);
					
					//rMatrix W(cell->RealArrayDV(tag_W).data(),NF,NF); //Matrix for gradient
					if( hybrid )
					{
						for(int k = 0; k < NF; ++k)
							DP(k,0) = (P[faces[k]] - P[cell]);
					}
					else
					{
						for(int k = 0; k < NF; ++k)
						{
							Cell n = cell.Neighbour(faces[k]);
							if( n.isValid() )
								DP(k,0) = (P[n] - P[cell]);
							else if( faces[k].HaveData(tag_BC) )
							{
								real_array BC = faces[k].RealArray(tag_BC);
								DP(k,0) = BC[2] - BC[0]*P[cell];
							}
							else DP(k,0) = 0;
						}
					}
					
					//std::cout << "cell " << cell->LocalID() << std::endl;
					for(int k = 0; k < NF; ++k) if( !faces[k].Boundary() )//loop over faces of current cell
					{
						var_array Q = faces[k]->VariableArray(tag_Q);
						int ind1,ind2;
						ind1 = faces[k].BackCell() == cell ? 0 : 1;
						ind2 = 1-ind1;
						/*
						if( Q[ind1]*Q[ind2] <= 0 )
							J(k,k) = soft_fabs(Q[ind2],1.0e-8)/(soft_fabs(Q[ind2],1.0e-8)+soft_fabs(Q[ind1],1.0e-8))*2;//(1-soft_sign(Q[ind1]*Q[ind2],1.0e-8)) + 1.0e-9;
						else
							J(k,k) = 0.0;
						*/
						//J(k,k) = variation(soft_fabs(Q[ind2],1.0e-8)/(soft_fabs(Q[ind2],1.0e-8)+soft_fabs(Q[ind1],1.0e-8))*(1-soft_sign(Q[ind1]*Q[ind2],1.0e-8)),0.5) + 1.0e-9;
						
						J(k,k) = variation((Q[ind2]*Q[ind2] - Q[ind1]*Q[ind1] - 2*Q[ind1]*Q[ind2] + 1.0e-20)/(2*(Q[ind2]*Q[ind2]+Q[ind1]*Q[ind1])+2.0e-20),0.75) + 0.5 ;
						//
						
						//J(k,k) = (-2*soft_fabs(Q[ind1]*Q[ind2],1.0e-8) + Q[ind1]*Q[ind2] + 2*Q[ind2]*soft_fabs(Q[ind2],1.0e-8)*soft_sign(Q[ind1],1.0e-8) - Q[ind1]*Q[ind1])/(Q[ind1]*Q[ind1]+4*soft_fabs(Q[ind1]*Q[ind2],1.0e-8)+Q[ind2]*Q[ind2]+1.0e-20) + 1;
						
						
						//limiter *= 0.5*J(k,k);
						//if( Q[ind1]*Q[ind2] < 0.0 )
						//	J(k,k) = variation(2*(fabs(Q[ind2])+1.0e-36)/(fabs(Q[ind1])+fabs(Q[ind2])+2.0e-36),0.25);
						//else J(k,k) = 0;
						//J(k,k) = (Q[ind2]*Q[ind2] - Q[ind1]*Q[ind1] - 2*Q[ind1]*Q[ind2])/(2*(Q[ind2]*Q[ind2]+Q[ind1]*Q[ind1])+1.0e-20) + 0.5;
						//std::cout << k << ": " << get_value(Q[ind1]) << " " << get_value(Q[ind2]) << std::endl;
						//std::cout << "Q1 " << get_value(Q[ind1]) << " Q2 " << get_value(Q[ind2]) << " " << get_value(J(k,k)) << std::endl;
					}
					else
					{
						//std::cout << k << ": bnd" << std::endl;
						J(k,k) = 1;
					}
					
					FLUX = (L + J*W)*DP;
					
					
					if( cell.GetStatus() != Element::Ghost )
					{
						for(int k = 0; k < NF; ++k) //loop over faces of current cell
							Resid[P.Index(cell)] += FLUX(k,0);
					}
					if( hybrid ) for(int k = 0; k < NF; ++k) //loop over faces of current cell
					{
						if( faces[k].GetStatus() == Element::Ghost ) continue;
						int index = P.Index(faces[k]);
						Locks.Lock(index);
						if( tag_BC.isValid() && faces[k].HaveData(tag_BC) )
						{
							real_array BC = faces[k].RealArray(tag_BC);
							Resid[index] -= BC[0]*P[faces[k]] + BC[1]*FLUX(k,0) - BC[2];
						}
						else
							Resid[index] -= FLUX(k,0);
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
						if( cell->HaveData(tag_F) ) Resid[P.Index(cell)] += cell->Real(tag_F)*cell->Volume();
					}
				}
				
				std::cout << "assembled in " << Timer() - tttt << "\t\t\t" << std::endl;
				
				
				Resid.Rescale();
				//Resid.GetJacobian().Save("jacobian.mtx",&Text);
				//Resid.GetResidual().Save("residual.mtx");
				
				
				std::cout << "Nonlinear residual: " << Resid.Norm() << "\t\t" << std::endl;
				
				if( Resid.Norm() < 1.0e-5 ) break;
				
				tttt = Timer();
				//Solver S(Solver::INNER_ILU2);
				Solver S(Solver::INNER_MPTILUC);
				//Solver S(Solver::SUPERLU);
				S.SetParameterReal("relative_tolerance", 1.0e-12);
				S.SetParameterReal("absolute_tolerance", 1.0e-9);
				S.SetParameterReal("drop_tolerance", 5.0e-2);
				S.SetParameterReal("reuse_tolerance", 1.0e-3);
				
				
				
				
				S.SetMatrix(Resid.GetJacobian());
				//std::fill(Update.Begin(),Update.End(),0.0);
				if( S.Solve(Resid.GetResidual(),Update) )
				{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
					for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
					{
						Cell cell = m->CellByLocalID(q);
						cell->Real(tag_P) -= Update[P.Index(cell)];
					}
					if( hybrid )
					{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
						{
							Face face = m->FaceByLocalID(q);
							face->Real(tag_P) -= Update[P.Index(face)];
						}
					}
					m->ExchangeData(tag_P, CELL|(hybrid?FACE:NONE), 0);
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
				std::cout << "solved in " << Timer() - tttt << "\t\t\t" << std::endl;
				++nit;
			} while( Resid.Norm() > 1.0e-5 && nit < 100); //check the residual norm
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
			L2 = sqrt(L2/volume);
			std::cout << "Error on cells, C-norm " << C << " L2-norm " << L2 << std::endl;
			C = L2 = volume = 0.0;
			if( false )
			{
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
					L2 = sqrt(L2/volume);
					std::cout << "Error on faces, C-norm " << C << " L2-norm " << L2 << std::endl;
				}
				else std::cout << "Reference solution was not defined on faces" << std::endl;
			}
		}
		
		if( m->GetProcessorsNumber() == 1 )
			m->Save("out.vtk");
		else
			m->Save("out.pvtk");
		m->Save("out.xml");
		
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
