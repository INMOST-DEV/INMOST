#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "checks.h"
#include "save_mesh.h"
#include "conv_diff.h"
using namespace INMOST;

#ifndef M_PI
#define M_PI 3.141592653589
#endif

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

//Solver::Type linear_solver_type = "inner_ilu2";
Solver::Type linear_solver_type = "inner_mptiluc";
//Solver::Type linear_solver_type = Solver::SUPERLU;

//shortcuts
typedef Storage::bulk            bulk;
typedef Storage::real            real;
typedef Storage::integer         integer;
typedef Storage::enumerator      enumerator;
typedef Storage::real_array      real_array;
typedef Storage::var_array       var_array;
typedef Storage::reference_array ref_array;


const bool use_reference_solution = false;

SchemeType scheme_type = MPFA;
//SchemeType scheme_type = NMPFA_QUADRATIC_CORRECTED;
//SchemeType scheme_type = NMPFA_VAN_ALBADA_CORRECTED;
//SchemeType scheme_type = NMPFA_VAN_ALBADA_CORRECTED;
//SchemeType scheme_type = NMPFA_VAN_ALBADA;
//SchemeType scheme_type = NMPFA_VAN_LEER;
//SchemeType scheme_type = NTPFA_PICARD;
//SchemeType scheme_type = NTPFA;
OutputMeshFormat output_format = PMF; //current choice of mesh format

//scheme behavior
bool split_diffusion = true; //diffusion is represented by two-point part plus correction
bool joint_advection_diffusion = true; //nonlinear weighting is performed for advection and diffusion simultaneously
//time stepping constants
//#define FRACTIONAL //tvd-limited fractional time-stepping scheme, bacwards euler if commented
bool bdf2_stepping = false; //Activate limited time derivative corrector that potentially yields bdf2 scheme

//checks for M-matrix property (LED property)
bool check_fluxes = false; //check individual fluxes
bool check_matrix = false; //check the jacobian (will not check coefficient properties for right hand side)
bool check_div = false; //check that the velocity is divergence free.

//parameters for nonlinear solver
double nonlinear_iterations                = 100; //maximal number of iterations
double nonlinear_abs_tolerance             = 1.0e-6; //absolute tolerance
double nonlinear_rel_tolerance             = 1.0e-5; //relative tolerance
double regularization                      = 1.0e-32; //eps regularization in |x| = sqrt(x*x+eps)
double degenerate_diffusion_regularization = 1.0e-32; //eps regularization when diffusion goes to zero

//constants for stencils
int max_layers = 5; //maximum number of layers to be considered in interpolation
ElementType bridge_layers = FACE; //used as bridge elements for layers


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
		(void) repartition;
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
			{ // Compute mesh partitioning.
				ttt = Timer();
				Partitioner p(m); //Create Partitioning object
				//p.SetMethod(Partitioner::Inner_RCM,repartition ? Partitioner::Repartition : Partitioner::Partition); // Specify the partitioner
				p.SetMethod(Partitioner::INNER_KMEANS,Partitioner::Partition); // Specify the partitioner
				//p.SetMethod(Partitioner::Zoltan_PHG,repartition ? Partitioner::Repartition : Partitioner::Partition); // Specify the partitioner
				p.Evaluate(); // Compute the partitioner and store new processor ID in the mesh
				BARRIER
				if( m->GetProcessorRank() == 0 ) std::cout << "Evaluate: " << Timer()-ttt << std::endl;
			}
			
			{ // Distribute the mesh.
				ttt = Timer();
				m->Redistribute(); // Redistribute the mesh data
				m->ReorderEmpty(CELL|FACE|EDGE|NODE); // Clean the data after reordring
				BARRIER
				if( m->GetProcessorRank() == 0 ) std::cout << "Redistribute: " << Timer()-ttt << std::endl;
			}
		}
#endif
		
		{ // Prepare geometrical data on the mesh.
			ttt = Timer();
			Mesh::GeomParam table;
			table[CENTROID]    = CELL | FACE; //Compute averaged center of mass
			table[NORMAL]      = FACE;        //Compute normals
			table[ORIENTATION] = FACE;        //Check and fix normal orientation
			table[MEASURE]     = CELL | FACE; //Compute volumes and areas
			table[BARYCENTER]  = CELL | FACE; //Compute volumetric center of mass
			m->RemoveGeometricData(table);
			m->PrepareGeometricData(table); //Ask to precompute the data
			BARRIER
			if( m->GetProcessorRank() == 0 ) std::cout << "Prepare geometric data: " << Timer()-ttt << std::endl;
		}
		
		// data tags for
		Tag tag_C;  // Concentration
		Tag tag_C_Update; // Concentration update
		Tag tag_C0; // Concentration on previous time step for non-steady-state problem
		Tag tag_C1; // for BDF2
		
		Tag tag_U;  // Normal velocity vector on faces
		Tag tag_K;  // Diffusion tensor
		Tag tag_F;  // Forcing term
		Tag tag_R;  // Reaction term of the form R*C
		Tag tag_BC; // Boundary conditions
		
		Tag tag_FLUX;       // Store computed flux
		Tag tag_PORO;       // Porocity of rock for time derivative

		
		MarkerType boundary_face = m->CreateMarker(); //defines boundary faces in parallel

		m->MarkBoundaryFaces(boundary_face);
		
		
		bool perform_correction_convection = true; //add nonlinear correction to spu
		bool perform_correction_diffusion = true; //add nonlinear correction to spu
		bool steady_state = true; //is there time derivative?
		real TimeStep = 1, TotalTime = 1, CurrentTime = 0;
		
		tag_FLUX = m->CreateTag("FLUX",DATA_REAL,FACE,NONE,1);
		
		if( argc > 2 && atoi(argv[2]) == 0 )
		{
			 if( m->GetProcessorRank() == 0 )
				std::cout << "Switched to single point upwind scheme" << std::endl;
			perform_correction_convection = false;
		}
		
		if( argc > 4 || m->HaveTag("TIME_INFO") ) //command line parameters or mesh information
		{
			if( argc > 4 )
			{
				TimeStep = atof(argv[3]);
				TotalTime = atof(argv[4]);
			}
			else if( m->HaveTag("TIME_INFO") )
			{
				real_array time = m->self()->RealArray(m->GetTag("TIME_INFO")); //first dt then T
				TimeStep = time[0];
				TotalTime = time[1];
			}
			steady_state = false;
			 if( m->GetProcessorRank() == 0 )
				std::cout << "Solving problem in time with time step: " << TimeStep << " total time: " << TotalTime << "." << std::endl;
		}
		else if( m->GetProcessorRank() == 0 ) std::cout << "Solving steady-state problem." << std::endl;
		
		
		if( m->GetProcessorsNumber() > 1 ) //skip for one processor job
		{ // Exchange ghost cells
			ttt = Timer();
			m->ExchangeGhost(2,NODE); // Produce layer of ghost cells
			BARRIER
			if( m->GetProcessorRank() == 0 ) std::cout << "Exchange ghost: " << Timer()-ttt << std::endl;
		}
		
		{ //initialize data
			if( m->HaveTag("PORO") ) //is porosity defined on the mesh?
				tag_PORO = m->GetTag("PORO"); //get the porosity

			if( m->HaveTag("PERM") ) // is diffusion tensor already defined on the mesh? (PERM from permeability)
				tag_K = m->GetTag("PERM"); // get the diffusion tensor
			
			if( m->HaveTag("CONCENTRATION") ) //Is there a concentration on the mesh?
				tag_C = m->GetTag("CONCENTRATION"); //Get the concentration
			else //Create the concentration
				tag_C = m->CreateTag("CONCENTRATION",DATA_REAL,CELL,NONE,1);

			tag_C_Update= m->CreateTag("CONCENTRATION_UPDATE",DATA_REAL,CELL,NONE,1);
			
			if( !steady_state ) //storage for previous time step
			{
				tag_C0 = m->CreateTag("CONCENTRATION0",DATA_REAL,CELL,NONE,1);
				if( bdf2_stepping ) tag_C1 = m->CreateTag("CONCENTRATION1",DATA_REAL,CELL,NONE,1); //for BDF2
			}
			
			if( m->HaveTag("BOUNDARY_CONDITION") ) //Is there boundary condition on the mesh?
			{
				tag_BC = m->GetTag("BOUNDARY_CONDITION");
				assert(tag_BC.isDefined(FACE)); //defined on faces
				assert(tag_BC.GetSize() == 3); //should be 3 entries, bc[0] C + bc[1] dC/dn = bc[2]
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
			
			
			if( m->HaveTag("FORCE") ) //Is there force on the mesh?
			{
				tag_F = m->GetTag("FORCE"); //retrive force
				assert(tag_F.isDefined(CELL)); //assuming it was defined on cells
				assert(tag_F.GetSize() == 1); //only one component
			} // end of force
			
			if( m->HaveTag("REACTION") )
			{
				tag_R = m->GetTag("REACTION"); //reaction
				assert(tag_R.isDefined(CELL)); //should be on cell
				assert(tag_R.GetSize() == 1); //only one component
			}
		} //end of initialize data

		ttt = Timer();
		
		if( m->GetProcessorRank() == 0 )
			std::cout << "Create problem" << std::endl;
		ConvectionDiffusion * Problem = new ConvectionDiffusion(m,tag_U,tag_K,tag_BC,boundary_face,perform_correction_convection,perform_correction_diffusion);
		if( m->GetProcessorRank() == 0 )
			std::cout << "Done create problem" << std::endl;

		if( m->GetProcessorRank() == 0 )
			std::cout << "Precompute gradients: " << Timer() - ttt << std::endl;

		if( Problem->Failure() ) 
		{
			std::string filename = SaveMesh(m,"failed",output_format);
			if( m->GetProcessorRank() == 0 )
				std::cout << "There was error during initialization, check " << filename << std::endl;
			exit(-1);
		}
		
		if( m->GetProcessorRank() == 0 )
			std::cout << "Initialization done" << std::endl;
		
		
		integer	nit = 0, //number of nonlinear iterations performed
				step = 0;  //number of time steps performed
		real	norm, //norm of the residual of nonlinear problem
				norm0 //norm of the residual at first iteration
				;
		real Cmax, Cmin, Cint;

		// test problem with reference solution as initial guess
		
		if( m->HaveTag("REFERENCE_SOLUTION") && use_reference_solution )
		{
			Tag r = m->GetTag("REFERENCE_SOLUTION");
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				it->Real(tag_C) = it->Real(r);
		}
		
		
		{ //Main loop for problem solution
			Automatizator aut; // declare class to help manage unknowns
			Automatizator::MakeCurrent(&aut);
			dynamic_variable C(aut,aut.RegisterTag(tag_C,CELL)); //register concentration as primary unknown
			static_variable C0(tag_C0);
			static_variable C1(tag_C1); //for BDF2
			aut.EnumerateEntries(); //enumerate all primary variables
			
			if( m->GetProcessorRank() == 0 )
				std::cout << "Enumeration done" << std::endl;
			
			
			Tag slope_limiter;
			Tag tag_residual;//, tag_residual_no_force;
			tag_residual = m->CreateTag("RESIDUAL",DATA_REAL,CELL,NONE,1);
			//tag_residual_no_force = m->CreateTag("RESIDUAL_NO_FORCE",DATA_REAL,CELL,NONE,1);
			if( !steady_state && bdf2_stepping )
				slope_limiter = m->CreateTag("BDF2_SLOPE_LIMITER",DATA_REAL,CELL,NONE,1);
			
			
			Residual R("",aut.GetFirstIndex(),aut.GetLastIndex());
			R.InitLocks();
			Sparse::Vector Update  ("",aut.GetFirstIndex(),aut.GetLastIndex()); //vector for update
			
			if( !steady_state )
				SaveMesh(m,"step",0,output_format);
			
			while(CurrentTime < TotalTime)
			{
				
				ttt = Timer();
				
				if( !steady_state )
				{
					//remember concentration
#if defined(USE_OMP)
#pragma omp parallel for
#endif
					for(int q = 0; q < m->CellLastLocalID(); ++q) if( m->isValidCell(q) )
					{
						Cell c = m->CellByLocalID(q);
						if( bdf2_stepping )
							c->Real(tag_C1) = c->Real(tag_C0); //for BDF2
						c->Real(tag_C0) = c->Real(tag_C);
					}
				}
				nit = 0;
				do
				{
					R.Clear(); //clean up the residual
					//double tttt = Timer();
#if defined(USE_OMP)
#pragma omp parallel
#endif
					{
						//shared variables
						variable fluxK, fluxL, fluxT; //two fluxes for Picard's method and single point upwind flux
						variable corrK, corrL;
						variable corrK_cK, corrL_cL;
						Face fKL;
						INMOST_DATA_ENUM_TYPE iK,iL;
						Cell cK, cL;
						//reset flux
#if defined(USE_OMP)
#pragma omp for
#endif
						for(integer q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
							m->FaceByLocalID(q)->Real(tag_FLUX) = 0.0;
#if defined(USE_OMP)
#pragma omp for
#endif
						for(integer q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) ) //loop over faces
						{
							fKL = m->FaceByLocalID(q);
							
							if( !Problem->BuildFlux(fKL) ) continue;

							cK = fKL->BackCell();
							cL = fKL->FrontCell();
							iK = C.Index(cK);
							iL = cL.isValid() ? C.Index(cL) : ENUMUNDEF;
							
							//zero out variables
							fluxT = 0.0;
							corrK = 0.0;
							corrL = 0.0;
							corrK_cK = 0.0;
							corrL_cL = 0.0;

							//advection part
							Problem->AdvectionFluxes(C,fKL,fluxT,corrK,corrL,corrK_cK,corrL_cL);

							if( !joint_advection_diffusion )
							{
								Problem->Averaging(scheme_type,regularization,fluxT,corrK,corrL,corrK_cK,corrL_cL,fluxK,fluxL, iL == ENUMUNDEF);
								if( iL != ENUMUNDEF )
								{
									//check coefficients in the flux
									if( check_fluxes )
									{
										check_flux_properties(iK,fluxK,"fluxK");
										check_flux_properties(iL,-fluxL,"fluxL");
									}
									//add to back cell equation
									if( cK.GetStatus() != Element::Ghost )
									{
										R.Lock(iK);
										R[iK] += fluxK*TimeStep;
										R.UnLock(iK);
									}
									
									//add to front cell equation
									if( cL.GetStatus() != Element::Ghost )
									{
										R.Lock(iL);
										R[iL] -= fluxL*TimeStep;
										R.UnLock(iL);
									}
									//accumulate flux on interfaces
									fKL.Real(tag_FLUX) += get_value(fluxK+fluxL)*0.5*TimeStep;
								}
								else
								{
									//check coefficients in the flux
									if( check_fluxes ) check_flux_properties(iK,fluxK,"fluxK(bnd)");

									if( cK.GetStatus() != Element::Ghost )
									{
										R.Lock(iK);
										R[iK] += fluxK*TimeStep;
										R.UnLock(iK);
									}
									//accumulate flux on interfaces
									fKL.Real(tag_FLUX) += get_value(fluxK)*TimeStep;
								}
								fluxT = 0.0;
								corrK = 0.0;
								corrL = 0.0;
								corrK_cK = 0.0;
								corrL_cL = 0.0;
							}

							//diffusion part
							Problem->DiffusionFluxes(C,fKL,fluxT,corrK,corrL,corrK_cK,corrL_cL);

							Problem->Averaging(scheme_type,regularization,fluxT,corrK,corrL,corrK_cK,corrL_cL,fluxK,fluxL, iL == ENUMUNDEF);			
							if( iL != ENUMUNDEF )
							{
								//check coefficients in the flux
								if( check_fluxes )
								{
									check_flux_properties(iK,fluxK,"fluxK");
									check_flux_properties(iL,-fluxL,"fluxL");
								}
								//add to back cell equation
								if( cK.GetStatus() != Element::Ghost )
								{
									R.Lock(iK);
									R[iK] += fluxK*TimeStep;
									R.UnLock(iK);
								}
									
								//add to front cell equation
								if( cL.GetStatus() != Element::Ghost )
								{
									R.Lock(iL);
									R[iL] -= fluxL*TimeStep;
									R.UnLock(iL);
								}
								//accumulate flux on interfaces
								fKL.Real(tag_FLUX) += get_value(fluxK+fluxL)*0.5*TimeStep;
							}
							else
							{			
								//check coefficients in the flux
								if( check_fluxes ) check_flux_properties(iK,fluxK,"fluxK(bnd)");

								if( cK.GetStatus() != Element::Ghost )
								{
									R.Lock(iK);
									R[iK] += fluxK*TimeStep;
									R.UnLock(iK);
								}
								//accumulate flux on interfaces
								fKL.Real(tag_FLUX) += get_value(fluxK)*TimeStep;
							}
						} //end of loop over faces
					}
					
					if( !steady_state ) //add time derivative and explicit part of the step
					{
						
#if defined(USE_OMP)
#pragma omp parallel
#endif
						{
							variable	flux,  //time flux
										corr;  //high order correction
							real slope, phi;
#if defined(USE_OMP)
#pragma omp for
#endif
							for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
							{
								Cell cell = m->CellByLocalID(q);
								if( cell.GetStatus() == Element::Ghost ) 
									continue;

								real porocity = 1.0;

								if( tag_PORO.isValid() ) porocity = cell->Real(tag_PORO);
								
								flux = C[cell] - C0[cell];
								
								if( bdf2_stepping && step > 0 )//BDF2 time-stepping with limiter
								{
									corr = (C[cell] - 2.0*C0[cell] + C1[cell])*0.5; // hi-lo
									if( fabs(C0[cell] - C1[cell]) > 1.0e-14 )
									{
										slope = get_value(flux/(C0[cell] - C1[cell]));
										phi = (slope+fabs(slope))/(1.0+fabs(slope)); //van-leer
										//phi = (slope*slope+slope)/(slope*slope+1); //van-albada
										//phi = std::max(0.0,std::min(1.0,slope)); //minmod
										//phi = std::max(0.0,std::max(std::min(2.0*slope,1.0), std::min(slope,2.0))); //superbee
									}
									else if( fabs(flux) > 1.0e-14 && flux*(C0[cell] - C1[cell]) > 0.0 ) // slope = +inf
									{
										phi = 2.0; //van-leer,superbee maximum
										//phi = 1.0; //van-albada,minmod maximum
									}
									else phi = 0;
									
									cell->Real(slope_limiter) = phi;
									
									flux += phi*corr; //lo + phi*(hi-lo)
								}
								
								
								R[C.Index(cell)] += flux*cell->Volume()*porocity;
							}
						}
					}
					
					
					if( tag_R.isValid() ) //add reaction term
					{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
						{
							Cell cell = m->CellByLocalID(q);
							if( cell.GetStatus() == Element::Ghost ) 
								continue;
							if( cell->HaveData(tag_R) ) R[C.Index(cell)] += cell->Real(tag_R)*C[cell]*cell->Volume()*TimeStep;
						}
					}
					/*
					//write residual without force into mesh
#if defined(USE_OMP)
#pragma omp parallel for
#endif
					for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
					{
						Cell cell = m->CellByLocalID(q);
						if( cell.GetStatus() == Element::Ghost ) 
								continue;
						cell->Real(tag_residual_no_force) = get_value(R[C.Index(cell)])/cell->Volume();
					}
					m->ExchangeData(tag_residual_no_force,CELL,0);
					*/
					if( tag_F.isValid() ) //add forcing term
					{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
						{
							Cell cell = m->CellByLocalID(q);
							if( cell.GetStatus() == Element::Ghost ) 
								continue;
							if( cell->HaveData(tag_F) ) R[C.Index(cell)] -= cell->Real(tag_F)*cell->Volume()*TimeStep;
						}
					}
					
					//write residual with force into mesh
#if defined(USE_OMP)
#pragma omp parallel for
#endif
					for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
					{
						Cell cell = m->CellByLocalID(q);
						if( cell.GetStatus() == Element::Ghost ) 
								continue;
						cell->Real(tag_residual) = get_value(R[C.Index(cell)])/cell->Volume();
					}
					m->ExchangeData(tag_residual,CELL,0);
					
					//std::cout << "assembled in " << Timer() - tttt << "\t\t\t" << std::endl;
					
					norm = R.Norm();
					if( nit == 0 ) norm0 = norm;
					
					if( m->GetProcessorRank() == 0 )
					{
						std::cout << "Nonlinear residual: " << std::setw(12) << norm << " relative: " << std::setw(12) << norm/norm0 << "\r";
						std::cout.flush();
					}
					
					if( norm < nonlinear_abs_tolerance || norm/norm0 < nonlinear_rel_tolerance )
					{
						if( m->GetProcessorRank() == 0 )
							std::cout << std::endl;
						//std::cout << "Processor " << m->GetProcessorRank() << " exit iterations" << std::endl;
						break;
					}
					
					if( check_matrix )
						check_matrix_properties(R.GetJacobian());
					
					//R.Rescale();
					//R.GetJacobian().Save("jacobian.mtx");
					//R.GetResidual().Save("residual.mtx");
					
					Solver S(linear_solver_type);
					
					S.SetParameterReal("relative_tolerance", 1.0e-18);
					S.SetParameterReal("absolute_tolerance", 1.0e-15);
					S.SetParameterReal("drop_tolerance", 1.0e-2);
					S.SetParameterReal("reuse_tolerance", 1.0e-3);
					S.SetMatrix(R.GetJacobian());
					std::fill(Update.Begin(),Update.End(),0.0);
					if( S.Solve(R.GetResidual(),Update) )
					{
						
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
						{
							Cell cell = m->CellByLocalID(q);
							if( cell.GetStatus() == Element::Ghost ) 
								continue;
							cell->Real(tag_C_Update) = Update[C.Index(cell)];
						}
						m->ExchangeData(tag_C_Update, CELL, 0);
						INMOST_DATA_REAL_TYPE alpha = 1;


						
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
						{
							Cell cell = m->CellByLocalID(q);
							cell->Real(tag_C) -= alpha*cell->Real(tag_C_Update);
						}
						 
						/*
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
						{
							Cell cell = m->CellByLocalID(q);
							if( cell.GetStatus() == Element::Ghost )
								continue;
							cell->Real(tag_C) -= Update[C.Index(cell)];
						}
						m->ExchangeData(tag_C, CELL, 0);
						 */
						
						Cmax = -1.0e20;
						Cmin = 1.0e20;
						Cint = 0;
						
						for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
						{
							Cell cell = m->CellByLocalID(q);
							if( cell->GetStatus() == Element::Ghost )
								continue;
							real C = cell->Real(tag_C);
							if( C > Cmax ) Cmax = C;
							if( C < Cmin ) Cmin = C;
							Cint += C*cell->Volume();
						}

						Cmax = m->AggregateMax(Cmax);
						Cmin = m->AggregateMin(Cmin);
						Cint = m->Integrate(Cint);
						
						if( m->GetProcessorRank() == 0 )
						{
							std::cout << "Nonlinear residual: " << std::setw(12) << norm;
							std::cout << " relative: " << std::setw(12) << norm/norm0;
							std::cout << " C in [" << std::setw(10) << Cmin << "," << std::setw(10) << Cmax << "]";
							std::cout << " integral: " << std::setw(10) << Cint << std::endl;
						}
						
						//if( steady_state )
							SaveMesh(m,"iter",nit,output_format);
					}
					else
					{
						if( m->GetProcessorRank() == 0 )
							std::cout << " Unable to solve: " << S.GetReason() << std::endl;
						break;
					}
					++nit;
				} while(nit < nonlinear_iterations); //check the residual norm
				
				CurrentTime += TimeStep;
				step++;
				
				if( m->GetProcessorRank() == 0 )
					std::cout << "Time: " << std::setw(10) << CurrentTime << " Step: " << step << " Solved problem in " << std::setw(10) << Timer() - ttt << " seconds with " << std::setw(4) << nit << " iterations " << std::endl;
				
				if( !steady_state )
					SaveMesh(m,"step",step,output_format);
			}
		}
		
		Storage::integer check_reference_solution = m->HaveTag("REFERENCE_SOLUTION") ? 1 : 0;
		check_reference_solution = m->AggregateMax(check_reference_solution);
		
		if( check_reference_solution )
		{
			Tag tag_E = m->CreateTag("ERROR",DATA_REAL,CELL,NONE,1);
			Tag tag_Q = m->HaveTag("REFERENCE_SOLUTION") ? m->GetTag("REFERENCE_SOLUTION") : Tag();
			real C, L2, volume;
			C = L2 = volume = 0.0;
			for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
			{
				Cell cell = m->CellByLocalID(q);
				if( cell.GetStatus() == Element::Ghost )
					continue;
				real err = cell.Real(tag_C) - cell.Real(tag_Q);
				real vol = cell.Volume();
				if( C < fabs(err) ) C = fabs(err);
				L2 += err*err*vol;
				volume += vol;
				cell.Real(tag_E) = err;
			}
			C = m->AggregateMax(C);
			L2 = sqrt(m->Integrate(L2)/m->Integrate(volume));
			if( m->GetProcessorRank() == 0 )
				std::cout << "Error of solution, C-norm " << C << " L2-norm " << L2 << std::endl;
		}
		
		Storage::integer check_reference_flux = m->HaveTag("REFERENCE_FLUX") ? 1 : 0;
		check_reference_flux = m->AggregateMax(check_reference_flux);
		
		if( check_reference_flux )
		{
			Tag tag_E = m->CreateTag("ERROR_FLUX",DATA_REAL,FACE,NONE,1);
			Tag tag_Q = m->HaveTag("REFERENCE_FLUX") ? m->GetTag("REFERENCE_FLUX") : Tag();
			real C, L2, volume;
			C = L2 = volume = 0.0;
			for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
			{
				Face face = m->FaceByLocalID(q);
				if( !Problem->BuildFlux(face) ) continue;
				real err = face.Real(tag_FLUX)/face.Area() - face.Real(tag_Q);
				real vol = 0;
				vol += face.BackCell()->Volume() / face.BackCell().nbAdjElements(FACE);
				vol += face.FrontCell().isValid() ? face.FrontCell().Volume() / face.FrontCell().nbAdjElements(FACE) : 0.0;
				if( C < fabs(err) ) C = fabs(err);
				L2 += err*err*vol;
				volume += vol;
				face.Real(tag_E) = err;
			}
			C = m->AggregateMax(C);
			L2 = sqrt(m->Integrate(L2)/m->Integrate(volume));
			if( m->GetProcessorRank() == 0 )
				std::cout << "Error of flux, C-norm " << C << " L2-norm " << L2 << std::endl;

		}

		m->ReleaseMarker(boundary_face,FACE);
		
		SaveMesh(m,"out",output_format);

		delete Problem; //cleanup problem
		delete m; //clean up the mesh
	}
	else
	{
		std::cout << argv[0] << " mesh_file scheme<0:spu,1:nonlinear>" << std::endl;
	}
	
#if defined(USE_PARTITIONER)
	Partitioner::Finalize(); // Finalize the partitioner activity
#endif
	Solver::Finalize(); // Finalize solver and close MPI activity
	return 0;
}
