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

#undef USE_OMP
//#define OPTIMIZATION

int main(int argc,char ** argv)
{
    Solver::Initialize(&argc,&argv,""); // Initialize the solver and MPI activity
#if defined(USE_PARTITIONER)
    Partitioner::Initialize(&argc,&argv); // Initialize the partitioner activity
#endif
    if( argc > 1 )
    {
		std::string force_name = "FORCE";
		int force_comp = 0;
        double ttt; // Variable used to measure timing
        bool repartition = false; // Is it required to redistribute the mesh?
        Mesh * m = new Mesh(); // Create an empty mesh
        { // Load the mesh
            ttt = Timer();
            //~ m->SetParallelFileStrategy(0);
			//~ m->SetParallelStrategy(1);
            //~ m->SetCommunicator(INMOST_MPI_COMM_WORLD); // Set the MPI communicator for the mesh
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

		if( argc > 2 ) force_name = std::string(argv[2]);
		if( argc > 3 ) force_comp = atoi(argv[3]);
		
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
            table[BARYCENTER]  = CELL | FACE; //Compute volumetric center of mass
            m->PrepareGeometricData(table); //Ask to precompute the data
            BARRIER
            if( m->GetProcessorRank() == 0 ) std::cout << "Prepare geometric data: " << Timer()-ttt << std::endl;
        }

        // data tags for
        Tag tag_P;  // Pressure
        Tag tag_K;  // Diffusion tensor
        Tag tag_F;  // Forcing term
        Tag tag_BC; // Boundary conditions
        TagRealArray     tag_WA; // Approximation mimetic matrix
        TagRealArray     tag_WS; // Stability mimetic matrix
        TagVariableArray tag_WF; // Approximation of fluxes
        
		
		
		//~ TagRealArray tag_PG; // Pressure gradient
		//~ TagRealArray tag_WG; // matrix to reconstruct gradient

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
                
            //~ tag_PG = m->CreateTag("PRESSURE_GRADIENT",DATA_REAL,CELL,NONE,3);
            //~ tag_WG = m->CreateTag("WGRAD",DATA_REAL,CELL,NONE);

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
            tag_WA = m->CreateTag("WA",DATA_REAL,CELL,NONE);
            tag_WS = m->CreateTag("WS",DATA_REAL,CELL,NONE);
            tag_WF = m->CreateTag("WF",DATA_VARIABLE,FACE,NONE,2);
			
            ttt = Timer();
            //Assemble gradient matrix W on cells
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				rMatrix N, L, R, W;
				rMatrix x(1,3), xf(1,3), n(1,3);
				double area, dist; //area of the face
#if defined(USE_OMP)
#pragma omp for
#endif
				for( integer q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
				{
					Cell cell = m->CellByLocalID(q);
					ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
					int NF = (int)faces.size(); //number of faces;
					cell->Barycenter(x.data());
					//get permeability for the cell
					rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),
													cell->RealArrayDF(tag_K).size());
					N.Resize(NF,3); //co-normals
					R.Resize(NF,3); //directions
					L.Resize(NF,NF); //two-point transmissibility coefficient
					L.Zero();
					tag_WS[cell].resize(NF);
					for(int k = 0; k < NF; ++k) //loop over faces
					{
						area = faces[k].Area();
						faces[k].Barycenter(xf.data());
						faces[k].OrientedUnitNormal(cell->self(),n.data());
						dist = (n.DotProduct(xf-x));
						if( dist < 0 ) std::cout << "negative dist! " << dist << " cell " << q << " face " << k << std::endl;
						// assemble matrix of directions
						R(k,k+1,0,3) = (xf-x);
						// assemble matrix of co-normals
						N(k,k+1,0,3) = n*K*area;
						L(k,k) =  n.DotProduct(n*K)*area/dist;
						tag_WS[cell][k] = L(k,k);
					} //end of loop over faces
					 
			 		//~ double ca = 0, ba = ca; //parameters
			 		tag_WA[cell].resize(NF*NF);
			 		
			 		//~ tag_WG[cell].resize(3*NF);
			 		//tag_WA(cell,NF,NF) = (N+ca*L*R)*((N+ca*L*R).Transpose()*R).Invert()*(N+ca*L*R).Transpose();
			 		//tag_WS(cell,NF,NF) = L - (1+ba)*(L*R)*((L*R).Transpose()*R).Invert()*(L*R).Transpose();
			 		tag_WA(cell,NF,NF) = N*(N.Transpose()*R).Invert()*N.Transpose() - (L*R)*((L*R).Transpose()*R).Invert()*(L*R).Transpose();
			 		/*
			 		#pragma omp critical
			 		{
			 		std::cout << "cell " << cell.LocalID() << std::endl;
			 		std::cout << "WA: " << std::endl;
			 		tag_WA(cell,NF,NF).Print();
			 		for(int k = 0; k < NF; ++k)
			 		{
						double rowsum = 0;
						for(int q = 0; q < NF; ++q)
							rowsum += tag_WA(cell,NF,NF)(k,q);
						std::cout << k << " T " << L(k,k) << " rowsum " << rowsum+L(k,k) << " no T " << rowsum << std::endl;
					}
					}
					*/
			 		//~ tag_WS(cell,NF,1) = L;
					//~ tag_WG(cell,3,NF) = (N.Transpose()*R).Invert()*N.Transpose();
				 } //end of loop over cells
			}
            std::cout << "Construct W matrix: " << Timer() - ttt << std::endl;
			 
            if( m->HaveTag(force_name) ) //Is there force on the mesh?
            {
                tag_F = m->GetTag(force_name); //initial force
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

            Residual Resid("",aut.GetFirstIndex(),aut.GetLastIndex());
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
				Resid.Clear(); //clean up the residual
				//~ double tttt = Timer();
#if defined(USE_OMP)
#pragma omp parallel
#endif
				{
					vMatrix pF; //vector of pressure differences on faces
					vMatrix FLUX; //computed flux on faces
					variable u, v;
					variable muu, muv;
					//~ double muu,muv;
					//~ double area; //area of the face
#if defined(USE_OMP)
#pragma omp for
#endif
					for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) ) //loop over cells
					{
						Cell cell = m->CellByLocalID(q);
						ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
						int NF = (int)faces.size();
						pF.Resize(NF,1);
						for(int k = 0; k < NF; ++k)
							pF(k,0) = (P(faces[k]) - P(cell));
						FLUX = tag_WA(cell,NF,NF)*pF; //approximation of fluxes on faces
						//~ tag_PG(cell,3,1) = tag_WG(cell,3,NF)*pF; //gradient of function
						for(int k = 0; k < NF; ++k) //loop over faces of current cell
						{
							int ind = (faces[k].BackCell() == cell ? 0 : 1);
							tag_WF[faces[k]][ind] = FLUX(k,0); //memorize fluxes
						}
					} //end of loop over cells
					
#if defined(USE_OMP)
#pragma omp for
#endif
					for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) ) //loop over cells
					{
						Cell cell = m->CellByLocalID(q);
						ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
						int NF = (int)faces.size();
						//get permeability for the cell
						/*
						cell->Barycenter(x.data());
						rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),
														cell->RealArrayDF(tag_K).size());
						N.Resize(NF,3); //co-normals
						R.Resize(NF,3); //directions
						coefs.Resize(NF,1);
						fluxes.Resize(NF,2);
						for(int k = 0; k < NF; ++k) //loop over faces
						{
							area = faces[k].Area();
							faces[k].Barycenter(xf.data());
							faces[k].OrientedUnitNormal(cell->self(),n.data());
							// assemble matrix of directions
							R(k,k+1,0,3) = (xf-x);
							// assemble matrix of co-normals
							N(k,k+1,0,3) = n*K*area;
							
							int ind = (faces[k].BackCell() == cell ? 0 : 1);
							FLUX(k,0); = tag_WF[faces[k]][ind] //restore fluxes
							if( !faces[k].Boundary() )
							{
								const double eps = 1.0e-6;
								
								// f = ((u^2+e) v + (v^2+e) u) / (u^2 + v^2 + 2e) = u * ( (u + e/u)*v + (v^2+e) ) / (u^2 + v^2 + 2e)
								//~ variable coef = variation((tag_WF[faces[k]][0]*tag_WF[faces[k]][1] + pow(tag_WF[faces[k]][1-ind],2)+eps)/(pow(tag_WF[faces[k]][0],2)+pow(tag_WF[faces[k]][1],2)+2*eps),0);
								variable coef = variation(2*sqrt(pow(tag_WF[faces[k]][1-ind],2)+eps*eps)/(sqrt(pow(tag_WF[faces[k]][0],2)+eps*eps)+sqrt(pow(tag_WF[faces[k]][1],2)+eps*eps)),0.0);
								coefs(k,0) = get_value(coef);
								N(k,k+1,0,3) *= coef;
								fluxes(k,0) = get_value(tag_WF[faces[k]][ind]);
								fluxes(k,1) = get_value(tag_WF[faces[k]][1-ind]);
							}
							else 
							{
								fluxes(k,0) = get_value(tag_WF[faces[k]][ind]);
								fluxes(k,1) = 0;
								coefs(k,0) = 1;
							}
						} //end of loop over faces
						*/
						for(int k = 0; k < NF; ++k)
						{
							int ind = (faces[k].BackCell() == cell ? 0 : 1);
							//FLUX(k,0) = tag_WF[faces[k]][ind];
							u = tag_WF[faces[k]][ind];
							if( !faces[k].Boundary() )
							{
								v = tag_WF[faces[k]][1-ind];
								//~ muu = fabs(get_value(v)) + 1.0e-5;
								//~ muv = fabs(get_value(u)) + 1.0e-5;
								//~ muu = get_value(v)*get_value(v) + 1.0e-5;
								//~ muv = get_value(u)*get_value(u) + 1.0e-5;
								//~ muu = 0.5;
								//~ muv = 0.5;
								muu = sqrt(v*v + 1.0e-4);
								muv = sqrt(u*u + 1.0e-4);
								//~ muu = variation(sqrt(v*v + 1.0e-10),0.25);
								//~ muv = variation(sqrt(u*u + 1.0e-10),0.25);
								//~ muu = v*v + 1.0e-2;
								//~ muv = u*u + 1.0e-2;
								//~ FLUX(k,0) = u*variation(muu/(muu+muv),0.1) - v*variation(muv/(muu+muv),0.1);
								FLUX(k,0) = (u*muu - v*muv)/(muu+muv);
								//~ FLUX(k,0) = 2*u*variation(muu/(muu+muv),0.5);
							}
							else FLUX(k,0) = u;
						}
						for(int k = 0; k < NF; ++k)
							FLUX(k,0) += tag_WS[cell][k]*(P(faces[k]) - P(cell));
						if( cell.GetStatus() != Element::Ghost )
						{
							for(int k = 0; k < NF; ++k) //loop over faces of current cell
								Resid[P.Index(cell)] += FLUX(k,0);
						}
						for(int k = 0; k < NF; ++k) //loop over faces of current cell
						{
							if( faces[k].GetStatus() == Element::Ghost ) continue;
							int index = P.Index(faces[k]);
							Locks.Lock(index);
							if( tag_BC.isValid() && faces[k].HaveData(tag_BC) )
							{
								real_array BC = faces[k].RealArray(tag_BC);
								Resid[index] -= BC[0]*P(faces[k]) + BC[1]*FLUX(k,0) - BC[2];
							}
							else
							{
								
								Resid[index] -= FLUX(k,0);
							}
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
							 if( cell->HaveData(tag_F) ) Resid[P.Index(cell)] += cell->RealArray(tag_F)[force_comp]*cell->Volume();
						 }
					 }
					 
					 //R[P.Index(m->BeginCell()->self())] = P[m->BeginCell()->self()];
				}
				//~ std::cout << "assembled in " << Timer() - tttt << "\t\t\t" << std::endl;

                std::cout << "Nonlinear residual: " << Resid.Norm() << "                 " << std::endl;

                if( Resid.Norm() < 1.0e-4 ) break;

				//~ Solver S(Solver::INNER_ILU2);
                Solver S(Solver::INNER_MPTILUC);
                //~ Solver S(Solver::K3BIILU2);
				//Solver S("superlu");
				S.SetParameter("verbosity","1");
                S.SetParameter("relative_tolerance", "1.0e-14");
                S.SetParameter("absolute_tolerance", "1.0e-12");
                S.SetParameter("drop_tolerance",  "1.0e-3");
                S.SetParameter("reuse_tolerance", "1.0e-4");

                S.SetMatrix(Resid.GetJacobian());
				
                if( S.Solve(Resid.GetResidual(),Update) )
                {
					double alpha = 1;
#if defined(USE_OMP)
#pragma omp parallel for
#endif
                    for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
                    {
                        Cell cell = m->CellByLocalID(q);
						if( cell->GetStatus() == Element::Ghost ) continue;
                        cell->Real(tag_P) -= Update[P.Index(cell)]*alpha;
                    }
#if defined(USE_OMP)
#pragma omp parallel for
#endif
                    for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
                    {
                        Face face = m->FaceByLocalID(q);
						if( face->GetStatus() == Element::Ghost ) continue;
                        face->Real(tag_P) -= Update[P.Index(face)]*alpha;
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
            } while( Resid.Norm() > 1.0e-4 && nit < 50); //check the residual norm
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
            
        m->Save("solution.pmf");
        m->Save("solution.vtu");

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
