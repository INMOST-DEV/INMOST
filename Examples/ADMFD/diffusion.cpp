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
bool output_matrix = false;
bool norne_wells = true;

//data for wells
Cell cc[3] = {InvalidCell(),InvalidCell(),InvalidCell()};
real WI[3] = {50000,50000,50000};
real pbhp[3] = {265,105,110};
real ccnt[3][3] =
{
	{4.567151e+05, 7.321079e+06, 2.767665e+03},
	{4.609346e+05, 7.323503e+06, 2.597767e+03},
	{4.595400e+05, 7.326078e+06, 2.803586e+03}
};


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

        //{ // prepare geometrical data on the mesh
		ttt = Timer();
		Mesh::GeomParam table;
		//~ table[CENTROID]    = CELL | FACE; //Compute averaged center of mass
		table[NORMAL]      = FACE;        //Compute normals
		table[ORIENTATION] = FACE;        //Check and fix normal orientation
		table[MEASURE]     = CELL | FACE; //Compute volumes and areas
		table[BARYCENTER]  = CELL | FACE; //Compute volumetric center of mass
		m->PrepareGeometricData(table); //Ask to precompute the data
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Prepare geometric data: " << Timer()-ttt << std::endl;
        //}

        // data tags for
        Tag tag_P;  // Pressure
        Tag tag_K;  // Diffusion tensor
        Tag tag_F;  // Forcing term
        Tag tag_W;  // Gradient matrix acting on harmonic points on faces and returning gradient on faces
        TagRealArray tag_BC; // Boundary conditions
        
		
		
		TagRealArray tag_PG; // Pressure gradient
		TagRealArray tag_WG; // matrix to reconstruct gradient

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
                tag_K = m->CreateTag("PERM",DATA_REAL,CELL,NONE,1); // create a new tag for symmetric diffusion tensor K
                for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) ) // loop over mesh cells
                {
                    Cell cell = m->CellByLocalID(q);
                    real_array K = cell->RealArray(tag_K);
                    // assign a symmetric positive definite tensor K
                    K[0] = 1.0; //XX
                    //~ K[1] = 0.0; //XY
                    //~ K[2] = 0.0; //XZ
                    //~ K[3] = 1.0; //YY
                    //~ K[4] = 0.0; //YZ
                    //~ K[5] = 1.0; //ZZ
                }

                m->ExchangeData(tag_K,CELL,0); //Exchange diffusion tensor
            }
            
            if( norne_wells )
            {
				for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) if( it->GetStatus() != Element::Ghost )
				{
					for(int q = 0; q < 3; ++q) if( it->Inside(ccnt[q]) )
					{
						cc[q] = it->self();
						std::cout << "proc " << m->GetProcessorRank() << " found c" << q << " " << cc[q].LocalID() << std::endl;
					}
				}
			}

            if( m->HaveTag("PRESSURE") ) //Is there a pressure on the mesh?
                tag_P = m->GetTag("PRESSURE"); //Get the pressure
                
            tag_PG = m->CreateTag("PRESSURE_GRADIENT",DATA_REAL,CELL,NONE,3);
            tag_WG = m->CreateTag("WGRAD",DATA_REAL,CELL,NONE);

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
			
            ttt = Timer();
            //Assemble gradient matrix W on cells
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				rMatrix NK, L, R, Areas;
				rMatrix x(1,3), xf(1,3), n(1,3);
				double area, dist; //area of the face
#if defined(USE_OMP)
#pragma omp for
#endif
				 for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
				 {
					 Cell cell = m->CellByLocalID(q);
					 ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
					 int NF = (int)faces.size(); //number of faces;
					 rMatrix W(NF,NF);
					 cell->Barycenter(x.data());
					 //get permeability for the cell
					 rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),
													 cell->RealArrayDF(tag_K).size());
					 NK.Resize(NF,3); //co-normals
					 R.Resize(NF,3); //directions
					 L.Resize(NF,NF);
					 Areas.Resize(NF,NF); //areas
					 Areas.Zero();
					 L.Zero();
					 for(int k = 0; k < NF; ++k) //loop over faces
					 {
						 area = faces[k].Area();
						 faces[k].Barycenter(xf.data());
						 faces[k].OrientedUnitNormal(cell->self(),n.data());
						 dist = n.DotProduct(xf-x);
						 // assemble matrix of directions
						 R(k,k+1,0,3) = (xf-x);
						 // assemble matrix of co-normals
						 NK(k,k+1,0,3) = n*K*area;
						 L(k,k) =  n.DotProduct(n*K)*area/dist;
						 Areas(k,k) = area;
					 } //end of loop over faces
					 //~ W = NK*(NK.Transpose()*R).PseudoInvert(1.0e-12)*NK.Transpose(); //stability part
					 //~ W+=(rMatrix::Unit(NF) - R*(R.Transpose()*R).CholeskyInvert()*R.Transpose())*
						//~ (2.0/(static_cast<real>(NF)*volume)*(NK*K.CholeskyInvert()*NK.Transpose()).Trace());
			 		 
			 		 
			 		 W = (NK)*((NK).Transpose()*R).Invert()*(NK).Transpose();
					 W+= L - (L*R)*((L*R).Transpose()*R).Invert()*(L*R).Transpose();
					 
					 //~ W = Areas*W*Areas;
					 //access data structure for gradient matrix in mesh
					 real_array store_W = cell->RealArrayDV(tag_W);
					 //resize the structure
					 store_W.resize(NF*NF);
					 //write down the gradient matrix
					 std::copy(W.data(),W.data()+NF*NF,store_W.data());
					 
					 tag_WG[cell].resize(3*NF);
					 tag_WG(cell,3,NF) = (NK.Transpose()*R).PseudoInvert(1.0e-12)*(NK).Transpose();
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
            
            MarkerType unk = m->CreateMarker();
            for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
				m->CellByLocalID(q).SetMarker(unk);
			for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
			{
				Face face = m->FaceByLocalID(q);
				face.SetMarker(unk);
				if( face.Boundary() )
				{
					double alpha = 0, beta = 1, gamma = 0;
					if( tag_BC.isValid() && face.HaveData(tag_BC) )
					{
						alpha = tag_BC[face][0];
						beta  = tag_BC[face][1];
						gamma = tag_BC[face][2];
					}
					if( alpha != 0 && beta == 0 )
					{
						face.RemMarker(unk);
						face.Real(tag_P) = gamma/alpha;
					}
				}
			}
            
            dynamic_variable P(aut,aut.RegisterTag(tag_P,CELL|FACE,unk)); //register pressure as primary unknown
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
                    if( face.GetStatus() != Element::Ghost && face.GetMarker(unk) )
                    {
                        if( face.Boundary() )
                        {
							double alpha = 0, beta = 1, gamma = 0;
							if( tag_BC.isValid() && face.HaveData(tag_BC) )
							{
								alpha = tag_BC[face][0];
								beta  = tag_BC[face][1];
								gamma = tag_BC[face][2];
							}
							std::stringstream str;
							str << "Pressure guided by boundary condition " << alpha << " " << beta << " " << gamma;
                            Text.SetAnnotation(P.Index(face),str.str());
						}
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
                //~ int total = 0, dmp = 0;
#if defined(USE_OMP)
#pragma omp parallel
#endif
				{
					vMatrix pF; //vector of pressure differences on faces
					vMatrix FLUX; //computed flux on faces
					rMatrix n(3,1);
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
							pF(k,0) = (P[faces[k]] - P[cell]);
						FLUX = W*pF; //fluxes on faces
						tag_PG(cell,3,1) = tag_WG(cell,3,NF)*pF;
						if( cell.GetStatus() != Element::Ghost )
						{
							for(int k = 0; k < NF; ++k) //loop over faces of current cell
								R[P.Index(cell)] += FLUX(k,0);//faces[k].Area();
						}
						for(int k = 0; k < NF; ++k) //loop over faces of current cell
						{
							if( faces[k].GetStatus() == Element::Ghost ) continue;
							if( !faces[k].GetMarker(unk) ) continue;
							int index = P.Index(faces[k]);
							Locks.Lock(index);
							
							if( faces[k].Boundary() )
							{
								double a = 0, b = 1, c = 0;
								//~ faces[k].UnitNormal(n.data());
								//~ double a = 1, b = 0, c = 0;
								//~ if( fabs(n(2,0)-1) < 1.0e-4 )
								//~ {
									//~ a = 0;
									//~ b = 1;
								//~ }
								
								if( tag_BC.isValid() && faces[k].HaveData(tag_BC) )
								{
									real_array BC = faces[k].RealArray(tag_BC);
									a = BC[0], b = BC[1], c = BC[2];
								}
								R[index] -= a*P(faces[k]) + b*FLUX(k,0) - c;
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
							 if( cell->HaveData(tag_F) ) R[P.Index(cell)] += cell->RealArray(tag_F)[force_comp]*cell->Volume();
						 }
					 }
					 
					 if( norne_wells )
					 {
						 for(int q = 0; q < 3; ++q) if( cc[q].isValid() )
							R[P.Index(cc[q])] += WI[q]*(pbhp[q] - P[cc[q]])*cc[q].Volume();
					 }
					 
					 //R[P.Index(m->BeginCell()->self())] = P[m->BeginCell()->self()];
				}
				std::cout << "assembled in " << Timer() - tttt << "\t\t\t" << std::endl;

                std::cout << "Nonlinear residual: " << R.Norm() << "\t\t" << std::endl;

                if( R.Norm() < 1.0e-4 ) break;

				//~ R.GetJacobian().Save("A.mtx");
				//~ R.GetResidual().Save("b.mtx");
				//Solver S(Solver::INNER_ILU2);
                Solver S(Solver::INNER_MPTILUC);
                //~ Solver S(Solver::INNER_MLMPTILUC);
                //~ Solver S(Solver::K3BIILU2);
				//Solver S("superlu");
				S.SetParameter("verbosity","3");
				S.SetParameter("rescale_iterations", "10");
                S.SetParameter("relative_tolerance", "1.0e-14");
                S.SetParameter("absolute_tolerance", "1.0e-12");
                //~ S.SetParameter("drop_tolerance", "5.0e-2");
                //~ S.SetParameter("reuse_tolerance", "2.5e-3");
                S.SetParameter("drop_tolerance", "1.0e-2");
                S.SetParameter("reuse_tolerance", "1.0e-4");
                S.SetParameter("pivot_condition", "20");
                double tset = Timer(), titr;

                S.SetMatrix(R.GetJacobian());
                
                tset = Timer() - tset;
                
				if( output_matrix )
				{
					std::cout << "write A.mtx" << std::endl;
					R.GetJacobian().Save("A.mtx");//,&Text);
					std::cout << "write b.mtx" << std::endl;
					R.GetResidual().Save("b.mtx");
					std::cout << "write done, solve" << std::endl;
				}
				
				titr = Timer();
				
				bool success =  S.Solve(R.GetResidual(),Update);
				
				titr = Timer() - titr;
				
                if( success )
                {
					std::cout << "Solved system in " << S.Iterations() << " iterations, time for setup " << tset << "s iterations " << titr << "s" << std::endl;
					if( output_matrix )
					{
						std::cout << "write x.mtx" << std::endl;
						Update.Save("x.mtx");
						//std::cout << "exit" << std::endl;
						//exit(-1);
					}
#if defined(USE_OMP)
#pragma omp parallel for
#endif
                    for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
                    {
                        Cell cell = m->CellByLocalID(q);
						if( cell->GetStatus() == Element::Ghost ) continue;
						if( !cell->GetMarker(unk) ) continue;
                        cell->Real(tag_P) -= Update[P.Index(cell)];
                    }
#if defined(USE_OMP)
#pragma omp parallel for
#endif
                    for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
                    {
                        Face face = m->FaceByLocalID(q);
						if( face->GetStatus() == Element::Ghost ) continue;
						if( !face->GetMarker(unk) ) continue;
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
        
        tag_W = m->DeleteTag(tag_W);
        tag_WG = m->DeleteTag(tag_WG);
        tag_K = m->DeleteTag(tag_K);
        tag_BC = m->DeleteTag(tag_BC);
        m->RemoveGeometricData(table);

        if( m->GetProcessorsNumber() == 1 )
            m->Save("out.vtk");
        else
            m->Save("out.pvtk");
            
        m->Save("solution.pmf");

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
