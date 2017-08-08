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

bool rt0 = false;

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
		if( argc > 2 )
		{
			if( std::string(argv[2]) == "MFD" )
				rt0 = false;
			else if ( std::string(argv[2]) == "MHFE" )
				rt0 = true;
		}

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

		if( rt0 )
			std::cout << "Running MHFE RT0" << std::endl;
		else
			std::cout << "Running MFD" << std::endl;

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
				Mesh * mesh = m;
                Cell cell = m->CellByLocalID(q);
				ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell

				int NF = (int)faces.size(); //number of faces;	
				rMatrix W(NF,NF);
				
				if( cell->GetGeometricType() == Element::Tet && rt0 ) // RT0 consturction of W matrix
				{
					double V = cell.Volume();
					double dN[12];
					double J[9];
					double J_T[9];
					double K_inv[9];
					double K_inv_ref1[9];
					double W_RT0[12];
					double B_RT0[16];
					double gauss_pt_xyz[12];
					double gauss_wei[4];
					double tetra[12];
					double xyz_i[3];
					double xyz_j[3];
					int i, j,k,Z;
					double w,l,m, J_det,x,y,z;

					rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),cell->RealArrayDF(tag_K).size()); //get permeability for the cell
					rMatrix K_inverse (3,3); // inverse of permeability 

					K_inverse= K.Invert(true).first;
					// gauss points for intgration
					gauss_pt_xyz[0]=(5.0+3.0*(sqrt(5.0)))/20.0;
					gauss_pt_xyz[1]=(5.0-(sqrt(5.0)))/20.0;
					gauss_pt_xyz[2]=(5.0-(sqrt(5.0)))/20.0;

					gauss_pt_xyz[3]=(5.0-(sqrt(5.0)))/20.0;
					gauss_pt_xyz[4]=(5.0+3.0*(sqrt(5.0)))/20.0;
					gauss_pt_xyz[5]=(5.0-(sqrt(5.0)))/20.0;

					gauss_pt_xyz[6]=(5.0-(sqrt(5.0)))/20.0;
					gauss_pt_xyz[7]=(5.0-(sqrt(5.0)))/20.0;
					gauss_pt_xyz[8]=(5.0+3.0*(sqrt(5.0)))/20.0;

					gauss_pt_xyz[9]=(5.0-(sqrt(5.0)))/20.0;
					gauss_pt_xyz[10]=(5.0-(sqrt(5.0)))/20.0;
					gauss_pt_xyz[11]=(5.0-(sqrt(5.0)))/20.0;

					// gauss weights
					gauss_wei[0]=  1.0/24.0;
					gauss_wei[1]=  1.0/24.0;
					gauss_wei[2]=  1.0/24.0;
					gauss_wei[3]= 1.0/24.0;

					// derivative of the finite element shape functions for tetraherdon (local) 
					dN[0]=-1; 
					dN[1]=-1;
					dN[2]=-1;

					dN[3]=1; 
					dN[4]=0;
					dN[5]=0;

					dN[6]=0; 
					dN[7]=1;
					dN[8]=0;

					dN[9]=0; 

					dN[10]=0;
					dN[11]=1;

					// local to reference tranformation
					ElementArray<Node> node(mesh);
					for(j = 0; j < 4; ++j)
					{
						ElementArray<Node> cell_nodes = cell->getNodes();
						cell_nodes.Subtract(faces[j].getNodes());
						node.Unite(cell_nodes);
					}

					for ( j = 0; j < 4; j++ )
					{
						for ( i = 0; i < 3; i++ )
							tetra[i+j*3]= node[j].Coords()[i]; // <-- X
					}

					J[0]=dN[0]*tetra[0]+dN[3]*tetra[3]+dN[6]*tetra[6]+dN[9]*tetra[9];
					J[1]=dN[0]*tetra[1]+dN[3]*tetra[4]+dN[6]*tetra[7]+dN[9]*tetra[10];
					J[2]=dN[0]*tetra[2]+dN[3]*tetra[5]+dN[6]*tetra[8]+dN[9]*tetra[11];

					J[3]=dN[1]*tetra[0]+dN[4]*tetra[3]+dN[7]*tetra[6]+dN[10]*tetra[9];
					J[4]=dN[1]*tetra[1]+dN[4]*tetra[4]+dN[7]*tetra[7]+dN[10]*tetra[10];
					J[5]=dN[1]*tetra[2]+dN[4]*tetra[5]+dN[7]*tetra[8]+dN[10]*tetra[11];

					J[6]=dN[2]*tetra[0]+dN[5]*tetra[3]+dN[8]*tetra[6]+dN[11]*tetra[9];
					J[7]=dN[2]*tetra[1]+dN[5]*tetra[4]+dN[8]*tetra[7]+dN[11]*tetra[10];
					J[8]=dN[2]*tetra[2]+dN[5]*tetra[5]+dN[8]*tetra[8]+dN[11]*tetra[11];

					w = J[0] * (J[4] * J[8] -J [5] * J[7]);
					l = J[1] * (J[3] * J[8] -J [5] * J[6]);
					m = J[2] * (J[3] * J[7] -J [4] * J[6]);
					J_det = w - l + m;
					// cout << " Determinent equals " << J_det << endl;
					if (J_det == 0.0)
					{
						std::cout << "As Determinent=0 so it is singular matrix and its inverse cannot exist for element "
							<< std::endl;

						exit(0);
					}
					// transpose of jacobian
					J_T[0]=J[0]; J_T[1]=J[3]; J_T[2]=J[6];
					J_T[3]=J[1]; J_T[4]=J[4]; J_T[5]=J[7];
					J_T[6]=J[2]; J_T[7]=J[5]; J_T[8]=J[8];

					// inverse of K tensor 
					for ( j = 0; j < 3; j++ )
					{

						for ( i = 0; i < 3; i++ )
						{
							K_inv[i+j*3]=K_inverse (i,j);
						}

					}

					for ( i = 0; i < 3; i++) { // row number of output
						for ( j = 0; j < 3; j++) { // column number of output
							K_inv_ref1[3*i+j] = 0.0;
							for ( k = 0; k < 3; k++) { // three elements are added for this output
								K_inv_ref1[3*i+j] += J[3*i+k] * K_inv[3*k+j];
							}
						}
					}

					for ( i = 0; i < 16; i++) 
					{
						B_RT0[i]=0.0;
					}

					for ( i = 0; i < 4; i++) { // row number of output

						for ( j = 0; j < 4; j++) { // column number of output

							for ( Z = 0; Z < 4; Z++) { // GAUSS POINT IN Z DIRECTION 
								x=gauss_pt_xyz[Z*3];
								y=gauss_pt_xyz[Z*3+1];
								z=gauss_pt_xyz[Z*3+2];



								W_RT0[0]=  2*x;///(sqrt(3.0));
								W_RT0[1]=  2*y;///(sqrt(3.0));
								W_RT0[2]=  2*z;///(sqrt(3.0));

								W_RT0[3]=  2*(-1+x);
								W_RT0[4]=  2*y;

								W_RT0[5]=  2*z;

								W_RT0[6]=  (2*x);///(sqrt(3.0));
								W_RT0[7]=  2*(-1+y);;///(sqrt(3.0));
								W_RT0[8]=  (2*z);///(sqrt(3.0));

								W_RT0[9]=  2*x;

								W_RT0[10]= 2*y;

								W_RT0[11]= 2*(-1+z);





								xyz_i[0]=W_RT0[i*3]*J[0]+W_RT0[i*3+1]*J[3]+W_RT0[i*3+2]*J[6];
								xyz_i[1]=W_RT0[i*3]*J[1]+W_RT0[i*3+1]*J[4]+W_RT0[i*3+2]*J[7];
								xyz_i[2]=W_RT0[i*3]*J[2]+W_RT0[i*3+1]*J[5]+W_RT0[i*3+2]*J[8];



								xyz_j[0]=W_RT0[j*3]*K_inv_ref1[0]+W_RT0[j*3+1]*K_inv_ref1[3]+W_RT0[j*3+2]*K_inv_ref1[6];
								xyz_j[1]=W_RT0[j*3]*K_inv_ref1[1]+W_RT0[j*3+1]*K_inv_ref1[4]+W_RT0[j*3+2]*K_inv_ref1[7];
								xyz_j[2]=W_RT0[j*3]*K_inv_ref1[2]+W_RT0[j*3+1]*K_inv_ref1[5]+W_RT0[j*3+2]*K_inv_ref1[8];


								B_RT0[i+j*4]+=gauss_wei[Z]*(xyz_i[0]*xyz_j[0]+ xyz_i[1]*xyz_j[1]+xyz_i[2]*xyz_j[2])/std::abs(J_det);


							}

						}
					}	

					for ( j = 0; j < 4; j++) { // row number of output
						double sum = 0;
						for ( i = 0; i < 4; i++) { // column number of output

							sum += W(i,j);
							W(i,j)= B_RT0[i+j*4];
						}
						//std::cout << "row " << j << " sum " << sum << std::endl;
					}


					//std::cout << "W" << std::endl;
					//W.Print();

					W = W.Invert(true).first ;


					for ( j = 0; j < 4; j++) { // row number of output
						double sum = 0;
						for ( i = 0; i < 4; i++) { // column number of output

							sum += W(i,j);
						}
						//std::cout << "row inverted " << j << " sum " << sum << std::endl;
					}

					//std::cout << "W inverted" << std::endl;

					//W.Print();

					//scanf("%*c");

				}
				else
				{
					real xP[3]; //center of the cell
					real yF[3]; //center of the face
					real nF[3]; //normal to the face
					real aF; //area of the face
					real vP = cell->Volume(); //volume of the cell
					cell->Centroid(xP);
					rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),cell->RealArrayDF(tag_K).size()); //get permeability for the cell
					//rMatrix U,S,V;
					//K0.SVD(U,S,V);
					//for(int k = 0; k < 3; ++k) S(k,k) = sqrt(S(k,k));
					//rMatrix K = U*S*V;
					rMatrix NK(NF,3), R(NF,3), D(NF,NF), U(NF,NF), Areas(NF,NF); //big gradient matrix, co-normals, directions
					Areas.Zero();
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

						Areas(k,k) = aF;
					} //end of loop over faces
					rMatrix SU,SS,SV;
					W = NK*(NK.Transpose()*R).Invert(true).first*NK.Transpose(); //stability part
					/*
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

					//std::cout << "UDtR" << std::endl;
					//(U*D.Transpose()*R).Print();

					W += D*U*D.Transpose();


					W = Areas*W*Areas;
				}
					//std::cout << "W: " << std::endl;
					//W.Print();
				bulk & isDMP = cell->Bulk(tag_DMP);
				isDMP = 1;
				for(int k = 0; k < NF; ++k)
				{
					real row_sum = 0;
					if( W(k,k) < 0.0 ) isDMP = 0;
					for(int j = 0; j < NF; ++j)
						row_sum += W(k,j);
					if( row_sum < 0.0 ) isDMP = 0;
					for(int j = k+1; j < NF; ++j)
						if( W(k,j) > 0.0 )
							isDMP = 0;
				}
				++total;
				if( isDMP ) ++dmp;
				real_array store_W = cell->RealArrayDV(tag_W); //access data structure for gradient matrix in mesh
				store_W.resize(NF*NF); //resize the structure
				std::copy(W.data(),W.data()+NF*NF,store_W.data()); //write down the gradient matrix
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
            Automatizator aut; // declare class to help manage unknowns
            Automatizator::MakeCurrent(&aut);
            dynamic_variable P(aut,aut.RegisterTag(tag_P,CELL|FACE)); //register pressure as primary unknown
            aut.EnumerateTags(); //enumerate all primary variables
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
                    vMatrix FLUX_MASS(NF,1); //computed flux on faces for mass balance on cell
					vMatrix POT (2,1); //vector of pressure differences on faces

					for(int k = 0; k < NF; ++k)
					{
					
					if (faces[k]->FrontCell().isValid())
					{
						Cell cell_n = cell.Neighbour(faces[k]); 

						ElementArray<Face> faces_n = cell_n->getFaces(); //obtain faces of the cell neighoubr
					    rMatrix nKGRAD_n(cell_n->RealArrayDV(tag_W).data(),NF,NF); //Matrix for gradient
						double B_1, B_2, L_1, L_2;
						int face_n_k;

                        L_1 =0.0; 

						B_1= nKGRAD(k,k) ;
						for(int j = 0; j < NF; ++j)
						L_1 +=nKGRAD(k,j) ;
						for(int j = 0; j < NF; ++j) if (faces_n[j] == faces[k])
						{
							face_n_k= j ;
						}

						L_2 =0.0; 

						for(int j = 0; j < NF; ++j)
						L_2 +=nKGRAD_n(face_n_k,j) ;

						B_2 = nKGRAD_n(face_n_k,face_n_k) ;

						POT(0,0)=L_1*P(cell) ;
						POT(1,0)=L_2*P(cell_n) ;

						for(int j = 0; j < NF; ++j)if (j != k)
						POT(0,0) -= nKGRAD(k,j)* P(faces[j]);

						for(int j = 0; j < NF; ++j) if (j != face_n_k)
						POT(1,0) -= nKGRAD_n(face_n_k,j)* P(faces_n[j]);

                        FLUX_MASS(k,0) = -(B_2*POT(0,0)-B_1 *POT(1,0))/(B_1+B_2) ;


					}
					else {
						FLUX_MASS(k,0)=0.0;

					for(int j = 0; j < NF; ++j)
						FLUX_MASS(k,0)+= nKGRAD(k,j)* (P(faces[j]) - P(cell)) ;
					}
					}

				//	for(int k = 0; k < NF; ++k)
				//	{
				//	FLUX_MASS(k,0)=0.0;
				////	Cell cell_l =faces[k].BackCell();
				//	Cell cell_l= cell.Neighbour(faces[k]);
				//	Cell cell_n= cell_l.Neighbour(faces[k]);
				//	for(int j = 0; j < NF; ++j)
				//	FLUX_MASS(k,0)+= nKGRAD(k,j)* (P(faces[j]) - P(cell_n)) ;
				//	}

	                    for(int k = 0; k < NF; ++k)
                    pF(k,0) = (P(faces[k]) - P(cell));//*faces[k].Area();
                    FLUX = nKGRAD*pF; //fluxes on faces
                    //Change matrix by nonlinear correction for DMP
                    ++total;
                    if( cell.GetStatus() != Element::Ghost )
                    {
                        for(int k = 0; k < NF; ++k) //loop over faces of current cell
                        { //R[P.Index(cell)] += FLUX(k,0);//faces[k].Area();
						   R[P.Index(cell)] +=  FLUX_MASS(k,0);//faces[k].Area();
						}
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

				//Solver S(Solver::INNER_ILU2);
                //Solver S(Solver::INNER_MPTILUC);
				Solver S("superlu");
                S.SetParameter("relative_tolerance", "1.0e-14");
                S.SetParameter("absolute_tolerance", "1.0e-12");
                S.SetParameter("drop_tolerance", "1.0e-1");
                S.SetParameter("reuse_tolerance", "1.0e-2");

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
                    std::cout << "Unable to solve: " << S.ReturnReason() << std::endl;
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