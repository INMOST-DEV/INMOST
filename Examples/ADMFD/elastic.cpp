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

//If true, then tag ELASTIC_TENSOR contains stiffness tensor C,
//determining strain to stress relation C*\varepsilon = \sigma, otherwise
//the tag contains compliance tensor S, determining stress to strain relation
//S*\sigma = \varepsilon. Here \varepsilon = \frac{G+G^T}{2},
//with G - gradient of displacement matrix
bool inverse_tensor = true;

//indices for tensor
const int IC11 = 0;
const int IC12 = 1;
const int IC13 = 2;
const int IC14 = 3;
const int IC15 = 4;
const int IC16 = 5;
const int IC22 = 6;
const int IC23 = 7;
const int IC24 = 8;
const int IC25 = 9;
const int IC26 = 10;
const int IC33 = 11;
const int IC34 = 12;
const int IC35 = 13;
const int IC36 = 14;
const int IC44 = 15;
const int IC45 = 16;
const int IC46 = 17;
const int IC55 = 18;
const int IC56 = 19;
const int IC66 = 20;


/*
 
 stored tensor
 C_1  C_2  C_3  C_4  C_5  C_6
      C_7  C_8  C_9  C_10 C_11
           C_12 C_13 C_14 C_15
                C_16 C_17 C_18
                     C_19 C_20
                          C_21
 
 
 voigt notation tensor
 
 c_11 c_12 c_13 c_14 c_15 c_16
 c_12 c_22 c_23 c_24 c_25 c_26
 c_13 c_23 c_33 c_34 c_35 c_36
 c_14 c_24 c_34 c_44 c_45 c_46
 c_15 c_25 c_35 c_45 c_55 c_56
 c_16 c_26 c_36 c_46 c_56 c_66
 
 co-normal tensors
 
     | K_1   \vec{n} \cdot \nabla u + K_4   \vec{n} \cdot \nabla v + K_5 \vec{n} \cdot \nabla w |
 F = | K_4^T \vec{n} \cdot \nabla u + K_2   \vec{n} \cdot \nabla v + K_6 \vec{n} \cdot \nabla w |
     | K_5^T \vec{n} \cdot \nabla u + K_6^T \vec{n} \cdot \nabla v + K_3 \vec{n} \cdot \nabla w |
 
 tensors representations
 
       | c_11 c_16 c_15 |
 K_1 = | c_16 c_66 c_56 |
       | c_15 c_56 c_55 |
 
       | c_66 c_26 c_46 |
 K_2 = | c_26 c_22 c_24 |
       | c_46 c_24 c_44 |
 
       | c_55 c_45 c_35 |
 K_3 = | c_45 c_44 c_34 |
       | c_35 c_34 c_33 |
 
       | c_16 c_66 c_56 |
 K_4 = | c_12 c_26 c_25 |
       | c_14 c_46 c_45 |
 
       | c_15 c_56 c_55 |
 K_5 = | c_14 c_46 c_45 |
       | c_13 c_36 c_35 |
 
       | c_56 c_25 c_45 |
 K_6 = | c_46 c_24 c_44 |
       | c_36 c_23 c_34 |
 
 two-point flux approximation
 
 F = T (\vec{u}_2 - \vec{u}_1)
 
 here the transimissibility matrix is
 
 T = C'_1 C'_2 (d_1 C'_2 + d_2 C'_1)^{-1}
 
      | k_1 k_4 k_5 |
 C' = | k_4 k_2 k_6 |
      | k_5 k_6 k_3 |
 
 with
 
 k_i = \vec{n} \cdot K_i \vec{n} = \vec{n} \cdot K_i^T \vec{n}
 
 */


void CTensor(const real_array & Cv, rMatrix & C)
{
	C = rMatrix::FromTensor(Cv.data(),Cv.size(),6);
}

void KTensor(const real_array & Cv, rMatrix & K)
{
	int si = 0, sj = 0;
	rMatrix C = rMatrix::FromTensor(Cv.data(),Cv.size(),6);
	
	K(0+si,0+sj) = C(0,0); K(0+si,1+sj) = C(0,5); K(0+si,2+sj) = C(0,4);
	K(1+si,0+sj) = C(0,5); K(1+si,1+sj) = C(5,5); K(1+si,2+sj) = C(4,5);
	K(2+si,0+sj) = C(0,4); K(2+si,1+sj) = C(4,5); K(2+si,2+sj) = C(4,4);
	
	si = sj = 3;
	
	K(0+si,0+sj) = C(5,5); K(0+si,1+sj) = C(1,5); K(0+si,2+sj) = C(3,5);
	K(1+si,0+sj) = C(1,5); K(1+si,1+sj) = C(1,1); K(1+si,2+sj) = C(1,3);
	K(2+si,0+sj) = C(3,5); K(2+si,1+sj) = C(1,3); K(2+si,2+sj) = C(3,3);
	
	si = sj = 6;
	
	K(0+si,0+sj) = C(4,4); K(0+si,1+sj) = C(3,4); K(0+si,2+sj) = C(2,4);
	K(1+si,0+sj) = C(3,4); K(1+si,1+sj) = C(3,3); K(1+si,2+sj) = C(2,3);
	K(2+si,0+sj) = C(2,4); K(2+si,1+sj) = C(2,3); K(2+si,2+sj) = C(2,2);
	
	si = 0, sj = 3;
	
	K(0+si,0+sj) = C(0,5); K(0+si,1+sj) = C(0,1); K(0+si,2+sj) = C(0,3);
	K(1+si,0+sj) = C(5,5); K(1+si,1+sj) = C(1,5); K(1+si,2+sj) = C(3,5);
	K(2+si,0+sj) = C(4,5); K(2+si,1+sj) = C(1,4); K(2+si,2+sj) = C(3,4);
	
	si = 3, sj = 0;
	
	K(0+si,0+sj) = C(0,5); K(0+si,1+sj) = C(5,5); K(0+si,2+sj) = C(4,5);
	K(1+si,0+sj) = C(0,1); K(1+si,1+sj) = C(1,5); K(1+si,2+sj) = C(1,4);
	K(2+si,0+sj) = C(0,3); K(2+si,1+sj) = C(3,5); K(2+si,2+sj) = C(3,4);
	
	si = 0, sj = 6;
	
	K(0+si,0+sj) = C(0,4); K(0+si,1+sj) = C(0,3); K(0+si,2+sj) = C(0,2);
	K(1+si,0+sj) = C(4,5); K(1+si,1+sj) = C(3,5); K(1+si,2+sj) = C(2,5);
	K(2+si,0+sj) = C(4,4); K(2+si,1+sj) = C(3,4); K(2+si,2+sj) = C(2,4);
	
	si = 6, sj = 0;
	
	K(0+si,0+sj) = C(0,4); K(0+si,1+sj) = C(4,5); K(0+si,2+sj) = C(4,4);
	K(1+si,0+sj) = C(0,3); K(1+si,1+sj) = C(3,5); K(1+si,2+sj) = C(3,4);
	K(2+si,0+sj) = C(0,2); K(2+si,1+sj) = C(2,5); K(2+si,2+sj) = C(2,4);
	
	si = 3, sj = 6;
	
	K(0+si,0+sj) = C(4,5); K(0+si,1+sj) = C(3,5); K(0+si,2+sj) = C(2,5);
	K(1+si,0+sj) = C(1,4); K(1+si,1+sj) = C(1,3); K(1+si,2+sj) = C(1,2);
	K(2+si,0+sj) = C(3,4); K(2+si,1+sj) = C(3,3); K(2+si,2+sj) = C(2,3);
	
	si = 6, sj = 3;
	
	K(0+si,0+sj) = C(4,5); K(0+si,1+sj) = C(1,4); K(0+si,2+sj) = C(3,4);
	K(1+si,0+sj) = C(3,5); K(1+si,1+sj) = C(1,3); K(1+si,2+sj) = C(3,3);
	K(2+si,0+sj) = C(2,5); K(2+si,1+sj) = C(1,2); K(2+si,2+sj) = C(2,3);
}

void GetBC(const real_array & bc, const rMatrix & n, rMatrix & Ra, rMatrix & Rb, rMatrix & r)
{
	const rMatrix I = rMatrix::Unit(3);
	double alpha_perp, alpha_parallel, beta_perp, beta_parallel;

	if( bc.size() == 6 )
	{
		double alpha, beta, proj;
		alpha = bc[0];
		beta = bc[1];
		proj = bc[2];
		r(0,0) = bc[3];
		r(1,0) = bc[4];
		r(2,0) = bc[5];
		if( proj )
		{
			alpha_perp = alpha;
			beta_parallel = beta;
			alpha_parallel = beta_perp = 0;
		}
		else
		{
			alpha_perp = alpha_parallel = alpha;
			beta_perp = beta_parallel = beta;
		}
	}
	else if( bc.size() == 7 )
	{
		alpha_perp = bc[0];
		beta_perp = bc[1];
		alpha_parallel = bc[2];
		beta_parallel = bc[3];
		r(0,0) = bc[4];
		r(1,0) = bc[5];
		r(2,0) = bc[6];
	}
	Ra = alpha_parallel*I + (alpha_perp-alpha_parallel)*n*n.Transpose();
	Rb = beta_parallel*I + (beta_perp-beta_parallel)*n*n.Transpose();
}

void PrintSV(const rMatrix & A)
{
	rMatrix U,S,V;
	A.SVD(U,S,V);
	std::cout << "singular values:";
	int cnt = 0;
	for(int k = 0; k < A.Rows(); ++k) if( fabs(S(k,k)) > 1.0e-10 )
	{
		std::cout << " " << S(k,k);
		cnt++;
	}
	else
		std::cout << " 0";
	std::cout << " count: " << cnt << "/" << A.Rows() << std::endl;
}

void PrintRS(const rMatrix & A)
{
	double sum;
	std::cout << "row sum:";
	for(int i = 0; i < A.Rows(); ++i)
	{
		sum = 0;
		for(int j = 0; j < A.Cols(); ++j)
			sum+= A(i,j);
		std::cout << " " << sum;
	}
	std::cout << std::endl;
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
        TagRealArray tag_UVW;// Displacement
        TagRealArray tag_C;  // Elasticity tensor in Voigt format
        TagRealArray tag_F;  // Forcing term
        TagRealArray tag_BC; // Boundary conditions
        TagRealArray tag_W;  // Gradient matrix
		TagRealArray tag_FLUX; // Flux (for error)

        if( m->GetProcessorsNumber() > 1 ) //skip for one processor job
        { // Exchange ghost cells
            ttt = Timer();
            m->ExchangeGhost(1,FACE); // Produce layer of ghost cells
            BARRIER
            if( m->GetProcessorRank() == 0 ) std::cout << "Exchange ghost: " << Timer()-ttt << std::endl;
        }
		
		std::cout << "Cells " << m->TotalNumberOf(CELL) << std::endl;
		std::cout << "Faces " << m->TotalNumberOf(FACE) << std::endl;
		std::cout << "Edges " << m->TotalNumberOf(EDGE) << std::endl;
		std::cout << "Nodes " << m->TotalNumberOf(NODE) << std::endl;
		
		const double vB[] =
		{
			1,0,0,0,0,0,
			0,0,0,0,0,1,
			0,0,0,0,1,0,
			0,0,0,0,0,1,
			0,1,0,0,0,0,
			0,0,0,1,0,0,
			0,0,0,0,1,0,
			0,0,0,1,0,0,
			0,0,1,0,0,0
		};
		const rMatrix B(vB,9,6);
		const rMatrix I = rMatrix::Unit(3);
		
		(B.Transpose()*B).Print();

        { //initialize data
            if( m->HaveTag("ELASTIC_TENSOR") ) // is elasticity tensor already defined on the mesh?
                tag_C = m->GetTag("ELASTIC_TENSOR"); // get the elasticity tensor
			else
			{
				std::cout << "No ELASTIC_TENSOR on mesh" << std::endl;
				return -1;
			}

			if( m->HaveTag("FORCE") ) //Is there force on the mesh?
			{
				tag_F = m->GetTag("FORCE"); //initial force
				assert(tag_F.isDefined(CELL)); //assuming it was defined on cells
			} // end of force

            if( m->HaveTag("BOUNDARY_CONDITION") ) //Is there boundary condition on the mesh?
            {
                tag_BC = m->GetTag("BOUNDARY_CONDITION");
                //initialize unknowns at boundary
            }
            tag_W = m->CreateTag("W",DATA_REAL,CELL,NONE);
			
            ttt = Timer();
            //Assemble gradient matrix W on cells
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				rMatrix N, R, L, T, K(9,9), C(6,6), U,S,V, W1, W2;
				rMatrix x(3,1), xf(3,1), n(3,1);
				double area; //area of the face
				double volume; //volume of the cell
				double dist; //distance from cell center to face
#if defined(USE_OMP)
#pragma omp for
#endif
				 for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
				 {
					 Mesh * mesh = m;
					 Cell cell = m->CellByLocalID(q);
					 ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
					 int NF = (int)faces.size(); //number of faces;
					 volume = cell->Volume(); //volume of the cell
					 cell->Centroid(x.data());
					 //get permeability for the cell
					 KTensor(tag_C[cell],K);
					 //CTensor(tag_C[cell],C);
					 if( !inverse_tensor ) K = K.Invert();
					 //K += rMatrix::Unit(9)*1.0e-6*K.FrobeniusNorm();
					 //PrintSV(K);
					 N.Resize(3*NF,9); //co-normals
					 //T.Resize(3*NF,9); //transversals
					 R.Resize(3*NF,9); //directions
					 L.Resize(3*NF,3*NF);
					 L.Zero();
					 //A.Resize(3*NF,3*NF);
					 //A.Zero();
					 for(int k = 0; k < NF; ++k) //loop over faces
					 {
						 area = faces[k].Area();
						 faces[k].Centroid(xf.data());
						 faces[k].OrientedUnitNormal(cell->self(),n.data());
						 dist = n.DotProduct(xf-x);
						 // assemble matrix of directions
						 R(3*k,3*(k+1),0,9) = I.Kronecker((xf-x).Transpose());
						 // assemble matrix of co-normals
						 
						 L(3*k,3*(k+1),3*k,3*(k+1)) = (area/dist)*I.Kronecker(n.Transpose())*K*I.Kronecker(n);
						 
						 N(3*k,3*(k+1),0,9) = area*I.Kronecker(n.Transpose());
						 
						 //std::cout << "I\otimes n^T" << std::endl;
						 //I.Kronecker(n.Transpose()).Print();
						 
						 //std::cout << "I\otimes n^T B" << std::endl;
						 
						 //(I.Kronecker(n.Transpose())*B).Print();
						 
						 //T(3*k,3*(k+1),0,9) = N(3*k,3*(k+1),0,9)*K - L(3*k,3*(k+1),3*k,3*(k+1))*I.Kronecker(n.Transpose());
						 //A(3*k,3*(k+1),3*k,3*(k+1)) = I*area;
						 //NK(3*k,3*(k+1),0,9) = area*I.Kronecker(n.Transpose());
					 } //end of loop over faces
					 tag_W[cell].resize(9*NF*NF);
					 //int ierr = -1;
					 //tag_W(cell,3*NF,3*NF) = N*(N.Transpose()*R).PseudoInvert(1.0e-9)*N.Transpose() //stability part
					 // + (rMatrix::Unit(3*NF) - R*(R.Transpose()*R).PseudoInvert(1.0e-9)*R.Transpose())*
					 //(2.0/(static_cast<real>(NF)*volume)*(N*K.PseudoInvert(1.0e-9)*N.Transpose()).Trace());
					 
					 //+  (rMatrix::Unit(3*NF) - R*(N.Transpose()*R).PseudoInvert(1.0e-12)*N.Transpose())*
					  //   (2.0/(static_cast<real>(NF)*volume)*(N*K.Invert()*N.Transpose()).Trace());
					 
					 //tag_W(cell,3*NF,3*NF) = N*(N.Transpose()*R).PseudoInvert(1.0e-9)*N.Transpose() + (L - L*R*((L*R).Transpose()*R).PseudoInvert(1.0e-9)*(L*R).Transpose());
					 
					 
					 //tag_W(cell,3*NF,3*NF) = L + (N*K-L*R)*((N*K-L*R).Transpose()*R).PseudoInvert(1.0e-9)*(N*K-L*R).Transpose() + (L - L*R*((L*R).Transpose()*R).PseudoInvert(1.0e-9)*(L*R).Transpose());
					 
					 //tag_W(cell,3*NF,3*NF) = (N*K)*((N*K).Transpose()*R).PseudoInvert(1.0e-9)*(N*K).Transpose() + (L - L*R*((L*R).Transpose()*R).PseudoInvert(1.0e-9)*(L*R).Transpose());
					 
					 //W1.Resize(3*NF,3*NF);
					 //W2.Resize(3*NF,3*NF);
					 
					 W1 = (N*K)*((N*K).Transpose()*R).PseudoInvert(1.0e-7)*(N*K).Transpose();
					 W2 = L - (L*R)*((L*R).Transpose()*R).PseudoInvert(1.0e-7)*(L*R).Transpose();
					 
					 tag_W(cell,3*NF,3*NF) = W1+W2;
					 //W2 = W1.Trace()*(rMatrix::Unit(3*NF) - R*(R.Transpose()*R).PseudoInvert(1.0e-7)*R.Transpose());
					 //W2 = W1.Trace()*(rMatrix::Unit(3*NF) - R*((N*K).Transpose()*R).PseudoInvert(1.0e-7)*(N*K).Transpose());
					 /*
					 std::cout << "L " << std::endl;
					 L.Print();
					 
					 std::cout << "N*K*N^T" << std::endl;
					 (N*K*N.Transpose()).Print();
					 
					 std::cout << "N*R^T" << std::endl;
					 (N*R.Transpose()).Print();
					 */
					 //L = N*K*N.Transpose()*(N*R.Transpose()).PseudoInvert(1.0e-9);
					 //std::cout << "L2 " << std::endl;
					 // L.Print();
					 
					 //PrintSV(W1);
					 //PrintSV(W2);
					 //PrintSV(W1+W2);
					 
					 //tag_W(cell,3*NF,3*NF) = W1+W2;
					 
					 //std::cout << "W1(" << W1.Rows() << "," << W1.Cols() << ")" << std::endl;
					 //W1.Print();
					 
					 //std::cout << "W2(" << W2.Rows() << "," << W2.Cols() << ")" << std::endl;
					 //W2.Print();
					 /*
					 tag_W(cell,3*NF,3*NF) =
					 L
					 //+ (N*K - L*R)*((N*K-L*R).Transpose()*R).PseudoInvert(1.0e-12)*(N*K-L*R).Transpose()
					 + (N*K)*((N*K).Transpose()*R).PseudoInvert(1.0e-7)*(N*K).Transpose()
					 //+ L.Trace()*(rMatrix::Unit(3*NF) - R*((R).Transpose()*R).PseudoInvert(1.0e-14)*(R).Transpose())
					 - (L*R)*((L*R).Transpose()*R).PseudoInvert(1.0e-7)*(L*R).Transpose()
					 ;
					 */
					 //(tag_W(cell,3*NF,3*NF) - W1-W2).Print();
					 //tag_W(cell,3*NF,3*NF) = W1+W2;
					 
					 
					 //std::cout << " (NK)^T*R ";
					 //PrintSV((N*K).Transpose()*R);
					 
					 //std::cout << " (LR)^T*R ";
					 //PrintSV((L*R).Transpose()*R);
					 
					 
					 //std::cout << " vol mat: " << std::endl;
					 //((N).Transpose()*R).PseudoInvert(1.0e-12).Print();
					 
					 //std::cout << " grad mat: ";// << std::endl;
					 //PrintSV(((N).Transpose()*R).PseudoInvert(1.0e-12)*(N).Transpose());
					 //(((N).Transpose()*R).PseudoInvert(1.0e-12)*(N).Transpose()).Print();
					 /*
					 std::cout << " W1 ";
					 PrintSV(L);
					 
					 L.Print();
					 
					 std::cout << " W2 ";
					 PrintSV((N*K - L*R)*((N*K).Transpose()*R).PseudoInvert(1.0e-12)*(N*K).Transpose());
					 
					 ((N*K - L*R)*((N*K).Transpose()*R).PseudoInvert(1.0e-12)*(N*K).Transpose()).Print();
					 
					 std::cout << " W ";
					 PrintSV(tag_W(cell,3*NF,3*NF));
					 tag_W(cell,3*NF,3*NF).Print();
					 */
					 //Convergent on tetra
					 //tag_W(cell,3*NF,3*NF) =
					 //L
					 //+ (N*K-L*R)*((N*K).Transpose()*R).PseudoInvert(1.0e-14)*(N*K).Transpose()
					 //+ (N*K-L*R)*((N*K-L*R).Transpose()*R).PseudoInvert(1.0e-14)*(N*K-L*R).Transpose()
					 //+ (L - L*R*((L*R).Transpose()*R).PseudoInvert(1.0e-14)*(L*R).Transpose())
					 //+ L.Trace()*(rMatrix::Unit(3*NF) - R*((R).Transpose()*R).PseudoInvert(1.0e-14)*(R).Transpose())
					 ;
					 
					 //tag_W(cell,3*NF,3*NF) = (N*K)*((N*K).Transpose()*R).PseudoInvert(1.0e-14)*(N*K).Transpose() + (L - L*R*((L*R).Transpose()*R).PseudoInvert(1.0e-14)*(L*R).Transpose());
					 
					 
					 //tag_W(cell,3*NF,3*NF) = N*K*((N*K).Transpose()*R).PseudoInvert(1.0e-9)*(N*K).Transpose() + (L - L*R*((L*R).Transpose()*R).PseudoInvert(1.0e-9)*(L*R).Transpose());
					 
					 //tag_W(cell,3*NF,3*NF) = L + (N*K-L*R)*(R.Transpose()*R).PseudoInvert(1.0e-9)*R.Transpose();
					 //tag_W(cell,3*NF,3*NF) = (N*K)*(N.Transpose()*R).PseudoInvert(1.0e-9)*N.Transpose();
					 
					 //std::cout << "K*V" << std::endl;
					 //(K).Print();
					 //std::cout << "N^T*R" << std::endl;
					 //(N.Transpose()*R/volume).Print();
					 
					 //std::cout << "R^T*R" << std::endl;
					 //(R.Transpose()*R).Print();
					 //std::cout << "(R^T*R)^{-1}" << std::endl;
					 //(R.Transpose()*R).PseudoInvert(1.0e-12).Print();
					 //std::cout << "W:" <<std::endl;
					 //tag_W(cell,3*NF,3*NF).Print();
					 //PrintRS(tag_W(cell,3*NF,3*NF));
					 if(false)if(!K.isSymmetric() )
					 {
						 std::cout << "K nonsymmetric: " << std::endl;
						 (K-K.Transpose()).Print();
					 }
					 
					 if( true ) if( !tag_W(cell,3*NF,3*NF).isSymmetric(1.0e-3) )
					 {
						 std::cout << "W nonsymmetric: " << std::endl;
						 (tag_W(cell,3*NF,3*NF)-tag_W(cell,3*NF,3*NF).Transpose()).Print();
					 }
					 
					 if( false ) if( (N*K - tag_W(cell,3*NF,3*NF)*R).FrobeniusNorm() > 1.0e-3 )
					 {
						 std::cout << "error: " << std::endl;
						 (N*K - tag_W(cell,3*NF,3*NF)*R).Print();
					 }
					 
					 //std::cout << "L ";
					 //PrintSV(L);
					 
					 //std::cout << "W2 ";
					 //PrintSV((N*K-L*R)*((N*K-L*R).Transpose()*R).PseudoInvert(1.0e-9)*(N*K-L*R).Transpose());
					 
					 //std::cout << "R ";
					 //PrintSV(rMatrix::Unit(3*NF) - R*(R.Transpose()*R).PseudoInvert(1.0e-9)*R.Transpose());
					 
					 //std::cout << "consistency ";
					 //PrintSV(NK*(NK.Transpose()*R).PseudoInvert(1.0e-12,&ierr)*NK.Transpose());
					 //std::cout << "stability   ";
					 //PrintSV((1.5/(static_cast<real>(NF)*volume)*(NK*K.Invert()*NK.Transpose()).Trace())*(rMatrix::Unit(3*NF) - R*(R.Transpose()*R).PseudoInvert(1.0e-12)*R.Transpose()));
					 
					 
					 //std::cout << "K ";
					 //PrintSV(K);
					 if( false )
					 {
						 std::cout << "total       ";
						 PrintSV(tag_W(cell,3*NF,3*NF));
					 }
					 
					 
					 //std::cout << "Check:" << std::endl;
					 //(NK - tag_W(cell,3*NF,3*NF)*R).Print();
					 //std::cout << "R^T*R:" << std::endl;
					 //(R.Transpose()*R).Print();
					 //std::cout << "NK^T*R:" << std::endl;
					 //(NK.Transpose()*R).Print();
					 /*
					 if( ierr )
					 {
						 std::cout << "K:" << std::endl;
						 K.Print();
						 std::cout << "NK:" << std::endl;
						 NK.Print();
						 std::cout << "R:" << std::endl;
						 R.Print();
						 std::cout << "NK^T*R:" << std::endl;
						 (NK.Transpose()*R).Print();
						 std::cout << "sym? " << tag_W(cell,3*NF,3*NF).isSymmetric() << std::endl;
						 std::cout << "ierr " << ierr << std::endl;
					 }
					  */
				 } //end of loop over cells
			}
            std::cout << "Construct W matrix: " << Timer() - ttt << std::endl;
			
			if( m->HaveTag("UVW") ) //Is there a displacement on the mesh?
				tag_UVW = m->GetTag("UVW"); //Get the pressure
			
			if( !tag_UVW.isValid() || !tag_UVW.isDefined(CELL) ) // Pressure was not initialized or was not defined on nodes
			{
				//srand(1); // Randomization
				tag_UVW = m->CreateTag("UVW",DATA_REAL,CELL|FACE,NONE,3); // Create a new tag for the displacement
				for(Mesh::iteratorElement e = m->BeginElement(CELL|FACE); e != m->EndElement(); ++e) //Loop over mesh cells
					tag_UVW(*e,3,1).Zero();//(rand()*1.0)/(RAND_MAX*1.0); // Prescribe random value in [0,1]
			}
			
			if( !tag_UVW.isDefined(FACE) )
			{
				tag_UVW = m->CreateTag("UVW",DATA_REAL,FACE,NONE,3);
				for(Mesh::iteratorElement e = m->BeginElement(FACE); e != m->EndElement(); ++e) //Loop over mesh cells
					tag_UVW(*e,3,1).Zero();//(rand()*1.0)/(RAND_MAX*1.0); // Prescribe random value in [0,1]
			}
			m->ExchangeData(tag_UVW,CELL|FACE,0); //Synchronize initial solution with boundary unknowns
			
			
			tag_FLUX = m->CreateTag("FLUX",DATA_REAL,FACE,NONE,3);
        } //end of initialize data
		

        std::cout << "Initialization done" << std::endl;


        integer nit = 0;
        ttt = Timer();

        { //Main loop for problem solution
            Automatizator aut; // declare class to help manage unknowns
            Automatizator::MakeCurrent(&aut);
			BlockEntry UVW(CELL|FACE);
			UVW.AddTag(tag_UVW);
			aut.RegisterEntry(UVW);
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
					{
                        Text.SetAnnotation(UVW.Index(cell,0),"Divergence, U");
						Text.SetAnnotation(UVW.Index(cell,1),"Divergence, V");
						Text.SetAnnotation(UVW.Index(cell,2),"Divergence, W");
					}
                }
                for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
                {
                    Face face = m->FaceByLocalID(q);
                    if( face.GetStatus() != Element::Ghost )
                    {
                        if( tag_BC.isValid() && face.HaveData(tag_BC) )
						{
                            Text.SetAnnotation(UVW.Index(face,0),"Boundary condition, U");
							Text.SetAnnotation(UVW.Index(face,1),"Boundary condition, V");
							Text.SetAnnotation(UVW.Index(face,2),"Boundary condition, W");
						}
                        else
						{
                            Text.SetAnnotation(UVW.Index(face,0),"Flux continuity, U");
							Text.SetAnnotation(UVW.Index(face,1),"Flux continuity, V");
							Text.SetAnnotation(UVW.Index(face,2),"Flux continuity, W");
						}
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
					
					vMatrix U,T; //vector of pressure differences and fluxes on faces
					rMatrix Ra(3,3), Rb(3,3), r(3,1), n(3,1);
#if defined(USE_OMP)
#pragma omp for
#endif
					for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
					{
						Face face = m->FaceByLocalID(q);
						tag_FLUX(face,3,1).Zero();
					}
#if defined(USE_OMP)
#pragma omp for
#endif
					for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) ) //loop over cells
					{
						Cell cell = m->CellByLocalID(q);
						ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
						int NF = (int)faces.size();
						
						raMatrix W = raMatrixMake(tag_W[cell].data(),3*NF,3*NF); //Matrix for gradient
						
						U.Resize(3*NF,1);
						T.Resize(3*NF,1);
						
						for(int k = 0; k < NF; ++k)
							U(3*k,3*(k+1),0,1) = (UVW[faces[k]] - UVW[cell]);
						T = W*U; //fluxes on faces
						if( cell.GetStatus() != Element::Ghost )
						{
							for(int k = 0; k < NF; ++k) //loop over faces of current cell
								R[UVW.Index(cell)] += T(3*k,3*(k+1),0,1);
						}
						for(int k = 0; k < NF; ++k) //loop over faces of current cell
						{
							if( faces[k].GetStatus() == Element::Ghost ) continue;
							int index = UVW.Index(faces[k],0);
							Locks.Lock(index);
							if( tag_BC.isValid() && faces[k].HaveData(tag_BC) )
							{
								faces[k]->UnitNormal(n.data());
								GetBC(faces[k].RealArray(tag_BC),n,Ra,Rb,r);
								if( Rb.FrobeniusNorm() ) Rb.Print();
								R[UVW.Index(faces[k])] -= Ra*UVW[faces[k]] + Rb*T(3*k,3*(k+1),0,1)/faces[k].Area() - r;
								
								
								tag_FLUX(faces[k],3,1) += rMatrix(T(3*k,3*(k+1),0,1))/faces[k].Area();
							}
							else
							{
								R[UVW.Index(faces[k])] -= T(3*k,3*(k+1),0,1);
								
								tag_FLUX(faces[k],3,1) += 0.5*rMatrix(T(3*k,3*(k+1),0,1))/faces[k].Area()*(faces[k].FaceOrientedOutside(cell)?1:-1);
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
							 if( cell->HaveData(tag_F) ) R[UVW.Index(cell)] += tag_F(cell,3,1)*cell->Volume();
						 }
					 }
				}
				std::cout << "assembled in " << Timer() - tttt << "\t\t\t" << std::endl;

                std::cout << "Nonlinear residual: " << R.Norm() << "\t\t" << std::endl;

                if( R.Norm() < 1.0e-4 ) break;
				tttt = Timer();

				Solver S(Solver::INNER_MPTILU2);
                //Solver S(Solver::INNER_MPTILUC);
				//Solver S("superlu");
                S.SetParameter("relative_tolerance", "1.0e-14");
                S.SetParameter("absolute_tolerance", "1.0e-12");
                S.SetParameter("drop_tolerance", "1.0e-2");
                S.SetParameter("reuse_tolerance", "1.0e-4");

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
                        tag_UVW(cell,3,1) -= Update[UVW.Index(cell)];
                    }
#if defined(USE_OMP)
#pragma omp parallel for
#endif
                    for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
                    {
                        Face face = m->FaceByLocalID(q);
						if( face->GetStatus() == Element::Ghost ) continue;
                        tag_UVW(face,3,1) -= Update[UVW.Index(face)];
                    }
                    m->ExchangeData(tag_UVW, CELL|FACE, 0);
					
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
				std::cout << "solved in " << Timer() - tttt << "\t\t\t" << std::endl;
                ++nit;
            } while( R.Norm() > 1.0e-4 && nit < 10); //check the residual norm
        }
        std::cout << "Solved problem in " << Timer() - ttt << " seconds with " << nit << " iterations " << std::endl;

        if( m->HaveTag("REFERENCE_SOLUTION") )
        {
            TagReal tag_E = m->CreateTag("ERROR",DATA_REAL,CELL,NONE,1);
            TagRealArray tag_R = m->GetTag("REFERENCE_SOLUTION");
            real C, L2, volume;
            C = L2 = volume = 0.0;
            for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
            {
                Cell cell = m->CellByLocalID(q);
                real err = (tag_UVW(cell,3,1) - tag_R(cell,3,1)).FrobeniusNorm();
                real vol = cell->Volume();
                if( C < fabs(err) ) C = fabs(err);
                L2 += err*err*vol;
                volume += vol;
                tag_E[cell] = err;
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
                    real err = (tag_UVW(face,3,1) - tag_R(face,3,1)).FrobeniusNorm();
                    real vol = (face->BackCell()->Volume() + (face->FrontCell().isValid() ? face->FrontCell()->Volume() : 0))*0.5;
                    if( C < fabs(err) ) C = fabs(err);
                    L2 += err*err*vol;
                    volume += vol;
                    tag_E[face] = err;
                }
				C = m->AggregateMax(C);
				L2 = m->Integrate(L2);
				volume = m->Integrate(volume);
                L2 = sqrt(L2/volume);
                if( m->GetProcessorRank() == 0 ) std::cout << "Error on faces, C-norm " << C << " L2-norm " << L2 << std::endl;
            }
            else std::cout << "Reference solution was not defined on faces" << std::endl;
        }
		
		if( m->HaveTag("REFERENCE_FLUX") )
		{
			TagReal tag_E = m->CreateTag("ERROR_FLUX",DATA_REAL,FACE,NONE,1);
			TagRealArray tag_R = m->GetTag("REFERENCE_FLUX");
			real C, L2, volume;
			C = L2 = volume = 0.0;
			for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
			{
				Face face = m->FaceByLocalID(q);
				real err = (tag_FLUX(face,3,1) - tag_R(face,3,1)).FrobeniusNorm();
				real vol = face->BackCell()->Volume()/(real)face->BackCell().nbAdjElements(FACE);
				if( face->FrontCell().isValid() )
					vol += face->FrontCell()->Volume()/(real)face->FrontCell().nbAdjElements(FACE);
				if( C < fabs(err) ) C = fabs(err);
				L2 += err*err*vol;
				volume += vol;
				tag_E[face] = err;
			}
			C = m->AggregateMax(C);
			L2 = m->Integrate(L2);
			volume = m->Integrate(volume);
			L2 = sqrt(L2/volume);
			if( m->GetProcessorRank() == 0 ) std::cout << "Error of fluxes, C-norm " << C << " L2-norm " << L2 << std::endl;
		}

        if( m->GetProcessorsNumber() == 1 )
            m->Save("out.vtk");
        else
            m->Save("out.pvtk");
		m->Save("out.pmf");

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
