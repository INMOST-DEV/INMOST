#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../../Source/Misc/utils.h"


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
bool reference_solution = false;
bool print_matrix = false;


template<typename T>
static void SVD2Eigen(const Matrix<T> & U, Matrix<T> & S, Matrix<T> & V)
{
	for (unsigned i = 0; i < V.Cols(); ++i)
	{
		Storage::real dot = 0.0;
		for (unsigned j = 0; j < V.Rows(); ++j)
			dot += get_value(U(j, i))*get_value(V(j, i));
		if (dot < 0.0)
		{
			S(i, i) *= -1;
			for (unsigned j = 0; j < V.Rows(); ++j)
				V(j, i) *= -1;
		}
	}
	//check
	//if ((U - V).FrobeniusNorm() > 1.0e-8)
	//	(U - V).Print();
}


class BlockRow
{
	std::map< unsigned, rMatrix > entries;
	unsigned height, width;
	rMatrix & GetAdd(unsigned col)
	{
		std::map< unsigned, rMatrix >::iterator s = entries.find(col/width);
		if( s == entries.end() )
			return (entries[col/width] = rMatrix(height,width,0.0));
		else
			return s->second;
	}
public:
	BlockRow(unsigned h, unsigned w) : height(h), width(w) {}
	rMatrix & operator [] (unsigned col) { return GetAdd(col); }
	Storage::real & operator () (unsigned i, unsigned j) { return GetAdd(j)(i%height,j%width);}
	Storage::real operator () (unsigned i, unsigned j) const
	{
		std::map< unsigned, rMatrix >::const_iterator s = entries.find(j/width);
		if( s != entries.end() )
			return s->second(i%height,j%width);
		throw -1;
	}
	std::map< unsigned, rMatrix >::iterator Begin() {return entries.begin();}
	std::map< unsigned, rMatrix >::iterator End() {return entries.end();}
};

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
	Storage::real 
		alpha_perp = 0, 
		alpha_parallel = 0, 
		beta_perp = 1, 
		beta_parallel = 1;

	if( bc.size() == 6 )
	{
		Storage::real alpha, beta, proj;
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
	for(unsigned k = 0; k < A.Cols(); ++k) if( fabs(S(k,k)) > 1.0e-10 )
	{
		std::cout << " " << S(k,k);
		cnt++;
	}
	else
		std::cout << " 0";
	std::cout << " count: " << cnt << "/" << A.Cols() << std::endl;
}

void PrintRS(const rMatrix & A)
{
	Storage::real sum;
	std::cout << "row sum:";
	for(unsigned i = 0; i < A.Rows(); ++i)
	{
		sum = 0;
		for(unsigned j = 0; j < A.Cols(); ++j)
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
            table[BARYCENTER]    = CELL | FACE; //Compute averaged center of mass
            table[NORMAL]      = FACE;        //Compute normals
            table[ORIENTATION] = FACE;        //Check and fix normal orientation
            table[MEASURE]     = CELL | FACE; //Compute volumes and areas
            //table[BARYCENTER]  = CELL | FACE; //Compute volumetric center of mass
			m->RemoveGeometricData(table);
            m->PrepareGeometricData(table); //Ask to precompute the data
            BARRIER
            if( m->GetProcessorRank() == 0 ) std::cout << "Prepare geometric data: " << Timer()-ttt << std::endl;
        }

        // data tags for
        TagRealArray tag_UVW;// Displacement
        TagRealArray tag_S;  // Stress (for error)
        TagRealArray tag_C;  // Elasticity tensor in Voigt format
        TagRealArray tag_F;  // Forcing term
        TagRealArray tag_BC; // Boundary conditions
        TagRealArray tag_W;  // Gradient matrix
        TagRealArray tag_Ws; // Matrix to reconstruct stress from fluxes
		TagRealArray tag_FLUX; // Flux (for error)
		TagRealArray tag_f_coefs, tag_i_coefs, tag_f_rhs, tag_i_rhs;
		TagReferenceArray tag_f_elems, tag_i_elems;
		TagIntegerArray tag_Row;

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
		
		const Storage::real vB[] =
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
		const Storage::real viBtB[] =
		{
			1,0,0,0.0,0.0,0.0,
			0,1,0,0.0,0.0,0.0,
			0,0,1,0.0,0.0,0.0,
			0,0,0,0.5,0.0,0.0,
			0,0,0,0.0,0.5,0.0,
			0,0,0,0.0,0.0,0.5
		};
		// |  0  -z   y | u   | w_y - v_z |
		// |  z   0  -x | v = | u_z - w_x |
		// | -y   x   0 | w   | v_x - u_y |
		//
		// 0  0  0  0  0 -z  0 +y  0
		// 0  0 +z  0  0  0 -x  0  0
		// 0 -y  0 +x  0  0  0  0  0
		//
		// 0  0  0  0  0  0  0  0  0
		// 0  0  0  0  0  0  0  1  0
		// 0  0  0  0  0 -1  0  0  0
		// 0  0  0  0  0  0 -1  0  0
		// 0  0  0  0  0  0  0  0  0
		// 0  0  1  0  0  0  0  0  0
		// 0  0  0  1  0  0  0  0  0
		// 0 -1  0  0  0  0  0  0  0
		// 0  0  0  0  0  0  0  0  0
		//
		const Storage::real vCurl[] =
		{
			0,  0,  0,  0,  0,  0,  0,  0,  0,
			0,  0,  0,  0,  0,  0,  0,  1,  0,
			0,  0,  0,  0,  0, -1,  0,  0,  0,
			0,  0,  0,  0,  0,  0, -1,  0,  0,
			0,  0,  0,  0,  0,  0,  0,  0,  0,
			0,  0,  1,  0,  0,  0,  0,  0,  0,
			0,  0,  0,  1,  0,  0,  0,  0,  0,
			0, -1,  0,  0,  0,  0,  0,  0,  0,
			0,  0,  0,  0,  0,  0,  0,  0,  0
		};
		const rMatrix B(vB,9,6);
		const rMatrix iBtB(viBtB,6,6);
		const rMatrix Curl(vCurl,9,9);
		const rMatrix I = rMatrix::Unit(3);
		const rMatrix I9 = rMatrix::Unit(9);
		
		//PrintSV(B);
		std::cout << "B^T*B" << std::endl;
		(B.Transpose()*B).Print();
		
		std::cout << "B*iBtB*B^T" << std::endl;
		(B*iBtB*B.Transpose()).Print();
		
		std::cout << "B^T*B*iBtB" << std::endl;
		(B.Transpose()*B*iBtB).Print();
		
		std::cout << "B*B^T" << std::endl;
		(B*B.Transpose()).Print();

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
			tag_i_elems = m->CreateTag("i_elems",DATA_REFERENCE,FACE,NONE);
			tag_f_elems = m->CreateTag("f_elems",DATA_REFERENCE,FACE,NONE);
			tag_i_coefs = m->CreateTag("i_coefs",DATA_REAL,FACE,NONE);
			tag_f_coefs = m->CreateTag("f_coefs",DATA_REAL,FACE,NONE);
			tag_i_rhs = m->CreateTag("i_rhs",DATA_REAL,FACE,NONE,3);
			tag_f_rhs = m->CreateTag("f_rhs",DATA_REAL,FACE,NONE,3);
            tag_W = m->CreateTag("W",DATA_REAL,CELL,NONE);
            tag_Ws = m->CreateTag("Ws",DATA_REAL,CELL,NONE);
			
            ttt = Timer();
            //Assemble gradient matrix W on cells
			
			tag_Row = m->CreateTag("row_index",DATA_INTEGER,FACE,NONE,2);
#if defined(USE_OMP)
//#pragma omp parallel
#endif
			{
				rMatrix N, R, L, M(9,9), T, K(9,9), iK(9,9), C(6,6), Q, W1, W2, W3, W3s, U,S,V, w, u, v;
				rMatrix x(3,1), xf(3,1), n(3,1);
				Storage::real area; //area of the face
				//double volume; //volume of the cell
				Storage::real dist; //distance from cell center to face
#if defined(USE_OMP)
#pragma omp for
#endif
				for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
				{
					Cell cell = m->CellByLocalID(q);
					ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
					int NF = (int)faces.size(); //number of faces;
					//~ volume = cell->Volume(); //volume of the cell
					cell->Barycenter(x.data());;
					//get permeability for the cell
					KTensor(tag_C[cell],K);
					CTensor(tag_C[cell],C);
					iK = K.PseudoInvert(1.0e-11);
					if( !inverse_tensor ) K = K.Invert();
					//K += rMatrix::Unit(9)*1.0e-6*K.FrobeniusNorm();
					//PrintSV(K);
					N.Resize(3*NF,9); //co-normals
					//NQ.Resize(3*NF,9);
					//T.Resize(3*NF,9); //transversals
					R.Resize(3*NF,9); //directions
					L.Resize(3*NF,3*NF);
					//M.Resize(3*NF,3*NF);
					L.Zero();
					
					for(int k = 0; k < NF; ++k) //loop over faces
					{
						area = faces[k].Area();
						faces[k].Barycenter(xf.data());
						faces[k].OrientedUnitNormal(cell->self(),n.data());
						dist = n.DotProduct(xf-x);
						// assemble matrix of directions
						R(3*k,3*(k+1),0,9) = I.Kronecker((xf-x).Transpose());
						// assemble matrix of co-normals
						
						L(3*k,3*(k+1),3*k,3*(k+1)) = (area/dist)*(I.Kronecker(n.Transpose())*K*I.Kronecker(n));
						
						N(3*k,3*(k+1),0,9) = area*I.Kronecker(n.Transpose());
#if defined(USE_OMP)
#pragma omp critical
#endif
						{
							if( faces[k].BackCell() == cell )
								tag_Row[faces[k]][0] = k; //detect my row
							else
								tag_Row[faces[k]][1] = k; //detect my row
						}
					} //end of loop over faces
					//R += N*Curl;
					tag_Ws[cell].resize(6*3*NF);
					tag_Ws(cell,6,3*NF) = ((R*B*iBtB).Transpose()*N*B).Invert()*(R*B*iBtB).Transpose();
					tag_W[cell].resize(9*NF*NF);
					
					if( true )
					{
						Storage::real alpha = 0;
						Storage::real beta = alpha;
						
						//M = B*iBtB*B.Transpose();
						//R = R*(I9 + M)*0.5;
						W1 = (N*K+alpha*L*R)*((N*K+alpha*L*R).Transpose()*R).PseudoInvert(1.0e-11)*(N*K+alpha*L*R).Transpose();
						W2 = L - (1+beta)*(L*R)*((L*R).Transpose()*R).PseudoInvert(1.0e-11)*(L*R).Transpose();
					}
					else
					{
						M = B*iBtB*B.Transpose();
						W1 = (N*B*C)*((N*B*C).Transpose()*R*B*iBtB).Invert()*(N*B*C).Transpose();
						W2 = L - (L*R*B*iBtB)*((L*R*B*iBtB).Transpose()*R*B*iBtB).Invert()*(L*R*B*iBtB).Transpose();
						
					}
					
					tag_W(cell,3*NF,3*NF) = W1+W2;// + rMatrix::Unit(3*NF)*volume;
					
					if(false)if(!K.isSymmetric() )
					{
						std::cout << "K nonsymmetric: " << std::endl;
						(K-K.Transpose()).Print();
					}
					
					//if( false )
					if( !tag_W(cell,3*NF,3*NF).isSymmetric(1.0e-5) )
					{
						std::cout << "W nonsymmetric: " << std::endl;
						(tag_W(cell,3*NF,3*NF)-tag_W(cell,3*NF,3*NF).Transpose()).Print();
						//std::cout << "(N*C)^T*R" << std::endl;
						//((N*C).Transpose()*R).Print();
					}
					
					if( false )
						if( (N*C - tag_W(cell,3*NF,3*NF)*R).FrobeniusNorm() > 1.0e-3 )
						{
							std::cout << "error: " << std::endl;
							(N*C - tag_W(cell,3*NF,3*NF)*R).Print();
						}
					
					
					
					if( false )
					{
						std::cout << "total       ";
						PrintSV(tag_W(cell,3*NF,3*NF));
					}
					
					
				} //end of loop over cells
				rMatrix C1(3,3),C2(3,3), iC12sum, C1sum(3,3), C2sum(3,3), CC(3,3), CCI(3,3), CCF(3,3), B_D(3,3), B_N(3,3), B_R(3,1);
				rMatrix W[2];
#if defined(USE_OMP)
#pragma omp for
#endif
				for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
				{
					Face f = m->FaceByLocalID(q);
					
					Storage::integer_array rows = f->IntegerArray(tag_Row);
					ElementArray<Face> faces[2];
					ElementArray<Cell> cells = f->getCells();
					//std::cout << "face " << f->LocalID() << " NC " << cells.size() << std::endl;
					for(int q = 0; q < (int)cells.size(); ++q)
					{
						Cell cell = cells[q];
						faces[q] = cell->getFaces(); //obtain faces of the cell
						W[q] = raMatrixMake(cell->RealArray(tag_W).data(),3*faces[q].size(),3*faces[q].size());
					}
					Storage::reference_array f_elems = f->ReferenceArray(tag_f_elems);
					Storage::reference_array i_elems = f->ReferenceArray(tag_i_elems);
					Storage::real_array f_coefs = f->RealArray(tag_f_coefs);
					Storage::real_array i_coefs = f->RealArray(tag_i_coefs);
					tag_f_rhs(f,3,1).Zero();
					tag_i_rhs(f,3,1).Zero();
					real aF = 1;
					//f_rhs = i_rhs = 0;
					real multi = 1, multf = -aF;//, multfbnd = -aF;
					if( cells.size() == 2 ) //internal face
					{
						C1 = W[0](rows[0]*3,rows[0]*3+3,rows[0]*3,rows[0]*3+3)*aF;
						C2 = W[1](rows[1]*3,rows[1]*3+3,rows[1]*3,rows[1]*3+3)*aF;
						iC12sum = (C1+C2).Invert();
						C1sum.Zero();
						C2sum.Zero();
						//add back cell faces and compute sum
						for(int q = 0; q < (int)faces[0].size(); ++q)
						{
							CC = W[0](rows[0]*3,rows[0]*3+3,q*3,q*3+3);
							if( CC.FrobeniusNorm() > 1.0e-12 )
							{
								C1sum += CC;
								CC = iC12sum*CC;
								if( faces[0][q] != f )
								{
									CCI = -CC*multi;
									CCF = -CC*C2*multf;
									i_coefs.insert(i_coefs.end(),CCI.data(),CCI.data()+9);
									f_coefs.insert(f_coefs.end(),CCF.data(),CCF.data()+9);
									i_elems.push_back(faces[0][q]);
									f_elems.push_back(faces[0][q]);
									
								}
							}
							
						}
						//add front cell faces and compute sum
						for(int q = 0; q < (int)faces[1].size(); ++q)
						{
							CC = W[1](rows[1]*3,rows[1]*3+3,q*3,q*3+3);
							if( CC.FrobeniusNorm() > 1.0e-12 )
							{
								C2sum += CC;
								CC = iC12sum*CC;
								if( faces[1][q] != f )
								{
									CCI = -CC*multi;
									CCF = CC*C1*multf;
									i_coefs.insert(i_coefs.end(),CCI.data(),CCI.data()+9);
									f_coefs.insert(f_coefs.end(),CCF.data(),CCF.data()+9);
									i_elems.push_back(faces[1][q]);
									f_elems.push_back(faces[1][q]);
								}
							}
						}
						//already have back cell and front cell
						raMatrixMake(f_coefs.data()+0,3,3) = iC12sum*C1sum*C2*multf; //BackCell
						assert(f_elems[0] == f->BackCell());
						raMatrixMake(f_coefs.data()+9,3,3) = -iC12sum*C1*C2sum*multf; //FrontCell
						assert(f_elems[1] == f->FrontCell());
						//i_coefs[0] = iC12sum*C1sum*multi;
						//i_coefs[1] = iC12sum*C2sum*multi;
					}
					else //Boundary face
					{
						B_D.Zero();
						B_N = I;
						B_R.Zero();
						if( tag_BC.isValid() && f.HaveData(tag_BC) )
						{
							f.UnitNormal(n.data());
							GetBC(f.RealArray(tag_BC),n,B_D,B_N,B_R);
						}
						C1 = W[0](rows[0]*3,rows[0]*3+3,rows[0]*3,rows[0]*3+3)*aF;
						C1sum.Zero();
						//TODO
						//iC12sum = (alpha*I+beta*C1).Invert();
						
						tag_f_rhs(f,3,1).Zero();
						tag_i_rhs(f,3,1) = iC12sum*B_R;
						for(int q = 0; q < (int)faces[0].size(); ++q) if( fabs(W[0](rows[0],q)) > 1.0e-12 )
						{
							CC = W[0](rows[0]*3,rows[0]*3+3,q*3,q*3+3);
							if( CC.FrobeniusNorm() > 1.0e-12 )
							{
								//TODO
								//f_elems.push_back(faces[0][q]);
								//f_coefs.push_back(-CC*multfbnd);
								if( faces[0][q] != f )
								{
									//i_elems.push_back(faces[0][q]);
									//i_coefs.push_back(-beta*iC12sum*CC);
								}
								C1sum += CC;
							}
						}
						//TODO
						//f_coefs[0] = C1sum*multfbnd;
						//assert(f_elems[0] == f->BackCell());
						//i_coefs[0] = beta*C1sum*iC12sum;
					}
				} //end of loop over faces
				
			}
            std::cout << "Construct W matrix: " << Timer() - ttt << std::endl;
			
			/*
			ttt = Timer();
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				
#if defined(USE_OMP)
#pragma omp for
#endif
				for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
				{
					Face face = m->FaceByLocalID(q);
				}
			}
			std::cout << "Construct fluxes and interpolations: " << Timer() - ttt << std::endl;
			*/
			if( m->HaveTag("UVW") ) //Is there a displacement on the mesh?
				tag_UVW = m->GetTag("UVW"); //Get the pressure
			
			if( !tag_UVW.isValid() || !tag_UVW.isDefined(CELL) ) // Pressure was not initialized or was not defined on nodes
			{
				//srand(1); // Randomization
				tag_UVW = m->CreateTag("UVW",DATA_REAL,CELL|FACE,NONE,3); // Create a new tag for the displacement
				for(Mesh::iteratorElement e = m->BeginElement(CELL|FACE); e != m->EndElement(); ++e) //Loop over mesh cells
				{
					//tag_UVW(*e,3,1).Zero();//(rand()*1.0)/(RAND_MAX*1.0); // Prescribe random value in [0,1]
					tag_UVW[*e][0] = (rand()*1.0)/(RAND_MAX*1.0);
					tag_UVW[*e][1] = (rand()*1.0)/(RAND_MAX*1.0);
					tag_UVW[*e][2] = (rand()*1.0)/(RAND_MAX*1.0);
				}
			}
			
			if( !tag_UVW.isDefined(FACE) )
			{
				tag_UVW = m->CreateTag("UVW",DATA_REAL,FACE,NONE,3);
				for(Mesh::iteratorElement e = m->BeginElement(FACE); e != m->EndElement(); ++e) //Loop over mesh cells
				{
					//tag_UVW(*e,3,1).Zero();//(rand()*1.0)/(RAND_MAX*1.0); // Prescribe random value in [0,1]
					tag_UVW[*e][0] = (rand()*1.0)/(RAND_MAX*1.0);
					tag_UVW[*e][1] = (rand()*1.0)/(RAND_MAX*1.0);
					tag_UVW[*e][2] = (rand()*1.0)/(RAND_MAX*1.0);
				}
			}
			m->ExchangeData(tag_UVW,CELL|FACE,0); //Synchronize initial solution with boundary unknowns
			
			if( reference_solution && m->HaveTag("REFERENCE_SOLUTION") )
			{
				TagRealArray tag_R = m->GetTag("REFERENCE_SOLUTION");
				for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE); it != m->EndElement(); ++it)
					tag_UVW(*it,3,1) = tag_R(*it,3,1);
			}
			
			tag_FLUX = m->CreateTag("FLUX",DATA_REAL,FACE,NONE,3);
			tag_S = m->CreateTag("STRESS",DATA_REAL,CELL,NONE,6);
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
                        Text.SetAnnotation(UVW.Index(cell,0),"Divergence, U " + to_string(UVW.Index(cell,0)));
						Text.SetAnnotation(UVW.Index(cell,1),"Divergence, V " + to_string(UVW.Index(cell,1)));
						Text.SetAnnotation(UVW.Index(cell,2),"Divergence, W " + to_string(UVW.Index(cell,2)));
					}
                }
                for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
                {
                    Face face = m->FaceByLocalID(q);
                    if( face.GetStatus() != Element::Ghost )
                    {
                        if( tag_BC.isValid() && face.HaveData(tag_BC) )
						{
                            Text.SetAnnotation(UVW.Index(face,0),"Boundary condition, U " + to_string(UVW.Index(face,0)));
							Text.SetAnnotation(UVW.Index(face,1),"Boundary condition, V " + to_string(UVW.Index(face,1)));
							Text.SetAnnotation(UVW.Index(face,2),"Boundary condition, W " + to_string(UVW.Index(face,2)));
						}
                        else
						{
                            Text.SetAnnotation(UVW.Index(face,0),"Flux continuity, U " + to_string(UVW.Index(face,0)));
							Text.SetAnnotation(UVW.Index(face,1),"Flux continuity, V " + to_string(UVW.Index(face,1)));
							Text.SetAnnotation(UVW.Index(face,2),"Flux continuity, W " + to_string(UVW.Index(face,2)));
						}
                    }
                }
            }

            std::cout << "Matrix was annotated" << std::endl;
			Storage::real condest = 0;
            do
			{
                R.Clear(); //clean up the residual
                double tttt = Timer();
                //~ int total = 0, dmp = 0;
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
						tag_S(cell,6,1) = tag_Ws(cell,6,3*NF)*T;
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
								
								tag_FLUX(faces[k],3,1) += 0.5*rMatrix(T(3*k,3*(k+1),0,1))/faces[k].Area()*(faces[k].FaceOrientedOutside(cell)?1.:-1.);
							}
							Locks.UnLock(index);
						}
					} //end of loop over cells
					
					
					if( print_matrix )
					{
						std::fstream ff;
						ff.open("mat.txt",std::ios::out);
						for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE); it != m->EndElement(); ++it)
						{
							
							ff << "Element " << ElementTypeName(it->GetElementType()) << ":" << it->LocalID();
							if( it->GetElementType() == FACE && it->getAsFace()->Boundary() )
								ff << " boundary";
							ff << std::endl;
							const int bs = 3;
							std::vector<Storage::real> VR(3);
							BlockRow BR(bs,bs);
							
							for(int q = 0; q < bs; ++q)
							{
								VR[q] = R[UVW.Index(it->self(),q)].GetValue();
								Sparse::Row & r = R[UVW.Index(it->self(),q)].GetRow();
								//r.Print();
								for(INMOST_DATA_ENUM_TYPE l = 0; l < r.Size(); ++l) if( r.GetValue(l) )
									BR(q,r.GetIndex(l)) = r.GetValue(l);
							}
							
							ff << std::setprecision(5);
							ff << std::scientific;
							std::map< unsigned, rMatrix >::iterator jt;
							for(int q = 0; q < bs; ++q)
							{
								
								ff << std::setw(13) << VR[q] << ": ";
								for(jt = BR.Begin(); jt != BR.End(); ++jt)
								{
									for(int l = 0; l < bs; ++l)
										ff << std::setw(13) << jt->second(q,l) << " ";
									ff << "|";
								}
								ff << std::endl;
							}
							
							rMatrix Sum(bs,bs);
							Sum.Zero();
							for(jt = BR.Begin(); jt != BR.End(); ++jt)
								Sum += jt->second;
							ff << "Sum: " << std::endl;
							for(int q = 0; q < bs; ++q)
							{
								for(int l = 0; l < bs; ++l)
									ff << std::setw(13) << Sum(q,l) << " ";
								ff << std::endl;
							}
							ff << "diagonal block: " << UVW.Index(it->self(),0)/bs << std::endl;
							rMatrix Diag = BR[UVW.Index(it->self(),0)/bs*bs];
							rMatrix iDiag = Diag.PseudoInvert(1.0e-13);
							ff << "Diag: " << std::endl;
							for(int q = 0; q < bs; ++q)
							{
								for(int l = 0; l < bs; ++l)
									ff << std::setw(13) << Diag(q,l) << " ";
								ff << std::endl;
							}
							ff << "iDiag: " << std::endl;
							for(int q = 0; q < bs; ++q)
							{
								for(int l = 0; l < bs; ++l)
									ff << std::setw(13) << iDiag(q,l) << " ";
								ff << std::endl;
							}
							
							rMatrix W(bs,bs), U(bs,bs), S(bs,bs), V(bs,bs);
							Sum.Zero();
							bool diag = false;
							for(jt = BR.Begin(); jt != BR.End(); ++jt)
							{
								int nz = 0;
								bool curdiag = false;
								W = jt->second;
								W.SVD(U,S,V);
								SVD2Eigen(U,S,V);
								ff << "Block " << jt->first << std::endl;
								if( jt->first == UVW.Index(it->self(),0)/bs )
								{
									ff << "diagonal" << std::endl;
									diag = true;
									curdiag = true;
								}
								for(int k = 0; k < bs; ++k) if( fabs(S(k,k)) < 1.0e-4 )
									nz++;
								if( curdiag && nz ) ff << "problem: " << bs-nz << "/" << bs << std::endl;
								ff << "A | U | S | V" << std::endl;
								for(int k = 0; k < bs; ++k)
								{
									for(int l = 0; l < bs; ++l)
										ff << std::setw(13) << W(k,l) << " ";
									ff << "|";
									for(int l = 0; l < bs; ++l)
										ff << std::setw(13) << U(k,l) << " ";
									ff << "|";
									for(int l = 0; l < bs; ++l)
										ff << std::setw(13) << S(k,l) << " ";
									ff << "|";
									for(int l = 0; l < bs; ++l)
										ff << std::setw(13) << V(k,l) << " ";
									ff << std::endl;
								}
							}
							if( !diag ) std::cout << "no diagonal block!" << std::endl;
						}
						ff.close();
					}


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
				
				//R.GetJacobian().Save("A.mtx",&Text);

				//Solver S(Solver::INNER_MPTILU2);
                Solver S(Solver::INNER_MPTILUC);
				//Solver S("superlu");
                S.SetParameter("relative_tolerance", "1.0e-14");
                S.SetParameter("absolute_tolerance", "1.0e-12");
                S.SetParameter("drop_tolerance", "1.0e-4");
                S.SetParameter("reuse_tolerance", "1.0e-6");

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
                    condest = 0;//S.Condest(1.0e-12,500);
                    std::cout << "condition number " << condest << "            " << std::endl;
                }
                else
                {
                    std::cout << "Unable to solve: " << S.ReturnReason() << std::endl;
                    break;
                }
				std::cout << "solved in " << Timer() - tttt << "\t\t\t\t" << std::endl;
                ++nit;
            } while( R.Norm() > 1.0e-4 && nit < 10); //check the residual norm
        }
        std::cout << "Solved problem in " << Timer() - ttt << " seconds with " << nit << " iterations " << std::endl;
        
        
        std::fstream stat;
		
		if( m->GetProcessorRank() == 0 )
		{
			stat.open("stat.csv",std::ios::in);
			if( stat.fail() )
			{
				stat.clear();
				stat.open("stat.csv",std::ios::out);
				stat << "cells; faces; ";
				stat << "uvw_C; uvw_L2; ";
				stat << "uvw_f_C; uvw_f_L2; ";
				stat << "F_C; F_L2; ";
				stat << "S_C; S_L2; ";
				stat << "energy; grid; problem;" << std::endl;
				stat.close();
			}
			else stat.close();
			
			stat.open("stat.csv",std::ios::app);
			
			stat << m->TotalNumberOf(CELL) << "; ";
			stat << m->TotalNumberOf(FACE) << "; ";
		}
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
            if( m->GetProcessorRank() == 0 ) 
            {
				std::cout << "Error on cells, C-norm " << C << " L2-norm " << L2 << std::endl;
				stat << C << "; " << L2 << "; ";
			}
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
                if( m->GetProcessorRank() == 0 ) 
                {
					std::cout << "Error on faces, C-norm " << C << " L2-norm " << L2 << std::endl;
					stat << C << "; " << L2 << "; ";
				}
            }
            else 
            {
				std::cout << "Reference solution was not defined on faces" << std::endl;
				stat << "NA; NA; ";
			}
        }
        else stat << "NA; NA; NA; NA; ";
		
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
			if( m->GetProcessorRank() == 0 ) 
			{
				std::cout << "Error of fluxes, C-norm " << C << " L2-norm " << L2 << std::endl;
				stat << C << "; " << L2 << "; ";
			}
		}
		else stat << "NA; NA; ";
		if( m->HaveTag("REFERENCE_STRESS") )
		{
			TagReal tag_E = m->CreateTag("ERROR_STRESS",DATA_REAL,CELL,NONE,1);
			TagRealArray tag_R = m->GetTag("REFERENCE_STRESS");
			real C, L2, volume;
			C = L2 = volume = 0.0;
			for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
            {
                Cell cell = m->CellByLocalID(q);
                real err = (tag_S(cell,6,1) - tag_R(cell,6,1)).FrobeniusNorm();
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
			if( m->GetProcessorRank() == 0 ) 
			{
				std::cout << "Error of stress C-norm " << C << " L2-norm " << L2 << std::endl;
				stat << C << "; " << L2 << "; ";
			}
		}
		else stat << "NA; NA; ";
		
		//compute energy by applying inverse of elastic tensor on stress to get strain
		{
			real energy = 0;
#if defined(USE_OMP)
#pragma omp parallel reduction(+:energy)
#endif
			{
				rMatrix s(6,1), e(6,1);
#if defined(USE_OMP)
#pragma omp for
#endif
				for(int i = 0; i < m->CellLastLocalID(); ++i) if( m->isValidCell(i) )
				{
					Cell c1 = m->CellByLocalID(i);
					//extract stress
					s = rMatrix::FromVector(tag_S[c1].data(),6);
					e = rMatrix::FromTensor(tag_C[c1].data(),21,6).Solve(s);
					//calculate energy
					energy += e.DotProduct(s)*c1.Volume();
					
				}
			}
			energy *= 0.5;
			std::cout << "Energy: " << energy << std::endl;
			stat << energy << "; ";
		}
		
		if( m->HaveTag("GRIDNAME") )
		{
			Tag tagname = m->GetTag("GRIDNAME");
			Storage::bulk_array name = m->self().BulkArray(tagname);
			stat << std::string(name.begin(),name.end()) << ";";
		}
		else stat << "NA;";
		if( m->HaveTag("PROBLEMNAME") )
		{
			Tag tagname = m->GetTag("PROBLEMNAME");
			Storage::bulk_array name = m->self().BulkArray(tagname);
			stat << std::string(name.begin(),name.end()) << ";";
		}
		else stat << "NA;";
		stat << std::endl;
		
		stat.close();
		

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
