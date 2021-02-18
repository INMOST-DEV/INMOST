#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

using namespace INMOST;

#define MESH_SIZE 32

#define V_ID(x, y, z) ((x-localstart[0])*(localsize[1]+1)*(localsize[2]+1) + (y-localstart[1])*(localsize[2]+1) + (z-localstart[2]))

/*      (4)*-------*(6)
          /|      /|
         /       / |
        /  |    /  |
    (5)*-------*(7)|
       |   |   |   |
       |       |   |
       |   |   |   |
       |(0)*- -|- -*(2)
       |  /    |  /
       |       | /
       |/      |/
    (1)*-------*(3)      */
void CreateCubeElement(Mesh *m, ElementArray<Node> verts)
{
	// Define six cube faces assuming verts are numerated in the way presented above
	const INMOST_DATA_INTEGER_TYPE face_nodes[24] = {0,4,6,2, 1,3,7,5, 0,1,5,4, 2,6,7,3, 0,2,3,1, 4,5,7,6};
	const INMOST_DATA_INTEGER_TYPE num_nodes[6]   = {4,       4,       4,       4,       4,       4};

	m->CreateCell(verts,face_nodes,num_nodes,6); // Create the cubic cell in the mesh
}

/*      (4)*-------*(6)
          /|  __/ /|
         / __/   / |
        / /|    /  |
    (5)*-------*(7)|
       |   |   |   |
       |       |   |
       |   |   |   |
       |(0)*- -|- -*(2)
       |  /    |/ /
       |    _/ | /
       |/_/    |/
    (1)*-------*(3)      */
void CreateNEPrismElements(Mesh *m, ElementArray<Node> verts)
{
	// Define prism faces assuming verts are numerated in the way presented above
	const INMOST_DATA_INTEGER_TYPE nw_face_nodes[18] = {0,4,6,2, 2,6,5,1, 1,5,4,0, 0,2,1, 4,5,6};
	const INMOST_DATA_INTEGER_TYPE nw_num_nodes[5]   = {4,       4,       4,       3,     3};
	const INMOST_DATA_INTEGER_TYPE se_face_nodes[18] = {1,5,6,2, 1,3,7,5, 7,3,2,6, 1,2,3, 6,5,7};
	const INMOST_DATA_INTEGER_TYPE se_num_nodes[5]   = {4,       4,       4,       3,     3};

	m->CreateCell(verts,nw_face_nodes,nw_num_nodes,5); // Create north-west prismatic cell
	m->CreateCell(verts,se_face_nodes,se_num_nodes,5); // Create south-east prismatic cell
}

/*      (4)*-------*(6)
          /|\     /|
         /   \   / |
        /  |  \ /  |
    (5)*-------*(7)|
       |   |   |   |
       |       |   |
       |   |   |   |
       |(0)*- -|- -*(2)
       |  / \  |  /
       |       | /
       |/     \|/
    (1)*-------*(3)      */
void CreateNWPrismElements(Mesh *m, ElementArray<Node> verts)
{
	// Define prism faces assuming verts are numerated in the way presented above
	const INMOST_DATA_INTEGER_TYPE ne_face_nodes[18] = {0,4,6,2, 0,3,7,4, 2,6,7,3, 0,2,3, 4,7,6};
	const INMOST_DATA_INTEGER_TYPE ne_num_nodes[5]   = {4,       4,       4,       3,     3};
	const INMOST_DATA_INTEGER_TYPE sw_face_nodes[18] = {0,4,7,3, 1,3,7,5, 0,1,5,4, 0,3,1, 4,5,7};
	const INMOST_DATA_INTEGER_TYPE sw_num_nodes[5]   = {4,       4,       4,       3,     3};

	m->CreateCell(verts,ne_face_nodes,ne_num_nodes,5); // Create north-east prismatic cell
	m->CreateCell(verts,sw_face_nodes,sw_num_nodes,5); // Create south-west prismatic cell
}


/*      (4)*-------*(6)
          /|\     /|
         /   \   / |
        /  |  \ /  |
    (5)*-------*(7)|
       |   |   |   |
       |       |   |
       |   |   |   |
       |(0)*- -|- -*(2)
       |  / \  |  /
       |       | /
       |/     \|/
    (1)*-------*(3)      */
void CreateNWTetElements(Mesh *m, ElementArray<Node> verts)
{
	// Define prism faces assuming verts are numerated in the way presented above
	const INMOST_DATA_INTEGER_TYPE ne_face_nodes1[12] = {0,1,5,  5,1,3,   1,0,3, 3,0,5};
	const INMOST_DATA_INTEGER_TYPE ne_num_nodes1[4]   = {3,3,3,3};

	const INMOST_DATA_INTEGER_TYPE ne_face_nodes2[12] = {0,3,5, 0,7,3, 5,3,7, 0,5,7};
	const INMOST_DATA_INTEGER_TYPE ne_num_nodes2[4]   = {3,3,3,3};

	const INMOST_DATA_INTEGER_TYPE ne_face_nodes3[12] = {0,7,5, 4,5,7, 0,5,4, 0,4,7};
	const INMOST_DATA_INTEGER_TYPE ne_num_nodes3[4]   = {3,3,3,3};


	const INMOST_DATA_INTEGER_TYPE sw_face_nodes1[12] = {0,3,7, 2,7,3, 0,7,2, 0,2,3};
	const INMOST_DATA_INTEGER_TYPE sw_num_nodes1[4]   = {3,3,3,3};

	const INMOST_DATA_INTEGER_TYPE sw_face_nodes2[12] = {0,7,4, 0,2,7, 2,4,7, 0,4,2};
	const INMOST_DATA_INTEGER_TYPE sw_num_nodes2[4]   = {3,3,3,3};

	const INMOST_DATA_INTEGER_TYPE sw_face_nodes3[12] = {4,6,2, 6,7,2, 4,7,6, 4,2,7};
	const INMOST_DATA_INTEGER_TYPE sw_num_nodes3[4]   = {3,3,3,3};

	m->CreateCell(verts,ne_face_nodes1,ne_num_nodes1,4); // Create north-east prismatic cell
	m->CreateCell(verts,ne_face_nodes2,ne_num_nodes2,4); // Create north-east prismatic cell
	m->CreateCell(verts,ne_face_nodes3,ne_num_nodes3,4); // Create north-east prismatic cell
	m->CreateCell(verts,sw_face_nodes1,sw_num_nodes1,4); // Create south-west prismatic cell
	m->CreateCell(verts,sw_face_nodes2,sw_num_nodes2,4); // Create south-west prismatic cell
	m->CreateCell(verts,sw_face_nodes3,sw_num_nodes3,4); // Create south-west prismatic cell
}

Mesh * ParallelGenerator(INMOST_MPI_Comm comm, int ng, int nx, int ny, int nz)
{
	int procs_per_axis[3] = {1,1,1};
	int sizes[3] = {nx,ny,nz};
	int rank,size;
	Mesh * m = new Mesh(); // Create a mesh to be constructed

	m->SetCommunicator(comm); // Set the MPI communicator, usually MPI_COMM_WORLD

#if defined(USE_MPI)
	MPI_Comm_set_errhandler(comm,MPI_ERRORS_RETURN);
	//MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
#endif

	rank = m->GetProcessorRank(); // Get the rank of the current process
	size = m->GetProcessorsNumber(); // Get the number of processors used in communicator comm

	// Compute the configuration of processes connection
	{
		int divsize = size;
		std::vector<int> divs;
		while( divsize > 1 )
		{
			for(int k = 2; k <= divsize; k++)
				if( divsize % k == 0 )
				{
					divs.push_back(k);
					divsize /= k;
					break;
				}
		}
		int elements_per_procs[3] = {sizes[0],sizes[1],sizes[2]};
		for(std::vector<int>::reverse_iterator it = divs.rbegin(); it != divs.rend(); it++)
		{
			int * max = std::max_element(elements_per_procs+0,elements_per_procs+3);
			procs_per_axis[max-elements_per_procs] *= *it;
			(*max) /= *it;
		}
	}

	//rank = proc_coords[2] * procs_per_axis[0] *procs_per_axis[1] + proc_coords[1] * procs_per_axis[0] + proc_coords[0];
	int proc_coords[3] = {rank % procs_per_axis[0] , rank / procs_per_axis[0] % procs_per_axis[1], rank / (procs_per_axis[0] *procs_per_axis[1]) };

	int localsize[3], localstart[3], localend[3];
	int avgsize[3] =
			{
					(int)ceil((double)sizes[0]/procs_per_axis[0]),
					(int)ceil((double)sizes[1]/procs_per_axis[1]),
					(int)ceil((double)sizes[2]/procs_per_axis[2])
			};

	for(int j = 0; j < 3; j++)
	{
		localstart[j] = avgsize[j] * proc_coords[j];
		if( proc_coords[j] == procs_per_axis[j] - 1 )
			localsize[j] = sizes[j] - avgsize[j] * (procs_per_axis[j]-1);
		else localsize[j] = avgsize[j];
		localend[j] = localstart[j] + localsize[j];
	}

	// Create i-j-k structure of nodes
	ElementArray<Node> newverts(m);
	newverts.reserve(localsize[0]*localsize[1]*localsize[2]);

	for(int i = localstart[0]; i <= localend[0]; i++)
		for(int j = localstart[1]; j <= localend[1]; j++)
			for(int k = localstart[2]; k <= localend[2]; k++)
			{
				Storage::real xyz[3];
				xyz[0] = i * 1.0 / (sizes[0]);
				xyz[1] = j * 1.0 / (sizes[1]);
				xyz[2] = k * 1.0 / (sizes[2]);
				newverts.push_back(m->CreateNode(xyz)); // Create node in the mesh with index V_ID(i,j,k)
			}

	// Create i-j-k structure of elements
	for(int i = localstart[0]+1; i <= localend[0]; i++)
		for(int j = localstart[1]+1; j <= localend[1]; j++)
			for(int k = localstart[2]+1; k <= localend[2]; k++)
			{
				// Create local array of eight nodes                           /*      (4)*-------*(6)  */
				// using representation on the right figure                    /*        /|      /|     */
				ElementArray<Node> verts(m);                                   /*       /       / |     */
				verts.push_back(newverts[V_ID(i - 1,j - 1, k - 1)]); // 0      /*      /  |    /  |     */
				verts.push_back(newverts[V_ID(i - 0,j - 1, k - 1)]); // 1      /*  (5)*-------*(7)|     */
				verts.push_back(newverts[V_ID(i - 1,j - 0, k - 1)]); // 2      /*     |   |   |   |     */
				verts.push_back(newverts[V_ID(i - 0,j - 0, k - 1)]); // 3      /*     |       |   |     */
				verts.push_back(newverts[V_ID(i - 1,j - 1, k - 0)]); // 4      /*     |   |   |   |     */
				verts.push_back(newverts[V_ID(i - 0,j - 1, k - 0)]); // 5      /*     |(0)*- -|- -*(2)  */
				verts.push_back(newverts[V_ID(i - 1,j - 0, k - 0)]); // 6      /*     |  /    |  /      */
				verts.push_back(newverts[V_ID(i - 0,j - 0, k - 0)]); // 7      /*     |       | /       */
				/*     |/      |/        */
				// Create cells based on parameter ng                          /*  (1)*-------*(3)      */
				if (ng == 5) // Create tetrahedral grid
				{
					CreateNWTetElements(m,verts);
				}
				else if (ng == 4) // Create cubic cell
				{
					CreateCubeElement(m,verts);
				}
				else if ((i + j) % 2 == 0 || ng == 6) // Create two prism cells
				{
					CreateNWPrismElements(m,verts);
				}
				else // Two prism cells with different diagonal direction
				{
					CreateNEPrismElements(m,verts);
				}
			}

	m->ResolveShared(); // Resolve duplicate nodes

	m->ExchangeGhost(2,NODE); // Construct Ghost cells in 2 layers connected via nodes
	return m;
}

int main(int argc, char *argv[])
{
	std::string mesh_name;
	Mesh * mesh;
	Mesh::Initialize(&argc,&argv);

	// Specify mesh configuration
	int ng = 3, nx = MESH_SIZE, ny = MESH_SIZE, nz = MESH_SIZE, args = 0;
	if( argc > 4 )
	{
		ng = atoi(argv[1]); // 3 - Prismatic, 4 - Cubic generator
		nx = atoi(argv[2]);
		ny = atoi(argv[3]);
		nz = atoi(argv[4]);
		args = 1;
	}

	// Construct a mesh in parallel
	//MPI_Barrier(mesh->GetCommunicator());
	Storage::real tt = Timer();
	mesh = ParallelGenerator(INMOST_MPI_COMM_WORLD,ng,nx,ny,nz);
	tt = Timer() - tt;
	tt = mesh->Integrate(tt)/mesh->GetProcessorsNumber(); //???
	if( mesh->GetProcessorRank() == 0 )
	{
		if( args == 0 ) std::cout << "Usage: " << argv[0] << " ng nx ny nz [output.[p]vtk]" << std::endl << "ng - 3 for prismatic mesh, 4 for cubic mesh" << std::endl;
		if( args == 0 ) std::cout << "Default ";
		if( ng == 5 )
		{
			std::cout << "Tetrahedral ";
			std::stringstream str;
			str << "TETRAHEDRAL_" << nx << "x" << ny << "x" << nz;
			mesh_name = str.str();
		}
		else if( ng == 4 )
		{
			std::cout << "Cubic ";
			std::stringstream str;
			str << "CUBIC_" << nx << "x" << ny << "x" << nz;
			mesh_name = str.str();
		}
		else
		{
			std::cout << "Prismatic ";
			std::stringstream str;
			str << "PRISMATIC_" << nx << "x" << ny << "x" << nz;
			mesh_name = str.str();
		}
		std::cout << "Grid: " << nx << " x " << ny << " x " << nz << std::endl;
		std::cout << "Processors: " <<mesh->GetProcessorsNumber() << std::endl;
		std::cout << "Mesh generator time: " << tt << std::endl;
	}

	std::string filename;
	if (argc > 5)
	{
		filename = argv[5];
	}
	else
	{
		filename = "grid";
		if( mesh->GetProcessorsNumber() == 1 )
			filename += ".vtk";
		else
			filename += ".pvtk";
	}

	Storage::bulk_array name = mesh->self()->BulkArray(mesh->CreateTag("GRIDNAME",DATA_BULK,MESH,NONE));
	name.replace(name.begin(),name.end(),mesh_name.begin(),mesh_name.end());

#if defined(USE_MPI)
	MPI_Barrier(mesh->GetCommunicator());
#endif
	tt = Timer();
	mesh->Save(filename); // Save constructed mesh to the file
#if defined(USE_MPI)
	MPI_Barrier(mesh->GetCommunicator());
#endif
	tt = Timer() - tt;
	if( mesh->GetProcessorRank() == 0 ) std::cout << "Save to file \"" << filename << "\" time: " << tt << std::endl;

	delete mesh;
	Mesh::Finalize();
	return 0;
}
