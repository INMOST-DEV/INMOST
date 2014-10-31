#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "../../inmost.h"
#include <string.h>

using namespace INMOST;
#define MESH_SIZE 64

#define V_ID(x, y, z) ((x-localstart[0])*(localsize[1]+1)*(localsize[2]+1) + (y-localstart[1])*(localsize[2]+1) + (z-localstart[2]))

Mesh * ParallelCubeGenerator(INMOST_MPI_Comm comm, int nx, int ny, int nz)
{
	int procs_per_axis[3] = {1,1,1};
	int sizes[3] = {nx,ny,nz};
	int rank,size;
	Mesh * m = new Mesh();

	m->SetCommunicator(comm);
	
#if defined(USE_MPI)
	MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
#endif
	
	rank = m->GetProcessorRank();
	size = m->GetProcessorsNumber();
	
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
	
	//???
	int localsize[3], localstart[3], localend[3];
	int avgsize[3] = 
	{
		(int)ceil(sizes[0]/procs_per_axis[0]),
		(int)ceil(sizes[1]/procs_per_axis[1]),
		(int)ceil(sizes[2]/procs_per_axis[2])
	};
	
	for(int j = 0; j < 3; j++)
	{
		localstart[j] = avgsize[j] * proc_coords[j];
		if( proc_coords[j] == procs_per_axis[j] - 1 )
			localsize[j] = sizes[j] - avgsize[j] * (procs_per_axis[j]-1);
		else localsize[j] = avgsize[j];
		localend[j] = localstart[j] + localsize[j];
	}
	
	std::vector<Node *> newverts;
	newverts.reserve(localsize[0]*localsize[1]*localsize[2]);
	
	for(int i = localstart[0]; i <= localend[0]; i++)
		for(int j = localstart[1]; j <= localend[1]; j++)
			for(int k = localstart[2]; k <= localend[2]; k++)
			{
				Storage::real xyz[3];
				xyz[0] = i * 1.0 / (sizes[0]);
				xyz[1] = j * 1.0 / (sizes[1]);
				xyz[2] = k * 1.0 / (sizes[2]);
				newverts.push_back(m->CreateNode(xyz));
				if( ((int)newverts.size()-1) != V_ID(i,j,k))
					printf("v_id = %ld, [%d,%d,%d] = %d\n", newverts.size()-1,i,j,k, V_ID(i, j, k));
			}
			
	std::vector< std::vector<Node *> > c_f_nodes(6);
	for (int i = 0; i < 6; i++)
		c_f_nodes[i].resize(4);
	
	for(int i = localstart[0]+1; i <= localend[0]; i++)
		for(int j = localstart[1]+1; j <= localend[1]; j++)
			for(int k = localstart[2]+1; k <= localend[2]; k++)
			{
				const INMOST_DATA_ENUM_TYPE nvf[24] = {0,4,6,2,1,3,7,5,0,1,5,4,2,6,7,3,0,2,3,1,4,5,7,6};
				const INMOST_DATA_ENUM_TYPE numnodes[6] = { 4, 4, 4, 4, 4, 4 };
				Node * verts[8]; 
				verts[0] = newverts[V_ID(i - 1,j - 1, k - 1)];
				verts[1] = newverts[V_ID(i - 0,j - 1, k - 1)];
				verts[2] = newverts[V_ID(i - 1,j - 0, k - 1)];
				verts[3] = newverts[V_ID(i - 0,j - 0, k - 1)];
				verts[4] = newverts[V_ID(i - 1,j - 1, k - 0)];
				verts[5] = newverts[V_ID(i - 0,j - 1, k - 0)];
				verts[6] = newverts[V_ID(i - 1,j - 0, k - 0)];
				verts[7] = newverts[V_ID(i - 0,j - 0, k - 0)];
				
				m->CreateCell(verts,nvf,numnodes,6).first;
			}
			
	m->ResolveShared();
	
	m->ExchangeGhost(2,NODE);
	//m->Save("test.pvtk");
	return m;
}

Mesh * ParallelCubePrismGenerator(INMOST_MPI_Comm comm, int nx, int ny, int nz)
{
	int procs_per_axis[3] = {1,1,1};
	int sizes[3] = {nx,ny,nz};
	int rank,size;
	Mesh * m = new Mesh();

	m->SetCommunicator(comm);
	
#if defined(USE_MPI)
	MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
#endif
	
	rank = m->GetProcessorRank();
	size = m->GetProcessorsNumber();
	
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
	
	//???
	int localsize[3], localstart[3], localend[3];
	int avgsize[3] = 
	{
		(int)ceil(sizes[0]/procs_per_axis[0]),
		(int)ceil(sizes[1]/procs_per_axis[1]),
		(int)ceil(sizes[2]/procs_per_axis[2])
	};
	
	for(int j = 0; j < 3; j++)
	{
		localstart[j] = avgsize[j] * proc_coords[j];
		if( proc_coords[j] == procs_per_axis[j] - 1 )
			localsize[j] = sizes[j] - avgsize[j] * (procs_per_axis[j]-1);
		else localsize[j] = avgsize[j];
		localend[j] = localstart[j] + localsize[j];
	}
	
	std::vector<Node *> newverts;
	newverts.reserve(localsize[0]*localsize[1]*localsize[2]);
	
	for(int i = localstart[0]; i <= localend[0]; i++)
		for(int j = localstart[1]; j <= localend[1]; j++)
			for(int k = localstart[2]; k <= localend[2]; k++)
			{
				Storage::real xyz[3];
				xyz[0] = i * 1.0 / (sizes[0]);
				xyz[1] = j * 1.0 / (sizes[1]);
				xyz[2] = k * 1.0 / (sizes[2]);
				newverts.push_back(m->CreateNode(xyz));
				if( ((int)newverts.size()-1) != V_ID(i,j,k))
					printf("v_id = %ld, [%d,%d,%d] = %d\n", newverts.size()-1,i,j,k, V_ID(i, j, k));
			}
			
	std::vector< std::vector<Node *> > c_f_nodes(5);
	for (int i = 0; i < 3; i++)
		c_f_nodes[i].resize(4);
	for (int i = 3; i < 5; i++)
		c_f_nodes[i].resize(3);
	
	for(int i = localstart[0]+1; i <= localend[0]; i++)
		for(int j = localstart[1]+1; j <= localend[1]; j++)
			for(int k = localstart[2]+1; k <= localend[2]; k++)
			{
				const INMOST_DATA_ENUM_TYPE NE_nvf1[18] = {0,4,6,2,0,3,7,4,2,6,7,3,0,2,3,4,7,6};
				const INMOST_DATA_ENUM_TYPE NE_nvf2[18] = {0,4,7,3,1,3,7,5,0,1,5,4,0,3,1,4,5,7};
				
				const INMOST_DATA_ENUM_TYPE NE_nvf3[18] = {0,4,6,2,2,6,5,1,1,5,4,0,0,2,1,4,5,6};
				const INMOST_DATA_ENUM_TYPE NE_nvf4[18] = {1,5,6,2,1,3,7,5,7,3,2,6,1,2,3,6,5,7};
				
				const INMOST_DATA_ENUM_TYPE numnodes[5] = {4,4,4,3,3};
				
				Node * verts[8]; 
				verts[0] = newverts[V_ID(i - 1,j - 1, k - 1)];
				verts[1] = newverts[V_ID(i - 0,j - 1, k - 1)];
				verts[2] = newverts[V_ID(i - 1,j - 0, k - 1)];
				verts[3] = newverts[V_ID(i - 0,j - 0, k - 1)];
				verts[4] = newverts[V_ID(i - 1,j - 1, k - 0)];
				verts[5] = newverts[V_ID(i - 0,j - 1, k - 0)];
				verts[6] = newverts[V_ID(i - 1,j - 0, k - 0)];
				verts[7] = newverts[V_ID(i - 0,j - 0, k - 0)];
				
				if ((i + j) % 2 == 0)
				{
					m->CreateCell(verts,NE_nvf1,numnodes,5).first;
					m->CreateCell(verts,NE_nvf2,numnodes,5).first;
				}
				
				else
				{
					m->CreateCell(verts,NE_nvf3,numnodes,5).first;
					m->CreateCell(verts,NE_nvf4,numnodes,5).first;
				}
			}
			
	m->ResolveShared();
	
	m->ExchangeGhost(2,NODE);
	//m->Save("test.pvtk");
	return m;
}


/****************************************************************************/
int main(int argc, char *argv[]) 
{
	Mesh * mesh;
	Mesh::Initialize(&argc,&argv);

	double tt = Timer();
	int nx = MESH_SIZE, ny = MESH_SIZE, nz = 1;
	
	if( argc > 3 )
	{
		nx = atoi(argv[1]);
		ny = atoi(argv[2]);
		nz = atoi(argv[3]);
	}
	else
		if( mesh->GetProcessorRank() == 0 ) std::cout << "Default grid: " << nx << " x " << ny << " x " << nz << std::endl;
	
// 	mesh = ParallelCubeGenerator(INMOST_MPI_COMM_WORLD,nx,ny,nz);
	mesh = ParallelCubePrismGenerator(INMOST_MPI_COMM_WORLD,nx,ny,nz);
	tt = Timer() - tt;
	tt = mesh->Integrate(tt)/mesh->GetProcessorsNumber();
	if( mesh->GetProcessorRank() == 0 ) std::cout << "Gen time: " << tt << std::endl;
	if( mesh->GetProcessorRank() == 0 ) std::cout << "Procs: " <<mesh->GetProcessorsNumber() << std::endl;	
	
	if (mesh->GetProcessorsNumber() == 1 )
		mesh->Save("grid.vtk");
	else
		mesh->Save("grid.pvtk");
	
	MPI_Barrier(mesh->GetCommunicator());
	if( mesh->GetProcessorRank() == 0 )
	{
//		printf ("Initialization time: %.2f sec. Time step = %.2f\nComputation time: %.2f sec, total iterations = %d (linear = %d)\n",
//				tinit, timestep, tcomp, totnit, totlit);
	}
	
	delete mesh;
	Mesh::Finalize();
	return 0;
}
