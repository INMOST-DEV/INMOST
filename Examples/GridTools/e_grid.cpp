#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


std::string mesh_name = "E-GRID";

using namespace INMOST;
const Storage::real pi = 3.1415926535897932384626433832795;
#define V_ID(x, y, z) ((x)*n*nz + (y)*nz + (z))
#define V_ID2(x, y, z) ((x)*n*nz + (y)*nz + (z))
#define V_ID3(x, y, z, m) (((x)*(n-1)*nz + (y)*nz + (z))*8 + m)

int main(int argc, char *argv[]) 
{
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " n [nz = 1]" << std::endl;
    return -1;
  }
  Mesh * mesh;
  int nz = 2;
  int n = 16;
  if (argc > 1) 
    n = atoi(argv[1])+1;
  if( argc > 2 )
    nz = atoi(argv[2])+1;
	
	double h = 0.2;
	double w = 0.8;
  
  mesh = new Mesh();
  
  array<HandleType> primary_nodes(n*n*nz);
  array<HandleType> secondary_nodes((n-1)*n*nz); // n-1 by x, n by y
  array<HandleType> internal_nodes((n-1)*(n-1)*nz*8); // n-1 by x, n-1 by y
  

  Tag mat = mesh->CreateTag("MATERIAL",DATA_INTEGER,CELL|NODE,NONE,1);

  for (int i = 0; i < n; i++) 
  {
    for (int j = 0; j < n; j++) 
    {
      for (int k = 0; k < nz; k++) 
      {
        Storage::real xyz[3];
        xyz[0] = i * 1.0 / (n-1);
        xyz[1] = j * 1.0 / (n-1);
        xyz[2] = k * 1.0 / (nz - 1);
        Node c = mesh->CreateNode(xyz);
        primary_nodes[V_ID(i,j,k)] = c.GetHandle();
        c.Integer(mat) = 0;
      }
    }
  }

  //go through X lines
  for(int i = 0; i < n-1; i++)
  {
    for(int j = 0; j < n; j++)
    {
      for (int k = 0; k < nz; k++) 
      {
        Storage::real xyz[3];
        xyz[0] = (i+w) * 1.0 / (n-1);
        xyz[1] = j * 1.0 / (n-1);
        xyz[2] = k * 1.0 / (nz - 1);
        Node c = mesh->CreateNode(xyz);
        c.Integer(mat) = 1;
        secondary_nodes[V_ID2(i,j,k)] = c.GetHandle();
      }
    }
  }

  //go through cells
  for(int i = 0; i < n-1; i++)
  {
    for(int j = 0; j < n-1; j++)
    {
      for (int k = 0; k < nz; k++) 
      {
		  for(int m = 0; m < 8; ++m)
		  {
			Storage::real xyz[3];
			xyz[0] = (i+(m % 2 == 1 ? 1-w:w)) * 1.0 / (n-1);
			xyz[1] = (j+h*(m/2+1)) * 1.0 / (n-1);
			xyz[2] = k * 1.0 / (nz - 1);
			Node c = mesh->CreateNode(xyz);
			c.Integer(mat) = 2;
			internal_nodes[V_ID3(i,j,k,m)] = c.GetHandle();
		  }
      }
    }
  }
  

  //internals of the mesh
  {
	  const INMOST_DATA_INTEGER_TYPE nvf1[] =
	  {
		  11,10,9,8,7,6,5,4,3,2,1,0,
		  0,1,13,12,
		  1,2,14,13,
		  2,3,15,14,
		  3,4,16,15,
		  4,5,17,16,
		  5,6,18,17,
		  6,7,19,18,
		  7,8,20,19,
		  8,9,21,20,
		  9,10,22,21,
		  10,11,23,22,
		  11,0,12,23,
		  12,13,14,15,16,17,18,19,20,21,22,23
	  };
	  const INMOST_DATA_INTEGER_TYPE numnodes1[] =
	  { 12,4,4,4,4,4,4,4,4,4,4,4,4,12 };
	  ElementArray<Node> verts(mesh,24);
	  {
		  for(int i = 0; i < n-1; ++i)
		  {
			  for(int j = 0; j < n-1; ++j)
			  {
				  for(int k = 0; k < nz-1; ++k)
				  {
					  
					  Cell c;
					  
					  verts[0] = Node(mesh,primary_nodes[V_ID(i,j,k)]);
					  verts[1] = Node(mesh,secondary_nodes[V_ID2(i,j,k)]);
					  verts[2] = Node(mesh,internal_nodes[V_ID3(i,j,k,0)]);
					  verts[3] = Node(mesh,internal_nodes[V_ID3(i,j,k,1)]);
					  verts[4] = Node(mesh,internal_nodes[V_ID3(i,j,k,3)]);
					  verts[5] = Node(mesh,internal_nodes[V_ID3(i,j,k,2)]);
					  verts[6] = Node(mesh,internal_nodes[V_ID3(i,j,k,4)]);
					  verts[7] = Node(mesh,internal_nodes[V_ID3(i,j,k,5)]);
					  verts[8] = Node(mesh,internal_nodes[V_ID3(i,j,k,7)]);
					  verts[9] = Node(mesh,internal_nodes[V_ID3(i,j,k,6)]);
					  verts[10] = Node(mesh,secondary_nodes[V_ID2(i,j+1,k)]);
					  verts[11] = Node(mesh,primary_nodes[V_ID(i,j+1,k)]);
					  
					  verts[12] = Node(mesh,primary_nodes[V_ID(i,j,k+1)]);
					  verts[13] = Node(mesh,secondary_nodes[V_ID2(i,j,k+1)]);
					  verts[14] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,0)]);
					  verts[15] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,1)]);
					  verts[16] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,3)]);
					  verts[17] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,2)]);
					  verts[18] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,4)]);
					  verts[19] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,5)]);
					  verts[20] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,7)]);
					  verts[21] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,6)]);
					  verts[22] = Node(mesh,secondary_nodes[V_ID2(i,j+1,k+1)]);
					  verts[23] = Node(mesh,primary_nodes[V_ID(i,j+1,k+1)]);
					  
					  c = mesh->CreateCell(verts,nvf1,numnodes1,14).first;
					  c.Integer(mat) = 0;
					  
					  verts[0] = Node(mesh,secondary_nodes[V_ID2(i,j,k)]);
					  verts[1] = Node(mesh,primary_nodes[V_ID(i+1,j,k)]);
					  verts[2] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k)]);
					  verts[3] = Node(mesh,secondary_nodes[V_ID2(i,j+1,k)]);
					  verts[4] = Node(mesh,internal_nodes[V_ID3(i,j,k,6)]);
					  verts[5] = Node(mesh,internal_nodes[V_ID3(i,j,k,7)]);
					  verts[6] = Node(mesh,internal_nodes[V_ID3(i,j,k,5)]);
					  verts[7] = Node(mesh,internal_nodes[V_ID3(i,j,k,4)]);
					  verts[8] = Node(mesh,internal_nodes[V_ID3(i,j,k,2)]);
					  verts[9] = Node(mesh,internal_nodes[V_ID3(i,j,k,3)]);
					  verts[10] = Node(mesh,internal_nodes[V_ID3(i,j,k,1)]);
					  verts[11] = Node(mesh,internal_nodes[V_ID3(i,j,k,0)]);
					  
					  verts[12] = Node(mesh,secondary_nodes[V_ID2(i,j,k+1)]);
					  verts[13] = Node(mesh,primary_nodes[V_ID(i+1,j,k+1)]);
					  verts[14] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k+1)]);
					  verts[15] = Node(mesh,secondary_nodes[V_ID2(i,j+1,k+1)]);
					  verts[16] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,6)]);
					  verts[17] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,7)]);
					  verts[18] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,5)]);
					  verts[19] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,4)]);
					  verts[20] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,2)]);
					  verts[21] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,3)]);
					  verts[22] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,1)]);
					  verts[23] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,0)]);
					  
					  c = mesh->CreateCell(verts,nvf1,numnodes1,14).first;
					  c.Integer(mat) = 1;
					  
				  }
			  }
		  }
	  }
  }
	

  Storage::bulk_array name = mesh->self()->BulkArray(mesh->CreateTag("GRIDNAME",DATA_BULK,MESH,NONE));
  std::stringstream str;
  str << mesh_name << "_" << n << "x" << n << "x" << nz;
  mesh_name = str.str();
  name.replace(name.begin(),name.end(),mesh_name.begin(),mesh_name.end());

  printf("I'm ready!\n");

  mesh->Save("grid.vtk");
  mesh->Save("grid.pmf");
  mesh->Save("grid.gmv");

  printf("File written!\n");


  delete mesh;
  return 0;
}
