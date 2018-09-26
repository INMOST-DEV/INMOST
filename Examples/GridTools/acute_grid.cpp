#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "inmost.h"

std::string mesh_name = "ACUTE_TRIANGULAR";

using namespace INMOST;
const Storage::real pi = 3.1415926535897932384626433832795;
#define V_ID(x, y, z) ((x)*n*nz + (y)*nz + (z))
#define V_ID2_X(x, y, z) ((x)*(n-1)*nz + (y)*nz + (z))
#define V_ID2_Y(x, y, z) ((x)*(n)*nz + (y)*nz + (z))
#define V_ID3(x, y, z, q) (((x)*(n-1)*nz + (y)*nz + (z))*4+(q))

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
  
  mesh = new Mesh();
  
  array<HandleType> primary_nodes(n*n*nz);
  array<HandleType> secondary_nodes_y((n-1)*n*nz); // n-1 by x, n by y
  array<HandleType> secondary_nodes_x((n-1)*n*nz); // n by x, n-1 by y
  array<HandleType> internal_nodes((n-1)*(n-1)*nz*4);

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
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n-1; j++)
    {
      for (int k = 0; k < nz; k++) 
      {
        Storage::real xyz[3];
        xyz[0] = i * 1.0 / (n-1);
        xyz[1] = (j+0.5) * 1.0 / (n-1);
        xyz[2] = k * 1.0 / (nz - 1);
        Node c = mesh->CreateNode(xyz);
        c.Integer(mat) = 1;
        secondary_nodes_x[V_ID2_X(i,j,k)] = c.GetHandle();
      }
    }
  }

  //go through Y lines
  for(int i = 0; i < n-1; i++)
  {
    for(int j = 0; j < n; j++)
    {
      for (int k = 0; k < nz; k++) 
      {
        Storage::real xyz[3];
        xyz[0] = (i+0.5) * 1.0 / (n-1);
        xyz[1] = j * 1.0 / (n-1);
        xyz[2] = k * 1.0 / (nz - 1);
        Node c = mesh->CreateNode(xyz);
        c.Integer(mat) = 2;
        secondary_nodes_y[V_ID2_Y(i,j,k)] = c.GetHandle();
      }
    }
  }

  //go through internals
  for(int i = 0; i < n-1; i++)
  {
    for(int j = 0; j < n-1; j++)
    {
      for (int k = 0; k < nz; k++) 
      {
        Node c;
        Storage::real xyz[3];
        xyz[0] = (i+0.35) * 1.0 / (n-1);
        xyz[1] = (j+0.35) * 1.0 / (n-1);
        xyz[2] = k * 1.0 / (nz - 1);
        c = mesh->CreateNode(xyz);
        c.Integer(mat) = 3;
        internal_nodes[V_ID3(i,j,k,0)] = c.GetHandle();


        xyz[0] = (i+0.3) * 1.0 / (n-1);
        xyz[1] = (j+0.7) * 1.0 / (n-1);
        xyz[2] = k * 1.0 / (nz - 1);
        c = mesh->CreateNode(xyz);
        c.Integer(mat) = 3;
        internal_nodes[V_ID3(i,j,k,1)] = c.GetHandle();


        xyz[0] = (i+0.7) * 1.0 / (n-1);
        xyz[1] = (j+0.3) * 1.0 / (n-1);
        xyz[2] = k * 1.0 / (nz - 1);
        c = mesh->CreateNode(xyz);
        c.Integer(mat) = 3;
        internal_nodes[V_ID3(i,j,k,2)] = c.GetHandle();

        xyz[0] = (i+0.65) * 1.0 / (n-1);
        xyz[1] = (j+0.65) * 1.0 / (n-1);
        xyz[2] = k * 1.0 / (nz - 1);
        c = mesh->CreateNode(xyz);
        c.Integer(mat) = 3;
        internal_nodes[V_ID3(i,j,k,3)] = c.GetHandle();
      }
    }
  }

  

  const INMOST_DATA_INTEGER_TYPE nvf[18] = { 0,2,1, 3,4,5, 0,1,4,3, 1,2,5,4, 0,3,5,2 };
  const INMOST_DATA_INTEGER_TYPE numnodes[5] = { 3, 3, 4, 4, 4 };
  for (int i = 0; i < n-1; i++) 
  {
    for (int j = 0; j < n-1; j++) 
    {
      for (int k = 0; k < nz-1; k++) 
      {
        Cell c;
        ElementArray<Node> verts(mesh,6); 
        verts[0] = Node(mesh,primary_nodes[V_ID(i,j,k)]);
        verts[1] = Node(mesh,secondary_nodes_x[V_ID2_X(i,j,k)]);
        verts[2] = Node(mesh,internal_nodes[V_ID3(i,j,k,0)]);
        verts[3] = Node(mesh,primary_nodes[V_ID(i,j,k+1)]);
        verts[4] = Node(mesh,secondary_nodes_x[V_ID2_X(i,j,k+1)]);
        verts[5] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,0)]);
        
        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 0;



        verts[0] = Node(mesh,primary_nodes[V_ID(i,j,k)]);
        verts[1] = Node(mesh,internal_nodes[V_ID3(i,j,k,0)]);
        verts[2] = Node(mesh,secondary_nodes_y[V_ID2_Y(i,j,k)]);
        verts[3] = Node(mesh,primary_nodes[V_ID(i,j,k+1)]);
        verts[4] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,0)]);
        verts[5] = Node(mesh,secondary_nodes_y[V_ID2_Y(i,j,k+1)]);

        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 1;

        verts[0] = Node(mesh,secondary_nodes_y[V_ID2_Y(i,j,k)]);
        verts[1] = Node(mesh,internal_nodes[V_ID3(i,j,k,0)]);
        verts[2] = Node(mesh,internal_nodes[V_ID3(i,j,k,2)]);
        verts[3] = Node(mesh,secondary_nodes_y[V_ID2_Y(i,j,k+1)]);
        verts[4] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,0)]);
        verts[5] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,2)]);

        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 2;

        verts[0] = Node(mesh,secondary_nodes_y[V_ID2_Y(i,j,k)]);
        verts[1] = Node(mesh,internal_nodes[V_ID3(i,j,k,2)]);
        verts[2] = Node(mesh,primary_nodes[V_ID(i+1,j,k)]);
        verts[3] = Node(mesh,secondary_nodes_y[V_ID2_Y(i,j,k+1)]);
        verts[4] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,2)]);
        verts[5] = Node(mesh,primary_nodes[V_ID(i+1,j,k+1)]);
        
        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 3;

        verts[0] = Node(mesh,primary_nodes[V_ID(i+1,j,k)]);
        verts[1] = Node(mesh,internal_nodes[V_ID3(i,j,k,2)]);
        verts[2] = Node(mesh,secondary_nodes_x[V_ID2_X(i+1,j,k)]);
        verts[3] = Node(mesh,primary_nodes[V_ID(i+1,j,k+1)]);
        verts[4] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,2)]);
        verts[5] = Node(mesh,secondary_nodes_x[V_ID2_X(i+1,j,k+1)]);
        
        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 4;

        verts[0] = Node(mesh,secondary_nodes_x[V_ID2_X(i,j,k)]);
        verts[1] = Node(mesh,internal_nodes[V_ID3(i,j,k,1)]);
        verts[2] = Node(mesh,internal_nodes[V_ID3(i,j,k,0)]);
        verts[3] = Node(mesh,secondary_nodes_x[V_ID2_X(i,j,k+1)]);
        verts[4] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,1)]);
        verts[5] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,0)]);
        
        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 5;

        verts[0] = Node(mesh,internal_nodes[V_ID3(i,j,k,0)]);
        verts[1] = Node(mesh,internal_nodes[V_ID3(i,j,k,1)]);
        verts[2] = Node(mesh,internal_nodes[V_ID3(i,j,k,3)]);
        verts[3] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,0)]);
        verts[4] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,1)]);
        verts[5] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,3)]);
        
        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 6;

        verts[0] = Node(mesh,internal_nodes[V_ID3(i,j,k,0)]);
        verts[1] = Node(mesh,internal_nodes[V_ID3(i,j,k,3)]);
        verts[2] = Node(mesh,internal_nodes[V_ID3(i,j,k,2)]);
        verts[3] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,0)]);
        verts[4] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,3)]);
        verts[5] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,2)]);
        
        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 7;

        verts[0] = Node(mesh,internal_nodes[V_ID3(i,j,k,2)]);
        verts[1] = Node(mesh,internal_nodes[V_ID3(i,j,k,3)]);
        verts[2] = Node(mesh,secondary_nodes_x[V_ID2_X(i+1,j,k)]);
        verts[3] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,2)]);
        verts[4] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,3)]);
        verts[5] = Node(mesh,secondary_nodes_x[V_ID2_X(i+1,j,k+1)]);
        
        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 8;

        verts[0] = Node(mesh,secondary_nodes_x[V_ID2_X(i,j,k)]);
        verts[1] = Node(mesh,primary_nodes[V_ID(i,j+1,k)]);
        verts[2] = Node(mesh,internal_nodes[V_ID3(i,j,k,1)]);
        verts[3] = Node(mesh,secondary_nodes_x[V_ID2_X(i,j,k+1)]);
        verts[4] = Node(mesh,primary_nodes[V_ID(i,j+1,k+1)]);
        verts[5] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,1)]);
        
        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 9;

        verts[0] = Node(mesh,internal_nodes[V_ID3(i,j,k,1)]);
        verts[1] = Node(mesh,primary_nodes[V_ID(i,j+1,k)]);
        verts[2] = Node(mesh,secondary_nodes_y[V_ID2_Y(i,j+1,k)]);
        verts[3] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,1)]);
        verts[4] = Node(mesh,primary_nodes[V_ID(i,j+1,k+1)]);
        verts[5] = Node(mesh,secondary_nodes_y[V_ID2_Y(i,j+1,k+1)]);
        
        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 10;

        verts[0] = Node(mesh,internal_nodes[V_ID3(i,j,k,1)]);
        verts[1] = Node(mesh,secondary_nodes_y[V_ID2_Y(i,j+1,k)]);
        verts[2] = Node(mesh,internal_nodes[V_ID3(i,j,k,3)]);
        verts[3] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,1)]);
        verts[4] = Node(mesh,secondary_nodes_y[V_ID2_Y(i,j+1,k+1)]);
        verts[5] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,3)]);
        
        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 11;

        verts[0] = Node(mesh,internal_nodes[V_ID3(i,j,k,3)]);
        verts[1] = Node(mesh,secondary_nodes_y[V_ID2_Y(i,j+1,k)]);
        verts[2] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k)]);
        verts[3] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,3)]);
        verts[4] = Node(mesh,secondary_nodes_y[V_ID2_Y(i,j+1,k+1)]);
        verts[5] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k+1)]);
        
        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 12;

        verts[0] = Node(mesh,internal_nodes[V_ID3(i,j,k,3)]);
        verts[1] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k)]);
        verts[2] = Node(mesh,secondary_nodes_x[V_ID2_X(i+1,j,k)]);
        verts[3] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,3)]);
        verts[4] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k+1)]);
        verts[5] = Node(mesh,secondary_nodes_x[V_ID2_X(i+1,j,k+1)]);
        
        c = mesh->CreateCell(verts,nvf,numnodes,5).first;
        c.Integer(mat) = 13;

      }
    }
  }


  {
	Storage::bulk_array name = mesh->self()->BulkArray(mesh->CreateTag("GRIDNAME",DATA_BULK,MESH,NONE));
	std::stringstream str;
	str << mesh_name << "_" << n << "_" << nz;
	mesh_name = str.str();
	name.replace(name.begin(),name.end(),mesh_name.begin(),mesh_name.end());
  }

  printf("I'm ready!\n");

  mesh->Save("grid.vtk");
  mesh->Save("grid.pmf");
  mesh->Save("grid.gmv");

  printf("File written!\n");

  delete mesh;
  return 0;
}
