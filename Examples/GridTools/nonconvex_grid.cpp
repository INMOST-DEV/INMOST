#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "inmost.h"

std::string mesh_name = "NONCONVEX";

using namespace INMOST;
const Storage::real pi = 3.1415926535897932384626433832795;
#define V_ID(x, y, z) ((x)*n*nz + (y)*nz + (z))
#define V_ID2(x, y, z) ((x)*(n-1)*nz + (y)*nz + (z))

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
  array<HandleType> secondary_nodes_y((n-1)*(n-1)*nz); // n-1 by x, n-1 by y
  array<HandleType> secondary_nodes_x((n-1)*(n-1)*nz); // n-1 by x, n-1 by y
  

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
    for(int j = 0; j < n-2; j++)
    {
      for (int k = 0; k < nz; k++) 
      {
        Storage::real xyz[3];
        xyz[0] = (i+0.8) * 1.0 / (n-1);
        xyz[1] = (j+1.2) * 1.0 / (n-1);
        xyz[2] = k * 1.0 / (nz - 1);
        Node c = mesh->CreateNode(xyz);
        c.Integer(mat) = 1;
        secondary_nodes_x[V_ID2(i,j,k)] = c.GetHandle();
      }
    }
  }

  //go through Y lines
  for(int i = 0; i < n-2; i++)
  {
    for(int j = 0; j < n-1; j++)
    {
      for (int k = 0; k < nz; k++) 
      {
        Storage::real xyz[3];
        xyz[0] = (i+1.2) * 1.0 / (n-1);
        xyz[1] = (j+0.8) * 1.0 / (n-1);
        xyz[2] = k * 1.0 / (nz - 1);
        Node c = mesh->CreateNode(xyz);
        c.Integer(mat) = 2;
        secondary_nodes_y[V_ID2(i,j,k)] = c.GetHandle();
      }
    }
  }
  

  //lower-left corner of the mesh
  {
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,5,4,3,2,1,  6,7,8,9,10,11,  0,1,7,6,  1,2,8,7,  2,3,9,8,  3,4,10,9,   4,5,11,10,  5,0,6,11 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 6, 6, 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,12); 
    {
      int i = 0, j = 0;
      for(int k = 0; k < nz-1; ++k)
      {
        verts[0] = Node(mesh,primary_nodes[V_ID(i,j,k)]);
        verts[1] = Node(mesh,primary_nodes[V_ID(i+1,j,k)]);
        verts[2] = Node(mesh,secondary_nodes_y[V_ID2(i,j,k)]);
        verts[3] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k)]);
        verts[4] = Node(mesh,secondary_nodes_x[V_ID2(i,j,k)]);
        verts[5] = Node(mesh,primary_nodes[V_ID(i,j+1,k)]);

        verts[6] = Node(mesh,primary_nodes[V_ID(i,j,k+1)]);
        verts[7] = Node(mesh,primary_nodes[V_ID(i+1,j,k+1)]);
        verts[8] = Node(mesh,secondary_nodes_y[V_ID2(i,j,k+1)]);
        verts[9] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k+1)]);
        verts[10] = Node(mesh,secondary_nodes_x[V_ID2(i,j,k+1)]);
        verts[11] = Node(mesh,primary_nodes[V_ID(i,j+1,k+1)]);

        Cell c = mesh->CreateCell(verts,nvf,numnodes,8).first;
        c.Integer(mat) = 0;
      }
    }
  }
  //left boundary of the mesh
  {
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,6,5,4,3,2,1,  7,8,9,10,11,12,13,  0,1,8,7,  1,2,9,8,  2,3,10,9,  3,4,11,10,   4,5,12,11,  5,6,13,12,  6,0,7,13 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 7, 7, 4, 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,14); 
    {
      int i = 0;
      for(int j = 1; j < n-2; ++j)
      {
        for(int k = 0; k < nz-1; ++k)
        {
          verts[0] = Node(mesh,primary_nodes[V_ID(i,j,k)]);
          verts[1] = Node(mesh,secondary_nodes_x[V_ID2(i,j-1,k)]);
          verts[2] = Node(mesh,primary_nodes[V_ID(i+1,j,k)]);
          verts[3] = Node(mesh,secondary_nodes_y[V_ID2(i,j,k)]);
          verts[4] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k)]);
          verts[5] = Node(mesh,secondary_nodes_x[V_ID2(i,j,k)]);
          verts[6] = Node(mesh,primary_nodes[V_ID(i,j+1,k)]);

          verts[7] = Node(mesh,primary_nodes[V_ID(i,j,k+1)]);
          verts[8] = Node(mesh,secondary_nodes_x[V_ID2(i,j-1,k+1)]);
          verts[9] = Node(mesh,primary_nodes[V_ID(i+1,j,k+1)]);
          verts[10] = Node(mesh,secondary_nodes_y[V_ID2(i,j,k+1)]);
          verts[11] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k+1)]);
          verts[12] = Node(mesh,secondary_nodes_x[V_ID2(i,j,k+1)]);
          verts[13] = Node(mesh,primary_nodes[V_ID(i,j+1,k+1)]);

          Cell c = mesh->CreateCell(verts,nvf,numnodes,9).first;
          c.Integer(mat) = 1;
        }
      }
    }
  }
  //left-top corner of the mesh
  {
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,5,4,3,2,1,  6,7,8,9,10,11,  0,1,7,6,  1,2,8,7,  2,3,9,8,  3,4,10,9,   4,5,11,10,  5,0,6,11 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 6, 6, 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,12); 
    {
      int i = 0, j = n-2;
      for(int k = 0; k < nz-1; ++k)
      {
        verts[0] = Node(mesh,primary_nodes[V_ID(i,j,k)]);
        verts[1] = Node(mesh,secondary_nodes_x[V_ID2(i,j-1,k)]);
        verts[2] = Node(mesh,primary_nodes[V_ID(i+1,j,k)]);
        verts[3] = Node(mesh,secondary_nodes_y[V_ID2(i,j,k)]);
        verts[4] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k)]);
        verts[5] = Node(mesh,primary_nodes[V_ID(i,j+1,k)]);

        verts[6] = Node(mesh,primary_nodes[V_ID(i,j,k+1)]);
        verts[7] = Node(mesh,secondary_nodes_x[V_ID2(i,j-1,k+1)]);
        verts[8] = Node(mesh,primary_nodes[V_ID(i+1,j,k+1)]);
        verts[9] = Node(mesh,secondary_nodes_y[V_ID2(i,j,k+1)]);
        verts[10] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k+1)]);
        verts[11] = Node(mesh,primary_nodes[V_ID(i,j+1,k+1)]);

        Cell c = mesh->CreateCell(verts,nvf,numnodes,8).first;
        c.Integer(mat) = 2;
      }
    }
  }
  //bottom boundary of the mesh
  {
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,6,5,4,3,2,1,  7,8,9,10,11,12,13,  0,1,8,7,  1,2,9,8,  2,3,10,9,  3,4,11,10,   4,5,12,11,  5,6,13,12,  6,0,7,13 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 7, 7, 4, 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,14); 
    {
      int j = 0;
      for(int i = 1; i < n-2; ++i)
      {
        for(int k = 0; k < nz-1; ++k)
        {
          verts[0] = Node(mesh,primary_nodes[V_ID(i,j,k)]);
          verts[1] = Node(mesh,primary_nodes[V_ID(i+1,j,k)]);
          verts[2] = Node(mesh,secondary_nodes_y[V_ID2(i,j,k)]);
          verts[3] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k)]);
          verts[4] = Node(mesh,secondary_nodes_x[V_ID2(i,j,k)]);
          verts[5] = Node(mesh,primary_nodes[V_ID(i,j+1,k)]);
          verts[6] = Node(mesh,secondary_nodes_y[V_ID2(i-1,j,k)]);

          verts[7] = Node(mesh,primary_nodes[V_ID(i,j,k+1)]);
          verts[8] = Node(mesh,primary_nodes[V_ID(i+1,j,k+1)]);
          verts[9] = Node(mesh,secondary_nodes_y[V_ID2(i,j,k+1)]);
          verts[10] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k+1)]);
          verts[11] = Node(mesh,secondary_nodes_x[V_ID2(i,j,k+1)]);
          verts[12] = Node(mesh,primary_nodes[V_ID(i,j+1,k+1)]);
          verts[13] = Node(mesh,secondary_nodes_y[V_ID2(i-1,j,k+1)]);
          
          Cell c = mesh->CreateCell(verts,nvf,numnodes,9).first;
          c.Integer(mat) = 3;
        }
      }
    }
  }
  //right-bottom corner of the mesh
  {
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,5,4,3,2,1,  6,7,8,9,10,11,  0,1,7,6,  1,2,8,7,  2,3,9,8,  3,4,10,9,   4,5,11,10,  5,0,6,11 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 6, 6, 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,12); 
    {
      int i = n-2, j = 0;
      for(int k = 0; k < nz-1; ++k)
      {
        verts[0] = Node(mesh,primary_nodes[V_ID(i,j,k)]);
        verts[1] = Node(mesh,primary_nodes[V_ID(i+1,j,k)]);
        verts[2] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k)]);
        verts[3] = Node(mesh,secondary_nodes_x[V_ID2(i,j,k)]);
        verts[4] = Node(mesh,primary_nodes[V_ID(i,j+1,k)]);
        verts[5] = Node(mesh,secondary_nodes_y[V_ID2(i-1,j,k)]);
        
        verts[6] = Node(mesh,primary_nodes[V_ID(i,j,k+1)]);
        verts[7] = Node(mesh,primary_nodes[V_ID(i+1,j,k+1)]);
        verts[8] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k+1)]);
        verts[9] = Node(mesh,secondary_nodes_x[V_ID2(i,j,k+1)]);
        verts[10] = Node(mesh,primary_nodes[V_ID(i,j+1,k+1)]);
        verts[11] = Node(mesh,secondary_nodes_y[V_ID2(i-1,j,k+1)]);

        Cell c = mesh->CreateCell(verts,nvf,numnodes,8).first;
        c.Integer(mat) = 4;
      }
    }
  }
  //internals of the mesh
  {
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,7,6,5,4,3,2,1,  8,9,10,11,12,13,14,15,  0,1,9,8,  1,2,10,9,  2,3,11,10,  3,4,12,11,   4,5,13,12,  5,6,14,13,  6,7,15,14, 7,0,8,15 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 8, 8, 4, 4, 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,16); 
    {
      for(int i = 1; i < n-2; ++i)
      {
        for(int j = 1; j < n-2; ++j)
        {
          for(int k = 0; k < nz-1; ++k)
          {
            verts[0] = Node(mesh,primary_nodes[V_ID(i,j,k)]);
            verts[1] = Node(mesh,secondary_nodes_x[V_ID2(i,j-1,k)]);
            verts[2] = Node(mesh,primary_nodes[V_ID(i+1,j,k)]);
            verts[3] = Node(mesh,secondary_nodes_y[V_ID2(i,j,k)]);
            verts[4] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k)]);
            verts[5] = Node(mesh,secondary_nodes_x[V_ID2(i,j,k)]);
            verts[6] = Node(mesh,primary_nodes[V_ID(i,j+1,k)]);
            verts[7] = Node(mesh,secondary_nodes_y[V_ID2(i-1,j,k)]);

            verts[8] = Node(mesh,primary_nodes[V_ID(i,j,k+1)]);
            verts[9] = Node(mesh,secondary_nodes_x[V_ID2(i,j-1,k+1)]);
            verts[10] = Node(mesh,primary_nodes[V_ID(i+1,j,k+1)]);
            verts[11] = Node(mesh,secondary_nodes_y[V_ID2(i,j,k+1)]);
            verts[12] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k+1)]);
            verts[13] = Node(mesh,secondary_nodes_x[V_ID2(i,j,k+1)]);
            verts[14] = Node(mesh,primary_nodes[V_ID(i,j+1,k+1)]);
            verts[15] = Node(mesh,secondary_nodes_y[V_ID2(i-1,j,k+1)]);
            
            Cell c = mesh->CreateCell(verts,nvf,numnodes,10).first;
            c.Integer(mat) = 5;
          }
        }
      }
    }
  }
  //top-right corner of the mesh
  {
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,5,4,3,2,1,  6,7,8,9,10,11,  0,1,7,6,  1,2,8,7,  2,3,9,8,  3,4,10,9,   4,5,11,10,  5,0,6,11 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 6, 6, 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,12); 
    {
      int i = n-2, j = n-2;
      for(int k = 0; k < nz-1; ++k)
      {
        verts[0] = Node(mesh,primary_nodes[V_ID(i,j,k)]);
        verts[1] = Node(mesh,secondary_nodes_x[V_ID2(i,j-1,k)]);
        verts[2] = Node(mesh,primary_nodes[V_ID(i+1,j,k)]);
        verts[3] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k)]);
        verts[4] = Node(mesh,primary_nodes[V_ID(i,j+1,k)]);
        verts[5] = Node(mesh,secondary_nodes_y[V_ID2(i-1,j,k)]);
        

        verts[6] = Node(mesh,primary_nodes[V_ID(i,j,k+1)]);
        verts[7] = Node(mesh,secondary_nodes_x[V_ID2(i,j-1,k+1)]);
        verts[8] = Node(mesh,primary_nodes[V_ID(i+1,j,k+1)]);
        verts[9] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k+1)]);
        verts[10] = Node(mesh,primary_nodes[V_ID(i,j+1,k+1)]);
        verts[11] = Node(mesh,secondary_nodes_y[V_ID2(i-1,j,k+1)]);

        Cell c = mesh->CreateCell(verts,nvf,numnodes,8).first;
        c.Integer(mat) = 6;
      }
    }
  }
  //top boundary of the mesh
  {
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,6,5,4,3,2,1,  7,8,9,10,11,12,13,  0,1,8,7,  1,2,9,8,  2,3,10,9,  3,4,11,10,   4,5,12,11,  5,6,13,12,  6,0,7,13 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 7, 7, 4, 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,14); 
    {
      int j = n-2;
      for(int i = 1; i < n-2; ++i)
      {
        for(int k = 0; k < nz-1; ++k)
        {
          verts[0] = Node(mesh,primary_nodes[V_ID(i,j,k)]);
          verts[1] = Node(mesh,secondary_nodes_x[V_ID2(i,j-1,k)]);
          verts[2] = Node(mesh,primary_nodes[V_ID(i+1,j,k)]);
          verts[3] = Node(mesh,secondary_nodes_y[V_ID2(i,j,k)]);
          verts[4] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k)]);
          verts[5] = Node(mesh,primary_nodes[V_ID(i,j+1,k)]);
          verts[6] = Node(mesh,secondary_nodes_y[V_ID2(i-1,j,k)]);

          verts[7] = Node(mesh,primary_nodes[V_ID(i,j,k+1)]);
          verts[8] = Node(mesh,secondary_nodes_x[V_ID2(i,j-1,k+1)]);
          verts[9] = Node(mesh,primary_nodes[V_ID(i+1,j,k+1)]);
          verts[10] = Node(mesh,secondary_nodes_y[V_ID2(i,j,k+1)]);
          verts[11] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k+1)]);
          verts[12] = Node(mesh,primary_nodes[V_ID(i,j+1,k+1)]);
          verts[13] = Node(mesh,secondary_nodes_y[V_ID2(i-1,j,k+1)]);

          Cell c = mesh->CreateCell(verts,nvf,numnodes,9).first;
          c.Integer(mat) = 7;
        }
      }
    }
  }
  //right boundary of the mesh
  {
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,6,5,4,3,2,1,  7,8,9,10,11,12,13,  0,1,8,7,  1,2,9,8,  2,3,10,9,  3,4,11,10,   4,5,12,11,  5,6,13,12,  6,0,7,13 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 7, 7, 4, 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,14); 
    {
      int i = n-2;
      for(int j = 1; j < n-2; ++j)
      {
        for(int k = 0; k < nz-1; ++k)
        {
          verts[0] = Node(mesh,primary_nodes[V_ID(i,j,k)]);
          verts[1] = Node(mesh,secondary_nodes_x[V_ID2(i,j-1,k)]);
          verts[2] = Node(mesh,primary_nodes[V_ID(i+1,j,k)]);
          verts[3] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k)]);
          verts[4] = Node(mesh,secondary_nodes_x[V_ID2(i,j,k)]);
          verts[5] = Node(mesh,primary_nodes[V_ID(i,j+1,k)]);
          verts[6] = Node(mesh,secondary_nodes_y[V_ID2(i-1,j,k)]);

          verts[7] = Node(mesh,primary_nodes[V_ID(i,j,k+1)]);
          verts[8] = Node(mesh,secondary_nodes_x[V_ID2(i,j-1,k+1)]);
          verts[9] = Node(mesh,primary_nodes[V_ID(i+1,j,k+1)]);
          verts[10] = Node(mesh,primary_nodes[V_ID(i+1,j+1,k+1)]);
          verts[11] = Node(mesh,secondary_nodes_x[V_ID2(i,j,k+1)]);
          verts[12] = Node(mesh,primary_nodes[V_ID(i,j+1,k+1)]);
          verts[13] = Node(mesh,secondary_nodes_y[V_ID2(i-1,j,k+1)]);

          Cell c = mesh->CreateCell(verts,nvf,numnodes,9).first;
          c.Integer(mat) = 8;
        }
      }
    }
  }
  /*
  int orphaned_nodes = 0;
  for(Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
  {
    if( it->nbAdjElements(CELL) == 0 ) orphaned_nodes++;
  }
  std::cout << "orphan nodes: " << orphaned_nodes << std::endl;
  */

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
