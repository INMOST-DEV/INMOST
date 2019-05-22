#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

std::string mesh_name = "HEXAGONAL";

using namespace INMOST;
const Storage::real pi = 3.1415926535897932384626433832795;
#define V_ID_B(y,z) ((y)*nz + (z))
#define V_ID_P(x, y, z) ((x)*(n+2)*nz + (y)*nz + (z))
#define V_ID_S(x, y, z) ((x)*(n+1)*nz + (y)*nz + (z))



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
    n = atoi(argv[1]);
  if( argc > 2 )
    nz = atoi(argv[2]);
  
  mesh = new Mesh();
  
  array<HandleType> lb((n+1)*nz), rb((n+1)*nz);
  array<HandleType> primary_nodes((n)*(n+2)*nz); // n by x, n+1 by y
  array<HandleType> secondary_nodes((n)*(n+1)*nz);  // n by x, n by y
  
  Storage::real h = 1.0/n;
  //Storage::real hb = 1.65*h; //(hb+hs = 2.5h)
  //Storage::real hs = 0.85*h;
  Storage::real hb = 1.2*h; //(hb+hs = 2h)
  Storage::real hs = 0.8*h;


  Tag mat = mesh->CreateTag("MATERIAL",DATA_INTEGER,CELL|NODE,NONE,1);

  Storage::real xyz[3];

  //lb
  xyz[0] = 0.0;
  xyz[1] = 0.0;
  for (int j = 0; j < n+1; j++) 
  {
    for (int k = 0; k < nz; k++) 
    {
      xyz[2] = k * 1.0 / (nz-1);
      Node c = mesh->CreateNode(xyz);
      lb[V_ID_B(j,k)] = c.GetHandle();
      c.Integer(mat) = 0;
    }
    xyz[1] += h;
  }
  //rb
  if( n%2 == 0 )
  {
    xyz[0] = 1.0;
    xyz[1] = 0.0;
    for (int j = 0; j < n+1; j++) 
    {
      for (int k = 0; k < nz; k++) 
      {
        xyz[2] = k * 1.0 / (nz-1);
        Node c = mesh->CreateNode(xyz);
        rb[V_ID_B(j,k)] = c.GetHandle();
        c.Integer(mat) = 1;
      }
      xyz[1] += h;
    }
  }
  else
  {
    xyz[0] = 1.0;
    xyz[1] = 0.5*h;
    for (int j = 0; j < n; j++) 
    {
      for (int k = 0; k < nz; k++) 
      {
        xyz[2] = k * 1.0 / (nz-1);
        Node c = mesh->CreateNode(xyz);
        rb[V_ID_B(j,k)] = c.GetHandle();
        c.Integer(mat) = 1;
      }
      xyz[1] += h;
    }
  }

  //primary nodes

  xyz[0] = 0.5*hs;
  for (int i = 0; i < n; i++) 
  {
    xyz[1] = 0.0;
    for (int j = 0; j < n+1; j++) 
    {
      for (int k = 0; k < nz; k++) 
      {
        xyz[2] = k * 1.0 / (nz - 1);
        Node c = mesh->CreateNode(xyz);
        primary_nodes[V_ID_P(i,j,k)] = c.GetHandle();
        c.Integer(mat) = 2;
      }
      xyz[1] += h;
    }
    xyz[0] += (i%2 ? hs : hb);
  }

  //secondary nodes
  xyz[0] = 0.5*hb;
  for (int i = 0; i < n; i++) 
  {
    xyz[1] = 0.5*h;
    for (int j = 0; j < n; j++) 
    {
      for (int k = 0; k < nz; k++) 
      {
        xyz[2] = k * 1.0 / (nz - 1);
        Node c = mesh->CreateNode(xyz);
        secondary_nodes[V_ID_S(i,j,k)] = c.GetHandle();
        c.Integer(mat) = 3;
      }
      xyz[1] += h;
    }
    xyz[0] += (i%2 ? hb : hs);
  }

  {
    //even internal hexahedra's
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,5,4,3,2,1,  6,7,8,9,10,11,  0,1,7,6,  1,2,8,7,  2,3,9,8,  3,4,10,9,   4,5,11,10,  5,0,6,11 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 6, 6, 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,12);
    for(int i = 0; i < n/2; ++i)
    {
      for(int j = 0; j < n-1; ++j)
      {
        for(int k = 0; k < nz-1; ++k)
        {
          verts[0] = Node(mesh,primary_nodes[V_ID_P(i*2,j+1,k)]);
          verts[1] = Node(mesh,secondary_nodes[V_ID_S(i*2,j+1,k)]);
          verts[2] = Node(mesh,secondary_nodes[V_ID_S(i*2+1,j+1,k)]);
          verts[3] = Node(mesh,primary_nodes[V_ID_P(i*2+1,j+1,k)]);
          verts[4] = Node(mesh,secondary_nodes[V_ID_S(i*2+1,j,k)]);
          verts[5] = Node(mesh,secondary_nodes[V_ID_S(i*2,j,k)]);

          verts[6] = Node(mesh,primary_nodes[V_ID_P(i*2,j+1,k+1)]);
          verts[7] = Node(mesh,secondary_nodes[V_ID_S(i*2,j+1,k+1)]);
          verts[8] = Node(mesh,secondary_nodes[V_ID_S(i*2+1,j+1,k+1)]);
          verts[9] = Node(mesh,primary_nodes[V_ID_P(i*2+1,j+1,k+1)]);
          verts[10] = Node(mesh,secondary_nodes[V_ID_S(i*2+1,j,k+1)]);
          verts[11] = Node(mesh,secondary_nodes[V_ID_S(i*2,j,k+1)]);

          Cell c = mesh->CreateCell(verts,nvf,numnodes,8).first;
          c.Integer(mat) = 0;
        }
      }
    }
  }
  {
    //odd internal hexahedra
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,5,4,3,2,1,  6,7,8,9,10,11,  0,1,7,6,  1,2,8,7,  2,3,9,8,  3,4,10,9,   4,5,11,10,  5,0,6,11 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 6, 6, 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,12);
    for(int i = 0; i < n/2-(n%2?0:1); ++i)
    {
      for(int j = 0; j < n; ++j)
      {
        for(int k = 0; k < nz-1; ++k)
        {
          verts[0] = Node(mesh,secondary_nodes[V_ID_S(i*2+1,j,k)]);
          verts[1] = Node(mesh,primary_nodes[V_ID_P(i*2+1,j+1,k)]);
          verts[2] = Node(mesh,primary_nodes[V_ID_P(i*2+2,j+1,k)]);
          verts[3] = Node(mesh,secondary_nodes[V_ID_S(i*2+2,j,k)]);
          verts[4] = Node(mesh,primary_nodes[V_ID_P(i*2+2,j,k)]);
          verts[5] = Node(mesh,primary_nodes[V_ID_P(i*2+1,j,k)]);

          
          verts[6] = Node(mesh,secondary_nodes[V_ID_S(i*2+1,j,k+1)]);
          verts[7] = Node(mesh,primary_nodes[V_ID_P(i*2+1,j+1,k+1)]);
          verts[8] = Node(mesh,primary_nodes[V_ID_P(i*2+2,j+1,k+1)]);
          verts[9] = Node(mesh,secondary_nodes[V_ID_S(i*2+2,j,k+1)]);
          verts[10] = Node(mesh,primary_nodes[V_ID_P(i*2+2,j,k+1)]);
          verts[11] = Node(mesh,primary_nodes[V_ID_P(i*2+1,j,k+1)]);

          Cell c = mesh->CreateCell(verts,nvf,numnodes,8).first;
          c.Integer(mat) = 1;
        }
      }
    }
  }
  {
    //left boundary
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,4,3,2,1,  5,6,7,8,9,  0,1,6,5,  1,2,7,6,  2,3,8,7,  3,4,9,8,   4,0,5,9 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 5, 5, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,10);
    for(int j = 0; j < n; ++j)
    {
      for(int k = 0; k < nz-1; ++k)
      {
        verts[0] = Node(mesh,lb[V_ID_B(j,k)]);
        verts[1] = Node(mesh,primary_nodes[V_ID_P(0,j,k)]);
        verts[2] = Node(mesh,secondary_nodes[V_ID_S(0,j,k)]);
        verts[3] = Node(mesh,primary_nodes[V_ID_P(0,j+1,k)]);
        verts[4] = Node(mesh,lb[V_ID_B(j+1,k)]);

        verts[5] = Node(mesh,lb[V_ID_B(j,k+1)]);
        verts[6] = Node(mesh,primary_nodes[V_ID_P(0,j,k+1)]);
        verts[7] = Node(mesh,secondary_nodes[V_ID_S(0,j,k+1)]);
        verts[8] = Node(mesh,primary_nodes[V_ID_P(0,j+1,k+1)]);
        verts[9] = Node(mesh,lb[V_ID_B(j+1,k+1)]);
        

        Cell c = mesh->CreateCell(verts,nvf,numnodes,7).first;
        c.Integer(mat) = 2;
      }
    }
  }

  {
    //bottom boundary
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,3,2,1,  4,5,6,7,  0,1,5,4,  1,2,6,5,  2,3,7,6,  3,0,4,7 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,8);
    for(int i = 0; i < n/2; ++i)
    {
      for(int k = 0; k < nz-1; ++k)
      {
        verts[0] = Node(mesh,primary_nodes[V_ID_P(i*2,0,k)]);
        verts[1] = Node(mesh,secondary_nodes[V_ID_S(i*2,0,k)]);
        verts[2] = Node(mesh,secondary_nodes[V_ID_S(i*2+1,0,k)]);
        verts[3] = Node(mesh,primary_nodes[V_ID_P(i*2+1,0,k)]);

        verts[4] = Node(mesh,primary_nodes[V_ID_P(i*2,0,k+1)]);
        verts[5] = Node(mesh,secondary_nodes[V_ID_S(i*2,0,k+1)]);
        verts[6] = Node(mesh,secondary_nodes[V_ID_S(i*2+1,0,k+1)]);
        verts[7] = Node(mesh,primary_nodes[V_ID_P(i*2+1,0,k+1)]);

        Cell c = mesh->CreateCell(verts,nvf,numnodes,6).first;
        c.Integer(mat) = 3;
      }
    }
  }

  {
    //top boundary
    const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,3,2,1,  4,5,6,7,  0,1,5,4,  1,2,6,5,  2,3,7,6,  3,0,4,7 };
    const INMOST_DATA_INTEGER_TYPE numnodes[] = { 4, 4, 4, 4, 4, 4 };
    ElementArray<Node> verts(mesh,8);
    for(int i = 0; i < n/2; ++i)
    {
      for(int k = 0; k < nz-1; ++k)
      {
        verts[0] = Node(mesh,primary_nodes[V_ID_P(i*2,n,k)]);
        verts[1] = Node(mesh,secondary_nodes[V_ID_S(i*2,n-1,k)]);
        verts[2] = Node(mesh,secondary_nodes[V_ID_S(i*2+1,n-1,k)]);
        verts[3] = Node(mesh,primary_nodes[V_ID_P(i*2+1,n,k)]);

        verts[4] = Node(mesh,primary_nodes[V_ID_P(i*2,n,k+1)]);
        verts[5] = Node(mesh,secondary_nodes[V_ID_S(i*2,n-1,k+1)]);
        verts[6] = Node(mesh,secondary_nodes[V_ID_S(i*2+1,n-1,k+1)]);
        verts[7] = Node(mesh,primary_nodes[V_ID_P(i*2+1,n,k+1)]);

        Cell c = mesh->CreateCell(verts,nvf,numnodes,6).first;
        c.Integer(mat) = 5;
      }
    }
  }

  if( n%2 == 0 )
  {
    {
      //right boundary
      const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,4,3,2,1,  5,6,7,8,9,  0,1,6,5,  1,2,7,6,  2,3,8,7,  3,4,9,8,   4,0,5,9 };
      const INMOST_DATA_INTEGER_TYPE numnodes[] = { 5, 5, 4, 4, 4, 4, 4 };
      ElementArray<Node> verts(mesh,10);
      for(int j = 0; j < n; ++j)
      {
        for(int k = 0; k < nz-1; ++k)
        {
          verts[0] = Node(mesh,rb[V_ID_B(j,k)]);
          verts[1] = Node(mesh,rb[V_ID_B(j+1,k)]);
          verts[2] = Node(mesh,primary_nodes[V_ID_P(n-1,j+1,k)]);
          verts[3] = Node(mesh,secondary_nodes[V_ID_S(n-1,j,k)]);
          verts[4] = Node(mesh,primary_nodes[V_ID_P(n-1,j,k)]);
        

          verts[5] = Node(mesh,rb[V_ID_B(j,k+1)]);
          verts[6] = Node(mesh,rb[V_ID_B(j+1,k+1)]);
          verts[7] = Node(mesh,primary_nodes[V_ID_P(n-1,j+1,k+1)]);
          verts[8] = Node(mesh,secondary_nodes[V_ID_S(n-1,j,k+1)]);
          verts[9] = Node(mesh,primary_nodes[V_ID_P(n-1,j,k+1)]);
        

          Cell c = mesh->CreateCell(verts,nvf,numnodes,7).first;
          c.Integer(mat) = 4;
        }
      }
    }
  }
  else
  {
    {
      //right boundary
      const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,4,3,2,1,  5,6,7,8,9,  0,1,6,5,  1,2,7,6,  2,3,8,7,  3,4,9,8,   4,0,5,9 };
      const INMOST_DATA_INTEGER_TYPE numnodes[] = { 5, 5, 4, 4, 4, 4, 4 };
      ElementArray<Node> verts(mesh,10);
      for(int j = 0; j < n-1; ++j)
      {
        for(int k = 0; k < nz-1; ++k)
        {
          verts[0] = Node(mesh,rb[V_ID_B(j,k)]);
          verts[1] = Node(mesh,rb[V_ID_B(j+1,k)]);
          verts[2] = Node(mesh,secondary_nodes[V_ID_S(n-1,j+1,k)]);
          verts[3] = Node(mesh,primary_nodes[V_ID_P(n-1,j+1,k)]);
          verts[4] = Node(mesh,secondary_nodes[V_ID_S(n-1,j,k)]);
        

          verts[5] = Node(mesh,rb[V_ID_B(j,k+1)]);
          verts[6] = Node(mesh,rb[V_ID_B(j+1,k+1)]);
          verts[7] = Node(mesh,secondary_nodes[V_ID_S(n-1,j+1,k+1)]);
          verts[8] = Node(mesh,primary_nodes[V_ID_P(n-1,j+1,k+1)]);
          verts[9] = Node(mesh,secondary_nodes[V_ID_S(n-1,j,k+1)]);
        

          Cell c = mesh->CreateCell(verts,nvf,numnodes,7).first;
          c.Integer(mat) = 4;
        }
      }
    }
    {
      //right-bottom corner
      array<HandleType> rbc(nz);
      xyz[0] = 1;
      xyz[1] = 0;
      for(int k = 0; k < nz; ++k)
      {
        xyz[2] = k * 1.0 / (nz-1);
        Node c = mesh->CreateNode(xyz);
        rbc[k] = c.GetHandle();
        c.Integer(mat) = 3;
      }
      const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,3,2,1,  4,5,6,7,  0,1,5,4,  1,2,6,5,  2,3,7,6,  3,0,4,7 };
      const INMOST_DATA_INTEGER_TYPE numnodes[] = { 4, 4, 4, 4, 4, 4 };
      ElementArray<Node> verts(mesh,8);
      for(int k = 0; k < nz-1; ++k)
      {
        verts[0] = Node(mesh,rbc[k]);
        verts[1] = Node(mesh,rb[V_ID_B(0,k)]);
        verts[2] = Node(mesh,secondary_nodes[V_ID_S(n-1,0,k)]);
        verts[3] = Node(mesh,primary_nodes[V_ID_P(n-1,0,k)]);

        verts[4] = Node(mesh,rbc[k+1]);
        verts[5] = Node(mesh,rb[V_ID_B(0,k+1)]);
        verts[6] = Node(mesh,secondary_nodes[V_ID_S(n-1,0,k+1)]);
        verts[7] = Node(mesh,primary_nodes[V_ID_P(n-1,0,k+1)]);

        Cell c = mesh->CreateCell(verts,nvf,numnodes,6).first;
        c.Integer(mat) = 6;
      }
    }
    {
      //right-top corner
      array<HandleType> rtc(nz);
      xyz[0] = 1;
      xyz[1] = 1;
      for(int k = 0; k < nz; ++k)
      {
        xyz[2] = k * 1.0 / (nz-1);
        Node c = mesh->CreateNode(xyz);
        rtc[k] = c.GetHandle();
        c.Integer(mat) = 4;
      }
      const INMOST_DATA_INTEGER_TYPE nvf[] = { 0,3,2,1,  4,5,6,7,  0,1,5,4,  1,2,6,5,  2,3,7,6,  3,0,4,7 };
      const INMOST_DATA_INTEGER_TYPE numnodes[] = { 4, 4, 4, 4, 4, 4 };
      ElementArray<Node> verts(mesh,8);
      for(int k = 0; k < nz-1; ++k)
      {
        verts[0] = Node(mesh,rtc[k]);
        verts[1] = Node(mesh,rb[V_ID_B(n-1,k)]);
        verts[2] = Node(mesh,secondary_nodes[V_ID_S(n-1,n-1,k)]);
        verts[3] = Node(mesh,primary_nodes[V_ID_P(n-1,n,k)]);

        verts[4] = Node(mesh,rtc[k+1]);
        verts[5] = Node(mesh,rb[V_ID_B(n-1,k+1)]);
        verts[6] = Node(mesh,secondary_nodes[V_ID_S(n-1,n-1,k+1)]);
        verts[7] = Node(mesh,primary_nodes[V_ID_P(n-1,n,k+1)]);

        Cell c = mesh->CreateCell(verts,nvf,numnodes,6).first;
        c.Integer(mat) = 7;
      }
    }
  }

  
  int orphaned_nodes = 0;
  for(Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
  {
    if( it->nbAdjElements(CELL) == 0 ) orphaned_nodes++;
  }
  std::cout << "orphan nodes: " << orphaned_nodes << std::endl;
  
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
