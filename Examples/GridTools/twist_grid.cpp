#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


std::string mesh_name = "TWIST-GRID";

using namespace INMOST;
const Storage::real pi = 3.1415926535897932384626433832795;
#define V_ID(x, y, z) ((x)*n*nz + (y)*nz + (z))
#define V_ID3(x, y, z, m) (((x)*(n-1)*nz + (y)*nz + (z))*(np*4+1) + m)

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
	
  int np = 10;
  double w = 0.5/(double)(np+1);
  
  mesh = new Mesh();
  
  array<HandleType> primary_nodes(n*n*nz);
  array<HandleType> internal_nodes((n-1)*(n-1)*nz*(np*4+1)); // n-1 by x, n-1 by y
  

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


  //go through cells
  for(int i = 0; i < n-1; i++)
  {
    for(int j = 0; j < n-1; j++)
    {
      for (int k = 0; k < nz; k++) 
      {
		  Node c;
		  Storage::real xyz[3];
		  xyz[0] = (i+0.5) * 1.0 / (n-1);
		  xyz[1] = (j+0.5) * 1.0 / (n-1);
		  xyz[2] = k * 1.0 / (nz - 1);
		  c = mesh->CreateNode(xyz);
		  c.Integer(mat) = 1;
		  internal_nodes[V_ID3(i,j,k,0)] = c.GetHandle();
		  for(int m = 0; m < 4; ++m)
		  {
			  int pp[4][2] = {{0,0},{1,0},{1,1},{0,1}};
			  for(int l = 0; l < np; ++l)
			  {
				 // std::cout << "l " << l << " m " << m << " x " << (i + pp[m][0] + (l+1)*w*(1-2*pp[m][0]) ) * 1.0 / (n-1) << " y " << (j + pp[m][1] + (l+1)*w*(1-2*pp[m][1]) ) * 1.0 / (n-1) << std::endl;
				xyz[0] = (i + pp[m][0] + (l+1)*w*(1-2*pp[m][0]) ) * 1.0 / (n-1);
				xyz[1] = (j + pp[m][1] + (l+1)*w*(1-2*pp[m][1]) ) * 1.0 / (n-1);
				xyz[2] = k * 1.0 / (nz - 1);
				c = mesh->CreateNode(xyz);
				c.Integer(mat) = m*np+l+1;
				internal_nodes[V_ID3(i,j,k,m*np+l+1)] = c.GetHandle();
			  }
		  }
      }
    }
  }
  

  //internals of the mesh
  {
	  std::vector<INMOST_DATA_INTEGER_TYPE> nvf1( (3+2*np)*2 + (3+2*np)*4 );
	  std::vector<INMOST_DATA_INTEGER_TYPE> numnodes1(3+2*np + 2);
	  ElementArray<Node> verts(mesh,(3+2*np)*2);
	  numnodes1[0] = 3+2*np; //bottom
	  numnodes1[1] = 3+2*np; //top
	  for(int k = 0; k < 3 + 2*np; ++k)
	  {
		  numnodes1[k+2] = 4; // sides
		  
		  nvf1[k] = 3 + 2*np - 1 - k; //bottom reverse
		  nvf1[k + 3 + 2*np] = k + 3 + 2*np; //top forward
		  
		  nvf1[k*4 + 0 + (3+2*np)*2] = k;
		  nvf1[k*4 + 1 + (3+2*np)*2] = (k+1)%(3+2*np);
		  nvf1[k*4 + 2 + (3+2*np)*2] = 3+2*np + (k+1)%(3+2*np);
		  nvf1[k*4 + 3 + (3+2*np)*2] = 3+2*np + k;
	  }
	  
	  
	  {
		  for(int i = 0; i < n-1; ++i)
		  {
			  for(int j = 0; j < n-1; ++j)
			  {
				  for(int k = 0; k < nz-1; ++k)
				  {
					  
					  Cell c;
					  
					  for(int m = 0; m < 4; ++m)
					  {
						  int pp[4][2] = {{0,0},{1,0},{1,1},{0,1}};
						  
						  
						  verts[0] = Node(mesh,primary_nodes[V_ID(i+pp[(m+0)%4][0],j+pp[(m+0)%4][1],k)]);
						  verts[1] = Node(mesh,primary_nodes[V_ID(i+pp[(m+1)%4][0],j+pp[(m+1)%4][1],k)]);
						  
						  verts[(3+2*np) + 0] = Node(mesh,primary_nodes[V_ID(i+pp[(m+0)%4][0],j+pp[(m+0)%4][1],k+1)]);
						  verts[(3+2*np) + 1] = Node(mesh,primary_nodes[V_ID(i+pp[(m+1)%4][0],j+pp[(m+1)%4][1],k+1)]);
						  
						  for(int l = 0; l < np; ++l)
						  {
							  verts[2+l] = Node(mesh,internal_nodes[V_ID3(i,j,k,1 + ((m+2+l)%4)*np + l)]);
							  verts[(3+2*np) + 2+l] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,1 + ((m+2+l)%4)*np + l)]);
						  }
						  
						  verts[2+np] = Node(mesh,internal_nodes[V_ID3(i,j,k,0)]);
						  verts[(3+2*np) + 2+np] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,0)]);
						  
						  for(int l = 0; l < np; ++l)
						  {
							  verts[3+np+l] = Node(mesh,internal_nodes[V_ID3(i,j,k,1 + ((m+np-l)%4)*np + (np-l-1))]);
							  verts[(3+2*np) + 3+np+l] = Node(mesh,internal_nodes[V_ID3(i,j,k+1,1 + ((m+np-l)%4)*np + (np-l-1))]);
						  }
						  
						  c = mesh->CreateCell(verts,&nvf1[0],&numnodes1[0],3+2*np + 2).first;
						  c.Integer(mat) = m;
					  }
					  
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
