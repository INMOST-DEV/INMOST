#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


std::string mesh_name = "SINUSOIDAL";

using namespace INMOST;
const Storage::real pi = 3.1415926535897932384626433832795;
#define V_ID(x, y, z) ((x)*n*nz + (y)*nz + (z))
#define RZIND(x,y) ((x)*n + (y))

int main(int argc, char *argv[]) 
{
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " n [nz = 1] [eta=1/(4pi)]" << std::endl;
    return -1;
  }
  Mesh * mesh;
  double eta = 1.0/(2.0*pi), L = 1;
  int nz = 2;
  int n = 16;
  if (argc > 1) 
    n = atoi(argv[1])+1;

  if( argc > 2 )
    nz = atoi(argv[2])+1;
  if( argc > 3 )
    eta = atof(argv[3]);

  if( eta > 1.0/(2*pi) )
  {
    printf("Eta %g is too large, should be 0 <= eta <= %g\n",eta,1.0/(2.0*pi));
  }

  double * pr = new double[n*n*2];
  double * pz = pr+n*n;
  memset(pr,0,sizeof(double)*n*n*2);
  shell<double> r(pr,n*n);
  shell<double> z(pz,n*n);
  
  for(int i = 0; i < n; ++i)
  {
    for(int j = 0; j < n; ++j)
    {
      double x = i * 1.0 / (n - 1);
      double y = j * 1.0 / (n - 1);

      r[RZIND(i,j)] = x + eta * sin(2*pi*x) * sin(2*pi*y);
      z[RZIND(i,j)] = y + eta * sin(2*pi*x) * sin(2*pi*y);
    }
  }
  
  mesh = new Mesh();
  

  for (int i = 0; i < n; i++) 
  {
    for (int j = 0; j < n; j++) 
    {
      for (int k = 0; k < nz; k++) 
      {
        Storage::real xyz[3];
        xyz[0] = r[i*n+j];
        xyz[1] = z[i*n+j];
        xyz[2] = k * 1.0 / (nz - 1);
        Node c = mesh->CreateNode(xyz);
        if (c->LocalID() != V_ID(i, j, k)) printf("v_id = %d, [i,j,k] = %d\n", c->LocalID(), V_ID(i, j, k));
      }
    }
  }

  delete [] pr;


  const INMOST_DATA_INTEGER_TYPE nvf[24] = { 2, 3, 1, 0, 4, 5, 7, 6, 0, 1, 5, 4, 3, 2, 6, 7, 2, 0, 4, 6, 1, 3, 7, 5 };
  const INMOST_DATA_INTEGER_TYPE numnodes[6] = { 4, 4, 4, 4, 4, 4 };
  for (int i = 1; i < n; i++) 
  {
    for (int j = 1; j < n; j++) 
    {
      for (int k = 1; k < nz; k++) 
      {

        ElementArray<Node> verts(mesh,8); 
        verts[0] = mesh->NodeByLocalID(V_ID(i - 1, j - 1, k - 1));
        verts[1] = mesh->NodeByLocalID(V_ID(i - 0, j - 1, k - 1));
        verts[2] = mesh->NodeByLocalID(V_ID(i - 1, j - 0, k - 1));
        verts[3] = mesh->NodeByLocalID(V_ID(i - 0, j - 0, k - 1));
        verts[4] = mesh->NodeByLocalID(V_ID(i - 1, j - 1, k - 0));
        verts[5] = mesh->NodeByLocalID(V_ID(i - 0, j - 1, k - 0));
        verts[6] = mesh->NodeByLocalID(V_ID(i - 1, j - 0, k - 0));
        verts[7] = mesh->NodeByLocalID(V_ID(i - 0, j - 0, k - 0));

        mesh->CreateCell(verts,nvf,numnodes,6);
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
