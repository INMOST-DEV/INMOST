#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


std::string mesh_name = "KERSHAW";

using namespace INMOST;
const Storage::real eps = 1.0e-9;
#define V_ID(x, y, z) ((x)*n*nz + (y)*nz + (z))
#define RZIND(x,y) ((x)*n + (y))

int main(int argc, char *argv[]) 
{
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " n [nz = 1] [alpha=0.5]" << std::endl;
    return -1;
  }
  Mesh * mesh;
    double alpha = 0.1;//, L = 1;//, rm = L, zm = L;
  int nz = 2;
  int n = 20;
  if (argc > 1) 
    n = atoi(argv[1])+1;
  if( argc > 2 )
    nz = atoi(argv[2])+1;
  if( argc > 3 )
    alpha = atof(argv[3]);

  if( alpha < eps || alpha > 1-eps )
  {
    std::cout << "Bad alpha " << alpha << ", should be " << eps << " < alpha < 1-" << eps << "\n";
    return -1;
  }

  int nsub[5], nmult;

  if( (n-1)%20 == 0 )
  {
    nsub[0] = 4;
    nsub[1] = 7;
    nsub[2] = 13;
    nsub[3] = 16;
    nsub[4] = 20;
    nmult = n/20;
  }
  else
  {
    std::cout << "n should be a multiply of 20 but n is " << n << "\n";
    return -1;
  }

  double slope1 = (1-2*alpha)/0.15;
  double slope2 = (2*alpha-1)/0.30;

  for(int i = 0; i < 5; ++i) nsub[i] = nsub[i]*nmult;

  double * pr = new double[n*n*2];
  double * pz = pr+n*n;
  memset(pr,0,sizeof(double)*n*n*2);
  shell<double> r(pr,n*n);
  shell<double> z(pz,n*n);
  
  for(int i = 0; i < n; ++i)
  {
    double x = i * 1.0 / (n - 1);
    double ymid = 0;
    if( i <= nsub[0] )
      ymid = alpha;
    else if( i <= nsub[1] )
      ymid = slope1 * (x - 0.2) + alpha;
    else if( i <= nsub[2] )
      ymid = slope2 * (x - 0.35) + (1-alpha);
    else if( i <= nsub[3] )
      ymid = slope1 * (x - 0.65) + alpha;
    else if( i <= nsub[4] )
      ymid = 1-alpha;
    else std::cout << __FILE__ << ":" << __LINE__ << " oops!\n";

    for(int j = 0; j < n/2; ++j)
    {
      double y = j * 1.0 / (n - 1); //from 0 to 0.5

      r[RZIND(i,j)] = x;
      z[RZIND(i,j)] = ymid * y * 2;
    }

    for(int j = n/2; j < n; ++j)
    {
      double y = j * 1.0 / (n - 1); //from 0.5 to 1.0

      r[RZIND(i,j)] = x;
      z[RZIND(i,j)] = ymid + (1.0-ymid) * (y-0.5) * 2;
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
        if (c->LocalID() != V_ID(i, j, k)) std::cout << "v_id = " << c->LocalID() << ", [i,j,k] = " << V_ID(i,j,k) << "\n";
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

  std::cout << "I'm ready!\n";

  mesh->Save("grid.vtk");
  mesh->Save("grid.pmf");
  mesh->Save("grid.gmv");

  std::cout << "File written!\n";

  delete mesh;
  return 0;
}
