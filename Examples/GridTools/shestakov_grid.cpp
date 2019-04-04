#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

std::string mesh_name = "SHESTAKOV";

using namespace INMOST;
#define V_ID(x, y, z) ((x)*nm*nz + (y)*nz + (z))
#define RZIND(x,y) (((x)-1)*nm + ((y)-1))

int main(int argc, char *argv[]) 
{
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " nc [nz = 1] [alpha=0.25]" << std::endl;
    return -1;
  }
  Mesh * mesh;
  double alpha=0.35, L = 1;
  int nz = 2;
  int nc = 4;
  if (argc > 1) 
    nc = atoi(argv[1])+1;
  if( argc > 2 )
    nz = atoi(argv[2])+1;
  if( argc > 3 )
    alpha = atof(argv[3]);

  if( nc >= 10 )
  {
    printf("You are about to create a grid of size %dx%dx%d. Are you sure?\n",static_cast<int>(pow(2.0,nc)),static_cast<int>(pow(2.0,nc)),nz);
    char ans = getchar();
    if( towlower(ans) != 'y' ) return 0;
  }
  srand(0);
  int nm = static_cast<int>(pow(2.0,nc))+1;
  int ind = nm-1;

  double * pr = new double[nm*nm*2];
  double * pz = pr+nm*nm;
  memset(pr,0,sizeof(double)*nm*nm*2);
  shell<double> r(pr,nm*nm);
  shell<double> z(pz,nm*nm);

  for(int i = 0; i <= 1; ++i)
  {
    int k = 1 + i*ind;
    for(int j = 0; j <= 1; ++j)
    {
      int l = 1 + j*ind;
      r[RZIND(k,l)] = L*(1-i);
      z[RZIND(k,l)] = L*j;
    }
  }

  for(int nl = 0; nl <= nc-1; ++nl)
  {
    int nn = static_cast<int>(pow(2.0,nl));
    int inc = ind/nn;
    int inch = inc/2;
    for(int k = 1; k <= nn; ++k)
    {
      int k1 = 1 + (k-1)*inc;
      int k2 = k1 + inc;
      int k3 = k1 + inch;
      for(int l = 1; l <= nn; ++l)
      {
        int l1 = 1 + (l-1)*inc;
        int l2 = l1 + inc;
        int l3 = l1 + inch;
        double ar,br,r1,r2,r3,r4,z1,z2,z3,z4;
        if( l == 1 )
        {
          ar = alpha + static_cast<double>(rand())/static_cast<double>(RAND_MAX)*(1-2*alpha);
          r[RZIND(k3,l1)] = ar*r[RZIND(k1,l1)] + (1-ar)*r[RZIND(k2,l1)];
          z[RZIND(k3,l1)] = ar*z[RZIND(k1,l1)] + (1-ar)*z[RZIND(k2,l1)];
        }
        ar = alpha + static_cast<double>(rand())/static_cast<double>(RAND_MAX)*(1-2*alpha);
        r[RZIND(k2,l3)] = ar*r[RZIND(k2,l1)] + (1-ar)*r[RZIND(k2,l2)];
        z[RZIND(k2,l3)] = ar*z[RZIND(k2,l1)] + (1-ar)*z[RZIND(k2,l2)];
        ar = alpha + static_cast<double>(rand())/static_cast<double>(RAND_MAX)*(1-2*alpha);
        r[RZIND(k3,l2)] = ar*r[RZIND(k1,l2)] + (1-ar)*r[RZIND(k2,l2)];
        z[RZIND(k3,l2)] = ar*z[RZIND(k1,l2)] + (1-ar)*z[RZIND(k2,l2)];
        if( k == 1 )
        {
          ar = alpha + static_cast<double>(rand())/static_cast<double>(RAND_MAX)*(1-2*alpha);
          r[RZIND(k1,l3)] = ar*r[RZIND(k1,l1)] + (1-ar)*r[RZIND(k1,l2)];
          z[RZIND(k1,l3)] = ar*z[RZIND(k1,l1)] + (1-ar)*z[RZIND(k1,l2)];
        }
        ar = alpha + static_cast<double>(rand())/static_cast<double>(RAND_MAX)*(1-2*alpha);
        br = alpha + static_cast<double>(rand())/static_cast<double>(RAND_MAX)*(1-2*alpha);
        r1 = r[RZIND(k1,l1)];
        r2 = r[RZIND(k2,l1)];
        r3 = r[RZIND(k2,l2)];
        r4 = r[RZIND(k1,l2)];
        z1 = z[RZIND(k1,l1)];
        z2 = z[RZIND(k2,l1)];
        z3 = z[RZIND(k2,l2)];
        z4 = z[RZIND(k1,l2)];
        // check for boomerang zones
        double det2,det3,det4,det1,d,r3p,z3p,r1p,z1p,r4p,z4p,r2p,z2p;
        det2 = (r2-r1)*(z3-z2)-(r3-r2)*(z2-z1);
        det3 = (r3-r2)*(z4-z3)-(r4-r3)*(z3-z2);
        det4 = (r4-r3)*(z1-z4)-(r1-r4)*(z4-z3);
        det1 = (r1-r4)*(z2-z1)-(r2-r1)*(z1-z4);
        if (det2>0)
        {
          d = (r4-r3)*(z2-z1)-(r2-r1)*(z4-z3);
          r3p = ((r2-r1)*(r4*z3-r3*z4)-(r4-r3)*(r2*z1-r1*z2))/d;
          z3p = ((z2-z1)*(r4*z3-r3*z4)-(z4-z3)*(r2*z1-r1*z2))/d;
          d = (r4-r1)*(z2-z3)-(r2-r3)*(z4-z1);
          r1p = ((r2-r3)*(r4*z1-r1*z4)-(r4-r1)*(r2*z3-r3*z2))/d;
          z1p = ((z2-z3)*(r4*z1-r1*z4)-(z4-z1)*(r2*z3-r3*z2))/d;
          r3 = r3p;
          z3 = z3p;
          r1 = r1p;
          z1 = z1p;
        }
        else if (det3>0)
        {
          d = (r1-r4)*(z3-z2)-(r3-r2)*(z1-z4);
          r4p = ((r3-r2)*(r1*z4-r4*z1)-(r1-r4)*(r3*z2-r2*z3))/d;
          z4p = ((z3-z2)*(r1*z4-r4*z1)-(z1-z4)*(r3*z2-r2*z3))/d;
          d = (r1-r2)*(z3-z4)-(r3-r4)*(z1-z2);
          r2p = ((r3-r4)*(r1*z2-r2*z1)-(r1-r2)*(r3*z4-r4*z3))/d;
          z2p = ((z3-z4)*(r1*z2-r2*z1)-(z1-z2)*(r3*z4-r4*z3))/d;
          r4 = r4p;
          z4 = z4p;
          r2 = r2p;
          z2 = z2p;
        }
        else if (det4>0)
        {
          d = (r2-r1)*(z4-z3)-(r4-r3)*(z2-z1);
          r1p = ((r4-r3)*(r2*z1-r1*z2)-(r2-r1)*(r4*z3-r3*z4))/d;
          z1p = ((z4-z3)*(r2*z1-r1*z2)-(z2-z1)*(r4*z3-r3*z4))/d;
          d = (r2-r3)*(z4-z1)-(r4-r1)*(z2-z3);
          r3p = ((r4-r1)*(r2*z3-r3*z2)-(r2-r3)*(r4*z1-r1*z4))/d;
          z3p = ((z4-z1)*(r2*z3-r3*z2)-(z2-z3)*(r4*z1-r1*z4))/d;
          r1 = r1p;
          z1 = z1p;
          r3 = r3p;
          z3 = z3p;
        }
        else if (det1>0)
        {
          d = (r3-r2)*(z1-z4)-(r1-r4)*(z3-z2);
          r2p = ((r1-r4)*(r3*z2-r2*z3)-(r3-r2)*(r1*z4-r4*z1))/d;
          z2p = ((z1-z4)*(r3*z2-r2*z3)-(z3-z2)*(r1*z4-r4*z1))/d;
          d = (r3-r4)*(z1-z2)-(r1-r2)*(z3-z4);
          r4p = ((r1-r2)*(r3*z4-r4*z3)-(r3-r4)*(r1*z2-r2*z1))/d;
          z4p = ((z1-z2)*(r3*z4-r4*z3)-(z3-z4)*(r1*z2-r2*z1))/d;
          r2 = r2p;
          z2 = z2p;
          r4 = r4p;
          z4 = z4p;
        }
        r[RZIND(k3,l3)] = ar*br*r1 + ar*(1-br)*r4 + (1-ar)*br*r2 + (1-ar)*(1-br)*r3;
        z[RZIND(k3,l3)] = ar*br*z1 + ar*(1-br)*z4 + (1-ar)*br*z2 + (1-ar)*(1-br)*z3;
      }
    }
  }

  
  mesh = new Mesh();
  

  for (int i = 0; i < nm; i++) 
  {
    for (int j = 0; j < nm; j++) 
    {
      for (int k = 0; k < nz; k++) 
      {
        Storage::real xyz[3];
        xyz[0] = r[i*nm+j];
        xyz[1] = z[i*nm+j];
        xyz[2] = k * 1.0 / (nz - 1);
        Node c = mesh->CreateNode(xyz);
        if (c->LocalID() != V_ID(i, j, k)) printf("v_id = %d, [i,j,k] = %d\n", c->LocalID(), V_ID(i, j, k));
      }
    }
  }

  delete [] pr;


  const INMOST_DATA_INTEGER_TYPE nvf[24] = { 2, 3, 1, 0, 4, 5, 7, 6, 0, 1, 5, 4, 3, 2, 6, 7, 2, 0, 4, 6, 1, 3, 7, 5 };
  const INMOST_DATA_INTEGER_TYPE numnodes[6] = { 4, 4, 4, 4, 4, 4 };
  for (int i = 1; i < nm; i++) 
  {
    for (int j = 1; j < nm; j++) 
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
  str << mesh_name << "_" << pow(2.0,nc) << "x" << pow(2.0,nc) << "x" << nz;
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
