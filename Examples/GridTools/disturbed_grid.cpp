#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "inmost.h"

std::string mesh_name = "DISTURBED";

#define SLICE_CELLS

using namespace INMOST;
INMOST_DATA_REAL_TYPE kx = 1000, ky = 10, kz = 10, zone1_angle = 45.0, zone2_angle = 0.0, zone3_angle = 90.0, zone4_angle = 45.0;
#define V_ID(x, y, z) ((x)*n*nz + (y)*nz + (z))

inline INMOST::Storage::real vec_dot_product(INMOST::Storage::real * vecin1, INMOST::Storage::real * vecin2, unsigned int size)
{
  INMOST::Storage::real ret = 0;
  for(unsigned int i = 0; i < size; i++)
    ret += vecin1[i]*vecin2[i];
  return ret;
}

inline void vec_diff(INMOST::Storage::real * vecin1, INMOST::Storage::real * vecin2, INMOST::Storage::real * vecout, unsigned int size)
{
  for(unsigned int i = 0; i < size; i++)
    vecout[i] = vecin1[i] - vecin2[i];
}

Storage::integer prepare_tensor(Storage::real t[9], Storage::real xyz[3])
{
  Storage::integer mat;
  Storage::real k[3] = { kx, ky, kz };
  Storage::real angle = zone2_angle, _sin, _cos, _sincos, _sin2, _cos2;
  if ((xyz[0] + xyz[1] < 0.5) )
  {
    angle = zone1_angle;
    mat = 0;
  }
  else if((xyz[0] + xyz[1] > 1.5))
  {
    angle = zone4_angle;
    mat = 3;
  }
  else if (xyz[0] + xyz[1] > 1)
  {
    angle = zone3_angle;
    mat = 2;
  }
  else mat = 1;
  angle = angle * 3.1415926535897932384626433832795 / 180.0;
  //angle = 0;
  _sin = sin(angle);
  _cos = cos(angle);
  _sincos = _sin*_cos;
  _sin2 = _sin*_sin;
  _cos2 = _cos*_cos;
  t[0] = k[0] * _cos2 + k[1] * _sin2;
  t[1] = (k[0] - k[1]) * _sincos;
  t[2] = 0.0;
  t[3] = t[1];
  t[4] = k[0] * _sin2 + k[1] * _cos2;
  t[5] = 0.0;
  t[6] = 0.0;
  t[7] = 0.0;
  t[8] = k[2];
  return mat;
}

int main(int argc, char *argv[]) 
{
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " nx|ny [nz = 1] [alpha=0.0] [populate_tensor=0] [kx=1000] [ky=10] [kz=10] [zone1_angle=45] [zone2_angle=0] [zone3_angle=90] [zone4_angle=45]" << std::endl;
    return -1;
  }
  Mesh * mesh;
  double alpha=0.0;
  double h;
  int nz, prep = 0;
  int n = 8 + 1;
  nz = 1 + 1;
  if (argc > 1) 
    n = atoi(argv[1])+1;
  if( argc > 2 )
    nz = atoi(argv[2]) +1;
  if( argc > 3 )
    alpha = atof(argv[3]);
  if (argc > 4)
    prep = atoi(argv[4]);

  if (argc > 5) kx = atoi(argv[5]);
  if (argc > 6) ky = atoi(argv[6]);
  if (argc > 7) kz = atoi(argv[7]);
  if (argc > 8) zone1_angle = atoi(argv[8]);
  if (argc > 9) zone2_angle = atoi(argv[9]);
  if (argc > 10) zone3_angle = atoi(argv[10]);
  if (argc > 11) zone4_angle = atoi(argv[11]);

  mesh = new Mesh();

  h = 1.0 / (n-1);

  //srand(time(NULL));
  srand(0);


  int numx = (n-1) / 5;
  int numy = (n-1) / 5;
  numx = numx/2+1;
  numy = numy/2+1;

  MarkerType marker = mesh->CreateMarker();
  int nmarked = 0;
  for (int i = 0; i < n; i++) 
  {
    for (int j = 0; j < n; j++) 
    {
      for (int k = 0; k < nz; k++) 
      {
        Storage::real xyz[3];
        bool mark = false;
        xyz[0] = i * 1.0 / (n - 1);
        xyz[1] = j * 1.0 / (n - 1);
        xyz[2] = k * 1.0 / (nz - 1);
        if ((fabs(xyz[0] + xyz[1] - 0.5) < 1.0e-9 || 
          fabs(xyz[0] + xyz[1] - 1.5) < 1.0e-9 ||
          fabs(xyz[0] + xyz[1] - 1.0) < 1.0e-9) )
        {
          nmarked++;
          mark = true;
        }

        if (fabs(alpha) > 1e-8)
        {
          if((!((i == numx) && (j == numy)) ||
            ((i == numx) && (j == numy-1)) ||
            ((i == numx-1) && (j == numy)) ||
            ((i == numx-1) && (j == numy-1)) ||
            ((i == n - numx) && (j == n - numy)) ||
            ((i == n - numx) && (j == n - numy - 1)) ||
            ((i == n - numx - 1) && (j == n - numy)) ||
            ((i == n - numx - 1) && (j == n - numy - 1))) && !mark)
          {


            if ((i>0) && (i<n-1))  xyz[0] += alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);
            if ((j>0) && (j<n-1))  xyz[1] += alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);
            if ((k>0) && (k<nz-1))  xyz[2] += alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);

          }
        }
        Node c = mesh->CreateNode(xyz);
        if( mark ) c->SetMarker(marker);
        if (c->LocalID() != V_ID(i, j, k)) printf("v_id = %d, [i,j,k] = %d\n", c->LocalID(), V_ID(i, j, k));
      }
    }
  }


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

  //slice cells
#if defined(SLICE_CELLS)
  mesh->BeginModification();
  for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it) if(!it->New())
  {
    int numnodes = 0;
    ElementArray<Node> nodes = it->getNodes();
    ElementArray<Node> quad_nodes(mesh,4);
    ElementArray<Face> f(mesh,1);
    for(ElementArray<Node>::iterator jt = nodes.begin(); jt != nodes.end(); ++jt)
    {
      if( jt->GetMarker(marker) )
      {
        quad_nodes[numnodes++] = jt->self();
      }
    }
    if(numnodes == 4) 
    {

      ElementArray<Face> face_sets[4];
      for(int i = 0; i < 4; i++)
      {
        face_sets[i].SetMeshLink(mesh);
        face_sets[i] = quad_nodes[i]->getFaces();
      }
      bool done[4] = {false,false,false,false};
      for(int i = 0; i < 4; i++) if( !done[i] )
      {
        for(int j = i+1; j < 4; j++) if( !done[j] )
        {
          ElementArray<Face> test = face_sets[i];
          test.Intersect(face_sets[j]);
          if( test.size() == 1)
          {
            ElementArray<Edge> e(mesh,1);
            ElementArray<Node> edge_nodes(mesh,2);
            edge_nodes[0] = quad_nodes[i];
            edge_nodes[1] = quad_nodes[j];
            done[i] = true;
            done[j] = true;
            e[0] = mesh->CreateEdge(edge_nodes).first;
            Face::SplitFace(test[0],e,0);
            break;
          }
        }
      }
      /*
      Edge * e; 
      adjacent<Face> faces;
      e = mesh->CreateEdge(quad_nodes,2).first;
      faces = quad_nodes[0]->getFaces();
      faces.intersect(quad_nodes[1]->getFaces());
      assert(faces.size()==1);
      Face::SplitFace(&faces[0],&e,1,0);
      e = mesh->CreateEdge(quad_nodes+2,2).first;
      faces = quad_nodes[2]->getFaces();
      faces.intersect(quad_nodes[3]->getFaces());
      assert(faces.size()==1);
      Face::SplitFace(&faces[0],&e,1,0);

      */

      Node temp = quad_nodes[2];
      quad_nodes[2] = quad_nodes[3];
      quad_nodes[3] = temp;

      f[0] = mesh->CreateFace(quad_nodes).first;


      Cell::SplitCell(it->self(),f,0);

    }
  }
	mesh->ApplyModification();
  mesh->EndModification();
#endif
  printf("nodes: %d edges: %d faces: %d cells: %d\n", mesh->NumberOfNodes(), mesh->NumberOfEdges(), mesh->NumberOfFaces(), mesh->NumberOfCells());

  if (prep)
  {
    std::cout << "put tensor" << std::endl;
    Tag mat = mesh->CreateTag("MATERIAL",DATA_INTEGER,CELL,NONE,1);
    Storage::real cnt[3];
    Tag tensor = mesh->CreateTag("K", DATA_REAL, CELL, NONE, 9);
    for (Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
    {
      it->Centroid(cnt);
      it->IntegerDF(mat) = prepare_tensor(&it->RealArrayDF(tensor)[0], cnt);
    }
  }

  Storage::bulk_array name = mesh->self()->BulkArray(mesh->CreateTag("GRIDNAME",DATA_BULK,MESH,NONE));
  name.replace(name.begin(),name.end(),mesh_name.begin(),mesh_name.end());

  printf("I'm ready!\n");

  mesh->Save("grid.vtk");
  mesh->Save("grid.pmf");
  mesh->Save("grid.gmv");

  printf("File written!\n");

  delete mesh;
  return 0;
}
