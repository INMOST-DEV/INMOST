#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


std::string mesh_name = "TRANSFORMED";

using namespace INMOST;



void normalize(Storage::real nrm[3])
{
	double l = 0.0;
	for(int i = 0; i < 3; ++i) l += nrm[i]*nrm[i];
	l = sqrt(l);
	if( l ) for(int i = 0; i < 3; ++i) nrm[i] /= l;
}

const Storage::real pi = 3.14;

void transform(Storage::real xyz[3], Storage::real & ztop, Storage::real & zbottom, Storage::real nrmtop[3], Storage::real nrmbottom[3])
{
	unknown x(xyz[0],0), y = unknown(xyz[1],1);
	//zbottom = x*x*4-y*y*4-sin(6*x)*8 - cos(4*y)*4 - x*15;
	//ztop = x*x*4-y*y*4-sin(6*x)*8 - cos(4*y)*4 - y*15 + 15;
	variable zb = sin(y*pi *2)*0.2;// ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5))*0.4;
	variable zt = sin(y*pi *2)*0.2+1;// 1 + cos(4 * y)*0.2 + ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5))*0.4 + 1 + cos(4 * y)*0.5;
	
	//variable zb = 0;// ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))*0.4;
	//variable zt = 1 + 1*y;//cos(4*y)*0.2;// + ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))*0.4 + 1 + cos(4*y)*0.5;
	
	//variable zb = 0;// ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))*0.4;
	//variable zt = 1 + cos(4*y)*0.2;// + ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))*0.4 + 1 + cos(4*y)*0.5;
	
	
	zbottom = zb.GetValue();
	ztop = zt.GetValue();
	nrmtop[0] = zt.GetRow()[0];
	nrmtop[1] = zt.GetRow()[1];
	nrmtop[2] = -1;
	normalize(nrmtop);
	nrmbottom[0] = zb.GetRow()[0];
	nrmbottom[1] = zb.GetRow()[1];
	nrmbottom[2] = -1;
	normalize(nrmbottom);
}

struct quat
{
	Storage::real vec[3];
	Storage::real w;
};

Storage::real quatnorm(const quat & a)
{
	return a.vec[0]*a.vec[0] + a.vec[1]*a.vec[1]
		 + a.vec[2]*a.vec[2] + a.w*a.w;
}

quat quatconj(const quat & a)
{
	quat ret = a;
	ret.vec[0] = -ret.vec[0];
	ret.vec[1] = -ret.vec[1];
	ret.vec[2] = -ret.vec[2];
	return ret;
}

void rotate_tensor(Storage::real nrm[3], const rMatrix & Kin, rMatrix & Kout, Storage::real phi[3])
{
	//~ double qmat[16], qrmat[16], sclmat[16];
	Storage::real qnrm;
	struct quat q, qr;
	rMatrix N = rMatrix::FromVector(nrm,3);
	rMatrix NX = rMatrix::CrossProductMatrix(nrm);
	rMatrix U(3,1), NXU;
	U(0,0) = 0;
	U(1,0) = 0;
	U(2,0) = -1.0;
	NXU = NX*U;
	
	q.vec[0] = NXU(0,0);
	q.vec[1] = NXU(1,0);
	q.vec[2] = NXU(2,0);
	q.w = 1.0 + N.DotProduct(U);

	qnrm = quatnorm(q);
	q.vec[0] /= qnrm;
	q.vec[1] /= qnrm;
	q.w /= qnrm;
	qr = quatconj(q);

	rMatrix Q(3,3);
	
	Q(0,0) = (q.w*q.w + q.vec[0]*q.vec[0] - q.vec[1]*q.vec[1]) *qnrm;
	Q(0,1) = 2.*(q.vec[0]*q.vec[1] ) *qnrm                          ;
	Q(0,2) = 2.*(q.w*q.vec[1]) *qnrm                                ;
	
	Q(1,0) = 2.*(q.vec[0]*q.vec[1]) *qnrm                           ;
	Q(1,1) = (q.w*q.w - q.vec[0]*q.vec[0] + q.vec[1]*q.vec[1]) *qnrm;
	Q(1,2) = 2.*( - q.w*q.vec[0])*qnrm                              ;
	
	Q(2,0) = 2.*( - q.w*q.vec[1])*qnrm                              ;
	Q(2,1) = 2.*( + q.w*q.vec[0])*qnrm                              ;
	Q(2,2) = (q.w*q.w - q.vec[0]*q.vec[0] - q.vec[1]*q.vec[1]) *qnrm;


	Kout = Q.Transpose()*Kin*Q;

	Storage::real t0 = -2.0 * (q.vec[1] * q.vec[1] + q.vec[2] * q.vec[2]) + 1.0;
	Storage::real t1 = +2.0 * (q.vec[0] * q.vec[1] - q.w * q.vec[2]);
	Storage::real t2 = -2.0 * (q.vec[0] * q.vec[2] + q.w * q.vec[1]);
	Storage::real t3 = +2.0 * (q.vec[1] * q.vec[2] - q.w * q.vec[0]);
	Storage::real t4 = -2.0 * (q.vec[0] * q.vec[0] + q.vec[1] * q.vec[1]) + 1.0;

	t2 = t2 > 1.0 ? 1.0 : t2;
	t2 = t2 < -1.0 ? -1.0 : t2;

	phi[0] = std::atan2(t3, t4) * 180 / 3.1415926535897932384626433832795;
	phi[1] = std::asin(t2)*180/3.1415926535897932384626433832795;
	phi[2] = std::atan2(t1, t0)*180/3.1415926535897932384626433832795;

}


int main(int argc, char *argv[]) 
{
  Mesh * mesh = new Mesh();
  //~ int prep = 0;
  if( argc < 2 )
  {
    std::cout << "Usage: " << argv[0] << " mesh_file [input_perm]" << std::endl;
    return -1;
  }

  mesh->Load(argv[1]);

  Tag nrm_tag = mesh->CreateTag("NORMAL",DATA_REAL,CELL,NONE,3);
  Tag degreez_tag = mesh->CreateTag("DEGREEZ",DATA_REAL,CELL,NONE,1);
  Tag degreex_tag = mesh->CreateTag("DEGREEX",DATA_REAL,CELL,NONE,1);
  Tag degreey_tag = mesh->CreateTag("DEGREEY",DATA_REAL,CELL,NONE,1);
  Tag perm_tag;
  if( argc > 2 )
  {
	  if( mesh->HaveTag(argv[2]) )
		perm_tag = mesh->GetTag(argv[2]);
	  else std::cout << "Do not have tag " << argv[2] << " on mesh" << std::endl;
  }
  if( !perm_tag.isValid() )
  {
	  if (mesh->HaveTag("PERM")) mesh->DeleteTag(mesh->GetTag("PERM"));
	  perm_tag = mesh->CreateTag("PERM",DATA_REAL,CELL,NONE,3);
	  for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
	  {
		  Storage::real_array perm = it->RealArray(perm_tag);
		  perm[0] = 1000;// + 50*rand()/(1.*RAND_MAX);
		  perm[1] = 1000;// + 50*rand()/(1.*RAND_MAX);
		  perm[2] = 10;// + 50*rand()/(1.*RAND_MAX);
	  }
  }
  Tag outperm_tag = mesh->CreateTag("PERM_TRANSFORM",DATA_REAL,CELL,NONE,9);

  Storage::real min[3] = {1.0e+20,1.0e+20,1.0e+20}, max[3] = {-1.0e+20,-1.0e+20,-1.0e+20};

  for(Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
  {
	  Storage::real_array c = it->Coords();
	  for(int i = 0; i < 3; ++i)
	  {
		  if( c[i] < min[i] ) min[i] = c[i];
		  if( c[i] > max[i] ) max[i] = c[i];
	  }
  }

  for(Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
  {
	  Storage::real_array c = it->Coords();
	  for(int i = 0; i < 3; ++i)
		  c[i] = (c[i]-min[i])/(max[i]-min[i]);
  }


  for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
  {
	  Storage::real cnt[3], zb, zt, nrmb[3], nrmt[3], phi[3];
	  it->Centroid(cnt);
	  transform(cnt,zt,zb,nrmt,nrmb);
	  Storage::real_array nrm = it->RealArray(nrm_tag);
	  for(int i = 0; i < 3; ++i)
		nrm[i] = (nrmt[i]-nrmb[i])*cnt[2] + nrmb[i];
	  normalize(nrm.data());


	  
	  
	  rMatrix Kperm = rMatrix::FromTensor(it->RealArray(perm_tag).data(),it->RealArray(perm_tag).size());
	  rMatrix Kout;
	  rotate_tensor(nrm.data(),Kperm,Kout,phi);

	  it->Real(degreex_tag) = phi[0];
	  it->Real(degreey_tag) = phi[1];
	  it->Real(degreez_tag) = phi[2];

	  Storage::real_array K = it->RealArray(outperm_tag);
	  std::copy(Kout.data(),Kout.data()+9,K.data());
  }

  for(Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
  {
	  Storage::real ztop, zbottom, nrmb[3],nrmt[3];
	  Storage::real_array c = it->Coords();
	  transform(c.data(), ztop, zbottom,nrmt,nrmb);
	  c[2] = (ztop-zbottom)*c[2] + zbottom;
  }

  for(Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
  {
	  Storage::real_array c = it->Coords();
	  for(int i = 0; i < 3; ++i)
		  c[i] = c[i]*(max[i]-min[i]) + min[i];
  }

  Storage::bulk_array name = mesh->self()->BulkArray(mesh->CreateTag("GRIDNAME",DATA_BULK,MESH,NONE));
  name.replace(name.begin(),name.end(),mesh_name.begin(),mesh_name.end());

  std::cout << "nodes: " << mesh->NumberOfNodes() << " edges: " << mesh->NumberOfEdges() << " faces: " << mesh->NumberOfFaces() << " cells: " << mesh->NumberOfCells() << std::endl;
  std::cout << "I'm ready!" << std::endl;
  mesh->Save("grid_out.xml");
  
  delete mesh;
  return 0;
}
