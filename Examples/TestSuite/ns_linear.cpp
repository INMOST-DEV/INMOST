#include "inmost.h"

using namespace INMOST;
//setup on mesh:
//
// FORCE - 3 entries
//
// BOUNDARY_CONDITION_VELOCITY - 7 entries
// r = (a5,a6,a7)
// n.(a1*u + a2*t) = n.r
// (I - nn.)(a3*u + a4*t) = (I-nn.)r
//
// BOUNDARY_CONDITION_PRESSURE - 1 entry
//
// REFERENCE_VELOCITY - 3 entries
//
// REFERENCE_PRESSURE - 1 entry
//

typedef Storage::real real;


vMatrix refUVW(real vx, real vy, real vz, real t, real nu)
{
	(void) t,(void)nu;
	unknown x(vx,0);
	unknown y(vy,1);
	unknown z(vz,2);
	vMatrix UVW(3,1);
	UVW(0,0) = 1*x+2*y+3*z;
	UVW(1,0) = 4*x+5*y+6*z;
	UVW(2,0) = 7*x+8*y+9*z;
	return UVW;
}

variable refP(real x, real y, real z, real t, real nu)
{
	(void)x,(void)y,(void)z,(void)t,(void)nu;
	return 0;
}


real refq(real vx, real vy, real vz, real t, real nu, const rMatrix & n)
{
	rMatrix UVW = refUVW(vx,vy,vz,t,nu);
	return n.DotProduct(UVW);
}

rMatrix reft(real vx, real vy, real vz, real t, real nu, real rho, const rMatrix & n)
{
	static const rMatrix I = rMatrix::Unit(3);
	rMatrix N;
	//if( fullstress)
	//	N = 0.5*(I.Kronecker(n.Transpose())+n.Transpose().Kronecker(I));
	//else
		N = I.Kronecker(n.Transpose());
	vMatrix UVW = refUVW(vx,vy,vz,t,nu);
	rMatrix vUVW = UVW;
	rMatrix G(9,1);
	//real p = (addpressure || addpressure_analytic) ? get_value(refP(vx,vy,vz,t,nu)) : 0.0;
	real p = get_value(refP(vx,vy,vz,t,nu));
	for(int k = 0; k < 3; ++k)
		for(int l = 0; l < 3; ++l)
			G(l+3*k,0) = UVW(k,0).GetRow()[l];
	rMatrix ret = rho*UVW*UVW.Transpose()*n - nu*N*G + p*n;
	return ret;
}

int main(int argc, char ** argv)
{
	
	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " mesh time=1 nu=0.1 rho=1 mesh_out=grid_out.pmf" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try
	{
		m->Load(argv[1]);
	}
	catch(...)
	{
		std::cout << "Cannot load the mesh " << argv[1] << std::endl;
		return -1;
	}
	
	std::string fout = "grid_out.pmf";
	double time = 1;
	double nu = 0.1;
	double rho = 1;
	if( argc > 2 ) time = atof(argv[2]);
	if( argc > 3 ) nu = atof(argv[3]);
	if( argc > 4 ) rho = atof(argv[4]);
	if( argc > 5 ) fout = std::string(argv[5]);
	
	TagRealArray force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,3);
	TagRealArray bc = m->CreateTag("BOUNDARY_CONDITION_VELOCITY",DATA_REAL,FACE,FACE,7);
	//TagReal bcp = m->CreateTag("BOUNDARY_CONDITION_PRESSURE",DATA_REAL,FACE,FACE,1);
	TagRealArray uvw = m->CreateTag("REFERENCE_VELOCITY",DATA_REAL,CELL,NONE,3);
	TagReal p = m->CreateTag("REFERENCE_PRESSURE",DATA_REAL,CELL,NONE,1);
	
	rMatrix x(3,1), n(3,1);
	
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		it->Centroid(x.data());
		uvw(*it,3,1) = refUVW(x(0,0),x(1,0),x(2,0),time,nu);
		p[*it] = get_value(refP(x(0,0),x(1,0),x(2,0),time,nu));
		
		ElementArray<Face> faces = it->getFaces();
		force(*it,3,1).Zero();
		for(ElementArray<Face>::iterator f = faces.begin(); f != faces.end(); ++f)
		{
			f->Barycenter(x.data());
			f->OrientedUnitNormal(it->self(),n.data());
			force(*it,3,1) += reft(x(0,0),x(1,0),x(2,0),time,nu,rho,n)*f->Area();
		}
		force(*it,3,1) /= it->Volume();
	}
	
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it)
	{
		it->Centroid(x.data());
		it->UnitNormal(n.data());
		if( it->Boundary() )
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 1;
			bc[*it][3] = 0;
			bc(*it,7,1)(4,7,0,1) = refUVW(x(0,0),x(1,0),x(2,0),time,nu);
		}
	}
	
	std::cout << "Saving output to " << fout << std::endl;
	
	m->Save(fout);
	
	return 0;
}
