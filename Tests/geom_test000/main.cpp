#include "inmost.h"
#include <cstdio>
#include <cmath>

using namespace INMOST;

typedef Storage::real_array real_array;
typedef Storage::real real;


__INLINE static void vec_cross_product(const Storage::real * vecin1, const Storage::real * vecin2, Storage::real * vecout)
{
	Storage::real temp[3];
	temp[0] = vecin1[1]*vecin2[2] - vecin1[2]*vecin2[1];
	temp[1] = vecin1[2]*vecin2[0] - vecin1[0]*vecin2[2];
	temp[2] = vecin1[0]*vecin2[1] - vecin1[1]*vecin2[0];
	vecout[0] = temp[0];
	vecout[1] = temp[1];
	vecout[2] = temp[2];
}


int main(int argc,char ** argv)
{
	int errors = 0;
#if defined(USE_FP64)
	const Storage::real tol = 1.0e-9;
#else
	const Storage::real tol = 1.0e-5;
#endif
	Mesh::Initialize(&argc,&argv);
	{
		Mesh m;
		m.Load((argc>1)?argv[1]:"d4.pmf");
		rMatrix c(3,1), n(3,1), x(3,1), y(3,1), z(3,1), ct(3,1), nt(3,1);
		rMatrix e1(3,1,0.0), e2(3,1,0.0), e3(3,1,0.0);
		rMatrix x0(3,1),y0(3,1),z0(3,1);
		rMatrix Et(3,3), Ef(3,3), d(3,1), W(9,3);
		rMatrix nsum(3,1), ntsum(3,1);
		e1(0,0) = 1.0;
		e2(1,0) = 1.0;
		e3(2,0) = 1.0;
		Mesh::GeomParam t;
		t[BARYCENTER] = CELL|FACE;
		t[ORIENTATION] = FACE;
		t[NORMAL] = FACE;
		m.RemoveGeometricData(t);
		m.PrepareGeometricData(t);
		for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
			it->FixNormalOrientation();
		int np = 0;
		for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
			if( !it->Planarity() ) np++;
		std::cout << "nonplanar: " << np << std::endl;
		
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		{
			ElementArray<Face> faces = it->getFaces();
			//face center of mass, face normal, face area, cell volume
			x.Zero(), y.Zero(), z.Zero();
			//~ real vol = 0;
			for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
			{
				jt->Barycenter(c.data());
				jt->OrientedNormal(it->self(),n.data());
				x += c(0,0) * n;
				y += c(1,0) * n;
				z += c(2,0) * n;
			}
			
			
			x /= it->Volume();
			y /= it->Volume();
			z /= it->Volume();

			if( (x-e1).FrobeniusNorm() > tol || x.CheckNans() )
			{
				std::cout << "On CELL:" << it->LocalID() << " dx: " << x(0,0) << " " << x(1,0) << " " << x(2,0) << std::endl;
				x.Print();
				e1.Print();
				std::cout << (x-e1).FrobeniusNorm() << std::endl;
				errors++;
			}
			if( (y-e2).FrobeniusNorm() > tol || y.CheckNans() )
			{
				std::cout << "On CELL:" << it->LocalID() << " dy: " << y(0,0) << " " << y(1,0) << " " << y(2,0) << std::endl;
				errors++;
			}
			if( (z-e3).FrobeniusNorm() > tol || z.CheckNans() )
			{
				std::cout << "On CELL:" << it->LocalID() << " dz: " << z(0,0) << " " << z(1,0) << " " << z(2,0) << std::endl;
				errors++;
			}
		}
	}
	Mesh::Finalize();
	if( errors )
	{
		std::cout << "Errors: " << errors << std::endl;
		return 1;
	}
	return 0;
}
