#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
typedef Storage::real real;
typedef Storage::real_array real_array;

int main(int argc, char ** argv)
{
	if( argc < 2 )
	{
		printf("Usage: %s input_mesh [output_mesh]\n",argv[0]);
		return -1;
	}
	
	Mesh A("");
	A.Load(argv[1]);
	
	Mesh::GeomParam table;
	table[ORIENTATION] = FACE;
	table[BARYCENTER] = FACE | CELL | EDGE;
	A.PrepareGeometricData(table);
	
	rMatrix n(3,1), x(3,1), E(3,3), I(3,3), cx(3,1);
	I = rMatrix::Unit(3);
	
	double vol;
	int collapsed = 0;
	//MarkerType rev = A.CreatePrivateMarker();
	MarkerType set = A.CreateMarker();
	for(Mesh::iteratorCell it = A.BeginCell(); it != A.EndCell(); ++it)
	{
		ElementArray<Face> faces = it->getFaces();
		// sum(x*n^T - (n^T*x)/3.0*I) = 0
		double vol = 0,vol0 = 0;
		E.Zero();
		it->Centroid(cx.data());
		//A.FacesOrientation(faces,rev);
		
		for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
		{
			jt->Barycenter(x.data());
			//jt->Normal(n.data());
			//n*= jt->GetPrivateMarker(rev) ? -1.0:1.0;
			jt->OrientedNormal(it->self(),n.data());
			E += (x-cx)*n.Transpose();
			vol += n.DotProduct(x-cx)/3.0;
			vol0 += n.DotProduct(cx)/3.0;
		}
		
		//faces.RemPrivateMarker(rev);
		E = E/vol-I;
		if( E.FrobeniusNorm() > 1.0e-5 || vol0 > 1.0e-5)
		{
			it->SetMarker(set);
			collapsed ++;
		}
	}
	//A.ReleasePrivateMarker(rev);
	for(Mesh::iteratorCell it = A.BeginCell(); it != A.EndCell(); ++it) if( it->GetMarker(set) )
		it->Delete();
	A.ReleaseMarker(set);
	std::cout << "total collapsed: " << collapsed << std::endl;
	
	for(ElementType et = FACE; et > NONE; et = PrevElementType(et) )
		for(Mesh::iteratorElement it = A.BeginElement(et); it != A.EndElement(); ++it)
			if( it->nbAdjElements(CELL) == 0 ) it->Delete();
	
	
	if( A.HaveTag("GRIDNAME") )
	{
		Storage::bulk_array nameA = A.self().BulkArray(A.GetTag("GRIDNAME"));
		std::string ins = "_collapse_degenerate";
		nameA.insert(nameA.end(),ins.c_str(),ins.c_str()+ins.size());
	}
	
	if( argc > 2 )
	{
		std::cout << "Save to " << argv[2] << std::endl;
		A.Save(argv[2]);
	}
	else
	{
		std::cout << "Save to out.vtk" << std::endl;
		A.Save("out.vtk");
	}

	return 0;
}
