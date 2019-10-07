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



int main(int argc, char ** argv)
{
	
	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " mesh is2D=1 mesh_out=grid_out.pmf" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try
	{
		m->Load(argv[1]);
		Mesh::GeomParam table;
		table[ORIENTATION] = FACE;
		m->PrepareGeometricData(table);
	}
	catch(...)
	{
		std::cout << "Cannot load the mesh " << argv[1] << std::endl;
		return -1;
	}
	
	std::string fout = "grid_out.pmf";
	int is2D = 1;
	double skew = 0;
	if( argc > 2 ) is2D = atoi(argv[2]);
	if( argc > 3 ) fout = std::string(argv[3]);
	if( argc > 4 ) skew = atof(argv[4]);
	
	if( skew > 0)
	{
		double cmax[3] = {-1.0e20,-1.0e20,-1.0e20}, cmin[3] = {1.0e20,1.0e20,1.0e20};
		for(Mesh::iteratorNode n = m->BeginNode(); n != m->EndNode(); ++n)
		{
			Storage::real_array c = n->Coords();
			for(int k = 0; k < 3; ++k)
			{
				if( cmax[k] < c[k] ) cmax[k] = c[k];
				if( cmin[k] > c[k] ) cmin[k] = c[k];
			}
		}
		for(Mesh::iteratorNode n = m->BeginNode(); n != m->EndNode(); ++n)
		{
			Storage::real_array c = n->Coords();
			for(int k = 0; k < 3; ++k)
			{
				double t = (c[k]-cmin[k])/(cmax[k]-cmin[k]);
				double a = 0.5 + 0.5*(t-0.5)/sqrt(pow(t-0.5,2)+pow(skew,2))*sqrt(0.25+pow(skew,2))*2;
				c[k] = cmin[k] + (cmax[k]-cmin[k])*a;
			}
		}
	}
	
	//TagRealArray force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,3);
	TagRealArray bc = m->CreateTag("BOUNDARY_CONDITION_VELOCITY",DATA_REAL,FACE,FACE,7);
	//TagReal bcp = m->CreateTag("BOUNDARY_CONDITION_PRESSURE",DATA_REAL,FACE,FACE,1);
	TagRealArray uvw = m->CreateTag("UVW",DATA_REAL,CELL,NONE,3);
	TagReal p = m->CreateTag("P",DATA_REAL,CELL,NONE,1);
	
	rMatrix n(3,1);
	
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		uvw(*it,3,1).Zero();
		p[*it] = 0;
	}
	
	if( is2D ) std::cout << "2D setup" << std::endl;
	else std::cout << "3D setup" << std::endl;
	
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
	{
		it->UnitNormal(n.data());
		//std::cout << "face " << it->LocalID() << std::endl;
		//std::cout << std::scientific;
		//std::cout <<  "n " << (n(2,0)-1) << " " << (fabs(n(2,0)-1) < 1.0e-3) << " " << (n(2,0)+1) << " " << (fabs(n(2,0)+1) < 1.0e-3) << std::endl;
		//n.Transpose().Print();
		if( is2D && (fabs(n(2,0)-1) < 1.0e-3 || fabs(n(2,0)+1) < 1.0e-3) ) //roller bc
		{
			//std::cout << "roller bc" << std::endl;
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 0;
			bc[*it][3] = 1;
			bc[*it][4] = 0;
			bc[*it][5] = 0;
			bc[*it][6] = 0;
		}
		else
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 1;
			bc[*it][3] = 0;
			if( fabs(n(1,0)-1) < 1.0e-3 ) //top
			{
				//std::cout << "dirichlet bc top" << std::endl;
				bc[*it][4] = 1;
				bc[*it][5] = 0;
				bc[*it][6] = 0;
			}
			else
			{
				//std::cout << "dirichlet bc" << std::endl;
				bc[*it][4] = 0;
				bc[*it][5] = 0;
				bc[*it][6] = 0;
			}
		}
	}
	
	std::cout << "Saving output to " << fout << std::endl;
	
	m->Save(fout);
	
	return 0;
}
