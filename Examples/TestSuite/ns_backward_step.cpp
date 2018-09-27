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
		std::cout << "Usage: " << argv[0] << " mesh is2D=1 height=0.6 length=0.2 mesh_out=grid_out.pmf" << std::endl;
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
	int is2D = 1;
	double height = 0.7;
	double length = 0.2;
	if( argc > 2 ) is2D = atoi(argv[2]);
	if( argc > 3 ) height = atof(argv[3]);
	if( argc > 4 ) length = atof(argv[4]);
	if( argc > 5 ) fout = std::string(argv[5]);
	
	//TagRealArray force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,3);
	TagRealArray bc = m->CreateTag("BOUNDARY_CONDITION_VELOCITY",DATA_REAL,FACE,FACE,7);
	TagReal bcp = m->CreateTag("BOUNDARY_CONDITION_PRESSURE",DATA_REAL,FACE,FACE,1);
	TagRealArray uvw = m->CreateTag("UVW",DATA_REAL,CELL,NONE,3);
	TagReal p = m->CreateTag("P",DATA_REAL,CELL,NONE,1);
	
	//cut part of the grid
	
	rMatrix n(3,1), x(3,1);
	
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		it->Centroid(x.data());
		if( x(0,0) < length && x(1,0) < height )
			it->Delete();
	}
	
	
	//this should not be needed?
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
		it->FixNormalOrientation();
	
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		uvw(*it,3,1).Zero();
		p[*it] = 0;
	}
	
	if( is2D ) std::cout << "2D setup" << std::endl;
	else std::cout << "3D setup" << std::endl;
	std::cout << "height=" << height << " length=" << length << std::endl;
	
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
	{
		it->Centroid(x.data());
		it->UnitNormal(n.data());
		if( is2D && (fabs(n(2,0)-1) < 1.0-3 || fabs(n(2,0)+1) < 1.0e-3 ) ) //roller bc
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 0;
			bc[*it][3] = 1;
			bc[*it][4] = 0;
			bc[*it][5] = 0;
			bc[*it][6] = 0;
		}
		else if( fabs(n(0,0)+1) < 1.0e-3 && x(0,0) < 1.0e-3 ) //inlet
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 1;
			bc[*it][3] = 0;
			bc[*it][4] = 1;
			bc[*it][5] = 0;
			bc[*it][6] = 0;
		}
		else if( fabs(n(0,0)-1) < 1.0e-3 ) //outlet
			bcp[*it] = 0;
		else //no-slip
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 1;
			bc[*it][3] = 0;
			bc[*it][4] = 0;
			bc[*it][5] = 0;
			bc[*it][6] = 0;
		}
	}
	
	std::cout << "Saving output to " << fout << std::endl;
	
	m->Save(fout);
	
	return 0;
}
