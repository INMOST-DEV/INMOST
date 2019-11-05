#include "inmost.h"

using namespace INMOST;

// use ADMFD or ADVDIFF to solve for phi

//setup BC
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
// BOUNDARY_CONDTION - 3 entries (for calculation of phi and measures of drag and lift)


typedef Storage::real real;
const double eps = 1.0e-5;


int main(int argc, char ** argv)
{
	
	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " mesh [mesh_out=grid_out.pmf] [Umax=2.25]" << std::endl;
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
	double Umax = 2.25;
	if( argc > 2 ) fout = std::string(argv[2]);
	if( argc > 3 ) Umax = atof(argv[3]);
	
	
	
	
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
	for(int d = 0; d < 3; ++d)
		std::cout << d << " " << cmin[d] << ":" << cmax[d] << std::endl;
			
	
	
	//TagRealArray force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,3);
	TagRealArray  bc      = m->CreateTag("BOUNDARY_CONDITION_VELOCITY",DATA_REAL,FACE,FACE,7);
	TagRealArray  bcphi   = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,3);
	TagReal       bcp     = m->CreateTag("BOUNDARY_CONDITION_PRESSURE",DATA_REAL,FACE,FACE,1);
	TagRealArray  uvw     = m->CreateTag("UVW",DATA_REAL,CELL,NONE,3);
	TagReal       p       = m->CreateTag("P",DATA_REAL,CELL,NONE,1);
	m->self().Real(m->CreateTag("Umax",DATA_REAL,MESH,NONE,1)) = Umax;
	
	
	//this should not be needed?
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
		it->FixNormalOrientation();
	
	
	
	
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
	{
		double  n[3], c[3];
		it->UnitNormal(n);
		it->Barycenter(c);
		bcphi[*it][0] = 1;
		bcphi[*it][1] = 0;
		bcphi[*it][2] = 0;
		if( fabs(n[0]-1) < 1.0e-3 && c[0] > 2.5-eps) // outflow
		{
			bcp[*it] = 0;
		}
		else if(  fabs(n[0]+1) < 1.0e-3 && c[0] < 0.0+eps) //inflow
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 1;
			bc[*it][3] = 0;
			bc[*it][4] = 16*Umax*c[1]*c[2]*(0.41-c[1])*(0.41-c[2])/pow(0.41,4);
			bc[*it][5] = 0;
			bc[*it][6] = 0;
		}
		else //no-slip walls
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 1;
			bc[*it][3] = 0;
			bc[*it][4] = 0;
			bc[*it][5] = 0;
			bc[*it][6] = 0;
			if( fabs(n[2]) < 0.7 && c[0] > 0.45-eps && c[0] < 0.55+eps && c[1] > 0.15-eps && c[1] < 0.25+eps )
				bcphi[*it][2] = 1;
		}
	}
	
	std::cout << "Saving output to " << fout << std::endl;
	
	m->Save(fout);
	
	return 0;
}
