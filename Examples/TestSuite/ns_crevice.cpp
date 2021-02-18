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
		std::cout << "Usage: " << argv[0] << " mesh [mesh_out=grid_out.pmf]" << std::endl;
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
	if( argc > 2 ) fout = std::string(argv[2]);
	
	
	
	
	
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
	TagReal       bcp     = m->CreateTag("BOUNDARY_CONDITION_PRESSURE",DATA_REAL,FACE,FACE,1);
	TagInteger    bb      = m->CreateTag("BOUNDARY_BLOOD",DATA_INTEGER,FACE,FACE,1);
	TagRealArray  uvw     = m->CreateTag("UVW",DATA_REAL,CELL,NONE,3);
	TagReal       p       = m->CreateTag("P",DATA_REAL,CELL,NONE,1);
	
	
	
	//this should not be needed?
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
		it->FixNormalOrientation();
	
	{ // prepare geometrical data on the mesh
		Mesh::GeomParam table;
		table[CENTROID]    = CELL | FACE; //Compute averaged center of mass
		table[NORMAL]      = FACE;        //Compute normals
		table[MEASURE]     = CELL | FACE; //Compute volumes and areas
		table[BARYCENTER]  = CELL | FACE; //Compute volumetric center of mass
		m->RemoveGeometricData(table); //Ask to precompute the data
	}
	//~ bool onep = false;
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
	{
		Storage::real  n[3], c[3];
		it->UnitNormal(n);
		it->Centroid(c);
		bb[*it] = 0;
		if( fabs(n[0]-1) < 1.0e-3 && c[0] > 1.8+eps) // outflow
		{
			//~ bcp[*it] = 0;
			//~ if( !onep )
			{
				//~ onep = true;
				bc[*it][0] = 0;
				bc[*it][1] = 1;
				bc[*it][2] = 0;
				bc[*it][3] = 1;
				bc[*it][4] = 0;
				bc[*it][5] = 0;
				bc[*it][6] = 0;
			}
			//~ else
			//~ {
				//~ bc[*it][0] = 1;
				//~ bc[*it][1] = 0;
				//~ bc[*it][2] = 1;
				//~ bc[*it][3] = 0;
				//~ bc[*it][4] = 17.3 * c[1]*(3.0-c[1])/9.0 * 4.0 * (3.0 / 2.0); //parabolic profile with average of 17.3
				//~ bc[*it][5] = 0;
				//~ bc[*it][6] = 0;
			//~ }
		}
		else if(  fabs(n[0]+1) < 1.0e-3 && c[0] < 0.0+eps) //inflow
		{
			bb[*it] = 1;
			//~ bcp[*it] = 0.13;
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 1;
			bc[*it][3] = 0;
			bc[*it][4] = 17.3 * c[1]*(3.0-c[1])/9.0 * 4.0 * (3.0 / 2.0); //parabolic profile with average of 17.3
			bc[*it][5] = 0;
			bc[*it][6] = 0;
		}
		else if( (fabs(n[1]-1) < 1.0e-3 && c[1] > 1.5-eps) || (fabs(n[2]-1) < 1.0e-3 && c[2] > 0.1-eps) || (fabs(n[2]+1) < 1.0e-3 && c[2] < 0.0+eps) ) //slip wall
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 0;
			bc[*it][3] = 1;
			bc[*it][4] = 0;
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
		}
	}
	
	std::cout << "Saving output to " << fout << std::endl;
	
	m->Save(fout);
	
	return 0;
}
