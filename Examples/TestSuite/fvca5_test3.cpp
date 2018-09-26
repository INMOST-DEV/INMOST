#include "inmost.h"

using namespace INMOST;


const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "oblique_flow";


int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out [sigma:0.0001]" << std::endl;
		return 0;
	}

	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try{m->Load(argv[1]);} catch(...) { std::cout << "Cannot load the mesh " << argv[1] << std::endl; return -1;}

	Mesh::GeomParam t;
	t[MEASURE] = CELL | FACE;
	t[ORIENTATION] = FACE;
	t[NORMAL] = FACE;
	//t[BARYCENTER] = CELL;
	t[CENTROID] = CELL | FACE;
	m->AssignGlobalID(CELL|FACE);
	m->PrepareGeometricData(t);

	double max[3] = {-1.0e20, -1.0e20, -1.0e20}, min[3] = {1.0e20,1.0e20,1.0e20};
	Storage::real c[3] = {0,0,0};
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		it->Centroid(c);
		if( c[0] > max[0] ) max[0] = c[0];
		if( c[1] > max[1] ) max[1] = c[1];
		if( c[2] > max[2] ) max[2] = c[2];
		if( c[0] < min[0] ) min[0] = c[0];
		if( c[1] < min[1] ) min[1] = c[1];
		if( c[2] < min[2] ) min[2] = c[2];
	}

	if( max[0] <= min[0] ) {std::cout << "strange X " << min[0] << ":" << max[0] << std::endl; return -1;}
	if( max[1] <= min[1] ) {std::cout << "strange Y " << min[1] << ":" << max[1] << std::endl; return -1;}
	if( max[2] <= min[2] ) 
	{
		//2d mesh
		if( m->GetDimensions() == 3 )
		{
			//offset from z-plane
			min[2] -= 0.0001;
			max[2] += 0.0001; 
		}
		else
		{
			min[2] = -0.0001;
			max[2] = +0.0001;
		}
	}

	std::cout << "Mesh bounds: " << min[0] << ":" << max[0] << " " << min[1] << ":" << max[1] << " " << min[2] << ":" << max[2] << std::endl;

	Tag material;
	if( m->HaveTag("PERM") ) m->DeleteTag(m->GetTag("PERM"));

	Storage::real sigma = 0.0001;

	if( argc > 3 ) sigma = atof(argv[3]);

	//Tag force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,1);
	Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,6);
	//Tag solution_val = m->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,1);
	//Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,FACE,1);
	Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,3);

	{
		Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}

	Storage::real phi = 40.0/180.0*pi;
	Storage::real cphi = cos(phi), sphi = sin(phi);
	Storage::real scphi = sphi*cphi, s2phi = sphi*sphi, c2phi = cphi*cphi; 

	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		it->Centroid(c);


		if( it->GetElementType() == CELL )
		{
			Storage::real_array perm = it->RealArray(tensor);
			perm[0] = 1.0 * c2phi + sigma * s2phi;
			perm[1] = (1.0-sigma) * scphi;
			perm[2] = 0.0;
			perm[3] = 1.0 * s2phi + sigma * c2phi;
			perm[4] = 0.0;
			perm[5] = 1.0;
		}
	}

	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it)
	{
		if( it->getAsFace()->Boundary() )
		{

			it->Centroid(c);
			Storage::real x = c[0];//(c[0]-min[0])/(max[0]-min[0]);
			Storage::real y = c[1];//(c[1]-min[1])/(max[1]-min[1]);
			Storage::real z = c[2];//(c[2]-min[2])/(max[2]-min[2]);

			if( !(z < eps || z > 1-eps ) )
			{

				Storage::real_array bc = it->RealArray(bndcond);
				bc[0] = 1.0;
				bc[1] = 0.0;

				if( (y < eps && (x >= 0.0 && x <= 0.2)) || (x < eps && (y >= 0.0 && y <= 0.2)) )
					bc[2] = 1;
				else if( (y > 1-eps && (x >= 0.8 && x <= 1.0)) || (x > 1 - eps && (y >= 0.8 && y <= 1.0)) )
					bc[2] = 0;
				else if( (y < eps && (x >= 0.3 && x <= 1.0)) || (x < eps && (y >= 0.3 && y <= 1.0)) )
					bc[2] = 0.5;
				else if( (y > 1-eps && (x >= 0.0 && x <= 0.7)) || (x > 1 - eps && (y >= 0.0 && y <= 0.7)) )
					bc[2] = 0.5;
				else if( y > 1-eps && (x >= 0.7 && x <= 0.8 ) )
					bc[2] = 0.5 - 5.0*(x-0.7);
				else if( x > 1-eps && (y >= 0.7 && y <= 0.8 ) )
					bc[2] = 0.5 - 5.0*(y-0.7);
				else if( (y < eps && (x >= 0.2 && x <= 0.3)) )
					bc[2] = 1.0 - 10.0*(x-0.2);
				else if( (x < eps && (y >= 0.2 && y <= 0.3)) )
					bc[2] = 1.0 - 10.0*(y-0.2);
				else printf("%s:%d ooops boundary at (%g,%g)\n",__FILE__,__LINE__,x,y);
			}
		}

	}

	std::cout << "Saving output to " << argv[2] << std::endl;

	m->Save(argv[2]);

	delete m;
}


