#include "inmost.h"

using namespace INMOST;


const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "perturbed parallelogram";

static void GetBox(Element e, Storage::real min[3], Storage::real max[3])
{
	min[0] = min[1] = min[2] = 1.0e20;
	max[0] = max[1] = max[2] = -1.0e20;
	ElementArray<Node> nodes = e->getNodes();
	for (ElementArray<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		Storage::real_array c = it->Coords();
		for (int i = 0; i < (int)c.size(); i++) 
		{
			if (max[i] < c[i]) max[i] = c[i]; //max
			if (min[i] > c[i]) min[i] = c[i]; //min
		}
	}
	for(int i = 0; i < 3; ++i)
	{
		if( max[i] < min[i] )
		{
			max[i] = 0.0001;
			min[i] = -0.0001;
		}
		else if( max[i] == min[i] )
		{
			max[i] += 0.0001;
			min[i] += -0.0001;
		}
	}
}

int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out [incline:1 X:1 Y:1/30 angle:30 perturbation:0.01]" << std::endl;
		return 0;
	}

	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try{m->Load(argv[1]);} catch(...) { std::cout << "Cannot load the mesh " << argv[1] << std::endl; return -1;}

	Storage::real angle = 30, X = 1.0, Y = 1.0/30.0, tanangle, perturbation = 0.01;

	if( argc > 4 ) X = atof(argv[4]);
	if( argc > 5 ) Y = atof(argv[5]);
	if( argc > 6 ) angle = atof(argv[6]);
	if( argc > 7 ) perturbation = atof(argv[7]);

	angle = angle/180.0*pi;
	tanangle = tan(angle);

	//for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); ++t)
	//  if( t->GetTagName().substr(0,9) == "GEOM_UTIL" ) m->DeleteTag(*t);

	if( argc <= 3 || (argc > 3 && atoi(argv[3])) )
	{
		double max[3] = {-1.0e20, -1.0e20, -1.0e20}, min[3] = {1.0e20,1.0e20,1.0e20};
		for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
		{
			Storage::real_array c = it->Coords();
			if( c[0] > max[0] ) max[0] = c[0];
			if( c[1] > max[1] ) max[1] = c[1];
			if( c[2] > max[2] ) max[2] = c[2];
			if( c[0] < min[0] ) min[0] = c[0];
			if( c[1] < min[1] ) min[1] = c[1];
			if( c[2] < min[2] ) min[2] = c[2];
		}

		srand(0);

		std::cout << "You asked to incline the grid for parallelogram test" << std::endl;
		std::cout << "and I assume that the initial grid is rectangular." << std::endl;
		for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
		{

			Storage::real_array c = it->Coords();
			Storage::real x = c[0], dx = 0;
			Storage::real y = c[1], dy = 0;
			//Storage::real d = 0.0;


			if( x > 0.0 && x < 1.0 && y > 0.0 && y < 1.0 )
			{
				dx = 2.0*(static_cast<Storage::real>(rand())/static_cast<Storage::real>(RAND_MAX)-0.5)*perturbation*X;
				dy = 2.0*(static_cast<Storage::real>(rand())/static_cast<Storage::real>(RAND_MAX)-0.5)*perturbation*Y;
			}


			x = (x-min[0])/(max[0]-min[0])*X;
			y = (y-min[1])/(max[1]-min[1])*Y;

			x -= tanangle*y;

			c[0] = x+dx;
			c[1] = y+dy;
		}
	}

	Mesh::GeomParam t;
	t[MEASURE] = CELL | FACE;
	t[ORIENTATION] = FACE;
	t[NORMAL] = FACE;
	//t[BARYCENTER] = CELL;
	t[CENTROID] = CELL | FACE | EDGE | NODE;
	m->AssignGlobalID(CELL|FACE);
	//m->RemoveGeometricData(t);
	m->PrepareGeometricData(t);

	double max[3] = {-1.0e20, -1.0e20, -1.0e20}, min[3] = {1.0e20,1.0e20,1.0e20};
	Storage::real c[3] = {0,0,0};//, nrm[3] = {0,0,0};
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


	Tag force = m->CreateTag("FORCE",DATA_REAL,CELL,CELL,1);
	Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,1);
	Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,3);

	{
		Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}

	for(Mesh::iteratorCell c = m->BeginCell(); c != m->EndCell(); ++c)
		c->Real(tensor) = 1.0;

	Storage::real vol = 0.0;
	Storage::real source_x = 0.5*(X-tanangle*Y), source_y = 0.5*Y;
	ElementArray<Cell> source_cells(m);

	for(Mesh::iteratorCell c = m->BeginCell(); c != m->EndCell(); ++c)
	{
		Storage::real max[3],min[3];
		GetBox(c->self(),min,max);
		if( min[0] < source_x && source_x < max[0] && min[1] < source_y && source_y < max[1] )
		{
			source_cells.push_back(c->self());
			vol += c->Volume();
		}
	}

	if( source_cells.empty() ) printf("No source cells found!\n");

	for(ElementArray<Cell>::iterator c = source_cells.begin(); c != source_cells.end(); ++c)
	{
		c->Real(force) = 1.0/vol;
	}


	for(Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
	{
		if( f->Boundary() )
		{
			f->Centroid(c);
			if( !(c[2] < eps && c[2] > 1-eps) )
			{
				Storage::real_array bc = f->RealArray(bndcond);
				bc[0] = 1.0;
				bc[1] = 0.0;
				bc[2] = 0.0;
			}
		}
	}

	std::cout << "Saving output to " << argv[2] << std::endl;

	m->Save(argv[2]);

	delete m;
}


