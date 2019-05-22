#include "inmost.h"

using namespace INMOST;


const Storage::real pi = 3.1415926535897932384626433832795;
//const Storage::real eps = 1.0e-6;
std::string problem_name = "two_wells";
  
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
		std::cout << "Usage: " << argv[0] << " mesh mesh_out [angle:67.5 perturbation:0.01 anisotropy_ratio:1000]" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try{m->Load(argv[1]);} catch(...) { std::cout << "Cannot load the mesh " << argv[1] << std::endl; return -1;}

  Storage::real angle = 67.5, perturbation = 0.0, ratio = 1000;

  if( argc > 3 ) angle = atof(argv[3]);
  if( argc > 4 ) perturbation = atof(argv[4]);
  if( argc > 5 ) ratio = atof(argv[5]);

  angle = angle/180.0*pi;

  
  if( perturbation ) 
  {
    
    std::cout << "You asked to perturb the grid with parameter " << perturbation << std::endl;
    for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
    {
      
      Storage::real_array c = it->Coords();
      Storage::real x = c[0], dx = 0;
      Storage::real y = c[1], dy = 0;

      
      if( x > 0.0 && x < 1.0 && y > 0.0 && y < 1.0 )
      {
        dx = 2.0*(static_cast<Storage::real>(rand())/static_cast<Storage::real>(RAND_MAX)-0.5)*perturbation;
        dy = 2.0*(static_cast<Storage::real>(rand())/static_cast<Storage::real>(RAND_MAX)-0.5)*perturbation;
      }
      
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

	if( m->HaveTag("PERM") ) m->DeleteTag(m->GetTag("PERM"));
  
  
  Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,6);
  Tag bndcond = m->CreateTag("BOUNDARY_CONDITION_CELL",DATA_REAL,CELL,CELL,1);
    
    {
        Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		std::stringstream str;
		str << problem_name << "_ratio_" << ratio;
		problem_name = str.str();
        name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
    }


  Storage::real cphi = cos(angle), sphi = sin(angle);
  Storage::real scphi = sphi*cphi, s2phi = sphi*sphi, c2phi = cphi*cphi; 
  Storage::real alpha = ratio, beta = 1;
  Storage::real source_0_x = 3.5/11.0, source_0_y = 0.5, source_1_x = 7.5/11.0, source_1_y = 0.5;

  int numcells = 0;
  
  for(Mesh::iteratorCell c = m->BeginCell(); c != m->EndCell(); ++c)
  {
    Storage::real_array perm = c->RealArray(tensor);
    perm[0] = alpha * c2phi + beta * s2phi;
    perm[1] = (alpha-beta) * scphi;
    perm[2] = 0.0;
    perm[3] = alpha * s2phi + beta * c2phi;
    perm[4] = 0;
    perm[5] = 1;

    Storage::real max[3], min[3];
    GetBox(c->self(),min,max);

    if( min[0] < source_0_x && source_0_x < max[0] && min[1] < source_0_y && source_0_y < max[1] )
    {
      c->Real(bndcond) = 0.0;
      numcells++;
    }

    if( min[0] < source_1_x && source_1_x < max[0] && min[1] < source_1_y && source_1_y < max[1] )
    {
      c->Real(bndcond) = 1.0;
      numcells++;
    }

  }
  std::cout << "Fixed values for cells " << numcells << std::endl;

  std::cout << "Saving output to " << argv[2] << std::endl;

	m->Save(argv[2]);

	delete m;
}


