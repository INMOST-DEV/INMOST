#include "inmost.h"

using namespace INMOST;


//const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "mild_anisotropy_2";

int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try{m->Load(argv[1]);} catch(...) { std::cout << "Cannot load the mesh " << argv[1] << std::endl; return -1;}

	

	double max[3] = {-1.0e20, -1.0e20, -1.0e20}, min[3] = {1.0e20,1.0e20,1.0e20};
  Storage::real c[3] = {0,0,0}, nrm[3] = {0,0,0};
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
	
	std::cout << "Original mesh bounds: " << min[0] << ":" << max[0] << " " << min[1] << ":" << max[1] << " " << min[2] << ":" << max[2] << std::endl;

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
	
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		for(int k = 0; k < 3; ++k)
			it->Coords()[k] = (it->Coords()[k]-min[k])/(max[k]-min[k]);
	}
	
	max[0] = max[1] = max[2] = -1.0e20;
	min[0] = min[1] = min[2] = 1.0e20;
	
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
	
	std::cout << "New mesh bounds: " << min[0] << ":" << max[0] << " " << min[1] << ":" << max[1] << " " << min[2] << ":" << max[2] << std::endl;
	
	Mesh::GeomParam t;
	t[MEASURE] = CELL | FACE;
	t[ORIENTATION] = FACE;
	t[NORMAL] = FACE;
	//t[BARYCENTER] = CELL;
	t[CENTROID] = CELL | FACE;
	m->AssignGlobalID(CELL|FACE);
	m->PrepareGeometricData(t);
	
	

	Tag material;
	if( m->HaveTag("PERM") ) m->DeleteTag(m->GetTag("PERM"));

	

  Tag force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,1);
  Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,6);
  Tag solution_val = m->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,1);
  Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,FACE,1);
  Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE|EDGE,FACE|EDGE,3);

    {
        Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
        name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
    }
  
	for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE|EDGE|NODE); it != m->EndElement(); ++it)
	{
    it->Centroid(c);
    
    if( it->GetElementType() == CELL )
    {
      Storage::real_array perm = it->RealArray(tensor);
      perm[0] = 1.5;
      perm[1] = 0.5;
      perm[2] = 0.0;
      perm[3] = 1.5;
      perm[4] = 0.0;
      perm[5] = 1.5;
    }

		

    
    Storage::real alpha = c[0];//(c[0]-min[0])/(max[0]-min[0]);
		Storage::real beta  = c[1];//(c[1]-min[1])/(max[1]-min[1]);
    Storage::real gamma = c[2];//(c[2]-min[2])/(max[2]-min[2]);
    Storage::real sol = 0.5*(sin((1.0-alpha)*(1.0-beta))/sin(1.0) + pow(1-alpha,3)*pow(1-beta,2));
    Storage::real flux;

    Storage::real dsolx = 0.5*(-(1-beta)*cos((1.0-alpha)*(1.0-beta))/sin(1.0) - 3*pow(1-alpha,2)*pow(1-beta,2));
    Storage::real dsoly = 0.5*(-(1-alpha)*cos((1.0-alpha)*(1.0-beta))/sin(1.0) - 2*pow(1-alpha,3)*pow(1-beta,1));

    if( it->GetElementType() == CELL )
    it->Real(force) = 
      0.5*((3*pow(beta,2)+(2*alpha-8)*beta+3*pow(alpha,2)-8*alpha+8)*sin((alpha-1)*beta-alpha+1)-2*cos((alpha-1)*beta-alpha+1)+(18*sin(1.0)*alpha-18*sin(1.0))*pow(beta,2)+(12*sin(1.0)*pow(alpha,2)-60*sin(1.0)*alpha+48*sin(1.0))*beta+6*sin(1.0)*pow(alpha,3)-30*sin(1.0)*pow(alpha,2)+60*sin(1.0)*alpha-36*sin(1.0))/(2*sin(1.0));

    it->Real(solution_val) = sol;
    
    if( it->GetElementType() == FACE )
    {
      it->getAsFace()->UnitNormal(nrm);
      flux = (1.5*nrm[0]+0.5*nrm[1])*dsolx + (1.5*nrm[1]+0.5*nrm[0])*dsoly;
      
      if( !(gamma < eps || gamma > 1.0-eps ) )
      {
        it->Real(solution_flux) = -flux;

        if( it->getAsFace()->Boundary() )
        {
          Storage::real_array bc = it->RealArray(bndcond);
          bc[0] = 1.0;
          bc[1] = 0.0;
          bc[2] = sol;
        }
      }
    }
    if( it->GetElementType() == EDGE )
    {
      if( it->Boundary() )
      {
        Storage::real_array bc = it->RealArray(bndcond);
        bc[0] = 1.0;
        bc[1] = 0.0;
        bc[2] = sol;
      }
    }
  }

  std::cout << "Saving output to " << argv[2] << std::endl;

	m->Save(argv[2]);

	delete m;
}


