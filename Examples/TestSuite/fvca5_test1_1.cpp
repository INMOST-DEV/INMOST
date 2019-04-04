#include "inmost.h"

using namespace INMOST;


const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "mild_anisotropy_1";
  

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

  Mesh::GeomParam t;
	t[MEASURE] = CELL | FACE;
	t[ORIENTATION] = FACE;
	t[NORMAL] = FACE;
	//t[BARYCENTER] = CELL;
	t[CENTROID] = CELL | FACE;
	m->AssignGlobalID(CELL|FACE);
	m->PrepareGeometricData(t);

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

	

	Tag vel = m->CreateTag("REFERENCE_VELOCITY",DATA_REAL,CELL,NONE,3);
  Tag force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,1);
  Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,6);
  Tag solution_val = m->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,1);
  Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,FACE,1);
	Tag saddle_val = m->CreateTag("SADDLE_SOLUTION",DATA_REAL,CELL,NONE,4);
  Tag saddle_flux = m->CreateTag("SADDLE_FLUX",DATA_REAL,FACE,NONE,4);
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
      perm[5] = 1.0;
    }

		
    

    
    Storage::real x = c[0];//(c[0]-min[0])/(max[0]-min[0]);
		Storage::real y = c[1];//(c[1]-min[1])/(max[1]-min[1]);
    Storage::real z = c[2];//(c[2]-min[2])/(max[2]-min[2]);
    Storage::real sol = 16.0*(1-x)*y*(1-y);
    Storage::real flux;

    Storage::real dsolx = 16*(y-1)*y;
    Storage::real dsoly = 16*(x-1)*(2*y-1);

    if( it->GetElementType() == CELL )
      it->Real(force) = -16*(2*y+3*x-4);

    it->Real(solution_val) = sol;
		
		if( it->GetElementType() == CELL )
		{
		it->RealArray(vel)[0] = 1.5*dsolx + 0.5*dsoly;
		it->RealArray(vel)[1] = 0.5*dsolx + 1.5*dsoly;
		it->RealArray(vel)[2] = 0;
			
			it->RealArray(saddle_val)[0] = 1.5*dsolx + 0.5*dsoly;
			it->RealArray(saddle_val)[1] = 0.5*dsolx + 1.5*dsoly;
			it->RealArray(saddle_val)[2] = 0;
			it->RealArray(saddle_val)[3] = sol;
		}
    
    if( it->GetElementType() == FACE )
    {
      it->getAsFace()->UnitNormal(nrm);
      flux = (1.5*nrm[0]+0.5*nrm[1])*dsolx + (0.5*nrm[0]+1.5*nrm[1])*dsoly;
      
      if( !(z < eps || z > 1.0-eps ) )
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
		
		it->RealArray(saddle_flux)[0] = sol*(1.5*nrm[0]+0.5*nrm[1]);
		it->RealArray(saddle_flux)[1] = sol*(0.5*nrm[0]+1.5*nrm[1]);
		it->RealArray(saddle_flux)[2] = sol*nrm[2];
		it->RealArray(saddle_flux)[3] = flux;
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


