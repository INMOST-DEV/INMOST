#include "inmost.h"

using namespace INMOST;


const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "wugao_test2";
  

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

	

  Tag force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,1);
  Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,6);
  Tag solution_val = m->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,1);
  Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,FACE,1);
  Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,3);
    
    {
        Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
        name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
    }


  
	for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE|EDGE|NODE); it != m->EndElement(); ++it)
	{
    it->Centroid(c);
    Storage::real x = c[0];//(c[0]-min[0])/(max[0]-min[0]);
		Storage::real y = c[1];//(c[1]-min[1])/(max[1]-min[1]);
    Storage::real z = c[2];//(c[2]-min[2])/(max[2]-min[2]);
    

    if( it->GetElementType() == CELL )
    {
      Storage::real_array perm = it->RealArray(tensor);

      if( x <= 0.5 )
      {
        perm[0] = 10;
        perm[1] = 2;
        perm[2] = 0.0;
        perm[3] = 5;
        perm[4] = 0.0;
        perm[5] = 1;
      }
      else
      {
        perm[0] = 1;
        perm[1] = 0.0;
        perm[2] = 0.0;
        perm[3] = 1;
        perm[4] = 0.0;
        perm[5] = 1;
      }
    }

		
    

    
    Storage::real sol, dsolx, dsoly;
    Storage::real flux;

    if( x <= 0.5 )
    {
      sol = (1+(x-0.5)*(0.1+8*pi*(y-0.5)))*exp(-20*pi*(y-0.5)*(y-0.5));
      dsolx = ((40*pi*(2*y-1)+1)*exp(-20*pi*y*y+20*pi*y-5*pi))/10;
      dsoly = -pi*(2*y*(2*(80*pi*x*y-40*pi*y-80*pi*x+x)+80*pi+19)+10*(8*pi-1)*x-40*pi-15)*exp(-20*pi*y*y+20*pi*y-5*pi);
      if( it->GetElementType() == CELL )
        it->Real(force) = 
          -2*pi*(8*y*(5*pi*(y*(10*(80*pi*x*y-40*pi*y-120*pi*x+x)+600*pi+79)+10*(60*pi-7)*x-300*pi-49)-1)-10*(400*pi*pi-130*pi+1)*x+5*(400*pi*pi+38*pi-15))*exp(-20*pi*y*y+20*pi*y-5*pi);
      
    }
    else
    {
      sol = exp(x-0.5)*exp(-20*pi*(y-0.5)*(y-0.5));
      dsolx = exp(-20*pi*y*y+20*pi*y+x-5*pi-0.5);
      dsoly = -20*pi*(2*y-1)*exp(-20*pi*y*y+20*pi*y+x-5*pi-0.5);
      if( it->GetElementType() == CELL )
        it->Real(force) = -(40*pi*(10*pi*(2*y-1)*(2*y-1)-1)+1)*exp(-20*pi*y*y+20*pi*y+x-5*pi-0.5);

      
    }

    
    
    it->Real(solution_val) = sol;
    
    if( it->GetElementType() == FACE )
    {
      it->getAsFace()->UnitNormal(nrm);
      if( x <= 0.5 )
        flux = (10*nrm[0]+2*nrm[1])*dsolx + (2*nrm[0]+5*nrm[1])*dsoly;
      else
        flux = (1*nrm[0]+0*nrm[1])*dsolx + (0*nrm[0]+1*nrm[1])*dsoly;
      
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
    }
  }

  std::cout << "Saving output to " << argv[2] << std::endl;

	m->Save(argv[2]);

	delete m;
}


