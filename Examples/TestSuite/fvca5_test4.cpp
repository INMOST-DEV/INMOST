#include "inmost.h"

using namespace INMOST;



const Storage::real eps = 1.0e-6;
  std::string problem_name = "vertical_fault";

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

  //Tag force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,1);
  Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,3);
  //Tag solution_val = m->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,1);
  //Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,FACE,1);
  Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,3);
    
    {
        Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
        name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
    }

  for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
    it->Centroid(c);
    Storage::real x = c[0];//(c[0]-min[0])/(max[0]-min[0]);
		Storage::real y = c[1];//(c[1]-min[1])/(max[1]-min[1]);
    //Storage::real z = c[2];//(c[2]-min[2])/(max[2]-min[2]);

    int zone = 2;

    if( x <= 0.5 && static_cast<int>(floor((y-0.05)/0.1))%2 == 0 )
      zone = 1;

    if( x >= 0.5 && static_cast<int>(floor(y/0.1))%2 == 0 )
      zone = 1;
    

    if( it->GetElementType() == CELL )
    {
      Storage::real_array perm = it->RealArray(tensor);
      if( zone == 1)
      {
        perm[0] = 100;
        perm[1] = 10;
      }
      else if( zone == 2 )
      {
        perm[0] = 0.01;
        perm[1] = 0.001;
      }
    }
  }

  for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it)
  {
    if( it->getAsFace()->Boundary() )
    {

      it->Centroid(c);
      Storage::real x = c[0];//(c[0]-min[0])/(max[0]-min[0]);
		//Storage::real y = c[1];//(c[1]-min[1])/(max[1]-min[1]);
      Storage::real z = c[2];//(c[2]-min[2])/(max[2]-min[2]);
    
      if( !(z < eps || z > 1-eps ) )
      {

        Storage::real_array bc = it->RealArray(bndcond);
        bc[0] = 1.0;
        bc[1] = 0.0;
        bc[2] = 1-x;
      }
    }

  }

  std::cout << "Saving output to " << argv[2] << std::endl;

	m->Save(argv[2]);

	delete m;
}


