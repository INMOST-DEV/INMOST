#include "inmost.h"

using namespace INMOST;


//#define MAKE_DIRICHLET_BC_FOR_NEUMANN
const Storage::real sol_shift = 1;
const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;

std::string problem_name = "parametric_locking";
  

int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out [parameter:0.01 boundary:A flip:0 alpha:0]" << std::endl;
    std::cout << "bounday: A - dirichlet, B - half dirichlet, half neumann, C - almost neumann except 2 faces" << std::endl;
    std::cout << "flip: 0 - as in paper by Manzini, 1 - as in FVCA test case." << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try{m->Load(argv[1]);} catch(...) { std::cout << "Cannot load the mesh " << argv[1] << std::endl; return -1;}
  std::cout << std::endl;

  double alpha = 0.0;
  if( argc > 6 ) alpha = atof(argv[6]);
  if( alpha )
  {
    std::cout << "Disturbing mesh with the parameter " << alpha << std::endl;
    Storage::real h = 0;
    for(Mesh::iteratorCell e = m->BeginCell(); e != m->EndCell(); ++e)
    {
		Storage::real maxmin[6];
		maxmin[0] = -1e20;
	    maxmin[1] = 1e20;
	    maxmin[2] = -1e20;
	    maxmin[3] = 1e20;
	    maxmin[4] = -1e20;
	    maxmin[5] = 1e20;
	    ElementArray<Node> nodes = e->getNodes();
	    for (ElementArray<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
	    {
		    Storage::real_array c = it->Coords();
		    for (int i = 0; i < (int)c.size(); i++) {
			    if (maxmin[2 * i + 0] < c[i]) maxmin[2 * i + 0] = c[i]; //max
			    if (maxmin[2 * i + 1] > c[i]) maxmin[2 * i + 1] = c[i]; //min
		    }
		    if( c.size() < 3 )
		    {
			    for(int i = c.size(); i < 3; i++)
			    {
				    maxmin[2*i+0] = maxmin[2*i+1] = 0;
			    }
		    }
	    }
      h = 0;
      h = std::max(h,maxmin[0]-maxmin[1]);
      h = std::max(h,maxmin[2]-maxmin[3]);
      //h = std::max(h,maxmin[4]-maxmin[5]);
    }
    srand(0);
    const Storage::real eps = 1.0e-6;
    Storage::real cnt[3], dx, dy;
    for( Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it)
    {
      it->Centroid(cnt);
      if( cnt[0] > eps && cnt[0] < 1.0-eps && cnt[1] > eps && cnt[1] < 1.0-eps && cnt[2] > eps && cnt[2] < 1.0-eps )
      {
        dx = alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);
        dy = alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);
        it->getBeg()->Coords()[0] += dx;
        it->getBeg()->Coords()[1] += dy;
        it->getEnd()->Coords()[0] += dx;
        it->getEnd()->Coords()[1] += dy;
      }
    }
  }


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
	if( m->HaveTag("MATERIAL") ) material = m->GetTag("MATERIAL");
	if( m->HaveTag("PERM") ) m->DeleteTag(m->GetTag("PERM"));

	Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,3);

  
  Storage::real parameter = 0.01;
  Storage::real vol = 0;
  char bc_type = 'a';
  int flip = 0;

  if( argc > 3 ) parameter = atof(argv[3]);
  if( argc > 4 ) bc_type = tolower(argv[4][0]);
  if( argc > 5 ) flip = atoi(argv[5]);
  
  std::cout << "Setting parametric locking permeability with " << parameter << std::endl;

	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		it->Centroid(c);
    vol += it->Volume();
    Storage::real_array perm = it->RealArray(tensor);
    perm[0] = 1.0;
    perm[1] = parameter;
    perm[2] = 1.0;
	}
    
    {
        Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
        std::stringstream str;
        str << problem_name << "_param_" << parameter << "_flip_" << flip;
        problem_name = str.str();
        name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
    }

  m->self()->Real(m->CreateTag("INTEGRAL_CONSTRAINT",DATA_REAL,MESH,NONE,1)) = vol*sol_shift;

  std::cout << "Setting custom boundary conditions " << (char)toupper(bc_type) << std::endl;

  std::cout << "Problem flip: " << flip << std::endl;


  Storage::real mesh_radius = 0.0;

  Tag solution_val = m->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL | FACE | EDGE | NODE,NONE,1);
  Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,FACE,1);
	Tag rvel = m->CreateTag("REFERENCE_VELOCITY",DATA_REAL,CELL,NONE,3);

  for(Mesh::iteratorCell e = m->BeginCell(); e != m->EndCell(); ++e)
  {
    Storage::real maxmin[6];
    maxmin[0] = -1e20;
	  maxmin[1] = 1e20;
	  maxmin[2] = -1e20;
	  maxmin[3] = 1e20;
	  maxmin[4] = -1e20;
	  maxmin[5] = 1e20;
	  ElementArray<Node> nodes = e->getNodes();
	  for (ElementArray<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
	  {
		  Storage::real_array c = it->Coords();
		  for (int i = 0; i < (int)c.size(); i++) {
			  if (maxmin[2 * i + 0] < c[i]) maxmin[2 * i + 0] = c[i]; //max
			  if (maxmin[2 * i + 1] > c[i]) maxmin[2 * i + 1] = c[i]; //min
		  }
		  if( c.size() < 3 )
		  {
			  for(int i = c.size(); i < 3; i++)
			  {
				  maxmin[2*i+0] = maxmin[2*i+1] = 0;
			  }
		  }
	  }
    mesh_radius = std::max(mesh_radius,maxmin[0]-maxmin[1]);
    mesh_radius = std::max(mesh_radius,maxmin[2]-maxmin[3]);
    //mesh_radius = std::max(mesh_radius,maxmin[4]-maxmin[5]);
  }
  mesh_radius *= 1;

  std::cout << "XY mesh element radius " << mesh_radius << std::endl;

  Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE|EDGE|NODE,FACE|EDGE|NODE,3);

  

  int total_faces = 0, boundary_faces = 0;
  int dirichlet_bc = 0, neumann_bc = 0;

  for(Mesh::iteratorCell e = m->BeginCell(); e != m->EndCell(); ++e)
  {
    e->Centroid(c);
    //project to [0,1]
	  double G, dGdx, dGdy;
    Storage::real alpha = c[0];//(c[0]-min[0])/(max[0]-min[0]);
    Storage::real beta  = c[1];//(c[1]-min[1])/(max[1]-min[1]);
    //solution is exp(-2 pi sqrt(parameter) x) sin(2 pi y)
	  
	  if( flip )
	  {
		  G = exp(-2*pi*sqrt(1.0/parameter)*beta)*sin(2*pi*alpha)+sol_shift; //solution
		  dGdx = 2*pi*cos(2*pi*alpha)*exp(-2*pi*sqrt(1.0/parameter)*beta);
		  dGdy = -2*pi*sqrt(1.0/parameter)*sin(2*pi*alpha)*exp(-2*pi*sqrt(1.0/parameter)*beta);
	  }
	  else
	  {
		  G = exp(-2*pi*sqrt(parameter)*alpha)*sin(2*pi*beta)+sol_shift; //solution
		  dGdx = -2*pi*sqrt(parameter)*exp(-2*pi*sqrt(parameter)*alpha)*sin(2*pi*beta); //x derivative
		  dGdy = 2*pi*exp(-2*pi*sqrt(parameter)*alpha)*cos(2*pi*beta);
	  }
	  
      e->Real(solution_val) = G;
	  
	  Storage::real_array perm = e->RealArray(tensor);
	  e->RealArray(rvel)[0] = perm[0]*dGdx;
	  e->RealArray(rvel)[1] = perm[1]*dGdy;
	  e->RealArray(rvel)[2] = perm[2]*0;

  }
  
  for(Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
  {
    f->Centroid(c);
    //project to [0,1]
    Storage::real alpha = c[0];//(c[0]-min[0])/(max[0]-min[0]);
		Storage::real beta  = c[1];//(c[1]-min[1])/(max[1]-min[1]);
    Storage::real gamma = c[2];//(c[2]-min[2])/(max[2]-min[2]);
    //normal direction is defined by construction, normal in calculation should match
    Storage::real nrm[3];
    f->OrientedUnitNormal(f->BackCell(),nrm);
    Storage::real G, dGdx, dGdy;
    if( flip )
    {
      G = exp(-2*pi*sqrt(1.0/parameter)*beta)*sin(2*pi*alpha)+sol_shift; //solution
      dGdx = 2*pi*cos(2*pi*alpha)*exp(-2*pi*sqrt(1.0/parameter)*beta);
      dGdy = -2*pi*sqrt(1.0/parameter)*sin(2*pi*alpha)*exp(-2*pi*sqrt(1.0/parameter)*beta);
    }
    else
    {
      G = exp(-2*pi*sqrt(parameter)*alpha)*sin(2*pi*beta)+sol_shift; //solution
      dGdx = -2*pi*sqrt(parameter)*exp(-2*pi*sqrt(parameter)*alpha)*sin(2*pi*beta); //x derivative
      dGdy = 2*pi*exp(-2*pi*sqrt(parameter)*alpha)*cos(2*pi*beta);
    }

    Storage::real Flux = (dGdx*nrm[0] + dGdy*parameter*nrm[1]);

    if( !(gamma < eps || gamma > 1.0-eps ) ) //do not write reference flux to top or bottom
    {
      f->Real(solution_flux) = Flux;
    }
    f->Real(solution_val) = G;

    total_faces++;
      
    if( f->Boundary() 
      && (alpha < eps || alpha > 1.0-eps || beta < eps || beta > 1.0 - eps)  ) //restrict top and bottom to be pure neumann bc
    {
      // bnd[0]*p + bnd[1]*K dp/dn = bnd[2]

      Storage::real gD = G; //solution is exp(-2 pi sqrt(parameter) x) sin(2 pi y)
      Storage::real gN = Flux; //boundary normals should always direct outside

      Storage::real_array bnd = f->RealArray(bndcond);
      if( bc_type == 'a' )
      {
        bnd[0] = 1.0;
        bnd[1] = 0.0;
        bnd[2] = gD;
        dirichlet_bc++;
      }
      else if( bc_type == 'b' )
      {
        if( alpha < eps || beta < eps )
        {
          bnd[0] = 1.0;
          bnd[1] = 0.0;
          bnd[2] = gD;
          dirichlet_bc++;
        }
        else
        {
          bnd[0] = 0.0;
          bnd[1] = 1.0;
          bnd[2] = gN;
          neumann_bc++;
        }
      }
      else if( bc_type == 'c' )
      {
#if defined(MAKE_DIRICHLET_BC_FOR_NEUMANN)
        if( (alpha >= 1.0-mesh_radius-eps && alpha <= 1.0+eps && beta >= 1.0-eps) ||
            (alpha <= mesh_radius + eps && alpha >= -eps && beta <= eps) )
        {
          bnd[0] = 1.0;
          bnd[1] = 0.0;
          bnd[2] = gD;
          dirichlet_bc++;
        }
        else if(  beta >= 1.0 - mesh_radius - eps && beta <= 1.0 + eps && alpha >= 1.0 - eps ||
                  beta <= mesh_radius + eps && beta >= -eps && alpha <= eps)
        {
          bnd[0] = 1.0;
          bnd[1] = 0.0;
          bnd[2] = gD;
          dirichlet_bc++;
        }
        else
#endif
        {
          bnd[0] = 0.0;
          bnd[1] = 1.0;
          bnd[2] = gN;
          neumann_bc++;
        }
      }

      boundary_faces++;
    }
  }

  for(Mesh::iteratorElement f = m->BeginElement(EDGE|NODE); f != m->EndElement(); ++f)
  {
    f->Centroid(c);
    Storage::real alpha = c[0];//(c[0]-min[0])/(max[0]-min[0]);
	Storage::real beta  = c[1];//(c[1]-min[1])/(max[1]-min[1]);
    //~ Storage::real gamma = c[2];//(c[2]-min[2])/(max[2]-min[2]);
    Storage::real G = 0;
    if( flip )
    {
      G = exp(-2*pi*sqrt(1.0/parameter)*beta)*sin(2*pi*alpha)+sol_shift; //solution
      //dGdx = 2*pi*cos(2*pi*alpha)*exp(-2*pi*sqrt(1.0/parameter)*beta);
      //dGdy = -2*pi*sqrt(1.0/parameter)*sin(2*pi*alpha)*exp(-2*pi*sqrt(1.0/parameter)*beta);
    }
    else
    {
      G = exp(-2*pi*sqrt(parameter)*alpha)*sin(2*pi*beta)+sol_shift; //solution
      //dGdx = -2*pi*sqrt(parameter)*exp(-2*pi*sqrt(parameter)*alpha)*sin(2*pi*beta); //x derivative
      //dGdy = 2*pi*exp(-2*pi*sqrt(parameter)*alpha)*cos(2*pi*beta);
    }
    //Storage::real G = exp(-2*pi*sqrt(parameter)*alpha)*sin(2*pi*beta); //solution
    f->Real(solution_val) = G;

    if( f->Boundary() 
      && (alpha < eps || alpha > 1.0-eps || beta < eps || beta > 1.0 - eps)  ) //restrict top and bottom to be pure neumann bc
    {
      if( bc_type == 'a' )
      {
        Storage::real_array bnd = f->RealArray(bndcond);
        bnd[0] = 1.0;
        bnd[1] = 0.0;
        bnd[2] = G;
      }
      else if( bc_type == 'b' )
      {
        if( alpha < eps || beta < eps )
        {
          Storage::real_array bnd = f->RealArray(bndcond);
          bnd[0] = 1.0;
          bnd[1] = 0.0;
          bnd[2] = G;
          dirichlet_bc++;
        }
      }
#if defined(MAKE_DIRICHLET_BC_FOR_NEUMANN)
      else if( bc_type == 'c' )
      {
        if( alpha >= 1.0-mesh_radius-eps && alpha <= 1.0+eps && beta >= 1.0-eps )
        {
          Storage::real_array bnd = f->RealArray(bndcond);
          bnd[0] = 1.0;
          bnd[1] = 0.0;
          bnd[2] = G;
          dirichlet_bc++;
        }
        else if( beta >= 1.0 - mesh_radius - eps && beta <= 1.0 + eps && alpha >= 1.0 - eps )
        {
          Storage::real_array bnd = f->RealArray(bndcond);
          bnd[0] = 1.0;
          bnd[1] = 0.0;
          bnd[2] = G;
          dirichlet_bc++;
        }
      }
#endif
    }
  }

  std::cout << "Total faces " << total_faces << " boundary faces " << boundary_faces << std::endl;
  std::cout << "Dirichlet bc " << dirichlet_bc << " neumann bc " << neumann_bc << std::endl;
  


  std::cout << "Saving output to " << argv[2] << std::endl;

	m->Save(argv[2]);

	delete m;
}


