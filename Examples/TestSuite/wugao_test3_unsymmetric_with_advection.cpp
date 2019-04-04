#include "inmost.h"

using namespace INMOST;


const Storage::real eps = 1.0e-6;
std::string problem_name = "wugao_test3_unsymmetric";

int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out [a:1 b:2]" << std::endl;
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

	Storage::real a = 1, b = 2;

	if( argc > 3 ) a = atof(argv[3]);
	if( argc > 4 ) b = atof(argv[4]);

	
	Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,9);
	Tag solution_val = m->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,1);
	Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,FACE,1);
	Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,3);
	Tag velocity, force;

	{
		Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}

	double velmult = 10.0;

	if( velmult )
	{
		force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,1);
		velocity = m->CreateTag("VELOCITY",DATA_REAL,FACE,NONE,1);
	}


	for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE|EDGE|NODE); it != m->EndElement(); ++it)
	{
		it->Centroid(c);
		Storage::real x = c[0];//(c[0]-min[0])/(max[0]-min[0]);
		Storage::real y = c[1];//(c[1]-min[1])/(max[1]-min[1]);
		Storage::real z = c[2];//(c[2]-min[2])/(max[2]-min[2]);

		Storage::real vel[3] = {-(y-0.5)*velmult,(x-0.5)*velmult,0.0};

		if( it->GetElementType() == CELL )
		{
			Storage::real_array perm = it->RealArray(tensor);
			perm[0] = a;
			perm[1] = b*y;
			perm[2] = 0.0;
			perm[3] = -b*y;
			perm[4] = a;
			perm[5] = 0.0;
			perm[6] = 0.0;
			perm[7] = 0.0;
			perm[8] = 1.0;
		}





		Storage::real sol, dsolx, dsoly = 0;
		Storage::real flux;

		sol = (exp(b/a*x)-1)/(exp(b/a)-1);
		dsolx = (b*exp((b*x)/a))/(a*(exp(b/a)-1));

		it->Real(solution_val) = sol;

		if( it->GetElementType() == CELL && velmult )
		{
			it->Real(force) = (dsolx*vel[0] + dsoly*vel[1]);
		}

		if( it->GetElementType() == FACE )
		{
			it->getAsFace()->UnitNormal(nrm);

			double adv_flux = 0;

			if( velmult )
			{
				double nvel = vel[0]*nrm[0] + vel[1]*nrm[1] + vel[2]*nrm[2];
				it->Real(velocity) = nvel;
				adv_flux = sol*nvel;
			}

			flux = -((a*nrm[0]-b*y*nrm[1])*dsolx + (+b*y*nrm[0]+a*nrm[1])*dsoly) + adv_flux;

			if( !(z < eps || z > 1.0-eps ) )
			{
				it->Real(solution_flux) = flux;

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


