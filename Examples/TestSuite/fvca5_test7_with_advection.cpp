#include "inmost.h"

using namespace INMOST;


//const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "oblique_barrier";

int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out [sigma:0.2 incline:1 y_low:0.25 y_high:0.5]" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try{m->Load(argv[1]);} catch(...) { std::cout << "Cannot load the mesh " << argv[1] << std::endl; return -1;}
	
	Storage::real sigma = 0.2, y_low = 0.25, y_high = 0.5;
	
	if( argc > 3 ) sigma = atof(argv[3]);
	
	//for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); ++t)
	//  if( t->GetTagName().substr(0,9) == "GEOM_UTIL" ) m->DeleteTag(*t);
	
	if( argc <= 4 || (argc > 4 && atoi(argv[4])) )
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

		std::cout << "You asked to incline the grid to honour discontinuities" << std::endl;
		std::cout << "and I assume that the initial grid is rectangular." << std::endl;
		for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
		{
			
			Storage::real_array c = it->Coords();
			Storage::real x = c[0];
			Storage::real y = c[1];
			Storage::real d = 0.0;
			Storage::real y1 = sigma*(x-0.5)+0.475;
			Storage::real y2 = sigma*(x-0.5)+0.475+0.05;
			if( y <= y_low ) d = y/y_low*(y1-y_low);
			else if( y >= y_low && y <= y_high ) d = (y1-y_low) + (y-y_low)/(y_high-y_low)*((y2-y_high)-(y1-y_low));
			else d = (1-(y-y_high)/(1-y_high))*(y2-y_high);
			c[1] += d;

			
		}


		//introduce skewness
		/*
		for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
		{
			Storage::real_array c = it->Coords();
			Storage::real x = c[0];
			Storage::real y = c[1];
			Storage::real y1 = sigma*(x-0.5)+0.475;
			Storage::real y2 = sigma*(x-0.5)+0.475+0.05;
			Storage::real rndx = (rand()*1.0/RAND_MAX*2-1);
			Storage::real rndy = (rand()*1.0/RAND_MAX*2-1);
			Storage::real h = 1.0e20;

			ElementArray<Cell> it_cells = it->getCells();

			for(int k = 0; k < it_cells.size(); ++k)
			{
				Storage::real maxmin[6];
				maxmin[0] = -1e20;
				maxmin[1] = 1e20;
				maxmin[2] = -1e20;
				maxmin[3] = 1e20;
				ElementArray<Node> nodes = it_cells[k].getNodes();
				for (ElementArray<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
				{
					Storage::real_array cc = it->Coords();
					for (int i = 0; i < (int)cc.size(); i++) 
					{
						if (maxmin[2 * i + 0] < cc[i]) maxmin[2 * i + 0] = cc[i]; //max
						if (maxmin[2 * i + 1] > cc[i]) maxmin[2 * i + 1] = cc[i]; //min
					}
				}
				h = std::min(h,sqrt((maxmin[2]-maxmin[3])*(maxmin[2]-maxmin[3])+(maxmin[0]-maxmin[1])*(maxmin[0]-maxmin[1])));
			}


			if( fabs(y-y1) > 1.0e-9 && fabs(y-y2) > 1.0e-9 && x > 0 && x < 1 && y > 0 && y < 1 )
			{
				c[0] += rndx*h*alpha;
				c[1] += rndy*h*alpha;
			}
		}
		*/
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
	
	if( m->HaveTag("PERM") ) m->DeleteTag(m->GetTag("PERM"));
	if( m->HaveTag("MATERIAL") ) m->DeleteTag(m->GetTag("MATERIAL"));
	
	
	Tag material = m->CreateTag("MATERIAL",DATA_INTEGER,CELL|EDGE|FACE|NODE,NONE,1);
	Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,3);
	Tag solution_val = m->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,1);
	Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,FACE,1);
	Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE|EDGE|NODE,FACE|EDGE|NODE,3);
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
			//~ Storage::real r = 1;//sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
		//Storage::real vel[3] = {-(y-0.5)/(r*r)*velmult,(x-0.5)/(r*r)*velmult,0.0};
		Storage::real vel[3] = {-0.5*velmult,velmult,0.0};

		
		
		int zone = 0;
		
		double phi1 = y - sigma*(x-0.5) - 0.475;
		double phi2 = phi1 - 0.05;
		double alpha;
		
		if( phi1 < 0 ) zone = 1;
		else if( phi1 >= 0 && phi2 < 0 ) zone = 2;
		else zone = 3;
		//else printf("%s:%d oops!\n",__FILE__,__LINE__);
		
		it->IntegerDF(material) = zone;
		
		Storage::real sol, dsolx, dsoly;
		Storage::real flux;
		
		if( zone == 2 )
		{
			alpha = 0.01;
		}
		else
		{
			alpha = 1;
		}
		sol = 0;
		if( zone == 1 )
		{
			sol = -phi1;
			dsolx = sigma;
			dsoly = -1;
		}
		else if( zone == 2 )
		{
			sol = -phi1/0.01;
			dsolx = sigma/0.01;
			dsoly = -1/0.01;
		}
		else if( zone == 3 )
		{
			sol = -phi2 - 0.05/0.01;
			dsolx = sigma;
			dsoly = -1;
		}
		
		sol += 6; //pick solution up for positivity
		
		Storage::real K[6] =
		{
			alpha,
			0.0,
			0.0,
			alpha,
			0,
			1
		};
		
		
		if( it->GetElementType() == CELL )
		{
			Storage::real_array perm = it->RealArray(tensor);
			perm[0] = alpha;
			perm[1] = alpha;
			perm[2] = 1;
		}
		
		
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
			
			flux = -((K[0]*nrm[0]+K[1]*nrm[1])*dsolx + (K[1]*nrm[0]+K[3]*nrm[1])*dsoly) + adv_flux;
			
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
		if( it->GetElementType() & (EDGE | NODE) )
		{
			if( (x < eps || x > 1.0-eps || y < eps || y > 1.0-eps) )
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
	}
	
	std::cout << "Saving output to " << argv[2] << std::endl;
	
	m->Save(argv[2]);
	
	delete m;
}


