#include "inmost.h"

using namespace INMOST;


//const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "discontinuous_4zones";

int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out [lambda=100]" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try{m->Load(argv[1]);} catch(...) { std::cout << "Cannot load the mesh " << argv[1] << std::endl; return -1;}
	
	
	
	Mesh::GeomParam t;
	t[MEASURE] = FACE;
	t[ORIENTATION] = FACE;
	t[NORMAL] = FACE;
	t[CENTROID] = CELL | FACE | EDGE | NODE;
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
	
	double lambda = 100;
	if( argc > 3 ) lambda = atof(argv[3]);
	std::cout << "Lambda is: " << lambda << std::endl;
	
	if( m->HaveTag("PERM") ) m->DeleteTag(m->GetTag("PERM"));
	if( m->HaveTag("MATERIAL") ) m->DeleteTag(m->GetTag("MATERIAL"));
	
	
	Tag material = m->CreateTag("MATERIAL",DATA_INTEGER,CELL,NONE,1);
	Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,6);
	Tag solution_val = m->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL|FACE,NONE,1);
	Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,NONE,1);
	Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,3);

	{
		Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}
	
	
	
	
	for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE); it != m->EndElement(); ++it)
	{
		it->Centroid(c);
		unknown x(c[0],0);
		unknown y(c[1],1);
		unknown z(c[2],2);
		variable sol;
		rMatrix dsol(3,1), Kgrad(3,1);
		rMatrix K(3,3);
		int zone = -1;
		
		if( x < 0.5 )
		{
			if( y < 0.5 )
			{
				zone = 1;
				
				sol = 2*x+5*y+2*z;
				
				K(0,0) = 10;
				K(0,1) = K(1,0) = 2;
				K(1,1) = 20;
				K(0,2) = K(2,0) = 5;
				K(1,2) = K(2,1) = 4;
				K(2,2) = 20;
			}
			else
			{
				zone = 3;
				
				sol = 2*x + ((28-3*lambda)/(5.0*lambda))*y + 2*z + 14/5.0*(1-1/lambda);
				
				K(0,0) = 10*lambda;
				K(0,1) = K(1,0) = 2*lambda;
				K(1,1) = 20*lambda;
				K(0,2) = K(2,0) = 5*lambda;
				K(1,2) = K(2,1) = 4*lambda;
				K(2,2) = 20*lambda;
			}
		}
		else
		{
			if( y < 0.5 )
			{
				zone = 2;
				
				sol = 2*(2/lambda-1)*x + 5*y + 2*z + 2*(1-1/lambda);
				
				K(0,0) = 10*lambda;
				K(0,1) = K(1,0) = 2*lambda;
				K(1,1) = 20*lambda;
				K(0,2) = K(2,0) = 5*lambda;
				K(1,2) = K(2,1) = 4*lambda;
				K(2,2) = 20*lambda;
			}
			else
			{
				zone = 4;
				
				sol = 2*(2/lambda-1)*x + ((28-3*lambda)/(5.0*lambda))*y + 2*z + 24/5.0*(1-1/lambda);
				
				K(0,0) = 100*lambda;
				K(0,1) = K(1,0) = 25*lambda;
				K(1,1) = 200*lambda;
				K(0,2) = K(2,0) = (1219*lambda-2644)/10.0;
				K(1,2) = K(2,1) = 137*lambda - 606;
				K(2,2) = 200*lambda;
			}
		}
		/*
		std::cout << "sol: " << sol.GetValue() << std::endl;
		std::cout << "Jacobian: " << std::endl;
		sol.GetRow().Print();
		std::cout << "Hessian: " << std::endl;
		sol.GetHessianRow().Print();
		*/
		dsol(0,0) = sol.GetRow()[0];
		dsol(1,0) = sol.GetRow()[1];
		dsol(2,0) = sol.GetRow()[2];
		
		
		
		
		
		Kgrad = K*dsol;
		
		
		if( it->GetElementType() == CELL )
		{
			Storage::real_array perm = it->RealArray(tensor);
			perm[0] = K(0,0);
			perm[1] = K(0,1);
			perm[2] = K(0,2);
			perm[3] = K(1,1);
			perm[4] = K(1,2);
			perm[5] = K(2,2);
			
			it->Integer(material) = zone;
		}
		
		
		it->Real(solution_val) = sol.GetValue();

		if( it->GetElementType() == FACE )
		{
			//it->getAsFace()->OrientedUnitNormal(it->getAsFace()->BackCell(),nrm);
			it->getAsFace()->UnitNormal(nrm);
			
			double ff = -(Kgrad(0,0)*nrm[0] + Kgrad(1,0)*nrm[1] + Kgrad(2,0)*nrm[2]);
			
			
			it->Real(solution_flux) = ff;
				
			if( it->getAsFace()->Boundary() )
			{
				Storage::real_array bc = it->RealArray(bndcond);
				bc[0] = 1.0;
				bc[1] = 0.0;
				bc[2] = sol.GetValue();
			}
			
		}
	}
	
	std::cout << "Saving output to " << argv[2] << std::endl;
	
	m->Save(argv[2]);
	
	delete m;
}


