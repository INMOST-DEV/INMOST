#include "inmost.h"

using namespace INMOST;


//const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "linear_solution";

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
	
	if( m->HaveTag("PERM") ) m->DeleteTag(m->GetTag("PERM"));
	if( m->HaveTag("MATERIAL") ) m->DeleteTag(m->GetTag("MATERIAL"));
	
	
	Tag material = m->CreateTag("MATERIAL",DATA_INTEGER,CELL|EDGE|FACE|NODE,NONE,1);
	Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,6);
	Tag solution_val = m->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,1);
	Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,FACE,1);
	Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE|EDGE|NODE,FACE|EDGE|NODE,3);
	Tag force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,1);

	{
		Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}
	
	
	
	
	for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE|EDGE|NODE); it != m->EndElement(); ++it)
	{
		it->Centroid(c);
		unknown x(c[0],0);
		unknown y(c[1],1);
		unknown z(c[2],2);
		hessian_variable sol;
		vMatrix dsol(3,1), Kgrad(3,1);
		rMatrix K(3,3);
		
		sol = 1 + 2*x + 3*y + 4*z;
		K(0,0) = 1;
		K(1,1) = 2;
		K(2,2) = 3;
		K(0,1) = K(1,0) = -0.1;
		K(0,2) = K(2,0) = -0.5;
		K(1,2) = K(2,1) = -0.2;
		/*
		std::cout << "sol: " << sol.GetValue() << std::endl;
		std::cout << "Jacobian: " << std::endl;
		sol.GetRow().Print();
		std::cout << "Hessian: " << std::endl;
		sol.GetHessianRow().Print();
		*/
		dsol(0,0) = sol.GetVariable(0);
		dsol(1,0) = sol.GetVariable(1);
		dsol(2,0) = sol.GetVariable(2);
		
		
		
		
		
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
			
			it->Real(force) = -(Kgrad(0,0).GetRow()[0] + Kgrad(1,0).GetRow()[1] + Kgrad(2,0).GetRow()[2]);
		}
		
		
		it->Real(solution_val) = sol.GetValue();

		if( it->GetElementType() == FACE )
		{
			//it->getAsFace()->OrientedUnitNormal(it->getAsFace()->BackCell(),nrm);
			it->getAsFace()->UnitNormal(nrm);
			
			double ff = -(Kgrad(0,0).GetValue()*nrm[0] + Kgrad(1,0).GetValue()*nrm[1] + Kgrad(2,0).GetValue()*nrm[2]);
			
			
			it->Real(solution_flux) = ff;
				
			if( it->getAsFace()->Boundary() )
			{
				Storage::real_array bc = it->RealArray(bndcond);
				bc[0] = 1.0;
				bc[1] = 0.0;
				bc[2] = sol.GetValue();
			}
			
		}
		if( it->GetElementType() & (EDGE | NODE) )
		{
			if( it->Boundary() )
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


