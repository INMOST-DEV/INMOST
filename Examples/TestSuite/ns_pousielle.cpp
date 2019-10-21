#include "inmost.h"

using namespace INMOST;

//TODO:
//attach scale to rescale grid
//attach slice grid to cut cylinder
//option to set shear rate / pressure drop / maximum velocity
//setup BC
//cmake compilation

//setup on mesh:
//
// FORCE - 3 entries
//
// BOUNDARY_CONDITION_VELOCITY - 7 entries
// r = (a5,a6,a7)
// n.(a1*u + a2*t) = n.r
// (I - nn.)(a3*u + a4*t) = (I-nn.)r
//
// BOUNDARY_CONDITION_PRESSURE - 1 entry
//
// REFERENCE_VELOCITY - 3 entries
//
// REFERENCE_PRESSURE - 1 entry
//



typedef Storage::real real;



int main(int argc, char ** argv)
{
	
	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " mesh [axis=2 (0:x,1:y,2:z)] [control=0 (0:shear rate,1:pressure drop,2:maximal velocity)] [control_value=10] [visc=1.0e-5] [mesh_out=grid_out.pmf]" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try
	{
		m->Load(argv[1]);
	}
	catch(...)
	{
		std::cout << "Cannot load the mesh " << argv[1] << std::endl;
		return -1;
	}
	
	std::string fout = "grid_out.pmf";
	int axis = 2;
	double visc = 1.0e-5;
	int control = 0;
	double control_value = 10;
	double Len = 8;
	double Diam = 1;
	double shear = 10;//0; //shear rate
	double dp = 0;//36e-6;
	double vmax = 0;

	if( argc > 2 ) axis = atoi(argv[2]);
	if( argc > 3 ) control = atoi(argv[3]);
	if( argc > 4 ) control_value = atof(argv[4]);
	if( argc > 5 ) visc = atof(argv[5]);
	if( argc > 6 ) fout = std::string(argv[6]);
	
	
	
	if( axis < 0 || axis > 2 )
	{
		std::cout << "bad axis: " << axis << " should be 0 - x, 1 - y, 2 - z" << std::endl;
		return -1;
	}
	
	if( control < 0 || control > 2 )
	{
		std::cout << "bad control: " << control << " should be 0 - shear rate, 1 - pressure drop, 2 - maximal velocity" << std::endl;
		return -1;
	}
	
	if( visc <= 0 )
	{
		std::cout << "bad viscosity: " << visc << std::endl;
		return -1;
	}
	
	if( control == 0 )
	{
		shear = control_value;
		dp = vmax = 0;
	}
	else if( control == 1 )
	{
		dp = control_value;
		shear = vmax = 0;
	}
	else if( control == 2 )
	{
		vmax = control_value;
		dp = shear = 0;
	}
	
	
	double cmax[3] = {-1.0e20,-1.0e20,-1.0e20}, cmin[3] = {1.0e20,1.0e20,1.0e20};
	for(Mesh::iteratorNode n = m->BeginNode(); n != m->EndNode(); ++n)
	{
		Storage::real_array c = n->Coords();
		for(int k = 0; k < 3; ++k)
		{
			if( cmax[k] < c[k] ) cmax[k] = c[k];
			if( cmin[k] > c[k] ) cmin[k] = c[k];
		}
	}
	Len = cmax[axis] - cmin[axis];
	Diam = std::max(cmax[(axis+1)%3]-cmin[(axis+1)%3],cmax[(axis+2)%3]-cmin[(axis+2)%3]);
			
	
	if( shear )
	{
		dp = 4*Len/Diam*shear*visc;
		vmax = Diam*Diam*dp/(16*visc*Len);
	}
	else if( dp )
	{
		vmax = Diam*Diam*dp/(16*visc*Len);
		shear = 4*vmax/Diam;
	}
	else if( vmax )
	{
		shear = 4*vmax/Diam;
		dp = 4*Len/Diam*shear*visc;
	}
	
	std::cout << "viscosity " << visc << " shear " << shear << " dp " << dp << " vmax " << vmax << " length " << Len << " diameter " << Diam << std::endl;
	
	//TagRealArray force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,3);
	TagRealArray  bc      = m->CreateTag("BOUNDARY_CONDITION_VELOCITY",DATA_REAL,FACE,FACE,7);
	TagReal       bcp     = m->CreateTag("BOUNDARY_CONDITION_PRESSURE",DATA_REAL,FACE,FACE,1);
	TagRealArray  uvw     = m->CreateTag("UVW",DATA_REAL,CELL,NONE,3);
	TagRealArray  uvw_ref = m->CreateTag("REFERENCE_VELOCITY",DATA_REAL,CELL,NONE,3); //depends on viscosity
	TagReal       p_ref   = m->CreateTag("REFERENCE_PRESSURE",DATA_REAL,CELL,NONE,1);
	TagReal       p       = m->CreateTag("P",DATA_REAL,CELL,NONE,1);
	
	
	
	//this should not be needed?
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
		it->FixNormalOrientation();
	
	double cnt0[3];
	cnt0[axis] = 0;
	cnt0[(axis+1)%2] = (cmax[(axis+1)%2] + cmin[(axis+1)%2])*0.5;
	cnt0[(axis+2)%2] = (cmax[(axis+2)%2] + cmin[(axis+2)%2])*0.5;
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		double cnt[3];
		it->Barycenter(cnt);
		double r = sqrt(pow(cnt[(axis+1)%2]-cnt0[(axis+1)%2],2)+pow(cnt[(axis+2)%2]-cnt0[(axis+2)%2],2));
		uvw(*it,3,1).Zero();
		p[*it] = 0;
		uvw_ref(*it,3,1).Zero();
		uvw_ref(*it,3,1)(axis,0) = (Diam*Diam/4.0-r*r)*dp/(4.0*visc*Len);
		p_ref[*it] = dp*(cnt[axis]-cmin[axis])/Len;
	}
	
	
	
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
	{
		double  n[3];
		it->UnitNormal(n);
		if( fabs(n[axis]-1) < 1.0e-3 ) // outflow
		{
			bcp[*it] = 0+10;
		}
		else if(  fabs(n[axis]+1) < 1.0e-3 ) //inflow
		{
			bcp[*it] = dp+10;
		}
		else //no-slip walls
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 1;
			bc[*it][3] = 0;
			bc[*it][4] = 0;
			bc[*it][5] = 0;
			bc[*it][6] = 0;
		}
	}
	
	std::cout << "Saving output to " << fout << std::endl;
	
	m->Save(fout);
	
	return 0;
}
