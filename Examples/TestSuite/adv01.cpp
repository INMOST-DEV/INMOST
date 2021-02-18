#include "inmost.h"

using namespace INMOST;

typedef Storage::real real;
typedef Storage::real_array real_array;
typedef Storage::bulk_array bulk_array;


const real pi = 3.1415926535897932384626433832795;
const real eps = 1.0e-6;
std::string problem_name = "zalesak_disk";

real max[3] = {-1.0e20, -1.0e20, -1.0e20}, min[3] = {1.0e20,1.0e20,1.0e20};
real disc_cx = 0.5;
real disc_cy = 0.75;
real disc_rad = 15.0/100.0;
real slot_top = disc_cy + 0.1;
real slot_w = 5.0/100.0;
real slot_l = disc_cx - slot_w*0.5;
real slot_r = disc_cx + slot_w*0.5;

real zalesak_shape(real * coord, real time)
{
	(void)time;
	real x = (coord[0]-min[0])/(max[0]-min[0]); //normalize x into [0,1]
	real y = (coord[1]-min[1])/(max[1]-min[1]); //normalize y into [0,1]

	if( sqrt((x-disc_cx)*(x-disc_cx) + (y-disc_cy)*(y-disc_cy)) < disc_rad )
	{
		if( y < slot_top && x > slot_l && x < slot_r )
			return 1.0;
		return 3.0;
	}
	return 1.0;
}


int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out cfl" << std::endl;
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

	
	real c[3] = {0,0,0}, nrm[3] = {0,0,0};
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

	real mesh_radius = 0;
	for(Mesh::iteratorCell e = m->BeginCell(); e != m->EndCell(); ++e)
	{
		real maxmin[6];
		maxmin[0] = -1e20;
		maxmin[1] = 1e20;
		maxmin[2] = -1e20;
		maxmin[3] = 1e20;
		maxmin[4] = -1e20;
		maxmin[5] = 1e20;
		ElementArray<Node> nodes = e->getNodes();
		for (ElementArray<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
		{
			real_array c = it->Coords();
			for (int i = 0; i < (int)c.size(); i++) 
			{
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

	std::cout << "Mesh bounds: " << min[0] << ":" << max[0] << " " << min[1] << ":" << max[1] << " " << min[2] << ":" << max[2] << " radius " << mesh_radius << std::endl;

	Tag material;
	//if( m->HaveTag("PERM") ) m->DeleteTag(m->GetTag("PERM"));


	Tag time = m->CreateTag("TIME_INFO",DATA_REAL,MESH,NONE,2); //first entry - time step, second entry - total time
	Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,3); //boundary conditions
	Tag velocity = m->CreateTag("VELOCITY",DATA_REAL,FACE,NONE,1); //normal component of the velocity
	Tag conc = m->CreateTag("CONCENTRATION",DATA_REAL,CELL,NONE,1); //initial concentration

	{
		bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}

	real_array time_info = m->self()->RealArray(time);

	time_info[0] = 1.0/628.0;// * 14.277739148561303614695220241265 / 2.5; // cfl 10 on cartesian mesh 1/100
	time_info[1] = 1;
	real dt = time_info[0], cfl = 0;

	real origin_x = (max[0]+min[0])*0.5;
	real origin_y = (max[1]+min[1])*0.5;
	real omega = 2*pi; //one circle in one second

	

	std::cout << "Origin axis: (" << origin_x << "," << origin_y << ") angular velocity " << omega << std::endl;

	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it)
	{
		it->Centroid(c);
		it->UnitNormal(nrm);

		real v = omega*(c[0] - origin_x);
		real u = omega*(origin_y - c[1]);

		//fill in normal component of the velocity
		it->Real(velocity) = nrm[0]*u + nrm[1]*v;
		
		if( it->Boundary() )
		{
			real_array bc = it->RealArray(bndcond); //create bc
			//dirichlet bc with concentration 1
			bc[0] = 1;
			bc[1] = 0;
			bc[2] = 1;
		}
	}

	Tag velrec = m->CreateTag("VELREC",DATA_REAL,CELL,NONE,3);
	real Vmax = 0.0;
	for(Mesh::iteratorCell it  = m->BeginCell(); it != m->EndCell(); ++it)
	{
		real_array V = it->RealArray(velrec);
		real xK[3], xF[3], sgn, A, U, Vol = it->Volume();
		ElementArray<Face> faces = it->getFaces();
		it->Centroid(xK);
		for(int k = 0; k < (int)faces.size(); ++k)
		{
			faces[k].Centroid(xF);
			A = faces[k].Area();
			U = faces[k].Real(velocity);
			sgn = (faces[k].FaceOrientedOutside(it->self()) ? 1 : -1);
			V[0] += sgn*A*U*(xF[0]-xK[0]);
			V[1] += sgn*A*U*(xF[1]-xK[1]);
			V[2] += sgn*A*U*(xF[2]-xK[2]);
		}
		V[0] /= Vol;
		V[1] /= Vol;
		V[2] /= Vol;
		real mag = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
		Vmax = std::max(mag,Vmax);
	}
	if( argc > 3 )
	{
		cfl = atof(argv[3]);
		dt = cfl*mesh_radius/Vmax;
		time_info[0] = dt;
		std::cout << "Selected time step: " << dt << " based on CFL" << std::endl;
	}
	else cfl = Vmax*dt/mesh_radius;
	std::cout << "CFL: " << cfl << std::endl;

	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		it->Real(conc) = it->Mean(zalesak_shape,1.0); //initial concentration
	
	std::cout << "Saving output to " << argv[2] << std::endl;

	//m->Save("test.vtk");
	m->Save(argv[2]);

	delete m;
}


