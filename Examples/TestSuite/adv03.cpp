#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


using namespace INMOST;
std::string problem_name = "advection_degenerate_diffusion_reaction";

const Storage::real pi = 3.1415926535897932384626433832795;
#define V_ID(x, y, z) ((x)*n*n + (y)*n + (z))

int main(int argc, char *argv[]) 
{

	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " nx|grid_name [alpha=0.4] [outer_boundary_pressure=0.0] [inner_boundary_pressure=1.0] [cut_grid=1] [is2d=1]" << std::endl;
		return -1;
	}

	Mesh * mesh;
	double alpha=0.2;
	double h;
	double outer_boundary_pressure = 0.0;
	double inner_boundary_pressure = 1.0;
	int cut_grid = 1;
	int is2d = 0;
	int n = 20 + 1;

	mesh = new Mesh();

	if (argc > 1) 
	{
		if( atoi(argv[1]) )
		{
			n = atoi(argv[1])+1;
			h = 1.0 / (n-1);
			std::cout << "Mesh: cube " << n << std::endl;
		}
		else
		{
			mesh->Load(argv[1]);
			std::cout << "Mesh: " << argv[1] << std::endl;
			n = 0;
		}

	}

	std::cout << "Mesh radius: " << h << std::endl;
	if( argc > 2 )
		alpha = atof(argv[2]);

	std::cout << "Alpha: " << alpha << std::endl;

	if( argc > 3 )
		outer_boundary_pressure = atof(argv[3]);

	if( argc > 4 )
		inner_boundary_pressure = atof(argv[4]);

	std::cout << "Boundaries: " << outer_boundary_pressure << " " << inner_boundary_pressure << std::endl;

	if( argc > 5 )
		cut_grid = atoi(argv[5]);

	if( cut_grid )
		std::cout << "Cutting center of the grid." << std::endl;

	if( argc > 6 )
		is2d = atoi(argv[6]);

	if( is2d )
		std::cout << "Using 2d setup" << std::endl;
	else
		std::cout << "Using 3d setup" << std::endl;

	double ratio = 1000;
	if( argc > 7 )
		ratio = atof(argv[7]);

	std::cout << "Anisotropy ratio: " << ratio << std::endl;


	srand(0);//(unsigned)time(NULL));

	if( n )
	{
		for (int i = 0; i < n; i++) 
		{
			for (int j = 0; j < n; j++) 
			{
				for (int k = 0; k < n; k++) 
				{
					Storage::real xyz[3];
					xyz[0] = i * 1.0 / (n - 1);
					xyz[1] = j * 1.0 / (n - 1);
					xyz[2] = k * 1.0 / (n - 1);
					Node c = mesh->CreateNode(xyz);
					if (c->LocalID() != V_ID(i, j, k)) std::cout << "v_id = " << c->LocalID() << ", [i,j,k] = " << V_ID(i,j,k) << std::endl;
				}
			}
		}


		const INMOST_DATA_INTEGER_TYPE nvf[24] = { 2, 3, 1, 0, 4, 5, 7, 6, 0, 1, 5, 4, 3, 2, 6, 7, 2, 0, 4, 6, 1, 3, 7, 5 };
		const INMOST_DATA_INTEGER_TYPE numnodes[6] = { 4, 4, 4, 4, 4, 4 };
		for (int i = 1; i < n; i++) 
		{
			for (int j = 1; j < n; j++) 
			{
				for (int k = 1; k < n; k++) 
				{

					ElementArray<Node> verts(mesh,8); 
					verts[0] = mesh->NodeByLocalID(V_ID(i - 1, j - 1, k - 1));
					verts[1] = mesh->NodeByLocalID(V_ID(i - 0, j - 1, k - 1));
					verts[2] = mesh->NodeByLocalID(V_ID(i - 1, j - 0, k - 1));
					verts[3] = mesh->NodeByLocalID(V_ID(i - 0, j - 0, k - 1));
					verts[4] = mesh->NodeByLocalID(V_ID(i - 1, j - 1, k - 0));
					verts[5] = mesh->NodeByLocalID(V_ID(i - 0, j - 1, k - 0));
					verts[6] = mesh->NodeByLocalID(V_ID(i - 1, j - 0, k - 0));
					verts[7] = mesh->NodeByLocalID(V_ID(i - 0, j - 0, k - 0));

					mesh->CreateCell(verts,nvf,numnodes,6);
				}
			}
		}
	}


	if( cut_grid )
	{
		for (Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
		{
			Storage::real cnt[3];
			it->Centroid(cnt);
			if( cnt[0] > 0.25 && cnt[0] < 0.75 && cnt[1] > 0.25 && cnt[1] < 0.75 )
				it->Delete();
			
		}

		for (Mesh::iteratorElement it = mesh->BeginElement(FACE|EDGE|NODE); it != mesh->EndElement(); ++it)
			if( it->nbAdjElements(CELL) == 0 ) it->Delete();
	}

	if( alpha > 0.0 ) //skewness
	{
		const double eps = 1.0e-6;
		for( Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
		{
			Storage::real_array coords = it->Coords();
			Storage::real h = 1.0e20;

			ElementArray<Cell> it_cells = it->getCells();

			for(INMOST_DATA_ENUM_TYPE k = 0; k < it_cells.size(); ++k)
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
				h = std::min<Storage::real>(h,sqrt((maxmin[2]-maxmin[3])*(maxmin[2]-maxmin[3])+(maxmin[0]-maxmin[1])*(maxmin[0]-maxmin[1])));
			}
			if( coords[0] > eps && coords[0] < 1.0-eps && coords[1] > eps && coords[1] < 1.0-eps )
			{
				if( !(coords[0] > 4.0/9.0-eps && coords[0] < 5.0/9.0+eps && coords[1] > 4.0/9.0-eps && coords[1] < 5.0/9.0+eps) )
				{
					coords[0] += alpha*h*(2.0*rand()/RAND_MAX-1.0);
					coords[1] += alpha*h*(2.0*rand()/RAND_MAX-1.0);
				}
			}
		}
	}

	for( Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
	{
		Storage::real_array coords = it->Coords();
		coords[0] = 2*coords[0]-1;
		coords[1] = 2*coords[1]-1;
	}

	mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);

	{ // Prepare geometrical data on the mesh.
		Mesh::GeomParam table;
		table[CENTROID]    = CELL | FACE; //Compute averaged center of mass
		table[NORMAL]      = FACE;        //Compute normals
		table[ORIENTATION] = FACE;        //Check and fix normal orientation
		table[MEASURE]     = CELL | FACE; //Compute volumes and areas
		//table[BARYCENTER]  = CELL | FACE; //Compute volumetric center of mass
		mesh->RemoveGeometricData(table);
		mesh->PrepareGeometricData(table); //Ask to precompute the data
	}

	std::cout << "nodes: " << mesh->NumberOfNodes() << " edges: " << mesh->NumberOfEdges() << " faces: " << mesh->NumberOfFaces() << " cells: " << mesh->NumberOfCells() << std::endl;


	{
		Storage::bulk_array name = mesh->self()->BulkArray(mesh->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}

	Tag bndcond = mesh->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,3);
	Tag velocity = mesh->CreateTag("VELOCITY",DATA_REAL,FACE,NONE,1);
	Tag reaction = mesh->CreateTag("REACTION",DATA_REAL,CELL,NONE,1);
	Tag tensor = mesh->CreateTag("PERM", DATA_REAL, CELL, NONE, 1);
	Tag force = mesh->CreateTag("FORCE",DATA_REAL,CELL,NONE,1);
	Tag vel3 = mesh->CreateTag("VELVEC",DATA_REAL,CELL,NONE,3);
	Tag refsol = mesh->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL,NONE,1);
	Tag refflux = mesh->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,NONE,1);
	Tag zone = mesh->CreateTag("ZONE",DATA_INTEGER,CELL,NONE,1);
	int numinner = 0, numouter = 0;
	int numinnern = 0, numoutern = 0;
	const Storage::real eps = 1.0e-6;
	const Storage::real mu = 1.0e-3;
	const Storage::real velmult = 1;
	for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); ++it) 
	{
		Storage::real cnt[3], nrm[3];
		it->Centroid(cnt);
		it->UnitNormal(nrm);
		Storage::real x = cnt[0];
		Storage::real y = cnt[1];
		Storage::real z = cnt[2];
		Storage::real r = sqrt(x*x+y*y);
		Storage::real theta = atan2(y,x);//+pi;
		Storage::real sol, diff, flux, velnrm, dsolx, dsoly;
		Storage::real vel[3] = {-y/(r*r)*velmult,x/(r*r)*velmult,0.0};

		if( theta < 0 ) theta = 2*pi + theta;
		
		velnrm = nrm[0]*vel[0] + nrm[1]*vel[1] + nrm[2]*vel[2];

		if( theta < pi )
		{
			diff = pi;
			sol = (theta-pi)*(theta-pi);
			dsolx = -2*y*(theta-pi)/(r*r);
			dsoly = 2*x*(theta-pi)/(r*r);
			flux = velnrm*sol - diff*(dsolx*nrm[0] + dsoly*nrm[1]);
		}
		else
		{
			diff = 0;
			sol = 3*pi*(theta-pi);
			flux = velnrm*sol;
		}
		
		it->Real(refflux) = flux;
		it->Real(velocity) = velnrm;

		if( it->Boundary() && z > 0+eps && z < 1-eps)
		{
			Storage::real_array bc = it->RealArray(bndcond);
			bc[0] = 1;
			bc[1] = 0;
			bc[2] = sol;
		}
	}
	std::cout << "Outer faces " << numouter << " inner faces " << numinner << std::endl;
	std::cout << "Outer nodes " << numoutern << " inner nodes " << numinnern << std::endl;

	
	
	for (Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
	{
		Storage::real cnt[3];
		it->Centroid(cnt);
		Storage::real f = 0;
		Storage::real x = cnt[0];
		Storage::real y = cnt[1];
		Storage::real r = sqrt(x*x+y*y);
		Storage::real theta = atan2(y,x);//+pi;
		Storage::real sol, diff;
		Storage::real vel[3] = {-y/(r*r)*velmult,x/(r*r)*velmult,0.0};

		if( theta < 0 ) theta = 2*pi + theta;
		
		it->Centroid(cnt);
		if( theta < pi )
		{
			diff = pi;
			//diff = 0;
			sol = (theta-pi)*(theta-pi);
			f = mu*sol + 2*(theta-pi)/(r*r)*velmult - diff*2/(r*r);
			
			it->Integer(zone) = 0;
		}
		else 
		{
			diff = 0;
			sol = 3*pi*(theta-pi);
			f = mu*sol + 3*pi/(r*r)*velmult;
			
			it->Integer(zone) = 1;
		}

		it->Real(force) = f;
		it->Real(reaction) = mu;
		it->Real(refsol) = sol;
		it->Real(tensor) = diff;

		Storage::real_array svel = it->RealArray(vel3);
		svel[0] = vel[0];
		svel[1] = vel[1];
		svel[2] = vel[2];
	}

	std::cout << "I'm ready!" << std::endl;

	//mesh->Save("grid.vtk");
	mesh->Save("grid_out.pmf");
	//mesh->Save("grid.gmv");

	std::cout << "File written!" << std::endl;

	delete mesh;
	return 0;
}
