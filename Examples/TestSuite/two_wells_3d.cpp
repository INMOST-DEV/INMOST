#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


using namespace INMOST;
std::string problem_name = "two_wells_3d";

#define V_ID(x, y, z) ((x)*n*n + (y)*n + (z))

void matmul(Storage::real * a, Storage::real * b, Storage::real * out)
{
	int i,j,k;
	double ret[9] = {0,0,0,0,0,0,0,0,0};
	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			for(k = 0; k < 3; k++)
			{
				ret[i*3+j] += a[i*3+k]*b[k*3+j];
			}
	for(i = 0; i < 9; i++) out[i] = ret[i];
}

void multangle(Storage::real t[9], Storage::real anglex, Storage::real angley, Storage::real anglez)
{
	Storage::real rot[9];
	rot[0] = cos(anglex);	rot[1] = -sin(anglex); rot[2] = 0.0;
	rot[3] = -rot[1];     rot[4] = rot[0];      rot[5] = 0.0;
	rot[6] = 0.0;         rot[7] = 0.0;         rot[8] = 1.0;
	
	matmul(rot,t,t);
	
	rot[1] = -rot[1];
	rot[3] = -rot[3];
	
	matmul(t,rot,t);
	
	
	rot[0] = cos(angley); rot[1] = 0.0;	rot[2] = -sin(angley);
	rot[3] = 0.0;         rot[4] = 1.0; rot[5] = 0.0;
	rot[6] = -rot[2];     rot[7] = 0.0;	rot[8] = rot[0];
	
	matmul(rot,t,t);
	
	rot[2] = -rot[2];
	rot[6] = -rot[6];
	
	matmul(t,rot,t);
	
	
	rot[0] = 1.0; rot[1] = 0.0;	        rot[2] = 0.0;
	rot[3] = 0.0;	rot[4] = cos(anglez);	rot[5] = -sin(anglez);
	rot[6] = 0.0;	rot[7] = -rot[5];	    rot[8] = rot[4];
	
	matmul(rot,t,t);
	
	rot[5] = -rot[5];
	rot[7] = -rot[7];
	
	matmul(t,rot,t);
}

void prepare_tensor(Storage::real t[9], Storage::real ratio)
{
	Storage::real k[9] = { ratio, 0, 0, 0, 15, 0, 0, 0, 1 };
	Storage::real anglex = 17, angley = 35.5, anglez = 67.5;
	anglex = anglex * 3.1415926535897932384626433832795 / 180.0;
	angley = angley * 3.1415926535897932384626433832795 / 180.0;
	anglez = anglez * 3.1415926535897932384626433832795 / 180.0;
	multangle(k,anglex,angley,anglez);
	memcpy(t,k,sizeof(Storage::real)*9);
	
}


void prepare_tensor_2d(Storage::real t[9], Storage::real ratio)
{
	Storage::real k[9] = { ratio, 0, 0, 0, 1, 0, 0, 0, 1 };
	Storage::real anglex = 67.5, angley = 0, anglez = 0;
	anglex = anglex * 3.1415926535897932384626433832795 / 180.0;
	multangle(k,anglex,angley,anglez);
	memcpy(t,k,sizeof(Storage::real)*9);
	
}

int main(int argc, char *argv[])
{
	
	if (argc < 2)
	{
  std::cout << "Usage: " << argv[0] << " nx|grid_name [alpha=0.4] [left_well_pressure=0.0] [right_well_pressure=1.0] [cut_grid=1] [is2d=1] [anisotropy=1000]" << std::endl;
  return -1;
	}
	
	Mesh * mesh;
	double alpha=0.4;
	Storage::real h;
	double left_well = 0.0;
	double right_well = 1.0;
	int cut_grid = 1;
	int is2d = 1;
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
			h = 0;
			for(Mesh::iteratorCell e = mesh->BeginCell(); e != mesh->EndCell(); ++e)
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
				h = std::max(h,maxmin[0]-maxmin[1]);
				h = std::max(h,maxmin[2]-maxmin[3]);
				//h = std::max(h,maxmin[4]-maxmin[5]);
			}
		}
		
	}
	
	std::cout << "Mesh radius: " << h << std::endl;
	if( argc > 2 )
		alpha = atof(argv[2]);
	
	std::cout << "Alpha: " << alpha << std::endl;
	
	if( argc > 3 )
		left_well = atof(argv[3]);
	
	if( argc > 4 )
		right_well = atof(argv[4]);
	
	std::cout << "Well pressures: " << left_well << " " << right_well << std::endl;
	
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
					if (c->LocalID() != V_ID(i, j, k)) std::cout << "v_id = " << c->LocalID() << ", [i,j,k] = " << V_ID(i, j, k) << std::endl;
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
			if( is2d )
			{
				if( (cnt[0] > 3.0/11.0 && cnt[0] < 4.0/11.0 && cnt[1] > 5.0/11.0 && cnt[1] < 6.0/11.0) ||
				    (cnt[0] > 7.0/11.0 && cnt[0] < 8.0/11.0 && cnt[1] > 5.0/11.0 && cnt[1] < 6.0/11.0) )
					it->Delete();
			}
			else
			{
				if( (cnt[0] > 3.0/11.0 && cnt[0] < 4.0/11.0 &&
					 cnt[1] > 5.0/11.0 && cnt[1] < 6.0/11.0 &&
					 cnt[2] > 5.0/11.0 && cnt[2] < 6.0/11.0) ||
					(cnt[0] > 7.0/11.0 && cnt[0] < 8.0/11.0 &&
					 cnt[1] > 5.0/11.0 && cnt[1] < 6.0/11.0 &&
					 cnt[2] > 5.0/11.0 && cnt[2] < 6.0/11.0) )
					it->Delete();
			}
		}
		
		for (Mesh::iteratorElement it = mesh->BeginElement(FACE|EDGE|NODE); it != mesh->EndElement(); ++it)
			if( it->nbAdjElements(CELL) == 0 ) it->Delete();
	}
	
	if( alpha > 0.0 )
	{
		if( is2d )
		{
			const double eps = 1.0e-6;
			for( Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
			{
				Storage::real_array coords = it->Coords();
				if( coords[0] > eps && coords[0] < 1.0-eps && coords[1] > eps && coords[1] < 1.0-eps )
				{
					if( !(((coords[0] > 3.0/11.0-eps && coords[0] < 4.0/11.0+eps) || (coords[0] > 7.0/11.0-eps && coords[0] < 8.0/11.0+eps)) && (coords[1] > 5.0/11.0-eps && coords[1] < 6.0/11.0+eps)) )
					{
						coords[0] += alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);
						coords[1] += alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);
					}
				}
			}
		}
		else
		{
			for( Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it) if( !it->Boundary() )
			{
				it->Coords()[0] += alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);
				it->Coords()[1] += alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);
				it->Coords()[2] += alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);
			}
		}
		{
			std::stringstream str;
			str << problem_name << "_disturb_" << alpha << "_ratio_" << ratio;
			problem_name = str.str();
		}
	}
	
	mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
	
	std::cout << "nodes: " << mesh->NumberOfNodes() << " edges: " << mesh->NumberOfEdges() << " faces: " << mesh->NumberOfFaces() << " cells: " << mesh->NumberOfCells() << std::endl;
	
	
	{
		Storage::bulk_array name = mesh->self()->BulkArray(mesh->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}
	
	Tag bndcond = mesh->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE|NODE,FACE|NODE,3);
	
	Storage::integer numinner = 0, numouter = 0;
	Storage::integer numinnern = 0, numoutern = 0;
	const double eps = 1.0e-6;
	for(Mesh::iteratorElement it = mesh->BeginElement(FACE|NODE); it != mesh->EndElement(); ++it) if( it->Boundary() )
	{
		Storage::real cnt[3];
		it->Centroid(cnt);
		//if( is2d && (cnt[2] < 0.0+1.0e-9 || cnt[2] > 1.0-1.0e-9) && it->GetElementType() == FACE) continue;
		
		
		//internal part bounary
		if( (cnt[0] > 3.0/11.0-eps && cnt[0] < 4.0/11.0+eps) && (cnt[1] > 5.0/11.0-eps && cnt[1] < 6.0/11.0+eps) && (is2d || (cnt[2] > 5.0/11.0-eps && cnt[2] < 6.0/11.0+eps)) )
		{
			Storage::real_array bnd = it->RealArray(bndcond);
			bnd[0] = 1.0; //dirichlet
			bnd[1] = 0.0;
			bnd[2] = left_well;
			if( it->GetElementType() == FACE ) numinner++;
			if( it->GetElementType() == NODE ) numinnern++;
		}
		else if( (cnt[0] > 7.0/11.0-eps && cnt[0] < 8.0/11.0+eps) && (cnt[1] > 5.0/11.0-eps && cnt[1] < 6.0/11.0+eps) && (is2d || (cnt[2] > 5.0/11.0-eps && cnt[2] < 6.0/11.0+eps)) )
		{
			//std::cout << ElementTypeName(it->GetElementType()) << " " << cnt[0] << " " << cnt[1] << " " << cnt[2] << " " << (cnt[0] < eps || cnt[0] > 1.0-eps || cnt[1] < eps || cnt[1] > eps) << std::endl;
			Storage::real_array bnd = it->RealArray(bndcond);
			bnd[0] = 1.0; //dirichlet
			bnd[1] = 0.0;
			bnd[2] = right_well;
			if( it->GetElementType() == FACE )  numouter++;
			if( it->GetElementType() == NODE )  numoutern++;
		}
	}
	std::cout << "Left well faces " << numouter << " right well faces " << numinner << std::endl;
	std::cout << "Left well nodes " << numoutern << " right well nodes " << numinnern << std::endl;
	
	Storage::real mat[9];
	if( !is2d )
	{
		prepare_tensor(mat,ratio);
		std::cout << "3d tensor: " << std::endl;
	}
	else
	{
		prepare_tensor_2d(mat,ratio);
		std::cout << "2d tensor: " << std::endl;
	}
	for(int k = 0; k < 9; ++k) std::cout << mat[k] << " " ;
	std::cout << std::endl;
	
	Tag tensor = mesh->CreateTag("PERM", DATA_REAL, CELL, NONE, 9);
	for (Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
		memcpy(&it->RealArrayDF(tensor)[0],mat,sizeof(Storage::real)*9);
	
	std::cout << "I'm ready!" << std::endl;
	
	//mesh->Save("grid.vtk");
	mesh->Save("grid_out.pmf");
	//mesh->Save("grid.gmv");
	
	std::cout << "File written!" << std::endl;
	
	delete mesh;
	return 0;
}
