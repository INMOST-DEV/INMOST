#include "inmost.h"

using namespace INMOST;

//todo: want to separate all faces into those that share edges
//see todo: inside text


double dot(const double a[3], const double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}



struct coords
{
	double v[3];
	double operator [] (int i) const {return v[i];}
	double & operator [] (int i) {return v[i];}
	int compare(const double bv[3]) const
	{
		const double eps = 1.0e-11;//pow(10,-14);
		for (int i = 0; i < 3; i++)
		{
			if (fabs(v[i] - bv[i]) > eps)
			{
				if (v[i] > bv[i])
					return 1;
				else
					return -1;
			}
		}
		return 0;
	}
	bool operator <(const coords & b) const { return compare(b.v) < 0; }
	coords() {v[0] = v[1] = v[2] = 0;}
	coords(double x, double y, double z)
	{
		v[0] = x;
		v[1] = y;
		v[2] = z;
	}
	coords(double _v[3]) { memcpy(v, _v, sizeof(double)* 3); }
	coords(const coords &b) { memcpy(v, b.v, sizeof(double)* 3); }
	coords & operator =(coords const & b) { memmove(v, b.v, sizeof(double)* 3); return *this; }
	coords operator -(const coords & b) const { return coords(v[0] - b.v[0], v[1] - b.v[1], v[2] - b.v[2]); }
	coords operator +(const coords & b) const { return coords(v[0] + b.v[0], v[1] + b.v[1], v[2] + b.v[2]); }
	coords operator /(double l) const {return coords(v[0]/l, v[1]/l, v[2]/l);}
};

double len2(const coords & p)
{
	return p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
}

struct Point
{
	coords c;
	int cluster;
	int position;
	int node;
	void set_cluster(int c, int p)
	{
		cluster = c;
		position = p;
	}
	coords & get_center()
	{
		return c;
	}
};

struct Cluster
{
	std::vector<int> points;
	std::vector<int> emptyp;
	coords center;
	int add_position(int p)
	{
		if( emptyp.empty() )
		{
			points.push_back(p);
			return points.size()-1;
		}
		else
		{
			int pos = emptyp.back();
			points[pos] = p;
			emptyp.pop_back();
			return pos;
		}
	}
	void rem_position(int k)
	{
		points[k] = -1;
		emptyp.push_back(k);
	}
	void set_center(coords c)
	{
		center = c;
	}
	int num_points()
	{
		return (int)(points.size() - emptyp.size());
	}
	coords & get_center() {return center;}
};

int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " mesh clusters [clusters=10] [max_iterations=10] [mesh_out=grid.vtk]" << std::endl;
		return -1;
	}
	
	int K = 10;
	if( argc > 2 ) K = atoi(argv[2]);
	int max_iterations = 10;
	if( argc > 3 ) max_iterations = atoi(argv[3]);
	
	if( K == 0 )
		std::cout << "Problem with number of clusters argument " << argv[2] << std::endl;
	
	if( max_iterations == 0 )
		std::cout << "Problem with max iterations argument " << argv[3] << std::endl;

	std::string grid_out = "grid.vtk";
	if (argc > 4) grid_out = std::string(argv[4]);


	Mesh m;
	m.Load(argv[1]);

	std::cout << "Start:" << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;

	double tt = Timer();
	//There seems to be a problem with the mesh
	int fixed = 0;
	for (Mesh::iteratorFace f = m.BeginFace(); f != m.EndFace(); ++f)
	{
		if (f->FixNormalOrientation())
			fixed++;
	}
	std::cout << "Time to fix normals: " << Timer() - tt << std::endl;
	std::cout << "Total face normals fixed: " << fixed << std::endl;
	
	int total_points = m.NumberOfCells();
	
	std::vector< Cluster > clusters(K);
	std::vector< Point > points(total_points);
	
	int k = 0;
	double cnt[3];
	for(Mesh::iteratorCell n = m.BeginCell(); n != m.EndCell(); ++n)
	{
		n->Centroid(cnt);
		points[k].c = coords(cnt);
		points[k].cluster = -1;
		points[k].position = -1;
		points[k].node = n->LocalID();
		k++;
	}
	
	
	// choose K distinct values for the centers of the clusters
	{
		std::vector<int> prohibited_indexes;
		for(int i = 0; i < K; i++)
		{
			while(true)
			{
				int index_point = rand() % total_points;
				
				if(find(prohibited_indexes.begin(), prohibited_indexes.end(),
						index_point) == prohibited_indexes.end())
				{
					prohibited_indexes.push_back(index_point);
					int pos = clusters[i].add_position(index_point);
					points[index_point].set_cluster(i,pos);
					clusters[i].set_center(points[index_point].c);
					break;
				}
			}
		}
	}
	
	int iter = 1;
	
	while(true)
	{
		bool done = true;
		
		// associates each point to the nearest center
		for(int i = 0; i < total_points; i++)
		{
			int id_old_cluster = points[i].cluster;
			int id_nearest_center = -1;
			double lmin = 1.0e+100;
			
			for(int j = 0; j < K; ++j)
			{
				double l = len2(points[i].get_center() - clusters[j].get_center());
				if( l < lmin )
				{
					lmin = l;
					id_nearest_center = j;
				}
			}
			
			if(id_old_cluster != id_nearest_center)
			{
				if(id_old_cluster != -1)
					clusters[id_old_cluster].rem_position(points[i].position);
				
				int pos = clusters[id_nearest_center].add_position(i);
				points[i].set_cluster(id_nearest_center,pos);
				done = false;
			}
		}
		
		// recalculating the center of each cluster
		for(int i = 0; i < K; i++)
		{
			int total_points_cluster = clusters[i].num_points();
			coords sum(0,0,0);
			
			if(total_points_cluster > 0)
			{
				for(int p = 0; p < (int)clusters[i].points.size(); p++)
					if( clusters[i].points[p] != -1 )
						sum = sum + points[clusters[i].points[p]].get_center();
				clusters[i].set_center(sum / (double) total_points_cluster);
			}
		}
		
		if(done == true || iter >= max_iterations)
		{
			std::cout << "Break in iteration " << iter << "\n\n";
			break;
		}
		else std::cout << "Iteration " << iter << std::endl;
		iter++;
	}
	
	// shows elements of clusters
	TagInteger mat = m.CreateTag("MATERIAL",DATA_INTEGER,CELL,NONE,1);
	for(int i = 0; i < K; i++)
	{
		for(int j = 0; j < (int)clusters[i].points.size(); ++j)
			if( clusters[i].points[j] != -1 )
				mat[m.CellByLocalID(points[clusters[i].points[j]].node)] = i+1;
	}
	

	m.Save(grid_out);
	return 0;
}
