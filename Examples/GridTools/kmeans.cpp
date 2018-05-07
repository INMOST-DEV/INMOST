#include "inmost.h"

using namespace INMOST;

//todo: want to separate all faces into those that share edges
//see todo: inside text

struct kdtree
{
	struct entry
	{
		double v[3];
		int pos;
	};
	struct compare
	{
		int comp;
		compare(int c) : comp(c) {}
		compare(const compare & b) : comp(b.comp) {}
		bool operator () (const entry & a, const entry & b) { return a.v[comp] < b.v[comp]; }
	};
	entry * set;
	kdtree * l, * r;
	int size;
	double split;
	kdtree()
	{
		size = 0;
		set = NULL;
		l = NULL;
		r = NULL;
	}
	~kdtree()
	{
		if( l != NULL ) delete l;
		if( r != NULL ) delete r;
		set = NULL;
		size = 0;
	}
	void resize(int set_size)
	{
		clear();
		size = set_size;
		set = new entry[size];
	}
	void clear()
	{
		if( set != NULL ) delete [] set;
		size = 0;
	}
	void build_tree(int comp)
	{
		if( l != NULL ) {delete l; l = NULL;}
		if( r != NULL ) {delete r; r = NULL;}
		if( size > 1 )
		{
			std::sort(set,set+size,compare(comp));
			int middle = size/2;
			l = new kdtree;
			l->set = set;
			l->size = middle;
			r = new kdtree;
			r->set = set + middle;
			r->size = size - middle;
			split = set[middle].v[comp];
			l->build_tree((comp+1)%3);
			r->build_tree((comp+1)%3);
		}
	}
	void build(double * coords)
	{
		for(int k = 0; k < size; ++k)
		{
			set[k].v[0] = coords[k*3+0];
			set[k].v[1] = coords[k*3+1];
			set[k].v[2] = coords[k*3+2];
			set[k].pos = k;
		}
		build_tree(0);
	}
	void closest_tree(double p[3], int & pos, double & min_dist, int comp)
	{
		if( size == 0 )
		{
			std::cout << "size is " << size << std::endl;
		}
		else if ( size == 1 )
		{
			double v[3];
			v[0] = (set[0].v[0] - p[0]);
			v[1] = (set[0].v[1] - p[1]);
			v[2] = (set[0].v[2] - p[2]);
			double dist = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
			if (dist < min_dist)
			{
				min_dist = dist;
				pos = set[0].pos;
			}
		}
		else
		{
			if (p[comp] < split)
			{
				// search left first
				l->closest_tree(p, pos, min_dist, (comp+1)%3);
				if (p[comp] + min_dist >= split)
					r->closest_tree(p, pos, min_dist, (comp+1)%3);
			}
			else
			{
				// search right first
				r->closest_tree(p, pos, min_dist, (comp+1)%3);
				if (p[comp] - min_dist <= split)
					l->closest_tree(p, pos, min_dist, (comp+1)%3);
			}
		}
	}
	int closest(double p[3])
	{
		double min_dist = std::numeric_limits<double>::max();
		int ret;
		closest_tree(p,ret,min_dist,0);
		return ret;
	}
};

int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " mesh clusters [clusters=10] [max_iterations=10] [mesh_out=grid.vtk]" << std::endl;
		return -1;
	}
	
	Mesh::Initialize(&argc,&argv);
	
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
	m.SetCommunicator(INMOST_MPI_COMM_WORLD);
	m.Load(argv[1]);

	std::cout << "Start:" << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;

	double tt = Timer();
	//There seems to be a problem with the mesh
	/*
	int fixed = 0;
	for (Mesh::iteratorFace f = m.BeginFace(); f != m.EndFace(); ++f)
	{
		if (f->FixNormalOrientation())
			fixed++;
	}
	std::cout << "Time to fix normals: " << Timer() - tt << std::endl;
	std::cout << "Total face normals fixed: " << fixed << std::endl;
	*/
	int total_points = 0;
	
	
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:total_points)
#endif
	for(int q = 0; q < m.CellLastLocalID(); ++q) if( m.isValidCell(q) )
	{
		Cell n = m.CellByLocalID(q);
		if( n->GetStatus() != Element::Ghost ) total_points++;
	}
	std::vector< int > points_node(total_points);
	std::vector< double > points_center(total_points*3);
	std::vector< int > points_cluster(total_points,-1);
	
	int k = 0;
	for(int q = 0; q < m.CellLastLocalID(); ++q) if( m.isValidCell(q) )
	{
		Cell n = m.CellByLocalID(q);
		if( n->GetStatus() != Element::Ghost ) points_node[k++] = n->LocalID();
	}

#if defined(USE_OMP)
#pragma omp paralell for
#endif
	for(int q = 0; q < total_points; ++q)
	{
		Cell n = m.CellByLocalID(points_node[q]);
		double cnt[3];
		n->Centroid(cnt);
		points_center[q*3+0] = cnt[0];
		points_center[q*3+1] = cnt[1];
		points_center[q*3+2] = cnt[2];
	}
	
#if defined(USE_MPI)
	std::vector<int> displs(m.GetProcessorsNumber());
	std::vector<int> counts(m.GetProcessorsNumber());
	std::vector<int> npoints(m.GetProcessorsNumber());
	int total_global_points = 0;
	int total_local_points = 0;
	bool balance = false;
	MPI_Allgather(&total_points,1,MPI_INT,&npoints[0],1,MPI_INT,MPI_COMM_WORLD);
	double imbalance = 1;
	for(int k = 0; k < (int) m.GetProcessorsNumber(); ++k)
	{
		for(int l = k+1; l < (int) m.GetProcessorsNumber(); ++l)
		{
			imbalance = std::max(imbalance,npoints[k]/(double)npoints[l]);
			imbalance = std::max(imbalance,npoints[l]/(double)npoints[k]);
		}
	}
	if( m.GetProcessorRank() == 0 )
		std::cout << "Imbalance is " << imbalance << std::endl;
	if( imbalance > 1.2 )
		balance = true;
	if( balance )//redistribute points
	{
		for(int k = 0; k < (int) m.GetProcessorsNumber(); ++k)
			total_global_points += npoints[k];
		total_local_points = (int)ceil((double)total_global_points/(double)m.GetProcessorsNumber());
		
		std::vector<double> points_center_global(m.GetProcessorRank() == 0 ? total_global_points*3 : 1);
		displs[0] = 0;
		counts[0] = npoints[0]*3;
		for(int k = 1; k < (int) m.GetProcessorsNumber(); ++k)
		{
			displs[k] = displs[k-1] + counts[k-1];
			counts[k] = npoints[k]*3;
		}
		MPI_Gatherv(&points_center[0],total_points*3,MPI_DOUBLE,&points_center_global[0],&counts[0],&displs[0],MPI_DOUBLE,0,MPI_COMM_WORLD);
		displs[0] = 0;
		counts[0] = total_local_points*3;
		for(int k = 1; k < (int) m.GetProcessorsNumber(); ++k)
		{
			displs[k] = displs[k-1] + counts[k-1];
			counts[k] = total_local_points*3;
		}
		counts.back() = total_global_points*3 - displs.back();
		total_points = counts[m.GetProcessorRank()]/3;
		points_center.resize(total_points*3);
		MPI_Scatterv(&points_center_global[0],&counts[0],&displs[0],MPI_DOUBLE,&points_center[0],total_points*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
		points_cluster.resize(total_points,-1);
	}
#endif
	
	
	
	
	std::vector< double > cluster_center(K*3);
	std::vector< int > cluster_npoints(K);
#if defined(USE_MPI)
	std::vector< double > cluster_center_tmp(K*3);
	std::vector< int > cluster_npoints_tmp(K);
#endif
	//kdtree tree;
	//tree.resize(K);
	
	
	if( m.GetProcessorRank() == 0 )
		std::cout << "Init clusters" << std::endl;
	
	// choose K distinct values for the centers of the clusters
	{
		std::vector<int> prohibited_indexes;
		int Kpart = (int)ceil((double)K/(double)m.GetProcessorsNumber());
		int Kstart = Kpart * m.GetProcessorRank();
		int Kend = Kstart + Kpart;
		if( Kend > K ) Kend = K;
		for(int i = Kstart; i < Kend; i++)
		{
			while(true)
			{
				int index_point = rand() % total_points;
				
				if(find(prohibited_indexes.begin(), prohibited_indexes.end(),
						index_point) == prohibited_indexes.end())
				{
					prohibited_indexes.push_back(index_point);
					cluster_center[i*3+0] = points_center[index_point*3+0];
					cluster_center[i*3+1] = points_center[index_point*3+1];
					cluster_center[i*3+2] = points_center[index_point*3+2];
					points_cluster[index_point] = i;
					break;
				}
			}
		}
		//guter cluster centers over all processes
#if defined(USE_MPI)
		{
			//std::vector<int> counts(m.GetProcessorsNumber()), displs(m.GetProcessorsNumber());
			displs[0] = 0;
			counts[0] = Kpart*3;
			for(int i = 1; i < (int)m.GetProcessorsNumber(); ++i)
			{
				displs[i] = displs[i-1] + counts[i-1];
				counts[i] = Kpart*3;
			}
			counts.back() = K*3 - displs.back();
			for(int i = Kstart; i < Kend; ++i)
			{
				cluster_center_tmp[i*3+0] = cluster_center[i*3+0];
				cluster_center_tmp[i*3+1] = cluster_center[i*3+1];
				cluster_center_tmp[i*3+2] = cluster_center[i*3+2];
			}
			MPI_Allgatherv(&cluster_center_tmp[Kstart*3],(Kend-Kstart)*3,MPI_DOUBLE,&cluster_center[0],&counts[0],&displs[0],MPI_DOUBLE,MPI_COMM_WORLD);
		}
#endif
	}
	
	if( m.GetProcessorRank() == 0 )
		std::cout << "Start clustering" << std::endl;
	
	int iter = 1;
	double t = Timer();
	while(true)
	{
		
		
		int changed = 0;
		
		//tree.build(&cluster_center[0]);
		
		// associates each point to the nearest center
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:changed)
#endif
		for(int i = 0; i < total_points; i++)
		{
			int id_old_cluster = points_cluster[i];
			int id_nearest_center = -1;
			
			double lmin = 1.0e+100;
			
			for(int j = 0; j < K; ++j)
			{
				double v[3];
				v[0] = (points_center[i*3+0] - cluster_center[j*3+0]);
				v[1] = (points_center[i*3+1] - cluster_center[j*3+1]);
				v[2] = (points_center[i*3+2] - cluster_center[j*3+2]);
				
				double l = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				if( l < lmin )
				{
					lmin = l;
					id_nearest_center = j;
				}
			}
			
			//id_nearest_center = tree.closest(&points_center[i*3]);
			
			if(id_old_cluster != id_nearest_center)
			{
				points_cluster[i] = id_nearest_center;
				changed++;
			}
		}
		
#if defined(USE_MPI)
		int tmp = changed;
		MPI_Allreduce(&tmp,&changed,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
		
		if(changed == 0 || iter >= max_iterations)
		{
			if( m.GetProcessorRank() == 0 )
				std::cout << "Break in iteration " << iter << std::endl;
			break;
		}
		
		for(int i = 0; i < K; i++)
		{
			cluster_center[i*3+0] = 0;
			cluster_center[i*3+1] = 0;
			cluster_center[i*3+2] = 0;
			cluster_npoints[i] = 0;
		}
		// recalculating the center of each cluster
#if defined(USE_OMP)
#pragma omp parallel
#endif
		{
			std::vector< double > local_sum(K*3,0.0);
			std::vector< int > local_npoints(K,0);
			for(int j = 0; j < total_points; ++j)
			{
				local_sum[points_cluster[j]*3+0] += points_center[j*3+0];
				local_sum[points_cluster[j]*3+1] += points_center[j*3+1];
				local_sum[points_cluster[j]*3+2] += points_center[j*3+2];
				local_npoints[points_cluster[j]]++;
			}
#if defined(USE_OMP)
#pragma omp critical
#endif
			{
				for(int i = 0; i < K; ++i)
				{
					cluster_center[i*3+0] += local_sum[i*3+0];
					cluster_center[i*3+1] += local_sum[i*3+1];
					cluster_center[i*3+2] += local_sum[i*3+2];
					cluster_npoints[i] += local_npoints[i];
				}
			}
		}
#if defined(USE_MPI)
		MPI_Allreduce(&cluster_center[0],&cluster_center_tmp[0],K*3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(&cluster_npoints[0],&cluster_npoints_tmp[0],K,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		cluster_center.swap(cluster_center_tmp);
		cluster_npoints.swap(cluster_npoints_tmp);
#endif
		for(int i = 0; i < K; i++)
		{
			cluster_center[i*3+0] /= (double) cluster_npoints[i];
			cluster_center[i*3+1] /= (double) cluster_npoints[i];
			cluster_center[i*3+2] /= (double) cluster_npoints[i];
		}
		
		if( m.GetProcessorRank() == 0 )
			std::cout << "Iteration " << iter << " changed " << changed << std::endl;
		iter++;
	}
	//tree.clear();
	std::cout << "Clustering in " << Timer() - t << " secs " << std::endl;
#if defined(USE_MPI)
	if( balance )
	{
		std::vector<int> points_cluster_global(total_global_points);
		displs[0] = 0;
		counts[0] = total_local_points;
		for(int k = 1; k < (int) m.GetProcessorsNumber(); ++k)
		{
			displs[k] = displs[k-1] + counts[k-1];
			counts[k] = total_local_points;
		}
		counts.back() = total_global_points - displs.back();
			MPI_Gatherv(&points_cluster[0],total_points,MPI_INT,&points_cluster_global[0],&counts[0],&displs[0],MPI_INT,0,MPI_COMM_WORLD);
		displs[0] = 0;
		counts[0] = npoints[0];
		for(int k = 1; k < (int) m.GetProcessorsNumber(); ++k)
		{
			displs[k] = displs[k-1] + counts[k-1];
			counts[k] = npoints[k];
		}
		total_points = counts[m.GetProcessorRank()];
		points_cluster.resize(total_points);
		MPI_Scatterv(&points_cluster_global[0],&counts[0],&displs[0],MPI_INT,&points_cluster[0],total_points,MPI_INT,0,MPI_COMM_WORLD);
	}
#endif
	// shows elements of clusters
	TagInteger mat = m.CreateTag("MATERIAL",DATA_INTEGER,CELL,NONE,1);
#if defined(USE_OMP)
#pragma omp parallel for
#endif
	for(int j = 0; j < total_points; ++j)
		mat[m.CellByLocalID(points_node[j])] = points_cluster[j]+1;
	m.ExchangeData(mat,CELL,0);

	m.Save(grid_out);
	
	Mesh::Finalize();
	return 0;
}
