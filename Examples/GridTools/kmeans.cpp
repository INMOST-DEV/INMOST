#include "inmost.h"

using namespace INMOST;

//todo: want to separate all faces into those that share edges
//see todo: inside text




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
	int total_points = 0, total_global_points = 0;
	
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:total_points)
#endif
	for(int q = 0; q < m.CellLastLocalID(); ++q) if( m.isValidCell(q) )
	{
		Cell n = m.CellByLocalID(q);
		if( n->GetStatus() != Element::Ghost ) total_points++;
	}
	
#if defined(USE_MPI)
	MPI_Allreduce(&total_global_points,&total,points,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD)
#else
	total_global_points = total_points;
#endif
	
	
	std::vector< double > points_center(total_points*3);
	std::vector< int > points_node(total_points);
	std::vector< int > points_cluster(total_points,-1);
	std::vector< double > cluster_center(K*3);
	std::vector< int > cluster_npoints(K);
#if defined(USE_MPI)
	std::vector< double > cluster_center_tmp(K*3);
	std::vector< int > cluster_npoints_tmp(K);
#endif
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
	
	std::cout << m.GetProcessorRank() << " init clusters" << std::endl;
	
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
			std::vector<int> recvcounts(m.GetProcessorsNumber()), displs(m.GetProcessorsNumber());
			displs[0] = 0;
			recvcounts[0] = Kpart*3;
			for(int i = 1; i < (int)m.GetProcessorsNumber(); ++i)
			{
				displs[i] = displs[i-1] + recvcounts[i-1];
				recvcounts[i] = Kpart*3;
			}
			recvcounts.back() = K*3 - displs.back();
			std::cout << m.GetProcessorRank() << " " << Kstart << " " << Kend << std::endl;
			for(int i = 0; i < (int)m.GetProcessorsNumber(); ++i)
			{
				std::cout << m.GetProcessorRank() << " " << i << " " << displs[i] << " " << recvcounts[i] << std::endl;
			}
			for(int i = Kstart; i < Kend; ++i)
			{
				cluster_center_tmp[i*3+0] = cluster_center[i*3+0];
				cluster_center_tmp[i*3+1] = cluster_center[i*3+1];
				cluster_center_tmp[i*3+2] = cluster_center[i*3+2];
			}
			MPI_Allgatherv(&cluster_center_tmp[Kstart*3],(Kend-Kstart)*3,MPI_DOUBLE,&cluster_center[0],&recvcounts[0],&displs[0],MPI_DOUBLE,MPI_COMM_WORLD);
		}
#endif
	}
	
	std::cout << m.GetProcessorRank() << " start clustering" << std::endl;
	
	int iter = 1;
	double t = Timer();
	while(true)
	{
		int changed = 0;
		
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
			
			if(id_old_cluster != id_nearest_center)
			{
				points_cluster[i] = id_nearest_center;
				changed++;
			}
		}
		
#if defined(USE_MPI)
		int tmp = changed;
		MPI_Allreduce(&tmp,&changed,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
		
		if(changed == 0 || iter >= max_iterations)
		{
			std::cout << m.GetProcessorRank() << " break in iteration " << iter << std::endl;
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
		
		std::cout << "Iteration " << iter << std::endl;
		iter++;
	}
	std::cout << "Clustering in " << Timer() - t << " secs " << std::endl;
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
