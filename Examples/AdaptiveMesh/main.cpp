#include "amesh.h"

using namespace INMOST;

bool output_file = false;
bool balance_mesh = true;
bool balance_mesh_refine = true;
bool balance_mesh_coarse = false;
std::string file_format = ".pmf";

int main(int argc, char ** argv)
{
	Mesh::Initialize(&argc,&argv);
	
	if( argc > 1 )
	{
		Mesh m;
		
		m.SetCommunicator(INMOST_MPI_COMM_WORLD);
		if( m.isParallelFileFormat(argv[1]) )
			m.Load(argv[1]);
		else if( m.GetProcessorRank() == 0 )
			m.Load(argv[1]);
	
		//m.SetTopologyCheck(PROHIBIT_MULTIPOLYGON | PROHIBIT_MULTILINE
		//				   | DEGENERATE_CELL | DEGENERATE_EDGE | DEGENERATE_FACE
		//				   | PRINT_NOTIFY | TRIPLE_SHARED_FACE | FLATTENED_CELL
		//				   | INTERLEAVED_FACES | NEED_TEST_CLOSURE
		//				   | DUPLICATE_EDGE | DUPLICATE_FACE | DUPLICATE_CELL
		//				   | ADJACENT_DUPLICATE | ADJACENT_HIDDEN | ADJACENT_VALID | ADJACENT_DIMENSION);
		//m.RemTopologyCheck(THROW_EXCEPTION);
#if defined(USE_PARTITIONER)
		std::vector<Storage::integer> nc(m.GetProcessorsNumber());
		Partitioner p(&m);
		if( true )
		{
			std::fill(nc.begin(),nc.end(),0); nc[m.GetProcessorRank()] = m.NumberOfCells(); m.Integrate(&nc[0],nc.size()); if( !m.GetProcessorRank() ) {std::cout << "init before "; for(unsigned q = 0; q < nc.size(); ++q) std::cout << nc[q] << " "; std::cout << std::endl;}
			//p.SetMethod(Partitioner::INNER_KMEANS,Partitioner::Partition);
			p.SetMethod(Partitioner::Parmetis,Partitioner::Partition);
			//p.SetMethod(Partitioner::Parmetis,Partitioner::Refine);
			p.Evaluate();
			m.Redistribute();
			m.ReorderEmpty(CELL|FACE|EDGE|NODE);
			
			std::fill(nc.begin(),nc.end(),0); nc[m.GetProcessorRank()] = m.NumberOfCells(); m.Integrate(&nc[0],nc.size()); if( !m.GetProcessorRank() ) {std::cout << "init after "; for(unsigned q = 0; q < nc.size(); ++q) std::cout << nc[q] << " "; std::cout << std::endl;}
			m.Barrier();
			p.SetMethod(Partitioner::Parmetis,Partitioner::Repartition);
		}
#endif
		m.ExchangeGhost(1,FACE);
		AdaptiveMesh am(m);
		//m.SetTopologyCheck(NEED_TEST_CLOSURE);
		//m.SetTopologyCheck(PROHIBIT_MULTILINE);
		//m.SetTopologyCheck(PROHIBIT_MULTIPOLYGON);
		TagInteger indicator = m.CreateTag("INDICATOR",DATA_INTEGER,CELL,NONE,1);
		TagReal wgt = m.CreateTag("WEIGHT",DATA_REAL,CELL,NONE,1);
#if defined(USE_PARTITIONER)
		p.SetWeight(wgt);
#endif
		/*
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
			indicator[*it] = 1;
		if( !m.Refine(indicator) ) std::cout << "refine failed" << std::endl;
		else std::cout << "refine ok!" << std::endl;
		m.Save("grid.pmf");
		return 0;
		*/
		int max_levels = 2;
		if( argc > 2 ) max_levels = atoi(argv[2]);
		//bounding box around mesh
		Storage::real cmax[3] = {-1.0e20,-1.0e20,-1.0e20};
		Storage::real cmin[3] = {+1.0e20,+1.0e20,+1.0e20};
		Storage::real xyz[3], r, q, cnt[3], cnt0[3], r0;
		//find bounding box around mesh
		for(Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it)
		{
			for(int d = 0; d < 3; ++d)
			{
				if( it->Coords()[d] > cmax[d] ) cmax[d] = it->Coords()[d];
				if( it->Coords()[d] < cmin[d] ) cmin[d] = it->Coords()[d];
			}
		}
		m.AggregateMax(cmax,3);
		m.AggregateMin(cmin,3);
		r0 = 1;
		for(int d = 0; d < 3; ++d)
		{
			r0 *= cmax[d]-cmin[d];
			cnt[d] = (cmax[d]+cmin[d])*0.5;
			cnt0[d] = cnt[d];
		}
		//r = pow(r,1.0/3.0)/20.0;
		r0 = pow(r0,1.0/3.0);
		r = r0/8.0;
		
		int ncells, nfaces, nedges, nnodes;

		double time = Timer();
		
		for(int k = 0; k < 80; ++k)
		{


			double step_time = Timer();

			cnt[0] = cnt0[0] + 0.25*r0*sin(k/20.0*M_PI);
			cnt[1] = cnt0[1] + 0.25*r0*cos(k/20.0*M_PI);

			m.ClearFile();
			
			//std::fill(nc.begin(),nc.end(),0); nc[m.GetProcessorRank()] = m.NumberOfCells(); m.Integrate(&nc[0],nc.size()); if( !m.GetProcessorRank() ) {std::cout << "start "; for(unsigned q = 0; q < nc.size(); ++q) std::cout << nc[q] << " "; std::cout << std::endl;}
			
			
			int numref;
			int refcnt = 0;
			do
			{
				numref = 0;
				for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) 
				{
					if( am.GetLevel(it->self()) < max_levels )
					{
						it->Barycenter(xyz);
						//refine a circle
						q = 0;
						for(int d = 0; d < 2; ++d)
							q += (xyz[d]-cnt[d])*(xyz[d]-cnt[d]);
						q = sqrt(q);
						//if( q < r*(k+1) && q > r*k)
						if( q < r*1.4 && q > r)
						{
							indicator[*it] = 1;
							numref++;
						}
					}
					else indicator[*it] = 0;
				}
				numref = m.Integrate(numref);
				if( numref )
				{
#if defined(USE_PARTITIONER)
					if( balance_mesh_refine && refcnt == 0)
					{
						//m.Barrier();
						am.ComputeWeightRefine(indicator,wgt);
						p.SetWeight(wgt);
						std::fill(nc.begin(),nc.end(),0); nc[m.GetProcessorRank()] = m.NumberOfCells(); m.Integrate(&nc[0],nc.size()); if( !m.GetProcessorRank() ) {std::cout << "refine before "; for(unsigned q = 0; q < nc.size(); ++q) std::cout << nc[q] << " "; std::cout << std::endl;}
						
						p.Evaluate();
						m.Redistribute();
						
						std::fill(nc.begin(),nc.end(),0); nc[m.GetProcessorRank()] = m.NumberOfCells(); m.Integrate(&nc[0],nc.size()); if( !m.GetProcessorRank() ) {std::cout << "refine after "; for(unsigned q = 0; q < nc.size(); ++q) std::cout << nc[q] << " "; std::cout << std::endl;}
						//m.Barrier();
					}
#endif
					ncells = m.TotalNumberOf(CELL);
					nfaces = m.TotalNumberOf(FACE);
					nedges = m.TotalNumberOf(EDGE);
					nnodes = m.TotalNumberOf(NODE);
					if( m.GetProcessorRank() == 0 )
						std::cout << "beg k " << k << " refcnt " << refcnt << " cells " << ncells << " faces " << nfaces << " nedges " << nedges << " nnodes " << nnodes << std::endl;
					//m.BeginSequentialCode();
					//std::cout << m.GetProcessorRank() << " cells " << m.NumberOfCells() << std::endl;
					//m.EndSequentialCode();

					if (!am.Refine(indicator)) break;
					
					
					ncells = m.TotalNumberOf(CELL);
					nfaces = m.TotalNumberOf(FACE);
					nedges = m.TotalNumberOf(EDGE);
					nnodes = m.TotalNumberOf(NODE);
					if( m.GetProcessorRank() == 0 )
						std::cout << "end k " << k << " refcnt " << refcnt << " cells " << ncells << " faces " << nfaces << " nedges " << nedges << " nnodes " << nnodes << std::endl;
					
					if( false )
					{
						std::stringstream file;
						file << "ref_" << k << "_" << refcnt << file_format;
						m.Save(file.str());
						if( m.GetProcessorRank() == 0 )
							std::cout << "Save " << file.str() << std::endl;
					}
					
				}
				refcnt++;
			}
			while(numref);
			refcnt = 0;
			do
			{
				numref = 0;
				for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
				{
					if( am.GetLevel(it->self()) > 0 )
					{
						it->Barycenter(xyz);
						//refine a circle
						q = 0;
						for(int d = 0; d < 2; ++d)
							q += (xyz[d]-cnt[d])*(xyz[d]-cnt[d]);
						q = sqrt(q);
						//if( q < r*k)
						if( !(q < r*1.4 && q > r) )
						{
							indicator[*it] = 1;
							numref++;
						}
					}
					else indicator[*it] = 0;
				}
				numref = m.Integrate(numref);
				if( numref )
				{
#if defined(USE_PARTITIONER)
					if( balance_mesh_coarse && refcnt == 0)
					{
						m.Barrier();
						am.ComputeWeightCoarse(indicator,wgt);
						p.SetWeight(wgt);
						std::fill(nc.begin(),nc.end(),0); nc[m.GetProcessorRank()] = m.NumberOfCells(); m.Integrate(&nc[0],nc.size()); if( !m.GetProcessorRank() ) {std::cout << "coarse before "; for(unsigned q = 0; q < nc.size(); ++q) std::cout << nc[q] << " "; std::cout << std::endl;}
						p.Evaluate();
						m.Redistribute();
						std::fill(nc.begin(),nc.end(),0); nc[m.GetProcessorRank()] = m.NumberOfCells(); m.Integrate(&nc[0],nc.size()); if( !m.GetProcessorRank() ) {std::cout << "coarse after "; for(unsigned q = 0; q < nc.size(); ++q) std::cout << nc[q] << " "; std::cout << std::endl;}
						m.Barrier();
					}
#endif
					
					ncells = m.TotalNumberOf(CELL);
					nfaces = m.TotalNumberOf(FACE);
					nedges = m.TotalNumberOf(EDGE);
					nnodes = m.TotalNumberOf(NODE);
					if( m.GetProcessorRank() == 0 )
						std::cout << ":beg k " << k << " crscnt " << refcnt << " cells " << ncells << " faces " << nfaces << " nedges " << nedges << " nnodes " << nnodes << std::endl;
					//m.BeginSequentialCode();
					//std::cout << m.GetProcessorRank() << " cells " << m.NumberOfCells() << std::endl;
					//m.EndSequentialCode();
					
					if( !am.Coarse(indicator) ) break;
					
					ncells = m.TotalNumberOf(CELL);
					nfaces = m.TotalNumberOf(FACE);
					nedges = m.TotalNumberOf(EDGE);
					nnodes = m.TotalNumberOf(NODE);
					if( m.GetProcessorRank() == 0 )
						std::cout << ":end k " << k << " crscnt " << refcnt << " cells " << ncells << " faces " << nfaces << " nedges " << nedges << " nnodes " << nnodes << std::endl;
					
					if( false ) 
					{
						std::stringstream file;
						file << "crs_" << k << "_" << refcnt << file_format;
						m.Save(file.str());
						if( m.GetProcessorRank() == 0 )
							std::cout << "Save " << file.str() << std::endl;
					}
					
				}
				
				refcnt++;
			}
			while(numref);
			
			
#if defined(USE_PARTITIONER)
			if( balance_mesh )
			{
				//m.Barrier();
				p.SetWeight(Tag());
				//std::fill(nc.begin(),nc.end(),0); nc[m.GetProcessorRank()] = m.NumberOfCells(); m.Integrate(&nc[0],nc.size()); if( !m.GetProcessorRank() ) {std::cout << "finish before "; for(unsigned q = 0; q < nc.size(); ++q) std::cout << nc[q] << " "; std::cout << std::endl;}
				p.Evaluate();
				m.Redistribute();
				//std::fill(nc.begin(),nc.end(),0); nc[m.GetProcessorRank()] = m.NumberOfCells(); m.Integrate(&nc[0],nc.size()); if( !m.GetProcessorRank() ) {std::cout << "finish after "; for(unsigned q = 0; q < nc.size(); ++q) std::cout << nc[q] << " "; std::cout << std::endl;}
				//m.Barrier();
			}
#endif
			
			std::fill(nc.begin(),nc.end(),0); nc[m.GetProcessorRank()] = m.NumberOfCells(); m.Integrate(&nc[0],nc.size()); if( !m.GetProcessorRank() ) {std::cout << "finish "; for(unsigned q = 0; q < nc.size(); ++q) std::cout << nc[q] << " "; std::cout << std::endl;}
			 
			 

			
			if( output_file )
			{
				TagInteger tag_owner = m.CreateTag("OWN",DATA_INTEGER,CELL,NONE,1);
				TagInteger tag_owner0 = m.GetTag("OWNER_PROCESSOR");
				TagInteger tag_stat = m.CreateTag("STAT",DATA_INTEGER,CELL,NONE,1);
				for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
				{
					tag_owner[*it] = tag_owner0[*it];
					tag_stat[*it] = it->GetStatus();
				}
				std::stringstream file;
				file << "step_" << k << file_format;
				m.Save(file.str());
				if( m.GetProcessorRank() == 0 )
					std::cout << "Save " << file.str() << std::endl;
			}
			else if( m.GetProcessorRank() == 0 ) std::cout << "step " << k << std::endl;
			
			step_time = Timer() - step_time;
			if( m.GetProcessorRank() == 0 ) std::cout << "step time " << step_time << std::endl;
		}


		time = Timer() - time;
		if( m.GetProcessorRank() == 0 ) std::cout << "total time: " << time << " processors " << m.GetProcessorsNumber() << std::endl;
	}
	else std::cout << "Usage: " << argv[0] << " mesh_file [max_levels=2]" << std::endl;
	
	Mesh::Finalize();
	return 0;
}
