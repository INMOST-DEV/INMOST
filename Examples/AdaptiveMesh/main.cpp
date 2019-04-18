#include "amesh.h"

using namespace INMOST;

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
	
		//m.SetTopologyCheck(PROHIBIT_MULTIPOLYGON | PROHIBIT_MULTILINE | DEGENERATE_CELL | DEGENERATE_EDGE | DEGENERATE_FACE | PRINT_NOTIFY | TRIPLE_SHARED_FACE | FLATTENED_CELL | INTERLEAVED_FACES | NEED_TEST_CLOSURE);
		//m.RemTopologyCheck(THROW_EXCEPTION);
#if defined(USE_PARTITIONER)
		Partitioner p(&m);
		if( true )
		{
			std::cout << "before on " << m.GetProcessorRank() << " " << m.NumberOfCells() << std::endl;
			p.SetMethod(Partitioner::INNER_KMEANS,Partitioner::Partition);
			p.Evaluate();
			m.Redistribute();
			m.ReorderEmpty(CELL|FACE|EDGE|NODE);
			std::cout << "after on " << m.GetProcessorRank() << " " << m.NumberOfCells() << std::endl;
		}
#endif
		m.ExchangeGhost(2,FACE);
		AdaptiveMesh am(m);
		//m.SetTopologyCheck(NEED_TEST_CLOSURE);
		//m.SetTopologyCheck(PROHIBIT_MULTILINE);
		//m.SetTopologyCheck(PROHIBIT_MULTIPOLYGON);
		TagInteger indicator = m.CreateTag("INDICATOR",DATA_INTEGER,CELL,NONE,1);
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
		Storage::real xyz[3], r, q, cnt[3];
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
		r = 1;
		for(int d = 0; d < 3; ++d)
		{
			r *= cmax[d]-cmin[d];
			cnt[d] = (cmax[d]+cmin[d])*0.5;
		}
		r = pow(r,1.0/3.0)/20.0;
		
		for(int k = 0; k < 15; ++k)
		{

			m.ClearFile();
			
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
						for(int d = 0; d < 3; ++d)
							q += (xyz[d]-cnt[d])*(xyz[d]-cnt[d]);
						q = sqrt(q);
						if( q < r*(k+1) && q > r*k)
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
					int ncells = m.TotalNumberOf(CELL);
					if( m.GetProcessorRank() == 0 )
						std::cout << "k " << k << " refcnt " << refcnt << " " <<  r*k << " < r < " << r*(k+1) << " cells " << ncells << std::endl;
					//m.BeginSequentialCode();
					//std::cout << m.GetProcessorRank() << " cells " << m.NumberOfCells() << std::endl;
					//m.EndSequentialCode();

					if (!am.Refine(indicator)) break;
					
					if( false )
					{
						std::stringstream file;
						file << "ref_" << k << "_" << refcnt << ".pvtk";
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
						for(int d = 0; d < 3; ++d)
							q += (xyz[d]-cnt[d])*(xyz[d]-cnt[d]);
						q = sqrt(q);
						if( q < r*k)
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
					int ncells = m.TotalNumberOf(CELL);
					if( m.GetProcessorRank() == 0 )
						std::cout << ": k " << k << " crscnt " << refcnt << " " << r*k << " < r < " << r*(k+1) << " cells " << ncells <<  std::endl;
					//m.BeginSequentialCode();
					//std::cout << m.GetProcessorRank() << " cells " << m.NumberOfCells() << std::endl;
					//m.EndSequentialCode();
					
					if( !am.Coarse(indicator) ) break;
					
					if( false ) 
					{
						std::stringstream file;
						file << "crs_" << k << "_" << refcnt << ".pvtk";
						m.Save(file.str());
						if( m.GetProcessorRank() == 0 )
							std::cout << "Save " << file.str() << std::endl;
					}
					
				}
				
				refcnt++;
			}
			while(numref);
			
			
#if defined(USE_PARTITIONER)
			if( false )
			{
				p.Evaluate();
				am.Repartition();
			}
#endif
			 
			 

			
			if( true )
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
				file << "step_" << k << ".pvtk";
				m.Save(file.str());
				if( m.GetProcessorRank() == 0 )
					std::cout << "Save " << file.str() << std::endl;
			}
			else if( m.GetProcessorRank() == 0 ) std::cout << "step " << k << std::endl;
		}
	}
	else std::cout << "Usage: " << argv[0] << " mesh_file [max_levels=2]" << std::endl;
	
	Mesh::Finalize();
	return 0;
}
