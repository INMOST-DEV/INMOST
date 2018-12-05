#include "amesh.h"

using namespace INMOST;

int main(int argc, char ** argv)
{
	Mesh::Initialize(&argc,&argv);
	
	if( argc > 1 )
	{
		Mesh m;
		m.SetCommunicator(INMOST_MPI_COMM_WORLD);
		m.Load(argv[1]);
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
		
		for(int k = 0; k < 16; ++k)
		{
			
			int numref;
			int refcnt = 0;
			do
			{
				numref = 0;
				for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) if( am.GetLevel(it->self()) < max_levels )
				{
					it->Barycenter(xyz);
					//refine a circle
					q = 0;
					for(int d = 0; d < 3; ++d)
						q += (xyz[d]-cnt[d])*(xyz[d]-cnt[d]);
					q = sqrt(q);
					if( q < r*(k+1) && q > r*k)
					//if( q < 0.02 )
					{
						indicator[it->self()] = 1;
						numref++;
					}
				}
				numref = m.Integrate(numref);
				if( numref )
				{
					std::cout << "k " << k << " refcnt " << refcnt << " " <<  r*k << " < r < " << r*(k+1) << " numref " << numref << " cells " << m.NumberOfCells() << std::endl;
					/*
					{
						std::stringstream file;
						file << "indicator_" << k << "_" << refcnt << "r.pmf";
						m.Save(file.str());
					}
					*/
					if( !am.Refine(indicator) ) break;
					
					{
						std::stringstream file;
						file << "dump_" << k << "_" << refcnt << "r.pmf";
						m.Save(file.str());
					}
					
					//cleanup indicator
					for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) indicator[it->self()] = 0;
				}
				refcnt++;
			}
			while(numref);
			refcnt = 0;
			do
			{
				numref = 0;
				for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) if( am.GetLevel(it->self()) > 0 )
				{
					it->Barycenter(xyz);
					//refine a circle
					q = 0;
					for(int d = 0; d < 3; ++d)
						q += (xyz[d]-cnt[d])*(xyz[d]-cnt[d]);
					q = sqrt(q);
					if( q < r*k)
					{
						indicator[it->self()] = 1;
						numref++;
					}
				}
				numref = m.Integrate(numref);
				if( numref )
				{
					std::cout << "k " << k << " crscnt " << refcnt << " " << r*k << " < r < " << r*(k+1) << " numcrs " << numref << " cells " << m.NumberOfCells() <<  std::endl;
					/*
					{
						std::stringstream file;
						file << "indicator_" << k << "_" << refcnt << "c.pmf";
						m.Save(file.str());
					}
					 */
					if( !am.Coarse(indicator) ) break;
					
					{
						std::stringstream file;
						file << "dump_" << k << "_" << refcnt << "c.pmf";
						m.Save(file.str());
					}
					
					
					
					//cleanup indicator
					for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) indicator[it->self()] = 0;
				}
				//return -1;
				refcnt++;
			}
			while(numref);
			
			{
				std::stringstream file;
				file << "step_" << k << ".pvtk";
				m.Save(file.str());
			}
			
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
				file << "step_" << k << ".pmf";
				m.Save(file.str());
			}
		}
	}
	else std::cout << "Usage: " << argv[0] << " mesh_file [max_levels=2]" << std::endl;
	
	Mesh::Finalize();
	return 0;
}
