#include "inmost.h"



#if defined(USE_PARTITIONER)

#if defined(USE_PARTITIONER_ZOLTAN)
#include <zoltan.h>
#endif
#if defined(USE_PARTITIONER_PARMETIS)
#include <metis.h>
#include <parmetis.h>
#endif
#include <sstream>
#include <deque>
#define ZOLTAN_CHKERR(x) if( x != ZOLTAN_OK )  {throw ErrorInPartitioner;}


#if defined(USE_PARALLEL_WRITE_TIME)
#define REPORT_MPI(x) {m->WriteTab(m->GetStream()) << "<MPI><![CDATA[" << #x << "]]></MPI>" << std::endl; x;}
#define REPORT_STR(x) {m->WriteTab(m->GetStream()) << "<TEXT><![CDATA[" << x << "]]></TEXT>" << std::endl;}
#define REPORT_VAL(str,x) {m->WriteTab(m->GetStream()) << "<VALUE name=\"" << str << "\"> <CONTENT><![CDATA[" << x << "]]> </CONTENT><CODE><![CDATA[" << #x << "]]></CODE></VALUE>" << std::endl;}
#define ENTER_FUNC() long double all_time = Timer(); m->WriteTab(m->GetStream()) << "<FUNCTION name=\"" << __FUNCTION__ << "\" id=\"func" << m->GetFuncID()++ << "\">" << std::endl; m->Enter();
#define EXIT_FUNC() m->WriteTab(m->GetStream()) << "<TIME>" << Timer() - all_time << "</TIME>" << std::endl; m->Exit(); m->WriteTab(m->GetStream()) << "</FUNCTION>" << std::endl;
#else
#define REPORT_MPI(x) x
#define REPORT_STR(x) {}
#define REPORT_VAL(str,x) {}
#define ENTER_FUNC() {}
#define EXIT_FUNC() {}
#endif

namespace INMOST
{
#if defined(USE_PARTITIONER_ZOLTAN)
	static int get_number_of_objects(void *data, int *ierr)
	{
		Partitioner * p = static_cast<Partitioner *>(data);
		Mesh * m = p->GetMesh();
		*ierr = ZOLTAN_OK;
		int num = 0;
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++)
			if( it->GetStatus() != Element::Ghost ) num++;
		return num;
	}
	
	static void get_object_list(void *data, int sizeGID, int sizeLID,
								ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
								int wgt_dim, float *obj_wgts, int *ierr)
	{
		Partitioner * p = static_cast<Partitioner *>(data);
		Mesh * m = p->GetMesh();
		*ierr = ZOLTAN_OK;
		int i = 0;	 
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++)
			if( it->GetStatus() != Element::Ghost )
			{
				globalID[i] = it->GlobalID();
				localID[i] = i;
				for(int j = 0; j < wgt_dim; j++)
					 obj_wgts[i] = static_cast<float>(it->RealArray(p->GetWeight())[j]);
				i++;
			}
	}
	
	static int get_num_geometry(void *data, int *ierr)
	{
		Partitioner * p = static_cast<Partitioner *>(data);
		Mesh * m = p->GetMesh();
		*ierr = ZOLTAN_OK;
		return m->GetDimensions();
	}
	
	static void get_geometry_list(void *data, int sizeGID, int sizeLID,
								  int num_obj,
								  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
								  int num_dim, double *geom_vec, int *ierr)
	{
		Partitioner * p = static_cast<Partitioner *>(data);
		Mesh * m = p->GetMesh();
		int dim = m->GetDimensions();
		if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != dim))
		{
			*ierr = ZOLTAN_FATAL;
			return;
		}
		
		*ierr = ZOLTAN_OK;
		
		
		int i = 0;
		std::vector<Storage::real> coords(dim);
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++)
			if( it->GetStatus() != Element::Ghost )
			{
				it->Centroid(&coords[0]);
				for(int j = 0; j < dim; j++)
					geom_vec[dim*i + j] = coords[j];
				i++;
			}
		return;
	}
	
	static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
								   int num_obj,
								   ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
								   int *numEdges, int *ierr)
	{
		Partitioner * p = static_cast<Partitioner *>(data);
		Mesh * m = p->GetMesh();
		
		if ( (sizeGID != 1) || (sizeLID != 1))
		{
			*ierr = ZOLTAN_FATAL;
			return;
		}
		
		int i = 0;	
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++)
			if( it->GetStatus() != Element::Ghost )
			{
				numEdges[i] = 0;
				ElementArray<Face> faces = it->getFaces();
				for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); jt++)
				{
					Cell n = it->Neighbour(jt->self());
					if( n.isValid() ) numEdges[i]++;
				}		
				i++;
			}
		
		*ierr = ZOLTAN_OK;
		return;
	}
	
	static void get_edge_list(void *data, int sizeGID, int sizeLID,
							  int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
							  int *num_edges,
							  ZOLTAN_ID_PTR nborGID, int *nborProc,
							  int wgt_dim, float *ewgts, int *ierr)
	{
		int *nextProc;
		float * nextWgt;
		ZOLTAN_ID_TYPE *nextNbor;
		
		Partitioner * p = static_cast<Partitioner *>(data);
		Mesh * m = p->GetMesh();
		*ierr = ZOLTAN_OK;
		
		if ( (sizeGID != 1) || (sizeLID != 1))
		{
			*ierr = ZOLTAN_FATAL;
			return;
		}
		
		nextNbor = nborGID;
		nextProc = nborProc;
		nextWgt = ewgts;
		
		int i = 0;
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++)
			if( it->GetStatus() != Element::Ghost )
			{
				ElementArray<Face> faces = it->getFaces();
				for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); jt++)
				{
					Cell n = it->Neighbour(jt->self());
					if( n.isValid() )
					{
						*nextNbor++ = n->GlobalID();
						*nextProc++ = n->Integer(m->OwnerTag());
						for(int j = 0; j < wgt_dim; j++)
							*nextWgt++ = static_cast<float>(jt->RealArray(p->GetWeight())[j]);
					}
				}		
				i++;
			}
		return;
	}
#endif
	
	Partitioner::Partitioner(Mesh * _m)
	{
		m = _m;
		if( m->GetMeshState() == Mesh::Serial ) m->ResolveShared();
		pzz = NULL;
	}
	
	Partitioner::Partitioner(const Partitioner & other)
	{
#if defined(USE_PARTITIONER_ZOLTAN)
		if( pt == Zoltan_Parmetis || pt == Zoltan_Scotch || pt == Zoltan_PHG || pt == Zoltan_RCB || pt == Zoltan_RIB || pt == Zoltan_HSFC )
		{
			if( other.pzz != NULL ) pzz = static_cast<void *>(Zoltan_Copy(static_cast<struct Zoltan_Struct *>(other.pzz)));
			else pzz = NULL;
		}
#endif
		m = other.m;
		weight_tag = other.weight_tag;
	}
	
	Partitioner & Partitioner::operator = (Partitioner const & other)
	{
#if defined(USE_PARTITIONER_ZOLTAN)
		if( pt == Zoltan_Parmetis || pt == Zoltan_Scotch || pt == Zoltan_PHG || pt == Zoltan_RCB || pt == Zoltan_RIB || pt == Zoltan_HSFC )
		{
			if( other.pzz != NULL ) 
			{
				if( pzz == NULL ) pzz = static_cast<void *>(Zoltan_Copy(static_cast<struct Zoltan_Struct *>(other.pzz)));
				else Zoltan_Copy_To(static_cast<struct Zoltan_Struct *>(pzz),static_cast<struct Zoltan_Struct *>(other.pzz));
			}
			else pzz = NULL;
		}
#endif
		m = other.m;
		weight_tag = other.weight_tag;
		return *this;
	}
	
	Partitioner::~Partitioner()
	{
#if defined(USE_PARTITIONER_ZOLTAN)
		if( pt == Zoltan_Parmetis || pt == Zoltan_Scotch || pt == Zoltan_PHG || pt == Zoltan_RCB || pt == Zoltan_RIB || pt == Zoltan_HSFC )
		{
			if( pzz != NULL )
			{
				struct Zoltan_Struct * zz = static_cast<struct Zoltan_Struct *>(pzz);
				Zoltan_Destroy(&zz);
				pzz = static_cast<void *>(zz);
			}
		}
#endif
	}
	
	void Partitioner::Initialize(int * argc, char *** argv)
	{
		(void) argc;
		(void) argv;
#if defined(USE_PARTITIONER_ZOLTAN)
		float ver;
		Zoltan_Initialize(*argc,*argv,&ver);
#endif
	}
	
	void Partitioner::Finalize()
	{
	}
	
	
	void Partitioner::Evaluate()
	{
		ENTER_FUNC();
		int package = 0;
		switch(pt)
		{
			case Inner_RCM:
				package = 0;
				break;
			case Zoltan_Parmetis:
			case Zoltan_PHG:
			case Zoltan_RIB:
			case Zoltan_RCB:
			case Zoltan_HSFC:
			case Zoltan_Scotch:
				package = 1;
				break;
			case Parmetis:
				package = 2;
				break;
		}
		if( package == 0 )
		{
#if defined(USE_MPI)
			if( m->GetMeshState() == Mesh::Serial ) m->ResolveShared();
			
			unsigned mpisize = m->GetProcessorsNumber();
			
//~ #if defined(NDEBUG)
			//Check!!
			unsigned ncell = m->NumberOfCells();
			std::vector<unsigned> ncells(mpisize);
			REPORT_MPI(MPI_Allgather(&ncell,1,MPI_UNSIGNED,&ncells[0],1,MPI_UNSIGNED,m->GetCommunicator()));
			unsigned nhave = 0;
			for(unsigned k = 0; k < mpisize; k++) if( ncells[k] ) nhave++;
			if( nhave > 1 )
			{
				std::cout << "currently only non-distributed meshes are supported by reverse cuthill-mckee algorithm" << std::endl;
				throw NotImplemented;
			}
//~ #endf
			
			
			Tag index = m->CreateTag("CM_INDEX",DATA_INTEGER,CELL,NONE,1);
			Tag order = m->CreateTag("ORDER",DATA_INTEGER,CELL,NONE,1);
			
			
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			{
				ElementArray<Face> f = it->getFaces();
				for(ElementArray<Face>::iterator jt = f.begin(); jt != f.end(); ++jt)
					if( !jt->Boundary() ) it->Integer(order)++;
			}
			
			Storage::integer visited = 1;
			Storage::integer chunk = m->NumberOfCells() / mpisize + 1;
			std::deque<HandleType> queue;
			ElementArray<Element> around(m);
			HandleType select;
			while(visited != static_cast<Storage::integer>(m->NumberOfCells()+1) )
			{
				select = InvalidHandle();
				Mesh::iteratorCell halt_it;
				for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				{
					if( it->Integer(index) == 0 )
					{
						select = *it;
						halt_it = it;
						break;
					}
				}
				if( select == InvalidHandle() ) break;
				Storage::integer sorder = m->Integer(select,order);
				for(Mesh::iteratorCell it = halt_it; it != m->EndCell(); ++it)
				{
					if( it->Integer(index) == 0 )
					{
						Storage::integer norder = it->Integer(order);
						if( norder < sorder ) 
						{
							select = *it;
							sorder = norder;
						}
					}
				}
				queue.push_back(select);
				m->Integer(queue.back(),index) = visited++;
				while(!queue.empty())
				{
					around = Cell(m,queue.front())->BridgeAdjacencies(FACE,CELL);
					queue.pop_front();
					std::sort(around.data(),around.data()+around.size(),Mesh::IntegerComparator(m,order));
					for(ElementArray<Element>::iterator it = around.begin(); it != around.end(); ++it)
					{
						if( it->Integer(index) == 0 )
						{
							queue.push_back(*it);
							m->Integer(queue.back(),index) = visited++;

							if( visited % chunk == 0 ) 
							{
								queue.clear();
								break;
							}

						}
					}
					//this may not fire properly
					//if( visited % chunk == 0 ) queue.clear();
				}
			}
			m->DeleteTag(order);
			//index may be used
			
			//reverse
			//for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			//	it->Integer(index) = m->NumberOfCells() - it->Integer(index);
				
			
			
			Tag redist = m->RedistributeTag();
			
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			{
				it->Integer(redist) = it->Integer(index)/chunk;
				//~ std::cout << it->Integer(index) << " -> " <<  it->Integer(redist) << std::endl;
			}

			m->DeleteTag(index);
#else
			throw NotImplemented;
#endif
		}
		if( package == 1 )
		{
#if defined(USE_PARTITIONER_ZOLTAN)
			int changes, numGidEntries, numLidEntries, numImport, numExport;
			ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
			int *importProcs, *importToPart, *exportProcs, *exportToPart;
			int rc;
			if( m->Integer(m->GetHandle(),m->LayersTag()) == 0 ) m->ExchangeGhost(1,FACE);
			
			rc = Zoltan_LB_Partition(static_cast<struct Zoltan_Struct *>(pzz), /* input (all remaining fields are output) */
									 &changes,        /* 1 if partitioning was changed, 0 otherwise */ 
									 &numGidEntries,  /* Number of integers used for a global ID */
									 &numLidEntries,  /* Number of integers used for a local ID */
									 &numImport,      /* Number of vertices to be sent to me */
									 &importGlobalGids,  /* Global IDs of vertices to be sent to me */
									 &importLocalGids,   /* Local IDs of vertices to be sent to me */
									 &importProcs,    /* Process rank for source of each incoming vertex */
									 &importToPart,   /* New partition for each incoming vertex */
									 &numExport,      /* Number of vertices I must send to other processes*/
									 &exportGlobalGids,  /* Global IDs of the vertices I must send */
									 &exportLocalGids,   /* Local IDs of the vertices I must send */
									 &exportProcs,    /* Process to which I send each of the vertices */
									 &exportToPart);  /* Partition to which each vertex will belong */
			ZOLTAN_CHKERR(rc);		
			
			if( changes )
			{
				std::vector<int> new_distribution;
				for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++)
					if( it->GetStatus() != Element::Ghost )
						new_distribution.push_back(it->Integer(m->OwnerTag()));	
				for(int i = 0; i < numExport; i++)
					new_distribution[exportLocalGids[i]] = exportToPart[i];
				
				Tag redist = m->RedistributeTag();
				int i = 0;
				for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++)
					if( it->GetStatus() != Element::Ghost )
						it->Integer(redist) = new_distribution[i++];
			}
			Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
								&importProcs, &importToPart);
			Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
								&exportProcs, &exportToPart);
#endif
		}
		if( package == 2 )
		{
#if defined(USE_PARTITIONER_PARMETIS)
			INMOST_DATA_ENUM_TYPE dim = m->GetDimensions();
			double time;
			bool debug_output = false;
			bool time_output = false;
			int result;
			int rank = m->GetProcessorRank();
			int size = m->GetProcessorsNumber();
			int have_vwgt = 0, have_adjwgt = 0, have_mwgt = 0;
			int nreserve = 32;
			idx_t nparts;
			//Array storing parts indexes
			std::vector<unsigned int> xparts;
			//These arrays determine local sparse matrix
			std::vector<idx_t> vtxdist;
			std::vector<idx_t> xadj;
			std::vector<idx_t> adjncy;
			//These arrays determine weights defined for vertices
			std::vector<idx_t> vwgt; //unused yet
			std::vector<idx_t> adjwgt; //unused yet
			idx_t wgtflag = 0;
			idx_t numflag = 0;
			idx_t ndims = m->GetDimensions();
			std::vector<real_t> xyz;
			//This sets weight for every part set in respect to according weight and tolerance
			idx_t ncon = 1;
			std::vector<real_t> tpwgts;
			std::vector<real_t> ubvec;
			//set options
			std::vector<idx_t> options(3,0);
			//number of edgecuts
			idx_t edgecut;
			//output data - which entity should belong to which part
			std::vector<idx_t> part;
			real_t itr = 1000.0;
			//send data to processors with zero nodes
			std::vector<unsigned int> send;
			idx_t mysize;
			//to estimate memory used by entities
			std::vector<idx_t> vsize; //unused yet
			
			//compute part ranges for every processor
			//and total number of parts
			xparts.resize(size+1);
			xparts[0] = 0;
			
			for(unsigned int i = 1; i < xparts.size(); i++)
				xparts[i] = i;
				
			
			nparts = xparts[xparts.size()-1];
			
			
			//set graph size
			vtxdist.resize(size+1);
			
			


			//set options
			if( debug_output )
			{	
				options[0] = 1;
				options[1] = 1;
				options[2] = 15;
			}
			
			
			if( time_output ) time = Timer();
			//set ghost boundaries so we certanly know we can access adjacencies
			
			if( m->GetMeshState() != Mesh::Parallel )  
				m->ResolveShared();
			if( m->Integer(m->GetHandle(),m->LayersTag()) == 0 ) 
				m->ExchangeGhost(1,FACE);
			//~ if( !m->GlobalIDTag().isValid() || !m->GlobalIDTag().isDefined(CELL) )
				m->AssignGlobalID(CELL);
			

			
			
			if( time_output ) time = Timer();
			
			
			mysize = 0;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) if( it->GetStatus() != Element::Ghost ) mysize++;
			vtxdist[0] = 0;
			REPORT_MPI(result = MPI_Allgather(&mysize,1,IDX_T,&vtxdist[1],1,IDX_T,m->GetCommunicator()));
			if( result != MPI_SUCCESS ) throw Impossible;
			for(unsigned int i = 2; i < vtxdist.size(); i++)
				vtxdist[i] += vtxdist[i-1];
				
				
			if( time_output )
			{
				time = Timer()-time;
				REPORT_STR("time to compute graph information");
				REPORT_VAL("time",time);
			}


			{	
				xadj.resize(mysize+1);
				adjncy.reserve(mysize*nreserve); //preallocate array
				xadj[0] = 0;
			}
			
			xyz.resize(mysize*dim);
			part.resize(mysize);
									
			
			if( time_output ) time = Timer();
			
			//compute itr parameter
			{
				itr = 1000; //communication time to distribution time
			}


			
			
			//estimate storage space
			if( pa == Repartition ) 
			{
				vsize.resize(mysize,1);
				int k = 0;
				for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) if( it->GetStatus() != Element::Ghost )
				{
					vsize[k++] = m->MemoryUsage(*it);
				}
			}
			
			
			
			
			{
				if( GetWeight().isValid() && GetWeight().GetSize() != ENUMUNDEF )
				{
					if( GetWeight().isDefined(MESH) ) have_mwgt = 1;
					if( GetWeight().isDefined(CELL) && !GetWeight().isSparse(CELL) ) have_vwgt = 1;
					if( GetWeight().isDefined(FACE) && !GetWeight().isSparse(FACE) ) have_adjwgt = 1;
					if( have_vwgt + have_adjwgt ) ncon = GetWeight().GetSize();
				}
				
				
				if( have_vwgt ) vwgt.resize(mysize*ncon);
				if( have_adjwgt) adjwgt.reserve(mysize*ncon*nreserve);
				if( !have_adjwgt && !have_vwgt ) ncon = 1;
				//set constraint tolerances
				const real_t v = static_cast<real_t>(1.05);
				ubvec.resize(ncon,v);
				if( have_mwgt )
				{
					std::vector<real_t> l_tpwgts(ncon);
					for(idx_t q = 0; q < ncon; q++) l_tpwgts[q] = static_cast<real_t>(m->RealArray(m->GetHandle(),GetWeight())[q]);
					REPORT_MPI(MPI_Allgather(&l_tpwgts[0],ncon,IDX_T,&tpwgts[0],ncon,IDX_T,m->GetCommunicator()));
					
					for(idx_t j = 0; j < ncon; j++)
					{
						real_t sum_tpwgts = 0;
						for(unsigned i = 0; i < tpwgts.size()/ncon; i++) sum_tpwgts += tpwgts[i*ncon+j];
						for(unsigned i = 0; i < tpwgts.size()/ncon; i++) tpwgts[i*ncon+j] /= sum_tpwgts;
					}
				}
				else
				{
					real_t c = static_cast<real_t>(1.0/static_cast<real_t>(nparts));
					tpwgts.resize(nparts*ncon,c);
				}
			}			
			
			if( time_output ) time = Timer();
			
			//get part tag
			
			if( pa == Repartition )
			{
				part.resize(mysize,rank);
			}
			
			if( time_output )
			{
				time = Timer()-time;
				REPORT_STR("time to compute weights");
				REPORT_VAL("time",time);
			}
			
			REPORT_STR("start making graph, number of local elements:");
			REPORT_VAL("local_size",mysize);
			
			
			if( time_output ) time = Timer();
			
			
			//make graph
			int k = 0;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			if( it->GetStatus() != Element::Ghost )
			{
				idx_t sum = 0;
				if( pa == Partition )
				{
					Storage::real xyz_temp[3] = {0.0,0.0,0.0};
					it->Centroid(xyz_temp);
					for(int q = 0; q < m->GetDimensions(); ++q)
						xyz[k*dim+q] = static_cast<real_t>(xyz_temp[q]);
				}
				{
					ElementArray<Face> faces = it->getFaces();
					for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
					{
						Cell n = it->Neighbour(jt->self());
						if( n.isValid() )
						{
							adjncy.push_back(n->GlobalID());
							if( have_adjwgt ) 
							{
								for(idx_t q = 0; q < ncon; q++) 
									adjwgt.push_back(jt->RealArray(GetWeight())[q]); //why adjwgt is not real???
							}
							sum++;
						}
					}
					
					xadj[k+1] = xadj[k] + sum;
					
					//Fill weight arrays here
					if( have_vwgt )
					{
						for(idx_t q = 0; q < ncon; q++)
							vwgt[k*ncon+q] = it->RealArrayDF(GetWeight())[q]; //why vwgt is not real???
					}
				}
				k++;
			}
			
			REPORT_STR("graph made, number of adjacencies: ");
			REPORT_VAL("adjacencies_size",adjncy.size());
			
			
			if( time_output )
			{
				time = Timer()-time;
				REPORT_STR("time to compute graph");
				REPORT_VAL("time",time);
			}
			
			
			wgtflag = have_vwgt + have_adjwgt*2;
			
			if( debug_output )
			{
				m->BeginSequentialCode();
				std::cout << "on processor " << rank << std::endl;
				std::cout << "vtxdist: " << std::endl;
				for(unsigned int i = 0; i < vtxdist.size(); i++)
				{
					std::cout << vtxdist[i] << " ";
					if( (i+1)%30 == 0 ) std::cout << std::endl;
				}
				std::cout << std::endl;
				
				{
					std::cout << "xadj size " << xadj.size() << ":" << std::endl;
					for(unsigned int i = 0; i < xadj.size(); i++)
					{
						std::cout << xadj[i] << " ";
						if( (i+1)%30 == 0 ) std::cout << std::endl;
					}
					std::cout << std::endl;
					std::cout << "adjncy size " << adjncy.size() << ":" << std::endl;
					for(unsigned int i = 0; i < adjncy.size(); i++)
					{
						std::cout << adjncy[i] << " ";
						if( (i+1)%30 == 0 ) std::cout << std::endl;
					}
					std::cout << std::endl;
					if( have_vwgt )
					{
						std::cout << "vwgt: " << std::endl;
						for(unsigned int i = 0; i < vwgt.size(); i++)
						{
							std::cout << vwgt[i] << " ";
							if( (i+1)%30 == 0 ) std::cout << std::endl;
						}
						std::cout << std::endl;
					}
					if( have_adjwgt )
					{
						std::cout << "adjwgt: " << std::endl;
						for(unsigned int i = 0; i < adjwgt.size(); i++)
						{
							std::cout << adjwgt[i] << " ";
							if( (i+1)%30 == 0 ) std::cout << std::endl;
						}
						std::cout << std::endl;
					}
					if( pa == Repartition )
					{
						std::cout << "vsize: " << std::endl;
						for(unsigned int i = 0; i < vsize.size(); i++)
						{
							std::cout << vsize[i] << " ";
							if( (i+1)%30 == 0 ) std::cout << std::endl;
						}
						std::cout << std::endl;
					}			
					std::cout << "wgtflag = " << wgtflag << std::endl;
					std::cout << "numflag = " << numflag << std::endl;
				}
				if( pa == Partition )
				{
					std::cout << "ndims = " << ndims << std::endl;
					std::cout << "xyz: " << std::endl;
					for(unsigned int i = 0; i < xyz.size(); i++)
					{
						std::cout << xyz[i] << " ";
						if( (i+1)%30 == 0 ) std::cout << std::endl;
					}
					std::cout << std::endl;
				}
				std::cout << "ncon = " << ncon << std::endl;
				std::cout << "nparts = " << nparts << std::endl;
				std::cout << "tpwgts: " << std::endl;
				for(unsigned int i = 0; i < tpwgts.size(); i++)
				{
					std::cout << tpwgts[i] << " ";
					if( (i+1)%30 == 0 ) std::cout << std::endl;
				}
				std::cout << std::endl;
				std::cout << "ubvec: " << std::endl;
				for(unsigned int i = 0; i < ubvec.size(); i++)
				{
					std::cout << ubvec[i] << " ";
					if( (i+1)%30 == 0 ) std::cout << std::endl;
				}
				std::cout << std::endl;
				std::cout << "options: " << std::endl;
				for(unsigned int i = 0; i < options.size(); i++)
				{
					std::cout << options[i] << " ";
					if( (i+1)%30 == 0 ) std::cout << std::endl;
				}
				std::cout << std::endl;
				m->EndSequentialCode();
			}
			
			//redistribute entities;
			if( time_output ) time = Timer();
			
			
			
			
			{
				int it,jt,kt,flag;
				unsigned int fill[3];
				
				send.reserve(3*nparts*2);
				for(it = size; it >= 1; it--)
				{
					if( vtxdist[it]-vtxdist[it-1] == 0 )
					{
						fill[0] = it - 1; //to
						flag = 0;
						for(jt = it-1; jt >= 1; jt--)
						{
							if( vtxdist[jt]-vtxdist[jt-1] > 0 )
							{
								fill[1] = jt - 1; //from
								fill[2] = vtxdist[jt]-1-vtxdist[rank]; //id
								send.insert(send.end(),fill,fill+3);
								vtxdist[jt]--;
								for(kt = jt+1; kt < it; kt++)
									vtxdist[kt] = vtxdist[jt];
								flag = 1;
								break;
							}
						}
						if( !flag ) throw Impossible;
					}
				}
				
				REPORT_STR("redistribute graph:");
				REPORT_VAL("send_size",send.size()/3);
				
				
				for(unsigned int i = 0; i < send.size()/3; i++)
				{
					idx_t trans_adjncy_size;
					real_t trans_xyz[3];
					std::vector<idx_t> trans_adjncy;
					if( send[i*3+0] == send[i*3+1] ) throw Impossible;
					if( send[i*3+1] == static_cast<unsigned int>(rank) )
					{	
						{
							trans_adjncy_size = xadj[send[i*3+2]+1]-xadj[send[i*3+2]];
							for(jt = xadj[send[i*3+2]]; jt < xadj[send[i*3+2]+1]; jt++)
								trans_adjncy.push_back(adjncy[jt]);
							if( have_adjwgt )
								for(jt = xadj[send[i*3+2]]; jt < xadj[send[i*3+2]+1]; jt++)
									trans_adjncy.push_back(adjwgt[jt]);
							if( have_vwgt )
								trans_adjncy.push_back(vwgt[send[i*3+2]]);
							if( pa == Repartition )
								trans_adjncy.push_back(part[send[i*3+2]]);
							REPORT_VAL("to",send[i*3+0]);
							REPORT_VAL("send_size",(trans_adjncy_size*(1+have_adjwgt)));
							REPORT_MPI(result = MPI_Send(&trans_adjncy[0],trans_adjncy_size*(1+have_adjwgt),IDX_T,send[i*3+0],2,m->GetCommunicator()));
							if( result != MPI_SUCCESS ) 
								throw Impossible;
						}
						if( pa == Partition )
						{
							for(jt = 0; jt < static_cast<int>(dim); jt++)
								trans_xyz[jt] = xyz[send[i*3+2]*dim+jt];
							REPORT_VAL("to",send[i*3+0]);
							REPORT_MPI(result = MPI_Send(trans_xyz,dim,REAL_T,send[i*3+0],3,m->GetCommunicator()));
							if( result != MPI_SUCCESS ) throw Impossible;
						}
					}
					else if( send[i*3+0] == static_cast<unsigned int>(rank) )
					{
						MPI_Status stat;
						int msgsize;
						
						{
							REPORT_VAL("from",send[i*3+1]);
							REPORT_MPI(result = MPI_Probe(send[i*3+1],2,m->GetCommunicator(),&stat));
							if( result != MPI_SUCCESS ) throw Impossible;
							result = MPI_Get_count(&stat,IDX_T,&msgsize);
							if( result != MPI_SUCCESS ) throw Impossible;
							trans_adjncy_size = (msgsize-have_vwgt-((pa == Repartition)?1:0))/(1+have_adjwgt);
							trans_adjncy.resize(msgsize);
							REPORT_VAL("from",send[i*3+1]);
							REPORT_MPI(result = MPI_Recv(&trans_adjncy[0],msgsize,IDX_T,send[i*3+1],2,m->GetCommunicator(),&stat));
							if( result != MPI_SUCCESS ) throw Impossible;
							xadj.resize(2);
							xadj[0] = 0;
							xadj[1] = trans_adjncy_size;
							for(jt = 0; jt < trans_adjncy_size; jt++)
								adjncy.push_back(trans_adjncy[jt]);
							if( have_adjwgt )
								for(jt = 0; jt < trans_adjncy_size; jt++)
									adjwgt.push_back(trans_adjncy[jt+trans_adjncy_size]);
							if( have_vwgt )
								vwgt.push_back(trans_adjncy[trans_adjncy_size*2]);
							if( pa == Repartition )
								part.push_back(trans_adjncy[trans_adjncy_size*2+1]);
						}
						if( pa == Partition )
						{
							REPORT_VAL("from",send[i*3+1]);
							REPORT_MPI(result = MPI_Recv(&trans_xyz,dim,REAL_T,send[i*3+1],3,m->GetCommunicator(),&stat));
							if( result != MPI_SUCCESS ) throw Impossible;
							xyz.resize(dim);
							for(jt = 0; jt < static_cast<int>(dim); jt++)
								xyz[jt] = trans_xyz[jt];
						}
						part.resize(1);
					}
				}
				
				REPORT_STR("graph redistributed");
			}
			
			if( time_output )
			{
				time = Timer()-time;
				REPORT_STR("time to redistribute graph");
				REPORT_VAL("time",time);
			}
			
			REPORT_STR("BEGIN PARMETIS");
			
			if( time_output ) time = Timer();
			
			MPI_Comm comm = m->GetCommunicator();
			
			idx_t  * link_vtxdist = vtxdist.empty() ? NULL : &vtxdist[0];
			idx_t  * link_xadj    = xadj.empty()    ? NULL : &xadj[0];
			idx_t  * link_adjncy  = adjncy.empty()  ? NULL : &adjncy[0];
			idx_t  * link_vwgt    = vwgt.empty()    ? NULL : &vwgt[0];
			idx_t  * link_adjwgt  = adjwgt.empty()  ? NULL : &adjwgt[0];
			real_t * link_tpwgts  = tpwgts.empty()  ? NULL : &tpwgts[0];
			real_t * link_ubvec   = ubvec.empty()   ? NULL : &ubvec[0];
			idx_t  * link_options = options.empty() ? NULL : &options[0];
			idx_t  * link_part    = part.empty()    ? NULL : &part[0];
			real_t * link_xyz     = xyz.empty()     ? NULL : &xyz[0];
			idx_t  * link_vsize   = vsize.empty()   ? NULL : &vsize[0];
			
			REPORT_VAL("vtxdist",link_vtxdist);
			REPORT_VAL("xadj",link_xadj);
			REPORT_VAL("adjncy",link_adjncy);
			REPORT_VAL("vwgt",link_vwgt);
			REPORT_VAL("adjwgt",link_adjwgt);
			REPORT_VAL("tpwgts",link_tpwgts);
			REPORT_VAL("ubvec",link_ubvec);
			REPORT_VAL("options",link_options);
			REPORT_VAL("part",link_part);
			REPORT_VAL("xyz",link_xyz);
			REPORT_VAL("vsize",link_vsize);

			switch(pa)
			{
				//~ case Partition:
				//~ result = ParMETIS_V3_PartKway(&vtxdist[0],&xadj[0],&adjncy[0],&vwgt[0],&adjwgt[0],
											  //~ &wgtflag,&numflag,&ncon,&nparts,&tpwgts[0],
											  //~ &ubvec[0],&options[0],&edgecut,&part[0],&comm);
				//~ break;
				case Refine:
				REPORT_STR("run refine");
				result = ParMETIS_V3_RefineKway(link_vtxdist,link_xadj,link_adjncy,link_vwgt,link_adjwgt,
												&wgtflag,&numflag,&ncon,&nparts,link_tpwgts,
												link_ubvec,link_options,&edgecut,link_part,&comm);	
				break;
				case Partition:
				REPORT_STR("run partition");
				result = ParMETIS_V3_PartGeomKway(link_vtxdist,link_xadj,link_adjncy,link_vwgt,link_adjwgt,
												  &wgtflag,&numflag,&ndims,link_xyz,&ncon,&nparts,link_tpwgts,
												  link_ubvec,link_options,&edgecut,link_part,&comm);
				break;
				//~ case P_PartGeom:
				//~ result = ParMETIS_V3_PartGeom(&vtxdist[0],&ndims,&xyz[0],&part[0],&comm);
				//~ break;
				case Repartition:
				REPORT_STR("run repartition");
				result = ParMETIS_V3_AdaptiveRepart(link_vtxdist,link_xadj,link_adjncy,link_vwgt,link_vsize,
													link_adjwgt,&wgtflag,&numflag,&ncon,&nparts,link_tpwgts,
													link_ubvec,&itr,link_options,&edgecut,link_part,&comm);
				break;
			}
			
			
			
			if( time_output )
			{
				time = Timer()-time;
				REPORT_STR("time in parmetis");
				REPORT_VAL("time",time);
			}
			
			REPORT_STR("END PARMETIS");

			if( time_output ) time = Timer();

			//distribute computed parts back
			{
				REPORT_VAL("send_size",send.size()/3);
				for(unsigned int i = 0; i < send.size()/3; i++)
				{
					if( send[i*3+0] == static_cast<unsigned int>(rank) )
					{
						REPORT_VAL("send_part",part[0]);
						REPORT_MPI(result = MPI_Send(&part[0],1,IDX_T,send[i*3+1],4,m->GetCommunicator()));
						if( result != MPI_SUCCESS ) throw Impossible;
						part.clear();
					}
					else if( send[i*3+1] == static_cast<unsigned int>(rank) )
					{ 
						REPORT_MPI(result = MPI_Recv(&part[send[i*3+2]],1,IDX_T,send[i*3+0],4,m->GetCommunicator(),MPI_STATUS_IGNORE));
						REPORT_VAL("recv_part",part[send[i*3+2]]);
						REPORT_VAL("from",part[send[i*3+0]]);
						if( result != MPI_SUCCESS ) throw Impossible;
					}
				}
			}

			if( time_output )
			{
				time = Timer()-time;
				REPORT_STR("redistribute graph back");
				REPORT_VAL("time",time);
			}
									 
			//debug
			if( debug_output )
			{
				m->BeginSequentialCode();
				std::cout << "distribution on " << rank << std::endl;
				for(unsigned int i = 0; i < part.size(); i++)
				{
					std::cout << part[i] << " ";
					if( (i+1)%30 == 0 ) std::cout << std::endl;
				}
				std::cout << std::endl;
				m->EndSequentialCode();
			}

			//assign entities to part sets
			{
				
				Tag redist = m->RedistributeTag();
				
				int k = 0;
				for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) 
					if( it->GetStatus() != Element::Ghost ) 
						it->Integer(redist) = part[k++];
				
				
			}
#endif
		}
		EXIT_FUNC();
	}
	
	void Partitioner::SetMethod(enum Type t, enum Action a)
	{
		pt = t;
		pa = a;
		int package = 0;
		switch(pt)
		{
			case Inner_RCM:
				package = 0;
				break;
			case Zoltan_Parmetis:
			case Zoltan_PHG:
			case Zoltan_RIB:
			case Zoltan_RCB:
			case Zoltan_HSFC:
			case Zoltan_Scotch:
				package = 1;
				break;
			case Parmetis:
				package = 2;
				break;
		}
		if( package == 1 )
		{
#if defined(USE_PARTITIONER_ZOLTAN)
			if(pzz == NULL)
			{
				struct Zoltan_Struct * zz = Zoltan_Create(m->GetCommunicator());
				Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, this);
				Zoltan_Set_Obj_List_Fn(zz, get_object_list, this);
				Zoltan_Set_Num_Geom_Fn(zz, get_num_geometry, this);
				Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list, this);
				Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, this);
				Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, this);
				
				Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");		
				Zoltan_Set_Param(zz, "LB_METHOD", "RIB");
				Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
				Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
				Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
				Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");
				Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
				pzz = static_cast<void *>(zz);
			}
			struct Zoltan_Struct * zz = static_cast<struct Zoltan_Struct *>(pzz);
			switch(t)
			{
				case Zoltan_Parmetis:
				{
					Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
					Zoltan_Set_Param(zz,"GRAPH_PACKAGE","PARMETIS");
					Zoltan_Set_Param(zz,"PARMETIS_METHOD","AdaptiveRepart");
					Zoltan_Set_Param(zz,"SCATTER_GRAPH","2");
					break;
				}
				case Zoltan_Scotch:
				{
					Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
					Zoltan_Set_Param(zz,"GRAPH_PACKAGE","SCOTCH");
					Zoltan_Set_Param(zz,"SCATTER_GRAPH","2");
					break;
				}
				case Zoltan_RCB:
				{
					Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
					break;
				}
				case Zoltan_RIB:
				{
					Zoltan_Set_Param(zz, "LB_METHOD", "RIB");
					break;
				}
				case Zoltan_HSFC:
				{
					Zoltan_Set_Param(zz, "LB_METHOD", "HSFC");
					
					break;
				}
				case Zoltan_PHG:
				{
					Zoltan_Set_Param(zz,"LB_METHOD", "GRAPH");
					Zoltan_Set_Param(zz,"GRAPH_PACKAGE","PHG");
					Zoltan_Set_Param(zz,"SCATTER_GRAPH","0");
					break;
				}
			}
			switch(a)
			{
				case Partition:
					Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
					break;
				case Repartition:
					Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");
					break;
				case Refine:	
					Zoltan_Set_Param(zz, "LB_APPROACH", "REFINE");
					break;
			}
#else
			throw NotImplemented;
#endif
		}
		if( package == 2 )
		{
#if defined(USE_PARTITIONER_PARMETIS)
#endif
		}
	}
	
	void Partitioner::SetWeight(Tag weight)
	{
		(void) weight;
		int package = 0;
		switch(pt)
		{
			case Inner_RCM:
				package = 0;
				break;
			case Zoltan_Parmetis:
			case Zoltan_PHG:
			case Zoltan_RIB:
			case Zoltan_RCB:
			case Zoltan_HSFC:
			case Zoltan_Scotch:
				package = 1;
				break;
			case Parmetis:
				package = 2;
				break;
		}
		if( package == 1 )
		{
#if defined(USE_PARTITIONER_ZOLTAN)
			struct Zoltan_Struct * zz = static_cast<struct Zoltan_Struct *>(pzz);
			std::stringstream Wobj, Wedge;
			weight_tag = weight;
			if( weight.isValid() )
			{
				if( weight.GetSize() == ENUMUNDEF ) throw UnknownWeightSize;
				if( weight.isDefined(CELL) )
					Wobj << weight.GetSize();
				else
					Wobj << 0;
				if( weight.isDefined(FACE) )
					Wedge << weight.GetSize();
				else
					Wedge << 0;
			}
			else
			{
				Wobj << 0;
				Wedge << 0;
			}
			Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", Wobj.str().c_str());
			Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", Wedge.str().c_str());
#else
			throw NotImplemented;
#endif
		}
		if( package == 2 )
		{
#if defined(USE_PARTITIONER_PARMETIS)
			weight_tag = weight;
#endif
		}
	}
	void Partitioner::ResetWeight()
	{
		int package = 0;
		switch(pt)
		{
			case Inner_RCM:
				package = 0;
				break;
			case Zoltan_Parmetis:
			case Zoltan_PHG:
			case Zoltan_RIB:
			case Zoltan_RCB:
			case Zoltan_HSFC:
			case Zoltan_Scotch:
				package = 1;
				break;
			case Parmetis:
				package = 2;
				break;
		}
		if( package == 1 )
		{
#if defined(USE_PARTITIONER_ZOLTAN)
			struct Zoltan_Struct * zz = static_cast<struct Zoltan_Struct *>(pzz);
			Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
			Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");
#else
			throw NotImplemented;
#endif
		}
		if( package == 2 )
		{
#if defined(USE_PARTITIONER_PARMETIS)
#endif
		}
		weight_tag == Tag();
	}
	Mesh * Partitioner::GetMesh()
	{
		return m;
	}
	Tag Partitioner::GetWeight()
	{
		return weight_tag;
	}
}
#endif
