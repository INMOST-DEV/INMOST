#include "inmost.h"


#if defined(USE_MESH)
#include <iostream>
#include <fstream>
#include <sstream>

//#define DEBUG_REMOVE_GHOST
//#define DEBUG_COMPUTE_SHARED_SKIN_SET

#if defined(USE_MPI)
static INMOST_DATA_BIG_ENUM_TYPE pmid = 0;
#endif

#if defined(USE_PARALLEL_WRITE_TIME)
#define REPORT_MPI(x) {WriteTab(out_time) << __FUNCTION__ << " call " << #x << std::endl; x;}
#define REPORT_STR(x) {WriteTab(out_time) << x << std::endl;}
#define REPORT_VAL(x) {WriteTab(out_time) << #x << " = " << x << std::endl;}
#else
#define REPORT_MPI(x) x
#define REPORT_STR(x) 
#define REPORT_VAL(x)
#endif

namespace INMOST
{
	Storage::integer Mesh::TotalNumberOf(ElementType mask)
	{
		Storage::integer number = 0, ret = 0;
		for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			if( it->GetStatus() != Element::Ghost ) number++;
		ret = number;
#if defined(USE_MPI)
		MPI_Allreduce(&number,&ret,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm);
#endif
		return ret;
	}
	Storage::integer Mesh::EnumerateSet(ElementSet * set, Tag num_tag, Storage::integer start)
	{
		Storage::integer shift = 0, ret = 0;
		ElementType mask = CELL | FACE | EDGE | NODE;
#if defined(USE_MPI)
		Storage::integer number = 0;
		for(ElementSet::iterator it = set->begin(); it != set->end(); it++)
			if( it->GetStatus() != Element::Ghost )
				number++;
		MPI_Scan(&number,&shift,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm);
		shift -= number;
#endif
		shift += start;
		for(ElementSet::iterator it = set->begin(); it != set->end(); it++)
			if( it->GetStatus() != Element::Ghost )
				it->Integer(num_tag) = shift++;
		ExchangeData(num_tag,mask);
		ret = shift;
#if defined(USE_MPI)
		MPI_Bcast(&ret,1,INMOST_MPI_DATA_INTEGER_TYPE,GetProcessorsNumber()-1,comm);
		//MPI_Allreduce(&shift,&ret,1,INMOST_DATA_INTEGER_TYPE,MPI_MAX,comm);
#endif
		return ret;
	}
	Storage::integer Mesh::Enumerate(std::vector<Element *> set, Tag num_tag, Storage::integer start)
	{
		Storage::integer shift = 0, ret = 0;
		ElementType mask = CELL | FACE | EDGE | NODE;
#if defined(USE_MPI)
		Storage::integer number = 0;
		for(std::vector<Element *>::iterator it = set.begin(); it != set.end(); it++)
			if( (*it)->GetStatus() != Element::Ghost )
				number++;
		MPI_Scan(&number,&shift,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm);
		shift -= number;
#endif
		shift += start;
		for(std::vector<Element *>::iterator it = set.begin(); it != set.end(); it++)
			if( (*it)->GetStatus() != Element::Ghost )
				(*it)->Integer(num_tag) = shift++;
		ExchangeData(num_tag,mask);
		ret = shift;
#if defined(USE_MPI)
		MPI_Bcast(&ret,1,INMOST_MPI_DATA_INTEGER_TYPE,GetProcessorsNumber()-1,comm);
		//MPI_Allreduce(&shift,&ret,1,INMOST_DATA_INTEGER_TYPE,MPI_MAX,comm);
#endif
		return ret;
	}
	Storage::integer Mesh::Enumerate(ElementType mask, Tag num_tag, Storage::integer start)
	{
		Storage::integer shift = 0, ret = 0;
#if defined(USE_MPI)
		Storage::integer number = 0;
		for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			if( it->GetStatus() != Element::Ghost )
				number++;
		MPI_Scan(&number,&shift,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm);
		shift -= number;
#endif
		shift += start;
		for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			if( it->GetStatus() != Element::Ghost )
				it->Integer(num_tag) = shift++;
		ExchangeData(num_tag,mask);
		ret = shift;
#if defined(USE_MPI)
		MPI_Bcast(&ret,1,INMOST_MPI_DATA_INTEGER_TYPE,GetProcessorsNumber()-1,comm);
		//MPI_Allreduce(&shift,&ret,1,INMOST_DATA_INTEGER_TYPE,MPI_MAX,comm);
#endif
		return ret;
	}
	
	Storage::real Mesh::Integrate(Storage::real input)
	{
		Storage::real output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_REAL_TYPE,MPI_SUM,comm);
#endif
		return output;
	}
	
	Storage::integer Mesh::Integrate(Storage::integer input)
	{
		Storage::integer output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm);
#endif
		return output;
	}
	
	Storage::integer Mesh::ExclusiveSum(Storage::integer input)
	{
		Storage::integer output = 0;
#if defined(USE_MPI)
		MPI_Scan(&input,&output,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm);
		output -= input;
#endif
		return output;
	}
	
	Storage::real Mesh::Integrate(Tag t, ElementType mask)
	{
		Storage::real output = 0, input = 0;
		for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			if( it->GetStatus() != Element::Ghost )
				input += it->Real(t);
		output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_REAL_TYPE,MPI_SUM,comm);
#endif
		return output;
	}
	
	
	
	
	void DefaultUnpack(Tag tag, Element * element, INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size == 0 )
		{
			if( tag.isDefined(element->GetElementType()) ) 
			{
				if( tag.isSparse(element->GetElementType()) )
				{
					if( element->HaveData(tag) )
						element->DelData(tag); 
				}
				else throw Impossible;
			}
			return;
		}
		if( !element->HaveData(tag) )
			element->SetDataSize(tag,size);
		else if( size != element->GetDataSize(tag) )
		{
			if( tag.GetSize() == ENUMUNDEF )
				element->SetDataSize(tag,size);
			else throw Impossible;
		}
		element->SetData(tag,0,size,data);
	}

	void ArraySetMarker(std::vector<Element *> & arr, MIDType marker)
	{
		for(std::vector<Element *>::iterator it = arr.begin(); it != arr.end(); it++)
			(*it)->SetMarker(marker);
	}
	
	void ArrayRemMarker(std::vector<Element *> & arr, MIDType marker)
	{
		for(std::vector<Element *>::iterator it = arr.begin(); it != arr.end(); it++)
			(*it)->RemMarker(marker);
	}


	


	INMOST_MPI_Comm Mesh::GetCommunicator()
	{
#if defined(USE_MPI)
		return comm;
#else
		return 0;
#endif
	}
	
	int Mesh::GetProcessorRank()
	{
#if defined(USE_MPI)
		int rank;
		MPI_Comm_rank(comm,&rank);
		return rank;
#else
		return 0;
#endif
	}
	
	int Mesh::GetProcessorsNumber()
	{
#if defined(USE_MPI)
		int size;
		MPI_Comm_size(comm,&size);
		return size;
#else
		return 1;
#endif
	}
	
	void Mesh::SetParallelStrategy(int strategy)
	{
		if( strategy < 0 || strategy > 3 ) throw NotImplemented;
		parallel_strategy = strategy;
	}
	
	int Mesh::GetParallelStrategy()
	{
		return parallel_strategy;
	}
	
	void Mesh::Initialize(int * argc, char *** argv)
	{
#if defined(USE_MPI)
		int test;
		MPI_Initialized(&test);
        if( test == 0 ) MPI_Init(argc,argv);
#endif
	}
	void Mesh::Finalize()
	{
#if defined(USE_MPI)
		MPI_Finalize();
#endif
	}
	
	
	
	void Mesh::SetCommunicator(INMOST_MPI_Comm _comm)
	{

		tag_shared = CreateTag("PROTECTED_P_SHARED",DATA_BULK,CELL | FACE | EDGE | NODE,NONE,1);
		tag_owner = CreateTag("P_OWNER_PROCESSOR",DATA_INTEGER, CELL | FACE | EDGE | NODE,NONE,1);
		tag_processors = CreateTag("P_PROCESSORS_LIST",DATA_INTEGER, MESH | NODE | EDGE | FACE | CELL,NONE);
		tag_layers = CreateTag("P_LAYERS",DATA_INTEGER,MESH,NONE,1);
		tag_bridge = CreateTag("P_BRIDGE",DATA_INTEGER,MESH,NONE,1);
		tag_sendto = CreateTag("PROTECTED_P_SENDTO",DATA_INTEGER, CELL | FACE | EDGE | NODE, CELL | FACE | EDGE | NODE);

#if defined(USE_MPI)
		
		
		parallel_strategy = 1;

		

		Mesh::Initialize(NULL,NULL);
		//~ MPI_Comm_dup(_comm,&comm);
		comm = _comm;
		{
			INMOST_DATA_BIG_ENUM_TYPE t = pmid;
			REPORT_MPI(MPI_Allreduce(&t,&parallel_mesh_unique_id,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,MPI_MAX,MPI_COMM_WORLD));
			pmid = parallel_mesh_unique_id+1;
		}
		m_state = Mesh::Parallel;
#if defined(USE_PARALLEL_WRITE_TIME)
		std::stringstream temp;
		temp << "time_" << GetProcessorRank() << ".txt";
		out_time.open(temp.str().c_str(),std::ios::out);
		tab = 0;
#endif

#if defined(USE_MPI2)
		int err;
		REPORT_MPI(err = MPI_Alloc_mem(GetProcessorsNumber()*sizeof(unsigned)*2,MPI_INFO_NULL,&shared_space));
		if( err )
		{
			int errclass;
			MPI_Error_class(err,&errclass);
			std::cout << "Cannot allocate shared space of size " << GetProcessorsNumber()*sizeof(unsigned)*2 << " reason is "; 
			switch(err)
			{
				case MPI_SUCCESS: std::cout << "success"; break;
				case MPI_ERR_INFO: std::cout << "bad info"; break;
				case MPI_ERR_ARG: std::cout << "bad argument"; break;
				case MPI_ERR_NO_MEM: std::cout << "no memory"; break;
			}
			std::cout << std::endl;
			MPI_Abort(comm,err);
		}
		REPORT_MPI(MPI_Win_create(shared_space,sizeof(unsigned)*GetProcessorsNumber()*2,sizeof(unsigned),MPI_INFO_NULL,comm,&window));
#endif
#endif
	}
	
#if defined(USE_PARALLEL_WRITE_TIME)
	void Mesh::Enter() { tab++; }
	void Mesh::Exit() {tab--; }
	std::fstream & Mesh::WriteTab(std::fstream & f)
	{
		for(int i = 0; i < tab; i++)
			f << "   ";
		return f;
	}
	std::fstream & Mesh::GetStream()
	{
		return out_time;
	}
#endif
	
	int comp_mapping(const void * pa, const void * pb)
	{
		return static_cast<const int *>(pa)[0] - static_cast<const int *>(pb)[0];
	}
	
	int comp_mapping_search(const void * pa,void * udata)
	{
		return static_cast<const int *>(pa)[0] - static_cast<int *>(udata)[0];
	}
		
	void determine_my_procs_low(adjacent<Element> & subelements,Storage::integer_array & my_procs,Tag & procs_tag)
	{
		if( subelements.empty() ) return;
		adjacent<Element>::iterator i = subelements.begin();
		std::vector<Storage::integer> result;
		Storage::integer_array p = i->IntegerArray(procs_tag);
		result.insert(result.begin(),p.begin(),p.end());
		i++;
		while(i != subelements.end())
		{
			Storage::integer_array q = i->IntegerArray(procs_tag);
			std::vector<Storage::integer> intersection(result.size());
			std::vector<Storage::integer>::iterator qt = std::set_intersection(result.begin(),result.end(),q.begin(),q.end(),intersection.begin());
			intersection.resize(qt-intersection.begin());
			result.swap(intersection);
			i++;
		}
		my_procs.clear();
		my_procs.insert(my_procs.begin(),result.begin(),result.end());
	}
	
	void determine_my_procs_high(adjacent<Element> & overelements,Storage::integer_array & my_procs,Tag & procs_tag)
	{
		if( overelements.empty() ) return;
		adjacent<Element>::iterator i = overelements.begin();
		std::vector<Storage::integer> result;
		Storage::integer_array p = i->IntegerArray(procs_tag);
		result.insert(result.begin(),p.begin(),p.end());
		i++;
		while(i != overelements.end())
		{
			Storage::integer_array q = i->IntegerArray(procs_tag);
			std::vector<Storage::integer> intersection(result.size()+q.size());
			std::vector<Storage::integer>::iterator qt = std::set_union(result.begin(),result.end(),q.begin(),q.end(),intersection.begin());
			intersection.resize(qt-intersection.begin());
			result.swap(intersection);
			i++;
		}
		my_procs.clear();
		my_procs.insert(my_procs.begin(),result.begin(),result.end());
	}
	
	
	void Mesh::ResolveShared()
	{
#if defined(USE_MPI)
		if( m_state == Mesh::Serial ) SetCommunicator(INMOST_MPI_COMM_WORLD);
		unsigned int dim = GetDimensions();
		int sendsize;
		int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber();
#if defined(USE_PARALLEL_STORAGE)
		shared_elements.clear();
		ghost_elements.clear();
#endif

#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
		std::vector<Storage::real> bbox(dim*2);
		std::vector<Storage::real> bboxs(mpisize*dim*2);
		std::vector<unsigned> hashs(mpisize);
		
		std::map<GeometricData,ElementType> table;
		for(ElementType etype = EDGE; etype <= CELL; etype = etype << 1)
			if( !HaveGeometricData(CENTROID,etype) )
				table[CENTROID] |= etype;
		PrepareGeometricData(table);
		
		for(unsigned k = 0; k < dim; k++)
		{
			bbox[k] = 1e20;
			bbox[k+dim] = -1e20;
		}
		
		
		for(iteratorNode it = BeginNode(); it != EndNode(); it++)
		{
			Storage::integer_array arr = it->IntegerArrayDV(tag_processors);
			arr.resize(1);
			arr[0] = mpirank;
		}
		

		std::vector< INMOST_DATA_BULK_TYPE > exch_data;
		std::vector< INMOST_DATA_REAL_TYPE > pack_real;
		std::vector< INMOST_DATA_REAL_TYPE > unpack_real;
		pack_real.reserve(NumberOfNodes()*dim);
		
		std::vector<Element *> sorted_nodes;
		sorted_nodes.reserve(NumberOfNodes());
		for(iteratorElement n = BeginElement(NODE); n != EndElement(); n++)
			sorted_nodes.push_back(&*n);
		
		qsort(&sorted_nodes[0],sorted_nodes.size(),sizeof(Element *),CompareElementsCCentroid);



		for(std::vector<Element *>::iterator it = sorted_nodes.begin(); it != sorted_nodes.end(); it++)
		{
			Storage::real_array arr = (*it)->getAsNode()->Coords();
			pack_real.insert(pack_real.end(),arr.begin(),arr.end());
			for(unsigned k = 0; k < dim; k++)
			{
				if( arr[k] < bbox[k] ) bbox[k] = arr[k];
				if( arr[k] > bbox[k+dim] ) bbox[k+dim] = arr[k];
			}
		}
		
		MPI_Pack_size(pack_real.size(),INMOST_MPI_DATA_REAL_TYPE,comm,&sendsize);
		exch_data.resize(sendsize);
		int position = 0;
		if( sendsize > 0 ) MPI_Pack(&pack_real[0],pack_real.size(),INMOST_MPI_DATA_REAL_TYPE,&exch_data[0],exch_data.size(),&position,comm);
		//pack_real.clear();
		
		unsigned hash = hashfunci32(reinterpret_cast<char *>(&exch_data[0]),exch_data.size());
		
		REPORT_MPI(MPI_Allgather(&hash,1,MPI_UNSIGNED,&hashs[0],1,MPI_UNSIGNED,comm));
		
		bool flag = true;
		for(int k = 1; k < mpisize; k++)
			if( hashs[k] != hashs[0] )
			{
				flag = false;
				break;
			}
		REPORT_VAL(flag);
		if( flag ) // this grid is duplicated over all processors - using fast algorithm
		{
			REPORT_VAL(mpisize);
			for(Mesh::iteratorElement it = BeginElement(CELL | EDGE | FACE | NODE); it != EndElement(); it++)
			{
				Storage::integer_array arr = it->IntegerArrayDV(tag_processors);
				arr.resize(mpisize);
				for(int k = 0; k < mpisize; k++)
					arr[k] = k;
				it->IntegerDF(tag_owner) = 0;
				if( mpirank == 0 )
					it->BulkDF(tag_shared) = Element::Shared;
				else
					it->BulkDF(tag_shared) = Element::Ghost;
			}
			ComputeSharedProcs();
			RecomputeParallelStorage(CELL | EDGE | FACE | NODE);
			AssignGlobalID(CELL | EDGE | FACE | NODE);
#if defined(USE_PARALLEL_STORAGE)
			for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
				for(int i = 0; i < 4; i++)
					qsort(&it->second[i][0],it->second[i].size(),sizeof(Element *),CompareElementsCGID);
			for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)			
				for(int i = 0; i < 4; i++)
					qsort(&it->second[i][0],it->second[i].size(),sizeof(Element *),CompareElementsCGID);
#endif
		}
		else
		
		{
		
			
			{
				Storage::real epsilon = GetEpsilon();
				
				REPORT_MPI(MPI_Allgather(&bbox[0],dim*2,INMOST_MPI_DATA_REAL_TYPE,&bboxs[0],dim*2,INMOST_MPI_DATA_REAL_TYPE,comm));
				
				//determine which bboxes i intersect
				std::vector<int> procs;
				
				
				for(int k = 0; k < mpisize; k++)
					if( k != mpirank )
					{
						bool flag = true;
						for(unsigned q = 0; q < dim; q++)
							flag &= !((bbox[q] > bboxs[k*dim*2+q+dim]) || (bbox[dim+q] < bboxs[k*dim*2+q]));
						if( flag ) procs.push_back(k);
					}
				
				
				std::vector<int> sendsizeall(mpisize*2);
				int sendsize2[2] = {sendsize,pack_real.size()};
				REPORT_MPI(MPI_Allgather(sendsize2,2,MPI_INT,&sendsizeall[0],2,MPI_INT,comm));
				
				{
					std::vector< std::pair< int, std::vector<INMOST_DATA_BULK_TYPE> > > send_buffs(procs.size());
					std::vector< std::pair< int, std::vector<INMOST_DATA_BULK_TYPE> > > recv_buffs(procs.size());
					std::vector< MPI_Request > send_reqs, recv_reqs;
					
					for(unsigned k = 0; k < procs.size(); k++)
					{
						send_buffs[k].first = procs[k];
						send_buffs[k].second = exch_data;
						recv_buffs[k].first = procs[k];
						recv_buffs[k].second.resize(sendsizeall[procs[k]*2]);
					}
					
					//PrepareReceiveInner(send_buffs,recv_buffs);
					ExchangeBuffersInner(send_buffs,recv_buffs,send_reqs,recv_reqs);
					
					
					std::vector<int> done;
					while( !(done = FinishRequests(recv_reqs)).empty() )
					{
						for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
						{
							int position = 0;
							unpack_real.resize(sendsizeall[recv_buffs[*qt].first*2+1]);
							MPI_Unpack(&recv_buffs[*qt].second[0],recv_buffs[*qt].second.size(),&position,&unpack_real[0],unpack_real.size(),INMOST_MPI_DATA_REAL_TYPE,comm);
							std::vector<Storage::real>::iterator it1 = pack_real.begin() , it2 = unpack_real.begin();
							while(it1 != pack_real.end() && it2 != unpack_real.end() )
							{
								int res = 0;
								for(unsigned k = 0; k < dim; k++)
									if( fabs((*(it1+k))-(*(it2+k))) > epsilon )
									{
										if( (*(it1+k)) < (*(it2+k)) ) res = -1;
										else res = 1;
										break;
									}
								if( res < 0 ) 
									it1 += dim;
								else if( res > 0 ) 
									it2 += dim;
								else
								{
									sorted_nodes[(it1-pack_real.begin())/dim]->IntegerArrayDV(tag_processors).push_back(recv_buffs[*qt].first);
									it1 += dim;
									it2 += dim;
								}
							}
						}
					}
					REPORT_MPI(MPI_Waitall(send_reqs.size(),&send_reqs[0],MPI_STATUSES_IGNORE));
				}			
				
				for(Mesh::iteratorElement it = BeginElement(NODE); it != EndElement(); it++)
				{
					int owner;
					Storage::integer_array v = it->IntegerArrayDV(tag_processors);
					std::sort(v.begin(),v.end());
					if( v.empty() )
					{
						owner = mpirank;
						v.push_back(mpirank);
					}
					else
						owner = std::min(mpirank,v[0]);
					
					it->IntegerDF(tag_owner) = owner;
					
					if( mpirank == owner )
					{
						if( v.size() == 1 )
							it->BulkDF(tag_shared) = Element::Owned;
						else
							it->BulkDF(tag_shared) = Element::Shared;
					}
					else
						it->BulkDF(tag_shared) = Element::Ghost;
				}
				
				ComputeSharedProcs();
				RecomputeParallelStorage(NODE);
				AssignGlobalID(NODE);
				
#if defined(USE_PARALLEL_STORAGE)
				for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
					qsort(&it->second[0][0],it->second[0].size(),sizeof(Element *),CompareElementsCGID);
				for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)			
					qsort(&it->second[0][0],it->second[0].size(),sizeof(Element *),CompareElementsCGID);
#endif
			}
			
			for(ElementType current_mask = EDGE; current_mask <= CELL; current_mask = current_mask << 1 )
			{
				
				int owner;
				Element::Status estat;
				
				//Determine what processors potentially share the element
				for(Mesh::iteratorElement it = BeginElement(current_mask); it != EndElement(); it++)
				{
					adjacent<Element> sub = it->getAdjElements(current_mask >> 1);
					Storage::integer_array p = it->IntegerArrayDV(tag_processors);
					determine_my_procs_low(sub,p,tag_processors);
					if( p.empty() ) p.push_back(mpirank);
				}
				
				//Initialize mapping that helps get local id by global id
				std::vector<int> mapping;
				for(Mesh::iteratorElement it = BeginElement(current_mask >> 1); it != EndElement(); it++)
				{
					mapping.push_back(it->GlobalID());
					mapping.push_back(it->LocalID());
				}
				qsort(&mapping[0],mapping.size()/2,sizeof(int)*2,comp_mapping);
				//Initialize arrays
				Storage::integer_array procs = IntegerArrayDV(tag_processors);
				Storage::integer_array::iterator p = procs.begin();
				std::vector<int> message_send;
				std::vector< std::vector<int> > message_recv(procs.size());
				std::vector< std::vector<Element *> > elements(procs.size());
				
				std::vector< MPI_Request > send_reqs,recv_reqs;
				std::vector< std::pair<int, std::vector<INMOST_DATA_BULK_TYPE> > > send_buffs(procs.size()), recv_buffs(procs.size());
				
				//Gather all possible shared elements and send global ids of their connectivity to the neighbouring proccessors
				for(p = procs.begin(); p != procs.end(); p++)
				{
					int m = p-procs.begin();
					{
						message_send.clear();
						message_send.push_back(0);
						message_send.push_back(0);
						for(Mesh::iteratorElement it = BeginElement(current_mask); it != EndElement(); it++)
						{
							Storage::integer_array pr = it->IntegerArrayDV(tag_processors);
							if( std::binary_search(pr.begin(),pr.end(),*p) )
							{
								adjacent<Element> sub = it->getAdjElements(current_mask >> 1);
								if( sub.size() == 0 ) throw Impossible;
								message_send.push_back(sub.size());
								for(adjacent<Element>::iterator kt = sub.begin(); kt != sub.end(); kt++)
									message_send.push_back(kt->GlobalID());
								message_send[1]++;
								elements[m].push_back(&*it);
							}
						}
						message_send[0] = message_send.size();
						MPI_Pack_size(message_send.size(),MPI_INT,comm,&sendsize);
						send_buffs[m].first = *p;
						send_buffs[m].second.resize(sendsize);
						int position = 0;
						MPI_Pack(&message_send[0],message_send.size(),MPI_INT,&send_buffs[m].second[0],send_buffs[m].second.size(),&position,comm);
						send_buffs[m].second.resize(position);
						recv_buffs[m].first = *p;
					}
				}
				
				PrepareReceiveInner(UnknownSize,send_buffs,recv_buffs);
				ExchangeBuffersInner(send_buffs,recv_buffs,send_reqs,recv_reqs);
				
				
				std::vector<int> done;
				
				while( !(done = FinishRequests(recv_reqs)).empty() )
				{
					for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
					{
						int position = 0;
						int size;
						int pos = -1;
						for(p = procs.begin(); p != procs.end(); p++)
							if( *p == recv_buffs[*qt].first )
							{
								pos = p - procs.begin();
								break;
							}
						if( pos == -1 ) throw Impossible;
						MPI_Unpack(&recv_buffs[*qt].second[0],recv_buffs[*qt].second.size(),&position,&size,1,MPI_INT,comm);
						message_recv[pos].resize(size-1);
						MPI_Unpack(&recv_buffs[*qt].second[0],recv_buffs[*qt].second.size(),&position,&message_recv[pos][0],message_recv[pos].size(),MPI_INT,comm);
					}
				}
				//Now find the difference of local elements with given processor number and remote elements
				for(p = procs.begin(); p != procs.end(); p++)
				{
					int m = p-procs.begin();
					std::vector<Element *> remote_elements;
					int pos = 0;
					if( message_recv[m].empty() ) continue;
					int num_remote_elements = message_recv[m][pos++];
					for(int i = 0; i < num_remote_elements; i++)
					{
						int conn_size = message_recv[m][pos++], flag = 1;
						
						std::vector<Element *> sub_elements;
						for(int j = 0; j < conn_size; j++)
						{
							int global_id = message_recv[m][pos++];
							int find = binary_search(&mapping[0],mapping.size()/2,sizeof(int)*2,comp_mapping_search,&global_id);
							if( find == -1 ) 
							{
								flag = 0;
								pos += conn_size-j-1;
								break;
							}
							int find_local_id = mapping[2*find+1];
							switch(current_mask)
							{
								case EDGE: sub_elements.push_back(nodes[find_local_id]); break;
								case FACE: sub_elements.push_back(edges[find_local_id]); break;
								case CELL: sub_elements.push_back(faces[find_local_id]); break;
								default: throw Impossible;
							}
							
						}
						if( flag )
						{
							Element * e = FindSharedAdjacency(&sub_elements[0],sub_elements.size());
							if( e == NULL ) continue;
							remote_elements.push_back(e);
						}
					}
					std::sort(remote_elements.begin(),remote_elements.end());
					std::sort(elements[m].begin(),elements[m].end());
					std::vector<Element *> result;
					std::vector<Element *>::iterator set_end;
					
					result.resize(elements[m].size());
					set_end = std::set_difference(elements[m].begin(),elements[m].end(),remote_elements.begin(),remote_elements.end(), result.begin());
					result.resize(set_end-result.begin());				
					
					//elements in result are wrongly marked as ghost
					for(std::vector<Element *>::iterator qt = result.begin(); qt != result.end(); qt++)
					{
						Storage::integer_array pr = (*qt)->IntegerArrayDV(tag_processors);
						Storage::integer_array::iterator find = std::lower_bound(pr.begin(),pr.end(),*p);
						pr.erase(find);
					}
				}
				//Now mark all the processors status
				for(Mesh::iteratorElement it = BeginElement(current_mask); it != EndElement(); it++)
				{
					Storage::integer_array pr = it->IntegerArrayDV(tag_processors);
					if( pr.empty() )
					{
						owner = mpirank;
						pr.push_back(mpirank);
					}
					else
						owner = std::min(mpirank,pr[0]);
					it->IntegerDF(tag_owner) = owner;
					if( mpirank == owner )
					{
						if( pr.size() == 1 )
							estat = Element::Owned;
						else
							estat = Element::Shared;
					}
					else
						estat = Element::Ghost;
					it->BulkDF(tag_shared) = estat;
				}
				RecomputeParallelStorage(current_mask);
				REPORT_MPI(MPI_Waitall(send_reqs.size(),&send_reqs[0],MPI_STATUSES_IGNORE));
				AssignGlobalID(current_mask);
				
#if defined(USE_PARALLEL_STORAGE)
				for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
					qsort(&it->second[ElementNum(current_mask)][0],it->second[ElementNum(current_mask)].size(),sizeof(Element *),CompareElementsCGID);
				for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)			
					qsort(&it->second[ElementNum(current_mask)][0],it->second[ElementNum(current_mask)].size(),sizeof(Element *),CompareElementsCGID);
#endif
			}	
		}
		//RemoveGhost();
		RemoveGeometricData(table);
#else
		AssignGlobalID(CELL | FACE | EDGE | NODE);
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}
	
	void DeleteUnpack(Tag tag,Element * e,INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size ) 
		{
			int old_size;
			if( e->HaveData(tag) ) old_size = e->GetDataSize(tag);
			else old_size = 0;
			e->SetDataSize(tag,old_size+size);
			e->SetData(tag,old_size,size,data);
		}
	}
	
	void Mesh::RemoveGhost()
	{
		if( m_state == Mesh::Serial ) return;
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
		std::vector<Element *> delete_elements;
#if defined(USE_PARALLEL_STORAGE)
		std::map< int, std::vector< Element *> > del_shared, del_ghost;
#endif
		int mpirank = GetProcessorRank();
		Tag tag_delete = CreateTag("TEMPORARY_DELETE_GHOST_TAG",DATA_INTEGER,FACE | EDGE | NODE,FACE | EDGE | NODE);
		
		for(cells_container::iterator it = cells.begin(); it != cells.end(); it++) if( (*it) != NULL )
		{
			Element::Status estat = (*it)->BulkDF(tag_shared);
			if( estat == Element::Ghost ) 
			{
#if defined(USE_PARALLEL_STORAGE)
				del_ghost[(*it)->IntegerDF(tag_owner)].push_back(*it);
#endif
				delete_elements.push_back(*it);
			}
			else if( estat == Element::Shared )
			{
				(*it)->Bulk(tag_shared) = Element::Owned;
				Storage::integer_array it_procs = (*it)->IntegerArrayDV(tag_processors);
				
#if defined(USE_PARALLEL_STORAGE)
				for(Storage::integer_array::iterator vit = it_procs.begin(); vit != it_procs.end(); vit++)
					if( *vit != mpirank )
						del_shared[*vit].push_back(*it);
#endif
				
				it_procs.resize(1);
				it_procs[0] = mpirank;
			}
		}
		
#if defined(USE_PARALLEL_STORAGE)
		for(std::map< int, std::vector<Element *> >::iterator it = del_ghost.begin(); it != del_ghost.end(); it++)
		{
			std::vector<Element * > & ref = ghost_elements[it->first][ElementNum(CELL)];
			qsort(&it->second[0],it->second.size(),sizeof(Element *),CompareElementsCGID);
			std::vector<Element *> result(ref.size());
			std::vector<Element *>::iterator end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),LessElementsGID);
			result.resize(end-result.begin());
			ref.swap(result);
		}
		del_ghost.clear();
		for(std::map< int, std::vector<Element *> >::iterator it = del_shared.begin(); it != del_shared.end(); it++)
		{
			std::vector<Element * > & ref = shared_elements[it->first][ElementNum(CELL)];
			qsort(&it->second[0],it->second.size(),sizeof(Element *),CompareElementsCGID);
			std::vector<Element *> result(ref.size());
			std::vector<Element *>::iterator end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),LessElementsGID);
			result.resize(end-result.begin());
			ref.swap(result);
		}
		del_shared.clear();
#endif

		for(std::vector<Element *>::iterator it = delete_elements.begin(); it != delete_elements.end(); it++) delete *it;
		delete_elements.clear();
		
		for(ElementType mask = FACE; mask >= NODE; mask = mask >> 1)
		{
			for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			{
				Element::Status estat = it->BulkDF(tag_shared);
				if( estat == Element::Owned || estat == Element::Shared ) continue;
				if( it->nbAdjElements(mask<<1) == 0 ) it->Integer(tag_delete) = mpirank;
			}
			ReduceData(tag_delete,mask,DeleteUnpack);
			ExchangeData(tag_delete,mask);
			for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			{
				Element::Status estat = it->BulkDF(tag_shared);
				if( estat == Element::Owned ) continue;
				if( estat == Element::Ghost && it->nbAdjElements(mask<<1) == 0 ) 
				{
#if defined(USE_PARALLEL_STORAGE)
					del_ghost[it->IntegerDF(tag_owner)].push_back(&*it);
#endif
					delete_elements.push_back(&*it);
				}
				else if( it->HaveData(tag_delete) )
				{
					Storage::integer_array procs = it->IntegerArrayDV(tag_processors);
					Storage::integer_array del_procs = it->IntegerArray(tag_delete);
					std::vector<Storage::integer> result(procs.size());
					std::sort(del_procs.begin(),del_procs.end());								
#if defined(USE_PARALLEL_STORAGE)				
					for(Storage::integer_array::iterator vit = del_procs.begin(); vit != del_procs.end(); vit++)
						del_shared[*vit].push_back(&*it);
#endif						
					std::vector<Storage::integer>::iterator end = std::set_difference(procs.begin(),procs.end(),del_procs.begin(),del_procs.end(),result.begin());
					result.resize(end-result.begin());
					procs.clear();
					procs.insert(procs.begin(),result.begin(),result.end());
					
					if( procs.size() == 1 && procs[0] == mpirank )
						it->Bulk(tag_shared) = Element::Owned;
				}
			}
			
#if defined(USE_PARALLEL_STORAGE)
			for(std::map< int, std::vector<Element *> >::iterator it = del_ghost.begin(); it != del_ghost.end(); it++)
			{
				std::vector<Element * > & ref = ghost_elements[it->first][ElementNum(mask)];
				qsort(&it->second[0],it->second.size(),sizeof(Element *),CompareElementsCGID);
				std::vector<Element *> result(ref.size());
				std::vector<Element *>::iterator end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),LessElementsGID);
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_ghost.clear();
			for(std::map< int, std::vector<Element *> >::iterator it = del_shared.begin(); it != del_shared.end(); it++)
			{
				std::vector<Element * > & ref = shared_elements[it->first][ElementNum(mask)];
				qsort(&it->second[0],it->second.size(),sizeof(Element *),CompareElementsCGID);
				std::vector<Element *> result(ref.size());
				std::vector<Element *>::iterator end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),LessElementsGID);
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_shared.clear();
#endif
			for(std::vector<Element *>::iterator it = delete_elements.begin(); it != delete_elements.end(); it++) delete *it;
			delete_elements.clear();
		}		
		DeleteTag(tag_delete);
		ComputeSharedProcs();
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}
	
	void Mesh::RemoveGhostElements(std::vector<Element *> ghost)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
		int mpirank = GetProcessorRank();
		Tag tag_delete = CreateTag("TEMPORARY_DELETE_GHOST_ELEMENTS_TAG",DATA_INTEGER,CELL | FACE | EDGE | NODE,CELL | FACE | EDGE | NODE);
		std::vector<Element *> delete_elements;
#if defined(USE_PARALLEL_STORAGE)
		std::map< int, std::vector< Element *> > del_shared, del_ghost;
#endif
		
		for(ElementType mask = CELL; mask >= NODE; mask = mask >> 1)
		{
			if( mask & CELL )
			{
				for(std::vector<Element *>::iterator it = ghost.begin(); it != ghost.end(); it++) if( (*it) != NULL )
					if( (*it)->Bulk(tag_shared) == Element::Ghost ) (*it)->Integer(tag_delete) = mpirank;
			}
			else 
			{
				for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
				{
					Element::Status estat = it->BulkDF(tag_shared);
					if( estat == Element::Owned || estat == Element::Shared ) continue;
					if( it->nbAdjElements(mask<<1) == 0 ) it->Integer(tag_delete) = mpirank;
				}
			}
			ReduceData(tag_delete,mask,DeleteUnpack);
			ExchangeData(tag_delete,mask);
			for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			{
				Element::Status estat = it->BulkDF(tag_shared);
				if( estat == Element::Owned ) continue;
				if( it->HaveData(tag_delete) )
				{
						Storage::integer_array del_procs = it->IntegerArray(tag_delete);
						std::sort(del_procs.begin(),del_procs.end());
						
						if( estat == Element::Ghost && std::binary_search(del_procs.begin(),del_procs.end(),mpirank) )
						{
#if defined(USE_PARALLEL_STORAGE)
							del_ghost[it->IntegerDF(tag_owner)].push_back(&*it);
#endif
							delete_elements.push_back(&*it);
						}
						else
						{
							Storage::integer_array procs = it->IntegerArrayDV(tag_processors);
							std::vector<Storage::integer> result(procs.size());
#if defined(USE_PARALLEL_STORAGE)
							if( estat == Element::Shared )
							{
								for(Storage::integer_array::iterator vit = del_procs.begin(); vit != del_procs.end(); vit++)
									del_shared[*vit].push_back(&*it);
							}
#endif						
							std::vector<Storage::integer>::iterator end = std::set_difference(procs.begin(),procs.end(),del_procs.begin(),del_procs.end(),result.begin());
							result.resize(end-result.begin());
							procs.clear();
							procs.insert(procs.begin(),result.begin(),result.end());
							
							if( procs.size() == 1 && procs[0] == mpirank )
								it->Bulk(tag_shared) = Element::Owned;						
						}
				}
			}
#if defined(USE_PARALLEL_STORAGE)
			for(std::map< int, std::vector<Element *> >::iterator it = del_ghost.begin(); it != del_ghost.end(); it++)
			{
				std::vector<Element * > & ref = ghost_elements[it->first][ElementNum(mask)];
				qsort(&it->second[0],it->second.size(),sizeof(Element *),CompareElementsCGID);
				std::vector<Element *> result(ref.size());
				std::vector<Element *>::iterator end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),LessElementsGID);
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_ghost.clear();
			for(std::map< int, std::vector<Element *> >::iterator it = del_shared.begin(); it != del_shared.end(); it++)
			{
				std::vector<Element * > & ref = shared_elements[it->first][ElementNum(mask)];
				qsort(&it->second[0],it->second.size(),sizeof(Element *),CompareElementsCGID);
				std::vector<Element *> result(ref.size());
				std::vector<Element *>::iterator end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),LessElementsGID);
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_shared.clear();
#endif
			for(std::vector<Element *>::iterator it = delete_elements.begin(); it != delete_elements.end(); it++) delete *it;
			delete_elements.clear();			
		}
		DeleteTag(tag_delete);
		ComputeSharedProcs();
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}

	
	void Mesh::AssignGlobalID(ElementType mask)
	{
		
		for(ElementType currenttype = NODE; currenttype <= CELL; currenttype = currenttype << 1 )
			if( (currenttype & mask) && !(currenttype & have_global_id ) )
			{
				tag_global_id = CreateTag("GLOBAL_ID",DATA_INTEGER, currenttype, NONE,1);
				have_global_id |= currenttype;
			}
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
		if( m_state == Parallel )
		{
			INMOST_DATA_BIG_ENUM_TYPE number,shift[4], shift_recv[4], local_shift;
			//int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber();
			int num;
			//std::vector<INMOST_DATA_BIG_ENUM_TYPE> shifts(mpisize*4);
			memset(shift,0,sizeof(INMOST_DATA_BIG_ENUM_TYPE));
			for(ElementType currenttype = NODE; currenttype <= CELL; currenttype = currenttype << 1 )
			{
				if( mask & currenttype )
				{
					num = ElementNum(currenttype);
					number = 0;
					for(Mesh::iteratorElement it = BeginElement(currenttype); it != EndElement(); it++)
					{
						if( it->GetStatus() != Element::Ghost )
							number++;
					}
					shift[num] = number;
				}
			}
			{
				int ierr;
				REPORT_MPI(ierr = MPI_Scan(shift,shift_recv,4,INMOST_MPI_DATA_BIG_ENUM_TYPE,MPI_SUM,GetCommunicator()));
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			}
			for(ElementType currenttype = NODE; currenttype <= CELL; currenttype = currenttype << 1 )
			{
				if( mask & currenttype )
				{
					num = ElementNum(currenttype);
					local_shift = shift_recv[num] - shift[num];
					for(Mesh::iteratorElement it = BeginElement(currenttype); it != EndElement(); it++)
					{
						if( it->GetStatus() != Element::Ghost )
							it->Integer(tag_global_id) = local_shift++;
					}	
				}
			}
			if( tag_global_id.isValid() )ExchangeData(tag_global_id,mask);
		}
		else
		{
#endif
			INMOST_DATA_BIG_ENUM_TYPE number;
			for(ElementType currenttype = NODE; currenttype <= CELL; currenttype = currenttype << 1 )
				if( mask & currenttype )
				{
					number = 0;
					for(Mesh::iteratorElement it = BeginElement(currenttype); it != EndElement(); it++)
						it->Integer(tag_global_id) = number++;
				}
#if defined(USE_MPI)
		}
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif

	}
		
	void Mesh::ComputeSharedProcs()
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
		std::set<int> shared_procs;
		int mpirank = GetProcessorRank();
		for(nodes_container::iterator it = nodes.begin(); it != nodes.end(); it++) if( (*it) != NULL )
		{
			Storage::integer_array p = (*it)->IntegerArrayDV(tag_processors);
			for(Storage::integer_array::iterator kt = p.begin(); kt != p.end(); kt++) shared_procs.insert(*kt);
		}
		std::set<int>::iterator ir = shared_procs.find(mpirank);
		if( ir != shared_procs.end() ) shared_procs.erase(ir);
		Storage::integer_array procs = IntegerArrayDV(tag_processors);
		procs.clear();
		procs.insert(procs.begin(),shared_procs.begin(),shared_procs.end());
		
		REPORT_VAL(procs.size());
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}
	
	void UnpackSkin(Tag tag, Element * e, INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size )
		{
			bool flag = true;
			Storage::integer * recv = static_cast<Storage::integer *>(static_cast<void *>(data));
			Storage::integer_array arr = e->IntegerArray(tag);
			for(Storage::integer_array::iterator it = arr.begin(); it != arr.end(); it++)
				if( *it == recv[0] )
				{
					flag = false;
					break;
				}
			if( flag ) 
			{
				arr.push_back(recv[0]);
				arr.push_back(recv[1]);
			}
		}
	}
	
	void UnpackOnSkin(Tag tag, Element * e, INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size )
		{
			Storage::integer * recv = static_cast<Storage::integer *>(static_cast<void *>(data));
			Storage::integer_array arr = e->IntegerArray(tag);
			arr.push_back(recv[0]);
		}
	}
	
	std::map< int , std::vector<Element *> > Mesh::ComputeSharedSkinSet(ElementType bridge_type)
	{
#if defined(USE_MPI)
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif

		int mpirank = GetProcessorRank();
		Tag tag_skin = CreateTag("TEMPORARY_COMPUTE_SHARED_SKIN_SET",DATA_INTEGER,FACE,FACE);
				
		for(iteratorFace it = BeginFace(); it != EndFace(); it++)
		{
			Element::Status estat = it->BulkDF(tag_shared);
			if( estat == Element::Ghost || estat == Element::Shared )
			{
				Storage::integer_array arr = it->IntegerArray(tag_skin);
				adjacent<Element> adj = it->getAdjElements(CELL);
				for(adjacent<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
				{
					arr.push_back(jt->GlobalID());
					arr.push_back(jt->IntegerDF(tag_owner));
				}
			}
		}
		ReduceData(tag_skin,FACE,UnpackSkin);
		ExchangeData(tag_skin,FACE);

		std::map< int, std::vector<Element *> > skin_faces;
		for(iteratorFace it = BeginFace(); it != EndFace(); it++)
		{
			bool flag = false;
			Storage::integer_array skin_data = it->IntegerArray(tag_skin);
			Element::Status estat = it->BulkDF(tag_shared), estat2;
			if( estat == Element::Owned ) continue;
			if( skin_data.size()/2 < 2 ) continue;
			if( skin_data.size()/2 > 2 ) throw Impossible;
			Cell * c1 = it->BackCell(), * c2 = it->FrontCell();
			if( c1 == NULL && c2 == NULL ) continue;
			else if( c2 == NULL )
			{
				estat = c1->Bulk(tag_shared);
				if( estat == Element::Owned || estat == Element::Shared )
					flag = true;
			}
			else
			{
				estat  = c1->Bulk(tag_shared);
				estat2 = c2->Bulk(tag_shared);
				if( (estat  == Element::Ghost && estat2 == Element::Shared) ||
					(estat2 == Element::Ghost && estat  == Element::Shared) ||
					(estat  == Element::Owned && estat2 == Element::Ghost)  ||
					(estat2 == Element::Owned && estat  == Element::Ghost))
					flag = true;
			}
			
			if( flag )
			{
				
				for(int i = 0; i < 2; i++)
				{
					Storage::integer owner = skin_data[i*2+1];
					if( owner != mpirank )
					{
						skin_faces[owner].push_back(&*it);
						break;
					}
				}	
			}
		}
		DeleteTag(tag_skin);
		
		
#if defined(DEBUG_COMPUTE_SHARED_SKIN_SET)
		{
			Storage::integer keynum;
			FILE * file;
			char keyword[2048];
			if( GetProcessorRank() == 0)
			{
				file = fopen("skin.gmv","wb");
				sprintf(keyword,"gmvinput"); fwrite(keyword,1,8,file);
				sprintf(keyword,"ieee");
				if( sizeof(Storage::real) != 8 || sizeof(Storage::integer) != 8 )
					sprintf(keyword,"ieeei%ldr%ld",sizeof(Storage::integer),sizeof(Storage::real));
				fwrite(keyword,1,8,file);
				sprintf(keyword,"polygons"); fwrite(keyword,1,8,file);
				fclose(file);
			}
			for(int i = 0; i < GetProcessorRank(); i++) MPI_Barrier(comm);
			file = fopen("skin.gmv","ab");
			for(std::map< int, std::vector<Element *> >::iterator it = skin_faces.begin(); it != skin_faces.end(); it++)
			for(std::vector<Element *>::iterator f = it->second.begin(); f != it->second.end(); f++)
			{
				adjacent<Node> fnodes = (*f)->getNodes();
				keynum = GetProcessorRank()*GetProcessorsNumber()+(it->first+1); 
				fwrite(&keynum,sizeof(Storage::integer),1,file);
				keynum = static_cast<Storage::integer>(fnodes.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
				for(adjacent<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
					fwrite(&fn->Coords()[0],sizeof(Storage::real),1,file);
				for(adjacent<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
					fwrite(&fn->Coords()[1],sizeof(Storage::real),1,file);
				for(adjacent<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
					fwrite(&fn->Coords()[2],sizeof(Storage::real),1,file);
			}
			fclose(file);
			for(int i = GetProcessorRank(); i < GetProcessorsNumber(); i++) MPI_Barrier(comm);
			if( GetProcessorRank() == 0 )
			{
				file = fopen("skin.gmv","ab");
				sprintf(keyword,"endpoly"); fwrite(keyword,1,8,file);
				sprintf(keyword,"endgmv"); fwrite(keyword,1,8,file);
				fclose(file);
			}
		}
#endif
		

		
		
		std::map< int , std::vector<Element *> > bridge;
		if( bridge_type & FACE ) bridge.swap(skin_faces);
		else
		{
			std::vector<Element *> all_visited;
			Tag on_skin = CreateTag("TEMPORARY_ON_SKIN",DATA_INTEGER,bridge_type,bridge_type);
			MIDType busy = CreateMarker();
			for(std::map< int , std::vector<Element *> >::iterator kt = skin_faces.begin(); kt != skin_faces.end(); kt++)
			{
				for(std::vector<Element *>::iterator it = kt->second.begin(); it != kt->second.end(); it++)
				{
					adjacent<Element> adj = (*it)->getAdjElements(bridge_type);
					for(adjacent<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
					{
						if( !jt->GetMarker(busy) )
						{
							jt->IntegerArray(on_skin).push_back(mpirank);
							jt->SetMarker(busy);
							all_visited.push_back(&*jt);
						}
					}
				}
			}
			ArrayRemMarker(all_visited,busy);
			ReleaseMarker(busy);
			
			ReduceData(on_skin,bridge_type,UnpackOnSkin);
			ExchangeData(on_skin,bridge_type);
			
			for(std::vector<Element *>::iterator it = all_visited.begin(); it != all_visited.end(); it++)
			{
				Storage::integer_array os = (*it)->IntegerArray(on_skin);
				for(Storage::integer_array::iterator p = os.begin(); p != os.end(); p++)
					if( *p != mpirank )
						bridge[*p].push_back(*it);
			}
			
						
			ReleaseMarker(busy);
		}		
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time  << std::endl;
#endif
		return bridge;		
#else
		return std::map< int, std::vector<Element *> >();
#endif
	}
	
	
	void Mesh::PackTagData(Tag tag, std::vector< std::vector<Element *> > & elements, ElementType mask, std::vector<INMOST_DATA_BULK_TYPE> & buffer)
	{
		if( tag.GetDataType() == DATA_REFERENCE ) return; //NOT IMPLEMENTED
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
		std::vector<Element *>::iterator eit;
		std::vector<INMOST_DATA_BULK_TYPE> array_data_send;
		std::vector<INMOST_DATA_ENUM_TYPE> array_size_send(2);
		array_data_send.reserve(4096);
		array_size_send.reserve(4096);
		unsigned int size = tag.GetSize();
		for(int i = 0; i < 4; i++) if( (mask & (1 << i)) && tag.isDefined(1 << i) )
		{
			if( tag.isSparse(1 << i) )
			{
				unsigned int count = array_size_send.size();
				array_size_send.push_back(0);
				for(eit = elements[i].begin(); eit != elements[i].end(); eit++)
					if( (*eit)->HaveData(tag) )
					{
						array_size_send.push_back(eit-elements[i].begin());
						array_size_send[count]++;
						size_t s = (*eit)->GetDataSize(tag);
						size_t had_s = array_data_send.size();
						array_data_send.resize(had_s+s*tag.GetBytesSize());
						(*eit)->GetData(tag,0,s,&array_data_send[had_s]);
						if( size == ENUMUNDEF ) array_size_send.push_back(s);						
					}
			}
			else
			{
				for(eit = elements[i].begin(); eit != elements[i].end(); eit++)
				{
					size_t s = (*eit)->GetDataSize(tag);
					size_t had_s = array_data_send.size();
					array_data_send.resize(had_s+s*tag.GetBytesSize());
					(*eit)->GetData(tag,0,s,&array_data_send[had_s]);
					if( size == ENUMUNDEF ) array_size_send.push_back(s);
				}
			}
		}
		array_size_send[0] = array_size_send.size()-2;
		array_size_send[1] = array_data_send.size();
		int buffer_size = 0,position = buffer.size(),temp;
		MPI_Pack_size(array_size_send.size(),INMOST_MPI_DATA_ENUM_TYPE,comm,&temp);
		buffer_size+= temp;
		MPI_Pack_size(array_data_send.size()/tag.GetBytesSize(),tag.GetBulkDataType(),comm,&temp);
		buffer_size+= temp;
		buffer.resize(position+buffer_size);
		MPI_Pack(&array_size_send[0],array_size_send.size(),INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],buffer.size(),&position,comm);
		MPI_Pack(&array_data_send[0],array_data_send.size()/tag.GetBytesSize(),tag.GetBulkDataType(),&buffer[0],buffer.size(),&position,comm);
		buffer.resize(position);		
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << " tag " << tag.GetTagName() << std::endl;
#endif

	}
	
	void Mesh::UnpackTagData(Tag tag, std::vector< std::vector<Element *> > & elements, ElementType mask, std::vector<INMOST_DATA_BULK_TYPE> & buffer, int & position, void (*Operation)(Tag tag, Element * s,INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size))
	{
		if( tag.GetDataType() == DATA_REFERENCE ) return; //NOT IMPLEMENTED
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
		int pos = 0, k = 0;
		INMOST_DATA_ENUM_TYPE data_recv, size_recv;
		std::vector<Element *>::iterator eit;
		unsigned int size = tag.GetSize();
		std::vector<INMOST_DATA_BULK_TYPE> array_data_recv;
		std::vector<INMOST_DATA_ENUM_TYPE> array_size_recv;
		
		MPI_Unpack(&buffer[0],buffer.size(),&position,&size_recv,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
		MPI_Unpack(&buffer[0],buffer.size(),&position,&data_recv,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
		array_size_recv.resize(size_recv);
		array_data_recv.resize(data_recv);
		MPI_Unpack(&buffer[0],buffer.size(),&position,&array_size_recv[0],array_size_recv.size(),INMOST_MPI_DATA_ENUM_TYPE,comm);
		MPI_Unpack(&buffer[0],buffer.size(),&position,&array_data_recv[0],array_data_recv.size()/tag.GetBytesSize(),tag.GetBulkDataType(),comm);
		for(int i = 0; i < 4; i++) if( (mask & (1<<i)) && tag.isDefined(1 << i) )
		{
			if( tag.isSparse(1 << i) )
			{
				unsigned int count = array_size_recv[k++];
				if( size == ENUMUNDEF )
				{
					for(unsigned int j = 0; j < count; j++)
					{
						eit = elements[i].begin() + array_size_recv[k++];
						Operation(tag,*eit,&array_data_recv[pos],array_size_recv[k]);
						pos += array_size_recv[k]*tag.GetBytesSize();
						k++;
					}
				}
				else
				{
					for(unsigned int j = 0; j < count; j++)
					{
						eit = elements[i].begin() + array_size_recv[k++];
						Operation(tag,*eit,&array_data_recv[pos],size);
						pos += size*tag.GetBytesSize();
					}
				}
			}
			else
			{
				if( size == ENUMUNDEF )
				{
					for(eit = elements[i].begin(); eit != elements[i].end(); eit++)
					{
						Operation(tag,*eit,&array_data_recv[pos],array_size_recv[k]);
						pos += array_size_recv[k]*tag.GetBytesSize();
						k++;
					}
				}
				else
				{
					for(eit = elements[i].begin(); eit != elements[i].end(); eit++)
					{
						Operation(tag,*eit,&array_data_recv[pos],size);
						pos += size*tag.GetBytesSize();
					}
				}
			}
		}
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif

	}
	
	
	void Mesh::ExchangeDataInnerBegin(std::vector<Tag> tags, parallel_storage & from, parallel_storage & to, ElementType mask, exchange_data & storage)
	{
		{ // checks for bad input
			if( mask == NONE ) return;
			std::vector<Tag>::iterator it = tags.begin();
			while(it != tags.end() )
				if( !it->isValid() ) it = tags.erase(it);
				else ++it;
			if( tags.empty() ) return;
		}
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
		if( m_state == Serial ) return;
#if defined(USE_MPI)
		//int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber();
		Storage::integer_array procs = IntegerArrayDV(tag_processors);
		Storage::integer_array::iterator p = procs.begin();
		storage.send_buffers.resize(procs.size());
		storage.recv_buffers.resize(procs.size());
		std::vector<int> done;
		parallel_storage::iterator find;
		std::vector<unsigned int> send_size(procs.size(),0), recv_size(procs.size(),0);
		
		bool unknown_size = false;
		for(unsigned int k = 0; k < tags.size(); k++) if( tags[k].GetSize() == ENUMUNDEF ) unknown_size = true;
		
		//precompute sizes
		for(p = procs.begin(); p != procs.end(); p++ )
		{
			int pos = p-procs.begin();
			find = from.find(*p);
			if( find != from.end() )
				for(int i = 0; i < 4; i++)  if( mask & (1 << i) )
					send_size[pos] += find->second[i].size();
			find = to.find(*p);
			if( find != to.end() )
				for(int i = 0; i < 4; i++)  if( mask & (1 << i) )
					recv_size[pos] += find->second[i].size();
		}
		
		int num_send = 0, num_recv = 0;
		for(p = procs.begin(); p != procs.end(); p++ )
		{
			if( send_size[p-procs.begin()] )
			{
				for(unsigned int k = 0; k < tags.size(); k++)
					PackTagData(tags[k],from[*p].get_container(),mask,storage.send_buffers[num_send].second);
				storage.send_buffers[num_send].first = *p;
				num_send++;
			}
			if( recv_size[p-procs.begin()] )
			{
				if( !unknown_size )
				{
					int buffer_size = 0,n = recv_size[p-procs.begin()];
					for(unsigned int k = 0; k < tags.size(); k++)
					{
						int temp;
						MPI_Pack_size(3+n*2,INMOST_MPI_DATA_ENUM_TYPE,comm,&temp); buffer_size += temp;
						MPI_Pack_size(n*tags[k].GetSize(),tags[k].GetBulkDataType(),comm,&temp); buffer_size += temp;
					}
					storage.recv_buffers[num_recv].second.resize(buffer_size);
				}
				storage.recv_buffers[num_recv].first = *p;
				num_recv++;
			}
		}
		storage.send_buffers.resize(num_send);
		storage.recv_buffers.resize(num_recv);
		

		if( unknown_size && parallel_strategy != 0 ) PrepareReceiveInner(UnknownSize,storage.send_buffers,storage.recv_buffers);
		
		ExchangeBuffersInner(storage.send_buffers,storage.recv_buffers,storage.send_reqs,storage.recv_reqs);
#endif

#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}
	
	
	void Mesh::ExchangeDataInnerEnd(std::vector<Tag> tags, parallel_storage & from, parallel_storage & to, ElementType mask, void (*Operation)(Tag tag,Element * element,INMOST_DATA_BULK_TYPE * recv_data,INMOST_DATA_ENUM_TYPE recv_size), exchange_data & storage)
	{
		{ // checks for bad input
			if( mask == NONE ) return;
			std::vector<Tag>::iterator it = tags.begin();
			while(it != tags.end() )
				if( !it->isValid() ) it = tags.erase(it);
				else ++it;
			if( tags.empty() ) return;
		}
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
		if( m_state == Serial ) return;
#if defined(USE_MPI)
		std::vector<int> done;
		while( !(done = FinishRequests(storage.recv_reqs)).empty() )
		{
			for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
			{
				int position = 0;
				for(unsigned int k = 0; k < tags.size(); k++)
					UnpackTagData(tags[k],to[storage.recv_buffers[*qt].first].get_container(),mask,storage.recv_buffers[*qt].second,position,Operation);
			}
		}
		REPORT_MPI(MPI_Waitall(storage.send_reqs.size(),&storage.send_reqs[0],MPI_STATUSES_IGNORE));
#endif

#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}
	

	void Mesh::GatherParallelStorage(parallel_storage & ghost, parallel_storage & shared, ElementType mask)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
#if defined(USE_PARALLEL_STORAGE)
		int mpirank = GetProcessorRank();
		for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); it++)
		{
			Element::Status estat = it->GetStatus();
			if( estat == Element::Shared ) 
			{ 
				Storage::integer_array v = it->IntegerArrayDV(tag_processors);
				for(Storage::integer_array::iterator vit = v.begin(); vit != v.end(); vit++)
					if( *vit != mpirank )
						shared[*vit][ElementNum(it->GetElementType())].push_back(&*it);
			}
			else if( estat == Element::Ghost )
				ghost[it->IntegerDF(tag_owner)][ElementNum(it->GetElementType())].push_back(&*it);
		}
		for(parallel_storage::iterator it = shared.begin(); it != shared.end(); it++)
			for(int i = 0; i < 4; i++) if( mask & (1 << i) )
			{
				qsort(&it->second[i][0],it->second[i].size(),sizeof(Element *),(have_global_id & (1<<i))? CompareElementsCGID : CompareElementsCCentroid);
				REPORT_STR("shared elements size");
				REPORT_STR(ElementTypeName(mask & (1 << i)));
				REPORT_VAL(it->second[i].size());
			}
		for(parallel_storage::iterator it = ghost.begin(); it != ghost.end(); it++)
			for(int i = 0; i < 4; i++) if( mask & (1 << i) )
			{
				qsort(&it->second[i][0],it->second[i].size(),sizeof(Element *),(have_global_id & (1<<i))? CompareElementsCGID : CompareElementsCCentroid);	
				REPORT_STR("ghost elements size");
				REPORT_STR(ElementTypeName(mask & (1 << i)));
				REPORT_VAL(it->second[i].size());
			}
#endif
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}
	
	void Mesh::ExchangeData(Tag tag, ElementType mask)
	{
		ExchangeData(std::vector<Tag>(1,tag),mask);
	}
	
	void Mesh::ExchangeData(std::vector<Tag> tags, ElementType mask)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
		if( m_state == Serial ) return;
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
		exchange_data storage;
		ExchangeDataInnerBegin(tags,shared_elements,ghost_elements,mask,storage);
		ExchangeDataInnerEnd(tags,shared_elements,ghost_elements,mask,DefaultUnpack,storage);
#endif

#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}


	void Mesh::ExchangeDataBegin(Tag tag, ElementType mask, exchange_data & storage)
	{
		ExchangeDataBegin(std::vector<Tag>(1,tag),mask,storage);
	}

	void Mesh::ExchangeDataEnd(Tag tag, ElementType mask, exchange_data & storage)
	{
		ExchangeDataEnd(std::vector<Tag>(1,tag),mask,storage);
	}	
	
	void Mesh::ExchangeDataBegin(std::vector<Tag> tags, ElementType mask, exchange_data & storage)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
		if( m_state == Serial ) return;
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
		ExchangeDataInnerBegin(tags,shared_elements,ghost_elements,mask,storage);
#endif

#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}
	
	void Mesh::ExchangeDataEnd(std::vector<Tag> tags, ElementType mask, exchange_data & storage)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
		if( m_state == Serial ) return;
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
		ExchangeDataInnerEnd(tags,shared_elements,ghost_elements,mask,DefaultUnpack,storage);
#endif

#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}
	
	void Mesh::ReduceData(Tag tag, ElementType mask, void (*Operation)(Tag tag,Element * element,INMOST_DATA_BULK_TYPE * recv_data,INMOST_DATA_ENUM_TYPE recv_size))
	{
		ReduceData(std::vector<Tag>(1,tag),mask,Operation);
	}
	
	void Mesh::ReduceData(std::vector<Tag> tags, ElementType mask, void (*Operation)(Tag tag,Element * element,INMOST_DATA_BULK_TYPE * recv_data,INMOST_DATA_ENUM_TYPE recv_size))
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif

		if( m_state == Serial ) return;
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
		exchange_data storage;
		ExchangeDataInnerBegin(tags,ghost_elements,shared_elements,mask,storage);
		ExchangeDataInnerEnd(tags,ghost_elements,shared_elements,mask,Operation,storage);
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time  << std::endl;
#endif
	}
	
	void Mesh::ReduceDataBegin(Tag tag, ElementType mask, exchange_data & storage)
	{
		ReduceDataBegin(std::vector<Tag>(1,tag),mask,storage);
	}
	
	void Mesh::ReduceDataBegin(std::vector<Tag> tags, ElementType mask, exchange_data & storage)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif

		if( m_state == Serial ) return;
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
		ExchangeDataInnerBegin(tags,ghost_elements,shared_elements,mask,storage);
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time <<  std::endl;
#endif
	}
	
	void Mesh::ReduceDataEnd(Tag tag, ElementType mask, void (*Operation)(Tag tag,Element * element,INMOST_DATA_BULK_TYPE * recv_data,INMOST_DATA_ENUM_TYPE recv_size),exchange_data & storage)
	{
		ReduceDataEnd(std::vector<Tag>(1,tag),mask,Operation,storage);
	}
	
	void Mesh::ReduceDataEnd(std::vector<Tag> tags, ElementType mask, void (*Operation)(Tag tag,Element * element,INMOST_DATA_BULK_TYPE * recv_data,INMOST_DATA_ENUM_TYPE recv_size),exchange_data & storage)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif

		if( m_state == Serial ) return;
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
		ExchangeDataInnerEnd(tags,ghost_elements,shared_elements,mask,Operation,storage);
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}
	
	void Mesh::PackElementsData(std::vector<Element *> & all, std::vector<INMOST_DATA_BULK_TYPE> & buffer, int destination, std::vector<std::string> tag_list)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
		INMOST_DATA_ENUM_TYPE num;
		int temp, position,new_size,k, mpirank = GetProcessorRank();
		std::vector<INMOST_DATA_BULK_TYPE> output;
		std::vector<Element *> scells, snodes, sedges, sfaces;
		std::vector<std::vector<Element *> > pack_tags(4);
		{
			for(std::vector<Element *>::iterator it = all.begin(); it != all.end(); it++)
			{
				switch( (*it)->GetElementType() )
				{
					case CELL: scells.push_back(*it); break;
					case FACE: sfaces.push_back(*it); break;
					case EDGE: sedges.push_back(*it); break;
					case NODE: snodes.push_back(*it); break;
				}
			}
			all.clear();
		}
		{
			adjacent<Element> adj;
			MIDType busy = CreateMarker();
			ArraySetMarker(snodes,busy);
			ArraySetMarker(sedges,busy);
			ArraySetMarker(sfaces,busy);
			ArraySetMarker(scells,busy);
			for(std::vector<Element *>::iterator it = scells.begin(); it != scells.end(); it++)
			{
				adj = (*it)->getAdjElements(NODE);
				for(adjacent<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
					if( !jt->GetMarker(busy) )
					{
						snodes.push_back(&*jt);
						jt->SetMarker(busy);
					}
				adj = (*it)->getAdjElements(EDGE);
				for(adjacent<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
					if( !jt->GetMarker(busy) )
					{
						sedges.push_back(&*jt);
						jt->SetMarker(busy);
					}
				adj = (*it)->getAdjElements(FACE);
				for(adjacent<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
					if( !jt->GetMarker(busy) )
					{
						sfaces.push_back(&*jt);
						jt->SetMarker(busy);
					}
			}
			for(std::vector<Element *>::iterator it = sfaces.begin(); it != sfaces.end(); it++)
			{
				adj = (*it)->getAdjElements(NODE);
				for(adjacent<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
					if( !jt->GetMarker(busy) )
					{
						snodes.push_back(&*jt);
						jt->SetMarker(busy);
					}
				adj = (*it)->getAdjElements(EDGE);
				for(adjacent<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
					if( !jt->GetMarker(busy) )
					{
						sedges.push_back(&*jt);
						jt->SetMarker(busy);
					}
			}
			for(std::vector<Element *>::iterator it = sedges.begin(); it != sedges.end(); it++)
			{
				adj = (*it)->getAdjElements(NODE);
				for(adjacent<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
					if( !jt->GetMarker(busy) )
					{
						snodes.push_back(&*jt);
						jt->SetMarker(busy);
					}
			}
			ArrayRemMarker(snodes,busy);
			ArrayRemMarker(sedges,busy);
			ArrayRemMarker(sfaces,busy);
			ArrayRemMarker(scells,busy);
			ReleaseMarker(busy);
			
			std::sort(snodes.begin(),snodes.end());
			std::sort(sedges.begin(),sedges.end());
			std::sort(sfaces.begin(),sfaces.end());
			std::sort(scells.begin(),scells.end());
		}
		//pack nodes coords
		{
			std::vector<Storage::integer> global_ids;
			position = buffer.size();
			new_size = 0;
			MPI_Pack_size(1								,INMOST_MPI_DATA_ENUM_TYPE		,comm,&temp); new_size += temp;
			MPI_Pack_size(snodes.size()*GetDimensions()	,INMOST_MPI_DATA_REAL_TYPE		,comm,&temp); new_size += temp;
			if( have_global_id & NODE )
				MPI_Pack_size(snodes.size()					,INMOST_MPI_DATA_INTEGER_TYPE		,comm,&temp); new_size += temp;
			buffer.resize(position+new_size);
			num = snodes.size();
			MPI_Pack(&num,1,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],buffer.size(),&position,comm);
			{
				k = 0;
				std::vector<Storage::real> coords(snodes.size()*GetDimensions());
				for(std::vector<Element*>::iterator it = snodes.begin(); it != snodes.end(); it++)
				{
					Storage::real_array c = (*it)->RealArray(tag_coords);
					for(unsigned int j = 0; j < GetDimensions(); j++) coords[k*GetDimensions()+j] = c[j];
					Storage::integer & owner = (*it)->IntegerDF(tag_owner);
					{
						if( owner == mpirank )
						{
							Storage::integer_array proc = (*it)->IntegerArrayDV(tag_processors);
							Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),destination);
							if( ip == proc.end() || (*ip) != destination ) proc.insert(ip,destination);
							(*it)->BulkDF(tag_shared) = Element::Shared;
						}
					}
					if( owner != destination ) 
						pack_tags[0].push_back(*it);
					if( have_global_id & NODE ) global_ids.push_back((*it)->Integer(GlobalIDTag()));
					k++;
				}
				MPI_Pack(&coords[0],snodes.size()*GetDimensions(),INMOST_MPI_DATA_REAL_TYPE,&buffer[0],buffer.size(),&position,comm);
				if(have_global_id & NODE ) 
					MPI_Pack(&global_ids[0],snodes.size(),INMOST_MPI_DATA_INTEGER_TYPE,&buffer[0],buffer.size(),&position,comm);
			}
			buffer.resize(position);
		}
		//pack edges
		{
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size(sedges.size());
			std::vector<Storage::integer> low_conn_nums;
			position = buffer.size();
			new_size = 0;
			num = 0; k = 0;
			for(std::vector<Element *>::iterator it = sedges.begin(); it != sedges.end(); it++) 
			{
				//num += (*it)->low_conn.size();
				low_conn_size[k] = 0;//(*it)->low_conn.size();
				for(Element::adj_iterator jt = (*it)->low_conn.begin(); jt != (*it)->low_conn.end(); jt++) if( !(*jt)->Hidden() )
				{
					std::vector<Element *>::iterator find = std::lower_bound(snodes.begin(),snodes.end(),(*jt));
					if( find == snodes.end() || (*find) != (*jt) ) throw Failure;
					low_conn_nums.push_back(static_cast<Storage::integer>(find - snodes.begin()));
					low_conn_size[k]++;
					num++;
				}
				Storage::integer & owner = (*it)->IntegerDF(tag_owner);
				{
					if( owner == mpirank )
					{
						Storage::integer_array proc = (*it)->IntegerArrayDV(tag_processors);
						Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),destination);
						if( ip == proc.end() || (*ip) != destination ) proc.insert(ip,destination);
						(*it)->BulkDF(tag_shared) = Element::Shared;
					}
				}
				if( owner != destination ) 
					pack_tags[1].push_back(*it);
				k++;
			}
			MPI_Pack_size(1+sedges.size()	,INMOST_MPI_DATA_ENUM_TYPE		,comm,&temp); new_size += temp;
			MPI_Pack_size(num				,INMOST_MPI_DATA_INTEGER_TYPE		,comm,&temp); new_size += temp;
			MPI_Pack_size(sedges.size()		,INMOST_MPI_DATA_INTEGER_TYPE		,comm,&temp); new_size += temp;
			buffer.resize(position+new_size);
			temp = sedges.size();
			MPI_Pack(&temp				,1				,INMOST_MPI_DATA_ENUM_TYPE		,&buffer[0],buffer.size(),&position,comm);
			MPI_Pack(&low_conn_size[0]	,sedges.size()	,INMOST_MPI_DATA_ENUM_TYPE		,&buffer[0],buffer.size(),&position,comm);
			MPI_Pack(&low_conn_nums[0]	,num			,INMOST_MPI_DATA_INTEGER_TYPE		,&buffer[0],buffer.size(),&position,comm);
			buffer.resize(position);
		}
		//pack faces
		{
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size(sfaces.size());
			std::vector<Storage::integer> low_conn_nums;
			position = buffer.size();
			new_size = 0;
			num = 0; k = 0;
			for(std::vector<Element *>::iterator it = sfaces.begin(); it != sfaces.end(); it++) 
			{
				//num += (*it)->low_conn.size();
				low_conn_size[k] = 0;//(*it)->low_conn.size();
				for(Element::adj_iterator jt = (*it)->low_conn.begin(); jt != (*it)->low_conn.end(); jt++) if( !(*jt)->Hidden() )
				{
					std::vector<Element *>::iterator find = std::lower_bound(sedges.begin(),sedges.end(),(*jt));
					if( find == sedges.end() || (*find) != (*jt) ) throw Failure;
					low_conn_nums.push_back(static_cast<Storage::integer>(find - sedges.begin()));
					low_conn_size[k]++;
					num++;
				}
				Storage::integer & owner = (*it)->IntegerDF(tag_owner);
				{
					if( owner == mpirank )
					{
						Storage::integer_array proc = (*it)->IntegerArrayDV(tag_processors);
						Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),destination);
						if( ip == proc.end() || (*ip) != destination ) proc.insert(ip,destination);
						(*it)->BulkDF(tag_shared) = Element::Shared;
					}
				}
				if( owner != destination ) 
					pack_tags[2].push_back(*it);
				k++;
			}
			MPI_Pack_size(1+sfaces.size()	,INMOST_MPI_DATA_ENUM_TYPE		,comm,&temp); new_size += temp;
			MPI_Pack_size(num				,INMOST_MPI_DATA_INTEGER_TYPE		,comm,&temp); new_size += temp;
			MPI_Pack_size(sfaces.size()		,INMOST_MPI_DATA_INTEGER_TYPE		,comm,&temp); new_size += temp;
			buffer.resize(position+new_size);
			temp = sfaces.size();
			MPI_Pack(&temp				,1				,INMOST_MPI_DATA_ENUM_TYPE		,&buffer[0],buffer.size(),&position,comm);
			MPI_Pack(&low_conn_size[0]	,sfaces.size()	,INMOST_MPI_DATA_ENUM_TYPE		,&buffer[0],buffer.size(),&position,comm);
			MPI_Pack(&low_conn_nums[0]	,num			,INMOST_MPI_DATA_INTEGER_TYPE		,&buffer[0],buffer.size(),&position,comm);
			buffer.resize(position);
		}
		//pack cells
		{
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size(scells.size()), high_conn_size(scells.size());
			std::vector<Storage::integer> low_conn_nums;
			std::vector<Storage::integer> high_conn_nums;
			INMOST_DATA_ENUM_TYPE num_high = 0;
			position = buffer.size();
			new_size = 0;
			num = 0; k = 0;
			for(std::vector<Element *>::iterator it = scells.begin(); it != scells.end(); it++) 
			{
				//num += (*it)->low_conn.size();
				low_conn_size[k] = 0;//(*it)->low_conn.size();				
				for(Element::adj_iterator jt = (*it)->low_conn.begin(); jt != (*it)->low_conn.end(); jt++) if( !(*jt)->Hidden() )
				{
					std::vector<Element *>::iterator find = std::lower_bound(sfaces.begin(),sfaces.end(),(*jt));
					if( find == sfaces.end() || (*find) != (*jt) ) throw Failure;
					low_conn_nums.push_back(static_cast<Storage::integer>(find - sfaces.begin()));
					low_conn_size[k]++;
					num++;
				}
				//num_high += (*it)->high_conn.size();
				high_conn_size[k] = 0;//(*it)->high_conn.size();
				for(Element::adj_iterator jt = (*it)->high_conn.begin(); jt != (*it)->high_conn.end(); jt++) if( !(*jt)->Hidden() )
				{
					std::vector<Element *>::iterator find = std::lower_bound(snodes.begin(),snodes.end(),(*jt));
					if( find == snodes.end() || (*find) != (*jt) ) throw Failure;
					high_conn_nums.push_back(static_cast<Storage::integer>(find - snodes.begin()));
					high_conn_size[k]++;
					num_high++;
				}
				Storage::integer & owner = (*it)->IntegerDF(tag_owner);
				{
					if( owner == mpirank )
					{
						Storage::integer_array proc = (*it)->IntegerArrayDV(tag_processors);
						Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),destination);
						if( ip == proc.end() || (*ip) != destination ) proc.insert(ip,destination);
						(*it)->BulkDF(tag_shared) = Element::Shared;
					}
				}
				if( owner != destination ) 
					pack_tags[3].push_back(*it);
				k++;
			}
			MPI_Pack_size(1					,INMOST_MPI_DATA_ENUM_TYPE	,comm,&temp); new_size += temp;
			MPI_Pack_size(scells.size()		,INMOST_MPI_DATA_ENUM_TYPE	,comm,&temp); new_size += temp;
			MPI_Pack_size(num				,INMOST_MPI_DATA_INTEGER_TYPE	,comm,&temp); new_size += temp;
			MPI_Pack_size(scells.size()		,INMOST_MPI_DATA_ENUM_TYPE	,comm,&temp); new_size += temp;
			MPI_Pack_size(num_high			,INMOST_MPI_DATA_INTEGER_TYPE	,comm,&temp); new_size += temp;
			MPI_Pack_size(scells.size()		,INMOST_MPI_DATA_INTEGER_TYPE	,comm,&temp); new_size += temp;
			buffer.resize(position+new_size);
			temp = scells.size();
			MPI_Pack(&temp				,1				,INMOST_MPI_DATA_ENUM_TYPE	,&buffer[0],buffer.size(),&position,comm);
			MPI_Pack(&low_conn_size[0]	,scells.size()	,INMOST_MPI_DATA_ENUM_TYPE	,&buffer[0],buffer.size(),&position,comm);
			MPI_Pack(&low_conn_nums[0]	,num			,INMOST_MPI_DATA_INTEGER_TYPE	,&buffer[0],buffer.size(),&position,comm);
			MPI_Pack(&high_conn_size[0]	,scells.size()	,INMOST_MPI_DATA_ENUM_TYPE	,&buffer[0],buffer.size(),&position,comm);
			MPI_Pack(&high_conn_nums[0]	,num_high		,INMOST_MPI_DATA_INTEGER_TYPE	,&buffer[0],buffer.size(),&position,comm);
			buffer.resize(position);
		}
		all.insert(all.end(),scells.begin(),scells.end());
		all.insert(all.end(),sfaces.begin(),sfaces.end());
		all.insert(all.end(),sedges.begin(),sedges.end());
		all.insert(all.end(),snodes.begin(),snodes.end());
		Storage::integer_array proc = IntegerArrayDV(tag_processors);
		Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),destination);
		if( ip == proc.end() || (*ip) != destination ) proc.insert(ip,destination);		
		//pack number of tags
		{
			position = buffer.size();
			num = tag_list.size();
			MPI_Pack_size(1,INMOST_MPI_DATA_ENUM_TYPE,comm,&new_size);
			buffer.resize(position+new_size);
			MPI_Pack(&num,1,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],buffer.size(),&position,comm);
			buffer.resize(position);
		}
				
		for(INMOST_DATA_ENUM_TYPE i = 0; i < tag_list.size(); i++)
		{
			{
				//pack tag name
				{
					Tag tag = GetTag(tag_list[i]);
					INMOST_DATA_ENUM_TYPE defined, sparse, datatype, datasize;
					datatype = static_cast<INMOST_DATA_ENUM_TYPE>(tag.GetDataType());
					datasize = tag.GetSize();
					defined = NONE;
					sparse = NONE;
					for(ElementType mask = NODE; mask <= MESH; mask = mask << 1)
					{
						if( tag.isDefined(mask) ) defined |= mask;
						if( tag.isSparse(mask) ) sparse |= mask;
					}
					
					
					
					position = buffer.size();
					new_size = 0;
					num = tag_list[i].size();
					MPI_Pack_size(5		,INMOST_MPI_DATA_ENUM_TYPE,comm,&temp); new_size += temp;
					MPI_Pack_size(num	,MPI_CHAR				,comm,&temp); new_size += temp;
					buffer.resize(position+new_size);
					MPI_Pack(&datatype			,1	,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],buffer.size(),&position,comm);
					MPI_Pack(&defined			,1	,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],buffer.size(),&position,comm);
					MPI_Pack(&sparse			,1	,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],buffer.size(),&position,comm);
					MPI_Pack(&datasize			,1	,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],buffer.size(),&position,comm);
					MPI_Pack(&num				,1	,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],buffer.size(),&position,comm);
					MPI_Pack(&tag_list[i][0]	,num,MPI_CHAR				,&buffer[0],buffer.size(),&position,comm);
					
					
					
					buffer.resize(position);
				}
				PackTagData(GetTag(tag_list[i]),pack_tags,NODE | EDGE | FACE | CELL | ESET,buffer);
			}
		}		
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << " for " << destination << std::endl;
#endif
	}
	
	
	void Mesh::UnpackElementsData(std::vector<Element *> & all, std::vector<INMOST_DATA_BULK_TYPE> & buffer, int source, std::vector<std::string> & tag_recv)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
		int mpirank = GetProcessorRank();
		INMOST_DATA_ENUM_TYPE num, temp;
		INMOST_DATA_ENUM_TYPE shift = 0;
		int position = 0;
		std::vector<Element *> scells, snodes, sedges, sfaces;
		std::vector<std::vector<Element *> > unpack_tags(4);
		std::vector<Node *> old_nodes;
		all.clear();
		for(nodes_container::iterator it = nodes.begin(); it != nodes.end(); it++) if( *it != NULL )
			old_nodes.push_back(*it);		
		
		if( have_global_id & NODE )
			qsort(&old_nodes[0],old_nodes.size(),sizeof(Node *),CompareElementsCGID);
		else
			qsort(&old_nodes[0],old_nodes.size(),sizeof(Node *),CompareElementsCCentroid);
		//sort(&old_nodes[0],old_nodes.size(),sizeof(Node *),CompareElement,SwapPointer,NULL);
		//std::sort(old_nodes.begin(),old_nodes.end(),LessElementsCentroid);
		//unpack nodes
		{
			unsigned int dim = GetDimensions();
			std::vector<Storage::real> coords;
			MPI_Unpack(&buffer[0],buffer.size(),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			coords.resize(num*dim);
			MPI_Unpack(&buffer[0],buffer.size(),&position,&coords[0],num*dim,INMOST_MPI_DATA_REAL_TYPE,comm);
			
			std::vector<Storage::integer> global_ids;
			if( have_global_id & NODE )
			{
				global_ids.resize(num);
				MPI_Unpack(&buffer[0],buffer.size(),&position,&global_ids[0],num,INMOST_MPI_DATA_INTEGER_TYPE,comm);
			}
			
			
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				Node * new_node;
				int find;
				if( have_global_id & NODE )
					find = binary_search(&old_nodes[0],old_nodes.size(),sizeof(Node*),CompareGIDSearch,&global_ids[i]);
				else
					find = binary_search(&old_nodes[0],old_nodes.size(),sizeof(Node*),CompareCoordSearch,&coords[i*dim]);
				if( find != -1 ) 
				{
					new_node = old_nodes[find];
					if( new_node->IntegerDF(tag_owner) != mpirank ) 
						unpack_tags[0].push_back(static_cast<Element *>(new_node));
				}
				else
				{
					new_node = CreateNode(&coords[i*dim]);
					{
						new_node->BulkDF(tag_shared) = Element::Ghost;
						new_node->IntegerDF(tag_owner) = -1;
						unpack_tags[0].push_back(static_cast<Element *>(new_node));
					}
				}
				snodes.push_back(static_cast<Element *>(new_node));
			}
		}
		//unpack edges
		{
			std::vector<Node *> e_nodes;
			shift = 0;
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size;
			std::vector<Storage::integer> low_conn_nums;
			MPI_Unpack(&buffer[0],buffer.size(),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			low_conn_size.resize(num);
			MPI_Unpack(&buffer[0],buffer.size(),&position,&low_conn_size[0],num,INMOST_MPI_DATA_ENUM_TYPE,comm);
			temp = 0;
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += low_conn_size[i];
			low_conn_nums.resize(temp);
			MPI_Unpack(&buffer[0],buffer.size(),&position,&low_conn_nums[0],temp,INMOST_MPI_DATA_INTEGER_TYPE,comm);
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				e_nodes.resize(low_conn_size[i]);
				for(INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
					e_nodes[j] = static_cast<Node *>(snodes[low_conn_nums[shift+j]]);
				Edge * new_edge = static_cast<Edge *>(FindSharedAdjacency(reinterpret_cast<Element **>(&e_nodes[0]),e_nodes.size()));
				if( new_edge == NULL )
				{
					new_edge = CreateEdge(e_nodes);
					{
						new_edge->BulkDF(tag_shared) = Element::Ghost;
						new_edge->IntegerDF(tag_owner) = -1;
						unpack_tags[1].push_back(static_cast<Element *>(new_edge));
					}
				}
				else if( new_edge->IntegerDF(tag_owner) != mpirank ) unpack_tags[1].push_back(static_cast<Element *>(new_edge));
				sedges.push_back(static_cast<Element *>(new_edge));
				shift+= low_conn_size[i];
			}
		}
		//unpack faces
		{
			std::vector<Edge *> f_edges;
			shift = 0;
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size;
			std::vector<Storage::integer> low_conn_nums;
			MPI_Unpack(&buffer[0],buffer.size(),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			low_conn_size.resize(num);
			MPI_Unpack(&buffer[0],buffer.size(),&position,&low_conn_size[0],num,INMOST_MPI_DATA_ENUM_TYPE,comm);
			temp = 0;
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += low_conn_size[i];
			low_conn_nums.resize(temp);
			MPI_Unpack(&buffer[0],buffer.size(),&position,&low_conn_nums[0],temp,INMOST_MPI_DATA_INTEGER_TYPE,comm);
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				f_edges.resize(low_conn_size[i]);
				for(INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
					f_edges[j] = static_cast<Edge *>(sedges[low_conn_nums[shift+j]]);
				Face * new_face = static_cast<Face *>(FindSharedAdjacency(reinterpret_cast<Element **>(&f_edges[0]),f_edges.size()));
				if( new_face == NULL )
				{
					new_face = CreateFace(f_edges);
					{
						new_face->BulkDF(tag_shared) = Element::Ghost;
						new_face->IntegerDF(tag_owner) = -1;
						unpack_tags[2].push_back(static_cast<Element *>(new_face));
					}
				} else if( new_face->IntegerDF(tag_owner) != mpirank ) unpack_tags[2].push_back(static_cast<Element *>(new_face));
				sfaces.push_back(static_cast<Element *>(new_face));
				shift+= low_conn_size[i];
			}
		}
		//unpack cells
		{
			std::vector<Face *> c_faces;
			std::vector<Node *> c_nodes;
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size;
			std::vector<INMOST_DATA_ENUM_TYPE> high_conn_size;
			std::vector<Storage::integer> low_conn_nums;
			std::vector<Storage::integer> high_conn_nums;
			INMOST_DATA_ENUM_TYPE shift_high = 0;
			shift = 0;
			MPI_Unpack(&buffer[0],buffer.size(),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			{
				low_conn_size.resize(num);
				MPI_Unpack(&buffer[0],buffer.size(),&position,&low_conn_size[0],num,INMOST_MPI_DATA_ENUM_TYPE,comm);
				temp = 0;
				for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += low_conn_size[i];
				low_conn_nums.resize(temp);
				MPI_Unpack(&buffer[0],buffer.size(),&position,&low_conn_nums[0],temp,INMOST_MPI_DATA_INTEGER_TYPE,comm);
			}
			{
				high_conn_size.resize(num);
				MPI_Unpack(&buffer[0],buffer.size(),&position,&high_conn_size[0],num,INMOST_MPI_DATA_ENUM_TYPE,comm);
				temp = 0;
				for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += high_conn_size[i];
				high_conn_nums.resize(temp);
				MPI_Unpack(&buffer[0],buffer.size(),&position,&high_conn_nums[0],temp,INMOST_MPI_DATA_INTEGER_TYPE,comm);
			}
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				c_faces.resize(low_conn_size[i]);
				for(INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
					c_faces[j] = static_cast<Face *>(sfaces[low_conn_nums[shift+j]]);
				c_nodes.resize(high_conn_size[i]);
				for(INMOST_DATA_ENUM_TYPE j = 0; j < high_conn_size[i]; j++)
					c_nodes[j] = static_cast<Node *>(snodes[high_conn_nums[shift_high+j]]);
				Cell * new_cell = static_cast<Cell *>(FindSharedAdjacency(reinterpret_cast<Element **>(&c_faces[0]),c_faces.size()));
				if( new_cell == NULL )
				{
					new_cell = CreateCell(c_faces,c_nodes);
					{
						new_cell->BulkDF(tag_shared) = Element::Ghost;
						new_cell->IntegerDF(tag_owner) = -1;
						unpack_tags[3].push_back(static_cast<Element *>(new_cell));
					}
				} else if( new_cell->IntegerDF(tag_owner) != mpirank ) unpack_tags[3].push_back(static_cast<Element *>(new_cell));
				scells.push_back(static_cast<Element *>(new_cell));
				shift += low_conn_size[i];
				shift_high += high_conn_size[i];
			}
		}
		all.insert(all.end(),scells.begin(),scells.end());
		all.insert(all.end(),sfaces.begin(),sfaces.end());
		all.insert(all.end(),sedges.begin(),sedges.end());
		all.insert(all.end(),snodes.begin(),snodes.end());
		Storage::integer_array proc = IntegerArrayDV(tag_processors);
		Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),source);
		if( ip == proc.end() || (*ip) != source ) proc.insert(ip,source);
		//unpack tags
		{	
			MPI_Unpack(&buffer[0],buffer.size(),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				INMOST_DATA_ENUM_TYPE defined, sparse, datatype, datasize;
				std::string tag_name;
				MPI_Unpack(&buffer[0],buffer.size(),&position,&datatype,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
				MPI_Unpack(&buffer[0],buffer.size(),&position,&defined,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
				MPI_Unpack(&buffer[0],buffer.size(),&position,&sparse,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
				MPI_Unpack(&buffer[0],buffer.size(),&position,&datasize,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
				MPI_Unpack(&buffer[0],buffer.size(),&position,&temp,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
				tag_name.resize(temp);
				MPI_Unpack(&buffer[0],buffer.size(),&position,&tag_name[0],temp,MPI_CHAR,comm);
				tag_recv.push_back(tag_name);
				Tag tag = CreateTag(tag_name,static_cast<enum DataType>(datatype),static_cast<ElementType>(defined),static_cast<ElementType>(sparse),datasize);
				UnpackTagData(tag,unpack_tags,NODE | EDGE | FACE | CELL | ESET, buffer,position,DefaultUnpack);
			}
		}
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << " from " << source << std::endl;
#endif
	}
	
	
	void Mesh::ExchangeBuffersInner(std::vector< std::pair<int, std::vector<INMOST_DATA_BULK_TYPE> > > & send_bufs, 
									std::vector< std::pair<int, std::vector<INMOST_DATA_BULK_TYPE> > > & recv_bufs,
									std::vector<INMOST_MPI_Request> & send_reqs, std::vector<INMOST_MPI_Request> & recv_reqs)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
		unsigned i;
		int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber(), rand_num = randomizer.Number()+1;
		int mpi_tag;
		int max_tag, flag;
		void * p_max_tag;
		MPI_Comm_get_attr(comm,MPI_TAG_UB,&p_max_tag,&flag);
		max_tag = *static_cast<int *>(p_max_tag);
		if( !flag ) throw Impossible;
		recv_reqs.resize(recv_bufs.size());
		send_reqs.resize(send_bufs.size());
		if( parallel_strategy == 0 )
		{
			for(i = 0; i < send_bufs.size(); i++)
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (send_bufs[i].first+mpisize+rand_num))%max_tag;
				REPORT_MPI(MPI_Isend(&send_bufs[i].second[0],send_bufs[i].second.size(),MPI_PACKED,send_bufs[i].first,mpi_tag,comm,&send_reqs[i]));	
			}	
			for(i = 0; i < recv_bufs.size(); i++)
			{
				int buffer_size;
				MPI_Status recv_stat;
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (mpirank+mpisize+rand_num))%max_tag;
				REPORT_MPI(MPI_Probe(MPI_ANY_SOURCE,mpi_tag,comm,&recv_stat));
				REPORT_MPI(MPI_Get_count(&recv_stat,MPI_PACKED,&buffer_size));
				recv_bufs[i].second.resize(buffer_size);
				recv_bufs[i].first = recv_stat.MPI_SOURCE;
				REPORT_MPI(MPI_Irecv(&recv_bufs[i].second[0],recv_bufs[i].second.size(),MPI_PACKED,recv_stat.MPI_SOURCE,mpi_tag,comm,&recv_reqs[i]));
			}
		}
		else if( parallel_strategy == 1 )
		{
			for(i = 0; i < recv_bufs.size(); i++) if( !recv_bufs[i].second.empty() )
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (mpirank+mpisize+rand_num))%max_tag;
				//mpi_tag = parallel_mesh_unique_id*mpisize*mpisize+recv_bufs[i].first*mpisize+mpirank;
				REPORT_MPI(MPI_Irecv(&recv_bufs[i].second[0],recv_bufs[i].second.size(),MPI_PACKED,recv_bufs[i].first,mpi_tag,comm,&recv_reqs[i]));
			}
			for(i = 0; i < send_bufs.size(); i++) if( !send_bufs[i].second.empty() )
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (send_bufs[i].first+mpisize+rand_num))%max_tag;
				//mpi_tag = parallel_mesh_unique_id*mpisize*mpisize+mpirank*mpisize+send_bufs[i].first;
				REPORT_MPI(MPI_Isend(&send_bufs[i].second[0],send_bufs[i].second.size(),MPI_PACKED,send_bufs[i].first,mpi_tag,comm,&send_reqs[i]));	
			}
		}
		else if( parallel_strategy == 2 )
		{
			for(i = 0; i < recv_bufs.size(); i++) if( !recv_bufs[i].second.empty() )
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (mpirank+mpisize+rand_num))%max_tag;
				REPORT_MPI(MPI_Irecv(&recv_bufs[i].second[0],recv_bufs[i].second.size(),MPI_PACKED,recv_bufs[i].first,mpi_tag,comm,&recv_reqs[i]));
			}
			REPORT_MPI(MPI_Barrier(comm));
			for(i = 0; i < send_bufs.size(); i++) if( !send_bufs[i].second.empty() )
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (send_bufs[i].first+mpisize+rand_num))%max_tag;
				REPORT_MPI(MPI_Irsend(&send_bufs[i].second[0],send_bufs[i].second.size(),MPI_PACKED,send_bufs[i].first,mpi_tag,comm,&send_reqs[i]));	
			}
		}
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}
	
	void Mesh::PrepareReceiveInner(Prepare todo, std::vector< std::pair<int, std::vector<INMOST_DATA_BULK_TYPE> > > & send_bufs, std::vector< std::pair<int, std::vector<INMOST_DATA_BULK_TYPE> > > & recv_bufs)
	{
		if( parallel_strategy == 0 && todo == UnknownSize ) return; //in this case we know all we need
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
#if defined(USE_MPI2)
		int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber();
		unsigned i, end = send_bufs.size();
		memset(shared_space,0,sizeof(unsigned)*mpisize); //zero bits where we receive data
		REPORT_MPI(MPI_Win_fence(MPI_MODE_NOPRECEDE,window)); //start exchange session
		for(i = 0; i < end; i++) shared_space[mpisize+i] = send_bufs[i].second.size()+1; //put data to special part of the memory
		for(i = 0; i < end; i++) REPORT_MPI(MPI_Put(&shared_space[mpisize+i],1,MPI_UNSIGNED,send_bufs[i].first,mpirank,1,MPI_UNSIGNED,window)); //request rdma
		REPORT_MPI(MPI_Win_fence(MPI_MODE_NOSUCCEED,window)); //end exchange session
		if( parallel_strategy == 0 )
		{
			unsigned num = 0;
			for(int ii = 0; ii < mpisize; ii++) if( shared_space[ii] > 0 ) num++;
			recv_bufs.resize(num);
		}
		else if( todo == UnknownSize )
		{
			end = recv_bufs.size();
			for(i = 0; i < end; i++)
				recv_bufs[i].second.resize(shared_space[recv_bufs[i].first]-1);
		}
		else if( todo == UnknownSource )
		{
			recv_bufs.clear();
			for(int ii = 0; ii < mpisize; ii++)
				if( shared_space[ii] > 0 )
					recv_bufs.push_back(std::pair<int,std::vector<INMOST_DATA_BULK_TYPE> >(ii,std::vector<INMOST_DATA_BULK_TYPE>(shared_space[ii]-1))); // this call would be optimized by compiler
		}
#else
		if( todo == UnknownSize )
		{
			int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber(), rand_num = randomizer.Number()+1;
			int mpi_tag;
			int max_tag, flag;
			void * p_max_tag;
			MPI_Comm_get_attr(comm,MPI_TAG_UB,&p_max_tag,&flag);
			max_tag = *static_cast<int *>(p_max_tag);
			std::vector<int> send_recv_size(send_bufs.size()+recv_bufs.size());
			std::vector<INMOST_MPI_Request> reqs(send_bufs.size()+recv_bufs.size());
			for(i = 0; i < send_bufs.size(); i++) send_recv_size[i+recv_bufs.size()] = send_bufs[i].second.size();

			for(i = 0; i < recv_bufs.size(); i++)
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (mpirank+mpisize+rand_num))%max_tag;
				//mpi_tag = parallel_mesh_unique_id*mpisize*mpisize+recv_bufs[i].first*mpisize+mpirank;
				REPORT_MPI(MPI_Irecv(&send_recv_size[i],1,MPI_INT,recv_bufs[i].first,mpi_tag,comm,&reqs[i]));
			}
			for(i = 0; i < send_bufs.size(); i++)
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (send_bufs[i].first+mpisize+rand_num))%max_tag;
				//mpi_tag = parallel_mesh_unique_id*mpisize*mpisize+mpirank*mpisize+send_bufs[i].first;
				REPORT_MPI(MPI_Isend(&send_recv_size[i+recv_bufs.size()],1,MPI_INT,send_bufs[i].first,mpi_tag,comm,&reqs[i+recv_bufs.size()]));	
			}
			REPORT_MPI(MPI_Waitall(recv_bufs.size(),&reqs[0],MPI_STATUSES_IGNORE));
			for(i = 0; i < recv_bufs.size(); i++) recv_bufs[i].second.resize(send_recv_size[i]);
			REPORT_MPI(MPI_Waitall(send_bufs.size(),&reqs[recv_bufs.size()],MPI_STATUSES_IGNORE));
		}
		else if( todo == UnknownSource )
		{
			int mpisize = GetProcessorsNumber(),mpirank = GetProcessorRank();
			
			if( parallel_strategy == 0 )
			{
				std::vector<unsigned> recvs_at_proc_temp(mpisize,0);
				std::vector<unsigned> recvs_at_proc(mpisize,0);
				for(unsigned i = 0; i < send_bufs.size(); i++)
					recvs_at_proc_temp[send_bufs[i].first]++;
				REPORT_MPI(MPI_Allreduce(&recvs_at_proc_temp[0],&recvs_at_proc[0],mpisize,MPI_UNSIGNED,MPI_SUM,comm));
				recv_bufs.resize(recvs_at_proc[mpirank]);
			}
			else if( parallel_strategy == 1 || parallel_strategy == 2 )
			{
				std::vector< unsigned > sends_dest_and_size;
				sends_dest_and_size.reserve(send_bufs.size()*2);
				std::vector< std::pair<int, std::vector<INMOST_DATA_BULK_TYPE> > >::iterator kt;
				for(kt = send_bufs.begin(); kt != send_bufs.end(); kt++)
				{
					sends_dest_and_size.push_back(kt->first);
					sends_dest_and_size.push_back(kt->second.size());
				}
				unsigned recvsize = 0;
				int k,j;
				std::vector<int> allsize(mpisize);
				int size = sends_dest_and_size.size();
				REPORT_MPI(MPI_Allgather(&size,1,MPI_INT,&allsize[0],1,MPI_INT,comm));
				std::vector<int> displs(mpisize+1,0);
				for(k = 0; k < mpisize; k++)
					recvsize += allsize[k];
				for(k = 1; k < mpisize+1; k++)
					displs[k] = displs[k-1]+allsize[k-1];
				std::vector<unsigned> recvs_dest_and_size(recvsize);
				REPORT_MPI(MPI_Allgatherv(&sends_dest_and_size[0],sends_dest_and_size.size(),MPI_UNSIGNED,&recvs_dest_and_size[0],&allsize[0],&displs[0],MPI_UNSIGNED,comm));
				recv_bufs.clear();
				for(k = 0; k < mpisize; k++)
					for(j = displs[k]; j < displs[k+1]; j+=2)
					{
						//if( mpirank == 0 ) std::cout << "proc " << k << " sends " << recvs_dest_and_size[j] << " size " << recvs_dest_and_size[j+1] << std::endl;
						if(  static_cast<int>(recvs_dest_and_size[j]) == mpirank )
						{
							recv_bufs.push_back(std::pair<int, std::vector<INMOST_DATA_BULK_TYPE> >(k,std::vector<INMOST_DATA_BULK_TYPE>()));
							recv_bufs.back().second.resize(recvs_dest_and_size[j+1]);
						}
					}

			}
		}
#endif
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}
	
	std::vector<int> Mesh::FinishRequests(std::vector<INMOST_MPI_Request> & recv_reqs)
	{
        std::vector<int> ret(recv_reqs.size(),-1);
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
		int outcount = 0;
		REPORT_VAL(recv_reqs.size());
		if( recv_reqs.empty() ) 
			ret.clear();
		else
		{
			REPORT_MPI(MPI_Waitsome(recv_reqs.size(),&recv_reqs[0],&outcount,&ret[0],MPI_STATUSES_IGNORE));
			if( outcount == MPI_UNDEFINED ) ret.clear();
			else ret.resize(outcount);
		}
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
        return ret;
	}
	
	
	void Mesh::ExchangeMarked(enum Action action)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
		if( m_state == Serial ) return;
#if defined(USE_MPI)
		INMOST_DATA_BIG_ENUM_TYPE num_wait;
		int mpirank = GetProcessorRank();
		std::vector<MPI_Request> send_reqs, recv_reqs;
		std::vector<std::string> tag_list, tag_list_recv;
		std::vector<std::string> tag_list_empty;
		std::vector< std::pair<int, std::vector<INMOST_DATA_BULK_TYPE> > > send_bufs;
		std::vector< std::pair<int, std::vector<INMOST_DATA_BULK_TYPE> > > recv_bufs;
		std::map< int , std::vector<Element *> > send_elements;
		std::vector<int> done;
		
		ListTagNames(tag_list);
		{
			std::vector<std::string>::iterator it = tag_list.begin();
			while( it != tag_list.end() )
			{
				if( it->substr(0,9) == "PROTECTED" ) 
					it = tag_list.erase(it);
				else it++;
			}
		}
		{
			for(Mesh::iteratorElement it = BeginElement(CELL | FACE | EDGE | NODE); it != EndElement(); ++it) if( it->HaveData(tag_sendto) )
			{
				Storage::integer_array mark = it->IntegerArray(tag_sendto);
				for(Storage::integer_array::iterator kt = mark.begin(); kt != mark.end(); kt++)
					send_elements[*kt].push_back(&*it);
				it->DelData(tag_sendto);
			}
		}
		num_wait = 0;
		send_bufs.resize(send_elements.size());
		for(std::map< int , std::vector<Element *> >::iterator it = send_elements.begin(); it != send_elements.end(); it++)
			if( !it->second.empty() )
			{
				PackElementsData(it->second,send_bufs[num_wait].second,it->first,tag_list);
				send_bufs[num_wait].first = it->first;
				num_wait++;
			}
		send_bufs.resize(num_wait);	

		//NEW ALGORITHM
		if( false )
		{ //in most cases the exchange of elements will effect only nearest neighbours, we want to detect this
			Storage::integer_array p = IntegerArrayDV(tag_processors);
			int halo_exchange = 0, test = 1;
			for(unsigned k = 0; k < send_bufs.size(); k++)
				if( !std::binary_search(p.begin(),p.end(),send_bufs[k].first) ) test = 0;
			//we should exchange, because other processor may want to send something to us
			REPORT_MPI(MPI_Allreduce(&test,&halo_exchange,1,MPI_INT,MPI_MIN,GetCommunicator())); 
			if( halo_exchange )
			{
				//prepare missing send buffers
				{
					std::vector<int> present(send_bufs.size()), missing(p.size());
					for(unsigned k = 0; k < send_bufs.size(); k++) present[k] = send_bufs[k].first;
					missing.resize(std::set_difference(p.begin(),p.end(),present.begin(),present.end(),missing.begin())-missing.begin());
					for(unsigned k = 0; k < missing.size(); k++) send_bufs.push_back(std::pair<int,std::vector<INMOST_DATA_BULK_TYPE> >(missing[k],std::vector<INMOST_DATA_BULK_TYPE>()));
				}
				//prepare recv buffers
				{
					recv_bufs.resize(p.size());
					for(unsigned k = 0; k < p.size(); k++) recv_bufs[k].first = p[k];
				}
				PrepareReceiveInner(UnknownSize, send_bufs,recv_bufs);
			}
			else PrepareReceiveInner(UnknownSource, send_bufs,recv_bufs);
		}
		else PrepareReceiveInner(UnknownSource, send_bufs,recv_bufs);
		ExchangeBuffersInner(send_bufs,recv_bufs,send_reqs,recv_reqs);
			
		if( action == AMigrate ) // delete packed elements
		{
			Tag tag_new_processors = GetTag("TEMPORARY_NEW_PROCESSORS");
			{
				std::vector< std::vector<Element *> > delete_elements(4);
				MIDType delete_marker = CreateMarker();
				
				for(std::map< int , std::vector<Element *> >::iterator it = send_elements.begin(); it != send_elements.end(); it++)
					if( !it->second.empty() )
					{
						//Have packed all entities, now delete those entities that i should not own	
						for(std::vector<Element *>::iterator kt = it->second.begin(); kt != it->second.end(); kt++)
							if( !(*kt)->GetMarker(delete_marker) )
							{
								Storage::integer_array procs = (*kt)->IntegerArray(tag_new_processors);
								if( !std::binary_search(procs.begin(),procs.end(),mpirank) ) 
								{
									delete_elements[ElementNum((*kt)->GetElementType())].push_back(*kt);
									(*kt)->SetMarker(delete_marker);
								}
							}
					}
				
				
				//Should delete elements in order by CELL..NODE
				for(int j = 3; j >= 0; j--)
				{
					ArrayRemMarker(delete_elements[j],delete_marker);
					for(std::vector<Element *>::iterator kt = delete_elements[j].begin(); kt != delete_elements[j].end(); kt++)
						delete (*kt);
				}
				
				ReleaseMarker(delete_marker);
			}
			
			//now delete elements that we have not yet deleted
			for(iteratorElement it = BeginElement(CELL | FACE | EDGE | NODE); it != EndElement(); it++)
			{
				Storage::integer_array procs = it->IntegerArray(tag_new_processors);
				if( !std::binary_search(procs.begin(),procs.end(),mpirank) ) delete(&*it);
			}
		}
		

		send_elements.clear();
		
		while( !(done = FinishRequests(recv_reqs)).empty() )
		{
			for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
			{
				std::vector<Element *> recv_elements;
				UnpackElementsData(recv_elements,recv_bufs[*qt].second,recv_bufs[*qt].first,tag_list_recv);
				if( action == AGhost )
				{
					for(std::vector<Element *>::iterator it = recv_elements.begin(); it != recv_elements.end(); it++)
					{
						Storage::integer owner = (*it)->IntegerDF(tag_owner);
						if( owner == mpirank ) continue;
						Storage::integer_array proc = (*it)->IntegerArrayDV(tag_processors);	
						Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),mpirank);
						if( ip == proc.end() || (*ip) != mpirank ) 
						{
							proc.insert(ip,mpirank);
							send_elements[owner].push_back(*it);
						}
					}
				}
			}
		}
				
		if( action == AGhost ) //second round to inform owner about ghosted elements
		{
			REPORT_MPI(MPI_Waitall(send_reqs.size(),&send_reqs[0],MPI_STATUSES_IGNORE));
			//Inform owners about the fact that we have received their elements
			tag_list.clear();
			tag_list_recv.clear();
			
			send_bufs.clear();
			send_reqs.clear();
			recv_bufs.clear();
			recv_reqs.clear();
			send_bufs.resize(send_elements.size());
			
			num_wait = 0;
			for(std::map< int , std::vector<Element *> >::iterator it = send_elements.begin(); it != send_elements.end(); it++)
				if( !it->second.empty() )
				{
					PackElementsData(it->second,send_bufs[num_wait].second,it->first,tag_list);
					send_bufs[num_wait].first = it->first;
					num_wait++;
				}
			send_bufs.resize(num_wait);
				
			if( false )
			{ //in most cases the exchange of elements will effect only nearest neighbours, we want to detect this
				Storage::integer_array p = IntegerArrayDV(tag_processors);
				int halo_exchange = 0, test = 1;
				for(unsigned k = 0; k < send_bufs.size(); k++)
					if( !std::binary_search(p.begin(),p.end(),send_bufs[k].first) ) test = 0;
				//we should exchange, because other processor may want to send something to us
				REPORT_MPI(MPI_Allreduce(&test,&halo_exchange,1,MPI_INT,MPI_MIN,GetCommunicator())); 
				if( halo_exchange )
				{
					//prepare missing send buffers
					{
						std::vector<int> present(send_bufs.size()), missing(p.size());
						for(unsigned k = 0; k < send_bufs.size(); k++) present[k] = send_bufs[k].first;
						missing.resize(std::set_difference(p.begin(),p.end(),present.begin(),present.end(),missing.begin())-missing.begin());
						for(unsigned k = 0; k < missing.size(); k++) send_bufs.push_back(std::pair<int,std::vector<INMOST_DATA_BULK_TYPE> >(missing[k],std::vector<INMOST_DATA_BULK_TYPE>()));
					}
					//prepare recv buffers
					{
						recv_bufs.resize(p.size());
						for(unsigned k = 0; k < p.size(); k++) recv_bufs[k].first = p[k];
					}
					PrepareReceiveInner(UnknownSize, send_bufs,recv_bufs);
				}
				else PrepareReceiveInner(UnknownSource, send_bufs,recv_bufs);
			}
			else PrepareReceiveInner(UnknownSource, send_bufs,recv_bufs);
			ExchangeBuffersInner(send_bufs,recv_bufs,send_reqs,recv_reqs);
			
			send_elements.clear();
			
			while( !(done = FinishRequests(recv_reqs)).empty() )
			{
				for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
				{
					std::vector<Element *> recv_elements;
					UnpackElementsData(recv_elements,recv_bufs[*qt].second,recv_bufs[*qt].first,tag_list_recv);	
					for(std::vector<Element *>::iterator it = recv_elements.begin(); it != recv_elements.end(); it++)
					{
						Storage::integer_array proc = (*it)->IntegerArrayDV(tag_processors);
						Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),recv_bufs[*qt].first);
						if( ip == proc.end() || (*ip) != recv_bufs[*qt].first ) proc.insert(ip,recv_bufs[*qt].first);
						if( (*it)->IntegerDF(tag_owner) == mpirank )
							(*it)->BulkDF(tag_shared) = Element::Shared;
						else
							(*it)->BulkDF(tag_shared) = Element::Ghost;
					}
				}
			}
		}
		
		if( action == AMigrate ) //Compute new values
		{
			Tag tag_new_owner = GetTag("TEMPORARY_NEW_OWNER");
			Tag tag_new_processors = GetTag("TEMPORARY_NEW_PROCESSORS");	
			for(iteratorElement it = BeginElement(CELL | FACE | EDGE | NODE); it != EndElement(); it++)
			{
				Storage::integer new_owner = it->Integer(tag_new_owner);
				it->IntegerDF(tag_owner) = new_owner;
				Storage::integer_array old_procs = it->IntegerArrayDV(tag_processors);
				Storage::integer_array new_procs = it->IntegerArray(tag_new_processors);
				old_procs.clear();
				old_procs.insert(old_procs.end(),new_procs.begin(),new_procs.end());
				if( new_owner == mpirank )
				{
					if( old_procs.size() == 1 ) //there is only one processors and that's me
						it->BulkDF(tag_shared) = Element::Owned;
					else
						it->BulkDF(tag_shared) = Element::Shared;
				}
				else it->BulkDF(tag_shared) = Element::Ghost;
			}
		}
		
		RecomputeParallelStorage(CELL | FACE | EDGE | NODE);

		if( action == AGhost ) 
		{
			//wait for second round of communication
			REPORT_MPI(MPI_Waitall(num_wait,&send_reqs[0],MPI_STATUSES_IGNORE));
			//Probably now owner should send processors_tag data
			ExchangeData(tag_processors,CELL | FACE | EDGE | NODE);
		}
		ComputeSharedProcs();
		
		if( action == AMigrate ) REPORT_MPI(MPI_Waitall(send_reqs.size(),&send_reqs[0],MPI_STATUSES_IGNORE));
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}
	
	void Mesh::RecomputeParallelStorage(ElementType mask)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
#if defined(USE_MPI)
#if defined(USE_PARALLEL_STORAGE)
		for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
			for(int i = 0; i < 4; i++) if( mask & (1 << i) )
				it->second[i].clear();
		for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)			
			for(int i = 0; i < 4; i++) if( mask & (1 << i) )
				it->second[i].clear();
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif
	}

	
	void UnpackLayersMarker(Tag tag, Element * e, INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size )
		{
			int old_size;
			if( e->HaveData(tag) ) old_size = e->GetDataSize(tag);
			else old_size = 0;
			e->SetDataSize(tag,old_size+size);
			e->SetData(tag,old_size,size,data);
		}
	}
	
	void Mesh::ExchangeGhost(Storage::integer layers, ElementType bridge)
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
		if( m_state == Serial ) return;
#if defined(USE_MPI)
		int mpirank = GetProcessorRank();
		bool delete_ghost = false;
		//if( layers == Integer(tag_layers) && bridge == Integer(tag_bridge) ) return;
		if( layers < Integer(tag_layers) ) delete_ghost = true;
		else if( layers == Integer(tag_layers) && bridge < Integer(tag_bridge) ) delete_ghost = true;
		int test_bridge = 0;
		
		if( (bridge & MESH) || (bridge & ESET) || (bridge & CELL) ) throw Impossible;
		for(ElementType mask = NODE; mask <= FACE; mask = mask << 1)
			test_bridge += (bridge & mask)? 1 : 0;
		if( test_bridge == 0 || test_bridge > 1 ) throw Impossible;
		
		//RemoveGhost();
		Tag layers_marker = CreateTag("TEMPORARY_LAYERS_MARKER",DATA_INTEGER,CELL,CELL);
		Integer(tag_layers) = layers;
		Integer(tag_bridge) = bridge;
		Storage::integer_array procs = IntegerArrayDV(tag_processors);
		std::map< int, std::vector<Element *> > old_layers;
		std::map< int, std::vector<Element *> > current_layers;
		std::vector<Element *> all_visited;
		{
			std::map< int , std::vector<Element *> > shared_skin = ComputeSharedSkinSet(bridge);
			for(Storage::integer_array::iterator p = procs.begin(); p != procs.end(); p++)
			{
				MIDType busy = CreateMarker();
				all_visited.clear();
				for(std::vector<Element *>::iterator it = shared_skin[*p].begin(); it != shared_skin[*p].end(); it++)
				{
					adjacent<Element> adj = (*it)->getAdjElements(CELL);
					for(adjacent<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
						if( !jt->GetMarker(busy) )
						{
							//if( jt->IntegerDF(tag_owner) != *p )
							if( jt->IntegerDF(tag_owner) == mpirank )
							{
								current_layers[*p].push_back(&*jt);
								if( layers > 0 )
								{
									Storage::integer_array adj_procs = jt->IntegerArrayDV(tag_processors);
									if( !std::binary_search(adj_procs.begin(),adj_procs.end(),*p) )
										jt->IntegerArray(tag_sendto).push_back(*p);
									if( delete_ghost ) jt->IntegerArray(layers_marker).push_back(*p);
								}
							}
							jt->SetMarker(busy);
							all_visited.push_back(&*jt);
						}
				}
				ArrayRemMarker(all_visited,busy);
				ReleaseMarker(busy);
			}
		}
		for(Storage::integer k = layers-1; k >= 0; k--)
		{
			ExchangeMarked();
			old_layers.swap(current_layers);
			current_layers.clear();
			if( k > 0 ) 
			{
				for(Storage::integer_array::iterator p = procs.begin(); p != procs.end(); p++)
				{
					std::vector<Element *> & ref_cur = current_layers[*p];
					std::vector<Element *> & ref_old = old_layers[*p];
					MIDType busy = CreateMarker();
					all_visited.clear();
					ArraySetMarker(ref_old,busy);
					for(std::vector<Element *>::iterator it = ref_old.begin(); it != ref_old.end(); it++)
					{
						adjacent<Element> adj_bridge = (*it)->getAdjElements(bridge);
						for(adjacent<Element>::iterator jt = adj_bridge.begin(); jt != adj_bridge.end(); jt++)
							if( !jt->GetMarker(busy) )
							{
								adjacent<Element> adj = jt->getAdjElements(CELL);
								for(adjacent<Element>::iterator kt = adj.begin(); kt != adj.end(); kt++)
									if( !kt->GetMarker(busy) )
									{
										if( jt->IntegerDF(tag_owner) != *p )
										{
											ref_cur.push_back(&*kt);
											Storage::integer_array adj_procs = kt->IntegerArrayDV(tag_processors);
											if( !std::binary_search(adj_procs.begin(),adj_procs.end(),*p) )
												kt->IntegerArray(tag_sendto).push_back(*p);
											if( delete_ghost ) jt->IntegerArray(layers_marker).push_back(*p);
										}
										kt->SetMarker(busy);
										all_visited.push_back(&*kt);
									}
								jt->SetMarker(busy);
								all_visited.push_back(&*jt);
							}
					}
					ArrayRemMarker(ref_old,busy);
					ArrayRemMarker(all_visited,busy);
					ReleaseMarker(busy);
				}
			}
		}
		if( delete_ghost )
		{
			ReduceData(layers_marker,CELL,UnpackLayersMarker);
			ExchangeData(layers_marker,CELL);
			std::vector<Element *> del_ghost;
			for(iteratorElement it = BeginElement(CELL); it != EndElement(); it++)
			{
				if( it->BulkDF(tag_shared) == Element::Ghost )
				{
					if( !it->HaveData(layers_marker) )
						del_ghost.push_back(&*it);
					else
					{
						bool delete_element = true;
						Storage::integer_array arr = it->IntegerArray(layers_marker);
						for(Storage::integer_array::iterator kt = arr.begin(); kt != arr.end(); kt++)
							if( *kt == mpirank )
							{
								delete_element = false;
								break;
							}
						if( delete_element ) 
							del_ghost.push_back(&*it);
					}
				}
			}
			RemoveGhostElements(del_ghost);
		}
		DeleteTag(layers_marker);
		//throw NotImplemented;
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif

	}
	
	
	void RedistUnpack(Tag tag,Element * e,INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size ) 
		{
			Storage::integer_array p1 = e->IntegerArray(tag);
			Storage::integer * p2 = static_cast<Storage::integer *>(static_cast<void *>(data));
			std::vector<Storage::integer> result(p1.size()+size);
			std::vector<Storage::integer>::iterator end;
			end = std::set_union(p1.begin(),p1.end(),p2,p2+size,result.begin());
			result.resize(end-result.begin());
			p1.clear();
			p1.insert(p1.end(),result.begin(),result.end());
		}
	}

	
	void Mesh::Redistribute()
	{
#if defined(USE_PARALLEL_WRITE_TIME)
		long double all_time = Timer();
		WriteTab(out_time) << __FUNCTION__ << std::endl;
		Enter();
#endif
		if( m_state == Serial ) return;
#if defined(USE_MPI)
		int mpirank = GetProcessorRank();
		Tag tag_new_owner = CreateTag("TEMPORARY_NEW_OWNER",DATA_INTEGER,FACE | EDGE | NODE, NONE,1);
		if( !tag_new_owner.isDefined(CELL) ) throw DistributionTagWasNotFilled;
		Tag tag_new_processors = CreateTag("TEMPORARY_NEW_PROCESSORS",DATA_INTEGER,CELL | FACE | EDGE | NODE, NONE);
		std::map<int, std::vector<Element *> > redistribute_skin, skin_faces;
		std::vector<Element *> all_visited;
		ElementType bridge = Integer(tag_bridge);
		Storage::integer layers = Integer(tag_layers);
		
		ExchangeData(tag_new_owner,CELL);
		
		
		for(iteratorElement it = BeginElement(CELL); it != EndElement(); it++)
			it->IntegerArray(tag_new_processors).push_back(it->Integer(tag_new_owner));
		
		
		//determine tag_new_processors for FACEs to calculate shared skin
		for(ElementType mask = FACE; mask >= NODE; mask = mask >> 1)
			for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			{
				Storage::integer_array procs = it->IntegerArray(tag_new_processors);
				adjacent<Element> over = it->getAdjElements(mask << 1);
				determine_my_procs_high(over,procs,tag_new_processors);
				if( !procs.empty() )
					it->Integer(tag_new_owner) = procs[0];
				else
				{
					it->Integer(tag_new_owner) = mpirank;
					procs.push_back(mpirank);
				}
			}
			
		ExchangeData(tag_new_owner,FACE | EDGE | NODE);

		if( bridge != NONE && layers != 0 )
		{			
			ReduceData(tag_new_processors,FACE ,RedistUnpack);
			ExchangeData(tag_new_processors,FACE);
			//compute skin around migrated cells to compute ghost layer
			for(iteratorElement it = BeginElement(FACE); it != EndElement(); it++)
			{
				Storage::integer_array procs = it->IntegerArray(tag_new_processors);
				if( procs.size() == 2 ) //there may be only 2 procs for every face max, because we have unique redistribution for every cell at the beginning
				{
					skin_faces[procs[0]].push_back(&*it);
					skin_faces[procs[1]].push_back(&*it);
				}
			}
			if( bridge == FACE ) redistribute_skin.swap(skin_faces);
			else
			{
				MIDType busy = CreateMarker();
				// Here should first find all adjacent elements of skin_faces of given type,
				// then distribute them between sets according to new processors
				for(std::map<int, std::vector<Element *> >::iterator it = skin_faces.begin(); it != skin_faces.end(); it++)
				{
					std::vector<Element *> & ref = redistribute_skin[it->first];
					for(std::vector<Element *>::iterator jt = it->second.begin(); jt != it->second.end(); jt++)
					{
						adjacent<Element> adj = (*jt)->getAdjElements(bridge);
						for(adjacent<Element>::iterator kt = adj.begin(); kt != adj.end(); kt++)
							if( !kt->GetMarker(busy) )
							{
								ref.push_back(&*kt);
								kt->SetMarker(busy);
							}
					}
					ArrayRemMarker(ref,busy);
				}
				ReleaseMarker(busy);
			}
			skin_faces.clear();
			//now mark all cell layers 
			{
				std::map<int, std::vector<Element *> > current_layers,old_layers;
				MIDType busy = CreateMarker();
				for(std::map<int, std::vector<Element *> >::iterator it = redistribute_skin.begin(); it != redistribute_skin.end(); it++)
				{
					all_visited.clear();
					for(std::vector<Element *>::iterator jt = it->second.begin(); jt != it->second.end(); jt++)
					{
						std::vector<Element *> & ref = current_layers[it->first];
						adjacent<Element> adj = (*jt)->getAdjElements(CELL);
						for(adjacent<Element>::iterator kt = adj.begin(); kt != adj.end(); kt++)
							if( !kt->GetMarker(busy) )
							{
								Storage::integer_array procs = kt->IntegerArray(tag_new_processors);
								Storage::integer_array::iterator find = std::lower_bound(procs.begin(),procs.end(),it->first);
								if( find == procs.end() || *find != it->first ) 
								{
									procs.insert(find,it->first);
									ref.push_back(&*kt);
								}
								all_visited.push_back(&*kt);
								kt->SetMarker(busy);	
							}
					}
					ArrayRemMarker(all_visited,busy);
				}
				ReleaseMarker(busy);
				for(Storage::integer k = layers-1; k > 0; k--)
				{
					old_layers.swap(current_layers);
					current_layers.clear();
					for(std::map< int, std::vector<Element *> >::iterator qt = old_layers.begin(); qt != old_layers.end(); qt++)
					{
						MIDType busy = CreateMarker();
						all_visited.clear();
						ArraySetMarker(qt->second,busy);
						for(std::vector<Element *>::iterator it = qt->second.begin(); it != qt->second.end(); it++)
						{
							std::vector<Element *> & ref = current_layers[qt->first];
							adjacent<Element> adj_bridge = (*it)->getAdjElements(bridge);
							for(adjacent<Element>::iterator jt = adj_bridge.begin(); jt != adj_bridge.end(); jt++)
								if( !jt->GetMarker(busy) )
								{
									adjacent<Element> adj = jt->getAdjElements(CELL);
									for(adjacent<Element>::iterator kt = adj.begin(); kt != adj.end(); kt++)
										if( !kt->GetMarker(busy) )
										{
											Storage::integer_array procs = kt->IntegerArray(tag_new_processors);
											Storage::integer_array::iterator find = std::lower_bound(procs.begin(),procs.end(),qt->first);
											if( find == procs.end() || *find != qt->first )
											{
												procs.insert(find,qt->first);
												ref.push_back(&*kt);
											}
											kt->SetMarker(busy);
											all_visited.push_back(&*kt);
										}
									jt->SetMarker(busy);
									all_visited.push_back(&*jt);
								}
						}
						ArrayRemMarker(qt->second,busy);
						ArrayRemMarker(all_visited,busy);
						ReleaseMarker(busy);
					}
				}
			}
			
			ReduceData(tag_new_processors,CELL,RedistUnpack);
			ExchangeData(tag_new_processors,CELL);
			
			for(ElementType mask = FACE; mask >= NODE; mask = mask >> 1)
			{
				for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
				{
					Storage::integer_array procs = it->IntegerArray(tag_new_processors);
					adjacent<Element> over = it->getAdjElements(mask << 1);
					determine_my_procs_high(over,procs,tag_new_processors);
					if( procs.empty() ) procs.push_back(mpirank);
				}
			}
			
			ReduceData(tag_new_processors,FACE| EDGE| NODE,RedistUnpack);
			ExchangeData(tag_new_processors,FACE|EDGE|NODE);
		}
		else 
		{
			ReduceData(tag_new_processors,FACE | EDGE | NODE,RedistUnpack);
			ExchangeData(tag_new_processors,FACE | EDGE | NODE);
		}
		


		for(iteratorElement it = BeginElement(CELL | FACE | EDGE | NODE); it != EndElement(); it++)
		{
			Storage::integer_array new_procs = it->IntegerArray(tag_new_processors);
			if( it->IntegerDF(tag_owner) == mpirank ) // deal with my entities, others will deal with theirs
			{
				//compute processors thay should have the entity but they have not
				Storage::integer_array old_procs = it->IntegerArrayDV(tag_processors);
				std::vector<Storage::integer> result(new_procs.size());
				std::vector<Storage::integer>::iterator end;
				end = std::set_difference(new_procs.begin(),new_procs.end(),old_procs.begin(),old_procs.end(),result.begin());
				result.resize(end-result.begin());
				//mark to send entity to processors that don't have it
				for(std::vector<Storage::integer>::iterator qt = result.begin(); qt != result.end(); qt++)
					it->IntegerArray(tag_sendto).push_back(*qt);
			}
		}


		ExchangeMarked(AMigrate);
		DeleteTag(tag_new_owner);
		DeleteTag(tag_new_processors);
		//throw NotImplemented;
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		Exit();
		WriteTab(out_time) << __FUNCTION__ << " time: " << Timer() - all_time << std::endl;
#endif

	}
	void Mesh::BeginSequentialCode() 
	{
#if defined(USE_MPI)
		for(int i = 0; i < GetProcessorRank(); i++) 
			MPI_Barrier(GetCommunicator());
#endif
	}
	void Mesh::EndSequentialCode() 
	{
#if defined(USE_MPI)
		for(int i = GetProcessorRank(); i < GetProcessorsNumber(); i++) 
			MPI_Barrier(GetCommunicator());
#endif
	}
}

#endif
