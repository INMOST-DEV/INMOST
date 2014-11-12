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
#define REPORT_MPI(x) {WriteTab(out_time) << "<MPI><![CDATA[" << #x << "]]></MPI>" << std::endl; x;}
#define REPORT_STR(x) {WriteTab(out_time) << "<TEXT><![CDATA[" << x << "]]></TEXT>" << std::endl;}
#define REPORT_VAL(str,x) {WriteTab(out_time) << "<VALUE name=\"" << str << "\"> <CONTENT><![CDATA[" << x << "]]></CONTENT> <CODE><![CDATA[" << #x << "]]></CODE></VALUE>" << std::endl;}
#define ENTER_FUNC() long double all_time = Timer(); WriteTab(out_time) << "<FUNCTION name=\"" << __FUNCTION__ << "\" id=\"func" << func_id++ << "\">" << std::endl; Enter();
#define EXIT_FUNC() WriteTab(out_time) << "<TIME>" << Timer() - all_time << "</TIME>" << std::endl; Exit(); WriteTab(out_time) << "</FUNCTION>" << std::endl;
#else
#define REPORT_MPI(x) x
#define REPORT_STR(x)
#define REPORT_VAL(str,x)
#define ENTER_FUNC()
#define EXIT_FUNC() 
#endif

namespace INMOST
{

	void UnpackSyncMarkerOR(Tag tag, Element * element, INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		element->Bulk(tag) |= *data;
	}

	void UnpackSyncMarkerXOR(Tag tag, Element * element, INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		element->Bulk(tag) ^= *data;
	}

	
	void UnpackSyncMarkerAND(Tag tag, Element * element, INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		element->Bulk(tag) &= *data;
	}

	void Mesh::SynchronizeMarker(MarkerType marker, ElementType mask, SyncBitOp op)
	{
#if defined(USE_MPI)
		if( m_state == Mesh::Parallel )
		{
			Tag t = CreateTag("TEMP_SYNC_MARKER",DATA_BULK,mask,mask,1);

			//workaround for old gcc compiler
			const Element::Status SGhost = Element::Ghost;
			const Element::Status SAny = Element::Any;
			Element::Status Expr = (Element::Shared | ((op != SYNC_BIT_NEW) ? SGhost : SAny));

			for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); ++it)
				if( it->GetMarker(marker) && (it->GetStatus() & Expr) )
					it->Bulk(t) = 1;

			
			
			switch(op)
			{
			case SYNC_BIT_NEW:
				ExchangeData(t,mask);
				for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); ++it)
				{
					if( it->GetStatus() == Element::Ghost )
					{
						if( it->HaveData(t) ) it->SetMarker(marker); else it->RemMarker(marker);
					}
				}
				break;
			case SYNC_BIT_OR:
				ReduceData(t,mask,UnpackSyncMarkerOR);
				for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); ++it)
				{
					if( it->GetStatus() & (Element::Ghost | Element::Shared) )
					{
						if( !it->GetMarker(marker) && it->HaveData(t) ) it->SetMarker(marker);
					}
				}
				break;
			case SYNC_BIT_AND:
				ReduceData(t,mask,UnpackSyncMarkerAND);
				for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); ++it)
				{
					if( it->GetStatus() & (Element::Ghost | Element::Shared))
					{
						if( it->HaveData(t) && it->Bulk(t) ) it->SetMarker(marker); else it->RemMarker(marker);
					}
				}
				break;
			case SYNC_BIT_XOR:
				ReduceData(t,mask,UnpackSyncMarkerXOR);
				for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); ++it)
				{
					if( it->GetStatus() & (Element::Ghost | Element::Shared))
					{
						if( it->HaveData(t) && it->Bulk(t) ) it->SetMarker(marker); else it->RemMarker(marker);
					}
				}
				break;
			}
			

			DeleteTag(t,mask);
		}
#endif
	}

	ElementType Mesh::SynchronizeElementType(ElementType etype)
	{
		ElementType etypeout = etype;
#if defined(USE_MPI)
		MPI_Allreduce(&etype,&etypeout,1,INMOST_MPI_DATA_BULK_TYPE,MPI_BOR,GetCommunicator());
#endif
		return etypeout;
	}
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
			if( it->GetStatus() == Element::Owned )
				it->Integer(num_tag) = shift++;
		for(ElementSet::iterator it = set->begin(); it != set->end(); it++)
			if( it->GetStatus() == Element::Shared )
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
			if( (*it)->GetStatus() != Element::Owned )
				(*it)->Integer(num_tag) = shift++;
		for(std::vector<Element *>::iterator it = set.begin(); it != set.end(); it++)
			if( (*it)->GetStatus() != Element::Shared )
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
			if( it->GetStatus() != Element::Owned )
				it->Integer(num_tag) = shift++;
		for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			if( it->GetStatus() != Element::Shared )
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
	
	Storage::real Mesh::AggregateMax(Storage::real input)
	{
		Storage::real output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_REAL_TYPE,MPI_MAX,comm);
#endif
		return output;
	}
	
	Storage::integer Mesh::AggregateMax(Storage::integer input)
	{
		Storage::integer output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_MAX,comm);
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
				else if( tag.GetSize() == ENUMUNDEF )
					element->SetDataSize(tag,size);
				else
				{
					std::cerr << "received zero size for dense nonzero tag" << std::endl;
					throw Impossible;
				}
			}
			return;
		}
		if( !element->HaveData(tag) )
			element->SetDataSize(tag,size);
		else if( size != element->GetDataSize(tag) )
		{
			if( tag.GetSize() == ENUMUNDEF )
				element->SetDataSize(tag,size);
			else
			{
				throw Impossible;
			}
		}
		element->SetData(tag,0,size,data);
	}

	void ArraySetMarker(std::vector<Element *> & arr, MarkerType marker)
	{
		for(std::vector<Element *>::iterator it = arr.begin(); it != arr.end(); it++)
			(*it)->SetMarker(marker);
	}
	
	void ArrayRemMarker(std::vector<Element *> & arr, MarkerType marker)
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
		int test = 0;
		MPI_Finalized(&test);
		if( !test )
			MPI_Finalize();
#endif
	}
	
	
	
	void Mesh::SetCommunicator(INMOST_MPI_Comm _comm)
	{
		ENTER_FUNC();
		tag_shared = CreateTag("PROTECTED_P_SHARED",DATA_BULK,CELL | FACE | EDGE | NODE,NONE,1);
		tag_owner = CreateTag("P_OWNER_PROCESSOR",DATA_INTEGER, CELL | FACE | EDGE | NODE,NONE,1);
		tag_processors = CreateTag("P_PROCESSORS_LIST",DATA_INTEGER, MESH | NODE | EDGE | FACE | CELL,NONE);
		tag_layers = CreateTag("P_LAYERS",DATA_INTEGER,MESH,NONE,1);
		tag_bridge = CreateTag("P_BRIDGE",DATA_INTEGER,MESH,NONE,1);
		tag_sendto = CreateTag("PROTECTED_P_SENDTO",DATA_INTEGER, CELL | FACE | EDGE | NODE, CELL | FACE | EDGE | NODE);

#if defined(USE_MPI)
		
		
		parallel_strategy = 1;
		parallel_file_strategy = 1;

		

		Mesh::Initialize(NULL,NULL);
		//~ MPI_Comm_dup(_comm,&comm);
		comm = _comm;
		{
			INMOST_DATA_BIG_ENUM_TYPE t = pmid;
			REPORT_MPI(MPI_Allreduce(&t,&parallel_mesh_unique_id,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,MPI_MAX,MPI_COMM_WORLD));
			pmid = parallel_mesh_unique_id+1;
		}
		m_state = Mesh::Parallel;


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
		EXIT_FUNC();
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
	
	void determine_my_procs_low(adjacent<Element> & subelements,Tag & procs_tag, dynarray<Storage::integer,64> & result, dynarray<Storage::integer,64> & intersection)
	{
		if( subelements.empty() ) return;
		adjacent<Element>::iterator i = subelements.begin();
		Storage::integer_array p = i->IntegerArrayDV(procs_tag);
		result.clear();
		result.insert(result.end(),p.begin(),p.end());
		i++;
		while(i != subelements.end())
		{
			Storage::integer_array q = i->IntegerArrayDV(procs_tag);
			intersection.resize(std::max(result.size(),q.size()));
			dynarray<Storage::integer,64>::iterator qt = std::set_intersection(result.begin(),result.end(),q.begin(),q.end(),intersection.begin());
			intersection.resize(qt-intersection.begin());
			result.swap(intersection);
			i++;
		}
	}
	
	void determine_my_procs_high(adjacent<Element> & overelements,Tag & procs_tag, dynarray<Storage::integer,64> & result, dynarray<Storage::integer,64> & intersection)
	{
		if( overelements.empty() ) return;
		adjacent<Element>::iterator i = overelements.begin();
		result.clear();
		while(i != overelements.end())
		{
			Storage::integer_array q = i->IntegerArrayDV(procs_tag);
			intersection.resize(result.size()+q.size());
			dynarray<Storage::integer,64>::iterator qt = std::set_union(result.begin(),result.end(),q.begin(),q.end(),intersection.begin());
			intersection.resize(qt-intersection.begin());
			result.swap(intersection);
			i++;
		}
	}
	
	/*
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
	*/
	/*
	void Mesh::MarkShared(ElementType mask)
	{
		for(ElementType etype = NODE; etype <= CELL; etype = etype << 1 ) if( etype & mask )
		for(Mesh::iteratorElement it = BeginElement(etype); it != EndElement(); it++)
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
		RecomputeParallelStorage(mask);
		AssignGlobalID(mask);
		
		//have to do it for all types
#if defined(USE_PARALLEL_STORAGE)
		for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
			qsort(&it->second[0][0],it->second[0].size(),sizeof(Element *),CompareElementsCGID);
		for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)			
			qsort(&it->second[0][0],it->second[0].size(),sizeof(Element *),CompareElementsCGID);
#endif
	}
	*/
	
	
	__INLINE bool point_in_bbox(Storage::real * p, Storage::real bbox[6], unsigned dim)
	{
		bool ret = true;
		for(unsigned k = 0; k < dim; k++) 
			ret &= (p[k] >= bbox[k] && p[k] <= bbox[dim+k]);
		return ret;
	}
	
	__INLINE bool compare_coord(Storage::real * a, Storage::real * b, INMOST_DATA_ENUM_TYPE dim, Storage::real eps)
	{
		for(INMOST_DATA_ENUM_TYPE i = 0; i <  dim; i++)
			if( fabs(a[i]-b[i]) > eps ) return a[i] <= b[i];
		return true;
	}
	
	void Mesh::ResolveShared()
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		if( m_state == Mesh::Serial ) SetCommunicator(INMOST_MPI_COMM_WORLD);
		unsigned int dim = GetDimensions();
		int sendsize;
		int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber();
#if defined(USE_PARALLEL_STORAGE)
		shared_elements.clear();
		ghost_elements.clear();
#endif
		//determine which bboxes i intersect
		dynarray<int,64> procs;
		Storage::real bbox[6];
		dynarray<Storage::real,384> bboxs(mpisize*6);
		
		for(unsigned k = 0; k < dim; k++)
		{
			bbox[k] = 1e20;
			bbox[k+dim] = -1e20;
		}
		
		for(iteratorNode it = BeginNode(); it != EndNode(); it++)
		{
			Storage::real_array arr = it->Coords();
			for(unsigned k = 0; k < dim; k++)
			{
				if( arr[k] < bbox[k] ) bbox[k] = arr[k];
				if( arr[k] > bbox[k+dim] ) bbox[k+dim] = arr[k];
			}
		}
		for(unsigned k = 0; k < dim; k++)
		{
			REPORT_VAL("min",bbox[k]);
			REPORT_VAL("max",bbox[dim+k]);
		}
		REPORT_MPI(MPI_Allgather(&bbox[0],dim*2,INMOST_MPI_DATA_REAL_TYPE,&bboxs[0],dim*2,INMOST_MPI_DATA_REAL_TYPE,comm));
		for(int k = 0; k < mpisize; k++)
			if( k != mpirank )
			{
				bool flag = true;
				for(unsigned q = 0; q < dim; q++)
					flag &= !((bbox[q] > bboxs[k*dim*2+q+dim]) || (bbox[dim+q] < bboxs[k*dim*2+q]));
				if( flag ) procs.push_back(k);
			}
		REPORT_VAL("neighbour processors",procs.size());
		
		//~ if( procs.empty() )
		//~ {
			//~ REPORT_STR("no processors around - all elements are owned");
			//~ for(Mesh::iteratorElement it = BeginElement(CELL | EDGE | FACE | NODE); it != EndElement(); it++)
			//~ {
				//~ it->IntegerArrayDV(tag_processors).resize(1);
				//~ it->IntegerDF(tag_owner) = it->IntegerDV(tag_processors) = mpirank;
				//~ it->BulkDF(tag_shared) = Element::Owned;
			//~ }
			//~ ComputeSharedProcs();
			//~ RecomputeParallelStorage(CELL | EDGE | FACE | NODE);
			//~ AssignGlobalID(CELL | EDGE | FACE | NODE);
		//~ }
		//~ else
		{
			bool same_boxes = true, same_box;
			for(int k = 0; k < mpisize && same_boxes; k++)
			{
				same_box = true;
				for(unsigned j = 0; j < dim*2; j++)
					same_box &= fabs(bbox[j] - bboxs[k*dim*2+j]) < epsilon;
				same_boxes &= same_box;
			}
			
			if( same_boxes )
			{
				REPORT_STR("All bounding boxes are the same - assuming that mesh is replicated over all nodes");
				for(Mesh::iteratorElement it = BeginElement(CELL | EDGE | FACE | NODE); it != EndElement(); it++)
				{
					Storage::integer_array arr = it->IntegerArrayDV(tag_processors);
					arr.resize(mpisize);
					for(int k = 0; k < mpisize; k++) arr[k] = k;
					it->IntegerDF(tag_owner) = 0;
					if( mpirank == 0 ) it->BulkDF(tag_shared) = Element::Shared;
					else it->BulkDF(tag_shared) = Element::Ghost;
				}
				ComputeSharedProcs();
				RecomputeParallelStorage(CELL | EDGE | FACE | NODE);
				AssignGlobalID(CELL | EDGE | FACE | NODE);
#if defined(USE_PARALLEL_STORAGE)
				for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
					for(int i = 0; i < 4; i++)
					{
						if( !it->second[i].empty() )
							qsort(&it->second[i][0],it->second[i].size(),sizeof(Element *),CompareElementsCGID);
					}
				for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)			
					for(int i = 0; i < 4; i++)
					{
						if( !it->second[i].empty() )
							qsort(&it->second[i][0],it->second[i].size(),sizeof(Element *),CompareElementsCGID);
					}
#endif
			}
			else
			{
				double time = Timer();
				Storage::real epsilon = GetEpsilon();
			
				std::map<GeometricData,ElementType> table;
				for(ElementType etype = EDGE; etype <= CELL; etype = etype << 1)
					if( !HaveGeometricData(CENTROID,etype) )
						table[CENTROID] |= etype;
				PrepareGeometricData(table);
				
				time = Timer() - time;
			
				REPORT_STR("Prepare geometric data");
				REPORT_VAL("time",time);
			
				
			
				for(iteratorNode it = BeginNode(); it != EndNode(); it++)
				{
					Storage::integer_array arr = it->IntegerArrayDV(tag_processors);
					arr.resize(1);
					arr[0] = mpirank;
				}
			

				std::vector< INMOST_DATA_BULK_TYPE > exch_data;
				std::vector< INMOST_DATA_REAL_TYPE > unpack_real;
				std::vector< INMOST_DATA_REAL_TYPE > pack_real;
				std::vector< Node * > sorted_nodes;
				
				sorted_nodes.reserve(NumberOfNodes());
				
				
				time = Timer();
				
				
				for(iteratorNode n = BeginNode(); n != EndNode(); n++)
				{
					Storage::real_array c = n->Coords();
					for(unsigned k = 0; k < procs.size(); k++)
						if( point_in_bbox(&c[0],&bboxs[procs[k]*dim*2],dim) )
						{
							sorted_nodes.push_back(&*n);
							break;
						}
				}
				
				time = Timer() - time;
				REPORT_STR("Prepare array of nodes");
				REPORT_VAL("time",time);
				REPORT_VAL("share nodes", sorted_nodes.size());
				REPORT_VAL("total nodes", NumberOfNodes());
				
				
				
				time = Timer();
				if( !sorted_nodes.empty() )
					qsort(&sorted_nodes[0],sorted_nodes.size(),sizeof(Element *),CompareElementsCCentroid);
				time = Timer() - time;
				REPORT_STR("Sort nodes");
				REPORT_VAL("time",time);
				
				
				pack_real.reserve(sorted_nodes.size()*dim);
				
				
				time = Timer();
				for(std::vector<Node *>::iterator it = sorted_nodes.begin(); it != sorted_nodes.end(); it++)
				{
					Storage::real_array arr = (*it)->Coords();
					pack_real.insert(pack_real.end(),arr.begin(),arr.end());
				}
				time = Timer() - time;
				REPORT_STR("Gather coordinates");
				REPORT_VAL("time",time);
				
				
				time = Timer();
			
				MPI_Pack_size(pack_real.size(),INMOST_MPI_DATA_REAL_TYPE,comm,&sendsize);
				exch_data.resize(sendsize);
				int position = 0;
				if( sendsize > 0 ) MPI_Pack(&pack_real[0],pack_real.size(),INMOST_MPI_DATA_REAL_TYPE,&exch_data[0],exch_data.size(),&position,comm);
				
				
				time = Timer() - time;
				REPORT_STR("Pack coordinates");
				REPORT_VAL("time",time);

				
				{
					std::vector< MPI_Request > send_reqs, recv_reqs;
					exch_buffer_type send_buffs(procs.size()), recv_buffs(procs.size());
					std::vector<int> done;
//~ #if defined(USE_MPI2)
					//~ unsigned * sendsizeall = shared_space;
					//~ unsigned usend[2] = {sendsize,pack_real.size()};
					//~ REPORT_MPI(MPI_Win_fence(MPI_MODE_NOPRECEDE,window)); //start exchange session
					//~ for(unsigned k = 0; k < procs.size(); k++)
						//~ REPORT_MPI(MPI_Put(usend,2,MPI_UNSIGNED,procs[k],mpirank*2,2,MPI_UNSIGNED,window));
					//~ REPORT_MPI(MPI_Win_fence(MPI_MODE_NOSUCCEED,window)); //end exchange session
//~ #else
					std::vector<unsigned> sendsizeall(mpisize*2);
					int pack_size2 = 0;
					unsigned long usend[2] = {(unsigned)sendsize,pack_real.size()};
					MPI_Pack_size(2,MPI_UNSIGNED,comm,&pack_size2);
					
					for(unsigned k = 0; k < procs.size(); k++)
					{
						send_buffs[k].first = procs[k];
						send_buffs[k].second.resize(pack_size2);
						position = 0;
						MPI_Pack(usend,2,MPI_UNSIGNED,&send_buffs[k].second[0],send_buffs[k].second.size(),&position,comm);
						recv_buffs[k].first = procs[k];
						recv_buffs[k].second.resize(pack_size2);
					}
					ExchangeBuffersInner(send_buffs,recv_buffs,send_reqs,recv_reqs);
					while( !(done = FinishRequests(recv_reqs)).empty() )
					{
						for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
						{
							position = 0;
							MPI_Unpack(&recv_buffs[*qt].second[0],recv_buffs[*qt].second.size(),&position,&sendsizeall[procs[*qt]*2],2,MPI_UNSIGNED,comm);
						}
					}
					if( !send_reqs.empty() )
					{
						REPORT_MPI(MPI_Waitall(send_reqs.size(),&send_reqs[0],MPI_STATUSES_IGNORE));
					}
					//~ REPORT_MPI(MPI_Allgather(usend,2,MPI_UNSIGNED,&sendsizeall[0],2,MPI_UNSIGNED,comm));
//~ #endif
					double time2 = Timer();
					{
						
						
						for(unsigned k = 0; k < procs.size(); k++)
						{
							send_buffs[k].first = procs[k];
							send_buffs[k].second = exch_data;
							recv_buffs[k].first = procs[k];
							recv_buffs[k].second.resize(sendsizeall[procs[k]*2]);
						}
						
						//PrepareReceiveInner(send_buffs,recv_buffs);
						ExchangeBuffersInner(send_buffs,recv_buffs,send_reqs,recv_reqs);
						
						
						
						while( !(done = FinishRequests(recv_reqs)).empty() )
						{
							for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
							{
								time = Timer();
								REPORT_STR("receive node coordinates");
								REPORT_VAL("processor",recv_buffs[*qt].first);
								int count = 0;
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
										count++;
										sorted_nodes[(it1-pack_real.begin())/dim]->IntegerArrayDV(tag_processors).push_back(recv_buffs[*qt].first);
										it1 += dim;
										it2 += dim;
									}
								}
								REPORT_VAL("intersected coords",count);
								time = Timer() - time;
								REPORT_STR("Intersect coordinates");
								REPORT_VAL("time",time);
							}
						}
						if( !send_reqs.empty() )
						{
							REPORT_MPI(MPI_Waitall(send_reqs.size(),&send_reqs[0],MPI_STATUSES_IGNORE));
						}
					}			
					
					time2 = Timer() - time2;
					REPORT_STR("Intersect all coordinates");
					REPORT_VAL("time",time2);
					
					time = Timer();
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
					{
						if( !it->second[0].empty() )
							qsort(&it->second[0][0],it->second[0].size(),sizeof(Element *),CompareElementsCGID);
					}
					for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)
					{
						if( !it->second[0].empty() )
							qsort(&it->second[0][0],it->second[0].size(),sizeof(Element *),CompareElementsCGID);
					}
	#endif
					
					
					time = Timer() - time;
					REPORT_STR("Set parallel info for nodes");
					REPORT_VAL("time",time);
					
					
					dynarray<Storage::integer,64> result, intersection;
					
					for(ElementType current_mask = EDGE; current_mask <= CELL; current_mask = current_mask << 1 )
					{
						REPORT_STR("Set parallel info for");
						REPORT_VAL("type",ElementTypeName(current_mask));
						
						
						int owner;
						Element::Status estat;
						
						time = Timer();
						//Determine what processors potentially share the element
						for(Mesh::iteratorElement it = BeginElement(current_mask); it != EndElement(); it++)
						{
							adjacent<Element> sub = it->getAdjElements(current_mask >> 1);
							
							determine_my_procs_low(sub,tag_processors, result, intersection);
							Storage::integer_array p = it->IntegerArrayDV(tag_processors);
							if( result.empty() ) 
							{
								p.clear();
								p.push_back(mpirank);
							}
							else p.replace(p.begin(),p.end(),result.begin(),result.end());
						}
						time = Timer() - time;
						REPORT_STR("Predict processors for elements");
						REPORT_VAL("time",time);
						
						
						time = Timer();
						//Initialize mapping that helps get local id by global id
						std::vector<int> mapping;
						for(Mesh::iteratorElement it = BeginElement(current_mask >> 1); it != EndElement(); it++)
						{
							mapping.push_back(it->GlobalID());
							mapping.push_back(it->LocalID());
						}
						if( !mapping.empty() ) 
							qsort(&mapping[0],mapping.size()/2,sizeof(int)*2,comp_mapping);
						time = Timer() - time;
						REPORT_STR("Compute global to local indexes mapping");
						REPORT_VAL("time",time);
						//Initialize arrays
						
						time = Timer();
						Storage::integer_array procs = IntegerArrayDV(tag_processors);
						Storage::integer_array::iterator p = procs.begin();
						std::vector<int> message_send;
						std::vector< std::vector<int> > message_recv(procs.size());
						std::vector< std::vector<Element *> > elements(procs.size());
						
						std::vector< MPI_Request > send_reqs,recv_reqs;
						exch_buffer_type send_buffs(procs.size()), recv_buffs(procs.size());
						
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
						
						time = Timer() - time;
						REPORT_STR("Pack messages for other processors");
						REPORT_VAL("time",time);
						
						
						time = Timer();
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
						
						time = Timer() - time;
						REPORT_STR("Exchange messages");
						REPORT_VAL("time",time);
						
						time = Timer();
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
						
						time = Timer() - time;
						REPORT_STR("Determine true shared based on remote info");
						REPORT_VAL("time",time);
						
						time = Timer();
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
						if( !send_reqs.empty() )
						{
							REPORT_MPI(MPI_Waitall(send_reqs.size(),&send_reqs[0],MPI_STATUSES_IGNORE));
						}
						AssignGlobalID(current_mask);					
#if defined(USE_PARALLEL_STORAGE)
						for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
						{
							if( !it->second[ElementNum(current_mask)].empty() )
								qsort(&it->second[ElementNum(current_mask)][0],it->second[ElementNum(current_mask)].size(),sizeof(Element *),CompareElementsCGID);
						}
						for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)
						{
							if( !it->second[ElementNum(current_mask)].empty() )
								qsort(&it->second[ElementNum(current_mask)][0],it->second[ElementNum(current_mask)].size(),sizeof(Element *),CompareElementsCGID);
						}
#endif
	
						time = Timer() - time;
						REPORT_STR("Set parallel info");
						REPORT_VAL("time",time);
					}
				}
				RemoveGeometricData(table);	
			}
		}
		//RemoveGhost();
		
#else
		AssignGlobalID(CELL | FACE | EDGE | NODE);
#endif
		EXIT_FUNC();
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
		ENTER_FUNC()
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
			if( !it->second.empty() ) qsort(&it->second[0],it->second.size(),sizeof(Element *),have_global_id & CELL ? CompareElementsCGID : CompareElementsCCentroid );
			std::vector<Element *> result(ref.size());
			std::vector<Element *>::iterator end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),have_global_id & CELL ? LessElementsGID : LessElementsCentroid);
			result.resize(end-result.begin());
			ref.swap(result);
		}
		del_ghost.clear();
		for(std::map< int, std::vector<Element *> >::iterator it = del_shared.begin(); it != del_shared.end(); it++)
		{
			std::vector<Element * > & ref = shared_elements[it->first][ElementNum(CELL)];
			if( !it->second.empty() ) qsort(&it->second[0],it->second.size(),sizeof(Element *),have_global_id & CELL ? CompareElementsCGID : CompareElementsCCentroid);
			std::vector<Element *> result(ref.size());
			std::vector<Element *>::iterator end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),have_global_id & CELL ? LessElementsGID : LessElementsCentroid);
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
				if( !it->second.empty() ) qsort(&it->second[0],it->second.size(),sizeof(Element *),have_global_id & mask ? CompareElementsCGID : CompareElementsCCentroid);
				std::vector<Element *> result(ref.size());
				std::vector<Element *>::iterator end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),have_global_id & mask ? LessElementsGID : LessElementsCentroid);
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_ghost.clear();
			for(std::map< int, std::vector<Element *> >::iterator it = del_shared.begin(); it != del_shared.end(); it++)
			{
				std::vector<Element * > & ref = shared_elements[it->first][ElementNum(mask)];
				if( !it->second.empty() ) qsort(&it->second[0],it->second.size(),sizeof(Element *),have_global_id & mask ? CompareElementsCGID : CompareElementsCCentroid);
				std::vector<Element *> result(ref.size());
				std::vector<Element *>::iterator end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),have_global_id & mask ? LessElementsGID : LessElementsCentroid);
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
		Integer(tag_layers) = 0;
		Integer(tag_bridge) = NONE;
#endif
		EXIT_FUNC();
	}
	
	void Mesh::RemoveGhostElements(std::vector<Element *> ghost)
	{
		ENTER_FUNC();
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
				if( !it->second.empty() ) qsort(&it->second[0],it->second.size(),sizeof(Element *),have_global_id & mask ? CompareElementsCGID : CompareElementsCCentroid);
				std::vector<Element *> result(ref.size());
				std::vector<Element *>::iterator end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),have_global_id & mask ? LessElementsGID : LessElementsCentroid);
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_ghost.clear();
			for(std::map< int, std::vector<Element *> >::iterator it = del_shared.begin(); it != del_shared.end(); it++)
			{
				std::vector<Element * > & ref = shared_elements[it->first][ElementNum(mask)];
				if( !it->second.empty() ) qsort(&it->second[0],it->second.size(),sizeof(Element *),have_global_id & mask ? CompareElementsCGID : CompareElementsCCentroid);
				std::vector<Element *> result(ref.size());
				std::vector<Element *>::iterator end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),have_global_id & mask ? LessElementsGID : LessElementsCentroid);
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
		EXIT_FUNC();
	}

	
	void Mesh::AssignGlobalID(ElementType mask)
	{
		
		for(ElementType currenttype = NODE; currenttype <= CELL; currenttype = currenttype << 1 )
			if( (currenttype & mask) && !(currenttype & have_global_id ) )
			{
				tag_global_id = CreateTag("GLOBAL_ID",DATA_INTEGER, currenttype, NONE,1);
				have_global_id |= currenttype;
			}
		ENTER_FUNC();

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
					//~ for(Mesh::iteratorElement it = BeginElement(currenttype); it != EndElement(); it++)
					//~ {
						//~ if( it->GetStatus() == Element::Shared )
							//~ it->Integer(tag_global_id) = local_shift++;
					//~ }
				}
			}
			if( tag_global_id.isValid() )
				ExchangeData(tag_global_id,mask);
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
		EXIT_FUNC();
	}
		
	void Mesh::ComputeSharedProcs()
	{
		ENTER_FUNC();
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
		
		REPORT_VAL("processors",procs.size());
#endif
		EXIT_FUNC();
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
		ENTER_FUNC();
#if defined(USE_MPI)
		if( !(have_global_id & CELL) ) AssignGlobalID(CELL);

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
			MarkerType busy = CreateMarker();
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
		EXIT_FUNC();		
		return bridge;		
#else
		EXIT_FUNC();
		return std::map< int, std::vector<Element *> >();
#endif
	}
	
	
	void Mesh::PackTagData(Tag tag, std::vector< std::vector<Element *> > & elements, ElementType mask, std::vector<INMOST_DATA_BULK_TYPE> & buffer)
	{
		if( tag.GetDataType() == DATA_REFERENCE ) return; //NOT IMPLEMENTED
		ENTER_FUNC();
#if defined(USE_MPI)
		ElementType pack_types[2] = {NONE,NONE};
		std::vector<Element *>::iterator eit;
		std::vector<INMOST_DATA_BULK_TYPE> array_data_send;
		std::vector<INMOST_DATA_ENUM_TYPE> array_size_send(2);
		array_data_send.reserve(4096);
		array_size_send.reserve(4096);
		unsigned int size = tag.GetSize();
		for(int i = 0; i < 4; i++) if( (mask & (1 << i)) && tag.isDefined(1 << i) )
		{
			pack_types[0] |= 1 << i;
			REPORT_VAL(ElementTypeName(mask & (1<<i)),elements[i].size());
			if( tag.isSparse(1 << i) )
			{
				pack_types[1] |= 1 << i;
				unsigned int count = array_size_send.size();
				array_size_send.push_back(0);
				for(eit = elements[i].begin(); eit != elements[i].end(); eit++)
					if( (*eit)->HaveData(tag) )
					{
						array_size_send.push_back(eit-elements[i].begin());
						array_size_send[count]++;
						INMOST_DATA_ENUM_TYPE s = (*eit)->GetDataSize(tag);
						INMOST_DATA_ENUM_TYPE had_s = array_data_send.size();
						array_data_send.resize(had_s+s*tag.GetBytesSize());
						(*eit)->GetData(tag,0,s,&array_data_send[had_s]);
						if( size == ENUMUNDEF ) array_size_send.push_back(s);						
					}
			}
			else
			{
				for(eit = elements[i].begin(); eit != elements[i].end(); eit++)
				{
					INMOST_DATA_ENUM_TYPE s = (*eit)->GetDataSize(tag);
					INMOST_DATA_ENUM_TYPE had_s = array_data_send.size();
					array_data_send.resize(had_s+s*tag.GetBytesSize());
					(*eit)->GetData(tag,0,s,&array_data_send[had_s]);
					if( size == ENUMUNDEF ) array_size_send.push_back(s);
				}
			}
		}
		array_size_send[0] = array_size_send.size()-2;
		array_size_send[1] = array_data_send.size();
		REPORT_VAL("size_size",array_size_send[0]);
		REPORT_VAL("data_size",array_size_send[1]);
		int buffer_size = 0,position = buffer.size(),temp;
		MPI_Pack_size(2,INMOST_MPI_DATA_BULK_TYPE,comm,&temp);
		buffer_size+= temp;
		MPI_Pack_size(array_size_send.size(),INMOST_MPI_DATA_ENUM_TYPE,comm,&temp);
		buffer_size+= temp;
		MPI_Pack_size(array_data_send.size()/tag.GetBytesSize(),tag.GetBulkDataType(),comm,&temp);
		buffer_size+= temp;
		buffer.resize(position+buffer_size);
		MPI_Pack(pack_types,2,INMOST_MPI_DATA_BULK_TYPE,&buffer[0],buffer.size(),&position,comm);
		if( !array_size_send.empty() ) MPI_Pack(&array_size_send[0],array_size_send.size(),INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],buffer.size(),&position,comm);
		if( !array_data_send.empty() ) MPI_Pack(&array_data_send[0],array_data_send.size()/tag.GetBytesSize(),tag.GetBulkDataType(),&buffer[0],buffer.size(),&position,comm);
		buffer.resize(position);		
#endif
		REPORT_VAL("TagName",tag.GetTagName());
		EXIT_FUNC();
	}
	
	void Mesh::UnpackTagData(Tag tag, std::vector< std::vector<Element *> > & elements, ElementType mask, std::vector<INMOST_DATA_BULK_TYPE> & buffer, int & position, void (*Operation)(Tag tag, Element * s,INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size))
	{
		if( tag.GetDataType() == DATA_REFERENCE ) return; //NOT IMPLEMENTED
		ENTER_FUNC();
		REPORT_VAL("TagName",tag.GetTagName());
		
#if defined(USE_MPI)
		if( !buffer.empty() )
		{
			int pos = 0, k = 0;
			ElementType recv_mask[2];
			INMOST_DATA_ENUM_TYPE data_recv, size_recv;
			std::vector<Element *>::iterator eit;
			unsigned int size = tag.GetSize();
			std::vector<INMOST_DATA_BULK_TYPE> array_data_recv;
			std::vector<INMOST_DATA_ENUM_TYPE> array_size_recv;
			MPI_Unpack(&buffer[0],1,&position,recv_mask,2,INMOST_MPI_DATA_BULK_TYPE,comm);
			MPI_Unpack(&buffer[0],buffer.size(),&position,&size_recv,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			MPI_Unpack(&buffer[0],buffer.size(),&position,&data_recv,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			array_size_recv.resize(size_recv);
			array_data_recv.resize(data_recv);
			if( !array_size_recv.empty() ) MPI_Unpack(&buffer[0],buffer.size(),&position,&array_size_recv[0],array_size_recv.size(),INMOST_MPI_DATA_ENUM_TYPE,comm);
			if( !array_data_recv.empty() ) MPI_Unpack(&buffer[0],buffer.size(),&position,&array_data_recv[0],array_data_recv.size()/tag.GetBytesSize(),tag.GetBulkDataType(),comm);
			for(int i = 0; i < 4; i++) if( (recv_mask[0] & (1<<i)) )//&& tag.isDefined(1 << i) )
			{
				REPORT_VAL("etype",ElementTypeName(recv_mask[0] & (1<<i)));
				REPORT_VAL("elements",elements[i].size());
				if( !tag.isDefined(1 << i) ) 
				{
					tag = CreateTag(tag.GetTagName(),tag.GetDataType(),1<<i,recv_mask[1] & (1<<i),size);
				}

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
		}
#endif
		EXIT_FUNC();
	}
	
	
	void Mesh::ExchangeDataInnerBegin(std::vector<Tag> tags, parallel_storage & from, parallel_storage & to, ElementType mask, exchange_data & storage)
	{
		/*
		{ // checks for bad input
			if( mask == NONE ) return;
			std::vector<Tag>::iterator it = tags.begin();
			while(it != tags.end() )
				if( !it->isValid() ) it = tags.erase(it);
				else ++it;
			if( tags.empty() ) return;
		}
		*/
		if( m_state == Serial ) return;
		ENTER_FUNC();
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
				{
					REPORT_VAL("processor",*p);
					PackTagData(tags[k],from[*p].get_container(),mask,storage.send_buffers[num_send].second);
				}
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
		EXIT_FUNC();
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
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
		std::vector<int> done;
		while( !(done = FinishRequests(storage.recv_reqs)).empty() )
		{
			for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
			{
				int position = 0;
				for(unsigned int k = 0; k < tags.size(); k++)
				{
					REPORT_VAL("processor",storage.recv_buffers[*qt].first);
					UnpackTagData(tags[k],to[storage.recv_buffers[*qt].first].get_container(),mask,storage.recv_buffers[*qt].second,position,Operation);
				}
			}
		}
		if( !storage.send_reqs.empty() )
		{
			REPORT_MPI(MPI_Waitall(storage.send_reqs.size(),&storage.send_reqs[0],MPI_STATUSES_IGNORE));
		}
#endif
		EXIT_FUNC();
	}
	

	void Mesh::GatherParallelStorage(parallel_storage & ghost, parallel_storage & shared, ElementType mask)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
#if defined(USE_PARALLEL_STORAGE)
		int mpirank = GetProcessorRank();
		REPORT_STR("gather ghost and shared elements")
		double time = Timer();
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
		time = Timer() - time;
		REPORT_VAL("time",time);
		REPORT_STR("sort shared elements")
		time = Timer();
		for(parallel_storage::iterator it = shared.begin(); it != shared.end(); it++)
		{
			REPORT_VAL("processor",it->first);
			for(int i = 0; i < 4; i++) if( mask & (1 << i) )
			{
				if( !it->second[i].empty() )qsort(&it->second[i][0],it->second[i].size(),sizeof(Element *),(have_global_id & (1<<i))? CompareElementsCGID : CompareElementsCCentroid);
				REPORT_VAL(ElementTypeName(mask & (1 << i)),it->second[i].size());
			}
		}
		time = Timer() - time;
		REPORT_VAL("time",time);
		REPORT_STR("sort ghost elements")
		time = Timer();
		for(parallel_storage::iterator it = ghost.begin(); it != ghost.end(); it++)
		{
			REPORT_VAL("processor",it->first);
			for(int i = 0; i < 4; i++) if( mask & (1 << i) )
			{
				if( !it->second[i].empty() ) qsort(&it->second[i][0],it->second[i].size(),sizeof(Element *),(have_global_id & (1<<i))? CompareElementsCGID : CompareElementsCCentroid);	
				REPORT_VAL(ElementTypeName(mask & (1 << i)),it->second[i].size());
			}
		}
		time = Timer() - time;
		REPORT_VAL("time",time);
#endif
#endif
		EXIT_FUNC();
	}
	
	void Mesh::ExchangeData(Tag tag, ElementType mask)
	{
		ExchangeData(std::vector<Tag>(1,tag),mask);
	}
	
	void Mesh::ExchangeData(std::vector<Tag> tags, ElementType mask)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
		exchange_data storage;
		ExchangeDataInnerBegin(tags,shared_elements,ghost_elements,mask,storage);
		ExchangeDataInnerEnd(tags,shared_elements,ghost_elements,mask,DefaultUnpack,storage);
#endif
		EXIT_FUNC();
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
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
		ExchangeDataInnerBegin(tags,shared_elements,ghost_elements,mask,storage);
#endif
		EXIT_FUNC();
	}
	
	void Mesh::ExchangeDataEnd(std::vector<Tag> tags, ElementType mask, exchange_data & storage)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
		ExchangeDataInnerEnd(tags,shared_elements,ghost_elements,mask,DefaultUnpack,storage);
#endif
		EXIT_FUNC();
	}
	
	void Mesh::ReduceData(Tag tag, ElementType mask, void (*Operation)(Tag tag,Element * element,INMOST_DATA_BULK_TYPE * recv_data,INMOST_DATA_ENUM_TYPE recv_size))
	{
		ReduceData(std::vector<Tag>(1,tag),mask,Operation);
	}
	
	void Mesh::ReduceData(std::vector<Tag> tags, ElementType mask, void (*Operation)(Tag tag,Element * element,INMOST_DATA_BULK_TYPE * recv_data,INMOST_DATA_ENUM_TYPE recv_size))
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
		exchange_data storage;
		ExchangeDataInnerBegin(tags,ghost_elements,shared_elements,mask,storage);
		ExchangeDataInnerEnd(tags,ghost_elements,shared_elements,mask,Operation,storage);
#endif
		EXIT_FUNC();
	}
	
	void Mesh::ReduceDataBegin(Tag tag, ElementType mask, exchange_data & storage)
	{
		ReduceDataBegin(std::vector<Tag>(1,tag),mask,storage);
	}
	
	void Mesh::ReduceDataBegin(std::vector<Tag> tags, ElementType mask, exchange_data & storage)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
		ExchangeDataInnerBegin(tags,ghost_elements,shared_elements,mask,storage);
#endif
		EXIT_FUNC();
	}
	
	void Mesh::ReduceDataEnd(Tag tag, ElementType mask, void (*Operation)(Tag tag,Element * element,INMOST_DATA_BULK_TYPE * recv_data,INMOST_DATA_ENUM_TYPE recv_size),exchange_data & storage)
	{
		ReduceDataEnd(std::vector<Tag>(1,tag),mask,Operation,storage);
	}
	
	void Mesh::ReduceDataEnd(std::vector<Tag> tags, ElementType mask, void (*Operation)(Tag tag,Element * element,INMOST_DATA_BULK_TYPE * recv_data,INMOST_DATA_ENUM_TYPE recv_size),exchange_data & storage)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif
		ExchangeDataInnerEnd(tags,ghost_elements,shared_elements,mask,Operation,storage);
#endif
		EXIT_FUNC();
	}
	
	void Mesh::PackElementsData(std::vector<Element *> & all, std::vector<INMOST_DATA_BULK_TYPE> & buffer, int destination, std::vector<std::string> tag_list)
	{
		ENTER_FUNC();
		REPORT_VAL("dest",destination);
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
			MarkerType busy = CreateMarker();
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
			
			REPORT_VAL("NODE",snodes.size());
			REPORT_VAL("EDGE",sedges.size());
			REPORT_VAL("FACE",sfaces.size());
			REPORT_VAL("CELL",scells.size());
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
				if( !coords.empty() ) MPI_Pack(&coords[0],snodes.size()*GetDimensions(),INMOST_MPI_DATA_REAL_TYPE,&buffer[0],buffer.size(),&position,comm);
				if(have_global_id & NODE ) 
					if( !global_ids.empty() ) MPI_Pack(&global_ids[0],snodes.size(),INMOST_MPI_DATA_INTEGER_TYPE,&buffer[0],buffer.size(),&position,comm);
			}
			buffer.resize(position);
		}
		//pack edges
		{
			dynarray<INMOST_DATA_ENUM_TYPE,64> low_conn_size(sedges.size());
			dynarray<Storage::integer,64> low_conn_nums;
			position = buffer.size();
			new_size = 0;
			num = 0; k = 0;
			for(std::vector<Element *>::iterator it = sedges.begin(); it != sedges.end(); it++) 
			{
				//num += (*it)->low_conn.size();
				low_conn_size[k] = 0;//(*it)->low_conn.size();
				for(Element::adj_iterator jt = (*it)->LowConn().begin(); jt != (*it)->LowConn().end(); jt++) if( !(*jt)->Hidden() )
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
			if( !low_conn_size.empty() ) MPI_Pack(low_conn_size.data()	,sedges.size()	,INMOST_MPI_DATA_ENUM_TYPE		,&buffer[0],buffer.size(),&position,comm);
			if( !low_conn_nums.empty() ) MPI_Pack(low_conn_nums.data()	,num			,INMOST_MPI_DATA_INTEGER_TYPE	,&buffer[0],buffer.size(),&position,comm);
			buffer.resize(position);
		}
		//pack faces
		{
			dynarray<INMOST_DATA_ENUM_TYPE,64> low_conn_size(sfaces.size());
			dynarray<Storage::integer,64> low_conn_nums;
			position = buffer.size();
			new_size = 0;
			num = 0; k = 0;
			for(std::vector<Element *>::iterator it = sfaces.begin(); it != sfaces.end(); it++) 
			{
				//num += (*it)->low_conn.size();
				low_conn_size[k] = 0;//(*it)->low_conn.size();
				for(Element::adj_iterator jt = (*it)->LowConn().begin(); jt != (*it)->LowConn().end(); jt++) if( !(*jt)->Hidden() )
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
			if( !low_conn_size.empty() ) MPI_Pack(low_conn_size.data()	,sfaces.size()	,INMOST_MPI_DATA_ENUM_TYPE		,&buffer[0],buffer.size(),&position,comm);
			if( !low_conn_nums.empty() ) MPI_Pack(low_conn_nums.data()	,num			,INMOST_MPI_DATA_INTEGER_TYPE	,&buffer[0],buffer.size(),&position,comm);
			buffer.resize(position);
		}
		//pack cells
		{
			dynarray<INMOST_DATA_ENUM_TYPE,64> low_conn_size(scells.size()), high_conn_size(scells.size());
			dynarray<Storage::integer,64> low_conn_nums;
			dynarray<Storage::integer,64> high_conn_nums;
			INMOST_DATA_ENUM_TYPE num_high = 0;
			position = buffer.size();
			new_size = 0;
			num = 0; k = 0;
			for(std::vector<Element *>::iterator it = scells.begin(); it != scells.end(); it++) 
			{
				//num += (*it)->low_conn.size();
				low_conn_size[k] = 0;//(*it)->low_conn.size();				
				for(Element::adj_iterator jt = (*it)->LowConn().begin(); jt != (*it)->LowConn().end(); jt++) if( !(*jt)->Hidden() )
				{
					std::vector<Element *>::iterator find = std::lower_bound(sfaces.begin(),sfaces.end(),(*jt));
					if( find == sfaces.end() || (*find) != (*jt) ) throw Failure;
					low_conn_nums.push_back(static_cast<Storage::integer>(find - sfaces.begin()));
					low_conn_size[k]++;
					num++;
				}
				//num_high += (*it)->high_conn.size();
				high_conn_size[k] = 0;//(*it)->high_conn.size();
				for(Element::adj_iterator jt = (*it)->HighConn().begin(); jt != (*it)->HighConn().end(); jt++) if( !(*jt)->Hidden() )
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
			if( !low_conn_size.empty() ) MPI_Pack(low_conn_size.data()	,scells.size()	,INMOST_MPI_DATA_ENUM_TYPE	,&buffer[0],buffer.size(),&position,comm);
			if( !low_conn_nums.empty() ) MPI_Pack(low_conn_nums.data()	,num			,INMOST_MPI_DATA_INTEGER_TYPE	,&buffer[0],buffer.size(),&position,comm);
			if( !high_conn_size.empty() ) MPI_Pack(high_conn_size.data()	,scells.size()	,INMOST_MPI_DATA_ENUM_TYPE	,&buffer[0],buffer.size(),&position,comm);
			if( !high_conn_nums.empty() ) MPI_Pack(high_conn_nums.data()	,num_high		,INMOST_MPI_DATA_INTEGER_TYPE	,&buffer[0],buffer.size(),&position,comm);
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
		REPORT_VAL("destination",destination);
		EXIT_FUNC();
	}
	
	
	void Mesh::UnpackElementsData(std::vector<Element *> & all, std::vector<INMOST_DATA_BULK_TYPE> & buffer, int source, std::vector<std::string> & tag_recv)
	{
		ENTER_FUNC();
		REPORT_VAL("source",source);
#if defined(USE_MPI)
		int mpirank = GetProcessorRank();
		INMOST_DATA_ENUM_TYPE num, temp;
		INMOST_DATA_ENUM_TYPE shift = 0;
		int position = 0;
		std::vector<Element *> scells, snodes, sedges, sfaces;
		std::vector<std::vector<Element *> > unpack_tags(4);
		std::vector<Node *> old_nodes;
		all.clear();
		double time = Timer();
		
		
		for(nodes_container::iterator it = nodes.begin(); it != nodes.end(); it++) if( *it != NULL )
			old_nodes.push_back(*it);		
		
		if( !old_nodes.empty() )
		{
			if( have_global_id & NODE ) qsort(&old_nodes[0],old_nodes.size(),sizeof(Node *),CompareElementsCGID);
			else qsort(&old_nodes[0],old_nodes.size(),sizeof(Node *),CompareElementsCCentroid);
		}
			
		time = Timer() - time;
		REPORT_STR("gather and sort old nodes");
		REPORT_VAL("time", time);	
		
		time = Timer();
		//sort(&old_nodes[0],old_nodes.size(),sizeof(Node *),CompareElement,SwapPointer,NULL);
		//std::sort(old_nodes.begin(),old_nodes.end(),LessElementsCentroid);
		//unpack nodes
		{
			unsigned int dim = GetDimensions();
			std::vector<Storage::real> coords;
			if( !buffer.empty() ) MPI_Unpack(&buffer[0],buffer.size(),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			coords.resize(num*dim);
			if( !buffer.empty() ) MPI_Unpack(&buffer[0],buffer.size(),&position,&coords[0],num*dim,INMOST_MPI_DATA_REAL_TYPE,comm);
			
			std::vector<Storage::integer> global_ids;
			if( have_global_id & NODE )
			{
				global_ids.resize(num);
				if( !buffer.empty() ) MPI_Unpack(&buffer[0],buffer.size(),&position,&global_ids[0],num,INMOST_MPI_DATA_INTEGER_TYPE,comm);
			}
			
			
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				Node * new_node;
				int find = -1;
				if( !old_nodes.empty() )
				{
					if( have_global_id & NODE ) find = binary_search(&old_nodes[0],old_nodes.size(),sizeof(Node*),CompareGIDSearch,&global_ids[i]);
					else find = binary_search(&old_nodes[0],old_nodes.size(),sizeof(Node*),CompareCoordSearch,&coords[i*dim]);
				}
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
		time = Timer() - time;
		REPORT_STR("unpack nodes");
		REPORT_VAL("time", time);	
		
		time = Timer();
		//unpack edges
		{
			dynarray<Node *,2> e_nodes;
			shift = 0;
			dynarray<INMOST_DATA_ENUM_TYPE,64> low_conn_size;
			dynarray<Storage::integer,64> low_conn_nums;
			if( !buffer.empty() ) MPI_Unpack(&buffer[0],buffer.size(),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			if( num > 0 )
			{
				low_conn_size.resize(num);
				MPI_Unpack(&buffer[0],buffer.size(),&position,low_conn_size.data(),num,INMOST_MPI_DATA_ENUM_TYPE,comm);
				temp = 0;
				for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += low_conn_size[i];
				if( temp > 0 )
				{
					low_conn_nums.resize(temp);
					MPI_Unpack(&buffer[0],buffer.size(),&position,low_conn_nums.data(),temp,INMOST_MPI_DATA_INTEGER_TYPE,comm);
				}
			}
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				e_nodes.resize(low_conn_size[i]);
				for(INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
					e_nodes[j] = static_cast<Node *>(snodes[low_conn_nums[shift+j]]);
				Edge * new_edge = NULL;
				if( !e_nodes.empty() ) 
					new_edge = static_cast<Edge *>(FindSharedAdjacency(reinterpret_cast<Element **>(&e_nodes[0]),e_nodes.size()));
				if( new_edge == NULL )
				{
					new_edge = CreateEdge(&e_nodes[0],e_nodes.size()).first;
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
		time = Timer() - time;
		REPORT_STR("unpack edges");
		REPORT_VAL("time", time);	
		//unpack faces
		{
			dynarray<Edge *,64> f_edges;
			shift = 0;
			dynarray<INMOST_DATA_ENUM_TYPE,64> low_conn_size;
			dynarray<Storage::integer,128> low_conn_nums;
			if( !buffer.empty() ) MPI_Unpack(&buffer[0],buffer.size(),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			if( num > 0 )
			{
				low_conn_size.resize(num);
				if( !buffer.empty() && !low_conn_size.empty() ) MPI_Unpack(&buffer[0],buffer.size(),&position,low_conn_size.data(),num,INMOST_MPI_DATA_ENUM_TYPE,comm);
				temp = 0;
				for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += low_conn_size[i];
				if( temp > 0 )
				{
					low_conn_nums.resize(temp);
					MPI_Unpack(&buffer[0],buffer.size(),&position,low_conn_nums.data(),temp,INMOST_MPI_DATA_INTEGER_TYPE,comm);
				}
			}
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				f_edges.resize(low_conn_size[i]);
				for(INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
					f_edges[j] = static_cast<Edge *>(sedges[low_conn_nums[shift+j]]);
				Face * new_face = NULL;
				if( !f_edges.empty() ) 
					new_face = static_cast<Face *>(FindSharedAdjacency(reinterpret_cast<Element **>(&f_edges[0]),f_edges.size()));
				if( new_face == NULL )
				{
					new_face = CreateFace(&f_edges[0], f_edges.size()).first;
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
		time = Timer();
		//unpack cells
		{
			dynarray<Face *,64> c_faces;
			dynarray<Node *,128> c_nodes;
			dynarray<INMOST_DATA_ENUM_TYPE,64> low_conn_size;
			dynarray<INMOST_DATA_ENUM_TYPE,128> high_conn_size;
			dynarray<Storage::integer,128> low_conn_nums;
			dynarray<Storage::integer,256> high_conn_nums;
			INMOST_DATA_ENUM_TYPE shift_high = 0;
			shift = 0;
			MPI_Unpack(&buffer[0],buffer.size(),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			if( num > 0 )
			{
				low_conn_size.resize(num);
				MPI_Unpack(&buffer[0],buffer.size(),&position,low_conn_size.data(),num,INMOST_MPI_DATA_ENUM_TYPE,comm);
				temp = 0;
				for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += low_conn_size[i];
				if( temp > 0 )
				{
					low_conn_nums.resize(temp);
					MPI_Unpack(&buffer[0],buffer.size(),&position,low_conn_nums.data(),temp,INMOST_MPI_DATA_INTEGER_TYPE,comm);
				}
				high_conn_size.resize(num);
				MPI_Unpack(&buffer[0],buffer.size(),&position,high_conn_size.data(),num,INMOST_MPI_DATA_ENUM_TYPE,comm);
				temp = 0;
				for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += high_conn_size[i];
				if( temp > 0 )
				{
					high_conn_nums.resize(temp);
					MPI_Unpack(&buffer[0],buffer.size(),&position,high_conn_nums.data(),temp,INMOST_MPI_DATA_INTEGER_TYPE,comm);
				}
			}
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				c_faces.resize(low_conn_size[i]);
				for(INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
					c_faces[j] = static_cast<Face *>(sfaces[low_conn_nums[shift+j]]);
				c_nodes.resize(high_conn_size[i]);
				for(INMOST_DATA_ENUM_TYPE j = 0; j < high_conn_size[i]; j++)
					c_nodes[j] = static_cast<Node *>(snodes[high_conn_nums[shift_high+j]]);
				Cell * new_cell = NULL;
				if( !c_faces.empty() ) 
					new_cell = static_cast<Cell *>(FindSharedAdjacency(reinterpret_cast<Element **>(&c_faces[0]),c_faces.size()));
				if( new_cell == NULL )
				{
					new_cell = CreateCell(&c_faces[0], c_faces.size(), &c_nodes[0], c_nodes.size()).first;
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
		time = Timer() - time;
		REPORT_STR("unpack cells");
		REPORT_VAL("time", time);
		
		
		time = Timer();
		all.reserve(snodes.size()+sedges.size()+sfaces.size()+scells.size());
		all.insert(all.end(),scells.begin(),scells.end());
		all.insert(all.end(),sfaces.begin(),sfaces.end());
		all.insert(all.end(),sedges.begin(),sedges.end());
		all.insert(all.end(),snodes.begin(),snodes.end());
		Storage::integer_array proc = IntegerArrayDV(tag_processors);
		Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),source);
		if( ip == proc.end() || (*ip) != source ) proc.insert(ip,source);
		
		time = Timer() - time;
		REPORT_STR("prepare elements to unpack tag data");
		REPORT_VAL("time", time);	
		
		
		time = Timer();
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
		time = Timer() - time;
		REPORT_STR("unpack tag data");
		REPORT_VAL("time", time);	
#endif
		REPORT_VAL("source",source);
		EXIT_FUNC();
	}
	
	
	void Mesh::ExchangeBuffersInner(exch_buffer_type & send_bufs, exch_buffer_type & recv_bufs,
									std::vector<INMOST_MPI_Request> & send_reqs, std::vector<INMOST_MPI_Request> & recv_reqs)
	{
		ENTER_FUNC();
		REPORT_VAL("exchange number", ++num_exchanges);
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
		EXIT_FUNC();
	}
	
	void Mesh::PrepareReceiveInner(Prepare todo, exch_buffer_type & send_bufs, exch_buffer_type & recv_bufs)
	{
		if( parallel_strategy == 0 && todo == UnknownSize ) return; //in this case we know all we need
		ENTER_FUNC();
#if defined(USE_MPI)
/*
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
					recv_bufs.push_back(proc_buffer_type(ii,std::vector<INMOST_DATA_BULK_TYPE>(shared_space[ii]-1))); // this call would be optimized by compiler
		}
#else
*/
		if( todo == UnknownSize )
		{
			if( parallel_strategy != 0 )
			{
				REPORT_VAL("exchange number", ++num_exchanges);
				unsigned i;
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
				if( !recv_bufs.empty() )
				{
					REPORT_MPI(MPI_Waitall(recv_bufs.size(),&reqs[0],MPI_STATUSES_IGNORE));
				}
				for(i = 0; i < recv_bufs.size(); i++) recv_bufs[i].second.resize(send_recv_size[i]);
				if( !send_bufs.empty() )
				{
					REPORT_MPI(MPI_Waitall(send_bufs.size(),&reqs[recv_bufs.size()],MPI_STATUSES_IGNORE));
				}
			}
		}
		else if( todo == UnknownSource )
		{
#if defined(USE_MPI2)
			int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber();
			unsigned i, end = send_bufs.size();
			memset(shared_space,0,sizeof(unsigned)*mpisize); //zero bits where we receive data
			REPORT_MPI(MPI_Win_fence(MPI_MODE_NOPRECEDE,window)); //start exchange session
			for(i = 0; i < end; i++) shared_space[mpisize+i] = send_bufs[i].second.size()+1; //put data to special part of the memory
			for(i = 0; i < end; i++) 
			{
				REPORT_MPI(MPI_Put(&shared_space[mpisize+i],1,MPI_UNSIGNED,send_bufs[i].first,mpirank,1,MPI_UNSIGNED,window)); //request rdma
			}
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
						recv_bufs.push_back(proc_buffer_type(ii,std::vector<INMOST_DATA_BULK_TYPE>(shared_space[ii]-1))); // this call would be optimized by compiler
			}
#else
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
				std::vector< unsigned > sends_dest_and_size(send_bufs.size()*2+1);
				for(int i = 0; i < send_bufs.size(); i++)
				{
					sends_dest_and_size[i*2+0] = send_bufs[i].first;
					sends_dest_and_size[i*2+1] = send_bufs[i].second.size();
				}
				unsigned recvsize = 0;
				int k,j;
				std::vector<int> allsize(mpisize);
				int size = send_bufs.size()*2;
				REPORT_MPI(MPI_Allgather(&size,1,MPI_INT,&allsize[0],1,MPI_INT,comm));
				std::vector<int> displs(mpisize+1,0);
				for(k = 0; k < mpisize; k++)
					recvsize += allsize[k];
				for(k = 1; k < mpisize+1; k++)
					displs[k] = displs[k-1]+allsize[k-1];
				std::vector<unsigned> recvs_dest_and_size(recvsize+1);
				REPORT_MPI(MPI_Allgatherv(&sends_dest_and_size[0],send_bufs.size()*2,MPI_UNSIGNED,&recvs_dest_and_size[0],&allsize[0],&displs[0],MPI_UNSIGNED,comm));
				recv_bufs.clear();
				for(k = 0; k < mpisize; k++)
					for(j = displs[k]; j < displs[k+1]; j+=2)
					{
						//if( mpirank == 0 ) std::cout << "proc " << k << " sends " << recvs_dest_and_size[j] << " size " << recvs_dest_and_size[j+1] << std::endl;
						if(  static_cast<int>(recvs_dest_and_size[j]) == mpirank )
						{
							recv_bufs.push_back(proc_buffer_type(k,std::vector<INMOST_DATA_BULK_TYPE>()));
							recv_bufs.back().second.resize(recvs_dest_and_size[j+1]);
						}
					}

			}
#endif
		}
//~ #endif
#endif
		EXIT_FUNC();
	}
	
	std::vector<int> Mesh::FinishRequests(std::vector<INMOST_MPI_Request> & recv_reqs)
	{
        std::vector<int> ret(recv_reqs.size(),-1);
		ENTER_FUNC();
#if defined(USE_MPI)
		int outcount = 0;
		REPORT_VAL("requests",recv_reqs.size());
		if( recv_reqs.empty() ) 
			ret.clear();
		else
		{
			if( !recv_reqs.empty() )
			{
				REPORT_MPI(MPI_Waitsome(recv_reqs.size(),&recv_reqs[0],&outcount,&ret[0],MPI_STATUSES_IGNORE));
			}
			else outcount = MPI_UNDEFINED;
			if( outcount == MPI_UNDEFINED ) ret.clear();
			else ret.resize(outcount);
		}
#endif
		EXIT_FUNC();
        return ret;
	}
	
	
	void Mesh::ExchangeMarked(enum Action action)
	{
		ENTER_FUNC();
		if( m_state == Serial ) return;
#if defined(USE_MPI)
		INMOST_DATA_BIG_ENUM_TYPE num_wait;
		int mpirank = GetProcessorRank();
		std::vector<MPI_Request> send_reqs, recv_reqs;
		std::vector<std::string> tag_list, tag_list_recv;
		std::vector<std::string> tag_list_empty;
		exch_buffer_type send_bufs;
		exch_buffer_type recv_bufs;
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
					if( *kt != mpirank )
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
					for(unsigned k = 0; k < missing.size(); k++) send_bufs.push_back(proc_buffer_type(missing[k],std::vector<INMOST_DATA_BULK_TYPE>()));
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
				MarkerType delete_marker = CreateMarker();
				
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
			if( !send_reqs.empty() )
			{
				REPORT_MPI(MPI_Waitall(send_reqs.size(),&send_reqs[0],MPI_STATUSES_IGNORE));
			}
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
						for(unsigned k = 0; k < missing.size(); k++) send_bufs.push_back(proc_buffer_type(missing[k],std::vector<INMOST_DATA_BULK_TYPE>()));
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
			if( !send_reqs.empty() )
			{
				REPORT_MPI(MPI_Waitall(num_wait,&send_reqs[0],MPI_STATUSES_IGNORE));
			}
			//Probably now owner should send processors_tag data
			ExchangeData(tag_processors,CELL | FACE | EDGE | NODE);
		}
		ComputeSharedProcs();
		
		if( action == AMigrate ) 
		{
			if( !send_reqs.empty() )
			{
				REPORT_MPI(MPI_Waitall(send_reqs.size(),&send_reqs[0],MPI_STATUSES_IGNORE));
			}
		}
#endif
		EXIT_FUNC();
	}
	
	void Mesh::RecomputeParallelStorage(ElementType mask)
	{
		ENTER_FUNC();
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
		EXIT_FUNC();
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
		if( m_state == Serial ) return;
		ENTER_FUNC();
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
		double time;
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
			time = Timer();
			for(Storage::integer_array::iterator p = procs.begin(); p != procs.end(); p++)
			{
				MarkerType busy = CreateMarker();
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
			time = Timer() - time;
			REPORT_STR("Mark first layer");
			REPORT_VAL("time",time);
		}
		for(Storage::integer k = layers-1; k >= 0; k--)
		{
			ExchangeMarked();
			old_layers.swap(current_layers);
			current_layers.clear();
			if( k > 0 ) 
			{
				time = Timer();
				for(Storage::integer_array::iterator p = procs.begin(); p != procs.end(); p++)
				{
					std::vector<Element *> & ref_cur = current_layers[*p];
					std::vector<Element *> & ref_old = old_layers[*p];
					MarkerType busy = CreateMarker();
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
				time = Timer() - time;
				REPORT_STR("Mark layer");
				REPORT_VAL("layer",k);
				REPORT_VAL("time",time);
			}
		}
		if( delete_ghost )
		{
			time = Timer();
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
			time = Timer() - time;
			REPORT_STR("Select ghost elements to remove");
			REPORT_VAL("time",time);
			RemoveGhostElements(del_ghost);
		}
		DeleteTag(layers_marker);
		//throw NotImplemented;
#endif
		EXIT_FUNC();

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
		if( m_state == Serial ) return;
		ENTER_FUNC();
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
		
		
		double time = Timer();
		
		
		for(iteratorElement it = BeginElement(CELL); it != EndElement(); it++)
			it->IntegerArrayDV(tag_new_processors).push_back(it->Integer(tag_new_owner));
		
		dynarray<Storage::integer,64> result,intersection;
		//determine tag_new_processors for FACEs to calculate shared skin
		for(ElementType mask = FACE; mask >= NODE; mask = mask >> 1)
			for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			{
				Storage::integer_array procs = it->IntegerArrayDV(tag_new_processors);
				adjacent<Element> over = it->getAdjElements(mask << 1);
				determine_my_procs_high(over,tag_new_processors,result,intersection);
				if( !result.empty() )
				{
					it->Integer(tag_new_owner) = result[0];
					procs.replace(procs.begin(),procs.end(),result.begin(),result.end());
				}
				else
				{
					it->Integer(tag_new_owner) = mpirank;
					procs.clear();
					procs.push_back(mpirank);
				}
			}
			
		time = Timer() - time;
		REPORT_STR("Determine new processors");
		REPORT_VAL("time",time);
			
		ExchangeData(tag_new_owner,FACE | EDGE | NODE);
		
		
		
		
		
		time = Timer();

		if( bridge != NONE && layers != 0 )
		{
			ReduceData(tag_new_processors,FACE ,RedistUnpack);
			ExchangeData(tag_new_processors,FACE);
			
			double time2 = Timer();
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
				MarkerType busy = CreateMarker();
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
			time2 = Timer() - time2;
			REPORT_STR("Compute shared skin");
			REPORT_VAL("time",time2);
			
			double time3 = Timer();
			//now mark all cell layers 
			{
				time2 = Timer();
				std::map<int, std::vector<Element *> > current_layers,old_layers;
				MarkerType busy = CreateMarker();
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
				time2 = Timer() - time2;
				REPORT_STR("Mark first layer");
				REPORT_VAL("time",time2);
				
				
				
				for(Storage::integer k = layers-1; k > 0; k--)
				{
					time2 = Timer();
					old_layers.swap(current_layers);
					current_layers.clear();
					for(std::map< int, std::vector<Element *> >::iterator qt = old_layers.begin(); qt != old_layers.end(); qt++)
					{
						MarkerType busy = CreateMarker();
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
					time2 = Timer() - time2;
					REPORT_STR("Mark layer");
					REPORT_VAL("layer",k);
					REPORT_VAL("time",time2);
				}
			}
			
			ReduceData(tag_new_processors,CELL,RedistUnpack);
			ExchangeData(tag_new_processors,CELL);
			
			time3 = Timer() - time3;
			REPORT_STR("Mark all layers");
			REPORT_VAL("time",time3);
			
			time2 = Timer();
			
			for(ElementType mask = FACE; mask >= NODE; mask = mask >> 1)
			{
				for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
				{
					Storage::integer_array procs = it->IntegerArrayDV(tag_new_processors);
					adjacent<Element> over = it->getAdjElements(mask << 1);
					determine_my_procs_high(over,tag_new_processors,result,intersection);
					if( result.empty() ) 
					{
						procs.clear();
						procs.push_back(mpirank);
					}
					else procs.replace(procs.begin(),procs.end(),result.begin(),result.end());
				}
			}
			
			
			
			ReduceData(tag_new_processors,FACE| EDGE| NODE,RedistUnpack);
			ExchangeData(tag_new_processors,FACE|EDGE|NODE);
			time2 = Timer() - time2;
			REPORT_STR("Detect processors");
			REPORT_VAL("time",time2);
		}
		else 
		{
			ReduceData(tag_new_processors,FACE | EDGE | NODE,RedistUnpack);
			ExchangeData(tag_new_processors,FACE | EDGE | NODE);
		}
		
		time = Timer() - time;
		REPORT_STR("Determine new processors for layers");
		REPORT_VAL("time",time);
		
		
		time = Timer();


		
		for(iteratorElement it = BeginElement(CELL | FACE | EDGE | NODE); it != EndElement(); it++)
		{
			Storage::integer_array new_procs = it->IntegerArray(tag_new_processors);
			if( it->IntegerDF(tag_owner) == mpirank ) // deal with my entities, others will deal with theirs
			{
				//compute processors that should have the entity but they have not
				Storage::integer_array old_procs = it->IntegerArrayDV(tag_processors);
				result.resize(new_procs.size());
				dynarray<Storage::integer,64>::iterator end = std::set_difference(new_procs.begin(),new_procs.end(),old_procs.begin(),old_procs.end(),result.begin());
				result.resize(end-result.begin());
				//mark to send entity to processors that don't have it
				Storage::integer_array sendto = it->IntegerArray(tag_sendto);
				sendto.insert(sendto.end(),result.begin(),result.end());
			}
		}
		
		time = Timer() - time;
		REPORT_STR("Determine local entities to send");
		REPORT_VAL("time",time);

		time = Timer();
		ExchangeMarked(AMigrate);
		time = Timer() - time;
		REPORT_STR("Migrate elements");
		REPORT_VAL("time",time);
		
		
		
		time = Timer();
		DeleteTag(tag_new_owner);
		DeleteTag(tag_new_processors);
		time = Timer() - time;
		REPORT_STR("Delete auxilarry tags");
		REPORT_VAL("time",time);
		//throw NotImplemented;
#endif
		EXIT_FUNC();

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
