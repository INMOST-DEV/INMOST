#ifdef _MSC_VER //kill some warnings
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#endif

#include "inmost.h"


#if defined(USE_MESH)
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

#if defined(USE_MPI)
static INMOST_DATA_BIG_ENUM_TYPE pmid = 0;
#endif

#if defined(USE_PARALLEL_WRITE_TIME)
#define REPORT_MPI(x) {WriteTab(out_time) << "<MPI><![CDATA[" << #x << "]]></MPI>" << std::endl; x;}
#define REPORT_STR(x) {WriteTab(out_time) << "<TEXT><![CDATA[" << x << "]]></TEXT>" << std::endl;}
#define REPORT_VAL(str,x) {WriteTab(out_time) << "<VALUE name=\"" << str << "\"> <CONTENT><![CDATA[" << x << "]]></CONTENT> <CODE><![CDATA[" << #x << "]]></CODE></VALUE>" << std::endl;}
#define ENTER_FUNC() double all_time = Timer(); {WriteTab(out_time) << "<FUNCTION name=\"" << __FUNCTION__ << "\" id=\"func" << func_id++ << "\">" << std::endl; Enter();}
#define EXIT_FUNC() {WriteTab(out_time) << "<TIME>" << Timer() - all_time << "</TIME>" << std::endl; Exit(); WriteTab(out_time) << "</FUNCTION>" << std::endl;}
#define EXIT_FUNC_DIE() {WriteTab(out_time) << "<TIME>" << -1 << "</TIME>" << std::endl; Exit(); WriteTab(out_time) << "</FUNCTION>" << std::endl;}
#else
#define REPORT_MPI(x) x
#define REPORT_STR(x) {}
#define REPORT_VAL(str,x) {}
#define ENTER_FUNC() {}
#define EXIT_FUNC() {}
#define EXIT_FUNC_DIE()  {}
#endif


#if defined(__LINUX__) || defined(__linux__) || defined(__APPLE__)
#include <unistd.h>
#define PROCESSID getpid()
#else
#define PROCESSID -1
#endif

namespace INMOST
{
    static int block_recursion = 0;

    std::string ro()
    {
        int rank = 0;
#ifdef USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
        std::stringstream ss;
        for (int i = 0; i < rank; i++)
            ss << "   ";
        return ss.str();
    }
    
	//////////////////////////////
	/// REDUCTION FUNCTIONS    ///
	//////////////////////////////
	void DefaultUnpack(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
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
					assert(false);
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
				assert(false);
			}
		}
		element->SetData(tag,0,size,data);
	}


	void UnpackOnSkin(const Tag & tag, const Element & e, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size )
		{
      const Storage::integer * recv = static_cast<const Storage::integer *>(static_cast<const void *>(data));
			Storage::integer_array arr = e->IntegerArray(tag);
			arr.push_back(recv[0]);
		}
	}


	void UnpackSkin(const Tag & tag, const Element & e, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size )
		{
      for(INMOST_DATA_ENUM_TYPE k = 0; k < size; k+=2 )
      {
			  bool flag = true;
			  const Storage::integer * recv = static_cast<const Storage::integer *>(static_cast<const void *>(data))+k;
			  Storage::integer_array arr = e->IntegerArray(tag);
      
        for(Storage::integer_array::iterator it = arr.begin(); it != arr.end(); it+=2)
        {
          if( *it == recv[0] )
				  {
					  flag = false;
					  break;
				  }
        }
			  if( flag ) 
			  {
				  arr.push_back(recv[0]);
				  arr.push_back(recv[1]);
			  }
      }
		}
	}
	


	void DeleteUnpack(const Tag & tag, const Element & e, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
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
	
	void RedistUnpack(const Tag & tag,const Element & e, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size ) 
		{
			Storage::integer_array p1 = e->IntegerArray(tag);
			const Storage::integer * p2 = static_cast<const Storage::integer *>(static_cast<const void *>(data));
			dynarray<Storage::integer,64> result(p1.size()+size);
			dynarray<Storage::integer,64>::iterator end;
			end = std::set_union(p1.begin(),p1.end(),p2,p2+size,result.begin());
			result.resize(end-result.begin());
			p1.clear();
			p1.insert(p1.end(),result.begin(),result.end());
		}
	}

	void UnpackLayersMarker(const Tag & tag, const Element & e, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
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
	

	void UnpackSyncMarkerOR(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		(void) size;
		element->Bulk(tag) |= *data;
	}

	void UnpackSyncMarkerXOR(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		(void) size;
		element->Bulk(tag) ^= *data;
	}

	
	void UnpackSyncMarkerAND(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		(void) size;
		element->Bulk(tag) &= *data;
	}

	//////////////////////////////
	/// REDUCTION FUNCTIONS END///
	//////////////////////////////
	bool Mesh::HaveGlobalID(ElementType type) const
	{
		assert(OneType(type));
		Tag test;
		if( GlobalIDTag().isValid() ) test = GlobalIDTag();
		//else if( HaveTag("GLOBAL_ID") ) test = GetTag("GLOBAL_ID") );
		if( test.isValid() )
		{
			//bool ret = true;
			//for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype) ) if( etype & types )
			//	ret &= test.isDefined(etype);
			//return ret;
			return test.isDefined(type);
		}
		else return false;
	}


	void Mesh::SynchronizeMarker(MarkerType marker, ElementType mask, SyncBitOp op)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		if( m_state == Mesh::Parallel )
		{
			Tag t = CreateTag("TEMP_SYNC_MARKER",DATA_BULK,mask,mask,1);
			
			//workaround for old gcc compiler
			const Element::Status SGhost = Element::Ghost;
			const Element::Status SAny = Element::Any;
			Element::Status Expr = (Element::Shared | ((op != SYNC_BIT_SET) ? SGhost : SAny));
			
			for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); ++it)
				if( it->GetMarker(marker) && (it->GetStatus() & Expr) )
					it->Bulk(t) = 1;
			
			
			
			switch(op)
			{
				case SYNC_BIT_SET:
					ExchangeData(t,mask,0);
					for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); ++it)
					{
						if( it->GetStatus() == Element::Ghost )
						{
							if( it->HaveData(t) ) it->SetMarker(marker); else it->RemMarker(marker);
						}
					}
					break;
				case SYNC_BIT_OR:
					ReduceData(t,mask,0,UnpackSyncMarkerOR);
					for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); ++it)
					{
						if( it->GetStatus() & (Element::Ghost | Element::Shared) )
						{
							if( !it->GetMarker(marker) && it->HaveData(t) ) it->SetMarker(marker);
						}
					}
					break;
				case SYNC_BIT_AND:
					ReduceData(t,mask,0,UnpackSyncMarkerAND);
					for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); ++it)
					{
						if( it->GetStatus() & (Element::Ghost | Element::Shared))
						{
							if( it->HaveData(t) && it->Bulk(t) ) it->SetMarker(marker); else it->RemMarker(marker);
						}
					}
					break;
				case SYNC_BIT_XOR:
					ReduceData(t,mask,0,UnpackSyncMarkerXOR);
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
#else//USE_MPI
		(void) marker;
		(void) mask;
		(void) op;
#endif//USE_MPI
		EXIT_FUNC();
	}
	
	ElementType Mesh::SynchronizeElementType(ElementType etype)
	{
		ENTER_FUNC();
		ElementType etypeout = etype;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Allreduce(&etype,&etypeout,1,INMOST_MPI_DATA_BULK_TYPE,MPI_BOR,GetCommunicator()));
#endif//USE_MPI
		return etypeout;
	}
	Storage::integer Mesh::TotalNumberOf(ElementType mask)
	{
		ENTER_FUNC();
		Storage::integer number = 0, ret = 0;
		for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			if( it->GetStatus() != Element::Ghost ) number++;
		ret = number;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Allreduce(&number,&ret,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm));
#endif//USE_MPI
		EXIT_FUNC();
		return ret;
	}
	Storage::integer Mesh::EnumerateSet(const ElementSet & set, const Tag & num_tag, Storage::integer start, bool define_sparse)
	{
		ENTER_FUNC();
		Storage::integer shift = 0, ret = 0;
		ElementType mask = CELL | FACE | EDGE | NODE;
#if defined(USE_MPI)
		Storage::integer number = 0;
		for(ElementSet::iterator it = set.Begin(); it != set.End(); it++)
			if( it->GetStatus() != Element::Ghost && (define_sparse || it->HaveData(num_tag)) )
				number++;
		REPORT_MPI(MPI_Scan(&number,&shift,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm));
		shift -= number;
#endif//USE_MPI
		shift += start;
		for(ElementSet::iterator it = set.Begin(); it != set.End(); it++)
			if( it->GetStatus() == Element::Owned && (define_sparse || it->HaveData(num_tag)) )
				it->Integer(num_tag) = shift++;
		for(ElementSet::iterator it = set.Begin(); it != set.End(); it++)
			if( it->GetStatus() == Element::Shared && (define_sparse || it->HaveData(num_tag)) )
				it->Integer(num_tag) = shift++;
		ExchangeData(num_tag,mask,0);
		ret = shift;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Bcast(&ret,1,INMOST_MPI_DATA_INTEGER_TYPE,GetProcessorsNumber()-1,comm));
		//MPI_Allreduce(&shift,&ret,1,INMOST_DATA_INTEGER_TYPE,MPI_MAX,comm);
#endif//USE_MPI
		EXIT_FUNC();
		return ret;
	}
	Storage::integer Mesh::Enumerate(const HandleType * set, enumerator n, const Tag & num_tag, Storage::integer start, bool define_sparse)
	{
		ENTER_FUNC();
		Storage::integer shift = 0, ret = 0;
		ElementType mask = CELL | FACE | EDGE | NODE;
#if defined(USE_MPI)
		Storage::integer number = 0;
		for(const HandleType * it = set; it != set+n; ++it)
			if( GetStatus(*it) != Element::Ghost && (define_sparse || HaveData(*it,num_tag))) number++;
		REPORT_MPI(MPI_Scan(&number,&shift,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm));
		shift -= number;
#endif//USE_MPI
		shift += start;
		for(const HandleType * it = set; it != set+n; ++it)
			if( GetStatus(*it) == Element::Owned && (define_sparse || HaveData(*it,num_tag))) Integer(*it,num_tag) = shift++;
		for(const HandleType * it = set; it != set+n; ++it)
			if( GetStatus(*it) == Element::Shared && (define_sparse || HaveData(*it,num_tag))) Integer(*it,num_tag) = shift++;
		ExchangeData(num_tag,mask,0);
		ret = shift;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Bcast(&ret,1,INMOST_MPI_DATA_INTEGER_TYPE,GetProcessorsNumber()-1,comm));
		//MPI_Allreduce(&shift,&ret,1,INMOST_DATA_INTEGER_TYPE,MPI_MAX,comm);
#endif//USE_MPI
		EXIT_FUNC();
		return ret;
	}
	Storage::integer Mesh::Enumerate(ElementType mask, Tag num_tag, Storage::integer start, bool define_sparse)
	{
		ENTER_FUNC();
		Storage::integer shift = 0, ret = 0;
#if defined(USE_MPI)
		Storage::integer number = 0;
		for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			if( it->GetStatus() != Element::Ghost && (define_sparse || it->HaveData(num_tag)) )
				number++;
		REPORT_MPI(MPI_Scan(&number,&shift,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm));
		shift -= number;
#endif//USE_MPI
		shift += start;
		for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			if( it->GetStatus() == Element::Owned && (define_sparse || it->HaveData(num_tag)) )
				it->Integer(num_tag) = shift++;
		for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			if( it->GetStatus() == Element::Shared && (define_sparse || it->HaveData(num_tag)) )
				it->Integer(num_tag) = shift++;
		ExchangeData(num_tag,mask,0);
		ret = shift;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Bcast(&ret,1,INMOST_MPI_DATA_INTEGER_TYPE,GetProcessorsNumber()-1,comm));
		//MPI_Allreduce(&shift,&ret,1,INMOST_DATA_INTEGER_TYPE,MPI_MAX,comm);
#endif//USE_MPI
		EXIT_FUNC();
		return ret;
	}
	
	Storage::real Mesh::Integrate(Storage::real input)
	{
		ENTER_FUNC();
		Storage::real output = input;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_REAL_TYPE,MPI_SUM,comm));
#else//USE_MPI
		(void) input;
#endif//USE_MPI
		EXIT_FUNC();
		return output;
	}
	
	
	
	
	Storage::integer Mesh::Integrate(Storage::integer input)
	{
		ENTER_FUNC();
		Storage::integer output = input;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm));
#else//USE_MPI
		(void) input;
#endif//USE_MPI
		EXIT_FUNC();
		return output;
	}
	
	void Mesh::Integrate(Storage::real * input, Storage::integer size)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		static dynarray<Storage::real,64> temp;
		temp.resize(size);
		memcpy(temp.data(),input,sizeof(Storage::real)*size);
		REPORT_MPI(MPI_Allreduce(temp.data(),input,size,INMOST_MPI_DATA_REAL_TYPE,MPI_SUM,comm));
#else//USE_MPI
		(void) input;
		(void) size;
#endif//USE_MPI
		EXIT_FUNC();
	}
	
	void Mesh::Integrate(Storage::integer * input, Storage::integer size)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		static dynarray<Storage::integer,64> temp;
		temp.resize(size);
		memcpy(temp.data(),input,sizeof(Storage::integer)*size);
		REPORT_MPI(MPI_Allreduce(temp.data(),input,size,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm));
#else//USE_MPI
		(void) input;
		(void) size;
#endif//USE_MPI
		EXIT_FUNC();
	}
	
	Storage::integer Mesh::ExclusiveSum(Storage::integer input)
	{
		ENTER_FUNC();
		Storage::integer output = 0;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Scan(&input,&output,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm));
		output -= input;
#else//USE_MPI
		(void) input;
#endif//USE_MPI
		EXIT_FUNC();
		return output;
	}
	
	Storage::real Mesh::Integrate(const Tag & t, enumerator entry, ElementType mask)
	{
		ENTER_FUNC();
		Storage::real output = 0, input = 0;
		for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			if( GetStatus(*it) != Element::Ghost && HaveData(*it,t) )
			{
				real_array arr = RealArray(*it,t);
				if( arr.size() > entry ) input += arr[entry];
			}
		output = input;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_REAL_TYPE,MPI_SUM,comm));
#endif//USE_MPI
		EXIT_FUNC();
		return output;
	}
	
	Storage::real Mesh::AggregateMax(Storage::real input)
	{
		ENTER_FUNC();
		Storage::real output = input;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_REAL_TYPE,MPI_MAX,comm));
#else //USE_MPI
		(void) input;
#endif //USE_MPI
		EXIT_FUNC();
		return output;
	}
	
	Storage::integer Mesh::AggregateMax(Storage::integer input)
	{
		ENTER_FUNC();
		Storage::integer output = input;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_MAX,comm));
#else //USE_MPI
		(void) input;
#endif //USE_MPI
		EXIT_FUNC();
		return output;
	}
	
	void Mesh::AggregateMax(Storage::real * input, Storage::integer size)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		static dynarray<Storage::real,64> temp;
		temp.resize(size);
		memcpy(temp.data(),input,sizeof(Storage::real)*size);
		REPORT_MPI(MPI_Allreduce(temp.data(),input,size,INMOST_MPI_DATA_REAL_TYPE,MPI_MAX,comm));
#else//USE_MPI
		(void) input;
		(void) size;
#endif//USE_MPI
		EXIT_FUNC();
	}
	
	void Mesh::AggregateMax(Storage::integer * input, Storage::integer size)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		static dynarray<Storage::integer,64> temp;
		temp.resize(size);
		memcpy(temp.data(),input,sizeof(Storage::integer)*size);
		REPORT_MPI(MPI_Allreduce(temp.data(),input,size,INMOST_MPI_DATA_INTEGER_TYPE,MPI_MAX,comm));
#else//USE_MPI
		(void) input;
		(void) size;
#endif//USE_MPI
		EXIT_FUNC();
	}
	
	Storage::real Mesh::AggregateMin(Storage::real input)
	{
		ENTER_FUNC();
		Storage::real output = input;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_REAL_TYPE,MPI_MIN,comm));
#else //USE_MPI
		(void) input;
#endif //USE_MPI
		EXIT_FUNC();
		return output;
	}
	
	Storage::integer Mesh::AggregateMin(Storage::integer input)
	{
		ENTER_FUNC();
		Storage::integer output = input;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_INTEGER_TYPE,MPI_MIN,comm));
#else //USE_MPI
		(void) input;
#endif //USE_MPI
		EXIT_FUNC();
		return output;
	}
	
	void Mesh::AggregateMin(Storage::real * input, Storage::integer size)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		static dynarray<Storage::real,64> temp;
		temp.resize(size);
		memcpy(temp.data(),input,sizeof(Storage::real)*size);
		REPORT_MPI(MPI_Allreduce(temp.data(),input,size,INMOST_MPI_DATA_REAL_TYPE,MPI_MIN,comm));
#else//USE_MPI
		(void) input;
		(void) size;
#endif//USE_MPI
		EXIT_FUNC();
	}
	
	void Mesh::AggregateMin(Storage::integer * input, Storage::integer size)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		static dynarray<Storage::integer,64> temp;
		temp.resize(size);
		memcpy(temp.data(),input,sizeof(Storage::integer)*size);
		REPORT_MPI(MPI_Allreduce(temp.data(),input,size,INMOST_MPI_DATA_INTEGER_TYPE,MPI_MIN,comm));
#else//USE_MPI
		(void) input;
		(void) size;
#endif//USE_MPI
		EXIT_FUNC();
	}
	
	
	INMOST_MPI_Comm Mesh::GetCommunicator() const
	{
#if defined(USE_MPI)
		return comm;
#else //USE_MPI
		return 0;
#endif //USE_MPI
	}
	
	int Mesh::GetProcessorRank() const
	{
#if defined(USE_MPI)
		int rank;
		MPI_Comm_rank(comm,&rank);
		return rank;
#else //USE_MPI
		return 0;
#endif //USE_MPI
	}
	
	int Mesh::GetProcessorsNumber() const
	{
#if defined(USE_MPI)
		int size;
		MPI_Comm_size(comm,&size);
		return size;
#else //USE_MPI
		return 1;
#endif //USE_MPI
	}	
	
	void Mesh::Initialize(int * argc, char *** argv)
	{
#if defined(USE_MPI)
		int test;
		MPI_Initialized(&test);
        if( test == 0 ) MPI_Init(argc,argv);
#else //USE_MPI
		(void) argc;
		(void) argv;
#endif //USE_MPI
#if defined(USE_PARALLEL_WRITE_TIME)
		atexit(Mesh::AtExit);
#endif
	}
	void Mesh::Finalize()
	{
#if defined(USE_MPI)
		int test = 0;
		MPI_Finalized(&test);
		if( !test )
			MPI_Finalize();
#endif //USE_MPI
	}
	
    INMOST_MPI_Group Mesh::GetGroup() const
    {
        INMOST_MPI_Group ret = INMOST_MPI_GROUP_EMPTY;
#if defined(USE_MPI)
        MPI_Comm_group(GetCommunicator(), &ret);
#endif
        return ret;
    }
	
	void Mesh::SetCommunicator(INMOST_MPI_Comm _comm)
	{
		ENTER_FUNC();
		tag_shared = CreateTag("PROTECTED_STATUS",DATA_BULK,  ESET |CELL | FACE | EDGE | NODE,NONE,1);
		tag_owner = CreateTag("OWNER_PROCESSOR",DATA_INTEGER, ESET | CELL | FACE | EDGE | NODE,NONE,1);
		tag_processors = CreateTag("PROCESSORS_LIST",DATA_INTEGER, ESET|  MESH | NODE | EDGE | FACE | CELL,NONE);
		tag_layers = CreateTag("LAYERS",DATA_INTEGER,MESH,NONE,1);
		tag_bridge = CreateTag("BRIDGE",DATA_INTEGER,MESH,NONE,1);
		tag_sendto = CreateTag("PROTECTED_SENDTO",DATA_INTEGER, ESET | CELL | FACE | EDGE | NODE, ESET | CELL | FACE | EDGE | NODE);
		tag_an_id = CreateTag("PROTECTED_ID",DATA_INTEGER, CELL | FACE | EDGE | NODE, CELL | FACE | EDGE | NODE);

#if defined(USE_MPI)
		randomizer = Random();
		
		parallel_strategy = 1;
		parallel_file_strategy = 1;

		

		Mesh::Initialize(NULL,NULL);
		//~ MPI_Comm_dup(_comm,&comm);
		comm = _comm;
		{
			INMOST_DATA_BIG_ENUM_TYPE t = pmid;
			REPORT_MPI(MPI_Allreduce(&t,&parallel_mesh_unique_id,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,MPI_MAX,comm));
			pmid = parallel_mesh_unique_id+1;
		}
		m_state = Mesh::Parallel;


#if defined(USE_MPI_P2P)
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
#endif //USE_MPI_P2P
#else //USE_MPI
		(void) _comm;
#endif //USE_MPI
		EXIT_FUNC();
	}
	
#if defined(USE_PARALLEL_WRITE_TIME)
	void Mesh::Enter() { tab++; }
	void Mesh::Exit() {tab--; }
	std::ostream & Mesh::WriteTab(std::ostream & f)
	{
		for(int i = 0; i < tab; i++)
			f << "   ";
		return f;
	}
	std::fstream & Mesh::GetStream()
	{
		return out_time;
	}
	void Mesh::FinalizeFile()
	{
    //std::stringstream str;
    if( tab > 1 )
		{
      out_time << "<TEXT><![CDATA[Died!]]></TEXT>\n";// << std::endl;
		}
    while(tab > 1)
		{
      out_time << "<TIME>-1</TIME>\n</FUNCTION>\n";// << std::endl; 
      Exit(); 
		}
		out_time << "</Debug>" << std::endl;
    //out_time << str;
    //out_time.flush();
    out_time.close();
	}
#endif //USE_PARALLEL_WRITE_TIME
	
	
	
	void determine_my_procs_low(Mesh * m, HandleType h, dynarray<Storage::integer,64> & result, dynarray<Storage::integer,64> & intersection)
	{
		Element::adj_type const & subelements = m->LowConn(h);
		Element::adj_type::const_iterator i = subelements.begin();
		Storage::integer_array p = m->IntegerArrayDV(*i,m->ProcessorsTag());
		result.clear();
		result.insert(result.end(),p.begin(),p.end());
		i++;
		while(i != subelements.end())
		{
			Storage::integer_array q = m->IntegerArrayDV(*i,m->ProcessorsTag());
			intersection.resize(std::max(static_cast<unsigned>(result.size()),static_cast<unsigned>(q.size())));
			dynarray<Storage::integer,64>::iterator qt = std::set_intersection(result.begin(),result.end(),q.begin(),q.end(),intersection.begin());
			intersection.resize(qt-intersection.begin());
			result.swap(intersection);
			i++;
		}
	}
	
	void determine_my_procs_high(Mesh * m, HandleType h, const Tag & procs, dynarray<Storage::integer,64> & result, dynarray<Storage::integer,64> & intersection)
	{
		Element::adj_type const & overelements = m->HighConn(h);
		if( overelements.empty() ) return;
		Element::adj_type::const_iterator i = overelements.begin();
		result.clear();
		while(i != overelements.end())
		{
			Storage::integer_array q = m->IntegerArrayDV(*i,procs);
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
	
	
	__INLINE bool point_in_bbox(Storage::real * p, Storage::real bbox[6], unsigned dim, Storage::real eps)
	{
		bool ret = true;
		for(unsigned k = 0; k < dim; k++) 
			ret &= (p[k] > bbox[k]-eps && p[k] < bbox[dim+k]+eps);
		return ret;
	}
	
	__INLINE bool compare_coord(Storage::real * a, Storage::real * b, INMOST_DATA_ENUM_TYPE dim, Storage::real eps)
	{
		for(INMOST_DATA_ENUM_TYPE i = 0; i <  dim; i++)
			if( ::fabs(a[i]-b[i]) > eps ) return a[i] <= b[i];
		return true;
	}

	class MappingComparator
	{
	public: bool operator () (const std::pair<int,int> & a, const std::pair<int,int> & b) {return a.first < b.first;}
	};
    
	
	
	void Mesh::ResolveShared(bool only_new)
	{
		ENTER_FUNC();
		
#if defined(USE_MPI)
		if( m_state == Mesh::Serial ) SetCommunicator(INMOST_MPI_COMM_WORLD);
		if( tag_global_id.isValid() ) tag_global_id = DeleteTag(tag_global_id,CELL | EDGE | FACE | NODE);
		integer dim = GetDimensions();
		int sendsize;
		int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber();
        int rank = mpirank;
#if defined(USE_PARALLEL_STORAGE)
		shared_elements.clear();
		ghost_elements.clear();
#endif //USE_PARALLEL_STORAGE
		//determine which bboxes i intersect
		dynarray<int,64> procs;
		Storage::real bbox[6]; //local bounding box
		std::vector<Storage::real> bboxs(mpisize*6);
		//Compute local bounding box containing nodes.
		//Will be more convinient to compute (or store)
		//and communicate local octree over all the nodes.
		for(integer k = 0; k < dim; k++)
		{
			bbox[k] = 1e20;
			bbox[k+dim] = -1e20;
		}
#if defined(USE_OMP)
#pragma omp parallel
#endif
		{
			real bbox0[6];
			for(integer k = 0; k < dim; k++)
			{
				bbox0[k] = 1e20;
				bbox0[k+dim] = -1e20;
			}
#if defined(USE_OMP)
#pragma omp for
#endif
			for(integer nit = 0; nit < NodeLastLocalID(); ++nit) if( isValidNode(nit) )
			{
				Node it = NodeByLocalID(nit);
				Storage::real_array arr = it->Coords();
				for(integer k = 0; k < dim; k++)
				{
					if( arr[k] < bbox0[k] ) bbox0[k] = arr[k];
					if( arr[k] > bbox0[k+dim] ) bbox0[k+dim] = arr[k];
				}
			}
#if defined(USE_OMP)
#pragma omp critical
#endif
			{
				for(integer k = 0; k < dim; k++)
				{
					if( bbox0[k] < bbox[k] ) bbox[k] = bbox0[k];
					if( bbox0[k+dim] > bbox[k+dim] ) bbox[k+dim] = bbox0[k+dim];
				}
			}
		}
		REPORT_VAL("dim",dim);
		assert(dim*2 <= 6);
		// write down bounding boxes
		for(integer k = 0; k < dim; k++)
		{
			REPORT_VAL("min",bbox[k]);
			REPORT_VAL("max",bbox[dim+k]);
		}
		// communicate bounding boxes
		REPORT_MPI(MPI_Allgather(&bbox[0],dim*2,INMOST_MPI_DATA_REAL_TYPE,&bboxs[0],dim*2,INMOST_MPI_DATA_REAL_TYPE,comm));
		// find all processors that i communicate with
		for(int k = 0; k < mpisize; k++)
			if( k != mpirank )
			{
				bool flag = true;
				for(integer q = 0; q < dim; q++)
					flag &= !((bbox[q]-GetEpsilon() > bboxs[k*dim*2+q+dim]) || (bbox[dim+q]+GetEpsilon() < bboxs[k*dim*2+q]));
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
				for(integer j = 0; j < dim*2; j++)
					same_box &= ::fabs(bbox[j] - bboxs[k*dim*2+j]) < GetEpsilon();
				same_boxes &= same_box;
			}
			
			if( same_boxes )
			{
				REPORT_STR("All bounding boxes are the same - assuming that mesh is replicated over all nodes");
				//for(Mesh::iteratorElement it = BeginElement(CELL | EDGE | FACE | NODE); it != EndElement(); it++)
				for(ElementType etype = NODE; etype <= CELL; etype = NextElementType(etype) )
				{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
					for(integer eit = 0; eit < LastLocalID(etype); ++eit) if( isValidElement(etype,eit) )
					{
						Element it = ElementByLocalID(etype,eit);
                        if (only_new && GetMarker(it->GetHandle(),NewMarker()) == false) continue;
                        
						integer_array arr = it->IntegerArrayDV(tag_processors);
						arr.resize(mpisize);
						for(int k = 0; k < mpisize; k++)
							arr[k] = k;
						it->IntegerDF(tag_owner) = 0;
						if( mpirank == 0 )
							SetStatus(it->GetHandle(),Element::Shared);
						else
							SetStatus(it->GetHandle(),Element::Ghost);
					}
				}
				ComputeSharedProcs();
				RecomputeParallelStorage(CELL | EDGE | FACE | NODE);
				AssignGlobalID(CELL | EDGE | FACE | NODE);
#if defined(USE_PARALLEL_STORAGE)
				for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
					for(int i = 0; i < 4; i++)
					{
						if( !it->second[i].empty() )
							std::sort(it->second[i].begin(),it->second[i].end(),GlobalIDComparator(this));
							//qsort(&it->second[i][0],it->second[i].size(),sizeof(Element *),CompareElementsCGID);
					}
				for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)
					for(int i = 0; i < 4; i++)
					{
						if( !it->second[i].empty() )
							std::sort(it->second[i].begin(),it->second[i].end(),GlobalIDComparator(this));
							//qsort(&it->second[i][0],it->second[i].size(),sizeof(Element *),CompareElementsCGID);
					}
#endif //USE_PARALLEL_STORAGE
			}
			else
			{
                
                if (only_new)
                {
		            for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); it++) if (GetMarker(*it,NewMarker()))
                    {
                        ElementArray<Node> nodes = it->getNodes();
                        for (int i = 0; i < nodes.size(); i++)
                        {
                            nodes[i].SetMarker(NewMarker());
                        }
                    }
                }

				double time = Timer();
				Storage::real epsilon = GetEpsilon();
			
				GeomParam table;
				for(ElementType etype = EDGE; etype <= CELL; etype = NextElementType(etype))
					if( !HaveGeometricData(CENTROID,etype) )
						table[CENTROID] |= etype;
				PrepareGeometricData(table);
				
				time = Timer() - time;
			
				REPORT_STR("Prepare geometric data");
				REPORT_VAL("time",time);
			
				//for(iteratorNode it = BeginNode(); it != EndNode(); it++)
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(integer nit = 0; nit < NodeLastLocalID(); ++nit) if( isValidNode(nit) )
				{
					Node it = NodeByLocalID(nit);
					Storage::integer_array arr = it->IntegerArrayDV(tag_processors);
					arr.resize(1);
					arr[0] = mpirank;
				}
			

				buffer_type exch_data;
				std::vector< INMOST_DATA_REAL_TYPE > unpack_real;
				std::vector< INMOST_DATA_REAL_TYPE > pack_real;
				element_set sorted_nodes;
				
				sorted_nodes.reserve(NumberOfNodes());
				
				
				time = Timer();
				
				
				for(iteratorNode n = BeginNode(); n != EndNode(); n++)
				{
					real_array c = n->Coords();
					for(real_array::size_type k = 0; k < procs.size(); k++)
						if( point_in_bbox(c.data(),bboxs.data()+procs[k]*dim*2,dim,GetEpsilon()) )
						{
							sorted_nodes.push_back(*n);
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
					std::sort(sorted_nodes.begin(),sorted_nodes.end(),CentroidComparator(this));
					//qsort(&sorted_nodes[0],sorted_nodes.size(),sizeof(Element *),CompareElementsCCentroid);
				time = Timer() - time;
				REPORT_STR("Sort nodes");
				REPORT_VAL("time",time);
				
				
				pack_real.reserve(sorted_nodes.size()*dim);
				
				
				time = Timer();
				for(element_set::iterator it = sorted_nodes.begin(); it != sorted_nodes.end(); it++)
				{
					Storage::real_array arr = RealArrayDF(*it,CoordsTag());
					pack_real.insert(pack_real.end(),arr.begin(),arr.end());
				}
				time = Timer() - time;
				REPORT_STR("Gather coordinates");
				REPORT_VAL("time",time);
				
				
				time = Timer();
			
				MPI_Pack_size(static_cast<int>(pack_real.size()),INMOST_MPI_DATA_REAL_TYPE,comm,&sendsize);
				exch_data.resize(sendsize);
				int position = 0;
				if( sendsize > 0 ) MPI_Pack(&pack_real[0],static_cast<INMOST_MPI_SIZE>(pack_real.size()),INMOST_MPI_DATA_REAL_TYPE,&exch_data[0],static_cast<INMOST_MPI_SIZE>(exch_data.size()),&position,comm);
				
				
				time = Timer() - time;
				REPORT_STR("Pack coordinates");
				REPORT_VAL("time",time);

				
				
				{
					std::vector< MPI_Request > send_reqs, recv_reqs;
					exch_buffer_type send_buffs(procs.size()), recv_buffs(procs.size());
					std::vector<int> done;
//~ #if defined(USE_MPI_P2P)
					//~ unsigned * sendsizeall = shared_space;
					//~ unsigned usend[2] = {sendsize,pack_real.size()};
					//~ REPORT_MPI(MPI_Win_fence(0,window)); //start exchange session
					//~ for(unsigned k = 0; k < procs.size(); k++)
						//~ REPORT_MPI(MPI_Put(usend,2,MPI_UNSIGNED,procs[k],mpirank*2,2,MPI_UNSIGNED,window));
					//~ REPORT_MPI(MPI_Win_fence(MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED,window)); //end exchange session
//~ #else
					std::vector<unsigned> sendsizeall(mpisize*2);
					int pack_size2 = 0;
					unsigned usend[2] = {static_cast<unsigned>(sendsize),static_cast<unsigned>(pack_real.size())};
					MPI_Pack_size(2,MPI_UNSIGNED,comm,&pack_size2);
					for(dynarray<integer,64>::size_type k = 0; k < procs.size(); k++)
					{
						send_buffs[k].first = procs[k];
						send_buffs[k].second.resize(pack_size2);
						position = 0;
						MPI_Pack(usend,2,MPI_UNSIGNED,&send_buffs[k].second[0],static_cast<INMOST_MPI_SIZE>(send_buffs[k].second.size()),&position,comm);
						recv_buffs[k].first = procs[k];
						recv_buffs[k].second.resize(pack_size2);
					}
					ExchangeBuffersInner(send_buffs,recv_buffs,send_reqs,recv_reqs);
					while( !(done = FinishRequests(recv_reqs)).empty() )
					{
						for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
						{
							position = 0;
							MPI_Unpack(&recv_buffs[*qt].second[0],static_cast<INMOST_MPI_SIZE>(recv_buffs[*qt].second.size()),&position,&sendsizeall[procs[*qt]*2],2,MPI_UNSIGNED,comm);
						}
					}
					if( !send_reqs.empty() )
					{
						REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(send_reqs.size()),&send_reqs[0],MPI_STATUSES_IGNORE));
					}
					//~ REPORT_MPI(MPI_Allgather(usend,2,MPI_UNSIGNED,&sendsizeall[0],2,MPI_UNSIGNED,comm));
//~ #endif
					double time2 = Timer();
					{
						
						
						for(dynarray<integer,64>::size_type k = 0; k < procs.size(); k++)
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
								MPI_Unpack(&recv_buffs[*qt].second[0],static_cast<INMOST_MPI_SIZE>(recv_buffs[*qt].second.size()),&position,&unpack_real[0],static_cast<INMOST_MPI_SIZE>(unpack_real.size()),INMOST_MPI_DATA_REAL_TYPE,comm);
								std::vector<Storage::real>::iterator it1 = pack_real.begin() , it2 = unpack_real.begin();
								while(it1 != pack_real.end() && it2 != unpack_real.end() )
								{
									int res = 0;
									for(integer k = 0; k < dim; k++)
										if( ::fabs((*(it1+k))-(*(it2+k))) > epsilon )
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
										IntegerArrayDV(sorted_nodes[(it1-pack_real.begin())/dim],tag_processors).push_back(recv_buffs[*qt].first);
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
							REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(send_reqs.size()),&send_reqs[0],MPI_STATUSES_IGNORE));
						}
					}			
					
					time2 = Timer() - time2;
					REPORT_STR("Intersect all coordinates");
					REPORT_VAL("time",time2);
					
					time = Timer();
					Element::Status estat;
					for(Mesh::iteratorElement it = BeginElement(NODE); it != EndElement(); it++)
					{
                        if (only_new && GetMarker(*it,NewMarker()) == false) continue;
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
								estat = Element::Owned;
							else 
								estat = Element::Shared;
						}
						else
							estat = Element::Ghost;
						SetStatus(*it,estat);
					}
					
					ComputeSharedProcs();
					RecomputeParallelStorage(NODE);
					AssignGlobalID(NODE);
					
#if defined(USE_PARALLEL_STORAGE)
					for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
					{
						if( !it->second[0].empty() )
							std::sort(it->second[0].begin(),it->second[0].end(),GlobalIDComparator(this));
							//qsort(&it->second[0][0],it->second[0].size(),sizeof(Element *),CompareElementsCGID);
					}
					for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)
					{
						if( !it->second[0].empty() )
							std::sort(it->second[0].begin(),it->second[0].end(),GlobalIDComparator(this));
							//qsort(&it->second[0][0],it->second[0].size(),sizeof(Element *),CompareElementsCGID);
					}
#endif //USE_PARALLEL_STORAGE
					
					
					time = Timer() - time;
					REPORT_STR("Set parallel info for nodes");
					REPORT_VAL("time",time);
					
					
					
					
					for(ElementType current_mask = EDGE; current_mask <= CELL; current_mask = NextElementType(current_mask))
					{
						REPORT_STR("Set parallel info for");
						REPORT_VAL("type",ElementTypeName(current_mask));
						
						
						//int owned_elems = 0;
						//int shared_elems = 0;
						int owner;
						Element::Status estat;
						
						time = Timer();
						//Determine what processors potentially share the element
#if defined(USE_OMP)
#pragma omp parallel
#endif
						{
							dynarray<Storage::integer,64> result, intersection;
#if defined(USE_OMP)
#pragma omp for
#endif
							for(integer eit = 0; eit < LastLocalID(current_mask); ++eit)
							{
								if( isValidElement(current_mask,eit) )
								{
									Element it = ElementByLocalID(current_mask,eit);
                                    if (only_new && GetMarker(it.GetHandle(),NewMarker()) == false) continue;

									determine_my_procs_low(this,it->GetHandle(), result, intersection);
									Storage::integer_array p = it->IntegerArrayDV(tag_processors);
									if( result.empty() )
									{
										p.clear();
										p.push_back(mpirank);
										//++owned_elems;
									}
									else
									{
										p.replace(p.begin(),p.end(),result.begin(),result.end());
										//if( result.size() == 1 && result[0] == mpirank )
										//  ++owned_elems;
										//else ++shared_elems;
									}
								}
							}
						}
						time = Timer() - time;
						//REPORT_VAL("predicted owned elements",owned_elems);
						//REPORT_VAL("predicted shared elements",shared_elems);
						REPORT_STR("Predict processors for elements");
						REPORT_VAL("time",time);
						
						
						time = Timer();
						//Initialize mapping that helps get local id by global id
						std::vector<std::pair<int,int> > mapping;
						REPORT_VAL("mapping type",ElementTypeName(PrevElementType(current_mask)));
						for(Mesh::iteratorElement it = BeginElement(PrevElementType(current_mask)); it != EndElement(); it++)
						{
							mapping.push_back(std::make_pair(it->GlobalID(),it->LocalID()));
						}
						if( !mapping.empty() ) 
							std::sort(mapping.begin(),mapping.end(),MappingComparator());
						/*
#if defined(USE_PARALLEL_WRITE_TIME)
						for(std::vector<std::pair<int,int> >::iterator it = mapping.begin(); it != mapping.end(); ++it)
							REPORT_STR("global " << it->first << " local " << it->second << Element::StatusName(ElementByLocalID(PrevElementType(current_mask),it->second)->GetStatus()));
#endif
						 */
						time = Timer() - time;
						REPORT_VAL("mapping size",mapping.size())
						REPORT_STR("Compute global to local indexes mapping");
						REPORT_VAL("time",time);
						//Initialize arrays
						
						time = Timer();
						Storage::integer_array procs = IntegerArrayDV(GetHandle(),tag_processors);
						Storage::integer_array::iterator p = procs.begin();
						std::vector<int> message_send;
						std::vector< std::vector<int> > message_recv(procs.size());
						std::vector< element_set > elements(procs.size());
						
						std::vector< MPI_Request > send_reqs,recv_reqs;
						exch_buffer_type send_buffs(procs.size()), recv_buffs(procs.size());
						
						//Gather all possible shared elements and send global ids of their connectivity to the neighbouring proccessors
						for(p = procs.begin(); p != procs.end(); p++)
						{
							REPORT_VAL("for processor",*p);
							int m = static_cast<int>(p-procs.begin());
							{
								message_send.clear();
								message_send.push_back(0);
								message_send.push_back(0);
								for(Mesh::iteratorElement it = BeginElement(current_mask); it != EndElement(); it++)
								{
                                    if (only_new && GetMarker(*it,NewMarker()) == false) continue;
									Storage::integer_array pr = it->IntegerArrayDV(tag_processors);
									if( std::binary_search(pr.begin(),pr.end(),*p) )
									{
										Element::adj_type & sub = LowConn(*it);
										if( sub.size() == 0 ) throw Impossible;
										message_send.push_back(static_cast<int>(sub.size()));
										REPORT_VAL("number of connections",sub.size());
										for(Element::adj_type::iterator kt = sub.begin(); kt != sub.end(); kt++)
										{
											message_send.push_back(GlobalID(*kt));
											REPORT_STR("global id " << GlobalID(*kt) << " local id " << GetHandleID(*kt) << " " << Element::StatusName(Element(this,*kt)->GetStatus()));
										}
										message_send[1]++;
										elements[m].push_back(*it);
									}
								}
								
								REPORT_VAL("gathered elements",elements[m].size());
								message_send[0] = static_cast<int>(message_send.size());
								MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(message_send.size()),MPI_INT,comm,&sendsize);
								send_buffs[m].first = *p;
								send_buffs[m].second.resize(sendsize);
								int position = 0;
								MPI_Pack(&message_send[0],static_cast<INMOST_MPI_SIZE>(message_send.size()),MPI_INT,&send_buffs[m].second[0],static_cast<INMOST_MPI_SIZE>(send_buffs[m].second.size()),&position,comm);
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
										pos = static_cast<int>(p - procs.begin());
										break;
									}
								if( pos == -1 ) throw Impossible;
								MPI_Unpack(&recv_buffs[*qt].second[0],static_cast<INMOST_MPI_SIZE>(recv_buffs[*qt].second.size()),&position,&size,1,MPI_INT,comm);
								REPORT_VAL("unpacked message size",size-1);
								message_recv[pos].resize(size-1);
								MPI_Unpack(&recv_buffs[*qt].second[0],static_cast<INMOST_MPI_SIZE>(recv_buffs[*qt].second.size()),&position,&message_recv[pos][0],static_cast<INMOST_MPI_SIZE>(message_recv[pos].size()),MPI_INT,comm);
							}
						}
						
						time = Timer() - time;
						REPORT_STR("Exchange messages");
						REPORT_VAL("time",time);
						
						time = Timer();
						//Now find the difference of local elements with given processor number and remote elements
						for(p = procs.begin(); p != procs.end(); p++)
						{
							REPORT_VAL("on processor",*p);
							int m = static_cast<int>(p-procs.begin());
							element_set remote_elements;
							int pos = 0;
							if( message_recv[m].empty() ) continue;
							int num_remote_elements = message_recv[m][pos++];
							REPORT_VAL("number of remote elements",num_remote_elements);
							for(int i = 0; i < num_remote_elements; i++)
							{
								int conn_size = message_recv[m][pos++], flag = 1;
								REPORT_VAL("number of connections",conn_size);
								dynarray<HandleType,64> sub_elements;
								for(int j = 0; j < conn_size; j++)
								{
									int global_id = message_recv[m][pos++];
									int find = -1;
									std::vector<std::pair<int,int> >::iterator it = std::lower_bound(mapping.begin(),mapping.end(),std::make_pair(global_id,0),MappingComparator());
									if( it != mapping.end() && it->first == global_id) 
										find = static_cast<int>(it-mapping.begin());
									if( find == -1 ) 
									{
										REPORT_STR("global id " << global_id << " local id -1");
										flag = 0;
										pos += conn_size-j-1;
										break;
									}
									int find_local_id = mapping[find].second;
									REPORT_STR("global id " << global_id << " local id " << find_local_id << " " << Element::StatusName(ElementByLocalID(PrevElementType(current_mask), find_local_id)->GetStatus()));
									sub_elements.push_back(ComposeHandle(PrevElementType(current_mask), find_local_id));
									
								}
								REPORT_VAL("is element found",flag);
								if( flag )
								{
									HandleType e = FindSharedAdjacency(sub_elements.data(),static_cast<enumerator>(sub_elements.size()));
									REPORT_VAL("found element",e);
									if( e == InvalidHandle() ) continue;
									remote_elements.push_back(e);
								}
							}
							REPORT_VAL("number of unpacked remote elements",remote_elements.size());
							if( !remote_elements.empty() )
							{
								REPORT_VAL("first",remote_elements.front());
								REPORT_VAL("first type",ElementTypeName(GetHandleElementType(remote_elements.front())));
								REPORT_VAL("last",remote_elements.back());
								REPORT_VAL("last type",ElementTypeName(GetHandleElementType(remote_elements.back())));
							}
							std::sort(remote_elements.begin(),remote_elements.end());
							REPORT_VAL("original elements size",elements[m].size());
							if( !elements[m].empty() )
							{
								REPORT_VAL("first",elements[m].front());
								REPORT_VAL("first type",ElementTypeName(GetHandleElementType(elements[m].front())));
								REPORT_VAL("last",elements[m].back());
								REPORT_VAL("last type",ElementTypeName(GetHandleElementType(elements[m].back())));
							}
							std::sort(elements[m].begin(),elements[m].end());
							element_set result;
							element_set::iterator set_end;
							
							result.resize(elements[m].size());
							set_end = std::set_difference(elements[m].begin(),elements[m].end(),remote_elements.begin(),remote_elements.end(), result.begin());
							result.resize(set_end-result.begin());

							REPORT_VAL("set difference size",result.size());
							
							//elements in result are wrongly marked as ghost
							for(element_set::iterator qt = result.begin(); qt != result.end(); qt++)
							{
								integer_array pr = IntegerArrayDV(*qt,tag_processors);
								integer_array::iterator find = std::lower_bound(pr.begin(),pr.end(),*p);
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
                            if (only_new && (GetMarker(*it,NewMarker()) == false)) continue;
    
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
							SetStatus(*it, estat);
						}
						RecomputeParallelStorage(current_mask);
						if( !send_reqs.empty() )
						{
							REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(send_reqs.size()),&send_reqs[0],MPI_STATUSES_IGNORE));
						}
						AssignGlobalID(current_mask);					
						integer num = ElementNum(current_mask);
#if defined(USE_PARALLEL_STORAGE)
						for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
						{
							if( !it->second[num].empty() )
								std::sort(it->second[num].begin(),it->second[num].end(),GlobalIDComparator(this));
						}
						for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)
						{
							if( !it->second[num].empty() )
								std::sort(it->second[num].begin(),it->second[num].end(),GlobalIDComparator(this));
						}
#endif //USE_PARALLEL_STORAGE
	
						time = Timer() - time;
						REPORT_STR("Set parallel info");
						REPORT_VAL("time",time);
					}
				}
				RemoveGeometricData(table);	
			}
		}
		//RemoveGhost();
		
#else //USE_MPI
		AssignGlobalID(CELL | FACE | EDGE | NODE);
        (void)only_new;
#endif //USE_MPI
		EXIT_FUNC();
	}
	
	
	void Mesh::RemoveGhost()
	{
		if( m_state == Mesh::Serial ) return;
		ENTER_FUNC()
#if defined(USE_MPI)
		element_set delete_elements;
#if defined(USE_PARALLEL_STORAGE)
		proc_elements del_shared, del_ghost;
#endif //USE_PARALLEL_STORAGE
		int mpirank = GetProcessorRank();
		Tag tag_delete = CreateTag("TEMPORARY_DELETE_GHOST_TAG",DATA_INTEGER,FACE | EDGE | NODE,FACE | EDGE | NODE);
		
		for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); it++)
		{
			Element::Status estat = GetStatus(*it);
			if( estat == Element::Ghost ) 
			{
#if defined(USE_PARALLEL_STORAGE)
				del_ghost[it->IntegerDF(tag_owner)].push_back(*it);
#endif //USE_PARALLEL_STORAGE
				delete_elements.push_back(*it);
			}
			else if( estat == Element::Shared )
			{
				SetStatus(*it,Element::Owned);
				Storage::integer_array it_procs = it->IntegerArrayDV(tag_processors);
				
#if defined(USE_PARALLEL_STORAGE)
				for(Storage::integer_array::iterator vit = it_procs.begin(); vit != it_procs.end(); vit++)
					if( *vit != mpirank )
						del_shared[*vit].push_back(*it);
#endif //USE_PARALLEL_STORAGE
				
				it_procs.resize(1);
				it_procs[0] = mpirank;
			}
		}
		
#if defined(USE_PARALLEL_STORAGE)
		for(proc_elements::iterator it = del_ghost.begin(); it != del_ghost.end(); it++)
		{
			element_set & ref = ghost_elements[it->first][ElementNum(CELL)];
			if( !it->second.empty() ) 
			{
				if( HaveGlobalID(CELL) )
					std::sort(it->second.begin(),it->second.end(),GlobalIDComparator(this));
				else
					std::sort(it->second.begin(),it->second.end(),CentroidComparator(this));
			}
			element_set result(ref.size());
			element_set::iterator end;
			if( HaveGlobalID(CELL) )
				end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),GlobalIDComparator(this));
			else
				end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),CentroidComparator(this));
			result.resize(end-result.begin());
			ref.swap(result);
		}
		del_ghost.clear();
		for(proc_elements::iterator it = del_shared.begin(); it != del_shared.end(); it++)
		{
			element_set & ref = shared_elements[it->first][ElementNum(CELL)];
			if( !it->second.empty() ) 
			{
				if( HaveGlobalID(CELL) )
					std::sort(it->second.begin(),it->second.end(),GlobalIDComparator(this));
				else
					std::sort(it->second.begin(),it->second.end(),CentroidComparator(this));
			}
			element_set result(ref.size());
			element_set::iterator end;
			if( HaveGlobalID(CELL) )
				end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),GlobalIDComparator(this));
			else
				end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),CentroidComparator(this));
			result.resize(end-result.begin());
			ref.swap(result);
		}
		del_shared.clear();
#endif //USE_PARALLEL_STORAGE

		for(element_set::iterator it = delete_elements.begin(); it != delete_elements.end(); it++) Destroy(*it);
		
		delete_elements.clear();
		for(ElementType mask = FACE; mask >= NODE; mask = PrevElementType(mask))
		{
			for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			{
				Element::Status estat = GetStatus(*it);
				if( estat == Element::Owned || estat == Element::Shared ) continue;
				if( it->nbAdjElements(NextElementType(mask)) == 0 ) it->Integer(tag_delete) = mpirank;
			}
			ReduceData(tag_delete,mask,0,DeleteUnpack);
			ExchangeData(tag_delete,mask,0);
			for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			{
				Element::Status estat = GetStatus(*it);
				if( estat == Element::Owned ) continue;
				if( estat == Element::Ghost && it->nbAdjElements(NextElementType(mask)) == 0 )
				{
#if defined(USE_PARALLEL_STORAGE)
					del_ghost[it->IntegerDF(tag_owner)].push_back(*it);
#endif //USE_PARALLEL_STORAGE
					delete_elements.push_back(*it);
				}
				else if( it->HaveData(tag_delete) )
				{
					Storage::integer_array procs = it->IntegerArrayDV(tag_processors);
					Storage::integer_array del_procs = it->IntegerArray(tag_delete);
					std::vector<Storage::integer> result(procs.size());
					std::sort(del_procs.begin(),del_procs.end());								
#if defined(USE_PARALLEL_STORAGE)				
					for(Storage::integer_array::iterator vit = del_procs.begin(); vit != del_procs.end(); vit++)
						del_shared[*vit].push_back(*it);
#endif //USE_PARALLEL_STORAGE				
					std::vector<Storage::integer>::iterator end = std::set_difference(procs.begin(),procs.end(),del_procs.begin(),del_procs.end(),result.begin());
					result.resize(end-result.begin());
					procs.clear();
					procs.insert(procs.begin(),result.begin(),result.end());
					
					if( procs.size() == 1 && procs[0] == mpirank )
						SetStatus(*it,Element::Owned);
				}
			}
			
#if defined(USE_PARALLEL_STORAGE)
			for(proc_elements::iterator it = del_ghost.begin(); it != del_ghost.end(); it++)
			{
				element_set & ref = ghost_elements[it->first][ElementNum(mask)];
				if( !it->second.empty() ) 
				{
					if( HaveGlobalID(mask) )
						std::sort(it->second.begin(),it->second.end(),GlobalIDComparator(this));
					else
						std::sort(it->second.begin(),it->second.end(),CentroidComparator(this));
				}
				element_set result(ref.size());
				element_set::iterator end;
				if( HaveGlobalID(CELL) )
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),GlobalIDComparator(this));
				else
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),CentroidComparator(this));
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_ghost.clear();
			for(proc_elements::iterator it = del_shared.begin(); it != del_shared.end(); it++)
			{
				element_set & ref = shared_elements[it->first][ElementNum(mask)];
				if( !it->second.empty() ) 
				{
					if( HaveGlobalID(mask) )
						std::sort(it->second.begin(),it->second.end(),GlobalIDComparator(this));
					else
						std::sort(it->second.begin(),it->second.end(),CentroidComparator(this));
				}
				element_set result(ref.size());
				element_set::iterator end;
				if( HaveGlobalID(CELL) )
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),GlobalIDComparator(this));
				else
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),CentroidComparator(this));
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_shared.clear();
#endif //USE_PARALLEL_STORAGE
			for(element_set::iterator it = delete_elements.begin(); it != delete_elements.end(); it++) Destroy(*it);
			delete_elements.clear();
		}		
		DeleteTag(tag_delete);
		ComputeSharedProcs();
		Integer(GetHandle(),tag_layers) = 0;
		Integer(GetHandle(),tag_bridge) = NONE;
#endif //USE_MPI
		EXIT_FUNC();
	}
	
	void Mesh::RemoveGhostElements(const HandleType * ghost, enumerator num)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		int mpirank = GetProcessorRank();
		Tag tag_delete = CreateTag("TEMPORARY_DELETE_GHOST_ELEMENTS_TAG",DATA_INTEGER,CELL | FACE | EDGE | NODE,CELL | FACE | EDGE | NODE);
		element_set delete_elements;
#if defined(USE_PARALLEL_STORAGE)
		proc_elements del_shared, del_ghost;
#endif //USE_PARALLEL_STORAGE
		
		for(ElementType mask = CELL; mask >= NODE; mask = PrevElementType(mask))
		{
			if( mask & CELL )
			{
				for(const HandleType * it = ghost; it != ghost + num; it++)
					if( GetStatus(*it) == Element::Ghost ) Integer(*it,tag_delete) = mpirank;
			}
			else 
			{
				for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
				{
					Element::Status estat = GetStatus(*it);
					if( estat == Element::Owned || estat == Element::Shared ) continue;
					if( it->nbAdjElements(NextElementType(mask)) == 0 ) it->Integer(tag_delete) = mpirank;
				}
			}
			ReduceData(tag_delete,mask,0,DeleteUnpack);
			ExchangeData(tag_delete,mask,0);
			for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			{
				Element::Status estat = GetStatus(*it);
				if( estat == Element::Owned ) continue;
				if( it->HaveData(tag_delete) )
				{
						Storage::integer_array del_procs = it->IntegerArray(tag_delete);
						std::sort(del_procs.begin(),del_procs.end());
						
						if( estat == Element::Ghost && std::binary_search(del_procs.begin(),del_procs.end(),mpirank) )
						{
#if defined(USE_PARALLEL_STORAGE)
							del_ghost[it->IntegerDF(tag_owner)].push_back(*it);
#endif //USE_PARALLEL_STORAGE
							delete_elements.push_back(*it);
						}
						else
						{
							Storage::integer_array procs = it->IntegerArrayDV(tag_processors);
							std::vector<Storage::integer> result(procs.size());
#if defined(USE_PARALLEL_STORAGE)
							if( estat == Element::Shared )
							{
								for(Storage::integer_array::iterator vit = del_procs.begin(); vit != del_procs.end(); vit++)
									del_shared[*vit].push_back(*it);
							}
#endif	//USE_PARALLEL_STORAGE
							std::vector<Storage::integer>::iterator end = std::set_difference(procs.begin(),procs.end(),del_procs.begin(),del_procs.end(),result.begin());
							result.resize(end-result.begin());
							procs.clear();
							procs.insert(procs.begin(),result.begin(),result.end());
							
							if( procs.size() == 1 && procs[0] == mpirank )
								SetStatus(*it,Element::Owned);
						}
				}
			}
#if defined(USE_PARALLEL_STORAGE)
			for(proc_elements::iterator it = del_ghost.begin(); it != del_ghost.end(); it++)
			{
				element_set & ref = ghost_elements[it->first][ElementNum(mask)];
				if( !it->second.empty() ) 
				{
					if( HaveGlobalID(mask) )
						std::sort(it->second.begin(),it->second.end(),GlobalIDComparator(this));
					else
						std::sort(it->second.begin(),it->second.end(),CentroidComparator(this));
				}
				element_set result(ref.size());
				element_set::iterator end;
				if( HaveGlobalID(CELL) )
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),GlobalIDComparator(this));
				else
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),CentroidComparator(this));
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_ghost.clear();
			for(proc_elements::iterator it = del_shared.begin(); it != del_shared.end(); it++)
			{
				element_set & ref = shared_elements[it->first][ElementNum(mask)];
				if( !it->second.empty() ) 
				{
					if( HaveGlobalID(mask) )
						std::sort(it->second.begin(),it->second.end(),GlobalIDComparator(this));
					else
						std::sort(it->second.begin(),it->second.end(),CentroidComparator(this));
				}
				element_set result(ref.size());
				element_set::iterator end;
				if( HaveGlobalID(CELL) )
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),GlobalIDComparator(this));
				else
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),CentroidComparator(this));
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_shared.clear();
#endif //USE_PARALLEL_STORAGE
			for(element_set::iterator it = delete_elements.begin(); it != delete_elements.end(); it++) Destroy(*it);
			delete_elements.clear();			
		}
		DeleteTag(tag_delete);
		ComputeSharedProcs();
#else //USE_MPI
		(void) ghost;
		(void) num;
#endif //USE_MPI
		EXIT_FUNC();
	}

	
	void Mesh::AssignGlobalID(ElementType mask)
	{
		tag_global_id = CreateTag("GLOBAL_ID",DATA_INTEGER, mask, NONE,1);
		ENTER_FUNC();

#if defined(USE_MPI)
		if( m_state == Parallel )
		{
			INMOST_DATA_BIG_ENUM_TYPE number,shift[4], shift_recv[4], local_shift;
			//int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber();
			int num;
			//std::vector<INMOST_DATA_BIG_ENUM_TYPE> shifts(mpisize*4);
			memset(shift,0,sizeof(INMOST_DATA_BIG_ENUM_TYPE));
			for(ElementType currenttype = NODE; currenttype <= CELL; currenttype = NextElementType(currenttype) )
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
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
			}
			for(ElementType currenttype = NODE; currenttype <= CELL; currenttype = NextElementType(currenttype) )
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
			if( tag_global_id.isValid() )
			{
				Tag exchange = tag_global_id;
				//Protect ExchangeData from using global id to sort ghost/shared elements
				tag_global_id = Tag();
				ExchangeData(exchange,mask,0);
				tag_global_id = exchange;
			}
			SortParallelStorage(mask);
		}
		else
		{
#endif //USE_MPI
			INMOST_DATA_BIG_ENUM_TYPE number;
			for(ElementType currenttype = NODE; currenttype <= CELL; currenttype = NextElementType(currenttype) )
				if( mask & currenttype )
				{
					number = 0;
					for(Mesh::iteratorElement it = BeginElement(currenttype); it != EndElement(); it++)
						it->Integer(tag_global_id) = number++;
				}
#if defined(USE_MPI)
		}
#endif //USE_MPI
		EXIT_FUNC();
	}
		
	void Mesh::ComputeSharedProcs()
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		std::set<int> shared_procs;
		int mpirank = GetProcessorRank();
#if defined(USE_OMP)
#pragma omp parallel
#endif
		{
			std::set<int> shared_procs_local;
#if defined(USE_OMP)
#pragma omp for
#endif
			for(integer nit = 0; nit < NodeLastLocalID(); nit++) if( isValidNode(nit) )
			{
				Node it = NodeByLocalID(nit);
				integer_array p = it->IntegerArrayDV(tag_processors);
				for(integer_array::iterator kt = p.begin(); kt != p.end(); kt++)
					shared_procs_local.insert(*kt);
			}
#if defined(USE_OMP)
#pragma omp critical
#endif
			{
				shared_procs.insert(shared_procs_local.begin(),shared_procs_local.end());
			}
		}
		std::set<int>::iterator ir = shared_procs.find(mpirank);
		if( ir != shared_procs.end() ) shared_procs.erase(ir);
		Storage::integer_array procs = IntegerArrayDV(GetHandle(),tag_processors);
		procs.clear();
		procs.insert(procs.begin(),shared_procs.begin(),shared_procs.end());
		
		REPORT_VAL("processors",procs.size());
#endif
		EXIT_FUNC();
	}
	
	
	
	
	Mesh::proc_elements Mesh::ComputeSharedSkinSet(ElementType bridge_type)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		if( !HaveGlobalID(CELL) ) AssignGlobalID(CELL);

		int mpirank = GetProcessorRank();
		Tag tag_skin = CreateTag("TEMPORARY_COMPUTE_SHARED_SKIN_SET",DATA_INTEGER,FACE,FACE);

		REPORT_STR("filling neighbouring cell's global identificators for faces")
#if defined(USE_PARALLEL_WRITE_TIME)
		std::map<int,int> numfacesperproc;
#endif

		for(iteratorFace it = BeginFace(); it != EndFace(); it++)
		{
			Element::Status estat = GetStatus(*it);
			if( estat == Element::Ghost || estat == Element::Shared )
			{
				Storage::integer_array arr = it->IntegerArray(tag_skin);
				ElementArray<Element> adj = it->getAdjElements(CELL);
        //REPORT_STR("face " << it->LocalID() << " global " << it->GlobalID() << " type " << Element::StatusName(estat));
				for(ElementArray<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
				{
					arr.push_back(jt->GlobalID()); //identificator of the cell
					arr.push_back(jt->IntegerDF(tag_owner)); //owner of the cell
          //REPORT_STR("id " << jt->GlobalID() << " owner " << jt->IntegerDF(tag_owner));
				}
			}
#if defined(USE_PARALLEL_WRITE_TIME)
        ++numfacesperproc[it->IntegerDF(tag_owner)];
#endif
		}
		REPORT_STR("number of ghosted or shared faces");
#if defined(USE_PARALLEL_WRITE_TIME)
		for(std::map<int,int>::iterator it = numfacesperproc.begin(); it != numfacesperproc.end(); ++it)
		{
			REPORT_VAL("processor",it->first);
			REPORT_VAL("skin faces",it->second);
		}
#endif

		REPORT_STR("reducing cell's global identificators information")

		ReduceData(tag_skin,FACE,0,UnpackSkin);


    REPORT_STR("exchanging cell's global identificators information")

#if defined(USE_PARALLEL_WRITE_TIME)
    //for(iteratorFace it = BeginFace(); it != EndFace(); it++)
		//{
		//	Element::Status estat1 = GetStatus(*it), estat2;
		//	if( estat1 == Element::Owned ) continue;
    //  Storage::integer_array skin_data = it->IntegerArray(tag_skin);
    //  REPORT_STR("face " << it->LocalID() << " global " << it->GlobalID() << " type " << Element::StatusName(estat1));
    //  for(Storage::integer_array::iterator kt = skin_data.begin(); kt != skin_data.end(); kt+=2)
    //  {
    //    REPORT_STR("id " << *kt << " owner " << *(kt+1));
    //  }
    //}
#endif

		ExchangeData(tag_skin,FACE,0);


		REPORT_STR("synchornization done")

		proc_elements skin_faces;
#if defined(USE_PARALLEL_WRITE_TIME)
    numfacesperproc.clear();
#endif
		for(iteratorFace it = BeginFace(); it != EndFace(); it++)
		{
			bool flag = false;
			Element::Status estat1 = GetStatus(*it), estat2;
			if( estat1 == Element::Owned ) continue;
      //Storage::integer_array skin_data = it->IntegerArray(tag_skin);
      //REPORT_STR("face " << it->LocalID() << " global " << it->GlobalID() << " type " << Element::StatusName(estat1));
      //for(Storage::integer_array::iterator kt = skin_data.begin(); kt != skin_data.end(); kt+=2)
      //{
      //  REPORT_STR("id " << *kt << " owner " << *(kt+1));
      //}
			Cell c1 = it->BackCell();
			Cell c2 = it->FrontCell();
			if( !c1.isValid() && !c2.isValid() ) continue; //no cells per face - skip hanging face
			Storage::integer_array skin_data = it->IntegerArray(tag_skin);
			if( skin_data.size()/2 < 2 ) continue; //only one cell per face - skip boundary face
			assert( skin_data.size()/2 == 2 ); //more then 2 cells per face - invalid grid topology
			if( !c2.isValid() ) //other cell is not present on current processor
			{
				estat1 = c1->GetStatus();
				if( estat1 == Element::Owned || estat1 == Element::Shared )
					flag = true;
			}
			else //other cell is present on current processor
			{
				estat1 = c1->GetStatus();
				estat2 = c2->GetStatus();
				if( (estat1 == Element::Ghost && estat2 == Element::Shared) ||
					(estat2 == Element::Ghost && estat1 == Element::Shared) ||
					(estat1 == Element::Owned && estat2 == Element::Ghost)  ||
					(estat2 == Element::Owned && estat1 == Element::Ghost))
					flag = true;
			}
			
			if( flag )
			{
        //printf("hello!\n");
				for(int i = 0; i < 2; i++) //assert checks that there are two cells
				{
					Storage::integer owner = skin_data[i*2+1]; //cell owner
					if( owner != mpirank )
					{
						skin_faces[owner].push_back(*it);
#if defined(USE_PARALLEL_WRITE_TIME)
						++numfacesperproc[owner];
#endif
						break;
					}
				}	
			}
		}
		REPORT_STR("number of shared faces");
#if defined(USE_PARALLEL_WRITE_TIME)
		for(std::map<int,int>::iterator it = numfacesperproc.begin(); it != numfacesperproc.end(); ++it)
		{
			REPORT_VAL("processor",it->first);
			REPORT_VAL("skin faces",it->second);
		}
#endif
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
			for(proc_elements::iterator it = skin_faces.begin(); it != skin_faces.end(); it++)
			for(element_set::iterator f = it->second.begin(); f != it->second.end(); f++)
			{
				ElementArray<Node> fnodes = Element(this,*f)->getNodes();
				keynum = GetProcessorRank()*GetProcessorsNumber()+(it->first+1); 
				fwrite(&keynum,sizeof(Storage::integer),1,file);
				keynum = static_cast<Storage::integer>(fnodes.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
				for(ElementArray<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
					fwrite(&fn->Coords()[0],sizeof(Storage::real),1,file);
				for(ElementArray<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
					fwrite(&fn->Coords()[1],sizeof(Storage::real),1,file);
				for(ElementArray<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
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
#endif //DEBUG_COMPUTE_SHARED_SKIN_SET
		proc_elements bridge;
		if( bridge_type & FACE ) bridge.swap(skin_faces);
		else
		{
			element_set all_visited;
			Tag on_skin = CreateTag("TEMPORARY_ON_SKIN",DATA_INTEGER,bridge_type,bridge_type);

			REPORT_STR("gathering specified elements on skin");
			REPORT_VAL("type",ElementTypeName(bridge_type));

			MarkerType busy = CreateMarker();
			for(proc_elements::iterator kt = skin_faces.begin(); kt != skin_faces.end(); kt++)
			{
				for(element_set::iterator it = kt->second.begin(); it != kt->second.end(); it++)
				{
					ElementArray<Element> adj = Element(this,*it)->getAdjElements(bridge_type);
					for(ElementArray<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
					{
						if( !jt->GetMarker(busy) )
						{
							jt->IntegerArray(on_skin).push_back(mpirank); //put current rank for element
							jt->SetMarker(busy);
							all_visited.push_back(*jt);
						}
					}
				}
			}
			if( !all_visited.empty() ) RemMarkerArray(&all_visited[0],static_cast<enumerator>(all_visited.size()),busy);
			//for(element_set::iterator it = all_visited.begin(); it != all_visited.end(); ++it) RemMarker(*it,busy);
			ReleaseMarker(busy);
			
			REPORT_STR("exchange of data");
			ReduceData(on_skin,bridge_type,0,UnpackOnSkin); //on each processor each element will know
			ExchangeData(on_skin,bridge_type,0);            //what processors share this element
			
#if defined(USE_PARALLEL_WRITE_TIME)
			std::map<int,int> numelemsperproc;
#endif
			for(element_set::iterator it = all_visited.begin(); it != all_visited.end(); it++)
			{
				Storage::integer_array os = IntegerArray(*it,on_skin);
        for(Storage::integer_array::iterator p = os.begin(); p != os.end(); p++)
				{
					if( *p != mpirank )
					{
						bridge[*p].push_back(*it);
#if defined(USE_PARALLEL_WRITE_TIME)
						++numelemsperproc[*p];
#endif
					}
				}
			}
			REPORT_STR("number of shared elements")
#if defined(USE_PARALLEL_WRITE_TIME)
			for(std::map<int,int>::iterator it = numelemsperproc.begin(); it != numelemsperproc.end(); ++it)
			{
				REPORT_VAL("processor",it->first);
				REPORT_VAL("skin elements",it->second);
			}
#endif
		}
		EXIT_FUNC();		
		return bridge;		
#else //USE_MPI
		(void) bridge_type;
		EXIT_FUNC();
		return proc_elements();
#endif //USE_MPI
	}
	
	
	void Mesh::PackTagData(const Tag & tag, const elements_by_type & elements, int destination, ElementType mask, MarkerType select, buffer_type & buffer)
	{
		if( tag.GetDataType() == DATA_REMOTE_REFERENCE ) return; //NOT IMPLEMENTED TODO 14
		ENTER_FUNC();
#if defined(USE_MPI)
		REPORT_VAL("Processor",destination);
		REPORT_VAL("Buffer size before pack",buffer.size());
		ElementType pack_types[2] = {NONE,NONE};
		element_set::const_iterator eit;
		buffer_type array_data_send;
		std::vector<INMOST_DATA_ENUM_TYPE> array_size_send;
		array_data_send.reserve(4096);
		array_size_send.reserve(4096);
		unsigned int size = tag.GetSize();
		for(int i = ElementNum(NODE); i <= ElementNum(ESET); i++) if( (mask & ElementTypeFromDim(i)) && tag.isDefinedByDim(i) )
		{
			pack_types[0] |= ElementTypeFromDim(i);
			REPORT_VAL("select marker",select);
			REPORT_VAL("pack for type",ElementTypeName(ElementTypeFromDim(i)));
			REPORT_VAL("number of elements",elements[i].size());
			int total_packed = 0;
			if( tag.isSparseByDim(i) )
			{
				pack_types[1] |= ElementTypeFromDim(i);
				unsigned count = static_cast<unsigned>(array_size_send.size());
				array_size_send.push_back(0);
				for(eit = elements[i].begin(); eit != elements[i].end(); eit++)
					if( (!select || GetMarker(*eit,select)) && HaveData(*eit,tag) )
					{
            //REPORT_STR("element type " << ElementTypeName(GetHandleElementType(*eit)) << " global id " << Integer(*eit,GlobalIDTag()));
						REPORT_VAL("element index", static_cast<INMOST_DATA_ENUM_TYPE>(eit-elements[i].begin()));

						array_size_send.push_back(static_cast<INMOST_DATA_ENUM_TYPE>(eit-elements[i].begin()));
						array_size_send[count]++;
						INMOST_DATA_ENUM_TYPE s = GetDataSize(*eit,tag);
						INMOST_DATA_ENUM_TYPE had_s = static_cast<INMOST_DATA_ENUM_TYPE>(array_data_send.size());
						//array_data_send.resize(had_s+s*tag.GetBytesSize());
						if( s )
						{
#if defined(USE_AUTODIFF)
							if( tag.GetDataType() == DATA_VARIABLE )
							{
								REPORT_VAL("data size: ", s);
								REPORT_VAL("data capacity: ", GetDataCapacity(*eit,tag));
								REPORT_VAL("array size: ", had_s);
							}
#endif
							array_data_send.resize(had_s+GetDataCapacity(*eit,tag));
							GetData(*eit,tag,0,s,&array_data_send[had_s]);
            //REPORT_VAL("size",s);
            //for(int qq = 0; qq < s; ++qq)
            //{
            //  REPORT_VAL("value " << qq, (*(Storage::integer *)&array_data_send[had_s+qq*tag.GetBytesSize()]));
            //}
						}
						if( size == ENUMUNDEF ) array_size_send.push_back(s);
						++total_packed;
					}
				REPORT_VAL("count",array_size_send[count]);
				REPORT_VAL("index", count);
			}
			else
			{
				for(eit = elements[i].begin(); eit != elements[i].end(); eit++) if( !select || GetMarker(*eit,select) )
				{
         // REPORT_STR("element type " << ElementTypeName(*eit) << " global id " << Integer(*eit,GlobalIDTag()));
					INMOST_DATA_ENUM_TYPE s = GetDataSize(*eit,tag);
					INMOST_DATA_ENUM_TYPE had_s = static_cast<INMOST_DATA_ENUM_TYPE>(array_data_send.size());
					//array_data_send.resize(had_s+s*tag.GetBytesSize());
					if( s )
					{
#if defined(USE_AUTODIFF)
						if( tag.GetDataType() == DATA_VARIABLE )
						{
							REPORT_VAL("on element ",Element(this,*eit).GlobalID());
							REPORT_VAL("position: ", had_s);
							REPORT_VAL("data capacity: ", GetDataCapacity(*eit,tag));
							REPORT_VAL("size: ", s);
						}
#endif
						array_data_send.resize(had_s+GetDataCapacity(*eit,tag));
                        if (tag.GetDataType() == DATA_REFERENCE)
                        {
                            reference_array refs = ReferenceArray(*eit, tag);
                            int bytes = tag.GetBytesSize();
                            for(Storage::reference_array::size_type i = 0; i < refs.size(); ++i)
                            {
                                if (refs[i] == InvalidElement()) continue;
                                HandleType data = ComposeHandle(refs[i]->GetElementType(), refs[i]->GlobalID());
                                memcpy(&array_data_send[had_s+i*bytes],&data,sizeof(HandleType));
                            }
                        }
                        else
                        {
                            GetData(*eit,tag,0,s,&array_data_send[had_s]);
                        }
					}
					if( size == ENUMUNDEF ) array_size_send.push_back(s);
					++total_packed;
				}
			}
			REPORT_VAL("total packed records",total_packed);
		}
		//static INMOST_DATA_ENUM_TYPE tail = 1;
		//REPORT_VAL("add tail",tail);
		//array_size_send.push_back(tail++);
		INMOST_DATA_ENUM_TYPE size_send, data_send;
		size_send = static_cast<INMOST_DATA_ENUM_TYPE>(array_size_send.size());
		data_send = static_cast<INMOST_DATA_ENUM_TYPE>(array_data_send.size());
		REPORT_VAL("tag defined on",static_cast<int>(pack_types[0]));
		REPORT_VAL("tag sparse on",static_cast<int>(pack_types[1]));
		REPORT_VAL("size_send",size_send);
		REPORT_VAL("data_send",data_send);
		int buffer_size = 0,position = static_cast<int>(buffer.size()),temp,bytes;
		bytes = tag.GetPackedBytesSize();
		MPI_Pack_size(2                                                                      ,INMOST_MPI_DATA_BULK_TYPE,comm,&temp); buffer_size+= temp;
		MPI_Pack_size(1                                                                      ,INMOST_MPI_DATA_ENUM_TYPE,comm,&temp); buffer_size+= temp;
		MPI_Pack_size(1                                                                      ,INMOST_MPI_DATA_ENUM_TYPE,comm,&temp); buffer_size+= temp;
		MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(array_size_send.size())                   ,INMOST_MPI_DATA_ENUM_TYPE,comm,&temp); buffer_size+= temp;
		MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(array_data_send.size()/bytes),tag.GetBulkDataType()    ,comm,&temp); buffer_size+= temp;
		buffer.resize(position+buffer_size);
		MPI_Pack(pack_types,2,INMOST_MPI_DATA_BULK_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
		MPI_Pack(&size_send,1,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
		MPI_Pack(&data_send,1,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
		if( !array_size_send.empty() ) MPI_Pack(&array_size_send[0],static_cast<INMOST_MPI_SIZE>(array_size_send.size()),INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
		if( !array_data_send.empty() ) MPI_Pack(&array_data_send[0],static_cast<INMOST_MPI_SIZE>(array_data_send.size()/bytes),tag.GetBulkDataType(),&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
		buffer.resize(position);
		REPORT_VAL("Buffer size after pack",buffer.size());
		if( size_send < 6 )
		{
			for(INMOST_DATA_ENUM_TYPE qq = 0; qq < size_send; ++qq)
			{
				REPORT_VAL("array_size_send["<<qq<<"]",array_size_send[qq]);
			}
		}
#else
		(void) tag;
		(void) elements;
        (void) destination;
		(void) mask;
		(void) select;
		(void) buffer;
#endif
		REPORT_VAL("TagName",tag.GetTagName());
		EXIT_FUNC();
	}
	


    bool Mesh::FindSharedGhost(int global_id, INMOST_DATA_INTEGER_TYPE el_type_num,  HandleType& res)
    {
        int dim = el_type_num;
        for (parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
        {
            for(element_set::iterator p = it->second[dim].begin(); p != it->second[dim].end(); p++) 
            {
                if (GlobalID(*p) == global_id)
                {
                    res = *p;
                    return true;
                }
            }
        }
        for (parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)
        {
            for(element_set::iterator p = it->second[dim].begin(); p != it->second[dim].end(); p++) 
            {
                if (GlobalID(*p) == global_id)
                {
                    res = *p;
                    return true;
                }
            }
        }
        return false;
     }               


	void Mesh::UnpackTagData(const Tag & tag, const elements_by_type & elements, ElementType mask, MarkerType select, buffer_type & buffer, int & position, ReduceOperation op)
	{
        (void) mask;
		if( tag.GetDataType() == DATA_REMOTE_REFERENCE) return; //NOT IMPLEMENTED TODO 14
		ENTER_FUNC();
		REPORT_VAL("TagName",tag.GetTagName());
		REPORT_VAL("select marker",select);
#if defined(USE_MPI)
		REPORT_VAL("Position before unpack",position);
		if( !buffer.empty() )
		{
			int pos = 0, k = 0;
			ElementType recv_mask[2] = {NONE,NONE};
			INMOST_DATA_ENUM_TYPE data_recv, size_recv;
			element_set::const_iterator eit;
			unsigned int size = tag.GetSize();
			std::vector<INMOST_DATA_BULK_TYPE> array_data_recv;
			std::vector<INMOST_DATA_ENUM_TYPE> array_size_recv;
			MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,recv_mask ,2,INMOST_MPI_DATA_BULK_TYPE,comm);
			//assert(recv_mask[0] == mask);//Element types mask is not synchronized among processors, this may lead to nasty errors
			MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&size_recv,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&data_recv,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			REPORT_VAL("size of size array",size_recv);
			REPORT_VAL("size of data array",data_recv);
			array_size_recv.resize(size_recv);
			array_data_recv.resize(data_recv);
			if( !array_size_recv.empty() ) 
			{
				MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&array_size_recv[0],static_cast<INMOST_MPI_SIZE>(array_size_recv.size()),INMOST_MPI_DATA_ENUM_TYPE,comm);
				REPORT_VAL("size array last", array_size_recv.back());
			}

			if( !array_data_recv.empty() )
			{
				int bytes = tag.GetPackedBytesSize();
				REPORT_VAL("occupied by type",bytes);
				REPORT_VAL("bytes in entry",sizeof(Sparse::Row::entry));
				REPORT_VAL("stored in type",tag.GetBytesSize());
				REPORT_VAL("all size recv",array_data_recv.size());
				REPORT_VAL("incoming size of data",data_recv);
				REPORT_VAL("calculated size of data",array_data_recv.size()/bytes);
				REPORT_VAL("calculated size of data",array_data_recv.size()/sizeof(Sparse::Row::entry));
				MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&array_data_recv[0],static_cast<INMOST_MPI_SIZE>(array_data_recv.size()/bytes),tag.GetBulkDataType(),comm);
			}
			REPORT_VAL("Position after unpack",position);
			for(int i = ElementNum(NODE); i <= ElementNum(ESET); i++) if( (recv_mask[0] & ElementTypeFromDim(i)) )
			{
				REPORT_VAL("unpack for type",ElementTypeName(ElementTypeFromDim(i)));
				REPORT_VAL("number of elements",elements[i].size());
				if( !tag.isDefinedByDim(i) ) 
				{
					CreateTag(tag.GetTagName(),tag.GetDataType(),ElementTypeFromDim(i),recv_mask[1] & ElementTypeFromDim(i),size);
				}
				int total_unpacked = 0;
				int total_skipped = 0;
				if( tag.isSparseByDim(i) )
				{
					REPORT_VAL("sparse for type",ElementTypeName(ElementTypeFromDim(i)));
					REPORT_VAL("index", k);
					REPORT_VAL("count", array_size_recv[k]);
					REPORT_VAL("pos",pos);
					unsigned count = static_cast<unsigned>(array_size_recv[k++]);
					if( size == ENUMUNDEF )
					{
						REPORT_STR("variable size");
						for(unsigned j = 0; j < count; j++)
						{
							eit = elements[i].begin() + array_size_recv[k++];
							//REPORT_STR("element type " << ElementTypeName(GetHandleElementType(*eit)) << " global id " << Integer(*eit,GlobalIDTag()));
							assert( !select || GetMarker(*eit,select) ); //if fires then very likely that marker was not synchronized
							//REPORT_VAL("size",array_size_recv[k]);
							//for(int qq = 0; qq < array_size_recv[k]; ++qq)
							//{
							//  REPORT_VAL("value " << qq, (*(Storage::integer *)&array_data_recv[pos+qq*tag.GetBytesSize()]));
							//}
							op(tag,Element(this,*eit),&array_data_recv[pos],array_size_recv[k]);
							pos += GetDataCapacity(&array_data_recv[pos],array_size_recv[k],tag);
							//pos += array_size_recv[k]*tag.GetBytesSize();
							++k;
							++total_unpacked;
						}
					}
					else
					{
						REPORT_STR("fixed size");
						for(unsigned j = 0; j < count; j++)
						{
							REPORT_VAL("element index",array_size_recv[k]);
							REPORT_VAL("pos",pos);
							eit = elements[i].begin() + array_size_recv[k++];
							assert( !select || GetMarker(*eit,select) ); //if fires then very likely that marker was not synchronized
							if( pos >= array_data_recv.size() )
							{
								REPORT_STR("ERROR!");
								REPORT_VAL("pos", pos );
								REPORT_VAL("array_size", array_size_recv.size());
								REPORT_VAL("array_data", array_data_recv.size());
								REPORT_VAL("size_recv", size_recv);
								REPORT_VAL("data_recv", data_recv);
								REPORT_VAL("count", count);
								REPORT_VAL("buffer position",position);
								REPORT_VAL("size[0]", array_size_recv[0]);
								REPORT_VAL("size[1]", array_size_recv[1]);
								REPORT_VAL("size[2]", array_size_recv[2]);
								REPORT_VAL("k",k);
							}
							op(tag,Element(this,*eit),&array_data_recv[pos],size);
							pos += GetDataCapacity(&array_data_recv[pos],size,tag);
							//pos += size*tag.GetBytesSize();
							++total_unpacked;
						}
					}
				}
				else
				{
					REPORT_VAL("dense for type",ElementTypeName(ElementTypeFromDim(i)));
					if( size == ENUMUNDEF )
					{
						REPORT_STR("variable size");
						for(eit = elements[i].begin(); eit != elements[i].end(); eit++) 
						{
							if( !select || GetMarker(*eit,select) )
							{
                                if (tag.GetDataType() == DATA_REFERENCE)
                                {
                                    int bytes = tag.GetBytesSize();
                                    for (int i = 0; i < array_size_recv[k]; i++)
                                    {
                                        int global_id;
                                        HandleType data;
                                        memcpy(&data,&array_data_recv[pos + i*bytes],sizeof(HandleType));
                                        global_id = GetHandleID(data);
                                        int type = GetHandleElementNum(data);
                                        HandleType target;
                                        if (FindSharedGhost(global_id,GetHandleElementNum(data),target))
                                        {
                                            TagReferenceArray ref_tag = tag;
                                            ref_tag[*eit].push_back(target);
                                        }
                                    }
                                }
                                else
                                {
								    op(tag,Element(this,*eit),&array_data_recv[pos],array_size_recv[k]);
                                }

								pos += GetDataCapacity(&array_data_recv[pos],array_size_recv[k],tag);
								//pos += array_size_recv[k]*tag.GetBytesSize();
								++k;
								++total_unpacked;
							} 
							else ++total_skipped;
						}
					}
					else
					{
						REPORT_STR("fixed size");
						for(eit = elements[i].begin(); eit != elements[i].end(); eit++) 
						{
							if( !select || GetMarker(*eit,select) )
							{
#if defined(USE_AUTODIFF)
								if( tag.GetDataType() == DATA_VARIABLE )
								{
									REPORT_VAL("on element ",Element(this,*eit).GlobalID());
									REPORT_VAL("position ", pos);
									REPORT_VAL("capacity ", GetDataCapacity(&array_data_recv[pos],size,tag));
									REPORT_VAL("size ", size);
								}
#endif
								op(tag,Element(this,*eit),&array_data_recv[pos],size);
								pos += GetDataCapacity(&array_data_recv[pos],size,tag);
								//pos += size*tag.GetBytesSize();
								++total_unpacked;
							}
							else total_skipped++;
						}
					}
				}
				REPORT_VAL("total skipped",total_skipped);
				REPORT_VAL("total unpacked records",total_unpacked);
			}
		}
		REPORT_VAL("Position after unpack",position);
#else
		(void) tag;
		(void) elements;
		(void) mask;
		(void) select;
		(void) buffer;
		(void) position;
		(void) op;
#endif //USE_MPI
		EXIT_FUNC();
	}
	
	
	void Mesh::ExchangeDataInnerBegin(const tag_set & tags, const parallel_storage & from, const parallel_storage & to, ElementType mask, MarkerType select, exchange_data & storage)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
		//int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber();
		Storage::integer_array procs = IntegerArrayDV(GetHandle(),tag_processors);
		Storage::integer_array::iterator p = procs.begin();
		REPORT_VAL("processors",procs.size());
		storage.send_buffers.resize(procs.size());
		storage.recv_buffers.resize(procs.size());
		REPORT_VAL("send buffers size",storage.send_buffers.size());
		REPORT_VAL("recv buffers size",storage.recv_buffers.size());
		std::vector<int> done;
		parallel_storage::const_iterator find;
		std::vector<INMOST_DATA_ENUM_TYPE> send_size(procs.size(),0), recv_size(procs.size(),0);
		
		bool unknown_size = false;
		for(unsigned int k = 0; k < tags.size(); k++)
		{
			if( tags[k].GetSize() == ENUMUNDEF
#if defined(USE_AUTODIFF)
			   || tags[k].GetDataType() == DATA_VARIABLE
#endif
			   ) unknown_size = true;
			//for(int i = 0; i < 5; ++i)
			//	if( (mask & ElementTypeFromDim(i)) && tags[k].isSparseByDim(i) )
			//		unknown_size = true;
		}
    	int rank = GetProcessorRank();
		//precompute sizes
		for(p = procs.begin(); p != procs.end(); p++ )
		{
			int pos = static_cast<int>(p-procs.begin());
			if( select )
			{
				find = from.find(*p);
				if( find != from.end() )
					for(int i = 0; i < 5; i++)  if( mask & ElementTypeFromDim(i) )
						for(element_set::const_iterator it = find->second[i].begin(); it != find->second[i].end(); ++it)
							if( GetMarker(*it,select) ) send_size[pos]++;
				
				find = to.find(*p);
				if( find != to.end() )
					for(int i = 0; i < 5; i++)  if( mask & ElementTypeFromDim(i) )
						for(element_set::const_iterator it = find->second[i].begin(); it != find->second[i].end(); ++it)
							if( GetMarker(*it,select) ) recv_size[pos]++;
			}
			else
			{
				find = from.find(*p);
				if( find != from.end() )
                {
					for(int i = 0; i < 5; i++)  if( mask & ElementTypeFromDim(i) )
						send_size[pos] += static_cast<INMOST_DATA_ENUM_TYPE>(find->second[i].size());
                }
				find = to.find(*p);
				if( find != to.end() )
					for(int i = 0; i < 5; i++)  if( mask & ElementTypeFromDim(i) )
						recv_size[pos] += static_cast<INMOST_DATA_ENUM_TYPE>(find->second[i].size());
			}
		}
		
		int num_send = 0, num_recv = 0;
        ///////////
        if ( block_recursion == 0)
        {
			bool call_exchange = false;
            for(unsigned int k = 0; k < tags.size(); k++)
            {
                if(tags[k].GetDataType() == DATA_REFERENCE)
                {
                    for(p = procs.begin(); p != procs.end(); p++ )
                    {
                        const elements_by_type& elements = from.find(*p)->second;
                        for(int i = ElementNum(NODE); i <= ElementNum(ESET); i++) if( (mask & ElementTypeFromDim(i)) && tags[k].isDefinedByDim(i) )
                        {
                            for (int j = 0; j < elements[i].size(); j++)
                            {
                                if (!isValidHandleRange(elements[i][j])) continue;
                                reference_array refs = ReferenceArray(elements[i][j], tags[k]);
                                if (refs.size() == 0) continue;
                                if (tags[k] == HighConnTag()) 
                                {
                                    assert(i == ElementNum(ESET));
                                    ElementSet set(this,elements[i][j]);
                                    for(ElementSet child = set.GetChild(); child.isValid(); child = child.GetSibling())
                                    {
                                        child.IntegerArray(tag_sendto).push_back(*p);
                                        call_exchange = true;
                                    }
                                }
                                else
                                {
                                    for(Storage::reference_array::size_type i = 0; i < refs.size(); ++i)
                                    {
                                        if (refs[i] == InvalidElement()) continue;
                                        refs[i].IntegerArray(tag_sendto).push_back(*p);
                                        call_exchange = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if( call_exchange )
            {
		REPORT_STR("call to ExchangeMarked");
				block_recursion = 1;
				ExchangeMarked();
				block_recursion = 0;
			}
        }
        ///////////
		for(p = procs.begin(); p != procs.end(); p++ )
		{
			REPORT_VAL("for processor",*p);//p-procs.begin());
			REPORT_VAL("send size",send_size[p-procs.begin()]);
			REPORT_VAL("recv size",recv_size[p-procs.begin()]);
			if( send_size[p-procs.begin()] )
			{
				for(unsigned int k = 0; k < tags.size(); k++)
					PackTagData(tags[k],from.find(*p)->second,*p,mask,select,storage.send_buffers[num_send].second);
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
						MPI_Pack_size(2,INMOST_MPI_DATA_BULK_TYPE,comm,&temp); buffer_size += temp;
						MPI_Pack_size(1,INMOST_MPI_DATA_ENUM_TYPE,comm,&temp); buffer_size += temp;
						MPI_Pack_size(1,INMOST_MPI_DATA_ENUM_TYPE,comm,&temp); buffer_size += temp;
						MPI_Pack_size(6+n*2,INMOST_MPI_DATA_ENUM_TYPE,comm,&temp); buffer_size += temp;
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
		if( unknown_size && parallel_strategy != 0 ) 
		{
			REPORT_STR("prepare receive with unknown size");
			PrepareReceiveInner(UnknownSize,storage.send_buffers,storage.recv_buffers);
		}
		ExchangeBuffersInner(storage.send_buffers,storage.recv_buffers,storage.send_reqs,storage.recv_reqs);
#else
		(void) tags;
		(void) from;
		(void) to;
		(void) mask;
		(void) select;
		(void) storage;
#endif //USE_MPI
		EXIT_FUNC();
	}
	
	
	void Mesh::ExchangeDataInnerEnd(const tag_set & tags, const parallel_storage & from, const parallel_storage & to, ElementType mask, MarkerType select, ReduceOperation op, exchange_data & storage)
	{
		(void) from;
		{ // checks for bad input
			if( mask == NONE ) return;
			for(tag_set::const_iterator it = tags.begin(); it != tags.end(); ++it) assert( it->isValid() );
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
					UnpackTagData(tags[k],to.find(storage.recv_buffers[*qt].first)->second,mask,select,storage.recv_buffers[*qt].second,position,op);
				}
			}
		}
		if( !storage.send_reqs.empty() )
		{
			REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(storage.send_reqs.size()),&storage.send_reqs[0],MPI_STATUSES_IGNORE));
		}
#else //USE_MPI
		(void) tags;
		(void) from;
		(void) to;
		(void) mask;
		(void) select;
		(void) op;
		(void) storage;
#endif //USE_MPI
		EXIT_FUNC();
	}
	
	

	void Mesh::GatherParallelStorage(parallel_storage & ghost, parallel_storage & shared, ElementType mask)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
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
						shared[*vit][ElementNum(it->GetElementType())].push_back(*it);
			}
			else if( estat == Element::Ghost )
				ghost[it->IntegerDF(tag_owner)][ElementNum(it->GetElementType())].push_back(*it);
		}
		time = Timer() - time;
		REPORT_VAL("time",time);
		SortParallelStorage(ghost,shared,mask);
#else //USE_MPI
		(void) ghost;
		(void) shared;
		(void) mask;
#endif
		EXIT_FUNC();
	}
	
	void Mesh::ExchangeData(const Tag & tag, ElementType mask, MarkerType select)
	{
		ExchangeData(tag_set(1,tag),mask,select);
	}
	
	void Mesh::ExchangeData(const tag_set & tags, ElementType mask, MarkerType select)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif //USE_PARALLEL_STORAGE
		exchange_data storage;
		ExchangeDataInnerBegin(tags,shared_elements,ghost_elements,mask,select,storage);
		ExchangeDataInnerEnd(tags,shared_elements,ghost_elements,mask,select,DefaultUnpack,storage);
#else //USE_MPI
		(void) tags;
		(void) mask;
		(void) select;
#endif //USE_MPI
		EXIT_FUNC();
	}


	void Mesh::ExchangeDataBegin(const Tag & tag, ElementType mask, MarkerType select, exchange_data & storage)
	{
		ExchangeDataBegin(tag_set(1,tag),mask,select,storage);
	}

	void Mesh::ExchangeDataEnd(const Tag & tag, ElementType mask, MarkerType select, exchange_data & storage)
	{
		ExchangeDataEnd(tag_set(1,tag),mask,select,storage);
	}	
	
	void Mesh::ExchangeDataBegin(const tag_set & tags, ElementType mask, MarkerType select, exchange_data & storage)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif //USE_PARALLEL_STORAGE
		ExchangeDataInnerBegin(tags,shared_elements,ghost_elements,mask,select,storage);
#else //USE_MPI
		(void) tags;
		(void) mask;
		(void) select;
		(void) storage;
#endif //USE_MPI
		EXIT_FUNC();
	}
	
	void Mesh::ExchangeDataEnd(const tag_set & tags, ElementType mask, MarkerType select, exchange_data & storage)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif //USE_PARALLEL_STORAGE
		ExchangeDataInnerEnd(tags,shared_elements,ghost_elements,mask,select,DefaultUnpack,storage);
#else //USE_MPI
		(void) tags;
		(void) mask;
		(void) select;
		(void) storage;
#endif //USE_MPI
		EXIT_FUNC();
	}
	
	void Mesh::ReduceData(const Tag & tag, ElementType mask, MarkerType select, ReduceOperation op)
	{
		ReduceData(tag_set(1,tag),mask,select,op);
	}
	
	void Mesh::ReduceData(const tag_set & tags, ElementType mask, MarkerType select, ReduceOperation op)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif //USE_PARALLEL_STORAGE
		exchange_data storage;
		ExchangeDataInnerBegin(tags,ghost_elements,shared_elements,mask,select,storage);
		ExchangeDataInnerEnd(tags,ghost_elements,shared_elements,mask,select,op,storage);
#else //USE_MPI
		(void) tags;
		(void) mask;
		(void) select;
		(void) op;
#endif //USE_MPI
		EXIT_FUNC();
	}
	
	void Mesh::ReduceDataBegin(const Tag & tag, ElementType mask, MarkerType select,exchange_data & storage)
	{
		ReduceDataBegin(tag_set(1,tag),mask,select,storage);
	}
	
	void Mesh::ReduceDataBegin(const tag_set & tags, ElementType mask, MarkerType select, exchange_data & storage)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif //USE_PARALLEL_STORAGE
		ExchangeDataInnerBegin(tags,ghost_elements,shared_elements,mask,select,storage);
#else //USE_MPI
		(void) tags;
		(void) mask;
		(void) select;
		(void) storage;
#endif
		EXIT_FUNC();
	}
	
	void Mesh::ReduceDataEnd(const Tag & tag, ElementType mask, MarkerType select, ReduceOperation op, exchange_data & storage)
	{
		ReduceDataEnd(tag_set(1,tag),mask,select,op,storage);
	}
	
	void Mesh::ReduceDataEnd(const tag_set & tags, ElementType mask, MarkerType select, ReduceOperation op, exchange_data & storage)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#endif //USE_PARALLEL_STORAGE
		ExchangeDataInnerEnd(tags,ghost_elements,shared_elements,mask,select,op,storage);
#else //USE_MPI
		(void) tags;
		(void) mask;
		(void) select;
		(void) op;
		(void) storage;
#endif
		EXIT_FUNC();
	}
	
	void Mesh::PackElementsData(element_set & all, buffer_type & buffer, int destination, const std::vector<std::string> & tag_list)
	{
		ENTER_FUNC();
		REPORT_VAL("dest",destination);
#if defined(USE_MPI)
        int rank = GetProcessorRank();
        //std::cout << ro() << rank << " In pack elements data " << all.size() << std::endl;
        //std::cout << rank << " In pack elements Data" << std::endl;
        INMOST_DATA_ENUM_TYPE num;
		//assume hex mesh for forward allocation 
		//8 nodes per cell, 2 nodes per edge, 4 edges per face, 6 faces per cell
		const INMOST_DATA_ENUM_TYPE size_hint[5] = {8,2,4,6,1};
		INMOST_DATA_ENUM_TYPE prealloc[5] = {0,0,0,0,0};
		int position,new_size,k, mpirank = GetProcessorRank();
		int temp;
		elements_by_type selems;
		//TODO 46 old
		//elements_by_type pack_tags
		{
			for(element_set::iterator it = all.begin(); it != all.end(); it++)
            {
				selems[GetHandleElementNum(*it)].push_back(*it);
            }
			all.clear();
		}
		REPORT_STR("initial number of elements");
		REPORT_VAL("NODE",selems[0].size());
		REPORT_VAL("EDGE",selems[1].size());
		REPORT_VAL("FACE",selems[2].size());
		REPORT_VAL("CELL",selems[3].size());

		Tag arr_position = CreateTag("TEMP_ARRAY_POSITION",DATA_INTEGER,ESET|CELL|FACE|EDGE|NODE,NONE,1);
		//ElementArray<Element> adj(this);
		MarkerType busy = CreateMarker();
        
        // Recursevely looking for all high connections for ESETs
        for(element_set::iterator it = selems[4].begin(); it != selems[4].end(); ++it) SetMarker(*it,busy);
        int ind = 0;
        while (ind < selems[4].size())
        {
            ElementSet set(this,selems[4][ind]);
            // looking to child, sibling and parent
            for (int i = 0; i <= 2; i++)
            {
                ElementSet _set;
                switch (i)
                {
                    case 0: _set = set.GetChild();   break;
                    case 1: _set = set.GetSibling(); break;
                    case 2: _set = set.GetParent();  break;
                }
                if (_set.GetHandle() == InvalidHandle()) continue;
                if (!GetMarker(_set.GetHandle(),busy))
                {
                    selems[4].push_back(_set.GetHandle());
                    SetMarker(_set.GetHandle(),busy);
                }
            }
            ind++;
        }
        for(element_set::iterator it = selems[4].begin(); it != selems[4].end(); ++it) RemMarker(*it,busy);
        // Add low conns elements to selems array. Low conns for ESET can contain any element
        // Low conns for ESETs
        {
            // Mark all elems as busy
		    for(int etypenum = ElementNum(CELL); etypenum > ElementNum(NODE); --etypenum)
			    for(element_set::iterator it = selems[etypenum-1].begin(); it != selems[etypenum-1].end(); ++it) SetMarker(*it,busy);

			for(element_set::iterator it = selems[4].begin(); it != selems[4].end(); it++)
			{
               // GetHandleElementNum(*it) - 4:ESET, 3:CELL, ... , 0:NODE
				Element::adj_type const & adj = LowConn(*it);
				for(Element::adj_type::const_iterator jt = adj.begin(); jt != adj.end(); jt++)
				{
					if( !GetMarker(*jt,busy) )
					{
				        prealloc[GetHandleElementNum(*jt)] += 1;
						selems[GetHandleElementNum(*jt)].push_back(*jt);
						SetMarker(*jt,busy);
					}
				}
			}
		    for(int etypenum = ElementNum(CELL); etypenum >= ElementNum(NODE); --etypenum)
			    for(element_set::iterator it = selems[etypenum].begin(); it != selems[etypenum].end(); ++it) RemMarker(*it,busy);
		}

        // Low conns for CELLs, FACEs, EDGEs 
		for(int etypenum = ElementNum(CELL); etypenum > ElementNum(NODE); --etypenum)
		{
			selems[etypenum-1].reserve(size_hint[etypenum]*selems[etypenum].size());
			for(element_set::iterator it = selems[etypenum-1].begin(); it != selems[etypenum-1].end(); ++it) SetMarker(*it,busy);
			for(element_set::iterator it = selems[etypenum].begin(); it != selems[etypenum].end(); it++)
			{
				Element::adj_type const & adj = LowConn(*it);
				prealloc[etypenum] += static_cast<INMOST_DATA_ENUM_TYPE>(adj.size()); //compute how many connections should be sent
				for(Element::adj_type::const_iterator jt = adj.begin(); jt != adj.end(); jt++)
				{
					if( !GetMarker(*jt,busy) )
					{
						selems[etypenum-1].push_back(*jt);
						SetMarker(*jt,busy);
					}
				}
			}
			for(element_set::iterator it = selems[etypenum-1].begin(); it != selems[etypenum-1].end(); ++it) RemMarker(*it,busy);
		}
		ReleaseMarker(busy);
		//compute how many node connections should be sent
		for(element_set::iterator it = selems[3].begin(); it != selems[3].end(); it++)
			prealloc[0] += static_cast<INMOST_DATA_ENUM_TYPE>(HighConn(*it).size());
			
		for(int etypenum = ElementNum(NODE); etypenum <= ElementNum(ESET); ++etypenum)
		{
			int q = 0;
			for(element_set::iterator it = selems[etypenum].begin(); it != selems[etypenum].end(); it++)
			{
				Integer(*it,arr_position) = q++;
				assert(GetHandleElementNum(*it) == etypenum);
			}
			//TODO: 44 old
			//std::sort(selems[etypenum].begin(),selems[etypenum].end());
		}
		/*
        stringstream ss;
        ss << ro() << mpirank << ": to send: ";
        ss << "nodes: " << selems[0].size() << " | ";
        ss << "edge: "  << selems[1].size() << " | ";
        ss << "faces: " << selems[2].size() << " | ";
        ss << "cells: " << selems[3].size() << " | ";
        ss << "esets: " << selems[4].size() << endl;
        ss << "esets: ";
        for(element_set::iterator it = selems[4].begin(); it != selems[4].end(); it++)
        {
            ss << ElementSet(this,*it).GetName() << " ";
        }
        //cout << ss.str() << endl;
		*/
		REPORT_STR("final number of elements");
		REPORT_VAL("NODE",selems[0].size());
		REPORT_VAL("EDGE",selems[1].size());
		REPORT_VAL("FACE",selems[2].size());
		REPORT_VAL("CELL",selems[3].size());

		REPORT_VAL("prealloc cell nodes",prealloc[0]);
		REPORT_VAL("prealloc edge nodes",prealloc[1]);
		REPORT_VAL("prealloc face edges",prealloc[2]);
		REPORT_VAL("prealloc cell faces",prealloc[3]);
		MarkerType pack_tags_mrk = CreateMarker();
		//pack nodes coords
		{
			integer dim = GetDimensions();
			std::vector<Storage::integer> global_ids;
			position = static_cast<int>(buffer.size());
			new_size = 0;
			MPI_Pack_size(1                                                 ,INMOST_MPI_DATA_ENUM_TYPE   ,comm,&temp); new_size += temp;
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(selems[0].size()*dim),INMOST_MPI_DATA_REAL_TYPE   ,comm,&temp); new_size += temp;
			if( HaveGlobalID(NODE) )
				MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(selems[0].size()),INMOST_MPI_DATA_INTEGER_TYPE,comm,&temp); new_size += temp;
			buffer.resize(position+new_size);
			num = static_cast<INMOST_DATA_ENUM_TYPE>(selems[0].size());
			MPI_Pack(&num,1,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			int marked_for_data = 0, marked_shared = 0;
			{
				k = 0;
				std::vector<Storage::real> coords(selems[0].size()*dim);
				for(element_set::iterator it = selems[0].begin(); it != selems[0].end(); it++)
				{
					Storage::real_array c = RealArray(*it,tag_coords);
					//REPORT_VAL("packed coords",k << " " << c[0] << " " << c[1] << " " << c[2]);
					for(integer j = 0; j < dim; j++) coords[k*GetDimensions()+j] = c[j];
					Storage::integer & owner = IntegerDF(*it,tag_owner);
					{
						if( owner == mpirank )
						{
							Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);
							Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),destination);
							if( ip == proc.end() || (*ip) != destination ) proc.insert(ip,destination);
							if( GetStatus(*it) != Element::Shared) ++marked_shared;
							SetStatus(*it, Element::Shared);
						}
					}
					//TODO: 45
					if( owner != destination ) 
					{
						SetMarker(*it,pack_tags_mrk);
						//TODO 46 old
						//	pack_tags[0].push_back(*it);
						++marked_for_data;
					}
					if( HaveGlobalID(NODE) )
					{
						global_ids.push_back(Integer(*it,GlobalIDTag()));
						//REPORT_VAL("packed global_id",global_ids.back());
					}
					k++;
				}
				if( !coords.empty() ) 
					MPI_Pack(&coords[0]    ,static_cast<INMOST_MPI_SIZE>(selems[0].size()*dim),INMOST_MPI_DATA_REAL_TYPE   ,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
				if( HaveGlobalID(NODE) && !global_ids.empty() )
					MPI_Pack(&global_ids[0],static_cast<INMOST_MPI_SIZE>(selems[0].size())    ,INMOST_MPI_DATA_INTEGER_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			}
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[0].size());
			REPORT_VAL("total marked as shared", marked_shared << " / " << selems[0].size());
			buffer.resize(position);
		}
		//pack edges
		{
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size(selems[1].size());
			std::vector<Storage::integer> low_conn_nums;
			low_conn_nums.reserve(prealloc[1]);
			position = static_cast<int>(buffer.size());
			new_size = 0;
			num = 0; k = 0;
			int marked_for_data = 0, marked_shared = 0;
			for(element_set::iterator it = selems[1].begin(); it != selems[1].end(); it++) 
			{
				low_conn_size[k] = 0;
				Element::adj_type const & lc = LowConn(*it);
				//REPORT_VAL("edge number",k);
				for(Element::adj_type::const_iterator jt = lc.begin(); jt != lc.end(); jt++) if( !Hidden(*jt) )
				{
					//TODO 44 old
					//element_set::iterator find = std::lower_bound(selems[0].begin(),selems[0].end(),(*jt));
					//assert( find != snodes.end() );
					//assert( (*find) == (*jt) ); 
					//low_conn_nums.push_back(static_cast<Storage::integer>(find - selems[0].begin()));
					assert(Integer(*jt,arr_position) == Integer(selems[0][Integer(*jt,arr_position)],arr_position));
					//REPORT_VAL("pack node position",Integer(*jt,arr_position));
					low_conn_nums.push_back(Integer(*jt,arr_position));
					low_conn_size[k]++;
					num++;
				}
				//REPORT_VAL("pack edge nodes ",k << " " << low_conn_nums[low_conn_nums.size()-2] << " " << low_conn_nums[low_conn_nums.size()-1]);
				Storage::integer & owner = IntegerDF(*it,tag_owner);
				{
					if( owner == mpirank )
					{
						Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);
						Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),destination);
						if( ip == proc.end() || (*ip) != destination ) proc.insert(ip,destination);
						if( GetStatus(*it) != Element::Shared) ++marked_shared;
						SetStatus(*it,Element::Shared);
					}
				}
				if( owner != destination ) 
				{
					SetMarker(*it,pack_tags_mrk);
					//TODO 46 old
					//	pack_tags[1].push_back(*it);
					++marked_for_data;
				}
				k++;
			}
			REPORT_VAL("number of edge nodes",num);
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[1].size());
			REPORT_VAL("total marked as shared", marked_shared << " / " << selems[1].size());
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(1+selems[1].size()),INMOST_MPI_DATA_ENUM_TYPE   ,comm,&temp); new_size += temp;
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(num)               ,INMOST_MPI_DATA_INTEGER_TYPE,comm,&temp); new_size += temp;
			//MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(selems[1].size())  ,INMOST_MPI_DATA_INTEGER_TYPE,comm,&temp); new_size += temp;
			buffer.resize(position+new_size);
			temp = static_cast<INMOST_DATA_ENUM_TYPE>(selems[1].size());
			MPI_Pack(&temp,1,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			if( !low_conn_size.empty() ) MPI_Pack(&low_conn_size[0],static_cast<INMOST_MPI_SIZE>(selems[1].size()),INMOST_MPI_DATA_ENUM_TYPE   ,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			if( !low_conn_nums.empty() ) MPI_Pack(&low_conn_nums[0],static_cast<INMOST_MPI_SIZE>(num)             ,INMOST_MPI_DATA_INTEGER_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			buffer.resize(position);
		}
		//pack faces
		{
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size(selems[2].size());
			std::vector<Storage::integer> low_conn_nums;
			low_conn_nums.reserve(prealloc[2]); //assume quads (4 edges)
			position = static_cast<int>(buffer.size());
			new_size = 0;
			num = 0; k = 0;
			int marked_for_data = 0, marked_shared = 0;
			for(element_set::iterator it = selems[2].begin(); it != selems[2].end(); it++) 
			{
				low_conn_size[k] = 0;
				Element::adj_type const & lc = LowConn(*it);
				//REPORT_VAL("face number",k);
				for(Element::adj_type::const_iterator jt = lc.begin(); jt != lc.end(); jt++) if( !Hidden(*jt) )
				{
					//TODO 44 old
					//element_set::iterator find = std::lower_bound(selems[1].begin(),selems[1].end(),(*jt));
					//assert( find != selems[1].end() );
					//assert( (*find) == (*jt) );
					//low_conn_nums.push_back(static_cast<Storage::integer>(find - selems[1].begin()));
					assert(Integer(*jt,arr_position) == Integer(selems[1][Integer(*jt,arr_position)],arr_position));
					//REPORT_VAL("pack edge position",Integer(*jt,arr_position));
					low_conn_nums.push_back(Integer(*jt,arr_position));
					low_conn_size[k]++;
					num++;
				}
				//REPORT_VAL("pack size",low_conn_size[k]);
				Storage::integer & owner = IntegerDF(*it,tag_owner);
				{
					if( owner == mpirank )
					{
						Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);
						Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),destination);
						if( ip == proc.end() || (*ip) != destination ) proc.insert(ip,destination);
						if( GetStatus(*it) != Element::Shared) ++marked_shared;
						SetStatus(*it,Element::Shared);
					}
				}
				if( owner != destination ) 
				{
					SetMarker(*it, pack_tags_mrk);
					//TODO 46 old
					//	pack_tags[2].push_back(*it);
					++marked_for_data;
				}
				k++;
			}
			REPORT_VAL("number of face edges",num);
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[2].size());
			REPORT_VAL("total marked as shared", marked_shared << " / " << selems[2].size());
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(1+selems[2].size()),INMOST_MPI_DATA_ENUM_TYPE   ,comm,&temp); new_size += temp;
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(num)               ,INMOST_MPI_DATA_INTEGER_TYPE,comm,&temp); new_size += temp;
			//MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(selems[2].size())  ,INMOST_MPI_DATA_INTEGER_TYPE,comm,&temp); new_size += temp;
			buffer.resize(position+new_size);
			temp = static_cast<INMOST_DATA_ENUM_TYPE>(selems[2].size());
			MPI_Pack(&temp,1,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			if( !low_conn_size.empty() ) MPI_Pack(&low_conn_size[0],static_cast<INMOST_MPI_SIZE>(selems[2].size()),INMOST_MPI_DATA_ENUM_TYPE   ,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			if( !low_conn_nums.empty() ) MPI_Pack(&low_conn_nums[0],static_cast<INMOST_MPI_SIZE>(num)             ,INMOST_MPI_DATA_INTEGER_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			buffer.resize(position);
		}
		//pack cells
		{
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size(selems[3].size());
			std::vector<INMOST_DATA_ENUM_TYPE> high_conn_size(selems[3].size());
			std::vector<Storage::integer> low_conn_nums;
			std::vector<Storage::integer> high_conn_nums;
			low_conn_nums.reserve(prealloc[3]); //assume hexes (6 faces)
			high_conn_size.reserve(prealloc[0]); //assume hexes (8 nodes)
			INMOST_DATA_ENUM_TYPE num_high = 0;
			position = static_cast<int>(buffer.size());
			new_size = 0;
			num = 0; k = 0;
			int marked_for_data = 0, marked_shared = 0;
			for(element_set::iterator it = selems[3].begin(); it != selems[3].end(); it++) 
			{
				//REPORT_VAL("cell number",k);
				low_conn_size[k] = 0;
				Element::adj_type const & lc = LowConn(*it);
				for(Element::adj_type::const_iterator jt = lc.begin(); jt != lc.end(); jt++) if( !Hidden(*jt) )
				{
					//TODO 44 old
					//element_set::iterator find = std::lower_bound(selems[2].begin(),selems[2].end(),(*jt));
					//assert( find != sfaces.end());
					//assert( (*find) == (*jt) );
					//low_conn_nums.push_back(static_cast<Storage::integer>(find - selems[2].begin()));
					assert(Integer(*jt,arr_position) == Integer(selems[2][Integer(*jt,arr_position)],arr_position));
					//REPORT_VAL("pack face position",Integer(*jt,arr_position));
					low_conn_nums.push_back(Integer(*jt,arr_position));
					low_conn_size[k]++;
					num++;
				}
				//REPORT_VAL("pack low size",low_conn_size[k]);
				high_conn_size[k] = 0;
				Element::adj_type const & hc = HighConn(*it);
				for(Element::adj_type::const_iterator jt = hc.begin(); jt != hc.end(); jt++) if( !Hidden(*jt) )
				{
					//TODO 44 old
					//element_set::iterator find = std::lower_bound(selems[0].begin(),selems[0].end(),(*jt));
					//assert( find == selems[0].end());
					//assert((*find) == (*jt) );
					//high_conn_nums.push_back(static_cast<Storage::integer>(find - selems[0].begin()));
					assert(Integer(*jt,arr_position) == Integer(selems[0][Integer(*jt,arr_position)],arr_position));
					//REPORT_VAL("pack node position",Integer(*jt,arr_position));
					high_conn_nums.push_back(Integer(*jt,arr_position));
					high_conn_size[k]++;
					num_high++;
				}
				//REPORT_VAL("pack high size",high_conn_size[k]);
				Storage::integer & owner = IntegerDF(*it,tag_owner);
				{
					if( owner == mpirank )
					{
						Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);
						Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),destination);
						if( ip == proc.end() || (*ip) != destination ) proc.insert(ip,destination);
						if( GetStatus(*it) != Element::Shared) ++marked_shared;
						SetStatus(*it,Element::Shared);
					}
				}
				if( owner != destination ) 
				{
					SetMarker(*it, pack_tags_mrk);
					//TODO 46 old
					//	pack_tags[3].push_back(*it);
					++marked_for_data;
				}
				k++;
			}
			REPORT_VAL("number of cell faces",num);
			REPORT_VAL("number of cell nodes",num_high);
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[3].size());
			REPORT_VAL("total marked as shared", marked_shared << " / " << selems[3].size());
			MPI_Pack_size(1                                             ,INMOST_MPI_DATA_ENUM_TYPE   ,comm,&temp); new_size += temp;
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(selems[3].size()),INMOST_MPI_DATA_ENUM_TYPE   ,comm,&temp); new_size += temp;
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(num)             ,INMOST_MPI_DATA_INTEGER_TYPE,comm,&temp); new_size += temp;
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(selems[3].size()),INMOST_MPI_DATA_ENUM_TYPE   ,comm,&temp); new_size += temp;
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(num_high)        ,INMOST_MPI_DATA_INTEGER_TYPE,comm,&temp); new_size += temp;
			//MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(selems[3].size()),INMOST_MPI_DATA_INTEGER_TYPE,comm,&temp); new_size += temp;
			buffer.resize(position+new_size);
			temp = static_cast<INMOST_DATA_ENUM_TYPE>(selems[3].size());
			MPI_Pack(&temp,1,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			if( !low_conn_size.empty() )  MPI_Pack(&low_conn_size[0] ,static_cast<INMOST_MPI_SIZE>(selems[3].size()),INMOST_MPI_DATA_ENUM_TYPE   ,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			if( !low_conn_nums.empty() )  MPI_Pack(&low_conn_nums[0] ,static_cast<INMOST_MPI_SIZE>(num)             ,INMOST_MPI_DATA_INTEGER_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			if( !high_conn_size.empty() ) MPI_Pack(&high_conn_size[0],static_cast<INMOST_MPI_SIZE>(selems[3].size()),INMOST_MPI_DATA_ENUM_TYPE   ,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			if( !high_conn_nums.empty() ) MPI_Pack(&high_conn_nums[0],static_cast<INMOST_MPI_SIZE>(num_high)        ,INMOST_MPI_DATA_INTEGER_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			buffer.resize(position);
		}

        /////////////////////////////////////////
        // pack esets
        //if( false )
		{
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size(selems[4].size());
			std::vector<INMOST_DATA_ENUM_TYPE> high_conn_size(selems[4].size()); // TODO - 3
			std::vector<Storage::integer> low_conn_nums;  // array composed elements : ElementType and position in array 
			std::vector<int> high_conn_nums(selems[4].size() * 3);        // array of indexes of children, sibling, parent. -1 if has't
			INMOST_DATA_ENUM_TYPE num_high = 0;
			position = static_cast<int>(buffer.size());
			new_size = 0;
			num = 0; k = 0;
			int marked_for_data = 0, marked_shared = 0;
            int names_buff_size = 0; 

            // Pack sequence: 
            // 1) all names of sets
            // 2) low conn for sets (arr_pos + element_type)
            // 3) high conn for sets (3 arr_pos: child, sibling, parent.  -1 if has't)
            
            // Compute names_buff_size
			for(element_set::iterator it = selems[4].begin(); it != selems[4].end(); it++) names_buff_size += ElementSet(this,*it).GetName().size() + 1;
            //cout << ro() << mpirank << ": Names buff size = " << names_buff_size << endl;
            std::vector<char> names_buff(names_buff_size);
            int names_buff_pos = 0;

			for(element_set::iterator it = selems[4].begin(); it != selems[4].end(); it++) 
			{
                ElementSet set(this, *it);
                
                // Add set name to names_buffer
                strcpy(&names_buff[names_buff_pos],set.GetName().c_str());
                names_buff_pos += set.GetName().length() + 1;

                // Add all low conns to low_conn_nums
                stringstream ss;
                ss << ro() << mpirank << ": For set " << ElementSet(this,*it).GetName() << " low conns (";
				low_conn_size[k] = 0;
				Element::adj_type const & lc = LowConn(*it);
				for(Element::adj_type::const_iterator jt = lc.begin(); jt != lc.end(); jt++) if( !Hidden(*jt) )
				{
                    INMOST_DATA_INTEGER_TYPE el_num = GetHandleElementNum(*jt);
                    ss << "(" << el_num << "," <<  GlobalID(*jt) << ") ";
			    	assert(Integer(*jt,arr_position) == Integer(selems[el_num][Integer(*jt,arr_position)],arr_position));
                    HandleType composed = ComposeHandle(GetHandleElementType(*jt), Integer(*jt,arr_position));
                    low_conn_nums.push_back(composed);

					low_conn_size[k]++;
					num++;
				}
                ss << ")";
                //cout << ss.str() << endl;

                if (set.HaveChild())   high_conn_nums[k*3+0] = Integer(selems[4][Integer(set.GetChild().GetHandle(),  arr_position)],arr_position); else high_conn_nums[k*3 + 0] = -1;
                if (set.HaveSibling()) high_conn_nums[k*3+1] = Integer(selems[4][Integer(set.GetSibling().GetHandle(),arr_position)],arr_position); else high_conn_nums[k*3 + 1] = -1;
                if (set.HaveParent())  high_conn_nums[k*3+2] = Integer(selems[4][Integer(set.GetParent().GetHandle(), arr_position)],arr_position); else high_conn_nums[k*3 + 2] = -1;

                stringstream ss5;
                ss5 << ro() << mpirank << ": high_conn_nums for set " << set.GetName() << ": ";
                ss5 << high_conn_nums[k*3 + 0] << " " << high_conn_nums[k*3 + 1] << " " << high_conn_nums[k*3 + 2] << endl;
                cout << ss5.str();
            	
                k++;
			}
			/*
            stringstream s1;
            s1 << ro() << mpirank << ": Packed names: ";
            for (int i = 0; i < names_buff_size; i++)
                if (names_buff[i] == '\0')
                    s1 << "|";
                else 
                    s1 << names_buff[i];

            stringstream ss;
            ss << ro() << mpirank << ": packed low_conns_size array: ";
            for (int i = 0; i < num; i++) ss << low_conn_size[i] << " ";
			*/
			MPI_Pack_size(1                                               ,INMOST_MPI_DATA_ENUM_TYPE   ,comm,&temp); new_size += temp;   // count of sets
			MPI_Pack_size(1                                               ,INMOST_MPI_DATA_ENUM_TYPE   ,comm,&temp); new_size += temp;   // names_buff_size
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(names_buff_size)   ,INMOST_MPI_DATA_BULK_TYPE   ,comm,&temp); new_size += temp;   // names_buff
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(selems[4].size())  ,INMOST_MPI_DATA_ENUM_TYPE   ,comm,&temp); new_size += temp; // low_conn_sizes array
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(num)               ,INMOST_MPI_DATA_INTEGER_TYPE,comm,&temp); new_size += temp; // low_conn_nums array
			MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(selems[4].size()*3),INMOST_MPI_DATA_ENUM_TYPE   ,comm,&temp); new_size += temp; // high_conn_nums array
            
			buffer.resize(position+new_size);
			temp = static_cast<INMOST_DATA_ENUM_TYPE>(selems[4].size());
            
			                         MPI_Pack(&temp           ,1                                            ,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			                         MPI_Pack(&names_buff_size,1                                            ,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			if(names_buff_size > 0)  MPI_Pack(&names_buff[0]  ,static_cast<INMOST_MPI_SIZE>(names_buff_size),INMOST_MPI_DATA_BULK_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);

			if( !low_conn_size.empty() )  MPI_Pack(&low_conn_size[0] ,static_cast<INMOST_MPI_SIZE>(selems[4].size())  ,INMOST_MPI_DATA_ENUM_TYPE   ,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			if( !low_conn_nums.empty() )  MPI_Pack(&low_conn_nums[0] ,static_cast<INMOST_MPI_SIZE>(num)               ,INMOST_MPI_DATA_INTEGER_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			if( selems[4].size() > 0 )    MPI_Pack(&high_conn_nums[0],static_cast<INMOST_MPI_SIZE>(selems[4].size()*3),INMOST_MPI_DATA_INTEGER_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			buffer.resize(position);
		}
        /////////////////////////////////////////
		DeleteTag(arr_position);
		for(int i = 3; i >= 0; --i)
			all.insert(all.end(),selems[i].begin(),selems[i].end());
		REPORT_VAL("number of all elements",all.size());
		Storage::integer_array proc = IntegerArrayDV(GetHandle(),tag_processors);
		Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),destination);
		if( ip == proc.end() || (*ip) != destination ) proc.insert(ip,destination);		
		//pack number of tags
		REPORT_STR("prepare elements to pack tag data");
		{
			position = static_cast<int>(buffer.size());
			num = static_cast<INMOST_DATA_ENUM_TYPE>(tag_list.size());
			REPORT_VAL("number of tags",num);
			MPI_Pack_size(1,INMOST_MPI_DATA_ENUM_TYPE,comm,&new_size);
			buffer.resize(position+new_size);
			MPI_Pack(&num,1,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
			buffer.resize(position);
		}
		for(INMOST_DATA_ENUM_TYPE i = 0; i < tag_list.size(); i++)
		{
			{
				assert(HaveTag(tag_list[i]));
				Tag tag = GetTag(tag_list[i]);
				//pack tag name
				{
					std::string name = tag.GetTagName();
					INMOST_DATA_ENUM_TYPE defined, sparse, datatype, datasize;
					datatype = static_cast<INMOST_DATA_ENUM_TYPE>(tag.GetDataType());
					datasize = tag.GetSize();
					defined = NONE;
					sparse = NONE;
					for(ElementType mask = NODE; mask <= MESH; mask = NextElementType(mask))
					{
						if( tag.isDefined(mask) ) defined |= mask;
						if( tag.isSparse(mask) ) sparse |= mask;
					}
					position = static_cast<int>(buffer.size());
					new_size = 0;
					num = static_cast<INMOST_DATA_ENUM_TYPE>(name.size());
					REPORT_VAL("characters in name",num);
					REPORT_VAL("name of tag",name);
					REPORT_VAL("data type",DataTypeName(static_cast<DataType>(datatype)));
					REPORT_VAL("data size",datasize);
#if defined(USE_PARALLEL_WRITE_TIME)
					for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype) )
					{
						if( defined & etype ) REPORT_VAL("defined on",ElementTypeName(etype));
						if( sparse & etype ) REPORT_VAL("sparse on",ElementTypeName(etype));
					}
#endif
					MPI_Pack_size(5                                ,INMOST_MPI_DATA_ENUM_TYPE,comm,&temp); new_size += temp;
					MPI_Pack_size(static_cast<INMOST_MPI_SIZE>(num),MPI_CHAR                 ,comm,&temp); new_size += temp;
					buffer.resize(position+new_size);
					MPI_Pack(&datatype                       ,1                                ,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
					MPI_Pack(&defined                        ,1                                ,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
					MPI_Pack(&sparse                         ,1                                ,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
					MPI_Pack(&datasize                       ,1                                ,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
					MPI_Pack(&num                            ,1                                ,INMOST_MPI_DATA_ENUM_TYPE,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
					MPI_Pack(const_cast<char *>(name.c_str()),static_cast<INMOST_MPI_SIZE>(num),MPI_CHAR                 ,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,comm);
					buffer.resize(position);
				}
				//TODO 46 old
				//PackTagData(GetTag(tag_list[i]),pack_tags,NODE | EDGE | FACE | CELL | ESET,0,buffer);
				PackTagData(tag,selems,destination,NODE | EDGE | FACE | CELL | ESET,pack_tags_mrk,buffer);
				//PackTagData(tag,selems,NODE | EDGE | FACE | CELL | ESET,0,buffer);
                //std::cout << mpirank << " After pack_tag_data\n" << std::endl;

			}
		}
		for(integer i = ElementNum(NODE); i <= ElementNum(CELL); i++) 
			if( !selems[i].empty() ) RemMarkerArray(&selems[i][0],static_cast<enumerator>(selems[i].size()),pack_tags_mrk);
		ReleaseMarker(pack_tags_mrk);
#else
		(void) all;
		(void) buffer;
		(void) destination;
		(void) tag_list;
#endif
		REPORT_VAL("destination",destination);
		EXIT_FUNC();
	}
	
	
	void Mesh::UnpackElementsData(element_set & all, buffer_type & buffer, int source, std::vector<std::string> & tag_recv)
	{
		ENTER_FUNC();
		REPORT_VAL("source",source);
#if defined(USE_MPI)
		int mpirank = GetProcessorRank();
        int rank = mpirank;
		INMOST_DATA_ENUM_TYPE num, temp;
		INMOST_DATA_ENUM_TYPE shift = 0;
		int position = 0;
		elements_by_type selems;
		MarkerType unpack_tags_mrk = CreateMarker();
		//MarkerType swap_normal = 0;
		//if( HaveGeometricData(NORMAL,FACE) )
		//	swap_normal = CreateMarker();
		//TODO 46 old
		//elements_by_type unpack_tags;
		std::vector<HandleType> old_nodes(NumberOfNodes());
		all.clear();
		double time = Timer();
		//TODO 49
		int k = 0;
		for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); it++)
			old_nodes[k++] = *it;
		if( !old_nodes.empty() )
		{
			if( HaveGlobalID(NODE) )
			{
				REPORT_STR("sorted for global id");
				std::sort(old_nodes.begin(),old_nodes.end(),GlobalIDComparator(this));
			}
			else
			{
				REPORT_STR("sorted for centroids");
				std::sort(old_nodes.begin(),old_nodes.end(),CentroidComparator(this));
				//REPORT_VAL("number of old nodes", old_nodes.size());
				//for(std::vector<HandleType>::iterator it = old_nodes.begin(); it != old_nodes.end(); ++it)
				//{
				//	Storage::real_array itcoords = RealArrayDF(*it,CoordsTag());
				//	REPORT_VAL("coords",itcoords[0] << " " << itcoords[1] << " " << itcoords[2]);
				//}
			}
		}
		time = Timer() - time;
		REPORT_STR("gather and sort old nodes");
		REPORT_VAL("time", time);	
		
		time = Timer();
		//unpack nodes
		{
			std::vector<Storage::real> coords;
			std::vector<Storage::integer> global_ids;
			integer dim = GetDimensions();
			REPORT_STR("unpack number of nodes");
            num = 0;
			if( !buffer.empty() ) MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			REPORT_VAL("number of nodes",num);
			REPORT_STR("unpack nodes coordinates");
			coords.resize(num*dim);
			if( !buffer.empty() ) MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&coords[0],static_cast<INMOST_MPI_SIZE>(num*dim),INMOST_MPI_DATA_REAL_TYPE,comm);
			if( HaveGlobalID(NODE) )
			{
				REPORT_STR("unpack nodes global identificators");
				global_ids.resize(num);
				if( !buffer.empty() ) MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&global_ids[0],static_cast<INMOST_MPI_SIZE>(num),INMOST_MPI_DATA_INTEGER_TYPE,comm);
			}
			selems[0].reserve(num);
			//TODO 46 old
			//unpack_tags[0].reserve(num);
			REPORT_STR("creating nodes");
			int found = 0, marked_for_data = 0, marked_ghost = 0;
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				HandleType new_node;
				int find = -1;
				//TODO 49
				if( !old_nodes.empty() )
				{
					if( HaveGlobalID(NODE) )
					{
						//REPORT_STR("searching node by global id");
						//REPORT_VAL("global_id",global_ids[i]);
						std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),global_ids[i],GlobalIDComparator(this));
						if( it != old_nodes.end() ) 
						{
							//REPORT_VAL("lower bound found",GlobalID(*it));
							if( global_ids[i] == GlobalID(*it) )
							{
								find = static_cast<int>(it - old_nodes.begin());
								found++;
							}
						} //else REPORT_STR("lower bound points to end");
						//REPORT_VAL("found",find);
					}
					else 
					{
						//REPORT_STR("searching node by coords");
						//REPORT_VAL("xyz",coords[i*dim+0] << " " << coords[i*dim+1] << " " << coords[i*dim+2]);
						std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),&coords[i*dim],CentroidComparator(this));
						if( it != old_nodes.end() ) 
						{
							Storage::real_array itcoords = RealArrayDF(*it,CoordsTag());
							//REPORT_VAL("lower bound found",itcoords[0] << " " << itcoords[1] << " " << itcoords[2]);
							if( CentroidComparator(this).Compare(&coords[i*dim],itcoords.data()) == 0 )
							{
								find = static_cast<int>(it - old_nodes.begin());
								found++;
							}
						} //else REPORT_STR("lower bound points to end");
						//REPORT_VAL("found",find);
					}
				}
				if( find != -1 ) 
				{
					new_node = old_nodes[find];
					if( IntegerDF(new_node,tag_owner) != mpirank ) 
					{
						SetMarker(new_node,unpack_tags_mrk);
						++marked_for_data;
					}
					//TODO 46 old
					//unpack_tags[0].push_back(new_node);
				}
				else
				{
					new_node = CreateNode(&coords[i*dim])->GetHandle();
					{
						assert(GetStatus(new_node) == Element::Any);
						//REPORT_STR(Element::StatusName(GetStatus(new_node)));
						//REPORT_VAL("node handle",new_node);
						SetStatus(new_node,Element::Ghost);
						IntegerDF(new_node,tag_owner) = -1;
						//TODO 46 old
						//unpack_tags[0].push_back(new_node);
						SetMarker(new_node,unpack_tags_mrk);
						++marked_for_data;
						++marked_ghost;
					}
				}
				selems[0].push_back(new_node);
//#if defined(USE_PARALLEL_WRITE_TIME)
//				{
//					Storage::real_array c = RealArrayDF(new_node,CoordsTag());
//					REPORT_VAL("unpack coords",i << " " << c[0] << " " << c[1] << " " << c[2]);
//				}
//#endif
			}
			REPORT_VAL("total found", found << " / " << selems[0].size());
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[0].size());
			REPORT_VAL("total marked ghost", marked_ghost << " / " << selems[0].size());
		}
		time = Timer() - time;
		REPORT_STR("unpack nodes");
		REPORT_VAL("time", time);	
		
		time = Timer();
		//unpack edges
		{
			ElementArray<Node> e_nodes(this,2);
			shift = 0;
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size;
			std::vector<Storage::integer> low_conn_nums;
			if( !buffer.empty() ) MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			if( num > 0 )
			{
				low_conn_size.resize(num);
				if( num != 0 ) MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&low_conn_size[0],static_cast<INMOST_MPI_SIZE>(num),INMOST_MPI_DATA_ENUM_TYPE,comm);
				temp = 0;
				for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += low_conn_size[i];
				if( temp > 0 )
				{
					low_conn_nums.resize(temp);
					MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&low_conn_nums[0],static_cast<INMOST_MPI_SIZE>(temp),INMOST_MPI_DATA_INTEGER_TYPE,comm);
				}
			}
			selems[1].reserve(num);
			REPORT_VAL("number of edges",num);
			int found = 0, marked_for_data = 0, marked_ghost = 0;
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				e_nodes.resize(low_conn_size[i]);
				//REPORT_VAL("Unpack size",low_conn_size[i]);
				for(INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
				{
					//REPORT_VAL("Unpack node position",low_conn_nums[shift+j]);
					e_nodes.at(j) = selems[0][low_conn_nums[shift+j]];
				}
				//REPORT_VAL("unpack edge nodes ",i << " " << low_conn_nums[shift+0] << " " << low_conn_nums[shift+1]);
				HandleType new_edge = InvalidHandle();
				if( !e_nodes.empty() ) 
					new_edge = FindSharedAdjacency(e_nodes.data(),static_cast<enumerator>(e_nodes.size()));
				if( new_edge == InvalidHandle() )
				{
					new_edge = CreateEdge(e_nodes).first->GetHandle();
					assert(GetStatus(new_edge) == Element::Any);
					//REPORT_STR(Element::StatusName(GetStatus(new_edge)));
					//REPORT_VAL("edge handle",new_edge);
					SetStatus(new_edge,Element::Ghost);
					IntegerDF(new_edge,tag_owner) = -1;
					//TODO 46 old
					//unpack_tags[1].push_back(new_edge);
					SetMarker(new_edge,unpack_tags_mrk);
					++marked_for_data;
					++marked_ghost;

					assert(IntegerArrayDV(new_edge,tag_processors).size() == 0);
				}
				else 
				{
					if( IntegerDF(new_edge,tag_owner) != mpirank ) 
					{
						//TODO 46 old
						//unpack_tags[1].push_back(new_edge);
						SetMarker(new_edge,unpack_tags_mrk);
						++marked_for_data;
					}
					++found;
				}
				selems[1].push_back(new_edge);
				shift+= low_conn_size[i];
			}
			REPORT_VAL("total found", found << " / " << selems[1].size());
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[1].size());
			REPORT_VAL("total marked ghost", marked_ghost << " / " << selems[1].size());
		}
		time = Timer() - time;
		REPORT_STR("unpack edges");
		REPORT_VAL("time", time);
		time = Timer();
		//unpack faces
		{
			ElementArray<Edge> f_edges(this);
			shift = 0;
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size;
			std::vector<Storage::integer> low_conn_nums;
			if( !buffer.empty() ) MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			if( num > 0 )
			{
				low_conn_size.resize(num);
				if( !buffer.empty() && !low_conn_size.empty() ) MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&low_conn_size[0],static_cast<INMOST_MPI_SIZE>(num),INMOST_MPI_DATA_ENUM_TYPE,comm);
				temp = 0;
				for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += low_conn_size[i];
				if( temp > 0 )
				{
					low_conn_nums.resize(temp);
					MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&low_conn_nums[0],static_cast<INMOST_MPI_SIZE>(temp),INMOST_MPI_DATA_INTEGER_TYPE,comm);
				}
			}
			selems[2].reserve(num);
			REPORT_VAL("number of faces",num);
			int found = 0, marked_for_data = 0, marked_ghost = 0;
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				f_edges.resize(low_conn_size[i]);
				//REPORT_VAL("Unpack size",low_conn_size[i]);
				for(INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
				{
					//REPORT_VAL("Unpack edge position",low_conn_nums[shift+j]);
					f_edges.at(j) = selems[1][low_conn_nums[shift+j]];
				}
				HandleType new_face = InvalidHandle();
				if( !f_edges.empty() ) 
					new_face = FindSharedAdjacency(f_edges.data(),static_cast<enumerator>(f_edges.size()));
				if( new_face == InvalidHandle() )
				{
					new_face = CreateFace(f_edges).first->GetHandle();
					assert(GetStatus(new_face) == Element::Any);
					SetStatus(new_face,Element::Ghost);
					IntegerDF(new_face,tag_owner) = -1;
					//TODO 46 old
					//unpack_tags[2].push_back(new_face);
					SetMarker(new_face,unpack_tags_mrk);
					++marked_for_data;
					++marked_ghost;
				} 
				else 
				{
					//current solution - no communication of geometric data
					//here we should check tat remote face orientation match
					//existing face orientation, otherwise we would need
					//to swap normal if it is precomputed and comes from remote
					//processor the worong way
					if( IntegerDF(new_face,tag_owner) != mpirank ) 
					{
						//TODO 46 old
						//unpack_tags[2].push_back(new_face);
						SetMarker(new_face,unpack_tags_mrk);
						++marked_for_data;
						
//						if( swap_normal )//&& GetMarker(new_face,unpack_tags_mrk) )
//						{
//							const Element::adj_type & lc = LowConn(new_face);
//							assert(lc.size() > 1 );
//							assert(lc.size() == f_edges.size());
//							if(lc.size() != f_edges.size() )
//								std::cout << "This should not happen" << std::endl;
//							for(ElementArray<Edge>::size_type q = 0; q < f_edges.size(); ++q)
//							{
//								if( f_edges.at(q) == lc[0] )
//								{
//									if( f_edges.at((q+1)%f_edges.size()) != lc[1] )
//										SetMarker(new_face,swap_normal);
//									break;
//								}
//							}
//						}
						
					}
					++found;
				}
				
				selems[2].push_back(new_face);
				shift+= low_conn_size[i];
			}
			REPORT_VAL("total found", found << " / " << selems[2].size());
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[2].size());
			REPORT_VAL("total marked ghost", marked_ghost << " / " << selems[2].size());
		}
		time = Timer() - time;
		REPORT_STR("unpack faces");
		REPORT_VAL("time", time);
		time = Timer();
		//unpack cells
		{
			ElementArray<Face> c_faces(this);
			ElementArray<Node> c_nodes(this);
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size;
			std::vector<INMOST_DATA_ENUM_TYPE> high_conn_size;
			std::vector<Storage::integer> low_conn_nums;
			std::vector<Storage::integer> high_conn_nums;
			INMOST_DATA_ENUM_TYPE shift_high = 0;
			shift = 0;
			MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			if( num > 0 )
			{
				low_conn_size.resize(num);
				if( num != 0 ) MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&low_conn_size[0],static_cast<INMOST_MPI_SIZE>(num),INMOST_MPI_DATA_ENUM_TYPE,comm);
				temp = 0;
				for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += low_conn_size[i];
				if( temp > 0 )
				{
					low_conn_nums.resize(temp);
					MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&low_conn_nums[0],static_cast<INMOST_MPI_SIZE>(temp),INMOST_MPI_DATA_INTEGER_TYPE,comm);
				}
				high_conn_size.resize(num);
				if( num != 0 ) MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&high_conn_size[0],static_cast<INMOST_MPI_SIZE>(num),INMOST_MPI_DATA_ENUM_TYPE,comm);
				temp = 0;
				for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += high_conn_size[i];
				if( temp > 0 )
				{
					high_conn_nums.resize(temp);
					MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&high_conn_nums[0],static_cast<INMOST_MPI_SIZE>(temp),INMOST_MPI_DATA_INTEGER_TYPE,comm);
				}
			}
			selems[3].reserve(num);
			REPORT_VAL("number of cells",num);
			int found = 0, marked_for_data = 0, marked_ghost = 0;
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				//REPORT_VAL("Unpack low size",low_conn_size[i]);
				c_faces.resize(low_conn_size[i]);
				for(INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
				{
					//REPORT_VAL("Unpack face position",low_conn_nums[shift+j]);
					c_faces.at(j) = selems[2][low_conn_nums[shift+j]];
				}
				//REPORT_VAL("Unpack high size",low_conn_size[i]);
				c_nodes.resize(high_conn_size[i]);
				for(INMOST_DATA_ENUM_TYPE j = 0; j < high_conn_size[i]; j++)
				{
					//REPORT_VAL("Unpack node position",high_conn_nums[shift_high+j]);
					c_nodes.at(j) = selems[0][high_conn_nums[shift_high+j]];
				}
				HandleType new_cell = InvalidHandle();
				if( !c_faces.empty() ) 
					new_cell = FindSharedAdjacency(c_faces.data(),static_cast<enumerator>(c_faces.size()));
				if( new_cell == InvalidHandle() )
				{
					new_cell = CreateCell(c_faces,c_nodes).first->GetHandle();
					assert(GetStatus(new_cell) == Element::Any);
					SetStatus(new_cell,Element::Ghost);
					IntegerDF(new_cell,tag_owner) = -1;
					//TODO 46 old
					//unpack_tags[3].push_back(new_cell);
					SetMarker(new_cell,unpack_tags_mrk);
					++marked_for_data;
					++marked_ghost;
				} 
				else 
				{
					if( IntegerDF(new_cell,tag_owner) != mpirank ) 
					{
						//TODO 46 old
						//unpack_tags[3].push_back(new_cell);
						SetMarker(new_cell,unpack_tags_mrk);
						++marked_for_data;
					}
					++found;
				}
				selems[3].push_back(new_cell);
				shift += low_conn_size[i];
				shift_high += high_conn_size[i];
			}
			REPORT_VAL("total found", found << " / " << selems[3].size());
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[3].size());
			REPORT_VAL("total marked ghost", marked_ghost << " / " << selems[3].size());
		}
		time = Timer() - time;
		REPORT_STR("unpack cells");
		REPORT_VAL("time", time);

        /////////////////////////////////////////////////////////////
        //unpack esets
        //if( false )
		{
			ElementArray<Face> c_faces(this);
			ElementArray<Node> c_nodes(this);
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size;
			std::vector<Storage::integer> low_conn_nums;
			std::vector<int> high_conn_nums;
			INMOST_DATA_ENUM_TYPE shift_high = 0;
			INMOST_DATA_ENUM_TYPE names_buff_size  = 0;
			shift = 0;
            
            // Count of esets
			MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
        //    cout << ro() << rank << ": Unpack num - " << num << endl;
            
            high_conn_nums.resize(num*3);
            
            // Count of chars in names buffer
			MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&names_buff_size,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
          //  cout << ro() << rank << ": Unpack " << names_buff_size << " names buff size" << endl;
            std::vector<char> names_buff(names_buff_size);

            // Names buffer
	        if( names_buff_size != 0 ) MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&names_buff[0],static_cast<INMOST_MPI_SIZE>(names_buff_size),INMOST_MPI_DATA_BULK_TYPE,comm);

            // Gather sets names to array
            int pos = 0;
            std::vector<string> names;
            while (pos < names_buff_size)
            {
                names.push_back(string(&names_buff[pos]));
                pos += names[names.size() - 1].length() + 1;
            }
            /*   
            stringstream ss;
            ss << ro() << rank << ": unpacked names: ";
            for (int i = 0; i < names.size(); i++)
                ss << names[i] << " ";
            */
            //cout << ss.str() << endl;

            // Looking to all names and create the sets if it's needed
            for (int i = 0; i < names.size(); i++)
            {
                ElementSet set = GetSet(names[i]);
                if (set == InvalidElementSet())
                    set = CreateSetUnique(names[i]).first;
                selems[4].push_back(set.GetHandle());
            }

			if( num > 0 )
			{
				low_conn_size.resize(num);
                // Unpack low_conn_size array
				if( num != 0 ) MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&low_conn_size[0],static_cast<INMOST_MPI_SIZE>(num),INMOST_MPI_DATA_ENUM_TYPE,comm);
                
                stringstream ss;
                ss << ro() << rank << ": unpack low_conns_size array: ";
                for (int i = 0; i < num; i++) ss << low_conn_size[i] << " ";
              //  cout << ss.str() << endl;

                temp = 0;
				for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++) temp += low_conn_size[i];
                // Unpack low_conn_num array
				if( temp > 0 )
				{
					low_conn_nums.resize(temp);
					MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&low_conn_nums[0],static_cast<INMOST_MPI_SIZE>(temp),INMOST_MPI_DATA_INTEGER_TYPE,comm);
				}

                // Add low conns for sets
                stringstream ss1;
                ss1 << ro() << rank << ": unpack low_conns_nums array: " << endl;
                int ind = 0;
                for (int i = 0; i < num; i++) 
                {
                    ss1 << ro();
                    for (int j = 0; j < low_conn_size[i]; j++)
                    {
                        INMOST_DATA_ENUM_TYPE type = GetHandleElementNum(low_conn_nums[ind]);
                        int array_pos = GetHandleID(low_conn_nums[ind]);
                        Element elem = Element(this,selems[type][array_pos]);

                        assert(type == GetHandleElementNum(elem.GetHandle()));
                        ss1 << "(" << type  << "," << array_pos << ") ";

                        ElementSet set(this, selems[4][i]);
                        set.PutElement(elem);

                        ind++;
                    }
                    ss1 << " | ";
                }
                //cout << ss1.str() << endl;

                
                // Unpack high conn array
                if (num > 0)
                    MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&high_conn_nums[0],static_cast<INMOST_MPI_SIZE>(num*3),INMOST_MPI_DATA_INTEGER_TYPE,comm);


                stringstream ss6;
                ss6 << ro() << rank << ": unpack high conn nums: ";
                for (int i = 0; i < num*3; i++)
                {
                    if (i%3 == 0) ss6 << "|";
                    ss6 << high_conn_nums[i] << " ";
                }
                //cout << ss6.str() << endl;

                for (int i = 0; i < num; i++) if (high_conn_nums[i*3+0] != -1) ElementSet(this,selems[4][i]).AddChild(ElementSet(this,selems[4][high_conn_nums[i*3+0]]));
                for (int i = 0; i < num; i++) if (high_conn_nums[i*3+1] != -1) ElementSet(this,selems[4][i]).AddSibling(ElementSet(this,selems[4][high_conn_nums[i*3+1]]));
            }
		}
        /////////////////////////////////////////////////////////////
        //cout << rank << "UNPACK ESET COMPLETE" << endl;
		
		
		time = Timer();
		all.reserve(selems[0].size()+selems[1].size()+selems[2].size()+selems[3].size());
		for(int i = 3; i >= 0; --i)
			all.insert(all.end(),selems[i].begin(),selems[i].end());
		REPORT_VAL("number of all elements",all.size());
		Storage::integer_array proc = IntegerArrayDV(GetHandle(),tag_processors);
		Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),source);
		if( ip == proc.end() || (*ip) != source ) proc.insert(ip,source);
		time = Timer() - time;
		REPORT_STR("prepare elements to unpack tag data");
		REPORT_VAL("time", time);	
		time = Timer();
		//unpack tags
		{	
			MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&num,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
			REPORT_VAL("number of tags",num);
			for(INMOST_DATA_ENUM_TYPE i = 0; i < num; i++)
			{
				INMOST_DATA_ENUM_TYPE defined, sparse, datatype, datasize;
				std::string tag_name;
				MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&datatype,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
				MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&defined,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
				MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&sparse,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
				MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&datasize,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
				MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&temp,1,INMOST_MPI_DATA_ENUM_TYPE,comm);
				REPORT_VAL("characters in name",temp);
				tag_name.resize(temp);
				MPI_Unpack(&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),&position,&tag_name[0],temp,MPI_CHAR,comm);
				REPORT_VAL("name of tag",tag_name);
				REPORT_VAL("data type",DataTypeName(static_cast<DataType>(datatype)));
				REPORT_VAL("data size",datasize);
				for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype) )
				{
					if( defined & etype ) REPORT_VAL("defined on",ElementTypeName(etype));
					if( sparse & etype ) REPORT_VAL("sparse on",ElementTypeName(etype));
				}
				tag_recv.push_back(tag_name);
				Tag tag = CreateTag(tag_name,static_cast<enum DataType>(datatype),static_cast<ElementType>(defined),static_cast<ElementType>(sparse),datasize);
				//TODO 46 old
				//UnpackTagData(tag,unpack_tags,0,NODE | EDGE | FACE | CELL | ESET, buffer,position,DefaultUnpack);
				UnpackTagData(tag,selems,NODE | EDGE | FACE | CELL | ESET,unpack_tags_mrk, buffer,position,DefaultUnpack);
				//UnpackTagData(tag,selems,NODE | EDGE | FACE | CELL | ESET,0, buffer,position,DefaultUnpack);
			}
			
		}
		for(integer k = ElementNum(NODE); k <= ElementNum(CELL); ++k)
			if( !selems[k].empty() ) RemMarkerArray(&selems[k][0],static_cast<enumerator>(selems[k].size()),unpack_tags_mrk);
		ReleaseMarker(unpack_tags_mrk);
//		if( swap_normal )
//		{
//			for(enumerator q = 0; q < selems[2].size(); ++q)
//			{
//				if( GetMarker(selems[2][q],swap_normal) )
//				{
//					real_array nrm = RealArrayDF(selems[2][q],GetGeometricTag(NORMAL));
//					std::cout << "Swap normal face " << GetHandleID(selems[2][q]) << " rank " << GetProcessorRank() <<  std::endl;
//					for(real_array::size_type it = 0; it < nrm.size(); ++it)
//						nrm[it] = -nrm[it];
//					RemMarker(selems[2][q],swap_normal);
//				}
//				else
//					std::cout << "No swap normal face " << GetHandleID(selems[2][q]) << " rank " << GetProcessorRank() <<  std::endl;
//				
//			}
//			ReleaseMarker(swap_normal);
//		}
		time = Timer() - time;
		REPORT_STR("unpack tag data");
		REPORT_VAL("time", time);	
#else
		(void) all;
		(void) buffer;
		(void) source;
		(void) tag_recv;
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
		int max_tag = 32767;
		int flag = 0;
		int * p_max_tag;
#if defined(USE_MPI2)
		MPI_Comm_get_attr(comm,MPI_TAG_UB,&p_max_tag,&flag);
#else //USE_MPI2
		MPI_Attr_get(comm,MPI_TAG_UB,&p_max_tag,&flag);
#endif //USE_MPI2
		if( flag ) max_tag = *p_max_tag;
		recv_reqs.resize(recv_bufs.size());
		send_reqs.resize(send_bufs.size());
		REPORT_VAL("strategy",parallel_strategy);
		if( parallel_strategy == 0 )
		{
			for(i = 0; i < send_bufs.size(); i++)
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (send_bufs[i].first+mpisize+rand_num))%max_tag;
				REPORT_MPI(MPI_Isend(&send_bufs[i].second[0],static_cast<INMOST_MPI_SIZE>(send_bufs[i].second.size()),MPI_PACKED,send_bufs[i].first,mpi_tag,comm,&send_reqs[i]));	
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
				REPORT_MPI(MPI_Irecv(&recv_bufs[i].second[0],static_cast<INMOST_MPI_SIZE>(recv_bufs[i].second.size()),MPI_PACKED,recv_stat.MPI_SOURCE,mpi_tag,comm,&recv_reqs[i]));
			}
		}
		else if( parallel_strategy == 1 )
		{
			INMOST_DATA_BULK_TYPE stub;
			REPORT_VAL("recv bufs size",recv_bufs.size());
			for(i = 0; i < recv_bufs.size(); i++)// if( !recv_bufs[i].second.empty() )
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (mpirank+mpisize+rand_num))%max_tag;
				//mpi_tag = parallel_mesh_unique_id*mpisize*mpisize+recv_bufs[i].first*mpisize+mpirank;
				REPORT_VAL("mpi_tag",mpi_tag);
				REPORT_VAL("size",recv_bufs[i].second.size());
				REPORT_VAL("proc",recv_bufs[i].first);
				REPORT_VAL("empty",recv_bufs[i].second.empty());
				REPORT_MPI(MPI_Irecv(recv_bufs[i].second.empty()?&stub:&recv_bufs[i].second[0],static_cast<INMOST_MPI_SIZE>(recv_bufs[i].second.size()),MPI_PACKED,recv_bufs[i].first,mpi_tag,comm,&recv_reqs[i]));
			}
			REPORT_VAL("send bufs size",send_bufs.size());
			for(i = 0; i < send_bufs.size(); i++) //if( !send_bufs[i].second.empty() )
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (send_bufs[i].first+mpisize+rand_num))%max_tag;
				//mpi_tag = parallel_mesh_unique_id*mpisize*mpisize+mpirank*mpisize+send_bufs[i].first;
				REPORT_VAL("mpi_tag",mpi_tag);
				REPORT_VAL("size",send_bufs[i].second.size());
				REPORT_VAL("proc",send_bufs[i].first);
				REPORT_VAL("empty",send_bufs[i].second.empty());
				REPORT_MPI(MPI_Isend(send_bufs[i].second.empty()?&stub:&send_bufs[i].second[0],static_cast<INMOST_MPI_SIZE>(send_bufs[i].second.size()),MPI_PACKED,send_bufs[i].first,mpi_tag,comm,&send_reqs[i]));	
			}
		}
		else if( parallel_strategy == 2 )
		{
			INMOST_DATA_BULK_TYPE stub;
			REPORT_VAL("recv bufs size",recv_bufs.size());
			for(i = 0; i < recv_bufs.size(); i++) //if( !recv_bufs[i].second.empty() )
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (mpirank+mpisize+rand_num))%max_tag;
				REPORT_MPI(MPI_Irecv(recv_bufs[i].second.empty()? &stub : &recv_bufs[i].second[0],static_cast<INMOST_MPI_SIZE>(recv_bufs[i].second.size()),MPI_PACKED,recv_bufs[i].first,mpi_tag,comm,&recv_reqs[i]));
			}
			REPORT_MPI(MPI_Barrier(comm));
			REPORT_VAL("send bufs size",send_bufs.size());
			for(i = 0; i < send_bufs.size(); i++)// if( !send_bufs[i].second.empty() )
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (send_bufs[i].first+mpisize+rand_num))%max_tag;
				REPORT_MPI(MPI_Irsend(send_bufs[i].second.empty()? &stub : &send_bufs[i].second[0],static_cast<INMOST_MPI_SIZE>(send_bufs[i].second.size()),MPI_PACKED,send_bufs[i].first,mpi_tag,comm,&send_reqs[i]));	
			}
		}
#else //USE_MPI
		(void) send_bufs;
		(void) recv_bufs;
		(void) send_reqs;
		(void) recv_reqs;
#endif
		EXIT_FUNC();
	}
	
	void Mesh::PrepareReceiveInner(Prepare todo, exch_buffer_type & send_bufs, exch_buffer_type & recv_bufs)
	{

		if( parallel_strategy == 0 && todo == UnknownSize ) return; //in this case we know all we need
		ENTER_FUNC();
#if defined(USE_MPI)
        int mpirank = GetProcessorRank();
#if defined(USE_MPI_P2P) && defined(PREFFER_MPI_P2P)
		unsigned i, end = send_bufs.size();
        REPORT_MPI(MPI_Win_fence(MPI_MODE_NOPRECEDE,window)); //start exchange session
		memset(shared_space,0,sizeof(unsigned)*mpisize); //zero bits where we receive data
		REPORT_MPI(MPI_Win_fence(0,window)); //wait memset finish
		for(i = 0; i < end; i++) shared_space[mpisize+i] = send_bufs[i].second.size()+1; //put data to special part of the memory
		for(i = 0; i < end; i++) REPORT_MPI(MPI_Put(&shared_space[mpisize+i],1,MPI_UNSIGNED,send_bufs[i].first,mpirank,1,MPI_UNSIGNED,window)); //request rdma
		REPORT_MPI(MPI_Win_fence(MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED,window)); //end exchange session
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
#else //PREFFER_MPI_P2P
		if( todo == UnknownSize )
		{
			REPORT_STR("Unknown size");
			if( parallel_strategy != 0 )
			{
				REPORT_VAL("exchange number", ++num_exchanges);
				unsigned i;
				int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber(), rand_num = randomizer.Number()+1;
				int mpi_tag;
				int max_tag = 32767;
				int flag = 0;
				int * p_max_tag;
#if defined(USE_MPI2)
				MPI_Comm_get_attr(comm,MPI_TAG_UB,&p_max_tag,&flag);
#else //USE_MPI2
				MPI_Attr_get(comm,MPI_TAG_UB,&p_max_tag,&flag);
#endif //USE_MPI2
				if( flag ) max_tag = *p_max_tag;
				REPORT_VAL("max_tag",max_tag);
				std::vector<int> send_recv_size(send_bufs.size()+recv_bufs.size());
				std::vector<INMOST_MPI_Request> reqs(send_bufs.size()+recv_bufs.size());
				for(i = 0; i < send_bufs.size(); i++)
					send_recv_size[i+recv_bufs.size()] = static_cast<int>(send_bufs[i].second.size());
				REPORT_VAL("recv buffers size",recv_bufs.size());
				for(i = 0; i < recv_bufs.size(); i++)
				{
					mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (mpirank+mpisize+rand_num))%max_tag;
					REPORT_VAL("origin",recv_bufs[i].first);
					REPORT_VAL("mpi_tag",mpi_tag);
					//mpi_tag = parallel_mesh_unique_id*mpisize*mpisize+recv_bufs[i].first*mpisize+mpirank;
					REPORT_MPI(MPI_Irecv(&send_recv_size[i],1,MPI_INT,recv_bufs[i].first,mpi_tag,comm,&reqs[i]));
				}
				REPORT_VAL("send buffers size",send_bufs.size());
				for(i = 0; i < send_bufs.size(); i++)
				{
					mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + (send_bufs[i].first+mpisize+rand_num))%max_tag;
					REPORT_VAL("destination",send_bufs[i].first);
					REPORT_VAL("mpi_tag",mpi_tag);
					REPORT_VAL("size",send_recv_size[i+recv_bufs.size()]);
					//mpi_tag = parallel_mesh_unique_id*mpisize*mpisize+mpirank*mpisize+send_bufs[i].first;
					REPORT_MPI(MPI_Isend(&send_recv_size[i+recv_bufs.size()],1,MPI_INT,send_bufs[i].first,mpi_tag,comm,&reqs[i+recv_bufs.size()]));	
				}
				if( !recv_bufs.empty() )
				{
					REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(recv_bufs.size()),&reqs[0],MPI_STATUSES_IGNORE));
				}
				REPORT_VAL("recieved buffers size",recv_bufs.size());
				for(i = 0; i < recv_bufs.size(); i++)
				{
					REPORT_VAL("origin",recv_bufs[i].first);
					REPORT_VAL("size",send_recv_size[i]);
					recv_bufs[i].second.resize(send_recv_size[i]);
				}
				if( !send_bufs.empty() )
				{
					REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(send_bufs.size()),&reqs[recv_bufs.size()],MPI_STATUSES_IGNORE));
				}
			}
		}
		else if( todo == UnknownSource )
		{
			REPORT_STR("Unknown source");
#if defined(USE_MPI_P2P)
			int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber();
			unsigned i, end = static_cast<unsigned>(send_bufs.size());
            REPORT_MPI(MPI_Win_fence(MPI_MODE_NOPRECEDE,window)); //start exchange session
			memset(shared_space,0,sizeof(unsigned)*mpisize); //zero bits where we receive data
			//REPORT_MPI(MPI_Win_fence( MPI_MODE_NOPRECEDE,window)); //start exchange session
            REPORT_MPI(MPI_Win_fence( 0,window)); //wait memset finish
			for(i = 0; i < end; i++) shared_space[mpisize+i] = static_cast<unsigned>(send_bufs[i].second.size()+1); //put data to special part of the memory
			for(i = 0; i < end; i++) 
			{
                REPORT_VAL("put value", shared_space[mpisize+i]);
                REPORT_VAL("destination", send_bufs[i].first);
                REPORT_VAL("displacement", mpirank);
				REPORT_MPI(MPI_Put(&shared_space[mpisize+i],1,MPI_UNSIGNED,send_bufs[i].first,mpirank,1,MPI_UNSIGNED,window)); //request rdma to target processors for each value
			}
			REPORT_MPI(MPI_Win_fence(MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED,window)); //end exchange session
			if( parallel_strategy == 0 )
			{
				unsigned num = 0;
				for(int ii = 0; ii < mpisize; ii++) if( shared_space[ii] > 0 ) num++;
				recv_bufs.resize(num);
			}
			else if( todo == UnknownSize )
			{
				end = static_cast<unsigned>(recv_bufs.size());
				for(i = 0; i < end; i++)
					recv_bufs[i].second.resize(shared_space[recv_bufs[i].first]-1);
			}
			else if( todo == UnknownSource )
			{
				recv_bufs.clear();
				for(int ii = 0; ii < mpisize; ii++)
					if( shared_space[ii] > 0 )
                    {
                        REPORT_VAL("position", ii);
                        REPORT_VAL("value", shared_space[ii]);
						recv_bufs.push_back(proc_buffer_type(ii,std::vector<INMOST_DATA_BULK_TYPE>(shared_space[ii]-1))); // this call would be optimized by compiler
                    }
                REPORT_VAL("recvs",recv_bufs.size());
			}
#else //USE_MPI_P2P
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
				for(int i = 0; i < static_cast<int>(send_bufs.size()); i++)
				{
					sends_dest_and_size[i*2+0] = send_bufs[i].first;
					sends_dest_and_size[i*2+1] = static_cast<unsigned>(send_bufs[i].second.size());
				}
				unsigned recvsize = 0;
				int k,j;
				std::vector<int> allsize(mpisize);
				int size = static_cast<int>(send_bufs.size()*2);
				REPORT_MPI(MPI_Allgather(&size,1,MPI_INT,&allsize[0],1,MPI_INT,comm));
				std::vector<int> displs(mpisize+1,0);
				for(k = 0; k < mpisize; k++)
					recvsize += allsize[k];
				for(k = 1; k < mpisize+1; k++)
					displs[k] = displs[k-1]+allsize[k-1];
				std::vector<unsigned> recvs_dest_and_size(recvsize+1);
				REPORT_MPI(MPI_Allgatherv(&sends_dest_and_size[0],static_cast<INMOST_MPI_SIZE>(send_bufs.size()*2),MPI_UNSIGNED,&recvs_dest_and_size[0],&allsize[0],&displs[0],MPI_UNSIGNED,comm));
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
#endif //USE_MPI_P2P
		}
#endif //PREFFER_MPI_P2P
#else //USE_MPI
		(void) send_bufs;
		(void) recv_bufs;
#endif //USE_MPI
		EXIT_FUNC();
	}
	
	std::vector<int> Mesh::FinishRequests(std::vector<INMOST_MPI_Request> & recv_reqs)
	{
        std::vector<int> ret;
		ENTER_FUNC();
#if defined(USE_MPI)
		int outcount = 0;
		REPORT_VAL("requests",recv_reqs.size());
		if( recv_reqs.empty() ) 
			ret.clear();
		else
		{
			ret.resize(recv_reqs.size(),-1);
			if( !recv_reqs.empty() )
			{
				REPORT_MPI(MPI_Waitsome(static_cast<INMOST_MPI_SIZE>(recv_reqs.size()),&recv_reqs[0],&outcount,&ret[0],MPI_STATUSES_IGNORE));
			}
			else outcount = MPI_UNDEFINED;
			if( outcount == MPI_UNDEFINED ) ret.clear();
			else ret.resize(outcount);
		}
#else
		(void) recv_reqs;
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
		proc_elements send_elements;
		std::vector<int> done;

		if( action == AGhost ) REPORT_STR("Ghosting algorithm")
		else if( action == AMigrate ) REPORT_STR("Migration algorithm")
		
		ListTagNames(tag_list);
		{
			std::vector<std::string>::iterator it = tag_list.begin();
			while( it != tag_list.end() )
			{
				if( it->substr(0,9) == "PROTECTED" ) 
					it = tag_list.erase(it);
                else if(GetTag(*it).GetDataType() == DATA_REFERENCE)
					it = tag_list.erase(it);
				else it++;
			}
		}
		{
			REPORT_STR("Gathering elements to send");
			for(ElementType etype = NODE; etype <= ESET; etype = NextElementType(etype))
			for(iteratorElement it = BeginElement(etype); it != EndElement(); it++)
//			for(Mesh::iteratorElement it = BeginElement(CELL | FACE | EDGE | NODE); it != EndElement(); ++it) if( it->HaveData(tag_sendto) )
			{
				Storage::integer_array mark = it->IntegerArray(tag_sendto);
				for(Storage::integer_array::iterator kt = mark.begin(); kt != mark.end(); kt++)
					if( *kt != mpirank )
                    {
						send_elements[*kt].push_back(*it);
                    }
				it->DelData(tag_sendto);
			}
		}
		num_wait = 0;
		send_bufs.resize(send_elements.size());
        //std::cout << mpirank << ": Send size: " << send_elements.size() << std::endl;
		REPORT_STR("Packing elements to send");
		for(proc_elements::iterator it = send_elements.begin(); it != send_elements.end(); it++)
			if( !it->second.empty() )
			{
				REPORT_VAL("pack for", it->first);
				REPORT_VAL("number of provided elements",it->second.size());
//                std::cout << mpirank << "Number of provided els " << it->second.size() << std::endl;
				PackElementsData(it->second,send_bufs[num_wait].second,it->first,tag_list);
				REPORT_VAL("number of got elements",it->second.size());
				send_bufs[num_wait].first = it->first;
				num_wait++;
			}
		send_bufs.resize(num_wait);	

		//NEW ALGORITHM
		if( false )
		{ //in most cases the exchange of elements will effect only nearest neighbours, we want to detect this
			Storage::integer_array p = IntegerArrayDV(GetHandle(),tag_processors);
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
					for(integer_array::size_type k = 0; k < p.size(); k++) recv_bufs[k].first = p[k];
				}
				PrepareReceiveInner(UnknownSize, send_bufs,recv_bufs);
			}
			else PrepareReceiveInner(UnknownSource, send_bufs,recv_bufs);
		}
		else PrepareReceiveInner(UnknownSource, send_bufs,recv_bufs);
		ExchangeBuffersInner(send_bufs,recv_bufs,send_reqs,recv_reqs);
			
		if( action == AMigrate ) // delete packed elements
		{
			REPORT_STR("Gathering elements to delete");
			Tag tag_new_processors = GetTag("TEMPORARY_NEW_PROCESSORS");
			{
				elements_by_type delete_elements;
				MarkerType delete_marker = CreateMarker();
				for(proc_elements::iterator it = send_elements.begin(); it != send_elements.end(); it++)
				{
					REPORT_VAL("for processor",it->first);
					REPORT_VAL("number of elements",it->second.size());
					if( !it->second.empty() )
					{
						//Have packed all entities, now delete those entities that i should not own	
						for(element_set::iterator kt = it->second.begin(); kt != it->second.end(); kt++)
							if( !GetMarker(*kt,delete_marker) )
							{
								Storage::integer_array procs = IntegerArray(*kt,tag_new_processors);
								if( !std::binary_search(procs.begin(),procs.end(),mpirank) )
								{
									delete_elements[GetHandleElementNum(*kt)].push_back(*kt);
									SetMarker(*kt,delete_marker);
								}
							}
					}
				}
				//Should delete elements in order by CELL..NODE
				REPORT_STR("Deleting elements");
				REPORT_VAL("number of nodes to delete",delete_elements[0].size());
				REPORT_VAL("number of edges to delete",delete_elements[1].size());
				REPORT_VAL("number of faces to delete",delete_elements[2].size());
				REPORT_VAL("number of cells to delete",delete_elements[3].size());
				for(int j = ElementNum(CELL); j >= ElementNum(NODE); j--)
				{
					//this is not really needed because we destroy those elements
					//for(element_set::iterator kt = delete_elements[j].begin(); kt != delete_elements[j].end(); ++kt) RemMarker(*kt,delete_marker);
					for(element_set::iterator kt = delete_elements[j].begin(); kt != delete_elements[j].end(); ++kt) Destroy(*kt);
				}
				ReleaseMarker(delete_marker);
			}
			REPORT_STR("Deleting some other not owned elements");
			//now delete elements that we have not yet deleted
			int count = 0;
			for(ElementType etype = NODE; etype <= CELL; etype = NextElementType(etype))
			for(iteratorElement it = BeginElement(etype); it != EndElement(); it++)
			//for(iteratorElement it = BeginElement(CELL | FACE | EDGE | NODE); it != EndElement(); it++)
			{
				Storage::integer_array procs = it->IntegerArray(tag_new_processors);
				if( !std::binary_search(procs.begin(),procs.end(),mpirank) ) 
				{
					Destroy(*it);
					++count;
				}
			}
			REPORT_VAL("number of some other",count);
			REPORT_STR("Done deleting");
		}
		send_elements.clear();
		REPORT_STR("Unpacking received data");
		while( !(done = FinishRequests(recv_reqs)).empty() )
		{
			for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
			{
				element_set recv_elements;
				REPORT_STR("call unpack");
				UnpackElementsData(recv_elements,recv_bufs[*qt].second,recv_bufs[*qt].first,tag_list_recv);
				if( action == AGhost )
				{
					for(element_set::iterator it = recv_elements.begin(); it != recv_elements.end(); it++)
					{
						Storage::integer owner = IntegerDF(*it,tag_owner);
						if( owner == mpirank ) continue;
						Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);	
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

		REPORT_STR("Ended recv");
				
		if( action == AGhost ) //second round to inform owner about ghosted elements
		{
			REPORT_STR("Second round for ghost exchange");
			if( !send_reqs.empty() )
			{
				REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(send_reqs.size()),&send_reqs[0],MPI_STATUSES_IGNORE));
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
			for(proc_elements::iterator it = send_elements.begin(); it != send_elements.end(); it++)
				if( !it->second.empty() )
				{
					PackElementsData(it->second,send_bufs[num_wait].second,it->first,tag_list);
					send_bufs[num_wait].first = it->first;
					num_wait++;
				}
			send_bufs.resize(num_wait);
				
			if( false )
			{ //in most cases the exchange of elements will effect only nearest neighbours, we want to detect this
				Storage::integer_array p = IntegerArrayDV(GetHandle(),tag_processors);
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
						for(integer_array::size_type k = 0; k < p.size(); k++) recv_bufs[k].first = p[k];
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
					element_set recv_elements;
					UnpackElementsData(recv_elements,recv_bufs[*qt].second,recv_bufs[*qt].first,tag_list_recv);	
					for(element_set::iterator it = recv_elements.begin(); it != recv_elements.end(); it++)
					{
						Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);
						Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),recv_bufs[*qt].first);
						if( ip == proc.end() || (*ip) != recv_bufs[*qt].first ) proc.insert(ip,recv_bufs[*qt].first);
						if( IntegerDF(*it,tag_owner) == mpirank )
							SetStatus(*it,Element::Shared);
						else
							SetStatus(*it,Element::Ghost);
					}
				}
			}
		}
		
		if( action == AMigrate ) //Compute new values
		{
			REPORT_STR("Second round for elements migration");
			REPORT_STR("Computing new values");
			Tag tag_new_owner = GetTag("TEMPORARY_NEW_OWNER");
			Tag tag_new_processors = GetTag("TEMPORARY_NEW_PROCESSORS");	
			for(ElementType etype = NODE; etype <= CELL; etype = NextElementType(etype))
			for(iteratorElement it = BeginElement(etype); it != EndElement(); it++)
			//for(iteratorElement it = BeginElement(CELL | FACE | EDGE | NODE); it != EndElement(); it++)
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
						SetStatus(*it,Element::Owned);
					else
						SetStatus(*it,Element::Shared);
				}
				else SetStatus(*it,Element::Ghost);
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
			ExchangeData(tag_processors,CELL | FACE | EDGE | NODE,0);
		}
		ComputeSharedProcs();
		
		if( action == AMigrate ) 
		{
			if( !send_reqs.empty() )
			{
				REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(send_reqs.size()),&send_reqs[0],MPI_STATUSES_IGNORE));
			}
		}
#else
		(void) action;
#endif
		EXIT_FUNC();
	}
	
	void Mesh::RecomputeParallelStorage(ElementType mask)
	{
		ENTER_FUNC();
#if defined(USE_MPI) && defined(USE_PARALLEL_STORAGE)
		for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
			for(int i = 0; i < 4; i++) if( mask & ElementTypeFromDim(i) )
				it->second[i].clear();
		for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)			
			for(int i = 0; i < 4; i++) if( mask & ElementTypeFromDim(i) )
				it->second[i].clear();
		GatherParallelStorage(ghost_elements,shared_elements,mask);
#else //USE_MPI and USE_PARALLEL_STORAGE
		(void) mask;
#endif //USE_MPI and USE_PARALLEL_STORAGE
		EXIT_FUNC();
	}
	
	void Mesh::SortParallelStorage(ElementType mask)
	{
#if defined(USE_PARALLEL_STORAGE)
		SortParallelStorage(ghost_elements,shared_elements,mask);
#else
		(void)mask;
#endif
	}

	void Mesh::SortParallelStorage(parallel_storage & ghost, parallel_storage & shared, ElementType mask)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		REPORT_STR("sort shared elements")
		double time = Timer();
		for(parallel_storage::iterator it = shared.begin(); it != shared.end(); it++)
		{
			REPORT_VAL("processor",it->first);
			for(int i = 0; i < 4; i++) if( mask & ElementTypeFromDim(i) )
			{
				REPORT_STR(ElementTypeName(ElementTypeFromDim(i)) << " have global id " << (!tag_global_id.isValid() ? "invalid" : (tag_global_id.isDefined(ElementTypeFromDim(i))?"YES":"NO")));
				if( !it->second[i].empty() )
				{
					if(HaveGlobalID(ElementTypeFromDim(i)))
						std::sort(it->second[i].begin(),it->second[i].end(),GlobalIDComparator(this));
					else
						std::sort(it->second[i].begin(),it->second[i].end(),CentroidComparator(this));
				}
				REPORT_VAL(ElementTypeName(mask & ElementTypeFromDim(i)),it->second[i].size());
			}

            // ESET sort
			if (mask & ElementTypeFromDim(4))
            {
				if( !it->second[4].empty() )
                {
			        std::sort(it->second[4].begin(),it->second[4].end(), SetNameComparator(this));
                }
            }
		}
		time = Timer() - time;
		REPORT_VAL("time",time);
		REPORT_STR("sort ghost elements")
		time = Timer();
		for(parallel_storage::iterator it = ghost.begin(); it != ghost.end(); it++)
		{
			REPORT_VAL("processor",it->first);
			for(int i = 0; i < 4; i++) if( mask & ElementTypeFromDim(i) )
			{
				REPORT_STR(ElementTypeName(ElementTypeFromDim(i)) << " have global id " << (!tag_global_id.isValid() ? "invalid" : (tag_global_id.isDefined(ElementTypeFromDim(i))?"YES":"NO")));
				if( !it->second[i].empty() )
				{
					if( HaveGlobalID(ElementTypeFromDim(i)) )
						std::sort(it->second[i].begin(),it->second[i].end(),GlobalIDComparator(this));
					else
						std::sort(it->second[i].begin(),it->second[i].end(),CentroidComparator(this));
				}
				REPORT_VAL(ElementTypeName(mask & ElementTypeFromDim(i)),it->second[i].size());
			}

            // ESET sort
			if (mask & ElementTypeFromDim(4))
            {
				if( !it->second[4].empty() )
                {
			        std::sort(it->second[4].begin(),it->second[4].end(), SetNameComparator(this));
                }
            }
		}
		time = Timer() - time;
		REPORT_VAL("time",time);
#else //USE_MPI and USE_PARALLEL_STORAGE
		(void) ghost; (void) shared; (void) mask;
#endif //USE_MPI and USE_PARALLEL_STORAGE
		EXIT_FUNC();
	}
	
	void Mesh::ExchangeGhost(Storage::integer layers, ElementType bridge)
	{
    //printf("%d called exchange ghost with layers %d bridge %s\n",GetProcessorRank(), layers,ElementTypeName(bridge));
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
		int mpirank = GetProcessorRank();
		bool delete_ghost = false;
		//if( layers == Integer(tag_layers) && bridge == Integer(tag_bridge) ) return;
		if( layers < Integer(GetHandle(),tag_layers) ) delete_ghost = true;
		else if( layers == Integer(GetHandle(),tag_layers) && bridge < Integer(GetHandle(),tag_bridge) ) delete_ghost = true;
		int test_bridge = 0;
		
		if( (bridge & MESH) || (bridge & ESET) || (bridge & CELL) ) throw Impossible;
		for(ElementType mask = NODE; mask <= FACE; mask = NextElementType(mask))
			test_bridge += (bridge & mask)? 1 : 0;
		if( test_bridge == 0 || test_bridge > 1 ) throw Impossible;
		double time;
		//RemoveGhost();
		Tag layers_marker = CreateTag("TEMPORARY_LAYERS_MARKER",DATA_INTEGER,CELL,CELL);
		Integer(GetHandle(),tag_layers) = layers;
		Integer(GetHandle(),tag_bridge) = bridge;
		Storage::integer_array procs = IntegerArrayDV(GetHandle(),tag_processors);
		proc_elements old_layers;
		proc_elements current_layers;
		element_set all_visited;
		{
			
			proc_elements shared_skin = ComputeSharedSkinSet(bridge);
      //printf("%d shared skin size %d\n",GetProcessorRank(),shared_skin.size());
			time = Timer();
			for(Storage::integer_array::iterator p = procs.begin(); p != procs.end(); p++)
			{
				MarkerType busy = CreateMarker();
				all_visited.clear();
				for(element_set::iterator it = shared_skin[*p].begin(); it != shared_skin[*p].end(); it++)
				{
					ElementArray<Element> adj = Element(this,*it)->getAdjElements(CELL);
					for(ElementArray<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
						if( !jt->GetMarker(busy) )
						{
							//if( jt->IntegerDF(tag_owner) != *p )
							if( jt->IntegerDF(tag_owner) == mpirank )
							{
								current_layers[*p].push_back(*jt);
								if( layers > 0 )
								{
									Storage::integer_array adj_procs = jt->IntegerArrayDV(tag_processors);
									if( !std::binary_search(adj_procs.begin(),adj_procs.end(),*p) )
										jt->IntegerArray(tag_sendto).push_back(*p);
									if( delete_ghost ) jt->IntegerArray(layers_marker).push_back(*p);
								}
							}
							jt->SetMarker(busy);
							all_visited.push_back(*jt);
						}
				}
				if( !all_visited.empty() ) RemMarkerArray(&all_visited[0],static_cast<enumerator>(all_visited.size()),busy);
				//for(element_set::iterator it = all_visited.begin(); it != all_visited.end(); ++it) RemMarker(*it,busy);
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
					element_set & ref_cur = current_layers[*p];
					element_set & ref_old = old_layers[*p];
					MarkerType busy = CreateMarker();
					all_visited.clear();
					for(element_set::iterator it = ref_old.begin(); it != ref_old.end(); ++it) SetMarker(*it,busy);
					for(element_set::iterator it = ref_old.begin(); it != ref_old.end(); it++)
					{
						ElementArray<Element> adj_bridge = Element(this,*it)->getAdjElements(bridge);
						for(ElementArray<Element>::iterator jt = adj_bridge.begin(); jt != adj_bridge.end(); jt++)
							if( !jt->GetMarker(busy) )
							{
								ElementArray<Element> adj = jt->getAdjElements(CELL);
								for(ElementArray<Element>::iterator kt = adj.begin(); kt != adj.end(); kt++)
									if( !kt->GetMarker(busy) )
									{
										if( jt->IntegerDF(tag_owner) != *p )
										{
											ref_cur.push_back(*kt);
											Storage::integer_array adj_procs = kt->IntegerArrayDV(tag_processors);
											if( !std::binary_search(adj_procs.begin(),adj_procs.end(),*p) )
												kt->IntegerArray(tag_sendto).push_back(*p);
											if( delete_ghost ) jt->IntegerArray(layers_marker).push_back(*p);
										}
										kt->SetMarker(busy);
										all_visited.push_back(*kt);
									}
								jt->SetMarker(busy);
								all_visited.push_back(*jt);
							}
					}
					for(element_set::iterator it = ref_old.begin(); it != ref_old.end(); ++it) RemMarker(*it,busy);
					for(element_set::iterator it = all_visited.begin(); it != all_visited.end(); ++it) RemMarker(*it,busy);
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
			ReduceData(layers_marker,CELL,0,UnpackLayersMarker);
			ExchangeData(layers_marker,CELL,0);
			ElementArray<Element> del_ghost(this);
			for(iteratorElement it = BeginElement(CELL); it != EndElement(); it++)
			{
				if( GetStatus(*it) == Element::Ghost )
				{
					if( !it->HaveData(layers_marker) )
						del_ghost.push_back(*it);
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
							del_ghost.push_back(*it);
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
#else
		(void) layers;
		(void) bridge;
#endif
		EXIT_FUNC();

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
		ElementType bridge = Integer(GetHandle(),tag_bridge);
		Storage::integer layers = Integer(GetHandle(),tag_layers);
		
		ExchangeData(tag_new_owner,CELL,0);
		
		
		double time = Timer();
		
		
		for(iteratorElement it = BeginElement(CELL); it != EndElement(); it++)
			it->IntegerArrayDV(tag_new_processors).push_back(it->Integer(tag_new_owner));
		
		dynarray<Storage::integer,64> result,intersection;
		//determine tag_new_processors for FACEs to calculate shared skin
		for(ElementType mask = FACE; mask >= NODE; mask = PrevElementType(mask))
		{
#if defined(USE_PARALLEL_WRITE_TIME)
			std::map<int,int> numelems;
#endif
			for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
			{
				Storage::integer_array procs = it->IntegerArrayDV(tag_new_processors);
				determine_my_procs_high(this,*it,tag_new_processors,result,intersection);
				if( !result.empty() )
				{
					it->Integer(tag_new_owner) = result[0];
					procs.replace(procs.begin(),procs.end(),result.begin(),result.end());
#if defined(USE_PARALLEL_WRITE_TIME)
					numelems[result[0]]++;
#endif
				}
				else
				{
					it->Integer(tag_new_owner) = mpirank;
					procs.clear();
					procs.push_back(mpirank);
#if defined(USE_PARALLEL_WRITE_TIME)
					numelems[mpirank]++;
#endif
				}
			}
#if defined(USE_PARALLEL_WRITE_TIME)
			for(std::map<int,int>::iterator it = numelems.begin(); it != numelems.end(); ++it)
			{
				REPORT_VAL(ElementTypeName(mask),it->second);
				REPORT_STR("elements on lower ierarhy belong to");
				REPORT_VAL("processor",it->first);
			}
#endif
		}	
		time = Timer() - time;
		REPORT_STR("Determine new processors");
		REPORT_VAL("time",time);
		ExchangeData(tag_new_owner,FACE | EDGE | NODE,0);
		time = Timer();
		if( bridge != NONE && layers != 0 )
		{
			REPORT_STR("Multiple layers");
			proc_elements redistribute_skin, skin_faces;
			element_set all_visited;
		
			ReduceData(tag_new_processors,FACE,0,RedistUnpack);
			ExchangeData(tag_new_processors,FACE,0);
			
			double time2 = Timer();
			//compute skin around migrated cells to compute ghost layer
			for(iteratorElement it = BeginElement(FACE); it != EndElement(); it++)
			{
				Storage::integer_array procs = it->IntegerArray(tag_new_processors);
				if( procs.size() == 2 ) //there may be only 2 procs for every face max, because we have unique redistribution for every cell at the beginning
				{
					skin_faces[procs[0]].push_back(*it);
					skin_faces[procs[1]].push_back(*it);
				}
			}
			if( bridge == FACE ) redistribute_skin.swap(skin_faces);
			else
			{
				MarkerType busy = CreateMarker();
				// Here should first find all adjacent elements of skin_faces of given type,
				// then distribute them between sets according to new processors
				for(proc_elements::iterator it = skin_faces.begin(); it != skin_faces.end(); it++)
				{
					element_set & ref = redistribute_skin[it->first];
					for(element_set::iterator jt = it->second.begin(); jt != it->second.end(); jt++)
					{
						ElementArray<Element> adj = Element(this,*jt)->getAdjElements(bridge);
						for(ElementArray<Element>::iterator kt = adj.begin(); kt != adj.end(); kt++)
							if( !kt->GetMarker(busy) )
							{
								ref.push_back(*kt);
								kt->SetMarker(busy);
							}
					}
					for(element_set::iterator jt = ref.begin(); jt != ref.end(); ++jt) RemMarker(*jt,busy);
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
				proc_elements current_layers,old_layers;
				MarkerType busy = CreateMarker();
				for(proc_elements::iterator it = redistribute_skin.begin(); it != redistribute_skin.end(); it++)
				{
					all_visited.clear();
					for(element_set::iterator jt = it->second.begin(); jt != it->second.end(); jt++)
					{
						element_set & ref = current_layers[it->first];
						ElementArray<Element> adj = Element(this,*jt)->getAdjElements(CELL);
						for(ElementArray<Element>::iterator kt = adj.begin(); kt != adj.end(); kt++)
							if( !kt->GetMarker(busy) )
							{
								Storage::integer_array procs = kt->IntegerArray(tag_new_processors);
								Storage::integer_array::iterator find = std::lower_bound(procs.begin(),procs.end(),it->first);
								if( find == procs.end() || *find != it->first ) 
								{
									procs.insert(find,it->first);
									ref.push_back(*kt);
								}
								all_visited.push_back(*kt);
								kt->SetMarker(busy);	
							}
					}
					for(element_set::iterator it = all_visited.begin(); it != all_visited.end(); ++it) RemMarker(*it,busy);
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
					for(proc_elements::iterator qt = old_layers.begin(); qt != old_layers.end(); qt++)
					{
						MarkerType busy = CreateMarker();
						all_visited.clear();
						for(element_set::iterator it = qt->second.begin(); it != qt->second.end(); it++) SetMarker(*it,busy);
						for(element_set::iterator it = qt->second.begin(); it != qt->second.end(); it++)
						{
							element_set & ref = current_layers[qt->first];
							ElementArray<Element> adj_bridge = Element(this,*it)->getAdjElements(bridge);
							for(ElementArray<Element>::iterator jt = adj_bridge.begin(); jt != adj_bridge.end(); jt++)
								if( !jt->GetMarker(busy) )
								{
									ElementArray<Element> adj = jt->getAdjElements(CELL);
									for(ElementArray<Element>::iterator kt = adj.begin(); kt != adj.end(); kt++)
									{
										if( !kt->GetMarker(busy) )
										{
											Storage::integer_array procs = kt->IntegerArray(tag_new_processors);
											Storage::integer_array::iterator find = std::lower_bound(procs.begin(),procs.end(),qt->first);
											if( find == procs.end() || *find != qt->first )
											{
												procs.insert(find,qt->first);
												ref.push_back(*kt);
											}
											kt->SetMarker(busy);
											all_visited.push_back(*kt);
										}
									}
									jt->SetMarker(busy);
									all_visited.push_back(*jt);
								}
						}
						for(element_set::iterator it = qt->second.begin(); it != qt->second.end(); it++) RemMarker(*it,busy);
						for(element_set::iterator it = all_visited.begin(); it != all_visited.end(); it++) RemMarker(*it,busy);
						ReleaseMarker(busy);
					}
					time2 = Timer() - time2;
					REPORT_STR("Mark layer");
					REPORT_VAL("layer",k);
					REPORT_VAL("time",time2);
				}
			}
			
			ReduceData(tag_new_processors,CELL,0,RedistUnpack);
			ExchangeData(tag_new_processors,CELL,0);
			
			time3 = Timer() - time3;
			REPORT_STR("Mark all layers");
			REPORT_VAL("time",time3);
			
			time2 = Timer();
			
			for(ElementType mask = FACE; mask >= NODE; mask = PrevElementType(mask))
			{
				for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
				{
					Storage::integer_array procs = it->IntegerArrayDV(tag_new_processors);
					determine_my_procs_high(this,*it,tag_new_processors,result,intersection);
					if( result.empty() ) 
					{
						procs.clear();
						procs.push_back(mpirank);
					}
					else procs.replace(procs.begin(),procs.end(),result.begin(),result.end());
				}
			}
			
			
			
			ReduceData(tag_new_processors,FACE| EDGE| NODE,0,RedistUnpack);
			ExchangeData(tag_new_processors,FACE|EDGE|NODE,0);
			time2 = Timer() - time2;
			REPORT_STR("Detect processors");
			REPORT_VAL("time",time2);
		}
		else 
		{
			REPORT_STR("No layers");
			ReduceData(tag_new_processors,FACE | EDGE | NODE,0,RedistUnpack);
			ExchangeData(tag_new_processors,FACE | EDGE | NODE,0);
		}
		
		time = Timer() - time;
		REPORT_STR("Determine new processors for layers");
		REPORT_VAL("time",time);
		
		
		time = Timer();


		for(ElementType etype = NODE; etype <= CELL; etype = NextElementType(etype))
		for(iteratorElement it = BeginElement(etype); it != EndElement(); it++)
		//for(iteratorElement it = BeginElement(CELL | FACE | EDGE | NODE); it != EndElement(); it++)
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
		
		
		//std::stringstream str;
		//str << "before_migrate" << GetProcessorRank() << ".xml";
		//Save(str.str());
		
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
	
	
	void UnpackMarkNormalOrientation(const Tag & tag,
								 const Element & element,
								 const INMOST_DATA_BULK_TYPE * data,
								 INMOST_DATA_ENUM_TYPE size)
	{
		Storage::real_array    local_nrm = element.RealArrayDF(tag);
		const Storage::real * remote_nrm = (const Storage::real *) data;
		assert(size == local_nrm.size());
		Storage::real dot = 0;
		for(Storage::enumerator k = 0; k < local_nrm.size(); ++k)
			dot += local_nrm[k]*remote_nrm[k];
		local_nrm[0] = dot;
        (void)size;
	}
	
	void Mesh::MarkNormalOrientation(MarkerType mrk)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,FACE);
#endif //USE_PARALLEL_STORAGE
		TagRealArray tag_nrm = CreateTag("TEMPORARY_NORMAL",DATA_REAL,FACE,NONE,3);
		for(iteratorFace it = BeginFace(); it != EndFace(); ++it)
			if( it->GetStatus() & (Element::Ghost | Element::Shared) )
				it->UnitNormal(tag_nrm[it->self()].data());
		exchange_data storage;
		ExchangeDataInnerBegin(tag_set(1,tag_nrm),shared_elements,ghost_elements,FACE,0,storage);
		ExchangeDataInnerEnd(tag_set(1,tag_nrm),shared_elements,ghost_elements,FACE,0,UnpackMarkNormalOrientation,storage);
		for(iteratorFace it = BeginFace(); it != EndFace(); ++it)
			if( it->GetStatus() == Element::Ghost && tag_nrm[it->self()][0] < 0 )
				it->SetMarker(mrk);
		DeleteTag(tag_nrm);
#else //USE_MPI
		(void) mrk;
#endif //USE_MPI
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

    void Mesh::ResolveSets()
    {
#ifdef USE_MPI
        int mpirank = GetProcessorRank();
        int mpisize = GetProcessorsNumber();

        map<string, vector<int> > map_names; // key - set_name, value - array of processors ranks which has this set

        // Collect all set names to vector<string>
        //vector<string>* all_set_names = new vector<string>[mpisize];
        //std::vector<std::string>& set_names = all_set_names[mpirank];
        std::vector<std::string> set_names;
        int bytes_size = 0;
        for(Mesh::iteratorSet it = BeginSet(); it != EndSet(); ++it)
        {
            set_names.push_back(it->GetName());
            map_names[it->GetName()].push_back(mpirank);
            bytes_size += it->GetName().length() + 1; // +1 for terminator
        }

        // Gather size of buffers from other processors
        int* all_bytes_sizes = new int[mpisize];
        all_bytes_sizes[mpirank] = bytes_size;

        for (int p = 0; p < mpisize; p++)
        {
            if (mpirank == p)
                MPI_Bcast(&bytes_size,1,MPI_INTEGER,p,GetCommunicator());
            else
                MPI_Bcast(&all_bytes_sizes[p],1,MPI_INTEGER,p,GetCommunicator());
        }

/*
        stringstream ss;
        for (int i = 0; i < mpisize; i++)
        {
            ss << all_bytes_sizes[i] << " ";
        }
        cout << mpirank << ": " << ss.str() << endl;
*/

        // Configure buffer for send to other processors: name1/0 name2/0 ...
        char* buffer = new char[bytes_size];
        char* ptr = buffer;
        for (int i = 0; i < set_names.size(); i++)
        {
            strcpy(ptr,set_names[i].c_str());
            ptr += set_names[i].length() + 1;
        }

        // send/receive buffers to/from other processors
        char* recv_buffer;
        for (int p = 0; p < mpisize; p++)
        {
            if (mpirank == p)
            {
                MPI_Bcast(buffer,all_bytes_sizes[p],MPI_CHAR,p,GetCommunicator());
            }
            else
            {
                recv_buffer = new char[all_bytes_sizes[p]];
                MPI_Bcast(recv_buffer,all_bytes_sizes[p],MPI_CHAR,p,GetCommunicator());


                // after receive - configure map_names
                ptr = recv_buffer;
                while (ptr-recv_buffer < all_bytes_sizes[p])
                {
                    int length = strlen(ptr);
                    char* set_name = new char[length+1];
                    strcpy(set_name,ptr);
                    //all_set_names[p].push_back(std::string(set_name));
                    map_names[set_name].push_back(p);
                    ptr += length+1;
                }
                delete[] recv_buffer;
            }
        }


        /*
        stringstream sss;
        sss << mpirank << endl;
        for (int i = 0; i < mpisize; i++)
        {
            sss << "Proc " << i << ": ";
            for (int j = 0; j < all_set_names[i].size(); j++)
            {
                sss << all_set_names[i][j] << " ";
            }
            sss << endl;
        }
        //cout << sss.str() << endl;

        stringstream s3;
        for (map<string,vector<int> >::iterator name = map_names.begin(); name != map_names.end(); name++)
        {
            s3 << name->first << " ";
            for (int i = 0; i < name->second.size(); i++)
                s3 << name->second[i] << " ";
            s3 << endl;

        }

        //cout << sss.str() << endl << s3.str() << endl;
        */
        
       // Change status for self sets
       for(Mesh::iteratorSet set = BeginSet(); set != EndSet(); ++set)
       {
           string set_name = set->GetName();

           Storage::integer_array arr = set->IntegerArrayDV(tag_processors);
           arr.resize(map_names[set_name].size());
           for (int i = 0; i < map_names[set_name].size(); i++)
               arr[i] = map_names[set_name][i];   

           assert(map_names[set_name].size() > 0);

           if (map_names[set_name].size() == 1)
           {
               assert(map_names[set_name][0] == mpirank);
               SetStatus(set->GetHandle(), Element::Owned);
			   set->IntegerDF(tag_owner) = mpirank;
               continue;
           }

           int min = map_names[set_name][0];
           for (int i = 1; i < map_names[set_name].size(); i++)
               if (map_names[set_name][i] < min)
                   min = map_names[set_name][i];

		   set->IntegerDF(tag_owner) = min;
           if (min == mpirank)
               SetStatus(set->GetHandle(), Element::Shared);
           else
               SetStatus(set->GetHandle(), Element::Ghost);
       }

       GatherParallelStorage(ghost_elements, shared_elements, ESET);

        
       /* 
       if (mpirank == 0)
           for (parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
           {
               cout << it->first << ": Sets in parallel storage: ";
               for(element_set::iterator p = it->second[4].begin(); p != it->second[4].end(); p++) 
               {
                   cout << ElementSet(this,*p).GetName() << " ";
               }
               cout << endl;
           }
        */
#endif

    }
}


#endif
