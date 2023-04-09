#ifdef _MSC_VER //kill some warnings
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#endif

#include "inmost.h"

//#define DEBUG_UNPACK

#if defined(USE_MESH)
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
//using namespace std;

#if defined(USE_PARALLEL_STORAGE)
bool allow_pack_by_gids = true;
#else // USE_PARALLEL_STORAGE
bool allow_pack_by_gids = false;
#endif // USE_PARALLEL_STORAGE


#if defined(USE_MPI)
static INMOST_DATA_BIG_ENUM_TYPE pmid = 0;
#endif // USE_MPI
#if defined(USE_PARALLEL_WRITE_TIME)
__INLINE std::string NameSlash(std::string input)
{
	for(size_t l = input.size(); l > 0; --l)
		if( input[l-1] == '/' || input[l-1] == '\\' )
			return std::string(input.c_str() + l);
	return input;
}
#define REPORT_MPI(x) {WriteTab(out_time) << "<MPI><![CDATA[" << #x << "]]></MPI>" << std::endl; x;}
#define REPORT_STR(x) {WriteTab(out_time) << "<TEXT><![CDATA[" << x << "]]></TEXT>" << std::endl;}
#define REPORT_VAL(str,x) {WriteTab(out_time) << "<VALUE name=\"" << str << "\"> <CONTENT><![CDATA[" << x << "]]></CONTENT> <CODE><![CDATA[" << #x << "]]></CODE></VALUE>" << std::endl;}
#define ENTER_FUNC() double all_time = Timer(); {WriteTab(out_time) << "<FUNCTION name=\"" << __FUNCTION__ << "\" id=\"func" << func_id++ << "\">" << std::endl; Enter();}
#define ENTER_BLOCK() { double btime = Timer(); WriteTab(out_time) << "<FUNCTION name=\"" << __FUNCTION__ << ":" << NameSlash(__FILE__) << ":" << __LINE__ << "\" id=\"func" << GetFuncID()++ << "\">" << std::endl; Enter();
#define EXIT_BLOCK() WriteTab(out_time) << "<TIME>" << Timer() - btime << "</TIME>" << std::endl; Exit(); WriteTab(out_time) << "</FUNCTION>" << std::endl;}
#define EXIT_FUNC() {WriteTab(out_time) << "<TIME>" << Timer() - all_time << "</TIME>" << std::endl; Exit(); WriteTab(out_time) << "</FUNCTION>" << std::endl;}
#define EXIT_FUNC_DIE() {WriteTab(out_time) << "<TIME>" << -1 << "</TIME>" << std::endl; Exit(); WriteTab(out_time) << "</FUNCTION>" << std::endl;}
#else // USE_PARALLEL_WRITE_TIME
#define REPORT_MPI(x) x
#define REPORT_STR(x) {}
#define REPORT_VAL(str,x) {}
#define ENTER_FUNC() {}
#define ENTER_BLOCK()
#define EXIT_BLOCK()
#define EXIT_FUNC() {}
#define EXIT_FUNC_DIE()  {}
#endif // USE_PARALLEL_WRITE_TIME


#if defined(__LINUX__) || defined(__linux__) || defined(__APPLE__)
#include <unistd.h>
#define PROCESSID getpid()
#else
#define PROCESSID -1
#endif


#if 1
#define MPI_Pack_call(data,size,type,buf,buf_size,pos,comm)   MPI_Pack(data,size,type,buf,buf_size,pos,comm)
#define MPI_Pack_size_call(size,type,comm,ret)                MPI_Pack_size(size,type,comm,ret)
#define MPI_Unpack_call(buf,buf_size,pos,data,size,type,comm) MPI_Unpack(buf,buf_size,pos,data,size,type,comm)
#define POS_TYPE int
#define SEND_AS MPI_PACKED
#else
#define MPI_Pack_call(data,size,type,buf,buf_size,pos,comm)   MPI_Pack_external("external32",data,size,type,buf,buf_size,pos)
#define MPI_Pack_size_call(size,type,comm,ret)                MPI_Pack_external_size("external32",size,type,ret)
#define MPI_Unpack_call(buf,buf_size,pos,data,size,type,comm) MPI_Unpack_external("external32",buf,buf_size,pos,data,size,type)
#define POS_TYPE MPI_Aint
#define SEND_AS MPI_BYTE
#endif

#ifndef SIZE_MAX
# ifdef __SIZE_MAX__
#  define SIZE_MAX __SIZE_MAX__
# else
#  define SIZE_MAX (static_cast<size_t>(-1))
# endif
#endif

#if SIZE_MAX == UCHAR_MAX
#define MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
#error "Cannot detect size_t type"
#endif


namespace INMOST
{
#if defined(USE_MPI)
	template<typename T> struct MPIType;
	template<> struct MPIType<bool> {static MPI_Datatype Type() {return MPI_C_BOOL;}};
	template<> struct MPIType<char> {static MPI_Datatype Type() {return MPI_CHAR;}};
	template<> struct MPIType<const char> {static MPI_Datatype Type() {return MPI_CHAR;}};
	template<> struct MPIType<int> {static MPI_Datatype Type() {return MPI_INT;}};
	template<> struct MPIType<short> {static MPI_Datatype Type() {return MPI_SHORT;}};
	template<> struct MPIType<long> {static MPI_Datatype Type() {return MPI_LONG;}};
	template<> struct MPIType<long long> {static MPI_Datatype Type() {return MPI_LONG_LONG;}};
	template<> struct MPIType<unsigned char> {static MPI_Datatype Type() {return MPI_UNSIGNED_CHAR;}};
	template<> struct MPIType<unsigned short> {static MPI_Datatype Type() {return MPI_UNSIGNED_SHORT;}};
	template<> struct MPIType<unsigned> {static MPI_Datatype Type() {return MPI_UNSIGNED;}};
	template<> struct MPIType<unsigned long> {static MPI_Datatype Type() {return MPI_UNSIGNED_LONG;}};
	template<> struct MPIType<unsigned long long> {static MPI_Datatype Type() {return MPI_UNSIGNED_LONG_LONG;}};
	template<> struct MPIType<double> {static MPI_Datatype Type() {return MPI_DOUBLE;}};
	template<> struct MPIType<long double> {static MPI_Datatype Type() {return MPI_LONG_DOUBLE;}};
	template<> struct MPIType<float> {static MPI_Datatype Type() {return MPI_FLOAT;}};
	//template<> struct MPIType<size_t> {static MPI_Datatype Type() {return MPI_SIZE_T;}};
	template<typename T>
	void pack_data(Mesh::buffer_type & buf, const T & data, INMOST_MPI_Comm comm)
	{
		int ierr;
		POS_TYPE pack_size = 0, shift = 0;
		size_t write_pos = buf.size();
		ierr = MPI_Pack_size_call(1,MPIType<T>::Type(),comm,&pack_size); 
		if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
		buf.resize(write_pos+pack_size);
		ierr = MPI_Pack_call((void*)&data,1,MPIType<T>::Type(),&buf[write_pos],pack_size,&shift,comm);
		if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
		buf.resize(write_pos+shift);
	}
	
	template<typename T>
	size_t pack_data_size(INMOST_MPI_Comm comm)
	{
		int ierr;
		POS_TYPE pack_size = 0;//, shift = 0;
		size_t ret = 0;
		ierr = MPI_Pack_size_call(1,MPIType<T>::Type(),comm,&pack_size); 
		if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
		ret += pack_size;
		return ret;
	}
	
	
	template<typename T>
	void pack_data_array(Mesh::buffer_type & buf, const T * data, size_t N, INMOST_MPI_Comm comm)
	{
		int ierr;
		if( N )
		{
			POS_TYPE pack_size = 0;
			ierr = MPI_Pack_size_call(1,MPIType<T>::Type(),comm,&pack_size);
			if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
			size_t rec_bytes = static_cast<size_t>(pack_size);
			size_t max_bytes = rec_bytes*N;
			size_t chunk_bytes = std::min(max_bytes,static_cast<size_t>(INT_MAX));
			size_t chunk_size = chunk_bytes / rec_bytes;
			size_t offset = 0;
			while( offset != N )
			{
				size_t write_pos = buf.size();
				size_t chunk = std::min(chunk_size, N-offset);
				ierr = MPI_Pack_size_call(static_cast<POS_TYPE>(chunk),MPIType<T>::Type(),comm,&pack_size); //pack_size is expected to be within INT_MAX
				if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
				buf.resize(write_pos+pack_size);
				POS_TYPE shift = 0;
				ierr = MPI_Pack_call((void *)&data[offset],static_cast<POS_TYPE>(chunk),MPIType<T>::Type(),&buf[write_pos],pack_size,&shift,comm);
				if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
				buf.resize(write_pos+shift);
				offset += chunk;
			}
		}
	}
	
	
	template<typename T>
	void pack_data_vector(Mesh::buffer_type & buf, const std::vector<T> & data, INMOST_MPI_Comm comm)
	{
		//int ierr;
		pack_data(buf,data.size(),comm);
		if( !data.empty() ) 
			pack_data_array(buf,&data[0],data.size(),comm);
	}
	
	template<typename T>
	size_t pack_data_array_size(size_t N, INMOST_MPI_Comm comm)
	{
		int ierr;
		size_t ret = 0;
		if( N )
		{
			POS_TYPE pack_size = 0;
			ierr = MPI_Pack_size_call(1,MPIType<T>::Type(),comm,&pack_size);
			if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
			size_t rec_bytes = static_cast<size_t>(pack_size);
			size_t max_bytes = rec_bytes*N;
			size_t chunk_bytes = std::min(max_bytes,static_cast<size_t>(INT_MAX));
			size_t chunk_size = chunk_bytes / rec_bytes;
			size_t offset = 0;
			while( offset != N )
			{
				size_t chunk = std::min(chunk_size, N-offset);
				ierr = MPI_Pack_size_call(static_cast<POS_TYPE>(chunk),MPIType<T>::Type(),comm,&pack_size); //pack_size is expected to be within INT_MAX
				if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
				ret += pack_size;
				offset += chunk;
			}
		}
		return ret;
	}
	
	template<typename T>
	size_t pack_data_vector_size(size_t N, INMOST_MPI_Comm comm)
	{
		size_t ret = pack_data_size<size_t>(comm);
		if( N ) ret += pack_data_array_size<T>(N,comm);
		return ret;
	}
	
	
	template<typename T>
	void unpack_data(Mesh::buffer_type & buf, size_t & buf_pos, T & data, INMOST_MPI_Comm comm)
	{
		int ierr;
		POS_TYPE shift = 0, pack_size;
		ierr = MPI_Pack_size_call(1,MPIType<T>::Type(),comm,&pack_size); 
		if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
		ierr = MPI_Unpack_call((void*)&buf[buf_pos],pack_size,&shift,&data,1,MPIType<T>::Type(),comm);
		if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
		buf_pos += shift;
	}
	
	
	template<typename T>
	void unpack_data_array(Mesh::buffer_type & buf, size_t & buf_pos, T * data, size_t N, INMOST_MPI_Comm comm)
	{
		int ierr;
		if( N )
		{
			POS_TYPE pack_size = 0;
			ierr = MPI_Pack_size_call(1,MPIType<T>::Type(),comm,&pack_size);
			if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
			size_t rec_bytes = static_cast<size_t>(pack_size);
			size_t max_bytes = rec_bytes*N;
			size_t chunk_bytes = std::min(max_bytes,static_cast<size_t>(INT_MAX));
			size_t chunk_size = chunk_bytes / rec_bytes;
			size_t offset = 0;
			while( offset != N )
			{
				size_t chunk = std::min(chunk_size, N-offset);
				ierr = MPI_Pack_size_call(static_cast<POS_TYPE>(chunk),MPIType<T>::Type(),comm,&pack_size); //pack_size is expected to be within INT_MAX
				if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
				POS_TYPE shift = 0;
				ierr = MPI_Unpack_call((void *)&buf[buf_pos],pack_size,&shift,&data[offset],static_cast<POS_TYPE>(chunk),MPIType<T>::Type(),comm);
				if( ierr != MPI_SUCCESS ) MPI_Abort(comm,__LINE__);
				buf_pos += shift;
				offset += chunk;
			}
		}
	}
	
	template<typename T>
	void unpack_data_vector(Mesh::buffer_type & buf, size_t & buf_pos, std::vector<T> & data, INMOST_MPI_Comm comm)
	{
		//int ierr;
		size_t unpack_size;
		unpack_data(buf,buf_pos,unpack_size,comm);
		data.resize(unpack_size);
		if( !data.empty() )
			unpack_data_array(buf,buf_pos,&data[0],data.size(),comm);
	}
#endif
    
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
				else assert(false); //received zero size for dense nonzero tag
			}
			return;
		}
		if( !element->HaveData(tag) )
			element->SetDataSize(tag,size);
		else if( size != element->GetDataSize(tag) )
		{
			if( tag.GetSize() == ENUMUNDEF )
				element->SetDataSize(tag,size);
			else assert(false); //received different size for fixed-sized array
		}
		element->SetData(tag,0,size,data);
	}


	void UnpackOnSkin(const Tag & tag, const Element & e, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size )
		{
			const Storage::integer * recv = static_cast<const Storage::integer *>(static_cast<const void *>(data));
			Storage::integer_array arr = e->IntegerArray(tag);
			//for(int k = 0; k < size; ++k)
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
			INMOST_DATA_ENUM_TYPE old_size;
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
			std::vector<Storage::integer> result(static_cast<size_t>(p1.size())+size);
			result.resize(std::set_union(p1.begin(),p1.end(),p2,p2+size,result.begin())-result.begin());
			p1.replace(p1.begin(),p1.end(),result.begin(),result.end());
		}
	}

	void UnpackLayersMarker(const Tag & tag, const Element & e, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size )
		{
			INMOST_DATA_ENUM_TYPE old_size;
			if( e->HaveData(tag) ) old_size = e->GetDataSize(tag);
			else old_size = 0;
			e->SetDataSize(tag,old_size+size);
			e->SetData(tag,old_size,size,data);
		}
	}
/*
    void Mesh::CheckFaces()
    {
        //std::cout << "Check faces" << endl;

        for(Mesh::iteratorFace it = BeginFace(); it != EndFace(); ++it) 
        {
            std::set<int> set_nodes;
            ElementArray<Node> nodes = it->getNodes();
            bool suc = true;
            for (ElementArray<Node>::iterator node = nodes.begin(); node != nodes.end(); node++)
                if (set_nodes.find(node->LocalID()) != set_nodes.end())
                {
                    suc = false;
                    break;
                }
                else
                {
                    set_nodes.insert(node->LocalID());
                }
            //if (suc) std::cout << "-=== Good face: " << setw(2) << it->LocalID();
            //else std::cout << "-=== Error face: " << setw(2) << it->LocalID();

            //cout << ". Nodes = " << nodes.size() << ": ";
            //for (ElementArray<Node>::iterator node = nodes.begin(); node != nodes.end(); node++)
             //   std::cout << node->LocalID() << " ";
            //std::cout << endl;
        }
    }
*/

	static void OperationMinDistance(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		(void)size;
		assert(size == 2);
		INMOST_DATA_INTEGER_TYPE * idata = (INMOST_DATA_INTEGER_TYPE *)data;
		Storage::integer_array rdata = element->IntegerArray(tag);
		if (idata[1] < rdata[1])
		{
			rdata[0] = idata[0];
			rdata[1] = idata[1];
		}
		else if(idata[1] == rdata[1] && idata[0] < rdata[0])
		{
			rdata[0] = idata[0];
			rdata[1] = idata[1];
		}
	}
	
	static void OperationMinOwner(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		(void)size;
		assert(size == 1);
		INMOST_DATA_INTEGER_TYPE idata = *(INMOST_DATA_INTEGER_TYPE *)data;
		element->Integer(tag) = std::min(element->Integer(tag),idata);
	}
	
	void Mesh::EquilibrateGhost()//bool only_new)
	{
		static int counter = 0;
		if( GetProcessorsNumber() == 1 ) return;
		ENTER_FUNC();
		ENTER_BLOCK();
		CheckGhostSharedCount(__FILE__, __LINE__);
		CheckCentroids(__FILE__, __LINE__);
		CheckOwners(__FILE__, __LINE__);
		CheckGIDs(__FILE__, __LINE__);
		CheckProcessors();
		EXIT_BLOCK();
#if defined(USE_MPI)
		//std::cout << "before_equiv"+std::to_string(counter)+".pvtk" << std::endl;
		//Save("before_equiv"+std::to_string(counter)+".pvtk");
		ElementType bridge = FACE;
		if( tag_bridge.isValid() ) bridge = ElementType(Integer(GetHandle(),tag_bridge));
		if( bridge == NONE ) bridge = FACE; // tag can be valid and uninitialized (for example, after Mesh::SetCommunicator)
		assert(OneType(bridge)); //what if there is multiple types? (todo: the lowest type is required!)
#if 0
		TagIntegerArray tag = CreateTag("TEMPORARY_OWNER_DISTANCE_TAG", DATA_INTEGER, CELL, CELL, 2);
		std::queue<HandleType> cells_queue;
		double t, tadj = 0;
		ENTER_BLOCK();
		REPORT_STR("Collect cells into queue");
		//TODO: compute graph distance only in new elements!!!!
		// Push nearest to owned cell into queue
		for(integer k = 0; k < CellLastLocalID(); k++) if( isValidCell(k) )
		{
			Cell c = CellByLocalID(k);
			if (!c.Hidden() && c.GetStatus() != Element::Owned)
			{
				assert( c.GetStatus() == Element::Shared || c.GetStatus() == Element::Ghost );
				tag[c][0] = -1;
				tag[c][1] = INT_MAX; 
				//ElementArray<Cell> cells = c.NeighbouringCells();
				t = Timer();
				ElementArray<Cell> cells = c.BridgeAdjacencies2Cell(bridge);
				tadj += Timer()-t;
				for(INMOST_DATA_ENUM_TYPE l = 0; l < cells.size(); l++) if (cells[l].GetStatus() == Element::Owned)
				{
					cells_queue.push(c.GetHandle());
					tag[c][0] = cells[l].Integer(tag_owner);
					tag[c][1] = 1;
					break;
				}
			}
		}
		EXIT_BLOCK();
		int cur_dist = 1; // For assert
		ENTER_BLOCK();
		REPORT_STR("graph distance calculation");
		while (!cells_queue.empty())
		{
			Cell c(this,cells_queue.front());
			cells_queue.pop();
			if (cur_dist > tag[c][1]) assert(0); 
			if (cur_dist < tag[c][1]) cur_dist++;
			//ElementArray<Cell> cells = c.NeighbouringCells();
			t = Timer();
			ElementArray<Cell> cells = c.BridgeAdjacencies2Cell(bridge);
			tadj += Timer()-t;
			for(INMOST_DATA_ENUM_TYPE l = 0; l < cells.size(); l++) if (cells[l].GetStatus() != Element::Owned)
			{
				assert( cells[l].GetStatus() == Element::Shared || cells[l].GetStatus() == Element::Ghost );
				if (tag[cells[l]][1] > tag[c][1] + 1)
				{
					tag[cells[l]][0] = tag[c][0];
					tag[cells[l]][1] = tag[c][1] + 1;
					cells_queue.push(cells[l].GetHandle());
				}
			}
		}
		EXIT_BLOCK();
		REPORT_VAL("adjacency requests", tadj);
#else
		TagIntegerArray tag;
		ENTER_BLOCK();
		tag = CreateTag("TEMPORARY_OWNER_DISTANCE_TAG", DATA_INTEGER, CELL | bridge, NONE, 2);
		EXIT_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for schedule(dynamic,1)
#endif
		for (integer k = 0; k < CellLastLocalID(); k++) if (isValidCell(k))
		{
			Cell c = CellByLocalID(k);
			if (!c.Hidden() && c.GetStatus() != Element::Owned)
			{
				tag[c][0] = INT_MAX;
				tag[c][1] = INT_MAX;
			}
		}
		//init bridge elements
		ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for (integer k = 0; k < LastLocalID(bridge); k++) if (isValidElement(bridge,k))
		{
			Element e = ElementByLocalID(bridge, k);
			if (!e.Hidden())
			{
				tag[e][0] = INT_MAX;
				tag[e][1] = INT_MAX;
				ElementArray<Cell> cells = e.getCells();
				Element::Status stat = Element::Any;
				for (ElementArray<Cell>::iterator it = cells.begin(); it != cells.end(); ++it)
					stat |= it->GetStatus();
				if ((stat & Element::Owned) && (stat & (Element::Shared | Element::Ghost)))
				{
					int owner = INT_MAX;
					for (ElementArray<Cell>::iterator it = cells.begin(); it != cells.end(); ++it) if (it->GetStatus() == Element::Owned)
						owner = std::min(owner, it->IntegerDF(tag_owner));
					tag[e][0] = owner;
					tag[e][1] = 0;
				}
			}
		}
		EXIT_BLOCK();
		int cur_dist = 0, incdist = 0;
		ENTER_BLOCK();
		do
		{
			incdist = 0;
			REPORT_VAL("current distance", cur_dist);
			ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for schedule(dynamic,1) reduction(+:incdist)
#endif
			for (integer k = 0; k < CellLastLocalID(); k++) if (isValidCell(k))
			{
				Cell c = CellByLocalID(k);
				if (!c.Hidden() && c.GetStatus() != Element::Owned)
				{
					ElementArray<Element> adj = c.getAdjElements(bridge);
					std::pair<int, int> val(INT_MAX, INT_MAX);
					for (ElementArray<Element>::iterator it = adj.begin(); it != adj.end(); ++it)
						val = std::min(val, std::make_pair(tag[*it][1], tag[*it][0]));
					if (val.first == cur_dist)
					{
						tag[c][0] = val.second;
						tag[c][1] = val.first + 1;
						incdist++;
					}
				}
			}
			EXIT_BLOCK();
			REPORT_VAL("increased distance", incdist);
			cur_dist++;
			ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for (integer k = 0; k < LastLocalID(bridge); k++) if (isValidElement(bridge, k))
			{
				Element e = ElementByLocalID(bridge, k);
				if (!e.Hidden())
				{
					if (tag[e][1] == INT_MAX)
					{
						ElementArray<Cell> cells = e.getCells();
						std::pair<int, int> val(INT_MAX, INT_MAX);
						for (ElementArray<Cell>::iterator it = cells.begin(); it != cells.end(); ++it)
							val = std::min(val, std::make_pair(tag[*it][1], tag[*it][0]));
						tag[e][0] = val.second;
						tag[e][1] = val.first;
					}
				}
			}
			EXIT_BLOCK();
		} while (incdist);
		EXIT_BLOCK();
		ENTER_BLOCK();
		tag = DeleteTag(tag, bridge);
		EXIT_BLOCK();
		ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for schedule(dynamic,1)
#endif
		for (integer k = 0; k < CellLastLocalID(); k++) if (isValidCell(k))
		{
			Cell c = CellByLocalID(k);
			if (!c.Hidden() && c.GetStatus() != Element::Owned)
			{
				if (tag[c][0] == INT_MAX)
					tag[c][0] = -1;
			}
		}
		EXIT_BLOCK();
#endif
		//static int equiv = 0;
		//std::string file;
		//file = "aequiv"+std::to_string(equiv)+".pvtk";
		//std::cout << "Save " << file << std::endl;
		//Save(file);
		
		//CheckCentroids(__FILE__,__LINE__);
		ReduceData(tag, CELL, 0, OperationMinDistance);
		ExchangeData(tag, CELL, 0);
		
		//std::cout << "comput_equiv"+std::to_string(counter)+".pvtk" << std::endl;
		//Save("comput_equiv"+std::to_string(counter)+".pvtk");
		//file = "bequiv"+std::to_string(equiv)+".pvtk";
		//std::cout << "Save " << file << std::endl;
		//Save(file);
		//equiv++;
		ENTER_BLOCK();
		CheckGhostSharedCount(__FILE__, __LINE__);
		CheckCentroids(__FILE__, __LINE__);
		CheckOwners(__FILE__, __LINE__);
		CheckGIDs(__FILE__, __LINE__);
		CheckProcessors();
		EXIT_BLOCK();
		
		//TagInteger tag_new_owner = CreateTag("TEMPORARY_NEW_OWNER",DATA_INTEGER,CELL,NONE,1);
		//for(int k = 0; k < CellLastLocalID(); k++) if( isValidCell(k) )
		//{
		//	Cell c = CellByLocalID(k);
		//	if( !c.Hidden() && c.GetStatus() != Element::Owned )
		//	{
		//		assert( c.GetStatus() == Element::Shared || c.GetStatus() == Element::Ghost );
		//		tag_new_owner[c] = tag[c][0];
		//	}
		//}
		//ExchangeData(tag_new_owner,CELL,0);
		ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for(integer k = 0; k < CellLastLocalID(); k++) if( isValidCell(k) )
		{
			Cell c = CellByLocalID(k);
			if( !c.Hidden() && c.GetStatus() != Element::Owned )
			{
				integer new_owner = tag[c][0];
				//if( new_owner == -1 ) continue;
				if (new_owner == -1) continue;
#if !defined(NDEBUG)
				Storage::integer_array procs = c.IntegerArray(ProcessorsTag());
				if( !std::binary_search(procs.begin(),procs.end(),tag[c][0]) ) 
				{
					std::cout << "new owner " << tag[c][0] << " dist " << tag[c][1] << " is not among processors ";
					for(INMOST_DATA_ENUM_TYPE k = 0; k < procs.size();++k) std::cout << procs[k] << " ";
					std::cout << " for cell " << c.LocalID() << std::endl;
				}
#endif
				c.IntegerDF(tag_owner) = new_owner;
				if (GetProcessorRank() == new_owner)
					c.SetStatus(Element::Shared);
				else
					c.SetStatus(Element::Ghost);
			}
		}
		EXIT_BLOCK();
		//DeleteTag(tag_new_owner);
		//ExchangeData(tag_owner,CELL,0);
		ENTER_BLOCK();
		RecomputeParallelStorage(CELL);
		CheckGhostSharedCount(__FILE__,__LINE__,CELL);
		EXIT_BLOCK();
		DeleteTag(tag);
		
		
		
		ENTER_BLOCK();
		
		TagInteger new_owner = CreateTag("TEMPORARY_NEW_OWNER",DATA_INTEGER,NODE|EDGE|FACE,NONE,1);
		for(ElementType etype = FACE; etype >= NODE; etype = PrevElementType(etype))
		{
			REPORT_VAL("etype ", ElementTypeName(etype));
			ENTER_BLOCK();
			if (HideMarker())
			{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for (integer k = 0; k < LastLocalID(etype); k++) if (isValidElement(etype, k))
				{
					Element e = ElementByLocalID(etype, k);
					if (!e.Hidden() && e.GetStatus() != Element::Owned)
					{
						Element::adj_type& hc = HighConn(e.GetHandle());
						Storage::integer owner = INT_MAX;
						for (INMOST_DATA_ENUM_TYPE l = 0; l < hc.size(); ++l)
						{
							if (!GetMarker(hc[l], HideMarker()))
								owner = std::min(owner, IntegerDF(hc[l], tag_owner));
						}
						new_owner[e] = owner;
					}
				}
			}
			else
			{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for (integer k = 0; k < LastLocalID(etype); k++) if (isValidElement(etype, k))
				{
					Element e = ElementByLocalID(etype, k);
					if (e.GetStatus() != Element::Owned)
					{
						Element::adj_type& hc = HighConn(e.GetHandle());
						Storage::integer owner = INT_MAX;
						for (INMOST_DATA_ENUM_TYPE l = 0; l < hc.size(); ++l)
							owner = std::min(owner, IntegerDF(hc[l], tag_owner));
						new_owner[e] = owner;
					}
				}
			}
			EXIT_BLOCK();
			ReduceData(new_owner,etype,0,OperationMinOwner);
			ExchangeData(new_owner,etype);
			ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(integer k = 0; k < LastLocalID(etype); k++) if( isValidElement(etype,k) )
			{
				Element e = ElementByLocalID(etype,k);
				if( !e.Hidden() && e.GetStatus() != Element::Owned )
				{
					Storage::integer nowner = new_owner[e], & owner = e.IntegerDF(tag_owner);
					if( nowner != INT_MAX && owner != nowner )
					{
						owner = nowner;
						if (GetProcessorRank() == nowner)
							e.SetStatus(Element::Shared);
						else
							e.SetStatus(Element::Ghost);
					}
				}
			}
			EXIT_BLOCK();
			RecomputeParallelStorage(etype);
			ENTER_BLOCK();
			CheckGhostSharedCount(__FILE__,__LINE__,etype);
			EXIT_BLOCK();
		}
		DeleteTag(new_owner);
		//ComputeSharedProcs();
		//std::cout << "after_equiv"+std::to_string(counter)+".pvtk" << std::endl;
		//Save("after_equiv"+std::to_string(counter)+".pvtk");
		EXIT_BLOCK();
		//RecomputeParallelStorage(CELL | EDGE | FACE | NODE);
#endif
		EXIT_FUNC();
		counter++;
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
		//ENTER_FUNC();
		assert(OneType(type));
		//REPORT_VAL("test gid on ", ElementTypeName(type));
		bool ret = false;
		Tag test;
		if( GlobalIDTag().isValid() ) test = GlobalIDTag();
		else
		{
			//REPORT_STR("gid tag is not valid");
		}
		//else if( HaveTag("GLOBAL_ID") ) test = GetTag("GLOBAL_ID") );
		
		if( test.isValid() )
		{
			//bool ret = true;
			//for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype) ) if( etype & types )
			//	ret &= test.isDefined(etype);
			//return ret;
			ret = test.isDefined(type);
			//REPORT_STR("gid is valid and " << (ret ? "":"not") << " defined " << ElementTypeName(type));
		}
		//EXIT_FUNC();
		return ret;
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
					ExchangeData(t,mask,0);
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
					ExchangeData(t,mask,0);
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
					ExchangeData(t,mask,0);
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
		EXIT_FUNC();
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
			if( (it->GetStatus() == Element::Owned || it->GetStatus() == Element::Any) && (define_sparse || it->HaveData(num_tag)) )
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
			if( (GetStatus(*it) == Element::Owned || GetStatus(*it) == Element::Any) && (define_sparse || HaveData(*it,num_tag))) Integer(*it,num_tag) = shift++;
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
			if( (it->GetStatus() == Element::Owned || it->GetStatus() == Element::Any) && (define_sparse || it->HaveData(num_tag)) )
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
	
	Storage::enumerator Mesh::Integrate(Storage::enumerator input)
	{
		ENTER_FUNC();
		Storage::enumerator output = input;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,comm));
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
		std::vector<Storage::real> temp(size);
		if( !temp.empty() )
		{
			memcpy(&temp[0],input,sizeof(Storage::real)*size);
			REPORT_MPI(MPI_Allreduce(&temp[0],input,size,INMOST_MPI_DATA_REAL_TYPE,MPI_SUM,comm));
		}
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
		std::vector<Storage::integer> temp(size);
		if( !temp.empty() )
		{
			memcpy(&temp[0],input,sizeof(Storage::integer)*size);
			REPORT_MPI(MPI_Allreduce(&temp[0],input,size,INMOST_MPI_DATA_INTEGER_TYPE,MPI_SUM,comm));
		}
#else//USE_MPI
		(void) input;
		(void) size;
#endif//USE_MPI
		EXIT_FUNC();
	}
	
	void Mesh::Integrate(Storage::enumerator * input, Storage::integer size)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		std::vector<Storage::enumerator> temp(size);
		if( !temp.empty() )
		{
			memcpy(&temp[0],input,sizeof(Storage::enumerator)*size);
			REPORT_MPI(MPI_Allreduce(&temp[0],input,size,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,comm));
		}
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
	
	Storage::enumerator Mesh::ExclusiveSum(Storage::enumerator input)
	{
		ENTER_FUNC();
		Storage::enumerator output = 0;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Scan(&input,&output,1,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,comm));
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
	
	Storage::enumerator Mesh::AggregateMax(Storage::enumerator input)
	{
		ENTER_FUNC();
		Storage::enumerator output = input;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_ENUM_TYPE,MPI_MAX,comm));
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
		std::vector<Storage::real> temp(size);
		if( !temp.empty() )
		{
			memcpy(&temp[0],input,sizeof(Storage::real)*size);
			REPORT_MPI(MPI_Allreduce(&temp[0],input,size,INMOST_MPI_DATA_REAL_TYPE,MPI_MAX,comm));
		}
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
		std::vector<Storage::integer> temp(size);
		if( !temp.empty() )
		{
			memcpy(&temp[0],input,sizeof(Storage::integer)*size);
			REPORT_MPI(MPI_Allreduce(&temp[0],input,size,INMOST_MPI_DATA_INTEGER_TYPE,MPI_MAX,comm));
		}
#else//USE_MPI
		(void) input;
		(void) size;
#endif//USE_MPI
		EXIT_FUNC();
	}
	
	void Mesh::AggregateMax(Storage::enumerator * input, Storage::integer size)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		std::vector<Storage::enumerator> temp(size);
		if( !temp.empty() )
		{
			memcpy(&temp[0],input,sizeof(Storage::enumerator)*size);
			REPORT_MPI(MPI_Allreduce(&temp[0],input,size,INMOST_MPI_DATA_ENUM_TYPE,MPI_MAX,comm));
		}
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
	
	Storage::enumerator Mesh::AggregateMin(Storage::enumerator input)
	{
		ENTER_FUNC();
		Storage::enumerator output = input;
#if defined(USE_MPI)
		REPORT_MPI(MPI_Allreduce(&input,&output,1,INMOST_MPI_DATA_ENUM_TYPE,MPI_MIN,comm));
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
		std::vector<Storage::real> temp(size);
		if( !temp.empty() )
		{
			memcpy(&temp[0],input,sizeof(Storage::real)*size);
			REPORT_MPI(MPI_Allreduce(&temp[0],input,size,INMOST_MPI_DATA_REAL_TYPE,MPI_MIN,comm));
		}
#else//USE_MPI
		(void) input;
		(void) size;
#endif//USE_MPI
		EXIT_FUNC();
	}
	
	void Mesh::AggregateMin(Storage::enumerator * input, Storage::integer size)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		std::vector<Storage::enumerator> temp(size);
		if( !temp.empty() )
		{
			memcpy(&temp[0],input,sizeof(Storage::enumerator)*size);
			REPORT_MPI(MPI_Allreduce(&temp[0],input,size,INMOST_MPI_DATA_ENUM_TYPE,MPI_MIN,comm));
		}
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
		std::vector<Storage::integer> temp(size);
		if( !temp.empty() )
		{
			memcpy(&temp[0],input,sizeof(Storage::integer)*size);
			REPORT_MPI(MPI_Allreduce(&temp[0],input,size,INMOST_MPI_DATA_INTEGER_TYPE,MPI_MIN,comm));
		}
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
		tag_processors = CreateTag("PROCESSORS_LIST",DATA_INTEGER, ESET | NODE | EDGE | FACE | CELL,NONE);
		tag_layers = CreateTag("LAYERS",DATA_INTEGER,MESH,NONE,1);
		tag_bridge = CreateTag("BRIDGE",DATA_INTEGER,MESH,NONE,1);
		tag_sendto = CreateTag("PROTECTED_SENDTO",DATA_INTEGER, ESET | CELL | FACE | EDGE | NODE, ESET | CELL | FACE | EDGE | NODE);
		
		SyncDimensions();

#if defined(USE_MPI)
		randomizer = Random();
		
		//parallel_file_strategy = 1;
		parallel_file_strategy = 0;

		

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
		REPORT_MPI(err = MPI_Alloc_mem(GetProcessorsNumber()*sizeof(INMOST_DATA_BIG_ENUM_TYPE)*2,MPI_INFO_NULL,&shared_space));
		if( err )
		{
			int errclass;
			MPI_Error_class(err,&errclass);
			std::cout << "Cannot allocate shared space of size " << GetProcessorsNumber()*sizeof(INMOST_DATA_BIG_ENUM_TYPE)*2 << " reason is "; 
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
		REPORT_MPI(MPI_Win_create(shared_space,sizeof(INMOST_DATA_BIG_ENUM_TYPE)*GetProcessorsNumber()*2,sizeof(INMOST_DATA_BIG_ENUM_TYPE),MPI_INFO_NULL,comm,&window));
#endif //USE_MPI_P2P
#else //USE_MPI
		(void) _comm;
#endif //USE_MPI
		EXIT_FUNC();
	}
	
#if defined(USE_PARALLEL_WRITE_TIME)
	void Mesh::Enter() const { tab++; }
	void Mesh::Exit() const {tab--; }
	std::ostream & Mesh::WriteTab(std::ostream & f) const
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
	
	
	
	void determine_my_procs_low(Mesh * m, HandleType h, std::vector<Storage::integer> & result, std::vector<Storage::integer> & intersection)
	{
		MarkerType hm = m->HideMarker();
		Element::adj_type const & subelements = m->LowConn(h);
		Element::adj_type::const_iterator i = subelements.begin();
		result.clear();
		while(i != subelements.end())
		{
			if( hm && m->GetMarker(*i,hm) ) {i++; continue;}
			Storage::integer_array p = m->IntegerArrayDV(*i,m->ProcessorsTag());
			result.insert(result.end(),p.begin(),p.end());
			i++;
			break;
		}
		while(i != subelements.end())
		{
			if( hm && m->GetMarker(*i,hm) ) {i++; continue;}
			Storage::integer_array q = m->IntegerArrayDV(*i,m->ProcessorsTag());
			intersection.resize(std::max(static_cast<INMOST_DATA_ENUM_TYPE>(result.size()),static_cast<INMOST_DATA_ENUM_TYPE>(q.size())));
			std::vector<Storage::integer>::iterator qt = std::set_intersection(result.begin(),result.end(),q.begin(),q.end(),intersection.begin());
			intersection.resize(qt-intersection.begin());
			result.swap(intersection);
			i++;
		}
	}
	
	void determine_my_procs_high(Mesh * m, HandleType h, const Tag & procs, std::vector<Storage::integer> & result, std::vector<Storage::integer> & intersection)
	{
		MarkerType hm = m->HideMarker();
		Element::adj_type const & overelements = m->HighConn(h);
		if( overelements.empty() ) return;
		Element::adj_type::const_iterator i = overelements.begin();
		result.clear();
		while(i != overelements.end())
		{
			if( hm && m->GetMarker(*i,hm) ) {i++; continue;}
			Storage::integer_array q = m->IntegerArrayDV(*i,procs);
			intersection.resize(result.size()+q.size());
			std::vector<Storage::integer>::iterator qt = std::set_union(result.begin(),result.end(),q.begin(),q.end(),intersection.begin());
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
	public: bool operator () (const std::pair<Storage::integer,Storage::integer> & a, const std::pair<Storage::integer,Storage::integer> & b) {return a.first < b.first;}
	};
    

	void Mesh::SyncDimensions()
	{
#if defined(USE_MPI)
		MPI_Bcast(&dim, 1, INMOST_MPI_DATA_INTEGER_TYPE, 0, GetCommunicator());
#endif //USE_MPI
	}
	
	
	void Mesh::ResolveShared(bool only_new)
	{
		//only_new = false;
		ENTER_FUNC();
#if defined(USE_MPI)
		if( m_state == Mesh::Serial ) SetCommunicator(INMOST_MPI_COMM_WORLD);
		if( tag_global_id.isValid() ) tag_global_id = DeleteTag(tag_global_id,CELL | EDGE | FACE | NODE);
		integer dim = GetDimensions();
		int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber();
		REPORT_STR(__FUNCTION__ << " " << (only_new ? "only new" : "all"));
#if defined(USE_PARALLEL_STORAGE)
		shared_elements.clear();
		ghost_elements.clear();
#endif //USE_PARALLEL_STORAGE
		ReportParallelStorage();
		//determine which bboxes i intersect
		std::vector<int> procs;
		Storage::real bbox[6]; //local bounding box
		std::vector<Storage::real> bboxs((size_t)mpisize*6);
		//Compute local bounding box containing nodes.
		//Will be more convinient to compute (or store)
		//and communicate local octree over all the nodes.
		ENTER_BLOCK();
		for(integer k = 0; k < dim; k++)
		{
			bbox[k] = 1e20;
			bbox[k+dim] = -1e20;
		}
		ENTER_BLOCK();
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
		EXIT_BLOCK();
		for (integer k = 0; k < dim; k++)
		{
			bbox[k] -= GetEpsilon() * 2;
			bbox[k + dim] += GetEpsilon() * 2;
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
		ENTER_BLOCK();
		REPORT_MPI(MPI_Allgather(&bbox[0],dim*2,INMOST_MPI_DATA_REAL_TYPE,&bboxs[0],dim*2,INMOST_MPI_DATA_REAL_TYPE,comm));
		EXIT_BLOCK();
		// find all processors that i communicate with
		ENTER_BLOCK();
		for(int k = 0; k < mpisize; k++)
			if( k != mpirank )
			{
				bool flag = true;
				for(integer q = 0; q < dim; q++)
					flag &= !((bbox[q] > bboxs[(size_t)k*dim*2+q+dim]) || (bbox[dim+q] < bboxs[(size_t)k*dim*2+q]));
				if( flag ) procs.push_back(k);
			}
		EXIT_BLOCK();

#if !defined(NDEBUG)
		ENTER_BLOCK();
		{
			REPORT_STR("check exchange");
			Storage::real* mbox = &bboxs[mpirank * dim * 2];
			REPORT_STR("old bbox " << bbox[0] << ":" << bbox[3] << " " << bbox[1] << ":" << bbox[4] << " " << bbox[2] << ":" << bbox[5]);
			REPORT_STR("exh bbox " << mbox[0] << ":" << mbox[3] << " " << mbox[1] << ":" << mbox[4] << " " << mbox[2] << ":" << mbox[5]);
			REPORT_STR("dif bbox " << mbox[0]-bbox[0] << " " << mbox[3]-bbox[3] << " " << mbox[1]-bbox[1] << " " << mbox[4]-bbox[4] << " " << mbox[2]-bbox[2] << " " << mbox[5]-bbox[5]);
			std::vector<char> ploc(mpisize, 0), pglob(mpisize * mpisize, 0);
			for (std::vector<int>::iterator pt = procs.begin(); pt != procs.end(); ++pt)
				ploc[*pt] = 1; // fill my line
			REPORT_MPI(MPI_Allgather(&ploc[0], mpisize, MPI_CHAR, &pglob[0], mpisize, MPI_CHAR, comm));
			std::stringstream str;
			for (int i = 0; i < mpisize; ++i)
			{
				for (int j = 0; j < mpisize; ++j)
				{
					str << (pglob[i * mpisize + j] == 1 ? '1' : '0') << ' ';
					if (pglob[i * mpisize + j] != pglob[j * mpisize + i])
					{
						if (!mpirank) std::cout << __FILE__ << ":" << __LINE__ << ": no symmetry! " << i << " and " << j << std::endl;
						Storage::real* ibox = &bboxs[i * dim * 2];
						Storage::real* jbox = &bboxs[j * dim * 2];
						REPORT_STR("no symmetry i " << i << " j " << j);
						REPORT_STR("ith bbox " << ibox[0] << ":" << ibox[3] << " " << ibox[1] << ":" << ibox[4] << " " << ibox[2] << ":" << ibox[5]);
						REPORT_STR("jth bbox " << jbox[0] << ":" << jbox[3] << " " << jbox[1] << ":" << jbox[4] << " " << jbox[2] << ":" << jbox[5]);
						REPORT_STR("dif bbox " << ibox[0] - jbox[0] << " " << ibox[3] - jbox[3] << " " << ibox[1] - jbox[1] << " " << ibox[4] - jbox[4] << " " << ibox[2] - jbox[2] << " " << ibox[5] - jbox[5]);
						bool flag1 = true, flag2 = true;
						for (integer q = 0; q < dim; q++)
						{
							flag1 &= !((ibox[q] > jbox[q + dim]) || (ibox[dim + q] < jbox[q]));
							flag2 &= !((jbox[q] > ibox[q + dim]) || (jbox[dim + q] < ibox[q]));
						}
						REPORT_STR("check intersection: i<->j " << (flag1?"yes":"no") << " j<->i " << (flag2 ? "yes" : "no"));
					}
				}
				str << std::endl;
			}
			REPORT_STR(str.str());
		}
		EXIT_BLOCK();
#endif

		REPORT_VAL("neighbour processors",procs.size());
		EXIT_BLOCK();
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
		ReportParallelStorage();
		ENTER_BLOCK();
		{
			bool same_boxes = true, same_box;
			ENTER_BLOCK();
			for(int k = 0; k < mpisize && same_boxes; k++)
			{
				same_box = true;
				for(integer j = 0; j < dim*2; j++)
					same_box &= ::fabs(bbox[j] - bboxs[(size_t)k*dim*2+j]) < GetEpsilon();
				same_boxes &= same_box;
			}
			EXIT_BLOCK();
			
			//if( same_boxes )
			if( false )
			{
				REPORT_STR("All bounding boxes are the same - assuming that mesh is replicated over all nodes");
				//for(Mesh::iteratorElement it = BeginElement(CELL | EDGE | FACE | NODE); it != EndElement(); it++)
				ENTER_BLOCK();
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
				EXIT_BLOCK();
				//ComputeSharedProcs();
				RecomputeParallelStorage(CELL | EDGE | FACE | NODE);
				AssignGlobalID(CELL | EDGE | FACE | NODE);
			}
			else
			{
                /*
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
                */

				
				Storage::real epsilon = GetEpsilon();
				
				const bool precompute_geom = true;
				GeomParam table;	
				ENTER_BLOCK();
				if (precompute_geom)
				{
					for (ElementType etype = EDGE; etype <= CELL; etype = NextElementType(etype))
						if (!HaveGeometricData(CENTROID, etype))
							table[CENTROID] |= etype;
					PrepareGeometricData(table);

					REPORT_STR("Prepare geometric data");
				}
				EXIT_BLOCK();

			
				//for(iteratorNode it = BeginNode(); it != EndNode(); it++)
				ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(integer nit = 0; nit < NodeLastLocalID(); ++nit) if( isValidNode(nit) )
				{
					Node it = NodeByLocalID(nit);
					if (it->Hidden()) continue;
					if (only_new && !GetMarker(it->GetHandle(),NewMarker())) continue;
					Storage::integer_array arr = it->IntegerArrayDV(tag_processors);
					arr.resize(1);
					arr[0] = mpirank;
				}
				EXIT_BLOCK();
			

				buffer_type exch_data;
				std::vector< INMOST_DATA_REAL_TYPE > unpack_real;
				std::vector< INMOST_DATA_REAL_TYPE > pack_real;
				element_set sorted_nodes;
				
				sorted_nodes.reserve(NumberOfNodes());
				
				ENTER_BLOCK();
				//for(iteratorNode n = BeginNode(); n != EndNode(); n++)
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for (integer nit = 0; nit < NodeLastLocalID(); ++nit) if (isValidNode(nit))
				{
					Node n = NodeByLocalID(nit);
					if (n->Hidden()) continue;
					if (only_new && !n->GetMarker(NewMarker())) continue;
					real_array c = n->Coords();
					for(real_array::size_type k = 0; k < procs.size(); k++)
						if( point_in_bbox(c.data(),bboxs.data()+(size_t)procs[k]*dim*2,dim,GetEpsilon()) )
						{
#if defined(USE_OMP)
#pragma omp critical
#endif
							sorted_nodes.push_back(n->GetHandle());
							break;
						}
				}
				
				REPORT_STR("Prepare array of nodes");
				REPORT_VAL("share nodes", sorted_nodes.size());
				REPORT_VAL("total nodes", NumberOfNodes());
				EXIT_BLOCK();
				
				
				ENTER_BLOCK();
				{
					CentroidComparator cmp(this);
					if (!sorted_nodes.empty())
						std::sort(sorted_nodes.begin(), sorted_nodes.end(), cmp);
					//qsort(&sorted_nodes[0],sorted_nodes.size(),sizeof(Element *),CompareElementsCCentroid);
#if !defined(NDEBUG)
					bool have_gid = HaveGlobalID(NODE);
					integer bad = 0;
					if (!sorted_nodes.empty())
					{
						element_set::iterator next = sorted_nodes.begin(), cur = next++;
						while (next < sorted_nodes.end())
						{
							if (cmp.Compare(*next, *cur) == 0)
							{
								real_array c1 = Node(this, *cur).Coords();
								real_array c2 = Node(this, *next).Coords();
								real dist = sqrt(pow(c1[0] - c2[0], 2) + pow(c1[1] - c2[1], 2) + pow(c1[2] - c2[2], 2));
								std::cout << __FILE__ << ":" << __LINE__ << " same "
									<< c1[0] << " " << c1[1] << " " << c1[2] << (have_gid ? GlobalID(*cur) : -1)
									<< " and "
									<< c2[0] << " " << c2[1] << " " << c2[2] << (have_gid ? GlobalID(*cur) : -1)
									<< " dist " << dist << std::endl;
								REPORT_STR("same "
									<< c1[0] << " " << c1[1] << " " << c1[2] << (have_gid ? GlobalID(*cur) : -1)
									<< " and "
									<< c2[0] << " " << c2[1] << " " << c2[2] << (have_gid ? GlobalID(*cur) : -1)
									<< " dist " << dist);
								bad++;
							}
							cur = next++;
						}
						if (bad) std::cout << __FILE__ << ":" << __LINE__ << " " << GetProcessorRank() << " duplicate coords " << bad << std::endl;
					}
#endif
				}
				REPORT_STR("Sort nodes");
				EXIT_BLOCK();
				
				ENTER_BLOCK();
				pack_real.reserve(sorted_nodes.size()*dim);
				EXIT_BLOCK();
				
				ENTER_BLOCK();
				for(element_set::iterator it = sorted_nodes.begin(); it != sorted_nodes.end(); it++)
				{
					Storage::real_array arr = RealArrayDF(*it,CoordsTag());
					pack_real.insert(pack_real.end(),arr.begin(),arr.end());
				}
				REPORT_STR("Gather coordinates");
				EXIT_BLOCK();
				
				
				//~ int position = 0;

				
				ENTER_BLOCK();
				pack_data_vector(exch_data,pack_real,GetCommunicator());
				EXIT_BLOCK();
				
				
				REPORT_STR("Pack coordinates");
				

				
				ENTER_BLOCK();
				{
					exch_reqs_type send_reqs;
					exch_recv_reqs_type recv_reqs;
					exch_buffer_type send_buffs(procs.size()), recv_buffs(procs.size());
					std::vector<int> done;

					std::vector<size_t> sendsizeall(mpisize);
					ENTER_BLOCK();
					
					for(size_t k = 0; k < procs.size(); k++)
					{
						send_buffs[k].first = procs[k];
						pack_data(send_buffs[k].second,exch_data.size(),GetCommunicator());
						recv_buffs[k].first = procs[k];
						recv_buffs[k].second.resize(send_buffs[k].second.size());
					}
					EXIT_BLOCK();
					ExchangeBuffersInner(send_buffs,recv_buffs,send_reqs,recv_reqs);
					ENTER_BLOCK();
					while( !(done = FinishRequests(recv_reqs)).empty() )
					{
						for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
						{
							//~ position = 0;
							//~ MPI_Unpack(&recv_buffs[*qt].second[0],static_cast<INMOST_MPI_SIZE>(recv_buffs[*qt].second.size()),&position,&sendsizeall[procs[*qt]*2],2,MPI_UNSIGNED,comm);
							ENTER_BLOCK();
							REPORT_VAL("processor",procs[*qt]);
							size_t buf_pos = 0;
							unpack_data(recv_buffs[*qt].second,buf_pos,sendsizeall[procs[*qt]],GetCommunicator());
							EXIT_BLOCK();
						}
					}
					EXIT_BLOCK();
					ENTER_BLOCK();
					if( !send_reqs.empty() )
					{
						REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(send_reqs.size()),&send_reqs[0],MPI_STATUSES_IGNORE));
					}
					EXIT_BLOCK();

					ENTER_BLOCK();
					{
						
						ENTER_BLOCK();
						for(size_t k = 0; k < procs.size(); k++)
						{
							send_buffs[k].first = procs[k];
							send_buffs[k].second = exch_data;
							recv_buffs[k].first = procs[k];
							recv_buffs[k].second.resize(sendsizeall[procs[k]]);
						}
						EXIT_BLOCK();
						//PrepareReceiveInner(send_buffs,recv_buffs);
						ExchangeBuffersInner(send_buffs,recv_buffs,send_reqs,recv_reqs);
						
						
						ENTER_BLOCK();
						while( !(done = FinishRequests(recv_reqs)).empty() )
						{
							for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
							{
								REPORT_STR("receive node coordinates");
								REPORT_VAL("processor",recv_buffs[*qt].first);
								size_t count = 0;
								//~ int position = 0;
								//~ unpack_real.resize(sendsizeall[recv_buffs[*qt].first*2+1]);
								
								//TODO: overflow
								//~ MPI_Unpack(&recv_buffs[*qt].second[0],static_cast<INMOST_MPI_SIZE>(recv_buffs[*qt].second.size()),&position,&unpack_real[0],static_cast<INMOST_MPI_SIZE>(unpack_real.size()),INMOST_MPI_DATA_REAL_TYPE,comm);
								size_t buf_pos = 0;
								unpack_data_vector(recv_buffs[*qt].second,buf_pos,unpack_real,GetCommunicator());
								std::vector<Storage::real>::iterator it1 = pack_real.begin() , it2 = unpack_real.begin();
								CentroidComparator cmp(this);
								while(it1 != pack_real.end() && it2 != unpack_real.end() )
								{
									int res = cmp.Compare(&*it1,&*it2);
									/*
									for(integer k = 0; k < dim; k++)
										if( ::fabs((*(it1+k))-(*(it2+k))) > epsilon )
										{
											if( (*(it1+k)) < (*(it2+k)) ) res = -1;
											else res = 1;
											break;
										}
									*/
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
								REPORT_STR("Intersect coordinates");
							}
						}
						EXIT_BLOCK();
						ENTER_BLOCK();
						if( !send_reqs.empty() )
						{
							REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(send_reqs.size()),&send_reqs[0],MPI_STATUSES_IGNORE));
						}
						EXIT_BLOCK();
					}
					REPORT_STR("Intersect all coordinates");
					EXIT_BLOCK();
					
					
					Element::Status estat;
					ENTER_BLOCK();
					for(Mesh::iteratorElement it = BeginElement(NODE); it != EndElement(); it++)
					{
						if (only_new && it->GetMarker(NewMarker()) == false) continue;
						integer owner;
						Storage::integer_array v = it->IntegerArrayDV(tag_processors);
						std::sort(v.begin(),v.end());
						if( v.empty() )
						{
							owner = mpirank;
							v.push_back(mpirank);
						}
						else
							owner = std::min<Storage::integer>(mpirank,v[0]);
						
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
						it->SetStatus(estat);
					}
					EXIT_BLOCK();
					//ComputeSharedProcs();
					RecomputeParallelStorage(NODE);
					CheckGhostSharedCount(__FILE__, __LINE__, NODE);
					AssignGlobalID(NODE);
					ReportParallelStorage();
					REPORT_STR("Set parallel info for nodes");
					
					ENTER_BLOCK();
					CheckGhostSharedCount(__FILE__, __LINE__,NODE);
					CheckOwners(__FILE__, __LINE__);
					CheckGIDs(__FILE__, __LINE__,NODE);
					CheckProcessors();
					CheckCentroids(__FILE__, __LINE__);
					EXIT_BLOCK();
					
					MarkerType new_lc = 0;
					if( only_new ) new_lc = CreateMarker();
					//TagInteger lc_mrk = CreateTag("LC_MARKER",DATA_INTEGER,CELL|FACE|EDGE,NONE,1);
					MarkerType hm = HideMarker();
					MarkerType nm = NewMarker();
					
					//for (iteratorElement it = BeginElement(CELL|FACE|EDGE); it != EndElement(); it++) lc_mrk[*it] = 0;
					ENTER_BLOCK();
					for(ElementType current_mask = EDGE; current_mask <= CELL; current_mask = NextElementType(current_mask))
					{
						REPORT_STR("Set parallel info for");
						REPORT_VAL("type",ElementTypeName(current_mask));
						
						
						//int owned_elems = 0;
						//int shared_elems = 0;
						integer owner;
						Element::Status estat;
						
						ENTER_BLOCK();
						if(only_new)
						{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
							for(integer eit = 0; eit < LastLocalID(current_mask); ++eit) if( isValidElement(current_mask,eit) )
							{
								Element it = ElementByLocalID(current_mask,eit);
								if( it->GetMarker(hm) ) continue;
								if( it->GetMarker(nm) ) 
								{
									it->SetMarker(new_lc);
									//lc_mrk[ComposeHandle(current_mask,eit)] = 1;
									continue;
								}
								Element::adj_type low_conn = LowConn(ComposeHandle(current_mask,eit));
								for (Element::adj_type::iterator jt = low_conn.begin(); jt != low_conn.end(); jt++)
									if( GetMarker(*jt,nm)) 
									{
										it->SetMarker(new_lc);
										//lc_mrk[ComposeHandle(current_mask,eit)] = 1;
										break;
									}
							}
						}
						EXIT_BLOCK();
						//stringstream ss;
						//ss << "file_" << GetProcessorRank() << ".txt";
						//				ofstream ofs(ss.str().c_str());
						//Determine what processors potentially share the element
						ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel
#endif
						{
							std::vector<Storage::integer> result, intersection;
#if defined(USE_OMP)
#pragma omp for
#endif
							for(integer eit = 0; eit < LastLocalID(current_mask); ++eit)
							{
								if( isValidElement(current_mask,eit) )
								{
									Element it = ElementByLocalID(current_mask,eit);
									if( it->GetMarker(hm) ) continue;
									if (only_new && !it->GetMarker(new_lc)) continue;
									determine_my_procs_low(this,it->GetHandle(), result, intersection);
									Storage::integer_array p = it->IntegerArrayDV(tag_processors);
									if( result.empty() )
									{
										p.clear();
										p.push_back(mpirank);
									}
									else p.replace(p.begin(),p.end(),result.begin(),result.end());
								}
							}
						}
						
						//REPORT_VAL("predicted owned elements",owned_elems);
						//REPORT_VAL("predicted shared elements",shared_elems);
						REPORT_STR("Predict processors for elements");
						EXIT_BLOCK();
						
						
						//Initialize mapping that helps get local id by global id
						std::vector<std::pair<integer,integer> > mapping;
						ENTER_BLOCK();
						REPORT_VAL("mapping type",ElementTypeName(PrevElementType(current_mask)));
						//for(Mesh::iteratorElement it = BeginElement(PrevElementType(current_mask)); it != EndElement(); it++)
						for(integer i = 0; i < LastLocalID(PrevElementType(current_mask)); ++i) if( isValidElement(PrevElementType(current_mask),i))
						{
							Element it = ElementByLocalID(PrevElementType(current_mask), i);
							if (it->Hidden()) continue;
							mapping.push_back(std::make_pair(it->GlobalID(),it->LocalID()));
						}
						if( !mapping.empty() ) 
							std::sort(mapping.begin(),mapping.end(),MappingComparator());
						
						REPORT_VAL("mapping size",mapping.size())
						REPORT_STR("Compute global to local indexes mapping");
						EXIT_BLOCK();
						//Initialize arrays
						
						
						std::vector<int> procs = ComputeSharedProcs(PrevElementType(current_mask));
						std::vector< std::vector<integer> > message_recv(procs.size());
						std::vector< element_set > elements(procs.size());
						
						exch_reqs_type send_reqs;
						exch_recv_reqs_type recv_reqs;
						exch_buffer_type send_buffs(procs.size()), recv_buffs(procs.size());

						ENTER_BLOCK();
						//Gather all possible shared elements and send global ids of their connectivity to the neighbouring proccessors
//#if defined(USE_OMP)
//#pragma omp parallel
//#endif
						{
							std::vector<integer> message_send;
//#if defined(USE_OMP)
//#pragma omp for schedule(dynamic,1)
//#endif
							for (std::vector<int>::iterator p = procs.begin(); p != procs.end(); p++)
							//for (int m = 0; m < (int)procs.size(); ++m)
							{
								//std::vector<int>::iterator p = procs.begin() + m;
								int m = static_cast<int>(p - procs.begin());
								REPORT_VAL("for processor", *p);
								{
									message_send.clear();
									message_send.push_back(0);
									ENTER_BLOCK();
									//for(Mesh::iteratorElement it = BeginElement(current_mask); it != EndElement(); it++)
									for (int i = 0; i < LastLocalID(current_mask); ++i) if (isValidElement(current_mask, i))
									{
										Element it = ElementByLocalID(current_mask, i);
										if (it->Hidden()) continue;
										if (only_new && !it->GetMarker(new_lc)) continue;
										Storage::integer_array pr = it->IntegerArrayDV(tag_processors);
										if (std::binary_search(pr.begin(), pr.end(), *p))
										{
											Element::adj_type& sub = LowConn(it->GetHandle());
											if (sub.size() == 0) throw Impossible;
											integer message_size_pos = (integer)message_send.size();
											message_send.push_back(0);
											REPORT_VAL("number of connections",sub.size());
											REPORT_STR("element " << ElementTypeName(current_mask) << ":" << it->LocalID());
											for (Element::adj_type::iterator kt = sub.begin(); kt != sub.end(); kt++) if (!hm || !GetMarker(*kt, hm))
											{
												message_send.push_back(GlobalID(*kt));
												message_send[message_size_pos]++;
#if defined(USE_PARALLEL_WRITE_TIME)
												INMOST_DATA_REAL_TYPE cnt[3] = { 0,0,0 };
												ElementByLocalID(PrevElementType(current_mask), GetHandleID(*kt))->Centroid(cnt);
												REPORT_STR("global id " << GlobalID(*kt) << " local id " << GetHandleID(*kt) << " " << Element::StatusName(Element(this,*kt)->GetStatus()) << " cnt " << cnt[0] << " " << cnt[1] << " " << cnt[2]);
#endif
											}
											//~ message_send[1]++;
											message_send[0]++;
											elements[m].push_back(it->GetHandle());
										}
									}
									EXIT_BLOCK();
									REPORT_VAL("gathered elements", elements[m].size());

									//ENTER_BLOCK();
									send_buffs[m].first = *p;
									pack_data_vector(send_buffs[m].second, message_send, GetCommunicator());
									recv_buffs[m].first = *p;
									//EXIT_BLOCK();
								}
							}
						}
						
						
						REPORT_STR("Pack messages for other processors");
						EXIT_BLOCK();
						
						PrepareReceiveInner(UnknownSize,send_buffs,recv_buffs);
						ExchangeBuffersInner(send_buffs,recv_buffs,send_reqs,recv_reqs);
						
						
						std::vector<int> done;
						ENTER_BLOCK();
						while( !(done = FinishRequests(recv_reqs)).empty() )
						{
							for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
							{
								//~ int position = 0;
								//~ int size;
								int pos = -1;
								for(std::vector<int>::iterator p = procs.begin(); p != procs.end(); p++)
									if( *p == recv_buffs[*qt].first )
									{
										pos = static_cast<int>(p - procs.begin());
										break;
									}
								if( pos == -1 ) throw Impossible;
								size_t buf_pos = 0;
								unpack_data_vector(recv_buffs[*qt].second,buf_pos,message_recv[pos],GetCommunicator());
							}
						}
						
						
						REPORT_STR("Exchange messages");
						EXIT_BLOCK();
						
						ENTER_BLOCK();
						//Now find the difference of local elements with given processor number and remote elements
						//for(p = procs.begin(); p != procs.end(); p++)
//#if defined(USE_OMP)
//#pragma omp parallel for schedule(dynamic,1)
//#endif
//						for(int m = 0; m < (int)procs.size(); ++m)
						for (std::vector<int>::iterator p = procs.begin(); p != procs.end(); p++)
						{
							//std::vector<int>::iterator p = procs.begin() + m;
							int m = static_cast<int>(p - procs.begin());
							//REPORT_VAL("on processor",*p);
							element_set remote_elements;
							int pos = 0;
							if( message_recv[m].empty() ) continue;
							integer num_remote_elements = message_recv[m][pos++];
							REPORT_VAL("number of remote elements",num_remote_elements);
							ENTER_BLOCK();
							for(integer i = 0; i < num_remote_elements; i++)
							{
								integer conn_size = message_recv[m][pos++], flag = 1;
								REPORT_VAL("number of connxections",conn_size);
								std::vector<HandleType> sub_elements;
								for(integer j = 0; j < conn_size; j++)
								{
									integer global_id = message_recv[m][pos++];
									integer find = -1;
									std::vector<std::pair<integer,integer> >::iterator it = std::lower_bound(mapping.begin(),mapping.end(),std::make_pair(global_id,0),MappingComparator());
									if( it != mapping.end() && it->first == global_id) 
										find = static_cast<integer>(it-mapping.begin());
									if( find == -1 ) 
									{
										REPORT_STR("global id " << global_id << " local id -1");
										flag = 0;
										pos += conn_size-j-1;
										break;
									}
									integer find_local_id = mapping[find].second;
#if defined(USE_PARALLEL_WRITE_TIME)
									real cnt[3] = { 0,0,0 };
									ElementByLocalID(PrevElementType(current_mask), find_local_id).Centroid(cnt);
									REPORT_STR("global id " << global_id 
										<< " local id " << find_local_id 
										<< " " << Element::StatusName(ElementByLocalID(PrevElementType(current_mask), find_local_id)->GetStatus())
										<< " " << cnt[0] << " " << cnt[1] << " " << cnt[2]);
#endif
									if( GetMarker(ComposeHandle(PrevElementType(current_mask), find_local_id),hm) ) std::cout << "Found hidden element" << std::endl;
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
							EXIT_BLOCK();
							if( !remote_elements.empty() )
							{
								REPORT_VAL("first",remote_elements.front());
								REPORT_VAL("first type",ElementTypeName(GetHandleElementType(remote_elements.front())));
								REPORT_VAL("last",remote_elements.back());
								REPORT_VAL("last type",ElementTypeName(GetHandleElementType(remote_elements.back())));
							}
							ENTER_BLOCK();
							std::sort(remote_elements.begin(),remote_elements.end());
							EXIT_BLOCK();
							REPORT_VAL("original elements size",elements[m].size());
							if( !elements[m].empty() )
							{
								REPORT_VAL("first",elements[m].front());
								REPORT_VAL("first type",ElementTypeName(GetHandleElementType(elements[m].front())));
								REPORT_VAL("last",elements[m].back());
								REPORT_VAL("last type",ElementTypeName(GetHandleElementType(elements[m].back())));
							}
							ENTER_BLOCK();
							std::sort(elements[m].begin(),elements[m].end());
							EXIT_BLOCK();
							element_set result;
							element_set::iterator set_end;
							ENTER_BLOCK();
							result.resize(elements[m].size());
							set_end = std::set_difference(elements[m].begin(),elements[m].end(),remote_elements.begin(),remote_elements.end(), result.begin());
							result.resize(set_end-result.begin());
							EXIT_BLOCK();
							REPORT_VAL("set difference size",result.size());
							ENTER_BLOCK();
							//elements in result are wrongly marked as ghost
#if defined(USE_OMP)
#pragma omp critical
#endif
							for(element_set::iterator qt = result.begin(); qt != result.end(); qt++)
							{
								integer_array pr = IntegerArrayDV(*qt,tag_processors);
								integer_array::iterator find = std::lower_bound(pr.begin(),pr.end(),*p);
								pr.erase(find);
							}
							EXIT_BLOCK();
						}
						
						
						REPORT_STR("Determine true shared based on remote info");
						EXIT_BLOCK();
						
						ENTER_BLOCK();
						//Now mark all the processors status
						//for(Mesh::iteratorElement it = BeginElement(current_mask); it != EndElement(); it++)
						for(int i = 0; i < LastLocalID(current_mask); ++i) if( isValidElement(current_mask,i))
						{
							Element it = ElementByLocalID(current_mask, i);
							if (only_new && !it->GetMarker(new_lc)) continue;
							Storage::integer_array pr = it->IntegerArrayDV(tag_processors);
							if( pr.empty() )
							{
								owner = mpirank;
								pr.push_back(mpirank);
							}
							else
								owner = std::min<Storage::integer>(mpirank,pr[0]);
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
							SetStatus(it->GetHandle(), estat);
						}
						EXIT_BLOCK();
						RecomputeParallelStorage(current_mask);
						ENTER_BLOCK();

						ENTER_BLOCK();
						CheckGhostSharedCount(__FILE__, __LINE__,current_mask);
						CheckOwners(__FILE__, __LINE__);
						CheckGIDs(__FILE__, __LINE__,current_mask);
						CheckProcessors();
						CheckCentroids(__FILE__, __LINE__);
						EXIT_BLOCK();

						if( !send_reqs.empty() )
						{
							REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(send_reqs.size()),&send_reqs[0],MPI_STATUSES_IGNORE));
						}
						EXIT_BLOCK();
						AssignGlobalID(current_mask);
						ReportParallelStorage();

						ENTER_BLOCK();
						CheckGhostSharedCount(__FILE__, __LINE__,current_mask);
						CheckOwners(__FILE__, __LINE__);
						CheckGIDs(__FILE__, __LINE__,current_mask);
						CheckProcessors();
						CheckCentroids(__FILE__, __LINE__);
						EXIT_BLOCK();
					}
					EXIT_BLOCK();
					if( only_new ) ReleaseMarker(new_lc,EDGE|FACE|CELL);
					
				}
				EXIT_BLOCK();
				//if( !only_new ) 
				if (precompute_geom)
				{
					ENTER_BLOCK();
					RemoveGeometricData(table);
					EXIT_BLOCK();
				}
			}
		}
		EXIT_BLOCK();
		ENTER_BLOCK();
		ReportParallelStorage();
		
		ENTER_BLOCK();
		CheckGhostSharedCount(__FILE__, __LINE__);
		CheckCentroids(__FILE__, __LINE__);
		CheckOwners(__FILE__, __LINE__);
		CheckGIDs(__FILE__, __LINE__);
		CheckProcessors();
		EXIT_BLOCK();


		EquilibrateGhost();

		ENTER_BLOCK();
		CheckGhostSharedCount(__FILE__, __LINE__);
		CheckCentroids(__FILE__, __LINE__);
		CheckOwners(__FILE__, __LINE__);
		CheckGIDs(__FILE__, __LINE__);
		CheckProcessors();
		EXIT_BLOCK();

		CheckProcsSorted(__FILE__,__LINE__);
		EXIT_BLOCK();
		//RemoveGhost();
#else //USE_MPI
		AssignGlobalID(CELL | FACE | EDGE | NODE);
        (void)only_new;
#endif //USE_MPI
		EXIT_FUNC();
	}
	
	
	void Mesh::RemoveGhost(MarkerType marker)
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
            if (marker &&  !it->GetMarker(marker)) continue;
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

		if( !delete_elements.empty())
		{
			MarkerType del = CreateMarker();
			SetMarkerArray(&delete_elements[0], (Storage::enumerator)delete_elements.size(), del);
			RemoveLinksToDeletedElements(del);
			RemMarkerArray(&delete_elements[0], (Storage::enumerator)delete_elements.size(), del);
			ReleaseMarker(del);
		}

		//for(element_set::iterator it = delete_elements.begin(); it != delete_elements.end(); it++) Destroy(*it);
		for(element_set::iterator it = delete_elements.begin(); it != delete_elements.end(); it++) Delete(*it);
		
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
				if( HaveGlobalID(mask) )
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
				if( HaveGlobalID(mask) )
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),GlobalIDComparator(this));
				else
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),CentroidComparator(this));
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_shared.clear();
#endif //USE_PARALLEL_STORAGE

			if (!delete_elements.empty())
			{
				MarkerType del = CreateMarker();
				SetMarkerArray(&delete_elements[0], (Storage::enumerator)delete_elements.size(), del);
				RemoveLinksToDeletedElements(del);
				RemMarkerArray(&delete_elements[0], (Storage::enumerator)delete_elements.size(), del);
				ReleaseMarker(del);
			}

			//for(element_set::iterator it = delete_elements.begin(); it != delete_elements.end(); it++) Destroy(*it);
			for(element_set::iterator it = delete_elements.begin(); it != delete_elements.end(); it++) Delete(*it);
			delete_elements.clear();
		}		
		DeleteTag(tag_delete);
		//ComputeSharedProcs();
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
		Tag tag_delete = CreateTag("TEMPORARY_DELETE_GHOST_ELEMENTS_TAG",DATA_INTEGER,ESET | CELL | FACE | EDGE | NODE, ESET | CELL | FACE | EDGE | NODE);
		element_set delete_elements;
#if defined(USE_PARALLEL_STORAGE)
		proc_elements del_shared, del_ghost;
#endif //USE_PARALLEL_STORAGE
		
		for(ElementType mask = ESET; mask >= NODE; mask = PrevElementType(mask))
		{
			INMOST_DATA_ENUM_TYPE cnt = 0;
			ENTER_BLOCK();
			for(const HandleType * it = ghost; it != ghost + num; it++)
				if( (GetHandleElementType(*it) & mask) && GetStatus(*it) == Element::Ghost ) 
				{
					Integer(*it,tag_delete) = mpirank;
					cnt++;
				}
			EXIT_BLOCK();
			if( mask & (ESET|CELL) )
			{
				cnt = Integrate(cnt);
				if( !cnt ) continue;
			}
			ENTER_BLOCK();
			if( mask & (FACE|EDGE|NODE) )
			{
				for(integer i = 0; i < LastLocalID(mask); ++i) if( isValidElement(mask,i) )
				{
					Element it = ElementByLocalID(mask, i);
					if (it->Hidden()) continue;
					Element::Status estat = it->GetStatus();
					if( estat == Element::Owned || estat == Element::Shared ) continue;
					if( it->nbAdjElements(NextElementType(mask)) == 0 ) it->Integer(tag_delete) = mpirank;
				}
			}
			EXIT_BLOCK();
			ReduceData(tag_delete,mask,0,DeleteUnpack);
			ExchangeData(tag_delete,mask,0);
			ENTER_BLOCK();
			for (integer i = 0; i < LastLocalID(mask); ++i) if (isValidElement(mask,i))
			{
				Element it = ElementByLocalID(mask, i);
				if (it->Hidden()) continue;
				Element::Status estat = it->GetStatus();
				if( estat == Element::Owned ) continue;
				if( it->HaveData(tag_delete) )
				{
					Storage::integer_array del_procs = it->IntegerArray(tag_delete);
					std::sort(del_procs.begin(),del_procs.end());
					
					if( estat == Element::Ghost && std::binary_search(del_procs.begin(),del_procs.end(),mpirank) )
					{
#if defined(USE_PARALLEL_STORAGE)
						del_ghost[it->IntegerDF(tag_owner)].push_back(it->GetHandle());
#endif //USE_PARALLEL_STORAGE
						delete_elements.push_back(it->GetHandle());
					}
					else
					{
						Storage::integer_array procs = it->IntegerArrayDV(tag_processors);
						std::vector<Storage::integer> result(procs.size());
#if defined(USE_PARALLEL_STORAGE)
						if( estat == Element::Shared )
						{
							for(Storage::integer_array::iterator vit = del_procs.begin(); vit != del_procs.end(); vit++)
								del_shared[*vit].push_back(it->GetHandle());
						}
#endif	//USE_PARALLEL_STORAGE
						std::vector<Storage::integer>::iterator end = std::set_difference(procs.begin(),procs.end(),del_procs.begin(),del_procs.end(),result.begin());
						result.resize(end-result.begin());
						procs.clear();
						procs.insert(procs.begin(),result.begin(),result.end());
						
						if( procs.size() == 1 && procs[0] == mpirank )
							it->SetStatus(Element::Owned);
					}
				}
			}
			EXIT_BLOCK();
			ENTER_BLOCK();
#if defined(USE_PARALLEL_STORAGE)
			for(proc_elements::iterator it = del_ghost.begin(); it != del_ghost.end(); it++)
			{
				//std::cout << GetProcessorRank() << " ghost delete size " << it->second.size() << std::endl;
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
				if( HaveGlobalID(mask) )
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),GlobalIDComparator(this));
				else
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),CentroidComparator(this));
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_ghost.clear();
			for(proc_elements::iterator it = del_shared.begin(); it != del_shared.end(); it++)
			{
				//std::cout << GetProcessorRank() << " shared delete size " << it->second.size() << std::endl;
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
				if( HaveGlobalID(mask) )
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),GlobalIDComparator(this));
				else
					end = std::set_difference(ref.begin(),ref.end(),it->second.begin(),it->second.end(),result.begin(),CentroidComparator(this));
				result.resize(end-result.begin());
				ref.swap(result);
			}
			del_shared.clear();
#endif //USE_PARALLEL_STORAGE
			EXIT_BLOCK();
			if (!delete_elements.empty())
			{
				ENTER_BLOCK();
				MarkerType del = CreateMarker();
				SetMarkerArray(&delete_elements[0], (Storage::enumerator)delete_elements.size(), del);
				RemoveLinksToDeletedElements(del);
				RemMarkerArray(&delete_elements[0], (Storage::enumerator)delete_elements.size(), del);
				ReleaseMarker(del);
				EXIT_BLOCK();
			}
			ENTER_BLOCK();
			//for(element_set::iterator it = delete_elements.begin(); it != delete_elements.end(); it++) Destroy(*it);
			for(element_set::iterator it = delete_elements.begin(); it != delete_elements.end(); it++) Delete(*it);
			EXIT_BLOCK();
			delete_elements.clear();
			//std::cout << GetProcessorRank() << " finish " << ElementTypeName(mask) << std::endl;
		}
		DeleteTag(tag_delete);
		//ComputeSharedProcs();
#else //USE_MPI
		(void) ghost;
		(void) num;
#endif //USE_MPI
		EXIT_FUNC();
		ReportParallelStorage();
	}

	
	void Mesh::AssignGlobalID(ElementType mask)
	{
		ENTER_FUNC();
		ENTER_BLOCK();
		tag_global_id = CreateTag("GLOBAL_ID", DATA_INTEGER, mask, NONE, 1);
		EXIT_BLOCK();
#if defined(USE_MPI)
		if( m_state == Parallel )
		{
			INMOST_DATA_ENUM_TYPE number, shift[5] = { 0,0,0,0,0 }, shift_recv[5] = { 0,0,0,0,0 };
			ENTER_BLOCK();
			//for(ElementType currenttype = NODE; currenttype <= ESET; currenttype = NextElementType(currenttype) )
			//num = ElementNum(currenttype);
			for (int num = ElementNum(NODE); num <= ElementNum(ESET); ++num) if (mask & ElementTypeFromDim(num))
			{
				number = 0;
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:number)
#endif
				for (integer i = 0; i < LastLocalIDNum(num); ++i) if (isValidElementNum(num, i))
				{
					Element it = ElementByLocalIDNum(num, i);
					if( it->GetStatus() != Element::Ghost )
						number++;
				}
				shift[num] = number;
				REPORT_STR("enumerate on " << ElementTypeName(ElementTypeFromDim(num)) << " number " << number);
			}
			EXIT_BLOCK();
			ENTER_BLOCK();
			{
				int ierr;
				REPORT_MPI(ierr = MPI_Scan(shift,shift_recv,5,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,GetCommunicator()));
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
			}
			EXIT_BLOCK();
			ENTER_BLOCK();
			//for(ElementType currenttype = NODE; currenttype <= ESET; currenttype = NextElementType(currenttype) )
			//num = ElementNum(currenttype);
//#if defined(USE_OMP)
//#pragma omp parallel for schedule(dynamic,1)
//#endif
			for (int num = ElementNum(NODE); num <= ElementNum(ESET); ++num) if( mask & ElementTypeFromDim(num) )
			{
				INMOST_DATA_ENUM_TYPE local_shift = shift_recv[num] - shift[num];
				for(integer i = 0; i < LastLocalIDNum(num); ++i) if( isValidElementNum(num,i))
				{
					Element it = ElementByLocalIDNum(num, i);
					if( it->GetStatus() != Element::Ghost )
						it->IntegerDF(tag_global_id) = local_shift++;
				}
				REPORT_STR("shift on " << ElementTypeName(ElementTypeFromDim(num)) << " number " << local_shift);
			}
			EXIT_BLOCK();
			ENTER_BLOCK();
			if( tag_global_id.isValid() )
			{
				Tag exchange = tag_global_id;
				//Protect ExchangeData from using global id to sort ghost/shared elements
				tag_global_id = Tag();
				ExchangeData(exchange,mask,0);
				tag_global_id = exchange;
			}
			EXIT_BLOCK();
			ENTER_BLOCK();
			CheckGIDs(__FILE__, __LINE__,mask);
			SortParallelStorage(mask);
			CheckGIDs(__FILE__, __LINE__,mask);
			EXIT_BLOCK();
		}
		else
		{
#endif //USE_MPI
#if defined(USE_OMP)
#pragma omp parallel for schedule(dynamic,1)
#endif
			for (int num = ElementNum(NODE); num <= ElementNum(ESET); ++num) if (mask & ElementTypeFromDim(num))
			{
				INMOST_DATA_ENUM_TYPE number = 0;
				for (integer i = 0; i < LastLocalIDNum(num); ++i) if (isValidElementNum(num, i))
					ElementByLocalIDNum(num, i)->IntegerDF(tag_global_id) = number++;
			}
#if defined(USE_MPI)
		}
#endif //USE_MPI
		EXIT_FUNC();
	}
		
	std::vector<int> Mesh::ComputeSharedProcs(const parallel_storage & from, const parallel_storage & to)
	{
		std::vector<int> ret;
		ENTER_FUNC();
#if defined(USE_MPI)
		std::set<int> shared_procs;
		int mpirank = GetProcessorRank();
		for(parallel_storage::const_iterator it = from.begin(); it != from.end(); ++it) 
			shared_procs.insert(it->first);
		for(parallel_storage::const_iterator it = to.begin(); it != to.end(); ++it) 
			shared_procs.insert(it->first);
		if( !shared_procs.empty() )
		{
			std::set<int>::iterator ir = shared_procs.find(mpirank);
			if( ir != shared_procs.end() ) shared_procs.erase(ir);
			if( !shared_procs.empty() )
			{
				ret.insert(ret.end(),shared_procs.begin(),shared_procs.end());

				REPORT_STR("result:");
				for(std::set<int>::iterator it = shared_procs.begin(); it != shared_procs.end(); ++it)
					REPORT_VAL("proc ",*it);

				ir = shared_procs.find(-1);
				if( ir != shared_procs.end() )
				{
					std::cout << GetProcessorRank() << " found processor -1 " << std::endl;
					REPORT_STR(GetProcessorRank() << " found processor -1 ");
					exit(-1);
				}
			}
		}
		//Storage::integer_array procs = IntegerArrayDV(GetHandle(),tag_processors);
		//procs.clear();
		//procs.insert(procs.begin(),shared_procs.begin(),shared_procs.end());
		
		REPORT_VAL("number of processors",ret.size());
		for(size_t k = 0; k < ret.size(); ++k)
		{
			REPORT_VAL("proc ",ret[k])
		}
#if !defined(NDEBUG)
#endif //NDEBUG
#endif
		EXIT_FUNC();
		return ret;
	}

	std::vector<int> Mesh::ComputeSharedProcs(ElementType mask)
	{
		std::vector<int> ret;
		ENTER_FUNC();
#if defined(USE_MPI)
		std::set<int> shared_procs;
		int mpirank = GetProcessorRank();
		for(ElementType etype = NODE; etype <= ESET; etype = NextElementType(etype)) if( etype & mask )
		{
			REPORT_VAL("check type: ",ElementTypeName(etype));
#if defined(USE_OMP)
#pragma omp parallel 
#endif //USE_OMP
			{
				std::set<int> shared_procs_local;
#if defined(USE_OMP)
#pragma omp for
#endif //USE_OMP
				for(integer k = 0; k < LastLocalID(etype); ++k) if( isValidElement(etype,k))
				{
					Element e = ElementByLocalID(etype,k);
					Storage::integer_array procs = e.IntegerArray(ProcessorsTag());
					shared_procs_local.insert(procs.begin(),procs.end());
				}
#if defined(USE_OMP)
#pragma omp critical
#endif
				{
					shared_procs.insert(shared_procs_local.begin(),shared_procs_local.end());
				}
			}
		}
		REPORT_STR("result:");
		for(std::set<int>::iterator it = shared_procs.begin(); it != shared_procs.end(); ++it)
			REPORT_VAL("proc ",*it);
		if( !shared_procs.empty() ) 
		{
			std::set<int>::iterator ir = shared_procs.find(mpirank);
			if( ir != shared_procs.end() ) shared_procs.erase(ir);
			if( !shared_procs.empty() ) 
			{
				ret.insert(ret.end(),shared_procs.begin(),shared_procs.end());
				//Storage::integer_array procs = IntegerArrayDV(GetHandle(),tag_processors);
				//procs.clear();
				//procs.insert(procs.begin(),shared_procs.begin(),shared_procs.end());
#if !defined(NDEBUG)
				ir = shared_procs.find(-1);
				if( ir != shared_procs.end() )
				{
					std::cout << GetProcessorRank() << " found processor -1 " << std::endl;
					REPORT_STR(GetProcessorRank() << " found processor -1 ");
					exit(-1);
				}
#endif
			}
		}
		//REPORT_VAL("processors",procs.size());
		REPORT_VAL("number of processors",ret.size());
		for(size_t k = 0; k < ret.size(); ++k)
		{
			REPORT_VAL("proc ",ret[k])
		}
	
#endif
		EXIT_FUNC();
		return ret;
	}
	
	
	
	
	Mesh::proc_elements Mesh::ComputeSharedSkinSet(ElementType bridge_type, MarkerType marker)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		CheckGIDs(__FILE__, __LINE__);
		if (!HaveGlobalID(CELL))
		{
			REPORT_STR("assign cell gids");
			AssignGlobalID(CELL);
		}

		bool mesh_2d_in_3d = false;
		int mpirank = GetProcessorRank();
		Tag tag_skin = CreateTag("TEMPORARY_COMPUTE_SHARED_SKIN_SET",DATA_INTEGER,FACE,NONE);

		REPORT_STR("filling neighbouring cell's global identificators for faces")

		//for(iteratorFace it = BeginFace(); it != EndFace(); it++)
		ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for(int i = 0; i < FaceLastLocalID(); ++i) if( isValidFace(i) )
		{
			Face it = FaceByLocalID(i);
			Element::Status estat = GetStatus(it->GetHandle());
			if( estat == Element::Ghost || estat == Element::Shared )
			{
				Storage::integer_array arr = it->IntegerArray(tag_skin);
				ElementArray<Element> adj = it->getAdjElements(CELL);
        			for(ElementArray<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
				{
					assert(jt->IntegerDF(tag_owner) >= 0 && jt->IntegerDF(tag_owner) < GetProcessorsNumber());
					arr.push_back(jt->GlobalID()); //identificator of the cell
					arr.push_back(jt->IntegerDF(tag_owner)); //owner of the cell
        			}
			}
		}
		EXIT_BLOCK();

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
		//for(iteratorFace it = BeginFace(); it != EndFace(); it++)
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for (int i = 0; i < FaceLastLocalID(); ++i) if (isValidFace(i))
		{
			Face it = FaceByLocalID(i);
			bool flag = false;
			Element::Status estat1 = GetStatus(it->GetHandle()), estat2;
			if (marker && !it->GetMarker(marker)) continue;

			int ncells;
			Storage::integer_array skin_data = it->IntegerArray(tag_skin);
			if(it->GetElementDimension() == 2) 
			{
				if( estat1 == Element::Owned ) continue;
				Cell c1 = it->BackCell();
				Cell c2 = it->FrontCell();
				if( !c1.isValid() && !c2.isValid() ) continue; //no cells per face - skip hanging face
				if( skin_data.size()/2 < 2 ) continue; //only one cell per face - skip boundary face
				assert( skin_data.size()/2 == 2 ); //more then 2 cells per face - invalid grid topology
				ncells = 2;
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
			}
			else //if(it->GetElementDimension() == 1)
			{
				// attempt to generalize algorithm: more than two 2d cells neighbours on 1d face when multiple 2d meshes intersected by 1d faces.
				mesh_2d_in_3d = true;//GetDimensions() == 3;
				if( skin_data.size()/2 < 2 ) continue; //only one cell per face - skip boundary face
				ElementArray<Cell> cells = it->getCells();
				ncells = (int)cells.size();
				if( ncells != skin_data.size()/2 ) //some cells are not present on current processor
				{
					for(int k = 0; !flag && k < ncells; ++k) if(cells[k].isValid())
					{
						estat1 = cells[k].GetStatus();
						if( estat1 == Element::Owned || estat1 == Element::Shared )
							flag = true;
					}
				}
				else //all cells are present on current processor
				{
					for(int k = 0; !flag && k < ncells; ++k) if(cells[k].isValid())
					{
						estat1 = cells[k].GetStatus();
						for(int j = k+1; !flag && j < ncells; ++j) if(cells[j].isValid())
						{
							estat2 = cells[j].GetStatus();
							if( (estat1 == Element::Ghost && estat2 == Element::Shared) ||
								(estat2 == Element::Ghost && estat1 == Element::Shared) ||
								(estat1 == Element::Owned && estat2 == Element::Ghost)  ||
								(estat2 == Element::Owned && estat1 == Element::Ghost))
								flag = true;
						}
					}
				}
			}
			
			if( flag )
			{
				for(int i = 0; i < (int)skin_data.size()/2; i++) //assert checks that there are two cells
				{
					Storage::integer owner = skin_data[i*2+1]; //cell owner
					assert(owner >= 0 && owner < GetProcessorsNumber());
					if( owner != mpirank )
					{
#if defined(USE_OMP)
#pragma omp critical
#endif
						skin_faces[owner].push_back(it->GetHandle());
						break;
					}
				}	
			}
		}
		REPORT_STR("number of shared faces");
		tag_skin = DeleteTag(tag_skin);
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
			Tag on_skin = CreateTag("TEMPORARY_ON_SKIN_BRIDGE",DATA_INTEGER,bridge_type,bridge_type);

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
			
			for(element_set::iterator it = all_visited.begin(); it != all_visited.end(); it++)
			{
				Storage::integer_array os = IntegerArray(*it,on_skin);
				for(Storage::integer_array::iterator p = os.begin(); p != os.end(); p++)
				{
					assert(*p >= 0 && *p < GetProcessorsNumber());
					if( *p != mpirank )
						bridge[*p].push_back(*it);
				}
			}
			on_skin = DeleteTag(on_skin);
			REPORT_STR("number of shared elements");
			CheckGIDs(__FILE__, __LINE__);
		}
		EXIT_FUNC();
		return bridge;
#else //USE_MPI
		(void) bridge_type;
		EXIT_FUNC();
		return proc_elements();
#endif //USE_MPI
	}
	
	
	void Mesh::PackTagData(const Tag & tag, const elements_by_type & elements, int destination, ElementType mask, MarkerType select, buffer_type & buffer, TagInteger pack_position)
	{
		(void)destination;
		if( tag.GetDataType() == DATA_REMOTE_REFERENCE ) return; //NOT IMPLEMENTED TODO 14
		//ENTER_FUNC();
		//REPORT_VAL("TagName",tag.GetTagName());
#if defined(USE_MPI)
		//REPORT_VAL("Processor",destination);
		//REPORT_VAL("Buffer size before pack",buffer.size());
		ElementType pack_types[2] = {NONE,NONE};
		element_set::const_iterator eit;
		//elements_by_type pack_elements;
		//TagInteger pack_position;
		if( tag.GetDataType() == DATA_REFERENCE )
		{
		//	pack_position = GetTag("TEMP_ARRAY_POSITION");
			assert(pack_position.isValid());
		}
		buffer_type array_data_send;
		std::vector<INMOST_DATA_ENUM_TYPE> array_size_send;
		array_data_send.reserve(4096);
		array_size_send.reserve(4096);
		INMOST_DATA_ENUM_TYPE size = tag.GetSize();
		for(int i = ElementNum(NODE); i <= ElementNum(MESH); i++) if( (mask & ElementTypeFromDim(i)) && tag.isDefinedByDim(i) )
		{
			pack_types[0] |= ElementTypeFromDim(i);
			//REPORT_VAL("select marker",select);
			//REPORT_VAL("pack for type",ElementTypeName(ElementTypeFromDim(i)));
			//REPORT_VAL("number of elements",elements[i].size());
			int total_packed = 0;
			if( tag.isSparseByDim(i) )
			{
				//REPORT_VAL("sparse for type", ElementTypeName(ElementTypeFromDim(i)));
				pack_types[1] |= ElementTypeFromDim(i);
				INMOST_DATA_ENUM_TYPE count = (INMOST_DATA_ENUM_TYPE)array_size_send.size();
				array_size_send.push_back(0);
				for(eit = elements[i].begin(); eit != elements[i].end(); eit++)
				{
					if( (!select || GetAnyMarker(*eit,select)) && HaveData(*eit,tag) )
					{
						//REPORT_VAL("element index", static_cast<INMOST_DATA_ENUM_TYPE>(eit-elements[i].begin()));
						array_size_send.push_back(static_cast<INMOST_DATA_ENUM_TYPE>(eit-elements[i].begin()));
						array_size_send[count]++;
						INMOST_DATA_ENUM_TYPE s = GetDataSize(*eit,tag);
						size_t had_s = array_data_send.size();
						if( s )
						{
							array_data_send.resize(had_s+GetDataCapacity(*eit,tag));
							if( tag.GetDataType() == DATA_REFERENCE )
							{
								reference_array refs = ReferenceArray(*eit, tag);
								INMOST_DATA_ENUM_TYPE bytes = tag.GetBytesSize();
								for(Storage::reference_array::size_type i = 0; i < refs.size(); ++i)
								{
									HandleType data = InvalidHandle();
									if ( refs[i].isValid() && refs[i].HaveData(pack_position) ) 
										data = ComposeHandle(refs[i]->GetElementType(), pack_position[refs[i]]);
									memcpy(&array_data_send[had_s+(size_t)i*bytes],&data,sizeof(HandleType));
								}
							}
							else GetData(*eit,tag,0,s,&array_data_send[had_s]);
						}
						if( size == ENUMUNDEF ) array_size_send.push_back(s);
						++total_packed;
					}
				}
				//REPORT_VAL("index", count);
				//REPORT_VAL("count",array_size_send[count]);
			}
			else
			{
				//REPORT_VAL("dense for type", ElementTypeName(ElementTypeFromDim(i)));
				for(eit = elements[i].begin(); eit != elements[i].end(); eit++) if( !select || GetAnyMarker(*eit,select) )
				{
					INMOST_DATA_ENUM_TYPE s = GetDataSize(*eit,tag);
					size_t had_s = array_data_send.size();
					if( s )
					{
						array_data_send.resize(had_s+GetDataCapacity(*eit,tag));
						if (tag.GetDataType() == DATA_REFERENCE)
						{
							reference_array refs = ReferenceArray(*eit, tag);
							INMOST_DATA_ENUM_TYPE bytes = tag.GetBytesSize();
							for(Storage::reference_array::size_type i = 0; i < refs.size(); ++i)
							{
								HandleType data = InvalidHandle();
								if ( refs[i].isValid() && refs[i].HaveData(pack_position) ) 
									data = ComposeHandle(refs[i]->GetElementType(), pack_position[refs[i]]);
								memcpy(&array_data_send[had_s+(size_t)i*bytes],&data,sizeof(HandleType));
							}
						}
						else GetData(*eit,tag,0,s,&array_data_send[had_s]);
					}
					if( size == ENUMUNDEF ) array_size_send.push_back(s);
					++total_packed;
				}
			}
			//REPORT_VAL("total packed records",total_packed);
			//REPORT_VAL("Buffer size after pack for type", buffer.size());
		}
		//REPORT_VAL("tag defined on",static_cast<int>(pack_types[0]));
		//REPORT_VAL("tag sparse on",static_cast<int>(pack_types[1]));
		//REPORT_VAL("size_send",array_size_send.size());
		//REPORT_VAL("data_send",array_data_send.size());
		
		pack_data_array(buffer,pack_types,2,GetCommunicator());
		pack_data_vector(buffer,array_size_send,GetCommunicator());
		pack_data_vector(buffer,array_data_send,GetCommunicator());
		//REPORT_VAL("Buffer size after pack",buffer.size());
#else
		(void) tag;
		(void) elements;
        (void) destination;
		(void) mask;
		(void) select;
		(void) buffer;
#endif
		//REPORT_VAL("TagName",tag.GetTagName());
		//EXIT_FUNC();
	}
	

	
	HandleType Mesh::FindSharedGhost(ElementType etype, Storage::integer global_id, int source_proc, int owner_proc)
    {
#if defined(USE_PARALLEL_STORAGE)
		int k = ElementNum(etype), mpirank = GetProcessorRank();
		element_set::iterator search;
		if( owner_proc == mpirank )
		{
			element_set & shared = shared_elements[source_proc][k];
			//for(unsigned q = 0; q < shared.size(); ++q)
			//	if( GlobalID(shared[q]) == global_id )
			//	{
					//std::cout << "found in shared, proc " << source_proc << " " << ElementTypeName(etype) << ":" << GetHandleID(shared[q]) << " gid " << global_id << std::endl;
			//		return shared[q];
			//	}
			search = std::lower_bound(shared.begin(),shared.end(),global_id,Mesh::GlobalIDComparator(this));
			if( search != shared.end() && GlobalID(*search) == global_id )
				return *search;
		}
		else
		{
			element_set & ghost = ghost_elements[owner_proc][k];
			//for(unsigned q = 0; q < ghost.size(); ++q)
			//	if( GlobalID(ghost[q]) == global_id )
			//	{
					//std::cout << "found in ghost " << owner_proc << " " << ElementTypeName(etype) << ":" << GetHandleID(ghost[q]) << " gid " << global_id << std::endl;
			//		return ghost[q];
			//	}
			search = std::lower_bound(ghost.begin(),ghost.end(),global_id,Mesh::GlobalIDComparator(this));
			if( search != ghost.end() && GlobalID(*search) == global_id )
				return *search;
		}
#endif
		//std::cout << "not found source " << source_proc << " owner " << owner_proc << " " << ElementTypeName(etype) << " gid " << global_id << std::endl;
		return InvalidHandle();
     }
	


	void Mesh::UnpackTagData(const Tag & tag, const elements_by_type & elements, int source, ElementType mask, MarkerType select, buffer_type & buffer, size_t & buffer_position, ReduceOperation op, const elements_by_type & unpack_elements)//, proc_elements_by_type * send_elements)
	{
        (void) mask;
		//if( tag.GetDataType() == DATA_REMOTE_REFERENCE) return; //NOT IMPLEMENTED TODO 14
		ENTER_FUNC();
		REPORT_VAL("TagName",tag.GetTagName());
		REPORT_VAL("select marker",select);
#if defined(USE_MPI)
		REPORT_VAL("Position before unpack",buffer_position);
		if( !buffer.empty() )
		{
			INMOST_DATA_ENUM_TYPE pos = 0, k = 0;
			ElementType recv_mask[2] = {NONE,NONE};
			element_set::const_iterator eit;
			//elements_by_type unpack_elements;
			INMOST_DATA_ENUM_TYPE size = tag.GetSize();
			buffer_type array_data_recv;
			std::vector<INMOST_DATA_ENUM_TYPE> array_size_recv;
			unpack_data_array(buffer,buffer_position,recv_mask,2,GetCommunicator());
			unpack_data_vector(buffer,buffer_position,array_size_recv,GetCommunicator());
			unpack_data_vector(buffer,buffer_position,array_data_recv,GetCommunicator());
			REPORT_VAL("size of size array",array_size_recv.size());
			REPORT_VAL("size of data array",array_data_recv.size());
			REPORT_VAL("Position after unpack",buffer_position);
			for(int i = ElementNum(NODE); i <= ElementNum(MESH); i++) if( (recv_mask[0] & ElementTypeFromDim(i)) )
			{
				REPORT_VAL("unpack for type",ElementTypeName(ElementTypeFromDim(i)));
				REPORT_VAL("number of elements",elements[i].size());
				if( !tag.isDefinedByDim(i) ) 
				{
					CreateTag(tag.GetTagName(),tag.GetDataType(),ElementTypeFromDim(i),recv_mask[1] & ElementTypeFromDim(i),size);
				}
				int total_unpacked = 0;
				int total_skipped = 0;
				INMOST_DATA_ENUM_TYPE array_size;
				if( tag.isSparseByDim(i) )
				{
					REPORT_VAL("sparse for type",ElementTypeName(ElementTypeFromDim(i)));
					REPORT_VAL("index", k);
					REPORT_VAL("count", array_size_recv[k]);
					REPORT_VAL("pos",pos);
					INMOST_DATA_ENUM_TYPE count = array_size_recv[k++];
					for(INMOST_DATA_ENUM_TYPE j = 0; j < count; j++)
					{
						assert( k < array_size_recv.size() );
						assert( array_size_recv[k] < elements[i].size() );
						//REPORT_VAL("element index", array_size_recv[k]);
						eit = elements[i].begin() + array_size_recv[k++];
						assert( !select || GetAnyMarker(*eit,select) ); //if fires then very likely that marker was not synchronized
						array_size = size;
						if (array_size == ENUMUNDEF)
						{
							assert(k < array_size_recv.size());
							array_size = array_size_recv[k++];
						}
						if( tag.GetDataType() == DATA_REFERENCE ) //change received content with actual reference to element
						{
							for (INMOST_DATA_ENUM_TYPE i = 0; i < array_size; i++)
							{
								INMOST_DATA_ENUM_TYPE spos = pos + i * tag.GetBytesSize();
								INMOST_DATA_ENUM_TYPE dpos = ENUMUNDEF;
								HandleType * data = (HandleType*)(&array_data_recv[spos]);
								if( *data != InvalidHandle() ) 
								{
									dpos = GetHandleID(*data);
									*data = unpack_elements[GetHandleElementNum(*data)][dpos];
								}
							}
						}
						if (array_size)
						{
							INMOST_DATA_ENUM_TYPE data_size = GetDataCapacity(&array_data_recv[pos], array_size, tag);
							assert(pos + data_size <= array_data_recv.size());
							op(tag, Element(this, *eit), &array_data_recv[pos], array_size);
							pos += data_size;
						}
						++total_unpacked;
					}
				}
				else
				{
					REPORT_VAL("dense for type",ElementTypeName(ElementTypeFromDim(i)));
					for(eit = elements[i].begin(); eit != elements[i].end(); eit++) 
					{
						if( !select || GetAnyMarker(*eit,select) )
						{
							array_size = size;
							if (array_size == ENUMUNDEF)
							{
								assert(k < array_size_recv.size());
								array_size = array_size_recv[k++];
							}
							if (tag.GetDataType() == DATA_REFERENCE) //convert input offsets into real handles
							{
								for (INMOST_DATA_ENUM_TYPE i = 0; i < array_size; i++)
								{
									INMOST_DATA_ENUM_TYPE spos = pos + i * tag.GetBytesSize();
									INMOST_DATA_ENUM_TYPE dpos = ENUMUNDEF;
									HandleType * data = (HandleType*)(&array_data_recv[spos]);
									if( *data != InvalidHandle() ) 
									{
										dpos = GetHandleID(*data);
										*data = unpack_elements[GetHandleElementNum(*data)][dpos];
									}
								}
							}
							if (array_size)
							{
								INMOST_DATA_ENUM_TYPE data_size = GetDataCapacity(&array_data_recv[pos], array_size, tag);
								assert(pos + data_size <= array_data_recv.size());
								op(tag, Element(this, *eit), &array_data_recv[pos], array_size);
								pos += data_size;
							}
							++total_unpacked;
						} 
						else ++total_skipped;
					}
				}
				REPORT_VAL("total skipped",total_skipped);
				REPORT_VAL("total unpacked records",total_unpacked);
				REPORT_VAL("Position after unpack for type", buffer_position);
			}
		}
		REPORT_VAL("Position after unpack",buffer_position);
#else
		(void) tag;
		(void) elements;
		(void) mask;
		(void) select;
		(void) buffer;
		(void) buffer_position;
		(void) op;
#endif //USE_MPI
		REPORT_VAL("TagName",tag.GetTagName());
		EXIT_FUNC();
	}
	
	
	void Mesh::ExchangeDataInnerBegin(const tag_set & tags, const parallel_storage & from, const parallel_storage & to, ElementType mask, MarkerType select, exchange_data & storage)
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
		//int mpirank = GetProcessorRank(),mpisize = GetProcessorsNumber();
		//Storage::integer_array procs = IntegerArrayDV(GetHandle(),tag_processors);
		//Storage::integer_array::iterator p = procs.begin();
		std::vector<int> procs = ComputeSharedProcs(from,to);
		std::vector<int>::iterator p = procs.begin();
		REPORT_VAL("processors",procs.size());
		storage.send_buffers.resize(procs.size());
		storage.recv_buffers.resize(procs.size());
		REPORT_VAL("send buffers size",storage.send_buffers.size());
		REPORT_VAL("recv buffers size",storage.recv_buffers.size());
		std::vector<int> done;
		parallel_storage::const_iterator find;
		std::vector<size_t> send_size(procs.size(),0), recv_size(procs.size(),0);
		
		bool unknown_size = false;
		bool have_reference_tag = false;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < tags.size(); k++)
		{
			if( tags[k].GetSize() == ENUMUNDEF
#if defined(USE_AUTODIFF)
			   || tags[k].GetDataType() == DATA_VARIABLE
#endif
			   || tags[k].GetDataType() == DATA_REFERENCE
			   || tags[k].GetDataType() == DATA_REMOTE_REFERENCE
			   ) unknown_size = true;
			//for(int i = 0; i < 5; ++i)
			//	if( (mask & ElementTypeFromDim(i)) && tags[k].isSparseByDim(i) )
			//		unknown_size = true;
			if( tags[k].GetDataType() == DATA_REFERENCE )
				have_reference_tag = true;
		}
		//precompute sizes
		ENTER_BLOCK();
		REPORT_STR("precompute sizes");
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
		EXIT_BLOCK();
		//int num_send = 0, num_recv = 0;
		//TagInteger pack_position;
		tag_set tag_pack_list;
		if( have_reference_tag )
		{
			//ENTER_BLOCK();
			//REPORT_STR("Create tag");
			//pack_position = CreateTag("PROTECTED_TEMPORARY_PACK_POSITION",DATA_INTEGER,ESET|CELL|FACE|EDGE|NODE,ESET|CELL|FACE|EDGE|NODE,1);
			//EXIT_BLOCK();
			ENTER_BLOCK();
			REPORT_STR("List tags");
			ListTags(tag_pack_list);
			{
				tag_set::iterator it = tag_pack_list.begin();
				while( it != tag_pack_list.end() )
				{
					if( it->GetTagName().substr(0,9) == "PROTECTED" ) 
						it = tag_pack_list.erase(it);
					//else if(it->GetDataType() == DATA_REFERENCE)
					//	it = tag_pack_list.erase(it);
					else if(it->GetDataType() == DATA_REMOTE_REFERENCE)
						it = tag_pack_list.erase(it);
					else it++;
				}
			}
			EXIT_BLOCK();
		}
		ENTER_BLOCK();
#if defined(USE_OMP)
//#pragma omp parallel
#endif
		{
			TagInteger pack_position;
			if (have_reference_tag)
			{
				std::stringstream tag_name;
				tag_name << "PROTECTED_TEMPORARY_PACK_POSITION_" << GetLocalProcessorRank();
				pack_position = CreateTag(tag_name.str(), DATA_INTEGER, ESET | CELL | FACE | EDGE | NODE, ESET | CELL | FACE | EDGE | NODE, 1);
				//pack_position = CreateTag(tag_name.str(), DATA_INTEGER, ESET | CELL | FACE | EDGE | NODE, NONE, 1);
			}
#if defined(USE_OMP)
//#pragma omp for schedule(dynamic,1)
#endif
			for(int pos = 0; pos < (int)procs.size(); ++pos)
			{
				std::vector<int>::iterator p = procs.begin() + pos;
				if (send_size[pos])
				{
					elements_by_type selems;
					if (have_reference_tag)
					{
						PackElementsGather(selems, from.find(*p)->second, *p, mask, select, tags, true, true, true);
						PackElementsEnumerate(selems, pack_position);
						PackElementsData(selems, storage.send_buffers[pos].second, *p, tag_pack_list, pack_position, true);
					}
					for (INMOST_DATA_ENUM_TYPE k = 0; k < tags.size(); k++)
						PackTagData(tags[k], from.find(*p)->second, *p, mask, select, storage.send_buffers[pos].second, pack_position);
					if (have_reference_tag)
						PackElementsUnenumerate(selems, pack_position);
					storage.send_buffers[pos].first = *p;
				}
				else storage.send_buffers[pos].first = -1;
				if (recv_size[pos])
				{
					if (!unknown_size)
					{
						size_t buffer_size = 0, n = recv_size[pos];
						for (INMOST_DATA_ENUM_TYPE k = 0; k < tags.size(); k++)
						{
							buffer_size += pack_data_array_size<ElementType>(2, GetCommunicator());
							buffer_size += pack_data_vector_size<INMOST_DATA_ENUM_TYPE>(6 + 2 * n, GetCommunicator());
							buffer_size += pack_data_vector_size<INMOST_DATA_BULK_TYPE>(n * tags[k].GetSize() * tags[k].GetPackedBytesSize(), GetCommunicator());
						}
						storage.recv_buffers[pos].second.resize(buffer_size);
					}
					storage.recv_buffers[pos].first = *p;
				}
				else storage.recv_buffers[pos].first = -1;
			}
			if (have_reference_tag)
				DeleteTag(pack_position);
		}
		EXIT_BLOCK();
		ENTER_BLOCK();
		{
			exch_buffer_type::iterator it = storage.send_buffers.begin();
			while (it != storage.send_buffers.end())
			{
				if (it->first == -1)
					it = storage.send_buffers.erase(it);
				else ++it;
			}
		}
		EXIT_BLOCK();
		ENTER_BLOCK();
		{
			exch_buffer_type::iterator it = storage.recv_buffers.begin();
			while (it != storage.recv_buffers.end())
			{
				if (it->first == -1)
					it = storage.recv_buffers.erase(it);
				else ++it;
			}
		}
		EXIT_BLOCK();
		if( unknown_size ) 
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
	
	
	void Mesh::InformElementsOwners(proc_elements_by_type & send_elements, exchange_data & storage)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		int mpirank = GetProcessorRank();
		std::vector<int> done;
		storage.send_buffers.clear();
		storage.send_reqs.clear();
		storage.recv_buffers.clear();
		storage.recv_reqs.clear();
		storage.send_buffers.resize(send_elements.size());
		
		ENTER_BLOCK();
		
		ENTER_BLOCK();
		REPORT_STR("Pack elements");
#if defined(USE_OMP)
//#pragma omp parallel
#endif
		{
			TagInteger pack_position;
			std::stringstream tag_name;
			tag_name << "PROTECTED_TEMPORARY_PACK_POSITION_" << GetLocalProcessorRank();
			pack_position = CreateTag(tag_name.str(), DATA_INTEGER, ESET | CELL | FACE | EDGE | NODE, ESET | CELL | FACE | EDGE | NODE, 1);

#if defined(USE_OMP)
//#pragma omp for schedule(dynamic,1)
#endif
			for (int i = 0; i < (int)send_elements.size(); ++i)//for(proc_elements_by_type::iterator it = send_elements.begin(); it != send_elements.end(); it++)
			{
				proc_elements_by_type::iterator it = send_elements.begin();
				{ int k = i; while (k) ++it, --k; } //advance
				if (!it->second.empty())
				{
					PackElementsGather(it->second, it->second, it->first, NONE, 0, tag_set(), false, true, false);
					PackElementsEnumerate(it->second, pack_position);
					PackElementsData(it->second, storage.send_buffers[i].second, it->first, tag_set(), pack_position, false);
					PackElementsUnenumerate(it->second, pack_position);
					storage.send_buffers[i].first = it->first;
				}
				else storage.send_buffers[i].first = -1;
			}
			DeleteTag(pack_position);
		}
		EXIT_BLOCK();
		ENTER_BLOCK();
		{
			exch_buffer_type::iterator it = storage.send_buffers.begin();
			while (it != storage.send_buffers.end())
			{
				if (it->first == -1)
					it = storage.send_buffers.erase(it);
				else ++it;
			}
		}
		EXIT_BLOCK();
		EXIT_BLOCK();

		PrepareReceiveInner(UnknownSource, storage.send_buffers,storage.recv_buffers);
		ExchangeBuffersInner(storage.send_buffers,storage.recv_buffers,storage.send_reqs,storage.recv_reqs);
		
		send_elements.clear();
		
		while( !(done = FinishRequests(storage.recv_reqs)).empty() )
		{
			for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
			{
				tag_set tag_list_recv;
				elements_by_type recv_elements;
				size_t position = 0;
				UnpackElementsData(recv_elements,storage.recv_buffers[*qt].second,storage.recv_buffers[*qt].first,position,tag_list_recv);
				for(int k = 0; k < 5; ++k)
					for(element_set::iterator it = recv_elements[k].begin(); it != recv_elements[k].end(); it++)
					{
						if( IntegerDF(*it,tag_owner) == -1 )
						{
							std::cout << __FILE__ << ":" << __LINE__ << " proc " << GetProcessorRank() << " receive from " << storage.recv_buffers[*qt].first;
							std::cout << " element " << ElementTypeName(GetHandleElementType(*it)) << ":" << GetHandleID(*it);
							std::cout << std::endl;
							exit(-1);
						}
						Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);
						Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),storage.recv_buffers[*qt].first);
						if( ip == proc.end() || (*ip) != storage.recv_buffers[*qt].first ) proc.insert(ip,storage.recv_buffers[*qt].first);
						if( IntegerDF(*it,tag_owner) == mpirank )
							SetStatus(*it,Element::Shared);
						else
							SetStatus(*it,Element::Ghost);
					}
			}
		}
		//wait for second round of communication
		if( !storage.send_reqs.empty() )
		{
			REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(storage.send_reqs.size()),&storage.send_reqs[0],MPI_STATUSES_IGNORE));
		}
#endif
		EXIT_FUNC();
	}
	
	void Mesh::ExchangeDataInnerEnd(const tag_set & tags, const parallel_storage & from, const parallel_storage & to, ElementType mask, MarkerType select, ReduceOperation op, exchange_data & storage)
	{
		(void) from;
		bool have_reference_tag = false;
		{ // checks for bad input
			if( mask == NONE ) return;
			for(tag_set::const_iterator it = tags.begin(); it != tags.end(); ++it)
			{
				assert( it->isValid() );
				if( it->GetDataType() == DATA_REFERENCE )
					have_reference_tag = true;
			}
			if( tags.empty() ) return;
		}
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
		int mpirank = GetProcessorRank();
		std::vector<int> done;
		proc_elements_by_type send_elements; //for DATA_REFERENCE
		while( !(done = FinishRequests(storage.recv_reqs)).empty() )
		{
			for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
			{
				tag_set tag_list_recv;
				elements_by_type selems;
				size_t position = 0;
				if( have_reference_tag )
				{
					UnpackElementsData(selems,storage.recv_buffers[*qt].second,storage.recv_buffers[*qt].first,position,tag_list_recv);
					for(int k = 0; k < 5; ++k)
						for(element_set::iterator it = selems[k].begin(); it != selems[k].end(); it++)
						{
							Storage::integer owner = IntegerDF(*it,tag_owner);
							if( owner == mpirank ) continue;
							Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);
							Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),mpirank);
							if( ip == proc.end() || (*ip) != mpirank )
							{
								proc.insert(ip,mpirank);
								send_elements[owner][GetHandleElementNum(*it)].push_back(*it);
								//std::cout << "inform " << owner << " of " << ElementTypeName(GetHandleElementType(*it)) << " from " << mpirank << std::endl;
							}
						}
				}
				for(INMOST_DATA_ENUM_TYPE k = 0; k < tags.size(); k++)
				{
					REPORT_VAL("processor",storage.recv_buffers[*qt].first);
					UnpackTagData(tags[k],to.find(storage.recv_buffers[*qt].first)->second,storage.recv_buffers[*qt].first,mask,select,storage.recv_buffers[*qt].second,position,op, selems);
				}
			}
		}
		if( !storage.send_reqs.empty() )
		{
			REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(storage.send_reqs.size()),&storage.send_reqs[0],MPI_STATUSES_IGNORE));
		}
		//additional round of communication to tell the tag owner about exchange of elements
		if( have_reference_tag )
		{
			InformElementsOwners(send_elements,storage);
			CheckProcsSorted(__FILE__,__LINE__);
			
			RecomputeParallelStorage(ESET | CELL | FACE | EDGE | NODE);
			
			CheckGhostSharedCount(__FILE__,__LINE__);
			CheckCentroids(__FILE__,__LINE__);
			
			
			ExchangeData(tag_processors,ESET | CELL | FACE | EDGE | NODE,0);

			CheckGhostSharedCount(__FILE__, __LINE__);
			CheckCentroids(__FILE__, __LINE__);
			CheckOwners(__FILE__, __LINE__);
			CheckGIDs(__FILE__, __LINE__);
			CheckProcessors();
			//ComputeSharedProcs();
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
		int mpisize = GetProcessorsNumber();
		(void)mpisize;
		//int err = 0;
		REPORT_STR("gather ghost and shared elements")
		double time = Timer();
		for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); it++)
		{
			Element::Status estat = it->GetStatus();
			if( it->Hidden() )
			{
				// err++;
				REPORT_STR(GetProcessorRank() <<  " " << __FILE__ << ":" << __LINE__ << "Adding hidden element " << ElementTypeName(it->GetElementType()));
				std::cout << GetProcessorRank() <<  " " << __FILE__ << ":" << __LINE__ <<  "Adding hidden element " << ElementTypeName(it->GetElementType()) << std::endl;
			}
			if( estat == Element::Shared ) 
			{ 
				Storage::integer_array v = it->IntegerArrayDV(tag_processors);
				for(Storage::integer_array::iterator vit = v.begin(); vit != v.end(); vit++)
				{
					//~ if( !(*vit >= 0 && *vit < mpisize) )
					//~ {
						//~ // err++;
						//~ REPORT_STR(GetProcessorRank() <<  " " << __FILE__ << ":" << __LINE__ <<  " bad proc in list " << *vit << " " << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() );
						//~ std::cout << GetProcessorRank() << " " << __FILE__ << ":" << __LINE__ << " bad proc in list " << *vit << " " << ElementTypeName(it->GetElementType()) << ":" << it->LocalID()  << std::endl;
					//~ }
					assert(*vit >= 0 && *vit < mpisize);
					if( *vit != mpirank )
						shared[*vit][ElementNum(it->GetElementType())].push_back(*it);
				}
			}
			else if( estat == Element::Ghost )
			{
				int proc = it->IntegerDF(tag_owner);
				//~ if( !(proc >= 0 && proc < mpisize) )
				//~ {
					//~ // err++;
					//~ REPORT_STR(GetProcessorRank() <<  " " << __FILE__ << ":" << __LINE__ <<  " bad proc owner " << proc << " " << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() );
					//~ std::cout << GetProcessorRank() <<  " " << __FILE__ << ":" << __LINE__ <<  " bad proc owner " << proc << " " << ElementTypeName(it->GetElementType()) << ":" << it->LocalID()  << std::endl;
				//~ }
				assert(proc >= 0 && proc < mpisize);
				ghost[proc][ElementNum(it->GetElementType())].push_back(*it);
			}
		}
		// err = Integrate(err);
		// if( err )
		// {
		// 	std::cout << "crash" << std::endl;
		// 	exit(-1);
		// }
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

	void Mesh::ExchangeOrientedData(const Tag& tag, ElementType mask, MarkerType select, MarkerType orient)
	{
		ExchangeOrientedData(tag_set(1, tag), mask, select, orient);
	}

	void Mesh::ExchangeOrientedData(const tag_set& tags, ElementType mask, MarkerType select, MarkerType orient)
	{
		if (m_state == Serial) return;
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements, shared_elements, mask);
#endif //USE_PARALLEL_STORAGE
		{
			ReportParallelStorage();
			exchange_data storage;
			ExchangeDataInnerBegin(tags, shared_elements, ghost_elements, mask, select, storage);
			ExchangeDataInnerEnd(tags, shared_elements, ghost_elements, mask, select, DefaultUnpack, storage);
		}
		bool clear_marker = false;
		if (!orient)
		{
			orient = CreateMarker();
			MarkNormalOrientation(orient);
			clear_marker = true;
		}
		for (tag_set::const_iterator it = tags.begin(); it != tags.end(); ++it)
			for (tag_set::iterator jt = orient_tags.begin(); jt != orient_tags.end(); ++jt)
				if (*it == *jt)
				{
					for (integer kt = 0; kt < FaceLastLocalID(); ++kt) if (isValidFace(kt))
					{
						Face f = FaceByLocalID(kt);
						if (f.GetMarker(orient)) OrientTag(f, *it);
					}
					break;
				}
		if (clear_marker)
		{
			ReleaseMarker(orient, FACE);
			orient = 0;
		}
		EXIT_FUNC();
#else //USE_MPI
		(void)tags;
		(void)mask;
		(void)select;
		(void)orient;
#endif //USE_MPI
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
		ReportParallelStorage();
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
		ReportParallelStorage();
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

	void Mesh::PackElementsEnumerate(elements_by_type & selems, TagInteger pack_position)
	{
		for(int etypenum = ElementNum(NODE); etypenum <= ElementNum(ESET); ++etypenum)
		{
			INMOST_DATA_INTEGER_TYPE q = 0;
			for(element_set::iterator it = selems[etypenum].begin(); it != selems[etypenum].end(); it++)
			{
				pack_position[*it] = q++;
				assert(GetHandleElementNum(*it) == etypenum);
				assert(*it == selems[etypenum][pack_position[*it]]);
			}
		}
	}
	void Mesh::PackElementsUnenumerate(elements_by_type & selems, TagInteger pack_position)
	{
		for(int etypenum = ElementNum(NODE); etypenum <= ElementNum(ESET); ++etypenum)
		{
			for(element_set::iterator it = selems[etypenum].begin(); it != selems[etypenum].end(); it++)
				DelData(*it,pack_position);
		}
	}

	void Mesh::PackElementsGather(elements_by_type & selems, const elements_by_type & elements, int destination, ElementType mask, MarkerType select, const tag_set & tag_list, bool force_send, bool send_links_to_owner, bool pack_by_gids)
	{
		ENTER_FUNC();
		if( !allow_pack_by_gids ) pack_by_gids = false;

		MarkerType busy = CreatePrivateMarker();

		REPORT_STR("initial number of elements");
		REPORT_VAL("NODE",selems[0].size());
		REPORT_VAL("EDGE",selems[1].size());
		REPORT_VAL("FACE",selems[2].size());
		REPORT_VAL("CELL",selems[3].size());
		REPORT_VAL("ESET",selems[4].size());
		
		
		for(int etypenum = ElementNum(NODE); etypenum <= ElementNum(CELL); ++etypenum) if( !selems[etypenum].empty() )
			SetPrivateMarkerArray(&selems[etypenum][0],(INMOST_DATA_ENUM_TYPE)selems[etypenum].size(),busy);



		//Add all elements from tags
		if( send_links_to_owner )
		for(INMOST_DATA_ENUM_TYPE k = 0; k < tag_list.size(); ++k)
		{
			const Tag & tag = tag_list[k];
			if( tag.GetDataType() == DATA_REFERENCE )
			{
				for(int i = ElementNum(NODE); i <= ElementNum(ESET); i++) if( (mask & ElementTypeFromDim(i)) && tag.isDefinedByDim(i) )
				{
					INMOST_DATA_ENUM_TYPE last = (INMOST_DATA_ENUM_TYPE)elements[i].size(); //avoid resize if selems == elements
					for(INMOST_DATA_ENUM_TYPE q = 0; q < last; ++q) 
					{
						HandleType eit = elements[i][q];
						if( !select || GetMarker(eit,select) )
						{
							if( !tag.isSparseByDim(i) || HaveData(eit,tag) )
							{
								reference_array refs = ReferenceArray(eit, tag);
								for(Storage::reference_array::size_type i = 0; i < refs.size(); ++i) if ( refs[i].isValid() ) 
								{
									if( !refs[i].Hidden() && !refs[i].GetPrivateMarker(busy) )
									{
										bool add = force_send;
										if( !add )
										{
											integer_array procs = refs[i].IntegerArray(ProcessorsTag());
											if( std::binary_search(procs.begin(),procs.end(),destination) )
												add = true;
										}
										if( add )
										{
											selems[refs[i].GetElementNum()].push_back(refs[i].GetHandle());
											refs[i]->SetPrivateMarker(busy);
										}
									} //if
								} //for
							} //if
						} //if
					} //for
				} //for
			} //if
		} //for k

		// Add low conns elements to selems array.
		// Low conns for ESET can contain any element
		if( send_links_to_owner )
		for(element_set::iterator it = selems[4].begin(); it != selems[4].end(); it++)
		{
			Element::adj_type const & adj = LowConn(*it);
			for(Element::adj_type::const_iterator jt = adj.begin(); jt != adj.end(); jt++)
			{
				if( *jt != InvalidHandle() && !Hidden(*jt) && !GetPrivateMarker(*jt,busy) )
				{
					Storage::integer_array procs = IntegerArrayDV(*jt,ProcessorsTag());
					if( std::binary_search(procs.begin(),procs.end(),destination) )
					{
						selems[GetHandleElementNum(*jt)].push_back(*jt);
						SetPrivateMarker(*jt,busy);
					}
				}
			}
		}

		// Low conns for CELLs, FACEs, EDGEs
		for(int etypenum = ElementNum(CELL); etypenum > ElementNum(NODE); --etypenum)
		{
			//selems[etypenum-1].reserve(selems[etypenum-1].size()+size_hint[etypenum]*selems[etypenum].size());
			for(element_set::iterator it = selems[etypenum].begin(); it != selems[etypenum].end(); it++)
			{
				//if element is already found on remote processor, then don't have to add it's connections
				bool add_connections = true;
				if( pack_by_gids && HaveGlobalID(ElementTypeFromDim(etypenum)) )
				{
					Storage::integer_array procs = IntegerArray(*it,tag_processors);
					if( std::binary_search(procs.begin(),procs.end(),destination) )
						add_connections = false;
				}
				if( add_connections )
				{
					Element::adj_type const & adj = LowConn(*it);
					for(Element::adj_type::const_iterator jt = adj.begin(); jt != adj.end(); jt++)
					{
						if( !Hidden(*jt) && !GetPrivateMarker(*jt,busy) )
						{
							selems[etypenum-1].push_back(*jt);
							SetPrivateMarker(*jt,busy);
						}
					}
					if (etypenum == ElementNum(CELL) && HighConnTag().isDefined(CELL))
					{
						Element::adj_type const& adj = HighConn(*it);
						for (Element::adj_type::const_iterator jt = adj.begin(); jt != adj.end(); jt++)
						{
							if (!Hidden(*jt) && !GetPrivateMarker(*jt, busy))
							{
								selems[0].push_back(*jt);
								SetPrivateMarker(*jt, busy);
							}
						}
					}
				}
			}
		}

		// Recursevely looking for all high connections for ESETs
		if( !selems[4].empty() ) SetPrivateMarkerArray(&selems[4][0],(INMOST_DATA_ENUM_TYPE)selems[4].size(),busy);
		
		for(INMOST_DATA_ENUM_TYPE i = 0; i < selems[4].size(); ++i)
		{
			ElementSet set(this,selems[4][i]);
			if( set.HaveParent() )
			{
				//if set is on the remote processor, then don't have to reconstruct the tree
				if( !set.GetParent().GetPrivateMarker(busy) && !set.GetParent().Hidden() )
				{
					//Storage::integer_array procs = set.GetParent().IntegerArrayDV(ProcessorsTag());
					//if( std::binary_search(procs.begin(),procs.end(),destination) )
					{
						selems[4].push_back(set.GetParent().GetHandle());
						set.GetParent().SetPrivateMarker(busy);
					}
				}
			}
			//add all children
			if( send_links_to_owner )
			if( set.HaveChild() )
			{
				for(ElementSet jt = set.GetChild(); jt.isValid(); jt = jt.GetSibling() )
				{
					if( !jt.GetPrivateMarker(busy) && !jt.Hidden() )
					{
						Storage::integer_array procs = jt.IntegerArrayDV(ProcessorsTag());
						if( std::binary_search(procs.begin(),procs.end(),destination) )
						{
							selems[4].push_back(jt.GetHandle());
							jt.SetPrivateMarker(busy);
						}
					}
				}
			}
		}

		REPORT_STR("final number of elements");
		REPORT_VAL("NODE",selems[0].size());
		REPORT_VAL("EDGE",selems[1].size());
		REPORT_VAL("FACE",selems[2].size());
		REPORT_VAL("CELL",selems[3].size());
		REPORT_VAL("ESET",selems[4].size());
		
		
		
		for(int etypenum = ElementNum(NODE); etypenum <= ElementNum(ESET); ++etypenum) if( !selems[etypenum].empty() )
			RemPrivateMarkerArray(&selems[etypenum][0],(INMOST_DATA_ENUM_TYPE)selems[etypenum].size(),busy);

		ReleasePrivateMarker(busy);

		EXIT_FUNC();
	}
	
	void Mesh::PackElementsData(elements_by_type & selems, buffer_type & buffer, int destination, const tag_set & tag_list, TagInteger arr_position, bool pack_by_gids)
	{
		ENTER_FUNC();
		if( !allow_pack_by_gids ) pack_by_gids = false;
		REPORT_VAL("dest",destination);
#if defined(USE_MPI)
		REPORT_VAL("buffer size",buffer.size());
		
		
		MarkerType pack_tags_mrk = CreatePrivateMarker();
		char orient = false;
		if (!orient_tags.empty())
		{
			//check if there is any oriented tag to be exchanged
			bool found = false;
			for (tag_set::const_iterator it = tag_list.begin(); it != tag_list.end() && !found; ++it)
				for (tag_set::iterator jt = orient_tags.begin(); jt != orient_tags.end() && !found; ++jt)
					if (*it == *jt) found = true;
			if (found)
			{
				orient = true;
				REPORT_STR("track orientation data");
			}
		}
		pack_data(buffer, dim, GetCommunicator());
		pack_data(buffer, orient, GetCommunicator());
		//pack nodes coords
		ENTER_BLOCK();
		{
			integer dim = GetDimensions();
			char have_gids = HaveGlobalID(NODE);
			std::vector<Storage::integer> global_ids;
			if( have_gids )
			{
				REPORT_STR("pack nodes global id");
				global_ids.resize(selems[0].size());
			}
			size_t marked_for_data = 0, marked_shared = 0, k = 0;
			std::vector<Storage::real> coords(selems[0].size()*dim);
			std::vector<char> flags(selems[0].size(), 0);
			ENTER_BLOCK();
			for(element_set::iterator it = selems[0].begin(); it != selems[0].end(); it++)
			{
				assert(!Hidden(*it));
				Storage::real_array c = RealArray(*it,tag_coords);
				//REPORT_VAL("packed coords",k << " " << c[0] << " " << c[1] << " " << c[2]);
				for(integer j = 0; j < dim; j++) coords[k*GetDimensions()+j] = c[j];
				Storage::integer & owner = IntegerDF(*it,tag_owner);
				bool flag_shared = false, flag_data = false;
#if defined(USE_PARALLEL_WRITE_TIME) && defined(DEBUG_UNPACK)
				std::stringstream pstr;
				Storage::integer_array proc = IntegerArrayDV(*it, tag_processors);
				pstr << " stat: " << Element::StatusName(GetStatus(*it)) << " owner: " << owner << " procs: ";
				for (Storage::integer_array::iterator ip = proc.begin(); ip != proc.end(); ++ip) pstr << *ip << " ";
#endif
				{
					if( owner == GetProcessorRank() )
					{
						assert(GetStatus(*it) != Element::Ghost);
						Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);
#if defined(USE_OMP)
#pragma omp critical
#endif
						{
							Storage::integer_array::iterator ip = std::lower_bound(proc.begin(), proc.end(), destination);
							if (ip == proc.end() || (*ip) != destination) proc.insert(ip, destination);
							if (GetStatus(*it) != Element::Shared)
							{
								//if( GetStatus(*it) != Element::Owned )
								//	std::cout << "node " << GetHandleID(*it) << " " << Element::StatusName(GetStatus(*it)) << std::endl;
								//assert(GetStatus(*it) == Element::Owned);
								flag_shared = true;
								++marked_shared;
								SetStatus(*it, Element::Shared);
							}
						}
					}
				}
				//TODO: 45
				if( owner != destination ) 
				{
					flags[k] |= 1;
					flag_data = true;
					SetPrivateMarker(*it,pack_tags_mrk);
					//TODO 46 old
					//	pack_tags[0].push_back(*it);
					++marked_for_data;
				}
				if( have_gids )
				{
					global_ids[k] = GlobalID(*it);
					//REPORT_VAL("packed global_id",global_ids.back());
				}
#if defined(USE_PARALLEL_WRITE_TIME) && defined(DEBUG_UNPACK)
				REPORT_STR("pack node " << k << " coords: " << c[0] << " " << c[1] << " " << c[2] 
						<< " flag: " << (flag_shared ? "shared" : "") << " " << (flag_data ? "data" : "")
						<< pstr.str());
#endif
				k++;
			}
			EXIT_BLOCK();
			pack_data(buffer,selems[0].size(),GetCommunicator());
			pack_data(buffer,marked_for_data,GetCommunicator());
			pack_data(buffer,have_gids,GetCommunicator());
			pack_data_vector(buffer,coords,GetCommunicator());
			pack_data_vector(buffer,global_ids,GetCommunicator());
			pack_data_vector(buffer, flags, GetCommunicator());
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[0].size());
			REPORT_VAL("total marked as shared", marked_shared << " / " << selems[0].size());
			REPORT_VAL("buffer position",buffer.size());
		}
		EXIT_BLOCK();
		//pack edges
		ENTER_BLOCK();
		{
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size(selems[1].size());
			std::vector<Storage::integer> low_conn_nums;
			std::vector<char> flags(selems[1].size(),0);
			low_conn_nums.reserve(selems[1].size()*2);
			size_t num = 0, k = 0, shifti = 0;
			size_t marked_for_data = 0, marked_shared = 0, packed_only_gid = 0;
			//pack edges
			ENTER_BLOCK();
			for(element_set::iterator it = selems[1].begin(); it != selems[1].end(); it++) 
			{
				assert(!Hidden(*it));
				bool pack_gid = false;
				shifti = low_conn_nums.size();
				//check that element is present on remote processor and pack global_id only
				if( pack_by_gids && HaveGlobalID(EDGE) )
				{
					Storage::integer_array procs = IntegerArray(*it,tag_processors);
					if( std::binary_search(procs.begin(),procs.end(),destination) )
						pack_gid = true;
				}
				if( pack_gid )
				{
					low_conn_size[k] = ENUMUNDEF;
					/*
					if( FindSharedGhost(EDGE,GlobalID(*it),destination) == InvalidHandle() )
					{
						Storage::integer_array procs = IntegerArray(*it,tag_processors);
						std::cout << "on " << GetProcessorRank() << " edge " << GlobalID(*it) << " procs ";
						for(unsigned q = 0; q < procs.size(); ++q) std::cout << procs[q] << " ";
						std::cout << Element::StatusName(GetStatus(*it));
						std::cout << " owner " << Integer(*it,tag_owner) << std::endl;
					}
					assert(FindSharedGhost(EDGE,GlobalID(*it),destination) != InvalidHandle());*/
					low_conn_nums.push_back(Integer(*it,tag_owner));
					low_conn_nums.push_back(GlobalID(*it));
					num+=2;
					packed_only_gid++;
				}
				else
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
						//assert( GetMarker(*jt,busy) );
						assert( HaveData(*jt,arr_position) );
						assert( (size_t)arr_position[*jt] < selems[0].size() );
						assert( *jt == selems[0][arr_position[*jt]] );
						low_conn_nums.push_back(arr_position[*jt]);
						low_conn_size[k]++;
						num++;
					}
					/*
					if( low_conn_size[k] != 2 )
					{
						std::cout << "pack edge " << *it << " lid " << ElementTypeName(GetHandleElementType(*it)) << ":" << GetHandleID(*it) << (Hidden(*it) ? " hidden " : " ok ") << std::endl;
						std::cout << "lc[" << lc.size() << "]: ";
						for(INMOST_DATA_ENUM_TYPE q = 0; q < lc.size(); ++q)
							std::cout << ElementTypeName(GetHandleElementType(lc[q])) << ":" << GetHandleID(lc[q]) << (Hidden(lc[q]) ? " hidden " : " ok ");
						Element::adj_type const & hc = HighConn(*it);
						std::cout << "hc[" << hc.size() << "]: ";
						for(INMOST_DATA_ENUM_TYPE q = 0; q < hc.size(); ++q)
							std::cout << ElementTypeName(GetHandleElementType(hc[q])) << ":" << GetHandleID(hc[q]) << (Hidden(hc[q]) ? " hidden " : " ok ");
					}
					assert(low_conn_size[k] == 2);
					*/
				}
				//REPORT_VAL("pack edge nodes ",k << " " << low_conn_nums[low_conn_nums.size()-2] << " " << low_conn_nums[low_conn_nums.size()-1]);
				bool flag_shared = false, flag_data = false;
				Storage::integer & owner = IntegerDF(*it,tag_owner);
#if defined(USE_PARALLEL_WRITE_TIME) && defined(DEBUG_UNPACK)
				std::stringstream pstr;
				Storage::integer_array proc = IntegerArrayDV(*it, tag_processors);
				pstr << " stat: " << Element::StatusName(GetStatus(*it)) << " owner: " << owner << " procs: ";
				for (Storage::integer_array::iterator ip = proc.begin(); ip != proc.end(); ++ip) pstr << *ip << " ";
#endif
				//
				{
					if( owner == GetProcessorRank() )
					{
						assert(GetStatus(*it) != Element::Ghost);
						Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);
#if defined(USE_OMP)
#pragma omp critical
#endif
						{
							Storage::integer_array::iterator ip = std::lower_bound(proc.begin(), proc.end(), destination);
							if (ip == proc.end() || (*ip) != destination) proc.insert(ip, destination);
							if (GetStatus(*it) != Element::Shared)
							{
								flag_shared = true;
								//REPORT_STR("edge " << k << " status " << Element::StatusName(GetStatus(*it)) << " now shared");
								//assert(GetStatus(*it) == Element::Owned);
								++marked_shared;
								SetStatus(*it, Element::Shared);
							}
						}
					}
				}
				if( owner != destination ) 
				{
					flag_data = true;
					flags[k] |= 1;
					SetPrivateMarker(*it,pack_tags_mrk);
					//TODO 46 old
					//	pack_tags[1].push_back(*it);
					++marked_for_data;
				}
#if defined(USE_PARALLEL_WRITE_TIME) && defined(DEBUG_UNPACK)
				REPORT_STR("pack edge " << k << " shift: " << shifti
					<< " nn: " << low_conn_nums[shifti] << " " << low_conn_nums[shifti+1]
					<< " flag: " << (flag_shared ? "shared" : "") << " " << (flag_data ? "data" : "") 
					<< pstr.str());
#endif
				k++;
			}
			EXIT_BLOCK();
			assert(low_conn_nums.size() == num);
			REPORT_VAL("number of edge nodes",num);
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[1].size());
			REPORT_VAL("total marked as shared", marked_shared << " / " << selems[1].size());
			REPORT_VAL("total packed by globalid ", packed_only_gid << " / " << selems[1].size());
			pack_data(buffer,selems[1].size(),GetCommunicator());
			pack_data(buffer,marked_for_data,GetCommunicator());
			pack_data_vector(buffer,low_conn_size,GetCommunicator());
			pack_data_vector(buffer,low_conn_nums,GetCommunicator());
			pack_data_vector(buffer,flags, GetCommunicator());
			REPORT_VAL("buffer position",buffer.size());
		}
		EXIT_BLOCK();
		//pack faces
		ENTER_BLOCK();
		{
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size(selems[2].size());
			std::vector<Storage::integer> low_conn_nums;
			std::vector<Storage::real> normals;
			std::vector<char> flags(selems[2].size(),0);
			low_conn_nums.reserve(selems[2].size()*4);
			if (orient) normals.reserve(selems[2].size() * 3);
			size_t num = 0, k = 0;
			size_t marked_for_data = 0, marked_shared = 0, packed_only_gid = 0;
			for(element_set::iterator it = selems[2].begin(); it != selems[2].end(); it++) 
			{
				assert(!Hidden(*it));
				//check that element is present on remote processor and pack global_id only
				bool pack_gid = false;
				if( pack_by_gids && HaveGlobalID(FACE) )
				{
					Storage::integer_array procs = IntegerArray(*it,tag_processors);
					if( std::binary_search(procs.begin(),procs.end(),destination) )
						pack_gid = true;
				}
				if( pack_gid )
				{
					low_conn_size[k] = ENUMUNDEF;
					low_conn_nums.push_back(Integer(*it,tag_owner));
					low_conn_nums.push_back(GlobalID(*it));
					num+=2;
					packed_only_gid++;
				}
				else
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
						//assert( GetMarker(*jt,busy) );
						assert( HaveData(*jt,arr_position) );
						assert( (size_t)arr_position[*jt] < selems[1].size() );
						assert( *jt == selems[1][arr_position[*jt]] );
						//REPORT_VAL("pack edge position",Integer(*jt,arr_position));
						low_conn_nums.push_back(arr_position[*jt]);
						low_conn_size[k]++;
						num++;
					}
				}
				//REPORT_VAL("pack size",low_conn_size[k]);
				Storage::integer & owner = IntegerDF(*it,tag_owner);
				{
					if( owner == GetProcessorRank() )
					{
						assert(GetStatus(*it) != Element::Ghost);
						Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);
#if defined(USE_OMP)
#pragma omp critical
#endif
						{
							Storage::integer_array::iterator ip = std::lower_bound(proc.begin(), proc.end(), destination);
							if (ip == proc.end() || (*ip) != destination) proc.insert(ip, destination);
							if (GetStatus(*it) != Element::Shared)
							{
								//assert(GetStatus(*it) == Element::Owned);
								++marked_shared;
								SetStatus(*it, Element::Shared);
							}
						}
					}
				}
				if( owner != destination ) 
				{
					flags[k] |= 1;
					SetPrivateMarker(*it, pack_tags_mrk);
					//TODO 46 old
					//	pack_tags[2].push_back(*it);
					++marked_for_data;

					if (orient)
					{
						Storage::real nrm[3];
						GetGeometricData(*it, NORMAL, nrm);
						normals.insert(normals.end(), nrm, nrm + 3);
					}
				}
				k++;
			}
			REPORT_VAL("number of faces",selems[2].size());
			REPORT_VAL("number of face edges",num);
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[2].size());
			REPORT_VAL("total marked as shared", marked_shared << " / " << selems[2].size());
			REPORT_VAL("total packed by globalid ", packed_only_gid << " / " << selems[2].size());
			pack_data(buffer,selems[2].size(),GetCommunicator());
			pack_data(buffer,marked_for_data,GetCommunicator());
			pack_data_vector(buffer,low_conn_size,GetCommunicator());
			pack_data_vector(buffer,low_conn_nums,GetCommunicator());
			pack_data_vector(buffer,flags,GetCommunicator());
			if (orient)
			{
				pack_data_vector(buffer, normals, GetCommunicator());
				REPORT_VAL("normals array size", normals.size());
			}
			REPORT_VAL("buffer position",buffer.size());
		}
		EXIT_BLOCK();
		//pack cells
		ENTER_BLOCK();
		{
			bool node_conns = HighConnTag().isDefined(CELL);
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size(selems[3].size());
			std::vector<INMOST_DATA_ENUM_TYPE> high_conn_size(node_conns ? selems[3].size() : 0);
			std::vector<Storage::integer> low_conn_nums;
			std::vector<Storage::integer> high_conn_nums;
			std::vector<char> flags(selems[3].size(),0);
			low_conn_nums.reserve(selems[3].size() * 6);
			if( node_conns ) high_conn_nums.reserve(selems[3].size() * 8);
			size_t num = 0, k = 0, num_high = 0;
			size_t marked_for_data = 0, marked_shared = 0, packed_only_gid = 0;
			for(element_set::iterator it = selems[3].begin(); it != selems[3].end(); it++) 
			{
				assert(!Hidden(*it));
				bool pack_gid = false;
				if( pack_by_gids && HaveGlobalID(CELL) )
				{
					Storage::integer_array procs = IntegerArray(*it,tag_processors);
					if( std::binary_search(procs.begin(),procs.end(),destination) )
						pack_gid = true;
				}
				if( pack_gid )
				{
					low_conn_size[k] = ENUMUNDEF;
					if( node_conns ) high_conn_size[k] = 0;
					low_conn_nums.push_back(Integer(*it,tag_owner));
					low_conn_nums.push_back(GlobalID(*it));
					num+=2;
					packed_only_gid++;
				}
				else
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
						//assert( GetMarker(*jt,busy) );
						assert( HaveData(*jt,arr_position) );
						assert( (size_t)arr_position[*jt] < selems[2].size() );
						assert( *jt == selems[2][arr_position[*jt]] );
						//REPORT_VAL("pack face position",Integer(*jt,arr_position));
						low_conn_nums.push_back(arr_position[*jt]);
						low_conn_size[k]++;
						num++;
					}
					if (node_conns)
					{
						high_conn_size[k] = 0;
						Element::adj_type const& hc = HighConn(*it);
						for (Element::adj_type::const_iterator jt = hc.begin(); jt != hc.end(); jt++) if (!Hidden(*jt))
						{
							//TODO 44 old
							//element_set::iterator find = std::lower_bound(selems[0].begin(),selems[0].end(),(*jt));
							//assert( find == selems[0].end());
							//assert((*find) == (*jt) );
							//high_conn_nums.push_back(static_cast<Storage::integer>(find - selems[0].begin()));
							//assert( GetMarker(*jt,busy) );
							assert(HaveData(*jt, arr_position));
							assert((size_t)arr_position[*jt] < selems[0].size());
							assert(*jt == selems[0][arr_position[*jt]]);
							//REPORT_VAL("pack node position",Integer(*jt,arr_position));
							high_conn_nums.push_back(arr_position[*jt]);
							high_conn_size[k]++;
							num_high++;
						}
					}
				}
				Storage::integer & owner = IntegerDF(*it,tag_owner);
				{
					if( owner == GetProcessorRank() )
					{
						assert(GetStatus(*it) != Element::Ghost);
						Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);
#if defined(USE_OMP)
#pragma omp critical
#endif
						{
							Storage::integer_array::iterator ip = std::lower_bound(proc.begin(), proc.end(), destination);
							if (ip == proc.end() || (*ip) != destination) proc.insert(ip, destination);
							if (GetStatus(*it) != Element::Shared)
							{
								//assert(GetStatus(*it) == Element::Owned);
								++marked_shared;
								SetStatus(*it, Element::Shared);
							}
						}
					}
				}
				if( owner != destination ) 
				{
					flags[k] |= 1;
					SetPrivateMarker(*it, pack_tags_mrk);
					//TODO 46 old
					//	pack_tags[3].push_back(*it);
					++marked_for_data;
				}
				k++;
			}
			REPORT_VAL("number of cell faces",num);
			REPORT_VAL("number of cell nodes", num_high);
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[3].size());
			REPORT_VAL("total marked as shared", marked_shared << " / " << selems[3].size());
			REPORT_VAL("total packed by globalid ", packed_only_gid << " / " << selems[3].size());
			pack_data(buffer,selems[3].size(),GetCommunicator());
			pack_data(buffer,marked_for_data,GetCommunicator());
			pack_data_vector(buffer,low_conn_size,GetCommunicator());
			pack_data_vector(buffer,low_conn_nums,GetCommunicator());
			if (node_conns)
			{
				pack_data_vector(buffer, high_conn_size, GetCommunicator());
				pack_data_vector(buffer, high_conn_nums, GetCommunicator());
			}
			pack_data_vector(buffer, flags, GetCommunicator());
			REPORT_VAL("buffer position",buffer.size());
		}
		EXIT_BLOCK();
        /////////////////////////////////////////
        // pack esets
		ENTER_BLOCK();
		{
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size(selems[4].size()); //store number of elements in set
			std::vector<INMOST_DATA_ENUM_TYPE> high_conn_size(selems[4].size()); // store number of children in set
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_nums;  // array composed elements : ElementType and position in array 
			std::vector<Storage::integer> high_conn_nums;  // array of indexes of parent (-1 if no parent), then all children
			low_conn_nums.reserve(selems[4].size()*8);
			high_conn_nums.reserve(selems[4].size()*3);
			std::vector<char> flags(selems[4].size(), 0);
			std::vector<char> names_buff;
			size_t k = 0;
			size_t marked_for_data = 0, marked_shared = 0;//, packed_only_name = 0;
			//int names_buff_size = 0;
            // Pack sequence: 
            // 1) all names of sets
            // 2) low conn for sets (arr_pos + element_type)
            // 3) high conn for sets (3 arr_pos: child, sibling, parent.  -1 if has't)
            
            // Compute names_buff_size
			//for(element_set::iterator it = selems[4].begin(); it != selems[4].end(); it++) names_buff_size += ElementSet(this,*it).GetName().size() + 1;
			//names_buff.resize(names_buff_size);
			//int names_buff_pos = 0;
			//MarkerType pack_set = CreateMarker();
			//if( selems[4].size() ) SetMarkerArray(&selems[4][0],selems[4].size(),pack_set);
			for(element_set::iterator it = selems[4].begin(); it != selems[4].end(); it++) 
			{
				assert(!Hidden(*it));
				ElementSet set(this, *it);
				//bool print = (set.GetName() == "AM_R446C0") || (set.GetName() == "AM_R446");
				//REPORT_VAL("pack set ", set.GetName());
				//REPORT_VAL("set owner ", set.Integer(tag_owner));
				Storage::integer & owner = Integer(*it,tag_owner);
				//if( GetStatus(*it) == Element::Owned || GetStatus(*it) == Element::Any )
				//	owner = mpirank;
				REPORT_STR("set name \"" << set.GetName() << "\" owner " << owner);
				
				if( owner == GetProcessorRank() )
				{
					assert(GetStatus(*it) != Element::Ghost);
					Storage::integer_array proc = IntegerArray(*it, tag_processors);
#if defined(USE_OMP)
#pragma omp critical
#endif
					{
						Storage::integer_array::iterator ip = std::lower_bound(proc.begin(), proc.end(), destination);
						if (ip == proc.end() || (*ip) != destination) proc.insert(ip, destination);
						if (GetStatus(*it) != Element::Shared)
						{
							//assert(GetStatus(*it) == Element::Owned);
							++marked_shared;
							SetStatus(*it, Element::Shared);
						}
					}
				}

				
				
				if( owner != destination )
				{
					flags[k] |= 1;
					SetPrivateMarker(*it, pack_tags_mrk);
					//TODO 46 old
					//	pack_tags[3].push_back(*it);
					++marked_for_data;
				}
				/*
				if( (set.GetName() == "AM_ROOT_SET" || set.GetName() == "AM_R270") && GetProcessorRank() == 4 && destination == 0 )
				{
					std::cout << "Hello there!!!" << std::endl;
					Storage::integer_array arr = set.IntegerArray(tag_processors);
					std::cout << set.GetName() << " offset " << k;
					std::cout << " procs ";
					for(int r = 0; r < arr.size(); ++r) std::cout << arr[r] << " ";
					//std::cout << std::endl;
					std::cout << "owner " << owner << (set.GetMarker(pack_tags_mrk) ? " has marker " : " no marker ") << std::endl;
				}
				*/
				// Add set name to names_buffer
				std::string set_name = set.GetName();
				names_buff.insert(names_buff.end(),set_name.begin(),set_name.end());
				names_buff.push_back('\0');
				
				/*
				bool pack_data = true;
				{
					Storage::integer_array procs = IntegerArray(*it,tag_processors);
					if( std::binary_search(procs.begin(),procs.end(),destination) )
						pack_data = false;
				}
				if( !pack_data )
				{
					low_conn_size[k] = 0;
					high_conn_size[k] = 1; //reserve for parent
					high_conn_nums.push_back(-1);
					packed_only_name++;
				}
				else*/
				{
					//strcpy(&names_buff[names_buff_pos],set.GetName().c_str());
					//names_buff[names_buff_pos+set.GetName().length()] = '\0';
					//names_buff_pos += set.GetName().length() + 1;
					// Add all low conns to low_conn_nums
					// std::cout << "From " << GetProcessorRank() << " to " << destination;
					// std::cout << " for set " << set.GetName() << " send elements ";
					low_conn_size[k] = 0;
					Element::adj_type const & lc = LowConn(*it);
					for(Element::adj_type::const_iterator jt = lc.begin(); jt != lc.end(); jt++)
					{
						if( *jt != InvalidHandle() && !Hidden(*jt) && HaveData(*jt,arr_position) )
						{
							assert(HaveData(*jt,arr_position));
							assert((size_t)arr_position[*jt] < selems[GetHandleElementNum(*jt)].size() );
							assert(*jt == selems[GetHandleElementNum(*jt)][arr_position[*jt]]);
							// std::cout << " " << ElementTypeName(GetHandleElementType(*jt));
							// std::cout << ":" << GetHandleID(*jt);
							low_conn_nums.push_back(ComposeHandle(GetHandleElementType(*jt),arr_position[*jt]));
							low_conn_size[k]++;
						}
					}
					// std::cout << " total " << low_conn_size[k] << " elements " <<  std::endl;
					//if( print ) std::cout << GetProcessorRank() << " set " << set.GetName();
					high_conn_size[k] = 1; //reserve for parent
					if( set.HaveParent() && set.GetParent().HaveData(arr_position) )
					{
						assert(HaveData(set.GetParent().GetHandle(),arr_position));
						assert( (size_t)arr_position[set.GetParent()] < selems[4].size() );
						assert(set.GetParent().GetHandle() == selems[4][arr_position[set.GetParent()]]);
						high_conn_nums.push_back(arr_position[set.GetParent()]);
						//if( print ) std::cout << " pack parent " << set.GetParent().GetName();
					}
					else
					{
						high_conn_nums.push_back(-1);
						//if( print ) std::cout << " no pack parent " << (set.HaveChild() ? set.GetChild().GetName() : "none");
					}
					if( set.HaveChild() )
					{
						//if( print ) std::cout << " pack";
						for(ElementSet jt = set.GetChild(); jt.isValid(); jt = jt.GetSibling() )
						{
							if( jt.HaveData(arr_position) )
							{
								assert(HaveData(jt.GetHandle(),arr_position));
								assert( (size_t)arr_position[jt] < selems[4].size() );
								assert(jt.GetHandle() == selems[4][arr_position[jt]]);
								high_conn_nums.push_back(arr_position[jt]);
								high_conn_size[k]++;
								//if( print ) std::cout << " " << jt.GetName();
							}
							//else if( print ) std::cout << " no pack " << jt.GetName();
						}
					}
				}
				//else  if( print ) std::cout << " no children";
				//if( print ) std::cout << " for " << destination << std::endl;
				
				/*	
				if( set.HaveSibling() && set.GetSibling().GetMarker(busy) )
                {
					assert(set.GetSibling().GetHandle() == selems[4][Integer(set.GetSibling().GetHandle(),arr_position)]);
					high_conn_nums[k*3 + 2] = Integer(set.GetSibling().GetHandle(),arr_position);
				}
				else 
					high_conn_nums[k*3 + 2] = -1;
					*/
				k++;
			}

			//std::cout << "number of sets: " << selems[4].size() << std::endl;
			//if( selems[4].size() ) RemMarkerArray(&selems[4][0],selems[4].size(),pack_set);
			//ReleaseMarker(pack_set);
			pack_data(buffer,selems[4].size(),GetCommunicator());
			pack_data(buffer,marked_for_data,GetCommunicator());
			pack_data_vector(buffer,names_buff,GetCommunicator());
			pack_data_vector(buffer,low_conn_size,GetCommunicator());
			pack_data_vector(buffer,low_conn_nums,GetCommunicator());
			pack_data_vector(buffer,high_conn_size,GetCommunicator());
			pack_data_vector(buffer,high_conn_nums,GetCommunicator());
			pack_data_vector(buffer, flags, GetCommunicator());
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[4].size());
			REPORT_VAL("total marked as shared", marked_shared << " / " << selems[4].size());
			//REPORT_VAL("total packed only name ", packed_only_name << " / " << selems[4].size());
			REPORT_VAL("buffer position",buffer.size());
		}
		EXIT_BLOCK();
		//for(int etypenum = ElementNum(NODE); etypenum <= ElementNum(ESET); ++etypenum) if( !selems[etypenum].empty() )
		//	RemMarkerArray(&selems[etypenum][0],selems[etypenum].size(),busy);
			//for(element_set::iterator it = selems[etypenum].begin(); it != selems[etypenum].end(); ++it) RemMarker(*it,busy);
		
		//ReleaseMarker(busy);
		
        /////////////////////////////////////////
		SetPrivateMarker(GetHandle(), pack_tags_mrk);
		selems[5].push_back(GetHandle()); //data on mesh
		
		//DeleteTag(arr_position);
		//for(int i = 3; i >= 0; --i)
		//	all.insert(all.end(),selems[i].begin(),selems[i].end());
		//REPORT_VAL("number of all elements",all.size());
		//Storage::integer_array proc = IntegerArrayDV(GetHandle(),tag_processors);
		//Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),destination);
		//if( ip == proc.end() || (*ip) != destination ) proc.insert(ip,destination);		
		//pack number of tags
		ENTER_BLOCK();
		REPORT_STR("prepare elements to pack tag data");
		{
			REPORT_VAL("number of tags",tag_list.size());
			pack_data(buffer,tag_list.size(),GetCommunicator());
			REPORT_VAL("buffer position",buffer.size());
		}
		for(INMOST_DATA_ENUM_TYPE i = 0; i < tag_list.size(); i++)
		{
			//assert(HaveTag(tag_list[i]));
			const Tag & tag = tag_list[i];
			//pack tag name
			ENTER_BLOCK();
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
				//REPORT_VAL("characters in name",name.size());
				REPORT_VAL("name of tag",name);
				REPORT_VAL("data type",DataTypeName(static_cast<DataType>(datatype)));
				REPORT_VAL("data size",datasize);
//#if defined(USE_PARALLEL_WRITE_TIME)
//				for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype) )
//				{
//					if( defined & etype ) REPORT_VAL("defined on",ElementTypeName(etype));
//					if( sparse & etype ) REPORT_VAL("sparse on",ElementTypeName(etype));
//				}
//#endif
				pack_data(buffer,datatype,GetCommunicator());
				pack_data(buffer,defined,GetCommunicator());
				pack_data(buffer,sparse,GetCommunicator());
				pack_data(buffer,datasize,GetCommunicator());
				pack_data(buffer,name.size(),GetCommunicator());
				pack_data_array(buffer,name.c_str(),name.size(),GetCommunicator());
				REPORT_VAL("buffer position",buffer.size());
			}
			EXIT_BLOCK();
			//TODO 46 old
			//PackTagData(GetTag(tag_list[i]),pack_tags,NODE | EDGE | FACE | CELL | ESET,0,buffer);
			PackTagData(tag,selems,destination,NODE | EDGE | FACE | CELL | ESET | MESH,pack_tags_mrk,buffer,arr_position);
			//PackTagData(tag,selems,NODE | EDGE | FACE | CELL | ESET,0,buffer);
			//std::cout << mpirank << " After pack_tag_data\n" << std::endl;
			REPORT_VAL("buffer position",buffer.size());
		}
		EXIT_BLOCK();
		
		for(integer i = ElementNum(NODE); i <= ElementNum(MESH); i++) if( !selems[i].empty() )
			RemPrivateMarkerArray(&selems[i][0],static_cast<enumerator>(selems[i].size()),pack_tags_mrk);
		ReleasePrivateMarker(pack_tags_mrk);

#else
		(void) selems;
		(void) buffer;
		(void) destination;
		(void) tag_list;
#endif
		REPORT_VAL("destination",destination);
		EXIT_FUNC();
	}
	
	
	void Mesh::UnpackElementsData(elements_by_type & selems, buffer_type & buffer, int source, size_t & buffer_position, tag_set & tag_recv)
	{
		ENTER_FUNC();
		REPORT_VAL("source",source);
#if defined(USE_MPI)
		MarkerType unpack_tags_mrk = CreateMarker();
		MarkerType orient = 0;
		integer remote_dim;
		char orient_flag = false;
		unpack_data(buffer, buffer_position, remote_dim, GetCommunicator());
		unpack_data(buffer,buffer_position,orient_flag, GetCommunicator());
		if (remote_dim != dim) //change number of dimensions
		{
			if( NumberOfNodes() == 0)
				dim = remote_dim;
			else throw Impossible; // Remote has different number of dimensions and local has nodes!!!
		}
		if (orient_flag)
		{
			orient = CreateMarker();
			REPORT_STR("track orientation data");
		}
		//TODO 46 old
		//elements_by_type unpack_tags;
		std::vector<HandleType> old_nodes(NumberOfNodes());
		//all.clear();
		double time = Timer();
		ENTER_BLOCK();
		REPORT_VAL("buffer size",buffer.size());
		REPORT_VAL("buffer position",buffer_position);
		for(int k = ElementNum(NODE); k <= ElementNum(ESET); ++k) selems[k].clear();
		//TODO 49
		{
			int k = 0;
			for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); it++)
				old_nodes[k++] = *it;
		}
		if( !old_nodes.empty() )
		{
			bool have_node_gid = HaveGlobalID(NODE);
			if( have_node_gid )
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
#if defined(USE_PARALLEL_WRITE_TIME) && defined(DEBUG_UNPACK)
			ENTER_BLOCK();
			for (std::vector<HandleType>::iterator it = old_nodes.begin(); it != old_nodes.end(); ++it)
			{
				real_array c = RealArray(*it, tag_coords);
				integer_array proc = IntegerArray(*it, tag_processors);
				std::stringstream pstr;
				pstr << " stat: " << Element::StatusName(GetStatus(*it)) << " owner: " << IntegerDF(*it, tag_owner) << " procs: ";
				for (integer_array::iterator ip = proc.begin(); ip != proc.end(); ++ip) pstr << *ip << " ";
				REPORT_STR("node: " << *it << " coords " << c[0] << " " << c[1] << " " << c[2] << " gid " << (have_node_gid ? GlobalID(*it) : -1) << pstr.str());
			}
			EXIT_BLOCK();
#endif
		}
		time = Timer() - time;
		REPORT_STR("gather and sort old nodes");
		REPORT_VAL("time", time);	
		EXIT_BLOCK();
		
		ENTER_BLOCK();
		time = Timer();
		//unpack nodes
		{
			char have_gids = 0;
			std::vector<Storage::real> coords;
			std::vector<Storage::integer> global_ids;
			std::vector<char> flags;
			integer dim = GetDimensions();
			REPORT_STR("unpack number of nodes");
            size_t num = 0, marked_remote = 0;
            unpack_data(buffer,buffer_position,num,GetCommunicator());
            unpack_data(buffer,buffer_position,marked_remote,GetCommunicator());
            unpack_data(buffer,buffer_position,have_gids,GetCommunicator());
            REPORT_VAL("number of nodes",num);
            REPORT_STR("unpack nodes coordinates");
            unpack_data_vector(buffer,buffer_position,coords,GetCommunicator());
            REPORT_STR("unpack nodes global identificators");
            unpack_data_vector(buffer,buffer_position,global_ids,GetCommunicator());
			REPORT_STR("unpack flags");
			unpack_data_vector(buffer, buffer_position, flags, GetCommunicator());
			REPORT_VAL("buffer position",buffer_position);
			selems[0].resize(num);
			//TODO 46 old
			//unpack_tags[0].reserve(num);
			REPORT_STR("creating nodes");
			size_t found = 0, marked_for_data = 0, marked_ghost = 0;
			ENTER_BLOCK();
			for(size_t i = 0; i < num; i++)
			{
				HandleType new_node = InvalidHandle();
				size_t find = static_cast<size_t>(-1);
				bool flag_found = false, flag_data = false;
				//TODO 49
				if( !old_nodes.empty() )
				{
					if( have_gids && HaveGlobalID(NODE) )
					{
						//REPORT_STR("searching node by global id");
						//REPORT_VAL("global_id",global_ids[i]);
						std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),global_ids[i],GlobalIDComparator(this));
						if( it != old_nodes.end() ) 
						{
							//REPORT_VAL("lower bound found",GlobalID(*it));
							if( global_ids[i] == GlobalID(*it) )
							{
								found = true;
								find = it - old_nodes.begin();
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
								found = true;
								find = it - old_nodes.begin();
								found++;
							}
						} //else REPORT_STR("lower bound points to end");
						//REPORT_VAL("found",find);
					}
				}
				if( find != static_cast<size_t>(-1) ) 
				{
					new_node = old_nodes[find];
					if( IntegerDF(new_node,tag_owner) != GetProcessorRank() ) 
					{
						flag_data = true;
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
						flag_data = true;
						SetMarker(new_node,unpack_tags_mrk);
						++marked_for_data;
						++marked_ghost;
						//std::cout << "New ghost node at " << GetProcessorRank() << " arrived from " << source << std::endl;
					}
				}
				if ((flags[i] & 1 ? true : false) != GetMarker(new_node, unpack_tags_mrk))
				{
					std::cout << __FILE__ << ":" << __LINE__ << " on " << GetProcessorRank() << " remote has " << (flags[i] & 1 ? "data" : "no data")
						<< " for node " << i << " and local expects " << (GetMarker(new_node, unpack_tags_mrk) ? "data" : "no data") << std::endl;
					REPORT_STR(__FILE__ << ":" << __LINE__ << " on " << GetProcessorRank() << " remote has " << (flags[i] & 1 ? "data" : "no data")
						<< " for node " << i << " and local expects " << (GetMarker(new_node, unpack_tags_mrk) ? "data" : "no data"));
				}
				assert(new_node != InvalidHandle());
				selems[0][i] = new_node;
#if defined(USE_PARALLEL_WRITE_TIME) && defined(DEBUG_UNPACK)
				std::stringstream pstr;
				Storage::integer_array proc = IntegerArrayDV(new_node, tag_processors);
				pstr << " stat: " << Element::StatusName(GetStatus(new_node)) << " owner: " << IntegerDF(new_node, tag_owner) << " procs: ";
				for (Storage::integer_array::iterator ip = proc.begin(); ip != proc.end(); ++ip) pstr << *ip << " ";
				REPORT_STR("unpack node" << i << " coords in " 
					<< coords[i * dim + 0] << " " << coords[i * dim + 1] << " " << coords[i * dim + 2]
					<< " got " << RealArray(new_node,tag_coords)[0] << " " << RealArray(new_node, tag_coords)[1] << " " << RealArray(new_node, tag_coords)[2]
					<< " gid in " << (have_gids ? global_ids[i] : -1) << " got " << (have_gids && HaveGlobalID(NODE) ? GlobalID(new_node) : -1)
					<< " flag: " << (flag_found ? "found" : "created") << " " << (flag_data ? "data" : "") 
					<< pstr.str());
#endif
			}
			EXIT_BLOCK();
			REPORT_VAL("total found", found << " / " << selems[0].size());
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[0].size());
			REPORT_VAL("total marked ghost", marked_ghost << " / " << selems[0].size());
			REPORT_VAL("buffer position",buffer_position);
			if( marked_for_data != marked_remote )
			{
				std::cout << __FILE__ << ":" << __LINE__ << " marked for data " << marked_for_data << " remote " << marked_remote;
				std::cout << " on " << GetProcessorRank() << " from " << source;
				std::cout << " found " << found << " marked ghost " << marked_ghost << " total " << selems[0].size();
				std::cout << std::endl;
			}
			assert(marked_for_data == marked_remote);
		}
		time = Timer() - time;
		REPORT_STR("unpack nodes");
		REPORT_VAL("time", time);	
		EXIT_BLOCK();
		
		ENTER_BLOCK();
		time = Timer();
		//unpack edges
		{
			ElementArray<Node> e_nodes(this,2);
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size;
			std::vector<Storage::integer> low_conn_nums;
			std::vector<char> flags;
			size_t num = 0, marked_remote = 0;
			unpack_data(buffer,buffer_position,num,GetCommunicator());
			unpack_data(buffer,buffer_position,marked_remote,GetCommunicator());
			REPORT_VAL("number of edges",num);
			unpack_data_vector(buffer,buffer_position,low_conn_size,GetCommunicator());
			unpack_data_vector(buffer,buffer_position,low_conn_nums,GetCommunicator());
			unpack_data_vector(buffer, buffer_position, flags, GetCommunicator());
			REPORT_VAL("buffer position",buffer_position);
			selems[1].reserve(num);
			
			size_t found = 0, marked_for_data = 0, marked_ghost = 0, unpacked_by_gid = 0, shift = 0, shifti = 0;
			ENTER_BLOCK();
			for(size_t i = 0; i < num; i++)
			{
				HandleType new_edge = InvalidHandle();
				bool flag_found = false, flag_data = false;
				shifti = shift;
#if defined(USE_PARALLEL_WRITE_TIME) && defined(DEBUG_UNPACK)
				std::stringstream istr;
#endif
				if( low_conn_size[i] == ENUMUNDEF )
				{
					if( !HaveGlobalID(EDGE) ) std::cout << GetProcessorRank() << " should find edge by global id but it is not defined " << std::endl;
					assert(HaveGlobalID(EDGE));
					Storage::integer owner = low_conn_nums[shift++];
					Storage::integer gid = low_conn_nums[shift++];
					
					// I have this element and remote processor also has this element, so it should be either in ghost or shared elements
					new_edge = FindSharedGhost(EDGE,gid,source,owner);
					
					if( new_edge == InvalidHandle() ) std::cout << GetProcessorRank() << " edge not found by global id " << gid << std::endl;
					assert(new_edge != InvalidHandle());
					unpacked_by_gid++;
				}
				else
				{
					//~ assert(low_conn_size[i] == 2);
					e_nodes.resize(low_conn_size[i]);
					//REPORT_VAL("Unpack size",low_conn_size[i]);
					for(INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
					{
						assert((size_t)low_conn_nums[shift+j] < selems[0].size());
						//REPORT_VAL("Unpack node position",low_conn_nums[shift+j]);
						e_nodes.at(j) = selems[0][low_conn_nums[shift+j]];
#if defined(USE_PARALLEL_WRITE_TIME) && defined(DEBUG_UNPACK)
						istr << " node ( " << low_conn_nums[shift + j] << " id " << e_nodes[j].LocalID() << " "
							<< e_nodes[j].Coords()[0] << " " << e_nodes[j].Coords()[1] << " " << e_nodes[j].Coords()[2]
							<< std::endl << "\tedges: ";
						ElementArray<Edge> edges = e_nodes[j].getEdges();
						for (ElementArray<Edge>::iterator kt = edges.begin(); kt != edges.end(); ++kt)
							istr << "\t\t" << kt->LocalID() << " " << kt->getBeg()->LocalID() << "<->" << kt->getEnd()->LocalID() << " stat " << Element::StatusName(kt->GetStatus()) << std::endl;
						istr << ")";
#endif
					}
					shift+= low_conn_size[i];
					if( !e_nodes.empty() )
						new_edge = FindSharedAdjacency(e_nodes.data(),static_cast<enumerator>(e_nodes.size()));
				}
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
					flag_data = true;
					SetMarker(new_edge,unpack_tags_mrk);
					++marked_for_data;
					++marked_ghost;
					
					assert(IntegerArrayDV(new_edge,tag_processors).size() == 0);
					//std::cout << "New ghost edge at " << GetProcessorRank() << " arrived from " << source << std::endl;
				}
				else
				{
					flag_found = true;
					if( IntegerDF(new_edge,tag_owner) != GetProcessorRank() )
					{
						//TODO 46 old
						//unpack_tags[1].push_back(new_edge);
						flag_data = true;
						SetMarker(new_edge,unpack_tags_mrk);
						++marked_for_data;
					}
					++found;
				}
				if ((flags[i] & 1 ? true : false) != GetMarker(new_edge, unpack_tags_mrk))
				{
					std::cout << __FILE__ << ":" << __LINE__ << " on " << GetProcessorRank() << " remote has " << (flags[i] & 1 ? "data" : "no data")
						<< " for edge " << i << " and local expects " << (GetMarker(new_edge, unpack_tags_mrk) ? "data" : "no data") << std::endl;
					REPORT_STR(__FILE__ << ":" << __LINE__ << " on " << GetProcessorRank() << " remote has " << (flags[i] & 1 ? "data" : "no data")
						<< " for edge " << i << " and local expects " << (GetMarker(new_edge, unpack_tags_mrk) ? "data" : "no data"));
				}
#if defined(USE_PARALLEL_WRITE_TIME) && defined(DEBUG_UNPACK)
				std::stringstream pstr;
				Storage::integer_array proc = IntegerArrayDV(new_edge, tag_processors);
				pstr << " stat: " << Element::StatusName(GetStatus(new_edge)) << " owner: " << IntegerDF(new_edge, tag_owner) << " procs: ";
				for (Storage::integer_array::iterator ip = proc.begin(); ip != proc.end(); ++ip) pstr << *ip << " ";
				REPORT_STR("unpack edge: " << i << " shift: " << shifti
					<< " nn: " << low_conn_nums[shifti] << " " << low_conn_nums[shifti + 1]
					<< " flag: " << (flag_found ? "found" : "created") << " " << (flag_data ? "data" : "") 
					<< pstr.str() << " info " << istr.str());
#endif
				selems[1].push_back(new_edge);
			}
			EXIT_BLOCK();
			REPORT_VAL("total found", found << " / " << selems[1].size());
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[1].size());
			REPORT_VAL("total marked ghost", marked_ghost << " / " << selems[1].size());
			REPORT_VAL("total unpacked by gid", unpacked_by_gid << " / " << selems[1].size());
			REPORT_VAL("buffer position",buffer_position);
			if( marked_for_data != marked_remote )
			{
				std::cout << __FILE__ << ":" << __LINE__ << " marked for data " << marked_for_data << " remote " << marked_remote;
				std::cout << " on " << GetProcessorRank() << " from " << source << std::endl;
			}
			assert(marked_for_data == marked_remote);
		}
		time = Timer() - time;
		REPORT_STR("unpack edges");
		REPORT_VAL("time", time);
		EXIT_BLOCK();
		
		ENTER_BLOCK();
		time = Timer();
		//unpack faces
		{
			ElementArray<Edge> f_edges(this);
			//~ shift = 0;
			std::vector<Storage::real> normals;
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size;
			std::vector<Storage::integer> low_conn_nums;
			std::vector<char> flags;
			size_t num = 0, marked_remote = 0;
			unpack_data(buffer,buffer_position,num,GetCommunicator());
			unpack_data(buffer,buffer_position,marked_remote,GetCommunicator());
			REPORT_VAL("number of faces",num);
			unpack_data_vector(buffer,buffer_position,low_conn_size,GetCommunicator());
			unpack_data_vector(buffer,buffer_position,low_conn_nums,GetCommunicator());
			unpack_data_vector(buffer,buffer_position,flags,GetCommunicator());
			if (orient)
			{
				unpack_data_vector(buffer, buffer_position, normals, GetCommunicator());
				REPORT_VAL("normals array size", normals.size());
			}
			REPORT_VAL("buffer position",buffer_position);
			selems[2].reserve(num);
			
			size_t found = 0, marked_for_data = 0, marked_ghost = 0, shift = 0, ncnt = 0;
			for(size_t i = 0; i < num; i++)
			{
				HandleType new_face = InvalidHandle();
				if( low_conn_size[i] == ENUMUNDEF )
				{
					if( !HaveGlobalID(FACE) ) std::cout << GetProcessorRank() << " should find face by global id but it is not defined " << std::endl;
					assert(HaveGlobalID(FACE));
					Storage::integer owner = low_conn_nums[shift++];
					Storage::integer gid = low_conn_nums[shift++];
					
					// I have this element and remote processor also has this element, so it should be either in ghost or shared elements
					new_face = FindSharedGhost(FACE,gid,source,owner);
					if( new_face == InvalidHandle() ) std::cout << GetProcessorRank() << " face not found by global id " << gid << std::endl;
					assert(new_face != InvalidHandle());
				}
				else
				{
					f_edges.resize(low_conn_size[i]);
					//REPORT_VAL("Unpack size",low_conn_size[i]);
					for(INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
					{
						//REPORT_VAL("Unpack edge position",low_conn_nums[shift+j]);
						assert((size_t)low_conn_nums[shift+j] < selems[1].size());
						f_edges.at(j) = selems[1][low_conn_nums[shift+j]];
					}
					shift+= low_conn_size[i];
					if( !f_edges.empty() )
						new_face = FindSharedAdjacency(f_edges.data(),static_cast<enumerator>(f_edges.size()));
				}
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

					//std::cout << "New ghost face at " << GetProcessorRank() << " arrived from " << source << std::endl;
				} 
				else 
				{
					//current solution - no communication of geometric data
					//here we should check that remote face orientation match
					//existing face orientation, otherwise we would need
					//to swap normal if it is precomputed and comes from remote
					//processor the wrong way
					if( IntegerDF(new_face,tag_owner) != GetProcessorRank() ) 
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
				//if we unpack data for this face,
				//and there are orientable tags,
				//then compare local and remote normals and
				//later change the sign accordingly
				if (orient && GetMarker(new_face,unpack_tags_mrk))
				{
					Storage::real lnrm[3], * rnrm = &normals[ncnt * 3];
					GetGeometricData(new_face, NORMAL, lnrm);
					if (lnrm[0] * rnrm[0] + lnrm[1] * rnrm[1] + lnrm[2] * rnrm[2] < 0.0)
						SetMarker(new_face, orient);
					ncnt++;
				}

				// if( !Face(this,new_face).CheckNormalOrientation() )
				// {
				// 	std::cout << __FILE__ << ":" << __LINE__ << " rank " << GetProcessorRank() << " face " << new_face << " " << (Face(this,new_face).Boundary()?"bnd":"int") << " " << Element::StatusName(Face(this,new_face).GetStatus()) << std::endl;
				// }
				if ((flags[i] & 1 ? true : false) != GetMarker(new_face, unpack_tags_mrk))
				{
					std::cout << __FILE__ << ":" << __LINE__ << " on " << GetProcessorRank() << " remote has " << (flags[i] & 1 ? "data" : "no data")
						<< " for face " << i << " and local expects " << (GetMarker(new_face, unpack_tags_mrk) ? "data" : "no data") << std::endl;
					REPORT_STR(__FILE__ << ":" << __LINE__ << " on " << GetProcessorRank() << " remote has " << (flags[i] & 1 ? "data" : "no data")
						<< " for face " << i << " and local expects " << (GetMarker(new_face, unpack_tags_mrk) ? "data" : "no data"));
				}

				selems[2].push_back(new_face);
			}
			REPORT_VAL("total found", found << " / " << selems[2].size());
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[2].size());
			REPORT_VAL("total marked ghost", marked_ghost << " / " << selems[2].size());
			REPORT_VAL("buffer position",buffer_position);
			if( marked_for_data != marked_remote )
			{
				std::cout << __FILE__ << ":" << __LINE__ << " marked for data " << marked_for_data << " remote " << marked_remote;
				std::cout << " on " << GetProcessorRank() << " from " << source << std::endl;
			}
			assert(marked_for_data == marked_remote);
		}
		time = Timer() - time;
		REPORT_STR("unpack faces");
		REPORT_VAL("time", time);
		EXIT_BLOCK();
		
		ENTER_BLOCK();
		time = Timer();
		//unpack cells
		{
			bool node_conns = HighConnTag().isDefined(CELL);
			ElementArray<Face> c_faces(this);
			ElementArray<Node> c_nodes(this);
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size;
			std::vector<INMOST_DATA_ENUM_TYPE> high_conn_size;
			std::vector<Storage::integer> low_conn_nums;
			std::vector<Storage::integer> high_conn_nums;
			std::vector<char> flags;
			size_t num = 0, marked_remote = 0;
			unpack_data(buffer,buffer_position,num,GetCommunicator());
			unpack_data(buffer,buffer_position,marked_remote,GetCommunicator());
			REPORT_VAL("number of cells",num);
			unpack_data_vector(buffer,buffer_position,low_conn_size,GetCommunicator());
			unpack_data_vector(buffer,buffer_position,low_conn_nums,GetCommunicator());
			if (node_conns)
			{
				unpack_data_vector(buffer, buffer_position, high_conn_size, GetCommunicator());
				unpack_data_vector(buffer, buffer_position, high_conn_nums, GetCommunicator());
			}
			unpack_data_vector(buffer, buffer_position, flags, GetCommunicator());
			REPORT_VAL("buffer position",buffer_position);
			selems[3].reserve(num);
			
			size_t found = 0, marked_for_data = 0, marked_ghost = 0, shift = 0, shift_high = 0;
			for(size_t i = 0; i < num; i++)
			{
				HandleType new_cell = InvalidHandle();
				if( low_conn_size[i] == ENUMUNDEF )
				{
					if( !HaveGlobalID(CELL) ) std::cout << GetProcessorRank() << " should find cell by global id but it is not defined " << std::endl;
					assert(HaveGlobalID(CELL));
					Storage::integer owner = low_conn_nums[shift++];
					Storage::integer gid = low_conn_nums[shift++];
					
					// I have this element and remote processor also has this element, so it should be either in ghost or shared elements
					new_cell = FindSharedGhost(CELL,gid,source,owner);
					if( new_cell == InvalidHandle() ) std::cout << GetProcessorRank() << " cell not found by global id " << std::endl;
					assert(new_cell != InvalidHandle());
				}
				else
				{
					//REPORT_VAL("Unpack low size",low_conn_size[i]);
					c_faces.resize(low_conn_size[i]);
					for(INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
					{
						assert((size_t)low_conn_nums[shift+j] < selems[2].size());
						//REPORT_VAL("Unpack face position",low_conn_nums[shift+j]);
						c_faces.at(j) = selems[2][low_conn_nums[shift+j]];
					}
					shift += low_conn_size[i];
					if (node_conns)
					{
						c_nodes.resize(high_conn_size[i]);
						for (INMOST_DATA_ENUM_TYPE j = 0; j < high_conn_size[i]; j++)
						{
							assert((size_t)high_conn_nums[shift_high + j] < selems[0].size());
							//REPORT_VAL("Unpack node position",high_conn_nums[shift_high+j]);
							c_nodes.at(j) = selems[0][high_conn_nums[shift_high + j]];
						}
						shift_high += high_conn_size[i];
					}
					if( !c_faces.empty() )
						new_cell = FindSharedAdjacency(c_faces.data(),static_cast<enumerator>(c_faces.size()));
				}
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
					//std::cout << "New ghost cell at " << GetProcessorRank() << " arrived from " << source << std::endl;

				} 
				else 
				{
					if( IntegerDF(new_cell,tag_owner) != GetProcessorRank() ) 
					{
						//TODO 46 old
						//unpack_tags[3].push_back(new_cell);
						SetMarker(new_cell,unpack_tags_mrk);
						++marked_for_data;
					}
					++found;
				}

				if ((flags[i] & 1 ? true : false) != GetMarker(new_cell, unpack_tags_mrk))
				{
					std::cout << __FILE__ << ":" << __LINE__ << " on " << GetProcessorRank() << " remote has " << (flags[i] & 1 ? "data" : "no data")
						<< " for cell " << i << " and local expects " << (GetMarker(new_cell, unpack_tags_mrk) ? "data" : "no data") << std::endl;
					REPORT_STR(__FILE__ << ":" << __LINE__ << " on " << GetProcessorRank() << " remote has " << (flags[i] & 1 ? "data" : "no data")
						<< " for cell " << i << " and local expects " << (GetMarker(new_cell, unpack_tags_mrk) ? "data" : "no data"));
				}


				selems[3].push_back(new_cell);
			}
			REPORT_VAL("total found", found << " / " << selems[3].size());
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[3].size());
			REPORT_VAL("total marked ghost", marked_ghost << " / " << selems[3].size());
			REPORT_VAL("buffer position",buffer_position);
			if( marked_for_data != marked_remote )
			{
				std::cout << __FILE__ << ":" << __LINE__ << " marked for data " << marked_for_data << " remote " << marked_remote;
				std::cout << " on " << GetProcessorRank() << " from " << source << std::endl;
			}
			assert(marked_for_data == marked_remote);
		}
		time = Timer() - time;
		REPORT_STR("unpack cells");
		REPORT_VAL("time", time);
		EXIT_BLOCK();
        /////////////////////////////////////////////////////////////
        //unpack esets
		MarkerType mrk_chld = CreateMarker();
		ENTER_BLOCK();
		{
			std::vector<char> names_buff;
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_size;
			std::vector<INMOST_DATA_ENUM_TYPE> high_conn_size;
			std::vector<INMOST_DATA_ENUM_TYPE> low_conn_nums;
			std::vector<Storage::integer> high_conn_nums;
			std::vector<char> flags;
			size_t num = 0, marked_remote = 0;
			unpack_data(buffer,buffer_position,num,GetCommunicator());
			unpack_data(buffer,buffer_position,marked_remote,GetCommunicator());
			unpack_data_vector(buffer,buffer_position,names_buff,GetCommunicator());
			REPORT_VAL("number of sets ",num);
			REPORT_VAL("length of buffer ",names_buff.size());
			unpack_data_vector(buffer,buffer_position,low_conn_size,GetCommunicator());
			unpack_data_vector(buffer,buffer_position,low_conn_nums,GetCommunicator());
			unpack_data_vector(buffer,buffer_position,high_conn_size,GetCommunicator());
			unpack_data_vector(buffer,buffer_position,high_conn_nums,GetCommunicator());
			unpack_data_vector(buffer, buffer_position, flags, GetCommunicator());
			REPORT_VAL("buffer position",buffer_position);
			
			std::vector<std::string> names;
			size_t pos = 0;
			while (pos < names_buff.size())
			{
				names.push_back(std::string(&names_buff[pos]));
				pos += names.back().length() + 1;
			}
			assert(names.size() == num);
			
			//std::cout << "number of sets: " << num << std::endl;
			
			// Looking to all names and create the sets if it's needed
			size_t found = 0, marked_for_data = 0, marked_ghost = 0;
			for (size_t i = 0; i < names.size(); i++)
			{
				//REPORT_VAL("set name", names[i]);
				//REPORT_STR("set name \"" << names[i] << "\"");
				ElementSet set = GetSet(names[i]);
				if (set == InvalidElementSet())
				{
					// std::cout << "Create set " << names[i] << std::endl;
					//REPORT_VAL("create new set", names[i]);
					set = CreateSetUnique(names[i]).first;
					//REPORT_VAL("create new ",set->GetHandle());
					
					assert(set.GetStatus() == Element::Any);
					set.SetStatus(Element::Ghost);
					set.IntegerDF(tag_owner) = -1;
					//TODO 46 old
					//unpack_tags[3].push_back(new_cell);
					set.SetMarker(unpack_tags_mrk);
					++marked_for_data;
					++marked_ghost;
					

					//std::cout << "New ghost eset at " << GetProcessorRank() << " arrived from " << source << std::endl;
				}
				else
				{
					// std::cout << "Found set " << names[i] << std::endl;
					//REPORT_VAL("found old ",set->GetHandle());
					//REPORT_VAL("owner ",IntegerDF(set.GetHandle(),tag_owner));
					if( set.IntegerDF(tag_owner) != GetProcessorRank() )
					{
						//TODO 46 old
						//unpack_tags[3].push_back(new_cell);
						set.SetMarker(unpack_tags_mrk);
						++marked_for_data;
					}
					++found;
				}
				
				//else
				//{
				//	REPORT_VAL("found  old ",set->GetHandle());
				//}

				if ((flags[i] & 1 ? true : false) != set.GetMarker(unpack_tags_mrk))
				{
					std::cout << __FILE__ << ":" << __LINE__ << " on " << GetProcessorRank() << " remote has " << (flags[i] & 1 ? "data" : "no data")
						<< " for set " << i << " and local expects " << (set.GetMarker(unpack_tags_mrk) ? "data" : "no data") << std::endl;
					REPORT_STR(__FILE__ << ":" << __LINE__ << " on " << GetProcessorRank() << " remote has " << (flags[i] & 1 ? "data" : "no data")
						<< " for set " << i << " and local expects " << (set.GetMarker(unpack_tags_mrk) ? "data" : "no data"));
				}

				
				selems[4].push_back(set.GetHandle());
			}
			if( num )
			{
				// Add low conns for sets
				size_t ind = 0, indh = 0;
				for (size_t i = 0; i < num; i++) 
				{
					ElementSet set(this, selems[4][i]);
					// std::cout << "Set " << set.GetName() << " receives " << low_conn_size[i] << " elements ";
					// std::cout << "from " << source;
					//unpack elements
					for (INMOST_DATA_ENUM_TYPE j = 0; j < low_conn_size[i]; j++)
					{
						Storage::integer etype = GetHandleElementNum(low_conn_nums[ind+j]);
						Storage::integer array_pos = GetHandleID(low_conn_nums[ind+j]);
						assert(etype == GetHandleElementNum(selems[etype][array_pos]));
						low_conn_nums[ind+j] = selems[etype][array_pos];
						// std::cout << " " << ElementTypeName(GetHandleElementType(low_conn_nums[ind+j]));
						// std::cout << ":" << GetHandleID(low_conn_nums[ind+j]);
					}
					// std::cout << std::endl;
					if( low_conn_size[i] ) set.AddElements(&low_conn_nums[ind],low_conn_size[i]);
					ind += low_conn_size[i];
					//unpack parent
					if( high_conn_nums[indh] != -1 )
					{
						assert(indh < high_conn_nums.size());
						assert((size_t)high_conn_nums[indh] < selems[4].size());
						ElementSet own_set = ElementSet(this,selems[4][i]);
						ElementSet par_set = ElementSet(this,selems[4][high_conn_nums[indh]]);
						if( own_set->HaveParent() && own_set->GetParent() != par_set  )
						{
							std::cout << __FILE__ << ":" << __LINE__ << "PROBLEM!!!!";
							REPORT_VAL("own set name", own_set->GetName());
							REPORT_VAL("own set handle", own_set->GetHandle());
							REPORT_VAL("par set name", par_set->GetName());
							REPORT_VAL("par set handle", par_set->GetHandle());
							REPORT_VAL("own parent set", own_set->GetParent()->GetName());
							REPORT_VAL("own parent set", own_set->GetParent()->GetHandle());
							exit(-1);
						}
						assert( !own_set->HaveParent() || own_set->GetParent() == par_set );
						if( !own_set->HaveParent() ) par_set->AddChild(own_set);
					}
					//unpack children
					for(ElementSet jt = set.GetChild(); jt.isValid(); jt = jt.GetSibling()) jt.SetMarker(mrk_chld);
					for(INMOST_DATA_ENUM_TYPE j = 1; j < high_conn_size[i]; j++)
					{
						assert(indh+j < high_conn_nums.size());
						assert((size_t)high_conn_nums[indh+j] < selems[4].size());
						ElementSet chld_set = ElementSet(this,selems[4][high_conn_nums[indh+j]]);
						if( !chld_set.GetMarker(mrk_chld) )
						{
							if( chld_set.HaveParent() )
							{
								std::cout << "child " << chld_set.GetName() << " for " << set.GetName() << " has parent " << chld_set.GetParent().GetName() << std::endl;
							}
							set.AddChild(chld_set);
							chld_set.SetMarker(mrk_chld);
						}
					}
					for(ElementSet jt = set.GetChild(); jt.isValid(); jt = jt.GetSibling()) jt.RemMarker(mrk_chld);
					indh += high_conn_size[i];
				}
				// Unpack high conn array
				/*
				for (int i = 0; i < num; i++) if (high_conn_nums[i*3+0] != -1) 
				{
					ElementSet par_set = ElementSet(this,selems[4][i]);
					ElementSet chd_set = ElementSet(this,selems[4][high_conn_nums[i*3+0]]);
					bool chd_exists = false;
					for(ElementSet it = par_set->GetChild(); it->isValid(); it = it->GetSibling())
					{
						if( chd_set == it ) 
						{ 
							chd_exists = true; 
							break; 
						}
					}
					if (!chd_exists) 
						par_set.AddChild(chd_set);
				}
				for (int i = 0; i < num; i++) if (high_conn_nums[i*3+1] != -1) 
				{
					ElementSet own_set = ElementSet(this,selems[4][i]);
					ElementSet par_set = ElementSet(this,selems[4][high_conn_nums[i*3+1]]);
					if( own_set->HaveParent() && own_set->GetParent() != par_set  )
					{
						REPORT_STR("PROBLEM!!!!");
						REPORT_VAL("own set name", own_set->GetName());
						REPORT_VAL("own set handle", own_set->GetHandle());
						REPORT_VAL("par set name", par_set->GetName());
						REPORT_VAL("par set handle", par_set->GetHandle());
						REPORT_VAL("own parent set", own_set->GetParent()->GetName());
						REPORT_VAL("own parent set", own_set->GetParent()->GetHandle());
						exit(-1);
					}
					assert( !own_set->HaveParent() || own_set->GetParent() == par_set );
					if( !own_set->HaveParent() ) par_set->AddChild(own_set);
					assert( own_set->GetParent().GetElementType() == ESET );
				}
				for (int i = 0; i < num; i++) if (high_conn_nums[i*3+2] != -1) 
				{
					ElementSet own_set = ElementSet(this,selems[4][i]);
					if( own_set->HaveParent() ) //if no parent then no sibling!!!
					{
						ElementSet sib_set = ElementSet(this,selems[4][high_conn_nums[i*3+2]]);
						bool sib_exists = false;
						//sibling may be already connected to other set???
						for(ElementSet it = own_set->GetParent()->GetChild(); it->isValid(); it = it->GetSibling())
						{
							if( sib_set == it )
							{
								sib_exists = true;
								break;
							}
						}
						if (!sib_exists) 
							own_set.AddSibling(sib_set);
					}
                }
				*/
            }
			REPORT_VAL("total found", found << " / " << selems[4].size());
			REPORT_VAL("total marked for data", marked_for_data << " / " << selems[4].size());
			REPORT_VAL("total marked ghost", marked_ghost << " / " << selems[4].size());
			REPORT_VAL("buffer position",buffer_position);
			if( marked_for_data != marked_remote )
			{
				std::cout << __FILE__ << ":" << __LINE__ << " marked for data " << marked_for_data << " remote " << marked_remote;
				std::cout << " on " << GetProcessorRank() << " from " << source << std::endl;
			}
			assert(marked_for_data == marked_remote);
		}
		EXIT_BLOCK();
		ReleaseMarker(mrk_chld);
        /////////////////////////////////////////////////////////////
		SetMarker(GetHandle(), unpack_tags_mrk);
		selems[5].push_back(GetHandle()); //data on mesh
		
		//all.reserve(selems[0].size()+selems[1].size()+selems[2].size()+selems[3].size());
		//for(int i = 3; i >= 0; --i)
		//	all.insert(all.end(),selems[i].begin(),selems[i].end());
		//REPORT_VAL("number of all elements",all.size());
		//Storage::integer_array proc = IntegerArrayDV(GetHandle(),tag_processors);
		//Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),source);
		//if( ip == proc.end() || (*ip) != source ) proc.insert(ip,source);
		
		ENTER_BLOCK();
		time = Timer();
		//unpack tags
		{	
			size_t num = 0, namesize;
			unpack_data(buffer,buffer_position,num,GetCommunicator());
			REPORT_VAL("number of tags",num);
			REPORT_VAL("buffer position",buffer_position);
			for(size_t i = 0; i < num; i++)
			{
				INMOST_DATA_ENUM_TYPE defined, sparse, datatype, datasize;
				std::string tag_name;
				Tag tag;
				ENTER_BLOCK();
				unpack_data(buffer,buffer_position,datatype,GetCommunicator());
				unpack_data(buffer,buffer_position,defined,GetCommunicator());
				unpack_data(buffer,buffer_position,sparse,GetCommunicator());
				unpack_data(buffer,buffer_position,datasize,GetCommunicator());
				unpack_data(buffer,buffer_position,namesize,GetCommunicator());
				REPORT_VAL("characters in name",namesize);
				tag_name.resize(namesize);
				unpack_data_array(buffer,buffer_position,&tag_name[0],namesize,GetCommunicator());
				REPORT_VAL("name of tag",tag_name);
				REPORT_VAL("data type",DataTypeName(static_cast<DataType>(datatype)));
				REPORT_VAL("data size",datasize);
				for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype) )
				{
					if( defined & etype ) REPORT_VAL("defined on",ElementTypeName(etype));
					if( sparse & etype ) REPORT_VAL("sparse on",ElementTypeName(etype));
				}
				REPORT_VAL("buffer position",buffer_position);
				tag = CreateTag(tag_name,static_cast<enum DataType>(datatype),static_cast<ElementType>(defined),static_cast<ElementType>(sparse),datasize);
				tag_recv.push_back(tag);
				EXIT_BLOCK();
				//TODO 46 old
				//UnpackTagData(tag,unpack_tags,0,NODE | EDGE | FACE | CELL | ESET, buffer,position,DefaultUnpack);
				UnpackTagData(tag,selems,source,NODE | EDGE | FACE | CELL | ESET | MESH,unpack_tags_mrk, buffer,buffer_position,DefaultUnpack,selems);
				//UnpackTagData(tag,selems,NODE | EDGE | FACE | CELL | ESET,0, buffer,position,DefaultUnpack);
				REPORT_VAL("buffer position",buffer_position);
			}
			
		}
		if (orient)
		{
			for (tag_set::iterator it = tag_recv.begin(); it != tag_recv.end(); ++it)
				for (tag_set::iterator jt = orient_tags.begin(); jt != orient_tags.end(); ++jt)
					if (*it == *jt)
					{
						for (enumerator q = 0; q < selems[2].size(); ++q)
							if (GetMarker(selems[2][q], orient))
							{
								assert(GetMarker(selems[2][q],unpack_tags_mrk));//it should also have data marker
								OrientTag(Face(this, selems[2][q]), *it);
							}
					}
			RemMarkerArray(&selems[2][0], static_cast<enumerator>(selems[2].size()), orient);
			ReleaseMarker(orient);
		}
		for(integer k = ElementNum(NODE); k <= ElementNum(MESH); ++k)
			if( !selems[k].empty() ) RemMarkerArray(&selems[k][0],static_cast<enumerator>(selems[k].size()),unpack_tags_mrk);
		ReleaseMarker(unpack_tags_mrk);
		time = Timer() - time;
		REPORT_STR("unpack tag data");
		REPORT_VAL("time", time);
		EXIT_BLOCK();
#else
		(void) selems;
		(void) buffer;
		(void) source;
		(void) tag_recv;
#endif
		REPORT_VAL("source",source);
		EXIT_FUNC();
	}
	
	
	void Mesh::ExchangeBuffersInner(exch_buffer_type & send_bufs, exch_buffer_type & recv_bufs,
									exch_reqs_type & send_reqs, exch_recv_reqs_type & recv_reqs)
	{
		ENTER_FUNC();
		REPORT_VAL("exchange number", ++num_exchanges);
#if defined(USE_MPI)
		INMOST_DATA_ENUM_TYPE i;
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
		recv_reqs.clear();//resize(recv_bufs.size());
		send_reqs.clear();//resize(send_bufs.size());
		{
			INMOST_DATA_BULK_TYPE stub;
			REPORT_VAL("recv bufs size",recv_bufs.size());
			for(i = 0; i < recv_bufs.size(); i++)// if( !recv_bufs[i].second.empty() )
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + ((size_t)mpirank+mpisize+rand_num))%max_tag;
				//mpi_tag = parallel_mesh_unique_id*mpisize*mpisize+recv_bufs[i].first*mpisize+mpirank;
				INMOST_DATA_BIG_ENUM_TYPE shift = 0, chunk, datasize = recv_bufs[i].second.size();
				int it = 0, mpi_tag_it; // for mpi tag
				REPORT_VAL("mpi_tag",mpi_tag);
				REPORT_VAL("total size",datasize);
				recv_reqs.cnt.push_back(0);
				recv_reqs.buf.push_back(recv_reqs.buf.back());
				do
				{
					MPI_Request req;
					chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),datasize - shift);
					mpi_tag_it = (mpi_tag*1000 + it)%max_tag;
					REPORT_VAL("it",it);
					REPORT_VAL("using mpi_tag",mpi_tag_it);
					REPORT_VAL("size",chunk);
					REPORT_VAL("proc",recv_bufs[i].first);
					REPORT_VAL("empty",recv_bufs[i].second.empty());
					REPORT_MPI(MPI_Irecv(recv_bufs[i].second.empty()?&stub:&recv_bufs[i].second[shift],static_cast<INMOST_MPI_SIZE>(chunk),SEND_AS,recv_bufs[i].first,mpi_tag_it,comm,&req));
					recv_reqs.requests.push_back(req);
					recv_reqs.buf.back()++;
					recv_reqs.cnt.back()++;
					shift += chunk;
					it++;
					if( it >= 1000 )
						std::cout << __FILE__ << ":" << __LINE__ << " too many iterations!!! " << it << " datasize " << datasize << std::endl;
				} while( shift != datasize );
			}
			REPORT_VAL("send bufs size",send_bufs.size());
			for(i = 0; i < send_bufs.size(); i++) //if( !send_bufs[i].second.empty() )
			{
				mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + ((size_t)send_bufs[i].first+mpisize+rand_num))%max_tag;
				//mpi_tag = parallel_mesh_unique_id*mpisize*mpisize+mpirank*mpisize+send_bufs[i].first;
				INMOST_DATA_BIG_ENUM_TYPE shift = 0, chunk, datasize = send_bufs[i].second.size();
				int it = 0, mpi_tag_it; // for mpi tag
				REPORT_VAL("mpi_tag",mpi_tag);
				REPORT_VAL("total size",datasize);
				do
				{
					MPI_Request req;
					chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),datasize - shift);
					mpi_tag_it = (mpi_tag*1000 + it)%max_tag;
					REPORT_VAL("it",it);
					REPORT_VAL("using mpi_tag",mpi_tag_it);
					REPORT_VAL("size",chunk);
					REPORT_VAL("proc",send_bufs[i].first);
					REPORT_VAL("empty",send_bufs[i].second.empty());
					REPORT_MPI(MPI_Isend(send_bufs[i].second.empty()?&stub:&send_bufs[i].second[shift],static_cast<INMOST_MPI_SIZE>(chunk),SEND_AS,send_bufs[i].first,mpi_tag_it,comm,&req));	
					send_reqs.push_back(req);
					shift += chunk;
					it++;
					if( it >= 1000 )
						std::cout << __FILE__ << ":" << __LINE__ << " too many iterations!!! " << it << " datasize " << datasize << std::endl;
				} while( shift != datasize );
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
		ENTER_FUNC();
#if defined(USE_MPI)
        int mpirank = GetProcessorRank(), mpisize = GetProcessorsNumber();
#if defined(USE_MPI_P2P) && defined(PREFFER_MPI_P2P)
		INMOST_DATA_ENUM_TYPE i, end = send_bufs.size();
        REPORT_MPI(MPI_Win_fence(MPI_MODE_NOPRECEDE,window)); //start exchange session
		memset(shared_space,0,sizeof(INMOST_DATA_BIG_ENUM_TYPE)*mpisize); //zero bits where we receive data
		REPORT_MPI(MPI_Win_fence(0,window)); //wait memset finish
		for(i = 0; i < end; i++) shared_space[mpisize+i] = static_cast<INMOST_DATA_BIG_ENUM_TYPE>(send_bufs[i].second.size()+1); //put data to special part of the memory
		for(i = 0; i < end; i++) 
		{
			REPORT_VAL("put value", shared_space[mpisize+i]);
			REPORT_VAL("destination", send_bufs[i].first);
			REPORT_VAL("displacement", mpirank);
			assert( send_bufs[i].first >= 0 && send_bufs[i].first < GetProcessorsNumber() );
			REPORT_MPI(MPI_Put(&shared_space[mpisize+i],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,send_bufs[i].first,mpirank,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,window)); //request rdma
		}
		REPORT_MPI(MPI_Win_fence(MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED,window)); //end exchange session
		if( todo == UnknownSize )
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
			{
				REPORT_VAL("exchange number", ++num_exchanges);
				INMOST_DATA_ENUM_TYPE i;
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
				std::vector<INMOST_DATA_BIG_ENUM_TYPE> send_recv_size(send_bufs.size()+recv_bufs.size());
				std::vector<INMOST_MPI_Request> reqs(send_bufs.size()+recv_bufs.size());
				for(i = 0; i < send_bufs.size(); i++)
					send_recv_size[i+recv_bufs.size()] = static_cast<INMOST_DATA_BIG_ENUM_TYPE>(send_bufs[i].second.size());
				REPORT_VAL("recv buffers size",recv_bufs.size());
				for(i = 0; i < recv_bufs.size(); i++)
				{
					mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + ((size_t)mpirank+mpisize+rand_num))%max_tag;
					REPORT_VAL("origin",recv_bufs[i].first);
					REPORT_VAL("mpi_tag",mpi_tag);
					//mpi_tag = parallel_mesh_unique_id*mpisize*mpisize+recv_bufs[i].first*mpisize+mpirank;
					REPORT_MPI(MPI_Irecv(&send_recv_size[i],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,recv_bufs[i].first,mpi_tag,comm,&reqs[i]));
				}
				REPORT_VAL("send buffers size",send_bufs.size());
				for(i = 0; i < send_bufs.size(); i++)
				{
					mpi_tag = ((parallel_mesh_unique_id+1)*mpisize*mpisize + ((size_t)send_bufs[i].first+mpisize+rand_num))%max_tag;
					REPORT_VAL("destination",send_bufs[i].first);
					REPORT_VAL("mpi_tag",mpi_tag);
					REPORT_VAL("size",send_recv_size[i+recv_bufs.size()]);
					//mpi_tag = parallel_mesh_unique_id*mpisize*mpisize+mpirank*mpisize+send_bufs[i].first;
					REPORT_MPI(MPI_Isend(&send_recv_size[i+recv_bufs.size()],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,send_bufs[i].first,mpi_tag,comm,&reqs[i+recv_bufs.size()]));	
				}
				if( !recv_bufs.empty() )
				{
					REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(recv_bufs.size()),&reqs[0],MPI_STATUSES_IGNORE));
				}
				REPORT_VAL("received buffers size",recv_bufs.size());
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
			INMOST_DATA_ENUM_TYPE i, end = (INMOST_DATA_ENUM_TYPE)send_bufs.size();
            REPORT_MPI(MPI_Win_fence(MPI_MODE_NOPRECEDE,window)); //start exchange session
			memset(shared_space,0,sizeof(INMOST_DATA_BIG_ENUM_TYPE)*mpisize); //zero bits where we receive data
			REPORT_MPI(MPI_Win_fence( 0,window)); //wait memset finish
			for(i = 0; i < end; i++) shared_space[mpisize+i] = static_cast<INMOST_DATA_BIG_ENUM_TYPE>(send_bufs[i].second.size()+1); //put data to special part of the memory
			for(i = 0; i < end; i++) 
			{
                REPORT_VAL("put value", shared_space[mpisize+i]);
                REPORT_VAL("destination", send_bufs[i].first);
                REPORT_VAL("displacement", mpirank);
                assert( send_bufs[i].first >= 0 && send_bufs[i].first < GetProcessorsNumber() );
				REPORT_MPI(MPI_Put(&shared_space[mpisize+i],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,send_bufs[i].first,mpirank,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,window)); //request rdma to target processors for each value
			}
			REPORT_MPI(MPI_Win_fence(MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED,window)); //end exchange session
			recv_bufs.clear();
			for(int ii = 0; ii < mpisize; ii++)
				if( shared_space[ii] > 0 )
				{
					REPORT_VAL("position", ii);
					REPORT_VAL("value", shared_space[ii]);
					recv_bufs.push_back(proc_buffer_type(ii,std::vector<INMOST_DATA_BULK_TYPE>(shared_space[ii]-1))); // this call would be optimized by compiler
				}
			REPORT_VAL("recvs",recv_bufs.size());
#else //USE_MPI_P2P
			int mpisize = GetProcessorsNumber(),mpirank = GetProcessorRank();
			{
				std::vector< INMOST_DATA_BIG_ENUM_TYPE > sends_dest_and_size(send_bufs.size()*2+1);
				for(int i = 0; i < static_cast<int>(send_bufs.size()); i++)
				{
					sends_dest_and_size[i*2+0] = static_cast<INMOST_DATA_BIG_ENUM_TYPE>(send_bufs[i].first);
					sends_dest_and_size[i*2+1] = static_cast<INMOST_DATA_BIG_ENUM_TYPE>(send_bufs[i].second.size());
				}
				INMOST_DATA_ENUM_TYPE recvsize = 0;
				int k,j;
				std::vector<int> allsize(mpisize);
				int size = static_cast<int>(send_bufs.size()*2);
				REPORT_MPI(MPI_Allgather(&size,1,MPI_INT,&allsize[0],1,MPI_INT,comm));
				std::vector<int> displs(mpisize+1,0);
				for(k = 0; k < mpisize; k++)
					recvsize += allsize[k];
				for(k = 1; k < mpisize+1; k++)
					displs[k] = displs[k-1]+allsize[k-1];
				std::vector<INMOST_DATA_BIG_ENUM_TYPE> recvs_dest_and_size(recvsize+1);
				REPORT_MPI(MPI_Allgatherv(&sends_dest_and_size[0],static_cast<INMOST_MPI_SIZE>(send_bufs.size()*2),INMOST_MPI_DATA_BIG_ENUM_TYPE,&recvs_dest_and_size[0],&allsize[0],&displs[0],INMOST_MPI_DATA_BIG_ENUM_TYPE,comm));
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
	
	std::vector<int> Mesh::FinishRequests(exch_recv_reqs_type & recv_reqs)
	{
        std::vector<int> ret;
		ENTER_FUNC();
#if defined(USE_MPI)
		int outcount = 0, ierr;
		REPORT_VAL("requests",recv_reqs.requests.size());
		if( recv_reqs.requests.empty() ) 
			ret.clear();
		else
		{
			std::vector<int> retreq;
			while( ret.empty() && recv_reqs.count() )
			{
				ENTER_BLOCK();
				retreq.resize(recv_reqs.requests.size(),-1);
				EXIT_BLOCK();
				ENTER_BLOCK();
#if defined(USE_PARALLEL_WRITE_TIME)
				std::stringstream waitfrom;
				for (size_t k = 0; k < recv_reqs.requests.size(); ++k)
					if (recv_reqs.cnt[k])
					{
						MPI_Status stat;
						int flag = 0;
						MPI_Request_get_status(recv_reqs.requests[k], &flag, &stat);
						waitfrom << " src " << stat.MPI_SOURCE << " err " << stat.MPI_ERROR << " tag " << stat.MPI_TAG << " " << (flag ? "completed" : "pending") << std::endl;
					}
				REPORT_STR(waitfrom.str());
#endif
				REPORT_MPI(ierr = MPI_Waitsome(static_cast<INMOST_MPI_SIZE>(recv_reqs.requests.size()),&recv_reqs.requests[0],&outcount,&retreq[0],MPI_STATUSES_IGNORE));
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				if( outcount == MPI_UNDEFINED ) 
				{
					std::cout << __FILE__ << ":" << __LINE__ << " no more requests, error? " << std::endl;
					retreq.clear();
					break;
				}
				else retreq.resize(outcount);
				EXIT_BLOCK();
				ENTER_BLOCK();
				for(size_t i = 0; i < retreq.size(); ++i)
				{
					ENTER_BLOCK();
					size_t p = retreq[i]; //corresponds to closed request
					//determine number of buffer corresponding to the request
					size_t k = std::lower_bound(recv_reqs.buf.begin(),recv_reqs.buf.end(),p) - recv_reqs.buf.begin();
					//decrase amount of messages we wait for
					recv_reqs.cnt[k]--;
					if( recv_reqs.cnt[k] == 0 ) // all messages received
						ret.push_back((int)k);
					EXIT_BLOCK();
				}
				EXIT_BLOCK();
			}
		}
#else
		(void) recv_reqs;
#endif
		EXIT_FUNC();
        return ret;
	}
	

	void Mesh::Barrier()
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		MPI_Barrier(GetCommunicator());
#endif
		EXIT_FUNC();
	}

	void Mesh::ListTags(tag_set & list)
	{
		for(iteratorTag t = BeginTag(); t != EndTag(); ++t) list.push_back(*t);
	}
	
	void Mesh::ExchangeMarked(enum Action action)
	{
		ENTER_FUNC();
		bool dely_delete = false;
		if( m_state == Serial ) return;
#if defined(USE_MPI)
        int mpirank = GetProcessorRank();
		tag_set tag_list, tag_list_recv;
		tag_set tag_list_empty;
		//std::vector<MPI_Request> send_reqs, recv_reqs;
		//exch_buffer_type send_bufs;
		//exch_buffer_type recv_bufs;
		exchange_data storage;
		proc_elements_by_type send_elements;
		std::vector<int> done;

		if( action == AGhost ) REPORT_STR("Ghosting algorithm")
		else if( action == AMigrate ) REPORT_STR("Migration algorithm")
		
		ListTags(tag_list);
		{
			tag_set::iterator it = tag_list.begin();
			while( it != tag_list.end() )
			{
				if( it->GetTagName().substr(0,9) == "PROTECTED" ) 
					it = tag_list.erase(it);
            	//else if(it->GetDataType() == DATA_REFERENCE)
				//	it = tag_list.erase(it);
				else if(it->GetDataType() == DATA_REMOTE_REFERENCE)
					it = tag_list.erase(it);
				else it++;
			}
		}
		{
			REPORT_STR("Gathering elements to send");
			for (ElementType etype = NODE; etype <= ESET; etype = NextElementType(etype))
			{
				int itype = ElementNum(etype);
				//for (iteratorElement it = BeginElement(etype); it != EndElement(); it++)
				for(integer i = 0; i < LastLocalID(etype); ++i) if( isValidElement(etype,i))
				{
					Element it = ElementByLocalID(etype, i);
					Storage::integer_array mark = it->IntegerArray(tag_sendto);
					for (Storage::integer_array::iterator kt = mark.begin(); kt != mark.end(); kt++)
					{
						assert(*kt >= 0 && *kt < GetProcessorsNumber());
						if (*kt != mpirank) send_elements[*kt][itype].push_back(it->GetHandle());
					}
					it->DelData(tag_sendto);
				}
			}
		}
		storage.send_buffers.resize(send_elements.size());
        //std::cout << mpirank << ": Send size: " << send_elements.size() << std::endl;
		ENTER_BLOCK();
		REPORT_STR("Packing elements to send");
		
		ENTER_BLOCK();
		REPORT_STR("pack elements");
#if defined(USE_OMP)
//#pragma omp parallel
#endif
		{
			TagInteger pack_position;
			std::stringstream tag_name;
			tag_name << "PROTECTED_TEMPORARY_PACK_POSITION_" << GetLocalProcessorRank();
			pack_position = CreateTag(tag_name.str(), DATA_INTEGER, ESET | CELL | FACE | EDGE | NODE, ESET | CELL | FACE | EDGE | NODE, 1);
#if defined(USE_OMP)
//#pragma omp for
#endif
			for (int i = 0; i < (int)send_elements.size(); ++i) //for(proc_elements_by_type::iterator it = send_elements.begin(); it != send_elements.end(); it++)
			{
				proc_elements_by_type::iterator it = send_elements.begin();
				{ int k = i; while (k) ++it, --k; } //advance
				if (!it->second.empty())
				{
					PackElementsGather(it->second, it->second, it->first, ESET | CELL | FACE | EDGE | NODE, 0, tag_list, false, action == AGhost || dely_delete, action == AGhost || dely_delete);
					PackElementsEnumerate(it->second, pack_position);
					PackElementsData(it->second, storage.send_buffers[i].second, it->first, tag_list, pack_position, action == AGhost || dely_delete);
					PackElementsUnenumerate(it->second, pack_position);
					storage.send_buffers[i].first = it->first;
				}
				else storage.send_buffers[i].first = -1;
			}
			DeleteTag(pack_position);
		}
		EXIT_BLOCK();
		ENTER_BLOCK();
		{
			exch_buffer_type::iterator it = storage.send_buffers.begin();
			while (it != storage.send_buffers.end())
			{
				if (it->first == -1)
					it = storage.send_buffers.erase(it);
				else ++it;
			}
		}
		EXIT_BLOCK();
		EXIT_BLOCK();

		int num_wait = (int)storage.send_buffers.size();

		PrepareReceiveInner(UnknownSource, storage.send_buffers,storage.recv_buffers);
		ExchangeBuffersInner(storage.send_buffers,storage.recv_buffers,storage.send_reqs,storage.recv_reqs);
			
		if( !dely_delete && action == AMigrate ) // delete packed elements
		{
			Tag tag_new_processors = GetTag("TEMPORARY_NEW_PROCESSORS");
			/*
			ENTER_BLOCK();
			REPORT_STR("Gathering elements to delete");
			{
				elements_by_type delete_elements;
				MarkerType delete_marker = CreateMarker();
				for(proc_elements_by_type::iterator it = send_elements.begin(); it != send_elements.end(); it++)
				{
					REPORT_VAL("for processor",it->first);
					REPORT_VAL("number of elements",it->second.size());
					if( !it->second.empty() )
					{
						//Have packed all entities, now delete those entities that i should not own	
						for(ElementType etype = NODE; etype <= ESET; etype = NextElementType(etype)) if( tag_new_processors.isDefined(etype) )
						{
							int k = ElementNum(etype);
							for(element_set::iterator kt = it->second[k].begin(); kt != it->second[k].end(); kt++) if( !GetMarker(*kt,delete_marker) )
							{
								Storage::integer_array procs = IntegerArray(*kt,tag_new_processors);
								if( procs.empty() ) continue; //don't delete elements without processors
								if( !std::binary_search(procs.begin(),procs.end(),mpirank) )
								{
									delete_elements[GetHandleElementNum(*kt)].push_back(*kt);
									SetMarker(*kt,delete_marker);
								}
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
				REPORT_VAL("number of esets to delete",delete_elements[4].size());
				for(int j = ElementNum(ESET); j >= ElementNum(NODE); j--)
				{
					//this is not really needed because we destroy those elements
					for(element_set::iterator kt = delete_elements[j].begin(); kt != delete_elements[j].end(); ++kt) RemMarker(*kt,delete_marker);
					//for(element_set::iterator kt = delete_elements[j].begin(); kt != delete_elements[j].end(); ++kt) Destroy(*kt);
					for(element_set::iterator kt = delete_elements[j].begin(); kt != delete_elements[j].end(); ++kt) Delete(*kt);
				}
				ReleaseMarker(delete_marker);
			}
			EXIT_BLOCK();
			*/
			ENTER_BLOCK();
			/*
			//This is not needed, we delete links in sets in RemoveLinksToDeletedElements later
			REPORT_STR("Deleting links to deleted elements from sets");
			for(iteratorSet it = BeginSet(); it != EndSet(); ++it)
			{
				ElementSet::iterator jt = it->Begin();
				while(jt != it->End())
				{
					Storage::integer_array procs = jt->IntegerArrayDV(tag_new_processors);
					if( !std::binary_search(procs.begin(),procs.end(),mpirank) )
						jt = it->Erase(jt);
					else jt++;
				}
			}
			*/
			REPORT_STR("Deleting not owned elements");
			//now delete elements that we have not yet deleted
			MarkerType del = CreateMarker();
			int count = 0;
			for(ElementType etype = ESET; etype >= NODE; etype = PrevElementType(etype)) if( tag_new_processors.isDefined(etype) )
			{
				int ecount = 0;
				//for(iteratorElement it = BeginElement(etype); it != EndElement(); it++)
				for(integer i = 0; i < LastLocalID(etype); ++i) if(isValidElement(etype,i))
				{
					Element it = ElementByLocalID(etype, i);
					Storage::integer_array procs = it->IntegerArray(tag_new_processors);
					if( procs.empty() ) continue;//don't delete elements without processors
					if( !std::binary_search(procs.begin(),procs.end(),mpirank) ) 
					{
						it->SetMarker(del);
						++ecount;
					}
				}
				REPORT_STR("number of " << ElementTypeName(etype) << " to delete " << ecount);
				count += ecount;
			}

			RemoveLinksToDeletedElements(del);

			for (ElementType etype = ESET; etype >= NODE; etype = PrevElementType(etype))
			{
				//for (iteratorElement it = BeginElement(etype); it != EndElement(); it++)
				for (integer i = 0; i < LastLocalID(etype); ++i) if (isValidElement(etype, i))
				{
					Element it = ElementByLocalID(etype, i);
					if (it->GetMarker(del))
					{
						if (etype == ESET)
							it->getAsSet().DeleteSet();
						else it->Delete();
					}
				}
			}

			ReleaseMarker(del);
			REPORT_VAL("number of some other",count);
			REPORT_STR("Done deleting");
			EXIT_BLOCK();
		}
		send_elements.clear();
		
		ENTER_BLOCK();
		REPORT_STR("Unpacking received data");
		while( !(done = FinishRequests(storage.recv_reqs)).empty() )
		{
			for(std::vector<int>::iterator qt = done.begin(); qt != done.end(); qt++)
			{
				elements_by_type recv_elements;
				REPORT_STR("call unpack");
				size_t position = 0;
				UnpackElementsData(recv_elements,storage.recv_buffers[*qt].second,storage.recv_buffers[*qt].first,position,tag_list_recv);
				if( action == AGhost )
				{
					for(int k = 0; k < 5; ++k)
						for(element_set::iterator it = recv_elements[k].begin(); it != recv_elements[k].end(); it++)
						{
							Storage::integer owner = IntegerDF(*it,tag_owner);
							if( owner == mpirank ) continue;
							Storage::integer_array proc = IntegerArrayDV(*it,tag_processors);	
							Storage::integer_array::iterator ip = std::lower_bound(proc.begin(),proc.end(),mpirank);
							if( ip == proc.end() || (*ip) != mpirank ) 
							{
								proc.insert(ip,mpirank);
								send_elements[owner][GetHandleElementNum(*it)].push_back(*it);
							}
						}
				}
			}
		}
		REPORT_STR("Ended recv");
		EXIT_BLOCK();
				
		if( action == AGhost ) //second round to inform owner about ghost elements
		{
			ENTER_BLOCK();
			REPORT_STR("Second round for ghost exchange");
			if( !storage.send_reqs.empty() )
			{
				REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(storage.send_reqs.size()),&storage.send_reqs[0],MPI_STATUSES_IGNORE));
			}
			//Inform owners about the fact that we have received their elements
			InformElementsOwners(send_elements,storage);
			EXIT_BLOCK();
		}
		
		if( action == AMigrate ) //Compute new values
		{
			Tag tag_new_owner = GetTag("TEMPORARY_NEW_OWNER");
			Tag tag_new_processors = GetTag("TEMPORARY_NEW_PROCESSORS");
			/*
			for(iteratorTag t = BeginTag(); t != EndTag(); t++)
			{
				if( t->GetTagName().substr(0,9) == "PROTECTED" ) continue;
				if( t->GetDataType() == DATA_REFERENCE )
				{
					for(ElementType etype = NODE; etype <= ESET; etype = NextElementType(etype)) if( t->isDefined(etype) )
						for(INMOST_DATA_ENUM_TYPE eit = 0; eit < LastLocalID(etype); ++eit) if( isValidElement(etype,eit) )
						{
							Element it = ElementByLocalID(etype,eit);
							Storage::integer owner = it->Integer(tag_new_owner);
							Storage::integer_array procs = it->IntegerArray(tag_new_processors);
							if( std::binary_search(procs.begin(),procs.end(),owner) )
							{
								Storage::reference_array refs = it->ReferenceArray(*t);
								for(Storage::reference_array::iterator jt = refs.begin(); jt != refs.end(); ++jt) if( jt->isValid() )
								{
									procs = jt->IntegerArray(tag_new_processors);
									Storage::integer_array::iterator find = std::lower_bound(procs.begin(),procs.end(),owner);
									if( find == procs.end() || *find != owner ) 
									{
										std::cout << t->GetTagName() << " reference from " << ElementTypeName(it->GetElementType()) << ":" << it->LocalID();
										std::cout << " to " << ElementTypeName(jt->GetElementType()) << ":" << jt->LocalID() << std::endl;
									}
								}
							}
						}
				}
			}
			 */
			
			CheckSetLinks(__FILE__,__LINE__);

			ENTER_BLOCK();
			REPORT_STR("Second round for elements migration");
			REPORT_STR("Computing new values");
			/*
			//this is not needed (and there is a bug!)
			if( dely_delete )
			{
				for(iteratorSet it = BeginSet(); it != EndSet(); ++it)
				{
					ElementSet::iterator jt = it->Begin();
					while(jt != it->End())
					{
						Storage::integer_array procs = jt->IntegerArrayDV(tag_new_processors);
						if( !std::binary_search(procs.begin(),procs.end(),mpirank) )
							jt = it->Erase(jt);
						else jt++;
					}
				}
			}
			*/
			MarkerType del = CreateMarker();
			for(ElementType etype = ESET; etype >= NODE; etype = PrevElementType(etype))
			{
				if( tag_new_owner.isDefined(etype) && tag_new_processors.isDefined(etype) )
				{
					//for(iteratorElement it = BeginElement(etype); it != EndElement(); it++)
					for (integer i = 0; i < LastLocalID(etype); ++i) if (isValidElement(etype, i))
					{
						Element it = ElementByLocalID(etype, i);
						integer_array new_procs = it->IntegerArrayDV(tag_new_processors);
						integer_array old_procs = it->IntegerArrayDV(tag_processors);
						if( new_procs.empty() ) continue;
						integer new_owner = it->IntegerDF(tag_new_owner);
						it->IntegerDF(tag_owner) = new_owner;
						old_procs.replace(old_procs.begin(),old_procs.end(),new_procs.begin(),new_procs.end());
						if( new_owner == mpirank )
						{
							if( old_procs.size() == 1 ) //there is only one processors and that's me
								it->SetStatus(Element::Owned);
							else
								it->SetStatus(Element::Shared);
						}
						else it->SetStatus(Element::Ghost);
						if( dely_delete && !std::binary_search(old_procs.begin(),old_procs.end(),mpirank) )
							it->SetMarker(del);
					}
				}
			}

			RemoveLinksToDeletedElements(del);

			for (ElementType etype = ESET; etype >= NODE; etype = PrevElementType(etype))
			{
				//for (iteratorElement it = BeginElement(etype); it != EndElement(); it++)
				for (integer i = 0; i < LastLocalID(etype); ++i) if (isValidElement(etype, i))
				{
					Element it = ElementByLocalID(etype, i);
					if (it->GetMarker(del))
					{
						if (etype == ESET)
							it->getAsSet().DeleteSet();
						else it->Delete();
					}
				}
			}

			ReleaseMarker(del);
			EXIT_BLOCK();
			CheckSetLinks(__FILE__,__LINE__);
		}
		
		CheckProcsSorted(__FILE__,__LINE__);
		
		RecomputeParallelStorage(ESET | CELL | FACE | EDGE | NODE);
		
		CheckGhostSharedCount(__FILE__,__LINE__);
		CheckOwners(__FILE__, __LINE__);
		CheckGIDs(__FILE__, __LINE__);
		CheckCentroids(__FILE__,__LINE__);

		if (action == AGhost)
		{
			//Probably now owner should send processors_tag data
			ExchangeData(tag_processors, ESET | CELL | FACE | EDGE | NODE, 0);
		}
		CheckCentroids(__FILE__, __LINE__);
		CheckOwners(__FILE__, __LINE__);
		CheckGIDs(__FILE__, __LINE__);
		CheckProcessors();
		
		//ComputeSharedProcs();
		
		if( action == AMigrate ) 
		{
			if( !storage.send_reqs.empty() )
			{
				REPORT_MPI(MPI_Waitall(static_cast<INMOST_MPI_SIZE>(storage.send_reqs.size()),&storage.send_reqs[0],MPI_STATUSES_IGNORE));
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
			for(int i = 0; i < 5; i++) if( mask & ElementTypeFromDim(i) )
				it->second[i].clear();
		for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)			
			for(int i = 0; i < 5; i++) if( mask & ElementTypeFromDim(i) )
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

	void Mesh::CheckGhostSharedCount(std::string file, int line, ElementType etype)
	{
		(void)file,(void)line,(void)etype; //no warning in release
#if !defined(NDEBUG)
		ENTER_FUNC();
#if defined(USE_MPI)
#if !defined(USE_PARALLEL_STORAGE)
		parallel_storage ghost_elements, shared_elements;
		GatherParallelStorage(ghost_elements,shared_elements,etype);
#endif //USE_PARALLEL_STORAGE
		int size = GetProcessorsNumber();
		int rank = GetProcessorRank();
		std::vector<int> send_counts(size*5*2,0);
		std::vector<int> recv_counts(size*size*5*2);
		for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); ++it)
		{
			for(int k = 0; k < 5; ++k) if( ElementTypeFromDim(k) & etype )
				send_counts[it->first*5+k] = (int)it->second[k].size();
		}
		for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); ++it)
		{
			for(int k = 0; k < 5; ++k) if( ElementTypeFromDim(k) & etype )
				send_counts[it->first*5+k + size*5] = (int)it->second[k].size();
		}
		REPORT_MPI(MPI_Allgather(&send_counts[0],size*5*2,MPI_INT,&recv_counts[0],size*5*2,MPI_INT,GetCommunicator()));
		//check only my against others
		//if I have N shared elements with processor I, then procesoor I should have N ghost elements
		Storage::integer err = 0;
		for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); ++it)
		{
			for(int k = 0; k < 5; ++k) if( ElementTypeFromDim(k) & etype && static_cast<size_t>(recv_counts[it->first*size*5*2 + rank*5+k + size*5]) != it->second[k].size() )
			{
				std::cout << "processor " << GetProcessorRank() << " has " << it->second[k].size() << " of shared " << ElementTypeName(ElementTypeFromDim(k));
				std::cout  << " for " << it->first << " but there are " << recv_counts[it->first*size*5*2 + rank*5+k + size*5] << " ghost elements on " << it->first << std::endl;
				REPORT_STR("processor " << GetProcessorRank() << " has " << it->second[k].size() << " of shared " << ElementTypeName(ElementTypeFromDim(k)) <<
					" for " << it->first << " but there are " << recv_counts[it->first*size*5*2 + rank*5+k + size*5] << " ghost elements on " << it->first);
				err++;
			}
			/*
			else if( ElementTypeFromDim(k)& etype )
			{
				std::cout << "processor " << GetProcessorRank() << " has " << it->second[k].size() << " of shared " << ElementTypeName(ElementTypeFromDim(k));
				std::cout << " for " << it->first << " and there are " << recv_counts[it->first * size * 5 * 2 + rank * 5 + k + size * 5] << " ghost elements on " << it->first;
				std::cout << " from " << file << ":" << line << std::endl;
			}
			*/
		}
		for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); ++it)
		{
			for(int k = 0; k < 5; ++k) if( ElementTypeFromDim(k) & etype && static_cast<size_t>(recv_counts[it->first*size*5*2 + rank*5+k]) != it->second[k].size() )
			{
				std::cout << "processor " << GetProcessorRank() << " has " << it->second[k].size() << " of ghost " << ElementTypeName(ElementTypeFromDim(k));
				std::cout  << " for " << it->first << " but there are " << recv_counts[it->first*size*5*2 + rank*5+k] << " shared elements on " << it->first << std::endl;
				REPORT_STR("processor " << GetProcessorRank() << " has " << it->second[k].size() << " of ghost " << ElementTypeName(ElementTypeFromDim(k)) <<
					" for " << it->first << " but there are " << recv_counts[it->first*size*5*2 + rank*5+k] << " shared elements on " << it->first);
				err++;
			}
			/*
			else if (ElementTypeFromDim(k) & etype)
			{
				std::cout << "processor " << GetProcessorRank() << " has " << it->second[k].size() << " of ghost " << ElementTypeName(ElementTypeFromDim(k));
				std::cout << " for " << it->first << " and there are " << recv_counts[it->first * size * 5 * 2 + rank * 5 + k] << " shared elements on " << it->first;
				std::cout << " from " << file << ":" << line << std::endl;
			}
			*/
		}
		err = Integrate(err);
#endif //USE_MPI
		EXIT_FUNC();
#if defined(USE_MPI)
		if( err ) 
		{
			std::cout << "crash from " << file << ":" << line << std::endl;
			REPORT_STR("crash from " << file << ":" << line);
			std::cout.flush();
			assert(false);
			exit(-1);
		}
#endif //USE_MPI
#endif //NDEBUG
	}

	void Mesh::CheckProcsSorted(std::string file, int line)
	{
		(void)file,(void)line;
#if !defined(NDEBUG)
		ENTER_FUNC();
		int err = 0;
		if( ProcessorsTag().isValid() )
		{
			for(ElementType etype = NODE; etype <= ESET; etype = NextElementType(etype))
			{
				//std::cout << GetProcessorRank() << " start " << ElementTypeName(etype) << std::endl;
				for(int k = 0; k < LastLocalID(etype); ++k) if( isValidElement(etype,k))
				{
					Element e = ElementByLocalID(etype,k);
					Storage::integer_array arr = e.IntegerArray(ProcessorsTag());
					//std::cout << ElementTypeName(etype) << " " << k << " " << arr.size() << std::endl;
					for(INMOST_DATA_ENUM_TYPE q = 1; q < arr.size(); q++)
						if( arr[q-1] >= arr[q] ) err++;
				}
				//std::cout << GetProcessorRank() << " finish " << ElementTypeName(etype) << std::endl;
			}
		}
		if( err ) std::cout << "error on " << GetProcessorRank() << std::endl;
		//err = Integrate(err);
		if( err ) 
		{
			std::cout << "crash from " << file << ":" << line << std::endl;
			REPORT_STR("crash from " << file << ":" << line);
			std::cout.flush();
			exit(-1);
		}
		EXIT_FUNC();
#endif //NDEBUG
	}
	
	void Mesh::CheckSetLinks(std::string file, int line)
	{
		(void)file,(void)line;
		ENTER_FUNC();
#if !defined(NDEBUG)
		Storage::integer err = 0;
		for(iteratorSet it = BeginSet(); it != EndSet(); ++it)
		{
			const Element::adj_type & lc = LowConn(*it);
			for(Element::adj_type::const_iterator jt = lc.begin(); jt != lc.end(); ++jt)
				if( *jt != InvalidHandle() && !isValidElement(*jt) )
				{
					std::cout << "set " << it->GetName() << " has bad link to " << ElementTypeName(GetHandleElementType(*jt)) << ":" << GetHandleID(*jt) << std::endl;
					REPORT_STR("set " << it->GetName() << " has bad link to " << ElementTypeName(GetHandleElementType(*jt)) << ":" << GetHandleID(*jt));
					err++;
				}
		}
		err = Integrate(err);
		if( err )
		{
			std::cout << file << ":" <<line <<  " " << err << " invalid links in sets" << std::endl;
			REPORT_STR(file << ":" <<line <<  " " << err << " invalid links in sets");
			std::cout.flush();
			exit(-1);
		}
#endif //NDEBUG
		EXIT_FUNC();
	}
	
	void Mesh::CheckCentroids(std::string file, int line)
	{
		(void)file,(void)line; //no warning
#if !defined(NDEBUG)
		ENTER_FUNC();
		TagRealArray checker = CreateTag("CHECK_CENTROIDS",DATA_REAL,CELL|FACE|EDGE|NODE,CELL|FACE|EDGE|NODE,3);
		for(Mesh::iteratorElement it = BeginElement(CELL|FACE|EDGE|NODE); it != EndElement(); ++it) if( it->GetStatus() & (Element::Ghost | Element::Shared) )
			it->Centroid(checker[*it].data());
		ExchangeData(checker,CELL|FACE|EDGE|NODE,0);
		int bad[4] = {0,0,0,0}, total[4] = {0,0,0,0};
		bool have_gid[4] = { HaveGlobalID(NODE), HaveGlobalID(EDGE), HaveGlobalID(FACE), HaveGlobalID(CELL) };
		for(Mesh::iteratorElement it = BeginElement(CELL|FACE|EDGE|NODE); it != EndElement(); ++it) if( it->GetStatus() & (Element::Ghost | Element::Shared) )
		{
			INMOST_DATA_REAL_TYPE cnt[3] = { 0,0,0 };
			it->Centroid(cnt);
			bool problem = false;
			double dist = 0;
			for (int k = 0; k < 3; ++k)
			{
				if (fabs(cnt[k] - checker[*it][k]) > 1.0e-9)
					problem = true;
				dist += pow(cnt[k] - checker[*it][k], 2);
			}
			if( problem )
			{
				dist = sqrt(dist);
				bad[it->GetElementNum()]++;
				std::cout << "bad centroid on " << GetProcessorRank() << " " << ElementTypeName(it->GetElementType()) 
					<< ":" << it->LocalID() << " gid " << (have_gid[it->GetElementNum()]?it->GlobalID():-1) << " "
					<< cnt[0] << " " << cnt[1] << " " << cnt[2] << " != " 
					<< checker[*it][0] << " " << checker[*it][1] << " " << checker[*it][2] 
					<< " diff " << cnt[0] - checker[*it][0] << " " << cnt[1] - checker[*it][1] << " " << cnt[2] - checker[*it][2]
					<< " dist " << dist << std::endl;
				REPORT_STR("bad centroid on " << GetProcessorRank() << " " << ElementTypeName(it->GetElementType())
					<< ":" << it->LocalID() << " gid " << (have_gid[it->GetElementNum()] ? it->GlobalID() : -1) << " "
					<< cnt[0] << " " << cnt[1] << " " << cnt[2] << " != "
					<< checker[*it][0] << " " << checker[*it][1] << " " << checker[*it][2]
					<< " diff " << cnt[0] - checker[*it][0] << " " << cnt[1] - checker[*it][1] << " " << cnt[2] - checker[*it][2] 
					<< " dist " << dist);
			}
			total[it->GetElementNum()]++;
		}
		DeleteTag(checker);
		Storage::integer crash = 0;
		for(int k = 0; k < 4; ++k) if( bad[k] )
		{
			std::cout << "processor " << GetProcessorRank() << " on " <<  ElementTypeName(ElementTypeFromDim(k)) << " bad centroids " << bad[k] << "/" << total[k] << " " << (HaveGlobalID(ElementTypeFromDim(k)) ? "have gid" : "no gid") << std::endl;
			REPORT_STR("processor " << GetProcessorRank() << " on " <<  ElementTypeName(ElementTypeFromDim(k)) << " bad centroids " << bad[k] << "/" << total[k] << " " << (HaveGlobalID(ElementTypeFromDim(k)) ? "have gid" : "no gid"));
			crash++;
		}
		crash = Integrate(crash);
		if( crash ) 
		{
			std::cout << "crash from " << file << ":" << line << std::endl;
			std::cout.flush();
			REPORT_STR("crash from " << file << ":" << line);
			std::cout.flush();
			exit(-1);
		}
		EXIT_FUNC();
#endif //NDEBUG
	}


	void UnpackCheckOwner(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( element->Integer(tag) != *(const Storage::integer *)data || size != 1 )
		{
			if( size != 1 ) 
				std::cout << __FILE__ << ":" << __LINE__ <<  " wrong size of array " << std::endl;
			std::cout << __FILE__ << ":" << __LINE__ <<  " wrong owner on " << tag.GetMeshLink()->GetProcessorRank();
			std::cout << " element " << ElementTypeName(element->GetElementType()) << ":" << element->LocalID();
			std::cout << " owner " << element->Integer(tag) << " remote " << *(const Storage::integer*)data;
			std::cout << std::endl;
			std::cout.flush();
			exit(-1);
		}
	}

	
	void Mesh::CheckOwners(std::string file, int line)
	{
#if !defined(NDEBUG)
		ENTER_FUNC();
		int bad = 0;
		for (Mesh::iteratorElement it = BeginElement(ESET | CELL | FACE | EDGE | NODE); it != EndElement(); ++it)
		{
			integer owner = Integer(*it, tag_owner);
			Element::Status s = GetStatus(*it);
			if (owner == GetProcessorRank() && s == Element::Ghost) bad++;
			if (owner != GetProcessorRank() && (s == Element::Owned || s == Element::Shared)) bad++;
		}
		Storage::integer crash = 0;
		if(bad)
		{
			std::cout << "processor " << GetProcessorRank() << " bad owner/status " << bad << std::endl;
			REPORT_STR("processor " << GetProcessorRank() << " bad owner/status " << bad);
			crash++;
		}
		crash = Integrate(crash);
		if (crash)
		{
			std::cout << "crash " << __FUNCTION__ << " " << file << ":" << line << std::endl;
			REPORT_STR("crash from " << __FUNCTION__);
			std::cout.flush();
			exit(-1);
		}
		ReduceData(tag_owner,ESET|CELL|FACE|EDGE|NODE,0,UnpackCheckOwner);
		EXIT_FUNC();
#endif
	}

	void UnpackCheckGID(const Tag& tag, const Element& element, const INMOST_DATA_BULK_TYPE* data, INMOST_DATA_ENUM_TYPE size)
	{
		if (element->Integer(tag) != *(const Storage::integer*)data || size != 1)
		{
			if (size != 1)
				std::cout << __FILE__ << ":" << __LINE__ << " wrong size of array " << std::endl;
			std::cout << __FILE__ << ":" << __LINE__ << " wrong GID on " << tag.GetMeshLink()->GetProcessorRank();
			std::cout << " element " << ElementTypeName(element->GetElementType()) << ":" << element->LocalID();
			std::cout << " GID " << element->Integer(tag) << " remote " << *(const Storage::integer*)data;
			std::cout << std::endl;
			std::cout.flush();
			exit(-1);
		}
	}

	void Mesh::CheckGIDs(std::string file, int line, ElementType mask)
	{
#if !defined(NDEBUG)
		ENTER_FUNC();
		ElementType exch = NONE;
		int bad = 0;
		for (ElementType etype = NODE; etype <= ESET; etype = NextElementType(etype))
		{
			if ((etype & mask) && HaveGlobalID(etype))
			{
				exch |= etype;
				std::vector<integer> gids;
				gids.reserve(NumberOf(etype));
				for (Mesh::iteratorElement it = BeginElement(etype); it != EndElement(); ++it)
					gids.push_back(it->GlobalID());
				size_t size0 = gids.size();
				std::sort(gids.begin(), gids.end());
				gids.resize(std::unique(gids.begin(), gids.end()) - gids.begin());
				if (size0 != gids.size())
				{
					bad++;
					std::cout << ElementTypeName(etype) << " on " << GetProcessorRank() << " has duplicate global ids! " << size0 << " != " << gids.size() << std::endl;
					REPORT_STR("processor " << GetProcessorRank() << " etype " << ElementTypeName(etype) << " on " << GetProcessorRank() << " has duplicate global ids!" << size0 << " != " << gids.size());
				}
			}
		}
		Storage::integer crash = 0;
		if (bad)
		{
			std::cout << "processor " << GetProcessorRank() << " duplicate gids in element types " << bad << std::endl;
			REPORT_STR("processor " << GetProcessorRank() << " duplicate gids in element types  " << bad);
			crash++;
		}
		crash = Integrate(crash);
		if (crash)
		{
			std::cout << "crash " << __FUNCTION__ << " " << file << ":" << line << std::endl;
			REPORT_STR("crash from " << __FUNCTION__);
			std::cout.flush();
			exit(-1);
		}
		if( exch != NONE ) ReduceData(tag_global_id, exch, 0, UnpackCheckGID);
		EXIT_FUNC();
#endif
	}

	bool crash_status = 0;

	void UnpackCheckProcessors(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		Storage::integer_array procs = element->IntegerArray(tag);
		bool compare = true;
		if( size == procs.size() )
		{
			const INMOST_DATA_INTEGER_TYPE * rdata = (const INMOST_DATA_INTEGER_TYPE *)data;
			int skip = 0;
#if defined(USE_PARALLEL_STORAGE)
			skip = 1;
#endif
			for(INMOST_DATA_ENUM_TYPE k = 0; k < size-skip; ++k)
				if( procs[k] != rdata[k] )
					compare = false;
		}
		else compare = false;
		if( !compare )
		{
			if( size != procs.size() ) 
				std::cout << __FILE__ << ":" << __LINE__ <<  " wrong size of array " << std::endl;
			std::cout << __FILE__ << ":" << __LINE__ <<  " on " << tag.GetMeshLink()->GetProcessorRank();
			INMOST_DATA_REAL_TYPE cnt[3] = { 0,0,0 };
			element->Centroid(cnt);
			std::cout << " element " << ElementTypeName(element->GetElementType()) << ":" << element->LocalID();
			std::cout << " procs";
			for(INMOST_DATA_ENUM_TYPE k = 0; k < procs.size(); ++k) std::cout << " " << procs[k];
			std::cout << " remote";
			for(INMOST_DATA_ENUM_TYPE k = 0; k < size; ++k) std::cout << " " << ((const INMOST_DATA_INTEGER_TYPE *)data)[k];
			std::cout << " cnt " << cnt[0] << " " << cnt[1] << " " << cnt[2];
			if (element.GetMeshLink()->HaveGlobalID(element->GetElementType()))
				std::cout << " gid " << element->GlobalID();
			std::cout << (element->New() ? " new" : " old");
			std::cout << (element->Hidden() ? " hidden" : "");
			std::cout << std::endl;
			std::cout.flush();
#if defined(USE_PARALLEL_WRITE_TIME)	
			Mesh* m = element.GetMeshLink();
			m->WriteTab(m->GetStream()) << " element " << ElementTypeName(element->GetElementType()) << ":" << element->LocalID();
			m->WriteTab(m->GetStream()) << " procs";
			for (INMOST_DATA_ENUM_TYPE k = 0; k < procs.size(); ++k) m->WriteTab(m->GetStream()) << " " << procs[k];
			m->WriteTab(m->GetStream()) << " remote";
			for (INMOST_DATA_ENUM_TYPE k = 0; k < size; ++k) m->WriteTab(m->GetStream()) << " " << ((const INMOST_DATA_INTEGER_TYPE*)data)[k];
			m->WriteTab(m->GetStream()) << " cnt " << cnt[0] << " " << cnt[1] << " " << cnt[2];
			if (element.GetMeshLink()->HaveGlobalID(element->GetElementType()))
				m->WriteTab(m->GetStream()) << " gid " << element->GlobalID();
			m->WriteTab(m->GetStream()) << (element->New() ? " new" : " old");
			m->WriteTab(m->GetStream()) << (element->Hidden() ? " hidden" : "");
			m->WriteTab(m->GetStream()) << std::endl;
			m->GetStream().flush();
#endif
			crash_status = 1;
			//exit(-1);
		}
	}

	void Mesh::CheckProcessors()
	{
#if !defined(NDEBUG)
		ENTER_FUNC();
		int bad = 0;
		for (Mesh::iteratorElement it = BeginElement(ESET | CELL | FACE | EDGE | NODE); it != EndElement(); ++it)
		{
			integer_array procs = it->IntegerArray(tag_processors);
			if (!procs.empty())
			{
				integer_array::iterator first, next;
				first = next = procs.begin();
				while (++next != procs.end()) 
				{
					if (*next <= *first)
					{
						bad++;
						break;
					}
					++first;
				}
			}
		}
			
		if (NewMarker()) //add new/old status as the last processor
		{
#if defined(USE_PARALLEL_STORAGE)
			for (Mesh::iteratorElement it = BeginElement(ESET | CELL | FACE | EDGE | NODE); it != EndElement(); ++it)
				it->IntegerArray(tag_processors).push_back(it->New());
#endif
		}
		
		ReduceData(tag_processors,ESET|CELL|FACE|EDGE|NODE,0,UnpackCheckProcessors);
		
		if (NewMarker()) //rem new/old status as the last processor
		{
#if defined(USE_PARALLEL_STORAGE)
			for (Mesh::iteratorElement it = BeginElement(ESET | CELL | FACE | EDGE | NODE); it != EndElement(); ++it)
				it->IntegerArray(tag_processors).pop_back();
#endif
		}
		
		Storage::integer crash = 0;
		if (bad)
		{
			std::cout << "processor " << GetProcessorRank() << " bad processors order " << bad << std::endl;
			REPORT_STR("processor " << GetProcessorRank() << " bad processors order " << bad);
			crash++;
		}
		crash += crash_status;
		crash = Integrate(crash);
		if (crash)
		{
			std::cout << "crash " << __FUNCTION__ << std::endl;
			REPORT_STR("crash from " << __FUNCTION__);
			std::cout.flush();
			exit(-1);
		}
		EXIT_FUNC();
#endif
	}

	/*
	template<typename Pred>
	void mergeSortRecursive(HandleType * v, Storage::enumerator left, Storage::enumerator right, const Pred & Pr)
	{
		if (left < right) 
		{
			if (right - left >= 32)
			{
				Storage::enumerator mid = (left + right) / 2;
				mergeSortRecursive(v, left, mid,Pr);
				mergeSortRecursive(v, mid + 1, right,Pr);
				std::inplace_merge(v + left, v + mid + 1, v + right + 1, Pr);
			}
			else std::sort(v + left, v + right + 1, Pr);
		}
	}

	template<typename Pred>
	void mergeSort(HandleType* v, Storage::enumerator size, const Pred & Pr)
	{
		mergeSortRecursive(v, 0, size - 1, Pr);
	}
	*/

	void Mesh::SortParallelStorage(parallel_storage & ghost, parallel_storage & shared, ElementType mask)
	{
		ENTER_FUNC();
#if defined(USE_MPI)
		REPORT_STR("sort shared elements")
		ENTER_BLOCK();
#if defined(USE_PARALLEL_WRITE_TIME)
		for (int j = 0; j < 5 * (int)shared.size(); ++j)
		{
			int i = j % 5;
			if (mask & ElementTypeFromDim(i))
			{
				parallel_storage::iterator it = shared.begin();
				int k = j / 5;
				while (k) ++it, --k;
				REPORT_VAL("processor", it->first);
				REPORT_STR(ElementTypeName(ElementTypeFromDim(i)) << " i " << i << " have global id " << (!tag_global_id.isValid() ? "invalid" : (tag_global_id.isDefined(ElementTypeFromDim(i)) ? "YES" : "NO")));
				REPORT_VAL(ElementTypeName(mask & ElementTypeFromDim(i)), it->second[i].size());
			}
		}
#endif
		
//#if defined(USE_OMP)
//#pragma omp parallel for schedule(dynamic,1)
//#endif
		//for(int j = 0; j < 5 * (int)shared.size(); ++j)
		for (parallel_storage::iterator it = shared.begin(); it != shared.end(); it++)
		{
			//int i = j % 5;
			for(int i = 0; i < 5; ++i ) if( mask & ElementTypeFromDim(i) )
			{
				//parallel_storage::iterator it = shared.begin();
				//int k = j / 5;
				//while (k) ++it, --k;
				//REPORT_VAL("processor", it->first);
				//ENTER_BLOCK();
				//REPORT_STR(ElementTypeName(ElementTypeFromDim(i)) << " have global id " << (!tag_global_id.isValid() ? "invalid" : (tag_global_id.isDefined(ElementTypeFromDim(i))?"YES":"NO")));
				if( !it->second[i].empty() )
				{
					if (HaveGlobalID(ElementTypeFromDim(i)))
					{
						std::sort(it->second[i].begin(), it->second[i].end(), GlobalIDComparator(this));
#if !defined(NDEBUG)
						integer bad = 0;
						element_set::iterator next = it->second[i].begin(), cur = next++;
						while(next < it->second[i].end())
						{
							if (GlobalID(*next) == GlobalID(*cur)) bad++;
							cur = next++;
						}
						if (bad) std::cout << __FILE__ << ":" << __LINE__ << " " << GetProcessorRank() << " duplicate global id " << bad << std::endl;
#endif
					}
					else
					{
						if (i < 4)
						{
							CentroidComparator cmp(this);
							std::sort(it->second[i].begin(), it->second[i].end(), cmp);
#if !defined(NDEBUG)
							integer bad = 0;
							element_set::iterator next = it->second[i].begin(), cur = next++;
							while (next < it->second[i].end())
							{
								if (cmp.Compare(*next,*cur) == 0) bad++;
								cur = next++;
							}
							if (bad) std::cout << __FILE__ << ":" << __LINE__ << " " << GetProcessorRank() << " duplicate coords " << bad << std::endl;
#endif
						}
						else
						{
							SetNameComparator cmp(this);
							std::sort(it->second[4].begin(), it->second[4].end(), cmp);
#if !defined(NDEBUG)
							integer bad = 0;
							element_set::iterator next = it->second[i].begin(), cur = next++;
							while (next < it->second[i].end())
							{
								if (cmp.Compare(*next, *cur) == 0) bad++;
								cur = next++;
							}
							if (bad) std::cout << __FILE__ << ":" << __LINE__ << " " << GetProcessorRank() << " duplicate set names " << bad << std::endl;
#endif
						}
					}

				}
				//REPORT_VAL(ElementTypeName(mask & ElementTypeFromDim(i)),it->second[i].size());
				//EXIT_BLOCK();
			}
		}
		EXIT_BLOCK();
		REPORT_STR("sort ghost elements")
		ENTER_BLOCK();
#if defined(USE_PARALLEL_WRITE_TIME)
		for (int j = 0; j < 5 * (int)ghost.size(); ++j)
		{
			int i = j % 5;
			if (mask & ElementTypeFromDim(i))
			{
				parallel_storage::iterator it = ghost.begin();
				int k = j / 5;
				while (k) ++it, --k;
				REPORT_VAL("processor", it->first);
				REPORT_STR(ElementTypeName(ElementTypeFromDim(i)) << " i " << i << " have global id " << (!tag_global_id.isValid() ? "invalid" : (tag_global_id.isDefined(ElementTypeFromDim(i))?"YES":"NO")));
				REPORT_VAL(ElementTypeName(mask & ElementTypeFromDim(i)), it->second[i].size());
			}
		}
#endif
		//for(parallel_storage::iterator it = ghost.begin(); it != ghost.end(); it++)
//#if defined(USE_OMP)
//#pragma omp parallel for schedule(dynamic,1)
//#endif
//		for (int j = 0; j < 5 * (int)ghost.size(); ++j)
		for (parallel_storage::iterator it = ghost.begin(); it != ghost.end(); it++)
		{
			//int i = j % 5;
			for(int i = 0; i < 5; ++i) if( mask & ElementTypeFromDim(i) )
			{
				//parallel_storage::iterator it = ghost.begin();
				//int k = j / 5;
				//while (k) ++it, --k;
				//REPORT_VAL("processor", it->first);
				//ENTER_BLOCK();
				//REPORT_STR(ElementTypeName(ElementTypeFromDim(i)) << " have global id " << (!tag_global_id.isValid() ? "invalid" : (tag_global_id.isDefined(ElementTypeFromDim(i))?"YES":"NO")));
				if( !it->second[i].empty() )
				{
					if (HaveGlobalID(ElementTypeFromDim(i)))
					{
						std::sort(it->second[i].begin(), it->second[i].end(), GlobalIDComparator(this));
#if !defined(NDEBUG)
						integer bad = 0;
						element_set::iterator next = it->second[i].begin(), cur = next++;
						while (next < it->second[i].end())
						{
							if (GlobalID(*next) == GlobalID(*cur)) bad++;
							cur = next++;
						}
						if (bad) std::cout << __FILE__ << ":" << __LINE__ << " " << GetProcessorRank() << " duplicate global id " << bad << std::endl;
#endif
					}
					else
					{
						if (i < 4)
						{
							CentroidComparator cmp(this);
							std::sort(it->second[i].begin(), it->second[i].end(), cmp);
#if !defined(NDEBUG)
							integer bad = 0;
							element_set::iterator next = it->second[i].begin(), cur = next++;
							while (next < it->second[i].end())
							{
								if (cmp.Compare(*next, *cur) == 0) bad++;
								cur = next++;
							}
							if (bad) std::cout << __FILE__ << ":" << __LINE__ << " " << GetProcessorRank() << " duplicate coords " << bad << std::endl;
#endif
						}
						else
						{
							SetNameComparator cmp(this);
							std::sort(it->second[4].begin(), it->second[4].end(), cmp);
#if !defined(NDEBUG)
							integer bad = 0;
							element_set::iterator next = it->second[i].begin(), cur = next++;
							while (next < it->second[i].end())
							{
								if (cmp.Compare(*next, *cur) == 0) bad++;
								cur = next++;
							}
							if (bad) std::cout << __FILE__ << ":" << __LINE__ << " " << GetProcessorRank() << " duplicate set names " << bad << std::endl;
#endif
						}
					}
				}
				//REPORT_VAL(ElementTypeName(mask & ElementTypeFromDim(i)),it->second[i].size());
				//EXIT_BLOCK();
			}
		}
		EXIT_BLOCK();
#else //USE_MPI and USE_PARALLEL_STORAGE
		(void) ghost; (void) shared; (void) mask;
#endif //USE_MPI and USE_PARALLEL_STORAGE
		EXIT_FUNC();
	}

	void Mesh::ReportParallelStorage()
	{
#if !defined(NDEBUG)
		ENTER_FUNC();
#if defined(USE_PARALLEL_STORAGE)
		REPORT_STR("shared");
		for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
		{
			REPORT_VAL("processor",it->first);
			for(int i = 0; i < 4; i++) //if( mask & ElementTypeFromDim(i) )
			{
				REPORT_STR(ElementTypeName(ElementTypeFromDim(i)) << " have global id " << (!tag_global_id.isValid() ? "invalid" : (tag_global_id.isDefined(ElementTypeFromDim(i))?"YES":"NO")));
				REPORT_VAL(ElementTypeName(ElementTypeFromDim(i)),it->second[i].size());
			}

            // ESET sort
			if (ElementTypeFromDim(4))
            {
				REPORT_VAL(ElementTypeName(ElementTypeFromDim(4)),it->second[4].size());
            }
		}
		REPORT_STR("ghost");
		for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)
		{
			REPORT_VAL("processor",it->first);
			for(int i = 0; i < 4; i++) //if( mask & ElementTypeFromDim(i) )
			{
				REPORT_STR(ElementTypeName(ElementTypeFromDim(i)) << " have global id " << (!tag_global_id.isValid() ? "invalid" : (tag_global_id.isDefined(ElementTypeFromDim(i))?"YES":"NO")));
				REPORT_VAL(ElementTypeName(ElementTypeFromDim(i)),it->second[i].size());
			}

            // ESET sort
			if (ElementTypeFromDim(4))
            {
				REPORT_VAL(ElementTypeName(ElementTypeFromDim(4)),it->second[4].size());
            }
		}
#endif //USE_PARALLEL_STORAGE
		EXIT_FUNC();
#endif //NDEBUG
	}
	
	void Mesh::ExchangeGhost(Storage::integer layers, ElementType bridge, MarkerType marker, bool delete_ghost)
	{
		//printf("%d called exchange ghost with layers %d bridge %s\n",GetProcessorRank(), layers,ElementTypeName(bridge));
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
		int mpirank = GetProcessorRank();
		// bool delete_ghost = false;
		//if( layers == Integer(tag_layers) && bridge == Integer(tag_bridge) ) return;
		//cout << "Check " << layers << " " << Integer(GetHandle(),tag_layers) << endl;
		//~ if( layers <= Integer(GetHandle(),tag_layers) ) delete_ghost = true;
		//~ else if( layers == Integer(GetHandle(),tag_layers) && bridge < Integer(GetHandle(),tag_bridge) ) delete_ghost = true;
		//if (marker != 0)
		// delete_ghost = true;
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
		CheckSetLinks(__FILE__,__LINE__);
		CheckGhostSharedCount(__FILE__, __LINE__);
		CheckOwners(__FILE__,__LINE__);
		CheckGIDs(__FILE__, __LINE__);
		CheckProcessors();
		CheckCentroids(__FILE__, __LINE__);
		//Storage::integer_array procs = IntegerArrayDV(GetHandle(),tag_processors);
		proc_elements old_layers;
		proc_elements current_layers;
		element_set all_visited;
		ENTER_BLOCK();
		{
			proc_elements shared_skin = ComputeSharedSkinSet(bridge, marker);
			//printf("%d shared skin size %d\n",GetProcessorRank(),shared_skin.size());
			time = Timer();
			for(proc_elements::iterator p = shared_skin.begin(); p != shared_skin.end(); p++)
			{
				MarkerType busy = CreateMarker();
				all_visited.clear();
				for(element_set::iterator it = p->second.begin(); it != p->second.end(); it++)
				{
					ElementArray<Element> adj = Element(this,*it)->getAdjElements(CELL);
					for(ElementArray<Element>::iterator jt = adj.begin(); jt != adj.end(); jt++)
						if( !jt->GetMarker(busy) )
						{
							if( jt->IntegerDF(tag_owner) != p->first )
							//if( jt->IntegerDF(tag_owner) == mpirank )
							{
								current_layers[p->first].push_back(*jt);
								if( layers > 0 )
								{
									Storage::integer_array adj_procs = jt->IntegerArrayDV(tag_processors);
									if( !std::binary_search(adj_procs.begin(),adj_procs.end(),p->first) )
									{
										assert(p->first >= 0 && p->first <= GetProcessorsNumber());
										jt->IntegerArray(tag_sendto).push_back(p->first);
									}
									if( delete_ghost ) jt->IntegerArray(layers_marker).push_back(p->first);
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
		EXIT_BLOCK();
		for(Storage::integer k = layers-1; k >= 0; k--)
		{
			ENTER_BLOCK();
			CheckSetLinks(__FILE__,__LINE__);
			CheckGhostSharedCount(__FILE__,__LINE__);
			CheckOwners(__FILE__, __LINE__);
			CheckGIDs(__FILE__, __LINE__);
			CheckProcessors();
			CheckCentroids(__FILE__, __LINE__);
			ExchangeMarked();
			old_layers.swap(current_layers);
			current_layers.clear();
			EXIT_BLOCK();
			if( k > 0 ) 
			{
				time = Timer();
				//!!!
				ENTER_BLOCK();
				for(proc_elements::iterator p = old_layers.begin(); p != old_layers.end(); p++)
				{
					element_set & ref_cur = current_layers[p->first];
					element_set & ref_old = old_layers[p->first];
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
										if( jt->IntegerDF(tag_owner) != p->first )
										{
											ref_cur.push_back(*kt);
											Storage::integer_array adj_procs = kt->IntegerArrayDV(tag_processors);
											if( !std::binary_search(adj_procs.begin(),adj_procs.end(),p->first) )
											{
												assert(p->first >= 0 && p->first <= GetProcessorsNumber());
												kt->IntegerArray(tag_sendto).push_back(p->first);
											}
											if( delete_ghost ) kt->IntegerArray(layers_marker).push_back(p->first);
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
				EXIT_BLOCK();
				time = Timer() - time;
				REPORT_STR("Mark layer");
				REPORT_VAL("layer",k);
				REPORT_VAL("time",time);
			}
		}
		if( delete_ghost )
		{
			ENTER_BLOCK();
			CheckSetLinks(__FILE__,__LINE__);
			CheckGhostSharedCount(__FILE__, __LINE__);
			CheckOwners(__FILE__, __LINE__);
			CheckGIDs(__FILE__, __LINE__);
			CheckProcessors();
			CheckCentroids(__FILE__, __LINE__);
			EXIT_BLOCK();
			time = Timer();
			ReduceData(layers_marker,CELL,0,UnpackLayersMarker);
			ExchangeData(layers_marker,CELL,0);
			ElementArray<Element> del_ghost(this);
			ENTER_BLOCK();
			for(iteratorElement it = BeginElement(CELL); it != EndElement(); it++)
			{
				if( GetStatus(*it) == Element::Ghost )
				{
					if( marker && !it->GetMarker(marker) ) continue;
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
			EXIT_BLOCK();
			//TagInteger del = CreateTag("DELETE_GHOST",DATA_INTEGER,CELL,NONE,1);
			//for(ElementArray<Element>::iterator it = del_ghost.begin(); it != del_ghost.end(); ++it) del[*it] = 1;
			//Save("before_delete_ghost.pvtk");
			time = Timer() - time;
			REPORT_STR("Select ghost elements to remove");
			REPORT_VAL("time",time);
			RemoveGhostElements(del_ghost);
			CheckSetLinks(__FILE__,__LINE__);
			CheckGhostSharedCount(__FILE__, __LINE__);
			CheckOwners(__FILE__, __LINE__);
			CheckGIDs(__FILE__, __LINE__);
			CheckProcessors();
			CheckCentroids(__FILE__, __LINE__);
			//Save("after_delete_ghost.pvtk");
			//exit(-1);
		}
		layers_marker = DeleteTag(layers_marker);
		//throw NotImplemented;
#else
		(void) layers;
		(void) bridge;
#endif
		EXIT_FUNC();

	}
	
	std::set<Storage::integer> CollectSetProcessors(ElementSet root, TagIntegerArray tag_processors)
	{
		Storage::integer_array parr;
		std::set<Storage::integer> ret_procs;
		if( root.HaveChild() )
		{
			for(ElementSet it = root.GetChild(); it != InvalidElementSet(); it = it.GetSibling())
			{
				std::set<Storage::integer> fill_procs = CollectSetProcessors(it,tag_processors);
				ret_procs.insert(fill_procs.begin(),fill_procs.end());
			}
			parr = root.IntegerArray(tag_processors);
			ret_procs.insert(parr.begin(),parr.end());
			parr.replace(parr.begin(),parr.end(),ret_procs.begin(),ret_procs.end());
		}
		else
		{
			parr = root.IntegerArray(tag_processors);
			ret_procs.insert(parr.begin(),parr.end());
		}
		return ret_procs;
	}
	/*
	void ReportSetProcs(ElementSet root, TagIntegerArray tag_processors, int tab, std::fstream & fout)
	{
		if( root.HaveChild() )
		{
			for(ElementSet it = root.GetChild(); it != InvalidElementSet(); it = it.GetSibling())
				ReportSetProcs(it,tag_processors,tab+1,fout);
		}
		Storage::integer_array parr = root.IntegerArray(tag_processors);
		for(int q = 0; q < tab; ++q) fout << "  ";
		fout << "set " << root.GetName() << " procs ";
		for(INMOST_DATA_ENUM_TYPE q = 0; q < parr.size(); ++q) fout << parr[q] << " ";
		fout << " owner " << root.Integer(root.GetMeshLink()->OwnerTag());
		fout << " kids " << root.Size();
		fout << " status " << Element::StatusName(root.GetStatus()) << std::endl;
	}
	*/
	void Mesh::Redistribute()
	{
		if( m_state == Serial ) return;
		ENTER_FUNC();
#if defined(USE_MPI)
		int mpirank = GetProcessorRank();
		Tag tag_new_owner = CreateTag("TEMPORARY_NEW_OWNER",DATA_INTEGER,ESET | FACE | EDGE | NODE, NONE,1);
		if( !tag_new_owner.isDefined(CELL) ) throw DistributionTagWasNotFilled;
		Tag tag_new_processors = CreateTag("TEMPORARY_NEW_PROCESSORS",DATA_INTEGER,ESET | CELL | FACE | EDGE | NODE, NONE);
		ElementType bridge = Integer(GetHandle(),tag_bridge);
		Storage::integer layers = Integer(GetHandle(),tag_layers);
		
		ExchangeData(tag_new_owner,CELL,0);
		
		
		ENTER_BLOCK();
		REPORT_STR("Determine new processors");
		//for(iteratorElement it = BeginElement(CELL); it != EndElement(); it++)
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for(integer i = 0; i < CellLastLocalID(); ++i) if( isValidCell(i))
			CellByLocalID(i)->IntegerArrayDV(tag_new_processors).push_back(CellByLocalID(i)->IntegerDF(tag_new_owner));
		//determine tag_new_processors for FACEs to calculate shared skin
		for(ElementType mask = FACE; mask >= NODE; mask = PrevElementType(mask))
		{
			REPORT_VAL("for ", ElementTypeName(mask));
			ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				std::vector<integer> result, intersection;
#if defined(USE_OMP)
#pragma omp for
#endif
				for (integer i = 0; i < LastLocalID(mask); ++i) if (isValidElement(mask, i))
				{
					Element it = ElementByLocalID(mask, i);
					Storage::integer_array procs = it->IntegerArrayDV(tag_new_processors);
					determine_my_procs_high(this, it->GetHandle(), tag_new_processors, result, intersection);
					if (!result.empty())
					{
						it->IntegerDF(tag_new_owner) = result[0];
						procs.replace(procs.begin(), procs.end(), result.begin(), result.end());
					}
					else
					{
						it->IntegerDF(tag_new_owner) = mpirank;
						procs.clear();
						procs.push_back(mpirank);
					}
				}
			}
			EXIT_BLOCK();
		}	
		EXIT_BLOCK();
		ExchangeData(tag_new_owner,FACE | EDGE | NODE,0);
		ENTER_BLOCK();
		REPORT_STR("Determine new processors for layers");
		if( bridge != NONE && layers != 0 )
		{
			REPORT_STR("Multiple layers");
			proc_elements redistribute_skin, skin_faces;
			element_set all_visited;
		
			ReduceData(tag_new_processors,FACE,0,RedistUnpack);
			ExchangeData(tag_new_processors,FACE,0);
			
			ENTER_BLOCK();
			REPORT_STR("Compute shared skin");
			//compute skin around migrated cells to compute ghost layer
			//for(iteratorElement it = BeginElement(FACE); it != EndElement(); it++)
			for(integer i = 0; i < FaceLastLocalID(); ++i) if( isValidFace(i) )
			{
				Face it = FaceByLocalID(i);
				Storage::integer_array procs = it->IntegerArrayDV(tag_new_processors);
				if( procs.size() == 2 ) //there may be only 2 procs for every face max, because we have unique redistribution for every cell at the beginning
				{
					skin_faces[procs[0]].push_back(it->GetHandle());
					skin_faces[procs[1]].push_back(it->GetHandle());
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
			EXIT_BLOCK();
			
			ENTER_BLOCK();
			REPORT_STR("Mark all layers");
			//now mark all cell layers 
			{
				proc_elements current_layers,old_layers;
				ENTER_BLOCK();
				REPORT_STR("Mark first layer");
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
				EXIT_BLOCK();
				
				for(Storage::integer k = layers-1; k > 0; k--)
				{
					ENTER_BLOCK();
					REPORT_STR("Mark layer");
					REPORT_VAL("layer",k);
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
					EXIT_BLOCK();
				}
			}
			EXIT_BLOCK();
			ReduceData(tag_new_processors,CELL,0,RedistUnpack);
			ExchangeData(tag_new_processors,CELL,0);
			
			ENTER_BLOCK();
			REPORT_STR("Extend processors list to lower level adjacencies");
			for(ElementType mask = FACE; mask >= NODE; mask = PrevElementType(mask))
			{
				//for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
#if defined(USE_OMP)
#pragma omp parallel
#endif
				{
					std::vector<integer> result, intersection;
#if defined(USE_OMP)
#pragma omp for
#endif
					for (integer i = 0; i < LastLocalID(mask); ++i) if (isValidElement(mask,i))
					{
						Element it = ElementByLocalID(mask, i);
						Storage::integer_array procs = it->IntegerArrayDV(tag_new_processors);
						determine_my_procs_high(this, it->GetHandle(), tag_new_processors, result, intersection);
						if (result.empty())
						{
							procs.clear();
							procs.push_back(mpirank);
						}
						else procs.replace(procs.begin(), procs.end(), result.begin(), result.end());
					}
				}
			}
			
			EXIT_BLOCK();
			ReduceData(tag_new_processors,FACE| EDGE| NODE,0,RedistUnpack);
			ExchangeData(tag_new_processors,FACE|EDGE|NODE,0);
			REPORT_STR("Detect processors");
		}
		else 
		{
			REPORT_STR("No layers");
			ReduceData(tag_new_processors,FACE | EDGE | NODE,0,RedistUnpack);
			ExchangeData(tag_new_processors,FACE | EDGE | NODE,0);
		}
		EXIT_BLOCK();

		MarkerType reffered = CreateMarker();
		ENTER_BLOCK();
		REPORT_STR("Determine entities to store links to elements");
		for(iteratorTag t = BeginTag(); t != EndTag(); t++)
		{
			if( t->GetTagName().substr(0,9) == "PROTECTED" ) continue;
			if( t->GetDataType() == DATA_REFERENCE )
			{
				for (ElementType etype = NODE; etype <= ESET; etype = NextElementType(etype)) if (t->isDefined(etype))
				{
					for (integer eit = 0; eit < LastLocalID(etype); ++eit) if (isValidElement(etype, eit))
					{
						Element it = ElementByLocalID(etype, eit);
						if (it->HaveData(*t))
						{
							Storage::integer owner = it->Integer(tag_new_owner);
							Storage::integer_array procs;
							Storage::reference_array refs = it->ReferenceArray(*t);
							for (Storage::reference_array::iterator jt = refs.begin(); jt != refs.end(); ++jt) if (jt->isValid())
							{
								procs = jt->IntegerArrayDV(tag_new_processors);
								Storage::integer_array::iterator find = std::lower_bound(procs.begin(), procs.end(), owner);
								if (find == procs.end() || *find != owner) procs.insert(find, owner);
								jt->SetMarker(reffered);
							}
						}
					}
				}
			}
		}
		ReduceData(tag_new_processors,ESET|CELL|FACE|EDGE|NODE,0,RedistUnpack);
		ExchangeData(tag_new_processors,ESET|CELL|FACE|EDGE|NODE,0);
		EXIT_BLOCK();

		ENTER_BLOCK();
		REPORT_STR("Extend processors list to lower level adjacencies");
		//OPTIMIZE!!!
		for(ElementType mask = FACE; mask >= NODE; mask = PrevElementType(mask))
		{
			//for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				std::vector<integer> result, intersection;
#if defined(USE_OMP)
#pragma omp for
#endif
				for (integer i = 0; i < LastLocalID(mask); ++i) if (isValidElement(mask, i))
				{
					Element it = ElementByLocalID(mask, i);
					Storage::integer_array procs = it->IntegerArrayDV(tag_new_processors);
					determine_my_procs_high(this, it->GetHandle(), tag_new_processors, result, intersection);
					if (result.empty())
					{
						procs.clear();
						procs.push_back(mpirank);
					}
					else procs.replace(procs.begin(), procs.end(), result.begin(), result.end());
				}
			}
		}
		ReduceData(tag_new_processors,FACE| EDGE| NODE,0,RedistUnpack);
		ExchangeData(tag_new_processors,FACE|EDGE|NODE,0);
		EXIT_BLOCK();

		ENTER_BLOCK();
		REPORT_STR("Determine new processors for sets");
		ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for(integer eit = 0; eit < EsetLastLocalID(); ++eit) if( isValidElementSet(eit) )
		{
			ElementSet it = EsetByLocalID(eit);
			if( !it.Empty() )
			{
				std::set<Storage::integer> elem_procs;
				for(ElementSet::iterator jt = it.Begin(); jt != it.End(); ++jt)
				{
					Storage::integer_array new_procs = jt->IntegerArrayDV(tag_new_processors);
					elem_procs.insert(new_procs.begin(),new_procs.end());
				}
				Storage::integer_array set_procs = it->IntegerArrayDV(tag_new_processors);
				elem_procs.insert(set_procs.begin(),set_procs.end());
				set_procs.replace(set_procs.begin(),set_procs.end(),elem_procs.begin(),elem_procs.end());
			}
		}
		EXIT_BLOCK();
		
		ENTER_BLOCK();
		
		for(integer eit = 0; eit < EsetLastLocalID(); ++eit) if( isValidElementSet(eit) )
		{
			ElementSet it = EsetByLocalID(eit);
			if( !it.HaveParent() )  CollectSetProcessors(it,tag_new_processors);
		}
		EXIT_BLOCK();
		ReduceData(tag_new_processors,ESET,0,RedistUnpack);
		ExchangeData(tag_new_processors,ESET,0);

		for(integer eit = 0; eit < EsetLastLocalID(); ++eit) if( isValidElementSet(eit) )
		{
			ElementSet it = EsetByLocalID(eit);
			Storage::integer_array new_procs = it->IntegerArray(tag_new_processors);
			if( !new_procs.empty() ) 
				it->IntegerDF(tag_new_owner) = new_procs[0];
		}
		
		ExchangeData(tag_new_owner,ESET,0);
		/*
		{
			std::fstream fout("file"+std::to_string(GetProcessorRank())+".txt",std::ios::out);
			for(INMOST_DATA_ENUM_TYPE eit = 0; eit < EsetLastLocalID(); ++eit) if( isValidElementSet(eit) )
			{
				ElementSet it = EsetByLocalID(eit);
				if( !it.HaveParent() )  ReportSetProcs(it,tag_new_processors,0,fout);
			}
		}
		*/
		EXIT_BLOCK();

		
		
		
		ENTER_BLOCK();
		REPORT_STR("Determine local entities to send");
		for(ElementType etype = NODE; etype <= ESET; etype = NextElementType(etype))
			for (integer eit = 0; eit < LastLocalID(etype); ++eit) if (isValidElement(etype, eit))
			{
				Element it = ElementByLocalID(etype, eit);
				//Storage::integer_array new_procs = it->IntegerArray(tag_new_processors);
				if (it->GetMarker(reffered))
				{
					Storage::integer_array sendto = it->IntegerArray(SendtoTag());
					Storage::integer_array procs = it->IntegerArray(tag_new_processors);
					sendto.replace(sendto.begin(), sendto.end(), procs.begin(), procs.end());
				}
				else if (etype == ESET || it->IntegerDF(tag_owner) == mpirank) // deal with my entities, others will deal with theirs
					it->SendTo(it->IntegerArray(tag_new_processors));
				/*
				{
					//compute processors that should have the entity but they have not
					Storage::integer_array old_procs = it->IntegerArrayDV(tag_processors);
					result.resize(new_procs.size());
					std::vector<Storage::integer>::iterator end = std::set_difference(new_procs.begin(),new_procs.end(),old_procs.begin(),old_procs.end(),result.begin());
					result.resize(end-result.begin());
					//mark to send entity to processors that don't have it
					Storage::integer_array sendto = it->IntegerArray(tag_sendto);
					sendto.insert(sendto.end(),result.begin(),result.end());
				}
				*/
			}
		EXIT_BLOCK();

		ReleaseMarker(reffered,ESET|CELL|FACE|EDGE|NODE);
		/*
		{
			std::fstream fout("sendto"+std::to_string(GetProcessorRank())+".txt",std::ios::out);
			for(INMOST_DATA_ENUM_TYPE eit = 0; eit < EsetLastLocalID(); ++eit) if( isValidElementSet(eit) )
			{
				ElementSet it = EsetByLocalID(eit);
				if( !it.HaveParent() )  ReportSetProcs(it,SendtoTag(),0,fout);
			}
		}*/

		
		ENTER_BLOCK();
		REPORT_STR("Migrate elements");
		ExchangeMarked(AMigrate);
		EXIT_BLOCK();
		
		ENTER_BLOCK();
		REPORT_STR("Delete auxilarry tags");
		DeleteTag(tag_new_owner);
		DeleteTag(tag_new_processors);
		EXIT_BLOCK();


		CheckGhostSharedCount(__FILE__, __LINE__);
		CheckOwners(__FILE__, __LINE__);
		CheckGIDs(__FILE__, __LINE__);
		CheckProcessors();
		CheckCentroids(__FILE__, __LINE__);

		
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
		for (iteratorFace it = BeginFace(); it != EndFace(); ++it) it->RemMarker(mrk);
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
		ENTER_FUNC();
#ifdef USE_MPI
		if( tag_global_id.isValid() ) tag_global_id = DeleteTag(tag_global_id,ESET);
        int mpirank = GetProcessorRank();
        int mpisize = GetProcessorsNumber();

        std::map<std::string, std::vector<int> > map_names; // key - set_name, value - array of processors ranks which has this set

		std::vector<char> set_names_send;
		std::vector<char> set_names_recv;
		std::vector<int> size_recv(mpisize);
		int size_send, size_recv_all = 0, prealloc = 0;
		ENTER_BLOCK();
		for(std::map<std::string,HandleType>::iterator it = set_search.begin(); it != set_search.end(); ++it) if( !Hidden(it->second) )
			prealloc += (int)it->first.size()+1;
		set_names_send.reserve(prealloc);
		for(std::map<std::string,HandleType>::iterator it = set_search.begin(); it != set_search.end(); ++it) if( !Hidden(it->second) )
		{
			set_names_send.insert(set_names_send.end(),it->first.begin(),it->first.end());
			set_names_send.push_back('\0');
		}
		size_send = (int)set_names_send.size();
		EXIT_BLOCK();
		ENTER_BLOCK();
		REPORT_MPI(MPI_Allgather(&size_send,1,MPI_INT,&size_recv[0],1,MPI_INT,comm));
		EXIT_BLOCK();
		
		std::vector<int> displs((size_t)mpisize+1,0);
		ENTER_BLOCK();
		for(int k = 0; k < mpisize; k++)
		{
			size_recv_all += size_recv[k];
			displs[(size_t)k+1] = displs[k]+size_recv[k];
		}
		EXIT_BLOCK();
		set_names_recv.resize(size_recv_all);
		char * send_ptr = set_names_send.empty() ? NULL : &set_names_send[0];
		char * recv_ptr = set_names_recv.empty() ? NULL : &set_names_recv[0];
		ENTER_BLOCK();
		REPORT_MPI(MPI_Allgatherv(send_ptr,size_send,MPI_CHAR,recv_ptr,&size_recv[0],&displs[0],MPI_CHAR,comm));
		EXIT_BLOCK();
		ENTER_BLOCK();
		for(int k = 0; k < mpisize; ++k)
		{
			char * beg_ptr = recv_ptr + displs[k];
			char * end_ptr = recv_ptr + displs[k] + size_recv[k];
			while( beg_ptr < end_ptr )
			{
				std::string name(beg_ptr);
//				REPORT_STR("\"" << name << "\" size " << name.size() << " proc " << k);
				map_names[name].push_back(k);
				beg_ptr += name.size()+1;
			}
		}
		EXIT_BLOCK();
		
		// Change status for self sets
		//for(Mesh::iteratorSet set = BeginSet(); set != EndSet(); ++set)
		/*
		ENTER_BLOCK();
		for(std::map<std::string,std::vector<int> >::iterator it = map_names.begin(); it != map_names.end(); ++it)
		{
			std::map<std::string,HandleType>::iterator find = set_search.find(it->first);
			REPORT_STR("\"" << it->first << "\" size " << it->first.size() << " procs " << it->second.size());
			if( find == set_search.end() )
			{REPORT_STR("handle not found");}
			else
			{
				REPORT_VAL("handle ", find->second);
				REPORT_VAL("status ",Element::StatusName(GetStatus(find->second)));
				REPORT_VAL("owner ",Integer(find->second,tag_owner));
			}
			for(std::vector<int>::iterator jt = it->second.begin(); jt != it->second.end(); ++jt) REPORT_STR("proc " << *jt);
		}
		EXIT_BLOCK();
		 */
		ENTER_BLOCK();
		REPORT_VAL("number of sets in set_search: ", set_search.size());
		REPORT_VAL("number of sets in mesh: ", NumberOfSets());
		REPORT_VAL("number of sets in map_names: ", map_names.size());
		for(std::map<std::string,HandleType>::iterator it = set_search.begin(); it != set_search.end(); ++it)  if( !Hidden(it->second) )
		{
			//REPORT_STR("\"" << it->first << "\" size " << it->first.size() << " handle " << it->second);
			std::vector<int> & procs = map_names[it->first];
			//if( procs.empty() )
			//	std::cout << "no procs for " << it->first << std::endl;
			//REPORT_VAL("procs ",procs.size());
			//std::string set_name = set->GetName();
			//REPORT_VAL("set name: ", set_name);
			Storage::integer_array arr = IntegerArrayDV(it->second,tag_processors);
			//if( !std::is_sorted(procs.begin(),procs.end()))
			//	std::cout << __FILE__ << ":" << __LINE__ << " not sorted!" << std::endl;
			//std::sort(map_names[set_name].begin(),map_names[set_name].end());
			//map_names[set_name].resize(std::unique(map_names[set_name].begin(),map_names[set_name].end())-map_names[set_name].begin());
			arr.resize((INMOST_DATA_ENUM_TYPE)procs.size());
			for (size_t i = 0; i < procs.size(); i++) arr[(INMOST_DATA_ENUM_TYPE)i] = procs[i];
			assert(procs.size() > 0);
			if (procs.size() == 1)
			{
				assert(procs[0] == mpirank);
				SetStatus(it->second, Element::Owned);
				IntegerDF(it->second, tag_owner) = mpirank;
			}
			else
			{
				//int min = map_names[set_name][0];
				//for (int i = 1; i < map_names[set_name].size(); i++)
				//	if (map_names[set_name][i] < min)
				//		min = map_names[set_name][i];
				IntegerDF(it->second,tag_owner) = procs[0];
				if (procs[0] == mpirank)
					SetStatus(it->second, Element::Shared);
				else
					SetStatus(it->second, Element::Ghost);
			}
			/*
			std::stringstream pstr;
			pstr << "name " << it->first;
			pstr << " handle " << it->second;
			pstr << " owner " << IntegerDF(it->second, tag_owner);
			pstr << " status " << Element::StatusName(GetStatus(it->second));
			pstr << " procs";
			if( Hidden(it->second) ) pstr << " hidden";
			else pstr << " not hidden";
			for (int i = 0; i < arr.size(); i++) pstr << " " << arr[i];
			REPORT_STR(pstr.str());*/
		}
		EXIT_BLOCK();
		//CheckGhostSharedCount(__FILE__,__LINE__,ESET);
		ENTER_BLOCK();
		RecomputeParallelStorage(ESET);
		EXIT_BLOCK();
		//ComputeSharedProcs();
		ENTER_BLOCK();
		CheckGhostSharedCount(__FILE__,__LINE__,ESET);
		EXIT_BLOCK();
		AssignGlobalID(ESET);
		ENTER_BLOCK();
		CheckGhostSharedCount(__FILE__,__LINE__,ESET);
		CheckProcsSorted(__FILE__,__LINE__);
		EXIT_BLOCK();
#endif
		EXIT_FUNC();
    }
}


#endif
