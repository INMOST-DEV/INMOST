#pragma once
#ifndef INMOST_MESH_H_INCLUDED
#define INMOST_MESH_H_INCLUDED

#include "inmost_common.h"

#if defined(USE_MESH)

#define __NDT 3
#define __NET 6

#define NEW_MARKERS
#define NEW_CONNECTIONS
#define NEW_SPARSE

namespace INMOST
{


	class Mesh;
	class Storage;
	class Element;
	class TagManager;
	class Node;
	class Edge;
	class Face;
	class Cell;

	typedef INMOST_DATA_BULK_TYPE ElementType;
	static const ElementType NONE = 0x00;
	static const ElementType NODE = 0x01;
	static const ElementType EDGE = 0x02;
	static const ElementType FACE = 0x04;
	static const ElementType CELL = 0x08;
	static const ElementType ESET = 0x10;
	static const ElementType MESH = 0x20;
	
	
	typedef INMOST_DATA_BULK_TYPE GeometricData;
	static const GeometricData CENTROID     = 0;
	static const GeometricData NORMAL       = 1;
	static const GeometricData ORIENTATION  = 2;
	static const GeometricData MEASURE      = 3;
	static const GeometricData BARYCENTER   = 4;
	
	typedef INMOST_DATA_BULK_TYPE SyncBitOp; //< This type is used for marker synchronization
	static const SyncBitOp SYNC_BIT_NEW = 0;
	static const SyncBitOp SYNC_BIT_OR  = 1;
	static const SyncBitOp SYNC_BIT_XOR = 2;
	static const SyncBitOp SYNC_BIT_AND = 3;

	//Use topolgy checking for debug purposes
	typedef INMOST_DATA_ENUM_TYPE TopologyCheck; 
	static const TopologyCheck THROW_EXCEPTION        = 0x00000001; //done//throw TopologyError exception on error
	static const TopologyCheck PRINT_NOTIFY           = 0x00000002; //done//print topology notify to std::cerr
	static const TopologyCheck DELETE_ON_ERROR        = 0x00000004; //done//element should be deleted if there is an error in EndTopologyCheck
	static const TopologyCheck MARK_ON_ERROR          = 0x00000008; //done//bitwise OR of errors is recorded to sparse integer tag of element (if not deleted), availible by TopologyErrorTag()
	static const TopologyCheck DUPLICATE_EDGE         = 0x00000010; //done//check that edge already exists, on occurance CreateEdge silently returns existant edge
	static const TopologyCheck DUPLICATE_FACE         = 0x00000020; //done//check that face already exists, on occurance CreateFace silently returns existant face
	static const TopologyCheck DUPLICATE_CELL         = 0x00000040; //done//check that cell already exists, on occurance CreateCell silently returns existant cell
	static const TopologyCheck DEGENERATE_EDGE        = 0x00000080; //done//check that edge have more then two nodes (no support for this kind of object)
	static const TopologyCheck DEGENERATE_FACE        = 0x00000100; //done//produce error if face consists of less then 3 edges
	static const TopologyCheck DEGENERATE_CELL        = 0x00000200; //done//produce error if cell consists of less then 4 faces
	static const TopologyCheck FACE_ORIENTATION       = 0x00000400; //done//produce error if face have wrong orientation of edges, for star-shaped elements only
	static const TopologyCheck FACE_PLANARITY         = 0x00000800; //done//produce error if face is non planar
	static const TopologyCheck INTERLEAVED_FACES      = 0x00001000; //done//produce error if there is another face that already uses same nodes
	static const TopologyCheck TRIPLE_SHARED_FACE     = 0x00002000; //done//check that every face have exectly two neighbours, if DUPLICATE_CELL is activated, then checks for cell duplication
	static const TopologyCheck FLATTENED_CELL         = 0x00004000; //done//produce error if one of the faces of the cell contains all the nodes of the cell
	static const TopologyCheck ADJACENT_DUPLICATE     = 0x00008000; //done//produce error if provided array of elements for creation contain duplications
	static const TopologyCheck ADJACENT_HIDDEN        = 0x00010000; //done//hidden elements should not be used when new elements are created
	static const TopologyCheck ADJACENT_ALIEN         = 0x00020000; //done//elements that belong to another mesh should not be used
	static const TopologyCheck ADJACENT_DIMENSION     = 0x00040000; //done//produce error if provided array of elements have wrong geometric dimension
	static const TopologyCheck PROHIBIT_MULTILINE     = 0x00080000; //done//(needs NEED_TEST_CLOSURE) produce error if edges of faces are not closed, or face is folded
	static const TopologyCheck PROHIBIT_POLYGON       = 0x00100000; //done//allow only known types of elements: Tet,Quad, Line
	static const TopologyCheck PROHIBIT_MULTIPOLYGON  = 0x00200000; //done//(needs NEED_TEST_CLOSURE) produce error if faces of cell are not closed, or cell is folded
	static const TopologyCheck PROHIBIT_POLYHEDRON    = 0x00400000; //done//allow only known types of elements: Tet,Hex,Prism,Pyramid
	static const TopologyCheck FACE_EDGES_ORDER       = 0x00800000; //not implemented(implement CheckEdgeOrder,FixEdgeOrder, class Cell,Face)//edges of the face should form one closed loop
	static const TopologyCheck PROHIBIT_CONCAVE_FACE  = 0x01000000; //not implemented//don't allow concave face
	static const TopologyCheck PROHIBIT_CONCAVE_CELL  = 0x02000000; //not implemented//don't allow concave cell
	static const TopologyCheck PROHIBIT_NONSTAR_FACE  = 0x04000000; //not implemented//don't allow non-star shaped face
	static const TopologyCheck PROHIBIT_NONSTAR_CELL  = 0x08000000; //not implemented//don't allow non-star shaped concave cell
	static const TopologyCheck FACE_SELF_INTERSECTION = 0x10000000; //not implemented//edges of the face don't cross each other
	static const TopologyCheck CELL_SELF_INTERSECTION = 0x20000000; //not implemented//faces of the cell don't cross each other
	static const TopologyCheck NEED_TEST_CLOSURE      = 0x40000000; //done//silent, test's for closure in ComputeGeometricType, needed to detect MultiLine and MultiPolygon
	static const TopologyCheck DISABLE_2D             = 0x80000000; //done//don't allow 2d grids, where edges appear to be vertexes, faces are edges and cells are faces
	static const TopologyCheck GRID_CONFORMITY        = NEED_TEST_CLOSURE | PROHIBIT_MULTILINE | PROHIBIT_MULTIPOLYGON  | INTERLEAVED_FACES | TRIPLE_SHARED_FACE;
	static const TopologyCheck DEFAULT_CHECK          = THROW_EXCEPTION | DUPLICATE_EDGE | DUPLICATE_FACE | PRINT_NOTIFY;
	
	__INLINE static const char * TopologyCheckNotifyString(TopologyCheck c)
	{
		switch(c)
		{
			case THROW_EXCEPTION:        return "exception thrown";
			case PRINT_NOTIFY:           return "print notify";
			case DELETE_ON_ERROR:        return "element should be deleted on error";
			case MARK_ON_ERROR:          return "element is marked on error";
			case DUPLICATE_EDGE:         return "search for duplicate edge";
			case DUPLICATE_FACE:         return "search for duplicate face";
			case DUPLICATE_CELL:         return "search for duplicate cell";
			case DEGENERATE_EDGE:        return "TOPOLOGY ERROR: curvilinear edge found"; 
			case DEGENERATE_FACE:        return "TOPOLOGY ERROR: degenerate face found";
			case DEGENERATE_CELL:        return "TOPOLOGY ERROR: degenerate cell found";
			case FACE_ORIENTATION:       return "TOPOLOGY ERROR: bad face orientation";
			case FACE_PLANARITY:         return "TOPOLOGY ERROR: non-planar face found";
			case INTERLEAVED_FACES:      return "TOPOLOGY ERROR: interleaving faces found";
			case TRIPLE_SHARED_FACE:     return "TOPOLOGY ERROR: face have more then two neighbours"; 
			case FLATTENED_CELL:         return "TOPOLOGY ERROR: flattened cell found"; 
			case ADJACENT_DUPLICATE:     return "TOPOLOGY ERROR: duplicates in adjacent elements";
			case ADJACENT_HIDDEN:        return "TOPOLOGY ERROR: hidden element is used as adjacent"; 
			case ADJACENT_ALIEN:         return "TOPOLOGY ERROR: element from other mesh is used as adjacent"; 
			case ADJACENT_DIMENSION:     return "TOPOLOGY ERROR: wrong dimension of adjacent elements";
			case PROHIBIT_MULTILINE:     return "TOPOLOGY ERROR: multiline is prohibited"; 
			case PROHIBIT_POLYGON:       return "TOPOLOGY ERROR: polygon is prohibited"; 
			case PROHIBIT_MULTIPOLYGON:  return "TOPOLOGY ERROR: multipolygon is prohibited"; 
			case PROHIBIT_POLYHEDRON:    return "TOPOLOGY ERROR: polyhedron is prohibited"; 
			case FACE_EDGES_ORDER:       return "TOPOLOGY ERROR: no order in face edges"; 
			case PROHIBIT_CONCAVE_FACE:  return "TOPOLOGY ERROR: concave faces are prohibited"; 
			case PROHIBIT_CONCAVE_CELL:  return "TOPOLOGY ERROR: concave cells are prohibited"; 
			case PROHIBIT_NONSTAR_FACE:  return "TOPOLOGY ERROR: non star-shaped faces are prohibited"; 
			case PROHIBIT_NONSTAR_CELL:  return "TOPOLOGY ERROR: non star-shpaed cells are prohibited"; 
			case FACE_SELF_INTERSECTION: return "TOPOLOGY ERROR: self intersection of face edges detected"; 
			case CELL_SELF_INTERSECTION: return "TOPOLOGY ERROR: self intersection of cell faces detected"; 
			case DISABLE_2D:             return "TOPOLOGY ERROR: 2d mesh support is disabled"; 
			default: return "unknown";
		}
	}
	
	/// This function helps to determine whether one type is chosen or multiple
	__INLINE static bool OneType(ElementType t) {return t > 0 && (t & (t-1)) == 0;}

	__INLINE static int ElementNum(ElementType t)
	{
		unsigned int v = static_cast<unsigned int>(t);  // 32-bit value to find the log2 of 
		static const unsigned int b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 0xFF00FF00, 0xFFFF0000};
		register unsigned int r = (v & b[0]) != 0;
		r |= ((v & b[4]) != 0) << 4;
		r |= ((v & b[3]) != 0) << 3;
		r |= ((v & b[2]) != 0) << 2;
		r |= ((v & b[1]) != 0) << 1;
		return static_cast<int>(r);
	}
	
	
	__INLINE static const char * ElementTypeName(ElementType t)
	{
		switch(t)
		{
			case NONE: return "NONE";
			case NODE: return "NODE";
			case EDGE: return "EDGE";
			case FACE: return "FACE";
			case CELL: return "CELL";
			case ESET: return "ESET";
			case MESH: return "MESH";
		}
		return "UNKNOWN";
	}
	
	
	
	typedef array<INMOST_DATA_REAL_TYPE>    inner_real_array;
	typedef array<INMOST_DATA_INTEGER_TYPE> inner_integer_array;
	typedef array<INMOST_DATA_BULK_TYPE>    inner_bulk_array;
	typedef array<Element *>                inner_reference_array;
		
	enum DataType
	{
		DATA_REAL      = 0, 
		DATA_INTEGER   = 1, 
		DATA_BULK      = 2,
		DATA_REFERENCE = 3
	};
	
	__INLINE static INMOST_DATA_ENUM_TYPE DataTypeBytesSize(DataType t)
	{
		switch(t)
		{
			case DATA_BULK:      return sizeof(INMOST_DATA_BULK_TYPE);
			case DATA_INTEGER:   return sizeof(INMOST_DATA_INTEGER_TYPE);
			case DATA_REAL:      return sizeof(INMOST_DATA_REAL_TYPE);
			case DATA_REFERENCE: return sizeof(Element *);
		}
		return 0;
	}
	
	__INLINE static INMOST_DATA_ENUM_TYPE VariableDataSize(DataType t)
	{
		switch(t)
		{
			case DATA_REAL:      return sizeof(inner_real_array);
			case DATA_INTEGER:   return sizeof(inner_integer_array);
			case DATA_BULK:      return sizeof(inner_bulk_array);
			case DATA_REFERENCE: return sizeof(inner_reference_array);
		}
		return 0;
	}
	
	__INLINE static const char * DataTypeName(DataType t)
	{
		switch(t)
		{
			case DATA_REAL:      return "REAL";
			case DATA_INTEGER:   return "INTEGER";
			case DATA_BULK:      return "BULK";
			case DATA_REFERENCE: return "REFERENCE";
		}
		return "UNKNOWN";
	}


	typedef INMOST_DATA_ENUM_TYPE MarkerType; // low 8 bits - marker mask, rest high bits - position of marker
	static const INMOST_DATA_ENUM_TYPE MarkerFields = 16;   // number of chars to hold all markers, total number (MarkerFields * bits_per_char)
	static const INMOST_DATA_BULK_TYPE MarkerMask   = static_cast<INMOST_DATA_BULK_TYPE>(-1); // bit mask to obtain marker mask within MarkerType
	static const INMOST_DATA_BULK_TYPE MarkerShift  = sizeof(INMOST_DATA_BULK_TYPE)*8;    // sizeof(char) * bits_per_char


	class TagMemory //implemented in tag.cpp
	{
	public:
		~TagMemory();
		TagMemory(Mesh * m, const TagMemory & other);
		TagMemory & operator =(TagMemory const & other);
	private:
		TagMemory();
		INMOST_DATA_ENUM_TYPE  pos[__NET]; 
		DataType               dtype;                 
		std::string            tagname;            
		INMOST_MPI_Type        bulk_data_type;   
		INMOST_DATA_ENUM_TYPE  size;       
		bool                   sparse[__NET];
		INMOST_DATA_ENUM_TYPE  record_size;
		TagManager *           tag_manager;
		Mesh *                 m_link;
		friend class Tag;
		friend class Storage; //for debug
	};

	class Tag //implemented in tag.cpp
	{
	private:
		TagMemory * mem;
		Tag(Mesh * m, std::string name, DataType _dtype, INMOST_DATA_ENUM_TYPE size); //Create me through TagManager::CreateTag
		__INLINE INMOST_DATA_ENUM_TYPE  GetRecordSize() const {return mem->record_size;}
		__INLINE void                   SetSize(INMOST_DATA_ENUM_TYPE size) {mem->size = size;}
		__INLINE void                   SetPosition(INMOST_DATA_ENUM_TYPE pos, ElementType type) {mem->pos[ElementNum(type)] = pos;}
		__INLINE INMOST_DATA_ENUM_TYPE  GetPosition(ElementType type) const {assert(mem != NULL); return mem->pos[ElementNum(type)];}
		__INLINE void                   SetSparse(ElementType type) {mem->sparse[ElementNum(type)] = true;}
		__INLINE void                   SetTagManager(TagManager * tm) {mem->tag_manager = tm;}
		__INLINE TagManager *           GetTagManager() const {return mem->tag_manager;}
		void                            ReallocateData(ElementType t,INMOST_DATA_ENUM_TYPE new_size);
		__INLINE INMOST_DATA_ENUM_TYPE  GetPositionNum(INMOST_DATA_ENUM_TYPE typenum) const {assert(mem != NULL); return mem->pos[typenum];}
		__INLINE void                   SetBulkDataType(INMOST_MPI_Type type){assert(mem!=NULL); if( mem->dtype == DATA_BULK ) mem->bulk_data_type = type;}
	public:
		~Tag() {mem = NULL;}
		Tag() {mem = NULL;}
		Tag(const Tag & other) {mem = other.mem;}
		bool operator <(const Tag & other) const {return mem < other.mem;}
		bool operator >(const Tag & other) const {return mem > other.mem;}
		bool operator ==(const Tag & other) const {return mem == other.mem;}
		bool operator !=(const Tag & other) const {return mem != other.mem;}
		Tag & operator =(Tag const & other) {mem = other.mem; return *this;	}
		__INLINE DataType               GetDataType() const {assert(mem!=NULL); return mem->dtype;}
		__INLINE INMOST_MPI_Type        GetBulkDataType() const {assert(mem!=NULL); return mem->bulk_data_type;}
		__INLINE INMOST_DATA_ENUM_TYPE  GetBytesSize() const {assert(mem!=NULL); return DataTypeBytesSize(mem->dtype);}
		__INLINE INMOST_DATA_ENUM_TYPE  GetSize() const {assert(mem!=NULL); return mem->size;}
		__INLINE std::string            GetTagName() const {assert(mem!=NULL); return mem->tagname;}
		__INLINE bool                   isDefined(ElementType type) const {assert(mem!=NULL && OneType(type)); return GetPosition(type) != ENUMUNDEF;}
		__INLINE bool                   isSparse(ElementType type) const {assert(mem!=NULL && OneType(type)); return mem->sparse[ElementNum(type)];}
		__INLINE bool                   isValid() const {return mem != NULL;}
		__INLINE Mesh *                 GetMeshLink() const {assert(mem!=NULL); return mem->m_link;}		
		__INLINE bool                   isSparseNum(INMOST_DATA_ENUM_TYPE typenum) const {assert(mem!=NULL); return mem->sparse[typenum];}
		friend class TagManager;
		friend class Storage;
		friend class Mesh;
	};
	
	class TagManager //implemented in tag.cpp
	{
	protected:
		TagManager();
		TagManager(Mesh * m, const TagManager & other);
		TagManager & assign(Mesh * m, TagManager const & other);
	public:
		//typedef chunk_array<INMOST_DATA_ENUM_TYPE>   empty_data;
		//typedef chunk_array<Tag>                     tag_array_type;
		typedef std::vector<INMOST_DATA_ENUM_TYPE> 	   empty_data;
		typedef std::vector<Tag>                       tag_array_type;
		typedef tag_array_type::iterator               iteratorTag;
		typedef chunk_array<INMOST_DATA_BULK_TYPE,15>  dense_sub_type;
		typedef chunk_array< dense_sub_type,8>         dense_data_array_type;
		typedef struct{void * tag, * rec;}             sparse_sub_record;
		typedef array< sparse_sub_record >             sparse_sub_type;
		typedef chunk_array< sparse_sub_type, 15>      sparse_data_array_type;
		virtual ~TagManager();
		
		bool                            HaveTag            (std::string name) const;
		Tag                             GetTag             (std::string name) const;
		void                            ListTagNames       (std::vector<std::string> & list) const;
		Tag                             CreateTag          (Mesh * m, std::string name, DataType dtype, ElementType etype, ElementType sparse, INMOST_DATA_ENUM_TYPE size = ENUMUNDEF); 
		virtual Tag                     DeleteTag          (Tag tag, ElementType mask); 
		__INLINE INMOST_DATA_ENUM_TYPE  NumberOfTags       () const { return static_cast<INMOST_DATA_ENUM_TYPE>(tags.size()); }
		__INLINE iteratorTag            BeginTag           () {return tags.begin(); }
		__INLINE iteratorTag            EndTag             () {return tags.end(); }
		bool                            ElementDefined     (Tag const & tag, ElementType etype) const;
	protected:
		void                            ReallocateData     (ElementType etype);
	private:
		typedef tag_array_type::iterator            tag_iterator;
		typedef tag_array_type::const_iterator      tag_const_iterator;
		tag_array_type                              tags;
		empty_data                                  empty_dense_data;
		dense_data_array_type                       dense_data;
#if defined(NEW_SPARSE)
		sparse_data_array_type                      sparse_data[6];
		__INLINE sparse_sub_type const & GetSparseData      (int etypenum, int local_id) const {return sparse_data[etypenum][local_id];}
		__INLINE sparse_sub_type &       GetSparseData      (int etypenum, int local_id) {return sparse_data[etypenum][local_id];}
#endif
		__INLINE dense_sub_type const &  GetDenseData       (int pos) const {return dense_data[pos];}
		__INLINE dense_sub_type &        GetDenseData       (int pos) {return dense_data[pos];}
		friend class Storage;
		friend class Tag;
		friend class Mesh;
	};
	
	/// Base class for Mesh, Element, and ElementSet classes.
	/// This base class is used for the Mesh class, as well as Element classes, and ElementSet class.
	/// 
	/// Storage class is used for representing different data objects in memory.
	/// Each data object is associated with corresponding Tag.
	class Storage //implemented in storage.cpp
	{
	protected:
		Storage(const Storage & other);           // never use this
		Storage(Mesh * m, INMOST_DATA_ENUM_TYPE lid, const Storage & other); // use this instead
		Storage & assign(Mesh * m, INMOST_DATA_ENUM_TYPE lid, Storage const & other);
		Storage(Mesh *m, ElementType _etype);
#if !defined(NEW_SPARSE)
		typedef std::pair<Tag, void *>    sparse_sub_type;
		typedef array< sparse_sub_type >  sparse_data_array_type;
#else
		typedef TagManager::sparse_sub_type sparse_type;
		typedef TagManager::sparse_sub_record sparse_rec;
		__INLINE static sparse_rec mkrec(const Tag & t) {sparse_rec ret; ret.tag = t.mem; ret.rec = NULL; return ret;}
#endif
		static void                       CopyData(Tag t, void * adata, void * bdata);
		static void                       DestroyVariableData(Tag t, void * adata);
		void                              MoveData(INMOST_DATA_INTEGER_TYPE new_local_id);
		INMOST_DATA_ENUM_TYPE             etypenum;
		INMOST_DATA_ENUM_TYPE             local_id;
		Mesh *                            m_link;
#if !defined(NEW_MARKERS)
		MarkerType                        markers;
#endif
#if !defined(NEW_SPARSE)
		sparse_data_array_type            inner_data;
		__INLINE void * &                 GetSparseLink(const Tag & t) {for(sparse_data_array_type::iterator it = inner_data.begin(); it != inner_data.end(); ++it) if( it->first == t ) return it->second; inner_data.push_back(sparse_sub_type(t,(void *)NULL)); return inner_data.back().second;}
		__INLINE void *                   GetSparseLink(const Tag & t) const {for(sparse_data_array_type::const_iterator it = inner_data.begin(); it != inner_data.end(); ++it) if( it->first == t ) return it->second; return NULL;}
#else
		__INLINE sparse_type const &      SLink() const {return reinterpret_cast<TagManager *>(GetMeshLink())->GetSparseData(etypenum,LocalID());}
		__INLINE sparse_type &            SLink() {return reinterpret_cast<TagManager *>(GetMeshLink())->GetSparseData(etypenum,LocalID());}
		__INLINE void * &                 GetSparseLink(const Tag & t) {sparse_type & s = SLink(); for(int i = 0; i < s.size(); ++i) if( s[i].tag == t.mem ) return s[i].rec; s.push_back(mkrec(t)); return s.back().rec;}
		__INLINE void *                   GetSparseLink(const Tag & t) const {sparse_type const & s = SLink(); for(int i = 0; i < s.size(); ++i) if( s[i].tag == t.mem ) return s[i].rec; return NULL;}
#endif
		__INLINE void *                   GetDenseLink (const Tag & t) const {return &(t.GetTagManager()->GetDenseData(t.GetPositionNum(etypenum))[LocalID()]);}
		//__INLINE void *                   GetDenseLink (const Tag & t) const {return &(t.mem->tag_manager->dense_data[t.mem->pos[etypenum]][local_id*t.mem->record_size]);}
		__INLINE void *                   GetLink      (const Tag & t) {void * p; if( !t.isSparseNum(etypenum) ) p = GetDenseLink(t); else {void * & q = GetSparseLink(t); if( q == NULL ) q = calloc(1,t.GetRecordSize()); p = q;} return p;}
		__INLINE void *                   GetLink      (const Tag & t) const {void * p; if( !t.isSparseNum(etypenum) ) p = GetDenseLink(t); else p = GetSparseLink(t); return p;}
	public:
		/// Storage type for representing real values.
		typedef INMOST_DATA_REAL_TYPE     real; 
		/// Storage type for representing integer values.
		typedef INMOST_DATA_INTEGER_TYPE  integer;
		/// Storage type for representing one byte of abstact data.
		typedef INMOST_DATA_BULK_TYPE     bulk;
		/// Storage type for representing references to Element.
		typedef Element *                 reference;
		/// Storage type for representing arrays of real values.
		typedef shell<real>               real_array;
		/// Storage type for representing arrays of integer values.
		typedef shell<integer>            integer_array;
		/// Storage type for representing abstact data as a series of bytes.
		typedef shell<bulk>               bulk_array;
		/// Storage type for representing arrays of Element references.
		typedef shell<reference>          reference_array;
		virtual ~Storage();
		/// Retrieve real value associated with Tag.
		__INLINE real      &              Real            (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_REAL     ); void * p = GetLink(tag); if( tag.GetSize() != ENUMUNDEF ) return static_cast<Storage::real     *>(p)[0]; else return static_cast<inner_real_array     *>(p)->at_safe(0);}
		/// Retrieve integer value associated with Tag.
		__INLINE integer   &              Integer         (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_INTEGER  ); void * p = GetLink(tag); if( tag.GetSize() != ENUMUNDEF ) return static_cast<Storage::integer  *>(p)[0]; else return static_cast<inner_integer_array  *>(p)->at_safe(0);}
		/// Retrieve one byte of abstract data associated with Tag.
		__INLINE bulk      &              Bulk            (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_BULK     ); void * p = GetLink(tag); if( tag.GetSize() != ENUMUNDEF ) return static_cast<Storage::bulk     *>(p)[0]; else return static_cast<inner_bulk_array     *>(p)->at_safe(0);}
		/// Retrieve Element reference associated with Tag.
		__INLINE reference &              Reference       (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_REFERENCE); void * p = GetLink(tag); if( tag.GetSize() != ENUMUNDEF ) return static_cast<Storage::reference*>(p)[0]; else return static_cast<inner_reference_array*>(p)->at_safe(0);}
		/// Retrieve array of real values associated with Tag.
		__INLINE real_array               RealArray       (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_REAL     ); void * p = GetLink(tag); if( tag.GetSize() == ENUMUNDEF ) return Storage::real_array     (*static_cast<inner_real_array     *>(p)); else return Storage::real_array     (static_cast<Storage::real      *>(p),tag.GetSize());}
		/// Retrieve array of integer values associated with Tag.
		__INLINE integer_array            IntegerArray    (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_INTEGER  ); void * p = GetLink(tag); if( tag.GetSize() == ENUMUNDEF ) return Storage::integer_array  (*static_cast<inner_integer_array  *>(p)); else return Storage::integer_array  (static_cast<Storage::integer   *>(p),tag.GetSize());}
		/// Retrieve abstract data associated with Tag as a series of bytes.
		__INLINE bulk_array               BulkArray       (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_BULK     ); void * p = GetLink(tag); if( tag.GetSize() == ENUMUNDEF ) return Storage::bulk_array     (*static_cast<inner_bulk_array     *>(p)); else return Storage::bulk_array     (static_cast<Storage::bulk      *>(p),tag.GetSize());}
		/// Retrieve array of Element references associated with Tag.
		__INLINE reference_array          ReferenceArray  (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_REFERENCE); void * p = GetLink(tag); if( tag.GetSize() == ENUMUNDEF ) return Storage::reference_array(*static_cast<inner_reference_array*>(p)); else return Storage::reference_array(static_cast<Storage::reference *>(p),tag.GetSize());}
		
		//optimized data requests for dense data with fixed size
		__INLINE real_array               RealArrayDF     (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_REAL      && tag.GetSize() != ENUMUNDEF && !tag.isSparse(GetElementType())); return Storage::real_array     (static_cast<Storage::real     *>(GetDenseLink(tag)),tag.GetSize());}
		__INLINE integer_array            IntegerArrayDF  (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_INTEGER   && tag.GetSize() != ENUMUNDEF && !tag.isSparse(GetElementType())); return Storage::integer_array  (static_cast<Storage::integer  *>(GetDenseLink(tag)),tag.GetSize());}
		__INLINE bulk_array               BulkArrayDF     (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_BULK      && tag.GetSize() != ENUMUNDEF && !tag.isSparse(GetElementType())); return Storage::bulk_array     (static_cast<Storage::bulk     *>(GetDenseLink(tag)),tag.GetSize());}
		__INLINE reference_array          ReferenceArrayDF(const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_REFERENCE && tag.GetSize() != ENUMUNDEF && !tag.isSparse(GetElementType())); return Storage::reference_array(static_cast<Storage::reference*>(GetDenseLink(tag)),tag.GetSize());}
		__INLINE real      &              RealDF          (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_REAL      && tag.GetSize() != ENUMUNDEF && !tag.isSparse(GetElementType()));  return static_cast<Storage::real     *>(GetDenseLink(tag))[0];}
		__INLINE integer   &              IntegerDF       (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_INTEGER   && tag.GetSize() != ENUMUNDEF && !tag.isSparse(GetElementType()));  return static_cast<Storage::integer  *>(GetDenseLink(tag))[0];}
		__INLINE bulk      &              BulkDF          (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_BULK      && tag.GetSize() != ENUMUNDEF && !tag.isSparse(GetElementType()));  return static_cast<Storage::bulk     *>(GetDenseLink(tag))[0];}
		__INLINE reference &              ReferenceDF     (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_REFERENCE && tag.GetSize() != ENUMUNDEF && !tag.isSparse(GetElementType()));  return static_cast<Storage::reference*>(GetDenseLink(tag))[0];}
		
		//optimized data requests for dense data with variable size
		__INLINE real_array               RealArrayDV     (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_REAL      && tag.GetSize() == ENUMUNDEF && !tag.isSparse(GetElementType())); return Storage::real_array     (*static_cast<inner_real_array     *>(GetDenseLink(tag)));}
		__INLINE integer_array            IntegerArrayDV  (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_INTEGER   && tag.GetSize() == ENUMUNDEF && !tag.isSparse(GetElementType())); return Storage::integer_array  (*static_cast<inner_integer_array  *>(GetDenseLink(tag)));}
		__INLINE bulk_array               BulkArrayDV     (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_BULK      && tag.GetSize() == ENUMUNDEF && !tag.isSparse(GetElementType())); return Storage::bulk_array     (*static_cast<inner_bulk_array     *>(GetDenseLink(tag)));}
		__INLINE reference_array          ReferenceArrayDV(const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_REFERENCE && tag.GetSize() == ENUMUNDEF && !tag.isSparse(GetElementType())); return Storage::reference_array(*static_cast<inner_reference_array*>(GetDenseLink(tag)));}
		__INLINE real      &              RealDV          (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_REAL      && tag.GetSize() == ENUMUNDEF && !tag.isSparse(GetElementType()));  return static_cast<inner_real_array     *>(GetDenseLink(tag))->at_safe(0);}
		__INLINE integer   &              IntegerDV       (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_INTEGER   && tag.GetSize() == ENUMUNDEF && !tag.isSparse(GetElementType()));  return static_cast<inner_integer_array  *>(GetDenseLink(tag))->at_safe(0);}
		__INLINE bulk      &              BulkDV          (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_BULK      && tag.GetSize() == ENUMUNDEF && !tag.isSparse(GetElementType()));  return static_cast<inner_bulk_array     *>(GetDenseLink(tag))->at_safe(0);}
		__INLINE reference &              ReferenceDV     (const Tag & tag) {assert(tag.isValid() && GetMeshLink() == tag.GetMeshLink() && tag.GetDataType() == DATA_REFERENCE && tag.GetSize() == ENUMUNDEF && !tag.isSparse(GetElementType()));  return static_cast<inner_reference_array*>(GetDenseLink(tag))->at_safe(0);}
		
		/// Return the data length associated with Tag.
		/// For abstract data return the number of bytes, otherwise return the length of associated array. 
		/// @see Storage::SetDataSize
		INMOST_DATA_ENUM_TYPE            GetDataSize      (const Tag & tag) const; //For DATA_BULK return number of bytes, otherwise return the length of array
		/// Set the length of  data associated with Tag.
		/// @param tag Identifying Tag.
		/// @param new_size The number of bytes for abstract data, otherwise the length of the array.
		/// @see Storage::GetDataSize
		void                             SetDataSize      (const Tag & tag,INMOST_DATA_ENUM_TYPE new_size);
		/// Extract part of the data associated with Tag.
		/// Copy part of the associated array or data to the destination memory.
		/// @param tag Identifying Tag.
		/// @param shift Starting position of the copied data.
		/// For abstact data – number of bytes to skip, otherwise number of values to skip.
		/// @param size Number of elements to copy.
		/// For abstact data – number of bytes to copy, otherwise number of values to copy.
		/// @param data Destination position to copy data to.
		/// @see Storage::SetData
		void                                   GetData        (const Tag & tag,INMOST_DATA_ENUM_TYPE shift, INMOST_DATA_ENUM_TYPE size, void * data) const;
		void                                   SetData        (const Tag & tag,INMOST_DATA_ENUM_TYPE shift, INMOST_DATA_ENUM_TYPE size, void * data);
		void                                   DelData        (const Tag & tag);
		/// Check if any data is associated with Tag.
		bool                                   HaveData       (const Tag & tag) const {if(tag.isSparseNum(etypenum)) { if( GetSparseLink(tag) != NULL ) return true; return false; } else {if( tag.GetPositionNum(etypenum) != ENUMUNDEF ) return true; return false;}}
		/// Swap dense data with the element of the same type and given local_id 
		void                                   SwapDenseData  (INMOST_DATA_INTEGER_TYPE local_id);
		/// Swap sparse data with the element of the same type and given local_id 
		void                                   SwapSparseData (INMOST_DATA_INTEGER_TYPE local_id);
		__INLINE ElementType                   GetElementType () const {return 1 << etypenum;}
		__INLINE INMOST_DATA_ENUM_TYPE         GetElementNum  () const {return etypenum;}
#if defined(NEW_MARKERS)
		void                                   SetMarker      (MarkerType n);
		bool                                   GetMarker      (MarkerType n) const;
		void                                   RemMarker      (MarkerType n) ;
		void                                   ClearMarkerSpace();
		void                                   GetMarkerSpace(Storage::bulk copy[MarkerFields]) const;
		void                                   SetMarkerSpace(Storage::bulk source[MarkerFields]);
		static INMOST_DATA_ENUM_TYPE           MaxMarker() {return MarkerShift * MarkerFields;}
#else
		__INLINE void                          SetMarker       (MarkerType n)  {markers |= n;}
		__INLINE bool                          GetMarker       (MarkerType n) const  {return (markers & n) != 0;}
		__INLINE void                          RemMarker       (MarkerType n) {markers &= ~n;}
		__INLINE void                          ClearMarkerSpace() {markers = 0;}
		__INLINE MarkerType                    GetMarkerSpace  () const {return markers;}
		__INLINE void                          SetMarkerSpace  (MarkerType _markers) {markers = _markers;}
		__INLINE static INMOST_DATA_ENUM_TYPE  MaxMarker() {return sizeof(MarkerType)*8;}
#endif
		__INLINE INMOST_DATA_ENUM_TYPE         LocalID         () const {return local_id;}
		__INLINE Mesh *                        GetMeshLink     () const {return m_link;}
		friend void                            SwapElement     (void * pa, void * pb, void * udata);
		friend class Mesh;
	};

	
	class ElementSet : public Storage //implemented in eset.cpp
	{
	private:
		typedef std::set<Element *,bool(*)(Element *,Element *)> element_set_type;
		bool ordered;
	public:
		class iterator : public element_set_type::iterator
		{
		public:
			iterator() {}
			iterator(const element_set_type::iterator & other ) : element_set_type::iterator(other) { }
			iterator(const iterator & other ) : element_set_type::iterator(other) { }
			iterator & operator =(iterator const & other) { element_set_type::iterator::operator=(static_cast<element_set_type::iterator const &>(other)); return *this; }
			Element & operator *() { return *element_set_type::iterator::operator *(); }
			Element * operator->() { return element_set_type::iterator::operator *(); }
		};
		class reverse_iterator : public element_set_type::reverse_iterator
		{
		public:
			reverse_iterator() {};
			reverse_iterator(const element_set_type::reverse_iterator & other ) : element_set_type::reverse_iterator(other) { }
			reverse_iterator(const reverse_iterator & other ) : element_set_type::reverse_iterator(other) { }
			reverse_iterator & operator =(reverse_iterator const & other) { element_set_type::reverse_iterator::operator=(static_cast<element_set_type::reverse_iterator const &>(other)); return *this; }
			Element & operator *() { return *element_set_type::reverse_iterator::operator *(); }
			Element * operator->() { return element_set_type::reverse_iterator::operator *(); }
		};
		ElementSet(Mesh *m,  bool ordered);
		ElementSet(bool ordered = false);
		ElementSet(const ElementSet & other);
		ElementSet(Mesh *m, INMOST_DATA_ENUM_TYPE lid,const ElementSet & other);
		ElementSet & operator =(ElementSet const & other);
		~ElementSet();
		std::pair< ElementSet::iterator, bool > Insert(const Element * e);
		template<class InputIterator>
		void Insert(InputIterator first, InputIterator last) {isInputForwardIterators<Element *, InputIterator>(); eset.insert(first,last);}
		void Insert(const ElementSet & e);
		void Insert(std::vector<Element *> e);
		void Insert(array<Element *> e);
		bool Erase(Element * e);
		void Erase(iterator e);
		void Intersection(ElementSet other);
		void Union(ElementSet other);
		void Difference(ElementSet other);
		iterator find(Element * e);
		iterator begin();
		iterator end();
		reverse_iterator rbegin();
		reverse_iterator rend();
		INMOST_DATA_ENUM_TYPE size() const;
		void clear();
		bool empty() const;
		void SetElementsMarker(MarkerType marker);
		void RemElementsMarker(MarkerType marker);
		bool isOrdered(){return ordered;}
	private:
		element_set_type eset;
	};
	
	
	template <typename AdjacentType>
	class adjacent
	{
	public:
		//typedef typename std::vector<Element *> container_type;
		//typedef container_type::size_t enumerator;
		typedef dynarray<Element *,64> container_type;
		typedef container_type::enumerator enumerator;
	private:
		container_type container;
		adjacent(const container_type & other) : container(other) {}
		container_type & get_container() {return container;};
	public:
		adjacent() {}
		adjacent(enumerator n) : container(n) {}
		template<class InputIterator>
		adjacent(InputIterator first, InputIterator last) :container(first,last) {isInputForwardIterators<Element *, InputIterator>();}
		adjacent(const adjacent & other ) {container = other.container;}
		~adjacent() {container.clear();}
		adjacent & operator=(adjacent const & other) {container = other.container; return *this;}
		class iterator : public container_type::iterator
		{
		public:
			iterator(const container_type::iterator & other ) : container_type::iterator(other) { }
			iterator(const iterator & other ) : container_type::iterator(other) { }
			iterator & operator =(iterator const & other) { container_type::iterator::operator=(static_cast<container_type::iterator const &>(other)); return *this; }
			AdjacentType & operator *() { return static_cast<AdjacentType &>(*(container_type::iterator::operator *())); }
			AdjacentType * operator->() { return static_cast<AdjacentType *>(container_type::iterator::operator *()); }
		};
		class reverse_iterator : public container_type::reverse_iterator
		{
		public:
			reverse_iterator(const container_type::reverse_iterator & other ) : container_type::reverse_iterator(other) { }
			reverse_iterator(const reverse_iterator & other ) : container_type::reverse_iterator(other) { }
			reverse_iterator & operator =(reverse_iterator const & other) { container_type::reverse_iterator::operator=(static_cast<container_type::reverse_iterator const &>(other)); return *this; }
			AdjacentType & operator *() { return static_cast<AdjacentType &>(*container_type::reverse_iterator::operator *()); }
			AdjacentType * operator->() { return static_cast<AdjacentType *>(container_type::reverse_iterator::operator *()); }
		};
		template<class InputIterator>
		void insert(iterator pos,InputIterator pbeg, InputIterator pend) {container.insert(pos,pbeg,pend);}
		iterator erase(iterator pos) {return container.erase(pos);}
		__INLINE iterator begin() { return iterator(container.begin()); }
		__INLINE iterator end() { return iterator(container.end()); }
		__INLINE reverse_iterator rbegin() { return reverse_iterator(container.rbegin()); }
		__INLINE reverse_iterator rend() { return reverse_iterator(container.rend()); }
		__INLINE AdjacentType & operator [] (enumerator n) {return static_cast<AdjacentType &>(*(container[n]));}
		__INLINE AdjacentType & front() { return static_cast<AdjacentType &>(*(container.front())); }
		__INLINE AdjacentType & back() { return static_cast<AdjacentType &>(*(container.back())); }
		__INLINE AdjacentType & at(enumerator n) { return static_cast<AdjacentType &>(*(container.at(n))); }
		void swap(adjacent<AdjacentType> & other) {container.swap(other.container);}
		__INLINE void push_back(Element & x) {container.push_back(&x);}
		__INLINE void push_back(Element * x) {container.push_back(x);}
		void resize(enumerator n) {container.resize(n);}
		__INLINE bool empty() {return container.empty();}
		void clear() {container.clear();}
		void reserve(enumerator n) {container.reserve(n);}
		__INLINE enumerator size() const { return container.size(); }
		__INLINE AdjacentType ** data() {return reinterpret_cast<AdjacentType **>(container.data());}
		void unite(const adjacent<AdjacentType>  & other);
		void substract(const adjacent<AdjacentType>  & other);
		void intersect(const adjacent<AdjacentType>  & other);
	};
	
			
	class Element : public Storage //implemented in element.cpp
	{
	public:
		typedef INMOST_DATA_BULK_TYPE GeometricType;
		static const GeometricType Unset        = 0;
		static const GeometricType Vertex       = 1;
		static const GeometricType Line         = 2;
		static const GeometricType MultiLine    = 3;
		static const GeometricType Tri          = 4;
		static const GeometricType Quad         = 5;
		static const GeometricType Polygon      = 6;
		static const GeometricType MultiPolygon = 7;
		static const GeometricType Tet          = 8;
		static const GeometricType Hex          = 9;
		static const GeometricType Prism        = 10;
		static const GeometricType Pyramid      = 11;
		static const GeometricType Polyhedron   = 12;
		static const GeometricType Set          = 100;
		//enum GeometricType {Unset,Vertex,Line,MultiLine,Tri,Quad,Polygon,MultiPolygon,Tet,Hex,Prism,Pyramid,Polyhedron,Set};
		static const char *       GeometricTypeName(GeometricType t);
		static unsigned int       GetGeometricDimension(GeometricType m_type);
		typedef INMOST_DATA_BULK_TYPE Status;
		static const Status Owned  = 1;
		static const Status Shared = 2;
		static const Status Ghost  = 4;
		static const Status Any    = 0;
		static const char * StatusName(Status s);
	public:
		typedef inner_reference_array            adj_type;
		typedef adj_type::iterator               adj_iterator;
		typedef adj_type::const_iterator         const_adj_iterator;
		typedef adj_type::reverse_iterator       adj_reverse_iterator;
		typedef adj_type::const_reverse_iterator const_adj_reverse_iterator;
	protected:
#if defined(NEW_CONNECTIONS)
		adj_type &              HighConn();
		adj_type &              LowConn ();
		adj_type const &        HighConn() const;
		adj_type const &        LowConn () const;
#else
		adj_type            high_conn_;
		adj_type            low_conn_;
		__INLINE adj_type &              HighConn() {return high_conn_;}
		__INLINE adj_type &              LowConn () {return low_conn_;}
		__INLINE adj_type const &        HighConn() const {return high_conn_;}
		__INLINE adj_type const &        LowConn () const {return low_conn_;}
#endif
		void                    SetGeometricType(GeometricType t);
		friend class Mesh;
		friend class Node;
		friend class Edge;
		friend class Face;
		friend class Cell;
		friend int              CompareElementsUnique(Element * a,Element * b);
		friend int              CompareElementsCentroid(Element * a,Element * b);
	protected:
		Element(Mesh * m, ElementType _etype);
		Element(const Element & other);            // don't use this
		Element(Mesh * m, INMOST_DATA_ENUM_TYPE lid,const Element & other);  // use this instead
		Element & operator =(Element const & other);
	public:
		virtual ~Element();
		INMOST_DATA_ENUM_TYPE   nbAdjElements(ElementType _etype) const;
		adjacent<Element>       getAdjElements(ElementType _etype) const;  //unordered
		INMOST_DATA_ENUM_TYPE   nbAdjElements(ElementType _etype, MarkerType mask, bool invert_mask = false) const;
		adjacent<Element>       getAdjElements(ElementType _etype, MarkerType mask, bool invert_mask = false) const;  //unordered
		adjacent<Element>       BridgeAdjacencies(ElementType Bridge, ElementType Dest, MarkerType mask = 0, bool invert_mask = false);
		adjacent<Node>          BridgeAdjacencies2Node(ElementType Bridge, MarkerType mask = 0, bool invert_mask = false);
		adjacent<Edge>          BridgeAdjacencies2Edge(ElementType Bridge, MarkerType mask = 0, bool invert_mask = false);
		adjacent<Face>          BridgeAdjacencies2Face(ElementType Bridge, MarkerType mask = 0, bool invert_mask = false);
		adjacent<Cell>          BridgeAdjacencies2Cell(ElementType Bridge, MarkerType mask = 0, bool invert_mask = false);
		virtual adjacent<Node>  getNodes(); //unordered
		virtual adjacent<Edge>  getEdges(); //unordered
		virtual adjacent<Face>  getFaces(); //unordered
		virtual adjacent<Cell>  getCells(); //unordered
		virtual adjacent<Node>  getNodes(MarkerType mask,bool invert_mask = false); //unordered
		virtual adjacent<Edge>  getEdges(MarkerType mask,bool invert_mask = false); //unordered
		virtual adjacent<Face>  getFaces(MarkerType mask,bool invert_mask = false); //unordered
		virtual adjacent<Cell>  getCells(MarkerType mask,bool invert_mask = false); //unordered
		Node *                  getAsNode(); //does dynamic conversation for you, if not a node returns NULL
		Edge *                  getAsEdge(); //does dynamic conversation for you, if not an edge returns NULL
		Face *                  getAsFace(); //does dynamic conversation for you, if not a face returns NULL
		Cell *                  getAsCell(); //does dynamic conversation for you, if not a cell returns NULL
		GeometricType           GetGeometricType() const;
		unsigned int            GetElementDimension() const {return GetGeometricDimension(GetGeometricType());}
		Status                  GetStatus();
		void                    SetStatus(Status status);
		Storage::integer &      GlobalID();
		bool                    CheckElementConnectivity();
		static bool             CheckConnectivity(Mesh * m);
		//implemented in geometry.cpp
		void                    CastRay(Storage::real * pos, Storage::real * dir, dynarray< std::pair<Element *, Storage::real> , 16 > & hits);
		void                    ComputeGeometricType();
		void                    Centroid(Storage::real * cnt);
		void                    Barycenter(Storage::real * cnt);
		Storage::real           Mean(Storage::real (*func)(Storage::real* x,Storage::real t),Storage::real time);
		bool                    Boundary();
		bool                    Planarity(); // check that all nodes lay on one plane
		//implemented in modify.cpp
		bool                    Hide(); // if true then element was hidden, works only inside BeginModification and EndModification, on EndModification all Hidden elements are deleted
		bool                    Show(); // if true then element was recovered
		bool                    Delete(); // if true then element was deleted, otherwise it was hidden
		bool                    Hidden();
		bool                    New();
		void                    Disconnect(bool delete_upper_adjacent); //disconnect all elements, delete upper dependent
		/// Disconnects nodes from this edge, edges from this face, faces from this cell, cannot disconnect cells from this node;
		/// Disconnects edges from this node, faces from this edge, cells from this face, cannot disconnect nodes from this cell;
		/// Updates geometric data and cell nodes automatically.
		void                    Disconnect(Element ** adjacent, INMOST_DATA_ENUM_TYPE num);
		/// Connects lower adjacencies to current element, 
		/// geometric data and cell nodes are updated automatically.
		/// TODO:
		///		1. asserts in this function should be replaced by Topography checks;
		///		2. this function should be used for creation of elements instead of current implementation.
		///     3. should correctly account for order of edges (may be implemented through CheckEdgeOrder, FixEdgeOrder)
		void                    Connect(Element ** adjacent, INMOST_DATA_ENUM_TYPE num); 
		/// Update geometric data for element, calls RecomputeGeometricData from Mesh.
		void                    UpdateGeometricData(); 
	};
	
	class Node : public Element //implemented in node.cpp
	{
	private:
		Node(Mesh * m);
		Node(const Node & other);
		Node(Mesh * m, INMOST_DATA_ENUM_TYPE lid, const Node & other);
		Node & operator =(Node const & other);
	public:
		~Node();
		
		adjacent<Edge>              getEdges(); //unordered
		adjacent<Face>              getFaces(); //unordered
		adjacent<Cell>              getCells(); //unordered

		adjacent<Edge>              getEdges(MarkerType mask,bool invert_mask = false); //unordered
		adjacent<Face>              getFaces(MarkerType mask,bool invert_mask = false); //unordered
		adjacent<Cell>              getCells(MarkerType mask,bool invert_mask = false); //unordered

		Storage::real_array         Coords(); 
		friend class Mesh;
	};
	
	class Edge : public Element //implemented in edge.cpp
	{
	private:
		Edge(Mesh * m);
		Edge(const Edge & other);
		Edge(Mesh * m, INMOST_DATA_ENUM_TYPE lid, const Edge & other);
		Edge & operator =(Edge const & other);
	public:
		~Edge();
		
		adjacent<Node>              getNodes(); //ordered
		adjacent<Face>              getFaces(); //unordered
		adjacent<Cell>              getCells(); //unordered

		adjacent<Node>              getNodes(MarkerType mask,bool invert_mask = false); //ordered
		adjacent<Face>              getFaces(MarkerType mask,bool invert_mask = false); //unordered
		adjacent<Cell>              getCells(MarkerType mask,bool invert_mask = false); //unordered

		Node *                      getBeg() const;
		Node *                      getEnd() const;
		//implemented in modify.cpp
		static Edge *               UniteEdges    (Edge ** edges, INMOST_DATA_ENUM_TYPE nedges, MarkerType del_protect);
		static bool                 TestUniteEdges(Edge ** edges, INMOST_DATA_ENUM_TYPE nedges, MarkerType del_protect);
		static dynarray<Edge *,32>  SplitEdge     (Edge * e, Node ** nodes, INMOST_DATA_ENUM_TYPE nnodes, MarkerType del_protect); //provide ordered array of nodes, that lay between former nodes of the edge
		static bool                 TestSplitEdge (Edge * e, Node ** nodes, INMOST_DATA_ENUM_TYPE nnodes, MarkerType del_protect);
		//implemented in geometry.cpp
		Storage::real               Length();
		friend class Mesh;
	};
	
	class Face : public Element //implemented in face.cpp
	{
	private:
		Face(Mesh * m);
		Face(const Face & other);
		Face(Mesh * m, INMOST_DATA_ENUM_TYPE lid, const Face & other);
		Face & operator =(Face const & other);
	public:
		~Face();
		
		adjacent<Node>              getNodes(); //ordered
		adjacent<Edge>              getEdges(); //ordered
		adjacent<Cell>              getCells(); //unordered

		adjacent<Node>              getNodes(MarkerType mask,bool invert_mask = false); //ordered
		adjacent<Edge>              getEdges(MarkerType mask,bool invert_mask = false); //ordered
		adjacent<Cell>              getCells(MarkerType mask,bool invert_mask = false); //unordered

		//this is for 2d case when the face is represented by segment
		Node *                      getBeg() const;
		Node *                      getEnd() const;

		Cell *                      BackCell() const;
		Cell *                      FrontCell() const;
		bool                        FaceOrientedOutside(Cell * c) const;
		void                        ReorderEdges();
		bool                        CheckEdgeOrder(); //not implemented// returns true if edges of face form an ordered closed loop
		bool                        FixEdgeOrder(); //not implemented// returns true if edges were successfully reordered to form a closed loop
		//implemented in modify.cpp
		static Face *               UniteFaces    (Face ** faces, INMOST_DATA_ENUM_TYPE nfaces, MarkerType del_protect);
		static bool                 TestUniteFaces(Face ** faces, INMOST_DATA_ENUM_TYPE nfaces,  MarkerType del_protect);
		static dynarray<Face *,32>  SplitFace     (Face * face, Edge ** edges, INMOST_DATA_ENUM_TYPE nedges, MarkerType del_protect); //provide all edges that lay inside face
		static bool                 TestSplitFace (Face * face, Edge ** edges, INMOST_DATA_ENUM_TYPE nedges, MarkerType del_protect);	
		void                        SwapCells(); //swap back cell and front cell
		//implemented in geometry.cpp
		Storage::real               Area();
		void                        Normal(Storage::real * nrm);
		void                        UnitNormal(Storage::real * nrm);
		void                        OrientedNormal(Cell * c, Storage::real * nrm);
		void                        OrientedUnitNormal(Cell * c, Storage::real * nrm);
		bool                        FixNormalOrientation();  //returns true if orientation was corrected, otherwise returns false
		bool                        CheckNormalOrientation(); //returns true if orientation is correct, otherwise returns false
		bool                        Closure(); // test integrity of polygon
		friend class Mesh;
	};
	
	class Cell : public Element //implemented in cell.cpp
	{
	private:
		Cell(Mesh * m);
		Cell(const Cell & other);
		Cell(Mesh * m, INMOST_DATA_ENUM_TYPE lid, const Cell & other);
		Cell & operator =(Cell const & other);
	public:
		~Cell();
		
		adjacent<Node>              getNodes(); //ordered (for known geometric types only)
		adjacent<Edge>              getEdges(); //unordered
		adjacent<Face>              getFaces(); //ordered (keeps order it was created in)

		adjacent<Node>              getNodes(MarkerType mask,bool invert_mask = false); //ordered (for known geometric types only)
		adjacent<Edge>              getEdges(MarkerType mask,bool invert_mask = false); //unordered
		adjacent<Face>              getFaces(MarkerType mask,bool invert_mask = false); //ordered (keeps order it was created in)
		
		
		bool                        CheckEdgeOrder(); //not implemented//2D only, returns true if edges of face form an ordered closed loop
		bool                        FixEdgeOrder(); //not implemented//2D only, returns true if edges were successfully reordered to form a closed loop
		//implemented in modify.cpp
		static Cell *               UniteCells    (Cell ** cells, INMOST_DATA_ENUM_TYPE ncells, MarkerType del_protect);
		static bool                 TestUniteCells(Cell ** cells, INMOST_DATA_ENUM_TYPE ncells, MarkerType del_protect);
		static dynarray<Cell *,32>  SplitCell     (Cell * cell, Face ** faces, INMOST_DATA_ENUM_TYPE nfaces, MarkerType del_protect); //provide all faces, that lay inside cell
		static bool                 TestSplitCell (Cell * cell, Face ** faces, INMOST_DATA_ENUM_TYPE nfaces, MarkerType del_protect);
		//implemented in geometry.cpp
		Cell *                      Neighbour(Face * f);
		adjacent<Cell>              NeighbouringCells(); // get all cells that share any face with current
		bool                        Inside(Storage::real * point); //is point inside cell, check for 2d case
		Storage::real               Volume();
		bool                        Closure(); // test integrity of cell
		friend class Mesh;
	};
	
	class Mesh : public TagManager, public Storage //implemented in mesh.cpp
	{
	public:
		enum MeshState {Serial, Parallel};
		typedef chunk_array<Node *,15>                 nodes_container;
		typedef chunk_array<Edge *,15>                 edges_container;
		typedef chunk_array<Face *,15>                 faces_container;
		typedef chunk_array<Cell *,15>                 cells_container;
		typedef chunk_array<ElementSet *,7>            sets_container;
		typedef std::vector<INMOST_DATA_ENUM_TYPE>     empty_container;

		//typedef std::vector<Node *>                 nodes_container;
		//typedef std::vector<Edge *>                 edges_container;
		//typedef std::vector<Face *>                 faces_container;
		//typedef std::vector<Cell *>                 cells_container;
		//typedef std::vector<ElementSet *>           sets_container;
		//typedef std::vector<INMOST_DATA_ENUM_TYPE>  empty_container;
	private:
		Storage::real                               epsilon;
		cells_container                             cells;
		empty_container                             empty_cells;
		faces_container                             faces;
		empty_container                             empty_faces;
		edges_container                             edges;
		empty_container                             empty_edges;
		nodes_container                             nodes;
		empty_container                             empty_nodes;
		sets_container                              sets;
		empty_container                             empty_sets;
		Tag                                         tag_global_id;
		Tag                                         tag_coords;
		Tag                                         tag_low_conn;
		Tag                                         tag_high_conn;
		Tag                                         tag_markers;
		Tag                                         tag_geom_type;
		MeshState                                   m_state;
		unsigned int                                dim;
		Element *                                   last_created_element;
	public:
		Mesh();
		Mesh(const Mesh & other);
		~Mesh();
		Mesh & operator =(Mesh const & other);
		MarkerType                   CreateMarker();
		void                         ReleaseMarker(MarkerType n);
		__INLINE void                SetEpsilon(Storage::real e) {epsilon = e;}
		__INLINE Storage::real       GetEpsilon() const {return epsilon;}
		void                         SetDimensions(unsigned int dim);
		__INLINE unsigned int        GetDimensions() const {return dim;}
		__INLINE MeshState           GetMeshState() const {return m_state;}
		__INLINE const Tag &         GlobalIDTag() const {return tag_global_id;}
		__INLINE const Tag &         CoordsTag() const {return tag_coords;}
		__INLINE const Tag &         LowConnTag() const {return tag_low_conn;}
		__INLINE const Tag &         HighConnTag() const {return tag_high_conn;}
		__INLINE const Tag &         MarkersTag() const {return tag_markers;}
		__INLINE const Tag &         GeomTypeTag() const {return tag_geom_type;}
		Tag                          CreateTag(std::string name, DataType dtype, ElementType etype,ElementType sparse, INMOST_DATA_ENUM_TYPE size = ENUMUNDEF);
		Tag                          DeleteTag(Tag tag, ElementType mask = NODE | EDGE | FACE | CELL | ESET | MESH);
		Node *                       CreateNode(Storage::real * coords);
		std::pair<Edge *,bool>       CreateEdge(Node ** e_nodes, INMOST_DATA_ENUM_TYPE num_e_nodes);
		std::pair<Face *,bool>       CreateFace(Edge ** f_edges, INMOST_DATA_ENUM_TYPE num_f_edges);
		std::pair<Face *,bool>       CreateFace(Node ** f_nodes, INMOST_DATA_ENUM_TYPE num_f_nodes);
		std::pair<Cell *,bool>       CreateCell(Face ** c_faces, INMOST_DATA_ENUM_TYPE num_c_faces, 
												Node ** suggest_nodes_order = NULL, INMOST_DATA_ENUM_TYPE numsuggest_nodes_order = 0);
		std::pair<Cell *,bool>       CreateCell(Node ** c_f_nodes, const INMOST_DATA_ENUM_TYPE * c_f_numnodes, INMOST_DATA_ENUM_TYPE num_c_faces, 
												Node ** suggest_nodes_order = NULL, INMOST_DATA_ENUM_TYPE numsuggest_nodes_order = 0);
		std::pair<Cell *,bool>       CreateCell(Node ** c_nodes, const INMOST_DATA_ENUM_TYPE * c_f_nodeinds, const INMOST_DATA_ENUM_TYPE * c_f_numnodes, INMOST_DATA_ENUM_TYPE num_c_faces, 
												Node ** suggest_nodes_order = NULL, INMOST_DATA_ENUM_TYPE numsuggest_nodes_order = 0);
		Element *                    ElementByLocalID(ElementType etype, INMOST_DATA_INTEGER_TYPE lid)
		{
			assert(OneType(etype));
			switch(etype)
			{
			case NODE: return nodes[lid];
			case EDGE: return edges[lid];
			case FACE: return faces[lid];
			case CELL: return cells[lid];
			}
			return NULL;
		}
		Node *                       NodeByLocalID(INMOST_DATA_INTEGER_TYPE lid) { return nodes[lid]; }
		Edge *                       EdgeByLocalID(INMOST_DATA_INTEGER_TYPE lid) { return edges[lid]; }
		Face *                       FaceByLocalID(INMOST_DATA_INTEGER_TYPE lid) {return  faces[lid];}
		Cell *                       CellByLocalID(INMOST_DATA_INTEGER_TYPE lid) { return cells[lid]; }
		ElementSet *                 EsetByLocalID(INMOST_DATA_INTEGER_TYPE lid) { return sets[lid]; }

		INMOST_DATA_INTEGER_TYPE     MaxLocalIDNODE() const {return static_cast<INMOST_DATA_INTEGER_TYPE>(nodes.size());}
		INMOST_DATA_INTEGER_TYPE     MaxLocalIDEDGE() const {return static_cast<INMOST_DATA_INTEGER_TYPE>(edges.size());}
		INMOST_DATA_INTEGER_TYPE     MaxLocalIDFACE() const {return static_cast<INMOST_DATA_INTEGER_TYPE>(faces.size());}
		INMOST_DATA_INTEGER_TYPE     MaxLocalIDCELL() const {return static_cast<INMOST_DATA_INTEGER_TYPE>(cells.size());}
		INMOST_DATA_INTEGER_TYPE     MaxLocalIDESET() const {return static_cast<INMOST_DATA_INTEGER_TYPE>(sets.size());}
		INMOST_DATA_INTEGER_TYPE     MaxLocalID(ElementType etype) const
		{
			assert(OneType(etype));
			switch(etype)
			{
			case NODE: return static_cast<INMOST_DATA_INTEGER_TYPE>(nodes.size());
			case EDGE: return static_cast<INMOST_DATA_INTEGER_TYPE>(edges.size());
			case FACE: return static_cast<INMOST_DATA_INTEGER_TYPE>(faces.size());
			case CELL: return static_cast<INMOST_DATA_INTEGER_TYPE>(cells.size());
			}
			return 0;
		}
		ElementSet * CreateSet();
		ElementSet * CreateOrderedSet();
		Element * FindSharedAdjacency(Element * const * arr, unsigned num) ;
		void ReorderEmpty(ElementType reordertypes);
		void ReorderApply(Tag index, ElementType mask);
		bool isOriginal(Element * e); //don't know now what this function was for, should detect the copy-constructed element, but copy-construction is prohibited
		INMOST_DATA_ENUM_TYPE GetArrayCapacity(ElementType etype); //This function is needed by TagManager, may be made private in future
		void RestoreCellNodes(Cell * c, dynarray<Node *,64> & ret);
	private:
		void MoveStorage(Storage * e, int new_local_id);
		void UntieElement(Storage * e);
		void TieElement(Storage * e);
		//implemented in mesh_parallel.cpp
	public:
		enum Action {AGhost, AMigrate};
		enum Prepare {UnknownSize, UnknownSource};
		typedef std::pair<int, std::vector<INMOST_DATA_BULK_TYPE> > proc_buffer_type;
		typedef std::vector< proc_buffer_type > exch_buffer_type;
		class exchange_data
		{
		public:
			std::vector<INMOST_MPI_Request> send_reqs, recv_reqs;
			exch_buffer_type send_buffers, recv_buffers;
		};
	private:
		ElementType have_global_id;
		INMOST_DATA_BIG_ENUM_TYPE parallel_mesh_unique_id;
		INMOST_MPI_Comm comm;
		Tag tag_shared, tag_owner, tag_processors, tag_layers, tag_sendto, tag_bridge, tag_redistribute;
		void ComputeSharedProcs();
		std::map< int , std::vector<Element *> > ComputeSharedSkinSet(ElementType bridge);
		void PackTagData(Tag tag, std::vector< std::vector<Element *> > & elements, ElementType mask, std::vector<INMOST_DATA_BULK_TYPE> & buffer);
		void UnpackTagData(Tag tag, std::vector< std::vector<Element *> > & elements, ElementType mask, std::vector<INMOST_DATA_BULK_TYPE> & buffer, int & position, void (*Operation)(Tag tag, Element * element,INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size));
		void PackElementsData(std::vector<Element *> & input, std::vector<INMOST_DATA_BULK_TYPE> & buffer, int destination, std::vector<std::string> tag_list);
		void UnpackElementsData(std::vector<Element *> & output, std::vector<INMOST_DATA_BULK_TYPE> & buffer, int source, std::vector<std::string> & tag_list);
		class elements_by_type
		{
		private:
			std::vector< std::vector<Element *> > container;
		public:
			elements_by_type() :container(4) {}
			elements_by_type(const elements_by_type & other){container = other.container;}
			~elements_by_type(){}
			std::vector<Element *> & operator [](int i){ return container[i]; }
			std::vector<std::vector<Element *> >::iterator begin(){ return container.begin(); }
			std::vector<std::vector<Element *> >::iterator end() { return container.end(); }
			std::vector<std::vector<Element *> > & get_container() {return container;}
		};
		typedef std::map<int, elements_by_type > parallel_storage;
#if defined(USE_PARALLEL_STORAGE)
		parallel_storage shared_elements, ghost_elements;
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		int num_exchanges;
		std::fstream out_time;
		int tab;
		int func_id;
#endif
#if defined(USE_MPI2)
		INMOST_MPI_Win window;
		unsigned * shared_space;
#endif
		int parallel_strategy, parallel_file_strategy;
		void PrepareReceiveInner(Prepare todo, exch_buffer_type & send_bufs, exch_buffer_type & recv_bufs);
		void ExchangeDataInnerBegin(std::vector<Tag> tag, parallel_storage & from, parallel_storage & to, ElementType mask, exchange_data & storage);
		void ExchangeDataInnerEnd(std::vector<Tag> tag, parallel_storage & from, parallel_storage & to, ElementType mask, void (*Operation)(Tag tag,Element * element,INMOST_DATA_BULK_TYPE * recv_data,INMOST_DATA_ENUM_TYPE recv_size), exchange_data & storage);
		void ExchangeBuffersInner(exch_buffer_type & send_bufs, exch_buffer_type & recv_bufs,std::vector<INMOST_MPI_Request> & send_reqs, std::vector<INMOST_MPI_Request> & recv_reqs);
		std::vector<int> FinishRequests(std::vector<INMOST_MPI_Request> & recv_reqs);
		void GatherParallelStorage(parallel_storage & ghost, parallel_storage & shared, ElementType mask);
		class Random // random generator to provide tag for communication
		{
		private: unsigned int n,a,c,m;
		public:
			Random(unsigned int seed = 50);
			Random(const Random & other);
			unsigned int Number();
		} randomizer;
	public:
#if defined(USE_PARALLEL_WRITE_TIME)	
		//this part is needed to test parallel performance
		void Enter();
		void Exit();
		int & GetFuncID() {return func_id;}
		std::fstream & GetStream();
		std::fstream & WriteTab(std::fstream & f);
#endif
		static void Initialize(int * argc, char *** argv);
		static void Finalize();
		void SetParallelStrategy(int strategy){if( strategy < 0 || strategy > 3 ) throw NotImplemented; parallel_strategy = strategy;}
		int GetParallelStrategy() {return parallel_strategy;}
		void SetParallelFileStrategy(int strategy){if( strategy < 0 || strategy > 1 ) throw NotImplemented; parallel_file_strategy = strategy;}
		int GetParallelFileStrategy() {return parallel_file_strategy;}
		int GetProcessorRank();
		int GetProcessorsNumber();
		INMOST_MPI_Comm GetCommunicator();
		void SetCommunicator(INMOST_MPI_Comm _comm);
		void ResolveShared();
		void RemoveGhost();
		void RemoveGhostElements(std::vector<Element *> ghost);
		void AssignGlobalID(ElementType mask);
		//exchange for tag of type DATA_REFERENCE is not implemented!
		void ExchangeData(Tag tag, ElementType mask);
		void ExchangeDataBegin(Tag tags, ElementType mask, exchange_data & storage);
		void ExchangeDataEnd(Tag tags, ElementType mask, exchange_data & storage);
		void ExchangeData(std::vector<Tag> tags, ElementType mask);
		void ExchangeDataBegin(std::vector<Tag> tags, ElementType mask, exchange_data & storage);
		void ExchangeDataEnd(std::vector<Tag> tags, ElementType mask, exchange_data & storage);
		void ReduceData(Tag tag, ElementType mask, void (*Operation)(Tag tag, Element * element, INMOST_DATA_BULK_TYPE * recv_data, INMOST_DATA_ENUM_TYPE recv_size) );
		void ReduceDataBegin(Tag tag, ElementType mask, exchange_data & storage);
		void ReduceDataEnd(Tag tag, ElementType mask, void (*Operation)(Tag tag, Element * element, INMOST_DATA_BULK_TYPE * recv_data, INMOST_DATA_ENUM_TYPE recv_size), exchange_data & storage );
		void ReduceData(std::vector<Tag> tags, ElementType mask, void (*Operation)(Tag tag, Element * element, INMOST_DATA_BULK_TYPE * recv_data, INMOST_DATA_ENUM_TYPE recv_size) );
		void ReduceDataBegin(std::vector<Tag> tags, ElementType mask, exchange_data & storage);
		void ReduceDataEnd(std::vector<Tag> tags, ElementType mask, void (*Operation)(Tag tag, Element * element, INMOST_DATA_BULK_TYPE * recv_data, INMOST_DATA_ENUM_TYPE recv_size), exchange_data & storage );
		void ExchangeMarked(enum Action action = AGhost);
		void ExchangeGhost(Storage::integer layers, ElementType bridge);
		void Redistribute();
		Storage::integer Enumerate(ElementType mask, Tag num_tag, Storage::integer start = 0);
		Storage::integer Enumerate(std::vector<Element *> elements, Tag num_tag, Storage::integer start = 0);
		Storage::integer EnumerateSet(ElementSet * set, Tag num_tag, Storage::integer start = 0);
		Storage::integer TotalNumberOf(ElementType mask);
		Storage::real Integrate(Storage::real input);
		Storage::integer Integrate(Storage::integer input);
		Storage::integer ExclusiveSum(Storage::integer input); 
		Storage::real Integrate(Tag t,ElementType mask);
		Storage::real AggregateMax(Storage::real input);
		Storage::integer AggregateMax(Storage::integer input);
		void RecomputeParallelStorage(ElementType mask);
		__INLINE const Tag SendtoTag() const {return tag_sendto;}
		__INLINE const Tag SharedTag() const {return tag_shared;}
		__INLINE const Tag OwnerTag() const {return tag_owner;}
		__INLINE const Tag LayersTag() const {return tag_layers;}
		__INLINE const Tag ProcessorsTag() const {return tag_processors;}
		__INLINE Tag RedistributeTag() {return CreateTag("TEMPORARY_NEW_OWNER",DATA_INTEGER,CELL,NONE,1);}
		
		ElementType SynchronizeElementType(ElementType etype);
		void SynchronizeMarker(MarkerType marker, ElementType mask, SyncBitOp op);
		
		//for debug
		void                 BeginSequentialCode();
		void                   EndSequentialCode();
		//iterator.cpp::::::::::::::::::::::::::::::::::::::::::::::::::
	public:
		class base_iterator;
		class iteratorElement;
		class iteratorSet;
		class iteratorCell;
		class iteratorFace;
		class iteratorEdge;
		class iteratorNode;
		
		__INLINE INMOST_DATA_ENUM_TYPE NumberOfCells()   {return static_cast<INMOST_DATA_ENUM_TYPE>(cells.size() - empty_cells.size());}
		__INLINE INMOST_DATA_ENUM_TYPE NumberOfFaces()   { return static_cast<INMOST_DATA_ENUM_TYPE>(faces.size() - empty_faces.size()); }
		__INLINE INMOST_DATA_ENUM_TYPE NumberOfEdges()   { return static_cast<INMOST_DATA_ENUM_TYPE>(edges.size() - empty_edges.size()); }
		__INLINE INMOST_DATA_ENUM_TYPE NumberOfNodes()   { return static_cast<INMOST_DATA_ENUM_TYPE>(nodes.size() - empty_nodes.size()); }
		__INLINE INMOST_DATA_ENUM_TYPE NumberOfSets()    { return static_cast<INMOST_DATA_ENUM_TYPE>(sets.size() - empty_sets.size()); }
		__INLINE INMOST_DATA_ENUM_TYPE NumberOfElements(){ return NumberOfCells() + NumberOfFaces() + NumberOfEdges() + NumberOfNodes(); }
		__INLINE INMOST_DATA_ENUM_TYPE NumberOfAll()     { return NumberOfSets() + NumberOfElements(); }
		INMOST_DATA_ENUM_TYPE NumberOf(ElementType t);
		base_iterator Begin(ElementType Types);
		base_iterator End();
		iteratorElement BeginElement(ElementType Types);
		iteratorElement EndElement();
		iteratorSet  BeginSet();
		iteratorSet  EndSet();
		iteratorCell BeginCell();
		iteratorCell EndCell();
		iteratorFace BeginFace();
		iteratorFace EndFace();
		iteratorEdge BeginEdge();
		iteratorEdge EndEdge();
		iteratorNode BeginNode();
		iteratorNode EndNode();
		class base_iterator
		{
		protected:
			typedef int itenum;
			Mesh * m;
			itenum Number;
			ElementType CurrentType;
			ElementType types;
			base_iterator(ElementType Types, Mesh * mesh, bool last);
			base_iterator(Mesh * mesh) {m = mesh; CurrentType = NONE; types = NONE; Number = -1;}
		public:
			typedef Storage * pointer;
			typedef Storage & reference;
			typedef Storage value_type;
			typedef ptrdiff_t difference_type;
			typedef std::bidirectional_iterator_tag iterator_category;
			base_iterator() {Number = -1; CurrentType = NONE; m = NULL; types = NONE;}
			base_iterator(const base_iterator & other) {m = other.m; Number = other.Number;	types = other.types; CurrentType = other.CurrentType;}
			virtual ~base_iterator() {}
			base_iterator & operator ++();
			__INLINE base_iterator operator ++(int) {Mesh::base_iterator ret(*this); operator++(); return ret;}
			base_iterator & operator --();
			__INLINE base_iterator operator --(int) {Mesh::base_iterator ret(*this); operator--(); return ret;}
			__INLINE virtual reference operator *()
			{
				switch(CurrentType)
				{
					case NODE: return static_cast<Storage &>(*m->nodes[Number]); break;
					case EDGE: return static_cast<Storage &>(*m->edges[Number]); break;
					case FACE: return static_cast<Storage &>(*m->faces[Number]); break;
					case CELL: return static_cast<Storage &>(*m->cells[Number]); break;
					case ESET: return static_cast<Storage &>(*m->sets[Number]); break;
				}
				return *m;
			}
			__INLINE virtual pointer operator ->()
			{
				switch(CurrentType)
				{
					case NODE: return static_cast<Storage *>(m->nodes[Number]); break;
					case EDGE: return static_cast<Storage *>(m->edges[Number]); break;
					case FACE: return static_cast<Storage *>(m->faces[Number]); break;
					case CELL: return static_cast<Storage *>(m->cells[Number]); break;
					case ESET: return static_cast<Storage *>(m->sets[Number]); break;
				}
				return NULL;
			}
			__INLINE base_iterator & operator =(base_iterator const & other) {m = other.m; Number = other.Number; types = other.types; CurrentType = other.CurrentType; return *this;}
			__INLINE bool operator ==(const base_iterator & other) const { assert( m == other.m ); if( Number == other.Number && CurrentType == other.CurrentType ) return true; return false;}
			__INLINE bool operator !=(const base_iterator & other) const { assert( m == other.m ); if( Number != other.Number || CurrentType != other.CurrentType ) return true; return false;}
			__INLINE bool operator <(const base_iterator & other) const { assert( m == other.m ); if( (CurrentType < other.CurrentType) || (CurrentType == other.CurrentType && Number < other.Number) ) return true; return false;}
			__INLINE bool operator >(const base_iterator & other) const { assert( m == other.m ); if( (CurrentType > other.CurrentType) || (CurrentType == other.CurrentType && Number > other.Number) ) return true; return false;}
			__INLINE bool operator <=(const base_iterator & other) const {assert( m != other.m ); if( (CurrentType < other.CurrentType) || (CurrentType == other.CurrentType && Number <= other.Number) ) return true; return false;}
			__INLINE bool operator >=(const base_iterator & other) const {assert( m == other.m ); if( (CurrentType > other.CurrentType) || (CurrentType == other.CurrentType && Number >= other.Number) ) return true; return false;}
			void Print();
			friend base_iterator Mesh::Begin(ElementType Types);
			friend base_iterator Mesh::End();
		};
		class iteratorElement : public base_iterator
		{
		private:
			iteratorElement(const base_iterator & other) :base_iterator(other) {}
		public:
			__INLINE Element & operator  *() {return static_cast<Element &>(Mesh::base_iterator::operator *());}
			__INLINE Element * operator ->() {return static_cast<Element *>(Mesh::base_iterator::operator ->());}
			friend iteratorElement Mesh::BeginElement(ElementType Types);
			friend iteratorElement Mesh::EndElement();
		};
		class iteratorSet : public base_iterator
		{
		private:
			iteratorSet(const base_iterator & other) :base_iterator(other) {}
		public:
			__INLINE ElementSet & operator  *() {return *m->sets[Number];}
			__INLINE ElementSet * operator ->() {return m->sets[Number];}
			friend iteratorSet Mesh::BeginSet();
			friend iteratorSet Mesh::EndSet();
		};
		class iteratorCell : public base_iterator
		{
		private:
			iteratorCell(const base_iterator & other) :base_iterator(other) {}
		public:
			__INLINE Cell & operator  *() {return *m->cells[Number];}
			__INLINE Cell * operator ->() {return  m->cells[Number];}
			friend iteratorCell Mesh::BeginCell();
			friend iteratorCell Mesh::EndCell();
		};
		class iteratorFace : public base_iterator
		{
		private:
			iteratorFace(const base_iterator & other) :base_iterator(other) {}
		public:
			__INLINE Face & operator  *() {return *m->faces[Number];}
			__INLINE Face * operator ->() {return  m->faces[Number];}
			friend iteratorFace Mesh::BeginFace();
			friend iteratorFace Mesh::EndFace();
		};
		class iteratorEdge : public base_iterator
		{
		private:
			iteratorEdge(const base_iterator & other) :base_iterator(other) {}
		public:
			__INLINE Edge & operator  *() {return *m->edges[Number];} 
			__INLINE Edge * operator ->() {return  m->edges[Number];}
			friend iteratorEdge Mesh::BeginEdge();
			friend iteratorEdge Mesh::EndEdge();
		};
		class iteratorNode : public base_iterator
		{
		private:
			iteratorNode(const base_iterator & other) :base_iterator(other) {}
		public:
			__INLINE Node & operator  *() {return *m->nodes[Number];}
			__INLINE Node * operator ->() {return  m->nodes[Number];}
			friend iteratorNode Mesh::BeginNode();
			friend iteratorNode Mesh::EndNode();
		};
	private:
		int          GetTypeEnd(ElementType t);
	private:
		std::vector< std::pair<std::string, std::string> > file_options;
		//implemented in io.hpp
		io_converter<INMOST_DATA_INTEGER_TYPE,INMOST_DATA_REAL_TYPE> iconv;
		io_converter<INMOST_DATA_ENUM_TYPE   ,INMOST_DATA_REAL_TYPE> uconv;
	private:
		//implemented in mesh_file.cpp
		void         WriteTag  (std::ostream & output, Tag X);
		Tag          ReadTag   (std::istream & input);
		void         ReadData  (std::istream & input, Tag t, Storage * e,std::vector<Node *> & new_nodes,std::vector<Edge *> & new_edges,std::vector<Face *> & new_faces,std::vector<Cell *> & new_cells);
		void         WriteData (std::ostream & output, Tag t, Storage & e);
		void         WriteElementSet(std::ostream & output, ElementSet * X);
		ElementSet * ReadElementSet(std::istream & input,std::vector<Node *> & new_nodes,std::vector<Edge *> & new_edges,std::vector<Face *> & new_faces,std::vector<Cell *> & new_cells);
		Node *       ReadNode(std::istream & input, std::vector<Node *> & old_nodes);
		void         WriteNode(std::ostream & output,Node * X);
		Edge *       ReadEdge(std::istream & input, std::vector<Node *> & nodes);
		void         WriteEdge(std::ostream & output,Edge * X);
		Face *       ReadFace(std::istream & input, std::vector<Edge *> & edges);
		void         WriteFace(std::ostream & output,Face * X);
		Cell *       ReadCell(std::istream & input,std::vector<Face *> & faces, std::vector<Node *> & nodes);
		void         WriteCell(std::ostream & output,Cell * X);
	public:
		/// Current availible file options:
		/// "VTK_GRID_DIMS" - set "2" for two-dimensional vtk grids, "3" for three-dimensional vtk grids
		/// "VERBOSITY"     - set "2" for progress messages, "1" for reports, "0" for silence
		void         SetFileOption(std::string,std::string);
		std::string  GetFileOption(std::string);
		void         Load(std::string File); // .vtk, .pvtk, .pmf
		/// Remeber: .pmf stores all references to elements. If reference are broken due to mesh modification,
		///          saving or loading such a mesh may lead to seagfault. To automatically maintain correct
		///          references modify mesh using BeginModification, ApplyModification, EndModification
		void         Save(std::string File); // .vtk, .pvtk, .gmv, .pmf
		bool         isParallelFileFormat(std::string File);
	public:
		
		//implemented in geometry.cpp
	private:
		void RestoreGeometricTags();
		Tag          measure_tag;
		Tag          centroid_tag;
		Tag          normal_tag;
		Tag          barycenter_tag;
		Tag          dot_tag;
		bool         remember[5][3];
		bool         HideGeometricData(GeometricData type, ElementType mask) {return remember[type][ElementNum(mask)-1] = false;}
		bool         ShowGeometricData(GeometricData type, ElementType mask) {return remember[type][ElementNum(mask)-1] = true;}
		//~ MarkerType reorient;
	public:
		typedef std::map<GeometricData, ElementType> GeomParam;
		// types for MEASURE:     EDGE | FACE | CELL   (length, area, volume)
		// types for CENTROID:    EDGE | FACE | CELL
		// types for BARYCENTER:  EDGE | FACE | CELL
		// types for NORMAL:      FACE | CELL          (may precompute normal for cells in 2d case)
		// types for ORIENTATION: FACE
		void PrepareGeometricData(GeomParam table);
		void RemoveGeometricData(GeomParam table);
		bool HaveGeometricData(GeometricData type, ElementType mask) const {return remember[type][ElementNum(mask)-1];} // requests to only one geometric and element type allowed
		void GetGeometricData(Element * e, GeometricData type, Storage::real * ret);
		void RecomputeGeometricData(Element * e); // Update all stored geometric data, runs automatically on element creation
		Tag GetGeometricTag(GeometricData type) {switch(type) {case MEASURE: return measure_tag; case CENTROID: return centroid_tag; case BARYCENTER: return barycenter_tag; case NORMAL: return normal_tag;} return Tag();}
		
		
		bool TestClosure(Element ** elements, unsigned num);

		std::vector<Face *> GatherBoundaryFaces();
		std::vector<Face *> GatherInteriorFaces();
		Storage::integer CountBoundaryFaces();
		Storage::integer CountInteriorFaces();
		Element::GeometricType ComputeGeometricType(ElementType element_type, Element ** lower_adjacent, INMOST_DATA_ENUM_TYPE lower_adjacent_size);
		//implemented in modify.cpp
	private:
		MarkerType hide_element, new_element;
	public:
		bool isMeshModified() const {return new_element != 0;} //In case mesh is modified, on element creation TieElements will always place elements to the end
		MarkerType HideMarker() const {return hide_element;}
		MarkerType NewMarker() const {return new_element;}
		void SwapModification(); // swap hidden and new elements, so that old mesh is recovered
		void BeginModification();  //allow elements to be hidden
		void ApplyModification();  //modify DATA_REFERENCE tags so that links to hidden elements are converted to NULL and removed from sets
		void ResolveModification(); //resolve parallel state of newly created elements, restore ghost layers; not implemented
		void EndModification();    //delete hidden elements
		
		
		//Elements and data is only moved, copy implementation from ReorderEmpty
		//~ void BeginMoveModification();
		//~ void EndMoveModification();
		
		static INMOST_DATA_ENUM_TYPE getNext(Element * const * arr, INMOST_DATA_ENUM_TYPE size, INMOST_DATA_ENUM_TYPE k, MarkerType marker);
		static INMOST_DATA_ENUM_TYPE Count(Element * const * arr, INMOST_DATA_ENUM_TYPE size, MarkerType marker);
		//implemented in mesh.cpp
	private:
		Tag            tag_topologyerror;
		TopologyCheck  checkset;
		TopologyCheck  errorset;
		TopologyCheck  BeginTopologyCheck(ElementType etype, Element ** adj, INMOST_DATA_ENUM_TYPE num); //check provided elements
		TopologyCheck  EndTopologyCheck(Element * e); //check created element
	public:
		Tag            TopologyErrorTag() const {return tag_topologyerror;}
		TopologyCheck  GetTopologyCheck(TopologyCheck mask = ENUMUNDEF) const {return checkset & mask;}
		void           SetTopologyCheck(TopologyCheck mask) {checkset = checkset | mask;}
		void           RemTopologyCheck(TopologyCheck mask) {checkset = checkset & ~mask;}

		void           SetTopologyError(TopologyCheck mask) {errorset = errorset | mask;}
		TopologyCheck  GetTopologyError(TopologyCheck mask = ENUMUNDEF) const {return errorset & mask;}
		void           ClearTopologyError(TopologyCheck mask = ENUMUNDEF) {errorset = errorset & ~mask;}

		friend class Storage;
	};

	template <typename AdjacentType>
	void adjacent<AdjacentType>::unite(const adjacent<AdjacentType> & other)
	{
		if( empty() ) insert(container.end(),other.container.begin(),other.container.end());
		else
		{
			INMOST_DATA_ENUM_TYPE s = size();
			MarkerType mrk = container[0]->GetMeshLink()->CreateMarker();
			for(INMOST_DATA_ENUM_TYPE it = 0; it < size(); it++) container[it]->SetMarker(mrk);
			for(INMOST_DATA_ENUM_TYPE it = 0; it < other.size(); it++) if( !other.container[it]->GetMarker(mrk) ) container.push_back(other.container[it]);
			for(INMOST_DATA_ENUM_TYPE it = 0; it < s; it++) container[it]->RemMarker(mrk);
			container[0]->GetMeshLink()->ReleaseMarker(mrk);
		}
	}
	template <typename AdjacentType>
	void adjacent<AdjacentType>::substract(const adjacent<AdjacentType> & other)
	{
		if( !empty() )
		{
			Mesh * m = container[0]->GetMeshLink();
			MarkerType mrk = m->CreateMarker();
			for(INMOST_DATA_ENUM_TYPE it = 0; it < other.size(); it++) other.container[it]->SetMarker(mrk);
			{
				INMOST_DATA_ENUM_TYPE m = 0, n = 0;
				while( m < size() ) 
				{
					if( !container[m]->GetMarker(mrk) )
						container[n++] = container[m];
					m++;
				}
				container.resize(n);
			}
			for(INMOST_DATA_ENUM_TYPE it = 0; it < other.size(); it++) other.container[it]->RemMarker(mrk);
			m->ReleaseMarker(mrk);
		}
	}
	template <typename AdjacentType>
	void adjacent<AdjacentType>::intersect(const adjacent<AdjacentType> & other)
	{
		if( !empty() )
		{
			Mesh * m = container[0]->GetMeshLink();
			MarkerType mrk = m->CreateMarker();
			for(INMOST_DATA_ENUM_TYPE it = 0; it < other.size(); it++) other.container[it]->SetMarker(mrk);
			{
				INMOST_DATA_ENUM_TYPE m = 0, n = 0;
				while( m < size() ) 
				{
					if( container[m]->GetMarker(mrk) )
						container[n++] = container[m];
					m++;
				}
				container.resize(n);
			}
			for(INMOST_DATA_ENUM_TYPE it = 0; it < other.size(); it++) other.container[it]->RemMarker(mrk);
			m->ReleaseMarker(mrk);
		}
	}
	
	bool LessElements(Element * a,Element * b); //function for std::sort and std::set based on USE_COMPARE
	int CompareElements(Element * a,Element * b); //function based on USE_COMPARE
	int CompareElementsC(const void * a,const void * b); //function for qsort based on USE_COMPARE 

	int CompareElementsGID(Element * a,Element * b);
	int CompareElementsPointer(Element * a,Element * b);
	int CompareElementsUnique(Element * a,Element * b);
	int CompareElementsCentroid(Element * a,Element * b);
	int CompareElementsHybrid(Element * a,Element * b);

	bool LessElementsGID(Element * a,Element * b);	
	bool LessElementsUnique(Element * a,Element * b);
	bool LessElementsPointer(Element * a,Element * b);
	bool LessElementsCentroid(Element * a,Element * b);
	
	int CompareElementsCGID(const void * pa, const void * pb);
	int CompareElementsCPointer(const void * pa, const void * pb);
	int CompareElementsCUnique(const void * pa, const void * pb);
	int CompareElementsCCentroid(const void * pa, const void * pb);
	
	
	//functions for algorithm cpp
	void SwapElement(void * pa, void * pb, void * udata);
	void SwapPointer(void * pa, void * pb, void * udata);
	int CompareElement(const void * pa, const void * pb, void * udata);
	int CompareElementPointer(const void * pa, const void * pb, void * udata);
	int CompareCoordSearch(const void * pa, void * udata);
	int CompareGIDSearch(const void * pa, void * udata);
	
	//algorithm.cpp
	void sort(void *pbase, unsigned num, unsigned width, int (*comp)(const void *, const void *,void *), void (*swap)(void *, void *, void *), void * user_data);
	int binary_search(void * data, unsigned num, unsigned width, int comp(const void *,void *), void * user_data);
	size_t hashfunci32(const char * pos, int n);



}

#endif

#endif // INMOST_MESH_H_INCLUDED
