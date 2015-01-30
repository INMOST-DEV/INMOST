
#ifndef INMOST_MESH_H_INCLUDED
#define INMOST_MESH_H_INCLUDED

#include "inmost_common.h"

#if defined(USE_MESH)

#define __NDT 3
#define __NET 6


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
	class ElementSet;

	typedef INMOST_DATA_BULK_TYPE ElementType;
	static const ElementType NONE = 0x00;
	static const ElementType NODE = 0x01;
	static const ElementType EDGE = 0x02;
	static const ElementType FACE = 0x04;
	static const ElementType CELL = 0x08;
	static const ElementType ESET = 0x10;
	static const ElementType MESH = 0x20;

	typedef INMOST_DATA_ENUM_TYPE HandleType;
	static const INMOST_DATA_ENUM_TYPE handle_etype_bits  = 3;
	static const INMOST_DATA_ENUM_TYPE handle_etype_shift = sizeof(HandleType)*8-handle_etype_bits;
	static const INMOST_DATA_ENUM_TYPE handle_id_mask     = (1 << handle_etype_shift)-1;
	
	typedef INMOST_DATA_BULK_TYPE GeometricData;
	static const GeometricData CENTROID     = 0;
	static const GeometricData NORMAL       = 1;
	static const GeometricData ORIENTATION  = 2;
	static const GeometricData MEASURE      = 3;
	static const GeometricData BARYCENTER   = 4;
	
	typedef INMOST_DATA_BULK_TYPE SyncBitOp; //< This type is used for marker synchronization
	static const SyncBitOp SYNC_BIT_SET = 0;
	static const SyncBitOp SYNC_BIT_OR  = 1;
	static const SyncBitOp SYNC_BIT_XOR = 2;
	static const SyncBitOp SYNC_BIT_AND = 3;

	typedef array<INMOST_DATA_REAL_TYPE>    inner_real_array;
	typedef array<INMOST_DATA_INTEGER_TYPE> inner_integer_array;
	typedef array<INMOST_DATA_BULK_TYPE>    inner_bulk_array;
	typedef array<HandleType>               inner_reference_array;

	enum DataType
	{
		DATA_REAL      = 0, 
		DATA_INTEGER   = 1, 
		DATA_BULK      = 2,
		DATA_REFERENCE = 3
	};

	typedef INMOST_DATA_ENUM_TYPE MarkerType; // low 8 bits - marker mask, rest high bits - position of marker
	static const INMOST_DATA_ENUM_TYPE MarkerFields = 16;   // number of chars to hold all markers, total number (MarkerFields * bits_per_char)
	static const INMOST_DATA_ENUM_TYPE MarkerMask   = static_cast<INMOST_DATA_BULK_TYPE>(-1); // bit mask to obtain marker mask within MarkerType
	static const INMOST_DATA_ENUM_TYPE MarkerShift  = sizeof(INMOST_DATA_BULK_TYPE)*8;    // sizeof(char) * bits_per_char
	__INLINE static MarkerType InvalidMarker() {return ENUMUNDEF;}


	static const INMOST_DATA_ENUM_TYPE chunk_bits_elems       = 13;
	static const INMOST_DATA_ENUM_TYPE chunk_bits_empty       = 8;
	static const INMOST_DATA_ENUM_TYPE chunk_bits_tags        = 6;
	static const INMOST_DATA_ENUM_TYPE chunk_bits_dense_array = 6;

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
	static const TopologyCheck ADJACENT_VALID         = 0x00020000; //done//check that all handles are valid
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
	static const TopologyCheck DEFAULT_CHECK          = THROW_EXCEPTION | DUPLICATE_EDGE | DUPLICATE_FACE | DUPLICATE_CELL | PRINT_NOTIFY;
	const char *                             TopologyCheckNotifyString(TopologyCheck c); //mesh.cpp
	const char *                             DataTypeName         (DataType t); //tag.cpp
	__INLINE bool                            OneType              (ElementType t) {return t > 0 && (t & (t-1)) == 0;}
	__INLINE ElementType                     FirstElementType     () {return NODE;}
	__INLINE ElementType                     LastElementType      () {return MESH << 1;}
	__INLINE ElementType                     NextElementType      (ElementType etype) {return etype << 1;}
	__INLINE ElementType                     PrevElementType      (ElementType etype) {return etype >> 1;}
	__INLINE ElementType                     ElementTypeFromDim   (INMOST_DATA_INTEGER_TYPE dim) {return 1 << dim;}
	const char *                             ElementTypeName      (ElementType t); //mesh.cpp
	__INLINE INMOST_DATA_INTEGER_TYPE        ElementNum           (ElementType t)
	{
		unsigned int v = static_cast<unsigned int>(t);  // 32-bit value to find the log2 of 
		//static const unsigned int b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 0xFF00FF00, 0xFFFF0000};
		static const unsigned int b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0};
		register unsigned int r = (v & b[0]) != 0;
		//r |= ((v & b[4]) != 0) << 4;
		//r |= ((v & b[3]) != 0) << 3;
		r |= ((v & b[2]) != 0) << 2;
		r |= ((v & b[1]) != 0) << 1;
		return static_cast<INMOST_DATA_INTEGER_TYPE>(r);
	}
	__INLINE HandleType                      InvalidHandle        () {return 0;}
	__INLINE INMOST_DATA_INTEGER_TYPE        GetHandleID          (HandleType h) {return (h & handle_id_mask)-1;}
	__INLINE INMOST_DATA_INTEGER_TYPE        GetHandleElementNum  (HandleType h) {return h >> handle_etype_shift;}
	__INLINE ElementType                     GetHandleElementType (HandleType h) {return 1 << GetHandleElementNum(h);}
	__INLINE HandleType                      ComposeHandle        (ElementType etype, INMOST_DATA_INTEGER_TYPE ID) {return ID == -1 ? InvalidHandle() : ((ElementNum(etype) << handle_etype_shift) + (1+ID));}
	__INLINE HandleType                      ComposeHandle        (INMOST_DATA_INTEGER_TYPE etypenum, INMOST_DATA_INTEGER_TYPE ID) {return ID == -1 ? InvalidHandle() : ((etypenum << handle_etype_shift) + (1+ID));}
	__INLINE bool                            isValidHandle        (HandleType h) {return h != 0;}
	

	

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
		INMOST_DATA_ENUM_TYPE  size, bytes_size;
		bool                   sparse[__NET];
		INMOST_DATA_ENUM_TYPE  record_size;
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
		__INLINE INMOST_DATA_ENUM_TYPE  GetPositionByDim(INMOST_DATA_ENUM_TYPE typenum) const {return mem->pos[typenum];}
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
		__INLINE INMOST_DATA_ENUM_TYPE  GetBytesSize() const {assert(mem!=NULL); return mem->bytes_size;}
		__INLINE INMOST_DATA_ENUM_TYPE  GetSize() const {assert(mem!=NULL); return mem->size;}
		__INLINE std::string            GetTagName() const {assert(mem!=NULL); return mem->tagname;}
		__INLINE bool                   isDefined(ElementType type) const {assert(mem!=NULL && OneType(type)); return GetPosition(type) != ENUMUNDEF;}
		__INLINE bool                   isSparse(ElementType type) const {assert(mem!=NULL && OneType(type)); return mem->sparse[ElementNum(type)];}
		__INLINE bool                   isValid() const {return mem != NULL;}
		__INLINE Mesh *                 GetMeshLink() const {assert(mem!=NULL); return mem->m_link;}		
		__INLINE bool                   isSparseByDim(INMOST_DATA_INTEGER_TYPE typenum) const {assert(mem!=NULL); return mem->sparse[typenum];}
		__INLINE bool                   isDefinedByDim(INMOST_DATA_INTEGER_TYPE typenum) const {assert(mem!=NULL); return GetPositionByDim(typenum) != ENUMUNDEF;}
		__INLINE void                   SetBulkDataType(INMOST_MPI_Type type) {assert(mem!=NULL && mem->dtype == DATA_BULK ); mem->bulk_data_type = type;}
		friend class TagManager;
		friend class Storage;
		friend class Mesh;
	};
	
	class TagManager //implemented in tag.cpp
	{
	protected:
		TagManager();
		TagManager(const TagManager & other);
		TagManager & operator = (TagManager const & other);
		typedef chunk_array<INMOST_DATA_ENUM_TYPE,chunk_bits_empty>    empty_data;
		typedef chunk_array<Tag, chunk_bits_tags>                      tag_array_type;
		typedef chunk_bulk_array<chunk_bits_elems>                     dense_sub_type;
		typedef chunk_array<dense_sub_type,chunk_bits_dense_array>     dense_data_array_type;
		typedef struct{void * tag, * rec;}                             sparse_sub_record;
		typedef array< sparse_sub_record >                             sparse_sub_type;
		typedef chunk_array< sparse_sub_type,chunk_bits_elems>         sparse_data_array_type;
		typedef chunk_array<INMOST_DATA_INTEGER_TYPE,chunk_bits_elems> back_links_type;
		//typedef std::vector<INMOST_DATA_INTEGER_TYPE>                  back_links_type;
	public:
		typedef tag_array_type::iterator                               iteratorTag;
	public:
		virtual ~TagManager();
		bool                             HaveTag            (std::string name) const;
		Tag                              GetTag             (std::string name) const;
		void                             ListTagNames       (std::vector<std::string> & list) const;
		Tag                              CreateTag          (Mesh * m, std::string name, DataType dtype, ElementType etype, ElementType sparse, INMOST_DATA_ENUM_TYPE size = ENUMUNDEF); 
		virtual Tag                      DeleteTag          (Tag tag, ElementType mask); 
		bool                             ElementDefined     (Tag const & tag, ElementType etype) const;
	protected:
		void                             ReallocateData     (const Tag & t, INMOST_DATA_INTEGER_TYPE etypenum,INMOST_DATA_ENUM_TYPE new_size);
		void                             ReallocateData     (INMOST_DATA_INTEGER_TYPE etypenum, INMOST_DATA_ENUM_TYPE new_size);
		__INLINE sparse_sub_type const & GetSparseData      (int etypenum, int local_id) const {return sparse_data[etypenum][local_id];}
		__INLINE sparse_sub_type &       GetSparseData      (int etypenum, int local_id) {return sparse_data[etypenum][local_id];}
		__INLINE dense_sub_type const &  GetDenseData       (int pos) const {return dense_data[pos];}
		__INLINE dense_sub_type &        GetDenseData       (int pos) {return dense_data[pos];}
		static void                      CopyData           (const Tag & t, void * adata, const void * bdata);
		static void                      DestroyVariableData(const Tag & t, void * adata);
	protected:
		typedef tag_array_type::iterator            tag_iterator;
		typedef tag_array_type::const_iterator      tag_const_iterator;
	protected:
		tag_array_type                              tags;
		empty_data                                  empty_dense_data;
		dense_data_array_type                       dense_data;
		sparse_data_array_type                      sparse_data[6];
		back_links_type                             back_links[6];
	};
	
	/// Base class for Mesh, Element, and ElementSet classes.
	/// This base class is used for the Mesh class, as well as Element classes, and ElementSet class.
	/// 
	/// Storage class is used for representing different data objects in memory.
	/// Each data object is associated with corresponding Tag.
	class Storage //implemented in storage.cpp
	{
	public:
		/// Storage type for representing real values.
		typedef INMOST_DATA_REAL_TYPE     real; 
		/// Storage type for representing integer values.
		typedef INMOST_DATA_INTEGER_TYPE  integer;
		/// Storage type for representing one byte of abstact data.
		typedef INMOST_DATA_BULK_TYPE     bulk;
		/// type for representing unsigned integer values.
		typedef INMOST_DATA_ENUM_TYPE     enumerator;
		/// Storage type for representing references to Element.
		typedef HandleType                reference;
		/// Storage type for representing arrays of real values.
		typedef shell<real>               real_array;
		/// Storage type for representing arrays of integer values.
		typedef shell<integer>            integer_array;
		/// Storage type for representing abstact data as a series of bytes.
		typedef shell<bulk>               bulk_array;
		/// Storage type for representing arrays of Element references.
		class reference_array : public shell<HandleType>
		{
			Mesh * m;
		public:
			reference_array() : m(NULL) {}
			reference_array(Mesh * m, inner_reference_array & arr) : shell<reference>(arr), m(m) {} 
			reference_array(Mesh * m, reference * arr, size_type size) : shell<reference>(arr,size), m(m) {}
			~reference_array() {}
			void push_back(const Storage & elem);
			void push_back(HandleType h) {shell<reference>::push_back(h);} //is it needed?
			Element operator[] (size_type n);
			Element operator[] (size_type n) const;
			class iterator : public shell<HandleType>::iterator
			{
				Mesh * m;
			public:
				iterator(Mesh * m, const shell<HandleType>::iterator & other) : shell<HandleType>::iterator(other), m(m) {}
				iterator(const iterator & other) : shell<HandleType>::iterator(other), m(other.m) {}
				iterator & operator =(iterator const & other) {m = other.m; shell<HandleType>::iterator::operator=(other); return *this;}
				iterator & operator ++() {shell<HandleType>::iterator::operator++(); return *this;}
				iterator   operator ++(int) {iterator ret(*this); shell<HandleType>::iterator::operator++(); return ret;}
				iterator & operator --() {shell<HandleType>::iterator::operator--(); return *this;}
				iterator   operator --(int) {iterator ret(*this); shell<HandleType>::iterator::operator--(); return ret;}
				Element operator->();
			};
			class const_iterator : public shell<HandleType>::const_iterator
			{
				Mesh * m;
			public:
				const_iterator(Mesh * m, const shell<HandleType>::const_iterator & other) : shell<HandleType>::const_iterator(other) , m(m) {}
				const_iterator(const const_iterator & other) : shell<HandleType>::const_iterator(other), m(other.m) {}
				const_iterator & operator =(const_iterator const & other) {m = other.m; shell<HandleType>::const_iterator::operator=(other); return *this;}
				const_iterator & operator ++() {shell<HandleType>::const_iterator::operator++(); return *this;}
				const_iterator   operator ++(int) {const_iterator ret(*this); shell<HandleType>::const_iterator::operator++(); return ret;}
				const_iterator & operator --() {shell<HandleType>::const_iterator::operator--(); return *this;}
				const_iterator   operator --(int) {const_iterator ret(*this); shell<HandleType>::const_iterator::operator--(); return ret;}
				Element operator->();
			};
			class reverse_iterator : public shell<HandleType>::reverse_iterator
			{
				Mesh * m;
			public:
				reverse_iterator(Mesh * m, const shell<HandleType>::reverse_iterator & other) : shell<HandleType>::reverse_iterator(other), m(m) {}
				reverse_iterator(const reverse_iterator & other) : shell<HandleType>::reverse_iterator(other), m(other.m)  {}
				reverse_iterator & operator =(reverse_iterator const & other) {m = other.m; shell<HandleType>::reverse_iterator::operator=(other); return *this;}
				reverse_iterator & operator ++() {shell<HandleType>::reverse_iterator::operator++(); return *this;}
				reverse_iterator   operator ++(int) {reverse_iterator ret(*this); shell<HandleType>::reverse_iterator::operator++(); return ret;}
				reverse_iterator & operator --() {shell<HandleType>::reverse_iterator::operator--(); return *this;}
				reverse_iterator   operator --(int) {reverse_iterator ret(*this); shell<HandleType>::reverse_iterator::operator--(); return ret;}
				Element operator->();
			};
			class const_reverse_iterator : public shell<HandleType>::const_reverse_iterator
			{
				Mesh * m;
			public:
				const_reverse_iterator(Mesh * m, const shell<HandleType>::const_reverse_iterator & other) : shell<HandleType>::const_reverse_iterator(other), m(m) {}
				const_reverse_iterator(const const_reverse_iterator & other) : shell<HandleType>::const_reverse_iterator(other), m(other.m) {}
				const_reverse_iterator & operator =(const_reverse_iterator const & other) {m = other.m; shell<HandleType>::const_reverse_iterator::operator=(other); return *this;}
				const_reverse_iterator & operator ++() {shell<HandleType>::const_reverse_iterator::operator++(); return *this;}
				const_reverse_iterator   operator ++(int) {const_reverse_iterator ret(*this); shell<HandleType>::const_reverse_iterator::operator++(); return ret;}
				const_reverse_iterator & operator --() {shell<HandleType>::const_reverse_iterator::operator--(); return *this;}
				const_reverse_iterator   operator --(int) {const_reverse_iterator ret(*this); shell<HandleType>::const_reverse_iterator::operator--(); return ret;}
				Element operator->();
			};
			iterator begin() {return iterator(m,shell<HandleType>::begin());}
			const_iterator begin() const {return const_iterator(m,shell<HandleType>::begin());}
			reverse_iterator rbegin() {return reverse_iterator(m,shell<HandleType>::rbegin());}
			const_reverse_iterator rbegin() const {return const_reverse_iterator(m,shell<HandleType>::rbegin());}
			iterator end() {return iterator(m,shell<HandleType>::end());}
			const_iterator end() const {return const_iterator(m,shell<HandleType>::end());}
			reverse_iterator rend() {return reverse_iterator(m,shell<HandleType>::rend());}
			const_reverse_iterator rend() const {return const_reverse_iterator(m,shell<HandleType>::rend());}
		};
		//typedef shell<reference>          reference_array;
	protected:
		HandleType                        handle;
		HandleType *                      handle_link;
	private:
		Mesh *                            m_link;
	public:
		Storage(const Storage & other) : handle(other.handle), handle_link(other.handle_link), m_link(other.m_link) {}
		Storage(Mesh * mesh, HandleType handle) : handle(handle), handle_link(NULL), m_link(mesh) {}
		/// This constructor allows for remote handle modification
		Storage(Mesh * mesh, HandleType * handle) : handle(*handle), handle_link(handle), m_link(mesh) {}
		/// If there is a link to handle provided (automatically by ElementArray and reference_array),
		/// then remote handle value will be modified
		Storage & operator =(Storage const & other); 
		bool      operator ==(const Storage & other) {return handle == other.handle;} //m_link may be NULL if Invalidxxx is used
		bool      operator !=(const Storage & other) {return handle != other.handle;} //m_link may be NULL if Invalidxxx is used
		Storage * operator ->() {return this;}
		const Storage * operator->() const {return this;}
		Storage & self() {return *this;}
		const Storage & self() const {return *this;}
		virtual ~Storage() {}
	public:
		/// Retrieve real value associated with Tag.
		real      &              Real            (const Tag & tag) const;
		/// Retrieve integer value associated with Tag.
		integer   &              Integer         (const Tag & tag) const;
		/// Retrieve one byte of abstract data associated with Tag.
		bulk      &              Bulk            (const Tag & tag) const;
		/// Retrieve Element reference associated with Tag.
		reference &              Reference       (const Tag & tag) const;
		/// Retrieve array of real values associated with Tag.
		real_array               RealArray       (const Tag & tag) const;
		/// Retrieve array of integer values associated with Tag.
		integer_array            IntegerArray    (const Tag & tag) const;
		/// Retrieve abstract data associated with Tag as a series of bytes.
		bulk_array               BulkArray       (const Tag & tag) const;
		/// Retrieve array of Element references associated with Tag.
		reference_array          ReferenceArray  (const Tag & tag) const;
		
		//optimized data requests for dense data with fixed size
		real_array               RealArrayDF     (const Tag & tag) const;
		integer_array            IntegerArrayDF  (const Tag & tag) const;
		bulk_array               BulkArrayDF     (const Tag & tag) const;
		reference_array          ReferenceArrayDF(const Tag & tag) const;
		real      &              RealDF          (const Tag & tag) const;
		integer   &              IntegerDF       (const Tag & tag) const;
		bulk      &              BulkDF          (const Tag & tag) const;
		reference &              ReferenceDF     (const Tag & tag) const;
		
		//optimized data requests for dense data with variable size
		real_array               RealArrayDV     (const Tag & tag) const;
		integer_array            IntegerArrayDV  (const Tag & tag) const;
		bulk_array               BulkArrayDV     (const Tag & tag) const;
		reference_array          ReferenceArrayDV(const Tag & tag) const;
		real      &              RealDV          (const Tag & tag) const;
		integer   &              IntegerDV       (const Tag & tag) const;
		bulk      &              BulkDV          (const Tag & tag) const;
		reference &              ReferenceDV     (const Tag & tag) const;
		
		/// Return the data length associated with Tag.
		/// For abstract data return the number of bytes, otherwise return the length of associated array. 
		/// @see Storage::SetDataSize
		INMOST_DATA_ENUM_TYPE    GetDataSize     (const Tag & tag) const; //For DATA_BULK return number of bytes, otherwise return the length of array
		/// Set the length of  data associated with Tag.
		/// @param tag Identifying Tag.
		/// @param new_size The number of bytes for abstract data, otherwise the length of the array.
		/// @see Storage::GetDataSize
		void                     SetDataSize     (const Tag & tag,INMOST_DATA_ENUM_TYPE new_size) const;
		/// Extract part of the data associated with Tag.
		/// Copy part of the associated array or data to the destination memory.
		/// @param tag Identifying Tag.
		/// @param shift Starting position of the copied data.
		/// For abstact data – number of bytes to skip, otherwise number of values to skip.
		/// @param size Number of elements to copy.
		/// For abstact data – number of bytes to copy, otherwise number of values to copy.
		/// @param data Destination position to copy data to.
		/// @see Storage::SetData
		void                     GetData         (const Tag & tag, INMOST_DATA_ENUM_TYPE shift, INMOST_DATA_ENUM_TYPE size, void * data) const;
		void                     SetData         (const Tag & tag, INMOST_DATA_ENUM_TYPE shift, INMOST_DATA_ENUM_TYPE size, const void * data) const;
		void                     DelData         (const Tag & tag) const;
		/// Deallocates space allocated for sparse data, frees variable array if necessary
		void                     DelSparseData   (const Tag & tag) const;
		/// Frees variable array or fills field with zeroes
		void                     DelDenseData    (const Tag & tag) const;
		/// Check if any data is associated with Tag.
		bool                     HaveData        (const Tag & tag) const;
		__INLINE ElementType     GetElementType  () const {return GetHandleElementType(handle);}
		__INLINE integer         GetElementNum   () const {return GetHandleElementNum(handle);}
		void                     SetMarker       (MarkerType n) const;
		bool                     GetMarker       (MarkerType n) const;
		void                     RemMarker       (MarkerType n) const;
		void                     ClearMarkerSpace() const;
		void                     GetMarkerSpace  (bulk copy[MarkerFields]) const;
		void                     SetMarkerSpace  (bulk source[MarkerFields]) const;
		__INLINE integer         LocalID         () const {return GetHandleID(handle);}
		/// This number is guaranteed to be between 0 and Mesh::NumberOf(type of element)
		/// after Mesh::ReorderEmpty
		integer                  DataLocalID     () const;
		bool                     isValid         () const;
		__INLINE Mesh *          GetMeshLink     () const {return m_link;}
		__INLINE HandleType      GetHandle       () const {return handle;}
		friend class Mesh;
	};
	
	template <typename StorageType>
	class ElementArray
	{
	public:
		typedef dynarray<HandleType,64>      cont_t;
		typedef cont_t::size_type            size_type;
	private:
		Mesh *                               m_link;
		cont_t                               container;
	public:
		const cont_t & GetContainer() {return container;}
		~ElementArray() {}
		ElementArray() : m_link(NULL) {}
		ElementArray(Mesh * m_link) : m_link(m_link) {}
		ElementArray(Mesh * m_link, size_type n, HandleType h = InvalidHandle())  : m_link(m_link), container(n,h) {}
		ElementArray(Mesh * m_link, const cont_t & c)  : m_link(m_link), container(c) {}
		ElementArray(const ElementArray & other) : m_link(other.m_link), container(other.container) {}
		template<class InputIterator>
		ElementArray(Mesh * m_link, InputIterator first, InputIterator last) :m_link(m_link), container(first,last) {}
		ElementArray & operator=(ElementArray const & other) {m_link = other.m_link; container = other.container; return *this;}
	public:
		class iterator : public cont_t::iterator
		{
			Mesh * m_link;
		public:
			iterator(Mesh * m, const cont_t::iterator & other ) : cont_t::iterator(other) , m_link(m){assert(m_link);}
			iterator(const iterator & other ) : cont_t::iterator(other), m_link(other.m_link) {assert(m_link);}
			iterator &    operator ++() {cont_t::iterator::operator++(); return *this;}
			iterator      operator ++(int) {iterator ret(*this); cont_t::iterator::operator++(); return ret;}
			iterator &    operator --() {cont_t::iterator::operator--(); return *this;}
			iterator      operator --(int) {iterator ret(*this); cont_t::iterator::operator--(); return ret;}
			iterator &    operator =(iterator const & other) {m_link = other.m_link; cont_t::iterator::operator=(static_cast<cont_t::iterator const &>(other)); return *this; }
			HandleType &  operator *() { return cont_t::iterator::operator *(); }
			StorageType   operator->() { return StorageType(m_link,&cont_t::iterator::operator *()); }
		};
		class reverse_iterator : public cont_t::reverse_iterator
		{
			Mesh * m_link;
		public:
			reverse_iterator(Mesh * m, const cont_t::reverse_iterator & other ) : cont_t::reverse_iterator(other), m_link(m) {assert(m_link);}
			reverse_iterator(const reverse_iterator & other ) : cont_t::reverse_iterator(other), m_link(other.m_link) {assert(m_link);}
			reverse_iterator &    operator ++() {cont_t::reverse_iterator::operator++(); return *this;}
			reverse_iterator      operator ++(int) {reverse_iterator ret(*this); cont_t::reverse_iterator::operator++(); return ret;}
			reverse_iterator &    operator --() {cont_t::reverse_iterator::operator--(); return *this;}
			reverse_iterator      operator --(int) {reverse_iterator ret(*this); cont_t::reverse_iterator::operator--(); return ret;}
			reverse_iterator & operator =(reverse_iterator const & other) {m_link = other.m_link; cont_t::reverse_iterator::operator=(static_cast<cont_t::reverse_iterator const &>(other)); return *this; }
			HandleType &       operator *() { return cont_t::reverse_iterator::operator *(); }
			StorageType        operator->() { return StorageType(m_link,&cont_t::reverse_iterator::operator *()); }
		};
		class const_iterator : public cont_t::const_iterator
		{
			Mesh * m_link;
		public:
			const_iterator(Mesh * m, const cont_t::const_iterator & other ) : cont_t::const_iterator(other), m_link(m) {assert(m_link);}
			const_iterator(const const_iterator & other ) : cont_t::const_iterator(other), m_link(other.m_link) {assert(m_link);}
			const_iterator &    operator ++() {cont_t::const_iterator::operator++(); return *this;}
			const_iterator      operator ++(int) {const_iterator ret(*this); cont_t::const_iterator::operator++(); return ret;}
			const_iterator &    operator --() {cont_t::const_iterator::operator--(); return *this;}
			const_iterator      operator --(int) {const_iterator ret(*this); cont_t::const_iterator::operator--(); return ret;}
			const_iterator &    operator =(const_iterator const & other) {m_link = other.m_link; cont_t::const_iterator::operator=(static_cast<cont_t::const_iterator const &>(other)); return *this; }
			const HandleType &  operator *() { return cont_t::const_iterator::operator *(); }
			StorageType         operator->() { return StorageType(m_link,cont_t::const_iterator::operator *()); }
		};
		class const_reverse_iterator : public cont_t::const_reverse_iterator
		{
			Mesh * m_link;
		public:
			const_reverse_iterator(Mesh * m, const cont_t::const_reverse_iterator & other) : cont_t::const_reverse_iterator(other), m_link(m) {assert(m_link);}
			const_reverse_iterator(const const_reverse_iterator & other ) : cont_t::const_reverse_iterator(other), m_link(other.m_link) {assert(m_link);}
			const_reverse_iterator &    operator ++() {cont_t::const_reverse_iterator::operator++(); return *this;}
			const_reverse_iterator      operator ++(int) {const_reverse_iterator ret(*this); cont_t::const_reverse_iterator::operator++(); return ret;}
			const_reverse_iterator &    operator --() {cont_t::const_reverse_iterator::operator--(); return *this;}
			const_reverse_iterator      operator --(int) {const_reverse_iterator ret(*this); cont_t::const_reverse_iterator::operator--(); return ret;}
			const_reverse_iterator & operator =(const_reverse_iterator const & other) { cont_t::const_reverse_iterator::operator=(static_cast<cont_t::const_reverse_iterator const &>(other)); return *this; }
			const HandleType &       operator *() { return cont_t::const_reverse_iterator::operator *(); }
			StorageType              operator->() { return StorageType(m_link,cont_t::const_reverse_iterator::operator *()); }
		};
	public:
		template<class InputIterator>
		__INLINE void             insert      (iterator pos,InputIterator pbeg, InputIterator pend) {container.insert(static_cast<cont_t::iterator>(pos),pbeg,pend);}
		__INLINE iterator         erase       (iterator pos) {return iterator(m_link,container.erase(static_cast<cont_t::iterator>(pos)));}

		__INLINE iterator         begin       () { return iterator(m_link,container.begin()); }
		__INLINE iterator         end         () { return iterator(m_link,container.end()); }
		__INLINE reverse_iterator rbegin      () { return reverse_iterator(m_link,container.rbegin()); }
		__INLINE reverse_iterator rend        () { return reverse_iterator(m_link,container.rend()); }
		__INLINE const_iterator         begin () const { return const_iterator(m_link,container.begin()); }
		__INLINE const_iterator         end   () const { return const_iterator(m_link,container.end()); }
		__INLINE const_reverse_iterator rbegin() const { return const_reverse_iterator(m_link,container.rbegin()); }
		__INLINE const_reverse_iterator rend  () const { return const_reverse_iterator(m_link,container.rend()); }
		
		__INLINE StorageType      operator [] (size_type n) {assert(m_link); return StorageType(m_link,&container[n]);}
		__INLINE StorageType      operator [] (size_type n) const {assert(m_link); return StorageType(m_link,container[n]);}
		__INLINE StorageType      front       () {assert(m_link); return StorageType(m_link,&container.front()); }
		__INLINE StorageType      front       () const {assert(m_link); return StorageType(m_link,container.front()); }
		__INLINE StorageType      back        ()  {assert(m_link); return StorageType(m_link,&container.back()); }
		__INLINE StorageType      back        () const {assert(m_link); return StorageType(m_link,container.back()); }
		__INLINE HandleType       atfront     () const { return container.front(); }
		__INLINE HandleType       atback      () const { return container.back(); }
		__INLINE HandleType &     atfront     () { return container.front(); }
		__INLINE HandleType &     atback      () { return container.back(); }
		__INLINE HandleType &     at          (size_type n) { return container.at(n); }
		__INLINE HandleType       at          (size_type n) const { return container.at(n); }
		__INLINE void             swap        (ElementArray<StorageType> & other) {Mesh * t = m_link; m_link = other.m_link; other.m_link = t; container.swap(other.container);}
		__INLINE void             push_back   (const Storage & x) {if( m_link == NULL ) m_link = x.GetMeshLink(); assert(x.GetMeshLink() == m_link); container.push_back(x.GetHandle());}
		//__INLINE void             push_back   (const StorageType & x) {container.push_back(x.GetHandle());}
		__INLINE void             push_back   (HandleType x) {assert(m_link != NULL); container.push_back(x);}
		__INLINE void             pop_back    () {container.pop_back();}
		__INLINE void             resize      (size_type n) {container.resize(n);}
		__INLINE bool             empty       () const {return container.empty();}
		__INLINE void             clear       () {container.clear();}
		__INLINE void             reserve     (size_type n) {container.reserve(n);}
		__INLINE size_type        size        () const {return container.size(); }
		__INLINE HandleType *     data        () {return container.data();}
		__INLINE const HandleType*data        () const {return container.data();}
		__INLINE Mesh *           GetMeshLink () const {assert(m_link); return m_link;}
		__INLINE void             SetMeshLink (Mesh * m) {m_link = m;}
		//implemented in mesh_array.cpp
		void                      Unite       (const HandleType * h, INMOST_DATA_ENUM_TYPE num);
		void                      Subtract    (const HandleType * h, INMOST_DATA_ENUM_TYPE num);
		void                      Intersect   (const HandleType * h, INMOST_DATA_ENUM_TYPE num);
		template<typename EType>
		void                      Unite       (const ElementArray<EType> & other) {Unite(other.data(),static_cast<INMOST_DATA_ENUM_TYPE>(other.size()));}
		template<typename EType>
		void                      Subtract    (const ElementArray<EType> & other) {Subtract(other.data(),static_cast<INMOST_DATA_ENUM_TYPE>(other.size()));}
		template<typename EType>
		void                      Intersect   (const ElementArray<EType> & other) {Intersect(other.data(),static_cast<INMOST_DATA_ENUM_TYPE>(other.size()));}
		void                      SetMarker   (MarkerType n) const;
		void                      RemMarker   (MarkerType n) const;
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
		static const GeometricType Set          = 253;
		static const GeometricType MeshPart     = 254;
		static const GeometricType MaxType      = 255;
		//enum GeometricType {Unset,Vertex,Line,MultiLine,Tri,Quad,Polygon,MultiPolygon,Tet,Hex,Prism,Pyramid,Polyhedron,Set};
		static const char *       GeometricTypeName(GeometricType t);
		static integer            GetGeometricDimension(GeometricType m_type);
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
	public:
		//adj_type &                  HighConn                () const;
		//adj_type &                  LowConn                 () const;
	protected:
		void                        SetGeometricType        (GeometricType t);
	public:
		Element() : Storage(NULL,InvalidHandle()) {}
		Element(Mesh * m, HandleType h) : Storage(m,h) {}
		Element(Mesh * m, HandleType * h) : Storage(m,h) {}
		Element(const Element & other) : Storage(other) {}
		Element & operator =(Element const & other) {Storage::operator =(other); return *this;}
		Element * operator ->() {return this;}
		const Element * operator ->() const {return this;}
		Element & self() {return *this;}
		const Element & self() const {return *this;}
		virtual ~Element() {}
	public:
		/// Retrive number of adjacent elements.
		/// For etype you can either pass one type as CELL,
		/// or several types as bitwise mask: NODE | CELL.
		/// @param etype bitwise mask of element types
		/// @return number of adjacent elements.
		virtual enumerator                nbAdjElements           (ElementType etype) const;
		/// Retrive unordered array of adjacent elements.
		/// If you care about orderness of nodes for face you should you Face::getNodes() instead.
		/// If you want faster access you may use direct access to handles stored in memory
		/// through Mesh::HighConn for upper adjacencies and Mesh::LowConn for lower adjacencies.
		/// @param etype bitwise mask of element types
		/// @return array of elements
		virtual ElementArray<Element>     getAdjElements          (ElementType etype) const;  //unordered
		/// Retrive number of adjacent elements with marker.
		/// As etype you can either pass one type as CELL,
		/// or several types as bitwise mask: NODE | CELL
		/// @param etype bitwise mask of element types
		/// @param mask marker to be set 
		/// @param invert_mask if true then those are selected on wich marker is not set
		/// @return number of adjacent elements.
		virtual enumerator                nbAdjElements           (ElementType etype, MarkerType mask, bool invert_mask = false) const;
		/// Retrive unordered array of adjacent elements with marker.
		/// @param etype bitwise mask of element types
		/// @param mask marker to be set 
		/// @param invert_mask if true then those are selected on wich marker is not set
		/// @return array of elements
		virtual ElementArray<Element>     getAdjElements          (ElementType etype, MarkerType mask, bool invert_mask = false) const;  //unordered
		ElementArray<Element>       BridgeAdjacencies       (ElementType Bridge, ElementType Dest, MarkerType mask = 0, bool invert_mask = false) const;
		ElementArray<Node>          BridgeAdjacencies2Node  (ElementType Bridge, MarkerType mask = 0, bool invert_mask = false) const;
		ElementArray<Edge>          BridgeAdjacencies2Edge  (ElementType Bridge, MarkerType mask = 0, bool invert_mask = false) const;
		ElementArray<Face>          BridgeAdjacencies2Face  (ElementType Bridge, MarkerType mask = 0, bool invert_mask = false) const;
		ElementArray<Cell>          BridgeAdjacencies2Cell  (ElementType Bridge, MarkerType mask = 0, bool invert_mask = false) const;
		/// Retrive all the nodes of the element.
		/// For a node returns itself.
		/// For an edge returns ordered pair of nodes. The order of nodes in the edge is preserved from the first creation.
		/// For a face returns ordered set of nodes. In the case Face::CheckNormalOrientation returns true the
		/// order of nodes will follow right hand side rule with respect to normal vector, otherwise it follows
		/// left hand side rule with respect to normal vector.
		/// For a cell returns the same order that was provided through suggest_nodes_oreder in Mesh::CreateCell.
		/// In the case suggest_nodes_order was not provided, the order of nodes follows VTK format for known types
		/// of elements such as Element::Tet, Element::Prism, Element::Hex, Element::Pyramid. For a general polyhedron
		/// the order is unspecified.
		/// @return array of elements
		/// @see Face::CheckNormalOrientation
		/// @see Face::FaceOrientedOutside
		virtual ElementArray<Node>  getNodes                () const; //unordered
		/// Retrive all the edges of the element.
		/// For a node returns unordered set of edges.
		/// For an edge returns itself.
		/// For a face returns ordered set of edges.
		/// For a cell returns unordered set of edges.
		/// @return array of elements
		virtual ElementArray<Edge>  getEdges                () const; //unordered
		/// Retrive all the faces of the element.
		/// For a node returns unordered set of faces.
		/// For an edge returns unordered set of faces.
		/// For a face returns itself.
		/// For a cell return ordered set of faces. The order of faces in the cell is preserved from the first creation.
		/// @return array of elements
		virtual ElementArray<Face>  getFaces                () const; //unordered
		/// Return all the cells of the element.
		/// For a node returns unordered set of cells.
		/// For an edge returns unordered set of cells.
		/// For a face returns a pair of cells. In the case Face::CheckNormalOrientation returns true
		/// then the normal points from the first cell to the second and in oppisite direction otherwise.
		/// For a cell returns itself.
		/// @return array of elements
		/// @see Face::FaceOrientedOutside
		virtual ElementArray<Cell>  getCells                () const; //unordered
		virtual ElementArray<Node>  getNodes                (MarkerType mask,bool invert_mask = false) const; //unordered
		virtual ElementArray<Edge>  getEdges                (MarkerType mask,bool invert_mask = false) const; //unordered
		virtual ElementArray<Face>  getFaces                (MarkerType mask,bool invert_mask = false) const; //unordered
		virtual ElementArray<Cell>  getCells                (MarkerType mask,bool invert_mask = false) const; //unordered
		Node                        getAsNode               () const;
		Edge                        getAsEdge               () const;
		Face                        getAsFace               () const;
		Cell                        getAsCell               () const;
		ElementSet                  getAsSet                () const;
		GeometricType               GetGeometricType        () const;
		unsigned int                GetElementDimension     () const {return GetGeometricDimension(GetGeometricType());}
		Status                      GetStatus               () const;
		void                        SetStatus               (Status status) const;
		Storage::integer &          GlobalID                () const;
		bool                        CheckElementConnectivity() const;
		void                        PrintElementConnectivity() const;
		static bool                 CheckConnectivity       (Mesh * m);
		//implemented in geometry.cpp
		void                        CastRay                 (real * pos, real * dir, dynarray< std::pair<Element, real> , 16 > & hits) const;
		void                        ComputeGeometricType    () const;
		void                        Centroid                (real * cnt) const;
		void                        Barycenter              (real * cnt) const;
		Storage::real               Mean                    (real (*func)(real* x,real t),real time) const;
		bool                        Boundary                () const;
		bool                        Planarity               () const; // check that all nodes lay on one plane
		//implemented in modify.cpp
		bool                        Hide                    () const; // if true then element was hidden, works only inside BeginModification and EndModification, on EndModification all Hidden elements are deleted
		bool                        Show                    () const; // if true then element was recovered
		bool                        Delete                  (); // if true then element was deleted, otherwise it was hidden
		bool                        Hidden                  () const;
		bool                        New                     () const;
		void                        Disconnect              (bool delete_upper_adjacent) const; //disconnect all elements, delete upper dependent
		/// Disconnects nodes from this edge, edges from this face, faces from this cell, cannot disconnect cells from this node; \n
		/// Disconnects edges from this node, faces from this edge, cells from this face, cannot disconnect nodes from this cell; \n
		/// Updates geometric data and cell nodes automatically.
		void                        Disconnect              (const HandleType * adjacent, INMOST_DATA_ENUM_TYPE num) const;
		/// Connects lower adjacencies to current element, 
		/// geometric data and cell nodes are updated automatically. \n
		/// TODO:
		///		1. asserts in this function should be replaced by Topography checks; \n
		///		2. this function should be used for creation of elements instead of current implementation. \n
		///   3. should correctly account for order of edges (may be implemented through CheckEdgeOrder, FixEdgeOrder)
		void                        Connect                 (const HandleType * adjacent, INMOST_DATA_ENUM_TYPE num) const; 
		/// Update geometric data for element, calls RecomputeGeometricData from Mesh.
		void                        UpdateGeometricData     () const; 
	};
	
	__INLINE const Element & InvalidElement() {static Element ret(NULL,InvalidHandle()); return ret;}
	
	class Node : public Element //implemented in node.cpp
	{
	public:
		Node() : Element() {}
		Node(const Node & other) : Element(other) {}
		Node(Mesh * m, HandleType h) : Element(m,h) {}
		Node(Mesh * m, HandleType * h) : Element(m,h) {}
		Node & operator =(Node const & other) {Element::operator =(other); return *this;}
		Node * operator ->() {return this;}
		const Node * operator ->() const {return this;}
		Node & self() {return *this;}
		const Node & self() const {return *this;}
	public:
		ElementArray<Edge>          getEdges                () const; //unordered
		ElementArray<Face>          getFaces                () const; //unordered
		ElementArray<Cell>          getCells                () const; //unordered

		ElementArray<Edge>          getEdges                (MarkerType mask,bool invert_mask = false) const; //unordered
		ElementArray<Face>          getFaces                (MarkerType mask,bool invert_mask = false) const; //unordered
		ElementArray<Cell>          getCells                (MarkerType mask,bool invert_mask = false) const; //unordered

		Storage::real_array         Coords                  () const; 
	};

	__INLINE const Node & InvalidNode() {static Node ret(NULL,InvalidHandle()); return ret;}
	
	class Edge : public Element //implemented in edge.cpp
	{
	public:
		Edge() : Element() {}
		Edge(const Edge & other) : Element(other) {}
		Edge(Mesh * m, HandleType h) : Element(m,h) {}
		Edge(Mesh * m, HandleType * h) : Element(m,h) {}
		Edge & operator =(Edge const & other) {Element::operator =(other); return *this;}
		Edge * operator ->() {return this;}
		const Edge * operator ->() const {return this;}
		Edge & self() {return *this;}
		const Edge & self() const {return *this;}
	public:
		
		ElementArray<Node>          getNodes                () const; //ordered
		ElementArray<Face>          getFaces                () const; //unordered
		ElementArray<Cell>          getCells                () const; //unordered

		ElementArray<Node>          getNodes                (MarkerType mask,bool invert_mask = false) const; //ordered
		ElementArray<Face>          getFaces                (MarkerType mask,bool invert_mask = false) const; //unordered
		ElementArray<Cell>          getCells                (MarkerType mask,bool invert_mask = false) const; //unordered

		Node                        getBeg                  () const;
		Node                        getEnd                  () const;
		//implemented in modify.cpp
		static Edge                 UniteEdges              (ElementArray<Edge> & edges, MarkerType del_protect);
		static bool                 TestUniteEdges          (const ElementArray<Edge> & edges, MarkerType del_protect);
		static ElementArray<Edge>   SplitEdge               (Edge e, const ElementArray<Node> & nodes, MarkerType del_protect); //provide ordered array of nodes, that lay between former nodes of the edge
		static bool                 TestSplitEdge           (Edge e, const ElementArray<Node> & nodes, MarkerType del_protect);
		//implemented in geometry.cpp
		Storage::real               Length                  () const;
	};

	__INLINE const Edge & InvalidEdge() {static Edge ret(NULL,InvalidHandle()); return ret;}
	
	class Face : public Element //implemented in face.cpp
	{
	public:
		Face() : Element() {}
		Face(const Face & other) : Element(other) {}
		Face(Mesh * m, HandleType h) : Element(m,h) {}
		Face(Mesh * m, HandleType * h) : Element(m,h) {}
		Face & operator =(Face const & other) {Element::operator =(other); return *this;}
		Face * operator ->() {return this;}
		const Face * operator ->() const {return this;}
		Face & self() {return *this;}
		const Face & self() const {return *this;}
	public:
		
		ElementArray<Node>          getNodes                () const; //ordered
		ElementArray<Edge>          getEdges                () const; //ordered
		ElementArray<Cell>          getCells                () const; //unordered

		ElementArray<Node>          getNodes                (MarkerType mask,bool invert_mask = false) const; //ordered
		ElementArray<Edge>          getEdges                (MarkerType mask,bool invert_mask = false) const; //ordered
		ElementArray<Cell>          getCells                (MarkerType mask,bool invert_mask = false) const; //unordered

		//this is for 2d case when the face is represented by segment
		Node                        getBeg                  () const;
		Node                        getEnd                  () const;
		/// Retrive the cell for which the normal points outwards.
		/// Depending on the grid construction the normal may incorrectly point inwards.
		/// You can resolve this situation by Face::FixNormalOrientation.
		/// @return the cell for which normal points outwards.
		/// @see Face::FixNormalOrientation
		Cell                        BackCell                () const;
		/// Retrive the cell for which the normal points inwards.
		/// Depending on the grid construction the normal may incorrectly point outwards.
		/// You can resolve this situation by Face::FixNormalOrientation.
		/// @return the cell for which normal points inwards.
		/// @see Face::FixNormalOrientation
		Cell                        FrontCell               () const;
		bool                        FaceOrientedOutside     (Cell c) const;
		void                        ReorderEdges            () const;
		bool                        CheckEdgeOrder          () const; //not implemented// returns true if edges of face form an ordered closed loop
		bool                        FixEdgeOrder            () const; //not implemented// returns true if edges were successfully reordered to form a closed loop
		//implemented in modify.cpp
		static Face                 UniteFaces              (ElementArray<Face> & faces, MarkerType del_protect);
		static bool                 TestUniteFaces          (const ElementArray<Face> & faces,  MarkerType del_protect);
		static ElementArray<Face>   SplitFace               (Face face, const ElementArray<Edge> & edges, MarkerType del_protect); //provide all edges that lay inside face
		static bool                 TestSplitFace           (Face face, const ElementArray<Edge> & edges, MarkerType del_protect);	
		void                        SwapCells               (); //swap back cell and front cell
		//implemented in geometry.cpp
		Storage::real               Area                    () const;
		void                        Normal                  (real * nrm) const;
		void                        UnitNormal              (real * nrm) const;
		void                        OrientedNormal          (Cell c, real * nrm) const;
		void                        OrientedUnitNormal      (Cell c, real * nrm) const;
		bool                        FixNormalOrientation    () const;  //returns true if orientation was corrected, otherwise returns false
		bool                        CheckNormalOrientation  () const; //returns true if orientation is correct, otherwise returns false
		bool                        Closure                 () const; // test integrity of polygon
	};

	__INLINE const Face & InvalidFace() {static Face ret(NULL,InvalidHandle()); return ret;}
	
	class Cell : public Element //implemented in cell.cpp
	{
	public:
		Cell() : Element() {}
		Cell(Mesh * m, HandleType h) : Element(m,h) {}
		Cell(Mesh * m, HandleType * h) : Element(m,h) {}
		Cell(const Cell & other) : Element(other) {}
		Cell & operator =(Cell const & other) {Element::operator =(other); return *this;}
		Cell * operator ->() {return this;}
		const Cell * operator ->() const {return this;}
		Cell & self() {return *this;}
		const Cell & self() const {return *this;}
	public:
		ElementArray<Node>          getNodes                () const; //ordered (for known geometric types only)
		ElementArray<Edge>          getEdges                () const; //unordered
		ElementArray<Face>          getFaces                () const; //ordered (keeps order it was created in)

		ElementArray<Node>          getNodes                (MarkerType mask,bool invert_mask = false) const; //ordered (for known geometric types only)
		ElementArray<Edge>          getEdges                (MarkerType mask,bool invert_mask = false) const; //unordered
		ElementArray<Face>          getFaces                (MarkerType mask,bool invert_mask = false) const; //ordered (keeps order it was created in)
		
		
		bool                        CheckEdgeOrder          () const; //not implemented//2D only, returns true if edges of face form an ordered closed loop
		bool                        FixEdgeOrder            () const; //not implemented//2D only, returns true if edges were successfully reordered to form a closed loop
		//implemented in modify.cpp
		static Cell                 UniteCells              (ElementArray<Cell> & cells, MarkerType del_protect);
		static bool                 TestUniteCells          (const ElementArray<Cell> & cells, MarkerType del_protect);
		static ElementArray<Cell>   SplitCell               (Cell cell, const ElementArray<Face> & faces, MarkerType del_protect); //provide all faces, that lay inside cell
		static bool                 TestSplitCell           (Cell cell, const ElementArray<Face> & faces, MarkerType del_protect);
		//implemented in geometry.cpp
		Cell                        Neighbour               (Face f) const;
		ElementArray<Cell>          NeighbouringCells       () const; // get all cells that share any face with current
		/// This function can be used to determine, if points lies inside element.
		/// Now it works only for 3-dimensional elements, at future it will be
		/// extended to support polygons and segments.
		/// @param point coordinates of the point, it is assumed that number of the
		/// coordinates is the same as the number of dimensions in the mesh.
		/// @return returns true if points inside element, false otherwise
		/// @see Mesh::GetDimensions
		bool                        Inside                  (real * point) const; //is point inside cell, check for 2d case
		real                        Volume                  () const;
		bool                        Closure                 () const; // test integrity of cell
	};

	__INLINE const Cell & InvalidCell() {static Cell ret(NULL,InvalidHandle()); return ret;}

	class ElementSet : public Element //implemented in eset.cpp
	{
	public:
		static const enumerator high_conn_reserved = 4; ///< number of reserved positions in HighConn array.
                                                    ///< The first position is the handle to parent set.
                                                    ///< The second position is handle to sibling set.
                                                    ///< The third position is handle to child set.
                                                    ///< The fourth position is number of sorted elements in the set.
                                                    ///< All the rest are positions of deleted elements.
		typedef INMOST_DATA_BULK_TYPE ComparatorType;
		static const ComparatorType UNSORTED_COMPARATOR = 0;
		static const ComparatorType GLOBALID_COMPARATOR = 1;
		static const ComparatorType CENTROID_COMPARATOR = 2;
		static const ComparatorType  IERARHY_COMPARATOR = 3;
		static const ComparatorType   HANDLE_COMPARATOR = 4;
		ElementSet() : Element() {}
		ElementSet(Mesh * m, HandleType h) : Element(m,h) {}
		ElementSet(Mesh * m, HandleType * h) : Element(m,h) {}
		ElementSet(const Cell & other) : Element(other) {}
		ElementSet & operator =(ElementSet const & other) {Element::operator =(other); return *this;}
		ElementSet * operator ->() {return this;}
		const ElementSet * operator ->() const {return this;}
		ElementSet & self() {return *this;}
		const ElementSet & self() const {return *this;}
	public:

		/// Get name of the set
		std::string                 GetName() const;
		/// Retrive parent of the set
		ElementSet                  GetParent() const;
		/// Retrive sibling set of the set, this will be next child for the parent
		ElementSet                  GetSibling() const;
		/// Retrive child set of the set
		/// Following example shows how to iterate through children of current set:
		/// for(ElementSet it = set->GetChild(); it->isValid(); it = it->GetSibling()) ...
		ElementSet                  GetChild() const;
		
		/// This will create new child for the parent
		void                        AddSibling(const ElementSet & sibling) const;
		/// Add child to current set
		void                        AddChild(const ElementSet & child) const;

		/// This will erase sibling or parent's child
		void                        RemSibling(const ElementSet & sibling) const;
		/// This will erase my child
		void                        RemChild(const ElementSet & child) const;
		/// How many there are siblings to the right of me including me
		enumerator                  CountSiblings() const;
		/// How many children I have 
		enumerator                  CountChildren() const;

		bool                        HaveSibling() const;
		bool                        HaveParent() const;
		bool                        HaveChild() const;

		/// Direct raw access to stored elements, no copy involved.
		/// The actual pointer may change when you add new elements due to reallocation of memory.
		/// Modify at your own risk.
		HandleType *                getHandles() const;
		/// Retrive number of stored handles, including invalid
		/// if you want to get number of valid elements use ElementSet::Size
		enumerator                  nbHandles() const;
		/// Retrive position after last sorted element, this does not correctly
		/// represent total number of sorted elements, since some of them may be deleted
		enumerator                  nbSorted() const;
		/// Retrive all elements by type
		enumerator                  nbAdjElements(ElementType etype) const;
		enumerator                  nbAdjElements(ElementType etype, MarkerType select, bool invert = false) const;
		/// Retrive all elements by type
		ElementArray<Element>       getAdjElements(ElementType etype) const;
		ElementArray<Element>       getAdjElements(ElementType etype, MarkerType select, bool invert = false) const;
		
		/// Retrive only nodes
		ElementArray<Node>          getNodes() const;
		ElementArray<Node>          getNodes(MarkerType select, bool invert = false) const;
		/// Retrive only edges
		ElementArray<Edge>          getEdges() const;
		ElementArray<Edge>          getEdges(MarkerType select, bool invert = false) const;
		/// Retrive only faces
		ElementArray<Face>          getFaces() const;
		ElementArray<Face>          getFaces(MarkerType select, bool invert = false) const;
		/// Retrive only cells
		ElementArray<Cell>          getCells() const;
		ElementArray<Cell>          getCells(MarkerType select, bool invert = false) const;
		/// Put one element without checking of the existance of duplicate
		/// For sorted set new element will appear at unsorted part
		void                        PutElement(HandleType e) const;
		/// Put one element without checking of the existance of duplicate
		/// For sorted set new element will appear at unsorted part
		void                        PutElement(const Storage & e) const {PutElement(e->GetHandle());}
		/// Put multiple handles without checking of the existance of duplicate
		void                        PutElements(const HandleType * handles, enumerator num) const;
		/// Put multiple handles of the other set without checking of the existance of duplicate
		void                        PutElements(const ElementSet & other) const {PutElements(other->getHandles(),other->nbHandles());}
		/// Put multiple handles without checking
		template<typename EType>
		void                        PutElements(const ElementArray<EType> & elems) const {PutElements(elems.data(),static_cast<enumerator>(elems.size()));}
		/// Put one element with checking of the existance of duplicate
		/// preserves order for sorted set, thus may be expensive
		void                        AddElement(HandleType e) const;
		/// Put one element with checking of the existance of duplicate
		/// preserves order for sorted set, thus may be expensive
		void                        AddElement(const Storage & e) const {AddElement(e->GetHandle());}
		/// Add multiple elements with checking of the existance of duplicate,
		/// preserves order for sorted set, thus may be expensive.
		/// This will also remove any duplicates in unsorted part of the set.
		/// If you inserted duplicated elements through PutElements into previously sorted array
		/// then this operation does not guarantee that those duplications will be stored.
		void                        AddElements(const HandleType * handles, enumerator num) const;
		/// Add elements of other set
		void                        AddElements(const ElementSet & other) {Unite(other);}
		/// Add multiple elements with checking of the existance of duplicate
		template<typename EType>
		void                        AddElements(const ElementArray<EType> & elems) const {AddElements(elems.data(),static_cast<enumerator>(elems.size()));}

		void                        RemoveElement(const Storage & e) const;
		void                        RemoveElements(const HandleType * handles, enumerator num) const;
		/// Remove multiple elements from the set
		template<typename EType>
		void                        RemoveElements(const ElementArray<EType> & elems) const {RemoveElements(elems.data(),static_cast<enumerator>(elems.size()));}

		/// Compute and return union with other set.
		/// Result is unordered.
		/// All elements will be unique.
		ElementArray<Element> Union(const ElementSet & other) const;
		/// Compute and return union with raw handles.
		/// Result is unordered.
		/// All elements will be unique.
		ElementArray<Element> Union(const HandleType * handles, enumerator num) const;
		/// Compute and return union with elements
		template<typename EType>
		ElementArray<Element> Union(const ElementArray<EType> & elems) const {return Union(elems.data(),static_cast<enumerator>(elems.size()));}

		ElementArray<Element> Difference(const ElementSet & other) const;
		/// Compute and return difference with raw handles.
		/// If initial set was ordered, result will preserve the order.
		/// All elements will be unique.
		ElementArray<Element> Difference(const HandleType * handles, enumerator num) const;
		/// Compute and return difference with elements
		template<typename EType>
		ElementArray<Element> Difference(const ElementArray<EType> & elems) const {return Difference(elems.data(),static_cast<enumerator>(elems.size()));}

		ElementArray<Element> Intersection(const ElementSet & other) const;
		/// Compute and return intersection with raw handles.
		/// If initial set was ordered, result will preserve the order.
		/// All elements will be unique.
		ElementArray<Element> Intersection(const HandleType * handles, enumerator num) const;
		/// Compute and return intersection with elements.
		template<typename EType>
		ElementArray<Element> Intersection(const ElementArray<EType> & elems) const {return Intersection(elems.data(),static_cast<enumerator>(elems.size()));}
		
		/// Compute and store union with raw handles.
		void Unite(const ElementSet & other) const;
		/// Compute and store union with raw handles.
		void Unite(const HandleType * handles, enumerator num) const;
		/// Compute and store union with elements.
		template<typename EType>
		void Unite(const ElementArray<EType> & elems) const {Unite(elems.data(),static_cast<enumerator>(elems.size()));}

		/// TODO
		///   If other and current sets are sorted in same way, may perform narrowing traversal by retriving
		///   mutual lower_bound/higher_bound O(log(n)) operations for detecting common subsets in sorted sets.
		///   May work good when deleting handles by small chunks, ApplyModification may greatly benefit.
		void Subtract(const ElementSet & other) const;
		/// Compute and store difference with raw handles
		void Subtract(const HandleType * handles, enumerator num) const;
		/// Compute and store difference with elements
		template<typename EType>
		void Subtract(const ElementArray<EType> & elems) const {Subtract(elems.data(),static_cast<enumerator>(elems.size()));}


		/// Compute and store intersection with raw handles
		void Intersect(const ElementSet & other) const;
		/// Compute and store intersection with raw handles
		void Intersect(const HandleType * handles, enumerator num) const;
		/// Compute and store intersection with elements
		template<typename EType>
		void Intersect(const ElementArray<EType> & elems) const {Intersect(elems.data(),static_cast<enumerator>(elems.size()));}
		/// Performs sort of the set of elements, if the set was previously sorted but
		/// have unsorted part, then unsorted part will be sorted and two parts will be merged.
		/// If you need all the set to be resorted (for example in the case global ids were changed)
		/// then invoke SortSet with UNSORTED_COMPARATOR first and then with needed comparator.
		/// Internally it uses:
		/// std::sort with Mesh::CentroidComparator for CENTROID_COMPARATOR
		/// std::sort with Mesh::IerarhyComparator  for  IERARHY_COMPARATOR
		/// radix sort starting from certain size for   GLOBALID_COMPARATOR
		/// radix sort starting from certain size for     HANDLE_COMPARATOR
		/// only changes state, no sort performed for   UNSORTED_COMPARATOR 
		/// After the set was sorted all the invalid handles should appear at the end of the set
		/// and then removed from array, so it will implicitly work as ReorderEmpty function.
		/// No checks that elements are hidden performed (Maybe this checks should be done in comparators)
		/// In the case you formed the set by running over all mesh elements from NODE to ESET in 
		/// increasing order then your set will be automatically sorted by handles, in this case
		/// you are encouraged to override Mesh::SetComparatorTag with HANDLE_COMPARATOR on the set
		/// without invoking SortSet, so that SortSet does not redundantly perform any work.
		/// You are encouraged to do so even if you are not going to use this information -
		/// some internal algorithms may benefit from it.
		/// !TODO 52 - check radix sort on big endian computer
		/// @param comp one of the comparators from description
		/// @see Mesh::SetComparatorTag
		void SortSet(ComparatorType comp) const;
		/// Perform binary search by global id in set sorted with GLOBALID_COMPARATOR in O(log(n)) time
		/// otherwise search needs O(n) comparisons
		Element FindElementByGlobalID(integer global_id) const;
		/// Perform binary search by centroid in set sorted with CENTROID_COMPARATOR in O(log(n)) time
		/// otherwise search needs O(n) comparisons
		Element FindElementByCentroid(real * centroid) const;
		/// Performs linear search in unsorted set 
		/// binary search by handle in set sorted with      HANDLE_COMPARATOR
		/// binary search by centroid in set sorted with  CENTROID_COMPARATOR
		/// binary search by global id in set sorted with GLOBALID_COMPARATOR
		/// binary search by ierarhy in set sorted with    IERARHY_COMPARATOR (may be too expensive)
		/// If you have a lot of elements to test against set, then you'd better first put marker
		/// on elements of the set by SetMarkerElements then test which of your elements have marker
		/// and then remove markers by RemMarkerElements. This will consume only O(n+m) operations
		/// where n is the number of my elements and m is the size of the test set.
		/// @param h handle of element
		/// @param use_comparator use information about current comparator or perform linear search
		/// @return true if element is present
		/// @see ElementSet::SetMarkerElements
		/// @see ElementSet::RemMarkerElements
		bool FindHandle(HandleType h, bool use_comparator) const;
		/// Set markers on all the elements of given type
		void SetMarkerElements(MarkerType m, ElementType etype = ESET|CELL|FACE|EDGE|NODE) const;
		/// Remove markers from all the elements of given type
		void RemMarkerElements(MarkerType m, ElementType etype = ESET|CELL|FACE|EDGE|NODE) const;
		class iterator
		{
			Mesh * m;
			Element::adj_type const * ptr;
			Element::adj_type::size_type pos;
		public:
			typedef std::forward_iterator_tag iterator_category;
			iterator() : m(NULL), ptr(NULL), pos(0) {}
			iterator(const iterator & other) : m(other.m), ptr(other.ptr), pos(other.pos) {}
			iterator(Mesh * m, Element::adj_type const * ptr, Element::adj_type::size_type pos) : m(m), ptr(ptr), pos(pos) {}
			iterator & operator = (iterator const & other) {m = other.m; ptr = other.ptr; pos = other.pos; return *this;}
			iterator & operator ++();
			iterator & operator ++(int) {iterator ret(*this); operator++(); return *this;}
			bool       operator ==(const iterator & other) const {assert(ptr == other.ptr); return pos == other.pos;}
			bool       operator !=(const iterator & other) const {assert(ptr == other.ptr); return pos != other.pos;}
			bool       operator < (const iterator & other) const {assert(ptr == other.ptr); return pos < other.pos;}
			bool       operator > (const iterator & other) const {assert(ptr == other.ptr); return pos > other.pos;}
			bool       operator <=(const iterator & other) const {assert(ptr == other.ptr); return pos <= other.pos;}
			bool       operator >=(const iterator & other) const {assert(ptr == other.ptr); return pos >= other.pos;}
			const HandleType & operator *() const {return ptr->at(pos);}
			Element operator->() const {return Element(m,ptr->at(pos));}
			friend class ElementSet;
		};
		/// Provides forward iterator that skips deleted and hidden elements within set.
		/// The iterator will be valid during removal or insertion of elements.
		/// While that may seem usefull, many functions like AddElements, SortSet, Union will internally
		/// reorder handles thus making iterators useless.
		/// To correctly resolve situation when size of array shrink due to removal of elements
		/// use less then (<) instead of not equal (!=) operator to check for end of the loop.
		/// Note that iterators would not let you change underlaying handles, you can use getHandles
		/// for that instead.
		/// @return iterator that points to the first valid element
		iterator Begin() const;
		/// Provides end for forward iterator to stop the loop.
		/// @return iterator that points to the position after last element
		iterator End() const;
		/// Provides iterator that points to element located after the last element that belong to
		/// presorted part of the set as well as on the first element of unsorted part
		/// @return see description
		iterator EndSorted() const;
		/// Erase one element pointed by iterator and return next valid element
		iterator Erase(iterator pos) const;
		/// Erase set of elements pointed by iterators
		void Erase(iterator beg, iterator end) const;
		/// Retrive current set comparator
		ComparatorType GetComparator() const;
		/// Compact holes in inner representation
		void ReorderEmpty() const;
		/// Is there any elements in the set
		bool Empty() const;
		/// Get total number of elements
		enumerator Size() const;
		/// Remove all elements, clear all data, removes sorted marker
		void Clear() const;
	};

	__INLINE const ElementSet & InvalidElementSet() {static ElementSet ret(NULL,InvalidHandle()); return ret;}
	
	class Mesh : public TagManager, public Storage //implemented in mesh.cpp
	{
	public:
		enum MeshState {Serial, Parallel};
		typedef chunk_array<integer,chunk_bits_empty>               empty_container;
		typedef chunk_array<integer,chunk_bits_elems>               links_container;
		//typedef std::vector<integer>                                empty_container;
		//typedef std::vector<integer>                                links_container;
		typedef TagManager::sparse_sub_type                         sparse_type;
		typedef TagManager::sparse_sub_record                       sparse_rec;
		typedef sparse_type::size_type                              senum;
	private:
		Storage::real                               epsilon;
		empty_container                             empty_space[6];
		empty_container                             empty_links[6];
		links_container                             links[6];
		Tag                                         tag_global_id;
		Tag                                         tag_coords;
		Tag                                         tag_low_conn;
		Tag                                         tag_high_conn;
		Tag                                         tag_markers;
		Tag                                         tag_geom_type;
		Tag                                         tag_setname;
		Tag                                         tag_setcomparator;
		MeshState                                   m_state;
		integer                                     dim;
		HandleType                                  last_created;
	private:
		ElementType                                 have_global_id;
		INMOST_DATA_BIG_ENUM_TYPE                   parallel_mesh_unique_id;
		INMOST_MPI_Comm                             comm;
		Tag                                         tag_shared;
		Tag                                         tag_owner;
		Tag                                         tag_processors;
		Tag                                         tag_layers;
		Tag                                         tag_sendto;
		Tag                                         tag_bridge;
		Tag                                         tag_redistribute;
	private:
		__INLINE static sparse_rec        mkrec              (const Tag & t) {sparse_rec ret; ret.tag = t.mem; ret.rec = NULL; return ret;}
		__INLINE sparse_type const &      MGetSparseLink     (integer etypenum, integer ID) const {return GetSparseData(etypenum,links[etypenum][ID]);}
		__INLINE sparse_type &            MGetSparseLink     (integer etypenum, integer ID) {return GetSparseData(etypenum,links[etypenum][ID]);}
		__INLINE sparse_type const &      MGetSparseLink     (HandleType h) const {return MGetSparseLink(GetHandleElementNum(h),GetHandleID(h));}
		__INLINE sparse_type &            MGetSparseLink     (HandleType h) {return MGetSparseLink(GetHandleElementNum(h),GetHandleID(h));}
		__INLINE const void *             MGetSparseLink     (HandleType h, const Tag & t) const {sparse_type const & s = MGetSparseLink(GetHandleElementNum(h),GetHandleID(h)); for(senum i = 0; i < s.size(); ++i) if( s[i].tag == t.mem ) return s[i].rec; return NULL;}
		__INLINE void * &                 MGetSparseLink     (HandleType h, const Tag & t) {sparse_type & s = MGetSparseLink(GetHandleElementNum(h),GetHandleID(h)); for(senum i = 0; i < s.size(); ++i) if( s[i].tag == t.mem ) return s[i].rec; s.push_back(mkrec(t)); return s.back().rec;}
		__INLINE const void *             MGetDenseLink      (integer n, integer id, const Tag & t) const {return &(GetDenseData(t.GetPositionByDim(n))[links[n][id]]);}
		__INLINE void *                   MGetDenseLink      (integer n, integer id, const Tag & t) {return &(GetDenseData(t.GetPositionByDim(n))[links[n][id]]);}
		__INLINE const void *             MGetDenseLink      (HandleType h, const Tag & t) const {return MGetDenseLink(GetHandleElementNum(h),GetHandleID(h),t);}
		__INLINE void *                   MGetDenseLink      (HandleType h, const Tag & t) {return MGetDenseLink(GetHandleElementNum(h),GetHandleID(h),t);}
		__INLINE const void *             MGetLink           (HandleType h, const Tag & t) const {if( !t.isSparseByDim(GetHandleElementNum(h)) ) return MGetDenseLink(h,t); else return MGetSparseLink(h,t);}
		__INLINE void *                   MGetLink           (HandleType h, const Tag & t) {if( !t.isSparseByDim(GetHandleElementNum(h)) ) return MGetDenseLink(h,t); else {void * & q = MGetSparseLink(h,t); if( q == NULL ) q = calloc(1,t.GetRecordSize()); return q;}}
	public:
		/// For debug purposes
		integer                           HandleDataPos      (HandleType h) {return links[GetHandleElementNum(h)][GetHandleID(h)];}
		/// For parmetis
		/// return total number in bytes of occupied memory by element and its data
		enumerator                        MemoryUsage        (HandleType h);
		Mesh();
		Mesh(const Mesh & other);
		Mesh & operator =(Mesh const & other);
		~Mesh();
		/// Allocate new marker.
		/// Assert will fire in debug mode (NDEBUG not set) if you run out of space for markers, in this case you either
		/// use too many markers, then you can just increase MarkerFields variable (increasing by 1 gives you 8 more markers)
		/// or you forget to release markers after you use them.
		/// In release mode (NDEBUG is set) if you run out of space for markers function will return InvalidMarker()
		/// @return new marker or InvalidMarker(), see description
		MarkerType                        CreateMarker       ();
		/// Release marker back for reuse.
		/// This function will only notice mesh that the marker is free to be reused, this will not clear markers
		/// from elements. Before releasing the marker you should ensure that all the marker is removed from all the elements.
		/// Since many inner algorithms use markers, not properly clearing markers from elements may lead to undefined behavior.
		/// Since it is expensive to check asserts will fire in debug mode (NDEBUG not set) only if you define CHECKS_MARKERS in inmost_common.h,
		/// no checks will be performed in release mode(NDEBUG is set).
		/// @param n byte position and byte bit mask
		void                              ReleaseMarker      (MarkerType n);
		/// Set tolerance for coordinates comparison, this tolerance is used in comparators 
		/// when two meshes are merged during loading, in ResolveShared to check that nodes on different processors match 
		/// and in UnpackElementsData
		/// @param e small real value
		__INLINE void                     SetEpsilon         (real e) {epsilon = e;}
		/// Retrive tolerance for coordinates comparison
		/// @return real value
		__INLINE real                     GetEpsilon         () const {return epsilon;}
		/// Set number of dimensions for mesh, it sets the size for number of internally stored coordinates.
		/// Size of the array returned by Node::Coords will match this number.
		/// When you change number of dimensions and there are nodes with bigger number of dimensions then
		/// first dim coordinates are copied and the rest are dropped.
		/// @see Node::Coords
		void                              SetDimensions      (integer dim);
		/// Get number of dimensions of mesh. 
		/// Size of the array returned by Node::Coords will match this number.
		/// @see Node::Coords
		__INLINE integer                  GetDimensions      () const {return dim;}
		/// Get parallel state of the mesh
		/// @return either Mesh::Serial or Mesh::Parallel
		__INLINE MeshState                GetMeshState       () const {return m_state;}
		__INLINE const Tag &              GlobalIDTag        () const {return tag_global_id;}
		__INLINE const Tag &              CoordsTag          () const {return tag_coords;}
		__INLINE const Tag &              LowConnTag         () const {return tag_low_conn;}
		__INLINE const Tag &              HighConnTag        () const {return tag_high_conn;}
		__INLINE const Tag &              MarkersTag         () const {return tag_markers;}
		__INLINE const Tag &              GeomTypeTag        () const {return tag_geom_type;}
		__INLINE const Tag &              SendtoTag          () const {return tag_sendto;}
		__INLINE const Tag &              SharedTag          () const {return tag_shared;}
		__INLINE const Tag &              OwnerTag           () const {return tag_owner;}
		__INLINE const Tag &              LayersTag          () const {return tag_layers;}
		__INLINE const Tag &              BridgeTag          () const {return tag_bridge;}
		__INLINE const Tag &              ProcessorsTag      () const {return tag_processors;}
		__INLINE const Tag &              SetNameTag         () const {return tag_setname;}
		__INLINE const Tag &              SetComparatorTag   () const {return tag_setcomparator;}
		/// Don't put this shortcut to any function directly, as it creates tag inside
		/// assign to other object of type Tag and put this object to functions.
		__INLINE Tag                      RedistributeTag    () {return CreateTag("TEMPORARY_NEW_OWNER",DATA_INTEGER,CELL,NONE,1);}
		/// Create the tag by name, type and size.
		/// You cannot create the tag with the same name and different type or size.
		/// You may make subsequent calls to the function with the same name, same type and same size 
		/// (or undefined size, then it will be deduced) but different selection of element.
		/// The following recomendation is for future:
		///   When you create data of variable size, every array of data for every element will be empty at first, so
		///   before accessing it through mechanism for single-valued data (Real, Integer, Bulk, Reference) 
		///   you should first resize arrays otherwise your program very likely to be halted by segmentation fault. 
		///   For this case arrays are not allocated automatically from performance considerations.
		/// @param name name of the tag
		/// @param dtype type of the tag
		/// @param etype the selection of elements on which the data of the tag is defined, you may use bitwise or operations 
		/// to define tag on multiple types of elements, example CELL | FACE
		/// @param sparse the selection of elements from etype on which the tag is sparse, for example, if you know that the data is used 
		/// on all cells and only on boundary faces, then you may should set etype = CELL | FACE and sparse = FACE
		/// @return returns the tag that represents the data
		Tag                               CreateTag          (std::string name, DataType dtype, ElementType etype,ElementType sparse, INMOST_DATA_ENUM_TYPE size = ENUMUNDEF);
		/// Remove the data that is represented by the tag from elements of selected type
		/// @param tag tag that indicates the data
		/// @param mask the selection of the elements on which the data should be removed, may be set by bitwise or operation
		/// @return returns tag that returns false on isValid() request if all the data was removed and the tag is no more occupied, otherwise returns the tag
		/// @see Tag::isValid
		Tag                               DeleteTag          (Tag tag, ElementType mask = NODE | EDGE | FACE | CELL | ESET | MESH);
		/// Returns the first tag defined on the mesh
		/// For safety while iterating through tags you should check for validity of the tag
		/// @return first tag
		__INLINE iteratorTag              BeginTag           () {return tags.begin(); }
		/// Retrive the total number of valid tags
		/// @return returns the number of valid tags
		__INLINE enumerator               NumberOfTags       () const { return static_cast<INMOST_DATA_ENUM_TYPE>(tags.size()); }
		/// Returns the indicator for loop to end iteration over tags
		/// @return the inexistant tag that is located the one position after the last tag
		__INLINE iteratorTag              EndTag             () {return tags.end(); }
		/// Create node by given coordinates.
		/// This operation would not involve searching existing nodes for node with the same coordinates.
		/// It is potentially dangerous to have nodes whose coordinates differ by less then GetEpsilon since ResolveShared
		/// algorithm would not know how to resolve parallel statuses of the elements. However this may be subject to change,
		/// if ResolveShared algorithm will be rewritten to resolve cells first by their centroids and only then resolve all the rest
		/// elements by adjacency information.
		/// @param coords array of coordinates at least of size GetDimensions()
		/// @return interface to created node
		Node                              CreateNode         (const real * coords);
		std::pair<Edge,bool>              CreateEdge         (const ElementArray<Node> & nodes);
		std::pair<Face,bool>              CreateFace         (const ElementArray<Edge> & edges);
		std::pair<Face,bool>              CreateFace         (const ElementArray<Node> & nodes);
		std::pair<Cell,bool>              CreateCell         (const ElementArray<Face> & faces, const ElementArray<Node> & suggest_nodes_order = ElementArray<Node>(NULL));
		std::pair<Cell,bool>              CreateCell         (const ElementArray<Node> & c_f_nodes, const integer * c_f_numnodes, integer num_c_faces, 
		                                                      const ElementArray<Node> & suggest_nodes_order = ElementArray<Node>(NULL));
		std::pair<Cell,bool>              CreateCell         (const ElementArray<Node> & c_nodes, const integer * c_f_nodeinds, const integer * c_f_numnodes, integer num_c_faces, 
		                                                      const ElementArray<Node> & suggest_nodes_order = ElementArray<Node>(NULL));
		std::pair<ElementSet,bool>        CreateSet          (std::string name);
		/// Retrive set by name.
		/// @param name set name
		/// @return set whose name match or InvalidHandle()
		ElementSet                        GetSet             (std::string name);
		/// Retrive all the sets whose names start with given prefix
		ElementArray<ElementSet>          GetSetsByPrefix    (std::string prefix);
		HandleType                        LastCreated        () const {return last_created;}

		bool                              isValidHandleRange (HandleType h) const; //for asserts
		bool                              isValidElement     (integer etypenum, integer lid) const {return links[etypenum][lid] != -1;}
		bool                              isValidElement     (ElementType etype, integer lid) const {return links[ElementNum(etype)][lid] != -1;}
		bool                              isValidElement     (HandleType h) const {return isValidHandle(h) && isValidElement(GetHandleElementNum(h),GetHandleID(h));}
		/// Retrive upper adjacent that is shared by multiple lower adjacencies
		/// @return handle of found element or InvalidHandle()
		HandleType                        FindSharedAdjacency(const HandleType * arr, enumerator num) const;
		void                              ReorderEmpty       (ElementType reordertypes);
		//Bug inside: sort would not work on chunk_array, because it is not contiguous in memory
		void                              ReorderApply       (Tag index, ElementType mask);
		void                              RestoreCellNodes   (HandleType hc, ElementArray<Node> & ret);
	private:
		//those functions contain asserts for debug purposes, in release mode (NDEBUG is set) they are empty and call should be optimized out in worst case by linker
		void                              Asserts            (HandleType h, const Tag & tag, DataType expected) const;
		void                              AssertsDF          (HandleType h, const Tag & tag, DataType expected) const;
		void                              AssertsDV          (HandleType h, const Tag & tag, DataType expected) const;
	public:
		/// Returns a reference to inner memory location of the first element of the array of real values.
		/// Future recomendation:
		///   If variable size array was not allocated then this function will generate segmentation fault.
		/// If you know that data is certanly dense and fixed or variable on elements you access then it is
		/// faster to use specialized variants of this function.
		/// Reference to the data is guaranteed to be valid during mesh modification.
		/// @param h element handle
		/// @param tag tag that represents data
		/// @see Mesh::RealDF
		/// @see Mesh::RealDV
		real      &                       Real               (HandleType h, const Tag & tag);
		/// Returns a reference to inner memory location of the first element of the array of integer values.
		/// Future recomendation:
		///   If variable size array was not allocated then this function will generate segmentation fault.
		/// If you know that data is certanly dense and fixed or variable on elements you access then it is
		/// faster to use specialized variants of this function.
		/// Reference to the data is guaranteed to be valid during mesh modification.
		/// @param h element handle
		/// @param tag tag that represents data
		integer   &                       Integer            (HandleType h, const Tag & tag);
		/// Returns a reference in inner representation to the first element of array of bytes.
		/// Future recomendation:
		///   If variable size array was not allocated then this function will generate segmentation fault.
		/// If you know that data is certanly dense and fixed or variable on elements you access then it is
		/// faster to use specialized variants of this function.
		/// Reference to the data is guaranteed to be valid during mesh modification.
		/// @param h element handle
		/// @param tag tag that represents data
		bulk      &                       Bulk               (HandleType h, const Tag & tag);
		/// Returns a reference in inner representation to the first element of array of element handles.
		/// Future recomendation:
		///   If variable size array was not allocated then this function will generate segmentation fault.
		/// If you know that data is certanly dense and fixed or variable on elements you access then it is
		/// faster to use specialized variants of this function.
		/// Reference to the data is guaranteed to be valid during mesh modification.
		/// Using handle you can construct objects of type Storage, Element, Node, Edge,
		/// Face, Cell, ElementSet by calling their constructor with pointer to mesh and handle
		/// as arguments.
		/// @param h element handle
		/// @param tag tag that represents data
		reference &                       Reference          (HandleType h, const Tag & tag);
		/// Returns an array of real values.
		/// If you know that data is certanly dense on elements you access then it is faster to use
		/// variants of this function with hint data structure.
		/// Array data structure is guaranteed to be valid during mesh modification.
		/// @param h element handle
		/// @param tag tag that represents data
		real_array                        RealArray          (HandleType h, const Tag & tag);
		/// Returns an array of integer values.
		/// If you know that data is certanly dense and fixed or variable on elements you access then it is
		/// faster to use specialized variants of this function.
		/// variants of this function with hint data structure.
		/// Array data structure is guaranteed to be valid during mesh modification.
		/// @param h element handle
		/// @param tag tag that represents data
		integer_array                     IntegerArray       (HandleType h, const Tag & tag);
		/// Returns an array of bytes.
		/// If variable size array was not allocated then this function will generate segmentation fault.
		/// If you know that data is certanly sparse or dense on elements you access then it is faster to use
		/// variants of this function with hint data structure.
		/// Array data structure is guaranteed to be valid during mesh modification.
		/// @param h element handle
		/// @param tag tag that represents data
		bulk_array                        BulkArray          (HandleType h, const Tag & tag);
		/// Returns an array of element handles.
		/// If variable size array was not allocated then this function will generate segmentation fault.
		/// If you know that data is certanly sparse or dense on elements you access then it is faster to use
		/// variants of this function with hint data structure.
		/// The class reference_array that is used to represent array of elements stores handles inside but
		/// accessing them through square scopes [] or by arrow -> in iterator will automatically
		/// form object of type Element. If you are not sure that stored handles are valid, you
		/// should either check that by Element::isValid (involves deep check) or test handle
		/// against InvalidHandle (simple check). To obtain handle you may use reference_array::at function or
		/// dereference * operator for iterator. If you need custom object like Node, Edge, Face, Cell
		/// or ElementSet you may use Element::getAsNode, Element::getAsEdge, Element::getAsFace,
		/// Element::getAsCell and Element::getAsSet functions.
		/// Array data structure is guaranteed to be valid during mesh modification. If you delete
		/// elements by Mesh::Delete or Element::Delete all the references are also will be valid
		/// and reverted to InvalidHandle on Mesh::ApplyModification. If you use Mesh::Destroy to
		/// delete mesh elements or delete elements not within modification state then references 
		/// may become either invalid but not testable against InvalidHandle (situation may be tested 
		/// by Element::isValid or Mesh::isValidHandle) or reference may be reused by another element 
		/// if you mix deletion and creation of elements and then is no way to resolve this situation,
		/// except if you have created only one element the it may be retrived by Mesh::LastHandle.
		/// @param h element handle
		/// @param tag tag that represents data
		/// @see InvalidHandle
		/// @see Element::isValid
		/// @see Element::getAsNode
		/// @see Element::getAsEdge
		/// @see Element::getAsFace
		/// @see Element::getAsCell
		/// @see Element::getAsSet
		reference_array                   ReferenceArray     (HandleType h, const Tag & tag);
		/// Returns a reference to inner memory location of the first element of the array of real values.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// @param h element handle
		/// @param tag tag that represents data
		real      &                       RealDF             (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_REAL     ); return static_cast<real     *>(MGetDenseLink(h,tag))[0];}
		/// Returns a reference to inner memory location of the first element of the array of integer values.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		integer   &                       IntegerDF          (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_INTEGER  ); return static_cast<integer  *>(MGetDenseLink(h,tag))[0];}
		/// Returns a reference in dense array to the first element of constant size array of bytes.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		bulk      &                       BulkDF             (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_BULK     ); return static_cast<bulk     *>(MGetDenseLink(h,tag))[0];}
		/// Returns a reference in dense array to the first element of constant size array of element handles.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// Using handle you can construct objects of type Storage, Element, Node, Edge,
		/// Face, Cell, ElementSet by calling their constructor with pointer to mesh and handle
		/// as arguments.
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		reference &                       ReferenceDF        (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_REFERENCE); return static_cast<reference*>(MGetDenseLink(h,tag))[0];}
		/// Returns an array of real values in dense array.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// Note that as array is fixed you shouldn't use any functions that alter size of the array
		/// as resize, erase, insert, you may use replace if initial and final size will match,
		/// in debug mode assert will fire if you try to do this in release (NDEBUG is set) it will lead to segfault.
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		real_array                        RealArrayDF        (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_REAL     ); return real_array     (static_cast<real     *>(MGetDenseLink(h,tag)),tag.GetSize());}
		/// Returns an array of integer values in dense array.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set), likely to result in segfault.
		/// Note that as array is fixed you shouldn't use any functions that alter size of the array
		/// as resize, erase, insert, you may use replace if initial and final size will match,
		/// in debug mode assert will fire if you try to do this in release (NDEBUG is set) it will lead to segfault.
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		integer_array                     IntegerArrayDF     (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_INTEGER  ); return integer_array  (static_cast<integer  *>(MGetDenseLink(h,tag)),tag.GetSize());}
		/// Returns an array of bytes in dense array.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set), likely to result in segfault.
		/// Note that as array is fixed you shouldn't use any functions that alter size of the array
		/// as resize, erase, insert, you may use replace if initial and final size will match,
		/// in debug mode assert will fire if you try to do this in release (NDEBUG is set) it will lead to segfault.
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		bulk_array                        BulkArrayDF        (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_BULK     ); return bulk_array     (static_cast<bulk     *>(MGetDenseLink(h,tag)),tag.GetSize());}
		/// Returns an array of element handles in dense array.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set), likely to result in segfault.
		/// Note that as array is fixed you shouldn't use any functions that alter size of the array
		/// as resize, erase, insert, you may use replace if initial and final size will match,
		/// in debug mode assert will fire if you try to do this in release (NDEBUG is set) it will lead to segfault.
		/// The class reference_array that is used to represent array of elements stores handles inside but
		/// accessing them through square scopes [] or by arrow -> in iterator will automatically
		/// form object of type Element. If you are not sure that stored handles are valid, you
		/// should either check that by Element::isValid (involves deep check) or test handle
		/// against InvalidHandle (simple check). To obtain handle you may use reference_array::at function or
		/// dereference * operator for iterator. If you need custom object like Node, Edge, Face, Cell
		/// or ElementSet you may use Element::getAsNode, Element::getAsEdge, Element::getAsFace,
		/// Element::getAsCell and Element::getAsSet functions.
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		/// @see InvalidHandle
		/// @see Element::isValid
		/// @see Element::getAsNode
		/// @see Element::getAsEdge
		/// @see Element::getAsFace
		/// @see Element::getAsCell
		/// @see Element::getAsSet
		reference_array                   ReferenceArrayDF   (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_REFERENCE); return reference_array(this,static_cast<reference*>(MGetDenseLink(h,tag)),tag.GetSize());}
		/// Returns a reference in dense array to the first element of variable size array of real values.
		/// Future recomendation:
		///    If array was not allocated (resized) then this function will generate segmentation fault.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// @param h element handle
		/// @param tag tag that represents data
		real      &                       RealDV             (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_REAL     ); return static_cast<inner_real_array     *>(MGetDenseLink(h,tag))->at_safe(0);}
		/// Returns a reference in dense array to the first element of variable size array of integer values.
		/// Future recomendation:
		///    If array was not allocated then this function will generate segmentation fault.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// @param h element handle
		/// @param tag tag that represents data
		integer   &                       IntegerDV          (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_INTEGER  ); return static_cast<inner_integer_array  *>(MGetDenseLink(h,tag))->at_safe(0);}
		/// Returns a reference in dense array to the first element of variable size array of bytes.
		/// Future recomendation:
		///    If array was not allocated then this function will generate segmentation fault.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// @param h element handle
		/// @param tag tag that represents data
		/// @see TagDenseVariable
		bulk      &                       BulkDV             (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_BULK     ); return static_cast<inner_bulk_array     *>(MGetDenseLink(h,tag))->at_safe(0);}
		/// Returns a reference in dense array to the first element of variable size array of element handles.
		/// Future recomendation:
		///    If array was not allocated then this function will generate segmentation fault.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// Using handle you can construct objects of type Storage, Element, Node, Edge,
		/// Face, Cell, ElementSet by calling their constructor with pointer to mesh and handle
		/// as arguments.
		/// @param h element handle
		/// @param tag tag that represents data
		reference &                       ReferenceDV        (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_REFERENCE); return static_cast<inner_reference_array*>(MGetDenseLink(h,tag))->at_safe(0);}
		/// Returns an array of real values in dense array of variable size.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// @param h element handle
		/// @param tag tag that represents data
		real_array                        RealArrayDV        (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_REAL     ); return real_array     (*static_cast<inner_real_array     *>(MGetDenseLink(h,tag)));}
		/// Returns an array of integer values in dense array of variable size.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// @param h element handle
		/// @param tag tag that represents data
		integer_array                     IntegerArrayDV     (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_INTEGER  ); return integer_array  (*static_cast<inner_integer_array  *>(MGetDenseLink(h,tag)));}
		/// Returns an array of bytes in dense array of variable size.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// @param h element handle
		/// @param tag tag that represents data
		bulk_array                        BulkArrayDV        (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_BULK     ); return bulk_array     (*static_cast<inner_bulk_array     *>(MGetDenseLink(h,tag)));}
		/// Returns an array of element handles in dense array of variable size.
		/// If you don't know any hint information about tag data you should not use this function.
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		/// The class reference_array that is used to represent array of elements stores handles inside but
		/// accessing them through square scopes [] or by arrow -> in iterator will automatically
		/// form object of type Element. If you are not sure that stored handles are valid, you
		/// should either check that by Element::isValid (involves deep check) or test handle
		/// against InvalidHandle (simple check). To obtain handle you may use reference_array::at function or
		/// dereference * operator for iterator. If you need custom object like Node, Edge, Face, Cell
		/// or ElementSet you may use Element::getAsNode, Element::getAsEdge, Element::getAsFace,
		/// Element::getAsCell and Element::getAsSet functions.
		/// @param h element handle
		/// @param tag tag that represents data
		/// @see InvalidHandle
		/// @see Element::isValid
		/// @see Element::getAsNode
		/// @see Element::getAsEdge
		/// @see Element::getAsFace
		/// @see Element::getAsCell
		/// @see Element::getAsSet
		reference_array                   ReferenceArrayDV   (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_REFERENCE); return reference_array(this,*static_cast<inner_reference_array*>(MGetDenseLink(h,tag)));}
		/// Set a marker on the element represented by handle
		/// @param h element handle
		/// @param mask stores byte number and byte bit mask that represent marker
		void                              SetMarker          (HandleType h,MarkerType n) {static_cast<bulk *>(MGetDenseLink(h,MarkersTag()))[n >> MarkerShift] |= static_cast<bulk>(n & MarkerMask);}
		/// Set a marker on the set of handles
		/// @param h set of handles
		/// @param n number of handles
		/// @param m stores byte number and byte bit mask that represent marker
		/// @see Mesh::SetMarker
		void                              SetMarkerArray     (const HandleType * h, enumerator n, MarkerType m) {for(enumerator i = 0; i < n; ++i) if( h[i] != InvalidHandle() )SetMarker(h[i],m);}
		/// Check weather the marker is set one the element
		/// @param h element handle
		/// @param mask stores byte number and byte bit mask that represent marker
		bool                              GetMarker          (HandleType h,MarkerType n) const {return (static_cast<const bulk *>(MGetDenseLink(h,MarkersTag()))[n >> MarkerShift] & static_cast<bulk>(n & MarkerMask)) != 0;}
		/// Remove the marker from the element
		/// @param h element handle
		/// @param mask stores byte number and byte bit mask that represent marker
		void                              RemMarker          (HandleType h,MarkerType n) {static_cast<bulk *>(MGetDenseLink(h,MarkersTag()))[n >> MarkerShift] &= ~static_cast<bulk>(n & MarkerMask);}
		/// Remove the marker from the set of handles
		/// @param h set of handles
		/// @param n number of handles
		/// @param m stores byte number and byte bit mask that represent marker
		/// @see Mesh::RemMarker
		void                              RemMarkerArray     (const HandleType * h, enumerator n, MarkerType m) {for(enumerator i = 0; i < n; ++i) if( h[i] != InvalidHandle() ) RemMarker(h[i],m);}
		/// Remove all the markers from the element
		void                              ClearMarkerSpace   (HandleType h);
		/// Get a copy of the bytes that store markers on the element
		/// @param h element handle
		/// @param copy storage to put data to
		void                              GetMarkerSpace     (HandleType h,bulk copy  [MarkerFields]) const;
		/// Overwrite the bytes that store markers on the element
		/// @param h element handle
		/// @param source storage to get data from
		void                              SetMarkerSpace     (HandleType h,bulk source[MarkerFields]);
		/// Access directly higher order adjacencies of current element with right of modification.
		/// This function gives direct access to elements and allows you to overwrite handles
		/// which is not recomended. You bypass topology checks, correct connectivity estableshment and geometric data updates.
		/// If you do so then check connectivity afterwards by Element::CheckElementConnectivity for debugging purposes.
		/// Then for correct function of geometrical algorithms in order from lower modified adjacenices to upper modified adjacencies
		/// (if some element have lower adjacencies updated then it should be updated otherwise it shouldn't)
		/// call ComputeGeometricType to deduce geometric representation or force it by SetGeometricType,
		/// then for element of type CELL call RestoreCellNodes, clear current mutual connection to the nodes and 
		/// establish new mutual connections from nodes returned by RestoreCellNodes;
		/// then for all elements call RecomputeGeometricData to automatically recompute all geometric quantities.
		/// Don't forget that edges of the face should be ordered for correct retrival of nodes, otherwise
		/// Face::getNodes and Element::getNodes for FACE in 3 dimensions or
		/// Cell::getNodes and Element::getNodes for CELL in 2 dimensions will return garbage
		///
		/// For NODE this returns edges that are connected to this node;
		/// For EDGE this returns faces that are connected to this edge;
		/// For FACE this returns cells that are connected to this face
		/// For CELL this returns nodes of the cell
		/// For ESET first three entries are parent, sibling, child then records represent empty positions in LowConn
		/// @param h handle of the element
		/// @return see description
		/// @see Element::CheckElementConnectivity
		/// @see Element::Connect
		/// @see Element::Disconnect
		/// @see Mesh::SetGeometricType
		/// @see Mesh::ComputeGeometricType
		/// @see Mesh::RestoreCellNodes
		/// @see Mesh::RecomputeGeometricData
		/// @see Face::FixNormalOrientation
		/// @see Mesh::LowConn
		Element::adj_type &               HighConn           (HandleType h) {return *static_cast<inner_reference_array*>(MGetDenseLink(h,HighConnTag()));}
		/// Access directly higher order adjacencies of current element without right of modification.
		Element::adj_type const&          HighConn           (HandleType h) const {return *static_cast<const inner_reference_array*>(MGetDenseLink(h,HighConnTag()));}
		/// Access directly lower order adjacencies of current element with right of modification.
		/// This function gives direct access to elements and allows you to overwrite handles.
		/// If you are going to overwrite them then read recomendations in description for HighConn function.
		/// For NODE this returns cells that are connected to this node;
		/// For EDGE this returns nodes of the edge
		/// For FACE this returns edges of the face
		/// For CELL this returns faces of the cell
		/// For ESET handles of the elements that this set contain
		/// @param h handle of the element
		/// @return see description
		/// @see Mesh::HighConn
		Element::adj_type &               LowConn            (HandleType h) {return *static_cast<inner_reference_array*>(MGetDenseLink(h,LowConnTag()));}
		/// Access directly lower order adjacencies of current element without right of modification.
		Element::adj_type const&          LowConn            (HandleType h) const {return *static_cast<const inner_reference_array*>(MGetDenseLink(h,LowConnTag()));}
		/// Return the size of the array.
		/// For variable size arrays returns current size of the array.
		/// For constant size array returns the same value that may be obtained through GetSize. 
		/// @param h handle of element
		/// @param tag tag that represents the data
		/// @see Mesh::GetSize
		INMOST_DATA_ENUM_TYPE             GetDataSize        (HandleType h,const Tag & tag) const; //For DATA_BULK return number of bytes, otherwise return the length of array
		/// Sets the size of the array for data of variable size.
		/// If you try to change size of data of constant size then if size is
		/// different from current then in debug mode (NDEBUG not set) assertion will fail,
		/// in release mode nothing happens.
		/// @param h handle of element
		/// @param tag tag that represents the data
		/// @param new_size new size of the array
		void                              SetDataSize        (HandleType h,const Tag & tag, enumerator new_size);
		/// Copy inner data array of size elements to provided array begining from shift element.
		/// It is assumed that user-provided array don't overlap inner data.
		/// @param h handle of element
		/// @param tag tag that represents data
		/// @param shift for which element to start to copy
		/// @param size how many elements to copy
		/// @param data user-provided array where data should be copied
		void                              GetData            (HandleType h,const Tag & tag, enumerator shift, enumerator size, void * data) const;
		/// Copy into inner data array of size elements from provided array begining from shift element.
		/// @param h handle of element
		/// @param tag tag that represents data
		/// @param shift for which element to start to copy
		/// @param size how many elements to copy
		/// @param data user-provided array where data should be copied
		void                              SetData            (HandleType h,const Tag & tag, enumerator shift, enumerator size, const void * data);
		/// Remove tag data from given element.
		/// Removes data of variable size and sparse tag data.
		/// Clears to zero data of fixed size.
		/// @param h handle to the element
		/// @param tag tag that indicates the data
		void                              DelData            (HandleType h,const Tag & tag);
		/// Removes data of variable size, clears to zero data of fixed size.
		/// @param h handle to the element
		/// @param tag tag that indicates the data
		void                              DelDenseData       (HandleType h,const Tag & tag);
		/// Removes data of variable size and sparse tag data.
		/// @param h handle to the element
		/// @param tag tag that indicates the data
		void                              DelSparseData      (HandleType h,const Tag & tag);
		/// Check whether data is present on given element.
		/// Always returns true for dense tag data.
		/// @param h handle to the element
		/// @param tag tag that indicates data
		/// @return true, if data exists otherwise false
		bool                              HaveData           (HandleType h,const Tag & tag) const;
		
		Element::GeometricType            GetGeometricType   (HandleType h) const {return static_cast<const bulk *>(MGetDenseLink(h,GeomTypeTag()))[0];}
		void                              SetGeometricType   (HandleType h, Element::GeometricType type) {static_cast<bulk *>(MGetDenseLink(h,GeomTypeTag()))[0] = type;}
		/// Retrive global id of the element with right of modification (dangerous to modify).
		/// Run AssignGlobalID so that tag is automatically allocated and shortcut is set within mesh,
		/// otherwise tag is not created and call will fail.
		/// @param h handle of the element
		/// @return global id
		/// @see Mesh::AssignGlobalID
		integer &                         GlobalID           (HandleType h) {return static_cast<integer *>(MGetDenseLink(h,GlobalIDTag()))[0];}
		/// Retrive global id of the element without right of modification.
		/// Run AssignGlobalID so that tag is automatically allocated and shortcut is set within mesh,
		/// otherwise tag is not created and call will fail.
		/// @param h handle of the element
		/// @return global id
		/// @see Mesh::AssignGlobalID
		integer                           GlobalID           (HandleType h) const {return static_cast<const integer *>(MGetDenseLink(h,GlobalIDTag()))[0];}
		/// Retrive position of the data position of current element, after ReorderEmpty
		/// this number is guaranteed to be between 0 and NumberOf(type of element)
		/// @param h handle of the element
		/// @return local id of data
		integer                        DataLocalID           (HandleType h) const {return static_cast<integer>(links[GetHandleElementNum(h)][GetHandleID(h)]);}
		/// Retrive parallel status of the element.
		/// If mesh is in Serial state then call always returns Element::Owned.
		/// otherwise it will return:
		///    Element::Owned  if there is a single copy of the element on the current processor
		///    Element::Shared if the main copy of the element is located on the current processor
		///    Element::Ghost  if current processor stores dependent copy of the element
		/// @param h handle of the element
		/// @return Element::Status, see function description
		Element::Status                   GetStatus          (HandleType h) const { if( SharedTag().isValid() ) return static_cast<const bulk *>(MGetDenseLink(h,SharedTag()))[0]; return Element::Owned;}
		/// Set parallel status of the element.
		/// If mesh is in Serial state then call will fire asserts in debug mode and segfault in relese mode.
		/// Parallel status controls how exchange of data between elements will be performed,
		/// it is expected that number of elements and copies of elements do match between processors.
		/// This kind of check is never performed and if you fail to setup statuses correctly,
		/// you may end up with data being copied to incorrect elements.
		/// If you modify statuses on your own don't forget to call RecomputeParallelStorage,
		/// otherwise exchange will be performed by old status guidelines.
		/// @param h handle of the element
		/// @param s new status of the element
		void                              SetStatus          (HandleType h, Element::Status s) {BulkDF(h,SharedTag()) = s;}
		//implemented in modify.
		/// Completely destroy element from mesh.
		/// This function bypass check that mesh is in modification state and will remove element immediatly.
		/// It will disconnect element from lower adjacencies and delete all the upper adjacencies that depend on
		/// current element. If you don't want upper adjacencies to be deleted you should first use 
		/// Element::Disconnect function to explicitly disconnect current element and then destroy it.
		/// @param h handle of the element
		/// @see Element::Disconnect
		void                              Destroy            (HandleType h);
		/// Shortcut for typed elements
		void                              Destroy            (const Storage & e) {Destroy(e->GetHandle());}
		/// Hide element from mesh. 
		/// All the functions (except direct access like LowConn,HighConn or ElementSet::getElementsHandles) involving adjacencies
		/// retrival would not return this element.
		/// Works only inside BeginModification and EndModification, on EndModification all Hidden elements are destroyed.
		/// @param h handle of the element
		/// @return if true then element was hidden
		/// @see Mesh::BeginModification
		/// @see Mesh::EndModification
		bool                              Hide               (HandleType h); 
		/// Show element from mesh
		/// @param h handle of the element
		/// @return true then element was recovered
		/// @see Mesh::Hide
		bool                              Show               (HandleType h);
		/// This function will hide element in modification state (between BeginModification and EndModification)
		/// or call Destroy in non-modification state
		/// @param h handle of the element
		/// @return  if true then element was deleted, otherwise it was hidden
		/// @see Mesh::Hide
		/// @see Mesh::BeginModification
		/// @see Mesh::EndModification
		/// @see Mesh::Destroy
		bool                              Delete             (HandleType h);
		/// Check whether element is hidden
		/// @return true if hidden
		/// @see Mesh::Hide
		bool                              Hidden             (HandleType h) const;
		/// Check whether element is new.
		/// Works only in modification state (between BeginModification and EndModification), when you create elements all of them are marked.
		/// @return true if new
		/// @see Mesh::BeginModification
		/// @see Mesh::EndModification
		bool                              New                (HandleType h) const;
		//geometry.cpp
		/// Recompute geometrical type of current element and set it to element.
		/// @param h handle of element
		void                              ComputeGeometricType(HandleType h);
		//mesh.cpp
		/// This function is needed by TagManager, may be made private in future
		/// follows definition of chunk_array to estimate current occupancy of arrays
		INMOST_DATA_ENUM_TYPE             GetArrayCapacity   (integer etypenum);
	private:
		/// Move data position to new location
		void                              MoveStorage        (integer etypenum, integer old_addr, integer new_addr);
		/// Remove data and link position for destroyed element.
		void                              UntieElement       (integer etypenum, integer ID);
		/// Find free data and link positions for new element
		integer                           TieElement         (integer etypenum);
	public:
		//implemented in mesh_parallel.cpp
		enum Action  {AGhost, AMigrate};
		enum Prepare {UnknownSize, UnknownSource};
		typedef void (*ReduceOperation)(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * recv_data, INMOST_DATA_ENUM_TYPE recv_size);
		typedef std::vector<Tag>                                             tag_set;
		typedef std::vector<HandleType>                                  element_set;
		typedef std::vector<INMOST_DATA_BULK_TYPE>                       buffer_type;
		typedef std::map<int, element_set >                            proc_elements;
		typedef std::pair<int, buffer_type >                        proc_buffer_type;
		typedef std::vector< proc_buffer_type >                     exch_buffer_type;
		class exchange_data
		{
		public:
			std::vector<INMOST_MPI_Request> send_reqs, recv_reqs;
			exch_buffer_type send_buffers, recv_buffers;
		};
	private:
		class Random // random generator to provide tag for communication
		{
		private: unsigned int n,a,c,m;
		public:
			Random(unsigned int seed = 50);
			Random(const Random & other);
			unsigned int Number();
		} randomizer;
		class elements_by_type
		{
		private:
			element_set container[4];
		public:
			elements_by_type() {}
			elements_by_type(const elements_by_type & other) {for(int i = 0; i < 4; i++) container[i] = other.container[i];}
			~elements_by_type(){}
			element_set & operator [](int i){ return container[i]; }
			const element_set & operator [](int i) const { return container[i]; }
		};
		typedef std::map<int, elements_by_type > parallel_storage;
	private:
#if defined(USE_PARALLEL_STORAGE)
		parallel_storage                    shared_elements;
		parallel_storage                    ghost_elements;
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		int                                 num_exchanges;
		std::fstream                        out_time;
		int                                 tab;
		int                                 func_id;
#endif
#if defined(USE_MPI_P2P)
		INMOST_MPI_Win                      window;
		unsigned *                          shared_space;
#endif
		int                                 parallel_strategy;
		int                                 parallel_file_strategy;
	private:
		void                              ComputeSharedProcs ();
		proc_elements                     ComputeSharedSkinSet(ElementType bridge);
		void                              PackTagData        (const Tag & tag, const elements_by_type & elements, ElementType mask, MarkerType select, buffer_type & buffer);
		void                              UnpackTagData      (const Tag & tag, const elements_by_type & elements, ElementType mask, MarkerType select, buffer_type & buffer, int & position, ReduceOperation op);
		void                              PackElementsData   (element_set & input, buffer_type & buffer, int destination, const std::vector<std::string> & tag_list);
		void                              UnpackElementsData (element_set & output, buffer_type & buffer, int source, std::vector<std::string> & tag_list);
		void                              PrepareReceiveInner(Prepare todo, exch_buffer_type & send_bufs, exch_buffer_type & recv_bufs);
		void                              ExchangeDataInnerBegin(const tag_set & tag, const parallel_storage & from, const parallel_storage & to, ElementType mask, MarkerType select, exchange_data & storage);
		void                              ExchangeDataInnerEnd(const tag_set & tag, const parallel_storage & from, const parallel_storage & to, ElementType mask, MarkerType select, ReduceOperation op, exchange_data & storage);
		void                              ExchangeBuffersInner(exch_buffer_type & send_bufs, exch_buffer_type & recv_bufs,std::vector<INMOST_MPI_Request> & send_reqs, std::vector<INMOST_MPI_Request> & recv_reqs);
		std::vector<int>                  FinishRequests     (std::vector<INMOST_MPI_Request> & recv_reqs);
		void                              GatherParallelStorage(parallel_storage & ghost, parallel_storage & shared, ElementType mask);
	public:
#if defined(USE_PARALLEL_WRITE_TIME)	
		//this part is needed to test parallel performance
		void                              Enter              ();
		void                              Exit               ();
		int &                             GetFuncID          () {return func_id;}
		std::fstream &                    GetStream          ();
		std::fstream &                    WriteTab           (std::fstream & f);
		void                              FinalizeFile       ();
		static void                       AtExit             (void);
#endif
		/// Initial initialization, calls MPI_Initialize, if MPI was not initialized
		/// it is necessery to invoke this function if you plan to use any parallel algorithms
		/// Accepts arguments passed to console aplication or NULL
		/// @param argc number of arguments for command line
		/// @param argv strings of arguments of command line
		static void                       Initialize         (int * argc, char *** argv);
		/// Finalizes operation with MPI, recomended to call, otherwise MPI may produce warnings
		static void                       Finalize           ();
		/// Set parallel strategy for inner communications
		/// There are three possible scenaries in parallel communication associated in
		/// accordance to enum Prepare structre:
		///   1) The communicating processors and sizes of the messages are known apriori
		///   2) UnknownSize: Communicating processors are known but sizes are unknown
		///   3) unknownSource: Communicationg processors are unknown
		/// Currently with UnknownSize it will run following algorithm
		/// none for strategy 0, following for strategies 1 and 2:
		///   1) Post asynchronous receive with MPI_Irecv of size of buffer to be sent
		///   2) Post asynchronous send with MPI_Isend for required corresponding receive buffer size
		///   3) Wait for all asynchronous operations by MPI_Waitall
		/// With UnknownSource there are two options depending from the USE_MPI_P2P define
		/// If USE_MPI_P2P is defined then MPI-2 api operations will be used
		///   1) MPI_Win_fence to start operations
		///   2) MPI_Put from specially allocated memory location to remote processor 
		///      of array size is performed
		///   3) MPI_Win_fence to stop operations
		/// if USE_MPI_P2P not set then MPI-1 api will be used
		///   1) MPI_Allgather for number of processors to which current processor wants to send data
		///   2) MPI_Allgatherv for sizes and destinations for each processors
		/// Initially it was intended to mainly use MPI-2 functionality for both scenarios but
		/// it was realized that there is no availible hardware on which MPI-2 functionalty
		/// performs match better then MPI-1 counterparts, especially in the case of UnknownSize.
		/// Probably this happens due to lack of support of RDMA operations.
		/// If you believe it will be better to use MPI-2 in both cases you are free to uncomment
		/// definition of PREFFER_MPI_P2P in inmost_common.h then MPI-2 will be used in both scenaries.
		/// These algorthms above are implemented in Mesh::ExchangeBufferInner.
		/// After first problem was resolved following strategies are availible for main communication:
		/// strategy = 0
		///   1) Asyncronous send of data by MPI_Isend;
		///   2) Check incoming messages by MPI_Probe;
		///   3) Check incoming message size by MPI_Get_size;
		///   4) Allocation of buffers of required size;
		///   5) Asynchronous receive of data
		///   6) MPI_Waitsome to copy received results to buffers
		/// This strategy shows to be the fastest on mac with intel core 2 duo
		/// it appears to be independend on apriori knowledge of sizes of incoming
		/// messages and skips the step of determining sizes in all the cases but still
		/// it requires establishing knowledge of communicating processors
		/// Asynchronous sending and receiving may be performed by breaking
		/// the steps 1) and 2-5) but should be considered bad since it
		/// will be performed without appropriate receive buffers posted for sends, as a result
		/// messages will stuck in network pipeline and would be repeatedly rejected
		/// resulting in bad networking performance especially if processors have small memory.
		/// As a result non-asynchronous communication is realized with this stategy breaking
		/// steps 1-5) and 6) when you as for asynchronous communication
		/// startegy = 1
		///   1) Post asynchronous receive of data by MPI_Irecv
		///   2) Post asynchronous send of data by MPI_Isend
		///   3) MPI_Waitsome for any received data
		/// True asynchronous behavior is reproduced by breaking 1-2) and 3)
		/// strategy = 2
		///   1) Post asynchronous receive of data by MPI_Irecv
		///   2) Set MPI_Barrier to ensure that all the receives were properly set by the time
		///   2) Perform direct send of data by MPI_Irsend
		///   3) MPI_Waitsome for any received data
		/// For asynchronous communication algorithm is broken into 1-3) and 4) which is fairly
		/// asynchronous. The only provisional benefit it may have on machines with small memory 
		/// since it should bypass any allocation of buffers for sent and received data by MPI and 
		/// probably not perform any randezvous communication to ensure data allocation.
		/// But MPI_Barrier looks like elephant here.
		/// Algorithms above are implemented in Mesh::ExchangeBuffersInner
		/// @see Mesh::PrepareReceiveInner
		/// @see Mesh::ExchangeBuffersInner
		void                              SetParallelStrategy(int strategy){assert( !(strategy < 0 || strategy > 3) ); parallel_strategy = strategy;}
		/// Retrive currently set parallel strategy
		/// @see Mesh::SetParallelStrategy
		int                               GetParallelStrategy() {return parallel_strategy;}
		/// This strategy correspond only to internal ".pmf" mesh format
		/// there are two availible strategies for ".pmf" files loading and saving:
		/// strategy = 0
		///  on save 
		///   1) every processor gather local data to some buffer
		///   2) MPI_Gather to obtain sizes of data among processors
		///   3) MPI_Gatherv to obtain the whole data on zeros processor
		///   4) the first processor writes all the data to disk by std::fstream
		///  on load
		///   1) first processors reads the whole file by std::fstream
		///   2) MPI_Scatter distributes block sizes among processors
		///   3) MPI_Scatterv distributes blocks among processors
		///   4) Each processor parses it's block
		/// This strategy requires one processor to hold all the data, which
		/// is quite bad for large files. New strategy may be created from this one in future
		/// when each processors consequently obtain access to the file using std::fstream and
		/// writes the data.
		/// strategy = 1
		///  on save it will perform:
		///   1) MPI_Gather to obtain sizes of data among processors on zeroth processor
		///   2) MPI_File_open to get parallel handle for the file
		///   3) MPI_File_write_shared called by processor with zeroth rank to write header
		///   4) MPI_File_write_ordered to write contents of individual data
		///   5) MPI_File_close to close parallel file handle
		///  on load it will perform
		///   1) MPI_File_open to open the file in parallel
		///   2) MPI_File_read_shared to get contents of header on zeroth processor
		///   3) MPI_Scatter to distribute block sizes among processors
		///   4) MPI_File_read_ordered to obtain contents
		///   5) MPI_File_close to close parallel file handle
		///  Availible only when USE_MPI_P2P is set because it rely on MPI-2 api that begins with MPI_File_xxx
		///  some MPI-1 standards contain this api as extension.
		/// The strategy 1 appeared to be considerably slower on INM cluster then strategy 0, this may
		/// happen due to lack of read-write devices that able to work in parallel. On IBM Bluegene/p
		/// strategy 1 was not working due to same old problem with shared file pointers in their MPI realization
		void                              SetParallelFileStrategy(int strategy){assert( !(strategy < 0 || strategy > 1) ); parallel_file_strategy = strategy;}
		/// Retrive currently set parallel strategy for ".pmf" files
		/// @see Mesh::GetParallelStrategy
		int                               GetParallelFileStrategy() {return parallel_file_strategy;}
		/// Get rank of current processor
		int                               GetProcessorRank   ();
		/// Get number of processors
		int                               GetProcessorsNumber();
		/// Retrive MPI communicator
		INMOST_MPI_Comm                   GetCommunicator    ();
		/// Set MPI communicator
		void                              SetCommunicator    (INMOST_MPI_Comm _comm);
		void                              ResolveShared      ();
		/// Delete all the ghost cells.
		void                              RemoveGhost        ();
		/// Delete some ghost cells provided in array.
		/// Non-ghost elements will also be deleted.
		/// This algorithm will properly communicate between processors so that 
		/// parallel states of deleted elements properly updated on remote processors.
		/// Internally function invokes Destroy function, not Delete, which hides elements during modification state,
		/// currently it is not expected that any parallel algorithms will be performed between BeginModification and EndModification
		/// since modification may break parallel state though it is never checked whether the mesh is in the modification state,
		/// so you are free to experiment. This behavior may change in future.
		/// Collective point-2-point.
		/// @param ghost array of handles
		/// @param num number of handles
		/// TODO
		///  1) Currently request for deletion of elements of lower level then cell will be simply ignored, ensure
		///     in future that algorithm will properly rise deletion data from lower to upper adjacencies 
		///     to delete all the upper adjacencies that depend on deleted lower adjacencies
		void                              RemoveGhostElements(const HandleType * ghost, enumerator num);
		template<typename EType>
		void                              RemoveGhostElements(const ElementArray<EType> & ghost) {RemoveGhostElements(ghost.data(),static_cast<enumerator>(ghost.size()));}
		void                              RemoveGhostElements(const ElementSet & ghost) {RemoveGhostElements(ghost.getHandles(),ghost.nbHandles());}
		/// Assign unique numbers to elements, internally this will create Mesh::GlobalIDTag and 
		/// make call to Element::GlobalID and Mesh::GlobalID functions valid. Internally this 
		/// will also set have_global_id variable that will indicate that all the comparisons 
		/// in parallel algorithms should be performed using global identificators instead 
		/// of centroids which is faster.
		/// TODO
		///     1) invoking function before loading mesh will not renew global identificators after load
		///        but would not unset have_global_id either. There are probably too many places when
		///        global ids may become invalid but no flag will be set. It may be benefitial to set
		///        such flags along with updating geometrical data which seems to be maintained fairly well
		///        during mesh modification
		void                              AssignGlobalID     (ElementType mask);
		/// Update data from Shared elements to Ghost elements, for backward direction please see Mesh::ReduceData
		/// If you have a tag of DATA_BULK type and you store your own custom data structure in it, it is highly
		/// recomended that you provide MPI information about your structure through Tag::SetBulkDataType,
		/// this would not do any difference on homogeneous architecture, but it may help you save a lot of 
		/// time and nerves in heterogeneous parallel environment.
		/// TODO: see TODO in Mesh::ReduceData
		/// Blocking, Collective point-2-point
		/// @param tag tag that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @see Mesh::ReduceData
		void                              ExchangeData       (const Tag & tag, ElementType mask, MarkerType select);
		/// Start asynchronous synchronization of data.
		/// You should define object of type exchange_data that will hold temporary buffers for data.
		/// every Mesh::ExchangeDataBegin should be matched with Mesh::ExchangeDataEnd with the same
		/// exchange_data object. After matching Mesh::ExchangeDataEnd the exchange_data object may be reused
		/// If you will go out of the scope where exchange_data object was defined it will be deallocated
		/// and may result in segmentation fault.
		/// You should also never put the same exchange_data object to any other Mesh::ExchangeDataBegin or
		/// Mesh::ReduceDataBegin, until matching Mesh::ExchangeDataEnd because it may override or reallocate 
		/// buffers, internally used by MPI and remote processor will receive garbage instead of data.
		/// Nonblocking, Collective point-2-point
		/// @param tag tag that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		/// @see Mesh::ExchangeDataEnd
		void                              ExchangeDataBegin  (const Tag & tag, ElementType mask, MarkerType select, exchange_data & storage);
		/// Complete asynchronous synchronization of data.
		/// Blocking
		/// @param tag tag that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		/// @see Mesh::ExchangeDataBegin
		void                              ExchangeDataEnd    (const Tag & tag, ElementType mask, MarkerType select, exchange_data & storage);
		/// This function will perform exchange of multiple data tags.
		/// Blocking, Collective point-2-point
		/// @param tags multiple tags that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		void                              ExchangeData       (const tag_set & tags, ElementType mask, MarkerType select);
		/// This function will initialize exchange of multiple data tags, 
		/// using this function may lead to good overlapping between communication and computation.
		/// Nonblocking, Collective point-2-point
		/// @param tags multiple tags that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		void                              ExchangeDataBegin  (const tag_set & tags, ElementType mask, MarkerType select, exchange_data & storage);
		/// This function will finalize exchange of multiple data tags.
		/// Blocking
		/// @param tags multiple tags that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		void                              ExchangeDataEnd    (const tag_set & tags, ElementType mask, MarkerType select, exchange_data & storage);
		/// Accumulation of data from ghost elements to shared elements, accumulation is performed 
		/// based on user-provided function. When processor - owner of the element receives data addressed
		/// to this element then it calls user-defined op function. Since there may be multiple ghost
		/// elements per one shared element, function may be called multiple times.
		/// Several examples of reduction functions may be found within the mesh_parallel.cpp source.
		/// Remember that the result will be up to date only on the owner of the processor. You will have
		/// to run Mesh::ExchangeData to make the data up to date among all of the processors.
		/// If you have a tag of DATA_BULK type and you store your own custom data structure in it, it is highly
		/// recomended that you provide MPI information about your structure through Tag::SetBulkDataType,
		/// this would not do any difference on homogeneous architecture, but it may help you save a lot of 
		/// time and nerves in heterogeneous parallel environment.
		/// Exchanging tags of DATA_REFERNCE is not implemented, TODO 14.
		/// TODO:
		///    1) Exchanging DATA_REFERENCE tags not implemented, this is due to the absence of any conclusion
		///       on how it should behave:
		///         either only search within elements owned by the other processor and
		///         establish references and set InvalidHandle() to elements that are not found (fairly easy,
		///         will involve search operations to check against owned elements for similar entry, efficient
		///         implementation will require bounding search trees (see TODO 49);
		///       or:
		///         send all the referenced elements through PackElementsData and establish all the links within
		///         elements reproduced by UnpackElementsData (UnpackElementsData calls UnpackTagData with set
		///         of unpacked elements using which it will be vary comfortable to establish references on
		///         remote processor)
		/// Blocking, Collective point-2-point
		/// @param tag tag that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param op user-defined operation on received data
		/// @see Mesh::ExchangeData
		void                              ReduceData         (const Tag & tag, ElementType mask, MarkerType select, ReduceOperation op );
		/// This function intializes data reduction.
		/// Read recomendations about exchange_storage object in Mesh::ExchangeDataBegin
		/// Nonblocking, Collective point-2-point
		/// @param tag tag that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		/// @see Mesh::ReduceData
		/// @see Mesh::ExchangeDataBegin
		void                              ReduceDataBegin    (const Tag & tag, ElementType mask, MarkerType select, exchange_data & storage);
		/// This function completes data reduction/
		/// Read Mesh::ReduceData for information about op function
		/// Blocking
		/// @param tag tag that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		/// @param op user-defined operation on received data
		/// @see Mesh::ReduceData
		void                              ReduceDataEnd      (const Tag & tag, ElementType mask, MarkerType select, ReduceOperation op, exchange_data & storage );
		/// This function will perform reduction of multiple data tags.
		/// Note that function will be the same for all tags, you can differentiate behavior of function
		/// depending on tag name (may be expensive).
		/// Blocking, collective point-2-point
		/// @param tag multiple tags that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param op user-defined operation on received data
		void                              ReduceData         (const tag_set & tags, ElementType mask, MarkerType select, ReduceOperation op );
		/// This function will initialize reduction of multiple data tags.
		/// Nonblocking, collective point-2-point
		/// @param tags multiple tags that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		void                              ReduceDataBegin    (const tag_set & tags, ElementType mask, MarkerType select, exchange_data & storage);
		/// This function will finalize exchange of multiple data tags.
		/// Blocking
		/// @param tag multiple tags that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		/// @param op user-defined operation on received data
		void                              ReduceDataEnd      (const tag_set & tags, ElementType mask, MarkerType select, ReduceOperation op, exchange_data & storage );
		/// This function realizes two algorithms: ghosting of elements and migration of elements
		/// ghosting:
		///  Creates ghosted elements at other processors prescribed in SendtoTag() performs
		///  sending and receiving of elements and all operations to keep parallel state of the mesh.
		///   Given that all the data was up to date among processors all the data at the end of the
		///   algorithm will be also up to data
		/// migration:
		///  To correctly perform migration of elements it is necessery to set up tags besides SendToTag(),
		///  which indicates to where every element should be sent:
		/// tag "TEMPORARY_NEW_PROCESSORS", of type DATA_INTEGER with variable size, tells which processors will have copy of the element after migration;
		/// tag "TEMPORARY_NEW_OWNER", of type DATA_INTEGER of size 1, tells which processor will have main copy of the element.
		/// if there is no current processor in "TEMPORARY_NEW_PROCESSORS", then current processor will remove copy of the element.
		/// All this actions are performed automatically by Mesh::Redistribute based on information provided in Mesh::RedistributeTag
		/// which effectively contains new owner.
		///   Given that all the data was up to date among processors all the data at the end of the algorithm will be also up to data
		/// TODO 
		///     1) test halo exchange algorithm (if used then change collective point-2-point to collective)
		///     2) see TODO 2 in Mesh::Redistribute
		/// Collective point-2-point.
		/// @param action
		/// @see Mesh::SendtoTag
		/// @see Mesh::Redistribute
		void                              ExchangeMarked     (enum Action action = AGhost);
		/// Form several layers of ghosted cells that are adjacent through bridge elements to current cells.
		/// This function acceptes any mesh topology, failure on some mesh should be considered a bug and
		/// sample example should be provided for testing purposes.
		/// This function internally calculates layer by layer and invokes ExchangeMarked for each layer,
		/// you can either reproduce the algorithm on your own if you want to bypass the function and
		/// call Mesh::ExchangeMarked directly, but then you will lose optimization in Mesh::Redistribute,
		/// that expects that layers are formed the same way they are formed in Mesh::ExchangeGhost.
		/// Internally it sets up LayersTag and BridgeTag for the mesh with provided values, which
		/// are used by Mesh::Redistribute but you are discouraged to override these tags since using
		/// non-matching algorithms is not tested and should be considered dangerous.
		/// Nevertheless you can use this function first for layers then request any additional ghosted elements
		/// by ExchangeMarked.
		/// Collective point-2-point.
		/// @param layers number of required layers of ghosted elements
		/// @param bridge bitwise mask of elements for which neighbouring cells should be considered a layer
		/// @see Mesh::ExchangeMarked
		/// @see Mesh::Redistribute
		void                              ExchangeGhost      (integer layers, ElementType bridge);
		/// Migrate all the elements to the new owners prescribed in data corresponding to RedistributeTag.
		/// This will perform all the actions to send mesh elements and data and reproduce new mesh partitions
		/// on remote elements and correctly resolve parallel state of the mesh. If you have priviously
		/// prescribed number of layers through ExchangeGhost, then minimal number of actions will be performed
		/// to reproduce layers of ghosted elements whithout involving removal of all ghosted elements.
		/// Internally function sets up following data on elements using provided information:
		/// "TEMPORARY_NEW_PROCESSORS" - new set processors that contain copy of the element
		/// "TEMPORARY_NEW_OWNER"      - new owner for each processor (effectively RedistributeTag)
		///
		/// Action of this function regarding restoration of layers of ghosted elements in the case you have 
		/// modified mesh without involving Mesh::ResolveModification is yet to be tested and should be
		/// considered dangerous.
		///
		/// If you have output from Zoltan or ParMetis  for cells of the mesh
		/// then just write this output to RedistributeTag and call Mesh::Redistribute.
		///
		/// TODO: 1)introduce "TEMPORARY_KEEP_GHOSTED" tag that will store processors on which copy of element
		///         should be kept, internally just merge it with "TEMPORARY_NEW_PROCESSORS" tag
		///         this will allow user to control ghosting of certain elements and not to invoke ExchangeMarked
		///         every time after Redistribute.
		///         This is probably already done using Mesh::SendtoTag, because function fills it without
		///         clearing and ExchangeMarked performs initial action based on SendtoTag, it is due to
		///         check that SendtoTag is properly merged with "TEMPORARY_NEW_PROCESSORS" before call to ExchangeMarked
		///         and received elements are not deleted by accident.
		///       2)let user provide any integer tag as input without involving RedistributeTag
		/// Collective point-2-point.
		/// @see Mesh::RedistributeTag
		/// @see Mesh::ExchangeGhost
		void                              Redistribute       ();
		/// Enumerate all elements begining with start and put numeration to data associated with num_tag for all elements with given type mask
		/// Collective operation.
		/// @param mask bitwise or of types of elements
		/// @param num_tag a tag that is associated with the data
		/// @param start starting value for enumeration
		/// @param define_sparse if true then function will define sparse data on elements that don't have it, otherwise it will skip those elements
		/// @return last value on all the processors
		integer                           Enumerate          (ElementType mask, Tag num_tag, integer start = 0, bool define_sparse = false);
		/// Enumerate all elements begining with start and put numeration to data associated with num_tag 
		/// Collective operation.
		/// @param h array of handles
		/// @param num number of handles
		/// @param num_tag a tag that is associated with the data
		/// @param start starting value for enumeration
		/// @param define_sparse if true then function will define sparse data on elements that don't have it, otherwise it will skip those elements
		/// @return last value on all the processors
		integer                           Enumerate          (const HandleType * h, enumerator num, const Tag & num_tag, integer start = 0, bool define_sparse = true);
		/// Enumerate all elements begining with start and put numeration to data associated with num_tag 
		/// Collective operation.
		/// @param elements array of elements
		/// @param num_tag a tag that is associated with the data
		/// @param start starting value for enumeration
		/// @param define_sparse if true then function will define sparse data on elements that don't have it, otherwise it will skip those elements
		/// @return last value on all the processors
		template<typename EType>
		integer                           Enumerate          (const ElementArray<EType> & elements, const Tag & num_tag, integer start = 0, bool define_sparse = true) {return Enumerate(elements.data(),static_cast<enumerator>(elements.size()),num_tag,start,define_sparse);}
		/// Enumerate all elements in the set
		/// Collective operation.
		/// @param set handle of the set
		/// @param num_tag a tag that is associated with the data
		/// @param start starting value for enumeration
		/// @param define_sparse if true then function will define sparse data on elements that don't have it, otherwise it will skip those elements
		/// @return last value on all the processors
		integer                           EnumerateSet       (const ElementSet & set, const Tag & num_tag, integer start = 0, bool define_sparse = true);
		/// Sum of all physical elements, it excludes ghosted copies.
		/// To compute total number including ghosted copies run Integrate(NumberOf(mask))
		/// Collective operation.
		/// @param mask bitwise mask of element types, example: NODE | CELL
		/// @return sum of elements
		/// @see Mesh::NumberOf
		integer                           TotalNumberOf      (ElementType mask);
		/// Integrate real value over all processors
		/// Collective operation.
		/// @param input value on current processor
		/// @return sum over all processors
		real                              Integrate          (real input);
		/// Integrate integer value over all processors.
		/// Collective operation.
		/// @param value on current processor
		/// @return sum over all processors
		integer                           Integrate          (integer input);
		/// Integrate data corresponding to tag between all processors.
		/// Elements without the data defined on them or when entry not present will be skipped.
		/// Collective operation.
		/// @param t tag that correspond to data to be integrated
		/// @param entry in the array of data
		/// @param mask bitwise or of types of elements on which to integrate
		/// @return sum between all processors
		real                              Integrate          (const Tag & t, enumerator entry, ElementType mask);
		/// Compute sum of integer values for all processors with rank lower then current, excluding current processor
		/// Collective operation.
		/// @param input value on current processor
		/// @return described sum
		integer                           ExclusiveSum       (integer input);
		real                              AggregateMax       (real input);
		integer                           AggregateMax       (integer input);
		/// Regather ghosted and shared element sets for data exchange, 
		/// this function will be quite useful if you change statuses of elements
		/// or modify mesh on your own bypassing internal algorithms.
		/// No action will be performed if USE_PARALLEL_STORAGE is not set in inmost_common.h,
		/// since all the elements are computed during exchange phase.
		/// Generally this is not needed if you use high-level algorithms for mesh modification
		/// or mesh redistribution.
		/// @param bitwise type mask
		/// @see Mesh::BeginModification
		/// @see Mesh::EndModification
		/// @see Mesh::ExchangeMarked
		/// @see Mesh::RemoveGhostElements
		void                              RecomputeParallelStorage(ElementType mask);
		/// Synchronize bitwise mask of element types between processors
		/// Collective operation
		/// @param bitwise type mask
		/// @return bitwise result among processors
		ElementType                       SynchronizeElementType(ElementType etype);
		/// Syncronize marker on elements between processors using provided operation
		/// Depending on requested operation following action is performed:
		/// SYNC_BIT_SET - value on ghost elements is set by value on corresponding shared processors;
		/// SYNC_BIT_OR  - bitwise OR between values in ghosted and shared elements;
		/// SYNC_BIT_AND - bitwise AND between values in ghosted and shared elements;
		/// SYNC_BIT_XOR - bitwise XOR between values in ghosted and shared elements.
		/// @param marker marker to be synchronized
		/// @param mask bitwise or type mask
		/// @param op operation, one of SYNC_BIT_SET, SYNC_BIT_OR, SYNC_BIT_XOR, SYNC_BIT_AND
		void                              SynchronizeMarker  (MarkerType marker, ElementType mask, SyncBitOp op);	
		//for debug
		void                              BeginSequentialCode();
		void                              EndSequentialCode  ();
		//iterator.cpp::::::::::::::::::::::::::::::::::::::::::::::::::
	public:
		Element                           ElementByLocalID   (integer etypenum, integer lid) {assert(etypenum < 5 && (lid >= 0 && lid < static_cast<integer>(links[etypenum].size())) || (etypenum == 5 && lid == 0)); return Element(this,ComposeHandle(etypenum,lid));}
		Element                           ElementByLocalID   (ElementType etype, integer lid) {return ElementByLocalID(ElementNum(etype),lid);}
		Element                           ElementByHandle    (HandleType h) {return Element(this,h);}
		
		HandleType                        NextHandle         (HandleType h) const;
		HandleType                        PrevHandle         (HandleType h) const; //returns InvalidHandle() when go beyond first element
		HandleType                        NextHandle         (HandleType h, ElementType mask) const;
		HandleType                        PrevHandle         (HandleType h, ElementType mask) const; //returns InvalidHandle() when go beyond first element
		HandleType                        FirstHandle        () const {return ComposeHandle(ElementNum(NODE),0);}
		HandleType                        LastHandle         () const {return ComposeHandle(ElementNum(MESH),1);} 
		HandleType                        FirstHandle        (ElementType etype) const {return ComposeHandle(ElementNum(etype),0);}
		HandleType                        LastHandle         (ElementType etype) const  {integer num = ElementNum(etype); return ComposeHandle(num,static_cast<integer>(links[num].size()));}

		Node                              NodeByLocalID      (integer lid) { assert(lid >= 0 && lid < static_cast<integer>(links[0].size())); return Node(this,ComposeHandle(0,lid)); }
		Edge                              EdgeByLocalID      (integer lid) { assert(lid >= 0 && lid < static_cast<integer>(links[1].size())); return Edge(this,ComposeHandle(1,lid)); }
		Face                              FaceByLocalID      (integer lid) { assert(lid >= 0 && lid < static_cast<integer>(links[2].size())); return Face(this,ComposeHandle(2,lid));}
		Cell                              CellByLocalID      (integer lid) { assert(lid >= 0 && lid < static_cast<integer>(links[3].size())); return Cell(this,ComposeHandle(3,lid)); }
		ElementSet                        EsetByLocalID      (integer lid) { assert(lid >= 0 && lid < static_cast<integer>(links[4].size())); return ElementSet(this,ComposeHandle(4,lid)); }
		
		integer                           NodeNextLocalID    (integer lid) const {++lid; while(lid < static_cast<integer>(links[0].size()) && links[0][lid] == -1) ++lid; return lid;}
		integer                           EdgeNextLocalID    (integer lid) const {++lid; while(lid < static_cast<integer>(links[1].size()) && links[1][lid] == -1) ++lid; return lid;}
		integer                           FaceNextLocalID    (integer lid) const {++lid; while(lid < static_cast<integer>(links[2].size()) && links[2][lid] == -1) ++lid; return lid;}
		integer                           CellNextLocalID    (integer lid) const {++lid; while(lid < static_cast<integer>(links[3].size()) && links[3][lid] == -1) ++lid; return lid;}
		integer                           EsetNextLocalID    (integer lid) const {++lid; while(lid < static_cast<integer>(links[4].size()) && links[4][lid] == -1) ++lid; return lid;}
		
		integer                           NodePrevLocalID    (integer lid) const {--lid; while(lid >= 0 && links[0][lid] == -1) --lid; return lid;}
		integer                           EdgePrevLocalID    (integer lid) const {--lid; while(lid >= 0 && links[1][lid] == -1) --lid; return lid;}
		integer                           FacePrevLocalID    (integer lid) const {--lid; while(lid >= 0 && links[2][lid] == -1) --lid; return lid;}
		integer                           CellPrevLocalID    (integer lid) const {--lid; while(lid >= 0 && links[3][lid] == -1) --lid; return lid;}
		integer                           EsetPrevLocalID    (integer lid) const {--lid; while(lid >= 0 && links[4][lid] == -1) --lid; return lid;}
		
		__INLINE integer                  NodeLastLocalID    () const {return static_cast<integer>(links[0].size());}
		__INLINE integer                  EdgeLastLocalID    () const {return static_cast<integer>(links[1].size());}
		__INLINE integer                  FaceLastLocalID    () const {return static_cast<integer>(links[2].size());}
		__INLINE integer                  CellLastLocalID    () const {return static_cast<integer>(links[3].size());}
		__INLINE integer                  EsetLastLocalID    () const {return static_cast<integer>(links[4].size());}
		integer                           NextLocalID        (ElementType etype, integer lid) const {integer q = ElementNum(etype); ++lid; while(lid < static_cast<integer>(links[q].size()) && links[q][lid] == -1) ++lid; return lid;}
		integer                           PrevLocalID        (ElementType etype, integer lid) const {integer q = ElementNum(etype); --lid; while(lid > 0 && links[q][lid] == -1) --lid; return lid;}
		integer                           FirstLocalID       (ElementType etype) const;
		integer                           LastLocalID        (integer n) const {assert(n >= 0 && n < 6); return n == 5 ? 1 : static_cast<integer>(links[n].size());}
		integer                           LastLocalID        (ElementType etype) const {assert(OneType(etype)); return LastLocalID(ElementNum(etype));}


		__INLINE integer                  NumberOfSets       () const { return static_cast<integer>(links[4].size() - empty_links[4].size()); }
		__INLINE integer                  NumberOfCells      () const { return static_cast<integer>(links[3].size() - empty_links[3].size());}
		__INLINE integer                  NumberOfFaces      () const { return static_cast<integer>(links[2].size() - empty_links[2].size()); }
		__INLINE integer                  NumberOfEdges      () const { return static_cast<integer>(links[1].size() - empty_links[1].size()); }
		__INLINE integer                  NumberOfNodes      () const { return static_cast<integer>(links[0].size() - empty_links[0].size()); }
		__INLINE integer                  NumberOfElements   () const { return NumberOfCells() + NumberOfFaces() + NumberOfEdges() + NumberOfNodes(); }
		__INLINE integer                  NumberOfAll        () const { return NumberOfSets() + NumberOfElements(); }
		integer                           NumberOf           (ElementType t) const;

		template<typename EType> class base_iterator;
		typedef base_iterator<Storage>    iteratorStorage;
		typedef base_iterator<Element>    iteratorElement;
		typedef base_iterator<ElementSet> iteratorSet;
		typedef base_iterator<Cell>       iteratorCell;
		typedef base_iterator<Face>       iteratorFace;
		typedef base_iterator<Edge>       iteratorEdge;
		typedef base_iterator<Node>       iteratorNode;
		/// These iterators skip invalid elements but don't skip modified elements.
		iteratorStorage                   Begin              (ElementType Types);
		iteratorStorage                   End                ();
		iteratorElement                   BeginElement       (ElementType Types);
		iteratorElement                   EndElement         ();
		iteratorSet                       BeginSet           ();
		iteratorSet                       EndSet             ();
		iteratorCell                      BeginCell          ();
		iteratorCell                      EndCell            ();
		iteratorFace                      BeginFace          ();
		iteratorFace                      EndFace            ();
		iteratorEdge                      BeginEdge          ();
		iteratorEdge                      EndEdge            ();
		iteratorNode                      BeginNode          ();
		iteratorNode                      EndNode            ();


		template<typename EType>
		class base_iterator
		{
		protected:
			Mesh *               m;
			Storage::integer     lid;
			ElementType          etype;
			ElementType          types;
			base_iterator(ElementType Types, Mesh * mesh, bool last);
			base_iterator(Mesh * mesh) {m = mesh; etype = NONE; types = NONE; lid = -1;}
		public:
			typedef HandleType *                     pointer;
			typedef HandleType &                     reference;
			typedef HandleType                       value_type;
			typedef ptrdiff_t                        difference_type;
			typedef std::bidirectional_iterator_tag  iterator_category;
			base_iterator() {m = NULL; lid = -1; etype = types = NONE;}
			base_iterator(const base_iterator & other) {m = other.m; lid = other.lid; types = other.types; etype = other.etype;}
			virtual ~base_iterator() {}
			base_iterator &             operator ++();
			__INLINE base_iterator      operator ++(int) {Mesh::base_iterator<EType> ret(*this); operator++(); return ret;}
			base_iterator &             operator --();
			__INLINE base_iterator      operator --(int) {Mesh::base_iterator<EType> ret(*this); operator--(); return ret;}
			__INLINE value_type         operator * () {return ComposeHandle(etype,lid);}
			__INLINE EType              operator ->() {return EType(m,ComposeHandle(etype,lid));}
			__INLINE base_iterator &    operator = (base_iterator const & other) {m = other.m; lid = other.lid; types = other.types; etype = other.etype; return *this;}
			__INLINE bool               operator ==(const base_iterator & other) const {return lid == other.lid && etype == other.etype;}
			__INLINE bool               operator !=(const base_iterator & other) const {return lid != other.lid || etype != other.etype;}
			__INLINE bool               operator < (const base_iterator & other) const {return (etype < other.etype) || (etype == other.etype && lid <  other.lid);}
			__INLINE bool               operator > (const base_iterator & other) const {return (etype > other.etype) || (etype == other.etype && lid >  other.lid);}
			__INLINE bool               operator <=(const base_iterator & other) const {return (etype < other.etype) || (etype == other.etype && lid <= other.lid);}
			__INLINE bool               operator >=(const base_iterator & other) const {return (etype > other.etype) || (etype == other.etype && lid >= other.lid);}
			void Print();
			friend iteratorStorage Mesh::Begin(ElementType Types);
			friend iteratorStorage Mesh::End();
			friend iteratorElement Mesh::BeginElement(ElementType Types);
			friend iteratorElement Mesh::EndElement();
			friend iteratorSet     Mesh::BeginSet();
			friend iteratorSet     Mesh::EndSet();
			friend iteratorCell    Mesh::BeginCell();
			friend iteratorCell    Mesh::EndCell();
			friend iteratorFace    Mesh::BeginFace();
			friend iteratorFace    Mesh::EndFace();
			friend iteratorEdge    Mesh::BeginEdge();
			friend iteratorEdge    Mesh::EndEdge();
			friend iteratorNode    Mesh::BeginNode();
			friend iteratorNode    Mesh::EndNode();
		};
	private:
		std::vector< std::pair<std::string, std::string> > file_options;
	public:
		/// Current availible file options:
		/// "VTK_GRID_DIMS" - set "2" for two-dimensional vtk grids, "3" for three-dimensional vtk grids
		/// "VERBOSITY"     - set "2" for progress messages, "1" for reports, "0" for silence
		/// TODO:
		///      introduce "SET_TAGS_LOAD", "SET_TAGS_SAVE" to explicitly provide set of tags to write
		///      or "SKIP_TAGS_LOAD", "SKIP_TAGS_SAVE" tags to skip
		void         SetFileOption(std::string,std::string);
		/// Get current option corresponding to key
		/// @param key options for which options should be retrieven
		std::string  GetFileOption(std::string key);
		/// Acceptable file formats for reading
		/// ".vtk"    - legacy vtk format for unstructured grid
		/// ".pvtk"   - legacy parallel vtk format
		/// ".gmv"    - format acceptable by general mesh viewer
		/// ".msh"    - gmsh generator format
		/// ".grdecl" - eclipse format (under construction)
		/// ".grid"   - mesh format by Mohammad Karimi-Fard
		/// ".pmf"    - internal parallel portable binary format, saves all features
		/// @param File path to the file
		void         Load(std::string File); 
		/// Acceptable file formats for writing
		/// ".vtk"  - legacy vtk format for unstructured grid
		/// ".pvtk" - legacy parallel vtk format
		/// ".gmv"  - format acceptable by general mesh viewer
		/// ".pmf"  - internal parallel portable binary format, saves all features
		/// Remeber: .pmf stores all references to elements. If reference are broken due to mesh modification,
		///          saving or loading such a mesh may lead to seagfault. To automatically maintain correct
		///          references modify mesh using BeginModification, ApplyModification, EndModification
		/// @param File path to the file
		void         Save(std::string File);
		bool         isParallelFileFormat(std::string File);
	public:
		
		//implemented in geometry.cpp
	private:
		Tag          measure_tag;
		Tag          centroid_tag;
		Tag          normal_tag;
		Tag          barycenter_tag;
		Tag          boundary_tag;
		bool         remember[5][3];
	private:
		void                              RestoreGeometricTags();
		bool                              HideGeometricData  (GeometricData type, ElementType mask) {return remember[type][ElementNum(mask)-1] = false;}
		bool                              ShowGeometricData  (GeometricData type, ElementType mask) {return remember[type][ElementNum(mask)-1] = true;}
	public:
		typedef tiny_map<GeometricData, ElementType,5> GeomParam;
		// types for MEASURE:     EDGE | FACE | CELL   (length, area, volume)
		// types for CENTROID:    EDGE | FACE | CELL
		// types for BARYCENTER:  EDGE | FACE | CELL
		// types for NORMAL:      FACE | CELL          (may precompute normal for cells in 2d case)
		// types for ORIENTATION: FACE
		void                              PrepareGeometricData(GeomParam table);
		void                              RemoveGeometricData(GeomParam table);
		bool                              HaveGeometricData  (GeometricData type, ElementType mask) const {return remember[type][ElementNum(mask)-1];} // requests to only one geometric and element type allowed
		void                              GetGeometricData   (HandleType e, GeometricData type, real * ret);
		const Tag &                       GetGeometricTag    (GeometricData type) const;
		bool                              TestClosure        (const HandleType * elements, integer num) const;
		ElementArray<Face>                GatherBoundaryFaces();
		ElementArray<Face>                GatherInteriorFaces();
		integer                           CountBoundaryFaces ();
		integer                           CountInteriorFaces ();
		void                              RecomputeGeometricData(HandleType e); // Update all stored geometric data, runs automatically on element creation
		Element::GeometricType            ComputeGeometricType(ElementType element_type, const HandleType * lower_adjacent, INMOST_DATA_ENUM_TYPE lower_adjacent_size) const;
		//implemented in modify.cpp
	private:
		MarkerType hide_element, new_element;
	public:
		/// Check weather code runs between Mesh::BeginModification, Mesh::EndModification scope
		/// In case mesh is modified, on element creation Mesh::TieElements will always place elements 
		/// to the end of the array as a result all the newly created elements will be iterated after current
		/// or hidden elements.
		bool                              isMeshModified     () const {return new_element != 0;} 
		MarkerType                        HideMarker         () const {return hide_element;}
		MarkerType                        NewMarker          () const {return new_element;}
		void                              SwapModification   (); // swap hidden and new elements, so that old mesh is recovered
		void                              BeginModification  ();  //allow elements to be hidden
		/// After this function any link to deleted element will be replaced by InvalidHandle()
		/// TODO:
		///      1) maybe instead of forming set of deleted elements and subtracting set from other sets it is better
		///         to remove each modified element
		///         (done, check and compare)
		///      2) parent/child elements in set would not be replaced or reconnected, this may lead to wrong behavior
		///         (done, check and compare)
		void                              ApplyModification  ();  //modify DATA_REFERENCE tags so that links to hidden elements are converted to NULL and removed from sets
		/// This function is not yet implemented. It should correctly resolve parallel state of 
		/// newly created elements, provide them valid global identificators, resolve owners of
		/// the elements potentially optimized using information from BridgeTag and LayersTag
		/// May use ResolveShared function as basis but instead the whole mesh run the same algorithm for subset.
		void                              ResolveModification(); //resolve parallel state of newly created elements, restore ghost layers; not implemented, resuse ResolveShared code
		void                              EndModification    ();    //delete hidden elements
		enumerator                        getNext            (const HandleType * arr, enumerator size, enumerator k, MarkerType marker) const;
		enumerator                        Count              (const HandleType * arr, enumerator size, MarkerType marker) const;
		//implemented in mesh.cpp
	private:
		Tag                   tag_topologyerror;
		TopologyCheck         checkset;
		TopologyCheck         errorset;
	public:
		/// This function allows you to perform some topologycal checks before you create an element.
		/// Function is used internally by CreateEdge, CreateFace, CreateCell functions
		/// If you perform topological checks on your own, you'd probably better turn off checks before calling CreateXXX
		/// functions. Note that check for duplicates within mesh is performed by Mesh::FindSharedAdjacency.
		/// TODO: list checks performed inside in description
		TopologyCheck                     BeginTopologyCheck (ElementType etype, const HandleType * adj, enumerator num);
		/// This function performs some topologycal checks after creation of element.
		/// Function is used internally by CreateEdge, CreateFace, CreateCell functions.
		/// TODO: list checks performed inside in description.
		TopologyCheck                     EndTopologyCheck   (HandleType e); //check created element
		/// This will return tag by which you can retrive error mark to any element on which topogy check failed.
		/// As this is sparse tag you can check presence of error by Element::HaveData or Mesh::HaveData check.
		/// This tag will be valid only if you pass MARK_ON_ERROR to Mesh::GetTopologyCheck
		/// and will be deleted if you pass MARK_ON_ERROR to Mesh::RemTopologyCheck
		Tag                               TopologyErrorTag   () const {return tag_topologyerror;}
		/// Retrive currently set topology checks
		TopologyCheck                     GetTopologyCheck   (TopologyCheck mask = ENUMUNDEF) const {return checkset & mask;}
		/// Set topology checks
		void                              SetTopologyCheck   (TopologyCheck mask);
		/// Remove topology checks
		void                              RemTopologyCheck   (TopologyCheck mask);
		/// This will turn mesh into the state indicating that some topology error occured
		void                              SetTopologyError   (TopologyCheck mask) {errorset = errorset | mask;}
		/// Retrive topology error state, this indicates that some error have occured
		TopologyCheck                     GetTopologyError   (TopologyCheck mask = ENUMUNDEF) const {return errorset & mask;}
		/// Revert mesh to clean topology error state
		void                              ClearTopologyError (TopologyCheck mask = ENUMUNDEF) {errorset = errorset & ~mask;}
		//implemented in comparator.cpp
	public:
		class CentroidComparator
		{
			Mesh * m;
		public:
			CentroidComparator(Mesh * m) :m(m) {}
			CentroidComparator(const CentroidComparator & other) :m(other.m){}
			CentroidComparator & operator = (CentroidComparator const & other) { m = other.m; return *this;}
			int Compare(const real * a, const real * b);
			bool operator() (HandleType a, HandleType b);
			bool operator() (HandleType a, const real * b);
		};

		class GlobalIDComparator
		{
			Mesh * m;
		public:
			GlobalIDComparator(Mesh * m) :m(m) {}
			GlobalIDComparator(const GlobalIDComparator & other) :m(other.m){}
			GlobalIDComparator & operator = (GlobalIDComparator const & other) { m = other.m; return *this;}
			bool operator() (HandleType a, HandleType b) {if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->GlobalID(a) < m->GlobalID(b);}
			bool operator() (HandleType a, integer gid) {if( a == InvalidHandle() ) return false; return m->GlobalID(a) < gid;}
		};

		class IerarhyComparator
		{
			Mesh * m;
		public:
			IerarhyComparator(Mesh * m) :m(m) {}
			IerarhyComparator(const IerarhyComparator & other) :m(other.m){}
			IerarhyComparator & operator = (IerarhyComparator const & other) { m = other.m; return *this;}
			int CompareNodes(HandleType a, HandleType b);
			int CompareElements(HandleType a, HandleType b);
			bool operator() (HandleType a, HandleType b) {if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return CompareElements(a,b) < 0;}
		};

		class RealComparator
		{
			Mesh * m; Tag t;
		public:
			RealComparator(Mesh * m, Tag t) :m(m), t(t) {}
			RealComparator(const RealComparator & other) :m(other.m), t(other.t){}
			RealComparator & operator = (RealComparator const & other) { m = other.m; t = other.t; return *this;}
			bool operator() (HandleType a, HandleType b){if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->Real(a,t) < m->Real(b,t);}
			bool operator() (HandleType a, real b) {if( a == InvalidHandle() ) return true; return m->Real(a,t) < b;}
		};

		class IntegerComparator
		{
			Mesh * m; Tag t;
		public:
			IntegerComparator(Mesh * m, Tag t) :m(m), t(t) {}
			IntegerComparator(const IntegerComparator & other) :m(other.m), t(other.t){}
			IntegerComparator & operator = (IntegerComparator const & other) { m = other.m; t = other.t; return *this;}
			bool operator() (HandleType a, HandleType b){if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->Integer(a,t) < m->Integer(b,t);}
			bool operator() (HandleType a, integer b) {if( a == InvalidHandle() ) return true; return m->Integer(a,t) < b;}
		};

		class BulkComparator
		{
			Mesh * m; Tag t;
		public:
			BulkComparator(Mesh * m, Tag t) :m(m), t(t) {}
			BulkComparator(const BulkComparator & other) :m(other.m), t(other.t){}
			BulkComparator & operator = (BulkComparator const & other) { m = other.m; t = other.t; return *this;}
			bool operator() (HandleType a, HandleType b){if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->Bulk(a,t) < m->Bulk(b,t);}
			bool operator() (HandleType a, bulk b) {if( a == InvalidHandle() ) return true; return m->Bulk(a,t) < b;}
		};
		
		class RealDFComparator
		{
			Mesh * m; Tag t;
		public:
			RealDFComparator(Mesh * m, Tag t) :m(m), t(t) {}
			RealDFComparator(const RealDFComparator & other) :m(other.m), t(other.t){}
			RealDFComparator & operator = (RealDFComparator const & other) { m = other.m; t = other.t; return *this;}
			bool operator() (HandleType a, HandleType b){if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->RealDF(a,t) < m->RealDF(b,t);}
			bool operator() (HandleType a, real b) {if( a == InvalidHandle() ) return true; return m->RealDF(a,t) < b;}
		};

		class IntegerDFComparator
		{
			Mesh * m; Tag t;
		public:
			IntegerDFComparator(Mesh * m, Tag t) :m(m), t(t) {}
			IntegerDFComparator(const IntegerDFComparator & other) :m(other.m), t(other.t){}
			IntegerDFComparator & operator = (IntegerDFComparator const & other) { m = other.m; t = other.t; return *this;}
			bool operator() (HandleType a, HandleType b){if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->IntegerDF(a,t) < m->IntegerDF(b,t);}
			bool operator() (HandleType a, integer b) {if( a == InvalidHandle() ) return true; return m->IntegerDF(a,t) < b;}
		};

		class BulkDFComparator
		{
			Mesh * m; Tag t;
		public:
			BulkDFComparator(Mesh * m, Tag t) :m(m), t(t) {}
			BulkDFComparator(const BulkDFComparator & other) :m(other.m), t(other.t){}
			BulkDFComparator & operator = (BulkDFComparator const & other) { m = other.m; t = other.t; return *this;}
			bool operator() (HandleType a, HandleType b){if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->BulkDF(a,t) < m->BulkDF(b,t);}
			bool operator() (HandleType a, bulk b) {if( a == InvalidHandle() ) return true; return m->BulkDF(a,t) < b;}
		};

		void SortHandles(HandleType * h, enumerator num);
		/// TODO 53 check that putting global ids to array will be faster
		void SortByGlobalID(HandleType * h, enumerator num);
	};
}

#endif

#endif // INMOST_MESH_H_INCLUDED
