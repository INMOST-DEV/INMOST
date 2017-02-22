#ifndef INMOST_DATA_H_INCLUDED
#define INMOST_DATA_H_INCLUDED

#include "inmost_common.h"
#if defined(USE_AUTODIFF)
#include "inmost_expression.h"
#endif


#if defined(USE_MESH)


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

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// ElementType
  
	typedef INMOST_DATA_BULK_TYPE ElementType;
	static const ElementType NONE = 0x00;
	static const ElementType NODE = 0x01;
	static const ElementType EDGE = 0x02;
	static const ElementType FACE = 0x04;
	static const ElementType CELL = 0x08;
	static const ElementType ESET = 0x10;
	static const ElementType MESH = 0x20;
#define NUM_ELEMENT_TYPS 6

	__INLINE bool                         OneType             (ElementType t) {return t > 0 && (t & (t-1)) == 0;}
	__INLINE ElementType                  FirstElementType    () {return NODE;}
	__INLINE ElementType                  LastElementType     () {return MESH << 1;}
	__INLINE ElementType                  NextElementType     (ElementType etype) {return etype << 1;}
	__INLINE ElementType                  PrevElementType     (ElementType etype) {return etype >> 1;}
	__INLINE ElementType                  ElementTypeFromDim  (INMOST_DATA_INTEGER_TYPE dim) {return 1 << dim;}
	const char *                          ElementTypeName     (ElementType t); //mesh.cpp
	__INLINE INMOST_DATA_INTEGER_TYPE     ElementNum          (ElementType t)
	{
		unsigned int v = static_cast<unsigned int>(t);  // 32-bit value to find the log2 of 
		//static const unsigned int b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 0xFF00FF00, 0xFFFF0000};
		static const unsigned int b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0};
		unsigned int r = (v & b[0]) != 0;
		//r |= ((v & b[4]) != 0) << 4;
		//r |= ((v & b[3]) != 0) << 3;
		r |= ((v & b[2]) != 0) << 2;
		r |= ((v & b[1]) != 0) << 1;
		return static_cast<INMOST_DATA_INTEGER_TYPE>(r);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// MarkerType
  
	/// Low 8 bits - marker mask, rest high bits - position of marker.
	typedef INMOST_DATA_ENUM_TYPE         MarkerType; 
	/// Number of chars to hold all markers, total number (MarkerFields * bits_per_char).
	static const INMOST_DATA_ENUM_TYPE    MarkerFields        = 16;   
	/// Number of chars to hold all private markers, total number (MarkerFields * bits_per_char).
	static const INMOST_DATA_ENUM_TYPE    MarkerFieldsPrivate = 4;
	/// Last bit indicate whether the marker is private.
	static const INMOST_DATA_ENUM_TYPE    MarkerPrivateBit    = 1 << (sizeof(INMOST_DATA_ENUM_TYPE)*8-1); 
	/// Bit mask to obtain marker mask within MarkerType.
	static const INMOST_DATA_ENUM_TYPE    MarkerMask          = static_cast<INMOST_DATA_BULK_TYPE>(-1); 
	/// sizeof(char) * bits_per_char.
	static const INMOST_DATA_ENUM_TYPE    MarkerShift         = sizeof(INMOST_DATA_BULK_TYPE)*8;    
	__INLINE static bool                  isPrivate           (MarkerType n) {return (n & MarkerPrivateBit) == MarkerPrivateBit;}
	__INLINE static MarkerType            InvalidMarker       (){return ENUMUNDEF;}


	//////////////////////////////////////////////////////////////////////////////////////////////////
	// HandleType

	typedef INMOST_DATA_ENUM_TYPE         HandleType;
	static const INMOST_DATA_ENUM_TYPE    handle_etype_bits   = 3;
	static const INMOST_DATA_ENUM_TYPE    handle_etype_shift  = sizeof(HandleType)*8-handle_etype_bits;
	static const INMOST_DATA_ENUM_TYPE    handle_id_mask      = (1 << handle_etype_shift)-1;

	static const INMOST_DATA_ENUM_TYPE    chunk_bits_elems    = 13;
	static const INMOST_DATA_ENUM_TYPE    chunk_bits_empty    = 8;
	static const INMOST_DATA_ENUM_TYPE    chunk_bits_tags     = 6;
	static const INMOST_DATA_ENUM_TYPE    chunk_bits_dense    = 6;

	__INLINE HandleType                   InvalidHandle       () {return 0;}
	__INLINE INMOST_DATA_INTEGER_TYPE     GetHandleID         (HandleType h) {return (h & handle_id_mask)-1;}
	__INLINE INMOST_DATA_INTEGER_TYPE     GetHandleElementNum (HandleType h) {return h >> handle_etype_shift;}
	__INLINE ElementType                  GetHandleElementType(HandleType h) {return 1 << GetHandleElementNum(h);}
	__INLINE HandleType                   ComposeHandle       (ElementType etype, INMOST_DATA_INTEGER_TYPE ID) {return ID == -1 ? InvalidHandle() : ((ElementNum(etype) << handle_etype_shift) + (1+ID));}
	__INLINE HandleType                   ComposeCellHandle   (INMOST_DATA_INTEGER_TYPE ID) {return ID == -1 ? InvalidHandle() : ((ElementNum(CELL) << handle_etype_shift) + (1+ID));}
	__INLINE HandleType                   ComposeFaceHandle   (INMOST_DATA_INTEGER_TYPE ID) {return ID == -1 ? InvalidHandle() : ((ElementNum(FACE) << handle_etype_shift) + (1+ID));}
	__INLINE HandleType                   ComposeEdgeHandle   (INMOST_DATA_INTEGER_TYPE ID) {return ID == -1 ? InvalidHandle() : ((ElementNum(EDGE) << handle_etype_shift) + (1+ID));}
	__INLINE HandleType                   ComposeNodeHandle   (INMOST_DATA_INTEGER_TYPE ID) {return ID == -1 ? InvalidHandle() : ((ElementNum(NODE) << handle_etype_shift) + (1+ID));}
	__INLINE HandleType                   ComposeSetHandle    (INMOST_DATA_INTEGER_TYPE ID) {return ID == -1 ? InvalidHandle() : ((ElementNum(ESET) << handle_etype_shift) + (1+ID));}
	__INLINE HandleType                   ComposeHandleNum    (INMOST_DATA_INTEGER_TYPE etypenum, INMOST_DATA_INTEGER_TYPE ID) {return ID == -1 ? InvalidHandle() : ((etypenum << handle_etype_shift) + (1+ID));}
	__INLINE bool                         isValidHandle       (HandleType h) {return h != 0;}

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// RemoteElementType

	typedef std::pair<Mesh*,HandleType>   RemoteHandleType;
	/// Construct an object of type Element, hande cannot be modified.
	Element                               MakeElement         (const RemoteHandleType & rh); //storage.cpp
	/// Construct an object of type Element, hande can be modified.
	Element                               MakeElementRef      (RemoteHandleType & rh); //storage.cpp
	


#if defined(USE_AUTODIFF)
	typedef array<variable>               inner_variable_array;
#endif
	typedef array<INMOST_DATA_REAL_TYPE>  inner_real_array;
	typedef array<INMOST_DATA_INTEGER_TYPE> 
		                                  inner_integer_array;
	typedef array<INMOST_DATA_BULK_TYPE>  inner_bulk_array;
	typedef array<HandleType>             inner_reference_array;
	typedef array<RemoteHandleType>       inner_remote_reference_array;

	enum DataType
	{
		DATA_REAL                           = 0, 
		DATA_INTEGER                        = 1, 
		DATA_BULK                           = 2,
		DATA_REFERENCE                      = 3,
		DATA_REMOTE_REFERENCE               = 4,
#if defined(USE_AUTODIFF)
		DATA_VARIABLE                       = 5
#endif
	};

	///Returns a name of the data type as a string.
	const char *                          DataTypeName        (DataType t);


	///This class is a data container for class Tag, contains all the necessery information to access
	/// mesh data in class TagManager. It is never exposed to the user directly, should be accessed 
	/// through interface class Tag. Functions are implemented in tag.cpp.
	class TagMemory 
	{
	public:
		///Destructor should not do anything.
		~TagMemory() {}
		///Copy constructor, copies all the data except for m_link. Main purpose is to create an exact 
		/// copy for different mesh, whenever another mesh is created.
		TagMemory(Mesh * m, const TagMemory & other);
		///Assignment operator should not be ever used, but is here for convinience.
		TagMemory & operator =(TagMemory const & other);
	private:
		///Common constructor shouldn't be called from outside.
		TagMemory();
		///Position of data in memory for each type of element.
		INMOST_DATA_ENUM_TYPE pos[NUM_ELEMENT_TYPS]; 
		///Type of represented data.
		DataType dtype; 
		///Specified name for the data.
		std::string tagname; 
		///Associated MPI-type for communication.
		INMOST_MPI_Type bulk_data_type; 
		///Number of entries in each record of the type.
		///May be set to ENUMUNDEF to represent variable-sized arrays.
		INMOST_DATA_ENUM_TYPE size;
		///Number of bytes used to represent data type in memory.
		INMOST_DATA_ENUM_TYPE bytes_size;
		///Indicates whether the data is represented as sparse
		/// on certain elements of the mesh.
		bool sparse[NUM_ELEMENT_TYPS];
		///Number of bytes used to store data for one element. It is size times bytes_size for data of 
		// fixed size or number of bytes for the structure used to represent data of variable size.
		INMOST_DATA_ENUM_TYPE record_size;
		///Link to the mesh.
		Mesh * m_link;
		/// Provide access to interface.
		friend class Tag;
		/// For debug purposes only.
		friend class Storage;
	};

  
  ///This class provides the access to the individual mesh datum and general information about it. 
	class Tag //implemented in tag.cpp
	{
	private:
		///A link to the data for the current instance of the class Tag.
		TagMemory * mem;
		///The use of this privite constructor is reserved for TagManager::CreateTag function.
		Tag (Mesh * m,std::string name, DataType _dtype, INMOST_DATA_ENUM_TYPE size); 
		__INLINE INMOST_DATA_ENUM_TYPE GetRecordSize() const;
		__INLINE void SetSize(INMOST_DATA_ENUM_TYPE size);
		__INLINE void SetPosition(INMOST_DATA_ENUM_TYPE pos, ElementType type);
		__INLINE INMOST_DATA_ENUM_TYPE GetPosition(ElementType type) const;
		__INLINE void SetSparse(ElementType type);
		__INLINE INMOST_DATA_ENUM_TYPE GetPositionByDim(INMOST_DATA_ENUM_TYPE typenum) const;
	public:
		~Tag();
		Tag();
		Tag(const Tag & other);
		__INLINE bool operator <(const Tag & other) const;
		__INLINE bool operator >(const Tag & other) const;
		__INLINE bool operator ==(const Tag & other) const;
		__INLINE bool operator !=(const Tag & other) const;
		__INLINE Tag & operator =(Tag const & other);
		__INLINE DataType GetDataType() const;
		__INLINE INMOST_MPI_Type GetBulkDataType() const;
		__INLINE INMOST_DATA_ENUM_TYPE GetBytesSize() const;
		__INLINE INMOST_DATA_ENUM_TYPE GetSize() const;
		__INLINE std::string GetTagName() const;
		__INLINE bool isDefined(ElementType type) const;
		__INLINE bool isSparse(ElementType type) const;
		__INLINE bool isValid() const;
		__INLINE Mesh * GetMeshLink() const;
		__INLINE bool isSparseByDim(INMOST_DATA_INTEGER_TYPE typenum)const;
		__INLINE bool isDefinedByDim(INMOST_DATA_INTEGER_TYPE typenum)const;
		__INLINE void SetBulkDataType(INMOST_MPI_Type type);
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
		typedef chunk_array<dense_sub_type,chunk_bits_dense>           dense_data_array_type;
		typedef struct{void * tag, * rec;}                             sparse_sub_record;
		typedef array< sparse_sub_record >                             sparse_sub_type;
		typedef chunk_array< sparse_sub_type,chunk_bits_elems>         sparse_data_array_type;
		typedef chunk_array<INMOST_DATA_INTEGER_TYPE,chunk_bits_elems> back_links_type;
		//typedef std::vector<INMOST_DATA_INTEGER_TYPE>                  back_links_type;
	public:
		typedef tag_array_type::iterator iteratorTag;
	public:
		virtual ~TagManager();
		/// Check existance of a data tag by it's name.
		bool HaveTag(std::string name) const;
		/// Retrive a data tag by it's name.
		Tag GetTag(std::string name) const;
		/// Retrive names for all the tags present on the mesh.
		void ListTagNames(std::vector<std::string> & list) const;
		/// Create tag with prescribed attributes.
		Tag CreateTag(Mesh * m, std::string name, DataType dtype, ElementType etype, ElementType sparse, INMOST_DATA_ENUM_TYPE size = ENUMUNDEF); 
		/// Delete tag from certain elements.
		virtual Tag DeleteTag(Tag tag, ElementType mask); 
		/// Check that the tag was defined on certain elements.
		bool ElementDefined(Tag const & tag, ElementType etype) const;
	protected:
		/// Shrink or enlarge arrays for a dense data.
		void ReallocateData(const Tag & t, INMOST_DATA_INTEGER_TYPE etypenum,INMOST_DATA_ENUM_TYPE new_size);
		/// Reallocate all the data for all the tags.
		void ReallocateData(INMOST_DATA_INTEGER_TYPE etypenum, INMOST_DATA_ENUM_TYPE new_size);
		///Retrive substructure for representation of the sparse data without permission for modification.
		__INLINE sparse_sub_type const & GetSparseData(int etypenum, int local_id) const {return sparse_data[etypenum][local_id];}
		///Retrive substructure for representation of the sparse data.
		__INLINE sparse_sub_type & GetSparseData(int etypenum, int local_id) {return sparse_data[etypenum][local_id];}
		///Retrive substructure for representation of the dense data without permission for modification.
		__INLINE dense_sub_type const & GetDenseData(int pos) const {return dense_data[pos];}
		///Retrive substructure for representation of the dense data.
		__INLINE dense_sub_type & GetDenseData(int pos) {return dense_data[pos];}
		///Copy data from one element to another.
		static void CopyData(const Tag & t, void * adata, const void * bdata);
		///Destroy data that represents array of variable size.
		static void DestroyVariableData(const Tag & t, void * adata);
	protected:
		typedef tag_array_type::iterator       tag_iterator; //< Use this type to iterate over tags of the mesh.
		typedef tag_array_type::const_iterator tag_const_iterator; //< Use this type to iterate over tags of the mesh without right for modification.
	protected:
		tag_array_type         tags;
		empty_data             empty_dense_data;
		dense_data_array_type  dense_data;
		sparse_data_array_type sparse_data[NUM_ELEMENT_TYPS];
		back_links_type        back_links[NUM_ELEMENT_TYPS];
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
		typedef INMOST_DATA_REAL_TYPE       real; 
		/// Storage type for representing integer values.
		typedef INMOST_DATA_INTEGER_TYPE    integer;
		/// Storage type for representing one byte of abstact data.
		typedef INMOST_DATA_BULK_TYPE       bulk;
		/// type for representing unsigned integer values.
		typedef INMOST_DATA_ENUM_TYPE       enumerator;
		/// Storage type for representing references to Element.
		typedef HandleType                  reference;
		/// Storage type for representing references to Element in another Mesh.
		typedef RemoteHandleType            remote_reference;
		/// Storage type for representing arrays of real values.
		typedef shell<real>                 real_array;
		/// Storage type for representing arrays of integer values.
		typedef shell<integer>              integer_array;
		/// Storage type for representing abstact data as a series of bytes.
		typedef shell<bulk>                 bulk_array;
#if defined(USE_AUTODIFF)
		/// Storage type for representing real value with vector of variations
		typedef variable                    var;
		/// Storage type for representing array of values with vectors of variations
		typedef shell<variable>             var_array;
#endif
		/// Storage type for representing arrays of Element references.
		class reference_array : public shell<reference>
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
				iterator() :shell<HandleType>::iterator() {}
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
				const_iterator() :shell<HandleType>::const_iterator() {}
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
				reverse_iterator() :shell<HandleType>::reverse_iterator() {}
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
				const_reverse_iterator() :shell<HandleType>::const_reverse_iterator() {}
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
			void SetMeshLink(Mesh * new_m) {m = new_m;}
			Element back();
			Element back() const;
			Element front();
			Element front() const;
		};
		/// Storage type for representing arrays of Element references on another Mesh.
		class remote_reference_array : public shell<remote_reference>
		{
		public:
			remote_reference_array() :shell<remote_reference>() {}
			remote_reference_array(const shell<remote_reference> & other) :shell<remote_reference>(other) {}
			remote_reference_array(const remote_reference_array & other) :shell<remote_reference>(other) {}
			remote_reference_array(remote_reference * pntr, shell<remote_reference>::size_type psize) : shell<remote_reference>(pntr,psize) {}
			void push_back(const Storage & elem);
			void push_back(Mesh * m, HandleType h) {shell<remote_reference>::push_back(RemoteHandleType(m,h));} //is it needed?
			Element operator[] (size_type n);
			Element operator[] (size_type n) const;
			class iterator : public shell<remote_reference>::iterator 
			{
			public:
				iterator() :shell<remote_reference>::iterator() {}
				iterator(const shell<remote_reference>::iterator & other) : shell<remote_reference>::iterator(other) {}
				iterator(const iterator & other) : shell<remote_reference>::iterator(other) {}
				Element operator->();
			};
			class const_iterator : public shell<remote_reference>::const_iterator 
			{
			public:
				const_iterator() :shell<remote_reference>::const_iterator() {}
				const_iterator(const shell<remote_reference>::const_iterator & other) : shell<remote_reference>::const_iterator(other) {}
				const_iterator(const const_iterator & other) : shell<remote_reference>::const_iterator(other) {}
				Element operator->();
			};
			class reverse_iterator : public shell<remote_reference>::reverse_iterator 
			{
			public:
				reverse_iterator() :shell<remote_reference>::reverse_iterator() {}
				reverse_iterator(const shell<remote_reference>::reverse_iterator & other) : shell<remote_reference>::reverse_iterator(other) {}
				reverse_iterator(const reverse_iterator & other) : shell<remote_reference>::reverse_iterator(other) {}
				Element operator->();
			};
			class const_reverse_iterator : public shell<remote_reference>::const_reverse_iterator 
			{
			public:
				const_reverse_iterator() :shell<remote_reference>::const_reverse_iterator() {}
				const_reverse_iterator(const shell<remote_reference>::const_reverse_iterator & other) : shell<remote_reference>::const_reverse_iterator(other) {}
				const_reverse_iterator(const const_reverse_iterator & other) : shell<remote_reference>::const_reverse_iterator(other) {}
				Element operator->();
			};
			Element back();
			Element back() const;
			Element front();
			Element front() const;
		};
		//typedef shell<reference>          reference_array;
	protected:
		HandleType                          handle;
		HandleType *                        handle_link;
	private:
		Mesh *                              m_link;
	public:
		Storage(const Storage & other) : handle(other.handle), handle_link(other.handle_link), m_link(other.m_link) {}
		Storage(Mesh * mesh, HandleType handle) : handle(handle), handle_link(NULL), m_link(mesh) {}
		/// This constructor allows for remote handle modification
		Storage(Mesh * mesh, HandleType * handle) : handle(*handle), handle_link(handle), m_link(mesh) {}
		/// If there is a link to handle provided (automatically by ElementArray and reference_array),
		/// then remote handle value will be modified
		Storage &                           operator =          (Storage const & other); 
		__INLINE bool                       operator <          (const Storage & other) const {return handle < other.handle;}
		__INLINE bool                       operator >          (const Storage & other) const {return handle > other.handle;}
		__INLINE bool                       operator <=         (const Storage & other) const {return handle <= other.handle;}
		__INLINE bool                       operator >=         (const Storage & other) const {return handle >= other.handle;}
		__INLINE bool                       operator ==         (const Storage & other) const {return handle == other.handle;}
		__INLINE bool                       operator !=         (const Storage & other) const {return handle != other.handle;}
		__INLINE Storage *                  operator->          () {return this;}
		__INLINE const Storage *            operator->          () const {return this;}
		__INLINE Storage &                  self                () {return *this;}
		__INLINE const Storage &            self                () const {return *this;}
		virtual ~Storage() {}
	public:
		/// Retrieve real value associated with Tag. Implemented in inmost_mesh.h.
		/// @param tag instance of class Tag that represent asking datum. 
		__INLINE real &                     Real                (const Tag & tag) const;
		/// Retrieve integer value associated with Tag. Implemented in inmost_mesh.h.
		/// @param tag instance of class Tag that represent asking datum. 
		__INLINE integer &                  Integer             (const Tag & tag) const;
		/// Retrieve one byte of abstract data associated with Tag. Implemented in inmost_mesh.h.
		/// @param tag instance of class Tag that represent asking datum. 
		__INLINE bulk &                     Bulk                (const Tag & tag) const;
		/// Retrieve Element reference associated with Tag. Implemented in inmost_mesh.h.
		/// @param tag instance of class Tag that represent asking datum. 
		__INLINE reference &                Reference           (const Tag & tag) const;
		/// Retrieve remote Element reference associated with Tag. Implemented in inmost_mesh.h.
		/// @param tag instance of class Tag that represent asking datum. 
		__INLINE remote_reference &         RemoteReference     (const Tag & tag) const;
		/// Retrieve array of real values associated with Tag. Implemented in inmost_mesh.h.
		/// @param tag instance of class Tag that represent asking datum. 
		__INLINE real_array                 RealArray           (const Tag & tag) const;
		/// Retrieve array of integer values associated with Tag. Implemented in inmost_mesh.h.
		/// @param tag instance of class Tag that represent asking datum. 
		__INLINE integer_array              IntegerArray        (const Tag & tag) const;
		/// Retrieve abstract data associated with Tag as a series of bytes. 
		/// Implemented in inmost_mesh.h.
		/// @param tag instance of class Tag that represent asking datum. 
		__INLINE bulk_array                 BulkArray           (const Tag & tag) const;
		/// Retrieve array of Element references associated with Tag. 
		/// Implemented in inmost_mesh.h.
		/// @param tag instance of class Tag that represent asking datum. 
		__INLINE reference_array            ReferenceArray      (const Tag & tag) const;
		/// Retrieve array of Element references associated with Tag. 
		/// Implemented in inmost_mesh.h.
		/// @param tag instance of class Tag that represent asking datum. 
		__INLINE remote_reference_array     RemoteReferenceArray(const Tag & tag) const;
    
		//optimized data requests for dense data with fixed size
		__INLINE real_array                 RealArrayDF         (const Tag & tag) const;
		__INLINE integer_array              IntegerArrayDF      (const Tag & tag) const;
		__INLINE bulk_array                 BulkArrayDF         (const Tag & tag) const;
		__INLINE reference_array            ReferenceArrayDF    (const Tag & tag) const;
		__INLINE remote_reference_array     RemoteReferenceArrayDF(const Tag & tag) const;
		__INLINE real      &                RealDF              (const Tag & tag) const;
		__INLINE integer   &                IntegerDF           (const Tag & tag) const;
		__INLINE bulk      &                BulkDF              (const Tag & tag) const;
		__INLINE reference &                ReferenceDF         (const Tag & tag) const;
		__INLINE remote_reference &         RemoteReferenceDF   (const Tag & tag) const;
		
		//optimized data requests for dense data with variable size
		__INLINE real_array                 RealArrayDV         (const Tag & tag) const;
		__INLINE integer_array              IntegerArrayDV      (const Tag & tag) const;
		__INLINE bulk_array                 BulkArrayDV         (const Tag & tag) const;
		__INLINE reference_array            ReferenceArrayDV    (const Tag & tag) const;
		__INLINE remote_reference_array     RemoteReferenceArrayDV(const Tag & tag) const;
		__INLINE real      &                RealDV              (const Tag & tag) const;
		__INLINE integer   &                IntegerDV           (const Tag & tag) const;
		__INLINE bulk      &                BulkDV              (const Tag & tag) const;
		__INLINE reference &                ReferenceDV         (const Tag & tag) const;
		__INLINE remote_reference &         RemoteReferenceDV   (const Tag & tag) const;
#if defined(USE_AUTODIFF)
		/// Retrieve variable reference associated with Tag.
		__INLINE var &                      Variable            (const Tag & tag) const;
		__INLINE var &                      VariableDF          (const Tag & tag) const;
		__INLINE var &                      VariableDV          (const Tag & tag) const;
		/// Retrieve array of variables associated with Tag.
		__INLINE var_array                  VariableArray       (const Tag & tag) const;
		__INLINE var_array                  VariableArrayDF     (const Tag & tag) const;
		__INLINE var_array                  VariableArrayDV     (const Tag & tag) const;
#endif
		
		/// Return the data length associated with Tag.
		/// For abstract data return the number of bytes, otherwise return the length of associated array. 
		/// @param tag tag that represents the data
		/// @see Storage::SetDataSize
		/// @see Mesh::GetDataSize
		__INLINE INMOST_DATA_ENUM_TYPE      GetDataSize         (const Tag & tag) const;
		/// Return the size of the structure required to represent the data on current element.
		/// This is equal to GetDataSize times Tag::GetBytesSize for all the data types,
		/// except for DATA_VARIABLE, that requires a larger structure to accomodate derivatives.
		/// @param tag tag that represents the data
		/// @see Mesh::GetDataCapacity
		__INLINE INMOST_DATA_ENUM_TYPE      GetDataCapacity     (const Tag & tag) const;
		/// Set the length of  data associated with Tag.
		/// @param tag Identifying Tag.
		/// @param new_size The number of bytes for abstract data, otherwise the length of the array.
		/// @see Storage::GetDataSize
		/// @see Mesh::SetDataSize
		__INLINE void                       SetDataSize         (const Tag & tag,
                                                             INMOST_DATA_ENUM_TYPE new_size) const;
		/// Extract part of the data associated with Tag.
		/// Copy part of the associated array or data to the destination memory.
		/// @param tag Identifying Tag.
		/// @param shift Starting position of the copied data.
        /// For abstact data - number of bytes to skip, otherwise number of values to skip.
		/// @param size Number of elements to copy.
        /// For abstact data - number of bytes to copy, otherwise number of values to copy.
		/// @param data Destination position to copy data to.
		/// @see Storage::SetData
		__INLINE void                       GetData             (const Tag & tag, 
                                                             INMOST_DATA_ENUM_TYPE shift, 
                                                             INMOST_DATA_ENUM_TYPE size, 
                                                             void * data) const;
		__INLINE void                       SetData             (const Tag & tag, 
                                                             INMOST_DATA_ENUM_TYPE shift, 
                                                             INMOST_DATA_ENUM_TYPE size, 
                                                             const void * data) const;
		__INLINE void                       DelData             (const Tag & tag) const;
		/// Deallocates space allocated for sparse data, frees variable array if necessary
		__INLINE void                       DelSparseData       (const Tag & tag) const;
		/// Frees variable array or fills field with zeroes
		__INLINE void                       DelDenseData        (const Tag & tag) const;
		/// Check if any data is associated with Tag.
		__INLINE bool                       HaveData            (const Tag & tag) const;
		__INLINE ElementType                GetElementType      () const;
		__INLINE integer                    GetElementNum       () const;
		__INLINE void                       SetMarker           (MarkerType n) const;
		__INLINE bool                       GetMarker           (MarkerType n) const;
		__INLINE void                       RemMarker           (MarkerType n) const;
		__INLINE void                       SetPrivateMarker    (MarkerType n) const;
		__INLINE bool                       GetPrivateMarker    (MarkerType n) const;
		__INLINE void                       RemPrivateMarker    (MarkerType n) const;
		__INLINE void                       ClearMarkerSpace    () const;
		__INLINE void                       GetMarkerSpace      (bulk copy[MarkerFields]) const;
		__INLINE void                       SetMarkerSpace      (bulk source[MarkerFields]) const;
		__INLINE integer                    LocalID             () const;
		/// This number is guaranteed to be between 0 and Mesh::NumberOf(type of element)
		/// after Mesh::ReorderEmpty
		__INLINE integer                    DataLocalID         () const;
		__INLINE bool                       isValid             () const;
		__INLINE Mesh *                     GetMeshLink         () const;
		__INLINE HandleType                 GetHandle           () const;
	    __INLINE Element                    getAsElement        () const;
		__INLINE Node                       getAsNode           () const;
		__INLINE Edge                       getAsEdge           () const;
		__INLINE Face                       getAsFace           () const;
		__INLINE Cell                       getAsCell           () const;
		__INLINE ElementSet                 getAsSet            () const;
		friend class Mesh;
	};
	
	//////////////////////////////////////////////////////////////////////
	/// Helper classes for class Tag                                    //
	//////////////////////////////////////////////////////////////////////
	
	class TagReal : public Tag
	{
	public:
		TagReal() : Tag() {}
		TagReal(const TagReal & b) : Tag(b) {}
		TagReal(const Tag & b) : Tag(b) {}
		TagReal & operator = (TagReal const & b) {Tag::operator =(b); return *this;}
		TagReal & operator = (Tag const & b) {Tag::operator =(b); return *this;}
		Storage::real & operator [](const Storage & arg) const {return arg.Real(*static_cast<const Tag*>(this));}
	};
	
	class TagInteger : public Tag
	{
	public:
		TagInteger() : Tag() {}
		TagInteger(const TagInteger & b) : Tag(b) {}
		TagInteger(const Tag & b) : Tag(b) {}
		TagInteger & operator = (TagInteger const & b) {Tag::operator =(b); return *this;}
		TagInteger & operator = (Tag const & b) {Tag::operator =(b); return *this;}
		Storage::integer & operator [](const Storage & arg) const {return arg.Integer(*static_cast<const Tag*>(this));}
	};
	
	class TagBulk : public Tag
	{
	public:
		TagBulk() : Tag() {}
		TagBulk(const TagBulk & b) : Tag(b) {}
		TagBulk(const Tag & b) : Tag(b) {}
		TagBulk & operator = (TagBulk const & b) {Tag::operator =(b); return *this;}
		TagBulk & operator = (Tag const & b) {Tag::operator =(b); return *this;}
		Storage::bulk & operator [](const Storage & arg) const {return arg.Bulk(*static_cast<const Tag*>(this));}
	};
	
	class TagReference : public Tag
	{
	public:
		TagReference() : Tag() {}
		TagReference(const TagReference & b) : Tag(b) {}
		TagReference(const Tag & b) : Tag(b) {}
		TagReference & operator = (TagReference const & b) {Tag::operator =(b); return *this;}
		TagReference & operator = (Tag const & b) {Tag::operator =(b); return *this;}
		Storage::reference & operator [](const Storage & arg) const {return arg.Reference(*static_cast<const Tag*>(this));}
	};
	
	
	class TagRealArray : public Tag
	{
	public:
		TagRealArray() : Tag() {}
		TagRealArray(const TagRealArray & b) : Tag(b) {}
		TagRealArray(const Tag & b) : Tag(b) {}
		TagRealArray & operator = (TagRealArray const & b) {Tag::operator =(b); return *this;}
		TagRealArray & operator = (Tag const & b) {Tag::operator =(b); return *this;}
		Storage::real_array operator [](const Storage & arg) const {return arg.RealArray(*static_cast<const Tag*>(this));}
	};
	
	class TagIntegerArray : public Tag
	{
	public:
		TagIntegerArray() : Tag() {}
		TagIntegerArray(const TagIntegerArray & b) : Tag(b) {}
		TagIntegerArray(const Tag & b) : Tag(b) {}
		TagIntegerArray & operator = (TagIntegerArray const & b) {Tag::operator =(b); return *this;}
		TagIntegerArray & operator = (Tag const & b) {Tag::operator =(b); return *this;}
		Storage::integer_array operator [](const Storage & arg) const {return arg.IntegerArray(*static_cast<const Tag*>(this));}
	};
	
	class TagBulkArray : public Tag
	{
	public:
		TagBulkArray() : Tag() {}
		TagBulkArray(const TagBulkArray & b) : Tag(b) {}
		TagBulkArray(const Tag & b) : Tag(b) {}
		TagBulkArray & operator = (TagBulkArray const & b) {Tag::operator =(b); return *this;}
		TagBulkArray & operator = (Tag const & b) {Tag::operator =(b); return *this;}
		Storage::bulk_array operator [](const Storage & arg) const {return arg.BulkArray(*static_cast<const Tag*>(this));}
	};
	
	class TagReferenceArray : public Tag
	{
	public:
		TagReferenceArray() : Tag() {}
		TagReferenceArray(const TagReferenceArray & b) : Tag(b) {}
		TagReferenceArray(const Tag & b) : Tag(b) {}
		TagReferenceArray & operator = (TagReferenceArray const & b) {Tag::operator =(b); return *this;}
		TagReferenceArray & operator = (Tag const & b) {Tag::operator =(b); return *this;}
		Storage::reference_array operator [](const Storage & arg) const {return arg.ReferenceArray(*static_cast<const Tag*>(this));}
	};
	
#if defined(USE_AUTODIFF)
	class TagVariable : public Tag
	{
	public:
		TagVariable() : Tag() {}
		TagVariable(const TagVariable & b) : Tag(b) {}
		TagVariable(const Tag & b) : Tag(b) {}
		TagVariable & operator = (TagVariable const & b) {Tag::operator =(b); return *this;}
		TagVariable & operator = (Tag const & b) {Tag::operator =(b); return *this;}
		Storage::var & operator [](const Storage & arg) const {return arg.Variable(*static_cast<const Tag*>(this));}
	};
	
	class TagVariableArray : public Tag
	{
	public:
		TagVariableArray() : Tag() {}
		TagVariableArray(const TagVariableArray & b) : Tag(b) {}
		TagVariableArray(const Tag & b) : Tag(b) {}
		TagVariableArray & operator = (TagVariableArray const & b) {Tag::operator =(b); return *this;}
		TagVariableArray & operator = (Tag const & b) {Tag::operator =(b); return *this;}
		Storage::var_array operator [](const Storage & arg) const {return arg.VariableArray(*static_cast<const Tag*>(this));}
	};
#endif //USE_AUTODIFF

	//////////////////////////////////////////////////////////////////////
	/// Inline functions for class Tag                                  //
	//////////////////////////////////////////////////////////////////////

	__INLINE INMOST_DATA_ENUM_TYPE  Tag::GetRecordSize() const 
	{
		return mem->record_size;
	}

	__INLINE void Tag::SetSize(INMOST_DATA_ENUM_TYPE size) 
	{
		mem->size = size;
	}

	__INLINE void Tag::SetPosition(INMOST_DATA_ENUM_TYPE pos, 
		ElementType type) 
	{
		mem->pos[ElementNum(type)] = pos;
	}

	__INLINE INMOST_DATA_ENUM_TYPE Tag::GetPosition(ElementType type) const 
	{
		assert(mem != NULL); 
		return mem->pos[ElementNum(type)];
	}

	__INLINE void Tag::SetSparse(ElementType type) 
	{
		mem->sparse[ElementNum(type)] = true;
	}

	__INLINE INMOST_DATA_ENUM_TYPE Tag::GetPositionByDim(INMOST_DATA_ENUM_TYPE typenum) const 
	{
		return mem->pos[typenum];
	}

	__INLINE bool Tag::operator <(const Tag & other) const 
	{
		return mem < other.mem;
	}

	__INLINE bool Tag::operator >(const Tag & other) const 
	{
		return mem > other.mem;
	}

	__INLINE bool Tag::operator ==(const Tag & other) const 
	{
		return mem == other.mem;
	}

	__INLINE bool Tag::operator !=(const Tag & other) const 
	{
		return mem != other.mem;
	}

	__INLINE Tag & Tag::operator =(Tag const & other) 
	{
		mem = other.mem; 
		return *this;	
	}

	__INLINE DataType Tag::GetDataType() const 
	{
		assert(mem!=NULL); 
		return mem->dtype;
	}

	__INLINE INMOST_MPI_Type Tag::GetBulkDataType() const 
	{
		assert(mem!=NULL); 
		return mem->bulk_data_type;
	}

	__INLINE INMOST_DATA_ENUM_TYPE Tag::GetBytesSize() const 
	{
		assert(mem!=NULL); 
		return mem->bytes_size;
	}
	__INLINE INMOST_DATA_ENUM_TYPE Tag::GetSize() const 
	{
		assert(mem!=NULL); 
		return mem->size;
	}
	__INLINE std::string Tag::GetTagName() const 
	{
		assert(mem!=NULL); 
		return mem->tagname;
	}
	__INLINE bool Tag::isDefined(ElementType type) const 
	{
		assert(mem!=NULL);
		assert(OneType(type)); 
		return GetPosition(type) != ENUMUNDEF;
	}
	__INLINE bool Tag::isSparse(ElementType type) const 
	{
		assert(mem!=NULL);
		assert(OneType(type)); 
		return mem->sparse[ElementNum(type)];
	}
	__INLINE bool Tag::isValid() const 
	{
		return mem != NULL;
	}
	__INLINE Mesh * Tag::GetMeshLink() const 
	{
		assert(mem!=NULL); 
		return mem->m_link;
	}		
	__INLINE bool Tag::isSparseByDim(INMOST_DATA_INTEGER_TYPE typenum) const 
	{
		assert(mem!=NULL); 
		return mem->sparse[typenum];
	}
	__INLINE bool Tag::isDefinedByDim(INMOST_DATA_INTEGER_TYPE typenum) const 
	{
		assert(mem!=NULL); 
		return GetPositionByDim(typenum) != ENUMUNDEF;
	}
	__INLINE void Tag::SetBulkDataType(INMOST_MPI_Type type) 
	{
		assert(mem!=NULL);
		assert(mem->dtype == DATA_BULK ); 
		mem->bulk_data_type = type;
	}

}

//Implementation of inlined functions
//#include "Source/Data/tag_inline.hpp"


#endif



#endif //INMOST_DATA_H_INCLUDED
