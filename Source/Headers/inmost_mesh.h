
#ifndef INMOST_MESH_H_INCLUDED
#define INMOST_MESH_H_INCLUDED

#include "inmost_common.h"
#include "inmost_data.h"

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


	
	typedef INMOST_DATA_BULK_TYPE GeometricData;
	static const GeometricData CENTROID    = 0;
	static const GeometricData NORMAL      = 1;
	static const GeometricData ORIENTATION = 2;
	static const GeometricData MEASURE     = 3;
	static const GeometricData BARYCENTER  = 4;
	
	typedef INMOST_DATA_BULK_TYPE SyncBitOp; //< This type is used for marker synchronization
	static const SyncBitOp SYNC_BIT_SET = 0;
	static const SyncBitOp SYNC_BIT_OR  = 1;
	static const SyncBitOp SYNC_BIT_XOR = 2;
	static const SyncBitOp SYNC_BIT_AND = 3;

    

	///Definition for data type of topology error or event. This type may be extended later to 64 bits
	/// to accomodate more topological errors or events.
	typedef INMOST_DATA_ENUM_TYPE         TopologyCheck;
	///Throw TopologyError exception on error.
	static const TopologyCheck            THROW_EXCEPTION        = 0x00000001;
	///Print topology notify to std::cerr.
	static const TopologyCheck            PRINT_NOTIFY           = 0x00000002;
	///Element should be deleted if there is an error in EndTopologyCheck.
	static const TopologyCheck            DELETE_ON_ERROR        = 0x00000004;
	///Bitwise OR of errors is recorded to sparse integer tag of element (if not deleted), availible
	/// by TopologyErrorTag().
	static const TopologyCheck            MARK_ON_ERROR          = 0x00000008;
	///Check that edge already exists, on occurance CreateEdge silently returns existant edge.
	static const TopologyCheck            DUPLICATE_EDGE         = 0x00000010;
	///Check that face already exists, on occurance CreateFace silently returns existant face.
	static const TopologyCheck            DUPLICATE_FACE         = 0x00000020;
	///Check that cell already exists, on occurance CreateCell silently returns existant cell.
	static const TopologyCheck            DUPLICATE_CELL         = 0x00000040;
	///Check that edge have more then two nodes (no support for this kind of object).
	static const TopologyCheck            DEGENERATE_EDGE        = 0x00000080;
	///Produce error if face consists of less then 3 edges.
	static const TopologyCheck            DEGENERATE_FACE        = 0x00000100;
	///Produce error if cell consists of less then 4 faces.
	static const TopologyCheck            DEGENERATE_CELL        = 0x00000200;
	///Produce error if face have wrong orientation of edges, for star-shaped elements only.
	static const TopologyCheck            FACE_ORIENTATION       = 0x00000400;
	///Produce error if face is non planar.
	static const TopologyCheck            FACE_PLANARITY         = 0x00000800;
	///Produce error if there is another face that already uses same nodes.
	static const TopologyCheck            INTERLEAVED_FACES      = 0x00001000;
	///Check that every face have exectly two neighbours, if DUPLICATE_CELL is activated, then checks
	/// for cell duplication.
	static const TopologyCheck            TRIPLE_SHARED_FACE     = 0x00002000;
	///Produce error if one of the faces of the cell contains all the nodes of the cell.
	static const TopologyCheck            FLATTENED_CELL         = 0x00004000;
	///Produce error if provided array of elements for construction contain duplications.
	static const TopologyCheck            ADJACENT_DUPLICATE     = 0x00008000;
	///Hidden elements should not be used when new elements are created.
	static const TopologyCheck            ADJACENT_HIDDEN        = 0x00010000;
	///Check that all handles are valid.
	static const TopologyCheck            ADJACENT_VALID         = 0x00020000;
	///Produce error if provided array of elements have wrong geometric dimension.
	static const TopologyCheck            ADJACENT_DIMENSION     = 0x00040000;
	///Produce error if edges of faces are not closed, or face is folded.
	/// \warning This check needs NEED_TEST_CLOSURE to be activated.
	static const TopologyCheck            PROHIBIT_MULTILINE     = 0x00080000;
	/// Allow only known types of elements: Tet, Quad, Line.
	static const TopologyCheck            PROHIBIT_POLYGON       = 0x00100000;
	///Produce error if faces of cell are not closed, or cell is folded.
	/// \warning This check needs NEED_TEST_CLOSURE to be activated.
	static const TopologyCheck            PROHIBIT_MULTIPOLYGON  = 0x00200000;
	///Allow only known types of elements: Tet,Hex,Prism,Pyramid.
	static const TopologyCheck            PROHIBIT_POLYHEDRON    = 0x00400000;
	///Edges of the face should form one closed loop.
	static const TopologyCheck            FACE_EDGES_ORDER       = 0x00800000;
	///Don't allow concave face.
	/// \warning Not implemented.
	static const TopologyCheck            PROHIBIT_CONCAVE_FACE  = 0x01000000;
	///Don't allow concave cell.
	/// \warning Not implemented.
	static const TopologyCheck            PROHIBIT_CONCAVE_CELL  = 0x02000000;
	///Don't allow non-star shaped face.
	/// \warning Not implemented.
	static const TopologyCheck            PROHIBIT_NONSTAR_FACE  = 0x04000000;
	///Don't allow non-star shaped concave cell.
	/// \warning Not implemented.
	static const TopologyCheck            PROHIBIT_NONSTAR_CELL  = 0x08000000;
	///Edges of the face don't cross each other.
	/// \warning Not implemented.
	static const TopologyCheck            FACE_SELF_INTERSECTION = 0x10000000;
	///Faces of the cell don't cross each other.
	/// \warning Not implemented.
	static const TopologyCheck            CELL_SELF_INTERSECTION = 0x20000000;
	///Silent, test's for closure in ComputeGeometricType, needed to detect MultiLine and
	/// MultiPolygon.
	static const TopologyCheck            NEED_TEST_CLOSURE      = 0x40000000;
	///Don't allow 2d grids, where edges appear to be vertexes, faces are edges and cells are faces
	static const TopologyCheck            DISABLE_2D             = 0x80000000;
	///A shortcut to test for grid conformity.
	static const TopologyCheck            GRID_CONFORMITY        = NEED_TEST_CLOSURE 
                                                               | PROHIBIT_MULTILINE 
                                                               | PROHIBIT_MULTIPOLYGON  
    //                                                           | INTERLEAVED_FACES 
                                                               | TRIPLE_SHARED_FACE;
	///Default set of options.
	static const TopologyCheck            DEFAULT_CHECK          = THROW_EXCEPTION 
                                                               | DUPLICATE_EDGE 
                                                               | DUPLICATE_FACE 
                                                               | DUPLICATE_CELL 
                                                               | PRINT_NOTIFY;
	///Returns a string explaining each topology check. Implemented in mesh.cpp.
	/// @param c Single topology check.
	const char * TopologyCheckNotifyString(TopologyCheck c);
  
 
	
	template <typename StorageType>
	class ElementArray
	{
	public:
		typedef std::vector<HandleType>      cont_t;
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
			ptrdiff_t     operator -(const iterator & other) const {return static_cast<const cont_t::iterator>(*this)-static_cast<const cont_t::iterator>(other);}
			iterator      operator +(size_t n) const{return iterator(m_link,cont_t::iterator::operator +(n));}
			iterator      operator -(size_t n) const{return iterator(m_link,cont_t::iterator::operator -(n));}
			iterator &    operator ++() {cont_t::iterator::operator++(); return *this;}
			iterator      operator ++(int) {iterator ret(*this); cont_t::iterator::operator++(); return ret;}
			iterator &    operator --() {cont_t::iterator::operator--(); return *this;}
			iterator      operator --(int) {iterator ret(*this); cont_t::iterator::operator--(); return ret;}
			iterator &    operator =(iterator const & other) {m_link = other.m_link; cont_t::iterator::operator=(static_cast<cont_t::iterator const &>(other)); return *this; }
			HandleType &  operator *() const { return cont_t::iterator::operator *(); }
			StorageType   operator->() const { return StorageType(m_link,&cont_t::iterator::operator *()); }
		};
		class reverse_iterator : public cont_t::reverse_iterator
		{
			Mesh * m_link;
		public:
			reverse_iterator(Mesh * m, const cont_t::reverse_iterator & other ) : cont_t::reverse_iterator(other), m_link(m) {assert(m_link);}
			reverse_iterator(const reverse_iterator & other ) : cont_t::reverse_iterator(other), m_link(other.m_link) {assert(m_link);}
			ptrdiff_t             operator -(const reverse_iterator & other) const {return static_cast<const cont_t::reverse_iterator>(*this)-static_cast<const cont_t::reverse_iterator>(other);}
			reverse_iterator      operator +(size_t n) const {return reverse_iterator(m_link,cont_t::reverse_iterator::operator+(n));}
			reverse_iterator      operator -(size_t n) const {return reverse_iterator(m_link,cont_t::reverse_iterator::operator-(n));}
			reverse_iterator &    operator ++() {cont_t::reverse_iterator::operator++(); return *this;}
			reverse_iterator      operator ++(int) {reverse_iterator ret(*this); cont_t::reverse_iterator::operator++(); return ret;}
			reverse_iterator &    operator --() {cont_t::reverse_iterator::operator--(); return *this;}
			reverse_iterator      operator --(int) {reverse_iterator ret(*this); cont_t::reverse_iterator::operator--(); return ret;}
			reverse_iterator & operator =(reverse_iterator const & other) {m_link = other.m_link; cont_t::reverse_iterator::operator=(static_cast<cont_t::reverse_iterator const &>(other)); return *this; }
			HandleType &       operator *() const { return cont_t::reverse_iterator::operator *(); }
			StorageType        operator->() const { return StorageType(m_link,&cont_t::reverse_iterator::operator *()); }
		};
		class const_iterator : public cont_t::const_iterator
		{
			Mesh * m_link;
		public:
			const_iterator(Mesh * m, const cont_t::const_iterator & other ) : cont_t::const_iterator(other), m_link(m) {assert(m_link);}
			const_iterator(const const_iterator & other ) : cont_t::const_iterator(other), m_link(other.m_link) {assert(m_link);}
			ptrdiff_t     operator -(const const_iterator & other) const {return static_cast<const cont_t::const_iterator>(*this)-static_cast<const cont_t::const_iterator>(other);}
			const_iterator &    operator ++() {cont_t::const_iterator::operator++(); return *this;}
			const_iterator      operator ++(int) {const_iterator ret(*this); cont_t::const_iterator::operator++(); return ret;}
			const_iterator &    operator --() {cont_t::const_iterator::operator--(); return *this;}
			const_iterator      operator --(int) {const_iterator ret(*this); cont_t::const_iterator::operator--(); return ret;}
			const_iterator &    operator =(const_iterator const & other) {m_link = other.m_link; cont_t::const_iterator::operator=(static_cast<cont_t::const_iterator const &>(other)); return *this; }
			const HandleType &  operator *() const { return cont_t::const_iterator::operator *(); }
			StorageType         operator->() const { return StorageType(m_link,cont_t::const_iterator::operator *()); }
		};
		class const_reverse_iterator : public cont_t::const_reverse_iterator
		{
			Mesh * m_link;
		public:
			const_reverse_iterator(Mesh * m, const cont_t::const_reverse_iterator & other) : cont_t::const_reverse_iterator(other), m_link(m) {assert(m_link);}
			const_reverse_iterator(const const_reverse_iterator & other ) : cont_t::const_reverse_iterator(other), m_link(other.m_link) {assert(m_link);}
			ptrdiff_t     operator -(const const_reverse_iterator & other) const {return static_cast<const cont_t::const_reverse_iterator>(*this)-static_cast<const cont_t::const_reverse_iterator>(other);}
			const_reverse_iterator &    operator ++() {cont_t::const_reverse_iterator::operator++(); return *this;}
			const_reverse_iterator      operator ++(int) {const_reverse_iterator ret(*this); cont_t::const_reverse_iterator::operator++(); return ret;}
			const_reverse_iterator &    operator --() {cont_t::const_reverse_iterator::operator--(); return *this;}
			const_reverse_iterator      operator --(int) {const_reverse_iterator ret(*this); cont_t::const_reverse_iterator::operator--(); return ret;}
			const_reverse_iterator & operator =(const_reverse_iterator const & other) { cont_t::const_reverse_iterator::operator=(static_cast<cont_t::const_reverse_iterator const &>(other)); return *this; }
			const HandleType &       operator *() const { return cont_t::const_reverse_iterator::operator *(); }
			StorageType              operator->() const { return StorageType(m_link,cont_t::const_reverse_iterator::operator *()); }
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
		void                      SetPrivateMarker   (MarkerType n) const;
		void                      RemPrivateMarker   (MarkerType n) const;
		template<typename Etype>
		ElementArray<Etype>       Convert() {return ElementArray<Etype>(m_link,container);}
	};
	
	// Abstract function for calculation of mean value on element
	struct MeanFunc { virtual Storage::real func(Storage::real* x, Storage::real t) const = 0; };
	struct MeanFuncRaw : public MeanFunc
	{ 
		Storage::real(*pfunc)(Storage::real* x, Storage::real t);
		MeanFuncRaw(Storage::real(*pfunc)(Storage::real* x, Storage::real t)) : pfunc(pfunc) {}
		Storage::real func(Storage::real* x, Storage::real t) const { return pfunc(x, t); }
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
		/// Retrieve number of adjacent elements.
		/// For etype you can either pass one type as CELL,
		/// or several types as bitwise mask: NODE | CELL.
		/// @param etype bitwise mask of element types
		/// @return number of adjacent elements.
		virtual enumerator                nbAdjElements           (ElementType etype) const;
		/// Retrieve unordered array of adjacent elements.
		/// If you care about orderness of nodes for face you should you Face::getNodes() instead.
		/// If you want faster access you may use direct access to handles stored in memory
		/// through Mesh::HighConn for upper adjacencies and Mesh::LowConn for lower adjacencies.
		/// @param etype bitwise mask of element types
		/// @return array of elements
		virtual ElementArray<Element>     getAdjElements          (ElementType etype) const;  //unordered
		/// Retrieve number of adjacent elements with marker.
		/// As etype you can either pass one type as CELL,
		/// or several types as bitwise mask: NODE | CELL
		/// @param etype bitwise mask of element types
		/// @param mask marker to be set 
		/// @param invert_mask if true then those are selected on wich marker is not set
		/// @return number of adjacent elements.
		virtual enumerator                nbAdjElements           (ElementType etype, MarkerType mask, bool invert_mask = false) const;
		/// Retrieve unordered array of adjacent elements with marker.
		/// @param etype bitwise mask of element types
		/// @param mask marker to be set 
		/// @param invert_mask if true then those are selected on wich marker is not set
		/// @return array of elements
		virtual ElementArray<Element>     getAdjElements          (ElementType etype, MarkerType mask, bool invert_mask = false) const;  //unordered
		ElementArray<Element>       BridgeAdjacencies       (ElementType Bridge, ElementType Dest, MarkerType bridge_mask = 0, bool bridge_invert = false, MarkerType target_mask = 0, bool target_invert = false) const;
		ElementArray<Node>          BridgeAdjacencies2Node  (ElementType Bridge, MarkerType bridge_mask = 0, bool bridge_invert = false, MarkerType target_mask = 0, bool target_invert = false) const;
		ElementArray<Edge>          BridgeAdjacencies2Edge  (ElementType Bridge, MarkerType bridge_mask = 0, bool bridge_invert = false, MarkerType target_mask = 0, bool target_invert = false) const;
		ElementArray<Face>          BridgeAdjacencies2Face  (ElementType Bridge, MarkerType bridge_mask = 0, bool bridge_invert = false, MarkerType target_mask = 0, bool target_invert = false) const;
		ElementArray<Cell>          BridgeAdjacencies2Cell  (ElementType Bridge, MarkerType bridge_mask = 0, bool bridge_invert = false, MarkerType target_mask = 0, bool target_invert = false) const;
		/// Retrieve all the nodes of the element.
		///
		/// For a node returns itself.
		///
		/// For an edge returns ordered pair of nodes. The order of nodes in the edge is preserved from the first construction.
		///
		/// For a face returns ordered set of nodes. In the case Face::CheckNormalOrientation returns true the
		/// order of nodes will follow right hand side rule with respect to normal vector, otherwise it follows
		/// left hand side rule with respect to normal vector.
		///
		/// For a cell returns the same order that was provided through suggest_nodes_oreder in Mesh::CreateCell.
		/// In the case suggest_nodes_order was not provided, the order of nodes follows VTK format for known types
		/// of elements such as Element::Tet, Element::Prism, Element::Hex, Element::Pyramid. For a general polyhedron
		/// the order is unspecified.
		///
		/// @return array of elements
		/// @see Face::CheckNormalOrientation
		/// @see Face::FaceOrientedOutside
		virtual ElementArray<Node>  getNodes                () const; //unordered
		/// Retrieve all the edges of the element.
		///
		/// For a node returns unordered set of edges.
		///
		/// For an edge returns itself.
		///
		/// For a face returns ordered set of edges.
		///
		/// For a cell returns unordered set of edges.
		///
		/// @return array of elements
		virtual ElementArray<Edge>  getEdges                () const; //unordered
		/// Retrieve all the faces of the element.
		///
		/// For a node returns unordered set of faces.
		///
		/// For an edge returns unordered set of faces.
		///
		/// For a face returns itself.
		///
		/// For a cell return ordered set of faces. The order of faces in the cell is preserved from the first construction.
		///
		/// @return array of elements
		virtual ElementArray<Face>  getFaces                () const; //unordered
		/// Return all the cells of the element.
		///
		/// For a node returns unordered set of cells.
		///
		/// For an edge returns unordered set of cells.
		///
		/// For a face returns a pair of cells. In the case Face::CheckNormalOrientation returns true
		/// then the normal points from the first cell to the second and in oppisite direction otherwise.
		///
		/// For a cell returns itself.
		///
		/// @return array of elements
		/// @see Face::FaceOrientedOutside
		virtual ElementArray<Cell>  getCells                () const; //unordered
		virtual ElementArray<Node>  getNodes                (MarkerType mask,bool invert_mask = false) const; //unordered
		virtual ElementArray<Edge>  getEdges                (MarkerType mask,bool invert_mask = false) const; //unordered
		virtual ElementArray<Face>  getFaces                (MarkerType mask,bool invert_mask = false) const; //unordered
		virtual ElementArray<Cell>  getCells                (MarkerType mask,bool invert_mask = false) const; //unordered
		GeometricType               GetGeometricType        () const;
		unsigned int                GetElementDimension     () const {return GetGeometricDimension(GetGeometricType());}
		Status                      GetStatus               () const;
		void                        SetStatus               (Status status) const;
		Storage::integer &          GlobalID                () const;
		bool                        CheckElementConnectivity() const;
		void                        PrintElementConnectivity() const;
		static bool                 CheckConnectivity       (Mesh * m);
		//implemented in geometry.cpp
		void                        CastRay                 (const real * pos, const real * dir, std::map<HandleType, real> & hits) const;
		
		void                        ComputeGeometricType    () const;
		void                        Centroid                (real * cnt) const;
		void                        Barycenter              (real * cnt) const;
		Storage::real               Mean                    (const MeanFunc & f, real time) const;
		Storage::real               Mean                    (real(*func)(real* x, real t), real time) const { return Mean(MeanFuncRaw(func), time); }
		/// Determine that the element is on the boundary.
		/// \warning This function does not involve communication
		/// for distributed mesh and works only for local partition
		/// of the mesh. It cannot correctly determine boundary elements
		/// of the global mesh. For this purpose you should use
		/// Mesh::MarkBoundaryFaces
		/// @return True if the element is on the boundary, otherwise false.
		/// @see Mesh::MarkBoundaryFaces
		bool                        Boundary                () const;
		bool                        Planarity               () const; // check that all nodes lay on one plane
		//implemented in modify.cpp
		/// If the function returns true then element was hidden,
		/// works only inside BeginModification and EndModification,
		/// on EndModification all Hidden elements are deleted.
		/// \todo Probably have to check normal orientation after
		/// hiding a back cell for a face.
		bool                        Hide                    () const;
		/// If the function returns true then element was recovered
		/// from hidden state, works only inside BeginModification
		/// and EndModification.
		/// \todo Probably have to check normal orientation after
		/// recovering a back cell for a face.
		bool                        Show                    () const; // if true then element was recovered
		/// Remove element from mesh.
		/// If you call this function inside modification phase, see Mesh::BeginModification and Mesh::EndModification,
		/// and the element was not created during modification phase (not marked as Element::New),
		/// then the element will not be actually destroyed but hidden. 
		/// You can restore all the hidden elements by using Mesh::ToggleModification.
		/// \warning This function will not resolve an ierarchical strucutre of ElementSet, use ElemnetSet::DeleteSet instead.
		/// @return Returns true if the element was actually destroyed. Returns false if the element was hidden.
		bool                        Delete                  (); // if true then element was deleted, otherwise it was hidden
		bool                        Hidden                  () const;
		bool                        New                     () const;
		void                        Disconnect              (bool delete_upper_adjacent) const; //disconnect all elements, delete upper dependent
		/// Disconnect element.
		/// Disconnects nodes from this edge, edges from this face, faces from this cell, cannot disconnect cells from this node; \n
		/// Disconnects edges from this node, faces from this edge, cells from this face, cannot disconnect nodes from this cell; \n
		/// Updates geometric data and cell nodes automatically.
		void                        Disconnect              (const HandleType * adjacent, INMOST_DATA_ENUM_TYPE num) const;
		/// \brief Connects lower adjacencies to current element.
		/// Geometric data and cell nodes are updated automatically.
		///
		/// \todo
		///	  1. Asserts in this function should be replaced by Topography checks.
		///	  2. This function should be used for construction of elements instead of current implementation.
		///	  3. Should correctly account for order of edges (may be implemented through CheckEdgeOrder, FixEdgeOrder).
		void                        Connect                 (const HandleType * adjacent, INMOST_DATA_ENUM_TYPE num) const; 
		/// Update geometric data for element, calls RecomputeGeometricData from Mesh.
		//void                        UpdateGeometricData     () const;
		/// Marks element to be sent to remote processors that current processor don't belong to.
		/// Call Mesh::ExchangeMarked to perform the exchange.
		void                        SendTo                  (std::set<Storage::integer> & procs) const; 
		void                        SendTo                  (std::vector<Storage::integer> & procs) const; 
		void                        SendTo                  (Storage::integer_array procs) const; 
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
		static Edge                 UniteEdges              (const ElementArray<Edge> & edges, MarkerType del_protect);
		static bool                 TestUniteEdges          (const ElementArray<Edge> & edges, MarkerType del_protect);
		static ElementArray<Edge>   SplitEdge               (Edge e, const ElementArray<Node> & nodes, MarkerType del_protect); //provide ordered array of nodes, that lay between former nodes of the edge
		static bool                 TestSplitEdge           (Edge e, const ElementArray<Node> & nodes, MarkerType del_protect);
		static Node                 CollapseEdge            (Edge e);
		//implemented in geometry.cpp
		Storage::real               Length                  () const;
		///Swap positions of first node and last node
		void                        SwapEnds                ();
	};

	__INLINE const Edge & InvalidEdge() {static Edge ret(NULL,InvalidHandle()); return ret;}
	

	/// \brief An interface for elements of type FACE.
	///
	/// This interface carry the link to the mesh and the element's handle that
	/// represents position of the element's data in the mesh. 
	///
	/// Interface provides some operations that can be done uniquely on the object of class 
	/// Element for which Element::GetElementType retruns FACE.
	///
	/// For the basic set of operations on all of the elements check class Element.
	///
	/// For the basic set of operations on the data of the element check class Storage.
	///
	/// You can obtain object of class Face from Mesh::iteratorFace,
	/// in this case obtained object is always valid.
	/// Also you can get it through Mesh::FaceByLocalID, check with Element::isValid to see whether you
	/// have obtained a valid object.
	/// You can convert an object of class Element into an object of class Face by Element::getAsFace. In debug
	/// mode it will internally check that the element is of type FACE.
	/// You may compose an object of class Face by using constructor and specifing pointer to the mesh
	/// and providing element's handle. You can make a handle with ComposeHandle(FACE,local_id) function.
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
		/// \brief Retrieve the cell for which the normal points outwards.
		///
		/// \warning Depending on the grid construction the normal may incorrectly point inwards.
		/// You can resolve this situation by Face::FixNormalOrientation.
		///
		/// @return the cell for which normal points outwards.
		/// @see Face::FixNormalOrientation
		Cell                        BackCell                () const;
		/// \brief Retrieve the cell for which the normal points inwards.
		///
		/// \warning Depending on the grid construction the normal may incorrectly point outwards.
		/// You can resolve this situation by Face::FixNormalOrientation.
		///
		/// @return the cell for which normal points inwards.
		/// @see Face::FixNormalOrientation
		Cell                        FrontCell               () const;
		bool                        FaceOrientedOutside     (Cell c) const;
		void                        ReorderEdges            () const;
		bool                        CheckEdgeOrder          () const; //not implemented// returns true if edges of face form an ordered closed loop
		bool                        FixEdgeOrder            () const; //not implemented// returns true if edges were successfully reordered to form a closed loop
		//implemented in modify.cpp
		static Face                 UniteFaces              (const ElementArray<Face> & faces, MarkerType del_protect);
		static bool                 TestUniteFaces          (const ElementArray<Face> & faces,  MarkerType del_protect);
		static ElementArray<Face>   SplitFace               (Face face, const ElementArray<Edge> & edges, MarkerType del_protect); //provide all edges that lay inside face
		static bool                 TestSplitFace           (Face face, const ElementArray<Edge> & edges, MarkerType del_protect);
		///This function changes Face::BackCell() with Face::FrontCell().
		/// \warning Function does not modify normal orientation. You may have
		/// to call Face::FixNormalOrientation().
		void                        SwapCells               () const; 
		//implemented in geometry.cpp
		Storage::real               Area                    () const;
		void                        Normal                  (real * nrm) const;
		void                        UnitNormal              (real * nrm) const;
		void                        OrientedNormal          (Cell c, real * nrm) const;
		void                        OrientedUnitNormal      (Cell c, real * nrm) const;
		bool                        FixNormalOrientation    (bool allow_swap = true) const;  //returns true if orientation was corrected, otherwise returns false
		bool                        CheckNormalOrientation  () const; //returns true if orientation is correct, otherwise returns false
		bool                        Closure                 () const; // test integrity of polygon
		bool                        Inside                  (const real * point) const; //is point inside face
	};

	__INLINE const Face & InvalidFace() {static Face ret(NULL,InvalidHandle()); return ret;}
	
	//implemented in cell.cpp
	/// \brief An interface for elements of type CELL.
	///
	/// This interface carry the link to the mesh and the element's handle that
	/// represents position of the element's data in the mesh. 
	///
	/// Interface provides some operations that can be done uniquely on the object of class 
	/// Element for which Element::GetElementType retruns CELL.
	///
	/// For the basic set of operations on all of the elements check class Element.
	///
	/// For the basic set of operations on the data of the element check class Storage.
	///
	/// You can obtain object of class Cell from Mesh::iteratorCell,
	/// in this case obtained object is always valid.
	/// Also you can get it through Mesh::CellByLocalID, check with Element::isValid to see whether you
	/// have obtained valid object.
	/// You can convert an object of class Element into an object of class Cell by Element::getAsCell. In debug
	/// mode it will internally check that the element is of type CELL.
	/// You may compose an object of class Cell by using constructor and specifing pointer to the mesh
	/// and providing element's handle. You can make a handle with ComposeHandle(CELL,local_id) function.
	class Cell : public Element
	{
	public:
		/// \brief Basic constructor.
		///
		/// Constructs invalid cell.
		Cell() : Element() {} 
		/// \brief Basic constructor with fixed handle.
		///
		/// When constructed this way, only handle within object may change.
		///
		/// @param m Pointer to the mesh to which the element belongs.
		/// @param h Unique handle that describes position of the element within the mesh.
		Cell(Mesh * m, HandleType h) : Element(m,h) {}
		/// \brief Basic constructor with an assignable handle.
		///
		/// When constructed this way, the provided memory location for handle
		/// will be modified on assignment.
		///
		/// The purpose of this function is to be used in various non-constant iterators
		/// of containers that allow underlaying contents to be changed.
		///
		/// @param m Pointer to the mesh to which the element belongs.
		/// @param h Pointer to the handle that describes position of the element within mesh.
		Cell(Mesh * m, HandleType * h) : Element(m,h) {}
		/// \brief Copy constructor.
		///
		/// New object will inherit the assignment behavior of the initial object.
		///
		/// @param other Object of type Cell to be duplicated.
		Cell(const Cell & other) : Element(other) {}
		/// \brief Assignment operator
		///
		/// Assigned object will inherit the asignment behavior of the initial object only
		/// in case it was constructed without pointer to a handle. Otherwise it will
		/// only modify the memory location to which it points.
		///
		/// @param other Object of type Cell to be duplicated.
		/// @return Reference to the current object.
		Cell & operator =(Cell const & other) {Element::operator =(other); return *this;} ///< Assignment Operator.
		/// \brief Operator of dereference to pointer.
		///
		/// This is needed for iterators to work properly.
		///
		/// @return Pointer to the current object.
		Cell * operator ->() {return this;} 
		/// \brief Operator of dereference to constant pointer. 
		///
		/// This is needed for const_iterators to work properly.
		///
		/// @return Constant pointer to the current object.
		const Cell * operator ->() const {return this;}
		/// \brief Get self-reference.
		///
		/// Main purpose is to convert iterators into elements and to pass them as arguments of functions.
		///
		/// @return Constant reference to the current object.
		Cell & self() {return *this;} 
		/// \brief Get constant self-reference.
		///
		/// Main purpose is to convert iterators into elements and to pass them as arguments of functions.
		///
		/// @return Reference to the current object.
		const Cell & self() const {return *this;}
	public:
		/// \brief Get all the nodes of the current cell.
		///
		/// This function traverses up the adjacency graph by one level.
		///
		/// The order of nodes will be preserved from the moment of the construction of the cell.
		/// When suggest_nodes array is not supplied into the Mesh::CreateCell functions, then
		/// for known types the order of nodes follows VTK convention. For polyhedral cells
		/// the order is unspecified. When suggest_nodes_order was provided into Mesh::CreateCell
		/// then the order follows provided order.
		///
		/// @return Set of nodes that compose current cell.
		ElementArray<Node>          getNodes                () const;
		/// \brief Get all the edges of the current cell.
		///
		/// This function traverses down the adjacency graph by two levels.
		///
		/// The order of the edges is unspecified.
		///
		/// @return Set of edges that compose current cell.
		///
		/// \todo
		/// One should thoroughly check three scenarios of function execution in shared parallel environment
		/// for different types of cells (simple tet/hex cells as well as complex polyhedral cells) and
		/// draw a conclusion on the best scenario for each condition. One of the development versions
		/// contains all the algorithms, ask for the files.
		///  1. Use of markers (current wariant).
		///  2. Put all elements into array with duplications, then run std::sort and std::unique.
		///  3. Put all elements into array, check for duplication by running through array.
		///
		/// \warning Note that this function uses markers to check for duplication of edges
		/// in output array. This allows for faster algorithm, but will lead to penalties
		/// in shared parallel execution, since operations on markers are declared atomic.
		/// For atomic declaration to work you have to define USE_OMP during CMake configuration
		/// or in inmost_common.h. If you attempt to use this function in shared parallel execution
		/// without USE_OMP you may expect side effects.
		ElementArray<Edge>          getEdges                () const;
		/// \brief Get all the faces of the current cell.
		///
		/// This function traverses down the adjacency graph by one level.
		///
		/// The order of the faces returned by this function is preserved from 
		/// the moment of the construction of the cell.
		///
		/// @return Set of faces that compose current cell.
		ElementArray<Face>          getFaces                () const;
		/// \brief Get the subset of the nodes of the current cell that are (not) marked by provided marker.
		///
		/// This function traverses up the adjacency graph by one level.
		///
		/// The order of nodes will be preserved from the moment of the construction of the cell.
		/// When suggest_nodes array is not supplied into the Mesh::CreateCell functions, then
		/// for known types the order of nodes follows VTK convention. For polyhedral cells
		/// the order is unspecified. When suggest_nodes_order was provided into Mesh::CreateCell
		/// then the order follows provided order.
		///
		/// @param mask Marker that should be used to filter elements.
		/// @param invert_mask If false then those elements that are marked will be taken, otherwise
		///                    elements that are not marked will be taken.
		/// @return Set of the nodes that compose current cell and (not) marked by marker.
		///
		/// \warning To work correctly in shared parallel environment this function requires
		/// USE_OMP to be defined in CMake or in inmost_common.h.
		ElementArray<Node>          getNodes                (MarkerType mask,bool invert_mask = false) const;
		/// \brief Get the subset of the edges of the current cell that are (not) marked by provided marker.
		///
		/// This function traverses down the adjacency graph by two levels.
		///
		/// The order of the edges is unspecified.
		///
		/// @param mask Marker that should be used to filter elements.
		/// @param invert_mask If false then those elements that are marked will be taken, otherwise
		///                    elements that are not marked will be taken.
		/// @return Set of the edges that compose current cell and (not) marked by marker.
		///
		/// \warning To work correctly in shared parallel environment this function requires
		/// USE_OMP to be defined in CMake or in inmost_common.h.
		ElementArray<Edge>          getEdges                (MarkerType mask,bool invert_mask = false) const;
		/// \brief Get the subset of the faces of the current cell that are (not) marked by provided marker.
		///
		/// This function traverses down the adjacency graph by one level.
		///
		/// The order of the faces returned by this function is preserved from 
		/// the moment of the construction of the cell.
		///
		/// @param mask Marker that should be used to filter elements.
		/// @param invert_mask If false then those elements that are marked will be taken, otherwise
		///                    elements that are not marked will be taken.
		/// @return Set of the faces that compose current cell and (not) marked by marker.
		///
		/// \warning To work correctly in shared parallel environment this function requires
		/// USE_OMP to be defined in CMake or in inmost_common.h.
		ElementArray<Face>          getFaces                (MarkerType mask,bool invert_mask = false) const;
		/// \brief Check that sequence of edges form a closed loop and each edge have a node that matches one of the nodes at the next edge.
		///
		/// This function works for cells for which Element::GetElementDimension returns 2.
		///
		/// @return True if edges form the correct closed loop.
		///
		/// \todo
		///   1. Implement.
		///   2. Use in topology check algorithms.
		///   3. Use in Element::Connect.
		bool                        CheckEdgeOrder          () const; //not implemented//2D only, returns true if edges of face form an ordered closed loop
		/// \brief Repair the sequence of edges so that each edge have node that matches one of the nodes at the next edge.
		///
		/// This function works for cells for which Element::GetGeometricDimension returns 2.
		///
		/// This function may fail if all the edges does not form a closed loop.
		///
		/// @return True if edges were successfully reordered to form a closed loop.
		///
		/// \todo
		///   1. Implement.
		///   2. Use in topology check algorithms.
		///   3. Use in Element::Connect.
		bool                        FixEdgeOrder            () const; //not implemented//2D only, returns true if edges were successfully reordered to form a closed loop
		//implemented in modify.cpp
		/// \brief Unite a set of given cells into one cell.
		///
		/// This will create a cell whose faces are formed by symmetric difference
		/// of faces of given cells.
		/// If you specify a nonzero marker then the procedure will fail 
		/// if any marked element have to be deleted during union.
		/// @param cells A set of cells to be united.
		/// @param del_protect A marker that protects elements from deletion. Zero means no check.
		/// @return A new cell.
		static Cell                 UniteCells              (const ElementArray<Cell> & cells, MarkerType del_protect);
		/// \brief Test that no marked element will be deleted during union of given cells.
		/// @param cells A set of cells to be united.
		/// @param del_protect A marker that protects elements from deletion. Zero means no check.
		/// @return True if no element deleted, false otherwise.
		static bool                 TestUniteCells          (const ElementArray<Cell> & cells, MarkerType del_protect);
		/// \brief Separate a cell according to the set of provided faces.
		///
		/// You should first separate all edges with new nodes by Edge::SplitEdge, then faces with new edges by Face::SplitFace
		/// and then create faces using new edges from faces of initial cell and optionally edges from inside of the cell. 
		/// All faces are supposed to be internal with respect to the current cell. 
		/// Internally this function will resolve geometry of new cells inside of the 
		/// current cell by running through adjacency graph.
		///
		/// @param cell A cell to be split.
		/// @param faces A set of faces, internal for current cell, that will be used for construction of new cells.
		/// @param del_protect Marker that may be used to protect some elements from deletion by algorithm. 
		///        Zero means no check.
		/// @return A set of new cells.
		/// @see Edge::SplitEdge
		/// @see Face::SplitFace
		///
		/// \todo
		///   1. The algorithm inside is minimizing the size of the adjacency graph for each new cell.
		///      The correct behavior is to calculate volume of the cell for each adjacency graph and choose the graph with minimal volume.
		///      This requires calculation of volume for non-convex cells. For correct calculation of volume on non-convex cells one should 
		///      find one face for which normal orientation can be clearly determined and then orient all edges of the cell with respect to
		///      the orientation of edges of this face and establish normals for all faces. Once the algorithm is implemented here it should 
		///      be implemented in geometrical services or vice verse.
		///   2. Probably the algorithm should minimize the volume and adjacency graph size altogether. Between the cells with smallest volume
		///      within some tolerance select those that have smallest adjacency graph.
		static ElementArray<Cell>   SplitCell               (Cell cell, const ElementArray<Face> & faces, MarkerType del_protect); //provide all faces, that lay inside cell
		/// \brief This functions checks is it possible to split the cell by the given set of faces without deleting marked elements.
		///
		/// @param cell A cell to be split.
		/// @param faces Set of faces, internal for current cell, that will be used for new cells.
		/// @param del_protect Marker that may be used to protect some elements from deletion by algorithm. Zero means no check.
		/// @return True if possible, otherwise false.
		static bool                 TestSplitCell           (Cell cell, const ElementArray<Face> & faces, MarkerType del_protect);
		//implemented in geometry.cpp
		/// \brief Get a cell that share a face with the current cell.
		///
		/// Don't forget to check that the returned cell is valid.
		/// It would return invalid cell for boundary face.
		///
		/// @param face A face of current cell for which neighbouring cell is determined.
		/// @return A cell that shares a given face with the current cell.
		Cell                        Neighbour               (Face face) const;
		/// \brief Get all cells that share the face with the current cell.
		/// @return Set of cells that share a face with the current cell.
		ElementArray<Cell>          NeighbouringCells       () const; // get all cells that share any face with current
		/// \brief Determine, if point lies inside element.
		///
		/// Now it works only for 3-dimensional elements, at future it will be
		/// extended to support polygons and segments.
		/// @param point coordinates of the point, it is assumed that number of the
		/// coordinates is the same as the number of dimensions in the mesh.
		/// @return returns true if points inside element, false otherwise
		/// @see Mesh::GetDimensions
		///
		/// \todo
		///   1. Should be checked or extended for 2d cells. (done, testing)
		bool                        Inside                  (const real * point) const; //is point inside cell, check for 2d case
		/// \brief Return volume of the cell.
		///
		/// Note that currently the volume for non-convex cells may be calculated incorrectly.
		/// \todo
		///   1. Geometric services should correctly resolve volume for non-convex cells.
		real                        Volume                  () const;
		/// \brief Test that faces of the cell form the closed set.
		///
		/// This is automatically checked for if you activate NEED_TEST_CLOSURE
		/// in Mesh::SetTopologyCheck.
		/// @return True if faces of the cell form the closed set, false otherwise.
		bool                        Closure                 () const;
		/// For each adjacent cell make me a front cell and fix normal orientation accordingly.
		void                        SwapBackCell            () const;

		bool                        CheckConvexity          () const;
	};

	///
	__INLINE const Cell & InvalidCell() {static Cell ret(NULL,InvalidHandle()); return ret;}

	class ElementSet : public Element //implemented in eset.cpp
	{
	public:
		///Number of reserved positions in HighConn array.
		///The first position is the handle to parent set.
		///The second position is handle to sibling set.
		///The third position is handle to child set.
		///The fourth position is number of sorted elements in the set.
		///All the rest are positions of deleted elements.
		static const enumerator             high_conn_reserved  = 4;
		__INLINE static HandleType &        hParent             (Element::adj_type & arr) {return arr[0];}
		__INLINE static HandleType &        hSibling            (Element::adj_type & arr) {return arr[1];}
		__INLINE static HandleType &        hChild              (Element::adj_type & arr) {return arr[2];}
		__INLINE static HandleType &        hSorted             (Element::adj_type & arr) {return arr[3];}
		typedef INMOST_DATA_BULK_TYPE       ComparatorType;
		static const ComparatorType         UNSORTED_COMPARATOR = 0;
		static const ComparatorType         GLOBALID_COMPARATOR = 1;
		static const ComparatorType         CENTROID_COMPARATOR = 2;
		static const ComparatorType         HIERARCHY_COMPARATOR= 3;
		static const ComparatorType         HANDLE_COMPARATOR   = 4;
		//typedef INMOST_DATA_BULK_TYPE       ExchangeType;
		//static const ExchangeType           SYNC_ELEMENTS_NONE  = 0; //elements are not synced in parallel set
		//static const ExchangeType           SYNC_ELEMENTS_ALL   = 1; //all elements are synchronized in parallel set
		//static const ExchangeType           SYNC_ELEMENTS_SHARED= 2; //only shared elements are present
		                                    ElementSet          () : Element() {}
		                                    ElementSet          (Mesh * m, HandleType h) : Element(m,h) {}
		                                    ElementSet          (Mesh * m, HandleType * h) : Element(m,h) {}
		                                    ElementSet          (const ElementSet & other) : Element(other) {}
		__INLINE ElementSet &               operator =          (ElementSet const & other) {Element::operator =(other); return *this;}
		__INLINE ElementSet *               operator ->         () {return this;}
		__INLINE const ElementSet *         operator ->         () const {return this;}
		__INLINE ElementSet &               self                ()              {return *this;}
		__INLINE const ElementSet &         self                () const        {return *this;}
	private:
		/// Compute union of processors on elements.
		/// @param procs Set of processors to fill.
		/// @param dir Directon for tree traversal, 1 - downwards, 2, upwards, 3 - both directions
		void                                CollectProcessors   (std::set<Storage::integer> & procs, char dir) const;
		/// Fill Mesh::SendtoTag() to send sets to processors from procs set that currently don't have this set
		/// @param procs Set of processors to send the set to.
		/// @param dir Directon for tree traversal, 1 - downwards, 2, upwards, 3 - both directions
		void                                SetSendTo(std::set<Storage::integer> & procs, char dir) const;
		void                                SetSendTo(std::vector<Storage::integer> & procs, char dir) const;
		void                                SetSendTo(Storage::integer_array procs, char dir) const;
	public:
		/// Get name of the set
		std::string                 GetName() const;
		/// Retrieve parent of the set
		ElementSet                  GetParent() const;
		/// Retrieve sibling set of the set, this will be next child for the parent
		ElementSet                  GetSibling() const;
		/// Retrieve child set of the set.
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
		/// Retrieve number of stored handles, including invalid.
		/// If you want to get number of valid elements use ElementSet::Size
		enumerator                  nbHandles() const;
		/// Retrieve position after last sorted element. This does not correctly
		/// represent total number of sorted elements, since some of them may be deleted
		enumerator                  nbSorted() const;
		/// Retrieve all elements by type
		enumerator                  nbAdjElements(ElementType etype) const;
		enumerator                  nbAdjElements(ElementType etype, MarkerType select, bool invert = false) const;
		/// Retrieve all elements by type
		ElementArray<Element>       getAdjElements(ElementType etype) const;
		ElementArray<Element>       getAdjElements(ElementType etype, MarkerType select, bool invert = false) const;
		
		/// Retrieve only nodes
		ElementArray<Node>          getNodes() const;
		ElementArray<Node>          getNodes(MarkerType select, bool invert = false) const;
		/// Retrieve only edges
		ElementArray<Edge>          getEdges() const;
		ElementArray<Edge>          getEdges(MarkerType select, bool invert = false) const;
		/// Retrieve only faces
		ElementArray<Face>          getFaces() const;
		ElementArray<Face>          getFaces(MarkerType select, bool invert = false) const;
		/// Retrieve only cells
		ElementArray<Cell>          getCells() const;
		ElementArray<Cell>          getCells(MarkerType select, bool invert = false) const;
		/// Put one element without checking of the existance of duplicate.
		/// For sorted set new element will appear at unsorted part.
		void                        PutElement(HandleType e) const;
		/// Put one element without checking of the existance of duplicate.
		/// For sorted set new element will appear at unsorted part
		void                        PutElement(const Storage & e) const {PutElement(e->GetHandle());}
		/// Put multiple handles without checking of the existance of duplicate
		void                        PutElements(const HandleType * handles, enumerator num) const;
		/// Put multiple handles of the other set without checking of the existance of duplicate
		void                        PutElements(const ElementSet & other) const {PutElements(other->getHandles(),other->nbHandles());}
		/// Put multiple handles without checking
		template<typename EType>
		void                        PutElements(const ElementArray<EType> & elems) const {PutElements(elems.data(),static_cast<enumerator>(elems.size()));}
		/// Put one element with checking of the existance of duplicate.
		/// Preserves order for sorted set, thus may be expensive.
		void                        AddElement(HandleType e) const;
		/// Put one element with checking of the existance of duplicate.
		/// Preserves order for sorted set, thus may be expensive.
		void                        AddElement(const Storage & e) const {AddElement(e->GetHandle());}
		/// Add multiple elements with checking of the existance of duplicate.
		/// Preserves order for sorted set, thus may be expensive.
		/// This will also remove any duplicates in unsorted part of the set.
		/// If you inserted duplicated elements through PutElements into previously sorted array
		/// then this operation does not guarantee that those duplications will be stored.
		/// \todo Recheck usage of markers.
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

		/// Compute and store difference with raw handles.
		/// \todo
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
		/// Performs sort of the set of elements. If the set was previously sorted but
		/// have unsorted part, then unsorted part will be sorted and two parts will be merged.
		/// If you need all the set to be resorted (for example in the case global ids were changed)
		/// then invoke SortSet with UNSORTED_COMPARATOR first and then with needed comparator.
		///
		/// Internally it uses:
		/// - std::sort with Mesh::CentroidComparator for CENTROID_COMPARATOR
		/// - std::sort with Mesh::IerarhyComparator  for  HIERARCHY_COMPARATOR
		/// - radix sort starting from certain size for   GLOBALID_COMPARATOR
		/// - radix sort starting from certain size for     HANDLE_COMPARATOR
		/// - only changes state, no sort performed for   UNSORTED_COMPARATOR 
		///
		/// After the set was sorted all the invalid handles should appear at the end of the set
		/// and then removed from array, so it will implicitly work as ReorderEmpty function.
		/// No checks that elements are hidden performed (Maybe this checks should be done in comparators)
		/// In the case you formed the set by running over all mesh elements from NODE to ESET in 
		/// increasing order then your set will be automatically sorted by handles, in this case
		/// you are encouraged to override Mesh::SetComparatorTag with HANDLE_COMPARATOR on the set
		/// without invoking SortSet, so that SortSet does not redundantly perform any work.
		/// You are encouraged to do so even if you are not going to use this information -
		/// some internal algorithms may benefit from it.
		///
		/// \todo
		/// !TODO 52 - check radix sort on big endian computer
		/// @param comp one of the comparators from description.
		/// @see Mesh::SetComparatorTag
		void SortSet(ComparatorType comp) const;
		/// Sets the synchronization regime for set elements
		/// @param comp one of the synchronization types from description.
		//void SetExchange(ExchangeType comp) const;
		/// Perform binary search by global id. In set sorted with GLOBALID_COMPARATOR in O(log(n)) time
		/// otherwise search needs O(n) comparisons
		Element FindElementByGlobalID(integer global_id) const;
		/// Perform binary search by centroid. In set sorted with CENTROID_COMPARATOR in O(log(n)) time
		/// otherwise search needs O(n) comparisons
		Element FindElementByCentroid(real * centroid) const;
		/// Performs linear search in unsorted set. 
		/// - binary search by handle in set sorted with      HANDLE_COMPARATOR
		/// - binary search by centroid in set sorted with  CENTROID_COMPARATOR
		/// - binary search by global id in set sorted with GLOBALID_COMPARATOR
		/// - binary search by ierarhy in set sorted with    HIERARCHY_COMPARATOR (may be too expensive)
		///
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
		void SetPrivateMarkerElements(MarkerType m, ElementType etype = ESET|CELL|FACE|EDGE|NODE) const;
		/// Remove markers from all the elements of given type
		void RemMarkerElements(MarkerType m, ElementType etype = ESET|CELL|FACE|EDGE|NODE) const;
		void RemPrivateMarkerElements(MarkerType m, ElementType etype = ESET|CELL|FACE|EDGE|NODE) const;
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
		/// Retrieve current set comparator
		ComparatorType GetComparator() const;
		/// Retrieve current set exchange type
		//ExchangeType GetExchange() const;
		/// Compact holes in inner representation
		void ReorderEmpty() const;
		/// Is there any elements in the set
		bool Empty() const;
		/// Get total number of elements
		enumerator Size() const;
		/// Remove all elements, clear all data, removes sorted marker
		void Clear();
		/// Remove the set and resolve it's ierarchical structure.
		/// This will not remove childrens of the tree.
		/// To remove set as a tree, see ElementSet::DeleteSetTree.
		/// @return Same as Element::Delete.
		bool DeleteSet();
		/// Remove the set and all it's children.
		/// @return Same as Element::Delete.
		bool DeleteSetTree();
		/// Asks all the elements to be sent to other processors.
		/// Call ExchangeMarked afterwards.
		/// @see Mesh::ExchangeMarked
		void SynchronizeSetElements();
		/// Asks all the elements of ghost sets to be sent to the owner processors.
		/// Call ExchangeMarked afterwards.
		/// @see Mesh::ExchangeMarked
		void SynchronizeSetElementsWithOwner();
		/// Asks all the children to be sent to other processors.
		/// Call ExchangeMarked afterwards.
		/// @see Mesh::ExchangeMarked
		void SynchronizeSetChildren();
		/// Asks all the parents upwards to be sent to other processors.
		/// This function does not ask children of the parents to be synchronized,
		/// for this traverse to the upppermost parent and call ElementSet::SynchronizeSetChildren
		/// Call ExchangeMarked afterwards.
		/// @see Mesh::ExchangeMarked
		void SynchronizeSetParents();
	};

	__INLINE const ElementSet & InvalidElementSet() {static ElementSet ret(NULL,InvalidHandle()); return ret;}
	
	class Mesh : public TagManager, public Storage //implemented in mesh.cpp
	{
	public:
#if defined(CHECKS_MARKERS)
		bool check_shared_mrk, check_private_mrk;
#endif
		enum MeshState {Serial, Parallel};
		typedef chunk_array<integer,chunk_bits_empty>               
											empty_container;
		typedef chunk_array<integer,chunk_bits_elems>               
											links_container;
		//typedef std::vector<integer>                                empty_container;
		//typedef std::vector<integer>                                links_container;
		typedef TagManager::sparse_sub_type                         
											sparse_type;
		typedef TagManager::sparse_sub_record                       
											sparse_rec;
		typedef sparse_type::size_type                              
											senum;
	private:
		std::string                         name;
		std::map<std::string,HandleType>    set_search;
		real                                epsilon;
		empty_container                     empty_space[6];
		empty_container                     empty_links[6];
		links_container                     links[6];
		Tag                                 tag_global_id;
		Tag                                 tag_coords;
		Tag                                 tag_low_conn;
		Tag                                 tag_high_conn;
		Tag                                 tag_markers;
		Tag *                               tag_private_markers;
		Tag                                 tag_geom_type;
		Tag                                 tag_setname;
		Tag                                 tag_setcomparator;
		//Tag                                 tag_setexchange;
		MeshState                           m_state;
		integer                             dim;
		HandleType                          last_created;
		integer								hidden_count[6];
		integer								hidden_count_zero[6];
	private:
		INMOST_DATA_BIG_ENUM_TYPE           parallel_mesh_unique_id;
		INMOST_MPI_Comm                     comm;
        //INMOST_MPI_Group                    group;
		Tag                                 tag_shared;
		Tag                                 tag_owner;
		Tag                                 tag_layers;
		Tag                                 tag_bridge;
		Tag                                 tag_redistribute;
	private:
        //double dist(Cell a, Cell b);
		void AllocatePrivateMarkers();
		void DeallocatePrivateMarkers();
		__INLINE static sparse_rec          mkrec               (const Tag & t) {sparse_rec ret; ret.tag = t.mem; ret.rec = NULL; return ret;}
		__INLINE sparse_type const &        MGetSparseLink      (integer etypenum, integer ID) const {assert(links[etypenum][ID] != -1); return GetSparseData(etypenum,links[etypenum][ID]);}
		__INLINE sparse_type &              MGetSparseLink      (integer etypenum, integer ID) {assert(links[etypenum][ID] != -1); return GetSparseData(etypenum,links[etypenum][ID]);}
		__INLINE sparse_type const &        MGetSparseLink      (HandleType h) const {return MGetSparseLink(GetHandleElementNum(h),GetHandleID(h));}
		__INLINE sparse_type &              MGetSparseLink      (HandleType h) {return MGetSparseLink(GetHandleElementNum(h),GetHandleID(h));}
		__INLINE const void *               MGetSparseLink      (HandleType h, const Tag & t) const {sparse_type const & s = MGetSparseLink(GetHandleElementNum(h),GetHandleID(h)); for(senum i = 0; i < s.size(); ++i) if( s[i].tag == t.mem ) return s[i].rec; return NULL;}
		__INLINE void * &                   MGetSparseLink      (HandleType h, const Tag & t) {sparse_type & s = MGetSparseLink(GetHandleElementNum(h),GetHandleID(h)); for(senum i = 0; i < s.size(); ++i) if( s[i].tag == t.mem ) return s[i].rec; s.push_back(mkrec(t)); return s.back().rec;}
		__INLINE const void *               MGetDenseLink       (integer n, integer ID, const Tag & t) const {assert(links[n][ID] != -1); return &(GetDenseData(t.GetPositionByDim(n))[links[n][ID]]);}
		__INLINE void *                     MGetDenseLink       (integer n, integer ID, const Tag & t) {assert(links[n][ID] != -1); return &(GetDenseData(t.GetPositionByDim(n))[links[n][ID]]);}
		__INLINE const void *               MGetDenseLink       (HandleType h, const Tag & t) const {return MGetDenseLink(GetHandleElementNum(h),GetHandleID(h),t);}
		__INLINE void *                     MGetDenseLink       (HandleType h, const Tag & t) {return MGetDenseLink(GetHandleElementNum(h),GetHandleID(h),t);}
		__INLINE const void *               MGetLink            (HandleType h, const Tag & t) const {if( !t.isSparseByDim(GetHandleElementNum(h)) ) return MGetDenseLink(h,t); else return MGetSparseLink(h,t);}
		__INLINE void *                     MGetLink            (HandleType h, const Tag & t) {if( !t.isSparseByDim(GetHandleElementNum(h)) ) return MGetDenseLink(h,t); else {void * & q = MGetSparseLink(h,t); if( q == NULL ) AllocateSparseData(q,t); return q;}}
		void                                AllocateSparseData  (void * & q, const Tag & t);
		void                                Init                (std::string name);
	public:
		Tag                                 tag_sendto;
		Tag                                 tag_processors;
		/// Go through all elements and detect presence of prescribed element in
		/// any reference data tag.
		void                                ReportConnection(HandleType h);
		/// Test whether global identificator was set on certain type of elements.
		/// This function does not validate correctness of global identificators.
		/// @param type Single type of elements on which to test presence of global identificators.
		/// @return Returns true if global identificators are present on provided type of elements.
		bool                                HaveGlobalID        (ElementType types) const;
		/// Remove all data and all elements from the mesh
		/// Reset geometry service and topology check flags
		void                                Clear               ();
		/// For parmetis
		/// return total number in bytes of occupied memory by element and its data
		enumerator                          MemoryUsage         (HandleType h);
		                                    Mesh                ();
											Mesh                (std::string name);
		                                    Mesh                (const Mesh & other);
		Mesh &                              operator =          (Mesh const & other);
		virtual                             ~Mesh               ();
		/// Allocate a new marker.
		/// Assert will fire in debug mode (NDEBUG not set) if you run out of space for markers, in this case you either
		/// use too many markers, then you can just increase MarkerFields variable (increasing by 1 gives you 8 more markers)
		/// or you forget to release markers after you use them.
		///
		/// In release mode (NDEBUG is set) if you run out of space for markers function will return InvalidMarker()
		/// @return New marker or InvalidMarker(), see description.
		MarkerType                          CreateMarker         ();
		MarkerType                          CreatePrivateMarker  ();
		/// Release marker back for reuse.
		/// This function will only notice mesh that the marker is free to be reused, this will not clear markers
		/// from elements. Before releasing the marker you should ensure that all the marker is removed from all the elements.
		/// Since many inner algorithms use markers, not properly clearing markers from elements may lead to undefined behavior.
		///
		/// Since it is expensive to check asserts will fire in debug mode (NDEBUG not set) only if you define CHECKS_MARKERS in inmost_common.h,
		/// no checks will be performed in release mode(NDEBUG is set).
		/// @param n Byte position and byte bit mask.
		/// @param cleanup Elements on which marker should be removed.
		void                                ReleaseMarker       (MarkerType n, ElementType cleanup = NONE);
		void                                ReleasePrivateMarker(MarkerType n, ElementType cleanup = NONE);
		/// Set tolerance for coordinates comparison. This tolerance is used in comparators 
		/// when two meshes are merged during loading, in ResolveShared to check that nodes on different processors match 
		/// and in UnpackElementsData
		/// @param e small real value
		__INLINE void                       SetEpsilon          (real e) {epsilon = e;}
		/// Retrieve tolerance for coordinates comparison.
		/// @return real value
		__INLINE real                       GetEpsilon          () const {return epsilon;}
		/// Set number of dimensions for mesh. It sets the size for number of internally stored coordinates.
		/// Size of the array returned by Node::Coords will match this number.
		/// When you change number of dimensions and there are nodes with bigger number of dimensions then
		/// first dim coordinates are copied and the rest are dropped.
		/// @see Node::Coords
		void                                SetDimensions       (integer dim);
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
		//__INLINE const Tag &              SetExchangeTag     () const {return tag_setexchange;}
		/// Don't put this shortcut to any function directly, as it creates tag inside
		/// assign to other object of type Tag and put this object to functions.
		__INLINE Tag                      RedistributeTag    () {return CreateTag("TEMPORARY_NEW_OWNER",DATA_INTEGER,CELL,NONE,1);}
		/// Create the tag by name, type and size.
		/// You cannot create the tag with the same name and different type or size.
		///
		/// You may make subsequent calls to the function with the same name, same type and same size 
		/// (or undefined size, then it will be deduced) but different selection of element.
		///
		/// The following recomendation is for future:
		///   When you create data of variable size, every array of data for every element will be empty at first, so
		///   before accessing it through mechanism for single-valued data (Real, Integer, Bulk, Reference) 
		///   you should first resize arrays otherwise your program very likely to be halted by segmentation fault. 
		///   For this case arrays are not allocated automatically from performance considerations.
		///
		/// @param name name of the tag
		/// @param dtype type of the tag
		/// @param etype the selection of elements on which the data of the tag is defined, you may use bitwise or operations 
		/// to define tag on multiple types of elements, example CELL | FACE
		/// @param sparse the selection of elements from etype on which the tag is sparse, for example, if you know that the data is used 
		/// on all cells and only on boundary faces, then you may should set etype = CELL | FACE and sparse = FACE
		/// @param size size of associated data
		/// @return returns the tag that represents the data
		Tag                               CreateTag          (std::string name, DataType dtype, ElementType etype,ElementType sparse, INMOST_DATA_ENUM_TYPE size = ENUMUNDEF);
		/// Remove the data that is represented by the tag from elements of selected type.
		/// @param tag tag that indicates the data
		/// @param mask the selection of the elements on which the data should be removed, may be set by bitwise or operation
		/// @return returns tag that returns false on isValid() request if all the data was removed and the tag is no more occupied, otherwise returns the tag
		/// @see Tag::isValid
		Tag                               DeleteTag          (Tag tag, ElementType mask = NODE | EDGE | FACE | CELL | ESET | MESH);
		/// Returns the first tag defined on the mesh.
		/// For safety while iterating through tags you should check for validity of the tag
		/// @return first tag
		__INLINE iteratorTag              BeginTag           () {return tags.begin(); }
		/// Retrieve the total number of valid tags.
		/// @return returns the number of valid tags
		__INLINE enumerator               NumberOfTags       () const { return static_cast<INMOST_DATA_ENUM_TYPE>(tags.size()); }
		/// Returns the indicator for loop to end iteration over tags.
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
		std::pair<Cell,bool>              CreateCell         (const ElementArray<Face> & faces);//, const ElementArray<Node> & suggest_nodes_order = ElementArray<Node>(NULL));
		std::pair<Cell,bool>              CreateCell         (const ElementArray<Node> & c_f_nodes, const integer * c_f_numnodes, integer num_c_faces);//, 
		                                                      //const ElementArray<Node> & suggest_nodes_order = ElementArray<Node>(NULL));
		std::pair<Cell,bool>              CreateCell         (const ElementArray<Node> & c_nodes, const integer * c_f_nodeinds, const integer * c_f_numnodes, integer num_c_faces);//, 
		                                                      //const ElementArray<Node> & suggest_nodes_order = ElementArray<Node>(NULL));
		std::pair<ElementSet,bool>        CreateSet          (std::string name);
		/// Same as Mesh::CreateSet without checking existance of the set
		std::pair<ElementSet,bool>        CreateSetUnique    (std::string name);
		/// Retrieve set by name.
		/// @param name set name
		/// @return set whose name match or InvalidHandle()
		ElementSet                        GetSet             (std::string name);
		/// Retrieve all the sets whose names start with given prefix
		ElementArray<ElementSet>          GetSetsByPrefix    (std::string prefix);
		HandleType                        LastCreated        () const {return last_created;}

		bool                              isValidHandleRange (HandleType h) const; //for asserts
		bool                              isValidElementNum  (integer etypenum, integer lid) const {return links[etypenum][lid] != -1;}
		bool                              isValidElement     (ElementType etype, integer lid) const {return links[ElementNum(etype)][lid] != -1;}
		bool                              isValidCell        (integer lid) const {return links[ElementNum(CELL)][lid] != -1;}
		bool                              isValidFace        (integer lid) const {return links[ElementNum(FACE)][lid] != -1;}
		bool                              isValidEdge        (integer lid) const {return links[ElementNum(EDGE)][lid] != -1;}
		bool                              isValidNode        (integer lid) const {return links[ElementNum(NODE)][lid] != -1;}
		bool                              isValidElementSet  (integer lid) const {return links[ElementNum(ESET)][lid] != -1;}
		bool                              isValidElement     (HandleType h) const {return  isValidHandleRange(h) && isValidElementNum(GetHandleElementNum(h),GetHandleID(h));}
		/// Retrieve upper adjacent that is shared by multiple lower adjacencies.
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
		///
		/// If you know that data is certanly dense and fixed or variable on elements you access then it is
		/// faster to use specialized variants of this function.
		///
		/// Reference to the data is guaranteed to be valid during mesh modification.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		/// @see Mesh::RealDF
		/// @see Mesh::RealDV
		real      &                         Real                (HandleType h, const Tag & tag);
		/// Returns a reference to inner memory location of the first element of the array of integer values.
		/// Future recomendation:
		///   If variable size array was not allocated then this function will generate segmentation fault.
		///
		/// If you know that data is certanly dense and fixed or variable on elements you access then it is
		/// faster to use specialized variants of this function.
		///
		/// Reference to the data is guaranteed to be valid during mesh modification.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		integer   &                         Integer             (HandleType h, const Tag & tag);
		/// Returns a reference in inner representation to the first element of array of bytes.
		/// Future recomendation:
		///   If variable size array was not allocated then this function will generate segmentation fault.
		///
		/// If you know that data is certanly dense and fixed or variable on elements you access then it is
		/// faster to use specialized variants of this function.
		///
		/// Reference to the data is guaranteed to be valid during mesh modification.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		bulk      &                         Bulk                (HandleType h, const Tag & tag);
		/// Returns a reference in an inner representation to the first element of array of element handles.
		/// Future recomendation:
		///   If variable size array was not allocated then this function will generate segmentation fault.
		///
		/// If you know that data is certanly dense and fixed or variable on elements you access then it is
		/// faster to use specialized variants of this function.
		///
		/// Reference to the data is guaranteed to be valid during mesh modification.
		///
		/// Using handle you can construct objects of type Storage, Element, Node, Edge,
		/// Face, Cell, ElementSet by calling their constructor with pointer to mesh and handle
		/// as arguments.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		reference &                         Reference           (HandleType h, const Tag & tag);
    /// Returns a reference in an inner representation to the first element of array of element remote handles.
		/// Future recomendation:
		///   If variable size array was not allocated then this function will generate segmentation fault.
		///
		/// If you know that data is certanly dense and fixed or variable on elements you access then it is
		/// faster to use specialized variants of this function.
		///
		/// Reference to the data is guaranteed to be valid during mesh modification.
		///
		/// Using remote handle you can construct objects of type Element with the function
		/// MakeElement or MakeElementRef.
		///
		/// @param h element handle
		/// @param tag tag that represents data
    /// @see MakeElement
    /// @see MakeElementRef
		remote_reference &                  RemoteReference     (HandleType h, const Tag & tag);
		/// Returns an array of real values.
		/// If you know that data is certanly dense on elements you access then it is faster to use
		/// variants of this function with hint data structure.
		///
		/// Array data structure is guaranteed to be valid during mesh modification.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		real_array                          RealArray           (HandleType h, const Tag & tag);
		/// Returns an array of integer values.
		/// If you know that data is certanly dense and fixed or variable on elements you access then it is
		/// faster to use specialized variants of this function.
		/// variants of this function with hint data structure.
		///
		/// Array data structure is guaranteed to be valid during mesh modification.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		integer_array                       IntegerArray        (HandleType h, const Tag & tag);
		/// Returns an array of bytes.
		///
		/// If you know that data is certanly sparse or dense on elements you access then it is faster to use
		/// variants of this function with hint data structure.
		///
		/// Array data structure is guaranteed to be valid during mesh modification.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		bulk_array                          BulkArray           (HandleType h, const Tag & tag);
		/// Returns an array of element handles.
		///
		/// If you know that data is certanly sparse or dense on elements you access then it is faster to use
		/// variants of this function with hint data structure.
		///
		/// The class reference_array that is used to represent array of elements stores handles inside but
		/// accessing them through square scopes [] or by arrow -> in iterator will automatically
		/// form object of type Element. If you are not sure that stored handles are valid, you
		/// should either check that by Element::isValid (involves deep check) or test handle
		/// against InvalidHandle (simple check). To obtain handle you may use reference_array::at function or
		/// dereference * operator for iterator. If you need custom object like Node, Edge, Face, Cell
		/// or ElementSet you may use Element::getAsNode, Element::getAsEdge, Element::getAsFace,
		/// Element::getAsCell and Element::getAsSet functions.
		///
		/// Array data structure is guaranteed to be valid during mesh modification. If you delete
		/// elements by Mesh::Delete or Element::Delete all the references are also will be valid
		/// and reverted to InvalidHandle on Mesh::ApplyModification. If you use Mesh::Destroy to
		/// delete mesh elements or you delete elements not within modification state then references 
		/// may become either invalid but not testable against InvalidHandle (situation may be tested 
		/// by Element::isValid or Mesh::isValidHandle) or reference may be reused by another element. 
		/// If you mix deletion and construction of elements then there is no way to resolve this situation,
		/// except if you have created only one element, then it may be retrieved by Mesh::LastHandle.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		/// @see InvalidHandle
		/// @see Element::isValid
		/// @see Element::getAsNode
		/// @see Element::getAsEdge
		/// @see Element::getAsFace
		/// @see Element::getAsCell
		/// @see Element::getAsSet
		reference_array                     ReferenceArray      (HandleType h, const Tag & tag);
    /// Returns an array of element remote handles.
		///
		/// If you know that data is certanly sparse or dense on elements you access then it is faster to use
		/// variants of this function with hint data structure.
		///
		/// The class reference_array that is used to represent array of elements stores remote handles 
    /// inside but accessing them through square scopes [] or by arrow -> in iterator will 
		/// automatically form an object of type Element. If you are not sure that stored handles are valid, you
		/// should either check that by Element::isValid (involves deep check) or test handle
		/// against InvalidHandle (simple check). To obtain remote handle you may use remote_reference_array::at 
    /// function or dereference * operator for iterator. If you need custom object like Node, Edge, Face, Cell
		/// or ElementSet you may use Element::getAsNode, Element::getAsEdge, Element::getAsFace,
		/// Element::getAsCell and Element::getAsSet functions.
		///
		/// Array data structure is guaranteed to be valid during mesh modification. If you delete
		/// elements by Mesh::Delete or Element::Delete all the references are also will be valid
		/// and reverted to InvalidHandle on Mesh::ApplyModification. If you use Mesh::Destroy to
		/// delete mesh elements or you delete elements not within modification state then references 
		/// may become either invalid but not testable against InvalidHandle (situation may be tested 
		/// by Element::isValid or Mesh::isValidHandle) or reference may be reused by another element. 
		/// If you mix deletion and construction of elements then there is no way to resolve this situation,
		/// except if you have created only one element, then it may be retrieved by Mesh::LastHandle.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		/// @see InvalidHandle
		/// @see Element::isValid
		/// @see Element::getAsNode
		/// @see Element::getAsEdge
		/// @see Element::getAsFace
		/// @see Element::getAsCell
		/// @see Element::getAsSet
		remote_reference_array              RemoteReferenceArray(HandleType h, const Tag & tag);
		/// Returns a reference to inner memory location of the first element of the array of real values.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// @param h element handle
		/// @param tag tag that represents data
		real      &                         RealDF              (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_REAL     ); return static_cast<real     *>(MGetDenseLink(h,tag))[0];}
		/// Returns a reference to inner memory location of the first element of the array of integer values.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		integer   &                         IntegerDF           (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_INTEGER  ); return static_cast<integer  *>(MGetDenseLink(h,tag))[0];}
		/// Returns a reference in dense array to the first element of constant size array of bytes.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		bulk      &                         BulkDF              (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_BULK     ); return static_cast<bulk     *>(MGetDenseLink(h,tag))[0];}
		/// Returns a reference in dense array to the first element of constant size array of element handles.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// Using handle you can construct objects of type Storage, Element, Node, Edge,
		/// Face, Cell, ElementSet by calling their constructor with pointer to mesh and handle
		/// as arguments.
		///
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		reference &                         ReferenceDF         (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_REFERENCE); return static_cast<reference*>(MGetDenseLink(h,tag))[0];}
    /// Returns a reference in dense array to the first element of constant size array of element remote handles.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// Using handle you can construct objects of type Storage, Element, Node, Edge,
		/// Face, Cell, ElementSet by calling their constructor with pointer to mesh and handle
		/// as arguments.
		///
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		remote_reference &                  RemoteReferenceDF   (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_REMOTE_REFERENCE); return static_cast<remote_reference*>(MGetDenseLink(h,tag))[0];}
		/// Returns an array of real values in dense array.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// Note that as array is fixed you shouldn't use any functions that alter size of the array
		/// as resize, erase, insert, you may use replace if initial and final size will match,
		/// in debug mode assert will fire if you try to do this in release (NDEBUG is set) it will lead to segfault.
		///
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		real_array                          RealArrayDF         (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_REAL     ); return real_array     (static_cast<real     *>(MGetDenseLink(h,tag)),tag.GetSize());}
		/// Returns an array of integer values in dense array.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set), likely to result in segfault.
		///
		/// Note that as array is fixed you shouldn't use any functions that alter size of the array
		/// as resize, erase, insert, you may use replace if initial and final size will match,
		/// in debug mode assert will fire if you try to do this in release (NDEBUG is set) it will lead to segfault.
		///
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		integer_array                       IntegerArrayDF      (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_INTEGER  ); return integer_array  (static_cast<integer  *>(MGetDenseLink(h,tag)),tag.GetSize());}
		/// Returns an array of bytes in dense array.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set), likely to result in segfault.
		///
		/// Note that as array is fixed you shouldn't use any functions that alter size of the array
		/// as resize, erase, insert, you may use replace if initial and final size will match,
		/// in debug mode assert will fire if you try to do this in release (NDEBUG is set) it will lead to segfault.
		///
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		bulk_array                          BulkArrayDF         (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_BULK     ); return bulk_array     (static_cast<bulk     *>(MGetDenseLink(h,tag)),tag.GetSize());}
		/// Returns an array of element handles in dense array.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set), likely to result in segfault.
		///
		/// Note that as array is fixed you shouldn't use any functions that alter size of the array
		/// as resize, erase, insert, you may use replace if initial and final size will match,
		/// in debug mode assert will fire if you try to do this in release (NDEBUG is set) it will lead to segfault.
		///
		/// The class reference_array that is used to represent array of elements stores handles inside but
		/// accessing them through square scopes [] or by arrow -> in iterator will automatically
		/// form object of type Element. If you are not sure that stored handles are valid, you
		/// should either check that by Element::isValid (involves deep check) or test handle
		/// against InvalidHandle (simple check). To obtain handle you may use reference_array::at function or
		/// dereference * operator for iterator. If you need custom object like Node, Edge, Face, Cell
		/// or ElementSet you may use Element::getAsNode, Element::getAsEdge, Element::getAsFace,
		/// Element::getAsCell and Element::getAsSet functions.
		///
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		/// @see InvalidHandle
		/// @see Element::isValid
		/// @see Element::getAsNode
		/// @see Element::getAsEdge
		/// @see Element::getAsFace
		/// @see Element::getAsCell
		/// @see Element::getAsSet
		reference_array                     ReferenceArrayDF    (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_REFERENCE); return reference_array(this,static_cast<reference*>(MGetDenseLink(h,tag)),tag.GetSize());}
    /// Returns an array of element remote handles in dense array.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and fixed is incorrect,
		/// no checks performed in release mode (NDEBUG is set), likely to result in segfault.
		///
		/// Note that as array is fixed you shouldn't use any functions that alter size of the array
		/// as resize, erase, insert, you may use replace if initial and final size will match,
		/// in debug mode assert will fire if you try to do this in release (NDEBUG is set) it will lead to segfault.
		///
		/// The class remote_reference_array that is used to represent array of elements stores handles inside but
		/// accessing them through square scopes [] or by arrow -> in iterator will automatically
		/// form object of type Element. If you are not sure that stored handles are valid, you
		/// should either check that by Element::isValid (involves deep check) or test handle
		/// against InvalidHandle (simple check). To obtain handle you may use reference_array::at function or
		/// dereference * operator for iterator. If you need custom object like Node, Edge, Face, Cell
		/// or ElementSet you may use Element::getAsNode, Element::getAsEdge, Element::getAsFace,
		/// Element::getAsCell and Element::getAsSet functions.
		///
		/// @param h element handle
		/// @param tag tag that represents dense data of fixed size on given handle
		/// @see InvalidHandle
		/// @see Element::isValid
		/// @see Element::getAsNode
		/// @see Element::getAsEdge
		/// @see Element::getAsFace
		/// @see Element::getAsCell
		/// @see Element::getAsSet
		remote_reference_array              RemoteReferenceArrayDF(HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_REMOTE_REFERENCE); return remote_reference_array(static_cast<remote_reference*>(MGetDenseLink(h,tag)),tag.GetSize());}
		/// Returns a reference in dense array to the first element of variable size array of real values.
		/// Future recomendation:
		///    If array was not allocated (resized) then this function will generate segmentation fault.
		///
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// @param h element handle
		/// @param tag tag that represents data
		real      &                         RealDV              (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_REAL     ); return static_cast<inner_real_array     *>(MGetDenseLink(h,tag))->at_safe(0);}
		/// Returns a reference in dense array to the first element of variable size array of integer values.
		/// Future recomendation:
		///    If array was not allocated then this function will generate segmentation fault.
		///
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// @param h element handle
		/// @param tag tag that represents data
		integer   &                         IntegerDV           (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_INTEGER  ); return static_cast<inner_integer_array  *>(MGetDenseLink(h,tag))->at_safe(0);}
		/// Returns a reference in dense array to the first element of variable size array of bytes.
		/// Future recomendation:
		///    If array was not allocated then this function will generate segmentation fault.
		///
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// @param h element handle
		/// @param tag tag that represents data
		/// @see TagDenseVariable
		bulk      &                         BulkDV              (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_BULK     ); return static_cast<inner_bulk_array     *>(MGetDenseLink(h,tag))->at_safe(0);}
		/// Returns a reference in dense array to the first element of variable size array of element handles.
		/// Future recomendation:
		///    If array was not allocated then this function will generate segmentation fault.
		///
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// Using handle you can construct objects of type Storage, Element, Node, Edge,
		/// Face, Cell, ElementSet by calling their constructor with pointer to mesh and handle
		/// as arguments.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		reference &                         ReferenceDV         (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_REFERENCE); return static_cast<inner_reference_array*>(MGetDenseLink(h,tag))->at_safe(0);}
    /// Returns a reference in dense array to the first element of variable size array of element remote handles.
		/// Future recomendation:
		///    If array was not allocated then this function will generate segmentation fault.
		///
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// Using handle you can construct objects of type Storage, Element, Node, Edge,
		/// Face, Cell, ElementSet by calling their constructor with pointer to mesh and handle
		/// as arguments.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		remote_reference &                  RemoteReferenceDV   (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_REMOTE_REFERENCE); return static_cast<inner_remote_reference_array*>(MGetDenseLink(h,tag))->at_safe(0);}
		/// Returns an array of real values in dense array of variable size.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// @param h element handle
		/// @param tag tag that represents data
		real_array                          RealArrayDV         (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_REAL     ); return real_array     (*static_cast<inner_real_array     *>(MGetDenseLink(h,tag)));}
		/// Returns an array of integer values in dense array of variable size.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// @param h element handle
		/// @param tag tag that represents data
		integer_array                       IntegerArrayDV      (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_INTEGER  ); return integer_array  (*static_cast<inner_integer_array  *>(MGetDenseLink(h,tag)));}
		/// Returns an array of bytes in dense array of variable size.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// @param h element handle
		/// @param tag tag that represents data
		bulk_array                          BulkArrayDV         (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_BULK     ); return bulk_array     (*static_cast<inner_bulk_array     *>(MGetDenseLink(h,tag)));}
		/// Returns an array of element handles in dense array of variable size.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// The class reference_array that is used to represent array of elements stores handles inside but
		/// accessing them through square scopes [] or by arrow -> in iterator will automatically
		/// form object of type Element. If you are not sure that stored handles are valid, you
		/// should either check that by Element::isValid (involves deep check) or test handle
		/// against InvalidHandle (simple check). To obtain handle you may use reference_array::at function or
		/// dereference * operator for iterator. If you need custom object like Node, Edge, Face, Cell
		/// or ElementSet you may use Element::getAsNode, Element::getAsEdge, Element::getAsFace,
		/// Element::getAsCell and Element::getAsSet functions.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		/// @see InvalidHandle
		/// @see Element::isValid
		/// @see Element::getAsNode
		/// @see Element::getAsEdge
		/// @see Element::getAsFace
		/// @see Element::getAsCell
		/// @see Element::getAsSet
		reference_array                     ReferenceArrayDV    (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_REFERENCE); return reference_array(this,*static_cast<inner_reference_array*>(MGetDenseLink(h,tag)));}
		/// Returns an array of element remote handles in dense array of variable size.
		/// If you don't know any hint information about tag data you should not use this function.
		///
		/// Asserts will fire in debug mode if assumption that data is dense and variable is incorrect,
		/// no checks performed in release mode (NDEBUG is set).
		///
		/// The class remote_reference_array that is used to represent array of elements stores remote handles inside but
		/// accessing them through square scopes [] or by arrow -> in iterator will automatically
		/// form object of type Element. If you are not sure that stored handles are valid, you
		/// should either check that by Element::isValid (involves deep check) or test handle
		/// against InvalidHandle (simple check). To obtain handle you may use reference_array::at function or
		/// dereference * operator for iterator. If you need custom object like Node, Edge, Face, Cell
		/// or ElementSet you may use Element::getAsNode, Element::getAsEdge, Element::getAsFace,
		/// Element::getAsCell and Element::getAsSet functions.
		///
		/// @param h element handle
		/// @param tag tag that represents data
		/// @see InvalidHandle
		/// @see Element::isValid
		/// @see Element::getAsNode
		/// @see Element::getAsEdge
		/// @see Element::getAsFace
		/// @see Element::getAsCell
		/// @see Element::getAsSet
		remote_reference_array              RemoteReferenceArrayDV(HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_REMOTE_REFERENCE); return remote_reference_array(*static_cast<inner_remote_reference_array*>(MGetDenseLink(h,tag)));}

#if defined(USE_AUTODIFF)
		var      &                          Variable            (HandleType h, const Tag & tag);
		var      &                          VariableDF          (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_VARIABLE); return static_cast<var     *>(MGetDenseLink(h,tag))[0];}
		var      &                          VariableDV          (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_VARIABLE); return static_cast<inner_variable_array     *>(MGetDenseLink(h,tag))->at_safe(0);}
		var_array                           VariableArray       (HandleType h, const Tag & tag);
		var_array                           VariableArrayDF     (HandleType h, const Tag & tag) {AssertsDF(h,tag,DATA_VARIABLE); return var_array(static_cast<var *>(MGetDenseLink(h,tag)),tag.GetSize());}
		var_array                           VariableArrayDV     (HandleType h, const Tag & tag) {AssertsDV(h,tag,DATA_VARIABLE); return var_array(*static_cast<inner_variable_array*>(MGetDenseLink(h,tag)));}
#endif
		/// Set a marker on the element represented by handle.
		/// @param h element handle
		/// @param n stores byte number and byte bit mask that represent marker
		void                              SetMarker          (HandleType h,MarkerType n)  {assert(!isPrivate(n)); static_cast<bulk *>(MGetDenseLink(h,MarkersTag()))[n >> MarkerShift] |= static_cast<bulk>(n & MarkerMask);}
		void                              SetPrivateMarker   (HandleType h,MarkerType n);
		/// Set a marker on the set of handles.
		/// @param h set of handles
		/// @param n number of handles
		/// @param m stores byte number and byte bit mask that represent marker
		/// @see Mesh::SetMarker
		void                              SetMarkerArray     (const HandleType * h, enumerator n, MarkerType m) {for(enumerator i = 0; i < n; ++i) if( h[i] != InvalidHandle() )SetMarker(h[i],m);}
		void                              SetPrivateMarkerArray     (const HandleType * h, enumerator n, MarkerType m) {for(enumerator i = 0; i < n; ++i) if( h[i] != InvalidHandle() )SetPrivateMarker(h[i],m);}
		/// Check whether the marker is set one the element.
		/// @param h element handle
		/// @param n stores byte number and byte bit mask that represent marker
		bool                              GetMarker          (HandleType h,MarkerType n) const {assert(!isPrivate(n)); return (static_cast<const bulk *>(MGetDenseLink(h,MarkersTag()))[n >> MarkerShift] & static_cast<bulk>(n & MarkerMask)) != 0;}
		bool                              GetPrivateMarker   (HandleType h,MarkerType n) const;
		/// Remove the marker from the element.
		/// @param h element handle
		/// @param n stores byte number and byte bit mask that represent marker
		void                              RemMarker          (HandleType h,MarkerType n) {assert(!isPrivate(n)); static_cast<bulk *>(MGetDenseLink(h,MarkersTag()))[n >> MarkerShift] &= ~static_cast<bulk>(n & MarkerMask);}
		void                              RemPrivateMarker   (HandleType h,MarkerType n);
		/// Remove the marker from the set of handles.
		/// @param h set of handles
		/// @param n number of handles
		/// @param m stores byte number and byte bit mask that represent marker
		/// @see Mesh::RemMarker
		void                              RemMarkerArray     (const HandleType * h, enumerator n, MarkerType m) {for(enumerator i = 0; i < n; ++i) if( h[i] != InvalidHandle() ) RemMarker(h[i],m);}
		void                              RemPrivateMarkerArray     (const HandleType * h, enumerator n, MarkerType m) {for(enumerator i = 0; i < n; ++i) if( h[i] != InvalidHandle() ) RemPrivateMarker(h[i],m);}
		/// Remove all the markers from the element
		void                              ClearMarkerSpace   (HandleType h);
		/// Get a copy of the bytes that store markers on the element.
		/// @param h element handle
		/// @param copy storage to put data to
		void                              GetMarkerSpace     (HandleType h,bulk copy  [MarkerFields]) const;
		/// Overwrite the bytes that store markers on the element.
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
		/// - For NODE this returns edges that are connected to this node;
		/// - For EDGE this returns faces that are connected to this edge;
		/// - For FACE this returns cells that are connected to this face
		/// - For CELL this returns nodes of the cell
		/// - For ESET first three entries are parent, sibling, child then records represent empty positions in LowConn
		///
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
		/// Check that upper adjacencies are stored
		ElementType                       HaveUpperAdjacencies() const;
		/// Delete all upper adjacencies, access to HighConn should fire assertion and retrival of upper adjacencies is no longer valid
		void                              RemoveUpperAdjacencies(ElementType mask = (NODE|EDGE|FACE));
		/// Restore all upper adjacencies
		void                              RestoreUpperAdjacencies(ElementType mask = (NODE|EDGE|FACE));
		/// Check that upper adjacencies are stored
		ElementType                       HaveLowerAdjacencies() const;
		/// Delete all upper adjacencies, access to HighConn should fire assertion and retrival of upper adjacencies is no longer valid
		void                              RemoveLowerAdjacencies(ElementType mask = (EDGE|FACE|CELL));
		/// Restore all upper adjacencies
		void                              RestoreLowerAdjacencies(ElementType mask = (EDGE|FACE|CELL));
		/// Access directly higher order adjacencies of current element without right of modification.
		Element::adj_type const&          HighConn           (HandleType h) const {return *static_cast<const inner_reference_array*>(MGetDenseLink(h,HighConnTag()));}
		/// Access directly lower order adjacencies of current element with right of modification.
		/// This function gives direct access to elements and allows you to overwrite handles.
		/// If you are going to overwrite them then read recomendations in description for HighConn function.
		/// - For NODE this returns cells that are connected to this node;
		/// - For EDGE this returns nodes of the edge
		/// - For FACE this returns edges of the face
		/// - For CELL this returns faces of the cell
		/// - For ESET handles of the elements that this set contain
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
		/// @see Tag::GetSize
		INMOST_DATA_ENUM_TYPE             GetDataSize        (HandleType h,const Tag & tag) const; //For DATA_BULK return number of bytes, otherwise return the length of array
		/// Return the size of the structure in bytes required to represent the data on current element.
		/// This is equal to GetDataSize times Tag::GetBytesSize for all the data types,
		/// except for DATA_VARIABLE, that requires a larger structure to accomodate derivatives.
		/// @param h handle of element
		/// @param tag tag that represents the data
		INMOST_DATA_ENUM_TYPE             GetDataCapacity    (HandleType h,const Tag & tag) const;
		/// Returns the number of bytes in data used for given type of tag.
		/// Trivial for all the types except DATA_VARIABLE.
		INMOST_DATA_ENUM_TYPE             GetDataCapacity    (const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size, const Tag & tag) const;
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
		/// Removes data of variable size, clears to zero data of fixed size.
		/// @param data link to data
		/// @param tag tag that indicates the data
		void                              DelDenseData       (void * data,const Tag & tag);
		/// Removes data of variable size and sparse tag data.
		/// @param h handle to the element
		/// @param tag tag that indicates the data
		bool                              DelSparseData      (HandleType h,const Tag & tag);
		/// Check whether data is present on given element.
		/// Always returns true for dense tag data.
		/// @param h handle to the element
		/// @param tag tag that indicates data
		/// @return true, if data exists otherwise false
		bool                              HaveData           (HandleType h,const Tag & tag) const;
		
		Element::GeometricType            GetGeometricType   (HandleType h) const {return static_cast<const bulk *>(MGetDenseLink(h,GeomTypeTag()))[0];}
		void                              SetGeometricType   (HandleType h, Element::GeometricType type) {static_cast<bulk *>(MGetDenseLink(h,GeomTypeTag()))[0] = type;}
		/// Retrieve global id of the element with right of modification (dangerous to modify).
		/// Run AssignGlobalID so that tag is automatically allocated and shortcut is set within mesh,
		/// otherwise tag is not created and call will fail.
		/// @param h handle of the element
		/// @return global id
		/// @see Mesh::AssignGlobalID
		integer &                         GlobalID           (HandleType h) {return static_cast<integer *>(MGetDenseLink(h,GlobalIDTag()))[0];}
		/// Retrieve global id of the element without right of modification.
		/// Run AssignGlobalID so that tag is automatically allocated and shortcut is set within mesh,
		/// otherwise tag is not created and call will fail.
		/// @param h handle of the element
		/// @return global id
		/// @see Mesh::AssignGlobalID
		integer                           GlobalID           (HandleType h) const {return static_cast<const integer *>(MGetDenseLink(h,GlobalIDTag()))[0];}
		/// Retrieve position of the data position of current element. After ReorderEmpty
		/// this number is guaranteed to be between 0 and NumberOf(type of element)
		/// @param h handle of the element
		/// @return local id of data
		integer                        DataLocalID           (HandleType h) const {return static_cast<integer>(links[GetHandleElementNum(h)][GetHandleID(h)]);}
		/// Retrieve parallel status of the element.
		/// If mesh is in Serial state then call always returns Element::Owned.
		/// otherwise it will return:
		/// -  Element::Owned  if there is a single copy of the element on the current processor
		/// -  Element::Shared if the main copy of the element is located on the current processor
		/// -  Element::Ghost  if current processor stores dependent copy of the element
		/// @param h handle of the element
		/// @return Element::Status, see function description
		Element::Status                   GetStatus          (HandleType h) const { if( SharedTag().isValid() ) return static_cast<const bulk *>(MGetDenseLink(h,SharedTag()))[0]; return Element::Owned;}
		/// Set parallel status of the element.
		/// If mesh is in Serial state then call will fire asserts in debug mode and segfault in relese mode.
		///
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
		/// Show element from mesh.
		/// @param h handle of the element
		/// @return true then element was recovered
		/// @see Mesh::Hide
		bool                              Show               (HandleType h);
		/// This function will hide element in modification state (between BeginModification and EndModification)
		/// or call Destroy in non-modification state.
		/// @param h handle of the element
		/// @return  if true then element was deleted, otherwise it was hidden
		/// @see Mesh::Hide
		/// @see Mesh::BeginModification
		/// @see Mesh::EndModification
		/// @see Mesh::Destroy
		bool                              Delete             (HandleType h);
		/// Check whether element is hidden.
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
		typedef std::vector<INMOST_MPI_Request>                       exch_reqs_type;
		typedef struct exch_recv_reqs_t
		{
			std::vector<unsigned> buf;
			std::vector<unsigned> cnt;
			std::vector<INMOST_MPI_Request> requests;
			unsigned count()
			{
				unsigned ret = 0;
				for(size_t i = 0; i < cnt.size(); ++i)
					ret += cnt[i];
				return ret;
			}
			void clear()
			{
				buf.clear();
				buf.push_back(0);
				cnt.clear();
				requests.clear();
			}
		} exch_recv_reqs_type;
		
		class exchange_data
		{
		public:
			exch_reqs_type send_reqs;
			exch_recv_reqs_type recv_reqs;
			exch_buffer_type send_buffers, recv_buffers;
		};
	private:
		class Random // random generator to provide tag for communication
		{
		private: unsigned int n,a,c,m;
		public:
			Random(unsigned int seed = 50);
			//~ Random(const Random & other);
			unsigned int Number();
		} randomizer;
	public:
		class elements_by_type
		{
		private:
			element_set container[6];
		public:
			elements_by_type() {}
			elements_by_type(const elements_by_type & other) {for(int i = 0; i < 6; i++) container[i] = other.container[i];}
			~elements_by_type(){}
			element_set & operator [](int i){ return container[i]; }
			const element_set & operator [](int i) const { return container[i]; }
			bool empty() {bool ret = true; for(int i = 0; i < 6 && ret; i++) ret &= container[i].empty(); return ret;}
			unsigned size() {unsigned ret = 0; for(int i = 0; i < 6; ++i) ret += (unsigned)container[i].size(); return ret;}
		};
		typedef std::map<int, elements_by_type > parallel_storage;
	
		typedef std::map<int, elements_by_type >               proc_elements_by_type;
	private:
#if defined(USE_PARALLEL_STORAGE)
		parallel_storage                    shared_elements;
		parallel_storage                    ghost_elements;
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		int                                 num_exchanges;
	public:
		std::fstream                        out_time;
	private:
		int                                 tab;
		int                                 func_id;
#endif
#if defined(USE_MPI_P2P)
		INMOST_MPI_Win                      window;
		INMOST_DATA_BIG_ENUM_TYPE *         shared_space;
#endif
		int                                 parallel_file_strategy;
	private:
		void                              ListTags           (tag_set & list);
		std::vector<int>                  ComputeSharedProcs (const parallel_storage & from, const parallel_storage & to);
		std::vector<int>                  ComputeSharedProcs (ElementType etype);
		proc_elements                     ComputeSharedSkinSet(ElementType bridge, MarkerType marker = 0);
		void                              PackElementsGather (elements_by_type & selems, const elements_by_type & elements, int destination, ElementType mask, MarkerType select, const tag_set & tag_list, bool force_send, bool send_links_to_owner, bool pack_by_gids);
		void                              PackElementsEnumerate   (elements_by_type & selems, TagInteger pack_position);
		void                              PackElementsUnenumerate (elements_by_type & selems, TagInteger pack_position);
		void                              PackTagData        (const Tag & tag, const elements_by_type & elements, int destination, ElementType mask, MarkerType select, buffer_type & buffer, TagInteger pack_position);
		void                              UnpackTagData      (const Tag & tag, const elements_by_type & elements, int source, ElementType mask, MarkerType select, buffer_type & buffer, size_t & buffer_position, ReduceOperation op, const elements_by_type & unpack_elements);//, proc_elements_by_type * send_elements = NULL);
		void                              PackElementsData   (elements_by_type & input, buffer_type & buffer, int destination, const tag_set & tag_list,TagInteger pack_position, bool pack_gids);
		void                              UnpackElementsData (elements_by_type & output, buffer_type & buffer, int source, size_t & buffer_position, tag_set & tag_list);
		void                              PrepareReceiveInner(Prepare todo, exch_buffer_type & send_bufs, exch_buffer_type & recv_bufs);
		void                              ExchangeDataInnerBegin(const tag_set & tag, const parallel_storage & from, const parallel_storage & to, ElementType mask, MarkerType select, exchange_data & storage);
		void                              ExchangeDataInnerEnd(const tag_set & tag, const parallel_storage & from, const parallel_storage & to, ElementType mask, MarkerType select, ReduceOperation op, exchange_data & storage);
		void                              ExchangeBuffersInner(exch_buffer_type & send_bufs, exch_buffer_type & recv_bufs, exch_reqs_type & send_reqs, exch_recv_reqs_type & recv_reqs);
		std::vector<int>                  FinishRequests     (exch_recv_reqs_type & recv_reqs);
		void                              SortParallelStorage(parallel_storage & ghost, parallel_storage & shared,ElementType mask);
		void                              ReportParallelStorage();
		void                              GatherParallelStorage(parallel_storage & ghost, parallel_storage & shared, ElementType mask);
		void                              InformElementsOwners(proc_elements_by_type & send_elements, exchange_data & storage);
	public:
		HandleType                        FindSharedGhost(ElementType etype, Storage::integer global_id, int source_proc, int owner_proc);
#if defined(USE_PARALLEL_WRITE_TIME)	
		//this part is needed to test parallel performance
		void                              Enter              ();
		void                              Exit               ();
		int &                             GetFuncID          () {return func_id;}
		std::fstream &                    GetStream          ();
		std::ostream &                    WriteTab           (std::ostream & f);
		void                              FinalizeFile       ();
		static void                       AtExit             (void);
#endif
		void                              ClearFile          ();
		/// Initial initialization. Calls MPI_Initialize, if MPI was not initialized
		/// it is necessery to invoke this function if you plan to use any parallel algorithms
		/// Accepts arguments passed to console aplication or NULL
		/// @param argc number of arguments for command line
		/// @param argv strings of arguments of command line
		static void                       Initialize         (int * argc, char *** argv);
		/// Finalizes operation with MPI, recomended to call, otherwise MPI may produce warnings
		static void                       Finalize           ();
		/// Set parallel strategy for inner communications.
		/// There are three possible scenaries in parallel communication associated in
		/// accordance to enum Prepare structre:
		///   1. The communicating processors and sizes of the messages are known apriori
		///   2. UnknownSize: Communicating processors are known but sizes are unknown
		///   3. unknownSource: Communicationg processors are unknown
		///
		/// Currently with UnknownSize it will run following algorithm
		/// none for strategy 0, following for strategies 1 and 2:
		///   1. Post asynchronous receive with MPI_Irecv of size of buffer to be sent
		///   2. Post asynchronous send with MPI_Isend for required corresponding receive buffer size
		///   3. Wait for all asynchronous operations by MPI_Waitall
		///
		/// With UnknownSource there are two options depending from the USE_MPI_P2P define.
		/// If USE_MPI_P2P is defined then MPI-2 api operations will be used
		///   1. MPI_Win_fence to start operations
		///   2. MPI_Put from specially allocated memory location to remote processor 
		///      of array size is performed
		///   3. MPI_Win_fence to stop operations
		///
		/// if USE_MPI_P2P not set then MPI-1 api will be used
		///   1. MPI_Allgather for number of processors to which current processor wants to send data
		///   2. MPI_Allgatherv for sizes and destinations for each processors
		///
		/// Initially it was intended to mainly use MPI-2 functionality for both scenarios but
		/// it was realized that there is no availible hardware on which MPI-2 functionalty
		/// performs match better then MPI-1 counterparts, especially in the case of UnknownSize.
		/// Probably this happens due to lack of support of RDMA operations.
		/// If you believe it will be better to use MPI-2 in both cases you are free to uncomment
		/// definition of PREFFER_MPI_P2P in inmost_common.h then MPI-2 will be used in both scenaries.
		/// These algorthms above are implemented in Mesh::ExchangeBufferInner.
		///
		/// After first problem was resolved following strategies are availible for main communication:
		/// strategy = 0
		///   1. Asyncronous send of data by MPI_Isend;
		///   2. Check incoming messages by MPI_Probe;
		///   3. Check incoming message size by MPI_Get_size;
		///   4. Allocation of buffers of required size;
		///   5. Asynchronous receive of data
		///   6. MPI_Waitsome to copy received results to buffers
		///
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
		///   1. Post asynchronous receive of data by MPI_Irecv
		///   2. Post asynchronous send of data by MPI_Isend
		///   3. MPI_Waitsome for any received data
		///
		/// True asynchronous behavior is reproduced by breaking 1-2) and 3)
		/// strategy = 2
		///   1. Post asynchronous receive of data by MPI_Irecv
		///   2. Set MPI_Barrier to ensure that all the receives were properly set by the time
		///   2. Perform direct send of data by MPI_Irsend
		///   3. MPI_Waitsome for any received data
		///
		/// For asynchronous communication algorithm is broken into 1-3) and 4) which is fairly
		/// asynchronous. The only provisional benefit it may have on machines with small memory 
		/// since it should bypass any allocation of buffers for sent and received data by MPI and 
		/// probably not perform any randezvous communication to ensure data allocation.
		/// But MPI_Barrier looks like elephant here.
		///
		/// Algorithms above are implemented in Mesh::ExchangeBuffersInner
		/// @see Mesh::PrepareReceiveInner
		/// @see Mesh::ExchangeBuffersInner
		//no longer used
		//void                              SetParallelStrategy(int strategy){assert( !(strategy < 0 || strategy > 3) ); parallel_strategy = strategy;}
		/// Retrieve currently set parallel strategy.
		/// @see Mesh::SetParallelStrategy
		//no longer used
		//int                               GetParallelStrategy() {return parallel_strategy;}
		/// This strategy correspond only to internal ".pmf" mesh format.
		/// There are two availible strategies for ".pmf" files loading and saving:
		///
		/// strategy = 0
		/// - on save 
		///   1. every processor gather local data to some buffer
		///   2. MPI_Gather to obtain sizes of data among processors
		///   3. MPI_Gatherv to obtain the whole data on zeros processor
		///   4. the first processor writes all the data to disk by std::fstream
		/// - on load
		///   1. first processors reads the whole file by std::fstream
		///   2. MPI_Scatter distributes block sizes among processors
		///   3. MPI_Scatterv distributes blocks among processors
		///   4. Each processor parses it's block
		///
		/// This strategy requires one processor to hold all the data, which
		/// is quite bad for large files. New strategy may be created from this one in future
		/// when each processors consequently obtain access to the file using std::fstream and
		/// writes the data.
		///
		/// strategy = 1
		/// - on save it will perform:
		///   1. MPI_Gather to obtain sizes of data among processors on zeroth processor
		///   2. MPI_File_open to get parallel handle for the file
		///   3. MPI_File_write_shared called by processor with zeroth rank to write header
		///   4. MPI_File_write_ordered to write contents of individual data
		///   5. MPI_File_close to close parallel file handle
		/// - on load it will perform
		///   1. MPI_File_open to open the file in parallel
		///   2. MPI_File_read_shared to get contents of header on zeroth processor
		///   3. MPI_Scatter to distribute block sizes among processors
		///   4. MPI_File_read_ordered to obtain contents
		///   5. MPI_File_close to close parallel file handle
		///
		///  Availible only when USE_MPI_P2P is set because it rely on MPI-2 api that begins with MPI_File_xxx
		///  some MPI-1 standards contain this api as extension.
		///
		/// The strategy 1 appeared to be considerably slower on INM cluster then strategy 0, this may
		/// happen due to lack of read-write devices that able to work in parallel. On IBM Bluegene/p
		/// strategy 1 was not working due to same old problem with shared file pointers in their MPI realization
		void                              SetParallelFileStrategy(int strategy){assert( !(strategy < 0 || strategy > 1) ); parallel_file_strategy = strategy;}
		/// Retrieve currently set parallel strategy for ".pmf" files
		/// @see Mesh::GetParallelStrategy
		int                               GetParallelFileStrategy() const {return parallel_file_strategy;}
		/// Get rank of current processor
		int                               GetProcessorRank   () const;
		/// Get number of processors
		int                               GetProcessorsNumber() const;
		/// Get rank of current processor in shared environment (OpenMP)
		int                               GetLocalProcessorRank() const;
		/// Get number of processors in shared environment (OpenMP)
		int                               GetLocalProcessorNumber() const;
		/// Retrieve MPI communicator
		INMOST_MPI_Comm                   GetCommunicator    () const;
        /// Retrieve MPI group corresponding to the communicator
        INMOST_MPI_Group                  GetGroup           () const;
		/// Set MPI communicator
		void                              SetCommunicator    (INMOST_MPI_Comm _comm);
		/// Find elements that are common between processors.
		void                              ResolveShared      (bool only_new = false);
		/// Find sets that are common between processors.
		void                              ResolveSets        ();
		/// Delete all the ghost cells.
		void                              RemoveGhost        (MarkerType marker = 0);
		/// Delete some ghost cells provided in array.
		/// Non-ghost elements will also be deleted.
		///
		/// This algorithm will properly communicate between processors so that 
		/// parallel states of deleted elements properly updated on remote processors.
		/// Internally function invokes Destroy function, not Delete, which hides elements during modification state,
		/// currently it is not expected that any parallel algorithms will be performed between BeginModification and EndModification
		/// since modification may break parallel state though it is never checked whether the mesh is in the modification state,
		/// so you are free to experiment. This behavior may change in future.
		///
		/// Collective point-2-point.
		///
		/// @param ghost array of handles
		/// @param num number of handles
		/// \todo
		///  1. Currently request for deletion of elements of lower level then cell will be simply ignored, ensure
		///     in future that algorithm will properly rise deletion data from lower to upper adjacencies 
		///     to delete all the upper adjacencies that depend on deleted lower adjacencies
		void                              RemoveGhostElements(const HandleType * ghost, enumerator num);
		template<typename EType>
		void                              RemoveGhostElements(const ElementArray<EType> & ghost) {RemoveGhostElements(ghost.data(),static_cast<enumerator>(ghost.size()));}
		void                              RemoveGhostElements(const ElementSet & ghost) {RemoveGhostElements(ghost.getHandles(),ghost.nbHandles());}
		/// Assign unique numbers to elements. Internally this will create Mesh::GlobalIDTag and 
		/// make call to Element::GlobalID and Mesh::GlobalID functions valid. Internally this 
		/// will also set have_global_id variable that will indicate that all the comparisons 
		/// in parallel algorithms should be performed using global identificators instead 
		/// of centroids which is faster.
		///
		/// \todo
		///     1. invoking function before loading mesh will not renew global identificators after load
		///        but would not unset have_global_id either. There are probably too many places when
		///        global ids may become invalid but no flag will be set. It may be benefitial to set
		///        such flags along with updating geometrical data which seems to be maintained fairly well
		///        during mesh modification
		void                              AssignGlobalID     (ElementType mask);
		/// Update data from Shared elements to Ghost elements. For backward direction please see Mesh::ReduceData.
		/// If you have a tag of DATA_BULK type and you store your own custom data structure in it, it is highly
		/// recomended that you provide MPI information about your structure through Tag::SetBulkDataType,
		/// this would not do any difference on homogeneous architecture, but it may help you save a lot of 
		/// time and nerves in heterogeneous parallel environment.
		///
		/// \todo
		/// see TODO in Mesh::ReduceData
		///
		/// Blocking, Collective point-2-point
		///
		/// @param tag tag that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @see Mesh::ReduceData
		void                              ExchangeData       (const Tag & tag, ElementType mask, MarkerType select = 0);
		/// Start asynchronous synchronization of data.
		/// You should define object of type exchange_data that will hold temporary buffers for data.
		/// every Mesh::ExchangeDataBegin should be matched with Mesh::ExchangeDataEnd with the same
		/// exchange_data object. After matching Mesh::ExchangeDataEnd the exchange_data object may be reused
		/// If you will go out of the scope where exchange_data object was defined it will be deallocated
		/// and may result in segmentation fault.
		///
		/// You should also never put the same exchange_data object to any other Mesh::ExchangeDataBegin or
		/// Mesh::ReduceDataBegin, until matching Mesh::ExchangeDataEnd because it may override or reallocate 
		/// buffers, internally used by MPI and remote processor will receive garbage instead of data.
		///
		/// Nonblocking, Collective point-2-point
		///
		/// @param tag tag that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		/// @see Mesh::ExchangeDataEnd
		void                              ExchangeDataBegin  (const Tag & tag, ElementType mask, MarkerType select, exchange_data & storage);
		/// Complete asynchronous synchronization of data.
		///
		/// Blocking
		///
		/// @param tag tag that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		/// @see Mesh::ExchangeDataBegin
		void                              ExchangeDataEnd    (const Tag & tag, ElementType mask, MarkerType select, exchange_data & storage);
		/// This function will perform exchange of multiple data tags.
		///
		/// Blocking, Collective point-2-point
		///
		/// @param tags multiple tags that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		void                              ExchangeData       (const tag_set & tags, ElementType mask, MarkerType select = 0);
		/// This function will initialize exchange of multiple data tags.
		/// Using this function may lead to good overlapping between communication and computation.
		///
		/// Nonblocking, Collective point-2-point
		///
		/// @param tags multiple tags that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		void                              ExchangeDataBegin  (const tag_set & tags, ElementType mask, MarkerType select, exchange_data & storage);
		/// This function will finalize exchange of multiple data tags.
		///
		/// Blocking
		///
		/// @param tags multiple tags that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		void                              ExchangeDataEnd    (const tag_set & tags, ElementType mask, MarkerType select, exchange_data & storage);
		/// \brief Accumulation of data from ghost elements to shared elements.
		///
		/// Accumulation is performed based on user-provided function. 
		/// When processor - owner of the element receives data addressed
		/// to this element then it calls user-defined op function. Since there may be multiple ghost
		/// elements per one shared element, function may be called multiple times.
		///
		/// Several examples of reduction functions may be found within the mesh_parallel.cpp source.
		///
		/// Remember that the result will be up to date only on the owner processor of the element. You will have
		/// to run Mesh::ExchangeData to make the data up to date among all of the processors.
		///
		/// If you have a tag of DATA_BULK type and you store your own custom data structure in it, it is highly
		/// recomended that you provide MPI information about your structure through Tag::SetBulkDataType,
		/// this would not do any difference on homogeneous architecture, but it may help you save a lot of 
		/// time and nerves in heterogeneous parallel environment.
		///
		/// Exchanging tags of DATA_REFERNCE is not implemented, TODO 14.
		/// \todo
		///    1. Exchanging DATA_REFERENCE,DATA_REMOTE_REFERENCE tags not implemented, this is due to the absence of any conclusion
		///    -  on how it should behave:
		///         either only search within elements owned by the other processor and
		///         establish references and set InvalidHandle() to elements that are not found (fairly easy,
		///         will involve search operations to check against owned elements for similar entry, efficient
		///         implementation will require bounding search trees (see TODO 49);
		///    -  or:
		///         send all the referenced elements through PackElementsData and establish all the links within
		///         elements reproduced by UnpackElementsData (UnpackElementsData calls UnpackTagData with set
		///         of unpacked elements using which it will be very comfortable to establish references on
		///         remote processor). Drawback is that exchanging laplacian operator in such a manner should
		///         result in the whole grid being shared among all the processors.
		///
		/// Blocking, Collective point-2-point
		///
		/// @param tag tag that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param op user-defined operation on received data
		/// @see Mesh::ExchangeData
		void                              ReduceData         (const Tag & tag, ElementType mask, MarkerType select, ReduceOperation op );
		/// This function intializes data reduction.
		/// Read recomendations about exchange_storage object in Mesh::ExchangeDataBegin
		///
		/// Nonblocking, Collective point-2-point
		///
		/// @param tag tag that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		/// @see Mesh::ReduceData
		/// @see Mesh::ExchangeDataBegin
		void                              ReduceDataBegin    (const Tag & tag, ElementType mask, MarkerType select, exchange_data & storage);
		/// This function completes data reduction.
		/// Read Mesh::ReduceData for information about op function
		///
		/// Blocking
		///
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
		///
		/// Blocking, collective point-2-point
		///
		/// @param tags multiple tags that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param op user-defined operation on received data
		void                              ReduceData         (const tag_set & tags, ElementType mask, MarkerType select, ReduceOperation op );
		/// This function will initialize reduction of multiple data tags.
		///
		/// Nonblocking, collective point-2-point
		///
		/// @param tags multiple tags that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		void                              ReduceDataBegin    (const tag_set & tags, ElementType mask, MarkerType select, exchange_data & storage);
		/// This function will finalize exchange of multiple data tags.
		///
		/// Blocking
		///
		/// @param tags multiple tags that represents data
		/// @param mask bitwise or of element types
		/// @param select set the marker to filter elements that perform operation, set 0 to select all elements
		/// @param storage buffer that will temporary hold sended data
		/// @param op user-defined operation on received data
		void                              ReduceDataEnd      (const tag_set & tags, ElementType mask, MarkerType select, ReduceOperation op, exchange_data & storage );
		/// This function realizes two algorithms: ghosting of elements and migration of elements.
		/// ghosting:
		///
		///  Creates ghosted elements at other processors prescribed in SendtoTag() performs
		///  sending and receiving of elements and all operations to keep parallel state of the mesh.
		///
		///   Given that all the data was up to date among processors all the data at the end of the
		///   algorithm will be also up to data
		///
		/// migration:
		///
		///  To correctly perform migration of elements it is necessery to set up tags besides SendToTag(),
		///  which indicates to where every element should be sent:
		/// - tag "TEMPORARY_NEW_PROCESSORS", of type DATA_INTEGER with variable size, tells which processors will have copy of the element after migration;
		/// - tag "TEMPORARY_NEW_OWNER", of type DATA_INTEGER of size 1, tells which processor will have main copy of the element.
		///
		/// if there is no current processor in "TEMPORARY_NEW_PROCESSORS", then current processor will remove copy of the element.
		/// All this actions are performed automatically by Mesh::Redistribute based on information provided in Mesh::RedistributeTag
		/// which effectively contains new owner.
		///
		///   Given that all the data was up to date among processors all the data at the end of the algorithm will be also up to data
		///
		/// \todo
		///     1. test halo exchange algorithm (if used then change collective point-2-point to collective)
		///     2. see TODO 2 in Mesh::Redistribute
		///
		/// Collective point-2-point.
		///
		/// @param action
		/// @see Mesh::SendtoTag
		/// @see Mesh::Redistribute
		void                              ExchangeMarked     (enum Action action = AGhost);
		/// Form several layers of ghosted cells that are adjacent through bridge elements to current cells.
		/// This function acceptes any mesh topology, failure on some mesh should be considered a bug and
		/// sample example should be provided for testing purposes.
		///
		/// This function internally calculates layer by layer and invokes ExchangeMarked for each layer,
		/// you can either reproduce the algorithm on your own if you want to bypass the function and
		/// call Mesh::ExchangeMarked directly, but then you will lose optimization in Mesh::Redistribute,
		/// that expects that layers are formed the same way they are formed in Mesh::ExchangeGhost.
		/// Internally it sets up LayersTag and BridgeTag for the mesh with provided values, which
		/// are used by Mesh::Redistribute but you are discouraged to override these tags since using
		/// non-matching algorithms is not tested and should be considered dangerous.
		///
		/// Nevertheless you can use this function first for layers then request any additional ghosted elements
		/// by ExchangeMarked.
		///
		/// Collective point-2-point.
		///
		/// @param layers number of required layers of ghosted elements
		/// @param bridge bitwise mask of elements for which neighbouring cells should be considered a layer
		/// @see Mesh::ExchangeMarked
		/// @see Mesh::Redistribute
		void                              ExchangeGhost      (integer layers, ElementType bridge, MarkerType select = 0);
		/// Migrate all the elements to the new owners prescribed in data corresponding to RedistributeTag.
		/// This will perform all the actions to send mesh elements and data and reproduce new mesh partitions
		/// on remote elements and correctly resolve parallel state of the mesh. If you have priviously
		/// prescribed number of layers through ExchangeGhost, then minimal number of actions will be performed
		/// to reproduce layers of ghosted elements whithout involving removal of all ghosted elements.
		///
		/// Internally function sets up following data on elements using provided information:
		/// - "TEMPORARY_NEW_PROCESSORS" - new set processors that contain copy of the element
		/// - "TEMPORARY_NEW_OWNER"      - new owner for each processor (effectively RedistributeTag)
		///
		/// Action of this function regarding restoration of layers of ghosted elements in the case you have 
		/// modified mesh without involving Mesh::ResolveModification is yet to be tested and should be
		/// considered dangerous.
		///
		/// If you have output from Zoltan or ParMetis  for cells of the mesh
		/// then just write this output to RedistributeTag and call Mesh::Redistribute.
		///
		/// \todo
		///      1. introduce "TEMPORARY_KEEP_GHOSTED" tag that will store processors on which copy of element
		///         should be kept, internally just merge it with "TEMPORARY_NEW_PROCESSORS" tag
		///         this will allow user to control ghosting of certain elements and not to invoke ExchangeMarked
		///         every time after Redistribute.
		///         This is probably already done using Mesh::SendtoTag, because function fills it without
		///         clearing and ExchangeMarked performs initial action based on SendtoTag, it is due to
		///         check that SendtoTag is properly merged with "TEMPORARY_NEW_PROCESSORS" before call to ExchangeMarked
		///         and received elements are not deleted by accident.
		///      2. let user provide any integer tag as input without involving RedistributeTag
		///
		/// Collective point-2-point.
		///
		/// @see Mesh::RedistributeTag
		/// @see Mesh::ExchangeGhost
		void                              Redistribute       ();
		/// Enumerate all elements begining with start and put numeration to data associated with num_tag for all elements with given type mask.
		///
		/// Collective operation.
		///
		/// @param mask bitwise or of types of elements
		/// @param num_tag a tag that is associated with the data
		/// @param start starting value for enumeration
		/// @param define_sparse if true then function will define sparse data on elements that don't have it, otherwise it will skip those elements
		/// @return last value on all the processors
		integer                           Enumerate          (ElementType mask, Tag num_tag, integer start = 0, bool define_sparse = false);
		/// Enumerate all elements begining with start and put numeration to data associated with num_tag.
		///
		/// Collective operation.
		///
		/// @param h array of handles
		/// @param num number of handles
		/// @param num_tag a tag that is associated with the data
		/// @param start starting value for enumeration
		/// @param define_sparse if true then function will define sparse data on elements that don't have it, otherwise it will skip those elements
		/// @return last value on all the processors
		integer                           Enumerate          (const HandleType * h, enumerator num, const Tag & num_tag, integer start = 0, bool define_sparse = true);
		/// Enumerate all elements begining with start and put numeration to data associated with num_tag.
		///
		/// Collective operation.
		///
		/// @param elements array of elements
		/// @param num_tag a tag that is associated with the data
		/// @param start starting value for enumeration
		/// @param define_sparse if true then function will define sparse data on elements that don't have it, otherwise it will skip those elements
		/// @return last value on all the processors
		template<typename EType>
		integer                           Enumerate          (const ElementArray<EType> & elements, const Tag & num_tag, integer start = 0, bool define_sparse = true) {return Enumerate(elements.data(),static_cast<enumerator>(elements.size()),num_tag,start,define_sparse);}
		/// Enumerate all elements in the set.
		///
		/// Collective operation.
		///
		/// @param set handle of the set
		/// @param num_tag a tag that is associated with the data
		/// @param start starting value for enumeration
		/// @param define_sparse if true then function will define sparse data on elements that don't have it, otherwise it will skip those elements
		/// @return last value on all the processors
		integer                           EnumerateSet       (const ElementSet & set, const Tag & num_tag, integer start = 0, bool define_sparse = true);
		/// Sum of all physical elements, it excludes ghosted copies.
		/// To compute total number including ghosted copies run Integrate(NumberOf(mask))
		///
		/// Collective operation.
		///
		/// @param mask bitwise mask of element types, example: NODE | CELL
		/// @return sum of elements
		/// @see Mesh::NumberOf
		integer                           TotalNumberOf      (ElementType mask);
		/// Integrate real value over all processors.
		///
		/// Collective operation.
		///
		/// @param input Value on current processor
		/// @return Sum over all processors
		real                              Integrate          (real input);
		/// Integrate unsigned integer value over all processors.
		///
		/// Collective operation.
		///
		/// @param input Value on current processor
		/// @return Sum over all processors
		enumerator                        Integrate          (enumerator input);
		/// Integrate integer value over all processors.
		///
		/// Collective operation.
		///
		/// @param input Value on current processor
		/// @return Sum over all processors
		integer                           Integrate          (integer input);
		/// Integrate an array of real values over all processors.
		///
		/// Collective operation.
		///
		/// @param input An array of values on current processor.
		/// @param size Size of the array.
		/// @return sum over all processors.
		void                              Integrate          (real * input, integer size);
		/// Integrate an array of unsigned integer values over all processors.
		///
		/// Collective operation.
		///
		/// @param input An array of values on current processor.
		/// @param size Size of the array.
		/// @return Sum over all processors.
		void                              Integrate          (enumerator * input, integer size);
		/// Integrate an array of integer values over all processors.
		///
		/// Collective operation.
		///
		/// @param input An array of values on current processor.
		/// @param size Size of the array.
		/// @return Sum over all processors.
		void                              Integrate          (integer * input, integer size);
		/// Integrate data corresponding to tag between all processors.
		/// Elements without the data defined on them or when entry not present will be skipped.
		///
		/// Collective operation.
		///
		/// @param t tag that correspond to data to be integrated
		/// @param entry in the array of data
		/// @param mask bitwise or of types of elements on which to integrate
		/// @return sum between all processors
		real                              Integrate          (const Tag & t, enumerator entry, ElementType mask);
		/// Compute sum of integer values for all processors with rank lower then current, excluding current processor.
		///
		/// Collective operation.
		///
		/// @param input value on current processor
		/// @return described sum
		integer                           ExclusiveSum       (integer input);
		enumerator                        ExclusiveSum       (enumerator input);
		
		real                              AggregateMax       (real input);
		integer                           AggregateMax       (integer input);
		enumerator                        AggregateMax       (enumerator input);
		void                              AggregateMax       (real * input, integer size);
		void                              AggregateMax       (integer * input, integer size);
		void                              AggregateMax       (enumerator * input, integer size);

		real                              AggregateMin       (real input);
		integer                           AggregateMin       (integer input);
		enumerator                        AggregateMin       (enumerator input);
		void                              AggregateMin       (real * input, integer size);
		void                              AggregateMin       (integer * input, integer size);
		void                              AggregateMin       (enumerator * input, integer size);
		/// Regather ghosted and shared element sets for data exchange.
		/// This function will be quite useful if you change statuses of elements
		/// or modify mesh on your own bypassing internal algorithms.
		///
		/// No action will be performed if USE_PARALLEL_STORAGE is not set in inmost_common.h,
		/// since all the elements are computed during exchange phase.
		///
		/// Generally this is not needed if you use high-level algorithms for mesh modification
		/// or mesh redistribution.
		///
		/// @param mask Bitwise mask of element types for which to recompute the parallel storage.
		/// @see Mesh::BeginModification
		/// @see Mesh::EndModification
		/// @see Mesh::ExchangeMarked
		/// @see Mesh::RemoveGhostElements
		void                              RecomputeParallelStorage(ElementType mask);
		/// Sort parallel storage. Parallel storage is sorted according to
		/// global identificators or centroids of elements if global identificators
		/// are not availible. If you manually change global identificators or
		/// if global identificators are not availible and coordinates of nodes
		/// change, then you should invoke this function.
		/// You can check presence of global identificators on single
		/// type of elements using function Mesh::HaveGlobalID.
		/// This function is called automatically inside Mesh::AssignGlobalID.
		/// No action will be performed if USE_PARALLEL_STORAGE is not set in inmost_common.h.
		/// @param mask Bitwise mask of element types for which to sort the parallel storage.
		void                              SortParallelStorage(ElementType mask);
		/// Outputs parallel storage into xml log files.
		/// USE_PARALLEL_STORAGE and USE_PARALLEL_WRITE_TIME should be activated in inmost_common.h.
		/// @param mask Bitwise mask of element types for which to log the parallel storage.
		void                              RecordParallelStorage(ElementType mask);
		/// Synchronize bitwise mask of element types between processors.
		///
		/// Collective operation
		///
		/// @param etype bitwise type mask
		/// @return bitwise result among processors
		ElementType                       SynchronizeElementType(ElementType etype);
		/// Syncronize marker on elements between processors using provided operation.
		/// Depending on requested operation following action is performed:
		/// - SYNC_BIT_SET - value on ghost elements is set by value on corresponding shared processors;
		/// - SYNC_BIT_OR  - bitwise OR between values in ghosted and shared elements;
		/// - SYNC_BIT_AND - bitwise AND between values in ghosted and shared elements;
		/// - SYNC_BIT_XOR - bitwise XOR between values in ghosted and shared elements.
		/// @param marker marker to be synchronized
		/// @param mask bitwise or type mask
		/// @param op operation, one of SYNC_BIT_SET, SYNC_BIT_OR, SYNC_BIT_XOR, SYNC_BIT_AND
		void                              SynchronizeMarker  (MarkerType marker, ElementType mask, SyncBitOp op);	
		//for debug
		void                              Barrier            ();
		void                              BeginSequentialCode();
		void                              EndSequentialCode  ();
		//iterator.cpp::::::::::::::::::::::::::::::::::::::::::::::::::
	public:
		Element                           ElementByLocalIDNum(integer etypenum, integer lid) {assert((etypenum < 5 && (lid >= 0 && lid < static_cast<integer>(links[etypenum].size()))) || (etypenum == 5 && lid == 0)); return Element(this,ComposeHandleNum(etypenum,lid));}
		Element                           ElementByLocalID   (ElementType etype, integer lid) {return ElementByLocalIDNum(ElementNum(etype),lid);}
		Element                           ElementByHandle    (HandleType h) {return Element(this,h);}
		
		HandleType                        NextHandle         (HandleType h) const;
		HandleType                        PrevHandle         (HandleType h) const; //returns InvalidHandle() when go beyond first element
		HandleType                        NextHandle         (HandleType h, ElementType mask) const;
		HandleType                        PrevHandle         (HandleType h, ElementType mask) const; //returns InvalidHandle() when go beyond first element
		HandleType                        FirstHandle        () const {return ComposeHandleNum(ElementNum(NODE),0);}
		HandleType                        LastHandle         () const {return ComposeHandleNum(ElementNum(MESH),1);} 
		HandleType                        FirstHandle        (ElementType etype) const {return ComposeHandleNum(ElementNum(etype),0);}
		HandleType                        LastHandle         (ElementType etype) const  {integer num = ElementNum(etype); return ComposeHandleNum(num,static_cast<integer>(links[num].size()));}

		Node                              NodeByLocalID      (integer lid) { assert(lid >= 0 && lid < static_cast<integer>(links[0].size())); return Node(this,ComposeHandleNum(0,lid)); }
		Edge                              EdgeByLocalID      (integer lid) { assert(lid >= 0 && lid < static_cast<integer>(links[1].size())); return Edge(this,ComposeHandleNum(1,lid)); }
		Face                              FaceByLocalID      (integer lid) { assert(lid >= 0 && lid < static_cast<integer>(links[2].size())); return Face(this,ComposeHandleNum(2,lid));}
		Cell                              CellByLocalID      (integer lid) { assert(lid >= 0 && lid < static_cast<integer>(links[3].size())); return Cell(this,ComposeHandleNum(3,lid)); }
		ElementSet                        EsetByLocalID      (integer lid) { assert(lid >= 0 && lid < static_cast<integer>(links[4].size())); return ElementSet(this,ComposeHandleNum(4,lid)); }
		
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

		integer                           LastLocalIDNum     (integer n) const {assert(n >= 0 && n < 6); return n == 5 ? 1 : static_cast<integer>(links[n].size());}

		integer                           NextLocalID        (ElementType etype, integer lid) const {integer q = ElementNum(etype); ++lid; while(lid < static_cast<integer>(links[q].size()) && links[q][lid] == -1) ++lid; return lid;}
		integer                           PrevLocalID        (ElementType etype, integer lid) const {integer q = ElementNum(etype); --lid; while(lid > 0 && links[q][lid] == -1) --lid; return lid;}
		integer                           FirstLocalID       (ElementType etype) const;
		integer                           LastLocalID        (ElementType etype) const {assert(OneType(etype)); return LastLocalIDNum(ElementNum(etype));}

		integer                           NextLocalIDIter    (ElementType etype, integer lid) const;
		integer                           PrevLocalIDIter    (ElementType etype, integer lid) const;
		integer                           FirstLocalIDIter   (ElementType etype) const;
		integer                           LastLocalIDIter    (ElementType etype) const;

		__INLINE integer                  NumberOfSets       () const { return static_cast<integer>(links[4].size() - empty_links[4].size()) - hidden_count[4]; }
		__INLINE integer                  NumberOfCells      () const { return static_cast<integer>(links[3].size() - empty_links[3].size()) - hidden_count[3];}
		__INLINE integer                  NumberOfFaces      () const { return static_cast<integer>(links[2].size() - empty_links[2].size()) - hidden_count[2]; }
		__INLINE integer                  NumberOfEdges      () const { return static_cast<integer>(links[1].size() - empty_links[1].size()) - hidden_count[1]; }
		__INLINE integer                  NumberOfNodes      () const { return static_cast<integer>(links[0].size() - empty_links[0].size()) - hidden_count[0]; }
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
		/// Set file option.
		/// Current availible file options:
		/// - "VERBOSITY"        - Set "2" for progress messages, "1" for reports, "0" for silence
		///
		/// - "VTK_GRID_DIMS"    - Set "2" for two-dimensional vtk grids, "3" for three-dimensional vtk grids
		///                        or "AUTO" for automatic detection.
		/// - "VTK_OUTPUT_FACES" - Set "1" for vtk output of values on faces
		///
		/// - "ECL_SPLIT_GLUED"  - Set "TRUE" to triangulate faces of the blocks that degenerate on three pillars.
		/// - "ECL_PROJECT_PERM" - Set "TRUE" to project permeability tensor from grid block coordinates
		///                        into global coordinates. Otherwise the tensor is considered to be
		///                        defined on global coordinates. Default "FALSE".
		/// - "ECL_COMPUTE_TRAN" - compute and store transmissibilities on the faces using NEWTRAN approach in eclipse 100
		/// - "ECL_DEGENERATE"   - In GRDECL format some active grid block may have zero volume, which
		///                        means there is a fault. Set "TRUE" to introduce a gap between blocks
		///                        that share degenerate active grid block, set to "TRANM" to introduce
		///                        zero transmissibility multiplier and keep the grid connected and "FALSE" 
		///                        to simply keep the grid connected. Default: "TRUE".
		/// - "ECL_TOPOLOGY"     - If "TRUE" checks topology of the grid for errors, this may provide useful
		///                        warnings if layers of the mesh enter each other and the grid cannot be
		///                        considered conformal. Default: "FALSE".
		///   "ECL_PARALLEL_READ"- If "TRUE" then each processor loads part of the eclipse mesh, requires some synchronization.
		///                        Otherwise if "FALSE" then each processor loads entire mesh. Default: "TRUE".
		///   "Tag:TAGNAME"      - Set comma-separated rules for tag with the name TAGNAME, the rules list is:
		///                        nosave - do not save the tag data into files;
		///                        noload - do not load the tag data from files;
		///                        noderivatives - do not save/load the derivatives for data with type DATA_VARIABLE (for .xml and .pmf);
		///                        loadonly - this creates an exclusive list for data to be loaded from files, all other tag names will be ignored;
		///                        saveonly - this creates an exclusive list for data to be saved to files, all other tag names will be ignored.
		///                        Example: mesh->SetFileOption("Tag:PressureGradient","noload,noderivatives");
		///                        the tag with the name PressureGradient will not be loaded from files and 
		///                        when recording the derivaives data will be not saved.
		void         SetFileOption(std::string,std::string);
		/// Get current option corresponding to key.
		/// @param key options for which options should be retrieven
		std::string  GetFileOption(std::string key) const;
		/// Collect file options realated to records Tag:TAGNAME.
		/// @param given option name, such as nosave, noload, noderivatives, loadonly, saveonly
		/// @return a set of tags that has given option
		std::set<std::string> TagOptions(std::string name) const;
		/// Check if tag loading should be skipped.
		bool CheckLoadSkip(std::string name, const std::set<std::string> & noload, const std::set<std::string> & loadonly) const;
		/// Check if tag saving should be skipped.
		bool CheckSaveSkip(std::string name, const std::set<std::string> & noload, const std::set<std::string> & loadonly) const;
		/// Acceptable file formats for reading
		/// - ".vtk"    - legacy vtk format for unstructured grid
		/// - ".pvtk"   - legacy parallel vtk format
		/// - ".gmv"    - format acceptable by general mesh viewer
		/// - ".msh"    - gmsh generator format
		/// - ".grdecl" - eclipse format (under construction)
		/// - ".grid"   - mesh format by Mohammad Karimi-Fard
		/// - ".pmf"    - internal parallel portable binary format, saves all features
		///
		/// @param File path to the file
		/// \todo
		/// 1. When loading mesh with the same tag name but different type or size, load will fail.
		/// 2. When loading tags in internal format should remember definition and sparsity masks
		///    for subsequent data loading. This will cure the case when tags were already priviously defined
		///    on mesh with different masks and data will be red incorrectly.
		void         Load(std::string File);
		void         LoadMSH(std::string File);
		void         LoadECL(std::string File);
		void         LoadXML(std::string File);
		void         LoadPMF(std::string File); 
		void         LoadVTK(std::string File); 
		void         LoadVTU(std::string File);
		void         LoadPVTK(std::string File); 
		void         LoadPVTU(std::string File); 
		void         LoadMKF(std::string File);
		/// Acceptable file formats for writing
		/// - ".vtk"  - legacy vtk format for unstructured grid
		/// - ".pvtk" - legacy parallel vtk format
		/// - ".gmv"  - format acceptable by general mesh viewer
		/// - ".pmf"  - internal parallel portable binary format, saves all features
		///
		/// Remeber: .pmf stores all references to elements. If reference are broken due to mesh modification,
		///          saving or loading such a mesh may lead to seagfault. To automatically maintain correct
		///          references modify mesh using BeginModification, ApplyModification, EndModification
		///
		/// @param File path to the file
		/// \todo
		/// 1. Markers are not saved in internal format due to possible conflict during load.
		void         Save(std::string File);
		void         SaveXML(std::string File);
		void         SavePMF(std::string File);
		void         SaveVTK(std::string File);
		void         SaveVTU(std::string File);
		void         SavePVTK(std::string File);
		void         SavePVTU(std::string File);
		void         SaveGMV(std::string File);
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
	public:
		void                              RepairGeometricTags();
	public:
		bool                              HideGeometricData(GeometricData type, ElementType mask) { remember[type][ElementNum(mask) - 1] = false;  return  remember[type][ElementNum(mask) - 1]; }
		bool                              ShowGeometricData(GeometricData type, ElementType mask) { remember[type][ElementNum(mask) - 1] = true;  return  remember[type][ElementNum(mask) - 1]; }
	public:
		typedef std::map<GeometricData, ElementType> GeomParam;
		// types for MEASURE:     EDGE | FACE | CELL   (length, area, volume)
		// types for CENTROID:    EDGE | FACE | CELL
		// types for BARYCENTER:  EDGE | FACE | CELL
		// types for NORMAL:      FACE | CELL          (may precompute normal for cells in 2d case)
		// types for ORIENTATION: FACE
		/// Marks face with the orientation direction by marker.
		/// If marker is set then face is reversed.
		/// Then all faces are oriented either inside or outside of the cell.
		void                              FacesOrientation(ElementArray<Face> & faces, MarkerType rev);
		bool                              CheckConvexity(const ElementArray<Face> & faces) const;
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
		void                              RecomputeGeometricData(HandleType e); // Update all stored geometric data, runs automatically on element construction
		Element::GeometricType            ComputeGeometricType(ElementType element_type, const HandleType * lower_adjacent, INMOST_DATA_ENUM_TYPE lower_adjacent_size);
		/// Compute node-centered interpolation on 2d face for point.
		/// Point should be inside face or on its boundary.
		/// @param x Interpolation point
		/// @param f Face (should be 2d) which contains point
		/// @param nodes_stencil Vector of pairs (node handle, coefficient) to store interpolation data
		void                              WachspressInterpolation2D           (const real * x, const Face & f, std::map<HandleType,real> & nodes_stencil) const;
		/// Compute node-centered interpolation on 3d cell for point
		/// Point should be inside cell or on its boundary.
		/// @param x Interpolation point
		/// @param c Cell (should be 3d) which contains point
		/// @param nodes_stencil Vector of pairs (node handle, coefficient) to store interpolation data
		void                              WachspressInterpolation3D           (const real * x, const Cell & c, std::map<HandleType,real> & nodes_stencil) const;
		/// Sets marker for all the faces that have only one neighbouring cell, works correctly in parallel environment.
		/// @param boundary_marker Non-private marker that will indicate boundary faces.
		void                              MarkBoundaryFaces(MarkerType boundary_marker);
		/// This function should be used to detect normal inversion on ghost interfaces
		/// with respect to normal orientation on owner of the interface.
		/// Due to automatic control over normal orientation in the grid, it may
		/// happen that some ghost faces have different orientation rather then
		/// face on owner processor.
		/// It may happen that some data depends on normal orientation, then
		/// one should be aware on a local processor orientation may be different
		/// from owner value, then the data may have incorrect sign.
		/// @param mrk Non-private marker that will indicate inverted normals
		void MarkNormalOrientation(MarkerType mrk);
		//implemented in modify.cpp
	private:
		MarkerType hide_element, new_element, update_geometry, temp_hide_element;
	public:
		/// Check whether code runs between Mesh::BeginModification, Mesh::EndModification scope.
		/// In case mesh is modified, on element construction Mesh::TieElements will always place elements 
		/// to the end of the array as a result all the newly created elements will be iterated after current
		/// or hidden elements.
		bool                              isMeshModified     () const {return new_element != 0;} 
		MarkerType                        HideMarker         () const {return hide_element;}
		MarkerType                        NewMarker          () const {return new_element;}
		MarkerType                        UpdateGeometryMarker() const {return update_geometry;}
		void                              SwapModification   (bool recompute_geometry); // swap hidden and new elements, so that old mesh is recovered
		void                              BeginModification  ();  //allow elements to be hidden
		/// After this function any link to deleted element will be replaced by InvalidHandle().
		/// This will modify DATA_REFERENCE tags and contents of sets, so that all deleted elements are not referenced anymore.
		/// If you have any tags of type DATA_REMOTE_REFERENCE on current mesh linking to the elements of the current mesh
		/// or there are other meshes that posses tags of type DATA_REMOTE_REFERENCE and link elements on the current mesh,
		/// you should check that there are no links to deleted elements manually with Element::Old().
		/// \todo
		///      1. maybe instead of forming set of deleted elements and subtracting set from other sets it is better
		///         to remove each modified element
		///         (done, check and compare)
		///      2. parent/child elements in set would not be replaced or reconnected, this may lead to wrong behavior
		///         (done, check and compare)
		/// @see Element::Old
		void                              ApplyModification  ();  //modify DATA_REFERENCE, tags so that links to hidden elements are converted to NULL and removed from sets
		/// This function is not yet implemented. It should correctly resolve parallel state of 
		/// newly created elements, provide them valid global identificators, resolve owners of
		/// the elements potentially optimized using information from BridgeTag and LayersTag
		/// May use ResolveShared function as basis but instead the whole mesh run the same algorithm for subset.
		void                              ResolveModification(); //resolve parallel state of newly created elements, restore ghost layers; not implemented, resuse ResolveShared code
		void                              EndModification    ();    //delete hidden elements
		enumerator                        getNext            (const HandleType * arr, enumerator size, enumerator k, MarkerType marker) const;
		enumerator                        Count              (const HandleType * arr, enumerator size, MarkerType marker) const;
		void                              EquilibrateGhost   ();//bool only_new = false); //Use in ResolveShared
		//void CheckFaces();
		/// Check that centroids of ghost and shared elements match to each other.
		/// Exits if does not match.
		void                              CheckCentroids     (std::string file, int line);
		/// Check that processors are sorted on every element
		void                              CheckProcsSorted   (std::string file, int line);
		/// Check that number of ghost and shared elements match to each other.
		/// Exits if does not match.
		void                              CheckGhostSharedCount(std::string file, int line, ElementType etype = ESET | CELL | FACE | EDGE | NODE);
		/// Let ghost elements send processors list to master elements and see if they match
		void                              CheckOwners        ();
		/// Let ghost elements send owner processor to master elements and see if they match
		void                              CheckProcessors    ();
		/// Checks that there are no invalid links in sets
		void                              CheckSetLinks      (std::string file, int line);
		/// Copy all the data from b to a (a = b).
		/// Except for protected data.
		/// Non-private markers are copied.
		/// Elements should be of same type.
		static void                              CopyData(Element a, Element b);
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
		///
		/// \todo
		/// list checks performed inside in description
		TopologyCheck                     BeginTopologyCheck (ElementType etype, const HandleType * adj, enumerator num);
		/// This function performs some topologycal checks after construction of element.
		/// Function is used internally by CreateEdge, CreateFace, CreateCell functions.
		///
		/// \todo
		/// list checks performed inside in description.
		bool                              EndTopologyCheck   (HandleType e, TopologyCheck begin_check); //check created element
		/// This will return tag by which you can retrieve error mark to any element on which topogy check failed.
		/// As this is sparse tag you can check presence of error by Element::HaveData or Mesh::HaveData check.
		/// This tag will be valid only if you pass MARK_ON_ERROR to Mesh::GetTopologyCheck
		/// and will be deleted if you pass MARK_ON_ERROR to Mesh::RemTopologyCheck
		Tag                               TopologyErrorTag   () const {return tag_topologyerror;}
		/// Retrieve currently set topology checks
		TopologyCheck                     GetTopologyCheck   (TopologyCheck mask = ENUMUNDEF) const {return checkset & mask;}
		/// Set topology checks
		void                              SetTopologyCheck   (TopologyCheck mask);
		/// Remove topology checks
		void                              RemTopologyCheck   (TopologyCheck mask);
		/// This will turn mesh into the state indicating that some topology error occured
		void                              SetTopologyError   (TopologyCheck mask) {errorset = errorset | mask;}
		/// Retrieve topology error state, this indicates that some error have occured
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
			int Compare(const real * a, const real * b) const;
			bool operator() (HandleType a, HandleType b) const;
			bool operator() (HandleType a, const real * b) const;
			bool operator() (const real * a, HandleType b) const;
		};

		class GlobalIDComparator
		{
			Mesh * m;
		public:
			GlobalIDComparator(Mesh * m) :m(m) {}
			GlobalIDComparator(const GlobalIDComparator & other) :m(other.m){}
			GlobalIDComparator & operator = (GlobalIDComparator const & other) { m = other.m; return *this;}
			bool operator() (HandleType a, HandleType b) const {if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->GlobalID(a) < m->GlobalID(b);}
			bool operator() (HandleType a, integer gid) const {if( a == InvalidHandle() ) return false; return m->GlobalID(a) < gid;}
		};

        class SetNameComparator
        {
			Mesh * m;
            public:
			    SetNameComparator(Mesh * m) :m(m) {}
                bool operator()(HandleType a, HandleType b) const
			    {
                    if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; 

                    return ElementSet(m,a).GetName() < ElementSet(m,b).GetName();
                }  

        };

		class HierarchyComparator
		{
			Mesh * m;
		public:
			HierarchyComparator(Mesh * m) :m(m) {}
			HierarchyComparator(const HierarchyComparator & other) :m(other.m){}
			HierarchyComparator & operator = (HierarchyComparator const & other) { m = other.m; return *this;}
			int CompareNodes(HandleType a, HandleType b) const;
			int CompareElements(HandleType a, HandleType b) const;
			bool operator() (HandleType a, HandleType b) const {if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return CompareElements(a,b) < 0;}
		};

		class RealComparator
		{
			Mesh * m; Tag t;
		public:
			RealComparator(Mesh * m, Tag t) :m(m), t(t) {}
			RealComparator(const RealComparator & other) :m(other.m), t(other.t){}
			RealComparator & operator = (RealComparator const & other) { m = other.m; t = other.t; return *this;}
			bool operator() (HandleType a, HandleType b) const {if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->Real(a,t) < m->Real(b,t);}
			bool operator() (HandleType a, real b) const {if( a == InvalidHandle() ) return true; return m->Real(a,t) < b;}
		};

		class IntegerComparator
		{
			Mesh * m; Tag t;
		public:
			IntegerComparator(Mesh * m, Tag t) :m(m), t(t) {}
			IntegerComparator(const IntegerComparator & other) :m(other.m), t(other.t){}
			IntegerComparator & operator = (IntegerComparator const & other) { m = other.m; t = other.t; return *this;}
			bool operator() (HandleType a, HandleType b) const {if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->Integer(a,t) < m->Integer(b,t);}
			bool operator() (HandleType a, integer b) const {if( a == InvalidHandle() ) return true; return m->Integer(a,t) < b;}
		};

		class BulkComparator
		{
			Mesh * m; Tag t;
		public:
			BulkComparator(Mesh * m, Tag t) :m(m), t(t) {}
			BulkComparator(const BulkComparator & other) :m(other.m), t(other.t){}
			BulkComparator & operator = (BulkComparator const & other) { m = other.m; t = other.t; return *this;}
			bool operator() (HandleType a, HandleType b) const {if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->Bulk(a,t) < m->Bulk(b,t);}
			bool operator() (HandleType a, bulk b) const {if( a == InvalidHandle() ) return true; return m->Bulk(a,t) < b;}
		};
		
		class RealDFComparator
		{
			Mesh * m; Tag t;
		public:
			RealDFComparator(Mesh * m, Tag t) :m(m), t(t) {}
			RealDFComparator(const RealDFComparator & other) :m(other.m), t(other.t){}
			RealDFComparator & operator = (RealDFComparator const & other) { m = other.m; t = other.t; return *this;}
			bool operator() (HandleType a, HandleType b) const {if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->RealDF(a,t) < m->RealDF(b,t);}
			bool operator() (HandleType a, real b) const {if( a == InvalidHandle() ) return true; return m->RealDF(a,t) < b;}
		};

		class IntegerDFComparator
		{
			Mesh * m; Tag t;
		public:
			IntegerDFComparator(Mesh * m, Tag t) :m(m), t(t) {}
			IntegerDFComparator(const IntegerDFComparator & other) :m(other.m), t(other.t){}
			IntegerDFComparator & operator = (IntegerDFComparator const & other) { m = other.m; t = other.t; return *this;}
			bool operator() (HandleType a, HandleType b) const {if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->IntegerDF(a,t) < m->IntegerDF(b,t);}
			bool operator() (HandleType a, integer b) const {if( a == InvalidHandle() ) return true; return m->IntegerDF(a,t) < b;}
		};

		class BulkDFComparator
		{
			Mesh * m; Tag t;
		public:
			BulkDFComparator(Mesh * m, Tag t) :m(m), t(t) {}
			BulkDFComparator(const BulkDFComparator & other) :m(other.m), t(other.t){}
			BulkDFComparator & operator = (BulkDFComparator const & other) { m = other.m; t = other.t; return *this;}
			bool operator() (HandleType a, HandleType b) const {if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return m->BulkDF(a,t) < m->BulkDF(b,t);}
			bool operator() (HandleType a, bulk b) const {if( a == InvalidHandle() ) return true; return m->BulkDF(a,t) < b;}
		};
		
		class MeasureComparator
		{
			Mesh * m;
		public:
			MeasureComparator(Mesh * m) :m(m) {}
			MeasureComparator(const MeasureComparator & other) :m(other.m) {}
			MeasureComparator & operator = (MeasureComparator const & other) { m = other.m; return *this;}
			bool operator() (HandleType a, HandleType b) const
			{
				if( a == InvalidHandle() || b == InvalidHandle() )
					return a > b;
				INMOST_DATA_REAL_TYPE ma, mb;
				m->GetGeometricData(a,MEASURE,&ma);
				m->GetGeometricData(b,MEASURE,&mb);
				return ma < mb;
			}
			bool operator() (HandleType a, bulk b) const
			{
				if( a == InvalidHandle() )
					return true;
				INMOST_DATA_REAL_TYPE ma;
				m->GetGeometricData(a,MEASURE,&ma);
				return ma < b;
			}
		};
		
		class MarkerComparator
		{
			Mesh * m; MarkerType mrk; bool inverse;
		public:
			MarkerComparator(Mesh * m, MarkerType mrk, bool inverse = false) :m(m), mrk(mrk), inverse(inverse) {assert(!isPrivate(mrk));}
			MarkerComparator(const MarkerComparator & other) :m(other.m), mrk(other.mrk), inverse(other.inverse){}
			MarkerComparator & operator = (MarkerComparator const & other) { m = other.m; mrk = other.mrk; inverse = other.inverse; return *this;}
			bool operator() (HandleType a, HandleType b) const {if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return ((inverse ^ m->GetMarker(a,mrk))? 1 : 0) < ((inverse ^ m->GetMarker(b,mrk)) ? 1 : 0);}
			bool operator() (HandleType a, bool b) const {if( a == InvalidHandle() ) return true; return ((inverse ^ m->GetMarker(a,mrk))? 1 : 0) < (b ? 1 : 0);}
		};
		
		class PrivateMarkerComparator
		{
			Mesh * m; MarkerType mrk; bool inverse;
		public:
			PrivateMarkerComparator(Mesh * m, MarkerType mrk, bool inverse = false) :m(m), mrk(mrk), inverse(inverse) {assert(isPrivate(mrk));}
			PrivateMarkerComparator(const PrivateMarkerComparator & other) :m(other.m), mrk(other.mrk), inverse(other.inverse){}
			PrivateMarkerComparator & operator = (PrivateMarkerComparator const & other) { m = other.m; mrk = other.mrk; inverse = other.inverse; return *this;}
			bool operator() (HandleType a, HandleType b) const {if( a == InvalidHandle() || b == InvalidHandle() ) return a > b; return ((inverse ^ m->GetPrivateMarker(a,mrk))? 1 : 0) < ((inverse ^ m->GetPrivateMarker(b,mrk))? 1 : 0);}
			bool operator() (HandleType a, bool b) const {if( a == InvalidHandle() ) return true; return ((inverse ^ m->GetPrivateMarker(a,mrk))? 1 : 0) < (b ? 1 : 0);}
		};


		void SortHandles(HandleType * h, enumerator num);
		/// \todo
		/// TODO 53 check that putting global ids to array will be faster
		void SortByGlobalID(HandleType * h, enumerator num);

		/// Retrive the name of the current mesh.
		std::string GetMeshName();
		/// Be careful changing mesh name if you have already established remote links.
		void SetMeshName(std::string new_name);
		/// Find mesh by name.
		static Mesh * GetMesh(std::string name);
	};
	
	
	
	
	///This structure is a helper structure to aid with search of cells by position.
	///Currently the structure is very specific to the step of mesh modification,
	///as it performs search over old elements of the mesh.
	
	class SearchKDTree
	{
	public:
		
		inline static bool cell_point(const Cell & c, const Storage::real p[3]) {return c.Inside(p);}
		template<typename bbox_type>
		inline static int bbox_point(const Storage::real p[3], const bbox_type bbox[6]);
		template<typename bbox_type>
		inline static int bbox_point_print(const Storage::real p[3], const bbox_type bbox[6]);
		template<typename bbox_type>
		inline static void bbox_closest_point(const Storage::real p[3], const bbox_type bbox[6], Storage::real pout[3]);
		template<typename bbox_type>
		inline static int bbox_sphere(const Storage::real p[3], Storage::real r, const bbox_type bbox[6]);
		inline static Storage::real segment_distance(const Storage::real x1[3], const Storage::real x2[3], const Storage::real p[3]);
		inline static Storage::real triangle_distance(const Storage::real x1[3], const Storage::real x2[3], const Storage::real x3[3], const Storage::real p[3]);
	private:
		struct entry
		{
			HandleType e;
			float xyz[3];
			struct entry & operator =(const struct entry & other)
			{
				e = other.e;
				xyz[0] = other.xyz[0];
				xyz[1] = other.xyz[1];
				xyz[2] = other.xyz[2];
				return *this;
			}
		} * set;
		Mesh * m;
		INMOST_DATA_ENUM_TYPE size;
		float bbox[6];
		SearchKDTree * children;
		static int cmpElements0(const void * a,const void * b);
		static int cmpElements1(const void * a,const void * b);
		static int cmpElements2(const void * a,const void * b);
		inline static unsigned int flip(const unsigned int * fp);
		void radix_sort(int dim, struct entry * temp);
		void kdtree_build(int dim, int & done, int total, struct entry * temp);
		SearchKDTree() : set(NULL), m(NULL), size(0), bbox(), children(NULL) {}
		
		Cell SubSearchCell(const Storage::real p[3], bool print) const;
		void clear_children();

		inline int ray_bbox(double pos[3], double ray[3], double closest) const;
		inline int  sphere_bbox(const Storage::real p[3], Storage::real r) const;
		inline int  segment_bbox(const Storage::real p1[3], const Storage::real p2[3]) const;
		inline int  segment_tri(const Storage::real tri[3][3], const Storage::real p1[3], const Storage::real p2[3]) const;
		inline bool segment_face(const Element & f, const Storage::real p1[3], const Storage::real p2[3]) const;
		inline bool segment_cell(const Element & c, const Storage::real p1[3], const Storage::real p2[3]) const;
		inline int  sphere_tri(const Storage::real tri[3][3], const Storage::real p[3], Storage::real r) const;
		inline bool sphere_face(const Element& f, const Storage::real p[3], Storage::real r) const;
		inline bool sphere_cell(const Element& c, const Storage::real p[3], Storage::real r) const;
		void sub_intersect_segment(ElementArray<Element> & hits, MarkerType mrk, const Storage::real p1[3], const Storage::real p2[3]) const;
		void sub_intersect_sphere(ElementArray<Element>& hits, MarkerType mrk, const Storage::real p[3], Storage::real r) const;
	public:
		SearchKDTree(Mesh * m);
		SearchKDTree(Mesh * m, HandleType * _set, unsigned set_size);
		~SearchKDTree();
		Cell SearchCell(const Storage::real * point, bool print = false) const;
		void IntersectSphere(ElementArray<Cell>& cells, const Storage::real p[3], Storage::real r) const;
		void IntersectSphere(ElementArray<Face>& faces, const Storage::real p[3], Storage::real r) const;
		void IntersectSegment(ElementArray<Cell>& cells, const Storage::real p1[3], const Storage::real p2[3]) const;
		void IntersectSegment(ElementArray<Face>& faces, const Storage::real p1[3], const Storage::real p2[3]) const;
	};
	

	//////////////////////////////////////////////////////////////////////
	/// Inline functions for class Storage                              //
	//////////////////////////////////////////////////////////////////////
	__INLINE Storage::real & Storage::Real(const Tag & tag) const
	{
		return GetMeshLink()->Real(GetHandle(),tag);
	}
	__INLINE Storage::integer & Storage::Integer(const Tag & tag)  const
	{
		return GetMeshLink()->Integer(GetHandle(),tag);
	}
	__INLINE Storage::bulk & Storage::Bulk(const Tag & tag)  const
	{
		return GetMeshLink()->Bulk(GetHandle(),tag);
	}
	__INLINE Storage::reference & Storage::Reference(const Tag & tag)  const
	{
		return GetMeshLink()->Reference(GetHandle(),tag);
	}
	__INLINE Storage::remote_reference & Storage::RemoteReference(const Tag & tag)  const
	{
		return GetMeshLink()->RemoteReference(GetHandle(),tag);
	}
	__INLINE Storage::real_array Storage::RealArray(const Tag & tag)  const
	{
		return GetMeshLink()->RealArray(GetHandle(),tag);
	}
	__INLINE Storage::integer_array Storage::IntegerArray(const Tag & tag)  const
	{
		return GetMeshLink()->IntegerArray(GetHandle(),tag);
	}
	__INLINE Storage::bulk_array Storage::BulkArray(const Tag & tag)  const
	{
		return GetMeshLink()->BulkArray(GetHandle(),tag);
	}
	__INLINE Storage::reference_array Storage::ReferenceArray(const Tag & tag)  const
	{
		return GetMeshLink()->ReferenceArray(GetHandle(),tag);
	}
	__INLINE Storage::remote_reference_array Storage::RemoteReferenceArray(const Tag & tag)  const
	{
		return GetMeshLink()->RemoteReferenceArray(GetHandle(),tag);
	}
	__INLINE Storage::real_array Storage::RealArrayDF(const Tag & tag)  const
	{
		return GetMeshLink()->RealArrayDF(GetHandle(),tag);
	}
	__INLINE Storage::integer_array Storage::IntegerArrayDF(const Tag & tag)  const
	{
		return GetMeshLink()->IntegerArrayDF(GetHandle(),tag);
	}
	__INLINE Storage::bulk_array Storage::BulkArrayDF(const Tag & tag)  const
	{
		return GetMeshLink()->BulkArrayDF(GetHandle(),tag);
	}
	__INLINE Storage::reference_array Storage::ReferenceArrayDF(const Tag & tag)  const
	{
		return GetMeshLink()->ReferenceArrayDF(GetHandle(),tag);
	}
	__INLINE Storage::remote_reference_array Storage::RemoteReferenceArrayDF(const Tag & tag)  const
	{
		return GetMeshLink()->RemoteReferenceArrayDF(GetHandle(),tag);
	}
	__INLINE Storage::real & Storage::RealDF(const Tag & tag)  const
	{
		return GetMeshLink()->RealDF(GetHandle(),tag);
	}
	__INLINE Storage::integer & Storage::IntegerDF(const Tag & tag)  const
	{
		return GetMeshLink()->IntegerDF(GetHandle(),tag);
	}
	__INLINE Storage::bulk & Storage::BulkDF(const Tag & tag)  const
	{
		return GetMeshLink()->BulkDF(GetHandle(),tag);
	}
	__INLINE Storage::reference & Storage::ReferenceDF(const Tag & tag)  const
	{
		return GetMeshLink()->ReferenceDF(GetHandle(),tag);
	}
	__INLINE Storage::remote_reference & Storage::RemoteReferenceDF(const Tag & tag)  const
	{
		return GetMeshLink()->RemoteReferenceDF(GetHandle(),tag);
	}
	__INLINE Storage::real_array Storage::RealArrayDV(const Tag & tag)  const
	{
		return GetMeshLink()->RealArrayDV(GetHandle(),tag);
	}
	__INLINE Storage::integer_array Storage::IntegerArrayDV(const Tag & tag)  const
	{
		return GetMeshLink()->IntegerArrayDV(GetHandle(),tag);
	}
	__INLINE Storage::bulk_array Storage::BulkArrayDV(const Tag & tag)  const
	{
		return GetMeshLink()->BulkArrayDV(GetHandle(),tag);
	}
	__INLINE Storage::reference_array Storage::ReferenceArrayDV(const Tag & tag)  const
	{
		return GetMeshLink()->ReferenceArrayDV(GetHandle(),tag);
	}
	__INLINE Storage::remote_reference_array Storage::RemoteReferenceArrayDV(const Tag & tag)  const
	{
		return GetMeshLink()->RemoteReferenceArrayDV(GetHandle(),tag);
	}
	__INLINE Storage::real & Storage::RealDV(const Tag & tag)  const
	{
		return GetMeshLink()->RealDV(GetHandle(),tag);
	}
	__INLINE Storage::integer & Storage::IntegerDV(const Tag & tag)  const
	{
		return GetMeshLink()->IntegerDV(GetHandle(),tag);
	}
	__INLINE Storage::bulk & Storage::BulkDV(const Tag & tag)  const
	{
		return GetMeshLink()->BulkDV(GetHandle(),tag);
	}
	__INLINE Storage::reference & Storage::ReferenceDV(const Tag & tag)  const
	{
		return GetMeshLink()->ReferenceDV(GetHandle(),tag);
	}
	__INLINE Storage::remote_reference & Storage::RemoteReferenceDV(const Tag & tag)  const
	{
		return GetMeshLink()->RemoteReferenceDV(GetHandle(),tag);
	}
	__INLINE Storage::real & TagReal::operator [](HandleType h) const
	{
		return GetMeshLink()->Real(h,*static_cast<const Tag*>(this));
	}
	__INLINE Storage::integer & TagInteger::operator [](HandleType h) const
	{
		return GetMeshLink()->Integer(h,*static_cast<const Tag*>(this));
	}
	__INLINE Storage::bulk & TagBulk::operator [](HandleType h) const
	{
		return GetMeshLink()->Bulk(h,*static_cast<const Tag*>(this));
	}
	__INLINE Storage::reference & TagReference::operator [](HandleType h) const
	{
		return GetMeshLink()->Reference(h,*static_cast<const Tag*>(this));
	}
	__INLINE Storage::real_array TagRealArray::operator [](HandleType h) const
	{
		return GetMeshLink()->RealArray(h,*static_cast<const Tag*>(this));
	}
	__INLINE Matrix<Storage::real,Storage::real_array> TagRealArray::operator ()(HandleType h, int n, int m) const
	{
		Storage::real_array data = GetMeshLink()->RealArray(h,*static_cast<const Tag*>(this));
		assert((int)data.size() == n*m);
		return Matrix<Storage::real,Storage::real_array>(data,n,m);
	}
	__INLINE Storage::integer_array TagIntegerArray::operator [](HandleType h) const
	{
		return GetMeshLink()->IntegerArray(h,*static_cast<const Tag*>(this));
	}
	__INLINE Storage::bulk_array TagBulkArray::operator [](HandleType h) const
	{
		return GetMeshLink()->BulkArray(h,*static_cast<const Tag*>(this));
	}
	__INLINE Storage::reference_array TagReferenceArray::operator [](HandleType h) const
	{
		return GetMeshLink()->ReferenceArray(h,*static_cast<const Tag*>(this));
	}
#if defined(USE_AUTODIFF)
	__INLINE Storage::var & Storage::Variable(const Tag & tag) const
	{
		return GetMeshLink()->Variable(GetHandle(),tag);
	}
	__INLINE Storage::var & Storage::VariableDF(const Tag & tag) const
	{
		return GetMeshLink()->VariableDF(GetHandle(),tag);
	}
	__INLINE Storage::var & Storage::VariableDV(const Tag & tag) const
	{
		return GetMeshLink()->VariableDV(GetHandle(),tag);
	}
	__INLINE Storage::var_array Storage::VariableArray(const Tag & tag)  const
	{
		return GetMeshLink()->VariableArray(GetHandle(),tag);
	}
	__INLINE Storage::var_array Storage::VariableArrayDF(const Tag & tag)  const
	{
		return GetMeshLink()->VariableArrayDF(GetHandle(),tag);
	}
	__INLINE Storage::var_array Storage::VariableArrayDV(const Tag & tag)  const
	{
		return GetMeshLink()->VariableArrayDV(GetHandle(),tag);
	}
	__INLINE Storage::var & TagVariable::operator [](HandleType h) const
	{
		return GetMeshLink()->Variable(h,*static_cast<const Tag*>(this));
	}
	__INLINE Storage::var_array TagVariableArray::operator [](HandleType h) const
	{
		return GetMeshLink()->VariableArray(h,*static_cast<const Tag*>(this));
	}
	__INLINE Matrix<Storage::var,Storage::var_array> TagVariableArray::operator ()(HandleType h, int n, int m) const
	{
		Storage::var_array data = GetMeshLink()->VariableArray(h,*static_cast<const Tag*>(this));
		assert((int)data.size() == n*m);
		return Matrix<Storage::var,Storage::var_array>(data,n,m);
	}
	
#endif
	__INLINE bool Storage::HaveData(const Tag & tag) const
	{
		assert(isValid());
		return GetMeshLink()->HaveData(GetHandle(),tag);
	}
	__INLINE INMOST_DATA_ENUM_TYPE Storage::GetDataSize(const Tag & tag) const
	{
		return GetMeshLink()->GetDataSize(GetHandle(),tag);
	}
	__INLINE INMOST_DATA_ENUM_TYPE Storage::GetDataCapacity(const Tag & tag) const
	{
		return GetMeshLink()->GetDataCapacity(GetHandle(),tag);
	}
	__INLINE void Storage::SetDataSize(const Tag & tag,INMOST_DATA_ENUM_TYPE new_size) const
	{
		GetMeshLink()->SetDataSize(GetHandle(),tag,new_size);
	}
	__INLINE ElementType Storage::GetElementType() const
	{
		return GetHandleElementType(handle);
	}
	__INLINE Storage::integer Storage::GetElementNum   () const
	{
		return GetHandleElementNum(handle);
	}
	__INLINE void Storage::GetData(const Tag & tag,INMOST_DATA_ENUM_TYPE shift, INMOST_DATA_ENUM_TYPE size, void * data_out) const
	{
		GetMeshLink()->GetData(GetHandle(),tag,shift,size,data_out);
	}
	__INLINE void Storage::SetData(const Tag & tag,INMOST_DATA_ENUM_TYPE shift, INMOST_DATA_ENUM_TYPE size, const void * data_in) const
	{
		GetMeshLink()->SetData(GetHandle(),tag,shift,size,data_in);
	}
	__INLINE void Storage::DelData(const Tag & tag) const
	{
		GetMeshLink()->DelData(GetHandle(),tag);
	}
	__INLINE void Storage::DelDenseData(const Tag & tag) const
	{
		GetMeshLink()->DelDenseData(GetHandle(),tag);
	}
	__INLINE bool Storage::DelSparseData(const Tag & tag) const
	{
		return GetMeshLink()->DelSparseData(GetHandle(),tag);
	}
	__INLINE void Storage::SetMarker(MarkerType n)  const
	{
		assert( isValid() );
		GetMeshLink()->SetMarker(GetHandle(),n);
	}
	__INLINE bool Storage::GetMarker(MarkerType n) const
	{
		assert( isValid() );
		return GetMeshLink()->GetMarker(GetHandle(),n);
	}
	__INLINE void Storage::RemMarker(MarkerType n)  const
	{
		assert( isValid() );
		GetMeshLink()->RemMarker(GetHandle(),n);
	}
	__INLINE void Storage::SetPrivateMarker(MarkerType n)  const
	{
		assert( isValid() );
		GetMeshLink()->SetPrivateMarker(GetHandle(),n);
	}
	__INLINE bool Storage::GetPrivateMarker(MarkerType n) const
	{
		assert( isValid() );
		return GetMeshLink()->GetPrivateMarker(GetHandle(),n);
	}
	__INLINE void Storage::RemPrivateMarker(MarkerType n)  const
	{
		assert( isValid() );
		GetMeshLink()->RemPrivateMarker(GetHandle(),n);
	}
	__INLINE void Storage::ClearMarkerSpace()  const
	{
		GetMeshLink()->ClearMarkerSpace(GetHandle());
	}
	__INLINE void Storage::GetMarkerSpace(Storage::bulk copy[MarkerFields]) const
	{
		GetMeshLink()->GetMarkerSpace(GetHandle(),copy);
	}
	__INLINE void Storage::SetMarkerSpace(Storage::bulk source[MarkerFields])  const
	{
		GetMeshLink()->SetMarkerSpace(GetHandle(),source);
	}
	__INLINE bool Storage::isValid() const
	{
		return handle != InvalidHandle() && GetMeshLink() != NULL && GetMeshLink()->isValidElement(handle);
	}
	__INLINE Storage::integer Storage::LocalID() const
	{
		return GetHandleID(handle);
	}
	__INLINE Storage::integer Storage::DataLocalID() const
	{
		return GetMeshLink()->DataLocalID(GetHandle());
	}
	__INLINE Element Storage::getAsElement() const
	{
		assert(GetElementType() & (NODE | EDGE | FACE | CELL | ESET) );
		return Element(GetMeshLink(), GetHandle());
	}
	__INLINE Node Storage::getAsNode() const
	{
		assert(GetElementType() == NODE);
		return Node(GetMeshLink(),GetHandle());
	}
	__INLINE Edge Storage::getAsEdge() const
	{
		assert(GetElementType() == EDGE);
		return Edge(GetMeshLink(),GetHandle());
	}
	__INLINE Face Storage::getAsFace() const
	{
		assert(GetElementType() == FACE);
		return Face(GetMeshLink(),GetHandle());
	}
	__INLINE Cell Storage::getAsCell() const
	{
		assert(GetElementType() == CELL);
		return Cell(GetMeshLink(),GetHandle());
	}
	__INLINE ElementSet Storage::getAsSet() const
	{
		assert(GetElementType() == ESET);
		return ElementSet(GetMeshLink(),GetHandle());
	}
	__INLINE Mesh * Storage::GetMeshLink() const
	{
		return m_link;
	}
	__INLINE HandleType Storage::GetHandle() const
	{
		return handle;
	}
	
	

}

  //Implementation of inlined functions
//#include "../Data/storage_inline.hpp"

#endif

#endif // INMOST_MESH_H_INCLUDED
