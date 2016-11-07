#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "inmost.h"
#include "../Mesh/incident_matrix.hpp"
#include <string>
#if defined(USE_MESH)

// coords/zcorn algorithm
// 1) put all block nodes into pillars, sort each pillar nodes by actual depth
// 2) create edges down along pillars and mark them according to blocks
// 3) consider all pairs of pillars in Oxz and Oyz planes and create uncut block edges, mark them in block numbers
// 4) use line-sweep algorithm to intersect and cut all the edges between each pair of pillars, mark them in block numbers
// 5) mark all the nodes with block numbers that they belong by considering union of block numbers on adjacent edges
// 5) use incident_matrix algorithm to create all the faces between pair of pillars from pillar edges and inter-pillar edges
// 6) from intersection of block numbers on nodes figure out which blocks the face belongs
// 7) add top and bottom interface

//todo
// 1. (ok) remove edges on pillars that do not get any associated block
// 2. (ok) do not create cells with actnum = 0
// 2. fix fast intersect() algorithm, check against intersect_naive()
// 3. CheckEdgeOrder reports bad order of edges for certain interfaces - investigate
// 4. when encounter cells with one face with actnum = 1 should disconnect adjacent cells
// 5. (ok) populate keywords data, poro, perm, actnum, satnum, etc...
// 5.1 add saturation, pressure
// 6. local grid refinement
// 6.1 radial local grid refinement
// 6.2 nested local grid refinement
// 7. (ok) omp-parallel loading
// 8. mpi-parallel loading
// 9. (ok) do not convert arrays xyz and zcorn
//10. read wells into sets
//10.1 (ok) read compdat
//10.2 read properties
//11. optimize:
//11.1 detect that pillars have no faults, skip intersection in this case
//11.2 skip incident_matrix algorithm when no inner edges found

//eclipse states
#define ECLSTRCMP(x,y) strncmp(x,y,8)


//controls when to consider two nodes on pillar to be the same
#define ECL_PILLAR_DEPTH_EPS 1.0e-8
#define ECL_POINT_EPS 1.0e-8
#define ECL_INTERSECT_EPS 1.0e-12

#define ECL_NEVER -1
#define ECL_NONE 0
#define ECL_SKIP_SECTION 1
#define ECL_INCLUDE 2
#define ECL_DIMENS 3
#define ECL_DX 4
#define ECL_DY 5
#define ECL_DZ 6
#define ECL_TOPS 7
#define ECL_PERMX 8
#define ECL_PERMY 9
#define ECL_PERMZ 10
#define ECL_PORO 11
#define ECL_MAPAXIS 12
#define ECL_INRAD 13
#define ECL_COORDS 14
#define ECL_ZCORN 15
#define ECL_ACTNUM 16
#define ECL_SATNUM 17
#define ECL_COMPDAT 18

#define ECL_GTYPE_NONE 0
#define ECL_GTYPE_TOPS 1
#define ECL_GTYPE_ZCORN 2
#define ECL_GTYPE_RADIAL 3
#define ECL_GTYPE_CARTESIAN 4


#define ECL_VAR_NONE 0
#define ECL_VAR_REAL 1
#define ECL_VAR_INT 2

#define ECL_IJK_DATA(i,j,k) (i + ((j)+(k)*dims[1])*dims[0])
#define ECL_IJK_COORDS(i,j,l,coord) ((((i) + (j)*(dims[0]+1))*2+(l))*3+(coord))
#define ECL_IJK_ZCORN(i,j,k,l) ((k)*dims[1]*dims[0]*8 + (1-(l)/4)*dims[1]*dims[0]*4 + (j)*dims[0]*4 + ((l)/2%2)*dims[0]*2 + (i)*2 + (l)%2)

//line sweep events
#define SEG_START 0
#define SEG_END 1



namespace INMOST
{
	static std::string GetFolder(std::string file)
	{
		size_t found = file.find_last_of("/\\");
		if( found == std::string::npos )
			return "";
		else return file.substr(0,found);
	}
	struct compdat_entry
	{
		int i, j, k1, k2;
		bool open;
	};
	typedef std::vector<compdat_entry> compdat_entries;
	typedef std::map<std::string, compdat_entries > compdat_wells;
	//2d point with comparison operator for line sweep algorithm
	struct Point
	{
		HandleType node;
		Storage::real x, y;
		Point & operator = (Point const & b) { node = b.node; x = b.x; y = b.y; return *this; }
		Point(const Point & b) : node(b.node), x(b.x), y(b.y) {}
		Point(HandleType _node, Storage::real _x, Storage::real _y) : node(_node), x(_x), y(_y) {}
		Point(HandleType _node, Storage::real v[2]) : node(_node), x(v[0]), y(v[1]) {}
		bool operator <(const Point & b) const
		{
			if (y < b.y - ECL_POINT_EPS) return true;
			else if (y > b.y + ECL_POINT_EPS) return false;
			else if (x < b.x - ECL_POINT_EPS) return true;
			else return false;
		}
		bool operator ==(const Point & b) const
		{
			return std::abs(y - b.y) < ECL_POINT_EPS && std::abs(x - b.x) < ECL_POINT_EPS;
		}
		bool operator !=(const Point & b) const
		{
			return std::abs(y - b.y) > ECL_POINT_EPS || std::abs(x - b.x) > ECL_POINT_EPS;
		}
	};
	//class that returns a Point into 3d space based on coords of two pillars
	class Unproject
	{
		Storage::real a0[3], a1[3], b0[3], b1[3];
	public:
		Unproject(Storage::real _a0[3], Storage::real _a1[3], Storage::real _b0[3], Storage::real _b1[3])
		{
			a0[0] = _a0[0];
			a0[1] = _a0[1];
			a0[2] = _a0[2];
			a1[0] = _a1[0];
			a1[1] = _a1[1];
			a1[2] = _a1[2];
			b0[0] = _b0[0];
			b0[1] = _b0[1];
			b0[2] = _b0[2];
			b1[0] = _b1[0];
			b1[1] = _b1[1];
			b1[2] = _b1[2];
		}
		Unproject(const Unproject & b)
		{
			a0[0] = b.a0[0];
			a0[1] = b.a0[1];
			a0[2] = b.a0[2];
			a1[0] = b.a1[0];
			a1[1] = b.a1[1];
			a1[2] = b.a1[2];
			b0[0] = b.b0[0];
			b0[1] = b.b0[1];
			b0[2] = b.b0[2];
			b1[0] = b.b1[0];
			b1[1] = b.b1[1];
			b1[2] = b.b1[2];
		}
		Unproject & operator =(Unproject const & b)
		{
			a0[0] = b.a0[0];
			a0[1] = b.a0[1];
			a0[2] = b.a0[2];
			a1[0] = b.a1[0];
			a1[1] = b.a1[1];
			a1[2] = b.a1[2];
			b0[0] = b.b0[0];
			b0[1] = b.b0[1];
			b0[2] = b.b0[2];
			b1[0] = b.b1[0];
			b1[1] = b.b1[1];
			b1[2] = b.b1[2];
			return *this;
		}
		
		void ActTri(const Point & p, Storage::real ret[3]) const
		{
			//Storage::real ret1[3],ret2[3];
			Storage::real px = p.y, py = 1-p.x;
			
			// a1         b1
			//          /
			//         /
			//        /
			//       /
			//      /
			//     /
			//    /
			// a0         b0
			if( py > px ) //tri a0,a1,b1
			{
				//a = (0,0) (a0) u
				//b = (0,1) (a1) v
				//c = (1,1) (b1) w
				//Vector v0 = b-a = (0,1), v1 = c-a =(1,1), v2 = p-a;
				const Storage::real d00 = 1; //Dot(v0, v0);
				const Storage::real d01 = 1; //Dot(v0, v1);
				const Storage::real d11 = 2; //Dot(v1, v1);
				const Storage::real denom = d00 * d11 - d01 * d01;
				Storage::real d20 = py; //Dot(v2, v0);
				Storage::real d21 = px+py; //Dot(v2, v1);
				Storage::real v = (d11 * d20 - d01 * d21) / denom; //for a1
				Storage::real w = (d00 * d21 - d01 * d20) / denom; //for b1
				Storage::real u = 1.0 - v - w; //for a0
				for(int k = 0; k < 3; ++k)
					ret[k] = a0[k]*u + a1[k]*v + b1[k]*w;
			}
			else //tri a0,b0,b1
			{
				//a = (0,0) (a0) u
				//b = (1,0) (b0) v
				//c = (1,1) (b1) w
				//Vector v0 = b-a = (1,0), v1 = c-a = (1,1), v2 = p-a;
				const Storage::real d00 = 1; //Dot(v0, v0);
				const Storage::real d01 = 1; //Dot(v0, v1);
				const Storage::real d11 = 2; //Dot(v1, v1);
				const Storage::real denom = d00 * d11 - d01 * d01;
				Storage::real d20 = px; //Dot(v2, v0);
				Storage::real d21 = px+py; //Dot(v2, v1);
				Storage::real v = (d11 * d20 - d01 * d21) / denom; //for b0
				Storage::real w = (d00 * d21 - d01 * d20) / denom; //for b1
				Storage::real u = 1.0 - v - w; //for a0
				for(int k = 0; k < 3; ++k)
					ret[k] = a0[k]*u + b0[k]*v + b1[k]*w;
			}
			
			/*
			// a1         b1
			//    \
			//     \
			//      \
			//       \
			//        \
			//         \
			//          \
			// a0         b0
			if( py < 1-px ) //tri a0,a1,b0
			{
				//a = (0,0) (a0) u
				//b = (0,1) (a1) v
				//c = (1,0) (b0) w
				//Vector v0 = b-a = (0,1), v1 = c-a = (1,0), v2 = p-a;
				const Storage::real d00 = 1; //Dot(v0, v0);
				const Storage::real d01 = 0; //Dot(v0, v1);
				const Storage::real d11 = 1; //Dot(v1, v1);
				const Storage::real denom = d00 * d11 - d01 * d01;
				Storage::real d20 = py; //Dot(v2, v0);
				Storage::real d21 = px; //Dot(v2, v1);
				Storage::real v = (d11 * d20 - d01 * d21) / denom; //for a1
				Storage::real w = (d00 * d21 - d01 * d20) / denom; //for b0
				Storage::real u = 1.0 - v - w; //for a0
				for(int k = 0; k < 3; ++k)
					ret[k] = a0[k]*u + a1[k]*v + b0[k]*w;
			}
			else //tri a0,b0,b1
			{
				//a = (0,1) (a1) u
				//b = (1,0) (b0) v
				//c = (1,1) (b1) w
				//Vector v0 = b-a = (1,-1), v1 = c-a = (1,0), v2 = p-a = (px,py-1);
				const Storage::real d00 = 2; //Dot(v0, v0);
				const Storage::real d01 = 1; //Dot(v0, v1);
				const Storage::real d11 = 1; //Dot(v1, v1);
				const Storage::real denom = d00 * d11 - d01 * d01;
				Storage::real d20 = px-py+1; //Dot(v2, v0);
				Storage::real d21 = px; //Dot(v2, v1);
				Storage::real v = (d11 * d20 - d01 * d21) / denom; //for b0
				Storage::real w = (d00 * d21 - d01 * d20) / denom; //for b1
				Storage::real u = 1.0 - v - w; //for a1
				for(int k = 0; k < 3; ++k)
					ret[k] = a1[k]*u + b0[k]*v + b1[k]*w;
			}
			*/
			//for(int k = 0; k < 3; ++k) ret[k] = (ret1[k]+ret2[k])*0.5;
			
		}
		
		 
		
		void Act(const Point & p, Storage::real v[3]) const
		{
			Storage::real a, b;
			for(int k = 0; k < 3; ++k)
			{
				a = (p.x)*a0[k] + (1.0-p.x)*a1[k];
				b = (p.x)*b0[k] + (1.0-p.x)*b1[k];
				v[k] = (1.0-p.y)*a + p.y*b;
			}
		}
		
	};
	// Used to sort arrays by indices
	class index_comparator
	{
		Storage::integer_array & data;
	public:
		index_comparator(Storage::integer_array & data) : data(data) {}
		index_comparator(const index_comparator & b) : data(b.data) {}
		index_comparator & operator =(index_comparator const & b) {data = b.data; return *this;}
		bool operator ()(int a, int b) const {return data[a] < data[b];}
	};
	// Used to sort arrays by depth
	class depth_comparator
	{
		const Storage::real * data;
	public:
		depth_comparator(const Storage::real * data) : data(data) {}
		depth_comparator(const depth_comparator & b) : data(b.data) {}
		depth_comparator & operator =(depth_comparator const & b) {data = b.data; return *this;}
		bool operator ()(int a, int b) const {return data[a] < data[b];}
	};
	
	template<typename T>
	int count_duplicates(ElementArray<T> & array)
	{
		Mesh * m = array.GetMeshLink();
		MarkerType mrk = m->CreatePrivateMarker();
		int dups = 0;
		for(int k = 0; k < array.size(); ++k)
		{
			if( array[k].GetPrivateMarker(mrk) )
				dups++;
			array[k].SetPrivateMarker(mrk);
		}
		array.RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return dups;
	}
	
	template<typename T>
	void make_unique(ElementArray<T> & array)
	{
		Mesh * m = array.GetMeshLink();
		MarkerType mrk = m->CreatePrivateMarker();
		int i = 0, j = 0;
		while(j < array.size())
		{
			if( array[j].GetPrivateMarker(mrk) )
				++j;
			else
			{
				array[j].SetPrivateMarker(mrk);
				array[i++] = array[j++];
			}
		}
		array.RemPrivateMarker(mrk);
		array.resize(i);
		m->ReleasePrivateMarker(mrk);
	}
	
	Point make_point(Node n, Tag pnt)
	{
		return Point(n->GetHandle(),n->RealArray(pnt).data());
	}

	std::pair<Storage::real,Storage::real> intersect_pairs(const Point & pabeg, const Point & paend, const Point & pbbeg, const Point & pbend)
	{
		const Storage::real eps = ECL_INTERSECT_EPS;
		Point pfind(InvalidHandle(),0,0);
		Storage::real t1,t2;
		Storage::real div = (pabeg.x - paend.x)*(pbbeg.y - pbend.y) - (pabeg.y - paend.y)*(pbbeg.x - pbend.x);
		if (std::abs(div) < 1.0e-50)
			return std::make_pair(-1,-1);
		pfind.x = ((pabeg.x*paend.y - pabeg.y*paend.x)*(pbbeg.x - pbend.x) - (pabeg.x - paend.x)*(pbbeg.x*pbend.y - pbbeg.y*pbend.x)) / div;
		pfind.y = ((pabeg.x*paend.y - pabeg.y*paend.x)*(pbbeg.y - pbend.y) - (pabeg.y - paend.y)*(pbbeg.x*pbend.y - pbbeg.y*pbend.x)) / div;
		//optimization! - should belong unit cube
		if( pfind.x < 0 || pfind.x > 1 || pfind.y < 0 || pfind.y > 1 ) return std::make_pair(-1,-1);
		if (std::abs(paend.x - pabeg.x) > eps)
		{
			t1 = (pfind.x - pabeg.x) / (paend.x - pabeg.x);
			if (t1 < eps || t1 > 1.0 - eps)  { return std::make_pair(-1,-1); }
		}
		if (std::abs(paend.y - pabeg.y) > eps)
		{
			t1 = (pfind.y - pabeg.y) / (paend.y - pabeg.y);
			if (t1 < eps || t1 > 1.0 - eps)  { return std::make_pair(-1,-1); }
		}
		if (std::abs(pbend.x - pbbeg.x) > eps)
		{
			t2 = (pfind.x - pbbeg.x) / (pbend.x - pbbeg.x);
			if (t2 < eps || t2 > 1.0 - eps)  { return std::make_pair(-1,-1); }
		}
		if (std::abs(pbend.y - pbbeg.y) > eps)
		{
			t2 = (pfind.y - pbbeg.y) / (pbend.y - pbbeg.y);
			if (t2 < eps || t2 > 1.0 - eps)  { return std::make_pair(-1,-1); }
		}
		return std::make_pair(t1,t2);
	}
	
	std::pair<bool,Point> intersect_segments(Mesh * m, const Edge & a, const Edge & b, std::set<Point> & intersections, Tag pnt, const Unproject & unp, bool print)
	{
		const Storage::real eps = ECL_INTERSECT_EPS;
		Point pfind(InvalidHandle(),0,0);
		if( a->getBeg() == b->getBeg() ||
		    a->getBeg() == b->getEnd() ||
		    a->getEnd() == b->getEnd() ||
		    a->getEnd() == b->getBeg() )
			return std::make_pair(false,pfind);
		Point pabeg = make_point(a->getBeg(),pnt);
		Point paend = make_point(a->getEnd(),pnt);
		Point pbbeg = make_point(b->getBeg(),pnt);
		Point pbend = make_point(b->getEnd(),pnt);
		Storage::real div = (pabeg.x - paend.x)*(pbbeg.y - pbend.y) - (pabeg.y - paend.y)*(pbbeg.x - pbend.x);
		Storage::real t1,t2;
		if (std::abs(div) < 1.0e-50)
		{
			if (print) std::cout << "divisor is zero" << std::endl;
			return std::make_pair(false,pfind);
		}
		pfind.x = ((pabeg.x*paend.y - pabeg.y*paend.x)*(pbbeg.x - pbend.x) - (pabeg.x - paend.x)*(pbbeg.x*pbend.y - pbbeg.y*pbend.x)) / div;
		pfind.y = ((pabeg.x*paend.y - pabeg.y*paend.x)*(pbbeg.y - pbend.y) - (pabeg.y - paend.y)*(pbbeg.x*pbend.y - pbbeg.y*pbend.x)) / div;
		//optimization uses information that we stay in unit cube
		if( pfind.x < 0 || pfind.x > 1 || pfind.y < 0 || pfind.y > 1 )
			return std::make_pair(false,pfind);
		if (print) std::cout << "found ("<< pfind.x << ", " << pfind.y << ") for edges " << a->GetHandle() << " and " << b->GetHandle() << std::endl;
		//probably some of these tests are redundant
		if (std::abs(paend.x - pabeg.x) > eps)
		{
			t1 = (pfind.x - pabeg.x) / (paend.x - pabeg.x);
			if (t1 < eps || t1 > 1.0 - eps)  { if (print) std::cout << "out of bound: " << t1 << std::endl; return std::make_pair(false,pfind); }
		}
		if (std::abs(paend.y - pabeg.y) > eps)
		{
			t1 = (pfind.y - pabeg.y) / (paend.y - pabeg.y);
			if (t1 < eps || t1 > 1.0 - eps)  { if (print) std::cout << "out of bound: " << t1 << std::endl; return std::make_pair(false,pfind); }
		}
		if (std::abs(pbend.x - pbbeg.x) > eps)
		{
			t2 = (pfind.x - pbbeg.x) / (pbend.x - pbbeg.x);
			if (t2 < eps || t2 > 1.0 - eps)  { if (print) std::cout << "out of bound: " << t2 << std::endl; return std::make_pair(false,pfind); }
		}
		if (std::abs(pbend.y - pbbeg.y) > eps)
		{
			t2 = (pfind.y - pbbeg.y) / (pbend.y - pbbeg.y);
			if (t2 < eps || t2 > 1.0 - eps)  { if (print) std::cout << "out of bound: " << t2 << std::endl; return std::make_pair(false,pfind); }
		}
		std::set<Point>::iterator search = intersections.find(pfind);
		//check whether intersection already exists
		if( search != intersections.end() ) //there is a node!
		{
			pfind = *search;
			if( pfind.node == a->getBeg()->GetHandle() ||
			    pfind.node == b->getBeg()->GetHandle() ||
			    pfind.node == a->getEnd()->GetHandle() ||
			    pfind.node == b->getEnd()->GetHandle())
				return std::make_pair(false,pfind);
		}
		else //no node, create one
		{
			Storage::real find[3];
			//restore coordinate
			//this is treaky, we do not want self intersection in 3d space
			unp.Act(pfind,find);
			//find[0] = (a->getEnd()->Coords()[0]*t1 + a->getBeg()->Coords()[0]*(1-t1)+b->getEnd()->Coords()[0]*t2 + b->getBeg()->Coords()[0]*(1-t2))*0.5;
			//find[1] = (a->getEnd()->Coords()[1]*t1 + a->getBeg()->Coords()[1]*(1-t1)+b->getEnd()->Coords()[1]*t2 + b->getBeg()->Coords()[1]*(1-t2))*0.5;
			//find[2] = (a->getEnd()->Coords()[2]*t1 + a->getBeg()->Coords()[2]*(1-t1)+b->getEnd()->Coords()[2]*t2 + b->getBeg()->Coords()[2]*(1-t2))*0.5;
			pfind.node = m->CreateNode(find)->GetHandle();
			if (print) std::cout << "intersection accepted (" << find[0] << "," << find[1] << "," << find[2] << ") t1 " << t1 << " t2 " << t2 << " new node " << pfind.node << std::endl;
			Storage::real_array _pfind = m->RealArray(pfind.node,pnt);
			_pfind[0] = pfind.x;
			_pfind[1] = pfind.y;
			intersections.insert(pfind);
		}
		return std::make_pair(true,pfind);
	}

	void split_edge(Mesh * m, Node I, Edge a, ElementArray<Edge> & splitted_a, std::vector<Tag> & transfer, bool print)
	{
		std::vector< std::vector<char> > copy(transfer.size());
		//memorize data
		for(int k = 0; k < transfer.size(); ++k) //tags
		{
			int size = a->GetDataSize(transfer[k]);
				copy[k].resize(transfer[k].GetBytesSize()*size);
			if( !copy.empty() ) a->GetData(transfer[k],0,size,&copy[k][0]);
		}
		
		if( print ) std::cout << "split a " << a->GetHandle() << " " << a->getBeg()->GetHandle() << " <-> " << a->getEnd()->GetHandle() << ":-> ";
		splitted_a = Edge::SplitEdge(a,ElementArray<Node>(m,1,I->GetHandle()),0);
		if( print ) std::cout << splitted_a[0]->GetHandle() << " " << splitted_a[0]->getBeg()->GetHandle() << " <-> " << splitted_a[0]->getEnd()->GetHandle() << " and " << splitted_a[1]->GetHandle() << " " << splitted_a[1]->getBeg()->GetHandle() << " <-> " << splitted_a[1]->getEnd()->GetHandle() << std::endl;
		//duplicate data
		for(int k = 0; k < transfer.size();++k)
		{
			int size = (int)copy[k].size() / transfer[k].GetBytesSize();
			if( size ) for(int l = 0; l < 2; ++l) //two parts
			{
				splitted_a[l].SetDataSize(transfer[k],size);
				splitted_a[l].SetData(transfer[k],0,size,&copy[k][0]);
			}
		}
	}

	void split_edges(Mesh * m, Node I, Edge a, Edge b, ElementArray<Edge> & splitted_a, ElementArray<Edge> & splitted_b,std::vector<Tag> & transfer, bool print)
	{
		split_edge(m,I,a,splitted_a,transfer,print);
		split_edge(m,I,b,splitted_b,transfer,print);
		/*
		//storage for data
		std::vector< std::vector<char> > copy(transfer.size()*2);
		//memorize data
		{
			const Edge * s[2] = {&a,&b};
			for(int q = 0; q < 2; ++q) //segments
				for(int k = 0; k < transfer.size(); ++k) //tags
				{
					int size = s[q]->GetDataSize(transfer[k]);
					copy[k + q*transfer.size()].resize(transfer[k].GetBytesSize()*size);
					if( !copy.empty() ) s[q]->GetData(transfer[k],0,size,&copy[k + q*transfer.size()][0]);
				}
		}
		if( print ) std::cout << "split a " << a->GetHandle() << " " << a->getBeg()->GetHandle() << " <-> " << a->getEnd()->GetHandle() << ":-> ";
		splitted_a = Edge::SplitEdge(a,ElementArray<Node>(m,1,I->GetHandle()),0);
		if( print ) std::cout << splitted_a[0]->GetHandle() << " " << splitted_a[0]->getBeg()->GetHandle() << " <-> " << splitted_a[0]->getEnd()->GetHandle() << " and " << splitted_a[1]->GetHandle() << " " << splitted_a[1]->getBeg()->GetHandle() << " <-> " << splitted_a[1]->getEnd()->GetHandle() << std::endl;
		if( print ) std::cout << "split b " << b->GetHandle() << " " << b->getBeg()->GetHandle() << " <-> " << b->getEnd()->GetHandle() << ":-> ";
		splitted_b = Edge::SplitEdge(b,ElementArray<Node>(m,1,I->GetHandle()),0);
		if( print ) std::cout << splitted_b[0]->GetHandle() << " " << splitted_b[0]->getBeg()->GetHandle() << " <-> " << splitted_b[0]->getEnd()->GetHandle() << " and " << splitted_b[1]->GetHandle() << " " << splitted_b[1]->getBeg()->GetHandle() << " <-> " << splitted_b[1]->getEnd()->GetHandle() << std::endl;
		//duplicate data
		{
			const Edge splitted[2][2] =
			{
				{splitted_a[0],splitted_a[1]},
				{splitted_b[0],splitted_b[1]}
			};
			for(int q = 0; q < 2;++q) //segments
				for(int k = 0; k < transfer.size();++k)
				{
					int size = (int)copy[k + q*transfer.size()].size() / transfer[k].GetBytesSize();
					if( size ) for(int l = 0; l < 2; ++l) //two parts
					{
						splitted[q][l].SetDataSize(transfer[k],size);
						splitted[q][l].SetData(transfer[k],0,size,&copy[k + q*transfer.size()][0]);
					}
				}
	
		}
		*/
	}
	


	void intersect_naive(Mesh * m, ElementArray<Edge> & segments, ElementArray<Node> & nodes, std::vector<Tag> & transfer, Tag pnt, const Unproject & unp, bool print)
	{
		//Tag pnt = m->CreateTag("PROJ_PNT"+m->GetLocalProcessorRank(),DATA_REAL,NODE,NODE,2);
		std::set<Point> intersections;
		std::vector<HandleType> initials(segments.size()*2);
		MarkerType initial = m->CreatePrivateMarker();
		for (int k = 0; k < (int)segments.size(); ++k)
		{
			initials[k*2+0] = segments[k]->getBeg()->GetHandle();
			initials[k*2+1] = segments[k]->getEnd()->GetHandle();
			segments[k]->getBeg()->SetPrivateMarker(initial);
			segments[k]->getEnd()->SetPrivateMarker(initial);
			intersections.insert(make_point(segments[k]->getBeg(),pnt));
			intersections.insert(make_point(segments[k]->getEnd(),pnt));
		}
		for(int i = 0; i < (int)segments.size(); ++i)
		{
			for(int j = i+1; j < (int)segments.size(); ++j)
			{
				std::pair<bool,Point> I = intersect_segments(m,segments[i],segments[j],intersections,pnt,unp,print);
				if( I.first )
				{
					ElementArray<Edge> splitted_a, splitted_b;
					split_edges(m,Node(m,I.second.node),segments[i],segments[j],splitted_a,splitted_b,transfer,print);
					segments[i] = splitted_a[0];
					segments[j] = splitted_b[0];
					segments.push_back(splitted_a[1]);
					segments.push_back(splitted_b[1]);
				}
			}
		}
		//nodes.clear();
		for(std::set<Point>::iterator it = intersections.begin(); it != intersections.end(); ++it)
		{
			if( !m->GetPrivateMarker(it->node,initial) )
				nodes.push_back(it->node);
		}
		for(int k = 0; k < (int)initials.size(); ++k)
			m->RemPrivateMarker(initials[k],initial);
		m->ReleasePrivateMarker(initial);
	}

	void split_segments(Mesh * m, ElementArray<Edge> & segments, ElementArray<Node> & nodes, std::vector<Tag> & transfer, Tag pnt, const Unproject & unp, MarkerType mrk, bool print)
	{
		int i = 0;
		//MarkerType mrk = m->CreateMarker();
		while( i < segments.size() )
		{
			
			//skip edges along pillars
			if( !segments[i]->GetPrivateMarker(mrk) ) 
			{
				i++;
				continue;
			}
			
			//get projected center point of segment
			//vl - point on the middle
			//vu - unprojected middle point onto arc
			Storage::real vu[3], vl[3];
			for(int k = 0; k < 3; ++k)
				vl[k] = (segments[i]->getBeg()->Coords()[k]+segments[i]->getEnd()->Coords()[k])*0.5;
			Point pb = make_point(segments[i]->getBeg(),pnt);
			Point pe = make_point(segments[i]->getEnd(),pnt);
			Point pm(InvalidHandle(),(pb.x+pe.x)*0.5,(pb.y+pe.y)*0.5);
			unp.Act(pm,vu);
			// l - length of the segment
			// h - distance to the arc
			Storage::real l = segments[i]->Length();
			Storage::real h = sqrt((vl[0]-vu[0])*(vl[0]-vu[0])+(vl[1]-vu[1])*(vl[1]-vu[1])+(vl[2]-vu[2])*(vl[2]-vu[2]));
			if( h/l > 0.05 ) // it is rather curvy
			{
				Node n = m->CreateNode(vu);
				n->RealArray(pnt)[0] = pm.x;
				n->RealArray(pnt)[1] = pm.y;
				ElementArray<Edge> splitted_a;
				split_edge(m,n,segments[i],splitted_a,transfer,print);
				splitted_a[0]->SetPrivateMarker(mrk);
				splitted_a[1]->SetPrivateMarker(mrk);
				segments[i] = splitted_a[0];
				segments.push_back(splitted_a[1]);
				nodes.push_back(n);
			}
			else i++;
		}
		//segments.RemMarker(mrk);
		//m->ReleaseMarker(mrk);
	}


	void split_segments_tri(Mesh * m, ElementArray<Edge> & segments, ElementArray<Node> & nodes, std::vector<Tag> & transfer, Tag pnt, const Unproject & unp, MarkerType mrk, bool print)
	{
		int i = 0;
		while( i < segments.size() )
		{
			
			//skip edges along pillars
			if( !segments[i]->GetPrivateMarker(mrk) ) 
			{
				i++;
				continue;
			}
			
			//get projected center point of segment
			//vl - point on the middle
			//vu - unprojected middle point onto arc
			Storage::real vu[3], vl[3];//, vt[3];
			Point pb = make_point(segments[i]->getBeg(),pnt);
			Point pe = make_point(segments[i]->getEnd(),pnt);
			std::pair<Storage::real,Storage::real> tt = intersect_pairs(pb,pe,Point(InvalidHandle(),0,0),Point(InvalidHandle(),1,1));
			Storage::real t = tt.first;
			if( t > -0.5 )
			{
				Point pm(InvalidHandle(),pb.x*t+pe.x*(1-t),pb.y*t+pe.y*(1-t));
				unp.ActTri(pm,vu);
				for(int k = 0; k < 3; ++k)
					vl[k] = segments[i]->getBeg()->Coords()[k]*t+segments[i]->getEnd()->Coords()[k]*(1-t);
				// l - length of the segment
				// h - distance to the arc
				Storage::real l = segments[i]->Length();
				Storage::real h = sqrt((vl[0]-vu[0])*(vl[0]-vu[0])+(vl[1]-vu[1])*(vl[1]-vu[1])+(vl[2]-vu[2])*(vl[2]-vu[2]));
				if( h/l > 1.0e-2 ) // it is not on the same plane
				{
					//unp.ActTri(pm,vt);
					Node n = m->CreateNode(vu);
					n->RealArray(pnt)[0] = pm.x;
					n->RealArray(pnt)[1] = pm.y;
					ElementArray<Edge> splitted_a;
					split_edge(m,n,segments[i],splitted_a,transfer,print);
					splitted_a[0]->SetPrivateMarker(mrk);
					splitted_a[1]->SetPrivateMarker(mrk);
					segments[i] = splitted_a[0];
					segments.push_back(splitted_a[1]);
					nodes.push_back(n);
				}
				else i++;
			}
			else i++;
		}
	}

	
	void block_number_union(Element n, ElementArray<Element> & adj, Tag block, Tag write)
	{
		Storage::integer_array bn = n->IntegerArray(write);
		if( !adj.empty() )
		{
			std::vector<int> uni, tmp;
			for(ElementArray<Edge>::size_type k = 0; k < adj.size(); ++k)
			{
				Storage::integer_array be = adj[k]->IntegerArray(block);
				tmp.resize(uni.size()+be.size());
				tmp.resize(std::set_union(uni.begin(),uni.end(),be.begin(),be.end(),tmp.begin())-tmp.begin());
				uni.swap(tmp);
			}
			bn.replace(bn.begin(),bn.end(),uni.begin(),uni.end());
		}
		else bn.clear();
	}
	
	void block_number_union_merge(Element n, ElementArray<Element> & adj, Tag block, Tag write)
	{
		Storage::integer_array bn = n->IntegerArray(write);
		if( !adj.empty() )
		{
			std::vector<int> uni(bn.begin(),bn.end()), tmp(uni.size());
			for(ElementArray<Edge>::size_type k = 0; k < adj.size(); ++k)
			{
				Storage::integer_array be = adj[k]->IntegerArray(block);
				tmp.resize(uni.size()+be.size());
				tmp.resize(std::set_union(uni.begin(),uni.end(),be.begin(),be.end(),tmp.begin())-tmp.begin());
				uni.swap(tmp);
			}
			bn.replace(bn.begin(),bn.end(),uni.begin(),uni.end());
		}
	}

	void block_number_intersection(Element n,const ElementArray<Element> & adj, Tag block, Tag write)
	{
		Storage::integer_array bn = n->IntegerArray(write);
		if( !adj.empty() )
		{
			Storage::integer_array be = adj[0]->IntegerArray(block);
			std::vector<int> inter(be.begin(),be.end()), tmp(inter.size());
			for(ElementArray<Edge>::size_type k = 1; k < adj.size(); ++k)
			{
				be = adj[k]->IntegerArray(block);
				tmp.resize(std::set_intersection(inter.begin(),inter.end(),be.begin(),be.end(),tmp.begin())-tmp.begin());
				inter.swap(tmp);
			}
			bn.replace(bn.begin(),bn.end(),inter.begin(),inter.end());
		}
		else bn.clear();
		//Storage::integer_array bn = n->IntegerArray(block);
		//bn.replace(bn.begin(),bn.end(),inter.begin(),inter.end());
	}

	
	void block_number_intersection(const ElementArray<Element> & adj, Tag block, std::vector<int> & out)
	{
		out.clear();
		if( !adj.empty() )
		{
			Storage::integer_array be = adj[0]->IntegerArray(block);
			std::vector<int> inter(be.begin(),be.end()), tmp(inter.size());
			for(ElementArray<Edge>::size_type k = 1; k < adj.size(); ++k)
			{
				be = adj[k]->IntegerArray(block);
				tmp.resize(std::set_intersection(inter.begin(),inter.end(),be.begin(),be.end(),tmp.begin())-tmp.begin());
				inter.swap(tmp);
			}
			out.insert(out.end(),inter.begin(),inter.end());
		}
		//Storage::integer_array bn = n->IntegerArray(block);
		//bn.replace(bn.begin(),bn.end(),inter.begin(),inter.end());
	}
	

	void Mesh::LoadECL(std::string File)
	{
		std::cout << std::scientific;
		bool perform_splitting = false;
		bool curvilinear_edges = true;
		bool triangulated_edges = false; //this is not working, yet
		bool project_perm = true;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < file_options.size(); ++k)
		{
			if( file_options[k].first == "ECL_SPLIT_GLUED" )
			{
				if( file_options[k].second == "TRUE" )
					perform_splitting = true;
				else
					perform_splitting = false;
			}
			if( file_options[k].first == "ECL_CURVILINEAR" )
			{
				if( file_options[k].second == "TRUE" )
					curvilinear_edges = true;
				else
					curvilinear_edges = false;
			}
			if( file_options[k].first == "ECL_PROJECT_PERM" )
			{
				if( file_options[k].second == "TRUE" )
					project_perm = true;
				else
					project_perm = false;
			}
		}
		std::vector<HandleType> old_nodes(NumberOfNodes());
		{
			unsigned qq = 0;
			for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it)
				old_nodes[qq++] = *it;
		}
		if( !old_nodes.empty() ) 
			std::sort(old_nodes.begin(),old_nodes.end(),CentroidComparator(this));

		FILE * f = fopen(File.c_str(),"r");
		if( f == NULL )
		{
			std::cout << __FILE__ << ":" << __LINE__ << " cannot open file " << File << std::endl;
			throw BadFileName;
		}
		std::vector< std::pair< std::pair<FILE *,std::string>, int> > fs(1,std::make_pair(std::make_pair(f,File),0));
		char readline[2048], *p, *pend, rec[2048];
		int text_end, text_start, state = ECL_NONE, nchars;
		int waitlines = 0;
		int have_dimens = 0, totread, downread, numrecs, offset;
		int gtype = ECL_GTYPE_NONE;
		int argtype = ECL_VAR_NONE;
		int radial = ECL_GTYPE_NONE;
		Storage::real * read_arrayf = NULL;
		Storage::integer * read_arrayi = NULL;
		Storage::integer dims[3], mapaxis[6] = {0,1,0,0,1,0};
		Storage::real inrad = 0;
		std::vector<Storage::real> xyz,perm,poro, tops,zcorn;
		std::vector<Storage::integer> actnum, satnum;
		compdat_wells wells;
		while(!fs.empty())
		{
			while(fgets(readline,2048,fs.back().first.first) != NULL)
			{
				fs.back().second++; //line number
				{
					if( readline[strlen(readline)-1] == '\n' ) readline[strlen(readline)-1] = '\0';
					text_end = static_cast<int>(strlen(readline));
					for(text_start = 0; isspace(readline[text_start]) && text_start < text_end; text_start++);
					if( text_start == text_end ) continue;
					for(text_end = text_end-1; isspace(readline[text_end]) && text_end > text_start; text_end--);
					readline[text_end+1] = '\0';
					p = readline + text_start;
					pend = readline + text_end + 1;
					for(char * q = p; q < pend; q++) *q = toupper(*q);
				}
				if( p[0] == '-' && p[1] == '-' ) continue; //skip comment
				if(waitlines) {waitlines--; continue;} //skip meaningful lines
				switch(state)
				{
				case ECL_NONE:
					if( !ECLSTRCMP(p,"END") ) //end of data - don't read beyond
					{
						goto ecl_exit_loop;
					}
					else if( !ECLSTRCMP(p,"INCLUDE") ) state = ECL_INCLUDE;
					else if( !ECLSTRCMP(p,"DIMENS") || !ECLSTRCMP(p,"SPECGRID") )
					{
						read_arrayi = dims;
						numrecs = 1;
						downread = totread = 3;
						argtype = ECL_VAR_INT;
						offset = state = ECL_DIMENS;
						have_dimens = 1;
					}
					else if( !ECLSTRCMP(p,"MAPAXIS") )
					{
						read_arrayi = mapaxis;
						numrecs = 1;
						downread = totread = 6;
						argtype = ECL_VAR_INT;
						offset = state = ECL_MAPAXIS;
					}
					else if( !ECLSTRCMP(p,"DX") )
					{
						assert(have_dimens);
						if( xyz.empty() ) xyz.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = xyz.empty()? NULL : &xyz[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = state = ECL_DX;
					}
					else if( !ECLSTRCMP(p,"DY") )
					{
						assert(have_dimens);
						if( xyz.empty() ) xyz.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = xyz.empty() ? NULL : &xyz[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = ECL_DX;
						state = ECL_DY;
					}
					else if( !ECLSTRCMP(p,"DZ") )
					{
						assert(have_dimens);
						if( xyz.empty() ) xyz.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = xyz.empty() ? NULL : &xyz[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = ECL_DX;
						state = ECL_DZ;
					}
					else if( !ECLSTRCMP(p,"COORD") )
					{
						assert(have_dimens);
						if( xyz.empty() ) xyz.resize(3*2*(dims[0]+1)*(dims[1]+1));
						read_arrayf = xyz.empty() ? NULL : &xyz[0];
						numrecs = 1;
						downread = totread = 3*2*(dims[0]+1)*(dims[1]+1);
						argtype = ECL_VAR_REAL;
						offset = state = ECL_COORDS;
						gtype = ECL_GTYPE_ZCORN;
					}
					else if( !ECLSTRCMP(p,"ZCORN") )
					{
						assert(have_dimens);
						if( zcorn.empty() ) zcorn.resize(dims[0]*dims[1]*dims[2]*8);
						read_arrayf = zcorn.empty() ? NULL : &zcorn[0];
						numrecs = 1;
						downread = totread = dims[0]*dims[1]*dims[2]*8;
						argtype = ECL_VAR_REAL;
						state = offset = ECL_ZCORN;
					}
					else if( !ECLSTRCMP(p,"TOPS") )
					{
						assert(have_dimens);
						tops.resize(dims[0]*dims[1]*dims[2]);
						read_arrayf = tops.empty() ? NULL : &tops[0];
						numrecs = 1;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = state = ECL_TOPS;
						gtype = ECL_GTYPE_TOPS;
					}
					else if( !ECLSTRCMP(p,"PERMX") )
					{
						assert(have_dimens);
						if( perm.empty() ) perm.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = perm.empty() ? NULL : &perm[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = state = ECL_PERMX;
					}
					else if( !ECLSTRCMP(p,"PERMY") )
					{
						assert(have_dimens);
						if( perm.empty() ) perm.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = perm.empty() ? NULL : &perm[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = ECL_PERMX;
						state = ECL_PERMY;
					}
					else if( !ECLSTRCMP(p,"PERMZ") )
					{
						assert(have_dimens);
						if( perm.empty() ) perm.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = perm.empty() ? NULL : &perm[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = ECL_PERMX;
						state = ECL_PERMZ;
					}
					else if( !ECLSTRCMP(p,"PORO") )
					{
						assert(have_dimens);
						poro.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = poro.empty() ? NULL : &poro[0];
						numrecs = 1;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = state = ECL_PORO;
					}
					else if( !ECLSTRCMP(p,"ACTNUM") )
					{
						assert(have_dimens);
						actnum.resize(dims[0]*dims[1]*dims[2]);
						read_arrayi = actnum.empty() ? NULL : &actnum[0];
						numrecs = 1;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_INT;
						offset = state = ECL_ACTNUM;
					}
					else if( !ECLSTRCMP(p,"SATNUM") )
					{
						assert(have_dimens);
						satnum.resize(dims[0]*dims[1]*dims[2]);
						read_arrayi = satnum.empty() ? NULL : &satnum[0];
						numrecs = 1;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_INT;
						offset = state = ECL_SATNUM;
					}
					else if( !ECLSTRCMP(p,"RADIAL") )
					{
						radial = ECL_GTYPE_RADIAL;
					}
					else if( !ECLSTRCMP(p,"CART") )
					{
						radial = ECL_GTYPE_CARTESIAN;
					}
					else if( !ECLSTRCMP(p,"COMPDAT") ) //some well data
					{
						state = ECL_COMPDAT;
					}
					else if( !ECLSTRCMP(p,"INRAD") )
					{
						if( radial != ECL_GTYPE_RADIAL ) 
						{
							std::cout << __FILE__ << ":" << __LINE__ << " inner radius specified for cartesian grid ";
							std::cout << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
						}
						if( radial == ECL_GTYPE_NONE ) radial = ECL_GTYPE_RADIAL;
						state = ECL_INRAD;
					}
					else
					{
						//std::cout << __FILE__ << ":" << __LINE__ << " skipped " << p << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
					}
					break;
				case ECL_SKIP_SECTION:
					if( *(pend-1) == '/' || *p == '/' ) state = ECL_NONE;
					break;
				case ECL_INCLUDE:
					if( 1 == sscanf(p," %s",rec) )
					{
						int shift_one = 0;
						if( (rec[0] == '\'' || rec[0] == '"') && rec[0] == rec[strlen(rec)-1] ) //remove quotes
						{
							rec[strlen(rec)-1] = '\0';
							shift_one = 1;
						}
						f = fopen((GetFolder(fs.back().first.second) + "/" + std::string(rec+shift_one)).c_str(),"r");
						if( f == NULL )
						{
							std::cout << __FILE__ << ":" << __LINE__ << " cannot open file " << (GetFolder(fs.back().first.second) + "/" + std::string(rec+shift_one)) << " included from "<< fs.back().first.second << " line " << fs.back().second << std::endl;
							throw BadFileName;
						}
						fs.push_back(std::make_pair(std::make_pair(f,GetFolder(fs.back().first.second) + "/"+std::string(rec+shift_one)),0));
						if( *(pend-1) == '/' ) state = ECL_NONE; else state = ECL_SKIP_SECTION;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " cannot read file name, string " << p << " in " << fs.back().first.second << " line " << fs.back().second << std::endl;
						throw BadFile;
					}
					break;
				case ECL_COMPDAT:
					{
						char wellname[4096];
						char openshut[4096];
						int ii,jj,kk1,kk2;
						if( 6 == sscanf(p,"%s %d %d %d %d %s",wellname,&ii,&jj,&kk1,&kk2,openshut) )
						{
							compdat_entry entry;
							entry.i = ii;
							entry.j = jj;
							entry.k1 = kk1;
							entry.k2 = kk2;
							if( std::string(openshut) == "OPEN" )
								entry.open = true;
							else
								entry.open = false;
							wells[std::string(wellname)].push_back(entry);
						}
					}
					break;
				case ECL_ZCORN:
				case ECL_COORDS:
				case ECL_MAPAXIS:
				case ECL_DIMENS:
				case ECL_DX:
				case ECL_DY:
				case ECL_DZ:
				case ECL_PERMX:
				case ECL_PERMY:
				case ECL_PERMZ:
				case ECL_PORO:
				case ECL_ACTNUM:
				case ECL_SATNUM:
				case ECL_TOPS:
					while( downread > 0 && p < pend )
					{
						if( 1 == sscanf(p,"%s%n",rec,&nchars) )
						{
							p += nchars;
							while(isspace(*p) && p < pend) ++p;
							int count = 1;
							int begval = 0;
							for(int q = 0; q < static_cast<int>(strlen(rec)); ++q)
								if( rec[q] == '*' )
								{
									begval = q+1;
									rec[q] = '\0';
									break;
								}
							if( begval > 0 ) count = atoi(rec);
							if( argtype == ECL_VAR_REAL )
							{
								Storage::real val = atof(rec+begval);
								while(count)
								{
									read_arrayf[numrecs*(totread-(downread--))+(state-offset)] = val;
									count--;
								}
							}
							else if( argtype == ECL_VAR_INT )
							{
								Storage::integer val = atoi(rec+begval);
								while(count)
								{
									read_arrayi[numrecs*(totread-(downread--))+(state-offset)] = val;
									count--;
								}
							}
							else
							{
								std::cout << __FILE__ << ":" << __LINE__ << " probably forgot to set up argument type to read " << std::endl;
								throw Impossible;
							}
						}
						else
						{
							std::cout << __FILE__ << ":" << __LINE__ << " cannot read data " << p << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
							throw BadFile;
						}
						if( *p == '/' ) 
						{
							if( downread > 0 )
							{
								std::cout << __FILE__ << ":" << __LINE__ << " early data termination, read " << totread-downread << " of " << totread << " records in " << fs.back().first.second << ":" << fs.back().second << std::endl;
								throw BadFile;
							}
						}
					}
					if( *(pend-1) == '/' || *p == '/' ) state = ECL_NONE; else if( downread == 0 ) state = ECL_SKIP_SECTION;
					break;
				case ECL_INRAD:
					if( 1 == sscanf(p,"%lf%n",&inrad,&nchars) )
					{
						p += nchars;
						while(isspace(*p) && p < pend) ++p;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " cannot read data " << p << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
						throw BadFile;
					}
					if( *(pend-1) == '/' || *p == '/' ) state = ECL_NONE; else state = ECL_SKIP_SECTION;
					break;
				}
			}
ecl_exit_loop:
			fclose(fs.back().first.first);
			fs.pop_back();
		}
		if( radial == ECL_GTYPE_RADIAL )
		{
			std::cout << __FILE__ << ":" << __LINE__ << " radial grids not supported yet " << std::endl;
		}
		int beg_dims[2], end_dims[2];
		beg_dims[0] = 0;
		beg_dims[1] = 0;
		end_dims[0] = dims[0];
		end_dims[1] = dims[1];
		std::vector<HandleType> blocks(dims[0]*dims[1]*dims[2],InvalidHandle()); // for wells
		if( gtype == ECL_GTYPE_TOPS )
		{
			
			std::vector<HandleType> newnodes((dims[0]+1)*(dims[1]+1)*(dims[2]+1));
			Storage::real x, y, z, node_xyz[3];
			x = 0.0;
			int numnode = 0;
			for(int i = 0; i < dims[0]+1; i++)
			{
				Storage::integer pif = std::min(dims[0]-1,i), pib = std::max(i-1,0);
				y = 0.0;
				for(int j = 0; j < dims[1]+1; j++)
				{
					Storage::integer pjf = std::min(dims[1]-1,j), pjb = std::max(j-1,0);
					z = (
							tops[ECL_IJK_DATA(pib,pjb,0)]+
							tops[ECL_IJK_DATA(pib,pjf,0)]+
							tops[ECL_IJK_DATA(pif,pjb,0)]+
							tops[ECL_IJK_DATA(pif,pjf,0)]
						  )*0.25;
					z -= (
							xyz[3*ECL_IJK_DATA(pib,pjb,0)+2]+
							xyz[3*ECL_IJK_DATA(pib,pjf,0)+2]+
							xyz[3*ECL_IJK_DATA(pif,pjb,0)+2]+
							xyz[3*ECL_IJK_DATA(pif,pjf,0)+2]
							)*0.25;
					for(int k = 0; k < dims[2]+1; k++)
					{
						Storage::integer pkf = std::min(dims[2]-1,k), pkb = std::max(k-1,0);
						bool create = true;
						if( !actnum.empty() )
						{
							create = false;
							for(int l = 0; l < 8; ++l)
							{
								int bi = i - (l%2);
								int bj = j - (l/2%2);
								int bk = k - (l/4%2);
								if( (bi >= 0 && bi < dims[0]) && 
									(bj >= 0 && bj < dims[1]) && 
									(bk >= 0 && bk < dims[2]) )
									create |= (actnum[ECL_IJK_DATA(bi,bj,bk)] != 0);
							}
						}
						if( create )
						{
							
							node_xyz[0] = x;
							node_xyz[1] = y;
							node_xyz[2] = z;
							int find = -1;
							if( !old_nodes.empty() )
							{
								std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),node_xyz,CentroidComparator(this));
								if( it != old_nodes.end() ) 
								{
									Storage::real_array c = RealArrayDF(*it,CoordsTag());
									if( CentroidComparator(this).Compare(node_xyz,c.data()) == 0 )
										find = static_cast<int>(it - old_nodes.begin());
								}
							}
							if( find == -1 ) 
								newnodes[numnode++] = CreateNode(node_xyz)->GetHandle();
							else newnodes[numnode++] = old_nodes[find];
						}
						else numnode++;
						//std::cout << i << " " << j << " " << k << " ( " << x << " , " << y << " , " << z << ") " << newnodes.back()->LocalID() << std::endl; 
						x += (
								(
								  xyz[3*ECL_IJK_DATA(pib,pjb,pkf)+0]+
								  xyz[3*ECL_IJK_DATA(pib,pjf,pkf)+0]+
								  xyz[3*ECL_IJK_DATA(pif,pjb,pkf)+0]+
								  xyz[3*ECL_IJK_DATA(pif,pjf,pkf)+0]
								)
							-
								(
								  xyz[3*ECL_IJK_DATA(pib,pjb,pkb)+0]+
								  xyz[3*ECL_IJK_DATA(pib,pjf,pkb)+0]+
								  xyz[3*ECL_IJK_DATA(pif,pjb,pkb)+0]+
								  xyz[3*ECL_IJK_DATA(pif,pjf,pkb)+0]
								)
								)*0.25;
						y += (
								(
							      xyz[3*ECL_IJK_DATA(pib,pjb,pkf)+1]+
								  xyz[3*ECL_IJK_DATA(pib,pjf,pkf)+1]+
								  xyz[3*ECL_IJK_DATA(pif,pjb,pkf)+1]+
								  xyz[3*ECL_IJK_DATA(pif,pjf,pkf)+1]
								)
							-
								(
								  xyz[3*ECL_IJK_DATA(pib,pjb,pkb)+1]+
								  xyz[3*ECL_IJK_DATA(pib,pjf,pkb)+1]+
								  xyz[3*ECL_IJK_DATA(pif,pjb,pkb)+1]+
								  xyz[3*ECL_IJK_DATA(pif,pjf,pkb)+1]
								)
								)*0.25;
						z += (
								xyz[3*ECL_IJK_DATA(pib,pjb,pkb)+2]+
								xyz[3*ECL_IJK_DATA(pib,pjf,pkb)+2]+
								xyz[3*ECL_IJK_DATA(pif,pjb,pkb)+2]+
								xyz[3*ECL_IJK_DATA(pif,pjf,pkb)+2]+
							    xyz[3*ECL_IJK_DATA(pib,pjb,pkf)+2]+
								xyz[3*ECL_IJK_DATA(pib,pjf,pkf)+2]+
								xyz[3*ECL_IJK_DATA(pif,pjb,pkf)+2]+
								xyz[3*ECL_IJK_DATA(pif,pjf,pkf)+2]
								)*0.125;
					}
					y += (
							xyz[3*ECL_IJK_DATA(pib,pjb,0)+1]+
							xyz[3*ECL_IJK_DATA(pif,pjb,0)+1]+
							xyz[3*ECL_IJK_DATA(pib,pjf,0)+1]+
							xyz[3*ECL_IJK_DATA(pif,pjf,0)+1]
							)*0.25; 
				}
				x += (
						xyz[3*ECL_IJK_DATA(pib,0,0)+0]+
						xyz[3*ECL_IJK_DATA(pif,0,0)+0]
						)*0.5; 
			}
			Tag tagporo,tagperm, tagsatnum;
			Tag block_index = CreateTag("BLOCK_IJK",DATA_INTEGER,CELL,NONE,3);
			if( !poro.empty() ) tagporo = CreateTag("PORO",DATA_REAL,CELL,NONE,1);
			if( !perm.empty() ) tagperm = CreateTag("PERM",DATA_REAL,CELL,NONE,project_perm ? 6 : 3);
			if( !satnum.empty() ) tagsatnum = CreateTag("SATNUM",DATA_INTEGER,CELL,NONE,1);
			const Storage::integer nvf[24] = 
			{ 
			2, 3, 1, 0, //bottom face z-
			4, 5, 7, 6, //top face    z+
			0, 1, 5, 4, //front face  y-
			3, 2, 6, 7, //back face   y+ 
			2, 0, 4, 6, //left face   x-
			1, 3, 7, 5  //right face  x+
			};
			const Storage::integer numnodes[6] = { 4, 4, 4, 4, 4, 4 };
			rMatrix UVW(3,3), U(3,3), S(3,3), Sinv(3,3), V(3,3), K(3,3);
			for(int i = 0; i < dims[0]; i++)
			{
				for(int j = 0; j < dims[1]; j++)
					for(int k = 0; k < dims[2]; k++) if(actnum.empty() || actnum[ECL_IJK_DATA(i,j,k)])
					{
						HandleType verts[8];
						verts[0] = newnodes[((i+0)*(dims[1]+1)+(j+0))*(dims[2]+1)+(k+0)];
						verts[1] = newnodes[((i+1)*(dims[1]+1)+(j+0))*(dims[2]+1)+(k+0)];
						verts[2] = newnodes[((i+0)*(dims[1]+1)+(j+1))*(dims[2]+1)+(k+0)];
						verts[3] = newnodes[((i+1)*(dims[1]+1)+(j+1))*(dims[2]+1)+(k+0)];
						verts[4] = newnodes[((i+0)*(dims[1]+1)+(j+0))*(dims[2]+1)+(k+1)];
						verts[5] = newnodes[((i+1)*(dims[1]+1)+(j+0))*(dims[2]+1)+(k+1)];
						verts[6] = newnodes[((i+0)*(dims[1]+1)+(j+1))*(dims[2]+1)+(k+1)];
						verts[7] = newnodes[((i+1)*(dims[1]+1)+(j+1))*(dims[2]+1)+(k+1)];
						//for(int q = 0; q < 8; q++)
						//	std::cout << verts[q]->Coords()[0] << " " << verts[q]->Coords()[1] << " " << verts[q]->Coords()[2] << " " << verts[q]->LocalID() << std::endl;
						Cell c = CreateCell(ElementArray<Node>(this,verts,verts+8),nvf,numnodes,6).first;
						blocks[ECL_IJK_DATA(i,j,k)] = c->GetHandle();
						integer_array index = c->IntegerArray(block_index);
						index[0] = i;
						index[1] = j;
						index[2] = k;
						if( !poro.empty() ) c->RealDF(tagporo) = poro[ECL_IJK_DATA(i,j,k)];
						if( !satnum.empty() ) c->IntegerDF(tagsatnum) = satnum[ECL_IJK_DATA(i,j,k)];
						if( !perm.empty() )
						{
							Storage::real_array arr_perm = c->RealArrayDF(tagperm);
							if( project_perm )
							{
								ElementArray<Face> faces = c->getFaces();
								Storage::real cntp[3], cntn[3], u[3], L;
								//compute grid block coordinates
								for(int r = 0; r < 3; ++r)
								{
									faces[5-r*2]->Centroid(cntp);
									faces[4-r*2]->Centroid(cntn);
									u[0] = cntp[0] - cntn[0];
									u[1] = cntp[1] - cntn[1];
									u[2] = cntp[2] - cntn[2];
									L = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
									if( L ) u[0] /= L, u[1] /= L, u[2] /= L;
									//fill row
									UVW(r,0) = u[0];
									UVW(r,1) = u[1];
									UVW(r,2) = u[2];

								}
								UVW.SVD(U,S,V);
								for(int r = 0; r < 3; ++r)
								{
									if( S(r,r) > 0.0 )
										Sinv(r,r) = 1.0/S(r,r);
									else 
										Sinv(r,r) = 0.0;
								}
								K.Zero();
								K(0,0) = perm[3*(ECL_IJK_DATA(i,j,k))+0];
								K(1,1) = perm[3*(ECL_IJK_DATA(i,j,k))+1];
								K(2,2) = perm[3*(ECL_IJK_DATA(i,j,k))+2];
								
								K = V*Sinv*U.Transpose()*K*U*Sinv*V.Transpose();

								arr_perm[0] = K(0,0);
								arr_perm[1] = K(0,1);
								arr_perm[2] = K(0,2);
								arr_perm[3] = K(1,1);
								arr_perm[4] = K(1,2);
								arr_perm[5] = K(2,2);
							}
							else
							{
								arr_perm[0] = perm[3*(ECL_IJK_DATA(i,j,k))+0];
								arr_perm[1] = perm[3*(ECL_IJK_DATA(i,j,k))+1];
								arr_perm[2] = perm[3*(ECL_IJK_DATA(i,j,k))+2];
							}
						}
					}
			}
		}
		else if( gtype == ECL_GTYPE_ZCORN )
		{
#if defined(USE_MPI)
			{
				std::cout << "x: " << beg_dims[0] << " " << end_dims[0] << std::endl;
				std::cout << "y: " << beg_dims[1] << " " << end_dims[1] << std::endl;
				int mpi_dims[2] = {0,0};
				MPI_Dims_create(GetProcessorsNumber(),2,mpi_dims);
				std::cout << "mpi: " << mpi_dims[0] << " " << mpi_dims[1] << std::endl;
				int rank = GetProcessorRank();
				//my position on the grid
				int rankx = rank % mpi_dims[0];
				int ranky = rank / mpi_dims[0];
				//my subgrid
				int numx = dims[0]/mpi_dims[0];
				int numy = dims[1]/mpi_dims[1];
				std::cout << "num: " << numx << " " << numy << std::endl;
				beg_dims[0] = numx*rankx;
				beg_dims[1] = numy*ranky;
				end_dims[0] = numx*(rankx+1);
				end_dims[1] = numy*(ranky+1);
				if( rankx == mpi_dims[0]-1 ) end_dims[0] = dims[0];
				if( ranky == mpi_dims[1]-1 ) end_dims[1] = dims[1];
				std::cout << "x: " << beg_dims[0] << " " << end_dims[0] << std::endl;
				std::cout << "y: " << beg_dims[1] << " " << end_dims[1] << std::endl;
			}
#endif //USE_MPI


			SetTopologyCheck(PRINT_NOTIFY | NEED_TEST_CLOSURE | PROHIBIT_MULTIPOLYGON | PROHIBIT_MULTILINE | MARK_ON_ERROR);
			SetTopologyCheck(DEGENERATE_EDGE | DEGENERATE_FACE | DEGENERATE_CELL);
			SetTopologyCheck(TRIPLE_SHARED_FACE | FLATTENED_CELL | INTERLEAVED_FACES);
			SetTopologyCheck(DUPLICATE_EDGE | DUPLICATE_FACE | DUPLICATE_CELL);
			SetTopologyCheck(ADJACENT_DUPLICATE | ADJACENT_DIMENSION);
			RemTopologyCheck(THROW_EXCEPTION);
			//actnum.clear();
			if( zcorn.empty() )
			{
				std::cout << "ZCORN was not provided, cannot construct grid" << std::endl;
				throw BadFile;
			}
			if( xyz.empty() )
			{
				std::cout << "COORD was not provided, cannot construct grid" << std::endl;
				throw BadFile;
			}
			Tag cell_number = CreateTag("CELL_NUMBER",DATA_INTEGER,CELL,NONE,1);
			Tag block_index = CreateTag("BLOCK_IJK",DATA_INTEGER,CELL,NONE,3);
			Tag edge_number = CreateTag("EDGE_NUMBER",DATA_INTEGER,EDGE,NONE);
			Tag block_intersection_number = CreateTag("BLOCK_INTERSECTION_NUMBER",DATA_INTEGER,EDGE,NONE);
			Tag block_number = CreateTag("BLOCK_NUMBER",DATA_INTEGER,EDGE | NODE,NONE);
			Tag face_origin = CreateTag("FACE_ORIGIN",DATA_INTEGER,FACE,1);
			Tag pillar_num = CreateTag("PILLAR_NUM",DATA_INTEGER,NODE,NONE,2);
			Tag node_blocks = CreateTag("NODE_BLOCK",DATA_INTEGER,NODE,NONE);
			Tag node_number = CreateTag("NODE_NUMBER",DATA_INTEGER,NODE,NONE);
			Tag node_alpha = CreateTag("NODE_ALPHA",DATA_REAL,NODE,NONE,1);
			Tag node_pos = CreateTag("NODE_POSITION",DATA_INTEGER,NODE,NONE,1);
			Tag block_pair = CreateTag("BLOCK_PAIR",DATA_INTEGER,FACE,NONE,3);
			//typedef std::map<Storage::real,Node,pillar_less> pillar;
			//std::vector< pillar > pillars((dims[0]+1)*(dims[1]+1));
			std::vector< HandleType > block_nodes(dims[0]*dims[1]*dims[2]*8,InvalidHandle());
			//all edges along pillars
			std::vector< ElementArray<Edge> > pillar_edges((dims[0]+1)*(dims[1]+1),ElementArray<Edge>(this));
			//this variant goes over pillars
			printf("started creating nodes and edges on pillars\n");
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				//z-positions for all the nodes
				std::vector< real > pillar_node_depth(dims[2]*8);
				//block number that posses
				std::vector< integer > pillar_node_block(dims[2]*8);
				//position of node relative to block
				std::vector< integer > pillar_node_number(dims[2]*8);
				//position of each block's node in pillar_nodes array
				//std::vector< int > pillar_node_pos(dims[2]*8);
				//number of depths in array
				int num;
				//arrays to sort all the data fields
				std::vector< int > indices_sort(dims[2]*8);
				std::vector< integer > itemporary(dims[2]*8);
				std::vector< real > rtemporary(dims[2]*8);
				//nodes along pillar
				std::vector< HandleType > pillar_nodes(dims[2]*8);
				//number of unique nodes along pillar
				int num_nodes;
				//structure to create an edge from pair of nodes
				ElementArray<Node> edge_nodes(this,2);
				//compute union of block number that belong to edges
				std::vector<int> block_union;
				//pillar min and max positions
				real zmin, zmax, xmin, xmax, ymin, ymax;
#if defined(USE_OMP)
#pragma omp for
#endif
				for(int i = beg_dims[0]; i < end_dims[0]+1; i++)
				{
					for(int j = beg_dims[1]; j < end_dims[1]+1; j++)
					{
						//reset number of nodes
						num_nodes = num = 0;
						//loop over block corners on pillar
						for(int k = 0; k < dims[2]; k++)
						{
							//loop over 8 blocks around node
							for(int l = 0; l < 4; ++l)
							{
								//block number
								int bi = i + (l%2) - 1;
								int bj = j + ((l/2)%2) - 1;
								if( bi >= 0 && bj >= 0 && bi < dims[0] && bj < dims[1] && (actnum.empty() || actnum[ECL_IJK_DATA(bi,bj,k)]) )
								{
									for(int q = 0; q < 2; ++q)
									{
										pillar_node_block[num] = ECL_IJK_DATA(bi,bj,k);
										pillar_node_depth[num] = zcorn[ECL_IJK_ZCORN(bi,bj,k,3-l + q*4)];
										pillar_node_number[num] = 3-l + q*4;
										num++;
									}
								}
							}
						}
						//sort data
						for(int l = 0; l < num; ++l) indices_sort[l] = l;
						std::sort(&indices_sort[0],&indices_sort[num],depth_comparator(&pillar_node_depth[0]));
						for(int l = 0; l < num; ++l) itemporary[l] = pillar_node_block[l];
						for(int l = 0; l < num; ++l) pillar_node_block[l] = itemporary[indices_sort[l]];
						for(int l = 0; l < num; ++l) itemporary[l] = pillar_node_number[l];
						for(int l = 0; l < num; ++l) pillar_node_number[l] = itemporary[indices_sort[l]];
						for(int l = 0; l < num; ++l) rtemporary[l] = pillar_node_depth[l];
						for(int l = 0; l < num; ++l) pillar_node_depth[l] = rtemporary[indices_sort[l]];
						//assert(std::is_sorted(&pillar_node_depth[0],&pillar_node_depth[num]));
						//retrive pillar info
						xmin = xyz[ECL_IJK_COORDS(i,j,0,0)];
						xmax = xyz[ECL_IJK_COORDS(i,j,1,0)];
						ymin = xyz[ECL_IJK_COORDS(i,j,0,1)];
						ymax = xyz[ECL_IJK_COORDS(i,j,1,1)];
						zmin = xyz[ECL_IJK_COORDS(i,j,0,2)];
						zmax = xyz[ECL_IJK_COORDS(i,j,1,2)];
						//create nodes
						for(int l = 0; l < num; ++l)
						{
							real node_xyz[3];
							real mean_depth = pillar_node_depth[l];
							real mean_alpha = (mean_depth-zmin)/(zmax-zmin);
							real mean_num = 1;
							//compute node position
							node_xyz[0] = xmin + mean_alpha * (xmax - xmin);
							node_xyz[1] = ymin + mean_alpha * (ymax - ymin);
							node_xyz[2] = mean_depth;
							//search node among priviously existed nodes, create if not found
							int find = -1;
							if( !old_nodes.empty() )
							{
								std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),node_xyz,CentroidComparator(this));
								if( it != old_nodes.end() )
								{
									Storage::real_array c = RealArrayDF(*it,CoordsTag());
									if( CentroidComparator(this).Compare(node_xyz,c.data()) == 0 )
										find = static_cast<int>(it - old_nodes.begin());
								}
							}
							Node current = (find == -1 ? CreateNode(node_xyz) : Node(this,old_nodes[find]));
							//write down pillar node information (debug)
							integer_array pn = current->IntegerArray(pillar_num);
							pn[0] = i;
							pn[1] = j;
							//to detect edge possession
							integer_array nb = current->IntegerArray(node_blocks);
							// (debug)
							integer_array nn = current->IntegerArray(node_number);
							//remember which block the node belongs
							assert(block_nodes[pillar_node_block[l]*8 + pillar_node_number[l]] == InvalidHandle());
							block_nodes[pillar_node_block[l]*8 + pillar_node_number[l]] = current->GetHandle();
							nb.push_back(pillar_node_block[l]);
							nn.push_back(pillar_node_number[l]);
							//check for nodes sharing the same position
							for(int k = l+1; k < num; ++k)
							{
								if( std::abs(pillar_node_depth[l]-pillar_node_depth[k]) < ECL_PILLAR_DEPTH_EPS )
								{
									mean_depth += pillar_node_depth[l];
									block_nodes[pillar_node_block[k]*8 + pillar_node_number[k]] = current->GetHandle();
									nb.push_back(pillar_node_block[k]);
									nn.push_back(pillar_node_number[k]);
									mean_num++;
									l = k;
								}
								else break;
							}
							//write alpha position
							mean_depth /= mean_num;
							mean_alpha = (mean_depth-zmin)/(zmax-zmin);
							//correct position due to averaging
							if( mean_num > 1 )
							{
								current->Real(node_alpha) = mean_alpha;
								real_array coords = current->Coords();
								coords[0] = xmin + mean_alpha * (xmax - xmin);
								coords[1] = ymin + mean_alpha * (ymax - ymin);
								coords[2] = mean_depth;
							}
							else current->Real(node_alpha) = mean_alpha;
							//remember position of node
							current->Integer(node_pos) = num_nodes;
							//write nodes into array
							pillar_nodes[num_nodes] = current->GetHandle();
							num_nodes++;
						} //l
						//sort node blocks information on each node
						for(int h = 0; h < num_nodes; ++h)
						{
							Storage::integer_array b = IntegerArray(pillar_nodes[h],node_blocks);
							Storage::integer_array e = IntegerArray(pillar_nodes[h],node_number);
							assert(e.size() == b.size());
							//sort indices according to b
							for(int l = 0; l < (int)b.size(); ++l) indices_sort[l] = l;
							std::sort(&indices_sort[0],&indices_sort[b.size()],index_comparator(b));
							//arrange data in b and e arrays according to indices_sort
							//first b array
							for(int l = 0; l < (int)b.size(); ++l) itemporary[l] = b[l];
							for(int l = 0; l < (int)b.size(); ++l) b[l] = itemporary[indices_sort[l]];
							//then e array
							for(int l = 0; l < (int)e.size(); ++l) itemporary[l] = e[l];
							for(int l = 0; l < (int)e.size(); ++l) e[l] = itemporary[indices_sort[l]];
							//assert(std::is_sorted(b.begin(),b.end()));
						} //l
						//create edges
						ElementArray<Edge> & p_edges = pillar_edges[i*(dims[1]+1)+j];
						for(int l = 0; l < num_nodes-1; ++l)
						{
							edge_nodes.at(0) = pillar_nodes[l];
							edge_nodes.at(1) = pillar_nodes[l+1];
							p_edges.push_back(CreateEdge(edge_nodes).first);
						} //l
						//mark edges along pillars (vertical edges, 8, 9, 10, 11)
						for(int k = 0; k < dims[2]; k++)
						{
							//loop over 4 blocks
							for(int l = 0; l < 4; ++l)
							{
								//block number
								int bi = i + (l%2) - 1;
								int bj = j + (l/2) - 1;
								if( bi >= 0 && bj >= 0 && bi < dims[0] && bj < dims[1] && (actnum.empty() || actnum[ECL_IJK_DATA(bi,bj,k)]) )
								{
									assert(block_nodes[ECL_IJK_DATA(bi,bj,k)*8 + (3-l)+0] != InvalidHandle());
									assert(block_nodes[ECL_IJK_DATA(bi,bj,k)*8 + (3-l)+4] != InvalidHandle());
									int pos0 = Integer(block_nodes[ECL_IJK_DATA(bi,bj,k)*8 + (3-l)+0],node_pos); //bottom
									int pos1 = Integer(block_nodes[ECL_IJK_DATA(bi,bj,k)*8 + (3-l)+4],node_pos); //top
									if( pos0 > pos1 ) std::swap(pos0,pos1);
									//edge 11-l
									for(int q = pos0; q < pos1; ++q)
									{
										p_edges[q]->IntegerArray(block_number).push_back(ECL_IJK_DATA(bi,bj,k));
										p_edges[q]->IntegerArray(edge_number).push_back(11-l);
									} //q
								}
							} //l
						} //k
						//erase edges on pillars that are not used in any block
						ElementArray<Edge>::iterator it = p_edges.begin();
						while(it != p_edges.end())
						{
							if( it->IntegerArray(block_number).empty() )
							{
								it->Delete();
								it = p_edges.erase(it);
							}
							else ++it;
						} //it
						//sort block number information on edges
						for(int h = 0; h < (int)p_edges.size(); ++h)
						{
							Storage::integer_array b = p_edges[h].IntegerArray(block_number);
							Storage::integer_array e = p_edges[h].IntegerArray(edge_number);
							assert(e.size() == b.size());
							//sort indices according to b
							for(int l = 0; l < (int)b.size(); ++l) indices_sort[l] = l;
							std::sort(&indices_sort[0],&indices_sort[b.size()],index_comparator(b));
							//arrange data in b and e arrays according to indices_sort
							//first b array
							for(int l = 0; l < (int)b.size(); ++l) itemporary[l] = b[l];
							for(int l = 0; l < (int)b.size(); ++l) b[l] = itemporary[indices_sort[l]];
							//then e array
							for(int l = 0; l < (int)e.size(); ++l) itemporary[l] = e[l];
							for(int l = 0; l < (int)e.size(); ++l) e[l] = itemporary[indices_sort[l]];
							//assert(std::is_sorted(b.begin(),b.end()));
						}
						//find out who shares each node along the pillar
						for(int l = 0; l < num_nodes; ++l)
						{
							Node n(this,pillar_nodes[l]);
							ElementArray<Element> adj = n->getAdjElements(EDGE);
							Storage::integer_array nb = n->IntegerArray(node_blocks);
							Storage::integer_array bn = n->IntegerArray(block_number);
							bn.insert(bn.end(),nb.begin(),nb.end());
							block_number_union_merge(n,adj,block_number,block_number);
						}
					} //i
				} //j
			} //omp parallel
			printf("finished creating nodes and edges on pillars\n");
			DeleteTag(node_pos);
				
				
			/*      (6)*<<------[3]--------*(7)
						/|                  /|
					[6]                 [7]|
					/  |                /  |
				(4)*--------[2]------>>*(5)|
					|  [10]             |  [11]
					|                   |   |
					|   |               |   |
					|                   |   |
					[8]  |              [9]  |
					|(2)*- - - - [1]- - |->>*(3)
					|  /                |  /
					|[4]                |[5]
					|/                  |/
				(0)*<<------[0]--------*(1)
			*/
			//Tag pillar_mark = CreateTag("PILLAR_MARK",DATA_INTEGER,EDGE,NONE,3);
			//tags to be transfered
			std::vector<Tag> transfer(2), transfer2(3);
			transfer[0] = edge_number;
			transfer[1] = block_number;
			transfer2[0] = edge_number;
			transfer2[1] = block_number;
			transfer2[2] = block_intersection_number;
			//transfer[2] = pillar_mark;
			//store faces for each block to assemble block cells later
			std::vector< ElementArray<Face> > block_faces(dims[0]*dims[1]*dims[2],ElementArray<Face>(this));
			//store top-bottom edges to assemble top and bottom faces later
			std::vector< ElementArray<Edge> > block_edges(dims[0]*dims[1]*dims[2],ElementArray<Edge>(this));
			//print out intersection algorithm
			bool print_inter = false;
			//print out block number information
			bool print_bn = false;
			//print out block edges information
			bool print_bedges = false;
			//some intermediate info
			bool print_info = false;
//#if defined(USE_OMP)
//				std::vector<omp_lock_t> block_locks(dims[0]*dims[1]*dims[2]);
//				for(int q = 0; q < dims[0]*dims[1]*dims[2]; ++q) omp_init_lock(&block_locks[q]);
//#endif
			//check_shared_mrk = false;
			//check_private_mrk = false;

			//go over nx pairs of pillars, then ny pairs of pillars
			for(int q = 0; q < 2; ++q)
			{
				printf("started creating faces for pairs of pillar along %s\n",q ? "ny":"nx");
#if defined(USE_OMP)
#pragma omp parallel
#endif
				{
					//structure to create an edge from pair of nodes
					ElementArray<Node> edge_nodes(this,2);
					//set of lines to be intersected
					ElementArray<Edge> edges(this);
					ElementArray<Node> intersections(this);
					//for sorting a pair of data
					std::vector<int> indices_sort, temporary;
					//array to obtain block number
					//std::vector<int> block_number_inter, node_number_inter;
					//store edges along pillar to assemble faces along pillars, on back side of the pillar
					std::vector< ElementArray<Edge> > pillar_block_edges_back(dims[2],ElementArray<Edge>(this));
					//on front side of the pillar
					std::vector< ElementArray<Edge> > pillar_block_edges_front(dims[2],ElementArray<Edge>(this));
					//mark edges on pillar for adjacency retrival
					MarkerType mrk = CreatePrivateMarker();
					//mark original edges of each block face, so that we know outer boundary on constructed interface
					MarkerType outer = CreatePrivateMarker();
					//remember visited edges
					MarkerType visited = CreatePrivateMarker();
					//mark intersection nodes
					MarkerType internodes = CreatePrivateMarker();
					//create tag that will be used for intersection of segments in projected space of [0,1]^2 quad
					Tag pnt;
					{
						std::stringstream tag_name;
						tag_name << "PROJ_PNT_PROC_" << GetLocalProcessorRank();
						pnt = CreateTag(tag_name.str(),DATA_REAL,NODE,NONE,2);
					}
					for(int uneven = 0; uneven < 2; uneven++)
#if defined(USE_OMP)
#pragma omp for
#endif
					//for(int i = 0; i < dims[0]+q; i++)
					for(int jt = uneven+beg_dims[1-q]; jt < end_dims[1-q]+1; jt+=2)
					{
						//printf("%s %6.2f%%\r",q ? "ny":"nx", ((Storage::real)i)/((Storage::real)dims[0]+q-1)*100);
						//fflush(stdout);
						//for(int j = 0; j < dims[1]+(1-q); ++j)
						for(int it = beg_dims[q]; it < end_dims[q]; it++)
						{
							int i,j;
							if( q )
							{
								i = jt;
								j = it;
							}
							else
							{
								i = it;
								j = jt;
							}
							//if( i == 63 && j == 129 && q == 1 ) print_inter = true;
							//if( i == 9 && j == 91 && q == 1 ) 
							//if( i == 7 && j == 78 && q == 1 ) 
							if( false )
							{
								print_inter = true;
								print_info = true;
							}
							if( print_info )
							{
								std::cout << "working on " << (q ? "ny" : "nx") << " pair of pillars: " << i << "," << j << " and " << i-q << "," << j-1+q << std::endl;
							}
							//p0 is at i,j
							//p1 is at i+1,j, when q = 0 and at i,j+1, when q = 1
							//pillar & p0 = pillars[(i+  0)*(dims[1]+1)+j+0];
							//pillar & p1 = pillars[(i+1-q)*(dims[1]+1)+j+q];
							//if( p0.empty() || p1.empty() ) continue;
							//add projected point information
							//i << "," << j << " and " << i+1-q << "," << j+q <<
							Unproject unp(&xyz[ECL_IJK_COORDS(i,j,1,0)],&xyz[ECL_IJK_COORDS(i,j,0,0)],&xyz[ECL_IJK_COORDS(i+1-q,j+q,1,0)],&xyz[ECL_IJK_COORDS(i+1-q,j+q,0,0)]);
							//Unproject unp(p0.begin()->second->Coords().data(),p0.rbegin()->second.Coords().data(),p1.begin()->second->Coords().data(),p1.rbegin()->second.Coords().data());
							//preallocate array for edges
							//edges.reserve(2*std::max(p0.size(),p1.size()));
							//add far faces of j-1 block and near faces of j block
							for(int m = 0; m < 2; ++m)
							{
								if( ((1-q)*j + q*i - dims[1-q]*m)*(1-2*m) > 0 ) //test j when q = 0 and i when q = 1
								{
									int bi = i-(q)*(1-m);
									int bj = j-(1-q)*(1-m);
									if( print_info ) std::cout << (m ?"front":"back") << " cell: " << bi << "," << bj << std::endl;
									for(int k = 0; k < dims[2]; ++k) if( actnum.empty() || actnum[ECL_IJK_DATA(bi,bj,k)] )
									{
										for(int l = 0; l < 2; ++l) //top-bottom
										{
											//q = 0 - test 2,3 + 4*l (edge 1 and 3)
											//q = 1 - test 1,3 + 4*l (edge 5 and 7)
											assert(block_nodes[ECL_IJK_DATA(bi,bj,k)*8 + 2 - 2*m - q*(1-m) + 4*l] != InvalidHandle());
											assert(block_nodes[ECL_IJK_DATA(bi,bj,k)*8 + 3 - 2*m + q*(m)   + 4*l] != InvalidHandle());
											edge_nodes.at(0) = block_nodes[ECL_IJK_DATA(bi,bj,k)*8 + 2 - 2*m - q*(1-m) + 4*l];
											edge_nodes.at(1) = block_nodes[ECL_IJK_DATA(bi,bj,k)*8 + 3 - 2*m + q*(m)   + 4*l];
											if( edge_nodes.at(0) != edge_nodes.at(1) )
											{
												//edge_nodes.SetMarker(used);
												Edge e = CreateEdge(edge_nodes).first;
												if( !e->GetPrivateMarker(visited) )
												{
													edges.push_back(e);
													e->SetPrivateMarker(visited);
												}
												e->IntegerArray(block_number).push_back(ECL_IJK_DATA(bi,bj,k));
												e->IntegerArray(edge_number).push_back(1 - m + 2*l + 4*q);
												real_array pbeg = edge_nodes[0]->RealArray(pnt);
												real_array pend = edge_nodes[1]->RealArray(pnt);
												pbeg[0] = edge_nodes[0]->Real(node_alpha);
												pbeg[1] = 0;
												pend[0] = edge_nodes[1]->Real(node_alpha);
												pend[1] = 1;
												//Storage::integer_array mark = e->IntegerArray(pillar_mark);
												//back block indices
												//mark[0] = i;
												//mark[1] = j;
												//nx or ny
												//mark[2] = q;
											}
										}
									}
								}
							}
							edges.RemPrivateMarker(visited);
							//produce intersected edges
							if( print_inter )
							{
								std::cout << "input edges: " << edges.size() << std::endl;
								if( true ) for(int k = 0; k < edges.size(); ++k)
								{
									Storage::integer_array bn = edges[k]->IntegerArray(block_number);
									Storage::integer_array en = edges[k]->IntegerArray(edge_number);
									std::cout << "edge " << k << " " << edges[k]->GetHandle() << " " << edges[k]->getBeg()->GetHandle() << "<->" << edges[k]->getEnd()->GetHandle() << " blocks: ";
									for(int l = 0; l < (int)bn.size(); ++l)
										std::cout << bn[l] << "(" << bn[l]%dims[0] << "," << bn[l]/dims[0]%dims[1] << "," << bn[l]/dims[0]/dims[1] << "):" << en[l] << " ";
									std::cout << std::endl;
								}

								for(int k = 0; k < edges.size(); ++k)
								{
									std::cout << "(" << edges[k]->getBeg()->Coords()[0] << "," << edges[k]->getBeg()->Coords()[1] << "," << edges[k]->getBeg()->Coords()[2] << ") <-> (" << edges[k]->getEnd()->Coords()[0] << "," << edges[k]->getEnd()->Coords()[1] << "," << edges[k]->getEnd()->Coords()[2] << ")" << std::endl;
								}

								std::cout << "Projected: " << std::endl;

								for(int k = 0; k < edges.size(); ++k)
								{
									std::cout << "(" << edges[k]->getBeg()->RealArray(pnt)[0] << "," << edges[k]->getBeg()->RealArray(pnt)[1] << ",0) <-> (" << edges[k]->getEnd()->RealArray(pnt)[0] << "," << edges[k]->getEnd()->RealArray(pnt)[1] << ",0)" << std::endl;
								}
								/*
								std::cout << "Unprojected: " << std::endl;

								for(int k = 0; k < edges.size(); ++k)
								{
									real v0[3],v1[3];
									Point p0(edges[k]->getBeg()->GetHandle(),edges[k]->getBeg()->RealArray(pnt).data());
									Point p1(edges[k]->getEnd()->GetHandle(),edges[k]->getEnd()->RealArray(pnt).data());
									unp.Act(p0,v0);
									unp.Act(p1,v1);
									std::cout << "(" << v0[0] << "," << v0[1] << "," << v0[2] << ") <-> (" << v1[0] << "," << v1[1] << "," << v1[2] << ")" << std::endl;
								}


								std::cout << "Unprojected segmented: " << std::endl;

								for(int k = 0; k < edges.size(); ++k)
								{
									real v0[3],v1[3];
									Point p0(edges[k]->getBeg()->GetHandle(),edges[k]->getBeg()->RealArray(pnt).data());
									Point p1(edges[k]->getEnd()->GetHandle(),edges[k]->getEnd()->RealArray(pnt).data());
									int nsegments = 5;
									for(int q = 0; q < nsegments; ++q)
									{
										Point pb(InvalidHandle(), q*(p1.x-p0.x)/nsegments+p0.x,q*(p1.y-p0.y)/nsegments+p0.y);
										Point pe(InvalidHandle(), (q+1)*(p1.x-p0.x)/nsegments+p0.x,(q+1)*(p1.y-p0.y)/nsegments+p0.y);
										unp.Act2(pb,v0);
										unp.Act2(pe,v1);
										std::cout << "(" << v0[0] << "," << v0[1] << "," << v0[2] << ") <-> (" << v1[0] << "," << v1[1] << "," << v1[2] << ")" << std::endl;
									}	
								}
								*/
							}
							assert(count_duplicates(edges) == 0);
								
							intersections.clear();
							//intersect(this,edges,intersections,transfer,pnt,unp,print_inter);
							
							intersect_naive(this,edges,intersections,transfer,pnt,unp,print_inter);
							
							assert(count_duplicates(edges) == 0);
							if( print_inter ) std::cout << "intersections: " << intersections.size() << std::endl;
							if( print_inter )
							{

								std::cout << "output edges: " << edges.size() << std::endl;
									
								if( true ) for(int k = 0; k < edges.size(); ++k)
								{
									Storage::integer_array bn = edges[k]->IntegerArray(block_number);
									Storage::integer_array en = edges[k]->IntegerArray(edge_number);
									std::cout << "edge " << k << " " << edges[k]->GetHandle() << " " << edges[k]->getBeg()->GetHandle() << "<->" << edges[k]->getEnd()->GetHandle() << " blocks: ";
									for(int l = 0; l < (int)bn.size(); ++l)
										std::cout << bn[l] << "(" << bn[l]%dims[0] << "," << bn[l]/dims[0]%dims[1] << "," << bn[l]/dims[0]/dims[1] << "):" << en[l] << " ";
									std::cout << std::endl;
								}

									
								for(int k = 0; k < edges.size(); ++k)
								{
									//std::cout << "edge " << k << " " << edges[k]->getBeg()->GetHandle() << "<->" << edges[k]->getEnd()->GetHandle() << std::endl;
									std::cout << "(" << edges[k]->getBeg()->Coords()[0] << "," << edges[k]->getBeg()->Coords()[1] << "," << edges[k]->getBeg()->Coords()[2] << ") <-> (" << edges[k]->getEnd()->Coords()[0] << "," << edges[k]->getEnd()->Coords()[1] << "," << edges[k]->getEnd()->Coords()[2] << ")" << std::endl;

								}
								std::cout << "Projected: " << std::endl;

								for(int k = 0; k < edges.size(); ++k)
								{
									std::cout << "(" << edges[k]->getBeg()->RealArray(pnt)[0] << "," << edges[k]->getBeg()->RealArray(pnt)[1] << ",0) <-> (" << edges[k]->getEnd()->RealArray(pnt)[0] << "," << edges[k]->getEnd()->RealArray(pnt)[1] << ",0)" << std::endl;
								}
								/*
								std::cout << "Unprojected: " << std::endl;

								for(int k = 0; k < edges.size(); ++k)
								{
									real v0[3],v1[3];
									Point p0(edges[k]->getBeg()->GetHandle(),edges[k]->getBeg()->RealArray(pnt).data());
									Point p1(edges[k]->getEnd()->GetHandle(),edges[k]->getEnd()->RealArray(pnt).data());
									unp.Act2(p0,v0);
									unp.Act2(p1,v1);
									std::cout << "(" << v0[0] << "," << v0[1] << "," << v0[2] << ") <-> (" << v1[0] << "," << v1[1] << "," << v1[2] << ")" << std::endl;
								}
								
								std::cout << "Unprojected segmented: " << std::endl;

								for(int k = 0; k < edges.size(); ++k)
								{
									real v0[3],v1[3];
									Point p0(edges[k]->getBeg()->GetHandle(),edges[k]->getBeg()->RealArray(pnt).data());
									Point p1(edges[k]->getEnd()->GetHandle(),edges[k]->getEnd()->RealArray(pnt).data());
									int nsegments = 5;
									for(int q = 0; q < nsegments; ++q)
									{
										Point pb(InvalidHandle(), q*(p1.x-p0.x)/nsegments+p0.x,q*(p1.y-p0.y)/nsegments+p0.y);
										Point pe(InvalidHandle(), (q+1)*(p1.x-p0.x)/nsegments+p0.x,(q+1)*(p1.y-p0.y)/nsegments+p0.y);
										unp.Act2(pb,v0);
										unp.Act2(pe,v1);
										std::cout << "(" << v0[0] << "," << v0[1] << "," << v0[2] << ") <-> (" << v1[0] << "," << v1[1] << "," << v1[2] << ")" << std::endl;
									}	
								}
								*/
							}

							
							//sort block numbers on edges
							for(int k = 0; k < (int)edges.size(); ++k)
							{
								Storage::integer_array b = edges[k].IntegerArray(block_number);
								Storage::integer_array e = edges[k].IntegerArray(edge_number);
								assert(e.size() == b.size());
								//sort indices according to b
								indices_sort.resize(b.size());
								for(int l = 0; l < indices_sort.size(); ++l) indices_sort[l] = l;
								std::sort(indices_sort.begin(),indices_sort.end(),index_comparator(b));
								//arrange data in b and e arrays according to indices_sort
								temporary.resize(b.size());
								//first b array
								for(int l = 0; l < (int)b.size(); ++l) temporary[l] = b[l];
								for(int l = 0; l < (int)b.size(); ++l) b[l] = temporary[indices_sort[l]];
								//then e array
								for(int l = 0; l < (int)e.size(); ++l) temporary[l] = e[l];
								for(int l = 0; l < (int)e.size(); ++l) e[l] = temporary[indices_sort[l]];
								//assert(std::is_sorted(b.begin(),b.end()));
							}
							//mark inter-pillar edges
							edges.SetPrivateMarker(visited);
							//add vertical edges along pillars
							ElementArray<Edge> & p0_edges = pillar_edges[(i+  0)*(dims[1]+1)+j+0];
							ElementArray<Edge> & p1_edges = pillar_edges[(i+1-q)*(dims[1]+1)+j+q];
							edges.reserve(edges.size()+p0_edges.size()+p1_edges.size());
							edges.insert(edges.end(),p0_edges.begin(),p0_edges.end());
							edges.insert(edges.end(),p1_edges.begin(),p1_edges.end());
							assert(count_duplicates(edges) == 0);
							//put block numbers to all nodes involved in pillars, so that we can figure out block numbers for constructed faces
							edges.SetPrivateMarker(mrk);
							//report block numbers on edges
							if(print_bn)
							{
								for(int k = 0; k < (int)edges.size(); ++k)
								{
									Storage::integer_array b = edges[k]->IntegerArray(block_number);
									Storage::integer_array e = edges[k]->IntegerArray(edge_number);
									std::cout << "edge " << k << " " << edges[k]->GetHandle() << " " << edges[k]->getBeg()->GetHandle() << " <-> " << edges[k]->getEnd()->GetHandle() << " blocks [" << b.size() << "]: ";
									for(int l = 0; l < (int)b.size()-1; ++l) std::cout << b[l] << "(" << b[l]%dims[0] << "," << b[l]/dims[0]%dims[1] << "," << b[l]/dims[0]/dims[1] << "):" << e[l] << " ";
									std::cout << b.back() << "(" << b.back()%dims[0] << "," << b.back()/dims[0]%dims[1] << "," << b.back()/dims[0]/dims[1] << "):" << e.back() << std::endl;
								}
							}
							if( print_bn ) std::cout << "intersection nodes: " << std::endl;
							for(int k = 0; k < (int)intersections.size(); ++k)
							{
								intersections[k]->SetPrivateMarker(internodes);
								ElementArray<Element> nedges = intersections[k]->getAdjElements(EDGE,mrk);
								block_number_union(intersections[k].getAsElement(),nedges,block_number,block_number);
								block_number_union(intersections[k].getAsElement(),nedges,block_number,node_blocks);
								intersections[k].IntegerArray(node_number).resize(intersections[k].IntegerArray(node_blocks).size(),-1);
								if( print_bn )
								{
									Storage::integer_array b = intersections[k]->IntegerArray(block_number);
									std::cout << "node " << k << " " << intersections[k]->GetHandle() << " blocks [" << b.size() << "]: ";
									for(int l = 0; l < (int)b.size()-1; ++l) std::cout << b[l] << "(" << b[l]%dims[0] << "," << b[l]/dims[0]%dims[1] << "," << b[l]/dims[0]/dims[1] << ") ";
									std::cout << b.back() << "(" << b.back()%dims[0] << "," << b.back()/dims[0]%dims[1] << "," << b.back()/dims[0]/dims[1] << ")" << std::endl;
								}
							}
							//compute block number for each edge from node's block numbers
							for(int k = 0; k < (int)edges.size(); ++k)
								block_number_intersection(edges[k],edges[k]->getAdjElements(NODE),block_number,block_intersection_number);
							//make segments curvy when necessery
							if( triangulated_edges )
							{
								//re-project intersection nodes onto triangulated face
								for(int k = 0; k < (int)intersections.size(); ++k)
								{
									Storage::real vt[3];
									Point p = make_point(intersections[k],pnt);
									unp.ActTri(p,vt);
									//change point
									intersections[k]->Coords()[0] = vt[0];
									intersections[k]->Coords()[1] = vt[1];
									intersections[k]->Coords()[2] = vt[2];
								}
								split_segments_tri(this,edges,intersections,transfer2,pnt,unp,visited,print_inter);
							}
							else if( curvilinear_edges )
								split_segments(this,edges,intersections,transfer2,pnt,unp,visited,print_inter);
							//distribute all the edges among blocks
							//add intersected edges into blocks, so that we can reconstruct top and bottom faces of each block
							for(int k = 0; k < (int)edges.size(); ++k) if( edges[k]->GetPrivateMarker(visited) ) //do not visit edges along pillars
							{
								Storage::integer_array b = edges[k].IntegerArray(block_number);
								for(int r = 0; r < (int)b.size(); ++r)
								{
//#if defined(USE_OMP)
//										omp_set_lock(&block_locks[b[r]]);
//#endif
									block_edges[b[r]].push_back(edges[k]);
//#if defined(USE_OMP)
//										omp_unset_lock(&block_locks[b[r]]);
//#endif
								}
							}
							edges.RemPrivateMarker(visited);
							//distribute edges to front and back blocks, so that we can assemble faces
							for(int k = 0; k < (int)edges.size(); ++k)
							{
								//intersect block numbers on nodes to detect to which blocks this edge belongs
								//ElementArray<Element> nodes = edges[k]->getAdjElements(NODE);
								//block_number_intersection(nodes,block_number,block_number_inter);
								//block_number_intersection(nodes,node_blocks,block_number_inter);
								Storage::integer_array block_number_inter = edges[k]->IntegerArray(block_intersection_number);
								if( print_bn || block_number_inter.empty())
								{
									std::cout << "edge " << k << " " << edges[k]->GetHandle() << " " << edges[k]->getBeg()->GetHandle() << " <-> " << edges[k]->getEnd()->GetHandle() << std::endl;
									{
										std::cout << " blocks [" << block_number_inter.size() << "]: ";
										for(int m = 0; m < (int)block_number_inter.size(); ++m)
										{
											int bi = block_number_inter[m] % dims[0];
											int bj = block_number_inter[m] / dims[0] % dims[1];
											int bk = block_number_inter[m] / dims[0] / dims[1];
											std::cout << block_number_inter[m] << "(" << bi << "," << bj << "," << bk << ") ";
										}
										std::cout << std::endl;
									}
									{
										Storage::integer_array bn = edges[k]->IntegerArray(block_number);
										Storage::integer_array nn = edges[k]->IntegerArray(edge_number);
										std::cout << " created by blocks [" << bn.size() << "]: ";
										for(int m = 0; m < (int)bn.size(); ++m)
										{
											int bi = bn[m] % dims[0];
											int bj = bn[m] / dims[0] % dims[1];
											int bk = bn[m] / dims[0] / dims[1];
											std::cout << bn[m] << "(" << bi << "," << bj << "," << bk << ")," << nn[m] << " ";
										}
										std::cout << std::endl;
									}
									{
										Storage::integer_array bn = edges[k]->getBeg()->IntegerArray(block_number);
										std::cout << " node " << edges[k]->getBeg()->GetHandle() << " blocks [" << bn.size() << "]: ";
										for(int m = 0; m < (int)bn.size(); ++m)
										{
											int bi = bn[m] % dims[0];
											int bj = bn[m] / dims[0] % dims[1];
											int bk = bn[m] / dims[0] / dims[1];
											std::cout << bn[m] << "(" << bi << "," << bj << "," << bk << ") ";
										}
										std::cout << std::endl;
									}
									{
										Storage::integer_array bn = edges[k]->getEnd()->IntegerArray(block_number);
										std::cout << " node " << edges[k]->getEnd()->GetHandle() << " blocks [" << bn.size() << "]: ";
										for(int m = 0; m < (int)bn.size(); ++m)
										{
											int bi = bn[m] % dims[0];
											int bj = bn[m] / dims[0] % dims[1];
											int bk = bn[m] / dims[0] / dims[1];
											std::cout << bn[m] << "(" << bi << "," << bj << "," << bk << ") ";
										}
										std::cout << std::endl;
									}
								}
								for(int m = 0; m < (int)block_number_inter.size(); ++m)
								{
									int bi = block_number_inter[m] % dims[0];
									int bj = block_number_inter[m] / dims[0] % dims[1];
									int bk = block_number_inter[m] / dims[0] / dims[1];
									assert(bi >= 0 && bi < dims[0]);
									assert(bj >= 0 && bj < dims[1]);
									if( bi == i-q && bj == j-(1-q) ) //back blocks
									{
										pillar_block_edges_back[bk].push_back(edges[k]);
										if( print_bn ) std::cout << "add edge " << k << " " << edges[k]->GetHandle() << " to back  block " << block_number_inter[m] << " (" << bi << "," << bj << "," << bk << ")" << std::endl;
									}
									else if(bi == i && bj == j) //front blocks
									{
										pillar_block_edges_front[bk].push_back(edges[k]);
										if( print_bn ) std::cout << "add edge " << k << " " << edges[k]->GetHandle() << " to front block " << block_number_inter[m] << " (" << bi << "," << bj << "," << bk << ")" << std::endl;
									}
									//skip non-current blocks
								}
							}
							std::vector< ElementArray<Edge> > * pillar_block_edges[2] = {&pillar_block_edges_back,&pillar_block_edges_front};
							int blocki[2] = {i-q,i};
							int blockj[2] = {j-(1-q),j};
							if( print_info )
							{
								std::cout << "back  block " << blocki[0] << "," << blockj[0] << std::endl;
								std::cout << "front block " << blocki[1] << "," << blockj[1] << std::endl;
								std::cout << "dims " << dims[0] << "," << dims[1] << std::endl;
							}
							//construct interfaces (some computations will be redundant)
							for(int m = 0; m < 2; ++m) //back and front
							{
								if( blocki[m] == -1 || blockj[m] == -1 ) {if(print_info) std::cout << "skip " << (m?"front ":"back ") << std::endl; continue;}
								if( blocki[m] == dims[0] || blockj[m] == dims[1] ) {if(print_info) std::cout << "skip " << (m?"front ":"back ") << std::endl; continue;}
								if( print_bedges ) std::cout << (m?"front ":"back ") << " column of blocks: " << blocki[m] << "," << blockj[m] << std::endl;
								for(int k = 0; k < dims[2]; ++k)  if( actnum.empty() || actnum[ECL_IJK_DATA(blocki[m],blockj[m],k)] ) //go down the piller
								{
									//retrive edge for the side
									ElementArray<Edge> & bedges = pillar_block_edges[m]->at(k);
									if( bedges.empty() ) continue;
									//remove any duplicates
									make_unique(bedges);
									int num_outer = 0;
									std::set<int> outer_edge_number;
									for(int l = 0; l < bedges.size(); ++l) //loop through edges of the block
									{
										//retrive block numbers of edges
										Storage::integer_array bn = bedges[l]->IntegerArray(block_number);
										Storage::integer_array en = bedges[l]->IntegerArray(edge_number);
										for(int r = 0; r < (int)bn.size(); ++r)
										{
											if( bn[r] == ECL_IJK_DATA(blocki[m],blockj[m],k) ) //this edge originally created by the block
											{
												bedges[l]->SetPrivateMarker(outer); //mark edge
												num_outer++;
												outer_edge_number.insert(en[r]);
											}
										}
									}
										
									//here we want to test that outer edges form a closed non-degenerate loop
									if( outer_edge_number.size() > 2 )//&& num_outer > 2 )
									{
										//move marked edges to the end
										std::sort(bedges.begin(),bedges.end(),Mesh::PrivateMarkerComparator(this,outer));
										//assert(std::is_sorted(bedges.begin(),bedges.end(),Mesh::PrivateMarkerComparator(this,outer)));
										if( print_bedges )
										{
											std::cout << (m?"front ":"back ") << "depth " << k << " block " << ECL_IJK_DATA(blocki[m],blockj[m],k) << " edges [" << bedges.size() << "]:" << std::endl;
											for(int l = 0; l < bedges.size(); ++l)
											{
												Storage::integer_array bn = bedges[l]->IntegerArray(block_number);
												Storage::integer_array en = bedges[l]->IntegerArray(edge_number);
												std::cout << "edge " << l << " " << bedges[l]->GetHandle() << " " << (bedges[l]->GetPrivateMarker(outer) ? "outer":"inner");
												std::cout << " blocks ";
												for(int r = 0; r < (int)bn.size(); ++r)
													std::cout << bn[r] << ":" << en[r] << " ";
												std::cout << std::endl;
											}
											assert(count_duplicates(bedges) == 0);
										}
										//there is enough edges to handle triangle
										if( bedges.size() > 2 )
										{
											//form faces out of edges
											incident_matrix<Edge> matrix(this,bedges.data(),bedges.data()+bedges.size(),bedges.size()-num_outer);
											//incident_matrix<Edge> matrix(this,bedges.data(),bedges.data()+bedges.size(),bedges.size()-num_outer, TagCoords(pnt));
											//collect all faces
											ElementArray<Edge> loop(this);
											while(matrix.find_shortest_loop(loop))
											{
												if( loop.size() > 2 ) //at least triangle
												{
													if( print_bedges )
													{
														std::cout << "Found loop of " << loop.size() << " edges:" <<std::endl;
														for(int g = 0; g < loop.size(); ++g)
														{
															std::cout << "edge " << g << " " << loop[g]->GetHandle() << " ";
															Storage::integer_array bn = loop[g]->IntegerArray(block_number);
															Storage::integer_array en = loop[g]->IntegerArray(edge_number);
															std::cout << " blocks ";
															for(int r = 0; r < (int)bn.size(); ++r)
																std::cout << bn[r] << ":" << en[r] << " ";
															std::cout << std::endl;
														}

													}
													//make face
													Face f = CreateFace(loop).first;
													f->IntegerArray(block_pair)[0] = i;
													f->IntegerArray(block_pair)[1] = j;
													f->IntegerArray(block_pair)[2] = q;
													if( TopologyErrorTag().isValid() && f->HaveData(TopologyErrorTag()) ) std::cout << __FILE__ <<":"<<__LINE__ << " topology error on face " << f->GetHandle() << std::endl;
													f->Integer(face_origin) = q;
													if(!f->CheckEdgeOrder()) std::cout << __FILE__ << ":" << __LINE__ << " bad edge order, edges " << loop.size() << std::endl;
													//detect block number
//#if defined(USE_OMP)
//														omp_set_lock(&block_locks[ECL_IJK_DATA(blocki[m],blockj[m],k)]);
//#endif
													block_faces[ECL_IJK_DATA(blocki[m],blockj[m],k)].push_back(f);
//#if defined(USE_OMP)
//														omp_unset_lock(&block_locks[ECL_IJK_DATA(blocki[m],blockj[m],k)]);
//#endif
													/*
													ElementArray<Element> nodes = f->getAdjElements(NODE);
													//block_number_intersection(nodes,block_number,block_number_inter);
													block_number_intersection(nodes,block_number,block_number_inter);
													if( block_number_inter.size() > 2 || block_number_inter.empty() )
													{
														std::cout << "face have " << block_number_inter.size() << " neighbours q is " << q << " block is " << ECL_IJK_DATA(blocki[m],blockj[m],k) << std::endl;
														print_bedges = true;
													}
														
													if( print_bedges )
													{
														std::cout << "nodes of the face " << f->GetHandle() << "(" << f->LocalID() << "): " << std::endl;
														for(int g = 0; g < nodes.size(); ++g)
														{
															std::cout << "node " << g << " " << nodes[g]->GetHandle() << " pillar " << nodes[g]->IntegerArray(pillar_num)[0] << "," << nodes[g]->IntegerArray(pillar_num)[1] << " " << (nodes[g]->GetPrivateMarker(internodes) ? "intersection" : "initial") << " pnt " << nodes[g]->RealArray(pnt)[0] << "," << nodes[g]->RealArray(pnt)[1] << " " << std::endl;
															Storage::integer_array bn = nodes[g]->IntegerArray(block_number);
															std::cout << " blocks [" << bn.size() << "]: ";
															for(int r = 0; r < (int)bn.size(); ++r)
																std::cout << bn[r] << "(" << bn[r]%dims[0] << "," << bn[r]/dims[0]%dims[1] << "," << bn[r]/dims[0]/dims[1] << ") ";
															std::cout << std::endl;
															Storage::integer_array cbn = nodes[g]->IntegerArray(node_blocks);
															Storage::integer_array cnn = nodes[g]->IntegerArray(node_number);
															std::cout << " created by [" << cbn.size() << "]: ";
															for(int r = 0; r < (int)cbn.size(); ++r)
																std::cout << cbn[r] << "(" << cbn[r]%dims[0] << "," << cbn[r]/dims[0]%dims[1] << "," << cbn[r]/dims[0]/dims[1] << ")," << cnn[r] << " ";
															std::cout << std::endl;
															ElementArray<Edge> edges = nodes[g]->getEdges(mrk);
															std::cout << " adjacent edges [" << edges.size() << "]: " << std::endl;
															for(int r = 0; r < (int)edges.size(); ++r)
															{
																Storage::integer_array ebn = edges[r]->IntegerArray(block_number);
																Storage::integer_array een = edges[r]->IntegerArray(edge_number);
																std::cout << "  " << edges[r]->GetHandle() << " " << edges[r]->getBeg()->GetHandle() << " <-> " << edges[r]->getEnd()-> GetHandle() << " blocks [" << ebn.size() << "]: ";
																for(int h = 0; h < (int)ebn.size(); ++h)
																	std::cout << ebn[h] << "(" << ebn[h]%dims[0] << "," << ebn[h]/dims[0]%dims[1] << "," << ebn[h]/dims[0]/dims[1] << ")," << een[h] << " ";
																std::cout << std::endl;
																	
															}

														}
														ElementArray<Edge> edges = f->getEdges();
														std::cout << "edges of the face: " << std::endl;
														for(int g = 0; g < edges.size(); ++g)
														{
															std::cout << "edge " << g << " " << edges[g]->GetHandle() << " " << edges[g]->getBeg()->GetHandle() << " <-> " << edges[g]->getEnd()->GetHandle() << " length " << edges[g]->Length() << " " << (edges[g]->GetPrivateMarker(outer) ? "outer":"inner") << " ";
															Storage::integer_array bn = edges[g]->IntegerArray(block_number);
															Storage::integer_array en = edges[g]->IntegerArray(edge_number);
																
															std::cout << " blocks ";
															for(int r = 0; r < (int)bn.size(); ++r)
																std::cout << bn[r] << "(" << bn[r]%dims[0] <<"," <<bn[r]/dims[0]%dims[1] << "," << bn[r]/dims[0]/dims[1] << ")," << en[r] << " ";
															std::cout << std::endl;

														}
														block_number_intersection(nodes,node_blocks,node_number_inter);
														std::cout << "node block number intersection: ";
														for(int r = 0; r < (int)node_number_inter.size(); ++r) std::cout << node_number_inter[r] << " (" << node_number_inter[r]%dims[0] << "," << node_number_inter[r]/dims[0] % dims[1] << "," << node_number_inter[r]/dims[0]/dims[1] <<") ";
														std::cout << std::endl;

														std::cout << "intersection: ";
													}
													for(int r = 0; r < (int)block_number_inter.size(); ++r)
													{
														if( print_bedges ) std::cout << block_number_inter[r] << " (" << block_number_inter[r]%dims[0] << "," << block_number_inter[r]/dims[0] % dims[1] << "," << block_number_inter[r]/dims[0]/dims[1] <<") ";
														block_faces[block_number_inter[r]].push_back(f);
													}
													if( print_bedges ) std::cout << std::endl;
													//if( block_number_inter.size() > 2 )
													print_bedges = false;
													*/
												}
												else if( !loop.empty() ) std::cout << __FILE__ << ":" << __LINE__ << " skip degenerate face " << loop.size() << std::endl;
											}
											if( !matrix.all_visited() )
											{
												std::cout << "Not all edges were visited, matrix: " << std::endl;
												matrix.print_matrix();
												std::cout << "Current block " << ECL_IJK_DATA(blocki[m],blockj[m],k) << std::endl;
												std::cout << "pillars: " << i << "," << j << "," << q << std::endl;
												for(int l = 0; l < bedges.size(); ++l)
												{
													Storage::integer_array bn = bedges[l]->IntegerArray(block_number);
													Storage::integer_array en = bedges[l]->IntegerArray(edge_number);
													std::cout << "edge " << l << " " << bedges[l]->GetHandle() << " " << bedges[l]->getBeg()->GetHandle() << " <-> " << bedges[l]->getEnd()->GetHandle() << " " << (bedges[l]->GetPrivateMarker(outer) ? "outer":"inner");
													std::cout << " blocks ";
													for(int r = 0; r < (int)bn.size(); ++r)
														std::cout << bn[r] << ":" << en[r] << " ";
													std::cout << std::endl;
												}
												std::cout << "edges projected: " << std::endl;
												for(int l = 0; l < bedges.size(); ++l)
												{
													std::cout << "(" << bedges[l]->getBeg()->RealArray(pnt)[0] << "," << bedges[l]->getBeg()->RealArray(pnt)[1] << ",0) <-> (" << bedges[l]->getEnd()->RealArray(pnt)[0] << "," << bedges[l]->getEnd()->RealArray(pnt)[1] << ",0)" << std::endl;
												}

												std::cout << "Outer edges [" << outer_edge_number.size() << "]:" << std::endl;
												for(std::set<int>::iterator it = outer_edge_number.begin(); it != outer_edge_number.end(); ++it)
													std::cout << *it << " ";
												std::cout << std::endl;
												std::cout << " num outer " << num_outer << std::endl;

												//form faces out of edgesm
												incident_matrix<Edge> matrix(this,bedges.data(),bedges.data()+bedges.size(),bedges.size()-num_outer,GridCoords(),true);
												//incident_matrix<Edge> matrix(this,bedges.data(),bedges.data()+bedges.size(),bedges.size()-num_outer,TagCoords(pnt),true);
												//collect all faces
												ElementArray<Edge> loop2(this);
												while(matrix.find_shortest_loop(loop2));
													
												std::cout << "block i " << i << " j " << j << " q " << q << std::endl;
												std::cout << "All edges: " << std::endl;
												for(int k = 0; k < (int)edges.size(); ++k)
												{
													std::cout << "(" << edges[k]->getBeg()->Coords()[0] << "," << edges[k]->getBeg()->Coords()[1] << "," << edges[k]->getBeg()->Coords()[2] << ") <-> (" << edges[k]->getEnd()->Coords()[0] << "," << edges[k]->getEnd()->Coords()[1] << "," << edges[k]->getEnd()->Coords()[2] << ")" << std::endl;
												}
											}
										}
									}
									//remove marker
									bedges.RemPrivateMarker(outer);
									//cleanup structure for reuse
									bedges.clear();
										
									//print_bn = false;
								}
							}
							//unmark edges
							edges.RemPrivateMarker(mrk);
							//clean-up structures
							intersections.RemPrivateMarker(internodes);
							edges.clear();
							intersections.clear();


							print_info=false;
							print_inter=false;
						} //j
					} //i
					DeleteTag(pnt);
					ReleasePrivateMarker(visited);
					ReleasePrivateMarker(mrk);
					ReleasePrivateMarker(outer);
					ReleasePrivateMarker(internodes);
				} //omp parallel
				//printf("\n");
				printf("finished creating faces for pairs of pillar along %s\n",q ? "ny":"nx");
			} //q
//#if defined(USE_OMP)
//				for(int q = 0; q < dims[0]*dims[1]*dims[2]; ++q) omp_destroy_lock(&block_locks[q]);
//#endif
			//do not need this tag on nodes
			//DeleteTag(block_number, NODE);
			//now construct top and bottom interfaces
			printf("started tops/bottoms/cells\n");
			Tag split_face;
			if( perform_splitting )
				split_face = CreateTag("SPLIT_FACE",DATA_REFERENCE,FACE,NONE,1); //points to edge that should be used to split face
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				ElementArray<Node> edge_nodes(this,2);
				for(int uneven = 0; uneven < 2; ++uneven) //uneven cell construction for OpenMP
				{
#if defined(USE_OMP)
#pragma omp for
#endif
					for(int i = uneven+beg_dims[0]; i < end_dims[0]; i+=2)
					{
						//printf("top/bottom/cells %6.2f%%\r", ((Storage::real)i)/((Storage::real)dims[0]-1)*100);
						//fflush(stdout);
						for(int j = beg_dims[1]; j < end_dims[1]; ++j)
						{
							for(int k = 0; k < dims[2]; ++k) if( actnum.empty() || actnum[ECL_IJK_DATA(i,j,k)] )
							{
								//current block number
								int cur = ECL_IJK_DATA(i,j,k);
								//bottom face -> 0,4,1,5
								//top face -> 2,7,3,6
								const int map[8][2] =
								{
									{0,0}, //edge 0 -> first position, bottom face
									{2,0}, //edge 1 -> third position, bottom face
									{0,1}, //edge 2 -> first position, top face
									{2,1}, //edge 3 -> third position, top face
									{1,0}, //edge 4 -> second position, bottom face
									{3,0}, //edge 5 -> fourth position, bottom face
									{3,1}, //edge 6 -> fourth position, top face
									{1,1} //edge 7 -> second position, top face
								};
								//which edges should be considered in reverse order
								const bool rev[2][4] =
								{
									{true,false,false,true},
									{false,false,true,true}
								};
								//array of ordered edges, first dimension - top or bottom, second dimension - edge number
								std::vector<HandleType> edges[2][4];
								//retrive set of edges of the block
								ElementArray<Edge> & cbe = block_edges[cur];
								for(int q = 0; q < cbe.size(); ++q)
								{
									Storage::integer_array bn = cbe[q].IntegerArray(block_number);
									Storage::integer_array en = cbe[q].IntegerArray(edge_number);
									assert(bn.size() == en.size());
									for(int r = 0; r < (int)bn.size(); ++r) if( bn[r] == cur ) //how the edge is represented on current block
										edges[map[en[r]][1]][map[en[r]][0]].push_back(cbe[q].GetHandle());
								}
								ElementArray<Edge> face_edges(this);
								//collect edges and create face
								for(int q = 0; q < 2; ++q) //bottom and top
								{
									for(int r = 0; r < 4; ++r) //all sides
									{
										if( rev[q][r] )
											face_edges.insert(face_edges.end(),edges[q][r].rbegin(),edges[q][r].rend());
										else
											face_edges.insert(face_edges.end(),edges[q][r].begin(),edges[q][r].end());
									}
									make_unique(face_edges); //somehow duplicate appears there
									if( face_edges.size() > 2 )
									{
										Face f = CreateFace(face_edges).first;
										if( TopologyErrorTag().isValid() && f->HaveData(TopologyErrorTag()) ) std::cout << __FILE__ <<":"<<__LINE__ << " topology error on face " << f->GetHandle() << std::endl;
											
										if( perform_splitting )
										{
											int split = 0;
											bool senw = block_nodes[cur*8+0] == block_nodes[cur*8+4] && block_nodes[cur*8+3] == block_nodes[cur*8+7];
											bool swne = block_nodes[cur*8+1] == block_nodes[cur*8+5] &&	block_nodes[cur*8+2] == block_nodes[cur*8+6];
											if( senw && swne )
												split = 0;
											else if( senw )
												split = 1; //SE-NW diagonal
											else if( swne )
												split = 2; //SW-NE diagonal
												
												
											HandleType diag = f->Reference(split_face);
											//split up quad faces into triagnles if they degenerate
											if( split )
											{
												if( diag == InvalidHandle() )
												{
													if( split == 1 ) //SE-NW diagonal
													{
														edge_nodes.at(0) = block_nodes[cur*8+1+q*4];
														edge_nodes.at(1) = block_nodes[cur*8+2+q*4];
														f->Reference(split_face) = CreateEdge(edge_nodes).first.GetHandle();
													}
													else if( split == 2 ) //SW-NE diagonal
													{
														edge_nodes.at(0) = block_nodes[cur*8+0+q*4];
														edge_nodes.at(1) = block_nodes[cur*8+3+q*4];
														f->Reference(split_face) = CreateEdge(edge_nodes).first.GetHandle();
													}
												}
												else //the face was already split
												{
													int was_split = 0;
													Edge e(this,diag);
													if( e->getBeg()->GetHandle() == block_nodes[cur*8+1+q*4] &&
														e->getEnd()->GetHandle() == block_nodes[cur*8+2+q*4])
														was_split = 1; //was SE-NW
													else if( e->getBeg()->GetHandle() == block_nodes[cur*8+0+q*4] &&
																e->getEnd()->GetHandle() == block_nodes[cur*8+3+q*4])
														was_split = 2; //was SW-NE
														
													if( !was_split )
													{
														std::cout << "Cannot detect how the face was priviously split" << std::endl;
														std::cout << "edge:  " << e->getBeg()->GetHandle() << " <-> " << e->getEnd()->GetHandle() << std::endl;
														std::cout << "SE-NW: " << block_nodes[cur*8+0+q*4] << " <-> " << block_nodes[cur*8+3+q*4] << std::endl;
														std::cout << "SW-NE: " << block_nodes[cur*8+1+q*4] << " <-> " << block_nodes[cur*8+2+q*4] << std::endl;
													}
														
													if( split != was_split )
													{
														//replace with node
														e->Delete();
														Storage::real xyz[3] = {0,0,0};
														for(int l = 0; l < 4; ++l)
														{
															Storage::real_array c = Node(this,block_nodes[cur*8+l+q*4]).Coords();
															xyz[0] += c[0]*0.25;
															xyz[1] += c[1]*0.25;
															xyz[2] += c[2]*0.25;
														}
														f->Reference(split_face) = CreateNode(xyz)->GetHandle();
														/*
														std::cout << "Got different split direction" << std::endl;
														std::cout << "was " << e->getBeg()->GetHandle() << " <-> " << e->getEnd()->GetHandle() << " " << (was_split == 1 ? "SE-NW" : "SW-NE") << std::endl;
														std::cout << (q? "top" : "bottom") << " nodes: " << std::endl;
														for(int l = 0; l < 4; ++l)
															std::cout << l+q*4 << ":" << block_nodes[cur*8+l+q*4] << " ";
														std::cout << std::endl;
														std::cout << (1-q? "top" : "bottom") << " nodes: " << std::endl;
														for(int l = 0; l < 4; ++l)
															std::cout << l+(1-q)*4 << ":" << block_nodes[cur*8+l+(1-q)*4] << " ";
														std::cout << std::endl;
														std::cout << "SE-NW: " << block_nodes[cur*8+0+q*4] << " <-> " << block_nodes[cur*8+3+q*4] << std::endl;
														std::cout << "SW-NE: " << block_nodes[cur*8+1+q*4] << " <-> " << block_nodes[cur*8+2+q*4] << std::endl;
														std::cout << "want direction " << (split == 1? "SE-NW" : "SW-NE") << std::endl;
															*/
													}
												}
											}
										}


										f->Integer(face_origin) = 2;
										if(!f->FixEdgeOrder())
											//if(!f->CheckEdgeOrder())
										{
											std::cout << __FILE__ << ":" << __LINE__ << " bad edge order, edges " << face_edges.size() << std::endl;
											std::cout << "block: " << cur << " (" << i << "," << j << "," << k << ") " << (q ? "top" : "bottom") << std::endl;
											for(int l = 0; l < face_edges.size(); ++l)
											{
												Storage::integer_array bn = face_edges[l]->IntegerArray(block_number);
												Storage::integer_array en = face_edges[l]->IntegerArray(edge_number);
												std::cout << "edge " << l << " " << face_edges[l]->GetHandle() << " ";
												std::cout << face_edges[l]->getBeg()->GetHandle() << "<->" << face_edges[l]->getEnd()->GetHandle() << " ";
												std::cout << "blocks: ";
												for(int p = 0; p < (int)bn.size(); ++p)
													std::cout << bn[p] << "(" << bn[p]%dims[0] << "," << bn[p]/dims[0]%dims[1] << "," << bn[p]/dims[0]/dims[1] << "):" << en[p] << " ";
												std::cout << std::endl;
											}
										}
										block_faces[cur].push_back(f);
									}
									else if( !face_edges.empty() ) std::cout << __FILE__ << ":" << __LINE__ << " skip degenerate face " << face_edges.size() << " " << (q?"top":"bottom") << " of block " << cur << " (" << i << "," << j << "," << k << ") actnum " << actnum[ECL_IJK_DATA(i,j,k)] << std::endl;
									face_edges.clear();
								}
								make_unique(block_faces[cur]); //some faces may be added twice?
								if( block_faces[cur].size() > 3 )
								{
									Cell c = CreateCell(block_faces[cur]).first;
									blocks[cur] = c->GetHandle();
									c->Integer(cell_number) = cur+1;
									c->IntegerArray(block_index)[0] = i;
									c->IntegerArray(block_index)[1] = j;
									c->IntegerArray(block_index)[2] = k;

									//c->IntegerArray(cell_number)[1] = i;
									//c->IntegerArray(cell_number)[2] = j;
									//c->IntegerArray(cell_number)[3] = k;
									if( TopologyErrorTag().isValid() )
									{
										if( c->Integer(TopologyErrorTag()) & PROHIBIT_MULTIPOLYGON )
										{
											std::cout << "Cell " << c->GetHandle() << " id " << c->LocalID() << " " << Element::GeometricTypeName(c->GetGeometricType()) << " block " << i << "," << j << "," << k << "("<<cur<<") got topology error" << std::endl;
											std::cout << "input faces [" << block_faces[cur].size() << "]: " << std::endl;
											std::map<HandleType,int>  visits;
											for(int r = 0; r < (int)block_faces[cur].size(); ++r)
											{
												ElementArray<Cell> fcells = block_faces[cur][r].getCells();
												ElementArray<Node> fnodes = block_faces[cur][r].getNodes();
												ElementArray<Edge> fedges = block_faces[cur][r].getEdges();
												std::cout << " face " << block_faces[cur][r].GetHandle() << " " << Element::GeometricTypeName(block_faces[cur][r]->GetGeometricType()) << " origin " << block_faces[cur][r]->Integer(face_origin) << std::endl;
												std::cout << "  cells [" << fcells.size() << "]: " << std::endl;
												for(int l = 0; l < (int)fcells.size(); ++l)
													std::cout << "  " << fcells[l]->GetHandle() << " " << Element::GeometricTypeName(fcells[l]->GetGeometricType()) << " block " << fcells[l]->Integer(cell_number) << std::endl;

												std::cout << "  edges [" << fedges.size() << "]: " << std::endl;
												for(int l = 0; l < (int)fedges.size(); ++l)
												{
													Storage::integer_array bn = fedges[l]->IntegerArray(block_number);
													Storage::integer_array en = fedges[l]->IntegerArray(edge_number);
													std::cout << "  " << fedges[l]->GetHandle() << " nnodes " << fedges[l]->nbAdjElements(NODE) << " " << fedges[l]->getBeg()->GetHandle() << " <-> " << fedges[l]->getEnd()->GetHandle() << " block_numbers [" << bn.size() << "]: ";
													for(int v = 0; v < (int)bn.size(); ++v) std::cout << bn[v] << ":" << en[v] << " ";
													std::cout << std::endl;
													visits[fedges[l]->GetHandle()]++;
												}
												std::cout << "  nodes [" << fnodes.size() << "]: ";
												for(int l = 0; l < (int)fnodes.size(); ++l)
													std::cout << fnodes[l]->GetHandle() << " ";
												std::cout << std::endl;
											}
											for(std::map<HandleType,int>::iterator kt = visits.begin(); kt != visits.end(); ++kt)
												std::cout << " edge " << kt->first << " visited " << kt->second << std::endl;
											std::cout << " total cell nodes: " << c->nbAdjElements(NODE) << std::endl;
											ElementArray<Node> cnodes = c->getNodes();
											for(int r = 0; r < (int)cnodes.size(); ++r)
												std::cout << "  " << cnodes[r]->GetHandle() << std::endl;
										}
										else if( c->Integer(TopologyErrorTag()) & TRIPLE_SHARED_FACE )
										{
											ElementArray<Face> faces = c->getFaces();
											for(int r = 0; r < (int)faces.size(); ++r)
												if( faces[r]->nbAdjElements(CELL) > 2 )
												{
													ElementArray<Cell> fcells = faces[r]->getCells();
													std::cout << "face " << faces[r]->GetHandle() << "(" << faces[r]->LocalID() << ") " << Element::GeometricTypeName(faces[r]->GetGeometricType()) << " area " << faces[r]->Area() << " cells[" << fcells.size() << "]: ";
													for(int l = 0; l < (int)fcells.size(); ++l)
														std::cout << fcells[l]->LocalID() << " ";
													std::cout << std::endl;
												}
										}
									}
								}
								else if( !block_faces[cur].empty() )
								{
									//depending on actnum, mark faces that are not connected
									//std::cout << __FILE__ << ":" << __LINE__ << " skip degenerate cell " << block_faces[cur].size() << " block " << cur << " (" << i << "," << j << "," << k << ") actnum " << actnum[ECL_IJK_DATA(i,j,k)] << std::endl;
								}
							} //k
						} //j
					} //i
				} //q
			}
			printf("finished tops/bottoms/cells\n");
			if( perform_splitting )
			{
				printf("started splitting faces of degenerate cells\n");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(integer it = 0; it < FaceLastLocalID(); ++it) if( isValidFace(it) )
				{
					Face f = FaceByLocalID(it);
					HandleType diag = f->Reference(split_face);
					if( diag != InvalidHandle() )
					{
						if( GetHandleElementType(diag) == NODE )
						{
							ElementArray<Edge> split_edges(this);
							ElementArray<Node> edge_nodes(this,2);
							ElementArray<Node> face_nodes = f->getNodes();
							edge_nodes.at(0) = diag;
							for(int k = 0; k < (int)face_nodes.size(); ++k)
							{
								edge_nodes[1] = face_nodes[k];
								split_edges.push_back(CreateEdge(edge_nodes).first);
							}
							Face::SplitFace(f,split_edges,0);
						}
						else
						{
							if( diag != InvalidHandle() )
								Face::SplitFace(f,ElementArray<Edge>(this,1,diag),0);
						}
					}
				}
				printf("finished splitting faces of degenerate cells\n");
			}
			//printf("\n");
			//cleanup data
			//DeleteTag(split_face);
			DeleteTag(edge_number);
			//DeleteTag(block_number);
			//crack up the mesh along degenerate active cells
			//populate properties to blocks
			Tag tagporo, tagsatnum, tagperm;
			if( !poro.empty() ) tagporo = CreateTag("PORO",DATA_REAL,CELL,NONE,1);
			if( !perm.empty() ) tagperm = CreateTag("PERM",DATA_REAL,CELL,NONE,project_perm?6:3);
			if( !satnum.empty() ) tagsatnum = CreateTag("SATNUM",DATA_INTEGER,CELL,NONE,1);
			const int face_node_nums[3][2][4]=
			{
				
				{{0,2,4,6},{1,3,5,7}}, //u-direction face nodes n,p
				{{0,1,4,5},{2,3,6,7}}, //v-direction face nodes n,p
				{{0,1,2,3},{4,5,6,7}}  //w-direction face nodes n,p
			};
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				rMatrix UVW(3,3), U(3,3), S(3,3), Sinv(3,3), V(3,3), K(3,3);
#if defined(USE_OMP)
#pragma omp for
#endif
				for(integer it = 0; it < CellLastLocalID(); ++it) if( isValidCell(it) )
				{
					Cell c = CellByLocalID(it);
					if( c->Integer(cell_number) > 0 ) //maybe this cell existed before
					{
						integer q = c->Integer(cell_number)-1;
						if( q < 0 ) continue;
						if( !poro.empty() ) c->Real(tagporo) = poro[q];
						if( !satnum.empty() ) c->Integer(tagsatnum) = satnum[q];
						if( !perm.empty() )
						{
							Storage::real_array arr_perm = c->RealArray(tagperm);
							if( project_perm )
							{
								ElementArray<Face> faces = c->getFaces();
								Storage::real u[3], L;
								//compute grid block coordinates
								for(int r = 0; r < 3; ++r)
								{
									for(int l = 0; l < 4; ++l)
									{
										u[0] -= Node(this,block_nodes[q*8+face_node_nums[r][0][l]])->Coords()[0];
										u[1] -= Node(this,block_nodes[q*8+face_node_nums[r][0][l]])->Coords()[1];
										u[2] -= Node(this,block_nodes[q*8+face_node_nums[r][0][l]])->Coords()[2];
										u[0] += Node(this,block_nodes[q*8+face_node_nums[r][1][l]])->Coords()[0];
										u[1] += Node(this,block_nodes[q*8+face_node_nums[r][1][l]])->Coords()[1];
										u[2] += Node(this,block_nodes[q*8+face_node_nums[r][1][l]])->Coords()[2];
									}
									L = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
									if( L ) u[0] /= L, u[1] /= L, u[2] /= L;
									//fill row
									UVW(r,0) = u[0];
									UVW(r,1) = u[1];
									UVW(r,2) = u[2];
								}
								UVW.SVD(U,S,V);
								for(int r = 0; r < 3; ++r)
								{
									if( S(r,r) > 0.0 )
										Sinv(r,r) = 1.0/S(r,r);
									else 
										Sinv(r,r) = 0.0;
								}
								K.Zero();
								K(0,0) = perm[3*q+0];
								K(1,1) = perm[3*q+1];
								K(2,2) = perm[3*q+2];
								
								K = V*Sinv*U.Transpose()*K*U*Sinv*V.Transpose();

								arr_perm[0] = K(0,0);
								arr_perm[1] = K(0,1);
								arr_perm[2] = K(0,2);
								arr_perm[3] = K(1,1);
								arr_perm[4] = K(1,2);
								arr_perm[5] = K(2,2);
							}
							else
							{
								arr_perm[0] = perm[q*3+0];
								arr_perm[1] = perm[q*3+1];
								arr_perm[2] = perm[q*3+2];
							}
						}
					}
				}
			}
			//compute cell centers that lay inside
			if( !curvilinear_edges )
			{
				GeomParam table;
				table[CENTROID] = CELL;
				PrepareGeometricData(table);
				//overwrite centroid info
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(integer it = 0; it < CellLastLocalID(); ++it) if( isValidCell(it) )
				{
					Cell c = CellByLocalID(it);
					integer bnum = c->Integer(cell_number)-1;
					if( bnum >= 0 ) //maybe this cell existed before
					{
						real ctop[3] = {0.0,0.0,0.0}, cbottom[3] = {0.0,0.0,0.0};
						for(int l = 0; l < 4; ++l)
						{
							real_array bc = Node(this,block_nodes[bnum*8+l+0]).Coords();
							real_array tc = Node(this,block_nodes[bnum*8+l+4]).Coords();
							cbottom[0] += bc[0]*0.25;
							cbottom[1] += bc[1]*0.25;
							cbottom[2] += bc[2]*0.25;
							ctop[0] += tc[0]*0.25;
							ctop[1] += tc[1]*0.25;
							ctop[2] += tc[2]*0.25;
						}
						real_array cnt = c->RealArray(centroid_tag);
						cnt[0] = (cbottom[0]+ctop[0])*0.5;
						cnt[1] = (cbottom[1]+ctop[1])*0.5;
						cnt[2] = (cbottom[2]+ctop[2])*0.5;
						
					}
				}
			}
			int num_outside = 0, num_total = 0;
			for(integer it = 0; it < CellLastLocalID(); ++it) if( isValidCell(it) )
			{
				Cell c = CellByLocalID(it);
				integer bnum = c->Integer(cell_number)-1;
				if( bnum >= 0 ) //maybe this cell existed before
				{
					Storage::real cnt[3];
					c->Centroid(cnt);
					if(!c->Inside(cnt) )
					{
						//std::cout << "Centroid is outside of cell " << bnum << std::endl;
						num_outside++;
					}
					num_total++;
				}
			}

				
			std::cout << "Centroid is outside for " << num_outside << " cells out of " << num_total << " cells " << std::endl;
				
			if( false )
			for(integer it = 0; it < NodeLastLocalID(); ++it) if( isValidNode(it) )
			{
				Node n = NodeByLocalID(it);
				if( n->nbAdjElements(CELL) == 0 )
				{
					//bool mrk[MarkerFields], mrkprv[MarkerFields];
					Element::adj_type & hc = HighConn(n->GetHandle());
					Element::adj_type & lc = LowConn(n->GetHandle());
					Storage::bulk_array mrk = n->BulkArray(tag_markers);
					Storage::bulk_array private_mrk = n->BulkArray(tag_private_markers[0]);
					//n->GetMarkerSpace(mrk);
					//ElementArray<Edge> nedges = n->getEdges();
					ElementArray<Face> nfaces = n->getFaces();
					Storage::integer_array bn =  n->IntegerArray(node_blocks);
					std::cout << "orphan node " << n->GetHandle() << " lc " << lc.size() << " hc " << hc.size() << " blocks numbers [" << bn.size() << "]: ";
					for(int k = 0; k < (int)bn.size(); ++k) std::cout << bn[k] << " ";
					std::cout << std::endl;
					std::cout << " adjacent faces [" << nfaces.size() << "]: " << std::endl;
					for(int k = 0; k < nfaces.size(); ++k)
					{
						std::cout << "  " << nfaces[k]->GetHandle() << " " << Element::GeometricTypeName(nfaces[k]->GetGeometricType()) << std::endl;
						ElementArray<Node> fnodes = nfaces[k]->getNodes();
						std::cout << "  nodes [" << fnodes.size() << "]: ";
						for(int q = 0; q < fnodes.size(); ++q)
						{
							std::cout << fnodes[q]->GetHandle() << " ";
						}
						std::cout << std::endl;
						std::cout << "  cells [" << nfaces[k]->nbAdjElements(CELL) << "]: " << std::endl;
						ElementArray<Cell> fcells = nfaces[k]->getCells();
						for(int q = 0; q < fcells.size(); ++q)
						{
							ElementArray<Node> cnodes = fcells[q]->getNodes();
							std::cout << "  " <<fcells[q]->GetHandle() << " " << Element::GeometricTypeName(fcells[q]->GetGeometricType()) << " (";
							for(int r = 0; r < cnodes.size()-1; ++r) std::cout << cnodes[r]->GetHandle() << ",";
							std::cout << cnodes.back()->GetHandle() << ") " << std::endl;
							ElementArray<Face> cfaces = fcells[q]->getFaces();
							std::cout << "   cell faces [" << cfaces.size() << "]:" << std::endl;
							for(int r = 0; r < cfaces.size(); ++r)
							{
								ElementArray<Node> ffnodes = cfaces[r]->getNodes();
								std::cout << "   " << Element::GeometricTypeName(cfaces[r]->GetGeometricType()) << " (";
								for(int l = 0; l < ffnodes.size()-1; ++l) std::cout << ffnodes[l]->GetHandle() << ",";
								std::cout << ffnodes.back()->GetHandle() << ") " << std::endl;
							}
								
						}
						std::cout << std::endl;
					}
					std::cout << std::endl;
					std::cout << "markers: ";
					for(int k = 0; k < (int)mrk.size(); ++k) std::cout << (int)mrk[k] << " ";
					std::cout << std::endl;
					std::cout << "private markers: ";
					for(int k = 0; k < (int)private_mrk.size(); ++k) std::cout << (int)private_mrk[k] << " ";
					std::cout << std::endl;
						
						
				}//if
			}//it
		} //if gtype
		
		//std::cout << "Total wells: " << wells.size() << std::endl;

		Tag tagnwell = CreateTag("NWELL",DATA_INTEGER,CELL,CELL,1);
		std::vector<HandleType> well_sets;
		int nwell = 0;
		for(compdat_wells::iterator it = wells.begin(); it != wells.end(); ++it)
		{
			//std::cout << "Well " << nwell << " name " << it->first << " elements ";
			ElementSet set = CreateSet(it->first).first;
			for(compdat_entries::iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
			{
				if( jt->i >= beg_dims[0] && jt->i < end_dims[0] && jt->j >= beg_dims[1] && jt->j < end_dims[1] )
				{
					for(int k = jt->k1; k <= jt->k2; ++k) if( k >= 0 && k < dims[2] )
					{
						HandleType h = blocks[ECL_IJK_DATA(jt->i,jt->j,k)];
						if( h != InvalidHandle() )
						{
							set.PutElement(h);
							Integer(h,tagnwell) = nwell;
						}
					}
				}
			}
			//std::cout << set.Size() << std::endl;
			for(ElementSet::iterator it = set.Begin(); it != set.End(); ++it)
			{
				real c[3];
				it->Centroid(c);
				//std::cout << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << " " << c[0] << " " << c[1] << " " << c[2] << std::endl;
			}
			well_sets.push_back(set->GetHandle());
			nwell++;
		}

		
		
		/*
		{
			std::fstream fout("wells.txt",std::ios::out);
			fout << nwell << std::endl;
			for(int i = 0; i < nwell; ++i)
			{
				ElementSet set(this,well_sets[i]);
				fout << set->GetName() << " 0 bhp 100 injector peacman_vertical" << std::endl;
			}
			fout << nwell*2 << std::endl;
			for(int i = 0; i < nwell; ++i)
			{
				ElementSet set(this,well_sets[i]);
				HandleType * handles = set.getHandles();
				int nhandles = set.nbHandles();
				Cell beg(this,handles[0]), end(this,handles[nhandles-1]);
				real cbeg[3], cend[3];
				beg.Centroid(cbeg);
				end.Centroid(cend);
				fout << cbeg[0] << " " << cbeg[1] << " " << cbeg[2] << std::endl;
				fout << cend[0] << " " << cend[1] << " " << cend[2] << std::endl;
			}
			for(int i = 0; i < nwell; ++i)
			{
				fout << 1 << std::endl;
				fout << "1 " << i*2+0 << " " << i*2+1 << " 0.0001 0 0" << std::endl;
			}
			fout.close();
		}
		*/

		//if( GetProcessorsNumber() )
		ResolveShared();
	} //LoadECL
} //namespace

#endif
