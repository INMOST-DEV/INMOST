#include "inmost.h"
#include <stack>
#if defined(USE_MESH)

#define MEASURE_EPS 1.0e-11
//#define NEW_ALGORITHM
//#define RECORD_PATH
//#define DEEP_RECURSE
//TODO:
// incident_matrix class should measure for minimal volume,
// possibly check and update from projects/OctreeCutcell/octgrid.cpp

namespace INMOST
{

	__INLINE void make_vec(const Storage::real p1[3], const Storage::real p2[3], Storage::real out[3])
	{
		out[0] = p1[0] - p2[0];
		out[1] = p1[1] - p2[1];
		out[2] = p1[2] - p2[2];
	}

	__INLINE void cross_prod(const Storage::real v1[3],const Storage::real v2[3], Storage::real out[3])
	{
		out[0] = v1[1]*v2[2] - v1[2]*v2[1];
		out[1] = v1[2]*v2[0] - v1[0]*v2[2];
		out[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}

	__INLINE Storage::real __det3d(Storage::real a, Storage::real b, Storage::real c,
	                             Storage::real d, Storage::real e, Storage::real f,
	                             Storage::real g, Storage::real h, Storage::real i ) 
	{
		return a*e*i - c*e*g + b*f*g - a*f*h + c*d*h - b*d*i;
	}
	
	__INLINE Storage::real __det3v(const Storage::real * x,const Storage::real * y,const Storage::real * z) 
	{
		return __det3d(x[0], x[1], x[2],  y[0], y[1], y[2],  z[0], z[1], z[2]);
	}

	__INLINE Storage::real dot_prod(const Storage::real v1[3],const Storage::real v2[3])
	{
		return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
	}
	
	typedef struct orient_face_t
	{
		Edge bridge;
		Node first;
		Face face;
		orient_face_t(Edge _bridge, Node _first, Face _face)
		:bridge(_bridge),first(_first),face(_face)
		{
		}
	} orient_face;
	
	class AbstractCoords
	{
	public:
		virtual Storage::real_array Get(const Node & s) const = 0;
		virtual AbstractCoords * Copy() const = 0;
		virtual ~AbstractCoords() {}
	};
	
	class GridCoords : public AbstractCoords
	{
	public:
		Storage::real_array Get(const Node & s) const
		{
			return s.Coords();
		}
		AbstractCoords * Copy() const {return new GridCoords(*this);}
	};
	
	class TagCoords : public AbstractCoords
	{
		Tag t;
	public:
		TagCoords(Tag _t) : t(_t) {}
		TagCoords & operator =(TagCoords const & b) {t = b.t; return *this;}
		TagCoords(const TagCoords & b) : t(b.t) {}
		Storage::real_array Get(const Node & s) const
		{
			return s.RealArray(t);
		}
		AbstractCoords * Copy() const {return new TagCoords(*this);}
	};
	
#if !defined(NEW_ALGORITHM)
	template<class T>
	class incident_matrix
	{
		Mesh * mesh;
		dynarray< unsigned char, 4096 > matrix;
		dynarray< char ,256 > visits;
		dynarray< char ,256 > visits0;
		dynarray<HandleType, 256> head_column;
		dynarray<HandleType, 256> head_row;
		dynarray<unsigned char ,256> head_row_count;
		dynarray<unsigned, 256> insert_order;
		bool exit_recurse;
		ElementArray<T> min_loop, temp_loop; //used as return
		dynarray< char , 256 > hide_column;
		dynarray< char , 256 > hide_row;
		dynarray< char , 256 > stub_row;
#if defined(RECORD_PATH)
		std::vector< std::pair<std::vector<int>,double> > remember;
#endif
		double min_loop_measure;
		bool print;
		AbstractCoords * coords;
		
		bool do_hide_row(unsigned k)
		{
			if( hide_column[k] == 0 )
			{
				hide_column[k] = 1;
				for(unsigned i = 0; i < head_row_count.size(); i++)
					if( matrix[k*head_row_count.size()+i] == 1 )
					{
						head_row_count[i] -= 1;
						if( head_row_count[i] == 0 )
						{
							hide_row[i] = 1;
							stub_row[i] = 0;
						}
					}
				insert_order.pop_back();
			}
			return true;
		}
		
		bool do_show_row(unsigned k)
		{
			if( hide_column[k] == 1 )
			{
				hide_column[k] = 0;
				
				bool success = true;
				for(unsigned i = 0; i < head_row_count.size(); i++)
					if( matrix[k*head_row_count.size()+i] == 1 )
					{
						head_row_count[i] += 1;
						if( head_row_count[i] > 0 ) hide_row[i] = 0;
						if( head_row_count[i] > 2 ) success = false;
					}
				if( print )
				{
					std::cout << (success? "ok":"fail") << std::endl;
					std::cout << "row count:";
					for(unsigned i = 0; i < head_row_count.size(); i++)
						std::cout << " " << (int)head_row_count[i];
					std::cout << std::endl;
					//std::cout << " success? " << (success?"yes":"no") << std::endl;
				}
				insert_order.push_back(k);
				if( !success ) do_hide_row(k);
				return success;
				
			} else return true;
		}
		bool test_success()
		{
			bool success = true;
			for(unsigned j = 0; j < head_row_count.size(); j++)
			{
				if( head_row_count[j] == 1 )
				{
					success = false;
					break;
				}
			}
			return success;
		}
		Storage::real compute_measure(ElementArray<T> & data)
		{
			Storage::real measure = 0, tmp;
			for(typename ElementArray<T>::size_type k = 1; k < data.size(); ++k)
			{
				mesh->GetGeometricData(data[k].GetHandle(),MEASURE,&tmp);
				measure += tmp;
			}
			return measure;
		}
		Storage::real compute_measure_vol(ElementArray<T> & data)
		{
			Storage::real measure = 0;
			if( data[0]->GetElementDimension() == 1 ) //this is edge //use geometric dimension here for 2d compatibility
			{
				//calculate area
				//int mdim = data[0]->GetMeshLink()->GetDimensions();
				ElementArray<Node> nodes,n1,n2;
				n1 = data[0]->getNodes();
				n2 = data[1]->getNodes();
				if( n1[0] == n2[0] || n1[0] == n2[1])
				{
					nodes.push_back(n1[1]);
					nodes.push_back(n1[0]);
				}
				else
				{
					nodes.push_back(n1[0]);
					nodes.push_back(n1[1]);
				}
				for(typename ElementArray<T>::size_type j = 1; j < data.size(); j++)
				{
					n1 = data[j]->getNodes();
					assert((nodes.back() == n1[0] || nodes.back() == n1[1]));
					if( nodes.back() == n1[0] )
						nodes.push_back(n1[1]);
					else
						nodes.push_back(n1[0]);
				}
				
				Storage::real x[3] = {0,0,0};
				Storage::real_array x0 = coords->Get(nodes[0]);
				for(unsigned i = 1; i < nodes.size()-1; i++)
				{
					Storage::real_array v1 = coords->Get(nodes[i]);
					Storage::real_array v2 = coords->Get(nodes[i+1]);
					if( v1.size() == 3 && v2.size() == 3 )
					{
						x[0] += (v1[1]-x0[1])*(v2[2]-x0[2]) - (v1[2]-x0[2])*(v2[1]-x0[1]);
						x[1] += (v1[2]-x0[2])*(v2[0]-x0[0]) - (v1[0]-x0[0])*(v2[2]-x0[2]);
					}
					x[2] += (v1[0]-x0[0])*(v2[1]-x0[1]) - (v1[1]-x0[1])*(v2[0]-x0[0]);
				}
				measure = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])*0.5;
				
			}
			else //this is 3d face
			{
				//firstly, have to figure out orientation of each face
				//mark all faces, so that we can perform adjacency retrival
				MarkerType mrk = mesh->CreatePrivateMarker();
				MarkerType rev = mesh->CreatePrivateMarker(); //reverse orientation
				for(int k = 1; k < data.size(); ++k)
					data[k]->SetPrivateMarker(mrk); //0-th face orientation is default
				Node n1,n2; //to retrive edge
				bool reverse = false; //reverse orientation in considered face
				std::deque< orient_face > stack; //edge and first node and face for visiting
				//todo: can do faster by retriving edges and going over their nodes
				//should not use FindSharedAdjacency
				ElementArray<Edge> edges = data[0]->getEdges();
				do
				{
					//figure out starting node order
					if( edges[0]->getBeg() == edges[1]->getBeg() ||
					    edges[0]->getBeg() == edges[1]->getEnd() )
					{
						n1 = edges[0]->getEnd();
						n2 = edges[0]->getBeg();
					}
					else
					{
						n1 = edges[0]->getBeg();
						n2 = edges[0]->getEnd();
					}
					//schedule unvisited adjacent faces
					for(typename ElementArray<Edge>::size_type j = 0; j < edges.size(); j++)
					{
						//schedule face adjacent to considered edge
						ElementArray<Face> adjacent = edges[j]->getFaces(mrk);
						assert(adjacent.size() <= 1);
						if( !adjacent.empty() && adjacent[0].GetPrivateMarker(mrk))
						{
							adjacent.RemPrivateMarker(mrk);
							stack.push_back(orient_face(edges[j],reverse ? n2 : n1,adjacent[0]));
						}
						//update edge nodes
						n1 = n2; //current end is new begin
						//find new end
						if( n2 == edges[(j+1)%edges.size()]->getBeg() )
							n2 = edges[(j+1)%edges.size()]->getEnd();
						else
							n2 = edges[(j+1)%edges.size()]->getBeg();
					}
					if( stack.empty() ) break;
					//get entry from stack
					orient_face r = stack.front();
					//remove face from stack
					stack.pop_front();
					//retrive edges for new face
					edges = r.face->getEdges();
					reverse = false;
					//figure out starting node order
					if( edges[0]->getBeg() == edges[1]->getBeg() ||
					    edges[0]->getBeg() == edges[1]->getEnd() )
					{
						n1 = edges[0]->getEnd();
						n2 = edges[0]->getBeg();
					}
					else
					{
						n1 = edges[0]->getBeg();
						n2 = edges[0]->getEnd();
					}
					//find out common edge orientation
					for(typename ElementArray<Node>::size_type j = 0; j < edges.size(); j++)
					{
						if( edges[j] == r.bridge ) //found the edge
						{
							//reverse ordering on this face
							if( r.first == n1 )
							{
								r.face->SetPrivateMarker(rev);
								reverse = true;
							}
							break;
						}
						//update edge nodes
						n1 = n2; //current end is new begin
						//find new end
						if( n2 == edges[(j+1)%edges.size()]->getBeg() )
							n2 = edges[(j+1)%edges.size()]->getEnd();
						else
							n2 = edges[(j+1)%edges.size()]->getBeg();
					}
				} while(true);
				data.RemPrivateMarker(mrk);
				Storage::real fnrm[3], fcnt[3], ccnt[3] = {0,0,0}, nnodes = 0;
				Storage::real_array x0,a,b;
				//find cell centroid
				for(typename ElementArray<T>::size_type j = 0; j < data.size(); j++)
				{
					ElementArray<Node> nodes = data[j]->getAsFace()->getNodes(mrk,true);
					nodes.SetPrivateMarker(mrk);
					for(ElementArray<Node>::iterator qt = nodes.begin(); qt != nodes.end(); ++qt)
					{
						x0 = coords->Get(qt->self());
						ccnt[0] += x0[0];
						ccnt[1] += x0[1];
						ccnt[2] += x0[2];
						nnodes += 1;
					}
				}
				ccnt[0] /= nnodes;
				ccnt[1] /= nnodes;
				ccnt[2] /= nnodes;
				for(typename ElementArray<T>::size_type j = 0; j < data.size(); j++)
				{
					//compute normal to face
					ElementArray<Node> nodes = data[j]->getAsFace()->getNodes();
					nodes.RemPrivateMarker(mrk);
					fnrm[0] = fnrm[1] = fnrm[2] = 0;
					fcnt[0] = fcnt[1] = fcnt[2] = 0;
					nnodes = 0;
					x0 = coords->Get(nodes[0]);
					a = x0;
					for(unsigned i = 0; i < nodes.size(); i++)
					{
						b = coords->Get(nodes[(i+1)%nodes.size()]);
						fnrm[0] += (a[1]-x0[1])*(b[2]-x0[2]) - (a[2]-x0[2])*(b[1]-x0[1]);
						fnrm[1] += (a[2]-x0[2])*(b[0]-x0[0]) - (a[0]-x0[0])*(b[2]-x0[2]);
						fnrm[2] += (a[0]-x0[0])*(b[1]-x0[1]) - (a[1]-x0[1])*(b[0]-x0[0]);
						
						fcnt[0] += a[0];
						fcnt[1] += a[1];
						fcnt[2] += a[2];
						
						a.swap(b);
						
						nnodes += 1;
					}
					fnrm[0] *= 0.5;
					fnrm[1] *= 0.5;
					fnrm[2] *= 0.5;
					fcnt[0] /= nnodes;
					fcnt[1] /= nnodes;
					fcnt[2] /= nnodes;
					for(int r = 0; r < 3; ++r)
						fcnt[r] = fcnt[r]-ccnt[r];
					measure += (data[j]->GetPrivateMarker(rev) ? -1.0 : 1.0)*dot_prod(fcnt,fnrm);
				}
				mesh->ReleasePrivateMarker(mrk);
				data.RemPrivateMarker(rev);
				mesh->ReleasePrivateMarker(rev);
				measure /= 3.0;
				measure = fabs(measure);
			}
			return measure;
		}
		void recursive_find(unsigned node, unsigned length)
		{
#if !defined(DEEP_RECURSE)
			if( !min_loop.empty() && length > min_loop.size() ) return;
#endif
			bool success = false;
			
			if( print )
			{
				std::cout << "current loop:";
				for(unsigned j = 0; j < insert_order.size(); j++)
					std::cout << " " << ElementTypeName(GetHandleElementType(head_column[insert_order[j]])) << ":" << GetHandleID(head_column[insert_order[j]]);
				std::cout << " try add " << ElementTypeName(GetHandleElementType(head_column[node])) << ":" << GetHandleID(head_column[node]) << " (" <<(int)hide_column[node] << ") ";
			}
			if( do_show_row(node) )
			{
				//if( print ) std::cout << "ok" << std::endl;
				success = test_success();
				
				if( success )
				{
					//if( min_loop.empty() || min_loop.size() >= length )
					{
						
						
						temp_loop.resize(length);
						for(unsigned j = 0; j < insert_order.size(); j++)
							temp_loop.at(j) = head_column[insert_order[j]];
						Storage::real measure = compute_measure(temp_loop);
						
						if( print )
						{
							std::cout << "found loop [" << temp_loop.size() <<"]: ";
							for(unsigned j = 0; j < temp_loop.size(); ++j)
								std::cout << ElementTypeName(GetHandleElementType(head_column[insert_order[j]])) << ":" << GetHandleID(head_column[insert_order[j]]) << " ";
							std::cout << "measure " << measure << std::endl;
						}
						
						if( min_loop.empty() || min_loop_measure + MEASURE_EPS >= measure )
						{
							if( !(fabs(min_loop_measure-measure) < MEASURE_EPS && min_loop.size() < temp_loop.size()) )
							{
								min_loop.swap(temp_loop);
								min_loop_measure = measure;
								if( print ) std::cout << "selected as current loop" << std::endl;
							}
							//~ if( min_loop.size() == head_column.size() ) // all elements were visited
							//~ {
							//~ unsigned num = 0;
							//~ for(unsigned j = 0; j < head_row.size(); j++) //check that all bridge elements were visited - we don't have any other loop then
							//~ num += hide_row[j];
							//~ if( num == head_row.size() ) exit_recurse = true; //exit recursive loop
							//~ }
						}
					}
				}
				else
				{
					bool stub = false;
					for(dynarray<unsigned char,256>::size_type j = 0; j < head_row_count.size() && !exit_recurse; j++) //first try follow the order
					{
						if( stub_row[j] == 0 && matrix[node*head_row_count.size()+j] == 1 && head_row_count[j] == 1 )
						{
							for(dynarray<HandleType,256>::size_type q = 0; q < head_column.size() && !exit_recurse; q++)
							{
								if( visits[q] > 0 && matrix[q*head_row_count.size()+j] == 1 && hide_column[q] == 1 )
								{
									recursive_find(static_cast<unsigned>(q),length+1);
								}
							}
							if( head_row_count[j] == 1 )
							{
								stub_row[j] = 1;
								stub = true;
								break; //this is a stub path
							}
						}
					}
					
					if( !stub ) for(dynarray<unsigned char,256>::size_type j = 0; j < head_row_count.size() && !exit_recurse; j++)
					{
						if( stub_row[j] == 0 && matrix[node*head_row_count.size()+j] == 0 && head_row_count[j] == 1 )
						{
							for(dynarray<HandleType,256>::size_type q = 0; q < head_column.size() && !exit_recurse; q++)
							{
								if( visits[q] > 0 && matrix[q*head_row_count.size()+j] == 1 && hide_column[q] == 1 )
								{
									recursive_find(static_cast<unsigned>(q),length+1);
								}
							}
							if( head_row_count[j] == 1 )
							{
								stub_row[j] = 1;
								stub = true;
								break; //this is a stub path
							}
						}
					}
					
				}
				do_hide_row(node);
			}
			//else if( print ) std::cout << "fail" << std::endl;
			if( length == 1 )
			{
				for(dynarray<HandleType,256>::size_type j = 0; j < head_row.size(); j++)
					stub_row[j] = 0;
			}
		}
	public:
		bool all_visited()
		{
			for(dynarray<char,256>::size_type k = 0; k < visits.size(); k++)
				if( visits[k] != 0 ) return false;
			return true;
		}
		void print_matrix()
		{
			Storage::real cnt[3];
			for(dynarray<HandleType,256>::size_type k = 0; k < head_column.size(); k++)
			{
				for(dynarray<HandleType,256>::size_type j = 0; j < head_row.size(); j++)
					std::cout << static_cast<int>(matrix[k*head_row.size()+ j]);
				std::cout << " " << (int)visits[k] << " " << (int)visits0[k];
				Element(mesh,head_column[k])->Centroid(cnt);
				std::cout << " " << cnt[0] << " " << cnt[1] << " " << cnt[2];
				std::cout << " " << ElementTypeName(GetHandleElementType(head_column[k]));
				std::cout << ":" << GetHandleID(head_column[k]);
				std::cout << std::endl;
			}
#if defined(RECORD_PATH)
			std::cout << "loops [" << remember.size() << "]:" << std::endl;
			for(size_t k = 0; k < remember.size(); ++k)
			{
				std::cout << k << " size " << remember[k].first.size() << ": ";
				for(size_t l = 0; l < remember[k].first.size()-1; ++l)
					std::cout << ElementTypeName(GetHandleElementType(head_column[remember[k].first[l]])) << ":" << GetHandleID(head_column[remember[k].first[l]]) << ", ";
				std::cout << ElementTypeName(GetHandleElementType(head_column[remember[k].first.back()])) << ":" << GetHandleID(head_column[remember[k].first.back()]) << " measure " << remember[k].second << std::endl;
			}
#endif
			
			if( GetHandleElementType(head_column[0]) == EDGE )
			{
				std::cout << "edges [" << head_column.size() << "]" << std::endl;
				for(int k = 0; k < (int)head_column.size(); ++k)
				{
					Edge e(mesh,head_column[k]);
					std::cout << "(" << e->getBeg()->Coords()[0] << "," << e->getBeg()->Coords()[1] << "," << e->getBeg()->Coords()[2] << ") <-> (" << e->getEnd()->Coords()[0] << "," << e->getEnd()->Coords()[1] << "," << e->getEnd()->Coords()[2] << ")" << std::endl;
				}
				std::cout << "edges [" << head_column.size() << "]" << std::endl;
				for(int k = 0; k < (int)head_column.size(); ++k)
				{
					Edge e(mesh,head_column[k]);
					std::cout << e->getBeg()->GetHandle() << " <-> " << e->getEnd()->GetHandle() << std::endl;
				}
			}
		}
		template<typename InputIterator>
		incident_matrix(Mesh * mesh, InputIterator beg, InputIterator end, typename ElementArray<T>::size_type num_inner, const AbstractCoords & _coords = GridCoords(), bool print = false)
		: mesh(mesh), head_column(beg,end), min_loop(), print(print), coords(_coords.Copy())
		{
			min_loop.SetMeshLink(mesh);
			temp_loop.SetMeshLink(mesh);
			//isInputForwardIterators<T,InputIterator>();
			if( !head_column.empty() )
			{
				MarkerType hide_marker = mesh->CreatePrivateMarker();
				
				visits.resize(head_column.size());
				visits0.resize(head_column.size());
				/*
				for(typename dynarray<HandleType, 256>::size_type it = 0; it < head_column.size(); ++it)
				{
					Element::adj_type const & sub = mesh->LowConn(head_column[it]);
					for(Element::adj_type::size_type jt = 0; jt < sub.size(); ++jt)
						if( mesh->GetPrivateMarker(sub[jt],hide_marker) )
							printf("element %s:%d already have marker %d\n",ElementTypeName(GetHandleElementType(sub[jt])),GetHandleID(sub[jt]),hide_marker);
				}
				 */
				for(typename dynarray<HandleType, 256>::size_type it = 0; it < head_column.size(); ++it)
				{
					visits[it] = it < num_inner ? 2 : 1;
					visits0[it] = visits[it];
					Element::adj_type const & sub = mesh->LowConn(head_column[it]);
					for(Element::adj_type::size_type jt = 0; jt < sub.size(); ++jt)
					{
						if( !mesh->GetPrivateMarker(sub[jt],hide_marker) )
						{
							head_row.push_back(sub[jt]);
							mesh->SetPrivateMarker(sub[jt],hide_marker);
						}
					}
				}
				std::map<HandleType,int> mat_num;
				for(dynarray<HandleType,256>::size_type it = 0; it < head_row.size(); ++it)
				{
					mesh->RemPrivateMarker(head_row[it],hide_marker);
					mat_num[head_row[it]] = static_cast<int>(it);
				}
				
				mesh->ReleasePrivateMarker(hide_marker);
				
				matrix.resize(head_row.size()*head_column.size(),0);
				
				
				
				for(typename dynarray<HandleType,256>::size_type it = 0; it < head_column.size(); ++it)
				{
					Element::adj_type const & sub = mesh->LowConn(head_column[it]);
					for(Element::adj_type::size_type jt = 0; jt < sub.size(); ++jt)
					{
						matrix[(it)*head_row.size()+mat_num[sub[jt]]] = 1;
					}
				}
				head_row_count.resize(head_row.size(),0);
				stub_row.resize(head_row.size(),0);
				hide_row.resize(head_row.size(),1);
				hide_column.resize(head_column.size(),1);
			}
		}
		incident_matrix(const incident_matrix & other)
		: mesh(other.mesh), matrix(other.matrix), head_column(other.head_column), head_row(other.head_row),
		head_row_count(other.head_row_count), min_loop(other.min_loop),
		hide_row(other.hide_row), hide_column(other.hide_column),
		stub_row(other.stub_row), coords(other.coords->Copy())
		{
		}
		incident_matrix & operator =(const incident_matrix & other)
		{
			mesh = other.mesh;
			matrix = other.matrix;
			head_column = other.head_column;
			head_row = other.head_row;
			head_row_count = other.head_row_count;
			min_loop = other.min_loop;
			hide_row = other.hide_row;
			hide_column = other.hide_column;
			stub_row = other.stub_row;
			coords = other.coords->Copy();
			return *this;
		}
		~incident_matrix()
		{
			delete coords;
		}
		bool find_shortest_loop(ElementArray<T> & ret)
		{
			ret.clear();
			exit_recurse = false;
			unsigned first = UINT_MAX;
			do
			{
				first = UINT_MAX;
				for(unsigned q = 0; q < head_column.size(); q++)
					if( visits[q] == 1 )
					{
						first = q;
						break;
					}
				if( first != UINT_MAX )
				{
					min_loop_measure = 1.0e+20;
					recursive_find(first,1);
					if( min_loop.empty() )
					{
						if( print ) std::cout << "abandon " << first << std::endl;
#if defined(RECORD_PATH)
						remember.push_back(std::make_pair(std::vector<int>(1,first),-1.0));
#endif
						visits[first]--; //don't start again from this element
					}
				}
			} while( min_loop.empty() && first != UINT_MAX );
			

			
			ret.insert(ret.end(),min_loop.begin(),min_loop.end());
			min_loop.clear();
			
			if( !ret.empty() )
			{
#if defined(RECORD_PATH)
				std::vector<int> add_remember;
				add_remember.reserve(min_loop.size());
#endif
				MarkerType hide_marker = mesh->CreatePrivateMarker();
				for(typename ElementArray<T>::size_type k = 0; k < ret.size(); k++) mesh->SetPrivateMarker(ret.at(k),hide_marker);
				if( print ) std::cout << "return loop [" << ret.size() << "]:";
				for(dynarray<HandleType,256>::size_type k = 0; k < head_column.size(); k++)
					if( mesh->GetPrivateMarker(head_column[k],hide_marker) )
					{
						visits[k]--;
#if defined(RECORD_PATH)
						add_remember.push_back((int)k);
#endif
						if( print ) std::cout << ElementTypeName(GetHandleElementType(head_column[k])) << ":" << GetHandleID(head_column[k]) << " ";
					}
				if( print ) std::cout << std::endl;
				for(typename ElementArray<T>::size_type k = 0; k < ret.size(); k++) mesh->RemPrivateMarker(ret.at(k),hide_marker);
				mesh->ReleasePrivateMarker(hide_marker);
#if defined(RECORD_PATH)
				remember.push_back(std::make_pair(add_remember,min_loop_measure));
#endif
				return true;
			}
			return false;
		}
		
		Element get_element(unsigned k) {return Element(mesh,head_column[k]);}
		void set_visit(unsigned k, char vis ) { visits[k] = vis; }
	};

#else // NEW_ALGORITHM
	template<class T>
	class incident_matrix
	{
		Mesh * mesh;
		dynarray< unsigned char, 4096 > matrix;
		dynarray< char ,256 > visits;
		dynarray< char ,256 > visits0;
		dynarray<HandleType, 256> head_column;
		dynarray<HandleType, 256> head_row;
		dynarray<unsigned char ,256> head_row_count;
		dynarray<unsigned, 256> insert_order;
		bool exit_recurse;
		ElementArray<T> min_loop, temp_loop; //used as return
		dynarray< char , 256 > hide_column;
		dynarray< char , 256 > hide_row;
		dynarray< char , 256 > stub_row;
#if defined(RECORD_PATH)
		std::vector< std::pair<std::vector<int>,double> > remember;
#endif
		double min_loop_measure;
		bool print;
		AbstractCoords * coords;
		
		bool do_hide_row(unsigned k)
		{
			if( hide_column[k] == 0 )
			{
				hide_column[k] = 1;
				for(unsigned i = 0; i < head_row_count.size(); i++)
					if( matrix[k*head_row_count.size()+i] == 1 )
					{
						head_row_count[i] -= 1;
						if( head_row_count[i] == 0 )
						{
							hide_row[i] = 1;
							stub_row[i] = 0;
						}
					}
				insert_order.pop_back();
			}
			return true;
		}
		
		bool do_show_row(unsigned k)
		{
			if( hide_column[k] == 1 )
			{
				hide_column[k] = 0;
				
				bool success = true;
				for(unsigned i = 0; i < head_row_count.size(); i++)
					if( matrix[k*head_row_count.size()+i] == 1 )
					{
						head_row_count[i] += 1;
						if( head_row_count[i] > 0 ) hide_row[i] = 0;
						if( head_row_count[i] > 2 ) success = false;
					}
				insert_order.push_back(k);
				if( !success ) do_hide_row(k);
				return success;
				
			} else return true;
		}
		bool test_success()
		{
			bool success = true;
			for(unsigned j = 0; j < head_row_count.size(); j++)
			{
				if( head_row_count[j] == 1 )
				{
					success = false;
					break;
				}
			}
			return success;
		}
		Storage::real compute_measure(ElementArray<T> & data)
		{
			Storage::real measure = 0;
			if( data[0]->GetElementDimension() == 1 ) //this is edge //use geometric dimension here for 2d compatibility
			{
				//calculate area
				//int mdim = data[0]->GetMeshLink()->GetDimensions();
				ElementArray<Node> nodes,n1,n2;
				n1 = data[0]->getNodes();
				n2 = data[1]->getNodes();
				if( n1[0] == n2[0] || n1[0] == n2[1])
				{
					nodes.push_back(n1[1]);
					nodes.push_back(n1[0]);
				}
				else
				{
					nodes.push_back(n1[0]);
					nodes.push_back(n1[1]);
				}
				for(typename ElementArray<T>::size_type j = 1; j < data.size(); j++)
				{
					n1 = data[j]->getNodes();
					assert((nodes.back() == n1[0] || nodes.back() == n1[1]));
					if( nodes.back() == n1[0] )
						nodes.push_back(n1[1]);
					else
						nodes.push_back(n1[0]);
				}
				
				Storage::real x[3] = {0,0,0};
				Storage::real_array x0 = coords->Get(nodes[0]);
				for(unsigned i = 1; i < nodes.size()-1; i++)
				{
					Storage::real_array v1 = coords->Get(nodes[i]);
					Storage::real_array v2 = coords->Get(nodes[i+1]);
					if( v1.size() == 3 && v2.size() == 3 )
					{
						x[0] += (v1[1]-x0[1])*(v2[2]-x0[2]) - (v1[2]-x0[2])*(v2[1]-x0[1]);
						x[1] += (v1[2]-x0[2])*(v2[0]-x0[0]) - (v1[0]-x0[0])*(v2[2]-x0[2]);
					}
					x[2] += (v1[0]-x0[0])*(v2[1]-x0[1]) - (v1[1]-x0[1])*(v2[0]-x0[0]);
				}
				measure = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])*0.5;
				
			}
			else //this is 3d face
			{
				//firstly, have to figure out orientation of each face
				//mark all faces, so that we can perform adjacency retrival
				MarkerType mrk = mesh->CreatePrivateMarker();
				MarkerType rev = mesh->CreatePrivateMarker(); //reverse orientation
				for(int k = 1; k < data.size(); ++k)
					data[k]->SetPrivateMarker(mrk); //0-th face orientation is default
				Node n1,n2; //to retrive edge
				bool reverse = false; //reverse orientation in considered face
				std::deque< orient_face > stack; //edge and first node and face for visiting
				//todo: can do faster by retriving edges and going over their nodes
				//should not use FindSharedAdjacency
				ElementArray<Edge> edges = data[0]->getEdges();
				do
				{
					//figure out starting node order
					if( edges[0]->getBeg() == edges[1]->getBeg() ||
					   edges[0]->getBeg() == edges[1]->getEnd() )
					{
						n1 = edges[0]->getEnd();
						n2 = edges[0]->getBeg();
					}
					else
					{
						n1 = edges[0]->getBeg();
						n2 = edges[0]->getEnd();
					}
					//schedule unvisited adjacent faces
					for(typename ElementArray<Edge>::size_type j = 0; j < edges.size(); j++)
					{
						//schedule face adjacent to considered edge
						ElementArray<Face> adjacent = edges[j]->getFaces(mrk);
						assert(adjacent.size() <= 1);
						if( !adjacent.empty() && adjacent[0].GetPrivateMarker(mrk))
						{
							adjacent.RemPrivateMarker(mrk);
							stack.push_back(orient_face(edges[j],reverse ? n2 : n1,adjacent[0]));
						}
						//update edge nodes
						n1 = n2; //current end is new begin
						//find new end
						if( n2 == edges[(j+1)%edges.size()]->getBeg() )
							n2 = edges[(j+1)%edges.size()]->getEnd();
						else
							n2 = edges[(j+1)%edges.size()]->getBeg();
					}
					if( stack.empty() ) break;
					//get entry from stack
					orient_face r = stack.front();
					//remove face from stack
					stack.pop_front();
					//retrive edges for new face
					edges = r.face->getEdges();
					reverse = false;
					//figure out starting node order
					if( edges[0]->getBeg() == edges[1]->getBeg() ||
					   edges[0]->getBeg() == edges[1]->getEnd() )
					{
						n1 = edges[0]->getEnd();
						n2 = edges[0]->getBeg();
					}
					else
					{
						n1 = edges[0]->getBeg();
						n2 = edges[0]->getEnd();
					}
					//find out common edge orientation
					for(typename ElementArray<Node>::size_type j = 0; j < edges.size(); j++)
					{
						if( edges[j] == r.bridge ) //found the edge
						{
							//reverse ordering on this face
							if( r.first == n1 )
							{
								r.face->SetPrivateMarker(rev);
								reverse = true;
							}
							break;
						}
						//update edge nodes
						n1 = n2; //current end is new begin
						//find new end
						if( n2 == edges[(j+1)%edges.size()]->getBeg() )
							n2 = edges[(j+1)%edges.size()]->getEnd();
						else
							n2 = edges[(j+1)%edges.size()]->getBeg();
					}
				} while(true);
				data.RemPrivateMarker(mrk);
				Storage::real fnrm[3], fcnt[3], ccnt[3] = {0,0,0}, nnodes = 0;
				Storage::real_array x0,a,b;
				//find cell centroid
				for(typename ElementArray<T>::size_type j = 0; j < data.size(); j++)
				{
					ElementArray<Node> nodes = data[j]->getAsFace()->getNodes(mrk,true);
					nodes.SetPrivateMarker(mrk);
					for(ElementArray<Node>::iterator qt = nodes.begin(); qt != nodes.end(); ++qt)
					{
						x0 = coords->Get(qt->self());
						ccnt[0] += x0[0];
						ccnt[1] += x0[1];
						ccnt[2] += x0[2];
						nnodes += 1;
					}
				}
				ccnt[0] /= nnodes;
				ccnt[1] /= nnodes;
				ccnt[2] /= nnodes;
				for(typename ElementArray<T>::size_type j = 0; j < data.size(); j++)
				{
					//compute normal to face
					ElementArray<Node> nodes = data[j]->getAsFace()->getNodes();
					nodes.RemPrivateMarker(mrk);
					fnrm[0] = fnrm[1] = fnrm[2] = 0;
					fcnt[0] = fcnt[1] = fcnt[2] = 0;
					nnodes = 0;
					x0 = coords->Get(nodes[0]);
					a = x0;
					for(unsigned i = 0; i < nodes.size(); i++)
					{
						b = coords->Get(nodes[(i+1)%nodes.size()]);
						fnrm[0] += (a[1]-x0[1])*(b[2]-x0[2]) - (a[2]-x0[2])*(b[1]-x0[1]);
						fnrm[1] += (a[2]-x0[2])*(b[0]-x0[0]) - (a[0]-x0[0])*(b[2]-x0[2]);
						fnrm[2] += (a[0]-x0[0])*(b[1]-x0[1]) - (a[1]-x0[1])*(b[0]-x0[0]);
						
						fcnt[0] += a[0];
						fcnt[1] += a[1];
						fcnt[2] += a[2];
						
						a.swap(b);
						
						nnodes += 1;
					}
					fnrm[0] *= 0.5;
					fnrm[1] *= 0.5;
					fnrm[2] *= 0.5;
					fcnt[0] /= nnodes;
					fcnt[1] /= nnodes;
					fcnt[2] /= nnodes;
					for(int r = 0; r < 3; ++r)
						fcnt[r] = fcnt[r]-ccnt[r];
					measure += (data[j]->GetPrivateMarker(rev) ? -1.0 : 1.0)*dot_prod(fcnt,fnrm);
				}
				mesh->ReleasePrivateMarker(mrk);
				data.RemPrivateMarker(rev);
				mesh->ReleasePrivateMarker(rev);
				measure /= 3.0;
				measure = fabs(measure);
			}
			return measure;
		}
		void recursive_find(unsigned node, unsigned length)
		{
			if( !min_loop.empty() && length > min_loop.size() ) return;
			bool success = false;
			if( do_show_row(node) )
			{
				success = test_success();
				
				if( success )
				{
					//if( min_loop.empty() || min_loop.size() >= length )
					{
						
						
						temp_loop.resize(length);
						for(unsigned j = 0; j < insert_order.size(); j++)
							temp_loop.at(j) = head_column[insert_order[j]];
						Storage::real measure = compute_measure(temp_loop);
						
						if( print )
						{
							std::cout << "found loop [" << temp_loop.size() <<"]: ";
							for(unsigned j = 0; j < temp_loop.size(); ++j)
								std::cout << ElementTypeName(GetHandleElementType(head_column[insert_order[j]])) << ":" << GetHandleID(head_column[insert_order[j]]) << " ";
							std::cout << "measure " << measure << std::endl;
						}
						
						if( min_loop.empty() || min_loop_measure + MEASURE_EPS >= measure )
						{
							if( !(fabs(min_loop_measure-measure) < MEASURE_EPS && min_loop.size() < temp_loop.size()) )
							{
								min_loop.swap(temp_loop);
								min_loop_measure = measure;
								if( print ) std::cout << "selected as current loop" << std::endl;
							}
							//~ if( min_loop.size() == head_column.size() ) // all elements were visited
							//~ {
							//~ unsigned num = 0;
							//~ for(unsigned j = 0; j < head_row.size(); j++) //check that all bridge elements were visited - we don't have any other loop then
							//~ num += hide_row[j];
							//~ if( num == head_row.size() ) exit_recurse = true; //exit recursive loop
							//~ }
						}
					}
				}
				else
				{
					bool stub = false;
					for(dynarray<unsigned char,256>::size_type j = 0; j < head_row_count.size() && !exit_recurse; j++) //first try follow the order
					{
						if( stub_row[j] == 0 && matrix[node*head_row_count.size()+j] == 1 && head_row_count[j] == 1 )
						{
							for(dynarray<HandleType,256>::size_type q = 0; q < head_column.size() && !exit_recurse; q++)
							{
								if( visits[q] > 0 && matrix[q*head_row_count.size()+j] == 1 && hide_column[q] == 1 )
								{
									recursive_find(static_cast<unsigned>(q),length+1);
								}
							}
							if( head_row_count[j] == 1 )
							{
								stub_row[j] = 1;
								stub = true;
								break; //this is a stub path
							}
						}
					}
					
					if( !stub ) for(dynarray<unsigned char,256>::size_type j = 0; j < head_row_count.size() && !exit_recurse; j++)
					{
						if( stub_row[j] == 0 && matrix[node*head_row_count.size()+j] == 0 && head_row_count[j] == 1 )
						{
							for(dynarray<HandleType,256>::size_type q = 0; q < head_column.size() && !exit_recurse; q++)
							{
								if( visits[q] > 0 && matrix[q*head_row_count.size()+j] == 1 && hide_column[q] == 1 )
								{
									recursive_find(static_cast<unsigned>(q),length+1);
								}
							}
							if( head_row_count[j] == 1 )
							{
								stub_row[j] = 1;
								stub = true;
								break; //this is a stub path
							}
						}
					}
					
				}
				do_hide_row(node);
			}
			if( length == 1 )
			{
				for(dynarray<HandleType,256>::size_type j = 0; j < head_row.size(); j++)
					stub_row[j] = 0;
			}
		}
	public:
		bool all_visited()
		{
			for(dynarray<char,256>::size_type k = 0; k < visits.size(); k++)
				if( visits[k] != 0 ) return false;
			return true;
		}
		void print_matrix()
		{
			Storage::real cnt[3];
			for(dynarray<HandleType,256>::size_type k = 0; k < head_column.size(); k++)
			{
				for(dynarray<HandleType,256>::size_type j = 0; j < head_row.size(); j++)
					std::cout << static_cast<int>(matrix[k*head_row.size()+ j]);
				std::cout << " " << (int)visits[k] << " " << (int)visits0[k];
				Element(mesh,head_column[k])->Centroid(cnt);
				std::cout << " " << cnt[0] << " " << cnt[1] << " " << cnt[2];
				std::cout << std::endl;
			}
#if defined(RECORD_PATH)
			std::cout << "loops [" << remember.size() << "]:" << std::endl;
			for(size_t k = 0; k < remember.size(); ++k)
			{
				std::cout << k << " size " << remember[k].first.size() << ": ";
				for(size_t l = 0; l < remember[k].first.size()-1; ++l)
					std::cout << remember[k].first[l] << ", ";
				std::cout << remember[k].first.back() << " measure " << remember[k].second << std::endl;
			}
#endif
			
			if( GetHandleElementType(head_column[0]) == EDGE )
			{
				std::cout << "edges [" << head_column.size() << "]" << std::endl;
				for(int k = 0; k < (int)head_column.size(); ++k)
				{
					Edge e(mesh,head_column[k]);
					std::cout << "(" << e->getBeg()->Coords()[0] << "," << e->getBeg()->Coords()[1] << "," << e->getBeg()->Coords()[2] << ") <-> (" << e->getEnd()->Coords()[0] << "," << e->getEnd()->Coords()[1] << "," << e->getEnd()->Coords()[2] << ")" << std::endl;
				}
				std::cout << "edges [" << head_column.size() << "]" << std::endl;
				for(int k = 0; k < (int)head_column.size(); ++k)
				{
					Edge e(mesh,head_column[k]);
					std::cout << e->getBeg()->GetHandle() << " <-> " << e->getEnd()->GetHandle() << std::endl;
				}
			}
		}
		void MapData(const dynarray<HandleType,256> & set)
		{
			
		}
		template<typename InputIterator>
		incident_matrix(Mesh * mesh, InputIterator beg, InputIterator end, typename ElementArray<T>::size_type num_inner, const AbstractCoords & _coords = GridCoords(), bool print = false)
		: mesh(mesh), head_column(beg,end), min_loop(), print(print), coords(_coords.Copy())
		{
			min_loop.SetMeshLink(mesh);
			temp_loop.SetMeshLink(mesh);
			//isInputForwardIterators<T,InputIterator>();
			if( !head_column.empty() )
			{
				MarkerType hide_marker = mesh->CreatePrivateMarker();
				
				visits.resize(head_column.size());
				visits0.resize(head_column.size());
				/*
				 for(typename dynarray<HandleType, 256>::size_type it = 0; it < head_column.size(); ++it)
				 {
					Element::adj_type const & sub = mesh->LowConn(head_column[it]);
					for(Element::adj_type::size_type jt = 0; jt < sub.size(); ++jt)
				 if( mesh->GetPrivateMarker(sub[jt],hide_marker) )
				 printf("element %s:%d already have marker %d\n",ElementTypeName(GetHandleElementType(sub[jt])),GetHandleID(sub[jt]),hide_marker);
				 }
				 */
				for(typename dynarray<HandleType, 256>::size_type it = 0; it < head_column.size(); ++it)
				{
					visits[it] = it < num_inner ? 2 : 1;
					visits0[it] = visits[it];
					Element::adj_type const & sub = mesh->LowConn(head_column[it]);
					for(Element::adj_type::size_type jt = 0; jt < sub.size(); ++jt)
					{
						if( !mesh->GetPrivateMarker(sub[jt],hide_marker) )
						{
							head_row.push_back(sub[jt]);
							mesh->SetPrivateMarker(sub[jt],hide_marker);
						}
					}
				}
				std::map<HandleType,int> mat_num;
				for(dynarray<HandleType,256>::size_type it = 0; it < head_row.size(); ++it)
				{
					mesh->RemPrivateMarker(head_row[it],hide_marker);
					mat_num[head_row[it]] = static_cast<int>(it);
				}
				
				mesh->ReleasePrivateMarker(hide_marker);
				
				matrix.resize(head_row.size()*head_column.size(),0);
				
				
				
				for(typename dynarray<HandleType,256>::size_type it = 0; it < head_column.size(); ++it)
				{
					Element::adj_type const & sub = mesh->LowConn(head_column[it]);
					for(Element::adj_type::size_type jt = 0; jt < sub.size(); ++jt)
					{
						matrix[(it)*head_row.size()+mat_num[sub[jt]]] = 1;
					}
				}
				head_row_count.resize(head_row.size(),0);
				stub_row.resize(head_row.size(),0);
				hide_row.resize(head_row.size(),1);
				hide_column.resize(head_column.size(),1);
			}
		}
		incident_matrix(const incident_matrix & other)
		: mesh(other.mesh), matrix(other.matrix), head_column(other.head_column), head_row(other.head_row),
		head_row_count(other.head_row_count), min_loop(other.min_loop),
		hide_row(other.hide_row), hide_column(other.hide_column),
		stub_row(other.stub_row), coords(other.coords.Copy())
		{
		}
		incident_matrix & operator =(const incident_matrix & other)
		{
			mesh = other.mesh;
			matrix = other.matrix;
			head_column = other.head_column;
			head_row = other.head_row;
			head_row_count = other.head_row_count;
			min_loop = other.min_loop;
			hide_row = other.hide_row;
			hide_column = other.hide_column;
			stub_row = other.stub_row;
			coords = other.coords.Copy();
			return *this;
		}
		~incident_matrix()
		{
			delete coords;
		}
		bool find_shortest_loop(ElementArray<T> & ret)
		{
			ret.clear();
			exit_recurse = false;
			unsigned first = UINT_MAX;
			do
			{
				first = UINT_MAX;
				for(unsigned q = 0; q < head_column.size(); q++)
					if( visits[q] == 1 )
					{
						first = q;
						break;
					}
				if( first != UINT_MAX )
				{
					min_loop_measure = 1.0e+20;
					recursive_find(first,1);
					if( min_loop.empty() )
					{
						if( print ) std::cout << "abandon " << first << std::endl;
#if defined(RECORD_PATH)
						remember.push_back(std::make_pair(std::vector<int>(1,first),-1.0));
#endif
						visits[first]--; //don't start again from this element
					}
				}
			} while( min_loop.empty() && first != UINT_MAX );
			
			
			
			ret.insert(ret.end(),min_loop.begin(),min_loop.end());
			min_loop.clear();
			
			if( !ret.empty() )
			{
#if defined(RECORD_PATH)
				std::vector<int> add_remember;
				add_remember.reserve(min_loop.size());
#endif
				MarkerType hide_marker = mesh->CreatePrivateMarker();
				for(typename ElementArray<T>::size_type k = 0; k < ret.size(); k++) mesh->SetPrivateMarker(ret.at(k),hide_marker);
				if( print ) std::cout << "return loop [" << ret.size() << "]:";
				for(dynarray<HandleType,256>::size_type k = 0; k < head_column.size(); k++)
					if( mesh->GetPrivateMarker(head_column[k],hide_marker) )
					{
						visits[k]--;
#if defined(RECORD_PATH)
						add_remember.push_back((int)k);
#endif
						if( print ) std::cout << k << " ";
					}
				if( print ) std::cout << std::endl;
				for(typename ElementArray<T>::size_type k = 0; k < ret.size(); k++) mesh->RemPrivateMarker(ret.at(k),hide_marker);
				mesh->ReleasePrivateMarker(hide_marker);
#if defined(RECORD_PATH)
				remember.push_back(std::make_pair(add_remember,min_loop_measure));
#endif
				return true;
			}
			return false;
		}
		
		Element get_element(unsigned k) {return Element(mesh,head_column[k]);}
		void set_visit(unsigned k, char vis ) { visits[k] = vis; }
	};
	
	
#endif //NEW_ALGORITHM
}




#endif
