#include "inmost.h"
#include <stack>
#if defined(USE_MESH)

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
	
	template<class T>
	class incident_matrix
	{
		Mesh * mesh;
		dynarray< unsigned char, 4096 > matrix;
		dynarray< char ,256 > visits;
		dynarray<HandleType, 256> head_column;
		dynarray<HandleType, 256> head_row;
		dynarray<unsigned char ,256> head_row_count;
		dynarray<unsigned, 256> insert_order;
		bool exit_recurse;
		ElementArray<T> min_loop, temp_loop; //used as return
		dynarray< char , 256 > hide_column;
		dynarray< char , 256 > hide_row;
		dynarray< char , 256 > stub_row;
		dynarray< double, 192 > centroids, normals;
		double min_loop_measure;
		
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
				int mdim = data[0]->GetMeshLink()->GetDimensions();
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
					if( nodes.back() == n1[0] )
						nodes.push_back(n1[1]);
					else
						nodes.push_back(n1[0]);
				}
				
				Storage::real x[3] = {0,0,0};
				Storage::real_array x0 = nodes[0].Coords();
				for(unsigned i = 1; i < nodes.size()-1; i++)
				{
					Storage::real_array v1 = nodes[i].Coords();
					Storage::real_array v2 = nodes[i+1].Coords();
					if( mdim == 3 )
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
				data.SetPrivateMarker(mrk);
				data[0]->RemPrivateMarker(mrk); //0-th face orientation is default
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
						if( n2 == data[(j+1)%edges.size()]->getBeg() )
							n2 = data[(j+1)%edges.size()]->getEnd();
						else
							n2 = data[(j+1)%edges.size()]->getBeg();
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
						if( n2 == data[(j+1)%edges.size()]->getBeg() )
							n2 = data[(j+1)%edges.size()]->getEnd();
						else
							n2 = data[(j+1)%edges.size()]->getBeg();
					}
				} while(true);
				data.RemPrivateMarker(mrk);
				mesh->ReleasePrivateMarker(mrk);
				Storage::real v[3], d;
				for(typename ElementArray<T>::size_type j = 0; j < data.size(); j++)
				{
					d = 0;
					ElementArray<Node> nodes = data[j]->getNodes();
					if( !nodes.empty() )
					{
						Storage::real_array a = nodes[0].Coords();
						for(typename ElementArray<Node>::size_type j = 1; j < nodes.size()-1; j++)
						{
							Storage::real_array b = nodes[j].Coords();
							Storage::real_array c = nodes[j+1].Coords();
							d += __det3v(&a[0],&b[0],&c[0]);
						}
					}
					measure += (data[j]->GetPrivateMarker(rev) ? -1.0 : 1.0)*d;
				}
				data.RemPrivateMarker(rev);
				mesh->ReleasePrivateMarker(rev);
				measure /= 6.0;
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
					if( min_loop.empty() || min_loop.size() >= length )
					{
						
						
						temp_loop.resize(length);
						for(unsigned j = 0; j < insert_order.size(); j++)
							temp_loop.at(j) = head_column[insert_order[j]];
						/*
						 // TODO:
						 // MUST CORRECTLY COMPUTE VOLUMES HERE FOR CONCAVE POLYGONS/POLYHEDRONS!
						 // this should work instead of min_loop
						 bool ok = false;
						 if( temp_loop.size() <= min_loop.size() )
						 {
						 Storage::real measure = compute_measure(temp_loop);
						 
						 if( measure > 0 && measure < min_loop_measure )
						 {
							ok = true;
							min_loop_measure = measure;
						 }
						 }
						 else
						 {
						 Storage::real measure = compute_measure(temp_loop);
						 if( measure > 0 )
						 {
							ok = true;
							min_loop_measure = measure;
						 }
						 }
						 
						 if( ok )
						 */
						{
							min_loop.swap(temp_loop);
							
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
				std::cout << " " << (int)visits[k];
				Element(mesh,head_column[k])->Centroid(cnt);
				std::cout << " " << cnt[0] << " " << cnt[1] << " " << cnt[2];
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
		template<typename InputIterator>
		incident_matrix(Mesh * mesh, InputIterator beg, InputIterator end, typename ElementArray<T>::size_type num_inner)
		: mesh(mesh), head_column(beg,end), min_loop()
		{
			min_loop.SetMeshLink(mesh);
			temp_loop.SetMeshLink(mesh);
			//isInputForwardIterators<T,InputIterator>();
			if( !head_column.empty() )
			{
				MarkerType hide_marker = mesh->CreateMarker();
				
				visits.resize(head_column.size());
				for(typename dynarray<HandleType, 256>::size_type it = 0; it < head_column.size(); ++it)
				{
					visits[it] = it < num_inner ? 2 : 1;
					Element::adj_type const & sub = mesh->LowConn(head_column[it]);
					for(Element::adj_type::size_type jt = 0; jt < sub.size(); ++jt)
					{
						if( !mesh->GetMarker(sub[jt],hide_marker) )
						{
							head_row.push_back(sub[jt]);
							mesh->SetMarker(sub[jt],hide_marker);
						}
					}
				}
				tiny_map<HandleType,int,256> mat_num;
				for(dynarray<HandleType,256>::size_type it = 0; it < head_row.size(); ++it)
				{
					mesh->RemMarker(head_row[it],hide_marker);
					mat_num[head_row[it]] = static_cast<int>(it);
				}
				
				mesh->ReleaseMarker(hide_marker);
				
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
		stub_row(other.stub_row)
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
			return *this;
		}
		~incident_matrix()
		{
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
					recursive_find(first,1);
					if( min_loop.empty() )
						visits[first]--; //don't start again from this element
				}
			} while( min_loop.empty() && first != UINT_MAX );
			
			ret.insert(ret.end(),min_loop.begin(),min_loop.end());
			min_loop.clear();
			
			if( !ret.empty() )
			{
				MarkerType hide_marker = mesh->CreateMarker();
				for(typename ElementArray<T>::size_type k = 0; k < ret.size(); k++) mesh->SetMarker(ret.at(k),hide_marker);
				for(dynarray<HandleType,256>::size_type k = 0; k < head_column.size(); k++)
					if( mesh->GetMarker(head_column[k],hide_marker) ) visits[k]--;
				for(typename ElementArray<T>::size_type k = 0; k < ret.size(); k++) mesh->RemMarker(ret.at(k),hide_marker);
				mesh->ReleaseMarker(hide_marker);
				return true;
			}
			return false;
		}
		
		Element get_element(unsigned k) {return Element(mesh,head_column[k]);}
		void set_visit(unsigned k, char vis ) { visits[k] = vis; }
	};
}




#endif
