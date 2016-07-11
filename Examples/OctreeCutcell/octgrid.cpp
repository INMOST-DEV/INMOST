
#include "octgrid.h"
#include <math.h>
#include <new>
#include <deque>
extern bool allow_coarse;
extern bool allow_refine;

extern bool global_problem;
const bool unite_faces = false;
const INMOST_DATA_ENUM_TYPE multiparent = ENUMUNDEF; //set ENUMUNDEF to enable union of small cells
bool remove_orphan_elements = true;
extern int global_test_number;
#define MINVOLFRAC 0.001

void default_transformation(double xyz[3]) { (void) xyz; }
int default_cell_should_unite(struct grid * g, int cell) { (void) g; (void) cell; return 0; }
int default_cell_should_split(struct grid * g, int cell) { (void) g; (void) cell; return 0; }
void default_cell_unite_data(struct grid * g, int cell) { (void) g; (void) cell; }
void default_cell_split_data(struct grid * g, int cell) { (void) g; (void) cell; }
void default_vert_interpolate_data(struct grid * g,int big_cell, int nvert, int * verts, int * isnew) { (void) g; (void) big_cell; (void ) nvert; (void ) verts; (void ) isnew;  }
void default_vert_init_data(struct grid * g, int vert) {(void) g; (void) vert;};
void default_cell_init_data(struct grid * g, int cell) {(void) g; (void) cell;};
void default_vert_destroy_data(struct grid * g, int vert) {(void) g; (void) vert;};
void default_cell_destroy_data(struct grid * g, int cell) {(void) g; (void) cell;};
void default_vert_to_INMOST(struct grid * g, int vert, Node v) {(void) g; (void) vert; (void) v;};
void default_cell_to_INMOST(struct grid * g, int cell, Cell r) {(void) g; (void) cell; (void) r;};
void default_init_mesh(struct grid * g) {(void ) g;} ;


void make_vec(double p1[3], double p2[3], double out[3])
{
	out[0] = p1[0] - p2[0];
	out[1] = p1[1] - p2[1];
	out[2] = p1[2] - p2[2];
}

void cross_prod(double v1[3], double v2[3], double out[3])
{
	out[0] = v1[1]*v2[2] - v1[2]*v2[1];
	out[1] = v1[2]*v2[0] - v1[0]*v2[2];
	out[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

Storage::real __det3d(Storage::real a, Storage::real b, Storage::real c,
	                         Storage::real d, Storage::real e, Storage::real f,
	                         Storage::real g, Storage::real h, Storage::real i ) 
{
	return a*e*i - c*e*g + b*f*g - a*f*h + c*d*h - b*d*i;
}
	
Storage::real __det3v(const Storage::real * x,const Storage::real * y,const Storage::real * z) 
{
	return __det3d(x[0], x[1], x[2],  y[0], y[1], y[2],  z[0], z[1], z[2]);
}

double dot_prod(double v1[3],double v2[3])
{
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

int invert_side(int side)
{
	return side/2*2 + (side%2+1)%2;
}

int add_element(int * elements, int * nelements, int newe)
{
	int i;
	for(i = 0; i < (*nelements); i++)
		if( newe == elements[i] )
		{
			//printf("already have %d at %d\n",new,i);
			return i;
		}
	//printf("elements[%d] = %d\n",nelements,new);
	elements[(*nelements)++] = newe;
	return (*nelements)-1;
}



int cellSearch(struct grid * g, int m, double v[3])
{
	int c,i;
	if( !g->cells[m].leaf )
	{
		c = 0;
		for(i = 0; i < 3; i++)
			if( v[i] > g->cells[m].center[i] )
				c += 1 << i;
		return cellSearch(g,g->cells[m].children[c],v);
	}
	return m;
}

int gridSearch(struct grid * g, double v[3])
{
	int i,j,k,l, flag;
	for(i = 0; i < g->n[0]; i++)
	for(j = 0; j < g->n[1]; j++)
	for(k = 0; k < g->n[2]; k++)
	{
		flag = 1;
		for(l = 0; l < 3 && flag; l++)
			if(!(g->cells[g->mgrid[i][j][k]].center[l] - g->cells[g->mgrid[i][j][k]].side[l]*0.5 <= v[l] &&
				 g->cells[g->mgrid[i][j][k]].center[l] + g->cells[g->mgrid[i][j][k]].side[l]*0.5 >= v[l]))
				 flag = 0;
		if( flag ) return cellSearch(g,g->mgrid[i][j][k],v);
	}
	return -1;
}


void vertGetCoord(struct grid * g, int v, double coord[3])
{
	int i;
	if( v == -1 ) return;
	if( !g->verts[v].busy ) return;
	for(i = 0; i < 1<<DIM; i++)
	{
		if( g->verts[v].env[((1<<DIM)-1)-i] != -1 && g->cells[g->verts[v].env[((1<<DIM)-1)-i]].vertexes[i] == v )
		{
			coord[0] = g->cells[g->verts[v].env[((1<<DIM)-1)-i]].center[0] + ((i & 1) * 2 - 1) * g->cells[g->verts[v].env[((1<<DIM)-1)-i]].side[0]*0.5;
			coord[1] = g->cells[g->verts[v].env[((1<<DIM)-1)-i]].center[1] + ((i & 2)     - 1) * g->cells[g->verts[v].env[((1<<DIM)-1)-i]].side[1]*0.5;
			coord[2] = g->cells[g->verts[v].env[((1<<DIM)-1)-i]].center[2] + ((i & 4) / 2 - 1) * g->cells[g->verts[v].env[((1<<DIM)-1)-i]].side[2]*0.5;
			return;
		}
	}
}

int cellAround(struct grid * g, int m, int side, int neighbours[1<<(DIM-1)])
{
	int k,v,c,q,ret = 0;
	const int vert[3][4] = {{0,2,6,4},{0,1,5,4},{0,1,3,2}};
	const int chck[3][4] = {{1,3,7,5},{2,3,7,6},{4,5,7,6}};
	int add = (side%2)*(1<<side/2);
	for(k = 0; k < 1<<(DIM-1); k++)
	{
		v = g->cells[m].vertexes[vert[side/2][k]+add];
		if( v == -1 ) TSNH;
		c = g->verts[v].env[vert[side/2][(k+DIM-1)%(1<<(DIM-1))]+add];
		if( c != -1 )
		{
			q = g->cells[c].vertexes[chck[side/2][k]-add];
			if( q != v ) c = -1;
		}
		if( c != -1 ) add_element(neighbours,&ret,c);
	}
	return ret;
}




void vertDestroyINMOST(struct grid *g, int m)
{
	//printf("%s\n",__FUNCTION__);
	if( g->verts[m].mv != InvalidHandle() )
	{
		g->mesh->Delete(g->verts[m].mv);
	}
}

void vertCreateINMOST(struct grid * g, int m)
{
	Storage::real xyz[3];
	if( g->verts[m].mv != InvalidHandle() ) return;
	vertGetCoord(g,m,xyz);
	g->transformation(xyz);
	g->verts[m].mv = g->mesh->CreateNode(xyz)->GetHandle();
	g->mesh->SetMarker(g->verts[m].mv,g->octree_node);
	//g->verts[m].mv->Integer(g->new_marker) = 1;
	g->vert_to_INMOST(g,m,Node(g->mesh,g->verts[m].mv));
}

void cellDestroyINMOST(struct grid * g, int m)
{

	if( !g->cells[m].mr->empty() )
	{
		std::vector<Storage::integer> other_del;
		for(int k = 0; k < g->cells[m].mr->size(); k++)
		{
			Cell ck = Cell(g->mesh,(*g->cells[m].mr)[k]);
			Storage::integer_array p = ck->IntegerArray(g->parent);
			for(int j = 0; j < p.size(); j++) if( p[j] != m )
			{
				other_del.push_back(p[j]); //this cell is united, should delete it's other parents
				for(int q = 0; q < g->cells[p[j]].mr->size(); q++)
					if( (*g->cells[m].mr)[k] == (*g->cells[p[j]].mr)[q] ) //remove the cell from other parent, so we don't delete it twice
					{
						g->cells[p[j]].mr->erase(g->cells[p[j]].mr->begin()+q);
						break;
					}
			}
			ElementArray<Face> faces = ck->getFaces();
			ElementArray<Edge> edges = ck->getEdges();
			ElementArray<Node> nodes = ck->getNodes();
			ck->Delete();
			for(ElementArray<Face>::iterator f = faces.begin(); f != faces.end(); f++)
			{
				if( f->nbAdjElements(CELL) == 0 ) 
					f->Delete();
			}
			for(ElementArray<Edge>::iterator e = edges.begin(); e != edges.end(); e++)
			{
				if( e->nbAdjElements(FACE) == 0 ) 
					e->Delete();
			}
			for(ElementArray<Node>::iterator e = nodes.begin(); e != nodes.end(); e++)
			{
				if( !e->GetMarker(g->octree_node) && e->nbAdjElements(EDGE) == 0)
					e->Delete();
			}
		}
		g->cells[m].mr->clear();
		std::sort(other_del.begin(),other_del.end());
		other_del.resize(std::unique(other_del.begin(),other_del.end())-other_del.begin());
		for(int k = 0; k < other_del.size(); k++)
			cellDestroyINMOST(g,other_del[k]);
	}
}

void reverse_face(int * face, bool * mid)
{
	int i;
	for(i = 0; i < face[0]/2; i++)
	{
		int temp = face[i+1];
		face[i+1] = face[face[0]-i];
		face[face[0]-i] = temp;
		
		bool btemp = mid[i+1];
		mid[i+1] = mid[face[0]-i];
		mid[face[0]-i] = btemp;
	}
}

void reverse_face2(int * face, bool * mid)
{
	int i;
	for(i = 0; i < (face[0]-1)/2; i++)
	{
		int temp = face[i+2];
		face[i+2] = face[face[0]-i];
		face[face[0]-i] = temp;
		
		bool btemp = mid[i+2];
		mid[i+2] = mid[face[0]-i];
		mid[face[0]-i] = btemp;
	}
}

int vertGetMiddle(struct grid * g, int m, int side, int edge)
{
	int i,j,current_vert,opposit_vert,middle_vert,adj_cell;
	const int nvf[6][4] = {{0,4,6,2},{1,3,7,5},{0,1,5,4},{2,6,7,3},{0,2,3,1},{4,5,7,6}};
	const int nve[6][4][2][2] = 
	{
		{{{5,6},{1,2}},{{2,7},{0,5}},{{0,3},{4,7}},{{1,4},{3,6}}},
		{{{2,7},{0,5}},{{5,6},{1,2}},{{1,4},{3,6}},{{0,3},{4,7}}},
		{{{3,5},{2,4}},{{4,7},{0,3}},{{0,6},{1,7}},{{1,2},{5,6}}},
		{{{4,7},{0,3}},{{3,5},{2,4}},{{1,2},{5,6}},{{0,6},{1,7}}},
		{{{3,6},{1,4}},{{1,7},{0,6}},{{0,5},{2,7}},{{2,4},{3,5}}},
		{{{1,7},{0,6}},{{3,6},{1,4}},{{2,4},{3,5}},{{0,5},{2,7}}}
	};
	for(i = 0; i < 2; i++)
	{
		current_vert = g->cells[m].vertexes[nvf[side][(edge+i)%4]];
		opposit_vert = g->cells[m].vertexes[nvf[side][(edge+(i+1)%2)%4]];
		for(j = 0; j < 2; j++)
		{
			adj_cell = g->verts[current_vert].env[nve[side][edge][i][j]];
			if( adj_cell != -1 && g->cells[m].level < g->cells[adj_cell].level )
			{
				middle_vert = g->cells[adj_cell].vertexes[nve[side][edge][i][1-j]];
				if( middle_vert == opposit_vert ) continue;
				return middle_vert;
			}
		}
		adj_cell = g->verts[current_vert].env[nvf[side][(edge+(i+1)%2)%4]];
		if( adj_cell != -1 && g->cells[m].level < g->cells[adj_cell].level )
		{
			middle_vert = g->cells[adj_cell].vertexes[7-nvf[side][(edge+i)%4]];
			if( middle_vert == opposit_vert ) continue;
			return middle_vert;
		}
	}
	return -1;
}

void cellGetFaceVerts(struct grid * g, int m, int side, int * nverts, Node verts[54], bool * mid,int * faces,  int reverse)
{
	int i,middle;
	const int nvf[6][4] = {{0,4,6,2},{1,3,7,5},{0,1,5,4},{2,6,7,3},{0,2,3,1},{4,5,7,6}};
	faces[0] = 0;
	mid[0] = false;
	for(i = 0; i < 4; i++)
	{
		mid[1+faces[0]] = false;
		verts[(faces[1+faces[0]] = (*nverts))] = Node(g->mesh,g->verts[g->cells[m].vertexes[nvf[side][i]]].mv);
		faces[0]++;
		(*nverts)++;
		if( (middle = vertGetMiddle(g,m,side,i)) != -1)
		{
			mid[1+faces[0]] = true;
			verts[(faces[1+faces[0]] = (*nverts))] = Node(g->mesh,g->verts[middle].mv);
			faces[0]++;
			(*nverts)++;
			
		}
	}
	if( reverse ) 
	{
		reverse_face2(faces,mid);
	}
}



std::vector<Edge> traverse_edges_sub(Edge start, Edge current, MarkerType edgeset, MarkerType visited_bridge, MarkerType visited_edge)
{
	//~ if( current == start ) return std::vector<Edge *> (1,start);
	std::vector< std::vector<Edge> > paths;
	ElementArray<Node> n = current->getNodes();
	for(int j = 0; j < n.size(); j++)
	{
		if( !n[j].GetMarker(visited_bridge) )
		{
			n[j].SetMarker(visited_bridge);
			ElementArray<Edge> e = n[j].getEdges();
			for(int i = 0; i < e.size(); i++) 
			{
				if( e[i] == start ) 
				{
					n[j].RemMarker(visited_bridge);
					return std::vector<Edge> (1,start);
				}
				if( e[i].GetMarker(edgeset) && !e[i].GetMarker(visited_edge) )
				{
					e[i].SetMarker(visited_edge);
					std::vector<Edge> ret = traverse_edges_sub(start,e[i],edgeset,visited_bridge,visited_edge);
					e[i].RemMarker(visited_edge);
					if( !ret.empty() )  
					{
						ret.push_back(e[i]);
						paths.push_back(ret);
					}
				}
			}
			n[j].RemMarker(visited_bridge);
		}
	}
	if( !paths.empty() )
	{
		int min = 0;
		for(int j = 1; j < paths.size(); j++)
		{
			if( paths[j].size() < paths[min].size() )
				min = j;
		}
		return paths[min];
	}
	return std::vector<Edge>();
}
//This function may be slow, because we collect all the arrays
//should detect shortest path here, then collect one array with shortest path
std::vector<Edge> traverse_edges(Edge start, MarkerType edgeset, MarkerType visited_bridge, MarkerType visited_edge)
{
	std::vector< std::vector<Edge> > paths;
	ElementArray<Node> n = start->getNodes();
	start->SetMarker(visited_edge);
	//~ std::cout << "start edge " << start << std::endl;
	for(int j = 0; j < n.size(); j++)
	{
		if( !n[j].GetMarker(visited_bridge) )
		{
			//~ std::cout  << "enter to the bridge " << &n[j] << " " << n[j].Coords()[0] << "," << n[j].Coords()[1] << "," << n[j].Coords()[2] << std::endl;
			n[j].SetMarker(visited_bridge);
			ElementArray<Edge> e = n[j].getEdges();
			for(int i = 0; i < e.size(); i++) 
				if( e[i].GetMarker(edgeset) && !e[i].GetMarker(visited_edge) )
				{
					e[i].SetMarker(visited_edge);
					std::vector<Edge> ret = traverse_edges_sub(start,e[i],edgeset,visited_bridge,visited_edge);
					e[i].RemMarker(visited_edge);
					if( !ret.empty() )  
					{
						ret.push_back(e[i]);
						paths.push_back(ret);
					}
				}
			n[j].RemMarker(visited_bridge);
		}
	}
	start->RemMarker(visited_edge);
	if( !paths.empty() )
	{
		int min = 0;
		for(int j = 1; j < paths.size(); j++)
		{
			if( paths[j].size() < paths[min].size() )
				min = j;
		}
		return paths[min];
	}
	return std::vector<Edge>();
}

class matcenter
{
	Storage::real xyz[3];
	mat_ret_type mat;
	int visit_count;
public:
	matcenter() :mat()
	{
		xyz[0] = xyz[1] = xyz[2] = 0;
		visit_count = 0;
	}
	matcenter(Storage::integer_array _mat, Storage::real * _xyz) : mat(_mat.begin(),_mat.end())
	{
		xyz[0] = _xyz[0];
		xyz[1] = _xyz[1];
		xyz[2] = _xyz[2];
		visit_count = 0;
	}
	matcenter(mat_ret_type const & _mat, Storage::real * _xyz, int count) : mat(_mat)
	{
		xyz[0] = _xyz[0];
		xyz[1] = _xyz[1];
		xyz[2] = _xyz[2];
		visit_count = count;
	}
	matcenter(const matcenter & other) : mat(other.mat)
	{
		xyz[0] = other.xyz[0];
		xyz[1] = other.xyz[1];
		xyz[2] = other.xyz[2];
		visit_count = other.visit_count;
	}
	matcenter & operator =(matcenter const & other)
	{
		xyz[0] = other.xyz[0];
		xyz[1] = other.xyz[1];
		xyz[2] = other.xyz[2];
		mat = other.mat;
		visit_count = other.visit_count;
		return *this;
	}
	int get_count() {return visit_count;}
	Storage::real * get_center() {return xyz;}
	mat_ret_type & get_mat() {return mat;}
	bool contain_mat(Storage::integer m) const {return std::binary_search(mat.begin(),mat.end(),m);}
	mat_ret_type intersect(mat_ret_type & mats) const
	{
		mat_ret_type intersection(std::min(mat.size(),mats.size()));
		intersection.resize(std::set_intersection(mat.begin(),mat.end(),mats.begin(),mats.end(),intersection.begin())-intersection.begin());
		return intersection;
	}
	mat_ret_type unite(mat_ret_type & mats) const
	{
		mat_ret_type intersection(mat.size()+mats.size());
		intersection.resize(std::set_union(mat.begin(),mat.end(),mats.begin(),mats.end(),intersection.begin())-intersection.begin());
		return intersection;
	}
};

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
	dynarray< unsigned char, 4096 > matrix;
	dynarray< char ,256 > visits;
	dynarray< T , 256> head_column;
	dynarray<Element, 256> head_row;
	dynarray<unsigned char ,256> head_row_count;
	dynarray<unsigned, 256> insert_order;
	bool exit_recurse;
	dynarray<T,64> min_loop, temp_loop; //used as return
	dynarray< char , 256 > hide_column;
	dynarray< char , 256 > hide_row;
	dynarray< char , 256 > stub_row;
	dynarray< double, 192 > centroids, normals;
	double min_loop_measure;
	Mesh * mesh;
	
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
	Storage::real compute_measure(dynarray<T,64> & data)
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
						adjacent[0].RemPrivateMarker(mrk);
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
			for(int k = 0; k < data.size(); ++k)
				data[k].RemPrivateMarker(mrk);
			mesh->ReleasePrivateMarker(mrk);
			Storage::real d;
			for(typename ElementArray<T>::size_type j = 0; j < data.size(); j++)
			{
				d = 0;
				ElementArray<Node> nodes = data[j]->getNodes();
				if( !nodes.empty() )
				{
					if( data[j]->GetPrivateMarker(rev) )
					{
						Storage::real_array a = nodes.back().Coords();
						for(typename ElementArray<Node>::size_type j = nodes.size()-2; j > 1; j--)
						{
							Storage::real_array b = nodes[j].Coords();
							Storage::real_array c = nodes[j-1].Coords();
							d += __det3v(&a[0],&b[0],&c[0]);
						}
					}
					else
					{
						Storage::real_array a = nodes[0].Coords();
						for(typename ElementArray<Node>::size_type j = 1; j < nodes.size()-1; j++)
						{
							Storage::real_array b = nodes[j].Coords();
							Storage::real_array c = nodes[j+1].Coords();
							d += __det3v(&a[0],&b[0],&c[0]);
						}
					}
				}
				//measure += (data[j]->GetPrivateMarker(rev) ? -1.0 : 1.0)*d;
				measure += d;
			}
			for(int k = 0; k < data.size(); ++k)
				data[k].RemPrivateMarker(rev);
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
						temp_loop[j] = head_column[insert_order[j]];
					//Storage::real measure = compute_measure(temp_loop);
					
					
					//if( min_loop.empty() || min_loop_measure >= measure )
					{
						min_loop.swap(temp_loop);
						//min_loop_measure = measure;
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
				for(unsigned j = 0; j < head_row_count.size() && !exit_recurse; j++) //first try follow the order
				{
					if( stub_row[j] == 0 && matrix[node*head_row_count.size()+j] == 1 && head_row_count[j] == 1 )
					{
						for(unsigned q = 0; q < head_column.size() && !exit_recurse; q++)
						{
							if( visits[q] > 0 && matrix[q*head_row_count.size()+j] == 1 && hide_column[q] == 1 ) 
							{
								recursive_find(q,length+1);
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
			
				if( !stub ) for(unsigned j = 0; j < head_row_count.size() && !exit_recurse; j++)
				{
					if( stub_row[j] == 0 && matrix[node*head_row_count.size()+j] == 0 && head_row_count[j] == 1 )
					{
						for(unsigned q = 0; q < head_column.size() && !exit_recurse; q++)
						{
							if( visits[q] > 0 && matrix[q*head_row_count.size()+j] == 1 && hide_column[q] == 1 ) 
							{
								recursive_find(q,length+1);
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
			for(unsigned j = 0; j < head_row.size(); j++)
				stub_row[j] = 0;
		}
	}
public:
	bool all_visited()
	{
		for(unsigned k = 0; k < visits.size(); k++)
			if( visits[k] != 0 ) return false;
		return true;
	}
	void print_matrix()
	{
		Storage::real cnt[3];
		for(unsigned k = 0; k < head_column.size(); k++)
		{
			for(unsigned j = 0; j < head_row.size(); j++)
				std::cout << static_cast<int>(matrix[k*head_row.size()+ j]);
			std::cout << " " << (int)visits[k];
			head_column[k]->Centroid(cnt);
			std::cout << " " << cnt[0] << " " << cnt[1] << " " << cnt[2];
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	template<typename InputIterator>
	incident_matrix(InputIterator beg, InputIterator end, unsigned num_inner)
	: head_column(beg,end), min_loop()
	{
		min_loop_measure = 1.0e20;
		//isInputForwardIterators<T,InputIterator>();
		if( !head_column.empty() )
		{
			Mesh * m = head_column[0]->GetMeshLink();
			mesh = m;
			MarkerType hide_marker = m->CreateMarker();

			visits.resize(head_column.size());
			for(typename dynarray<T, 256>::iterator it = head_column.begin(); it != head_column.end(); ++it)
			{
				unsigned k = it-head_column.begin();
				visits[k] = k < num_inner ? 2 : 1;
				ElementArray<Element> sub = (*it)->getAdjElements((*it)->GetElementType() >> 1);
				for(ElementArray<Element>::iterator jt = sub.begin(); jt != sub.end(); ++jt)
					if( !jt->GetMarker(hide_marker) )
					{
						head_row.push_back(jt->self());
						jt->SetMarker(hide_marker);
					}
			}
			
			
			
			tiny_map<Element,int,256> mat_num;
			
			for(dynarray<Element,256>::iterator it = head_row.begin(); it != head_row.end(); ++it)
			{
				(*it)->RemMarker(hide_marker);
				mat_num[*it] = it-head_row.begin();
			}	
			
			m->ReleaseMarker(hide_marker);
				
			matrix.resize(head_row.size()*head_column.size(),0);
			
			
			
			for(typename dynarray<T,256>::iterator it = head_column.begin(); it != head_column.end(); ++it)
			{
				ElementArray<Element> sub = (*it)->getAdjElements((*it)->GetElementType() >> 1);
				for(ElementArray<Element>::iterator jt = sub.begin(); jt != sub.end(); ++jt)
				{
					matrix[(it-head_column.begin())*head_row.size()+mat_num[jt->self()]] = 1;
				}
			}
			
			head_row_count.resize(head_row.size(),0);
			
			stub_row.resize(head_row.size(),0);
			
			hide_row.resize(head_row.size(),1);
			
			hide_column.resize(head_column.size(),1);
		}
	}
	incident_matrix(const incident_matrix & other) 
	: matrix(other.matrix), head_column(other.head_column), head_row(other.head_row), 
	  head_row_count(other.head_row_count), min_loop(other.min_loop), 
	  hide_row(other.hide_row), hide_column(other.hide_column), 
	  stub_row(other.stub_row) 
	{
	}
	incident_matrix & operator =(const incident_matrix & other) 
	{
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
		min_loop_measure = 1.0e20;
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
		
		for(typename dynarray<T,64>::iterator it = min_loop.begin(); it != min_loop.end(); ++it)
			ret.push_back(it->self());
		//ret.insert(ret.end(),min_loop.begin(),min_loop.end());
		min_loop.clear();
		
		if( !ret.empty() )
		{
			Mesh * m = ret[0]->GetMeshLink();
			MarkerType hide_marker = m->CreateMarker();
			for(unsigned k = 0; k < ret.size(); k++) ret[k]->SetMarker(hide_marker);
			for(unsigned k = 0; k < head_column.size(); k++)
				if( head_column[k]->GetMarker(hide_marker) ) visits[k]--;
			for(unsigned k = 0; k < ret.size(); k++) ret[k]->RemMarker(hide_marker);
			m->ReleaseMarker(hide_marker);
			return true;
		}
		return false;
	}
	
	Element * get_element(unsigned k) {return head_column[k];}
	void set_visit(unsigned k, char vis ) { visits[k] = vis; }
};


class edge_Comparator
	{
	private: MarkerType medge;
	public:
		edge_Comparator(MarkerType medge):medge(medge){}
		bool operator()(Edge a, Edge b){return a->GetMarker(medge) < b->GetMarker(medge);}
	};



void cellCreateINMOST(struct grid * g, int m, bool print = false)
{
	
	//if( m == 337 ) print = true;
	//if( m == 12743 ) print = true;
	//print = true;
	//if( m == 7852 || m == 7853 ) print = true;
	//if( m == 9919 ) print = true;
	//~ if( m == 3950 ) print = true;
	//~ if( (m == 63 || m == 64) && global_test_number == 7 ) print = true;
	//~ if( m == 26 && global_test_number == 24 ) print = true;

	//if( m == 113 ) print = true;
	//~ if( m == 13 ) print = true;
	//~ if( m == 190 ) print = true;
	
	int trigger_problem = 0;
	//~ if( m == 37 ) 
	//~ {
		//~ print = true;
		//~ trigger_problem = 100;
	//~ }
	ElementArray<Node> edge_nodes(g->mesh,2);
	const bool check = false;
	const int nvf[6][4] = {{0,4,6,2},{1,3,7,5},{0,1,5,4},{2,6,7,3},{0,2,3,1},{4,5,7,6}};
	int i,j,k,l,nn, neighbours[1<<(DIM-1)];
	int faces[24][9];
	bool mid[24][9];
	int sides[24];
	int dirs[24];
	Node verts[216];
	if( !g->cells[m].mr->empty() ) return;
	k = l = 0;
	for(j = 0; j < 6; j++)
	{
		if( print ) std::cout << "side " << j << std::endl;
		nn = cellAround(g,m,j,neighbours);
		if( nn == 4 )
		{
			for(i = 0; i < nn; i++)
			{
				sides[l] = j;
				dirs[l] = (j+1)%2;
				cellGetFaceVerts(g,neighbours[i],invert_side(j),&k,verts,mid[l],faces[l],dirs[l]);
				if( print ) 
				{
					std::cout << "face " << l << std::endl;
					for(unsigned jj = 0; jj < faces[l][0]; jj++) std::cout << jj << " " <<faces[l][jj+1] << " mid " << mid[l][jj+1] << " vert " << verts[faces[l][jj+1]]->LocalID() << std::endl;
				}
				
				
				l++;
			}
		}
		else 
		{
			sides[l] = j;
			dirs[l] = j%2;
			cellGetFaceVerts(g,m,j,&k,verts,mid[l],faces[l],dirs[l]);
			if( print ) 
			{
				std::cout << "face " << l << std::endl;
				for(unsigned jj = 0; jj < faces[l][0]; jj++) std::cout << jj << " " << faces[l][jj+1] << " mid " << mid[l][jj+1] << " vert " << verts[faces[l][jj+1]]->LocalID() << std::endl;
			}
			l++;
		}
	}
	
	
	//detect that cutcell is needed
	bool cutcell = false;
	Storage::integer_array mat = verts[faces[0][1]]->IntegerArray(g->materials);
	if( mat.size() != 1 ) cutcell = true;
	for(i = 0; i < l && !cutcell; i++)
	{
		for(j = 0; j < faces[i][0]; j++)
		{
			Storage::integer_array mats2 = verts[faces[i][1+j]]->IntegerArray( g->materials );
			if( mats2.size() > 1 || mats2[0] != mat[0] )
			{
				cutcell = true;
				break;
			}
		}
	}
	
	g->cells[m].vol = 0;

	for(i = 0; i < l; i++)
	{
		Storage::real_array x0 = verts[faces[i][1]]->Coords();
		Storage::real x[3];
		Storage::real y[3];
		x[0] = x[1] = x[2] = 0;
		y[0] = x0[0];
		y[1] = x0[1];
		y[2] = x0[2];
		//std::cout << "0 " << x0[0] << " " << x0[1] << " " << x0[2] << std::endl;
		for(j = 1; j < faces[i][0]-1; j++)
		{
			Storage::real_array v1 = verts[faces[i][1+j]]->Coords();
			Storage::real_array v2 = verts[faces[i][1+(j+1)%faces[i][0]]]->Coords();
			x[0] += (v1[1]-x0[1])*(v2[2]-x0[2]) - (v1[2]-x0[2])*(v2[1]-x0[1]);
			x[1] += (v1[2]-x0[2])*(v2[0]-x0[0]) - (v1[0]-x0[0])*(v2[2]-x0[2]);
			x[2] += (v1[0]-x0[0])*(v2[1]-x0[1]) - (v1[1]-x0[1])*(v2[0]-x0[0]);
			y[0] += v1[0];
			y[1] += v1[1];
			y[2] += v1[2];
			//std::cout << j << " " << v1[0] << " " << v1[1] << " " << v1[2] << std::endl;
		}
		x[0] *= 0.5;
		x[1] *= 0.5;
		x[2] *= 0.5;
		x0 = verts[faces[i][faces[i][0]]]->Coords();
		//std::cout << faces[i][0]-1 << " " << x0[0] << " " << x0[1] << " " << x0[2] << std::endl;
		y[0] = (y[0] + x0[0])/faces[i][0];
		y[1] = (y[1] + x0[1])/faces[i][0];
		y[2] = (y[2] + x0[2])/faces[i][0];
		
		//if( !dirs[i] )
		if( sides[i]%2 )
		{
			x[0] = -x[0];
			x[1] = -x[1];
			x[2] = -x[2];
		}
		//std::cout << "x " << x[0] << " " << x[1] << " " << x[2] << " y " << y[0] << " " << y[1] << " " << y[2] << " dirs " << dirs[i] << std::endl;
		g->cells[m].vol += x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
	}

	g->cells[m].vol /= 3.0;
	
	if( !cutcell )
	{
		ElementArray<Face> c_faces;
		ElementArray<Edge> f_edges(g->mesh);
		for(i = 0; i < l; i++)
		{
			f_edges.resize(faces[i][0]);
			for(j = 0; j < faces[i][0]; j++)
			{
				
				edge_nodes[0] = verts[faces[i][1+j]];
				edge_nodes[1] = verts[faces[i][1+(j+1)%faces[i][0]]];
				f_edges[j] = g->mesh->CreateEdge(edge_nodes).first;
			}
			Face new_face = g->mesh->CreateFace(f_edges).first;
			c_faces.push_back(new_face);
		}
		Cell c = g->mesh->CreateCell(c_faces).first;
		if( print ) 
		{
			std::cout << __FILE__ << ":" << __LINE__ << " new cell " << c->GetHandle() << " id " << c->LocalID() << " type " << Element::GeometricTypeName(c->GetGeometricType()) << " nodes " << c->nbAdjElements(NODE) << std::endl;
			ElementArray<Element> nodes = c->getAdjElements(NODE);
			for(j = 0; j < nodes.size(); j++)
				std::cout << "node[" << j << "]: " << nodes[j].LocalID() << std::endl;
		}
		c->Integer(g->parent) = m;
		g->cells[m].mr->push_back(c->GetHandle());
		c->Integer(g->cell_material) = mat[0];
	}
	
	else
	{
		MarkerType face_on_face = g->mesh->CreateMarker();
		MarkerType edge_on_face = g->mesh->CreateMarker();
		MarkerType edge_on_edge = g->mesh->CreateMarker();
		MarkerType multi_edge   = g->mesh->CreateMarker();
		
		dynarray<Element,16> en1;
		tiny_map<int,int, 64> edges_mat;
		dynarray<Edge,128> face_edges;
		dynarray<Edge,128> potential_edges;
		dynarray<Edge,128> skipped_edges;
		dynarray<dynarray<Node,128> ,24> edge_cut_nodes(l);
		dynarray<dynarray<Node,128> ,24> edge_cut_nodes2(l);
		dynarray<Element,128> face_elements; // collect nodes and edges here in good order (normal outside)
		dynarray<Storage::integer,64> mat_intersection, mat_union, matse0,matse1;
		mat_ret_type mat0d,mat2d;
		dynarray<Node,24> node_in_face(l,InvalidNode());
		dynarray<HandleType,128> faces_on_face;
		dynarray<Face,128> inner_faces;
		dynarray<int,64> can_skip, cannot_skip, skip_mats;
		tiny_map<int,int,64> mat_cuts_on_edge[4];
		tiny_map<int, dynarray<Edge,128> ,64> edges_by_material; // edges that lay inside inital octree face
		ElementArray<Node> split_node(g->mesh,1);
		if( print ) std::cout << "calculate cutcell" << std::endl;

	
		
		for(i = 0; i < l; i++)
		{
			if( print ) std::cout << "face " << i << "/" << l << std::endl;
			
			
			
			bool is_create_center_node = false;
			face_elements.clear();
			edges_mat.clear();
			mat_cuts_on_edge[0].clear();
			mat_cuts_on_edge[1].clear();
			mat_cuts_on_edge[2].clear();
			mat_cuts_on_edge[3].clear();
			int edge_side = -1;
			int num_good_cuts = 0;


			for(j = 0; j < faces[i][0]; j++)
			{
				bool cut_edge = false;
				if( !mid[i][1+j] ) edge_side++;
				if( edge_side > 3 ) throw -1;
				Node n0 = verts[faces[i][1+j]];
				Node n1 = InvalidNode();
				Node n2 = verts[faces[i][1+(j+1)%faces[i][0]]];
				Storage::integer_array mat0 = n0->IntegerArray( g->materials );
				Storage::integer_array mat2 = n2->IntegerArray( g->materials );
				
				if( mat0.size() > 1 ) edge_cut_nodes2[i].push_back(n0);
				if( mat2.size() > 1 ) edge_cut_nodes2[i].push_back(n2);
				
				mat_intersection.resize(std::min(mat0.size(),mat2.size()));
				mat_intersection.resize(std::set_intersection(mat0.begin(),mat0.end(),mat2.begin(),mat2.end(),mat_intersection.begin())-mat_intersection.begin());
				
				if( mat_intersection.empty() )
				{
					cut_edge = true;
					mat0d.clear();
					mat0d.insert(mat0d.end(),mat0.begin(),mat0.end());
					mat2d.clear();
					mat2d.insert(mat2d.end(),mat2.begin(),mat2.end());
				}
				else
				{
					mat0d.resize(mat0.size());
					mat0d.resize(std::set_difference(mat0.begin(),mat0.end(),mat_intersection.begin(),mat_intersection.end(),mat0d.begin())-mat0d.begin());
					mat2d.resize(mat2.size());
					mat2d.resize(std::set_difference(mat2.begin(),mat2.end(),mat_intersection.begin(),mat_intersection.end(),mat2d.begin())-mat2d.begin());
					if( !mat0d.empty() && !mat2d.empty() )
					{
						cut_edge = true;
						//~ trigger_problem = 100;
					}
				}
				
				//~ if( mat0.size() == mat2.size() )
				//~ {
					//~ cut_edge = false;
					//~ for(k = 0; k < mat0.size(); k++)
						//~ if( mat0[k] != mat2[k] )
							//~ cut_edge = true;
				//~ }
				//~ else cut_edge = true;
				
				if( j == 0 ) face_elements.push_back(n0);
				
				if( cut_edge ) //cut edge!
				{
					if( print ) std::cout << "edge " << j << " is cut!" << std::endl;
					en1.clear();
					ElementArray<Element> nodesn0 = n0->BridgeAdjacencies(EDGE,NODE);
					ElementArray<Element> nodesn2 = n2->BridgeAdjacencies(EDGE,NODE);
					MarkerType inter = g->mesh->CreateMarker();
					
					if( inter == 0 ) throw -1;
					
					for(k = 0; k < nodesn0.size(); k++) nodesn0[k].SetMarker(inter);
					for(k = 0; k < nodesn2.size(); k++) if( nodesn2[k].GetMarker(inter) )	en1.push_back(nodesn2[k]);
					for(k = 0; k < nodesn0.size(); k++)	nodesn0[k].RemMarker(inter);
					g->mesh->ReleaseMarker(inter);
					
					
					if( print ) std::cout << "candidates: " << en1.size() << std::endl;

					{
						n1 = InvalidNode();
						double v1[3], v2[3];
						make_vec(&n0->Coords()[0],&n2->Coords()[0],v1);
						double l1 = sqrt(dot_prod(v1,v1)), l2;
						v1[0] /= l1;
						v1[1] /= l1;
						v1[2] /= l1;
						if( print ) std::cout << "l1: " << l1 << std::endl;
						for(k = 0; k < en1.size(); k++)
						{
							if( print ) std::cout << "cut edge candidate " << k << "/" << en1.size() << " " << en1[k]->LocalID() << std::endl; 
							make_vec(&en1[k]->getAsNode()->Coords()[0],&n2->Coords()[0],v2);
							l2 = sqrt(dot_prod(v2,v2));
							if( print ) std::cout << "l2: " << l2 << std::endl;
							if( l2 > 0 )
							{
								v2[0] /= l2;
								v2[1] /= l2;
								v2[2] /= l2;
								if( print )
								{
									double cp[3];
									cross_prod(v1,v2,cp);
									std::cout << "cross_product squared norm: " << dot_prod(cp,cp) << " dot product " << dot_prod(v1,v2) << " lengths " << l1 << " " << l2 << std::endl;
								}
								//if( dot_prod(cp,cp) < 1e-9 && l2 < l1)
								if( fabs(dot_prod(v1,v2)-1.0) < 1e-9 && l2 < l1 )
								{
									if( print ) std::cout << "choosen this node! " << std::endl;
									n1 = en1[k]->getAsNode();
									break;
								}
							}
						}	
					}
					
					if( n1 == InvalidNode() ) // node not found, create
					{
						if( print ) std::cout << "node not found - create one" << std::endl;
						Storage::real coord[3], coord2[3];
						matcenter centers[2];
						centers[0] = matcenter(mat0d,&n0->Coords()[0],0);
						centers[1] = matcenter(mat2d,&n2->Coords()[0],0);
						//make iterations to find correct position of the node
						int max = 16;
						int found = 0;
						for(int q = 0; q < max; q++)
						{
							//calculate mean center
							coord[0] = 0;
							coord[1] = 0;
							coord[2] = 0;
							for(int r = 0; r < 2; r++)
							{
								coord[0] += centers[r].get_center()[0];
								coord[1] += centers[r].get_center()[1];
								coord[2] += centers[r].get_center()[2];
							}
							coord[0] *= 0.5;
							coord[1] *= 0.5;
							coord[2] *= 0.5;
							
							mat_ret_type mats = g->get_material_types(coord);
							
							bool visited = false;
							for(int r = 0; r < 2; r++)
							{
								mat_ret_type inter = centers[r].intersect(mats);
								if( !inter.empty() )
								{
									coord2[0] = (coord[0] + centers[r].get_center()[0])*0.5;
									coord2[1] = (coord[1] + centers[r].get_center()[1])*0.5;
									coord2[2] = (coord[2] + centers[r].get_center()[2])*0.5;
									centers[r] = matcenter(inter,coord2, centers[r].get_count()+1);
									visited = true;
									found++;
								}
							}
							if( !visited ) 
								break;
						}
						
						if( found == 0 && !mat_intersection.empty()) 
						{
							cut_edge = false;
						}
						
						
						
						if( cut_edge )
						{
							
							mat_union.swap(mat_intersection);
							coord[0] = 0;
							coord[1] = 0;
							coord[2] = 0;
							for(int r = 0; r < 2; r++) 
							{
								coord[0] += centers[r].get_center()[0];
								coord[1] += centers[r].get_center()[1];
								coord[2] += centers[r].get_center()[2];
								
								{
									mat_intersection.resize(mat_union.size()+centers[r].get_mat().size());
									mat_intersection.resize(std::set_union(mat_union.begin(),mat_union.end(),centers[r].get_mat().begin(),centers[r].get_mat().end(),mat_intersection.begin())-mat_intersection.begin());
									mat_union.swap(mat_intersection);
								}
							}
							coord[0] *= 0.5;
							coord[1] *= 0.5;
							coord[2] *= 0.5;
							
							
							
							n1 = g->mesh->CreateNode(coord);
							//~ n1->SetMarker(g->cut_cell_node);
							Storage::integer_array mat1 = n1->IntegerArray(g->materials);
							
							
							mat1.resize(mat_union.size());
							for(int r = 0; r < mat_union.size(); ++r)
								mat1[r] = mat_union[r];

							if( print ) std::cout << "new node " << n1->LocalID() << std::endl;
						}
					}
					
					if( cut_edge )
					{
					
						Storage::integer_array mat1 = n1->IntegerArray(g->materials);
						Edge e0, e2;
						
						edge_cut_nodes[i].push_back(n1);
						num_good_cuts++;
						
						edge_nodes[0] = n0;
						edge_nodes[1] = n1;
						e0 = g->mesh->CreateEdge(edge_nodes).first;
						Storage::integer_array emat0 = e0->IntegerArray(g->materials);

						if( print ) std::cout << "new edge " << e0->LocalID() << std::endl;
						
						if( emat0.empty() )
						{
							mat_intersection.resize(std::min(mat0.size(),mat1.size()));
							mat_intersection.resize(std::set_intersection(mat0.begin(),mat0.end(),mat1.begin(),mat1.end(),mat_intersection.begin())-mat_intersection.begin());
							emat0.resize(mat_intersection.size());
							for(int r = 0; r < emat0.size(); r++) emat0[r] = mat_intersection[r];
						}

						e0->SetMarker(edge_on_edge);
						
						for(int r = 0; r < emat0.size(); r++)
						{
							mat_cuts_on_edge[edge_side][emat0[r]]++;
							edges_mat[emat0[r]] |= 1 << edge_side;
						}
						
						potential_edges.push_back(e0);
						
						e0->Integer(g->edge_face_number) = edge_side;
						
						
						
						face_elements.push_back(e0);
						
						face_elements.push_back(n1);
						
						
						edge_nodes[0] = n1;
						edge_nodes[1] = n2;
						e2 = g->mesh->CreateEdge(edge_nodes).first;
						Storage::integer_array emat2 = e2->IntegerArray(g->materials);

						if( print ) std::cout << "new edge " << e2->LocalID() << std::endl;
						
						if( emat2.empty() )
						{
							mat_intersection.resize(std::min(mat1.size(),mat2.size()));
							mat_intersection.resize(std::set_intersection(mat1.begin(),mat1.end(),mat2.begin(),mat2.end(),mat_intersection.begin())-mat_intersection.begin());
							emat2.resize(mat_intersection.size());
							for(int r = 0; r < emat2.size(); r++) emat2[r] = mat_intersection[r];
						}

						e2->SetMarker(edge_on_edge);
						
						
						for(int r = 0; r < emat2.size(); r++)
						{
							mat_cuts_on_edge[edge_side][emat2[r]]++;
							edges_mat[emat2[r]] |= 1 << edge_side;
						}
						potential_edges.push_back(e2);
							
						e2->Integer(g->edge_face_number) = edge_side;
							
						face_elements.push_back(e2);
						
						
							
						face_elements.push_back(n2);

					}
				}
				
				if( !cut_edge ) // regular edge
				{
					if( print ) std::cout << "edge " << j << " have no cut" << std::endl;
					Edge e0;
					edge_nodes[0] = n0;
					edge_nodes[1] = n2;
					e0 = g->mesh->CreateEdge(edge_nodes).first;
					Storage::integer_array emat0 = e0->IntegerArray(g->materials);
					
					if( emat0.empty() )
					{
						mat_intersection.resize(std::min(mat0.size(),mat2.size()));
						mat_intersection.resize(std::set_intersection(mat0.begin(),mat0.end(),mat2.begin(),mat2.end(),mat_intersection.begin())-mat_intersection.begin());
						emat0.insert(emat0.end(),mat_intersection.begin(),mat_intersection.end());
					}

					for(int r = 0; r < emat0.size(); r++)
						edges_mat[emat0[r]] |= 1 << edge_side;
					potential_edges.push_back(e0);
					e0->SetMarker(edge_on_edge);
					
					
					e0->Integer(g->edge_face_number) = edge_side;
					
					face_elements.push_back(e0);
					
					face_elements.push_back(n2);
				}
				
			}
			face_elements.pop_back(); //remove last node, that is equal to first
			
			skipped_edges.clear();
			face_edges.clear();
			
			//~ print = 1;
			//~ if( m == 636 ) print = 1;
			
			if( print ) for(k = 0; k < face_elements.size(); k++)
			{
				std::cout << k << " " << ElementTypeName(face_elements[k]->GetElementType()) << " " << face_elements[k]->LocalID() << " materials: ";
				Storage::integer_array mats = face_elements[k]->IntegerArray(g->materials);
				for(j = 0; j < mats.size(); j++) std::cout << mats[j] << " "; std::cout << std::endl;
			}
			
			{		
				cannot_skip.clear();
				can_skip.clear();
				
				for(dynarray<Element, 128>::iterator it = face_elements.begin(); it != face_elements.end(); ++it)
					if( (*it)->GetElementType() == EDGE )
					{
						Storage::integer_array a = (*it)->IntegerArray(g->materials);
						if( a.size() == 1 )
							cannot_skip.push_back(a[0]);
					}
				for(tiny_map<int,int,64>::iterator it = edges_mat.begin(); it != edges_mat.end() && !is_create_center_node; ++it) 
				{
					if( (it->second & 1) + ((it->second & 2) >> 1) + ((it->second & 4) >> 2) + ((it->second & 8) >> 3) <= 1 )
						can_skip.push_back(it->first);
					
				}
				std::sort(can_skip.begin(),can_skip.end());
				can_skip.resize(std::unique(can_skip.begin(),can_skip.end())-can_skip.begin());
				std::sort(cannot_skip.begin(),cannot_skip.end());
				cannot_skip.resize(std::unique(cannot_skip.begin(),cannot_skip.end())-cannot_skip.begin());
				skip_mats.resize(can_skip.size());
				skip_mats.resize(std::set_difference(can_skip.begin(),can_skip.end(),cannot_skip.begin(),cannot_skip.end(),skip_mats.begin())-skip_mats.begin());
				
				if( print )
				{
					std::cout << "cannot skip materials (" << cannot_skip.size() << "): ";
					for(j = 0; j < cannot_skip.size(); j++)
						std::cout << cannot_skip[j] << " ";
					std::cout << std::endl;
					std::cout << "can skip materials (" << can_skip.size() << "): ";
					for(j = 0; j < can_skip.size(); j++)
						std::cout << can_skip[j] << " ";
					std::cout << std::endl;
					std::cout << "skip materials (" << skip_mats.size() << "): ";
					for(j = 0; j < skip_mats.size(); j++)
						std::cout << skip_mats[j] << " ";
					std::cout << std::endl;
				}
			}
			{
				unsigned found = 0;
				Storage::integer_array mats0, mats1, matse0a, matse1a;
				tiny_map<int,int,5> edge_sides;
				bool stuck = false;
				while( !stuck )
				{
					tiny_map<int,int,64> nummats;
					tiny_map<int,int,64>::iterator minmat, qt;
					int local_found = 0, last_q = -1, last_e = -1, last_j = -1, cutnodes = 0;
					Storage::integer curmat;
					for(j = 1; j < face_elements.size(); j+=2) //traverse all edges
					{
						mats0 = face_elements[j]->IntegerArray(g->materials);
						for(Storage::integer_array::iterator kt = mats0.begin(); kt != mats0.end(); ++kt)
							if( !std::binary_search(skip_mats.begin(),skip_mats.end(),*kt) )
								nummats[*kt]++;
					}
					
					if( print )
					{
						std::cout << "nummats (" << nummats.size() << "): ";
						for(qt = nummats.begin(); qt != nummats.end(); qt++)
							std::cout << "(" << qt->first << "," << qt->second << ") ";
						std::cout << std::endl;
					}
					
					if( nummats.empty() ) break;
					
					bool success = false;
					while( !success )
					{
						minmat = nummats.begin();
						if( print ) std::cout << "select min: (" << minmat->first << "," << minmat->second << ")" << std::endl;
						for(qt = nummats.begin()+1; qt != nummats.end(); qt++)
							if( qt->second < minmat->second )
							{
								minmat = qt;
								if( print ) std::cout << "select min: (" << minmat->first << "," << minmat->second << ")" << std::endl;
							}
						curmat = minmat->first;
						
						if( print ) std::cout << "current material " << curmat << std::endl;
						
						dynarray< std::pair< std::pair<int,int> , int > , 64 > maxpath;
						found = 0;
						while( found < face_elements.size() )
						{
							int numedges = 0;
							local_found = 0;
							last_q = -1;
							last_e = -1;
							last_j = -1;
							cutnodes = 0;
							for(j = found; j < face_elements.size(); j+=2) //traverse nodes, find cuts on edges or nodes with multiple materials
							{
								mats0 = face_elements[j]->IntegerArray(g->materials);
								if( (!face_elements[j]->GetMarker(g->octree_node) || mats0.size() > 1) && std::binary_search(mats0.begin(),mats0.end(),curmat) ) //this is edge cut
								{
									found = j;
									matse0a = face_elements[found+1]->IntegerArray(g->materials);
									if( std::binary_search(matse0a.begin(),matse0a.end(),curmat) )
									{
										local_found = 1;
										matse0.resize(matse0a.size());
										matse0.resize(std::set_difference(matse0a.begin(),matse0a.end(),skip_mats.begin(),skip_mats.end(),matse0.begin())-matse0.begin());
										assert(!matse0.empty());
										if( !face_elements[found]->GetMarker(g->octree_node) ) cutnodes = 1;
										break;
									}
								}
							}
							if( !local_found ) 
							{
								if( print ) std::cout << "cut node for material not found" << std::endl;
								found = face_elements.size();
							}
							else
							{
								//scanf("%*c");
								if( print ) 
								{
									std::cout << "found node " << found << " materials: ";
									for(unsigned jjj = 0; jjj < mats0.size(); jjj++) std::cout << mats0[jjj] << " "; 
									std::cout << std::endl;
									std::cout << "found edge materials: ";
									for(unsigned jjj = 0; jjj < matse0.size(); jjj++) std::cout << matse0[jjj] << " "; 
									std::wcout << std::endl;
								}
								edge_sides.clear();
								for(j = 2; j < face_elements.size()+1; j+=2) // move node by node, skip found node
								{
									int e = (found+j-1)%face_elements.size(); // actual edge
									int q = (found+j)%face_elements.size(); // actual node
									mats1 = face_elements[q]->IntegerArray(g->materials);
									matse1a = face_elements[e]->IntegerArray(g->materials);
									matse1.resize(matse1a.size());
									matse1.resize(std::set_difference(matse1a.begin(),matse1a.end(),skip_mats.begin(),skip_mats.end(),matse1.begin())-matse1.begin());
									if( print )
									{
										std::cout << "next node " << q << " materials: ";
										for(unsigned jjj = 0; jjj < mats1.size(); jjj++) std::cout << mats1[jjj] << " "; std::cout << std::endl;
										
										std::cout << "prev edge " << e << " materials: ";
										for(unsigned jjj = 0; jjj < matse1.size(); jjj++) std::cout << matse1[jjj] << " "; std::cout << std::endl;
									}
									
									//enclosed edges should contain all the same materials
									//if they don't, we may cut out some other material
									bool ok = true;
									if( matse1.size() == matse0.size() )
									{
										
										for(unsigned q = 0; q < matse0.size(); q++)
											if( matse0[q] != matse1[q] )
												ok = false;
										
										
									}
									else ok = false;
									
									
									if( ok )
									{
										if( !face_elements[q]->GetMarker(g->octree_node) ) cutnodes++;
										//node may have less materials, but not any new material
										//~ mat_intersection.resize(mats1.size());
										//~ mat_intersection.resize(std::set_difference(mats1.begin(),mats1.end(),mats0.begin(),mats0.end(),mat_intersection.begin())-mat_intersection.begin());
										//~ if( mat_intersection.empty() ) //there is no other materials on this node
										{
											//check that this will not be degenerate face
											edge_sides[face_elements[e]->Integer(g->edge_face_number)]++;
											numedges++;
											//this node have all the same materials
											
											
											{
												if( print ) std::cout << "number of edge sides: " << edge_sides.size() << std::endl;
												if( edge_sides.size() > 1 )
												{
													if( print ) std::cout << "remember " << q << " " << e << " cutnodes " << cutnodes << std::endl;
													last_q = q;
													last_e = e;
													last_j = j;
													
												}
												
											}
											
											
											//we can't step through the node that have some other material
											//~ mat_intersection.resize(mats1.size());
											//~ mat_intersection.resize(std::set_difference(mats1.begin(),mats1.end(),matse0.begin(),matse0.end(),mat_intersection.begin())-mat_intersection.begin());
											//~ 
											//~ if( print )
											//~ {
												//~ std::cout << "difference of materials of first edge and current node: ";
												//~ if( mat_intersection.empty() )
													//~ std::cout << "is empty" << std::endl;
												//~ else
													//~ for(unsigned jjj = 0; jjj < mat_intersection.size(); jjj++) std::cout << mat_intersection[jjj] << " "; std::cout << std::endl;
											//~ }
											//~ 
											//~ 
											//~ if( !mat_intersection.empty() ) 
												//~ break;
										}
										//~ else 
										//~ {
											//~ if( print ) std::cout << " another material! " << std::endl;
											//~ break; //this node have another material - exit!
										//~ }
									} 
									else break;
								}
								if( last_q != -1 )
								{
									if( print ) std::cout << "loop " << found << "," << last_q << " have length " << numedges << " cutnodes " << cutnodes << std::endl; 
									//if( cutnodes > 0 || last_q == found) 
										maxpath.push_back( std::pair< std::pair<int,int >, int >( std::pair<int,int>(found,last_q) , numedges ) );
								}
								found += 2;
							}
						}
						if( !maxpath.empty() )
						{
							int maxpathj = 0;
							if( print ) std::cout << "select path length " << maxpath[maxpathj].second << std::endl;
							for(j = 1; j < maxpath.size(); j++)
								if( maxpath[j].second > maxpath[maxpathj].second )
								{
									maxpathj = j;
									if( print ) std::cout << "select path length " << maxpath[maxpathj].second << std::endl;
								}
							found = maxpath[maxpathj].first.first;
							last_q = maxpath[maxpathj].first.second;
							int numedges = maxpath[maxpathj].second;
							int q = last_q;
							if( print )
							{
								std::cout << "restored " << found << "," << last_q  << std::endl;
							}
							Edge new_edge = InvalidEdge();
							if( q != found )
							{
								edge_nodes[0] = face_elements[found]->getAsNode();
								edge_nodes[1] = face_elements[q]->getAsNode();
								if( face_elements.size()/2-numedges == 2 ) //check that we don't connect the edge over some node
								{
									Element middle_node = face_elements[(q+2)%face_elements.size()];
									Storage::real cnt0[3], cnt1[3], cntm[3], v1[3], v2[3], l1,l2;
									edge_nodes[0]->Centroid(cnt0);
									edge_nodes[1]->Centroid(cnt1);
									middle_node->Centroid(cntm);
									make_vec(cnt1,cnt0,v1);
									make_vec(cntm,cnt0,v2);
									l1 = sqrt(dot_prod(v1,v1));
									v1[0] /= l1;
									v1[1] /= l1;
									v1[2] /= l1;
									l2 = sqrt(dot_prod(v2,v2));
									v2[0] /= l2;
									v2[1] /= l2;
									v2[2] /= l2;
									if( print )
									{
										double cp[3];
										cross_prod(v1,v2,cp);									
										std::cout << "cross_product squared length: " << dot_prod(cp,cp) << " dot product: " << dot_prod(v1,v2) << " lengths " << l1 << " " << l2 << std::endl;
									}
									//if( dot_prod(cp,cp) < 1e-9 && l2 < l1) 
									if( fabs(dot_prod(v1,v2)-1.0) < 1e-9 && l2 < l1)
									{
										if( print ) std::cout << "prevent bad edge!" << std::endl;
										goto exit_work;
									}
									//check that they are not on one line
								}
								new_edge = g->mesh->CreateEdge(edge_nodes).first;
								//~ if( new_edge->LocalID() == 11801 ) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
								if( print ) std::cout << "connected! " << new_edge->LocalID() << std::endl;
								if( !new_edge->GetMarker(edge_on_edge) )
								{
									if( print ) std::cout << "new edge on face!" << std::endl;
									new_edge->SetMarker(edge_on_face);
									new_edge->Integer(g->edge_face_number) = sides[i]; //this algorithm shouldn't ever touch this new edge
									
									
									
									face_edges.push_back(new_edge);
									if( print ) std::cout << __FILE__ << ":" << __LINE__ << " put " << new_edge->GetHandle() << " to face_edges " << face_edges.size() << std::endl;
									
									
									Storage::integer_array mat = new_edge->IntegerArray(g->materials);
									if( mat.empty() )
									{
										Storage::integer_array maten0 = edge_nodes[0]->IntegerArray(g->materials);
										Storage::integer_array maten1 = edge_nodes[1]->IntegerArray(g->materials);
										mat_intersection.resize(maten0.size()+maten1.size());
										mat_intersection.resize(std::set_intersection(maten0.begin(),maten0.end(),maten1.begin(),maten1.end(),mat_intersection.begin())-mat_intersection.begin());
										mat.insert(mat.end(),mat_intersection.begin(),mat_intersection.end());
									}
									
									if( print )
									{
										std::cout << "materials: ";
										for(unsigned jj = 0; jj < mat.size(); jj++) std::cout << mat[jj] << " ";
										std::cout << std::endl;
									}
									
									for(unsigned qq = 0; qq < mat.size(); qq++)
										edges_by_material[mat[qq]].push_back(new_edge);
								}
								else 
								{
									if( print ) std::cout << "found existed edge on edge!" << std::endl;
									
									face_edges.push_back(new_edge); //the only variant this happen - we have found entire face, thus we need to add this edge, or it will be skipped
									if( print ) std::cout << __FILE__ << ":" << __LINE__ << " put " << new_edge->GetHandle() << " to face_edges " << face_edges.size() << std::endl;
								}
							}
							
							
							
								
							{
								//mark elements to skip
								mat_union.clear();
								MarkerType skip = g->mesh->CreateMarker();
								for(unsigned jj = 1; jj < face_elements.size()+1; jj++) //start from edge
								{
									int de = (found+jj)%face_elements.size(); // actual element
									if( de == q ) break;
									if( de % 2 == 0 ) //this is node
										face_elements[de]->SetMarker(skip);
									else
									{
										face_elements[de]->SetMarker(skip);
										skipped_edges.push_back(face_elements[de]->getAsEdge());
										mats1 = skipped_edges.back()->IntegerArray(g->materials);
										mat_intersection.resize(mats1.size()+mat_union.size());
										mat_intersection.resize(std::set_union(mats1.begin(),mats1.end(),mat_union.begin(),mat_union.end(),mat_intersection.begin())-mat_intersection.begin());
										mat_union.swap(mat_intersection);
									}
									if( print ) std::cout << "marked to skip: " << de << std::endl;
								}
								
								for(unsigned jj = 1; jj < face_elements.size(); jj+=2)
								{
									if( !face_elements[jj]->GetMarker(skip) )
									{
										mats1 = face_elements[jj]->IntegerArray(g->materials);
										mat_intersection.resize(mat_union.size());
										mat_intersection.resize(std::set_difference(mat_union.begin(),mat_union.end(),mats1.begin(),mats1.end(),mat_intersection.begin())-mat_intersection.begin());
										mat_union.swap(mat_intersection);
									}
								}
								
								//remove materials
								if( print )
								{
									std::cout << "erase materials:" << std::endl;
									for(unsigned jjj = 0; jjj < mat_union.size(); jjj++) std::cout << mat_union[jjj] << " "; std::cout << std::endl;
								}
								if(!mat_union.empty())
								{
									mat_intersection.resize(mat_union.size()+skip_mats.size());
									mat_intersection.resize(std::set_union(mat_union.begin(),mat_union.end(),skip_mats.begin(),skip_mats.end(),mat_intersection.begin())-mat_intersection.begin());
									skip_mats.swap(mat_intersection);
									if( print )
									{
										std::cout << "skip materials (" << skip_mats.size() << "): ";
										for(j = 0; j < skip_mats.size(); j++)
											std::cout << skip_mats[j] << " ";
										std::cout << std::endl;
									}
									//~ for(unsigned jjj = 0; jjj < mat_union.size(); jjj++)
										//~ for(tiny_map<int,int,64>::iterator it = edges_mat.begin(); it != edges_mat.end(); ++it) 
											//~ if( it->first == mat_union[jjj] )
											//~ {
												//~ if( print ) std::cout << "erase material: " << it->first << std::endl;
												//~ edges_mat.erase(it);
												//~ break;
											//~ }
									//~ if( print ) std::cout << " number of materials left: " << edges_mat.size() << std::endl;
								}
								
								{ //find new face_elements array
									dynarray<Element,128> replace;
									int started_good = 0, pushed = 0;
									for(unsigned jj = 0; jj < face_elements.size(); jj++)
									{
										if( face_elements[jj]->GetMarker(skip) )
										{
											if( started_good != pushed )
											{
												pushed++;
												found = replace.size()-1;
												replace.push_back(new_edge);
												if( print ) std::cout << "pushed new: " << (!new_edge.isValid() ? -1 : new_edge->LocalID()) << std::endl;
											}
											
											if( print ) std::cout << "skipped: " << face_elements[jj]->LocalID() << std::endl;
										}
										else
										{
											if( started_good == 0 ) started_good++;
											replace.push_back(face_elements[jj]);
											if( print ) std::cout << "pushed: " << face_elements[jj]->LocalID() << std::endl;
										}
										face_elements[jj]->RemMarker(skip);
									}
									
									
									if( !pushed ) 
									{
										found = replace.size()-1;
										replace.push_back(new_edge);
									}
									
									if( replace.size() == 4 && replace[1] == replace[3]  || replace.size() < 4) 
									{
										//if( new_edge != NULL ) face_edges.push_back(new_edge);
										replace.clear(); // there is only one node and one edge left - avoid this edge duplication
										stuck = true;
									}
									
									if( print ) std::cout << "replace size: " << replace.size() << std::endl;
									
									face_elements.swap(replace);
									
									if( print ) std::cout << "face_elements size: " << face_elements.size() << std::endl;
									
									if( print ) 
									{
										std::cout << "new array: " << std::endl;
										for(unsigned kk = 0; kk < face_elements.size(); kk++)
										{
											std::cout << kk << " " << ElementTypeName(face_elements[kk]->GetElementType()) << " " << face_elements[kk]->LocalID() << " materials: ";
											Storage::integer_array mats = face_elements[kk]->IntegerArray(g->materials);
											for(unsigned kkk = 0; kkk < mats.size(); kkk++) std::cout << mats[kkk] << " "; std::cout << std::endl;
										}
									}
									if( print ) std::cout << "done!" << std::endl;
								}
								g->mesh->ReleaseMarker(skip);
							}

							success = true;
						}
						
						
exit_work:				if( !success ) 
						{
							if( print ) std::cout << "material " << curmat << " not successful" << std::endl;
							nummats.erase(minmat);
							if( nummats.empty() ) 
							{
								if( print ) std::cout << "shutdown" << std::endl;
								stuck = true;
								break;
							}
						}
					}
					
					
				}
			}
			
			
			//~ if( m == 636 ) print = 0;
			
			if( print ) 
			{
				std::cout << "array: " << std::endl;
				for(unsigned kk = 0; kk < face_elements.size(); kk++)
				{
					std::cout << kk << " " << face_elements[kk]->GetHandle() << " " << ElementTypeName(face_elements[kk]->GetElementType()) << " " << face_elements[kk]->LocalID() << " materials: ";
					Storage::integer_array mats = face_elements[kk]->IntegerArray(g->materials);
					for(unsigned kkk = 0; kkk < mats.size(); kkk++) std::cout << mats[kkk] << " "; std::cout << std::endl;
				}
			}
	
			if( !face_elements.empty() )
			{
				Storage::integer_array mats = face_elements[1]->IntegerArray(g->materials);
				mat_intersection.clear();
				mat_intersection.insert(mat_intersection.end(),mats.begin(),mats.end());
				mat_union.resize(mat_intersection.size());
				for(j = 3; j < face_elements.size(); j+=2)
				{
					mats = face_elements[j]->IntegerArray(g->materials);
					mat_union.resize(std::set_intersection(mat_intersection.begin(),mat_intersection.end(),mats.begin(),mats.end(),mat_union.begin())-mat_union.begin());
					mat_intersection.swap(mat_union);
					if( mat_intersection.empty() ) break;
				}
				
				
				is_create_center_node = mat_intersection.empty();
				
				
				if( is_create_center_node )
				{
					//try to find the node
					ElementArray<Element> result = face_elements[0]->BridgeAdjacencies(EDGE,NODE);
					for(k = 2; k < face_elements.size(); k+=2) //all edges
						result.Intersect(face_elements[k]->BridgeAdjacencies(EDGE,NODE));
					
					
					if( !result.empty() ) 
					{
						for(k = 0; k < result.size(); k++)
							if( result[k].Integer(g->face_center_node) == 1 ) //(sides[i]/3+1) )
								node_in_face[i] = result[k].getAsNode();
					}
					if( !node_in_face[i].isValid() )
					{
						
						Storage::real coord[3],coord2[3];
						dynarray<matcenter,32> centers;
						for(j = 0; j < face_elements.size(); j+=2)
							centers.push_back(matcenter(face_elements[j]->IntegerArray(g->materials),&face_elements[j]->getAsNode()->Coords()[0]));
						int max = centers.size();
						Storage::real div = 1.0/(Storage::real)centers.size();
						for(int q = 0; q < max; q++)
						{
							//calculate mean center
							coord[0] = 0;
							coord[1] = 0;
							coord[2] = 0;
							for(j = 0; j < centers.size(); j++)
							{
								coord[0] += centers[j].get_center()[0];
								coord[1] += centers[j].get_center()[1];
								coord[2] += centers[j].get_center()[2];
							}
							coord[0] *= div;
							coord[1] *= div;
							coord[2] *= div;
							
							
							
							mat_ret_type mats = g->get_material_types(coord);
								
							bool visited = false;
							for(j = 0; j < centers.size(); j++)
							{
								mat_ret_type inter = centers[j].intersect(mats);
								if( !inter.empty() )
								{
									coord2[0] = (coord[0] + centers[j].get_center()[0])*0.5;
									coord2[1] = (coord[1] + centers[j].get_center()[1])*0.5;
									coord2[2] = (coord[2] + centers[j].get_center()[2])*0.5;
									centers[j] = matcenter(centers[j].unite(mats),coord2, centers[j].get_count()+1);
									visited = true;
								}
							}
							if( !visited ) break;
							
						
						}
						mat_union.clear();
						coord[0] = 0;
						coord[1] = 0;
						coord[2] = 0;
						div = 0;
						for(j = 0; j < centers.size(); j++)
						{
							div++;
							coord[0] += centers[j].get_center()[0];
							coord[1] += centers[j].get_center()[1];
							coord[2] += centers[j].get_center()[2];
							mat_intersection.resize(mat_union.size()+centers[j].get_mat().size());
							mat_intersection.resize(std::set_union(mat_union.begin(),mat_union.end(),centers[j].get_mat().begin(),centers[j].get_mat().end(),mat_intersection.begin())-mat_intersection.begin());
							mat_union.swap(mat_intersection);
						}
						div = 1.0/div;
						coord[0] *= div;
						coord[1] *= div;
						coord[2] *= div;
						
						
						node_in_face[i] = g->mesh->CreateNode(coord);
						Storage::integer_array mat = node_in_face[i]->IntegerArray(g->materials);
						mat.insert(mat.end(),mat_union.begin(),mat_union.end());
					}
					
					node_in_face[i]->Integer(g->face_center_node) = 1;//(sides[i]/3+1);
					
					if( print ) std::cout << "created face centered node!" << std::endl;
					
					MarkerType del = g->mesh->CreateMarker();
					for(j = 1; j < face_elements.size(); j+=2) if( !face_elements[j]->GetMarker(edge_on_edge) ) face_elements[j]->SetMarker(del);

					dynarray<Edge,128>::iterator it = face_edges.begin();
					while(it != face_edges.end())
						if((*it)->GetMarker(del) ) 
						{
							if( print ) std::cout << __FILE__ << ":" << __LINE__ << " erase " << it->GetHandle();
							it = face_edges.erase(it); 
							if( print ) std::cout << " to face_edges " << face_edges.size() << std::endl;
						}
						else it++;
						
					
					for(tiny_map<int, dynarray<Edge,128> ,64>::iterator jt = edges_by_material.begin();
						jt != edges_by_material.end(); ++jt)
					{
						it = jt->second.begin();
						while(it != jt->second.end())
							if((*it)->GetMarker(del) ) 
								it = jt->second.erase(it); 
							else it++;
					}


					for(j = 1; j < face_elements.size(); j+=2) face_elements[j]->RemMarker(del);
					g->mesh->ReleaseMarker(del);
					
					split_node[0] = node_in_face[i];
					edge_nodes[0] = node_in_face[i];
					for(j = 1; j < face_elements.size(); j+=2)
					{
						if( !face_elements[j]->GetMarker(edge_on_edge) )
						{
							if( print ) std::cout << "EDGE " << face_elements[j]->LocalID() << " is edge on edge, split" << std::endl;
							ElementArray<Edge> result = Edge::SplitEdge(face_elements[j]->getAsEdge(),split_node,0);

							for(k = 0; k < result.size(); k++)
							{
						
								result[k]->SetMarker(edge_on_face);
								result[k]->Integer(g->edge_face_number) = sides[i];
						
								
								face_edges.push_back(result[k]);
								if( print ) std::cout << __FILE__ << ":" << __LINE__ << " put " << result[k]->LocalID() << " to face_edges " << face_edges.size() << std::endl;
						
								Storage::integer_array mat = result[k]->IntegerArray(g->materials);
						
								if( mat.empty() )
								{
									ElementArray<Node> nodes = result[k]->getNodes();
									Storage::integer_array mat0 = nodes[0].IntegerArray(g->materials);
									Storage::integer_array mat1 = nodes[1].IntegerArray(g->materials);
									mat.resize(std::min(mat0.size(),mat1.size()));
									mat.resize(std::set_intersection(mat0.begin(),mat0.end(),mat1.begin(),mat1.end(),mat.begin())-mat.begin());
								}
						
								for(unsigned q = 0; q < mat.size(); q++)
									edges_by_material[mat[q]].push_back(result[k]);
							}
						}
						else
						{
							if( print ) std::cout << "EDGE " << face_elements[j]->LocalID() << " is edge on face, connect" << std::endl;
							for(k = 0; k < 2; k++)
							{
								edge_nodes[1] = face_elements[(j-1+k*2)%face_elements.size()]->getAsNode();
								Edge new_edge = g->mesh->CreateEdge(edge_nodes).first;
								new_edge->SetMarker(edge_on_face);
								new_edge->Integer(g->edge_face_number) = sides[i];
								
								
								face_edges.push_back(new_edge);
								if( print ) std::cout << __FILE__ << ":" << __LINE__ << " put " << new_edge->LocalID() << " to face_edges " << face_edges.size() << std::endl;
								Storage::integer_array mat = new_edge->IntegerArray(g->materials);
						
								if( mat.empty() )
								{
									Storage::integer_array mat0 = edge_nodes[0]->IntegerArray(g->materials);
									Storage::integer_array mat1 = edge_nodes[1]->IntegerArray(g->materials);
									mat.resize(std::min(mat0.size(),mat1.size()));
									mat.resize(std::set_intersection(mat0.begin(),mat0.end(),mat1.begin(),mat1.end(),mat.begin())-mat.begin());
								}
						
								for(unsigned q = 0; q < mat.size(); q++)
									edges_by_material[mat[q]].push_back(new_edge);
							}
							face_edges.push_back(face_elements[j]->getAsEdge());
							if( print ) std::cout << __FILE__ << ":" << __LINE__ << " put " << face_elements[j]->LocalID() << " to face_edges " << face_edges.size() << std::endl;
						}
					}
					face_elements.clear();
					std::sort(face_edges.begin(),face_edges.end());
					face_edges.resize(std::unique(face_edges.begin(),face_edges.end())-face_edges.begin());
					//try to adjust the position of the face center node
					/*
					{
						dynarray<Element *,64> other_node;
						for(j = 0; j < face_edges.size(); j++)
						{
							adjacent<Element> n = face_edges[j]->getAdjElements(NODE);
							if( &n[0] == node_in_face[i] ) other_node.push_back(&n[1]);
							else other_node.push_back(&n[0]);
						}
						
					}
					*/
				}
			}
			
			if( print ) for(j = 0; j < face_edges.size(); j++) std::cout << "face_edges [" << j << "]: " << face_edges[j]->LocalID() << std::endl;
			
			if( print ) 
			{
				std::cout << "array: " << std::endl;
				for(unsigned kk = 0; kk < face_elements.size(); kk++)
				{
					std::cout << kk << " " << face_elements[kk]->GetHandle() << " " << ElementTypeName(face_elements[kk]->GetElementType()) << " " << face_elements[kk]->LocalID() << " materials: ";
					Storage::integer_array mats = face_elements[kk]->IntegerArray(g->materials);
					for(unsigned kkk = 0; kkk < mats.size(); kkk++) std::cout << mats[kkk] << " "; std::cout << std::endl;
				}
			}
			
			for(j = 1; j < face_elements.size(); j+=2)
			{
				if( !face_elements[j]->GetMarker(edge_on_face) )
				{
					
					face_edges.push_back(face_elements[j]->getAsEdge());
					if( print ) std::cout << __FILE__ << ":" << __LINE__ << " put " << face_elements[j]->GetHandle() << " to face_edges " << face_edges.size() << std::endl;
				}
			}

			if( print ) 
			{
				std::cout << "on start" << std::endl;
				for(j = 0; j < face_edges.size(); j++) std::cout << "face_edges [" << j << "]: " << face_edges[j]->LocalID() << std::endl;
			}
				
			face_edges.insert(face_edges.end(),skipped_edges.begin(),skipped_edges.end());
			
			if( print ) 
			{
				std::cout << "insert skipped" << std::endl;
				for(j = 0; j < face_edges.size(); j++) std::cout << "face_edges [" << j << "]: " << face_edges[j]->LocalID() << std::endl;
			}
			
			std::sort(face_edges.begin(),face_edges.end());
			face_edges.resize(std::unique(face_edges.begin(),face_edges.end())-face_edges.begin());
			if( print ) 
			{
				std::cout << "remove unique" << std::endl;
				for(j = 0; j < face_edges.size(); j++) std::cout << "face_edges [" << j << "]: " << face_edges[j]->LocalID() << std::endl;
			}
			std::sort(face_edges.begin(),face_edges.end(),edge_Comparator(edge_on_edge));
			
			int num_inner = 0;
			for(k = 0; k < face_edges.size(); k++)
				if( !face_edges[k]->GetMarker(edge_on_edge) )
					num_inner++;
				else break;
			
			
			if( print )
			{
				std::cout << "num inner: " << num_inner << std::endl;
				for(k = 0; k < face_edges.size(); k++)
				{
					std::cout << "edge " << k << " id " << face_edges[k]->LocalID();
					if( face_edges[k]->GetMarker(edge_on_face) ) std::cout << " (inner)";
					if( face_edges[k]->GetMarker(edge_on_edge) )std::cout << " (outer)";
					std::cout << " material:";
					Storage::integer_array mat = face_edges[k]->IntegerArray(g->materials);
					for(j = 0; j < mat.size(); j++) std::cout << " " << mat[j];
					std::cout << " nodes: " << face_edges[k]->getNodes()[0].LocalID() << " , " << face_edges[k]->getNodes()[1].LocalID();
					std::cout << std::endl;
				}
			}
				
			incident_matrix<Edge> matrix(face_edges.begin(),face_edges.end(),num_inner);
			
			
			ElementArray<Edge> loop(g->mesh);
			while(matrix.find_shortest_loop(loop))
			{ 
				if( print ) std::cout << "found loop " << loop.size() << std::endl;
				if( loop.size() > 2 )
				{
					if( print )
						for(k = 0; k < loop.size(); k++)
						{
							std::cout << k << " " << loop[k]->LocalID() << " material";
							Storage::integer_array mat = loop[k]->IntegerArray(g->materials);
							for(j = 0; j < mat.size(); j++) std::cout << " " << mat[j];
							for(j = 0; j < face_edges.size(); j++)
								if( loop[k] == face_edges[j] ) 
								{
									std::cout << " edge " << j;
									break;
								}
							std::cout << std::endl;
						}

						Face new_face = g->mesh->CreateFace(loop).first;
					
					
					
					
					
					new_face->SetMarker(face_on_face);
						
					Storage::integer_array mat = new_face->IntegerArray(g->materials);
						

					if( mat.empty() )
					{
						Storage::integer_array mt = loop[0]->IntegerArray(g->materials);
						mat_intersection.clear();
						mat_intersection.insert(mat_intersection.end(),mt.begin(),mt.end());
						mat_union.resize(mat_intersection.size());
						for(k = 1; k < loop.size(); k++)
						{
							mt = loop[k]->IntegerArray(g->materials);
							mat_union.resize(std::set_intersection(mat_intersection.begin(),mat_intersection.end(),mt.begin(),mt.end(),mat_union.begin())-mat_union.begin());
							mat_intersection.swap(mat_union);
						}
						mat.insert(mat.end(),mat_intersection.begin(),mat_intersection.end());
					}
					
					
						
					for(unsigned q = 0; q < mat.size(); q++)
					{
						if( print ) 
						{
							std::cout << __FILE__ << ":" << __LINE__ << " "  << mat[q] << " added face " << new_face->LocalID() << " neighbours " << new_face->nbAdjElements(CELL) << " ";
							if( new_face->nbAdjElements(CELL) == 2 )
								std::cout << new_face->BackCell()->GetHandle() << " " <<  new_face->BackCell()->LocalID() << " " << new_face->FrontCell()->GetHandle() << " " << new_face->FrontCell()->LocalID() << std::endl;
							else std::cout << std::endl;
						}
						//~ faces_by_material[mat[q]].push_back(new_face);
						
					}
					faces_on_face.push_back(new_face->GetHandle());
				}
				else if( !loop.empty() ) std::cout << __FILE__ << ":" << __LINE__ << " skip degenerate face " << loop.size() << std::endl;
			}
			//~ if( !matrix.all_visited() ) trigger_problem = 1;
		}
		
		// Create one cell with many faces, each face for each material
		
		if( false )
		{
			dynarray<HandleType,128> faces_on_face_copy(faces_on_face);
			std::sort(faces_on_face_copy.begin(),faces_on_face_copy.end());
			int old_size = faces_on_face_copy.size();
			faces_on_face_copy.resize(std::unique(faces_on_face_copy.begin(),faces_on_face_copy.end())-faces_on_face_copy.begin());

			ElementArray<Face> c_faces(g->mesh,faces_on_face.size());
			for(i = 0; i < faces_on_face.size(); ++i) c_faces.at(i) = faces_on_face[i];
			//for(std::map<int,std::vector<Face *> >::iterator it = faces_by_material.begin(); it != faces_by_material.end(); ++it)
			//	c_faces.insert(c_faces.end(),it->second.begin(),it->second.end());
			Cell c = g->mesh->CreateCell(c_faces).first;
			if( print ) std::cout << __FILE__ << ":" << __LINE__ << " new cell " << c->GetHandle() << " id " << c->LocalID() << std::endl;
			(*g->cells[m].mr).push_back(c->GetHandle());
			if( c->GetGeometricType() == Element::MultiPolygon ) 
			{
				
				std::cout << "Problem" << std::endl;
				std::cout << "faces " << l << std::endl;
				std::map<HandleType, int> edge_visit;
				for(i = 0; i < c_faces.size(); i++)
				{
					ElementArray<Edge> edges = c_faces[i]->getEdges();
					for(ElementArray<Edge>::iterator it = edges.begin(); it != edges.end(); it++)
						edge_visit[it->GetHandle()]++;
				}
				for(i = 0; i < l; i++)
				{
					std::cout << " face " << i << std::endl;
					std::cout << " original face " << std::endl;
					for(j = 0; j < faces[i][0]; j++)
					{
						Storage::real_array c = verts[faces[i][1+j]]->Coords();
						std::cout << " (" << c[0] << "," << c[1] << "," << c[2] << ")";
					}
					std::cout << std::endl;
				}
				
				std::cout << "edge_cut_nodes " << std::endl;
				
				for(i = 0; i < l; i++)
				{
					std::cout << " face " << i << std::endl;
					for(j = 0; j < edge_cut_nodes[i].size(); j++)
					{
						Storage::real_array c = edge_cut_nodes[i][j]->Coords();
						std::cout << j << " " << edge_cut_nodes[i][j]->GetHandle() << " ( " << c[0] << "," << c[1] << "," << c[2] << ") ";
						std::cout << std::endl;
					}
					
				}
				
				
				//(*g->cells[m].mr).back()->Integer(g->problem) = 1;

				//~ throw -1;
			}
			Cell(g->mesh,(*g->cells[m].mr).back())->Integer(g->cell_material) = c_faces[0]->Integer(g->materials);
		}
		else
		
		{
			
			// create a node inside cell if there are more then two nodes on faces
			
			Node cell_center_node = InvalidNode();
			
			int num_face_center_nodes = 0;
			for(j = 0; j < l; j++)
				if( node_in_face[j].isValid() )
					num_face_center_nodes++;
			
			//if( m == 2106 ) print = 1;		
			
			if( print ) std::cout << "num_face_center_nodes: " << num_face_center_nodes << std::endl;
			
			
			
			if( num_face_center_nodes > 2 )
			{
				if( print ) std::cout << "added center node" << std::endl;
				Storage::real coords[3] = {0,0,0};
				//CALCULATE CORRECT LOCATION HERE!
				//for the cell-centered node
				/*
				for(j = 0; j < 8; j++)
				{
					Storage::real_array c = g->verts[g->cells[m].vertexes[j]].mv->Coords();
					coords[0] += c[0];
					coords[1] += c[1];
					coords[2] += c[2];
				}
				coords[0] *= 0.125;
				coords[1] *= 0.125;
				coords[2] *= 0.125;
				*/
				
				Storage::real coord[3],coord2[3];
				dynarray<matcenter,24> centers;
				for(j = 0; j < l; j++) if( node_in_face[j].isValid() )
					centers.push_back(matcenter(node_in_face[j]->IntegerArray(g->materials),&node_in_face[j]->Coords()[0]));
				//make iterations to find correct position of the node
				int max = 8;
				Storage::real div = 1.0/centers.size();
				for(int q = 0; q < max; q++)
				{
					//calculate mean center
					coord[0] = 0;
					coord[1] = 0;
					coord[2] = 0;
					for(j = 0; j < centers.size(); j++)
					{
						coord[0] += centers[j].get_center()[0];
						coord[1] += centers[j].get_center()[1];
						coord[2] += centers[j].get_center()[2];
					}
					coord[0] *= div;
					coord[1] *= div;
					coord[2] *= div;

					mat_ret_type mats = g->get_material_types(coord);
							
					bool visited = false;
					for(j = 0; j < centers.size(); j++)
					{
						mat_ret_type inter = centers[j].intersect(mats);
						if( !inter.empty() )
						{
							coord2[0] = (coord[0] + centers[j].get_center()[0])*0.5;
							coord2[1] = (coord[1] + centers[j].get_center()[1])*0.5;
							coord2[2] = (coord[2] + centers[j].get_center()[2])*0.5;
							centers[j] = matcenter(inter,coord2, centers[j].get_count()+1);
							visited = true;
						}
					}
					if( !visited ) break;
				}
				//calculate mean center the last time
				coord[0] = 0;
				coord[1] = 0;
				coord[2] = 0;
				for(j = 0; j < centers.size(); j++)
				{
					coord[0] += centers[j].get_center()[0];
					coord[1] += centers[j].get_center()[1];
					coord[2] += centers[j].get_center()[2];
				}
				coord[0] *= div;
				coord[1] *= div;
				coord[2] *= div;
				
				cell_center_node = g->mesh->CreateNode(coord);
				
				//~ cell_center_node->SetMarker(g->cut_cell_node);
				Storage::integer_array mat = cell_center_node->IntegerArray(g->materials);
				mat_union.clear();
				for(j = 0; j < l; j++)
				{
					if( node_in_face[j].isValid() )
					{
						Storage::integer_array matj = node_in_face[j]->IntegerArray(g->materials);
						mat_intersection.resize(mat_union.size()+matj.size());
						mat_intersection.resize(std::set_union(mat_union.begin(),mat_union.end(),matj.begin(),matj.end(),mat_intersection.begin())-mat_intersection.begin());
						
						mat_union.swap(mat_intersection);
						
						
						//now connect nodes in face centers to cell center
						
						edge_nodes[0] = node_in_face[j];
						edge_nodes[1] = cell_center_node;
						Edge new_edge = g->mesh->CreateEdge(edge_nodes).first;
						//~ if( new_edge->LocalID() == 11801 ) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
						//~ new_edge->SetMarker(hyper_edge);
						Storage::integer_array mate = new_edge->IntegerArray(g->materials);
						if( mate.empty() ) mate.insert(mate.end(),matj.begin(),matj.end());
						
						//add obtained edges to array material by material
						for(k = 0; k < mate.size(); k++)
						{
							edges_by_material[mate[k]].push_back(new_edge);
							if( print ) std::cout << __FILE__ << ":" << __LINE__ << " " << mate[k] << " added edge " << new_edge->GetHandle() << std::endl;
						}
					}
				}
				mat.insert(mat.end(),mat_union.begin(),mat_union.end());
				
				//now connect central node with cuts on edges
				MarkerType mrk = g->mesh->CreateMarker();
				for(j = 0; j < l; j++)
				{
					//edge_cut_nodes[j].insert(edge_cut_nodes[j].end(),edge_cut_nodes2[j].begin(),edge_cut_nodes2[j].end());// do we need to connect degenerate cuts?
					for(k = 0; k < edge_cut_nodes[j].size(); k++) 
					if( !edge_cut_nodes[j][k]->GetMarker(mrk) )
					{
						edge_nodes[0] = cell_center_node;
						edge_nodes[1] = edge_cut_nodes[j][k];
						
						Edge new_edge = g->mesh->CreateEdge(edge_nodes).first;
						//~ if( new_edge->LocalID() == 11801 ) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
						//~ new_edge->SetMarker(hyper_edge);
						Storage::integer_array matjk = edge_cut_nodes[j][k]->IntegerArray(g->materials);
						Storage::integer_array mate = new_edge->IntegerArray(g->materials);
						if( mate.empty() ) 
						{
							mat_intersection.resize(mat.size()+matjk.size());
							mat_intersection.resize(std::set_intersection(mat.begin(),mat.end(),matjk.begin(),matjk.end(),mat_intersection.begin())-mat_intersection.begin());
							mate.insert(mate.end(),mat_intersection.begin(),mat_intersection.end());
						
							//add obtained edges to array material by material
							for(int q = 0; q < mate.size(); q++)
							{
								edges_by_material[mate[q]].push_back(new_edge);
								if( print ) std::cout << __FILE__ << ":" << __LINE__ << " " << mate[q] << " added edge " << new_edge->GetHandle() << std::endl;
							}
						}
						edge_cut_nodes[j][k]->SetMarker(mrk);
					}
				}	
				for(j = 0; j < l; j++)
					for(k = 0; k < edge_cut_nodes[j].size(); k++)
						edge_cut_nodes[j][k]->RemMarker(mrk);
				g->mesh->ReleaseMarker(mrk);
				
			}
			else if( num_face_center_nodes == 2 )
			{
				
				int eside[2] = {-1,-1};
				Storage::real c1[2], c2;
				Storage::integer_array mat1[2],mat2;
				if( print ) std::cout << "connect center nodes" << std::endl;
				//only one case may happen when nodes are on opposite site
				//just connect them
				k = 0;
				for(j = 0; j < l; j++)
					if( node_in_face[j].isValid() )
					{
						eside[k] = sides[j];
						c1[k] = node_in_face[j]->Coords()[eside[k]/2];
						mat1[k] = node_in_face[j]->IntegerArray(g->materials);
						edge_nodes[k++] = node_in_face[j];
						if( k == 2 ) break;
					}
						
				Storage::integer_array mat0 = edge_nodes[0]->IntegerArray(g->materials);
				
				
				
				Edge new_edge = g->mesh->CreateEdge(edge_nodes).first;
				//~ if( new_edge->LocalID() == 11801 ) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
				//~ new_edge->SetMarker(hyper_edge);
				Storage::integer_array mate = new_edge->IntegerArray(g->materials);
				if( mate.empty() ) mate.insert(mate.begin(),mat0.begin(),mat0.end());
				
				
				//add obtained edges to array material by material
				for(k = 0; k < mate.size(); k++)
					edges_by_material[mate[k]].push_back(new_edge);
					
					
				//Node * remember[2] = {edge_nodes[0],edge_nodes[1]};
				//MarkerType mrk = g->mesh->CreateMarker();
				/*
				for(j = 0; j < l; j++)
				{
					edge_cut_nodes[j].insert(edge_cut_nodes[j].end(),edge_cut_nodes2[j].begin(),edge_cut_nodes2[j].end());// do we need to connect degenerate cuts?
				}
				*/
				/*
				for(unsigned qq = 0; qq < 2; qq++)
				{
					edge_nodes[0] = remember[qq];
					for(j = 0; j < l; j++)
					{
						for(k = 0; k < edge_cut_nodes[j].size(); k++) if( !edge_cut_nodes[j][k]->GetMarker(mrk) )
						{
							c2 = edge_cut_nodes[j][k]->Coords()[eside[qq]/2];
							if( fabs(c1[qq]-c2) > 1e-5 ) //not on one plain
							if( eside[qq] != sides[j] )
							{
								mat2 = edge_cut_nodes[j][k]->IntegerArray(g->materials);
								bool connect = true;
								
								mat_intersection.resize(mat2.size()+mat1[qq].size());
								mat_intersection.resize(std::set_difference(mat2.begin(),mat2.end(),mat1[qq].begin(),mat1[qq].end(),mat_intersection.begin())-mat_intersection.begin());
								if( mat_intersection.empty() ) connect = true;
								
								//if( mat1[qq].size() == mat2.size() )
								//{
								//	for(unsigned qqq = 0; qqq < mat1[qq].size(); qqq++)
								//		if( mat1[qq][qqq] != mat2[qqq] )
								//			connect = false;
								//}
								//else connect = false;
								
								if( connect )
								{
									
									edge_nodes[1] = edge_cut_nodes[j][k];
									
									Edge * new_edge = g->mesh->CreateEdge(edge_nodes);
									//~ if( new_edge->LocalID() == 11801 ) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
									Storage::integer_array mate = new_edge->IntegerArray(g->materials);
									if( mate.empty() ) 
									{
										mat_intersection.resize(mat1[qq].size()+mat2.size());
										mat_intersection.resize(std::set_intersection(mat2.begin(),mat2.end(),mat1[qq].begin(),mat1[qq].end(),mat_intersection.begin())-mat_intersection.begin());
										mate.insert(mate.begin(),mat_intersection.begin(),mat_intersection.end());
										for(unsigned q = 0; q < mate.size(); q++)
											edges_by_material[mate[q]].push_back(new_edge);
									}
									edge_cut_nodes[j][k]->SetMarker(mrk);
								}
								
							}
						}
					}
					
					
				}
				for(j = 0; j < l; j++)
					for(k = 0; k < edge_cut_nodes[j].size(); k++)
						edge_cut_nodes[j][k]->RemMarker(mrk);
				*/
				//g->mesh->ReleaseMarker(mrk);
			}
			
			else if( num_face_center_nodes == 1 )
			{
				if( print ) std::cout << " connect one face node to edge cuts " << std::endl; 
				int side = -1;
				Storage::real c1,c2;
				Storage::integer_array mat1,mat2;
				k = 0;
				for(j = 0; j < l; j++)
					if( node_in_face[j].isValid() )
					{
						edge_nodes[k++] = node_in_face[j];
						side = sides[j];
						c1 = node_in_face[j]->Coords()[side/2];
						mat1 = node_in_face[j]->IntegerArray(g->materials);
						break;
					}
				MarkerType mrk = g->mesh->CreateMarker();
				for(j = 0; j < l; j++)
				{
					edge_cut_nodes[j].insert(edge_cut_nodes[j].end(),edge_cut_nodes2[j].begin(),edge_cut_nodes2[j].end());// do we need to connect degenerate cuts?
					for(k = 0; k < edge_cut_nodes[j].size(); k++) if( !edge_cut_nodes[j][k]->GetMarker(mrk) )
					{
						c2 = edge_cut_nodes[j][k]->Coords()[side/2];
						if( fabs(c1-c2) > 1e-5 ) //not on one plain
						{
							mat2 = edge_cut_nodes[j][k]->IntegerArray(g->materials);
							
							bool connect = true;
							if( mat1.size() == mat2.size() )
							{
								for(unsigned qq = 0; qq < mat1.size(); qq++)
									if( mat1[qq] != mat2[qq] )
										connect = false;
							}
							else connect = false;
							//~ 
							//~ mat_intersection.resize(mat1.size()+mat2.size());
							//~ mat_intersection.resize(std::set_intersection(mat1.begin(),mat1.end(),mat2.begin(),mat2.end(),mat_intersection.begin())-mat_intersection.begin());
							//~ if( !mat_intersection.empty() )
							if( connect )
							{
								
								edge_nodes[1] = edge_cut_nodes[j][k];
								
								Edge new_edge = g->mesh->CreateEdge(edge_nodes).first;
								new_edge->SetMarker(multi_edge);
								//~ if( new_edge->LocalID() == 11801 ) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
								Storage::integer_array mate = new_edge->IntegerArray(g->materials);
								if( mate.empty() )
								{
									mate.insert(mate.begin(),mat1.begin(),mat1.end());
									//~ mate.insert(mate.begin(),mat_intersection.begin(),mat_intersection.end());
									for(unsigned q = 0; q < mate.size(); q++)
										edges_by_material[mate[q]].push_back(new_edge);
								}
								edge_cut_nodes[j][k]->SetMarker(mrk);
							}
							
						}
					}
					/*
					for(k = 0; k < edge_cut_nodes2[j].size(); k++) if( !edge_cut_nodes2[j][k]->GetMarker(mrk) )
					{
						c2 = edge_cut_nodes2[j][k]->Coords()[side/2];
						if( fabs(c1-c2) > 1e-5 ) //not on one plain
						{
							mat2 = edge_cut_nodes2[j][k]->IntegerArray(g->materials);
							
							bool connect = true;
							if( mat1.size() == mat2.size() )
							{
								for(unsigned qq = 0; qq < mat1.size(); qq++)
									if( mat1[qq] != mat2[qq] )
										connect = false;
							}
							else connect = false;
							
							//~ mat_intersection.resize(mat1.size()+mat2.size());
							//~ mat_intersection.resize(std::set_intersection(mat1.begin(),mat1.end(),mat2.begin(),mat2.end(),mat_intersection.begin())-mat_intersection.begin());
							//~ if( !mat_intersection.empty() )
							if( connect )
							{
								
								edge_nodes[1] = edge_cut_nodes2[j][k];
								
								Edge * new_edge = g->mesh->CreateEdge(edge_nodes).first;
								new_edge->SetMarker(hyper_edge);
								//~ if( new_edge->LocalID() == 11801 ) std::cout << __FILE__ << ":" << __LINE__ << std::endl;
								Storage::integer_array mate = new_edge->IntegerArray(g->materials);
								if( mate.empty() )
								{
									mate.insert(mate.begin(),mat1.begin(),mat1.end());
									//~ mate.insert(mate.begin(),mat_intersection.begin(),mat_intersection.end());
									for(unsigned q = 0; q < mate.size(); q++)
										edges_by_material[mate[q]].push_back(new_edge);
								}
								edge_cut_nodes2[j][k]->SetMarker(mrk);
							}
							
						}
					}
					*/
				}
				for(j = 0; j < l; j++)
				{
					for(k = 0; k < edge_cut_nodes[j].size(); k++)
						edge_cut_nodes[j][k]->RemMarker(mrk);
					//~ for(k = 0; k < edge_cut_nodes2[j].size(); k++)
						//~ edge_cut_nodes2[j][k]->RemMarker(mrk);
				}
				g->mesh->ReleaseMarker(mrk);
			}
			
			else 
			{
				//~ if( global_test_number == 10 && m == 27 ) print = true;
				if( print ) std::cout << " connections between edgecuts " << std::endl; 
				Storage::real c1,c2;
				Storage::integer_array mats1,mats2;
				MarkerType mrk = g->mesh->CreateMarker();
				for(j = 0; j < l; j++)
				{
					edge_cut_nodes[j].insert(edge_cut_nodes[j].end(),edge_cut_nodes2[j].begin(),edge_cut_nodes2[j].end());
				}
				for(j = 0; j < l; j++)
				{
					int s = sides[j];
					for(unsigned q1 = 0; q1 < edge_cut_nodes[j].size(); q1++) if( !edge_cut_nodes[j][q1]->GetMarker(mrk) )
					{
						c1 = edge_cut_nodes[j][q1]->Coords()[s/2];
						edge_cut_nodes[j][q1]->SetMarker(mrk);
						mats1 = edge_cut_nodes[j][q1]->IntegerArray(g->materials);
						
						if( print ) 
						{
							std::cout << "from face " << j << " side " << s << " cut " << q1 << " " << edge_cut_nodes[j][q1]->LocalID() << " coord " << c1 << " materials: ";
							for(unsigned qq = 0; qq < mats1.size(); qq++) std::cout << mats1[qq] << " ";
							std::cout << std::endl;
						}
						
						if( mats1.size() < 3 ) continue;
						edge_nodes[0] = edge_cut_nodes[j][q1];
						for(k = 0; k < l; k++) 
						{
							if( print ) std::cout << k << " is side " << sides[k] << " num cuts " << edge_cut_nodes[k].size() << " deteriorate cuts " << edge_cut_nodes2[k].size()  << std::endl;
							if( sides[k] != s )
							{
								for(unsigned q2 = 0; q2 < edge_cut_nodes[k].size(); q2++) 
								{
									if( print ) std::cout << "   "  << q2 << " is  " << edge_cut_nodes[k][q2]->LocalID() << " mrk: " << edge_cut_nodes[k][q2]->GetMarker(mrk) << std::endl;
									if( !edge_cut_nodes[k][q2]->GetMarker(mrk) )
									{
										c2 = edge_cut_nodes[k][q2]->Coords()[s/2];
										
										if( print ) std::cout << "to face " << k << " side " << sides[k] << " cut " << q2 << " " << edge_cut_nodes[k][q2]->LocalID() << " coord " << c2;
										
										
										
										if( fabs(c1-c2) > 1e-5 )
										{
											mats2 = edge_cut_nodes[k][q2]->IntegerArray(g->materials);
											if( print ) 
											{
												std::cout  << " materials: ";
												for(unsigned qq = 0; qq < mats2.size(); qq++) std::cout << mats2[qq] << " ";
												std::cout << std::endl;
											}
											bool connect = true;
											if( mats1.size() == mats2.size() )
											{
												for(unsigned qq = 0; qq < mats1.size(); qq++)
													if( mats1[qq] != mats2[qq] )
														connect = false;
											}
											else connect = false;
											if( connect )
											{
												
												edge_nodes[1] = edge_cut_nodes[k][q2];
												Edge new_edge = g->mesh->CreateEdge(edge_nodes).first;
												if( print ) std::cout << "connect! " << new_edge->LocalID() << std::endl;
												
												//~ new_edge->SetMarker(hyper_edge);
												Storage::integer_array mate = new_edge->IntegerArray(g->materials);
												if( mate.empty() ) 
												{
													mate.insert(mate.begin(),mats1.begin(),mats1.end());
													for(unsigned q = 0; q < mate.size(); q++)
														edges_by_material[mate[q]].push_back(new_edge);
												}	
												edge_cut_nodes[k][q2]->SetMarker(mrk);
											}
										} 
										else if( print ) std::cout << std::endl;
									}
								}
							}
						}
					}
				}
				for(j = 0; j < l; j++)
					for(unsigned q1 = 0; q1 < edge_cut_nodes[j].size(); q1++) edge_cut_nodes[j][q1]->RemMarker(mrk);
				g->mesh->ReleaseMarker(mrk);
				//~ if( global_test_number == 10 && m == 27 ) print = false;
			}
			
			//it may happen that cuts exist on edges inside refined faces, then we should connect them
			//may leave it to the algorithm that creates faces inside the cell
			
			//if ( m == 7833 ) print = true;
			
			//add all potential edges
			std::sort(potential_edges.begin(),potential_edges.end());
			potential_edges.resize(std::unique(potential_edges.begin(),potential_edges.end())-potential_edges.begin());
			if( print ) std::cout << "adding potential edges candidates " << potential_edges.size() << std::endl;
			for(j = 0; j < potential_edges.size(); j++)
			{
				ElementArray<Face> f = potential_edges[j]->getFaces();
				if( print ) std::cout << "potential edge " << potential_edges[j]->LocalID() << " faces ";
				int num_faces = 0;
				bool meet_first = false;
				mat_intersection.clear();
				for(ElementArray<Face>::iterator it = f.begin(); it != f.end(); ++it)
				{
					if( it->GetMarker(face_on_face) )
					{
						Storage::integer_array mats = it->IntegerArray(g->materials);
						if( !meet_first )
						{
							mat_intersection.insert(mat_intersection.end(),mats.begin(),mats.end());
							meet_first = true;
						}
						else
						{
							mat_union.resize(mat_intersection.size());
							mat_union.resize(std::set_intersection(mat_intersection.begin(),mat_intersection.end(),mats.begin(),mats.end(),mat_union.begin())-mat_union.begin());
							mat_intersection.swap(mat_union);
						}
						num_faces++;
						if( print ) 
						{
							std::cout << it->LocalID() << " (" << mats.size() << "): ";
							for(int jj = 0; jj < mats.size(); jj++) std::cout << mats[jj] << " ";
						}
					}
				}
				if( print ) std::cout << " materials " << mat_intersection.size() << std::endl;
				if( mat_intersection.empty() )
				{
					Storage::integer_array mats = potential_edges[j]->IntegerArray(g->materials);
					if( print ) std::cout << "adding edge " << potential_edges[j]->LocalID() << " to materials ";
					for(Storage::integer_array::iterator kt = mats.begin(); kt != mats.end(); ++kt)
					{
						edges_by_material[*kt].push_back(potential_edges[j]);
						if( print ) std::cout << *kt << " ";
					}
					if( print ) std::cout << std::endl;
				}
				
			}
		
			//Now create all the faces that lay inside the cell
			
			if( print ) std::cout << "create inner faces " << m << std::endl;
			for(tiny_map<int, dynarray<Edge,128>, 64 >::iterator it = edges_by_material.begin(); it != edges_by_material.end(); ++it) // iterate over materials
			{
				std::sort(it->second.begin(),it->second.end());
				it->second.resize(std::unique(it->second.begin(),it->second.end())-it->second.begin());
				
				
				if( print ) 
				{
					std::cout << "material " << it->first << " total edges " << it->second.size() << std::endl;
					for(k = 0; k < it->second.size(); k++)
					{
						ElementArray<Node> n = it->second[k]->getNodes();
						std::cout << k << " " << it->second[k]->LocalID();
						std::cout << " don't start " << !it->second[k]->GetMarker(edge_on_face);
						{
							std::cout << " materials ";
							Storage::integer_array mat = it->second[k]->IntegerArray(g->materials);
							for(j = 0; j < mat.size(); j++) std::cout << mat[j] << " ";
						}
						std::cout << n[0].LocalID() << " " << n[0].Coords()[0] << "," << n[0].Coords()[1] << "," << n[0].Coords()[2] << " ";
						std::cout << n[1].LocalID() << " " << n[1].Coords()[0] << "," << n[1].Coords()[1] << "," << n[1].Coords()[2] << " ";
						std::cout << std::endl;
					}
					
					std::cout << "start gather" << std::endl;
				}

				
				unsigned num_inner_edges = 0, num = 0;
				dynarray<Edge,128> edges(it->second.size());
				for(k = 0; k < it->second.size(); k++)
					if( !it->second[k]->GetMarker(edge_on_face) )
					{
						edges[num++] = it->second[k]; 
						num_inner_edges++;
					}
				for(k = 0; k < it->second.size(); k++)
					if( it->second[k]->GetMarker(edge_on_face) )
						edges[num++] = it->second[k]; 
				
				//std::vector<Edge *> & edges = it->second;
				
				incident_matrix<Edge> matrix(edges.begin(),edges.end(),num_inner_edges);
				
				
				/*
				for(unsigned q = 0; q < edges.size(); q++)
					if( matrix.get_element(q)->GetMarker(multi_edge) )
					{
						//TODO: correctly detect number of visits to edge
						adjacent<Element> e = matrix.get_element(q)->BridgeAdjacencies(NODE,EDGE);
						mat_intersection.clear();
						bool first = false;
						for(adjacent<Element>::iterator qt = e.begin(); qt != e.end(); ++qt)
						{
							if( qt->GetMarker(edge_on_face) )
							{
								Storage::integer_array matse = qt->IntegerArrayDF(g->materials);
								if( std::binary_search(matse.begin(),matse.end(),it->first) )
								{
									if( first )
										mat_intersection.insert(mat_intersection.end(),matse.begin(),matse.end());
									else
									{
										mat_union.resize(std::max(mat_intersection.size(),matse.size()));
										mat_union.resize(std::set_intersection(mat_intersection.begin(),mat_intersection.end(),matse.begin(),matse.end(),mat_union.begin())-mat_union.begin());
										mat_intersection.swap(mat_union);
									}
								}
							}
						}
						if( print ) std::cout << "set visit count for edge " << matrix.get_element(q)->LocalID() << " " << mat_intersection.size() << std::endl;
						matrix.set_visit(q,mat_intersection.size());
					}
				*/
				
				if( print ) matrix.print_matrix();
				
				MarkerType mrk = g->mesh->CreateMarker();
				int num_loops = 0;
				
				ElementArray<Edge> loop(g->mesh);
				while(matrix.find_shortest_loop(loop))
				{ 
					
					int num_edge_on_face = 0, num_edge_on_edge = 0;
					tiny_map<int,int,128> num_faces;
					for(k = 0; k < loop.size(); k++)
					{
						if( loop[k]->GetMarker(edge_on_edge) ) num_edge_on_edge++;
						if( loop[k]->GetMarker(edge_on_face) ) 
						{
							num_edge_on_face++;
							num_faces[loop[k]->Integer(g->edge_face_number)]++;
						}
					}
					
					num_loops++;
					for(k = 0; k < loop.size(); k++) loop[k]->SetMarker(mrk);
					
					if( num_edge_on_edge == loop.size() ) continue;
					if( num_edge_on_edge + num_edge_on_face == loop.size() )
					{
						if( num_faces.size() == 1 ) continue;
					}
					//probably skip other variant
					
					if( print ) std::cout << "found loop " << loop.size() << std::endl;
					if( loop.size() > 2 )
					{
						
						
						if( print )
							for(k = 0; k < loop.size(); k++)
							{
								std::cout << k << " " << loop[k]->GetHandle() << " material";
								Storage::integer_array mat = loop[k]->IntegerArray(g->materials);
								for(j = 0; j < mat.size(); j++) std::cout << " " << mat[j];
								for(j = 0; j < edges.size(); j++)
									if( loop[k] == edges[j] ) 
									{
										std::cout << " edge " << j;
										break;
									}
								if( loop[k]->GetMarker(edge_on_face) )
									std::cout << " on face ";
								else
									std::cout << " inside ";
								std::cout << std::endl;
							}

							Face new_face = g->mesh->CreateFace(loop).first;

						if( !new_face->GetMarker(face_on_face) )
						{
							Storage::integer_array mat = new_face->IntegerArray(g->materials);
								

							if( mat.empty() )
							{
								Storage::integer_array mt = loop[0]->IntegerArray(g->materials);
								mat_intersection.clear();
								mat_intersection.insert(mat_intersection.end(),mt.begin(),mt.end());
								mat_union.resize(mat_intersection.size());
								for(k = 1; k < loop.size(); k++)
								{
									mt = loop[k]->IntegerArray(g->materials);
									mat_union.resize(std::set_intersection(mat_intersection.begin(),mat_intersection.end(),mt.begin(),mt.end(),mat_union.begin())-mat_union.begin());
									mat_intersection.swap(mat_union);
								}
								mat.insert(mat.end(),mat_union.begin(),mat_union.end());
							}
								
							if( !std::binary_search(mat.begin(),mat.end(),it->first) ) throw -1;
								
							for(unsigned q = 0; q < mat.size(); q++)
							{
								if( print ) 
								{
									std::cout << __FILE__ << ":" << __LINE__ << " "  << mat[q] << " added face " << new_face->GetHandle() << " neighbours " << new_face->nbAdjElements(CELL) << " ";
									if( new_face->nbAdjElements(CELL) == 2 )
										std::cout << new_face->BackCell()->GetHandle() << " " <<  new_face->BackCell()->LocalID() << " " << new_face->FrontCell()->GetHandle() << " " << new_face->FrontCell()->LocalID() << std::endl;
									else std::cout << std::endl;
								}
								//faces_by_material[mat[q]].push_back(new_face);
								
							}
							inner_faces.push_back(new_face);
						}
					}
					//else if( !loop.empty() ) std::cout << __FILE__ << ":" << __LINE__ << " skip degenerate face " << loop.size() << std::endl;
				} 
				
				
				int unused = static_cast<int>(edges.size());
				for(unsigned q = 0; q < edges.size(); q++)
				{
					//~ if( !edges[q]->GetMarker(mrk) )
						//~ std::cout << "unused edge: " << edges[q]->LocalID() << " faces " << edges[q]->nbAdjElements(FACE) << " cells " << edges[q]->nbAdjElements(CELL) << std::endl;
					unused -= edges[q]->GetMarker(mrk);
					edges[q]->RemMarker(mrk);
				}
				
				g->mesh->ReleaseMarker(mrk);
				
				
				
				if( unused ) 
				{
					//~ std::cout << "total unused: " << unused << " total edges :" << edges.size() << " total loops: " << num_loops << std::endl;
					//trigger_problem = 6;
				}
			}
			
			
			std::sort(inner_faces.begin(),inner_faces.end());
			inner_faces.resize(std::unique(inner_faces.begin(),inner_faces.end())-inner_faces.begin());
			std::sort(faces_on_face.begin(),faces_on_face.end());
			faces_on_face.resize(std::unique(faces_on_face.begin(),faces_on_face.end())-faces_on_face.begin());
			unsigned num_inner_faces = inner_faces.size();

			for(dynarray<HandleType,128>::iterator it = faces_on_face.begin(); it != faces_on_face.end(); ++it)
				inner_faces.push_back(Face(g->mesh,*it));
			
			incident_matrix<Face> matrix(inner_faces.begin(),inner_faces.end(),num_inner_faces);
				
			if( print ) matrix.print_matrix();
			
			int loops_found = 0;
			ElementArray<Face> loop;
			while(matrix.find_shortest_loop(loop))
			{ 
				

				if( print ) std::cout << "found loop " << loop.size() << std::endl;
				
				if( print )
					for(k = 0; k < loop.size(); k++)
					{
						std::cout << k << " " << loop[k]->LocalID() << " material";
						Storage::integer_array mat = loop[k]->IntegerArray(g->materials);
						for(j = 0; j < mat.size(); j++) std::cout << " " << mat[j];
						for(j = 0; j < inner_faces.size(); j++)
							if( loop[k] == inner_faces[j] ) 
							{
								std::cout << " face " << j;
								break;
							}
						if( loop[k]->GetMarker(face_on_face) )
							std::cout << " on face ";
						else
							std::cout << " inside ";
						std::cout << std::endl;
					}
				
				
				if( loop.size() > 2 )
				{
					loops_found++;
					//try
					{
						Cell new_cell = g->mesh->CreateCell(loop).first;
						if( print ) std::cout << __FILE__ << ":" << __LINE__ << " new cell " << new_cell->GetHandle() << " id " << new_cell->LocalID() << std::endl;
						
						
						//~ if( new_cell->HaveData(g->mesh->TopologyErrorTag() ) )
							//~ global_problem = 1;
						
						if( print )
							std::cout << "added cell " << new_cell->LocalID() << std::endl;
						
						if( loop.size() < 4 ) new_cell->Integer(g->problem) = 10;
						
						{
							Storage::integer_array fmat = loop[0]->IntegerArray(g->materials);
							mat_intersection.clear();
							mat_intersection.insert(mat_intersection.end(),fmat.begin(),fmat.end());
							for(unsigned q = 1; q < loop.size(); q++)
							{
								fmat = loop[q]->IntegerArray(g->materials);
								mat_union.resize(mat_intersection.size());
								mat_union.resize(std::set_intersection(mat_intersection.begin(),mat_intersection.end(),fmat.begin(),fmat.end(),mat_union.begin())-mat_union.begin());
								mat_intersection.swap(mat_union);
							}
							/*
							if( !mat_intersection.empty() )
								new_cell->Integer(g->cell_material) = mat_intersection[0];
							else
							{
								new_cell->Integer(g->cell_material) = 6;
								new_cell->Integer(g->problem) = 6;
								Storage::real cnt[3], cntf[3];
								new_cell->Centroid(cnt);
								mat_ret_type mats = g->get_material_types(cnt);
								new_cell->Integer(g->cell_material) = mats[0];
							}
							*/
							
							bool need_fix = false;
							Storage::real cnt[3], cntf[3];
							new_cell->Centroid(cnt);
							mat_ret_type mats = g->get_material_types(cnt);
							//new_cell->Integer(g->cell_material) = mats[0];
														
							if( mat_intersection.size() == 1 )
							{
								new_cell->Integer(g->cell_material) = mat_intersection[0];
								/*
								if( mats.size() == 1 && mats[0] == mat_intersection[0] )
									new_cell->Integer(g->cell_material) = mat_intersection[0];//it->first;
								
								else if( mats.size() > 1 )
								{
									mat_union.resize(std::max(mat_intersection.size(),mats.size()));
									mat_union.resize(std::set_intersection(mats.begin(),mats.end(),mat_intersection.begin(),mat_intersection.end(),mat_union.begin())-mat_union.begin());
									if( mat_union.size() == 1 )
										new_cell->Integer(g->cell_material) = mat_union[0];
									else need_fix = true;
								}
								
								else
									need_fix = true;
								*/
							}
							else if( mat_intersection.size() > 1 )
							{
								
								if( mats.size() == 1 )
									new_cell->Integer(g->cell_material) = mats[0];//it->first;
								else
								{
									mat_union.resize(std::max(mat_intersection.size(),mats.size()));
									mat_union.resize(std::set_intersection(mats.begin(),mats.end(),mat_intersection.begin(),mat_intersection.end(),mat_union.begin())-mat_union.begin());
									if( mat_union.size() == 1 )
									{
										new_cell->Integer(g->cell_material) = mat_union[0];
										//new_cell->Integer(g->problem) = 6;
									}
									else need_fix = true;
								}
							}
							else 
								need_fix = true;
							
							if( need_fix )
							{
								double minarea = 1e20;
								int mink = -1;
								for( k = 0; k < loop.size(); k++) //search smallest face
								{
									Storage::integer_array fmat = loop[k]->IntegerArray(g->materials);
									if( fmat.size() > 1 )
									{
										double area = loop[k]->Area();
										if( area < minarea )
										{
											minarea = area;
											mink = k;
										}
									}
								}
								
								
								if( mink == -1 )
									new_cell->Integer(g->cell_material) = mats[0];
								else
								{
									loop[mink]->Centroid(cntf);
									for(k = 0; k < 3; k++) cnt[k] = cnt[k]*0.25 + cntf[k]*0.75;
									mats = g->get_material_types(cnt);
									new_cell->Integer(g->cell_material) = mats[0];
								}
								
								//new_cell->Integer(g->problem) = 6;
							}
							
						}
						new_cell->Integer(g->parent) = m;
						if( new_cell->GetGeometricType() == Element::MultiPolygon )
						{
							trigger_problem = 3;
							throw -1;
						}
						
						unsigned num_nodes = new_cell->nbAdjElements(NODE);
						
						for(k = 0; k < loop.size(); k++)
							if( loop[k]->nbAdjElements(NODE) == num_nodes )
							{
								new_cell->Integer(g->problem) = 10;
								break;
							}
						
						g->cells[m].mr->push_back(new_cell->GetHandle());
					}
					//catch(...) 
//					{
//						bool found = false;
//						for(unsigned q = 0; q < loop.size(); q++) if( loop[q]->nbAdjElements(CELL) == 2 )
//						{
//							adjacent<Cell> nb = loop[q]->getCells();
//							for(unsigned qq = 0; qq < nb.size(); qq++)
//								nb[qq].Integer(g->problem) = 1;
//							found = true;
//						}
//						if( found ) std::cout << "error found!" << std::endl;
//						else std::cout << "error not found!" << std::endl;
//					}
				}
				else if( !loop.empty() ) 
					std::cout  << __FILE__ << ":" << __LINE__ << " skip degenerate cell " << loop.size() << std::endl;
			} 
			
			if( !matrix.all_visited() ) 
			{
				//matrix.print_matrix();
				//std::cout << "not all faces used, total cells " << g->cells[m].mr->size() << " parent " << m << " marked as problem 4" << std::endl;
				global_problem = true;
				trigger_problem = 4;
				//~ for(unsigned q = 0; q < g->cells[m].mr->size(); q++)
					//~ (*g->cells[m].mr)[q]->Integer(g->problem) = 5;
			}
		}
		
		
		
		
		for(tiny_map<int, dynarray<Edge,128> ,64>::iterator it = edges_by_material.begin(); it != edges_by_material.end(); ++it) // iterate over materials
			for(dynarray<Edge,128>::iterator jt = it->second.begin(); jt != it->second.end(); jt++)
			{
				(*jt)->RemMarker(edge_on_face);
				(*jt)->RemMarker(edge_on_edge);
				(*jt)->RemMarker(multi_edge);
			}
			
		for(j = 0; j < potential_edges.size(); j++)
			potential_edges[j]->RemMarker(edge_on_edge);
		
		g->mesh->ReleaseMarker(edge_on_edge);
			
		g->mesh->ReleaseMarker(edge_on_face);
		g->mesh->ReleaseMarker(multi_edge);
		
		for(dynarray<HandleType,128>::iterator it = faces_on_face.begin(); it != faces_on_face.end(); ++it) // iterate over materials
			Face(g->mesh,*it)->RemMarker(face_on_face);
		
		//~ for(std::map<int, std::vector<Face *> >::iterator it = faces_by_material.begin(); it != faces_by_material.end(); ++it) // iterate over materials
			//~ for(std::vector<Face *>::iterator jt = it->second.begin(); jt != it->second.end(); jt++)
				//~ (*jt)->RemMarker(face_on_face);
		g->mesh->ReleaseMarker(face_on_face);
		
		for(k = 0; k < g->cells[m].mr->size(); k++)
		{
			if( Cell(g->mesh,(*g->cells[m].mr)[k])->Integer(g->problem) == 10 )
			{
				Cell ck = Cell(g->mesh,(*g->cells[m].mr)[k]);
				unsigned nv = ck->nbAdjElements(NODE);
				ElementArray<Face> faces = ck->getFaces();
				for(j = 0; j < faces.size(); j++)
				{
					if( faces[j].nbAdjElements(NODE) == nv )
					{
						Cell n = ck->Neighbour(faces[j]);
						if( n.isValid() )
						{
							Storage::integer mat = n->Integer(g->cell_material);
							ElementArray<Cell> unite(g->mesh,2);
							unite[0] = ck;
							unite[1] = n;
							
							{
								unsigned q = 0;
								while(q < g->cells[m].mr->size())
								{
									if( (*g->cells[m].mr)[q] == unite[0]->GetHandle() || (*g->cells[m].mr)[q] == unite[1]->GetHandle() )
										g->cells[m].mr->erase(g->cells[m].mr->begin()+q);
									else q++;
								}
							}
							Cell out = Cell::UniteCells(unite,g->octree_node);
							//~ std::cout << "cells united! " << out << " id " << out->LocalID() << " parent " << m <<  std::endl;
							if( out.isValid() )
							{
								out->Integer(g->cell_material) = mat;
								out->Integer(g->parent) = m;
								if( out->GetGeometricType() == Element::MultiPolygon )
								{
									trigger_problem = 3;
									throw -1;
								}
								g->cells[m].mr->push_back(out->GetHandle());
							}
							else
							{
								std::cout << "unite failed" << std::endl;
								g->cells[m].mr->push_back(unite[0]->GetHandle());
								g->cells[m].mr->push_back(unite[1]->GetHandle());
							}
						}
					}
				}
			}
		}
		
	}
	
	
	Storage::real vsum = 0;
	for(int k = 0; k < g->cells[m].mr->size(); k++)
	{
		Cell ck = Cell(g->mesh,(*g->cells[m].mr)[k]);
		g->cell_to_INMOST(g,m,ck);
		vsum += ck->Volume();
		if( ck->Volume() <= 0 ) trigger_problem = 1;
	}
	if( vsum < 0 ) trigger_problem = 1;
	if( fabs(g->cells[m].vol - vsum) > 1e-3 ) 
	{
		std::cout << "octree cell vol " << g->cells[m].vol << " cells sum " << vsum << " diff " << g->cells[m].vol - vsum << " l " << l << std::endl;
		trigger_problem = 14;
	}
	if( trigger_problem )
		for(int k = 0; k < g->cells[m].mr->size(); k++)
			Cell(g->mesh,(*g->cells[m].mr)[k])->Integer(g->problem) = trigger_problem;
	g->cell_init_data(g,m);
}

void vertDelete(struct grid * g, int m)
{
	int i;
	if( m == -1 ) return;
	vertDestroyINMOST(g,m);
	g->vert_destroy_data(g,m);
	g->verts[m].busy = 0;
	for (i = 0; i < 1<<DIM; i++)
		g->verts[m].env[i] = -1;
	C_ADD_ARR(g,vempty,m);
	C_INC_ARR(g,vempty,int,VBUFFER);
}



int vertCheck(struct grid * g, int m)
{
	const int lookdir[6][4] = {{0,2,6,4},{1,3,5,7},{0,1,4,5},{2,3,6,7},{0,1,2,3},{4,5,6,7}};
	int i, j, k, max = 5, d = 0,flag;
	int cell[1<<DIM];
	if( m == -1 ) return 0;
	for(i = 0; i < DIM*2; i++)
	{
		int flag = 1;
		for(k = 0; k < 1<<(DIM-1); k++)
			if( g->verts[m].env[lookdir[i][k]] != -1 )
			{
				flag = 0;
				break;
			}
		if( flag )
			d++;
	}
	
	if( d == 0 ) max = 5;
	else if( d == 1 ) max = 3;
	else if( d == 2 ) max = 2;
	else if( d == 3 ) max = 1;
	else return 1;
	
	k = 0;
	for (i = 0; i < 1<<DIM; i++) if( g->verts[m].env[i] != -1 )
	{
		flag = 1;
		for (j = 0; j < k; j++)
			if (g->verts[m].env[i] == cell[j])
			{
				flag = 0;
				break;
			}
		if (flag)
			cell[k++] = g->verts[m].env[i];
	}
	if( k < max ) return 1;
	return 0;
}

int vertNextFree(struct grid * g)
{
	int ret;
	if( C_CNT_ARR(g,vempty) )
	{
		ret = C_POP_ARR(g,vempty);
		g->verts[ret].busy = 1;
		return ret;
	}
	ret = C_CNT_ARR(g,verts);
	g->verts[ret].busy = 1;
	C_CNT_ARR(g,verts)++;
	C_INC_ARR(g,verts,struct vert,VBUFFER);
	return ret;
}




void cellDelete(struct grid * g, int m)
{
	int i;
	if( m == -1 ) return;
	cellDestroyINMOST(g,m);
	g->cell_destroy_data(g,m);
	g->cells[m].busy = 0;
	if( !g->cells[m].leaf )
	{
		for (i = 0; i < 1<<DIM; i++)
			cellDelete(g,g->cells[m].children[i]);
	}
	C_ADD_ARR(g,cempty,m);
	C_INC_ARR(g,cempty,int,CBUFFER);
	delete g->cells[m].mr;
	delete g->cells[m].data;
}

int cellNextFree(struct grid * g)
{
	int ret;
	if( C_CNT_ARR(g,cempty) )
	{
		ret = C_POP_ARR(g,cempty);
		g->cells[ret].busy = 1;
		return ret;
	}
	ret = C_CNT_ARR(g,cells);
	g->cells[ret].busy = 1;
	C_CNT_ARR(g,cells)++;	
	C_INC_ARR(g,cells,struct cell,CBUFFER);
	return ret;
}

void cellInit(struct grid * g, int cell, int parent)
{
	int i;
	g->cells[cell].mr = new cell_vector();
	g->cells[cell].data = new data_storage();
	g->cells[cell].leaf = 1;
	g->cells[cell].level = 0;
	g->cells[cell].parent = parent;
	for(i = 0; i < 1<<DIM; i++)
	{
		g->cells[cell].children[i] = -1;
		g->cells[cell].vertexes[i] = -1;
	}
	//data
	g->cells[cell].side[0] = g->cells[cell].side[1] = g->cells[cell].side[2] = 0.0;
	g->cells[cell].center[0] = g->cells[cell].center[1] = g->cells[cell].center[2] = 0.0;
}

void vertInit(struct grid * g, int vert)
{
	int i;
	g->verts[vert].mv = InvalidHandle();
	for(i = 0; i < 1<<DIM; i++)
		g->verts[vert].env[i] = -1;
}


void gridRecreateINMOST(struct grid * g)
{
	int i;
	//~ for(i = 0; i < C_CNT_ARR(g,cells); i++)
		//~ if(!g->cells[i].mr.empty() )
			//~ cellDestroyINMOST(g,i);
	for(i = 0; i < C_CNT_ARR(g,cells); i++)
	{
		if( g->cells[i].busy )
		{
			if( g->cells[i].leaf && g->cells[i].mr->empty() )
			{
				//~ try
				//~ {
					cellCreateINMOST(g,i);
				//~ }
				//~ catch(...)
				//~ {
					//~ cellCreateINMOST(g,i,true);
				//~ }
				//~ g->cell_init_data(g,i);
			}
		}
	}
	for(i = 0; i < C_CNT_ARR(g,cells); i++)
	{
		if( g->cells[i].busy )
		{
			if( !g->cells[i].leaf && !g->cells[i].mr->empty() )
			{
				cellDestroyINMOST(g,i);
				//~ g->cell_init_data(g,i);
			}
		}
	}
	
	if( remove_orphan_elements )
	{
	for(Mesh::iteratorFace f = g->mesh->BeginFace(); f != g->mesh->EndFace(); f++)
		if( f->nbAdjElements(CELL) == 0 ) f->Delete();
		
	for(Mesh::iteratorEdge f = g->mesh->BeginEdge(); f != g->mesh->EndEdge(); f++)
		if( f->nbAdjElements(CELL) == 0 ) f->Delete();
	}
}


void vertSplitData(struct grid * g, int big_cell, int nverts, int * verts, int * isnew)
{
	g->vert_interpolate_data(g,big_cell,nverts,verts,isnew);
}


int cellSplit(struct grid * g, int m)
{
#if (DIM == 2)
	const int ndirs = 4;
	const int nface = 4;
	const int nsubface = 1;
	const int kjl[4][1][4]= 
	{
		{{0,0,2,3}},
		{{1,0,1,3}},
		{{2,3,1,0}},
		{{3,3,2,0}}
	};
	const int qjl[4][8] =
	{
		{0,0,2,2,  2,3,0,1},
		{0,1,0,1,  1,0,3,2},
		{1,1,3,3,  2,3,0,1},
		{2,3,2,3,  1,0,3,2}
	};
#elif (DIM == 3)
	const int ndirs = 18;
	const int nface = 6;
	const int nsubface = 5;
	const int nedge = 12;
	const int kjl[6][5][4] = 
	{
		{{0,0,6,7},{10,0,6,3},{14,0,6,5},{16,6,0,3},{12,6,0,5}},
		{{1,0,5,7},{6,0,5,3},{14,0,5,6},{15,5,0,3},{8,5,0,6}},
		{{2,0,3,7},{10,0,3,6},{6,0,3,5},{7,3,0,6},{11,3,0,5}},
		{{3,7,1,0},{13,7,1,4},{17,7,1,2},{15,1,7,4},{11,1,7,2}},
		{{4,7,2,0},{9,7,2,4},{17,7,2,1},{16,2,7,4},{7,2,7,1}},
		{{5,7,4,0},{9,7,4,2},{13,7,4,1},{12,4,7,2},{8,4,7,1}}
	};
	const int mjl[12][4] =
	{
		{6,0,1,7} ,{7,3,2,4} ,{8,5,4,2} ,{9,6,7,1} ,
		{10,0,2,7},{11,3,1,4},{12,6,4,1},{13,5,7,2},
		{14,0,4,7},{15,5,1,2},{16,6,2,1},{17,3,7,4}
	};
	const int qjl[18][16] =
	{
		{0,0,2,2,4,4,6,6,  6,7,4,5,2,3,0,1},
		{0,1,0,1,4,5,4,5,  5,4,7,6,1,0,3,2},
		{0,1,2,3,0,1,2,3,  3,2,1,0,7,6,5,4},
		{1,1,3,3,5,5,7,7,  6,7,4,5,2,3,0,1},
		{2,3,2,3,6,7,6,7,  5,4,7,6,1,0,3,2},
		{4,5,6,7,4,5,6,7,  3,2,1,0,7,6,5,4},
		{0,1,0,1,0,1,0,1,  1,0,3,2,5,4,7,6},
		{2,3,2,3,2,3,2,3,  1,0,3,2,5,4,7,6},
		{4,5,4,5,4,5,4,5,  1,0,3,2,5,4,7,6},
		{6,7,6,7,6,7,6,7,  1,0,3,2,5,4,7,6},
		{0,0,2,2,0,0,2,2,  2,3,0,1,6,7,4,5},
		{1,1,3,3,1,1,3,3,  2,3,0,1,6,7,4,5},
		{4,4,6,6,4,4,6,6,  2,3,0,1,6,7,4,5},
		{5,5,7,7,5,5,7,7,  2,3,0,1,6,7,4,5},
		{0,0,0,0,4,4,4,4,  4,5,6,7,0,1,2,3},
		{1,1,1,1,5,5,5,5,  4,5,6,7,0,1,2,3},
		{2,2,2,2,6,6,6,6,  4,5,6,7,0,1,2,3},
		{3,3,3,3,7,7,7,7,  4,5,6,7,0,1,2,3}
	};
#endif
	int i,j, chld[(1 << DIM)], dirs[ndirs], repl[ndirs], v,c;
	if( m == -1 ) return 0;
	if( !g->cells[m].busy ) return 0;
	if( g->cells[m].leaf )
	{
		for(i = 0; i < ndirs; i++)
		{
			dirs[i] = -1;
			repl[i] = 0;
		}
		
		
		cellDestroyINMOST(g,m);
		for(i = 0; i < (1 << DIM); i++)
		{
			v = g->cells[m].vertexes[i];
			for(j = 0; j < (1 << DIM); j++)
				if( i != j && i != ((1<<DIM)-1)-j)
				{
					c = g->verts[v].env[j];
					if( c != -1 && g->cells[c].level <= g->cells[m].level )
						cellDestroyINMOST(g,c);
				}
		}
		
		for(i = 0; i < (1 << DIM); i++)
		{
			chld[i] = cellNextFree(g);
			cellInit(g,chld[i],m);
			g->cells[chld[i]].level = g->cells[m].level+1;
			g->cells[m].children[i] = chld[i];

			
			
			g->cells[chld[i]].side[0] = g->cells[m].side[0]*0.5;
			g->cells[chld[i]].side[1] = g->cells[m].side[1]*0.5;
			g->cells[chld[i]].side[2] = g->cells[m].side[2]*0.5;
			
			g->cells[chld[i]].center[0] = g->cells[m].center[0] + (1*(i%2)-0.5)*g->cells[chld[i]].side[0];
			g->cells[chld[i]].center[1] = g->cells[m].center[1] + (1*(i/2%2)-0.5)*g->cells[chld[i]].side[1];
			g->cells[chld[i]].center[2] = g->cells[m].center[2] + (1*(i/4)-0.5)*g->cells[chld[i]].side[2];
		}
		
		cellSplitData(g,m);
		
		for(i = 0; i < nface; i++)
		{
			v = g->cells[m].vertexes[ kjl[i][0][1] ];
			c = g->verts[v].env[ kjl[i][0][2] ];
			if( c == -1 ) continue;
			if( g->cells[c].level > g->cells[m].level )
			{
				for(j = 0; j < nsubface; j++ )
				{
					v = g->cells[m].vertexes[ kjl[i][j][1] ];
					c = g->verts[v].env[ kjl[i][j][2] ];
					if( c == -1 ) continue;
					dirs[ kjl[i][j][0] ] = g->cells[c].vertexes[ kjl[i][j][3] ];
				}
			}		
			else if( g->cells[c].level < g->cells[m].level )
				cellSplit(g, c);
		}
		for(i = 0; i < nedge; i++)
		{
			v = g->cells[m].vertexes[ mjl[i][1] ];
			c = g->verts[v].env[ mjl[i][2] ];
			if( c == -1 ) continue;
			if( dirs[mjl[i][0]] == -1 && g->cells[c].level > g->cells[m].level )
				dirs[mjl[i][0]] = g->cells[c].vertexes[ mjl[i][3] ];
			else if( g->cells[c].level < g->cells[m].level )
				cellSplit(g, c);
		}
		v = vertNextFree(g);
		vertInit(g,v);
		//reinterpolate data to vertex
		for(i = 0; i < (1<<DIM); i++ )
		{
			g->cells[ chld[i] ].vertexes[i] = g->cells[m].vertexes[i];
			g->verts[ g->cells[m].vertexes[i] ].env[((1<<DIM)-1)-i] = chld[i];
			
			g->cells[ chld[i] ].vertexes[((1<<DIM)-1)-i] = v;
			g->verts[v].env[i] = chld[i];
		}
		for(i = 0; i < ndirs; i++ )
		{
			if( dirs[i] == -1 )
			{
				dirs[i] = vertNextFree(g);
				vertInit(g,dirs[i]);
				repl[i] = 1;
			}
		}
		for(i = 0; i < ndirs; i++)
		{
			if( dirs[i] != -1)
			{
				for(j = 0; j < (1<<DIM); j++)
					g->verts[dirs[i]].env[j] = g->verts[g->cells[m].vertexes[qjl[i][j]]].env[qjl[i][j+(1<<DIM)]];
			}
			else TSNH
		}
		for(i = 0; i < ndirs; i++)
		{
			if( dirs[i] != -1 ) for(j = 0; j < (1<<DIM); j++)
			{
				c = g->verts[dirs[i]].env[j];
				if( c != -1 && g->cells[c].level > g->cells[m].level )
					g->cells[c].vertexes[((1<<DIM)-1)-j] = dirs[i];
			}
		}
		g->cells[m].leaf = 0;
		vertCreateINMOST(g,v);
		g->vert_init_data(g,v);
		for(i = 0; i < ndirs; i++) if( repl[i] == 1 )
		{
			g->vert_init_data(g,dirs[i]);
			vertCreateINMOST(g,dirs[i]);
		}
		vertSplitData(g,m,ndirs,dirs,repl);
	}
	return 1;
}

int cellUnite(struct grid * g, int m)
{
	int i,j, c,v,flag;
	if( m == -1 ) return 0;
	if( !g->cells[m].busy ) return 0;
	if( !g->cells[m].leaf )
	{
		flag = 1;
		for(i = 0; i < (1 << DIM); i++)
			for(j = 0; j < (1 << DIM); j++)
				if( j != i && j != ((1 << DIM) - 1)-i )
				{
					v = g->cells[g->cells[m].children[i]].vertexes[i];
					c = g->verts[v].env[j];
					if( c != -1 && g->cells[c].level > g->cells[m].level+1 )
						flag = 0;
				}
		if( !flag ) return 0;
		
		
		//~ cellUniteData(g,m);	

		for(i = 0; i < (1 << DIM); i++)
			cellDestroyINMOST(g,g->cells[m].children[i]);
		for(i = 0; i < (1 << DIM); i++)
		{
			v = g->cells[m].vertexes[i];
			for(j = 0; j < (1 << DIM); j++)
				if( i != j && i != ((1<<DIM)-1)-j)
				{
					c = g->verts[v].env[j];
					if( c != -1 && g->cells[c].level <= g->cells[m].level )
						cellDestroyINMOST(g,c);
				}
		}

		
		
		for(i = 0; i < (1 << DIM); i++)
			cellUnite(g,g->cells[m].children[i]);		
		for(i = 0; i < (1 << DIM); i++)
			for(j = 0; j < (1 << DIM); j++)
				g->verts[g->cells[g->cells[m].children[i]].vertexes[j]].env[((1<<DIM)-1)-j] = m;
		i = 0; j = (1 << DIM)-1;
		vertDelete(g,g->cells[g->cells[m].children[i]].vertexes[j]);
		i = 0;
		for(j = 1; j < (1 << DIM)-1; j++)
			if( vertCheck(g,g->cells[g->cells[m].children[i]].vertexes[j]) )
				vertDelete(g,g->cells[g->cells[m].children[i]].vertexes[j]);
		i = (1 << DIM)-1;
		for(j = 1; j < (1 << DIM)-1; j++)
			if( vertCheck(g,g->cells[g->cells[m].children[i]].vertexes[j]) )
				vertDelete(g,g->cells[g->cells[m].children[i]].vertexes[j]);
			
			
#if (DIM==3)
		i = 1; j = 3;
		if( vertCheck(g,g->cells[g->cells[m].children[i]].vertexes[j]) )
			vertDelete(g,g->cells[g->cells[m].children[i]].vertexes[j]);
			
		i = 1; j = 5;
		if( vertCheck(g,g->cells[g->cells[m].children[i]].vertexes[j]) )
			vertDelete(g,g->cells[g->cells[m].children[i]].vertexes[j]);
			
			
		i = 2; j = 3;
		if( vertCheck(g,g->cells[g->cells[m].children[i]].vertexes[j]) )
			vertDelete(g,g->cells[g->cells[m].children[i]].vertexes[j]);
		i = 2; j = 6;
		if( vertCheck(g,g->cells[g->cells[m].children[i]].vertexes[j]) )
			vertDelete(g,g->cells[g->cells[m].children[i]].vertexes[j]);
			
			
		i = 4; j = 5;
		if( vertCheck(g,g->cells[g->cells[m].children[i]].vertexes[j]) )
			vertDelete(g,g->cells[g->cells[m].children[i]].vertexes[j]);
		i = 4; j = 6;
		if( vertCheck(g,g->cells[g->cells[m].children[i]].vertexes[j]) )
			vertDelete(g,g->cells[g->cells[m].children[i]].vertexes[j]);
#endif
		for(i = 0; i < (1 << DIM); i++)
		{
			cellDelete(g,g->cells[m].children[i]);
			g->cells[m].children[i] = -1;	
		}
		g->cells[m].leaf = 1;
	}
	return 1;
}





int cellAllAround(struct grid * g, int m, int around[(1<<DIM)*((1<<DIM)-1)])
{

	int i,j,v,c,k = 0;
	for(i = 0; i < (1<<DIM); i++)
	{
		v = g->cells[m].vertexes[i];
		for(j = 0; j < (1<<DIM); j++)
			if( i != (1<<DIM)-1-j && (c = g->verts[v].env[j]) != -1)
				add_element(around,&k,c);
	}
	return k;
}






void gridInit(struct grid * g, int n[3])
{
	double ttt = Timer();
	int i,j,k,l,m, c[3], v, r, env[1<<DIM], flag;
	double dx,dy,dz;
	
	C_INI_ARR(g,verts,struct vert,VBUFFER);
	C_INI_ARR(g,cells,struct cell,CBUFFER);
	C_INI_ARR(g,vempty,int,VBUFFER);
	C_INI_ARR(g,cempty,int,CBUFFER);
	if( DIM == 2 ) n[2] = 1;
	dx = 1.0/(double)n[0];
	dy = 1.0/(double)n[1];
	dz = 1.0/(double)n[2];
	g->mesh = new Mesh();
	g->mesh->SetDimensions(3);
	//g->new_marker = g->mesh->CreateTag("New",DATA_INTEGER,NODE | EDGE | FACE | CELL,NODE | EDGE | FACE | CELL,1);
	g->init_mesh(g);
	g->materials = g->mesh->CreateTag("MATERIALS",DATA_INTEGER,NODE | EDGE | FACE,NONE);
	
	g->cell_material = g->mesh->CreateTag("MATERIAL",DATA_INTEGER,CELL,NONE,1);
	g->parent = g->mesh->CreateTag("PARENT",DATA_INTEGER,CELL,NONE,multiparent);
	//g->united = g->mesh->CreateTag("UNITED",DATA_INTEGER,CELL,NONE,1);
	g->problem = g->mesh->CreateTag("PROBLEM",DATA_INTEGER,CELL,NONE,1);
	g->face_center_node = g->mesh->CreateTag("FCN",DATA_INTEGER,NODE,NONE,1);
	g->edge_face_number = g->mesh->CreateTag("EFN",DATA_INTEGER,EDGE,NONE,1);
	//g->mesh->CreateTag("test",DATA_REAL,CELL | NODE,NONE,10);
	g->K = g->mesh->CreateTag("tensorK",DATA_REAL,CELL|FACE,NONE,9);
	g->Kvec = g->mesh->CreateTag("Kvec",DATA_REAL,CELL|FACE,NONE,3);

	g->mesh->SetTopologyCheck(PRINT_NOTIFY);// | NEED_TEST_CLOSURE | PROHIBIT_MULTIPOLYGON | PROHIBIT_MULTILINE | MARK_ON_ERROR);
	//g->mesh->SetTopologyCheck(DEGENERATE_EDGE | DEGENERATE_FACE | DEGENERATE_CELL);
	//g->mesh->SetTopologyCheck(TRIPLE_SHARED_FACE | FLATTENED_CELL | INTERLEAVED_FACES);
	g->mesh->RemTopologyCheck(THROW_EXCEPTION);
	
	g->octree_node = g->mesh->CreateMarker();

	
	Mesh::GeomParam table;
	table[CENTROID] = FACE | CELL;
	table[NORMAL] = FACE;
	table[MEASURE] = FACE | CELL;
	table[ORIENTATION] = FACE;
	
	g->mesh->PrepareGeometricData(table);
	
	
	g->n[0] = n[0];
	g->n[1] = n[1];
	g->n[2] = n[2];
	g->mgrid = (int ***)malloc(sizeof(int **)*n[0]);
	for(i = 0; i < n[0]; i++)
	{
		g->mgrid[i] = (int **)malloc(sizeof(int *)*n[1]);
		for(j = 0; j < n[1]; j++)
		{
			g->mgrid[i][j] = (int *)malloc(sizeof(int)*n[2]);
			for(k = 0; k < n[2]; k++)
			{
				g->mgrid[i][j][k] = cellNextFree(g);
				cellInit(g,g->mgrid[i][j][k],-1);
				g->cells[g->mgrid[i][j][k]].center[0] = dx*(i+0.5);
				g->cells[g->mgrid[i][j][k]].center[1] = dy*(j+0.5);
				g->cells[g->mgrid[i][j][k]].center[2] = dz*(k+0.5);
				g->cells[g->mgrid[i][j][k]].side[0] = dx;
				g->cells[g->mgrid[i][j][k]].side[1] = dy;
				g->cells[g->mgrid[i][j][k]].side[2] = dz;
			}
		}
	}
	for(i = 0; i < n[0]; i++)
		for(j = 0; j < n[1]; j++)
			for(k = 0; k < n[2]; k++)
			{
				for(l = 0; l < 1<<DIM; l++)
				{
					flag = 0; v = -1;
					for(m = 0; m < 1<<DIM; m++)
					{
						c[0] = i+(m%2)+(l%2)-1;
						c[1] = j+(m/2%2)+(l/2%2)-1;
						if( DIM == 2 )
							c[2] = 0;
						else
							c[2] = k+(m/4)+(l/4)-1;
						
						if( c[0] >= 0 && c[1] >= 0 && c[2] >= 0 && c[0] < n[0] && c[1] < n[1] && c[2] < n[2])
						{
							env[m] = g->mgrid[c[0]][c[1]][c[2]];
							r = g->cells[env[m]].vertexes[((1<<DIM)-1)-m];
							if( r == -1 )
								flag = 1;
							else if( v == -1 ) 
								v = r;
							if( v != r ) TSNH
						}
						else env[m] = -1;
					}
					if( flag )
					{
						v = vertNextFree(g);
						vertInit(g,v);
						for(m = 0; m < 1<<DIM; m++)
						{
							
							g->verts[v].env[m] = env[m];
							if( env[m] != -1 )
							{
								 g->cells[env[m]].vertexes[((1<<DIM)-1)-m] = v;
							}
						}
					}
				}
			}
	for(i = 0; i < C_CNT_ARR(g,verts); i++) if( g->verts[i].busy )
		vertCreateINMOST(g,i);
	for(i = 0; i < C_CNT_ARR(g,verts); i++) if( g->verts[i].busy )
		g->vert_init_data(g,i);
		
		
	for(i = 0; i < C_CNT_ARR(g,cells); i++) if( g->cells[i].busy )
		cellCreateINMOST(g,i);

	if( remove_orphan_elements )
	{
		for(Mesh::iteratorFace f = g->mesh->BeginFace(); f != g->mesh->EndFace(); f++)
			if( f->nbAdjElements(CELL) == 0 ) f->Delete();
			
		for(Mesh::iteratorEdge f = g->mesh->BeginEdge(); f != g->mesh->EndEdge(); f++)
			if( f->nbAdjElements(CELL) == 0 ) f->Delete();
	}
	
	
	if( multiparent == ENUMUNDEF ) gridUniteSmallElements(g);

	std::cout  << Timer() - ttt << " faces " << g->mesh->NumberOfFaces() << " cells " << g->mesh->NumberOfCells() << std::endl;
}


void gridDelete(struct grid * g)
{
	int i,j,k;
	for(i = 0; i < g->n[0]; i++)
	{
		for(j = 0; j < g->n[1]; j++)
		{
			for(k = 0; k < g->n[2]; k++)
				cellDelete(g,g->mgrid[i][j][k]);
			free(g->mgrid[i][j]);
		}
		free(g->mgrid[i]);
	}
	free(g->mgrid);
	C_KIL_ARR(g,cells);
	C_KIL_ARR(g,verts);
	C_KIL_ARR(g,cempty);
	C_KIL_ARR(g,vempty);
	delete g->mesh;
}


int cellRefine(struct grid * g, int m)
{
	int i;
	const double r = 0.4;
	double r2 = (1.0/(double)g->n[0]*1.0/(double)g->n[0]+1.0/(double)g->n[1]*1.0/(double)g->n[1])/2.0;
	int test,test2;
	if( m == -1 ) return 0;
	if( !g->cells[m].busy ) return 0;
	if( g->cells[m].leaf )
	{
		if( g->cell_should_split(g,m) )
		{
			cellSplit(g,m);
			for(i = 0; i < 1<<DIM; i++)
				cellRefine(g,g->cells[m].children[i]);
			return 1;
		}
	}
	else for(i = 0; i < 1<<DIM; i++)
		cellRefine(g,g->cells[m].children[i]);
	return 0;
}

int cellCoarse(struct grid * g, int m)
{
	int i,j,flag,v,c;
	const double r = 0.4;
	double r2 = (1.0/(double)g->n[0]*1.0/(double)g->n[0]+1.0/(double)g->n[1]*1.0/(double)g->n[1])/2.0;
	int test;
	if( m == -1 ) return 0;
	if( !g->cells[m].busy ) return 0;
	if( !g->cells[m].leaf )
	{
		flag = 1;
		for(i = 0; i < 1<<DIM; i++)
		{
			cellCoarse(g,g->cells[m].children[i]);
			flag *= g->cells[g->cells[m].children[i]].leaf;
		}
		if( flag )
		{
			if( g->cell_should_unite(g,m) && g->cells[m].level > -1 )
			{
				for(i = 0; i < 1<<DIM; i++)
					for(j = 0; j < 1<<DIM; j++)
						if( j != i && j != ((1<<DIM)-1)-i )
						{
							v = g->cells[g->cells[m].children[i]].vertexes[i];
							if( v == -1 ) continue;
							c = g->verts[v].env[j];
							if( c == -1 ) continue;
							if( g->cells[c].level > g->cells[m].level+1 && g->cells[c].parent != -1)
								if( !cellCoarse(g,g->cells[c].parent) )
									return 0;
						}
				cellUnite(g,m);
				return 1;
			}
		}
		
	}
	return 0;
}

void gridRefine(struct grid * g)
{
	int i,j,k;
	
	for(i = 0; i < g->n[0]; i++)
		for(j = 0; j < g->n[1]; j++)
			for(k = 0; k < g->n[2]; k++)
				cellRefine(g,g->mgrid[i][j][k]);
				
	
}

void gridCoarse(struct grid * g)
{
	int i,j,k;
				
	
	for(i = 0; i < g->n[0]; i++)
		for(j = 0; j < g->n[1]; j++)
			for(k = 0; k < g->n[2]; k++)
				cellCoarse(g,g->mgrid[i][j][k]);	
}



void cellUniteData(struct grid * g, int m)
{
	int i;
	if( m == -1 ) return;
	if( !g->cells[m].busy ) return;
	if( !g->cells[m].leaf )
	{
		for(i = 0; i < (1<<DIM); i++)
			cellUniteData(g,g->cells[m].children[i]);
			
		g->cell_unite_data(g,m);
	}
}

void cellSplitData(struct grid * g, int m)
{
	int i;
	if( m == -1 ) return;
	if( !g->cells[m].busy ) return;
	for(i = 0; i < (1<<DIM); i++)
		g->cell_init_data(g,g->cells[m].children[i]);
	g->cell_split_data(g,m);
}

void gridUniteData(struct grid * g)
{
	int i,j,k;
	for(i = 0; i < g->n[0]; i++)
		for(j = 0; j < g->n[1]; j++)
			for(k = 0; k < g->n[2]; k++)
				cellUniteData(g,g->mgrid[i][j][k]);

}

double time_united, time_precalc_area, time_precalc_vol, time_precalc;

void cellUniteSmallElements(struct grid * g, int m)
{
	//~ if( !g->cells[m].busy || !g->cells[m].leaf || g->cells[m].mr->empty() ) return;
	Storage::real cell_vol = g->cells[m].vol;
	MarkerType cell_visited = g->mesh->CreateMarker();
	bool restart = false;
	do
	{
		restart = false;
		
		for(int i = 0; i < g->cells[m].mr->size(); i++) if( !Cell(g->mesh,(*g->cells[m].mr)[i])->GetMarker(cell_visited) )
		{
			Cell ci = Cell(g->mesh,(*g->cells[m].mr)[i]);
			ci->SetMarker(cell_visited);
			Storage::real cv = ci->Volume();
			/*
			if( cv <= 0 || cv > cell_vol*1.5 ) 
			{
				std::cout << "warning: bad volume for cell ";
				std::cout << (*g->cells[m].mr)[i]->LocalID() << " vol " << cv << " octree vol " << cell_vol;
				std::cout << " parents ";
				Storage::integer_array ps = (*g->cells[m].mr)[i]->IntegerArray(g->parent);
				for(int j = 0; j < ps.size(); j++) std::cout << ps[j] << " ";
				std::cout << std::endl;
				(*g->cells[m].mr)[i]->Integer(g->problem) = 13;
			}
			*/
			if( cv/cell_vol < MINVOLFRAC )
			{
				//std::cout << "unite cell! " << cv << " " << cell_vol << std::endl;
				double tttt = Timer();
				ElementArray<Cell> unite(g->mesh);
				unite.push_back(ci);
				Storage::real unite_vol = cv < 0 ? 0 : cv;
				Storage::integer mat = ci->Integer(g->cell_material);
				ElementArray<Face> around = ci->getFaces();
				std::map<Cell, Storage::real> around_cell;
				for(int j = 0; j < around.size(); j++)
				{
					Cell c = ci->Neighbour(around[j]);
					if( c.isValid() )
					{
						Storage::integer_array pp = c->IntegerArray(g->parent);
						//if( abs(g->cells[m].level - g->cells[pp[0]].level) > 1 ) continue;
						if( c->Integer(g->cell_material) == mat ) around_cell[c] += around[j].Area();
					}
				}
				
				//~ for(std::map<Cell *, Storage::real>::iterator j = around_cell.begin(); j != around_cell.end(); ++j)
					//~ std::cout << j->second << std::endl;
				MarkerType skip_current = g->mesh->CreateMarker();
				MarkerType visited = g->mesh->CreateMarker();
				while(unite_vol/cell_vol < MINVOLFRAC )
				{
					std::map<Cell, Storage::real>::iterator k = around_cell.end();
					double vol = 0;///*cv < 0? 0 :*/ 1e20;
					for(std::map<Cell, Storage::real>::iterator j = around_cell.begin(); j != around_cell.end(); ++j)
						if( !j->first->GetMarker(visited) && !j->first->GetMarker(skip_current) )
						{
							//if( /*cv < 0 ? j->second >= vol :*/ j->second <= vol )
							if( j->second > vol )
							{
								vol = j->second;
								k = j;
							}
						}
					if( k == around_cell.end() ) // not found - break;
						break;
					else
					{
						unite.push_back(k->first);
						if( Cell::TestUniteCells(unite,g->octree_node) )
						{
							unite_vol += k->first->Volume();
							k->first->SetMarker(visited);
						}
						else 
						{
							unite.pop_back();
							k->first->SetMarker(skip_current);
						}
					}
				}
				for(std::map<Cell, Storage::real>::iterator j = around_cell.begin(); j != around_cell.end(); ++j) 
				{
					j->first->RemMarker(visited);
					j->first->RemMarker(skip_current);
				}
				g->mesh->ReleaseMarker(visited);
				g->mesh->ReleaseMarker(skip_current);
				if( unite.size() > 1 )
				{
					//~ std::cout << "Cells to unite: ";
					//~ for(int j = 0; j < unite.size(); j++)
						//~ std::cout << unite[j] << " " << Element::GeometricTypeName(unite[j]->GetGeometricType()) << " ";
					//~ std::cout << std::endl;
					//remove them from array we look at
					MarkerType rem = g->mesh->CreateMarker();
					for(int j = 0; j < unite.size(); j++) unite[j]->SetMarker(rem);
					for(int j = 0; j < g->cells[m].mr->size(); j++)
					{
						Cell cj = Cell(g->mesh,(*g->cells[m].mr)[j]);
						if( cj->GetMarker(rem) )
							cj->SetMarker(cell_visited);
					}
					//get all the parents that link to united cells
					std::vector<int> parents;
					for(int j = 0; j < unite.size(); j++)
					{
						Storage::integer_array pp = unite[j]->IntegerArray(g->parent);
						parents.insert(parents.end(),pp.begin(),pp.end());
					}
					std::sort(parents.begin(),parents.end());
					parents.resize(std::unique(parents.begin(),parents.end())-parents.begin());
					for(int j =0 ; j < parents.size(); j++)
					{
						cell_vector::iterator it = g->cells[parents[j]].mr->begin();
						while(it != g->cells[parents[j]].mr->end())
						{
							if( g->mesh->GetMarker(*it,rem) )
								it = g->cells[parents[j]].mr->erase(it);
							else it++;
						} 
					}
					for(int j = 0; j < unite.size(); j++) unite[j]->RemMarker(rem);
					g->mesh->ReleaseMarker(rem);
					
					std::map<Tag,Storage::real> newdata = g->cell_small_unite(unite);
					
					//for(int j = 0; j < unite.size(); j++) std::cout << "(" << unite[j]->LocalID() << ", " << unite[j]->Volume() << ") ";
					//reconstruct cell by outer faces
					Cell new_cell = Cell::UniteCells(unite,g->octree_node);

					//std::cout << "-> (" << new_cell->LocalID() << ", " << new_cell->Volume() << ") " << std::endl;
					//~ std::cout << "cells united! " << new_cell << " id " << new_cell->LocalID() << " parent " << m <<  std::endl;
					//~ std::cout << "United cell: " << new_cell << " " << Element::GeometricTypeName(new_cell->GetGeometricType()) << std::endl;
					if( new_cell.isValid() )
					{
						/*
						if( new_cell->GetGeometricType() == Element::MultiPolygon )
							new_cell->Integer(g->united) = 2;
						else
							new_cell->Integer(g->united) = 1;
						*/
						for(std::map<Tag,Storage::real>::iterator it = newdata.begin(); it != newdata.end(); it++)
							new_cell->Real(it->first) = it->second;
						new_cell->Integer(g->cell_material) = mat;
						
						//unite faces
						if( unite_faces )
						{
							std::map<Cell, std::vector<Face> > face_list;
							ElementArray<Face> faces = new_cell->getFaces();
							for(ElementArray<Face>::iterator it = faces.begin(); it != faces.end(); it++)
							{
								Cell c = new_cell->Neighbour(it->self());
								face_list[c].push_back(it->self());
							}
							for(std::map<Cell, std::vector<Face> >::iterator it = face_list.begin(); it != face_list.end(); it++)
								if( it->second.size() > 1 )
								{
									//~ if( it->first == NULL ) //boundary faces, compare normals
									{
										MarkerType fv = g->mesh->CreateMarker();
										
										bool found = false;
										do
										{
											Storage::real nrm[3];
											ElementArray<Face> unite(g->mesh);
											found = false;
											for(unsigned i = 0; i < it->second.size(); i++)
												if( !it->second[i]->GetMarker(fv) )
												{
													it->second[i]->Normal(nrm);
													it->second[i]->SetMarker(fv);
													unite.push_back(it->second[i]);
													found = true;
													break;
												}
											if( found )
											{
												for(unsigned i = 0; i < it->second.size(); i++)
													if( !it->second[i]->GetMarker(fv) )
													{
														Storage::real nrmi[3];
														it->second[i]->Normal(nrmi);
														Storage::real scal = nrm[0]*nrmi[0] + nrm[1]*nrmi[1] + nrm[2]*nrmi[2];
														if( fabs(scal) > 0.99 ) 
														{
															unite.push_back(it->second[i]);
															if( Face::TestUniteFaces(unite,g->octree_node) )
																it->second[i]->SetMarker(fv);
															else 
															{
																//std::cout << "unite kills octree node " << std::endl;
																unite.pop_back();
															}
														}
														//~ if( fabs(scal) > 1 ) 
														//std::cout << "scalar product is " << scal << std::endl;
													}
												if( unite.size() > 1 )
												{
													Storage::real unite_area = 0;
													MarkerType rem = g->mesh->CreateMarker();
													for(unsigned i = 0; i < unite.size(); i++) 
													{
														unite_area += unite[i]->Area();
														unite[i]->SetMarker(rem);
														unite[i]->RemMarker(fv);
													}
													//std::cout << "unite " << unite.size() << " faces " << std::endl;
													std::vector<Face>::iterator qt = it->second.begin();
													while(qt != it->second.end())
														if( (*qt)->GetMarker(rem) )
															qt = it->second.erase(qt);
														else ++qt;
													for(unsigned i = 0; i < unite.size(); i++) unite[i]->RemMarker(rem);
													g->mesh->ReleaseMarker(rem);
													Face ret = Face::UniteFaces(unite,g->octree_node);
													
												}
												//else std::cout << " no unite faces " << std::endl;
											}
										} while( found );
										
										for(unsigned i = 0; i < it->second.size(); i++) it->second[i]->RemMarker(fv);
										
										g->mesh->ReleaseMarker(fv);
									}
								}
							
						}
						
						Storage::integer_array pp = new_cell->IntegerArray(g->parent);
						for(int j = 0; j < parents.size(); j++)
						{
							g->cells[parents[j]].mr->push_back(new_cell->GetHandle());
							pp.push_back(parents[j]);
						}
					}
					else throw -1;
				}
				else 
				{
					Cell(g->mesh,(*g->cells[m].mr)[i])->SetMarker(cell_visited);
				}
				restart = true;
				break;
				
				time_united += Timer()-tttt;
			}
		}
		
	} while( restart );
	for(int i = 0; i < g->cells[m].mr->size(); i++) Cell(g->mesh,(*g->cells[m].mr)[i])->RemMarker(cell_visited);
	g->mesh->ReleaseMarker(cell_visited);
}

void gridUniteSmallElements(struct grid * g)
{
	
	time_united = 0;
	
	
	for(int i = 0; i < C_CNT_ARR(g,cells); i++)
	{
		if( g->cells[i].busy && g->cells[i].leaf && !g->cells[i].mr->empty() )
			cellUniteSmallElements(g,i);
	}
	
}

void gridAMR(struct grid * g, int recreate)
{
	double t = Timer(), tt = Timer();
	//~ for(int i = 0; i < C_CNT_ARR(g,cells); i++) if( g->cells[i].busy ) g->cell_init_data(g,i);
	//~ std::cout << "cell init " << Timer() - t << std::endl;
	
	//~ t = Timer();
	gridUniteData(g);
	//~ std::cout << "unite data " << Timer() - t << std::endl;
	
	t = Timer();
	if( allow_coarse )gridCoarse(g);
	
	//~ std::cout << "coarse " << Timer() - t << std::endl;
	
	t = Timer();
	if( allow_refine )gridRefine(g);
	
	//~ std::cout << "refine " << Timer() - t << std::endl;
	
	if( recreate) 
	{
		/*
		if( unite_faces )
		{
			t = Timer();
			for(int i = 0; i < C_CNT_ARR(g,cells); i++)
				if( g->cells[i].busy && !g->cells[i].mr->empty())
					cellDestroyINMOST(g,i);
			//~ std::cout << "destroy " << Timer() - t << std::endl;
		}
		*/
		
		t = Timer();
		gridRecreateINMOST(g);
		//~ std::cout << "recreate " << Timer() - t << std::endl;
		//~ for(Mesh::iteratorFace it = g->mesh->BeginFace(); it != g->mesh->EndFace(); ++it)
			//~ if( it->GetMarker(g->octree_node) )
			//~ {
				//~ Geometry::FixNormaleOrientation(&*it);
				//~ it->RemMarker(g->octree_node);
			//~ }
		t = Timer();
		
		//~ if( !Element::CheckConnectivity(g->mesh) ) throw -1;
		
		if( multiparent == ENUMUNDEF ) gridUniteSmallElements(g);
		
		//~ if( !Element::CheckConnectivity(g->mesh) ) throw -1;
		//~ std::cout << "unite " << Timer() - t << " " << time_united << std::endl;
	}
	
	//~ t = Timer();
	//~ for(int i = 0; i < C_CNT_ARR(g,cells); i++) if( g->cells[i].busy ) g->cell_destroy_data(g,i);
	//~ std::cout << "destroy " << Timer() - t << std::endl;
	
	//~ for(int i = 0; i < C_CNT_ARR(g,cells); i++) if( g->cells[i].busy ) g->cell_init_data(g,i);
	
	t = Timer();
	g->mesh->ReorderEmpty(CELL | FACE | EDGE | NODE);
	//~ std::cout << "reorder empty: " << Timer() - t << std::endl;
	
	std::cout  << Timer() - tt << " faces " << g->mesh->NumberOfFaces() << " cells " << g->mesh->NumberOfCells() << std::endl;

	{
		int negvol = 0;
		for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); ++it) 
			if( it->Volume() < 0 ) 
			{
				negvol++;
				//it->Integer(g->problem) = 1;
			}
		if( negvol ) std::cout << "negative volume: " << negvol <<std::endl;
	}
}


