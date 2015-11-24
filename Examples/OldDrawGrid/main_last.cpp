//g++ main.cpp rotate.cpp -L/usr/X11R6/lib -lX11 -lXi -lXmu -lGL -lglut -lGLU ../../INMOST.a -O5
// press space - explode mesh to see connection 
#include "../../inmost.h"
#include "rotate.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include "my_glut.h"


using namespace INMOST;
Mesh * mesh;
int interactive = 0;
double zoom = 1;
int width = 800, height = 800;
double sleft = 1e20, sright = -1e20, sbottom = 1e20, stop = -1e20, sfar = -1e20, snear = 1e20;
double shift[3] = {0,0,0};
bool perspective = false;
bool pick = false, boundary = true, planecontrol = false, clipupdate = false, clipboxupdate = false;
std::map<GeometricData,ElementType> table;

#define CLIP_NONE 0
#define CLIP_NODE 1
#define CLIP_FULL 2
#define CLIP_ENDP 3

#define CLIP_FACE_NONE      0
#define CLIP_FACE_INSIDE    1
#define CLIP_FACE_OUTSIDE   2
#define CLIP_FACE_INTERSECT 3

Storage::real p[3] = {0,0,0}, n[3] = {0,0,1};
std::vector<Element *> boundary_faces;


class face2gl
{
	double dist;
	double c[4];
	double cnt[3];
	bool flag;
	ElementType etype;
	Storage::integer id;
	std::vector<double> verts;
public:
	face2gl():verts() {etype = NONE; id = 0; dist = 0; flag = false; memset(c,0,sizeof(double)*4);}
	face2gl(const face2gl & other) :verts(other.verts) 
	{
		etype = other.etype;
		id = other.id;
		dist = other.dist; 
		cnt[0] = other.cnt[0];
		cnt[1] = other.cnt[1];
		cnt[2] = other.cnt[2];
		c[0] = other.c[0];
		c[1] = other.c[1];
		c[2] = other.c[2];
		c[3] = other.c[3];
		flag = other.flag;
	}
	face2gl & operator =(face2gl const & other) 
	{ 
		etype = other.etype;
		id = other.id;
		verts = other.verts; 
		dist = other.dist; 
		cnt[0] = other.cnt[0];
		cnt[1] = other.cnt[1];
		cnt[2] = other.cnt[2];
		c[0] = other.c[0];
		c[1] = other.c[1];
		c[2] = other.c[2];
		c[3] = other.c[3];
		flag = other.flag;
		return *this;
	}
	~face2gl() {}
	void draw() const
	{
		glColor4dv(c); 
		for(unsigned k = 0; k < verts.size(); k+=3) 
		{
			glVertex3dv(cnt);
			glVertex3dv(&verts[k]);
			glVertex3dv(&verts[(k+3)%verts.size()]);
		}
	}
	void drawedges() const
	{
		for(unsigned k = 0; k < verts.size(); k+=3) 
		{
			glVertex3dv(&verts[k]); 
			glVertex3dv(&verts[(k+3)%verts.size()]); 
		}
	}
	bool operator <(const face2gl & other) const {return dist < other.dist;}
	void set_color(double r, double g, double b, double a) {c[0] = r; c[1] = g; c[2] = b; c[3] = a;}
	void add_vert(double x, double y, double z) {unsigned s = verts.size(); verts.resize(s+3); verts[s] = x; verts[s+1] = y; verts[s+2] = z;}
	void add_vert(double v[3]) {verts.insert(verts.end(),v,v+3);}
	void set_center(double _cnt[3])
	{
		cnt[0] = _cnt[0];
		cnt[1] = _cnt[1];
		cnt[2] = _cnt[2];
	}
	void compute_center()
	{
		cnt[0] = cnt[1] = cnt[2] = 0;
		for(size_t k = 0; k < verts.size(); k+=3)
		{
			cnt[0] += verts[k+0];
			cnt[1] += verts[k+1];
			cnt[2] += verts[k+2];
		}
		cnt[0] /= (verts.size()/3)*1.0;
		cnt[1] /= (verts.size()/3)*1.0;
		cnt[2] /= (verts.size()/3)*1.0;
	}
	void compute_dist(double cam[3])
	{
		dist = sqrt((cnt[0]-cam[0])*(cnt[0]-cam[0])+(cnt[1]-cam[1])*(cnt[1]-cam[1])+(cnt[2]-cam[2])*(cnt[2]-cam[2]));
	}
	void set_flag(bool set) { flag = set;}
	bool get_flag() {return flag;}
	void set_elem(ElementType _etype, Storage::integer _id) {etype = _etype; id = _id;}
};

std::vector<face2gl> all_boundary;
std::vector<face2gl> clip_boundary;
std::vector<face2gl> clip_slice;

void draw_faces(std::vector<face2gl> & set)
{
	glBegin(GL_TRIANGLES);
	for(size_t q = 0; q < set.size() ; q++) set[q].draw();
	glEnd();
}

void draw_edges(std::vector<face2gl> & set)
{
	glBegin(GL_LINES);
	for(size_t q = 0; q < set.size() ; q++) set[q].drawedges();
	glEnd();
}

void draw_faces_interactive(std::vector<face2gl> & set)
{
	glBegin(GL_TRIANGLES);
	for(size_t q = 0; q < set.size() ; q++) if( set[q].get_flag() ) set[q].draw();
	glEnd();
}

void draw_edges_interactive(std::vector<face2gl> & set)
{
	glBegin(GL_LINES);
	for(size_t q = 0; q < set.size() ; q++) if( set[q].get_flag() ) set[q].drawedges();
	glEnd();
}

Storage::integer clip_plane_edge(double sp0[3], double sp1[3], double p[3], double n[3], double node[3])
{
	Storage::real u[3], w[3], D, N, sI;
	u[0] = sp1[0] - sp0[0]; u[1] = sp1[1] - sp0[1]; u[2] = sp1[2] - sp0[2];
	w[0] = sp0[0] - p[0];   w[1] = sp0[1] - p[1];   w[2] = sp0[2] - p[2];
	D =  (n[0]*u[0] + n[1]*u[1] + n[2]*u[2]);
	N = -(n[0]*w[0] + n[1]*w[1] + n[2]*w[2]);
	if( fabs(D) < 1.0e-9 )
	{
		if( fabs(N) < 1.0e-9 ) return CLIP_FULL;
		else return CLIP_NONE;
	}
	else
	{
		sI = N/D;
		if( sI < 0-1.0e-9 || sI > 1+1.0e-9 ) return CLIP_NONE;
		else
		{
			node[0] = sp0[0] + sI * u[0];
			node[1] = sp0[1] + sI * u[1];
			node[2] = sp0[2] + sI * u[2];
			return (sI > 1.0e-9 && sI < 1.0-1.0e-9) ? CLIP_NODE : CLIP_ENDP;
		}
	}
}

class kdtree
{
	int marked;
	struct entry
	{
		Element * e;
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
	INMOST_DATA_ENUM_TYPE size;
	Storage::real bbox[6];
	kdtree * children;
	static int cmpElements0(const void * a,const void * b) 
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[0];
		float bd = eb->xyz[0];
		return (ad > bd) - (ad < bd);
	}
	static int cmpElements1(const void * a,const void * b) 
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[1];
		float bd = eb->xyz[1];
		return (ad > bd) - (ad < bd);
	}
	static int cmpElements2(const void * a,const void * b) 
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[2];
		float bd = eb->xyz[2];
		return (ad > bd) - (ad < bd);
	}
	inline static unsigned int flip(const unsigned int * fp)
	{
		unsigned int mask = -((int)(*fp >> 31)) | 0x80000000;
		return *fp ^ mask;
	}
#define _0(x)	(x & 0x7FF)
#define _1(x)	(x >> 11 & 0x7FF)
#define _2(x)	(x >> 22 )
	void radix_sort(int dim, struct entry * temp)
	{
		unsigned int i;
		const unsigned int kHist = 2048;
		unsigned int  b0[kHist * 3];
		unsigned int *b1 = b0 + kHist;
		unsigned int *b2 = b1 + kHist;
		memset(b0,0,sizeof(unsigned int)*kHist*3);
		for (i = 0; i < size; i++) 
		{
			unsigned int fi = flip((unsigned int *)&set[i].xyz[dim]);
			++b0[_0(fi)]; ++b1[_1(fi)]; ++b2[_2(fi)];
		}
		{
			unsigned int sum0 = 0, sum1 = 0, sum2 = 0;
			for (i = 0; i < kHist; i++) 
			{
				b0[kHist-1] = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = b0[kHist-1];
				b1[kHist-1] = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = b1[kHist-1];
				b2[kHist-1] = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = b2[kHist-1];
			}
		}
		for (i = 0; i < size; i++) temp[++b0[_0(flip((unsigned int *)&set[i].xyz[dim]))]] = set[i];
		for (i = 0; i < size; i++) set[++b1[_1(flip((unsigned int *)&temp[i].xyz[dim]))]] = temp[i];
		for (i = 0; i < size; i++) temp[++b2[_2(flip((unsigned int *)&set[i].xyz[dim]))]] = set[i];
		for (i = 0; i < size; i++) set[i] = temp[i];
	}
	void kdtree_build(int dim, int & done, int total, struct entry * temp)
	{
		if( size > 1 )
		{
			if( size > 128 ) radix_sort(dim,temp); else 
			switch(dim)
			{
			case 0: qsort(set,size,sizeof(entry),cmpElements0);break;
			case 1: qsort(set,size,sizeof(entry),cmpElements1);break;
			case 2: qsort(set,size,sizeof(entry),cmpElements2);break;
			}
			children = new kdtree[2];
			children[0].set = set;
			children[0].size = size/2;
			children[1].set = set+size/2;
			children[1].size = size - size/2;
			children[0].kdtree_build((dim+1)%3,done,total,temp);
			children[1].kdtree_build((dim+1)%3,done,total,temp);
			for(int k = 0; k < 3; k++)
			{
				bbox[0+2*k] = std::min(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
				bbox[1+2*k] = std::max(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
			}
		}
		else 
		{
			assert(size == 1);
			if( set[0].e->GetElementType() == EDGE )
			{
				Storage::real_array n1 = set[0].e->getAsEdge()->getBeg()->Coords();
				Storage::real_array n2 = set[0].e->getAsEdge()->getEnd()->Coords();
				for(int k = 0; k < 3; k++)
				{
					bbox[0+2*k] = std::min(n1[k],n2[k]);
					bbox[1+2*k] = std::max(n1[k],n2[k]);
				}
				done++;
				if( done%150 == 0 )
				{
					printf("%3.1f%%\r",(done*100.0)/(total*1.0));
					fflush(stdout);
				}
			}
			else
			{
				adjacent<Node> nodes = set[0].e->getNodes();
				bbox[0] = bbox[2] = bbox[4] = 1.0e20;
				bbox[1] = bbox[3] = bbox[5] = -1.0e20;
				for(size_t k = 0; k < nodes.size(); ++k)
				{
					Storage::real_array coords = nodes[k].Coords();
					for(size_t q = 0; q < 3; q++)
					{
						bbox[q*2+0] = std::min(bbox[q*2+0],coords[q]);
						bbox[q*2+1] = std::max(bbox[q*2+1],coords[q]);
					}
				}
			}
		}
	}
	kdtree() : set(NULL), size(0), children(NULL), marked(0) {}
	inline int plane_bbox(double p[3], double n[3])
	{
		Storage::real pv[3], nv[3];
		for(int k = 0; k < 3; ++k)
		{
			if( n[k] >= 0 ) 
			{ 
				pv[k] = bbox[1+2*k]; //max
				nv[k] = bbox[0+2*k]; //min
			} 
			else 
			{ 
				pv[k] = bbox[0+2*k]; //min
				nv[k] = bbox[1+2*k]; //max
			}
		}
		Storage::real pvD, nvD;
		pvD = n[0]*(pv[0]-p[0])+n[1]*(pv[1]-p[1])+n[2]*(pv[2]-p[2]);
		nvD = n[0]*(nv[0]-p[0])+n[1]*(nv[1]-p[1])+n[2]*(nv[2]-p[2]);
		if( nvD*pvD <= 0.0 )
			return 2;
		else if( nvD < 0.0 )
			return 1;
		else return 0;
	}
	bool sub_intersect_plane_edge(Tag clip_point, Tag clip_state, std::vector<Cell *> & cells, MIDType mrk, double p[3], double n[3])
	{
		if( size == 1 )
		{
			assert( set[0].e->GetElementType() == EDGE );
			Storage::real_array sp0 = set[0].e->getAsEdge()->getBeg()->Coords();
			Storage::real_array sp1 = set[0].e->getAsEdge()->getEnd()->Coords();
			Storage::integer & clip = set[0].e->IntegerDF(clip_state);
			clip = clip_plane_edge(&sp0[0],&sp1[0],p,n,&set[0].e->RealArrayDF(clip_point)[0]);
			if( clip )
			{
				adjacent<Cell> ecells = set[0].e->getCells();
				for(size_t k = 0; k < ecells.size(); ++k) if( !ecells[k].GetMarker(mrk) )
				{
					ecells[k].SetMarker(mrk);
					cells.push_back(&ecells[k]);
				}
				marked = 1;
			}
		}
		else if( plane_bbox(p,n) == 2 )
		{
			bool test1 = children[0].sub_intersect_plane_edge(clip_point,clip_state,cells,mrk,p,n);
			bool test2 = children[1].sub_intersect_plane_edge(clip_point,clip_state,cells,mrk,p,n);
			if( test1 || test2 ) marked = 1;
		}
		return marked;
	}
	void sub_intersect_plane_faces(std::vector<face2gl> & out, double p[3], double n[3])
	{
		if( size == 1 )
		{
			int state;
			assert( set[0].e->GetElementDimension() == 2 );
			adjacent<Node> nodes = set[0].e->getNodes();
			dynarray<bool,64> nodepos(nodes.size());
			Storage::real_array coords = nodes[0].Coords();
			Storage::real dot0 = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
			if( dot0 < 0.0 ) state = CLIP_FACE_INSIDE; else state = CLIP_FACE_OUTSIDE;
			nodepos[0] = dot0 < 1.0e-10;
			for(size_t k = 0; k < nodes.size(); k++)
			{
				coords = nodes[k].Coords();
				Storage::real dot = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
				nodepos[k] = dot < 1.0e-10;
				if( dot*dot0 < 0.0 ) state = CLIP_FACE_INTERSECT;
			}
			if( state == CLIP_FACE_INTERSECT )
			{
				face2gl f;
				f.set_color(0.6,0.6,0.6,1);
				for(size_t q = 0; q < nodes.size(); q++)
				{
					if( nodepos[q] ) f.add_vert(&nodes[q].Coords()[0]);
					if( nodepos[q] != nodepos[(q+1)%nodes.size()] )
					{
						Storage::real_array sp0 = nodes[q].Coords();
						Storage::real_array sp1 = nodes[(q+1)%nodes.size()].Coords();
						Storage::real node[3];
						if( clip_plane_edge(&sp0[0],&sp1[0],p,n,node) > CLIP_NONE) f.add_vert(node);
					}
				}
				f.compute_center();
				f.set_elem(set[0].e->GetElementType(),set[0].e->LocalID());
				out.push_back(f);
			}
			else if( state == CLIP_FACE_INSIDE )
			{
				face2gl f;
				f.set_color(0.6,0.6,0.6,1);
				for(size_t q = 0; q < nodes.size(); q++) f.add_vert(&nodes[q].Coords()[0]);
				f.compute_center();
				f.set_elem(set[0].e->GetElementType(),set[0].e->LocalID());
				out.push_back(f);
			}
		}
		else
		{
			marked = plane_bbox(p,n);
			if( marked == 1 )
			{
				for(INMOST_DATA_ENUM_TYPE k = 0; k < size; k++) 
				{
					adjacent<Node> nodes = set[k].e->getNodes();
					face2gl f;
					f.set_color(0.6,0.6,0.6,1);
					for(size_t q = 0; q < nodes.size(); q++) f.add_vert(&nodes[q].Coords()[0]);
					f.compute_center();
					f.set_elem(set[0].e->GetElementType(),set[0].e->LocalID());
					out.push_back(f);
				}
			}
			else if( marked == 2 )
			{
				children[0].sub_intersect_plane_faces(out,p,n);
				children[1].sub_intersect_plane_faces(out,p,n);
			}
		}
	}
	void unmark_old_edges(Tag clip_state)
	{
		if( size == 1 )
		{
			assert(set[0].e->GetElementType() == EDGE);
			marked = 0;
			if( set[0].e->GetElementType() == EDGE )
				set[0].e->IntegerDF(clip_state) = CLIP_NONE;
			else if( set[0].e->GetElementType() == FACE )
				set[0].e->IntegerDF(clip_state) = CLIP_FACE_NONE;
		}
		else if( children )
		{
			if(children[0].marked) {children[0].unmark_old_edges(clip_state); marked = 0;}
			if(children[1].marked) {children[1].unmark_old_edges(clip_state); marked = 0;}
		}
	}
	void unmark_old_faces(Tag clip_points, Tag clip_state)
	{
		if( size == 1 )
		{
			assert(set[0].e->GetElementType() == FACE);
			if( set[0].e->IntegerDF(clip_state) == CLIP_FACE_INTERSECT )
				set[0].e->DelData(clip_points);
			set[0].e->IntegerDF(clip_state) == CLIP_FACE_NONE;
		}
		else if( children )
		{
			if(children[0].marked == 2 ) children[0].unmark_old_faces(clip_points,clip_state);
			else for(INMOST_DATA_ENUM_TYPE k = 0; k < size; k++) set[k].e->IntegerDF(clip_state) = CLIP_FACE_NONE;

			if(children[1].marked == 2 ) children[1].unmark_old_faces(clip_points,clip_state);
			else for(INMOST_DATA_ENUM_TYPE k = 0; k < size; k++) set[k].e->IntegerDF(clip_state) = CLIP_FACE_NONE;
		}
	}
	void clear_children() { if( children ) {children[0].clear_children(); children[1].clear_children();}}
public:
	kdtree(Mesh * m) : children(NULL), marked(0)
	{
		double tt;
		size = m->NumberOfEdges();
		assert(size > 1);
		set = new entry[size];
		INMOST_DATA_ENUM_TYPE k = 0;
		tt = Timer();
		printf("Prepearing edge set.\n");
		for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it) 
		{
			set[k].e = &*it;
			set[k].xyz[0] = (it->getBeg()->Coords()[0] + it->getEnd()->Coords()[0])*0.5;
			set[k].xyz[1] = (it->getBeg()->Coords()[1] + it->getEnd()->Coords()[1])*0.5;
			set[k].xyz[2] = (it->getBeg()->Coords()[2] + it->getEnd()->Coords()[2])*0.5;
			k++;
			if( k%150 == 0 ) 
			{
				printf("%3.1f%%\r",(k*100.0)/(size*1.0));
				fflush(stdout);
			}
		}
		printf("Done. Time %g\n",Timer()-tt);
		int done = 0, total = size;
		printf("Building KD-tree.\n");
		tt = Timer();
		struct entry *  temp = new entry[size];
		kdtree_build(0,done,total,temp);
		delete [] temp;
		printf("Done. Time %g\n",Timer()-tt);
		for(int k = 0; k < 3; k++)
		{
			bbox[0+2*k] = std::min(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
			bbox[1+2*k] = std::max(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
		}
	}
	kdtree(Mesh * m, Element ** eset, INMOST_DATA_ENUM_TYPE size) : children(NULL), marked(0), size(size)
	{
		double tt;
		assert(size > 1);
		set = new entry[size];
		tt = Timer();
		printf("Prepearing elements set.\n");
		for(INMOST_DATA_ENUM_TYPE k = 0; k < size; k++) 
		{
			set[k].e = eset[k];
			Storage::real cnt[3];
			set[k].e->Centroid(cnt);
			set[k].xyz[0] = cnt[0];
			set[k].xyz[1] = cnt[1];
			set[k].xyz[2] = cnt[2];
			if( k%150 == 0 ) 
			{
				printf("%3.1f%%\r",(k*100.0)/(size*1.0));
				fflush(stdout);
			}
		}
		printf("Done. Time %g\n",Timer()-tt);
		int done = 0, total = size;
		printf("Building KD-tree.\n");
		tt = Timer();
		struct entry *  temp = new entry[size];
		kdtree_build(0,done,total,temp);
		delete [] temp;
		printf("Done. Time %g\n",Timer()-tt);
		for(int k = 0; k < 3; k++)
		{
			bbox[0+2*k] = std::min(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
			bbox[1+2*k] = std::max(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
		}
	}
	void intersect_plane_edge(Tag clip_point, Tag clip_state, std::vector<Cell*> & cells, MIDType mark_cells, double p[3], double n[3])
	{
		if( marked ) 
		{
			unmark_old_edges(clip_state);
			cells.clear();
		}
		sub_intersect_plane_edge(clip_point, clip_state, cells,mark_cells,p,n);
	}
	void intersect_plane_face(std::vector<face2gl> & out, double p[3], double n[3])
	{
		double t = Timer();
		t = Timer();
		sub_intersect_plane_faces(out, p,n);
		printf("intersect %g\n",Timer()-t);
	}
	~kdtree()
	{
		delete [] set;
		clear_children();
	}
};

class clipper
{
	struct edge_point
	{
		Storage::real xyz[3];
		Storage::integer edge;
		edge_point(){}
		edge_point(Storage::real _xyz[3], Storage::integer n)
		{
			xyz[0] = _xyz[0];
			xyz[1] = _xyz[1];
			xyz[2] = _xyz[2];
			edge = n;
		}
		bool operator ==(const edge_point& b) const
		{
			Storage::real temp = 0.0;
			for(int k = 0; k < 3; k++) temp += (xyz[k]-b.xyz[k])*(xyz[k]-b.xyz[k]);
			if( temp < 1.0e-8 ) return true; else return false;
		}
		bool operator !=(const edge_point& b) const {return !(operator ==(b));}
		void print() {printf("%g %g %g e %d\n",xyz[0],xyz[1],xyz[2],edge);}
	};
	Tag clip_point, clip_state;
	kdtree * tree;
	Tag clips;
	MIDType marker;
	std::vector<Cell *> cells;
	Mesh * mm;
public:
	~clipper() 
	{
		delete tree; 
		for(size_t k = 0; k < cells.size(); k++) cells[k]->RemMarker(marker);
		mm->ReleaseMarker(marker); 
		mm->DeleteTag(clips); 
		mm->DeleteTag(clip_point); 
		mm->DeleteTag(clip_state);
	}
	clipper(Mesh * m)
	{
		mm = m;
		tree = new kdtree(m);
		marker = m->CreateMarker();
		clips = m->CreateTag("CLIPS",DATA_REAL,CELL,CELL);
		clip_point = m->CreateTag("CLIP_POINT",DATA_REAL,EDGE,NONE,3);
		clip_state = m->CreateTag("CLIP_STATE",DATA_INTEGER,EDGE,NONE,1);
	}
	void clip_plane(Storage::real p[3], Storage::real n[3])
	{
		const bool print = false;
		double t1,t2;
		t1 = Timer();
		for(size_t k = 0; k < cells.size(); ++k) 
			if( cells[k]->GetMarker(marker) )
			{
				cells[k]->DelData(clips);
				cells[k]->RemMarker(marker);
			}
		tree->intersect_plane_edge(clip_point,clip_state,cells,marker,p,n);
		t1 = Timer() - t1;

		t2 = Timer();
		adjacent<Face> faces;
		adjacent<Edge> edges;
		dynarray<edge_point,128> clipcoords, loopcoords;
		std::vector<bool> closed;
		for(size_t k = 0; k < cells.size(); ++k)
		{
			//assuming faces are convex we will have at most one clipping edge per polygon
			//otherwise every pair of clipping nodes forming edge should appear consequently
			//as long as we go through face's edges in ordered way
			clipcoords.clear();
			//we will gather all the pairs of nodes, then form closed loop
			faces = cells[k]->getFaces();
			Face * full_face = NULL;
			int ntotpoints = 0, ntotedges = 0;
			for(size_t q = 0; q < faces.size(); ++q)
			{
				int last_edge_type = CLIP_NONE;
				int nfulledges = 0, npoints = 0, nstartedge = ntotedges;
				edges = faces[q].getEdges();
				for(size_t r = 0; r < edges.size(); ++r)
				{
					Storage::integer state = edges[r].IntegerDF(clip_state);
					if( state == CLIP_FULL )
					{
						nfulledges++;
						edge_point n1 = edge_point(&edges[r].getBeg()->Coords()[0],ntotedges);
						edge_point n2 = edge_point(&edges[r].getEnd()->Coords()[0],ntotedges);
						if( npoints % 2 == 0 ) //all privious edges are closed, just add this one
						{
							clipcoords.push_back(n1);
							clipcoords.push_back(n2);
							npoints+=2;
							ntotedges++;
							last_edge_type = CLIP_FULL;
						} 
						else if(n1 == clipcoords.back()) //this may be prolongation of one point that hit one edge
						{
							clipcoords.push_back(n2);
							npoints++;
							ntotedges++;
							last_edge_type = CLIP_FULL;
						}
						else if( n2 == clipcoords.back() )
						{
							clipcoords.push_back(n1);
							npoints++;
							ntotedges++;
							last_edge_type = CLIP_FULL;
						}
						else printf("%s:%d strange orphan node before me\n",__FILE__,__LINE__);
					}
					else if( state == CLIP_ENDP )
					{
						edge_point n = edge_point(&edges[r].RealArrayDF(clip_point)[0],ntotedges);
						bool add = true;
						if( last_edge_type == CLIP_ENDP )
						{
							if( n == clipcoords.back() )
								add = false;
						}
						else if( last_edge_type == CLIP_FULL )
						{
							if( n == clipcoords.back() || n == clipcoords[clipcoords.size()-2])
								add = false;
						}
						if( add ) //this one node should be prolongation of privious edge
						{
							if( print )
							{
								printf("added: ");
								n.print();
							}
							clipcoords.push_back(n);
							npoints++;
							if( npoints % 2 == 0 ) 
							{
								if( print ) printf("edge %d accepted\n",ntotedges);
								ntotedges++;
							}
							last_edge_type = CLIP_ENDP;
						}
						else if( print )
						{
							printf("ignored: ");
							n.print();
						}
					}
					else if( state == CLIP_NODE )
					{
						edge_point n = edge_point(&edges[r].RealArrayDF(clip_point)[0],ntotedges);
						if( print )
						{
							printf("added: ");
							n.print();
						}
						clipcoords.push_back(n);
						npoints++;
						if( npoints % 2 == 0 ) 
						{
							if( print ) printf("edge %d accepted\n",ntotedges);
							ntotedges++;
						}
						last_edge_type = CLIP_NODE;
					}
				}
				if( npoints % 2 != 0 ) 
				{
					if( print ) printf("edge %d not closed - remove\n",ntotedges);
					clipcoords.pop_back();
					npoints--;
					//printf("%s:%d this should not happen!\n",__FILE__,__LINE__);
				}
				
				if( nfulledges == edges.size() )
				{
					full_face = &faces[q];
					break;
				}
				if( print )
				{
					printf("nodes on face %d\n",faces[q].LocalID());
					for(int m = nstartedge*2; m < clipcoords.size(); m++) clipcoords[m].print();
				}
				ntotpoints += npoints;
			}
			if( full_face )
			{
				adjacent<Node> nodes = full_face->getNodes();
				Storage::real_array cl = cells[k]->RealArray(clips);
				cl.resize(3*nodes.size());
				for(size_t r = 0; r < nodes.size(); r++)
				{
					Storage::real_array p = nodes[r].Coords();
					cl[0+3*r] = p[0];
					cl[1+3*r] = p[1];
					cl[2+3*r] = p[2];
				}
				cells[k]->SetMarker(marker);
			}
			else if( ntotedges > 2 )
			{
				if( print )
				{
					printf("coords on cell %d\n",cells[k]->LocalID());
					for(int m = 0; m < clipcoords.size(); m++) clipcoords[m].print();
				}
				//Can make this faster using hash
				closed.resize(ntotedges);
				std::fill(closed.begin(),closed.end(),false);
				loopcoords.push_back(clipcoords[0]); //this is starting point
				loopcoords.push_back(clipcoords[1]); //this is next
				closed[0] = true;
				for(int r = 0; r < ntotedges-2; ++r) //we need to add this number of points
				{
					bool hit = false;
					for(int q = 0; q < ntotedges; ++q) if( !closed[q] )
					{
						//some end of q-th edge connects to current end point - connect it
						if( clipcoords[q*2+0] == loopcoords.back() )
						{
							loopcoords.push_back(clipcoords[q*2+1]);
							closed[q] = true;
							hit = true;
							break;
						}
						else if( clipcoords[q*2+1] == loopcoords.back() )
						{
							loopcoords.push_back(clipcoords[q*2+0]);
							closed[q] = true;
							hit = true;
							break;
						}
					}
					if( !hit ) printf("%s:%d cannot find end for edge! total edges %d current loop size %d\n",
						__FILE__,__LINE__,ntotedges,loopcoords.size());
				}
				Storage::real_array cl = cells[k]->RealArray(clips);
				cl.resize(3*loopcoords.size());
				for(size_t r = 0; r < loopcoords.size(); ++r)
				{
					cl[r*3+0] = loopcoords[r].xyz[0];
					cl[r*3+1] = loopcoords[r].xyz[1];
					cl[r*3+2] = loopcoords[r].xyz[2];
				}

				loopcoords.clear();
				clipcoords.clear();

				cells[k]->SetMarker(marker);
			}
		}
		t2 = Timer() - t2;
		//printf("intersect edges %20g compute faces %20g\n",t1,t2);
		//printf("Total cells: %d marked: %d\n",cells.size(),numcells);
	}
	void gen_clip(std::vector<face2gl> & out)
	{
		size_t pace = std::max<size_t>(1,std::min<size_t>(15,size()/100));
		for(size_t k = 0; k < cells.size(); k++) if( cells[k]->GetMarker(marker))
		{
			face2gl f;
			f.set_color(0.6,0.6,0.6,1);
			Storage::real_array cl = cells[k]->RealArray(clips);
			for(size_t q = 0; q < cl.size(); q+=3) f.add_vert(&cl[q]);
			f.compute_center();
			f.set_elem(cells[k]->GetElementType(),cells[k]->LocalID());
			if( k%pace == 0 ) f.set_flag(true);
			out.push_back(f);
		}
	}
	void draw_clip(size_t pace)
	{
		for(size_t k = 0; k < cells.size(); k+=pace) if( cells[k]->GetMarker(marker))
		{
			Storage::real_array cl = cells[k]->RealArray(clips);
			glBegin(GL_POLYGON);
			for(size_t q = 0; q < cl.size(); q+=3) glVertex3dv(&cl[q]);
			glEnd();
		}
	}
	void draw_clip_edges(size_t pace)
	{
		for(size_t k = 0; k < cells.size(); k+=pace) if( cells[k]->GetMarker(marker))
		{
			Storage::real_array cl = cells[k]->RealArray(clips);
			glBegin(GL_LINE_LOOP);
			for(size_t q = 0; q < cl.size(); q+=3) glVertex3dv(&cl[q]);
			glEnd();
		}
	}
	size_t size() {return cells.size();}
} * oclipper = NULL;

class bnd_clipper
{
	Tag clip_points, clip_state;
	kdtree * tree;
	Mesh * mm;
	Element ** faces;
	INMOST_DATA_ENUM_TYPE nfaces;
public:
	~bnd_clipper()
	{
		mm->DeleteTag(clip_points);
		mm->DeleteTag(clip_state);
		delete tree;
		delete [ ]faces;
	}
	bnd_clipper(Mesh * m , Element ** _faces, INMOST_DATA_ENUM_TYPE size)
	{
		mm = m;
		clip_points = m->CreateTag("CLIP_FACE_POINTS",DATA_REAL,FACE,FACE);
		clip_state = m->CreateTag("CLIP_FACE_STATE",DATA_INTEGER,FACE,NONE,1);
		nfaces = size;
		faces = new Element *[nfaces];
		for(INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k++) faces[k] = _faces[k];
		tree = new kdtree(mm,faces,nfaces);
	}
	void clip_plane(std::vector<face2gl> & out, Storage::real p[3], Storage::real n[3])
	{
		out.clear();
		tree->intersect_plane_face(out,p,n);
	}
	/*
	void draw_clip(size_t pace)
	{
		for(INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k+=pace)
		{
			int state = faces[k]->IntegerDF(clip_state);
			if( state == CLIP_FACE_INSIDE )
			{
				adjacent<Node> nodes = faces[k]->getNodes();
				glBegin(GL_POLYGON);
				for(size_t q = 0; q < nodes.size(); q++) glVertex3dv(&nodes[q].Coords()[0]);
				glEnd();
			}
			else if( state == CLIP_FACE_INTERSECT )
			{
				Storage::real_array arr = faces[k]->RealArray(clip_points);
				glBegin(GL_POLYGON);
				for(size_t q = 0; q < arr.size(); q+=3) glVertex3dv(&arr[q]);
				glEnd();
			}
		}
	}
	void draw_clip_edges(size_t pace)
	{
		for(INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k+=pace)
		{
			int state = faces[k]->IntegerDF(clip_state);
			if( state == CLIP_FACE_INSIDE )
			{
				adjacent<Node> nodes = faces[k]->getNodes();
				glBegin(GL_LINE_LOOP);
				for(size_t q = 0; q < nodes.size(); q++) glVertex3dv(&nodes[q].Coords()[0]);
				glEnd();
			}
			else if( state == CLIP_FACE_INTERSECT )
			{
				Storage::real_array arr = faces[k]->RealArray(clip_points);
				glBegin(GL_LINE_LOOP);
				for(size_t q = 0; q < arr.size(); q+=3) glVertex3dv(&arr[q]);
				glEnd();
			}
		}
	}
	*/
} * bclipper = NULL;



void set_matrix3d()
{
	double aspect = (double)width/(double)height;
	double side = std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
	//double center[3] = { (sleft+sright)*0.5, (sbottom+stop)*0.5, (sfar+snear)*0.5};
	const double sc = 2;
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	//~ glOrtho(center[0]-sc*side*zoom*aspect,center[0]+sc*zoom*side*aspect,
		//~ center[1]-sc*side*zoom,center[1]+sc*zoom*side,
		//~ center[2]-sc*side*100,center[2]+sc*side*100);
	if( !perspective )
	{
		glOrtho(-sc*side*zoom*aspect,sc*side*zoom*aspect,
				-sc*side*zoom,sc*side*zoom,
				-sc*side*100,sc*side*100);
	}
	else
		gluPerspective(60.0, aspect, 0.00001, 1000.0);
	glMatrixMode (GL_MODELVIEW);
}

void set_matrix2d()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1,1,-1,1,-1,1);
	glMatrixMode(GL_MODELVIEW);
}

void reshape(int w, int h)
{
	width = w;
	height = h;
	set_matrix3d();
	glViewport(0, 0, w, h);
}


int actionstate  = 0;
double mymx = 0;
double mymy = 0;

void myclickmotion(int nmx, int nmy) // Mouse
{
	double lmx = 2.*(nmx/(double)width - 0.5),lmy = 2.*(0.5 - nmy/(double)height), dmx = lmx-mymx, dmy = lmy - mymy;
	if( actionstate == 1 ) //middle button
	{
		double shiftmod[3] = {0,0,0};
		shiftmod[0] += dmx*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear);
		shiftmod[1] += dmy*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear);
		if( planecontrol )
		{
			
			rotatevector_from_stack((double*)shiftmod);
			p[0] += shiftmod[0];
			p[1] += shiftmod[1];
			p[2] += shiftmod[2];
			clipupdate = true;
		}
		else
		{
			//shiftmod[0] += dmx*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
			//shiftmod[1] += dmy*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
			rotatevector((double*)shiftmod);
			shift[0] += shiftmod[0];
			shift[1] += shiftmod[1];
			shift[2] += shiftmod[2];
		}
		glutPostRedisplay();
		mymx = lmx;
		mymy = lmy;
	}
	else if( actionstate == 2 ) //right button
	{
		if( planecontrol )
		{
			double shiftmod[3] = {0,0,0};
			shiftmod[2] -= dmy*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear);
			rotatevector_from_stack((double*)shiftmod);
			p[0] += shiftmod[0];
			p[1] += shiftmod[1];
			p[2] += shiftmod[2];
			clipupdate = true;
		}
		else
		{
			zoom *= expf(-dmy);
			reshape(width,height);
			double shiftmod[3] = {0,0,0};
			shiftmod[2] += dmx*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear);
			rotatevector((double*)shiftmod);
			shift[0] += shiftmod[0];
			shift[1] += shiftmod[1];
			shift[2] += shiftmod[2];
		}
		glutPostRedisplay();
		mymx = lmx;
		mymy = lmy;
	}
	else if( actionstate == 3 ) //left buttion
	{
		clickmotion(nmx,nmy);
		if( planecontrol )
		{
			rotatevector((double *)n);
			clipupdate = true;
			quatinit();
		}
	}
}
void mymotion(int nmx, int nmy) // Mouse
{
	motion(nmx,nmy);
	mymx = 2.*(nmx/(double)width - 0.5);
	mymy = 2.*(0.5 - nmy/(double)height);
}

void myclick(int b, int s, int nmx, int nmy) // Mouse
{
	if( b == GLUT_LEFT_BUTTON )
	{
		if( s == GLUT_DOWN )
		{
			actionstate = 3;
		}
		else
		{
			actionstate = 0;
		}
		click(b,s,nmx,nmy);
	}
	else if( b == GLUT_MIDDLE_BUTTON )
	{
		if( s == GLUT_DOWN )
		{
			actionstate = 1;
			interactive = true;
		}
		else
		{
			actionstate = 0;
			interactive = false;
		}
		mymx = 2.*(nmx/(double)width - 0.5);
		mymy = 2.*(0.5 - nmy/(double)height);
	}
	else if( b == GLUT_RIGHT_BUTTON )
	{
		if( s == GLUT_DOWN )
		{
			actionstate = 2;
			interactive = true;
		}
		else
		{
			actionstate = 0;
			interactive = false;
		}
		mymx = 2.*(nmx/(double)width - 0.5);
		mymy = 2.*(0.5 - nmy/(double)height);
	}
	glutPostRedisplay();
}



void keyboard(unsigned char key, int x, int y)
{
	if( key == 27 )
	{
		if( oclipper ) delete oclipper;
		if( bclipper ) delete bclipper;
		delete mesh;
		exit(-1);
	}
	else if( key == 'p' )
	{
		pick = !pick;
		glutPostRedisplay();
	}
	
	else if( key == '=' || key == '+')
	{
		zoom /= 1.1;
		reshape(width,height);
		glutPostRedisplay();
		interactive = true;
		//reset_timer = Timer();
	}
	else if( key == '_' || key == '-')
	{
		zoom *= 1.1;	
		reshape(width,height);
		glutPostRedisplay();
		interactive = true;
		//reset_timer = Timer();
	}
	else if( key == 'w' )
	{
		double shiftmod[3] = {0,0,0};
		shiftmod[1] -= 0.03f*expf(zoom-1)*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
		rotatevector((double*)shiftmod);
		shift[0] += shiftmod[0];
		shift[1] += shiftmod[1];
		shift[2] += shiftmod[2];
		glutPostRedisplay();
		interactive = true;
	}
	else if( key == 's' )
	{
		if( !planecontrol )
		{
			double shiftmod[3] = {0,0,0};
			shiftmod[1] += 0.03f*expf(zoom-1)*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
			rotatevector((double*)shiftmod);
			shift[0] += shiftmod[0];
			shift[1] += shiftmod[1];
			shift[2] += shiftmod[2];
			glutPostRedisplay();
			interactive = true;
		}
	}
	else if( key == 'a' )
	{
		if( !planecontrol )
		{
			double shiftmod[3] = {0,0,0};
			shiftmod[0] += 0.03f*expf(zoom-1)*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
			rotatevector((double*)shiftmod);
			shift[0] += shiftmod[0];
			shift[1] += shiftmod[1];
			shift[2] += shiftmod[2];
			glutPostRedisplay();
			interactive = true;
		}
	}
	else if( key == 'd' )
	{
		if( !planecontrol )
		{
			double shiftmod[3] = {0,0,0};
			shiftmod[0] -= 0.03f*expf(zoom-1)*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
			rotatevector((double*)shiftmod);
			shift[0] += shiftmod[0];
			shift[1] += shiftmod[1];
			shift[2] += shiftmod[2];
			glutPostRedisplay();
			interactive = true;
		}
	}
	else if( key == 'r' )
	{
		if( planecontrol )
		{
			p[0] = (sleft+sright)*0.5;
			p[1] = (sbottom+stop)*0.5;
			p[2] = (sfar+snear)*0.5;
			n[0] = 0;
			n[1] = 0;
			n[2] = 1;
			quatinit();
			clipupdate = true;
		}
		else
		{
			shift[0] = 0.0f;
			shift[1] = 0.0f;
			shift[2] = 0.0f;
			zoom = 1;
			quatinit();
		}
		glutPostRedisplay();
		interactive = true;
	}
	else if( key == 'z' )
	{
		printf("Enter point of plane (x,y,z):\n");
		scanf("%lf %lf %lf",p,p+1,p+2);
		printf("Enter normal of plane (x,y,z):\n");
		scanf("%lf %lf %lf",n,n+1,n+2);
		Storage::real l = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
		if( l )
		{
			n[0] /= l;
			n[1] /= l;
			n[2] /= l;
		}
		if( oclipper ) 
		{
			oclipper->clip_plane(p,n);
			if( bclipper ) bclipper->clip_plane(clip_boundary,p,n);
			glutPostRedisplay();
		}
	}
	else if( key == 'b' )
	{
		boundary = !boundary;
		glutPostRedisplay();
	}
	else if( key == 'l' )
	{
		planecontrol = !planecontrol;
		if( planecontrol ) 
		{
			quatpush();
			quatinit();
		}
		else quatpop();
	}
}

void keyboard2(unsigned char key, int x, int y)
{
	if( key == '=' || key == '+' ||  key == '_' || key == '-' || key == 'w' || key == 's' || key == 'a' || key == 'd' || key == 'r' || key == 'p' || key == 'z')
	{
		interactive = false;
		glutPostRedisplay();
	}
	
}

face2gl DrawFace(Element * f)
{
	double cnt[3];
	face2gl ret;
	f->Centroid(cnt);
	ret.set_center(cnt);
	adjacent<Node> nodes = f->getNodes();

	if( f->nbAdjElements(CELL) == 0 ) ret.set_color(1,0,0,0.1);
	else ret.set_color(0,1,0,0.1);
	
	for( adjacent<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
		ret.add_vert(&(kt->Coords()[0]));
	ret.set_elem(f->GetElementType(),f->LocalID());
	return ret;
}

void whereami(double & cx, double & cy, double & cz)
{
   // Get the viewing matrix
   GLdouble modelview[16],projection[16];
   GLint viewport[4] = {0,0,1,1};
   glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
   glGetDoublev(GL_PROJECTION_MATRIX, projection);
   
   GLdouble outx, outy, outz;  // Var's to save the answer in

   gluUnProject(0.5, 0.5, 0.,
               modelview, projection, viewport,
               &outx, &outy, &outz);

   // Return the result.
   cx = outx;
   cy = outy;
   cz = outz;
}


void pick_mouse(double origin[3], double direction[3])
{
   // Get the viewing matrix
   GLdouble modelview[16],projection[16];
   GLint viewport[4] = {0,0,1,1};
   glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
   glGetDoublev(GL_PROJECTION_MATRIX, projection);
   
   GLdouble outx, outy, outz;  // Var's to save the answer in

   gluUnProject((mymx+1.0)/2.0, (mymy+1.0)/2.0, 0.,
               modelview, projection, viewport,
               &outx, &outy, &outz);

   // Return the result.
   origin[0] = outx;
   origin[1] = outy;
   origin[2] = outz;
   
   gluUnProject((mymx+1.0)/2.0, (mymy+1.0)/2.0, 1.,
               modelview, projection, viewport,
               &outx, &outy, &outz);
   direction[0] = outx - origin[0];
   direction[1] = outy - origin[1];
   direction[2] = outz - origin[2];
}

void draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	set_matrix3d();

	Storage::real mult = zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5*0.1;
	if( perspective )
		glTranslated(0,0,-zoom*2*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5);
	if( planecontrol ) 
		rotate_from_stack(); 
	else 
		rotate();

	//axis
	{
		glLineWidth(3.0);
		glBegin(GL_LINES);
		glColor3f(1,0,0);
		glVertex3d(0,0,0);
		glVertex3d(mult,0,0);
		glColor3f(0,1,0);
		glVertex3d(0,0,0);
		glVertex3d(0,mult,0);
		glColor3f(0,0,1);
		glVertex3d(0,0,0);
		glVertex3d(0,0,mult);
		glEnd();
		glLineWidth(1.0);
	}
	//glPointSize(1);
	
	//glTranslated(-(sleft+sright)*0.5+shift[0],-(sbottom+stop)*0.5 + shift[1],-(snear+sfar)*0.5 + shift[2]);
	
	glTranslated(shift[0],shift[1],shift[2]);


	//if( planecontrol )
	{
		glColor3f(0.6,0.4,0.2);
		glPointSize(5);
		glBegin(GL_POINTS);
		glVertex3dv(p);
		glEnd();
		glPointSize(1);

		glColor3f(0.2,0.4,0.6);
		glLineWidth(3.0);
		glBegin(GL_LINES);
		glVertex3dv(p);
		glVertex3d(p[0]+n[0]*mult*2,p[1]+n[1]*mult*2,p[2]+n[2]*mult*2);
		glEnd();
		glLineWidth(1.0);
	}

	

	double campos[3] = {0.5,0.5,0};
	whereami(campos[0],campos[1],campos[2]);

	//glTranslated((l+r)*0.5,(b+t)*0.5,(near+far)*0.5);

	{
		if( interactive && !bclipper)
		{
			std::vector<bool> nodepos;
			std::vector<Storage::real> faceverts;
			size_t pace = interactive ? std::max<size_t>(1,std::min<size_t>(25,boundary_faces.size()/100)) : 1;
			for(size_t k = 0; k < boundary_faces.size(); k+=pace)
			{
				int state = 0; //0 - outside, 1 - inside, 2 - intersect
				adjacent<Node> nodes = boundary_faces[k]->getNodes();
				Storage::real_array coords = nodes[0].Coords();
				Storage::real dot0 = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
				if( dot0 < 0.0 ) state = 1; else state = 0;
				nodepos.resize(nodes.size());
				nodepos[0] = dot0 < 1.0e-10;
				for(size_t q = 1; q < nodes.size(); ++q)
				{
					coords = nodes[q].Coords();
					Storage::real dot = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
					nodepos[q] = dot < 1.0e-10;
					if( dot*dot0 <= 0.0 ) state = 2;
				}
				if( state == 1 )
				{
					glColor3f(0.6,0.6,0.6);
					glBegin(GL_POLYGON);
					for(size_t q = 0; q < nodes.size(); q++) glVertex3dv(&nodes[q].Coords()[0]);
					glEnd();
					glColor3f(0.2,0.2,0.2);
					glBegin(GL_LINE_LOOP);
					for(size_t q = 0; q < nodes.size(); q++) glVertex3dv(&nodes[q].Coords()[0]);
					glEnd();
				}
				else if( state == 2 )
				{
					for(size_t q = 0; q < nodes.size(); ++q)
					{
						if( nodepos[q] ) //current node is inside
						{
							coords = nodes[q].Coords();
							faceverts.insert(faceverts.end(),coords.begin(),coords.end());
						}
						if( nodepos[q] != nodepos[(q+1)%nodes.size()] ) //edge is crossed by the plane, calculate crossing point
						{
							Storage::real_array sp0 = nodes[q].Coords(), sp1 = nodes[(q+1)%nodes.size()].Coords();
							Storage::real node[3];
							if( clip_plane_edge(&sp0[0],&sp1[0],p,n,node) > CLIP_NONE)
								faceverts.insert(faceverts.end(),node,node+3);
						}
					}
					glColor3f(0.6,0.6,0.6);
					glBegin(GL_POLYGON);
					for(size_t q = 0; q < faceverts.size(); q+=3) glVertex3dv(&faceverts[q]);
					glEnd();
					glColor3f(0.2,0.2,0.2);
					glBegin(GL_LINE_LOOP);
					for(size_t q = 0; q < faceverts.size(); q+=3) glVertex3dv(&faceverts[q]);
					glEnd();
					faceverts.clear();
				}
			}
		}
		
		
		if( oclipper ) 
		{
			glDisable(GL_BLEND);
			
			if( clipupdate ) 
			{
				oclipper->clip_plane(p,n);
				clip_slice.clear();
				oclipper->gen_clip(clip_slice);
				if( bclipper ) bclipper->clip_plane(clip_boundary,p,n);
				clipupdate = false;
				clipboxupdate = true;
			}
		
			if( !bclipper && !interactive && clipboxupdate )
			{
				double tt = Timer();
				clip_boundary.clear();
				std::vector<bool> nodepos;
				std::vector<Storage::real> faceverts;
				size_t pace = std::max<size_t>(1,std::min<size_t>(15,boundary_faces.size()/100));
				for(size_t k = 0; k < boundary_faces.size(); k++)
				{
					int state = 0; //0 - outside, 1 - inside, 2 - intersect
					adjacent<Node> nodes = boundary_faces[k]->getNodes();
					Storage::real_array coords = nodes[0].Coords();
					Storage::real dot0 = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
					if( dot0 < 0.0 ) state = 1; else state = 0;
					nodepos.resize(nodes.size());
					nodepos[0] = dot0 < 1.0e-10;
					for(size_t q = 1; q < nodes.size(); ++q)
					{
						coords = nodes[q].Coords();
						Storage::real dot = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
						nodepos[q] = dot < 1.0e-10;
						if( dot*dot0 <= 0.0 ) state = 2;
					}
					if( state == 1 )
					{
						face2gl f;
						f.set_color(0.6,0.6,0.6,1);
						for(size_t q = 0; q < nodes.size(); q++) f.add_vert(&nodes[q].Coords()[0]);
						f.compute_center();
						if( k%pace == 0 ) f.set_flag(true);
						f.set_elem(boundary_faces[k]->GetElementType(),boundary_faces[k]->LocalID());
						clip_boundary.push_back(f);
					}
					else if( state == 2 )
					{
						face2gl f;
						f.set_color(0.6,0.6,0.6,1);
						for(size_t q = 0; q < nodes.size(); ++q)
						{
							if( nodepos[q] ) //current node is inside
							{
								coords = nodes[q].Coords();
								f.add_vert(&coords[0]);
							}
							if( nodepos[q] != nodepos[(q+1)%nodes.size()] ) //edge is crossed by the plane, calculate crossing point
							{
								Storage::real_array sp0 = nodes[q].Coords(), sp1 = nodes[(q+1)%nodes.size()].Coords();
								Storage::real node[3];
								if( clip_plane_edge(&sp0[0],&sp1[0],p,n,node) > CLIP_NONE)
									f.add_vert(node);
							}
						}
						f.compute_center();
						if( k%pace == 0 ) f.set_flag(true);
						f.set_elem(boundary_faces[k]->GetElementType(),boundary_faces[k]->LocalID());
						clip_boundary.push_back(f);
						faceverts.clear();
					}
				}
				clipboxupdate = false;
				printf("time %g\n",Timer()-tt);
			}
			
			double tt, to = 0, tb = 0;
			size_t pace = interactive && !planecontrol ? std::max<size_t>(1,std::min<size_t>(15,oclipper->size()/100)) : 1;
			/*
			if( interactive )
			{
				glColor3f(0.6,0.6,0.6);
				oclipper->draw_clip(pace);
				if( bclipper ) bclipper->draw_clip(pace);
				glColor3f(0.2,0.2,0.2);
				oclipper->draw_clip_edges(pace);
				if( bclipper ) bclipper->draw_clip_edges(pace);
			}
			else*/
			{
				if( interactive && !planecontrol ) draw_faces_interactive(clip_slice);
				else draw_faces(clip_slice);
				glColor4f(0,0,0,1); 
				if( interactive && !planecontrol ) draw_edges_interactive(clip_slice);
				else draw_edges(clip_slice);
			}

			if( bclipper )
			{
				double tt = Timer();
				draw_faces(clip_boundary);
				glColor4f(0,0,0,1); 
				draw_edges(clip_boundary);
				printf("draw %g\n",Timer()-tt);
			}
	
			if( !bclipper && !interactive )
			{
				draw_faces(clip_boundary);
				glColor4f(0,0,0,1); 
				draw_edges(clip_boundary);
			}
			
			glEnable(GL_BLEND);
		}
		if( boundary )
		{
			if( !interactive )
			{
				for(size_t q = 0; q < all_boundary.size() ; q++) 
					all_boundary[q].compute_dist(campos);
				std::sort(all_boundary.rbegin(),all_boundary.rend());
			}
			if( interactive ) draw_faces_interactive(all_boundary);
			else draw_faces(all_boundary);
			glColor4f(0,0,0,0.25); 

			if( interactive ) draw_edges_interactive(all_boundary);
			else draw_edges(all_boundary);
			
		}
	}

	
	
	
	glutSwapBuffers();
}



int main(int argc, char ** argv)
{
	double tt;
	table[CENTROID]    = CELL | FACE;
	//table[NORMAL]      = FACE;
	//table[ORIENTATION] = FACE;
	//table[MEASURE]     = CELL|FACE;
	mesh = new Mesh();
	printf("Started loading mesh.\n");
	tt = Timer();
	mesh->SetFileOption("VERBOSITY","2");
	if( argc < 2 )
	{
		printf("Usage: %s mesh_file\n",argv[0]);
		return 0;
	}
	try
	{
		mesh->Load(argv[1]);
	} 
	catch(...)
	{
		printf("Error during loading %s\n",argv[1]);
		return -1;
	}
	printf("Done. Time %g\n",Timer()-tt);
	int fixed = 0;

	for(ElementType mask = NODE; mask <= CELL; mask = mask << 1)	
		std::cout << ElementTypeName(mask) << " " << mesh->NumberOf(mask) << std::endl;
	
	//for(Mesh::iteratorFace f = mesh->BeginFace(); f != mesh->EndFace(); f++)
	//	if( Geometry::FixNormaleOrientation(&*f) ) fixed++;
	//printf("fixed: %d\n",fixed);
	printf("Computing geometric quantities\n");
	tt = Timer();
	mesh->PrepareGeometricData(table);
	printf("Done. %g\n",Timer()-tt);
	
	printf("Calculating bounds.\n");
	tt = Timer();
	for(Mesh::iteratorNode n = mesh->BeginNode(); n != mesh->EndNode(); n++)
	{
		Storage::real_array c = n->Coords();
		if( c[0] > sright ) sright = c[0];
		if( c[0] < sleft ) sleft = c[0];
		if( c[1] > stop ) stop = c[1];
		if( c[1] < sbottom ) sbottom = c[1];
		if( c[2] > sfar ) sfar = c[2];
		if( c[2] < snear ) snear = c[2];
	}
	printf("Done. Time %g\n",Timer()-tt);
	printf("%g:%g %g:%g %g:%g\n",sleft,sright,sbottom,stop,snear,sfar);

	printf("Gathering boundary faces.\n");
	tt = Timer();
	for(Mesh::iteratorCell f = mesh->BeginCell(); f != mesh->EndCell(); f++) if( f->GetElementDimension() == 2 )
		boundary_faces.push_back(&*f);
	for(Mesh::iteratorFace f = mesh->BeginFace(); f != mesh->EndFace(); f++) if( f->GetElementDimension() == 2 )
	{
		if( f->Boundary() )
			boundary_faces.push_back(&*f);
	}
	printf("Done. Time %g\n",Timer()-tt);

	printf("Prepearing set of boundary faces for drawing.\n");
	tt = Timer();
	size_t pace = std::max<size_t>(1,std::min<size_t>(15,boundary_faces.size()/100));
	for(size_t k = 0; k < boundary_faces.size(); k++) 
	{
		all_boundary.push_back(DrawFace(boundary_faces[k]));
		if( k%pace == 0 ) all_boundary.back().set_flag(true);
	}
	printf("Done. Time %g\n",Timer() - tt);

	printf("Prepearing interactive mesh clipper.\n");
	//tt = Timer();
	oclipper = new clipper(mesh);
	bclipper = new bnd_clipper(mesh,&boundary_faces[0],boundary_faces.size());
	clipupdate = true;
	//printf("Done. Time %g\n",Timer() - tt);

	
	
	shift[0] = -(sleft+sright)*0.5;
	shift[1] = -(sbottom+stop)*0.5;
	shift[2] =  -(sfar+snear)*0.5;


	p[0] = (sleft+sright)*0.5;
	p[1] = (sbottom+stop)*0.5;
	p[2] = (sfar+snear)*0.5;

	quatinit();
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutInitWindowPosition (100, 100);
	glutCreateWindow("Graph");
	
	glDepthFunc(GL_LEQUAL);
	glClearDepth(1.f);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	
	//glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);
	//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	
	glClearColor (1.0f, 1.0f, 1.0f, 1.f);
	glutDisplayFunc(draw);
	glutReshapeFunc(reshape);
	
	glutKeyboardFunc(keyboard);
	glutMouseFunc(myclick);
	glutMotionFunc(myclickmotion);
	glutPassiveMotionFunc(mymotion);
	//glutIdleFunc(idle);
	
	glutPostRedisplay();
	glutMainLoop();
}
