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
bool boundary = true, planecontrol = false, clipupdate = false, bndupdate = true, clipboxupdate = false;

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

void printtext(const char * fmt, ... )
{
	
	unsigned int i;
	char stext[1024];
	va_list ap;
	if ( fmt == NULL ) return;
	va_start(ap,fmt);
	vsprintf(stext,fmt,ap);
	va_end(ap);
	for(i=0;i<strlen(stext);i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, 
							stext[i]);
	}
}


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
	void draw_colour() const
	{
		glColor4dv(c); 
		for(unsigned k = 0; k < verts.size(); k+=3) 
		{
			glVertex3dv(cnt);
			glVertex3dv(&verts[k]);
			glVertex3dv(&verts[(k+3)%verts.size()]);
		}
	}
	void draw() const
	{
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
	double * get_vert(int k) {return &verts[k*3];}
	int size() {return verts.size()/3;}
	void set_center(double _cnt[3])
	{
		cnt[0] = _cnt[0];
		cnt[1] = _cnt[1];
		cnt[2] = _cnt[2];
	}
	void get_center(float _cnt[3])
	{
		_cnt[0] = cnt[0];
		_cnt[1] = cnt[1];
		_cnt[2] = cnt[2];
	}
	void get_center(double _cnt[3])
	{
		_cnt[0] = cnt[0];
		_cnt[1] = cnt[1];
		_cnt[2] = cnt[2];
	}
	double * get_center() {return cnt;}
	void compute_center()
	{
		cnt[0] = cnt[1] = cnt[2] = 0;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < verts.size(); k+=3)
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
	Element * get_elem(Mesh * m) {return m->ElementByLocalID(etype,id);}
};

std::vector<face2gl> all_boundary;
std::vector<face2gl> clip_boundary;

void draw_faces(std::vector<face2gl> & set, int highlight = -1)
{
	glBegin(GL_TRIANGLES);
	for(INMOST_DATA_ENUM_TYPE q = 0; q < set.size() ; q++) set[q].draw_colour();
	glEnd();
	if( highlight != -1 )
	{
		glColor4f(1,0,0,1);
		glBegin(GL_TRIANGLES);
		set[highlight].draw();
		glEnd();
	}
}

void draw_edges(std::vector<face2gl> & set, int highlight = -1)
{
	glBegin(GL_LINES);
	for(INMOST_DATA_ENUM_TYPE q = 0; q < set.size() ; q++) set[q].drawedges();
	glEnd();
	if( highlight != -1 )
	{
		glColor4f(0,1,0,1);
		glBegin(GL_LINES);
		set[highlight].drawedges();
		glEnd();
	}
	
}

void draw_faces_interactive(std::vector<face2gl> & set)
{
	glBegin(GL_TRIANGLES);
	for(INMOST_DATA_ENUM_TYPE q = 0; q < set.size() ; q++) if( set[q].get_flag() ) set[q].draw_colour();
	glEnd();
}

void draw_edges_interactive(std::vector<face2gl> & set)
{
	glBegin(GL_LINES);
	for(INMOST_DATA_ENUM_TYPE q = 0; q < set.size() ; q++) if( set[q].get_flag() ) set[q].drawedges();
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
	float bbox[6];
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
			children = static_cast<kdtree *>(malloc(sizeof(kdtree)*2));//new kdtree[2];
			children[0].marked = 0;
			children[0].children = NULL;
			children[0].set = set;
			children[0].size = size/2;
			children[1].marked = 0;
			children[1].children = NULL;
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
				for(INMOST_DATA_ENUM_TYPE k = 0; k < nodes.size(); ++k)
				{
					Storage::real_array coords = nodes[k].Coords();
					for(INMOST_DATA_ENUM_TYPE q = 0; q < 3; q++)
					{
						bbox[q*2+0] = std::min<float>(bbox[q*2+0],coords[q]);
						bbox[q*2+1] = std::max<float>(bbox[q*2+1],coords[q]);
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
	bool sub_intersect_plane_edge(Tag clip_point, Tag clip_state, std::vector<Cell *> & cells, MarkerType mrk, double p[3], double n[3])
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
				for(INMOST_DATA_ENUM_TYPE k = 0; k < ecells.size(); ++k) if( !ecells[k].GetMarker(mrk) )
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
	void sub_intersect_plane_faces(Tag clip_state, double p[3], double n[3])
	{
		if( size == 1 )
		{
			Storage::integer state;
			assert( set[0].e->GetElementDimension() == 2 );
			adjacent<Node> nodes = set[0].e->getNodes();
			Storage::real_array coords = nodes[0].Coords();
			Storage::real dot0 = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
			if( dot0 <= 0.0 ) state = CLIP_FACE_INSIDE; else state = CLIP_FACE_OUTSIDE;
			for(INMOST_DATA_ENUM_TYPE k = 1; k < nodes.size(); k++)
			{
				coords = nodes[k].Coords();
				Storage::real dot = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
				if( dot*dot0 <= 0.0 ) 
				{
					state = CLIP_FACE_INTERSECT;
					break;
				}
			}
			set[0].e->IntegerDF(clip_state) = state;
		}
		else
		{
			marked = plane_bbox(p,n);
			if( marked == 0 )
			{
				for(INMOST_DATA_ENUM_TYPE k = 0; k < size; k++) set[k].e->IntegerDF(clip_state) = CLIP_FACE_OUTSIDE;
			}
			else if( marked == 1 )
			{
				for(INMOST_DATA_ENUM_TYPE k = 0; k < size; k++) set[k].e->IntegerDF(clip_state) = CLIP_FACE_INSIDE;
			}
			else
			{
				children[0].sub_intersect_plane_faces(clip_state,p,n);
				children[1].sub_intersect_plane_faces(clip_state,p,n);
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
	void clear_children() { if( children ) {children[0].clear_children(); children[1].clear_children(); free(children);}}
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
	void intersect_plane_edge(Tag clip_point, Tag clip_state, std::vector<Cell*> & cells, MarkerType mark_cells, double p[3], double n[3])
	{
		if( marked ) 
		{
			unmark_old_edges(clip_state);
			cells.clear();
		}
		sub_intersect_plane_edge(clip_point, clip_state, cells,mark_cells,p,n);
	}
	void intersect_plane_face(Tag clip_state, double p[3], double n[3])
	{
		sub_intersect_plane_faces(clip_state, p,n);
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
	MarkerType marker;
	std::vector<Cell *> cells;
	Mesh * mm;
public:
	~clipper() 
	{
		delete tree; 
		for(INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k++) cells[k]->RemMarker(marker);
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
		double t;

		for(INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); ++k) 
			if( cells[k]->GetMarker(marker) )
			{
				cells[k]->RealArray(clips).clear();
				cells[k]->RemMarker(marker);
			}
		tree->intersect_plane_edge(clip_point,clip_state,cells,marker,p,n);
		adjacent<Face> faces;
		adjacent<Edge> edges;
		dynarray<edge_point,128> clipcoords, loopcoords;
		std::vector<bool> closed;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); ++k)
		{
			//assuming faces are convex we will have at most one clipping edge per polygon
			//otherwise every pair of clipping nodes forming edge should appear consequently
			//as long as we go through face's edges in ordered way
			clipcoords.clear();
			//we will gather all the pairs of nodes, then form closed loop
			faces = cells[k]->getFaces();
			Face * full_face = NULL;
			int ntotpoints = 0, ntotedges = 0;
			for(INMOST_DATA_ENUM_TYPE q = 0; q < faces.size(); ++q)
			{
				int last_edge_type = CLIP_NONE;
				int nfulledges = 0, npoints = 0, nstartedge = ntotedges;
				edges = faces[q].getEdges();
				for(INMOST_DATA_ENUM_TYPE r = 0; r < edges.size(); ++r)
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
				for(INMOST_DATA_ENUM_TYPE r = 0; r < nodes.size(); r++)
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
				for(INMOST_DATA_ENUM_TYPE r = 0; r < loopcoords.size(); ++r)
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
	}
	void gen_clip(std::vector<face2gl> & out)
	{
		INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,size()/100));
		for(INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k++) if( cells[k]->GetMarker(marker))
		{
			face2gl f;
			f.set_color(0.6,0.6,0.6,1);
			Storage::real_array cl = cells[k]->RealArray(clips);
			for(INMOST_DATA_ENUM_TYPE q = 0; q < cl.size(); q+=3) f.add_vert(&cl[q]);
			f.compute_center();
			f.set_elem(cells[k]->GetElementType(),cells[k]->LocalID());
			if( k%pace == 0 ) f.set_flag(true);
			out.push_back(f);
		}
	}
	void draw_clip(INMOST_DATA_ENUM_TYPE pace)
	{
		glBegin(GL_TRIANGLES);
		for(INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k+=pace) if( cells[k]->GetMarker(marker))
		{
			Storage::real_array cl = cells[k]->RealArray(clips);
			Storage::real cnt[3] = {0,0,0};
			for(INMOST_DATA_ENUM_TYPE q = 0; q < cl.size(); q+=3)
			{
				cnt[0] += cl[q+0];
				cnt[1] += cl[q+1];
				cnt[2] += cl[q+2];
			}
			cnt[0] /= (cl.size()/3);
			cnt[1] /= (cl.size()/3);
			cnt[2] /= (cl.size()/3);
			for(INMOST_DATA_ENUM_TYPE q = 0; q < cl.size(); q+=3) 
			{
				glVertex3dv(cnt);
				glVertex3dv(&cl[q]);
				glVertex3dv(&cl[(q+3)%cl.size()]);
			}
		}
		glEnd();
	}
	void draw_clip_edges(INMOST_DATA_ENUM_TYPE pace)
	{
		glBegin(GL_LINES);
		for(INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k+=pace) if( cells[k]->GetMarker(marker))
		{
			Storage::real_array cl = cells[k]->RealArray(clips);
			
			for(INMOST_DATA_ENUM_TYPE q = 0; q < cl.size(); q+=3) 
			{
				glVertex3dv(&cl[q]);
				glVertex3dv(&cl[(q+3)%cl.size()]);
			}
		}
		glEnd();
	}
	INMOST_DATA_ENUM_TYPE size() {return cells.size();}
} * oclipper = NULL;

class bnd_clipper
{
	Tag clip_state;
	kdtree * tree;
	Mesh * mm;
	Element ** faces;
	INMOST_DATA_ENUM_TYPE nfaces;
public:
	~bnd_clipper()
	{
		mm->DeleteTag(clip_state);
		delete tree;
		delete [ ]faces;
	}
	bnd_clipper(Mesh * m , Element ** _faces, INMOST_DATA_ENUM_TYPE size)
	{
		mm = m;
		clip_state = m->CreateTag("CLIP_FACE_STATE",DATA_INTEGER,FACE,NONE,1);
		nfaces = size;
		faces = new Element *[nfaces];
		for(INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k++) faces[k] = _faces[k];
		tree = new kdtree(mm,faces,nfaces);
	}
	void clip_plane(Storage::real p[3], Storage::real n[3])
	{
		tree->intersect_plane_face(clip_state,p,n);
	}
	void gen_clip(std::vector<face2gl> & out )
	{
		INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,nfaces/100));
		for(INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k++)
		{
			int state = faces[k]->IntegerDF(clip_state);
			if( state == CLIP_FACE_INSIDE )
			{
				adjacent<Node> nodes = faces[k]->getNodes();
				face2gl f;
				f.set_color(0.6,0.6,0.6,1);
				for(INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++) f.add_vert(&nodes[q].Coords()[0]);
				f.compute_center();
				f.set_elem(faces[k]->GetElementType(),faces[k]->LocalID());
				if( k%pace == 0 ) f.set_flag(true);
				out.push_back(f);
			}
			else if( state == CLIP_FACE_INTERSECT )
			{
				face2gl f;
				f.set_color(0.6,0.6,0.6,1);
				adjacent<Node> nodes = faces[k]->getNodes();
				dynarray<bool,64> nodepos(nodes.size());
				dynarray<Storage::real,64> faceverts;
				Storage::real_array coords = nodes[0].Coords();
				for(INMOST_DATA_ENUM_TYPE k = 0; k < nodes.size(); k++)
				{
					coords = nodes[k].Coords();
					Storage::real dot = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
					nodepos[k] = dot < 1.0e-10;
				}
				for(INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++)
				{
					if( nodepos[q] )
					{
						coords = nodes[q].Coords();
						f.add_vert(&coords[0]);
					}
					if( nodepos[q] != nodepos[(q+1)%nodes.size()] )
					{
						Storage::real_array sp0 = nodes[q].Coords();
						Storage::real_array sp1 = nodes[(q+1)%nodes.size()].Coords();
						Storage::real node[3];
						if( clip_plane_edge(&sp0[0],&sp1[0],p,n,node) > CLIP_NONE) f.add_vert(node);
					}
				}
				f.compute_center();
				f.set_elem(faces[k]->GetElementType(),faces[k]->LocalID());
				if( k%pace == 0 ) f.set_flag(true);
				out.push_back(f);
			}
		}
	}
	void draw_clip(INMOST_DATA_ENUM_TYPE pace)
	{
		for(INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k+=pace)
		{
			int state = CLIP_FACE_NONE;
			adjacent<Node> nodes = boundary_faces[k]->getNodes();
			Storage::real_array coords = nodes[0].Coords();
			Storage::real dot0 = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
			if( dot0 <= 0.0 ) state = CLIP_FACE_INSIDE; else state = CLIP_FACE_OUTSIDE;
			for(INMOST_DATA_ENUM_TYPE q = 1; q < nodes.size(); ++q)
			{
				coords = nodes[q].Coords();
				Storage::real dot = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
				if( dot*dot0 <= 0.0 ) 
				{
					state = CLIP_FACE_INTERSECT;
					break;
				}
			}
			if( state == CLIP_FACE_INSIDE )
			{
				adjacent<Node> nodes = faces[k]->getNodes();
				glBegin(GL_POLYGON);
				for(INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++) glVertex3dv(&nodes[q].Coords()[0]);
				glEnd();
			}
			faces[k]->IntegerDF(clip_state) = state;
		}
	}
	void draw_clip_edges(INMOST_DATA_ENUM_TYPE pace)
	{
		for(INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k+=pace)
		{
			if( faces[k]->IntegerDF(clip_state) == CLIP_FACE_INSIDE )
			{
				adjacent<Node> nodes = faces[k]->getNodes();
				glBegin(GL_LINE_LOOP);
				for(INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++) glVertex3dv(&nodes[q].Coords()[0]);
				glEnd();
			}
		}
	}
	INMOST_DATA_ENUM_TYPE size() {return nfaces;}
} * bclipper = NULL;

class kdtree_picker
{
	struct entry
	{
		int index;
		float xyz[3];
	} * set;
	INMOST_DATA_ENUM_TYPE size;
	Storage::real bbox[6];
	kdtree_picker * children;
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
	void kdtree_build(int dim, std::vector<face2gl> & in, struct entry * temp)
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
			children = static_cast<kdtree_picker *>(malloc(sizeof(kdtree_picker)*2));//new kdtree_picker[2];
			children[0].children = NULL;
			children[0].set = set;
			children[0].size = size/2;
			children[1].children = NULL;
			children[1].set = set+size/2;
			children[1].size = size - size/2;
			children[0].kdtree_build((dim+1)%3,in,temp);
			children[1].kdtree_build((dim+1)%3,in,temp);
			for(int k = 0; k < 3; k++)
			{
				bbox[0+2*k] = std::min(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
				bbox[1+2*k] = std::max(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
			}
		}
		else 
		{
			assert(size == 1);
			bbox[0] = bbox[2] = bbox[4] = 1.0e20;
			bbox[1] = bbox[3] = bbox[5] = -1.0e20;
			for(INMOST_DATA_ENUM_TYPE k = 0; k < in[set[0].index].size(); ++k)
			{
				double * coords = in[set[0].index].get_vert(k);
				for(INMOST_DATA_ENUM_TYPE q = 0; q < 3; q++)
				{
					bbox[q*2+0] = std::min(bbox[q*2+0],coords[q]);
					bbox[q*2+1] = std::max(bbox[q*2+1],coords[q]);
				}
			}
		}
	}
	void clear_children() { if( children ) {children[0].clear_children(); children[1].clear_children(); free(children);}}
	int raybox(double pos[3], double ray[3], double closest)
	{
		double tnear = -1.0e20, tfar = 1.0e20, t1,t2,c;
		for(int i = 0; i < 3; i++)
		{
			if( fabs(ray[i]) < 1.0e-15 )
			{
				if( pos[i] < bbox[i*2] || pos[i] > bbox[i*2+1] )
					return 0;
			}
			else
			{
				t1 = (bbox[i*2+0] - pos[i])/ray[i];
				t2 = (bbox[i*2+1] - pos[i])/ray[i];
				if( t1 > t2 ) {c = t1; t1 = t2; t2 = c;}
				if( t1 > tnear ) tnear = t1;
				if( t2 < tfar ) tfar = t2;
				if( tnear > closest ) return 0;
				if( tnear > tfar ) return 0;
				if( tfar < 0 ) return 0;
			}
		}
		return 1;
	}
	void sub_intersect_ray_faces(std::vector<face2gl> & in, double p[3], double dir[3], std::pair<double,int> & closest)
	{
		if( size == 1 )
		{
			face2gl & f = in[set[0].index];
			double * tri[3], btri[3][3], dot[3], prod[3][3], norm[3], d, proj[3];
			tri[0] = f.get_center();
			Storage::real maxdist = 0, dist;
			for(INMOST_DATA_ENUM_TYPE i = 0; i < f.size(); i++)
			{
				INMOST_DATA_ENUM_TYPE j = (i+1)%f.size();
				tri[1] = f.get_vert(i);
				tri[2] = f.get_vert(j);
				norm[0] = (tri[2][1]-tri[0][1])*(tri[1][2]-tri[0][2]) - (tri[2][2]-tri[0][2])*(tri[1][1]-tri[0][1]);
				norm[1] = (tri[2][2]-tri[0][2])*(tri[1][0]-tri[0][0]) - (tri[2][0]-tri[0][0])*(tri[1][2]-tri[0][2]);
				norm[2] = (tri[2][0]-tri[0][0])*(tri[1][1]-tri[0][1]) - (tri[2][1]-tri[0][1])*(tri[1][0]-tri[0][0]);
				d = norm[0]*(tri[0][0]-p[0])+norm[1]*(tri[0][1]-p[1])+norm[2]*(tri[0][2]-p[2]);
				d /= norm[0]*dir[0] + norm[1]*dir[1] + norm[2]*dir[2];
				proj[0] = p[0] + d*dir[0];
				proj[1] = p[1] + d*dir[1];
				proj[2] = p[2] + d*dir[2];

				if( d > closest.first ) break;

				for(int k = 0; k < 3; k++)
				{
					btri[k][0] = tri[k][0] - proj[0];
					btri[k][1] = tri[k][1] - proj[1];
					btri[k][2] = tri[k][2] - proj[2];
				}
				for(int k = 0; k < 3; k++)
				{
					int l = (k+1)%3;
					prod[k][0] = btri[k][1]*btri[l][2] - btri[k][2]*btri[l][1];
					prod[k][1] = btri[k][2]*btri[l][0] - btri[k][0]*btri[l][2];
					prod[k][2] = btri[k][0]*btri[l][1] - btri[k][1]*btri[l][0];
				}

				for(int k = 0; k < 3; k++)
				{
					int l = (k+1)%3;
					dot[k] = prod[k][0]*prod[l][0]+prod[k][1]*prod[l][1]+prod[k][2]*prod[l][2];
				}
				if( dot[0] >= 0 && dot[1] >= 0 && dot[2] >= 0 ) 
				{
					closest.first = d;
					closest.second = set[0].index;
					break; //don't expect anything better here
				}
			}
		}
		else
		{
			if( raybox(p,dir,closest.first) )
			{
				children[0].sub_intersect_ray_faces(in,p,dir,closest);
				children[1].sub_intersect_ray_faces(in,p,dir,closest);
			}
		}
	}
public:
	kdtree_picker() : set(NULL), size(0), children(NULL) {}
	~kdtree_picker() { clear_children(); delete [] set;}
	kdtree_picker(std::vector<face2gl> & in)
	{
		size = in.size();
		set = new struct entry[size];
		for(INMOST_DATA_ENUM_TYPE k = 0; k < in.size(); ++k)
		{
			set[k].index = k;
			in[k].get_center(set[k].xyz);
		}
		struct entry * temp = new struct entry[size];
		kdtree_build(0,in,temp);
		delete [] temp;
	}
	int ray_faces(std::vector<face2gl> & in, double p[3], double dir[3])
	{
		std::pair<double, int> closest(1.0e20,-1);
		sub_intersect_ray_faces(in,p,dir,closest);
		return closest.second;
	}
};

class picker
{
	kdtree_picker * tree;
	std::vector<face2gl> * faces;
public:
	~picker() {delete tree;}
	picker(std::vector<face2gl> & _faces)
	{
		faces = &_faces;
		tree = new kdtree_picker(*faces);
	}
	int select(double p[3], double ray[3]) {return tree->ray_faces(*faces,p,ray);};
} * current_picker = NULL;


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
	{
		const double pi = 3.1415926535897932384626433832795;
		const double znear = 0.01;
		const double zfar  = 10000.0;
		const double fH = znear * tan( 60.0 / 360.0 * pi);
		double fW = fH * aspect;
		glFrustum(-fW,fW,-fH,fH,znear,zfar);
	}
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();
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
			bndupdate = true;
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
			bndupdate = true;
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
		else bndupdate = true;
	}
}
void mymotion(int nmx, int nmy) // Mouse
{
	motion(nmx,nmy);
	mymx = 2.*(nmx/(double)width - 0.5);
	mymy = 2.*(0.5 - nmy/(double)height);
	if( current_picker != NULL ) glutPostRedisplay();
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
		if( current_picker ) delete current_picker;
		delete mesh;
		exit(-1);
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
		bndupdate = true;
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
			bndupdate = true;
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
			bndupdate = true;
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
			bndupdate = true;
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
			bndupdate = true;
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
			bclipper->clip_plane(p,n);
			clipboxupdate = true;
		}
			
		glutPostRedisplay();
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
	adjacent<Node> nodes = f->getNodes();

	if( f->nbAdjElements(CELL) == 0 ) ret.set_color(1,0,0,0.1);
	else ret.set_color(0,1,0,0.1);
	
	for( adjacent<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
		ret.add_vert(&(kt->Coords()[0]));
	ret.set_elem(f->GetElementType(),f->LocalID());
	ret.compute_center();
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
   double l = sqrt(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
   if( l ) 
   {
	   direction[0] /= l;
	   direction[1] /= l;
	   direction[2] /= l;
   }
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

	

	double campos[3] = {0.5,0.5,0}, pickp[3], pickd[3];
	whereami(campos[0],campos[1],campos[2]);
	int picked = -1;

	//glTranslated((l+r)*0.5,(b+t)*0.5,(near+far)*0.5);

	{
		
		if( oclipper ) 
		{
			
			if( clipupdate ) 
			{
				if( current_picker != NULL ) {delete current_picker; current_picker = NULL;}
				oclipper->clip_plane(p,n);
				clipupdate = false;
				clipboxupdate = true;
			}
		
			if( !interactive && clipboxupdate )
			{
				clip_boundary.clear();
				oclipper->gen_clip(clip_boundary);
				bclipper->clip_plane(p,n);
				bclipper->gen_clip(clip_boundary);
				clipboxupdate = false;

				if( current_picker != NULL ) {delete current_picker; current_picker = NULL;}
				current_picker = new picker(clip_boundary);
			}

			

			if( current_picker != NULL )
			{
				pick_mouse(pickp,pickd);
				picked = current_picker->select(pickp,pickd);
			}
	
			if( interactive && clipboxupdate )
			{
				//printf("draw1 %d %d\n",interactive, clipboxupdate);
				INMOST_DATA_ENUM_TYPE opace = !planecontrol ? std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,oclipper->size()/100)) : 1;
				INMOST_DATA_ENUM_TYPE bpace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,bclipper->size()/100));
				glColor4f(0.6,0.6,0.6,1);
				oclipper->draw_clip(opace);
				bclipper->draw_clip(bpace);
				glColor4f(0,0,0,1); 
				oclipper->draw_clip_edges(opace);
				bclipper->draw_clip_edges(bpace);
			}
			else
			{
				//printf("draw2 %d %d\n",interactive, clipboxupdate);
				if( interactive )
				{
					draw_faces_interactive(clip_boundary);
					glColor4f(0,0,0,1); 
					draw_edges_interactive(clip_boundary);
				}
				else
				{
					draw_faces(clip_boundary,picked);
					glColor4f(0,0,0,1); 
					draw_edges(clip_boundary,picked);
				}
			}
		}
		if( boundary )
		{
			glEnable(GL_BLEND);
			if( !interactive && bndupdate)
			{
				for(INMOST_DATA_ENUM_TYPE q = 0; q < all_boundary.size() ; q++) 
					all_boundary[q].compute_dist(campos);
				std::sort(all_boundary.rbegin(),all_boundary.rend());
				bndupdate = false;
			}
			glColor4f(0,0,0,0.25); 
			if( interactive ) draw_edges_interactive(all_boundary);
			else draw_edges(all_boundary);
			

			if( interactive ) draw_faces_interactive(all_boundary);
			else draw_faces(all_boundary);
			
			glDisable(GL_BLEND);
		}
	}

	if( picked != -1 )
	{
		glDisable(GL_DEPTH_TEST);
		glLoadIdentity();
		set_matrix2d();
		

		double top = 0.96, left = 0.11, interval = 0.04, bottom = top;
		
		Element * e = clip_boundary[picked].get_elem(mesh);
		for(Mesh::iteratorTag t = mesh->BeginTag(); t != mesh->EndTag(); ++t) if( t->isDefined(e->GetElementType()) )
			if( e->HaveData(*t) )
				bottom -= interval;

		glColor3f(1,1,1);
		glBegin(GL_QUADS);
		glVertex2f(left-0.01,bottom-0.01);
		glVertex2f(left-0.01,0.99);
		glVertex2f(0.99,0.99);
		glVertex2f(0.99,bottom-0.01);
		glEnd();
		glColor3f(0,0,0);
		glBegin(GL_LINE_LOOP);
		glVertex2f(left-0.01,bottom-0.01);
		glVertex2f(left-0.01,0.99);
		glVertex2f(0.99,0.99);
		glVertex2f(0.99,bottom-0.01);
		glEnd();

		
		glColor3f(0.2,0.2,0.2);
		glRasterPos2f(left,top);
		printtext("%s %d",ElementTypeName(e->GetElementType()),e->LocalID());
		top -= interval;
		glColor3f(0.2,0.2,0.2);
		for(Mesh::iteratorTag t = mesh->BeginTag(); t != mesh->EndTag(); ++t) if( t->isDefined(e->GetElementType()) )
		{
			if( e->HaveData(*t) )
			{
				char str[1024];
				sprintf(str,"%s %s",t->GetTagName().c_str(),DataTypeName(t->GetDataType()));
				switch(t->GetDataType())
				{
				case DATA_INTEGER:
					{
						Storage::integer_array arr = e->IntegerArray(*t);
						for(INMOST_DATA_ENUM_TYPE k = 0; k < arr.size(); k++) sprintf(str,"%s %d",str,arr[k]);
						break;
					}
				case DATA_REAL:
					{
						Storage::real_array arr = e->RealArray(*t);
						for(INMOST_DATA_ENUM_TYPE k = 0; k < arr.size(); k++) sprintf(str,"%s %lf",str,arr[k]);
						break;
					}
				case DATA_BULK:
					{
						Storage::bulk_array arr = e->BulkArray(*t);
						for(INMOST_DATA_ENUM_TYPE k = 0; k < arr.size(); k++) sprintf(str,"%s %c",str,arr[k]);
						break;
					}
				case DATA_REFERENCE:
					{
						Storage::reference_array arr = e->ReferenceArray(*t);
						for(INMOST_DATA_ENUM_TYPE k = 0; k < arr.size(); k++) sprintf(str,"%s %s:%d",str,ElementTypeName(arr[k]->GetElementType()),arr[k]->LocalID());
						break;
					}
				}
				glRasterPos2f(left,top);
				printtext(str);
				top -= interval;
			}
		}


		glEnable(GL_DEPTH_TEST);
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
		printf("Usage: %s mesh_file [dims]\n",argv[0]);
		return 0;
	}
	//try
	{
		if( argc > 2 )	mesh->SetFileOption("VTK_GRID_DIMS",argv[2]);
		mesh->Load(argv[1]);
	} 
	/*
	catch(...)
	{
		printf("Error during loading %s\n",argv[1]);
		return -1;
	}
	*/
	printf("Done. Time %g\n",Timer()-tt);
	int fixed = 0;
	

	for(ElementType mask = NODE; mask <= CELL; mask = mask << 1)	
		std::cout << ElementTypeName(mask) << " " << mesh->NumberOf(mask) << std::endl;
	
	//return 0;
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
	for(Mesh::iteratorCell f = mesh->BeginCell(); f != mesh->EndCell(); f++) 
	{
		if( f->GetElementDimension() == 2 )
			boundary_faces.push_back(&*f);
	}
	for(Mesh::iteratorFace f = mesh->BeginFace(); f != mesh->EndFace(); f++) if( f->GetElementDimension() == 2 )
	{
		if( f->Boundary() )
			boundary_faces.push_back(&*f);
	}
	printf("Done. Time %g\n",Timer()-tt);

	printf("Prepearing set of boundary faces for drawing.\n");
	tt = Timer();
	INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,boundary_faces.size()/100));
	for(INMOST_DATA_ENUM_TYPE k = 0; k < boundary_faces.size(); k++) 
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
