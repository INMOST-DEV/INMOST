#ifndef _OCTGRID_H
#define _OCTGRID_H
#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DIM 3 
//must be 2 or 3

using namespace INMOST;

#define TSNH {printf("(%s:%d) this should not happen!\n",__FILE__,__LINE__); exit(-1);}
#define WAITNL {int ret; char c; printf("press Enter\n"); ret = scanf("%c",&c);}

#define VBUFFER 131072
#define CBUFFER 131072
#define VBUFFER 131072
#define CBUFFER 131072

#define DEF_ARR(x,type) type * x; int n##x;
#define INI_ARR(x,type,BUF) {x = (type *)malloc(sizeof(type)*BUF); n##x = 0; if( x == NULL ) {printf("(%s:%d) out of mem\n",__FILE__,__LINE__); exit(-1);}}
#define ADD_ARR(x,val) {x[n##x++] = val;}
#define POP_ARR(x) x[--(n##x)]
#define INC_ARR(x,type,BUF) if( (n##x)%BUF == 0 ) { x = (type *)realloc(x,sizeof(type)*((n##x)/BUF+1)*BUF); if( x == NULL ) {printf("(%s:%d) out of mem\n",__FILE__,__LINE__); exit(-1);}}
#define KIL_ARR(x) {if( x != NULL ) free(x); x = NULL; n##x = 0;}
#define CNT_ARR(x) (n##x)
#define PRM_ARR(x,type) type * x, int n##x
#define PUT_ARR(x) (x),(n##x)
#define LAST_ARR(x) (x[n##x-1])

#define C_INI_ARR(c,x,type,BUF) {c->x = (type *)malloc(sizeof(type)*BUF); c->n##x = 0; if( c->x == NULL ) {printf("(%s:%d) out of mem\n",__FILE__,__LINE__); exit(-1);}}
#define C_ADD_ARR(c,x,val) {c->x[c->n##x++] = val;}
#define C_POP_ARR(c,x) c->x[--(c->n##x)]
#define C_INC_ARR(c,x,type,BUF) if( (c->n##x)%BUF == 0 ) { c->x = (type *)realloc(c->x,sizeof(type)*((c->n##x)/BUF+1)*BUF); if( c->x == NULL ) {printf("(%s:%d) out of mem\n",__FILE__,__LINE__); exit(-1);}}
#define C_KIL_ARR(c,x) {if( c->x != NULL ) free(c->x); c->x = NULL; c->n##x = 0;}
#define C_CNT_ARR(c,x) (c->n##x)
#define C_LAST_ARR(c,x) (c->x[(c->n##x)-1])
#define C_PUT_ARR(c,x) (c->x),(c->n##x)

typedef std::vector<HandleType> cell_vector;
typedef std::map<Tag, Storage::real > data_by_mat;
typedef std::pair< Storage::real, std::map<Tag, Storage::real> > vol_and_data_by_mat;
typedef std::map<Storage::integer, vol_and_data_by_mat > data_storage;
typedef std::vector<Storage::integer> mat_ret_type;

struct cell
{	
	int busy;
	int leaf;
	int level;
	int vertexes[1<<DIM];
	int children[1<<DIM];
	int parent;
	
	//data
	double center[3], side[3];
	double vol;
	data_storage * data;
	cell_vector * mr;
};


struct vert
{
	int busy;
	int env[1<<DIM];
	HandleType mv;
};

struct grid
{
	Mesh * mesh;
	MarkerType octree_node;
	Tag materials, cell_material, parent, united, problem,edge_face_number, face_center_node, Kvec, K;
	
	
	DEF_ARR(vempty,int);
	DEF_ARR(verts,struct vert);
	
	DEF_ARR(cempty,int);
	DEF_ARR(cells,struct cell);

	int n[3];
	int *** mgrid;
	
	void (*transformation)(double xyz[3]);
	int (*cell_should_unite)(struct grid * g, int cell);
	int (*cell_should_split)(struct grid * g, int cell);
	void (*cell_unite_data)(struct grid * g, int cell);
	void (*cell_split_data)(struct grid * g, int cell);
	void (*cell_init_data)(struct grid * g, int cell);
	void (*cell_destroy_data)(struct grid * g, int cell);
	void (*cell_to_INMOST)(struct grid * g, int cell,Cell r);
	void (*vert_interpolate_data)(struct grid * g, int big_cell, int nvert, int * verts, int * isnew);
	void (*vert_init_data)(struct grid * g, int vert);
	void (*vert_destroy_data)(struct grid * g, int vert);
	void (*vert_to_INMOST)(struct grid * g, int vert,Node v);
	void (*init_mesh)(struct grid *g);
	std::map<Tag,Storage::real> (*cell_small_unite)(ElementArray<Cell> & unite);
	mat_ret_type (*get_material_types)(double xyz[3]);
};

void vertGetCoord(struct grid * g, int v, double coord[3]);
void vertSplitData(struct grid * g, int big_cell, int nverts, int * verts, int * isnew);

int cellSearch(struct grid * g, int m, double v[3]);
int cellAround(struct grid * g, int m, int side, int neighbours[1<<(DIM-1)]);
int cellSplit(struct grid * g, int m);
int cellUnite(struct grid * g, int m);
int cellRefine(struct grid * g, int m);
int cellCoarse(struct grid * g, int m);
void cellUniteData(struct grid * g, int m); // rise data to coarser cells
void cellSplitData(struct grid * g, int m); // put data to finer cells

int cellAllAround(struct grid * g, int m, int around[(1<<DIM)*((1<<DIM)-1)]);

//point should be in coordinate system of transformation
//interpolation is provided on transformed coordinate system
// z-coordinate for DIM = 2 is 0.5
int cellInterpolationPattern(struct grid * g, int m, double point[3], 
							 int cells[(1<<DIM)*((1<<DIM)-1)+1], 
							 double coefs[(1<<DIM)*((1<<DIM)-1)+1]);
//point should be in coordinate system of transformation
//interpolation is provided on transformed coordinate system
// z-coordinate for DIM = 2 is 0.5
int cellGradientPattern(struct grid * g, int m, double point[3], 
						int cells[(1<<DIM)*((1<<DIM)-1)+1], 
						double coefs[(1<<DIM)*((1<<DIM)-1)+1][DIM]);

void gridInit(struct grid * g, int n[3]);
void gridRefine(struct grid * g);
void gridCoarse(struct grid * g);
void gridUniteData(struct grid * g);
void gridRecreateINMOST(struct grid * g);
void gridAMR(struct grid * g, int recreate);
void gridDelete(struct grid * g);
int gridSearch(struct grid * g, double v[3]);
void gridUniteSmallElements(struct grid * g);



void default_transformation(double xyz[3]);
int default_cell_should_unite(struct grid * g, int cell);
int default_cell_should_split(struct grid * g, int cell);
void default_cell_unite_data(struct grid * g, int cell);
void default_cell_split_data(struct grid * g, int cell);
void default_vert_interpolate_data(struct grid * g, int big_cell, int nvert, int * verts, int * isnew);
void default_vert_init_data(struct grid * g, int vert);
void default_cell_init_data(struct grid * g, int cell);
void default_vert_destroy_data(struct grid * g, int vert);
void default_cell_destroy_data(struct grid * g, int cell);
void default_vert_to_INMOST(struct grid * g, int vert, Node v);
void default_cell_to_INMOST(struct grid * g, int cell, Cell r);

#endif
