#ifndef _OCTGRID_H
#define _OCTGRID_H
#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DIM 3 

using namespace INMOST;

#define TSNH {printf("(%s:%d) this should not happen!\n",__FILE__,__LINE__); exit(-1);}
#define WAITNL {int ret; char c; printf("press Enter\n"); ret = scanf("%c",&c);}

#define VBUFFER 131072
#define CBUFFER 131072

/// Struct conatains tags for one cell
struct Cell_Tags
{
    /// Center of cell (3 coordinates)
    Tag center;  
    /// Sides of cell (3 doubles)
    Tag side;    
    /// Number of cell in z dimension (up/down). 
    Tag floor;    
    /// Only for child cell. Octree division means 8 children. This tag contains id among other children.
    /// Numeration: 1 2 above them 5 6
    ///             0 3            4 7 
    Tag chld_num; 
    /// Child id for parent
    Tag par_chld_nums; 
    /// Level of partition
    Tag level;
    /// Cell should be split
    Tag to_split; 

    Tag busy;
    Tag leaf;
    Tag vertexes;
    Tag children;
    Tag parent;
    Tag is_valid;
    /// Number of processor before redistribute
    Tag proc; 
    Tag i;        
    Tag base_id; 
};

/// Main object, contains all information about mesh
struct grid
{
    int last_base_id;

    /// Inmost mesh
	Mesh* mesh;
    Cell_Tags c_tags;
	
	void (*transformation)(double xyz[3]);
	void (*rev_transformation)(double xyz[3]);
	int (*cell_should_split)(struct grid * g, Cell cell);
	int (*cell_should_unite)(struct grid * g, Cell cell);
};

/// Often after redistribution brother cells splits to different processor
/// This function corrects this division. Puts all children to one processor
void old_correct_brothers(struct grid* g, int size, int rank, int type);
void correct_brothers(struct grid* g, int size, int rank, int type);
void pre_eval(struct grid* g, int size, int rank);
void init_mesh(struct grid* g);
void fill_proc_tag(grid* g);
int calc_sends(grid* g);

/// Debug. Print information about all cells
void print_all_about_cells(struct grid * g);
void print_cell_center(struct grid* g, Cell cell);
void print_node_center(struct grid* g, Node node);
void print_edge(struct grid* g, Edge edge);
void print_face_nodes(struct grid* g, Face face);
void print_face_edges(struct grid* g, Face face);
void print_cells_statistics(grid* g);
void command(grid* g);
void dump_to_vtk(grid* g, const char* suffix="");
void redistribute(grid* g,int type);

void print_redist_tag(struct grid* g,  int rank);
void ghost_to_shared(grid* g);
void remove_ghost_edges(grid* g);
void resolve_ghost(grid* g);

/// Initialize grid. Create tags and other data
void gridInit(struct grid * g, int n[3]);
/// Split cells using "cell_should_split" rule
void gridRefine(struct grid * g);
/// Unite cells using "cell_should_unite" rule
void gridCoarse(struct grid * g);
/// Unite and coarse cells
void gridAMR(struct grid * g, int action);

/// Default method of transformation
void default_transformation(double xyz[3]);
/// Default method of unite
int default_cell_should_unite(struct grid * g, int cell);
/// Default method of split
int default_cell_should_split(struct grid * g, int cell);
#endif
