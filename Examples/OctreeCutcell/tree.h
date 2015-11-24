/*********************************
 Implementation of kd-tree

 Performance
 Build: O(NlogNlogN)

 Functions
 kdtree_mem_alloc   - alloc memory for tree

 kdtree_mem_destroy - release memory

 kdtree_build       - build a kd-tree by set
                      of points

 Dependency: point.h, sort.h
   Standard: stdlib.h, stdio.h, math.h
 **********************************/
#include "point.h"


#ifndef _TREE_H
#define _TREE_H
#define KDMIN -10e50
#define KDMAX 10e50
#define LEAF  1
#define WIDTH 2
struct tree
{
	struct point * set;
	SCALE_TYPE size;
	COORD_TYPE side[DIM];
	COORD_TYPE center[DIM];
    //~ COORD_TYPE cm[DIM];
	struct tree * children;
};
void kdtree_mem_alloc(struct tree * t, SCALE_TYPE setsize);
void kdtree_mem_destroy(struct tree *t);
void kdtree_build(struct tree * t, struct point * set, SCALE_TYPE setsize, int dim);
#endif
