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
   Standard: stdlib.h, stdio.h math.h
 **********************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "tree.h"
//#include "sort.h"



void kdtree_mem_alloc(struct tree * t, SCALE_TYPE setsize)
{
	SCALE_TYPE i;
	if( t == NULL )
	{
		printf("Tree pointer is NULL. Can not alloc memory for tree.\n");
		return;
	}
//	printf("malloc %d\n",setsize);
	if( setsize > LEAF )
	{
		SCALE_TYPE subsize = setsize/WIDTH;
		if( setsize > 0 && subsize == 0 )
			subsize = 1;
		t->children = (struct tree *)malloc(sizeof(struct tree)*WIDTH);
		if( t->children == NULL )
		{
			printf("Cannot allocate memory!");
			return;
		}
		for(i = 0; i < WIDTH; i++)
		{
			if( i*subsize>setsize-1 )
				subsize = 0;
			if( i == WIDTH-1 && subsize != 0 )
				subsize = setsize-i*subsize;
			kdtree_mem_alloc(t->children+i,subsize);
		}
	}
	else
	{
		t->children = NULL;
	}
}
void kdtree_mem_destroy(struct tree * t)
{
	if( t->children != NULL )
	{
		SCALE_TYPE i;
		for(i = 0; i < WIDTH; i++)
			kdtree_mem_destroy(t->children+i);
		free(t->children);
		t->children = NULL;
	}
}

int cmp0(const void * pa, const void * pb)
{
	struct point * a = (struct point *)pa;
	struct point * b = (struct point *)pb;
	return (int)ceil(a->coord[0]-b->coord[0]);
}

int cmp1(const void * pa, const void * pb)
{
	struct point * a = (struct point *)pa;
	struct point * b = (struct point *)pb;
	return (int)ceil(a->coord[1]-b->coord[1]);
}

int cmp2(const void * pa, const void * pb)
{
	struct point * a = (struct point *)pa;
	struct point * b = (struct point *)pb;
	return (int)ceil(a->coord[2]-b->coord[2]);
}


void kdtree_build(struct tree * t, struct point * set, SCALE_TYPE setsize, int dim)
{
	SCALE_TYPE i,j;
	t->set = set;
	t->size = setsize;
	if( t == NULL )
	{
		printf("Tree pointer is NULL. Can not build tree.\n");
		return;
	}
	if( setsize <= LEAF )
	{
		if( setsize > 0 )
		{
			COORD_TYPE max, min;
//			printf("leaf size %d\n",setsize);
//		printf("Leaf\n");
			for(j = 0; j < DIM; j++)
			{
				max = set[0].coord[j]+set[0].radius[j];
				min = set[0].coord[j]-set[0].radius[j];
                //~ t->cm[j] = set[0].coord[j];
				for(i = 1; i < setsize; i++)
				{
					if( set[i].coord[j]+set[i].radius[j] > max )
						max = set[i].coord[j]+set[i].radius[j];
					if( set[i].coord[j]-set[i].radius[j] < min )
						min = set[i].coord[j]-set[i].radius[j];
					//~ t->cm[j] += set[i].coord[j];
				}
				t->side[j] = (max-min)/2.0;
				t->center[j] = (max+min)/2.0;
                //~ t->cm[j] /= (COORD_TYPE)setsize;
//				printf("[%d]: %f %f\n",j,t->center[j],t->side[j]);
			}
//			printf("\n");
		}
		kdtree_mem_destroy(t);
	}
	else
	{
		struct point * subset;
		COORD_TYPE max,min;
		SCALE_TYPE subsize = setsize/WIDTH;
//		printf("size %d\n",setsize);
		if( setsize > 0 && subsize == 0 )
			subsize = 1;
//		printf("sub %d size %d width %d\n",subsize, setsize, WIDTH);
		//if( setsize > 5000 && sizeof(COORD_TYPE) == 4) radix_sort(set,setsize,dim);
		//else
		switch(dim)
		{
			case 0: qsort(set,setsize,sizeof(struct point),cmp0); break;
			case 1: qsort(set,setsize,sizeof(struct point),cmp1); break;
			case 2: qsort(set,setsize,sizeof(struct point),cmp2); break;
		}
		//quick_sort(set,setsize,dim);
		for(i = 0; i < WIDTH; i++)
		{
			subset = set + i*subsize;
			if( i*subsize>setsize-1 )
				subsize = 0;
			if( i == WIDTH-1 && subsize != 0 )
				subsize = setsize-i*subsize;
			if( t->children == NULL )
				kdtree_mem_alloc(t,setsize);
			// почему-то раньше здесь последним аргументом стояло dim+1!!!!!
			kdtree_build(t->children+i,subset,subsize,(dim+1)%DIM);
		}
//		printf("Root\n");
		for(j = 0; j < DIM; j++)
		{
        	COORD_TYPE num = 0;
//			t->center[j] = t->children[0].center[j];
			max = -10e35; //t->children[0].center[j]+t->children[0].side[j]/2.0;
			min = 10e35; //t->children[0].center[j]-t->children[0].side[j]/2.0;
//			printf("%f %f\n",t->children[0].center[j],t->children[0].side[j]);
//			printf("%f\n",t->children[0].center[j]+t->children[i].side[j]/2.0);
//			printf("%f %f\n",max,min);
			for(i = 0; i < WIDTH; i++)
				if(t->children[i].size > 0 )
				{
					//t->center[j] += t->children[i].center[j];
					if( t->children[i].center[j]+t->children[i].side[j] > max )
						max = t->children[i].center[j]+t->children[i].side[j];
					if( t->children[i].center[j]-t->children[i].side[j] < min )
						min = t->children[i].center[j]-t->children[i].side[j];
                    //~ t->cm[j] += t->children[i].cm[j];
                    num = num + 1;
				}
			t->side[j] = (max-min)/2.0;
			t->center[j] = (max+min)/2.0;
            //~ t->cm[j] /= (COORD_TYPE)num;
//			printf("[%d]: %f %f %f %f\n",j,t->center[j],t->side[j],max,min);
		}
//		printf("\n");
	}
	return;
}
