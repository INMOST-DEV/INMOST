#include <math.h>
#include "triang.h"


double trieps = 1.0e-8;

void set_eps(double new_eps)
{
	trieps = new_eps;
}


int compare(XYZ a , XYZ b)
{
	if( fabs(a.x - b.x) < trieps && fabs(a.y-b.y) < trieps && fabs(a.z-b.z) < trieps )
		return 1;
	return 0;
}



/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/
XYZ VertexInterp (double isolevel, XYZ p1, XYZ p2, double valp1, double valp2)
{
	double mu;
	XYZ p;

	if (fabs (isolevel - valp1) < 1e-13)
		return (p1);
	if (fabs (isolevel - valp2) < 1e-13)
		return (p2);
	if (fabs (valp1 - valp2) < 1e-13)
		return (p1);
	mu = (isolevel - valp1) / (valp2 - valp1);
	p.x = p1.x + mu * (p2.x - p1.x);
	p.y = p1.y + mu * (p2.y - p1.y);
	p.z = p1.z + mu * (p2.z - p1.z);

	return (p);
}



XYZ VertexInterpCh (XYZ p1, XYZ p2,int f1,int f2, int depth,bool (*inside)(double x, double y, double z))
{
	int i, fm;
	int plus, minus;
	plus = minus = 0;
	XYZ pm;
	XYZ pp1, pp2;
	pp1 = p1;
	pp2 = p2;
	
	for (i = 0; i < depth; i++)
	{
		pm.x = (p1.x + p2.x) * 0.5;
		pm.y = (p1.y + p2.y) * 0.5;
		pm.z = (p1.z + p2.z) * 0.5;
		fm = inside (pm.x, pm.y, pm.z) ? 1 : -1;
		if (f1 == fm)
		{
			minus++;
			p1.x = pm.x;
			p1.y = pm.y;
			p1.z = pm.z;
		}
		else
		{
			plus++;
			p2.x = pm.x;
			p2.y = pm.y;
			p2.z = pm.z;
		}
	}
	if (plus == i) return (p1); //!!!!!!!!
	if (minus == i) return (p2);	
//	if (dist(pm,pp1) < dist(pp1,pp2) / 1024)
	return (pm);
}

/*
   Given a grid cell and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the array "triangles"
   will be loaded up with the vertices at most 5 triangular facets.
   0 will be returned if the grid cell is either totally above
   of totally below the isolevel.
*/


const int edgeTable[256] = {
0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };
const int triTable[256][19] =
	{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1, -1, -1, -1},
	{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1, -1, -1, -1},
	{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1, -1, -1, -1},
	{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1, -1, -1, -1},
	{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1, -1, -1, -1},
	{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1, -1, -1, -1},
	{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1, -1, -1, -1},
	{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1, -1, -1, -1},
	{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1, -1, -1, -1},
	{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1, -1, -1, -1},
	{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1, -1, -1, -1},
	{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1, -1, -1, -1},
	{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1, -1, -1, -1},
	{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1, -1, -1, -1},
	{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1, -1, -1, -1},
	{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1, -1, -1, -1},
	{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1, -1, -1, -1},
	{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1, -1, -1, -1},
	{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1, -1, -1, -1},
	{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1, -1, -1, -1},
	{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1, -1, -1, -1},
	{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1, -1, -1, -1},
	{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1, -1, -1, -1},
// 	{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},	// !!!!!!!!!!!!!! 0111 1101 (6, 18)
	{0, 9, 7, 11, 6, 1, 9, 1, 6, 6, 7, 9, 1, 0, 11, 7, 11, 0, -1},			// !!!!!!!!!!!!!! 0111 1101 (6, 26)
	{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1, -1, -1, -1},
	{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1, -1, -1, -1},
	{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1, -1, -1, -1},
	{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1, -1, -1, -1},
	{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1, -1, -1, -1},
	{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1, -1, -1, -1},
	{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1, -1, -1, -1},
	{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1, -1, -1, -1},
	{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1, -1, -1, -1},
	{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1, -1, -1, -1},
	{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1, -1, -1, -1},
	{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1, -1, -1, -1},
	{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1, -1, -1, -1},
	{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1, -1, -1, -1},
	{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1, -1, -1, -1},
	{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1, -1, -1, -1},
	{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1, -1, -1, -1},
// 	{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},	// !!!!!!!!!!!!!! 1011 1110 (6, 18)
	{3, 8, 6, 5, 6, 8, 8, 0, 5, 10, 5, 0, 0, 3, 10, 6, 10, 3, -1},			// !!!!!!!!!!!!!! 1011 1110 (6, 26)
	{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1, -1, -1, -1},
	{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1, -1, -1, -1},
	{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1, -1, -1, -1},
	{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1, -1, -1, -1},
	{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1, -1, -1, -1},
	{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1, -1, -1, -1},
// 	{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, 	// !!!!!!!!!!!!!! 1101 0111 (6, 18)
	{9, 4, 3, 2, 11, 5, 4, 5, 11, 11, 3, 4, 5, 9, 2, 3, 2, 9, -1}, 			// !!!!!!!!!!!!!! 1101 0111 (6, 26)
	{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1, -1, -1, -1},
	{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1, -1, -1, -1},
	{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1, -1, -1, -1},
	{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1, -1, -1, -1},
	{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1, -1, -1, -1},
	{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1, -1, -1, -1},
// 	{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, 	// !!!!!!!!!!!!!! 1110 1011 (6, 18)
	{1, 10, 4, 8, 7, 2, 10, 2, 7, 7, 4, 10, 2, 1, 8, 4, 8, 1, -1}, 			// !!!!!!!!!!!!!! 1110 1011 (6, 26)
	{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1, -1, -1, -1},
	{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1, -1, -1, -1},
	{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1, -1, -1, -1},
	{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};


int PolygoniseCh (GRIDCELL& grid, TRIANGLE *triangles, int depth, bool (*inside)(double x, double y, double z))
{
	int i,ntriang;
	int cubeindex;
	int havelist[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	XYZ vertlist[12];
	XYZ q[3];
	/*
	int k,j;
	int ind[24] = { 0,3,8,
			0,1,9,
			1,2,10,
			2,3,11,
			4,7,8,
			4,5,9,
			5,6,10,
			6,7,11};
	*/
	/*
	Determine the index into the edge table which
	tells us which vertices are inside of the surface
	*/
	cubeindex = 0;
		
	
	if (grid.val[0] <= 0) cubeindex |= 1;
	if (grid.val[1] <= 0) cubeindex |= 2;
	if (grid.val[2] <= 0) cubeindex |= 4;
	if (grid.val[3] <= 0) cubeindex |= 8;
	if (grid.val[4] <= 0) cubeindex |= 16;
	if (grid.val[5] <= 0) cubeindex |= 32;
	if (grid.val[6] <= 0) cubeindex |= 64;
	if (grid.val[7] <= 0) cubeindex |= 128;

	/* Cube is entirely in/out of the surface */
	if (edgeTable[cubeindex] == 0)
		return (0);

	/* Find the vertices where the surface intersects the cube */
	if (edgeTable[cubeindex] & 1)
	{
		vertlist[0] = VertexInterpCh(grid.p[0],grid.p[1],(int)grid.val[0],(int)grid.val[1],depth,inside);
		havelist[0] = 1;
	}
	if (edgeTable[cubeindex] & 2)
	{
		vertlist[1] = VertexInterpCh(grid.p[1],grid.p[2],(int)grid.val[1],(int)grid.val[2],depth,inside);
		havelist[1] = 1;
	}
	if (edgeTable[cubeindex] & 4)
	{
		vertlist[2] = VertexInterpCh(grid.p[2],grid.p[3],(int)grid.val[2],(int)grid.val[3],depth,inside);
		havelist[2] = 1;
	}
	if (edgeTable[cubeindex] & 8)
	{
		vertlist[3] = VertexInterpCh(grid.p[3],grid.p[0],(int)grid.val[3],(int)grid.val[0],depth,inside);
		havelist[3] = 1;
	}
	if (edgeTable[cubeindex] & 16)
	{
		vertlist[4] = VertexInterpCh(grid.p[4],grid.p[5],(int)grid.val[4],(int)grid.val[5],depth,inside);
		havelist[4] = 1;
	}
	if (edgeTable[cubeindex] & 32)
	{
		vertlist[5] = VertexInterpCh(grid.p[5],grid.p[6],(int)grid.val[5],(int)grid.val[6],depth,inside);
		havelist[5] = 1;
	}
	if (edgeTable[cubeindex] & 64)
	{
		vertlist[6] = VertexInterpCh(grid.p[6],grid.p[7],(int)grid.val[6],(int)grid.val[7],depth,inside);
		havelist[6] = 1;
	}
	if (edgeTable[cubeindex] & 128)
	{
		vertlist[7] = VertexInterpCh(grid.p[7],grid.p[4],(int)grid.val[7],(int)grid.val[4],depth,inside);
		havelist[7] = 1;
	}
	if (edgeTable[cubeindex] & 256)
	{
		vertlist[8] = VertexInterpCh(grid.p[0],grid.p[4],(int)grid.val[0],(int)grid.val[4],depth,inside);
		havelist[8] = 1;
	}
	if (edgeTable[cubeindex] & 512)
	{
		vertlist[9] = VertexInterpCh(grid.p[1],grid.p[5],(int)grid.val[1],(int)grid.val[5],depth,inside);
		havelist[9] = 1;
	}
	if (edgeTable[cubeindex] & 1024)
	{
		vertlist[10] = VertexInterpCh(grid.p[2],grid.p[6],(int)grid.val[2],(int)grid.val[6],depth,inside);
		havelist[10] = 1;
	}
	if (edgeTable[cubeindex] & 2048)
	{
		vertlist[11] = VertexInterpCh(grid.p[3],grid.p[7],(int)grid.val[3],(int)grid.val[7],depth,inside);
		havelist[11] = 1;
	}
	/*
	for (i = 0 ;i < 8; i++)
	{
		for (j = 0;j<3;j++)
		{
			if (compare(grid.p[i],vertlist[ind[3*i+j]]) )
			{
				for (k=0;k<3;k++)
				{
					if (!(havelist[ind[3*i+k]]))
						continue;
					else
						vertlist[ind[3*i+k]] = grid.p[i];
				}
				for (k=0;k<3;k++)
				{
					if (!(havelist[ind[3*i+k]]))
						break;
				}
				if (k == 3) grid.val[i] = 1.;
				break;
			}
		}
	}
	*/
	/* Create the triangle */
	ntriang = 0;
	for (i = 0; triTable[cubeindex][i] != -1; i += 3)
	{
		q[0] = vertlist[triTable[cubeindex][i  ]];
		q[1] = vertlist[triTable[cubeindex][i+1]];
		q[2] = vertlist[triTable[cubeindex][i+2]];

		if ( compare(q[0],q[1]) || compare(q[0],q[2]) || compare(q[1],q[2]) ) continue;
		
		triangles[ntriang].p[0] = q[0];
		triangles[ntriang].p[1] = q[1];
		triangles[ntriang].p[2] = q[2];
		ntriang++;
	}

	return (ntriang);
}



int Polygonise (GRIDCELL grid, double isolevel, TRIANGLE *triangles)
{
   int i,ntriang;
   int cubeindex;
   XYZ vertlist[12];



   /*
      Determine the index into the edge table which
      tells us which vertices are inside of the surface
   */
   cubeindex = 0;
   if (grid.val[0] < isolevel) cubeindex |= 1;
   if (grid.val[1] < isolevel) cubeindex |= 2;
   if (grid.val[2] < isolevel) cubeindex |= 4;
   if (grid.val[3] < isolevel) cubeindex |= 8;
   if (grid.val[4] < isolevel) cubeindex |= 16;
   if (grid.val[5] < isolevel) cubeindex |= 32;
   if (grid.val[6] < isolevel) cubeindex |= 64;
   if (grid.val[7] < isolevel) cubeindex |= 128;

   /* Cube is entirely in/out of the surface */
   if (edgeTable[cubeindex] == 0)
      return (0);

   /* Find the vertices where the surface intersects the cube */
   if (edgeTable[cubeindex] & 1)
      vertlist[0] =
         VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
   if (edgeTable[cubeindex] & 2)
      vertlist[1] =
         VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
   if (edgeTable[cubeindex] & 4)
      vertlist[2] =
         VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
   if (edgeTable[cubeindex] & 8)
      vertlist[3] =
         VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
   if (edgeTable[cubeindex] & 16)
      vertlist[4] =
         VertexInterp(isolevel,grid.p[4],grid.p[5],grid.val[4],grid.val[5]);
   if (edgeTable[cubeindex] & 32)
      vertlist[5] =
         VertexInterp(isolevel,grid.p[5],grid.p[6],grid.val[5],grid.val[6]);
   if (edgeTable[cubeindex] & 64)
      vertlist[6] =
         VertexInterp(isolevel,grid.p[6],grid.p[7],grid.val[6],grid.val[7]);
   if (edgeTable[cubeindex] & 128)
      vertlist[7] =
         VertexInterp(isolevel,grid.p[7],grid.p[4],grid.val[7],grid.val[4]);
   if (edgeTable[cubeindex] & 256)
      vertlist[8] =
         VertexInterp(isolevel,grid.p[0],grid.p[4],grid.val[0],grid.val[4]);
   if (edgeTable[cubeindex] & 512)
      vertlist[9] =
         VertexInterp(isolevel,grid.p[1],grid.p[5],grid.val[1],grid.val[5]);
   if (edgeTable[cubeindex] & 1024)
      vertlist[10] =
         VertexInterp(isolevel,grid.p[2],grid.p[6],grid.val[2],grid.val[6]);
   if (edgeTable[cubeindex] & 2048)
      vertlist[11] =
         VertexInterp(isolevel,grid.p[3],grid.p[7],grid.val[3],grid.val[7]);

   /* Create the triangle */
   ntriang = 0;
   for (i = 0; triTable[cubeindex][i] != -1; i += 3)
   {
      triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i  ]];
      triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i+1]];
      triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i+2]];
      ntriang++;
   }

   return (ntriang);
}


double DistTriPoint (TRIANGLE tr, XYZ p)
{
	double A, B, C, D, t, h, k, l;
	XYZ P, Q1, Q2, Q3;
	A = (tr.p[1].y - tr.p[0].y) * (tr.p[2].z - tr.p[0].z) - (tr.p[1].z - tr.p[0].z) * (tr.p[2].y - tr.p[0].y);
	B = (tr.p[1].z - tr.p[0].z) * (tr.p[2].x - tr.p[0].x) - (tr.p[1].x - tr.p[0].x) * (tr.p[2].z - tr.p[0].z);
	C = (tr.p[1].x - tr.p[0].x) * (tr.p[2].y - tr.p[0].y) - (tr.p[1].y - tr.p[0].y) * (tr.p[2].x - tr.p[0].x);
	D = - (A * tr.p[0].x + B * tr.p[0].y + C * tr.p[0].z);
	
	t = ((tr.p[0].x - p.x) * A + (tr.p[0].y - p.y) * B + (tr.p[0].z - p.z) * C) / (A * A + B * B + C * C);
	
	P.x = p.x + t * A;
	P.y = p.y + t * B;
	P.z = p.z + t * C;
	
	Q1.x = (P.y - tr.p[0].y) * (P.z - tr.p[1].z) - (P.y - tr.p[1].y) * (P.z - tr.p[0].z);
	Q1.y = (P.z - tr.p[0].z) * (P.x - tr.p[1].x) - (P.z - tr.p[1].z) * (P.x - tr.p[0].x);
	Q1.z = (P.x - tr.p[0].x) * (P.y - tr.p[1].y) - (P.x - tr.p[1].x) * (P.y - tr.p[0].y);
		
	Q2.x = (P.y - tr.p[1].y) * (P.z - tr.p[2].z) - (P.y - tr.p[2].y) * (P.z - tr.p[1].z);
	Q2.y = (P.z - tr.p[1].z) * (P.x - tr.p[2].x) - (P.z - tr.p[2].z) * (P.x - tr.p[1].x);
	Q2.z = (P.x - tr.p[1].x) * (P.y - tr.p[2].y) - (P.x - tr.p[2].x) * (P.y - tr.p[1].y);
	
	Q3.x = (P.y - tr.p[2].y) * (P.z - tr.p[0].z) - (P.y - tr.p[0].y) * (P.z - tr.p[2].z);
	Q3.y = (P.z - tr.p[2].z) * (P.x - tr.p[0].x) - (P.z - tr.p[0].z) * (P.x - tr.p[2].x);
	Q3.z = (P.x - tr.p[2].x) * (P.y - tr.p[0].y) - (P.x - tr.p[0].x) * (P.y - tr.p[2].y);
	
	if (Q1.x * Q2.x + Q1.y * Q2.y + Q1.z * Q2.z >= 0 &&
	    Q2.x * Q3.x + Q2.y * Q3.y + Q2.z * Q3.z >= 0 &&
	    Q1.x * Q3.x + Q1.y * Q3.y + Q1.z * Q3.z >= 0)
		return fabs ((A * p.x + B * p.y + C * p.z + D) / sqrt (A * A + B * B + C * C));


	t = ((tr.p[1].x - tr.p[0].x) * (p.x - tr.p[0].x) + 
	     (tr.p[1].y - tr.p[0].y) * (p.y - tr.p[0].y) + 
	     (tr.p[1].z - tr.p[0].z) * (p.z - tr.p[0].z)) / 
	    ((tr.p[1].x - tr.p[0].x) * (tr.p[1].x - tr.p[0].x) + 
	     (tr.p[1].y - tr.p[0].y) * (tr.p[1].y - tr.p[0].y) + 
	     (tr.p[1].z - tr.p[0].z) * (tr.p[1].z - tr.p[0].z));
	
	if (t >= 0.0 && t <= 1.0)
		h = sqrt ((tr.p[0].x - p.x + t * (tr.p[1].x - tr.p[0].x)) * (tr.p[0].x - p.x + t * (tr.p[1].x - tr.p[0].x))
			+ (tr.p[0].y - p.y + t * (tr.p[1].y - tr.p[0].y)) * (tr.p[0].y - p.y + t * (tr.p[1].y - tr.p[0].y))
			+ (tr.p[0].z - p.z + t * (tr.p[1].z - tr.p[0].z)) * (tr.p[0].z - p.z + t * (tr.p[1].z - tr.p[0].z)));
	else if (t < 0.0)
		h = sqrt ((tr.p[0].x - p.x) * (tr.p[0].x - p.x)
			+ (tr.p[0].y - p.y) * (tr.p[0].y - p.y)
			+ (tr.p[0].z - p.z) * (tr.p[0].z - p.z));
	else
		h = sqrt ((tr.p[1].x - p.x) * (tr.p[1].x - p.x)
			+ (tr.p[1].y - p.y) * (tr.p[1].y - p.y)
			+ (tr.p[1].z - p.z) * (tr.p[1].z - p.z));
	
	t = ((tr.p[2].x - tr.p[0].x) * (p.x - tr.p[0].x) + 
	     (tr.p[2].y - tr.p[0].y) * (p.y - tr.p[0].y) + 
	     (tr.p[2].z - tr.p[0].z) * (p.z - tr.p[0].z)) / 
	    ((tr.p[2].x - tr.p[0].x) * (tr.p[2].x - tr.p[0].x) + 
	     (tr.p[2].y - tr.p[0].y) * (tr.p[2].y - tr.p[0].y) + 
	     (tr.p[2].z - tr.p[0].z) * (tr.p[2].z - tr.p[0].z));
		
	if (t >= 0.0 && t <= 1.0)
		k = sqrt ((tr.p[0].x - p.x + t * (tr.p[2].x - tr.p[0].x)) * (tr.p[0].x - p.x + t * (tr.p[2].x - tr.p[0].x))
			+ (tr.p[0].y - p.y + t * (tr.p[2].y - tr.p[0].y)) * (tr.p[0].y - p.y + t * (tr.p[2].y - tr.p[0].y))
			+ (tr.p[0].z - p.z + t * (tr.p[2].z - tr.p[0].z)) * (tr.p[0].z - p.z + t * (tr.p[2].z - tr.p[0].z)));
	else if (t < 0.0)
		k = sqrt ((tr.p[0].x - p.x) * (tr.p[0].x - p.x)
			+ (tr.p[0].y - p.y) * (tr.p[0].y - p.y)
			+ (tr.p[0].z - p.z) * (tr.p[0].z - p.z));
	else
		k = sqrt ((tr.p[2].x - p.x) * (tr.p[2].x - p.x)
			+ (tr.p[2].y - p.y) * (tr.p[2].y - p.y)
			+ (tr.p[2].z - p.z) * (tr.p[2].z - p.z));
	
	t = ((tr.p[2].x - tr.p[1].x) * (p.x - tr.p[1].x) + 
	     (tr.p[2].y - tr.p[1].y) * (p.y - tr.p[1].y) + 
	     (tr.p[2].z - tr.p[1].z) * (p.z - tr.p[1].z)) / 
	    ((tr.p[2].x - tr.p[1].x) * (tr.p[2].x - tr.p[1].x) + 
	     (tr.p[2].y - tr.p[1].y) * (tr.p[2].y - tr.p[1].y) + 
	     (tr.p[2].z - tr.p[1].z) * (tr.p[2].z - tr.p[1].z));	
		
	if (t >= 0.0 && t <= 1.0)
		l = sqrt ((tr.p[1].x - p.x + t * (tr.p[2].x - tr.p[1].x)) * (tr.p[1].x - p.x + t * (tr.p[2].x - tr.p[1].x))
			+ (tr.p[1].y - p.y + t * (tr.p[2].y - tr.p[1].y)) * (tr.p[1].y - p.y + t * (tr.p[2].y - tr.p[1].y))
			+ (tr.p[1].z - p.z + t * (tr.p[2].z - tr.p[1].z)) * (tr.p[1].z - p.z + t * (tr.p[2].z - tr.p[1].z)));
	else if (t < 0.0)
		l = sqrt ((tr.p[1].x - p.x) * (tr.p[1].x - p.x)
			+ (tr.p[1].y - p.y) * (tr.p[1].y - p.y)
			+ (tr.p[1].z - p.z) * (tr.p[1].z - p.z));
	else
		l = sqrt ((tr.p[2].x - p.x) * (tr.p[2].x - p.x)
			+ (tr.p[2].y - p.y) * (tr.p[2].y - p.y)
			+ (tr.p[2].z - p.z) * (tr.p[2].z - p.z));

	if (h <= k && h <= l)
		return h;
		
	if (k <= h && k <= l)
		return k;
		
	return l;
}
#define SQR(x) ((double)(x)*(x))

double DistPointPoint(XYZ p1, XYZ p2)
{
	return sqrt(SQR(p1.x-p2.x)+SQR(p1.y-p2.y)+SQR(p1.z-p2.z));
}

double TriArea(TRIANGLE tr)
{
	double a,b,c,s;
    a = DistPointPoint(tr.p[0],tr.p[1]);
    b = DistPointPoint(tr.p[1],tr.p[2]);
    c = DistPointPoint(tr.p[2],tr.p[0]);
    s = (a+b+c)/2.0f;
    return sqrt(s*(s-a)*(s-b)*(s-c));
}

int PolygoniseTri(GRIDCELL g,double iso,
   TRIANGLE *tri,int v0,int v1,int v2,int v3)
{
   int ntri = 0;
   int triindex;

   /*
      Determine which of the 16 cases we have given which vertices
      are above or below the isosurface
   */
   triindex = 0;
   if (g.val[v0] < iso) triindex |= 1;
   if (g.val[v1] < iso) triindex |= 2;
   if (g.val[v2] < iso) triindex |= 4;
   if (g.val[v3] < iso) triindex |= 8;

   /* Form the vertices of the triangles for each case */
   switch (triindex) {
   case 0x00:
   case 0x0F:
      break;
   case 0x0E:
   case 0x01:
      tri[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      tri[0].p[1] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      tri[0].p[2] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      ntri++;
      break;
   case 0x0D:
   case 0x02:
      tri[0].p[0] = VertexInterp(iso,g.p[v1],g.p[v0],g.val[v1],g.val[v0]);
      tri[0].p[1] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      tri[0].p[2] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      ntri++;
      break;
   case 0x0C:
   case 0x03:
      tri[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      tri[0].p[1] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      tri[0].p[2] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      ntri++;
      tri[1].p[0] = tri[0].p[2];
      tri[1].p[1] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      tri[1].p[2] = tri[0].p[1];
      ntri++;
      break;
   case 0x0B:
   case 0x04:
      tri[0].p[0] = VertexInterp(iso,g.p[v2],g.p[v0],g.val[v2],g.val[v0]);
      tri[0].p[1] = VertexInterp(iso,g.p[v2],g.p[v1],g.val[v2],g.val[v1]);
      tri[0].p[2] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      ntri++;
      break;
   case 0x0A:
   case 0x05:
      tri[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      tri[0].p[1] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      tri[0].p[2] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      ntri++;
      tri[1].p[0] = tri[0].p[0];
      tri[1].p[1] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      tri[1].p[2] = tri[0].p[1];
      ntri++;
      break;
   case 0x09:
   case 0x06:
      tri[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      tri[0].p[1] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      tri[0].p[2] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      ntri++;
      tri[1].p[0] = tri[0].p[0];
      tri[1].p[1] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      tri[1].p[2] = tri[0].p[2];
      ntri++;
      break;
   case 0x07:
   case 0x08:
      tri[0].p[0] = VertexInterp(iso,g.p[v3],g.p[v0],g.val[v3],g.val[v0]);
      tri[0].p[1] = VertexInterp(iso,g.p[v3],g.p[v2],g.val[v3],g.val[v2]);
      tri[0].p[2] = VertexInterp(iso,g.p[v3],g.p[v1],g.val[v3],g.val[v1]);
      ntri++;
      break;
   }

   return(ntri);
}

double PolygoniseTriVolume(GRIDCELL g,double iso,int v0,int v1,int v2,int v3)
{
   TRIANGLE tri[4];
   int triindex;
   int q1,q2;
   double vol = 0.0, V;
   /*
      Determine which of the 16 cases we have given which vertices
      are above or below the isosurface
   */
   triindex = 0;
   if (g.val[v0] < iso) triindex |= 1;
   if (g.val[v1] < iso) triindex |= 2;
   if (g.val[v2] < iso) triindex |= 4;
   if (g.val[v3] < iso) triindex |= 8;
   
   tri[3].p[0] = g.p[v1];
   tri[3].p[1] = g.p[v2];
   tri[3].p[2] = g.p[v3];
   
   V = TriArea(tri[3])*DistTriPoint(tri[3],g.p[v0])/3.0;

   /* Form the vertices of the triangles for each case */
   switch (triindex) {
   case 0x00:
   case 0x0F:
      break;
   case 0x0E:
   case 0x01:
      tri[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      tri[0].p[1] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      tri[0].p[2] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
	  if( g.val[v0] < iso ) vol = TriArea(tri[0])*DistTriPoint(tri[0],g.p[v0])/3.0;
	  else vol = V - TriArea(tri[0])*DistTriPoint(tri[0],g.p[v0])/3.0;
      break;
   case 0x0D:
   case 0x02:
      tri[0].p[0] = VertexInterp(iso,g.p[v1],g.p[v0],g.val[v1],g.val[v0]);
      tri[0].p[1] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      tri[0].p[2] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
	  if( g.val[v1] < iso ) vol = TriArea(tri[0])*DistTriPoint(tri[0],g.p[v1])/3.0;
	  else vol = V - TriArea(tri[0])*DistTriPoint(tri[0],g.p[v1])/3.0;
      break;
   case 0x0C:
   case 0x03:
      tri[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      tri[0].p[1] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      tri[0].p[2] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      tri[1].p[0] = tri[0].p[2];
      tri[1].p[1] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      tri[1].p[2] = tri[0].p[1];
	  if( g.val[v0] < iso ) {q1 = v0; q2 = v1;}
	  else  {q1 = v3; q2 = v2;}
	  tri[2].p[0] = g.p[q1];
	  tri[2].p[1] = g.p[q2];
	  tri[2].p[2] = tri[0].p[1];
	  vol  = DistTriPoint(tri[0],g.p[q1])*TriArea(tri[0])/3.0;
	  vol += DistTriPoint(tri[1],g.p[q2])*TriArea(tri[1])/3.0;
	  vol += DistTriPoint(tri[2],tri[0].p[2])*TriArea(tri[2])/3.0;
      break;
   case 0x0B:
   case 0x04:
      tri[0].p[0] = VertexInterp(iso,g.p[v2],g.p[v0],g.val[v2],g.val[v0]);
      tri[0].p[1] = VertexInterp(iso,g.p[v2],g.p[v1],g.val[v2],g.val[v1]);
      tri[0].p[2] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
	  if( g.val[v2] < iso ) vol = TriArea(tri[0])*DistTriPoint(tri[0],g.p[v2])/3.0;
	  else vol = V - TriArea(tri[0])*DistTriPoint(tri[0],g.p[v2])/3.0;
      break;
   case 0x0A:
   case 0x05:
      tri[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      tri[0].p[1] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      tri[0].p[2] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      tri[1].p[0] = tri[0].p[0];
      tri[1].p[1] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      tri[1].p[2] = tri[0].p[1];
	  if( g.val[v0] < iso ) {q1 = v0; q2 = v2;}
	  else  {q1 = v3; q2 = v1;}
	  tri[2].p[0] = g.p[q1];
	  tri[2].p[1] = g.p[q2];
	  tri[2].p[2] = tri[0].p[0];
	  vol  = DistTriPoint(tri[0],g.p[q1])*TriArea(tri[0])/3.0;
	  vol += DistTriPoint(tri[1],g.p[q2])*TriArea(tri[1])/3.0;
	  vol += DistTriPoint(tri[2],tri[0].p[1])*TriArea(tri[2])/3.0;
      break;
   case 0x09:
   case 0x06:
      tri[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      tri[0].p[1] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      tri[0].p[2] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      tri[1].p[0] = tri[0].p[1];
      tri[1].p[1] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      tri[1].p[2] = tri[0].p[2];
	  if( g.val[v0] < iso ) {q1 = v0; q2 = v3;}
	  else  {q1 = v1; q2 = v2;}
	  tri[2].p[0] = g.p[q1];
	  tri[2].p[1] = g.p[q2];
	  tri[2].p[2] = tri[0].p[1];
	  vol  = DistTriPoint(tri[0],g.p[q1])*TriArea(tri[0])/3.0;
	  vol += DistTriPoint(tri[1],g.p[q2])*TriArea(tri[1])/3.0;
	  vol += DistTriPoint(tri[2],tri[0].p[2])*TriArea(tri[2])/3.0;
      break;
   case 0x07:
   case 0x08:
      tri[0].p[0] = VertexInterp(iso,g.p[v3],g.p[v0],g.val[v3],g.val[v0]);
      tri[0].p[1] = VertexInterp(iso,g.p[v3],g.p[v2],g.val[v3],g.val[v2]);
      tri[0].p[2] = VertexInterp(iso,g.p[v3],g.p[v1],g.val[v3],g.val[v1]);
	  if( g.val[v3] < iso ) vol = TriArea(tri[0])*DistTriPoint(tri[0],g.p[v3])/3.0;
	  else vol = V - TriArea(tri[0])*DistTriPoint(tri[0],g.p[v3])/3.0;
      break;
   }

   return vol;
}
