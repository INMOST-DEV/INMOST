#if !defined(_TRIANG_H_)
#define _TRIANG_H_

typedef struct
{
	double x;
	double y;
	double z;
} XYZ;

typedef struct
{
	XYZ p[3];
	int n;
} TRIANGLE;

typedef struct
{
	XYZ p[8];
	double val[8];
} GRIDCELL;

XYZ VertexInterp (double isolevel, XYZ p1, XYZ p2, double valp1, double valp2);
int Polygonise (GRIDCELL grid, double isolevel, TRIANGLE *triangles);
double DistTriPoint (TRIANGLE tr, XYZ p);
double TriArea(TRIANGLE tr);
int PolygoniseTri(GRIDCELL g,double iso,TRIANGLE *tri,int v0,int v1,int v2,int v3);
double PolygoniseTriVolume(GRIDCELL g,double iso,int v0,int v1,int v2,int v3);



XYZ VertexInterpCh (XYZ p1, XYZ p2, int depth,bool (*inside)(double x, double y, double z));
int PolygoniseCh (GRIDCELL& grid, TRIANGLE *triangles, int depth, bool (*inside)(double x, double y, double z));

#endif // _TRIANG_H_
