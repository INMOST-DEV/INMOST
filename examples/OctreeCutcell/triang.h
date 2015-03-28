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


double DistTriPoint (TRIANGLE tr, XYZ p);
double TriArea(TRIANGLE tr);


#endif // _TRIANG_H_
