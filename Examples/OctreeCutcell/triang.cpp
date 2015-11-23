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
