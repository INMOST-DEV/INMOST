#include "coord.h"



INMOST_DATA_REAL_TYPE abs(const coord & p)
{
	return sqrt(p^p);
}


void get_matrix(const coord & a, const coord & b, INMOST_DATA_REAL_TYPE matrix[16])
{
	INMOST_DATA_REAL_TYPE d;
	coord z = (b - a) / sqrt((b - a) ^ (b - a));
	coord y;
	coord x;
	y = coord(z[1], -z[2], 0);
	d = sqrt(y^y);
	if (d < 1e-5)
	{
		y = coord(-z[2], 0, z[0]);
		d = sqrt(y^y);
	}
	y = y / d;
	x = y*z;
	x = x / sqrt(x^x);
	y = x*z;
	matrix[0] = x[0];
	matrix[1] = x[1];
	matrix[2] = x[2];
	matrix[3] = 0;
	matrix[4] = y[0];
	matrix[5] = y[1];
	matrix[6] = y[2];
	matrix[7] = 0;
	matrix[8] = z[0];
	matrix[9] = z[1];
	matrix[10] = z[2];
	matrix[11] = 0;
	matrix[12] = 0;
	matrix[13] = 0;
	matrix[14] = 0;
	matrix[15] = 1;
}
