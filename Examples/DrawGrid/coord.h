#ifndef _COORD_H
#define _COORD_H

#include <cmath>
#include "inmost_common.h"
class coord
{
	INMOST_DATA_REAL_TYPE p[3];
public:
	coord() { p[0] = p[1] = p[2] = 0; }
	coord(INMOST_DATA_REAL_TYPE xyz[3]) { p[0] = xyz[0]; p[1] = xyz[1]; p[2] = xyz[2]; }
	coord(INMOST_DATA_REAL_TYPE x, INMOST_DATA_REAL_TYPE y, INMOST_DATA_REAL_TYPE z) { p[0] = x; p[1] = y; p[2] = z; }
	coord(const coord & other) { p[0] = other.p[0]; p[1] = other.p[1]; p[2] = other.p[2]; }
	coord & operator = (coord const & other) { p[0] = other.p[0]; p[1] = other.p[1]; p[2] = other.p[2]; return *this; }
	coord & operator +=(const coord & other) { p[0] += other.p[0]; p[1] += other.p[1]; p[2] += other.p[2]; return *this; }
	coord & operator -=(const coord & other) { p[0] -= other.p[0]; p[1] -= other.p[1]; p[2] -= other.p[2]; return *this; }
	coord & operator *=(const coord & other)
	{
		INMOST_DATA_REAL_TYPE tmp[3] = { p[1] * other.p[2] - p[2] * other.p[1], p[2] * other.p[0] - p[0] * other.p[2], p[0] * other.p[1] - p[1] * other.p[0] };
		p[0] = tmp[0]; p[1] = tmp[1]; p[2] = tmp[2];
		return *this;
	}
	coord & operator *=(INMOST_DATA_REAL_TYPE other) { p[0] *= other; p[1] *= other; p[2] *= other; return *this; }
	coord & operator /=(INMOST_DATA_REAL_TYPE other) { p[0] /= other; p[1] /= other; p[2] /= other; return *this; }
	coord operator -(const coord & other) const { return coord(p[0] - other.p[0], p[1] - other.p[1], p[2] - other.p[2]); }
	coord operator +(const coord & other) const { return coord(p[0] + other.p[0], p[1] + other.p[1], p[2] + other.p[2]); }
	coord operator *(const coord & other) const { return coord(p[1] * other.p[2] - p[2] * other.p[1], p[2] * other.p[0] - p[0] * other.p[2], p[0] * other.p[1] - p[1] * other.p[0]); }
	coord operator /(INMOST_DATA_REAL_TYPE other) const { return coord(p[0] / other, p[1] / other, p[2] / other); }
	coord operator *(INMOST_DATA_REAL_TYPE other) const { return coord(p[0] * other, p[1] * other, p[2] * other); }
	INMOST_DATA_REAL_TYPE operator ^(const coord & other) const { return p[0] * other.p[0] + p[1] * other.p[1] + p[2] * other.p[2]; }
	~coord() {}
	INMOST_DATA_REAL_TYPE length() const { return sqrt((*this) ^ (*this)); }
	INMOST_DATA_REAL_TYPE & operator [](int i) { return p[i]; }
	INMOST_DATA_REAL_TYPE operator [](int i) const { return p[i]; }
	INMOST_DATA_REAL_TYPE * data() { return p; }
	const INMOST_DATA_REAL_TYPE * data()const { return p; }
};

INMOST_DATA_REAL_TYPE abs(const coord & p);
void get_matrix(const coord & a, const coord & b, INMOST_DATA_REAL_TYPE matrix[16]);


#endif
