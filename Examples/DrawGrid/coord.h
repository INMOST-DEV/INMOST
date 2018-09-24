#ifndef _COORD_H
#define _COORD_H

#include <cmath>

class coord
{
	double p[3];
public:
	coord() { p[0] = p[1] = p[2] = 0; }
	coord(double xyz[3]) { p[0] = xyz[0]; p[1] = xyz[1]; p[2] = xyz[2]; }
	coord(double x, double y, double z) { p[0] = x; p[1] = y; p[2] = z; }
	coord(const coord & other) { p[0] = other.p[0]; p[1] = other.p[1]; p[2] = other.p[2]; }
	coord & operator = (coord const & other) { p[0] = other.p[0]; p[1] = other.p[1]; p[2] = other.p[2]; return *this; }
	coord & operator +=(const coord & other) { p[0] += other.p[0]; p[1] += other.p[1]; p[2] += other.p[2]; return *this; }
	coord & operator -=(const coord & other) { p[0] -= other.p[0]; p[1] -= other.p[1]; p[2] -= other.p[2]; return *this; }
	coord & operator *=(const coord & other)
	{
		double tmp[3] = { p[1] * other.p[2] - p[2] * other.p[1], p[2] * other.p[0] - p[0] * other.p[2], p[0] * other.p[1] - p[1] * other.p[0] };
		p[0] = tmp[0]; p[1] = tmp[1]; p[2] = tmp[2];
		return *this;
	}
	coord & operator *=(double other) { p[0] *= other; p[1] *= other; p[2] *= other; return *this; }
	coord & operator /=(double other) { p[0] /= other; p[1] /= other; p[2] /= other; return *this; }
	coord operator -(const coord & other) const { return coord(p[0] - other.p[0], p[1] - other.p[1], p[2] - other.p[2]); }
	coord operator +(const coord & other) const { return coord(p[0] + other.p[0], p[1] + other.p[1], p[2] + other.p[2]); }
	coord operator *(const coord & other) const { return coord(p[1] * other.p[2] - p[2] * other.p[1], p[2] * other.p[0] - p[0] * other.p[2], p[0] * other.p[1] - p[1] * other.p[0]); }
	coord operator /(double other) const { return coord(p[0] / other, p[1] / other, p[2] / other); }
	coord operator *(double other) const { return coord(p[0] * other, p[1] * other, p[2] * other); }
	double operator ^(const coord & other) const { return p[0] * other.p[0] + p[1] * other.p[1] + p[2] * other.p[2]; }
	~coord() {}
	double length() const { return sqrt((*this) ^ (*this)); }
	double & operator [](int i) { return p[i]; }
	double operator [](int i) const { return p[i]; }
	double * data() { return p; }
	const double * data()const { return p; }
};

double abs(const coord & p);
void get_matrix(const coord & a, const coord & b, double matrix[16]);


#endif