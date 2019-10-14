#ifndef __SLICE__H
#define __SLICE__H
#include "inmost.h"

using namespace INMOST;

class Slice
{
	double epsf, epsl;
	int maxits, maxits_zero;
	Mesh * m;
	double Search(double r0, double r1, double c0[3], double c1[3], double p[3], bool binary = false) const;
	double SearchZero(double r0, double r1, double c0[3], double c1[3], double p[3]) const;
public:
	Slice(double epsf = 1.0e-6, double epsl = 1.0e-3, int maxits = 100, int maxits_zero = 20) : epsf(epsf), epsl(epsl), maxits(maxits), maxits_zero(maxits_zero) {}
	Slice(const Slice & b) : epsf(b.epsf), epsl(b.epsl), maxits(b.maxits), maxits_zero(b.maxits_zero) {}
	Slice & operator =(Slice const & b) {epsf = b.epsf; epsl = b.epsl; maxits = b.maxits; maxits_zero = maxits_zero; return *this;}
	virtual ~Slice() {}
	virtual double LevelFunction(double p[3]) const = 0;
	void SliceMesh(Mesh & m, bool remove_material_zero = true);
};
#endif // __SLICE__H


