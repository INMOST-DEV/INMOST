#ifndef __SLICE__H
#define __SLICE__H
#include "inmost.h"

using namespace INMOST;

class Slice
{
	double epsf, epsl;
	int maxits;
	Mesh * m;
	double Search(double r0, double r1, double c0[3], double c1[3], double p[3], bool binary = false) const;
	double SearchZero(double r0, double r1, double c0[3], double c1[3], double p[3]) const;
public:
	Slice(double epsf = 1.0e-6, double epsl = 1.0e-3, int maxits = 100) : epsf(epsf), epsl(epsl), maxits(maxits) {}
	Slice(const Slice & b) : epsf(b.epsf), epsl(b.epsl), maxits(b.maxits) {}
	Slice & operator =(Slice const & b) {epsf = b.epsf; epsl = b.epsl; maxits = b.maxits; return *this;}
	virtual ~Slice();
	virtual double LevelFunction(double p[3]) const = 0;
	void SliceMesh(Mesh & m);
};
#endif // __SLICE__H


