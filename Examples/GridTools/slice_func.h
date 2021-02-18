#ifndef __SLICE__H
#define __SLICE__H
#include "inmost.h"

using namespace INMOST;

class Slice
{
	Storage::real epsf, epsl;
	int maxits, maxits_zero;
	Mesh * m;
	Storage::real Search(Storage::real r0, Storage::real r1, Storage::real c0[3], Storage::real c1[3], Storage::real p[3], bool binary = false) const;
	Storage::real SearchZero(Storage::real r0, Storage::real r1, Storage::real c0[3], Storage::real c1[3], Storage::real p[3]) const;
public:
	Slice(Storage::real epsf = 1.0e-6, Storage::real epsl = 1.0e-3, int maxits = 100, int maxits_zero = 20) : epsf(epsf), epsl(epsl), maxits(maxits), maxits_zero(maxits_zero) {}
	Slice(const Slice & b) : epsf(b.epsf), epsl(b.epsl), maxits(b.maxits), maxits_zero(b.maxits_zero) {}
	Slice & operator =(Slice const & b) {epsf = b.epsf; epsl = b.epsl; maxits = b.maxits; maxits_zero = maxits_zero; return *this;}
	virtual ~Slice() {}
	virtual Storage::real LevelFunction(Storage::real p[3]) const = 0;
	void SliceMesh(Mesh & m, bool remove_material_zero = true);
};
#endif // __SLICE__H


