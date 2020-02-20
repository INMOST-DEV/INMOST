#ifndef _FIX_FAULTS_H
#define _FIX_FAULTS_H

#include "inmost.h"

using namespace INMOST;

class FixFaults
{
	Mesh & m;
public:
	FixFaults(Mesh & m) : m(m) {}
	FixFaults(const FixFaults & b) : m(b.m) {}
	FixFaults & operator =(FixFaults const & b) {m = b.m; return *this;}
	~FixFaults() {}
	
	void FixMeshFaults();
};

#endif //_FIX_FAULTS_H
