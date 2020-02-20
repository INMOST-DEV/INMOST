#ifndef _SET_LAYERS_H
#define _SET_LAYERS_H

#include "inmost.h"

using namespace INMOST;

class SetLayers
{
	Mesh & mesh;
	int rand_seed;
	
	std::vector<double> layers_z;
public:
	SetLayers(Mesh & m, int seed = 0) :mesh(m), rand_seed(seed) {srand(rand_seed);}
	SetLayers(SetLayers const & b) :mesh(b.mesh) {}
	SetLayers & operator =(const SetLayers & b) {mesh = b.mesh; return *this;}
	~SetLayers() {}
	
	
	void RandomLayersZ(int num_layers);
	// @param set_layers_z contains ascending values from 0 to 1.
	void SetLayersZ(std::vector<double> set_layers_z) {layers_z = set_layers_z;}
	
	void DeformLayers(double coef);
	
	void DeformHeight(double height);
};

#endif //_SET_LAYERS_H
