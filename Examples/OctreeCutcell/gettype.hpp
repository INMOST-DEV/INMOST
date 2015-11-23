#ifndef _GETTYPE_H
#define _GETTYPE_H
#define GEPS 5.0e-10
#define GEPS2 0.25e-9
#define GUNDEF 1e20
#include <vector>
#include <map>
#include "my_glut.h"
#include "obj.h"


class GetType
{
private:
	std::vector<int> objs;
	std::vector<int> min_layers;
	std::vector< std::pair<double,int> > tlayers;
public:
	GetType() {}
	GetType(const GetType & other) : objs(other.objs) {}
	GetType & operator =(GetType const & other) {objs = other.objs; return *this;}
	~GetType() {}
	void ReadLayers(std::vector<std::string> & input)
	{
		for(int i = 0; i < input.size(); i++)
			objs.push_back(ReadObj((char *)input[i].c_str()));
		min_layers.resize(objs.size());
		tlayers.resize(objs.size());
	}
	int Size() { return objs.size(); }
	int Layers(double xyz[3], int * layers) 
	{
		int num_ret = 0, num_min = 0;
		double pos[3];
		pos[0] = xyz[0];
		pos[1] = xyz[1];
		pos[2] = xyz[2];//-GEPS*0.5;
		if( pos[0] < 0.0+GEPS ) pos[0] = 0.0+GEPS;
		if( pos[0] > 1.0-GEPS ) pos[0] = 1.0-GEPS;
		if( pos[1] < 0.0+GEPS ) pos[1] = 0.0+GEPS;
		if( pos[1] > 1.0-GEPS ) pos[1] = 1.0-GEPS;
		double ray[3];
		ray[0] = 0;
		ray[1] = 0;
		ray[2] = 1;
		double current_dist = GUNDEF;
		for(int i = 0; i < objs.size(); i++)
		{
			
			double dist = RayObjIntersection(pos,ray,objs[i],current_dist+GEPS*2);
			if( fabs(dist - GUNDEF) > GEPS && (dist < current_dist-GEPS*0.25)) 
			{ 
				if( dist < GEPS2 )
					min_layers[num_min++] = i;
				else 
				{
					num_ret = 1;
					layers[0] = i;
					current_dist = dist;
				}
			}
		}
		
		//if( num_min > 0 ) num_ret = 0;
		
		{
			ray[2] = -1;
			pos[2] = xyz[2];// + GEPS*0.5;
			int kmin = objs.size();
			if( num_ret+num_min > 0 )
				current_dist = GEPS;
			for(int i = 0; i < objs.size(); i++)
			{
				double dist = RayObjIntersection(pos,ray,objs[i],current_dist);
				if( dist < current_dist )
				{
					current_dist = dist;
					kmin = i;
				}
			}
			if( kmin != objs.size() && (((num_ret+num_min) == 0) || (current_dist < GEPS)) )
				layers[num_ret++] = kmin;
				
			
		}
		
		for(int i = 0; i < num_min; i++)
			layers[num_ret++] = min_layers[i];
			
		if( num_ret == 0 ) throw -1;
		
		return num_ret;
	}
	bool Normal(double xyz[3], double nrm[3])
	{
		double lnrm[3], unrm[3], lh,uh;
		double pos[3], ray[3];
		pos[0] = xyz[0];
		pos[1] = xyz[1];
		pos[2] = -1;
		if( pos[0] < 0.0+GEPS ) pos[0] = 0.0+GEPS;
		if( pos[0] > 1.0-GEPS ) pos[0] = 1.0-GEPS;
		if( pos[1] < 0.0+GEPS ) pos[1] = 0.0+GEPS;
		if( pos[1] > 1.0-GEPS ) pos[1] = 1.0-GEPS;
		ray[0] = 0;
		ray[1] = 0;
		ray[2] = 1;
		double current_dist = GUNDEF;
		int k = 0;
		for(int i = 0; i < objs.size(); i++)
		{	
			double dist = RayObjIntersection(pos,ray,objs[i],GUNDEF);
			if( fabs(dist - GUNDEF) > GEPS ) 
				tlayers[k++] = std::make_pair(dist-1,i);
		}
		std::sort(tlayers.begin(),tlayers.begin()+k);
		std::vector< std::pair<double, int> >::iterator u = std::upper_bound(tlayers.begin(),tlayers.begin()+k, std::make_pair(xyz[2]-4*GEPS,0)), l;
		if( u == tlayers.begin()+k )
		{
			if( fabs(xyz[2] - (tlayers.begin()+k-1)->first) < 1e-4 )
				u = tlayers.begin()+k-1;
		}
		if( u == tlayers.begin()+k )
			return false;
		if( u == tlayers.begin() )
		{
			lnrm[0] = lnrm[1] = 0.0;
			lnrm[2] = 1.0;
			RayObjNormal(pos,ray,objs[u->second],unrm,GUNDEF);
			uh = u->first - xyz[2];
			lh = xyz[2];
		}
		else
		{
			l = u - 1;
			RayObjNormal(pos,ray,objs[u->second],unrm,GUNDEF);
			RayObjNormal(pos,ray,objs[l->second],lnrm,GUNDEF);
			uh = u->first - xyz[2];
			lh = xyz[2] - l->first;
		}
		nrm[0] = (lnrm[0] * uh + unrm[0] * lh) / (uh + lh);
		nrm[1] = (lnrm[1] * uh + unrm[1] * lh) / (uh + lh);
		nrm[2] = (lnrm[2] * uh + unrm[2] * lh) / (uh + lh);
		return true;
	}
#if defined(__GRAPHICS__)
	void Draw()
	{
		for(int i = 0; i < objs.size(); i++)
		{
			double alpha = i/(double)objs.size();
			glColor3f(alpha,1.0-alpha,0.5);
			DrawObj(objs[i],2,0);
		}
		
	}
#endif
};



#endif
