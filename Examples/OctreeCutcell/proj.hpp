#ifndef _PROJECTION_HPP
#define _PROJECTION_HPP
#define PEPS 2.5e-10
#include "obj.h"
#include <vector>
#include <string>
#include <math.h>
#include "my_glut.h"
#include <stdio.h>
#include <string.h>

#define WAITNL {int ret; char c; printf("press Enter\n"); ret = scanf("%c",&c);}

class Projection
{
	std::vector<int> Gbound;
	double * Gcoord;
	double * Gcoef;
	int Gstep;
public:
	Projection()
	{
		Gcoef = NULL;
		Gcoord = NULL;
		Gbound.clear();
		Gstep = 1;
	}
	~Projection()
	{
		if( Gcoef != NULL ) delete [] Gcoef;
		if( Gcoord != NULL ) delete [] Gcoord;
		Gbound.clear();
	}
	Projection(const Projection & other)
	{
		Gbound = other.Gbound;
		Gcoef = new double[Gbound.size()];
		memcpy(Gcoef,other.Gcoef,sizeof(double)*Gbound.size());
		Gcoord = new double[Gbound.size()+1];
		memcpy(Gcoord,other.Gcoord,sizeof(double)*(Gbound.size()+1));
		Gstep = other.Gstep;
	}
	
	Projection & operator =(Projection const & other)
	{
		Gbound = other.Gbound;
		Gcoef = new double[Gbound.size()];
		memcpy(Gcoef,other.Gcoef,sizeof(double)*Gbound.size());
		Gcoord = new double[Gbound.size()+1];
		memcpy(Gcoord,other.Gcoord,sizeof(double)*(Gbound.size()+1));
		Gstep = other.Gstep;
		return *this;
	}
	void Project(double in[3])
	{
		double pos[3];
		double floor = 0.0;
		double ray[3];
		ray[0] = 0;
		ray[1] = 0;
		ray[2] = 1;
		/*
		if((in[0] < -PEPS || in[0] > 1+PEPS) ||
		   (in[1] < -PEPS || in[1] > 1+PEPS) ||
		   (in[2] < -PEPS || in[2] > 1+PEPS) )
		{
			return;
		}
		*/
		//printf("in %g %g %g\n",in[0], in[1], in[2]);
		pos[2] = in[2];
		if( pos[2] < 0.0 ) pos[2] = 0;
		if( pos[2] > 1.0 ) pos[2] = 1.0;
		int k = Gbound.size();
		for(int i = 0; i < Gbound.size(); i++)
			if( pos[2] <= Gcoord[i+1]+PEPS )
			{
				k = i;
				break;
			}
		if( k == Gbound.size() ) throw 0;
		pos[0] = in[0];
		pos[1] = in[1];
		if( pos[0] < 0.0+PEPS ) pos[0] = 0.0+PEPS;
		if( pos[0] > 1.0-PEPS ) pos[0] = 1.0-PEPS;
		if( pos[1] < 0.0+PEPS ) pos[1] = 0.0+PEPS;
		if( pos[1] > 1.0-PEPS ) pos[1] = 1.0-PEPS;
		pos[2] = 0;
		double dist = RayObjIntersection(pos,ray,Gbound[k],1.0e+25);
		if( k > 0 )
			floor = RayObjIntersection(pos,ray,Gbound[k-1],1.0e+25);
		double alpha = (in[2]-Gcoord[k])/Gcoef[k];
		in[2] = (floor+(dist-floor)*alpha);
		//printf("k %d fl %g dist %g floor %g alpha %g Gcoord[k] %g Gcoef[k] %g\n",
		//	k,fl,dist,floor,alpha,Gcoord[k],Gcoef[k]);
		//printf("out %g %g %g\n",in[0], in[1], in[2]);
		//WAITNL;
	}

	int Layer(double in[3])
	{
		double distmin = 1.0e+25;
		int k = Gbound.size();
		double pos[3];
		pos[0] = in[0];
		pos[1] = in[1];
		pos[2] = in[2];
		double ray[3];
		ray[0] = 0;
		ray[1] = 0;
		ray[2] = 1;
		for(int i = 0; i < Gbound.size(); i++)
		{
			double dist = RayObjIntersection(pos,ray,Gbound[i],distmin);
			if( dist < distmin )
			{
				distmin = dist;
				k = i;
			}
		}
		return k;
	}
	int Size() {return Gbound.size();}
	void ReadLayers(std::vector<std::string> & input)
	{
		int i;
		double h = 1.0/(double)Gstep;
		double * Middle = new double[input.size()],norm = 0;
		for(i = 0; i < input.size(); i++)
		{
			Gbound.push_back(ReadObj((char *)input[i].c_str()));
			Middle[i] = MiddleZ(Gbound[i]);
			//printf("middle %d = %e\n",i,Middle[i]);
		}
		Gcoef = new double[Gbound.size()];
		Gcoord = new double[Gbound.size()+1];
		Gcoef[0] = Middle[0];
		for(i = 1; i < Gbound.size(); i++) Gcoef[i] = (Middle[i]-Middle[i-1]);
		for(i = 0; i < Gbound.size(); i++) norm += Gcoef[i];
		for(i = 0; i < Gbound.size(); i++) Gcoef[i] /= norm;
		Gcoord[0] = 0.0;
		for(i = 1; i < Gbound.size()+1; i++) Gcoord[i] = Gcoord[i-1]+Gcoef[i-1];
		if( Gstep > 1 )
		{
			for(i = 1; i < Gbound.size()+1; i++)
			{
				Gcoord[i] = ceil(Gcoord[i]/h-0.5)*h;
				if( i == 1 && Gcoord[i] < PEPS ) Gcoord[i] = h;
				if( Gcoord[i] <= Gcoord[i-1] ) Gcoord[i] = Gcoord[i-1] + h;
				//printf("Gcoord %d = %e\n", i, Gcoord[i]);
			}
			for(i = 0; i < Gbound.size(); i++)
			{
				Gcoef[i] = Gcoord[i+1]-Gcoord[i];
				//printf("Gcoef %d = %e\n",i, Gcoef[i]);
			}
		}
		delete [] Middle;
	}
	void SetGridStep(int n)
	{	
		Gstep = n;
	}
#if defined(__GRAPHICS__)
	void Draw()
	{
		for(int i = 0; i < Gbound.size(); i++)
		{
			glColor3f(0,0,0);
			DrawObj(Gbound[i],2,0);
		}
	}
#endif
};

#endif //_PROJECTION_HPP
