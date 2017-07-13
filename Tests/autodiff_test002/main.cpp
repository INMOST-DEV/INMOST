#include <cstdio>
#include <cmath>
#include "inmost.h"
using namespace INMOST;

int main(int argc,char ** argv)
{
	int test = 0;
	if (argc > 1)  test = atoi(argv[1]);

	const double pi = 3.1415926535897932384626433832795;
	double x = 0.51;
	double y = 0.25;
	double z = 0.15;
	double t = 0.1;
	double dx, dy, dz, dt, dxdx, dydy, dzdz, dtdt, dxdy, dxdz, dydz, dxdt, dydt, dzdt;
	unknown vx(x,0), vy(y,1), vz(z,2), vt(t,3);
	hessian_variable f;
	variable f2;


	if( test == 0 )
	{
		dx = exp(-t)*(sin(2*pi*x*y)*sin(2*pi*z)*sin(2*pi*x*z) + 2 * pi*x*y*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x*sin(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) + 3 * x*x * y*y * z);
		dy = exp(-t)*(2 * pi*x *x * cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * x *x*x * y*z);
		dz = exp(-t)*(2 * pi*x*sin(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x *x * sin(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z) + x *x*x * y*y);
		dt = -exp(-t)*(x*sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + x *x*x * y *y * z);
		dxdx = exp(-t)*(-4 * pi *pi * x*sin(2 * pi*x*y)*z *z * sin(2 * pi*z)*sin(2 * pi*x*z) - 4 * pi *pi * x*y *y * sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*
			y*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*sin(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) + 8 * pi *pi * x*y*cos(2 * pi*x*y)*z*
			sin(2 * pi*z)*cos(2 * pi*x*z) + 6 * x*y *y * z);
		dxdy = exp(-t)*(-4 * pi *pi * x *x * y*sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*x*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi *pi * x *x *
			cos(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) + 6 * x *x * y*z);
		dxdz = exp(-t)*(-4 * pi *pi * x *x * sin(2 * pi*x*y)*z*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*sin(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi *pi * x*y*
			cos(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*x*sin(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z) + 4 * pi *pi * x *x * y*cos(2 * pi*x*y)*sin(2 * pi*z)*
			cos(2 * pi*x*z) + 4 * pi *pi * x*sin(2 * pi*x*y)*z*cos(2 * pi*z)*cos(2 * pi*x*z) + 3 * x *x * y *y);
		dydy = exp(-t)*(2 * x *x*x * z - 4 * pi *pi * x *x*x * sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z));
		dydz = exp(-t)*(4 * pi *pi * x *x * cos(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi *pi * x *x*x * cos(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z) + 2 * x *x*x * y);
		dzdz = exp(-t)*(-4 * pi *pi * x *x*x * sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) - 4 * pi *pi * x*sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 8 * pi *pi * x *x *
			sin(2 * pi*x*y)*cos(2 * pi*z)*cos(2 * pi*x*z));
		dxdt = -exp(-t)*(sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x*y*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x*
			sin(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) + 3 * x *x * y *y * z);
		dydt = -exp(-t)*(2 * pi*x *x * cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * x *x*x * y*z);
		dzdt = -exp(-t)*(2 * pi*x*sin(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x *x * sin(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z) + x *x*x * y *y);
		dtdt = exp(-t)*(x*sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + x *x*x * y *y * z);
		
		f = (vx*vx*vx*vy*vy*vz + vx*sin(2 * pi*vx*vz)*sin(2 * pi*vx*vy)*sin(2 * pi*vz))*exp(-vt);
		f2 = (vx*vx*vx*vy*vy*vz + vx*sin(2 * pi*vx*vz)*sin(2 * pi*vx*vy)*sin(2 * pi*vz))*exp(-vt);
	}
	else if (test == 1)
	{
		dx = t*(sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x*y*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x*sin(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) + 3 * x*x * y*y * z);
		dy = t*(2 * pi*x *x * cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * x *x*x * y*z);
		dz = t*(2 * pi*x*sin(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x *x * sin(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z) + x *x*x * y*y);
		dt = (x*sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + x *x*x * y *y * z);
		dxdx = t*(-4 * pi *pi * x*sin(2 * pi*x*y)*z *z * sin(2 * pi*z)*sin(2 * pi*x*z) - 4 * pi *pi * x*y *y * sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*
			y*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*sin(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) + 8 * pi *pi * x*y*cos(2 * pi*x*y)*z*
			sin(2 * pi*z)*cos(2 * pi*x*z) + 6 * x*y *y * z);
		dxdy = t*(-4 * pi *pi * x *x * y*sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*x*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi *pi * x *x *
			cos(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) + 6 * x *x * y*z);
		dxdz = t*(-4 * pi *pi * x *x * sin(2 * pi*x*y)*z*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*sin(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi *pi * x*y*
			cos(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*x*sin(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z) + 4 * pi *pi * x *x * y*cos(2 * pi*x*y)*sin(2 * pi*z)*
			cos(2 * pi*x*z) + 4 * pi *pi * x*sin(2 * pi*x*y)*z*cos(2 * pi*z)*cos(2 * pi*x*z) + 3 * x *x * y *y);
		dydy = t*(2 * x *x*x * z - 4 * pi *pi * x *x*x * sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z));
		dydz = t*(4 * pi *pi * x *x * cos(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi *pi * x *x*x * cos(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z) + 2 * x *x*x * y);
		dzdz = t*(-4 * pi *pi * x *x*x * sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) - 4 * pi *pi * x*sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 8 * pi *pi * x *x *
			sin(2 * pi*x*y)*cos(2 * pi*z)*cos(2 * pi*x*z));
		dxdt = (sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x*y*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x*
			sin(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) + 3 * x *x * y *y * z);
		dydt = (2 * pi*x *x * cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * x *x*x * y*z);
		dzdt = (2 * pi*x*sin(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x *x * sin(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z) + x *x*x * y *y);
		dtdt = 0;

		f = (vx*vx*vx*vy*vy*vz + vx*sin(2 * pi*vx*vz)*sin(2 * pi*vx*vy)*sin(2 * pi*vz))*vt;
		f2 = (vx*vx*vx*vy*vy*vz + vx*sin(2 * pi*vx*vz)*sin(2 * pi*vx*vy)*sin(2 * pi*vz))*vt;
	}
	else if (test == 2)
	{
		dx = (sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x*y*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x*sin(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) + 3 * x*x * y*y * z);
		dy = (2 * pi*x *x * cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * x *x*x * y*z);
		dz = (2 * pi*x*sin(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x *x * sin(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z) + x *x*x * y*y);
		dt = 0;
		dxdx = (-4 * pi *pi * x*sin(2 * pi*x*y)*z *z * sin(2 * pi*z)*sin(2 * pi*x*z) - 4 * pi *pi * x*y *y * sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*
			y*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*sin(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) + 8 * pi *pi * x*y*cos(2 * pi*x*y)*z*
			sin(2 * pi*z)*cos(2 * pi*x*z) + 6 * x*y *y * z);
		dxdy = (-4 * pi *pi * x *x * y*sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*x*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi *pi * x *x *
			cos(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) + 6 * x *x * y*z);
		dxdz = (-4 * pi *pi * x *x * sin(2 * pi*x*y)*z*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*sin(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi *pi * x*y*
			cos(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*x*sin(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z) + 4 * pi *pi * x *x * y*cos(2 * pi*x*y)*sin(2 * pi*z)*
			cos(2 * pi*x*z) + 4 * pi *pi * x*sin(2 * pi*x*y)*z*cos(2 * pi*z)*cos(2 * pi*x*z) + 3 * x *x * y *y);
		dydy = (2 * x *x*x * z - 4 * pi *pi * x *x*x * sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z));
		dydz = (4 * pi *pi * x *x * cos(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi *pi * x *x*x * cos(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z) + 2 * x *x*x * y);
		dzdz = (-4 * pi *pi * x *x*x * sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) - 4 * pi *pi * x*sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 8 * pi *pi * x *x *
			sin(2 * pi*x*y)*cos(2 * pi*z)*cos(2 * pi*x*z));
		dxdt = 0;
		dydt = 0;
		dzdt = 0;
		dtdt = 0;

		f = (vx*vx*vx*vy*vy*vz + vx*sin(2 * pi*vx*vz)*sin(2 * pi*vx*vy)*sin(2 * pi*vz));
		f2 = (vx*vx*vx*vy*vy*vz + vx*sin(2 * pi*vx*vz)*sin(2 * pi*vx*vy)*sin(2 * pi*vz));
	}
	else if (test == 3)
	{
		dx = (sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x*y*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x*sin(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) );
		dy = (2 * pi*x *x * cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z));
		dz = (2 * pi*x*sin(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*x *x * sin(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z));
		dt = 0;
		dxdx = (-4 * pi *pi * x*sin(2 * pi*x*y)*z *z * sin(2 * pi*z)*sin(2 * pi*x*z) - 4 * pi *pi * x*y *y * sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*
			y*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*sin(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z) + 8 * pi *pi * x*y*cos(2 * pi*x*y)*z*
			sin(2 * pi*z)*cos(2 * pi*x*z));
		dxdy = (-4 * pi *pi * x *x * y*sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*x*cos(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi *pi * x *x *
			cos(2 * pi*x*y)*z*sin(2 * pi*z)*cos(2 * pi*x*z));
		dxdz = (-4 * pi *pi * x *x * sin(2 * pi*x*y)*z*sin(2 * pi*z)*sin(2 * pi*x*z) + 2 * pi*sin(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi *pi * x*y*
			cos(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi*x*sin(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z) + 4 * pi *pi * x *x * y*cos(2 * pi*x*y)*sin(2 * pi*z)*
			cos(2 * pi*x*z) + 4 * pi *pi * x*sin(2 * pi*x*y)*z*cos(2 * pi*z)*cos(2 * pi*x*z));
		dydy = ( - 4 * pi *pi * x *x*x * sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z));
		dydz = (4 * pi *pi * x *x * cos(2 * pi*x*y)*cos(2 * pi*z)*sin(2 * pi*x*z) + 4 * pi *pi * x *x*x * cos(2 * pi*x*y)*sin(2 * pi*z)*cos(2 * pi*x*z));
		dzdz = (-4 * pi *pi * x *x*x * sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) - 4 * pi *pi * x*sin(2 * pi*x*y)*sin(2 * pi*z)*sin(2 * pi*x*z) + 8 * pi *pi * x *x *
			sin(2 * pi*x*y)*cos(2 * pi*z)*cos(2 * pi*x*z));
		dxdt = 0;
		dydt = 0;
		dzdt = 0;
		dtdt = 0;

		f = (vx*sin(2 * pi*vx*vz)*sin(2 * pi*vx*vy)*sin(2 * pi*vz));
		f2 = (vx*sin(2 * pi*vx*vz)*sin(2 * pi*vx*vy)*sin(2 * pi*vz));
	}
	else if (test == 4)
	{
		dx = sin(2 * pi*z);
		dy = 0;
		dz = 2 * pi*x*cos(2 * pi*z);
		dt = 0;
		dxdx = 0;
		dxdy = 0;
		dxdz = 2 * pi*cos(2 * pi*z);
		dydy = 0;
		dydz = 0;
		dzdz = -4 * pi *pi * x*sin(2 * pi*z);
		dxdt = 0;
		dydt = 0;
		dzdt = 0;
		dtdt = 0;

		f = (vx*sin(2 * pi*vz));
		f2 = (vx*sin(2 * pi*vz));
	}
	else if (test == 5)
	{
		dx = sin(2 * pi*x*z) + 2 * pi*x*z*cos(2 * pi*x*z);
		dy = 0;
		dz = 2 * pi*x *x * cos(2 * pi*x*z);
		dt = 0;
		dxdx = 4 * pi*z*cos(2 * pi*x*z) - 4 * pi *pi * x*z *z * sin(2 * pi*x*z);
		dxdy = 0;
		dxdz = 4 * pi*x*cos(2 * pi*x*z) - 4 * pi *pi * x *x * z*sin(2 * pi*x*z);
		dydy = 0;
		dydz = 0;
		dzdz = -4 * pi *pi * x *x*x * sin(2 * pi*x*z);
		dxdt = 0;
		dydt = 0;
		dzdt = 0;
		dtdt = 0;

		f = (vx*sin(2 * pi*vx*vz));
		f2 = (vx*sin(2 * pi*vx*vz));
	}
	else if (test == 6)
	{
		dx = sin(2 * pi*x*y)*sin(2 * pi*z) + 2 * pi*x*y*cos(2 * pi*x*y)*sin(2 * pi*z);
		dy = 2 * pi*x *x * cos(2 * pi*x*y)*sin(2 * pi*z);
		dz = 2 * pi*x*sin(2 * pi*x*y)*cos(2 * pi*z);
		dt = 0;
		dxdx = 4 * pi*y*cos(2 * pi*x*y)*sin(2 * pi*z) - 4 * pi *pi * x*y *y * sin(2 * pi*x*y)*sin(2 * pi*z);
		dxdy = 4 * pi*x*cos(2 * pi*x*y)*sin(2 * pi*z) - 4 * pi *pi * x *x * y*sin(2 * pi*x*y)*sin(2 * pi*z);
		dxdz = 2 * pi*sin(2 * pi*x*y)*cos(2 * pi*z) + 4 * pi *pi * x*y*cos(2 * pi*x*y)*cos(2 * pi*z);
		dydy = -4 * pi *pi * x *x*x * sin(2 * pi*x*y)*sin(2 * pi*z);
		dydz = 4 * pi *pi * x *x * cos(2 * pi*x*y)*cos(2 * pi*z);
		dzdz = -4 * pi *pi * x*sin(2 * pi*x*y)*sin(2 * pi*z);
		dxdt = 0;
		dydt = 0;
		dzdt = 0;
		dtdt = 0;

		f = (vx*sin(2 * pi*vx*vy)*sin(2 * pi*vz));
		f2 = (vx*sin(2 * pi*vx*vy)*sin(2 * pi*vz));
	}
	else if( test == 7 )
	{
		dx = 16*exp(-t)*(x-0.5)*cos(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)));
		dy = 16*exp(-t)*(y-0.5)*cos(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)));
		dz = 16*exp(-t)*cos(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)))*(z-0.5);
		dt = -exp(-t)*(sin(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)))+1);
		dxdx = 16*exp(-t)*cos(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)))-256*exp(-t)*(x-0.5)*(x-0.5)*sin(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)));
		dxdy = -256*exp(-t)*(x-0.5)*(y-0.5)*sin(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)));
		dxdz = -256*exp(-t)*(x-0.5)*sin(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)))*(z-0.5);
		dydy = 16*exp(-t)*cos(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)))-256*exp(-t)*(y-0.5)*(y-0.5)*sin(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)));
		dydz = -256*exp(-t)*(y-0.5)*sin(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)))*(z-0.5);
		dzdz = 16*exp(-t)*cos(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)))-256*exp(-t)*sin(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)))*(z-0.5)*(z-0.5);
		dxdt = -16*exp(-t)*(x-0.5)*cos(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)));
		dydt = -16*exp(-t)*(y-0.5)*cos(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)));
		dzdt = -16*exp(-t)*cos(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)))*(z-0.5);
		dtdt = exp(-t)*(sin(8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5)))+1);
		f = (sin(((vx-0.5)*(vx-0.5) + (vy-0.5)*(vy-0.5) + (vz-0.5)*(vz-0.5))*8)+1)*exp(-vt);
		f2 = (sin(((vx-0.5)*(vx-0.5) + (vy-0.5)*(vy-0.5) + (vz-0.5)*(vz-0.5))*8)+1)*exp(-vt);
	}
	//mixed derivative computed twice: dxdy and dydx
	dxdy *= 2;
	dxdz *= 2;
	dydz *= 2;
	dxdt *= 2;
	dydt *= 2;
	dzdt *= 2;
	
	double vdx = f.GetRow()[0];
	double vdy = f.GetRow()[1];
	double vdz = f.GetRow()[2];
	double vdt = f.GetRow()[3];
	double vdxdx = f.GetHessianRow()[Sparse::HessianRow::make_index(0,0)];
	double vdydy = f.GetHessianRow()[Sparse::HessianRow::make_index(1,1)];
	double vdzdz = f.GetHessianRow()[Sparse::HessianRow::make_index(2,2)];
	double vdtdt = f.GetHessianRow()[Sparse::HessianRow::make_index(3,3)];
	double vdxdy = f.GetHessianRow()[Sparse::HessianRow::make_index(0,1)];
	double vdxdz = f.GetHessianRow()[Sparse::HessianRow::make_index(0,2)];
	double vdydz = f.GetHessianRow()[Sparse::HessianRow::make_index(1,2)];
	double vdxdt = f.GetHessianRow()[Sparse::HessianRow::make_index(0,3)];
	double vdydt = f.GetHessianRow()[Sparse::HessianRow::make_index(1,3)];
	double vdzdt = f.GetHessianRow()[Sparse::HessianRow::make_index(2,3)];
	

	bool error = false;
	std::cout << "GetHessian:" << std::endl;
	std::cout << std::setw(10) << "derivative " << std::setw(10) << "original " << std::setw(10) << "computed" << std::endl;
	std::cout << std::setw(10) << "dx " << std::setw(10) << dx << std::setw(10) << vdx << std::endl;
	std::cout << std::setw(10) << "dy " << std::setw(10) << dy << std::setw(10) << vdy << std::endl;
	std::cout << std::setw(10) << "dz " << std::setw(10) << dz << std::setw(10) << vdz << std::endl;
	std::cout << std::setw(10) << "dt " << std::setw(10) << dt << std::setw(10) << vdt << std::endl;
	std::cout << std::setw(10) << "dxdx " << std::setw(10) << dxdx << std::setw(10) << vdxdx << std::endl;
	std::cout << std::setw(10) << "dxdy " << std::setw(10) << dxdy << std::setw(10) << vdxdy << std::endl;
	std::cout << std::setw(10) << "dxdz " << std::setw(10) << dxdz << std::setw(10) << vdxdz << std::endl;
	std::cout << std::setw(10) << "dxdt " << std::setw(10) << dxdt << std::setw(10) << vdxdt << std::endl;
	std::cout << std::setw(10) << "dydy " << std::setw(10) << dydy << std::setw(10) << vdydy << std::endl;
	std::cout << std::setw(10) << "dydz " << std::setw(10) << dydz << std::setw(10) << vdydz << std::endl;
	std::cout << std::setw(10) << "dydt " << std::setw(10) << dydt << std::setw(10) << vdydt << std::endl;
	std::cout << std::setw(10) << "dzdz " << std::setw(10) << dzdz << std::setw(10) << vdzdz << std::endl;
	std::cout << std::setw(10) << "dzdt " << std::setw(10) << dzdt << std::setw(10) << vdzdt << std::endl;
	std::cout << std::setw(10) << "dtdt " << std::setw(10) << dtdt << std::setw(10) << vdtdt << std::endl;
	if( std::abs(dx-vdx) > 1.0e-9 ) error = true, std::cout << "Error in dx: " << std::abs(dx-vdx) << " original " << dx << " computed " << vdx << std::endl;
	if( std::abs(dy-vdy) > 1.0e-9 ) error = true, std::cout << "Error in dy: " << std::abs(dy-vdy) << " original " << dy << " computed " << vdy << std::endl;
	if( std::abs(dz-vdz) > 1.0e-9 ) error = true, std::cout << "Error in dz: " << std::abs(dz-vdz) << " original " << dz << " computed " << vdz << std::endl;
	if (std::abs(dt - vdt) > 1.0e-9) error = true, std::cout << "Error in dt: " << std::abs(dt - vdt) << " original " << dt << " computed " << vdt << std::endl;
	if( std::abs(dxdx-vdxdx) > 1.0e-9 ) error = true, std::cout << "Error in dxdx: " << std::abs(dxdx-vdxdx) << " original " << dxdx << " computed " << vdxdx << std::endl;
	if( std::abs(dxdy-vdxdy) > 1.0e-9 ) error = true, std::cout << "Error in dxdy: " << std::abs(dxdy-vdxdy) << " original " << dxdy << " computed " << vdxdy << std::endl;
	if( std::abs(dxdz-vdxdz) > 1.0e-9 ) error = true, std::cout << "Error in dxdz: " << std::abs(dxdz-vdxdz) << " original " << dxdz << " computed " << vdxdz << std::endl;
	if (std::abs(dxdt - vdxdt) > 1.0e-9) error = true, std::cout << "Error in dxdt: " << std::abs(dxdt - vdxdt) << " original " << dxdt << " computed " << vdxdt << std::endl;
	if( std::abs(dydy-vdydy) > 1.0e-9 ) error = true, std::cout << "Error in dydy: " << std::abs(dydy-vdydy) << " original " << dydy << " computed " << vdydy << std::endl;
	if( std::abs(dydz-vdydz) > 1.0e-9 ) error = true, std::cout << "Error in dydz: " << std::abs(dydz-vdydz) << " original " << dydz << " computed " << vdydz << std::endl;
	if (std::abs(dydt - vdydt) > 1.0e-9) error = true, std::cout << "Error in dydt: " << std::abs(dydt - vdydt) << " original " << dydt << " computed " << vdydt << std::endl;
	if( std::abs(dzdz-vdzdz) > 1.0e-9 ) error = true, std::cout << "Error in dzdz: " << std::abs(dzdz-vdzdz) << " original " << dzdz << " computed " << vdzdz << std::endl;
	if (std::abs(dzdt - vdzdt) > 1.0e-9) error = true, std::cout << "Error in dzdt: " << std::abs(dzdt - vdzdt) << " original " << dzdt << " computed " << vdzdt << std::endl;
	if (std::abs(dtdt - vdtdt) > 1.0e-9) error = true, std::cout << "Error in dtdt: " << std::abs(dtdt - vdtdt) << " original " << dtdt << " computed " << vdtdt << std::endl;
	//if( error ) return -1;
	vdx = f2.GetRow()[0];
	vdy = f2.GetRow()[1];
	vdz = f2.GetRow()[2];
	vdt = f2.GetRow()[3];
	std::cout << "GetJacobian:" << std::endl;
	std::cout << std::setw(10) << "derivative " << std::setw(10) << "original " << std::setw(10) << "computed" << std::endl;
	std::cout << std::setw(10) << "dx " << std::setw(10) << dx << std::setw(10) << vdx << std::endl;
	std::cout << std::setw(10) << "dy " << std::setw(10) << dy << std::setw(10) << vdy << std::endl;
	std::cout << std::setw(10) << "dz " << std::setw(10) << dz << std::setw(10) << vdz << std::endl;
	std::cout << std::setw(10) << "dt " << std::setw(10) << dt << std::setw(10) << vdt << std::endl;
	if (std::abs(dx - vdx) > 1.0e-9) error = true, std::cout << "Error in dx: " << std::abs(dx - vdx) << " original " << dx << " computed " << vdx << std::endl;
	if (std::abs(dy - vdy) > 1.0e-9) error = true, std::cout << "Error in dy: " << std::abs(dy - vdy) << " original " << dy << " computed " << vdy << std::endl;
	if (std::abs(dz - vdz) > 1.0e-9) error = true, std::cout << "Error in dz: " << std::abs(dz - vdz) << " original " << dz << " computed " << vdz << std::endl;
	if (std::abs(dt - vdt) > 1.0e-9) error = true, std::cout << "Error in dt: " << std::abs(dt - vdt) << " original " << dt << " computed " << vdt << std::endl;
	if (error) return -1;

	return 0;
}
