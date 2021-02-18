#include "inmost.h"
#include <cstdio>
#include <cmath>
using namespace INMOST;

int main(int argc,char ** argv)
{
	int test = 0;
	if (argc > 1)  test = atoi(argv[1]);

	INMOST_DATA_REAL_TYPE _x = 0.5;
	INMOST_DATA_REAL_TYPE _y = 0.25;
	INMOST_DATA_REAL_TYPE _dx, _dy, _dxdx, _dydy, _dxdy;
	unknown x(_x,0), y(_y,1);
	hessian_variable f;
	variable f2;
#if defined(USE_FP64)
	INMOST_DATA_REAL_TYPE tol = 1.0e-9;
#else
	INMOST_DATA_REAL_TYPE tol = 1.0e-6;
#endif

	if( test == 0 ) //check derivative and hessian of sin(x*x+y*y)
	{
		// (sin(x*x+y*y))'' = (2cos(x*x+y*y)(xdx+ydy))' =
		// 2cos(x*x+y*y)(dxdx+dydy) - 4sin(x*x+y*y)(x*x dxdx+y*y dydy + 2x*y dxdy)(xdx+ydy)
		//
		// Expected Jacobian
		// 0: 2cos(x*x+y*y)x dx
		// 1: 2cos(x*x+y*y)y dy
		//
		// Expected Hessian
		// 0,0: 2*cos(x*x+y*y)-4*sin(x*x+y*y) x*x dxdx
		// 0,1: -8*sin(x*x+y*y) x*y dx dy
		// 1,1: 2*cos(x*x+y*y)-4*sin(x*x+y*y) y*y dydy
		_dx = 2*cos(_x*_x+_y*_y)*_x;
		_dy = 2*cos(_x*_x+_y*_y)*_y;
		_dxdx = 2*cos(_x*_x+_y*_y)-4*sin(_x*_x+_y*_y)*_x*_x;
		_dydy = 2*cos(_x*_x+_y*_y)-4*sin(_x*_x+_y*_y)*_y*_y;
		_dxdy = -8*sin(_x*_x+_y*_y)*_x*_y;
		f = sin(x*x+y*y);
		f2 = sin(x*x + y*y);
	}
	else if( test == 1 )
	{
		// (sin(x*x+y*y+x*y))'' = (cos(x*x+y*y+x*y)(2xdx+2ydy+xdy+ydx))' =
		// 2cos(x*x+y*y+x*y)(dxdx+dydy+dxdy) - sin(x*x+y*y+x*y)((2x+y)dx+(2y+x)dy)*((2x+y)dx+(2y+x)dy) =
		// 2cos(x*x+y*y+x*y)(dxdx+dydy+dxdy) - sin(x*x+y*y+x*y)((2x+y)(2x+y)dxdx + 2(2x+y)(2y+x)dxdy + (2y+x)(2y+x)dydy)
		//
		// Expected Jacobian
		// 0: cos(x*x+y*y+x*y)(2x+y) dx
		// 1: cos(x*x+y*y+x*y)(2y+x) dy
		//
		// Expected Hessian
		// 0,0: 2cos(x*x+y*y+x*y) -  sin(x*x+y*y+x*y) (2x+y)(2x+y) dxdx
		// 0,1: 2cos(x*x+y*y+x*y) - 2sin(x*x+y*y+x*y) (2x+y)(2y+x) dxdy
		// 1,1: 2cos(x*x+y*y+x*y) -  sin(x*x+y*y+x*y) (x+2y)(x+2y) dydy
		_dx = cos(_x*_x+_y*_y+_x*_y)*(2*_x+_y);
		_dy = cos(_x*_x+_y*_y+_x*_y)*(_x+2*_y);
		_dxdx = 2*cos(_x*_x+_y*_y+_x*_y)-sin(_x*_x+_y*_y+_x*_y)*(2*_x+_y)*(2*_x+_y);
		_dydy = 2*cos(_x*_x+_y*_y+_x*_y)-sin(_x*_x+_y*_y+_x*_y)*(_x+2*_y)*(_x+2*_y);
		_dxdy = 2*cos(_x*_x+_y*_y+_x*_y)-2*sin(_x*_x+_y*_y+_x*_y)*(2*_x+_y)*(_x+2*_y);
		f = sin(x*x+y*y+x*y);
		f2 = sin(x*x + y*y + x*y);
	}
	else if( test == 2 )
	{
		// (cos(x*x+y*y+x*y))'' = (-sin(x*x+y*y+x*y)(2xdx+2ydy+xdy+ydx))' =
		// -2sin(x*x+y*y+x*y)(dxdx+dydy+dxdy) - cos(x*x+y*y+x*y)((2x+y)dx+(2y+x)dy)*((2x+y)dx+(2y+x)dy) =
		// -2sin(x*x+y*y+x*y)(dxdx+dydy+dxdy) - cos(x*x+y*y+x*y)((2x+y)(2x+y)dxdx + 2(2x+y)(2y+x)dxdy + (2y+x)(2y+x)dydy)
		//
		// Expected Jacobian
		// 0: -sin(x*x+y*y+x*y)(2x+y) dx
		// 1: -sin(x*x+y*y+x*y)(2y+x) dy
		//
		// Expected Hessian
		// 0,0: -2sin(x*x+y*y+x*y) -  cos(x*x+y*y+x*y) (2x+y)(2x+y) dxdx
		// 0,1: -2sin(x*x+y*y+x*y) - 2cos(x*x+y*y+x*y) (2x+y)(2y+x) dxdy
		// 1,1: -2sin(x*x+y*y+x*y) -  cos(x*x+y*y+x*y) (x+2y)(x+2y) dydy
		_dx = -sin(_x*_x+_y*_y+_x*_y)*(2*_x+_y);
		_dy = -sin(_x*_x+_y*_y+_x*_y)*(_x+2*_y);
		_dxdx = -2*sin(_x*_x+_y*_y+_x*_y)-cos(_x*_x+_y*_y+_x*_y)*(2*_x+_y)*(2*_x+_y);
		_dydy = -2*sin(_x*_x+_y*_y+_x*_y)-cos(_x*_x+_y*_y+_x*_y)*(_x+2*_y)*(_x+2*_y);
		_dxdy = -2*sin(_x*_x+_y*_y+_x*_y)-2*cos(_x*_x+_y*_y+_x*_y)*(2*_x+_y)*(_x+2*_y);
		f = cos(x*x+y*y+x*y);
		f2 = cos(x*x + y*y + x*y);
	}
	else if( test == 3 )
	{
		//f = sin(x)*y
		//Expected Jacobian
		// 0: y cos(x)
		// 1: sin(x)
		//Expected Hessian
		// 0,0: -y*sin(x)
		// 0,1: 2*cos(x)
		// 1,1: 0
		_dx = _y*cos(_x);
		_dy = sin(_x);
		_dxdx = -_y*sin(_x);
		_dxdy = 2*cos(_x);
		_dydy = 0;
		f = sin(x)*y;
		f2 = sin(x)*y;
	}
	else if( test == 4 )
	{
		//f = cos(x)*y
		//Expected Jacobian
		// 0: -y sin(x)
		// 1: cos(x)
		//Expected Hessian
		// 0,0: -y*cos(x)
		// 0,1: -2*sin(x)
		// 1,1: 0
		_dx = -_y*sin(_x);
		_dy = cos(_x);
		_dxdx = -_y*cos(_x);
		_dxdy = -2*sin(_x);
		_dydy = 0;
		f = cos(x)*y;
		f2 = cos(x)*y;
	}
	else if( test == 5 )
	{
		// f = (sqrt(x*x+y*y+x*y))'
		// Expected Jacobian
		// 0: (2x+y)/(2sqrt(x*x+x*y+y*y)
		// 1: (x+2y)/(2sqrt(x*x+x*y+y*y)
		// Expected Hessian
		// 0,0: 3*y*y/(4(x*x+x*y+y*y)^(3/2))
		// 0,1: -6*x*y/(4*(x*x+x*y+y*y)^(3/2))
		// 1,1: 3*x*x/(4(x*x+x*y+y*y)^(3/2))
		_dx = (2*_x+_y)/(2*sqrt(_x*_x+_x*_y+_y*_y));
		_dy = (2*_y+_x)/(2*sqrt(_x*_x+_x*_y+_y*_y));
		_dxdx = 3*_y*_y/(4*pow(_x*_x+_x*_y+_y*_y,1.5));
		_dxdy = -6*_x*_y/(4*pow(_x*_x+_x*_y+_y*_y,1.5));
		_dydy = 3*_x*_x/(4*pow(_x*_x+_x*_y+_y*_y,1.5));
		f = sqrt(x*x+y*y+x*y);
		f2 = sqrt(x*x + y*y + x*y);
	}
	else if( test == 6 )
	{
		_dx = 4*_x + 4*_y;
		_dy = 6*_y + 4*_x;
		_dxdx = 4;
		_dxdy = 8;
		_dydy = 6;
		f = 2*x*x+3*y*y+4*x*y;
		f2 = 2 * x*x + 3 * y*y + 4 * x*y;
	}
	else if( test == 7 )
	{
		_dx = 4 *(_x - 0.5)*cos(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5)));
		_dy = 4 *(_y - 0.5)*cos(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5)));
		_dxdx = 4*cos(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))) - 16*(_x - 0.5)*(_x - 0.5)*sin(2.0*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5)));
		_dxdy = -32*(_x-0.5)*(_y-0.5)*sin(2*((_x-0.5)*(_x-0.5)+(_y-0.5)*(_y-0.5)));
		_dydy = 4*cos(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))) - 16*(_y - 0.5)*(_y - 0.5)*sin(2.0*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5)));
		f = sin(2*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
		f2 = sin(2 * ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5)));
	}
	else if( test == 8 )
	{
		_dx = -(4 *(_x - 0.5)*cos(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))));
		_dy = -(4 *(_y - 0.5)*cos(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))));
		_dxdx = -(4*cos(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))) - 16*(_x - 0.5)*(_x - 0.5)*sin(2.0*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))));
		_dxdy = -(-32*(_x-0.5)*(_y-0.5)*sin(2*((_x-0.5)*(_x-0.5)+(_y-0.5)*(_y-0.5))));
		_dydy = -(4*cos(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))) - 16*(_y - 0.5)*(_y - 0.5)*sin(2.0*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))));
		f = 1.0-sin(2*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
		f2 = 1.0 - sin(2 * ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5)));
	}
	else if( test == 9 )
	{
		_dx = (4 *(_x - 0.5)*sin(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))));
		_dy = (4 *(_y - 0.5)*sin(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))));
		_dxdx = 4*sin(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))) + 16*(_x - 0.5)*(_x - 0.5)*cos(2.0*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5)));
		_dxdy = 32*(_x-0.5)*(_y-0.5)*cos(2*((_x-0.5)*(_x-0.5)+(_y-0.5)*(_y-0.5)));
		_dydy = 4*sin(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))) + 16*(_y - 0.5)*(_y - 0.5)*cos(2.0*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5)));
		f = 1.0-cos(2*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
		f2 = 1.0 - cos(2 * ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5)));
	}
	else if( test == 10 )
	{
		_dx = sin(2*((_x-0.5)*(_x-0.5)+(_y-0.5)*(_y-0.5)))+4 *(_x - 0.5)*_x*cos(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5)));
		_dy = 4 *_x*(_y - 0.5)*cos(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5)));
		_dxdx = 8*(_x-0.5)*cos(2.0*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))) + 4*_x*cos(2*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5))) - 16*_x*(_x - 0.5)*(_x - 0.5)*sin(2.0*((_x - 0.5)*(_x - 0.5) + (_y - 0.5)*(_y - 0.5)));
		_dxdy = 8*(_y-0.5)*cos(2*((_x-0.5)*(_x-0.5)+(_y-0.5)*(_y-0.5)))-4*(_x-0.5)*_x*sin(2*((_x-0.5)*(_x-0.5)+(_y-0.5)*(_y-0.5)));
		_dydy = 4*_x*cos(2*((_x-0.5)*(_x-0.5)+(_y-0.5)*(_y-0.5))) - 16*_x*(_y-0.5)*(_y-0.5)*sin(2*((_x-0.5)*(_x-0.5)+(_y-0.5)*(_y-0.5)));
		f = sin(2*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)))*x;
		f2 = sin(2 * ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5)))*x;
	}
	else if( test == 11 )
	{
		//f = sin(x)*y
		//Expected Jacobian
		// 0: y cos(x)
		// 1: sin(x)
		//Expected Hessian
		// 0,0: -y*sin(x)
		// 0,1: 2*cos(x)
		// 1,1: 0
		_dx = sin(_x*_x+_y) + 2*_x*_x*cos(_x*_x+_y);
		_dy = _x*cos(_x*_x+_y);
		_dxdx = 6*_x*cos(_x*_x+_y) - 4*_x*_x*_x*sin(_x*_x+_y);
		_dxdy = 2*cos(_x*_x+_y) - 4*_x*_x*sin(_x*_x+_y);
		_dydy = -_x*sin(_x*_x+_y);
		f = sin(x*x+y)*x;
		f2 = sin(x*x + y)*x;
	}
	else if( test == 12 )
	{
		//f = sin(x)*y
		//Expected Jacobian
		// 0: y cos(x)
		// 1: sin(x)
		//Expected Hessian
		// 0,0: -y*sin(x)
		// 0,1: 2*cos(x)
		// 1,1: 0
		_dx = sin(_x*_x+_y*_y) + 2*_x*_x*cos(_x*_x+_y*_y);
		_dy = 2*_x*_y*cos(_x*_x+_y*_y);
		_dxdx = 6*_x*cos(_x*_x+_y*_y) - 4*_x*_x*_x*sin(_x*_x+_y*_y);
		_dxdy = 4*_y*(cos(_x*_x+_y*_y) - 2*_x*_x*sin(_x*_x+_y*_y));
		_dydy = 2*_x*(cos(_x*_x+_y*_y)-2*_y*_y*sin(_x*_x+_y*_y));
		f = sin(x*x+y*y)*x;
		f2 = sin(x*x + y*y)*x;
	}
	INMOST_DATA_REAL_TYPE dx = f.GetRow()[0];
	INMOST_DATA_REAL_TYPE dy = f.GetRow()[1];
	INMOST_DATA_REAL_TYPE dxdx = f.GetHessianRow()[Sparse::HessianRow::make_index(0,0)];
	INMOST_DATA_REAL_TYPE dxdy = f.GetHessianRow()[Sparse::HessianRow::make_index(0,1)];
	INMOST_DATA_REAL_TYPE dydy = f.GetHessianRow()[Sparse::HessianRow::make_index(1,1)];

	bool error = false;
	std::cout << "For GetHessian:" << std::endl;
	std::cout << std::setw(10) << "derivative " << std::setw(10) << "original " << std::setw(10) << "computed" << std::endl;
	std::cout << std::setw(10) << "dx " << std::setw(10) << _dx << std::setw(10) << dx << std::endl;
	std::cout << std::setw(10) << "dy " << std::setw(10) << _dy << std::setw(10) << dy << std::endl;
	std::cout << std::setw(10) << "dxdx " << std::setw(10) << _dxdx << std::setw(10) << dxdx << std::endl;
	std::cout << std::setw(10) << "dxdy " << std::setw(10) << _dxdy << std::setw(10) << dxdy << std::endl;
	std::cout << std::setw(10) << "dydy " << std::setw(10) << _dydy << std::setw(10) << dydy << std::endl;
	if( std::abs(dx-_dx) > tol ) error = true, std::cout << "Error in dx: " << std::abs(dx-_dx) << " original " << _dx << " computed " << dx << std::endl;
	if( std::abs(dy-_dy) > tol ) error = true, std::cout << "Error in dy: " << std::abs(dy-_dy) << " original " << _dy << " computed " << dy << std::endl;
	if( std::abs(dxdx-_dxdx) > tol ) error = true, std::cout << "Error in dxdx: " << std::abs(dxdx-_dxdx) << " original " << _dxdx << " computed " << dxdx << std::endl;
	if( std::abs(dxdy-_dxdy) > tol ) error = true, std::cout << "Error in dxdy: " << std::abs(dxdy-_dxdy) << " original " << _dxdy << " computed " << dxdy << std::endl;
	if( std::abs(dydy-_dydy) > tol ) error = true, std::cout << "Error in dydy: " << std::abs(dydy-_dydy) << " original " << _dydy << " computed " << dydy << std::endl;

	

	dx = f2.GetRow()[0];
	dy = f2.GetRow()[1];
	std::cout << "For GetJacobian:" << std::endl;
	std::cout << std::setw(10) << "derivative " << std::setw(10) << "original " << std::setw(10) << "computed" << std::endl;
	std::cout << std::setw(10) << "dx " << std::setw(10) << _dx << std::setw(10) << dx << std::endl;
	std::cout << std::setw(10) << "dy " << std::setw(10) << _dy << std::setw(10) << dy << std::endl;

	if( error ) return -1;
	return 0;
}
