#include <cstdio>
#include <cmath>
#include "inmost.h"
using namespace INMOST;

int main(int argc,char ** argv)
{
	int test = 0;
	if (argc > 1)  test = atoi(argv[1]);


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
		double _x = 0.5;
		double _y = 0.25;
		double _dx = 2*cos(_x*_x+_y*_y)*_x;
		double _dy = 2*cos(_x*_x+_y*_y)*_y;
		double _dxdx = 2*cos(_x*_x+_y*_y)-4*sin(_x*_x+_y*_y)*_x*_x;
		double _dydy = 2*cos(_x*_x+_y*_y)-4*sin(_x*_x+_y*_y)*_y*_y;
		double _dxdy = -8*sin(_x*_x+_y*_y)*_x*_y;
		unknown x(_x,0), y(_y,1);
		hessian_variable f = sin(x*x+y*y);
		double dx = f.GetRow()[0];
		double dy = f.GetRow()[1];
		double dxdx = f.GetHessianRow()[Sparse::HessianRow::make_index(0,0)];
		double dxdy = f.GetHessianRow()[Sparse::HessianRow::make_index(0,1)];
		double dydy = f.GetHessianRow()[Sparse::HessianRow::make_index(1,1)];

		bool error = false;
		std::cout << std::setw(10) << "derivative " << std::setw(10) << "original " << std::setw(10) << "computed" << std::endl;
		std::cout << std::setw(10) << "dx " << std::setw(10) << _dx << std::setw(10) << dx << std::endl;
		std::cout << std::setw(10) << "dy " << std::setw(10) << _dy << std::setw(10) << dy << std::endl;
		std::cout << std::setw(10) << "dxdx " << std::setw(10) << _dxdx << std::setw(10) << dxdx << std::endl;
		std::cout << std::setw(10) << "dxdy " << std::setw(10) << _dxdy << std::setw(10) << dxdy << std::endl;
		std::cout << std::setw(10) << "dydy " << std::setw(10) << _dydy << std::setw(10) << dydy << std::endl;
		if( std::abs(dx-_dx) > 1.0e-9 ) error = true, std::cout << "Error in dx: " << std::abs(dx-_dx) << " original " << _dx << " computed " << dx << std::endl;
		if( std::abs(dy-_dy) > 1.0e-9 ) error = true, std::cout << "Error in dy: " << std::abs(dy-_dy) << " original " << _dy << " computed " << dy << std::endl;
		if( std::abs(dxdx-_dxdx) > 1.0e-9 ) error = true, std::cout << "Error in dxdx: " << std::abs(dxdx-_dxdx) << " original " << _dxdx << " computed " << dxdx << std::endl;
		if( std::abs(dxdy-_dxdy) > 1.0e-9 ) error = true, std::cout << "Error in dxdy: " << std::abs(dxdy-_dxdy) << " original " << _dxdy << " computed " << dxdy << std::endl;
		if( std::abs(dydy-_dydy) > 1.0e-9 ) error = true, std::cout << "Error in dydy: " << std::abs(dydy-_dydy) << " original " << _dydy << " computed " << dydy << std::endl;

		if( error ) return -1;
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
		double _x = 0.5;
		double _y = 0.25;
		double _dx = cos(_x*_x+_y*_y+_x*_y)*(2*_x+_y);
		double _dy = cos(_x*_x+_y*_y+_x*_y)*(_x+2*_y);
		double _dxdx = 2*cos(_x*_x+_y*_y+_x*_y)-sin(_x*_x+_y*_y+_x*_y)*(2*_x+_y)*(2*_x+_y);
		double _dydy = 2*cos(_x*_x+_y*_y+_x*_y)-sin(_x*_x+_y*_y+_x*_y)*(_x+2*_y)*(_x+2*_y);
		double _dxdy = 2*cos(_x*_x+_y*_y+_x*_y)-2*sin(_x*_x+_y*_y+_x*_y)*(2*_x+_y)*(_x+2*_y);
		unknown x(_x,0), y(_y,1);
		hessian_variable f = sin(x*x+y*y+x*y);
		double dx = f.GetRow()[0];
		double dy = f.GetRow()[1];
		double dxdx = f.GetHessianRow()[Sparse::HessianRow::make_index(0,0)];
		double dxdy = f.GetHessianRow()[Sparse::HessianRow::make_index(0,1)];
		double dydy = f.GetHessianRow()[Sparse::HessianRow::make_index(1,1)];

		bool error = false;
		std::cout << std::setw(10) << "derivative " << std::setw(10) << "original " << std::setw(10) << "computed" << std::endl;
		std::cout << std::setw(10) << "dx " << std::setw(10) << _dx << std::setw(10) << dx << std::endl;
		std::cout << std::setw(10) << "dy " << std::setw(10) << _dy << std::setw(10) << dy << std::endl;
		std::cout << std::setw(10) << "dxdx " << std::setw(10) << _dxdx << std::setw(10) << dxdx << std::endl;
		std::cout << std::setw(10) << "dxdy " << std::setw(10) << _dxdy << std::setw(10) << dxdy << std::endl;
		std::cout << std::setw(10) << "dydy " << std::setw(10) << _dydy << std::setw(10) << dydy << std::endl;
		if( std::abs(dx-_dx) > 1.0e-9 ) error = true, std::cout << "Error in dx: " << std::abs(dx-_dx) << " original " << _dx << " computed " << dx << std::endl;
		if( std::abs(dy-_dy) > 1.0e-9 ) error = true, std::cout << "Error in dy: " << std::abs(dy-_dy) << " original " << _dy << " computed " << dy << std::endl;
		if( std::abs(dxdx-_dxdx) > 1.0e-9 ) error = true, std::cout << "Error in dxdx: " << std::abs(dxdx-_dxdx) << " original " << _dxdx << " computed " << dxdx << std::endl;
		if( std::abs(dxdy-_dxdy) > 1.0e-9 ) error = true, std::cout << "Error in dxdy: " << std::abs(dxdy-_dxdy) << " original " << _dxdy << " computed " << dxdy << std::endl;
		if( std::abs(dydy-_dydy) > 1.0e-9 ) error = true, std::cout << "Error in dydy: " << std::abs(dydy-_dydy) << " original " << _dydy << " computed " << dydy << std::endl;

		if( error ) return -1;
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
		double _x = 0.5;
		double _y = 0.25;
		double _dx = -sin(_x*_x+_y*_y+_x*_y)*(2*_x+_y);
		double _dy = -sin(_x*_x+_y*_y+_x*_y)*(_x+2*_y);
		double _dxdx = -2*sin(_x*_x+_y*_y+_x*_y)-cos(_x*_x+_y*_y+_x*_y)*(2*_x+_y)*(2*_x+_y);
		double _dydy = -2*sin(_x*_x+_y*_y+_x*_y)-cos(_x*_x+_y*_y+_x*_y)*(_x+2*_y)*(_x+2*_y);
		double _dxdy = -2*sin(_x*_x+_y*_y+_x*_y)-2*cos(_x*_x+_y*_y+_x*_y)*(2*_x+_y)*(_x+2*_y);
		unknown x(_x,0), y(_y,1);
		hessian_variable f = cos(x*x+y*y+x*y);
		double dx = f.GetRow()[0];
		double dy = f.GetRow()[1];
		double dxdx = f.GetHessianRow()[Sparse::HessianRow::make_index(0,0)];
		double dxdy = f.GetHessianRow()[Sparse::HessianRow::make_index(0,1)];
		double dydy = f.GetHessianRow()[Sparse::HessianRow::make_index(1,1)];

		bool error = false;
		std::cout << std::setw(10) << "derivative " << std::setw(10) << "original " << std::setw(10) << "computed" << std::endl;
		std::cout << std::setw(10) << "dx " << std::setw(10) << _dx << std::setw(10) << dx << std::endl;
		std::cout << std::setw(10) << "dy " << std::setw(10) << _dy << std::setw(10) << dy << std::endl;
		std::cout << std::setw(10) << "dxdx " << std::setw(10) << _dxdx << std::setw(10) << dxdx << std::endl;
		std::cout << std::setw(10) << "dxdy " << std::setw(10) << _dxdy << std::setw(10) << dxdy << std::endl;
		std::cout << std::setw(10) << "dydy " << std::setw(10) << _dydy << std::setw(10) << dydy << std::endl;
		if( std::abs(dx-_dx) > 1.0e-9 ) error = true, std::cout << "Error in dx: " << std::abs(dx-_dx) << " original " << _dx << " computed " << dx << std::endl;
		if( std::abs(dy-_dy) > 1.0e-9 ) error = true, std::cout << "Error in dy: " << std::abs(dy-_dy) << " original " << _dy << " computed " << dy << std::endl;
		if( std::abs(dxdx-_dxdx) > 1.0e-9 ) error = true, std::cout << "Error in dxdx: " << std::abs(dxdx-_dxdx) << " original " << _dxdx << " computed " << dxdx << std::endl;
		if( std::abs(dxdy-_dxdy) > 1.0e-9 ) error = true, std::cout << "Error in dxdy: " << std::abs(dxdy-_dxdy) << " original " << _dxdy << " computed " << dxdy << std::endl;
		if( std::abs(dydy-_dydy) > 1.0e-9 ) error = true, std::cout << "Error in dydy: " << std::abs(dydy-_dydy) << " original " << _dydy << " computed " << dydy << std::endl;

		if( error ) return -1;
	}
	else if( test == 3 )
	{
		// (sqrt(x*x+y*y+x*y))'
	}

	return 0;
}
