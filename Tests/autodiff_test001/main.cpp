#include "inmost.h"
#include <cstdio>
#include <cmath>
using namespace INMOST;

int main(int argc,char ** argv)
{
	int test = 0;
	if (argc > 1)  test = atoi(argv[1]);

	double x = 0.5;
	double y = 0.25;
	double z = 0.15;
	double dx, dy, dz, dxdx, dydy, dzdz, dxdy, dxdz, dydz;
	unknown vx(x,0), vy(y,1), vz(z,2);
	hessian_variable f;
	variable f2;


	if( test == 0 ) //check derivative and hessian of sin(x*x+y*y)
	{
		double arg = 8*((z-0.5)*(z-0.5)+(y-0.5)*(y-0.5)+(x-0.5)*(x-0.5));
		dx = 2*sin(arg)+16*(x-0.5)*(2*x-1)*cos(arg);
		dy = 16*(2*x-1)*(y-0.5)*cos(arg);
		dz = 16*(2*x-1)*cos(arg)*(z-0.5);
		dxdx = -256*(x-0.5)*(x-0.5)*(2*x-1)*sin(arg)+16*(2*x-1)*cos(arg)+64*(x-0.5)*cos(arg);
		dxdy = 32*(y-0.5)*cos(arg)-256*(x-0.5)*(2*x-1)*(y-0.5)*sin(arg);
		dxdz = 32*cos(arg)*(z-0.5)-256*(x-0.5)*(2*x-1)*sin(arg)*(z-0.5);
		dydy = 16*(2*x-1)*cos(arg)-256*(2*x-1)*(y-0.5)*(y-0.5)*sin(arg);
		dydz = -256*(2*x-1)*(y-0.5)*sin(arg)*(z-0.5);
		dzdz = 16*(2*x-1)*cos(arg)-256*(2*x-1)*sin(arg)*(z-0.5)*(z-0.5);
		
		
		
		
		f = sin(8*((vx-0.5)*(vx-0.5)+(vy-0.5)*(vy-0.5)+(vz-0.5)*(vz-0.5)))*(2*vx-1);
		f2 = sin(8 * ((vx - 0.5)*(vx - 0.5) + (vy - 0.5)*(vy - 0.5) + (vz - 0.5)*(vz - 0.5)))*(2 * vx - 1);
	}
	//mixed derivative computed twice: dxdy and dydx
	dxdy *= 2;
	dxdz *= 2;
	dydz *= 2;
	
	double vdx = f.GetRow()[0];
	double vdy = f.GetRow()[1];
	double vdz = f.GetRow()[2];
	double vdxdx = f.GetHessianRow()[Sparse::HessianRow::make_index(0,0)];
	double vdydy = f.GetHessianRow()[Sparse::HessianRow::make_index(1,1)];
	double vdzdz = f.GetHessianRow()[Sparse::HessianRow::make_index(2,2)];
	double vdxdy = f.GetHessianRow()[Sparse::HessianRow::make_index(0,1)];
	double vdxdz = f.GetHessianRow()[Sparse::HessianRow::make_index(0,2)];
	double vdydz = f.GetHessianRow()[Sparse::HessianRow::make_index(1,2)];
	

	bool error = false;
	std::cout << "For GetHessian:" << std::endl;
	std::cout << std::setw(10) << "derivative " << std::setw(10) << "original " << std::setw(10) << "computed" << std::endl;
	std::cout << std::setw(10) << "dx " << std::setw(10) << dx << std::setw(10) << vdx << std::endl;
	std::cout << std::setw(10) << "dy " << std::setw(10) << dy << std::setw(10) << vdy << std::endl;
	std::cout << std::setw(10) << "dz " << std::setw(10) << dz << std::setw(10) << vdz << std::endl;
	std::cout << std::setw(10) << "dxdx " << std::setw(10) << dxdx << std::setw(10) << vdxdx << std::endl;
	std::cout << std::setw(10) << "dxdy " << std::setw(10) << dxdy << std::setw(10) << vdxdy << std::endl;
	std::cout << std::setw(10) << "dxdz " << std::setw(10) << dxdz << std::setw(10) << vdxdz << std::endl;
	std::cout << std::setw(10) << "dydy " << std::setw(10) << dydy << std::setw(10) << vdydy << std::endl;
	std::cout << std::setw(10) << "dydz " << std::setw(10) << dydz << std::setw(10) << vdydz << std::endl;
	std::cout << std::setw(10) << "dzdz " << std::setw(10) << dzdz << std::setw(10) << vdzdz << std::endl;
	if( std::abs(dx-vdx) > 1.0e-9 ) error = true, std::cout << "Error in dx: " << std::abs(dx-vdx) << " original " << dx << " computed " << vdx << std::endl;
	if( std::abs(dy-vdy) > 1.0e-9 ) error = true, std::cout << "Error in dy: " << std::abs(dy-vdy) << " original " << dy << " computed " << vdy << std::endl;
	if( std::abs(dz-vdz) > 1.0e-9 ) error = true, std::cout << "Error in dz: " << std::abs(dz-vdz) << " original " << dz << " computed " << vdz << std::endl;
	if( std::abs(dxdx-vdxdx) > 1.0e-9 ) error = true, std::cout << "Error in dxdx: " << std::abs(dxdx-vdxdx) << " original " << dxdx << " computed " << vdxdx << std::endl;
	if( std::abs(dxdy-vdxdy) > 1.0e-9 ) error = true, std::cout << "Error in dxdy: " << std::abs(dxdy-vdxdy) << " original " << dxdy << " computed " << vdxdy << std::endl;
	if( std::abs(dxdz-vdxdz) > 1.0e-9 ) error = true, std::cout << "Error in dxdz: " << std::abs(dxdz-vdxdz) << " original " << dxdz << " computed " << vdxdz << std::endl;
	if( std::abs(dydy-vdydy) > 1.0e-9 ) error = true, std::cout << "Error in dydy: " << std::abs(dydy-vdydy) << " original " << dydy << " computed " << vdydy << std::endl;
	if( std::abs(dydz-vdydz) > 1.0e-9 ) error = true, std::cout << "Error in dydz: " << std::abs(dydz-vdydz) << " original " << dydz << " computed " << vdydz << std::endl;
	if( std::abs(dzdz-vdzdz) > 1.0e-9 ) error = true, std::cout << "Error in dzdz: " << std::abs(dzdz-vdzdz) << " original " << dzdz << " computed " << vdzdz << std::endl;
	

	vdx = f2.GetRow()[0];
	vdy = f2.GetRow()[1];
	vdz = f2.GetRow()[2];
	std::cout << "For GetJacobian:" << std::endl;
	std::cout << std::setw(10) << "derivative " << std::setw(10) << "original " << std::setw(10) << "computed" << std::endl;
	std::cout << std::setw(10) << "dx " << std::setw(10) << dx << std::setw(10) << vdx << std::endl;
	std::cout << std::setw(10) << "dy " << std::setw(10) << dy << std::setw(10) << vdy << std::endl;
	std::cout << std::setw(10) << "dz " << std::setw(10) << dz << std::setw(10) << vdz << std::endl;

	if (std::abs(dx - vdx) > 1.0e-9) error = true, std::cout << "Error in dx: " << std::abs(dx - vdx) << " original " << dx << " computed " << vdx << std::endl;
	if (std::abs(dy - vdy) > 1.0e-9) error = true, std::cout << "Error in dy: " << std::abs(dy - vdy) << " original " << dy << " computed " << vdy << std::endl;
	if (std::abs(dz - vdz) > 1.0e-9) error = true, std::cout << "Error in dz: " << std::abs(dz - vdz) << " original " << dz << " computed " << vdz << std::endl;

	if (error) return -1;

	return 0;
}
