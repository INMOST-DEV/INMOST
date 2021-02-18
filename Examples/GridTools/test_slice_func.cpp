#include "inmost.h"
#include "slice_func.h"
using namespace INMOST;



class SliceTest : public Slice
{
	int testn;
	
	Storage::real func(Storage::real x, Storage::real y, Storage::real z, int n) const
	{
				
		if( n == 0 )
			return sqrt(x*x+y*y)-0.25;
		else if( n == 1 )
			return -(sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))-0.5);
		else if( n == 2 )
		{
			if( x > 0.5 )
				return -(sqrt((x-0.7)*(x-0.7)+(y-0.2)*(y-0.2)+(z-0.4)*(z-0.4))-0.4);
			else
				return sqrt((x-0.3)*(x-0.3)+(y-0.4)*(y-0.4)+(z-0.5)*(z-0.5))-0.3;
		}
		else if( n == 3 )
			return y-(4*x*x*x*x-2*x*x+0.5)+z*z*z;
		else if( n == 4 )
			return y-(10*x*x*x*x-8*x*x*x-5*x*x+0.2+4*x)+z*z*z;
		else if( n == 5 )
		{
			
			//double Lx = 0.2, Rx = 0.4, Ly = 0.1, Ry = 0.3;
			/*
			double Lx = 0., Rx = 2, Ly = 0.0, Ry = 4.5;
			if (x > Rx){
				if (y > Ry) return -sqrt( (x-Rx)*(x-Rx) + (y-Ry)*(y-Ry) );
				else return -(x-Rx);
			}
			if (x < Lx){
				if (y < Ly)  return -sqrt( (x-Lx)*(x-Lx) + (y-Ly)*(y-Ly) );
				else return (x - Lx);
			}
			if (y > Ry) return Ry - y;
			if (y < Ly) return y - Ly;
			return fmin( fmin(x-Lx,Rx-x), fmin(y-Ly, Ry-y) );
			*/
			
				double Lx = 0., Rx = 20, Ly = -2.5, Ry = 17.5;
			if (x > Rx)
			{
				if (y > Ry) return -sqrt( (x-Rx)*(x-Rx) + (y-Ry)*(y-Ry) );
				if (y < Ly) return -sqrt( (x-Rx)*(x-Rx) + (y-Ly)*(y-Ly) );
				return -(x-Rx);
			}
			if (x < Lx)
			{
				if (y < Ly)  return -sqrt( (x-Lx)*(x-Lx) + (y-Ly)*(y-Ly) );
				if (y > Ry)  return -sqrt( (x-Lx)*(x-Lx) + (y-Ry)*(y-Ry) );
				else return (x - Lx);
			}
			if (y > Ry) return Ry - y;
			if (y < Ly) return y - Ly;
			return fmin( fmin(x-Lx,Rx-x), fmin(y-Ly, Ry-y) );
			
		}
		else if( n == 6 )
		{
			return std::min(std::min(sqrt((x-0.48)*(x-0.48) + (z-1)*(z-1))-0.25,sqrt((y-0.53)*(y-0.53) + (z-3)*(z-3))-0.24),sqrt((x-0.6)*(x-0.6) + (z-5)*(z-5))-0.3);
		}
		return 1;
	}
public:
	SliceTest(int testn) :Slice(), testn(testn) {}
	SliceTest(const SliceTest &b) :Slice(b), testn(b.testn) {}
	SliceTest & operator =(SliceTest const & b) { Slice::operator =(b); testn = b.testn; return *this;}
	Storage::real LevelFunction(Storage::real p[3]) const {return func(p[0],p[1],p[2],testn);}
};



int main(int argc, char ** argv)
{
	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " mesh  [mesh_out=grid.pmf] [type=perforated_strip]" << std::endl;
		return -1;
	}
	
	std::string grid_out = "grid.pmf";

	int ntype = 0;
	if( argc > 2 ) grid_out = std::string(argv[2]);
	if( argc > 3 ) ntype = atoi(argv[3]);
	
	

	Mesh m;
	m.Load(argv[1]);
	//m.SetTopologyCheck(NEED_TEST_CLOSURE|PROHIBIT_MULTILINE|PROHIBIT_MULTIPOLYGON|GRID_CONFORMITY|DEGENERATE_EDGE|DEGENERATE_FACE|DEGENERATE_CELL | FACE_EDGES_ORDER | MARK_ON_ERROR | ADJACENT_DUPLICATE);
	//m.RemTopologyCheck(THROW_EXCEPTION);
	
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	
	SliceTest(ntype).SliceMesh(m);
	
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	    
	m.Save(grid_out);
	return 0;
}
