#include "inmost.h"

using namespace INMOST;


static bool MatchPoints(const Storage::real v1[3], const Storage::real v2[3])
{
	Storage::real l = 0, l2 = 0;
	for (int k = 0; k < 3; ++k)
	{
		l += (v1[k] - v2[k])*(v1[k] - v2[k]);
		l2 +=(v1[k]+v2[k])*(v1[k]+v2[k]);
	}
	return sqrt(l) < 1.0e-13*sqrt(l2/4.0);
}

static Storage::real DistPoints(const Storage::real v1[3], const Storage::real v2[3])
{
	Storage::real l = 0;
	for (int k = 0; k < 3; ++k)
		l += (v1[k] - v2[k])*(v1[k] - v2[k]);
	return sqrt(l);
}

int main(int argc, char ** argv)
{
	if( argc > 1 )
	{
		Mesh m;
		m.Load(argv[1]);
		std::string fout = "out.stl";
		if( argc > 2 ) fout = std::string(argv[2]);

		std::fstream f(fout.c_str(),std::ios::out);

		Storage::real min[3] = {1.0e+20,1.0e+20,1.0e+20}, max[3] = {-1.0e+20,-1.0e+20,-1.0e+20};

		for(Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it)
		{
			Storage::real_array c = it->Coords();
			for(int k = 0; k < 3; ++k)
			{
				if( c[k] > max[k] ) max[k] = c[k];
				if( c[k] < min[k] ) min[k] = c[k];
			}
		}

		std::cout << "mesh bounds: " << std::endl;
		const char xyz[3] = {'x','y','z'};
		for(int k = 0; k < 3; ++k)
			std::cout << xyz[k] << ": [" << min[k] << "," << max[k] << "]" << std::endl;


		f << "solid outer_space" << std::endl;

		f << std::scientific;
		f << std::setprecision(35);
		std::cout << std::scientific;
		std::cout << std::setprecision(35);
	

		//internal part
		int nontri = 0, q = 0;
		for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it) if( it->Boundary() )
		{
			Storage::real n[3];
			it->UnitNormal(n);
			
			ElementArray<Node> nodes = it->getNodes();
			
			if( it->Area() < 1.0e-8 ) std::cout << "Area: " << it->Area() << std::endl;
			
			if( it->GetGeometricType() == Element::Tri )
			{
				Storage::real_array c0 = nodes[2].Coords();
				Storage::real_array c1 = nodes[1].Coords();
				Storage::real_array c2 = nodes[0].Coords();
				
				if( MatchPoints(c0.data(),c1.data()) || MatchPoints(c0.data(),c2.data()) || MatchPoints(c1.data(),c2.data()) ) continue;
				
				f << "facet normal " << -n[0] << " " << -n[1] << " " << -n[2] << std::endl;
				f << "outer loop" << std::endl;
				
				f << "vertex " << c0[0] << " " << c0[1] << " " << c0[2] << std::endl;
				f << "vertex " << c1[0] << " " << c1[1] << " " << c1[2] << std::endl;
				f << "vertex " << c2[0] << " " << c2[1] << " " << c2[2] << std::endl;

				f << "endloop" << std::endl;
				f << "endfacet" << std::endl;
				
				
				
				
				q++;
				
				if( q == 34787 || q == 54053 )
				{
					std::cout << "vertex " << c0[0] << " " << c0[1] << " " << c0[2] << " c0c1 " << DistPoints(c0.data(),c1.data()) << " c0c2 " << DistPoints(c0.data(),c2.data()) << std::endl;
					std::cout << "vertex " << c1[0] << " " << c1[1] << " " << c1[2] << " c1c0 " << DistPoints(c0.data(),c1.data()) << " c1c2 " << DistPoints(c1.data(),c2.data()) << std::endl;
					std::cout << "vertex " << c2[0] << " " << c2[1] << " " << c2[2] << " c2c0 " << DistPoints(c0.data(),c2.data()) << " c2c1 " << DistPoints(c1.data(),c2.data()) <<  std::endl;
				}
				
				
			}
			else
			{
				std::cout << "non-tri begin: " << q << std::endl;
				nontri++;
				Storage::real c0[3];
				it->Centroid(c0);

				for(INMOST_DATA_ENUM_TYPE k = 0; k < nodes.size(); ++k)
				{
					
					//inverse normal
					f << "facet normal " << -n[0] << " " << -n[1] << " " << -n[2] << std::endl;
					f << "outer loop" << std::endl;

					Storage::real_array c2 = nodes[k].Coords();
					Storage::real_array c1 = nodes[(k+1)%nodes.size()].Coords();
					
					if( MatchPoints(c0,c1.data()) || MatchPoints(c0,c2.data()) || MatchPoints(c1.data(),c2.data())) continue;
					
					f << "vertex " << c0[0] << " " << c0[1] << " " << c0[2] << std::endl;
					f << "vertex " << c1[0] << " " << c1[1] << " " << c1[2] << std::endl;
					f << "vertex " << c2[0] << " " << c2[1] << " " << c2[2] << std::endl;

					f << "endloop" << std::endl;
					f << "endfacet" << std::endl;
					q++;
				}
				std::cout << "non-tri end: " << q << std::endl;
			}
		}
		std::cout << "Non-triangles: " << nontri << std::endl;
		
		
		//outer box
		for(int k = 0; k < 3; ++k)
		{
			Storage::real d = (max[k]-min[k]);
			min[k] -= d*0.3;
				max[k] += d*0.3;
		}
		
		//-Z

		f << "facet normal 0.0E+000 0.0E+000 -1.0E+000" << std::endl;
		f << "outer loop" << std::endl;
		f << "vertex " << min[0] << " " << min[1] << " " << min[2] << std::endl;
		f << "vertex " << min[0] << " " << max[1] << " " << min[2] << std::endl;
		f << "vertex " << max[0] << " " << max[1] << " " << min[2] << std::endl;
		f << "endloop" << std::endl;
		f << "endfacet" << std::endl;

		f << "facet normal 0.0E+000 0.0E+000 -1.0E+000" << std::endl;
		f << "outer loop" << std::endl;
		f << "vertex " << min[0] << " " << min[1] << " " << min[2] << std::endl;
		f << "vertex " << max[0] << " " << max[1] << " " << min[2] << std::endl;
		f << "vertex " << max[0] << " " << min[1] << " " << min[2] << std::endl;
		f << "endloop" << std::endl;
		f << "endfacet" << std::endl;
		
		//+Z
		
		f << "facet normal 0.0E+000 0.0E+000 1.0E+000" << std::endl;
		f << "outer loop" << std::endl;
		f << "vertex " << min[0] << " " << min[1] << " " << max[2] << std::endl;
		f << "vertex " << max[0] << " " << max[1] << " " << max[2] << std::endl;
		f << "vertex " << min[0] << " " << max[1] << " " << max[2] << std::endl;
		f << "endloop" << std::endl;
		f << "endfacet" << std::endl;

		f << "facet normal 0.0E+000 0.0E+000 1.0E+000" << std::endl;
		f << "outer loop" << std::endl;
		f << "vertex " << min[0] << " " << min[1] << " " << max[2] << std::endl;
		f << "vertex " << max[0] << " " << min[1] << " " << max[2] << std::endl;
		f << "vertex " << max[0] << " " << max[1] << " " << max[2] << std::endl;
		f << "endloop" << std::endl;
		f << "endfacet" << std::endl;
		
		//-X
		
		f << "facet normal -1.0E+000 0.0E+000 0.0E+000" << std::endl;
		f << "outer loop" << std::endl;
		f << "vertex " << min[0] << " " << min[1] << " " << min[2] << std::endl;
		f << "vertex " << min[0] << " " << min[1] << " " << max[2] << std::endl;
		f << "vertex " << min[0] << " " << max[1] << " " << max[2] << std::endl;
		f << "endloop" << std::endl;
		f << "endfacet" << std::endl;
		
		f << "facet normal -1.0E+000 0.0E+000 0.0E+000" << std::endl;
		f << "outer loop" << std::endl;
		f << "vertex " << min[0] << " " << min[1] << " " << min[2] << std::endl;
		f << "vertex " << min[0] << " " << max[1] << " " << max[2] << std::endl;
		f << "vertex " << min[0] << " " << max[1] << " " << min[2] << std::endl;
		f << "endloop" << std::endl;
		f << "endfacet" << std::endl;
		
		//+X
		
		f << "facet normal 1.0E+000 0.0E+000 0.0E+000" << std::endl;
		f << "outer loop" << std::endl;
		f << "vertex " << max[0] << " " << min[1] << " " << min[2] << std::endl;
		f << "vertex " << max[0] << " " << max[1] << " " << max[2] << std::endl;
		f << "vertex " << max[0] << " " << min[1] << " " << max[2] << std::endl;
		f << "endloop" << std::endl;
		f << "endfacet" << std::endl;
		
		f << "facet normal 1.0E+000 0.0E+000 0.0E+000" << std::endl;
		f << "outer loop" << std::endl;
		f << "vertex " << max[0] << " " << min[1] << " " << min[2] << std::endl;
		f << "vertex " << max[0] << " " << max[1] << " " << min[2] << std::endl;
		f << "vertex " << max[0] << " " << max[1] << " " << max[2] << std::endl;
		f << "endloop" << std::endl;
		f << "endfacet" << std::endl;
		
		//-Y
		
		f << "facet normal 0.0E+000 -1.0E+000 0.0E+000" << std::endl;
		f << "outer loop" << std::endl;
		f << "vertex " << min[0] << " " << min[1] << " " << min[2] << std::endl;
		f << "vertex " << max[0] << " " << min[1] << " " << max[2] << std::endl;
		f << "vertex " << min[0] << " " << min[1] << " " << max[2] << std::endl;
		f << "endloop" << std::endl;
		f << "endfacet" << std::endl;
		
		f << "facet normal 0.0E+000 -1.0E+000 0.0E+000" << std::endl;
		f << "outer loop" << std::endl;
		f << "vertex " << min[0] << " " << min[1] << " " << min[2] << std::endl;
		f << "vertex " << max[0] << " " << min[1] << " " << min[2] << std::endl;
		f << "vertex " << max[0] << " " << min[1] << " " << max[2] << std::endl;
		f << "endloop" << std::endl;
		f << "endfacet" << std::endl;
		
		//+Y
		
		f << "facet normal 0.0E+000 1.0E+000 0.0E+000" << std::endl;
		f << "outer loop" << std::endl;
		f << "vertex " << min[0] << " " << max[1] << " " << min[2] << std::endl;
		f << "vertex " << min[0] << " " << max[1] << " " << max[2] << std::endl;
		f << "vertex " << max[0] << " " << max[1] << " " << max[2] << std::endl;
		f << "endloop" << std::endl;
		f << "endfacet" << std::endl;
		
		f << "facet normal 0.0E+000 1.0E+000 0.0E+000" << std::endl;
		f << "outer loop" << std::endl;
		f << "vertex " << min[0] << " " << max[1] << " " << min[2] << std::endl;
		f << "vertex " << max[0] << " " << max[1] << " " << max[2] << std::endl;
		f << "vertex " << max[0] << " " << max[1] << " " << min[2] << std::endl;
		f << "endloop" << std::endl;
		f << "endfacet" << std::endl;


		f << "endsolid outer_space" << std::endl;

		f.close();
	}
	else std::cout << "Usage: " << argv[0] << " mesh [out.stl]" << std::endl;
}
