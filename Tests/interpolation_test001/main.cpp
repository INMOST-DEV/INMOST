#include "inmost.h"
#include "cstdlib"

using namespace INMOST;

int main(int argc, char ** argv)
{
	int test = 0;
	if( argc > 1 )	test = atoi(argv[1]);

	int err = 0, err2 = 0;

	Mesh m;

	if( test == 0 )
	{
		Storage::real coords[8][3] = {
			{0,0,0},
			{1,0,0},
			{0,1,0},
			{1,1,0},
			{0,0,1},
			{1,0,1},
			{0,1,1},
			{1,1,1}
		};
		ElementArray<Node> nodes(&m), enodes(&m);
		ElementArray<Edge> edges(&m), fedges(&m);
		ElementArray<Face> faces(&m);
		enodes.resize(2);
		fedges.resize(4);
		for(int k = 0; k < 8; k++)
			nodes.push_back(m.CreateNode(coords[k]));
		int eindices[12][2] = {{0,1},{2,3},{4,5},{6,7},{0,2},{1,3},{4,6},{5,7},{0,4},{1,5},{2,6},{3,7}};
		for(int k = 0; k < 12; k++)
		{
			enodes[0] = nodes[eindices[k][0]];
			enodes[1] = nodes[eindices[k][1]];
			Edge e = m.CreateEdge(enodes).first;
			edges.push_back(e);
		}
		int findices[6][4] = {{8,6,10,4},{11,7,9,5},{0,9,2,8},{3,11,1,10},{0,4,1,5},{3,6,2,7}};
		for(int k = 0; k < 6; k++)
		{
			for(int j = 0; j < 4; j++)	fedges[j] = edges[findices[k][j]];
			Face f = m.CreateFace(fedges).first;
			faces.push_back(f);
		}
		Cell c = m.CreateCell(faces).first;

		Storage::real midedges[12][3];
		for(int k = 0; k < 12; k++)
		{
			for(int j = 0; j < 3; j++)	midedges[k][j] = 0.5*(coords[eindices[k][0]][j] + coords[eindices[k][1]][j]);
		}
		Storage::real midfaces[6][3];
		for(int k = 0; k < 6; k++)
		{
			for(int j = 0; j < 3; j++)	midfaces[k][j] = 0.5*(midedges[findices[k][0]][j] + midedges[findices[k][2]][j]);
		}

		Storage::real midcell[3] = {0.5,0.5,0.5};

		std::map<HandleType,Storage::real> nodes_stencil;
		for(int j = 0; j < 8; j++)
		{
			if(!nodes_stencil.empty())	nodes_stencil.clear();
			m.WachspressInterpolation3D(coords[j]+0, c, nodes_stencil);
			if( nodes_stencil.size() != 2 ) err++;
			for(std::map<HandleType,Storage::real>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it)
				if( fabs(it->second-1.0) > 1e-8 && fabs(it->second) > 1e-8 )
				{
					std::cout << "Coef " << it->second << " expected 1.0 or 0.0 num " << j << std::endl;
					err2++;
				}
		}
		
		for(int j = 0; j < 12; j++)
		{
			if(!nodes_stencil.empty())	nodes_stencil.clear();
			m.WachspressInterpolation3D(midedges[j], c, nodes_stencil);
			if( nodes_stencil.size() != 2 ) err++;
			for(std::map<HandleType,Storage::real>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it)
				if( fabs(it->second-0.5) > 1e-8 )
				{
					std::cout << "Coef " << it->second << " expected 0.5 " << std::endl;
					err2++;
				}
		}

		for(int j = 0; j < 6; j++)
		{
			if(!nodes_stencil.empty())	nodes_stencil.clear();
			m.WachspressInterpolation3D(midfaces[j], c, nodes_stencil);
			if( nodes_stencil.size() != 4 ) err++;
			for(std::map<HandleType,Storage::real>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it)
				if( fabs(it->second-0.25) > 1e-8 )
				{
					std::cout << "Coef " << it->second << " expected 0.25 " << std::endl;
					err2++;
				}
		}

		if(!nodes_stencil.empty())	nodes_stencil.clear();
		m.WachspressInterpolation3D(midcell, c, nodes_stencil);
		if( nodes_stencil.size() != 8 ) err++;
		for(std::map<HandleType,Storage::real>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it)
			if( fabs(it->second-0.125) > 1e-8 )
			{
				std::cout << "Coef " << it->second << " expected 0.125 " << std::endl;
				err2++;
			}
		/*
		Storage::real point[3] = {0.75,0.25,0.0};
		Storage::real coefs[4] = {0.1875,0.0625,0.1875,0.5625};
		if(!nodes_stencil.empty())	nodes_stencil.clear();
		m.wachspress_2d(point, f, nodes_stencil);
		if( nodes_stencil.size() != 4 ) err++;
		int j = 0;
		for(std::map<HandleType,Storage::real>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it,++j)
			if( fabs(it->second-coefs[j]) > 1e-8 )
			{
				std::cout << "Coef " << it->second << " expected " << coefs[j] << std::endl;
				err2++;
			}
			*/
	}
	else if( test == 1 )
	{
		Storage::real coords[4][3] = {
			{1,0,0},
			{0,1,0},
			{0,0,1},
			{0,0,0}
		};
		ElementArray<Node> nodes(&m), enodes(&m);
		ElementArray<Edge> edges(&m), fedges(&m);
		ElementArray<Face> faces(&m);
		enodes.resize(2);
		fedges.resize(3);
		for(int k = 0; k < 8; k++)
			nodes.push_back(m.CreateNode(coords[k]));
		int eindices[6][2] = {{3,0},{3,1},{3,2},{0,1},{0,2},{1,2}};
		for(int k = 0; k < 6; k++)
		{
			enodes[0] = nodes[eindices[k][0]];
			enodes[1] = nodes[eindices[k][1]];
			Edge e = m.CreateEdge(enodes).first;
			edges.push_back(e);
		}
		int findices[4][3] = {{1,2,5},{0,4,2},{0,3,1},{3,5,4}};
		for(int k = 0; k < 4; k++)
		{
			for(int j = 0; j < 3; j++)	fedges[j] = edges[findices[k][j]];
			Face f = m.CreateFace(fedges).first;
			faces.push_back(f);
		}
		Cell c = m.CreateCell(faces).first;

		Storage::real point[3] = {0.25,0.2,0.3};
		Storage::real coefs[4][4], coefs_sum = 0.0;
		coefs[0][3] = 200.0/3, coefs[1][3] = 160.0/3, coefs[2][3] = 80.0, coefs[3][3] = 200.0/3;
		for(int k = 0; k < 4; k++)
		{
			for(int l = 0; l < 3; l++)	coefs[k][l] = coords[k][l];
			coefs_sum += coefs[k][3];
		}
		for(int k = 0; k < 4; k++)	coefs[k][3] /= coefs_sum;

		std::map<HandleType,Storage::real> nodes_stencil;

		if(!nodes_stencil.empty())	nodes_stencil.clear();
		m.WachspressInterpolation3D(point, c, nodes_stencil);
		if( nodes_stencil.size() != 4 ) err++;
		int k = 0;
		for(std::map<HandleType,Storage::real>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it, ++k)
		{
			Storage::real cnt[3];
			m.ElementByHandle(it->first).Centroid(cnt);
			int k = -1;
			for(int j = 0; j < 4 && k == -1; j++)
				if( sqrt( (cnt[0]-coords[j][0])*(cnt[0]-coords[j][0]) + (cnt[1]-coords[j][1])*(cnt[1]-coords[j][1]) + (cnt[2]-coords[j][2])*(cnt[2]-coords[j][2]) ) < 1e-8 )
				{
					k = j;
				}

			if( fabs(it->second-coefs[k][3]) > 1e-8 )
			{
				std::cout	<< "Coef " << it->second << " node " << cnt[0] << " " << cnt[1] << " " << cnt[2] 
								<< " expected " << coefs[k] << " num " << k << std::endl;
				err2++;
			}
		}
	}

	if( err > 0 || err2 > 0 )
		std::cout << "There are " << err << " wrong stencil sizes, " << err2 << " wrong coefs." << std::endl;
	else
		std::cout << "WachspressInterpolation3D OK" << std::endl;

	return 0;
}
