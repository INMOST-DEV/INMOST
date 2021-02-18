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
		double coords[4][3] = {
			{0,0,0},
			{0,1,0},
			{1,1,0},
			{1,0,0}
		};
		ElementArray<Node> nodes(&m), enodes(&m);
		ElementArray<Edge> edges(&m);
		enodes.resize(2);
		for(int k = 0; k < 4; k++)
			nodes.push_back(m.CreateNode(coords[k]));
		for(int k = 0; k < 4; k++)
		{
			enodes[0] = nodes[k];
			enodes[1] = nodes[(k+1)%4];
			Edge e = m.CreateEdge(enodes).first;
			edges.push_back(e);
		}
		Face f = m.CreateFace(edges).first;

		double midedges[4][3];
		for(int k = 0; k < 4; k++)
			for(int j = 0; j < 3; j++)	midedges[k][j] = 0.5*(coords[k][j] + coords[(k+1)%4][j]);
		double midface[3] = {0.5,0.5,0.5};


		std::map<HandleType,double> nodes_stencil;
		for(int j = 0; j < 4; j++)
		{
			if(!nodes_stencil.empty())	nodes_stencil.clear();
			m.WachspressInterpolation2D(coords[j]+0, f, nodes_stencil);
			if( nodes_stencil.size() != 2 ) err++;
			for(std::map<HandleType,double>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it)
				if( fabs(it->second-1.0) > 1e-8 && fabs(it->second) > 1e-8 )
				{
					std::cout << "Coef " << it->second << " expected 1.0 or 0.0 num " << j << std::endl;
					err2++;
				}
		}
		
		for(int j = 0; j < 4; j++)
		{
			if(!nodes_stencil.empty())	nodes_stencil.clear();
			m.WachspressInterpolation2D(midedges[j], f, nodes_stencil);
			if( nodes_stencil.size() != 2 ) err++;
			for(std::map<HandleType,double>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it)
				if( fabs(it->second-0.5) > 1e-8 )
				{
					std::cout << "Coef " << it->second << " expected 0.5 " << std::endl;
					err2++;
				}
		}

		if(!nodes_stencil.empty())	nodes_stencil.clear();
		m.WachspressInterpolation2D(midface, f, nodes_stencil);
		if( nodes_stencil.size() != 4 ) err++;
		for(std::map<HandleType,double>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it)
			if( fabs(it->second-0.25) > 1e-8 )
			{
				std::cout << "Coef " << it->second << " expected 0.25 " << std::endl;
				err2++;
			}
		
		double point[3] = {0.75,0.25,0.0};
		double coefs[4] = {0.1875,0.0625,0.1875,0.5625};
		if(!nodes_stencil.empty())	nodes_stencil.clear();
		m.WachspressInterpolation2D(point, f, nodes_stencil);
		if( nodes_stencil.size() != 4 ) err++;
		int j = 0;
		for(std::map<HandleType,double>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it,++j)
			if( fabs(it->second-coefs[j]) > 1e-8 )
			{
				std::cout << "Coef " << it->second << " expected " << coefs[j] << std::endl;
				err2++;
			}
	}
	else if( test == 1 )
	{
		double coords[5][3] = {
			{0,0,0},
			{0,1,0},
			{1,1,0},
			{1.5,0.5,0},
			{1,0,0}
		};
		ElementArray<Node> nodes(&m), enodes(&m);
		ElementArray<Edge> edges(&m);
		enodes.resize(2);
		for(int k = 0; k < 5; k++)
			nodes.push_back(m.CreateNode(coords[k]));
		for(int k = 0; k < 5; k++)
		{
			enodes[0] = nodes[k];
			enodes[1] = nodes[(k+1)%5];
			Edge e = m.CreateEdge(enodes).first;
			edges.push_back(e);
		}
		Face f = m.CreateFace(edges).first;

		std::map<HandleType,double> nodes_stencil;
		
		double point[3] = {0.4,0.6,0.0};
		double coefs[5] = {25.0/6, 25.0/4, 2.5, 5./3, 25./18}, sum = 0;
		for(int k = 0; k < 5; k++)	sum += coefs[k];
		for(int k = 0; k < 5; k++)	coefs[k] /= sum;
		if(!nodes_stencil.empty())	nodes_stencil.clear();
		m.WachspressInterpolation2D(point, f, nodes_stencil);
		if( nodes_stencil.size() != 5 ) err++;
		int k = 0;
		for(std::map<HandleType,double>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it, ++k)
			if( fabs(it->second-coefs[k]) > 1e-8 )
			{
				std::cout << "Coef " << it->second << " expected " << coefs[k] << std::endl;
				err2++;
			}
	}

	if( err > 0 || err2 > 0 )
		std::cout << "There are " << err << " wrong stencil sizes, " << err2 << " wrong coefs." << std::endl;
	else
		std::cout << "WachspressInterpolation2D OK" << std::endl;

	return 0;
}
