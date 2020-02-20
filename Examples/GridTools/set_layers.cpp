#include "set_layers.h"

#define ind(r,c) ((r)*N + (c))

void rand1d(double * arr, int N)
{
	if( N == 1 ) return;
	const double t = 0.15;
	double a = 0.5+(2*(rand()*1.0/RAND_MAX)-1)*t; // in 0.5+[-1,1]*t
	arr[N/2] = arr[0]*a + arr[N]*(1-a);
	rand1d(arr,N/2);
	rand1d(arr+N/2,N-N/2);
}

std::pair<int,double> layer1d(double * arr, int N, double z)
{
	for(int k = 0; k < N; ++k)
	{
		if( z >= arr[k] && z <= arr[k+1] )
			return std::make_pair(k, (z-arr[k])/(arr[k+1]-arr[k]));
	}
	return std::make_pair(-1,0.0);
}
const double noind = std::numeric_limits<double>::quiet_NaN();

void rand2d(double * arr, int N, int Nl, int Nr, int Nb, int Nt, double t)
{
	//std::cout << "Nl:Nr " << Nl <<":" << Nr << " Nb:Nt " << Nb << ":" << Nt << std::endl;
	if(  Nr - Nl < 2 && Nt - Nb < 2 ) 
	{
		//std::cout << "exit" << std::endl;
		return;
	}
	//const double t = 0.15;
	int Nk = (Nb+Nt)/2;
	int Nm = (Nl+Nr)/2;
	//std::cout << "Nk " << Nk << " Nm " << Nm << std::endl;
	double lb = arr[ind(Nb,Nl)];
	double rb = arr[ind(Nb,Nr)];
	double lt = arr[ind(Nt,Nl)];
	double rt = arr[ind(Nt,Nr)];
	if( lb != lb || rb != rb || lt != lt || rt != rt ) throw -1;
	if( arr[ind(Nk,Nl)] != arr[ind(Nk,Nl)] ) arr[ind(Nk,Nl)] = 0.5*(lb + lt) + (2*(rand()*1.0/RAND_MAX)-1)*t;
	if( arr[ind(Nk,Nr)] != arr[ind(Nk,Nr)] ) arr[ind(Nk,Nr)] = 0.5*(rb + rt) + (2*(rand()*1.0/RAND_MAX)-1)*t;
	if( arr[ind(Nb,Nm)] != arr[ind(Nb,Nm)] ) arr[ind(Nb,Nm)] = 0.5*(lb + rb) + (2*(rand()*1.0/RAND_MAX)-1)*t;
	if( arr[ind(Nt,Nm)] != arr[ind(Nt,Nm)] ) arr[ind(Nt,Nm)] = 0.5*(lt + rt) + (2*(rand()*1.0/RAND_MAX)-1)*t;
	arr[ind(Nk,Nm)] = 0.25*(lb+rb+lt+rt) + (2*(rand()*1.0/RAND_MAX)-1)*t;
	rand2d(arr,N,Nl,Nm,Nb,Nk,t*0.5); rand2d(arr,N,Nm,Nr,Nb,Nk,t*0.5);
	rand2d(arr,N,Nl,Nm,Nk,Nt,t*0.5); rand2d(arr,N,Nm,Nr,Nk,Nt,t*0.5);
}

void init2d(double * arr, int N, double mint, double maxt)
{
	for(int k = 0; k < N*N; ++k) arr[k] = noind;
	arr[ind(0  ,0  )] = mint+(rand()*1.0/RAND_MAX)*(maxt-mint);
	arr[ind(0  ,N-1)] = mint+(rand()*1.0/RAND_MAX)*(maxt-mint);
	arr[ind(N-1,0  )] = mint+(rand()*1.0/RAND_MAX)*(maxt-mint);
	arr[ind(N-1,N-1)] = mint+(rand()*1.0/RAND_MAX)*(maxt-mint);
}

double intrp2d(double * arr, int N, double x, double y)
{
	int n = ceil(x*(N-1));
	int m = ceil(y*(N-1));
	if( n == 0 ) n = 1;
	if( m == 0 ) m = 1;
	double dh = 1.0/(double)(N-1);
	double kx = (x-(n-1)*dh)/dh;
	double ky = (y-(m-1)*dh)/dh;
	//if( kx < 0 || kx > 1 ) std::cout << "bad kx: " << kx << " x is " << x << " n is " << n << " dh is " << dh << " N is " << N << std::endl;
	//if( ky < 0 || ky > 1 ) std::cout << "bad ky: " << ky << " y is " << y << " m is " << m << " dh is " << dh << " N is " << N << std::endl;
	double lb = arr[ind(m-1,n-1)];
	double rb = arr[ind(m-1,n+0)];
	double lt = arr[ind(m+0,n-1)];
	double rt = arr[ind(m+0,n+0)];
	if( lb != lb || rb != rb || lt != lt || rt != rt ) throw -1;
	return (1-ky)*(lb*(1-kx) + rb*kx) + ky*(lt*(1-kx) + rt*kx);
}

void SetLayers::RandomLayersZ(int num_layers)
{
	layers_z.resize(num_layers+1);
	layers_z[0] = 0;
	layers_z[num_layers] = 1;
	rand1d(&layers_z[0],num_layers);
}

void SetLayers::DeformHeight(double height)
{
	const int N = 128;
	std::vector<double> map(N*N,0.0);
	init2d(&map[0],N,0.0,height);
	rand2d(&map[0],N,0,N-1,0,N-1,height*0.5);
	
	for(int k = 0; k < N*N; ++k) if( map[k] != map[k] ) std::cout << "NaN at row " << k/N << " col " << k%N << std::endl;
	
	double cmin[3], cmax[3];
	for(int d = 0; d < 3; ++d)
	{
		cmin[d]=1.0e+20;
		cmax[d]=-1.0e+20;
	}
	
	for(Mesh::iteratorNode n = mesh.BeginNode(); n != mesh.EndNode(); ++n)
	{
		for(int d = 0; d < 3; ++d)
		{
			if( cmin[d] > n->Coords()[d] ) cmin[d] = n->Coords()[d];
			if( cmax[d] < n->Coords()[d] ) cmax[d] = n->Coords()[d];
		}
	}
	for(Mesh::iteratorNode n = mesh.BeginNode(); n != mesh.EndNode(); ++n)
	{
		double c[3] = {0,0,0};
		for(int d = 0; d < 3; ++d)
			c[d] = (n->Coords()[d]-cmin[d])/(cmax[d]-cmin[d]);
		n->Coords()[2] += intrp2d(&map[0],N,c[0],c[1]);
	}
}

void SetLayers::DeformLayers(double coef)
{
	const int N = 128;
	std::vector< std::vector<double> > map_layer(layers_z.size());
#if defined(USE_OMP)
#pragma omp parallel for
#endif
	for(int k = 0; k < map_layer.size(); ++k)
	{
		map_layer[k].resize(N*N);
		init2d(&map_layer[k][0],N,1-coef,1);
		rand2d(&map_layer[k][0],N,0,N-1,0,N-1,coef);
		for(int kk = 0; kk < N*N; ++kk) if( map_layer[k][kk] < 0.0 ) map_layer[k][kk] = 0;
	}
	double cmin[3], cmax[3];
	for(int d = 0; d < 3; ++d)
	{
		cmin[d]=1.0e+20;
		cmax[d]=-1.0e+20;
	}
	for(Mesh::iteratorNode n = mesh.BeginNode(); n != mesh.EndNode(); ++n)
	{
		for(int d = 0; d < 3; ++d)
		{
			if( cmin[d] > n->Coords()[d] ) cmin[d] = n->Coords()[d];
			if( cmax[d] < n->Coords()[d] ) cmax[d] = n->Coords()[d];
		}
	}
	TagInteger layer_tag = mesh.CreateTag("LAYER",DATA_INTEGER,CELL,NONE,1);
	TagReal    coef_tag  = mesh.CreateTag("LAYER_COEF",DATA_REAL,CELL,NONE,1);
	for(Mesh::iteratorCell c = mesh.BeginCell(); c != mesh.EndCell(); ++c)
	{
		double cnt[3], cc[3] = {0,0,0};
		c->Centroid(cnt);
		for(int d = 0; d < 3; ++d)
			cc[d] = (cnt[d]-cmin[d])/(cmax[d]-cmin[d]);
		std::pair<int,double> lc = layer1d(&layers_z[0],layers_z.size()-1,cc[2]);
		layer_tag[*c] = lc.first;
		coef_tag[*c] = lc.second;
		
	}
	for(Mesh::iteratorNode n = mesh.BeginNode(); n != mesh.EndNode(); ++n)
	{
		double c[3] = {0,0,0};
		for(int d = 0; d < 3; ++d)
			c[d] = (n->Coords()[d]-cmin[d])/(cmax[d]-cmin[d]);
		std::pair<int,double> lc = layer1d(&layers_z[0],layers_z.size()-1,c[2]);
		if( lc.first == -1 ) std::cout << "layer not found for " << c[0] << "," << c[1] << "," << c[2] << " node " << n->LocalID() << std::endl;
		
		double h = 0;
		int kk = lc.first;
		for(int k = 0; k < kk; ++k)
			h += (layers_z[k+1]-layers_z[k])*(cmax[2]-cmin[2])*intrp2d(&map_layer[k][0],N,c[0],c[1]);
		
		h += lc.second*(layers_z[kk+1]-layers_z[kk])*(cmax[2]-cmin[2])*intrp2d(&map_layer[kk][0],N,c[0],c[1]);
		
		n->Coords()[2] = h;
	}
	
	//TODO: fix zero edges
	int bad_edge = 0,collapse_edge =0;
	mesh.BeginModification();
	for(Mesh::iteratorEdge e = mesh.BeginEdge(); e != mesh.EndEdge(); ++e) if( !e->Hidden() )
	{
		if( e->Length() < 1.0e-9 ) 
		{
			if( !e->Hidden() )
			{
				Edge::CollapseEdge(e->self(),0);
				collapse_edge++;
			}
			bad_edge++;
		}
	}
	mesh.EndModification();
	std::cout << "Bad edges: " << bad_edge << " collapsed " << collapse_edge << std::endl;
}


