#include "inmost.h"
#include "slice_func.h"
#include "fracture.h"
#include "fix_faults.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class SliceFault : public Slice
{
	double v[3]; // vector for projection in z direction
	std::vector< std::pair<double,double> > curvxy; //a curve that defines a fault
public:
	SliceFault(double cmin[3], double cmax[3]) :Slice() 
	{
		v[0] = 0;//0.1*(rand()*1.0/RAND_MAX);
		v[1] = 0;//0.1*(rand()*1.0/RAND_MAX);
		v[2] = 1;
		double l = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
		v[0] /= l; v[1] /= l; v[2] /= l;
		double dx = (cmax[0]-cmin[0])/50.0;
		double dy = (cmax[1]-cmin[1])/50.0;
		double dl = sqrt(dx*dx+dy*dy);
		double t = 2.0*(rand()*1.0/RAND_MAX)-1.0;
		double a = M_PI*0.5 + t*M_PI*0.15;
		double y = cmin[1];
		double x = (0.2 + (rand()*1.0/RAND_MAX)*0.6)*(cmax[0]-cmin[0]) + cmin[0];
		std::cout << " starting point " << x << ", " << y << std::endl;
		curvxy.push_back(std::make_pair(x,y));
		bool ended = false;
		while( !ended )
		{
			x += cos(a)*dl;
			y += sin(a)*dl;
			if( x >= cmax[0] ) {x = cmax[0]; ended = true;}
			if( x <= cmin[0] ) {x = cmin[0]; ended = true;}
			if( y >= cmax[1] ) {y = cmax[1]; ended = true;}
			if( y <= cmin[1] ) {y = cmin[1]; ended = true;}
			curvxy.push_back(std::make_pair(x,y));
			t = 2.0*(rand()*1.0/RAND_MAX)-1.0;
			a = a + M_PI*0.05*t;
		}
		std::cout << " ending point " << x << ", " << y << std::endl;
		std::cout << " generated points " << curvxy.size() << std::endl;
	}
	SliceFault(const SliceFault &b) :Slice(b) { std::copy(b.v,b.v+3,v); }
	SliceFault & operator =(SliceFault const & b) { Slice::operator =(b); std::copy(b.v,b.v+3,v); return *this;}
	double LevelFunction(double p[3]) const 
	{
		//project p to p[2] = 0;
		double px = p[0] - v[0]/v[2]*p[2], py = p[1] - v[0]/v[2]*p[2];
		double d[2], r[2], l, lmin = 1.0e+20, smin = 0;
		for(size_t k = 0; k < curvxy.size()-1; ++k)
		{
			d[0] = curvxy[k+1].first  - curvxy[k].first;
			d[1] = curvxy[k+1].second - curvxy[k].second;
			r[0] = px - (curvxy[k+1].first  + curvxy[k].first )*0.5;
			r[1] = py - (curvxy[k+1].second + curvxy[k].second)*0.5;
			l = sqrt(r[0]*r[0]+r[1]*r[1]);
			if( l < lmin ) 
			{
				lmin = l;
				smin = d[0]*r[1] - d[1]*r[0];
			}
		}
		return smin < 0.0 ? -lmin : lmin;
	}
	
	void Getv(double vv[3])
	{
		vv[0] = v[0];
		vv[1] = v[1];
		vv[2] = v[2];
	}
};


int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " mesh [mesh_out=grid.vtk] [nfaults=1]" << std::endl;
		return -1;
	}

	std::string grid_out = "grid.vtk";
	if (argc > 2) grid_out = std::string(argv[2]);
	int nfaults = argc > 3 ? atoi(argv[3]) : 1;

	{
		Mesh m;
		m.Load(argv[1]);
		
		double cmin[3], cmax[3];
		for(int d = 0; d < 3; ++d)
		{
			cmin[d]=1.0e+20;
			cmax[d]=-1.0e+20;
		}
		
		for(Mesh::iteratorNode n = m.BeginNode(); n != m.EndNode(); ++n)
		{
			for(int d = 0; d < 3; ++d)
			{
				if( cmin[d] > n->Coords()[d] ) cmin[d] = n->Coords()[d];
				if( cmax[d] < n->Coords()[d] ) cmax[d] = n->Coords()[d];
			}
		}
		
		for(int kk = 0; kk < nfaults; ++kk)
		{
			SliceFault sf(cmin,cmax);
			sf.SliceMesh(m,false);
			//create fake aperture tag
			TagReal aperture = m.CreateTag("APERTURE",DATA_REAL,FACE,FACE,1);
			TagBulk sliced = m.GetTag("SLICED"); //this was created by Slice
			TagInteger mat = m.GetTag("MATERIAL");
		
			int nmrk = 0;
			for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it) if( it->HaveData(sliced) ) {aperture[*it] = 1; nmrk++;}
		
			std::cout << "marked for fracture " << nmrk << std::endl;
		
			Fracture fr(m);
			//fr.Open(aperture,false,0.5);
			fr.Open(aperture,false,1.0);
		
			//shift nodes
			
			
			double v[3];
			sf.Getv(v);
			
			double t = 0.1+0.9*rand()*1.0/RAND_MAX;
			double sc = t*(cmax[2]-cmin[2])*0.05;
			
			std::cout << "shift nodes by " << sc << std::endl;
			
			
			for(Mesh::iteratorNode n = m.BeginNode(); n != m.EndNode(); ++n)
			{
				ElementArray<Cell> adj_cells = n->getCells();
				int mats[3] = {0,0,0};
				for(ElementArray<Cell>::iterator kt = adj_cells.begin(); kt != adj_cells.end(); ++kt)
					mats[mat[*kt]]++;
				if( mats[1] == 0 )
				{
					n->Coords()[0] += v[0]*sc;
					n->Coords()[1] += v[1]*sc;
					n->Coords()[2] += v[2]*sc;
				}
			}
			
			
			FixFaults fix(m);	
			MarkerType fault = m.CreateMarker();
			std::cout << "fault marker: " << fault << std::endl;
			for(Mesh::iteratorFace f = m.BeginFace(); f != m.EndFace(); ++f)
				if( mat[*f] == 2 ) f->SetMarker(fault);
			fix.FixMeshFaults(fault);
			std::cout << "release fault marker: " << fault << std::endl;
			m.ReleaseMarker(fault,FACE);
			
			//delete to clear data
			m.DeleteTag(aperture);
			m.DeleteTag(sliced);
			m.DeleteTag(mat);
		}
		
		m.Save(grid_out);
	}
	return 0;
}
