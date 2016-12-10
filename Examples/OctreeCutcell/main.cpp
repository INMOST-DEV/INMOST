//Compilation g++ *.cpp -D__GRAPHICS__ -L/usr/X11R6/lib -lX11 -lXi -lXmu -lGL -lglut -lGLU ../../INMOST.a -O5
//run ./a.out 
//move cursor around the screen to see grid refined and field interpolated
//quadtree for default, For an octree compile with #define DIM 3 in file octgrid.h


// in obj.cpp z coordinate is not normalized from 0 to 1!
// may be a problem, lines 791-794 
#include "inmost.h"
#include "octgrid.h"
#include "my_glut.h"
#include "rotate.h"
#include "proj.hpp"
#include "gettype.hpp"
#include "obj.h"
#include <math.h>
#include <sstream>
GetType get_type;
Projection make_proj;
bool allow_coarse = 1;
bool allow_refine = 1;
const int draw_edges = 1;
const int draw_regions = 1;
const int draw_faces = 0;
const int draw_nodes = 1;
int recreate = 1;
int level_max = 2; 
int level_max_well = 0;
int draw_quantity = 0;
int allow_amr = 0;
int draw_scheme = 0;
double mx = 0, my = 0;
Tag pressure;
Tag saturation;
Tag gas;

struct grid thegrid;
int print_element = 0;
int initialize = 1;
bool wasmotion = false;
double last_time;
double gleft, gright, gbottom, gtop, gnear, gfar, zoom = 1;
bool transparent = false;

//#define MASAHIKO_MESH
//#define SPHERE_MESH
#define OBJ_READERS
#define INIT_PERM


void transformation(double xyz[3]) 
{
#if defined(SPHERE_MESH)
	return;
#elif defined(MASAHIKO_MESH)
	if( xyz[2] < 11.0/36.0 )
		xyz[2] = (xyz[2]-0)*11.9/11.0;
	else if( xyz[2] < 12.0/36.0 )
		xyz[2] = (xyz[2]-11.0/36.0)*0.1 + 11.9/36.0;
	else if( xyz[2] < 23.0/36.0 )
		xyz[2] = (xyz[2]-12.0/36.0)*11.9/11.0 + 12.0/36.0;
	else if( xyz[2] < 24.0/36.0 )
		xyz[2] = (xyz[2]-23.0/36.0)*0.1 + 23.9/36.0;
	return;
#endif
	make_proj.Project(xyz);
	/*
    double tmp[3];
	tmp[0] = xyz[0];
	tmp[1] = xyz[1];
	tmp[2] = xyz[2];
	tmp[0] -= 0.5;
	tmp[0] *= 100.0;
	tmp[1] -= 0.5;
	tmp[1] *= 100.0;
	xyz[0] =  tmp[0];
	xyz[1] =  tmp[1];
	xyz[2] = 4010.0 + 10.0 * (tmp[2]-0.5);
	*/
}


mat_ret_type get_material_types(double xyz[3])
{
	mat_ret_type ret1;
	ret1.push_back(1);
	return ret1;
#if defined(MASAHIKO_MESH)
	mat_ret_type ret;
	//ret.push_back(1);
	
	const double eps = 1.0e-5;
	double test = xyz[2] - (2*xyz[0]-0.5);
	if( test < eps )  ret.push_back(1);
	if( test > -eps ) ret.push_back(2);
	
	return ret;
#elif defined(SPHERE_MESH)
	mat_ret_type ret;
	const double eps = 1.0e-5;
	double sphere = 0.5 - sqrt((xyz[0]-0.5)*(xyz[0]-0.5) + (xyz[1]-0.5)*(xyz[1]-0.5) + (xyz[2]-0.5)*(xyz[2]-0.5));
	//double sphere = 1.0 - sqrt((xyz[0])*(xyz[0]) + (xyz[1])*(xyz[1]) + (xyz[2])*(xyz[2]));
	if( sphere < eps ) ret.push_back(-1);
	if( sphere > -eps ) ret.push_back(1);
	return ret;
#else
	mat_ret_type ret(get_type.Size());
	ret.resize(get_type.Layers(xyz,&ret[0]));
	std::sort(ret.begin(),ret.end());
	ret.resize(std::unique(ret.begin(),ret.end())-ret.begin());
	return ret;
#endif
}




int cell_should_unite(struct grid * g, int cell)
{
#if defined(MASAHIKO_MESH)
	return 0;
	//refinement along well
	double x,y,z;
	x = g->cells[cell].center[0];
	y = g->cells[cell].center[1];
	z = g->cells[cell].center[2];
	double plane = fabs(z - (1.25-2*x));
	double mid = fabs(y-0.5);
	return plane > 0.06 && mid > 0.06;

#endif
	const double r = 0.03;
	//const double r = 0.001;
	double r2 = (1.0/(double)g->n[0]*1.0/(double)g->n[0]+1.0/(double)g->n[1]*1.0/(double)g->n[1])/2.0;
	int test = 1;	
	//~ test &= (g->cells[cell].center[0]-0)*(g->cells[cell].center[0]-0)+(g->cells[cell].center[1]-0)*(g->cells[cell].center[1]-0) > r2/(double)pow(g->cells[cell].level+1,2.5);
	//~ test &= (g->cells[cell].center[0]-1)*(g->cells[cell].center[0]-1)+(g->cells[cell].center[1]-1)*(g->cells[cell].center[1]-1) > r2/(double)pow(g->cells[cell].level+1,2.5); 
	
	test &= 0.0//(g->cells[cell].center[2]-0.5)*(g->cells[cell].center[2]-0.5)
		+(g->cells[cell].center[0]-mx)*(g->cells[cell].center[0]-mx)+(g->cells[cell].center[1]-my)*(g->cells[cell].center[1]-my) > r;
	//if( saturation_grad(g,cell) > r ) test = 0;
	return test;
}
int cell_should_split(struct grid * g, int cell)
{
#if defined(MASAHIKO_MESH)
	return 0;
	//refinement along well
	double x,y,z;
	x = g->cells[cell].center[0];
	y = g->cells[cell].center[1];
	z = g->cells[cell].center[2];
	double plane = fabs(z - (1.25-2*x));
	double mid = fabs(y-0.5);
	return plane < 0.05 && mid < 0.05 && g->cells[cell].level < level_max;

#endif	
	const double r = 0.03;
	//const double r = 0.001;
	double r2 = (1.0/(double)g->n[0]*1.0/(double)g->n[0]+1.0/(double)g->n[1]*1.0/(double)g->n[1])/2.0;
	int test = 0, test2 = 0;
	//~ test2 |= (g->cells[cell].center[0]-0)*(g->cells[cell].center[0]-0)+(g->cells[cell].center[1]-0)*(g->cells[cell].center[1]-0) < r2/(double)pow(g->cells[cell].level+1,2.5);
	//~ test2 |= (g->cells[cell].center[0]-1)*(g->cells[cell].center[0]-1)+(g->cells[cell].center[1]-1)*(g->cells[cell].center[1]-1) < r2/(double)pow(g->cells[cell].level+1,2.5);

	test |= //(g->cells[cell].center[2]-0.5)*(g->cells[cell].center[2]-0.5) + 
		(g->cells[cell].center[0]-mx)*(g->cells[cell].center[0]-mx)+(g->cells[cell].center[1]-my)*(g->cells[cell].center[1]-my) < r;	
	//test |= saturation_grad(g,cell) > r;
	   
	if( (test && !test2 && g->cells[cell].level < level_max) 
		|| 
		(test2 && g->cells[cell].level < level_max_well) 
		)
		return 1;
	return 0;
	
	/*
	int test = 0;
	double x,y,z;
	x = g->cells[cell].center[0];
	y = g->cells[cell].center[1];
	z = g->cells[cell].center[2];

	test |= (x>0.1 && x < 0.5 && y>0.3 && y < 0.6 && z > 0.4 && z < 1.);
  
	if( test && g->cells[cell].level < level_max)
	return 1;
	return 0;
	*/
}


void cell_init_data(struct grid * g, int cell)
{
	g->cells[cell].data->clear();
	if( !g->cells[cell].mr->empty() )
	{
		const int ntags = 2;
		Tag tags[ntags];
		tags[0] = pressure;
		tags[1] = saturation;
		
		for(int i = 0; i < g->cells[cell].mr->size(); i++)
		{
			Storage::real vol = Cell(g->mesh,(*g->cells[cell].mr)[i])->Volume();
			Storage::integer mat = Cell(g->mesh,(*g->cells[cell].mr)[i])->Integer(g->cell_material);
			vol_and_data_by_mat & current_mat = (*g->cells[cell].data)[mat];
			for(int j = 0; j < ntags; j++)
			{
				Storage::real val = Cell(g->mesh,(*g->cells[cell].mr)[i])->Real(tags[j]);
				current_mat.second[tags[j]] += vol*val;
				//~ if( isnan(current_mat.second[tags[j]]) ) throw -1;
			}
			current_mat.first += vol;
		}
	}
}

void cell_destroy_data(struct grid * g, int cell)
{
	g->cells[cell].data->clear();
}

void cell_unite_data(struct grid * g, int cell)
{
	int i;
	if( !g->cells[cell].leaf )
	{
		const int ntags = 2;
		Tag tags[ntags];
		tags[0] = pressure;
		tags[1] = saturation;
		g->cells[cell].data->clear();
		for(i = 0; i < 1<<DIM; i++)
		{
			for(data_storage::iterator it = g->cells[g->cells[cell].children[i]].data->begin(); it != g->cells[g->cells[cell].children[i]].data->end(); ++it)
			{
				vol_and_data_by_mat & my_current_mat = (*g->cells[cell].data)[it->first];
				for(int j = 0; j < ntags; j++)
				{
					my_current_mat.second[tags[j]] += it->second.second[tags[j]];
					//~ if( isnan(my_current_mat.second[tags[j]]) ) throw -1; 
				}
				my_current_mat.first += it->second.first;
			}
		}
	}
}


void cell_split_data(struct grid * g, int cell)
{
	//~ g->cells[cell].data.clear();
}

Storage::real get_mean(std::pair< Storage::real, Storage::real > & data)
{
	return data.first/data.second;
}

Storage::real get_sat_by_mat(Storage::integer mat)
{
	switch(mat)
	{
		case 0: return 0.15; break;
		case 1: return 0.35; break;
		case 2: return 0.05; break;
		case 3: return 1; break;
		case 4: return 0.25; break;
		case 5: return 0.15; break;
	}
	return 0;
}
Storage::real get_press_by_mat(Storage::integer mat)
{
	return 100*mat/(Storage::real)get_type.Size();
}

void fill_K(Storage::real * center, Storage::real * Kvec, Storage::real * K)
{
	return;
	Storage::real Knrm[3];
	if( !get_type.Normal(center,Knrm) )
	{
		Knrm[0] = 0;
		Knrm[1] = 0;
		Knrm[2] = 1;
	}
	else
	{
		double l = sqrt(Knrm[0]*Knrm[0]+Knrm[1]*Knrm[1]+Knrm[2]*Knrm[2]);
		Knrm[0] /= l;
		Knrm[1] /= l;
		Knrm[2] /= l;
	}
	memcpy(Kvec,Knrm,sizeof(double)*3);
	double qmat[16], qrmat[16], sclmat[16], qnrm;
	double kx = 100, ky = 100, kz = 10;
	struct quat q, qr;
	q.vec[0] = - Knrm[1];
	q.vec[1] = Knrm[0];
	q.vec[2] = 0;
	q.w = 1.0 + Knrm[2];
	qnrm = quatnorm(q);
	q.vec[0] /= qnrm;
	q.vec[1] /= qnrm;
	q.w /= qnrm;
	qr = quatconj(q);
	
	qmat[0] = (q.w*q.w + q.vec[0]*q.vec[0] - q.vec[1]*q.vec[1]) *qnrm;
	qmat[1] = 2.*(q.vec[0]*q.vec[1] ) *qnrm                          ;
	qmat[2] = 2.*(q.w*q.vec[1]) *qnrm                                ;
	
	qmat[3] = 2.*(q.vec[0]*q.vec[1]) *qnrm                           ;
	qmat[4] = (q.w*q.w - q.vec[0]*q.vec[0] + q.vec[1]*q.vec[1]) *qnrm;
	qmat[5] = 2.*( - q.w*q.vec[0])*qnrm                              ;
	
	qmat[6] = 2.*( - q.w*q.vec[1])*qnrm                              ;
	qmat[7] = 2.*( + q.w*q.vec[0])*qnrm                              ;
	qmat[8] = (q.w*q.w - q.vec[0]*q.vec[0] - q.vec[1]*q.vec[1]) *qnrm;



	qrmat[0] = qmat[0];
	qrmat[1] = qmat[1];
	qrmat[2] = -qmat[2];
	
	qrmat[3] = qmat[3];
	qrmat[4] = qmat[4];
	qrmat[5] = -qmat[5];
	
	qrmat[6] = -qmat[6];
	qrmat[7] = -qmat[7];
	qrmat[8] = qmat[8];

	qmat[0] *= kx; qmat[3] *= kx; qmat[6] *= kx;
	qmat[1] *= ky; qmat[4] *= ky; qmat[7] *= ky;
	qmat[2] *= kz; qmat[5] *= kz; qmat[8] *= kz;

	memset(K,0,sizeof(Storage::real)*9);
	int i,j,k;
	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			for(k = 0; k < 3; k++)
				K[i*3+j] += qmat[i*3+k]*qrmat[k*3+j];
}

void cell_to_INMOST(struct grid * g, int cell, Cell r)
{
	Storage::integer mat = r->Integer(g->cell_material);
	Storage::real center[3];
	r->Centroid(center);
	if( initialize )
	{		
		double ar = 1.0/(1+(center[0]-0.5)*(center[0]-0.5)+(center[1]-0.5)*(center[1]-0.5)+(center[2]-0.5)*(center[2]-0.5));
		r->Real(pressure) = get_press_by_mat(mat)*ar;
		r->Real(saturation) = get_sat_by_mat(mat)*ar;
	}
	else
	{
		/*
		int find_cell = cell;
		while( g->cells[find_cell].data->empty() && g->cells[find_cell].parent != -1 ) find_cell = g->cells[find_cell].parent;
		if( !g->cells[find_cell].data->empty() )
		{
			if( g->cells[find_cell].data->find(mat) == g->cells[find_cell].data->end() )
			{
				//Search in cells around
				int neighbours[1<<(DIM-1)];
				bool done = false;
				for(int k = 0; k < 6 && !done; k++)
				{
					int ret = cellAround(g,find_cell,k,neighbours);
					for(int q = 0; q < ret; q++)
					{
						int find_cell2 = neighbours[q];
						while( g->cells[find_cell2].data->empty() && g->cells[find_cell2].parent != -1 ) find_cell2 = g->cells[find_cell2].parent;
						if( g->cells[find_cell2].data->find(mat) != g->cells[find_cell2].data->end() )
						{
							find_cell = find_cell2;
							done = true;
							break;
						}
					}
				}
				if( done == false ) throw -1;
			}
			r->Real(pressure) = (*g->cells[find_cell].data)[mat].second[pressure] / (*g->cells[find_cell].data)[mat].first;
			r->Real(saturation) = (*g->cells[find_cell].data)[mat].second[saturation] / (*g->cells[find_cell].data)[mat].first;
			//~ if( isnan(r->Real(pressure)) || isnan(r->Real(saturation) ) )
			//~ {
				//~ std::cout << "sat " << r->Real(saturation) << " have to be " << get_sat_by_mat(mat) << std::endl;
				//~ std::cout << g->cells[find_cell].data[mat].second[saturation] << " " << g->cells[find_cell].data[mat].first << std::endl;
				//~ std::cout << "pressure " << r->Real(pressure) << " have to be " << get_press_by_mat(mat) << std::endl;
				//~ std::cout << g->cells[find_cell].data[mat].second[pressure] << " " << g->cells[find_cell].data[mat].first << std::endl;
				//~ 
			//~ }
		}
		else throw -1;
		*/
	}
#if defined(INIT_PERM)
	fill_K(center,&r->RealArrayDF(g->Kvec)[0],&r->RealArrayDF(g->K)[0]);
#endif
}

void init_mesh(struct grid * g)
{
	pressure   = g->mesh->CreateTag("P", DATA_REAL, CELL,false,1);
	saturation = g->mesh->CreateTag("S", DATA_REAL, CELL,false,1);
	gas = g->mesh->CreateTag("Sg", DATA_REAL, CELL,false,1);
	
}

std::map<Tag,Storage::real> cell_small_unite(ElementArray<Cell> & unite)
{
	std::map<Tag,Storage::real> ret;
	Storage::real vol = 0;
	const int ntags = 2;
	Tag tags[ntags];
	tags[0] = pressure;
	tags[1] = saturation;
	for(int i = 0; i < unite.size(); i++)
	{
		double v = unite[i]->Volume();
		for(int j = 0; j < ntags; j++)
		{
			Storage::real val = unite[i]->Real(tags[j]);
			ret[tags[j]] += val*v;
			
		}
		vol += v;
	}
	for(std::map<Tag,Storage::real>::iterator it = ret.begin(); it != ret.end(); it++)
		it->second /= vol;
	return ret;
}

void onclose()
{
	gridDelete(&thegrid);
	DestroyObj();
}


#if defined( __GRAPHICS__)
int show_region = 0;
int width = 800, height = 600;
double maxS = 1,minS = 0,maxP = 100,minP = 0;


void motion(int nmx, int nmy) // Mouse
{
	mx = ((nmx/(double)(width))-0.5)*((double)width/(double)height)+0.5;
	my = (1.0 - nmy/(double)height);	
	wasmotion = true;
	
	
	//~ if( print_element )
	//~ {
		//~ double val;
		//~ double v[3] = {mx,my,0.5};
		//~ int c = gridSearch(&thegrid,v);
		//~ if( c != -1 )
		//~ {
			//~ val = thegrid.cells[c].mr->Real(saturation);
			//~ printf("saturation on octree %g on grid %g\n",((double *)thegrid.cells[c].cell_data)[0],val);
			//~ val = thegrid.cells[c].mr->Real(pressure);
			//~ printf("pressure on octree %g on grid %g\n",((double *)thegrid.cells[c].cell_data)[1],val);
		//~ }
	//~ }
	
	glutPostRedisplay();
}



class face2gl
{
	double dist;
	double c[4];
	std::vector<double> verts;
public:
	face2gl():verts() {dist = 0; memset(c,0,sizeof(double)*4);}
	face2gl(const face2gl & other) :verts(other.verts) {dist = other.dist; memcpy(c,other.c,sizeof(double)*4);}
	face2gl & operator =(face2gl const & other) { verts = other.verts; dist = other.dist; memmove(c,other.c,sizeof(double)*4); return *this;}
	~face2gl() {}
	void draw() const
	{
		glColor4dv(c); 
		glBegin(GL_POLYGON); 
		{
			for(unsigned k = 0; k < verts.size(); k+=3) 
				glVertex3dv(&verts[k]); 
		}
		glEnd();
		glColor3d(0,0,0);
		glBegin(GL_LINE_LOOP);
		{
			for(unsigned k = 0; k < verts.size(); k+=3) 
				glVertex3dv(&verts[k]); 
		}
		glEnd();
	}
	bool operator <(const face2gl & other) const {return dist < other.dist;}
	void set_color(double r, double g, double b, double a) {c[0] = r; c[1] = g; c[2] = b; c[3] = a;}
	void add_vert(double x, double y, double z) {unsigned s = verts.size(); verts.resize(s+3); verts[s] = x; verts[s+1] = y; verts[s+2] = z;}
	void add_vert(double v[3]) {verts.insert(verts.end(),v,v+3);}
	void set_center(double cnt[3], double cam[3])
	{
		dist = sqrt((cnt[0]-cam[0])*(cnt[0]-cam[0])+(cnt[1]-cam[1])*(cnt[1]-cam[1])+(cnt[2]-cam[2])*(cnt[2]-cam[2]));
	}
};



face2gl DrawFace(Element f, int mmat, double campos[3])
{
	double cnt[3];
	face2gl ret;
	f->Centroid(cnt);
	ret.set_center(cnt,campos);
	ElementArray<Node> nodes = f->getNodes();

	if( f->nbAdjElements(CELL) == 0 )
		ret.set_color(1,0,0,0.25);
	else if( mmat == -1 )
	{
		ret.set_color(0.2,0,1,0.25);
	}
	else
	{
		double r = mmat / (double)get_type.Size();
		ret.set_color(r,1-r,0.5,0.25);
	}
	
	for( ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
		ret.add_vert(&(kt->Coords()[0]));
	

	//~ glColor3f(0,0,0);
	//~ glBegin(GL_LINE_LOOP);
	//~ for( adjacent<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
		//~ glVertex3dv(&(kt->Coords()[0]));
	//~ glEnd();
	
	return ret;
}


void whereami(double * pos)
{
   // Get the viewing matrix
   GLdouble modelview[16],projection[16];
   GLint viewport[4] = {0,0,1,1};
   glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
   glGetDoublev(GL_PROJECTION_MATRIX, projection);
   
   GLdouble outx, outy, outz;  // Var's to save the answer in

   gluUnProject(0.5, 0.5, 0.,
               modelview, projection, viewport,
               &outx, &outy, &outz);

   // Return the result.
   pos[0] = outx;
   pos[1] = outy;
   pos[2] = outz;
}


void draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	
	glTranslated((gleft+gright)*0.5,(gbottom+gtop)*0.5,(gnear+gfar)*0.5);
	rotate();
	glTranslated(-(gleft+gright)*0.5,-(gbottom+gtop)*0.5,-(gnear+gfar)*0.5);
	if( draw_quantity == 0 )
	{
		double campos[3];
		whereami(campos);
		std::vector<face2gl> polygons;
		for(Mesh::iteratorFace it = thegrid.mesh->BeginFace(); it != thegrid.mesh->EndFace(); it++)
		{
			if( it->nbAdjElements(CELL) <= 1  )
			{
				int mmat = (!it->BackCell().isValid() ? -1 : it->BackCell()->Integer(thegrid.cell_material));
				polygons.push_back(DrawFace(it->self(),mmat,campos));
			}
		}
		std::sort(polygons.begin(),polygons.end());
		for(int q = polygons.size(); q > 0 ; q--) polygons[q-1].draw();
		if( draw_nodes )
		{
			glPointSize(5);
		
			glBegin(GL_POINTS);
			for(Mesh::iteratorNode f = thegrid.mesh->BeginNode(); f != thegrid.mesh->EndNode(); f++)
			{
				if( !f->GetMarker(thegrid.octree_node) )
				{
					if( f->IntegerArray(thegrid.materials).size() > 2 )
						glColor3f(0,0,1);
					else glColor3f(1,0,0);
					glVertex3dv(&f->Coords()[0]);
				}
			}
			glEnd();
			glPointSize(1);
		}
	}
	else if( draw_quantity == 1 )
	{
		get_type.Draw();
	}
	else if( draw_quantity == 2 )
		make_proj.Draw();
	
	glutSwapBuffers();
}


void draw2()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	
	glTranslated((gleft+gright)*0.5,(gbottom+gtop)*0.5,(gnear+gfar)*0.5);
	rotate();
	glTranslated(-(gleft+gright)*0.5,-(gbottom+gtop)*0.5,-(gnear+gfar)*0.5);
	
	
	if( draw_quantity < 3 )
	{
		if( draw_regions )
		{
			
			int m = 0;
			for(Mesh::iteratorCell it = thegrid.mesh->BeginCell(); it != thegrid.mesh->EndCell(); it++)
			{
				double value,t;
				
				if( draw_quantity == 1 ) //draw_pressure
				{
					value = it->Real(pressure);
					t = (value-minP)/(maxP-minP);
					glColor3f(t*0.65,0.65*(t < 0.5 ? t : 1.0-t),0.65*(1-t));
				}
				else if( draw_quantity == 2 )
				{
					value = it->Real(saturation);
					t = (value-minS)/(maxS-minS);
					glColor3f(t*0.65,0.65*(t < 0.5 ? t : 1.0-t),0.65*(1-t));
				}
				else if( draw_quantity == 0 ) //draw type
				{
					//Storage::real ret[3];
					//~ Geometry::Centroid(&*it,ret);
					//~ double alpha = get_type.Layer(ret)*1.0/get_type.Size();
					double alpha = it->Integer(thegrid.cell_material)*1.0/get_type.Size();
					glColor3f(alpha,1.0-alpha,0.5);
				}
			
				//if( m == show_region )
				{
					int i = 0;
					ElementArray<Face> faces = it->getFaces();
					for(ElementArray<Face>::iterator f = faces.begin(); f != faces.end(); f++) if( f->Boundary() )
					{
						ElementArray<Node> nodes = f->getNodes();
						glBegin(GL_POLYGON);
						for(ElementArray<Node>::iterator n = nodes.begin(); n != nodes.end(); n++)
							glVertex3dv(&n->RealArray(thegrid.mesh->CoordsTag())[0]);
						glEnd();
					}
				}
				m++;
			}
		}
		
		
		if( draw_edges )
		{
			glLineWidth(2);
			glColor3f(0.0f,0.0f,0.0f);
			{
				glBegin(GL_LINES);
				for(Mesh::iteratorEdge f = thegrid.mesh->BeginEdge(); f != thegrid.mesh->EndEdge(); f++)
				{
					ElementArray<Node> nodes = f->getNodes();
					glVertex3dv(&nodes[0].RealArray(thegrid.mesh->CoordsTag())[0]);
					glVertex3dv(&nodes[1].RealArray(thegrid.mesh->CoordsTag())[0]);
				}
				glEnd();
			}
			glLineWidth(1);
		}
		
		if( draw_faces )
		{
			for(Mesh::iteratorFace f = thegrid.mesh->BeginFace(); f != thegrid.mesh->EndFace(); f++)
			{
				ElementArray<Node> nodes = f->getNodes();
				glColor3f(0.0f,0.5f,0.0f);	
				glBegin(GL_POLYGON);
				for(ElementArray<Node>::iterator n = nodes.begin(); n != nodes.end(); n++)
					glVertex3dv(&n->RealArray(thegrid.mesh->CoordsTag())[0]);
				glEnd();
			}
		}
		
		if( draw_nodes )
		{
			glPointSize(5);
			
			glBegin(GL_POINTS);
			for(Mesh::iteratorNode f = thegrid.mesh->BeginNode(); f != thegrid.mesh->EndNode(); f++)
			{
				if( !f->GetMarker(thegrid.octree_node) )
				{
					if( f->IntegerArray(thegrid.materials).size() > 2 )
						glColor3f(0,0,1);
					else glColor3f(1,0,0);
					glVertex3dv(&f->Coords()[0]);
				}
			}
			glEnd();
			glPointSize(1);
		}
	}
	else if( draw_quantity == 3 )
	{
		get_type.Draw();
	}
	else if( draw_quantity == 4 )
		make_proj.Draw();
	
	
	
		
	
	
	glutSwapBuffers();

}


void reshape(int w, int h)
{
	/*
	double aspect = (double)w/(double)h;
	width = w;
	height = h;
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho(-50*aspect,50*aspect,-50,50,-5000,5000);
	glMatrixMode (GL_MODELVIEW);
	glViewport(0, 0, w, h);
	*/
	const double sc = 1.5;
	double aspect = (double)w/(double)h;
	double center[3] = { (gleft+gright)*0.5, (gbottom+gtop)*0.5, (gfar+gnear)*0.5};
	double side = std::max(std::max( gright-gleft, gtop-gbottom ), gfar-gnear)*0.5;
	width = w;
	height = h;
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho(center[0]-sc*side*zoom*aspect,center[0]+sc*zoom*side*aspect,
			center[1]-sc*side*zoom,center[1]+sc*zoom*side,
			center[2]-sc*side*100,center[2]+sc*side*100);
	glMatrixMode (GL_MODELVIEW);
	glViewport(0, 0, w, h);
}

int s = 0;


void keyboard(unsigned char key, int x, int y)
{
	(void) x, (void) y;
	if( key == 27 )
	{
		onclose();
		//~ glutLeaveMainLoop();
		exit(-1);
	}
	else if( key == ' ' ) 
	{
		std::cout << "mx: " << mx << "my: " << my << std::endl;
		gridAMR(&thegrid,recreate);
		glutPostRedisplay();
	}
	else if( key == 'd' || key == 'D' )
	{
		draw_quantity = (draw_quantity+1)%3;
		glutPostRedisplay();
	}
	else if( key == '+' || key == '=' )
	{
		/*
		int t = thegrid.mesh->NumberOfCells();
		show_region = (show_region+1+t)%t;
		printf("%d / %d \n",show_region,t);*/
		zoom /= 1.1;
		reshape(width,height);
		glutPostRedisplay();
	}
	else if( key == '-' || key == '_' )
	{
		/*
		int t = thegrid.mesh->NumberOfCells();
		show_region = (show_region-1+t)%t;
		printf("%d / %d \n",show_region,t);*/
		zoom *= 1.1;
		reshape(width,height);
		glutPostRedisplay();
	}
	else if( key == 'f' || key == 'F' )
	{
		thegrid.mesh->Save("out.gmv");
		thegrid.mesh->Save("out.vtk");
		thegrid.mesh->Save("out.pmf");
	}
	else if( key == 'A' || key == 'a' )
	{
		allow_amr = !allow_amr;
	}
	else if( key == 'S' || key == 's' )
	{
		draw_scheme = (draw_scheme+1)%3;
	}
	else if( key == 'e' || key == 'E' )
	{
		print_element = !print_element;
	}
	else if( key == 'r' || key == 'R' )
	{
		recreate = !recreate;
	}
	else if( key == 'c' || key == 'C' )
	{
		allow_coarse = !allow_coarse;
	}
	else if( key == 'v' || key == 'V' )
	{
		allow_refine = !allow_refine;
	}
	else if( key == 't' )
	{
		transparent = !transparent;
		if( !transparent ) 
		{
			glDisable(GL_BLEND);
			glLineWidth(2);
		}
		else 
		{
			glEnable(GL_BLEND);
			glLineWidth(1.0);
		}
		glutPostRedisplay();
		
	}
}




void idle(void)
{
	if( wasmotion && allow_amr == 1 && (Timer()-last_time) > 0.01)
	{
		last_time = Timer();
		gridAMR(&thegrid,recreate);
		glutPostRedisplay();
		wasmotion = false;
	}
}
#endif



void vert_to_INMOST(struct grid * g, int vert, Node v) 
{
	(void) g; (void) vert;
	Storage::real_array c = v->Coords();
	mat_ret_type ret = get_material_types(&c[0]);
	Storage::integer_array mats = v->IntegerArray(g->materials);
	mats.clear();
	mats.insert(mats.end(),ret.begin(),ret.end());
}


bool global_problem = false;
int global_test_number = -1;
int main(int argc, char ** argv)
{
	int i;
#if defined(MASAHIKO_MESH)
	int n[3] = {36,36,36};
#else
	int n[3] = {8,8,8};
#endif
	last_time = Timer();
	
	if( argc > 3 )
	{
		n[0] = atoi(argv[1]);
		n[1] = atoi(argv[2]);
		n[2] = atoi(argv[3]);
	}
	
#if defined(OBJ_READERS)
	make_proj.SetGridStep(n[2]);
	
	InitObj();
	{
		std::vector<std::string> layers;
        for(int i = 1; i <= 3; i++)
		{
				 std::stringstream name;
				//name << "Obj/bound" << i+1 << ".obj";
				 name << "oil_obj2/proj/layer" << i << ".obj";
				layers.push_back(name.str());
		}
		make_proj.ReadLayers(layers);
	}
	{
		
		std::vector<std::string> layers;
		for(int i = 1; i <= 1; i++)
		{
			 std::stringstream name;
				//name << "Obj/rt" << i+1 << ".obj";
				 name << "oil_obj2/mat/layer" << i << ".obj";
				layers.push_back(name.str());
		}
		 
		/*
		for(int i = 1; i <= 3; i++)
		{
				 std::stringstream name;
				//name << "Obj/rt" << i+1 << ".obj";
				 name << "oil_obj2/proj/layer" << i << ".obj";
				layers.push_back(name.str());
		}
		*/
		get_type.ReadLayers(layers);
	}
#endif
	
	
	thegrid.transformation = transformation;
	thegrid.cell_should_unite = cell_should_unite;
	thegrid.cell_should_split = cell_should_split;
	thegrid.cell_unite_data = cell_unite_data;
	thegrid.cell_split_data = cell_split_data;
	thegrid.cell_init_data = cell_init_data;
	thegrid.cell_destroy_data = cell_destroy_data;
	thegrid.cell_to_INMOST = cell_to_INMOST;
	thegrid.vert_init_data = default_vert_init_data;
	thegrid.vert_interpolate_data = default_vert_interpolate_data;
	thegrid.vert_destroy_data = default_vert_destroy_data;
	thegrid.vert_to_INMOST = vert_to_INMOST;
	thegrid.cell_small_unite = cell_small_unite;
	thegrid.init_mesh = init_mesh;
	thegrid.get_material_types = get_material_types;

	
	
	gridInit(&thegrid,n);
	
	
	
	//~ Tag boundary = thegrid.mesh->CreateTag("BOUNDARY",DATA_INTEGER,FACE,false,1);
	//~ for(Mesh::iteratorFace f = thegrid.mesh->BeginFace(); f != thegrid.mesh->EndFace(); f++)
		//~ if( Geometry::Boundary(&*f) )
			//~ f->Integer(boundary) = 1;

	initialize = 0;
		
	//mx = 0.675; my = 0.65;
	//mx = 0.4083333; my = 0.355;
	//gridAMR(&thegrid,recreate);
	//mx = 0.783333; my = 0.645;
	//gridAMR(&thegrid,recreate);
	/*
	//~ 
	thegrid.mesh->Save("out.pmf");
	thegrid.mesh->Save("out.gmv");
	gridDelete(&thegrid);
	DestroyObj();
	return 0;
	*/
	/*
	for(int i = 0; i < 500; i++)
	{
		global_test_number = i;
		srand(i*100);
		mx = rand()/(double)RAND_MAX; my = rand()/(double)RAND_MAX;
		std::cout << "i: " << i <<  " mx: " << mx << " my: " << my << std::endl;
		gridAMR(&thegrid,recreate);
		
		//~ if( global_problem || thegrid.mesh->GetTopologyError() )
		//~ {
			//~ std::stringstream sstrpmf, sstrgmv;
			//~ sstrpmf << "problem_grids/" << i << ".pmf";
			//~ thegrid.mesh->Save(sstrpmf.str());
			//~ sstrgmv << "problem_grids/" << i << ".gmv";
			//~ thegrid.mesh->Save(sstrgmv.str());
			//~ std::cout << "problem " << i << " files written " << std::endl;
			//~ global_problem = false;
			//~ thegrid.mesh->ClearTopologyError();
		//~ }
		
	}
	
	//thegrid.mesh->Save("out.pmf");
	//thegrid.mesh->Save("out.gmv");
	
	gridDelete(&thegrid);
	DestroyObj();
	return 0;
	*/
#if defined(INIT_PERM)
	for(Mesh::iteratorFace it = thegrid.mesh->BeginFace(); it != thegrid.mesh->EndFace(); ++it)
	{
		Storage::real cnt[3];
		it->Centroid(cnt);
		fill_K(cnt,&it->RealArrayDF(thegrid.Kvec)[0],&it->RealArrayDF(thegrid.K)[0]);
	}
#endif
	
#if defined(SPHERE_MESH)
	thegrid.mesh->BeginModification();
	int cnt = 0;
	for(Mesh::iteratorCell it = thegrid.mesh->BeginCell(); it != thegrid.mesh->EndCell(); ++it)
		if( it->Integer(thegrid.cell_material) < 0 )
		{
			it->Delete();
			cnt++;
		}

	std::cout << "deleted CELL " << cnt << std::endl;
	{
		ElementType etype = CELL;
		do
		{
			cnt = 0;
			etype = PrevElementType(etype);
			for(Mesh::iteratorElement it = thegrid.mesh->BeginElement(etype); it != thegrid.mesh->EndElement(); ++it)
				if( it->nbAdjElements(NextElementType(etype)) == 0 )
				{
					it->Delete();
					cnt++;
				}

			std::cout << "deleted " << ElementTypeName(etype) << " " << cnt << std::endl;
		}
		while(etype != NODE);
	}
	thegrid.mesh->Save("sphere.vtk");
	thegrid.mesh->SwapModification();
	thegrid.mesh->EndModification();
	onclose();
	return 0;
#endif


#if defined(MASAHIKO_MESH)

	gridAMR(&thegrid,recreate);
	//gridRefine(&thegrid);
	gridRecreateINMOST(&thegrid);
	{
		Storage::real x[3];
		Tag tensor_K = thegrid.mesh->CreateTag("PERM",DATA_REAL,CELL | FACE, FACE, 3);
		Tag poro = thegrid.mesh->CreateTag("PORO",DATA_REAL,CELL | FACE, FACE, 1);
		Tag aperture = thegrid.mesh->CreateTag("APERTURE",DATA_REAL,FACE,FACE,1);
		Tag mats = thegrid.mesh->GetTag("MATERIALS"); //face,edge,node
		Tag mat = thegrid.mesh->GetTag("MATERIAL");
		int fault_faces = 0;
		for(Mesh::iteratorFace f = thegrid.mesh->BeginFace(); f != thegrid.mesh->EndFace(); ++f)
		{
			Storage::real_array fmats = f->RealArray(mats);
			if( fmats.size() == 2 ) //fault
			{
				Storage::real_array K = f->RealArray(tensor_K);
				K[0] = K[1] = K[2] = 10000; // 10D perm
				f->Real(aperture) = 0.01;
				f->Real(poro) = 0.6;
				fault_faces++;
			}
		}
		std::cout << "fault faces: " << fault_faces << std::endl;
		for(Mesh::iteratorCell c = thegrid.mesh->BeginCell(); c != thegrid.mesh->EndCell(); ++c)
		{
			c->Centroid(x);
			Storage::real_array K = c->RealArray(tensor_K);
			if( x[2] <= 11.9/36.0 ) // bottom reservoir
			{
				c->Integer(mat) = 0;
				K[0] = (rand()*1.0)/static_cast<double>(RAND_MAX)*90+10;
				K[1] = (rand()*1.0)/static_cast<double>(RAND_MAX)*90+10;
				K[2] = (rand()*1.0)/static_cast<double>(RAND_MAX)*9+1;
				c->Real(poro) = sqrt((K[0]*K[0] + K[1]*K[1] + K[2]*K[2])/(100*100+100*100+10*10))*0.6;
			}
			else if( x[2] <= 12.0/36.0 ) //empty part between bottom and middle reservoir
			{
				c->Integer(mat) = 1;
				K[0] = K[1] = K[2] = 0.1;
				c->Real(poro) = 0.01;
			}
			else if( x[2] <= 23.9/36.0 ) //middle reservoir
			{
				c->Integer(mat) = 2;
				K[0] = (rand()*1.0)/static_cast<double>(RAND_MAX)*120+30;
				K[1] = (rand()*1.0)/static_cast<double>(RAND_MAX)*170+30;
				K[2] = (rand()*1.0)/static_cast<double>(RAND_MAX)*10+5;
				c->Real(poro) = sqrt((K[0]*K[0] + K[1]*K[1] + K[2]*K[2])/(150*150+200*200+15*15))*0.6;
			}
			else if( x[2] <= 24.0/36.0 ) //empty part between middle and top reservoir
			{
				c->Integer(mat) = 3;
				K[0] = K[1] = K[2] = 0.1;
				c->Real(poro) = 0.01;
			}
			else //top reservoir
			{
				c->Integer(mat) = 4;
				K[0] = (rand()*1.0)/static_cast<double>(RAND_MAX)*60+20;
				K[1] = (rand()*1.0)/static_cast<double>(RAND_MAX)*40+20;
				K[2] = (rand()*1.0)/static_cast<double>(RAND_MAX)*6+3;
				c->Real(poro) = sqrt((K[0]*K[0] + K[1]*K[1] + K[2]*K[2])/(80*80+60*60+9*9))*0.6;
			}
		}
		thegrid.mesh->Save("masahiko.xml");
		thegrid.mesh->Save("masahiko.vtk");
		onclose();
		return 0;
	}
#endif



#if defined(__GRAPHICS__)

	for(Mesh::iteratorNode n = thegrid.mesh->BeginNode(); n != thegrid.mesh->EndNode(); n++)
	{
		Storage::real_array c = n->Coords();
		if( c[0] > gright ) gright = c[0];
		if( c[0] < gleft ) gleft = c[0];
		if( c[1] > gtop ) gtop = c[1];
		if( c[1] < gbottom ) gbottom = c[1];
		if( c[2] > gfar ) gfar = c[2];
		if( c[2] < gnear ) gnear = c[2];
	}


	quatinit();
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutInitWindowPosition (100, 100);
	glutCreateWindow("Graph");

	
	glDepthFunc(GL_LEQUAL);
	glClearDepth(1.f);
	glEnable(GL_DEPTH_TEST);
	
	//glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	//~ glEnable(GL_POLYGON_SMOOTH);
	//~ glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);
	//~ glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	//~ glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
	//~ glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
	
	glClearColor (1.0f, 1.0f, 1.0f, 1.f);
	glutDisplayFunc(draw);
	glutReshapeFunc(reshape);
	glLineWidth(2);
	
	glutKeyboardFunc(keyboard);
	glutMotionFunc(rotate_clickmotion);
	//~ glutWMCloseFunc(onclose);
	glutMouseFunc(rotate_click);
	//glutPassiveMotionFunc(rotate_motion);
	glutPassiveMotionFunc(motion);
	glutIdleFunc(idle);
	

	glutPostRedisplay();
	glutMainLoop();
#endif
	//~ std::cout << "Hello!" << std::endl;
	onclose();
}
