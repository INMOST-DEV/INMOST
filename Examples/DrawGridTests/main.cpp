//g++ main.cpp rotate.cpp -L/usr/X11R6/lib -lX11 -lXi -lXmu -lGL -lglut -lGLU ../../inmost.a -O5
// press space - explode mesh to see connection 
//#define OCTREECUTCELL_DEBUG
#define DISCR_DEBUG

#include "inmost.h"
#include "my_glut.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include "rotate.h"
#include <stdarg.h>

#if defined(USE_FP64)
#define glVertex3v glVertex3dv
#define glVertex2v glVertex2dv
#define glRasterPos3v glRasterPos3dv
#else
#define glVertex3v glVertex3fv
#define glVertex2v glVertex2fv
#define glRasterPos3v glRasterPos3fv
#endif


using namespace INMOST;
Mesh * mesh;
int drawgrid = 0;
int interactive = 0;
double zoom = 1;
int width = 800, height = 800;
double sleft = 1e20, sright = -1e20, sbottom = 1e20, stop = -1e20, sfar = -1e20, snear = 1e20;
Tag mat,colort; int maxmat = -1, maxcolor = 0;
Mesh::GeomParam table;
double reset_timer = Timer();
double shift[3] = {0,0,0};
int boundary = 2;
int color = 0;
int edges = 1;
bool perspective = false;
int text = 4;
int display_orphans = 1;
int badfaces = 1;
ElementSet problem_set, orphan_edges, orphan_faces;


#if defined(DISCR_DEBUG)
Tag avg_coord, tensor;

static Storage::real dot_prod(Storage::real x[3], Storage::real y[3])  { return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]; }
static void tensor_prod3(Storage::real * K, Storage::real v[3], Storage::real out[3])
{
	out[0] = v[0] * K[0] + v[1] * K[1] + v[2] * K[2];
	out[1] = v[0] * K[3] + v[1] * K[4] + v[2] * K[5];
	out[2] = v[0] * K[6] + v[1] * K[7] + v[2] * K[8];
}

void FindHarmonicPoint(Face & fKL, Cell & cK, Cell & cL, Tag tensor, Storage::real xK[3], Storage::real xL[3], Storage::real y[3], Storage::real & coef)
{
	Storage::real yK[3], yL[3], xKL[3], nKL[3], dK, dL, lK, lL, lKs[3], lLs[3], D;
	Storage::real coefK, coefL, coefQ, coefDiv;
	Storage::real_array KK = cK->RealArray(tensor), KL = cL->RealArray(tensor);
	fKL->OrientedUnitNormal(cK,nKL);
	fKL->Centroid(xKL);
	tensor_prod3(&KK[0],nKL,lKs);
	tensor_prod3(&KL[0],nKL,lLs);
	lK = dot_prod(nKL,lKs);
	lL = dot_prod(nKL,lLs);
	lKs[0] -= nKL[0]*lK;
	lKs[1] -= nKL[1]*lK;
	lKs[2] -= nKL[2]*lK;
	lLs[0] -= nKL[0]*lL;
	lLs[1] -= nKL[1]*lL;
	lLs[2] -= nKL[2]*lL;
	D = -dot_prod(xKL,nKL);
	dK = fabs(dot_prod(xK,nKL)+D);
	dL = fabs(dot_prod(xL,nKL)+D);
	yK[0] = xK[0] + dK*nKL[0];
	yK[1] = xK[1] + dK*nKL[1];
	yK[2] = xK[2] + dK*nKL[2];
	yL[0] = xL[0] - dL*nKL[0];
	yL[1] = xL[1] - dL*nKL[1];
	yL[2] = xL[2] - dL*nKL[2];
	coefDiv = (lL*dK+lK*dL);
	coefL = lL*dK/coefDiv;
	coefK = lK*dL/coefDiv;
	coefQ = dK*dL/coefDiv;
	y[0] = yK[0]*coefK + yL[0]*coefL + (lKs[0]-lLs[0])*coefQ;
	y[1] = yK[1]*coefK + yL[1]*coefL + (lKs[1]-lLs[1])*coefQ;
	y[2] = yK[2]*coefK + yL[2]*coefL + (lKs[2]-lLs[2])*coefQ;
	coef = coefL;
}


void FindBoundaryPoint(Face & fK, Cell & cK, Tag tensor, Storage::real xK[3], Storage::real y[3])
{
	Storage::real nK[3], lKs[3], lK, xfK[3], t;
	Storage::real_array KK = cK->RealArray(tensor);
	fK->Centroid(xfK);
	fK->OrientedUnitNormal(cK,nK);
	tensor_prod3(&KK[0],nK,lKs);
	lK = dot_prod(nK,lKs);
	t = (dot_prod(nK,xfK) - dot_prod(nK,xK))/lK;
	y[0] = xK[0] + t*lKs[0];
	y[1] = xK[1] + t*lKs[1];
	y[2] = xK[2] + t*lKs[2];
}
#endif


void printtext(const char * fmt, ... )
{

	unsigned int i;
	char stext[1024];
	va_list ap;
	if ( fmt == NULL ) return;
	va_start(ap,fmt);
	vsprintf(stext,fmt,ap);
	va_end(ap);
	for(i=0;i<strlen(stext);i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,
							stext[i]);
	}
}

void set_matrix3d()
{
	double aspect = (double)width/(double)height;
	double side = std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
	//double center[3] = { (sleft+sright)*0.5, (sbottom+stop)*0.5, (sfar+snear)*0.5};
	const double sc = 2;
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	//~ glOrtho(center[0]-sc*side*zoom*aspect,center[0]+sc*zoom*side*aspect,
		//~ center[1]-sc*side*zoom,center[1]+sc*zoom*side,
		//~ center[2]-sc*side*100,center[2]+sc*side*100);
	if( !perspective )
	{
		glOrtho(-sc*side*zoom*aspect,sc*side*zoom*aspect,
				-sc*side*zoom,sc*side*zoom,
				-sc*side*100,sc*side*100);
	}
	else
		gluPerspective(60.0, aspect, 0.00001, 1000.0);
	glMatrixMode (GL_MODELVIEW);
}

void set_matrix2d()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1,1,-1,1,-1,1);
	glMatrixMode(GL_MODELVIEW);
}

void reshape(int w, int h)
{
	width = w;
	height = h;
	set_matrix3d();
	glViewport(0, 0, w, h);
}

class Input
{
public:
	enum InputType {Double,Integer};
private:
	std::string str;
	void * input_link;
	InputType type;
	bool done;


public:
	Input(int * val) {input_link = val; type = Integer; done = false; str = "";}
	Input(double * val) {input_link = val; type = Double; done = false; str = "";}
	Input(void * link, InputType type) : str(""), input_link(link), type(type), done(false) {}
	Input(const Input & other): str(other.str), input_link(other.input_link), type(other.type), done(other.done) {}
	Input & operator =(Input const & other) {input_link = other.input_link; str = other.str; type = other.type; done = other.done; return *this;}
	~Input() {}
	void KeyPress(char c)
	{
		if( c == 13 )
		{

			done = true;
			if( type == Double ) *((double *)input_link) = atof(str.c_str());
			else if( type == Integer ) *((int *)input_link) = atoi(str.c_str());
			glutPostRedisplay();
		}
		else if( c == 8 )
		{
			if( !str.empty() ) str.erase(str.size()-1);
			glutPostRedisplay();
		}
		else if( (c >= '0' && c <= '9') || c=='+' || c=='-' || c=='.' || c=='e' || c == 'E' )
		{
			str += c;
			glutPostRedisplay();
		}
		else if( c == 27 )
		{

			done = true;
			glutPostRedisplay();
		}
	}
	bool Done() {return done;}
	void Draw()
	{
		float h = 24.0f/(float)height;

		glColor3f(1,1,1);
		glBegin(GL_QUADS);
		glVertex3f(-0.99,-0.99,1);
		glVertex3f(-0.99,-0.99+h,1);
		glVertex3f(0.99,-0.99+h,1);
		glVertex3f(0.99,-0.99,1);
		glEnd();
		glColor3f(0,0,0);
		glBegin(GL_LINE_LOOP);
		glVertex3f(-0.99,-0.99,1);
		glVertex3f(-0.99,-0.99+h,1);
		glVertex3f(0.99,-0.99+h,1);
		glVertex3f(0.99,-0.99,1);
		glEnd();

		glColor4f(0,0,0,1);
		glRasterPos2d(-0.985,-0.985);
		//printtext(str.c_str());
		printtext("input number (%s): %s",type == Integer ? "integer":"double", str.c_str());
	}
} * CommonInput = NULL;





void keyboard2(unsigned char key, int x, int y)
{
	(void)x,(void)y;
	if( key == '=' || key == '+' ||  key == '_' || key == '-' || key == 'w' || key == 's' || key == 'a' || key == 'd' || key == 'r' || key == 'p' || key == 'z')
	{
		interactive = false;
		glutPostRedisplay();
	}

}

int matfilter = 0;
int show_cell = -1, show_face = -1;
int parentfilter = -1;
std::vector<int> show_cells, show_faces;

void keyboard(unsigned char key, int x, int y)
{
	(void)x,(void)y;
	if( CommonInput != NULL )
	{
		CommonInput->KeyPress(key);
		return;
	}
	if( key == 't' )
	{
		text = (text+1)%6;
		glutPostRedisplay();
	}
	else if( key == 'b' )
	{
		boundary = (boundary+1)%3;
		glutPostRedisplay();
	}
	else if( key == 'p' )
	{
		perspective = !perspective;
		reshape(width,height);
		glutPostRedisplay();
		interactive = true;
	}
	else if( key == 'w' )
	{
		double shiftmod[3] = {0,0,0};
		shiftmod[1] -= 0.03f*expf(zoom-1)*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
		rotatevector((double*)shiftmod);
		shift[0] += shiftmod[0];
		shift[1] += shiftmod[1];
		shift[2] += shiftmod[2];
		glutPostRedisplay();
		interactive = true;
	}
	else if( key == 's' )
	{
		double shiftmod[3] = {0,0,0};
		shiftmod[1] += 0.03f*expf(zoom-1)*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
		rotatevector((double*)shiftmod);
		shift[0] += shiftmod[0];
		shift[1] += shiftmod[1];
		shift[2] += shiftmod[2];
		glutPostRedisplay();
		interactive = true;
	}
	else if( key == 'a' )
	{
		double shiftmod[3] = {0,0,0};
		shiftmod[0] += 0.03f*expf(zoom-1)*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
		rotatevector((double*)shiftmod);
		shift[0] += shiftmod[0];
		shift[1] += shiftmod[1];
		shift[2] += shiftmod[2];
		glutPostRedisplay();
		interactive = true;
	}
	else if( key == 'd' )
	{
		double shiftmod[3] = {0,0,0};
		shiftmod[0] -= 0.03f*expf(zoom-1)*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
		rotatevector((double*)shiftmod);
		shift[0] += shiftmod[0];
		shift[1] += shiftmod[1];
		shift[2] += shiftmod[2];
		glutPostRedisplay();
		interactive = true;
	}
	else if( key == 'r' )
	{
		shift[0] = 0.0f;
		shift[1] = 0.0f;
		shift[2] = 0.0f;
		zoom = 1;
		quatinit();
		glutPostRedisplay();
		interactive = true;
	}
	else if( key == 'e' )
	{
		edges = !edges;
		glutPostRedisplay();
	}
	else if( key == 27 )
	{
		//delete mesh;
		exit(-1);
	}
	else if( key == '=' || key == '+')
	{
		zoom /= 1.1;
		reshape(width,height);
		glutPostRedisplay();
		interactive = true;
		//reset_timer = Timer();
	}
	else if( key == '_' || key == '-')
	{
		zoom *= 1.1;
		reshape(width,height);
		glutPostRedisplay();
		interactive = true;
		//reset_timer = Timer();
	}
	else if( key == ' ')
	{
		show_cells.clear();
		show_faces.clear();
		//mesh->Save("out.vtk");
		//drawgrid = (drawgrid+1)%2;
		//glutPostRedisplay();
	}
	else if( key == 'q' )
	{
		mesh->PrepareGeometricData(table);
	}
	else if( key == 'w' )
	{
		mesh->RemoveGeometricData(table);
	}
	else if( key-'0' >= 0 && key-'9' <= 0 )
	{
		matfilter = key-'0';
		glutPostRedisplay();
	}
	else if( key == 'm' )
	{
		if( parentfilter == -1 )
		{
			if( CommonInput == NULL ) CommonInput = new Input(&parentfilter);
		}
		else parentfilter = -1;
		glutPostRedisplay();
	}
	else if( key == 'c' )
	{
		if( CommonInput == NULL ) CommonInput = new Input(&show_cell);
		glutPostRedisplay();
	}
	else if( key == 'f' )
	{
		if( CommonInput == NULL ) CommonInput = new Input(&show_face);
		glutPostRedisplay();
	}
	else if( key == 'z' )
	{
		if( !colort.isValid() )
		{
			Storage::integer visited = 0;
			colort = mesh->CreateTag("COLOR",DATA_INTEGER,CELL,NONE,1);
			std::vector<Cell> queue;
			ElementArray<Element> around;
			queue.reserve(256);
			while(visited != mesh->NumberOfCells() )
			{
				maxcolor++;
				std::cout << "new part " << maxcolor << std::endl;
				for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
				{
					if( it->Integer(colort) == 0 )
					{
						queue.push_back(it->self());
						visited++;
						queue.back()->Integer(colort) = maxcolor;
						break;
					}
				}
				while(!queue.empty())
				{
					around = queue.back()->BridgeAdjacencies(FACE,CELL);

					queue.pop_back();
					for(ElementArray<Element>::iterator it = around.begin(); it != around.end(); ++it)
						if( it->Integer(colort) == 0 )
						{
							queue.push_back(it->getAsCell());
							visited++;
							queue.back()->Integer(colort) = maxcolor;
						}
				}
			}
		}
		color = !color;
	}
	else if( key == 'o' )
	{
		display_orphans = !display_orphans;
		glutPostRedisplay();
	}

	else if( key == 'l' )
	{
		badfaces = !badfaces;
		glutPostRedisplay();
	}
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
		if( boundary == 2 ) glColor3dv(c);
		else glColor4dv(c);
		glBegin(GL_POLYGON);
		for(unsigned k = 0; k < verts.size(); k+=3)
			glVertex3dv(&verts[k]);
		glEnd();
	}
	bool operator <(const face2gl & other) const {return dist < other.dist;}
	void set_color(double r, double g, double b, double a) {c[0] = r; c[1] = g; c[2] = b; c[3] = a;}
	void add_vert(double x, double y, double z) {unsigned s = verts.size(); verts.resize(s+3); verts[s] = x; verts[s+1] = y; verts[s+2] = z;}
	void add_vert(double v[3]) {verts.insert(verts.end(),v,v+3);}
	void add_vert(float v[3]) {verts.insert(verts.end(),v,v+3);}
	void set_center(Storage::real cnt[3], double cam[3])
	{
		dist = sqrt((cnt[0]-cam[0])*(cnt[0]-cam[0])+(cnt[1]-cam[1])*(cnt[1]-cam[1])+(cnt[2]-cam[2])*(cnt[2]-cam[2]));
	}
};



face2gl DrawFace(Element & f, int mmat, double campos[3])
{
	Storage::real cnt[3];
	face2gl ret;
	f->Centroid(cnt);
	ret.set_center(cnt,campos);
	ElementArray<Node> nodes = f->getNodes();

	if( f->nbAdjElements(CELL) == 0 )
		ret.set_color(1,0,0,0.25);
	else if( mmat == 0 )
	{
		ret.set_color(0.2,0,1,0.25);
	}
	else
	{
		double r = color ? mmat / (double)maxcolor : mmat / (double) maxmat;
		ret.set_color(1-r,r,1-r,0.25);
	}
	double v[3] = {0,0,0};
	for( ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
	{
		for(int k = 0; k < f->GetMeshLink()->GetDimensions(); ++k)
			v[k] = kt->Coords()[k];
		ret.add_vert(v);
	}


	if( edges )
	{
		//glLineWidth(2.0);
		if( boundary == 1 ) glColor3f(0,0,1); else glColor3f(0,0,0);
		glBegin(GL_LINE_LOOP);
		for( ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
		{
			if( f->GetMeshLink()->GetDimensions() == 2 )
				glVertex2v(&(kt->Coords()[0]));
			else 
				glVertex3v(&(kt->Coords()[0]));
		}
		glEnd();
		//glLineWidth(1.0);
	}
	return ret;
}

Tag problem_tag, mats_tag,parent_tag,cell_material_tag;

void fill_mat_str(Storage::integer id, Storage::integer_array mats, std::stringstream & str)
{
	str.str(std::string());
	str.clear();
	str << id;
	for(unsigned q = 0; q < mats.size(); q++)
		str << "," << mats[q]+1;
}


void fill_str(Storage::integer id, std::stringstream & str)
{
	str.str(std::string());
	str.clear();
	str << id;
}

void whereami(double & cx, double & cy, double & cz)
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
   cx = outx;
   cy = outy;
   cz = outz;
}

void draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	set_matrix3d();

	if( perspective )
		glTranslated(0,0,-zoom*2*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5);
	//glTranslated((sleft+sright)*0.5,(sbottom+stop)*0.5,(snear+sfar)*0.5);
	rotate();

	//glPointSize(5);
	//glColor3f(0,1,1);
	glLineWidth(3.0);
	glBegin(GL_LINES);
	glColor3f(1,0,0);
	glVertex3d(0,0,0);
	glVertex3d(zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5*0.1,0,0);
	glColor3f(0,1,0);
	glVertex3d(0,0,0);
	glVertex3d(0,zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5*0.1,0);
	glColor3f(0,0,1);
	glVertex3d(0,0,0);
	glVertex3d(0,0,zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5*0.1);
	glEnd();
	glLineWidth(1.0);
	//glPointSize(1);

	//glTranslated(-(sleft+sright)*0.5+shift[0],-(sbottom+stop)*0.5 + shift[1],-(snear+sfar)*0.5 + shift[2]);

	glTranslated(shift[0],shift[1],shift[2]);

	double campos[3] = {0.5,0.5,0};
	//~ double matrix[16], inverse[16];
	//~ glGetDoublev(GL_MODELVIEW_MATRIX,matrix);
	//~ inverse_v1(inverse,matrix);
	//~ campos[0] = inverse[12];
	//~ campos[1] = inverse[13];
	//~ campos[2] = inverse[14];

	//~ campos[0] = inverse[3];
	//~ campos[1] = inverse[7];
	//~ campos[2] = inverse[11];

	whereami(campos[0],campos[1],campos[2]);


	//glTranslated((l+r)*0.5,(b+t)*0.5,(near+far)*0.5);
	Storage::integer pacef = std::max<Storage::integer>(1,mesh->FaceLastLocalID()/10000);

	std::vector<face2gl> polygons;

	if( badfaces )
	{
		for(INMOST_DATA_INTEGER_TYPE it = 0; it < mesh->FaceLastLocalID(); it += (interactive ? pacef : 1)) if( mesh->isValidFace(it) )
		{
			Face f = mesh->FaceByLocalID(it);
			if( !f->CheckNormalOrientation() )
				polygons.push_back(DrawFace(f,1,campos));
		}
	}

	std::stringstream sstr;
	//~ char str[2048];

#if defined(DISCR_DEBUG)
	if(avg_coord.isValid())
	{
		glLineWidth(2.0);
		glColor3f(0,1,0);
		glBegin(GL_LINES);
		for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
		{
			Storage::real cnt[3];
			it->Centroid(cnt);
			ElementArray<Face> faces = it->getFaces();
			for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt) if(!jt->Boundary())
			{
				Storage::real_array cc = jt->RealArray(avg_coord);
				glVertex3v(cnt);
				glVertex3v(&cc[0]);
			}
		}
		glEnd();
		glLineWidth(1.0);
		glPointSize(4.0);
		glColor3f(1,0,0);
		glBegin(GL_POINTS);
		for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); ++it) if(!it->Boundary() )
		{
			Storage::real_array cc = it->RealArray(avg_coord);
			glVertex3v(&cc[0]);
		}
		glEnd();
		glPointSize(1.0);
	}

	if(tensor.isValid())
	{
		glColor3f(0,0,0);
		glBegin(GL_LINES);
		Storage::real h = 0, nh = 0;
		for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); ++it) if(!it->Boundary())
		{
			h += (it->BackCell()->Volume() + it->FrontCell()->Volume());
			nh += it->Area();
		}
		h /= nh;
		for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); ++it) if(!it->Boundary() )
		{
			Storage::real nrm[3], scale, f1[3], f2[3], l, cnt[3], v[3];
			Cell c1 = it->BackCell(), c2 = it->FrontCell();
			if( c1.isValid() && c2.isValid() )
			{
				it->Centroid(cnt);
				it->OrientedUnitNormal(c1,nrm);
				scale = h;
				Storage::real_array K1 = c1->RealArray(tensor), K2 = c2->RealArray(tensor);
				f1[0] = K1[0]*nrm[0] + K1[1]*nrm[1] + K1[2]*nrm[2];
				f1[1] = K1[3]*nrm[0] + K1[4]*nrm[1] + K1[5]*nrm[2];
				f1[2] = K1[6]*nrm[0] + K1[7]*nrm[1] + K1[8]*nrm[2];

				l = 0.15*scale/sqrt(f1[0]*f1[0]+f1[1]*f1[1]+f1[2]*f1[2]);

				f1[0] *= l;
				f1[1] *= l;
				f1[2] *= l;

				v[0] = f1[0] + cnt[0];
				v[1] = f1[1] + cnt[1];
				v[2] = f1[2] + cnt[2];

				glVertex3v(cnt);
				glVertex3v(v);

				v[0] = -f1[0] + cnt[0];
				v[1] = -f1[1] + cnt[1];
				v[2] = -f1[2] + cnt[2];

				glVertex3v(cnt);
				glVertex3v(v);

				f2[0] = K2[0]*nrm[0] + K2[1]*nrm[1] + K2[2]*nrm[2];
				f2[1] = K2[3]*nrm[0] + K2[4]*nrm[1] + K2[5]*nrm[2];
				f2[2] = K2[6]*nrm[0] + K2[7]*nrm[1] + K2[8]*nrm[2];



				l = 0.15*scale/sqrt(f2[0]*f2[0]+f2[1]*f2[1]+f2[2]*f2[2]);

				f2[0] *= l;
				f2[1] *= l;
				f2[2] *= l;

				v[0] = f2[0] + cnt[0];
				v[1] = f2[1] + cnt[1];
				v[2] = f2[2] + cnt[2];

				glVertex3v(cnt);
				glVertex3v(v);

				v[0] = -f2[0] + cnt[0];
				v[1] = -f2[1] + cnt[1];
				v[2] = -f2[2] + cnt[2];

				glVertex3v(cnt);
				glVertex3v(v);
			}
		}

		glEnd();

	}
#endif


	//glColor3f(0,0,0);

	if( !show_cells.empty() )
	{
		for(size_t qq = 0; qq < show_cells.size(); qq++)
		{
			if( show_cells[qq] < 0 || show_cells[qq] >= mesh->CellLastLocalID() ) continue;
			Element it = mesh->ElementByLocalID(CELL,show_cells[qq]);
			if( it.isValid() )
			{
				Storage::integer_array mats;
				ElementArray<Edge> edges = it->getEdges();
				Storage::real cnt[3];
				for(unsigned k = 0; k < edges.size(); k++)
				{
					edges[k].Centroid(cnt);
					if( mats_tag.isValid() )mats = edges[k].IntegerArray(mats_tag);

					ElementArray<Node> nodes = edges[k].getNodes();

					if( matfilter == 0 || (mats_tag.isValid() && std::binary_search(mats.begin(),mats.end(),matfilter-1)) )
					{
						if( mats_tag.isValid() ) fill_mat_str(edges[k].LocalID(), mats,sstr);
						else fill_str(edges[k].LocalID(),sstr);
						glColor3f(0,0,1);
						glRasterPos3v(cnt);
						if( text == 2 || text == 5 ) if(text) printtext(sstr.str().c_str());
						glColor3f(0,0,0);
						glBegin(GL_LINES);
						glVertex3v(&nodes[0].Coords()[0]);
						glVertex3v(&nodes[1].Coords()[0]);
						glEnd();
					}

					glColor3f(1,0,0);
					if(mats_tag.isValid())mats = nodes[0].IntegerArray(mats_tag);

					if( matfilter == 0 || (mats_tag.isValid() &&std::binary_search(mats.begin(),mats.end(),matfilter-1)) )
					{

						if( mats_tag.isValid() ) fill_mat_str(nodes[0].LocalID(), mats,sstr);
						else fill_str(nodes[0].LocalID(),sstr);
						glRasterPos3v(&nodes[0].Coords()[0]);
						if( text == 1 || text == 5 ) if(text) printtext(sstr.str().c_str());
					}

					if(mats_tag.isValid()) mats = nodes[1].IntegerArray(mats_tag);

					if( matfilter == 0 || (mats_tag.isValid() &&std::binary_search(mats.begin(),mats.end(),matfilter-1)) )
					{
						if( mats_tag.isValid() ) fill_mat_str(nodes[1].LocalID(),mats,sstr);
						else fill_str(nodes[1].LocalID(),sstr);
						glRasterPos3v(&nodes[1].Coords()[0]);
						if( text == 1 || text == 5 ) if(text) printtext(sstr.str().c_str());
					}
				}

				ElementArray<Face> faces = it->getFaces();

				for(unsigned k = 0; k < faces.size(); k++)
				{
					if(mats_tag.isValid())mats = faces[k].IntegerArray(mats_tag);
					if( matfilter == 0 || (mats_tag.isValid() &&std::binary_search(mats.begin(),mats.end(),matfilter-1)) )
					{
						if( faces[k].Boundary() )
							glColor3f(0.15,0.45,0);
						else glColor3f(0,1,0);
						faces[k].Centroid(cnt);
						if( mats_tag.isValid() ) fill_mat_str(faces[k].LocalID(),mats,sstr);
						else fill_str(faces[k].LocalID(), sstr);
						glRasterPos3v(cnt);
						if( text == 3 || text == 5 ) if(text) printtext(sstr.str().c_str());

						//~ glColor3f(0,0,0);
						//~ Storage::real nrm[3];
						//~ faces[k].OrientedNormal(it->getAsCell(),nrm);
						//~ glLineWidth(2.0);
						//~ glBegin(GL_LINES);
						//~ glVertex3v(cnt);
						//~ glVertex3d(cnt[0]+nrm[0]*20,cnt[1]+nrm[1]*20,cnt[2]+nrm[2]*20);
						//~ glEnd();
						//~ glLineWidth(1.0);
					}
					//polygons.push_back(DrawFace(&faces[k],0,campos));
				}


				if( matfilter == 0 || (mat.isValid() && it->Integer(mat) == matfilter-1) )
				{
					glColor3f(0,0,0);
					it->Centroid(cnt);
					fill_str(it->LocalID(),sstr);
					if(parent_tag.isValid()) sstr << ", " << it->Integer(parent_tag);
					if(cell_material_tag.isValid()) sstr << ", " << it->Integer(cell_material_tag);
					if( problem_tag.isValid() ) sstr << ", " << it->Integer(problem_tag);
					glRasterPos3v(cnt);
					if( text == 4 || text == 5 ) if(text) printtext(sstr.str().c_str());

					glPointSize(3);
					glColor3f(1,0,0);
					glBegin(GL_POINTS);
					glVertex3v(cnt);
					glEnd();
					glPointSize(1);
				}
			}
		}
	}
	std::vector<face2gl> inside_face;
	if( !show_faces.empty() )
	{
		for(size_t qq = 0; qq < show_faces.size(); qq++)
		{
			Element it = InvalidElement();
			if( show_faces[qq] < 0 || show_faces[qq] >= mesh->LastLocalID(FACE) ) continue;

			it = mesh->ElementByLocalID(FACE,show_faces[qq]);
			if( it.isValid() )
			{
				Storage::integer_array mats;
				Storage::real cnt[3], cntf[3];
				ElementArray<Edge> edges = it->getEdges();
				it->Centroid(cntf);
				for(unsigned k = 0; k < edges.size(); k++)
				//if( edges[k].nbAdjElements(CELL) == 0 )
				{
					if(mats_tag.isValid() ) mats = edges[k].IntegerArray(mats_tag);

					ElementArray<Node> nodes = edges[k].getNodes();

					if( matfilter == 0 || (mats_tag.isValid() &&std::binary_search(mats.begin(),mats.end(),matfilter-1)) )
					{
						edges[k].Centroid(cnt);
						if( mats_tag.isValid() ) fill_mat_str(edges[k].LocalID(), mats,sstr);
						else fill_str(edges[k].LocalID(),sstr);
						glColor3f(0,0,1);
						glRasterPos3v(cnt);
						if( text == 2 || text == 5 ) if(text) printtext(sstr.str().c_str());

						glColor3f(0,0,0);
						glBegin(GL_LINES);
						glVertex3v(&nodes[0].Coords()[0]);
						glVertex3v(&nodes[1].Coords()[0]);
						glEnd();

						face2gl f;
						f.set_color(0.6,0.6,0.6,0.1);
						f.add_vert(cntf);
						f.add_vert(&nodes[0].Coords()[0]);
						f.add_vert(&nodes[1].Coords()[0]);
						f.set_center(cntf,campos);
						inside_face.push_back(f);
					}

					glColor3f(1,0,0);
					 if( mats_tag.isValid() ) mats = nodes[0].IntegerArray(mats_tag);
					if( matfilter == 0 || (mats_tag.isValid() &&std::binary_search(mats.begin(),mats.end(),matfilter-1)) )
					{
						if( mats_tag.isValid() ) fill_mat_str(nodes[0].LocalID(), mats,sstr);
						else fill_str(nodes[0].LocalID(),sstr);
						glRasterPos3v(&nodes[0].Coords()[0]);
						if( text == 1 || text == 5 ) if(text) printtext(sstr.str().c_str());
					}
					if( mats_tag.isValid() ) mats = nodes[1].IntegerArray(mats_tag);
					if( matfilter == 0 || (mats_tag.isValid() &&std::binary_search(mats.begin(),mats.end(),matfilter-1)) )
					{

						if( mats_tag.isValid() ) fill_mat_str(nodes[1].LocalID(),mats,sstr);
						else fill_str(nodes[1].LocalID(),sstr);
						glRasterPos3v(&nodes[1].Coords()[0]);
						if( text == 1 || text == 5 ) if(text) printtext(sstr.str().c_str());
					}

				}

				if( matfilter == 0 || (mat.isValid() &&it->Integer(mat) == matfilter-1) )
				{
					glColor3f(0,1,0);
					if( mats_tag.isValid() ) fill_mat_str(it->LocalID(),it->IntegerArray(mats_tag),sstr);
					else fill_str(it->LocalID(),sstr);
					glRasterPos3v(cntf);
					if( text == 3 || text == 5 ) if(text) printtext(sstr.str().c_str());
				}
#if 0
				ElementArray<Cell> cells = it->getCells();

				for(ElementArray<Cell>::iterator jt = cells.begin(); jt != cells.end(); ++jt)
				{
					glColor3f(0,0,0);
					jt->Centroid(cnt);
					fill_str(jt->LocalID(),sstr.str().c_str());
					glRasterPos3v(cnt);
					if( text == 4 || text == 5 ) if(text) printtext(sstr.str().c_str());
				}

	#if defined(DISCR_DEBUG)
				Storage::real cnt0[3], cnt2[3], nrm[3], f1[3],f2[3],l;

				it->Centroid(cnt0);
				/*
				glColor3f(0,0,1);
				glBegin(GL_LINES);
				it->getAsFace()->UnitNormal(nrm);
				tensor_prod3(&it->getAsFace()->BackCell()->RealArrayDF(tensor)[0],
							 nrm,f1);
				l = dot_prod(f1,f1);
				glVertex3dv(cnt0);
				glVertex3d(cnt[0]-f1[0],cnt[1]-f1[1],cnt[2]-f1[2]);

				if(it->getAsFace()->FrontCell() != NULL )
				{
					tensor_prod3(&it->getAsFace()->BackCell()->RealArrayDF(tensor)[0],
								 nrm,f2);
					glVertex3dv(cnt0);
					glVertex3d(cnt[0]+f2[0],cnt[1]+f2[1],cnt[2]+f2[2]);
				}
				glEnd();
				*/
				for(ElementArray<Cell>::iterator jt = cells.begin(); jt != cells.end(); ++jt)
				{
					jt->Centroid(cnt);

					glColor3f(0,0,0);
					glBegin(GL_LINES);
					glVertex3v(cnt0);
					glVertex3v(cnt);
					glEnd();

					ElementArray<Face> faces = jt->getFaces();

					for(ElementArray<Face>::iterator qt = faces.begin(); qt != faces.end(); ++qt)
					{
						qt->Centroid(cnt);
						glColor3f(0,0,1);
						glBegin(GL_LINES);
						glVertex3v(cnt0);
						glVertex3v(cnt);
						glEnd();
						fill_str(qt->LocalID(),sstr);
						glRasterPos3v(cnt);
						printtext(sstr.str().c_str());

						if( !qt->Boundary() )
						{
							Storage::real stub;
							Cell * xL = qt->FrontCell();
							Cell * xK = qt->BackCell();
							Cell * xQ = xK == &*jt ? xL : xK;
							xQ->Centroid(cnt2);
							FindHarmonicPoint(&*qt,&*jt,xQ,tensor,cnt0,cnt2,cnt,stub);
							glColor3f(1,0,0);
							glBegin(GL_LINES);
							glVertex3v(cnt0);
							glVertex3v(cnt);
							glVertex3v(cnt);
							glVertex3v(cnt2);
							glEnd();

							fill_str(xQ->LocalID());
							glRasterPos3v(cnt2);
							printtext(sstr.str().c_str());
						}
						else
						{
							Cell * xK = qt->BackCell();
							xK->Centroid(cnt2);
							FindBoundaryPoint(&*qt,xK,tensor,cnt2,cnt);

							glColor3f(0,1,0);
							glBegin(GL_LINES);
							glVertex3v(cnt0);
							glVertex3v(cnt);
							glEnd();

							fill_str(qt->LocalID(),sstr);
							glRasterPos3v(cnt);
							printtext(sstr.str().c_str());
						}
					}
				}

				glEnd();
	#endif
#endif
				//mats = it->IntegerArray(mats_tag);
				//if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
				//	polygons.push_back(DrawFace(&*it,0,campos));
			}
		}
	}
	std::sort(inside_face.rbegin(),inside_face.rend());

	glEnable(GL_BLEND);

	for (size_t q = 0; q < inside_face.size(); q++) inside_face[q].draw();

	glDisable(GL_BLEND);

	//for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); it++)
#if defined(OCTREECUTCELL_DEBUG)
	if( show_cell == -1 && show_face == -1 )
	{
		int end = mesh->MaxLocalID(CELL);
		//for(int i = 0; i < end; ++i)
		if( problem_set != NULL ) for(ElementSet::iterator iit = problem_set->begin(); iit != problem_set->end(); ++iit)
		{
			Element * it = &*iit;//mesh->ElementByLocalID(CELL,i);
			//if( it == NULL ) continue;
			//if(it->HaveData(problem_tag) && it->Integer(problem_tag) != 0 )
			{
				if( parentfilter != -1 && parent_tag.isValid() && it->Integer(parent_tag) != parentfilter ) continue;
				Storage::integer_array mats;
				ElementArray<Edge> edges = it->getEdges();
				Storage::real cnt[3];
				for(unsigned k = 0; k < edges.size(); k++)
				{
					edges[k].Centroid(cnt);
					mats = edges[k].IntegerArray(mats_tag);

					ElementArray<Node> nodes = edges[k].getNodes();

					if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
					{
						fill_mat_str(edges[k].LocalID(), mats,sstr);
						glColor3f(0,0,1);
						glRasterPos3v(cnt);
						if( text == 2 || text == 5 ) if(text) printtext(sstr.str().c_str());
						glColor3f(0,0,0);
						glBegin(GL_LINES);
						glVertex3v(&nodes[0].Coords()[0]);
						glVertex3v(&nodes[1].Coords()[0]);
						glEnd();
					}

					glColor3f(1,0,0);
					mats = nodes[0].IntegerArray(mats_tag);

					if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
					{

						fill_mat_str(nodes[0].LocalID(), mats,sstr);
						glRasterPos3v(&nodes[0].Coords()[0]);
						if( text == 1 || text == 5 ) if(text) printtext(sstr.str().c_str());
					}

					mats = nodes[1].IntegerArray(mats_tag);

					if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
					{
						fill_mat_str(nodes[1].LocalID(),mats,sstr);
						glRasterPos3v(&nodes[1].Coords()[0]);
						if( text == 1 || text == 5 ) if(text) printtext(sstr.str().c_str());
					}
				}

				ElementArray<Face> faces = it->getFaces();

				for(unsigned k = 0; k < faces.size(); k++)
				{
					mats = faces[k].IntegerArray(mats_tag);
					if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
					{
						if( faces[k].Boundary() )
							glColor3f(0.15,0.45,0);
						else glColor3f(0,1,0);
						faces[k].Centroid(cnt);
						fill_mat_str(faces[k].LocalID(),mats,sstr);
						glRasterPos3v(cnt);
						if( text == 3 || text == 5 ) if(text) printtext(sstr.str().c_str());


						//~ glColor3f(0,0,0);
						//~ Storage::real nrm[3];
						//~ faces[k].OrientedNormal(it->getAsCell(),nrm);
						//~ glLineWidth(3.0);
						//~ glBegin(GL_LINES);
						//~ glVertex3v(cnt);
						//~ glVertex3d(cnt[0]+nrm[0]*20,cnt[1]+nrm[1]*20,cnt[2]+nrm[2]*20);
						//~ glEnd();
						//~ glLineWidth(1.0);
					}
					//polygons.push_back(DrawFace(&faces[k],0,campos));
				}


				if( matfilter == 0 || it->Integer(mat) == matfilter-1 )
				{
					glColor3f(0,0,0);
					it->Centroid(cnt);
					//sprintf(str,"%d, %d, %d",it->LocalID(),it->Integer(parent_tag),it->Integer(cell_material_tag)+1);
					fill_str(it->LocalID(),sstr);
					sstr << ", " << it->Integer(parent_tag), << ", " << it->Integer(cell_material_tag)+1;
					if( problem_tag.isValid() ) //sprintf(str,"%s, %d", str, it->Integer(problem_tag));
						sstr << ", " << it->Integer(problem_tag);
					glRasterPos3v(cnt);
					if( text == 4 || text == 5 ) if(text) printtext(sstr.str().c_str());

					glPointSize(3);
					glColor3f(1,0,0);
					glBegin(GL_POINTS);
					glVertex3v(cnt);
					glEnd();
					glPointSize(1);
				}






			}
		}
		if( display_orphans )
		{
			//end = mesh->MaxLocalID(EDGE);
			//for(Mesh::iteratorEdge it = mesh->BeginEdge(); it != mesh->EndEdge(); it++)
			//for(int i = 0; i < end; i++)
			if( orphan_edges != NULL ) for(ElementSet::iterator iit = orphan_edges.Begin(); iit != orphan_edges.End(); ++iit)
			{
				Element * it = &*iit;//mesh->ElementByLocalID(EDGE,i);
				//if( it == NULL ) continue;
				//if( it->nbAdjElements(CELL) == 0 )
				{
					Storage::integer_array mats = it->IntegerArray(mats_tag);
					Storage::real cnt[3];

					ElementArray<Node> nodes = it->getNodes();

					if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
					{
						it->Centroid(cnt);
						fill_mat_str(it->LocalID(), mats,sstr);
						glColor3f(1,0,1);
						glRasterPos3v(cnt);
						if( text == 2 || text == 5 ) if(text) printtext(sstr.str().c_str());


						glBegin(GL_LINES);
						glVertex3v(&nodes[0].Coords()[0]);
						glVertex3v(&nodes[1].Coords()[0]);
						glEnd();


						glColor3f(1,0,1);
						mats = nodes[0].IntegerArray(mats_tag);


					}


					//~ glColor3f(1,0,0);
					//~ mats = nodes[0].IntegerArray(mats_tag);
					//~ fill_mat_str(nodes[0].LocalID(), mats,str);
					//~ glRasterPos3v(&nodes[0].Coords()[0]);
					//~ printtext(str);
					//~
					//~ mats = nodes[1].IntegerArray(mats_tag);
					//~ fill_mat_str(nodes[1].LocalID(),mats,str);
					//~ glRasterPos3v(&nodes[1].Coords()[0]);
					//~ printtext(str);

				}
			}
			//end = mesh->MaxLocalID(FACE);
			//for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); it++)
			//for(int i = 0; i < end; i++)
			if( orphan_faces != NULL ) for(ElementSet::iterator iit = orphan_faces->begin(); iit != orphan_faces->end(); ++iit)
			{
				Element * it = &*iit;//mesh->ElementByLocalID(FACE,i);
				//if( it == NULL ) continue;
				//if( it->nbAdjElements(CELL) == 0 )
				{

					Storage::integer_array mats;
					Storage::real cnt[3];
					ElementArray<Edge> edges = it->getEdges();

					for(unsigned k = 0; k < edges.size(); k++)
					//if( edges[k].nbAdjElements(CELL) == 0 )
					{
						mats = edges[k].IntegerArray(mats_tag);

						ElementArray<Node> nodes = edges[k].getNodes();

						if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
						{
							edges[k].Centroid(cnt);
							fill_mat_str(edges[k].LocalID(), mats,sstr);
							glColor3f(0,0,1);
							glRasterPos3v(cnt);
							if( text == 2 || text == 5 ) if(text) printtext(sstr.str().c_str());

							glColor3f(0,0,0);
							glBegin(GL_LINES);
							glVertex3v(&nodes[0].Coords()[0]);
							glVertex3v(&nodes[1].Coords()[0]);
							glEnd();
						}

						glColor3f(1,0,0);
						mats = nodes[0].IntegerArray(mats_tag);
						if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
						{
							fill_mat_str(nodes[0].LocalID(), mats,sstr);
							glRasterPos3v(&nodes[0].Coords()[0]);
							if( text == 1 || text == 5 ) if(text) printtext(sstr.str().c_str());
						}

						if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
						{
							mats = nodes[1].IntegerArray(mats_tag);
							fill_mat_str(nodes[1].LocalID(),mats,sstr);
							glRasterPos3v(&nodes[1].Coords()[0]);
							if( text == 1 || text == 5 ) if(text) printtext(sstr.str().c_str());
						}
					}



					mats = it->IntegerArray(mats_tag);

					if( text == 3 || text == 5 )
					{
						it->Centroid(cnt);
						glColor3f(0,0.5,0.5);
						glRasterPos3v(cnt);
						fill_mat_str(it->LocalID(), mats,sstr);
						if(text) printtext(sstr.str().c_str());
					}

			//~ EDGE 647252
					if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
						polygons.push_back(DrawFace(&*it,0,campos));
				}
			}
		}
	}
#endif
	if( boundary )
	{
		for(ElementType etype = FACE; etype <= CELL; etype = NextElementType(etype) )
		{
			int end = mesh->LastLocalID(etype);
			//for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); it++)
			if( colort.isValid() )
			{
				for(int i = 0; i < end; i+= interactive ? pacef : 1)
				{
					Element it = mesh->ElementByLocalID(etype,i);
					if( !it.isValid() ) continue;
					if( it->Boundary() || (etype == CELL && Element::GetGeometricDimension(it->GetGeometricType()) == 2))
					{
						int mmat = color ? (!it->getAsFace()->BackCell().isValid() ? 0 : it->getAsFace()->BackCell()->Integer(colort)) : (!it->getAsFace()->BackCell().isValid() ? 0 : it->getAsFace()->BackCell()->Integer(mat));
						if( matfilter == 0 || matfilter-1 == mmat )
							polygons.push_back(DrawFace(it->self(),mmat,campos));
					}
				}
			}
			else
			{
				for(int i = 0; i < end; i+= interactive ? pacef : 1)
				{
					Element it = mesh->ElementByLocalID(etype,i);
					if( !it.isValid() ) continue;
					if (it->Boundary() || (etype == CELL && Element::GetGeometricDimension(it->GetGeometricType()) == 2))
						polygons.push_back(DrawFace(it->self(),0,campos));
				}
			}
		}
	}

	if( display_orphans )
	{
		if( orphan_edges.isValid()) for(ElementSet::iterator iit = orphan_edges.Begin(); iit != orphan_edges.End(); ++iit)
		{
			Element it = iit->self();//mesh->ElementByLocalID(EDGE,i);
			//if( it == NULL ) continue;
			//if( it->nbAdjElements(CELL) == 0 )
			{
				//~ Storage::real cnt[3];
				
				ElementArray<Node> nodes = it->getNodes();

				{
					glBegin(GL_LINES);
					glVertex3v(&nodes[0].Coords()[0]);
					glVertex3v(&nodes[1].Coords()[0]);
					glEnd();

				}
			}
		}
		if( orphan_faces.isValid() ) for(ElementSet::iterator iit = orphan_faces.Begin(); iit != orphan_faces.End(); ++iit)
		{
			Element it = iit->self();
			{

				Storage::integer_array mats;
				//~ Storage::real cnt[3];
				ElementArray<Edge> edges = it->getEdges();

				for(unsigned k = 0; k < edges.size(); k++)
				{
					ElementArray<Node> nodes = edges[k].getNodes();
					{
						glBegin(GL_LINES);
						glVertex3v(&nodes[0].Coords()[0]);
						glVertex3v(&nodes[1].Coords()[0]);
						glEnd();
					}
				}


				if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
					polygons.push_back(DrawFace(it->self(),0,campos));
			}
		}
	}


	std::sort(polygons.begin(),polygons.end());


	//glDisable(GL_DEPTH_TEST);
	//if( boundary == 1 ) glEnable(GL_BLEND);
	//else glDisable(GL_BLEND);
	glEnable(GL_BLEND);

	for(int q = polygons.size(); q > 0 ; q--) polygons[q-1].draw();

	glDisable(GL_BLEND);
//~
	//~ for(int q = 0; q < polygons.size() ; q++)
		//~ polygons[q].draw();

	//glEnable(GL_DEPTH_TEST);
	/*
	if( drawgrid)
	{
		const double scale = 2;
		double mc[3] = {(sleft+sright)*0.5,(stop+sbottom)*0.5,(snear+sfar)*0.5};
		Storage::real cc[3];
		for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); it++)
		{

			if( interactive )
			{
				int num = mesh->NumberOfCells()/1000;
				for(int k = 0; k < num; k++)
				{
					it++;
					if( it == mesh->EndCell() ) break;
				}
				if( it == mesh->EndCell() ) break;
			}


			//it->Centroid(cc);
			it->Barycenter(cc);
			double dir[3] = {cc[0] - mc[0],cc[1] - mc[1],cc[2] - mc[2]};
			{
				if( it->GetElementDimension() == 2 )
				{
					ElementArray<Node> nodes = it->getNodes();
					if( it->GetGeometricType() == MultiPolygon )
					{
						glColor3f(1,0,0);
					}
					else if( mat.isValid() )
					{
						double r = it->Integer(mat) / (double) maxmat;
						glColor3f(1-r,r,1-r);
					}
					else glColor3f(0,1,0);
					glBegin(GL_POLYGON);
					for( ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
						glVertex3d(kt->Coords()[0]+dir[0]*scale,kt->Coords()[1]+dir[1]*scale,kt->Coords()[2]+dir[2]*scale);
					glEnd();
				}
				ElementArray<Element> faces = it->getAdjElements(FACE);
				for(ElementArray<Element>::iterator f = faces.begin(); f != faces.end(); f++)// if( Geometry::Boundary(&*f) )
				{
					ElementArray<Node> nodes = f->getNodes();
					if( it->GetGeometricType() == MultiPolygon )
					{
						glColor3f(1,0,0);
					}
					else if( mat.isValid() )
					{
						double r = it->Integer(mat) / (double) maxmat;
						glColor3f(1-r,r,1-r);
					}
					else glColor3f(0,1,0);
					if( f->GetElementDimension() == 2 )
						glBegin(GL_POLYGON);
					else
						glBegin(GL_LINES);
					for(ElementArray<Node>::iterator n = nodes.begin(); n != nodes.end(); n++)
						glVertex3d(n->Coords()[0]+dir[0]*scale,n->Coords()[1]+dir[1]*scale,n->Coords()[2]+dir[2]*scale);
					glEnd();

					glColor3f(1,0,0);
					glBegin(GL_LINES);


					ElementArray<Cell> cells = f->getCells();
					for(ElementArray<Cell>::iterator c = cells.begin(); c != cells.end(); c++)
						if( &*c != &*it )
						{
							Storage::real cc2[3];
							//c->Centroid(cc2);
							c->Barycenter(cc2);
							double dir2[3] = {cc2[0] - mc[0],cc2[1] - mc[1],cc2[2] - mc[2]};
							glVertex3d(cc[0]+dir[0]*scale,cc[1]+dir[1]*scale,cc[2]+dir[2]*scale);
							glVertex3d(cc2[0]+dir2[0]*scale,cc2[1]+dir2[1]*scale,cc2[2]+dir2[2]*scale);
							//break;
						}
						glEnd();
				}
			}
			glColor3f(0,0,1);

			if( it->GetElementDimension() == 3 )
			{
				glBegin(GL_LINES);
				ElementArray<Edge> edges = it->getEdges();
				for(ElementArray<Edge>::iterator e = edges.begin(); e != edges.end(); e++)
				{
					ElementArray<Element> nodes = e->getAdjElements(NODE);
					for(ElementArray<Element>::iterator jt = nodes.begin(); jt != nodes.end(); jt++)
					{
						Storage::real_array nc = jt->getAsNode()->Coords();
						glVertex3d(nc[0]+dir[0]*scale,nc[1]+dir[1]*scale,nc[2]+dir[2]*scale);
					}
				}
				glEnd();
			}
		}
	}
	else
	{



		if( mesh->NumberOfCells() )
		{
			if( mesh->BeginCell()->GetElementDimension() == 2 )
			{
				for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); it++)
				{
					if( interactive )
					{
						int num = mesh->NumberOfCells()/1000;
						for(int k = 0; k < num; k++)
						{
							it++;
							if( it == mesh->EndFace() ) break;
						}
						if( it == mesh->EndFace() ) break;
					}
					DrawFace(&*it);
				}
			}
			else
			{
				for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); it++)
				{
					if( interactive )
					{
						int num = mesh->NumberOfCells()/1000;
						for(int k = 0; k < num; k++)
						{
							it++;
							if( it == mesh->EndFace() ) break;
						}
						if( it == mesh->EndFace() ) break;
					}
					if( it->nbAdjElements(CELL) <= 1 )
						DrawFace(&*it);
				}
			}
		}
		glColor3f(0,0,0);
		glPointSize(1.0);
		glBegin(GL_POINTS);
		for(Mesh::iteratorElement it = mesh->BeginElement(NODE); it != mesh->EndElement(); it++)
			glVertex3dv(&(it->getAsNode()->Coords()[0]));
		glEnd();
		glPointSize(1.0);



	}
	*/




	if( CommonInput != NULL )
	{
		glDisable(GL_DEPTH_TEST);
		glLoadIdentity();
		set_matrix2d();
		CommonInput->Draw();
		if( CommonInput->Done() )
		{
			delete CommonInput;
			CommonInput = NULL;
			Element it = InvalidElement();
			Storage::real cnt[3];
			if( show_cell != -1 )
			{
				it = mesh->ElementByLocalID(CELL,show_cell);

				if( !it.isValid() )
				{
					std::cout << "cell " << show_cell << " not accepted " << std::endl;
					show_cell = -1;
				}
				else
				{
					show_cells.push_back(show_cell);
					std::cout << "added cell " << show_cell << std::endl;
					show_cell = -1;
				}
			}
			if( show_face != -1 )
			{
				it = mesh->ElementByLocalID(FACE,show_face);

				if( !it.isValid() )
				{
					std::cout << "face " << show_face << " not accepted " << std::endl;
					show_face = -1;
				}
				else
				{
					show_faces.push_back(show_face);
					std::cout << "added face " << show_face << std::endl;
					show_face = -1;
				}
			}
			if( it.isValid() )
			{
				it->Centroid(cnt);
				shift[0] = -cnt[0];
				shift[1] = -cnt[1];
				shift[2] = -cnt[2];
			}
			glutPostRedisplay();
		}
		glEnable(GL_DEPTH_TEST);
	}

	glutSwapBuffers();
}

int actionstate  = 0;
double mymx = 0;
double mymy = 0;

void myclickmotion(int nmx, int nmy) // Mouse
{
	double lmx = 2.*(nmx/(double)width - 0.5),lmy = 2.*(0.5 - nmy/(double)height), dmx = lmx-mymx, dmy = lmy - mymy;
	if( actionstate == 1 )
	{
		double shiftmod[3] = {0,0,0};
		shiftmod[0] += dmx*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
		shiftmod[1] += dmy*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
		rotatevector((double*)shiftmod);
		shift[0] += shiftmod[0];
		shift[1] += shiftmod[1];
		shift[2] += shiftmod[2];
		glutPostRedisplay();
		mymx = lmx;
		mymy = lmy;
	}
	else if( actionstate == 2 )
	{
		zoom *= expf(-dmy);
		reshape(width,height);
		double shiftmod[3] = {0,0,0};
		shiftmod[2] += dmx*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
		rotatevector((double*)shiftmod);
		shift[0] += shiftmod[0];
		shift[1] += shiftmod[1];
		shift[2] += shiftmod[2];
		glutPostRedisplay();
		mymx = lmx;
		mymy = lmy;
	}
	else if( actionstate == 3 )
		clickmotion(nmx,nmy);
}
void mymotion(int nmx, int nmy) // Mouse
{
	motion(nmx,nmy);
	mymx = 2.*(nmx/(double)width - 0.5);
	mymy = 2.*(0.5 - nmy/(double)height);
}

void myclick(int b, int s, int nmx, int nmy) // Mouse
{
	if( b == GLUT_LEFT_BUTTON )
	{
		if( s == GLUT_DOWN )
		{
			actionstate = 3;
		}
		else
		{
			actionstate = 0;
		}
		click(b,s,nmx,nmy);
	}
	else if( b == GLUT_MIDDLE_BUTTON )
	{
		if( s == GLUT_DOWN )
		{
			actionstate = 1;
			interactive = true;
		}
		else
		{
			actionstate = 0;
			interactive = false;
		}
		mymx = 2.*(nmx/(double)width - 0.5);
		mymy = 2.*(0.5 - nmy/(double)height);
	}
	else if( b == GLUT_RIGHT_BUTTON )
	{
		if( s == GLUT_DOWN )
		{
			actionstate = 2;
			interactive = true;
		}
		else
		{
			actionstate = 0;
			interactive = false;
		}
		mymx = 2.*(nmx/(double)width - 0.5);
		mymy = 2.*(0.5 - nmy/(double)height);
	}
	glutPostRedisplay();
}


//void idle()
//{
//	if( interactive && reset_timer + 0.1 < Timer() )
//	{
//		interactive = false;
//		glutPostRedisplay();
//	}
//}

int main(int argc, char ** argv)
{
	table[CENTROID] = CELL | FACE;
	table[BARYCENTER] = CELL;
	table[NORMAL] = FACE;
	table[ORIENTATION] = FACE;

	Mesh * read = new Mesh();

	//~ read->PrepareGeometricData(table);
	//read->SetCommunicator(INMOST_MPI_COMM_WORLD);
	double t = Timer();
	read->Load(argv[1]);

	std::cout << "time to read: " << Timer() -t << std::endl;

	mesh = read;
	/*
	for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
		if(it->GetGeometricType() == MultiPolygon ) delete &*it;
	for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); ++it)
	{
		if(it->nbAdjElements(CELL) == 0 ) delete &*it;
	//	else it->FixNormalOrientation();
	}
	for(Mesh::iteratorEdge it = mesh->BeginEdge(); it != mesh->EndEdge(); ++it)
		if(it->nbAdjElements(FACE) == 0 ) delete &*it;
	for(Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
		if(it->nbAdjElements(NODE) == 0 ) delete &*it;
	*/



	if( argc <= 2 )
	{
		orphan_edges = mesh->CreateSet("ORPHAN_EDGES").first;
		orphan_faces = mesh->CreateSet("ORPHAN_FACES").first;
	}

	for(ElementType mask = FACE; mask >= NODE; mask = mask >> 1 )
	{
		int norphan = 0;
		for(Mesh::iteratorElement it = mesh->BeginElement(mask); it != mesh->EndElement(); ++it)
			if( it->nbAdjElements(mask << 1) == 0 )
			{
				norphan++;
				if( argc > 2 ) it->Delete();
				else if( mask == EDGE ) orphan_edges->PutElement(it->self());
				else if( mask == FACE ) orphan_faces->PutElement(it->self());
			}
		if( norphan > 0 ) std::cout << "orphan elements " << ElementTypeName(mask) << " " << norphan << std::endl;
	}

	//mesh = new Mesh(*read);

	int badorient = 0;
	std::map<Element::GeometricType,int> neighbours;
	std::map<Element::GeometricType,int> elemtypes;
	for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); ++it)
		if( !it->CheckNormalOrientation() )
		{
			badorient++;
			if( it->BackCell().isValid() ) neighbours[it->BackCell()->GetGeometricType()]++;
			if( it->FrontCell().isValid() ) neighbours[it->FrontCell()->GetGeometricType()]++;
			elemtypes[it->GetGeometricType()]++;
		}

	if( badorient )
	{
		std::cout << "Found " << badorient << " badly oriented faces" << std::endl;
		std::cout << "types of bad faces: " << std::endl;
		for(std::map<Element::GeometricType,int>::iterator it = elemtypes.begin(); it != elemtypes.end(); ++it)
			std::cout << Element::GeometricTypeName(it->first) << ": " << it->second << std::endl;
		std::cout << "neighbours of bad faces: " << std::endl;
		for(std::map<Element::GeometricType,int>::iterator it = neighbours.begin(); it != neighbours.end(); ++it)
			std::cout << Element::GeometricTypeName(it->first) << ": " << it->second << std::endl;
	}

	if( mesh->HaveTag("MATERIAL") ) mat = mesh->GetTag("MATERIAL");
	//~
	//~
#if defined(DISCR_DEBUG)
	if( mesh->HaveTag("AVGCOORD") )
		avg_coord = mesh->GetTag("AVGCOORD");
	if( mesh->HaveTag("K") )
		tensor = mesh->GetTag("K");
#endif
#if defined(OCTREECUTCELL_DEBUG)

	if( mesh->HaveTag("PROBLEM" ) )
		problem_tag = mesh->GetTag("PROBLEM");
	else
		problem_tag = mesh->GetTag("TOPOLOGY_ERROR_TAG");
	problem_set = mesh->CreateSet();
	for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
		if( it->HaveData(problem_tag) && it->Integer(problem_tag) > 0 )
			problem_set->Insert(&*it);
	parent_tag = mesh->GetTag("PARENT");
	cell_material_tag = mesh->GetTag("MATERIAL");
#endif
	if(mesh->HaveTag("MATERIALS"))
		mats_tag = mesh->GetTag("MATERIALS");

	if( mat.isValid() )
	{
		for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
			if( it->Integer(mat) > maxmat ) maxmat = it->Integer(mat);
		maxmat++;
	}

	if( mesh->HaveTag("GMSH_TAGS") )
	{
		Tag t = mesh->GetTag("GMSH_TAGS");
		for(int k = 0; k < mesh->FaceLastLocalID(); k++) if( mesh->isValidFace(k) )
		{
			Face f = mesh->FaceByLocalID(k);
			if( f->HaveData(t) && !f->Boundary() ) show_faces.push_back(k);
		}
	}


	for(ElementType mask = NODE; mask <= CELL; mask = mask << 1)
		std::cout << ElementTypeName(mask) << " " << mesh->NumberOf(mask) << std::endl;

	//for(Mesh::iteratorFace f = mesh->BeginFace(); f != mesh->EndFace(); f++)
	//	if( Geometry::FixNormaleOrientation(&*f) ) fixed++;
	//printf("fixed: %d\n",fixed);
	std::map<Element::GeometricType,int> num;
	for(Mesh::iteratorElement e = mesh->BeginElement(CELL/*|FACE|EDGE|NODE*/); e != mesh->EndElement(); ++e)
		num[e->GetGeometricType()]++;
	for(std::map<Element::GeometricType,int>::iterator it = num.begin(); it != num.end(); ++it)
		std::cout << Element::GeometricTypeName(it->first) << ": " << it->second << std::endl;

	if( mesh->GetDimensions() == 2 )
	{
		for(Mesh::iteratorNode n = mesh->BeginNode(); n != mesh->EndNode(); n++)
		{
			Storage::real_array c = n->Coords();
			if( c[0] > sright ) sright = c[0];
			if( c[0] < sleft ) sleft = c[0];
			if( c[1] > stop ) stop = c[1];
			if( c[1] < sbottom ) sbottom = c[1];
		}
		sfar = -1, snear = 1;
	}
	else
	{
		for(Mesh::iteratorNode n = mesh->BeginNode(); n != mesh->EndNode(); n++)
		{
			Storage::real_array c = n->Coords();
			if( c[0] > sright ) sright = c[0];
			if( c[0] < sleft ) sleft = c[0];
			if( c[1] > stop ) stop = c[1];
			if( c[1] < sbottom ) sbottom = c[1];
			if( c[2] > sfar ) sfar = c[2];
			if( c[2] < snear ) snear = c[2];
		}
	}
	printf("%g:%g %g:%g %g:%g\n",sleft,sright,sbottom,stop,snear,sfar);

	shift[0] = -(sleft+sright)*0.5;
	shift[1] = -(sbottom+stop)*0.5;
	shift[2] =  -(sfar+snear)*0.5;

	quatinit();
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutInitWindowPosition (100, 100);
	glutCreateWindow("Graph");

	glDepthFunc(GL_LEQUAL);
	glClearDepth(1.f);
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);


	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	glClearColor (1.0f, 1.0f, 1.0f, 1.f);
	glutDisplayFunc(draw);
	glutReshapeFunc(reshape);

	glutKeyboardFunc(keyboard);
	glutKeyboardUpFunc(keyboard2);
	glutMouseFunc(myclick);
	glutMotionFunc(myclickmotion);
	glutPassiveMotionFunc(mymotion);
	//glutIdleFunc(idle);

	glutPostRedisplay();
	glutMainLoop();
}
