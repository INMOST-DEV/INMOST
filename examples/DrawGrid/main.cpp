//g++ main.cpp rotate.cpp -L/usr/X11R6/lib -lX11 -lXi -lXmu -lGL -lglut -lGLU ../../inmost.a -O5
// press space - explode mesh to see connection 
//#define OCTREECUTCELL_DEBUG
#include "../../inmost.h"
#include "my_glut.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include "rotate.h"
#include <stdarg.h>

using namespace INMOST;
Mesh * mesh;
int drawgrid = 0;
int interactive = 0;
double zoom = 1;
int width = 800, height = 800;
double sleft = 1e20, sright = -1e20, sbottom = 1e20, stop = -1e20, sfar = -1e20, snear = 1e20;
Tag mat,colort; int maxmat = -1, maxcolor = 0;
std::map<GeometricData,ElementType> table;
double reset_timer = Timer();
double shift[3] = {0,0,0};
int boundary = 2;
int color = 0;
int edges = 1;
bool perspective = false;
int text = 4;
int display_orphans = 1;
ElementSet * problem_set, * orphan_edges, * orphan_faces;



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
	Input(void * link, InputType type) : input_link(link), type(type) {done = false; str = "";}
	Input(const Input & other):input_link(other.input_link), str(other.str), type(other.type), done(other.done) {}
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
	if( key == '=' || key == '+' ||  key == '_' || key == '-' || key == 'w' || key == 's' || key == 'a' || key == 'd' || key == 'r' || key == 'p' || key == 'z')
	{
		interactive = false;
		glutPostRedisplay();
	}
	
}

int matfilter = 0;
int show_cell = -1, show_face = -1;
int parentfilter = -1;

void keyboard(unsigned char key, int x, int y)
{
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
		
		drawgrid = (drawgrid+1)%2;
		glutPostRedisplay();
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
		if( show_cell == -1 )
		{
			if( CommonInput == NULL ) CommonInput = new Input(&show_cell);
		}
		else show_cell = -1;
		glutPostRedisplay();
	}
	else if( key == 'f' )
	{
		if( show_face == -1 )
		{
			if( CommonInput == NULL ) CommonInput = new Input(&show_face);
		}
		else show_face = -1;
		glutPostRedisplay();
	}
	else if( key == 'z' )
	{
		if( !colort.isValid() )
		{
			unsigned visited = 0;
			colort = mesh->CreateTag("COLOR",DATA_INTEGER,CELL,NONE,1);
			std::vector<Cell *> queue;
			adjacent<Element> around;
			queue.reserve(256);
			while(visited != mesh->NumberOfCells() )
			{
				maxcolor++;
				std::cout << "new part " << maxcolor << std::endl;
				for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
				{
					if( it->Integer(colort) == 0 )
					{
						queue.push_back(&*it);
						visited++;
						queue.back()->Integer(colort) = maxcolor;
						break;
					}
				}
				while(!queue.empty())
				{
					around = queue.back()->BridgeAdjacencies(FACE,CELL);
					
					queue.pop_back();
					for(adjacent<Element>::iterator it = around.begin(); it != around.end(); ++it)
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
	void set_center(double cnt[3], double cam[3])
	{
		dist = sqrt((cnt[0]-cam[0])*(cnt[0]-cam[0])+(cnt[1]-cam[1])*(cnt[1]-cam[1])+(cnt[2]-cam[2])*(cnt[2]-cam[2]));
	}
};



face2gl DrawFace(Element * f, int mmat, double campos[3])
{
	double cnt[3];
	face2gl ret;
	f->Centroid(cnt);
	ret.set_center(cnt,campos);
	adjacent<Node> nodes = f->getNodes();

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
	
	for( adjacent<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
		ret.add_vert(&(kt->Coords()[0]));
	

	if( edges )
	{
		glColor3f(0,0,1);
		glBegin(GL_LINE_LOOP);
		for( adjacent<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
			glVertex3dv(&(kt->Coords()[0]));
		glEnd();
	}
	return ret;
}

Tag problem_tag, mats_tag,parent_tag,cell_material_tag;

void fill_mat_str(Storage::integer id, Storage::integer_array mats, char str[2048])
{
	str[0] = '\0';
	sprintf(str,"%d",id);
	for(unsigned q = 0; q < mats.size(); q++)
		sprintf(str,"%s,%d",str,mats[q]+1);
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
	
	std::vector<face2gl> polygons;
	
	
	char str[2048];
	
	
	//glColor3f(0,0,0);
	
	if( show_cell != -1 )
	{
		Element * it = mesh->ElementByLocalID(CELL,show_cell);
		if( it == NULL )
			show_cell = -1;
		else
		{
			Storage::integer_array mats;
			adjacent<Edge> edges = it->getEdges();
			Storage::real cnt[3];
			for(unsigned k = 0; k < edges.size(); k++)
			{
				edges[k].Centroid(cnt);
				mats = edges[k].IntegerArray(mats_tag);
				
				adjacent<Node> nodes = edges[k].getNodes();
				
				
				if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
				{
					fill_mat_str(edges[k].LocalID(), mats,str);
					glColor3f(0,0,1);
					glRasterPos3dv(cnt);
					if( text == 2 || text == 5 ) if(text) printtext(str);				
					glColor3f(0,0,0);
					glBegin(GL_LINES);
					glVertex3dv(&nodes[0].Coords()[0]);
					glVertex3dv(&nodes[1].Coords()[0]);
					glEnd();
				}
				
				glColor3f(1,0,0);
				mats = nodes[0].IntegerArray(mats_tag);
				
				if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
				{
					
					fill_mat_str(nodes[0].LocalID(), mats,str);
					glRasterPos3dv(&nodes[0].Coords()[0]);
					if( text == 1 || text == 5 ) if(text) printtext(str);
				}
				
				mats = nodes[1].IntegerArray(mats_tag);
				
				if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
				{
					fill_mat_str(nodes[1].LocalID(),mats,str);
					glRasterPos3dv(&nodes[1].Coords()[0]);
					if( text == 1 || text == 5 ) if(text) printtext(str);
				}
			}
			
			adjacent<Face> faces = it->getFaces();
			
			for(unsigned k = 0; k < faces.size(); k++)
			{
				mats = faces[k].IntegerArray(mats_tag);
				if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
				{
					if( faces[k].Boundary() )
						glColor3f(0.15,0.45,0);
					else glColor3f(0,1,0);
					faces[k].Centroid(cnt);
					fill_mat_str(faces[k].LocalID(),mats,str);
					glRasterPos3dv(cnt);
					if( text == 3 || text == 5 ) if(text) printtext(str);
					
					//~ glColor3f(0,0,0);
					//~ Storage::real nrm[3];
					//~ faces[k].OrientedNormal(it->getAsCell(),nrm);
					//~ glLineWidth(2.0);
					//~ glBegin(GL_LINES);
					//~ glVertex3dv(cnt);
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
				sprintf(str,"%d, %d, %d",it->LocalID(),it->Integer(parent_tag),it->Integer(cell_material_tag));
				if( problem_tag.isValid() ) sprintf(str,"%s, %d", str, it->Integer(problem_tag));
				glRasterPos3dv(cnt);
				if( text == 4 || text == 5 ) if(text) printtext(str);
				
				glPointSize(3);
				glColor3f(1,0,0);
				glBegin(GL_POINTS);
				glVertex3dv(cnt);
				glEnd();
				glPointSize(1);
			}
		}
	}
	
	if( show_face != -1 )
	{
		Element * it = mesh->ElementByLocalID(FACE,show_face);
		if( it == NULL )
			show_face = -1;
		else
		{
			Storage::integer_array mats;
			Storage::real cnt[3];
			adjacent<Edge> edges = it->getEdges();
			
			for(unsigned k = 0; k < edges.size(); k++)
			//if( edges[k].nbAdjElements(CELL) == 0 )
			{
				mats = edges[k].IntegerArray(mats_tag);
				
				adjacent<Node> nodes = edges[k].getNodes();
				
				if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
				{
					edges[k].Centroid(cnt);
					fill_mat_str(edges[k].LocalID(), mats,str);
					glColor3f(0,0,1);
					glRasterPos3dv(cnt);
					if( text == 2 || text == 5 ) if(text) printtext(str);
					
					glColor3f(0,0,0);
					glBegin(GL_LINES);
					glVertex3dv(&nodes[0].Coords()[0]);
					glVertex3dv(&nodes[1].Coords()[0]);
					glEnd();
				}
				
				glColor3f(1,0,0);
				mats = nodes[0].IntegerArray(mats_tag);
				if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
				{
					fill_mat_str(nodes[0].LocalID(), mats,str);
					glRasterPos3dv(&nodes[0].Coords()[0]);
					if( text == 1 || text == 5 ) if(text) printtext(str);
				}
				
				if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
				{
					mats = nodes[1].IntegerArray(mats_tag);
					fill_mat_str(nodes[1].LocalID(),mats,str);
					glRasterPos3dv(&nodes[1].Coords()[0]);
					if( text == 1 || text == 5 ) if(text) printtext(str);
				}

			}

			if( matfilter == 0 || it->Integer(mat) == matfilter-1 )
			{
				glColor3f(0,1,0);
				it->Centroid(cnt);
				fill_mat_str(it->LocalID(),it->IntegerArray(mats_tag),str);
				glRasterPos3dv(cnt);
				if( text == 3 || text == 5 ) if(text) printtext(str);
			}
			
			//mats = it->IntegerArray(mats_tag);
			//if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
			//	polygons.push_back(DrawFace(&*it,0,campos));
		}
	}
	
	
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
				adjacent<Edge> edges = it->getEdges();
				Storage::real cnt[3];
				for(unsigned k = 0; k < edges.size(); k++)
				{
					edges[k].Centroid(cnt);
					mats = edges[k].IntegerArray(mats_tag);
					
					adjacent<Node> nodes = edges[k].getNodes();
					
					
					if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
					{
						fill_mat_str(edges[k].LocalID(), mats,str);
						glColor3f(0,0,1);
						glRasterPos3dv(cnt);
						if( text == 2 || text == 5 ) if(text) printtext(str);				
						glColor3f(0,0,0);
						glBegin(GL_LINES);
						glVertex3dv(&nodes[0].Coords()[0]);
						glVertex3dv(&nodes[1].Coords()[0]);
						glEnd();
					}
					
					glColor3f(1,0,0);
					mats = nodes[0].IntegerArray(mats_tag);
					
					if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
					{
						
						fill_mat_str(nodes[0].LocalID(), mats,str);
						glRasterPos3dv(&nodes[0].Coords()[0]);
						if( text == 1 || text == 5 ) if(text) printtext(str);
					}
					
					mats = nodes[1].IntegerArray(mats_tag);
					
					if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
					{
						fill_mat_str(nodes[1].LocalID(),mats,str);
						glRasterPos3dv(&nodes[1].Coords()[0]);
						if( text == 1 || text == 5 ) if(text) printtext(str);
					}
				}
				
				adjacent<Face> faces = it->getFaces();
				
				for(unsigned k = 0; k < faces.size(); k++)
				{
					mats = faces[k].IntegerArray(mats_tag);
					if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
					{
						if( faces[k].Boundary() )
							glColor3f(0.15,0.45,0);
						else glColor3f(0,1,0);
						faces[k].Centroid(cnt);
						fill_mat_str(faces[k].LocalID(),mats,str);
						glRasterPos3dv(cnt);
						if( text == 3 || text == 5 ) if(text) printtext(str);
						
						
						//~ glColor3f(0,0,0);
						//~ Storage::real nrm[3];
						//~ faces[k].OrientedNormal(it->getAsCell(),nrm);
						//~ glLineWidth(3.0);
						//~ glBegin(GL_LINES);
						//~ glVertex3dv(cnt);
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
					sprintf(str,"%d, %d, %d",it->LocalID(),it->Integer(parent_tag),it->Integer(cell_material_tag)+1);
					if( problem_tag.isValid() ) sprintf(str,"%s, %d", str, it->Integer(problem_tag));
					glRasterPos3dv(cnt);
					if( text == 4 || text == 5 ) if(text) printtext(str);
					
					glPointSize(3);
					glColor3f(1,0,0);
					glBegin(GL_POINTS);
					glVertex3dv(cnt);
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
			if( orphan_edges != NULL ) for(ElementSet::iterator iit = orphan_edges->begin(); iit != orphan_edges->end(); ++iit)
			{
				Element * it = &*iit;//mesh->ElementByLocalID(EDGE,i);
				//if( it == NULL ) continue;
				//if( it->nbAdjElements(CELL) == 0 )
				{
					Storage::integer_array mats = it->IntegerArray(mats_tag);
					Storage::real cnt[3];
				
					adjacent<Node> nodes = it->getNodes();
				
					if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
					{
						it->Centroid(cnt);
						fill_mat_str(it->LocalID(), mats,str);
						glColor3f(1,0,1);
						glRasterPos3dv(cnt);
						if( text == 2 || text == 5 ) if(text) printtext(str);
					
					
						glBegin(GL_LINES);
						glVertex3dv(&nodes[0].Coords()[0]);
						glVertex3dv(&nodes[1].Coords()[0]);
						glEnd();
					
					
						glColor3f(1,0,1);
						mats = nodes[0].IntegerArray(mats_tag);
					
					
					}
				
				
					//~ glColor3f(1,0,0);
					//~ mats = nodes[0].IntegerArray(mats_tag);
					//~ fill_mat_str(nodes[0].LocalID(), mats,str);
					//~ glRasterPos3dv(&nodes[0].Coords()[0]);
					//~ printtext(str);
					//~ 
					//~ mats = nodes[1].IntegerArray(mats_tag);
					//~ fill_mat_str(nodes[1].LocalID(),mats,str);
					//~ glRasterPos3dv(&nodes[1].Coords()[0]);
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
					adjacent<Edge> edges = it->getEdges();
				
					for(unsigned k = 0; k < edges.size(); k++)
					//if( edges[k].nbAdjElements(CELL) == 0 )
					{
						mats = edges[k].IntegerArray(mats_tag);
					
						adjacent<Node> nodes = edges[k].getNodes();
					
						if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
						{
							edges[k].Centroid(cnt);
							fill_mat_str(edges[k].LocalID(), mats,str);
							glColor3f(0,0,1);
							glRasterPos3dv(cnt);
							if( text == 2 || text == 5 ) if(text) printtext(str);
						
							glColor3f(0,0,0);
							glBegin(GL_LINES);
							glVertex3dv(&nodes[0].Coords()[0]);
							glVertex3dv(&nodes[1].Coords()[0]);
							glEnd();
						}
					
						glColor3f(1,0,0);
						mats = nodes[0].IntegerArray(mats_tag);
						if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
						{
							fill_mat_str(nodes[0].LocalID(), mats,str);
							glRasterPos3dv(&nodes[0].Coords()[0]);
							if( text == 1 || text == 5 ) if(text) printtext(str);
						}
					
						if( matfilter == 0 || std::binary_search(mats.begin(),mats.end(),matfilter-1) )
						{
							mats = nodes[1].IntegerArray(mats_tag);
							fill_mat_str(nodes[1].LocalID(),mats,str);
							glRasterPos3dv(&nodes[1].Coords()[0]);
							if( text == 1 || text == 5 ) if(text) printtext(str);
						}
					}
				
				

					mats = it->IntegerArray(mats_tag);

					if( text == 3 || text == 5 ) 
					{
						it->Centroid(cnt);
						glColor3f(0,0.5,0.5);
						glRasterPos3dv(cnt);
						fill_mat_str(it->LocalID(), mats,str);
						if(text) printtext(str);
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
		int end = mesh->MaxLocalID(FACE);
		//for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); it++)
		if( colort.isValid() )
		{
			for(int i = 0; i < end; i+= interactive ? 10 : 1)
			{
				Element * it = mesh->ElementByLocalID(FACE,i);
				if( it == NULL ) continue;
				if( it->Boundary() )
				{
					int mmat = color ? (it->getAsFace()->BackCell() == NULL ? 0 : it->getAsFace()->BackCell()->Integer(colort)) : (it->getAsFace()->BackCell() == NULL ? 0 : it->getAsFace()->BackCell()->Integer(mat));
					if( matfilter == 0 || matfilter-1 == mmat )
						polygons.push_back(DrawFace(&*it,mmat,campos));
				}
			}
		}
		else
		{
			for(int i = 0; i < end; i+= interactive ? 10 : 1)
			{
				Element * it = mesh->ElementByLocalID(FACE,i);
				if( it == NULL ) continue;
				if( it->Boundary() )
					polygons.push_back(DrawFace(&*it,0,campos));
			}
		}
	}
	
	std::sort(polygons.begin(),polygons.end());
	
	
	//glDisable(GL_DEPTH_TEST); 
	//if( boundary == 1 ) glEnable(GL_BLEND);
	//else glDisable(GL_BLEND);
	glEnable(GL_BLEND);

	for(int q = polygons.size(); q > 0 ; q--) polygons[q-1].draw();
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
					adjacent<Node> nodes = it->getNodes();
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
					for( adjacent<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
						glVertex3d(kt->Coords()[0]+dir[0]*scale,kt->Coords()[1]+dir[1]*scale,kt->Coords()[2]+dir[2]*scale);
					glEnd();
				}
				adjacent<Element> faces = it->getAdjElements(FACE);
				for(adjacent<Element>::iterator f = faces.begin(); f != faces.end(); f++)// if( Geometry::Boundary(&*f) )
				{
					adjacent<Node> nodes = f->getNodes();
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
					for(adjacent<Node>::iterator n = nodes.begin(); n != nodes.end(); n++)
						glVertex3d(n->Coords()[0]+dir[0]*scale,n->Coords()[1]+dir[1]*scale,n->Coords()[2]+dir[2]*scale);
					glEnd();

					glColor3f(1,0,0);
					glBegin(GL_LINES);


					adjacent<Cell> cells = f->getCells();
					for(adjacent<Cell>::iterator c = cells.begin(); c != cells.end(); c++)
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
				adjacent<Edge> edges = it->getEdges();
				for(adjacent<Edge>::iterator e = edges.begin(); e != edges.end(); e++)
				{
					adjacent<Element> nodes = e->getAdjElements(NODE);
					for(adjacent<Element>::iterator jt = nodes.begin(); jt != nodes.end(); jt++)
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
			Element * it = NULL;
			Storage::real cnt[3];
			if( show_cell != -1 )
			{
				it = mesh->ElementByLocalID(CELL,show_cell);
				if( it == NULL ) show_cell = -1;
			}
			else if( show_face != -1 )
			{
				it = mesh->ElementByLocalID(FACE,show_face);
				if( it == NULL ) show_face = -1;
			}
			if( it != NULL )
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
		orphan_edges = mesh->CreateSet();
		orphan_faces = mesh->CreateSet();
	}
	
	for(ElementType mask = FACE; mask >= NODE; mask = mask >> 1 )
	{
		int norphan = 0;
		for(Mesh::iteratorElement it = mesh->BeginElement(mask); it != mesh->EndElement(); ++it)
			if( it->nbAdjElements(mask << 1) == 0 ) 
			{
				norphan++;
				if( argc > 2 ) delete &*it;
				else if( mask == EDGE ) orphan_edges->Insert(&*it);
				else if( mask == FACE ) orphan_faces->Insert(&*it);
			}
		if( norphan > 0 ) std::cout << "orphan elements " << ElementTypeName(mask) << " " << norphan << std::endl;
	}
	
	//mesh = new Mesh(*read);

	
	
	if( mesh->HaveTag("MATERIAL") ) mat = mesh->GetTag("MATERIAL");
	//~ 
	//~ 
#if defined(OCTREECUTCELL_DEBUG)
	
	if( mesh->HaveTag("PROBLEM" ) )
		problem_tag = mesh->GetTag("PROBLEM");
	else 
		problem_tag = mesh->GetTag("TOPOLOGY_ERROR_TAG");
	problem_set = mesh->CreateSet();
	for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
		if( it->HaveData(problem_tag) && it->Integer(problem_tag) > 0 )
			problem_set->Insert(&*it);
	mats_tag = mesh->GetTag("MATERIALS");
	parent_tag = mesh->GetTag("PARENT");
	cell_material_tag = mesh->GetTag("MATERIAL");
#endif
	if( mat.isValid() )
	{
		for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
			if( it->Integer(mat) > maxmat ) maxmat = it->Integer(mat);
		maxmat++;
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
