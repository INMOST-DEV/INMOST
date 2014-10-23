//g++ main.cpp rotate.cpp -L/usr/X11R6/lib -lX11 -lXi -lXmu -lGL -lglut -lGLU ../../INMOST.a -O5
// press space - explode mesh to see connection 
#include "../../inmost.h"
#include "rotate.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include "my_glut.h"


using namespace INMOST;
Mesh * mesh;
int drawgrid = 0;
int interactive = 0;
double zoom = 1;
int width = 800, height = 800;
double sleft = 1e20, sright = -1e20, sbottom = 1e20, stop = -1e20, sfar = -1e20, snear = 1e20;
double shift[3] = {0,0,0};
bool perspective = false;
std::map<GeometricData,ElementType> table;

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


void clear_half()
{
	double tt = Timer();	
	int num = 0, maxnum = mesh->NumberOf(FACE) >> 1;
	if( maxnum > 0 )
	for(Mesh::iteratorElement h = mesh->BeginElement(FACE); h != mesh->EndElement(); ++h)
		//if( mesh->isValidHandle(h) )
		{
			//std::cout << "DELETE " << ElementTypeName(GetElementType(h)) << " " << LocalID(h) << std::endl;
			//mesh->DeleteElement(h,true);
			delete &*h;
			if( ++num == maxnum ) break;
		}
	/*
	before.reserve(mesh->NumberOf(CELL)+mesh->NumberOf(FACE));
	after.reserve(mesh->NumberOf(CELL)+mesh->NumberOf(FACE));
	integer i = 0;
	Tag temp = mesh->CreateTag("temp",DATA_INTEGER,CELL|FACE,NONE,1);
	for(ElementType t = FACE; t <= CELL; t = t << 1)
	for(Handle h = mesh->GetFirstElement(t); h != mesh->GetLastElement(t); ++h)
		if( mesh->isValidHandle(h) )
		{
			mesh->Integer(h,temp) = i;
			before.push_back(i++);
		}
		
	std::cout << "before CELLS: " << mesh->NumberOf(CELL) << std::endl;
	std::cout << "before FACES: " << mesh->NumberOf(FACE) << std::endl;
	*/
	double t = Timer();
	mesh->ReorderEmpty(CELL | FACE);
	std::cout << Timer()-t << std::endl;
	/*
	for(ElementType t = FACE; t <= CELL; t = t << 1)
	for(Handle h = mesh->GetFirstElement(t); h != mesh->GetLastElement(t); ++h)
	{
		assert(mesh->isValidHandle(h));
		after.push_back(mesh->Integer(h,temp));
	}
		
	assert(before.size() == after.size());
	std::sort(before.begin(),before.end());
	std::sort(after.begin(),after.end());
	for(unsigned k = 0; k < before.size(); k++)	
	{
		assert(before[k] == after[k]);
	}
	
	std::cout << "after CELLS: " << mesh->NumberOf(CELL) << std::endl;
	std::cout << "after FACES: " << mesh->NumberOf(FACE) << std::endl;
	*/
	std::cout << Timer()-tt << std::endl;
}


void keyboard(unsigned char key, int x, int y)
{
	if( key == 27 )
	{
		delete mesh;
		exit(-1);
	}
	else if( key == 'p' )
	{
		perspective = !perspective;
		reshape(width,height);
		glutPostRedisplay();
		interactive = true;
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
	if( key == ' ')
	{
		drawgrid = (drawgrid+1)%2;
		glutPostRedisplay();
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
	else if( key == 'q' )
	{
		mesh->PrepareGeometricData(table);
	}
	else if( key == 'w' )
	{
		mesh->RemoveGeometricData(table);
	}
}

void keyboard2(unsigned char key, int x, int y)
{
	if( key == '=' || key == '+' ||  key == '_' || key == '-' || key == 'w' || key == 's' || key == 'a' || key == 'd' || key == 'r' || key == 'p' || key == 'z')
	{
		interactive = false;
		glutPostRedisplay();
	}
	
}

void DrawFace(Element * f)
{
	adjacent<Node> nodes = f->getNodes();
	
	if( f->nbAdjElements(CELL) == 0 )
		glColor3f(1,0,0);
	else
		glColor3f(0,1,0);
	glBegin(GL_POLYGON);
	for( adjacent<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
		glVertex3dv(&(kt->Coords()[0]));
	glEnd();
	
	
	glColor3f(0,0,1);
	glBegin(GL_LINE_LOOP);
	for( adjacent<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
		glVertex3dv(&(kt->Coords()[0]));
	glEnd();
}

void draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	set_matrix3d();


	if( perspective )
		glTranslated(0,0,-zoom*2*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5);
	rotate();
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

	//glTranslated((l+r)*0.5,(b+t)*0.5,(near+far)*0.5);
	
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

		
			it->Centroid(cc);
			double dir[3] = {cc[0] - mc[0],cc[1] - mc[1],cc[2] - mc[2]};
			{
				if( it->GetElementDimension() == 2 )
				{
					adjacent<Node> nodes = it->getNodes();
					glColor3f(0,1,0);
					glBegin(GL_POLYGON);
					for( adjacent<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
						glVertex3d(kt->Coords()[0]+dir[0]*scale,kt->Coords()[1]+dir[1]*scale,kt->Coords()[2]+dir[2]*scale);
					glEnd();
				}
				adjacent<Element> faces = it->getAdjElements(FACE);
				for(adjacent<Element>::iterator f = faces.begin(); f != faces.end(); f++)// if( Geometry::Boundary(&*f) )
				{
					adjacent<Node> nodes = f->getNodes();
					glColor3f(0,1,0);
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
						c->Centroid(cc2);
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
	
	
	glutSwapBuffers();
}



int main(int argc, char ** argv)
{
	table[CENTROID] = CELL;
	mesh = new Mesh();
	mesh->Load(argv[1]);
	int fixed = 0;

	for(ElementType mask = NODE; mask <= CELL; mask = mask << 1)	
		std::cout << ElementTypeName(mask) << " " << mesh->NumberOf(mask) << std::endl;
	
	//for(Mesh::iteratorFace f = mesh->BeginFace(); f != mesh->EndFace(); f++)
	//	if( Geometry::FixNormaleOrientation(&*f) ) fixed++;
	printf("fixed: %d\n",fixed);
	
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
	
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	
	glClearColor (1.0f, 1.0f, 1.0f, 1.f);
	glutDisplayFunc(draw);
	glutReshapeFunc(reshape);
	
	glutKeyboardFunc(keyboard);
	glutMouseFunc(myclick);
	glutMotionFunc(myclickmotion);
	glutPassiveMotionFunc(mymotion);
	//glutIdleFunc(idle);
	
	glutPostRedisplay();
	glutMainLoop();
}
