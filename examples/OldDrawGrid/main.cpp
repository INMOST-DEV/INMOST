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
double l = 1e20, r = -1e20, b = 1e20, t = -1e20, zfar = -1e20, znear = 1e20;
std::map<GeometricData,ElementType> table;
void reshape(int w, int h)
{
	const double sc = 2;
	double aspect = (double)w/(double)h;
	double center[3] = { (l+r)*0.5, (b+t)*0.5, (zfar+znear)*0.5};
	double side = std::max(std::max( r-l, t-b ), zfar-znear)*0.5;
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
	if( key == '=' || key == '+')
	{
		zoom /= 1.1;
		reshape(width,height);
		glutPostRedisplay();
	}
	if( key == '_' || key == '-')
	{
		zoom *= 1.1;	
		reshape(width,height);
		glutPostRedisplay();
	}
	if( key == ' ')
	{
		drawgrid = (drawgrid+1)%2;
		glutPostRedisplay();
	}
	if( key == 'r' )
	{
		clear_half();
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

	glTranslated((l+r)*0.5,(b+t)*0.5,(znear+zfar)*0.5);
	rotate();
	glTranslated(-(l+r)*0.5,-(b+t)*0.5,-(znear+zfar)*0.5);

	//glTranslated((l+r)*0.5,(b+t)*0.5,(near+far)*0.5);
	
	if( drawgrid)
	{
		const double scale = 2;
		double mc[3] = {(l+r)*0.5,(t+b)*0.5,(znear+zfar)*0.5};
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
		if( c[0] > r ) r = c[0];
		if( c[0] < l ) l = c[0];
		if( c[1] > t ) t = c[1];
		if( c[1] < b ) b = c[1];
		if( c[2] > zfar ) zfar = c[2];
		if( c[2] < znear ) znear = c[2];
	}
	printf("%g:%g %g:%g %g:%g\n",l,r,b,t,znear,zfar);
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
	glutMouseFunc(click);
	glutMotionFunc(clickmotion);
	glutPassiveMotionFunc(motion);
	//glutIdleFunc(idle);
	
	glutPostRedisplay();
	glutMainLoop();
}
