#include "inmost.h"
#include "octgrid.h"
#include "my_glut.h"
#include "rotate.h"
#include <math.h>


#include <iomanip>
#include <iostream>

#define LOG(level,msg)  { if (log_level >= level) cout << msg << endl; }
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);

using namespace std;

struct grid thegrid;

// Variables for drawing
int draw_edges = 1;
int draw_faces = 1;
int draw_ghost = 1;
int draw_in_motion = 0; // redraw mode. true - always, false - by space press.
int draw_sem = 0;
int all_cells_count;

int cur_cell = 0; bool draw_all_cells = true;
int cur_face = 0; bool draw_all_faces = true;
int cur_edge = 0; bool draw_all_edges = true;

int from_file = 1;
int redist_after_amr = 0;

// Variables for refine and coarse
int refine_depth = 2;
double  mx = 0,  my = 0;
double rmx = 0, rmy = 0;
double base_radius = 0.01;
int action = 0;
int log_level = 0;

// Variables for MPI
int rank;
int size;
int counter = 0; // Count number of steps

// Variables for MPI and drawing
double recv_buff[2];
int wait_message = 0;
MPI_Request recv_req;
const int MAX_NODES_IN_FACE = 10;
const int MAX_PROCESSORS_COUNT = 128;
int current_proc_draw = -1;
void refresh_slaves_grid();
void send_coordinates_to_slaves(int action);

void redistribute_command();

struct drawing_face {
    int nodes_count;
    char stat;
    double** nodes;

    drawing_face()
    {
        stat = Element::Owned;
        nodes = new double*[MAX_NODES_IN_FACE];
        for (int i = 0; i < MAX_NODES_IN_FACE; i++)
            nodes[i] = new double[3];
    }
};

const int MAX_NODES_COUNT = 500000;
drawing_face** drawing_faces;
int*           drawing_faces_n;
double** sub_edges_nodes;
int*     sub_edges_nodes_n; 
char**    sub_edges_ghost;



void print_all_cells(grid* g)
{
    for(Mesh::iteratorCell it = thegrid.mesh->BeginCell(); it != thegrid.mesh->EndCell(); it++)
    {
        cout << ::rank << "-" << it->Integer(g->c_tags.base_id) << endl;
    }
}

/// Function provided to octgrid algorithm. Defines transformation from grid to grid for draw.
void transformation(double xyz[3]) 
{
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
	//xyz[2] = 4010.0 + 10.0 * (tmp[2]-0.5);
	xyz[2] = 4010.0 + 10.0 * (tmp[2]-0.5);
}

/// Function provided to octgrid algorithm. Defines transformation from grid for draw to grid.
void rev_transformation(double xyz[3]) 
{
    double tmp[3];
	tmp[0] = xyz[0];
	tmp[1] = xyz[1];
	tmp[2] = xyz[2];
	tmp[0] /= 100.0;
	tmp[0] += 0.5;
	tmp[1] /= 100.0;
	tmp[1] += 0.5;
	xyz[0] =  tmp[0];
	xyz[1] =  tmp[1];
	xyz[2] = 4010.0 + 10.0 * (tmp[2]-0.5);
	xyz[2] = (tmp[2] - 4010.0) / 10.0 + 0.5;
}


/// Function provided to octgrid algorithm. 
/// Defines that cell should be split. Returns 1 for split else returns 0.
int cell_should_split(struct grid * g, Cell cell)
{
	double r = base_radius;

    double x = cell.RealArrayDF(g->c_tags.center)[0]; 
    double y = cell.RealArrayDF(g->c_tags.center)[1];
    int c_level = cell.Integer(g->c_tags.level);

	if ((x-mx)*(x-mx)+(y-my)*(y-my) < r)
    {
        if (c_level < refine_depth) return 1;
    }

    for (int level = 2; level <= refine_depth; level++)
    {
	    if ((x-mx)*(x-mx)+(y-my)*(y-my) < r*5*(level-1))
            if (c_level < refine_depth - level + 1) 
                {
                    return 1;
                }
    }
    return 0;
}


/// Function provided to octgrid algorithm. 
/// Defines that cell should be unite. Returns 1 for unite else returns 0.
int cell_should_unite(struct grid * g, Cell cell)
{
//    return !cell_should_split(g,cell);
	const double r = base_radius;

    double x = cell.RealArrayDF(g->c_tags.center)[0];
    double y = cell.RealArrayDF(g->c_tags.center)[1];
    int c_level = cell.Integer(g->c_tags.level);

    if (c_level == refine_depth)
    {
        if ((x-mx)*(x-mx)+(y-my)*(y-my) > r) return 1;
    }
    else
    {
        double R = (refine_depth - c_level)*5*r;
        if ((x-mx)*(x-mx)+(y-my)*(y-my) > R) return 1;
    }
    return 0;
}

#if defined( __GRAPHICS__)
int show_region = 0;
int width = 800, height = 600;

/// Function calls while mouse is moving
void motion(int nmx, int nmy) // Mouse
{
    rmx = nmx; rmy = nmy;
	mx = ((nmx/(double)(width))-0.5)*((double)width/(double)height)+0.5;
	my = (1.0 - nmy/(double)height);	
	
	if(draw_sem == 0 && draw_in_motion)
	{
        if (size > 0) send_coordinates_to_slaves(0);
        gridAMR(&thegrid, 0); 
        if (redist_after_amr) redistribute_command();
        if (size > 0) refresh_slaves_grid();   
        draw_sem = 1;
	}

	glutPostRedisplay();
}

void set_color(int i)
{
    switch (i)
    {
        case 0: glColor3f(0.0, 0.0, 1.0); break;
        case 1: glColor3f(0.0, 1.0, 0.0); break;
        case 2: glColor3f(1.0, 0.0, 0.0); break;
        case 3: glColor3f(0.0, 1.0, 1.0); break;
        case 4: glColor3f(1.0, 0.0, 1.0); break;
        case 5: glColor3f(0.0, 0.5, 0.5); break;
        case 6: glColor3f(1.0, 1.0, 0.0); break;
        case 7: glColor3f(0.3, 0.8, 0.5); break;
        case 8: glColor3f(0.5, 0.0, 0.5); break;
        case 9: glColor3f(0.5, 0.5, 0.0); break;
        case 10: glColor3f(1.0, 1.0, 1.0); break;
    }
}

/// Main drawing function. MPI version draws cells of each processor with unique color.
void mpi_draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

    glTranslated(0,0,4010);
    rotate();
    glTranslated(0,0,-4010);

    if (draw_faces)
    {
        // Draw our (main proc) faces
        if (current_proc_draw == -1 || current_proc_draw == 0)
        {
            glColor3f(0.0,   0.0   ,255.0);
            int c = 0;
            int cf = 0;
            for(Mesh::iteratorCell it = thegrid.mesh->BeginCell(); it != thegrid.mesh->EndCell(); it++)
            {
                if (it->GetStatus() == Element::Ghost)
                {
                    if (draw_ghost == 0) continue;
                    glColor3f(1.0,   0.0   ,0.0);
                }
                else if (it->GetStatus() == Element::Shared) glColor3f(1.0,   1.0   ,0.0);
                else                                         glColor3f(0.0,   0.0   ,1.0);

                

                if (!draw_all_cells && cur_cell % thegrid.mesh->NumberOfCells() != c++) continue;
                
                // Draw faces
                ElementArray<Face> faces = it->getFaces();
                for (ElementArray<Face>::iterator f = faces.begin(); f != faces.end(); f++) 
                {

                    if (f->GetStatus() == Element::Ghost)        glColor3f(0.4,   0.4   ,0.0);
                    else if (f->GetStatus() == Element::Shared)  glColor3f(1.0,   1.0   ,1.0);
                    else                                         glColor3f(0.0,   0.0   ,1.0);

                    if (!draw_all_faces && cur_face % faces.size() != cf++) continue;
                    ElementArray<Node> nodes = f->getNodes();
                    glBegin(GL_POLYGON);
                    for (ElementArray<Node>::iterator n = nodes.begin(); n != nodes.end(); n++)
                    {
                        glVertex3dv(&n->RealArray(thegrid.mesh->CoordsTag())[0]);
                    }
                    glEnd();
                }
            }
        }

        // Draw other (slaves procs) faces
        for (int i = 1; i < size; i++)
        {
            if (current_proc_draw != -1 && current_proc_draw != i) continue;
            set_color(i);

            if (i > 10 && i < 13)
            {
                double colo = 40.0 + 215.0*(double(i - 11)/double(size - 11));
                colo /= 255;
                glColor3f(colo, 0.3,0.3 );
            }
            if (i >= 13)
            {
                double colo = 40.0 + 215.0*(double(i - 13)/double(size - 13));
                colo /= 255;
                glColor3f(0.3,0.3,colo);
            }

            int cf = 0;
            for (int j = 1; j < drawing_faces_n[i]; j++)
            {
                if (!draw_all_faces && cur_face % drawing_faces_n[i] != cf++) continue;
                glBegin(GL_POLYGON);
                for (int k = 0; k < drawing_faces[i][j].nodes_count; k++)
                {
                    if (drawing_faces[i][j].stat == Element::Ghost) 
                    {
                        if (draw_ghost == 0) continue;
                        glColor3f(0.4,0.4,0);
                    }
                    else
                        set_color(i);

                    glVertex3dv((drawing_faces[i][j].nodes[k]));
                }
                glEnd();
                
            }
        }       
    }

	if (draw_edges)
	{
        // Draw our edges
        glColor3f(0.0f,0.0f,0.0f);
        
        if (current_proc_draw == -1 || current_proc_draw == 0)
        {
            int c = 0;
            for (Mesh::iteratorEdge f = thegrid.mesh->BeginEdge(); f != thegrid.mesh->EndEdge(); f++)
            {
                if (!draw_all_edges && cur_edge % thegrid.mesh->NumberOfEdges() != c++)  {  glEnable(GL_LINE_STIPPLE); glLineStipple(1, 0x0FF0); }
                else                                                                     {  glDisable(GL_LINE_STIPPLE);                          }

                glBegin(GL_LINES);
                if (f->GetStatus() == Element::Ghost)        glColor3f(1.0,   0.0   ,0.0);
                else if (f->GetStatus() == Element::Shared)  glColor3f(0.8,   0.8   ,0.8);
                else                                         glColor3f(0.0,   0.0   ,0.0);
                ElementArray<Node> nodes = f->getNodes();
                glVertex3dv(&nodes[0].RealArray(thegrid.mesh->CoordsTag())[0]);
                glVertex3dv(&nodes[1].RealArray(thegrid.mesh->CoordsTag())[0]);
                glEnd();
            }
        }
        // Draw other edges
        glLineWidth(2);
        glColor3f(0.0,   0.0   ,0.0);
        for (int i = 1; i < size; i++)
        {
            if (current_proc_draw != -1 && current_proc_draw != i) continue;
            int gc = 0;
            int c = 0;
            for (int j = 0; j < sub_edges_nodes_n[i]; j+=6)
            {
                 if (!draw_all_edges && cur_edge != c++)  { glEnable(GL_LINE_STIPPLE); glLineStipple(1, 0x0FF0); }
                 else { glDisable(GL_LINE_STIPPLE);  }
                 
                 glBegin(GL_LINES);
                 if (sub_edges_ghost[i][gc] == Element::Shared)
                    glColor3f(0.8,   0.8   ,0.8);
                 else if (sub_edges_ghost[i][gc] == Element::Ghost)
                    glColor3f(0.1,   0.0   ,0.0);
                 else
                    glColor3f(0.0,   0.0   ,0.0);
                glVertex3dv(&(sub_edges_nodes[i][j]));
                glVertex3dv(&(sub_edges_nodes[i][j+3]));
                gc++;
                glEnd();
            }
        }       

    }	
    glEnable(GL_POINT_SMOOTH);
    glBegin(GL_POINTS);
    for(Mesh::iteratorNode it = thegrid.mesh->BeginNode(); it != thegrid.mesh->EndNode(); it++)
    {
        if (it->GetStatus() == Element::Ghost) glColor3f(1, 1, 1);
        else if (it->GetStatus() == Element::Shared) glColor3f(1, 0.0, 0.0);
        else glColor3f(0.0, 0.0, 0.0);
        glVertex3dv(&(it->getAsNode().RealArray(thegrid.mesh->CoordsTag())[0]));
    }
    glEnd();
    glDisable(GL_POINT_SMOOTH);

    glutSwapBuffers();
    draw_sem = 0;
}

/// Main drawing function. 
void draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glTranslated(0,0,4010);
	rotate();
	glTranslated(0,0,-4010);

	if (draw_faces)
	{
		int i = 0;
		for (Mesh::iteratorFace f = thegrid.mesh->BeginFace(); f != thegrid.mesh->EndFace(); f++)
		{
			double c[3];
			ElementArray<Node> nodes = f->getNodes();
			glColor3f(0.0f,0.0f,1.0f);	
			glBegin(GL_POLYGON);
			for (ElementArray<Node>::iterator n = nodes.begin(); n != nodes.end(); n++)
				glVertex3dv(&n->RealArray(thegrid.mesh->CoordsTag())[0]);
			glEnd();
		}
	}
	
	if (draw_edges)
	{
        glLineWidth(2);
		glColor3f(0.0f,0.0f,0.0f);
		{
			glBegin(GL_LINES);
			for (Mesh::iteratorEdge f = thegrid.mesh->BeginEdge(); f != thegrid.mesh->EndEdge(); f++)
			{
				ElementArray<Node> nodes = f->getNodes();
				glVertex3dv(&nodes[0].RealArray(thegrid.mesh->CoordsTag())[0]);
				glVertex3dv(&nodes[1].RealArray(thegrid.mesh->CoordsTag())[0]);
			}
			glEnd();
		}
        glLineWidth(1);
	}	
 
    glutSwapBuffers();
    draw_sem = 0;
}


void reshape(int w, int h)
{
	double aspect = (double)w/(double)h;
	width = w;
	height = h;
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho(-50*aspect,50*aspect,-50,50,-5000,5000);
	glMatrixMode (GL_MODELVIEW);
	glViewport(0, 0, w, h);
}

/// Send mpuse coordinates to slaves. Slaves receive the message 'm', execute "refine" with received coordinates.
void send_coordinates_to_slaves(int action)
{
    char buff[MAX_PROCESSORS_COUNT][10];
    MPI_Request req[MAX_PROCESSORS_COUNT];

    for (int i = 1; i < size; i++)
    {
        LOG(2,"Main process: send coordinates to " << i)
        buff[i][0] = 'm'; // Special key, means refine 
        *((double*)(buff[i] + 1)) = mx;
        *((double*)(buff[i] + 1 + sizeof(double))) = my;
        *((int*)(buff[i] + 1 + 2*sizeof(double))) = action;
        MPI_Isend(buff[i], 1 + 2 * sizeof(double) + sizeof(int), MPI_CHAR, i, 0, INMOST_MPI_COMM_WORLD, req + i);
    }
}

/// OpenGL drawing execute only on one processor. In our case it's process with rank = 0 (Master). Every time when grid should be redraw, master sends
/// special command 'm' to slaves. Slave receive the message 'm', sends his cells to master.
void refresh_slaves_grid()
{
	if (::rank == 0) { // Send command to slaves for they will sent grid to us
        char buff1[10][2];
        MPI_Request req[10];
        for (int i = 1; i < size; i++)
        {
            buff1[i][0] = 'r'; // Special key, means send mesh
            MPI_Isend(buff1[i], 1, MPI_CHAR, i, 0, INMOST_MPI_COMM_WORLD, req + i);
        }
    }

    // Waiting response from slaves
    char buff[MAX_NODES_COUNT * sizeof(double)];
    MPI_Status status;
    int count_e;
    int count_df = 0;
    char nodes_count;
    char stat;

    for (int i = 1; i < size; i++)
    {
        // First receive array of vertices of faces
        MPI_Recv(buff, MAX_NODES_COUNT * sizeof(double),MPI_CHAR,i,0,INMOST_MPI_COMM_WORLD, &status);
        count_df = *((int*)buff);
        LOG(2,"Main process: received faces from " << i << " with " << count_df << " faces")

        drawing_faces_n[i] = count_df;
        int offset = sizeof(int);
        for (int j = 0; j < count_df; j++)
        {
            nodes_count = *(buff + offset++); // Read count of vertex in face
            stat = *(buff + offset++);
            drawing_faces[i][j].nodes_count = nodes_count;
            drawing_faces[i][j].stat = stat;
            for (int k = 0; k < nodes_count; k++)
            {
                drawing_faces[i][j].nodes[k][0] = *((double*)(buff + offset)); 
                drawing_faces[i][j].nodes[k][1] = *((double*)(buff + offset + sizeof(double)  )); 
                drawing_faces[i][j].nodes[k][2] = *((double*)(buff + offset + sizeof(double)*2)); 
                offset += sizeof(double)*3;
            }
        }

        // Now receive array of vertices of edges
        MPI_Recv(buff, MAX_NODES_COUNT * sizeof(double),MPI_CHAR,i,0,INMOST_MPI_COMM_WORLD, &status);
        count_e = *((int*)buff);
        sub_edges_nodes_n[i] = (count_e/7)*6;
        LOG(2,"Main process: received edges from " << i << " with " << count_e << " doubles")

        int c = 0;
        int gc = 0;
        while (c+gc < count_e)
        {
            sub_edges_ghost[i][gc] = *((double*)(buff + sizeof(int) + sizeof(double)*(c+gc)));
            gc++;
            for (int j = 0; j < 6 ; j++)
            {
                sub_edges_nodes[i][c] = *((double*)(buff + sizeof(int) + sizeof(double)*(c+gc))); 
                c++;
            }
        }
    }
}

void prepare_to_correct_brothers()
{
	correct_brothers(&thegrid,size,::rank, 0);
    thegrid.mesh->RemoveGhost();
    thegrid.mesh->Redistribute(); 
    thegrid.mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
    thegrid.mesh->AssignGlobalID(CELL | EDGE | FACE | NODE);
}



/// Prepare to redistribute. Master sends special command to slaves which means redistribute
void redistribute_command()
{
    char buff[MAX_PROCESSORS_COUNT][10];
    MPI_Request req[MAX_PROCESSORS_COUNT];
    int type = 2;
    
    if (counter == 5)
    {
        counter = 0;
        type = 1;
    }
   
    for (int i = 1; i < size; i++)
    {
		LOG(3,"Master: send redistribute command to slave " << ::rank)
        buff[i][0] = 'x'; // Special key, means redistribute
        buff[i][1] = type + '0'; // Special key, means redistribute
        MPI_Isend(buff[i], 2, MPI_CHAR, i, 0, INMOST_MPI_COMM_WORLD, req + i);
    }       
    
    redistribute(&thegrid, type);
    counter++;

}


void keyboard(unsigned char key, int x, int y)
{
	(void) x, (void) y;
	if( key == 27 )
	{
		exit(-1);
	}
    if( key == '+')
    {
        cur_cell++;
		glutPostRedisplay();
    }
    if( key == '-')
    {
        cur_cell--;
		glutPostRedisplay();
    }
    if( key == 'q')
    {
        draw_all_cells = !draw_all_cells;
		glutPostRedisplay();
    }
    if( key == 'a')
    {
        draw_all_faces = !draw_all_faces;
		glutPostRedisplay();
    }
    if( key == 'w')
    {
        draw_all_edges = !draw_all_edges;
		glutPostRedisplay();
    }
    if( key == '6')
    {
        cur_edge++;
		glutPostRedisplay();
    }
    if( key == '4')
    {
        cur_edge--;
		glutPostRedisplay();
    }
    if( key == '>')
    {
        cur_face++;
		glutPostRedisplay();
    }
    if( key == '<')
    {
        cur_face--;
		glutPostRedisplay();
    }
    if( key == 'g' || key == 'G')
    {
        draw_ghost = !draw_ghost;
		glutPostRedisplay();
    }
    if( key == 'o' || key == 'O')
    {
        if (action == 0) action = 3;
        else action = 0;
		glutPostRedisplay();
    }
    if( key == 'p' || key == 'P')
    {
        redist_after_amr = !redist_after_amr; 
        draw_in_motion = !draw_in_motion;
		glutPostRedisplay();
    }
    if( key == 'c' || key == 'C')
    {
        char buff[10][1];
        MPI_Request req[10];
        for (int i = 1; i < size; i++)
        {
            buff[i][0] = 'c'; // Special key, means command
            MPI_Isend(buff[i], 1, MPI_CHAR, i, 0, INMOST_MPI_COMM_WORLD, req + i);
        } 
        command(&thegrid);
        
		glutPostRedisplay();
    }
    if( key == 'u' || key == 'U')
    {
        char buff[10][1];
        MPI_Request req[10];
        for (int i = 1; i < size; i++)
        {
            buff[i][0] = 'u'; // Special key, means remove_ghost
            MPI_Isend(buff[i], 1, MPI_CHAR, i, 0, INMOST_MPI_COMM_WORLD, req + i);
        } 
        thegrid.mesh->RemoveGhost();
		glutPostRedisplay();
    }
	if( key == ' ' ) 
	{
		if (::rank == 0)
        {
            send_coordinates_to_slaves(action);
		}

        gridAMR(&thegrid, action);
        //thegrid.mesh->AssignGlobalID(CELL | EDGE | FACE | NODE);
        if (redist_after_amr)
            redistribute_command();
		glutPostRedisplay();
	}
    if( key == '[' ) 
	{
		if (::rank == 0)
        {
            send_coordinates_to_slaves(action);
		}

        gridAMR(&thegrid, 1);
        //thegrid.mesh->AssignGlobalID(CELL | EDGE | FACE | NODE);
		glutPostRedisplay();
	}
    if( key == 'r' || key == 'R'|| key == ' ') 
    {
        refresh_slaves_grid();   
    }
    if( key == 'f' ) 
	{
		if (::rank == 0) { // Send dump command to other process
            dump_to_vtk(&thegrid);
            char buff[10][2];
            MPI_Request req[10];
            for (int i = 1; i < size; i++)
            {
                buff[i][0] = 'f'; // Special key, means dump file
                MPI_Isend(buff[i], 1, MPI_CHAR, i, 0, INMOST_MPI_COMM_WORLD, req + i);
            }
        }
		glutPostRedisplay();
	}
    if( key == 'x' || key == 'X' )
    {
        redistribute_command();
        refresh_slaves_grid();   
	glutPostRedisplay();
    }
    if( key == 'z' || key == 'Z' )
    {
        char buff[10][10];
        MPI_Request req[10];

        for (int i = 1; i < size; i++)
        {
            buff[i][0] = 'z'; // Special key, means correct_brothers
            MPI_Isend(buff[i], 1, MPI_CHAR, i, 0, INMOST_MPI_COMM_WORLD, req + i);
        }       
        prepare_to_correct_brothers();
    }
	if( key == '*' || key == '8' )
    {
        current_proc_draw++;
        if (current_proc_draw == size) current_proc_draw = -1;
	    glutPostRedisplay();
    }
	if( key == '/' || key == '?' )
    {
        current_proc_draw--;
        if (current_proc_draw == -2) current_proc_draw = size - 1;
	    glutPostRedisplay();
    }
	if( key == 'e' || key == 'E' )
	{
        draw_edges = !draw_edges;
	}
	if( key == 't' || key == 'T' )
	{
        draw_faces = !draw_faces;
	}
    if( key == 'i' || key == 'I' )
    {
        draw_in_motion = !draw_in_motion ;
    }
}

void idle(void)
{
}
#endif

/// Main loop for not main process
void NotMainProcess()
{
    char buff[256];
    MPI_Status status;
    int length;

    while (1) {
        MPI_Recv(buff, 256, MPI_CHAR, 0,0, INMOST_MPI_COMM_WORLD, &status);

		LOG(2, "Process " << ::rank << ": received message '" << buff[0] << "'")
        
        if (buff[0] == 'm') // Need to refine, mouse coordinates come
        {
            mx = *((double*)(buff + 1));
            my = *((double*)(buff + 1 + sizeof(double)));
            int action = *((int*)(buff + 1 + 2*sizeof(double)));

            gridAMR(&thegrid, action); 
            //thegrid.mesh->AssignGlobalID(CELL | EDGE | FACE | NODE);
        }
        if (buff[0] == 'x') // Need to redistribute 
        {
            char type_C = buff[1];
            int  type = type_C - '0';
            redistribute(&thegrid,type);
        }
        if (buff[0] == 'u') // Need remove ghosts
        {
            thegrid.mesh->RemoveGhost();
        }
        if (buff[0] == 'c') // Need command
        {
            command(&thegrid);
        }
        if (buff[0] == 'z') // Need to correct brothers 
        {
            prepare_to_correct_brothers();
        }
        if (buff[0] == 'f') // Dump mesh to file
        {
            dump_to_vtk(&thegrid);
        }
        if (buff[0] == 'r') // Master wants to redraw mesh. Now this procces need send own grid to master.
        {
            // First fill the buffer for sent
            char buff_f[MAX_NODES_COUNT * sizeof(double)];
            char buff_e[MAX_NODES_COUNT * sizeof(double)];
            int count_e = 0;
            int offset = sizeof(int);
            int count_df = 0; // Number of transmitted faces in array buff_f
           
            // Fill array buff_f. Strcuture of array:
            // [Faces count] [Vertices' count in face] [Status] [Points] [Vertices' count in face] [Status] [Points]...
            for(Mesh::iteratorCell it = thegrid.mesh->BeginCell(); it != thegrid.mesh->EndCell(); it++)
            {
                ElementArray<Face> faces = it->getFaces();
                for(ElementArray<Face>::iterator f = faces.begin(); f != faces.end(); f++) 
                {
                    count_df++;

                    ElementArray<Node> nodes = f->getNodes();
                    if (nodes.size() == 0) continue;

                    *(buff_f + offset++) = (char)nodes.size();
                    *(buff_f + offset++) = (char)it->GetStatus();

                    for(ElementArray<Node>::iterator n = nodes.begin(); n != nodes.end(); n++)
                    {
                        *((double*)(buff_f + offset                   )) = n->RealArray(thegrid.mesh->CoordsTag())[0];
                        *((double*)(buff_f + offset + sizeof(double)  )) = n->RealArray(thegrid.mesh->CoordsTag())[1];
                        *((double*)(buff_f + offset + sizeof(double)*2)) = n->RealArray(thegrid.mesh->CoordsTag())[2];
                        offset  += 3*sizeof(double);
                    }
                }
            }
            *((int*)buff_f) = count_df;

            // Fill array buff_e
            for(Mesh::iteratorEdge f = thegrid.mesh->BeginEdge(); f != thegrid.mesh->EndEdge(); f++)
			{
                ElementArray<Node> nodes = f->getNodes();
				*((double*)(buff_e + sizeof(int) + sizeof(double)*count_e++)) = (double)f->GetStatus();
				*((double*)(buff_e + sizeof(int) + sizeof(double)*count_e++)) = nodes[0].RealArray(thegrid.mesh->CoordsTag())[0];
				*((double*)(buff_e + sizeof(int) + sizeof(double)*count_e++)) = nodes[0].RealArray(thegrid.mesh->CoordsTag())[1];
				*((double*)(buff_e + sizeof(int) + sizeof(double)*count_e++)) = nodes[0].RealArray(thegrid.mesh->CoordsTag())[2];
				*((double*)(buff_e + sizeof(int) + sizeof(double)*count_e++)) = nodes[1].RealArray(thegrid.mesh->CoordsTag())[0];
				*((double*)(buff_e + sizeof(int) + sizeof(double)*count_e++)) = nodes[1].RealArray(thegrid.mesh->CoordsTag())[1];
				*((double*)(buff_e + sizeof(int) + sizeof(double)*count_e++)) = nodes[1].RealArray(thegrid.mesh->CoordsTag())[2];
			}
            *((int*)buff_e) = count_e;

            // Now we are ready to transmit the data. Transmit with 2 steps
			LOG(2, "Process " << ::rank << ": send buffer with faces = " << count_df << ". Edges " << count_e)
            MPI_Send(buff_f, offset                             ,MPI_CHAR,0,0,INMOST_MPI_COMM_WORLD);
            MPI_Send(buff_e,count_e*sizeof(double) + sizeof(int),MPI_CHAR,0,0,INMOST_MPI_COMM_WORLD);
        }
    }
}

void parse_arguments(int argc, char** argv, int* n, double* R, int* L, string& file_name)
{
  if (argc < 2) return;

  string str;
  string str1, str2;
  for (int i = 1; i < argc; i++)
  {
    str = argv[i];
    size_t pos = str.find('=');
    if (pos == string::npos)
    {
      cout << "Invalid argument: " << str << endl;
      continue;
    }

    str1 = str.substr(0, pos);
    str2 = str.substr(pos + 1);

    if (str1 == "-n")
    {
      for (int j = 0; j < 3; j++)
      {
        pos = str2.find('x');
        if (j < 2 && pos == string::npos) 
        {
          cout << "Invalid command: " << str2 << endl;
          break;
        }
        str1 = str2.substr(0, pos);
        n[j] = atoi(str1.c_str());
        str2 = str2.substr(pos + 1);
      }
    }
    else if (str1 == "-r")
    {

	    *R = atof(str2.c_str());
    }
    else if (str1 == "-l")
    {
	    *L = atoi(str2.c_str());
    }
    else if (str1 == "-f")
    {
        file_name = str2;
    }
    else
    {
      cout << "Invalid command: " << str1 << endl;
    }
  }
}

void print_help()
{
  cout << "Example of Octree refine on redistributed grid" << endl;
  cout << "Command arguments:" << endl;
  cout << "   -n=10x10x1     - grid size" << endl;
  cout << "   -r=0.01        - refine radius" << endl;
  cout << "   -l=2           - refine level"  << endl;
  cout << "   -f=grid_3.pvtk - takes grid from file" << endl;
  cout << endl;
  cout << "Hotkeys:" << endl;
  cout << "   Space - refine grid around mouse cursor" << endl;
  //cout << "       r - redraw grid" << endl;
  cout << "       f - dump grid to file (see grids folder)" << endl;
  cout << "       x - redistribute grid" << endl;
}


int main(int argc, char ** argv)
{
	int i;
	int n[3] = {10,10,1};


	thegrid.transformation = transformation;
	thegrid.rev_transformation = rev_transformation;
	thegrid.cell_should_split = cell_should_split;
	thegrid.cell_should_unite = cell_should_unite;
    Mesh::Initialize(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &::rank);

    if (::rank == 0) print_help();
    string file_name;
    parse_arguments(argc, argv, n, &base_radius,&refine_depth, file_name);
    all_cells_count = n[0]*n[1]*n[2] * 2;
    
    if (file_name.empty()) {
        gridInit(&thegrid,n);
    } else {   
        thegrid.mesh = new Mesh(); // Create an empty mesh
        thegrid.mesh->SetCommunicator(INMOST_MPI_COMM_WORLD); // Set the MPI communicator for the mesh
        file_name = "grids/" + file_name;
        thegrid.mesh->Load(file_name.c_str()); // Load input mesh
        //	thegrid.mesh->Load("dump.pvtk"); // Load input mesh
        init_mesh(&thegrid); 
    }

    Partitioner::Initialize(&argc,&argv);

	::size = thegrid.mesh->GetProcessorsNumber();
	::rank = thegrid.mesh->GetProcessorRank();

    //dump_to_vtk();

	if (::rank != 0)
    {
        NotMainProcess();
    }
    else
    {

        #if defined(USE_MPI)
        drawing_faces     = new drawing_face*[size-1];
        drawing_faces_n   = new int[size-1];
        sub_edges_nodes   = new double*[size-1];
        sub_edges_nodes_n = new int[size-1];
        sub_edges_ghost   = new char*[size-1];
        
        for (int i = 0; i < size; i++) {
            drawing_faces_n[i] = 0;
            drawing_faces[i] = new drawing_face[30000];
            sub_edges_nodes_n[i] = 0;
            sub_edges_ghost[i] = new char[25000];
            sub_edges_nodes[i] = new double[75000];
        }
        #endif

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

        glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);
        glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
        glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
        glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);

        glClearColor (1.0f, 1.0f, 1.0f, 1.f);
        if (size > 1) glutDisplayFunc(mpi_draw);
        else glutDisplayFunc(draw);
        glutReshapeFunc(reshape);

        glutKeyboardFunc(keyboard);
        glutMotionFunc(rotate_clickmotion);
        glutMouseFunc(rotate_click);
        //glutPassiveMotionFunc(rotate_motion);
        glutPassiveMotionFunc(motion);
        glutIdleFunc(idle);

        glutPostRedisplay();
        glutMainLoop();
    }
}
