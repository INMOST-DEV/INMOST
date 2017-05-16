
//make all && mpirun -np 1 ./Octree_test -log=3 -l=2 -n=50x50x2 
#include "octgrid.h"
#include <math.h>
#include "inmost.h"
#include <iomanip>
#include <iostream>
#include <sstream>

#define LOG(level,msg)  { if (log_level >= level) cout << msg << endl; }
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);

using namespace std;

const int MAX_PROCESSORS_COUNT = 10;

struct grid thegrid;
int refine_depth = 2;
double  mx = 0,  my = 0;
double rmx = 0, rmy = 0;
double base_radius = 0.01;
int action = 0;
int log_level = 0;
int all_cells_count;
int rank;
int size;
int counter = 0; // Count number of steps

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

void pre_redistribute(int type)
{
    redistribute(&thegrid, type);
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
            if (c_level < refine_depth - level + 1) return 1;
    }
    return 0;
}

/// Main loop for not main process
void NotMainProcess()
{

    char buff[256];
    MPI_Status status;
    int length;
    int iteration = 0;

    while (1) 
    {
        MPI_Recv(buff, 256, MPI_CHAR, 0,0, INMOST_MPI_COMM_WORLD, &status);

		LOG(2, "Process " << ::rank << ": received message '" << buff[0] << "'")
        
        if (buff[0] == 'q') // Need to refine, mouse coordinates come
		{
			break;

		}
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
            //redistribute(&thegrid,type);
            pre_redistribute(type);
            cout << ::rank << ": iteration " << iteration++ << " complete. Cells: " << thegrid.mesh->NumberOfCells() << endl;
        }
        if (buff[0] == 'u') // Need remove ghosts
        {
            thegrid.mesh->RemoveGhost();
        }
        if (buff[0] == 'c') // Need command
        {
            command(&thegrid);
        }
        if (buff[0] == 'f') // Dump mesh to file
        {
            dump_to_vtk(&thegrid);
        }
   }
}

void send_dump_command()
{
	char buff[10][2];
	MPI_Request req[10];
	for (int i = 1; i < size; i++)
	{
		buff[i][0] = 'f'; // Special key, means dump file
		MPI_Isend(buff[i], 1, MPI_CHAR, i, 0, INMOST_MPI_COMM_WORLD, req + i);
	}
}

void send_quit_command()
{
	char buff[10][2];
	MPI_Request req[10];
	for (int i = 1; i < size; i++)
	{
		buff[i][0] = 'q'; // Special key, means dump file
		MPI_Isend(buff[i], 1, MPI_CHAR, i, 0, INMOST_MPI_COMM_WORLD, req + i);
	}
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
    
//    redistribute(&thegrid, type);
    pre_redistribute(type);
    counter++;

}


void print_help()
{
  cout << "Example of Octree refine on redistributed grid" << endl;
  cout << "Command arguments:" << endl;
  cout << "   -n=10x10x1 - grid size" << endl;
  cout << "   -r=0.01    - refine radius" << endl;
  cout << "   -l=2       - refine level"  << endl;
  cout << "   -log=1     - log level"  << endl;
  cout << endl;
  cout << "Hotkeys:" << endl;
  cout << "   Space - refine grid around mouse cursor" << endl;
  //cout << "       r - redraw grid" << endl;
  cout << "       f - dump grid to file (see grids folder)" << endl;
  cout << "       x - redistribute grid" << endl;
}

void parse_arguments(int argc, char** argv, int* n, double* R, int* L, int* log, int* iters_count, double* h)
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
    else if (str1 == "-h")
    {
	    sscanf(str2.c_str(),"%lf",h);
    }
    else if (str1 == "-log")
    {
	    *log = atoi(str2.c_str());
    }
    else if (str1 == "-i")
    {
	    *iters_count = atoi(str2.c_str());
    }
    else
    {
      cout << "Invalid command: " << str1 << endl;
    }
  }
}


int main(int argc, char ** argv)
{
	int i;
	int n[3] = {10,10,1};
  int iters_count = 10;
  double h = 0.02;

	thegrid.transformation = transformation;
	thegrid.rev_transformation = rev_transformation;
	thegrid.cell_should_split = cell_should_split;
	thegrid.cell_should_unite = cell_should_unite;
    Mesh::Initialize(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &::rank);

    if (::rank == 0) print_help();
    parse_arguments(argc, argv, n, &base_radius,&refine_depth,&log_level,&iters_count,&h);
    all_cells_count = n[0]*n[1]*n[2] * 2;

    gridInit(&thegrid,n);
    Partitioner::Initialize(&argc,&argv);

	::size = thegrid.mesh->GetProcessorsNumber();
	::rank = thegrid.mesh->GetProcessorRank();


    //dump_to_vtk();
	if (::rank == 0)
  {
  	cout << "Test start: " << endl;
  	cout << "Size: " << n[0] << "x" << n[1] << "x" << n[2] << endl;
  	cout << "Radius: " << base_radius << endl;
  	cout << "Refine depth: " << refine_depth << endl;
  	cout << "Iters count: " << iters_count << endl;
  	cout << "Step: " << h << endl;
  	cout << "Processors: " << ::size << endl;
  } 

      {
		mx = 0.1;
		my = 0.6;
        // double h = 0.8 / iters_count;
//        double h = 0.03;
        int i = 0;
   		BARRIER
		double st = Timer();
        double ct , tt;
        double time_amr, time_red;
        double a_amr = 0;
        double a_red = 0;
        int sends = 0;
        bool forward = true;
        for (int iter = 0; iter < iters_count; iter++)
        {
    		BARRIER
            if (::rank == 0) LOG(1, "Iteration: " << i)

            ct = Timer();
            gridAMR(&thegrid,0);
    		BARRIER
            tt = Timer();
            time_amr = tt-ct;
            
            fill_proc_tag(&thegrid);
            
            ct = Timer();
            redistribute(&thegrid, 0);
    		BARRIER
            tt = Timer();
            time_red = tt-ct;
            
            int l_sends = calc_sends(&thegrid);
            sends += l_sends;

            ct = Timer();
            LOG(2, ::rank << ": iteration " << i << " complete. Cells: " << thegrid.mesh->NumberOfCells() << ". Sends: " << l_sends << ". All sends: " << sends)
            if (iter > 0) 
            {
                a_amr += time_amr;
                a_red += time_red;
            }
    		if (::rank == 0) LOG(1, "AMR time = " << time_amr);
    		if (::rank == 0) LOG(1, "Red time = " << time_red);
            stringstream suffix;
            suffix << "_" << iter;
//            dump_to_vtk(&thegrid,suffix.str().c_str());
    		BARRIER
    		if (::rank == 0) LOG(1, "===============");
			i++;
            if (forward) {
                mx += h;
                if (mx >= 0.9) forward = false;
            } else {
                mx -= h;
                if (mx <= 0.1) forward = true;
            }
        }

    	BARRIER
        LOG(2, ::rank << ": test completed. Average sends: " << double(sends)/double(iters_count))
        int* g_sends = new int[size];
        MPI_Allgather(&sends, 1, MPI_INTEGER, g_sends, 1, MPI_INTEGER, MPI_COMM_WORLD);

        int av_sends = 0;
        if (::rank == 0) 
            for (int i = 0; i < size; i++)
                av_sends += g_sends[i];
    

    	BARRIER
		tt = Timer() - st;
		if (::rank == 0) cout << "time = " << tt << endl;
		if (::rank == 0) cout << "Average AMR time = " << a_amr/(iters_count - 1) << endl;
		if (::rank == 0) cout << "Average RED time = " << a_red/(iters_count - 1) << endl;
		if (::rank == 0) cout << "Average sends = " << int(double(av_sends)/double(iters_count*size)) << endl;
		if (::rank == 0) cout << "time = " << tt << endl;

		dump_to_vtk(&thegrid);
//		send_dump_command();
//		send_quit_command();
    }

	Mesh::Finalize();
}
