#include "octgrid.h"
#include "my_glut.h"
#include "rotate.h"
#include <math.h>
#include "inmost.h"

#include <iomanip>
#include <iostream>

#define LOG(level,msg)  { if (log_level >= level) cout << msg << endl; }
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);

using namespace std;

struct grid thegrid;


// Variables for refine and coarse
int refine_depth = 1;
double  mx = 0.75,  my = 0.5;
double base_radius = 0.02;
int log_level = 3;

// Variables for MPI
int rank;
int size;
int counter = 0; // Count number of steps




/// Dump mesh to vtk file in folder "grids"
void dump_to_vtk()
{
	//thegrid.mesh->ResolveShared(); // Resolve duplicate nodes
	//thegrid.mesh->ExchangeGhost(2,NODE); // Construct Ghost cells in 2 layers connected via nodes

    std::stringstream filename;
    filename << "grids/grid_";
    filename << size;
    if( size == 1 )
        filename << ".vtk";
    else
        filename << ".pvtk";
    thegrid.mesh->Save(filename.str());
	cout << "Process " << ::rank << ": dumped mesh to file" << endl;
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
/// Defines that cell should be unite. Returns 1 for unite else returns 0.
int cell_should_unite(struct grid * g, Cell cell)
{
	const double r = base_radius;
	int test = 1;	

    double x = cell.RealArrayDF(g->c_tags.center)[0];
    double y = cell.RealArrayDF(g->c_tags.center)[1];

	test &= (x-mx)*(x-mx)+(y-my)*(y-my) > r;
	return test;
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


void prepare_to_correct_brothers()
{
	correct_brothers(&thegrid,size,::rank, 0);
    thegrid.mesh->RemoveGhost();
    thegrid.mesh->Redistribute(); 
    thegrid.mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
    //thegrid.mesh->AssignGlobalID(CELL | EDGE | FACE | NODE);
}

/// Redistribute grid by  partitioner
void redistribute(int type)
{
	//std::stringstream str1;
	//str1 << "startredist" << ::rank << ".xml";
	//thegrid.mesh->Save(str1.str());
	//thegrid.mesh->Save("startredist.pmf");
	thegrid.mesh->ResolveShared();
	//thegrid.mesh->Save("resolveshared.pmf");
	//std::stringstream str2;
	//str2 << "resolveshared" << ::rank << ".xml";
	//thegrid.mesh->Save(str2.str());
	LOG(2,"Process " << ::rank << ": redistribute. Cells: " << thegrid.mesh->NumberOfCells())
	if( type == 3 )
	{
		TagInteger r = thegrid.mesh->RedistributeTag();
		TagInteger o = thegrid.mesh->OwnerTag();
		for(Mesh::iteratorCell it = thegrid.mesh->BeginCell(); it != thegrid.mesh->EndCell(); ++it)
			r[it->self()] = (o[it->self()]+1)%size;
	}
	else
	{
    Partitioner * part = new Partitioner(thegrid.mesh);
    
    // Specify the partitioner
    if (type == 0) part->SetMethod(Partitioner::Parmetis, Partitioner::Partition);
    if (type == 1) part->SetMethod(Partitioner::Parmetis, Partitioner::Repartition);
    if (type == 2) part->SetMethod(Partitioner::Parmetis, Partitioner::Refine);
    try
    {
        part->Evaluate();
    }
    catch (INMOST::ErrorType er)
    {
        cout << "Exception: " << er << endl;
    }
    delete part;
	}

	//thegrid.mesh->Save("parmetis.pmf");

	correct_brothers(&thegrid,size,::rank, 2);

	//thegrid.mesh->Save("brothers.pmf");

    try
    {
        thegrid.mesh->Redistribute(); 
    }
    catch (INMOST::ErrorType er)
    {
        cout << "Exception: " << er << endl;
    }

	//thegrid.mesh->Save("redistribute.pmf");

    thegrid.mesh->RemoveGhost();

	//thegrid.mesh->Save("removeghost.pmf");

    thegrid.mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);

	//thegrid.mesh->Save("reorderempty.pmf");

    thegrid.mesh->AssignGlobalID(CELL | EDGE | FACE | NODE);

	//thegrid.mesh->Save("assigngid.pmf");

	LOG(2,"Process " << ::rank << ": redistribute completed")
}


int main(int argc, char ** argv)
{
	int i;
	int n[3] = {20,40,1};

	thegrid.transformation = transformation;
	thegrid.rev_transformation = rev_transformation;
	thegrid.cell_should_split = cell_should_split;
	thegrid.cell_should_unite = cell_should_unite;
    Mesh::Initialize(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &::rank);

    gridInit(&thegrid,n);

    Partitioner::Initialize(&argc,&argv);

	::size = thegrid.mesh->GetProcessorsNumber();
	::rank = thegrid.mesh->GetProcessorRank();
	std::cout << "rank " << ::rank << " size " << ::size << std::endl;
    //dump_to_vtk();
	std::cout << "gridAMR" << std::endl;
	gridAMR(&thegrid, 0);
	std::cout << "redistribute" << std::endl;
	redistribute(1);
	std::cout << "after redistribute " << thegrid.mesh->NumberOfCells() <<  std::endl;
	thegrid.mesh->Save("end.pvtk");
	Partitioner::Finalize();
	Mesh::Finalize();
	return 0;
	
}
