#include "inmost.h"
using namespace INMOST;


// Compute number of connected components
int components(Mesh *m, Tag t)
{
	if( m->NumberOfCells() == 0 )  return 0;

	// We use marker to mark unknown cells;
	MarkerType unknown = m->CreateMarker();
	if( unknown == InvalidMarker() )  return -1;

	// Last-In-First-Out structure for searching connected components
	ElementArray<Cell> lifo(m); 
	Cell lastcell;
	int comp = 0, color = m->GetProcessorRank();

	// Mark all cells
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		it->SetMarker(unknown);
		lastcell = it->getAsCell();
	}
	// Start with one cell
	lifo.push_back(lastcell);

	while( lifo.size() > 0 )
	{
		while( lifo.size() > 0 )
		{
			// Get last cell
			lastcell = lifo.back();
			lastcell->Integer(t) = color;
			lifo.pop_back();
			// Unmark cell
			lastcell->RemMarker(unknown);
			// Add all unknown neighbours
			ElementArray<Cell> neighbours = lastcell->BridgeAdjacencies2Cell(FACE,unknown);
			for(ElementArray<Cell>::iterator it = neighbours.begin(); it != neighbours.end(); ++it)
			{
				lifo.push_back(it->getAsCell());
			}
		}
		// Connected component is closed
		comp++;
		color+=m->GetProcessorsNumber();
		// Check if there are still unknown cells
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		{
			if( it->GetMarker(unknown) )
			{
				lastcell = it->getAsCell();
				break;
			}
		}
		if( lastcell->GetMarker(unknown) )  lifo.push_back(lastcell);

	}
	m->ReleaseMarker(unknown);

	return comp;
}

int main(int argc,char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " input.[p]vtk type [action]" << std::endl << std::endl
		          << "type: " << std::endl
		          << "0 - Inner_RCM" << std::endl
		          << "1 - Parmetis" << std::endl
		          << "2 - Zoltan_HSFC" << std::endl
		          << "3 - Zoltan_RIB" << std::endl
		          << "4 - Zoltan_RCB" << std::endl
		          << "5 - Zoltan_PHG" << std::endl
		          << "6 - Zoltan_Scotch" << std::endl
		          << "7 - Zoltan_Parmetis" << std::endl
		          << "-1 - No partitioner" << std::endl << std::endl
		          << "action (for Parmetis): " << std::endl
		          << "0 - Partition" << std::endl
		          << "1 - Repartition" << std::endl
		          << "2 - Refine" << std::endl;
		return 0;
	}

	Mesh::Initialize(&argc,&argv);
	Partitioner::Initialize(&argc,&argv);

	Mesh * m = new Mesh(); // Create an empty mesh
	m->SetCommunicator(INMOST_MPI_COMM_WORLD); // Set the MPI communicator for the mesh
	int rank = m->GetProcessorRank();//,  nproc = m->GetProcessorsNumber();

	if( rank == 0 || m->isParallelFileFormat(argv[1]) )  m->Load(argv[1]); // Load input mesh

	int itype = atoi(argv[2]), iaction = 0;
	
	if( argc > 3 )  iaction = atoi(argv[3]);
	
	Partitioner::Type type;
	switch(itype)
	{
		case 0: type = Partitioner::Inner_RCM; break;
		case 1: type = Partitioner::Parmetis; break;
		case 2: type = Partitioner::Zoltan_HSFC; break;
		case 3: type = Partitioner::Zoltan_RIB; break;
		case 4: type = Partitioner::Zoltan_RCB; break;
		case 5: type = Partitioner::Zoltan_PHG; break;
		case 6: type = Partitioner::Zoltan_Scotch; break;
		case 7: type = Partitioner::Zoltan_Parmetis; break;
		case 8: /*reserved*/ break;
	}

	Partitioner::Action action;
	switch(iaction)
	{
		case 0: action = Partitioner::Partition; break;
		case 1: action = Partitioner::Repartition; break;
		case 2: action = Partitioner::Refine; break;
	}

	if( itype == 8 ) //pass the whole local mesh on the current processor to the next one
	{
		if( m->OwnerTag().isValid() )
		{
			TagInteger r = m->RedistributeTag(); //new location of the element
			TagInteger o = m->OwnerTag(); //current owner
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) //for each cell
				r[it->self()] = (o[it->self()]+1)%m->GetProcessorsNumber(); //pass to the next
			m->Redistribute(); // Redistribute the mesh data
			m->ReorderEmpty(CELL|FACE|EDGE|NODE); // Clean the data after reordring
		}
		else std::cout << "Current mesh is not parallel for this type of distribution" << std::endl;
	}
	else if( itype >= 0 )
	{
		Partitioner * p = new Partitioner(m);
		p->SetMethod(type,action); // Specify the partitioner
		p->Evaluate(); // Compute the partitioner and store new processor ID in the mesh
		delete p;
		m->RemoveGhost(); // Delete all ghost cells
		m->Redistribute(); // Redistribute the mesh data
		m->ReorderEmpty(CELL|FACE|EDGE|NODE); // Clean the data after reordring
	}

	Tag c = m->CreateTag("Comp",DATA_INTEGER,CELL,NONE,1);

	std::cout << "Proc " << rank << ": " << m->NumberOfCells() << " cells, " << components(m,c) << " components" << std::endl;
	int maxcells = m->AggregateMax(m->NumberOfCells());
	if( rank == 0 )  std::cout << "Max number of cells per proc: " << maxcells << std::endl;

	// m->Save("dump.pvtk"); // Save mesh for visualization

	delete m;
	Partitioner::Finalize();
	Mesh::Finalize();
	return 0;
}
