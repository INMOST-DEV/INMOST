#include <cstdio>
#include <cmath>

#include "inmost.h"
using namespace INMOST;

int main(int argc,char ** argv)
{
	int errors = 0;
	Mesh::Initialize(&argc,&argv);
	Mesh * m = new Mesh(); // Create an empty mesh
	m->SetCommunicator(INMOST_MPI_COMM_WORLD); // Set the MPI communicator for the mesh
	int rank = m->GetProcessorRank(),  nproc = m->GetProcessorsNumber();

	/* Expected serial mesh with IxJxK box grid, you can try different prime values, like 13x17x19 */
	if( rank == 0 ) m->Load((argc>1)?argv[1]:"3x3x3.vtk");

	
	int total = m->TotalNumberOf(CELL);
	if( rank == 0 ) std::cout << "Total cells: " << total << std::endl;

	// Generate GlobalID for cells
	m->AssignGlobalID(CELL);
	// Copy them in order to check later
	Tag idcopy = m->CreateTag("ID_copy",DATA_INTEGER,CELL,NONE,1);
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
	    it->Integer(idcopy) = it->GlobalID();
	}

	// Create Tag to verify RedistributeTag
	Tag pccopy = m->CreateTag("PC_copy",DATA_INTEGER,CELL,NONE,1);
	for(int permut = 0; permut < nproc; permut++)
	{
		// Redistribute mesh in chess-board fashion, shifting each time by one proccess
		Tag redist = m->RedistributeTag();
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		{
			it->Integer(redist) = (it->GlobalID() + permut) % nproc;
			it->Integer(pccopy) = it->Integer(redist);
		}
		m->Redistribute(); // Redistribute the mesh data
		m->ReorderEmpty(CELL|FACE|EDGE|NODE); // Clean the data after reordring [optional]

		// Report number of local cells
		std::cout << "Permut: " << permut << ", proc: " << rank << ", cells: " << m->NumberOfCells() << std::endl;

		// Check that GlobalID and RedistributeTag are preserved 
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		{
			if( it->Integer(pccopy) != rank )
			{
				std::cout << "Permut: " << permut << ", proc: " << rank << ", cell: " << it->GlobalID() << ", pccopy: " << it->Integer(pccopy) << std::endl;
				errors++;
			}
			if( it->Integer(idcopy) != it->GlobalID() )
			{
				std::cout << "Permut: " << permut << ", proc: " << rank << ", cell: " << it->GlobalID() << ", idcopy: " << it->Integer(idcopy) << std::endl;
				errors++;
			}
		}
	}

	// Obtain 1 layer of ghost cells
	m->ExchangeGhost(1,FACE);
	std::cout << "Ghost: proc: " << rank << ", cells: " << m->NumberOfCells() << std::endl;

	if( nproc == 2 )  // In case np=2 and simple box 3x3x3 these should cover the whole mesh, let's check that!
	{
		int * mask = new int[total];
		for (int i=0; i<total; i++)  mask[i] = 0;
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		{
		    mask[it->GlobalID()]++;
		}
		for (int i=0; i<total; i++) if ( mask[i] != 1)
		{
			std::cout << "Check: proc: " << rank << ", cell: " << i << ((mask[i])? " duplicated" : " missing") << std::endl;
			errors++;
		}
		delete [] mask;
	}
  delete m;
  Mesh::Finalize();
  return (errors)? -1 : 0;
}
