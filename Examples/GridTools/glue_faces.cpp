#include "inmost.h"
#include <stdio.h>
#include <deque>
//#include "tetgen/tetgen.h"
using namespace INMOST;
typedef Storage::real real;
typedef Storage::real_array real_array;



int main(int argc, char ** argv)
{
	if( argc < 2 ) 
	{
		printf("Usage: %s input_mesh [output_mesh]\n",argv[0]);
		return -1;
	}

	Mesh m;
	m.SetFileOption("ECL_CURVILINEAR","FALSE");
	m.SetFileOption("ECL_SPLIT_GLUED","TRUE");
	m.Load(argv[1]);
	Tag sliced;
	if( m.HaveTag("SLICED") ) sliced = m.GetTag("SLICED");
	
	ElementArray<Cell> unite(&m,2);
	
	m.BeginModification();
	
	//MarkerType mrk = m.CreateMarker();
	//Face last, next, test(&m,ComposeHandle(FACE,382));
	for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
	{
		//std::cout << "test face " << test.LocalID() << " back " << test.BackCell().LocalID() << " front " << test.FrontCell().LocalID() << std::endl;
		/*
		int a = last.isValid() ? last.LocalID() : -1;
		if( last.isValid() && next.isValid() )
		{
			//std::cout << "unite " << last.BackCell().LocalID();
			//std::cout << " (" << Element::GeometricTypeName(last.BackCell().GetGeometricType()) << ":" << last.BackCell().nbAdjElements(FACE) << ")";
			//std::cout << " and " << last.FrontCell().LocalID();
			//std::cout << " (" << Element::GeometricTypeName(last.FrontCell().GetGeometricType()) << ":" << last.FrontCell().nbAdjElements(FACE) << ")";
			//std::cout << " for " << last.LocalID();
			//std::cout << " next " << next.LocalID() << std::endl;
			unite[0] = last->BackCell();
			unite[1] = last->FrontCell();
			Cell c = Cell::UniteCells(unite,0);
			c.SetMarker(mrk);
			last = next;
			next = InvalidFace();
		}
		*/
		if( !it->Boundary() && (!sliced.isValid() || !it->HaveData(sliced)) )
		{
			//if( !it->BackCell().GetMarker(mrk) && !it->FrontCell().GetMarker(mrk) )
			if( !it->BackCell().New() && !it->FrontCell().New() )
			{
				
				//std::cout << "candidate " << it->LocalID();
				//std::cout << " back " << it->BackCell()->LocalID();
				//std::cout << " (" << Element::GeometricTypeName(it->BackCell().GetGeometricType()) <<  ":" << it->BackCell().nbAdjElements(FACE) <<")";
				//std::cout << " front " << it->FrontCell()->LocalID();
				//std::cout << " (" << Element::GeometricTypeName(it->FrontCell().GetGeometricType()) <<  ":" << it->BackCell().nbAdjElements(FACE) <<")";
				//std::cout << std::endl;
				/*
				it->BackCell().SetMarker(mrk);
				it->FrontCell().SetMarker(mrk);
				if( !last.isValid() )
					last = it->self();
				else
					next = it->self();
				 */
				//break;
				
				unite[0] = it->BackCell();
				unite[1] = it->FrontCell();
				Cell c = Cell::UniteCells(unite,0);
				
			}
		}
	}
	
	
	m.EndModification();
	
	TagInteger tetra = m.CreateTag("TET",DATA_INTEGER,CELL,NONE,1);
	
	int ntet = 0;
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		if( it->GetGeometricType() == Element::Tet )
		{
			tetra[*it] = 1;
			ElementArray<Cell> cells = it->BridgeAdjacencies2Cell(FACE);
			for(ElementArray<Cell>::size_type k = 0; k < cells.size(); ++k)
				if( cells[k]->GetGeometricType() == Element::Tet )
					tetra[*it] = 2;
			//std::cout << "Cell:" << it->LocalID() << " " << it->nbAdjElements(FACE) << std::endl;
			ntet ++;
		}
	
	std::cout << "tetrahedras: " << ntet << std::endl;
	
							  
	if( argc > 2 )
		m.Save(argv[2]);
	else
		m.Save("out.vtk");

	return 0;
}
