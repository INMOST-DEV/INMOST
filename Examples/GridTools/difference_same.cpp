#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace INMOST;


int main(int argc, char *argv[]) 
{
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] ;
    std::cout << " mesh1 mesh2" << std::endl;
    return -1;
  }
  Mesh m1,m2;
	//m1.SetFileOption("VTK_GRID_DIMS","2");
	//m2.SetFileOption("VTK_GRID_DIMS","2");
	m1.SetFileOption("VERBOSITY","2");
	m2.SetFileOption("VERBOSITY","2");
  m1.Load(argv[1]);
  m2.Load(argv[2]);
	Tag volumes_tag = m1.CreateTag("VOLUMES_FOR_NORMS_CALCULATION",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,1);
	for(ElementType etype = NODE; etype <= CELL; etype = etype << 1)
	{
		for(int id = 0; id < m1.LastLocalID(etype); id++)
		{
			Element c1 = m1.ElementByLocalID(etype,id);
			if( c1.isValid() )
			{
				Storage::real vol = 0;
				ElementArray<Cell> cells = c1.getCells();
				for(ElementArray<Cell>::iterator it = cells.begin(); it != cells.end(); ++it)
				{
					vol += it->Volume() / static_cast<Storage::real>(it->nbAdjElements(etype));
				}
				c1->RealDF(volumes_tag) = vol;
			}
		}
	}

  for(Mesh::iteratorTag t = m1.BeginTag(); t != m1.EndTag(); t++)
  {
    if( *t == m1.CoordsTag() ) continue;
    if( t->GetSize() == ENUMUNDEF ) continue;
    if( t->GetDataType() != DATA_REAL ) continue;
    if( m2.HaveTag(t->GetTagName()) )
    {
      Tag t2 = m2.GetTag(t->GetTagName());

      if( t->GetSize() != t2.GetSize() ) continue;

      for(ElementType etype = NODE; etype <= CELL; etype = etype << 1)
      {
				if( m1.ElementDefined(*t,etype) && m2.ElementDefined(t2,etype) )
        {
					Storage::real Cnorm = 0, L1norm = 0, L2norm = 0, absval, Lvol = 0, vol;
          for(int id = 0; id < m1.LastLocalID(etype); id++)
          {
						Element c1 = m1.ElementByLocalID(etype,id);
            Element c2 = m2.ElementByLocalID(etype,id);
            if( c1.isValid() && c2.isValid() )
            {
							vol = c1->RealDF(volumes_tag);
              Storage::real_array arr1 = c1->RealArray(*t);
              Storage::real_array arr2 = c2->RealArray(t2);
              for(int k = 0; k < (int)arr1.size(); k++)
							{
								absval = fabs(arr1[k]-arr2[k]);

                arr1[k] = absval;
								

								if( Cnorm < absval ) Cnorm = absval;
								L1norm += absval*vol;
								L2norm += absval*absval*vol;
								Lvol += vol;
							}
            }
          }
					std::cout << "Data " << t->GetTagName() << " on " << ElementTypeName(etype) << " Cnorm " << Cnorm << " L1norm " << L1norm/Lvol << " L2norm " << sqrt(L2norm/Lvol) << std::endl;
        }
      }
    }
  }

	m1.DeleteTag(volumes_tag);
  m1.Save("diff.gmv");
  m1.Save("diff.vtk");
  m1.Save("diff.pmf");

  return 0;
}
