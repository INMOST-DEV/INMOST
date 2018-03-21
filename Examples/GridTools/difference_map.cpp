#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "inmost.h"

using namespace INMOST;

class kdtree_mapping
{
public:
	
	inline static bool cell_point(const Cell & c, const Storage::real p[3])
	{
		return c.Inside(p);
	}
	template<typename bbox_type>
	inline static int bbox_point(const Storage::real p[3], const bbox_type bbox[6])
	{
		for(int i = 0; i < 3; i++)
		{
			if( p[i] < bbox[i*2] || p[i] > bbox[i*2+1] )
				return 0;
		}
		return 1;
	}
private:
	struct entry
	{
		HandleType e;
		float xyz[3];
		struct entry & operator =(const struct entry & other)
		{
			e = other.e;
			xyz[0] = other.xyz[0];
			xyz[1] = other.xyz[1];
			xyz[2] = other.xyz[2];
			return *this;
		}
	} * set;
	Mesh * m;
	INMOST_DATA_ENUM_TYPE size;
	float bbox[6];
	kdtree_mapping * children;
	static int cmpElements0(const void * a,const void * b) 
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[0];
		float bd = eb->xyz[0];
		return (ad > bd) - (ad < bd);
	}
	static int cmpElements1(const void * a,const void * b) 
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[1];
		float bd = eb->xyz[1];
		return (ad > bd) - (ad < bd);
	}
	static int cmpElements2(const void * a,const void * b) 
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[2];
		float bd = eb->xyz[2];
		return (ad > bd) - (ad < bd);
	}
	inline static unsigned int flip(const unsigned int * fp)
	{
		unsigned int mask = -((int)(*fp >> 31)) | 0x80000000;
		return *fp ^ mask;
	}
#define _0(x)	(x & 0x7FF)
#define _1(x)	(x >> 11 & 0x7FF)
#define _2(x)	(x >> 22 )
	void radix_sort(int dim, struct entry * temp)
	{
		unsigned int i;
		const unsigned int kHist = 2048;
		unsigned int  b0[kHist * 3];
		unsigned int *b1 = b0 + kHist;
		unsigned int *b2 = b1 + kHist;
		memset(b0,0,sizeof(unsigned int)*kHist*3);
		for (i = 0; i < size; i++) 
		{
			unsigned int fi = flip((unsigned int *)&set[i].xyz[dim]);
			++b0[_0(fi)]; ++b1[_1(fi)]; ++b2[_2(fi)];
		}
		{
			unsigned int sum0 = 0, sum1 = 0, sum2 = 0;
			for (i = 0; i < kHist; i++) 
			{
				b0[kHist-1] = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = b0[kHist-1];
				b1[kHist-1] = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = b1[kHist-1];
				b2[kHist-1] = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = b2[kHist-1];
			}
		}
		for (i = 0; i < size; i++) temp[++b0[_0(flip((unsigned int *)&set[i].xyz[dim]))]] = set[i];
		for (i = 0; i < size; i++) set[++b1[_1(flip((unsigned int *)&temp[i].xyz[dim]))]] = temp[i];
		for (i = 0; i < size; i++) temp[++b2[_2(flip((unsigned int *)&set[i].xyz[dim]))]] = set[i];
		for (i = 0; i < size; i++) set[i] = temp[i];
	}
	void kdtree_build(int dim, int & done, int total, struct entry * temp)
	{
		if( size > 1 )
		{
			if( size > 128 ) radix_sort(dim,temp); else 
			switch(dim)
			{
			case 0: qsort(set,size,sizeof(entry),cmpElements0);break;
			case 1: qsort(set,size,sizeof(entry),cmpElements1);break;
			case 2: qsort(set,size,sizeof(entry),cmpElements2);break;
			}
			children = static_cast<kdtree_mapping *>(malloc(sizeof(kdtree_mapping)*2));//new kdtree[2];
			children[0].children = NULL;
			children[0].set = set;
			children[0].size = size/2;
			children[0].m = m;
			children[1].children = NULL;
			children[1].set = set+size/2;
			children[1].size = size - size/2;
			children[1].m = m;
			children[0].kdtree_build((dim+1)%3,done,total,temp);
			children[1].kdtree_build((dim+1)%3,done,total,temp);
			for(int k = 0; k < 3; k++)
			{
				bbox[0+2*k] = std::min<float>(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
				bbox[1+2*k] = std::max<float>(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
			}
		}
		else 
		{
			assert(size == 1);
			{
				ElementArray<Node> nodes = Element(m,set[0].e)->getNodes();
				bbox[0] = bbox[2] = bbox[4] = 1.0e20f;
				bbox[1] = bbox[3] = bbox[5] = -1.0e20f;
				for(ElementArray<Node>::size_type k = 0; k < nodes.size(); ++k)
				{
					Storage::real_array coords = nodes[k].Coords();
					for(Storage::real_array::size_type q = 0; q < coords.size(); q++)
					{
						bbox[q*2+0] = std::min<float>(bbox[q*2+0],(float)coords[q]);
						bbox[q*2+1] = std::max<float>(bbox[q*2+1],(float)coords[q]);
					}
				}
				for(int k = m->GetDimensions(); k < 3; ++k)
				{
					bbox[k*2+0] = -1.0e20f;
					bbox[k*2+1] = 1.0e20f;
				}
			}
		}
	}
	kdtree_mapping() : set(NULL), size(0), children(NULL) {}
	
	Cell sub_cell_point(const Storage::real p[3])
	{
		Cell ret = InvalidCell();
		if( size == 1 )
		{
			if( cell_point(Cell(m,set[0].e),p) )
				ret = Cell(m,set[0].e);
		}
		else 
		{
			assert(size > 1);
			if( bbox_point(p,bbox) )
			{
				ret = children[0].sub_cell_point(p);
				if( !ret.isValid() )
					ret = children[1].sub_cell_point(p);
			}
		}
		return ret;
	}
	void clear_children() { if( children ) {children[0].clear_children(); children[1].clear_children(); free(children);}}
public:
	kdtree_mapping(Mesh * m) : m(m), children(NULL)
	{
		size = m->NumberOfCells();
		assert(size > 1);
		set = new entry[size];
		INMOST_DATA_ENUM_TYPE k = 0;
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) 
		{
			set[k].e = it->GetHandle();
			Storage::real cnt[3] = {0,0,0};
			m->GetGeometricData(set[k].e,CENTROID,cnt);
			set[k].xyz[0] = (float)cnt[0];
			set[k].xyz[1] = (float)cnt[1];
			set[k].xyz[2] = (float)cnt[2];
			++k;
		}
		int done = 0, total = size;
		struct entry *  temp = new entry[size];
		kdtree_build(0,done,total,temp);
		delete [] temp;
		for(int k = 0; k < 3; k++)
		{
			bbox[0+2*k] = std::min<float>(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
			bbox[1+2*k] = std::max<float>(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
		}
	}
	
	~kdtree_mapping()
	{
		delete [] set;
		clear_children();
	}

	Cell cell_point(const Storage::real * point)
	{
		return sub_cell_point(point);
	}
};


int main(int argc, char *argv[]) 
{
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] ;
    std::cout << " mesh1 mesh2 [out_mesh1] [out_mesh2]" << std::endl;
    return -1;
  }
  Mesh m1,m2;
	//m1.SetFileOption("VTK_GRID_DIMS","2");
	//m2.SetFileOption("VTK_GRID_DIMS","2");
	m1.SetFileOption("VERBOSITY","2");
	m2.SetFileOption("VERBOSITY","2");
  m1.Load(argv[1]);
  m2.Load(argv[2]);

	Tag map_m2 = m1.CreateTag("MAPPING_TO_M2",DATA_REFERENCE,CELL,NONE);
	Tag map_m1 = m2.CreateTag("MAPPING_TO_M1",DATA_REFERENCE,CELL,NONE);
	//prepare mapping between grids
	{
		kdtree_mapping tm1(&m1), tm2(&m2);
		Storage::real cnt[3] = {0,0,0};
		int k = 0;
		for(Mesh::iteratorCell it = m1.BeginCell(); it != m1.EndCell(); ++it)
		{
			it->Centroid(cnt);
			Cell c = tm2.cell_point(cnt); // c belongs to m2
			if( c.isValid() )
			{
				c->ReferenceArrayDV(map_m1).push_back(it->GetHandle());
				it->ReferenceArrayDV(map_m2).push_back(c->GetHandle());
			}
			else
			{
				std::cout << "Cannot map cell " << k << " of first mesh at " << cnt[0] << " " << cnt[1] << " " << cnt[2] << std::endl;
			}
			++k;
		}
		k = 0;
		for(Mesh::iteratorCell it = m2.BeginCell(); it != m2.EndCell(); ++it)
		{
			it->Centroid(cnt);
			Cell c = tm1.cell_point(cnt); // c belongs to m1
			if( c.isValid() )
			{
				c->ReferenceArrayDV(map_m2).push_back(it->GetHandle());
				it->ReferenceArrayDV(map_m1).push_back(c->GetHandle());
			}
			else
			{
				std::cout << "Cannot map cell " << k << " of second mesh at " << cnt[0] << " " << cnt[1] << " " << cnt[2] << std::endl;
			}
			++k;
		}
		//remove duplicates
		k = 0;
		for(Mesh::iteratorCell it = m1.BeginCell(); it != m1.EndCell(); ++it)
		{
			Storage::reference_array arr = it->ReferenceArrayDV(map_m2);
			HandleType * beg = arr.data(), * end = arr.data()+arr.size(), * shr;
			std::sort(beg,end);
			//remove all invalid entries
			shr = beg;
			while(*shr == 0 && shr < end) ++shr;
			if( shr != beg ) 
			{
				memmove(beg,shr,(end-shr)*sizeof(HandleType));
				end = beg + (end-shr);
			}
			//remove duplicate entries
			arr.resize(std::unique(beg,end)-beg);
			if( arr.empty() ) std::cout << "No cells for cell " << k << " from first mesh to second" << std::endl;
			++k;
		}
		k = 0;
		for(Mesh::iteratorCell it = m2.BeginCell(); it != m2.EndCell(); ++it)
		{
			Storage::reference_array arr = it->ReferenceArrayDV(map_m1);
			HandleType * beg = arr.data(), * end = arr.data()+arr.size(), * shr;
			//sort handles
			std::sort(beg,end);
			//remove all invalid entries
			shr = beg;
			while(*shr == 0 && shr < end) ++shr;
			if( shr != beg ) 
			{
				memmove(beg,shr,(end-shr)*sizeof(HandleType));
				end = beg + (end-shr);
			}
			//remove duplicate entries
			arr.resize(std::unique(beg,end)-beg);
			if( arr.empty() ) std::cout << "No cells for cell " << k << " from second mesh to first" << std::endl;
			++k;
		}
	}
	//prepare volumes for norms
	Mesh::GeomParam params;
	params[ORIENTATION] = FACE;
	params[MEASURE] = CELL;
	params[NORMAL] = FACE;
	params[CENTROID] = CELL | FACE;
	params[BARYCENTER] = CELL | FACE;
	m1.RemoveGeometricData(params);
	m1.PrepareGeometricData(params);
	m2.RemoveGeometricData(params);
	m2.PrepareGeometricData(params);

	//for test
	Tag one_m1 = m1.CreateTag("SELF_TEST_ONES",DATA_REAL,CELL,NONE,1);
	Tag one_m2 = m2.CreateTag("SELF_TEST_ONES",DATA_REAL,CELL,NONE,1);

	for(Mesh::iteratorCell it = m1.BeginCell(); it != m1.EndCell(); ++it) it->RealDF(one_m1) = 1.0;
	for(Mesh::iteratorCell it = m2.BeginCell(); it != m2.EndCell(); ++it) it->RealDF(one_m2) = 1.0;

	std::vector<Storage::real> difference;
	
  for(Mesh::iteratorTag t = m1.BeginTag(); t != m1.EndTag(); t++)
  {
    if( *t == m1.CoordsTag() ) continue;
    if( t->GetSize() == ENUMUNDEF ) continue;
    if( t->GetDataType() != DATA_REAL ) continue;
		if( t->GetTagName() == "PROTECTED_GEOM_UTIL_MEASURE" ) continue;
    if( m2.HaveTag(t->GetTagName()) )
    {
			Tag t1 = *t;
      Tag t2 = m2.GetTag(t->GetTagName());

      if( t1.GetSize() != t2.GetSize() ) continue;

			if( t1.GetTagName() == "SELF_TEST_ONES" )
			{
				std::cout << "Data with name SELF_TEST_ONES have value of one on all the cells." << std::endl;
				std::cout << "It is used here to perform self test of mapping consistancy." << std::endl;
				std::cout << "All the norms of this data must be zero if cells of the mesh" << std::endl;
				std::cout << "ideally fit each other. Otherwise they should be close to zero." << std::endl;
			}

      if( m1.ElementDefined(t1,CELL) && m2.ElementDefined(t2,CELL) )
      {
				difference.resize(m1.NumberOfCells()*t1.GetSize());

				Storage::real Cnorm = 0, L1norm = 0, L2norm = 0, absval, Lvol = 0, vol;
				Storage::real mean_map, vol_map, val_map;
				int k = 0;
        for(int id = 0; id < m1.CellLastLocalID(); ++id)
        {
          Element c1 = m1.CellByLocalID(id);
					if( c1.isValid() )
					{
						Storage::real_array arr1 = c1.RealArray(t1);
						for(Storage::real_array::size_type entry = 0; entry < arr1.size(); ++entry)
						{
							Storage::reference_array map = c1->ReferenceArrayDV(map_m2);
							map.SetMeshLink(&m2);
							mean_map = vol_map = 0.0;
							for(Storage::reference_array::iterator it = map.begin(); it != map.end(); ++it)
							{
								vol = it->getAsCell()->Volume();
								mean_map += it->RealArray(t2)[entry]*vol;
								vol_map += vol;
							}
							if( vol_map > 0.0 )
								val_map = mean_map/vol_map;
							else
							{
								std::cout << __FILE__ << ":" << __LINE__ << std::endl;
								std::cout << "element " << c1.LocalID() << std::endl;
								std::cout << "vol is " << vol_map << std::endl;
								std::cout << "map size " << map.size() << std::endl;
								for(Storage::reference_array::iterator it = map.begin(); it != map.end(); ++it)
									std::cout << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << " (" << it->getAsCell()->Volume() << ") ";
								std::cout << std::endl;
								
								val_map = 0.0;
							}
							absval = fabs(arr1[entry] - val_map);
							difference[k*t1.GetSize() + entry] = absval;

							vol = c1.getAsCell().Volume();
							if( Cnorm < absval ) Cnorm = absval;
							L1norm += absval*vol;
							L2norm += absval*absval*vol;
							Lvol += vol;
						}
						++k;
					}
        }

				std::cout << "Mesh1 Data " << t1.GetTagName() << " Cnorm " << Cnorm << " L1norm " << L1norm/Lvol << " L2norm " << sqrt(L2norm/Lvol) << std::endl;

				Cnorm = 0;
				L1norm = 0;
				L2norm = 0;
				Lvol = 0;

				for(int id = 0; id < m2.CellLastLocalID(); ++id)
				{
					Element c2 = m2.CellByLocalID(id);
					if( c2.isValid() )
					{
						Storage::real_array arr2 = c2.RealArray(t2);
						for(Storage::real_array::size_type entry = 0; entry < arr2.size(); ++entry)
						{
							Storage::reference_array map = c2->ReferenceArrayDV(map_m1);
							map.SetMeshLink(&m1);
							mean_map = vol_map = 0.0;
							for(Storage::reference_array::iterator it = map.begin(); it != map.end(); ++it)
							{
								vol = it->getAsCell()->Volume();
								mean_map += it->RealArray(t1)[entry]*vol;
								vol_map += vol;
							}
							if( vol_map > 0.0 )
								val_map = mean_map/vol_map;
							else
							{
								std::cout << __FILE__ << ":" << __LINE__ << std::endl;
								std::cout << "element " << c2.LocalID() << std::endl;
								std::cout << "vol is " << vol_map << std::endl;
								std::cout << "map size " << map.size() << std::endl;
								for(Storage::reference_array::iterator it = map.begin(); it != map.end(); ++it)
									std::cout << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << " ";
								std::cout << std::endl;
								val_map = 0.0;
							}
							absval = fabs(arr2[entry] - val_map);
							arr2[entry] = absval;

							vol = c2.getAsCell().Volume();
							if( Cnorm < absval ) Cnorm = absval;
							L1norm += absval*vol;
							L2norm += absval*absval*vol;
							Lvol += vol;
						}
					}
				}

				std::cout << "Mesh2 Data " << t2.GetTagName() << " Cnorm " << Cnorm << " L1norm " << L1norm/Lvol << " L2norm " << sqrt(L2norm/Lvol) << std::endl;

				k = 0;
				for(int id = 0; id < m1.CellLastLocalID(); ++id)
        {
          Element c1 = m1.CellByLocalID(id);
					if( c1.isValid() )
					{
						Storage::real_array arr1 = c1.RealArray(t1);
						for(Storage::real_array::size_type entry = 0; entry < arr1.size(); ++entry)
							arr1[entry] = difference[k*t1.GetSize()+entry];

						++k;
					}
				}
      }
    }
  }

	if( argc > 3 ) m1.Save(argv[3]); else m1.Save("mesh1diff.vtk");
	if( argc > 4 ) m2.Save(argv[4]); else m2.Save("mesh2diff.vtk");

  return 0;
}
