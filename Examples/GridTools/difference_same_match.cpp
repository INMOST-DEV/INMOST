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
	//m1.SetFileOption("VERBOSITY","2");
	//m2.SetFileOption("VERBOSITY","2");
	m1.Load(argv[1]);
	m2.Load(argv[2]);
	
	if( m1.NumberOfNodes() != m2.NumberOfNodes() )
	{
		std::cout << "Unequal number of nodes!" << std::endl;
		return -1;
	}
	
	TagReal volumes_tag = m1.CreateTag("VOLUMES_FOR_NORMS_CALCULATION",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,1);
	
	Tag ref12 = m1.CreateTag("M1TOM2",DATA_REMOTE_REFERENCE,NODE|EDGE|FACE|CELL,NONE,1);
	Tag ref21 = m2.CreateTag("M2TOM1",DATA_REMOTE_REFERENCE,NODE|EDGE|FACE|CELL,NONE,1);
	
	int nnodes = (int)m1.NumberOfNodes();
	
	
	
	std::vector<HandleType> nodes1(nnodes), nodes2(nnodes);
	int k;
	k = 0;
	for(Mesh::iteratorNode it = m1.BeginNode(); it != m1.EndNode(); ++it)
		nodes1[k++] = it->GetHandle();
	std::sort(nodes1.begin(),nodes1.end(),Mesh::CentroidComparator(&m1));
	k = 0;
	for(Mesh::iteratorNode it = m2.BeginNode(); it != m2.EndNode(); ++it)
		nodes2[k++] = it->GetHandle();
	std::sort(nodes2.begin(),nodes2.end(),Mesh::CentroidComparator(&m2));
	
	TagRealArray coords1(m1.CoordsTag()), coords2(m2.CoordsTag());
	
	int nmiss = 0;
	for(k = 0; k < nnodes; ++k)
	{
		if( (coords1(nodes1[k],3,1) - coords2(nodes2[k],3,1)).FrobeniusNorm() < 1.0e-8 )
		{
			m1.RemoteReference(nodes1[k],ref12) = RemoteHandleType(&m2,nodes2[k]);
			m2.RemoteReference(nodes2[k],ref21) = RemoteHandleType(&m1,nodes1[k]);
		}
		else
		{
			std::cout << "cannot match node " << k << " coords" << std::endl;
			std::cout << "at m1 "; coords1(nodes1[k],1,3).Print();
			std::cout << "at m2 "; coords1(nodes1[k],1,3).Print();
			nmiss++;
		}
	}
	if( nmiss )
		std::cout << "cannot match " << nmiss << " of NODE" << std::endl;
	std::cout << "matching done for NODE" << std::endl;
	
	std::vector<HandleType> remote_adj;
	
	for(ElementType etype = EDGE; etype <= CELL; etype = NextElementType(etype) )
	{
		nmiss = 0;
		for(Mesh::iteratorElement it = m1.BeginElement(etype); it != m1.EndElement(); ++it)
		{
			ElementArray<Element> adj = it->getAdjElements(PrevElementType(etype));
			for(ElementArray<Element>::iterator jt = adj.begin(); jt != adj.end(); ++jt)
				remote_adj.push_back(jt->RemoteReference(ref12).second);
			HandleType remote = m2.FindSharedAdjacency(&remote_adj[0],remote_adj.size());
			if( remote != InvalidHandle() )
				it->RemoteReference(ref12) = RemoteHandleType(&m2,remote);
			else
			{
				std::cout << "cannot match " << ElementTypeName(etype) << ":" << it->LocalID() << std::endl;
				nmiss++;
			}
			remote_adj.clear();
		}
		if( nmiss )
			std::cout << "cannot match " << nmiss << " of " << ElementTypeName(etype) << " from m1 to m2 " << std::endl;
		nmiss = 0;
		for(Mesh::iteratorElement it = m2.BeginElement(etype); it != m2.EndElement(); ++it)
		{
			ElementArray<Element> adj = it->getAdjElements(PrevElementType(etype));
			for(ElementArray<Element>::iterator jt = adj.begin(); jt != adj.end(); ++jt)
				remote_adj.push_back(jt->RemoteReference(ref21).second);
			HandleType remote = m1.FindSharedAdjacency(&remote_adj[0],remote_adj.size());
			if( remote != InvalidHandle() )
				it->RemoteReference(ref21) = RemoteHandleType(&m1,remote);
			else
			{
				std::cout << "cannot match " << ElementTypeName(etype) << ":" << it->LocalID() << std::endl;
				nmiss++;
			}
			remote_adj.clear();
		}
		if( nmiss )
			std::cout << "cannot match " << nmiss << " of " << ElementTypeName(etype) << " from m2 to m1 " << std::endl;
		std::cout << "matching done for " << ElementTypeName(etype) << std::endl;
	}
		
	/*
	for(ElementType etype = NODE; etype <= CELL; etype = etype << 1)
	{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
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
		std::cout << "volumes precomputed for " << ElementTypeName(etype) << std::endl;
	}
	*/
	for(int id = 0; id < m1.CellLastLocalID(); id++)
	{
		Cell c1 = m1.CellByLocalID(id);
		if( c1.isValid() )
		{
			double vol = c1.Volume();
			volumes_tag[c1] = vol;
			
			for(ElementType etype = NODE; etype <= FACE; etype = NextElementType(etype) )
			{
				ElementArray<Element> adj = c1.getAdjElements(etype);
				for(ElementArray<Element>::iterator it = adj.begin(); it != adj.end(); ++it)
					volumes_tag[*it] += vol / adj.size();
			}
		}
	}
	std::cout << "volumes precomputed" << std::endl;
	
	for(Mesh::iteratorTag t = m1.BeginTag(); t != m1.EndTag(); t++)
	{
		//if( *t == m1.CoordsTag() ) continue;
		//if( t->GetSize() == ENUMUNDEF ) continue;
		if( t->GetDataType() != DATA_REAL && t->GetDataType() != DATA_VARIABLE ) continue;
		if( m2.HaveTag(t->GetTagName()) )
		{
			Tag t2 = m2.GetTag(t->GetTagName());
			
			if( t->GetSize() != t2.GetSize() ) continue;
			
			for(ElementType etype = NODE; etype <= CELL; etype = NextElementType(etype))
			{
				if( m1.ElementDefined(*t,etype) && m2.ElementDefined(t2,etype) )
				{
					Storage::real Cnorm = 0, L1norm = 0, L2norm = 0, absval, Lvol = 0, vol;
					Storage::real max_val = -1.0e+20, min_val = 1.0e+20;
					int err = 0, tot = 0, diff_size = 0, tot_size = 0;
					for(int id = 0; id < m1.LastLocalID(etype); id++)
					{
						Element c1 = m1.ElementByLocalID(etype,id);
						Element c2 = MakeElement(c1.RemoteReference(ref12));
						if( c1.isValid() && c2.isValid() )
						{
							vol = c1->RealDF(volumes_tag);
							if( t->GetDataType() == DATA_REAL )
							{
								Storage::real_array arr1 = c1->RealArray(*t);
								Storage::real_array arr2 = c2->RealArray(t2);
								tot_size++;
								if( arr1.size() != arr2.size() )
								{
									diff_size++;
									continue;
									std::cout << "tag " << t->GetTagName() << " position " << k << "/" << arr1.size();
									std::cout << " arrays of different size ";
									std::cout << " on MESH1:" << ElementTypeName(etype) << ":" << c1.LocalID() << " " << arr1.size();
									std::cout << " on MESH2:" << ElementTypeName(etype) << ":" << c2.LocalID() << " " << arr2.size();
									std::cout << std::endl;
									continue;
								}
								for(int k = 0; k < (int)arr1.size(); k++)
								{
									max_val = std::max(max_val,arr1[k]);
									max_val = std::max(max_val,arr2[k]);
									min_val = std::min(min_val,arr1[k]);
									min_val = std::min(min_val,arr2[k]);

									absval = fabs(arr1[k]-arr2[k]);
										
									if( absval > 1.0e-3 )
									{
										if( true )
										{
											std::cout << "tag " << t->GetTagName() << " position " << k << "/" << arr1.size();
											std::cout << " on MESH1:" << ElementTypeName(etype) << ":" << c1.LocalID() << " " << arr1[k];
											std::cout << " on MESH2:" << ElementTypeName(etype) << ":" << c2.LocalID() << " " << arr2[k];
											std::cout << std::endl;
										}
										err++;
									}
									
									if( *t != m1.CoordsTag() )
										arr1[k] = absval;
									
									tot++;
									
									if( Cnorm < absval ) Cnorm = absval;
									L1norm += absval*vol;
									L2norm += absval*absval*vol;
									Lvol += vol;
								}
							}
							else if( t->GetDataType() == DATA_VARIABLE )
							{
								Storage::var_array arr1 = c1->VariableArray(*t);
								Storage::var_array arr2 = c2->VariableArray(t2);
								tot_size++;
								if( arr1.size() != arr2.size() )
								{
									diff_size++;
									continue;
								}
								for(int k = 0; k < (int)arr1.size(); k++)
								{
									max_val = std::max(max_val,get_value(arr1[k]));
									max_val = std::max(max_val,get_value(arr2[k]));
									min_val = std::min(min_val,get_value(arr1[k]));
									min_val = std::min(min_val,get_value(arr2[k]));


									absval = fabs(get_value(arr1[k])-get_value(arr2[k]));
									
									arr1[k] = absval;
									
									if( absval > 1.0e-3 ) err++;
									
									if( Cnorm < absval ) Cnorm = absval;
									L1norm += absval*vol;
									L2norm += absval*absval*vol;
									Lvol += vol;
								}
							}
						}
					}
					if( diff_size ) std::cout << "Size is different on " << diff_size << " / " << tot_size << " of " << ElementTypeName(etype) << std::endl;
					if( err ) std::cout << "Reported error on " << err << " / " << tot << " of " << ElementTypeName(etype) << std::endl;
					std::cout << "Data " << t->GetTagName() << " on " << ElementTypeName(etype) << " Cnorm " << Cnorm << " L1norm " << L1norm/Lvol << " L2norm " << sqrt(L2norm/Lvol) << " Min " << min_val << " Max " << max_val << std::endl;
				}
			}
		}
	}
	m1.DeleteTag(ref12);
	m1.DeleteTag(volumes_tag);
	m1.Save("diff.gmv");
	m1.Save("diff.vtk");
	m1.Save("diff.pmf");
	
	return 0;
}
