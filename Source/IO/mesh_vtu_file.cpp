#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif


#include "inmost.h"
#include <cfloat>

#if defined(USE_MESH)

static int __isnan__(double x) { return x != x; }
//static int isinf(double x) { return !isnan(x) && isnan(x - x); }
static int __isinf__(double x) { return fabs(x) > DBL_MAX; }
static int __isbad(double x) { return __isnan__(x) || __isinf__(x); }



namespace INMOST
{

	
	void Mesh::LoadVTU(std::string File)
	{
		int verbosity = 0;
		for (INMOST_DATA_ENUM_TYPE k = 0; k < file_options.size(); ++k)
		{
			if (file_options[k].first == "VERBOSITY")
			{
				verbosity = atoi(file_options[k].second.c_str());
				if (verbosity < 0 || verbosity > 2)
				{
					printf("%s:%d Unknown verbosity option: %s\n", __FILE__, __LINE__, file_options[k].second.c_str());
					verbosity = 1;
				}
			}
		}

		MarkerType unused_marker = CreateMarker();
		int grid_is_2d = 2;
		for (INMOST_DATA_ENUM_TYPE k = 0; k < file_options.size(); ++k)
		{
			if (file_options[k].first == "VTK_GRID_DIMS")
			{
				if (file_options[k].second == "AUTO")
					grid_is_2d = 2;
				if (atoi(file_options[k].second.c_str()) == 2)
					grid_is_2d = 1;
				else if (atoi(file_options[k].second.c_str()) == 3)
					grid_is_2d = 0;
			}
		}

		//Determine whether there are already 3d elements so that the grid is 3d
		if (grid_is_2d == 2 && NumberOfCells())
		{
			for (Mesh::iteratorCell it = BeginCell(); it != EndCell() && grid_is_2d == 2; ++it)
			if (it->GetElementDimension() == 3)
				grid_is_2d = 0;
		}


		std::vector<HandleType> old_nodes(NumberOfNodes());
		{
			unsigned qq = 0;
			for (Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it)
				old_nodes[qq++] = *it;
		}
		
		
		
		if (!old_nodes.empty())
		{
			std::sort(old_nodes.begin(), old_nodes.end(), CentroidComparator(this));
			//for(std::vector<HandleType>::iterator it = old_nodes.begin(); it != old_nodes.end(); ++it)
			//{
			//	Storage::real_array c = RealArrayDF(*it,CoordsTag());
			//	REPORT_VAL("coord: ",c[0] << " " << c[1] << " " << c[2]);
			//}
		}

		std::vector<Tag> datatags;
		std::vector<HandleType> newnodes;
		std::vector<HandleType> newpolyh;
		std::vector<HandleType> newcells;
		

		std::fstream f(File.c_str(), std::ios::in);
		XMLReader r(File, f);
		XMLReader::XMLTree t = r.ReadXML();

		if (t.GetName() == "VTKFile")
		{
			assert(t.GetAttrib("type") == "UnstructuredGrid");
			if (t.FindAttrib("compression") != t.NumAttrib())
			{
				std::cout << "Compression is specified in " << File << " but not supported" << std::endl;
				throw BadFile;
			}
			const XMLReader::XMLTree * da, * pd;
			const XMLReader::XMLTree * v = t.GetChild("UnstructuredGrid")->GetChild("Piece");
			int nnodes, ncells, ncoords;
			nnodes = atoi(v->GetAttrib("NumberOfPoints").c_str());
			ncells = atoi(v->GetAttrib("NumberOfCells").c_str());
			//first read in all the nodes
			{
				da = v->GetChild("Points")->GetChild("DataArray");
				ncoords = atoi(da->GetAttrib("NumberOfComponents").c_str());
				std::stringstream readcoords(da->GetContents());
				Storage::real xyz[3] = { 0.0, 0.0, 0.0 };
				newnodes.reserve(nnodes);
				for (int q = 0; q < nnodes; ++q)
				{
					for (int l = 0; l < ncoords; ++l)
						readcoords >> xyz[l];
					int find = -1;
					if (!old_nodes.empty())
					{
						std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(), old_nodes.end(), xyz, CentroidComparator(this));
						if (it != old_nodes.end())
						{
							Storage::real_array c = RealArrayDF(*it, CoordsTag());
							if (CentroidComparator(this).Compare(c.data(), xyz) == 0)
								find = static_cast<int>(it - old_nodes.begin());
						}
					}
					if (find == -1)
					{
						newnodes.push_back(CreateNode(xyz)->GetHandle());
						SetMarker(newnodes.back(), unused_marker);
					}
					else
						newnodes.push_back(old_nodes[find]);
				}
			}
			//read all the faces
			da = v->GetChild("Cells")->GetChildWithAttrib("Name", "faces");
			if( da )
			{
				std::stringstream faces(v->GetChild("Cells")->GetChildWithAttrib("Name", "faces")->GetContents());
				std::stringstream faceoffsets(v->GetChild("Cells")->GetChildWithAttrib("Name", "faceoffsets")->GetContents());
				int cconn, coffset = 0, totread = 0, nread, nfaces, nfacenodes;
				ElementArray<Node> hnodes(this);
				ElementArray<Face> hfaces(this);
				while (!faceoffsets.eof())
				{
					faceoffsets >> coffset;
					nread = coffset - totread;
					faces >> nfaces;
					hfaces.resize(nfaces);
					totread++;
					for (int q = 0; q < nfaces; ++q)
					{
						faces >> nfacenodes;
						hnodes.resize(nfacenodes);
						totread++;
						for (int l = 0; l < nfacenodes; ++l)
						{
							faces >> cconn;
							hnodes.at(l) = newnodes[cconn];
							RemMarker(newnodes[cconn], unused_marker);
							totread++;
						}
						hfaces[q] = CreateFace(hnodes).first;
					}
					newpolyh.push_back(CreateCell(hfaces).first.GetHandle());
				}
			}
			//check grid type
			if( grid_is_2d == 2 && ncells) //detect grid type
			{
				std::stringstream type(v->GetChild("Cells")->GetChildWithAttrib("Name", "types")->GetContents());
				int ctype;
				//bool have_2d = false;
				for (int q = 0; q < ncells && grid_is_2d == 2; ++q)
				{
					type >> ctype;
					if( ctype > 9 ) grid_is_2d = 0;
				}
				if( grid_is_2d == 2 ) grid_is_2d = 1;
			}
			
			

			if (verbosity > 0)
			{
				switch (grid_is_2d)
				{
				case 0: std::cout << "Grid has three dimensions" << std::endl; break;
				case 1: std::cout << "Grid has two dimensions" << std::endl; break;
				case 2: std::cout << "Grid has undetermined dimension" << std::endl; break;
				}
			}
			
			if( grid_is_2d == 1 && old_nodes.empty() ) SetDimensions(2);

			bool have_faces = false; //some elements go as faces
			bool have_edges = false;
			bool have_nodes = false;
			//read all the cells
			{
				std::stringstream conn(v->GetChild("Cells")->GetChildWithAttrib("Name", "connectivity")->GetContents());
				std::stringstream type(v->GetChild("Cells")->GetChildWithAttrib("Name", "types")->GetContents());
				std::stringstream offset(v->GetChild("Cells")->GetChildWithAttrib("Name", "offsets")->GetContents());
				int ctype, coffset = 0, totread = 0, nread, cconn, npolyh = 0;
				ElementArray<Face> hfaces(this);
				ElementArray<Node> hnodes(this);
				newcells.resize(ncells);
				for (int q = 0; q < ncells; ++q)
				{
					type >> ctype;
					offset >> coffset;
					nread = coffset - totread;
					hnodes.resize(nread);
					for (int l = 0; l < nread; ++l)
					{
						conn >> cconn;
						hnodes.at(l) = newnodes[cconn];
						RemMarker(newnodes[cconn], unused_marker);
						totread++;
					}
					if (ctype == 1) //VTK_VERTEX
					{
						newcells[q] = hnodes.at(0);
						have_nodes = true;
					}
					else if (ctype == 2) //VTK_POLY_VERTEX
					{
						std::cout << __FILE__ << ":" << __LINE__ << " skipping VTK_POLY_VERTEX" << std::endl;
					}
					else if (ctype == 3) //VTK_LINE
					{
						if( grid_is_2d == 1 )
						{
							ElementArray<Edge> f_edges(this,hnodes.size());
							ElementArray<Node> e_nodes(this,1);
							for(int k = 0; k < (int)hnodes.size(); ++k)
							{
								e_nodes[0] = hnodes[k];
								f_edges[k] = CreateEdge(e_nodes).first;
							}
							newcells[q] = CreateFace(f_edges).first.GetHandle();
							have_faces = true;
						}
						else
						{
							newcells[q] = CreateEdge(hnodes).first.GetHandle();
							have_edges = true;
						}
					}
					else if (ctype == 4)
					{
						std::cout << __FILE__ << ":" << __LINE__ << " skipping VTK_POLY_LINE" << std::endl;
					}
					else if (ctype == 5) //VTK_TRIANGLE
					{
						if (grid_is_2d == 1)
						{
							ElementArray<Node> e_nodes(this, 1);
							ElementArray<Edge> f_edges(this, 2);
							ElementArray<Face> c_faces(this);
							for (unsigned int k = 0; k < 3; k++)
							{
								e_nodes.at(0) = hnodes.at(k);
								f_edges.at(0) = CreateEdge(e_nodes).first->GetHandle();
								e_nodes.at(0) = hnodes.at((k + 1) % 3);
								f_edges.at(1) = CreateEdge(e_nodes).first->GetHandle();
								c_faces.push_back(CreateFace(f_edges).first);
							}
							Cell c = CreateCell(c_faces,hnodes).first;
							newcells[q] = c->GetHandle();
						}
						else
						{
							newcells[q] = CreateFace(hnodes).first->GetHandle();
							have_faces = true;
						}
						break;
					}
					else if (ctype == 6) //VTK_TRIANGLE_STRIP
					{
						std::cout << __FILE__ << ":" << __LINE__ << " skipping VTK_TRIANGLE_STRIP" << std::endl;
					}
					else if (ctype == 7) //VTK_POLYGON
					{
						if (grid_is_2d == 1)
						{
							ElementArray<Node> e_nodes(this, 1);
							ElementArray<Edge> f_edges(this, 2);
							ElementArray<Face> c_faces(this);
							for (ElementArray<Node>::size_type k = 0; k < hnodes.size(); k++)
							{
								e_nodes.at(0) = hnodes.at(k);
								f_edges.at(0) = CreateEdge(e_nodes).first->GetHandle();
								e_nodes.at(0) = hnodes.at((k + 1) % hnodes.size());
								f_edges.at(1) = CreateEdge(e_nodes).first->GetHandle();
								c_faces.push_back(CreateFace(f_edges).first);
							}
							Cell c = CreateCell(c_faces, hnodes).first;
							newcells[q] = c->GetHandle();
						}
						else
						{
							newcells[q] = CreateFace(hnodes).first->GetHandle();
							have_faces = true;
						}
						break;
					}
					else if (ctype == 8) //VTK_PIXEL
					{
						HandleType temp = hnodes.at(2);
						hnodes.at(2) = hnodes.at(3);
						hnodes.at(3) = temp;
						if (grid_is_2d == 1)
						{
							ElementArray<Node> e_nodes(this, 1);
							ElementArray<Edge> f_edges(this, 2);
							ElementArray<Face> c_faces(this);
							e_nodes.resize(1);
							f_edges.resize(2);
							for (int k = 0; k < 4; k++)
							{
								e_nodes.at(0) = hnodes.at(k);
								f_edges.at(0) = CreateEdge(e_nodes).first->GetHandle();
								e_nodes.at(0) = hnodes.at((k + 1) % 4);
								f_edges.at(1) = CreateEdge(e_nodes).first->GetHandle();
								c_faces.push_back(CreateFace(f_edges).first);
							}
							Cell c = CreateCell(c_faces, hnodes).first;
							newcells[q] = c->GetHandle();
						}
						else
						{
							newcells[q] = CreateFace(hnodes).first->GetHandle();
							have_faces = true;
						}
					}
					else if (ctype == 9) //VTK_QUAD
					{
						if (grid_is_2d == 1)
						{
							ElementArray<Node> e_nodes(this,1);
							ElementArray<Edge> f_edges(this,2);
							ElementArray<Face> c_faces(this);
							for (int k = 0; k < 4; k++)
							{
								e_nodes.at(0) = hnodes.at(k);
								f_edges.at(0) = CreateEdge(e_nodes).first->GetHandle();
								e_nodes.at(0) = hnodes.at((k + 1) % 4);
								f_edges.at(1) = CreateEdge(e_nodes).first->GetHandle();
								c_faces.push_back(CreateFace(f_edges).first);
							}
							Cell c = CreateCell(c_faces, hnodes).first;
							newcells[q] = c->GetHandle();
						}
						else
						{
							newcells[q] = CreateFace(hnodes).first->GetHandle();
							have_faces = true;
						}
					}
					else if (ctype == 10) //VTK_TETRA
					{
						assert(nread == 4);
						const integer nodesnum[12] = { 0, 2, 1, 0, 1, 3, 1, 2, 3, 0, 3, 2 };
						const integer sizes[4] = { 3, 3, 3, 3 };
						newcells[q] = CreateCell(hnodes, nodesnum, sizes, 4).first.GetHandle();
					}
					else if (ctype == 12 || ctype == 11) // VTK_HEXAHEDRON or VTK_VOXEL
					{
						assert(nread == 8);
						if( ctype == 11 )
						{
							HandleType temp;
							temp = hnodes.at(2);
							hnodes.at(2) = hnodes.at(3);
							hnodes.at(3) = temp;
							temp = hnodes.at(6);
							hnodes.at(6) = hnodes.at(7);
							hnodes.at(7) = temp;
						}
						const integer nodesnum[24] = { 0, 4, 7, 3, 1, 2, 6, 5, 0, 1, 5, 4, 3, 7, 6, 2, 0, 3, 2, 1, 4, 5, 6, 7 };
						const integer sizes[6] = { 4, 4, 4, 4, 4, 4 };
						newcells[q] = CreateCell(hnodes, nodesnum, sizes, 6).first.GetHandle();
					}
					else if (ctype == 13) // VTK_WEDGE
					{
						assert(nread == 6);
						const integer nodesnum[18] = { 0, 2, 5, 3, 1, 4, 5, 2, 0, 3, 4, 1, 3, 5, 4, 0, 1, 2 };
						//const integer nodesnum[18] = { 0, 3, 5, 2, 0, 1, 4, 3, 1, 4, 5, 2, 3, 4, 5, 0, 2, 1 };
						const integer sizes[5] = { 4, 4, 4, 3, 3 };
						newcells[q] = CreateCell(hnodes, nodesnum, sizes, 5).first.GetHandle();
					}
					else if (ctype == 14) //VTK_PYRAMID
					{
						assert(nread == 5);
						const integer nodesnum[16] = { 0, 4, 3, 0, 1, 4, 1, 2, 4, 3, 4, 2, 0, 3, 2, 1 };
						const integer sizes[5] = { 3, 3, 3, 3, 4 };
						newcells[q] = CreateCell(hnodes, nodesnum, sizes, 5).first.GetHandle();
					}
					else if (ctype == 15) //VTK_PENTAGONAL_PRISM
					{
						assert(nread == 10);
						const integer nodesnum[30] =
						{
							0, 1, 6, 5,
							1, 2, 7, 6,
							2, 3, 8, 7,
							3, 4, 9, 8,
							4, 0, 5, 9,
							5, 6, 7, 8, 9,
							4, 3, 2, 1, 0
						};
						//5, 6, 7, 8, 9
						//0, 1, 2, 3, 4

						//9, 8, 7, 6, 5
						//0, 1, 2, 3, 4
						const integer sizes[7] = { 4, 4, 4, 4, 4, 5, 5 };
						newcells[q] = CreateCell(hnodes, nodesnum, sizes, 7).first.GetHandle();
					}
					else if (ctype == 16) //VTK_HEXAGONAL_PRISM
					{
						assert(nread == 12);
						const integer nodesnum[36] =
						{
							0, 1, 7, 6,
							1, 2, 8, 7,
							2, 3, 9, 8,
							3, 4, 10, 9,
							4, 5, 11, 10,
							5, 0, 6, 11,
							6, 7, 8, 9, 10, 11,
							5, 4, 3, 2, 1, 0
						};
						//6, 7, 8, 9, 10, 11
						//0, 1, 2, 3, 4, 5
						const integer sizes[8] = { 4, 4, 4, 4, 4, 4, 6, 6 };
						newcells[q] = CreateCell(hnodes, nodesnum, sizes, 8).first.GetHandle();
					}
					else if (ctype == 42) //VTK_POLYHEDRON
						newcells[q] = newpolyh[npolyh++];
					else std::cout << __FILE__ << ":" << __LINE__ << " Implement VTK format " << ctype << std::endl;
				}
			}
			//read data
			{
				std::string dname[2] = { "PointData", "CellData" };
				ElementType dtype[2] = { NODE, CELL };
				ElementType dsparse[2] = { NONE, NONE };
				if (have_faces)
				{
					dtype[1] |= FACE;
					dsparse[1] |= FACE;
				}
				if (have_edges)
				{
					dtype[1] |= EDGE;
					dsparse[1] |= EDGE;
				}
				if (have_nodes)
				{
					dtype[1] |= NODE;
					dsparse[1] |= NODE;
				}
				HandleType * darray[2] = { &newnodes[0], &newcells[0] };
				int dsize[2] = { nnodes, ncells };
				for (int j = 0; j < 2; ++j)
				{
					da = v->GetChild(dname[j]);
					if (da)
					{
						for (int k = 0; k < da->NumChildren(); ++k)
						{
							pd = &da->GetChild(k);
							if (pd->GetName() == "DataArray")
							{
								int ncomps = 1;
								int nca = pd->FindAttrib("NumberOfComponents");
								if (nca != pd->NumAttrib()) ncomps = atoi(pd->GetAttrib(nca).value.c_str());
								TagRealArray t = CreateTag(pd->GetAttrib("Name"), DATA_REAL, dtype[j], dsparse[j], ncomps);
								std::stringstream inp(pd->GetContents());
								for (int l = 0; l < dsize[j]; ++l)
								{
									for (INMOST_DATA_ENUM_TYPE q = 0; q < t.GetSize(); ++q)
										inp >> t[darray[j][l]][q];
								}
							}
							else std::cout << __FILE__ << ":" << __LINE__ << "I don't know yet what is " << pd->GetName() << " in point data" << std::endl;
						}
					}
				}
			}
		}
		else
		{
			std::cout << "Root tag is not VTKFile, " << File << " is not a valid .vtu file" << std::endl;
			throw BadFile;
		}

		ReleaseMarker(unused_marker,FACE|NODE);
		f.close();
	}
}

#endif
