#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "inmost.h"
//gmsh states

#if defined(USE_MESH)

#define GMSH_VERNONE 0
#define GMSH_VER1 1
#define GMSH_VER2 2

#define GMSH_NONE 0
#define GMSH_NODES_TOT 1
#define GMSH_NODES_NUM 2
#define GMSH_NODES_X 3
#define GMSH_NODES_Y 4
#define GMSH_NODES_Z 5
#define GMSH_ELEMENTS_TOT 6
#define GMSH_ELEMENTS_NUM 7
#define GMSH_ELEMENTS_TYPE 8
#define GMSH_ELEMENTS_REGPHYS 9
#define GMSH_ELEMENTS_REGELEM 10
#define GMSH_ELEMENTS_NUMTAGS 11
#define GMSH_ELEMENTS_TAGS 12
#define GMSH_ELEMENTS_NUMNODES 13
#define GMSH_ELEMENTS_NODELIST 14
#define GMSH_FORMAT 15
#define GMSH_SKIP_KEYWORD 16
#define GMSH_PHYSICAL_NAMES 17
#define GMSH_PHYSICAL_DIM 18
#define GMSH_PHYSICAL_TAG 19
#define GMSH_PHYSICAL_NAM 20

namespace INMOST
{
  void Mesh::LoadMSH(std::string File)
  {
    int verbosity = 0;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < file_options.size(); ++k)
		{
			if( file_options[k].first == "VERBOSITY" )
			{
				verbosity = atoi(file_options[k].second.c_str());
				if( verbosity < 0 || verbosity > 2 )
				{
					printf("%s:%d Unknown verbosity option: %s\n",__FILE__,__LINE__,file_options[k].second.c_str());
					verbosity = 1;
				}
			}
		}

		Tag gmsh_tags;
		char readline[2048], sname[2048], * p, * pend;
		std::vector<HandleType> newnodes;
		std::vector<HandleType> newelems;
		std::map< std::pair<int,int>, ElementSet> newsets; //tag,dim : set
		TagInteger set_tag, set_dim;
		double xyz[3];
		int text_start, text_end, state = GMSH_NONE, ver = GMSH_VERNONE;
		int nnodes, ncells, nchars, nnames, sdim, stag, nodenum, elemnum, elemtype, numtags, elemnodes, temp;
		int ascii, float_size, verlow, verhigh;
		char skip_keyword[2048];
		std::vector<int> elemtags, ntags, tags, elemtypes;
		std::vector<int> nodelist, nnods, nods;
		ElementArray<Node> c_nodes(this);
		FILE * f = fopen(File.c_str(),"r");
		if( f == NULL )
		{
			std::cout << __FILE__ << ":" << __LINE__ << " cannot open " << File << std::endl;
			throw BadFileName;
		}
		std::vector<HandleType> old_nodes(NumberOfNodes());
		{
			unsigned qq = 0;
			for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it)
				old_nodes[qq++] = *it;
		}
		if( !old_nodes.empty() ) 
			std::sort(old_nodes.begin(),old_nodes.end(),CentroidComparator(this));
		int line = 0;
		unsigned report_pace;
		while( fgets(readline,2048,f) != NULL )
		{
			line++;
			if( readline[strlen(readline)-1] == '\n' ) readline[strlen(readline)-1] = '\0';
			text_end = static_cast<int>(strlen(readline));
			for(text_start = 0; isspace(readline[text_start]) && text_start < text_end; text_start++);
			if( text_start == text_end ) continue;
			for(text_end = text_end-1; isspace(readline[text_end]) && text_end > text_start; text_end--);
			readline[text_end+1] = '\0';
			p = readline + text_start;
			pend = readline + text_end + 1;
			switch(state)
			{
			case GMSH_NONE:
				if(!strcmp(p,"$NOD"))
				{
					ver = GMSH_VER1;
					state = GMSH_NODES_TOT;
					elemtags.resize(2);
					gmsh_tags = CreateTag("GMSH_TAGS",DATA_INTEGER,CELL|FACE|EDGE|NODE,FACE|EDGE|NODE,2);
				}
				else if( !strcmp(p,"$ELM") )
				{
					if( ver == GMSH_NONE )
					{
						std::cout << __FILE__ << ":" << __LINE__ <<  " bad keyword sequence in " << File << " line " << line << std::endl;
						throw BadFile;
					}
					state = GMSH_ELEMENTS_TOT;
				}
				else if( !strcmp(p,"$MeshFormat") )
				{
					ver = GMSH_VER2;
					state = GMSH_FORMAT;
					gmsh_tags = CreateTag("GMSH_TAGS",DATA_INTEGER,CELL|FACE|EDGE|NODE,FACE|EDGE|NODE);
				}
				else if( !strcmp(p,"$Nodes") )
				{
					ver = GMSH_VER2;
					state = GMSH_NODES_TOT;
					gmsh_tags = CreateTag("GMSH_TAGS",DATA_INTEGER,CELL|FACE|EDGE|NODE,FACE|EDGE|NODE);
				}
				else if( !strcmp(p,"$Elements") )
				{
					if( ver == GMSH_NONE )
					{
						std::cout << __FILE__ << ":" << __LINE__ <<  " bad keyword sequence in " << File << " line " << line <<std::endl;
						throw BadFile;
					}
					state = GMSH_ELEMENTS_TOT;
				}
				else if (!strcmp(p, "$PhysicalNames"))
				{
					if (ver == GMSH_NONE)
					{
						std::cout << __FILE__ << ":" << __LINE__ << " bad keyword sequence in " << File << " line " << line << std::endl;
						throw BadFile;
					}
					state = GMSH_PHYSICAL_NAMES;
				}
				else
				{
					if( p[0] == '$' )
					{
						skip_keyword[0] = '$';
						skip_keyword[1] = 'E';
						skip_keyword[2] = 'n';
						skip_keyword[3] = 'd';
						for(int q = 1; !(p[q] == '\0' || isspace(p[q])); q++)
						{
							skip_keyword[q+3] = p[q];
							skip_keyword[q+4] = '\0';
						}
						state = GMSH_SKIP_KEYWORD;
						std::cout << __FILE__ << ":" << __LINE__ << " skipping unknown keyword " << skip_keyword+4 << " in " << File << " line " << line << std::endl;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ <<  " Unexpected line for " << File << " format: " << p << " line " << line <<std::endl;
						throw BadFile;
					}
				}
				break;
			case GMSH_SKIP_KEYWORD:
				if( !strncmp(p,skip_keyword,sizeof(skip_keyword) ) )
				{
					state = GMSH_NONE;
				}
				break;
			case GMSH_PHYSICAL_NAMES:
				if (1 == sscanf(p, "%d%n", &nnames, &nchars))
				{
					set_tag = CreateTag("GMSH_PHYSICAL_TAG", DATA_INTEGER, ESET, ESET, 1);
					set_dim = CreateTag("GMSH_PHYSICAL_DIM", DATA_INTEGER, ESET, ESET, 1);
					p += nchars;
					while (isspace(*p) && p < pend) ++p;
					state = GMSH_PHYSICAL_DIM;
					if (verbosity > 0) printf("Reading %d set names.\n", nnames);
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " error reading number of sets in " << File << std::endl;
					throw BadFile;
				}
				if (p >= pend) break;
				//FALLTHROUGH
			case GMSH_PHYSICAL_DIM:
				if (1 == sscanf(p, "%d%n", &sdim, &nchars))
				{
					p += nchars;
					while (isspace(*p) && p < pend) ++p;
					state = GMSH_PHYSICAL_TAG;
				}
				else
				{
					if (ver == GMSH_VER2 && !strcmp(p, "$EndPhysicalNames"))
					{
						if (nnames) std::cout << __FILE__ << ":" << __LINE__ << " number of element records mismatch with reported number of physical names " << nnames << " in " << File << " line " << line << std::endl;
						state = GMSH_NONE;
						break;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " error reading set dimension in " << File << " line " << line << std::endl;
						throw BadFile;
					}
				}
				if (p >= pend) break;
				//FALLTHROUGH
			case GMSH_PHYSICAL_TAG:
				if (1 == sscanf(p, "%d%n", &stag, &nchars))
				{
					p += nchars;
					while (isspace(*p) && p < pend) ++p;
					state = GMSH_PHYSICAL_NAM;
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " error reading set tag in " << File << " line " << line << std::endl;
					throw BadFile;
				}
				if (p >= pend) break;
				//FALLTHROUGH
			case GMSH_PHYSICAL_NAM:
				if (1 == sscanf(p, "%s%n", sname, &nchars))
				{
					nnames--; //count of remaining sets
					std::string strname = sname;
					if (strname[0] == '\"' && strname[strname.size() - 1] == '\"')
						strname = strname.substr(1, strname.size() - 2);
					ElementSet set = CreateSet(strname).first;
					set_tag[set] = stag;
					set_dim[set] = sdim;
					newsets[std::make_pair(stag, sdim)] = set;
					p += nchars;
					while (isspace(*p) && p < pend) ++p;
					state = GMSH_PHYSICAL_DIM;
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " error reading set dimension in " << File << " line " << line << std::endl;
					throw BadFile;
				}
				if (p >= pend) break;
				//FALLTHROUGH
			case GMSH_NODES_TOT:
				if( 1 == sscanf(p,"%d%n",&nnodes,&nchars) )
				{
					p += nchars;
					while(isspace(*p) && p < pend) ++p;
					newnodes.resize(nnodes);
					state = GMSH_NODES_NUM;
					report_pace = std::max<unsigned>(nnodes/250,1);
					if( verbosity > 0 ) printf("Reading %d nodes.\n",nnodes);
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " error reading number of nodes in " << File << std::endl;
					throw BadFile;
				}
				if( p >= pend ) break;
				//FALLTHROUGH
			case GMSH_NODES_NUM:
read_node_num_link:
				if( 1 == sscanf(p,"%d%n",&nodenum,&nchars) )
				{
					--nodenum;
					if( nodenum >= static_cast<int>(newnodes.size()) ) std::cout << __FILE__ << ":" << __LINE__ << " node number is bigger then total number of nodes " << nodenum << " / " << newnodes.size() << " line " << line <<std::endl;
					p += nchars;
					while(isspace(*p) && p < pend) ++p;
					state = GMSH_NODES_X;
				}
				else
				{
					if( ver == GMSH_VER1 && !strcmp(p,"$ENDNOD"))
					{
						state = GMSH_NONE;
						if( nnodes ) std::cout << __FILE__ << ":" << __LINE__ << " number of node records mismatch with reported number of nodes " << nnodes << " in " << File << " line " << line <<std::endl;
						break;
					}
					else if( ver == GMSH_VER2 && !strcmp(p,"$EndNodes") )
					{
						state = GMSH_NONE;
						if( nnodes ) std::cout << __FILE__ << ":" << __LINE__ << " number of node records mismatch with reported number of nodes " << nnodes << " in " << File << " line " << line <<std::endl;
						break;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " error reading node number in " << File << " got " << p << " instead " << " line " << line <<std::endl;
						throw BadFile;
					}
				}
				if( p >= pend ) break;
				//FALLTHROUGH
			case GMSH_NODES_X:
				if( 1 == sscanf(p,"%lf%n",xyz,&nchars) )
				{
					p += nchars;
					while(isspace(*p) && p < pend) ++p;
					state = GMSH_NODES_Y;
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " error reading X-coord in " << File << " line " << line <<std::endl;
					throw BadFile;
				}
				if( p >= pend ) break;
				//FALLTHROUGH
			case GMSH_NODES_Y:
				if( 1 == sscanf(p,"%lf%n",xyz+1,&nchars) )
				{
					p += nchars;
					while(isspace(*p) && p < pend) ++p;
					state = GMSH_NODES_Z;
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " error reading Y-coord in " << File << " line " << line <<std::endl;
					throw BadFile;
				}
				if( p >= pend ) break;
				//FALLTHROUGH
			case GMSH_NODES_Z:
				if( 1 == sscanf(p,"%lf%n",xyz+2,&nchars) )
				{
					p += nchars;
					while(isspace(*p) && p < pend) ++p;
					int find = -1;
					INMOST_DATA_REAL_TYPE rxyz[3];
					rxyz[0] = xyz[0], rxyz[1] = xyz[1], rxyz[2] = xyz[2];		
					if( !old_nodes.empty() )
					{
						std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),rxyz,CentroidComparator(this));
						if( it != old_nodes.end() ) 
						{
							Storage::real_array c = RealArrayDF(*it,CoordsTag());
							if( CentroidComparator(this).Compare(rxyz,c.data()) == 0 )
								find = static_cast<int>(it - old_nodes.begin());
						}
					}
					if( find == -1 ) newnodes[nodenum] = CreateNode(rxyz)->GetHandle();
					else newnodes[nodenum] = old_nodes[find];
					nnodes--;
					if( verbosity > 1 && (newnodes.size()-nnodes)%report_pace == 0 )
					{
						printf("%3.1f%%\r",((newnodes.size()-nnodes)*100.0)/(newnodes.size()*1.0));
						fflush(stdout);
					}
					state = GMSH_NODES_NUM;
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " error reading Z-coord in " << File << " line " << line <<std::endl;
					throw BadFile;
				}
				if( p < pend ) goto read_node_num_link; else break;
			case GMSH_ELEMENTS_TOT:
				if( 1 == sscanf(p,"%d%n",&ncells,&nchars) )
				{
					p += nchars;
					while(isspace(*p) && p < pend) ++p;
					newelems.resize(ncells);
					state = GMSH_ELEMENTS_NUM;
					report_pace = std::max<unsigned>(ncells/250,1);
					if(verbosity > 0) printf("Reading %d elements\n",ncells);
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " error reading total number of elements in " << File << " line " << line <<std::endl;
					throw BadFile;
				}
				if( p >= pend ) break;
				//FALLTHROUGH
			case GMSH_ELEMENTS_NUM:
read_elem_num_link:
				if( 1 == sscanf(p,"%d%n",&elemnum,&nchars) )
				{
					--elemnum;
					//std::cout << __FILE__ << ":" << __LINE__ << " element number: " << elemnum << " line " << line << std::endl;
					if( elemnum >= static_cast<int>(newelems.size()) ) std::cout << __FILE__ << ":" << __LINE__ << " element number is bigger then total number of elements " << elemnum << " / " << newelems.size() << " line " << line <<std::endl;
					p += nchars;
					while(isspace(*p) && p < pend) ++p;
					state = GMSH_ELEMENTS_TYPE;
				}
				else
				{
					if( ver == GMSH_VER1 && !strcmp(p,"$ENDELM"))
					{
						if( ncells ) std::cout << __FILE__ << ":" << __LINE__ << " number of element records mismatch with reported number of elements " << ncells << " in " << File << " line " << line <<std::endl;
						state = GMSH_NONE;
					}
					else if( ver == GMSH_VER2 && !strcmp(p,"$EndElements") )
					{
						if( ncells ) std::cout << __FILE__ << ":" << __LINE__ << " number of element records mismatch with reported number of elements " << ncells << " in " << File << " line " << line <<std::endl;
						state = GMSH_NONE;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " error reading element number in " << File << " got " << p << " instead " << " line " << line <<std::endl;
						throw BadFile;
					}
					if (state == GMSH_NONE) //End of data, process
					{
						bool is2d = true;
						for (size_t q = 0; q < elemtypes.size() && is2d; ++q)
						{
							elemtype = elemtypes[q];
							switch (elemtype)
							{
							case 1: //edge
							case 2: //triangle
							case 3: //quad
								break;
							case 4: //tetra
							case 5: //hex
							case 6: //prism
							case 7: //pyramid
								is2d = false; 
								break;
							case 8: //second order edge
							case 9: //second order tri with nodes on edges
							case 10: //second order quad with nodes on edges and in center
								break;
							case 11: // second order tet with nodes on edges
							case 12: // second order hex with nodes in centers of edges and faces and in center of volume
							case 13: // second order prism  with nodes in centers of edges and quad faces
							case 14: // second order pyramid with nodes in centers of edges and quad faces
								is2d = false;
								break;
							case 15: // vertex
							case 16: // second order quad with nodes on edges
								break;
							case 17: // second order hex with nodes on edges
							case 18: // second order prism with nodes on edges
							case 19: // second order pyramid with nodes on edges
								is2d = false;
								break;
							}
						}
						size_t qnod = 0, qtag = 0;
						for (size_t q = 0; q < elemtypes.size(); ++q)
						{
							elemnum = q;
							elemtype = elemtypes[q];
							//fill arrays
							for (size_t k = 0; k < nnods[q]; ++k)
								nodelist.push_back(nods[qnod + k]);
							qnod += nnods[q];
							for (size_t k = 0; k < ntags[q]; ++k)
								elemtags.push_back(tags[qtag + k]);
							qtag += ntags[q];

							for (size_t k = 0; k < nodelist.size(); ++k)
								c_nodes.push_back(newnodes[nodelist[k]]);

							//Create new element
							switch (elemtype)
							{
							case 1: //edge
								if( is2d )
									newelems[elemnum] = CreateFace(c_nodes).first->GetHandle();
								else newelems[elemnum] = CreateEdge(c_nodes).first->GetHandle(); 
								break; 
							case 2: //triangle
								if (is2d)
								{
									const integer nodesnum[6] = { 0,1,1,2,2,0 };
									const integer sizes[3] = { 2,2,2 };
									newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 3).first->GetHandle();
								}
								else newelems[elemnum] = CreateFace(c_nodes).first->GetHandle(); 
								break; 
							case 3: //quad
								if (is2d)
								{
									const integer nodesnum[8] = { 0,1,1,2,2,3,3,0 };
									const integer sizes[4] = { 2,2,2,2 };
									newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 4).first->GetHandle();
								}
								else newelems[elemnum] = CreateFace(c_nodes).first->GetHandle();
								break; 
							case 4: //tetra
							{
								const integer nodesnum[12] = { 0,2,1,0,1,3,1,2,3,0,3,2 };
								const integer sizes[4] = { 3,3,3,3 };
								newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 4).first->GetHandle();
							}
							break;
							case 5: //hex
							{
								//const integer nodesnum[24] = {0,4,7,3,1,2,6,5,0,1,5,4,3,7,6,2,0,3,2,1,4,5,6,7};
								const integer nodesnum[24] = { 0,3,2,1,4,5,6,7,0,4,7,3,1,2,6,5,0,1,5,4,2,3,7,6 };
								const integer sizes[6] = { 4,4,4,4,4,4 };
								newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 6).first->GetHandle();
							}
							break;
							case 6: //prism
							{
								//const integer nodesnum[18] = {0,2,5,3,1,4,5,2,0,3,4,1,3,5,4,0,1,2};
								const integer nodesnum[18] = { 0,3,5,2,1,2,5,4,0,1,4,3,3,4,5,0,2,1 };
								const integer sizes[5] = { 4,4,4,3,3 };
								newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 5).first->GetHandle();
							}
							break;
							case 7: //pyramid
							{
								const integer nodesnum[16] = { 0,4,3,0,1,4,1,2,4,3,4,2,0,3,2,1 };
								const integer sizes[5] = { 3,3,3,3,4 };
								newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 5).first->GetHandle();
							}
							break;
							case 8: //second order edge
							{
								//treat as just an edge
								if (is2d)
									newelems[elemnum] = CreateFace(ElementArray<Node>(this, c_nodes.data(), c_nodes.data() + 2)).first->GetHandle();
								else newelems[elemnum] = CreateEdge(ElementArray<Node>(this, c_nodes.data(), c_nodes.data() + 2)).first->GetHandle();
							}
							break;
							case 9: //second order tri with nodes on edges
							{
								if (is2d)
								{
									const integer nodesnum[6] = { 0,1,1,2,2,0 };
									const integer sizes[3] = { 2,2,2 };
									newelems[elemnum] = CreateCell(ElementArray<Node>(this, c_nodes.data(), c_nodes.data() + 3), nodesnum, sizes, 3).first->GetHandle();
								}
								else newelems[elemnum] = CreateFace(ElementArray<Node>(this, c_nodes.data(), c_nodes.data() + 3)).first->GetHandle();
							}
							break;
							case 10: //second order quad with nodes on edges and in center
							{
								if (is2d)
								{
									const integer nodesnum[8] = { 0,1,1,2,2,3,3,0 };
									const integer sizes[4] = { 2,2,2,2 };
									newelems[elemnum] = CreateCell(ElementArray<Node>(this, c_nodes.data(), c_nodes.data() + 4), nodesnum, sizes, 4).first->GetHandle();
								}
								else newelems[elemnum] = CreateFace(ElementArray<Node>(this, c_nodes.data(), c_nodes.data() + 4)).first->GetHandle();
							}
							break;
							case 11: // second order tet with nodes on edges
							{
								const integer nodesnum[12] = { 0,2,1,0,1,3,1,2,3,0,3,2 };
								const integer sizes[4] = { 3,3,3,3 };
								newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 4).first->GetHandle();
							}
							break;
							case 12: // second order hex with nodes in centers of edges and faces and in center of volume
							{
								const integer nodesnum[24] = { 0,4,7,3,1,2,6,5,0,1,5,4,3,7,6,2,0,3,2,1,4,5,6,7 };
								const integer sizes[6] = { 4,4,4,4,4,4 };
								newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 6).first->GetHandle();
							}
							break;
							case 13: // second order prism  with nodes in centers of edges and quad faces
							{
								const integer nodesnum[18] = { 0,3,5,2,1,2,5,4,0,1,4,3,3,4,5,0,2,1 };
								const integer sizes[5] = { 4,4,4,3,3 };
								newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 5).first->GetHandle();
							}
							break;
							case 14: // second order pyramid with nodes in centers of edges and quad faces
							{
								const integer nodesnum[16] = { 0,4,3,0,1,4,1,2,4,3,4,2,0,3,2,1 };
								const integer sizes[5] = { 3,3,3,3,4 };
								newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 5).first->GetHandle();
							}
							break;
							case 15: // vertex
							{
								newelems[elemnum] = c_nodes.at(0);
							}
							break;
							case 16: // second order quad with nodes on edges
							{
								if (is2d)
								{
									const integer nodesnum[8] = { 0,1,1,2,2,3,3,0 };
									const integer sizes[4] = { 2,2,2,2 };
									newelems[elemnum] = CreateCell(ElementArray<Node>(this, c_nodes.data(), c_nodes.data() + 4), nodesnum, sizes, 4).first->GetHandle();
								}
								else newelems[elemnum] = CreateFace(ElementArray<Node>(this, c_nodes.data(), c_nodes.data() + 4)).first->GetHandle();
							}
							break;
							case 17: // second order hex with nodes on edges
							{
								const integer nodesnum[24] = { 0,4,7,3,1,2,6,5,0,1,5,4,3,7,6,2,0,3,2,1,4,5,6,7 };
								const integer sizes[6] = { 4,4,4,4,4,4 };
								newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 6).first->GetHandle();
							}
							break;
							case 18: // second order prism with nodes on edges
							{
								const integer nodesnum[18] = { 0,3,5,2,1,2,5,4,0,1,4,3,3,4,5,0,2,1 };
								const integer sizes[5] = { 4,4,4,3,3 };
								newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 5).first->GetHandle();
							}
							break;
							case 19: // second order pyramid with nodes on edges
							{
								const integer nodesnum[16] = { 0,4,3,0,1,4,1,2,4,3,4,2,0,3,2,1 };
								const integer sizes[5] = { 3,3,3,3,4 };
								newelems[elemnum] = CreateCell(c_nodes, nodesnum, sizes, 5).first->GetHandle();
							}
							break;
							}


							//add to set
							if (!elemtags.empty())
							{
								stag = elemtags[0]; //assume first tag is physical tag
								sdim = Element::GetGeometricDimension(GetGeometricType(newelems[elemnum]));
								if (newsets.find(std::make_pair(stag, sdim)) != newsets.end())
									newsets.find(std::make_pair(stag, sdim))->second.AddElement(newelems[elemnum]);
							}


							Storage::integer_array ptags = IntegerArray(newelems[elemnum], gmsh_tags);

							if (ver == GMSH_VER2)
							{
								ptags.replace(ptags.begin(), ptags.end(), elemtags.begin(), elemtags.end());
								elemtags.clear();
							}
							else
							{
								ptags[0] = elemtags[0];
								ptags[1] = elemtags[1];
							}

							nodelist.clear();
							elemtags.clear();
							c_nodes.clear();
						}
						break;
					}
				}
				if( p >= pend ) break;
				//FALLTHROUGH
			case GMSH_ELEMENTS_TYPE:
				if( 1 == sscanf(p,"%d%n",&elemtype,&nchars) )
				{
					const int type_nodes[19] =
					{
						2,3,4,4,8,6,5,3,6,9,10,
						27,18,14,1,8,20,15,13
					};
					if( elemtype < 1 || elemtype > 19 )
					{
						std::cout << __FILE__ << ":" << __LINE__ << " bad element type " << elemtype << " expected in [1,19] range" << " line " << line <<std::endl;
						throw BadFile;
					}

					p += nchars;
					while(isspace(*p) && p < pend) ++p;
					if( ver == GMSH_VER1 )
						state = GMSH_ELEMENTS_REGPHYS;
					else if( ver == GMSH_VER2 )
						state = GMSH_ELEMENTS_NUMTAGS;

					elemnodes = type_nodes[elemtype-1];
					nodelist.resize(elemnodes);
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " error reading element type in " << File << " got " << p << " instead " << " line " << line <<std::endl;
					throw BadFile;
				}
				if( p >= pend ) break;
				//FALLTHROUGH
			case GMSH_ELEMENTS_REGPHYS:
				if( state == GMSH_ELEMENTS_REGPHYS )
				{
					if( 1 == sscanf(p,"%d%n",&elemtags[0],&nchars) )
					{
						//std::cout << __FILE__ << ":" << __LINE__ << " i should not be here " << std::endl;
						p += nchars;
						while(isspace(*p) && p < pend) ++p;
						state = GMSH_ELEMENTS_REGELEM;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " error reading physical tag in " << File << " got " << p << " instead " << " line " << line <<std::endl;
						throw BadFile;
					}
					if( p >= pend ) break;
				}
				//FALLTHROUGH
			case GMSH_ELEMENTS_REGELEM:
				if( state == GMSH_ELEMENTS_REGELEM )
				{
					if( 1 == sscanf(p,"%d%n",&elemtags[1],&nchars) )
					{
						//std::cout << __FILE__ << ":" << __LINE__ << " i should not be here " << std::endl;
						p += nchars;
						while(isspace(*p) && p < pend) ++p;
						state = GMSH_ELEMENTS_NUMNODES;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " error reading element tag in " << File << " got " << p << " instead " << " line " << line <<std::endl;
						throw BadFile;
					}
					if( p >= pend ) break;
				}
				//FALLTHROUGH
			case GMSH_ELEMENTS_NUMNODES:
				if( state == GMSH_ELEMENTS_NUMNODES )
				{
					if( 1 == sscanf(p,"%d%n",&temp,&nchars) )
					{
						//std::cout << __FILE__ << ":" << __LINE__ << " i should not be here " << std::endl;
						if( temp != elemnodes )
						{
							std::cout << __FILE__ << ":" << __LINE__ << " bad number of nodes " << temp << " for type " << elemtype << " expected " << elemnodes << " line " << line <<std::endl;
							throw BadFile;
						}
						p += nchars;
						while(isspace(*p) && p < pend) ++p;
						state = GMSH_ELEMENTS_NODELIST;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " error reading number of nodes of element in " << File << " got " << p << " instead " << " line " << line <<std::endl;
						throw BadFile;
					}
					if( p >= pend ) break;
				}
				//FALLTHROUGH
			case GMSH_ELEMENTS_NUMTAGS:
				if( state == GMSH_ELEMENTS_NUMTAGS )
				{
					if( 1 == sscanf(p,"%d%n",&numtags,&nchars) )
					{
						//std::cout << __FILE__ << ":" << __LINE__ << " numtags: " << numtags << " line " << line << std::endl;
						elemtags.resize(numtags);
						p += nchars;
						while(isspace(*p) && p < pend) ++p;
						state = GMSH_ELEMENTS_TAGS;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " error reading number of tags " << File << " got " << p << " instead " << " line " << line <<std::endl;
						throw BadFile;
					}
					if( p >= pend ) break;
				}
				//FALLTHROUGH
			case GMSH_ELEMENTS_TAGS:
				if( state == GMSH_ELEMENTS_TAGS )
				{
					while( numtags > 0 && p < pend )
					{
						if( 1 == sscanf(p,"%d%n",&elemtags[elemtags.size()-numtags],&nchars) )
						{
							p += nchars;
							while(isspace(*p) && p < pend) ++p;
							if( --numtags == 0 ) state = GMSH_ELEMENTS_NODELIST;
						}
						else
						{
							std::cout << __FILE__ << ":" << __LINE__ << " error reading element tag in " << File << " got " << p << " instead " << " line " << line <<std::endl;
							throw BadFile;
						}
					}
					if( p >= pend ) break;
				}
				//FALLTHROUGH
			case GMSH_ELEMENTS_NODELIST:
				while( elemnodes > 0 && p < pend )
				{
					if( 1 == sscanf(p,"%d%n",&nodelist[nodelist.size()-elemnodes],&nchars) )
					{
						--nodelist[nodelist.size()-elemnodes];
						if( nodelist[nodelist.size()-elemnodes] < 0 || nodelist[nodelist.size()-elemnodes] >= static_cast<int>(newnodes.size()) )
						{
							std::cout << __FILE__ << ":" << __LINE__ << " got node " << nodelist[nodelist.size()-elemnodes] << " but must be within [0:" << newnodes.size()-1 << "] in file " << File << " line " << line << std::endl;
							throw BadFile;
						}
						p += nchars;
						while(isspace(*p) && p < pend) ++p;

						//std::cout << __FILE__ << ":" << __LINE__ << " node[" << nodelist.size()-elemnodes << "] = " << nodelist[nodelist.size()-elemnodes] << " line " << line << std::endl;



						if( --elemnodes == 0 ) 
						{
							elemtypes.push_back(elemtype);
							nnods.push_back(nodelist.size());
							ntags.push_back(elemtags.size());
							nods.insert(nods.end(), nodelist.begin(), nodelist.end());
							tags.insert(tags.end(), elemtags.begin(), elemtags.end());
							
							nodelist.clear();
							if (ver == GMSH_VER2)
								elemtags.clear();

							state = GMSH_ELEMENTS_NUM;
							ncells--;
							if( verbosity > 1 && (elemtypes.size()-ncells)%report_pace == 0 )
							{
								printf("%3.1f%%\r",((elemtypes.size()-ncells)*100.0)/(elemtypes.size()*1.0));
								fflush(stdout);
							}
						}
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " error reading node list for element in " << File << " got " << p << " instead " << " line " << line <<std::endl;
						throw BadFile;
					}
				}
				if( p < pend ) goto read_elem_num_link; else break;
			case GMSH_FORMAT:
				if( 4 == sscanf(p, "%1d.%1d %d %d\n",&verlow,&verhigh,&ascii,&float_size) )
				{
					if( !(verlow == 2 && verhigh == 0) && verbosity > 0)
					{
						std::cout << __FILE__ << ":" << __LINE__ << " version of file " << File << " is " << verlow << "." << verhigh << " expected 2.0 ";
						std::cout << " in " << File << " line " << line << std::endl;
					}
					if( ascii != 0 )
					{
						std::cout << __FILE__ << ":" << __LINE__ << " file " << File << " is binary, gmsh binary is not currently supported " << " line " << line <<std::endl;
						throw BadFile;
					}
				}
				else
				{
					if( !strcmp(p,"$EndMeshFormat") )
					{
						state = GMSH_NONE;
						break;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " unexpected line content " << p << " in file " << File << " line " << line <<std::endl;
					}
				}
				break;
			}
		}

		if( state != GMSH_NONE )
		{
			std::cout << __FILE__ ":" << __LINE__ << " probably file " << File << " is not complitely written " << " line " << line <<std::endl;
		}

		fclose(f);
	}
}

#endif
