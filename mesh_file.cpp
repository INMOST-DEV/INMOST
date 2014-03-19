#include "inmost.h"

#if defined(USE_MESH)
#include <deque>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>

#define R_VERSION     0
#define R_USERDATA    1
#define R_DATATYPE    2
#define R_WAITDATA    3

#define R_SGRID       4
#define R_RGRID       5
#define R_UGRID       6
#define R_POLYDATA    7
#define R_FIELD       8
#define R_SPOINTS     9

#define R_ATTRIBUTES  10
#define R_ATTRDATA    11

#define R_QUIT		  100


template<typename T>
void ReadCoords(FILE * f,INMOST_DATA_REAL_TYPE c[3])
{
	T temp[3];
	if( fread(temp,sizeof(T),3,f) == 3 )
	{
		c[0] = temp[0];
		c[1] = temp[1];
		c[2] = temp[2];
	}
	else throw INMOST::BadFile;
}

#if defined(USE_PARALLEL_WRITE_TIME)
#define REPORT_MPI(x) {WriteTab(out_time) << "<MPI><![CDATA[" << #x << "]]></MPI>\n"; x;}
#define REPORT_STR(x) {WriteTab(out_time) << "<TEXT><![CDATA[" << x << "]]></TEXT>\n";}
#define REPORT_VAL(str,x) {WriteTab(out_time) << "<VALUE name=\"" << str << "\"> <CONTENT><![CDATA[" << x << "]]> </CONTENT><CODE><![CDATA[" << #x << "]]></CODE></VALUE>\n";}
#define ENTER_FUNC() long double all_time = Timer(); WriteTab(out_time) << "<FUNCTION name=\"" << __FUNCTION__ << "\" id=\"func" << func_id++ << "\">\n"; Enter();
#define EXIT_FUNC() WriteTab(out_time) << "<TIME>" << Timer() - all_time << "</TIME>\n"; Exit(); WriteTab(out_time) << "</FUNCTION>\n";
#else
#define REPORT_MPI(x) x
#define REPORT_STR(x) 
#define REPORT_VAL(str,x)
#define ENTER_FUNC()
#define EXIT_FUNC() 
#endif



#define GMV_HEADER 0
#define GMV_NODES 1
#define GMV

namespace INMOST
{


	typedef char HeaderType;
	const HeaderType EndOfData  = 0x01;
	const HeaderType NodeHeader = 0x02;
	const HeaderType EdgeHeader = 0x03;
	const HeaderType FaceHeader = 0x04;
	const HeaderType CellHeader = 0x05;
	const HeaderType ESetHeader = 0x06;
	const HeaderType TagsHeader = 0x07;
	const HeaderType MeshHeader = 0x08;
	const HeaderType EoMHeader  = 0x09;
	const HeaderType INMOSTFile   = 0x10;
	const HeaderType MeshDataHeader = 0x11;

	
	std::ostream & operator <<(std::ostream & out, HeaderType H)
	{
		out.put(H);
		return out;
	}
	
	std::istream & operator >>(std::istream & in, HeaderType &H)
	{
		in.get(H);
		return in;
	}

	void Mesh::WriteTag(std::ostream & out, Tag X)
	{
		std::string name = X.GetTagName();
		INMOST_DATA_ENUM_TYPE namesize = name.size();
		INMOST_DATA_BULK_TYPE datatype = static_cast<INMOST_DATA_BULK_TYPE>(X.GetDataType());
		ElementType sparsemask = NONE;
		ElementType definedmask = NONE;
		INMOST_DATA_ENUM_TYPE datalength = X.GetSize();
		for(ElementType current_type = NODE; current_type <= MESH; current_type = current_type << 1)
		{
			if( X.isSparse(current_type) ) sparsemask |= current_type;
			if( X.isDefined(current_type) ) definedmask |= current_type;
		}
		uconv.write_iValue(out,namesize);
		out.write(name.c_str(), name.size());
		out.put(datatype);
		out.put(sparsemask);
		out.put(definedmask);
		uconv.write_iValue(out,datalength);
	}
	
	Tag Mesh::ReadTag(std::istream & in)
	{
		INMOST_DATA_ENUM_TYPE namesize;
		std::string name;
		char datatype;
		char sparsemask,definedmask;
		INMOST_DATA_ENUM_TYPE datalength;
		uconv.read_iValue(in,namesize);
		name.resize(namesize);
		in.read(reinterpret_cast<char *>(&name[0]), namesize);
		in.get(datatype);
		in.get(sparsemask);
		in.get(definedmask);
		uconv.read_iValue(in,datalength);
		return CreateTag(name,static_cast<DataType>(datatype),static_cast<ElementType>(definedmask),static_cast<ElementType>(sparsemask),datalength);
	}
	void Mesh::ReadData(std::istream & in, Tag t, Storage * X,std::vector<Node *> & new_nodes,std::vector<Edge *> & new_edges,std::vector<Face *> & new_faces,std::vector<Cell *> & new_cells)
	{
		
		if( !t.isDefined(X->GetElementType()) ) return;
		//if( t.isSparse(X->GetElementType()) )
		//{
		//	char have_data;
		//	in.get(have_data);
		//	if(!have_data) return;
		//}
		INMOST_DATA_ENUM_TYPE size = t.GetSize(),k,lid;
		if( size == ENUMUNDEF )
		{
			uconv.read_iValue(in,size);
			X->SetDataSize(t,size);
		}
		switch(t.GetDataType())
		{
		case DATA_REAL:      
			{
			//	INMOST_DATA_REAL_TYPE stub;
				for(k = 0; k < size; k++) 
					iconv.read_fValue(in,X->RealArray(t)[k]); 
			//			iconv.read_fValue(in,stub); 
			}
			break;
		case DATA_INTEGER:   
			{
				for(k = 0; k < size; k++) 
					iconv.read_iValue(in,X->IntegerArray(t)[k]); 
			}
			break;
		case DATA_BULK:      
			{
				in.read(reinterpret_cast<char *>(&X->Bulk(t)),size); 
			}
			break;
		case DATA_REFERENCE: 
			{
				for(k = 0; k < size; k++)
				{
					char type;
					in.get(type);
					uconv.read_iValue(in,lid);
					switch(static_cast<ElementType>(type))
					{
					case NODE: X->ReferenceArray(t)[k] = new_nodes[lid]; break;
					case EDGE: X->ReferenceArray(t)[k] = new_edges[lid]; break;
					case FACE: X->ReferenceArray(t)[k] = new_faces[lid]; break;
					case CELL: X->ReferenceArray(t)[k] = new_cells[lid]; break;
					}
					
				}
			}
			break;
		}
	}
	void Mesh::WriteData(std::ostream & out, Tag t, Storage & X)
	{
		if( !t.isDefined(X.GetElementType()) ) return;
		//if( t.isSparse(X.GetElementType()) )
		//{
		//	char have_data = X.HaveData(t);
		//	out.put(have_data);
		//	if(!have_data) return;
		//}
		INMOST_DATA_ENUM_TYPE size = t.GetSize(),k, lid;
		if( size == ENUMUNDEF )
		{
			size = X.GetDataSize(t);
			uconv.write_iValue(out,size);
		}
		switch(t.GetDataType())
		{
		case DATA_REAL:      for(k = 0; k < size; k++) iconv.write_fValue(out,X.RealArray(t)[k]); break;
		case DATA_INTEGER:   for(k = 0; k < size; k++) iconv.write_iValue(out,X.IntegerArray(t)[k]); break;
		case DATA_BULK:      out.write(reinterpret_cast<char *>(&X.Bulk(t)),size); break;
		case DATA_REFERENCE: 
			{
				for(k = 0; k < size; k++)
				{
					Element * e = X.ReferenceArray(t)[k];
					char type = e->GetElementType();
					lid = e->LocalID();
					out.put(type);
					uconv.write_iValue(out,lid);
				}
			}
		break;
		}
	}
	Node * Mesh::ReadNode(std::istream & in, std::vector<Node *> & old_nodes)
	{
		Node * ret;
		unsigned int dim = GetDimensions();
		Storage::real coords[3];
		for(unsigned int i = 0; i < dim; i++)
			iconv.read_fValue(in,coords[i]);
		int find = -1;
		if( !old_nodes.empty() ) find = binary_search(&old_nodes[0],old_nodes.size(),sizeof(Node *),CompareCoordSearch,&coords[0]);
		if( find == -1 ) 
			ret = CreateNode(&coords[0]);
		else 
			ret = old_nodes[find];
		MIDType markers;
		uconv.read_iValue(in,markers);
		ret->SetMarkerSpace(markers);
		return ret;
	}
	void Mesh::WriteNode(std::ostream & out,Node * X)
	{
		Storage::real_array coords = X->Coords();
		for(Storage::real_array::iterator it = coords.begin(); it != coords.end(); it++)
			iconv.write_fValue(out,*it);
		uconv.write_iValue(out,X->GetMarkerSpace());
	}

	Edge * Mesh::ReadEdge(std::istream & in, std::vector<Node *> & nodes)
	{
		Edge * ret;
		INMOST_DATA_ENUM_TYPE nlow,lid,i;
		dynarray<Node *,64> sub_elements;
		uconv.read_iValue(in,nlow);
		//assert(nlow == 2);
		for(i = 0; i < nlow; i++)
		{
			uconv.read_iValue(in,lid);
			sub_elements.push_back(nodes[lid]);
		}
		ret = CreateEdge(&sub_elements[0], sub_elements.size()).first;
		MIDType markers;
		uconv.read_iValue(in,markers);
		ret->SetMarkerSpace(markers);
		return ret;	
	}
	void Mesh::WriteEdge(std::ostream & out,Edge * X)
	{
		adjacent<Node> sub = X->getNodes();
		INMOST_DATA_ENUM_TYPE nlow =sub.size(),lid;
		assert(nlow == 2 );
		uconv.write_iValue(out,nlow);
		for(adjacent<Node>::iterator it = sub.begin(); it != sub.end(); it++)
		{
			lid = it->LocalID();
			uconv.write_iValue(out,lid);
		}
		uconv.write_iValue(out,X->GetMarkerSpace());
	}
	
	Face * Mesh::ReadFace(std::istream & in, std::vector<Edge *> & edges)
	{
		Face * ret;
		INMOST_DATA_ENUM_TYPE nlow,lid,i;
		dynarray<Edge *,64> sub_elements;
		uconv.read_iValue(in,nlow);
		for(i = 0; i < nlow; i++)
		{
			uconv.read_iValue(in,lid);
			sub_elements.push_back(edges[lid]);
		}
		ret = CreateFace(&sub_elements[0], sub_elements.size()).first;
		MIDType markers;
		uconv.read_iValue(in,markers);
		ret->SetMarkerSpace(markers);
		return ret;	
	}
	void Mesh::WriteFace(std::ostream & out,Face * X)
	{
		adjacent<Edge> sub = X->getEdges();
		INMOST_DATA_ENUM_TYPE nlow =sub.size(),lid;
		uconv.write_iValue(out,nlow);
		for(adjacent<Node>::iterator it = sub.begin(); it != sub.end(); it++)
		{
			lid = it->LocalID();
			uconv.write_iValue(out,lid);
		}
		uconv.write_iValue(out,X->GetMarkerSpace());
	}

	Cell * Mesh::ReadCell(std::istream & in, std::vector<Face *> & faces, std::vector<Node *> & nodes)
	{
		Cell * ret = NULL;
		INMOST_DATA_ENUM_TYPE nlow, nhigh,lid,i;
		dynarray<Face *,64> sub_elements;
		dynarray<Node *,64> suggest_nodes;
		uconv.read_iValue(in,nlow);
		for(i = 0; i < nlow; i++)
		{
			uconv.read_iValue(in,lid);
			sub_elements.push_back(faces[lid]);
		}
		uconv.read_iValue(in,nhigh);
		for(i = 0; i < nhigh; i++)
		{
			uconv.read_iValue(in,lid);
			suggest_nodes.push_back(nodes[lid]);
		}
		ret = CreateCell(&sub_elements[0],sub_elements.size(),&suggest_nodes[0],suggest_nodes.size()).first;
		MIDType markers;
		uconv.read_iValue(in,markers);
		ret->SetMarkerSpace(markers);
		return ret;	
	}
	void Mesh::WriteCell(std::ostream & out,Cell * X)
	{
		adjacent<Node> nodes = X->getNodes();
		adjacent<Face> sub = X->getFaces();
		INMOST_DATA_ENUM_TYPE nlow = sub.size(), nhigh = nodes.size(),lid;
		uconv.write_iValue(out,nlow);
		for(adjacent<Face>::iterator it = sub.begin(); it != sub.end(); it++)
		{
			lid = it->LocalID();
			uconv.write_iValue(out,lid);
		}
		uconv.write_iValue(out,nhigh);
		for(adjacent<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
		{
			lid = it->LocalID();
			uconv.write_iValue(out,lid);
		}
		uconv.write_iValue(out,X->GetMarkerSpace());
	}
	
	void Mesh::WriteElementSet(std::ostream & out, ElementSet * X)
	{
		INMOST_DATA_ENUM_TYPE size = X->size(),lid;
		char ordered = X->isOrdered();
		uconv.write_iValue(out,size);
		out.put(ordered);
		for(ElementSet::iterator it = X->begin(); it != X->end(); it++)
		{
			char type = static_cast<char>(it->GetElementType());
			lid = it->LocalID();
			out.put(type);
			uconv.write_iValue(out,lid);
		}
		uconv.write_iValue(out,X->GetMarkerSpace());
	}
	

	
	ElementSet * Mesh::ReadElementSet(std::istream & in,std::vector<Node *> & new_nodes,std::vector<Edge *> & new_edges,std::vector<Face *> & new_faces,std::vector<Cell *> & new_cells)
	{
		INMOST_DATA_ENUM_TYPE size, lid, it;
		char ordered;
		uconv.read_iValue(in,size);
		in.get(ordered);
		ElementSet * ret = ordered? CreateOrderedSet() : CreateSet();
		for(it = 0; it < size; it++)
		{
			char type;
			in.get(type);
			uconv.read_iValue(in,lid);
			switch(static_cast<ElementType>(type))
			{
				case NODE: ret->Insert(new_nodes[lid]); break;
				case EDGE: ret->Insert(new_edges[lid]); break;
				case FACE: ret->Insert(new_faces[lid]); break;
				case CELL: ret->Insert(new_cells[lid]); break;
			}
		}
		MIDType markers;
		uconv.read_iValue(in,markers);
		ret->SetMarkerSpace(markers);
		return ret;
	}



	void Mesh::Load(std::string File)
	{
		ENTER_FUNC();
		REPORT_VAL("File",File);
		std::string LFile;
		LFile.resize(File.size());
		std::transform(File.begin(),File.end(),LFile.begin(),::tolower);
//		long double load_timer = Timer();
		/*
		if(File.find("gmv") != std::string::npos) //this is gmv
		{
			char readline[2048];
			std::vector<File *> fs;
			fs.push_back(fopen(File.c_str(),"r"));
			if( fs[0] == NULL ) throw BadFileName;
			while( !fs.empty() )
			{
				while( fgets(readline,2048,fs[fs.size()-1]) != NULL )
				{
					switch(state)
					{
						case 
					}
				}
				fclose(fs[fs.size()-1]);
				fs.pop();
			}
		}
		else 
		*/
		if(LFile.find(".grdecl") != std::string::npos) // this is eclipse grid
		{
		
			Storage::real origin[2], axisx[2], axisy[2];
			std::vector<Cell *> new_cells;
			Storage::integer nx,ny,nz, i , j , k, m, q;
			std::vector< std::vector< std::vector< std::vector<Storage::real> > > > coords, zcorn;
			std::vector< std::vector< std::vector<Storage::integer> > > pinch;
			size_t l;
			std::string str;
			std::fstream stream(File.c_str(),std::ios::in);
			while(stream.good())
			{
				getline(stream,str);
				l = str.find_first_not_of(" \t\r\n");
				if( l != std::string::npos ) str.substr(l,str.size()-l).swap(str);
				else str.clear();
				l = str.find_last_not_of(" \t\r\n");
				if( l != std::string::npos ) str.resize(l+1);
				else str.clear();
				
				if( str.empty() ) continue;
				if( str[0] == '-' && str[1] == '-' ) continue;
				
				if( str == "NOECHO" ) continue;
				else if( str == "PINCH" )
				{
					getline(stream,str,'/');
				}
				else if( str == "MAPUNITS" )
				{
					getline(stream,str,'/');
				}
				else if( str == "MAPAXES" )
				{
					getline(stream,str,'/');
					std::istringstream sstream(str);
					sstream >> origin[0] >> origin[1] >> axisx[0] >> axisx[1] >> axisy[0] >> axisy[1];
				}
				else if( str == "GRIDUNIT" )
				{
					getline(stream,str,'/');
				}
				else if( str == "SPECGRID" )
				{
					getline(stream,str,'/');
					std::istringstream sstream(str);
					sstream >> nx >> ny >> nz;
					
					
					coords.resize(nx+1);
					for(i = 0; i < nx+1; i++)
					{
						coords[i].resize(ny+1);
						for(j = 0; j < ny+1; j++)
						{
							coords[i][j].resize(2);
							for(l = 0; l < 2; l++)
								coords[i][j][l].resize(3);
						}
					}
						
					zcorn.resize(nx);
					for(i = 0; i < nx; i++)
					{
						zcorn[i].resize(ny);
						for(j = 0; j < ny; j++)
						{
							zcorn[i][j].resize(nz);
							for(k = 0; k < nz; k++)
								zcorn[i][j][k].resize(8);
						}
					}
					pinch.resize(nx+1);
					for(i = 0; i < nx+1; i++)
					{
						pinch[i].resize(ny+1);
						for(j = 0; j < ny+1; j++)
							pinch[i][j].resize(nz+1);
					}

				}
				else if( str == "COORDSYS" )
				{
					getline(stream,str,'/');
				}
				else if( str == "COORD" )
				{
					for(j = 0; j < ny+1; j++)
						for(i = 0; i < nx+1; i++)
							for(l = 0; l < 2; l++)
								for(k = 0; k < 3; k++)
									stream >> coords[i][j][l][k];
					getline(stream,str,'/');
				}
				else if( str == "ZCORN" )
				{
					for(k = 0; k < nz; k++)
					{
						for(q = 0; q < 2; q++)
							for(j = 0; j < ny; j++)
								for(m = 0; m < 2; m++)
									for(i = 0; i < nx; i++)
										for(l = 0; l < 2; l++)
											stream >> zcorn[i][j][k][l+m*2+(1-q)*4];
					}
					getline(stream,str,'/');
				}
				else if( str == "ACTNUM" )
				{
					pinch.resize(nx+1);
					for(k = 0; k < nz+1; k++)
						for(j = 0; j < ny+1; j++)
							for(i = 0; i < nx+1; i++)
								stream >> pinch[i][j][k];
					getline(stream,str,'/');
				}
				else
				{
					std::cout << "Ignoring keyword " << str << std::endl;
					getline(stream,str,'/');
				}
			}
			{
				std::vector< std::vector < std::map< Storage::real, Node * > > > nodes;
				nodes.resize(nx+1);
				for(i = 0; i < nx+1; i++) nodes[i].resize(ny+1);
				
				//nx = ny = 15;
				//nz = 80;
				
				//std::cout << "Creating grid" << std::endl;
				
				Node * c_nodes[8];
				int cbad = 0, cflat = 0;
				for(i = 0; i < nx; i++)
				{
					for(j = 0; j < ny; j++)
					{
						for(k = 0; k < nz; k++)
						{
							for(l = 0; l < 8; l++)
							{
								int i2 = i + l%2, j2 = j + (l/2)%2;
								std::map<Storage::real, Node * >::iterator search = nodes[i2][j2].find(zcorn[i][j][k][l]);
								if( search == nodes[i2][j2].end() )
								{
									Storage::real node_coords[3];
									Storage::real alpha = (zcorn[i][j][k][l]-coords[i2][j2][1][2])/(coords[i2][j2][0][2]-coords[i2][j2][1][2]);
									node_coords[0] = coords[i2][j2][1][0] + alpha * (coords[i2][j2][0][0] - coords[i2][j2][1][0]);
									node_coords[1] = coords[i2][j2][1][1] + alpha * (coords[i2][j2][0][1] - coords[i2][j2][1][1]);
									node_coords[2] = zcorn[i][j][k][l];
									search = nodes[i2][j2].insert(std::pair<Storage::real,Node *>(zcorn[i][j][k][l],CreateNode(node_coords))).first;
								}
								c_nodes[l] = search->second;
							}
							INMOST_DATA_ENUM_TYPE nodesnum[24] = {2,3,1,0,5,7,6,4,4,6,2,0,3,7,5,1,1,5,4,0,6,7,3,2};
							INMOST_DATA_ENUM_TYPE facessize[6] = {4,4,4,4,4,4};
							//const int nodesnum[6][4] = {{0,1,3,2},{4,6,7,5},{0,2,6,4},{1,5,7,3},{0,4,5,1},{2,3,7,6}};
							new_cells.push_back(CreateCell(c_nodes,nodesnum,facessize,6).first);
						}
					}
				}
				std::cout << "bad cells: " << cbad << " flat cells: " << cflat << std::endl;
				//std::cout << "Done" << std::endl;
			}
			stream.close();
		}
		else
		if(LFile.find(".pvtk") != std::string::npos) //this is legacy parallel vtk
		{
			int state = 0, np, nf = 0;
			size_t l,attrl, ql,qr;
			//~ if( m_state == Mesh::Serial ) SetCommunicator(INMOST_MPI_COMM_WORLD);
			std::string str, tag, attrval, path = "";
			std::vector<std::string> files;
			std::fstream stream(File.c_str(),std::ios::in);
			l = File.find_last_of("/\\");
			if( l != std::string::npos )
				path = File.substr(0,l+1);
			while(stream.good())
			{
				getline(stream,str,'>');
				l = str.find_first_not_of(" \t\r\n");
				if( l != std::string::npos ) str.substr(l,str.size()-l).swap(str);
				else str.clear();
				if( str.empty() ) continue;
				l = str.find_first_of(' ');
				tag = str.substr(1,l-1);
				//std::cout << " tag name " << tag << std::endl;
				if( state == 0 )
				{
					if( tag == "File" )
					{
						attrl = str.find("version=\"pvtk-1.0\"",l);
						if( attrl == std::string::npos ) throw BadFile;
						attrl = str.find("dataType=\"vtkUnstructuredGrid\"",l);
						if( attrl == std::string::npos ) throw BadFile;
						attrl = str.find("numberOfPieces");
						if( attrl == std::string::npos ) throw BadFile;
						ql = str.find_first_of('"',attrl);
						qr = str.find_first_of('"',ql+1);
						attrval = str.substr(ql+1,qr-ql-1);
						//std::cout << "number of peaces: " << attrval << " l: " << ql << " r: " << qr << std::endl;
						np = atoi(attrval.c_str());
						//std::cout << "np: " << np << std::endl;;
						files.resize(np);
						state = 1;
					}
					else throw BadFile;
				}
				else if( state == 1 )
				{
					if( tag == "Piece")
					{
						attrl = str.find("fileName");
						ql = str.find_first_of('"',attrl);
						qr = str.find_first_of('"',ql+1);
						attrval = str.substr(ql+1,qr-ql-1);
						//std::cout << "file" << nf << ":" << attrval << std::endl;
						files[nf++] = attrval;
						if( nf == np ) state = 2;
					}
					else throw BadFile;
				}
				else if( state == 2 )
				{
					if( tag != "/File" ) throw BadFile;
					else break;
				}
			}
			for(int i = GetProcessorRank(); i < static_cast<int>(files.size()); i+= GetProcessorsNumber()) 
			{
				//std::cout << "load " << i << ": " << path + files[i] << " by " << GetProcessorRank() << std::endl;
				Load(path + files[i]);
			}
			ResolveShared();
		}
		else if(LFile.find(".vtk") != std::string::npos) //this is legacy vtk
		{
			MIDType unused_marker = CreateMarker();
			std::vector<Node *> old_nodes;
			for(nodes_container::iterator it = nodes.begin(); it != nodes.end(); it++) if( *it != NULL )
				old_nodes.push_back(*it);

			if( !old_nodes.empty() ) qsort(&old_nodes[0],old_nodes.size(),sizeof(Node *),CompareElementsCCentroid);
			
			FILE * f = fopen(File.c_str(),"r");
			if( !f ) throw BadFileName;
			std::vector<Node *> newnodes;
			std::vector<Cell *> newcells;
			std::vector<int> cp;
			std::vector<int> ct;
			unsigned int state = R_VERSION;
			char readline[2048];
			int filled;
			bool binary = false;
			int read_into = NONE;
			while( state != R_QUIT )
			{
				
				if( fgets(readline,2048,f) == NULL )
				{
					state = R_QUIT;
					continue;
				}
				if( readline[strlen(readline)-1] == '\n' )
					readline[strlen(readline)-1] = '\0';
				if( strlen(readline) == 0 ) continue;
				switch( state )
				{
					case R_VERSION:
					{
						int h,l;
						filled = sscanf(readline,"# vtk DataFile Version %1d.%1d",&h,&l);
						if( filled == 0 ) throw BadFile;
						//file version checker
						state = R_USERDATA;
						break;
					}
					case R_USERDATA:
					{
						//printf("Attached info: %s\n",readline);
						state = R_DATATYPE;
						break;
					}
					case R_DATATYPE:
					{
						if( !strncmp(readline,"BINARY",6) )
							binary = 1;
						else if( !strncmp(readline,"ASCII",5) )
							state = R_WAITDATA;
						else throw BadFile;
						break;
					}
					case R_ATTRDATA:
					{
						char dataname[1024];
						char attrname[1024];
						char attrtype[1024];
						int nentries = 1;
						filled = sscanf(readline," %s ",dataname);
						if( filled != 1 ) throw BadFile;
						if( !strcmp(dataname,"SCALARS") )
						{
							DataType t = DATA_BULK;
							filled = sscanf(readline," SCALARS %s %s %d",attrname,attrtype,&nentries);
							if( filled < 2 ) {printf("found %d arguments to SCALARS field, must be >= 2\nline:\n%s\n",filled,readline); throw BadFile;}
							for(unsigned int i = 0; i < strlen(attrtype); i++) attrtype[i] = tolower(attrtype[i]);
							if( !strcmp(attrtype,"bit") || !strcmp(attrtype,"unsigned_char") || !strcmp(attrtype,"char") ||
							   !strcmp(attrtype,"unsigned_short") || !strcmp(attrtype,"short") || !strcmp(attrtype,"unsigned_int") ||
							   !strcmp(attrtype,"int") || !strcmp(attrtype,"unsigned_long") || !strcmp(attrtype,"long"))
								t = DATA_INTEGER;
							else if( !strcmp(attrtype,"float") || !strcmp(attrtype,"double") )
								t = DATA_REAL;
							else
								throw BadFile;
							Tag attr = CreateTag(attrname,t,read_into,NONE,nentries);
							if( fgets(readline,2048,f) == NULL ) throw BadFile; //LOOK_UP TABLE
							if( read_into == CELL )
							{
								for(unsigned int it = 0; it < newcells.size(); it++)
								{
									if( t == DATA_INTEGER )
									{
										if( newcells[it] != NULL )
										{
											Storage::integer_array attrdata = newcells[it]->IntegerArray(attr);
											for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
										}
										else for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
									}
									if( t == DATA_REAL )
									{
										if( newcells[it] != NULL )
										{
											Storage::real_array attrdata = newcells[it]->RealArray(attr);
											for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
										}
										else for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; }
									}
								}
								filled = fscanf(f,"\n");
							}
							if( read_into == NODE )
							{
								for(unsigned int it = 0; it < newnodes.size(); it++)
								{
									if( t == DATA_INTEGER )
									{
										Storage::integer_array attrdata = newnodes[it]->IntegerArray(attr);
										for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
									if( t == DATA_REAL )
									{
										Storage::real_array attrdata = newnodes[it]->RealArray(attr);
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
								}
								filled = fscanf(f,"\n");
							}
							break;
						}
						else if( !strcmp(dataname,"COLOR_SCALARS") )
						{	
							char attrname[1024];
							filled = sscanf(readline," COLOR_SCALARS %s %d",attrname,&nentries);
							Tag attr = CreateTag(attrname,DATA_REAL,read_into,NONE,nentries);
							if( read_into == CELL )
							{
								for(unsigned int it = 0; it < newcells.size(); it++)
								{
									if( newcells[it] != NULL )
									{
										Storage::real_array attrdata = newcells[it]->RealArray(attr);
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
									else for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
								}
								filled = fscanf(f,"\n");
							}
							if( read_into == NODE )
							{
								for(unsigned int it = 0; it < newnodes.size(); it++)
								{
									Storage::real_array attrdata = newnodes[it]->RealArray(attr);
									for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
								}
								filled = fscanf(f,"\n");
							}
							break;
						}
						else if( !strcmp(dataname," LOOKUP_TABLE") )
						{
							int nentries;
							filled = sscanf(readline," LOOKUP_TABLE %*s %d",&nentries);
							if( filled != 1 ) throw BadFile;
							for(int i = 0; i < nentries; i++)
							{
								filled = fscanf(f," %*f %*f %*f %*f");
							}
							filled = fscanf(f,"\n");
							break;
						}
						else if( !strcmp(dataname,"NORMALS") || !strcmp(dataname,"VECTORS") || !strcmp(dataname,"TENSORS"))
						{
							char attrname[1024];
							char attrtype[1024];
							int nentries = 3;
							if( !strcmp(dataname,"TENSORS") ) nentries = 9;
							DataType t = DATA_BULK;
							filled = sscanf(readline,"%*s %s %s",attrname,attrtype);
							if( filled != 2 ) throw BadFile;
							for(unsigned int i = 0; i < strlen(attrtype); i++) attrtype[i] = tolower(attrtype[i]);
							if( !strcmp(attrtype,"bit") || !strcmp(attrtype,"unsigned_char") || !strcmp(attrtype,"char") ||
							   !strcmp(attrtype,"unsigned_short") || !strcmp(attrtype,"short") || !strcmp(attrtype,"unsigned_int") ||
							   !strcmp(attrtype,"int") || !strcmp(attrtype,"unsigned_long") || !strcmp(attrtype,"long"))
								t = DATA_INTEGER;
							else if( !strcmp(attrtype,"float") || !strcmp(attrtype,"double") )
								t = DATA_REAL;
							else
								throw BadFile;
							Tag attr = CreateTag(attrname,t,read_into,NONE,nentries);
							if( fgets(readline,2048,f) == NULL ) throw BadFile; //LOOK_UP TABLE
							if( read_into == CELL )
							{
								for(unsigned int it = 0; it < newcells.size(); it++)
								{
									if( t == DATA_INTEGER )
									{
										if( newcells[it] != NULL )
										{
											Storage::integer_array attrdata = newcells[it]->IntegerArray(attr);
											for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
										}
										else for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
									}
									if( t == DATA_REAL )
									{
										if( newcells[it] != NULL )
										{
											Storage::real_array attrdata = newcells[it]->RealArray(attr);
											for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
										}
										else for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
									}
								}
								filled = fscanf(f,"\n");
							}
							if( read_into == NODE )
							{
								for(unsigned int it = 0; it < newnodes.size(); it++)
								{
									if( t == DATA_INTEGER )
									{
										Storage::integer_array attrdata = newnodes[it]->IntegerArray(attr);
										for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
									if( t == DATA_REAL )
									{
										Storage::real_array attrdata = newnodes[it]->RealArray(attr);
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
								}
								filled = fscanf(f,"\n");
							}
							break;
						}
						else if( !strcmp(dataname,"TEXTURE_COORDINATES") )
						{
							char attrname[1024];
							char attrtype[1024];
							int nentries;
							DataType t = DATA_BULK;
							filled = sscanf(readline," TEXTURE_COORDINATES %s %d %s",attrname,&nentries,attrtype);
							if( filled != 3 ) throw BadFile;
							for(unsigned int i = 0; i < strlen(attrtype); i++) attrtype[i] = tolower(attrtype[i]);
							if( !strcmp(attrtype,"bit") || !strcmp(attrtype,"unsigned_char") || !strcmp(attrtype,"char") ||
							   !strcmp(attrtype,"unsigned_short") || !strcmp(attrtype,"short") || !strcmp(attrtype,"unsigned_int") ||
							   !strcmp(attrtype,"int") || !strcmp(attrtype,"unsigned_long") || !strcmp(attrtype,"long"))
								t = DATA_INTEGER;
							else if( !strcmp(attrtype,"float") || !strcmp(attrtype,"double") )
								t = DATA_REAL;
							else
								throw BadFile;
							Tag attr = CreateTag(attrname,t,read_into,NONE,nentries);
							if( fgets(readline,2048,f) == NULL ) throw BadFile; //LOOK_UP TABLE
							if( read_into == CELL )
							{
								for(unsigned int it = 0; it < newcells.size(); it++)
								{
									if( t == DATA_INTEGER )
									{
										if( newcells[it] != NULL )
										{
											Storage::integer_array attrdata = newcells[it]->IntegerArray(attr);
											for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
										} else for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
									}
									if( t == DATA_REAL )
									{
										if( newcells[it] != NULL )
										{
											Storage::real_array attrdata = newcells[it]->RealArray(attr);
											for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
										} else for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
									}
								}
								filled = fscanf(f,"\n");
							}
							if( read_into == NODE )
							{
								for(unsigned int it = 0; it < newnodes.size(); it++)
								{
									if( t == DATA_INTEGER )
									{
										Storage::integer_array attrdata = newnodes[it]->IntegerArray(attr);
										for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
									if( t == DATA_REAL )
									{
										Storage::real_array attrdata = newnodes[it]->RealArray(attr);
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
								}
								filled = fscanf(f,"\n");
							}
							break;
						}
						else if( !strcmp(dataname,"FIELD") )
						{
							int nfields = 0;
							filled = sscanf(readline," FIELD %*s %d",&nfields);
							if( filled != 1 ) throw BadFile;
							for(int i = 0; i < nfields; i++)
							{
								char attrname[1024];
								char attrtype[1024];
								unsigned int nentries, ntuples;
								DataType t = DATA_BULK;
								filled = fscanf(f," %s %u %u %s\n",attrname,&nentries,&ntuples,attrtype);
								if( filled != 4 ) throw BadFile;
								for(unsigned int i = 0; i < strlen(attrtype); i++) attrtype[i] = tolower(attrtype[i]);
								if( !strcmp(attrtype,"bit") || !strcmp(attrtype,"unsigned_char") || !strcmp(attrtype,"char") ||
								   !strcmp(attrtype,"unsigned_short") || !strcmp(attrtype,"short") || !strcmp(attrtype,"unsigned_int") ||
								   !strcmp(attrtype,"int") || !strcmp(attrtype,"unsigned_long") || !strcmp(attrtype,"long"))
									t = DATA_INTEGER;
								else if( !strcmp(attrtype,"float") || !strcmp(attrtype,"double") )
									t = DATA_REAL;
								else
									throw BadFile;
								Tag attr = CreateTag(attrname,t,read_into,NONE,nentries);
								
								if( read_into == CELL )
								{
									if( ntuples != newcells.size() ) printf("number of tuples in field is not equal to number of cells\n");
									for(unsigned int it = 0; it < ntuples; it++)
									{
										if( t == DATA_INTEGER )
										{
											if( newcells[it] != NULL )
											{
												Storage::integer_array attrdata = newcells[it]->IntegerArray(attr);
												for(unsigned int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
											} else for(unsigned int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
										}
										if( t == DATA_REAL )
										{
											if( newcells[it] != NULL )
											{
												Storage::real_array attrdata = newcells[it]->RealArray(attr);
												for(unsigned int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
											} else for(unsigned int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
										}
									}
									filled = fscanf(f,"\n");
								}
								if( read_into == NODE )
								{
									if( ntuples != newcells.size() ) printf("number of tuples in field is not equal to number of nodes\n");
									for(unsigned int it = 0; it < ntuples; it++)
									{
										if( t == DATA_INTEGER )
										{
											Storage::integer_array attrdata = newnodes[it]->IntegerArray(attr);
											for(unsigned int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
										}
										if( t == DATA_REAL )
										{
											Storage::real_array attrdata = newnodes[it]->RealArray(attr);
											for(unsigned int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
										}
									}
									filled = fscanf(f,"\n");
								}
							}
							break;
						}
						else
						{
							state = R_ATTRIBUTES;
						}
					}
					case R_ATTRIBUTES:
					{
						char datatype[1024];
						unsigned int npoints;
						filled = sscanf(readline," %s %u",datatype,&npoints);
						if( filled != 2 ) throw BadFile;
						if( !strcmp(datatype,"CELL_DATA") )
						{
							if( npoints != newcells.size() ) printf("number of attributes is not equal to number of cells\n");
							state = R_ATTRDATA;
							read_into = CELL;
							break;
						}
						else if( !strcmp(datatype,"POINT_DATA") )
						{
							if( npoints != newnodes.size() ) printf("number of attributes is not equal to number of nodes\n");
							state = R_ATTRDATA;
							read_into = NODE;
							break;
						}
						else
						{
							printf("%s\n",readline);
							printf("Unknown type of attributes\n");
							state = R_WAITDATA;
						}
					}
					case R_WAITDATA:
					{
						char check[1024];
						filled = sscanf(readline," DATASET %s",check);
						if( filled != 1 ) throw BadFile;
						if( !strcmp(check,"STRUCTURED_POINTS") )
							state = R_SPOINTS;
						else if( !strcmp(check,"STRUCTURED_GRID" ) )
							state = R_SGRID;
						else if( !strcmp(check,"UNSTRUCTURED_GRID") )
							state = R_UGRID;
						else if( !strcmp(check,"POLYDATA") )
							state = R_POLYDATA;
						else if( !strcmp(check,"RECTLINEAR_GRID") )
							state = R_RGRID;
						else if( !strcmp(check,"FIELD") )
							state = R_FIELD;
						else throw BadFile;
						break;
					}
					case R_UGRID:
					{
						INMOST_DATA_REAL_TYPE coords[3] = {0,0,0};
						int npoints = 0, i = 0, ncells = 0, ncells2 = 0, nints = 0;
						char datatype[1024];
						filled = sscanf(readline,"%*s %d %s",&npoints,datatype);
						//printf("number of nodes: %d\n",npoints);
						if( filled != 2 ) throw BadFile;
						if( binary )
						{
							i = 0;
							while(i != npoints)
							{
								if( !strcmp(datatype,"double") ) ReadCoords<double>(f,coords);
								else if( !strcmp(datatype,"float") ) ReadCoords<float>(f,coords);
								else if( !strcmp(datatype,"unsigned char") ) ReadCoords<unsigned char>(f,coords);
								else if( !strcmp(datatype,"char") ) ReadCoords<char>(f,coords);
								else if( !strcmp(datatype,"unsigned short") ) ReadCoords<unsigned short>(f,coords);
								else if( !strcmp(datatype,"short") ) ReadCoords<short>(f,coords);
								else if( !strcmp(datatype,"unsigned int") ) ReadCoords<unsigned int>(f,coords);
								else if( !strcmp(datatype,"int") ) ReadCoords<int>(f,coords);
								else if( !strcmp(datatype,"unsigned long") ) ReadCoords<unsigned long>(f,coords);
								else if( !strcmp(datatype,"long") ) ReadCoords<long>(f,coords);
								else throw BadFile;
								int find = -1;
								if( !old_nodes.empty() )
									find = binary_search(&old_nodes[0],old_nodes.size(),sizeof(Element*),CompareCoordSearch,coords);
								if( find == -1 ) newnodes.push_back(CreateNode(coords));
								else newnodes.push_back(old_nodes[find]);
								newnodes.back()->SetMarker(unused_marker);
								i++;
							}
						}
						else
						{
							i = 0;
							while(i != npoints)
							{
								double c[3];
								filled = fscanf(f," %lg %lg %lg",&c[0],&c[1],&c[2]);
								coords[0] = c[0];
								coords[1] = c[1];
								coords[2] = c[2];
								if( filled != 3 ) throw BadFile;
								int find = -1;
								if( !old_nodes.empty() )
									find = binary_search(&old_nodes[0],old_nodes.size(),sizeof(Element*),CompareCoordSearch,coords);
								if( find == -1 ) newnodes.push_back(CreateNode(coords));
								else newnodes.push_back(old_nodes[find]);
								newnodes.back()->SetMarker(unused_marker);
								i++;
							}
							filled = fscanf(f,"\n");
						}
						if( fgets(readline,2048,f) == NULL ) throw BadFile;
						filled = sscanf(readline,"%*s %d %d\n",&ncells,&nints);
						//printf("number of cells: %d\nnumber of ints: %d\n",ncells,nints);
						if( filled != 2 ) throw BadFile;
						cp.resize(nints);
						ct.resize(ncells);
						if( binary )
							filled = fread(&cp[0],sizeof(int),nints,f);
						else
						{
							i = 0;
							while(i != nints)
							{
								filled = fscanf(f," %d",&cp[i]);
								if( filled != 1 ) throw BadFile;
								i++;
							}
							filled = fscanf(f,"\n");
						}
						if( fgets(readline,2048,f) == NULL ) throw BadFile;
						filled = sscanf(readline,"%*s %d\n",&ncells2);
						if( filled != 1 ) throw BadFile;
						if( ncells2 != ncells ) throw BadFile;
						if( binary )
							filled = fread(&ct[0],sizeof(int),nints,f);
						else
						{
							i = 0;
							while(i != ncells)
							{
								filled = fscanf(f," %d",&ct[i]);
								if( filled != 1 ) throw BadFile;
								i++;
							}
							filled = fscanf(f,"\n");
						}
						{
							int j = 0;
							dynarray<Node *,256> c_nodes;
							Node * e_nodes[2];
							dynarray<Edge *,64> f_edges;
							dynarray<Face *,64> c_faces;
							dynarray<Edge *,64> v_edges;
							for(i = 0; i < ncells; i++)
							{
								//printf("load progress: %20.2f%%\r",(float)(i+1)/(float)ncells*100.0f);
								fflush(stdin);
								c_nodes.clear();
								f_edges.clear();
								c_faces.clear();
								v_edges.clear();
								switch(ct[i])
								{
									case 1: // VTK_VERTEX
									{
										printf("VTK_VERTEX found, but 0d cell objects are not supported\n");
										newcells.push_back(NULL);
										break;
									}
									case 2: //VTK_POLY_VERTEX
									{
										printf("VTK_POLY_VERTEX cannot be represented\n");
										break;
									}
									case 3: //VTK_LINE
									{
										printf("VTK_LINE found, but 1d cell objects are not supported\n");
										newcells.push_back(NULL);
										break;
									}
									case 4: //VTK_POLY_LINE
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											newnodes[cp[k]]->RemMarker(unused_marker);
										}
										j = j + 1 + cp[j];
										if( c_nodes.size() < 2 ) throw BadFile;
										f_edges.resize(2);
										for(unsigned int k = 0; k < c_nodes.size()-1; k++)
										{
											e_nodes[0] = c_nodes[k];
											f_edges[0] = CreateEdge(e_nodes, 1).first;
											e_nodes[0] = c_nodes[(k+1)];
											f_edges[1] = CreateEdge(e_nodes, 1).first;
											c_faces.push_back(CreateFace(&f_edges[0],2).first);
										}
										Cell * c = CreateCell(&c_faces[0],c_faces.size(),&c_nodes[0],c_nodes.size()).first;
										c_faces.clear();
										newcells.push_back(c);
										break;
									}
									case 5: //VTK_TRIANGLE
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											newnodes[cp[k]]->RemMarker(unused_marker);
										}
										j = j + 1 + cp[j];
										if( c_nodes.size() != 3 ) throw BadFile;
										f_edges.resize(2);
										for(unsigned int k = 0; k < 3; k++)
										{
											e_nodes[0] = c_nodes[k];
											f_edges[0] = CreateEdge(e_nodes, 1).first;
											e_nodes[0] = c_nodes[(k+1)%3];
											f_edges[1] = CreateEdge(e_nodes, 1).first;
											c_faces.push_back(CreateFace(&f_edges[0],2).first);
										}
										Cell * c = CreateCell(&c_faces[0],2,&c_nodes[0],c_nodes.size()).first;
										c_faces.clear();
										newcells.push_back(c);
										break;
									}
									case 6: //VTK_TRIANGLE_STRIP
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											newnodes[cp[k]]->RemMarker(unused_marker);
										}
										j = j + 1 + cp[j];
										if( c_nodes.size() < 3 ) throw BadFile;
										std::deque<Node *> d_nodes;
										unsigned int l = 3;
										f_edges.resize(3);
										d_nodes.insert(d_nodes.begin(),c_nodes.begin(),c_nodes.begin()+3);
										do
										{
											for(int q = 0; q < 3; q++)
											{
												e_nodes[0] = d_nodes[q];
												e_nodes[1] = d_nodes[(q+1)%3];
												f_edges[q] = CreateEdge(e_nodes, 2).first;
											}
											if( l == c_nodes.size() ) break;
											c_faces.push_back(CreateFace(&f_edges[0],3).first);
											d_nodes.pop_front();
											d_nodes.push_back(c_nodes[l]);
											l++;
										}
										while(l <= c_nodes.size());
										Cell * c = CreateCell(&c_faces[0],c_faces.size(),&c_nodes[0],c_nodes.size()).first;
										newcells.push_back(c);
										break;
									}
									case 7: //VTK_POLYGON
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											newnodes[cp[k]]->RemMarker(unused_marker);
										}
										j = j + 1 + cp[j];
										f_edges.resize(2);
										for(unsigned int k = 0; k < c_nodes.size(); k++)
										{
											e_nodes[0] = c_nodes[k];
											f_edges[0] = CreateEdge(e_nodes, 1).first;
											e_nodes[0] = c_nodes[(k+1)%c_nodes.size()];
											f_edges[1] = CreateEdge(e_nodes, 1).first;
											c_faces.push_back(CreateFace(&f_edges[0],2).first);
										}
										Cell * c = CreateCell(&c_faces[0],c_faces.size(),&c_nodes[0],c_nodes.size()).first;
										c_faces.clear();
										newcells.push_back(c);
										break;
									}
									case 8: //VTK_PIXEL
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											newnodes[cp[k]]->RemMarker(unused_marker);
										}
										j = j + 1 + cp[j];
										if( c_nodes.size() != 4 ) throw BadFile;
										Node * temp = c_nodes[2];
										c_nodes[2] = c_nodes[3];
										c_nodes[3] = temp;
										f_edges.resize(2);
										for(int k = 0; k < 4; k++)
										{
											e_nodes[0] = c_nodes[k];
											f_edges[0] = CreateEdge(e_nodes, 1).first;
											e_nodes[0] = c_nodes[(k+1)%4];
											f_edges[1] = CreateEdge(e_nodes, 1).first;
											c_faces.push_back(CreateFace(&f_edges[0],2).first);
										}
										Cell * c = CreateCell(&c_faces[0],c_faces.size(),&c_nodes[0],4).first;
										newcells.push_back(c);
										break;
									}
									case 9: //VTK_QUAD
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											newnodes[cp[k]]->RemMarker(unused_marker);
										}
										j = j + 1 + cp[j];
										if( c_nodes.size() != 4 ) throw BadFile;
										f_edges.resize(2);
										for(int k = 0; k < 4; k++)
										{
											e_nodes[0] = c_nodes[k];
											f_edges[0] = CreateEdge(e_nodes, 1).first;
											e_nodes[0] = c_nodes[(k+1)%4];
											f_edges[1] = CreateEdge(e_nodes, 1).first;
											c_faces.push_back(CreateFace(&f_edges[0],2).first);
										}
										Cell * c = CreateCell(&c_faces[0],c_faces.size(),&c_nodes[0],4).first;
										c_faces.clear();
										newcells.push_back(c);
										break;
									}
									case 10: //VTK_TETRA
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											newnodes[cp[k]]->RemMarker(unused_marker);
										}
										j = j + 1 + cp[j];
										if( c_nodes.size() != 4 ) throw BadFile;
										INMOST_DATA_ENUM_TYPE nodesnum[12] = {0,2,1,0,1,3,1,2,3,0,3,2};
										INMOST_DATA_ENUM_TYPE sizes[4] = {3,3,3,3};
										Cell * c = CreateCell(&c_nodes[0],nodesnum,sizes,4).first;
										newcells.push_back(c);
										break;
									}
									case 11: //VTK_VOXEL
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											newnodes[cp[k]]->RemMarker(unused_marker);
										}
										j = j + 1 + cp[j];
										if (c_nodes.size() != 8) throw BadFile;
										{
											Node * temp;
											temp = c_nodes[2];
											c_nodes[2] = c_nodes[3];
											c_nodes[3] = temp;
											temp = c_nodes[6];
											c_nodes[6] = c_nodes[7];
											c_nodes[7] = temp;
										}
										//INMOST_DATA_ENUM_TYPE nodesnum[24] = {0,4,6,2,1,3,7,5,2,6,7,3,0,1,5,4,0,2,3,1,4,5,7,6};
										INMOST_DATA_ENUM_TYPE nodesnum[24] = {0,4,7,3,1,2,6,5,0,1,5,4,3,7,6,2,0,3,2,1,4,5,6,7};
										INMOST_DATA_ENUM_TYPE sizes[6] = {4,4,4,4,4,4};
										Cell * c = CreateCell(&c_nodes[0],nodesnum,sizes,6).first;
										newcells.push_back(c);
										break;
									}
									case 12: //VTK_HEXAHEDRON
									{

										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											newnodes[cp[k]]->RemMarker(unused_marker);
										}
										j = j + 1 + cp[j];
										if( c_nodes.size() != 8 ) throw BadFile;
										INMOST_DATA_ENUM_TYPE nodesnum[24] = {0,4,7,3,1,2,6,5,0,1,5,4,3,7,6,2,0,3,2,1,4,5,6,7};
										INMOST_DATA_ENUM_TYPE sizes[6] = {4,4,4,4,4,4};
										Cell * c = CreateCell(&c_nodes[0],nodesnum,sizes,6).first;
										newcells.push_back(c);
										break;
									}
									case 13: //VTK_WEDGE
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											newnodes[cp[k]]->RemMarker(unused_marker);
										}
										j = j + 1 + cp[j];
										if( c_nodes.size() != 6 ) throw BadFile;
										INMOST_DATA_ENUM_TYPE nodesnum[18] = {0,2,5,3,1,4,5,2,0,3,4,1,3,5,4,0,1,2};
										INMOST_DATA_ENUM_TYPE sizes[5] = {4,4,4,3,3};
										Cell * c = CreateCell(&c_nodes[0],nodesnum,sizes,5).first;
										newcells.push_back(c);
										break;
									}
									case 14: //VTK_PYRAMID
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											newnodes[cp[k]]->RemMarker(unused_marker);
										}
										j = j + 1 + cp[j];
										if( c_nodes.size() != 5 ) throw BadFile;
										INMOST_DATA_ENUM_TYPE nodesnum[16] = {0,4,3,0,1,4,1,2,4,3,4,2,0,3,2,1};
										INMOST_DATA_ENUM_TYPE sizes[5] = {3,3,3,3,4};
										Cell * c = CreateCell(&c_nodes[0],nodesnum,sizes,5).first;
										newcells.push_back(c);
										break;
									}
									case 15: //VTK_PENTAGONAL_PRISM
									{
										printf("no order description for VTK_PENTAGONAL_PRISM\n");
										break;
									}
									case 16: //VTK_HEXAGONAL_PRISM
									{
										printf("no order description for VTK_HEXAGONAL_PRISM\n");
										break;
									}
									case 21: //VTK_QUADRATIC_EDGE
									{
										printf("no order description for VTK_QUADRATIC_EDGE\n");
										break;
									}
									case 22: //VTK_QUADRATIC_TRIANGLE
									{
										printf("no order description for VTK_QUADRATIC_TRIANGLE\n");
										break;
									}
									case 23: //VTK_QUADRATIC_QUAD
									{
										printf("no order description for VTK_QUADRATIC_QUAD\n");
										break;
									}
									case 24: //VTK_QUADRATIC_TETRA
									{
										printf("no order description for VTK_QUADRATIC_TETRA\n");
										break;
									}
									case 25: //VTK_QUADRATIC_HEXAHEDRON
									{
										printf("no order description for VTK_QUADRATIC_HEXAHEDRON\n");
										break;
									}
									case 26: //VTK_QUADRATIC_WEDGE
									{
										printf("no order description for VTK_QUADRATIC_WEDGE\n");
										break;
									}
									case 27: //VTK_QUADRATIC_PYRAMID
									{
										printf("no order description for VTK_QUADRATIC_PYRAMID\n");
										break;
									}
									case 28: //VTK_BIQUADRATIC_QUAD
									{
										printf("no order description for VTK_BIQUADRATIC_QUAD\n");
										break;
									}
									case 29: //VTK_TRIQUADRATIC_HEXAHEDRON
									{
										printf("no order description for VTK_TRIQUADRATIC_HEXAHEDRON\n");
										break;
									}
									case 30: //VTK_QUADRATIC_LINEAR_QUAD
									{
										printf("no order description for VTK_QUADRATIC_LINEAR_QUAD\n");
										break;
									}
									case 31: //VTK_QUADRATIC_LINEAR_WEDGE
									{
										printf("no order description for VTK_QUADRATIC_LINEAR_WEDGE\n");
										break;
									}
									case 32: //VTK_BIQUADRATIC_QUADRATIC_WEDGE
									{
										printf("no order description for VTK_BIQUADRATIC_QUADRATIC_WEDGE\n");
										break;
									}
									case 33: //VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON
									{
										printf("no order description for VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON\n");
										break;
									}
									case 34: //VTK_BIQUADRATIC_TRIANGLE
									{
										printf("no order description for VTK_BIQUADRATIC_TRIANGLE\n");
										break;
									}
									case 35: //VTK_CUBIC_LINE
									{
										printf("no order description for VTK_CUBIC_LINE\n");
										break;
									}
									case 41: //VTK_CONVEX_POINT_SET
									{
										printf("no algorithm for VTK_CONVEX_POINT_SET yet\n");
										break;
									}
									case 42: //VTK_POLYHEDRON
									{
										int k = j;
										k++; //skip number of total number of integers
										dynarray<INMOST_DATA_ENUM_TYPE,64> sizes(cp[k++]);
										for(INMOST_DATA_ENUM_TYPE m = 0; m < sizes.size(); m++)
										{
											sizes[m] = cp[k++];
											for (INMOST_DATA_ENUM_TYPE l = 0; l < sizes[m]; l++)
											{
												c_nodes.push_back(newnodes[cp[k++]]);
												c_nodes.back()->RemMarker(unused_marker);
											}
										}
										j = k;
										Cell * c = CreateCell(&c_nodes[0],&sizes[0],sizes.size()).first;
										newcells.push_back(c);
										break;
									}
									case 51: //VTK_PARAMETRIC_CURVE
									{
										printf("no order description for VTK_PARAMETRIC_CURVE\n");
										break;
									}
									case 52: //VTK_PARAMETRIC_SURFACE
									{
										printf("no order description for VTK_PARAMETRIC_SURFACE\n");
										break;
									}
									case 53: //VTK_PARAMETRIC_TRI_SURFACE
									{
										printf("no order description for VTK_PARAMETRIC_TRI_SURFACE\n");
										break;
									}
									case 54: //VTK_PARAMETRIC_QUAD_SURFACE
									{
										printf("no order description for VTK_PARAMETRIC_QUAD_SURFACE\n");
										break;
									}
									case 55: //VTK_PARAMETRIC_TETRA_REGION
									{
										printf("no order description for VTK_PARAMETRIC_TETRA_REGION\n");
										break;
									}
									case 56: //VTK_PARAMETRIC_HEX_REGION
									{
										printf("no order description for VTK_PARAMETRIC_HEX_REGION\n");
										break;
									}
									case 60: //VTK_HIGHER_ORDER_EDGE
									{
										printf("no order description for VTK_HIGHER_ORDER_EDGE\n");
										break;
									}
									case 61: //VTK_HIGHER_ORDER_TRIANGLE
									{
										printf("no order description for VTK_HIGHER_ORDER_TRIANGLE\n");
										break;
									}
									case 62: //VTK_HIGHER_ORDER_QUAD
									{
										printf("no order description for VTK_HIGHER_ORDER_QUAD\n");
										break;
									}
									case 63: //VTK_HIGHER_ORDER_POLYGON
									{
										printf("no order description for VTK_HIGHER_ORDER_POLYGON\n");
										break;
									}
									case 64: //VTK_HIGHER_ORDER_TETRAHEDRON
									{
										printf("no order description for VTK_HIGHER_ORDER_TETRAHEDRON\n");
										break;
									}
									case 65: //VTK_HIGHER_ORDER_WEDGE
									{
										printf("no order description for VTK_HIGHER_ORDER_WEDGE\n");
										break;
									}
									case 66: //VTK_HIGHER_ORDER_PYRAMID
									{
										printf("no order description for VTK_HIGHER_ORDER_PYRAMID\n");
										break;
									}
									case 67: //VTK_HIGHER_ORDER_HEXAHEDRON
									{
										printf("no order description for VTK_HIGHER_ORDER_HEXAHEDRON\n");
										break;
									}
									default: //VTK_NUMBER_OF_CELL_TYPES
									{
										printf("Cell type %d is not known\n",ct[i]);
									}
								}
							}
							ReorderEmpty(FACE | EDGE);
							state = R_ATTRIBUTES;
						}
						break;
					}
					case R_SPOINTS:
					{
						printf("STRUCTURED POINTS is not yet implemented\n");
						break;
					}
					case R_POLYDATA:
					{
						printf("POLYDATA is not yet implemented\n");
						break;
					}
					case R_SGRID:
					{
						printf("STRUCTURED GRID is not yet implemented\n");
						break;
					}
					case R_RGRID:
					{
						printf("RECTLINEAR GRID is not yet implemented\n");
						break;
					}
					case R_FIELD:
					{
						printf("FIELD is not yet implemented\n");
						break;
					}
				}
			}
			{
				int count_unused = 0;
				for(unsigned q = 0; q < newnodes.size(); q++)
					if( newnodes[q]->GetMarker(unused_marker) )
					{
						count_unused++;
						delete newnodes[q];
						//~ newnodes[q]->RemMarker(unused_marker);
					}
				if( count_unused ) std::cout << __FILE__ << ":" << __LINE__ << " Warning: deleted " << count_unused << " unused nodes, file " << File << std::endl;
			}
			ReleaseMarker(unused_marker);
			fclose(f);
		} 
		else if(LFile.find(".pmf") != std::string::npos) //this is inner parallel/platform mesh format
		{
			REPORT_STR("start load pmf");
			std::vector<INMOST_DATA_ENUM_TYPE> myprocs;
			std::stringstream in(std::ios::in | std::ios::out | std::ios::binary);
			HeaderType token;
#if defined(USE_MPI)
			if( m_state == Mesh::Parallel )
			{
#if defined(USE_MPI2)
				if( parallel_file_strategy == 1 )
				{
					int ierr;
					std::string buffer;
					MPI_File fh;
					MPI_Status stat;
					ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(File.c_str()),MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
					INMOST_DATA_ENUM_TYPE numprocs = GetProcessorsNumber(), recvsize = 0, mpirank = GetProcessorRank();
					std::vector<INMOST_DATA_ENUM_TYPE> recvsizes;

					if( mpirank == 0 ) //the alternative is to read alltogether
					{
						INMOST_DATA_ENUM_TYPE datanum, chunk,pos,k,q;
						std::vector<INMOST_DATA_ENUM_TYPE> datasizes;
						std::stringstream header;
						
						buffer.resize(3);
						//ierr = MPI_File_read_all(fh,&buffer[0],3,MPI_CHAR,&stat);
						ierr = MPI_File_read_shared(fh,&buffer[0],3,MPI_CHAR,&stat);
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);

						if( static_cast<HeaderType>(buffer[0]) != INMOST::INMOSTFile ) throw BadFile;

						header << buffer.c_str()+1;
						uconv.read_iByteOrder(header);
						uconv.read_iByteSize(header);

						buffer.resize(uconv.get_source_iByteSize());
						//ierr = MPI_File_read_all(fh,&buffer[0],buffer.size(),MPI_CHAR,&stat);
						ierr = MPI_File_read_shared(fh,&buffer[0],buffer.size(),MPI_CHAR,&stat);
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);

						header << buffer;
						uconv.read_iValue(header,datanum);

						buffer.resize(datanum*uconv.get_source_iByteSize());
						datasizes.resize(datanum);
						//ierr = MPI_File_read_all(fh,&buffer[0],buffer.size(),MPI_CHAR,&stat);
						ierr = MPI_File_read_shared(fh,&buffer[0],buffer.size(),MPI_CHAR,&stat);
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);

						header << buffer;
						for(k = 0; k < datanum; k++) uconv.read_iValue(header,datasizes[k]);


						// use this commented code when all processors read the file alltogether through MPI_File_read_all
						//{
						//	MPI_Offset off;
						//	ierr = MPI_File_get_position(fh,&off);
						//	if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
						//	ierr = MPI_File_seek_shared( fh, off, MPI_SEEK_SET );
						//	if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
						//}
						//if( datanum <= numprocs )
						//{
						//	if( mpirank < datanum )
						//		recvsize = datasizes[mpirank];
						//}
						//else
						//{
						//	chunk = static_cast<INMOST_DATA_ENUM_TYPE>(floor(static_cast<double>(datanum)/static_cast<double>(numprocs)));
						//	if( mpirank < numprocs - 1)
						//	{
						//		for(k = chunk*(mpirank); k < chunk*(mpirank+1); k++)
						//			recvsize += datasizes[k];
						//	}
						//	else
						//	{
						//		for(k = chunk*(mpirank); k < datanum; k++)
						//			recvsize += datasizes[k];
						//	}
						//}

						recvsizes.resize(numprocs,0);
						if( datanum <= recvsizes.size() )
						{
							for(k = 0; k < datanum; k++)
								recvsizes[k] = datasizes[k];
						}
						else
						{
							chunk = static_cast<INMOST_DATA_ENUM_TYPE>(floor(static_cast<double>(datanum)/static_cast<double>(recvsizes.size())));
							pos = 0;
							for(k = 0; k < recvsizes.size()-1; k++)
								for(q = 0; q < chunk; q++)
									recvsizes[k] += datasizes[pos++];
							for(k = pos; k < datanum; k++)
								recvsizes[recvsizes.size()-1] += datasizes[k];
						}
					}
					else recvsizes.resize(1); //protect from dereferencing null

					ierr = MPI_Scatter(&recvsizes[0],1,INMOST_MPI_DATA_ENUM_TYPE,&recvsize,1,INMOST_MPI_DATA_ENUM_TYPE,0,GetCommunicator());
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);


					buffer.resize(std::max(1u,recvsize)); //protect from dereferencing null


					{
						ierr = MPI_File_read_ordered(fh,&buffer[0],recvsize,MPI_CHAR,&stat);
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
						in.write(&buffer[0],recvsize);
					}

					ierr = MPI_File_close(&fh);
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
				}
				else
#endif
				{
					int ierr;
					std::string buffer, local_buffer;
					INMOST_DATA_ENUM_TYPE recvsize;
					
					INMOST_DATA_ENUM_TYPE numprocs = GetProcessorsNumber();
					std::vector<INMOST_DATA_ENUM_TYPE> recvsizes(numprocs,0);
					std::vector<int> sendcnts(numprocs), displs(numprocs);
					if( GetProcessorRank() == 0 ) //zero reads everything
					{
						std::fstream fin(File.c_str(),std::ios::in | std::ios::binary);
						fin.get(token);
						if( token != INMOST::INMOSTFile ) throw BadFile;
						uconv.read_iByteOrder(fin);
						uconv.read_iByteSize(fin);
						INMOST_DATA_ENUM_TYPE datanum,k,q,datasum = 0,chunk,pos;
						uconv.read_iValue(fin,datanum);
						std::vector<INMOST_DATA_ENUM_TYPE> datasizes(datanum);
						for(k = 0; k < datanum; k++) 
						{
							uconv.read_iValue(fin,datasizes[k]);
							datasum += datasizes[k];
						}
						{
							buffer.resize(datasum);
							fin.read(&buffer[0],buffer.size());
						}
						fin.close();


						if( datanum <= recvsizes.size() )
						{
							for(k = 0; k < datanum; k++)
								recvsizes[k] = datasizes[k];
						}
						else
						{
							chunk = static_cast<INMOST_DATA_ENUM_TYPE>(floor(static_cast<double>(datanum)/static_cast<double>(recvsizes.size())));
							pos = 0;
							for(k = 0; k < recvsizes.size()-1; k++)
								for(q = 0; q < chunk; q++)
									recvsizes[k] += datasizes[pos++];
							for(k = pos; k < datanum; k++)
								recvsizes[recvsizes.size()-1] += datasizes[k];
						}

						displs[0] = 0;
						sendcnts[0] = recvsizes[0];
						for(k = 1; k < numprocs; k++)
						{
							sendcnts[k] = recvsizes[k];
							displs[k] = sendcnts[k-1]+displs[k-1];
						}
					}
					else 
					{
						//protect from dereferencing null
						buffer.resize(1);
					}
					ierr = MPI_Scatter(&recvsizes[0],1,INMOST_MPI_DATA_ENUM_TYPE,&recvsize,1,INMOST_MPI_DATA_ENUM_TYPE,0,GetCommunicator());
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
					local_buffer.resize(recvsize);
					ierr = MPI_Scatterv(&buffer[0],&sendcnts[0],&displs[0],MPI_CHAR,&local_buffer[0],recvsize,MPI_CHAR,0,GetCommunicator());
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
					in.write(&local_buffer[0],local_buffer.size());
				}

			}
			else
#endif
			{
				std::fstream fin(File.c_str(),std::ios::in | std::ios::binary);
				fin.get(token);
				if( token != INMOST::INMOSTFile ) throw BadFile;
				uconv.read_iByteOrder(fin);
				uconv.read_iByteSize(fin);
				INMOST_DATA_ENUM_TYPE datanum,k,datasum = 0;
				uconv.read_iValue(fin,datanum);
				std::vector<INMOST_DATA_ENUM_TYPE> datasizes(datanum);
				for(k = 0; k < datanum; k++) 
				{
					uconv.read_iValue(fin,datasizes[k]);
					datasum += datasizes[k];
				}
				{
					std::string buffer;
					buffer.resize(datasum);
					fin.read(&buffer[0],buffer.size());
					in.write(&buffer[0],buffer.size());
				}
				fin.close();
			}
			
			//std::fstream in(File.c_str(),std::ios::in | std::ios::binary);
			
			std::vector<Tag> tags;
			std::vector<Node *> old_nodes;
			std::vector<Node *> new_nodes;
			std::vector<Edge *> new_edges;
			std::vector<Face *> new_faces;
			std::vector<Cell *> new_cells;
			std::vector<ElementSet *> new_sets;
			INMOST_DATA_ENUM_TYPE size,i;
			TopologyCheck tmp;
			

			bool start = false;
			
			std::map<GeometricData,ElementType> table;
			
			BeginModification();
		
			while (!(in >> token).eof()) 
			{
				if( !start ) 
				{
					if( token != INMOST::INMOSTFile ) throw BadFile; //check that this is valid file
					else 
					{
						REPORT_STR("File chunk start read");
						//~ std::cout << "start read" << std::endl;
						tags.clear();
						old_nodes.clear();
						new_nodes.clear();
						new_edges.clear();
						new_faces.clear();
						new_cells.clear();
						new_sets.clear();
						old_nodes.reserve(NumberOfNodes());
						for(nodes_container::iterator it = nodes.begin(); it != nodes.end(); it++) if( *it != NULL )
							old_nodes.push_back(*it);
						if( !old_nodes.empty() )qsort(&old_nodes[0],old_nodes.size(),sizeof(Node *),CompareElementsCCentroid);
						if( old_nodes.empty() )
						{
							tmp = GetTopologyCheck(DUPLICATE_CELL | DUPLICATE_FACE | DUPLICATE_EDGE); //we expect not to have duplicates
							RemTopologyCheck(DUPLICATE_CELL | DUPLICATE_FACE | DUPLICATE_EDGE);
						}
						else
						{
							tmp = GetTopologyCheck(DUPLICATE_CELL | DUPLICATE_FACE | DUPLICATE_EDGE); //we expect to have duplicates
							SetTopologyCheck(DUPLICATE_CELL | DUPLICATE_FACE | DUPLICATE_EDGE);
						}
						
						start = true;
					}
				}
				else if (token == INMOST::EoMHeader)  
				{
					
					if( !start ) throw BadFile;
					REPORT_STR("File chunk end read");
					//~ std::cout << "end read" << std::endl;
					if( old_nodes.empty() )
					{
						SetTopologyCheck(tmp);
					}
					else
					{
						if( !(tmp & DUPLICATE_CELL) ) RemTopologyCheck(DUPLICATE_CELL);
						if( !(tmp & DUPLICATE_FACE) ) RemTopologyCheck(DUPLICATE_FACE);
						if( !(tmp & DUPLICATE_EDGE) ) RemTopologyCheck(DUPLICATE_EDGE);
					}
					REPORT_VAL("NODE",new_nodes.size());
					REPORT_VAL("EDGE",new_edges.size());
					REPORT_VAL("FACE",new_faces.size());
					REPORT_VAL("CELL",new_cells.size());
					REPORT_VAL("ESET",new_sets.size());
					REPORT_VAL("TAG",tags.size());
					start = false; //probably the next file is in the input
				}
				else if (token == INMOST::MeshHeader)
				{
					REPORT_STR("MeshHeader");
					uconv.read_iByteOrder(in);
					uconv.read_iByteSize(in);
					iconv.read_iByteOrder(in);
					iconv.read_iByteSize(in);
					iconv.read_fByteOrder(in);
					iconv.read_fByteSize(in);

					INMOST_DATA_ENUM_TYPE header[9],k;
					for(k = 0; k < 9; k++)
						uconv.read_iValue(in,header[k]);
					
					{
						
						char rtemp[5][3];
						in.read(reinterpret_cast<char *>(rtemp),sizeof(rtemp));
						for(GeometricData d = CENTROID; d <= BARYCENTER; d++)
							for(ElementType et = EDGE; et <= CELL; et = et << 1)
								if( rtemp[d][ElementNum(et)-1] ) table[d] |= et;
						
					}

					

					SetDimensions(header[0]);
					new_nodes.reserve(header[1]);
					new_edges.reserve(header[2]);
					new_faces.reserve(header[3]);
					new_cells.reserve(header[4]);
					new_sets.reserve(header[5]);
					tags.reserve(header[6]);
					//~ if( static_cast<Mesh::MeshState>(header[7]) == Mesh::Parallel && m_state != Mesh::Parallel)
						//~ SetCommunicator(INMOST_MPI_COMM_WORLD);
					myprocs.push_back(header[8]);
				}
				else if (token == INMOST::TagsHeader)
				{
					uconv.read_iValue(in,size);
					REPORT_STR("TagsHeader");
					REPORT_VAL("tag_size",size);
					for(i = 0; i < size; i++) tags.push_back(ReadTag(in));
				}
				else if (token == INMOST::NodeHeader)
				{
					uconv.read_iValue(in,size);
					REPORT_STR("NodeHeader");
					REPORT_VAL("node_size",size);
					for(i = 0; i < size; i++) new_nodes.push_back(ReadNode(in,old_nodes));
				}
				else if (token == INMOST::EdgeHeader)
				{
					uconv.read_iValue(in,size);
					REPORT_STR("EdgeHeader");
					REPORT_VAL("edge_size",size);
					for(i = 0; i < size; i++) new_edges.push_back(ReadEdge(in,new_nodes));
				}
				else if (token == INMOST::FaceHeader)
				{
					uconv.read_iValue(in,size);
					REPORT_STR("FaceHeader");
					REPORT_VAL("face_size",size);
					for(i = 0; i < size; i++) new_faces.push_back(ReadFace(in,new_edges));
				}
				else if (token == INMOST::CellHeader)
				{
					uconv.read_iValue(in,size);
					REPORT_STR("CellHeader");
					REPORT_VAL("cell_size",size);
					for(unsigned i = 0; i < size; i++) new_cells.push_back(ReadCell(in,new_faces,new_nodes));
				}
				else if (token == INMOST::ESetHeader)
				{
					uconv.read_iValue(in,size);
					REPORT_STR("EsetHeader");
					REPORT_VAL("eset_size",size);
					for(unsigned i = 0; i < size; i++) new_sets.push_back(ReadElementSet(in,new_nodes,new_edges,new_faces,new_cells));
				}
				else if (token == INMOST::MeshDataHeader)
				{
					REPORT_STR("MeshDataHeader");
					for(size_t j = 0; j < tags.size(); j++) 
					{
						REPORT_VAL("TagName",tags[j].GetTagName());
						Tag * jt = &tags[j];
						for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
						{
							if( jt->isDefined(etype) )
							{
								INMOST_DATA_ENUM_TYPE q, cycle_end;
								Storage * e;
								switch(etype)
								{
								case NODE: cycle_end = new_nodes.size(); break;
								case EDGE: cycle_end = new_edges.size(); break;
								case FACE: cycle_end = new_faces.size(); break;
								case CELL: cycle_end = new_cells.size(); break;
								case ESET: cycle_end = new_sets.size(); break;
								case MESH: cycle_end = 1; break;
								}
								if( jt->isSparse(etype) )
								{
									uconv.read_iValue(in,q);
									//~ std::cout << jt->GetTagName() << std::endl;
									while(q != cycle_end)
									{
										switch(etype)
										{
										case NODE: e = new_nodes[q]; break;
										case EDGE: e = new_edges[q]; break;
										case FACE: e = new_faces[q]; break;
										case CELL: e = new_cells[q]; break;
										case ESET: e = new_sets[q]; break;
										case MESH: e = this; break;
										}
										//~ std::cout << q << "/" << cycle_end << std::endl;
										ReadData(in,*jt,e,new_nodes,new_edges,new_faces,new_cells);
										uconv.read_iValue(in,q);
									}
								}
								else
								{
									for(q = 0; q < cycle_end; q++)
									{
										switch(etype)
										{
										case NODE: e = new_nodes[q]; break;
										case EDGE: e = new_edges[q]; break;
										case FACE: e = new_faces[q]; break;
										case CELL: e = new_cells[q]; break;
										case ESET: e = new_sets[q]; break;
										case MESH: e = this; break;
										}
										ReadData(in,*jt,e,new_nodes,new_edges,new_faces,new_cells);
									}
								}
							}
						}
					}
				}
				else throw BadFile;
			}
			
			
			
			if( m_state == Mesh::Parallel )
			{
#if defined(USE_MPI) 
				
				bool restore_state = false;
				INMOST_DATA_ENUM_TYPE numprocs = GetProcessorsNumber(), size = myprocs.size(),k, procs_sum = 0;
				std::vector<INMOST_DATA_ENUM_TYPE> procs_sizes(numprocs), procs;
				REPORT_MPI(MPI_Allgather(&size,1,INMOST_MPI_DATA_ENUM_TYPE,&procs_sizes[0],1,INMOST_MPI_DATA_ENUM_TYPE,GetCommunicator()));
				for(k = 0; k < numprocs; k++) if( procs_sizes[k] > 1 ) restore_state = true;
				REPORT_VAL("restore_state",restore_state);
				if( restore_state ) //we have to do something with parallel status data
				{
					int ierr;
					std::vector<int> recvcnts(numprocs), displs(numprocs);
					Storage::integer myrank = GetProcessorRank();
					std::sort(myprocs.begin(),myprocs.end());
					//have to allgatherv myprocs from all processors to detect the ownership of elements
					procs_sum = procs_sizes[0];
					recvcnts[0] = procs_sizes[0];
					displs[0] = 0;
					for(k = 1; k < numprocs; k++) 
					{
						procs_sum += procs_sizes[k];
						recvcnts[k] = procs_sizes[k];
						displs[k] = displs[k-1]+recvcnts[k-1];
					}
					myprocs.resize(std::max(1u,static_cast<unsigned>(myprocs.size())));
					procs.resize(procs_sum);
					REPORT_MPI(ierr = MPI_Allgatherv(&myprocs[0],procs_sizes[myrank],INMOST_MPI_DATA_ENUM_TYPE,&procs[0],&recvcnts[0],&displs[0],INMOST_MPI_DATA_ENUM_TYPE,GetCommunicator()));
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
					//we have to distinguish new elements and old elements
					//all new elements with owner in myprocs belong to me
					
					if( procs_sizes[myrank] > 0 )
					{
						for(Mesh::iteratorElement it = BeginElement(CELL | FACE | EDGE | NODE); it != EndElement(); ++it)
							if( it->New() )
							{
								Storage::integer & owner = it->IntegerDF(tag_owner);
								if( std::binary_search(myprocs.begin(),myprocs.end(),static_cast<INMOST_DATA_ENUM_TYPE>(owner)) ) 
									owner = myrank;
								else
								{
									for(k = 0; k < numprocs; k++)
										if( std::binary_search(procs.begin()+displs[k],procs.begin()+displs[k]+recvcnts[k],static_cast<INMOST_DATA_ENUM_TYPE>(owner)) )
										{
											owner = k;
											break;
										}
								}
								Storage::integer_array v = it->IntegerArrayDV(tag_processors);
								for(Storage::integer_array::iterator jt = v.begin(); jt != v.end(); ++jt)
									if( std::binary_search(myprocs.begin(),myprocs.end(),static_cast<INMOST_DATA_ENUM_TYPE>(*jt)) )
										*jt = myrank;
									else
									{
										for(k = 0; k < numprocs; k++)
											if( std::binary_search(procs.begin()+displs[k],procs.begin()+displs[k]+recvcnts[k],static_cast<INMOST_DATA_ENUM_TYPE>(*jt)) )
											{
												*jt = k;
												break;
											}
									}
								std::sort(v.begin(),v.end());
								v.resize(std::unique(v.begin(),v.end())-v.begin());
								if( myrank == owner )
								{
									if( v.size() == 1 )
										it->BulkDF(tag_shared) = Element::Owned;
									else
										it->BulkDF(tag_shared) = Element::Shared;
								}
								else
									it->BulkDF(tag_shared) = Element::Ghost;
							}
					}

				}
				ComputeSharedProcs();
				RecomputeParallelStorage(CELL | FACE | EDGE | NODE);
				
				//Share number of Layers
				REPORT_MPI(MPI_Bcast(&Integer(tag_layers),1,INMOST_MPI_DATA_INTEGER_TYPE,0,GetCommunicator()));
				REPORT_MPI(MPI_Bcast(&Integer(tag_bridge),1,INMOST_MPI_DATA_BULK_TYPE,0,GetCommunicator()));
#else // if there is no mpi, we don't care about statuses
				tag_shared = DeleteTag(tag_shared);
				tag_processors = DeleteTag(tag_processors);
				tag_owner = DeleteTag(tag_owner);
				tag_layers = DeleteTag(tag_layers);
				tag_bridge = DeleteTag(tag_bridge);
				tag_sendto = DeleteTag(tag_sendto);
				m_state = Mesh::Serial;
#endif
			}
			EndModification();
			
			//~ PrepareGeometricData(table);
			RestoreGeometricTags();
			
			if( HaveTag("GLOBAL_ID") )
			{
				tag_global_id = GetTag("GLOBAL_ID");
				for(ElementType etype = NODE; etype <= CELL; etype = etype << 1)
					if( tag_global_id.isDefined(etype) ) have_global_id |= etype;
			}
#if defined(USE_MPI)
			if( m_state == Mesh::Parallel )
			{
				int hgi = have_global_id, recvtype = 0;
				//~ int flag = 0, recvflag;
				REPORT_MPI(MPI_Allreduce(&hgi,&recvtype,1,MPI_INT,MPI_BOR,comm));
				REPORT_VAL("local  types", hgi);
				REPORT_VAL("global types", recvtype);
				for(ElementType etype = NODE; etype <= CELL; etype = etype << 1 )
				{
					REPORT_VAL("test global id tag type",ElementTypeName(etype));
					if( (etype & recvtype) && !(etype & have_global_id) ) 
					{
						tag_global_id = CreateTag("GLOBAL_ID",DATA_INTEGER, etype, NONE,1);
						have_global_id |= etype;
						//~ flag = 1;
						REPORT_VAL("created new global id tag",ElementTypeName(etype));
					}
				}
				//~ REPORT_MPI(MPI_Allreduce(&flag,&recvflag,1,MPI_INT,MPI_BOR,comm));
				//~ REPORT_VAL("flag",&recvflag);
				//~ if( recvflag ) 
				AssignGlobalID(recvtype);
			}
#endif
		}
		else throw NotImplemented;
		EXIT_FUNC();
	}
	
	int VtkElementType(ElementType t)
	{
		switch(t)
		{
			case Element::Tri: return 5;
			case Element::Quad: return 9;
			case Element::MultiLine: return 4;
			case Element::Polygon: return 7;
			case Element::Tet: return 10;
			case Element::Hex: return 12;
			case Element::Prism: return 13;
			case Element::Pyramid: return 14;
			case Element::Polyhedron: return 42;
			case Element::MultiPolygon: return 42;
		}
		return -1;
	}
	
	
	
	void Mesh::Save(std::string File)
	{
		ENTER_FUNC();
		REPORT_VAL("File",File);
		std::string LFile;
		LFile.resize(File.size());
		std::transform(File.begin(),File.end(),LFile.begin(),::tolower);
		if(LFile.find(".pvtk") != std::string::npos) //this is legacy parallel vtk
		{
			std::string name=File;
			int pos=name.rfind(".pvtk");
			name.erase(pos); 

			int l=name.find_last_of("/\\");
			std::string fname=name.substr(l+1,name.length());
			if(GetProcessorRank()==0)
			{//create pvtk file
				std::stringstream ss;
				std::string numproc;
				ss << GetProcessorsNumber();
				numproc=ss.str();
				FILE  *f=fopen(File.c_str(), "w");
				fprintf(f,"%s%s%s", "<File version=\"pvtk-1.0\" dataType=\"vtkUnstructuredGrid\" numberOfPieces=\"", numproc.c_str(),"\">\n");
				for (int i=0; i<GetProcessorsNumber(); i++)
				{
					std::stringstream ss;
					ss << i;//task: leading zeroes
					std::string end="_"+ss.str()+".vtk";

					std::string temp=fname;
					temp.append(end);
					fprintf(f, "%s%s%s", "<Piece fileName=\"" ,temp.c_str() ,"\"/>\n");
				}
				fprintf(f, "%s", "</File>");
				fclose(f);
			}
			
			std::stringstream ss;
			ss << GetProcessorRank();//task: leading zeroes
			std::string end="_"+ss.str()+".vtk";
			
			name.append(end);
			Save(name);
		}
		else if(LFile.find(".vtk") != std::string::npos) //this is legacy vtk
		{
			unsigned int dim = GetDimensions();
			if( dim > 3 )
			{
				printf("VTK file supports 3 dimensions max\n");
				return;
			}
			FILE * f = fopen(File.c_str(),"w");
			if( !f ) throw BadFileName;
			fprintf(f,"# vtk DataFile Version 3.0\n");
			fprintf(f,"file is written by INMOST\n");
			fprintf(f,"ASCII\n");
			fprintf(f,"DATASET UNSTRUCTURED_GRID\n");
			ReorderEmpty(CELL | NODE);
			fprintf(f,"POINTS %u double\n",nodes.size());
			for(Mesh::nodes_container::iterator it = nodes.begin(); it != nodes.end(); it++)
			{
				Storage::real_array coords = (*it)->RealArray(CoordsTag());
				for(unsigned int i = 0; i < dim; i++) 
				{
					double temp = coords[i];
					fprintf(f,"%.10f ",temp);
				}
				for(unsigned int i = dim; i < 3; i++)
					fprintf(f,"0 ");
				fprintf(f,"\n");
			}
			{
				std::vector<int> values;
				for(Mesh::cells_container::iterator it = cells.begin(); it != cells.end(); it++)
				{	
					switch((*it)->GetGeometricType())
					{
						case Element::Tri:
						case Element::Quad:
						case Element::MultiLine:
						case Element::Polygon:
						case Element::Tet:
						{
							adjacent<Node> nodes = (*it)->getNodes();
							values.push_back(nodes.size());
							for(adjacent<Node>::iterator jt = nodes.begin(); jt != nodes.end(); jt++)
								values.push_back(jt->LocalID());
							break;
						}
						case Element::Prism:
						{
							adjacent<Node> nodes = (*it)->getNodes();
							values.push_back(nodes.size());
							values.push_back(nodes[0].LocalID());
							values.push_back(nodes[2].LocalID());
							values.push_back(nodes[1].LocalID());
							values.push_back(nodes[3].LocalID());
							values.push_back(nodes[5].LocalID());
							values.push_back(nodes[4].LocalID());
							break;
						}
						case Element::Hex:
						{
							adjacent<Node> nodes = (*it)->getNodes();
							values.push_back(nodes.size());
							values.push_back(nodes[0].LocalID());
							values.push_back(nodes[3].LocalID());
							values.push_back(nodes[2].LocalID());
							values.push_back(nodes[1].LocalID());
							values.push_back(nodes[4].LocalID());
							values.push_back(nodes[7].LocalID());
							values.push_back(nodes[6].LocalID());
							values.push_back(nodes[5].LocalID());
							break;
						}
						case Element::Pyramid:
						{
							adjacent<Node> nodes = (*it)->getNodes();
							values.push_back(nodes.size());
							values.push_back(nodes[0].LocalID());
							values.push_back(nodes[3].LocalID());
							values.push_back(nodes[2].LocalID());
							values.push_back(nodes[1].LocalID());
							values.push_back(nodes[4].LocalID());
							break;
						}
						case Element::Polyhedron:
						case Element::MultiPolygon:
						{
							//printf("polyhedron!!!\n");
							adjacent<Face> faces = (*it)->getFaces();
                            int totalNum = 1 + faces.size();
                            for(adjacent<Face>::iterator jt = faces.begin(); jt != faces.end(); jt++)
                            {
                                adjacent<Node> nodes = jt->getNodes();
                                totalNum += nodes.size();

                            }
                            values.push_back(totalNum);
							values.push_back(faces.size());
							for(adjacent<Face>::iterator jt = faces.begin(); jt != faces.end(); jt++)
							{
								adjacent<Node> nodes = jt->getNodes();
								values.push_back(nodes.size());
								if( jt->FaceOrientedOutside(*it) )
									for(adjacent<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
										values.push_back(kt->LocalID());
								else
									for(adjacent<Node>::reverse_iterator kt = nodes.rbegin(); kt != nodes.rend(); kt++)
										values.push_back(kt->LocalID());
							}
							break;
						}
						default: printf("This should not happen %s\n",Element::GeometricTypeName((*it)->GetGeometricType()));
					}
				}
				fprintf(f,"CELLS %u %ld\n",cells.size(),values.size());
				for(unsigned int i = 0; i < values.size(); i++)
				{
					fprintf(f,"%d ",values[i]);
					if( (i+1) % 20 == 0) fprintf(f,"\n");
				}
				fprintf(f,"\n");
			}
			fprintf(f,"CELL_TYPES %u\n",cells.size());
			for(Mesh::cells_container::iterator it = cells.begin(); it != cells.end(); it++)
				fprintf(f,"%d\n",VtkElementType((*it)->GetGeometricType()));
			
			{
				std::vector<std::string> tag_names;
				std::vector<Tag> tags;
				ListTagNames(tag_names);
				for(unsigned int i = 0; i < tag_names.size(); i++)
				{
					Tag t = GetTag(tag_names[i]);
					//printf("%s %d %d %d\n",tag_names[i].c_str(),t.isDefined(CELL),!t.isSparse(CELL),t.GetDataType() != DATA_BULK);
					if( t.isDefined(CELL) && !t.isSparse(CELL) && t.GetDataType() != DATA_BULK && t.GetDataType() != DATA_REFERENCE &&
					   t != CoordsTag() && t != SharedTag() && t != SendtoTag() && t != ProcessorsTag())
					{
						//printf("added!\n");
						tags.push_back(t);
					}
				}
				
				if( !tags.empty() ) fprintf(f,"CELL_DATA %u\n",cells.size());
				for(unsigned int i = 0; i < tags.size(); i++)
				{
					unsigned int comps = tags[i].GetSize();
					if( comps == ENUMUNDEF )
					{
						//printf("Warning: vtk don't support arrays of variable size (tag name: %s)\n",tags[i].GetTagName().c_str());
						continue;
					}
					else
					{
						{
							fprintf(f,"SCALARS %s %s %d\n",tags[i].GetTagName().c_str(),(tags[i].GetDataType() == DATA_REAL ? "double" : "int"),comps);
							fprintf(f,"LOOKUP_TABLE default\n");
							for(Mesh::cells_container::iterator it = cells.begin(); it != cells.end(); it++)
							{
								switch( tags[i].GetDataType() )
								{
									case DATA_REAL:
									{
										Storage::real_array arr = (*it)->RealArray(tags[i]);
										for(unsigned int m = 0; m < comps; m++) fprintf(f,"%f ",arr[m]);
										fprintf(f,"\n");
									}
									break;
									case DATA_INTEGER:
									{
										Storage::integer_array arr = (*it)->IntegerArray(tags[i]);
										for(unsigned int m = 0; m < comps; m++) fprintf(f,"%d ",arr[m]);
										fprintf(f,"\n");
									}
									break;
									default: continue;
								}
							}
						}
					}
				}
			}
			{
				std::vector<std::string> tag_names;
				std::vector<Tag> tags;
				ListTagNames(tag_names);
				for(unsigned int i = 0; i < tag_names.size(); i++)
				{
					Tag t = GetTag(tag_names[i]);
					if( t.isDefined(NODE) && !t.isSparse(NODE) && t.GetDataType() != DATA_BULK && t.GetDataType() != DATA_REFERENCE &&
					   t != CoordsTag() && t != SharedTag() && t != SendtoTag() && t != ProcessorsTag())
						tags.push_back(t);
				}
				
				if( !tags.empty() ) fprintf(f,"POINT_DATA %u\n",nodes.size());
				for(unsigned int i = 0; i < tags.size(); i++)
				{
					unsigned int comps = tags[i].GetSize();
					if( comps == ENUMUNDEF )
					{
						//printf("Warning: vtk don't support arrays of variable size (tag name: %s)\n",tags[i].GetTagName().c_str());
						continue;
					}
					else
					{
						{
							fprintf(f,"SCALARS %s %s %d\n",tags[i].GetTagName().c_str(),(tags[i].GetDataType() == DATA_REAL ? "double" : "int"),comps);
							fprintf(f,"LOOKUP_TABLE default\n");
							for(Mesh::nodes_container::iterator it = nodes.begin(); it != nodes.end(); it++)
							{
								switch( tags[i].GetDataType() )
								{
									case DATA_REAL:
									{
										Storage::real_array arr = (*it)->RealArray(tags[i]);
										for(unsigned int m = 0; m < comps; m++) fprintf(f,"%f ",arr[m]);
										fprintf(f,"\n");
									}
									break;
									case DATA_INTEGER:
									{
										Storage::integer_array arr = (*it)->IntegerArray(tags[i]);
										for(unsigned int m = 0; m < comps; m++) fprintf(f,"%d ",arr[m]);
										fprintf(f,"\n");
									}
									break;
									default: continue;
								}
							}
						}
					}
				}
			}
			fclose(f);
			
		}
		else if( LFile.find(".gmv") != std::string::npos) //this is gmv file
		{
			Storage::integer keynum;
			Storage::real keyval;
			ReorderEmpty(CELL | FACE | NODE | ESET);
			FILE * file = fopen(File.c_str(),"wb");
			char keyword[2048];
			sprintf(keyword,"gmvinput"); fwrite(keyword,1,8,file);
			sprintf(keyword,"ieee");
			if( sizeof(Storage::real) != 8 || sizeof(Storage::integer) != 8 )
				sprintf(keyword,"ieeei%ldr%ld",sizeof(Storage::integer),sizeof(Storage::real));
			fwrite(keyword,1,8,file);
			sprintf(keyword,"nodev"); fwrite(keyword,1,8,file);
			keynum = static_cast<Storage::integer>(nodes.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
			for(Mesh::nodes_container::iterator n = nodes.begin(); n != nodes.end(); n++)
				fwrite(&(*n)->Coords()[0],sizeof(Storage::real),3,file);
			sprintf(keyword,"faces"); fwrite(keyword,1,8,file);
			keynum = static_cast<Storage::integer>(faces.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
			keynum = static_cast<Storage::integer>(cells.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
			for(Mesh::faces_container::iterator f = faces.begin(); f != faces.end(); f++)
			{
				adjacent<Node> fnodes = (*f)->getNodes();
				//sprintf(keyword,"nverts"); fwrite(keyword,1,8,file);
				keynum = static_cast<Storage::integer>(fnodes.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
				//sprintf(keyword,"vertex_ids"); fwrite(keyword,1,8,file);
				for(adjacent<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
				{
					keynum = fn->LocalID()+1;
					fwrite(&keynum,sizeof(Storage::integer),1,file);
				}
				adjacent<Cell> fcells = (*f)->getCells();
				//sprintf(keyword,"cellno1"); fwrite(keyword,1,8,file);
				keynum = fcells.size() > 0 ? fcells[0].LocalID()+1 : 0;
				fwrite(&keynum,sizeof(Storage::integer),1,file);
				//sprintf(keyword,"cellno2"); fwrite(keyword,1,8,file);
				keynum = fcells.size() > 1 ? fcells[1].LocalID()+1 : 0;
				fwrite(&keynum,sizeof(Storage::integer),1,file);
			}
			
			
			if( HaveTag("MATERIAL") )
			{
				Tag t = GetTag("MATERIAL");
				if( t.isDefined(CELL) )
				{
					std::set<Storage::integer> mat;
					for(Mesh::cells_container::iterator c = cells.begin(); c != cells.end(); c++)
						mat.insert((*c)->Integer(t)+1);
					sprintf(keyword,"material"); fwrite(keyword,1,8,file);
					keynum = static_cast<Storage::integer>(mat.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
					keynum = 0; fwrite(&keynum,sizeof(Storage::integer),1,file);
					for(std::set<Storage::integer>::iterator it = mat.begin(); it != mat.end(); it++)
					{
						sprintf(keyword,"mat%d",*it); fwrite(keyword,1,8,file);
					}
					for(Mesh::cells_container::iterator c = cells.begin(); c != cells.end(); c++)
					{
						keynum = (*c)->Integer(t)+1; fwrite(&keynum,sizeof(Storage::integer),1,file);
					}
				}
				
			}
			
			sprintf(keyword,"variable"); fwrite(keyword,1,8,file);
			for(ElementType etype = CELL; etype >= NODE; etype = etype >> 1)
			{
				if( etype == EDGE ) continue;
				for(Mesh::iteratorTag t = BeginTag(); t != EndTag(); t++) 
					if( t->isDefined(etype) && t->GetSize() == 1 && !t->isSparse(etype) && t->GetTagName().substr(0,9) != "PROTECTED" && t->GetDataType() != DATA_REFERENCE )
					{
						sprintf(keyword,"%s",t->GetTagName().substr(0,8).c_str());
						fwrite(keyword,1,8,file);
						switch(etype)
						{
							case NODE: keynum = 1; break;
							case FACE: keynum = 2; break;
							case CELL: keynum = 0; break;
							default: throw NotImplemented;
						}
						fwrite(&keynum,sizeof(Storage::integer),1,file);
						for(Mesh::iteratorElement e = BeginElement(etype); e != EndElement(); e++)
						{
							keyval = 1e+20;
							switch(t->GetDataType())
							{
								case DATA_INTEGER: keyval = static_cast<Storage::real>(e->Integer(*t)); break;
								case DATA_REAL: keyval = e->Real(*t); break;
								case DATA_BULK: keyval = static_cast<Storage::real>(e->Bulk(*t)); break;
								default: throw NotImplemented;
							}
							fwrite(&keyval,sizeof(Storage::real),1,file);
						}
					}
			}
			sprintf(keyword,"endvars"); fwrite(keyword,1,8,file);
			
			sprintf(keyword,"subvars"); fwrite(keyword,1,8,file);
			for(ElementType etype = CELL; etype >= NODE; etype = etype >> 1)
			{
				if( etype == EDGE ) continue;
				for(Mesh::iteratorTag t = BeginTag(); t != EndTag(); t++) 
					if( t->isDefined(etype) && t->GetSize() == 1 && t->isSparse(etype) && t->GetDataType() != DATA_REFERENCE && t->GetDataType() != DATA_REFERENCE)
					{
						Storage::integer temp;
						keynum = 0;
						for(Mesh::iteratorElement e = BeginElement(etype); e != EndElement(); e++)
							if( e->HaveData(*t) ) 
								keynum++;
						if( keynum )
						{
							sprintf(keyword,"%s",t->GetTagName().substr(0,8).c_str());
							fwrite(keyword,1,8,file);
							switch(etype)
							{
								case NODE: temp = 1; break;
								case FACE: temp = 2; break;
								case CELL: temp = 0; break;
								default: throw NotImplemented;
							}
							fwrite(&temp,sizeof(Storage::integer),1,file);
							fwrite(&keynum,sizeof(Storage::integer),1,file);
							for(Mesh::iteratorElement e = BeginElement(etype); e != EndElement(); e++)
								if( e->HaveData(*t) ) 
								{
									keynum = e->LocalID()+1;
									fwrite(&keynum,sizeof(Storage::integer),1,file);
								}
							for(Mesh::iteratorElement e = BeginElement(etype); e != EndElement(); e++)
								if( e->HaveData(*t) )
								{
									switch(t->GetDataType())
									{
										case DATA_INTEGER: keyval = static_cast<Storage::real>(e->Integer(*t)); break;
										case DATA_REAL: keyval = e->Real(*t); break;
										case DATA_BULK: keyval = static_cast<Storage::real>(e->Bulk(*t)); break;
										default: throw NotImplemented;
									}
									fwrite(&keyval,sizeof(Storage::real),1,file);
								}
						}
					}
			}
			sprintf(keyword,"endsubv"); fwrite(keyword,1,8,file);
			
			sprintf(keyword,"groups"); fwrite(keyword,1,8,file);
			for(Mesh::sets_container::iterator set = sets.begin(); set != sets.end(); set++)
			{

				for(ElementType etype = CELL; etype >= NODE; etype = etype >> 1)
				{
					if( etype == EDGE ) continue;
					keynum = 0;
					for(ElementSet::iterator it = (*set)->begin(); it != (*set)->end(); it++)
						if( it->GetElementType() == etype ) keynum++;
					if( keynum )
					{
						Storage::integer temp;
						switch(etype)
						{
							case NODE: temp = 1; break;
							case FACE: temp = 2; break;
							case CELL: temp = 0; break;
							default: throw NotImplemented;
						}
						sprintf(keyword,"set%d_%s\n",(*set)->LocalID()+1,ElementTypeName(etype));
						fwrite(keyword,1,8,file);
						fwrite(&temp,sizeof(Storage::integer),1,file);						
						fwrite(&keynum,sizeof(Storage::integer),1,file);
						for(ElementSet::iterator it = (*set)->begin(); it != (*set)->end(); it++)
							if( it->GetElementType() == etype )
							{
								keynum = it->LocalID()+1;
								fwrite(&keynum,sizeof(Storage::integer),1,file);
							} 
					}	
				}
			}
			sprintf(keyword,"endgrp"); fwrite(keyword,1,8,file);
			
			sprintf(keyword,"vectors"); fwrite(keyword,1,8,file);
			for(ElementType etype = CELL; etype >= NODE; etype = etype >> 1)
			{
				if( etype == EDGE ) continue;
				for(Mesh::iteratorTag t = BeginTag(); t != EndTag(); t++) 
					if( t->isDefined(etype) && t->GetSize() != 1 && t->GetSize() != ENUMUNDEF && !t->isSparse(etype) && t->GetDataType() != DATA_REFERENCE)
					{
						sprintf(keyword,"%s",t->GetTagName().substr(0,8).c_str());
						fwrite(keyword,1,8,file);
						switch(etype)
						{
							case NODE: keynum = 1; break;
							case FACE: keynum = 2; break;
							case CELL: keynum = 0; break;
							default: throw NotImplemented;
						}
						fwrite(&keynum,sizeof(Storage::integer),1,file);
						keynum = t->GetSize();
						fwrite(&keynum,sizeof(Storage::integer),1,file);
						keynum = 0;
						fwrite(&keynum,sizeof(Storage::integer),1,file);
						
						for(size_t q = 0; q < t->GetSize(); q++)
						{
							for(Mesh::iteratorElement e = BeginElement(etype); e != EndElement(); e++)
							{
								switch(t->GetDataType())
								{
									case DATA_INTEGER: keyval = static_cast<Storage::real>(e->IntegerArray(*t)[q]); break;
									case DATA_REAL: keyval = e->RealArray(*t)[q]; break;
									case DATA_BULK: keyval = static_cast<Storage::real>(e->BulkArray(*t)[q]); break;
									default: throw NotImplemented;
								}
								fwrite(&keyval,sizeof(Storage::real),1,file);
							}
						}
					}
			}
			sprintf(keyword,"endvect"); fwrite(keyword,1,8,file);
			
			sprintf(keyword,"polygons"); fwrite(keyword,1,8,file);
			for(Mesh::faces_container::iterator f = faces.begin(); f != faces.end(); f++) if( (*f)->Boundary() )
			{
				adjacent<Node> fnodes = (*f)->getNodes();
				keynum = 1; fwrite(&keynum,sizeof(Storage::integer),1,file);
				keynum = static_cast<Storage::integer>(fnodes.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
				for(adjacent<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
					fwrite(&fn->Coords()[0],sizeof(Storage::real),1,file);
				for(adjacent<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
					fwrite(&fn->Coords()[1],sizeof(Storage::real),1,file);
				for(adjacent<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
					fwrite(&fn->Coords()[2],sizeof(Storage::real),1,file);
			}
			sprintf(keyword,"endpoly"); fwrite(keyword,1,8,file);
			sprintf(keyword,"endgmv"); fwrite(keyword,1,8,file);
			fclose(file);
		}
		else if(LFile.find(".pmf") != std::string::npos) //this is inner parallel/platform mesh format
		{
			//~ if( m_state == Mesh::Serial ) SetCommunicator(INMOST_MPI_COMM_WORLD);
			std::stringstream out(std::ios::in | std::ios::out | std::ios::binary);
			
			out << INMOST::INMOSTFile;
			out << INMOST::MeshHeader;
			ReorderEmpty(NODE | EDGE | FACE | CELL | ESET);

			uconv.write_iByteOrder(out);
			uconv.write_iByteSize(out);
			iconv.write_iByteOrder(out);
			iconv.write_iByteSize(out);
			iconv.write_fByteOrder(out);
			iconv.write_fByteSize(out);

			INMOST_DATA_ENUM_TYPE header[9] = {GetDimensions(), NumberOfNodes(), NumberOfEdges(), NumberOfFaces(), NumberOfCells(), NumberOfSets(), NumberOfTags(), m_state, GetProcessorRank()},k;
			for(k = 0; k < 9; k++) uconv.write_iValue(out,header[k]);
			out.write(reinterpret_cast<char *>(remember),sizeof(remember));
		
			// Tags
			out << INMOST::TagsHeader;
			uconv.write_iValue(out,header[6]);
			REPORT_STR("TagsHeader");
			REPORT_VAL("tag_size",header[6]);
			for(Mesh::iteratorTag it = BeginTag(); it != EndTag(); it++)  WriteTag(out,*it);

			// Nodes 
			out << INMOST::NodeHeader;
			uconv.write_iValue(out,header[1]);
			REPORT_STR("NodeHeader");
			REPORT_VAL("node_size",header[1]);
			for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); it++) WriteNode(out,&*it);
		
			// Edges
			out << INMOST::EdgeHeader;
			uconv.write_iValue(out,header[2]);
			REPORT_STR("EdgeHeader");
			REPORT_VAL("edge_size",header[2]);
			for(Mesh::iteratorEdge it = BeginEdge(); it != EndEdge(); it++) WriteEdge(out,&*it);
		
			// Faces
			out << INMOST::FaceHeader;
			uconv.write_iValue(out,header[3]);
			REPORT_STR("FaceHeader");
			REPORT_VAL("face_size",header[3]);
			for(Mesh::iteratorFace it = BeginFace(); it != EndFace(); it++) WriteFace(out,&*it);
		
			// Cells
			out << INMOST::CellHeader;
			uconv.write_iValue(out,header[4]);
			REPORT_STR("CellHeader");
			REPORT_VAL("cell_size",header[4]);
			for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); it++)  WriteCell(out,&*it);
		
			// Element Sets
			out << INMOST::ESetHeader;
			REPORT_STR("ESetHeader");
			REPORT_VAL("eset_size",header[5]);
			uconv.write_iValue(out,header[5]);
			for(Mesh::iteratorSet it = BeginSet(); it != EndSet(); it++) WriteElementSet(out,&*it);
		
			out << INMOST::MeshDataHeader;
			REPORT_STR("MeshDataHeader");

			for(Mesh::iteratorTag jt = BeginTag(); jt != EndTag(); jt++)
			{
				REPORT_VAL("TagName",jt->GetTagName());
				for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
					if( jt->isDefined(etype) ) 
					{
						if( jt->isSparse(etype) )
						{
							INMOST_DATA_ENUM_TYPE q = 0;
							//~ std::cout << jt->GetTagName() << std::endl;
							for(Mesh::base_iterator it = Begin(etype); it != End(); it++) 
							{
								if( it->HaveData(*jt) )
								{
									uconv.write_iValue(out,q);
									//~ std::cout << q << std::endl;
									WriteData(out,*jt,*it);
								}
								q++;
							}
							//~ std::cout << q << std::endl;
							uconv.write_iValue(out,q);
						}
						else
						{
							for(Mesh::base_iterator it = Begin(etype); it != End(); it++)  
								WriteData(out,*jt,*it);
						}
					}
			}
		
		
			out << INMOST::EoMHeader;
#if defined(USE_MPI)
			if( m_state == Mesh::Parallel )
			{
				REPORT_STR("Parallel write");
				int ierr;
				INMOST_DATA_ENUM_TYPE numprocs = GetProcessorsNumber(), datasize = static_cast<INMOST_DATA_ENUM_TYPE>(out.tellp()),k;
				std::vector<INMOST_DATA_ENUM_TYPE> datasizes(numprocs,0);
				REPORT_VAL("local_write_file_size",datasize);
				REPORT_MPI(ierr = MPI_Gather(&datasize,1,INMOST_MPI_DATA_ENUM_TYPE,&datasizes[0],1,INMOST_MPI_DATA_ENUM_TYPE,0,GetCommunicator()));
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
#if defined(USE_MPI2) //We have access to MPI_File
				if( parallel_file_strategy == 1 )
				{
					MPI_File fh;
					MPI_Status stat;
					ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(File.c_str()),MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
					if( GetProcessorRank() == 0 )
					{
						std::stringstream header;
						header.put(INMOST::INMOSTFile);
						uconv.write_iByteOrder(header);
						uconv.write_iByteSize(header);
						uconv.write_iValue(header,numprocs);
						for(k = 0; k < numprocs; k++) uconv.write_iValue(header,datasizes[k]);

						std::string header_data(header.str());
						ierr = MPI_File_write_shared(fh,&header_data[0],header_data.size(),MPI_CHAR,&stat);
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
					}
					{
						std::string local_data(out.str());
						ierr = MPI_File_write_ordered(fh,&local_data[0],local_data.size(),MPI_CHAR,&stat);
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
					}

					ierr = MPI_File_close(&fh);
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
				}
				else
#endif
				{
					std::vector<int> displs(numprocs),recvcounts(numprocs);
					std::string file_contents;
					std::string local_data(out.str());
					if( GetProcessorRank() == 0 )
					{
						recvcounts[0] = datasizes[0];
						displs[0] = 0;
						int sizesum = recvcounts[0];
						for(k = 1; k < numprocs; k++)
						{
							recvcounts[k] = datasizes[k];
							displs[k] = displs[k-1]+recvcounts[k-1];
							sizesum += recvcounts[k];
						}
						file_contents.resize(sizesum);
					}
					else file_contents.resize(1); //protect from accessing bad pointer
					ierr = MPI_Gatherv(&local_data[0],local_data.size(),MPI_CHAR,&file_contents[0],&recvcounts[0],&displs[0],MPI_CHAR,0,GetCommunicator());
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
					if( GetProcessorRank() == 0 )
					{
						std::fstream fout(File.c_str(),std::ios::out | std::ios::binary);
						fout.put(INMOST::INMOSTFile);
						uconv.write_iByteOrder(fout);
						uconv.write_iByteSize(fout);
						uconv.write_iValue(fout,numprocs);
						for(k = 0; k < numprocs; k++) uconv.write_iValue(fout,datasizes[k]);
						fout.write(&file_contents[0],file_contents.size());
						fout.close();
					}
				}
			}
			else
#endif
			{
				REPORT_STR("Serial write");
				std::fstream fout(File.c_str(),std::ios::out | std::ios::binary);
				INMOST_DATA_ENUM_TYPE numprocs = 1, datasize = static_cast<INMOST_DATA_ENUM_TYPE>(out.tellp());
				REPORT_VAL("write_file_size",datasize);
				fout.put(INMOST::INMOSTFile);
				uconv.write_iByteOrder(fout);
				uconv.write_iByteSize(fout);
				uconv.write_iValue(fout,numprocs);
				uconv.write_iValue(fout,datasize);
				fout << out.rdbuf();
				fout.close();
			}
			
		}
		else throw NotImplemented;	
		EXIT_FUNC();
	}
	
	bool Mesh::isParallelFileFormat(std::string File)
	{
		std::string LFile;
		LFile.resize(File.size());
		std::transform(File.begin(),File.end(),LFile.begin(),::tolower);
		if(LFile.find(".grdecl") != std::string::npos) return false;
		else if(LFile.find(".vtk") != std::string::npos) return false;
		else if(LFile.find(".gmv") != std::string::npos) return false;
		else if(LFile.find(".pvtk") != std::string::npos) return true;
		else if(LFile.find(".pmf") != std::string::npos) return true;
		throw NotImplemented;
	}
}
#endif
