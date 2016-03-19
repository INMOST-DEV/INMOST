#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif


#include "inmost.h"

#if defined(USE_MESH)

//vtk states

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
	else 
	{
		std::cout << __FILE__ << ":" << __LINE__ << " Cannot read 3 coordinates from file " << std::endl;
		throw INMOST::BadFile;
	}
}

namespace INMOST
{
  
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
		assert(false);
		return -1;
	}
	
	INMOST_DATA_ENUM_TYPE VtkElementNodes(ElementType t)
	{
		switch(t)
		{
			case Element::Tri: return 3;
			case Element::Quad: return 4;
			case Element::MultiLine: return ENUMUNDEF;
			case Element::Polygon: return ENUMUNDEF;
			case Element::Tet: return 4;
			case Element::Hex: return 8;
			case Element::Prism: return 6;
			case Element::Pyramid: return 5;
			case Element::Polyhedron: return ENUMUNDEF;
			case Element::MultiPolygon: return ENUMUNDEF;
		}
		assert(false);
		return ENUMUNDEF;
	}
	

  void Mesh::SaveVTK(std::string File)
  {
		integer dim = GetDimensions();
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
		//ReorderEmpty(CELL | NODE);
		Tag set_id = CreateTag("TEMPORARY_ELEMENT_ID",DATA_INTEGER,CELL | NODE,NONE,1);
		Storage::integer cur_num = 0;
		for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it) it->IntegerDF(set_id) = cur_num++;
		cur_num = 0;
		for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); ++it) it->IntegerDF(set_id) = cur_num++;
		fprintf(f,"POINTS %u double\n",NumberOfNodes());
		for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); it++)
		{
			Storage::real_array coords = it->RealArray(CoordsTag());
			for(integer i = 0; i < dim; i++) 
			{
				double temp = coords[i];
				fprintf(f,"%.10f ",temp);
			}
			for(integer i = dim; i < 3; i++)
				fprintf(f,"0 ");
			fprintf(f,"\n");
		}
		{
			dynarray<int,64> values;
			for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); it++)
			{	
				switch(it->GetGeometricType())
				{
					case Element::Tri:
					case Element::Quad:
					case Element::MultiLine:
					case Element::Polygon:
					case Element::Tet:
					{
						ElementArray<Node> nodes = it->getNodes();
						values.push_back(static_cast<integer>(nodes.size()));
						for(ElementArray<Node>::iterator jt = nodes.begin(); jt != nodes.end(); jt++)
							values.push_back(jt->IntegerDF(set_id));
						break;
					}
					case Element::Prism:
					{
						ElementArray<Node> nodes = it->getNodes();
						if( nodes.size() != 6 ) goto safe_output;
						values.push_back(static_cast<integer>(nodes.size()));
						values.push_back(nodes[0].IntegerDF(set_id));
						values.push_back(nodes[2].IntegerDF(set_id));
						values.push_back(nodes[1].IntegerDF(set_id));
						values.push_back(nodes[3].IntegerDF(set_id));
						values.push_back(nodes[5].IntegerDF(set_id));
						values.push_back(nodes[4].IntegerDF(set_id));
						break;
					}
					case Element::Hex:
					{
						ElementArray<Node> nodes = it->getNodes();
						if( nodes.size() != 8 ) goto safe_output;
						values.push_back(static_cast<integer>(nodes.size()));
						values.push_back(nodes[0].IntegerDF(set_id));
						values.push_back(nodes[3].IntegerDF(set_id));
						values.push_back(nodes[2].IntegerDF(set_id));
						values.push_back(nodes[1].IntegerDF(set_id));
						values.push_back(nodes[4].IntegerDF(set_id));
						values.push_back(nodes[7].IntegerDF(set_id));
						values.push_back(nodes[6].IntegerDF(set_id));
						values.push_back(nodes[5].IntegerDF(set_id));
						break;
					}
					case Element::Pyramid:
					{
						ElementArray<Node> nodes = it->getNodes();
						if( nodes.size() != 5 ) goto safe_output;
						values.push_back(static_cast<integer>(nodes.size()));
						values.push_back(nodes[0].IntegerDF(set_id));
						values.push_back(nodes[3].IntegerDF(set_id));
						values.push_back(nodes[2].IntegerDF(set_id));
						values.push_back(nodes[1].IntegerDF(set_id));
						values.push_back(nodes[4].IntegerDF(set_id));
						break;
					}
					case Element::Polyhedron:
					case Element::MultiPolygon:
					{
safe_output:
						//printf("polyhedron!!!\n");
						ElementArray<Face> faces = it->getFaces();
                          integer totalNum = 1 + static_cast<integer>(faces.size());
                          for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); jt++)
                              totalNum += jt->nbAdjElements(NODE);

                          values.push_back(totalNum);
						values.push_back(static_cast<integer>(faces.size()));
						for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); jt++)
						{
							ElementArray<Node> nodes = jt->getNodes();
							values.push_back(static_cast<integer>(nodes.size()));
							if( jt->FaceOrientedOutside(it->self()) )
								for(ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
									values.push_back(kt->IntegerDF(set_id));
							else
								for(ElementArray<Node>::reverse_iterator kt = nodes.rbegin(); kt != nodes.rend(); kt++)
									values.push_back(kt->IntegerDF(set_id));
						}
						break;
					}
					default: printf("This should not happen %s\n",Element::GeometricTypeName(it->GetGeometricType()));
				}
			}
			fprintf(f,"CELLS %u %ld\n",NumberOfCells(),values.size());
			for(dynarray<Storage::integer,64>::size_type i = 0; i < values.size(); i++)
			{
				fprintf(f,"%d ",values[i]);
				if( (i+1) % 20 == 0) fprintf(f,"\n");
			}
			fprintf(f,"\n");
		}
		fprintf(f,"CELL_TYPES %u\n",NumberOfCells());
		for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); it++)
		{
			INMOST_DATA_ENUM_TYPE nnodes = VtkElementNodes(it->GetGeometricType());
			if( nnodes == ENUMUNDEF || nnodes == it->nbAdjElements(NODE) ) //nodes match - output correct type
				fprintf(f,"%d\n",VtkElementType(it->GetGeometricType()));
			else //number of nodes mismatch with expected - some topology checks must be off
				fprintf(f,"%d\n",VtkElementType(Element::MultiPolygon));
		}
		DeleteTag(set_id);
		{
			std::vector<std::string> tag_names;
			std::vector<Tag> tags;
			ListTagNames(tag_names);
			for(unsigned int i = 0; i < tag_names.size(); i++)
			{
				Tag t = GetTag(tag_names[i]);
				//printf("%s %d %d %d\n",tag_names[i].c_str(),t.isDefined(CELL),!t.isSparse(CELL),t.GetDataType() != DATA_BULK);
				if( t.isDefined(CELL) && 
            !t.isSparse(CELL) && 
            t.GetDataType() != DATA_BULK && 
            t.GetDataType() != DATA_REFERENCE &&
            t.GetDataType() != DATA_REMOTE_REFERENCE &&
					  t != CoordsTag() && 
            t != SharedTag() && 
            t != SendtoTag() && 
            t != ProcessorsTag())
				{
					//printf("added!\n");
					tags.push_back(t);
				}
			}
				
			if( !tags.empty() ) fprintf(f,"CELL_DATA %u\n",NumberOfCells());
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
            std::string type_str = "int";
            if(  tags[i].GetDataType() == DATA_REAL
#if defined(USE_AUTODIFF)
              || tags[i].GetDataType() == DATA_VARIABLE
#endif
              ) type_str = "double";
						fprintf(f,"SCALARS %s %s %d\n",tags[i].GetTagName().c_str(),type_str.c_str(),comps);
						fprintf(f,"LOOKUP_TABLE default\n");
						for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); it++)
						{
							switch( tags[i].GetDataType() )
							{
								case DATA_REAL:
								{
									Storage::real_array arr = it->RealArray(tags[i]);
									for(unsigned int m = 0; m < comps; m++) 
									{
										double val = static_cast<double>(arr[m]);
										fprintf(f,"%14e ",val != val ? -0.9999E30 : val);
									}
									fprintf(f,"\n");
								}
								break;
								case DATA_INTEGER:
								{
									Storage::integer_array arr = it->IntegerArray(tags[i]);
									for(unsigned int m = 0; m < comps; m++) fprintf(f,"%d ",arr[m]);
									fprintf(f,"\n");
								}
								break;
#if defined(USE_AUTODIFF)
                case DATA_VARIABLE:
								{
									Storage::var_array arr = it->VariableArray(tags[i]);
									for(unsigned int m = 0; m < comps; m++) 
									{
										double val = static_cast<double>(arr[m].GetValue());
										fprintf(f,"%14e ",val != val ? -0.9999E30 : val);
									}
									fprintf(f,"\n");
								}
                break;
#endif
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
				if( t.isDefined(NODE) && 
            !t.isSparse(NODE) && 
            t.GetDataType() != DATA_BULK && 
            t.GetDataType() != DATA_REFERENCE &&
            t.GetDataType() != DATA_REMOTE_REFERENCE &&
					  t != CoordsTag() && 
            t != SharedTag() && 
            t != SendtoTag() && 
            t != ProcessorsTag())
					tags.push_back(t);
			}
				
			if( !tags.empty() ) fprintf(f,"POINT_DATA %u\n",NumberOfNodes());
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
            std::string type_str = "int";
            if(  tags[i].GetDataType() == DATA_REAL
#if defined(USE_AUTODIFF)
              || tags[i].GetDataType() == DATA_VARIABLE
#endif
              ) type_str = "double";
            fprintf(f,"SCALARS %s %s %d\n",tags[i].GetTagName().c_str(),type_str.c_str(),comps);
						fprintf(f,"LOOKUP_TABLE default\n");
						for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); it++)
						{
							switch( tags[i].GetDataType() )
							{
								case DATA_REAL:
								{
									Storage::real_array arr = it->RealArray(tags[i]);
									for(unsigned int m = 0; m < comps; m++) 
									{
										double val = static_cast<double>(arr[m]);
										fprintf(f,"%14e ",(val != val ? -0.9999E30 : val));
									}
									fprintf(f,"\n");
								}
								break;
								case DATA_INTEGER:
								{
									Storage::integer_array arr = it->IntegerArray(tags[i]);
									for(unsigned int m = 0; m < comps; m++) fprintf(f,"%d ",arr[m]);
									fprintf(f,"\n");
								}
								break;
#if defined(USE_AUTODIFF)
                case DATA_VARIABLE:
								{
									Storage::var_array arr = it->VariableArray(tags[i]);
									for(unsigned int m = 0; m < comps; m++) 
									{
										double val = static_cast<double>(arr[m].GetValue());
										fprintf(f,"%14e ",(val != val ? -0.9999E30 : val));
									}
									fprintf(f,"\n");
								}
								break;
#endif
								default: continue;
							}
						}
					}
				}
			}
		}
		fclose(f);	
	}

  void Mesh::LoadVTK(std::string File)
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

		MarkerType unused_marker = CreateMarker();
		int grid_is_2d = 2;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < file_options.size(); ++k)
		{
			if( file_options[k].first == "VTK_GRID_DIMS" )
			{
				if( file_options[k].second == "AUTO" )
					grid_is_2d = 2;
				if( atoi(file_options[k].second.c_str()) == 2 )
					grid_is_2d = 1;
				else if( atoi(file_options[k].second.c_str()) == 3 )
					grid_is_2d = 0;
			}
		}

		//Determine whether there are already 3d elements so that the grid is 3d
		if( grid_is_2d == 2 && NumberOfCells() )
		{
			for(Mesh::iteratorCell it = BeginCell(); it != EndCell() && grid_is_2d == 2; ++it)
				if( it->GetElementDimension() == 3 )
					grid_is_2d = 0;
		}


		std::vector<HandleType> old_nodes(NumberOfNodes());
		{
			unsigned qq = 0;
			for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it)
				old_nodes[qq++] = *it;
		}
		if( !old_nodes.empty() ) 
		{
			std::sort(old_nodes.begin(),old_nodes.end(),CentroidComparator(this));
			//for(std::vector<HandleType>::iterator it = old_nodes.begin(); it != old_nodes.end(); ++it)
			//{
			//	Storage::real_array c = RealArrayDF(*it,CoordsTag());
			//	REPORT_VAL("coord: ",c[0] << " " << c[1] << " " << c[2]);
			//}
		}

			
		FILE * f = fopen(File.c_str(),"r");
		if( !f ) 
		{
			std::cout << __FILE__ << ":" << __LINE__ << " cannot open " << File << std::endl;
			throw BadFileName;
		}
		std::vector<Tag> datatags;
		std::vector<HandleType> newnodes;
		std::vector<HandleType> newcells;
		std::vector<int> cp;
		std::vector<int> ct;
		unsigned int state = R_VERSION;
		char readline[2048];
		int filled;
		bool binary = false;
		int read_into = NONE, read_into_cell = NONE;
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
					if( filled == 0 ) 
					{
						std::cout << __FILE__ << ":" << __LINE__ << " version information not found in " << File << std::endl;
						throw BadFile;
					}
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
					else 
					{
						std::cout << __FILE__ << ":" << __LINE__ << " unexpected data type " << readline << " in " << File << " expected BINARY or ASCII instead " << std::endl;
						throw BadFile;
					}
					break;
				}
				case R_ATTRDATA:
				{
					char dataname[1024];
					char attrname[1024];
					char attrtype[1024];
					int nentries = 1;
					filled = sscanf(readline," %s ",dataname);
					if( filled != 1 ) 
					{
						std::cout << __FILE__ << ":" << __LINE__ << " cannot read attribute name in " << File << std::endl;
						throw BadFile;
					}
					if( !strcmp(dataname,"SCALARS") )
					{
						DataType t = DATA_BULK;
						filled = sscanf(readline," SCALARS %s %s %d",attrname,attrtype,&nentries);
						if( verbosity > 0 ) printf("Reading attribute %s.\n",attrname);
						if( filled < 2 ) 
						{
							printf("%s:%d found %d arguments to SCALARS field, must be >= 2\nline:\n%s\n",__FILE__,__LINE__,filled,readline); 
							throw BadFile;
						}
						for(unsigned int i = 0; i < strlen(attrtype); i++) attrtype[i] = tolower(attrtype[i]);
						if( !strcmp(attrtype,"bit") || !strcmp(attrtype,"unsigned_char") || !strcmp(attrtype,"char") ||
							  !strcmp(attrtype,"unsigned_short") || !strcmp(attrtype,"short") || !strcmp(attrtype,"unsigned_int") ||
							  !strcmp(attrtype,"int") || !strcmp(attrtype,"unsigned_long") || !strcmp(attrtype,"long"))
							t = DATA_INTEGER;
						else if( !strcmp(attrtype,"float") || !strcmp(attrtype,"double") )
							t = DATA_REAL;
						else
						{
							std::cout << __FILE__ << ":" << __LINE__ << " unexpected data type " << attrtype << " in " << File << std::endl;
							std::cout << "expected one of: bit, unsigned_char, char, unsigned_short, short, unsigned_int, int, unsigned_long, long, float, double" << std::endl;
							throw BadFile;
						}
						if( fgets(readline,2048,f) == NULL ) 
						{
							std::cout << __FILE__ << ":" << __LINE__ << " expected LOOKUP_TABLE in " << File << std::endl;
							throw BadFile; //LOOK_UP TABLE
						}
						if( read_into == 2 )
						{
							Tag attr = CreateTag(attrname,t,read_into_cell,read_into_cell & ESET,nentries);
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newcells.size()/250),1);
							for(unsigned int it = 0; it < newcells.size(); it++)
							{

								if( t == DATA_INTEGER )
								{
									if( newcells[it] != InvalidHandle() )
									{
										Storage::integer_array attrdata = IntegerArray(newcells[it],attr);
										for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
									else for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
								}
								if( t == DATA_REAL )
								{
									if( newcells[it] != InvalidHandle() )
									{
										Storage::real_array attrdata = RealArray(newcells[it],attr);
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
									else for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; }
								}
								if( verbosity > 1 && it%report_pace == 0)
								{
									printf("data %3.1f%%\r",(it*100.0)/(1.0*newcells.size()));
									fflush(stdout);
								}
							}
							filled = fscanf(f,"\n");
							datatags.push_back(attr);
						}
						if( read_into == 1 )
						{
							Tag attr = CreateTag(attrname,t,NODE,NONE,nentries);
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newnodes.size())/250,1);
							for(unsigned int it = 0; it < newnodes.size(); it++)
							{
								if( t == DATA_INTEGER )
								{
									Storage::integer_array attrdata = IntegerArray(newnodes[it],attr);
									for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
								}
								if( t == DATA_REAL )
								{
									Storage::real_array attrdata = RealArray(newnodes[it],attr);
									for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
								}
								if( verbosity > 1 && it%report_pace == 0)
								{
									printf("data %3.1f%%\r",(it*100.0)/(1.0*newnodes.size()));
									fflush(stdout);
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
						if( verbosity > 0 ) printf("Reading attribute %s.\n",attrname);
						if( read_into == 2 )
						{
							Tag attr = CreateTag(attrname,DATA_REAL,read_into_cell,read_into_cell & ESET,nentries);
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newcells.size()/250),1);
							for(unsigned int it = 0; it < newcells.size(); it++)
							{
								if( newcells[it] != InvalidHandle() )
								{
									Storage::real_array attrdata = RealArray(newcells[it],attr);
									for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
								}
								else for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
								if( verbosity > 1 && it%100 == 0)
								{
									printf("data %3.1f%%\r",(it*report_pace)/(1.0*newcells.size()));
									fflush(stdout);
								}
							}
							filled = fscanf(f,"\n");
							datatags.push_back(attr);
						}
						if( read_into == 1 )
						{
							Tag attr = CreateTag(attrname,DATA_REAL,NODE,NONE,nentries);
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newnodes.size()/250),1);
							for(unsigned int it = 0; it < newnodes.size(); it++)
							{
								Storage::real_array attrdata = RealArray(newnodes[it],attr);
								for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
								if( verbosity > 1 && it%report_pace == 0)
								{
									printf("data %3.1f%%\r",(it*100.0)/(1.0*newnodes.size()));
									fflush(stdout);
								}
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
						if( verbosity > 0 ) printf("Reading attribute %s.\n",attrname);
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
						if( fgets(readline,2048,f) == NULL ) throw BadFile; //LOOK_UP TABLE
						if( read_into == 2 )
						{
							Tag attr = CreateTag(attrname,t,read_into_cell,read_into_cell & ESET,nentries);
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newcells.size()/250),1);
							for(unsigned int it = 0; it < newcells.size(); it++)
							{
								if( t == DATA_INTEGER )
								{
									if( newcells[it] != InvalidHandle() )
									{
										Storage::integer_array attrdata = IntegerArray(newcells[it],attr);
										for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
									else for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
								}
								if( t == DATA_REAL )
								{
									if( newcells[it] != InvalidHandle() )
									{
										Storage::real_array attrdata = RealArray(newcells[it],attr);
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
									else for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
								}
								if( verbosity > 1 && it%report_pace == 0)
								{
									printf("data %3.1f%%\r",(it*100.0)/(1.0*newcells.size()));
									fflush(stdout);
								}
							}
							filled = fscanf(f,"\n");
							datatags.push_back(attr);
						}
						if( read_into == 1 )
						{
							Tag attr = CreateTag(attrname,t,NODE,NONE,nentries);
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newnodes.size()/250),1);
							for(unsigned int it = 0; it < newnodes.size(); it++)
							{
								if( t == DATA_INTEGER )
								{
									Storage::integer_array attrdata = IntegerArray(newnodes[it],attr);
									for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
								}
								if( t == DATA_REAL )
								{
									Storage::real_array attrdata = RealArray(newnodes[it],attr);
									for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
								}
								if( verbosity > 1 && it%report_pace == 0)
								{
									printf("data %3.1f%%\r",(it*100.0)/(1.0*newnodes.size()));
									fflush(stdout);
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
						if( fgets(readline,2048,f) == NULL ) throw BadFile; //LOOK_UP TABLE
						if( read_into == 2 )
						{
							Tag attr = CreateTag(attrname,t,read_into_cell,read_into_cell & ESET,nentries);
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newcells.size()/250),1);
							for(unsigned int it = 0; it < newcells.size(); it++)
							{
								if( t == DATA_INTEGER )
								{
									if( newcells[it] != InvalidHandle() )
									{
										Storage::integer_array attrdata = IntegerArray(newcells[it],attr);
										for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									} else for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
								}
								if( t == DATA_REAL )
								{
									if( newcells[it] != InvalidHandle() )
									{
										Storage::real_array attrdata = RealArray(newcells[it],attr);
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									} else for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
								}
								if( verbosity > 1 && it%report_pace == 0)
								{
									printf("data %3.1f%%\r",(it*100.0)/(1.0*newcells.size()));
									fflush(stdout);
								}
							}
							filled = fscanf(f,"\n");
							datatags.push_back(attr);
						}
						if( read_into == 1 )
						{
							Tag attr = CreateTag(attrname,t,NODE,NONE,nentries);
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newnodes.size()/250),1);
							for(unsigned int it = 0; it < newnodes.size(); it++)
							{
								if( t == DATA_INTEGER )
								{
									Storage::integer_array attrdata = IntegerArray(newnodes[it],attr);
									for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
								}
								if( t == DATA_REAL )
								{
									Storage::real_array attrdata = RealArray(newnodes[it],attr);
									for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
								}
								if( verbosity > 1 && it%report_pace == 0)
								{
									printf("data %3.1f%%\r",(it*100.0)/(1.0*newnodes.size()));
									fflush(stdout);
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
							if( read_into == 2 )
							{
								Tag attr  = CreateTag(attrname,t,read_into_cell,read_into_cell & ESET,nentries);
								if( ntuples != newcells.size() ) printf("number of tuples in field is not equal to number of cells\n");
								unsigned report_pace = std::max<unsigned>(ntuples/250,1);
								for(unsigned int it = 0; it < ntuples; it++)
								{
									if( t == DATA_INTEGER )
									{
										if( newcells[it] != InvalidHandle() )
										{
											Storage::integer_array attrdata = IntegerArray(newcells[it],attr);
											for(unsigned int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
										} else for(unsigned int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
									}
									if( t == DATA_REAL )
									{
										if( newcells[it] != InvalidHandle() )
										{
											Storage::real_array attrdata = RealArray(newcells[it],attr);
											for(unsigned int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
										} else for(unsigned int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
									}
									if( verbosity > 1 && it%report_pace == 0)
									{
										printf("data %3.1f%%\r",(it*100.0)/(1.0*ntuples));
										fflush(stdout);
									}
								}
								filled = fscanf(f,"\n");
								datatags.push_back(attr);
							}
							if( read_into == 1 )
							{
								Tag attr = CreateTag(attrname,t,NODE,NONE,nentries);
								if( ntuples != newnodes.size() ) printf("number of tuples in field is not equal to number of nodes\n");
								unsigned report_pace = std::max<unsigned>(ntuples/250,1);
								for(unsigned int it = 0; it < ntuples; it++)
								{
									if( t == DATA_INTEGER )
									{
										Storage::integer_array attrdata = IntegerArray(newnodes[it],attr);
										for(unsigned int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
									if( t == DATA_REAL )
									{
										Storage::real_array attrdata = RealArray(newnodes[it],attr);
										for(unsigned int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
									if( verbosity > 1 && it%report_pace == 0)
									{
										printf("data %3.1f%%\r",(it*100.0)/(1.0*ntuples));
										fflush(stdout);
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
						read_into = 2;
						if( verbosity > 0 ) printf("Reading data for cells.\n");
						break;
					}
					else if( !strcmp(datatype,"POINT_DATA") )
					{
						if( npoints != newnodes.size() ) printf("number of attributes is not equal to number of nodes\n");
						state = R_ATTRDATA;
						read_into = 1;
						if( verbosity > 0 ) printf("Reading data for nodes.\n");
						break;
					}
					else
					{
						printf("%s\n",readline);
						printf("Unknown type of attributes\n");
						state = R_WAITDATA;
						read_into = 0;
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
					newnodes.resize(npoints);
					//printf("number of nodes: %d\n",npoints);
					if( verbosity > 0 ) printf("Reading %d nodes.\n",npoints);
					if( filled != 2 ) throw BadFile;
					unsigned report_pace = std::max<unsigned>(npoints/250,1);
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
							{
								std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),coords,CentroidComparator(this));
								if( it != old_nodes.end() ) 
								{
									Storage::real_array c = RealArrayDF(*it,CoordsTag());
									if( CentroidComparator(this).Compare(c.data(),coords) == 0 )
										find = static_cast<int>(it - old_nodes.begin());
								}
							}
							if( find == -1 ) 
							{
								newnodes[i] = CreateNode(coords)->GetHandle();
								SetMarker(newnodes[i],unused_marker);
							}
							else 
								newnodes[i] = old_nodes[find];
								
							if( verbosity > 1 && i%report_pace == 0)
							{
								printf("nodes %3.1f%%\r",(i*100.0)/(1.0*npoints));
								fflush(stdout);
							}
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
							{
								//REPORT_VAL("look up",coords[0] << " " << coords[1] << " " << coords[2]);
								std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),coords,CentroidComparator(this));
								if( it != old_nodes.end() ) 
								{
									Storage::real_array c = RealArrayDF(*it,CoordsTag());
									//REPORT_VAL("found",c[0] << " " << c[1] << " " << c[2]);
									if( CentroidComparator(this).Compare(c.data(),coords) == 0 )
									{
										find = static_cast<int>(it - old_nodes.begin());
										//REPORT_VAL("match",find << " " << old_nodes[find]);
									}
								}
							}
							if( find == -1 )
							{
								newnodes[i] = CreateNode(coords)->GetHandle();
								SetMarker(newnodes[i],unused_marker);
							}
							else newnodes[i] = old_nodes[find];
								
							if( verbosity > 1 && i%report_pace == 0)
							{
								printf("nodes %3.1f%%\r",(i*100.0)/(1.0*npoints));
								fflush(stdout);
							}
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
						filled = static_cast<int>(fread(&cp[0],sizeof(int),nints,f));
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
						filled = static_cast<int>(fread(&ct[0],sizeof(int),nints,f));
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

					if( grid_is_2d == 2 )
					{
						for(i = 0; i < ncells && grid_is_2d == 2; i++)
							if( ct[i] > 9 ) grid_is_2d = 0;

						grid_is_2d = 1;
					}

					if( verbosity > 0 )
					{
						switch( grid_is_2d )
						{
						case 0: std::cout << "Grid has three dimensions" << std::endl; break;
						case 1: std::cout << "Grid has two dimensions" << std::endl; break;
						case 2: std::cout << "Grid has undetermined dimension" << std::endl; break;
						}
					}

					{
						if( verbosity > 0 ) printf("Reading %d cells.\n",ncells);
						int j = 0;
						ElementArray<Node> c_nodes(this);
						ElementArray<Node> e_nodes(this);
						ElementArray<Edge> f_edges(this);
						ElementArray<Face> c_faces(this);
						ElementArray<Edge> v_edges(this);
						unsigned report_pace = std::max<unsigned>(ncells/250,1);
						newcells.resize(ncells);
						for(i = 0; i < ncells; i++)
						{
							//printf("load progress: %20.2f%%\r",(float)(i+1)/(float)ncells*100.0f);
							//fflush(stdin);
							e_nodes.clear();
							c_nodes.clear();
							f_edges.clear();
							c_faces.clear();
							v_edges.clear();
							switch(ct[i])
							{
								case 1: // VTK_VERTEX
								{
									for(int k = j+1; k < j+1+cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]],unused_marker);
									}
									j = j + 1 + cp[j];
									assert(c_nodes.size() == 1);
									newcells[i] = c_nodes.at(0);
									break;
								}
								case 2: //VTK_POLY_VERTEX
								{
									{
										std::stringstream setname;
										setname << "VTK_POLY_VERTEX_" << rand();
										ElementSet eset = CreateSet(setname.str()).first;
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											eset->PutElement(newnodes[cp[k]]);
											RemMarker(newnodes[cp[k]],unused_marker);
										}
										newcells[i] = eset->GetHandle();
									}
									j = j + 1 + cp[j];
									break;
								}
								case 3: //VTK_LINE
								{
									for(int k = j+1; k < j+1+cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]],unused_marker);
									}
									j = j + 1 + cp[j];
									assert(c_nodes.size() == 2);
									newcells[i] = CreateEdge(c_nodes).first->GetHandle();
									break;
								}
								case 4: //VTK_POLY_LINE
								{
									/*
#if defined(VTK_DEFINE_CELLS)
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
									newcells[i] = c;
#else
									*/
									{
										std::stringstream setname;
										setname << "VTK_POLY_LINE_" << rand();
										ElementSet eset = CreateSet(setname.str()).first;
										e_nodes.resize(2);
										for(int k = j+1; k < j+cp[j]; k++)
										{
											e_nodes.at(0) = newnodes[cp[k]];
											e_nodes.at(1) = newnodes[cp[k+1]];
											eset->PutElement(CreateEdge(e_nodes).first);
											RemMarker(newnodes[cp[k]],unused_marker);
										}
										RemMarker(newnodes[cp[j+cp[j]]],unused_marker);
										j = j + 1 + cp[j];
										newcells[i] = eset->GetHandle();
									}
//#endif
									break;
								}
								case 5: //VTK_TRIANGLE
								{
									for(int k = j+1; k < j+1+cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]],unused_marker);
									}
									j = j + 1 + cp[j];
									if( c_nodes.size() != 3 ) throw BadFile;
									if( grid_is_2d == 1 )
									{
										e_nodes.resize(1);
										f_edges.resize(2);
										for(unsigned int k = 0; k < 3; k++)
										{
											e_nodes.at(0) = c_nodes.at(k);
											f_edges.at(0) = CreateEdge(e_nodes).first->GetHandle();
											e_nodes.at(0) = c_nodes.at((k+1)%3);
											f_edges.at(1) = CreateEdge(e_nodes).first->GetHandle();
											c_faces.push_back(CreateFace(f_edges).first);
										}
										Cell c = CreateCell(c_faces).first;
										newcells[i] = c->GetHandle();
									}
									else newcells[i] = CreateFace(c_nodes).first->GetHandle();
									break;
								}
								case 6: //VTK_TRIANGLE_STRIP
								{
									for(int k = j+1; k < j+1+cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]],unused_marker);
									}
									j = j + 1 + cp[j];
									if( c_nodes.size() < 3 ) throw BadFile;
									/*
#if defined(VTK_DEFINE_CELLS)
									for(INMOST_DATA_ENUM_TYPE l = c_nodes.size(); l > 2; l--)
											c_faces.push_back(CreateFace(&c_nodes[l-3],3).first);
									Cell * c = CreateCell(&c_faces[0],c_faces.size(),&c_nodes[0],c_nodes.size()).first;
									newcells[i] = c;
#else //treat it as a set of faces
									*/
									{
										std::stringstream setname;
										setname << "VTK_TRIANGLE_STRIP_" << rand();
										ElementSet eset = CreateSet(setname.str()).first;
										for(ElementArray<Node>::size_type l = c_nodes.size(); l > 2; l--)
										{
											Face f = CreateFace(ElementArray<Node>(this,c_nodes.data()+l-3,c_nodes.data()+l)).first;
											eset->PutElement(f);
										}
									}
//#endif
									break;
								}
								case 7: //VTK_POLYGON
								{
									if( grid_is_2d == 1 )
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											RemMarker(newnodes[cp[k]],unused_marker);
										}
										j = j + 1 + cp[j];
										e_nodes.resize(1);
										f_edges.resize(2);
										for(ElementArray<Node>::size_type k = 0; k < c_nodes.size(); k++)
										{
											e_nodes.at(0) = c_nodes.at(k);
											f_edges.at(0) = CreateEdge(e_nodes).first->GetHandle();
											e_nodes.at(0) = c_nodes.at((k+1)%c_nodes.size());
											f_edges.at(1) = CreateEdge(e_nodes).first->GetHandle();
											c_faces.push_back(CreateFace(f_edges).first);
										}
										Cell c = CreateCell(c_faces,c_nodes).first;
										newcells[i] = c->GetHandle();
									}
									else
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											RemMarker(newnodes[cp[k]],unused_marker);
										}
										j = j + 1 + cp[j];
										newcells[i] = CreateFace(c_nodes).first->GetHandle();
									}
									break;
								}
								case 8: //VTK_PIXEL
								{
									for(int k = j+1; k < j+1+cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]],unused_marker);
									}
									j = j + 1 + cp[j];
									if( c_nodes.size() != 4 ) throw BadFile;
									HandleType temp = c_nodes.at(2);
									c_nodes.at(2) = c_nodes.at(3);
									c_nodes.at(3) = temp;
									if( grid_is_2d == 1 )
									{
										e_nodes.resize(1);
										f_edges.resize(2);
										for(int k = 0; k < 4; k++)
										{
											e_nodes.at(0) = c_nodes.at(k);
											f_edges.at(0) = CreateEdge(e_nodes).first->GetHandle();
											e_nodes.at(0) = c_nodes.at((k+1)%4);
											f_edges.at(1) = CreateEdge(e_nodes).first->GetHandle();
											c_faces.push_back(CreateFace(f_edges).first);
										}
										Cell c = CreateCell(c_faces,c_nodes).first;
										newcells[i] = c->GetHandle();
									}
									else newcells[i] = CreateFace(c_nodes).first->GetHandle();
									break;
								}
								case 9: //VTK_QUAD
								{
									for(int k = j+1; k < j+1+cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]],unused_marker);
									}
									j = j + 1 + cp[j];
									if( c_nodes.size() != 4 ) throw BadFile;
									if( grid_is_2d == 1 )
									{
										e_nodes.resize(1);
										f_edges.resize(2);
										for(int k = 0; k < 4; k++)
										{
											e_nodes.at(0) = c_nodes.at(k);
											f_edges.at(0) = CreateEdge(e_nodes).first->GetHandle();
											e_nodes.at(0) = c_nodes.at((k+1)%4);
											f_edges.at(1) = CreateEdge(e_nodes).first->GetHandle();
											c_faces.push_back(CreateFace(f_edges).first);
										}
										Cell c = CreateCell(c_faces,c_nodes).first;
										newcells[i] = c->GetHandle();
									}
									else newcells[i] = CreateFace(c_nodes).first->GetHandle();
									break;
								}
								case 10: //VTK_TETRA
								{
									for(int k = j+1; k < j+1+cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]],unused_marker);
									}
									j = j + 1 + cp[j];
									if( c_nodes.size() != 4 ) throw BadFile;
									const integer nodesnum[12] = {0,2,1,0,1,3,1,2,3,0,3,2};
									const integer sizes[4] = {3,3,3,3};
									Cell c = CreateCell(c_nodes,nodesnum,sizes,4).first;
									newcells[i] = c->GetHandle();
									break;
								}
								case 11: //VTK_VOXEL
								{
									for(int k = j+1; k < j+1+cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]],unused_marker);
									}
									j = j + 1 + cp[j];
									if (c_nodes.size() != 8) throw BadFile;
									{
										HandleType temp;
										temp = c_nodes.at(2);
										c_nodes.at(2) = c_nodes.at(3);
										c_nodes.at(3) = temp;
										temp = c_nodes.at(6);
										c_nodes.at(6) = c_nodes.at(7);
										c_nodes.at(7) = temp;
									}
									//INMOST_DATA_ENUM_TYPE nodesnum[24] = {0,4,6,2,1,3,7,5,2,6,7,3,0,1,5,4,0,2,3,1,4,5,7,6};
									const integer nodesnum[24] = {0,4,7,3,1,2,6,5,0,1,5,4,3,7,6,2,0,3,2,1,4,5,6,7};
									const integer sizes[6] = {4,4,4,4,4,4};
									Cell c = CreateCell(c_nodes,nodesnum,sizes,6).first;
									newcells[i] = c->GetHandle();
									break;
								}
								case 12: //VTK_HEXAHEDRON
								{

									for(int k = j+1; k < j+1+cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]],unused_marker);
									}
									j = j + 1 + cp[j];
									if( c_nodes.size() != 8 ) throw BadFile;
									const integer nodesnum[24] = {0,4,7,3,1,2,6,5,0,1,5,4,3,7,6,2,0,3,2,1,4,5,6,7};
									const integer sizes[6] = {4,4,4,4,4,4};
									Cell c = CreateCell(c_nodes,nodesnum,sizes,6).first;
									newcells[i] = c->GetHandle();
									break;
								}
								case 13: //VTK_WEDGE
								{
									for(int k = j+1; k < j+1+cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]],unused_marker);
									}
									j = j + 1 + cp[j];
									if( c_nodes.size() != 6 ) throw BadFile;
									const integer nodesnum[18] = {0,2,5,3,1,4,5,2,0,3,4,1,3,5,4,0,1,2};
									const integer sizes[5] = {4,4,4,3,3};
									Cell c = CreateCell(c_nodes,nodesnum,sizes,5).first;
									newcells[i] = c->GetHandle();
									break;
								}
								case 14: //VTK_PYRAMID
								{
									for(int k = j+1; k < j+1+cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]],unused_marker);
									}
									j = j + 1 + cp[j];
									if( c_nodes.size() != 5 ) throw BadFile;
									const integer nodesnum[16] = {0,4,3,0,1,4,1,2,4,3,4,2,0,3,2,1};
									const integer sizes[5] = {3,3,3,3,4};
									Cell c = CreateCell(c_nodes,nodesnum,sizes,5).first;
									newcells[i] = c->GetHandle();
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
									dynarray<integer,64> sizes(cp[k++]);
									for(dynarray<integer,64>::size_type m = 0; m < sizes.size(); m++)
									{
										sizes[m] = cp[k++];
										for (integer l = 0; l < sizes[m]; l++)
										{
											c_nodes.push_back(newnodes[cp[k++]]);
											c_nodes.back()->RemMarker(unused_marker);
										}
									}
									j = k;
									Cell c = CreateCell(c_nodes,sizes.data(),static_cast<integer>(sizes.size())).first;
									newcells[i] = c->GetHandle();
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
							if( verbosity > 1 && i%report_pace == 0)
							{
								printf("cells %3.1f%% cells %8d faces %8d edges %8d\r",(i*100.0)/(1.0*ncells),NumberOfCells(),NumberOfFaces(),NumberOfEdges());
								fflush(stdout);
							}
						}
						//ReorderEmpty(FACE | EDGE);
						state = R_ATTRIBUTES;
					}
					for(i = 0; i < ncells; i++)
					{
						read_into_cell |= GetHandleElementType(newcells[i]);
						if( GetHandleElementType(newcells[i]) == ESET )
						{
							ElementSet set = ElementSet(this,newcells[i]);
							for(ElementSet::iterator it = set->Begin(); it != set->End(); ++it)
								read_into_cell |= it->GetElementType();
						}
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
		//move data from sets to elements, delete sets
		for(std::vector<HandleType>::size_type i = 0; i < newcells.size(); i++)
		{
			if( GetHandleElementType(newcells[i]) & ESET )
			{
				ElementSet set = ElementSet(this,newcells[i]);
					
				for(ElementSet::iterator it = set->Begin(); it != set->End(); ++it)
				{
					for(std::vector<Tag>::iterator jt = datatags.begin(); jt != datatags.end(); ++jt)
					{
						if(jt->GetSize() == ENUMUNDEF ) 
						{
							std::cout << __FILE__ << ":" << __LINE__ << " there should be no variable-size data, tag " << jt->GetTagName() << std::endl;
							continue;
						}
						switch( jt->GetDataType() )
						{
						case DATA_INTEGER:
							{
								Storage::integer_array arra = it->IntegerArray(*jt);
								Storage::integer_array arrb = set->IntegerArray(*jt);
								for(enumerator k = 0; k < jt->GetSize(); k++) arra[k] = arrb[k];
							}
							break;
						case DATA_REAL:
							{
								Storage::real_array arra = it->RealArray(*jt);
								Storage::real_array arrb = set->RealArray(*jt);
								for(enumerator k = 0; k < jt->GetSize(); k++) arra[k] = arrb[k];
							}
							break;
						default:
							std::cout << __FILE__ << ":" << __LINE__ << " there should be no data of type " << DataTypeName(jt->GetDataType()) << " , tag " << jt->GetTagName() << std::endl;
							break;
						}
					}
				}
					
				Destroy(newcells[i]);
			}
		}


		{
			int count_unused = 0;
			for(unsigned q = 0; q < newnodes.size(); q++)
				if( GetMarker(newnodes[q],unused_marker) )
				{
					count_unused++;
					Destroy(newnodes[q]);
				}
			if( count_unused ) std::cout << __FILE__ << ":" << __LINE__ << " Warning: deleted " << count_unused << " unused nodes, file " << File << std::endl;
		}
		ReleaseMarker(unused_marker);
		fclose(f);
	}
}

#endif