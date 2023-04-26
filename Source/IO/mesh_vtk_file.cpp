#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif


#include "inmost.h"
#include <cfloat>

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

//static int __isnan__(double x) { return x != x; }
//static int isinf(double x) { return !isnan(x) && isnan(x - x); }
//static int __isinf__(double x) { return fabs(x) > DBL_MAX; }
//static int __isbad(double x) { return __isnan__(x) || __isinf__(x); }



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


	void skip_empty(char * readline,FILE * f)
	{
		bool empty = true;
		do
		{
			for(size_t k = 0; k < strlen(readline); ++k)
				if( !isspace(readline[k]) ) empty = false;
			if( empty )
			{
				if( fgets(readline,2048,f) == NULL ) 
					break;
			}
		} while(empty);
		
	}

	void skip_metadata(int comps, char * readline,FILE * f)
	{
		int dcnt,filled, scnt;
		//char ntype[256];
		if( !strncmp(readline,"METADATA",8) )
		{
			for(int k = 0; k < 2; ++k)
			{
				if( fgets(readline,2048,f) == NULL ) throw BadFile;
				if( !strncmp(readline,"INFORMATION",11) )
				{
					filled = sscanf(readline,"%*s %d",&dcnt); //INFORMATION %d
					if( filled != 1 ) throw BadFile;
					for(int k = 0; k < dcnt; ++k) //2 lines per data
					{
						if( fgets(readline,2048,f) == NULL )  //NAME
							throw BadFile;
						filled = fscanf(f,"%*s"); //DATA
						filled = fscanf(f,"%d",&scnt); //vector size
						if( filled != 1 ) throw BadFile;
						for(int q = 0; q < scnt; ++q) filled = fscanf(f,"%*s\n");
					}
				}
				else if( !strncmp(readline,"COMPONENT_NAMES",15) )
				{
					for(int k = 0; k < comps; ++k) //2 lines per data
					{
						if( fgets(readline,2048,f) == NULL ) 
							throw BadFile;
					}
				}
			}
		}
	}

	std::string Space2Underscore(const std::string & inp)
	{
		std::string ret = inp;
		for (size_t k = 0; k < ret.size(); ++k)
		{
			if (ret[k] == ' ') ret[k] = '_';
		}
		return ret;
	}
  
	int VtkElementType(ElementType t)
	{
		switch(t)
		{
			case Element::Line: return 3;
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
			case Element::Line: return 2;
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
			std::cout << "VTK file supports 3 dimensions max" << std::endl;
			return;
		}
		std::ofstream f(File.c_str());
		//~ FILE * f = fopen(File.c_str(),"w");
		//~ if( !f ) throw BadFileName;
		if( f.fail() ) throw BadFileName;
		f << "# vtk DataFile Version 3.0\n";
		if (this->GetFileOption("VTK_OUTPUT_FACES") == "1") f << "VTK_OUTPUT_FACES file is written by INMOST\n";
		else f << "file is written by INMOST\n";
		f << "ASCII\n";
		f << "DATASET UNSTRUCTURED_GRID\n";

		bool output_faces = false;
		bool keep_ghost = false;
		for (INMOST_DATA_ENUM_TYPE k = 0; k < file_options.size(); ++k)
		{
			if (file_options[k].first == "VTK_OUTPUT_FACES")
			{
				output_faces = true;
			}
			else if (file_options[k].first == "KEEP_GHOST")
			{
				keep_ghost = true;
			}
		}
		
		std::set< std::string > nosave, saveonly;
		
		nosave = TagOptions("nosave");
		saveonly = TagOptions("saveonly");
		
		//ReorderEmpty(CELL | NODE);
		Tag set_id = CreateTag("PROTECTED_TEMPORARY_ELEMENT_ID",DATA_INTEGER,CELL |FACE| NODE,NONE,1);
		integer num_cells = 0, num_faces = 0, num_nodes = 0;
		MarkerType used = CreateMarker();
		for (Mesh::iteratorCell it = BeginCell(); it != EndCell(); ++it)
			if (keep_ghost || it->GetStatus() != Element::Ghost)
			{
				it->IntegerDF(set_id) = num_cells++;
				it->getAdjElements(NODE).SetMarker(used);
			}
			else it->IntegerDF(set_id) = -1;
		if (output_faces)
		{
			for (Mesh::iteratorFace it = BeginFace(); it != EndFace(); ++it)
				if (keep_ghost || it->GetStatus() != Element::Ghost)
				{
					it->IntegerDF(set_id) = num_faces++;
					it->getAdjElements(NODE).SetMarker(used);
				}
				else it->IntegerDF(set_id) = -1;
		}
		for (Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it)
			it->IntegerDF(set_id) = it->GetMarker(used) ? num_nodes++ : -1;

		ReleaseMarker(used, NODE);
		
		f << "POINTS " << num_nodes << " double\n";
		for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); it++) if( it->IntegerDF(set_id) != -1 )
		{
			Storage::real_array coords = it->RealArray(CoordsTag());
			for(integer i = 0; i < dim; i++) 
			{
				double temp = coords[i];
				f << temp << " ";
			}
			for(integer i = dim; i < 3; i++)
				f << "0 ";
			f << "\n";
		}
		{
			std::vector<int> values;
			for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); it++) if( it->IntegerDF(set_id) != -1 )
			{	
				switch(it->GetGeometricType())
				{
					case Element::Tri:
					case Element::Quad:
					//case Element::MultiLine:
					//case Element::Polygon:
					case Element::Tet:
					case Element::Hex:
					case Element::Pyramid:
					case Element::Prism:
					{
						ElementArray<Node> nodes(this);// = it->getNodes();
						RestoreCellNodes(*it,nodes);
						if( nodes.size() != VtkElementNodes(it->GetGeometricType()) ) goto safe_output;
						values.push_back(static_cast<integer>(nodes.size()));
						for(ElementArray<Node>::iterator jt = nodes.begin(); jt != nodes.end(); jt++)
							values.push_back(jt->IntegerDF(set_id));
						break;
					}
					case Element::MultiLine:
					case Element::Polygon:
					{
						ElementArray<Node> nodes = it->getNodes();
						values.push_back(static_cast<integer>(nodes.size()));
						for(ElementArray<Node>::iterator jt = nodes.begin(); jt != nodes.end(); jt++)
							values.push_back(jt->IntegerDF(set_id));
						break;
					}
					case Element::Polyhedron:
					case Element::MultiPolygon:
					{
safe_output:
						
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
					default: std::cout << __FILE__ << ":" << __LINE__ << " This should not happen " << Element::GeometricTypeName(it->GetGeometricType()) << std::endl;
				}
			}
			if (output_faces)
			{
				for (Mesh::iteratorFace it = BeginFace(); it != EndFace(); it++) if( it->IntegerDF(set_id) != -1)
				{
					switch (it->GetGeometricType())
					{
					case Element::Line:
					{
						ElementArray<Node> nodes = it->getNodes();
						values.push_back(static_cast<integer>(nodes.size()));
						for (ElementArray<Node>::iterator jt = nodes.begin(); jt != nodes.end(); jt++)
							values.push_back(jt->IntegerDF(set_id));
						break;
					}
					case Element::Tri:
					case Element::Quad:
					case Element::MultiLine:
					case Element::Polygon:
					{
						ElementArray<Node> nodes = it->getNodes();
						values.push_back(static_cast<integer>(nodes.size()));
						for (ElementArray<Node>::iterator jt = nodes.begin(); jt != nodes.end(); jt++)
							values.push_back(jt->IntegerDF(set_id));
						break;
					}
					default: std::cout << __FILE__ << ":" << __LINE__ << " This should not happen " << Element::GeometricTypeName(it->GetGeometricType()) << std::endl;
					}
				}
				f << "CELLS " << num_cells + num_faces << " " << values.size() << "\n";
			}
			else
				f << "CELLS " << num_cells << " " << values.size() << "\n";

			for(size_t i = 0; i < values.size(); i++)
			{
				f << values[i] << " ";
				if( (i+1) % 20 == 0) f << "\n";
			}
			f << "\n";
		}
		if (output_faces)
			f << "CELL_TYPES " << num_cells + num_faces << "\n";
		else 			
			f << "CELL_TYPES " << num_cells << "\n";

		for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); it++) if( it->IntegerDF(set_id) != -1)
		{
			INMOST_DATA_ENUM_TYPE nnodes = VtkElementNodes(it->GetGeometricType());
			if( nnodes == ENUMUNDEF || nnodes == it->nbAdjElements(NODE) ) //nodes match - output correct type
				f << VtkElementType(it->GetGeometricType()) << "\n";
			else //number of nodes mismatch with expected - some topology checks must be off
				f << VtkElementType(Element::MultiPolygon) << "\n";
		}
		if (output_faces)
		{
			for (Mesh::iteratorFace it = BeginFace(); it != EndFace(); it++) if(it->IntegerDF(set_id) != -1)
			{
				INMOST_DATA_ENUM_TYPE nnodes = VtkElementNodes(it->GetGeometricType());
				if (nnodes == ENUMUNDEF || nnodes == it->nbAdjElements(NODE)) //nodes match - output correct type
					f << VtkElementType(it->GetGeometricType()) << "\n";
				else //number of nodes mismatch with expected - some topology checks must be off
					f << VtkElementType(Element::MultiPolygon) << "\n";
			}
		}
		
		{
			std::vector<std::string> tag_names;
			std::vector<Tag> tags;
			ListTagNames(tag_names);
			for(unsigned int i = 0; i < tag_names.size(); i++)
			{
				Tag t = GetTag(tag_names[i]);
				if (((t.isDefined(CELL) /*&& !t.isSparse(CELL)*/)
					|| (t.isDefined(FACE) && output_faces)) &&
						t.GetPrint() && //Temporary solution: @see Mesh::file_option
						t.GetDataType() != DATA_BULK && 
						t.GetDataType() != DATA_REFERENCE &&
						t.GetDataType() != DATA_REMOTE_REFERENCE &&
						t != CoordsTag() && 
						t != SharedTag() && 
						t != SendtoTag() && 
						t != ProcessorsTag())
				{
					bool skip = CheckSaveSkip(tag_names[i],nosave,saveonly);
					if( !skip )
						tags.push_back(t);
				}
			}
				
			if (!tags.empty() && output_faces) 
			{ 
				f << "CELL_DATA " <<num_cells + num_faces << "\n"; 
			}
			else if (!tags.empty()) f << "CELL_DATA " << num_cells << "\n";

			for(INMOST_DATA_ENUM_TYPE i = 0; i < tags.size(); i++)
			{
				INMOST_DATA_ENUM_TYPE comps = tags[i].GetSize();
				if( comps == ENUMUNDEF )
					continue;
				else
				{
					{
						std::string type_str = "int";
						if(  tags[i].GetDataType() == DATA_REAL
#if defined(USE_AUTODIFF)
						   || tags[i].GetDataType() == DATA_VARIABLE
#endif
						   ) type_str = "double";
						f << "SCALARS " << Space2Underscore(tags[i].GetTagName()) << " " << type_str << " " << comps << "\n";
						f << "LOOKUP_TABLE default\n";
						for (Mesh::iteratorCell it = BeginCell(); it != EndCell(); it++) if (it->IntegerDF(set_id) != -1)
						{
							switch (tags[i].GetDataType())
							{
							case DATA_REAL:
							{
								if (tags[i].isDefined(CELL) && it->HaveData(tags[i]))
								{
									Storage::real_array arr = it->RealArray(tags[i]);
									for (unsigned int m = 0; m < comps; m++)
									{
										double val = static_cast<double>(arr[m]);
										f << (__isbad(val) ? -0.9999E30 : val) << " ";
									}
								}
								else for (unsigned int m = 0; m < comps; m++) f << -0.9999E30 << " ";
								f << "\n";
							}
							break;
							case DATA_INTEGER:
							{
								if (tags[i].isDefined(CELL) && it->HaveData(tags[i]))
								{
									Storage::integer_array arr = it->IntegerArray(tags[i]);
									for (unsigned int m = 0; m < comps; m++) f << arr[m] << " ";
								}
								else for (unsigned int m = 0; m < comps; m++) f << INT_MIN << " ";
								f << "\n";
							}
							break;
#if defined(USE_AUTODIFF)
							case DATA_VARIABLE:
							{
								if (tags[i].isDefined(CELL) && it->HaveData(tags[i]))
								{
									Storage::var_array arr = it->VariableArray(tags[i]);
									for (unsigned int m = 0; m < comps; m++)
									{
										double val = static_cast<double>(arr[m].GetValue());
										f << (__isbad(val) ? -0.9999E30 : val) << " ";
									}
								}
								else for (unsigned int m = 0; m < comps; m++) f << -0.9999E30 << " ";
								f << "\n";
							}
							break;
#endif
							default: continue;
							}
						}
							


						if (output_faces)
						{
							for (Mesh::iteratorFace it = BeginFace(); it != EndFace(); it++) if (it->IntegerDF(set_id) != -1)
							{
								switch (tags[i].GetDataType())
								{
								case DATA_REAL:
								{
									if (tags[i].isDefined(FACE) && it->HaveData(tags[i]))
									{
										Storage::real_array arr = it->RealArray(tags[i]);
										for (unsigned int m = 0; m < comps; m++) f << arr[m] << " ";
									}
									else for (unsigned int m = 0; m < comps; m++) f << -0.9999E30 << " ";
									f << "\n";
								}
								break;
								case DATA_INTEGER:
								{
									if (tags[i].isDefined(FACE) && it->HaveData(tags[i]))
									{
										Storage::integer_array arr = it->IntegerArray(tags[i]);
										for (unsigned int m = 0; m < comps; m++) f << arr[m] << " ";
									}
									else for (unsigned int m = 0; m < comps; m++) f << INT_MIN << " ";
									f << "\n";
								}
								break;
#if defined(USE_AUTODIFF)
								case DATA_VARIABLE:
								{
									if (tags[i].isDefined(FACE) && it->HaveData(tags[i]))
									{
										Storage::var_array arr = it->VariableArray(tags[i]);
										for (unsigned int m = 0; m < comps; m++)
										{
											double val = static_cast<double>(arr[m].GetValue());
											f << (__isbad(val) ? -0.9999E30 : val) << " ";
										}
									}
									else for (unsigned int m = 0; m < comps; m++) f << -0.9999E30 << " ";
									f << "\n";
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
		}
		{
			std::vector<std::string> tag_names;
			std::vector<Tag> tags;
			ListTagNames(tag_names);
			for(unsigned int i = 0; i < tag_names.size(); i++)
			{
				Tag t = GetTag(tag_names[i]);
				if (t.isDefined(NODE) && /*!t.isSparse(NODE) &&*/
					t.GetDataType() != DATA_BULK &&
					t.GetDataType() != DATA_REFERENCE &&
					t.GetDataType() != DATA_REMOTE_REFERENCE &&
					t != CoordsTag() &&
					t != SharedTag() &&
					t != SendtoTag() &&
					t != ProcessorsTag())
				{
					bool skip = CheckSaveSkip(tag_names[i], nosave, saveonly);
					if (!skip)
						tags.push_back(t);
				}
			}
				
			if( !tags.empty() ) f << "POINT_DATA " << num_nodes << "\n";
			for(INMOST_DATA_ENUM_TYPE i = 0; i < tags.size(); i++)
			{
				INMOST_DATA_ENUM_TYPE comps = tags[i].GetSize();
				if( comps == ENUMUNDEF )
					continue;
				else
				{
					{
						std::string type_str = "int";
						if(  tags[i].GetDataType() == DATA_REAL
#if defined(USE_AUTODIFF)
						   || tags[i].GetDataType() == DATA_VARIABLE
#endif
						   ) type_str = "double";
						f << "SCALARS " << Space2Underscore(tags[i].GetTagName()) << " " << type_str << " " << comps << "\n";
						f << "LOOKUP_TABLE default\n";
						for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); it++) if (it->IntegerDF(set_id) != -1)
						{
							switch( tags[i].GetDataType() )
							{
								case DATA_REAL:
								{
									if (it->HaveData(tags[i]))
									{
										Storage::real_array arr = it->RealArray(tags[i]);
										for(unsigned int m = 0; m < comps; m++)
										{
											double val = static_cast<double>(arr[m]);
											f << (__isbad(val) ? -0.9999E30 : val) << " ";
										}
										f << "\n";
									}
									else for(unsigned int m = 0; m < comps; m++) f << -0.9999E30 << " ";
								}
								break;
								case DATA_INTEGER:
								{
									if (it->HaveData(tags[i]))
									{
										Storage::integer_array arr = it->IntegerArray(tags[i]);
										for(unsigned int m = 0; m < comps; m++) f << arr[m] << " ";
										f << "\n";
									}
									else for (unsigned int m = 0; m < comps; m++) f << INT_MIN << " ";
								}
								break;
#if defined(USE_AUTODIFF)
								case DATA_VARIABLE:
								{
									if (it->HaveData(tags[i]))
									{
										Storage::var_array arr = it->VariableArray(tags[i]);
										for(unsigned int m = 0; m < comps; m++)
										{
											double val = static_cast<double>(arr[m].GetValue());
											f << (__isbad(val) ? -0.9999E30 : val) << " ";
										}
										f << "\n";
									}
									else for(unsigned int m = 0; m < comps; m++) f << -0.9999E30 << " ";
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
		DeleteTag(set_id);
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
					std::cout << __FILE__ << ":" << __LINE__ << " Unknown verbosity option: " << file_options[k].second << std::endl;
					verbosity = 1;
				}
			}
		}
		
		std::set< std::string > noload, loadonly;
		
		noload = TagOptions("noload");
		loadonly = TagOptions("loadonly");
		

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
		bool read_next_line = true;
		while( state != R_QUIT )
		{
			
			if( read_next_line )
			{
				if( fgets(readline,2048,f) == NULL )
				{
					state = R_QUIT;
					continue;
				}
			} else read_next_line = true;
			//if( readline[strlen(readline)-1] == '\n' )
			//	readline[strlen(readline)-1] = '\0';
			while( strlen(readline) != 0 && isspace(readline[strlen(readline)-1]) )
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
					if (!strncmp(readline, "VTK_OUTPUT_FACES", 16)) this->SetFileOption("VTK_OUTPUT_FACES", "1");
					state = R_DATATYPE;
					break;
				}
				case R_DATATYPE:
				{
					if( !strncmp(readline,"BINARY",6) )
					{
						binary = 1;
						state = R_WAITDATA;
					}
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
						bool skip = CheckLoadSkip(attrname,noload,loadonly);
						if( verbosity > 0 ) std::cout << (skip?"Skipping":"Reading") << " scalar attribute " << attrname << std::endl;
						if( filled < 2 ) 
						{
							std::cout << __FILE__ << ":" << __LINE__ << " found " << filled << " arguments to SCALARS field, must be >= 2\nline:\n" << readline << "\n";
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
							std::cout << __FILE__ << ":" << __LINE__ << " unexpected data type " << attrtype << " in " << File << " for " << attrname << std::endl;
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
							bool sparse = (std::string(attrname).substr(0,9) != "PROTECTED");
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newcells.size()/250),1);
							if( skip )
							{
								for(unsigned int it = 0; it < newcells.size(); it++)
								{

									if( t == DATA_INTEGER )
										for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
									if( t == DATA_REAL )
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; }
									if( verbosity > 1 && it%report_pace == 0)
									{
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newcells.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
							else
							{
								Tag attr = CreateTag(attrname,t,read_into_cell,sparse ? read_into_cell & ~CELL : NONE,nentries);
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
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newcells.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
								datatags.push_back(attr);
							}
						}
						if( read_into == 1 )
						{
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newnodes.size())/250,1);
							if( skip )
							{
								for(unsigned int it = 0; it < newnodes.size(); it++)
								{
									if( t == DATA_INTEGER )
										for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
									if( t == DATA_REAL )
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
									if( verbosity > 1 && it%report_pace == 0)
									{
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newnodes.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
							else
							{
								Tag attr = CreateTag(attrname,t,NODE,NONE,nentries);
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
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newnodes.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
							//filled = fscanf(f,"\n");
						}
						if( fgets(readline,2048,f) == NULL ) throw BadFile;
						skip_empty(readline,f);
						skip_metadata(3,readline,f);
						skip_empty(readline,f);
						read_next_line = false;
						
						break;
					}
					else if( !strcmp(dataname,"COLOR_SCALARS") )
					{	
						char attrname[1024];
						filled = sscanf(readline," COLOR_SCALARS %s %d",attrname,&nentries);
						bool skip = CheckLoadSkip(attrname,noload,loadonly);
						if( verbosity > 0 ) std::cout << (skip?"Skipping":"Reading") << " color attribute " << attrname << "\n";
						if( read_into == 2 )
						{
							bool sparse = (std::string(attrname).substr(0,9) != "PROTECTED");
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newcells.size()/250),1);
							if( skip )
							{
								for(unsigned int it = 0; it < newcells.size(); it++)
								{
									for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
									if( verbosity > 1 && it%report_pace == 0)
									{
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newcells.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
							else
							{
								Tag attr = CreateTag(attrname,DATA_REAL,read_into_cell,sparse ? read_into_cell & ~CELL : NONE,nentries);
								for(unsigned int it = 0; it < newcells.size(); it++)
								{
									if( newcells[it] != InvalidHandle() )
									{
										Storage::real_array attrdata = RealArray(newcells[it],attr);
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									}
									else for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
									if( verbosity > 1 && it%report_pace == 0)
									{
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newcells.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
								datatags.push_back(attr);
							}
						}
						if( read_into == 1 )
						{
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newnodes.size()/250),1);
							if( skip ) 
							{
								for(unsigned int it = 0; it < newnodes.size(); it++)
								{
									for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
									if( verbosity > 1 && it%report_pace == 0)
									{
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newnodes.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
							else
							{
								Tag attr = CreateTag(attrname,DATA_REAL,NODE,NONE,nentries);
								for(unsigned int it = 0; it < newnodes.size(); it++)
								{
									Storage::real_array attrdata = RealArray(newnodes[it],attr);
									for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; attrdata[jt] = temp;}
									if( verbosity > 1 && it%report_pace == 0)
									{
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newnodes.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
							//filled = fscanf(f,"\n");
						}
						if( fgets(readline,2048,f) == NULL ) throw BadFile;
						skip_empty(readline,f);
						skip_metadata(3,readline,f);
						skip_empty(readline,f);
						read_next_line = false;
					
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
						bool skip = CheckLoadSkip(attrname,noload,loadonly);
						if( verbosity > 0 ) std::cout << (skip? "Skipping":"Reading") << " " << dataname << " attribute " << attrname << "n";
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
							bool sparse = (std::string(attrname).substr(0,9) != "PROTECTED");
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newcells.size()/250),1);
							if( skip )
							{
								for(unsigned int it = 0; it < newcells.size(); it++)
								{
									if( t == DATA_INTEGER )
										for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
									if( t == DATA_REAL )
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
									if( verbosity > 1 && it%report_pace == 0)
									{
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newcells.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
							else
							{
								Tag attr = CreateTag(attrname,t,read_into_cell,sparse ? read_into_cell & ~CELL : NONE,nentries);
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
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newcells.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
								datatags.push_back(attr);
							}
						}
						if( read_into == 1 )
						{
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newnodes.size()/250),1);
							if( skip )
							{
								for(unsigned int it = 0; it < newnodes.size(); it++)
								{
									if( t == DATA_INTEGER )
										for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; }
									if( t == DATA_REAL )
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; }
									if( verbosity > 1 && it%report_pace == 0)
									{
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newnodes.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
							else
							{
								Tag attr = CreateTag(attrname,t,NODE,NONE,nentries);
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
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newnodes.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
							//filled = fscanf(f,"\n");
						}
						
						if( fgets(readline,2048,f) == NULL ) throw BadFile;
						skip_empty(readline,f);
						skip_metadata(3,readline,f);
						skip_empty(readline,f);
						read_next_line = false;
						
						break;
					}
					else if( !strcmp(dataname,"TEXTURE_COORDINATES") )
					{
						char attrname[1024];
						char attrtype[1024];
						int nentries;
						DataType t = DATA_BULK;
						filled = sscanf(readline," TEXTURE_COORDINATES %s %d %s",attrname,&nentries,attrtype);
						bool skip = CheckLoadSkip(attrname,noload,loadonly);
						if( verbosity > 0 ) std::cout << (skip?"Skipping":"Reading") << " texture coordinates attribute " << attrname << "\n";
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
							bool sparse = (std::string(attrname).substr(0,9) != "PROTECTED");
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newcells.size()/250),1);
							if( skip )
							{
								for(unsigned int it = 0; it < newcells.size(); it++)
								{
									if( t == DATA_INTEGER )
										for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
									if( t == DATA_REAL )
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
									if( verbosity > 1 && it%report_pace == 0)
									{
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newcells.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
							else
							{
								Tag attr = CreateTag(attrname,t,read_into_cell,sparse ? read_into_cell & ~CELL : NONE,nentries);
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
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newcells.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
								datatags.push_back(attr);
							}
						}
						if( read_into == 1 )
						{
							unsigned report_pace = std::max<unsigned>(static_cast<unsigned>(newnodes.size()/250),1);
							if( skip )
							{
								for(unsigned int it = 0; it < newnodes.size(); it++)
								{
									if( t == DATA_INTEGER )
										for(int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile; }
									if( t == DATA_REAL )
										for(int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile; }
									if( verbosity > 1 && it%report_pace == 0)
									{
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newnodes.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
							else
							{
								Tag attr = CreateTag(attrname,t,NODE,NONE,nentries);
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
										std::ios save(NULL);
										save.copyfmt(std::cout);
										std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*newnodes.size()) << "%\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
							//filled = fscanf(f,"\n");
						}
						
						if( fgets(readline,2048,f) == NULL ) throw BadFile;
						skip_empty(readline,f);
						skip_metadata(3,readline,f);
						skip_empty(readline,f);
						read_next_line = false;
						
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
							if( read_next_line && fgets(readline,2048,f) == NULL ) throw BadFile;
							filled = sscanf(readline," %s %u %u %s\n",attrname,&nentries,&ntuples,attrtype);
							bool skip = CheckLoadSkip(attrname,noload,loadonly);
							if( verbosity > 0 ) std::cout << (skip?"Skipping":"Reading") << " field " << i << " / " << nfields << " attribute " << attrname << "\n";
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
								bool sparse = (std::string(attrname).substr(0,9) != "PROTECTED");
								if( ntuples != newcells.size() ) std::cout << __FILE__ << ":" << __LINE__ << " number of tuples in field is not equal to number of cells" << std::endl;
								unsigned report_pace = std::max<unsigned>(ntuples/250,1);
								if( skip )
								{
									for(unsigned int it = 0; it < ntuples; it++)
									{
										if( t == DATA_INTEGER )
											for(unsigned int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
										if( t == DATA_REAL )
											for(unsigned int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
										if( verbosity > 1 && it%report_pace == 0)
										{
											std::ios save(NULL);
											save.copyfmt(std::cout);
											std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*ntuples) << "%\r" << std::flush;
											std::cout.copyfmt(save);
										}
									}
								}
								else
								{
									Tag attr = CreateTag(attrname,t,read_into_cell,sparse ? read_into_cell & ~CELL : NONE,nentries);
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
											std::ios save(NULL);
											save.copyfmt(std::cout);
											std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*ntuples) << "%\r" << std::flush;
											std::cout.copyfmt(save);
										}
									}
									datatags.push_back(attr);
								}
							}
							if( read_into == 1 )
							{
								if( ntuples != newnodes.size() ) std::cout << __FILE__ << ":" << __LINE__ << " number of tuples in field is not equal to number of nodes" << std::endl;
								unsigned report_pace = std::max<unsigned>(ntuples/250,1);
								if( skip )
								{
									for(unsigned int it = 0; it < ntuples; it++)
									{
										if( t == DATA_INTEGER )
											for(unsigned int jt = 0; jt < nentries; jt++) {int temp; filled = fscanf(f," %d",&temp); if(filled != 1 ) throw BadFile;}
										if( t == DATA_REAL )
											for(unsigned int jt = 0; jt < nentries; jt++) {double temp; filled = fscanf(f," %lf",&temp); if(filled != 1 ) throw BadFile;}
										if( verbosity > 1 && it%report_pace == 0)
										{
											std::ios save(NULL);
											save.copyfmt(std::cout);
											std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*ntuples) << "%\r" << std::flush;
											std::cout.copyfmt(save);
										}
									}
								}
								else
								{
									Tag attr = CreateTag(attrname,t,NODE,NONE,nentries);
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
											std::ios save(NULL);
											save.copyfmt(std::cout);
											std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (it*100.0)/(1.0*ntuples) << "%\r" << std::flush;
											std::cout.copyfmt(save);
										}
									}
								}
								//filled = fscanf(f,"\n");
							}
							
							if( fgets(readline,2048,f) == NULL ) throw BadFile;
							skip_empty(readline,f);
							skip_metadata(3,readline,f);
							skip_empty(readline,f);
							read_next_line = false;
						}
						break;
					}
					else
					{
						state = R_ATTRIBUTES;
					}
				}
				//FALLTHROUGH
				case R_ATTRIBUTES:
				{
					char datatype[1024];
					unsigned int npoints;
					filled = sscanf(readline," %s %u",datatype,&npoints);
					if( filled != 2 ) throw BadFile;
					if( !strcmp(datatype,"CELL_DATA") )
					{
						if( npoints != newcells.size() ) std::cout << __FILE__ << ":" << __LINE__ << " number of entries in the attribute is not equal to number of cells" << std::endl;
						state = R_ATTRDATA;
						read_into = 2;
						if( verbosity > 0 ) std::cout << "Reading data for cells." << std::endl;
						break;
					}
					else if( !strcmp(datatype,"POINT_DATA") )
					{
						if( npoints != newnodes.size() ) std::cout << __FILE__ << ":" << __LINE__ << " number of entries in the attribute is not equal to number of nodes" << std::endl;
						state = R_ATTRDATA;
						read_into = 1;
						if( verbosity > 0 ) std::cout << "Reading data for nodes." << std::endl;
						break;
					}
					else
					{
						std::cout << readline << std::endl;
						std::cout << "Unknown type of attributes" << std::endl;
						state = R_WAITDATA;
						read_into = 0;
					}
				}
				//FALLTHROUGH
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
					
					while( !strncmp(readline,"FIELD",5) )
					{
						char dataname[1024], fieldname[1024];
						int nfields, ntuples, ncomps, nfields_done = 0;
						filled = sscanf(readline,"%*s %s %d",fieldname,&nfields);
						
						while( nfields_done < nfields )
						{
							if( fgets(readline,2048,f) == NULL )
							{
								state = R_QUIT;
								continue;
							}
							while( strlen(readline) != 0 && isspace(readline[strlen(readline)-1]) )
								readline[strlen(readline)-1] = '\0';
							if( strlen(readline) == 0 ) continue;
							
							sscanf(readline,"%s %d %d %s",dataname,&ncomps,&ntuples,datatype);
							DataType t;
							for(unsigned int i = 0; i < strlen(datatype); i++) datatype[i] = tolower(datatype[i]);
							if( !strcmp(datatype,"bit") || !strcmp(datatype,"unsigned_char") || !strcmp(datatype,"char") ||
								  !strcmp(datatype,"unsigned_short") || !strcmp(datatype,"short") || !strcmp(datatype,"unsigned_int") ||
								  !strcmp(datatype,"int") || !strcmp(datatype,"unsigned_long") || !strcmp(datatype,"long"))
								t = DATA_INTEGER;
							else if( !strcmp(datatype,"float") || !strcmp(datatype,"double") )
								t = DATA_REAL;
							else
								throw BadFile;
							Tag field_data = CreateTag(std::string(fieldname)+"_"+std::string(dataname),t,MESH,NONE,ncomps*ntuples);
							if( t == DATA_REAL )
							{
								for(int q = 0; q < ncomps*ntuples; ++q)
								{
									double tmp;
									filled = fscanf(f,"%lf ", &tmp);
									if( filled ) RealArray(GetHandle(),field_data)[q] = static_cast<Storage::real>(tmp);
								}
							}
							else
							{
								for(int q = 0; q < ncomps*ntuples; ++q)
								{
									int tmp;
									filled = fscanf(f,"%d ", &tmp);
									if( filled ) IntegerArray(GetHandle(),field_data)[q] = static_cast<Storage::integer>(tmp);
								}
							}

							
							nfields_done++;
						}
						do
						{
							if( fgets(readline,2048,f) == NULL )
							{
								state = R_QUIT;
								continue;
							}
							while( strlen(readline) != 0 && isspace(readline[strlen(readline)-1]) )
								readline[strlen(readline)-1] = '\0';
						}
						while( strlen(readline) == 0 );
					}
					
					
					filled = sscanf(readline,"%*s %d %s",&npoints,datatype);
					newnodes.resize(npoints);
					if( verbosity > 0 ) std::cout << "Reading " << npoints << " nodes." << std::endl;
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
								std::ios save(NULL);
								save.copyfmt(std::cout);
								std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (i*100.0)/(1.0*npoints) << "%\r" << std::flush;
								std::cout.copyfmt(save);
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
								std::ios save(NULL);
								save.copyfmt(std::cout);
								std::cout << "data " << std::setw(5) << std::fixed << std::setprecision(1) << (i*100.0)/(1.0*npoints) << "%\r" << std::flush;
								std::cout.copyfmt(save);
							}
							i++;
						}
						//filled = fscanf(f,"\n");
					}
					
					if( fgets(readline,2048,f) == NULL ) throw BadFile;
					skip_empty(readline,f);
					skip_metadata(3,readline,f);
					skip_empty(readline,f);
					
					filled = sscanf(readline,"%*s %d %d\n",&ncells,&nints);
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

					if( grid_is_2d == 2 && ncells )
					{
						for(i = 0; i < ncells && grid_is_2d == 2; i++)
							if( ct[i] > 9 ) grid_is_2d = 0;
						//this will create all 2d elements as cells
						//grid_is_2d = 1; 
						if (grid_is_2d == 2) grid_is_2d = 1;
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
					
					if( grid_is_2d && old_nodes.empty() ) SetDimensions(2);

					{
						if( verbosity > 0 ) std::cout << "Reading " << ncells << " cells." << std::endl;
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
									if( grid_is_2d )
									{
										c_nodes.resize(1);
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.at(0) = newnodes[cp[k]];
											RemMarker(newnodes[cp[k]],unused_marker);
											f_edges.push_back(CreateEdge(c_nodes).first);
										}
										j = j + 1 + cp[j];
										assert(f_edges.size() == 2);
										newcells[i] = CreateFace(f_edges).first->GetHandle();
									}
									else
									{
										for(int k = j+1; k < j+1+cp[j]; k++)
										{
											c_nodes.push_back(newnodes[cp[k]]);
											RemMarker(newnodes[cp[k]],unused_marker);
										}
										j = j + 1 + cp[j];
										assert(c_nodes.size() == 2);
										newcells[i] = CreateEdge(c_nodes).first->GetHandle();
									}
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
										//Cell c = CreateCell(c_faces).first;
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
										//Cell c = CreateCell(c_faces).first;
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
										//Cell c = CreateCell(c_faces).first;
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
									Cell c = CreateCell(c_nodes,nodesnum,sizes,4,c_nodes).first;
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
									Cell c = CreateCell(c_nodes,nodesnum,sizes,6,c_nodes).first;
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
									Cell c = CreateCell(c_nodes,nodesnum,sizes,6,c_nodes).first;
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
									//const integer nodesnum[18] = { 0, 3, 5, 2, 0, 1, 4, 3, 1, 4, 5, 2, 3, 4, 5, 0, 2, 1 };
									const integer sizes[5] = {4,4,4,3,3};
									Cell c = CreateCell(c_nodes,nodesnum,sizes,5,c_nodes).first;
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
									Cell c = CreateCell(c_nodes,nodesnum,sizes,5,c_nodes).first;
									newcells[i] = c->GetHandle();
									break;
								}
								case 15: //VTK_PENTAGONAL_PRISM
								{
									for (int k = j + 1; k < j + 1 + cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]], unused_marker);
									}
									j = j + 1 + cp[j];
									if (c_nodes.size() != 10) throw BadFile;
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
									const integer sizes[7] = { 4, 4, 4, 4, 4, 5, 5 };
									Cell c = CreateCell(c_nodes, nodesnum, sizes, 7,c_nodes).first;
									newcells[i] = c->GetHandle();
									break;
								}
								case 16: //VTK_HEXAGONAL_PRISM
								{
									for (int k = j + 1; k < j + 1 + cp[j]; k++)
									{
										c_nodes.push_back(newnodes[cp[k]]);
										RemMarker(newnodes[cp[k]], unused_marker);
									}
									j = j + 1 + cp[j];
									if (c_nodes.size() != 12) throw BadFile;
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
									Cell c = CreateCell(c_nodes, nodesnum, sizes, 8,c_nodes).first;
									newcells[i] = c->GetHandle();
									break;
								}
								case 21: //VTK_QUADRATIC_EDGE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_QUADRATIC_EDGE\n";
									break;
								}
								case 22: //VTK_QUADRATIC_TRIANGLE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_QUADRATIC_TRIANGLE\n";
									break;
								}
								case 23: //VTK_QUADRATIC_QUAD
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_QUADRATIC_QUAD\n";
									break;
								}
								case 24: //VTK_QUADRATIC_TETRA
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_QUADRATIC_TETRA\n";
									break;
								}
								case 25: //VTK_QUADRATIC_HEXAHEDRON
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_QUADRATIC_HEXAHEDRON\n";
									break;
								}
								case 26: //VTK_QUADRATIC_WEDGE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_QUADRATIC_WEDGE\n";
									break;
								}
								case 27: //VTK_QUADRATIC_PYRAMID
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_QUADRATIC_PYRAMID\n";
									break;
								}
								case 28: //VTK_BIQUADRATIC_QUAD
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_BIQUADRATIC_QUAD\n";
									break;
								}
								case 29: //VTK_TRIQUADRATIC_HEXAHEDRON
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_TRIQUADRATIC_HEXAHEDRON\n";
									break;
								}
								case 30: //VTK_QUADRATIC_LINEAR_QUAD
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_QUADRATIC_LINEAR_QUAD\n";
									break;
								}
								case 31: //VTK_QUADRATIC_LINEAR_WEDGE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_QUADRATIC_LINEAR_WEDGE\n";
									break;
								}
								case 32: //VTK_BIQUADRATIC_QUADRATIC_WEDGE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_BIQUADRATIC_QUADRATIC_WEDGE\n";
									break;
								}
								case 33: //VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON\n";
									break;
								}
								case 34: //VTK_BIQUADRATIC_TRIANGLE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_BIQUADRATIC_TRIANGLE\n";
									break;
								}
								case 35: //VTK_CUBIC_LINE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_CUBIC_LINE\n";
									break;
								}
								case 41: //VTK_CONVEX_POINT_SET
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no algorithm for VTK_CONVEX_POINT_SET yet\n";
									break;
								}
								case 42: //VTK_POLYHEDRON
								{
									int k = j;
									k++; //skip number of total number of integers
									std::vector<integer> sizes(cp[k++]);
									for(size_t m = 0; m < sizes.size(); m++)
									{
										sizes[m] = cp[k++];
										for (integer l = 0; l < sizes[m]; l++)
										{
											c_nodes.push_back(newnodes[cp[k++]]);
											c_nodes.back()->RemMarker(unused_marker);
										}
									}
									j = k;
									Cell c = CreateCell(c_nodes,sizes.data(),static_cast<integer>(sizes.size())/*, c_nodes*/).first;
									newcells[i] = c->GetHandle();
									break;
								}
								case 51: //VTK_PARAMETRIC_CURVE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_PARAMETRIC_CURVE\n";
									break;
								}
								case 52: //VTK_PARAMETRIC_SURFACE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_PARAMETRIC_SURFACE\n";
									break;
								}
								case 53: //VTK_PARAMETRIC_TRI_SURFACE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_PARAMETRIC_TRI_SURFACE\n";
									break;
								}
								case 54: //VTK_PARAMETRIC_QUAD_SURFACE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_PARAMETRIC_QUAD_SURFACE\n";
									break;
								}
								case 55: //VTK_PARAMETRIC_TETRA_REGION
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_PARAMETRIC_TETRA_REGION\n";
									break;
								}
								case 56: //VTK_PARAMETRIC_HEX_REGION
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_PARAMETRIC_HEX_REGION\n";
									break;
								}
								case 60: //VTK_HIGHER_ORDER_EDGE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_HIGHER_ORDER_EDGE\n";
									break;
								}
								case 61: //VTK_HIGHER_ORDER_TRIANGLE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_HIGHER_ORDER_TRIANGLE\n";
									break;
								}
								case 62: //VTK_HIGHER_ORDER_QUAD
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_HIGHER_ORDER_QUAD\n";
									break;
								}
								case 63: //VTK_HIGHER_ORDER_POLYGON
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_HIGHER_ORDER_POLYGON\n";
									break;
								}
								case 64: //VTK_HIGHER_ORDER_TETRAHEDRON
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_HIGHER_ORDER_TETRAHEDRON\n";
									break;
								}
								case 65: //VTK_HIGHER_ORDER_WEDGE
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_HIGHER_ORDER_WEDGE\n";
									break;
								}
								case 66: //VTK_HIGHER_ORDER_PYRAMID
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_HIGHER_ORDER_PYRAMID\n";
									break;
								}
								case 67: //VTK_HIGHER_ORDER_HEXAHEDRON
								{
									std::cout << __FILE__ << ":" << __LINE__ << " no order description for VTK_HIGHER_ORDER_HEXAHEDRON\n";
									break;
								}
								default: //VTK_NUMBER_OF_CELL_TYPES
								{
									std::cout << __FILE__ << ":" << __LINE__ << " Cell type " << ct[i] << " is not known\n";
								}
							}
							if( verbosity > 1 && i%report_pace == 0)
							{
								std::ios save(NULL);
								save.copyfmt(std::cout);
								std::cout << "cells " << std::setw(5) << std::fixed << std::setprecision(1) << (i*100.0)/(1.0*ncells) << "% ";
								std::cout << "cells " << std::setw(8) << NumberOfCells() << " faces " << std::setw(8) << NumberOfFaces() << " edges " << std::setw(8) << NumberOfEdges() << "\r" << std::flush;
								std::cout.copyfmt(save);
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
					if( verbosity > 1 )
					{
						std::cout << "\nCell types: ";
						for(ElementType etype = EDGE; etype <= MESH; etype = NextElementType(etype) ) if( etype & read_into_cell )
							std::cout << ElementTypeName(etype) << " ";
						std::cout << std::endl;
					}
					break;
				}
				case R_SPOINTS:
				{
					std::cout << "STRUCTURED POINTS is not yet implemented\n";
					break;
				}
				case R_POLYDATA:
				{
					std::cout << "POLYDATA is not yet implemented\n";
					break;
				}
				case R_SGRID:
				{
					std::cout << "STRUCTURED GRID is not yet implemented\n";
					break;
				}
				case R_RGRID:
				{
					std::cout << "RECTLINEAR GRID is not yet implemented\n";
					break;
				}
				case R_FIELD:
				{
					std::cout << "FIELD is not yet implemented\n";
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
		//this is a system tag that may have bad values leading to crash in parallel
		if (HaveTag("GLOBAL_ID")) DeleteTag(GetTag("GLOBAL_ID"));
	}
}

#endif
