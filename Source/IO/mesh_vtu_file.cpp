#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif


#include "inmost.h"
#include <cfloat>

#if defined(USE_MESH)

//static int __isnan__(double x) { return x != x; }
//static int isinf(double x) { return !isnan(x) && isnan(x - x); }
//static int __isinf__(double x) { return fabs(x) > DBL_MAX; }
//static int __isbad(double x) { return __isnan__(x) || __isinf__(x); }



namespace INMOST
{
	
	int VtkElementType(ElementType t); //from mesh_vtk_file.cpp
	
	static int getdigits(int val)
	{
		int digits = 0;
		do
		{
			digits++;
			val /= 10;
		}
		while(val);
		return digits;
	}

	
	void Mesh::SavePVTU(std::string file)
	{
		std::string name=file;
		std::string::size_type pos=name.rfind(".pvtu");
		name.erase(pos); 

		std::string::size_type l=name.find_last_of("/\\");
		std::string fname=name.substr(l+1,name.length());
		
		
		
		if( GetProcessorRank() == 0 )
		{
			std::set< std::string > nosave, saveonly;
			nosave = TagOptions("nosave");
			saveonly = TagOptions("saveonly");
			std::ofstream fh(file.c_str());
			
			std::map<DataType,std::string> type_name;
			type_name[DATA_INTEGER] = "Int32";
			type_name[DATA_BULK] = "Int8";
			type_name[DATA_REAL] = "Float64";
#if defined(USE_AUTODIFF)
			type_name[DATA_VARIABLE] = "Float64";
#endif
		
			fh << "<VTKFile type=\"PUnstructuredGrid\">" << std::endl;
			fh << "\t<PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
			fh << "\t\t<PPointData>" << std::endl;
			for(Mesh::iteratorTag it = BeginTag(); it != EndTag(); it++)
			{
				//SKIPHERE
				if( it->GetSize() == ENUMUNDEF ) continue;
				if( it->GetDataType() == DATA_REFERENCE ) continue;
				if( it->GetDataType() == DATA_REMOTE_REFERENCE ) continue;
				if( *it == MarkersTag() ) continue; 
				if( *it == HighConnTag() ) continue;
				if( *it == LowConnTag() ) continue;
				if( *it == CoordsTag() ) continue;
				if( *it == SetNameTag() ) continue;
				if( !it->isDefined(NODE) ) continue;
				if( it->GetTagName().substr(0,9) == "PROTECTED" ) continue;
				if( CheckSaveSkip(it->GetTagName(),nosave,saveonly) )
					continue;
				fh << "\t\t\t<PDataArray";
				fh << " type=\"" << type_name[it->GetDataType()] << "\"";
				fh << " Name=\"" << it->GetTagName() << "\"";
				fh << " Format=\"ascii\"";
				fh << " NumberOfComponents=\"" << it->GetSize() << "\"";
				fh << "/>" << std::endl;
			}
			fh << "\t\t</PPointData>" << std::endl;
			fh << "\t\t<PCellData>" << std::endl;
			for(Mesh::iteratorTag it = BeginTag(); it != EndTag(); it++)
			{
				//SKIPHERE
				if( it->GetSize() == ENUMUNDEF ) continue;
				if( it->GetDataType() == DATA_REFERENCE ) continue;
				if( it->GetDataType() == DATA_REMOTE_REFERENCE ) continue;
				if( *it == MarkersTag() ) continue; 
				if( *it == HighConnTag() ) continue;
				if( *it == LowConnTag() ) continue;
				if( *it == CoordsTag() ) continue;
				if( *it == SetNameTag() ) continue;
				if( !it->isDefined(CELL) ) continue;
				if( it->GetTagName().substr(0,9) == "PROTECTED" ) continue;
				if( CheckSaveSkip(it->GetTagName(),nosave,saveonly) )
					continue;
				fh << "\t\t\t<PDataArray";
				fh << " type=\"" << type_name[it->GetDataType()] << "\"";
				fh << " Name=\"" << it->GetTagName() << "\"";
				fh << " Format=\"ascii\"";
				fh << " NumberOfComponents=\"" << it->GetSize() << "\"";
				fh << "/>" << std::endl;
			}
			fh << "\t\t</PCellData>" << std::endl;
			fh << "\t\t<PPoints>" << std::endl;
			fh << "\t\t\t<PDataArray";
			fh << " type=\"Float64\"";
			fh << " NumberOfComponents=\"" << GetDimensions() << "\"";
			fh << " Format=\"ascii\"";
			fh << "/>" << std::endl;
			fh << "\t\t</PPoints>" << std::endl;
			for(int k = 0; k < GetProcessorsNumber(); ++k)
			{
				std::stringstream filename;
				filename << fname << 'p';
				for(int q = getdigits(k+1); q < getdigits(GetProcessorsNumber()); ++q) 
					filename << '0'; //leading zeros
				filename << (k+1);
				filename << ".vtu";
				fh << "\t\t<Piece Source=\"" << filename.str() << "\"/>" << std::endl; 
			}
			fh << "\t</PUnstructuredGrid>" << std::endl;
			fh << "</VTKFile>" << std::endl;
		}
		
		
		{
			std::stringstream filename;
			filename << name << 'p';
			for(int q = getdigits(GetProcessorRank()+1); q < getdigits(GetProcessorsNumber()); ++q) 
				filename << '0'; //leading zeros
			filename << (GetProcessorRank()+1);
			filename << ".vtu";
			Save(filename.str());
		}
	}
	
	void Mesh::LoadPVTU(std::string file)
	{
		std::vector<std::string> files;
		std::string path = "";
		
		size_t l = file.find_last_of("/\\");
		if( l != std::string::npos )
			path = file.substr(0,l+1);
			
		
		std::fstream f(file.c_str(), std::ios::in);
		XMLReader r(file, f);
		XMLReader::XMLTree t = r.ReadXML();

		if (t.GetName() == "VTKFile")
		{
			if(t.GetAttrib("type") != "PUnstructuredGrid")
			{
				std::cout << "Expected parallel unstructured grid type inside " << file << std::endl;
				throw BadFile;
			}
			//const XMLReader::XMLTree * v = t.GetChild("PUnstructuredGrid")->GetChild("Piece");
			const XMLReader::XMLTree * da = t.GetChild("PUnstructuredGrid");
			for (int k = 0; k < da->NumChildren(); ++k)
			{
				const XMLReader::XMLTree * pd = &da->GetChild(k);
				if (pd->GetName() == "Piece")
				{
					int nca = pd->FindAttrib("Source");
					if (nca != pd->NumAttrib()) 
					{
						std::string filename = pd->GetAttrib(nca).value;
						files.push_back(filename);
					}
				}
				//else if( GetProcessorRank() == 0 ) std::cout << __FILE__ << ":" << __LINE__ << "I don't use yet " << pd->GetName() << " in PUnstructuredGrid" << std::endl;
			}
		}
		for(int i = GetProcessorRank(); i < static_cast<int>(files.size()); i+= GetProcessorsNumber()) 
		{
			//std::cout << "load " << i << ": " << path + files[i] << " by " << GetProcessorRank() << std::endl;
			Load(path + files[i]);
		}
		ResolveShared();
	}
	
	void Mesh::SaveVTU(std::string File)
	{
		std::set< std::string > nosave, saveonly;
		nosave = TagOptions("nosave");
		saveonly = TagOptions("saveonly");
		std::ofstream f(File.c_str());
		std::map<DataType,std::string> type_name;
		type_name[DATA_INTEGER] = "Int32";
		type_name[DATA_BULK] = "Int8";
		type_name[DATA_REAL] = "Float64";
#if defined(USE_AUTODIFF)
		type_name[DATA_VARIABLE] = "Float64";
#endif
		std::map<DataType,std::string> type_undef;
		{
			std::stringstream s;
			s << INT_MIN;
			type_undef[DATA_INTEGER] = s.str();
		}
		type_undef[DATA_BULK] = "255";
		type_undef[DATA_REAL] = "-0.9999E30";
#if defined(USE_AUTODIFF)
		type_undef[DATA_VARIABLE] = "-0.9999E30";
#endif
		bool keep_ghost = false;
		for (INMOST_DATA_ENUM_TYPE k = 0; k < file_options.size(); ++k)
		{
			if (file_options[k].first == "KEEP_GHOST")
			{
				keep_ghost = true;
			}
		}
		int num_cells = NumberOfCells(), num_nodes = NumberOfNodes();
		//give id to nodes
		TagInteger nid = CreateTag("PROTECTED_TEMPORARY_NODE_ID", DATA_INTEGER, NODE, NONE, 1);
		if (!keep_ghost)
		{
			MarkerType used = CreateMarker();
			num_cells = 0;
			for(iteratorCell it = BeginCell(); it != EndCell(); ++it)
				if (it->GetStatus() != Element::Ghost)
				{
					num_cells++;
					it->getAdjElements(NODE).SetMarker(used);
				}
			num_nodes = 0;
			for (Mesh::iteratorNode jt = BeginNode(); jt != EndNode(); ++jt)
				if (jt->GetMarker(used)) nid[*jt] = num_nodes++; else nid[*jt] = -1;
			ReleaseMarker(used, NODE);
		}
		else
		{
			INMOST_DATA_INTEGER_TYPE q = 0;
			for (Mesh::iteratorNode jt = BeginNode(); jt != EndNode(); ++jt)
				nid[*jt] = q++;
		}

		
		f << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
		f << "\t<UnstructuredGrid>" << std::endl;
		//~ f << "\t\t<FieldData>" << std::endl;
		//~ f << "\t\t\t<DataArray>" << std::endl;
		//~ f << "\t\t\t</DataArray>" << std::endl;
		//~ f << "\t\t</FieldData>" << std::endl;
		f << "\t\t<Piece NumberOfPoints=\"" << num_nodes << "\" NumberOfCells=\"" << num_cells << "\">" << std::endl;
		f << "\t\t\t<PointData>" << std::endl;
		for(Mesh::iteratorTag it = BeginTag(); it != EndTag(); it++)
		{
			//SKIPHERE
			if( it->GetSize() == ENUMUNDEF ) continue;
			if( it->GetDataType() == DATA_REFERENCE ) continue;
			if( it->GetDataType() == DATA_REMOTE_REFERENCE ) continue;
			if( *it == MarkersTag() ) continue; 
			if( *it == HighConnTag() ) continue;
			if( *it == LowConnTag() ) continue;
			if( *it == CoordsTag() ) continue;
			if( *it == SetNameTag() ) continue;
			if( !it->isDefined(NODE) ) continue;
			if( it->GetTagName().substr(0,9) == "PROTECTED" ) continue;
			if( CheckSaveSkip(it->GetTagName(),nosave,saveonly) )
				continue;
			f << "\t\t\t\t<DataArray type=\"" << type_name[it->GetDataType()] << "\" Name=\"" << it->GetTagName() << "\" NumberOfComponents=\"" << it->GetSize() << "\" format=\"ascii\">" << std::endl;
			for(Mesh::iteratorNode jt = BeginNode(); jt != EndNode(); ++jt) if( nid[*jt] != -1 )
			{
				if( !jt->HaveData(*it) )
				{
					for(INMOST_DATA_ENUM_TYPE k = 0; k < it->GetSize(); ++k)
						f << type_undef[it->GetDataType()] << " ";
				}
				else if( it->GetDataType() == DATA_REAL )
				{
					for (INMOST_DATA_ENUM_TYPE k = 0; k < it->GetSize(); ++k)
					{
						INMOST_DATA_REAL_TYPE val = jt->RealArray(*it)[k];
						if (__isbad(val))
							f << type_undef[it->GetDataType()] << " ";
						else f << val << " ";
					}
				}
				else if( it->GetDataType() == DATA_INTEGER )
				{
					for(INMOST_DATA_ENUM_TYPE k = 0; k < it->GetSize(); ++k)
						f << jt->IntegerArray(*it)[k] << " ";
				}
				else if( it->GetDataType() == DATA_BULK )
				{
					for(INMOST_DATA_ENUM_TYPE k = 0; k < it->GetSize(); ++k)
						f << (int)jt->BulkArray(*it)[k] << " ";
				}
#if defined(USE_AUTODIFF)
				else if( it->GetDataType() == DATA_VARIABLE )
				{
					for (INMOST_DATA_ENUM_TYPE k = 0; k < it->GetSize(); ++k)
					{
						INMOST_DATA_REAL_TYPE val = jt->VariableArray(*it)[k].GetValue();
						if (__isbad(val))
							f << type_undef[it->GetDataType()] << " ";
						else f << val << " ";
					}
				}
#endif
				f << std::endl;
			}
			f << "\t\t\t\t</DataArray>" << std::endl;
		}
		f << "\t\t\t</PointData>" << std::endl;
		f << "\t\t\t<CellData>" << std::endl;
		for(Mesh::iteratorTag it = BeginTag(); it != EndTag(); it++)
		{
			//SKIPHERE
			if( it->GetSize() == ENUMUNDEF ) continue;
			if( it->GetDataType() == DATA_REFERENCE ) continue;
			if( it->GetDataType() == DATA_REMOTE_REFERENCE ) continue;
			if( *it == MarkersTag() ) continue; 
			if( *it == HighConnTag() ) continue;
			if( *it == LowConnTag() ) continue;
			if( *it == CoordsTag() ) continue;
			if( *it == SetNameTag() ) continue;
			if( !it->isDefined(CELL) ) continue;
			if( it->GetTagName().substr(0,9) == "PROTECTED" ) continue;
			if( CheckSaveSkip(it->GetTagName(),nosave,saveonly) )
				continue;
			f << "\t\t\t\t<DataArray type=\"" << type_name[it->GetDataType()] << "\" Name=\"" << it->GetTagName() << "\" NumberOfComponents=\"" << it->GetSize() << "\" format=\"ascii\">" << std::endl;
			for(Mesh::iteratorCell jt = BeginCell(); jt != EndCell(); ++jt) if( keep_ghost || jt->GetStatus() != Element::Ghost )
			{
				if( !jt->HaveData(*it) )
				{
					for(INMOST_DATA_ENUM_TYPE k = 0; k < it->GetSize(); ++k)
						f << type_undef[it->GetDataType()] << " ";
				}
				else if( it->GetDataType() == DATA_REAL )
				{
					for(INMOST_DATA_ENUM_TYPE k = 0; k < it->GetSize(); ++k)
					{
						INMOST_DATA_REAL_TYPE val = jt->RealArray(*it)[k];
						if (__isbad(val))
							f << type_undef[it->GetDataType()] << " ";
						else f << val << " ";
					}
				}
				else if( it->GetDataType() == DATA_INTEGER )
				{
					for(INMOST_DATA_ENUM_TYPE k = 0; k < it->GetSize(); ++k)
						f << jt->IntegerArray(*it)[k] << " ";
				}
				else if( it->GetDataType() == DATA_BULK )
				{
					for(INMOST_DATA_ENUM_TYPE k = 0; k < it->GetSize(); ++k)
						f << (int)jt->BulkArray(*it)[k] << " ";
				}
#if defined(USE_AUTODIFF)
				else if( it->GetDataType() == DATA_VARIABLE )
				{
					for(INMOST_DATA_ENUM_TYPE k = 0; k < it->GetSize(); ++k)
					{
						INMOST_DATA_REAL_TYPE val = jt->VariableArray(*it)[k].GetValue();
						if (__isbad(val))
							f << type_undef[it->GetDataType()] << " ";
						else f << val << " ";
					}
				}
#endif
				f << std::endl;
			}
			f << "\t\t\t\t</DataArray>" << std::endl;
		}
		f << "\t\t\t</CellData>" << std::endl;
		f << "\t\t\t<Points>" << std::endl;
		f << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"" << GetDimensions() << "\" format=\"ascii\">" << std::endl;
		for(Mesh::iteratorNode jt = BeginNode(); jt != EndNode(); ++jt) if( nid[*jt] != -1 )
		{
			for(int k = 0; k < GetDimensions(); ++k)
				f << jt->Coords()[k] << " ";
			f << std::endl;
		}
		f << "\t\t\t\t</DataArray>" << std::endl;
		f << "\t\t\t</Points>" << std::endl;
		f << "\t\t\t<Cells>" << std::endl;
		
		f << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
		for(Mesh::iteratorCell jt = BeginCell(); jt != EndCell(); ++jt) if( keep_ghost || jt->GetStatus() != Element::Ghost )
		{
			ElementArray<Node> nodes(this); //= jt->getNodes();
			RestoreCellNodes(*jt,nodes);
			assert(nodes.size() == jt->nbAdjElements(NODE));
			for(ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); ++kt)
				f << nid[*kt] << " ";
			f << std::endl;
		}
		f << "\t\t\t\t</DataArray>" << std::endl;
		f << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
		{
			size_t offset = 0;
			for(Mesh::iteratorCell jt = BeginCell(); jt != EndCell(); ++jt)  if (keep_ghost || jt->GetStatus() != Element::Ghost)
			{
				//ElementArray<Node> nodes = jt->getNodes();
				//offset += nodes.size();
				offset += jt->nbAdjElements(NODE);
				f << offset << std::endl;
			}
		}
		f << "\t\t\t\t</DataArray>" << std::endl;
		bool need_faces = false;
		f << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
		for(Mesh::iteratorCell jt = BeginCell(); jt != EndCell(); ++jt) if (keep_ghost || jt->GetStatus() != Element::Ghost)
		{
			int etype = VtkElementType(jt->GetGeometricType());
			f << etype << std::endl;
			if( etype == 42 ) need_faces = true;
		}
		f << "\t\t\t\t</DataArray>" << std::endl;
		if( need_faces )
		{
			f << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">" << std::endl;	
			for(Mesh::iteratorCell jt = BeginCell(); jt != EndCell(); ++jt)  if (keep_ghost || jt->GetStatus() != Element::Ghost)
			{
				int etype = VtkElementType(jt->GetGeometricType());
				if( etype == 42 ) //polyhedron
				{
					ElementArray<Face> faces = jt->getFaces();
					f << faces.size() << " ";
					for(ElementArray<Face>::iterator kt = faces.begin(); kt != faces.end(); ++kt)
					{
						ElementArray<Node> nodes = kt->getNodes();
						assert(nodes.size() == kt->nbAdjElements(NODE));
						f << nodes.size() << " ";
						for(ElementArray<Node>::iterator qt = nodes.begin(); qt != nodes.end(); ++qt)
							f << nid[*qt] << " ";
					}
					f << std::endl;
				}
				//else f << -1 << std::endl;
			}
			f << "\t\t\t\t</DataArray>" << std::endl;
			f << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;	
			{
				size_t offset = 0;
				for(Mesh::iteratorCell jt = BeginCell(); jt != EndCell(); ++jt) if (keep_ghost || jt->GetStatus() != Element::Ghost)
				{
					int etype = VtkElementType(jt->GetGeometricType());
					if( etype == 42 ) //polyhedron
					{
						offset++; //number of faces
						ElementArray<Face> faces = jt->getFaces();
						for(ElementArray<Face>::iterator kt = faces.begin(); kt != faces.end(); ++kt)
						{
							offset++; //number of nodes in face
							//offset+= kt->getNodes().size(); //node ids
							offset += kt->nbAdjElements(NODE);
						}
						f << offset << std::endl;
					}
					else f << -1 << std::endl;
				}
			}
			f << "\t\t\t\t</DataArray>" << std::endl;
		}
		f << "\t\t\t</Cells>" << std::endl;
		f << "\t\t</Piece>" << std::endl;
		f << "\t</UnstructuredGrid>" << std::endl;
		f << "</VTKFile>" << std::endl;
		DeleteTag(nid);
	}
	
	void Mesh::LoadVTU(std::string File)
	{
		std::set< std::string > noload, loadonly;		
		noload = TagOptions("noload");
		loadonly = TagOptions("loadonly");
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

		std::vector<HandleType> newnodes;
		std::vector<HandleType> newpolyh;
		std::vector<HandleType> newcells;
		

		std::fstream f(File.c_str(), std::ios::in);
		XMLReader r(File, f);
		XMLReader::XMLTree t = r.ReadXML();

		if (t.GetName() == "VTKFile")
		{
			if(t.GetAttrib("type") != "UnstructuredGrid")
			{
				std::cout << "Expected unstructured grid type inside " << File << std::endl;
				throw BadFile;
			}
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
				int cconn, coffset = 0, totread = 0, nfaces, nfacenodes;
				ElementArray<Node> hnodes(this);
				ElementArray<Face> hfaces(this);
				while (!faceoffsets.eof())
				{
					faceoffsets >> coffset;
					if( coffset == -1 ) continue;
					//nread = coffset - totread;
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
			
			

			if (verbosity > 0 && GetProcessorRank() == 0)
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
							//Cell c = CreateCell(c_faces).first;
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
							//Cell c = CreateCell(c_faces).first;
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
							//Cell c = CreateCell(c_faces).first;
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
							//Cell c = CreateCell(c_faces).first;
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
						newcells[q] = CreateCell(hnodes, nodesnum, sizes, 4, hnodes).first.GetHandle();
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
						newcells[q] = CreateCell(hnodes, nodesnum, sizes, 6, hnodes).first.GetHandle();
					}
					else if (ctype == 13) // VTK_WEDGE
					{
						assert(nread == 6);
						const integer nodesnum[18] = { 0, 2, 5, 3, 1, 4, 5, 2, 0, 3, 4, 1, 3, 5, 4, 0, 1, 2 };
						//const integer nodesnum[18] = { 0, 3, 5, 2, 0, 1, 4, 3, 1, 4, 5, 2, 3, 4, 5, 0, 2, 1 };
						const integer sizes[5] = { 4, 4, 4, 3, 3 };
						newcells[q] = CreateCell(hnodes, nodesnum, sizes, 5, hnodes).first.GetHandle();
					}
					else if (ctype == 14) //VTK_PYRAMID
					{
						assert(nread == 5);
						const integer nodesnum[16] = { 0, 4, 3, 0, 1, 4, 1, 2, 4, 3, 4, 2, 0, 3, 2, 1 };
						const integer sizes[5] = { 3, 3, 3, 3, 4 };
						newcells[q] = CreateCell(hnodes, nodesnum, sizes, 5, hnodes).first.GetHandle();
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
						newcells[q] = CreateCell(hnodes, nodesnum, sizes, 7, hnodes).first.GetHandle();
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
						newcells[q] = CreateCell(hnodes, nodesnum, sizes, 8, hnodes).first.GetHandle();
					}
					else if (ctype == 42) //VTK_POLYHEDRON
						newcells[q] = newpolyh[npolyh++];
					else std::cout << __FILE__ << ":" << __LINE__ << " Implement VTK format " << ctype << std::endl;
				}
			}
			//read data
			{
				std::string dname[2] = { "PointData", "CellData" };
				ElementType etype[2] = { NODE, CELL };
				ElementType dsparse[2] = { NONE, NONE };
				if (have_faces)
				{
					etype[1] |= FACE;
					dsparse[1] |= FACE;
				}
				if (have_edges)
				{
					etype[1] |= EDGE;
					dsparse[1] |= EDGE;
				}
				if (have_nodes)
				{
					etype[1] |= NODE;
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
								DataType dtype = DATA_REAL;
								int nca = pd->FindAttrib("NumberOfComponents");
								if (nca != pd->NumAttrib()) ncomps = atoi(pd->GetAttrib(nca).value.c_str());
								int nct = pd->FindAttrib("type");
								if( nct != pd->NumAttrib())
								{
									std::string t = pd->GetAttrib(nct).value;
									if( t == "UInt8" || t == "Int8" ) dtype = DATA_BULK;
									else if( t == "UInt16" || t == "Int16" ) dtype = DATA_INTEGER;
									else if( t == "UInt32" || t == "Int32" ) dtype = DATA_INTEGER;
									else if( t == "UInt64" || t == "Int64" ) dtype = DATA_INTEGER;
								}
								if( !CheckLoadSkip(pd->GetAttrib("Name"),noload,loadonly) )
								{
									if( dtype == DATA_REAL )
									{
										TagRealArray t = CreateTag(pd->GetAttrib("Name"), dtype, etype[j], dsparse[j], ncomps);
										std::stringstream inp(pd->GetContents());
										for (int l = 0; l < dsize[j]; ++l)
										{
											for (INMOST_DATA_ENUM_TYPE q = 0; q < t.GetSize(); ++q)
												inp >> t[darray[j][l]][q];
										}
									}
									else if( dtype == DATA_INTEGER )
									{
										TagIntegerArray t = CreateTag(pd->GetAttrib("Name"), dtype, etype[j], dsparse[j], ncomps);
										std::stringstream inp(pd->GetContents());
										for (int l = 0; l < dsize[j]; ++l)
										{
											for (INMOST_DATA_ENUM_TYPE q = 0; q < t.GetSize(); ++q)
												inp >> t[darray[j][l]][q];
										}
									}
									else if( dtype == DATA_BULK )
									{
										TagBulkArray t = CreateTag(pd->GetAttrib("Name"), dtype, etype[j], dsparse[j], ncomps);
										std::stringstream inp(pd->GetContents());
										for (int l = 0; l < dsize[j]; ++l)
										{
											for (INMOST_DATA_ENUM_TYPE q = 0; q < t.GetSize(); ++q)
												inp >> t[darray[j][l]][q];
										}
									}
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
		//this is a system tag that may have bad values leading to crash in parallel
		if (HaveTag("GLOBAL_ID")) DeleteTag(GetTag("GLOBAL_ID"));
	}
}

#endif
