
#include "inmost.h"
#include "io.hpp"

#if defined(USE_MESH)

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


namespace INMOST
{
	typedef char HeaderType;
	//const HeaderType EndOfData  = 0x01;
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
	
	void Mesh::SavePMF(std::string File)
	{
		io_converter<INMOST_DATA_INTEGER_TYPE ,INMOST_DATA_REAL_TYPE> iconv;
		io_converter<INMOST_DATA_ENUM_TYPE    ,INMOST_DATA_REAL_TYPE> uconv;
		io_converter<INMOST_DATA_BIG_ENUM_TYPE,INMOST_DATA_REAL_TYPE> buconv;
		INMOST_DATA_ENUM_TYPE nlow,nhigh, lid;
		char wetype;
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
		
		
		std::set< std::string > nosave, saveonly, noderivs;
		
		nosave = TagOptions("nosave");
		saveonly = TagOptions("saveonly");
		noderivs = TagOptions("noderivatives");
		
		INMOST_DATA_ENUM_TYPE ntags = 0;
		for(Mesh::iteratorTag it = BeginTag(); it != EndTag(); it++)
		{
			//SKIPHERE
			//don't forget to change header[6] if you skip more
			//should match with content after MeshDataHeader
			if( *it == MarkersTag() ) continue; //temporary fix to protect markers being loaded into hide_marker and then marked elements being destroyed by EndModification
			if( *it == HighConnTag() ) continue;
			if( *it == LowConnTag() ) continue;
			if( *it == CoordsTag() ) continue;
			if( *it == SetNameTag() ) continue;
			
			if( CheckSaveSkip(it->GetTagName(),nosave,saveonly) )
				continue;
			
			ntags++;
			
		}
		
		
		INMOST_DATA_ENUM_TYPE header[9] =
		{
			static_cast<INMOST_DATA_ENUM_TYPE>(GetDimensions()),
			static_cast<INMOST_DATA_ENUM_TYPE>(NumberOfNodes()),
			static_cast<INMOST_DATA_ENUM_TYPE>(NumberOfEdges()),
			static_cast<INMOST_DATA_ENUM_TYPE>(NumberOfFaces()),
			static_cast<INMOST_DATA_ENUM_TYPE>(NumberOfCells()),
			static_cast<INMOST_DATA_ENUM_TYPE>(NumberOfSets()),
			//NumberOfTags()-5, //add counter to skip unwanted tags here SKIPHERE, search by SKIPHERE for additional instructions
			ntags,
			static_cast<INMOST_DATA_ENUM_TYPE>(m_state),
			static_cast<INMOST_DATA_ENUM_TYPE>(GetProcessorRank())
		},k;
		for(k = 0; k < 9; k++) uconv.write_iValue(out,header[k]);
		out.write(reinterpret_cast<char *>(remember),sizeof(remember));
		
		
		
		// Tags
		out << INMOST::TagsHeader;
		uconv.write_iValue(out,header[6]);
		REPORT_STR("TagsHeader");
		REPORT_VAL("tag_size",header[6]);
		INMOST_DATA_ENUM_TYPE tags_written = 0;
		for(Mesh::iteratorTag it = BeginTag(); it != EndTag(); it++)
		{
			//SKIPHERE
			//don't forget to change header[6] if you skip more
			//should match with content after MeshDataHeader
			if( *it == MarkersTag() ) continue; //temporary fix to protect markers being loaded into hide_marker and then marked elements being destroyed by EndModification
			if( *it == HighConnTag() ) continue;
			if( *it == LowConnTag() ) continue;
			if( *it == CoordsTag() ) continue;
			if( *it == SetNameTag() ) continue;
			
			
			if( CheckSaveSkip(it->GetTagName(),nosave,saveonly) )
				continue;
			
			
			std::string name = it->GetTagName();
			INMOST_DATA_ENUM_TYPE namesize = static_cast<INMOST_DATA_BULK_TYPE>(name.size());
			INMOST_DATA_BULK_TYPE datatype = static_cast<INMOST_DATA_BULK_TYPE>(it->GetDataType());
			ElementType sparsemask = NONE;
			ElementType definedmask = NONE;
			INMOST_DATA_ENUM_TYPE datalength = it->GetSize();
			for(ElementType current_type = NODE; current_type <= MESH; current_type = current_type << 1)
			{
				if( it->isSparse (current_type) ) sparsemask  |= current_type;
				if( it->isDefined(current_type) ) definedmask |= current_type;
			}
			uconv.write_iValue(out,namesize);
			out.write(name.c_str(), name.size());
			out.put(datatype);
			out.put(sparsemask);
			out.put(definedmask);
			uconv.write_iValue(out,datalength);
			++tags_written;
		}
		assert(tags_written == header[6]);
		
		
		Tag set_id = CreateTag("TEMPORARY_ELEMENT_ID_PMF_WRITER",DATA_INTEGER,ESET|CELL|FACE|EDGE|NODE,NONE,1);
		{
			Storage::integer cur_num = 0;
			for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it) it->IntegerDF(set_id) = cur_num++;
			cur_num = 0;
			for(Mesh::iteratorEdge it = BeginEdge(); it != EndEdge(); ++it) it->IntegerDF(set_id) = cur_num++;
			cur_num = 0;
			for(Mesh::iteratorFace it = BeginFace(); it != EndFace(); ++it) it->IntegerDF(set_id) = cur_num++;
			cur_num = 0;
			for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); ++it) it->IntegerDF(set_id) = cur_num++;
			cur_num = 0;
			for(Mesh::iteratorSet it = BeginSet(); it != EndSet(); ++it) it->IntegerDF(set_id) = cur_num++;
		}
		
		
		// Nodes
		out << INMOST::NodeHeader;
		uconv.write_iValue(out,header[1]);
		REPORT_STR("NodeHeader");
		REPORT_VAL("node_size",header[1]);
		for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); it++)
		{
			Storage::real_array coords = it->Coords();
			for(Storage::real_array::size_type it = 0; it < coords.size(); it++)
				iconv.write_fValue(out,coords[it]);
		}
		
		// Edges
		out << INMOST::EdgeHeader;
		uconv.write_iValue(out,header[2]);
		REPORT_STR("EdgeHeader");
		REPORT_VAL("edge_size",header[2]);
		
		for(Mesh::iteratorEdge it = BeginEdge(); it != EndEdge(); it++)
		{
			Element::adj_type & lc = LowConn(*it);
			nlow = static_cast<INMOST_DATA_ENUM_TYPE>(lc.size());
			uconv.write_iValue(out,nlow);
			for(Element::adj_type::size_type kt = 0; kt < lc.size(); ++kt)
			{
				lid = IntegerDF(lc[kt],set_id);
				uconv.write_iValue(out,lid);
			}
		}
		
		// Faces
		out << INMOST::FaceHeader;
		uconv.write_iValue(out,header[3]);
		REPORT_STR("FaceHeader");
		REPORT_VAL("face_size",header[3]);
		for(Mesh::iteratorFace it = BeginFace(); it != EndFace(); it++)
		{
			Element::adj_type & lc = LowConn(*it);
			nlow = static_cast<INMOST_DATA_ENUM_TYPE>(lc.size());
			uconv.write_iValue(out,nlow);
			for(Element::adj_type::size_type kt = 0; kt < lc.size(); ++kt)
			{
				lid = IntegerDF(lc[kt],set_id);
				uconv.write_iValue(out,lid);
			}
		}
		
		// Cells
		out << INMOST::CellHeader;
		uconv.write_iValue(out,header[4]);
		REPORT_STR("CellHeader");
		REPORT_VAL("cell_size",header[4]);
		for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); it++)
		{
			Element::adj_type & lc = LowConn(*it);
			nlow = static_cast<INMOST_DATA_ENUM_TYPE>(lc.size());
			uconv.write_iValue(out,nlow);
			for(Element::adj_type::size_type kt = 0; kt < lc.size(); ++kt)
			{
				lid = IntegerDF(lc[kt],set_id);
				uconv.write_iValue(out,lid);
			}
			Element::adj_type & hc = HighConn(*it);
			nhigh = static_cast<INMOST_DATA_ENUM_TYPE>(hc.size());
			uconv.write_iValue(out,nhigh);
			for(Element::adj_type::size_type kt = 0; kt < hc.size(); ++kt)
			{
				lid = IntegerDF(hc[kt],set_id);
				uconv.write_iValue(out,lid);
			}
		}
		
		// Element Sets
		out << INMOST::ESetHeader;
		REPORT_STR("ESetHeader");
		REPORT_VAL("eset_size",header[5]);
		uconv.write_iValue(out,header[5]);
		for(Mesh::iteratorSet it = BeginSet(); it != EndSet(); ++it)
		{
			std::string name = it->GetName();
			INMOST_DATA_ENUM_TYPE name_size = static_cast<INMOST_DATA_ENUM_TYPE>(name.size());
			assert(name_size < 4096);
			uconv.write_iValue(out,name_size);
			out.write(name.c_str(),name.size());
			Element::adj_type & lc = LowConn(it->GetHandle());
			uconv.write_iValue(out,static_cast<enumerator>(lc.size()));
			for(Element::adj_type::iterator kt = lc.begin(); kt != lc.end(); ++kt)
			{
				if( *kt != InvalidHandle() )
				{
					wetype = GetHandleElementType(*kt);
					out.put(wetype);
					assert(wetype != NONE);
					lid = IntegerDF(*kt,set_id);
					uconv.write_iValue(out,lid);
				}
				else out.put(NONE);
			}
			Element::adj_type & hc = HighConn(it->GetHandle());
			//write tree information
			uconv.write_iValue(out,static_cast<enumerator>(hc.size()));
			for(Element::adj_type::iterator kt = hc.begin(); kt != hc.begin()+ElementSet::high_conn_reserved-1; ++kt)
			{
				if( *kt != InvalidHandle() )
				{
					wetype = GetHandleElementType(*kt);
					out.put(wetype);
					assert(wetype != NONE);
					lid = IntegerDF(*kt,set_id);
					uconv.write_iValue(out,lid);
				}
				else out.put(NONE);
			}
			//write additional information
			for(Element::adj_type::iterator kt = hc.begin()+ElementSet::high_conn_reserved-1; kt != hc.end(); ++kt)
			{
				uconv.write_iValue(out,*kt);
			}
		}
		
		out << INMOST::MeshDataHeader;
		REPORT_STR("MeshDataHeader");
		
		for(Mesh::iteratorTag jt = BeginTag(); jt != EndTag(); jt++)
		{
			std::string tagname = jt->GetTagName();
			//skipping should match with header[6] and content
			// after TagsHeader
			//SKIPHERE
			if( *jt == set_id ) continue;
			if( *jt == HighConnTag() ) continue;
			if( *jt == LowConnTag() ) continue;
			if( *jt == MarkersTag() ) continue;
			if( *jt == CoordsTag() ) continue;
			if( *jt == SetNameTag() ) continue;
			
			if( CheckSaveSkip(tagname,nosave,saveonly) )
				continue;
			
			
			REPORT_VAL("TagName",tagname);
			for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
				if( jt->isDefined(etype) )
				{
					INMOST_DATA_ENUM_TYPE q = 0;
					INMOST_DATA_ENUM_TYPE tagsize = jt->GetSize(), recsize = tagsize, lid, k;
					DataType data_type = jt->GetDataType();
					bool sparse = jt->isSparse(etype);
					for(Mesh::iteratorStorage it = Begin(etype); it != End(); it++)
					{
						if( !sparse || (sparse && it->HaveData(*jt)) )
						{
							if( sparse ) uconv.write_iValue(out,q);
							if( tagsize == ENUMUNDEF )
							{
								recsize = it->GetDataSize(*jt);
								uconv.write_iValue(out,recsize);
							}
							switch(data_type)
							{
								case DATA_REAL:
								{
									Storage::real_array arr = it->RealArray(*jt);
									for(k = 0; k < recsize; k++)
										uconv.write_fValue(out,arr[k]);
								} break;
								case DATA_INTEGER:
								{
									Storage::integer_array arr = it->IntegerArray(*jt);
									for(k = 0; k < recsize; k++)
										uconv.write_iValue(out,arr[k]);
								} break;
								case DATA_BULK:
								{
									out.write(reinterpret_cast<char *>(&it->Bulk(*jt)),recsize);
								} break;
								case DATA_REFERENCE:
								{
									Storage::reference_array arr = it->ReferenceArray(*jt);
									for(k = 0; k < recsize; k++)
									{
										if( arr[k].isValid() )
										{
											wetype = arr[k].GetElementType();
											out.put(wetype);
											lid = IntegerDF(arr[k]->GetHandle(),set_id);
											uconv.write_iValue(out,lid);
										}
										else out.put(NONE);
									}
								} break;
								case DATA_REMOTE_REFERENCE:
								{
									Storage::remote_reference_array arr = it->RemoteReferenceArray(*jt);
									for(k = 0; k < recsize; k++)
									{
										if( arr[k].isValid() )
										{
											uconv.write_iValue(out,static_cast<INMOST_DATA_ENUM_TYPE>(arr[k].GetMeshLink()->GetMeshName().size()));
											out.write(arr[k].GetMeshLink()->GetMeshName().c_str(),arr[k].GetMeshLink()->GetMeshName().size());
											wetype = arr[k].GetElementType();
											out.put(wetype);
											lid = IntegerDF(arr[k]->GetHandle(),set_id);
											uconv.write_iValue(out,lid);
										}
										else out.put(NONE);
									}
								} break;
#if defined(USE_AUTODIFF)
								case DATA_VARIABLE:
								{
									Storage::var_array arr = it->VariableArray(*jt);
									if( noderivs.find(jt->GetTagName()) != noderivs.end() )
									{
										//~ std::cout << "Saving " << jt->GetTagName() << " without derivatives" << std::endl;
										for(k = 0; k < recsize; k++)
										{
											uconv.write_fValue(out,arr[k].GetValue());
											uconv.write_iValue(out,0);
										}
									}
									else
									{
										//~ std::cout << "Saving " << jt->GetTagName() << " with derivatives" << std::endl;
										for(k = 0; k < recsize; k++)
										{
											const Sparse::Row & r = arr[k].GetRow();
											uconv.write_fValue(out,arr[k].GetValue());
											uconv.write_iValue(out,arr[k].GetRow().Size());
											for(int m = 0; m < (int)r.Size(); ++m)
											{
												uconv.write_fValue(out,r.GetValue(m));
												uconv.write_iValue(out,r.GetIndex(m));
											}
										}
									}
								} break;
#endif
							}
						}
						q++;
					}
					if( sparse ) uconv.write_iValue(out,q);
				}
		}
		DeleteTag(set_id);
		
		out << INMOST::EoMHeader;
#if defined(USE_MPI)
		if( m_state == Mesh::Parallel )
		{
			REPORT_STR("Parallel write");
			int ierr;
			INMOST_DATA_ENUM_TYPE numprocs = GetProcessorsNumber(),k;
			INMOST_DATA_BIG_ENUM_TYPE datasize, offset;
			std::vector<INMOST_DATA_BIG_ENUM_TYPE> datasizes(numprocs,0), offsets(numprocs,0);
			std::string local_data(out.str());
			datasize = local_data.size();
			REPORT_VAL("local_write_file_size",datasize);
			REPORT_MPI(ierr = MPI_Gather(&datasize,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,&datasizes[0],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,0,GetCommunicator()));
			if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
#if defined(USE_MPI_FILE) //We have access to MPI_File
			if( parallel_file_strategy == 1 )
			{
				MPI_File fh;
				MPI_Status stat;
				REPORT_MPI(ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(File.c_str()), MPI_MODE_CREATE | MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh));
				if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				REPORT_MPI(ierr = MPI_File_close(&fh));
				if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				REPORT_MPI(ierr = MPI_Barrier(GetCommunicator()));
				if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				REPORT_MPI(ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(File.c_str()),MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh));
				if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				if( GetProcessorRank() == 0 )
				{
					std::stringstream header;
					header.put(INMOST::INMOSTFile);
					buconv.write_iByteOrder(header);
					buconv.write_iByteSize(header);
					buconv.write_iValue(header,numprocs);
					for(k = 0; k < numprocs; k++) buconv.write_iValue(header,datasizes[k]);
					
					std::string header_data(header.str());
					REPORT_MPI(ierr = MPI_File_write(fh,&header_data[0],static_cast<INMOST_MPI_SIZE>(header_data.size()),MPI_CHAR,&stat));
					if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				
					MPI_Offset off;
					REPORT_MPI(ierr = MPI_File_get_position(fh,&off));
					if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
					offsets[0] = off;
					for(k = 1; k < numprocs; ++k)
						offsets[k] = offsets[k-1] + datasizes[k-1];
				}
				else offsets.resize(1,0);
				
				REPORT_MPI(ierr = MPI_Scatter(&offsets[0],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,&offset,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,0,GetCommunicator()));
				if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				REPORT_VAL("local write offset",offset);
				
				if( datasize )
				{
					INMOST_DATA_BIG_ENUM_TYPE shift = 0, chunk;
					//REPORT_MPI(ierr = MPI_File_write_ordered(fh,&local_data[0],static_cast<INMOST_MPI_SIZE>(local_data.size()),MPI_CHAR,&stat));
					while( shift != datasize )
					{
						chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),datasize-shift);
						REPORT_MPI(ierr = MPI_File_write_at(fh,offset+shift,local_data.c_str()+shift,static_cast<INMOST_MPI_SIZE>(chunk),MPI_CHAR,&stat));
						if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
						shift += chunk;
					}
				}
				
				REPORT_MPI(ierr = MPI_File_close(&fh));
				if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
			}
			else
#endif
			{
				std::string file_contents;
				std::string local_data(out.str());
				std::vector<MPI_Request> requests;
				if( GetProcessorRank() == 0 )
				{
					INMOST_DATA_BIG_ENUM_TYPE sizesum = datasizes[0];
					for(k = 1; k < numprocs; k++)
						sizesum += datasizes[k];
					file_contents.resize(sizesum);
					memcpy(&file_contents[0],&local_data[0],sizeof(char)*datasizes[0]);
					INMOST_DATA_BIG_ENUM_TYPE offset = datasizes[0];
					for(k = 1; k < numprocs; k++)
					{
						INMOST_DATA_BIG_ENUM_TYPE chunk, shift = 0;
						int it = 0; // for mpi tag
						while( shift != datasizes[k] )
						{
							MPI_Request req;
							chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),datasizes[k] - shift);
							REPORT_MPI(ierr = MPI_Irecv(&file_contents[offset+shift],chunk,MPI_CHAR,k, k*1000 + it, GetCommunicator(), &req) );
							if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
							requests.push_back(req);
							shift += chunk;
							it++;
							//TODO: remove temporary check
							if( it >= 1000 )
								std::cout << __FILE__ << ":" << __LINE__ << " too many iterations!!! " << it << " datasize " << datasizes[k] << std::endl;
						}
						offset += shift;
					}
				}
				else 
				{
					INMOST_DATA_BIG_ENUM_TYPE chunk, shift = 0;
					int it = 0; // for mpi tag
					while( shift != datasize )
					{
						MPI_Request req;
						chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),datasize - shift);
						REPORT_MPI(ierr = MPI_Isend(&local_data[shift],chunk,MPI_CHAR, 0, GetProcessorRank()*1000 + it, GetCommunicator(), &req) );
						if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
						requests.push_back(req);
						shift += chunk;
						it++;
						//TODO: remove temporary check
						if( it >= 1000 )
							std::cout << __FILE__ << ":" << __LINE__ << " too many iterations!!! " << it << " datasize " << datasize << std::endl;
					}
				}
				if( !requests.empty() )
				{
					REPORT_MPI(ierr = MPI_Waitall((int)requests.size(),&requests[0],MPI_STATUSES_IGNORE));
					if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				}
				//~ REPORT_MPI(ierr = MPI_Gatherv(&local_data[0],static_cast<INMOST_MPI_SIZE>(local_data.size()),MPI_CHAR,&file_contents[0],&recvcounts[0],&displs[0],MPI_CHAR,0,GetCommunicator()));
				//~ if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				if( GetProcessorRank() == 0 )
				{
					std::fstream fout(File.c_str(),std::ios::out | std::ios::binary);
					fout.put(INMOST::INMOSTFile);
					buconv.write_iByteOrder(fout);
					buconv.write_iByteSize(fout);
					buconv.write_iValue(fout,numprocs);
					for(k = 0; k < numprocs; k++) 
						buconv.write_iValue(fout,datasizes[k]);
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
			INMOST_DATA_ENUM_TYPE numprocs = 1;
			INMOST_DATA_BIG_ENUM_TYPE datasize = static_cast<INMOST_DATA_BIG_ENUM_TYPE>(out.tellp());
			REPORT_VAL("write_file_size",datasize);
			fout.put(INMOST::INMOSTFile);
			buconv.write_iByteOrder(fout);
			buconv.write_iByteSize(fout);
			buconv.write_iValue(fout,numprocs);
			buconv.write_iValue(fout,datasize);
			fout << out.rdbuf();
			fout.close();
		}
		
	}
	
	
	void Mesh::LoadPMF(std::string File)
	{
		int verbosity = 0;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < file_options.size(); ++k)
		{
			if( file_options[k].first == "VERBOSITY" )
			{
				verbosity = atoi(file_options[k].second.c_str());
				if( verbosity < 0 || verbosity > 2 )
				{
					std::cout << __FILE__ << ":" << __LINE__ << " Unknown verbosity option: " << file_options[k].second << "\n";
					verbosity = 1;
				}
			}
		}
		io_converter<INMOST_DATA_INTEGER_TYPE ,INMOST_DATA_REAL_TYPE> iconv;
		io_converter<INMOST_DATA_ENUM_TYPE    ,INMOST_DATA_REAL_TYPE> uconv;
		io_converter<INMOST_DATA_BIG_ENUM_TYPE,INMOST_DATA_REAL_TYPE> buconv;
		REPORT_STR("start load pmf");
		dynarray<INMOST_DATA_ENUM_TYPE,128> myprocs;
		std::stringstream in(std::ios::in | std::ios::out | std::ios::binary);
		HeaderType token;
		std::set< std::string > noload, loadonly, noderivs;
		
		noload = TagOptions("noload");
		loadonly = TagOptions("loadonly");
		noderivs = TagOptions("noderivatives");
		
		
#if defined(USE_MPI)
		if( m_state == Mesh::Parallel )
		{
#if defined(USE_MPI_FILE)
			if( parallel_file_strategy == 1 )
			{
				REPORT_STR("strategy 1");
				int ierr;
				std::vector<char> buffer;
				MPI_File fh;
				MPI_Status stat;
				REPORT_MPI(ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(File.c_str()),MPI_MODE_RDONLY,MPI_INFO_NULL,&fh));
				if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				INMOST_DATA_ENUM_TYPE numprocs = GetProcessorsNumber(), mpirank = GetProcessorRank();
				INMOST_DATA_BIG_ENUM_TYPE recvsize = 0, offset = 0;
				REPORT_VAL("number of processors",numprocs);
				REPORT_VAL("rank of processor", mpirank);
				std::vector<INMOST_DATA_BIG_ENUM_TYPE> recvsizes, offsets;
				if( mpirank == 0 ) //the alternative is to read alltogether
				{
					INMOST_DATA_ENUM_TYPE datanum, chunk,pos,k,q;\
					INMOST_DATA_BIG_ENUM_TYPE rdatanum;
					std::stringstream header;
					
					buffer.resize(3);
					//ierr = MPI_File_read_all(fh,&buffer[0],3,MPI_CHAR,&stat);
					REPORT_MPI(ierr = MPI_File_read(fh,&buffer[0],3,MPI_CHAR,&stat));
					if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
					
					if( static_cast<HeaderType>(buffer[0]) != INMOST::INMOSTFile ) throw BadFile;
					
					header.write(&buffer[1],2);
					buconv.read_iByteOrder(header);
					buconv.read_iByteSize(header);
					
					REPORT_VAL("integer_byte_order",buconv.str_iByteOrder(uconv.get_iByteOrder()));
					REPORT_VAL("integer_byte_size",(int)buconv.get_iByteSize());
					
					
					buffer.resize(buconv.get_source_iByteSize());
					//ierr = MPI_File_read_all(fh,&buffer[0],buffer.size(),MPI_CHAR,&stat);
					REPORT_MPI(ierr = MPI_File_read(fh,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),MPI_CHAR,&stat));
					if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
					
					header.write(&buffer[0],buffer.size());
					buconv.read_iValue(header,rdatanum);
					datanum = static_cast<INMOST_DATA_ENUM_TYPE>(rdatanum);
					
					REPORT_VAL("number of data entries",datanum);
					
					buffer.resize(static_cast<size_t>(datanum)*buconv.get_source_iByteSize());
					std::vector<INMOST_DATA_BIG_ENUM_TYPE> datasizes(datanum);
					//ierr = MPI_File_read_all(fh,&buffer[0],buffer.size(),MPI_CHAR,&stat);
					REPORT_MPI(ierr = MPI_File_read(fh,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),MPI_CHAR,&stat));
					if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
					
					INMOST_DATA_BIG_ENUM_TYPE datasum = 0;
					header.write(&buffer[0],buffer.size());
					for(k = 0; k < datanum; k++)
					{
						buconv.read_iValue(header,datasizes[k]);
						REPORT_VAL("size of data entry " << k,datasizes[k]);
						datasum += datasizes[k];
					}
					REPORT_VAL("total size",datasum);
					
					recvsizes.resize(numprocs,0);
					offsets.resize(numprocs,0);
					if( datanum <= recvsizes.size() )
					{
						REPORT_STR("number of processors is greater or equal then number of data entries");
						for(k = 0; k < datanum; k++)
							recvsizes[k] = datasizes[k];
					}
					else
					{
						REPORT_STR("number of processors is less then number of data entries - accumulating data");
						chunk = static_cast<INMOST_DATA_ENUM_TYPE>(floor(static_cast<double>(datanum)/static_cast<double>(recvsizes.size())));
						REPORT_VAL("chunk size" ,chunk);
						pos = 0;
						for(k = 0; k < recvsizes.size()-1; k++)
						{
							for(q = 0; q < chunk; q++)
								recvsizes[k] += datasizes[pos++];
							
							REPORT_VAL("recv on " << k, recvsizes[k]);
						}
						for(k = pos; k < datanum; k++)
							recvsizes[recvsizes.size()-1] += datasizes[k];
						REPORT_VAL("recv on " << recvsizes.size()-1, recvsizes[recvsizes.size()-1]);
					}
					MPI_Offset off;
					REPORT_MPI(ierr = MPI_File_get_position(fh,&off));
					if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
					
					//~ std::cout << "file offset: " << off << " size of big enum " << sizeof(INMOST_DATA_BIG_ENUM_TYPE) << std::endl;
					
					offsets[0] = off;
					for(k = 1; k < recvsizes.size(); k++)
						offsets[k] = offsets[k-1] + recvsizes[k-1];
				}
				else 
				{
					recvsizes.resize(1,0); //protect from dereferencing null
					offsets.resize(1,0);
				}
				
				REPORT_MPI(ierr = MPI_Scatter(&recvsizes[0],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,&recvsize,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,0,GetCommunicator()));
				if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				
				REPORT_MPI(ierr = MPI_Scatter(&offsets[0],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,&offset,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,0,GetCommunicator()));
				if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				
				REPORT_VAL("read on current processor",recvsize);
				REPORT_VAL("offset on current processor",offset);
				
				//~ std::cout << "read on " << GetProcessorRank() << " size " << recvsize << " offset " << offset << std::endl;
				//~ std::cout.flush();
				
				//~ {
					//~ REPORT_MPI(ierr = MPI_File_read_ordered(fh,&buffer[0],static_cast<INMOST_MPI_SIZE>(recvsize),MPI_CHAR,&stat));
					//~ if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
					//~ in.write(&buffer[0],recvsize);
				//~ }
				
				if( recvsize )
				{
					buffer.resize(recvsize);
					INMOST_DATA_BIG_ENUM_TYPE shift = 0, chunk;
					while( shift != recvsize )
					{
						chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),recvsize-shift);
						REPORT_MPI(ierr = MPI_File_read_at(fh,offset+shift,&buffer[shift],static_cast<INMOST_MPI_SIZE>(chunk),MPI_CHAR,&stat));
						if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
						shift += chunk;
					}
					in.write(&buffer[0],recvsize);
				}
				
				REPORT_MPI(ierr = MPI_File_close(&fh));
				if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
			}
			else
#endif
			{
				
				REPORT_STR("strategy 0");
				int ierr;
				std::vector<char> buffer, local_buffer;
				INMOST_DATA_BIG_ENUM_TYPE recvsize;
				
				INMOST_DATA_ENUM_TYPE numprocs = GetProcessorsNumber(),mpirank = GetProcessorRank();
				std::vector<INMOST_DATA_BIG_ENUM_TYPE> recvsizes(numprocs,0);
				std::vector<MPI_Request> requests;
				//~ std::vector<INMOST_MPI_SIZE> sendcnts(numprocs), displs(numprocs);
				REPORT_VAL("number of processors",numprocs);
				REPORT_VAL("rank of processor", mpirank);
				if( mpirank == 0 ) //zero reads everything
				{
					std::fstream fin(File.c_str(),std::ios::in | std::ios::binary);
					fin.get(token);
					if( token != INMOST::INMOSTFile ) throw BadFile;
					buconv.read_iByteOrder(fin);
					buconv.read_iByteSize(fin);
					
					REPORT_VAL("file position",fin.tellg());
					
					REPORT_VAL("integer_byte_order",buconv.str_iByteOrder(uconv.get_iByteOrder()));
					REPORT_VAL("integer_byte_size",(int)buconv.get_iByteSize());
					
					
					INMOST_DATA_ENUM_TYPE datanum,k,q,chunk,pos;
					INMOST_DATA_BIG_ENUM_TYPE datasum = 0, rdatanum;
					buconv.read_iValue(fin,rdatanum);
					datanum = static_cast<INMOST_DATA_ENUM_TYPE>(rdatanum);
					REPORT_VAL("number of data entries",datanum);
					
					REPORT_VAL("file position",fin.tellg());
					
					std::vector<INMOST_DATA_BIG_ENUM_TYPE> datasizes(datanum);
					for(k = 0; k < datanum; k++)
					{
						buconv.read_iValue(fin,datasizes[k]);
						REPORT_VAL("size of data entry " << k,datasizes[k]);
						datasum += datasizes[k];
					}
					
					REPORT_VAL("file position",fin.tellg());
					
					{
						buffer.resize(datasum);
						fin.read(&buffer[0],buffer.size());
					}
					REPORT_VAL("file position",fin.tellg());
					REPORT_VAL("total size",datasum);
					fin.close();
					
					
					if( datanum <= recvsizes.size() )
					{
						REPORT_STR("number of processors is greater or equal then number of data entries");
						for(k = 0; k < datanum; k++)
							recvsizes[k] = datasizes[k];
					}
					else
					{
						REPORT_STR("number of processors is less then number of data entries - accumulating data");
						chunk = static_cast<INMOST_DATA_ENUM_TYPE>(floor(static_cast<double>(datanum)/static_cast<double>(recvsizes.size())));
						REPORT_VAL("chunk size" ,chunk);
						pos = 0;
						for(k = 0; k < recvsizes.size()-1; k++)
						{
							for(q = 0; q < chunk; q++)
								recvsizes[k] += datasizes[pos++];
							REPORT_VAL("recv on " << k, recvsizes[k]);
						}
						for(k = pos; k < datanum; k++)
							recvsizes[recvsizes.size()-1] += datasizes[k];
						REPORT_VAL("recv on " << recvsizes.size()-1, recvsizes[recvsizes.size()-1]);
					}
					
					//~ displs[0] = 0;
					//~ sendcnts[0] = static_cast<INMOST_MPI_SIZE>(recvsizes[0]);
					//~ REPORT_VAL("disp on "<<0,displs[0]);
					//~ REPORT_VAL("send on "<<0,sendcnts[0]);
					
					//copy local
					local_buffer.resize(recvsizes[0]);
					memcpy(&local_buffer[0],&buffer[0],sizeof(char)*recvsizes[0]);
					in.write(&local_buffer[0],local_buffer.size());
					REPORT_VAL("output position",in.tellg());
					
					INMOST_DATA_BIG_ENUM_TYPE offset = recvsizes[0];
					for(k = 1; k < numprocs; k++)
					{
						//~ sendcnts[k] = static_cast<INMOST_MPI_SIZE>(recvsizes[k]);
						//~ displs[k] = sendcnts[k-1]+displs[k-1];
						//~ REPORT_VAL("disp on "<<k,displs[k]);
						//~ REPORT_VAL("send on "<<k,sendcnts[k]);
						
						INMOST_DATA_BIG_ENUM_TYPE chunk, shift = 0;
						int it = 0; // for mpi tag
						while( shift != recvsizes[k] )
						{
							MPI_Request req;
							chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),recvsizes[k] - shift);
							REPORT_MPI(ierr = MPI_Isend(&buffer[offset+shift],chunk,MPI_CHAR, k, k*1000 + it, GetCommunicator(), &req) );
							if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
							requests.push_back(req);
							shift += chunk;
							it++;
							//TODO: remove temporary check
							if( it >= 1000 )
								std::cout << __FILE__ << ":" << __LINE__ << " too many iterations!!! " << it << " datasize " << recvsizes[k] << std::endl;
						}
						offset += shift;
					}
				}
				
				REPORT_MPI(ierr = MPI_Scatter(&recvsizes[0],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,&recvsize,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,0,GetCommunicator()));
				if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				
				//~ std::cout << "on processor " << mpirank << " recieve " << recvsize << std::endl;
				
				REPORT_VAL("read on current processor",recvsize);
				
				if( mpirank && recvsize )
				{
					local_buffer.resize(recvsize);
					INMOST_DATA_BIG_ENUM_TYPE chunk, shift = 0;
					int it = 0; // for mpi tag
					while( shift != recvsize )
					{
						MPI_Request req;
						chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),recvsize - shift);
						REPORT_MPI(ierr = MPI_Irecv(&local_buffer[shift],chunk,MPI_CHAR, 0, GetProcessorRank()*1000 + it, GetCommunicator(), &req) );
						if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
						requests.push_back(req);
						shift += chunk;
						it++;
						//TODO: remove temporary check
						if( it >= 1000 )
							std::cout << __FILE__ << ":" << __LINE__ << " too many iterations!!! " << it << " datasize " << recvsize << std::endl;
					}
					in.write(&local_buffer[0],local_buffer.size());
					REPORT_VAL("output position",in.tellg());
				}
				
				if( !requests.empty() )
				{
					//~ std::cout << mpirank << " wait " << requests.size() << std::endl;
					REPORT_MPI(ierr = MPI_Waitall((int)requests.size(),&requests[0],MPI_STATUSES_IGNORE));
					if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				}
				//~ std::cout << mpirank << " is out " << std::endl;
				
				//~ REPORT_MPI(ierr = MPI_Scatterv(&buffer[0],&sendcnts[0],&displs[0],MPI_CHAR,&local_buffer[0],recvsize,MPI_CHAR,0,GetCommunicator()));
				//~ if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
				//~ in.write(&local_buffer[0],local_buffer.size());
				//~ REPORT_VAL("output position",in.tellg());
			}
			
		}
		else
#endif
		{
			std::fstream fin(File.c_str(),std::ios::in | std::ios::binary);
			fin.get(token);
			if( token != INMOST::INMOSTFile ) throw BadFile;
			buconv.read_iByteOrder(fin);
			buconv.read_iByteSize(fin);
			
			REPORT_VAL("file position",fin.tellg());
			
			REPORT_VAL("integer_byte_order",buconv.str_iByteOrder(uconv.get_iByteOrder()));
			REPORT_VAL("integer_byte_size",(int)buconv.get_iByteSize());
			
			INMOST_DATA_ENUM_TYPE datanum,k;
			INMOST_DATA_BIG_ENUM_TYPE datasum = 0, rdatanum;
			buconv.read_iValue(fin,rdatanum);
			datanum = static_cast<INMOST_DATA_ENUM_TYPE>(rdatanum);
			
			REPORT_VAL("file position",fin.tellg());
			
			REPORT_VAL("number of data entries",datanum);
			
			std::vector<INMOST_DATA_BIG_ENUM_TYPE> datasizes(datanum);
			for(k = 0; k < datanum; k++)
			{
				buconv.read_iValue(fin,datasizes[k]);
				REPORT_VAL("size of data entry " << k,datasizes[k]);
				datasum += datasizes[k];
			}
			
			REPORT_VAL("file position",fin.tellg());
			
			REPORT_VAL("total size",datasum);
			{
				std::vector<char> buffer;
				buffer.resize(datasum);
				fin.read(&buffer[0],buffer.size());
				in.write(&buffer[0],buffer.size());
				REPORT_VAL("output position",in.tellg());
			}
			
			REPORT_VAL("file position",fin.tellg());
			
			fin.close();
		}
		
		//std::fstream in(File.c_str(),std::ios::in | std::ios::binary);
		
		std::vector<Tag> tags;
		std::vector<bool> tags_skip;
		std::vector<DataType> tags_dtype;
		std::vector<std::string> tags_name;
		std::vector<ElementType> tags_defined;
		std::vector<ElementType> tags_sparse;
		std::vector<INMOST_DATA_ENUM_TYPE> tags_sizes;
		std::vector<HandleType> old_nodes;
		std::vector<HandleType> new_nodes;
		std::vector<HandleType> new_edges;
		std::vector<HandleType> new_faces;
		std::vector<HandleType> new_cells;
		std::vector<HandleType> new_sets;
		INMOST_DATA_ENUM_TYPE size,i,q;
		TopologyCheck tmp;
		INMOST_DATA_ENUM_TYPE current_dim = GetDimensions();
		
		bool start = false;
		
		std::map<GeometricData,ElementType> table;
		
		BeginModification();
		
		while (in >> token)
		{
			REPORT_VAL("output position, loop",in.tellg());
			if( !start )
			{
				if( token != INMOST::INMOSTFile ) throw BadFile; //check that this is valid file
				else
				{
					REPORT_STR("File chunk start read");
					//~ std::cout << "start read" << std::endl;
					tags.clear();
					tags_skip.clear();
					tags_dtype.clear();
					tags_name.clear();
					tags_defined.clear();
					tags_sparse.clear();
					tags_sizes.clear();
					old_nodes.clear();
					old_nodes.resize(NumberOfNodes());
					{
						unsigned qq = 0;
						for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it)
							old_nodes[qq++] = *it;
					}
					if( !old_nodes.empty() )
						std::sort(old_nodes.begin(),old_nodes.end(),CentroidComparator(this));
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
				
				
				current_dim = header[0];
				SetDimensions(header[0]);
				new_nodes.clear();
				new_nodes.resize(header[1]);
				new_edges.clear();
				new_edges.resize(header[2]);
				new_faces.clear();
				new_faces.resize(header[3]);
				new_cells.clear();
				new_cells.resize(header[4]);
				new_sets.clear();
				new_sets.resize(header[5]);
				tags.resize(header[6]);
				
				tags_skip.resize(header[6],false);
				tags_dtype.resize(header[6]);
				tags_name.resize(header[6]);
				tags_defined.resize(header[6]);
				tags_sparse.resize(header[6]);
				tags_sizes.resize(header[6]);
				//~ if( static_cast<Mesh::MeshState>(header[7]) == Mesh::Parallel && m_state != Mesh::Parallel)
				//~ SetCommunicator(INMOST_MPI_COMM_WORLD);
				myprocs.push_back(header[8]);
			}
			else if (token == INMOST::TagsHeader)
			{
				uconv.read_iValue(in,size);
				REPORT_STR("TagsHeader");
				REPORT_VAL("tag_size",size);
				for(i = 0; i < size; i++)
				{
					INMOST_DATA_ENUM_TYPE namesize;
					char name[4096];
					char datatype;
					char sparsemask,definedmask;
					uint64_t datalength;
					uconv.read_iValue(in,namesize);
					in.read(name, namesize);
					assert(namesize < 4096);
					name[namesize] = '\0';
					REPORT_VAL("tag name",name);
					in.get(datatype);
					REPORT_VAL("tag data type",DataTypeName(static_cast<DataType>(datatype)));
					in.get(sparsemask);
					in.get(definedmask);
					//for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype) )
					//{
					//  if( etype & definedmask ) REPORT_VAL("defined on",ElementTypeName(etype));
					//  if( etype & sparsemask ) REPORT_VAL("sparse on",ElementTypeName(etype));
					//}
					uconv.read_iValue_uint64_t(in,datalength);
					REPORT_VAL("length",datalength);
					
					
					{
						uint64_t undef;
						if( sizeof(uint8_t) == uconv.get_source_iByteSize() )
							undef = uint8_t(~uint8_t(0));
						else if( sizeof(uint16_t) == uconv.get_source_iByteSize() )
							undef = uint16_t(~uint16_t(0));
						else if( sizeof(uint32_t) == uconv.get_source_iByteSize() )
							undef = uint32_t(~uint32_t(0));
						if( sizeof(uint64_t) == uconv.get_source_iByteSize() )
							undef = uint64_t(~uint64_t(0));
						if( datalength == undef )
							datalength = ENUMUNDEF;
					}
					
					tags_skip[i] = false;
					tags_dtype[i] = static_cast<DataType>(datatype);
					tags_defined[i] = static_cast<ElementType>(definedmask);
					tags_sparse[i] = static_cast<ElementType>(sparsemask);
					tags_name[i] = std::string(name);
					tags_sizes[i] = (INMOST_DATA_ENUM_TYPE)datalength;
					
					tags_skip[i] = CheckLoadSkip(tags_name[i],noload,loadonly);
					
					if( !tags_skip[i] ) 
						tags[i] = CreateTag(std::string(name),static_cast<DataType>(datatype),
											static_cast<ElementType>(definedmask),
											static_cast<ElementType>(sparsemask),static_cast<INMOST_DATA_ENUM_TYPE>(datalength));
					
					REPORT_VAL("output position, tag " << i,in.tellg());
				}
			}
			else if (token == INMOST::NodeHeader)
			{
				uconv.read_iValue(in,size);
				assert(size == new_nodes.size());
				REPORT_STR("NodeHeader");
				REPORT_VAL("node_size",size);
				Storage::real coords[3] = {0,0,0};
				for(i = 0; i < size; i++)
				{
					for(unsigned int k = 0; k < current_dim; k++) iconv.read_fValue(in,coords[k]);
					int find = -1;
					if( !old_nodes.empty() )
					{
						std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),coords,CentroidComparator(this));
						if( it != old_nodes.end() )
						{
							Storage::real_array c = RealArrayDF(*it,CoordsTag());
							if( CentroidComparator(this).Compare(coords,c.data()) == 0 )
								find = static_cast<int>(it - old_nodes.begin());
						}
					}
					if( find == -1 ) new_nodes[i] = CreateNode(coords)->GetHandle();
					else  new_nodes[i] = old_nodes[find];
				}
				
				REPORT_VAL("output position, nodes",in.tellg());
			}
			else if (token == INMOST::EdgeHeader)
			{
				uconv.read_iValue(in,size);
				assert(size == new_edges.size());
				REPORT_STR("EdgeHeader");
				REPORT_VAL("edge_size",size);
				INMOST_DATA_ENUM_TYPE nlow, lid, i;
				ElementArray<Node> sub_elements(this);
				for(i = 0; i < size; i++)
				{
					uconv.read_iValue(in,nlow);
					for(q = 0; q < nlow; q++)
					{
						uconv.read_iValue(in,lid);
						sub_elements.push_back(new_nodes[lid]);
					}
					new_edges[i] = CreateEdge(sub_elements).first->GetHandle();
					sub_elements.clear();
				}
				
				REPORT_VAL("output position, edges",in.tellg());
			}
			else if (token == INMOST::FaceHeader)
			{
				uconv.read_iValue(in,size);
				assert(size == new_faces.size());
				REPORT_STR("FaceHeader");
				REPORT_VAL("face_size",size);
				INMOST_DATA_ENUM_TYPE nlow,lid,i;
				ElementArray<Edge> sub_elements(this);
				for(i = 0; i < size; i++)
				{
					uconv.read_iValue(in,nlow);
					for(q = 0; q < nlow; q++)
					{
						uconv.read_iValue(in,lid);
						sub_elements.push_back(new_edges[lid]);
					}
					new_faces[i] = CreateFace(sub_elements).first->GetHandle();
					sub_elements.clear();
				}
				
				REPORT_VAL("output position, faces",in.tellg());
			}
			else if (token == INMOST::CellHeader)
			{
				uconv.read_iValue(in,size);
				assert(size == new_cells.size());
				REPORT_STR("CellHeader");
				REPORT_VAL("cell_size",size);
				INMOST_DATA_ENUM_TYPE nlow, nhigh,lid;
				ElementArray<Face> sub_elements(this);
				ElementArray<Node> suggest_nodes(this);
				for(unsigned i = 0; i < size; i++)
				{
					uconv.read_iValue(in,nlow);
					for(q = 0; q < nlow; q++)
					{
						uconv.read_iValue(in,lid);
						sub_elements.push_back(new_faces[lid]);
					}
					uconv.read_iValue(in,nhigh);
					for(q = 0; q < nhigh; q++)
					{
						uconv.read_iValue(in,lid);
						suggest_nodes.push_back(new_nodes[lid]);
					}
					new_cells[i] = CreateCell(sub_elements, suggest_nodes).first->GetHandle();
					sub_elements.clear();
					suggest_nodes.clear();
				}
				
				REPORT_VAL("output position, cells",in.tellg());
			}
			else if (token == INMOST::ESetHeader)
			{
				uconv.read_iValue(in,size);
				assert(size == new_sets.size());
				REPORT_STR("EsetHeader");
				REPORT_VAL("eset_size",size);
				INMOST_DATA_ENUM_TYPE set_size, name_size, lid,val;
				char set_name[4096];
				HandleType * elem_links[4] =
				{
					new_nodes.empty() ? NULL : &new_nodes[0],
					new_edges.empty() ? NULL : &new_edges[0],
					new_faces.empty() ? NULL : &new_faces[0],
					new_cells.empty() ? NULL : &new_cells[0]
				};
				bool low_conn_have_sets = false;
				bool high_conn_have_sets = false;
				for(unsigned i = 0; i < size; i++)
				{
					uconv.read_iValue(in,name_size);
					in.read(set_name,name_size);
					assert(name_size < 4096);
					set_name[name_size] = '\0';
					new_sets[i] = CreateSet(std::string(set_name)).first->GetHandle();
					Element::adj_type & lc = LowConn(new_sets[i]);
					uconv.read_iValue(in,set_size);
					lc.resize(set_size);
					for(q = 0; q < set_size; ++q)
					{
						char type;
						in.get(type);
						if( type != 0 )
						{
							uconv.read_iValue(in,lid);
							if( static_cast<ElementType>(type) != ESET )
								lc[q] = elem_links[ElementNum(static_cast<ElementType>(type))][lid];
							else
							{
								lc[q] = ComposeHandle(static_cast<ElementType>(type),lid);
								low_conn_have_sets = true;
							}
						}
						else lc[q] = InvalidHandle();
					}
					Element::adj_type & hc = HighConn(new_sets[i]);
					//write tree information
					uconv.read_iValue(in,set_size);
					hc.resize(set_size);
					for(q = 0; q < ElementSet::high_conn_reserved-1; ++q)
					{
						char type;
						in.get(type);
						if( type != 0 )
						{
							uconv.read_iValue(in,lid);
							hc[q] = ComposeHandle(static_cast<ElementType>(type),static_cast<integer>(lid));
							high_conn_have_sets = true;
						}
						else hc[q] = InvalidHandle();
					}
					for(q = ElementSet::high_conn_reserved-1; q < set_size; ++q)
					{
						uconv.read_iValue(in,val);
						hc[q] = val;
					}
				}
				//convert handles to sets into links to new_sets
				if( high_conn_have_sets )
				{
					for(unsigned i = 0; i < size; i++)
					{
						Element::adj_type & hc = HighConn(new_sets[i]);
						for(enumerator j = 0; j < ElementSet::high_conn_reserved-1; ++j)
							if( hc[j] != InvalidHandle() )
								hc[j] = new_sets[GetHandleID(hc[j])];
					}
				}
				if( low_conn_have_sets ) //this may be expensive and redundant in some cases
				{
					for(unsigned i = 0; i < size; i++)
					{
						Element::adj_type & lc = LowConn(new_sets[i]);
						for(Element::adj_type::size_type j = 0; j < lc.size(); ++j)
							if( GetHandleElementType(lc[j]) == ESET )
								lc[j] = new_sets[GetHandleID(lc[j])];
					}
				}
				
				REPORT_VAL("output position, sets",in.tellg());
			}
			else if (token == INMOST::MeshDataHeader)
			{
				REPORT_STR("MeshDataHeader");
				HandleType m_storage = GetHandle();
				INMOST_DATA_ENUM_TYPE elem_sizes[6] =
				{
					static_cast<INMOST_DATA_ENUM_TYPE>(new_nodes.size()),
					static_cast<INMOST_DATA_ENUM_TYPE>(new_edges.size()),
					static_cast<INMOST_DATA_ENUM_TYPE>(new_faces.size()),
					static_cast<INMOST_DATA_ENUM_TYPE>(new_cells.size()),
					static_cast<INMOST_DATA_ENUM_TYPE>(new_sets.size()),
					1
				};
				HandleType * elem_links[6] =
				{
					new_nodes.empty() ? NULL : &new_nodes[0],
					new_edges.empty() ? NULL : &new_edges[0],
					new_faces.empty() ? NULL : &new_faces[0],
					new_cells.empty() ? NULL : &new_cells[0],
					new_sets.empty() ? NULL : &new_sets[0],
					&m_storage
				};
				for(INMOST_DATA_ENUM_TYPE j = 0; j < (INMOST_DATA_ENUM_TYPE)tags_skip.size(); j++)
				{
					REPORT_VAL("TagName",tags_name[j]);
					if( verbosity > 0 ) 
					{
						if( tags_skip[j] )
							std::cout << "Skipping " << tags_name[j] << std::endl;
						else
							std::cout << "Reading " << tags_name[j] << std::endl;
					}
					//Tag * jt = &tags[j];
					for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
					{
						if( etype & tags_defined[j] )
						{
							REPORT_VAL("defined on",ElementTypeName(etype));
							INMOST_DATA_ENUM_TYPE q, cycle_end, etypenum = ElementNum(etype);
							cycle_end = elem_sizes[etypenum];
							REPORT_VAL("cycle end",cycle_end);
							bool sparse = false;
							if( etype & tags_sparse[j] ) sparse = true;
							INMOST_DATA_ENUM_TYPE tagsize = tags_sizes[j], recsize = tagsize, lid;
							INMOST_DATA_ENUM_TYPE k;
							DataType data_type = tags_dtype[j];
							if( sparse )
							{
								REPORT_VAL("sparse on",ElementTypeName(etype));
								uconv.read_iValue(in,q);
								if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
							}
							else q = 0;
							
							REPORT_VAL("data type",DataTypeName(data_type));
							REPORT_VAL("tag size",tagsize);
							
							while(q != cycle_end)
							{
								HandleType he = elem_links[etypenum][q];
								if( tagsize == ENUMUNDEF )
								{
									uconv.read_iValue(in,recsize);
									if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
									if( !tags_skip[j] ) SetDataSize(he,tags[j],recsize);
								}
								switch(data_type)
								{
									case DATA_REAL:
									{
										if( tags_skip[j] )
										{
											Storage::real tmp;
											for(k = 0; k < recsize; k++)
											{
												iconv.read_fValue(in,tmp);
												if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
											}
										}
										else
										{
											Storage::real_array arr = RealArray(he,tags[j]);
											for(k = 0; k < recsize; k++)
											{
												iconv.read_fValue(in,arr[k]);
												if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
											}
										}
									} break;
									case DATA_INTEGER:
									{
										if( tags_skip[j] )
										{
											Storage::integer tmp;
											for(k = 0; k < recsize; k++)
											{
												iconv.read_iValue(in,tmp);
												if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
											}
										}
										else
										{
											Storage::integer_array arr = IntegerArray(he,tags[j]);
											for(k = 0; k < recsize; k++)
											{
												iconv.read_iValue(in,arr[k]);
												if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
											}
										}
									} break;
									case DATA_BULK:
									{
										if( tags_skip[j] )
										{
											std::vector<char> tmp(recsize);
											in.read(reinterpret_cast<char *>(&tmp[0]),recsize);
											if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
										}
										else
										{
											in.read(reinterpret_cast<char *>(&Bulk(he,tags[j])),recsize);
											if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
										}
									} break;
									case DATA_REFERENCE:
									{
										if( tags_skip[j] )
										{
											for(k = 0; k < recsize; k++)
											{
												char type;
												in.get(type);
												if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
												if (type != NONE)
												{
													uconv.read_iValue(in, lid);
													if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
												}
											}
										}
										else
										{
											Storage::reference_array arr = ReferenceArray(he,tags[j]);
											for(k = 0; k < recsize; k++)
											{
												char type;
												in.get(type);
												if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
												if (type != NONE)
												{
													uconv.read_iValue(in, lid);
													if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
													arr.at(k) = elem_links[ElementNum(type)][lid];
												}
												else arr.at(k) = InvalidHandle();
											}
										}
									} break;
									case DATA_REMOTE_REFERENCE:
									{
										Storage::remote_reference_array arr = RemoteReferenceArray(he,tags[j]);
										for(k = 0; k < recsize; k++)
										{
											if( tags_skip[j] )
											{
												INMOST_DATA_ENUM_TYPE size;
												std::vector<char> name;
												uconv.read_iValue(in,size);
												name.resize(size);
												if( !name.empty() ) in.read(&name[0],size);
												else std::cout << __FILE__ << ":" << __LINE__ << " Mesh of the name was not specified" << std::endl;
												char type;
												in.get(type);
												if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
												if (type != NONE)
												{
													uconv.read_iValue(in, lid);
													if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
												}
											}
											else
											{
												INMOST_DATA_ENUM_TYPE size;
												std::vector<char> name;
												uconv.read_iValue(in,size);
												name.resize(size);
												if( !name.empty() ) in.read(&name[0],size);
												else std::cout << __FILE__ << ":" << __LINE__ << " Mesh of the name was not specified" << std::endl;
												arr.at(k).first = GetMesh(std::string(name.begin(),name.end()));
												if( arr.at(k).first == NULL )
													std::cout << __FILE__ << ":" << __LINE__ << " Mesh with the name " << std::string(name.begin(),name.end()) << " do not exist, you should create the mesh with this name first" << std::endl;
												char type;
												in.get(type);
												if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
												if (type != NONE)
												{
													uconv.read_iValue(in, lid);
													if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
													arr.at(k).second = ComposeHandle(type,lid);
												}
												else arr.at(k).second = InvalidHandle();
											}
										}
									} break;
#if defined(USE_AUTODIFF)
									case DATA_VARIABLE:
									{
										if( tags_skip[j] )
										{
											Storage::real val;
											Storage::integer ival;
											for(k = 0; k < recsize; k++)
											{
												iconv.read_fValue(in,val);
												iconv.read_iValue(in,ival);
												int rend = ival;
												for(int l = 0; l < rend; ++l)
												{
													iconv.read_fValue(in,val);
													iconv.read_iValue(in,ival);
												}
											}
										}
										else
										{
											bool noder = false;
											if( noderivs.find(tags_name[j]) != noderivs.end() )
												noder = true;
											Storage::var_array arr = VariableArray(he,tags[j]);
											Storage::real val;
											Storage::integer ival;
											for(k = 0; k < recsize; k++)
											{
												Sparse::Row & r = arr[k].GetRow();
												iconv.read_fValue(in,val);
												arr[k].SetValue(val);
												iconv.read_iValue(in,ival);
												int rend = ival;
												if( noder )
												{
													for(int l = 0; l < rend; ++l)
													{
														iconv.read_fValue(in,val);
														iconv.read_iValue(in,ival);
													}
												}
												else
												{
													r.Resize(ival);
													for(int l = 0; l < rend; ++l)
													{
														iconv.read_fValue(in,val);
														iconv.read_iValue(in,ival);
														r.GetValue(l) = val;
														r.GetIndex(l) = ival;
													}
												}
											}
										}
									} break;
#endif
								}
								if( sparse )
								{
									uconv.read_iValue(in,q);
									if( in.eof() ) std::cout << __FILE__ << ":" << __LINE__ << " Unexpected end of file! " << tags_name[j] << " " << ElementTypeName(etype) << " " << (sparse? "sparse" : "dense") << std::endl;
								}
								else q++;
							}
						}
					}
					
					REPORT_VAL("output position, tag data " << j,in.tellg());
				}
				
				REPORT_VAL("output position, tag data",in.tellg());
				if( verbosity > 0 ) std::cout << "Finished reading data" << std::endl;
				REPORT_STR("EndOfData");
			}
			else
			{
				std::cout << "Unknown token on input" << std::endl;
				throw BadFile;
			}
		}
		
		
		
		if( m_state == Mesh::Parallel )
		{
#if defined(USE_MPI) 
			
			bool restore_state = false;
			INMOST_DATA_ENUM_TYPE numprocs = GetProcessorsNumber(), size = static_cast<INMOST_DATA_ENUM_TYPE>(myprocs.size()),k;
			INMOST_DATA_ENUM_TYPE procs_sum = 0;
			std::vector<INMOST_DATA_ENUM_TYPE> procs_sizes(numprocs);
			REPORT_MPI(MPI_Allgather(&size,1,INMOST_MPI_DATA_ENUM_TYPE,&procs_sizes[0],1,INMOST_MPI_DATA_ENUM_TYPE,GetCommunicator()));
			for(k = 0; k < numprocs; k++) if( procs_sizes[k] > 1 ) restore_state = true;
			REPORT_VAL("restore_state",restore_state);
			if( restore_state ) //we have to do something with parallel status data
			{
				int ierr;
				std::vector<INMOST_MPI_SIZE> recvcnts(numprocs), displs(numprocs);
				Storage::integer myrank = GetProcessorRank();
				std::sort(myprocs.begin(),myprocs.end());
				//have to allgatherv myprocs from all processors to detect the ownership of elements
				procs_sum = procs_sizes[0];
				recvcnts[0] = static_cast<INMOST_MPI_SIZE>(procs_sizes[0]);
				displs[0] = 0;
				for(k = 1; k < numprocs; k++)
				{
					procs_sum += procs_sizes[k];
					recvcnts[k] = static_cast<INMOST_MPI_SIZE>(procs_sizes[k]);
					displs[k] = displs[k-1]+recvcnts[k-1];
				}
				std::vector<INMOST_DATA_ENUM_TYPE> procs(procs_sum);
				REPORT_MPI(ierr = MPI_Allgatherv(myprocs.data(),procs_sizes[myrank],INMOST_MPI_DATA_ENUM_TYPE,&procs[0],&recvcnts[0],&displs[0],INMOST_MPI_DATA_ENUM_TYPE,GetCommunicator()));
				if( ierr != MPI_SUCCESS ) REPORT_MPI(MPI_Abort(GetCommunicator(),__LINE__));
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
							integer_array v = it->IntegerArrayDV(tag_processors);
							for(integer_array::iterator jt = v.begin(); jt != v.end(); ++jt)
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
							v.resize(static_cast<integer_array::size_type>(std::unique(v.begin(),v.end())-v.begin()));
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
			//ComputeSharedProcs();
			RecomputeParallelStorage(CELL | FACE | EDGE | NODE);
			
			//Share number of Layers
			REPORT_MPI(MPI_Bcast(&Integer(GetHandle(),tag_layers),1,INMOST_MPI_DATA_INTEGER_TYPE,0,GetCommunicator()));
			REPORT_MPI(MPI_Bcast(&Integer(GetHandle(),tag_bridge),1,INMOST_MPI_DATA_BULK_TYPE,0,GetCommunicator()));
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
			tag_global_id = GetTag("GLOBAL_ID");
		if( HaveTag("TOPOLOGY_ERROR_TAG") )
			tag_topologyerror = GetTag("TOPOLOGY_ERROR_TAG");
#if defined(USE_MPI)
		if( m_state == Mesh::Parallel )
		{
			ElementType have_global_id = NONE;
			if( GlobalIDTag().isValid() )
			{
				for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype))
					if( GlobalIDTag().isDefined(etype) ) have_global_id |= etype;
			}
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
}

#endif
