#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "inmost.h"

#if defined(USE_MESH)

namespace INMOST
{

  void Mesh::SavePVTK(std::string File)
  {
		std::string name=File;
		std::string::size_type pos=name.rfind(".pvtk");
		name.erase(pos); 

		std::string::size_type l=name.find_last_of("/\\");
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

  void Mesh::LoadPVTK(std::string File)
  {
		int state = 0, np = 0, nf = 0;
		std::string::size_type attrl, ql, qr, l;
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
					if( attrl == std::string::npos ) 
					{
						std::cout << __FILE__ << ":" << __LINE__ << " version information not found in " << File << std::endl;
						throw BadFile;
					}
					attrl = str.find("dataType=\"vtkUnstructuredGrid\"",l);
					if( attrl == std::string::npos ) 
					{
						std::cout << __FILE__ << ":" << __LINE__ << " data type information not found in " << File << std::endl;
						throw BadFile;
					}
					attrl = str.find("numberOfPieces");
					if( attrl == std::string::npos ) 
					{
						std::cout << __FILE__ << ":" << __LINE__ << " number of files not found in " << File << std::endl;
						throw BadFile;
					}
					ql = str.find_first_of('"',attrl);
					qr = str.find_first_of('"',ql+1);
					attrval = str.substr(ql+1,qr-ql-1);
					//std::cout << "number of peaces: " << attrval << " l: " << ql << " r: " << qr << std::endl;
					np = atoi(attrval.c_str());
					//std::cout << "np: " << np << std::endl;;
					files.resize(np);
					state = 1;
				}
				else 
				{
					std::cout << __FILE__ << ":" << __LINE__ << " unexpected xml tag " << tag << " in " << File << " expected File instead " << std::endl;
					throw BadFile;
				}
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
				else 
				{
					std::cout << __FILE__ << ":" << __LINE__ << " unexpected xml tag " << tag << " in " << File << " expected Piece instead " << std::endl;
					throw BadFile;
				}
			}
			else if( state == 2 )
			{
				if( tag != "/File" ) 
				{
					std::cout << __FILE__ << ":" << __LINE__ << " unexpected xml tag " << tag << " in " << File << " expected /File instead " << std::endl;
					throw BadFile;
				}
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
}

#endif
