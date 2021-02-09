#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "inmost.h"

#if defined(USE_MESH)
#include <deque>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include "io.hpp"
#include "../Misc/utils.h"










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

	

	void Mesh::SetFileOption(std::string key, std::string val)
	{
		for(INMOST_DATA_ENUM_TYPE k = 0; k < file_options.size(); k++)
		{
			if(file_options[k].first == key)
			{
				file_options[k].second = val;
			}
		}
		file_options.push_back(std::make_pair(key,val));
	}

	std::string Mesh::GetFileOption(std::string key) const
	{
		for(INMOST_DATA_ENUM_TYPE k = 0; k < file_options.size(); k++)
		{
			if(file_options[k].first == key)
			{
				return file_options[k].second;
			}
		}
		return "";
	}
	
	std::set<std::string> Mesh::TagOptions(std::string name) const
	{
		std::set<std::string> ret;
		std::vector<std::string> name_split;
		std::vector<std::string> value_split;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < file_options.size(); k++)
		{
			split_string(file_options[k].first,name_split,":");
			if(name_split[0] == "Tag")
			{
				//~ std::cout << "Tag option " << file_options[k].first << " value " << file_options[k].second << " tag " << name_split[1] << " searching " << name << std::endl;
				split_string(file_options[k].second,value_split,",");
				for(size_t q = 0; q < value_split.size(); ++q)
				{
					if( value_split[q] == name )
					{
						//~ std::cout << "found!" << std::endl;
						ret.insert(name_split[1]);
						break;
					}
				}
				value_split.clear();
			}
			name_split.clear();
		}
		return ret;
	}
	
	bool Mesh::CheckLoadSkip(std::string name, const std::set<std::string> & noload, const std::set<std::string> & loadonly) const
	{
		bool skip = false;
		if( !loadonly.empty() && loadonly.find(name) == loadonly.end() )
			skip = true;
		else if( noload.find(name) != noload.end() )
			skip = true;
		return skip;
	}
	
	bool Mesh::CheckSaveSkip(std::string name, const std::set<std::string> & nosave, const std::set<std::string> & saveonly) const
	{
		bool skip = false;
		if( !saveonly.empty() && saveonly.find(name) == saveonly.end() )
			skip = true;
		else if( nosave.find(name) != nosave.end() )
			skip = true;
		return skip;
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
			LoadECL(File);
		else if(LFile.find(".pvtk") != std::string::npos) //this is legacy parallel vtk
			LoadPVTK(File);
		else if(LFile.find(".pvtu") != std::string::npos) //this is new parallel vtk
			LoadPVTU(File);
		else if(LFile.find(".vtk") != std::string::npos) //this is legacy vtk
			LoadVTK(File);
		else if(LFile.find(".grid") != std::string::npos ) // mesh format by Mohammad Karimi-Fard
			LoadMKF(File);
		else if(LFile.find(".msh") != std::string::npos ) // GMSH file format
			LoadMSH(File);
		else if (LFile.find(".xml") != std::string::npos) //new mesh format 
			LoadXML(File);
		else if (LFile.find(".vtu") != std::string::npos) //this is new vtk
			LoadVTU(File);
		else if(LFile.find(".pmf") != std::string::npos) //this is inner parallel/platform mesh format
			LoadPMF(File);
		else 
		{
			std::cout << __FILE__ << ":" << __LINE__ << " File " << File << " cannot be loaded due to unknown extension." << std::endl;
			throw NotImplemented;
		}
		EXIT_FUNC();
	}
	
	void Mesh::Save(std::string File)
	{
		ENTER_FUNC();
		REPORT_VAL("File",File);
		std::string LFile;
		LFile.resize(File.size());
		std::transform(File.begin(),File.end(),LFile.begin(),::tolower);
		if(LFile.find(".pvtk") != std::string::npos) //this is legacy parallel vtk
			SavePVTK(File);
		else if(LFile.find(".pvtu") != std::string::npos) //this is new parallel vtk
			SavePVTU(File);
		else if(LFile.find(".xml") != std::string::npos)
			SaveXML(File);
		else if(LFile.find(".vtk") != std::string::npos) //this is legacy vtk
			SaveVTK(File);
		else if(LFile.find(".vtu") != std::string::npos) //this is new vtk
			SaveVTU(File);
		else if( LFile.find(".gmv") != std::string::npos) //this is gmv file
			SaveGMV(File);
		else if(LFile.find(".pmf") != std::string::npos) //this is inner parallel/platform mesh format
			SavePMF(File);
		else 
		{
			std::cout << __FILE__ << ":" << __LINE__ << " File " << File << " cannot be saved due to unknown extension." << std::endl;
			throw NotImplemented;	
		}
		EXIT_FUNC();
	}
	
	bool Mesh::isParallelFileFormat(std::string File)
	{
		std::string LFile;
		LFile.resize(File.size());
		std::transform(File.begin(),File.end(),LFile.begin(),::tolower);
		if(LFile.find(".grdecl") != std::string::npos) return true;
		if(LFile.find(".msh") != std::string::npos) return false;
		if(LFile.find(".grid") != std::string::npos) return false;
		else if(LFile.find(".vtk") != std::string::npos) return false;
		else if(LFile.find(".gmv") != std::string::npos) return false;
		else if(LFile.find(".pvtk") != std::string::npos) return true;
		else if(LFile.find(".pvtu") != std::string::npos) return true;
		else if(LFile.find(".pmf") != std::string::npos) return true;
		else if(LFile.find(".xml") != std::string::npos) return true;
		else if (LFile.find(".vtu") != std::string::npos) return false;
		else
		{
			std::cout << __FILE__ << ":" << __LINE__ << " File " << File << " format is unknown." << std::endl;
		}
		throw NotImplemented;
	}
}
#endif
