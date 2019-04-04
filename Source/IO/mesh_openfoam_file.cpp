#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif


#include "inmost.h"
#include <cfloat>

#if defined(USE_MESH)



namespace INMOST
{
	
	
	class stream : public std::istream
	{
		std::string src;
		int hadlinebreak, linebreak;
		int hadlinechar, linechar;
		int state;
		
		stream(std::string file)
		:src(file), linebreak(0), hadlinebreak(0), hadlinechar(0), linechar(0), state(0)
		{
			open(file.c_str(),std::ios::in);
			if( fail() ) 
			{
				Report("Failed to open %s",src.c_str());
				state = -1;
			}
		}
		
		char GetChar()
		{
			char c = '\0';
			get(c);
			hadlinebreak = linebreak;
			hadlinechar = linechar;
			if( c == '\n' ) 
			{
				++linebreak;
				linechar = 0;
			}
			else ++linechar;
			if( eof() ) 
			{
				state = -100;
			}
			if( fail() )
			{
				Report("Stream failed while getting the char, state %d",state);
				state = -1;
			}
			return c;
		}
		void RetChar()
		{
			linebreak = hadlinebreak;
			linechar = hadlinechar;
			unget();
			if( fail() ) 
			{
				Report("Stream failed while ungetting the char");
				state = -1;
			}
		}
		void Report(const char * fmt, ...) const
		{ 
			std::cout << src << ":row:" << linebreak << ":col:" << linechar << " ";
			{
				char stext[16384];
				va_list ap;
				if ( fmt == NULL ) {std::cout << std::endl; return;}
				va_start(ap,fmt);
				vsprintf(stext,fmt,ap);
				va_end(ap);
				std::cout << stext;
			}
			std::cout << std::endl;
		}
	}
	
	
	
	
	
	void Mesh::LoadOpenFoam(std::string Folder)
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

		std::vector<HandleType> old_nodes(NumberOfNodes());
		{
			unsigned qq = 0;
			for (Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it)
				old_nodes[qq++] = *it;
		}
		
		
		
		if (!old_nodes.empty())
			std::sort(old_nodes.begin(), old_nodes.end(), CentroidComparator(this));

		std::vector<Tag> datatags;
		std::vector<HandleType> newnodes;
		std::vector<HandleType> newpolyh;
		std::vector<HandleType> newcells;
		
		char line[4096];
		
		stream(Folder+"/nodes");
		
		
		
		f.close();
	}
}

#endif
