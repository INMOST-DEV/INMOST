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


//gmsh states

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

//eclipse states
#define ECLSTRCMP(x,y) strncmp(x,y,8)

#define ECL_NEVER -1
#define ECL_NONE 0
#define ECL_SKIP_SECTION 1
#define ECL_INCLUDE 2
#define ECL_DIMENS 3
#define ECL_DX 4
#define ECL_DY 5
#define ECL_DZ 6
#define ECL_TOPS 7
#define ECL_PERMX 8
#define ECL_PERMY 9
#define ECL_PERMZ 10
#define ECL_PORO 11
#define ECL_MAPAXIS 12
#define ECL_INRAD 13
#define ECL_COORDS 14
#define ECL_ZCORN 15

#define ECL_GTYPE_NONE 0
#define ECL_GTYPE_TOPS 1
#define ECL_GTYPE_ZCORN 2
#define ECL_GTYPE_RADIAL 3
#define ECL_GTYPE_CARTESIAN 4


#define ECL_VAR_NONE 0
#define ECL_VAR_REAL 1
#define ECL_VAR_INT 2

#define ECL_IJK_DATA(i,j,k) (i + (j+k*dims[1])*dims[0])



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

	std::string Mesh::GetFileOption(std::string key)
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

	
	void Mesh::Load(std::string File)
	{
		ENTER_FUNC();
		REPORT_VAL("File",File);
		std::string LFile;
		LFile.resize(File.size());
		std::transform(File.begin(),File.end(),LFile.begin(),::tolower);
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
			FILE * f = fopen(LFile.c_str(),"r");
			if( f == NULL )
			{
				std::cout << __FILE__ << ":" << __LINE__ << " cannot open file " << LFile << std::endl;
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
			std::vector< std::pair< std::pair<FILE *,std::string>, int> > fs(1,std::make_pair(std::make_pair(f,File),0));
			char readline[2048], *p, *pend, rec[2048];
			int text_end, text_start, state = ECL_NONE, nchars;
			int waitlines = 0;
			int have_dimens = 0, totread, downread, numrecs, offset;
			int gtype = ECL_GTYPE_NONE;
			int argtype = ECL_VAR_NONE;
			int radial = ECL_GTYPE_NONE;
			Storage::real * read_arrayf = NULL;
			Storage::integer * read_arrayi = NULL;
			Storage::integer dims[3], mapaxis[6] = {0,1,0,0,1,0};
			Storage::real inrad = 0;
			std::vector<Storage::real> xyz,perm,poro, tops,zcorn;
			std::vector<Storage::integer> actnum;
			while(!fs.empty())
			{
				while(fgets(readline,2048,fs.back().first.first) != NULL)
				{
					fs.back().second++; //line number
					{
						if( readline[strlen(readline)-1] == '\n' ) readline[strlen(readline)-1] = '\0';
						text_end = static_cast<int>(strlen(readline));
						for(text_start = 0; isspace(readline[text_start]) && text_start < text_end; text_start++);
						if( text_start == text_end ) continue;
						for(text_end = text_end-1; isspace(readline[text_end]) && text_end > text_start; text_end--);
						readline[text_end+1] = '\0';
						p = readline + text_start;
						pend = readline + text_end + 1;
						for(char * q = p; q < pend; q++) *q = toupper(*q);
					}
					if( p[0] == '-' && p[1] == '-' ) continue; //skip comment
					if(waitlines) {waitlines--; continue;} //skip meaningful lines
					switch(state)
					{
					case ECL_NONE:
						if( !ECLSTRCMP(p,"END") ) //end of data - don't read beyond
						{
							goto ecl_exit_loop;
						}
						else if( !ECLSTRCMP(p,"INCLUDE") ) state = ECL_INCLUDE;
						else if( !ECLSTRCMP(p,"DIMENS") || !ECLSTRCMP(p,"SPECGRID") )
						{
							read_arrayi = dims;
							numrecs = 1;
							downread = totread = 3;
							argtype = ECL_VAR_INT;
							offset = state = ECL_DIMENS;
							have_dimens = 1;
						}
						else if( !ECLSTRCMP(p,"MAPAXIS") )
						{
							read_arrayi = mapaxis;
							numrecs = 1;
							downread = totread = 6;
							argtype = ECL_VAR_INT;
							offset = state = ECL_MAPAXIS;
						}
						else if( !ECLSTRCMP(p,"DX") )
						{
							assert(have_dimens);
							if( xyz.empty() ) xyz.resize(3*dims[0]*dims[1]*dims[2]);
							read_arrayf = xyz.data();
							numrecs = 3;
							downread = totread = dims[0]*dims[1]*dims[2];
							argtype = ECL_VAR_REAL;
							offset = state = ECL_DX;
						}
						else if( !ECLSTRCMP(p,"DY") )
						{
							assert(have_dimens);
							if( xyz.empty() ) xyz.resize(3*dims[0]*dims[1]*dims[2]);
							read_arrayf = xyz.data();
							numrecs = 3;
							downread = totread = dims[0]*dims[1]*dims[2];
							argtype = ECL_VAR_REAL;
							offset = ECL_DX;
							state = ECL_DY;
						}
						else if( !ECLSTRCMP(p,"DZ") )
						{
							assert(have_dimens);
							if( xyz.empty() ) xyz.resize(3*dims[0]*dims[1]*dims[2]);
							read_arrayf = xyz.data();
							numrecs = 3;
							downread = totread = dims[0]*dims[1]*dims[2];
							argtype = ECL_VAR_REAL;
							offset = ECL_DX;
							state = ECL_DZ;
						}
						else if( !ECLSTRCMP(p,"COORD") )
						{
							assert(have_dimens);
							if( xyz.empty() ) xyz.resize(2*(dims[0]+1)*(dims[1]+1));
							read_arrayf = xyz.data();
							numrecs = 1;
							downread = totread = 2*(dims[0]+1)*(dims[1]+1);
							argtype = ECL_VAR_REAL;
							offset = state = ECL_COORDS;
							gtype = ECL_GTYPE_ZCORN;
						}
						else if( !ECLSTRCMP(p,"ZCORN") )
						{
							assert(have_dimens);
							if( zcorn.empty() ) zcorn.resize(dims[0]*dims[1]*dims[2]*8);
							read_arrayf = zcorn.data();
							numrecs = 1;
							downread = totread = dims[0]*dims[1]*dims[2]*8;
							argtype = ECL_VAR_REAL;
							state = offset = ECL_ZCORN;
						}
						else if( !ECLSTRCMP(p,"TOPS") )
						{
							assert(have_dimens);
							tops.resize(dims[0]*dims[1]*dims[2]);
							read_arrayf = tops.data();
							numrecs = 1;
							downread = totread = dims[0]*dims[1]*dims[2];
							argtype = ECL_VAR_REAL;
							offset = state = ECL_TOPS;
							gtype = ECL_GTYPE_TOPS;
						}
						else if( !ECLSTRCMP(p,"PERMX") )
						{
							assert(have_dimens);
							if( perm.empty() ) perm.resize(3*dims[0]*dims[1]*dims[2]);
							read_arrayf = perm.data();
							numrecs = 3;
							downread = totread = dims[0]*dims[1]*dims[2];
							argtype = ECL_VAR_REAL;
							offset = state = ECL_PERMX;
						}
						else if( !ECLSTRCMP(p,"PERMY") )
						{
							assert(have_dimens);
							if( perm.empty() ) perm.resize(3*dims[0]*dims[1]*dims[2]);
							read_arrayf = perm.data();
							numrecs = 3;
							downread = totread = dims[0]*dims[1]*dims[2];
							argtype = ECL_VAR_REAL;
							offset = ECL_PERMX;
							state = ECL_PERMY;
						}
						else if( !ECLSTRCMP(p,"PERMZ") )
						{
							assert(have_dimens);
							if( perm.empty() ) perm.resize(3*dims[0]*dims[1]*dims[2]);
							read_arrayf = perm.data();
							numrecs = 3;
							downread = totread = dims[0]*dims[1]*dims[2];
							argtype = ECL_VAR_REAL;
							offset = ECL_PERMX;
							state = ECL_PERMZ;
						}
						else if( !ECLSTRCMP(p,"PORO") )
						{
							assert(have_dimens);
							poro.resize(3*dims[0]*dims[1]*dims[2]);
							read_arrayf = poro.data();
							numrecs = 1;
							downread = totread = dims[0]*dims[1]*dims[2];
							argtype = ECL_VAR_REAL;
							offset = state = ECL_PORO;
						}
						else if( !ECLSTRCMP(p,"RADIAL") )
						{
							radial = ECL_GTYPE_RADIAL;
						}
						else if( !ECLSTRCMP(p,"CART") )
						{
							radial = ECL_GTYPE_CARTESIAN;
						}
						else if( !ECLSTRCMP(p,"INRAD") )
						{
							if( radial != ECL_GTYPE_RADIAL ) 
							{
								std::cout << __FILE__ << ":" << __LINE__ << " inner radius specified for cartesian grid ";
								std::cout << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
							}
							if( radial == ECL_GTYPE_NONE ) radial = ECL_GTYPE_RADIAL;
							state = ECL_INRAD;
						}
						else
						{
							std::cout << __FILE__ << ":" << __LINE__ << " skipped " << p << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
						}
						break;
					case ECL_SKIP_SECTION:
						if( *(pend-1) == '/' || *p == '/' ) state = ECL_NONE;
						break;
					case ECL_INCLUDE:
						if( 1 == sscanf(p," %s",rec) )
						{
							int shift_one = 0;
							if( (rec[0] == '\'' || rec[0] == '"') && rec[0] == rec[strlen(rec)-1] ) //remove quotes
							{
								rec[strlen(rec)-1] = '\0';
								shift_one = 1;
							}
							f = fopen(rec+shift_one,"r");
							if( f == NULL )
							{
								std::cout << __FILE__ << ":" << __LINE__ << " cannot open file " << rec+shift_one << " included from "<< fs.back().first.second << " line " << fs.back().second << std::endl;
								throw BadFileName;
							}
							fs.push_back(std::make_pair(std::make_pair(f,std::string(rec+shift_one)),0));
							if( *(pend-1) == '/' ) state = ECL_NONE; else state = ECL_SKIP_SECTION;
						}
						else
						{
							std::cout << __FILE__ << ":" << __LINE__ << " cannot read file name, string " << p << " in " << fs.back().first.second << " line " << fs.back().second << std::endl;
							throw BadFile;
						}
						break;
					case ECL_ZCORN:
					case ECL_COORDS:
					case ECL_MAPAXIS:
					case ECL_DIMENS:
					case ECL_DX:
					case ECL_DY:
					case ECL_DZ:
					case ECL_PERMX:
					case ECL_PERMY:
					case ECL_PERMZ:
					case ECL_PORO:
					case ECL_TOPS:
						while( downread > 0 && p < pend )
						{
							if( 1 == sscanf(p,"%s%n",rec,&nchars) )
							{
								p += nchars;
								while(isspace(*p) && p < pend) ++p;
								int count = 1;
								int begval = 0;
								for(int q = 0; q < static_cast<int>(strlen(rec)); ++q)
									if( rec[q] == '*' )
									{
										begval = q+1;
										rec[q] = '\0';
										break;
									}
								if( begval > 0 ) count = atoi(rec);
								if( argtype == ECL_VAR_REAL )
								{
									Storage::real val = atof(rec+begval);
									while(count)
									{
										read_arrayf[numrecs*(totread-(downread--))+(state-offset)] = val;
										count--;
									}
								}
								else if( argtype == ECL_VAR_INT )
								{
									Storage::integer val = atoi(rec+begval);
									while(count)
									{
										read_arrayi[numrecs*(totread-(downread--))+(state-offset)] = val;
										count--;
									}
								}
								else
								{
									std::cout << __FILE__ << ":" << __LINE__ << " probably forgot to set up argument type to read " << std::endl;
									throw Impossible;
								}
							}
							else
							{
								std::cout << __FILE__ << ":" << __LINE__ << " cannot read data " << p << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
								throw BadFile;
							}
							if( *p == '/' ) 
							{
								if( downread > 0 )
								{
									std::cout << __FILE__ << ":" << __LINE__ << " early data termination, read " << totread-downread << " of " << totread << " records in " << fs.back().first.second << ":" << fs.back().second << std::endl;
									throw BadFile;
								}
							}
						}
						if( *(pend-1) == '/' || *p == '/' ) state = ECL_NONE; else state = ECL_SKIP_SECTION;
						break;
					case ECL_INRAD:
						if( 1 == sscanf(p,"%lf%n",&inrad,&nchars) )
						{
							p += nchars;
							while(isspace(*p) && p < pend) ++p;
						}
						else
						{
							std::cout << __FILE__ << ":" << __LINE__ << " cannot read data " << p << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
							throw BadFile;
						}
						if( *(pend-1) == '/' || *p == '/' ) state = ECL_NONE; else state = ECL_SKIP_SECTION;
						break;
					}
				}
ecl_exit_loop:
				fclose(fs.back().first.first);
				fs.pop_back();
			}
			if( radial == ECL_GTYPE_RADIAL )
			{
				std::cout << __FILE__ << ":" << __LINE__ << " radial grids not supported yet " << std::endl;
			}
			if( gtype == ECL_GTYPE_TOPS )
			{
				std::vector<HandleType> newnodes((dims[0]+1)*(dims[1]+1)*(dims[2]+1));
				Storage::real x, y, z, node_xyz[3];
				x = 0.0;
				int numnode = 0;
				for(int i = 0; i < dims[0]+1; i++)
				{
					Storage::integer pif = std::min(dims[0]-1,i), pib = std::max(i-1,0);
					y = 0.0;
					for(int j = 0; j < dims[1]+1; j++)
					{
						Storage::integer pjf = std::min(dims[1]-1,j), pjb = std::max(j-1,0);
						z = (
							 tops[ECL_IJK_DATA(pib,pjb,0)]+
							 tops[ECL_IJK_DATA(pib,pjf,0)]+
							 tops[ECL_IJK_DATA(pif,pjb,0)]+
							 tops[ECL_IJK_DATA(pif,pjf,0)]
						    )*0.25;
						z -= (
							  xyz[3*ECL_IJK_DATA(pib,pjb,0)+2]+
							  xyz[3*ECL_IJK_DATA(pib,pjf,0)+2]+
							  xyz[3*ECL_IJK_DATA(pif,pjb,0)+2]+
							  xyz[3*ECL_IJK_DATA(pif,pjf,0)+2]
							 )*0.25;
						for(int k = 0; k < dims[2]+1; k++)
						{
							Storage::integer pkf = std::min(dims[2]-1,k), pkb = std::max(k-1,0);
							node_xyz[0] = x;
							node_xyz[1] = y;
							node_xyz[2] = z;
							int find = -1;
							if( !old_nodes.empty() )
							{
								std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),node_xyz,CentroidComparator(this));
								if( it != old_nodes.end() ) 
								{
									Storage::real_array c = RealArrayDF(*it,CoordsTag());
									if( CentroidComparator(this).Compare(node_xyz,c.data()) == 0 )
										find = static_cast<int>(it - old_nodes.begin());
								}
							}
							if( find == -1 ) newnodes[numnode++] = CreateNode(node_xyz)->GetHandle();
							else newnodes[numnode++] = old_nodes[find];
							//std::cout << i << " " << j << " " << k << " ( " << x << " , " << y << " , " << z << ") " << newnodes.back()->LocalID() << std::endl; 
							x += (
								  (
								   xyz[3*ECL_IJK_DATA(pib,pjb,pkf)+0]+
								   xyz[3*ECL_IJK_DATA(pib,pjf,pkf)+0]+
								   xyz[3*ECL_IJK_DATA(pif,pjb,pkf)+0]+
								   xyz[3*ECL_IJK_DATA(pif,pjf,pkf)+0]
								  )
								-
								  (
								   xyz[3*ECL_IJK_DATA(pib,pjb,pkb)+0]+
								   xyz[3*ECL_IJK_DATA(pib,pjf,pkb)+0]+
								   xyz[3*ECL_IJK_DATA(pif,pjb,pkb)+0]+
								   xyz[3*ECL_IJK_DATA(pif,pjf,pkb)+0]
								  )
								 )*0.25;
							y += (
								  (
							       xyz[3*ECL_IJK_DATA(pib,pjb,pkf)+1]+
								   xyz[3*ECL_IJK_DATA(pib,pjf,pkf)+1]+
								   xyz[3*ECL_IJK_DATA(pif,pjb,pkf)+1]+
								   xyz[3*ECL_IJK_DATA(pif,pjf,pkf)+1]
								  )
								-
								  (
								   xyz[3*ECL_IJK_DATA(pib,pjb,pkb)+1]+
								   xyz[3*ECL_IJK_DATA(pib,pjf,pkb)+1]+
								   xyz[3*ECL_IJK_DATA(pif,pjb,pkb)+1]+
								   xyz[3*ECL_IJK_DATA(pif,pjf,pkb)+1]
								  )
								 )*0.25;
							z += (
								  xyz[3*ECL_IJK_DATA(pib,pjb,pkb)+2]+
								  xyz[3*ECL_IJK_DATA(pib,pjf,pkb)+2]+
								  xyz[3*ECL_IJK_DATA(pif,pjb,pkb)+2]+
								  xyz[3*ECL_IJK_DATA(pif,pjf,pkb)+2]+
							      xyz[3*ECL_IJK_DATA(pib,pjb,pkf)+2]+
								  xyz[3*ECL_IJK_DATA(pib,pjf,pkf)+2]+
								  xyz[3*ECL_IJK_DATA(pif,pjb,pkf)+2]+
								  xyz[3*ECL_IJK_DATA(pif,pjf,pkf)+2]
								 )*0.125;
						}
						y += (
							  xyz[3*ECL_IJK_DATA(pib,pjb,0)+1]+
							  xyz[3*ECL_IJK_DATA(pif,pjb,0)+1]+
							  xyz[3*ECL_IJK_DATA(pib,pjf,0)+1]+
							  xyz[3*ECL_IJK_DATA(pif,pjf,0)+1]
							 )*0.25; 
					}
					x += (
						  xyz[3*ECL_IJK_DATA(pib,0,0)+0]+
						  xyz[3*ECL_IJK_DATA(pif,0,0)+0]
						 )*0.5; 
				}
				Tag tagporo,tagperm;
				if( !poro.empty() ) tagporo = CreateTag("PORO",DATA_REAL,CELL,NONE,1);
				if( !perm.empty() ) tagperm = CreateTag("PERM",DATA_REAL,CELL,NONE,3);

				const Storage::integer nvf[24] = { 2, 3, 1, 0, 4, 5, 7, 6, 0, 1, 5, 4, 3, 2, 6, 7, 2, 0, 4, 6, 1, 3, 7, 5 };
				const Storage::integer numnodes[6] = { 4, 4, 4, 4, 4, 4 };
				for(int i = 0; i < dims[0]; i++)
					for(int j = 0; j < dims[1]; j++)
						for(int k = 0; k < dims[2]; k++)
						{
							HandleType verts[8];
							verts[0] = newnodes[((i+0)*(dims[1]+1)+(j+0))*(dims[2]+1)+(k+0)];
							verts[1] = newnodes[((i+1)*(dims[1]+1)+(j+0))*(dims[2]+1)+(k+0)];
							verts[2] = newnodes[((i+0)*(dims[1]+1)+(j+1))*(dims[2]+1)+(k+0)];
							verts[3] = newnodes[((i+1)*(dims[1]+1)+(j+1))*(dims[2]+1)+(k+0)];
							verts[4] = newnodes[((i+0)*(dims[1]+1)+(j+0))*(dims[2]+1)+(k+1)];
							verts[5] = newnodes[((i+1)*(dims[1]+1)+(j+0))*(dims[2]+1)+(k+1)];
							verts[6] = newnodes[((i+0)*(dims[1]+1)+(j+1))*(dims[2]+1)+(k+1)];
							verts[7] = newnodes[((i+1)*(dims[1]+1)+(j+1))*(dims[2]+1)+(k+1)];
							//for(int q = 0; q < 8; q++)
							//	std::cout << verts[q]->Coords()[0] << " " << verts[q]->Coords()[1] << " " << verts[q]->Coords()[2] << " " << verts[q]->LocalID() << std::endl;
							Cell c = CreateCell(ElementArray<Node>(this,verts,verts+8),nvf,numnodes,6).first;
							if( !poro.empty() ) c->RealDF(tagporo) = poro[(i*dims[1]+j)*dims[2]+k];
							if( !perm.empty() )
							{
								Storage::real_array arr_perm = c->RealArrayDF(tagperm);
								arr_perm[0] = perm[3*((i*dims[1]+j)*dims[2]+k)+0];
								arr_perm[1] = perm[3*((i*dims[1]+j)*dims[2]+k)+1];
								arr_perm[2] = perm[3*((i*dims[1]+j)*dims[2]+k)+2];
							}
						}
			}
			else if( gtype == ECL_GTYPE_ZCORN )
			{
				//std::vector<HandleType> block_nodes(dims[0]*dims[1]*dims[2]*8);
				//Storage::integer block_zcorn[8];
			}
		}
		else
		if(LFile.find(".pvtk") != std::string::npos) //this is legacy parallel vtk
		{
			int state = 0, np, nf = 0;
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
		else if(LFile.find(".vtk") != std::string::npos) //this is legacy vtk
		{
			MarkerType unused_marker = CreateMarker();
			bool grid_is_2d = false;
			for(INMOST_DATA_ENUM_TYPE k = 0; k < file_options.size(); ++k)
			{
				if( file_options[k].first == "VTK_GRID_DIMS" )
				{
					if( atoi(file_options[k].second.c_str()) == 2 )
					{
						grid_is_2d = true;
						break;
					}
				}
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
									newnodes[i] = CreateNode(coords)->GetHandle();
								else 
									newnodes[i] = old_nodes[find];
								SetMarker(newnodes[i],unused_marker);
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
									std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),coords,CentroidComparator(this));
									if( it != old_nodes.end() ) 
									{
										Storage::real_array c = RealArrayDF(*it,CoordsTag());
										if( CentroidComparator(this).Compare(c.data(),coords) == 0 )
											find = static_cast<int>(it - old_nodes.begin());
									}
								}
								if( find == -1 ) newnodes[i] = CreateNode(coords)->GetHandle();
								else newnodes[i] = old_nodes[find];
								SetMarker(newnodes[i],unused_marker);
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
										if( grid_is_2d )
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
										if( grid_is_2d )
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
										if( grid_is_2d )
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
										if( grid_is_2d )
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
		else if(LFile.find(".grid") != std::string::npos ) // mesh format by Mohammad Karimi-Fard
		{
			Tag volume_factor, porosity, permiability, zone;
			zone = CreateTag("MATERIAL",DATA_INTEGER,CELL|FACE|ESET,ESET,1);
			volume_factor = CreateTag("VOLUME_FACTOR",DATA_REAL,CELL|FACE,FACE,1);
			porosity = CreateTag("PORO",DATA_REAL,CELL|FACE,FACE,1);
			permiability = CreateTag("PERM",DATA_REAL,CELL|FACE,FACE,9);
			std::vector<HandleType> old_nodes(NumberOfNodes());
			{
				unsigned qq = 0;
				for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it)
					old_nodes[qq++] = *it;
			}
			if( !old_nodes.empty() ) 
				std::sort(old_nodes.begin(),old_nodes.end(),CentroidComparator(this));

			FILE * f = fopen(File.c_str(),"r");
			if( f == NULL ) 
			{
				std::cout << __FILE__ << ":" << __LINE__ << " cannot open " << File << std::endl;
				throw BadFileName;
			}
			int nbnodes, nbpolygon, nbpolyhedra, nbzones, volcorr, nbfacenodes, nbpolyhedronfaces, num, nbK;
			int report_pace;
			Storage::real K[9],readK[9], poro, vfac;
			std::vector<HandleType> newnodes;
			std::vector<HandleType> newpolygon;
			std::vector<HandleType> newpolyhedron;
			ElementArray<ElementSet> newsets(this);
			ElementArray<Node> f_nodes(this);
			ElementArray<Face> c_faces(this);
			fscanf(f," %d %d %d %d %d",&nbnodes,&nbpolygon,&nbpolyhedra,&nbzones,&volcorr);
			newnodes.resize(nbnodes);
			newpolygon.resize(nbpolygon);
			newpolyhedron.resize(nbpolyhedra);
			newsets.resize(nbzones);
			if( verbosity > 0 ) printf("Creating %d sets for zones.\n",nbzones);
			report_pace = std::max<int>(nbzones/250,1);
			for(int i = 0; i < nbzones; i++)
			{
				std::stringstream str;
				str << "ZONE_" << i << "_SET";
				newsets[i] = CreateSet(str.str()).first;
				newsets[i]->Integer(zone) = i;
				if( verbosity > 1 &&  i % report_pace == 0 )
				{
					printf("sets %3.1f%%\r",(i*100.0)/(1.0*nbzones));
					fflush(stdout);
				}
			}
			if( verbosity > 0 ) printf("Reading %d nodes.\n",nbnodes);
			report_pace = std::max<int>(nbnodes/250,1);
			for(int i = 0; i < nbnodes; i++)
			{
				Storage::real xyz[3];
				if( 3 == fscanf(f," %lf %lf %lf",xyz,xyz+1,xyz+2) )
				{
					int find = -1;
					if( !old_nodes.empty() )
					{
						std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),xyz,CentroidComparator(this));
						if( it != old_nodes.end() ) 
						{
							Storage::real_array c = RealArrayDF(*it,CoordsTag());
							if( CentroidComparator(this).Compare(xyz,c.data()) == 0 )
								find = static_cast<int>(it - old_nodes.begin());
						}
					}
					if( find == -1 ) newnodes[i] = CreateNode(xyz)->GetHandle();
					else newnodes[i] = old_nodes[find];
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " cannot read coordinates from " << File << std::endl;
					throw BadFile;
				}
				if( verbosity > 1 &&  i % report_pace == 0 )
				{
					printf("nodes %3.1f%%\r",(i*100.0)/(1.0*nbnodes));
					fflush(stdout);
				}
			}
			if( verbosity > 0 ) printf("Reading %d faces.\n",nbpolygon);
			report_pace = std::max<int>(nbpolygon/250,1);
			for(int i = 0; i < nbpolygon; i++)
			{
				if( 1 != fscanf(f," %d",&nbfacenodes) )
				{
					std::cout << __FILE__ << ":" << __LINE__ << " cannot read number of nodes from " << File << std::endl;
					throw BadFile;
				}
				for(int j = 0; j < nbfacenodes; j++)
				{
					if( 1 != fscanf(f," %d", &num) )
					{
						std::cout << __FILE__ << ":" << __LINE__ << " cannot read node number from " << File << std::endl;
						throw BadFile;
					}
					if( num < 0 || num >= nbnodes )
					{
						std::cout << __FILE__ << ":" << __LINE__ << " node number is out of range: " << num << "/" << nbnodes << std::endl;
						throw BadFile;
					}
					f_nodes.push_back(newnodes[num]);
				}
				if( 1 != fscanf(f," %d",&num) )
				{
					std::cout << __FILE__ << ":" << __LINE__ << " cannot read zone number from " << File << std::endl;
					throw BadFile;
				}
				if( num >= nbzones )
				{
					std::cout << __FILE__ << ":" << __LINE__ << " zone number is out of range: " << num << "/" << nbzones << std::endl;
					throw BadFile;
				}
				Face f = CreateFace(f_nodes).first;
				f->Integer(zone) = num;
				if( num >= 0 ) newsets[num]->PutElement(f);
				newpolygon[i] = f->GetHandle();
				f_nodes.clear();

				if( verbosity > 1 &&  i % report_pace == 0 )
				{
					printf("faces %3.1f%%\r",(i*100.0)/(1.0*nbpolygon));
					fflush(stdout);
				}
			}
			if( verbosity > 0 ) printf("Reading %d cells.\n",nbpolyhedra);
			report_pace = std::max<int>(nbpolyhedra/250,1);
			for(int i = 0; i < nbpolyhedra; i++)
			{
				if( 1 != fscanf(f," %d",&nbpolyhedronfaces) )
				{
					std::cout << __FILE__ << ":" << __LINE__ << " cannot read number of faces from " << File << std::endl;
					throw BadFile;
				}
				for(int j = 0; j < nbpolyhedronfaces; j++)
				{
					if( 1 != fscanf(f," %d", &num) )
					{
						std::cout << __FILE__ << ":" << __LINE__ << " cannot read face number from " << File << std::endl;
						throw BadFile;
					}
					if( num < 0 || num >= nbpolygon )
					{
						std::cout << __FILE__ << ":" << __LINE__ << " face number is out of range: " << num << "/" << nbpolygon << std::endl;
						throw BadFile;
					}
					c_faces.push_back(newpolygon[num]);
				}
				if( 1 != fscanf(f," %d",&num) )
				{
					std::cout << __FILE__ << ":" << __LINE__ << " cannot read face number from " << File << std::endl;
					throw BadFile;
				}
				if( num >= nbzones )
				{
					std::cout << __FILE__ << ":" << __LINE__ << " zone number is out of range: " << num << "/" << nbzones << std::endl;
					throw BadFile;
				}
				Cell c = CreateCell(c_faces).first;
				c->Integer(zone) = num;
				newpolyhedron[i] = c->GetHandle();
				if( num >= 0 ) newsets[num]->PutElement(newpolyhedron[i]);
				c_faces.clear();

				if( verbosity > 1 &&  i % report_pace == 0 )
				{
					printf("cells %3.1f%%\r",(i*100.0)/(1.0*nbpolyhedra));
					fflush(stdout);
				}
			}
			report_pace = std::max<int>(nbzones/250,1);
			if( verbosity > 0 ) printf("Reading %d zones data.\n",nbzones);
			for(int i = 0; i < nbzones; i++)
			{
				if( 4 != fscanf(f," %d %lf %lf %d",&num,&vfac,&poro,&nbK) )
				{
					std::cout << __FILE__ << ":" << __LINE__ << " cannot read zone data from " << File << std::endl;
					throw BadFile;
				}
				if( num < 0 || num >= nbzones )
				{
					std::cout << __FILE__ << ":" << __LINE__ << " zone number out of range " << num << "/" << nbzones << " in file "  << File << std::endl;
					throw BadFile;
				}
				if( !(nbK == 1 || nbK == 2 || nbK == 3 || nbK == 6 || nbK == 9) )
				{
					std::cout << __FILE__ << ":" << __LINE__ << " strange size for permiability tensor: " << nbK << " in file "  << File << std::endl;
					throw BadFile;
				}
				for(int j = 0; j < nbK; j++)
				{
					if( 1 != fscanf(f," %lf",readK+j) )
					{
						std::cout << __FILE__ << ":" << __LINE__ << " cannot read permiability component " << j << "/" << nbK << " from " << File << std::endl;
						throw BadFile;
					}
				}
				memset(K,0,sizeof(Storage::real)*9);
				switch(nbK)
				{
				case 1:
					K[0] = K[4] = K[8] = readK[0];
					break;
				case 2:
					K[0] = readK[0];
					K[4] = readK[1];
					K[8] = 1; //eigenvalue should not be zero
				case 3:
					K[0] = readK[0];
					K[4] = readK[1];
					K[8] = readK[2];
					break;
				case 6:
					K[0] = readK[0];
					K[4] = readK[1];
					K[8] = readK[2];
					K[1] = K[3] = readK[3];
					K[2] = K[6] = readK[4];
					K[5] = K[7] = readK[5];
					break;
				case 9:
					//just copy
					memcpy(K,readK,sizeof(Storage::real)*9);
				}
				for(ElementSet::iterator it = newsets[num]->Begin(); it != newsets[num]->End(); ++it)
				{
					it->Real(volume_factor) = vfac;
					it->Real(porosity) = poro;
					memcpy(it->RealArray(permiability).data(),K,sizeof(Storage::real)*9);
				}
				if( verbosity > 1 &&  i % report_pace == 0 )
				{
					printf("data %3.1f%%\r",(i*100.0)/(1.0*nbzones));
					fflush(stdout);
				}
			}
			fclose(f);
		}
		else if(LFile.find(".msh") != std::string::npos ) // GMSH file format
		{
			Tag gmsh_tags;
			char readline[2048], * p, * pend;
			std::vector<HandleType> newnodes;
			std::vector<HandleType> newelems;
			Storage::real xyz[3];
			int text_start, text_end, state = GMSH_NONE, ver = GMSH_VERNONE;
			int nnodes, ncells, nchars, nodenum, elemnum, elemtype, numtags, elemnodes, temp;
			int ascii, float_size, verlow, verhigh;
			char skip_keyword[2048];
			dynarray<int,128> elemtags;
			dynarray<int,128> nodelist;
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
				case GMSH_NODES_Z:
					if( 1 == sscanf(p,"%lf%n",xyz+2,&nchars) )
					{
						p += nchars;
						while(isspace(*p) && p < pend) ++p;
						int find = -1;
						if( !old_nodes.empty() )
						{
							std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),xyz,CentroidComparator(this));
							if( it != old_nodes.end() ) 
							{
								Storage::real_array c = RealArrayDF(*it,CoordsTag());
								if( CentroidComparator(this).Compare(xyz,c.data()) == 0 )
									find = static_cast<int>(it - old_nodes.begin());
							}
						}
						if( find == -1 ) newnodes[nodenum] = CreateNode(xyz)->GetHandle();
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
						printf("Reading %d elements\n",ncells);
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " error reading total number of elements in " << File << " line " << line <<std::endl;
						throw BadFile;
					}
					if( p >= pend ) break;
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
							break;
						}
						else if( ver == GMSH_VER2 && !strcmp(p,"$EndElements") )
						{
							if( ncells ) std::cout << __FILE__ << ":" << __LINE__ << " number of element records mismatch with reported number of elements " << ncells << " in " << File << " line " << line <<std::endl;
							state = GMSH_NONE;
							break;
						}
						else
						{
							std::cout << __FILE__ << ":" << __LINE__ << " error reading element number in " << File << " got " << p << " instead " << " line " << line <<std::endl;
							throw BadFile;
						}
					}
					if( p >= pend ) break;
				case GMSH_ELEMENTS_TYPE:
					if( 1 == sscanf(p,"%d%n",&elemtype,&nchars) )
					{
						const int type_nodes[19] =
						{
							2,3,4,4,8,6,5,3,6,9,10,
							27,18,14,1,8,20,15,13
						};
						//std::cout << __FILE__ << ":" << __LINE__ << " element type: " << elemtype << " line " << line << std::endl;
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
				case GMSH_ELEMENTS_TAGS:
					if( state == GMSH_ELEMENTS_TAGS )
					{
						while( numtags > 0 && p < pend )
						{
							if( 1 == sscanf(p,"%d%n",&elemtags[elemtags.size()-numtags],&nchars) )
							{
								//std::cout << __FILE__ << ":" << __LINE__ << " tag[" << elemtags.size()-numtags << "] = " << elemtags[elemtags.size()-numtags] << " line " << line << std::endl;
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
								for(dynarray<int,128>::size_type k = 0; k < nodelist.size(); ++k)
									c_nodes.push_back(newnodes[nodelist[k]]);

								//Create new element
								switch(elemtype)
								{
								case 1: newelems[elemnum] = CreateEdge(c_nodes).first->GetHandle(); break; //edge
								case 2: newelems[elemnum] = CreateFace(c_nodes).first->GetHandle(); break; //triangle
								case 3: newelems[elemnum] = CreateFace(c_nodes).first->GetHandle(); break; //quad
								case 4: //tetra
									{
										const integer nodesnum[12] = {0,2,1,0,1,3,1,2,3,0,3,2};
										const integer sizes[4] = {3,3,3,3};
										newelems[elemnum] = CreateCell(c_nodes,nodesnum,sizes,4).first->GetHandle();
									}
									break;
								case 5: //hex
									{
										//const integer nodesnum[24] = {0,4,7,3,1,2,6,5,0,1,5,4,3,7,6,2,0,3,2,1,4,5,6,7};
										const integer nodesnum[24] = {0,3,2,1,4,5,6,7,0,4,7,3,1,2,6,5,0,1,5,4,2,3,7,6};
										const integer sizes[6] = {4,4,4,4,4,4};
										newelems[elemnum] = CreateCell(c_nodes,nodesnum,sizes,6).first->GetHandle();
									}
									break;
								case 6: //prism
									{
										//const integer nodesnum[18] = {0,2,5,3,1,4,5,2,0,3,4,1,3,5,4,0,1,2};
										const integer nodesnum[18] = {0,3,5,2,1,2,5,4,0,1,4,3,3,4,5,0,2,1};
										const integer sizes[5] = {4,4,4,3,3};
										newelems[elemnum] = CreateCell(c_nodes,nodesnum,sizes,5).first->GetHandle();
									}
									break;
								case 7: //pyramid
									{
										const integer nodesnum[16] = {0,4,3,0,1,4,1,2,4,3,4,2,0,3,2,1};
										const integer sizes[5] = {3,3,3,3,4};
										newelems[elemnum] = CreateCell(c_nodes,nodesnum,sizes,5).first->GetHandle();
									}
									break;
								case 8: //second order edge
									{
										//treat as just an edge
										newelems[elemnum] = CreateEdge(ElementArray<Node>(this,c_nodes.data(),c_nodes.data()+2)).first->GetHandle();
									}
									break;
								case 9: //second order tri with nodes on edges
									{
										newelems[elemnum] = CreateFace(ElementArray<Node>(this,c_nodes.data(),c_nodes.data()+3)).first->GetHandle();
									}
									break;
								case 10: //second order quad with nodes on edges and in center
									{
										newelems[elemnum] = CreateFace(ElementArray<Node>(this,c_nodes.data(),c_nodes.data()+4)).first->GetHandle();
									}
									break;
								case 11: // second order tet with nodes on edges
									{
										const integer nodesnum[12] = {0,2,1,0,1,3,1,2,3,0,3,2};
										const integer sizes[4] = {3,3,3,3};
										newelems[elemnum] = CreateCell(c_nodes,nodesnum,sizes,4).first->GetHandle();
									}
									break;
								case 12: // second order hex with nodes in centers of edges and faces and in center of volume
									{
										const integer nodesnum[24] = {0,4,7,3,1,2,6,5,0,1,5,4,3,7,6,2,0,3,2,1,4,5,6,7};
										const integer sizes[6] = {4,4,4,4,4,4};
										newelems[elemnum] = CreateCell(c_nodes,nodesnum,sizes,6).first->GetHandle();
									}
									break;
								case 13: // second order prism  with nodes in centers of edges and quad faces
									{
										const integer nodesnum[18] = {0,3,5,2,1,2,5,4,0,1,4,3,3,4,5,0,2,1};
										const integer sizes[5] = {4,4,4,3,3};
										newelems[elemnum] = CreateCell(c_nodes,nodesnum,sizes,5).first->GetHandle();
									}
									break;
								case 14: // second order pyramid with nodes in centers of edges and quad faces
									{
										const integer nodesnum[16] = {0,4,3,0,1,4,1,2,4,3,4,2,0,3,2,1};
										const integer sizes[5] = {3,3,3,3,4};
										newelems[elemnum] = CreateCell(c_nodes,nodesnum,sizes,5).first->GetHandle();
									}
									break;
								case 15: // vertex
									{
										newelems[elemnum] = c_nodes.at(0);
									}
									break;
								case 16: // second order quad with nodes on edges
									{
										newelems[elemnum] = CreateFace(ElementArray<Node>(this,c_nodes.data(),c_nodes.data()+4)).first->GetHandle();
									}
									break;
								case 17: // second order hex with nodes on edges
									{
										const integer nodesnum[24] = {0,4,7,3,1,2,6,5,0,1,5,4,3,7,6,2,0,3,2,1,4,5,6,7};
										const integer sizes[6] = {4,4,4,4,4,4};
										newelems[elemnum] = CreateCell(c_nodes,nodesnum,sizes,6).first->GetHandle();
									}
									break;
								case 18: // second order prism with nodes on edges
									{
										const integer nodesnum[18] = {0,3,5,2,1,2,5,4,0,1,4,3,3,4,5,0,2,1};
										const integer sizes[5] = {4,4,4,3,3};
										newelems[elemnum] = CreateCell(c_nodes,nodesnum,sizes,5).first->GetHandle();
									}
									break;
								case 19: // second order pyramid with nodes on edges
									{
										const integer nodesnum[16] = {0,4,3,0,1,4,1,2,4,3,4,2,0,3,2,1};
										const integer sizes[5] = {3,3,3,3,4};
										newelems[elemnum] = CreateCell(c_nodes,nodesnum,sizes,5).first->GetHandle();
									}
									break;
								}

								Storage::integer_array ptags = IntegerArray(newelems[elemnum],gmsh_tags);

								if( ver == GMSH_VER2 )
								{
									ptags.replace(ptags.begin(),ptags.end(),elemtags.begin(),elemtags.end());
									elemtags.clear();
								}
								else
								{
									ptags[0] = elemtags[0];
									ptags[1] = elemtags[1];
								}

								c_nodes.clear();
								state = GMSH_ELEMENTS_NUM;
								ncells--;
								if( verbosity > 1 && (newelems.size()-ncells)%report_pace == 0 )
								{
									printf("%3.1f%%\r",((newelems.size()-ncells)*100.0)/(newelems.size()*1.0));
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
		else if(LFile.find(".pmf") != std::string::npos) //this is inner parallel/platform mesh format
		{
			io_converter<INMOST_DATA_INTEGER_TYPE,INMOST_DATA_REAL_TYPE> iconv;
			io_converter<INMOST_DATA_ENUM_TYPE   ,INMOST_DATA_REAL_TYPE> uconv;
			REPORT_STR("start load pmf");
			dynarray<INMOST_DATA_ENUM_TYPE,128> myprocs;
			std::stringstream in(std::ios::in | std::ios::out | std::ios::binary);
			HeaderType token;
#if defined(USE_MPI)
			if( m_state == Mesh::Parallel )
			{
#if defined(USE_MPI_FILE)
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
						ierr = MPI_File_read_shared(fh,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),MPI_CHAR,&stat);
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);

						header << buffer;
						uconv.read_iValue(header,datanum);

						buffer.resize(datanum*uconv.get_source_iByteSize());
						std::vector<INMOST_DATA_ENUM_TYPE> datasizes(datanum);
						//ierr = MPI_File_read_all(fh,&buffer[0],buffer.size(),MPI_CHAR,&stat);
						ierr = MPI_File_read_shared(fh,&buffer[0],static_cast<INMOST_MPI_SIZE>(buffer.size()),MPI_CHAR,&stat);
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
					else recvsizes.resize(1,0); //protect from dereferencing null

					ierr = MPI_Scatter(&recvsizes[0],1,INMOST_MPI_DATA_ENUM_TYPE,&recvsize,1,INMOST_MPI_DATA_ENUM_TYPE,0,GetCommunicator());
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);


					buffer.resize(std::max(1u,recvsize)); //protect from dereferencing null


					{
						ierr = MPI_File_read_ordered(fh,&buffer[0],static_cast<INMOST_MPI_SIZE>(recvsize),MPI_CHAR,&stat);
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
					std::vector<INMOST_MPI_SIZE> sendcnts(numprocs), displs(numprocs);
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
						sendcnts[0] = static_cast<INMOST_MPI_SIZE>(recvsizes[0]);
						for(k = 1; k < numprocs; k++)
						{
							sendcnts[k] = static_cast<INMOST_MPI_SIZE>(recvsizes[k]);
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
						INMOST_DATA_ENUM_TYPE datalength;
						uconv.read_iValue(in,namesize);
						in.read(name, namesize);
						assert(namesize < 4096);
						name[namesize] = '\0';
						in.get(datatype);
						in.get(sparsemask);
						in.get(definedmask);
						uconv.read_iValue(in,datalength);
						tags[i] = CreateTag(std::string(name),static_cast<DataType>(datatype),
						                    static_cast<ElementType>(definedmask),
						                    static_cast<ElementType>(sparsemask),datalength);
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
						new_nodes.data(),
						new_edges.data(),
						new_faces.data(),
						new_cells.data()
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
						new_nodes.data(),
						new_edges.data(),
						new_faces.data(),
						new_cells.data(),
						new_sets.data(),
						&m_storage
					};
					for(INMOST_DATA_ENUM_TYPE j = 0; j < tags.size(); j++) 
					{
						REPORT_VAL("TagName",tags[j].GetTagName());
						Tag * jt = &tags[j];
						for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
						{
							if( jt->isDefined(etype) )
							{
								INMOST_DATA_ENUM_TYPE q, cycle_end, etypenum = ElementNum(etype);
								cycle_end = elem_sizes[etypenum];
								bool sparse = jt->isSparse(etype);
								INMOST_DATA_ENUM_TYPE tagsize = jt->GetSize(), recsize = tagsize, lid;
								INMOST_DATA_ENUM_TYPE k;
								DataType data_type = jt->GetDataType();
								if( sparse ) uconv.read_iValue(in,q); else q = 0;
								while(q != cycle_end)
								{
									HandleType he = elem_links[etypenum][q];
									if( tagsize == ENUMUNDEF ) 
									{
										uconv.read_iValue(in,recsize); 
										SetDataSize(he,*jt,recsize);
									}
									switch(data_type)
									{
									case DATA_REAL:
										{
											Storage::real_array arr = RealArray(he,*jt);
											for(k = 0; k < recsize; k++) 
												iconv.read_fValue(in,arr[k]); 
										} break;
									case DATA_INTEGER:   
										{
											Storage::integer_array arr = IntegerArray(he,*jt); 
											for(k = 0; k < recsize; k++) 
												iconv.read_iValue(in,arr[k]); 
										} break;
									case DATA_BULK:      
										{
											in.read(reinterpret_cast<char *>(&Bulk(he,*jt)),recsize);
										} break;
									case DATA_REFERENCE: 
										{
											Storage::reference_array arr = ReferenceArray(he,*jt);
											for(k = 0; k < recsize; k++)
											{
												char type;
												in.get(type);
												if (type != NONE)
												{
													uconv.read_iValue(in, lid);
													arr.at(k) = elem_links[ElementNum(type)][lid];
												}
												else arr.at(k) = InvalidHandle();
											}
										} break;
									}
									if( sparse ) uconv.read_iValue(in,q); else q++;
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
				INMOST_DATA_ENUM_TYPE numprocs = GetProcessorsNumber(), size = static_cast<INMOST_DATA_ENUM_TYPE>(myprocs.size()),k;
				INMOST_DATA_ENUM_TYPE procs_sum = 0;
				std::vector<INMOST_DATA_ENUM_TYPE> procs_sizes(numprocs);
				REPORT_MPI(MPI_Allgather(&size,1,INMOST_MPI_DATA_ENUM_TYPE,procs_sizes.data(),1,INMOST_MPI_DATA_ENUM_TYPE,GetCommunicator()));
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
					REPORT_MPI(ierr = MPI_Allgatherv(myprocs.data(),procs_sizes[myrank],INMOST_MPI_DATA_ENUM_TYPE,procs.data(),recvcnts.data(),displs.data(),INMOST_MPI_DATA_ENUM_TYPE,GetCommunicator()));
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
				ComputeSharedProcs();
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
		else if(LFile.find(".vtk") != std::string::npos) //this is legacy vtk
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
					if( t.isDefined(CELL) && !t.isSparse(CELL) && t.GetDataType() != DATA_BULK && t.GetDataType() != DATA_REFERENCE &&
					   t != CoordsTag() && t != SharedTag() && t != SendtoTag() && t != ProcessorsTag())
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
							fprintf(f,"SCALARS %s %s %d\n",tags[i].GetTagName().c_str(),(tags[i].GetDataType() == DATA_REAL ? "double" : "int"),comps);
							fprintf(f,"LOOKUP_TABLE default\n");
							for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); it++)
							{
								switch( tags[i].GetDataType() )
								{
									case DATA_REAL:
									{
										Storage::real_array arr = it->RealArray(tags[i]);
										for(unsigned int m = 0; m < comps; m++) fprintf(f,"%14e ",arr[m]);
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
							fprintf(f,"SCALARS %s %s %d\n",tags[i].GetTagName().c_str(),(tags[i].GetDataType() == DATA_REAL ? "double" : "int"),comps);
							fprintf(f,"LOOKUP_TABLE default\n");
							for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); it++)
							{
								switch( tags[i].GetDataType() )
								{
									case DATA_REAL:
									{
										Storage::real_array arr = it->RealArray(tags[i]);
										for(unsigned int m = 0; m < comps; m++) fprintf(f,"%14e ",arr[m]);
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
			//ReorderEmpty(CELL | FACE | NODE | ESET);
			FILE * file = fopen(File.c_str(),"wb");
			char keyword[2048];
			sprintf(keyword,"gmvinput"); fwrite(keyword,1,8,file);
			sprintf(keyword,"ieee");
			if( sizeof(Storage::real) != 8 || sizeof(Storage::integer) != 8 )
				sprintf(keyword,"ieeei%ldr%ld",sizeof(Storage::integer),sizeof(Storage::real));
			fwrite(keyword,1,8,file);
			sprintf(keyword,"nodev"); fwrite(keyword,1,8,file);
			keynum = static_cast<Storage::integer>(NumberOfNodes()); fwrite(&keynum,sizeof(Storage::integer),1,file);
			for(Mesh::iteratorNode n = BeginNode(); n != EndNode(); n++)
			{
				fwrite(n->Coords().data(),sizeof(Storage::real),GetDimensions(),file);
				if( GetDimensions() < 3 )
				{
					Storage::real zero = 0.0;
					fwrite(&zero,sizeof(Storage::real),3-GetDimensions(),file);
				}
			}
			sprintf(keyword,"faces"); fwrite(keyword,1,8,file);
			keynum = static_cast<Storage::integer>(NumberOfFaces()); fwrite(&keynum,sizeof(Storage::integer),1,file);
			keynum = static_cast<Storage::integer>(NumberOfCells()); fwrite(&keynum,sizeof(Storage::integer),1,file);
			Tag set_id = CreateTag("TEMPORARY_ELEMENT_ID",DATA_INTEGER,CELL|FACE|NODE,NONE,1);
			Storage::integer cur_num = 0;
			for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it) it->IntegerDF(set_id) = cur_num++;
			cur_num = 0;
			for(Mesh::iteratorFace it = BeginFace(); it != EndFace(); ++it) it->IntegerDF(set_id) = cur_num++;
			cur_num = 0;
			for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); ++it) it->IntegerDF(set_id) = cur_num++;
			for(Mesh::iteratorFace f = BeginFace(); f != EndFace(); f++)
			{
				ElementArray<Node> fnodes = f->getNodes();
				//sprintf(keyword,"nverts"); fwrite(keyword,1,8,file);
				keynum = static_cast<Storage::integer>(fnodes.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
				//sprintf(keyword,"vertex_ids"); fwrite(keyword,1,8,file);
				for(ElementArray<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
				{
					keynum = fn->IntegerDF(set_id)+1;
					fwrite(&keynum,sizeof(Storage::integer),1,file);
				}
				ElementArray<Cell> fcells = f->getCells();
				//sprintf(keyword,"cellno1"); fwrite(keyword,1,8,file);
				keynum = fcells.size() > 0 ? fcells[0].IntegerDF(set_id)+1 : 0;
				fwrite(&keynum,sizeof(Storage::integer),1,file);
				//sprintf(keyword,"cellno2"); fwrite(keyword,1,8,file);
				keynum = fcells.size() > 1 ? fcells[1].IntegerDF(set_id)+1 : 0;
				fwrite(&keynum,sizeof(Storage::integer),1,file);
			}
			
			
			if( HaveTag("MATERIAL") )
			{
				Tag t = GetTag("MATERIAL");
				if( t.isDefined(CELL) )
				{
					std::set<Storage::integer> mat;
					for(Mesh::iteratorCell c = BeginCell(); c != EndCell(); c++)
						mat.insert(c->Integer(t)+1);
					sprintf(keyword,"material"); fwrite(keyword,1,8,file);
					keynum = static_cast<Storage::integer>(mat.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
					keynum = 0; fwrite(&keynum,sizeof(Storage::integer),1,file);
					for(std::set<Storage::integer>::iterator it = mat.begin(); it != mat.end(); it++)
					{
						sprintf(keyword,"mat%d",*it); fwrite(keyword,1,8,file);
					}
					for(Mesh::iteratorCell c = BeginCell(); c != EndCell(); c++)
					{
						keynum = c->Integer(t)+1; fwrite(&keynum,sizeof(Storage::integer),1,file);
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
									keynum = e->IntegerDF(set_id)+1;
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
			for(Mesh::iteratorSet set = BeginSet(); set != EndSet(); set++)
			{

				for(ElementType etype = CELL; etype >= NODE; etype = etype >> 1)
				{
					if( etype == EDGE ) continue;
					keynum = 0;
					for(ElementSet::iterator it = set->Begin(); it != set->End(); it++)
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
						//sprintf(keyword,"set%d_%s\n",set->IntegerDF(set_id)+1,ElementTypeName(etype));
						sprintf(keyword,"%s",set->GetName().c_str());
						fwrite(keyword,1,8,file);
						fwrite(&temp,sizeof(Storage::integer),1,file);						
						fwrite(&keynum,sizeof(Storage::integer),1,file);
						for(ElementSet::iterator it = set->Begin(); it != set->End(); it++)
							if( it->GetElementType() == etype )
							{
								keynum = it->IntegerDF(set_id)+1;
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
						
						for(INMOST_DATA_ENUM_TYPE q = 0; q < t->GetSize(); q++)
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
			for(Mesh::iteratorFace f = BeginFace(); f != EndFace(); f++) if( f->Boundary() )
			{
				ElementArray<Node> fnodes = f->getNodes();
				keynum = 1; fwrite(&keynum,sizeof(Storage::integer),1,file);
				keynum = static_cast<Storage::integer>(fnodes.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
				for(ElementArray<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
					fwrite(&fn->Coords()[0],sizeof(Storage::real),1,file);
				for(ElementArray<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
					fwrite(&fn->Coords()[1],sizeof(Storage::real),1,file);
				if( GetDimensions() > 2 )
				{
					for(ElementArray<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
						fwrite(&fn->Coords()[2],sizeof(Storage::real),1,file);
				}
				else 
				{
					Storage::real zero = 0.0;
					fwrite(&zero,sizeof(Storage::real),fnodes.size(),file);
				}
			}
			sprintf(keyword,"endpoly"); fwrite(keyword,1,8,file);
			sprintf(keyword,"endgmv"); fwrite(keyword,1,8,file);
			fclose(file);
		}
		else if(LFile.find(".pmf") != std::string::npos) //this is inner parallel/platform mesh format
		{
			io_converter<INMOST_DATA_INTEGER_TYPE,INMOST_DATA_REAL_TYPE> iconv;
			io_converter<INMOST_DATA_ENUM_TYPE   ,INMOST_DATA_REAL_TYPE> uconv;
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

			INMOST_DATA_ENUM_TYPE header[9] = {GetDimensions(), NumberOfNodes(), NumberOfEdges(), NumberOfFaces(), NumberOfCells(), NumberOfSets(), NumberOfTags()-3, m_state, GetProcessorRank()},k;
			for(k = 0; k < 9; k++) uconv.write_iValue(out,header[k]);
			out.write(reinterpret_cast<char *>(remember),sizeof(remember));
		

			Tag set_id = CreateTag("TEMPORARY_ELEMENT_ID",DATA_INTEGER,ESET|CELL|FACE|EDGE|NODE,NONE,1);
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

			// Tags
			out << INMOST::TagsHeader;
			uconv.write_iValue(out,header[6]);
			REPORT_STR("TagsHeader");
			REPORT_VAL("tag_size",header[6]);
			for(Mesh::iteratorTag it = BeginTag(); it != EndTag(); it++)  
			{
				//don't forget to change header[6] if you skip more
				//should match with content after MeshDataHeader
				if( *it == set_id ) continue;
				if( *it == HighConnTag() ) continue;
				if( *it == LowConnTag() ) continue;

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
						lid = IntegerDF(*kt,set_id);
						uconv.write_iValue(out,lid);
					}
					else out.put(NONE);
				}
				//write additional information
				for(Element::adj_type::iterator kt = hc.begin()+ElementSet::high_conn_reserved-1; kt != hc.end(); ++kt)
				{
					if( *kt != InvalidHandle() )
					{
						wetype = GetHandleElementType(*kt);
						out.put(wetype);
						lid = IntegerDF(*kt,set_id);
						uconv.write_iValue(out,lid);
					}
					else out.put(NONE);
				}
			}
		
			out << INMOST::MeshDataHeader;
			REPORT_STR("MeshDataHeader");

			for(Mesh::iteratorTag jt = BeginTag(); jt != EndTag(); jt++) 
			{
				std::string tagname = jt->GetTagName();
				//skipping should match with header[6] and content
				// after TagsHeader
				if( *jt == set_id ) continue;
				if( *jt == HighConnTag() ) continue;
				if( *jt == LowConnTag() ) continue;
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
				INMOST_DATA_ENUM_TYPE numprocs = GetProcessorsNumber(), datasize = static_cast<INMOST_DATA_ENUM_TYPE>(out.tellp()),k;
				std::vector<INMOST_DATA_ENUM_TYPE> datasizes(numprocs,0);
				REPORT_VAL("local_write_file_size",datasize);
				REPORT_MPI(ierr = MPI_Gather(&datasize,1,INMOST_MPI_DATA_ENUM_TYPE,&datasizes[0],1,INMOST_MPI_DATA_ENUM_TYPE,0,GetCommunicator()));
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
#if defined(USE_MPI_FILE) //We have access to MPI_File
				if( parallel_file_strategy == 1 )
				{
					MPI_File fh;
					MPI_Status stat;
					REPORT_MPI(ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(File.c_str()),MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh));
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
						REPORT_MPI(ierr = MPI_File_write_shared(fh,&header_data[0],static_cast<INMOST_MPI_SIZE>(header_data.size()),MPI_CHAR,&stat));
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
					}
					{
						std::string local_data(out.str());
						REPORT_MPI(ierr = MPI_File_write_ordered(fh,&local_data[0],static_cast<INMOST_MPI_SIZE>(local_data.size()),MPI_CHAR,&stat));
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
					}

					REPORT_MPI(ierr = MPI_File_close(&fh));
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
				}
				else
#endif
				{
					std::vector<INMOST_MPI_SIZE> displs(numprocs),recvcounts(numprocs);
					std::string file_contents;
					std::string local_data(out.str());
					if( GetProcessorRank() == 0 )
					{
						recvcounts[0] = static_cast<INMOST_MPI_SIZE>(datasizes[0]);
						displs[0] = 0;
						int sizesum = recvcounts[0];
						for(k = 1; k < numprocs; k++)
						{
							recvcounts[k] = static_cast<INMOST_MPI_SIZE>(datasizes[k]);
							displs[k] = displs[k-1]+recvcounts[k-1];
							sizesum += recvcounts[k];
						}
						file_contents.resize(sizesum);
					}
					else file_contents.resize(1); //protect from accessing bad pointer
					REPORT_MPI(ierr = MPI_Gatherv(&local_data[0],static_cast<INMOST_MPI_SIZE>(local_data.size()),MPI_CHAR,&file_contents[0],&recvcounts[0],&displs[0],MPI_CHAR,0,GetCommunicator()));
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
		if(LFile.find(".msh") != std::string::npos) return false;
		if(LFile.find(".grid") != std::string::npos) return false;
		else if(LFile.find(".vtk") != std::string::npos) return false;
		else if(LFile.find(".gmv") != std::string::npos) return false;
		else if(LFile.find(".pvtk") != std::string::npos) return true;
		else if(LFile.find(".pmf") != std::string::npos) return true;
		throw NotImplemented;
	}
}
#endif
