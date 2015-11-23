#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "inmost.h"

#if defined(USE_MESH)

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


namespace INMOST
{
  void Mesh::LoadECL(std::string File)
  {
		FILE * f = fopen(File.c_str(),"r");
		if( f == NULL )
		{
			std::cout << __FILE__ << ":" << __LINE__ << " cannot open file " << File << std::endl;
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
						read_arrayf = xyz.empty()? NULL : &xyz[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = state = ECL_DX;
					}
					else if( !ECLSTRCMP(p,"DY") )
					{
						assert(have_dimens);
						if( xyz.empty() ) xyz.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = xyz.empty() ? NULL : &xyz[0];
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
						read_arrayf = xyz.empty() ? NULL : &xyz[0];
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
						read_arrayf = xyz.empty() ? NULL : &xyz[0];
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
						read_arrayf = zcorn.empty() ? NULL : &zcorn[0];
						numrecs = 1;
						downread = totread = dims[0]*dims[1]*dims[2]*8;
						argtype = ECL_VAR_REAL;
						state = offset = ECL_ZCORN;
					}
					else if( !ECLSTRCMP(p,"TOPS") )
					{
						assert(have_dimens);
						tops.resize(dims[0]*dims[1]*dims[2]);
						read_arrayf = tops.empty() ? NULL : &tops[0];
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
						read_arrayf = perm.empty() ? NULL : &perm[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = state = ECL_PERMX;
					}
					else if( !ECLSTRCMP(p,"PERMY") )
					{
						assert(have_dimens);
						if( perm.empty() ) perm.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = perm.empty() ? NULL : &perm[0];
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
						read_arrayf = perm.empty() ? NULL : &perm[0];
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
						read_arrayf = poro.empty() ? NULL : &poro[0];
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
}

#endif