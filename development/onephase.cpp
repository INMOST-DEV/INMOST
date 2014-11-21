#include "dlg.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <time.h>
#if defined(_WIN32)
#include <direct.h>
#else
#include <unistd.h>
#endif
#include <dirent.h> //for win32 - download from http://softagalleria.net/download/dirent/?C=M;O=D


const bool nmpfa_special = false;

using namespace INMOST;
INMOST_DATA_ENUM_TYPE wnum = 0;

void get_well_func(Storage * current_element, Automatizator::stencil_pairs & out_stencil,void * user_data)
{
	out_stencil.push_back(Automatizator::stencil_pair(((wells *)user_data)->GetWellElements(wnum),1.0));
}

void add_row(Solver::Row & Ar, Solver::Row & r, Storage::real mult)
{
	for(Solver::Row::iterator it = r.Begin(); it != r.End(); ++it) if( fabs(it->second) > 0.0 )
	{
		Ar[it->first] += mult*it->second;
	}
}

const int name_width = 26;
const int type_width = 14;
const int elems_width = 10;
const int sparse_width = 10;
const int length_width = 10;


void PrintTag(Tag t)
{
	std::cout << std::setw(name_width) << t.GetTagName() << std::setw(type_width) << DataTypeName(t.GetDataType());
	int num = 0;
	char elems[7] = "NEFCSM";
	std::string print = "";
	for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
	{
		if( t.isDefined(etype) )
		{
			print = print + elems[ElementNum(etype)];
			num++;
		}
	}
	std::cout << std::setw(elems_width) << print;
	print = "";
	for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
	{
		if( t.isSparse(etype) )
		{
			print = print + elems[ElementNum(etype)];
			num++;
		}
	}
	std::cout << std::setw(sparse_width) << print;
	if( t.GetSize() == ENUMUNDEF )
		std::cout << std::setw(length_width) << "VAR" << std::endl;
	else
		std::cout << std::setw(length_width) << t.GetSize() << std::endl;
}

void PrintTags(Mesh * m, ElementType etypes = MESH|ESET|CELL|FACE|EDGE|NODE)
{
	std::cout << std::setw(name_width) << "Name" << std::setw(type_width) << "Type" << std::setw(elems_width) << "Element" << std::setw(sparse_width) << "Sparse" << std::setw(length_width) << "Length" << std::endl;
	for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); ++t )
	{
		bool print = false;
		for(ElementType etype = NODE; etype <= MESH; etype = etype << 1) if( (etype&etypes) && t->isDefined(etype) ) {print = true; break;}
		if( print ) PrintTag(*t);
	}
}

void PrintRealTags(Mesh * m, ElementType etypes = MESH|ESET|CELL|FACE|EDGE|NODE)
{
	std::cout << std::setw(name_width) << "Name" << std::setw(type_width) << "Type" << std::setw(elems_width) << "Element" << std::setw(sparse_width) << "Sparse" << std::setw(length_width) << "Length" << std::endl;
	for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); ++t ) if( t->GetDataType() == DATA_REAL )
	{
		bool print = false;
		for(ElementType etype = NODE; etype <= MESH; etype = etype << 1) if( (etype&etypes) && t->isDefined(etype) ) {print = true; break;}
		if( print ) PrintTag(*t);
	}
}

void PrintNumericTags(Mesh * m, ElementType etypes = MESH|ESET|CELL|FACE|EDGE|NODE)
{
	std::cout << std::setw(name_width) << "Name" << std::setw(type_width) << "Type" << std::setw(elems_width) << "Element" << std::setw(sparse_width) << "Sparse" << std::setw(length_width) << "Length" << std::endl;
	for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); ++t ) if( t->GetDataType() != DATA_REFERENCE )
	{
		bool print = false;
		for(ElementType etype = NODE; etype <= MESH; etype = etype << 1) if( (etype&etypes) && t->isDefined(etype) ) {print = true; break;}
		if( print ) PrintTag(*t);
	}
}

std::string stolower(const std::string & input)
{
	std::string ret(input);
	for(size_t k = 0; k < ret.size(); k++) ret[k] = tolower(ret[k]);
	return ret;
}

bool question_yesno()
{
	while(true)
	{
		std::string ans;
		std::cin >> ans;
		std::string lans = stolower(ans);
		if( lans == "yes" || lans == "y" )
			return true;
		else if ( lans == "no" || lans == "n")
			return false;
		else 
		{
			std::cout << "I don't understand " << ans << "." << std::endl;
			std::cout << "Please enter \"yes\" or \"no\" again." << std::endl;
		}
	}
}

int question_keyword(const std::string keywords[], const int N)
{
	while(true)
	{
		std::string ans;
		std::cin >> ans;
		std::string lans = stolower(ans);
		for(int k = 0; k < N; k++)
		{
			std::string lkeyword = stolower(keywords[k]);
			if( lans == lkeyword )
				return k;
		}
		std::cout << "I don't understand " << ans << "." << std::endl;
		std::cout << "Please enter one of:" << std::endl;
		for(int k = 0; k < N; k++)
			std::cout << keywords[k] << ", ";
		std::cout << std::endl;
	}
}

void _list_directory(const char *dirname)
{
	DIR *dir;
	struct dirent *ent;
	dir = opendir (dirname);
	if (dir != NULL) 
	{
		const size_t maxw = 24, maxl = 86;
		size_t line = 0, width, recs;
		while ((ent = readdir (dir)) != NULL) 
		{
			width = strlen(ent->d_name);
			switch (ent->d_type) 
			{
			case DT_DIR: width+=2; break;
			default: width = 0;
			}

			recs = ceil(width/(maxw*1.0));
			if( line+recs*maxw > maxl )
			{
				std::cout << std::endl;
				line = 0;
			}
			line += recs*maxw;

			switch (ent->d_type) 
			{
			case DT_DIR: std::cout << std::setw(maxw*recs) << (std::string(ent->d_name)+"/"); break;
			}
		}
		rewinddir(dir);
		std::string formats[6]= {".vtk",".pvtk",".msh",".grid",".grdecl",".pmf"};
		while ((ent = readdir (dir)) != NULL) 
		{
			std::string fname(ent->d_name);
			for(size_t k = 0; k < fname.size(); k++) fname[k] = tolower(fname[k]);
			bool print = false;
			for(size_t k = 0; k < 6; k++)
			{
				if( fname.find(formats[k]) != std::string::npos )
				{
					print = true;
					break;
				}
			}
			if( !print ) continue;

			width = strlen(ent->d_name);
			switch (ent->d_type) 
			{
			case DT_DIR: width = 0; break;
			case DT_REG: width++; break;
			case DT_LNK: width+=2; break;
			default: width+=2;
			}

			recs = ceil(width/(maxw*1.0));
			if( line+recs*maxw > maxl )
			{
				std::cout << std::endl;
				line = 0;
			}
			line += recs*maxw;

			switch (ent->d_type) 
			{
			case DT_DIR: break;
			case DT_REG: std::cout << std::setw(recs*maxw) << ent->d_name; break;
			case DT_LNK: std::cout << std::setw(recs*maxw) << (std::string(ent->d_name)+"@"); break;
			default: std::cout << std::setw(recs*maxw) << (std::string(ent->d_name)+"*");
			}
		}
		std::cout << std::endl;
		closedir (dir);
	}
	else std::cout << "Cannot open directory" << dirname << std::endl;
}

struct intpair
{
	int vals[2];
	int & operator [](int k) {assert(k >= 0 && k < 2); return vals[k];}
	intpair & operator =(intpair const & b) {vals[0] = b.vals[0]; vals[1] = b.vals[1]; return *this;}
	intpair(const intpair & b){vals[0] = b.vals[0]; vals[1] = b.vals[1];}
	intpair() {vals[0] = vals[1] = 0;}
	intpair(int a, int b) {vals[0] = a; vals[1] = b;}
};

void mesh_load_help_message()
{
	char cwd[2048];
	getcwd(cwd,2048);
	std::cout << "Currently supported file formats: " << std::endl;
	std::cout << ".vtk    - legacy vtk file with unstructured mesh description;" << std::endl;
	std::cout << ".pvtk   - set of legacy vtk files with unstructured mesh" << std::endl;
	std::cout << "          description for parallel computation;" << std::endl;
	std::cout << ".msh    - gmsh mesh generator file format;" << std::endl;
	std::cout << ".grid   - Mohammad Karimi-Fard format of discrete feature modelling;" << std::endl;
	std::cout << ".grdecl - eclipse or petrel file format (not complete);" << std::endl;
	std::cout << ".pmf    - inner INMOST portable parallel binary file format." << std::endl;
	std::cout << "Type " << std::endl;
	std::cout << "  \"ls\" or \"dir\" to list dirrectory;" << std::endl;
	std::cout << "  \"cd\" to change current directory;" << std::endl;
	std::cout << "  \"help\" to reproduce this message." << std::endl;
	std::cout << "Current directory: " << std::endl << cwd << std::endl;
	std::cout << "Please enter the path to the mesh file in one of the supported formats. " << std::endl;
}

void permeability_setup_help_message()
{
	std::cout << "To generate discretization I need permeability tensor " << std::endl;
	std::cout << "to be defined on the grid." << std::endl;
	std::cout << "Availible options to define permeability: " << std::endl;
	std::cout << "data     - you have a data tag that define the tensor on the grid;" << std::endl;
	std::cout << "zones    - you have a data tag that define zones on the grid by" << std::endl;
	std::cout << "           which you can enter permeability tensor;" << std::endl;
	std::cout << "input    - you will input bounding boxes that cover portion of" << std::endl;
	std::cout << "           the grid and let you define permeability on this portion;" << std::endl;
	std::cout << "random2d - I will put the symmetric positive definitive tensor" << std::endl;
	std::cout << "           with 2D variation, Kz is fixed to 1;" << std::endl;
	std::cout << "random3d - I will put the symmetric positive definitive tensor" << std::endl;
	std::cout << "           with 3D variation." << std::endl;
	std::cout << "How would you like to define tensor on the grid? Please enter the keyword." << std::endl;
}

void scheme_setup_help_message()
{
	std::cout << "Availible discretization schemes: " << std::endl;
	std::cout << "TPFA    - conventional two-point flux approximation;" << std::endl;
	std::cout << "NTPFA-A - new nonlinear two-point scheme;" << std::endl;
	std::cout << "NMPFA-A - new nonlinear multi-point scheme," << std::endl;
	std::cout << "          respect discrete maximum principle;" << std::endl;
	std::cout << "MPFA-A  - new scheme derived from NMPFA-A by fixing weights;" << std::endl;
	std::cout << "NTPFA-B - improved nonlinear scheme from paper of " << std::endl;
	std::cout << "          Danilov,Vassilevsky, 2010;" << std::endl;
	std::cout << "NMPFA-B - nonlinear scheme from paper of Lipnikov," << std::endl;
	std::cout << "          Svyatsky, Vassilevsky, 2012;" << std::endl;
	std::cout << "MPFA-B  - derived from NMPFA-B by fixing weights;" << std::endl;
	std::cout << "MPFA-O  - hybrid variational finite volume formulation" <<std::endl;
	std::cout << "          from paper of Agelas, Masson, 2008;" << std::endl;
	std::cout << "MPFA-L  - realisation from paper of Aavatsmark, Eigestad," << std::endl;
	std::cout << "          Mallison, Nordbotten, Oian, 2005;" << std::endl;
	std::cout << "MPFA-G  - realisation from paper of Agelas, Di Pietro," << std::endl;
	std::cout << "          Droniou, 2008;" << std::endl;
	std::cout << "LDMP    - new scheme under construction." << std::endl;
	std::cout << "Please enter one of the scheme names." << std::endl;
}

int main(int argc, char ** argv)
{
	/*
	if (argc < 4)
	{
		std::cout << "Usage: " << argv[0] << " mesh scheme tensor [wells.txt=""] [allow_harmonic = 1] [allow_neumann = 1] [trans_scale = 1] [vol_scale = 1] [well_index_scale = 1] [x_scale = 1] [y_scale = 1] [z_scale = 1] [x_shift = 0] [y_shift = 0] [z_shift = 0]" << std::endl;
		return -1;
	}
	*/
	Mesh::Initialize(&argc,&argv);
	Solver::Initialize(&argc,&argv);
	std::string wellsfile = "";
	Storage::real trans_scale = 1.0, vol_scale = 1.0, wi_scale = 1.0;
	
	int Ktype = 0;
	int allow_harmonic = AVG_NONLINEAR | AVG_REVDIST;
	Storage::real Kxmin, Kxmax, Kymin, Kymax, Phimin, Phimax;
	int tensor_seed = 0;
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	Automatizator * aut = new Automatizator(m);
	discr_basic * d;

	std::stringstream record_dialog;

	//std::cout << argv[0] << std::endl;
	//system("echo %CD%");

	//m->SetCommunicator(INMOST_MPI_COMM_WORLD);

	char cwd[2048];
	getcwd(cwd,2048);

	record_dialog << "MESH" << std::endl;

	do
	{
		char meshfilestr[2048];
		
		mesh_load_help_message();
		
		bool loaded = false, asked2d = false;
		while( !loaded )
		{
			std::cin.getline(meshfilestr,2048);
			size_t linesize = strlen(meshfilestr);
			if( linesize == 0 ) continue;
			if( !strncmp(meshfilestr,"cd ",std::min<size_t>(linesize,3)) )
			{
				record_dialog << "cd " << (linesize > 3 ? meshfilestr+3 : ".") << std::endl;
				chdir(linesize > 3 ? meshfilestr+3 : ".");
			}
			else if( !strncmp(meshfilestr,"ls ",std::min<size_t>(linesize,3)) || 
				     !strncmp(meshfilestr,"dir ",std::min<size_t>(linesize,4)))
				_list_directory(linesize > 3 ? meshfilestr+3 : ".");
			else if( !strncmp(meshfilestr,"help ",std::min<size_t>(linesize,5)) )
				mesh_load_help_message();
			else
			{
				std::string meshfile(meshfilestr);
				int dims = 3;
				if( !asked2d && (meshfile.find(".vtk") != std::string::npos || meshfile.find(".pvtk") != std::string::npos) )
				{
					std::cout << "Should " << meshfile << " be treated as 2D or 3D mesh? Please enter number of dimensions. " << std::endl;
					std::string keywords[] = {"2","3"};
					dims = question_keyword(keywords,2);
					if( dims == 0 ) m->SetFileOption("VTK_GRID_DIMS","2");
					else m->SetFileOption("VTK_GRID_DIMS","3");
					asked2d = true;
				}
				try
				{
					if( m->isParallelFileFormat(meshfile) ) m->Load(meshfile);
					else if( m->GetProcessorRank() == 0 ) m->Load(meshfile);	
					loaded = true;
					record_dialog << meshfile;
					if( dims != 3 ) record_dialog << " " << dims;
					record_dialog << std::endl;
				}
				catch(...)
				{
					std::cout << "Something went wrong with: " << std::endl;
					std::cout << meshfile << std::endl;
					std::cout << "Please enter another path." << std::endl;
				}
			}
		}
		std::cout << "Would you like to load another mesh? " << std::endl;
		std::cout << "New mesh will be merged with existing if nodes match." << std::endl;
		std::cout << "Please enter \"yes\" or \"no\"." << std::endl;
	} while( question_yesno() );

	record_dialog << "/" << std::endl;

	chdir(cwd);

	record_dialog << "SHIFT" << std::endl;

	{
		std::cout << "Would you like to rescale or shift your mesh?" << std::endl;
		std::cout << "Please enter \"yes\" or \"no\"." << std::endl;
		if( question_yesno() )
		{
			record_dialog << "yes" << std::endl;
			Storage::real scale[3] = { 1.0, 1.0, 1.0 }, shift[3] = { 0, 0, 0 };
			Storage::real minmax[6] = {1.0e20, -1.0e20, 1.0e20, -1.0e20, 1.0e20, -1.0e20};
			for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
			{
				Storage::real_array c = it->Coords();
				if( c[0] < minmax[0] ) minmax[0] = c[0];
				if( c[0] > minmax[1] ) minmax[1] = c[0];
				if( c[1] < minmax[2] ) minmax[2] = c[1];
				if( c[1] > minmax[3] ) minmax[3] = c[1];
				if( c[2] < minmax[4] ) minmax[4] = c[2];
				if( c[2] > minmax[5] ) minmax[5] = c[2];
			}

			std::cout << "Your mesh is currently contained at the box: ";
			std::cout << minmax[0] << " < x < " << minmax[1] << ", " << minmax[2] << " < y < " << minmax[3] << ", " << minmax[4] << " < z < " <<  minmax[5] << std::endl;
			std::cout << "Please enter three numbers divided by space," << std::endl;
			std::cout << "shift by x, y and z." << std::endl;
			std::cin >> shift[0] >> shift[1] >> shift[2];

			record_dialog << shift[0] << " " << shift[1] << " " << shift[2] << std::endl;

			std::cout << "Please enter three numbers divided by space," << std::endl;
			std::cout << "rescale in x, y, and z directions." << std::endl;
			std::cin >> scale[0] >> scale[1] >> scale[2];

			record_dialog << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;

			std::cout << "Please wait while I rescale the mesh..." << std::endl;
			if (m->GetDimensions() == 3)
			{
				for (Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); it++)
				{
					Storage::real_array cnt = it->Coords();
					cnt[0] = cnt[0] * scale[0] + shift[0];
					cnt[1] = cnt[1] * scale[1] + shift[1];
					cnt[2] = cnt[2] * scale[2] + shift[2];
				}
			}
			else
			{
				for (Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); it++)
				{
					Storage::real_array cnt = it->Coords();
					cnt[0] = cnt[0] * scale[0] + shift[0];
					cnt[1] = cnt[1] * scale[1] + shift[1];
				}
			}
			std::cout << "Done." << std::endl;
		}
		else record_dialog << "no" << std::endl;
	}

	record_dialog << "/" << std::endl;

	{ //ask about partitioning?
#if defined(USE_PARTITIONER)
		if( m->GetProcessorsNumber() > 1 )
		{
			Partitioner * p = new Partitioner(m);
			p->SetMethod(Partitioner::Inner_RCM,Partitioner::Partition);
			p->Evaluate();
			delete p;
			m->Redistribute();
			m->ReorderEmpty(CELL|FACE|EDGE|NODE);
			m->ExchangeGhost(2,EDGE);
		}
#endif
	}

	std::string tensor_name = "";

	

	{
		record_dialog << "PERM" << std::endl;
restart_setup:
		std::cout << "Mesh has following data defined on it: " << std::endl;
		PrintTags(m);

		permeability_setup_help_message();

		const std::string keywords[] = {"data","zones","input","random2d","random3d"};
		switch(question_keyword(keywords,5))
		{
		case 0: //data
			{
restart_data:
				std::cout << "Following data tags are defined on grid: " << std::endl;
				PrintRealTags(m);
				std::cout << "Which of the data fields describe permeability tensor? Please enter the name." << std::endl;

				bool have_tensor = false;
				while(!have_tensor)
				{
					std::cin >> tensor_name;
					if( !m->HaveTag(tensor_name) )
					{
						std::cout << "There is no such data field: " << tensor_name << " Please try again. " << std::endl;
					}
					else 
					{
						std::cout << "Accepted for permiability tensor: " << std::endl;
						PrintTag(m->GetTag(tensor_name));

						bool tensor_ok = true;
						int Ktype = -1;
						std::cout << "Please wait while I check for consistancy of provided permiability tensor." << std::endl;

						Tag K = m->GetTag(tensor_name);

						if( K.GetSize() == ENUMUNDEF )
						{
							int bad_entries = 0, total = 0;
							for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
							{
								size_t s = it->GetDataSize(K) ;
								if( !(s == 1 || s == 3 || s == 6 || s == 9 ) ) bad_entries++;
								total++;
							}
							if( K.isDefined(FACE) )
							{
								for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->HaveData(K) )
								{
									size_t s = it->GetDataSize(K) ;
									if( !(s == 1 || s == 3 || s == 6 || s == 9 ) ) bad_entries++;
									total++;
								}
							}
							if( bad_entries ) 
							{
								std::cout << "There are " << bad_entries << " out of " << total << " entries that have unaccaptable record size." << std::endl;
								tensor_ok = false;
							}
						}
						else if( !(K.GetSize() == 1 || K.GetSize() == 3 || K.GetSize() == 6 || K.GetSize() == 9) )
						{
							std::cout << "Data field " << K.GetTagName() << " has unaccpatable record size: " << K.GetSize() << ". Acceptable sizes are: 1,3,6,9." << std::endl;
							tensor_ok = false;
						}
						else std::cout << "Your permeability data is ok." << std::endl;

						if( tensor_ok )
						{
							std::cout << "Would you like to check eigenvalues for your permeability tensor to be sure that it is positive definitive? Please enter \"yes\" or \"no\"" << std::endl;
							int bad_entries = 0, total = 0;
							if(question_yesno() )
							{
								for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
								{
									Storage::real mineigen;
									if( min_eigenvalue(4,it->RealArray(K),mineigen) || mineigen < 0.0) bad_entries++;
									total++;
								}
								if( K.isDefined(FACE) )
								{
									for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->HaveData(K) )
									{
										Storage::real mineigen;
										if( min_eigenvalue(4,it->RealArray(K),mineigen) || mineigen < 0.0) bad_entries++;
										total++;
									}
								}
								if( bad_entries ) 
								{
									std::cout << "There are " << bad_entries << " out of " << total << " entries that have negative or complex eigenvalues." << std::endl;
									tensor_ok = false;
								}
								else std::cout << "Your permeability is ok." << std::endl;
							}
							have_tensor = true;
						}
						if( !tensor_ok )
						{
							std::cout << "How would you like to proceed?" << std::endl;
							std::cout << "restart  - choose new permeability input method;" << std::endl;
							std::cout << "data     - select another data for permeability";
							if( have_tensor ) std::cout << ";" << std::endl << "continue - continue with current choice." << std::endl;
							else std::cout << "." << std::endl;
							std::cout << "Please enter one of the keywords." << std::endl;
							std::string keywords[] = {"restart","data","continue"};
							switch( question_keyword(keywords,2+(have_tensor?1:0)) )
							{
							case 0: goto restart_setup; break;
							case 1: goto restart_data; break;
							case 2: break;
							}
						}
					}
				}
				record_dialog << "data " << tensor_name << std::endl;
			}
			break;
		case 1: //zones
			{
				std::cout << "There are several options to describe permeability by number of entries." << std::endl;
				std::cout << "1 - your tensor is represented by single scalar; " << std::endl;
				std::cout << "3 - your tensor has distinctive diagonal entries; " << std::endl;
				std::cout << "6 - you have full symmetric tensor; " << std::endl;
				std::cout << "9 - you have full tensor. " << std::endl;
				std::cout << "Please enter the number of entries that describe your permeability." << std::endl;
				Tag Ktensor;
				const std::string keywords2[] = {"1","3","6","9"};
				int Ktype = question_keyword(keywords2,4);


restart_zonal:
				std::cout << "Following numeric data tags are defined on grid: " << std::endl;
				PrintNumericTags(m,CELL|FACE);
				std::cout << "Which of the data fields describe zones?" << std::endl;
				std::cout << "Please enter the name." << std::endl;
				bool have_tensor = false;

				tensor_name = "ZONAL_TENSOR";
				if( m->HaveTag(tensor_name) )
				{
					int attach = 1;
					std::string newname;
					do
					{
						std::stringstream str;
						str << tensor_name << attach++;
						newname = str.str();
					}
					while(m->HaveTag(newname));
					tensor_name = newname;
				}

				switch(Ktype)
				{
				case 0: Ktensor = m->CreateTag(tensor_name,DATA_REAL,CELL|FACE,FACE,1); break;
				case 1: Ktensor = m->CreateTag(tensor_name,DATA_REAL,CELL|FACE,FACE,3); break;
				case 2: Ktensor = m->CreateTag(tensor_name,DATA_REAL,CELL|FACE,FACE,6); break;
				case 3: Ktensor = m->CreateTag(tensor_name,DATA_REAL,CELL|FACE,FACE,9); break;
				}
				

				while(!have_tensor)
				{
					std::string zonal_name;
					std::cin >> zonal_name;
					if( !m->HaveTag(zonal_name) )
					{
						std::cout << "There is no such data field: " << tensor_name << std::endl;
						std::cout << "Please try again." << std::endl;
					}
					else 
					{
						int entry = 0, minentry = 10000000, maxentry = 0;
						Tag zone = m->GetTag(zonal_name);
						if( zone.GetSize() == ENUMUNDEF )
						{
							int types[2] = {0,0};
							for(Mesh::iteratorElement it = m->BeginElement(CELL | (zone.isDefined(FACE)? FACE : NONE)); it != m->EndElement(); ++it) if( it->HaveData(zone) )
							{
								size_t num = it->GetDataSize(zone);
								if( minentry > num ) minentry = num;
								if( maxentry < num ) maxentry = num;
								types[it->GetElementNum()-2]++;
							}
							std::cout << "Data field covers " << types[1] << " out of " << m->NumberOfCells() << " cells";
							if( types[0] > 0 ) std::cout << " and " << types[0] << " out of " << m->NumberOfFaces() << " faces";
							std::cout << std::endl;
							if( minentry != maxentry ) 
							{
								std::cout << "Data " << zone.GetTagName() << " have from " << minentry << " to " << maxentry << " entries per element." << std::endl;
								std::cout << "Element will be ignored if entry is absent on it. " << std::endl;
							}
							else if( maxentry > 1 )
							{
								std::cout << "Data " <<  zone.GetTagName() << " have " << maxentry << " entries per element." << std::endl;
							}
							
						}
						else minentry = maxentry = zone.GetSize();



						if( maxentry > 1 )
						{
restart_entry:
							std::cout << "Which data entry describes zone?" << std::endl;
							entry = 0;	
							while( entry == 0 )
							{
								std::cout << "Please enter the number starting from 1. " << std::endl;
								std::cin >> entry;
							}
							entry--;
						}

						if( zone.GetDataType() == DATA_INTEGER || zone.GetDataType() == DATA_BULK)
						{
							std::map<int,intpair> values;
							int types[2] = {0,0};
							std::cout << "Wait while I analyse " << zone.GetTagName() << "." << std::endl;
							for(Mesh::iteratorElement it = m->BeginElement(CELL | (zone.isDefined(FACE)? FACE : NONE)); it != m->EndElement(); ++it) if( it->HaveData(zone) )  
							{
								if( zone.GetDataType() == DATA_INTEGER )
								{
									Storage::integer_array arr = it->IntegerArray(zone);
									if( arr.size() > entry ) 
									{
										values[arr[entry]][it->GetElementNum()-2]++;
										types[it->GetElementNum()-2]++;
									}
								}
								else
								{
									Storage::bulk_array arr = it->BulkArray(zone);
									if( arr.size() > entry ) 
									{
										values[arr[entry]][it->GetElementNum()-2]++;
										types[it->GetElementNum()-2]++;
									}
								}
							}
							std::cout << "Done." << std::endl;
restart_message:
							std::cout << "Entry covers " << types[1] << " out of " << m->NumberOfCells() << " cells";
							if( types[0] > 0 ) std::cout << " and " << types[0] << " out of " << m->NumberOfFaces() << " faces";
							std::cout << std::endl;
							bool intervals = false;
							if( values.empty() )
							{
								std::cout << "There are no zones in your data." << std::endl;
								std::cout << "restart - choose new permeability input method;" << std::endl;
								std::cout << "zone    - select another zone for permeability;" << std::endl;
								std::cout << "entry   - select another zonal entry." << std::endl;
								std::cout << "Please enter one of the keywords." << std::endl;
								std::string keywords[] = {"restart","zone","entry"};
								switch( question_keyword(keywords,3) )
								{
								case 0: m->DeleteTag(Ktensor); goto restart_setup;
								case 1: m->DeleteTag(Ktensor); goto restart_zonal;
								case 2: goto restart_entry;
								}
							}
							else
							{

								if( values.size() == 1 )
									std::cout << "There is only one zone." << std::endl;
								else std::cout << "There are " << values.size() << " distinctive values in your data." << std::endl;
								std::cout << "How would you like to proceed?" << std::endl;
								std::cout << "restart  - choose new permeability input method;" << std::endl;
								std::cout << "zone     - select another zone for permeability;" << std::endl;
								std::cout << "entry    - select another zonal entry;" << std::endl;
								std::cout << "list     - list values and elements for zones;" << std::endl;
								std::cout << "continue - continue with current choice." << std::endl;
								std::cout << "Please enter one of the keywords." << std::endl;
								std::string keywords[] = {"restart","zone","entry","list","continue"};
								switch( question_keyword(keywords,5) )
								{
								case 0: m->DeleteTag(Ktensor); goto restart_setup;
								case 1: m->DeleteTag(Ktensor); goto restart_zonal;
								case 2: goto restart_entry;
								case 3:
									{
										int q = 0;
										int wn = std::max<int>(1,ceil(log((double)values.size())/log(10.0)))+1;
										int wv = std::max<int>(5,ceil(log((double)values.rbegin()->first)/log(10.0)))+1;
										std::cout << std::setw(wn) << " " << std::setw(wv) << "value" << std::setw(20) << "cells" << std::setw(20) << "faces" << std::endl;
										for(std::map<int,intpair>::iterator jt = values.begin(); jt != values.end(); ++jt)
										{
											std::cout << std::setw(wn) << q++ << std::setw(wv) << jt->first << std::setw(20) << jt->second[1] << std::setw(20) << jt->second[0] << std::endl;
										}
										goto restart_message;
									}
								}


								if( values.size() > 5 )
								{
									std::cout << "Would you like to indicate zone by intervals?" << std::endl;
									std::cout << "Otherwise you will have to enter data for zones for each value." << std::endl;
									std::cout << "Please enter \"yes\" or \"no\"." << std::endl;
									intervals = question_yesno();
								}

								record_dialog << "zones " << Ktype << " " << tensor_name << " " << zone.GetTagName() << " " << entry << std::endl;
							}
							
							
							if( intervals )
							{
								std::cout << "Values of your data have entries from " << values.begin()->first << " to " << values.rbegin()->first << std::endl;
								bool outputmsg = true;
								while(true)
								{
									if( outputmsg )
									{
										std::cout << "Please enter two numbers that indicate interval." << std::endl;
										std::cout << "If you are done write \"stop\"" << std::endl;
										outputmsg = false;
									}
									char intervalstr[2048];
									std::cin.getline(intervalstr,2048);
									if( strlen(intervalstr) == 0 ) continue;
									outputmsg = true;
									if( !strcmp(intervalstr,"stop") ) break;
									else
									{
										std::string interval(intervalstr);
										std::stringstream readit(interval);
										int ibeg, iend;
										readit >> ibeg >> iend;
										std::cout << "Please enter " << Ktensor.GetSize() << " tensor values for interval [" << ibeg << "," << iend << "]:" << std::endl;
										Storage::real K[9];
										record_dialog << ibeg << " " << iend;
										for(int q = 0; q < Ktensor.GetSize(); ++q) 
										{
											std::cin >> K[q];
											record_dialog << " " << K[q];
										}
										record_dialog << std::endl;
										int types[2] = {0,0};
										for(Mesh::iteratorElement it = m->BeginElement(CELL | (zone.isDefined(FACE)? FACE : NONE)); it != m->EndElement(); ++it) if( it->HaveData(zone) )  
										{
											if( zone.GetDataType() == DATA_INTEGER )
											{
												Storage::integer_array arr = it->IntegerArray(zone);
												if( arr.size() > entry && arr[entry] >= ibeg && arr[entry] <= iend) 
												{
													memcpy(&it->RealArray(Ktensor)[0],K,sizeof(Storage::real)*Ktensor.GetSize());
													types[it->GetElementNum()-2]++;
												}
											}
											else
											{
												Storage::bulk_array arr = it->BulkArray(zone);
												if( arr.size() > entry && arr[entry] >= ibeg && arr[entry] <= iend) 
												{
													memcpy(&it->RealArray(Ktensor)[0],K,sizeof(Storage::real)*Ktensor.GetSize());
													types[it->GetElementNum()-2]++;
												}
											}
										}
										std::cout << "Filled " << types[1] << " cells";
										if( types[0] ) std::cout << " and " << types[0] << " faces";
										std::cout << " on interval [" << ibeg << "," << iend << "]." << std::endl;

									}
								}
							}
							else
							{
								int nzone = 0;
								for(std::map<int,intpair>::iterator jt = values.begin(); jt != values.end(); ++jt)
								{
									std::cout << "Please enter " << Ktensor.GetSize() << " tensor values for zone " << jt->first << " cells " << jt->second[1] << " faces " << jt->second[0] << ":" << std::endl;
									Storage::real K[9];
									record_dialog << jt->first << " " << jt->first;
									for(int q = 0; q < Ktensor.GetSize(); ++q) 
									{
										std::cin >> K[q];
										record_dialog << " " << K[q];
									}
									record_dialog << std::endl;
									int types[2] = {0,0};
									for(Mesh::iteratorElement it = m->BeginElement(CELL | (zone.isDefined(FACE)? FACE : NONE)); it != m->EndElement(); ++it) if( it->HaveData(zone) )  
									{
										if( zone.GetDataType() == DATA_INTEGER )
										{
											Storage::integer_array arr = it->IntegerArray(zone);
											if( arr.size() > entry && arr[entry] == jt->first) 
											{
												memcpy(&it->RealArray(Ktensor)[0],K,sizeof(Storage::real)*Ktensor.GetSize());
												types[it->GetElementNum()-2]++;
											}
										}
										else
										{
											Storage::bulk_array arr = it->BulkArray(zone);
											if( arr.size() > entry && arr[entry] == jt->first) 
											{
												memcpy(&it->RealArray(Ktensor)[0],K,sizeof(Storage::real)*Ktensor.GetSize());
												types[it->GetElementNum()-2]++;
											}
										}
									}
									std::cout << "Filled " << types[1] << " cells";
									if( types[0] ) std::cout << " and " << types[0] << " faces";
									std::cout << " of zone " << jt->first << "." << std::endl;
									nzone++;
								}
							}

							record_dialog << "stop" << std::endl;
						}
						else if(zone.GetDataType() == DATA_REAL )
						{
							Storage::real vmin = 1.0e20, vmax = -1.0e20;
							std::cout << "Wait while I analyse " << zone.GetTagName() << "." << std::endl;
							int types[2] = {0,0};
							for(Mesh::iteratorElement it = m->BeginElement(CELL | (zone.isDefined(FACE)? FACE : NONE)); it != m->EndElement(); ++it) if( it->HaveData(zone) )  
							{
								Storage::real_array arr = it->RealArray(zone);
								if( arr.size() > entry )
								{
									if( arr[entry] < vmin ) vmin = arr[entry];
									if( arr[entry] > vmax ) vmax = arr[entry];
									types[it->GetElementNum()-2]++;
								}
							}
							std::cout << "Done." << std::endl;
							std::cout << "Entry covers " << types[1] << " out of " << m->NumberOfCells() << " cells";
							if( types[0] > 0 ) std::cout << " and " << types[0] << " out of " << m->NumberOfFaces() << " faces";
							std::cout << std::endl;
							std::cout << "Your data has values between " << vmin << " and " << vmax << std::endl;


							std::cout << "restart  - choose new permeability input method;" << std::endl;
							std::cout << "zone     - select another zone for permeability;" << std::endl;
							std::cout << "entry    - select another zonal entry;" << std::endl;
							std::cout << "continue - continue with current choice." << std::endl;
							std::cout << "Please enter one of the keywords." << std::endl;
							std::string keywords[] = {"restart","zone","entry","continue"};
							switch( question_keyword(keywords,4) )
							{
							case 0: m->DeleteTag(Ktensor); goto restart_setup;
							case 1: m->DeleteTag(Ktensor); goto restart_zonal;
							case 2: goto restart_entry;
							}

							record_dialog << "zones " << Ktype << " " << tensor_name << " " << zone.GetTagName() << " " << entry << std::endl;

							while(true)
							{
								std::cout << "Please enter two numbers that indicate interval." << std::endl;
								std::cout << "If you are done write \"stop\"" << std::endl;
								std::string interval;
								std::cin >> interval;
								if( interval == "stop" ) break;
								else
								{
									std::stringstream readit(interval);
									Storage::real ibeg, iend;
									readit >> ibeg >> iend;
									std::cout << "Please enter " << Ktensor.GetSize() << " tensor values for interval [" << ibeg << "," << iend << "]:" << std::endl;
									Storage::real K[9];
									record_dialog << ibeg << " " << iend;
									for(int q = 0; q < Ktensor.GetSize(); ++q) 
									{
										std::cin >> K[q];
										record_dialog << " " << K[q];
									}
									record_dialog << std::endl;
									int types[2] = {0,0};
									for(Mesh::iteratorElement it = m->BeginElement(CELL | (zone.isDefined(FACE)? FACE : NONE)); it != m->EndElement(); ++it) if( it->HaveData(zone) )  
									{
										if( zone.GetDataType() == DATA_INTEGER )
										{
											Storage::integer_array arr = it->IntegerArray(zone);
											if( arr.size() > entry && arr[entry] >= ibeg && arr[entry] <= iend) 
											{
												memcpy(&it->RealArray(Ktensor)[0],K,sizeof(Storage::real)*Ktensor.GetSize());
												types[it->GetElementNum()-2]++;
											}
										}
										else
										{
											Storage::bulk_array arr = it->BulkArray(zone);
											if( arr.size() > entry && arr[entry] >= ibeg && arr[entry] <= iend) 
											{
												memcpy(&it->RealArray(Ktensor)[0],K,sizeof(Storage::real)*Ktensor.GetSize());
												types[it->GetElementNum()-2]++;
											}
										}
									}
									std::cout << "Filled " << types[1] << " cells";
									if( types[0] ) std::cout << " and " << types[0] << " faces";
									std::cout << " on interval [" << ibeg << "," << iend << "]." << std::endl;
								}
							}

							record_dialog << "stop" << std::endl;
						}
						have_tensor = true;
					}
				}
			}
			break;
		case 2: //input
			{
				std::cout << "This brench is not finished yet" << std::endl;
			}
			break;
		case 3: //random2d
			{
				std::cout << "Please enter the number for the random seed, or \"time\" to get seed from current time." << std::endl;
				std::string seed;
				std::cin >> seed;
				if( seed == "time")
					tensor_seed = time(NULL);
				else tensor_seed = atoi(seed.c_str());
				std::cout << "Got seed " << tensor_seed << std::endl;
				std::cout << "The tensor for each cell will be represented by matrix multiplication R(phi) * diag(Kx,Ky) * R(-phi).";
				std::cout << "Where random Kx, Ky, phi will satisfy:" << std::endl; 
				std::cout << "Kx_min < Kx < Kx_max," << std::endl;
				std::cout << "Ky_min < Ky < Ky_max," << std::endl;
				std::cout << "phi_min < phi < phi_max." << std::endl;
				std::cout << "Please enter Kx_min and Kx_max." << std::endl;
				std::cin >> Kxmin >> Kxmax;
				std::cout << "Please enter Ky_min and Ky_max." << std::endl;
				std::cin >> Kymin >> Kymax;
				std::cout << "Please enter phi_min and phi_max." << std::endl;
				std::cin >> Phimin >> Phimax;
				std::cout << "Please wait while I put permeability field to cells of your mesh." << std::endl;
				const Storage::real pi = 3.1415926535;
				tensor_name = "RANDOM2D_TENSOR";
				Tag K = m->CreateTag(tensor_name,DATA_REAL,CELL,NONE,6);
				srand(tensor_seed);
				for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				{
					Storage::real_array pK = it->RealArrayDF(K);
					Storage::real Kx = rand() / (1.0*RAND_MAX) * 990.0 + 10.0;
					Storage::real Ky = rand() / (1.0*RAND_MAX) * 90.0 + 10.0;
					Storage::real Phi = (rand() / (1.0*RAND_MAX) * (pi-0.01)+0.01)/2.0;
					Storage::real cPhi = cos(Phi), sPhi = sin(Phi), cPhi2 = cPhi*cPhi, sPhi2 = sPhi*sPhi, csPhi = cPhi*sPhi;
					pK[0] = Kx*cPhi2 + Ky*sPhi2;
					pK[1] = (Kx-Ky)*csPhi;
					pK[2] = 0.0;
					pK[3] = Kx*sPhi2 + Ky*cPhi2;
					pK[4] = 0.0;
					pK[5] = 1.0;
				}
				std::cout << "Done." << std::endl;
				record_dialog << "random2d " << tensor_seed << " " << Kxmin << " " << Kxmax << " " << Kymin << " " << Kymax << " " << Phimin << " " << Phimax << std::endl;
			}
			break;
		case 4: //random3d
			{
				tensor_name = "random3d";
			}
			break;
		}
		record_dialog << "/" << std::endl;
	}

	record_dialog << "SCHEME" << std::endl;

	std::string str_scheme;
	{

		scheme_setup_help_message();

		bool have_scheme = false;
		
		int scheme = -1;
		while( !have_scheme)
		{
			std::string keywords[] = 
			{
				"TPFA_OLD","TPFA",
				"NTPFA-A","NMPFA-A","MPFA-A",
				"NTPFA-B","NMPFA-B","MPFA-B",
				"NTPFA-C",
				"MPFA-O","MPFA-L","MPFA-G",
				"LDMP"
			};
			scheme = question_keyword(keywords,13);
			if( scheme >= 0 ) 
			{
				have_scheme = true;
				if( scheme > 1 && scheme < 9 )
				{
					std::cout << "You use the scheme that have additional options." << std::endl;
					std::cout << "Harmonic averaging points allow for compact interpolation and" << std::endl;
					std::cout << "greatly increase efficiency, but may affect accuracy." << std::endl;
					std::cout << "Would you like to allow for scheme to use harmonic averaging points?" << std::endl;
					std::cout << "Please enter \"yes\" or \"no\"." << std::endl;
					if( question_yesno() ) allow_harmonic |= AVG_HARMONIC;
					std::cout << "Usign information about presence of neumann or dirichlet boundary" << std::endl;
					std::cout << "conditions will improve efficiency, but may affect accuracy." << std::endl;
					std::cout << "If other boundary conditions are present they will be treated correctly." << std::endl;
					std::cout << "Would you like to allow for scheme to use neumann or dirichlet" <<std::endl;
					std::cout << "boundary conditions information?" << std::endl;
					std::cout << "Please enter \"yes\" or \"no\"." << std::endl;
					if( question_yesno() ) allow_harmonic |= AVG_NEUMANN;
				}
				switch (scheme)
				{
				case 0: d = new tpfa(aut, m, tensor_name); break;
				case 1: d = new tpfa2(aut, m, tensor_name); break;
				case 2: d = new ntpfa_a(aut, m, tensor_name); break;
				case 3: d = new nmpfa_a(aut, m, tensor_name); break;
				case 4: d = new mpfa_a(aut, m, tensor_name); break;
				case 5: d = new ntpfa_b(aut, m, tensor_name); break;
				case 6: d = new nmpfa_b(aut, m, tensor_name); break;
				case 7: d = new mpfa_b(aut, m, tensor_name); break;
				case 8: d = new ntpfa_c(aut, m, tensor_name); break;
				case 9: d = new mpfa_o(aut, m, tensor_name); break;
				case 10: d = new mpfa_l(aut, m, tensor_name); break;
				case 11: d = new mpfa_g(aut, m, tensor_name); break;
				case 12: d = new ldmp(aut,m,tensor_name); break;
				default: d = NULL;  std::cout << "Sorry, scheme not implemented." << std::endl;
				}
	

				if (d == NULL)
				{
					std::cout << "Discretization was not generated due to some error." << std::endl;
					std::cout << "Would you like to try another? " << std::endl;
					if( question_yesno() ) have_scheme = false;
					else exit(0);
				}
				else
				{
					record_dialog << keywords[scheme];
					if( scheme > 1 && scheme < 9 )
					{
						if( allow_harmonic & AVG_HARMONIC ) record_dialog << " yes" << std::endl;
						else record_dialog << " no" << std::endl;
						if( allow_harmonic & AVG_NEUMANN ) record_dialog << " yes" << std::endl;
						else record_dialog << " no" << std::endl;
					}
					std::cout << std::endl;
				}
			}
		}
	}
	record_dialog << "/" << std::endl;

	record_dialog << "FRACTURES" << std::endl;
	{
		std::cout << "Would you like to indicate presence of fractuers?" << std::endl;
		std::cout << "Please enter \"yes\" or \"no\"." << std::endl;
		if( question_yesno() )
		{
			std::cout << "Mesh has following data defined on it: " << std::endl;
			PrintTags(m);
			std::cout << "You may indicate the presence of fractures in two different ways:" << std::endl;
			std::cout << "data   - you have a data tag which represents volume factors;" << std::endl;
			std::cout << "zones  - you have a data tag whose presence or value indicate fractures;" << std::endl;
			std::cout << "input  - you will input identificators and volume factors of faces." << std::endl;
			std::cout << "How would you like to distinguish fractures?" << std::endl;
			std::cout << "Please enter one of the keywords above." << std::endl;
			const std::string keywords[] = {"data","zones","input"};
			MIDType frac_marker = m->CreateMarker();
			Tag volume_factor;
			switch(question_keyword(keywords,2))
			{
			case 0://data
				{
					std::cout << "Following data is defined on faces:" << std::endl;
					PrintRealTags(m,FACE);


				}
				break;
			case 1: //zones
				{
				}
				break;
			case 2: //input
				{
					m->ReorderEmpty(FACE);
					bool outputmsg = true;
					record_dialog << "input" << std::endl;
					while(true)
					{
						volume_factor = m->CreateTag("VOLUME_FACTOR",DATA_REAL,FACE,FACE,1);
						if( outputmsg )
						{
							std::cout << "Please enter faces identificator from 1 up to " <<  m->MaxLocalIDFACE() << " or \"stop\" to finish." << std::endl;
							outputmsg = false;
						}
						char idstr[2048];
						std::cin.getline(idstr,2048);
						if( strlen(idstr) == 0 ) continue;
						outputmsg = true;
						if( !strcmp(idstr,"stop") ) break;
						else
						{
							int id = atoi(idstr)-1;
							if( id < 0 || id > m->MaxLocalIDFACE() ) std::cout << "Identificator " << id << " is out of bounds." << std::endl;
							else
							{
								Face * f = m->FaceByLocalID(id);
								f->SetMarker(frac_marker);
								double volfactor = (((f->BackCell() != NULL)? f->BackCell()->Volume() : 0.0)+((f->FrontCell() != NULL)? f->FrontCell()->Volume() : 0.0))/f->Area();
								std::cout << "Maximum reasonable volume factor for this face is " << volfactor << "." << std::endl;
								std::cout << "It is calculated from volume of neighbouring cells divided by area of the face." << std::endl;
								std::cout << "You are free to enter any volume factor you'd like." << std::endl;
								std::cout << "Please enter volume factor for current face." << std::endl;
								std::cin >> volfactor;
								f->Real(volume_factor) = volfactor;
								record_dialog << id << " " << volfactor << std::endl;
							}
						}
					}
					record_dialog << "stop" << std::endl;
				}
				break;
			}
			

			

			std::cout << "Availible options for fracture treatment: " << std::endl;
			std::cout << "NONE    - fractures are treated as no-flow boundary condition;" << std::endl;
			std::cout << "TPFA-MF - matrix-fracture and fracture-fracture connections are treated by two-point flux approximation;" << std::endl;
			std::cout << "TPFA-FF - fracture-fracture connections are treated by two-point flux approximation;" << std::endl;
			std::cout << "FULL    - fractures are fully incorporated into the scheme." << std::endl;
			int behavior = -1;
			while( behavior == -1 )
			{
				std::string entry;
				std::cin >> entry;
				if( entry == "NONE" ) behavior = FRAC_NONE;
				else if( entry == "TPFA-MF" ) behavior = FRAC_TPFA_FM;
				else if( entry == "TPFA-FF" ) behavior = FRAC_TPFA_FF;
				else if( entry == "FULL" ) behavior = FRAC_FULL;
				else std::cout << "I don't understand " << entry << ", please try again. " << std::endl;
			}
		}
		else std::cout << " OK, no fractures " << std::endl;
	}
	record_dialog << "/" << std::endl;

	record_dialog << "WELLS" << std::endl;
	{
	}
	record_dialog << "/" << std::endl;
	
	if (!m->HaveTag(tensor_name))
	{
		if( tensor_name == "random2d" )
		{
			const Storage::real pi = 3.1415926535;
			tensor_name = "K";
			Tag K = m->CreateTag("K",DATA_REAL,CELL,NONE,6);
			srand(12345);
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			{
				Storage::real_array pK = it->RealArrayDF(K);
				Storage::real Kx = rand() / (1.0*RAND_MAX) * 990.0 + 10.0;
				Storage::real Ky = rand() / (1.0*RAND_MAX) * 90.0 + 10.0;
				Storage::real Phi = (rand() / (1.0*RAND_MAX) * (pi-0.01)+0.01)/2.0;
				Storage::real cPhi = cos(Phi), sPhi = sin(Phi), cPhi2 = cPhi*cPhi, sPhi2 = sPhi*sPhi, csPhi = cPhi*sPhi;
				pK[0] = Kx*cPhi2 + Ky*sPhi2;
				pK[1] = (Kx-Ky)*csPhi;
				pK[2] = 0.0;
				pK[3] = Kx*sPhi2 + Ky*cPhi2;
				pK[4] = 0.0;
				pK[5] = 1.0;
			}
		}
		else
		{
			std::cout << "Permiability tensor " << tensor_name << " is not defined on the grid" << std::endl;
			return -1;
		}
	}

	


	if (argc > 7) trans_scale = atof(argv[7]);
	if (argc > 8) vol_scale = atof(argv[8]);
	if (argc > 9) wi_scale = atof(argv[9]);

	


	Mesh::GeomParam t;
	t[MEASURE] = CELL | FACE;
	t[ORIENTATION] = FACE;
	t[NORMAL] = FACE;
	//t[BARYCENTER] = CELL;
	t[CENTROID] = CELL | FACE;
	m->AssignGlobalID(CELL);
	m->PrepareGeometricData(t);

	
	
	


	if (d == NULL)
	{
		std::cout << "Discretization not generated" << std::endl;
	}
	else
	{
		discr_basic * bnd_d = d;

		d->AllowTypes(allow_harmonic);
		d->Init();
		

		//std::fstream fout("out.txt", std::ios::out);
		//d->Export(fout,trans_scale,vol_scale);
		//fout.close();

		if (argc > 4) wellsfile = std::string(argv[4]);

		wells * w = NULL;

		Tag WI, bhp;

		INMOST_DATA_REAL_TYPE init_pressure = 0, init_num = 0;

		if (wellsfile != "")
		{
			WI = m->CreateTag("WELL_INDEX", DATA_REAL, CELL, CELL, 1);
			w = new wells(m, m->GetTag(tensor_name), WI);
			std::fstream wf(wellsfile.c_str(), std::ios::in);
			w->Import(wf);
			wf.close();

			bhp = m->CreateTag("bhp",DATA_REAL,ESET,ESET,1);

			for(INMOST_DATA_ENUM_TYPE k = 0; k < w->GetWellNumber(); k++)
			{
				if( w->GetWellType(k) == well_data::Flux )
				{
					w->GetWellElements(k)->SetMarker(d->UnknownMarker());
					//w->GetWellElements(k)->Real(flux) = w->GetWell(k).GetWellControlValue();
				}
				else 
				{
					init_pressure += w->GetWell(k).GetWellControlValue();
					init_num++;
				}
			}
		}

		if( w != NULL )
		{
			for(INMOST_DATA_ENUM_TYPE k = 0; k < w->GetWellNumber(); k++)
			{
				if(  w->GetWellElements(k)->size() == 0 )
				{
					std::cout << m->GetProcessorRank() << " well " << k << " does not penetrate any cell " << std::endl;
				}
				else
				{
					std::cout << m->GetProcessorRank() << " well " << k << " penetrate " <<  w->GetWellElements(k)->size() <<  " cells " << std::endl;
				}
			}

			if( w->GetWellNumber() == 0 ) std::cout << m->GetProcessorRank() << " no wells imported" << std::endl;
		}

		if( init_num > 0 ) init_pressure /= init_num;
		
		
		INMOST_DATA_ENUM_TYPE types[5] = {0,0,0,0,0}, total[5];
		for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE|EDGE|NODE|ESET); it != m->EndElement(); ++it)
			if( it->GetStatus() != Element::Ghost && it->GetMarker(d->UnknownMarker()) )
			{
				types[it->GetElementNum()]++;
			}
		/*
		for(Mesh::iteratorElement it = m->BeginElement(FACE); it != m->EndElement(); ++it)
			if( it->GetMarker(d->UnknownMarker()) )
			{
				std::cout << "FACE " << it->GlobalID() << " on " << m->GetProcessorRank() << " status " << Element::StatusName(it->GetStatus()) << " owner " << it->Integer(m->OwnerTag()) << std::endl;
			}
		*/
		ElementType active_types = NONE;
		
		for(ElementType etype = NODE; etype <= ESET; etype = etype << 1)
		{
			if( types[ElementNum(etype)] ) 
			{
				active_types |= etype;
				//std::cout << ElementTypeName(etype) << " " << types[ElementNum(etype)] << std::endl;
			}
		}
		
#if defined(USE_MPI)
		MPI_Reduce(types,total,5,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,0,m->GetCommunicator());
#else
		memcpy(total,types,sizeof(INMOST_MPI_DATA_ENUM_TYPE)*5);
#endif

		if( m->GetProcessorRank() == 0 ) 
		{
			for(ElementType etype = NODE; etype <= ESET; etype = etype << 1)
			{
				if( total[ElementNum(etype)] > 0 )
				{
					std::cout << ElementTypeName(etype) << " " << total[ElementNum(etype)] << std::endl;
				}
			}
		}
		
		bool had_pressure = m->HaveTag("PRESSURE");
		active_types = m->SynchronizeElementType(active_types);

		ElementSet * additional = m->CreateSet(), * interpolation = m->CreateSet();

		//std::cout << "active types: " << std::endl;
		//for(ElementType etype = NODE; etype <= ESET; etype = etype << 1 )
		//	if( active_types & etype )
		//		std::cout << ElementTypeName(etype) << std::endl;
		Tag p = m->CreateTag("PRESSURE",DATA_REAL,active_types, active_types & ~CELL, 1);
		Tag pupdate, lambda_coef;//,update_coef;
		if( d->GetDiscrType() == discr_basic::NMPFA )
		{
			pupdate = m->CreateTag("PRESS_UPDATE",DATA_REAL,active_types, active_types & ~CELL,1);
			//update_coef = m->CreateTag("UPDATE_COEF",DATA_REAL,active_types, active_types & ~CELL,1);
			lambda_coef = m->CreateTag("LAMBDA_COEF",DATA_REAL,active_types, active_types & ~CELL,1);
		}

		if( !had_pressure )
		{
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				it->RealDF(p) = init_pressure + (rand() / (double)RAND_MAX)*10;
		}

		for(Mesh::iteratorElement it = m->BeginElement(active_types & ~CELL ); it != m->EndElement(); ++it)
			if( it->GetMarker(d->UnknownMarker()) )
			{
				if( it->GetStatus() != Element::Ghost )
				{
					if( !it->GetMarker(d->BoundaryMarker()) || it->GetElementType() != FACE)
					{
						interpolation->Insert(&*it);
					}
					additional->Insert(&*it);
				}
				it->Real(p) = init_pressure;// + (rand() / (double)RAND_MAX);
			}
		if( w != NULL )
		{
			for(INMOST_DATA_ENUM_TYPE k = 0; k < w->GetWellNumber(); k++)
			{
				if( w->GetWellType(k) == well_data::BHP )
					w->GetWellElements(k)->Real(bhp) = w->GetWell(k).GetWellControlValue();
			}
		}

		m->ExchangeData(p,active_types);

		INMOST_DATA_ENUM_TYPE tagp = aut->RegisterDynamicTag(p,active_types, d->UnknownMarker());
		INMOST_DATA_ENUM_TYPE welli, welle, wellp;
		if( w != NULL )
		{
			welli = aut->RegisterStaticTag(WI);
			wellp = aut->RegisterStaticTag(bhp);
			welle = aut->RegisterStencil("well",get_well_func);
		}
		else
		{
			welli = ENUMUNDEF;
			wellp = ENUMUNDEF;
			welle = ENUMUNDEF;
		}
		//m->Save("temp.pvtk");

		aut->EnumerateDynamicTags();

		Tag index_p = aut->GetDynamicIndexTag(tagp);

		//std::cout << "proc " << m->GetProcessorRank() << std::endl;

		

		expr pres = tagval(tagp);
		expr Kgradint = trans_scale*d->Grad(pres);
		expr Kgradintl, Kgradintr;
		if( nmpfa_special && d->GetDiscrType() == discr_basic::NMPFA )
		{
			//Kgradintl = trans_scale*nmpfa_d->Gradl(pres);
			//Kgradintr = trans_scale*nmpfa_d->Gradr(pres);
			Kgradintl = trans_scale*dynamic_cast<discr_nmpfa_basic *>(d)->Gradl(pres);
			Kgradintr = trans_scale*dynamic_cast<discr_nmpfa_basic *>(d)->Gradr(pres);
		}
		expr Kgradbnd = trans_scale*d->BoundaryGrad(pres);
		expr Interp[3] = {pres - d->Interp(NODE,pres),pres - d->Interp(EDGE,pres),pres - d->Interp(FACE,pres)};
		expr bhp_well = welli != ENUMUNDEF ? wi_scale*tagval(welli)*(pres - stencil(welle,tagval(wellp))) : 0.0;
		expr flux_well = welli != ENUMUNDEF ? wi_scale*tagval(welli)*(pres - stencil(welle,pres)) : 0.0;
		INMOST_DATA_REAL_TYPE err = 0.0, update_coef = 1.0;
		INMOST_DATA_ENUM_TYPE iter = 0;
		Solver::Row r;
		Solver::Matrix A;
		Solver::Vector rhs, sol;
		A.SetInterval(aut->GetFirstIndex(),aut->GetLastIndex());
		rhs.SetInterval(aut->GetFirstIndex(),aut->GetLastIndex());

		Storage::real trust_region_lambda = 0.0;//1.0e-5;

		Storage::real tnonlinear = Timer();

		do
		{
			INMOST_DATA_REAL_TYPE ret;
			std::fill(sol.Begin(),sol.End(),0);
			std::fill(rhs.Begin(),rhs.End(),0);
			for(Solver::Matrix::iterator it = A.Begin(); it != A.End(); ++it) it->Clear();

			
			INMOST_DATA_REAL_TYPE tmatrix = Timer(), tsolver, tupdate;
			
			for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); it++) if( d->needBuild(&*it) )
			{
				expr * Kgrad = &Kgradint;
				Element * r0 = it->BackCell();
				Element * r1 = it->FrontCell();
				INMOST_DATA_ENUM_TYPE i1 = r0->IntegerDF(index_p), i2;
				if( it->GetMarker(d->BoundaryMarker()) )
				{
					//special treatment for boundary face
					if( !it->GetMarker(d->UnknownMarker()) ) continue;
					else 
					{
						r1 = &*it;
						Kgrad = &Kgradbnd;
						i2 = r1->Integer(index_p);
					}
				} else i2 = r1->IntegerDF(index_p);
				Element::Status s0 = r0->GetStatus(), s1 = r1->GetStatus();
				if( s0 == Element::Ghost && s1 == Element::Ghost ) continue;
				if( nmpfa_special && d->GetDiscrType() == discr_basic::NMPFA && Kgrad == &Kgradint )
				{
					if( s0 != Element::Ghost ) 
					{
						r.Clear();
						ret = aut->Derivative(Kgradintl,&*it,r,1.0,NULL);
						rhs[i1] += ret; add_row(A[i1],r,-1.0);
					}
					if( s1 != Element::Ghost )
					{
						r.Clear();
						ret = aut->Derivative(Kgradintr,&*it,r,1.0,NULL);
						rhs[i2] -= ret; add_row(A[i2],r,+1.0);
					}
				}
				else
				{
					r.Clear();
					ret = aut->Derivative(*Kgrad,&*it,r,1.0,NULL);
					if( s0 != Element::Ghost ) {rhs[i1] += ret; add_row(A[i1],r,-1.0);}
					if( s1 != Element::Ghost ) {rhs[i2] -= ret; add_row(A[i2],r,+1.0);}
				}
			}
			for(ElementSet::iterator it = interpolation->begin(); it != interpolation->end(); ++it)
			{
				INMOST_DATA_ENUM_TYPE i1 = it->Integer(index_p);//aut->GetDynamicIndex(&*it,tagp);
				rhs[i1] -= aut->Derivative(Interp[it->GetElementNum()],&*it,A[i1],1.0,NULL);
			}

			if( w != NULL )
			{
				for(INMOST_DATA_ENUM_TYPE i = 0; i < w->GetWellNumber(); i++)
				{
					wnum = i;
					if( w->GetWellType(i) == well_data::BHP )
					{
						for(ElementSet::iterator it = w->GetWellElements(i)->begin(); it != w->GetWellElements(i)->end(); ++it)
						{
							INMOST_DATA_ENUM_TYPE i1 = it->IntegerDF(index_p);// aut->GetDynamicIndex(&*it,tagp);
							if( it->GetStatus() != Element::Ghost ) rhs[i1] -= aut->Derivative(bhp_well,&*it,A[i1],1.0,w);
						}
					}
					else if( w->GetWellType(i) == well_data::Flux )
					{
						INMOST_DATA_ENUM_TYPE i2 = w->GetWellElements(i)->Integer(index_p); //aut->GetDynamicIndex(w->GetWellElements(i),tagp);
						for(ElementSet::iterator it = w->GetWellElements(i)->begin(); it != w->GetWellElements(i)->end(); ++it)
						{
							INMOST_DATA_ENUM_TYPE i1 = it->IntegerDF(index_p);// aut->GetDynamicIndex(&*it,tagp);
							r.Clear();
							ret = aut->Derivative(flux_well,&*it,r,1.0,w);
							if( it->GetStatus() != Element::Ghost ) 
							{
								rhs[i1] += ret;
								add_row(A[i1],r,-1.0);
							}
							rhs[i2] -= ret;
							add_row(A[i2],r,1.0);
						}
						rhs[i2] += w->GetWell(i).GetWellControlValue();
					}
				}
			}
			


			//Mimic trust region methods
			/*
			if( d->GetDiscrType() == discr_basic::NMPFA )
			{
				trust_region_lambda = 0.0;
				Storage::real sum_p, sum, dist, temp;
				Tag lrtrp_tag[2] = {m->GetTag("TRP_L"), m->GetTag("TRP_R")};
				Tag lrcoef_tag[2] = {m->GetTag("COEF_L"), m->GetTag("COEF_R")};
				for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if(it->GetMarker(d->UnknownMarker()))
				{
					Storage::real_array lrcoef[2] = {it->RealArray(lrcoef_tag[0]), it->RealArray(lrcoef_tag[1])};
					Storage::reference_array lrtrp[2] = {it->ReferenceArray(lrtrp_tag[0]), it->ReferenceArray(lrtrp_tag[1])};
					for(int q = 0; q < 2; q++)
					{
						dist = 0.0;
						sum_p = sum = 0.0;
						for(int r = 0; r < 4; r++) if( lrtrp[q][r] != NULL )
						{
							sum_p += lrcoef[q][r] * lrtrp[q][r]->Real(p);
							sum += lrcoef[q][r]*lrcoef[q][r];
						}
						dist = fabs((sum_p) / (sqrt(sum)+1.0e-14));
						if( dist > 0.0 )
						{
							dist = 1.0/dist;
							if( trust_region_lambda < dist ) trust_region_lambda = dist;
						}
					}
				}
			}
			*/

			//different trust region sizes
			/*
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) 
			{
				if( it->GetMarker(d->UnknownMarker()) && it->GetStatus() != Element::Ghost )
				{
					INMOST_DATA_ENUM_TYPE i = it->IntegerDF(index_p);//aut->GetDynamicIndex(&*it,tagp);
					INMOST_DATA_REAL_TYPE & diag = A[i][i], & lambda = it->RealDF(lambda_coef);
					diag = (1+lambda)*diag;
					lambda = 0.0;

				}
			}

			for(ElementSet::iterator it = additional->begin(); it != additional->end(); ++it)
			{
				INMOST_DATA_ENUM_TYPE i = it->Integer(index_p);//aut->GetDynamicIndex(&*it,tagp);
				INMOST_DATA_REAL_TYPE & diag = A[i][i], & lambda = it->Real(lambda_coef);
				diag = (1+lambda)*diag;
				lambda = 0.0;
			}
			*/
			//uniform trust region size
			
			if( d->GetDiscrType() == discr_basic::NMPFA && trust_region_lambda > 0.0 )
			{
				for(INMOST_DATA_ENUM_TYPE i = aut->GetFirstIndex(); i != aut->GetLastIndex(); ++i) 
				{
					INMOST_DATA_REAL_TYPE & diag = A[i][i];
					diag = (1+trust_region_lambda)*diag;
				}
			}
			

			err = 0;
			for(INMOST_DATA_ENUM_TYPE i = aut->GetFirstIndex(); i != aut->GetLastIndex(); ++i) 	{err += rhs[i]*rhs[i];}

			

			tmatrix = Timer() - tmatrix;

			//std::cout << m->GetProcessorRank() << " " << aut->GetFirstIndex() << " " << aut->GetLastIndex() << " " << err << std::endl;

			err = sqrt(m->Integrate(err));
			
			if( m->GetProcessorRank() == 0 ) std::cout << iter << " residual: " << err << "\t\t\t\t\t" << std::endl;
			if(err < 1.0e-4 || iter > 10) break;

			tsolver = Timer();
			//Solver S(Solver::INNER_ILU2);
			A.Save(str_scheme+"_A.mtx");
			Solver S(Solver::INNER_MLILUC);
			S.SetParameterEnum("levels",2);
			S.SetParameterEnum("overlap",m->GetProcessorsNumber() == 1 ? 0 : 2);
			S.SetParameterReal("rtol",1.0e-12);
			S.SetParameterReal("atol",1.0e-7);
			//S.SetParameterEnum("maxits",100);
			INMOST_DATA_REAL_TYPE tau = 0.005;
			//INMOST_DATA_REAL_TYPE tau = 0.000005;
			S.SetParameterReal("tau",tau);
			S.SetParameterReal("tau2",tau*tau);
			S.SetParameterReal("ddpq_tau",0.75);
			S.SetParameterEnum("reorder_nnz",1);
			S.SetMatrix(A);
			S.Solve(rhs,sol);
			

			tsolver = Timer() - tsolver;

			


			tupdate = Timer();


			//Mimic line search methods
			
			//this line search method calculates minimal allowed update for each unknown and then each unknown has it's own restriction
			/*
			if( d->GetDiscrType() == discr_basic::NMPFA )
			{
				Storage::real sum_p , sum_pupdate, dist;
				Tag lrtrp_tag[2] = {m->GetTag("TRP_L"), m->GetTag("TRP_R")};
				Tag lrcoef_tag[2] = {m->GetTag("COEF_L"), m->GetTag("COEF_R")};
				Tag update_coef = m->CreateTag("UPDATE_COEF",DATA_REAL,active_types, active_types & ~CELL,1);
				for(Mesh::iteratorElement it = m->BeginElement(active_types); it != m->EndElement(); ++it) if( it->GetStatus() != Element::Ghost )
				{
					if( it->GetMarker(d->UnknownMarker()) )
					{
						it->Real(pupdate) = sol[aut->GetDynamicIndex(&*it,tagp)];
						it->Real(update_coef) = 1.0;
					}
				}
				for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it)
				{
					if(it->GetMarker(d->BoundaryMarker()) && !it->GetMarker(d->UnknownMarker())) continue;
					Storage::real threshold = 0.1;
					Storage::real coef = 1.0;
					Storage::real_array lrcoef[2] = {it->RealArray(lrcoef_tag[0]), it->RealArray(lrcoef_tag[1])};
					Storage::reference_array lrtrp[2] = {it->ReferenceArray(lrtrp_tag[0]), it->ReferenceArray(lrtrp_tag[1])};
					for(int q = 0; q < 2; q++)
					{
						sum_p = sum_pupdate = 0.0;
						for(int r = 0; r < 4; r++) if( lrtrp[q][r] != NULL )
						{
							sum_p += lrcoef[q][r] * lrtrp[q][r]->Real(p);
							sum_pupdate += lrcoef[q][r] *lrtrp[q][r]->Real(pupdate);
						}
						dist = -sum_p / sum_pupdate;

						if( dist > threshold && dist < coef) coef = dist;
					}
					
					for(int q = 0; q < 2; q++)
					{
						for(int q = 0; q < 2; q++)
						{
							for(int r = 0; r < 4; r++) if( lrtrp[q][r] != NULL )
							{
								Storage::real & updcoef = lrtrp[q][r]->Real(update_coef);
								if( updcoef > coef ) 
								{
									updcoef = coef;
									//it->Real(lambda_coef) = (1.0 - updcoef)*updcoef;
								}
							}
						}
					}
					
				}

				for(Mesh::iteratorElement it = m->BeginElement(active_types); it != m->EndElement(); ++it) if( it->GetStatus() != Element::Ghost )
				{
					if( it->GetMarker(d->UnknownMarker()) )
						it->Real(p) += it->Real(update_coef) * it->Real(pupdate);
				}
			}
			
			*/
			//this line search method calculates global minimal allowed update for all unknowns
			/*
			if( d->GetDiscrType() == discr_basic::NMPFA )
			{

				update_coef = 1.0e+20;
				Storage::real sum_p[2] , sum_pupdate[2], dist[2];
				Tag lrtrp_tag[2] = {m->GetTag("TRP_L"), m->GetTag("TRP_R")};
				Tag lrcoef_tag[2] = {m->GetTag("COEF_L"), m->GetTag("COEF_R")};

				for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				{
					if( it->GetMarker(d->UnknownMarker()) && it->GetStatus() != Element::Ghost )
						it->RealDF(pupdate) = sol[it->IntegerDF(index_p)];
				}

				for(ElementSet::iterator it = additional->begin(); it != additional->end(); ++it)
					it->Real(pupdate) = sol[it->Integer(index_p)];
				int count = 0;

				for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
				{
					if(it->GetMarker(d->BoundaryMarker()) && !it->GetMarker(d->UnknownMarker())) continue;
					//if(!it->GetMarker(d->UnknownMarker())) continue;
					count++;
					Storage::real local_update_coef = 1.0e+20;
					Storage::real_array lrcoef[2] = {it->RealArrayDF(lrcoef_tag[0]), it->RealArrayDF(lrcoef_tag[1])};
					Storage::reference_array lrtrp[2] = {it->ReferenceArrayDF(lrtrp_tag[0]), it->ReferenceArrayDF(lrtrp_tag[1])};
					for(int q = 0; q < 2; q++)
					{
						sum_p[q] = sum_pupdate[q] = 0.0;
						for(int r = 0; r < 4; r++) if( lrtrp[q][r] != NULL )
						{
							sum_p[q] += lrcoef[q][r] * lrtrp[q][r]->Real(p);
							sum_pupdate[q] += lrcoef[q][r] *lrtrp[q][r]->Real(pupdate);
						}
						dist[q] = -sum_p[q] / sum_pupdate[q] * 1.05;

					}

					INMOST_DATA_REAL_TYPE threshold = 0.05;

					if( sum_p[0] * sum_p[1] > 1.0e-4 && (sum_p[0]+sum_pupdate[0])*(sum_p[1]+sum_pupdate[1]) < 0.0 )
					{
						if(dist[0] > 1.0e-9 && dist[0] < local_update_coef) local_update_coef = std::max(threshold,dist[0]);
						if(dist[1] > 1.0e-9 && dist[1] < local_update_coef) local_update_coef = std::max(threshold,dist[1]);
					}

					if( local_update_coef < update_coef ) update_coef = local_update_coef;
				}
				
				std::cout << "line search update: " << update_coef << " trust region " << trust_region_lambda << " count " << count << std::endl;
				//update_coef = std::min(update_coef,1.0);
				INMOST_DATA_REAL_TYPE active_update_coef = std::max(std::min(update_coef,1.0),0.25);
				for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) 
				{
					if( it->GetStatus() != Element::Ghost && it->GetMarker(d->UnknownMarker()))
						it->RealDF(p) += 1.0* it->RealDF(pupdate);
				}

				for(ElementSet::iterator it = additional->begin(); it != additional->end(); ++it)
				{
					it->Real(p) += 1.0 * it->Real(pupdate);
				}
				//trust_region_lambda = (1.0 - update_coef)*update_coef;// + trust_region_lambda*0.01;
				
				
				trust_region_lambda = ((1.0 - std::min(update_coef,1.0))*0.01+ trust_region_lambda*0.05)/sqrt(update_coef);
				
				
				//if( update_coef > 0.95 ) trust_region_lambda /= 100.0;
				//else trust_region_lambda *= 100.0;
				//trust_region_lambda = std::min(1.0,trust_region_lambda);
				//trust_region_lambda = std::max(1.0e-9,trust_region_lambda);
			}
			
			else
			*/
			{
				for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) 
				{
					if( it->GetStatus() != Element::Ghost && it->GetMarker(d->UnknownMarker()))
						it->RealDF(p) += sol[it->IntegerDF(index_p)];
				}

				for(ElementSet::iterator it = additional->begin(); it != additional->end(); ++it)
				{
					it->Real(p) += sol[it->Integer(index_p)];
				}
			}
			m->ExchangeData(p,active_types);

			

			tupdate = Timer() - tupdate;

			//std::cout << "matrix time " << tmatrix << " solver time " << tsolver << " update time " << tupdate << std::endl;

			iter++;

		} while(err > 1.0e-4);

		tnonlinear = Timer() - tnonlinear;

		m->Integrate(tnonlinear);
		if(m->GetProcessorRank() == 0 ) std::cout << "total time: " << tnonlinear << std::endl;

		if( pupdate.isValid() ) pupdate = m->DeleteTag(pupdate);
		//if( update_coef.isValid() ) update_coef = m->DeleteTag(update_coef);
		//if( lambda_coef.isValid() ) lambda_coef = m->DeleteTag(lambda_coef);


		//Calculate fluxes
		Tag q = m->CreateTag("FACE_FLUX",DATA_REAL, FACE,NONE,1);
		Tag vel = m->CreateTag("VELOCITY",DATA_REAL, CELL,NONE,3);
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); it++) if( d->needBuild(&*it) )
		{
			expr * Kgrad = &Kgradint;
			Element * r0 = it->BackCell();
			Element * r1 = it->FrontCell();
			if( it->GetMarker(d->BoundaryMarker()) )
			{
				//special treatment for boundary face
				if( !it->GetMarker(d->UnknownMarker()) )
				{
					it->RealDF(q) = 0.0;
					continue;
				}
				else 
				{
					r1 = &*it;
					Kgrad = &Kgradbnd;
				}
			}
			Element::Status s0 = r0->GetStatus(), s1 = r1->GetStatus();
			if( s0 == Element::Ghost && s1 == Element::Ghost ) continue;
			it->RealDF(q) = aut->Evaluate(*Kgrad,&*it,NULL);
		} else it->RealDF(q) = 0.0;


		Storage::real M[9+3], * R = M + 9, qval, ret[3], area;
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it )
		{
			Storage::real nrmf[3];
			memset(M,0,sizeof(Storage::real)*(9+3));
			adjacent<Face> faces = it->getFaces();
			for(adjacent<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
			{
				jt->OrientedNormal(&*it,nrmf);
				area = sqrt(nrmf[0]*nrmf[0] + nrmf[1]*nrmf[1]+nrmf[2]*nrmf[2]);
				qval = jt->RealDF(q) / area * (jt->FaceOrientedOutside(&*it) ? -1 : 1);
				nrmf[0] /= area;
				nrmf[1] /= area;
				nrmf[2] /= area;
				M[0] += nrmf[0]*nrmf[0];
				M[1] += nrmf[1]*nrmf[0];
				M[2] += nrmf[2]*nrmf[0];
				M[3] += nrmf[0]*nrmf[1];
				M[4] += nrmf[1]*nrmf[1];
				M[5] += nrmf[2]*nrmf[1];
				M[6] += nrmf[0]*nrmf[2];
				M[7] += nrmf[1]*nrmf[2];
				M[8] += nrmf[2]*nrmf[2];
				R[0] += qval * nrmf[0];
				R[1] += qval * nrmf[1];
				R[2] += qval * nrmf[2];
			}
			solve3x3(M,ret,R);
			memcpy(&it->RealArrayDF(vel)[0],ret,sizeof(Storage::real)*3);
		}


		if ( err > 1.0e-4 )
			std::cout << "Stopped because iteration number is 100" << std::endl;


		{
			Storage::real xyz[3];
			std::fstream fout("out.m",std::ios::out);
			fout << "p = [" << std::endl;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			{
				fout << it->Real(p) << std::endl;
			}
			fout << "];" << std::endl;
			fout << "x = [" << std::endl;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			{
				it->Centroid(xyz);
				fout << xyz[0] << std::endl;
			}
			fout << "];" << std::endl;
			fout << "y = [" << std::endl;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			{
				it->Centroid(xyz);
				fout << xyz[1] << std::endl;
			}
			fout << "];" << std::endl;
			fout << "z = [" << std::endl;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			{
				it->Centroid(xyz);
				fout << xyz[2] << std::endl;
			}
			fout << "];" << std::endl;
			fout << "vel_x = [" << std::endl;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			{
				fout << it->RealArrayDF(vel)[0] << std::endl;
			}
			fout << "];" << std::endl;
			fout << "vel_y = [" << std::endl;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			{
				fout << it->RealArrayDF(vel)[1] << std::endl;
			}
			fout << "];" << std::endl;
			fout << "vel_z = [" << std::endl;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			{
				fout << it->RealArrayDF(vel)[2] << std::endl;
			}
			fout << "];" << std::endl;
			fout.close();
		}
	
		
	}

	if( m->GetProcessorsNumber() == 1 )
	{
		m->Save("out.vtk");
		m->Save("out.gmv");
		//m->Save("out.pmf");
	}
	else
	{
		m->Save("out.pvtk");
		//m->Save("out.pmf");
	}

	
	if( d!= NULL) delete d;
	delete aut;
	delete m;

#if defined(USE_PARTITIONER)
	Partitioner::Finalize();
#endif
	Solver::Finalize();
	Mesh::Finalize();
	return 0;
}