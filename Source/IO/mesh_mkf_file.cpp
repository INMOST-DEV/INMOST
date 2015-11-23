#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "inmost.h"

#if defined(USE_MESH)

namespace INMOST
{
  // mesh format by Mohammad Karimi-Fard
  void Mesh::LoadMKF(std::string File)
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

		//TEMPORARY
		//for(int i = 0; i < nbzones; i++) newsets[i]->Delete();

		fclose(f);
	}
}

#endif