#include "inmost.h"

using namespace INMOST;



int main(int argc, char ** argv)
{
	if (argc < 9)
	{
		std::cout << "Usage: " << argv[0] << " mesh x1 y1 z1 x2 y2 z2 data:component" << std::endl;
		return -1;
	}
	
	Mesh m;
	m.Load(argv[1]);

	double p1[3], p2[3];
	p1[0] = atof(argv[2]);
	p1[1] = atof(argv[3]);
	p1[2] = atof(argv[4]);
	
	p2[0] = atof(argv[5]);
	p2[1] = atof(argv[6]);
	p2[2] = atof(argv[7]);
	
	std::string odata = std::string(argv[8]);
	size_t pos = odata.find(":");
	int comp = -1; //output all components
	if( pos != std::string::npos )
	{
		comp = atoi(odata.substr(pos+1).c_str());
		odata = odata.substr(0,pos);
	}
	/*
	std::cout << "output ";
	if( comp == -1 )
		std::cout << "all components";
	else
		std::cout << "component " << comp;
	std::cout << " of " << odata << std::endl;
	*/
	if( !m.HaveTag(odata) )
	{
		std::cout << "data " << odata << " not found in mesh" << std::endl;
		return -1;
	}
	
	Tag otag = m.GetTag(odata);
	
	if( comp >= 0 && otag.GetSize() < comp )
	{
		std::cout << "data " << odata << " size is " << otag.GetSize() << " which is less then " << comp << std::endl;
		return -1;
	}
		
	SearchKDTree tree(&m);
	ElementArray<Cell> cells;
	tree.IntersectSegment(cells,p1,p2);
	if( cells.empty() ) 
	{
		std::cout << "no cells found" << std::endl;
		return -1;
	}
	std::cout << "x; y; z; ";
	if( comp != -1 )
		std::cout << odata << "[" << comp << "]; ";
	else if( otag.GetSize() != ENUMUNDEF )
		for(unsigned k = 0; k < otag.GetSize(); ++k) std::cout << odata << "[" << k << "]; ";
	else std::cout << odata << "; ";
	std::cout << std::endl;
	
	
	for(ElementArray<Cell>::iterator it = cells.begin(); it != cells.end(); ++it)
	{
		double cnt[3];
		it->Centroid(cnt);
		std::cout << cnt[0] << "; " << cnt[1] << "; " << cnt[2] << "; ";
		if( otag.GetDataType() == DATA_REAL )
		{
			Storage::real_array oarr = it->RealArray(otag);
			if( comp == -1 ) for(Storage::real_array::iterator jt = oarr.begin(); jt != oarr.end(); ++jt) std::cout << *jt << "; ";
			else if( comp < oarr.size() ) std::cout << oarr[comp] << "; ";
			else std::cout << "NAN; ";
		}
#if defined(USE_AUTODIFF)
		else if( otag.GetDataType() == DATA_VARIABLE )
		{
			Storage::var_array oarr = it->VariableArray(otag);
			if( comp == -1 ) for(Storage::var_array::iterator jt = oarr.begin(); jt != oarr.end(); ++jt) std::cout << get_value(*jt) << "; ";
			else if( comp < oarr.size() ) std::cout << get_value(oarr[comp]) << "; ";
			else std::cout << "NAN; ";
		}
#endif
		else if( otag.GetDataType() == DATA_INTEGER )
		{
			Storage::integer_array oarr = it->IntegerArray(otag);
			if( comp == -1 ) for(Storage::integer_array::iterator jt = oarr.begin(); jt != oarr.end(); ++jt) std::cout << *jt << "; ";
			else if( comp < oarr.size() ) std::cout << oarr[comp] << "; ";
			else std::cout << "NAN; ";
		}
		else if( otag.GetDataType() == DATA_BULK )
		{
			Storage::bulk_array oarr = it->BulkArray(otag);
			if( comp == -1 ) for(Storage::bulk_array::iterator jt = oarr.begin(); jt != oarr.end(); ++jt) std::cout << *jt << "; ";
			else if( comp < oarr.size() ) std::cout << oarr[comp] << "; ";
			else std::cout << "NAN; ";
		}
		std::cout << std::endl;
	}
	
	return 0;
}
