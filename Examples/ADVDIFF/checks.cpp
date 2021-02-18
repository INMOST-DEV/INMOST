#include "checks.h"

using namespace INMOST;


//shortcuts
typedef Storage::bulk       bulk;
typedef Storage::enumerator enumerator;
typedef Storage::real       real;

bool check_flux_properties(INMOST_DATA_ENUM_TYPE i, const variable & flux, const std::string & name, bool print)
{
	bool test = true;
	const Sparse::Row & r = flux.GetRow();
	for(int k = 0; k < (int)r.Size(); ++k)
	{
		if( r.GetIndex(k) == i ) //diagonal element
		{
			if( r.GetValue(k) < 0.0 )
			{
				if( print ) std::cout << "flux " << name << " gives negative diagonal contribution to " << i << "-th variable" << std::endl;
				test = false;
			}
		}
		else //off-diagonal element
		{
			if( r.GetValue(k) > 0.0 )
			{
				if( print ) std::cout << "flux " << name << " gives positive off-diagonal contribution to " << i << "-th variable" << std::endl;
				test = false;
			}
		}
	}
	if( !test && print ) 
	{
		r.Print();
		//scanf("%*c");
	}
	return test;
}


bool check_matrix_properties(const Sparse::Matrix & A, bool print)
{
	bool is_m_matrix = true;
	enumerator mbeg,mend;
	A.GetInterval(mbeg,mend);
	for(enumerator i = mbeg; i < mend; ++i)
	{
		const Sparse::Row & r = A[i];
		real rowsum = 0.0;
		for(enumerator k = 0; k < r.Size(); ++k)
		{
			rowsum += r.GetValue(k);
			if( r.GetIndex(k) != i ) 
			{
				if( r.GetValue(k) > 0.0 )
				{
					if( print ) std::cout << "Off-diagonal element (" << i << "," << r.GetIndex(k) << ") is positive: " << r.GetValue(k) << std::endl;
					is_m_matrix = false;
				}
			}
			else 
			{
				if( r.GetValue(k) < 0.0 )
				{
					if( print ) std::cout << "Diagonal element (" << i << "," << r.GetIndex(k) << ") is negative: " << r.GetValue(k) << std::endl;
					is_m_matrix = false;
				}
			}
		}
		if( rowsum < -1.0e-9 )
		{
			std::cout << "Row-sum is " << rowsum << std::endl;
			is_m_matrix = false;
		}
	}
	if( !is_m_matrix ) std::cout << "Not an M-matrix" << std::endl;
	return is_m_matrix;
}
