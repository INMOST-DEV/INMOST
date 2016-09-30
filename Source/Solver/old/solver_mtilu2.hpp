
#ifndef __SOLVER_MTILU2__
#define __SOLVER_MTILU2__
#include <iomanip>

#include "inmost_solver.h"
#include "solver_prototypes.hpp"
//#define REPORT_ILU
//#define REPORT_ILU_PROGRESS
//#define REPORT_ILU_SUMMARY
using namespace INMOST;



class MTILU2_preconditioner : public Method
{
private:
	
	Sparse::Matrix * Alink;
	Solver::OrderInfo * info;
	//Sparse::Matrix L,U;
	//Sparse::Vector div;
	std::vector<INMOST_DATA_REAL_TYPE> luv;
	std::vector<INMOST_DATA_ENUM_TYPE> lui;
	interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> ilu,iu;
	interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE > Perm;
	interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE > temporary; //used for reordering
	INMOST_DATA_ENUM_TYPE Lfill;
	INMOST_DATA_REAL_TYPE tau, tau2;
	INMOST_DATA_ENUM_TYPE nnz, sciters;
	bool init;
public:

	void DumpMatrix(interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & Address, 
									std::vector<Sparse::Row::entry> & Entries,
									INMOST_DATA_ENUM_TYPE wmbeg, INMOST_DATA_ENUM_TYPE wmend,
									std::string file_name);
	INMOST_DATA_REAL_TYPE & RealParameter(std::string name);
	INMOST_DATA_ENUM_TYPE & EnumParameter(std::string name);
	MTILU2_preconditioner(Solver::OrderInfo & info);
	bool Initialize();
	bool isInitialized();
	bool Finalize();
	bool isFinalized();
	~MTILU2_preconditioner();
	void Copy(const Method * other);
	MTILU2_preconditioner(const MTILU2_preconditioner & other);
	MTILU2_preconditioner & operator =(MTILU2_preconditioner const & other);
	bool Solve(Sparse::Vector & input, Sparse::Vector & output);
	bool ReplaceMAT(Sparse::Matrix & A) { if (isInitialized()) Finalize();  Alink = &A; return true; };
	bool ReplaceSOL(Sparse::Vector & x) {(void)x;return true;}
	bool ReplaceRHS(Sparse::Vector & b) {(void)b;return true;}
	Method * Duplicate() { return new MTILU2_preconditioner(*this); }
};




#endif //__SOLVER_ILU2__
