#include "inmost.h"

#if defined(USE_AUTODIFF)

namespace INMOST
{	
#if defined(USE_SOLVER)
	void Residual::GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const
	{
		start = residual.GetFirstIndex();
		end = residual.GetLastIndex();
	}
	void Residual::SetInterval(INMOST_DATA_ENUM_TYPE beg, INMOST_DATA_ENUM_TYPE end)
	{
		jacobian.SetInterval(beg,end);
		residual.SetInterval(beg,end);
	}
	void Residual::ClearResidual()
	{
		for(Sparse::Vector::iterator it = residual.Begin(); it != residual.End(); ++it) (*it) = 0.0;
	}
	void Residual::ClearJacobian()
	{
		for(Sparse::Matrix::iterator it = jacobian.Begin(); it != jacobian.End(); ++it) it->Clear();
	}
	void Residual::ClearHessian()
	{
		for(Sparse::HessianMatrix::iterator it = hessian.Begin(); it != hessian.End(); ++it) it->Clear();
	}
	void Residual::Clear()
	{
#if defined(USE_OMP)
#pragma omp for
#endif //USE_OMP
		for(int k = (int)GetFirstIndex(); k < (int)GetLastIndex(); ++k)
		{
			residual[k] = 0.0;
			if( !jacobian.Empty() ) jacobian[k].Clear();
			if( !hessian.Empty() ) hessian[k].Clear();
		}
	}
	INMOST_DATA_REAL_TYPE Residual::Norm()
	{
		INMOST_DATA_REAL_TYPE ret = 0;
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:ret)
#endif //USE_OMP
		for(int k = (int)GetFirstIndex(); k < (int)GetLastIndex(); ++k)
			ret += residual[k]*residual[k];
#if defined(USE_MPI)
		INMOST_DATA_REAL_TYPE tmp = ret;
		MPI_Allreduce(&tmp, &ret, 1, INMOST_MPI_DATA_REAL_TYPE, MPI_SUM, jacobian.GetCommunicator());
#endif
		return sqrt(ret);
	}
	Residual & Residual::operator =(Residual const & other)
	{
		hessian = other.hessian;
		jacobian = other.jacobian;
		residual = other.residual;
		return *this;
	}
	void Residual::Rescale(INMOST_DATA_ENUM_TYPE p)
	{
#if defined(USE_OMP)
#pragma omp parallel for
#endif //USE_OMP
		for(int k = (int)GetFirstIndex(); k < (int)GetLastIndex(); ++k)
		{
			INMOST_DATA_REAL_TYPE norm = 0.0;
			if( p == ENUMUNDEF ) //infinite norm
			{
				for(INMOST_DATA_ENUM_TYPE q = 0; q < jacobian[k].Size(); ++q)
					if( norm < fabs(jacobian[k].GetValue(q)) )
						norm = fabs(jacobian[k].GetValue(q));
			}
			else //p-norm
			{
				for(INMOST_DATA_ENUM_TYPE q = 0; q < jacobian[k].Size(); ++q)
					norm += pow(fabs(jacobian[k].GetValue(q)),p);
				norm = pow(norm,1.0/p);
			}
			if( norm )
			{
				norm = 1.0/norm;
				residual[k] *= norm;
				for(INMOST_DATA_ENUM_TYPE q = 0; q < jacobian[k].Size(); ++q)
					jacobian[k].GetValue(q) *= norm;
				if( !hessian.Empty() )
					for(INMOST_DATA_ENUM_TYPE q = 0; q < hessian[k].Size(); ++q)
						hessian[k].GetValue(q) *= norm;
			}
		}
	}
	Residual::Residual(std::string name, INMOST_DATA_ENUM_TYPE start, INMOST_DATA_ENUM_TYPE end, INMOST_MPI_Comm _comm)
	: hessian(name,0,0,_comm),jacobian(name,start,end,_comm),residual(name,start,end,_comm)
	{
	}
	Residual::Residual(const Residual & other)
	: hessian(other.hessian),jacobian(other.jacobian), residual(other.residual)
	{
	}
	
	Matrix<multivar_expression_reference> Residual::operator [](const AbstractMatrix<INMOST_DATA_INTEGER_TYPE> & rows)
	{
		Matrix<multivar_expression_reference> ret(rows.Rows(),rows.Cols());
		for(INMOST_DATA_ENUM_TYPE i = 0; i < rows.Rows(); ++i)
			for(INMOST_DATA_ENUM_TYPE j = 0; j < rows.Cols(); ++j)
				new (&ret(i,j)) multivar_expression_reference(residual[rows(i,j)],&jacobian[rows(i,j)]);
		return ret;
	}
	
	rpMatrix Residual::Value(const AbstractMatrix<INMOST_DATA_INTEGER_TYPE> & rows) const
	{
		rpMatrix ret(rows.Rows(),rows.Cols());
		for(INMOST_DATA_ENUM_TYPE i = 0; i < rows.Rows(); ++i)
			for(INMOST_DATA_ENUM_TYPE j = 0; j < rows.Cols(); ++j)
				ret(i,j) = residual[rows(i,j)];
		return ret;
	}
#endif //USE_SOLVER
} //namespace INMOST

#endif // USE_AUTODIFF
