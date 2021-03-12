#ifndef INMOST_RESIDUAL_INCLUDED
#define INMOST_RESIDUAL_INCLUDED

#include "inmost_common.h"
#include "inmost_sparse.h"

#if defined(USE_AUTODIFF) && defined(USE_SOLVER)

namespace INMOST
{
	/// The Residual class provides a representation for array of residuals of nonlinear equations.
	/// By working with the residual class you automatically assemble right hand side and
	/// the jacobian of a nonlinear system of equation.
	/// Jacobian matrix has a sparse representation.
	/// \todo
	///  1. Extend for hessian calculation.
	class Residual
	{
		Sparse::HessianMatrix hessian; ///< Hessian matrix
		Sparse::Matrix jacobian; ///< Jacobian matrix.
		Sparse::Vector residual; ///< Right hand side vector.
		Sparse::LockService locks; ///< Array of locks for openmp shared access.
	public:
		/// Constructor
		/// @param name Name for the matrix and right hand side. Can be used to set options for linear solver.
		/// @param start First index of the equation in the local partition. Use Automatizator::GetFirstIndex.
		/// @param end Last index of the equation in the local partition. Use Automatizator::GetLastIndex.
		/// @param _comm MPI Communicator.
		Residual(std::string name = "", INMOST_DATA_ENUM_TYPE start = 0, INMOST_DATA_ENUM_TYPE end = 0, INMOST_MPI_Comm _comm = INMOST_MPI_COMM_WORLD);
		/// Copy constructor.
		/// \warning May be expensive if matrices are large.
		Residual(const Residual & other);
		/// Assignment operator.
		/// \warning May be expensive if matrices are large.
		Residual & operator =(Residual const & other);
		/// Retrive the first index of the equations in the local partition.
		INMOST_DATA_ENUM_TYPE GetFirstIndex() const {return residual.GetFirstIndex();}
		/// Retrive the last index of the equations in the local partition.
		INMOST_DATA_ENUM_TYPE GetLastIndex() const {return residual.GetLastIndex();}
		/// Retrive the first and the last indices of the equations in the local partition
		/// @param start The first index of the equations will be recorded here.
		/// @param end The last index of the equations will be recorded here.
		void GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const;
		/// Assign the new first and last indices of the equations in the local partition.
		/// @param start The new first index of the equations.
		/// @param end The new last index of the equations.
		void SetInterval(INMOST_DATA_ENUM_TYPE beg, INMOST_DATA_ENUM_TYPE end);
		/// Retrive a residual value and a jacobian row corresponding to certain equation.
		/// @param row Equation number.
		/// @return A structure that can be used in or assigned an automatic differentiation expression.
		__INLINE multivar_expression_reference operator [](INMOST_DATA_ENUM_TYPE row)
		{return multivar_expression_reference(residual[row],&jacobian[row]);}
		/// Retrive a residual value corresponding to certain equation.
		/// @param row Equation number.
		/// @return A structure that can be used in or assigned an automatic differentiation expression.
		__INLINE double Value(INMOST_DATA_ENUM_TYPE row) const {return residual[row];}
		/// Retrive a vector of entries in residual, corresponding to a set of equations.
		/// @param rows A row-vector of equation numbers.
		/// @param A structure that can be used in or assigned an automatic differentiation matrix expression.
		Matrix<multivar_expression_reference> operator [](const AbstractMatrix<INMOST_DATA_INTEGER_TYPE> & rows);
		/// Retrive a vector of entries in residual, corresponding to a set of equations.
		/// @param rows A row-vector of equation numbers.
		/// @param A structure that can be used in or assigned an automatic differentiation matrix expression.
		rMatrix Value(const AbstractMatrix<INMOST_DATA_INTEGER_TYPE> & rows) const;
		/// Retrive hessian matrix. Use in nonlinear solver.
		Sparse::HessianMatrix & GetHessian() {return hessian;}
		/// Retrive hessian matrix without right of modificaiton.
		const Sparse::HessianMatrix & GetHessian() const {return hessian;}
		/// Retrive jacobian matrix. Use in Sparse::Solver::Solve function.
		Sparse::Matrix & GetJacobian() {return jacobian;}
		/// Retrive jacobian matrix without right of modificaiton.
		const Sparse::Matrix & GetJacobian() const {return jacobian;}
		/// Retrive right hand side vector. Use in Sparse::Solver::Solve function.
		Sparse::Vector & GetResidual() {return residual;}
		/// Retrive right hand side vector without right of modification.
		const Sparse::Vector & GetResidual() const {return residual;}
		/// Zero out right hand side vector.
		void ClearResidual();
		/// Remove all entries in jacobian matrix.
		void ClearJacobian();
		/// Remove all entries in hessian matrix.
		void ClearHessian();
		/// Zero out right hand side vector and remove all entries in jacobian matrix.
		void Clear();
		/// Compute the second norm of the right hand side vector over all of the processors.
		INMOST_DATA_REAL_TYPE Norm();
		/// Normalize jacobian rows to unit p-norms and scale right hand side accordingly.
		/// Use ENUMUNDEF as p to scale to infinite-norm.
		void Rescale(INMOST_DATA_ENUM_TYPE p = 2);
		/// Initialize openmp locks.
		void InitLocks() {locks.SetInterval(GetFirstIndex(),GetLastIndex());}
		/// Lock an equation to avoid simultaneous shared access.
		/// @param pos Equation number.
		__INLINE void Lock(INMOST_DATA_ENUM_TYPE pos) {if(!locks.Empty()) locks.Lock(pos);}
		/// UnLock an equation to allow simultaneous shared access.
		/// @param pos Equation number.
		__INLINE void UnLock(INMOST_DATA_ENUM_TYPE pos) {if(!locks.Empty()) locks.UnLock(pos);}
		/// Try to lock the equation.
		/// @param pos Equation number.
		/// @return True if equation was locked.
		__INLINE bool TestLock(INMOST_DATA_ENUM_TYPE pos) {if(!locks.Empty()) return locks.TestLock(pos); return false;}
	};
} //namespace INMOST

#endif //USE_SOLVER && USE_AUTODIFF

#endif //INMOST_RESIDUAL_INCLUDED
