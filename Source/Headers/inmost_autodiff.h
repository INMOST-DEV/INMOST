
#ifndef INMOST_AUTODIFF_H_INCLUDED
#define INMOST_AUTODIFF_H_INCLUDED
#include "inmost_common.h"
#include "inmost_mesh.h"
#include "inmost_solver.h"
#include "inmost_variable.h"

#if defined(USE_AUTODIFF)
#include <math.h>

namespace INMOST
{
	class Automatizator; //forward declaration

#if defined(USE_SOLVER)
	/// The Residual class provides a representation for array of residuals of nonlinear equations.
	/// By working with the residual class you automatically assemble right hand side and
	/// the jacobian of a nonlinear system of equation.
	/// Jacobian matrix has a sparse representation.
	/// \todo
	///  1. Extend for hessian calculation.
	class Residual
	{
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
		/// Retrive jacobian matrix. Use in Sparse::Solver::Solve function.
		Sparse::Matrix & GetJacobian() {return jacobian;}
		/// Retrive jacobian matrix without right of modificaiton.
		const Sparse::Matrix & GetJacobian() const {return jacobian;}
		/// Retrive right hand side vector. Use in Sparse::Solver::Solve function.
		Sparse::Vector & GetResidual() {return residual;}
		/// Retrive right hand side vector without right of modification.
		const Sparse::Vector & GetResidual() const {return residual;}
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
		__INLINE multivar_expression_reference operator [](INMOST_DATA_ENUM_TYPE row) {return multivar_expression_reference(residual[row],&jacobian[row]);}
		/// Zero out right hand side vector.
		void ClearResidual();
		/// Remove all entries in jacobian matrix.
		void ClearJacobian();
		/// Zero out right hand side vector and remove all entries in jacobian matrix.
		void Clear();
		/// Compute the second norm of the right hand side vector over all of the processors.
		INMOST_DATA_REAL_TYPE Norm();
		/// Normalize jacobian rows to unit second norms and scale right hand side accordingly.
		void Rescale();
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
#endif //USE_SOLVER
	
#if defined(USE_MESH)
	///
	/*
	class BlockEntry :public AbstractEntry
	{
		INMOST_DATA_ENUM_TYPE reg_index;
		std::vector<Tag> unknown_tags;
		std::vector<Tag> bitmask_tags;
		std::vector<INMOST_DATA_ENUM_TYPE> comp;
		std::vector<INMOST_DATA_ENUM_TYPE> eqn_pos;
		std::vector<INMOST_DATA_ENUM_TYPE> unk_pos;
		Tag index_tag;
		MarkerType mask;
		ElementType etype;
	public:
		BlockEntry() : reg_index(ENUMUNDEF) {}
		void AddTag(ElementType etype, Tag value, Tag bitmask = Tag(), MarkerType mask = 0, int eqn = -1, int unk = -1);
		rMatrix Value(const Storage & e) const;
		iMatrix Index(const Storage & e) const;
		vMatrix operator [](const Storage & e) const;
	};
	 */
	/// The Automatizator class helps in defining primary unknowns of the model and
	/// enhances user experience in interaction with automatic differentiation functionality.
	/// User may declare real-typed mesh data as primary unknowns on the mesh with the
	/// function RegisterTag. The returned integer code can be used with dynamic_variable
	/// class, that is part of automatic differentiation framework helping to assemble
	/// unknowns or expressions that could be evaluated into variable on the mesh element.
	/// In addition this class provides automatic differentiation with acceleration structures
	/// that help evaluate expressions faster.
	/// \todo
	/// 1. (test) UnRegisterTag.
	/// 2. (test) Copy constructor.
	class Automatizator
	{
	private:
		static Automatizator * CurrentAutomatizator; ///< Currently used automatizator for automatic differentiation acceleration structures.
#if defined(USE_OMP)
		std::vector<Sparse::RowMerger> merger; ///< Automatic differentiation acceleration structures.
#else
		Sparse::RowMerger merger; ///< Automatic differentiation acceleration structures.
#endif //USE_OMP
		typedef struct{ Tag t, indices; MarkerType domain_mask; bool active;} tagdata; ///< Pair of tag and it's mask marker.
		typedef std::vector<tagdata>               tag_enum; ///< A type for an array of registered tags.
		typedef std::vector<INMOST_DATA_ENUM_TYPE> del_enum; ///< A type for an array of deleted positions.
	private:
		std::string           name; ///< Name of the automatizator.
		del_enum              del_tags; ///< Array of deleted positions.
		tag_enum              reg_tags; ///< Array of registered tags.
		INMOST_DATA_ENUM_TYPE first_num; ///< First index in unknowns of locally owned elements.
		INMOST_DATA_ENUM_TYPE last_num; ///< Last index in unknowns of locally owned elements.
	public:
		/// Make a copy.
		/// \warning
		/// Calls Automatizator::EnumerateTags after copy for indices
		/// to be valid only if original was enumerated.
		Automatizator(const Automatizator & b);
		/// Performs assignment.
		/// \warning
		/// Calls Automatizator::EnumerateTags after copy for indices
		/// to be valid only if original was enumerated.
		Automatizator & operator =(Automatizator const & b);
		/// Makes an empty Automatizator.
		Automatizator(std::string name = "");
		/// Destructor for the automatizator, deletes all the tags corresponding to indices from
		/// respective meshes.
		~Automatizator();
		/// Retrive first index of unknowns, local to the processor.
		__INLINE INMOST_DATA_ENUM_TYPE GetFirstIndex() const { return first_num; }
		/// Retrive last index of unknowns, local to the processor.
		__INLINE INMOST_DATA_ENUM_TYPE GetLastIndex() const { return last_num; }
		/// Set data of tag t defined on domain_mask to be dynamic data.
		/// @param t Tag of DATA_REAL that represents independent data of the model.
		/// @param typemask Element types on which that data is independent.
		/// @param domain_mask Marker that may be used to mask indepndent data on certain elements.
		/// \warning
		/// 1. Don't register tag twice.
		/// 2. Have to call Automatizator::EnumerateTags to compute indices.
		/// \todo
		/// Read comments inside, change merger.Resize() behavior.
		INMOST_DATA_ENUM_TYPE RegisterTag(Tag t, ElementType typemask, MarkerType domain_mask = 0);
		/// Erase a registered tag.
		/// @param ind Integer returned from Automatizator::RegisterTag.
		/// \warning
		/// 1. Have to call Automatizator::EnumerateTags to recompute indices.
		void UnregisterTag(INMOST_DATA_ENUM_TYPE ind);
		/// Set index for every data entry of dynamic tag.
		void EnumerateTags();
		/// Check whether the tag is still registered.
		/// @param True if tag is still registered.
		__INLINE bool isRegistretedTag(INMOST_DATA_ENUM_TYPE ind) const {return reg_tags[ind].active;}
		/// Retrive the tag that is used to store values of the unknown on the mesh.
		/// This is the same tag as provided to Automatizator::RegisterTag.
		/// @param ind Integer returned from Automatizator::RegisterTag.
		/// @return Tag related to values.
		__INLINE Tag GetValueTag(INMOST_DATA_ENUM_TYPE ind) const { return reg_tags[ind].t; }
		/// Retrive the tag that is used to store indices of the unknown on the mesh.
		/// @param ind Integer returned from Automatizator::RegisterTag.
		/// @return Tag related to indices.
		__INLINE Tag GetIndexTag(INMOST_DATA_ENUM_TYPE ind) const { return reg_tags[ind].indices; }
		/// Retrive the mask marker of the registered tag.
		/// @param ind Integer returned from Automatizator::RegisterTag.
		/// @return Mask marker.
		__INLINE MarkerType GetMask(INMOST_DATA_ENUM_TYPE ind) const { return reg_tags[ind].domain_mask; }
		/// Retrive the type of elements for the registered tag.
		/// @param ind Integer returned from Automatizator::RegisterTag.
		/// @return Element types.
		ElementType GetElementType(INMOST_DATA_ENUM_TYPE ind) const;
		/// Retrive value of the unknown on provided mesh element.
		/// @param e Mesh element.
		/// @param ind Integer returned from Automatizator::RegisterTag.
		/// @return Value of the unknown on element.
		__INLINE INMOST_DATA_REAL_TYPE GetValue(const Storage & e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) const { return e->RealArray(GetValueTag(ind))[comp]; }
		/// Retrive index of the unknown on provided mesh element.
		/// @param e Mesh element.
		/// @param ind Integer returned from Automatizator::RegisterTag.
		/// @return Index of the unknown on element.
		__INLINE INMOST_DATA_ENUM_TYPE GetIndex(const Storage & e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) const { return e->IntegerArray(GetIndexTag(ind))[comp]; }
		/// Check that the data is defined and independent on provided mesh element.
		/// @param e Mesh element.
		/// @param ind Integer returned from Automatizator::RegisterTag.
		/// @return Returns true if the data is defined.
		__INLINE bool isValid(const Storage & e, INMOST_DATA_ENUM_TYPE ind) const { MarkerType mask = GetMask(ind); return mask == 0 || e->GetMarker(mask); }
		/// Retrive acceleration structure for automatic differentiation.
		/// This structure can be used for operations on the matrix.
		/// @return Acceleration structure.
		__INLINE Sparse::RowMerger & GetMerger()
		{
#if defined(USE_OMP)
			return merger[omp_get_thread_num()];
#else
			return merger;
#endif //USE_OMP
		}
		/// Remove global current automatizator used to set acceleration structures for automatic differentation.
		static void RemoveCurrent() {CurrentAutomatizator = NULL;}
		/// Set current global automatizator, so that variable will be optimized with row merger.
		static void MakeCurrent(Automatizator * aut) {CurrentAutomatizator = aut;}
		/// Check that there is an automatizator.
		static bool HaveCurrent() {return CurrentAutomatizator != NULL;}
		/// Retrive the automatizator.
		/// @return Currently set automatizator.
		static Automatizator * GetCurrent() {return CurrentAutomatizator;}
		/// Lists all the indices of registered tags.
		/// @return An array with indices corresponding to all registered tags.
		std::vector<INMOST_DATA_ENUM_TYPE> ListRegisteredTags() const;
	};
#endif //USE_MESH
} //namespace INMOST

#endif //USE_AUTODIFF
#endif //INMOST_AUTODIFF_H_INCLUDED
