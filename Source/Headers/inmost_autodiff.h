
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
	class Residual
	{
		Sparse::Matrix jacobian;
		Sparse::Vector residual;
		Sparse::LockService locks;
	public:
		Residual(std::string name = "", INMOST_DATA_ENUM_TYPE start = 0, INMOST_DATA_ENUM_TYPE end = 0, INMOST_MPI_Comm _comm = INMOST_MPI_COMM_WORLD)
			: jacobian(name,start,end,_comm),residual(name,start,end,_comm) {}
		Residual(const Residual & other)
			: jacobian(other.jacobian), residual(other.residual)
		{}
		Residual & operator =(Residual const & other)
		{
			jacobian = other.jacobian;
			residual = other.residual;
			return *this;
		}
		Sparse::Matrix & GetJacobian() {return jacobian;}
		const Sparse::Matrix & GetJacobian() const {return jacobian;}
		Sparse::Vector & GetResidual() {return residual;}
		const Sparse::Vector & GetResidual() const {return residual;}
		INMOST_DATA_ENUM_TYPE GetFirstIndex() const {return residual.GetFirstIndex();}
		INMOST_DATA_ENUM_TYPE GetLastIndex() const {return residual.GetLastIndex();}
		void GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const 
		{
			start = residual.GetFirstIndex(); 
			end = residual.GetLastIndex();
		}
		void SetInterval(INMOST_DATA_ENUM_TYPE beg, INMOST_DATA_ENUM_TYPE end)
		{
			jacobian.SetInterval(beg,end);
			residual.SetInterval(beg,end);
		}
		multivar_expression_reference operator [](INMOST_DATA_ENUM_TYPE row)
		{
			return multivar_expression_reference(residual[row],&jacobian[row]);
		}
		void ClearResidual()
		{
			for(Sparse::Vector::iterator it = residual.Begin(); it != residual.End(); ++it) (*it) = 0.0;
		}
		void ClearJacobian()
		{
			for(Sparse::Matrix::iterator it = jacobian.Begin(); it != jacobian.End(); ++it)
				it->Clear();
		}
		void Clear()
		{
#if defined(USE_OMP)
#pragma omp for
#endif
			for(int k = (int)GetFirstIndex(); k < (int)GetLastIndex(); ++k) 
			{
				residual[k] = 0.0;
				jacobian[k].Clear();
			}
		}
		INMOST_DATA_REAL_TYPE Norm()
		{
			INMOST_DATA_REAL_TYPE ret = 0;
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:ret)
#endif
			for(int k = (int)GetFirstIndex(); k < (int)GetLastIndex(); ++k) 
				ret += residual[k]*residual[k];
#if defined(USE_MPI)
			INMOST_DATA_REAL_TYPE tmp = ret;
			MPI_Allreduce(&tmp, &ret, 1, INMOST_MPI_DATA_REAL_TYPE, MPI_SUM, jacobian.GetCommunicator());
#endif
			return sqrt(ret);
		}
		/// Normalize entries in jacobian and right hand side
		void Rescale()
		{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(int k = (int)GetFirstIndex(); k < (int)GetLastIndex(); ++k)
			{
				INMOST_DATA_REAL_TYPE norm = 0.0;
				for(INMOST_DATA_ENUM_TYPE q = 0; q < jacobian[k].Size(); ++q)
					norm += jacobian[k].GetValue(q)*jacobian[k].GetValue(q);
				norm = sqrt(norm);
				if( norm )
				{
					norm = 1.0/norm;
					residual[k] *= norm;
					for(INMOST_DATA_ENUM_TYPE q = 0; q < jacobian[k].Size(); ++q)
						jacobian[k].GetValue(q) *= norm;
				}
			}
		}
		void InitLocks() {locks.SetInterval(GetFirstIndex(),GetLastIndex());}
		void Lock(INMOST_DATA_ENUM_TYPE pos) {if(!locks.Empty()) locks.Lock(pos);}
		void UnLock(INMOST_DATA_ENUM_TYPE pos) {if(!locks.Empty()) locks.UnLock(pos);}
		void TestLock(INMOST_DATA_ENUM_TYPE pos) {if(!locks.Empty()) locks.TestLock(pos);}
	};
#endif
	
#if defined(USE_MESH)
	/// The Automatizator class helps in defining primary unknowns of the model and
	/// enhances user experience in interaction with automatic differentiation functionality.
	/// User may declare real-typed mesh data as primary unknowns on the mesh with the
	/// function RegisterTag. The returned integer code can be used with dynamic_variable
	/// class, that is part of automatic differentiation framework helping to assemble
	/// unknowns or expressions that could be evaluated into variable on the mesh element.
	/// In addition this class provides automatic differentiation with acceleration structures
	/// that help evaluate expressions faster.
	/// \todo
	/// 1. UnRegisterTag.
	/// 2. (test) Copy constructor.
	class Automatizator
	{
	private:
		static Automatizator * CurrentAutomatizator; //< Currently used automatizator for automatic differentiation acceleration structures.
#if defined(USE_OMP)
		std::vector<Sparse::RowMerger> merger; //< Automatic differentiation acceleration structures.
#else
		Sparse::RowMerger merger; //< Automatic differentiation acceleration structures.
#endif
		typedef struct{ Tag t, indices; MarkerType domain_mask; bool active;} tagdata; //< Pair of tag and it's mask marker.
		typedef std::vector<tagdata>               tag_enum; //< A type for an array of registered tags.
		typedef std::vector<INMOST_DATA_ENUM_TYPE> del_enum; // A type for an array of deleted positions.
	private:
		std::string           name; //< Name of the automatizator.
		del_enum              del_tags; //< Array of deleted positions.
		tag_enum              reg_tags; //< Array of registered tags.
		INMOST_DATA_ENUM_TYPE first_num; //< First index in unknowns of locally owned elements.
		INMOST_DATA_ENUM_TYPE last_num; //< Last index in unknowns of locally owned elements.
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
		Sparse::RowMerger & GetMerger()
		{
#if defined(USE_OMP)
			return merger[omp_get_thread_num()];
#else
			return merger;
#endif
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
#endif
} //namespace INMOST

#endif //USE_AUTODIFF
#endif //INMOST_AUTODIFF_H_INCLUDED
