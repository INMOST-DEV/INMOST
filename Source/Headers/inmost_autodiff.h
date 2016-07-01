
#ifndef INMOST_AUTODIFF_H_INCLUDED
#define INMOST_AUTODIFF_H_INCLUDED
#include "inmost_common.h"
#include "inmost_mesh.h"
#include "inmost_solver.h"
#include "inmost_variable.h"
#include <sstream> //for debug

//#define NEW_VERSION

#if defined(USE_AUTODIFF) && (!defined(USE_MESH))
#warning "USE_AUTODIFF require USE_MESH"
#undef USE_AUTODIFF
#endif

//#define DPRNT



#if defined(USE_AUTODIFF)
#include <math.h>



namespace INMOST
{
	class Automatizator; //forward declaration


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

	class Automatizator
	{
	private:
		static Automatizator * CurrentAutomatizator;
#if defined(USE_OMP)
		std::vector<Sparse::RowMerger> merger;
#else
		Sparse::RowMerger merger;
#endif
		typedef struct{ Tag t; MarkerType domain_mask; } tagdomain;
		typedef struct{ tagdomain d; Tag indices; } tagpair;
		typedef dynarray<tagpair, 128> tagpairs_type;
		typedef std::vector<tagpair> index_enum;
	private:
		index_enum            index_tags;
		tagpairs_type         reg_tags;
		INMOST_DATA_ENUM_TYPE first_num;
		INMOST_DATA_ENUM_TYPE last_num;
		Mesh * m;
	public:
		Automatizator(Mesh * m);
		~Automatizator();
		__INLINE INMOST_DATA_ENUM_TYPE                                   GetFirstIndex() { return first_num; }
		__INLINE INMOST_DATA_ENUM_TYPE                                   GetLastIndex() { return last_num; }
		/// Set data of tag t defined on domain_mask to be dynamic data.
		/// \warning
		/// Don't register tag twice.
		/// \todo
		/// Read comments inside, change merger.Resize() behavior.
		INMOST_DATA_ENUM_TYPE                                            RegisterDynamicTag(Tag t, ElementType typemask, MarkerType domain_mask = 0);
		/// Set index for every data entry of dynamic tag.
		void                                                             EnumerateDynamicTags();
		__INLINE Tag                                                     GetDynamicValueTag(INMOST_DATA_ENUM_TYPE ind) { return reg_tags[ind].d.t; }
		__INLINE Tag                                                     GetDynamicIndexTag(INMOST_DATA_ENUM_TYPE ind) { return reg_tags[ind].indices; }
		__INLINE MarkerType                                              GetDynamicMask(INMOST_DATA_ENUM_TYPE ind) { return reg_tags[ind].d.domain_mask; }
		__INLINE INMOST_DATA_REAL_TYPE                                   GetDynamicValue(const Storage & e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) { return e->RealArray(GetDynamicValueTag(ind))[comp]; }
		__INLINE INMOST_DATA_ENUM_TYPE                                   GetDynamicIndex(const Storage & e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) { return e->IntegerArray(GetDynamicIndexTag(ind))[comp]; }
		__INLINE bool                                                    isDynamicValid(const Storage & e, INMOST_DATA_ENUM_TYPE ind) { MarkerType mask = GetDynamicMask(ind); return mask == 0 || e->GetMarker(mask); }
		__INLINE INMOST_DATA_REAL_TYPE                                   GetIndex(const Storage & e, INMOST_DATA_ENUM_TYPE tagind, INMOST_DATA_ENUM_TYPE comp = 0) { return e->IntegerArray(GetDynamicIndexTag(tagind))[comp]; }
		__INLINE INMOST_DATA_ENUM_TYPE                                   GetComponents(const Storage & e, INMOST_DATA_ENUM_TYPE tagind) { return static_cast<INMOST_DATA_ENUM_TYPE>(e->IntegerArray(GetDynamicIndexTag(tagind)).size()); }
		__INLINE Mesh *                                                  GetMesh() { return m; }
		Sparse::RowMerger &                                              GetMerger() 
		{
#if defined(USE_OMP)
			return merger[omp_get_thread_num()];
#else
			return merger;
#endif
		}
		/// Remove global current automatizator.
		static void RemoveCurrent() {CurrentAutomatizator = NULL;}
		/// Set current global automatizator, so that variable will be optimized with row merger.
		static void MakeCurrent(Automatizator * aut) {CurrentAutomatizator = aut;}
		/// Check that there is an automatizator.
		static bool HaveCurrent() {return CurrentAutomatizator != NULL;}
		/// Retrive the automatizator.
		static Automatizator * GetCurrent() {return CurrentAutomatizator;}
	};
} //namespace INMOST

#endif //USE_AUTODIFF
#endif //INMOST_AUTODIFF_H_INCLUDED
