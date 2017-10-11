
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

#if defined(USE_OMP)
#define OMP_THREAD omp_get_thread_num()
#define MAX_THREADS omp_get_max_threads()
#else //USE_OMP
#define OMP_THREAD 0
#define MAX_THREADS 1
#endif //USE_OMP

#if defined(USE_SOLVER)
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
		/// Retrive a vector of entries in residual, corresponding to a set of equations.
		/// @param rows A row-vector of equation numbers.
		/// @param A structure that can be used in or assigned an automatic differentiation matrix expression.
		Matrix<multivar_expression_reference> operator [](const AbstractMatrix<INMOST_DATA_INTEGER_TYPE> & rows);
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
#endif //USE_SOLVER
	
#if defined(USE_MESH)
	/// This class is used to organize unknowns in abstract way,
	/// it should be registered with and managed by class Automatizator.
	/// \todo
	/// 1. Is there a need for layout on how matrices are returned?
	/// 2. Is there a need for layout on how unknowns and equations are arrenged?
	class AbstractEntry
	{
		INMOST_DATA_ENUM_TYPE              reg_index;    ///< Index of block registry with Automatizator.
		TagInteger                         offset_tag;   ///< Starting index of the entry.
		MarkerType                         mask;         ///< Allows to enable or disable the entire block.
		ElementType                        etype;        ///< Type of elements on which unknowns are defined.
	public:
		AbstractEntry(ElementType etype = NONE, MarkerType mask = 0) : reg_index(ENUMUNDEF), offset_tag(Tag()), mask(mask), etype(etype) {}
		/// Retrive mask of the block.
		MarkerType GetMask() const {return mask;}
		/// Set mask for the block.
		void SetMask(MarkerType _mask) {mask = _mask;}
		/// Retrive element type of the block.
		ElementType GetElementType() const {return etype;}
		/// Set element type for the block.
		void SetElementType(ElementType _etype) {etype = _etype;}
		/// Retrive tag that stores enumeration offset on each element.
		TagInteger GetOffsetTag() const {return offset_tag;}
		/// Retrive tag that stores enumeration offset on each element.
		void SetOffsetTag(TagInteger tag) {offset_tag = tag;}
		/// Check that the block is valid on given element.
		bool isValid(const Storage & e) const {return (e.GetElementType() & etype) && (mask == 0 || e->GetMarker(mask));}
		/// Return value in vector of unknowns of the block at certain position.
		/// @param pos Position for which to extract the value, should be no larger the MatrixSize.
		virtual INMOST_DATA_REAL_TYPE Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const = 0;
		/// Return value in vector of unknowns of the block at certain position.
		/// @param pos Position for which to extract the value, should be no larger the MatrixSize.
		virtual INMOST_DATA_REAL_TYPE & Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) = 0;
		/// Return index in vector of indices of the block at certain position.
		/// The index may be ENUMUNDEF if the unknown is inactive.
		/// @param pos Position for which to extract the index, should be no larger the MatrixSize.
		virtual INMOST_DATA_ENUM_TYPE Index(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const = 0;
		/// Return unknown in vector of variables of the block at certain position.
		/// @param pos Position for which to extract the unknown, should be no larger the MatrixSize.
		virtual unknown Unknown(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const = 0;
		/// Return vector filled with values of unknowns of the block.
		virtual rMatrix Value(const Storage & e) const = 0;
		/// Return vector filled with indices of unknowns of the block.
		virtual iMatrix Index(const Storage & e) const = 0;
		/// Return vector filled with unknowns of the block with their derivatives.
		virtual uMatrix operator [](const Storage & e) const = 0;
		/// The intended size of the matrix for this entry.
		virtual INMOST_DATA_ENUM_TYPE MatrixSize(const Storage & e) const = 0;
		/// Number of tags in block.
		virtual INMOST_DATA_ENUM_TYPE Size() const = 0;
		/// Total number of entries that this block expands to on given element,
		/// count variable-sized entries.
		/// This is an important function since it determines offset of
		/// each block when enumerated by Automatizator::EnumerateEntries.
		/// See also MatrixSize.
		virtual INMOST_DATA_ENUM_TYPE Size(const Storage & e) const = 0;
		/// Retrive component of the tag related to unknown.
		virtual INMOST_DATA_ENUM_TYPE GetValueComp(INMOST_DATA_ENUM_TYPE unk) const = 0;
		/// Retrive tag related to unknown value.
		virtual TagRealArray GetValueTag(INMOST_DATA_ENUM_TYPE unk) const = 0;
		/// Retrive mesh pointer.
		virtual Mesh * GetMeshLink() const = 0;
		/// Make a copy of the object
		virtual AbstractEntry * Copy() const = 0;
		/// Retrive a registration index.
		INMOST_DATA_ENUM_TYPE GetRegistrationIndex() const {return reg_index;}
		/// Destructor.
		virtual ~AbstractEntry() {}
		
		friend class Automatizator; //provide registration index from inside of Automatizator
	};
	/// This class is used to organize unknowns into blocks,
	/// blocks enumeration are managed by class Automatizator.
	class BlockEntry : public AbstractEntry
	{
		std::vector<TagRealArray>          unknown_tags; ///< Tags for unknown values that enter the block.
		std::vector<INMOST_DATA_ENUM_TYPE> unknown_comp; ///< Component of the tag used as unknown.
	public:
		/// Default constructor.
		BlockEntry(ElementType etype = NONE, MarkerType mask = 0) : AbstractEntry(etype,mask) {}
		/// Add a component of the tag as unknown, by default all components of the tag are added.
		/// Adding all components of variable-sized tags is not supported.
		/// \warning Tags should belong to the same mesh and should be defined on the elements indicated in etype.
		void AddTag(Tag value, INMOST_DATA_ENUM_TYPE comp = ENUMUNDEF);
		/// Return value in vector of unknowns of the block at certain position.
		INMOST_DATA_REAL_TYPE Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {return unknown_tags[pos][e][unknown_comp[pos]];}
		/// Return value in vector of unknowns of the block at certain position.
		INMOST_DATA_REAL_TYPE & Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) {return unknown_tags[pos][e][unknown_comp[pos]];}
		/// Return index in vector of indices of the block at certain position.
		INMOST_DATA_ENUM_TYPE Index(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {return isValid(e) ? GetOffsetTag()[e]+pos : ENUMUNDEF;}
		/// Return unknown in vector of variables of the block at certain position.
		unknown Unknown(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {return unknown(Value(e,pos),Index(e,pos));}
		/// Return vector filled with values of unknowns of the block.
		rMatrix Value(const Storage & e) const {rMatrix ret(MatrixSize(e),1); for(unsigned k = 0; k < Size(); ++k) ret(k,0) = Value(e,k); return ret; }
		/// Return vector filled with indices of unknowns of the block.
		iMatrix Index(const Storage & e) const {iMatrix ret(MatrixSize(e),1); for(unsigned k = 0; k < Size(); ++k) ret(k,0) = Index(e,k); return ret; }
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix operator [](const Storage & e) const {uMatrix ret(MatrixSize(e),1); for(unsigned k = 0; k < Size(); ++k) ret(k,0) = Unknown(e,k); return ret; }
		/// The intended size of the matrix for this entry.
		INMOST_DATA_ENUM_TYPE MatrixSize(const Storage & e) const {return Size();}
		/// Number of tags in block.
		INMOST_DATA_ENUM_TYPE Size() const {return (INMOST_DATA_ENUM_TYPE)unknown_tags.size();}
		/// Number of entries for each tag in the block.
		INMOST_DATA_ENUM_TYPE Size(const Storage & e) const {return (INMOST_DATA_ENUM_TYPE)unknown_tags.size();}
		/// Retrive component of the tag related to unknown.
		INMOST_DATA_ENUM_TYPE GetValueComp(INMOST_DATA_ENUM_TYPE unk) const {return unknown_comp[unk];}
		/// Retrive tag related to unknown value.
		TagRealArray GetValueTag(INMOST_DATA_ENUM_TYPE unk) const {return unknown_tags[unk];}
		/// Retrive mesh pointer.
		Mesh * GetMeshLink() const {return unknown_tags.back().GetMeshLink();}
		/// Make a copy of the object
		AbstractEntry * Copy() const {BlockEntry * ret = new BlockEntry(GetElementType(),GetMask()); for(unsigned k = 0; k < Size(); ++k) ret->AddTag(unknown_tags[k],unknown_comp[k]); return ret; }
	};
	/// This class is used to organize a single unknown,
	class SingleEntry : public AbstractEntry
	{
		TagRealArray          unknown_tag;
		INMOST_DATA_ENUM_TYPE unknown_comp;
	public:
		///Default constructor.
		SingleEntry(ElementType etype = NONE, MarkerType mask = 0) : AbstractEntry(etype,mask) {}
		///Constructor with tag.
		SingleEntry(ElementType etype, MarkerType mask, Tag unknown_tag, INMOST_DATA_ENUM_TYPE unknown_comp = 0) : AbstractEntry(etype,mask), unknown_tag(unknown_tag), unknown_comp(unknown_comp) {}
		///Provide tag.
		void SetTag(Tag unknown_tag_in, INMOST_DATA_ENUM_TYPE unknown_comp_in = 0) {unknown_tag = unknown_tag_in; unknown_comp = unknown_comp_in;}
		/// Return value in vector of unknowns of the block at certain position.
		INMOST_DATA_REAL_TYPE Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {assert(pos==0); return unknown_tag[e][unknown_comp];}
		/// Return value in vector of unknowns of the block at certain position.
		INMOST_DATA_REAL_TYPE & Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) {assert(pos==0); return unknown_tag[e][unknown_comp];}
		/// Return index in vector of indices of the block at certain position.
		INMOST_DATA_ENUM_TYPE Index(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {assert(pos==0); return isValid(e) ? GetOffsetTag()[e] : ENUMUNDEF;}
		/// Return unknown in vector of variables of the block at certain position.
		unknown Unknown(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {assert(pos==0); return unknown(Value(e,pos),Index(e,pos));}
		/// Return vector filled with values of unknowns of the block.
		rMatrix Value(const Storage & e) const { rMatrix ret(1,1); ret(0,0) = Value(e,0); return ret; }
		/// Return vector filled with indices of unknowns of the block.
		iMatrix Index(const Storage & e) const { iMatrix ret(1,1); ret(0,0) = Index(e,0); return ret; }
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix operator [](const Storage & e) const { uMatrix ret(1,1); ret(0,0) = Unknown(e,0); return ret; }
		/// The intended size of the matrix for this entry.
		INMOST_DATA_ENUM_TYPE MatrixSize(const Storage & e) const {return 1;}
		/// Number of tags in block.
		INMOST_DATA_ENUM_TYPE Size() const {return 1;}
		/// Number of entries for each tag in the block.
		INMOST_DATA_ENUM_TYPE Size(const Storage & e) const {return 1;}
		/// Retrive component of the tag related to unknown.
		INMOST_DATA_ENUM_TYPE GetValueComp(INMOST_DATA_ENUM_TYPE unk) const {assert(unk == 0); return unknown_comp;}
		/// Retrive tag related to unknown value.
		TagRealArray GetValueTag(INMOST_DATA_ENUM_TYPE unk) const { assert(unk == 0); return unknown_tag;}
		/// Retrive mesh pointer.
		Mesh * GetMeshLink() const {return unknown_tag.GetMeshLink();}
		/// Make a copy of the object
		AbstractEntry * Copy() const {return new SingleEntry(GetElementType(),GetMask(),unknown_tag,unknown_comp);}
	};
	/// This class is used to organize multiple unknowns resided on single tag of variable or static size,
	class VectorEntry : public AbstractEntry
	{
		TagRealArray          unknown_tag;
	public:
		///Default constructor.
		VectorEntry(ElementType etype = NONE, MarkerType mask = 0) : AbstractEntry(etype,mask) {}
		///Constructor with tag.
		VectorEntry(ElementType etype, MarkerType mask, Tag unknown_tag) : AbstractEntry(etype,mask), unknown_tag(unknown_tag) {}
		///Provide tag.
		void SetTag(Tag unknown_tag_in) {unknown_tag = unknown_tag_in;}
		/// Return value in vector of unknowns of the block at certain position.
		INMOST_DATA_REAL_TYPE Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {assert(pos<unknown_tag[e].size()); return unknown_tag[e][pos];}
		/// Return value in vector of unknowns of the block at certain position.
		INMOST_DATA_REAL_TYPE & Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) {assert(pos<unknown_tag[e].size()); return unknown_tag[e][pos];}
		/// Return index in vector of indices of the block at certain position.
		INMOST_DATA_ENUM_TYPE Index(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {assert(pos<unknown_tag[e].size()); return isValid(e) ? GetOffsetTag()[e]+pos : ENUMUNDEF;}
		/// Return unknown in vector of variables of the block at certain position.
		unknown Unknown(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {assert(pos<unknown_tag[e].size()); return unknown(Value(e,pos),Index(e,pos));}
		/// Return vector filled with values of unknowns of the block.
		rMatrix Value(const Storage & e) const { rMatrix ret(MatrixSize(e),1); for(int k = 0; k < (int)unknown_tag[e].size(); ++k) ret(k,0) = Value(e,k); return ret; }
		/// Return vector filled with indices of unknowns of the block.
		iMatrix Index(const Storage & e) const { iMatrix ret(MatrixSize(e),1); for(int k = 0; k < (int)unknown_tag[e].size(); ++k) ret(k,0) = Index(e,k); return ret; }
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix operator [](const Storage & e) const { uMatrix ret(MatrixSize(e),1); for(int k = 0; k < (int)unknown_tag[e].size(); ++k) ret(0,0) = Unknown(e,k); return ret; }
		/// The intended size of the matrix for this entry.
		INMOST_DATA_ENUM_TYPE MatrixSize(const Storage & e) const {return (INMOST_DATA_ENUM_TYPE)unknown_tag[e].size();}
		/// Number of tags in block.
		INMOST_DATA_ENUM_TYPE Size() const {return 1;}
		/// Number of entries for each tag in the block.
		INMOST_DATA_ENUM_TYPE Size(const Storage & e) const {return (INMOST_DATA_ENUM_TYPE)unknown_tag[e].size();}
		/// Retrive component of the tag related to unknown.
		INMOST_DATA_ENUM_TYPE GetValueComp(INMOST_DATA_ENUM_TYPE unk) const {assert(unk==0); return ENUMUNDEF;}
		/// Retrive tag related to unknown value.
		TagRealArray GetValueTag(INMOST_DATA_ENUM_TYPE unk) const {assert(unk==0); return unknown_tag;}
		/// Retrive mesh pointer.
		Mesh * GetMeshLink() const {return unknown_tag.GetMeshLink();}
		/// Make a copy of the object
		AbstractEntry * Copy() const {return new VectorEntry(GetElementType(),GetMask(),unknown_tag);}
	};
	/// This class is used to organize unknowns into blocks and provides mechanisms to change activation statuses of individual unknowns,
	/// blocks enumeration are managed by class Automatizator. This is less efficient then BlockEntry for single status.
	/// \todo
	/// 1. Check it works
	class StatusBlockEntry : public AbstractEntry 
	{
		std::vector<TagRealArray>          unknown_tags; ///< Tags for unknown values that enter the block.
		std::vector<INMOST_DATA_ENUM_TYPE> unknown_comp; ///< Component of the tag used as unknown.
		TagInteger                         status_tag; ///< Integer tag that determines current status of the block, value should be no greater then number of entries in status_tbl array.
		std::vector< std::vector<bool> >   status_tbl; ///< Array of statuses for activation of unknowns, length should be equal to number of unknowns, provided by user.
	public:
		/// Default constructor.
		StatusBlockEntry(ElementType etype = NONE, MarkerType mask = 0) : AbstractEntry(etype,mask), status_tag(Tag()) {}
		/// Constructor with status tag.
		StatusBlockEntry(ElementType etype, MarkerType mask, TagInteger status_tag) : AbstractEntry(etype,mask), status_tag(status_tag) {}
		/// Constructor with status tag and status table.
		StatusBlockEntry(ElementType etype, MarkerType mask, TagInteger status_tag, const std::vector< std::vector<bool> >  & status_tbl) : AbstractEntry(etype,mask), status_tag(status_tag), status_tbl(status_tbl) {}
		/// Add a component of the tag as unknown, by default all components of the tag are added.
		/// Adding all components of variable-sized tags is not supported.
		/// \warning Tags should belong to the same mesh and should be defined on the elements indicated in etype.
		void AddTag(Tag value, INMOST_DATA_ENUM_TYPE comp = ENUMUNDEF);
		/// Add status into table.
		void AddStatus(const std::vector<bool> & stat) {status_tbl.push_back(stat);}
		/// Set status tag.
		void SetStatusTag(TagInteger input) {status_tag = input;}
		/// Return value in vector of unknowns of the block at certain position.
		INMOST_DATA_REAL_TYPE Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {return unknown_tags[pos][e][unknown_comp[pos]];}
		/// Return value in vector of unknowns of the block at certain position.
		INMOST_DATA_REAL_TYPE & Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) {return unknown_tags[pos][e][unknown_comp[pos]];}
		/// Return index in vector of indices of the block at certain position.
		INMOST_DATA_ENUM_TYPE Index(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {return isValid(e) && status_tbl[status_tag[e]][pos] ? GetOffsetTag()[e]+pos : ENUMUNDEF;}
		/// Return unknown in vector of variables of the block at certain position.
		unknown Unknown(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {return unknown(Value(e,pos),Index(e,pos));}
		/// Return vector filled with values of unknowns of the block.
		rMatrix Value(const Storage & e) const {rMatrix ret(MatrixSize(e),1); for(INMOST_DATA_ENUM_TYPE k = 0; k < Size(); ++k) ret(k,0) = Value(e,k); return ret; }
		/// Return vector filled with indices of unknowns of the block.
		iMatrix Index(const Storage & e) const {iMatrix ret(MatrixSize(e),1); for(INMOST_DATA_ENUM_TYPE k = 0; k < Size(); ++k) ret(k,0) = Index(e,k); return ret; }
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix operator [](const Storage & e) const {uMatrix ret(MatrixSize(e),1); for(INMOST_DATA_ENUM_TYPE k = 0; k < Size(); ++k) ret(k,0) = Unknown(e,k); return ret; }
		/// The intended size of the matrix for this entry.
		INMOST_DATA_ENUM_TYPE MatrixSize(const Storage & e) const {return Size();}
		/// Number of tags in block.
		INMOST_DATA_ENUM_TYPE Size() const {return (INMOST_DATA_ENUM_TYPE)unknown_tags.size();}
		/// Number of entries for each tag in the block.
		INMOST_DATA_ENUM_TYPE Size(const Storage & e) const;
		/// Retrive component of the tag related to unknown.
		INMOST_DATA_ENUM_TYPE GetValueComp(INMOST_DATA_ENUM_TYPE unk) const {return unknown_comp[unk];}
		/// Retrive tag related to unknown value.
		TagRealArray GetValueTag(INMOST_DATA_ENUM_TYPE unk) const {return unknown_tags[unk];}
		/// Retrive mesh pointer.
		Mesh * GetMeshLink() const {return unknown_tags.back().GetMeshLink();}
		/// Make a copy of the object
		AbstractEntry * Copy() const;
	};
	/// Stack together multiple objects of AbstractEntry class.
	/// This may help enumerate together variable-sized and constant-sized entries or blocks.
	/// Essentially can replace BlockEntry but may proove to be less efficient.
	/// \todo
	/// 1. Check it works
	/// 2. Check if it crashes with different combinations of entries on different element types or different masks.
	class MultiEntry : public AbstractEntry 
	{
		std::vector<AbstractEntry *> entries;
	public:
		///Default constructor.
		MultiEntry(ElementType etype = NONE, MarkerType mask = 0) : AbstractEntry(etype,mask) {}
		///Destructor.
		~MultiEntry() {for(unsigned k = 0; k < entries.size(); ++k) delete entries[k];}
		///Add entry into block of entries.
		void AddEntry(const AbstractEntry & entry);
		///Retirve entry from block of entries.
		AbstractEntry & GetEntry(INMOST_DATA_ENUM_TYPE k) {return *entries[k];}
		///Retirve entry from block of entries.
		const AbstractEntry & GetEntry(INMOST_DATA_ENUM_TYPE k) const {return *entries[k];}
		///Total number of entries.
		INMOST_DATA_ENUM_TYPE NumEntries() const {return (INMOST_DATA_ENUM_TYPE)entries.size();}
		/// Return value in vector of unknowns of the block at certain position.
		INMOST_DATA_REAL_TYPE Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const;
		/// Return value in vector of unknowns of the block at certain position.
		INMOST_DATA_REAL_TYPE & Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos);
		/// Return index in vector of indices of the block at certain position.
		INMOST_DATA_ENUM_TYPE Index(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const;
		/// Return unknown in vector of variables of the block at certain position.
		unknown Unknown(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const;
		/// Return vector filled with values of unknowns of the block.
		rMatrix Value(const Storage & e) const;
		/// Return vector filled with indices of unknowns of the block.
		iMatrix Index(const Storage & e) const;
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix operator [](const Storage & e) const;
		/// The intended size of the matrix for this entry.
		INMOST_DATA_ENUM_TYPE MatrixSize(const Storage & e) const;
		/// Number of tags in block.
		INMOST_DATA_ENUM_TYPE Size() const {INMOST_DATA_ENUM_TYPE ret = 0; for(unsigned k = 0; k < entries.size(); ++k) ret += entries[k]->Size(); return ret;}
		/// Number of entries for each tag in the block.
		INMOST_DATA_ENUM_TYPE Size(const Storage & e) const {INMOST_DATA_ENUM_TYPE ret = 0; for(unsigned k = 0; k < entries.size(); ++k) ret += entries[k]->Size(e); return ret;}
		/// Retrive component of the tag related to unknown.
		INMOST_DATA_ENUM_TYPE GetValueComp(INMOST_DATA_ENUM_TYPE unk) const;
		/// Retrive tag related to unknown value.
		TagRealArray GetValueTag(INMOST_DATA_ENUM_TYPE unk) const;
		/// Retrive mesh pointer.
		Mesh * GetMeshLink() const {assert(!entries.empty()); return entries.back()->GetMeshLink();}
		/// Make a copy of the object
		AbstractEntry * Copy() const;
	};
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
	/// 3. Rename
	class Automatizator
	{
	private:
		static Automatizator * CurrentAutomatizator; ///< Currently used automatizator for automatic differentiation acceleration structures.
		std::vector<Sparse::RowMerger> merger; ///< Automatic differentiation acceleration structures.
		typedef std::vector<AbstractEntry *>       blk_enum; ///< A type for an array of registered tags.
		typedef std::vector<INMOST_DATA_ENUM_TYPE> del_enum; ///< A type for an array of deleted positions.
		typedef std::vector<bool>                  act_enum; ///< A type for an array of deactivation flags.
	private:
		std::string           name; ///< Name of the automatizator.
		del_enum              del_blocks; ///< Array of deleted positions.
		blk_enum              reg_blocks; ///< Array of registered blocks.
		act_enum              act_blocks; ///< Deactivated blocks since they were deleted
		INMOST_DATA_ENUM_TYPE first_num; ///< First index in unknowns of locally owned elements.
		INMOST_DATA_ENUM_TYPE last_num; ///< Last index in unknowns of locally owned elements.
	public:
		/// Make a copy.
		/// \warning
		/// Calls Automatizator::EnumerateEntries after copy for indices
		/// to be valid only if original was enumerated.
		Automatizator(const Automatizator & b);
		/// Performs assignment.
		/// \warning
		/// Calls Automatizator::EnumerateEntries after copy for indices
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
		/// 2. Have to call Automatizator::EnumerateEntries to compute indices.
		INMOST_DATA_ENUM_TYPE RegisterTag(Tag t, ElementType typemask, MarkerType domain_mask = 0);
		/// Register block with the automatizator.
		/// Note that copy of entry is stored with Automatizator.
		/// @param b Entry that represents block of indepenedent unknowns of the model.
		/// \warning
		/// 1. If you create your entry with intention to use it after registration, then you should use the function as entry = aut.GetEntry(aut.RegisterEntry(entry));
		INMOST_DATA_ENUM_TYPE RegisterEntry(const AbstractEntry & e);
		/// Erase a registered tag.
		/// @param ind Integer returned from Automatizator::RegisterTag.
		/// \warning
		/// 1. Have to call Automatizator::EnumerateEntries to recompute indices.
		void UnregisterEntry(INMOST_DATA_ENUM_TYPE ind);
		/// Set index for every data entry of dynamic tag.
		void EnumerateEntries();
		/// Check whether the tag is still registered.
		/// @param True if tag is still registered.
		__INLINE bool isRegisteredEntry(INMOST_DATA_ENUM_TYPE ind) const {return act_blocks[ind];}
		/// Get index of the unknown associated with the entry on element.
		INMOST_DATA_ENUM_TYPE GetIndex(const Storage & e, INMOST_DATA_ENUM_TYPE reg_index, INMOST_DATA_ENUM_TYPE pos = 0) const {return GetEntry(reg_index).Index(e,pos);}
		/// Get value of the unknown associated with the entry on element.
		INMOST_DATA_REAL_TYPE GetValue(const Storage & e, INMOST_DATA_ENUM_TYPE reg_index, INMOST_DATA_ENUM_TYPE pos = 0) const {return GetEntry(reg_index).Value(e,pos);}
		/// Get unknown associated with the entry on element.
		unknown GetUnknown(const Storage & e, INMOST_DATA_ENUM_TYPE reg_index, INMOST_DATA_ENUM_TYPE pos = 0) const {return GetEntry(reg_index).Unknown(e,pos);}
		/// Retive the block from automatizator by index.
		AbstractEntry & GetEntry(INMOST_DATA_ENUM_TYPE ind) {return *reg_blocks[ind];}
		/// Retive the block from automatizator by index.
		const AbstractEntry & GetEntry(INMOST_DATA_ENUM_TYPE ind) const {return *reg_blocks[ind];}
		/// @return Acceleration structure.
		__INLINE Sparse::RowMerger & GetMerger() {return merger[OMP_THREAD];}
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
		std::vector<INMOST_DATA_ENUM_TYPE> ListRegisteredEntries() const;
	};
#endif //USE_MESH
} //namespace INMOST

#endif //USE_AUTODIFF
#endif //INMOST_AUTODIFF_H_INCLUDED
