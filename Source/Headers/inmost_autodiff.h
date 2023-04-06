
#ifndef INMOST_AUTODIFF_H_INCLUDED
#define INMOST_AUTODIFF_H_INCLUDED

#include "inmost_common.h"
#include "inmost_mesh.h"

#if defined(USE_AUTODIFF) && defined(USE_MESH)

namespace INMOST
{
	class Automatizator; //forward declaration
	
	//return specific type with specific template
	template<class T> struct Demote;
	template<> struct Demote<INMOST_DATA_INTEGER_TYPE> {typedef INMOST_DATA_INTEGER_TYPE type;};
	template<> struct Demote<INMOST_DATA_REAL_TYPE> {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Demote<unknown> {typedef unknown type;};
	template<> struct Demote<variable> {typedef unknown type;};
	template<> struct Demote<hessian_variable> {typedef unknown type;};
	
	/// This class is used to organize unknowns in abstract way,
	/// it should be registered with and managed by class Automatizator.
	/// \todo
	/// 1. Is there a need for layout on how matrices are returned?
	/// 2. Is there a need for layout on how unknowns and equations are arranged?
	/// 3. Function for update of variables.
	/// 4. Function for synchronization of variables.
	class AbstractEntry
	{
		INMOST_DATA_ENUM_TYPE              reg_index;    ///< Index of block registry with Automatizator.
		TagInteger                         offset_tag;   ///< Starting index of the entry.
		MarkerType                         mask;         ///< Allows to enable or disable the entire block.
		bool                               inverse_mask; ///< Invert marker mask
		ElementType                        etype;        ///< Type of elements on which unknowns are defined.
	public:
		AbstractEntry(ElementType etype = NONE, MarkerType mask = 0, bool inverse = false)
		: reg_index(ENUMUNDEF), offset_tag(Tag()), mask(mask), inverse_mask(inverse), etype(etype) {}
		/// Retrieve mask of the block.
		MarkerType GetMask() const {return mask;}
		/// Retrieve if the mask is inverted
		bool GetMaskInverse() const {return inverse_mask;}
		/// Set mask for the block.
		void SetMask(MarkerType _mask, bool inverse = false) {mask = _mask; inverse_mask = inverse;}
		/// Retrieve element type of the block.
		ElementType GetElementType() const {return etype;}
		/// Set element type for the block.
		void SetElementType(ElementType _etype) {etype = _etype;}
		/// Retrieve tag that stores enumeration offset on each element.
		TagInteger GetOffsetTag() const {return offset_tag;}
		/// Retrieve tag that stores enumeration offset on each element.
		void SetOffsetTag(TagInteger tag) {offset_tag = tag;}
		/// Check that the block is valid on given element.
		bool isValid(const Storage & e) const {return reg_index != ENUMUNDEF && offset_tag.isDefined(e.GetElementType()) && (e.GetElementType() & etype) && (mask == 0 || (e->GetMarker(mask) ^ inverse_mask));}
		/// Return value in vector of unknowns of the block at certain position.
		/// @param pos Position for which to extract the value, should be no larger then MatrixSize.
		virtual INMOST_DATA_REAL_TYPE Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const = 0;
		/// Return value in vector of unknowns of the block at certain position.
		/// @param pos Position for which to extract the value, should be no larger then MatrixSize.
		virtual INMOST_DATA_REAL_TYPE & Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) = 0;
		/// Return index in vector of indices of the block at certain position.
		/// The index may be ENUMUNDEF if the unknown is inactive.
		/// @param pos Position for which to extract the index, should be no larger then MatrixSize.
		virtual INMOST_DATA_ENUM_TYPE Index(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const = 0;
		/// Return unknown in vector of variables of the block at certain position.
		/// @param pos Position for which to extract the unknown, should be no larger then MatrixSize.
		virtual unknown Unknown(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const = 0;
		/// Return vector filled with references to values of unknowns of the block.
		virtual Matrix<value_reference> Value(const Storage& e) = 0;
		/// Return vector filled with values of unknowns of the block.
		virtual rMatrix Value(const Storage & e) const = 0;
		/// Return vector filled with indices of unknowns of the block.
		virtual iMatrix Index(const Storage & e) const = 0;
		/// Return vector filled with unknowns of the block with their derivatives.
		virtual uMatrix Unknown(const Storage & e) const = 0;
		/// Return vector filled with either values or indices or unknowns of the block,
		/// depending on the template parameter.
		template<typename T>
		Matrix<typename Demote<T>::type>
		Access(const Storage &e) const;
		/// Return either value or index or unknown at specified position of the block,
		/// depending on the template parameter.
		/// @param pos Position in the block, should be no larger then MatrixSize.
		template<typename T>
		typename Demote<T>::type
		Access(const Storage &e, INMOST_DATA_ENUM_TYPE pos) const;
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
		/// Retrieve component of the tag related to unknown.
		virtual INMOST_DATA_ENUM_TYPE GetValueComp(INMOST_DATA_ENUM_TYPE unk) const = 0;
		/// Retrieve tag related to unknown value.
		virtual TagRealArray GetValueTag(INMOST_DATA_ENUM_TYPE unk) const = 0;
		/// Retrieve mesh pointer.
		virtual Mesh * GetMeshLink() const = 0;
		/// Make a copy of the object
		virtual AbstractEntry * Copy() const = 0;
		/// Retrieve a registration index.
		INMOST_DATA_ENUM_TYPE GetRegistrationIndex() const {return reg_index;}
		/// Update variables  contained in block on ghost elements of the grid.
		/// For synchronization of data in all blocks see Automatizator::SynchronizeData.
		void SynchronizeData();
		/// Destructor.
		virtual ~AbstractEntry() {}
		
		friend class Automatizator; //provide registration index from inside of Automatizator
		friend class Model; //provide registration index from inside of Model
	};
	
	
	
	/// This class is used to organize unknowns into blocks,
	/// blocks enumeration are managed by class Automatizator.
	class BlockEntry : public AbstractEntry
	{
		std::vector<TagRealArray>          unknown_tags; ///< Tags for unknown values that enter the block.
		std::vector<INMOST_DATA_ENUM_TYPE> unknown_comp; ///< Component of the tag used as unknown.
	public:
		/// Default constructor.
		BlockEntry(ElementType etype = NONE, MarkerType mask = 0, bool inverse = false) : AbstractEntry(etype,mask,inverse) {}
		/// Remove all existing tags.
		void ClearTags() { unknown_tags.clear(); unknown_comp.clear(); }
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
		rMatrix Value(const Storage& e) const { rMatrix ret(MatrixSize(e), 1); for (unsigned k = 0; k < Size(); ++k) ret(k, 0) = Value(e, k); return ret; }
		/// Return vector filled with references to values of unknowns of the block.
		Matrix<value_reference> Value(const Storage & e) {Matrix<value_reference> ret(MatrixSize(e),1); for(unsigned k = 0; k < Size(); ++k) new (&ret(k,0)) value_reference(Value(e,k)); return ret; }
		/// Return vector filled with indices of unknowns of the block.
		iMatrix Index(const Storage & e) const {iMatrix ret(MatrixSize(e),1); for(unsigned k = 0; k < Size(); ++k) ret(k,0) = Index(e,k); return ret; }
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix Unknown(const Storage & e) const {return BlockEntry::operator [](e);}
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix operator [](const Storage & e) const {uMatrix ret(MatrixSize(e),1); for(unsigned k = 0; k < Size(); ++k) ret(k,0) = Unknown(e,k); return ret; }
		/// The intended size of the matrix for this entry.
        INMOST_DATA_ENUM_TYPE MatrixSize(const Storage & e) const {(void)e; return Size();}
		/// Number of tags in block.
		INMOST_DATA_ENUM_TYPE Size() const {return (INMOST_DATA_ENUM_TYPE)unknown_tags.size();}
		/// Number of entries for each tag in the block.
        INMOST_DATA_ENUM_TYPE Size(const Storage & e) const
		{
			INMOST_DATA_ENUM_TYPE ret = 0;
			for(unsigned k = 0; k < unknown_tags.size(); ++k)
				if( e.HaveData(unknown_tags[k]) )	ret++;
			return ret;
		}
		/// Retrieve component of the tag related to unknown.
		INMOST_DATA_ENUM_TYPE GetValueComp(INMOST_DATA_ENUM_TYPE unk) const {return unknown_comp[unk];}
		/// Retrieve tag related to unknown value.
		TagRealArray GetValueTag(INMOST_DATA_ENUM_TYPE unk) const {return unknown_tags[unk];}
		/// Retrieve mesh pointer.
		Mesh * GetMeshLink() const {return unknown_tags.back().GetMeshLink();}
		/// Make a copy of the object
		AbstractEntry * Copy() const {BlockEntry * ret = new BlockEntry(GetElementType(),GetMask(),GetMaskInverse()); for(unsigned k = 0; k < Size(); ++k) ret->AddTag(unknown_tags[k],unknown_comp[k]); return ret; }
	};
	/// This class is used to organize a single unknown
	class SingleEntry : public AbstractEntry
	{
		TagRealArray          unknown_tag;
		INMOST_DATA_ENUM_TYPE unknown_comp;
	public:
		///Default constructor.
		SingleEntry(ElementType etype = NONE, MarkerType mask = 0, bool inverse = false) : AbstractEntry(etype, mask, inverse) { unknown_comp = 0; }
		///Constructor with tag.
		SingleEntry(ElementType etype, MarkerType mask, bool inverse, Tag unknown_tag, INMOST_DATA_ENUM_TYPE unknown_comp = 0) : AbstractEntry(etype,mask,inverse), unknown_tag(unknown_tag), unknown_comp(unknown_comp) {}
		///Provide tag.
		void SetTag(Tag unknown_tag_in, INMOST_DATA_ENUM_TYPE unknown_comp_in = 0) {unknown_tag = unknown_tag_in; unknown_comp = unknown_comp_in;}
		/// Return value in vector of unknowns of the block at certain position.
        INMOST_DATA_REAL_TYPE Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {(void)pos; assert(pos==0); return unknown_tag[e][unknown_comp];}
		/// Return value in vector of unknowns of the block at certain position.
        INMOST_DATA_REAL_TYPE & Value(const Storage & e, INMOST_DATA_ENUM_TYPE pos) {(void)pos; assert(pos==0); return unknown_tag[e][unknown_comp];}
		/// Return index in vector of indices of the block at certain position.
        INMOST_DATA_ENUM_TYPE Index(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {(void)pos; assert(pos==0); return isValid(e) ? GetOffsetTag()[e] : ENUMUNDEF;}
		/// Return unknown in vector of variables of the block at certain position.
		unknown Unknown(const Storage & e, INMOST_DATA_ENUM_TYPE pos) const {assert(pos==0); return unknown(Value(e,pos),Index(e,pos));}
		/// Return vector filled with references to values of unknowns of the block.
		Matrix<value_reference> Value(const Storage& e) { Matrix<value_reference> ret(1, 1); new (&ret(0, 0)) value_reference(Value(e, 0)); return ret; }
		/// Return vector filled with values of unknowns of the block.
		rMatrix Value(const Storage & e) const { rMatrix ret(1,1); ret(0,0) = Value(e,0); return ret; }
		/// Return vector filled with indices of unknowns of the block.
		iMatrix Index(const Storage & e) const { iMatrix ret(1,1); ret(0,0) = Index(e,0); return ret; }
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix Unknown(const Storage & e) const {return SingleEntry::operator [](e);}
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix operator [](const Storage & e) const { uMatrix ret(1,1); ret(0,0) = Unknown(e,0); return ret; }
		/// The intended size of the matrix for this entry.
        INMOST_DATA_ENUM_TYPE MatrixSize(const Storage & e) const {(void)e; return 1;}
		/// Number of tags in block.
		INMOST_DATA_ENUM_TYPE Size() const {return 1;}
		/// Number of entries for each tag in the block.
        INMOST_DATA_ENUM_TYPE Size(const Storage & e) const {(void)e; return e.HaveData(unknown_tag) ? 1 : 0;}
		/// Retrieve component of the tag related to unknown.
        INMOST_DATA_ENUM_TYPE GetValueComp(INMOST_DATA_ENUM_TYPE unk) const {(void)unk; assert(unk == 0); return unknown_comp;}
		/// Retrieve tag related to unknown value.
        TagRealArray GetValueTag(INMOST_DATA_ENUM_TYPE unk) const {(void)unk; assert(unk == 0); return unknown_tag;}
		/// Retrieve mesh pointer.
		Mesh * GetMeshLink() const {return unknown_tag.GetMeshLink();}
		/// Make a copy of the object
		AbstractEntry * Copy() const {return new SingleEntry(GetElementType(),GetMask(),GetMaskInverse(),unknown_tag,unknown_comp);}
	};
	/// This class is used to organize multiple unknowns resided on single tag of variable or static size
	class VectorEntry : public AbstractEntry
	{
		TagRealArray          unknown_tag;
	public:
		///Default constructor.
		VectorEntry(ElementType etype = NONE, MarkerType mask = 0, bool inverse = false) : AbstractEntry(etype,mask,inverse) {}
		///Constructor with tag.
		VectorEntry(ElementType etype, MarkerType mask, bool inverse, Tag unknown_tag) : AbstractEntry(etype,mask,inverse), unknown_tag(unknown_tag) {}
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
		/// Return vector filled with references to values of unknowns of the block.
		Matrix<value_reference> Value(const Storage& e) { Matrix<value_reference> ret(MatrixSize(e), 1); for (int k = 0; k < (int)unknown_tag[e].size(); ++k) new (&ret(k, 0)) value_reference(Value(e, k)); return ret; }
		/// Return vector filled with values of unknowns of the block.
		rMatrix Value(const Storage & e) const { rMatrix ret(MatrixSize(e),1); for(int k = 0; k < (int)unknown_tag[e].size(); ++k) ret(k,0) = Value(e,k); return ret; }
		/// Return vector filled with indices of unknowns of the block.
		iMatrix Index(const Storage & e) const { iMatrix ret(MatrixSize(e),1); for(int k = 0; k < (int)unknown_tag[e].size(); ++k) ret(k,0) = Index(e,k); return ret; }
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix Unknown(const Storage & e) const {return VectorEntry::operator [](e);}
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix operator [](const Storage & e) const { uMatrix ret(MatrixSize(e),1); for(int k = 0; k < (int)unknown_tag[e].size(); ++k) ret(k,0) = Unknown(e,k); return ret; }
		/// The intended size of the matrix for this entry.
		INMOST_DATA_ENUM_TYPE MatrixSize(const Storage & e) const {return (INMOST_DATA_ENUM_TYPE)unknown_tag[e].size();}
		/// Number of tags in block.
		INMOST_DATA_ENUM_TYPE Size() const {return 1;}
		/// Number of entries for each tag in the block.
		INMOST_DATA_ENUM_TYPE Size(const Storage & e) const {return e.HaveData(unknown_tag) ? (INMOST_DATA_ENUM_TYPE)unknown_tag[e].size() : 0;}
		/// Retrieve component of the tag related to unknown.
        INMOST_DATA_ENUM_TYPE GetValueComp(INMOST_DATA_ENUM_TYPE unk) const {(void)unk; assert(unk==0); return ENUMUNDEF;}
		/// Retrieve tag related to unknown value.
        TagRealArray GetValueTag(INMOST_DATA_ENUM_TYPE unk) const {(void)unk; assert(unk==0); return unknown_tag;}
		/// Retrieve mesh pointer.
		Mesh * GetMeshLink() const {return unknown_tag.GetMeshLink();}
		/// Make a copy of the object
		AbstractEntry * Copy() const {return new VectorEntry(GetElementType(),GetMask(),GetMaskInverse(),unknown_tag);}
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
		StatusBlockEntry(ElementType etype = NONE, MarkerType mask = 0, bool inverse = false) : AbstractEntry(etype,mask,inverse), status_tag(Tag()) {}
		/// Constructor with status tag.
		StatusBlockEntry(ElementType etype, MarkerType mask, bool inverse, TagInteger status_tag) : AbstractEntry(etype,mask,inverse), status_tag(status_tag) {}
		/// Constructor with status tag and status table.
		StatusBlockEntry(ElementType etype, MarkerType mask, bool inverse, TagInteger status_tag, const std::vector< std::vector<bool> >  & status_tbl) : AbstractEntry(etype,mask,inverse), status_tag(status_tag), status_tbl(status_tbl) {}
		/// Remove all existing tags.
		void ClearTags() { unknown_tags.clear(); unknown_comp.clear(); }
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
		/// Return vector filled with references to values of unknowns of the block.
		Matrix<value_reference> Value(const Storage& e) { Matrix<value_reference> ret(MatrixSize(e), 1); for (INMOST_DATA_ENUM_TYPE k = 0; k < Size(); ++k) new (&ret(k, 0)) value_reference(Value(e, k)); return ret; }
		/// Return vector filled with values of unknowns of the block.
		rMatrix Value(const Storage & e) const {rMatrix ret(MatrixSize(e),1); for(INMOST_DATA_ENUM_TYPE k = 0; k < Size(); ++k) ret(k,0) = Value(e,k); return ret; }
		/// Return vector filled with indices of unknowns of the block.
		iMatrix Index(const Storage & e) const {iMatrix ret(MatrixSize(e),1); for(INMOST_DATA_ENUM_TYPE k = 0; k < Size(); ++k) ret(k,0) = Index(e,k); return ret; }
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix Unknown(const Storage & e) const {return StatusBlockEntry::operator [](e);}
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix operator [](const Storage & e) const {uMatrix ret(MatrixSize(e),1); for(INMOST_DATA_ENUM_TYPE k = 0; k < Size(); ++k) ret(k,0) = Unknown(e,k); return ret; }
		/// The intended size of the matrix for this entry.
        INMOST_DATA_ENUM_TYPE MatrixSize(const Storage & e) const {(void)e; return Size();}
		/// Number of tags in block.
		INMOST_DATA_ENUM_TYPE Size() const {return (INMOST_DATA_ENUM_TYPE)unknown_tags.size();}
		/// Number of entries for each tag in the block.
		INMOST_DATA_ENUM_TYPE Size(const Storage & e) const;
		/// Retrieve component of the tag related to unknown.
		INMOST_DATA_ENUM_TYPE GetValueComp(INMOST_DATA_ENUM_TYPE unk) const {return unknown_comp[unk];}
		/// Retrieve tag related to unknown value.
		TagRealArray GetValueTag(INMOST_DATA_ENUM_TYPE unk) const {return unknown_tags[unk];}
		/// Retrieve mesh pointer.
		Mesh * GetMeshLink() const {return unknown_tags.back().GetMeshLink();}
		/// Make a copy of the object
		AbstractEntry * Copy() const;
	};
	/// Stack together multiple objects of AbstractEntry class.
	/// This may help enumerate together variable-sized and constant-sized entries or blocks.
	/// Essentially can replace BlockEntry but it may be less efficient.
	/// \todo
	/// 1. Check it works
	/// 2. Check if it crashes with different combinations of entries on different element types or different masks.
	class MultiEntry : public AbstractEntry 
	{
		std::vector<AbstractEntry *> entries;
	public:
		///Default constructor.
		MultiEntry(ElementType etype = NONE, MarkerType mask = 0, bool inverse = false) : AbstractEntry(etype,mask,inverse) {}
		///Destructor.
		~MultiEntry() {for(unsigned k = 0; k < entries.size(); ++k) delete entries[k];}
		///Add entry into block of entries.
		void AddEntry(const AbstractEntry & entry);
		///Retrieve entry from block of entries.
		AbstractEntry & GetEntry(INMOST_DATA_ENUM_TYPE k) {return *entries[k];}
		///Retrieve entry from block of entries.
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
		/// Return vector filled with references to values of unknowns of the block.
		Matrix<value_reference> Value(const Storage& e);
		/// Return vector filled with values of unknowns of the block.
		rMatrix Value(const Storage & e) const;
		/// Return vector filled with indices of unknowns of the block.
		iMatrix Index(const Storage & e) const;
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix Unknown(const Storage & e) const {return MultiEntry::operator [](e);}
		/// Return vector filled with unknowns of the block with their derivatives.
		uMatrix operator [](const Storage & e) const;
		/// The intended size of the matrix for this entry.
		INMOST_DATA_ENUM_TYPE MatrixSize(const Storage & e) const;
		/// Number of tags in block.
		INMOST_DATA_ENUM_TYPE Size() const {INMOST_DATA_ENUM_TYPE ret = 0; for(unsigned k = 0; k < entries.size(); ++k) ret += entries[k]->Size(); return ret;}
		/// Number of entries for each tag in the block.
		INMOST_DATA_ENUM_TYPE Size(const Storage & e) const {INMOST_DATA_ENUM_TYPE ret = 0; for(unsigned k = 0; k < entries.size(); ++k) ret += entries[k]->Size(e); return ret;}
		/// Retrieve component of the tag related to unknown.
		INMOST_DATA_ENUM_TYPE GetValueComp(INMOST_DATA_ENUM_TYPE unk) const;
		/// Retrieve tag related to unknown value.
		TagRealArray GetValueTag(INMOST_DATA_ENUM_TYPE unk) const;
		/// Retrieve mesh pointer.
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
		/// Retrieve first index of unknowns, local to the processor.
		__INLINE INMOST_DATA_ENUM_TYPE GetFirstIndex() const { return first_num; }
		/// Retrieve last index of unknowns, local to the processor.
		__INLINE INMOST_DATA_ENUM_TYPE GetLastIndex() const { return last_num; }
		/// Set data of tag t defined on domain_mask to be dynamic data.
		/// @param t Tag of DATA_REAL that represents independent data of the model.
		/// @param typemask Element types on which that data is independent.
		/// @param domain_mask Marker that may be used to mask indepndent data on certain elements.
		/// \warning
		/// 1. Don't register tag twice.
		/// 2. Have to call Automatizator::EnumerateEntries to compute indices.
		INMOST_DATA_ENUM_TYPE RegisterTag(Tag t, ElementType typemask, MarkerType domain_mask = 0, bool inverse = false);
		/// Register block with the automatizator.
		/// Note that copy of entry is stored with Automatizator.
		/// @param b Entry that represents block of indepenedent unknowns of the model.
		/// \warning
		/// 1. If you create your entry with intention to use it after registration, then you should use the function as entry = aut.GetEntry(aut.RegisterEntry(entry));
		INMOST_DATA_ENUM_TYPE RegisterEntry(const AbstractEntry & e);
		INMOST_DATA_ENUM_TYPE RegisterEntry(AbstractEntry & e);
		/// Erase a registered tag.
		/// @param ind Integer returned from Automatizator::RegisterTag.
		/// \warning
		/// 1. Have to call Automatizator::EnumerateEntries to recompute indices.
		void UnregisterEntry(INMOST_DATA_ENUM_TYPE ind);
		/// Swith a registered tag to be non-active, in this case it's unknowns are considered to be constant.
		/// @param ind Integer returned from Automatizator::RegisterTag.
		/// \warning
		/// 1. Have to call Automatizator::EnumerateEntries to recompute indices.
		void DeactivateEntry(INMOST_DATA_ENUM_TYPE ind);
		/// Swith a registered tag to be active, in this case it's unknowns are considered to be variable.
		/// @param ind Integer returned from Automatizator::RegisterTag.
		/// \warning
		/// 1. Have to call Automatizator::EnumerateEntries to recompute indices.
		void ActivateEntry(INMOST_DATA_ENUM_TYPE ind);
		/// Set index for every data entry of dynamic tag.
		void EnumerateEntries(bool blocks = false);
		/// Check whether the tag is still registered.
		/// @param True if tag is still registered.
		__INLINE bool isRegisteredEntry(INMOST_DATA_ENUM_TYPE ind) const {return reg_blocks[ind] != NULL;}
		/// Get index of the unknown associated with the entry on element.
		INMOST_DATA_ENUM_TYPE GetIndex(const Storage & e, INMOST_DATA_ENUM_TYPE reg_index, INMOST_DATA_ENUM_TYPE pos = 0) const {return GetEntry(reg_index).Index(e,pos);}
		/// Get value of the unknown associated with the entry on element.
		INMOST_DATA_REAL_TYPE GetValue(const Storage & e, INMOST_DATA_ENUM_TYPE reg_index, INMOST_DATA_ENUM_TYPE pos = 0) const {return GetEntry(reg_index).Value(e,pos);}
		/// Get unknown associated with the entry on element.
		unknown GetUnknown(const Storage & e, INMOST_DATA_ENUM_TYPE reg_index, INMOST_DATA_ENUM_TYPE pos = 0) const {return GetEntry(reg_index).Unknown(e,pos);}
		/// Retrieve the block from automatizator by index.
		AbstractEntry & GetEntry(INMOST_DATA_ENUM_TYPE ind) {return *reg_blocks[ind];}
		/// Retrieve the block from automatizator by index.
		const AbstractEntry & GetEntry(INMOST_DATA_ENUM_TYPE ind) const {return *reg_blocks[ind];}
		/// Lists all the indices of registered tags.
		/// @return An array with indices corresponding to all registered tags.
		std::vector<INMOST_DATA_ENUM_TYPE> ListRegisteredEntries() const;
		/// Update variables  contained in all block of automatizator on ghost elements of the grid.
		/// For synchronization of data in individual blocks see AbstractEntry::SynchronizeData.
		void SynchronizeData();
	};
} //namespace INMOST

#endif //USE_AUTODIFF && USE_MESH

#endif //INMOST_AUTODIFF_H_INCLUDED
