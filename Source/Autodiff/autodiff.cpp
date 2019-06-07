#include "inmost.h"

#if defined(USE_AUTODIFF)

#if defined(USE_AUTODIFF_ASMJIT)
#include "asmjit.h"
#endif




namespace INMOST
{
	
	template<> Demote<INMOST_DATA_REAL_TYPE>::type    AbstractEntry::Access<INMOST_DATA_REAL_TYPE>   (const Storage& e, INMOST_DATA_ENUM_TYPE pos) const {return Value(e,pos);}
	template<> Demote<INMOST_DATA_INTEGER_TYPE>::type AbstractEntry::Access<INMOST_DATA_INTEGER_TYPE>(const Storage& e, INMOST_DATA_ENUM_TYPE pos) const {return Index(e,pos);}
	template<> Demote<unknown>::type                  AbstractEntry::Access<unknown>                 (const Storage& e, INMOST_DATA_ENUM_TYPE pos) const {return Unknown(e,pos);}
	template<> Demote<variable>::type                 AbstractEntry::Access<variable>                (const Storage& e, INMOST_DATA_ENUM_TYPE pos) const {return Unknown(e,pos);}
	template<> Demote<hessian_variable>::type         AbstractEntry::Access<hessian_variable>        (const Storage& e, INMOST_DATA_ENUM_TYPE pos) const {return Unknown(e,pos);}
	
	template<>
	Matrix<Demote<INMOST_DATA_REAL_TYPE>::type, pool_array_t<Demote<INMOST_DATA_REAL_TYPE>::type> >
	AbstractEntry::Access<INMOST_DATA_REAL_TYPE>   (const Storage& e) const {return Value(e);}
	template<>
	Matrix<Demote<INMOST_DATA_INTEGER_TYPE>::type, pool_array_t<Demote<INMOST_DATA_INTEGER_TYPE>::type> >
	AbstractEntry::Access<INMOST_DATA_INTEGER_TYPE>(const Storage& e) const {return Index(e);}
	template<>
	Matrix<Demote<unknown>::type, pool_array_t<Demote<unknown>::type> >
	AbstractEntry::Access<unknown>(const Storage& e) const {return Unknown(e);}
	template<>
	Matrix<Demote<variable>::type, pool_array_t<Demote<variable>::type> >
	AbstractEntry::Access<variable>(const Storage& e) const {return Unknown(e);}
	template<>
	Matrix<Demote<hessian_variable>::type,pool_array_t<Demote<hessian_variable>::type> >
	AbstractEntry::Access<hessian_variable>(const Storage& e) const {return Unknown(e);}
	
	
#if defined(USE_MESH)
	Automatizator * Automatizator::CurrentAutomatizator = NULL;
	bool print_ad_ctor = false;
	bool GetAutodiffPrint() {return print_ad_ctor;}
	void SetAutodiffPrint(bool set) {print_ad_ctor = set;}
	bool CheckCurrentAutomatizator() {return Automatizator::HaveCurrent();}

	void FromBasicExpression(Sparse::Row & entries, const basic_expression & expr)
	{
		Sparse::RowMerger & merger = Automatizator::GetCurrent()->GetMerger();
		expr.GetJacobian(1.0,merger);
		merger.RetriveRow(entries);
		merger.Clear();
	}

	void AddBasicExpression(Sparse::Row & entries, INMOST_DATA_REAL_TYPE multme, INMOST_DATA_REAL_TYPE multit, const basic_expression & expr)
	{
		Sparse::RowMerger & merger = Automatizator::GetCurrent()->GetMerger();
		merger.PushRow(multme,entries);
		expr.GetJacobian(multit,merger);
		merger.RetriveRow(entries);
		merger.Clear();
	}

	void FromGetJacobian(const basic_expression & expr, INMOST_DATA_REAL_TYPE mult, Sparse::Row & r)
	{
		Sparse::RowMerger & merger = Automatizator::GetCurrent()->GetMerger();
		expr.GetJacobian(mult,merger);
		merger.AddRow(1.0,r);
		merger.RetriveRow(r);
		merger.Clear();
	}
#else //USE_MESH
	bool CheckCurrentAutomatizator() {return false;}
	void FromBasicExpression(Sparse::Row & entries, const basic_expression & expr) {}
	void AddBasicExpression(Sparse::Row & entries, INMOST_DATA_REAL_TYPE multme, INMOST_DATA_REAL_TYPE multit, const basic_expression & expr) {}
	void FromGetJacobian(const basic_expression & expr, INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) {}
#endif //USE_MESH
	
#if defined(USE_MESH)
	Automatizator::Automatizator(const Automatizator & b) : name(b.name+"_copy")
	{
		std::vector<INMOST_DATA_ENUM_TYPE> regs = b.ListRegisteredEntries();
		for(std::vector<INMOST_DATA_ENUM_TYPE>::iterator kt = regs.begin(); kt != regs.end(); ++kt)
			RegisterEntry(b.GetEntry(*kt));
		if( b.last_num != 0 ) EnumerateEntries();
	}
	Automatizator & Automatizator::operator =(Automatizator const & b)
	{
		if( &b != this )
		{
			name = b.name+"_copy";
			for (unsigned k = 0; k < reg_blocks.size(); k++)
				if( isRegisteredEntry(k) ) UnregisterEntry(k);
			del_blocks.clear();
			reg_blocks.clear();
			act_blocks.clear();
			std::vector<INMOST_DATA_ENUM_TYPE> regs = b.ListRegisteredEntries();
			for(std::vector<INMOST_DATA_ENUM_TYPE>::iterator kt = regs.begin(); kt != regs.end(); ++kt)
				RegisterEntry(b.GetEntry(*kt));
			if( b.last_num != 0 ) EnumerateEntries();
		}
		return *this;
	}
	Automatizator::Automatizator(std::string _name) :name(_name), first_num(0), last_num(0) {}
	Automatizator::~Automatizator()
	{
		for (unsigned k = 0; k < reg_blocks.size(); k++) 
			if( isRegisteredEntry(k) ) UnregisterEntry(k);
		del_blocks.clear();
		act_blocks.clear();
		reg_blocks.clear();
	}
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterTag(Tag t, ElementType typemask, MarkerType domain_mask, bool inverse)
	{
		if( t.GetSize() == ENUMUNDEF )
			return RegisterEntry(VectorEntry(typemask,domain_mask,inverse,t));
		else if( t.GetSize() == 1 )
			return RegisterEntry(SingleEntry(typemask,domain_mask,inverse,t,0));
		else
		{
			BlockEntry b(typemask,domain_mask,inverse);
            for(INMOST_DATA_ENUM_TYPE k = 0; k < t.GetSize(); ++k)
				b.AddTag(t,k);
			return RegisterEntry(b);
		}
	}
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterEntry(const AbstractEntry & b)
	{
		Mesh * m = b.GetMeshLink();
		
		ElementType sparse = NONE;// b.GetElementType();
		for (ElementType q = NODE; q <= MESH; q = NextElementType(q)) if (q & b.GetElementType())
		{
			for(unsigned unk = 0; unk < b.Size(); ++unk)
				sparse |= b.GetValueTag(unk).isSparse(q) ? q : NONE;
		}
		
		INMOST_DATA_ENUM_TYPE ret = ENUMUNDEF;
		if( del_blocks.empty() )
		{
			ret = static_cast<INMOST_DATA_ENUM_TYPE>(reg_blocks.size());
			reg_blocks.push_back(b.Copy());
			act_blocks.push_back(true);
		}
		else
		{
			ret = del_blocks.back();
			assert(!act_blocks[ret]);
			del_blocks.pop_back();
			reg_blocks[ret] = b.Copy();
			act_blocks[ret] = true;
		}
		
		reg_blocks[ret]->reg_index = ret;
		//b.reg_index = ret;
		
		{
			std::stringstream tag_name;
			tag_name << name << "_BLK_" << ret << "_Offset";
			reg_blocks[ret]->SetOffsetTag(m->CreateTag(tag_name.str(),DATA_INTEGER,b.GetElementType() | MESH,sparse,1));
			//b.SetOffsetTag(reg_blocks[ret]->GetOffsetTag());
		}
					
		return ret;
	}
	
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterEntry(AbstractEntry & b)
	{
		INMOST_DATA_ENUM_TYPE ret = RegisterEntry(static_cast<const AbstractEntry &>(b));
		b.reg_index = reg_blocks[ret]->GetRegistrationIndex();
		b.SetOffsetTag(reg_blocks[ret]->GetOffsetTag());
		return ret;
	}
	
	void Automatizator::UnregisterEntry(INMOST_DATA_ENUM_TYPE ind)
	{
		assert(reg_blocks[ind]);
		if( reg_blocks[ind]->GetOffsetTag().isValid() ) 
			reg_blocks[ind]->GetOffsetTag().GetMeshLink()->DeleteTag(reg_blocks[ind]->GetOffsetTag());
		delete reg_blocks[ind];
		reg_blocks[ind] = NULL;
		del_blocks.push_back(ind);
		act_blocks[ind] = false;
	}
	
	void Automatizator::DeactivateEntry(INMOST_DATA_ENUM_TYPE ind)
	{
		assert(reg_blocks[ind] != NULL); ///This block was not deleted
		AbstractEntry & b = GetEntry(ind);
		Mesh * m = b.GetMeshLink();
		b.reg_index = ENUMUNDEF;
		//std::cout << "delete " << reg_blocks[ind]->GetOffsetTag().GetTagName() << " on " << ElementTypeName(reg_blocks[ind]->GetElementType()) << std::endl;
		m->DeleteTag(b.GetOffsetTag(),b.GetElementType());
		act_blocks[ind] = false;
	}
	
	void Automatizator::ActivateEntry(INMOST_DATA_ENUM_TYPE ind)
	{
		assert(reg_blocks[ind] != NULL); ///This block was not deleted
		AbstractEntry & b = GetEntry(ind);
		Mesh * m = b.GetMeshLink();
		b.reg_index = ind;
		{
			ElementType sparse = NONE;// b.GetElementType();
			for (ElementType q = NODE; q <= MESH; q = NextElementType(q)) if (q & b.GetElementType())
			{
				for(unsigned unk = 0; unk < b.Size(); ++unk)
					sparse |= b.GetValueTag(unk).isSparse(q) ? q : NONE;
			}
			std::stringstream tag_name;
			tag_name << name << "_BLK_" << ind << "_Offset";
			b.SetOffsetTag(m->CreateTag(tag_name.str(),DATA_INTEGER,b.GetElementType(),sparse,1));
		}
		act_blocks[ind] = true;
	}
	
	void Automatizator::EnumerateEntries()
	{
		first_num = last_num = 0;
		const ElementType paralleltypes = NODE | EDGE | FACE | CELL;
		
		for (unsigned it = 0; it < reg_blocks.size(); ++it) if( act_blocks[it] )
		{
			AbstractEntry & b = *reg_blocks[it];
			TagInteger offset_tag = b.GetOffsetTag();
			Mesh * m = offset_tag.GetMeshLink();
			for (ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype))
				if (offset_tag.isDefined(etype) && offset_tag.isSparse(etype))
				{
					for(int kt = 0; kt < m->LastLocalID(etype); ++kt) if( m->isValidElement(etype,kt) )
						m->ElementByLocalID(etype,kt).DelData(offset_tag);
				}
		}

		for (unsigned it = 0; it < reg_blocks.size(); ++it) if( act_blocks[it] )
		{
			AbstractEntry & b = *reg_blocks[it];
			TagInteger offset_tag = b.GetOffsetTag();
			Mesh * m = offset_tag.GetMeshLink();
			for (ElementType etype = MESH; etype >= NODE; etype = PrevElementType(etype)) if( b.GetElementType() & etype )
			{
				for(int kt = 0; kt < m->LastLocalID(etype); ++kt) if( m->isValidElement(etype,kt) )
				{
					Element jt = m->ElementByLocalID(etype,kt);
					if ((!(etype & paralleltypes) || (jt.GetStatus() != Element::Ghost)) && b.isValid(jt) && b.Size(jt))
					{
						offset_tag[jt] = last_num;
						last_num += b.Size(jt);
					}
				}
			}
		}
		std::set<INMOST_DATA_ENUM_TYPE> Pre, Post; //Nonlocal indices
#if defined(USE_MPI)
		int size;
		MPI_Comm_size(MPI_COMM_WORLD,&size);
		if (size > 1)
		{
			MPI_Scan(&last_num, &first_num, 1, INMOST_MPI_DATA_ENUM_TYPE, MPI_SUM, MPI_COMM_WORLD);
			first_num -= last_num;
			last_num += first_num;
			ElementType exch_mask = NONE;
			for (unsigned it = 0; it < reg_blocks.size(); ++it) if( act_blocks[it] )
			{
				AbstractEntry & b = *reg_blocks[it];
				TagInteger offset_tag = b.GetOffsetTag();
				Mesh * m = offset_tag.GetMeshLink();
				for (ElementType etype = MESH; etype >= NODE; etype = PrevElementType(etype)) if( b.GetElementType() & etype )
				{
					exch_mask |= etype;
					for(int kt = 0; kt < m->LastLocalID(etype); ++kt) if( m->isValidElement(etype,kt) )
					{
						Element jt = m->ElementByLocalID(etype,kt);
						if ((!(etype & paralleltypes) || (jt.GetStatus() != Element::Ghost)) && b.isValid(jt) && b.Size(jt))
							offset_tag[jt] += first_num;
					}
				}
			}
			//synchronize indices
			{
				std::map<Mesh *,std::vector<Tag> > exch_tags;
				for (unsigned it = 0; it < reg_blocks.size(); ++it) if( act_blocks[it] )
					exch_tags[reg_blocks[it]->GetOffsetTag().GetMeshLink()].push_back(reg_blocks[it]->GetOffsetTag());
				for(std::map<Mesh *,std::vector<Tag> >::iterator it = exch_tags.begin(); it != exch_tags.end(); ++it)
					it->first->ExchangeData(it->second, exch_mask,0);
			}
			//compute out-of-bounds indices
			for (unsigned it = 0; it < reg_blocks.size(); ++it) if( act_blocks[it] )
			{
				AbstractEntry & b = *reg_blocks[it];
				TagInteger offset_tag = b.GetOffsetTag();
				Mesh * m = offset_tag.GetMeshLink();
				for (ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype)) if( b.GetElementType() & etype )
				{
					for(int kt = 0; kt < m->LastLocalID(etype); ++kt) if( m->isValidElement(etype,kt) )
					{
						Element jt = m->ElementByLocalID(etype,kt);
						if ((!(etype & paralleltypes) || (jt.GetStatus() == Element::Ghost)) && b.isValid(jt) && b.Size(jt))
						{
							for(unsigned q = 0; q < b.MatrixSize(jt); ++q) 
							{
								INMOST_DATA_ENUM_TYPE ind =  b.Index(jt,q);
								if( ind != ENUMUNDEF ) 
								{
									if( ind < first_num ) Pre.insert(ind);
									if( ind >= last_num ) Post.insert(ind);
								}
							}
						}
					}
				} //etype
			} //it
			// after cycle
		}
#endif
		//INMOST_DATA_INTEGER_TYPE max_unknowns = m->AggregateMax(static_cast<INMOST_DATA_INTEGER_TYPE>(last_num));
		//std::cout << "Proc " << m->GetProcessorRank() << " size " << last_num-first_num <<  " pre " << Pre.size() << " post " << Post.size() << " max " << max_unknowns << std::endl;
#if defined(USE_OMP)
#pragma omp parallel
#endif //USE_OMP
		{
#if defined(USE_OMP)
#pragma omp single
#endif //USE_OMP
			{
				merger.resize(MAX_THREADS);
			}
			merger[OMP_THREAD].Resize(first_num,last_num,std::vector<INMOST_DATA_ENUM_TYPE>(Pre.begin(),Pre.end()),std::vector<INMOST_DATA_ENUM_TYPE>(Post.begin(),Post.end()),false);
		}
	}
	
	std::vector<INMOST_DATA_ENUM_TYPE> Automatizator::ListRegisteredEntries() const
	{
		std::vector<INMOST_DATA_ENUM_TYPE> ret;
        for(blk_enum::size_type it = 0; it < reg_blocks.size(); ++it)
            if( isRegisteredEntry(static_cast<INMOST_DATA_ENUM_TYPE>(it)) )
                ret.push_back(static_cast<INMOST_DATA_ENUM_TYPE>(it));
		return ret;
	}
	
	void BlockEntry::AddTag(Tag value, INMOST_DATA_ENUM_TYPE comp)
	{
		assert(unknown_tags.empty() || GetMeshLink() == value.GetMeshLink());
		if( comp == ENUMUNDEF )
		{
			if( value.GetSize() != ENUMUNDEF )
			{
				for(unsigned k = 0; k < value.GetSize(); ++k)
				{
					unknown_tags.push_back(value);
					unknown_comp.push_back(k);
				}
			}
			else throw "Cannot add variable-sized tag to block";
		}
		else
		{
			unknown_tags.push_back(value);
			unknown_comp.push_back(comp);
		}
	}
	
	
	INMOST_DATA_ENUM_TYPE MultiEntry::MatrixSize(const Storage & e) const 
	{
		INMOST_DATA_ENUM_TYPE ret = 0; 
		for(unsigned k = 0; k < entries.size(); ++k) if( entries[k]->isValid(e) )
			ret += entries[k]->MatrixSize(e); 
		return ret;
	}
	
	INMOST_DATA_REAL_TYPE MultiEntry::Value(const Storage & e, INMOST_DATA_ENUM_TYPE unk) const
	{
        unsigned pos = 0;
		for(unsigned k = 0; k < entries.size(); ++k) if( entries[k]->isValid(e) )
		{
			unsigned s = entries[k]->MatrixSize(e);
			if( pos + s > unk )
				return entries[k]->Value(e,unk-pos);
			else pos += s; 
		}
		throw Impossible;
	}
	
	INMOST_DATA_REAL_TYPE & MultiEntry::Value(const Storage & e, INMOST_DATA_ENUM_TYPE unk)
	{
        unsigned pos = 0;
		for(unsigned k = 0; k < entries.size(); ++k) if( entries[k]->isValid(e) )
		{
			unsigned s = entries[k]->MatrixSize(e);
			if( pos + s > unk )
				return entries[k]->Value(e,unk-pos);
			else pos += s; 
		}
		throw Impossible;
	}
	
	INMOST_DATA_ENUM_TYPE MultiEntry::Index(const Storage & e, INMOST_DATA_ENUM_TYPE unk) const
	{
        unsigned pos = 0;
		for(unsigned k = 0; k < entries.size(); ++k) if( entries[k]->isValid(e) )
		{
			unsigned s = entries[k]->MatrixSize(e);
			if( pos + s > unk )
				return entries[k]->Index(e,unk-pos);
			else pos += s; 
		}
		throw Impossible;
	}
	
	unknown MultiEntry::Unknown(const Storage & e, INMOST_DATA_ENUM_TYPE unk) const
	{
        unsigned pos = 0;
		for(unsigned k = 0; k < entries.size(); ++k) if( entries[k]->isValid(e) )
		{
			unsigned s = entries[k]->MatrixSize(e);
			if( pos + s > unk )
				return entries[k]->Unknown(e,unk-pos);
			else pos += s; 
		}
		throw Impossible;
	}
	
	rpMatrix MultiEntry::Value(const Storage & e) const
	{
		rpMatrix ret(MatrixSize(e),1);
		unsigned l = 0, r, t;
		for(unsigned k = 0; k < entries.size(); ++k) if( entries[k]->isValid(e) )
		{
			t = entries[k]->MatrixSize(e);
			for(r = 0; r < t; ++r)
				ret(l++,0) = entries[k]->Value(e,r);
		}
		return ret;
	}
	
	ipMatrix MultiEntry::Index(const Storage & e) const
	{
		ipMatrix ret(MatrixSize(e),1);
		unsigned l = 0, r, t;
		for(unsigned k = 0; k < entries.size(); ++k) if( entries[k]->isValid(e) )
		{
			t = entries[k]->MatrixSize(e);
			for(r = 0; r < t; ++r)
				ret(l++,0) = entries[k]->Index(e,r);
		}
		return ret;
	}
	
	upMatrix MultiEntry::operator [](const Storage & e) const
	{
		upMatrix ret(MatrixSize(e),1);
		unsigned l = 0, r, t;
		for(unsigned k = 0; k < entries.size(); ++k) if( entries[k]->isValid(e) )
		{
			t = entries[k]->MatrixSize(e);
			for(r = 0; r < t; ++r)
				ret(l++,0) = entries[k]->Unknown(e,r);
		}
		return ret;
	}

	INMOST_DATA_ENUM_TYPE MultiEntry::GetValueComp(INMOST_DATA_ENUM_TYPE unk) const
	{
		unsigned pos = 0, k = 0;
		while( pos + entries[k]->Size() <= unk ) pos += entries[k++]->Size();
		assert(k < entries.size());
		return entries[k]->GetValueComp(unk-pos);
	}
	
	TagRealArray MultiEntry::GetValueTag(INMOST_DATA_ENUM_TYPE unk) const
	{
		unsigned pos = 0, k = 0;
		while( pos + entries[k]->Size() <= unk ) pos += entries[k++]->Size();
		assert(k < entries.size());
		return entries[k]->GetValueTag(unk-pos);
	}
	
	AbstractEntry * MultiEntry::Copy() const
	{
		MultiEntry * ret = new MultiEntry(GetElementType(),GetMask(),GetMaskInverse());
		for(unsigned k = 0; k < entries.size(); ++k)
			ret->entries.push_back(entries[k]->Copy());
		return ret;
	}
	
	void MultiEntry::AddEntry(const AbstractEntry & entry)
	{
		assert(entries.empty() || (GetMeshLink() == entry.GetMeshLink())); 
		SetElementType(GetElementType() | entry.GetElementType());
		entries.push_back(entry.Copy());
	}
	
	
	void StatusBlockEntry::AddTag(Tag value, INMOST_DATA_ENUM_TYPE comp)
	{
		assert(unknown_tags.empty() || GetMeshLink() == value.GetMeshLink());
		if( comp == ENUMUNDEF )
		{
			if( value.GetSize() != ENUMUNDEF )
			{
				for(unsigned k = 0; k < value.GetSize(); ++k)
				{
					unknown_tags.push_back(value);
					unknown_comp.push_back(k);
				}
			}
			else throw "Cannot add variable-sized tag to block";
		}
		else
		{
			unknown_tags.push_back(value);
			unknown_comp.push_back(comp);
		}
	}
	
	INMOST_DATA_ENUM_TYPE StatusBlockEntry::Size(const Storage & e) const 
	{
		INMOST_DATA_ENUM_TYPE ret = 0; 
		int stat = status_tag[e];
		for(unsigned k = 0; k < unknown_tags.size(); ++k) 
			if( e.HaveData(unknown_tags[k]) && status_tbl[stat][k] ) ret++;
		return ret;
	}
	
	AbstractEntry * StatusBlockEntry::Copy() const 
	{
		StatusBlockEntry * ret = new StatusBlockEntry(GetElementType(),GetMask(),GetMaskInverse(),status_tag,status_tbl);
		for(unsigned k = 0; k < Size(); ++k) 
			ret->AddTag(unknown_tags[k],unknown_comp[k]); 
		return ret; 
	}
	
	void AbstractEntry::SynchronizeData()
	{
		//synchronize indices
		{
			Mesh * m = NULL;
			std::set<Tag> exch_tags;
			for(unsigned jt = 0; jt < Size(); ++jt)
			{
				if( m == NULL )
					m = GetValueTag(jt).GetMeshLink();
				else assert(m == GetValueTag(jt).GetMeshLink());
				exch_tags.insert(GetValueTag(jt));
			}
			m->ExchangeData(std::vector<Tag>(exch_tags.begin(),exch_tags.end()), GetElementType(),0);
		}
	}
	
	
	void Automatizator::SynchronizeData()
	{
		//TODO:
		// optimize std::map usage
		// don't sort Tag by itself, only by name, otherwise memory location may affect order and result
		//synchronize indices
		{
			std::map<Mesh *,std::map<std::string,Tag> > exch_tags;
			ElementType exch_mask = NONE;
			for (unsigned it = 0; it < reg_blocks.size(); ++it) if( act_blocks[it] )
			{
				exch_mask |= reg_blocks[it]->GetElementType();
				for(unsigned jt = 0; jt < reg_blocks[it]->Size(); ++jt)
					exch_tags[reg_blocks[it]->GetValueTag(jt).GetMeshLink()][reg_blocks[it]->GetValueTag(jt).GetTagName()] = reg_blocks[it]->GetValueTag(jt);
			}
			for(std::map<Mesh *,std::map<std::string,Tag> >::iterator it = exch_tags.begin(); it != exch_tags.end(); ++it)
			{
				std::vector<Tag> sync;
				//std::cout << "Synchronize tags: ";
				for(std::map<std::string,Tag>::iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
				{
					sync.push_back(jt->second);
					//std::cout << jt->first << " ";
				}
				//std::cout << std::endl;
				
				it->first->ExchangeData(sync, exch_mask,0);
			}
		}
	}
#endif //USE_MESH
} //namespace INMOST

#endif //USE_AUTODIFF
