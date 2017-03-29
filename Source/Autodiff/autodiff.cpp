#include "inmost.h"

#if defined(USE_AUTODIFF)

#if defined(USE_AUTODIFF_ASMJIT)
#include "asmjit.h"
#endif

#if defined(USE_AUTODIFF_OPENCL)
#include <CL/cl.h>
#endif


namespace INMOST
{
	
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
		std::vector<INMOST_DATA_ENUM_TYPE> regs = b.ListRegisteredTags();
		for(std::vector<INMOST_DATA_ENUM_TYPE>::iterator kt = regs.begin(); kt != regs.end(); ++kt)
			RegisterTag(b.GetValueTag(*kt),b.GetElementType(*kt),b.GetMask(*kt));
		if( b.last_num != 0 ) EnumerateTags();
	}
	Automatizator & Automatizator::operator =(Automatizator const & b)
	{
		if( &b != this )
		{
			name = b.name+"_copy";
			del_tags.clear();
			reg_tags.clear();
			std::vector<INMOST_DATA_ENUM_TYPE> regs = b.ListRegisteredTags();
			for(std::vector<INMOST_DATA_ENUM_TYPE>::iterator kt = regs.begin(); kt != regs.end(); ++kt)
				RegisterTag(b.GetValueTag(*kt),b.GetElementType(*kt),b.GetMask(*kt));
			if( b.last_num != 0 ) EnumerateTags();
		}
		return *this;
	}
	Automatizator::Automatizator(std::string _name) :name(_name), first_num(0), last_num(0) {}
	Automatizator::~Automatizator()
	{
		del_tags.clear();
		for (unsigned k = 0; k < reg_tags.size(); k++) if( reg_tags[k].active )
			reg_tags[k].indices = reg_tags[k].indices.GetMeshLink()->DeleteTag(reg_tags[k].indices);
	}
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterTag(Tag t, ElementType typemask, MarkerType domain_mask)
	{
		tagdata p;
		p.domain_mask = domain_mask;
		p.t = t;
		ElementType def = NONE, sparse = NONE;
		for (ElementType q = NODE; q <= MESH; q = q << 1) if (q & typemask)
		{
			if (t.isDefined(q)) def |= q;
			if (t.isSparse(q)) sparse |= q;
		}
		p.indices = t.GetMeshLink()->CreateTag(t.GetTagName() + "_index_" + name, DATA_INTEGER, def, sparse, t.GetSize());
		p.active = true;
		INMOST_DATA_ENUM_TYPE ret;
		
		if( del_tags.empty() )
		{
			ret = static_cast<INMOST_DATA_ENUM_TYPE>(reg_tags.size());
			reg_tags.push_back(p);
		}
		else
		{
			ret = del_tags.back();
			assert(!reg_tags[ret].active);
			del_tags.pop_back();
			reg_tags[ret] = p;
		}
		return ret;
	}
	void Automatizator::UnregisterTag(INMOST_DATA_ENUM_TYPE ind)
	{
		assert(reg_tags[ind].active);
		del_tags.push_back(ind);
		reg_tags[ind].active = false;
	}
	void Automatizator::EnumerateTags()
	{
		first_num = last_num = 0;
		const ElementType paralleltypes = NODE | EDGE | FACE | CELL;
		for (tag_enum::iterator it = reg_tags.begin(); it != reg_tags.end(); ++it) if( it->active )
		{
			Mesh * m = it->indices.GetMeshLink();
			for (ElementType etype = NODE; etype <= MESH; etype = etype << 1)
				if (it->indices.isDefined(etype) && it->indices.isSparse(etype))
				{
					for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
						jt->DelData(it->indices);
				}
		}


		for (tag_enum::iterator it = reg_tags.begin(); it != reg_tags.end(); ++it) if( it->active )
		{
			Mesh * m = it->indices.GetMeshLink();
			for (ElementType etype = MESH; etype >= NODE; etype = PrevElementType(etype))
			{
				if (it->indices.isDefined(etype))
				{
					if (it->indices.GetSize() == ENUMUNDEF)
					{
						if (!it->indices.isSparse(etype))
						{
							for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
							{
								if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (it->domain_mask == 0 || jt->GetMarker(it->domain_mask)))
								{
									Storage::integer_array indarr = jt->IntegerArray(it->indices);
									indarr.resize(jt->RealArray(it->t).size());
									for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
										*qt = last_num++;
								}
							}
						}
						else
						{
							for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
							{
								if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (jt->HaveData(it->t) && (it->domain_mask == 0 || jt->GetMarker(it->domain_mask))))
								{
									Storage::integer_array indarr = jt->IntegerArray(it->indices);
									indarr.resize(jt->RealArray(it->t).size());
									for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
										*qt = last_num++;
								}
							}
						}
					}
					else
					{
						if (!it->indices.isSparse(etype))
						{
							for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
							{
								if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (it->domain_mask == 0 || jt->GetMarker(it->domain_mask)))
								{
									Storage::integer_array indarr = jt->IntegerArray(it->indices);
									for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
										*qt = last_num++;
								}
							}
						}
						else
						{
							for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
							{
								if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (jt->HaveData(it->t) && (it->domain_mask == 0 || jt->GetMarker(it->domain_mask))))
								{
									Storage::integer_array indarr = jt->IntegerArray(it->indices);
									for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
										*qt = last_num++;
								}
							}
						}
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
			ElementType exch_mask = NONE;
			for (tag_enum::iterator it = reg_tags.begin(); it != reg_tags.end(); ++it) if( it->active )
			{
				Mesh * m = it->indices.GetMeshLink();
				for (ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype))
				{
					if (it->indices.isDefined(etype))
					{
						exch_mask |= etype;
						if (it->indices.GetSize() == ENUMUNDEF)
						{
							if (!it->indices.isSparse(etype))
							{
								for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
								{
									if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (it->domain_mask == 0 || jt->GetMarker(it->domain_mask)))
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
											*qt += first_num;
									}
								}
							}
							else
							{
								for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
								{
									if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (jt->HaveData(it->t) && (it->domain_mask == 0 || jt->GetMarker(it->domain_mask))))
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										indarr.resize(jt->RealArray(it->t).size());
										for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
											*qt += first_num;
									}
								}
							}
						}
						else
						{
							if (!it->indices.isSparse(etype))
							{
								for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
								{
									if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (it->domain_mask == 0 || jt->GetMarker(it->domain_mask)))
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
											*qt += first_num;
									}
								}
							}
							else
							{
								for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
								{
									if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (jt->HaveData(it->t) && (it->domain_mask == 0 || jt->GetMarker(it->domain_mask))))
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
											*qt += first_num;
									}
								}
							}
						}
					}
				}
			}
			//synchronize indices
			last_num += first_num;
			{
				std::map<Mesh *,std::vector<Tag> > exch_tags;
				for (tag_enum::iterator it = reg_tags.begin(); it != reg_tags.end(); ++it) if( it->active )
					exch_tags[it->indices.GetMeshLink()].push_back(it->indices);
				for(std::map<Mesh *,std::vector<Tag> >::iterator it = exch_tags.begin(); it != exch_tags.end(); ++it)
					it->first->ExchangeData(it->second, exch_mask,0);
			}
			//compute out-of-bounds indices
			for (tag_enum::iterator it = reg_tags.begin(); it != reg_tags.end(); ++it) if( it->active )
			{
				Mesh * m = it->indices.GetMeshLink();
				for (ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype))
				{
					if (it->indices.isDefined(etype))
					{
						exch_mask |= etype;
						if (it->indices.GetSize() == ENUMUNDEF)
						{
							if (!it->indices.isSparse(etype))
							{
								for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
								{
									if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() == Element::Ghost)) && (it->domain_mask == 0 || jt->GetMarker(it->domain_mask)))
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
										{
											if( static_cast<INMOST_DATA_ENUM_TYPE>(*qt) < first_num ) Pre.insert(*qt);
											else if( static_cast<INMOST_DATA_ENUM_TYPE>(*qt) >= last_num ) Post.insert(*qt);
										}
									}
								}
							}
							else
							{
								for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
								{
									if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() == Element::Ghost)) && (jt->HaveData(it->t) && (it->domain_mask == 0 || jt->GetMarker(it->domain_mask))))
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										indarr.resize(jt->RealArray(it->t).size());
										for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
										{
											if( static_cast<INMOST_DATA_ENUM_TYPE>(*qt) < first_num ) Pre.insert(*qt);
											else if( static_cast<INMOST_DATA_ENUM_TYPE>(*qt) >= last_num ) Post.insert(*qt);
										}
									}
								}
							}
						}
						else //getsize
						{
							if (!it->indices.isSparse(etype))
							{
								for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
								{
									if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() == Element::Ghost)) && (it->domain_mask == 0 || jt->GetMarker(it->domain_mask)))
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
										{
											if( static_cast<INMOST_DATA_ENUM_TYPE>(*qt) < first_num ) Pre.insert(*qt);
											else if( static_cast<INMOST_DATA_ENUM_TYPE>(*qt) >= last_num ) Post.insert(*qt);
										}
									}
								}
							}
							else
							{
								for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
								{
									if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() == Element::Ghost)) && (jt->HaveData(it->t) && (it->domain_mask == 0 || jt->GetMarker(it->domain_mask))))
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
										{
											if( static_cast<INMOST_DATA_ENUM_TYPE>(*qt) < first_num ) Pre.insert(*qt);
											else if( static_cast<INMOST_DATA_ENUM_TYPE>(*qt) >= last_num ) Post.insert(*qt);
										}
									}
								}
							}
						} //getsize
					} //isdefined
				} //etype
			} //it
			// after cycle
		}
#endif
		// this version will fail in parallel
		//merger.Resize(first_num,last_num,false);
		// use this version until there is a way to define multiple intervals in RowMerger
		//INMOST_DATA_INTEGER_TYPE max_unknowns = m->AggregateMax(static_cast<INMOST_DATA_INTEGER_TYPE>(last_num));
		//std::cout << "Proc " << m->GetProcessorRank() << " size " << last_num-first_num <<  " pre " << Pre.size() << " post " << Post.size() << " max " << max_unknowns << std::endl;
#if defined(USE_OMP)
#pragma omp parallel
		{
#pragma omp single
			{
				merger.resize(omp_get_num_procs());
			}
			merger[omp_get_thread_num()].Resize(first_num,last_num,std::vector<INMOST_DATA_ENUM_TYPE>(Pre.begin(),Pre.end()),std::vector<INMOST_DATA_ENUM_TYPE>(Post.begin(),Post.end()),false);
		}
#else
		merger.Resize(first_num,last_num,std::vector<INMOST_DATA_ENUM_TYPE>(Pre.begin(),Pre.end()),std::vector<INMOST_DATA_ENUM_TYPE>(Post.begin(),Post.end()),false);
#endif
	}
	
	std::vector<INMOST_DATA_ENUM_TYPE> Automatizator::ListRegisteredTags() const
	{
		std::vector<INMOST_DATA_ENUM_TYPE> ret;
		for(tag_enum::size_type it = 0; it < reg_tags.size(); ++it) if( reg_tags[it].active )
			ret.push_back(static_cast<INMOST_DATA_ENUM_TYPE>(it));
		return ret;
	}
	
	ElementType Automatizator::GetElementType(INMOST_DATA_ENUM_TYPE ind) const
	{
		Tag index = GetIndexTag(ind);
		ElementType ret = NONE;
		for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype))
			ret |= (index.isDefined(etype) ? etype : NONE);
			return ret;
	}
#endif //USE_MESH
};

#endif //USE_AUTODIFF
