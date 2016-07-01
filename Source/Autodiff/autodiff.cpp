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
	Automatizator::Automatizator(Mesh * m) :first_num(0), last_num(0), m(m) {}
	Automatizator::~Automatizator()
	{
		for (unsigned k = 0; k < index_tags.size(); k++)
			index_tags[k].indices = m->DeleteTag(index_tags[k].indices);
	}
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterDynamicTag(Tag t, ElementType typemask, MarkerType domain_mask)
	{
		tagpair p;
		p.d.domain_mask = domain_mask;
		p.d.t = t;
		ElementType def = NONE, sparse = NONE;
		for (ElementType q = NODE; q <= MESH; q = q << 1) if (q & typemask)
		{
			if (t.isDefined(q)) def |= q;
			if (t.isSparse(q)) sparse |= q;
		}
		p.indices = m->CreateTag(t.GetTagName() + "_index", DATA_INTEGER, def, sparse, t.GetSize());
		INMOST_DATA_ENUM_TYPE ret = static_cast<INMOST_DATA_ENUM_TYPE>(reg_tags.size());
		reg_tags.push_back(p);
		index_tags.push_back(p);
		return ret;
	}
	void Automatizator::EnumerateDynamicTags()
	{
		first_num = last_num = 0;
		const ElementType paralleltypes = NODE | EDGE | FACE | CELL;
		for (index_enum::iterator it = index_tags.begin(); it != index_tags.end(); ++it)
		{
			for (ElementType etype = NODE; etype <= MESH; etype = etype << 1)
				if (it->indices.isDefined(etype) && it->indices.isSparse(etype))
				{
					for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
						jt->DelData(it->indices);
				}
		}


		for (index_enum::iterator it = index_tags.begin(); it != index_tags.end(); ++it)
		{
			for (ElementType etype = NODE; etype <= MESH; etype = etype << 1)
			{
				if (it->indices.isDefined(etype))
				{
					if (it->indices.GetSize() == ENUMUNDEF)
					{
						if (!it->indices.isSparse(etype))
						{
							for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
							{
								if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask)))
								{
									Storage::integer_array indarr = jt->IntegerArray(it->indices);
									indarr.resize(jt->RealArray(it->d.t).size());
									for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
										*qt = last_num++;
								}
							}
						}
						else
						{
							for (Mesh::iteratorStorage jt = m->Begin(etype); jt != m->End(); ++jt)
							{
								if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (jt->HaveData(it->d.t) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask))))
								{
									Storage::integer_array indarr = jt->IntegerArray(it->indices);
									indarr.resize(jt->RealArray(it->d.t).size());
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
								if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask)))
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
								if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (jt->HaveData(it->d.t) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask))))
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
#if defined(USE_MPI)
		if (m->GetProcessorsNumber() > 1)
		{
			MPI_Scan(&last_num, &first_num, 1, INMOST_MPI_DATA_ENUM_TYPE, MPI_SUM, m->GetCommunicator());
			first_num -= last_num;
			ElementType exch_mask = NONE;
			for (index_enum::iterator it = index_tags.begin(); it != index_tags.end(); ++it)
			{
				for (ElementType etype = NODE; etype <= MESH; etype = etype << 1)
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
									if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask)))
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
									if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (jt->HaveData(it->d.t) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask))))
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										indarr.resize(jt->RealArray(it->d.t).size());
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
									if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask)))
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
									if ((!(etype & paralleltypes) || ((etype & paralleltypes) && jt->getAsElement()->GetStatus() != Element::Ghost)) && (jt->HaveData(it->d.t) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask))))
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
			last_num += first_num;
			{
				std::vector<Tag> exch_tags;
				for (index_enum::iterator it = index_tags.begin(); it != index_tags.end(); ++it) exch_tags.push_back(it->indices);
				m->ExchangeData(exch_tags, exch_mask,0);
			}
		}
#endif
		// this version will fail in parallel
		//merger.Resize(first_num,last_num,false);
		// use this version until there is a way to define multiple intervals in RowMerger
		INMOST_DATA_INTEGER_TYPE max_unknowns = m->AggregateMax(static_cast<INMOST_DATA_INTEGER_TYPE>(last_num));
#if defined(USE_OMP)
#pragma omp parallel
		{
#pragma omp single
			{
				merger.resize(omp_get_num_procs());
			}
			merger[omp_get_thread_num()].Resize(0,max_unknowns,false);
		}
#else
		merger.Resize(0,max_unknowns,false);
#endif
	}
};

#endif //USE_AUTODIFF
