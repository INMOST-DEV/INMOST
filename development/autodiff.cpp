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
#if !defined(EXPRESSION_TEMPLATES)
	INMOST_DATA_REAL_TYPE Automatizator::DerivativePrecompute(const expr & var, Storage * e, precomp_values_t & values, void * user_data)
	{
		assert(var.op != AD_NONE);
		INMOST_DATA_REAL_TYPE lval, rval, ret = 0.0;
		switch (var.op)
		{
		case AD_COND:
			lval = Evaluate(*var.left, e, user_data);
			rval = DerivativePrecompute(*(lval > 0.0 ? var.right->left : var.right->right), e, values, user_data);
			values.push_back(lval);
			return rval*var.coef;
		case AD_PLUS:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			rval = DerivativePrecompute(*var.right, e, values, user_data);
			return (lval + rval)*var.coef;
		case AD_MINUS:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			rval = DerivativePrecompute(*var.right, e, values, user_data);
			return (lval - rval)*var.coef;
		case AD_MULT:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			rval = DerivativePrecompute(*var.right, e, values, user_data);
			values.push_back(lval);
			values.push_back(rval);
			return (lval * rval)*var.coef;
		case AD_DIV:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			rval = DerivativePrecompute(*var.right, e, values, user_data);
			values.push_back(lval);
			values.push_back(rval);
			return (lval / rval)*var.coef;
		case AD_POW:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			rval = DerivativePrecompute(*var.right, e, values, user_data);
			values.push_back(lval);
			values.push_back(rval);
			ret = ::pow(lval, rval);
			values.push_back(ret);
			return ret*var.coef;
		case AD_INV:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			return var.coef / lval;
		case AD_ABS:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			return ::fabs(lval)*var.coef;
		case AD_EXP:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			ret = ::exp(lval);
			values.push_back(ret);
			return ret*var.coef;
		case AD_LOG:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			return ::log(lval)*var.coef;
		case AD_SIN:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			return ::sin(lval)*var.coef;
		case AD_COS:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			return ::cos(lval)*var.coef;
		case AD_CONST:
			return var.coef;
		case AD_MES:
			assert(!(e->GetElementType() & (ESET | MESH)));
			m->GetGeometricData(static_cast<Element *>(e), MEASURE, &ret);
			return ret*var.coef;
		}
		if (var.op >= AD_FUNC)
		{
			ret = reg_funcs[var.op].func(e, user_data);
			return ret*var.coef;
		}
		else if (var.op >= AD_TABLE)
		{
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			ret = reg_tables[var.op]->get_value(lval);
			return ret*var.coef;
		}
		else if (var.op >= AD_STNCL)
		{
			stencil_kind_domain st = reg_stencils[var.op];
			assert(st.domainmask == 0 || e->GetMarker(st.domainmask));
			if (st.kind == 0)
			{
				Storage::reference_array elems = e->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
				Storage::real_array coefs = e->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
				assert(elems.size() == coefs.size());
				for (INMOST_DATA_ENUM_TYPE k = 0; k < elems.size(); ++k)
				{
					lval = DerivativePrecompute(*var.left, elems[k], values, user_data);
					ret += lval * coefs[k];
				}
			}
			else if (st.kind == 1)
			{
				stencil_pairs get_st;
				reinterpret_cast<stencil_callback>(st.link)(e, get_st, user_data);
				for (INMOST_DATA_ENUM_TYPE k = 0; k < get_st.size(); ++k)
				{
					lval = DerivativePrecompute(*var.left, get_st[k].first, values, user_data);
					ret += lval * get_st[k].second;
				}
			}
			return ret*var.coef;
		}
		else if (var.op >= AD_CTAG)
		{
			INMOST_DATA_ENUM_TYPE comp = *(INMOST_DATA_ENUM_TYPE *)(&var.left);
			ret = GetStaticValue(e, var.op, comp);
			return ret * var.coef;
		}
		else if (var.op >= AD_TAG)
		{
			INMOST_DATA_ENUM_TYPE comp = *(INMOST_DATA_ENUM_TYPE *)(&var.left);
			ret = GetDynamicValue(e, var.op, comp);
			return ret*var.coef;
		}
		assert(false);
		return 0.0;
	}
	//! returns offset from the end of precomputed values
	void Automatizator::DerivativeFill(const expr & var, Storage * e, Solver::Row & entries, precomp_values_t & values, INMOST_DATA_REAL_TYPE multval, void * user_data)
	{
		assert(var.op != AD_NONE);
		INMOST_DATA_REAL_TYPE lval, rval, ret;
		switch (var.op)
		{
		case AD_COND:
			lval = values.back(); values.pop_back();
			DerivativeFill(*(lval > 0.0 ? var.right->left : var.right->right), e, entries, values, multval*var.coef, user_data);
			return;
		case AD_PLUS:
			DerivativeFill(*var.right, e, entries, values, multval*var.coef, user_data);
			DerivativeFill(*var.left, e, entries, values, multval*var.coef, user_data);
			return;
		case AD_MINUS:
			DerivativeFill(*var.right, e, entries, values, -multval*var.coef, user_data);
			DerivativeFill(*var.left, e, entries, values, multval*var.coef, user_data);
			return;
		case AD_MULT:
			rval = values.back(); values.pop_back();
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.right, e, entries, values, lval*multval*var.coef, user_data);
			DerivativeFill(*var.left, e, entries, values, rval*multval*var.coef, user_data);
			return;
		case AD_DIV:
			rval = values.back(); values.pop_back();
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.right, e, entries, values, -multval * lval / (rval*rval) * var.coef, user_data);
			DerivativeFill(*var.left, e, entries, values, multval / rval*var.coef, user_data);
			return;
		case AD_POW:
			ret = values.back(); values.pop_back();
			rval = values.back(); values.pop_back();
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.right, e, entries, values, multval * ret * ::log(lval) * var.coef, user_data);
			DerivativeFill(*var.left, e, entries, values, multval * ret * rval / lval * var.coef, user_data);
			return;
		case AD_INV:
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, -multval / (lval * lval) * var.coef, user_data);
			return;
		case AD_ABS:
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, multval * var.coef * (lval > 0 ? 1.0 : -1.0), user_data);
			return;
		case AD_EXP:
			ret = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, multval * ret * var.coef, user_data);
			return;
		case AD_LOG:
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, multval * var.coef / lval, user_data);
			return;
		case AD_SIN:
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, multval * var.coef * ::cos(lval), user_data);
			return;
		case AD_COS:
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, -multval * var.coef * ::sin(lval), user_data);
			return;
		case AD_CONST: return;
		case AD_MES: return;
		}
		if (var.op >= AD_FUNC)
		{
			return;
		}
		else if (var.op >= AD_TABLE)
		{
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, multval * var.coef * reg_tables[var.op]->get_derivative(lval), user_data);
			return;
		}
		else if (var.op >= AD_STNCL)
		{
			stencil_kind_domain st = reg_stencils[var.op];
			assert(st.domainmask == 0 || e->GetMarker(st.domainmask));
			if (st.kind == 0)
			{
				Storage::reference_array elems = e->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
				Storage::real_array coefs = e->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
				assert(elems.size() == coefs.size());
				for (INMOST_DATA_ENUM_TYPE k = elems.size(); k > 0; --k)
					DerivativeFill(*var.left, elems[k - 1], entries, values, var.coef * coefs[k - 1] * multval, user_data);
			}
			else if (st.kind == 1)
			{
				stencil_pairs get_st;
				reinterpret_cast<stencil_callback>(st.link)(e, get_st, user_data);
				for (INMOST_DATA_ENUM_TYPE k = get_st.size(); k > 0; --k)
					DerivativeFill(*var.left, get_st[k - 1].first, entries, values, var.coef * get_st[k - 1].second*multval, user_data);
			}
			return;
		}
		else if (var.op >= AD_CTAG) return;
		else if (var.op >= AD_TAG)
		{
			if (isDynamicValid(e, var.op))
			{
				INMOST_DATA_ENUM_TYPE ind = GetDynamicIndex(e, var.op, *(INMOST_DATA_ENUM_TYPE *)(&var.left));
				entries[ind] += multval * var.coef;
			}
			return;
		}
		assert(false);
		return;
	}
	INMOST_DATA_REAL_TYPE Automatizator::Evaluate(const expr & var, Storage * e, void * user_data)
	{
		assert(var.op != AD_NONE);
		switch (var.op)
		{
		case AD_COND:  return Evaluate(*(Evaluate(*var.left, e, user_data) > 0.0 ? var.right->left : var.right->right), e, user_data)*var.coef;
		case AD_PLUS:  return (Evaluate(*var.left, e, user_data) + Evaluate(*var.right, e, user_data))*var.coef;
		case AD_MINUS: return (Evaluate(*var.left, e, user_data) - Evaluate(*var.right, e, user_data))*var.coef;
		case AD_MULT:  return (Evaluate(*var.left, e, user_data) * Evaluate(*var.right, e, user_data))*var.coef;
		case AD_DIV:   return (Evaluate(*var.left, e, user_data) / Evaluate(*var.right, e, user_data))*var.coef;
		case AD_INV:   return var.coef / Evaluate(*var.left, e, user_data);
		case AD_POW:   return ::pow(Evaluate(*var.left, e, user_data), Evaluate(*var.right, e, user_data))*var.coef;
		case AD_ABS:   return ::fabs(Evaluate(*var.left, e, user_data))*var.coef;
		case AD_EXP:   return ::exp(Evaluate(*var.left, e, user_data))*var.coef;
		case AD_LOG:   return ::log(Evaluate(*var.left, e, user_data))*var.coef;
		case AD_SIN:   return ::sin(Evaluate(*var.left, e, user_data))*var.coef;
		case AD_COS:   return ::cos(Evaluate(*var.left, e, user_data))*var.coef;
		case AD_CONST: return var.coef;
		case AD_MES: assert(!(e->GetElementType() & (ESET | MESH))); Storage::real ret; m->GetGeometricData(static_cast<Element *>(e), MEASURE, &ret); return ret*var.coef;
		}
		if (var.op >= AD_FUNC) return reg_funcs[var.op].func(e, user_data);
		if (var.op >= AD_TABLE) return reg_tables[var.op]->get_value(Evaluate(*var.left, e, user_data))*var.coef;
		if (var.op >= AD_STNCL)
		{
			INMOST_DATA_REAL_TYPE ret = 0.0;
			stencil_kind_domain st = reg_stencils[var.op];
			assert(st.domainmask == 0 || e->GetMarker(st.domainmask));
			if (st.kind == 0)
			{
				Storage::reference_array elems = e->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
				Storage::real_array coefs = e->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
				assert(elems.size() == coefs.size());
				for (INMOST_DATA_ENUM_TYPE k = 0; k < elems.size(); ++k)
					ret += var.coef * Evaluate(*var.left, elems[k], user_data) * coefs[k];
			}
			else if (st.kind == 1)
			{
				stencil_pairs get_st;
				reinterpret_cast<stencil_callback>(st.link)(e, get_st, user_data);
				for (INMOST_DATA_ENUM_TYPE k = 0; k < get_st.size(); ++k)
					ret += var.coef * Evaluate(*var.left, get_st[k].first, user_data) * get_st[k].second;
			}
			return ret;
		}
		if (var.op >= AD_CTAG) return GetStaticValue(e, var.op, *(INMOST_DATA_ENUM_TYPE *)(&var.left))*var.coef;
		if (var.op >= AD_TAG) return  GetDynamicValue(e, var.op, *(INMOST_DATA_ENUM_TYPE *)(&var.left))*var.coef;
		assert(false);
		return 0.0;
	}

	INMOST_DATA_REAL_TYPE Automatizator::Derivative(const expr & var, Storage * e, Solver::Row & out, void * user_data)
	{
		INMOST_DATA_REAL_TYPE ret;
		precomp_values_t values;
		ret = DerivativePrecompute(var, e, values, user_data);
		DerivativeFill(var, e, out, values, 1.0, user_data);
		return ret;
	}
#endif
	Automatizator::Automatizator(Mesh * m) :first_num(0), last_num(0), m(m) {}
	Automatizator::~Automatizator()
	{
		for (unsigned k = 0; k < index_tags.size(); k++)
			index_tags[k].indices = m->DeleteTag(index_tags[k].indices);
		for (table_type::iterator it = reg_tables.begin(); it != reg_tables.end(); ++it)
		{
			delete[] it->second->args;
			delete[] it->second->vals;
			delete it->second;
		}
		for (stencil_type::iterator it = reg_stencils.begin(); it != reg_stencils.end(); ++it)
		if (it->second.kind == 0)
			delete static_cast<stencil_tag *>(it->second.link);
	}
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterFunc(std::string name, func_callback func)
	{
		INMOST_DATA_ENUM_TYPE ret = reg_funcs.size() + AD_FUNC;
		func_name_callback v;
		v.name = name;
		v.func = func;
		reg_funcs[ret] = v;
		return ret;
	}
	//register stencil that can be got from tags
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterStencil(std::string name, Tag elements_tag, Tag coefs_tag, MIDType domain_mask)
	{
		INMOST_DATA_ENUM_TYPE ret = reg_stencils.size() + AD_STNCL;
		stencil_kind_domain st;
		stencil_tag * save = new stencil_tag;
		st.name = name;
		save->coefs = coefs_tag;
		save->elements = elements_tag;
		st.kind = 0;
		st.link = static_cast<void *>(save);
		st.domainmask = domain_mask;
		reg_stencils[ret] = st;
		return ret;
	}
	//register stencil that can be got from function
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterStencil(std::string name, stencil_callback func, MIDType domain_mask)
	{
		INMOST_DATA_ENUM_TYPE ret = reg_stencils.size() + AD_STNCL;
		stencil_kind_domain st;
		st.name = name;
		st.kind = 1;
		st.link = reinterpret_cast<void *>(func);
		st.domainmask = domain_mask;
		reg_stencils[ret] = st;
		return ret;
	}
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterTable(std::string name, INMOST_DATA_REAL_TYPE * Arguments, INMOST_DATA_REAL_TYPE * Values, INMOST_DATA_ENUM_TYPE size)
	{
		INMOST_DATA_ENUM_TYPE ret = reg_tables.size() + AD_TABLE;
		table_ptr t = new table;
		t->name = name;
		t->args = new INMOST_DATA_REAL_TYPE[size];
		memcpy(t->args, Arguments, sizeof(INMOST_DATA_REAL_TYPE)*size);
		t->vals = new INMOST_DATA_REAL_TYPE[size];
		memcpy(t->vals, Values, sizeof(INMOST_DATA_REAL_TYPE)*size);
		t->size = size;
		reg_tables[ret] = t;
		return ret;
	}
	/// set data of tag t defined on domain_mask to be dynamic data
	/// don't register tag twice
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterDynamicTag(Tag t, ElementType typemask, MIDType domain_mask)
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
		INMOST_DATA_ENUM_TYPE ret = reg_tags.size() + AD_TAG;
		reg_tags[ret] = p;
		index_tags.push_back(p);
		return ret;
	}
	/// set index for every data entry of dynamic tag
	void Automatizator::EnumerateDynamicTags()
	{
		first_num = last_num = 0;
		const ElementType paralleltypes = NODE | EDGE | FACE | CELL;
		for (index_enum::iterator it = index_tags.begin(); it != index_tags.end(); ++it)
		for (ElementType etype = NODE; etype <= MESH; etype = etype << 1)
		if (it->indices.isDefined(etype) && it->indices.isSparse(etype))
		for (Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
			jt->DelData(it->indices);


		for (index_enum::iterator it = index_tags.begin(); it != index_tags.end(); ++it)
		{
			for (ElementType etype = NODE; etype <= MESH; etype = etype << 1)
			if (it->indices.isDefined(etype))
			{
				if (it->indices.GetSize() == ENUMUNDEF)
				{
					if (!it->indices.isSparse(etype))
					{
						for (Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
						if (((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask)))
						{
							Storage::integer_array indarr = jt->IntegerArray(it->indices);
							indarr.resize(jt->RealArray(it->d.t).size());
							for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
								*qt = last_num++;
						}
					}
					else
					{
						for (Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
						if (((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (jt->HaveData(it->d.t) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask))))
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
					if (!it->indices.isSparse(etype))
					{
						for (Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
						if (((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask)))
						{
							Storage::integer_array indarr = jt->IntegerArray(it->indices);
							for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
								*qt = last_num++;
						}
					}
					else
					{
						for (Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
						if (((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (jt->HaveData(it->d.t) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask))))
						{
							Storage::integer_array indarr = jt->IntegerArray(it->indices);
							for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
								*qt = last_num++;
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
			if (first_num > 0) for (index_enum::iterator it = index_tags.begin(); it != index_tags.end(); ++it)
			{
				for (ElementType etype = NODE; etype <= MESH; etype = etype << 1)
				if (it->indices.isDefined(etype))
				{
					exch_mask |= etype;
					if (it->indices.GetSize() == ENUMUNDEF)
					{
						if (!it->indices.isSparse(etype))
						{
							for (Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
							if (((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask)))
							{
								Storage::integer_array indarr = jt->IntegerArray(it->indices);
								for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
									*qt += first_num;
							}
						}
						else
						{
							for (Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
							if (((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (jt->HaveData(it->d.t) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask))))
							{
								Storage::integer_array indarr = jt->IntegerArray(it->indices);
								indarr.resize(jt->RealArray(it->d.t).size());
								for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
									*qt += first_num;
							}
						}
					}
					else
					{
						if (!it->indices.isSparse(etype))
						{
							for (Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
							if (((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask)))
							{
								Storage::integer_array indarr = jt->IntegerArray(it->indices);
								for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
									*qt += first_num;
							}
						}
						else
						{
							for (Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
							if (((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (jt->HaveData(it->d.t) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask))))
							{
								Storage::integer_array indarr = jt->IntegerArray(it->indices);
								for (Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
									*qt += first_num;
							}
						}
					}
				}
			}
			last_num += first_num;
			{
				std::vector<Tag> exch_tags;
				for (index_enum::iterator it = index_tags.begin(); it != index_tags.end(); ++it) exch_tags.push_back(it->indices);
				m->ExchangeData(exch_tags, exch_mask);
			}
		}
#endif

	}
	/// register tag, data for which don't change through iterations
	/// don't register tag twice
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterStaticTag(Tag t, MIDType domain_mask)
	{
		INMOST_DATA_ENUM_TYPE ret = reg_ctags.size() + AD_CTAG;
		tagdomain d;
		d.t = t;
		d.domain_mask = domain_mask;
		reg_ctags[ret] = d;
		return ret;
	}
#if defined(USE_AUTODIFF_ASMJIT)
	/*
	INMOST_DATA_REAL_TYPE Automatizator::DerivativePrecompute(const expr & var, Storage * e, precomp_values_t & values, void * user_data)
	{
		assert(var.op != AD_NONE);
		INMOST_DATA_REAL_TYPE lval, rval, ret = 0.0;
		switch (var.op)
		{
		case AD_COND:
			lval = Evaluate(*var.left, e, user_data);
			rval = DerivativePrecompute(*(lval > 0.0 ? var.right->left : var.right->right), e, values, user_data);
			values.push_back(lval);
			return rval*var.coef;
		case AD_PLUS:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			rval = DerivativePrecompute(*var.right, e, values, user_data);
			return (lval + rval)*var.coef;
		case AD_MINUS:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			rval = DerivativePrecompute(*var.right, e, values, user_data);
			return (lval - rval)*var.coef;
		case AD_MULT:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			rval = DerivativePrecompute(*var.right, e, values, user_data);
			values.push_back(lval);
			values.push_back(rval);
			return (lval * rval)*var.coef;
		case AD_DIV:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			rval = DerivativePrecompute(*var.right, e, values, user_data);
			values.push_back(lval);
			values.push_back(rval);
			return (lval / rval)*var.coef;
		case AD_POW:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			rval = DerivativePrecompute(*var.right, e, values, user_data);
			values.push_back(lval);
			values.push_back(rval);
			ret = ::pow(lval, rval);
			values.push_back(ret);
			return ret*var.coef;
		case AD_INV:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			return var.coef / lval;
		case AD_ABS:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			return ::fabs(lval)*var.coef;
		case AD_EXP:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			ret = ::exp(lval);
			values.push_back(ret);
			return ret*var.coef;
		case AD_LOG:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			return ::log(lval)*var.coef;
		case AD_SIN:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			return ::sin(lval)*var.coef;
		case AD_COS:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			return ::cos(lval)*var.coef;
		case AD_CONST:
			return var.coef;
		case AD_MES:
			assert(!(e->GetElementType() & (ESET | MESH)));
			m->GetGeometricData(static_cast<Element *>(e), MEASURE, &ret);
			return ret*var.coef;
		}
		if (var.op >= AD_FUNC)
		{
			ret = reg_funcs[var.op].func(e, user_data);
			return ret*var.coef;
		}
		else if (var.op >= AD_TABLE)
		{
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			ret = reg_tables[var.op]->get_value(lval);
			return ret*var.coef;
		}
		else if (var.op >= AD_STNCL)
		{
			stencil_kind_domain st = reg_stencils[var.op];
			assert(st.domainmask == 0 || e->GetMarker(st.domainmask));
			if (st.kind == 0)
			{
				Storage::reference_array elems = e->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
				Storage::real_array coefs = e->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
				assert(elems.size() == coefs.size());
				for (INMOST_DATA_ENUM_TYPE k = 0; k < elems.size(); ++k)
				{
					lval = DerivativePrecompute(*var.left, elems[k], values, user_data);
					ret += lval * coefs[k];
				}
			}
			else if (st.kind == 1)
			{
				stencil_pairs get_st;
				reinterpret_cast<stencil_callback>(st.link)(e, get_st, user_data);
				for (INMOST_DATA_ENUM_TYPE k = 0; k < get_st.size(); ++k)
				{
					lval = DerivativePrecompute(*var.left, get_st[k].first, values, user_data);
					ret += lval * get_st[k].second;
				}
			}
			return ret*var.coef;
		}
		else if (var.op >= AD_CTAG)
		{
			INMOST_DATA_ENUM_TYPE comp = *(INMOST_DATA_ENUM_TYPE *)(&var.left);
			ret = GetStaticValue(e, var.op, comp);
			return ret * var.coef;
		}
		else if (var.op >= AD_TAG)
		{
			INMOST_DATA_ENUM_TYPE comp = *(INMOST_DATA_ENUM_TYPE *)(&var.left);
			ret = GetDynamicValue(e, var.op, comp);
			return ret*var.coef;
		}
		assert(false);
		return 0.0;
	}
	//! returns offset from the end of precomputed values
	void Automatizator::DerivativeFill(const expr & var, Storage * e, Solver::Row & entries, precomp_values_t & values, INMOST_DATA_REAL_TYPE multval, void * user_data)
	{
		assert(var.op != AD_NONE);
		INMOST_DATA_REAL_TYPE lval, rval, ret;
		switch (var.op)
		{
		case AD_COND:
			lval = values.back(); values.pop_back();
			DerivativeFill(*(lval > 0.0 ? var.right->left : var.right->right), e, entries, values, multval*var.coef, user_data);
			return;
		case AD_PLUS:
			DerivativeFill(*var.right, e, entries, values, multval*var.coef, user_data);
			DerivativeFill(*var.left, e, entries, values, multval*var.coef, user_data);
			return;
		case AD_MINUS:
			DerivativeFill(*var.right, e, entries, values, -multval*var.coef, user_data);
			DerivativeFill(*var.left, e, entries, values, multval*var.coef, user_data);
			return;
		case AD_MULT:
			rval = values.back(); values.pop_back();
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.right, e, entries, values, lval*multval*var.coef, user_data);
			DerivativeFill(*var.left, e, entries, values, rval*multval*var.coef, user_data);
			return;
		case AD_DIV:
			rval = values.back(); values.pop_back();
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.right, e, entries, values, -multval * lval / (rval*rval) * var.coef, user_data);
			DerivativeFill(*var.left, e, entries, values, multval / rval*var.coef, user_data);
			return;
		case AD_POW:
			ret = values.back(); values.pop_back();
			rval = values.back(); values.pop_back();
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.right, e, entries, values, multval * ret * ::log(lval) * var.coef, user_data);
			DerivativeFill(*var.left, e, entries, values, multval * ret * rval / lval * var.coef, user_data);
			return;
		case AD_INV:
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, -multval / (lval * lval) * var.coef, user_data);
			return;
		case AD_ABS:
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, multval * var.coef * (lval > 0 ? 1.0 : -1.0), user_data);
			return;
		case AD_EXP:
			ret = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, multval * ret * var.coef, user_data);
			return;
		case AD_LOG:
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, multval * var.coef / lval, user_data);
			return;
		case AD_SIN:
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, multval * var.coef * ::cos(lval), user_data);
			return;
		case AD_COS:
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, -multval * var.coef * ::sin(lval), user_data);
			return;
		case AD_CONST: return;
		case AD_MES: return;
		}
		if (var.op >= AD_FUNC)
		{
			return;
		}
		else if (var.op >= AD_TABLE)
		{
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, multval * var.coef * reg_tables[var.op]->get_derivative(lval), user_data);
			return;
		}
		else if (var.op >= AD_STNCL)
		{
			stencil_kind_domain st = reg_stencils[var.op];
			assert(st.domainmask == 0 || e->GetMarker(st.domainmask));
			if (st.kind == 0)
			{
				Storage::reference_array elems = e->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
				Storage::real_array coefs = e->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
				assert(elems.size() == coefs.size());
				for (INMOST_DATA_ENUM_TYPE k = elems.size(); k > 0; --k)
					DerivativeFill(*var.left, elems[k - 1], entries, values, var.coef * coefs[k - 1] * multval, user_data);
			}
			else if (st.kind == 1)
			{
				stencil_pairs get_st;
				reinterpret_cast<stencil_callback>(st.link)(e, get_st, user_data);
				for (INMOST_DATA_ENUM_TYPE k = get_st.size(); k > 0; --k)
					DerivativeFill(*var.left, get_st[k - 1].first, entries, values, var.coef * get_st[k - 1].second*multval, user_data);
			}
			return;
		}
		else if (var.op >= AD_CTAG) return;
		else if (var.op >= AD_TAG)
		{
			if (isDynamicValid(e, var.op))
			{
				INMOST_DATA_ENUM_TYPE ind = GetDynamicIndex(e, var.op, *(INMOST_DATA_ENUM_TYPE *)(&var.left));
				entries[ind] += multval * var.coef;
			}
			return;
		}
		assert(false);
		return;
	}
	*/

	asmjit::kVarType GetEnumType()
	{
		switch (sizeof(INMOST_DATA_ENUM_TYPE))
		{
		case 1: return asmjit::kVarTypeUInt8;
		case 2: return asmjit::kVarTypeUInt16;
		case 4: return asmjit::kVarTypeUInt32;
		case 8: return asmjit::kVarTypeUInt64;
		}
		assert(0);
		return asmjit::kVarTypeInvalid;
	}
	asmjit::kVarType GetIntegerType()
	{
		switch (sizeof(INMOST_DATA_INTEGER_TYPE))
		{
		case 1: return asmjit::kVarTypeInt8;
		case 2: return asmjit::kVarTypeInt16;
		case 4: return asmjit::kVarTypeInt32;
		case 8: return asmjit::kVarTypeInt64;
		}
		assert(0);
		return asmjit::kVarTypeInvalid;
	}
	asmjit::kVarType GetRealType()
	{
		switch (sizeof(INMOST_DATA_REAL_TYPE))
		{
		case 4: return asmjit::kVarTypeFp32;
		case 8: return asmjit::kVarTypeFp64;
		case 10:
		case 16:
			return asmjit::kVarTypeFpEx;
		}
		assert(0);
		return asmjit::kVarTypeInvalid;
	}

	void * Automatizator::GenerateEvaluateASMJIT(const expr & var)
	{
		using namespace asmjit;
		using namespace asmjit::host;
		Compiler comp(&runtime);
		X86X64FuncNode * func = comp.addFunc(kFuncConvHostCDecl, FuncBuilder2<INMOST_DATA_REAL_TYPE, Storage *, void *>());
		GpVar a0(comp, kVarTypeUIntPtr);
		GpVar a1(comp, kVarTypeUIntPtr);
		GpVar ar(comp, GetRealType());
		comp.setArg(0, a0);
		comp.setArg(1, a1);
		XmmVar ret = GenerateEvaluateSubASMJIT(var, comp,a0,a1);
		comp.movq(ar, ret);
		comp.ret(ar);
		comp.endFunc();
		return comp.make();
		
	}

	
	void MultCoef(INMOST_DATA_REAL_TYPE coef, asmjit::host::Compiler & comp, asmjit::host::XmmVar & var)
	{
		using namespace asmjit;
		using namespace asmjit::host;
		if (fabs(coef - 1.0) > 1e-9)
		{
			XmmVar vconst(comp.newXmmVar());
			GpVar gpreg(comp,GetRealType());
			int64_t * i = reinterpret_cast<int64_t *>(&coef);
			comp.mov(gpreg, i[0]);
			if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
				comp.movd(vconst, gpreg);
			else
				comp.movq(vconst, gpreg);
			comp.unuse(gpreg);
			if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(var, vconst);
			else comp.mulsd(var, vconst);
			comp.unuse(vconst);
		}
	}
	
	asmjit::host::XmmVar Automatizator::GenerateEvaluateSubASMJIT(const expr & var, asmjit::host::Compiler & comp, asmjit::host::GpVar & elem, asmjit::host::GpVar & user_data)
	{
		using namespace asmjit;
		using namespace asmjit::host;
		assert(var.op != AD_NONE);
		
		switch (var.op)
		{
		case AD_COND:
		{
						XmmVar ret(comp.newXmmVar()), v0 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data);
						Label L_1(comp.newLabel()), L_2(comp.newLabel()), L_3(comp.newLabel());
						XmmVar zero(comp.newXmmVar());
						GpVar gpreg(comp.newGpVar());
						INMOST_DATA_REAL_TYPE zeroval = 0.0;
						uint64_t * i = reinterpret_cast<uint64_t *>(&zeroval);
						comp.mov(gpreg, i[0]);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						{
							comp.movd(zero, gpreg);
							comp.comiss(zero, v0);
						}
						else
						{
							comp.movq(zero, gpreg);
							comp.comisd(zero, v0);
						}
						comp.unuse(gpreg);
						comp.unuse(zero);
						comp.jae(L_1); // v0 >= 0.0
						comp.jb(L_2); // v0 < 0.0
						comp.bind(L_1);
						XmmVar v1 = GenerateEvaluateSubASMJIT(*var.right->left, comp, elem, user_data);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(ret, v1); else comp.movsd(ret, v1);
						comp.jmp(L_3);
						comp.bind(L_2);
						XmmVar v2 = GenerateEvaluateSubASMJIT(*var.right->right, comp, elem, user_data);
						comp.bind(L_3);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(ret, v2); else comp.movsd(ret, v2);
						MultCoef(var.coef, comp, ret);
						comp.unuse(v1);
						comp.unuse(v2);
						return ret;
		}
		case AD_PLUS:
		{
						XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data), v2 = GenerateEvaluateSubASMJIT(*var.right, comp, elem, user_data);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.addss(v1, v2); else comp.addsd(v1, v2);
						comp.unuse(v2);
						MultCoef(var.coef, comp, v1);
						return v1;
		}
		case AD_MINUS:
		{
						 XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data), v2 = GenerateEvaluateSubASMJIT(*var.right, comp, elem, user_data);
						 if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.subss(v1, v2); else comp.subsd(v1, v2);
						 comp.unuse(v2);
						 MultCoef(var.coef, comp, v1);
						 return v1;
		}
		case AD_MULT:
		{
						XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data), v2 = GenerateEvaluateSubASMJIT(*var.right, comp, elem, user_data);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(v1, v2); else comp.mulsd(v1, v2);
						comp.unuse(v2);
						MultCoef(var.coef, comp, v1);
						return v1;
		}
		case AD_DIV:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data), v2 = GenerateEvaluateSubASMJIT(*var.right, comp, elem, user_data);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.divss(v1, v2); else comp.divsd(v1, v2);
					   comp.unuse(v2);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_INV:
		{
					   
						XmmVar vconst(comp.newXmmVar());
						GpVar gpreg(comp.newGpVar());
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						{
							const uint32_t * i = reinterpret_cast<const uint32_t *>(&var.coef);
							comp.mov(gpreg, i[0]);
							comp.movd(vconst, gpreg);
						}
						else
						{
							const uint64_t * i = reinterpret_cast<const uint64_t *>(&var.coef);
							comp.mov(gpreg, i[0]);
							comp.movq(vconst, gpreg);
						}
						comp.unuse(gpreg);
						XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.divss(vconst, v1); else comp.divsd(vconst, v1);
					   comp.unuse(v1);
					   return vconst;
		}
		case AD_POW:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data), v2 = GenerateEvaluateSubASMJIT(*var.right, comp, elem, user_data);
					   GpVar address(comp.newGpVar());
					   X86X64CallNode *ctx;
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float, float)>(::pow))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double, double)>(::pow))));
					   ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder2<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setArg(0, v2);
					   ctx->setRet(0, v1);
					   comp.unuse(v2);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_ABS:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data);
					   GpVar address(comp.newGpVar());
					   X86X64CallNode *ctx;
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::fabs))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::fabs))));
					   ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_EXP:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data);
					   GpVar address(comp.newGpVar());
					   X86X64CallNode *ctx;
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::exp))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::exp))));
					   ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_LOG:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data);
					   GpVar address(comp.newGpVar());
					   X86X64CallNode *ctx;
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::log))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::log))));
					   ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_SIN:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data);
					   GpVar address(comp.newGpVar());
					   X86X64CallNode *ctx;
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::sin))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::sin))));
					   ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_COS:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data);
					   GpVar address(comp.newGpVar());
					   X86X64CallNode *ctx;
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::cos))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::cos))));
					   ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_CONST:
		{
						XmmVar vconst(comp.newXmmVar());
						GpVar gpreg(comp.newGpVar());
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						{
							const uint32_t * i = reinterpret_cast<const uint32_t *>(&var.coef);
							comp.mov(gpreg, i[0]);
							comp.movd(vconst, gpreg);
						}
						else
						{
							const uint64_t * i = reinterpret_cast<const uint64_t *>(&var.coef);
							comp.mov(gpreg, i[0]);
							comp.movq(vconst, gpreg);
						}
						comp.unuse(gpreg);
						 return vconst;
		}
		case AD_MES:
		{
					   XmmVar ret(comp.newXmmVar());
					   GpVar address(comp.newGpVar());
					   GpVar mesh(comp, kVarTypeUIntPtr);
					   GpVar fret(comp, GetRealType());
					   comp.mov(mesh, imm(reinterpret_cast<int64_t>(m)));
					   X86X64CallNode *ctx;
					   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Mesh *, Storage *)>(Automatizator::GetGeometricDataStatic))));
					   ctx = comp.call(address, kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, Mesh *, Storage *>());
					   ctx->setArg(0, mesh);
					   ctx->setArg(1, elem);
					   ctx->setRet(0, fret);
					   comp.movq(ret, fret);
					   comp.unuse(address);
					   comp.unuse(mesh);
					   comp.unuse(fret);
					   MultCoef(var.coef, comp, ret);
					   return ret;
		}
		}
		if (var.op >= AD_FUNC) 
		{
			XmmVar ret(comp.newXmmVar());
			GpVar address(comp.newGpVar());
			GpVar fret(comp, GetRealType());
			X86X64CallNode *ctx;
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Storage *, void *)>(reg_funcs[var.op].func))));
			ctx = comp.call(address, kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, Storage *, void *>());
			ctx->setArg(0, elem);
			ctx->setArg(1, user_data);
			ctx->setRet(0, fret);
			comp.movq(ret, fret);
			comp.unuse(address);
			comp.unuse(fret);
			MultCoef(var.coef, comp, ret);
			return ret;
		}
		if (var.op >= AD_TABLE)
		{
			XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data);
			GpVar address(comp.newGpVar());
			GpVar vtable(comp, kVarTypeUIntPtr);
			GpVar fin(comp, GetRealType()), fout(comp, GetRealType());
			comp.mov(vtable, imm_ptr(reinterpret_cast<void *>(reg_tables[var.op])));
			X86X64CallNode *ctx;
			comp.movq(fin, v1);
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(table_ptr,INMOST_DATA_REAL_TYPE)>(Automatizator::GetTableValueStatic))));
			ctx = comp.call(address, kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, table_ptr, INMOST_DATA_REAL_TYPE>());
			int64_t tableid = var.op;
			ctx->setArg(0, vtable);
			ctx->setArg(1, fin);
			ctx->setRet(0, fout);
			comp.movq(v1, fout);
			comp.unuse(address);
			comp.unuse(vtable);
			comp.unuse(fin);
			comp.unuse(fout);
			MultCoef(var.coef, comp, v1);
			return v1;
		}
		if (var.op >= AD_STNCL)
		{
			X86X64FuncNode * funcnew = comp.addFunc(kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, Storage *, void *>());
			GpVar a0(comp, kVarTypeUIntPtr);
			GpVar a1(comp, kVarTypeUIntPtr);
			GpVar fret(comp, GetRealType());
			comp.setArg(0, a0);
			comp.setArg(1, a1);
			XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, a0, a1);
			comp.movq(fret, v1);
			comp.ret(fret);
			comp.unuse(v1);
			comp.unuse(fret);
			comp.endFunc();
			Label L_1(comp.newLabel());
			XmmVar ret(comp.newXmmVar()), stencil_coef(comp.newXmmVar()), temp(comp.newXmmVar());

			GpVar address(comp.newGpVar()), address_coef(comp.newGpVar()), address_elem(comp.newGpVar());
			GpVar index(comp.newGpVar());
			GpVar stencil_size(comp, GetEnumType());
			GpVar stencil_elem(comp, kVarTypeUIntPtr);
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_ENUM_TYPE(*)(Automatizator *, INMOST_DATA_ENUM_TYPE, Storage *, void *)>(Automatizator::GetStencilStatic))));
			comp.mov(address_elem, imm(reinterpret_cast<int64_t>(static_cast<Storage *(*)(Automatizator *, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetStencilElemStatic))));
			comp.mov(address_coef, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE (*)(Automatizator *, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetStencilCoefStatic))));
			X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder4<INMOST_DATA_ENUM_TYPE, Automatizator *, INMOST_DATA_ENUM_TYPE, Storage *, void *>());
			GpVar self_ptr(comp, kVarTypeUIntPtr);
			GpVar stencil_id(comp, GetEnumType());
			comp.mov(stencil_id, imm(var.op));
			comp.mov(self_ptr, imm_ptr(this));
			ctx->setArg(0, self_ptr);
			ctx->setArg(1, stencil_id);
			ctx->setArg(2, elem);
			ctx->setArg(3, user_data);
			ctx->setRet(0, stencil_size);
			comp.mov(index, imm(0));
			comp.bind(L_1);

			ctx = comp.call(address_coef, kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, Automatizator *, INMOST_DATA_ENUM_TYPE>());
			ctx->setArg(0, imm(reinterpret_cast<int64_t>(this)));
			ctx->setArg(1, index);
			ctx->setRet(0, stencil_coef);

			ctx = comp.call(address_elem, kFuncConvHost, FuncBuilder2<Storage *, Automatizator *, INMOST_DATA_ENUM_TYPE>());
			ctx->setArg(0, imm(reinterpret_cast<int64_t>(this)));
			ctx->setArg(1, index);
			ctx->setRet(0, stencil_elem);

			ctx = comp.call(funcnew->getEntryLabel(), kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, Storage *, void *>());
			ctx->setArg(0, stencil_elem);
			ctx->setArg(1, user_data);
			ctx->setRet(0, temp);

			if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
			{
				comp.mulss(temp, stencil_coef);
				comp.addss(ret, temp);
			}
			else
			{
				comp.mulsd(temp, stencil_coef);
				comp.addsd(ret, temp);
			}



			comp.inc(index);
			comp.cmp(index, stencil_size);
			comp.jne(L_1);
			comp.unuse(address);
			comp.unuse(address_elem);
			comp.unuse(address_coef);
			comp.unuse(index);
			comp.unuse(stencil_size);
			comp.unuse(stencil_coef);
			comp.unuse(stencil_elem);
			comp.unuse(temp);

			MultCoef(var.coef, comp, ret);
			return ret;
		}
		if (var.op >= AD_CTAG)
		{
			XmmVar ret(comp.newXmmVar());
			GpVar address(comp.newGpVar());
			GpVar vtag(comp, kVarTypeUIntPtr);
			GpVar vcomp(comp, GetEnumType());
			GpVar fret(comp, GetRealType());
			comp.mov(vtag, imm_ptr(reinterpret_cast<void *>(&reg_ctags[var.op].t)));
			comp.mov(vcomp, imm(*(INMOST_DATA_ENUM_TYPE *)(&var.left)));
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Storage *, Tag *, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetValueStatic))));
			X86X64CallNode * ctx = comp.call(address, kFuncConvHost, FuncBuilder3<INMOST_DATA_REAL_TYPE, Storage *, Tag *, INMOST_DATA_ENUM_TYPE>());
			int64_t component = *(INMOST_DATA_ENUM_TYPE *)(&var.left);
			ctx->setArg(0, elem);
			ctx->setArg(1, vtag);
			ctx->setArg(2, vcomp);
			ctx->setRet(0,fret);
			comp.movq(ret, fret);
			comp.unuse(fret);
			comp.unuse(vtag);
			comp.unuse(vcomp);
			comp.unuse(address);
			MultCoef(var.coef, comp, ret);
			return ret;
		}
		if (var.op >= AD_TAG)
		{
			XmmVar ret(comp.newXmmVar());
			GpVar address(comp.newGpVar());
			GpVar vtag(comp, kVarTypeUIntPtr);
			GpVar vcomp(comp, GetEnumType());
			GpVar fret(comp, GetRealType());
			comp.mov(vtag, imm_ptr(reinterpret_cast<void *>(&reg_tags[var.op].d.t)));
			comp.mov(vcomp, imm(*(INMOST_DATA_ENUM_TYPE *)(&var.left)));
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Storage *, Tag *, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetValueStatic))));
			X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder3<INMOST_DATA_REAL_TYPE, Storage *, Tag *, INMOST_DATA_ENUM_TYPE>());
			ctx->setArg(0, elem);
			ctx->setArg(1, vtag);
			ctx->setArg(2, vcomp);
			ctx->setRet(0, fret);
			comp.movq(ret, fret);
			comp.unuse(fret);
			comp.unuse(vtag);
			comp.unuse(vcomp);
			comp.unuse(address);
			MultCoef(var.coef, comp, ret);
			return ret;
		}
		assert(false);
	}

#endif
};

#endif //USE_AUTODIFF
