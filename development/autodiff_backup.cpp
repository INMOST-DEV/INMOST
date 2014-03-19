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
				std::cout << ind << ", " << multval*var.coef << std::endl;
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
		case 4: return asmjit::kVarTypeInt32;//asmjit::kVarTypeFp32;
		case 8: return asmjit::kVarTypeInt64; //asmjit::kVarTypeFp64;
		case 10:
		case 16:
			return asmjit::kVarTypeFpEx;
		}
		assert(0);
		return asmjit::kVarTypeInvalid;
	}
	asmjit::host::kVarType GetRealTypeXmm()
	{
		switch (sizeof(INMOST_DATA_REAL_TYPE))
		{
		case 4: return asmjit::host::kVarTypeXmm;
		case 8: return asmjit::host::kVarTypeXmm;
		}
		assert(0);
		return asmjit::host::kVarTypeXmm;
	}
	void SetFloat(asmjit::host::Compiler & comp, asmjit::host::XmmVar & var, INMOST_DATA_REAL_TYPE value)
	{
		using namespace asmjit;
		using namespace asmjit::host;
		GpVar gpreg(comp.newGpVar(GetRealType()));
		const int64_t * i = reinterpret_cast<const int64_t *>(&value);
		comp.mov(gpreg, imm(i[0]));
		if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movd(var, gpreg); else comp.movq(var, gpreg);
		comp.unuse(gpreg);
	}
	void * Automatizator::GenerateEvaluateASMJIT(const expr & var)
	{
		using namespace asmjit;
		using namespace asmjit::host;
		Compiler comp(&runtime);
		X86X64FuncNode * func = comp.addFunc(kFuncConvHostCDecl, FuncBuilder2<INMOST_DATA_REAL_TYPE, Storage *, void *>());
		GpVar a0(comp.newGpVar(kVarTypeUIntPtr));
		GpVar a1(comp.newGpVar(kVarTypeUIntPtr));
		comp.setArg(0, a0);
		comp.setArg(1, a1);
		XmmVar ret = GenerateEvaluateSubASMJIT(var, comp,a0,a1, 0);
		comp.ret(ret);
		comp.endFunc();
		return comp.make();
	}

	void * Automatizator::GenerateDerivativeASMJIT(const expr & var)
	{
		using namespace asmjit;
		using namespace asmjit::host;
		Compiler comp(&runtime);
		X86X64FuncNode * func = comp.addFunc(kFuncConvHostCDecl, FuncBuilder3<INMOST_DATA_REAL_TYPE, Storage *, void *, Solver::Row *>());
		GpVar a0(comp.newGpVar(kVarTypeUIntPtr));
		GpVar a1(comp.newGpVar(kVarTypeUIntPtr));
		GpVar a2(comp.newGpVar(kVarTypeUIntPtr));
		comp.setArg(0, a0);
		comp.setArg(1, a1);
		comp.setArg(2, a2);
		XmmVar ret = GenerateDerivativePrecomputeSubASMJIT(var, comp, a0, a1, 0);
		XmmVar mult(comp.newXmmVar(GetRealTypeXmm()));
		SetFloat(comp, mult, 1.0);
		GenerateDerivativeFillSubASMJIT(var, comp, a0, a1, a2, mult, 0);
		comp.unuse(mult);
		comp.ret(ret);
		comp.endFunc();
		return comp.make();
	}

	void * Automatizator::GenerateCode(const expr & var)
	{
		typedef std::pair<void *, void *> ret_t;
		ret_t * code = new ret_t;
		code->first = GenerateEvaluateASMJIT(var);
		code->second = GenerateDerivativeASMJIT(var);
		return (void *)code;
	}

	
	
	void MultCoef(const INMOST_DATA_REAL_TYPE & coef, asmjit::host::Compiler & comp, asmjit::host::XmmVar & var)
	{
		using namespace asmjit;
		using namespace asmjit::host;
		if (fabs(coef - 1.0) > 1e-9)
		{
			XmmVar vconst(comp.newXmmVar(GetRealTypeXmm()));
			SetFloat(comp, vconst, coef);
			if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(var, vconst); else comp.mulsd(var, vconst);
			comp.unuse(vconst);
		}
	}

	INMOST_DATA_REAL_TYPE Automatizator::GetGeometricDataJIT(Mesh * m, Storage * e) 
	{ 
		Storage::real ret; 
		m->GetGeometricData(static_cast<Element *>(e), MEASURE, &ret);  
		return ret; 
	}
	INMOST_DATA_REAL_TYPE    Automatizator::GetTableValueJIT(INMOST_DATA_REAL_TYPE arg, table_ptr tab) { return tab->get_value(arg); }
	INMOST_DATA_REAL_TYPE    Automatizator::GetTableDerivativeJIT(INMOST_DATA_REAL_TYPE arg, table_ptr tab) { return tab->get_derivative(arg); }
	INMOST_DATA_REAL_TYPE    Automatizator::GetValueJIT(Storage * e, Tag * t, INMOST_DATA_ENUM_TYPE comp) {return e->RealArray(*t)[comp];}
	INMOST_DATA_ENUM_TYPE    Automatizator::GetIndexJIT(Storage * e, Tag * t, INMOST_DATA_ENUM_TYPE comp) { return static_cast<INMOST_DATA_ENUM_TYPE>(e->IntegerArray(*t)[comp]); }
	void                   Automatizator::PushJIT(INMOST_DATA_REAL_TYPE val, Automatizator * aut)
	{
		std::cout << "in: " << val << std::endl;
		aut->stacks[THREAD_NUM].push_back(val); 
		//return val;
	}
	INMOST_DATA_REAL_TYPE    Automatizator::PopJIT(Automatizator * aut)
	{ 
		INMOST_DATA_REAL_TYPE ret = aut->stacks[THREAD_NUM].back();
		std::cout << "out: " << ret << std::endl;
		aut->stacks[THREAD_NUM].pop_back();
		return ret;
	}
	INMOST_DATA_REAL_TYPE    Automatizator::InsertPairJIT(INMOST_DATA_REAL_TYPE val, INMOST_DATA_ENUM_TYPE ind, Solver::Row * row)
	{
		std::cout << ind << ", " << val << std::endl;
		(*row)[ind] += val;
		return val;
	}
	void                   Automatizator::ClearStencilJIT(Automatizator * aut, INMOST_DATA_ENUM_TYPE id, INMOST_DATA_ENUM_TYPE nested)
	{
		assert(nested < MAX_NESTED);
		aut->stencil_storage[THREAD_NUM][nested].pairs.clear();
	}
	INMOST_DATA_ENUM_TYPE    Automatizator::GetStencilJIT(Automatizator * aut, Storage * elem, void * user_data, INMOST_DATA_ENUM_TYPE id, INMOST_DATA_ENUM_TYPE nested)
	{
		assert(nested < MAX_NESTED);
		INMOST_DATA_ENUM_TYPE size = aut->GetStencil(id, elem, user_data, aut->stencil_storage[THREAD_NUM][nested].pairs);
		return size;
	}
	INMOST_DATA_REAL_TYPE    Automatizator::GetStencilCoefJIT(Automatizator * aut, INMOST_DATA_ENUM_TYPE k, INMOST_DATA_ENUM_TYPE nested) 
	{ 
		assert(nested < MAX_NESTED); 
		return aut->stencil_storage[THREAD_NUM][nested].pairs[k].second; 
	}
	Storage *              Automatizator::GetStencilElemJIT(Automatizator * aut, INMOST_DATA_ENUM_TYPE k, INMOST_DATA_ENUM_TYPE nested) 
	{ 
		assert(nested < MAX_NESTED); return aut->stencil_storage[THREAD_NUM][nested].pairs[k].first; 
	}
	
	asmjit::host::XmmVar   Automatizator::GenerateEvaluateSubASMJIT(const expr & var, asmjit::host::Compiler & comp, asmjit::host::GpVar & elem, asmjit::host::GpVar & user_data, INMOST_DATA_ENUM_TYPE nested)
	{
		using namespace asmjit;
		using namespace asmjit::host;
		assert(var.op != AD_NONE);
		
		switch (var.op)
		{
		case AD_COND:
		{
						XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
						XmmVar v0 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested);
						XmmVar zero(comp.newXmmVar(GetRealTypeXmm()));
						SetFloat(comp, zero, 0.0);
						Label Brench(comp.newLabel()), Exit(comp.newLabel());
						//GpVar vzero(comp.newGpVar(GetRealType()));
						//INMOST_DATA_REAL_TYPE zeroval = 0.0;
						//const int64_t * i = reinterpret_cast<const int64_t *>(&zeroval);
						//comp.mov(vzero, imm(i[0]));
						//comp.movq(zero, vzero);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.comiss(zero, v0); else comp.comisd(zero, v0); // zero -> zero <= v0
						//if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.cmpss(zero, v0, 2); else comp.cmpsd(zero, v0, 2); // zero -> zero <= v0
						//comp.movq(vzero, zero);
						//comp.cmp(vzero, imm(0));
						comp.unuse(v0);
						comp.ja(Brench);
						//comp.movq(ret, vzero); //redundant!!
						XmmVar v1 = GenerateEvaluateSubASMJIT(*var.right->left, comp, elem, user_data, nested);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(ret, v1); else comp.movsd(ret, v1);
						comp.unuse(v1);
						comp.jmp(Exit);
						comp.bind(Brench);
						//comp.movq(ret, vzero); //redundant!!
						XmmVar v2 = GenerateEvaluateSubASMJIT(*var.right->right, comp, elem, user_data, nested);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(ret, v2); else comp.movsd(ret, v2);
						comp.unuse(v2);
						comp.bind(Exit);
						MultCoef(var.coef, comp, ret);
						//comp.unuse(vzero);
						comp.unuse(zero);
						return ret;
		}
		case AD_PLUS:
		{
						XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested), v2 = GenerateEvaluateSubASMJIT(*var.right, comp, elem, user_data, nested);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.addss(v1, v2); else comp.addsd(v1, v2);
						comp.unuse(v2);
						MultCoef(var.coef, comp, v1);
						return v1;
		}
		case AD_MINUS:
		{
						 XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested), v2 = GenerateEvaluateSubASMJIT(*var.right, comp, elem, user_data, nested);
						 if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.subss(v1, v2); else comp.subsd(v1, v2);
						 comp.unuse(v2);
						 MultCoef(var.coef, comp, v1);	
						 return v1;
		}
		case AD_MULT:
		{
						XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested), v2 = GenerateEvaluateSubASMJIT(*var.right, comp, elem, user_data, nested);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(v1, v2); else comp.mulsd(v1, v2);
						comp.unuse(v2);
						MultCoef(var.coef, comp, v1);
						return v1;
		}
		case AD_DIV:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested), v2 = GenerateEvaluateSubASMJIT(*var.right, comp, elem, user_data, nested);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.divss(v1, v2); else comp.divsd(v1, v2);
					   comp.unuse(v2);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_INV:
		{
					   
					   XmmVar vconst(comp.newXmmVar(GetRealTypeXmm()));
						SetFloat(comp, vconst, var.coef);
						XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.divss(vconst, v1); else comp.divsd(vconst, v1);
					   comp.unuse(v1);
					   return vconst;
		}
		case AD_POW:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested), v2 = GenerateEvaluateSubASMJIT(*var.right, comp, elem, user_data, nested);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float, float)>(::pow))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double, double)>(::pow))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder2<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
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
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::fabs))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::fabs))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_EXP:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::exp))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::exp))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_LOG:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::log))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::log))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_SIN:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::sin))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::sin))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_COS:
		{
					   XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::cos))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::cos))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_CONST:
		{
						 XmmVar vconst(comp.newXmmVar(GetRealTypeXmm()));
						SetFloat(comp, vconst, var.coef);
						return vconst;
		}
		case AD_MES:
		{
					   XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
					   GpVar address(comp.newGpVar());
					   GpVar vmesh(comp.newGpVar(kVarTypeUIntPtr));
					   GpVar velem(comp.newGpVar(kVarTypeUIntPtr));
					   comp.mov(vmesh, imm(reinterpret_cast<int64_t>(m)));
					   comp.mov(velem, elem);
					   X86X64CallNode *ctx;
					   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Mesh *, Storage *)>(Automatizator::GetGeometricDataJIT))));
					   ctx = comp.call(address, kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, Mesh *, Storage *>());
					   ctx->setArg(0, vmesh);
					   ctx->setArg(1, velem);
					   ctx->setRet(0, ret);
					   comp.unuse(address);
					   comp.unuse(vmesh);
					   comp.unuse(velem);
					   MultCoef(var.coef, comp, ret);
					   return ret;
		}
		}
		if (var.op >= AD_FUNC) 
		{
			XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
			GpVar address(comp.newGpVar());
			GpVar velem(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vuser_data(comp.newGpVar(kVarTypeUIntPtr));
			comp.mov(velem, elem);
			comp.mov(vuser_data, user_data);
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Storage *, void *)>(reg_funcs[var.op].func))));
			X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, Storage *, void *>());
			ctx->setArg(0, velem);
			ctx->setArg(1, vuser_data);
			ctx->setRet(0, ret);
			comp.unuse(address);
			comp.unuse(velem);
			comp.unuse(vuser_data);
			MultCoef(var.coef, comp, ret);
			return ret;
		}
		if (var.op >= AD_TABLE)
		{
			XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested);
#if defined(DPRNT)
			{
				GpVar addrprnt(comp.newGpVar());
				comp.mov(addrprnt, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(INMOST_DATA_REAL_TYPE)>(Automatizator::PrintConst))));
				X86X64CallNode *ctx = comp.call(addrprnt, kFuncConvHost, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
				ctx->setArg(0, v1);
				ctx->setRet(0, v1);
			}
#endif
#if defined(DPRNT)
			{
				GpVar addrprnt(comp.newGpVar());
				GpVar opid(comp, GetEnumType());
				GpVar self(comp, kVarTypeUIntPtr);
				const int64_t vop = static_cast<const INMOST_DATA_ENUM_TYPE &>(var.op);
				comp.mov(opid, imm(vop));
				comp.mov(self, imm_ptr(this));
				comp.mov(addrprnt, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(INMOST_DATA_REAL_TYPE, Automatizator *, INMOST_DATA_ENUM_TYPE)>(Automatizator::PrintOp1))));
				X86X64CallNode *ctx = comp.call(addrprnt, kFuncConvHost, FuncBuilder3<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE, Automatizator *, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, v1);
				ctx->setArg(1, self);
				ctx->setArg(2, opid);
				ctx->setRet(0, v1);
			}
#endif

			GpVar address(comp.newGpVar());
			GpVar vtable(comp.newGpVar(kVarTypeUIntPtr));
			comp.mov(vtable, imm(reinterpret_cast<int64_t>(reg_tables[var.op])));
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(INMOST_DATA_REAL_TYPE,table_ptr)>(Automatizator::GetTableValueJIT))));
			X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE, table_ptr>());
			ctx->setArg(0, v1);
			ctx->setArg(1, vtable);
			ctx->setRet(0, v1);
			comp.unuse(address);
			comp.unuse(vtable);
			MultCoef(var.coef, comp, v1);
			return v1;
		}
		if (var.op >= AD_STNCL)
		{
			XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
			SetFloat(comp, ret, 0.0);
			
			GpVar index(comp.newGpVar(GetEnumType()));
			GpVar stencil_size(comp.newGpVar(GetEnumType()));
			GpVar stencil_elem(comp.newGpVar(kVarTypeUIntPtr));
			
			
			
			{
				const int64_t opid = static_cast<const INMOST_DATA_ENUM_TYPE &>(var.op);
				GpVar address(comp.newGpVar());
				GpVar self(comp.newGpVar(kVarTypeUIntPtr));
				GpVar velem(comp.newGpVar(kVarTypeUIntPtr));
				GpVar vuser_data(comp.newGpVar(kVarTypeUIntPtr));
				GpVar stencil_id(comp.newGpVar(GetEnumType()));
				GpVar vnested(comp.newGpVar(GetEnumType()));
				
				comp.mov(self, imm_ptr(this));
				comp.mov(velem, elem);
				comp.mov(vuser_data, user_data);
				comp.mov(stencil_id, imm(opid));
				comp.mov(vnested, imm(nested));
				comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_ENUM_TYPE(*)(Automatizator *, Storage *, void *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetStencilJIT))));
				X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder5<INMOST_DATA_ENUM_TYPE, Automatizator *, Storage *, void *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, self);
				ctx->setArg(1, velem);
				ctx->setArg(2, vuser_data);
				ctx->setArg(3, stencil_id);
				ctx->setArg(4, vnested);
				ctx->setRet(0, stencil_size);
				
				comp.unuse(address);
				comp.unuse(self);
				comp.unuse(velem);
				comp.unuse(vuser_data);
				comp.unuse(stencil_id);
				comp.unuse(vnested);
			}

			comp.mov(index, imm(0));
			Label L_1(comp.newLabel());
			comp.bind(L_1);
			{
				GpVar address_elem(comp.newGpVar());
				GpVar self(comp.newGpVar(kVarTypeUIntPtr));
				GpVar vindex(comp.newGpVar(GetEnumType()));
				GpVar vnested(comp.newGpVar(GetEnumType()));
				comp.mov(self, imm_ptr(this));
				comp.mov(vindex, index);
				comp.mov(vnested, imm(nested));
				comp.mov(address_elem, imm(reinterpret_cast<int64_t>(static_cast<Storage *(*)(Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetStencilElemJIT))));
				X86X64CallNode *ctx = comp.call(address_elem, kFuncConvHost, FuncBuilder3<Storage *, Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, self);
				ctx->setArg(1, vindex);
				ctx->setArg(2, vnested);
				ctx->setRet(0, stencil_elem);
				comp.unuse(address_elem);
				comp.unuse(vindex);
				comp.unuse(vnested);
			}
			//std::cout << "for " << reg_stencils[var.op].name << " nested " << nested << std::endl;
			XmmVar v1 = GenerateEvaluateSubASMJIT(*var.left, comp, stencil_elem, user_data, nested+1);
			{
				XmmVar stencil_coef(comp.newXmmVar(GetRealTypeXmm()));
				GpVar self(comp.newGpVar(kVarTypeUIntPtr));
				GpVar vindex(comp.newGpVar(GetEnumType()));
				GpVar vnested(comp.newGpVar(GetEnumType()));
				GpVar address_coef(comp.newGpVar());
				comp.mov(self, imm_ptr(this));
				comp.mov(vindex, index);
				comp.mov(vnested, imm(nested));
				comp.mov(address_coef, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetStencilCoefJIT))));
				X86X64CallNode * ctx = comp.call(address_coef, kFuncConvHost, FuncBuilder3<INMOST_DATA_REAL_TYPE, Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, self);
				ctx->setArg(1, vindex);
				ctx->setArg(2, vnested);
				ctx->setRet(0, stencil_coef);
				comp.unuse(vindex);
				comp.unuse(vnested);
				comp.unuse(address_coef);
				if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(v1, stencil_coef); else comp.mulsd(v1, stencil_coef);
				comp.unuse(stencil_coef);
			}
			
			if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.addss(ret, v1); else comp.addsd(ret, v1);
			comp.unuse(v1);

			comp.inc(index);
			comp.cmp(index, stencil_size);
			comp.jne(L_1);

			comp.unuse(index);
			comp.unuse(stencil_size);
			comp.unuse(stencil_elem);
			

			{
				const int64_t opid = static_cast<const INMOST_DATA_ENUM_TYPE &>(var.op);
				GpVar address(comp.newGpVar());
				GpVar self(comp.newGpVar(kVarTypeUIntPtr));
				GpVar stencil_id(comp.newGpVar(GetEnumType()));
				GpVar vnested(comp.newGpVar(GetEnumType()));
				comp.mov(self, imm_ptr(this));
				comp.mov(stencil_id, imm(opid));
				comp.mov(vnested, imm(nested));
				comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<void(*)(Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE)>(Automatizator::ClearStencilJIT))));
				X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder3<void, Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, self);
				ctx->setArg(1, stencil_id);
				ctx->setArg(2, vnested);
				comp.unuse(address);
				comp.unuse(self);
				comp.unuse(stencil_id);
				comp.unuse(vnested);
			}

			MultCoef(var.coef, comp, ret);
			return ret;
		}
		if (var.op >= AD_CTAG)
		{
			XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
			GpVar address(comp.newGpVar());
			GpVar velem(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vtag(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vcomp(comp.newGpVar(GetEnumType()));
			comp.mov(velem, elem);
			comp.mov(vtag, imm_ptr(reinterpret_cast<void *>(&reg_ctags[var.op].t)));
			comp.mov(vcomp, imm(*(INMOST_DATA_ENUM_TYPE *)(&var.left)));
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Storage *, Tag *, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetValueJIT))));
			X86X64CallNode * ctx = comp.call(address, kFuncConvHost, FuncBuilder3<INMOST_DATA_REAL_TYPE, Storage *, Tag *, INMOST_DATA_ENUM_TYPE>());
			ctx->setArg(0, velem);
			ctx->setArg(1, vtag);
			ctx->setArg(2, vcomp);
			ctx->setRet(0,ret);
			comp.unuse(velem);
			comp.unuse(vtag);
			comp.unuse(vcomp);
			comp.unuse(address);
			MultCoef(var.coef, comp, ret);
			return ret;
		}
		if (var.op >= AD_TAG)
		{
			XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
			GpVar address(comp.newGpVar());
			GpVar velem(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vtag(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vcomp(comp.newGpVar(GetEnumType()));
			comp.mov(velem, elem);
			comp.mov(vtag, imm_ptr(reinterpret_cast<void *>(&reg_tags[var.op].d.t)));
			comp.mov(vcomp, imm(*(INMOST_DATA_ENUM_TYPE *)(&var.left)));
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Storage *, Tag *, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetValueJIT))));
			X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder3<INMOST_DATA_REAL_TYPE, Storage *, Tag *, INMOST_DATA_ENUM_TYPE>());
			ctx->setArg(0, velem);
			ctx->setArg(1, vtag);
			ctx->setArg(2, vcomp);
			ctx->setRet(0, ret);
			comp.unuse(velem);
			comp.unuse(vtag);
			comp.unuse(vcomp);
			comp.unuse(address);
			MultCoef(var.coef, comp, ret);
			return ret;
		}
		assert(false);
	}

	void Automatizator::PushCode(asmjit::host::Compiler & comp, asmjit::host::XmmVar & var)
	{
		using namespace asmjit;
		using namespace asmjit::host;
		GpVar self(comp.newGpVar(kVarTypeUIntPtr));
		GpVar address(comp.newGpVar());
		//comp.spill(var);
		XmmVar tmp(comp.newXmmVar(GetRealTypeXmm()));
		if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(tmp, var); else comp.movsd(tmp, var);
		comp.mov(self, imm_ptr(this));
		comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<void(*)(INMOST_DATA_REAL_TYPE, Automatizator *)>(Automatizator::PushJIT))));
		X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder2<void, INMOST_DATA_REAL_TYPE, Automatizator *>());
		ctx->setArg(0, tmp);
		ctx->setArg(1, self);
		//ctx->setRet(0, tmp);
		//comp.alloc(var);
		comp.unuse(tmp);
		comp.unuse(self);
		comp.unuse(address);
	}


	void Automatizator::PopCode(asmjit::host::Compiler & comp, asmjit::host::XmmVar & var)
	{
		using namespace asmjit;
		using namespace asmjit::host;
		GpVar self(comp.newGpVar(kVarTypeUIntPtr));
		GpVar address(comp.newGpVar());
		XmmVar tmp(comp.newXmmVar(GetRealTypeXmm()));
		comp.mov(self, imm_ptr(this));
		comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Automatizator *)>(Automatizator::PopJIT))));
		X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder1<INMOST_DATA_REAL_TYPE, Automatizator *>());
		ctx->setArg(0, self);
		ctx->setRet(0, tmp);
		if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(var, tmp); else comp.movsd(var, tmp);
		comp.unuse(tmp);
		comp.unuse(self);
		comp.unuse(address);
	}

	asmjit::host::XmmVar   Automatizator::GenerateDerivativePrecomputeSubASMJIT(const expr & var, asmjit::host::Compiler & comp, asmjit::host::GpVar & elem, asmjit::host::GpVar & user_data, INMOST_DATA_ENUM_TYPE nested)
	{
		using namespace asmjit;
		using namespace asmjit::host;
		assert(var.op != AD_NONE);

		switch (var.op)
		{
		case AD_COND:
		{
						XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
						XmmVar v0 = GenerateEvaluateSubASMJIT(*var.left, comp, elem, user_data, nested);
						XmmVar zero(comp.newXmmVar(GetRealTypeXmm()));
						SetFloat(comp, zero, 0.0);
						Label Brench(comp.newLabel()), Exit(comp.newLabel());
						
						//GpVar vzero(comp.newGpVar(GetRealType()));
						//INMOST_DATA_REAL_TYPE zeroval = 0.0;
						//const int64_t * i = reinterpret_cast<const int64_t *>(&zeroval);
						//comp.mov(vzero, imm(i[0]));
						//comp.movq(zero, vzero);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.comiss(zero, v0); else comp.comisd(zero, v0); // zero -> zero <= v0
						//if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.cmpss(zero, v0, 2); else comp.cmpsd(zero, v0, 2); // zero -> zero <= v0
						//comp.movq(vzero, zero);
						//comp.cmp(vzero, imm(0));
						comp.ja(Brench);
						//comp.movq(ret, vzero); //redundant!!
						XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.right->left, comp, elem, user_data, nested);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(ret, v1); else comp.movsd(ret, v1);
						comp.unuse(v1);
						comp.jmp(Exit);
						comp.bind(Brench);
						//comp.movq(ret, vzero); //redundant!!
						XmmVar v2 = GenerateDerivativePrecomputeSubASMJIT(*var.right->right, comp, elem, user_data, nested);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(ret, v2); else comp.movsd(ret, v2);
						comp.unuse(v2);
						comp.bind(Exit);
						MultCoef(var.coef, comp, ret);
						//comp.unuse(vzero);
						comp.unuse(zero);
						PushCode(comp, v0);
						comp.unuse(v0);
						return ret;
		}
		case AD_PLUS:
		{
						XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, elem, user_data, nested), v2 = GenerateDerivativePrecomputeSubASMJIT(*var.right, comp, elem, user_data, nested);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.addss(v1, v2); else comp.addsd(v1, v2);
						comp.unuse(v2);
						MultCoef(var.coef, comp, v1);
						return v1;
		}
		case AD_MINUS:
		{
						 XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, elem, user_data, nested), v2 = GenerateDerivativePrecomputeSubASMJIT(*var.right, comp, elem, user_data, nested);
						 if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.subss(v1, v2); else comp.subsd(v1, v2);
						 comp.unuse(v2);
						 MultCoef(var.coef, comp, v1);
						 return v1;
		}
		case AD_MULT:
		{
						XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, elem, user_data, nested);
						comp.save(v1);
						XmmVar v2 = GenerateDerivativePrecomputeSubASMJIT(*var.right, comp, elem, user_data, nested);
						comp.save(v2);
						//comp.spill(v1);
						//comp.spill(v2);
						
						PushCode(comp, v1);
						PushCode(comp, v2);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(v1, v2); else comp.mulsd(v1, v2);
						MultCoef(var.coef, comp, v1);
						comp.unuse(v2);
						return v1;
		}
		case AD_DIV:
		{
					   XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, elem, user_data, nested), v2 = GenerateDerivativePrecomputeSubASMJIT(*var.right, comp, elem, user_data, nested);
					   PushCode(comp, v1);
					   PushCode(comp, v2);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.divss(v1, v2); else comp.divsd(v1, v2);
					   comp.unuse(v2);
					   return v1;
		}
		case AD_INV:
		{

					   XmmVar vconst(comp.newXmmVar(GetRealTypeXmm()));
					   SetFloat(comp, vconst, var.coef);
					   XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, elem, user_data, nested);
					   PushCode(comp, v1);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.divss(vconst, v1); else comp.divsd(vconst, v1);
					   comp.unuse(v1);
					   return vconst;
		}
		case AD_POW:
		{
					   XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, elem, user_data, nested), v2 = GenerateDerivativePrecomputeSubASMJIT(*var.right, comp, elem, user_data, nested);
					   PushCode(comp, v1);
					   PushCode(comp, v2);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float, float)>(::pow))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double, double)>(::pow))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder2<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setArg(0, v2);
					   ctx->setRet(0, v1);
					   PushCode(comp, v1);
					   comp.unuse(v2);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_ABS:
		{
					   XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, elem, user_data, nested);
					   PushCode(comp, v1);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::fabs))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::fabs))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_EXP:
		{
					   XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, elem, user_data, nested);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::exp))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::exp))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   PushCode(comp, v1);
					   return v1;
		}
		case AD_LOG:
		{
					   XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, elem, user_data, nested);
					   PushCode(comp, v1);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::log))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::log))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_SIN:
		{
					   XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, elem, user_data, nested);
					   PushCode(comp, v1);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::sin))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::sin))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_COS:
		{
					   XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, elem, user_data, nested);
					   PushCode(comp, v1);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::cos))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::cos))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   MultCoef(var.coef, comp, v1);
					   return v1;
		}
		case AD_CONST:
		{
						 XmmVar vconst(comp.newXmmVar(GetRealTypeXmm()));
						 SetFloat(comp,vconst, var.coef);
						 return vconst;
		}
		case AD_MES:
		{
					   XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
					   GpVar address(comp.newGpVar());
					   GpVar vmesh(comp.newGpVar(kVarTypeUIntPtr));
					   GpVar velem(comp.newGpVar(kVarTypeUIntPtr));
					   comp.mov(vmesh, imm(reinterpret_cast<int64_t>(m)));
					   comp.mov(velem, elem);
					   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Mesh *, Storage *)>(Automatizator::GetGeometricDataJIT))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, Mesh *, Storage *>());
					   ctx->setArg(0, vmesh);
					   ctx->setArg(1, velem);
					   ctx->setRet(0, ret);
					   comp.unuse(address);
					   comp.unuse(vmesh);
					   comp.unuse(velem);
					   MultCoef(var.coef, comp, ret);
					   return ret;
		}
		}
		if (var.op >= AD_FUNC)
		{
			XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
			GpVar address(comp.newGpVar());
			GpVar velem(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vuser_data(comp.newGpVar(kVarTypeUIntPtr));
			comp.mov(velem, elem);
			comp.mov(vuser_data, user_data);
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Storage *, void *)>(reg_funcs[var.op].func))));
			X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, Storage *, void *>());
			ctx->setArg(0, velem);
			ctx->setArg(1, vuser_data);
			ctx->setRet(0, ret);
			comp.unuse(address);
			comp.unuse(velem);
			comp.unuse(vuser_data);
			MultCoef(var.coef, comp, ret);
			return ret;
		}
		if (var.op >= AD_TABLE)
		{
			XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, elem, user_data, nested);
			PushCode(comp, v1);
			GpVar address(comp.newGpVar());
			GpVar vtable(comp.newGpVar(kVarTypeUIntPtr));
			comp.mov(vtable, imm(reinterpret_cast<int64_t>(reg_tables[var.op])));
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(INMOST_DATA_REAL_TYPE, table_ptr)>(Automatizator::GetTableValueJIT))));
			X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE, table_ptr>());
			ctx->setArg(0, v1);
			ctx->setArg(1, vtable);
			ctx->setRet(0, v1);
			comp.unuse(address);
			comp.unuse(vtable);
			MultCoef(var.coef, comp, v1);
			return v1;
		}
		if (var.op >= AD_STNCL)
		{
			XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
			SetFloat(comp, ret, 0.0);

			GpVar index(comp.newGpVar(GetEnumType()));
			GpVar stencil_size(comp.newGpVar(GetEnumType()));
			GpVar stencil_elem(comp.newGpVar(kVarTypeUIntPtr));



			{
				const int64_t opid = static_cast<const INMOST_DATA_ENUM_TYPE &>(var.op);
				GpVar address(comp.newGpVar());
				GpVar self(comp.newGpVar(kVarTypeUIntPtr));
				GpVar velem(comp.newGpVar(kVarTypeUIntPtr));
				GpVar vuser_data(comp.newGpVar(kVarTypeUIntPtr));
				GpVar stencil_id(comp.newGpVar(GetEnumType()));
				GpVar vnested(comp.newGpVar(GetEnumType()));

				comp.mov(self, imm_ptr(this));
				comp.mov(velem, elem);
				comp.mov(vuser_data, user_data);
				comp.mov(stencil_id, imm(opid));
				comp.mov(vnested, imm(nested));
				comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_ENUM_TYPE(*)(Automatizator *, Storage *, void *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetStencilJIT))));
				X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder5<INMOST_DATA_ENUM_TYPE, Automatizator *, Storage *, void *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, self);
				ctx->setArg(1, velem);
				ctx->setArg(2, vuser_data);
				ctx->setArg(3, stencil_id);
				ctx->setArg(4, vnested);
				ctx->setRet(0, stencil_size);

				comp.unuse(address);
				comp.unuse(self);
				comp.unuse(velem);
				comp.unuse(vuser_data);
				comp.unuse(stencil_id);
				comp.unuse(vnested);
			}
			comp.mov(index, imm(0));
			Label L_1(comp.newLabel());
			comp.bind(L_1);
			{
				GpVar address_elem(comp.newGpVar());
				GpVar self(comp.newGpVar(kVarTypeUIntPtr));
				GpVar vindex(comp.newGpVar(GetEnumType()));
				GpVar vnested(comp.newGpVar(GetEnumType()));
				comp.mov(self, imm_ptr(this));
				comp.mov(vindex, index);
				comp.mov(vnested, imm(nested));
				comp.mov(address_elem, imm(reinterpret_cast<int64_t>(static_cast<Storage *(*)(Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetStencilElemJIT))));
				X86X64CallNode *ctx = comp.call(address_elem, kFuncConvHost, FuncBuilder3<Storage *, Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, self);
				ctx->setArg(1, vindex);
				ctx->setArg(2, vnested);
				ctx->setRet(0, stencil_elem);
				comp.unuse(address_elem);
				comp.unuse(self);
				comp.unuse(vindex);
				comp.unuse(vnested);
			}
			//std::cout << "for " << reg_stencils[var.op].name << " nested " << nested << std::endl;
			XmmVar v1 = GenerateDerivativePrecomputeSubASMJIT(*var.left, comp, stencil_elem, user_data, nested + 1);
			{
				XmmVar stencil_coef(comp.newXmmVar(GetRealTypeXmm()));
				GpVar self(comp.newGpVar(kVarTypeUIntPtr));
				GpVar vindex(comp.newGpVar(GetEnumType()));
				GpVar vnested(comp.newGpVar(GetEnumType()));
				GpVar address_coef(comp.newGpVar());
				comp.mov(self, imm_ptr(this));
				comp.mov(vindex, index);
				comp.mov(vnested, imm(nested));
				comp.mov(address_coef, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetStencilCoefJIT))));
				X86X64CallNode * ctx = comp.call(address_coef, kFuncConvHost, FuncBuilder3<INMOST_DATA_REAL_TYPE, Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, self);
				ctx->setArg(1, vindex);
				ctx->setArg(2, vnested);
				ctx->setRet(0, stencil_coef);
				comp.unuse(vindex);
				comp.unuse(vnested);
				comp.unuse(address_coef);
				if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(v1, stencil_coef); else comp.mulsd(v1, stencil_coef);
				comp.unuse(stencil_coef);
			}

			if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.addss(ret, v1); else comp.addsd(ret, v1);
			comp.unuse(v1);

			comp.inc(index);
			comp.cmp(index, stencil_size);
			comp.jne(L_1);

			comp.unuse(index);
			comp.unuse(stencil_size);
			comp.unuse(stencil_elem);


			{
				const int64_t opid = static_cast<const INMOST_DATA_ENUM_TYPE &>(var.op);
				GpVar address(comp.newGpVar());
				GpVar self(comp.newGpVar(kVarTypeUIntPtr));
				GpVar stencil_id(comp.newGpVar(GetEnumType()));
				GpVar vnested(comp.newGpVar(GetEnumType()));
				comp.mov(self, imm_ptr(this));
				comp.mov(stencil_id, imm(opid));
				comp.mov(vnested, imm(nested));
				comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<void(*)(Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE)>(Automatizator::ClearStencilJIT))));
				X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder3<void, Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, self);
				ctx->setArg(1, stencil_id);
				ctx->setArg(2, vnested);
				comp.unuse(address);
				comp.unuse(self);
				comp.unuse(stencil_id);
				comp.unuse(vnested);
			}

			MultCoef(var.coef, comp, ret);
			return ret;
		}
		if (var.op >= AD_CTAG)
		{
			XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
			GpVar address(comp.newGpVar());
			GpVar velem(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vtag(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vcomp(comp.newGpVar(GetEnumType()));
			comp.mov(velem, elem);
			comp.mov(vtag, imm_ptr(reinterpret_cast<void *>(&reg_ctags[var.op].t)));
			comp.mov(vcomp, imm(*(INMOST_DATA_ENUM_TYPE *)(&var.left)));
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Storage *, Tag *, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetValueJIT))));
			X86X64CallNode * ctx = comp.call(address, kFuncConvHost, FuncBuilder3<INMOST_DATA_REAL_TYPE, Storage *, Tag *, INMOST_DATA_ENUM_TYPE>());
			ctx->setArg(0, velem);
			ctx->setArg(1, vtag);
			ctx->setArg(2, vcomp);
			ctx->setRet(0, ret);
			comp.unuse(velem);
			comp.unuse(vtag);
			comp.unuse(vcomp);
			comp.unuse(address);
			MultCoef(var.coef, comp, ret);
			return ret;
		}
		if (var.op >= AD_TAG)
		{
			XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
			GpVar address(comp.newGpVar());
			GpVar velem(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vtag(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vcomp(comp.newGpVar(GetEnumType()));
			comp.mov(velem, elem);
			comp.mov(vtag, imm_ptr(reinterpret_cast<void *>(&reg_tags[var.op].d.t)));
			comp.mov(vcomp, imm(*(INMOST_DATA_ENUM_TYPE *)(&var.left)));
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Storage *, Tag *, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetValueJIT))));
			X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder3<INMOST_DATA_REAL_TYPE, Storage *, Tag *, INMOST_DATA_ENUM_TYPE>());
			ctx->setArg(0, velem);
			ctx->setArg(1, vtag);
			ctx->setArg(2, vcomp);
			ctx->setRet(0, ret);
			comp.unuse(velem);
			comp.unuse(vtag);
			comp.unuse(vcomp);
			comp.unuse(address);
			MultCoef(var.coef, comp, ret);
			return ret;
		}
		assert(false);
	}


	void      Automatizator::GenerateDerivativeFillSubASMJIT(const expr & var, asmjit::host::Compiler & comp, asmjit::host::GpVar & elem, asmjit::host::GpVar & user_data, asmjit::host::GpVar & row, asmjit::host::XmmVar & multval, INMOST_DATA_ENUM_TYPE nested)
	{
		using namespace asmjit;
		using namespace asmjit::host;
		assert(var.op != AD_NONE);

		switch (var.op)
		{
		case AD_COND:
		{
						MultCoef(var.coef, comp, multval);
						XmmVar v0(comp.newXmmVar(GetRealTypeXmm()));
						PopCode(comp, v0);
						XmmVar zero(comp.newXmmVar(GetRealTypeXmm()));
						Label Brench(comp.newLabel()), Exit(comp.newLabel());
						SetFloat(comp, zero, 0.0);
						//GpVar vzero(comp.newGpVar(GetRealType()));
						//INMOST_DATA_REAL_TYPE zeroval = 0.0;
						//const int64_t * i = reinterpret_cast<const int64_t *>(&zeroval);
						//comp.mov(vzero, imm(i[0]));
						//comp.movq(zero, vzero);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.comiss(zero, v0); else comp.comisd(zero, v0); // zero -> zero <= v0
						//if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.cmpss(zero, v0, 2); else comp.cmpsd(zero, v0, 2); // zero -> zero <= v0
						//comp.movq(vzero, zero);
						//comp.cmp(vzero, imm(0));
						comp.ja(Brench);
						//comp.movq(v0, vzero); //redundant!!
						GenerateDerivativeFillSubASMJIT(*var.right->left, comp, elem, user_data,row, multval, nested);
						comp.jmp(Exit);
						comp.bind(Brench);
						//comp.movq(v0, vzero); //redundant!!
						GenerateDerivativeFillSubASMJIT(*var.right->right, comp, elem, user_data,row, multval, nested);
						comp.bind(Exit);
						comp.unuse(zero);
						comp.unuse(v0);
						return;
		}
		case AD_PLUS:
		{
						MultCoef(var.coef, comp, multval);
						GenerateDerivativeFillSubASMJIT(*var.right, comp, elem, user_data,row, multval, nested);
						GenerateDerivativeFillSubASMJIT(*var.left, comp, elem, user_data,row, multval, nested);
						return;
		}
		case AD_MINUS:
		{
						 MultCoef(-var.coef, comp, multval);
						 GenerateDerivativeFillSubASMJIT(*var.right, comp, elem, user_data,row, multval, nested);
						 MultCoef(-1.0, comp, multval);
						 GenerateDerivativeFillSubASMJIT(*var.left, comp, elem, user_data,row, multval, nested);
						 return;
		}
		case AD_MULT:
		{
						XmmVar v1(comp.newXmmVar(GetRealTypeXmm())), v2(comp.newXmmVar(GetRealTypeXmm()));
						PopCode(comp, v2);
						PopCode(comp, v1);
						MultCoef(var.coef, comp, multval);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(v1, multval); else comp.mulsd(v1, multval);
						if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(v2, multval); else comp.mulsd(v2, multval);
						GenerateDerivativeFillSubASMJIT(*var.right, comp, elem, user_data,row, v1, nested);
						comp.unuse(v1);
						GenerateDerivativeFillSubASMJIT(*var.left, comp, elem, user_data,row, v2, nested); 
						comp.unuse(v2);
						return;
		}
		case AD_DIV:
		{
					   XmmVar v1(comp.newXmmVar(GetRealTypeXmm())), v2(comp.newXmmVar(GetRealTypeXmm())), lv(comp.newXmmVar(GetRealTypeXmm())), rv0(comp.newXmmVar(GetRealTypeXmm())), rv(comp.newXmmVar(GetRealTypeXmm()));
					   PopCode(comp, v2);
					   PopCode(comp, v1);
					   MultCoef(var.coef, comp, multval);
					   SetFloat(comp, rv, 0.0);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(lv, multval); else comp.movsd(lv, multval);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.divss(lv, v2);      else comp.divsd(lv, v2);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(rv0, lv);     else comp.movsd(rv0, lv);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(rv0, v1);     else comp.mulsd(rv0, v1);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.divss(rv0, v2);     else comp.divsd(rv0, v2);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.subss(rv, rv0);     else comp.subsd(rv, rv0);
					   comp.unuse(v1);
					   comp.unuse(v2);
					   comp.unuse(rv0);
					   GenerateDerivativeFillSubASMJIT(*var.right, comp, elem, user_data,row, rv, nested);
					   comp.unuse(rv);
					   GenerateDerivativeFillSubASMJIT(*var.left, comp, elem, user_data,row, lv, nested);
					   comp.unuse(lv);
					   return;
		}
		case AD_INV:
		{
					   XmmVar v1(comp.newXmmVar(GetRealTypeXmm())), newmult(comp.newXmmVar(GetRealTypeXmm()));
					   PopCode(comp, v1);
					   SetFloat(comp,newmult, 0.0);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.subss(newmult, multval);     else comp.subsd(newmult, multval);
					   MultCoef(var.coef, comp, newmult);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.divss(newmult, v1);     else comp.divsd(newmult, v1);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.divss(newmult, v1);     else comp.divsd(newmult, v1);
					   comp.unuse(v1);
					   GenerateDerivativeFillSubASMJIT(*var.left, comp, elem, user_data,row, newmult, nested);
					   comp.unuse(newmult);
					   return;
		}
		case AD_POW:
		{
					   XmmVar v1(comp.newXmmVar(GetRealTypeXmm())), v2(comp.newXmmVar(GetRealTypeXmm())), ret(comp.newXmmVar(GetRealTypeXmm())), lv(comp.newXmmVar(GetRealTypeXmm())), rv(comp.newXmmVar(GetRealTypeXmm()));
					   PopCode(comp, ret);
					   PopCode(comp, v2);
					   PopCode(comp, v1);
					   MultCoef(var.coef, comp, multval);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(multval, ret); else comp.mulsd(multval,ret);
					   comp.unuse(ret);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(lv, multval); else comp.movsd(lv, multval);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(rv, multval); else comp.movsd(rv, multval);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(rv, v2);      else comp.mulsd(rv, v2);
					   comp.unuse(v2);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.divss(rv, v1);      else comp.divsd(rv, v1);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::log))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::log))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(lv, v1); else comp.mulsd(lv, v1);
					   comp.unuse(v1);
					   GenerateDerivativeFillSubASMJIT(*var.right, comp, elem, user_data,row, rv, nested);
					   comp.unuse(rv);
					   GenerateDerivativeFillSubASMJIT(*var.left, comp, elem, user_data,row, lv, nested);
					   comp.unuse(lv);
					   return;
		}
		case AD_ABS:
		{
					   XmmVar v1(comp.newXmmVar(GetRealTypeXmm()));
					   PopCode(comp, v1);
					   MultCoef(var.coef, comp, multval);
					   XmmVar zero(comp.newXmmVar(GetRealTypeXmm()));
					   SetFloat(comp, zero, 0.0);
					   Label Brench(comp.newLabel()), Exit(comp.newLabel());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.comiss(zero, v1); else comp.comisd(zero, v1); // zero -> zero <= v0
					   comp.unuse(v1);
					   comp.ja(Brench);
					   SetFloat(comp, zero, 1.0);
					   comp.jmp(Exit);
					   comp.bind(Brench);
					   SetFloat(comp, zero, -1.0);
					   comp.bind(Exit);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(multval, zero);      else comp.mulsd(multval, zero);
					   comp.unuse(zero);
					   GenerateDerivativeFillSubASMJIT(*var.left, comp, elem, user_data,row, multval, nested);
					   return;
		}
		case AD_EXP:
		{
					   XmmVar ret(comp.newXmmVar(GetRealTypeXmm()));
					   PopCode(comp, ret);
					   MultCoef(var.coef, comp, multval);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(multval, ret);      else comp.mulsd(multval, ret);
					   comp.unuse(ret);
					   GenerateDerivativeFillSubASMJIT(*var.left, comp, elem, user_data,row, multval, nested);
					   return;
		}
		case AD_LOG:
		{
					   XmmVar v1(comp.newXmmVar(GetRealTypeXmm()));
					   PopCode(comp, v1);
					   MultCoef(var.coef, comp, multval);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.divss(multval, v1);      else comp.divsd(multval, v1);
					   comp.unuse(v1);
					   GenerateDerivativeFillSubASMJIT(*var.left, comp, elem, user_data,row, multval, nested);
					   return;
		}
		case AD_SIN:
		{
					   XmmVar v1(comp.newXmmVar(GetRealTypeXmm()));
					   PopCode(comp, v1);
					   MultCoef(var.coef, comp, multval);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::cos))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::cos))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(multval, v1);      else comp.mulsd(multval, v1);
					   comp.unuse(v1);
					   GenerateDerivativeFillSubASMJIT(*var.left, comp, elem, user_data,row, multval, nested);
					   return;
		}
		case AD_COS:
		{
					   XmmVar v1(comp.newXmmVar(GetRealTypeXmm())), newmult(comp.newXmmVar(GetRealTypeXmm()));
					   SetFloat(comp, newmult, 0.0);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.subss(newmult, multval);      else comp.mulsd(newmult, multval);
					   MultCoef(var.coef, comp, newmult);
					   PopCode(comp, v1);
					   GpVar address(comp.newGpVar());
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4)
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<float(__cdecl*)(float)>(::sin))));
					   else
						   comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<double(__cdecl*)(double)>(::sin))));
					   X86X64CallNode *ctx = comp.call(address, kFuncConvHostCDecl, FuncBuilder1<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>());
					   ctx->setArg(0, v1);
					   ctx->setRet(0, v1);
					   comp.unuse(address);
					   if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(multval, v1);      else comp.mulsd(multval, v1);
					   comp.unuse(v1);
					   GenerateDerivativeFillSubASMJIT(*var.left, comp, elem, user_data,row, newmult, nested);
					   comp.unuse(newmult);
					   return;
		}
		case AD_CONST: return;
		case AD_MES: return;
		}
		if (var.op >= AD_FUNC)
		{
			return;
		}
		if (var.op >= AD_TABLE)
		{
			XmmVar v1(comp.newXmmVar(GetRealTypeXmm()));
			PopCode(comp, v1);
			GpVar address(comp.newGpVar());
			GpVar vtable(comp.newGpVar(kVarTypeUIntPtr));
			comp.mov(vtable, imm(reinterpret_cast<int64_t>(reg_tables[var.op])));
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(INMOST_DATA_REAL_TYPE, table_ptr)>(Automatizator::GetTableDerivativeJIT))));
			X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder2<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE, table_ptr>());
			ctx->setArg(0, v1);
			ctx->setArg(1, vtable);
			ctx->setRet(0, v1);
			comp.unuse(address);
			comp.unuse(vtable);
			if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(multval, v1);      else comp.mulsd(multval, v1);
			comp.unuse(v1);
			MultCoef(var.coef, comp, multval);
			GenerateDerivativeFillSubASMJIT(*var.left, comp, elem, user_data,row, multval, nested);
			return;
		}
		if (var.op >= AD_STNCL)
		{
			
			
			MultCoef(var.coef, comp, multval);
			GpVar index(comp.newGpVar(GetEnumType()));
			GpVar stencil_size(comp.newGpVar(GetEnumType()));
			GpVar stencil_elem(comp.newGpVar(kVarTypeUIntPtr));
			
			XmmVar newmult(comp.newXmmVar(GetRealTypeXmm()));

			
			{
				const int64_t opid = static_cast<const INMOST_DATA_ENUM_TYPE &>(var.op);
				GpVar address(comp.newGpVar());
				GpVar self(comp.newGpVar(kVarTypeUIntPtr));
				GpVar velem(comp.newGpVar(kVarTypeUIntPtr));
				GpVar vuser_data(comp.newGpVar(kVarTypeUIntPtr));
				GpVar stencil_id(comp.newGpVar(GetEnumType()));
				GpVar vnested(comp.newGpVar(GetEnumType()));

				comp.mov(self, imm_ptr(this));
				comp.mov(velem, elem);
				comp.mov(vuser_data, user_data);
				comp.mov(stencil_id, imm(opid));
				comp.mov(vnested, imm(nested));
				comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_ENUM_TYPE(*)(Automatizator *, Storage *, void *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetStencilJIT))));
				X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder5<INMOST_DATA_ENUM_TYPE, Automatizator *, Storage *, void *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, self);
				ctx->setArg(1, velem);
				ctx->setArg(2, vuser_data);
				ctx->setArg(3, stencil_id);
				ctx->setArg(4, vnested);
				ctx->setRet(0, stencil_size);

				comp.unuse(address);
				comp.unuse(self);
				comp.unuse(velem);
				comp.unuse(vuser_data);
				comp.unuse(stencil_id);
				comp.unuse(vnested);
			}
			
			comp.mov(index, stencil_size);
			//comp.mov(index, imm(0));
			
			
			Label L_1(comp.newLabel());
			comp.bind(L_1);
			
			comp.dec(index);
			XmmVar stencil_coef(comp.newXmmVar(GetRealTypeXmm()));
			{

				GpVar self(comp.newGpVar(kVarTypeUIntPtr));
				GpVar vindex(comp.newGpVar(GetEnumType()));
				GpVar vnested(comp.newGpVar(GetEnumType()));
				GpVar address_coef(comp.newGpVar());
				comp.mov(self, imm_ptr(this));
				comp.mov(vindex, index);
				comp.mov(vnested, imm(nested));
				comp.mov(address_coef, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_REAL_TYPE(*)(Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetStencilCoefJIT))));
				X86X64CallNode * ctx = comp.call(address_coef, kFuncConvHost, FuncBuilder3<INMOST_DATA_REAL_TYPE, Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, self);
				ctx->setArg(1, vindex);
				ctx->setArg(2, vnested);
				ctx->setRet(0, stencil_coef);
				comp.unuse(self);
				comp.unuse(vindex);
				comp.unuse(vnested);
				comp.unuse(address_coef);
			}
			if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(newmult, multval);      else comp.movsd(newmult, multval);
			if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.mulss(newmult, stencil_coef); else comp.mulsd(newmult, stencil_coef);
			comp.unuse(stencil_coef);
			

			

			{
				GpVar address_elem(comp.newGpVar());
				GpVar self(comp.newGpVar(kVarTypeUIntPtr));
				GpVar vindex(comp.newGpVar(GetEnumType()));
				GpVar vnested(comp.newGpVar(GetEnumType()));
				comp.mov(self, imm_ptr(this));
				comp.mov(vindex, index);
				comp.mov(vnested, imm(nested));
				comp.mov(address_elem, imm(reinterpret_cast<int64_t>(static_cast<Storage *(*)(Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetStencilElemJIT))));
				X86X64CallNode *ctx = comp.call(address_elem, kFuncConvHost, FuncBuilder3<Storage *, Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, self);
				ctx->setArg(1, vindex);
				ctx->setArg(2, vnested);
				ctx->setRet(0, stencil_elem);
				comp.unuse(address_elem);
				comp.unuse(self);
				comp.unuse(vindex);
				comp.unuse(vnested);
			}
			
			
			
			

			


			
			
			GenerateDerivativeFillSubASMJIT(*var.left, comp, stencil_elem, user_data,row,newmult, nested + 1);
			
			comp.cmp(index, 0);
			comp.jne(L_1);
			
			comp.unuse(index);
			comp.unuse(stencil_size);
			comp.unuse(stencil_elem);
			

			{
				const int64_t opid = static_cast<const INMOST_DATA_ENUM_TYPE &>(var.op);
				GpVar address(comp.newGpVar());
				GpVar self(comp.newGpVar(kVarTypeUIntPtr));
				GpVar stencil_id(comp.newGpVar(GetEnumType()));
				GpVar vnested(comp.newGpVar(GetEnumType()));
				comp.mov(self, imm_ptr(this));
				comp.mov(stencil_id, imm(opid));
				comp.mov(vnested, imm(nested));
				comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<void(*)(Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE)>(Automatizator::ClearStencilJIT))));
				X86X64CallNode *ctx = comp.call(address, kFuncConvHost, FuncBuilder3<void, Automatizator *, INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>());
				ctx->setArg(0, self);
				ctx->setArg(1, stencil_id);
				ctx->setArg(2, vnested);
				comp.unuse(address);
				comp.unuse(self);
				comp.unuse(stencil_id);
				comp.unuse(vnested);
			}
			
			
			return;
		}
		if (var.op >= AD_CTAG)
		{
			return;
		}
		if (var.op >= AD_TAG)
		{
			MultCoef(var.coef, comp, multval);
			GpVar address(comp.newGpVar());
			GpVar vrow(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vindex(comp.newGpVar(GetEnumType()));
			GpVar velem(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vtag(comp.newGpVar(kVarTypeUIntPtr));
			GpVar vcomp(comp.newGpVar(GetEnumType()));
			comp.mov(velem, elem);
			comp.mov(vtag, imm_ptr(reinterpret_cast<void *>(&reg_tags[var.op].indices)));
			comp.mov(vcomp, imm(*(INMOST_DATA_ENUM_TYPE *)(&var.left)));
			comp.mov(address, imm(reinterpret_cast<int64_t>(static_cast<INMOST_DATA_ENUM_TYPE(*)(Storage *, Tag *, INMOST_DATA_ENUM_TYPE)>(Automatizator::GetIndexJIT))));
			X86X64CallNode * ctx = comp.call(address, kFuncConvHost, FuncBuilder3<INMOST_DATA_ENUM_TYPE, Storage *, Tag *, INMOST_DATA_ENUM_TYPE>());
			ctx->setArg(0, velem);
			ctx->setArg(1, vtag);
			ctx->setArg(2, vcomp);
			ctx->setRet(0, vindex);
			XmmVar tmp(comp.newXmmVar(GetRealTypeXmm()));
			if (sizeof(INMOST_DATA_REAL_TYPE) == 4) comp.movss(tmp, multval); else comp.movsd(tmp, multval);
			comp.mov(vrow, row);
			comp.mov(address, imm_ptr((static_cast<INMOST_DATA_REAL_TYPE(*)(INMOST_DATA_REAL_TYPE, INMOST_DATA_ENUM_TYPE, Solver::Row *)>(Automatizator::InsertPairJIT))));
			ctx = comp.call(address, kFuncConvHost, FuncBuilder3<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE, INMOST_DATA_ENUM_TYPE, Solver::Row *>());
			ctx->setArg(0, tmp);
			ctx->setArg(1, vindex);
			ctx->setArg(2, vrow);
			ctx->setRet(0,tmp);
			return;
		}
		assert(false);
	}

	INMOST_DATA_REAL_TYPE Automatizator::Evaluate(void * prog, Storage * elem, void * user_data)
	{
		typedef std::pair<void *, void *> input_t;
		typedef INMOST_DATA_REAL_TYPE(*Func)(Storage *, void *);
		input_t * inp = (input_t *)prog;
		Func func = asmjit_cast<Func>(inp->first);
		return func(elem, user_data);
	}

	INMOST_DATA_REAL_TYPE Automatizator::Derivative(void * prog, Storage * elem, Solver::Row & row, void * user_data)
	{
		typedef std::pair<void *, void *> input_t;
		typedef INMOST_DATA_REAL_TYPE(*Func)(Storage *, void *, Solver::Row *);
		input_t * inp = (input_t *)prog;
		Func func = asmjit_cast<Func>(inp->second);
		return func(elem, user_data,&row);
	}
#endif

};

#endif //USE_AUTODIFF
