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
#if defined(NEW_VERSION)
	
	//! returns offset from the end of precomputed z
	void Automatizator::DerivativeFill(expr & var, INMOST_DATA_ENUM_TYPE element, INMOST_DATA_ENUM_TYPE parent, Sparse::RowMerger & entries, INMOST_DATA_REAL_TYPE multval, void * user_data)
	{
		INMOST_DATA_ENUM_TYPE voffset = var.values_offset(element), doffset = var.derivatives_offset(element);
		INMOST_DATA_ENUM_TYPE k = static_cast<INMOST_DATA_ENUM_TYPE>(var.data.size()-1);
		Storage e = Storage(m,var.current_stencil[element].first);
		INMOST_DATA_REAL_TYPE lval, rval, ret;
		var.values[doffset+k] = multval;
		expr::expr_data * arr = &var.data[0], *it;
		//for (expr::data_type::reverse_iterator it = var.data.rbegin(); it != var.data.rend(); ++it)
		do
		{
			it = arr+k;
			switch (it->op)
			{
			case AD_EXT: 
				assert(parent != ENUMUNDEF); 
				it->left.e->values[it->left.e->derivatives_offset(parent)+it->right.i] += var.values[doffset+k]; 
				break;
			case AD_COND:
				{
					expr & next = var.values[voffset+it->left.i] > 0.0 ? *it->right.q->left.e : *it->right.q->right.e;
					DerivativeFill(next, 0, element, entries,var.values[doffset+k], user_data);
				}
				break;
			case AD_PLUS:  
				var.values[doffset+it->left.i] += var.values[doffset+k];
				var.values[doffset+it->right.i] += var.values[doffset+k];
				break;
			case AD_MINUS:
				var.values[doffset+it->left.i] += var.values[doffset+k];
				var.values[doffset+it->right.i] -= var.values[doffset+k];
				break;
			case AD_MULT:
				rval = var.values[voffset+it->right.i];
				lval = var.values[voffset+it->left.i];
				var.values[doffset+it->left.i] += var.values[doffset+k]*rval;
				var.values[doffset+it->right.i] += var.values[doffset+k]*lval;
				break;
			case AD_DIV:
				rval = var.values[voffset+it->right.i];
				lval = var.values[voffset+it->left.i];
				var.values[doffset+it->left.i] += var.values[doffset+k]/rval;
				var.values[doffset+it->right.i] -= var.values[doffset+k]*lval/(rval*rval);
				break;
			case AD_POW:
				ret = var.values[voffset+k];
				rval = var.values[voffset+it->right.i];
				lval = var.values[voffset+it->left.i];
				var.values[doffset+it->right.i] += var.values[doffset+k] * ::log(lval);
				var.values[doffset+it->left.i] += var.values[doffset+k] * ret * rval / lval;
				break;
			case AD_SQRT:
				ret = var.values[voffset+k];
				var.values[doffset+it->left.i] += 0.5 * var.values[doffset+k] / ret;
				break;
			case AD_ABS:
				lval = var.values[voffset+it->left.i];
				var.values[doffset+it->left.i] += var.values[doffset+k] * (lval > 0.0 ? 1.0 : -1.0);
				break;
			case AD_EXP:
				ret = var.values[voffset+k];
				var.values[doffset+it->left.i] += var.values[doffset+k] * ret;
				break;
			case AD_LOG:   
				lval = var.values[voffset+it->left.i];
				var.values[doffset+it->left.i] += var.values[doffset+k] / lval;
				break;
			case AD_SIN:   
				lval = var.values[voffset+it->left.i];
				var.values[doffset+it->left.i] += var.values[doffset+k] * ::cos(lval);
				break;
			case AD_COS:   
				lval = var.values[voffset+it->left.i];
				var.values[doffset+it->left.i] -= var.values[doffset+k] * ::sin(lval);
				break;
			case AD_COND_MARK: break;
			case AD_COND_TYPE: break;
			case AD_CONST: break;
			case AD_MES: break;
			case AD_VAL: 
				{
					expr & next = *it->right.e;
					next.current_stencil.resize(1);
					next.current_stencil[0] = stencil_pair(e->GetHandle(),1.0);
					next.resize_for_stencil();
					var.values[doffset+it->left.i] += var.values[doffset+k] * EvaluateSub(next,0,element,user_data); 
					break;
				}
			default:
				if (it->op >= AD_FUNC) {}
				else if (it->op >= AD_TABLE) 
				{
					lval = var.values[voffset+it->left.i];
					var.values[doffset+it->left.i] += var.values[doffset+k] *  reg_tables[it->op-AD_TABLE]->get_derivative(lval);
				}
				else if (it->op >= AD_STNCL)
				{
					for (INMOST_DATA_ENUM_TYPE j = 0; j < it->left.e->current_stencil.size(); ++j) if(  it->left.e->current_stencil[j].first != NULL )
					{
						DerivativeFill(*it->left.e,  j, ENUMUNDEF,  entries, var.values[doffset+k] *  it->left.e->current_stencil[j].second, user_data);
					}
				}
				else if (it->op >= AD_CTAG) {}
				else if (it->op >= AD_TAG)
				{
					if (isDynamicValid(e, it->op))
					{
						entries[GetDynamicIndex(e,it->op,it->left.i)] += var.values[doffset+k];
					}
				}
				else assert(false);
			}
		} while(k-- != 0);
	}
	INMOST_DATA_REAL_TYPE Automatizator::EvaluateSub(expr & var, INMOST_DATA_ENUM_TYPE element, INMOST_DATA_ENUM_TYPE parent, void * user_data)
	{
		INMOST_DATA_ENUM_TYPE k = 0, offset = var.values_offset(element);
		Storage e = Storage(m,var.current_stencil[element].first);
		expr::expr_data * arr = &var.data[0], *it;
		//for (expr::data_type::iterator it = var.data.begin(); it != var.data.end(); ++it)
		do
		{
			it = arr+k;
			switch (it->op)
			{
			case AD_EXT: 
				assert(parent != ENUMUNDEF); 
				var.values[offset+k] = it->left.e->values[it->left.e->values_offset(parent)+it->right.i]; 
				break;
			case AD_COND:
				{
					expr & next = var.values[offset+it->left.i] > 0.0 ? *it->right.q->left.e : *it->right.q->right.e;
					next.current_stencil.resize(1);
					next.current_stencil[0] = stencil_pair(e->GetHandle(),1.0);
					next.resize_for_stencil();
					var.values[offset+k] = EvaluateSub(next, 0,element, user_data);
				}
				break;
			case AD_PLUS:  var.values[offset+k] = var.values[offset+it->left.i] + var.values[offset+it->right.i]; break;
			case AD_MINUS: var.values[offset+k] = var.values[offset+it->left.i] - var.values[offset+it->right.i]; break;
			case AD_MULT:  var.values[offset+k] = var.values[offset+it->left.i] * var.values[offset+it->right.i]; break;
			case AD_DIV:   var.values[offset+k] = var.values[offset+it->left.i] / var.values[offset+it->right.i]; break;
			case AD_POW:   var.values[offset+k] = ::pow(var.values[offset+it->left.i], var.values[offset+it->right.i]); break;
			case AD_SQRT:  var.values[offset+k] = ::sqrt(var.values[offset+it->left.i]); break;
			case AD_ABS:   var.values[offset+k] = ::fabs(var.values[offset+it->left.i]); break;
			case AD_EXP:   var.values[offset+k] = ::exp(var.values[offset+it->left.i]); break;
			case AD_LOG:   var.values[offset+k] = ::log(var.values[offset+it->left.i]); break;
			case AD_SIN:   var.values[offset+k] = ::sin(var.values[offset+it->left.i]); break;
			case AD_COS:   var.values[offset+k] = ::cos(var.values[offset+it->left.i]); break;
			case AD_CONST: var.values[offset+k] = it->left.r; break;
			case AD_COND_TYPE: var.values[offset+k] = ((e->GetElementType() & it->left.i)? 1.0 : -1.0); break;
			case AD_COND_MARK: var.values[offset+k] = e->GetMarker(it->left.i) ? 1.0 : -1.0; break;
			case AD_MES: assert(!(e->GetElementType() & (ESET | MESH))); m->GetGeometricData(e->GetHandle(), MEASURE, &var.values[offset+k]); break;
			case AD_VAL: var.values[offset+k] = var.values[offset+it->left.i]; break;
			default:
				if (it->op >= AD_FUNC) var.values[offset+k] = reg_funcs[it->op-AD_FUNC].func(e, user_data);
				else if (it->op >= AD_TABLE) var.values[offset+k] = reg_tables[it->op-AD_TABLE]->get_value(var.values[offset+it->left.i]);
				else if (it->op >= AD_STNCL)
				{
					it->left.e->current_stencil.clear();
					GetStencil(it->op,e,user_data,it->left.e->current_stencil);
					it->left.e->resize_for_stencil();
					var.values[offset+k] = 0.0;
					for (INMOST_DATA_ENUM_TYPE j = 0; j < it->left.e->current_stencil.size(); ++j) if( it->left.e->current_stencil[j].first != NULL )
						var.values[offset+k] += EvaluateSub(*it->left.e,j,ENUMUNDEF,user_data) * it->left.e->current_stencil[j].second;
				}
				else if (it->op >= AD_CTAG) var.values[offset+k] = GetStaticValue(e, it->op, it->left.i);
				else if (it->op >= AD_TAG) var.values[offset+k] = GetDynamicValue(e, it->op, it->left.i);
				else assert(false);
			}
			//k++;
		} while(++k != var.data.size());
		return var.values[offset+var.data.size()-1];
	}

	INMOST_DATA_REAL_TYPE Automatizator::Evaluate(expr & var, const Storage & e, void * user_data)
	{
		var.current_stencil.resize(1);
		var.current_stencil[0] = stencil_pair(e->GetHandle(),1.0);
		var.resize_for_stencil();
		return EvaluateSub(var,0,ENUMUNDEF,user_data);
	}
	
	INMOST_DATA_REAL_TYPE Automatizator::Derivative(expr & var, const Storage & e, Sparse::Row & out, Storage::real multiply, void * user_data)
	{
    Sparse::RowMerger & m = GetMerger();
		INMOST_DATA_REAL_TYPE ret;
		var.current_stencil.resize(1);
		var.current_stencil[0] = stencil_pair(e->GetHandle(),1.0);
		var.resize_for_stencil();
		ret = EvaluateSub(var,0,ENUMUNDEF,user_data);
    m.PushRow(1.0,out);
		DerivativeFill(var, 0, ENUMUNDEF, m, multiply, user_data);
    m.RetriveRow(out);
    m.Clear();
		return ret*multiply;
	}
#else

	INMOST_DATA_REAL_TYPE Automatizator::DerivativePrecompute(const expr & var, const Storage & e, precomp_values_t & values, void * user_data)
	{
		assert(var.op != AD_NONE);
		INMOST_DATA_REAL_TYPE lval, rval, ret = 0.0;
		switch (var.op)
		{
		case AD_COND_TYPE:
			lval = (e->GetElementType() & (*(ElementType*)&var.left))? 1.0 : -1.0;
			return lval;
		case AD_COND_MARK:
			lval = e->GetMarker(*(MarkerType *)&var.left)? 1.0 : -1.0;
			return lval;
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
		case AD_SQRT:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			ret = ::sqrt(lval);
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
			m->GetGeometricData(e->GetHandle(), MEASURE, &ret);
			return ret*var.coef;
		case AD_VAL:
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			rval = Evaluate(*var.right,e,user_data);
			values.push_back(rval);
			return lval*var.coef;
		}
		if (var.op >= AD_FUNC)
		{
			ret = reg_funcs[var.op-AD_FUNC].func(e, user_data);
			return ret*var.coef;
		}
		else if (var.op >= AD_TABLE)
		{
			lval = DerivativePrecompute(*var.left, e, values, user_data);
			values.push_back(lval);
			ret = reg_tables[var.op-AD_TABLE]->get_value(lval);
			return ret*var.coef;
		}
		else if (var.op >= AD_STNCL)
		{
			stencil_kind_domain & st = reg_stencils[var.op-AD_STNCL];
			assert(st.domainmask == 0 || e->GetMarker(st.domainmask));
			if (st.kind == 0)
			{
				Storage::reference_array elems = e->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
				Storage::real_array coefs = e->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
				assert(elems.size() == coefs.size());
				for (INMOST_DATA_ENUM_TYPE k = 0; k < elems.size(); ++k) if( elems.at(k) != InvalidHandle() )
				{
					lval = DerivativePrecompute(*var.left, elems[k], values, user_data);
					ret += lval * coefs[k];
				}
			}
			else if (st.kind == 1)
			{
				stencil_pairs get_st;
				reinterpret_cast<stencil_callback>(st.link)(e, get_st, user_data);
				for (INMOST_DATA_ENUM_TYPE k = 0; k < get_st.size(); ++k) if( get_st[k].first != InvalidHandle() )
				{
					lval = DerivativePrecompute(*var.left, Storage(m,get_st[k].first), values, user_data);
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
	void Automatizator::DerivativeFill(const expr & var, const Storage & e, Sparse::RowMerger & entries, precomp_values_t & values, INMOST_DATA_REAL_TYPE multval, void * user_data)
	{
		assert(var.op != AD_NONE);
		INMOST_DATA_REAL_TYPE lval, rval, ret;
		switch (var.op)
		{
		case AD_COND_MARK:
		case AD_COND_TYPE:
			return;
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
		case AD_SQRT:
			ret = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, 0.5 * multval / ret * var.coef, user_data);
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
		case AD_VAL: 
			rval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, multval * var.coef * rval, user_data);
			return;
		}
		if (var.op >= AD_FUNC)
		{
			return;
		}
		else if (var.op >= AD_TABLE)
		{
			lval = values.back(); values.pop_back();
			DerivativeFill(*var.left, e, entries, values, multval * var.coef * reg_tables[var.op-AD_TABLE]->get_derivative(lval), user_data);
			return;
		}
		else if (var.op >= AD_STNCL)
		{
			stencil_kind_domain & st = reg_stencils[var.op-AD_STNCL];
			assert(st.domainmask == 0 || e->GetMarker(st.domainmask));
			if (st.kind == 0)
			{
				Storage::reference_array elems = e->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
				Storage::real_array coefs = e->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
				assert(elems.size() == coefs.size());
				for (INMOST_DATA_ENUM_TYPE k = elems.size(); k > 0; --k) if( elems.at(k-1) != InvalidHandle() )
					DerivativeFill(*var.left, elems[k - 1], entries, values, var.coef * coefs[k - 1] * multval, user_data);
			}
			else if (st.kind == 1)
			{
				stencil_pairs get_st;
				reinterpret_cast<stencil_callback>(st.link)(e, get_st, user_data);
				for (INMOST_DATA_ENUM_TYPE k = static_cast<INMOST_DATA_ENUM_TYPE>(get_st.size()); k > 0; --k) if( get_st[k-1].first != InvalidHandle() )
					DerivativeFill(*var.left, Storage(m,get_st[k - 1].first), entries, values, var.coef * get_st[k - 1].second*multval, user_data);
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
	INMOST_DATA_REAL_TYPE Automatizator::Evaluate(const expr & var, const Storage & e, void * user_data)
	{
		assert(var.op != AD_NONE);
		switch (var.op)
		{
		case AD_COND_MARK: return e->GetMarker(*(MarkerType *)&var.left) ? 1.0 : -1.0;
		case AD_COND_TYPE: return (e->GetElementType() & (*(ElementType *)&var.left)) ? 1.0 : -1.0;
		case AD_COND:  return Evaluate(*(Evaluate(*var.left, e, user_data) > 0.0 ? var.right->left : var.right->right), e, user_data)*var.coef;
		case AD_PLUS:  return (Evaluate(*var.left, e, user_data) + Evaluate(*var.right, e, user_data))*var.coef;
		case AD_MINUS: return (Evaluate(*var.left, e, user_data) - Evaluate(*var.right, e, user_data))*var.coef;
		case AD_MULT:  return (Evaluate(*var.left, e, user_data) * Evaluate(*var.right, e, user_data))*var.coef;
		case AD_DIV:   return (Evaluate(*var.left, e, user_data) / Evaluate(*var.right, e, user_data))*var.coef;
		case AD_INV:   return var.coef / Evaluate(*var.left, e, user_data);
		case AD_POW:   return ::pow(Evaluate(*var.left, e, user_data), Evaluate(*var.right, e, user_data))*var.coef;
		case AD_SQRT:  return ::sqrt(Evaluate(*var.left, e, user_data))*var.coef;
		case AD_ABS:   return ::fabs(Evaluate(*var.left, e, user_data))*var.coef;
		case AD_EXP:   return ::exp(Evaluate(*var.left, e, user_data))*var.coef;
		case AD_LOG:   return ::log(Evaluate(*var.left, e, user_data))*var.coef;
		case AD_SIN:   return ::sin(Evaluate(*var.left, e, user_data))*var.coef;
		case AD_COS:   return ::cos(Evaluate(*var.left, e, user_data))*var.coef;
		case AD_CONST: return var.coef;
		case AD_MES: assert(!(e->GetElementType() & (ESET | MESH))); Storage::real ret; m->GetGeometricData(e->GetHandle(), MEASURE, &ret); return ret*var.coef;
		case AD_VAL:   return Evaluate(*var.left,e,user_data)*var.coef;
		}
		if (var.op >= AD_FUNC) return reg_funcs[var.op-AD_FUNC].func(e, user_data);
		if (var.op >= AD_TABLE) return reg_tables[var.op-AD_TABLE]->get_value(Evaluate(*var.left, e, user_data))*var.coef;
		if (var.op >= AD_STNCL)
		{
			INMOST_DATA_REAL_TYPE ret = 0.0;
			stencil_kind_domain & st = reg_stencils[var.op-AD_STNCL];
			assert(st.domainmask == 0 || e->GetMarker(st.domainmask));
			if (st.kind == 0)
			{
				Storage::reference_array elems = e->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
				Storage::real_array coefs = e->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
				assert(elems.size() == coefs.size());
				for (INMOST_DATA_ENUM_TYPE k = 0; k < elems.size(); ++k) if( elems.at(k) != InvalidHandle() )
					ret += var.coef * Evaluate(*var.left, elems[k], user_data) * coefs[k];
			}
			else if (st.kind == 1)
			{
				stencil_pairs get_st;
				reinterpret_cast<stencil_callback>(st.link)(e, get_st, user_data);
				for (INMOST_DATA_ENUM_TYPE k = 0; k < get_st.size(); ++k) if ( get_st[k].first != InvalidHandle() )
					ret += var.coef * Evaluate(*var.left, Storage(m,get_st[k].first), user_data) * get_st[k].second;
			}
			return ret;
		}
		if (var.op >= AD_CTAG) return GetStaticValue(e, var.op, *(INMOST_DATA_ENUM_TYPE *)(&var.left))*var.coef;
		if (var.op >= AD_TAG) return  GetDynamicValue(e, var.op, *(INMOST_DATA_ENUM_TYPE *)(&var.left))*var.coef;
		assert(false);
		return 0.0;
	}

	INMOST_DATA_REAL_TYPE Automatizator::Derivative(const expr & var, const Storage & e, Sparse::Row & out, Storage::real multiply, void * user_data)
	{
    Sparse::RowMerger & m = GetMerger();
		INMOST_DATA_REAL_TYPE ret;
		precomp_values_t values;
		ret = DerivativePrecompute(var, e, values, user_data);
    m.PushRow(1.0,out);
		DerivativeFill(var, e, m, values, multiply, user_data);
    m.RetriveRow(out);
    m.Clear();
		return ret*multiply;
	}
#endif
	Automatizator::Automatizator(Mesh * m) :first_num(0), last_num(0), m(m) {}
	Automatizator::~Automatizator()
	{
		for (unsigned k = 0; k < index_tags.size(); k++)
			index_tags[k].indices = m->DeleteTag(index_tags[k].indices);
		for (table_type::iterator it = reg_tables.begin(); it != reg_tables.end(); ++it)
		{
			delete[] (*it)->args;
			delete[] (*it)->vals;
			delete (*it);
		}
		for (stencil_type::iterator it = reg_stencils.begin(); it != reg_stencils.end(); ++it)
		if (it->kind == 0)
			delete static_cast<stencil_tag *>(it->link);
	}
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterFunc(std::string name, func_callback func)
	{
		INMOST_DATA_ENUM_TYPE ret = static_cast<INMOST_DATA_ENUM_TYPE>(reg_funcs.size()) + AD_FUNC;
		func_name_callback v;
		v.name = name;
		v.func = func;
		reg_funcs.push_back(v);
		return ret;
	}
	//register stencil that can be got from tags
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterStencil(std::string name, Tag elements_tag, Tag coefs_tag, MarkerType domain_mask)
	{
		INMOST_DATA_ENUM_TYPE ret = static_cast<INMOST_DATA_ENUM_TYPE>(reg_stencils.size()) + AD_STNCL;
		stencil_kind_domain st;
		stencil_tag * save = new stencil_tag;
		st.name = name;
		save->coefs = coefs_tag;
		save->elements = elements_tag;
		st.kind = 0;
		st.link = static_cast<void *>(save);
		st.domainmask = domain_mask;
		reg_stencils.push_back(st);
		return ret;
	}
	//register stencil that can be got from function
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterStencil(std::string name, stencil_callback func, MarkerType domain_mask)
	{
		INMOST_DATA_ENUM_TYPE ret = static_cast<INMOST_DATA_ENUM_TYPE>(reg_stencils.size()) + AD_STNCL;
		stencil_kind_domain st;
		st.name = name;
		st.kind = 1;
		st.link = reinterpret_cast<void *>(func);
		st.domainmask = domain_mask;
		reg_stencils.push_back(st);
		return ret;
	}
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterTable(std::string name, INMOST_DATA_REAL_TYPE * Arguments, INMOST_DATA_REAL_TYPE * Values, INMOST_DATA_ENUM_TYPE size)
	{
		INMOST_DATA_ENUM_TYPE ret = static_cast<INMOST_DATA_ENUM_TYPE>(reg_tables.size()) + AD_TABLE;
		table_ptr t = new table;
		t->name = name;
		t->args = new INMOST_DATA_REAL_TYPE[size];
		memcpy(t->args, Arguments, sizeof(INMOST_DATA_REAL_TYPE)*size);
		t->vals = new INMOST_DATA_REAL_TYPE[size];
		memcpy(t->vals, Values, sizeof(INMOST_DATA_REAL_TYPE)*size);
		t->size = size;
		reg_tables.push_back(t);
		return ret;
	}
	/// set data of tag t defined on domain_mask to be dynamic data
	/// don't register tag twice
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
		INMOST_DATA_ENUM_TYPE ret = static_cast<INMOST_DATA_ENUM_TYPE>(reg_tags.size()) + AD_TAG;
		reg_tags.push_back(p);
		index_tags.push_back(p);
		return ret;
	}
	/// set index for every data entry of dynamic tag
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
	/// register tag, data for which don't change through iterations
	/// don't register tag twice
	INMOST_DATA_ENUM_TYPE Automatizator::RegisterStaticTag(Tag t, MarkerType domain_mask)
	{
		INMOST_DATA_ENUM_TYPE ret = static_cast<INMOST_DATA_ENUM_TYPE>(reg_ctags.size()) + AD_CTAG;
		tagdomain d;
		d.t = t;
		d.domain_mask = domain_mask;
		reg_ctags.push_back(d);
		return ret;
	}

	INMOST_DATA_ENUM_TYPE Automatizator::GetStencil(INMOST_DATA_ENUM_TYPE stnclind, const Storage & elem, void * user_data, stencil_pairs & ret)
	{
		stencil_kind_domain & st = reg_stencils[stnclind-AD_STNCL];
		assert(st.domainmask == 0 || elem->GetMarker(st.domainmask));
		if (st.kind == 0)
		{
			Storage::reference_array elems = elem->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
			Storage::real_array coefs = elem->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
			assert(elems.size() == coefs.size());
			for (INMOST_DATA_ENUM_TYPE k = 0; k < elems.size(); ++k) 
				ret.push_back(std::make_pair(elems.at(k), coefs[k]));
		}
		else if (st.kind == 1) reinterpret_cast<stencil_callback>(st.link)(elem, ret, user_data);
		return static_cast<INMOST_DATA_ENUM_TYPE>(ret.size());
	}

	INMOST_DATA_ENUM_TYPE Automatizator::table::binary_search(INMOST_DATA_REAL_TYPE arg)
	{
		int l = 0, r = static_cast<int>(size)-1, mid = 0;
		while (r >= l)
		{
			mid = (l + r) / 2;
			if (args[mid] > arg) r = mid - 1;
			else if (args[mid] < arg) l = mid + 1;
			else return mid;
		}
		mid = (l + r) / 2;
		if (mid > static_cast<int>(size - 2)) mid = static_cast<int>(size - 2);
		return static_cast<INMOST_DATA_ENUM_TYPE>(mid);
	}
	INMOST_DATA_REAL_TYPE Automatizator::table::get_value(INMOST_DATA_REAL_TYPE arg)
	{
		if (arg < args[0]) return vals[0];
		INMOST_DATA_ENUM_TYPE i = binary_search(arg);
		return vals[i] + (vals[i + 1] - vals[i]) * (arg - args[i]) / (args[i + 1] - args[i]);
	}
	INMOST_DATA_REAL_TYPE Automatizator::table::get_derivative(INMOST_DATA_REAL_TYPE arg)
	{
		if (arg < args[0]) return 0.0;
		INMOST_DATA_ENUM_TYPE i = binary_search(arg);
		return (vals[i + 1] - vals[i]) / (args[i + 1] - args[i]);
	}
	std::pair<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> Automatizator::table::get_both(INMOST_DATA_REAL_TYPE arg)
	{
		if (arg < args[0]) return std::make_pair(vals[0], 0.0);
		INMOST_DATA_ENUM_TYPE i = binary_search(arg);
		INMOST_DATA_REAL_TYPE der = (vals[i + 1] - vals[i]) / (args[i + 1] - args[i]);
		return std::make_pair(vals[i] + der * (arg - args[i]), der);
	}
};

#endif //USE_AUTODIFF
