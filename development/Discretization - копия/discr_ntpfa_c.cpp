#include "discr.h"
#include <iomanip>
#define EPS 1e-12

#define ALLOWED_INTERPOLATION (AVG_HARMONIC | AVG_NONLINEAR | AVG_REVDIST)

typedef Storage::real vec[3];

static Storage::real lengthS(vec x) { return dot_prod(x, x); }
static Storage::real length(vec x) { return sqrt(lengthS(x)); }
static Storage::real normalize(Storage::real x[3]) { Storage::real r = length(x); if (r > 0.0)  x[0] /= r, x[1] /= r, x[2] /= r; return r; }

ntpfa_c::ntpfa_c(Automatizator * aut, Mesh * m, std::string tensor)
:discr_triplet_basic(aut, m, tensor)
{
	trpl1_elems = m->CreateTag("NTPFA_TRP_L1", DATA_REFERENCE, FACE, NONE, 3);
	trpr1_elems = m->CreateTag("NTPFA_TRP_R1", DATA_REFERENCE, FACE, NONE, 3);
	trpl1_coefs = m->CreateTag("NTPFA_COEF_L1", DATA_REAL, FACE, NONE, 3);
	trpr1_coefs = m->CreateTag("NTPFA_COEF_R1", DATA_REAL, FACE, NONE, 3);
	trpl1_add = m->CreateTag("NTPFA_ADD_L1", DATA_REAL, FACE, NONE, 2);
	trpr1_add = m->CreateTag("NTPFA_ADD_R1", DATA_REAL, FACE, NONE, 2);
	trpl2_elems = m->CreateTag("NTPFA_TRP_L2", DATA_REFERENCE, FACE, NONE, 3);
	trpr2_elems = m->CreateTag("NTPFA_TRP_R2", DATA_REFERENCE, FACE, NONE, 3);
	trpl2_coefs = m->CreateTag("NTPFA_COEF_L2", DATA_REAL, FACE, NONE, 3);
	trpr2_coefs = m->CreateTag("NTPFA_COEF_R2", DATA_REAL, FACE, NONE, 3);
	trpl2_add = m->CreateTag("NTPFA_ADD_L2", DATA_REAL, FACE, NONE, 2);
	trpr2_add = m->CreateTag("NTPFA_ADD_R2", DATA_REAL, FACE, NONE, 2);

	intrp_stncl = m->CreateTag("INTRP_STNCL", DATA_REFERENCE, EDGE | NODE, NONE);
	intrp_coefs = m->CreateTag("INTRP_COEFS", DATA_REAL,  EDGE | NODE, NONE);
	
	
	trpr1 = aut->RegisterStencil("ntpfa_right_triplet1", trpr1_elems, trpr1_coefs);
	trpl1 = aut->RegisterStencil("ntpfa_left_triplet1", trpl1_elems, trpl1_coefs);
	trpr2 = aut->RegisterStencil("ntpfa_right_triplet2", trpr2_elems, trpr2_coefs);
	trpl2 = aut->RegisterStencil("ntpfa_left_triplet2", trpl2_elems, trpl2_coefs);
	dwn = aut->RegisterStencil("back_cell", get_l);
	upw = aut->RegisterStencil("front_cell", get_r);
	addl1 = aut->RegisterStaticTag(trpl1_add);
	addr1 = aut->RegisterStaticTag(trpr1_add);
	addl2 = aut->RegisterStaticTag(trpl2_add);
	addr2 = aut->RegisterStaticTag(trpr2_add);

	intrp = aut->RegisterStencil("interp",intrp_stncl,intrp_coefs,add_markers);
}


void ntpfa_c::MarkRecursive(Element * f)
{
	Storage::reference_array lrtrp[4] = 
	{
		f->ReferenceArrayDF(trpl1_elems),
		f->ReferenceArrayDF(trpl2_elems),
		f->ReferenceArrayDF(trpr1_elems),
		f->ReferenceArrayDF(trpr2_elems)
	};
	int end = 2 + (f->GetMarker(bnd_markers)? 0 : 2);
	for(int k = 0; k < 2; k++)
	{
		for(Storage::reference_array::iterator it = lrtrp[k].begin(); it != lrtrp[k].end(); ++it) if( (*it) != NULL )
		{
			if( (*it)->GetElementType() != CELL && !(*it)->GetMarker(add_markers) )
			{
				(*it)->SetMarker(add_markers);
				(*it)->SetMarker(mark_next);
			}
		}
	}
	for(int k = 0; k < 2; k++)
	{
		for(Storage::reference_array::iterator it = lrtrp[k].begin(); it != lrtrp[k].end(); ++it) if( (*it) != NULL )
		{
			if( (*it)->GetMarker(mark_next) )
			{
				(*it)->RemMarker(mark_next);
				if((*it)->GetElementType() == FACE )// && (*it)->GetMarker(bnd_markers))
					MarkRecursive(*it);
				else
					MarkRecursiveEdge(*it);
			
			}
		}
	}
}

void ntpfa_c::MarkRecursiveEdge(Element * f)
{
	Storage::reference_array arr = f->ReferenceArrayDV(intrp_stncl);
	for(Storage::reference_array::iterator it = arr.begin(); it != arr.end(); ++it)
	{
		if( (*it)->GetElementType() != CELL && !(*it)->GetMarker(add_markers))
		{
			(*it)->SetMarker(add_markers);
			(*it)->SetMarker(mark_next);
		}
	}
	for(Storage::reference_array::iterator it = arr.begin(); it != arr.end(); ++it)
	{
		if( (*it)->GetMarker(mark_next) )
		{
			(*it)->RemMarker(mark_next);
			if((*it)->GetElementType() == FACE )// && (*it)->GetMarker(bnd_markers))
				MarkRecursive(*it);
			else
				MarkRecursiveEdge(*it);
			
		}
	}
}


void ntpfa_c::Init()
{
	if (!have_bnd) Boundary(0 , FACE | EDGE);
	for (INMOST_DATA_INTEGER_TYPE id = 0; id < m->MaxLocalID(FACE); ++id)
	{
		Face * f = m->FaceByLocalID(id);
		if (f != NULL) Build(f);
	}

	if( bnd_conds.isValid() )
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it)
			if( it->HaveData(bnd_conds) && it->GetMarker(bnd_markers) )
			{
				it->SetMarker(add_markers);
			}

	for(Mesh::iteratorElement it = m->BeginElement(EDGE | NODE); it != m->EndElement(); ++it)
	{
		Storage::real v[3], cnt[3], cnt0[3], l;
		it->Centroid(cnt0);
		Storage::real_array coefs = it->RealArray(intrp_coefs);
		Storage::reference_array stncl = it->ReferenceArray(intrp_stncl);
		adjacent<Element> cells = it->getAdjElements(it->GetMarker(bnd_markers) ? FACE : CELL);
		Storage::real sum = 0;
		for(int k = 0; k < cells.size(); ++k)
		if( !it->GetMarker(bnd_markers) || cells[k].GetMarker(bnd_markers) )
		{
			cells[k].Centroid(cnt);
			vec_diff(cnt0,cnt,v);
			l = 1.0/(dot_prod(v,v) + 1.0e-9);
			stncl.push_back(&cells[k]);
			coefs.push_back(l);
			sum += l;
		}
		for(int k = 0; k < coefs.size(); ++k) coefs[k] /= sum;
	}

	mark_next = m->CreateMarker();
	for(Mesh::iteratorElement it = m->BeginElement(FACE); it != m->EndElement(); ++it)
	{
		if( it->GetMarker(add_markers) )
		{
			MarkRecursive(&*it);
		}
	}
	for(Mesh::iteratorElement it = m->BeginElement(EDGE); it != m->EndElement(); ++it)
	{
		if( it->GetMarker(add_markers) )
		{
			MarkRecursiveEdge(&*it);
		}
	}
	m->ReleaseMarker(mark_next);


	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) it->SetMarker(add_markers);
	
	m->SynchronizeMarker(UnknownMarker(),FACE|EDGE|NODE,SYNC_BIT_OR);
}


void ntpfa_c::get_l(Storage * e, Automatizator::stencil_pairs & out, void * user_data)
{
	Face * f = static_cast<Face *>(e);
	Element * r = f->BackCell();
	out.push_back(Automatizator::stencil_pair(r == NULL ? f : r, 1.0));
}
void ntpfa_c::get_r(Storage * e, Automatizator::stencil_pairs & out, void * user_data)
{
	Face * f = static_cast<Face *>(e);
	Element * r = f->FrontCell();
	out.push_back(Automatizator::stencil_pair(r == NULL ? f : r, 1.0));
}


void ntpfa_c::Build(Face * f)
{
	if (needBuild(f))
	{
		Storage::real nrmf[3], area, area1, area2, f1[3], f2[3], cnt1[3], cnt2[3], cntf[3], deny;
		Element * r0 = f->BackCell();
		Element * r1 = f->FrontCell();

		
		r0->Centroid(cnt1);
		f->Centroid(cntf);
		f->UnitNormal(nrmf);
		area = f->Area();

		if( r1 != NULL )
		{
			r1->Centroid(cnt2);
			int ret1, ret2, ptypes[4][3];
			int ret11, ret12, ret21, ret22;
			ret1 = tensor_prod(Ktype, r0->RealArray(K), nrmf, f1);
			ret2 = tensor_prod(Ktype, (r1 == NULL ? r0 : r1)->RealArray(K), nrmf, f2);

			if (!ret1 || !ret2)	error = 1; //error in tensor

			Storage::reference_array ltrp1 = f->ReferenceArrayDF(trpl1_elems);
			Storage::reference_array ltrp2 = f->ReferenceArrayDF(trpl2_elems);
			Storage::reference_array rtrp1 = f->ReferenceArrayDF(trpr1_elems);
			Storage::reference_array rtrp2 = f->ReferenceArrayDF(trpr2_elems);
			Storage::real_array lcoef1 = f->RealArrayDF(trpl1_coefs);
			Storage::real_array lcoef2 = f->RealArrayDF(trpl2_coefs);
			Storage::real_array rcoef1 = f->RealArrayDF(trpr1_coefs);
			Storage::real_array rcoef2 = f->RealArrayDF(trpr2_coefs);
			Storage::real_array ltrp1_add = f->RealArrayDF(trpl1_add);
			Storage::real_array ltrp2_add = f->RealArrayDF(trpl2_add);
			Storage::real_array rtrp1_add = f->RealArrayDF(trpr1_add);
			Storage::real_array rtrp2_add = f->RealArrayDF(trpr2_add);

			area1 = area*normalize(f1);
			area2 = area*normalize(f2);
			/*
			ret11 = find_triplet(r0, &*f, NULL, f1, +1, area1, &ltrp1[0], &lcoef1[0],ptypes[0], allowed_types);//, AVG_NONLINEAR | AVG_REVDIST);
			ret21 = find_triplet(r1, &*f, NULL, f2, -1, area2, &rtrp1[0], &rcoef1[0],ptypes[2], allowed_types);//, AVG_NONLINEAR | AVG_REVDIST);
			ret12 = find_triplet(&*f, r0, NULL, f1, -1, area1, &ltrp2[0], &lcoef2[0],ptypes[1], allowed_types);//, AVG_NONLINEAR | AVG_REVDIST);
			ret22 = find_triplet(&*f, r1, NULL, f2, +1, area2, &rtrp2[0], &rcoef2[0],ptypes[3], allowed_types);//, AVG_NONLINEAR | AVG_REVDIST);
			*/
			ret11 = find_triplet(r0, &*f, NULL, f1, +1, area1, &ltrp1[0], &lcoef1[0],ptypes[0], allowed_types,true);//, AVG_NONLINEAR | AVG_REVDIST);
			ret21 = find_triplet(r1, &*f, NULL, f2, -1, area2, &rtrp1[0], &rcoef1[0],ptypes[2], allowed_types,true);//, AVG_NONLINEAR | AVG_REVDIST);
			ret12 = find_triplet(&*f, r0, NULL, f1, -1, area1, &ltrp2[0], &lcoef2[0],ptypes[1], allowed_types,true);//, AVG_NONLINEAR | AVG_REVDIST);
			ret22 = find_triplet(&*f, r1, NULL, f2, +1, area2, &rtrp2[0], &rcoef2[0],ptypes[3], allowed_types,true);//, AVG_NONLINEAR | AVG_REVDIST);


			if (ret11 && ret12 && ret21 && ret22)
			{
				ltrp1_add[0] = ltrp1_add[1] = ltrp2_add[0] = ltrp2_add[1] = 0.0;
				rtrp1_add[0] = rtrp1_add[1] = rtrp2_add[0] = rtrp2_add[1] = 0.0;

				for (int q = 0; q < 3; q++) 
				{
					assert( lcoef1[q] >= -1.0e-9 && lcoef2[q] >= -1.0e-9 && rcoef1[q] >= -1.0e-9 && rcoef2[q] >= -1.0e-9);

					ltrp1_add[1] += lcoef1[q];
					ltrp2_add[1] += lcoef2[q];
					rtrp1_add[1] += rcoef1[q];
					rtrp2_add[1] += rcoef2[q];
				}
				for (int q = 0; q < 3; q++) 
				{
					if( ltrp1[q] != &*f && ltrp1[q]->GetElementType() != CELL ) ltrp1[q]->SetMarker(add_markers);
					if( ltrp2[q] != &*f && ltrp2[q]->GetElementType() != CELL ) ltrp2[q]->SetMarker(add_markers);
					if( rtrp1[q] != &*f && rtrp1[q]->GetElementType() != CELL ) rtrp1[q]->SetMarker(add_markers);
					if( rtrp2[q] != &*f && rtrp2[q]->GetElementType() != CELL ) rtrp2[q]->SetMarker(add_markers);
				}
				//Flux from C_K to C_f
				deny = 0.0;
				for (int q = 0; q < 3; q++)
				{
				
					if (fabs(lcoef1[q] / ltrp1_add[1]) < 1.0e-8) //contribution is less then 0.5%
					{
						ltrp1[q] = NULL;
						deny += lcoef1[q];
						lcoef1[q] = 0.0;
					}
					else
					{
						if (ltrp1[q] == &*f)
						{
							ltrp2_add[0] += lcoef1[q]; //bmnc
							ltrp1[q] = NULL;
							lcoef1[q] = 0.0;
						}
					}
				}
				ltrp1_add[1] -= deny;
				//flux from C_f to C_K
				deny = 0.0;
				for (int q = 0; q < 3; q++)
				{
				
					if (fabs(lcoef2[q] / ltrp2_add[1]) < 1.0e-8) //contribution is less then 0.5%
					{
						ltrp2[q] = NULL;
						deny += lcoef2[q];
						lcoef2[q] = 0.0;
					}
					else
					{
						if (ltrp2[q] == r0)
						{
							ltrp1_add[0] += lcoef2[q]; //bmnc
							ltrp2[q] = NULL;
							lcoef2[q] = 0.0;
						}
					}
				}
				ltrp2_add[1] -= deny;
				//Flux from C_L to C_f
				deny = 0.0;
				for (int q = 0; q < 3; q++)
				{
				
					if (fabs(rcoef1[q] / rtrp1_add[1]) < 1.0e-8) //contribution is less then 0.5%
					{
						rtrp1[q] = NULL;
						deny += rcoef1[q];
						rcoef1[q] = 0.0;
					}
					else //TODO
					{
						if (rtrp1[q] == &*f)
						{
							rtrp2_add[0] += rcoef1[q]; //bmnc
							rtrp1[q] = NULL;
							rcoef1[q] = 0.0;
						}
					}
				}
				rtrp1_add[1] -= deny;
				//flux from C_f to C_L
				deny = 0.0;
				for (int q = 0; q < 3; q++)
				{
				
					if (fabs(rcoef2[q] / rtrp2_add[1]) < 1.0e-8) //contribution is less then 0.5%
					{
						rtrp2[q] = NULL;
						deny += rcoef2[q];
						rcoef2[q] = 0.0;
					}
					else //TODO
					{
						if (rtrp2[q] == r1)
						{
							rtrp1_add[0] += rcoef2[q]; //bmnc
							rtrp2[q] = NULL;
							rcoef2[q] = 0.0;
						}
					}
				}
				rtrp2_add[1] -= deny;
			
			}
			else
			{
				std::cout << "Failed to find triplet!" << std::endl;
				error = 1;
			}
		}
		else
		{
			int ret1, ret2, ptypes[4][3];
			ret1 = tensor_prod(Ktype, r0->RealArray(K), nrmf, f1);
			
			if (!ret1)	error = 1; //error in tensor

			Storage::reference_array ltrp = f->ReferenceArrayDF(trpl1_elems);
			Storage::reference_array rtrp = f->ReferenceArrayDF(trpl2_elems);
			Storage::real_array lcoef = f->RealArrayDF(trpl1_coefs);
			Storage::real_array rcoef = f->RealArrayDF(trpl2_coefs);
			Storage::real_array ltrp_add = f->RealArrayDF(trpl1_add);
			Storage::real_array rtrp_add = f->RealArrayDF(trpl2_add);

			area1 = area*normalize(f1);

		
			ret1 = find_triplet(r0, f, NULL, f1, +1, area1, &ltrp[0], &lcoef[0],ptypes[0],(allowed_types & (~AVG_NEUMANN)),true);
			ret2 = find_triplet(f, r0, NULL, f1, -1, area1, &rtrp[0], &rcoef[0],ptypes[2],(allowed_types & (~AVG_NEUMANN)),true);

			if (ret1 && ret2)
			{
				ltrp_add[0] = ltrp_add[1] =  0.0;
				rtrp_add[0] = rtrp_add[1] = 0.0;

				for (int q = 0; q < 3; q++) 
				{
				
					ltrp_add[1] += lcoef[q];
					rtrp_add[1] += rcoef[q];
				
				}
				//Flux from C_K to C_L
				deny = 0.0;
				for (int q = 0; q < 3; q++)
				{
				
					if (fabs(lcoef[q] / ltrp_add[1]) < 1.0e-8) //contribution is less then 0.5%
					{
						ltrp[q] = NULL;
						deny += lcoef[q];
						lcoef[q] = 0.0;
					}
					else
					{
						if (ltrp[q] == f)
						{
							rtrp_add[0] += lcoef[q]; //bmnc
							ltrp[q] = NULL;
							lcoef[q] = 0.0;
						}
					}
				}
				ltrp_add[1] -= deny;
				//flux from C_L to C_K
				deny = 0.0;
				for (int q = 0; q < 3; q++)
				{
				
					if (fabs(rcoef[q] / rtrp_add[1]) < 1.0e-8) //contribution is less then 0.5%
					{
						rtrp[q] = NULL;
						deny += rcoef[q];
						rcoef[q] = 0.0;
					}
					else
					{
						if (rtrp[q] == r0)
						{
							ltrp_add[0] += rcoef[q]; //bmpc
							rtrp[q] = NULL;
							rcoef[q] = 0.0;
						}
					}
				}
				rtrp_add[1] -= deny;
			}
			else
			{
				std::cout << "Failed to find triplet!" << std::endl;
				error = 1;
			}

		}
	}
}




expr ntpfa_c::Grad(const expr & param) const
{
	expr dp1 = 1.0e-9 + stencil(trpl1,param);
	expr dn1 = 1.0e-9 + stencil(trpl2,param); 
	expr dp2 = 1.0e-9 + stencil(trpr1,param);
	expr dn2 = 1.0e-9 + stencil(trpr2,param); 
	expr bmpc1 = tagval(addl1,0);
	expr tkp1 = tagval(addl1,1);
	expr bmnc1 = tagval(addl2,0);
	expr tkn1 = tagval(addl2,1);
	expr bmpc2 = tagval(addr1,0);
	expr tkp2 = tagval(addr1,1);
	expr bmnc2 = tagval(addr2,0);
	expr tkn2 = tagval(addr2,1);
	expr mp1 = dn1 / (dp1 + dn1);
	expr mn1 = dp1 / (dp1 + dn1);
	expr mp2 = dn2 / (dp2 + dn2);
	expr mn2 = dp2 / (dp2 + dn2);
	expr np = (mp1 * tkp1 + mn1 * bmpc1);
	expr npf = (mn1 * tkn1 + mp1 * bmnc1);
	expr nn = (mp2 * tkp2 + mn2 * bmpc2);
	expr nnf = (mn2 * tkn2 + mp2 * bmnc2);
	return (stencil(upw,param)*nn*npf-stencil(dwn,param)*np*nnf)/(npf+nnf);
}

expr ntpfa_c::BoundaryGrad(const expr & param) const
{
	expr dp = 1.0e-9 + stencil(trpl1, param);
	expr dn = 1.0e-9 + stencil(trpl2, param);
	expr bmpc = tagval(addl1, 0);
	expr tkp = tagval(addl1, 1);
	expr bmnc = tagval(addl2, 0);
	expr tkn = tagval(addl2, 1);
	expr mp = dn / (dn + dp);
	expr mn = dp / (dn + dp);
	expr bmp = (mp * tkp + mn * bmpc);
	expr bmn = (mn * tkn + mp * bmnc);
	return stencil(upw, param)*bmn-stencil(dwn,param)*bmp;
}

expr ntpfa_c::Interp(ElementType etype,const expr & param) const
{
	if( etype == FACE )
	{
		expr dp1 = 1.0e-9 + stencil(trpl1,param);
		expr dn1 = 1.0e-9 + stencil(trpl2,param); 
		expr dp2 = 1.0e-9 + stencil(trpr1,param);
		expr dn2 = 1.0e-9 + stencil(trpr2,param); 
		expr bmpc1 = tagval(addl1,0);
		expr tkp1 = tagval(addl1,1);
		expr bmnc1 = tagval(addl2,0);
		expr tkn1 = tagval(addl2,1);
		expr bmpc2 = tagval(addr1,0);
		expr tkp2 = tagval(addr1,1);
		expr bmnc2 = tagval(addr2,0);
		expr tkn2 = tagval(addr2,1);
		
		expr mp1 = dn1 / (dp1 + dn1);
		expr mn1 = dp1 / (dp1 + dn1);
		expr mp2 = dn2 / (dp2 + dn2);
		expr mn2 = dp2 / (dp2 + dn2);
		expr np = (mp1 * tkp1 + mn1 * bmpc1);
		expr npf = (mn1 * tkn1 + mp1 * bmnc1);
		expr nn = (mp2 * tkp2 + mn2 * bmpc2);
		expr nnf = (mn2 * tkn2 + mp2 * bmnc2);
		return (stencil(dwn,param)*np + stencil(upw,param)*nn)/(npf+nnf);
	}
	else return stencil(intrp,param);
}

void ntpfa_c::Export(std::ostream & fcout, Storage::real trans_scale, Storage::real vol_scale) const
{
	discr_basic::Export(fcout,trans_scale,vol_scale);
	Tag unknown_id = m->CreateTag("ID",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	Tag row = m->CreateTag("ROW",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	int cnt = 0, cntsupp = 0;
	INMOST_DATA_ENUM_TYPE idnum = 0;
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		it->IntegerDF(unknown_id) = idnum++;
	fcout << "-- block unknowns# " << idnum << std::endl;
	for(Mesh::iteratorElement it = m->BeginElement(FACE | EDGE | NODE); it != m->EndElement(); ++it)
		if(it->GetMarker(add_markers)) it->IntegerDF(unknown_id) = idnum+cntsupp++;
	fcout << "-- total unknowns# " << idnum+cntsupp << std::endl;
	fcout << std::setprecision(16);
	
	fcout << "NTPFACONNS" << std::endl;
	for (Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
	{
		if (needBuild(&*f)) 
		{
			if( f->GetMarker(bnd_markers) && !f->GetMarker(add_markers) ) continue;
			cnt++;
		}
	}
	fcout << cnt << std::endl;
	cnt = 0;
	MIDType marker = m->CreateMarker(), mrkunique = m->CreateMarker();
	dynarray<Element *,128> stack;
	Tag num = m->CreateTag("TEMP_POS",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	Tag coe = m->CreateTag("TEMP_COE",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,4);
	
	for (Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
	{
		if (needBuild(&*f))
		{
			//fcout << "-- face# " << f->LocalID() << " conn# " << cnt;
			//if( f->GetMarker(bnd_markers) ) fcout << " is boundary connection ";
			//fcout << std::endl;

			int end = 4;
			Element * r0 = f->BackCell(), *r1 = f->FrontCell();
			if( f->GetMarker(bnd_markers) )
			{
				if( !f->GetMarker(add_markers) ) continue;
				r1 = &*f;
				end = 2;
			}
			int abnP[4] = {0,0,0,0};
			stack.clear();
			stack.push_back(r0);
			stack.push_back(r1);
			r0->SetMarker(marker);
			r1->SetMarker(marker);
			r0->IntegerDF(num) = 0;
			r1->IntegerDF(num) = 1;
			Storage::reference_array lrtrp[4];
			lrtrp[0] = f->ReferenceArrayDF(trpl1_elems);
			lrtrp[1] = f->ReferenceArrayDF(trpl2_elems);
			lrtrp[2] = f->ReferenceArrayDF(trpr1_elems);
			lrtrp[3] = f->ReferenceArrayDF(trpr2_elems);
			Storage::real_array lrcoef[4];
			lrcoef[0] = f->RealArrayDF(trpl1_coefs);
			lrcoef[1] = f->RealArrayDF(trpl2_coefs);
			lrcoef[2] = f->RealArrayDF(trpr1_coefs);
			lrcoef[3] = f->RealArrayDF(trpr2_coefs);
			Storage::real_array ltrp1_add = f->RealArrayDF(trpl1_add);
			Storage::real_array ltrp2_add = f->RealArrayDF(trpl2_add);
			Storage::real_array rtrp1_add = f->RealArrayDF(trpr1_add);
			Storage::real_array rtrp2_add = f->RealArrayDF(trpr2_add);

			for(int l = 0; l < end; l++)
			{
				for (int q = 0; q < 3; q++) if (lrtrp[l][q] != NULL) 
				{
					if(!lrtrp[l][q]->GetMarker(mrkunique)) 
					{
						abnP[l]++;
						if(!lrtrp[l][q]->GetMarker(marker))
						{
							lrtrp[l][q]->IntegerDF(num) = stack.size();
							stack.push_back(lrtrp[l][q]);
							lrtrp[l][q]->SetMarker(marker);
						}
						lrtrp[l][q]->RealArrayDF(coe)[l] = lrcoef[l][q];
						lrtrp[l][q]->SetMarker(mrkunique);
					}
					else lrtrp[l][q]->RealArrayDF(coe)[l] += lrcoef[l][q];
				}
				for (int q = 0; q < 3; q++) if (lrtrp[l][q] != NULL) 
					lrtrp[l][q]->RemMarker(mrkunique);
			}


			
			fcout << stack.size() << " ";
			//fcout << std::endl;

			for (int k = 0; k < stack.size(); ++k) 
			{
				fcout << stack[k]->IntegerDF(unknown_id) << " ";
				stack[k]->RemMarker(marker);
			}

			fcout << 1+(end==2) << " "; //nonlinear type

			//fcout << std::endl;
			fcout << abnP[0] << " "; //number of elements that construct left coefficient
			fcout << abnP[1] << " "; //number of elements that construct right coefficient
			if( end > 2 )
			{
				fcout << abnP[2] << " "; //number of elements that construct left coefficient
				fcout << abnP[3] << " "; //number of elements that construct right coefficient
			}

			

			std::string words[4] = {"-- r0 to f","-- f to r0","-- r1 to f","-- f to r1"};

			
			//fcout << std::endl;
			fcout << std::scientific << ltrp1_add[0] * trans_scale << " "; //bmpc1
			fcout << std::scientific << ltrp1_add[1] * trans_scale << " "; //tkp1
			//fcout << "-- r0 to f";
			//fcout << std::endl;
			fcout << std::scientific << ltrp2_add[0] * trans_scale << " "; //bmnc1
			fcout << std::scientific << ltrp2_add[1] * trans_scale << " "; //tkn1
			//fcout << "-- f to r0";
			//fcout << std::endl;
			if( end > 2 )
			{
				fcout << std::scientific << rtrp1_add[0] * trans_scale << " "; //bmpc2
				fcout << std::scientific << rtrp1_add[1] * trans_scale << " "; //tkp2
				//fcout << "-- r1 to f";
				//fcout << std::endl;
				fcout << std::scientific << rtrp2_add[0] * trans_scale << " "; //bmnc2
				fcout << std::scientific << rtrp2_add[1] * trans_scale << " "; //tkn2
				//fcout << "-- f to r1";
				//fcout << std::endl;
			}
			
			for(int l = 0; l < end; l++)
			{
				for (int q = 0; q < 3; q++) if (lrtrp[l][q] != NULL && !lrtrp[l][q]->GetMarker(mrkunique)) 
				{
					fcout << lrtrp[l][q]->IntegerDF(num) << " " << std::scientific << lrtrp[l][q]->RealArrayDF(coe)[l]*trans_scale << " ";
					lrtrp[l][q]->SetMarker(mrkunique);
				}
				for (int q = 0; q < 3; q++) if (lrtrp[l][q] != NULL) lrtrp[l][q]->RemMarker(mrkunique);
				//fcout << words[l];
				//fcout << std::endl;
			}
			
			fcout << std::endl;
			f->IntegerDF(row) = cnt++;
		}
	}
	fcout << "/" << std::endl;
	if( cntsupp )
	{
		fcout << "SUPPCONNS" << std::endl;
		fcout << cntsupp << std::endl;
		for (Mesh::iteratorElement f = m->BeginElement(EDGE|NODE); f != m->EndElement(); ++f)
		{
			if( f->GetMarker(add_markers) )
			{
				//fcout << "--" << ElementTypeName(f->GetElementType()) << "# " << f->LocalID() << " supp# " << f->IntegerDF(unknown_id) << " ";
				//Storage::real cnt[3];
				//f->Centroid(cnt);
				//fcout << cnt[0] << " " << cnt[1] << " " << cnt[2] << std::endl;

				fcout << 0 << " "; //linear type
				Storage::reference_array stncl = f->ReferenceArrayDV(intrp_stncl);
				Storage::real_array coefs = f->RealArrayDV(intrp_coefs);
				fcout << stncl.size() << " ";
				for(int k = 0; k < stncl.size(); k++)
				{
					fcout << stncl[k]->IntegerDF(unknown_id) << " ";
				}
				for(int k = 0; k < stncl.size(); k++)
				{
					fcout << std::scientific << coefs[k] << " "; //revdist should not be scaled
				}
				fcout << std::endl;
			}
		}
		int recs = 0;
		for (Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
		{
			if( f->GetMarker(add_markers) )
			{
				if( !f->GetMarker(bnd_markers) )
				{
					//fcout << "--supp " << ElementTypeName(f->GetElementType()) << "# " << f->LocalID() << " supp# " << f->IntegerDF(unknown_id)<< std::endl;

					Cell * r0 = f->BackCell(), *r1 = f->FrontCell();
			
					int abnP[4] = {0,0,0,0};
					stack.clear();
					stack.push_back(r0);
					stack.push_back(r1);
					r0->SetMarker(marker);
					r1->SetMarker(marker);
					r0->IntegerDF(num) = 0;
					r1->IntegerDF(num) = 1;
					Storage::reference_array lrtrp[4];
					lrtrp[0] = f->ReferenceArrayDF(trpl1_elems);
					lrtrp[1] = f->ReferenceArrayDF(trpl2_elems);
					lrtrp[2] = f->ReferenceArrayDF(trpr1_elems);
					lrtrp[3] = f->ReferenceArrayDF(trpr2_elems);
					Storage::real_array lrcoef[4];
					lrcoef[0] = f->RealArrayDF(trpl1_coefs);
					lrcoef[1] = f->RealArrayDF(trpl2_coefs);
					lrcoef[2] = f->RealArrayDF(trpr1_coefs);
					lrcoef[3] = f->RealArrayDF(trpr2_coefs);
					Storage::real_array ltrp1_add = f->RealArrayDF(trpl1_add);
					Storage::real_array ltrp2_add = f->RealArrayDF(trpl2_add);
					Storage::real_array rtrp1_add = f->RealArrayDF(trpr1_add);
					Storage::real_array rtrp2_add = f->RealArrayDF(trpr2_add);

					for(int l = 0; l < 4; l++)
					{
						for (int q = 0; q < 3; q++) if (lrtrp[l][q] != NULL) 
						{
							if(!lrtrp[l][q]->GetMarker(mrkunique)) 
							{
								abnP[l]++;
								if(!lrtrp[l][q]->GetMarker(marker))
								{
									lrtrp[l][q]->IntegerDF(num) = stack.size();
									stack.push_back(lrtrp[l][q]);
									lrtrp[l][q]->SetMarker(marker);
								}
								lrtrp[l][q]->RealArrayDF(coe)[l] = lrcoef[l][q];
								lrtrp[l][q]->SetMarker(mrkunique);
							}
							else lrtrp[l][q]->RealArrayDF(coe)[l] += lrcoef[l][q];
						}
						for (int q = 0; q < 3; q++) if (lrtrp[l][q] != NULL) 
							lrtrp[l][q]->RemMarker(mrkunique);
					}


					//fcout << "-- face# " << f->LocalID() << std::endl;

					fcout << 1 << " "; //nonlinear type
					/*
					Storage::real cnt[3];
					f->Centroid(cnt);
					fcout << cnt[2] << " "; // depth
					*/
					fcout << stack.size() << " ";
					//fcout << std::endl;

					for (int k = 0; k < stack.size(); ++k) 
					{
						fcout << stack[k]->IntegerDF(unknown_id) << " ";
						stack[k]->RemMarker(marker);
					}

				

					//fcout << std::endl;
					fcout << abnP[0] << " "; //number of elements that construct left coefficient
					fcout << abnP[1] << " "; //number of elements that construct right coefficient
					fcout << abnP[2] << " "; //number of elements that construct left coefficient
					fcout << abnP[3] << " "; //number of elements that construct right coefficient

			

					std::string words[4] = {"-- r0 to f","-- f to r0","-- r1 to f","-- f to r1"};

			
					//fcout << std::endl;
					fcout << std::scientific << ltrp1_add[0] * trans_scale << " "; //bmpc1
					fcout << std::scientific << ltrp1_add[1] * trans_scale << " "; //tkp1
					//fcout << "-- r0 to f";
					//fcout << std::endl;
					fcout << std::scientific << ltrp2_add[0] * trans_scale << " "; //bmnc1
					fcout << std::scientific << ltrp2_add[1] * trans_scale << " "; //tkn1
					//fcout << "-- f to r0";
					//fcout << std::endl;
					fcout << std::scientific << rtrp1_add[0] * trans_scale << " "; //bmpc2
					fcout << std::scientific << rtrp1_add[1] * trans_scale << " "; //tkp2
					//fcout << "-- r1 to f";
					//fcout << std::endl;
					fcout << std::scientific << rtrp2_add[0] * trans_scale << " "; //bmnc2
					fcout << std::scientific << rtrp2_add[1] * trans_scale << " "; //tkn2
					//fcout << "-- f to r1";
					//fcout << std::endl;
			
					for(int l = 0; l < 4; l++)
					{
						for (int q = 0; q < 3; q++) if (lrtrp[l][q] != NULL && !lrtrp[l][q]->GetMarker(mrkunique)) 
						{
							fcout << lrtrp[l][q]->IntegerDF(num) << " " << std::scientific << lrtrp[l][q]->RealArrayDF(coe)[l]*trans_scale << " ";
							lrtrp[l][q]->SetMarker(mrkunique);
						}
						for (int q = 0; q < 3; q++) if (lrtrp[l][q] != NULL) lrtrp[l][q]->RemMarker(mrkunique);
						//fcout << words[l];
						//fcout << std::endl;
					}
			
					
				}
				else
				{
					//fcout << "--bnd " << ElementTypeName(f->GetElementType()) << "# " << f->LocalID() << " supp# " << f->IntegerDF(unknown_id) << " ";
					//Storage::real cnt[3];
					//f->Centroid(cnt);
					//fcout << cnt[0] << " " << cnt[1] << " " << cnt[2] << std::endl;
					//fcout << std::endl;
					fcout << 2 << " "; //boundary type
					fcout << f->IntegerDF(row) << " "; //corresponding entry in connection list
					// \alpha C + \beta dC/dn = \gamma
					if( bnd_conds.isValid() && f->HaveData(bnd_conds) )
					{
						Storage::real_array bnd = f->RealArray(bnd_conds);
						fcout << bnd[0] << " "; // \alpha
						fcout << bnd[1] << " "; // \beta
						fcout << bnd[2] << " "; // \gamma
					}
					else
					{
						fcout << 0.0 << " "; // \alpha
						fcout << 1.0 << " "; // \beta
						fcout << 0.0 << " "; // \gamma
					}
				}
				fcout << std::endl;
				recs++;
			}
		}
		fcout << "/" << std::endl;
	}
	m->DeleteTag(row);
	m->DeleteTag(coe);
	m->DeleteTag(num);
	m->ReleaseMarker(marker);
	m->ReleaseMarker(mrkunique);
	//fcout << "/" << std::endl;
}


void ntpfa_c::Update() //TODO
{
	//remember the stencil size, then analyze cells around, if some change - recompute
	for (INMOST_DATA_INTEGER_TYPE id = 0; id < m->MaxLocalID(FACE); ++id)
	{
		Face * f = m->FaceByLocalID(id);
		if (f != NULL)
		{
			f->RemMarker(add_markers);
			if (needUpdate(f))
			{
				Build(f); //always recompute
			}
		}
	}
}

ntpfa_c::~ntpfa_c()
{
	m->DeleteTag(intrp_stncl);
	m->DeleteTag(intrp_coefs);
	m->DeleteTag(trpl1_elems);
	m->DeleteTag(trpl2_elems);
	m->DeleteTag(trpr1_elems);
	m->DeleteTag(trpr2_elems);
	m->DeleteTag(trpl1_coefs);
	m->DeleteTag(trpl2_coefs);
	m->DeleteTag(trpr1_coefs);
	m->DeleteTag(trpr2_coefs);
	m->DeleteTag(trpl1_add);
	m->DeleteTag(trpl2_add);
	m->DeleteTag(trpr1_add);
	m->DeleteTag(trpr2_add);
}