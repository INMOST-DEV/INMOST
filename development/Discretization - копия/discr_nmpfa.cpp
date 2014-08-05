#include "discr.h"
#include <iomanip>
#define EPS 1e-12

typedef Storage::real vec[3];

static Storage::real lengthS(vec x) { return dot_prod(x, x); }
static Storage::real length(vec x) { return sqrt(lengthS(x)); }
static Storage::real normalize(Storage::real x[3]) { Storage::real r = length(x); if (r > 0.0)  x[0] /= r, x[1] /= r, x[2] /= r; return r; }

nmpfa::nmpfa(Automatizator * aut, Mesh * m, std::string tensor)
:discr_triplet_basic(aut, m, tensor)
{
	trpl_elems = m->CreateTag("NMPFA_TRP_L", DATA_REFERENCE, FACE, NONE, 4);
	trpr_elems = m->CreateTag("NMPFA_TRP_R", DATA_REFERENCE, FACE, NONE, 4);
	trpl_coefs = m->CreateTag("NMPFA_COEF_L", DATA_REAL, FACE, NONE, 4);
	trpr_coefs = m->CreateTag("NMPFA_COEF_R", DATA_REAL, FACE, NONE, 4);
	addval = m->CreateTag("NMPFA_ADD", DATA_REAL, FACE, NONE, 2);

	
	trpr = aut->RegisterStencil("nmpfa_right_triplet", trpr_elems, trpr_coefs);
	trpl = aut->RegisterStencil("nmpfa_left_triplet", trpl_elems, trpl_coefs);
	add = aut->RegisterStaticTag(addval);
	intrp_stncl = m->CreateTag("INTRP_STNCL", DATA_REFERENCE, FACE | EDGE | NODE, NONE);
	intrp_coefs = m->CreateTag("INTRP_COEFS", DATA_REAL, FACE | EDGE | NODE, NONE);
	
	
	dwn = aut->RegisterStencil("back_cell", get_l);
	upw = aut->RegisterStencil("front_cell", get_r);

	intrp = aut->RegisterStencil("interp",intrp_stncl,intrp_coefs,add_markers);
}


void nmpfa::MarkRecursive(Element * f)
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
			MarkRecursive(*it);
		}
	}
}


void nmpfa::Init()
{
	if (!have_bnd) Boundary(0, FACE | EDGE);
	for (INMOST_DATA_INTEGER_TYPE id = 0; id < m->MaxLocalID(FACE); ++id)
	{
		Face * f = m->FaceByLocalID(id);
		if (f != NULL) 
		{
			Build(f);
		}
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
	for (Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it)
	{
		if( it->GetMarker(add_markers) )
			MarkRecursive(&*it);
	}
	m->ReleaseMarker(mark_next);

	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) it->SetMarker(add_markers);
	
	m->SynchronizeMarker(UnknownMarker(),FACE|EDGE|NODE,SYNC_BIT_OR);
}


void nmpfa::get_l(Storage * e, Automatizator::stencil_pairs & out, void * user_data)
{
	Face * f = static_cast<Face *>(e);
	Element * r = f->BackCell();
	out.push_back(Automatizator::stencil_pair(r == NULL ? f : r, 1.0));
}
void nmpfa::get_r(Storage * e, Automatizator::stencil_pairs & out, void * user_data)
{
	Face * f = static_cast<Face *>(e);
	Element * r = f->FrontCell();
	out.push_back(Automatizator::stencil_pair(r == NULL ? f : r, 1.0));
}


void nmpfa::Build(Face * f)
{
	if (needBuild(f))
	{
		Storage::real nrmf[3], area, area1, area2, f1[3], f2[3], cnt1[3], cnt2[3], cntf[3], deny;
		Element * r0 = f->BackCell();
		Element * r1 = f->FrontCell();
		r0->Centroid(cnt1);
		(r1 == NULL ? f : r1)->Centroid(cnt2);
		f->Centroid(cntf);
		f->UnitNormal(nrmf);
		area = f->Area();

		if( !f->GetMarker(bnd_markers) )
		{
			int ret1, ret2, ptypes[4][3];
			ret1 = tensor_prod(Ktype, r0->RealArray(K), nrmf, f1);
			ret2 = tensor_prod(Ktype, r1->RealArray(K), nrmf, f2);

			if (!ret1 || !ret2)	error = 1; //error in tensor

			Storage::reference_array ltrp = f->ReferenceArrayDF(trpl_elems);
			Storage::reference_array rtrp = f->ReferenceArrayDF(trpr_elems);
			Storage::real_array lcoef = f->RealArrayDF(trpl_coefs);
			Storage::real_array rcoef = f->RealArrayDF(trpr_coefs);
			Storage::real_array addv = f->RealArrayDF(addval);
			area1 = area*normalize(f1);
			area2 = area*normalize(f2);

		
			ret1 = find_triplet(r0, NULL, f, f1, +1, area1, &ltrp[0], &lcoef[0],ptypes[0],allowed_types);
			ret2 = find_triplet(r1, NULL, f, f2, -1, area2, &rtrp[0], &rcoef[0],ptypes[2],allowed_types);


			if (ret1 && ret2)
			{
				Storage::real lsum = 0, rsum = 0;
				for (int q = 0; q < 3; q++) 
				{
					lsum += lcoef[q];
					rsum += rcoef[q];
				}
				
				//Flux from C_K to C_L
				deny = 0.0;
				for (int q = 0; q < 3; q++)
				{
				
					if (fabs(lcoef[q] / lsum) < 1.0e-8) //contribution is less then 0.5%
					{
						ltrp[q] = NULL;
						deny += lcoef[q];
						lcoef[q] = 0.0;
					}
					else
					{
						if (ltrp[q] == r1)
						{
							addv[0] += lcoef[q]; //bmpc
							ltrp[q] = NULL;
							deny += lcoef[q];
							lcoef[q] = 0.0;
							
						}
					}
				}
				ltrp[3] = r0;
				lcoef[3] = -(lsum-deny);
				//flux from C_L to C_K
				deny = 0.0;
				for (int q = 0; q < 3; q++)
				{
				
					if (fabs(rcoef[q] / rsum) < 1.0e-8) //contribution is less then 0.5%
					{
						rtrp[q] = NULL;
						deny += rcoef[q];
						rcoef[q] = 0.0;
					}
					else
					{
						if (rtrp[q] == r0)
						{
							addv[1] += rcoef[q]; //bmnc
							rtrp[q] = NULL;
							deny += rcoef[q];
							rcoef[q] = 0.0;
							
						}
					}
					rcoef[q] = -rcoef[q];
				}

				rtrp[3] = r1;
				rcoef[3] = (rsum-deny);

				for (int q = 0; q < 3; q++) 
				{
					if( ltrp[q] != NULL && ltrp[q]->GetElementType() != CELL ) ltrp[q]->SetMarker(add_markers);
					if( rtrp[q] != NULL && rtrp[q]->GetElementType() != CELL ) rtrp[q]->SetMarker(add_markers);
				}
			}
			else
			{
				std::cout << "Failed to find triplet!" << std::endl;
				error = 1;
			}
		

			Element * t[2][3];
			Storage::real c[2][3];
			ret1 = find_triplet(f, r0, NULL, f1, -1, area1, t[0], c[0],ptypes[0],allowed_types);
			ret2 = find_triplet(f, r1, NULL, f2, +1, area2, t[1], c[1],ptypes[1],allowed_types);

			if( ret1 && ret2 )
			{
				Storage::real sum = 0;
				Storage::reference_array stncl = f->ReferenceArrayDV(intrp_stncl);
				Storage::real_array coefs = f->RealArrayDV(intrp_coefs);

				for(int q = 0; q < 3; q++)
				{
					sum += c[0][q];
					sum += c[1][q];
				}

				for(int j = 0; j < 2; j++)
				for(int q = 0; q < 3; q++)
				{
					if(c[j][q]/sum > 1.0e-13 )
					{
						//if( ptypes[j][q] == AVG_NONLINEAR ) t[j][q]->SetMarker(add_markers);
						bool found = false;
						for(int k = 0; k < stncl.size(); k++)
						if( stncl[k] == t[j][q] )
						{
							coefs[k] += c[j][q]/sum;
							found = true;
							break;
						}
						if( !found )
						{
							stncl.push_back(t[j][q]);
							coefs.push_back(c[j][q]/sum);
						}
					}
				}

			}
			else
			{
				std::cout << "Failed to find triplet for interpolation" << std::endl;
				error = 1;
			}
		}

		else
		{
			/*
			int ptypes[2][3];
			int ret = tensor_prod(Ktype, r0->RealArray(K), nrmf, f1);
			if( ret )
			{
				Storage::reference_array ltrp = f->ReferenceArrayDF(trpl_elems);
				Storage::real_array lcoef = f->RealArrayDF(trpl_coefs);
				area1 = area*normalize(f1);

				ret = find_triplet(r0, f, NULL, f1, +1, area1, &ltrp[0], &lcoef[0],ptypes[0],allowed_types & ~AVG_NEUMANN, true);
				if( ret )
				{
					bool found_f = false;
					for (int q = 0; q < 3; q++) 
					{
						if( fabs(lcoef[q]) < 1.0e-9 ) 
						{
							ltrp[q] = NULL;
							lcoef[q] = 0.0;
						}

						if( ltrp[q] == f ) found_f = true;
					}

					if( !found_f ) 
					{
						std::cout << "face not found!" << std::endl;
					}
				}
				else
				{
					std::cout << "Failed to find triplet for interpolation" << std::endl;
					error = 1;
				}
			}
			*/

			int ptypes[3];
			int ret = tensor_prod(Ktype, r0->RealArray(K), nrmf, f1);
			if( ret )
			{
				Element * t[3];
				Storage::real c[3];
				Storage::reference_array astncl = f->ReferenceArrayDV(intrp_stncl);
				Storage::real_array acoefs = f->RealArrayDV(intrp_coefs);
				area1 = area*normalize(f1);

				astncl.push_back(r0);
				acoefs.push_back(0.0);
				astncl.push_back(f);
				acoefs.push_back(0.0);

				ret = find_triplet(r0, f, NULL, f1, 1, area1,t, c,ptypes,(allowed_types & (~AVG_NEUMANN)));
				if( ret )
				{
					for (int q = 0; q < 3; q++) 
					{
						if( fabs(c[q]) > 1.0e-9 ) 
						{
							acoefs[0] -= c[q]*0.5;
							bool found = false;
							for(int k = 0; k < astncl.size(); k++)
							{
								if(astncl[k] == t[q])
								{
									acoefs[k] += c[q]*0.5;
									found = true;
									break;
								}
							}
							if( !found )
							{
								astncl.push_back(t[q]);
								acoefs.push_back(c[q]*0.5);
							}
						}
					}
				}
				else
				{
					std::cout << "Failed to find triplet for interpolation" << std::endl;
					error = 1;
				}

				ret = find_triplet(f, r0, NULL, f1, -1, area1,t, c,ptypes,(allowed_types & (~AVG_NEUMANN)));
				if( ret )
				{
					for (int q = 0; q < 3; q++) 
					{
						if( fabs(c[q]) > 1.0e-9 ) 
						{
							acoefs[1] += c[q]*0.5;
							bool found = false;
							for(int k = 0; k < astncl.size(); k++)
							{
								if(astncl[k] == t[q])
								{
									acoefs[k] -= c[q]*0.5;
									found = true;
									break;
								}
							}
							if( !found )
							{
								astncl.push_back(t[q]);
								acoefs.push_back(-c[q]*0.5);
							}
						}
					}
				}
				else
				{
					std::cout << "Failed to find triplet for interpolation" << std::endl;
					error = 1;
				}
			}

		}
	}
}




expr nmpfa::Grad(const expr & param) const
{
	expr cp = tagval(add,0);
	expr cn = tagval(add,1);
	expr dp = stencil(trpl, param);
	expr dn = stencil(trpr, param);
	//expr dp2 = ad_abs(dp) + 1.0e-7;
	//expr dn2 = ad_abs(dn) + 1.0e-7;
	expr dp2 = dp*dp / ad_sqrt(dp*dp + 1.0e-5) + 1.0e-7;
	expr dn2 = dn*dn / ad_sqrt(dn*dn + 1.0e-5) + 1.0e-7;
	//expr mp = ad_val(dn2 / (dn2 + dp2));
	//expr mn = ad_val(dp2 / (dn2 + dp2));
	expr mp = dn2 / (dn2 + dp2);
	expr mn = dp2 / (dn2 + dp2);
	//expr mp2 = condition(dn2 - 1.0e-1,mp,ad_val(mp));
	//expr np2 = condition(dp2 - 1.0e-1,mn,ad_val(mn));

	return (mp*cp + mn*cn)*(stencil(upw, param)-stencil(dwn,param)) + mp*dp+mn*dn;
	//return (dp+dn)*0.5 + (cp+cn)*0.5*(stencil(upw, param)-stencil(dwn,param));
}

expr nmpfa::Interp(ElementType etype,const expr & param) const
{
	return stencil(intrp,param);
}

void nmpfa::Export(std::ostream & fcout, Storage::real trans_scale, Storage::real vol_scale) const
{
	discr_basic::Export(fcout,trans_scale,vol_scale);
	Tag unknown_id = m->CreateTag("ID",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	Tag row = m->CreateTag("ROW",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	int cnt = 0, cntsupp = 0;
	INMOST_DATA_ENUM_TYPE idnum = 0;
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		it->IntegerDF(unknown_id) = idnum++;
	for(Mesh::iteratorElement it = m->BeginElement(FACE | EDGE | NODE); it != m->EndElement(); ++it)
		if(it->GetMarker(add_markers)) it->IntegerDF(unknown_id) = idnum+cntsupp++;

	fcout << std::setprecision(16);
	fcout << "NMPFACONNS" << std::endl;
	
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
	for (Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
		if( f->GetMarker(marker) ) std::cout << "FACE " << f->LocalID() << std::endl;


	dynarray<Element *,128> stack;
	Tag num = m->CreateTag("TEMP_POS",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	Tag coe = m->CreateTag("TEMP_COE",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,4);
	
	for (Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
	{
		if (needBuild(&*f))
		{

			
			
			Cell * r0 = f->BackCell(), *r1 = f->FrontCell();
			
			int abnP[2] = {0,0};
			
			

			if( !f->GetMarker(bnd_markers) )
			{
				//fcout << "-- face# " << f->LocalID() << " conn# " << cnt << std::endl;
				
				stack.clear();
				stack.push_back(r0);
				r0->SetMarker(marker);
				r0->IntegerDF(num) = 0;
				stack.push_back(r1);
				r1->SetMarker(marker);
				r1->IntegerDF(num) = 1;
				Storage::reference_array lrtrp[4];
				lrtrp[0] = f->ReferenceArrayDF(trpl_elems);
				lrtrp[1] = f->ReferenceArrayDF(trpr_elems);
				Storage::real_array lrcoef[4];
				lrcoef[0] = f->RealArrayDF(trpl_coefs);
				lrcoef[1] = f->RealArrayDF(trpr_coefs);
				Storage::real_array addv = f->RealArrayDF(addval);

				for(int l = 0; l < 2; l++)
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

				fcout << 3 << " "; //nonlinear type

				//fcout << std::endl;
				fcout << abnP[0] << " "; //number of elements that construct left coefficient
				fcout << abnP[1] << " "; //number of elements that construct right coefficient
			
				std::string words[2] = {"-- r0 to r1","-- r1 to r0"};

			
				//fcout << std::endl;
				fcout << std::scientific << addv[0] * trans_scale << " "; //for mp
				fcout << std::scientific << addv[1] * trans_scale << " "; //for mn
				//fcout << "-- r0 to r1";
				//fcout << std::endl;
			
				for(int l = 0; l < 2; l++)
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
				if( !f->GetMarker(add_markers) ) continue;
				/*
				stack.clear();
				stack.push_back(r0);
				r0->SetMarker(marker);
				r0->IntegerDF(num) = 0;
				fcout << "-- face# " << f->LocalID() << " conn# " << cnt << " is boundary connection " << std::endl;
				Storage::reference_array ltrp = f->ReferenceArrayDF(trpl_elems);
				Storage::real_array lcoef = f->RealArrayDF(trpl_coefs);
				Storage::real sum = 0.0;
				stack[0]->RealDF(coe) = 0.0;

			
				for(int q = 0; q < 3; q++)
					if( ltrp[q] != NULL )
					{
						sum += lcoef[q];
						if( !ltrp[q]->GetMarker(marker) )
						{
							stack.push_back(ltrp[q]);
							ltrp[q]->SetMarker(marker);
							ltrp[q]->RealDF(coe) = lcoef[q];
						}
						else ltrp[q]->RealDF(coe) += lcoef[q]; 
					}
				fcout << stack.size() << " ";
				for(int q = 0; q < stack.size(); q++)
				{
					fcout << stack[q]->IntegerDF(unknown_id) << " ";
					stack[q]->RemMarker(marker);
				}
				fcout << 0 << " "; //linear stencil

				fcout << -sum*trans_scale << " ";
				for(int q = 1; q < stack.size(); q++)
					fcout << stack[q]->RealDF(coe) * trans_scale << " ";
				*/

				//fcout << "-- face# " << f->LocalID() << " conn# " << cnt << " is boundary connection " << std::endl;

				Storage::reference_array astncl = f->ReferenceArrayDV(intrp_stncl);
				Storage::real_array acoefs = f->RealArrayDV(intrp_coefs);
				fcout << astncl.size() << " ";
				for(int k = 0; k < astncl.size(); k++)
				{
					fcout << astncl[k]->IntegerDF(unknown_id) << " ";
				}

				fcout << 0 << " ";

				for(int k = 0; k < acoefs.size(); k++)
				{
					fcout << std::scientific << acoefs[k] * trans_scale << " ";
				}
			}
			f->IntegerDF(row) = cnt++;
			fcout << std::endl;
		}
	}
	fcout << "/" << std::endl;
	if( cntsupp )
	{
		fcout << "SUPPCONNS" << std::endl;
		fcout << cntsupp << std::endl;
		cntsupp = 0;
		for (Mesh::iteratorElement f = m->BeginElement(FACE|EDGE|NODE); f != m->EndElement(); ++f)
		{
			if( f->GetMarker(add_markers) )
			{
				if( !f->GetMarker(bnd_markers) || (f->GetElementType() & (EDGE | NODE)) )
				{
					//fcout << "--supp " << ElementTypeName(f->GetElementType()) << "# " << f->LocalID() << " supp# " << cntsupp  << std::endl;
					fcout << 0 << " "; //linear type
				
					/*
					Storage::real cnt[3];
					f->Centroid(cnt);
					fcout << cnt[2] << " "; // depth
					*/
					Storage::reference_array stncl = f->ReferenceArrayDV(intrp_stncl);
					Storage::real_array coefs = f->RealArrayDV(intrp_coefs);
					fcout << stncl.size() << " ";
					for(int k = 0; k < stncl.size(); k++)
					{
						fcout << stncl[k]->IntegerDF(unknown_id) << " ";
					}
					for(int k = 0; k < stncl.size(); k++)
					{
						fcout << std::scientific << coefs[k] << " ";
					}
				}
				else
				{
					//fcout << "--bnd " << ElementTypeName(f->GetElementType()) << "# " << f->LocalID() << " supp# " << cntsupp << std::endl;
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
				cntsupp++;
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


void nmpfa::Update() //TODO
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

nmpfa::~nmpfa()
{
	m->DeleteTag(addval);
	m->DeleteTag(intrp_stncl);
	m->DeleteTag(intrp_coefs);
	m->DeleteTag(trpl_elems);
	m->DeleteTag(trpr_elems);
	m->DeleteTag(trpl_coefs);
	m->DeleteTag(trpr_coefs);
}