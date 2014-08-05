#include "discr.h"
#include <iomanip>
#define EPS 1e-12



typedef Storage::real vec[3];

static Storage::real lengthS(vec x) { return dot_prod(x, x); }
static Storage::real length(vec x) { return sqrt(lengthS(x)); }
static Storage::real normalize(Storage::real x[3]) { Storage::real r = length(x); if (r > 0.0)  x[0] /= r, x[1] /= r, x[2] /= r; return r; }

mpfa_t::mpfa_t(Automatizator * aut, Mesh * m, std::string tensor)
:discr_triplet_basic(aut, m, tensor)
{
	stncl = m->CreateTag("SCHEME_STNCL", DATA_REFERENCE,  FACE , NONE);
	coefs = m->CreateTag("SCHEME_COEFS", DATA_REAL, FACE, NONE);
	interp_stncl = m->CreateTag("INTRP_STNCL", DATA_REFERENCE, FACE | EDGE | NODE, NONE);
	interp_coefs = m->CreateTag("INTRP_COEFS", DATA_REAL, FACE | EDGE | NODE, NONE);
	scheme = aut->RegisterStencil("scheme",stncl,coefs);
	interp = aut->RegisterStencil("interp",interp_stncl,interp_coefs);
}


void mpfa_t::MarkRecursive(Element * f)
{
	Storage::reference_array arr = f->ReferenceArrayDV((f->GetElementType() == FACE && f->GetMarker(bnd_markers)) ? stncl : interp_stncl);
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


void mpfa_t::Init()
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
		Storage::real_array coefs = it->RealArray(interp_coefs);
		Storage::reference_array stncl = it->ReferenceArray(interp_stncl);
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


void mpfa_t::Build(Face * f)
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

		if( !f->GetMarker(bnd_markers) )
		{
			r1->Centroid(cnt2);
			int ret1, ret2, ptypes[4][3];
			ret1 = tensor_prod(Ktype, r0->RealArray(K), nrmf, f1);
			ret2 = tensor_prod(Ktype, r1->RealArray(K), nrmf, f2);

			if (!ret1 || !ret2)	error = 1; //error in tensor

		

			area1 = area*normalize(f1);
			area2 = area*normalize(f2);
			Element * t[2][3];
			Storage::real c[2][3];
		
			ret1 = find_triplet(r0, NULL, f, f1, +1, area1, t[0],c[0],ptypes[0],allowed_types);
			ret2 = find_triplet(r1, NULL, f, f2, -1, area2, t[1],c[1],ptypes[1],allowed_types);


			if (ret1 && ret2)
			{
				Storage::reference_array astncl = f->ReferenceArrayDV(stncl);
				Storage::real_array acoefs = f->RealArrayDV(coefs);
				Storage::real sum[2] = {0,0}, sign[2] = {-0.5,0.5};

				for(int q = 0; q < 3; q++)
				{
					sum[0] += c[0][q];
					sum[1] += c[1][q];
				}

				astncl.push_back(r0);
				acoefs.push_back(sum[0]*sign[0]);
				astncl.push_back(r1);
				acoefs.push_back(sum[1]*sign[1]);

				for(int j = 0; j < 2; j++)
				for(int q = 0; q < 3; q++)
				{
					if(c[j][q]/sum[j] > 1.0e-13 )
					{
						if( t[j][q]->GetElementType() != CELL ) t[j][q]->SetMarker(add_markers);
					
						bool found = false;
						for(int k = 0; k < astncl.size(); k++)
						if( astncl[k] == t[j][q] )
						{
							acoefs[k] -= c[j][q]*sign[j];
							found = true;
							break;
						}
						if( !found )
						{
							astncl.push_back(t[j][q]);
							acoefs.push_back(c[j][q]*(-sign[j]));
						}
					}
				}
			}
			else
			{
				std::cout << "Failed to find triplet!" << std::endl;
				error = 1;
			}
		

		
			ret1 = find_triplet(f, r0, NULL, f1, -1, area1, t[0], c[0],ptypes[0], allowed_types);
			ret2 = find_triplet(f, r1, NULL, f2, +1, area2, t[1], c[1],ptypes[1], allowed_types);

			if( ret1 && ret2 )
			{
				Storage::reference_array astncl = f->ReferenceArrayDV(interp_stncl);
				Storage::real_array acoefs = f->RealArrayDV(interp_coefs);
				Storage::real sum = 0;

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
						for(int k = 0; k < astncl.size(); k++)
						if( astncl[k] == t[j][q] )
						{
							acoefs[k] += c[j][q]/sum;
							found = true;
							break;
						}
						if( !found )
						{
							astncl.push_back(t[j][q]);
							acoefs.push_back(c[j][q]/sum);
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
			int ptypes[3];
			int ret = tensor_prod(Ktype, r0->RealArray(K), nrmf, f1);
			if( ret )
			{
				Element * t[3];
				Storage::real c[3];
				Storage::reference_array astncl = f->ReferenceArrayDV(stncl);
				Storage::real_array acoefs = f->RealArrayDV(coefs);
				area1 = area*normalize(f1);

				astncl.push_back(r0);
				acoefs.push_back(0.0);
				astncl.push_back(f);
				acoefs.push_back(0.0);

				ret = find_triplet(r0, f, NULL, f1, 1, area1,t, c,ptypes,allowed_types);//(allowed_types & (~AVG_NEUMANN)));
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

				ret = find_triplet(f, r0, NULL, f1, -1, area1,t, c,ptypes,allowed_types);//(allowed_types & (~AVG_NEUMANN)));
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




expr mpfa_t::Grad(const expr & param) const
{
	return stencil(scheme,param);
}

expr mpfa_t::Interp(ElementType etype,const expr & param) const
{
	return stencil(interp,param);
}

void mpfa_t::Export(std::ostream & fcout, Storage::real trans_scale, Storage::real vol_scale) const
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
	fcout << "MPFACONNS" << std::endl;
	
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

			if( f->GetMarker(bnd_markers) && !f->GetMarker(add_markers) ) continue;

			Storage::reference_array astncl = f->ReferenceArrayDV(stncl);
			Storage::real_array acoefs = f->RealArrayDV(coefs);
			fcout << astncl.size() << " ";
			for(int k = 0; k < astncl.size(); k++)
			{
				fcout << astncl[k]->IntegerDF(unknown_id) << " ";
				fcout << std::scientific << acoefs[k] * trans_scale << " ";
			}
			
			fcout << std::endl;
			f->IntegerDF(row) = cnt++;
			
		}
		//else if (f->Boundary()) fcout << "-- face# " << f->LocalID() << " is boundary (appear as cell " << f->BackCell()->IntegerDF(unknown_id) << ")" << std::endl;
	}
	fcout << "/" << std::endl;
	if( cntsupp )
	{
		fcout << "SUPPCONNS" << std::endl;
		fcout << cntsupp << std::endl;
		for (Mesh::iteratorElement f = m->BeginElement(FACE|EDGE|NODE); f != m->EndElement(); ++f)
		{
			if( f->GetMarker(add_markers) )
			{
				if( !f->GetMarker(bnd_markers) || (f->GetElementType() & (EDGE | NODE)) )
				{
					//fcout << "--" << ElementTypeName(f->GetElementType()) << "# " << f->LocalID() << std::endl;

					fcout << 0 << " "; //linear type
				
					/*
					Storage::real cnt[3];
					f->Centroid(cnt);
					fcout << cnt[2] << " "; // depth
					*/
					Storage::reference_array astncl = f->ReferenceArrayDV(interp_stncl);
					Storage::real_array acoefs = f->RealArrayDV(interp_coefs);
					fcout << astncl.size() << " ";
					for(int k = 0; k < astncl.size(); k++)
					{
						fcout << astncl[k]->IntegerDF(unknown_id) << " ";
					}
					for(int k = 0; k < astncl.size(); k++)
					{
						fcout << std::scientific << acoefs[k] << " ";
					}
				}
				else
				{
					//fcout << "--bnd " << ElementTypeName(f->GetElementType()) << "# " << f->LocalID() << std::endl;
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


void mpfa_t::Update() //TODO
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

mpfa_t::~mpfa_t()
{
	m->DeleteTag(stncl);
	m->DeleteTag(coefs);
	m->DeleteTag(interp_stncl);
	m->DeleteTag(interp_coefs);
}