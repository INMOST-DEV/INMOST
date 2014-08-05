#include "discr.h"
#include <iomanip>
#define EPS 1e-12

typedef Storage::real vec[3];

static Storage::real lengthS(vec x) { return dot_prod(x, x); }
static Storage::real length(vec x) { return sqrt(lengthS(x)); }
static Storage::real normalize(Storage::real x[3]) { Storage::real r = length(x); if (r > 0.0)  x[0] /= r, x[1] /= r, x[2] /= r; return r; }

nmpfa_a::nmpfa_a(Automatizator * aut, Mesh * m, std::string tensor)
:basic_a(aut, m, tensor)
{
}



expr nmpfa_a::Grad(const expr & param) const
{
	expr T = tagval(add,0);
	expr dp = stencil(trpl, param);
	expr dn = stencil(trpr, param);
	
	expr dp2 = ad_abs(dp) + 1.0e-7;
	expr dn2 = ad_abs(dn) + 1.0e-7;
	expr mp = ad_val(dn2 / (dp2+dn2),0.85);
	expr mn = ad_val(dp2 / (dp2+dn2),0.85);
	
	/*
	expr dp2 = dp*dp / ad_sqrt(dp*dp + 1.0e-5) + 1.0e-7;
	expr dn2 = dn*dn / ad_sqrt(dn*dn + 1.0e-5) + 1.0e-7;
	expr mp = ad_val(dn2 / (dn2 + dp2),0.5);
	expr mn = ad_val(dp2 / (dn2 + dp2),0.5);
	*/
	return T*(stencil(upw,param) - stencil(dwn, param)) + dp*mp+dn*mn;
}

expr nmpfa_a::Interp(ElementType etype,const expr & param) const
{
	return stencil(intrp,param);
}

void nmpfa_a::Export(std::ostream & fcout, Storage::real trans_scale, Storage::real vol_scale) const
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

			
			
			Element * r0 = f->BackCell(), *r1 = f->FrontCell();
			
			int abnP[2] = {0,0};
			
			
			if( f->GetMarker(bnd_markers) && !f->GetMarker(add_markers) ) continue;


			//fcout << "-- face# " << f->LocalID() << " conn# " << cnt;
			//if( f->GetMarker(bnd_markers) ) fcout << " is boundary ";
			//fcout << std::endl;

			{
				
				
				stack.clear();
				stack.push_back(r0);
				r0->SetMarker(marker);
				r0->IntegerDF(num) = 0;
				if( r1 == NULL ) r1 = &*f;

				stack.push_back(r1);
				r1->SetMarker(marker);
				r1->IntegerDF(num) = 1;
				Storage::reference_array lrtrp[2];
				lrtrp[0] = f->ReferenceArrayDF(trpl_elems);
				lrtrp[1] = f->ReferenceArrayDF(trpr_elems);
				Storage::real_array lrcoef[2];
				lrcoef[0] = f->RealArrayDF(trpl_coefs);
				lrcoef[1] = f->RealArrayDF(trpr_coefs);
				Storage::real T = f->RealDF(trans);

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

				fcout << 4 << " "; //nonlinear type

				//fcout << std::endl;
				fcout << abnP[0] << " "; //number of elements that construct left coefficient
				fcout << abnP[1] << " "; //number of elements that construct right coefficient
			
				std::string words[2] = {"-- r0 to r1","-- r1 to r0"};

			
				//fcout << std::endl;
				fcout << std::scientific << T * trans_scale << " "; //for mp
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



nmpfa_a::~nmpfa_a()
{
}