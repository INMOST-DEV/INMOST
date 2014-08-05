#include "discr.h"
#include <iomanip>
#define EPS 1e-12

typedef Storage::real vec[3];

static Storage::real lengthS(vec x) { return dot_prod(x, x); }
static Storage::real length(vec x) { return sqrt(lengthS(x)); }
static Storage::real normalize(Storage::real x[3]) { Storage::real r = length(x); if (r > 0.0)  x[0] /= r, x[1] /= r, x[2] /= r; return r; }

basic_a::basic_a(Automatizator * aut, Mesh * m, std::string tensor)
:discr_triplet_basic(aut, m, tensor)
{
	trpl_elems = m->CreateTag("A_TRP_L", DATA_REFERENCE, FACE, NONE, 4);
	trpr_elems = m->CreateTag("A_TRP_R", DATA_REFERENCE, FACE, NONE, 4);
	trpl_coefs = m->CreateTag("A_COEF_L", DATA_REAL, FACE, NONE, 4);
	trpr_coefs = m->CreateTag("A_COEF_R", DATA_REAL, FACE, NONE, 4);
	trans = m->CreateTag("A_TRANS", DATA_REAL, FACE, NONE, 1);

	
	trpr = aut->RegisterStencil("a_right_triplet", trpr_elems, trpr_coefs);
	trpl = aut->RegisterStencil("a_left_triplet", trpl_elems, trpl_coefs);
	add = aut->RegisterStaticTag(trans);
	intrp_stncl = m->CreateTag("INTRP_STNCL", DATA_REFERENCE, FACE | EDGE | NODE, NONE);
	intrp_coefs = m->CreateTag("INTRP_COEFS", DATA_REAL, FACE | EDGE | NODE, NONE);
	
	
	dwn = aut->RegisterStencil("back_cell", get_l);
	upw = aut->RegisterStencil("front_cell", get_r);

	
	intrp = aut->RegisterStencil("interp",intrp_stncl,intrp_coefs,add_markers);
}


void basic_a::MarkRecursive(Element * f)
{
	if( f->GetElementType() == FACE && f->GetMarker(bnd_markers) )
	{
		Storage::reference_array lrtrp[2];
		lrtrp[0] = f->ReferenceArrayDF(trpl_elems);
		lrtrp[1] = f->ReferenceArrayDF(trpr_elems);
		for(int i = 0; i < 2; i++)
			for(Storage::reference_array::iterator it = lrtrp[i].begin(); it != lrtrp[i].end(); ++it) if( (*it) != NULL )
			{
				if( (*it)->GetElementType() != CELL && !(*it)->GetMarker(add_markers))
				{
					(*it)->SetMarker(add_markers);
					(*it)->SetMarker(mark_next);
				}
			}
		
		for(int i = 0; i < 2; i++)
			for(Storage::reference_array::iterator it = lrtrp[i].begin(); it != lrtrp[i].end(); ++it) if( (*it) != NULL )
			{
				if( (*it)->GetMarker(mark_next) )
				{
					(*it)->RemMarker(mark_next);
					MarkRecursive(*it);
				}
			}
	}
	else
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
}


void basic_a::Init()
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


void basic_a::get_l(Storage * e, Automatizator::stencil_pairs & out, void * user_data)
{
	Face * f = static_cast<Face *>(e);
	Element * r = f->BackCell();
	out.push_back(Automatizator::stencil_pair(r == NULL ? f : r, 1.0));
}
void basic_a::get_r(Storage * e, Automatizator::stencil_pairs & out, void * user_data)
{
	Face * f = static_cast<Face *>(e);
	Element * r = f->FrontCell();
	out.push_back(Automatizator::stencil_pair(r == NULL ? f : r, 1.0));
}


void basic_a::Build(Face * f)
{
	if (needBuild(f))
	{
		int ptypes[2][3];
		Storage::real & T = f->RealDF(trans);
		Storage::reference_array ltrp = f->ReferenceArrayDF(trpl_elems);
		Storage::reference_array rtrp = f->ReferenceArrayDF(trpr_elems);
		Storage::real_array lcoef = f->RealArrayDF(trpl_coefs);
		Storage::real_array rcoef = f->RealArrayDF(trpr_coefs);

		if( !f->GetMarker(bnd_markers) )
		{
			int ret1,ret2;
			
			Storage::real xK[3],xL[3], gamma[3];
			
			
			Cell * cK = f->BackCell(), * cL = f->FrontCell();



			Storage::real yK[3], yL[3], xKL[3], nKL[3], dK, dL, lKn, lLn, lKv[3], lLv[3],lKt[3], lLt[3], D;
			Storage::real coefK, coefL, coefY, coefDiv;
			Storage::real area = f->Area();
			Storage::real_array KK = cK->RealArrayDF(K), KL = cL->RealArrayDF(K);
			cK->Centroid(xK);
			cL->Centroid(xL);
			f->UnitNormal(nKL);
			f->Centroid(xKL);
			tensor_prod(Ktype,KK,nKL,lKv);
			tensor_prod(Ktype,KL,nKL,lLv);
			lKn = dot_prod(nKL,lKv);
			lLn = dot_prod(nKL,lLv);
			lKt[0] = lKv[0] - nKL[0]*lKn;
			lKt[1] = lKv[1] - nKL[1]*lKn;
			lKt[2] = lKv[2] - nKL[2]*lKn;
			lLt[0] = lLv[0] - nKL[0]*lLn;
			lLt[1] = lLv[1] - nKL[1]*lLn;
			lLt[2] = lLv[2] - nKL[2]*lLn;
			D = -dot_prod(xKL,nKL);
			dK = fabs(dot_prod(xK,nKL)+D);
			dL = fabs(dot_prod(xL,nKL)+D);
			yK[0] = xK[0] + dK*nKL[0];
			yK[1] = xK[1] + dK*nKL[1];
			yK[2] = xK[2] + dK*nKL[2];
			yL[0] = xL[0] - dL*nKL[0];
			yL[1] = xL[1] - dL*nKL[1];
			yL[2] = xL[2] - dL*nKL[2];
			coefDiv = (lLn*dK+lKn*dL);
			coefL = lLn*dK/coefDiv;
			coefK = lKn*dL/coefDiv;
			coefY = lKn*lLn/coefDiv;
			//uncomputed direction:
			gamma[0] = (yK[0] - yL[0])*coefY + coefL * lKt[0] + coefK * lLt[0];
			gamma[1] = (yK[1] - yL[1])*coefY + coefL * lKt[1] + coefK * lLt[1];
			gamma[2] = (yK[2] - yL[2])*coefY + coefL * lKt[2] + coefK * lLt[2];
			

			T = coefY * area; //two-point flux transmissibility

			Storage::real g = normalize(gamma);


			if( g > 1.0e-9 )
			{
				ret1 = find_triplet(cK,NULL, NULL, gamma, +1, area*g, &ltrp[0], &lcoef[0],ptypes[0],allowed_types);
				ret2 = find_triplet(cL,NULL, NULL, gamma, -1, area*g, &rtrp[0], &rcoef[0],ptypes[1],allowed_types);

				
				if( ret1 && ret2 )
				{
					Storage::real sum[2] = {0,0}, deny[2] = {0,0};
					for(int q = 0; q < 3; q++)
					{
						sum[0] += lcoef[q];
						sum[1] += rcoef[q];
					}

					for(int q = 0; q < 3; q++)
					{
						if( lcoef[q]/sum[0] < 1.0e-6 || lcoef[q] < 1.0e-11 )
						{
							deny[0] += lcoef[q];
							lcoef[q] = 0.0;
							ltrp[q] = NULL;
						}
						else 
						{
							if( ltrp[q]->GetElementType() != CELL ) ltrp[q]->SetMarker(add_markers);
						}
						if( rcoef[q]/sum[1] < 1.0e-6 || rcoef[q] < 1.0e-11 )
						{
							deny[1] += rcoef[q];
							rcoef[q] = 0.0;
							rtrp[q] = NULL;
						}
						else 
						{
							if( rtrp[q]->GetElementType() != CELL ) rtrp[q]->SetMarker(add_markers);
						}
						rcoef[q] = -rcoef[q];
					}
					ltrp[3] = cK;
					lcoef[3] = -(sum[0]-deny[0]);
					rtrp[3] = cL;
					rcoef[3] = (sum[1]-deny[1]);
				}
				else std::cout << "cannot find nonlinear correction to two-point flux, face " << f->LocalID() << std::endl;
			}
			else
			{
				for(int q = 0; q < 4; q++)
				{
					ltrp[q] = rtrp[q] = NULL;
					lcoef[q] = rcoef[q] = 0.0;
				}
			}
		

			Element * t[2][3];
			Storage::real c[2][3];
			ret1 = find_triplet(f, cK, NULL, lKv, -1, 1.0, t[0], c[0],ptypes[0],allowed_types);
			ret2 = find_triplet(f, cL, NULL, lLv, +1, 1.0, t[1], c[1],ptypes[1],allowed_types);

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
			
			int ret1,ret2;
			Storage::real nK[3],lKn, lKv[3], lKt[3], area = f->Area(), cntf[3], cntK[3], dK,D,g;
			Cell * cK = f->BackCell();
			f->UnitNormal(nK);
			f->Centroid(cntf);
			cK->Centroid(cntK);
			tensor_prod(Ktype,cK->RealArray(K),nK,lKv);
			lKn = dot_prod(nK,lKv);
			lKt[0] = lKv[0] - nK[0] * lKn;
			lKt[1] = lKv[1] - nK[1] * lKn;
			lKt[2] = lKv[2] - nK[2] * lKn;
			D = -dot_prod(cntf,nK);
			dK = fabs(dot_prod(cntK,nK)+D);
			T = lKn / dK * area;
			
			g = normalize(lKt);

			if( g > 1.0e-7 )
			{
				ret1 = find_triplet(cK,f,NULL,lKt,+1,area*g,&ltrp[0],&lcoef[0],ptypes[0],allowed_types);
				ret2 = find_triplet(f,cK,NULL,lKt,-1,area*g,&rtrp[0],&rcoef[0],ptypes[1],allowed_types);


				if( ret1 && ret2 )
				{
					Storage::real sum[2] = {0,0}, deny[2] = {0,0};
					for(int q = 0; q < 3; q++)
					{
						sum[0] += lcoef[q];
						sum[1] += rcoef[q];
					}

					for(int q = 0; q < 3; q++)
					{
						if( lcoef[q]/sum[0] < 1.0e-6 || lcoef[q] < 1.0e-11 )
						{
							deny[0] += lcoef[q];
							lcoef[q] = 0.0;
							ltrp[q] = NULL;
						}
						if( rcoef[q]/sum[1] < 1.0e-6 || rcoef[q] < 1.0e-11 )
						{
							deny[1] += rcoef[q];
							rcoef[q] = 0.0;
							rtrp[q] = NULL;
						}
						rcoef[q] = -rcoef[q];
					}
					ltrp[3] = cK;
					lcoef[3] = -(sum[0]-deny[0]);
					rtrp[3] = f;
					rcoef[3] = (sum[1]-deny[1]);
				}
			}
			else
			{
				for(int q = 0; q < 4; q++)
				{
					ltrp[q] = rtrp[q] = NULL;
					lcoef[q] = rcoef[q] = 0.0;
				}
			}
			
		}
	}
}



void basic_a::Update() //TODO
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

basic_a::~basic_a()
{
	m->DeleteTag(trans);
	m->DeleteTag(intrp_stncl);
	m->DeleteTag(intrp_coefs);
	m->DeleteTag(trpl_elems);
	m->DeleteTag(trpr_elems);
	m->DeleteTag(trpl_coefs);
	m->DeleteTag(trpr_coefs);
}