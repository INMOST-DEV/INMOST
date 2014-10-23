#include "discr.h"
#define EPS 1e-8

void tpfa2::Build(Face * f)
{
	if (needBuild(f))
	{
		Storage::real nrmf[3], cnt0[3], cnt1[3], cntf[3], d[3], d1[3], d2[3], f1[3], f2[3], T1, T2, area, T,dd1,dd2, ll1, ll2;
		Element * r0 = f->BackCell();
		Element * r1 = f->FrontCell();
		if (r1 == NULL) r1 = f;
		area = f->Area();
		r0->Barycenter(cnt0);
		r1->Barycenter(cnt1);
		f->Barycenter(cntf);
		f->UnitNormal(nrmf);

		vec_diff(cntf, cnt0, d1);
		vec_diff(cnt1, cntf, d2);
		vec_diff(cnt1, cnt0, d);

		int ret1 = tensor_prod(Ktype, r0->RealArray(K), nrmf, f1);
		int ret2 = r1 == f? memcpy(f2,f1,sizeof(Storage::real)*3),1 : tensor_prod(Ktype, r1->RealArray(K), nrmf, f2);

		if (!ret1 || !ret2)	error = 1; //error in tensor

		dd1 = sqrt(dot_prod(d1, d1));
		dd2 = sqrt(dot_prod(d2, d2));
		ll1 = dot_prod(nrmf, f1);
		ll2 = dot_prod(nrmf, f2);
		if (dd1 > EPS && dd2 > EPS)
			T = ll1*ll2 / (dd1*ll2 + dd2*ll1);
		else if (dd1 > EPS)
			T = ll1 / dd1;
		else if (dd2 > EPS)
			T = ll2 / dd2;
		else
			T = 0.0;

		if( T < 0.0 ) std::cout << "TPFA: negative trans" << std::endl;

		T *= area;
		//check for negative transmissibility?

		Storage::real_array & arr_coefs = f->RealArrayDF(trans);
		Storage::reference_array & arr_elems = f->ReferenceArrayDF(elems);

		arr_elems[0] = r0;
		arr_elems[1] = r1;

		arr_coefs[0] = -T;
		arr_coefs[1] = +T;
	}
}
tpfa2::tpfa2(Automatizator * out, Mesh *m, std::string tensor) :discr_basic(out,m, tensor)
{
	trans = m->CreateTag("TPFA2_TRANS_COEF", DATA_REAL, FACE, NONE, 2);
	elems = m->CreateTag("TPFA2_TRANS_ELEM", DATA_REFERENCE, FACE, NONE, 2);
	stncl = aut->RegisterStencil("tpfa2", elems, trans);
}

void tpfa2::Init()
{
	if (!have_bnd) Boundary(0, FACE);
	for (INMOST_DATA_INTEGER_TYPE id = 0; id < m->MaxLocalID(FACE); ++id)
	{
		Face * f = m->FaceByLocalID(id);
		if (f != NULL) Build(f);
	}
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) it->SetMarker(add_markers);
}
expr tpfa2::Grad(const expr & param) const
{
	return stencil(stncl, param);
}
expr tpfa2::Interp(ElementType etype,const expr & param) const
{
	return 0;
}
void tpfa2::Export(std::ostream & fdout, Storage::real trans_scale, Storage::real vol_scale) const
{
	discr_basic::Export(fdout,trans_scale,vol_scale);
	Tag unknown_id = m->CreateTag("ID",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	Tag row = m->CreateTag("ROW",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	
	int cnt = 0, cntsupp = 0;
	INMOST_DATA_ENUM_TYPE idnum = 0;
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		it->IntegerDF(unknown_id) = idnum++;
	}
	if( bnd_conds.isValid() )
		for(Mesh::iteratorElement it = m->BeginElement(FACE); it != m->EndElement(); ++it)
			if(it->HaveData(bnd_conds)) it->IntegerDF(unknown_id) = idnum+cntsupp++;

	for (Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
	{
		if (needBuild(&*f)) 
		{
			if( f->GetMarker(bnd_markers) && (!bnd_conds.isValid() || !f->HaveData(bnd_conds))) continue;
			cnt++;
		}
	}

	fdout << "TPFACONNS" << std::endl;
	fdout << cnt << std::endl;
	cnt = 0;
	
	for (Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
	{
		Cell * r0 = f->BackCell(), *r1 = f->FrontCell();
		if(needBuild(&*f))
		{
			if( f->GetMarker(bnd_markers) && (!bnd_conds.isValid() || !f->HaveData(bnd_conds))) continue;
			fdout << f->ReferenceArrayDF(elems)[0]->IntegerDF(unknown_id) << " ";
			fdout << f->ReferenceArrayDF(elems)[1]->IntegerDF(unknown_id) << " ";
			fdout << -f->RealArrayDF(trans)[0]*trans_scale << std::endl;

			f->IntegerDF(row) = cnt++;
		}
	}
	fdout << "/" << std::endl;

	if( cntsupp )
	{
		fdout << "SUPPCONNS" << std::endl;
		fdout << cntsupp << std::endl;
		for (Mesh::iteratorElement f = m->BeginElement(FACE); f != m->EndElement(); ++f)
		{
			if( f->GetMarker(bnd_markers) && f->HaveData(bnd_conds) )
			{
				//fdout << "--bnd " << ElementTypeName(f->GetElementType()) << "# " << f->LocalID() << std::endl;
				fdout << 2 << " "; //boundary type
				fdout << f->IntegerDF(row) << " "; //corresponding entry in connection list
				// \alpha C + \beta dC/dn = \gamma
				if( bnd_conds.isValid() && f->HaveData(bnd_conds) )
				{
					Storage::real_array bnd = f->RealArray(bnd_conds);
					fdout << bnd[0] << " "; // \alpha
					fdout << bnd[1] << " "; // \beta
					fdout << bnd[2] << " "; // \gamma
				}
				else
				{
					fdout << 0.0 << " "; // \alpha
					fdout << 1.0 << " "; // \beta
					fdout << 0.0 << " "; // \gamma
				}
				fdout << std::endl;
			}
		}
		fdout << "/" << std::endl;
	}
}
void tpfa2::Update()
{
	for (INMOST_DATA_INTEGER_TYPE id = 0; id < m->MaxLocalID(FACE); ++id)
	{
		Face * f = m->FaceByLocalID(id);
		if (f != NULL)
		{
			if (needBuild(f))
			{
				Cell * r0 = f->BackCell();
				Cell * r1 = f->FrontCell();
				if (r0->New() || r1->New()) Build(f);
			}
		}
	}
}
tpfa2::~tpfa2()
{
	m->DeleteTag(elems);
	m->DeleteTag(trans);
}