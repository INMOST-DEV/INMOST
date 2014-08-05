#include "discr.h"
#define EPS 1e-8

void mpfa_o::Build(Face * f)
{
	if (needBuild(f))
	{
	}
}
mpfa_o::mpfa_o(Automatizator * out, Mesh *m, std::string tensor) :discr_basic(out,m, tensor)
{
	trans = m->CreateTag("MPFA_TRANS_COEF", DATA_REAL, FACE, NONE, 2);
	elems = m->CreateTag("MPFA_TRANS_ELEM", DATA_REFERENCE, FACE, NONE, 2);
	stncl = aut->RegisterStencil("mpfa", elems, trans);
}

void mpfa_o::Init()
{
	if (!have_bnd) Boundary(0, FACE);
	

	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		
		adjacent<Face> faces = it->getFaces();
		Storage::real * mat = new Storage::real[faces.size()*faces.size()];
		int k = 0;
		for(adjacent<Face>::iterator f = faces.begin(); f != faces.end(); f++)
		{
			if( f->GetMarker(bnd_markers) )
			{
				
			}
			else
			{
			}
			k++;
		}
	}

}
expr mpfa_o::Grad(const expr & param) const
{
	return stencil(stncl, param);
}
expr mpfa_o::Interp(ElementType etype,const expr & param) const
{
	return 0;
}
void mpfa_o::Export(std::ostream & fdout, Storage::real trans_scale, Storage::real vol_scale) const
{
	discr_basic::Export(fdout,trans_scale,vol_scale);
	Tag unknown_id = m->CreateTag("ID",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	Tag row = m->CreateTag("ROW",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	
	int cnt = 0, cntsupp = 0;
	INMOST_DATA_ENUM_TYPE idnum = 0;
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		it->IntegerDF(unknown_id) = idnum++;
	if( bnd_conds.isValid() )
		for(Mesh::iteratorElement it = m->BeginElement(FACE); it != m->EndElement(); ++it)
			if(it->HaveData(bnd_conds)) it->IntegerDF(unknown_id) = idnum+cntsupp++;

	cnt = 0;
	fdout << "MPFACONNS" << std::endl;
	for (Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
	{
		Cell * r0 = f->BackCell(), *r1 = f->FrontCell();
		if(needBuild(&*f))
		{
			if( f->GetMarker(bnd_markers) && (!bnd_conds.isValid() || !f->HaveData(bnd_conds))) continue;
			fdout << f->ReferenceArrayDF(elems)[0]->IntegerDF(unknown_id) << " ";
			fdout << f->ReferenceArrayDF(elems)[1]->IntegerDF(unknown_id) << " ";
			fdout << f->RealArrayDF(trans)[0]*trans_scale << std::endl;

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
				fdout << "--bnd " << ElementTypeName(f->GetElementType()) << "# " << f->LocalID() << std::endl;
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
void mpfa_o::Update()
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
mpfa_o::~mpfa_o()
{
	m->DeleteTag(elems);
	m->DeleteTag(trans);
}