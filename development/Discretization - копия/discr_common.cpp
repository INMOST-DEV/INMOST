#include "discr.h"
#include <iomanip>
void tensor_prod1(Storage::real * K, Storage::real v[3], Storage::real out[3])
{
	out[0] = v[0] * K[0];
	out[1] = v[1] * K[1];
	out[2] = v[2] * K[2];
}
void tensor_prod2(Storage::real * K, Storage::real v[3], Storage::real out[3])
{
	out[0] = v[0] * K[0] + v[1] * K[1] + v[2] * K[2];
	out[1] = v[0] * K[1] + v[1] * K[3] + v[2] * K[4];
	out[2] = v[0] * K[2] + v[1] * K[4] + v[2] * K[5];
}
void tensor_prod3(Storage::real * K, Storage::real v[3], Storage::real out[3])
{
	out[0] = v[0] * K[0] + v[1] * K[1] + v[2] * K[2];
	out[1] = v[0] * K[3] + v[1] * K[4] + v[2] * K[5];
	out[2] = v[0] * K[6] + v[1] * K[7] + v[2] * K[8];
}
int tensor_prod(int Ktype, Storage::real_array & K, Storage::real v[3], Storage::real out[3])
{
	switch (Ktype)
	{
	case 1:
		tensor_prod1(&K[0], v, out);
		break;
	case 2:
		tensor_prod2(&K[0], v, out);
		break;
	case 3:
		tensor_prod3(&K[0], v, out);
		break;
	case 4:
		switch (K.size())
		{
		case 3: tensor_prod1(&K[0], v, out); break;
		case 6: tensor_prod2(&K[0], v, out); break;
		case 9: tensor_prod3(&K[0], v, out); break;
		default:
			std::cout << "Local permiability tensor for ";
			std::cout << " have " << K.size();
			std::cout << " entries and cannot be interpreted" << std::endl;
			return 0;
		}
		break;
	}
	return 1;
}

void tensor_add1(Storage::real * K, Storage::real Kout[9])
{
	Kout[0] += K[0]; 
	Kout[4] += K[1]; 
	Kout[8] += K[2];
}

void tensor_add2(Storage::real * K, Storage::real Kout[9])
{
	Kout[0] += K[0]; 
	Kout[4] += K[3]; 
	Kout[8] += K[5];
	Kout[1] += K[1]; 
	Kout[3] += K[1];
	Kout[2] += K[2];
	Kout[6] += K[2];
	Kout[5] += K[4];
	Kout[7] += K[4];
}

void tensor_add3(Storage::real * K, Storage::real Kout[9])
{
	for (int k = 0; k < 9; k++) Kout[k] += K[k];
}

int tensor_add(int Ktype, Storage::real_array & K, Storage::real Kout[9])
{
	switch (Ktype)
	{
	case 1: tensor_add1(&K[0], Kout); break;
	case 2: tensor_add2(&K[0], Kout); break;
	case 3: tensor_add3(&K[0], Kout); break;
	case 4:
		switch (K.size())
		{
		case 3: tensor_add1(&K[0], Kout); break;
		case 6: tensor_add2(&K[0], Kout); break;
		case 9: tensor_add3(&K[0], Kout); break;
		default:
			std::cout << "Local permiability tensor for ";
			std::cout << " have " << K.size();
			std::cout << " entries and cannot be interpreted" << std::endl;
			return 0;
		}
		break;
	}
	return 1;
}

int average_tensor(int Ktype, Storage::real_array & K1, Storage::real_array & K2, Storage::real Kout[9])
{
	int ret;
	memset(Kout, 0, sizeof(Storage::real) * 9);
	ret = tensor_add(Ktype, K1, Kout); if (!ret) return 0;
	ret = tensor_add(Ktype, K2, Kout); if (!ret) return 0;
	return 1;
}

void vec_diff(Storage::real x[3], Storage::real y[3], Storage::real out[3]) 
{ 
	out[0] = y[0] - x[0];
	out[1] = y[1] - x[1]; 
	out[2] = y[2] - x[2];
}

void cross_prod(Storage::real vecin1[3], Storage::real vecin2[3], Storage::real vecout[3])
{
	Storage::real temp[3];
	temp[0] = vecin1[1]*vecin2[2] - vecin1[2]*vecin2[1];
	temp[1] = vecin1[2]*vecin2[0] - vecin1[0]*vecin2[2];
	temp[2] = vecin1[0]*vecin2[1] - vecin1[1]*vecin2[0];
	vecout[0] = temp[0];
	vecout[1] = temp[1];
	vecout[2] = temp[2];
}
Storage::real dot_prod(Storage::real x[3], Storage::real y[3])  { return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]; }


discr_basic::discr_basic(Automatizator * aut, Mesh * m, std::string tensor_name)
: aut(aut), m(m)
{
	allowed_types = AVG_NONLINEAR | AVG_HARMONIC | AVG_REVDIST;
	bnd_types = NONE;
	rem_marker = 0;
	bnd_markers = 0;
	add_markers = m->CreateMarker();
	error = 0;
	have_bnd = 0;
	if (m->HaveTag(tensor_name))
	{
		K = m->GetTag(tensor_name);
		switch (K.GetSize())
		{
		case 3: Ktype = 1; break;
		case 6: Ktype = 2; break;
		case 9: Ktype = 3; break;
		case ENUMUNDEF: Ktype = 4; break;
		default:
			std::cout << "Permiability tensor have ";
			std::cout << K.GetSize() << " entries ";
			std::cout << "and cannot be interpreted ";
			std::cout << std::endl;
			Ktype = -1;
			error = 1;
		}
	}
	if( m->HaveTag("BOUNDARY_CONDITION") )
	{
		bnd_conds = m->GetTag("BOUNDARY_CONDITION");
		if( bnd_conds.GetSize() != 3 )
		{
			std::cout << "cannot understand BOUNDARY_CONDITIONS data, there should be 3 entries, but I see " << bnd_conds.GetSize() << std::endl;
			bnd_conds = Tag();
		}
	}
}

void discr_basic::Export(std::ostream & fdout, Storage::real trans_scale, Storage::real vol_scale) const
{
	m->AssignGlobalID(CELL);
	fdout << std::setprecision(14);
	fdout << "DIMENS" << std::endl << m->NumberOfCells() << " " << 1 << " " << 1 << "/" << std::endl << "/" << std::endl;
	fdout << "DEPTH" << std::endl;
	if (m->GetDimensions() == 3)
	{
		Storage::integer count = 0;
		Storage::real prev_depth = 1.0e+20, curr_depth, cnt[3];
		for (Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		{
			it->Centroid(cnt);
			curr_depth = cnt[2];
			//fdout << std::scientific << curr_depth << std::endl;
			
			if (fabs(prev_depth - curr_depth) < 1e-5)
				count++;
			else
			{
				if (count > 1)
					fdout << count << "*" << curr_depth << std::endl;
				else if (count != 0)
					fdout << curr_depth << std::endl;
				count = 1;
				prev_depth = curr_depth;
			}	
		}
		if( count > 0 ) fdout << count << "*" << curr_depth << std::endl;
	}
	else
	{
		fdout << m->NumberOfCells() << "*" << "0" << std::endl;
	}
	fdout << "/" << std::endl;
	fdout << "VOLUME" << std::endl;
	Storage::integer count = 0;
	Storage::real prev_vol = 1.0e+20, curr_vol;
	for (Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		curr_vol = it->Volume()*vol_scale;
		//fdout << std::scientific << curr_vol << std::endl;
		if (curr_vol <= 0.0)
		{
			std::cout << __FILE__ << ":" << __LINE__ << " Error: negative volume " << curr_vol << " for cell " << it->GlobalID() << std::endl;
		}
		
		if (fabs(prev_vol - curr_vol) < 1e-5)
			count++;
		else
		{
			if (count > 1)
				fdout << count << "*" << curr_vol << std::endl;
			else if (count != 0)
				fdout << curr_vol << std::endl;
			count = 1;
			prev_vol = curr_vol;
		}
	}
	if (count > 0) fdout << count << "*" << curr_vol << std::endl;
	fdout << "/" << std::endl;

	if( m->HaveTag("PORO") )
	{
		Tag poro = m->GetTag("PORO");
		fdout << "VOLUME" << std::endl;
		Storage::integer count = 0;
		Storage::real prev_poro = 1.0e+20, curr_poro;
		for (Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		{
			curr_poro = it->RealDF(poro);
			if (curr_poro <= 0.0)
			{
				std::cout << __FILE__ << ":" << __LINE__ << " Error: negative porosity " << curr_poro << " for cell " << it->GlobalID() << std::endl;
			}
		
			if (fabs(prev_poro - curr_poro) < 1e-5)
				count++;
			else
			{
				if (count > 1)
					fdout << count << "*" << curr_poro << std::endl;
				else if (count != 0)
					fdout << curr_poro << std::endl;
				count = 1;
				prev_poro = curr_poro;
			}
		}
		if (count > 0) fdout << count << "*" << curr_poro << std::endl;
		fdout << "/" << std::endl;
	}
}


static void UnpackBoundaryInner(Tag tag, Element * e, INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
{
	if (size)
	{
		Storage::integer * recv = static_cast<Storage::integer *>(static_cast<void *>(data));
		Storage::integer * self = &e->IntegerArrayDF(tag)[0];
		Storage::integer f1 = std::min(std::min(recv[0], recv[1]), std::min(self[0], self[1]));
		Storage::integer f2 = std::min(std::max(recv[0], recv[1]), std::max(self[0], self[1]));
		self[0] = f1;
		self[1] = f2;
	}
}

bool discr_basic::needBuild(Face * f) const
{
	Cell * c1 = f->BackCell();
	Cell * c2 = f->FrontCell();
	if (c1 == NULL) return false;
	if (c2 == NULL)
	{
		//return false; // TODO
		if (c1->GetStatus() == Element::Ghost && !c1->GetMarker(BoundaryMarker())) return false;
		return true;
	}
	return !(c1->GetStatus() == Element::Ghost && c2->GetStatus() == Element::Ghost);
}


bool discr_basic::needUpdate(Face * f) const
{
	Cell * c1 = f->BackCell();
	Cell * c2 = f->FrontCell();
	if (c1 == NULL) return false;
	if (c2 == NULL)
	{
		return false; // TODO
		if (c1->GetStatus() == Element::Ghost) return false;
		return c1->New();
	}
	return !(c1->GetStatus() == Element::Ghost && c2->GetStatus() == Element::Ghost) && (c1->New() || c2->New());
}

void discr_basic::Boundary(MIDType bnd_marker, ElementType types)
{
	bnd_types = types;
	if (bnd_marker == 0)
	{
		rem_marker = 1;
		bnd_markers = m->CreateMarker();
		Storage::integer maxval = std::numeric_limits<Storage::integer>::max();
		if (m->GetProcessorsNumber() > 1)
		{
			Tag tag_bnd = m->CreateTag("CALC_BOUNDARY_INNER", DATA_INTEGER, FACE, NONE, 2);
			for (Storage::integer it = 0; it < m->MaxLocalID(FACE); ++it)
			{
				Face * f = m->FaceByLocalID(it);
				if (f != NULL)
				{
					Cell * r0 = f->BackCell(), *r1 = f->FrontCell();
					Storage::integer_array arr = f->IntegerArrayDF(tag_bnd);
					if (r0 != NULL) arr[0] = r0->GlobalID(); else arr[0] = maxval;
					if (r1 != NULL) arr[1] = r1->GlobalID(); else arr[1] = maxval;
					if (arr[0] > arr[1]) //sort
					{
						Storage::integer temp = arr[0];
						arr[0] = arr[1];
						arr[1] = temp;
					}
				}
			}
			m->ReduceData(tag_bnd, FACE, UnpackBoundaryInner);
			m->ExchangeData(tag_bnd, FACE);
			for (Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
			{
				Storage::integer_array arr = f->IntegerArrayDF(tag_bnd);
				if (!(arr[0] != maxval && arr[1] != maxval)) f->SetMarker(bnd_markers);
				else f->RemMarker(bnd_markers);
			}
			m->DeleteTag(tag_bnd);
		}
		else
		{
			for (Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
				if (f->Boundary()) f->SetMarker(bnd_markers); else f->RemMarker(bnd_markers);
		}

		if (types & EDGE)
		{
			for (Mesh::iteratorEdge e = m->BeginEdge(); e != m->EndEdge(); ++e)
			{
				adjacent<Face> around = e->getFaces();
				for (adjacent<Face>::iterator f = around.begin(); f != around.end(); ++f)
					if (f->GetMarker(bnd_markers))
					{
						e->SetMarker(bnd_markers);
						break;
					}
			}
		}
	}
	else
	{
		rem_marker = 0;
		bnd_markers = bnd_marker;
	}
	have_bnd = 1;
}

discr_basic::~discr_basic()
{
	if (rem_marker)
	{
		for (Mesh::iteratorElement it = m->BeginElement(FACE | bnd_types); it != m->EndElement(); ++it)
		{
			it->RemMarker(bnd_markers);
		}
		m->RemMarker(bnd_markers);
	}
	for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE|EDGE|NODE|ESET); it != m->EndElement(); ++it)
		it->RemMarker(add_markers);
	m->ReleaseMarker(add_markers);
}

int discr_basic::isPointInside(Face * f, Storage::real x[3])
{
	Storage::real tri[3][3], btri[3][3], prod[3][3], nc[3];
	f->Centroid(tri[0]);
	adjacent<Node> nodes = f->getNodes();
	Storage::real maxdist = 0, dist;
	for(int i = 0; i < nodes.size(); i++)
	{
		int j = (i+1)%nodes.size();
		nodes[i].Centroid(tri[1]);
		nodes[j].Centroid(tri[2]);
		vec_diff(tri[0],tri[1],nc);
		dist = dot_prod(nc,nc);
		if( dist > maxdist ) maxdist = dist;

		vec_diff(x,tri[0],btri[0]);
		vec_diff(x,tri[1],btri[1]);
		vec_diff(x,tri[2],btri[2]);
		cross_prod(btri[0],btri[1],prod[0]);
		cross_prod(btri[1],btri[2],prod[1]);
		cross_prod(btri[2],btri[0],prod[2]);

		if( dot_prod(prod[0],prod[1]) >= 0 &&
			dot_prod(prod[1],prod[2]) >= 0 &&
			dot_prod(prod[0],prod[2]) >= 0 )
			return 0; //inside one of the triangles
	}
	vec_diff(tri[0],x,nc);
	if( sqrt(dot_prod(nc,nc) / maxdist ) > CRITICAL_DIST_RATIO )
		return 2; //the point is too far
	return 1;
}

int discr_basic::FindHarmonicPoint(Face * fKL, Cell * cK, Cell * cL, Tag tensor, int tensor_type, Storage::real xK[3], Storage::real xL[3], Storage::real y[3], Storage::real & coef)//, Storage::real gamma[3])
{
	Storage::real yK[3], yL[3], xKL[3], nKL[3], dK, dL, lK, lL, lKs[3], lLs[3], D, t;
	Storage::real coefK, coefL, coefQ, coefY, coefDiv;
	Storage::real_array KK = cK->RealArray(tensor), KL = cL->RealArray(tensor);
	fKL->OrientedUnitNormal(cK,nKL);
	fKL->Centroid(xKL);
	tensor_prod(tensor_type,KK,nKL,lKs);
	tensor_prod(tensor_type,KL,nKL,lLs);
	lK = dot_prod(nKL,lKs);
	lL = dot_prod(nKL,lLs);
	lKs[0] -= nKL[0]*lK;
	lKs[1] -= nKL[1]*lK;
	lKs[2] -= nKL[2]*lK;
	lLs[0] -= nKL[0]*lL;
	lLs[1] -= nKL[1]*lL;
	lLs[2] -= nKL[2]*lL;
	D = -dot_prod(xKL,nKL);
	dK = fabs(dot_prod(xK,nKL)+D);
	dL = fabs(dot_prod(xL,nKL)+D);
	yK[0] = xK[0] + dK*nKL[0];
	yK[1] = xK[1] + dK*nKL[1];
	yK[2] = xK[2] + dK*nKL[2];
	yL[0] = xL[0] - dL*nKL[0];
	yL[1] = xL[1] - dL*nKL[1];
	yL[2] = xL[2] - dL*nKL[2];
	coefDiv = (lL*dK+lK*dL);
	coefL = lL*dK/coefDiv;
	coefK = lK*dL/coefDiv;
	coefQ = dK*dL/coefDiv;
	coefY = lK*lL/coefDiv;
	//gamma[0] = (yK[0] - yL[0])*coefY + coefL * lKs[0] + coefK * lLs[0];
	//gamma[1] = (yK[1] - yL[1])*coefY + coefL * lKs[1] + coefK * lLs[1];
	//gamma[2] = (yK[2] - yL[2])*coefY + coefL * lKs[2] + coefK * lLs[2];
	y[0] = yK[0]*coefK + yL[0]*coefL + (lKs[0]-lLs[0])*coefQ;
	y[1] = yK[1]*coefK + yL[1]*coefL + (lKs[1]-lLs[1])*coefQ;
	y[2] = yK[2]*coefK + yL[2]*coefL + (lKs[2]-lLs[2])*coefQ;
	coef = coefL;
	return isPointInside(fKL,y);
}

int discr_basic::FindBoundaryPoint(Face * fK, Cell * cK, Tag tensor, int tensor_type, Storage::real xK[3], Storage::real y[3])
{
	Storage::real nK[3], lKs[3], lK, xfK[3], t;
	Storage::real_array KK = cK->RealArray(tensor);
	fK->Centroid(xfK);
	fK->OrientedUnitNormal(cK,nK);
	tensor_prod(tensor_type,KK,nK,lKs);
	lK = dot_prod(nK,lKs);
	t = (dot_prod(nK,xfK) - dot_prod(nK,xK))/lK;
	y[0] = xK[0] + t*lKs[0];
	y[1] = xK[1] + t*lKs[1];
	y[2] = xK[2] + t*lKs[2];
	return isPointInside(fK,y);
}
/*
void discr_basic::InitializeSupportPoints(int allowed_types)
{
	AvgCoord = m->CreateTag("AVGCOORD",DATA_REAL,FACE,NONE,3);
	AvgCoefs = m->CreateTag("AVGCOEFS",DATA_REAL,FACE,NONE,1);
	AvgStats = m->CreateTag("AVGSTATS",DATA_INTEGER,FACE,NONE,1);
	for(INMOST_DATA_INTEGER_TYPE k = 0; k < m->MaxLocalIDFACE(); ++k)
	{
		Face * f = m->FaceByLocalID(k);
		if(f != NULL)
		{
			Storage::integer & stat = f->Integer(AvgStats);
			Cell * cK = f->BackCell(), * cL = f->FrontCell();
			if(cK != NULL)
			{
				Storage::real_array coord = f->RealArrayDF(AvgCoord);
				Storage::real & coef = f->RealDF(AvgCoefs);
				if(cL != NULL) //interior cell
				{
					if( allowed_types & AVG_HARMONIC )
					{
						Storage::real xK[3], xL[3], coef;
						cK->Centroid(xK);
						cL->Centroid(xL);
						FindHarmonicPoint(f,cK,cL,K,Ktype,xK,xL,&coord[0],coef);
						//Check that the point is inside of the polygon?
						stat = AVG_HARMONIC;
					}
					else if( allowed_types & AVG_CENTROID )
					{
						Storage::real xK[3], xL[3], dK, dL, nKL[3], lKs[3], lLs[3], lK, lL;
						Storage::real_array KK = cK->RealArray(K), KL = cL->RealArray(K);
						f->Centroid(&coord[0]);
						f->OrientedUnitNormal(cK,nKL);
						cK->Centroid(xK);
						cL->Centroid(xL);
						tensor_prod(Ktype,KK,nKL,lKs);
						tensor_prod(Ktype,KL,nKL,lLs);
						lK = dot_prod(nKL,lKs);
						lL = dot_prod(nKL,lLs);
						xK[0] -= coord[0];
						xK[1] -= coord[1];
						xK[2] -= coord[2];
						xL[0] -= coord[0];
						xL[1] -= coord[1];
						xL[2] -= coord[2];
						dK = sqrt(dot_prod(xK,xK));
						dL = sqrt(dot_prod(xL,xL));
						coef = lL*dK/(lL*dK+lK*dL);
						stat = AVG_CENTROID;
					}
					else stat = AVG_NONE;
				}
				else
				{
					if( allowed_types & AVG_CENTROID )
					{
						f->Centroid(&coord[0]);
						coef = 1.0;
						stat = AVG_CENTROID;
					}
					else stat = AVG_NONE;
				}
			} else stat = AVG_NONE;
		}
	}
}
*/