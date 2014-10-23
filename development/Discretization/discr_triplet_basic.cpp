#include "discr.h"
#include <iomanip>
#define EPS 1e-12

typedef Storage::real vec[3];

static Storage::real lengthS(vec x) { return dot_prod(x, x); }
static Storage::real length(vec x) { return sqrt(lengthS(x)); }
static Storage::real normalize(Storage::real x[3]) { Storage::real r = length(x); if (r > 0.0)  x[0] /= r, x[1] /= r, x[2] /= r; return r; }

static int compute_triplet_coeffs(vec x, vec x0, vec x1, vec x2, Storage::real c[3])
{
	Storage::real A[3][3], b[3], temp;
	Storage::integer order[3], temp2;
	for(int i = 0; i < 3; i++)
	{
		A[i][0] = x0[i];
		A[i][1] = x1[i];
		A[i][2] = x2[i];
		b[i] = x[i];
		order[i] = i;
	}
	//Factorize matrix
	for(int i = 0; i < 3; i++)
	{
		int maxk = i, maxq = i;
		Storage::real max = fabs(A[maxk][maxq]);
		//Find best pivot
		for(int q = i; q < 3; q++) // over columns
		{
			for(int k = i; k < 3; k++) // over rows
			{
				if( fabs(A[k][q]) > max )
				{
					max = fabs(A[k][q]);
					maxk = k;
					maxq = q;
				}
			}
		}
		//Exchange rows
		if( maxk != i ) 
		{
			for(int q = 0; q < 3; q++)
			{
				temp = A[maxk][q];
				A[maxk][q] = A[i][q];
				A[i][q] = temp;
			}
			//exchange rhs
			{
				temp = b[maxk];
				b[maxk] = b[i];
				b[i] = temp;
			}
		}
		//Exchange columns
		if( maxq != i ) 
		{
			for(int k = 0; k < 3; k++)
			{
				temp = A[k][maxq];
				A[k][maxq] = A[k][i];
				A[k][i] = temp;
			}
			//remember order in sol
			{
				temp2 = order[maxq];
				order[maxq] = order[i];
				order[i] = temp2;
			}
		}
		if( fabs(A[i][i]) < 1.0e-12 )
		{
			if( fabs(b[i]/1.0e-12) > 1.0e-12 ) return -1;
			else A[i][i] = A[i][i] < 0.0 ? - 1.0e-12 : 1.0e-12;
		}
		for(int k = i+1; k < 3; k++)
		{
			A[i][k] /= A[i][i];
			A[k][i] /= A[i][i];
		}
		for(int k = i+1; k < 3; k++)
		for(int q = i+1; q < 3; q++)
		{
			A[k][q] -= A[k][i] * A[i][i] * A[i][q];
		}
		for(int j = i+1; j < 3; j++) //iterate over columns of L
		{
			b[j] -= b[i] * A[j][i];
		}
		b[i] /= A[i][i];
	}

	b[1] -= b[2] * A[1][2];
	b[0] -= b[2] * A[0][2];
	b[0] -= b[1] * A[0][1];
	/*
	for(int i = 2; i >= 0; i--) //iterate over rows of U
		for(int j = i+1; j < 3; j++) 
		{
			b[i] -= b[j] * A[i][j];
		}
	for(int i = 0; i < 3; i++)
		c[order[i]] = b[i];
	*/
	c[order[0]] = b[0];
	c[order[1]] = b[1];
	c[order[2]] = b[2];
	return 0;
}
// Check triplet x0, x1, x2 for vector x
static Storage::real check_triplet(vec x, vec x0, vec x1, vec x2, Storage::real c[3])
{
	Storage::real r = 0.0;
	if (compute_triplet_coeffs(x, x0, x1, x2, c))  return -1.0;
	assert( fabs((c[0]*x0[0]+c[1]*x1[0]+c[2]*x2[0] - x[0])/(fabs(c[0])+fabs(c[1])+fabs(c[2]))) < 1.0e-12 );
	assert( fabs((c[0]*x0[1]+c[1]*x1[1]+c[2]*x2[1] - x[1])/(fabs(c[0])+fabs(c[1])+fabs(c[2]))) < 1.0e-12 );
	assert( fabs((c[0]*x0[2]+c[1]*x1[2]+c[2]*x2[2] - x[2])/(fabs(c[0])+fabs(c[1])+fabs(c[2]))) < 1.0e-12 );
	if (c[0] > r)  r = c[0];
	if (c[1] > r)  r = c[1];
	if (c[2] > r)  r = c[2];
	return ((c[0] >= -1.0e-15) && (c[1] >= -1.0e-15) && (c[2] >= -1.0e-15)) ? r : -1.0; // Check it!!!
	//return ((c[0] >= 0.0) && (c[1] >= 0.0) && (c[2] >= 0.0)) ? r : -1.0;
}

static int find_triplet_sphere2(int n, vec *x, vec v, int ngood, int nnln, int nbad, int t[3], Storage::real c[3], bool necessery)
{
	std::pair<Storage::real,int> dist_pos[128];
	Storage::real r,cc[3], best = -1.0;
	vec tmp;
	int i, j, k, b;
	c[0] = c[1] = c[2] = -1;
	//int m;
	// compute distances
	for (i = 0; i<n; i++)
	{
		dist_pos[i].second = i;
		vec_diff(v, x[i], tmp);
		dist_pos[i].first = lengthS(tmp);
	}
	// sort according to distance
	std::sort(dist_pos+(necessery?1:0),dist_pos+ngood);
	std::sort(dist_pos+ngood,dist_pos+ngood+nnln);
	std::sort(dist_pos+ngood+nnln,dist_pos+ngood+nnln+nbad);

	//std::sort(dist_pos+(necessery?1:0),dist_pos+ngood+nnln+nbad);

	if( dist_pos[0].first < 1.0e-9 )
	{
		t[0] = dist_pos[0].second;
		t[1] = dist_pos[1].second;
		t[2] = dist_pos[2].second;
		c[0] = 1.0;
		c[1] = 0.0;
		c[2] = 0.0;
		return 1;
	}

	for (b = 1, k = 2; b && (k<n); k++)
	{ // order of for-loops is corrected 
		for (j = 1; b && (j<k); j++)
		{
			for (i = 0; (b && (i<j)) && (!necessery || i < 1); i++)
			{ // check points i,j,k 
				r = check_triplet(v, x[dist_pos[i].second], x[dist_pos[j].second], x[dist_pos[k].second],cc);
				if (r < 0.0 )  continue;
				//if( necessery && cc[0]/sqrt(cc[0]*cc[0]+cc[1]*cc[1]+cc[2]*cc[2]) < 1.0e-5 ) continue;
				//if (r > 1.0 ) continue;
				
				if( necessery )
				{
					if( best >= 0.0 && fabs(cc[0]-1) > fabs(best-1))  continue; // we seek those decompositions that are closer to 1
					//if( best >= 0.0 && fabs(cc[0]) > fabs(best))  continue;
					if( cc[0] < 1.0e-6 ) continue;
					
					best = cc[0];
				}
				else
				{
					if( best >= 0.0 )
					{
						//if( fabs(r-1) > fabs(best-1) ) continue;
						if( fabs(r-1) > fabs(best-1) ) continue;
						//if(	necessery && fabs(cc[0]) < fabs(c[0]) ) continue;
					}
					//if ( (best >= 0.0 && best >= r) && !(best <= 1.0 && r > 1.0) && !(best >= 1.0 && r < best) )  continue; // we seek those decompositions that are closer to 1
					best = r;
				}
				t[0] = dist_pos[i].second, t[1] = dist_pos[j].second, t[2] = dist_pos[k].second;
				c[0] = cc[0], c[1] = cc[1], c[2] = cc[2];
			}
		}
	}
	//if( necessery && best > 0 ) std::cout << c[0] << " " << c[1] << " " << c[2] << " " << best << std::endl;
	return (best < 0.0) ? 0 : 1;
}

int discr_triplet_basic::find_triplet(Element * c, Element * include, Face * harmonic_face,
									  Storage::real v[3], Storage::real dir, 
									  Storage::real area, Element * t[3], 
									  Storage::real coefs[3], int ptypes[3],
									  int allowed, bool f_necessery)
{
	Storage::real cc[3];
	int ngood, nnln, nbad;
	int res, k, tt[3];
	bool necessery = (include != NULL && f_necessery) || (harmonic_face != NULL);
	vec v0, cnt, x[128];
	Element * elem[128];
	Storage::real dist[128];
	Storage::real coef[128];
	Storage::real vd;
	//Storage::real gamma[3];
	v0[0] = dir * v[0], v0[1] = dir * v[1], v0[2] = dir * v[2];
	vd = normalize(v0); //for any case!
	t[0] = NULL, t[1] = NULL, t[2] = NULL;
	k = 0;
	c->Centroid(cnt);
	{
		if(c->GetElementType() == CELL && harmonic_face != NULL )
		{
			Cell * xL = harmonic_face->FrontCell();
			Cell * xK = harmonic_face->BackCell();
			elem[k] = xK == c ? xL : xK;
			elem[k]->Centroid(x[k]);
			FindHarmonicPoint(harmonic_face,c->getAsCell(),elem[k]->getAsCell(),K,Ktype,cnt,x[k],x[k],coef[k]);//,gamma);
			vec_diff(cnt,x[k],x[k]);
			dist[k] = normalize(x[k]);
			if( dist[k] > 1.0e-15 ) k++;
		}
		if(include != NULL)
		{
			include->Centroid(x[k]);
			vec_diff(cnt,x[k],x[k]);
			dist[k] = normalize(x[k]);
			coef[k] = 1.0;
			elem[k] = include;
			if( dist[k] > 1.0e-15 ) k++;
		}
		
		adjacent<Cell> cells = c->getCells();
		//Add central point
		for(adjacent<Cell>::iterator it = cells.begin(); it != cells.end(); ++it)
		{
			if( &*it != include )
			{
				it->Centroid(x[k]);
				vec_diff(cnt,x[k],x[k]);
				dist[k] = normalize(x[k]);
				coef[k] = 1.0;
				elem[k] = &*it;
				if( dist[k] > 1.0e-15 )// && dot_prod(v,x[k]) > 0.0)
					k++;
			}
		}
		for(adjacent<Cell>::iterator it = cells.begin(); it != cells.end(); ++it)
		{
			adjacent<Face> faces = it->getFaces();

			for(adjacent<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
			{
				int ret = 2;
				if(jt->GetMarker(bnd_markers) )
				{

					/*
					coef[k] = 1.0;
					elem[k] = jt->BackCell();
					jt->Centroid(x[k]);
					vec_diff(cnt,x[k],x[k]);
					dist[k] = normalize(x[k]);
					if( dist[k] > 1.0e-15 ) k++;
					*/
					bool not_neumann = bnd_conds.isValid();

					if( not_neumann )
					{
						if( !jt->HaveData(bnd_conds) )
							not_neumann = false;
						else
						{
							Storage::real_array bnd = jt->RealArray(bnd_conds);
							if( fabs(bnd[0]-0.0) < 1.0e-8 && fabs(bnd[1]-1.0) < 1.0e-8 && fabs(bnd[2]-0.0) < 1.0e-8 )
								not_neumann = false;
						}
					}

					if( (allowed & AVG_NEUMANN) && !not_neumann)
					{
						Cell * xK = jt->BackCell();
						elem[k] = xK;
						elem[k]->Centroid(x[k]);
						coef[k] = 1.0;
						ret = FindBoundaryPoint(&*jt,xK,K,Ktype,x[k],x[k]);
					}
					
				}
				else if( (allowed & AVG_HARMONIC) && (&*jt != harmonic_face))
				{
					Cell * xL = jt->FrontCell();
					Cell * xK = jt->BackCell();
					elem[k] = xK == &*it ? xL : xK;
					elem[k]->Centroid(x[k]);
					ret = FindHarmonicPoint(&*jt,&*it,elem[k]->getAsCell(),K,Ktype,cnt,x[k],x[k],coef[k]);//,gamma);
				}
				if( ret < 2 )
				{
					vec_diff(cnt,x[k],x[k]);
					dist[k] = normalize(x[k]);
					if( dist[k] > 1.0e-15 ) k++;
				}
			}
		}
	}
	
	ngood = k;
	res = k >= 3 ? find_triplet_sphere2(k,x,v0,ngood,0,0,tt,cc, necessery) : 0;
	if( res == 0 ) // failed to find triplet
	{
		if( allowed & AVG_NONLINEAR )
		{
			adjacent<Element> faces = c->BridgeAdjacencies(CELL,FACE);
			//Add face centroids to array
			//std::cout << "Warning: have to add points with nonlinear flux" << std::endl;
			//scanf("%*c");
			for(adjacent<Element>::iterator it = faces.begin(); it != faces.end(); ++it)
			{
				it->Centroid(x[k]);
				vec_diff(cnt,x[k],x[k]);
				dist[k] = normalize(x[k]);
				if( dist[k] > 1.0e-15 ) 
				{
					//elem[k] = it->GetMarker(bnd_markers) ? it->getAsFace()->BackCell() : &*it; //TODO Should not replace with cell
					elem[k] = &*it; //TODO Should not replace with cell
					coef[k] = 1.0;
					k++;
				}
			}
		}
		nnln = k-ngood;
		res = (k >= 3 && nnln > 0 ) ? find_triplet_sphere2(k,x,v0,ngood,nnln,0,tt,cc, necessery) : 0;
		if( res == 0 )
		{
			if( allowed & AVG_REVDIST )
			{
				//std::cout << "Warning: adding values on edges to approximate flux" << std::endl;
				//scanf("%*c");
				adjacent<Element> edges = c->BridgeAdjacencies(CELL | c->GetElementType(),EDGE);
				//std::cout << edges.size() << std::endl;
				for(adjacent<Edge>::iterator it = edges.begin(); it != edges.end(); ++it)
				{
					it->Centroid(x[k]);
					vec_diff(cnt,x[k],x[k]);
					dist[k] = normalize(x[k]);
					if( dist[k] > 1.0e-15 ) 
					{
						elem[k] = &*it;
						coef[k] = 1.0;
						k++;
					}
				}
				nbad = k - ngood - nnln;
				res = k >= 3 ? find_triplet_sphere2(k,x,v0,ngood,nnln,nbad,tt,cc, necessery) : 0;
				if( res == 0 )
				{
					//std::cout << "Warning: adding values on nodes to approximate flux" << std::endl;
					//scanf("%*c");
					adjacent<Element> nodes = c->BridgeAdjacencies(CELL | c->GetElementType(),NODE);
					//std::cout << nodes.size() << std::endl;
					for(adjacent<Node>::iterator it = nodes.begin(); it != nodes.end(); ++it)
					{
						it->Centroid(x[k]);
						vec_diff(cnt,x[k],x[k]);
						dist[k] = normalize(x[k]);
						if( dist[k] > 1.0e-15 && dot_prod(v,x[k]) > 0.0 ) 
						{
							elem[k] = &*it;
							coef[k] = 1.0;
							k++;
						}
					}
					nbad = k - ngood - nnln;
					res = k >= 3 ? find_triplet_sphere2(k,x,v0,ngood,nnln,nbad,tt,cc, necessery) : 0;
				}
			}
		}
		
	}

	if( res )
	{
		//compute_triplet_coeffs(v0,x[tt[0]],x[tt[1]],x[tt[2]],&coefs[0]);
		for(int i = 0; i < 3; ++i)
		{
			coefs[i] = vd * area / dist[tt[i]] * coef[tt[i]] * cc[i];
			t[i] = elem[tt[i]];
			if( tt[i] < ngood ) ptypes[i] = AVG_NONE;
			else if( tt[i] < ngood+nnln )
			{
				if( t[i]->GetElementType() == FACE )
					ptypes[i] = AVG_NONLINEAR;
				else
					ptypes[i] = AVG_NONE;
			}
			else ptypes[i] = AVG_REVDIST;

			//if( t[i]->GetMarker(bnd_markers) || (ptypes[i] == AVG_NONLINEAR && t[i]->GetElementType() == CELL))
			//{
			//	std::cout << "Used boundary" << std::endl;
			//}
		}
		/*
		std::cout << ElementTypeName(c->GetElementType()) << " at " << cnt[0] << " " << cnt[1] << " " << cnt[2] << std::endl;
		std::cout << "conormal: " << v0[0] << " " << v0[1] << " " << v0[2] << std::endl;
		std::cout << "chosen directions: ";
		for(int i = 0; i < 3; i++)
			std::cout << tt[i] << " ";
		std::cout << "coefficients: ";
		for(int i = 0; i < 3; i++)
			std::cout << coefs[i] << " ";
		std::cout << std::endl;
		std::cout << "all directions: " << std::endl;
		for(int i = 0; i < k; i++)
			std::cout << "vec " << i << " " << x[i][0] << " " << x[i][1] << " " << x[i][2] << " " << ElementTypeName(elem[i]->GetElementType()) << " " << coef[i] << std::endl;
		*/
	}
	else 
	{
		std::cout << "Error: cannot find triplet, total points gathered " << k << std::endl;
		std::cout << "conormal: " << v[0] << " " << v[1] << " " << v[2] << std::endl;
		for(int i = 0; i < k; i++)
			std::cout << "vec " << i << " " << x[i][0] << " " << x[i][1] << " " << x[i][2] << std::endl;
	}
	//return length to previous state
	//v[0] *= vd;
	//v[1] *= vd;
	//v[2] *= vd;

	return res;
}



discr_triplet_basic::~discr_triplet_basic()
{
}