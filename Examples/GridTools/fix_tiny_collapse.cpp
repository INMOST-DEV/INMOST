#include "inmost.h"

using namespace INMOST;

// cell-> node, edge, face
// face-> node, edge
// edge-> node
// based on eigenvalues in bounding ellipse

static void cross(const double a[3], const double b[3], double out[3])
{
	out[0] = a[1] * b[2] - a[2] * b[1];
	out[1] = a[2] * b[0] - a[0] * b[2];
	out[2] = a[0] * b[1] - a[1] * b[0];
}

static double dot(const double a[3], const double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void normalize(double v[3])
{
	double d = sqrt(dot(v, v));
	if (d) for (int k = 0; k < 3; ++k) v[k] /= d;
}


bool MatchPoints(const double v1[3], const double v2[3])
{
	double l = 0;
	for (int k = 0; k < 3; ++k)
		l += (v1[k] - v2[k])*(v1[k] - v2[k]);
	return sqrt(l) < 1.0e-7;
}


//should be similar to
//https://people.orie.cornell.edu/miketodd/TYKhach.pdf
//http://www.bnikolic.co.uk/blog/cpp-khachiyan-min-cov-ellipsoid.html
static bool BoundingEllipse(Element e, double eps, int iters, rMatrix & Q, rMatrix & c)
{
	const bool print = false;
	int d = 3;
	double orthx[3], orthy[3], nrm[3], cnt[3];
	ElementArray<Node> ss = e.getNodes();
	rMatrix A;
	if (e.Planarity())
	{
		double pcnt[3];
		d = 2;
		A.Resize(d, ss.size());
		if (e.GetElementType() == FACE)
			e.getAsFace().UnitNormal(nrm);
		else
		{
			Storage::real_array x0 = ss[0].Coords(), a = x0, b;
			for (ElementArray<Node>::size_type i = 0; i < ss.size(); i++)
			{
				b = ss[(i + 1) % ss.size()].Coords();
				nrm[0] += (a[1] - x0[1])*(b[2] - x0[2]) - (a[2] - x0[2])*(b[1] - x0[1]);
				nrm[1] += (a[2] - x0[2])*(b[0] - x0[0]) - (a[0] - x0[0])*(b[2] - x0[2]);
				nrm[2] += (a[0] - x0[0])*(b[1] - x0[1]) - (a[1] - x0[1])*(b[0] - x0[0]);
				a.swap(b);
			}
			normalize(nrm);
		}
		e.Centroid(cnt);
		//compute axises for projection
		ss[0].Centroid(pcnt);
		for (int k = 0; k < 3; ++k)
			orthx[k] = pcnt[k] - cnt[k];
		normalize(orthx);
		cross(orthx, nrm, orthy);
		for (int j = 0; j < (int)ss.size(); ++j)
		{
			
			double * v = ss[j].Coords().data();
			for (int k = 0; k < 3; ++k)
				pcnt[k] = v[k] - cnt[k];
			A(0, j) = dot(orthx, pcnt);
			A(1, j) = dot(orthy, pcnt);
		}
	}
	else
	{
		A.Resize(d, ss.size());
		for (int j = 0; j < (int)ss.size(); ++j)
		{
			for (size_t k = 0; k < 3; ++k)
				A(k, j) = ss[j].Coords()[k];
		}
	}
	// (x,1)^T (Q,c) (x,1) <= 1
	
	Q.Resize(d, d);
	c.Resize(d, 1);
	double init_p = 1.0 / static_cast<double>(A.Cols());
	//init_p = 1.0;
	rMatrix p(A.Cols(), 1, init_p);// m by 1
	rMatrix _p(A.Cols(), 1); // m by 1
	rMatrix A1 = A.ConcatRows(rMatrix::Row(A.Cols(), 1.0)); // 4 by m
	rMatrix D(A.Cols(), A.Cols()); // m by m
	rMatrix M(A.Cols(), A.Cols());
	
	
	int done_iters = 0;
	double n = d+1, ep , en, kp, kn, beta;
	int jp,jn;
	for (int i = 0; i < iters; ++i)
	{
		
		for (unsigned k = 0; k < p.Rows(); ++k) D(k, k) = p(k, 0);
		M = A1.Transpose()*(A1*D*A1.Transpose()).Solve(A1); // m by m
		//std::cout << "matrix M:" << std::endl;
		//M.Print();
		kp = -1.0e20;
		kn = 1.0e20;
		jp = -1;
		jn = -1;
		for (unsigned k = 0; k < M.Rows(); ++k)
		{
			if (M(k, k) > kp)
			{
				kp = M(k,k);
				jp = k;
			}
			if (M(k,k) < kn)
			{
				kn = M(k,k);
				jn = k;
			}
		}
		assert(jp != -1);
		assert(jn != -1);
		//algorithm 4.2
		ep = (kp / n) - 1;
		en = 1 - (kn / n);
		if (ep > en)
		{
			beta = (kp-n)/(n*(kp-1));
			_p = p*(1 - beta);
			_p(jp, 0) += beta;
		}
		else
		{
			beta = (n-kn)/(n*(kn-1));
			beta = std::min(beta,p(jn,0)/(1-p(jn,0)));
			_p = p*(1+beta);
			_p(jn,0) -= beta;
		}
		
		//algorithm  in nikolic
		//double beta = (kp - d - 1) / ((d + 1)*(kp - 1));
		//_p = p*(1 - beta);
		//_p(jp, 0) += beta;
		done_iters = i;
		if ((_p - p).FrobeniusNorm() < eps) break;
		p.Swap(_p);
		//std::cout << "step " << i << " ceps " << ceps << " max " << m << " norm " << M.FrobeniusNorm() <<  " step size " << l << std::endl;
		
	}
	for (unsigned k = 0; k < p.Rows(); ++k) D(k, k) = p(k, 0);
	Q = (A*(D - p*p.Transpose())*A.Transpose()).PseudoInvert() / static_cast<double>(d);
	c = A*p;
	//check
	if (d == 2)
	{
		if (print)
		{
			std::cout << "2d check " << std::endl;
			std::cout << "orthx: " << orthx[0] << " " << orthx[1] << " " << orthx[2] << std::endl;
			std::cout << "orthy: " << orthy[0] << " " << orthy[1] << " " << orthy[2] << std::endl;
			std::cout << "Q: " << std::endl;
			Q.Print();
			std::cout << "c: "; c.Transpose().Print();
			for (int j = 0; j < (int)ss.size(); ++j)
			{
				rMatrix x = A(0, 2, j, j + 1);
				double t = (x - c).DotProduct(Q*(x - c));
				std::cout << j << ": " << std::setw(14) << t << " " << (t <= 1.0 + eps * 2 ? "good" : "bad ");
				std::cout << " x: "; x.Transpose().Print();
			}
		}
		rMatrix U(2, 2), S(2, 2), V(2, 2), U3(3, 3), S3(3, 3), c3(3, 1);
		Q.SVD(U, S, V);
		//std::cout << "U: " << std::endl;
		//U.Print();
		//std::cout << "S: " << std::endl;
		//S.Print();
		U3.Zero();
		double u0[3], u1[3], u2[3];
		//std::cout << "orthx: " << orthx[0] << " " << orthx[1] << " " << orthx[2] << std::endl;
		//std::cout << "orthy: " << orthy[0] << " " << orthy[1] << " " << orthy[2] << std::endl;
		for (int k = 0; k < 3; ++k)
		{
			u0[k] = U(0, 0)*orthx[k] + U(0, 1)*orthy[k];
			u1[k] = U(1, 0)*orthx[k] + U(1, 1)*orthy[k];
		}
		cross(u0, u1, u2);
		//strange, but it works this way, vectors should go column-wise
		for (int k = 0; k < 3; ++k)
		{
			U3(k, 0) = u0[k];
			U3(k, 1) = u1[k];
			U3(k, 2) = u2[k];
		}
		//std::cout << "U3: " << std::endl;
		//U3.Print();
		S3.Zero();
		S3(0, 0) = S(0, 0);
		S3(1, 1) = S(1, 1);
		Q = U3*S3*U3.Transpose();
		
		
		for (int k = 0; k < 3; ++k)
			c3(k, 0) = c(0, 0)*orthx[k] + c(1, 0)*orthy[k] + cnt[k];
		c = c3;
	}
	if (print)
	{
		std::cout << "iters: " << done_iters << std::endl;
		std::cout << "d: " << d << std::endl;
		std::cout << "Q: " << std::endl;
		Q.Print();
		std::cout << "c: "; c.Transpose().Print();
		for (int j = 0; j < (int)ss.size(); ++j)
		{
			rMatrix x(ss[j].Coords().data(), 3, 1);
			double t = (x - c).DotProduct(Q*(x - c));
			std::cout << j << ": " << std::setw(14) << t << " " << (t <= 1.0 + eps * 2 ? "good" : "bad ");
			std::cout << " x: "; x.Transpose().Print();
		}
	}
	return true;
}

int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " mesh [mesh_out=grid.vtk]" << std::endl;
		return -1;
	}

	std::string grid_out = "grid.vtk";
	
	if (argc > 2) grid_out = std::string(argv[2]);


	Mesh m;
	m.Load(argv[1]);
	
	Mesh::GeomParam table;
	table[MEASURE] = CELL;
	m.PrepareGeometricData(table);
	

	std::cout << "Start:" << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;
	double tt = Timer();
	
	rMatrix Q(3,3), c(3,1), U(3,3), S(3,3), V(3,3);
	TagInteger cosm = m.CreateTag("fix", DATA_INTEGER, CELL | FACE | EDGE | NODE, NONE, 1);
	MarkerType collapse = m.CreateMarker();
	//mark cells that have to be collapsed
	//todo
	for (Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
	{
		ElementArray<Cell> around = it->BridgeAdjacencies2Cell(FACE);
		double vol0 = it->Volume(), vol = vol0;
		for(unsigned l = 0; l < around.size(); ++l)
			vol += around[l].Volume();
		if( vol0 / vol < 0.03 )
		{
			ElementArray<Edge> edges = it->getEdges();
			if (BoundingEllipse(it->self(), 5.0e-6, 1000, Q, c))
			{
				std::cout << "found ellipse" << std::endl;
				Q.Print();
				Q.SVD(U,S,V);
				std::cout << "eigenvalues" << std::endl;
				S.Print();
				double d[3] = {0,0,0};
				int nk = 0;
				for (unsigned k = 0; k < S.Rows(); ++k)
				{
					if (S(k, k) > 1.0e-8)
					{
						d[nk] = 2 * sqrt(1.0 / S(k, k));
						nk++;
					}
				}
				
				double l = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
				std::cout << "count: " << nk << " " << d[0] << " " << d[1] << " " << d[2] << " l " << l << std::endl;
				int nn = 0;
				for(int k = 0; k < 3; ++k)
					nn += (3 * d[k] / l > 1.0) ? 1 : 0;
				std::cout << "case " << nn << std::endl;
				cosm[it->self()] = nn;
				if( nn == 3 ) // collapse into node
				{
					for(unsigned k = 0; k < edges.size(); ++k)
						edges[k].SetMarker(collapse);
				}
				else if( nn == 2 ) // collapse into face, find largest face
				{
					Face f;
					double area = 0;
					ElementArray<Face> faces = it->getFaces();
					for(unsigned k = 0; k < faces.size(); ++k)
					{
						if( faces[k].Area() > area )
						{
							f = faces[k];
							area = faces[k].Area();
						}
					}
					edges.Subtract(f.getEdges());
					for(unsigned k = 0; k < edges.size(); ++k)
						edges[k].SetMarker(collapse);
				}
				else if( nn == 1 ) // collapse into edge, collapse all short edges
				{
					std::sort(edges.begin(),edges.end(),Mesh::MeasureComparator(&m));
					edges[0].SetMarker(collapse); //collapse the smallest
					for(unsigned k = 1; k < edges.size(); ++k)
					{
						if( edges[k].Length() > edges[k-1].Length()*1.5 )
							break;
						else edges[k].SetMarker(collapse);
					}
					for(unsigned k = 0; k < edges.size(); ++k)
						std::cout << "edge " << edges[k].GetHandle() << " length " << edges[k].Length() << " collapse " << (edges[k].GetMarker(collapse) ? 1 : 0) << std::endl;
				}
			}
			std::cout << it->LocalID() << " vol " << vol0 << " total " << vol << std::endl;
		}
	}
	//mark faces that have to be collapsed
	/*
	for (Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
	{
		if (cosm[it->self()] != 1)
		{
			ElementArray<Face> faces = it->getFaces();
			double surface = 0;
			for (ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
				surface += jt->Area();
			for (ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
			{
				if (jt->Area() / surface < 5.0e-3)
				{
					if (cosm[jt->self()] == 0)
					{
						cosm[jt->self()] = 1;
						std::cout << jt->LocalID() << " area " << jt->Area() << " total " << surface << std::endl;
					}
				}
				else cosm[jt->self()] = -1;
			}
		}
	}
	//mark edges that have to be collapsed
	for (Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
	{
		if (cosm[it->self()] != 1) //do not check for collapsed faces
		{
			ElementArray<Edge> edges = it->getEdges();
			double perimeter = 0;
			for (ElementArray<Edge>::iterator jt = edges.begin(); jt != edges.end(); ++jt)
				perimeter += jt->Length();
			for (ElementArray<Edge>::iterator jt = edges.begin(); jt != edges.end(); ++jt)
			{
				if (jt->Length() / perimeter < 5.0e-3)
				{
					if (cosm[jt->self()] == 0) //collapse was not prohibited
					{
						cosm[jt->self()] = 1;
						std::cout << jt->LocalID() << " length " << jt->Length() << " perimeter " << perimeter << std::endl;
					}
				}
				else cosm[jt->self()] = -1; //prohibit collapse
			}
		}
	}
	
	
#pragma omp parallel for reduction(+:nedges)	
	for(int k = 0; k < m.EdgeLastLocalID(); ++k) if( m.isValidEdge(k) )
	{
		Edge e = m.EdgeByLocalID(k);
		ElementArray<Element> around = e.getAdjElements(CELL|FACE|EDGE);
		for(unsigned l = 0; l < around.size(); ++l)
			if( cosm[around[l]] == 1 )
			{
				e.SetMarker(collapse);
				nedges++;
				break;
			}
	}

	std::cout << "collapse edges: " << nedges << std::endl;
	*/
	/*
	for(int k = 0; k < m.EdgeLastLocalID(); ++k) if( m.isValidEdge(k) )
	{
		Edge e = m.EdgeByLocalID(k);
		if( e.GetMarker(collapse) && e.Boundary() )
		{
			std::cout << "Unmark boundary edge " << e.LocalID() << std::endl;
			e.RemMarker(collapse);
		}
	}
	*/
	/*
	std::vector<HandleType> edges;
	for(int k = 0; k < m.EdgeLastLocalID(); ++k) if( m.isValidEdge(k) )
	{
		Edge e = m.EdgeByLocalID(k);
		if( e.GetMarker(collapse) )
			edges.push_back(e.GetHandle());
	}
	std::sort(edges.begin(),edges.end(),Mesh::MeasureComparator(&m));
	 */
	int ncollapsed = 0;
	m.Save("collapsed.pmf");
	//scanf("%*c");
	bool found_edge;
	do
	{
		found_edge = false;
		Edge e = InvalidEdge();
		for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it)
			if( it->GetMarker(collapse) )
			{
				if( e == InvalidEdge() || e.Length() > it->Length() )
				{
					e = it->self();
					found_edge = true;
				}
			}
		if( e.isValid() )
		{
			std::cout << "collapse edge " << e.LocalID() << std::endl;
			Node n = Edge::CollapseEdge(e);
			cosm[n] = 1;
			std::cout << "new node " << n.LocalID() << std::endl;
			ncollapsed++;
		}
			
		m.Save("collapsed.pmf");
		//scanf("%*c");
			
		if( !Element::CheckConnectivity(&m) )
		{
			std::cout << "connectivity is bad" << std::endl;
			throw -1;
		}
	}
	while( found_edge );
	
	
	std::cout << "Collapsed: " << ncollapsed << std::endl;
	std::cout << "Time to fix tiny:" << Timer() - tt << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;

	m.Save(grid_out);
	return 0;
}
