#include "inmost.h"

using namespace INMOST;

//todo: add algorithm for cells
// for cell find bounding ellipse, if all three dimensions of ellipse are almost equal replace with node,
// if 2 dimensions dominate, replace with face, project cell onto plane represented by face, detect 2 cells that use this as face, connect others to edges (think about it)
// if 1 dimension dominate, replace with edge
//todo: when collapsing multiple adjacent elements, figure out when to collapse them into one element

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


static bool MatchPoints(const double v1[3], const double v2[3])
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


	std::pair<rMatrix, bool> inv, iQ;
	inv.first.Resize(d + 1, d + 1);
	iQ.first.Resize(d, d);
	int done_iters = 0;
	double n = d+1, ep , en, kp, kn, beta;
	int jp,jn;
	for (int i = 0; i < iters; ++i)
	{

		for (int k = 0; k < p.Rows(); ++k) D(k, k) = p(k, 0);
		inv = (A1*D*A1.Transpose()).Invert(true); // d+1 by d+1
		if (!inv.second) return false;
		M = A1.Transpose()*inv.first*A1; // m by m
		//std::cout << "matrix M:" << std::endl;
		//M.Print();
		kp = -1.0e20;
		kn = 1.0e20;
		jp = -1;
		jn = -1;
		for (int k = 0; k < M.Rows(); ++k)
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
	for (int k = 0; k < p.Rows(); ++k) D(k, k) = p(k, 0);
	iQ = (A*(D - p*p.Transpose())*A.Transpose()).Invert(true);
	if (!iQ.second) return false;
	Q = iQ.first / static_cast<double>(d);
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

	std::cout << "Start:" << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;

	double tt = Timer();
	//There seems to be a problem with VTU reader
	int fixed = 0;
	for (Mesh::iteratorFace f = m.BeginFace(); f != m.EndFace(); ++f)
	{
		if (f->FixNormalOrientation())
			fixed++;
	}
	std::cout << "Time to fix normals: " << Timer() - tt << std::endl;
	std::cout << "Total face normals fixed: " << fixed << std::endl;
	//fix issue with tiny faces/edges and self-intersections
	{
		tt = Timer();
		std::cout << "This program removes tiny faces and edges," << std::endl;
		std::cout << "don't do this if you have legit resolved faults" << std::endl;
		TagInteger cosm = m.CreateTag("fix", DATA_INTEGER, CELL | FACE | EDGE, NONE, 1);

		//mark cells that have to be collapsed
		//todo
		for (Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) cosm[it->self()] = 0;
		//mark faces that have to be collapsed
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
		//report mark results
		std::cout << "faces: " << std::endl;
		for (Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
		if (cosm[it->self()] == 1) std::cout << it->LocalID() << std::endl;
		std::cout << " edges: " << std::endl;
		for (Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it)
		if (cosm[it->self()] == 1) std::cout << it->LocalID() << std::endl;

		TagReferenceArray disconnect = m.CreateTag("DISCONNECT", DATA_REFERENCE, CELL | FACE | EDGE, CELL | FACE | EDGE);
		TagReferenceArray connect = m.CreateTag("CONNECT", DATA_REFERENCE, CELL | FACE | EDGE, CELL | FACE | EDGE);
		TagReference replace = m.CreateTag("REPLACE", DATA_REFERENCE, EDGE | FACE, EDGE | FACE, 1);
		MarkerType collapse = m.CreateMarker();
		//collapsed edges should be replaced with nodes
		//reconnect all adjacent edges later
		for (Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it)
		{
			if (cosm[it->self()] == 1)
			{
				int eid = it->LocalID();
				double v[3];
				for (int k = 0; k < 3; ++k)
					v[k] = (it->getBeg().Coords()[k] + it->getEnd().Coords()[k])*0.5;
				Node n = m.CreateNode(v);
				//replacing edge with a node
				//1.disconnect all adjacent faces from current edge

				//2.in edges sharing nodes with current replace nodes to new node

				//3.disconnect all adjacent nodes from current edge

				//std::cout << "node instead of edge " << it->LocalID() << " ";
				//std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
				replace[it->self()] = n.GetHandle();

				ElementArray<Node> edge_nodes = it->getNodes();
				for (ElementArray<Node>::iterator jt = edge_nodes.begin(); jt != edge_nodes.end(); ++jt)
				{
					ElementArray<Edge> node_edges = jt->getEdges();
					for (ElementArray<Edge>::iterator kt = node_edges.begin(); kt != node_edges.end(); ++kt) if (kt->self() != it->self())
					{
						eid = kt->LocalID();
						disconnect[kt->self()].push_back(jt->self());
						connect[kt->self()].push_back(n);
					}
				}
				ElementArray<Face> edge_faces = it->getFaces();
				for (ElementArray<Face>::iterator jt = edge_faces.begin(); jt != edge_faces.end(); ++jt)
					disconnect[jt->self()].push_back(it->self());

				it->SetMarker(collapse);
			}
		}
		//collapsed faces should be replaced either with edges or nodes
		rMatrix Q(3, 3), c(3, 1), U(3, 3), S(3, 3), V(3, 3);
		MarkerType mrk = m.CreateMarker();
		for (Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
		{
			if (cosm[it->self()] == 1)
			{
				it->getEdges().SetMarker(collapse);
				it->SetMarker(collapse);
				if (BoundingEllipse(it->self(), 5.0e-6, 1000, Q, c))
				{
					//std::cout << "center: "; c.Transpose().Print();
					Q.SVD(U, S, V);
					double d[3], axis[3][3];
					int ek = -1;
					int nk = 0;
					for (int k = 0; k < S.Rows(); ++k)
					{
						if (S(k, k) > 1.0e-8)
						{
							d[nk] = 2 * sqrt(1.0 / S(k, k));
							axis[nk][0] = U(0, k);
							axis[nk][1] = U(1, k);
							axis[nk][2] = U(2, k);
							nk++;
							//std::cout << "axis: " << U(k,0) << " " << U(k,1) << " " << U(k,2);
							//std::cout << "axis: " << U(0,k) << " " << U(1,k) << " " << U(2,k);
							//std::cout << " diameter: " << 2 * sqrt(1.0 / S(k, k)) << " eigen " << S(k,k) << std::endl;
						}
					}
					assert(nk == 2);
					if (d[0] * 3 < d[1]) ek = 1; //create edge along axis 1
					else if (d[1] * 3 < d[0]) ek = 0; //create edge along axis 0
					else ek = -1; //create node
					if (ek == -1)
					{
						//it->Centroid(c.data());
						//replacing face with a node
						Node n = m.CreateNode(c.data());
						//std::cout << "node instead of face " << it->LocalID() << " ";
						//std::cout << c(0,0) << " " << c(1,0) << " " << c(2,0) << std::endl;

						replace[it->self()] = n->GetHandle();
						//1. reconnect edges to node
						//2. disconnect faces from former edges


						ElementArray<Element> face_elems = it->getAdjElements(NODE | EDGE);
						face_elems.SetMarker(mrk);
						ElementArray<Face> adj_faces = it->BridgeAdjacencies2Face(EDGE);
						for (ElementArray<Face>::iterator jt = adj_faces.begin(); jt != adj_faces.end(); ++jt)
						{
							ElementArray<Edge> disconnect_edges = jt->getEdges(mrk);
							disconnect[jt->self()].insert(disconnect[jt->self()].end(), disconnect_edges.begin(), disconnect_edges.end());
						}
						ElementArray<Edge> adj_edges = it->BridgeAdjacencies2Edge(NODE);
						for (ElementArray<Edge>::iterator jt = adj_edges.begin(); jt != adj_edges.end(); ++jt)
						{
							assert(!jt->GetMarker(mrk));
							int eid = jt->LocalID();
							if (jt->getBeg().GetMarker(mrk))
								disconnect[jt->self()].push_back(jt->getBeg());
							else
								disconnect[jt->self()].push_back(jt->getEnd());
							connect[jt->self()].push_back(n);
						}
						ElementArray<Cell> adj_cells = it->getCells();
						for (ElementArray<Cell>::iterator jt = adj_cells.begin(); jt != adj_cells.end(); ++jt)
							disconnect[jt->self()].push_back(it->self());

						face_elems.RemMarker(mrk);

					}
					else
					{
						//it->Centroid(c.data());
						//replacing face with an edge
						//intersect axis with the polygon, represented by current face
						//project both face nodes and axis onto same plane
						//find intersection nodes with the plane and create edge with them
						double vv[2][3]; // expect two coordinates
						int state[2] = { -1, -1 };
						ElementArray<Node> edge_nodes(&m, 2);
						Edge node_edges[2];
						int nvv = 0;
						double orthx[3], orthy[3];
						double nrm[3], cnt[3], pcnt[3];
						it->UnitNormal(nrm);
						it->Centroid(cnt);
						ElementArray<Edge> edges = it->getEdges();
						//compute axises for projection
						edges[0].getBeg().Centroid(pcnt);
						for (int k = 0; k < 3; ++k)
							orthx[k] = pcnt[k] - cnt[k];
						normalize(orthx);
						cross(orthx, nrm, orthy);
						//project l and axis
						double cx, cy, axisx, axisy;
						for (int k = 0; k < 3; ++k)
							pcnt[k] = c.data()[k] - cnt[k];
						cx = dot(pcnt, orthx);
						cy = dot(pcnt, orthy);
						axisx = dot(axis[ek], orthx);
						axisy = dot(axis[ek], orthy);

						//std::cout << "c " << cx << " " << cy << " axis " << axisx << " " << axisy << std::endl;
						//std::cout << "orthx: " << orthx[0] << " " << orthx[1] << " " << orthx[2] << std::endl;
						//std::cout << "orthy: " << orthy[0] << " " << orthy[1] << " " << orthy[2] << std::endl;
						//std::cout << "nrm:   " << nrm[0] << " " << nrm[1] << " " << nrm[2] << std::endl;
						//std::cout << "cnt:   " << cnt[0] << " " << cnt[1] << " " << cnt[2] << std::endl;
						for (int j = 0; j < (int)edges.size(); ++j)
						{
							double vx[3], vy[3], *v, a, b, q;
							v = edges[j].getBeg().Coords().data();
							for (int k = 0; k < 3; ++k)
								pcnt[k] = v[k] - cnt[k];
							vx[0] = dot(orthx, pcnt);
							vy[0] = dot(orthy, pcnt);
							v = edges[j].getEnd().Coords().data();
							for (int k = 0; k < 3; ++k)
								pcnt[k] = v[k] - cnt[k];
							vx[1] = dot(orthx, pcnt);
							vy[1] = dot(orthy, pcnt);

							//std::cout << "j " << j << " v0 " << vx[0] << " " << vy[0] << " v1 " << vx[1] << " " << vy[1];
							//intersect
							//cx + a*axisx
							//cy + a*axisy
							//with
							//vx[0]+b*(vx[1]-vx[0])
							//vy[0]+b*(vy[1]-vy[0])
							// cx + a*axisx = vx[0] + b*(vx[1]-vx[0])
							// cy + a*axisy = vy[0] + b*(vy[1]-vy[0]);
							// a*axisx + b*(vx[0]-vx[1]) = vx[0]-cx
							// a*axisy + b*(vy[0]-vy[1])= vy[0]-cy
							if (axisx)
							{
								q = vy[0] - vy[1] - (vx[0] - vx[1])*axisy / axisx;
								b = (vy[0] - cy - (vx[0] - cx)*axisy / axisx) / q;
							}
							else
							{
								q = vx[0] - vx[1] - (vy[0] - vy[1])*axisx / axisy;
								b = (vx[0] - cx - (vy[0] - cy)*axisx / axisy) / q;
							}
							if (fabs(q) > 1.0e-12) //not parallel
							{
								//line intersects segment
								if (b >= 0 - 1.0e-9 && b <= 1.0 + 1.0e-9)
								{
									vx[2] = vx[0] + b*(vx[1] - vx[0]);
									vy[2] = vy[0] + b*(vy[1] - vy[0]);
									//unproject point into 3d space
									double tvv[3];
									for (int k = 0; k < 3; ++k)
										tvv[k] = orthx[k] * vx[2] + orthy[k] * vy[2] + cnt[k];
									if (nvv >= 2) //expect only two points, otherwise face is non-convex
									{
										std::cout << __FILE__ << ":" << __LINE__ << " the face " << it->LocalID() << " is non-convex" << std::endl;
									}
									else
									{
										int eid = edges[j].LocalID();
										Node n;
										if (edges[j].HaveData(replace))
										{
											n = Node(&m, replace[edges[j]]);
											Storage::real_array c = n.Coords();
											if (!MatchPoints(c.data(), tvv))
											{
												for (int k = 0; k < 3; ++k)
													c[k] = (c[k] + tvv[k])*0.5;
											}
											state[nvv] = 1;
										}
										else
										{
											n = m.CreateNode(tvv);
											replace[edges[j]] = n.GetHandle();
											state[nvv] = 0;
										}
										edge_nodes[nvv] = n;
										node_edges[nvv] = edges[j];
										for (int k = 0; k < 3; ++k)
											vv[nvv][k] = tvv[k];
										nvv++;
									}
								}
							}
						}
						std::cout << "edge instead of face " << it->LocalID() << " ";
						std::cout << vv[0][0] << " " << vv[0][1] << " " << vv[0][2] << " " << state[0] << " ";
						std::cout << vv[1][0] << " " << vv[1][1] << " " << vv[1][2] << " " << state[1] << std::endl;
						Edge e = m.CreateEdge(edge_nodes).first;
						replace[it->self()] = e.GetHandle();

						//for node_edges 

						//for new edge
						ElementArray<Element> face_elems = it->getAdjElements(NODE | EDGE);
						face_elems.SetMarker(mrk);
						ElementArray<Face> adj_faces = it->BridgeAdjacencies2Face(EDGE);
						for (ElementArray<Face>::iterator jt = adj_faces.begin(); jt != adj_faces.end(); ++jt)
						{
							ElementArray<Edge> disconnect_edges = jt->getEdges(mrk);
							for (ElementArray<Edge>::iterator kt = disconnect_edges.begin(); kt != disconnect_edges.end(); ++kt)
							{
								disconnect[jt->self()].push_back(kt->self());
								if (!kt->HaveData(replace))
									connect[jt->self()].push_back(e);
							}
						}
						for (int k = 0; k < 2; ++k) if (state[k] == 0) //if state == 1 this action was already performed
						{
							assert(state[k] != -1);
							ElementArray<Edge> adj_edges = node_edges[k].BridgeAdjacencies2Edge(NODE);
							for (ElementArray<Edge>::iterator jt = adj_edges.begin(); jt != adj_edges.end(); ++jt) if (!jt->GetMarker(mrk))
							{
								int eid = jt->LocalID();
								if (jt->getBeg().GetMarker(mrk))
									disconnect[jt->self()].push_back(jt->getBeg());
								else
									disconnect[jt->self()].push_back(jt->getEnd());
								//std::cout << eid << " " << enid << " " << disid << std::endl;
								connect[jt->self()].push_back(edge_nodes[k]);
							}
						}
						ElementArray<Cell> adj_cells = it->getCells();
						for (ElementArray<Cell>::iterator jt = adj_cells.begin(); jt != adj_cells.end(); ++jt)
							disconnect[jt->self()].push_back(it->self());
						face_elems.RemMarker(mrk);

					}
				}
			}
		}
		m.ReleaseMarker(mrk);
		//collapsed cells should be replaced either with face or edge or node

		//first reconnect faces, then edges
		//for (ElementType etype = CELL; etype > NODE; etype = PrevElementType(etype))
		for (ElementType etype = EDGE; etype <= CELL; etype = NextElementType(etype))
		{
			std::cout << "Reconnect " << ElementTypeName(etype) << std::endl;
			for (Mesh::iteratorElement it = m.BeginElement(etype); it != m.EndElement(); ++it) if (!it->GetMarker(collapse))
			{
				int itid = it->LocalID();
				Storage::reference_array disc, conn;
				disc = disconnect[it->self()];
				conn = connect[it->self()];
				if (!disc.empty() || !conn.empty())
				{
					std::cout << "Element " << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << " ";
					if (!disc.empty())
					{
						std::cout << "disc[" << disc.size() << "]: ";
						for (int k = 0; k < disc.size(); ++k) std::cout << ElementTypeName(disc[k].GetElementType()) << ":" << disc[k].LocalID() << " ";
					}
					if (!conn.empty())
					{
						std::cout << "conn[" << conn.size() << "]: ";
						for (int k = 0; k < conn.size(); ++k) std::cout << ElementTypeName(conn[k].GetElementType()) << ":" << conn[k].LocalID() << " ";
					}
					std::cout << std::endl;
				}
				if (!disc.empty())
					it->Disconnect(disc.data(), disc.size());
				if (!conn.empty())
					it->Connect(conn.data(), conn.size());
			}
		}
		m.ReleaseMarker(collapse, EDGE | FACE);


		//remove unused disconnected elements
		for (Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) if (cosm[it->self()] == 1) it->Delete();
		for (ElementType etype = FACE; etype >= EDGE; etype = PrevElementType(etype))
		{
			for (Mesh::iteratorElement it = m.BeginElement(etype); it != m.EndElement(); ++it)
			{
				if (cosm[it->self()] == 1 && it->nbAdjElements(CELL) != 0)
				{
					ElementArray<Cell> cells = it->getCells();
					std::cout << __FILE__ << ":" << __LINE__ << " Element " << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << " is still connected to " << cells.size() << " cells:";
					for (int k = 0; k < cells.size(); ++k) std::cout << " " << cells[k].LocalID();
					std::cout << std::endl;
				}
				if (it->nbAdjElements(CELL) == 0) it->Delete();
			}
		}
		for (Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it) if (it->nbAdjElements(CELL) == 0) it->Delete();


		std::cout << "Time to fix tiny:" << Timer() - tt << std::endl;
		std::cout << "Cosmetics:" << std::endl;
		std::cout << "Cells: " << m.NumberOfCells() << std::endl;
		std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
		std::cout << "Edges: " << m.NumberOfEdges() << std::endl;


	}
	m.Save(grid_out);
	return 0;
}