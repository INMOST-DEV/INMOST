#include "streamline.h"
#include "inc_glut.h"
#include "color.h"
#include "svg_line.h"

namespace INMOST
{
	
	void GetVelocity(Element c, const Tag & velocity_tag, ElementType vel_adj, coord pnt, coord & ret)
	{
		coord cnt;
		const Storage::real eps = 1.0e-8;
		Storage::real dist = 0;
		ret[0] = ret[1] = ret[2] = 0;
		c->Centroid(cnt.data());
		if ((cnt - pnt).length() < 1.0e-5)
		{
			ret = coord(c->RealArray(velocity_tag).data());
		}
		else //inverse distance algorithm (consider wlsqr with linear basis)
		{
			ElementArray<Element> adj = c->BridgeAdjacencies(vel_adj == NODE ? CELL : NODE,vel_adj);
			adj.push_back(c);
			for (ElementArray<Element>::iterator it = adj.begin(); it != adj.end(); ++it)
			{
				it->Centroid(cnt.data());
				coord vel = coord(it->RealArray(velocity_tag).data());
				Storage::real l = (cnt - pnt).length() + eps;
				Storage::real omega = 1.0 / (l*l);
				ret += vel*omega;
				dist += omega;
			}
			ret /= dist;
		}
	}


	Storage::real GetSize(Cell c)
	{
		Storage::real bounds[3][2] = { { 1.0e20, -1.0e20 }, { 1.0e20, -1.0e20 }, { 1.0e20, -1.0e20 } };
		ElementArray<Node> nodes = c->getNodes();
		for (ElementArray<Node>::iterator n = nodes.begin(); n != nodes.end(); ++n)
		{
			Storage::real_array cnt = n->Coords();
			for (int k = 0; k < 3; ++k)
			{
				if (bounds[k][0] > cnt[k]) bounds[k][0] = cnt[k];
				if (bounds[k][1] < cnt[k]) bounds[k][1] = cnt[k];
			}
		}
		Storage::real ret = 1.0e+20;
		for (int k = 0; k < 3; ++k) if (bounds[k][1] - bounds[k][0])
			ret = std::min(ret, bounds[k][1] - bounds[k][0]);
		if (ret > 1.0e+19) std::cout << __FILE__ << ":" << __LINE__ << " oops" << std::endl;
		return ret;
	}
	
	void GetBbox(Element c, Storage::real bounds[3][2])
	{
		bounds[0][0] = 1.0e20;
		bounds[0][1] = -1.0e20;
		bounds[1][0] = 1.0e20;
		bounds[1][1] = -1.0e20;
		bounds[2][0] = 1.0e20;
		bounds[2][1] = -1.0e20;
		//bounds[3][2] = { { 1.0e20, -1.0e20 }, { 1.0e20, -1.0e20 }, { 1.0e20, -1.0e20 } };
		ElementArray<Node> nodes = c->getNodes();
		for (ElementArray<Node>::iterator n = nodes.begin(); n != nodes.end(); ++n)
		{
			Storage::real_array cnt = n->Coords();
			for (int k = 0; k < 3; ++k)
			{
				if (bounds[k][0] > cnt[k]) bounds[k][0] = cnt[k];
				if (bounds[k][1] < cnt[k]) bounds[k][1] = cnt[k];
			}
		}
		return;
		//Storage::real ret = 1.0e+20;
		//for (int k = 0; k < 3; ++k) if (bounds[k][1] - bounds[k][0])
		//	ret = std::min(ret, bounds[k][1] - bounds[k][0]);
		//if (ret > 1.0e+19) std::cout << __FILE__ << ":" << __LINE__ << " oops" << std::endl;
		//return ret;
	}

	INMOST_DATA_REAL_TYPE GetSizeProj(INMOST_DATA_REAL_TYPE bounds[3][2], const coord & v)
	{
		INMOST_DATA_REAL_TYPE vl = v.length();
		if( !vl ) return 0;
		coord uv = v/vl;
		INMOST_DATA_REAL_TYPE ret = 0;
		for(int k = 0; k < 3; ++k) ret += fabs(uv[k]*(bounds[k][1]-bounds[k][0]));
		return ret;
	}
	
	INMOST_DATA_REAL_TYPE GetSizeProj(Element c, const coord & v)
	{
		INMOST_DATA_REAL_TYPE bounds[3][2];
		GetBbox(c,bounds);
		return GetSizeProj(bounds,v);
	}


	Storage::real GetSize(Element n, const Tag & size_tag)
	{
		ElementArray<Cell> cells = n->getCells();
		Storage::real minsize = 1.0e+20, size;
		for (ElementArray<Cell>::iterator c = cells.begin(); c != cells.end(); ++c)
		{
			size = c->RealDF(size_tag);
			if (minsize > size) minsize = size;
		}
		if (minsize > 1.0e+19) std::cout << __FILE__ << ":" << __LINE__ << " oops" << std::endl;
		return minsize;
	}


	void BuildStreamlines(Mesh *mesh, Tag vel, ElementType vel_def, std::vector<Streamline> & output)
	{
		if (!vel.isDefined(vel_def))
		{
			std::cout << __FILE__ << ":" << __LINE__ << " Velocity was not defined on " << ElementTypeName(vel_def) << std::endl;
			return;
		}
		if (vel.GetDataType() != DATA_REAL)
		{
			std::cout << __FILE__ << ":" << __LINE__ << " Data type for velocity is not floating point" << std::endl;
			return;
		}
		if (vel.GetSize() != 3)
		{
			std::cout << __FILE__ << ":" << __LINE__ << " Expected 3 entries in velocity field for streamlines" << std::endl;
			return;
		}
		printf("preparing octree around mesh, was sets %d\n", (int)mesh->NumberOfSets());
		//~ SearchKDTree octsearch(mesh);
		Octree octsearch = Octree(mesh->CreateSet("octsearch").first);
		octsearch.Construct(vel_def, false); //auto-detect octree or quadtree
		printf("done, sets %d\n", (int)mesh->NumberOfSets());
		printf("building streamlines\n");
		Tag cell_size = mesh->CreateTag("STREAMLINES_TEMPORARY_CELL_SIZES", DATA_REAL, CELL, NONE, 1);
		Storage::real velmax = 0, velmin = 1.0e20, l;
		for (Mesh::iteratorCell c = mesh->BeginCell(); c != mesh->EndCell(); ++c)
		{
			coord velv(c->RealArray(vel).data());
			l = velv.length();
			if (l > velmax) velmax = l;
			if (l < velmin) velmin = l;
			c->RealDF(cell_size) = GetSize(c->self());
		}
		velmax = log(velmax + 1.0e-25);
		velmin = log(velmin + 1.0e-25);
		
		{
			MarkerType visited = mesh->CreateMarker();
			int tot = 0;
			int k = 0;
			
			printf("started building streamlines from boundary elements\n");
			
			for (Mesh::iteratorElement f = mesh->BeginElement(vel_def); f != mesh->EndElement(); ++f) if (f->Boundary())
				tot++;

			printf("total elements: %d\n",tot);
			
			for (Mesh::iteratorElement f = mesh->BeginElement(vel_def); f != mesh->EndElement(); ++f) if (f->Boundary())
			{
				if( k % 5 != 0 ) //every tenth face!
				{
					k++;
					continue;
				}
				coord v(f->RealArray(vel).data());
				if (v.length() > 1.0e-4)
				{
					coord nrm(0,0,0);
					if( f->GetElementType() == FACE )
						f->getAsFace().UnitNormal(nrm.data());
					else
					{
						coord ncur;
						int nn = 0;
						ElementArray<Face> faces = f->getFaces();
						for (ElementArray<Face>::iterator ff = faces.begin(); ff != faces.end(); ++ff) if (ff->Boundary())
						{
							ff->UnitNormal(ncur.data());
							nrm+=ncur;
							nn++;
						}
						if( nn )
							nrm /= (double)nn;
						else
							std::cout << __FILE__ << ":" << __LINE__ << " No boundary faces around boundary element" << std::endl;
					}
					double dir = nrm^v;
					if( dir > 0.0 ) //velocity points out of mesh
						dir = -1; 
					else dir = 1; //velocity points into mesh
					coord cntf;
					f->Centroid(cntf.data());
					output.push_back(Streamline(octsearch, cntf, vel, vel_def, cell_size, velmin, velmax, dir, visited));

					ElementArray<Edge> edges = f->getEdges();
					for (ElementArray<Edge>::iterator n = edges.begin(); n != edges.end(); ++n)
					{
						coord cntn;
						n->Centroid(cntn.data());
						if (cntn[2] == cntf[2])
						{
							const Storage::real coef[4] = { 0.4, 0.8 };
							for (int q = 0; q < 2; ++q)
								output.push_back(Streamline(octsearch, cntf*coef[q] + cntn*(1 - coef[q]), vel, vel_def, cell_size, velmin, velmax, 1.0, visited));
						}
					}
				}
				k++;
				if (k % 100 == 0)
				{
					printf("%5.2f%%\r", (double)k / (double)tot*100.0);
					fflush(stdout);
				}
			}

			printf("done from boundary faces, total streamlines = %lu\n", output.size());
			 
			
			printf("started building streamlines from unvisited cells\n");
			tot = 0;
			for (Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
				if (!it->GetMarker(visited)) tot++;
			/*

			printf("total elements: %d\n", tot);
			k = 0;
			for (Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
			{
				if (!it->GetMarker(visited))
				{
					coord cntc;
					it->Centroid(cntc.data());
					if (coord(it->RealArray(vel).data()).length() > 1.0e-4)
					{
						output.push_back(Streamline(octsearch, cntc, vel, vel_def, cell_size, velmin, velmax, 1.0, visited));
						output.push_back(Streamline(octsearch, cntc, vel, vel_def, cell_size, velmin, velmax, -1.0, visited));
					}
				}
				k++;
				if (k % 100 == 0)
				{
					printf("%5.2f%%\r", (double)k / (double)tot*100.0);
					fflush(stdout);
				}
			}
			printf("done from unvisited cells, total streamlines = %lu\n", output.size());
			*/
			mesh->ReleaseMarker(visited,vel_def);
			

		}
		mesh->DeleteTag(cell_size);
		printf("done, total streamlines = %lu\n", output.size());
		printf("killing octree, was sets %d\n", (int)mesh->NumberOfSets());
		octsearch.Destroy();
		printf("done, sets %d\n", (int)mesh->NumberOfSets());

	}


	Streamline::Streamline(const Octree & octsearch, coord pos, Tag velocity_tag, ElementType velocity_defined, Tag cell_size, Storage::real velocity_min, Storage::real velocity_max, Storage::real sign, MarkerType visited)
	{
		(void)cell_size;
		Storage::real coef, len, size;
		coord next = pos, vel;
		Element c;
		const int maxsteps = 25000;
		points.reserve(maxsteps / 2);
		velarr.reserve(maxsteps / 2);
		points.push_back(pos);
		velarr.push_back(0);
		while (points.size() < maxsteps)
		{
			c = octsearch.FindClosestCell(next.data());
			//~ c = octsearch.SearchCell(next.data());
			if (!c.isValid()) break;
			//if( !c.getAsCell().Inside(next.data()) ) break;
			//check we are inside mesh
			/*
			ElementArray<Cell> cells = c->BridgeAdjacencies2Cell(NODE);
			bool inside = false;
			for (ElementArray<Cell>::iterator it = cells.begin(); it != cells.end() && !inside; ++it)
				if( it->Inside(next.data()) ) inside = true;
			if( !inside ) break;
			*/
			c.SetMarker(visited);
			GetVelocity(c, velocity_tag, velocity_defined, next, vel);
			len = vel.length();
			if (len < 1.0e-8) break;
			//size = GetSize(c, cell_size);// c->RealDF(cell_size);
			size = GetSizeProj(c,vel);
			coef = 0.05*size / len;
			next += vel*coef*sign;
			points.push_back(next);
			velarr.push_back((log(len + 1.0e-25) - velocity_min) / (velocity_max - velocity_min));
		}
		//printf("%ld %ld\n",points.size(),velarr.size());
	}


	GLUquadric * cylqs = NULL;
	void drawcylinder(coord a, coord b, double width)
	{
		INMOST_DATA_REAL_TYPE matrix[16];
		if (cylqs == NULL)
		{
			cylqs = gluNewQuadric();
			gluQuadricNormals(cylqs, GLU_SMOOTH);
			gluQuadricOrientation(cylqs, GLU_OUTSIDE);
			gluQuadricDrawStyle(cylqs, GLU_FILL);//GLU_SILHOUETTE
		}
		glPushMatrix();
		glTranslated(a[0], a[1], a[2]);
		get_matrix(a, b, matrix);
#if defined(USE_FP64)
		glMultMatrixd(matrix);
#else
		glMultMatrixf(matrix);
#endif
		gluCylinder(cylqs, width, width, sqrt((b - a) ^ (b - a)), 4, 1);
		glPopMatrix();
	}



	void Streamline::Draw(int reduced)
	{
		if (reduced)
		{
			glBegin(GL_LINE_STRIP);
			for (unsigned int i = 0; i < points.size() - 1; i++)
			{
				glColor3f(velarr[i + 1] * 0.65, 0.65*(velarr[i + 1] < 0.5 ? velarr[i] : 1.0 - velarr[i]), 0.65*(1 - velarr[i + 1]));
				glVertex3d(points[i][0], points[i][1], points[i][2]);
			}
			glEnd();
		}
		else for (unsigned int i = 0; i < points.size() - 1; i++)
		{
			glColor3f(velarr[i + 1] * 0.65, 0.65*(velarr[i + 1] < 0.5 ? velarr[i] : 1.0 - velarr[i]), 0.65*(1 - velarr[i + 1]));
			drawcylinder(points[i], points[i + 1], 0.2*abs(points[i + 1] - points[i]));
		}
	}


	void Streamline::SVGDraw(std::ostream & file, double modelview[16], double projection[16], int viewport[4])
	{
		for (unsigned int i = 0; i < points.size() - 1; i++)
		{
			INMOST_DATA_REAL_TYPE * v0 = points[i].data();
			INMOST_DATA_REAL_TYPE * v1 = points[i + 1].data();
			color_t c(velarr[i + 1] * 0.65, 0.65*(velarr[i + 1] < 0.5 ? velarr[i] : 1.0 - velarr[i]), 0.65*(1 - velarr[i + 1]));
			file << "<g stroke=\"" << c.svg_rgb() << "\">" << std::endl;
			svg_line(file, v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], modelview, projection, viewport);
			file << "</g>" << std::endl;
		}
	}
}
