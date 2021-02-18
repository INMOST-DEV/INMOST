#include "inmost.h"

using namespace INMOST;

//todo: want to separate all faces into those that share edges
//see todo: inside text


double dot(const double a[3], const double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}



struct coords
{
	double v[3];
	int compare(const double bv[3]) const
	{
		const double eps = 1.0e-11;//pow(10,-14);
		for (int i = 0; i < 3; i++)
		{
			if (fabs(v[i] - bv[i]) > eps)
			{
				if (v[i] > bv[i])
					return 1;
				else
					return -1;
			}
		}
		return 0;
	}
	bool operator <(const coords & b) const { return compare(b.v) < 0; }
	coords(double _v[3]) { memcpy(v, _v, sizeof(double)* 3); }
	coords(const coords &b) { memcpy(v, b.v, sizeof(double)* 3); }
	coords & operator =(coords const & b) { memmove(v, b.v, sizeof(double)* 3); return *this; }
};





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
	//There seems to be a problem with the mesh
	int fixed = 0;
	for (Mesh::iteratorFace f = m.BeginFace(); f != m.EndFace(); ++f)
	{
		if (f->FixNormalOrientation())
		{
			//std::cout << "face " << f->LocalID() << " nodes " << f->nbAdjElements(NODE);
			//if( f->BackCell().isValid() ) std::cout << " back cell " << f->BackCell().LocalID() << " " << Element::GeometricTypeName(f->BackCell()->GetGeometricType()) << " nodes " << f->BackCell()->nbAdjElements(NODE) << " ";
			//if (f->FrontCell().isValid()) std::cout << " front cell " << f->FrontCell().LocalID() << " " << Element::GeometricTypeName(f->FrontCell()->GetGeometricType()) << " nodes " << f->FrontCell()->nbAdjElements(NODE) << " ";
			//std::cout << std::endl;
			fixed++;
		}
	}
	std::cout << "Time to fix normals: " << Timer() - tt << std::endl;
	std::cout << "Total face normals fixed: " << fixed << std::endl;
	
	//fix multiple flat faces
	{

		tt = Timer();
		m.BeginModification();
		MarkerType shareedge = m.CreateMarker();
		TagInteger united = m.CreateTag("united", DATA_INTEGER, FACE, NONE, 1);
		for (Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		{
			ElementArray<Face> faces = it->getFaces();
			//find out faces that share the same adjacent cells
			std::map<Cell, std::vector<HandleType> > cfaces;
			for (ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
			{
				Cell cK = jt->BackCell() == it->self() ? jt->FrontCell() : jt->BackCell();
				cfaces[cK].push_back(jt->GetHandle());
			}
			for (std::map<Cell, std::vector<HandleType> >::iterator jt = cfaces.begin(); jt != cfaces.end(); ++jt)
			{
				if (jt->second.size() > 1)
				{
					//group faces with the same normal
					std::vector< std::pair<coords, std::vector<HandleType> > > nfaces;
					for (std::vector<HandleType>::iterator kt = jt->second.begin(); kt != jt->second.end(); ++kt)
					{
						Face f(&m, *kt);
						double nrm[3];
						f->OrientedUnitNormal(it->self(), nrm);
						bool added = false;
						for (std::vector< std::pair< coords, std::vector<HandleType> > >::iterator qt = nfaces.begin(); qt != nfaces.end(); ++qt)
						{
							double * v = qt->first.v;
							if (dot(nrm, v) > 0.95)
							{
								for (int k = 0; k < 3; ++k)
									v[k] = (v[k] + nrm[k])*0.5;
								qt->second.push_back(*kt);
								added = true;
								break;
							}
						}
						if (!added) nfaces.push_back(std::make_pair(nrm, std::vector<HandleType>(1, *kt)));
					}
					//unite faces in the group
					for (std::vector< std::pair< coords, std::vector<HandleType> > >::iterator kt = nfaces.begin(); kt != nfaces.end(); ++kt)
					{
						if (kt->second.size() > 1)
						{
							//check that each face shares an edge with any other
							bool all_share = true;
							for (std::vector<HandleType>::iterator qt = kt->second.begin(); qt != kt->second.end(); ++qt)
							{
								bool qt_share = false;
								Face f(&m, *qt);
								ElementArray<Edge> fedges = f.getEdges();
								fedges.SetMarker(shareedge);
								for (std::vector<HandleType>::iterator vt = kt->second.begin(); vt != kt->second.end() && !qt_share; ++vt)
								if (Face(&m, *vt).nbAdjElements(EDGE, shareedge))
									qt_share = true;
								fedges.RemMarker(shareedge);
								all_share &= qt_share;
							}
							//todo: want to separate all faces into those that share edges
							//start from one face then add one that shares any edge
							//once there are no more faces to add, create another set and
							//add faces to it one by one
							if (all_share)
							{
								Face f = Face::UniteFaces(ElementArray<Face>(&m, kt->second.begin(), kt->second.end()), 0);
								united[f] = 1;
							}
						}
					}
				}
			}
		}
		m.ReleaseMarker(shareedge);
		m.EndModification();
		std::cout << "Time to unite faces:" << Timer() - tt << std::endl;
		std::cout << "Unite flat:" << std::endl;
		std::cout << "Cells: " << m.NumberOfCells() << std::endl;
		std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
		std::cout << "Edges: " << m.NumberOfEdges() << std::endl;
	}

	m.Save(grid_out);
	return 0;
}
