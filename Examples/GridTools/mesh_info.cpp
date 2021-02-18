#include "inmost.h"

using namespace INMOST;

int main(int argc, char *argv[]) 
{
	
	if( argc < 2 )
	{
		printf("Usage: %s mesh_file [output_file=]\n",argv[0]);
		return -1;
	}
	std::string infile(argv[1]);
	std::cout << " mesh_file=" << infile << std::endl;
	  
	
	Mesh mesh; 
	mesh.Load(infile);
	
	// prepare geometrical data on the mesh
	std::cout << "Prepare geometric data" << std::endl;
	Mesh::GeomParam table;
	table[NORMAL]      = FACE;        //Compute normals
	table[ORIENTATION] = FACE;        //Check and fix normal orientation
	table[MEASURE]     = CELL;        //Compute volumes
	table[BARYCENTER]  = CELL | FACE; //Compute volumetric center of mass
	mesh.RemoveGeometricData(table);
	mesh.PrepareGeometricData(table); //Ask to precompute the data
	std::cout << "Done with geometric data" << std::endl;



	std::cout << "Start mesh info" << std::endl;
	rMatrix I(3, 3), E(3, 3), n(3, 1), x(3, 1), xc(3,1), one(3, 1), g1(3, 1);
	one(0, 0) = 1;
	one(1, 0) = 1;
	one(2, 0) = 1;
	I = rMatrix::Unit(3);
	Storage::integer bad = 0;
	Storage::integer ccout = 0;
	Storage::integer bcout = 0;
	Storage::integer elog = 0;
	Storage::integer bad_div1 = 0, bad_divx = 0, bad_g1 = 0, bad_E = 0;
	Storage::real vmax = -1.0e+20, vmin = 1.0e+20;
	Storage::real divx, div1;
	Storage::integer histogram_E[10], histogram_divx[10], histogram_div1[10], histogram_g1[10];
	memset(histogram_E, 0, sizeof(int) * 10);
	memset(histogram_divx, 0, sizeof(int) * 10);
	memset(histogram_div1, 0, sizeof(int) * 10);
	memset(histogram_g1, 0, sizeof(int) * 10);
	TagInteger tbad_E = mesh.CreateTag("badE", DATA_INTEGER, CELL, NONE, 1);
	TagInteger tbad_C = mesh.CreateTag("badC", DATA_INTEGER, CELL, NONE, 1);
	for (Mesh::iteratorCell it = mesh.BeginCell(); it != mesh.EndCell(); ++it) if (it->GetStatus() != Element::Ghost)
	{
		ElementArray<Face> faces = it->getFaces();
		E.Zero();
		g1.Zero();
		divx = 0;
		div1 = 0;
		Storage::real vol = it->Volume();
		it->Barycenter(xc.data());
		for (ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
		{
			jt->OrientedNormal(it->self(), n.data());
			jt->Barycenter(x.data());
			//jt->Centroid(x.data());
			E += (x-xc) * n.Transpose();
			div1 += n.DotProduct(one);
			divx += n.DotProduct(x);
			g1 += n;
		}
		div1 = fabs(div1 / vol);
		divx = fabs(divx / vol - 3);
		//~ bool print = false;
		g1 /= vol;
		E = (E / vol - I);
		if (E.FrobeniusNorm() > 1.0e-4)
		{
			tbad_E[*it] = 1;
			bad_E++;
		}
		elog = -ceil(log(E.FrobeniusNorm() + 1.0e-33) / log(10));
		if (elog > 9) elog = 9;
		if (elog < 0) elog = 0;
		//~ if (elog == 0) print = true;
		histogram_E[elog]++;
		if (divx > 1.0e-4) bad_divx++;
		elog = -ceil(log(fabs(divx) + 1.0e-33) / log(10));
		if (elog > 9) elog = 9;
		if (elog < 0) elog = 0;
		//~ if (elog == 0) print = true;
		histogram_divx[elog]++;
		if (div1 > 1.0e-4) bad_div1++;
		elog = -ceil(log(fabs(div1) + 1.0e-33) / log(10));
		if (elog > 9) elog = 9;
		if (elog < 0) elog = 0;
		//~ if (elog == 0) print = true;
		histogram_div1[elog]++;
		if (g1.FrobeniusNorm() > 1.0e-4) bad_g1++;
		elog = -ceil(log(fabs(g1.FrobeniusNorm() + 1.0e-33)) / log(10));
		if (elog > 9) elog = 9;
		if (elog < 0) elog = 0;
		//~ if (elog == 0) print = true;
		histogram_g1[elog]++;
		//if( print )
		//	std::cout << std::setw(14) << E.FrobeniusNorm() << std::setw(14)  << divx << std::setw(14) <<  div1 << std::setw(14) <<  g1.FrobeniusNorm() << std::endl;
		if (E.FrobeniusNorm() > 1.0e-4 || vol < 0)
		{
			//std::cout << "Cell " << it->LocalID() << std::endl;
			//E.Print();
			bad++;
		}
		it->Barycenter(x.data());
		if (!it->Inside(x.data())) bcout++;
		it->Centroid(x.data());
		if (!it->Inside(x.data()))
		{
			tbad_C[*it] = 1;
			ccout++;
		}

		if (vmax < vol) vmax = vol;
		if (vmin > vol) vmin = vol;

	}
	bad_E = mesh.Integrate(bad_E);
	bad_g1 = mesh.Integrate(bad_g1);
	bad_divx = mesh.Integrate(bad_divx);
	bad_div1 = mesh.Integrate(bad_div1);
	mesh.Integrate(histogram_E,10);
	mesh.Integrate(histogram_g1,10);
	mesh.Integrate(histogram_divx,10);
	mesh.Integrate(histogram_div1,10);
	if( mesh.GetProcessorRank() == 0 )
	{
		std::cout << std::setw(4) << "k" << std::setw(10) << "Gx" << std::setw(10) << "G1" << std::setw(10) << "Divx" << std::setw(10) << "Div1" << std::endl;
		std::cout << std::setw(4) << " " << std::setw(10) << bad_E << std::setw(10) << bad_g1 << std::setw(10) << bad_divx << std::setw(10) << bad_div1 << std::endl;
		for (int k = 0; k < 10; ++k)
			std::cout << std::setw(4) << k << std::setw(10) << histogram_E[k] << std::setw(10) << histogram_g1[k] << std::setw(10) << histogram_divx[k] << std::setw(10) << histogram_div1[k] << std::endl;
	}
	bad = mesh.Integrate(bad);
	ccout = mesh.Integrate(ccout);
	int ncells = mesh.TotalNumberOf(CELL);
	vmin = mesh.AggregateMin(vmin);
	vmax = mesh.AggregateMax(vmax);
	if (mesh.GetProcessorRank() == 0)
		std::cout << "Cells with bad divergence: " << bad << "/" << ncells << " centroid outside " << ccout << " barycenter outside " << bcout << " volume " << vmin << ":" << vmax << std::endl;
	if( !bad_E )
		tbad_E = mesh.DeleteTag(tbad_E);
	if( !ccout )
		tbad_C = mesh.DeleteTag(tbad_C);
	Storage::integer npln = 0;
	Storage::integer ornt = 0;
	Storage::integer nfaces = mesh.TotalNumberOf(FACE);
	for (Mesh::iteratorFace it = mesh.BeginFace(); it != mesh.EndFace(); ++it)
	{
		if (!it->CheckNormalOrientation()) ornt++;
		if (!it->Planarity()) npln++;
	}
	ornt = mesh.Integrate(ornt);
	npln = mesh.Integrate(npln);
	if (mesh.GetProcessorRank() == 0)
		std::cout << "Bad face orientation " << ornt << "/" << nfaces << " no planarity " << npln << std::endl;
	Storage::real cmax[3] = { -1.0e20,-1.0e20,-1.0e20 };
	Storage::real cmin[3] = { 1.0e20, 1.0e20, 1.0e20 };
	for (Mesh::iteratorNode it = mesh.BeginNode(); it != mesh.EndNode(); ++it)
	{
		for (int k = 0; k < 3; ++k)
		{
			if (cmax[k] < it->Coords()[k]) cmax[k] = it->Coords()[k];
			if (cmin[k] > it->Coords()[k]) cmin[k] = it->Coords()[k];
		}
	}
	mesh.AggregateMax(cmax,3);
	mesh.AggregateMin(cmin,3);
	if( mesh.GetProcessorRank() == 0 )
	{
		std::cout << "x " << cmin[0] << ":" << cmax[0] << " ";
		std::cout << "y " << cmin[1] << ":" << cmax[1] << " ";
		std::cout << "z " << cmin[2] << ":" << cmax[2] << " ";
		std::cout << std::endl;
	}
	//mesh.Save("geom_info.vtk");
	
	Storage::real vol = 0;
	for (Mesh::iteratorCell it = mesh.BeginCell(); it != mesh.EndCell(); ++it) if (it->GetStatus() != Element::Ghost)
		vol += it->Volume();
	vol = mesh.Integrate(vol);
	if( mesh.GetProcessorRank() == 0 )
	{
		std::cout << "mesh volume " << vol  << std::endl;
	}
	std::cout << "Done with mesh info" << std::endl;
	if( argc > 2 )
	{
		std::cout << "Save marked mesh to " << std::string(argv[2]) << std::endl;
		mesh.Save(std::string(argv[2]));
	}
	return 0;
}
