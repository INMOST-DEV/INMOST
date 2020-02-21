#include "fracture.h"


INMOST_DATA_REAL_TYPE Fracture::Volume(Cell c) const
{
	if( c->GetMarker(isFracture()) )
		return c->Real(fracture_volume);
	else
		return c->Volume();
}

INMOST_DATA_REAL_TYPE Fracture::Area(Face f) const
{
	if( f->GetMarker(isFracture()) )
		return f->Real(fracture_area);
	else
		return f->Area();
}

void Fracture::FaceCenter(Face f, INMOST_DATA_REAL_TYPE cnt[3]) const
{
	if( f->GetMarker(matrix_fracture) )
	{
		Cell cK = f->BackCell();
		Cell cL = f->FrontCell();
		bool fcK = cK->GetMarker(isFracture());
		bool fcL = cL->GetMarker(isFracture());
		if( fcK && !fcL ) //cK is fracture cell
		{
			Storage::real cntK[3], nrm[3];
			cK->Centroid(cntK);
			f->OrientedUnitNormal(cK,nrm);
			Storage::real half_aperture = cK->Real(fracture_aperture)*0.5;
			cnt[0] = cntK[0] + nrm[0]*half_aperture;
			cnt[1] = cntK[1] + nrm[1]*half_aperture;
			cnt[2] = cntK[2] + nrm[2]*half_aperture;
		}
		else if( fcL && !fcK ) // cL is fracture cell
		{
			Storage::real cntL[3], nrm[3];
			cL->Centroid(cntL);
			f->OrientedUnitNormal(cL,nrm);
			Storage::real half_aperture = cL->Real(fracture_aperture)*0.5;
			cnt[0] = cntL[0] + nrm[0]*half_aperture;
			cnt[1] = cntL[1] + nrm[1]*half_aperture;
			cnt[2] = cntL[2] + nrm[2]*half_aperture;
		}
		else f->Centroid(cnt);
			}
	else if( f->GetMarker(multiple_fracture_joints) )
	{
		ElementArray<Edge> edges = f->getEdges(isFracture());
		assert(edges.size() == 1); //there must be only one such edge
		Storage::real cntf[3], cnte[3], vec[3],l;
		f->Centroid(cntf);
		edges[0]->Centroid(cnte);
		vec[0] = cntf[0] - cnte[0];
		vec[1] = cntf[1] - cnte[1];
		vec[2] = cntf[2] - cnte[2];
		l = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
		if( l )
		{
			vec[0] /= l;
			vec[1] /= l;
			vec[2] /= l;
		}
		Storage::real half_aperture = f->Real(fracture_aperture)*0.5;
		cnt[0] = cnte[0];// + vec[0]*half_aperture;
		cnt[1] = cnte[1];// + vec[1]*half_aperture;
		cnt[2] = cnte[2];// + vec[2]*half_aperture;
	}
	else f->Centroid(cnt);
}

void Fracture::Open(Tag aperture, bool fill_fracture, double gap_multiplier)
{
	fracture_aperture = aperture;
	std::cout << "Cells: " << m->NumberOfCells() << std::endl;
	std::cout << "Faces: " << m->NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m->NumberOfEdges() << std::endl;
	std::cout << "Nodes: " << m->NumberOfNodes() << std::endl;
	m->BeginModification();
	fracture_marker = m->CreateMarker();
	std::cout << "create marker fracture_marker " << fracture_marker << std::endl;
	for(Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
		if( f->GetMarker(fracture_marker) ) std::cout << "Face " << f->LocalID() << " already marked fracture " << fracture_marker << std::endl;
	for(Mesh::iteratorCell f = m->BeginCell(); f != m->EndCell(); ++f)
		if( f->GetMarker(fracture_marker) ) std::cout << "Cell " << f->LocalID() << " already marked fracture " << fracture_marker << std::endl;
	for(Mesh::iteratorFace f = m->BeginFace(); f != m->EndFace(); ++f)
		if( f->HaveData(aperture) ) f->SetMarker(fracture_marker);
	m->self()->Integer(m->CreateTag("FRACTURE_MARKER",DATA_INTEGER,MESH,NONE,1)) = fracture_marker;
	//break up all the faces with fractures into volumes?
	fracture_volume = m->CreateTag("VOLUME_FRACTURE",DATA_REAL,CELL,CELL,1);
	const Storage::real mult = gap_multiplier;
	Tag indexes = m->CreateTag("FRAC_TEMP_INDEXES",DATA_INTEGER,CELL,NONE,2);
	//associate a node or a pair of nodes to each non-fracture edge that ends with a fracture node
	Tag cell2node = m->CreateTag("CELL2NODE",DATA_REFERENCE,CELL,NONE);
	//new faces to reconnect
	Tag connface = m->CreateTag("CONNFACE",DATA_REFERENCE,CELL,NONE);
	//original position of the fracture nodes
	Tag origcoords = m->CreateTag("ORIGINAL_FRACTURE_COORDS",DATA_REAL,NODE,NODE,3);
	//std::cout << "Initial number of nodes: " << m->NumberOfNodes() << std::endl;
	
	std::cout << "mult is " << mult << std::endl;

	//Unmark any boundary faces
	int unmrk = 0, mrk = 0;
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it)
	{
		if( it->GetMarker(isFracture()) ) mrk++;
		if( !it->FrontCell().isValid() && it->GetMarker(isFracture()) )
		{
			it->RemMarker(isFracture());
			unmrk++;
		}
	}
	
	std::cout << "fracture marked " << mrk << " boundary unmarked " << unmrk << std::endl;

	std::vector<Tag> transfer_node_real_tags;
	std::vector<Tag> transfer_node_integer_tags;
	for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); ++t)
		if( t->isDefined(NODE) && *t != m->CoordsTag() )
		{
			if( t->GetDataType() == DATA_REAL )
				transfer_node_real_tags.push_back(*t);
			else if( t->GetDataType() == DATA_INTEGER )
				transfer_node_integer_tags.push_back(*t);
		}
	std::cout << "Transfer real tags:" << std::endl;
	for(int q = 0; q < (int)transfer_node_real_tags.size(); ++q)
		std::cout << transfer_node_real_tags[q].GetTagName() << std::endl;
	std::cout << "Transfer integer tags:" << std::endl;
	for(int q = 0; q < (int)transfer_node_integer_tags.size(); ++q)
		std::cout << transfer_node_integer_tags[q].GetTagName() << std::endl;
	//For each node create several new ones that will be used in open fracture
	std::cout << "Nodes in mesh " << m->NumberOfNodes() << std::endl;
	int nfrac_nodes =0, new_nodes=0;
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		ElementArray<Face> n_faces = it->getFaces(); //retrive all faces joining at the node
		INMOST_DATA_ENUM_TYPE num_frac_faces = 0;
		for(int k = 0; k < n_faces.size(); ++k) if( n_faces[k].GetMarker(isFracture()) ) num_frac_faces++;
		if( num_frac_faces )
		{
			ElementArray<Cell> n_cells = it->getCells(); //get all cells of the node
			for(int k = 0; k < n_cells.size(); ++k) //assign local enumeration to the cells
				n_cells[k].IntegerDF(indexes) = k;
			for(int k = 0; k < n_faces.size(); ++k) //stich together node's numbers corresponding to cells if no fracture separates them
			{
				if( !n_faces[k].GetMarker(isFracture()) && n_faces[k].FrontCell().isValid())
				{
					int bi = n_faces[k].BackCell()->IntegerDF(indexes);
					int fi = n_faces[k].FrontCell()->IntegerDF(indexes);
					if( bi != fi )
					{
						for(int q = 0; q < n_cells.size(); q++) if( n_cells[q].IntegerDF(indexes) == fi ) 
							n_cells[q].IntegerDF(indexes) = bi;
					}
				}
			}
			dynarray<int,64> nums(n_cells.size()); //store all numbers
			for(int k = 0; k < n_cells.size(); ++k) nums[k] = n_cells[k].IntegerDF(indexes);
			std::sort(nums.begin(),nums.end());
			nums.resize(std::unique(nums.begin(),nums.end()) - nums.begin());
			if( nums.size() > 1 ) //at least two distinctive nodes
			{
				nfrac_nodes++;
				for(int k = 0; k < nums.size(); k++)
				{
					Node image;
					Storage::real xyz[3] = {0,0,0}, cntc[3] = {0,0,0};
					int n_node_cells = 0;
					for(int q = 0; q < n_cells.size(); q++)
					{
						if( n_cells[q].IntegerDF(indexes) == nums[k] )
						{
							n_cells[q].Centroid(cntc);
							xyz[0] += cntc[0];
							xyz[1] += cntc[1];
							xyz[2] += cntc[2];
							n_node_cells++;
						}
					}
					xyz[0] /= static_cast<Storage::real>(n_node_cells);
					xyz[1] /= static_cast<Storage::real>(n_node_cells);
					xyz[2] /= static_cast<Storage::real>(n_node_cells);
					Storage::real_array coords = it->Coords();
					for( Storage::real_array::size_type q = 0; q < coords.size(); ++q)
						xyz[q] = coords[q]*mult + (1-mult)*xyz[q];
					image = m->CreateNode(xyz);
					new_nodes++;
					Storage::real_array save_coords = image->RealArray(origcoords);
					save_coords[0] = coords[0];
					save_coords[1] = coords[1];
					save_coords[2] = coords[2];
					for(int q = 0; q < (int)transfer_node_real_tags.size(); ++q)
						if( it->HaveData(transfer_node_real_tags[q]) )
						{
							Storage::real_array source = it->RealArray(transfer_node_real_tags[q]);
							Storage::real_array target = image->RealArray(transfer_node_real_tags[q]);
							if( target.size() != source.size() ) target.resize(source.size());
							for(int qq = 0; qq < source.size(); ++qq)
								target[qq] = source[qq];
						}
					
					for(int q = 0; q < (int)transfer_node_integer_tags.size(); ++q)
						if( it->HaveData(transfer_node_integer_tags[q]) )
						{
							Storage::integer_array source = it->IntegerArray(transfer_node_integer_tags[q]);
							Storage::integer_array target = image->IntegerArray(transfer_node_integer_tags[q]);
							if( target.size() != source.size() ) target.resize(source.size());
							for(int qq = 0; qq < source.size(); ++qq)
								target[qq] = source[qq];
						}
					
					for(int q = 0; q < n_cells.size(); q++)
					{
						if( n_cells[q].IntegerDF(indexes) == nums[k] )
						{
							Storage::reference_array fnodes = n_cells[q].ReferenceArray(cell2node);
							fnodes.push_back(*it); //TODO
							fnodes.push_back(image);
						}
					}
				}
				it->SetMarker(isFracture()); //mark node as fracture node
				//std::cout << "Number of nodes: " << m->NumberOfNodes() << std::endl;
			}
			
		}
	}
	std::cout << "Nodes in mesh " << m->NumberOfNodes() << " frac nodes " << nfrac_nodes << " new nodes " << new_nodes << std::endl;
	std::cout << "Cells: " << m->NumberOfCells() << std::endl;
	std::cout << "Faces: " << m->NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m->NumberOfEdges() << std::endl;
	std::cout << "Nodes: " << m->NumberOfNodes() << std::endl;
	m->DeleteTag(indexes);

	//this tag is used to transfer position of fracture edge center onto fracture face center
	Tag centroid_tag;
	if( m->HaveTag("GEOM_UTIL_CENTROID") )
		centroid_tag = m->GetTag("GEOM_UTIL_CENTROID");
	if( centroid_tag.isValid() && !centroid_tag.isDefined(FACE) )
		centroid_tag = Tag();
	////////
	std::vector<Tag> transfer_face_real_tags;
	std::vector<Tag> transfer_face_integer_tags;
	for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); ++t)
		if( t->isDefined(FACE) )
		{
			if( t->GetTagName().substr(0,10) == "GEOM_UTIL_" ) continue;
			if( t->GetDataType() == DATA_REAL )
				transfer_face_real_tags.push_back(*t);
			else if( t->GetDataType() == DATA_INTEGER )
				transfer_face_integer_tags.push_back(*t);
		}
	Tag fracture_edge_length;
	if( m->GetDimensions() == 3 )
	{
		fracture_edge_length = m->CreateTag("LENGTH_FRACTURE",DATA_REAL,FACE,FACE,1);
		//Mark edges that will be incorporated into fracture cells
		int nedges = 0;
		for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it)
		{
			if( it->nbAdjElements(FACE,isFracture()) > 2) 
			{
				it->SetMarker(isFracture());
				nedges++;
			}
		}
		
		std::cout << "fracture edges with multiple joints " << nedges << std::endl;
		
		
		

		//create new matrix-matrix faces, adjacent to fractures
		//for all non-fracture faces that have any fracture edge create new face
		int nfaces = 0;
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( !it->GetMarker(isFracture()) && it->nbAdjElements(NODE,isFracture()) )
		{
			ElementArray<Edge> edges = it->getEdges(); //go through edges of current face
			ElementArray<Edge> newedges(m);
			Storage::reference_array images = it->BackCell()->ReferenceArray(cell2node);
			//mark me and all my edges that they may contain images
			for(ElementArray<Edge>::iterator jt = edges.begin(); jt != edges.end(); ++jt)
			{
				if( jt->getBeg()->GetMarker(isFracture()) || jt->getEnd()->GetMarker(isFracture()) )
				{
					ElementArray<Node> enodes(m,2);
					if( jt->getBeg()->GetMarker(isFracture()) ) //search image of node between my elements
					{
						for(int k = 0; k < static_cast<int>(images.size()) && !enodes[0].isValid(); k+=2)
						{
							if( images[k] == jt->getBeg() ) // find my element between images
								enodes[0] = images[k+1]->getAsNode();
						}
						assert(enodes[0] != InvalidElement());
					}
					else enodes[0] = jt->getBeg();
					if( jt->getEnd()->GetMarker(isFracture()) )
					{
						for(int k = 0; k < static_cast<int>(images.size()) && !enodes[1].isValid(); k+=2)
						{
							if( images[k] == jt->getEnd() ) // find my element between images
								enodes[1] = images[k+1]->getAsNode();
						}
						assert(enodes[1] != InvalidElement());
					}
					else enodes[1] = jt->getEnd();
					if( enodes[0] != enodes[1] ) newedges.push_back(m->CreateEdge(enodes).first);
				}
				else newedges.push_back(*jt);
			}
			std::pair<Face,bool> f = m->CreateFace(newedges);
			nfaces++;
			if( !f.second ) std::cout << __FILE__ << ":" << __LINE__ << " Face already exists!!! source " << *it << " new " << f.first.GetHandle()  << std::endl;
			else
			{
				for(int q = 0; q < (int)transfer_face_real_tags.size(); ++q)
					if( it->HaveData(transfer_face_real_tags[q]) )
					{
						Storage::real_array source = it->RealArray(transfer_face_real_tags[q]);
						Storage::real_array target = f.first->RealArray(transfer_face_real_tags[q]);
						if( target.size() != source.size() ) target.resize(source.size());
						for(int qq = 0; qq < source.size(); ++qq)
							target[qq] = source[qq];
					}
				for(int q = 0; q < (int)transfer_face_integer_tags.size(); ++q)
				 if( it->HaveData(transfer_face_integer_tags[q]) )
				 {
					 Storage::integer_array source = it->IntegerArray(transfer_face_integer_tags[q]);
					 Storage::integer_array target = f.first->IntegerArray(transfer_face_integer_tags[q]);
					 if( target.size() != source.size() ) target.resize(source.size());
					 for(int qq = 0; qq < source.size(); ++qq)
						 target[qq] = source[qq];
				 }
				 
			}
			//indicate that this new face should be connected to respecting cells
			if(it->BackCell().isValid() ) it->BackCell()->ReferenceArray(connface).push_back(f.first);
			if(it->FrontCell().isValid() ) it->FrontCell()->ReferenceArray(connface).push_back(f.first);
		}
		std::cout << "new fracture faces " << nfaces << std::endl;
		if( !fracture_aperture.isDefined(CELL) ) // this is to extend volume factor
			m->CreateTag(fracture_aperture.GetTagName(),fracture_aperture.GetDataType(),CELL,CELL,fracture_aperture.GetSize());


		std::vector<Tag> transfer_real_tags;
		std::vector<Tag> transfer_integer_tags;
		for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); ++t)
		{
			if( t->isDefined(FACE) && (t->GetDataType() == DATA_REAL || t->GetDataType() == DATA_INTEGER) )
			{
				if( t->GetTagName().substr(0,10) == "GEOM_UTIL_" &&
				    t->GetTagName() != "GEOM_UTIL_CENTROID" &&
				    t->GetTagName() != "GEOM_UTIL_BARYCENTER") continue;
				Tag q;
				if( !t->isDefined(CELL) )
					q = m->CreateTag(t->GetTagName(),t->GetDataType(),CELL,CELL,t->GetSize());
				else 
					q = *t;
				switch(t->GetDataType() )
				{
				case DATA_REAL:    transfer_real_tags.push_back(q); break;
				case DATA_INTEGER: transfer_integer_tags.push_back(q); break;
					default: printf("%s:%d unexpected\n",__FILE__,__LINE__); break;
				}
			}
		}
		
		std::cout << "Copying following real-value tags from fracture faces to fracture cells:" << std::endl;
		for(std::vector<Tag>::iterator q = transfer_real_tags.begin(); q != transfer_real_tags.end(); ++q)
			std::cout << q->GetTagName() << std::endl;

		std::cout << "Copying following integer-value tags from fracture faces to fracture cells:" << std::endl;
		for(std::vector<Tag>::iterator q = transfer_integer_tags.begin(); q != transfer_integer_tags.end(); ++q)
			std::cout << q->GetTagName() << std::endl;
		
		
		
		//now create fracture-fracture control volumes, add fracture-matrix faces to matrix cells
		int kid = 0, new_edges = 0, new_faces = 0, new_cells = 0;
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->GetMarker(isFracture()) && it->nbAdjElements(NODE,isFracture()) )
		{
			//I need all elements adjacent to neighbouring cells of current face
			//and ending at any element of current face
			ElementArray<Face> fracfaces(m);
			ElementArray<Node> nodes = it->getNodes(), nodesb(m,nodes.size()), nodesf(m,nodes.size());
			ElementArray<Edge> edges = it->getEdges(), edgesb(m,nodes.size()), edgesf(m,nodes.size());
			INMOST_DATA_ENUM_TYPE N = (INMOST_DATA_ENUM_TYPE)nodes.size();
			Storage::reference_array images;
			{
				Cell c = it->BackCell();
				images = c->ReferenceArray(cell2node);
				//among images of fracture nodes select those highlighted by marker
				for(int q = 0; q < nodes.size(); ++q) if( nodes[q].GetMarker(isFracture()) )
				{
					Node n;
					for(int k = 0; k < (int)images.size() && !n.isValid(); k+=2)
						if( images[k] == nodes[q] ) n = images[k+1]->getAsNode();
					assert(n.isValid());
					nodesb[q] = n;
				}
				else nodesb[q] = nodes[q];
				//create edges
				for(int k = 0; k < edges.size(); k++)
				{
					if( nodes[k].GetMarker(isFracture()) || nodes[(k+1)%N].GetMarker(isFracture()) )
					{
						ElementArray<Node> enodes(m,2);
						enodes[0] = nodesb[k];
						enodes[1] = nodesb[(k+1)%N];
						edgesb[k] = m->CreateEdge(enodes).first;
						new_edges++;
					}
					else edgesb.data()[k] = edges.data()[k];
				}
				//This is matrix-fracture face
				std::pair<Face,bool> facesb = m->CreateFace(edgesb);
				new_faces++;
				if(!facesb.second ) std::cout << __FILE__ << ":" << __LINE__ << " Face already exists!!! " << *it << " " << facesb.first->GetHandle() << std::endl;
				else
				{
					for(int q = 0; q < (int)transfer_face_real_tags.size(); ++q)
						if( it->HaveData(transfer_face_real_tags[q]) )
						{
							Storage::real_array source = it->RealArray(transfer_face_real_tags[q]);
							Storage::real_array target = facesb.first->RealArray(transfer_face_real_tags[q]);
							if( target.size() != source.size() ) target.resize(source.size());
							for(int qq = 0; qq < source.size(); ++qq)
								target[qq] = source[qq];
						}
					for(int q = 0; q < (int)transfer_face_integer_tags.size(); ++q)
						if( it->HaveData(transfer_face_integer_tags[q]) )
						{
							Storage::integer_array source = it->IntegerArray(transfer_face_integer_tags[q]);
							Storage::integer_array target = facesb.first->IntegerArray(transfer_face_integer_tags[q]);
							if( target.size() != source.size() ) target.resize(source.size());
							for(int qq = 0; qq < source.size(); ++qq)
								target[qq] = source[qq];
						}
				}
				//add faces to indicate reconnection
				it->BackCell()->ReferenceArray(connface).push_back(facesb.first);
				fracfaces.push_back(facesb.first);
			}

			if( it->FrontCell().isValid() ) //there is another cell
			{
				Cell c = it->FrontCell();
				images = c->ReferenceArray(cell2node);
				//among images of fracture nodes select those highlited by marker
				for(int q = 0; q < nodes.size(); ++q) if( nodes[q].GetMarker(isFracture()) )
				{
					Node n;
					for(int k = 0; k < (int)images.size() && !n.isValid(); k+=2)
						if( images[k] == nodes[q] ) n = images[k+1]->getAsNode();
					assert(n.isValid());
					nodesf[q] = n;
				}
				else nodesf[q] = nodes[q];
				//create edges
				for(int k = 0; k < edges.size(); k++)
				{
					if( nodes[k].GetMarker(isFracture()) || nodes[(k+1)%N].GetMarker(isFracture()) )
					{
						ElementArray<Node> enodes(m,2);
						enodes[0] = nodesf[k];
						enodes[1] = nodesf[(k+1)%N];
						edgesf[k] = m->CreateEdge(enodes).first;
						new_edges++;
					}
					else edgesf[k] = edges[k];
				}
				//This is matrix-fracture face
				std::pair<Face,bool> facesf = m->CreateFace(edgesf);
				new_faces++;
				if( !facesf.second ) std::cout << __FILE__ << ":" << __LINE__ << " Face already exists!!! " << *it << " " << facesf.first->GetHandle() << std::endl;
				else
				{
					for(int q = 0; q < (int)transfer_face_real_tags.size(); ++q)
						if( it->HaveData(transfer_face_real_tags[q]) )
						{
							Storage::real_array source = it->RealArray(transfer_face_real_tags[q]);
							Storage::real_array target = facesf.first->RealArray(transfer_face_real_tags[q]);
							if( target.size() != source.size() ) target.resize(source.size());
							for(int qq = 0; qq < source.size(); ++qq)
								target[qq] = source[qq];
						}
					for(int q = 0; q < (int)transfer_face_integer_tags.size(); ++q)
						if( it->HaveData(transfer_face_integer_tags[q]) )
						{
							Storage::integer_array source = it->IntegerArray(transfer_face_integer_tags[q]);
							Storage::integer_array target = facesf.first->IntegerArray(transfer_face_integer_tags[q]);
							if( target.size() != source.size() ) target.resize(source.size());
							for(int qq = 0; qq < source.size(); ++qq)
								target[qq] = source[qq];
						}
				}
				//add faces to indicate reconnection
				it->FrontCell()->ReferenceArray(connface).push_back(facesf.first);
				fracfaces.push_back(facesf.first);
				//now create fracture-fracture faces
				//HERE NODES MAY BE SIMILAR
				ElementArray<Edge> f_edges(m);
				for(int k = 0; k < edges.size(); k++)
				{
					Storage::real ecnt[3];
					if( edges[k].nbAdjElements(FACE,isFracture()) > 2 ) // more then 2 faces merge at this edge
					{
						//incorporate edge into fracture-fracture interface
						//create 2 quads
						ElementArray<Node> enodes1(m,2);
						ElementArray<Node> enodes2(m,2);
						enodes1[0] = nodesb[k];
						enodes1[1] = nodes[k];
						enodes2[0] = nodesb[(k+1)%N];
						enodes2[1] = nodes[(k+1)%N];
						
						f_edges.push_back(edgesb[k]);
						if( enodes1[0] != enodes1[1] ) {f_edges.push_back(m->CreateEdge(enodes1).first); new_edges++;}
						f_edges.push_back(edges[k]);
						if( enodes2[0] != enodes2[1] ) {f_edges.push_back(m->CreateEdge(enodes2).first); new_edges++;}
						
						assert(f_edges.size() > 2 );
						Face f1 = m->CreateFace(f_edges).first;
						new_faces++;
						f1->Real(fracture_edge_length) = edges[k]->Length();
						if( centroid_tag.isValid() )
						{
							edges[k]->Centroid(ecnt);
							f1->RealArray(centroid_tag)[0] = ecnt[0];
							f1->RealArray(centroid_tag)[1] = ecnt[1];
							f1->RealArray(centroid_tag)[2] = ecnt[2];
						}
						fracfaces.push_back(f1);

						f_edges.clear();

						ElementArray<Node> enodes3(m,2);
						ElementArray<Node> enodes4(m,2);

						enodes3[0] = nodesf[k];
						enodes3[1] = nodes[k];
						enodes4[0] = nodesf[(k+1)%N];
						enodes4[1] = nodes[(k+1)%N];

						if( enodes3[0] != enodes3[1] ) {f_edges.push_back(m->CreateEdge(enodes3).first); new_edges++;}
						f_edges.push_back(edgesf[k]);
						if( enodes4[0] != enodes4[1] ) {f_edges.push_back(m->CreateEdge(enodes4).first); new_edges++;}
						f_edges.push_back(edges[k]);
						

						assert(f_edges.size() > 2 );

						Face f2 = m->CreateFace(f_edges).first;
						new_faces++;
						if( centroid_tag.isValid() )
						{
							edges[k]->Centroid(ecnt);
							f2->RealArray(centroid_tag)[0] = ecnt[0];
							f2->RealArray(centroid_tag)[1] = ecnt[1];
							f2->RealArray(centroid_tag)[2] = ecnt[2];
						}
						f2->Real(fracture_edge_length) = edges[k]->Length();
						
						fracfaces.push_back(f2);

						f_edges.clear();
					}
					else if( nodes[k].GetMarker(isFracture()) || nodes[(k+1)%N].GetMarker(isFracture()) )
					{
						if( edgesb[k] != edgesf[k] )
						{
							//create quad
							ElementArray<Node> enodes(m,2);
							if( nodes[k].nbAdjElements(EDGE,isFracture()) && nodes[k].GetMarker(isFracture()) )
							{
								enodes[0] = nodesb[k];
								enodes[1] = nodes[k];
								f_edges.push_back(m->CreateEdge(enodes).first);
								new_edges++;

								enodes[0] = nodes[k];
								enodes[1] = nodesf[k];
								f_edges.push_back(m->CreateEdge(enodes).first);
								new_edges++;
							}
							else if( nodes[k].GetMarker(isFracture()) )
							{
								enodes[0] = nodesb[k];
								enodes[1] = nodesf[k];
								f_edges.push_back(m->CreateEdge(enodes).first);
								new_edges++;
							}
							f_edges.push_back(edgesf[k]);

							if( nodes[(k+1)%N].nbAdjElements(EDGE,isFracture()) && nodes[(k+1)%N].GetMarker(isFracture()) )
							{
								enodes[0] = nodesf[(k+1)%N];
								enodes[1] = nodes[(k+1)%N];
								f_edges.push_back(m->CreateEdge(enodes).first);
								new_edges++;

								enodes[0] = nodes[(k+1)%N];
								enodes[1] = nodesb[(k+1)%N];
								f_edges.push_back(m->CreateEdge(enodes).first);
								new_edges++;
							}
							else if( nodes[(k+1)%N].GetMarker(isFracture()) )
							{
								enodes[0] = nodesf[(k+1)%N];
								enodes[1] = nodesb[(k+1)%N];
								f_edges.push_back(m->CreateEdge(enodes).first);
								new_edges++;
							}
							f_edges.push_back(edgesb[k]);


							assert(f_edges.size() > 2 );

							Face f = m->CreateFace(f_edges).first;
							new_faces++;
							if( centroid_tag.isValid() )
							{
								edges[k]->Centroid(ecnt);
								f->RealArray(centroid_tag)[0] = ecnt[0];
								f->RealArray(centroid_tag)[1] = ecnt[1];
								f->RealArray(centroid_tag)[2] = ecnt[2];
							}
							f->Real(fracture_edge_length) = edges[k]->Length();
							fracfaces.push_back(f);

							f_edges.clear();
						}
					}
				}
			}
			else //it is boundary cell - get old nodes for nodesf
			{
				Face facesf = it->self();
				//add faces to indicate reconnection
				fracfaces.push_back(facesf);
				//now create fracture-fracture faces
				ElementArray<Edge> f_edges(m);
				for(int k = 0; k < edges.size(); k++)
				{
					if( edges[k].nbAdjElements(FACE,isFracture()) > 2 || nodes[k].GetMarker(isFracture()) || nodes[(k+1)%N].GetMarker(isFracture()) ) // more then 2 faces merge at this edge
					{
						//incorporate edge into fracture-fracture interface
						//create 2 quads
						ElementArray<Node> enodes1(m,2), enodes2(m,2);
						enodes1[0] = nodesb[k];
						enodes1[1] = nodes[k];
						enodes2[0] = nodesb[(k+1)%N];
						enodes2[1] = nodes[(k+1)%N];
						
						if( enodes1[0] != enodes1[1] ) {f_edges.push_back(m->CreateEdge(enodes1).first); new_edges++;}
						f_edges.push_back(edgesb.data()[k]);
						if( enodes2[0] != enodes2[1] ) {f_edges.push_back(m->CreateEdge(enodes2).first); new_edges++;}
						f_edges.push_back(edges.data()[k]);
						
						assert(f_edges.size() > 2);

						Face f = m->CreateFace(f_edges).first;
						new_faces++;
						fracfaces.push_back(f);
						f->Real(fracture_edge_length) = edges[k]->Length();
						f_edges.clear();
					}
				}
			}
			if( fill_fracture )
			{
				assert(fracfaces.size() >= 4); //tetrahedra or more
				/*
				std::cout << "fracfaces.size() = " << fracfaces.size() << std::endl;
				std::fstream fout("edgelist.txt",std::ios::out);
				for(int kk = 0; kk < fracfaces.size(); ++kk)
				{
					ElementArray<Edge> edges = fracfaces[kk].getEdges();
					for(int qq = 0; qq < edges.size(); ++qq)
					{
						std::cout << "(" << edges[qq]->getBeg()->Coords()[0] << "," << edges[qq]->getBeg()->Coords()[1] << "," << edges[qq]->getBeg()->Coords()[2] << ") <-> (" << edges[qq]->getEnd()->Coords()[0] << "," << edges[qq]->getEnd()->Coords()[1] << "," << edges[qq]->getEnd()->Coords()[2] << ")" << std::endl;
					}
				}
				fout.close();
				*/
				Cell fcell = m->CreateCell(fracfaces).first;
				new_cells++;
				fcell->Real(fracture_volume) = it->Area()*it->Real(fracture_aperture);
				//move permiability and volume factor to cell
				/*
				fcell->Real(volume_factor) = it->Real(volume_factor)*it->Area();
				if( K.isDefined(FACE) ) //what to do if it is not?
				{
					Storage::real_array Kcell = fcell->RealArray(K);
					Storage::real_array Kface = it->RealArray(K);
					if( Kcell.size() != Kface.size() ) //must be resizable
						Kcell.resize(Kface.size());
					memcpy(&Kcell[0],&Kface[0],sizeof(Storage::real)*Kcell.size());
				}
				*/
				for(size_t q = 0; q < transfer_real_tags.size(); ++q)
				{
					Storage::real_array c_arr = fcell->RealArray(transfer_real_tags[q]);
					Storage::real_array f_arr = it->RealArray(transfer_real_tags[q]);
					c_arr.replace(c_arr.begin(),c_arr.end(),f_arr.begin(),f_arr.end());
				}

				for(size_t q = 0; q < transfer_integer_tags.size(); ++q)
				{
					Storage::integer_array c_arr = fcell->IntegerArray(transfer_integer_tags[q]);
					Storage::integer_array f_arr = it->IntegerArray(transfer_integer_tags[q]);
					c_arr.replace(c_arr.begin(),c_arr.end(),f_arr.begin(),f_arr.end());
				}

				fcell->SetMarker(isFracture());
			}
			kid++;
		}
		std::cout << "kid " << kid << " new edges " << new_edges << " new faces " << new_faces << " new_cells " << new_cells << std::endl;
		std::cout << "Cells: " << m->NumberOfCells() << std::endl;
		std::cout << "Faces: " << m->NumberOfFaces() << std::endl;
		std::cout << "Edges: " << m->NumberOfEdges() << std::endl;
		std::cout << "Nodes: " << m->NumberOfNodes() << std::endl;
		//now reconnect matrix cells to new faces
		int reconnect = 0;
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) if( !it->GetMarker(isFracture()) && it->nbAdjElements(NODE,isFracture()) )
		{
			ElementArray<Face> change(m);
			ElementArray<Face> faces = it->getFaces();
			for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); jt++)
			{
				if( jt->nbAdjElements(NODE,isFracture()) )
					change.push_back(*jt);
			}
			if( !change.empty() )
			{
				reconnect++;
				Storage::reference_array newfaces = it->ReferenceArray(connface);
				assert(change.size() == newfaces.size());
				it->Disconnect(change.data(),(Storage::enumerator)change.size());
				it->Connect(newfaces.data(),(Storage::enumerator)newfaces.size());
			}
		}
		std::cout << "reconnect " << reconnect << std::endl;
		std::cout << "Cells: " << m->NumberOfCells() << std::endl;
		std::cout << "Faces: " << m->NumberOfFaces() << std::endl;
		std::cout << "Edges: " << m->NumberOfEdges() << std::endl;
		std::cout << "Nodes: " << m->NumberOfNodes() << std::endl;
	}
	else //number of dimensions is 2
	{
		//std::cout << __FILE__ << ":" << __LINE__ << "This part should be reimplemented the same way as for 3D" << std::endl;
		//throw NotImplemented;
		
		//mark all edges (in fact nodes) that have multiple fracture-faces (in fact edges)
		for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it)
		{
			if( it->nbAdjElements(FACE,isFracture()) > 2) it->SetMarker(isFracture());
		}

		//create new matrix-matrix faces (in fact edges), adjacent to fractures
		//for all non-fracture faces (in fact edges) that have any fracture edge (in fact node) create new face
		ElementArray<Node> enodes(m,1);
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( !it->GetMarker(isFracture()) && it->nbAdjElements(NODE,isFracture()) )
		{
			Storage::reference_array images = it->BackCell()->ReferenceArray(cell2node);
			ElementArray<Edge> fnodes(m,2);
			if( it->getBeg()->GetMarker(isFracture()) ) //search image of node between my elements
			{
				for(int k = 0; k < static_cast<int>(images.size()) && !fnodes[0].isValid(); k+=2)
				{
					if( images[k] == it->getBeg() ) // find my element between images
					{
						enodes[0] = images[k+1]->getAsNode();
						fnodes[0] = m->CreateEdge(enodes).first;
					}
				}
				assert(enodes[0] != InvalidElement());
				assert(fnodes[0] != InvalidElement());
			}
			else 
			{
				fnodes[0] = it->getBeg()->getEdges()[0];
				assert(fnodes[0] != InvalidElement());
			}
			

			if( it->getEnd()->GetMarker(isFracture()) )
			{
				for(int k = 0; k < static_cast<int>(images.size()) && !fnodes[1].isValid(); k+=2)
				{
					if( images[k] == it->getEnd() ) // find my element between images
					{
						enodes[0] = images[k+1]->getAsNode();
						fnodes[1] = m->CreateEdge(enodes).first;
					}
				}
				assert(enodes[0] != InvalidElement());
				assert(fnodes[1] != InvalidElement());
			}
			else 
			{
				fnodes[1] = it->getEnd()->getEdges()[0];
				assert(fnodes[1] != InvalidElement());
			}
			

			std::pair<Face,bool> f = m->CreateFace(fnodes);

			if( !f.second ) std::cout << __FILE__ << ":" << __LINE__ << " Face already exists!!! source " << *it << " new " << f.first.GetHandle()  << std::endl;
			else
			{
				for(int q = 0; q < (int)transfer_face_real_tags.size(); ++q)
					if( it->HaveData(transfer_face_real_tags[q]) )
					{
						Storage::real_array source = it->RealArray(transfer_face_real_tags[q]);
						Storage::real_array target = f.first->RealArray(transfer_face_real_tags[q]);
						if( target.size() != source.size() ) target.resize(source.size());
						for(int qq = 0; qq < source.size(); ++qq)
							target[qq] = source[qq];
					}
				for(int q = 0; q < (int)transfer_face_integer_tags.size(); ++q)
				 if( it->HaveData(transfer_face_integer_tags[q]) )
				 {
					 Storage::integer_array source = it->IntegerArray(transfer_face_integer_tags[q]);
					 Storage::integer_array target = f.first->IntegerArray(transfer_face_integer_tags[q]);
					 if( target.size() != source.size() ) target.resize(source.size());
					 for(int qq = 0; qq < source.size(); ++qq)
						 target[qq] = source[qq];
				 }
			}

			if(it->BackCell().isValid() ) it->BackCell()->ReferenceArray(connface).push_back(f.first);
			if(it->FrontCell().isValid() ) it->FrontCell()->ReferenceArray(connface).push_back(f.first);
		}

		if( !fracture_aperture.isDefined(CELL) ) // this is to extend volume factor
			m->CreateTag(fracture_aperture.GetTagName(),fracture_aperture.GetDataType(),CELL,CELL,fracture_aperture.GetSize());


		std::vector<Tag> transfer_real_tags;
		std::vector<Tag> transfer_integer_tags;
		for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); ++t)
			if( t->isDefined(FACE) && (t->GetDataType() == DATA_REAL || t->GetDataType() == DATA_INTEGER) )
			{
				if( t->GetTagName().substr(0,10) == "GEOM_UTIL_" &&
				    t->GetTagName() != "GEOM_UTIL_CENTROID" &&
				    t->GetTagName() != "GEOM_UTIL_BARYCENTER") continue;
				Tag q;
				if( !t->isDefined(CELL) )
					q = m->CreateTag(t->GetTagName(),t->GetDataType(),CELL,CELL,t->GetSize());
				else 
					q = *t;
				switch(t->GetDataType() )
				{
				case DATA_REAL:    transfer_real_tags.push_back(q); break;
				case DATA_INTEGER: transfer_integer_tags.push_back(q); break;
					default: printf("%s:%d unexpected\n",__FILE__,__LINE__); break;
				}
			}

		

		Tag num_breaks = m->CreateTag("NUM_BREAKS",DATA_INTEGER,CELL,NONE,1);
		Tag frac_num = m->CreateTag("FRACTURE_NUMBER",DATA_INTEGER,CELL,NONE,1);
		//now create fracture-fracture control volumes, add fracture-matrix faces to matrix cells
		int nfracs = 0;
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->GetMarker(isFracture()) && it->nbAdjElements(NODE,isFracture()) )
		{
			int nbreaks = 0; //number of breaks of faces by multiple joints
			Face cur = it->self();
			//I need all elements adjacent to neighbouring cells of current face
			//and ending at any element of current face
			ElementArray<Face> fracfaces(m);
			ElementArray<Node> nodes = it->getNodes(), enodes(m,1);
			ElementArray<Edge> edges(m,2), edgesb(m,2), edgesf(m,2);
			INMOST_DATA_ENUM_TYPE N = (INMOST_DATA_ENUM_TYPE)nodes.size();
			Storage::reference_array images;
			edges[0] = nodes[0]->getEdges()[0];
			edges[1] = nodes[1]->getEdges()[0];
			{
				images = it->BackCell()->ReferenceArray(cell2node);
				//among images of fracture nodes select those highlighted by marker
				for(int q = 0; q < nodes.size(); ++q) if( nodes[q].GetMarker(isFracture()) )
				{
					Node n;
					for(int k = 0; k < (int)images.size() && !n.isValid(); k+=2)
						if( images[k] == nodes[q] ) n = images[k+1]->getAsNode();
					assert(n.isValid());
					enodes[0] = n;
					edgesb[q] = m->CreateEdge(enodes).first;
				}
				else edgesb[q] = edges[q];
				//This is matrix-fracture face
				std::pair<Face,bool> facesb = m->CreateFace(edgesb);
				if(!facesb.second ) std::cout << __FILE__ << ":" << __LINE__ << " Face already exists!!! " << *it << " " << facesb.first->GetHandle() << std::endl;
				else
				{
					for(int q = 0; q < (int)transfer_face_real_tags.size(); ++q)
						if( it->HaveData(transfer_face_real_tags[q]) )
						{
							Storage::real_array source = it->RealArray(transfer_face_real_tags[q]);
							Storage::real_array target = facesb.first->RealArray(transfer_face_real_tags[q]);
							if( target.size() != source.size() ) target.resize(source.size());
							for(int qq = 0; qq < source.size(); ++qq)
								target[qq] = source[qq];
						}
					for(int q = 0; q < (int)transfer_face_integer_tags.size(); ++q)
						if( it->HaveData(transfer_face_integer_tags[q]) )
						{
							Storage::integer_array source = it->IntegerArray(transfer_face_integer_tags[q]);
							Storage::integer_array target = facesb.first->IntegerArray(transfer_face_integer_tags[q]);
							if( target.size() != source.size() ) target.resize(source.size());
							for(int qq = 0; qq < source.size(); ++qq)
								target[qq] = source[qq];
						}
				}
				//add faces to indicate reconnection
				it->BackCell()->ReferenceArray(connface).push_back(facesb.first);
				fracfaces.push_back(facesb.first);
			}

			if( it->FrontCell().isValid() ) //there is another cell
			{
				images = it->FrontCell()->ReferenceArray(cell2node);
				//among images of fracture nodes select those highlited by marker
				for(int q = 0; q < nodes.size(); ++q) if( nodes[q].GetMarker(isFracture()) )
				{
					Node n;
					for(int k = 0; k < (int)images.size() && !n.isValid(); k+=2)
						if( images[k] == nodes[q] ) n = images[k+1]->getAsNode();
					assert(n.isValid());
					enodes[0] = n;
					edgesf[q] = m->CreateEdge(enodes).first;
				}
				else edgesf[q] = edges[q];
				//This is matrix-fracture face
				std::pair<Face,bool> facesf = m->CreateFace(edgesf);
				if( !facesf.second ) std::cout << __FILE__ << ":" << __LINE__ << " Face already exists!!! " << *it << " " << facesf.first->GetHandle() << std::endl;
				else
				{
					for(int q = 0; q < (int)transfer_face_real_tags.size(); ++q)
						if( it->HaveData(transfer_face_real_tags[q]) )
						{
							Storage::real_array source = it->RealArray(transfer_face_real_tags[q]);
							Storage::real_array target = facesf.first->RealArray(transfer_face_real_tags[q]);
							if( target.size() != source.size() ) target.resize(source.size());
							for(int qq = 0; qq < source.size(); ++qq)
								target[qq] = source[qq];
						}
					for(int q = 0; q < (int)transfer_face_integer_tags.size(); ++q)
						if( it->HaveData(transfer_face_integer_tags[q]) )
						{
							Storage::integer_array source = it->IntegerArray(transfer_face_integer_tags[q]);
							Storage::integer_array target = facesf.first->IntegerArray(transfer_face_integer_tags[q]);
							if( target.size() != source.size() ) target.resize(source.size());
							for(int qq = 0; qq < source.size(); ++qq)
								target[qq] = source[qq];
						}
				}
				//add faces to indicate reconnection
				it->FrontCell()->ReferenceArray(connface).push_back(facesf.first);
				
				//now create fracture-fracture faces
				//HERE NODES MAY BE SIMILAR

				//we already have back matrix-fracture face in fracfaces
				//before adding front matrix-fracture face we should
				//add intermediate fracture-fracture faces to
				//preserve the order

				//add first fracture-fracture face
				if( edges[0].nbAdjElements(FACE,isFracture()) > 2 ) //in this case we should incorporate the fracture node into fracture-fracture faces
				{
					nbreaks++;
					ElementArray<Edge> fedges(m,2);
					fedges[0] = edgesb[0];
					fedges[1] = edges[0];
					//it may appear that elements match, no need to add new face
					if( fedges[0] != fedges[1] ) fracfaces.push_back(m->CreateFace(fedges).first);

					fedges[0] = edges[0];
					fedges[1] = edgesf[0];
					//it may appear that elements match, no need to add new face
					if( fedges[0] != fedges[1] ) fracfaces.push_back(m->CreateFace(fedges).first);
				}
				else //we can skip the node
				{
					ElementArray<Edge> fedges(m,2);
					fedges[0] = edgesb[0];
					fedges[1] = edgesf[0];
					//it may appear that elements match, no need to add new face
					if( fedges[0] != fedges[1] ) fracfaces.push_back(m->CreateFace(fedges).first);
				}

				//Now add front matrix-fracture face
				fracfaces.push_back(facesf.first);

				//add second fracture-fracture face
				if( edges[1].nbAdjElements(FACE,isFracture()) > 2 ) //in this case we should incorporate the fracture node into fracture-fracture faces
				{
					nbreaks++;
					ElementArray<Edge> fedges(m,2);
					fedges[0] = edgesf[1];
					fedges[1] = edges[1];
					//it may appear that elements match, no need to add new face
					if( fedges[0] != fedges[1] ) fracfaces.push_back(m->CreateFace(fedges).first);

					fedges[0] = edges[1];
					fedges[1] = edgesb[1];
					//it may appear that elements match, no need to add new face
					if( fedges[0] != fedges[1] ) fracfaces.push_back(m->CreateFace(fedges).first);
				}
				else //we can skip the node
				{
					ElementArray<Edge> fedges(m,2);
					fedges[0] = edgesf[1];
					fedges[1] = edgesb[1];
					//it may appear that elements match, no need to add new face
					if( fedges[0] != fedges[1] ) fracfaces.push_back(m->CreateFace(fedges).first);
				}
			}
			else //it is boundary cell
			{
				//we already have matrix-fracture face in fracfaces
				//should preserve the order
				//this is boundary face
				Face facesf = it->self();
				//now create first fracture-fracture face
				{
					ElementArray<Edge> fedges(m,2);
					fedges[0] = edgesb[0];
					fedges[1] = nodes[0].getEdges()[0];
					//it may appear that elements match, no need to add new face
					if( fedges[0] != fedges[1] ) fracfaces.push_back(m->CreateFace(fedges).first);
				}
				//now add boundary face
				fracfaces.push_back(facesf);
				//now create second fracture-fracture face
				{
					ElementArray<Edge> fedges(m,2);
					fedges[0] = edgesb[1];
					fedges[1] = nodes[1].getEdges()[0];
					//it may appear that elements match, no need to add new face
					if( fedges[0] != fedges[1] ) fracfaces.push_back(m->CreateFace(fedges).first);
				}
			}

			if( fill_fracture )
			{
				assert(fracfaces.size() >= 3); //triangle or more
				Cell fcell = m->CreateCell(fracfaces).first;
				fcell->Integer(num_breaks) = nbreaks;
				fcell->Integer(frac_num) = nfracs++;
				Storage::real fvol = it->Area()*it->Real(fracture_aperture);
				fcell->Real(fracture_volume) = fvol;
				//move permiability and volume factor to cell
				/*
				fcell->Real(volume_factor) = it->Real(volume_factor)*it->Area();
				if( K.isDefined(FACE) ) //what to do if it is not?
				{
					Storage::real_array Kcell = fcell->RealArray(K);
					Storage::real_array Kface = it->RealArray(K);
					if( Kcell.size() != Kface.size() ) //must be resizable
						Kcell.resize(Kface.size());
					memcpy(&Kcell[0],&Kface[0],sizeof(Storage::real)*Kcell.size());
				}
				*/
				for(size_t q = 0; q < transfer_real_tags.size(); ++q)
				{
					Storage::real_array c_arr = fcell->RealArray(transfer_real_tags[q]);
					Storage::real_array f_arr = it->RealArray(transfer_real_tags[q]);
					c_arr.replace(c_arr.begin(),c_arr.end(),f_arr.begin(),f_arr.end());
				}

				for(size_t q = 0; q < transfer_integer_tags.size(); ++q)
				{
					Storage::integer_array c_arr = fcell->IntegerArray(transfer_integer_tags[q]);
					Storage::integer_array f_arr = it->IntegerArray(transfer_integer_tags[q]);
					c_arr.replace(c_arr.begin(),c_arr.end(),f_arr.begin(),f_arr.end());
				}

				fcell->SetMarker(isFracture());
			}
		}
		//m->Save("fracture_test0.vtk");

		//now reconnect matrix cells to new faces
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) if( !it->GetMarker(isFracture()) && it->nbAdjElements(NODE,isFracture()) )
		{
			ElementArray<Face> change(m);
			ElementArray<Face> faces = it->getFaces();
			for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); jt++)
			{
				if( jt->nbAdjElements(NODE,isFracture()) )
					change.push_back(*jt);
			}
			if( !change.empty() )
			{
				Storage::reference_array newfaces = it->ReferenceArray(connface);
				assert(change.size() == newfaces.size());
				it->Disconnect(change.data(),(Storage::enumerator)change.size());
				it->Connect(newfaces.data(),(Storage::enumerator)newfaces.size());
			}
		}
	}
	m->DeleteTag(connface);
	m->DeleteTag(cell2node);
	for(Mesh::iteratorElement it = m->BeginElement(FACE|EDGE|NODE); it != m->EndElement(); ++it)
	{
		it->RemMarker(isFracture());
	}
	//delete remenants
	for(ElementType etypei = CELL; etypei > NODE; etypei = etypei >> 1)
	{
		ElementType etype = etypei >> 1;
		for(Mesh::iteratorElement it = m->BeginElement(etype); it != m->EndElement(); ++it)
		{
			if( !it->nbAdjElements(etypei) )
				it->Delete();
		}
		/* //output of element types
		std::cout << m->NumberOf(etypei) << " " << ElementTypeName(etypei) << ":" << std::endl;
		std::map<Element::GeometricType,int> gnum;
		for(Mesh::iteratorElement it = m->BeginElement(etypei); it != m->EndElement(); ++it)
		{
			gnum[it->GetGeometricType()]++;
		}
		for(std::map<Element::GeometricType,int>::iterator it = gnum.begin(); it != gnum.end(); ++it)
			std::cout << Element::GeometricTypeName(it->first) << ": " << it->second << std::endl;
			*/
	}
	//Mark all nodes edges and faces that entirely belong to fracture
	matrix_fracture = m->CreateMarker();
	multiple_fracture_joints = m->CreateMarker();
	m->self()->Integer(m->CreateTag("MATRIX_FRACTURE_MARKER",DATA_INTEGER,MESH,NONE,1)) = matrix_fracture;
	m->self()->Integer(m->CreateTag("FRACTURE_JOINTS_MARKER",DATA_INTEGER,MESH,NONE,1)) = multiple_fracture_joints;
	for(Mesh::iteratorElement it = m->BeginElement(FACE|EDGE|NODE); it != m->EndElement(); ++it)
	{
		Storage::enumerator nf = it->nbAdjElements(CELL,isFracture());
		if( it->nbAdjElements(CELL) == nf )
		{
			it->SetMarker(isFracture());
		}
		else if( nf )
		{
			it->SetMarker(matrix_fracture);
		}
	}
	//unmark all edges that have any node outside fracture
	for(Mesh::iteratorElement it = m->BeginElement(EDGE); it != m->EndElement(); ++it) if( it->GetMarker(isFracture()) )
	{
		if( it->nbAdjElements(NODE,isFracture(),true) != 0 )
			it->RemMarker(isFracture());
	}
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->GetMarker(isFracture()) )
	{
		if( it->nbAdjElements(EDGE,isFracture()) )
			it->SetMarker(multiple_fracture_joints);
	}
	//precalculate area for fractures
	std::cout << "Cells: " << m->NumberOfCells() << std::endl;
	std::cout << "Faces: " << m->NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m->NumberOfEdges() << std::endl;
	std::cout << "Nodes: " << m->NumberOfNodes() << std::endl;
	fracture_area = m->CreateTag("AREA_FRACTURE",DATA_REAL,FACE,FACE,1);
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->GetMarker(isFracture()) )
	{
		ElementArray<Edge> edges = it->getEdges();
		Storage::real half_perimeter = 0.0;
		if( m->GetDimensions() == 3 )
		{
			half_perimeter = it->Real(fracture_edge_length);
			//for(ElementArray<Edge>::iterator jt = edges.begin(); jt != edges.end(); ++jt)
			//	half_perimeter += jt->Length();
			//half_perimeter *= 0.5; // expect two long edges, all others are negligable
		} else half_perimeter = 1;
		Storage::real mean_aperture = 0.0;
		ElementArray<Cell> cells = it->getCells(isFracture());
		if( it->GetMarker(multiple_fracture_joints) )
		{
			for(ElementArray<Cell>::iterator jt = cells.begin(); jt != cells.end(); ++jt)
				mean_aperture += jt->Real(fracture_aperture)*jt->Real(fracture_aperture);
			mean_aperture = sqrt(mean_aperture);
		}
		else
		{
			for(ElementArray<Cell>::iterator jt = cells.begin(); jt != cells.end(); ++jt)
				mean_aperture += jt->Real(fracture_aperture);
		}
		mean_aperture *= cells.size() == 2 ? 0.5 : 1;
		it->Real(fracture_aperture) = mean_aperture;
		it->Real(fracture_area) = mean_aperture * half_perimeter;
	}
	m->ApplyModification();
	m->EndModification();
	std::cout << "Cells: " << m->NumberOfCells() << std::endl;
	std::cout << "Faces: " << m->NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m->NumberOfEdges() << std::endl;
	std::cout << "Nodes: " << m->NumberOfNodes() << std::endl;
	//m->Save("fracture_test.vtk");
}
