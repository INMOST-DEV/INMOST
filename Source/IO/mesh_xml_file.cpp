#include "inmost.h"


#if defined(USE_MESH)

namespace INMOST
{

  std::string ToLower(std::string input)
  {
	  std::string ret(input);
	  for(size_t k = 0; k < ret.size(); ++k)
		  ret[k] = tolower(ret[k]);
	  return ret;
  }

  void Mesh::LoadXML(std::string File)
  {
    std::fstream infile(File.c_str(),std::ios::in);
    XMLReader reader(File,infile);
    std::vector<Tag> tags;
    std::vector<HandleType> new_nodes, new_edges, new_faces, new_cells, new_sets;
    std::vector<ElementSet::ComparatorType> set_comparators;
    XMLReader::XMLTag PassTag;
    bool pass_tag = false;
    int nmeshes = 1;

    XMLReader::XMLTag TagParallelMesh = reader.OpenTag();
    if( TagParallelMesh.name != "ParallelMesh" )
    {
      reader.Report("Incorrect XML tag %s expected ParallelMesh",TagParallelMesh.name.c_str());
      throw BadFile;
    }

    for(int q = 0; q < TagParallelMesh.NumAttrib(); ++q)
    {
      const XMLReader::XMLAttrib & attr = TagParallelMesh.GetAttrib(q);
      if( ToLower(attr.name) == "number" ) nmeshes = atoi(attr.value.c_str());
      else reader.Report("Unused attribute for ParallelMesh %s='%s'",attr.name.c_str(),attr.value.c_str());
    }
    
    for(XMLReader::XMLTag TagMesh = reader.OpenTag(); !TagMesh.Finalize() && TagMesh.name == "Mesh"; reader.CloseTag(TagMesh), TagMesh = reader.OpenTag())
    {
      bool repair_orientation = false;
      for(int q = 0; q < TagMesh.NumAttrib(); ++q)
      {
        const XMLReader::XMLAttrib & attr = TagMesh.GetAttrib(q);
        if( ToLower(attr.name) == "repairorientation" ) repair_orientation = reader.ParseBool(attr.value);
        else reader.Report("Unused attribute for Tags %s='%s'",attr.name.c_str(),attr.value.c_str());
      }

      { //Nodes
        int nnodes = 0, ndims = 3;
        XMLReader::XMLTag TagNodes;
        if( pass_tag )
        {
          TagNodes = PassTag;
          pass_tag = false;
        }
        else
          TagNodes = reader.OpenTag();

        if( TagNodes.name != "Nodes" )
        {
          reader.Report("Incorrect XML tag %s expected Nodes",TagNodes.name.c_str());
          throw BadFile;
        }
		  bool matchnnodes = false;

        for(int q = 0; q < TagNodes.NumAttrib(); ++q)
        {
          const XMLReader::XMLAttrib & attr = TagNodes.GetAttrib(q);
          if( ToLower(attr.name) == "number" )
		  {
			  nnodes = atoi(attr.value.c_str());
			  matchnnodes = true;
		  }
          else if( ToLower(attr.name) == "dimensions" )
          {
            ndims = atoi(attr.value.c_str());
            if( GetDimensions() != ndims ) SetDimensions(ndims);
          }
          else reader.Report("Unused attribute for Tags %s='%s'",attr.name.c_str(),attr.value.c_str());
        }

        new_nodes.reserve(nnodes);

        {
          if( reader.ReadOpenContents() )
          {
            std::vector<double> Vector;
            int Repeat;
            dynarray<Storage::real,3> xyz;
            for(std::string val = reader.GetContentsWord(); !reader.isContentsEnded(); val = reader.GetContentsWord() )
            {
              reader.ParseReal(val,Vector,Repeat,nnodes);
              for(int l = 0; l < Repeat; ++l)
              {
                for(int q = 0; q < (int)Vector.size(); ++q)
                {
                  xyz.push_back(Vector[q]);
                  if( xyz.size() == ndims )
                  {
                    new_nodes.push_back(CreateNode(xyz.data())->GetHandle());
                    xyz.clear();
                  }
                }
              }
            }
            reader.ReadCloseContents();
          }
          else 
          {
            reader.Report("Cannot find contents of XML tag");
            throw BadFile;
          }

          reader.CloseTag(TagNodes);
        }

        if( matchnnodes && new_nodes.size() != nnodes )
        {
          reader.Report("Number of records for nodes %d do not match specified number of coordinates %d",new_nodes.size(),nnodes);
        }
      }

      { //Edges, Faces, Cells
        bool repeat = false;
          
        do
        {
          int nelems = 0, ntotconns = 0;
          bool matchelems = false;
          std::vector<HandleType> * elems;
          ElementType curtype;


          XMLReader::XMLTag TagElems;
          
          if( pass_tag )
          {
            TagElems = PassTag;
            pass_tag = false;
          }
          else
            TagElems = reader.OpenTag();

          if( !(TagElems.name == "Cells" || TagElems.name == "Faces" || TagElems.name == "Edges" ) )
          {
            reader.Report("Unexpected tag %s while waiting for either Cells or Faces or Edges",TagElems.name.c_str());
            throw BadFile;
          }
          if( TagElems.name != "Cells" ) repeat = true; else repeat = false;
          if( TagElems.name == "Cells" ) 
          {
            elems = &new_cells;
            curtype = CELL;
          }
          else if( TagElems.name == "Faces" ) 
          {
            elems = &new_faces;
            curtype = FACE;
          }
          else if( TagElems.name == "Edges" ) 
          {
            elems = &new_edges;
            curtype = EDGE;
          }
          

          for(int q = 0; q < TagElems.NumAttrib(); ++q)
          {
            const XMLReader::XMLAttrib & attr = TagElems.GetAttrib(q);
            if( ToLower(attr.name) == "number" ) 
            {
                nelems = atoi(attr.value.c_str());
                matchelems = true;
                //std::cout << "nelems: " << nelems << " text " << attr.value.c_str() << std::endl;
            }
            else reader.Report("Unused attribute for %ss %s='%s'",TagElems.name.c_str(),attr.name.c_str(),attr.value.c_str());
          } 
          
          elems->reserve(nelems);

            
          HandleType * links[3] =
          {
            new_nodes.empty() ? NULL : &new_nodes[0],
            new_edges.empty() ? NULL : &new_edges[0],
            new_faces.empty() ? NULL : &new_faces[0]
          };

          
          XMLReader::XMLTag TagConns;
          for(TagConns = reader.OpenTag(); !TagConns.Finalize() && TagConns.name == "Connections"; reader.CloseTag(TagConns), TagConns = reader.OpenTag())
          {
			  bool matchnconns = false;
			  int nexpectconns = 0;
            int nconns = 0;
            int offset = 0;
			int dims = 3; //to distinguish 3d cells from 2d cells when they are created with nodes
            ElementType subtype = NODE;
            for(int q = 0; q < TagConns.NumAttrib(); ++q)
            {
              const XMLReader::XMLAttrib & attr = TagConns.GetAttrib(q);
              if( ToLower(attr.name) == "number" )
			  {
				  nexpectconns = atoi(attr.value.c_str());
				  matchnconns = true;
			  }
			  else if( ToLower(attr.name) == "dimensions" ) dims = atoi(attr.value.c_str());
              else if( ToLower(attr.name) == "type" ) subtype = reader.atoes(attr.value.c_str());
              else if( ToLower(attr.name) == "offset" ) offset = atoi(attr.value.c_str());
              else reader.Report("Unused attribute for %ss %s='%s'",TagConns.name.c_str(),attr.name.c_str(),attr.value.c_str());
            }
            
            //ntotconns += nconns;
            if( subtype >= curtype )
            {
              reader.Report("%ss cannot be constructed from %ss",ElementTypeName(curtype),ElementTypeName(subtype));
              throw BadFile;
            }
            reader.ReadOpenContents();
            ElementArray<Element> subarr(this);
            for(std::string val = reader.GetContentsWord(); !reader.isContentsEnded(); val = reader.GetContentsWord())
            {
              int num = atoi(val.c_str()), elem;
              subarr.clear();
              subarr.reserve(num);
              for(int q = 0; q < num && !reader.isContentsEnded(); ++q)
              {
                val = reader.GetContentsWord();
                elem = atoi(val.c_str());
                //subarr.push_back(ComposeHandle(subtype,elem-offset));
                subarr.push_back(links[ElementNum(subtype)][elem-offset]);
              }
              switch(curtype)
              {
              case EDGE:
                elems->push_back(CreateEdge(subarr.Convert<Node>()).first.GetHandle());
                break;
              case FACE:
                if( subtype == NODE )
				{
					Face f = CreateFace(subarr.Convert<Node>()).first;
					if( repair_orientation ) f.FixEdgeOrder();
					elems->push_back(f.GetHandle());
				}
                else if( subtype == EDGE )
				{
					Face f = CreateFace(subarr.Convert<Edge>()).first;
					if( repair_orientation ) f.FixEdgeOrder();
					elems->push_back(f.GetHandle());
				}
                break;
              case CELL:
                if( subtype == NODE )
                {
					if( dims == 3 )
					{
						switch(subarr.size())
						{
						case 4: //Tetrahedron
							{
								const integer nodesnum[12] = {0,2,1,0,1,3,1,2,3,0,3,2};
								const integer sizes[4] = {3,3,3,3};
								elems->push_back(CreateCell(subarr.Convert<Node>(),nodesnum,sizes,4).first.GetHandle());
							}
							break;
						case 5: //Pyramid
							{
								const integer nodesnum[16] = {0,4,3,0,1,4,1,2,4,3,4,2,0,3,2,1};
								const integer sizes[5] = {3,3,3,3,4};
								elems->push_back(CreateCell(subarr.Convert<Node>(),nodesnum,sizes,5).first.GetHandle());
							}
							break;
						case 6: //Wedge or Prism
							{
								//const integer nodesnum[18] = {0,2,5,3,1,4,5,2,0,3,4,1,3,5,4,0,1,2};
								const integer nodesnum[18] = { 0, 3, 5, 2, 0, 1, 4, 3, 1, 4, 5, 2, 3, 4, 5, 0, 2, 1 };
								const integer sizes[5] = {4,4,4,3,3};
								elems->push_back(CreateCell(subarr.Convert<Node>(),nodesnum,sizes,5).first.GetHandle());
							}
							break;
						case 8: //Hexahedron
							{
								const integer nodesnum[24] = {0,4,7,3,1,2,6,5,0,1,5,4,3,7,6,2,0,3,2,1,4,5,6,7};
								const integer sizes[6] = {4,4,4,4,4,4};
								elems->push_back(CreateCell(subarr.Convert<Node>(),nodesnum,sizes,6).first.GetHandle());
							}
							break;
						default: 
							reader.Report("Sorry, no rule to convert %d nodes into a cell",subarr.size());
							throw BadFile;
							break;
						}
					}
					else if( dims == 2 )
					{
						ElementArray<Node> e_nodes(this,1);
						ElementArray<Edge> f_edges(this,2);
						ElementArray<Face> c_faces;
						for(unsigned int k = 0; k < subarr.size(); k++)
						{
							e_nodes.at(0) = subarr.at(k);
							f_edges.at(0) = CreateEdge(e_nodes).first->GetHandle();
							e_nodes.at(0) = subarr.at((k+1)%subarr.size());
							f_edges.at(1) = CreateEdge(e_nodes).first->GetHandle();
							c_faces.push_back(CreateFace(f_edges).first);
						}
						Cell c = CreateCell(c_faces).first;
						if( repair_orientation ) c.FixEdgeOrder();
						elems->push_back(c.GetHandle());
					}
					else reader.Report("Sorry, cannot understand number of dimensions %d for a cell",dims);
                }
                else if( subtype == EDGE )
                {
                  reader.Report("Sorry, no rule to convert %d edges into a cell",subarr.size());
                  throw BadFile;
                }
				else if( subtype == FACE )
				{
					Cell c = CreateCell(subarr.Convert<Face>()).first;
					if( repair_orientation ) c.FixEdgeOrder();
					elems->push_back(c.GetHandle());
				}
                break;
              }
				nconns++;
            }
            reader.ReadCloseContents();
			
			  ntotconns += nconns;
			  if( matchnconns && nconns != nexpectconns )
				  reader.Report("Number %d of elements encountered do not match to the specified number %d",nconns,nexpectconns);
          }

          if( matchelems && nelems != ntotconns) reader.Report("Number %d of elements encountered do not match to the specified number %d",ntotconns,nelems);

          reader.CloseTag(TagElems);

          if( !TagConns.Finalize() )
          {
            PassTag = TagConns;
            pass_tag = true;
          }

        } while(repeat);
      }

      { //Sets and Data
        bool repeat = false;
        do
        {
          XMLReader::XMLTag TagSetsData;
          if( pass_tag )
          {
            TagSetsData = PassTag;
            pass_tag = false;
          }
          else TagSetsData = reader.OpenTag();
          
		  if( TagSetsData.name == "Tags" )
		  {
			  repeat = true;
			  int ntags = 0;
			  bool matchntags = false;
			  for(int q = 0; q < TagSetsData.NumAttrib(); ++q)
			  {
				  const XMLReader::XMLAttrib & attr = TagSetsData.GetAttrib(q);
				  if( ToLower(attr.name) == "number" ) 
				  {
					  ntags = atoi(attr.value.c_str());
					  matchntags = true;
				  }
				  else reader.Report("Unused attribute for Tags %s='%s'",attr.name.c_str(),attr.value.c_str());
			  }

			  tags.reserve(ntags);


			  XMLReader::XMLTag TagTag;
			  for(TagTag = reader.OpenTag(); !TagTag.Finalize() && TagTag.name == "Tag"; reader.CloseTag(TagTag), TagTag = reader.OpenTag())
			  {
				  std::vector<std::string> parsed;
				  std::string tagname = "";
				  enumerator size = ENUMUNDEF;
				  DataType type = DATA_REAL;
				  ElementType sparse = NONE;
				  ElementType defined = NONE;
				  for(int q = 0; q < TagTag.NumAttrib(); ++q)
				  {
					  const XMLReader::XMLAttrib & attr = TagTag.GetAttrib(q);
					  if( ToLower(attr.name) == "name" ) tagname = attr.value;
					  else if( ToLower(attr.name) == "size" ) 
					  {
						  if( ToLower(attr.value) == "variable" )
							  size = ENUMUNDEF;
						  else size = atoi(attr.value.c_str());
					  }
					  else if( ToLower(attr.name) == "type" )
					  {
						  if( ToLower(attr.value) == "real" ) type = DATA_REAL;
						  else if( ToLower(attr.value) == "integer" ) type = DATA_INTEGER;
						  else if( ToLower(attr.value) == "reference" ) type = DATA_REFERENCE;
						  else if( ToLower(attr.value) == "remoteReference" ) type = DATA_REMOTE_REFERENCE;
						  else if( ToLower(attr.value) == "bulk" ) type = DATA_BULK;
#if defined(USE_AUTODIFF)
						  else if( ToLower(attr.value) == "variable" ) type = DATA_VARIABLE;
#endif
					  }
					  else if( ToLower(attr.name) == "sparse" )
					  { 
						  reader.ParseCommaSeparated(attr.value,parsed);
						  for(int q = 0; q < (int)parsed.size(); ++q) sparse |= reader.atoes(parsed[q].c_str());
					  }
					  else if( ToLower(attr.name) == "definition" )
					  {
						  reader.ParseCommaSeparated(attr.value,parsed);
						  for(int q = 0; q < (int)parsed.size(); ++q) defined |= reader.atoes(parsed[q].c_str());
					  }
					  else reader.Report("Unused attribute for Tags %s='%s'",attr.name.c_str(),attr.value.c_str());
				  }
				  if( tagname == "" )
					  reader.Report("Tag name was not specified");
				  else if( defined == NONE )
					  reader.Report("Domain of definition for the tag was not specified");
				  tags.push_back(CreateTag(tagname,type,defined,sparse,size));
			  }
			  reader.CloseTag(TagSetsData);

			  if( !TagTag.Finalize() )
			  {
				  PassTag = TagTag;
				  pass_tag = true;
			  }

			  if( matchntags && ntags != tags.size() ) reader.Report("Number %d of XML tags Tag red do not match to the specified number %d",tags.size(),ntags);
		  }
          else if( TagSetsData.name == "Sets" )
          {
            repeat = true;
            HandleType * links[4] =
            {
              new_nodes.empty() ? NULL : &new_nodes[0],
              new_edges.empty() ? NULL : &new_edges[0],
              new_faces.empty() ? NULL : &new_faces[0],
              new_cells.empty() ? NULL : &new_cells[0]
            };
            int nsets = 0, sets_offset = 0;
            int nsets_read = 0;
            bool matchsets = false;
            for(int q = 0; q < TagSetsData.NumAttrib(); ++q)
            {
              const XMLReader::XMLAttrib & attr = TagSetsData.GetAttrib(q);
              if( ToLower(attr.name) == "number" ) 
              {
                nsets = atoi(attr.value.c_str());
                matchsets = true;
              }
              else if( ToLower(attr.name) == "offset" ) sets_offset = atoi(attr.value.c_str());
              else reader.Report("Unused attribute for %ss %s='%s'",TagSetsData.name.c_str(),attr.name.c_str(),attr.value.c_str());
            }
            new_sets.reserve(nsets);
            
            XMLReader::XMLTag Set;
            for(Set = reader.OpenTag(); !Set.Finalize() && Set.name == "Set"; reader.CloseTag(Set), Set = reader.OpenTag() )
            {
				int expectsize = 0;
				bool matchsize = false;
              int size = 0, offset = 0;
              std::string name;
              HandleType parent = InvalidHandle(), child = InvalidHandle(), sibling = InvalidHandle();
              ElementSet::ComparatorType comparator = ElementSet::UNSORTED_COMPARATOR;
              nsets_read++;
              for(int q = 0; q < Set.NumAttrib(); ++q)
              {
                const XMLReader::XMLAttrib & attr = Set.GetAttrib(q);
                if( ToLower(attr.name) == "size" )
				{
					expectsize = atoi(attr.value.c_str());
					matchsize = true;
				}
                else if( ToLower(attr.name) == "offset" ) offset = atoi(attr.value.c_str());
                else if( ToLower(attr.name) == "name" ) name = attr.value;
                else if( ToLower(attr.name) == "parent" ) 
                {
                  if( ToLower(attr.value) != "unset" )
                    parent = ComposeHandle(ESET,atoi(attr.value.c_str()));
                }
                else if( ToLower(attr.name) == "sibling" ) 
                {
                  if( ToLower(attr.value) != "unset" )
                    sibling = ComposeHandle(ESET,atoi(attr.value.c_str()));
                }
                else if( ToLower(attr.name) == "child" ) 
                {
                  if( ToLower(attr.value) != "unset" )
                    child = ComposeHandle(ESET,atoi(attr.value.c_str()));
                }
                else if( ToLower(attr.name) == "comparator" )
                {
                  if( ToLower(attr.value) == "unsorted" ) comparator = ElementSet::UNSORTED_COMPARATOR;
                  else if( ToLower(attr.value) == "identificator" ) comparator = ElementSet::GLOBALID_COMPARATOR;
                  else if( ToLower(attr.value) == "centroid" ) comparator = ElementSet::CENTROID_COMPARATOR;
                  else if( ToLower(attr.value) == "hierarchy" ) comparator = ElementSet::HIERARCHY_COMPARATOR;
                  else if( ToLower(attr.value) == "handle" ) comparator = ElementSet::HANDLE_COMPARATOR;
                  else reader.Report("Unexpected comparator type %s for attribute Comparator, expected Unsorted,Identificator,Centroid,Hierarchy,Handle",attr.value.c_str());
                }
                else reader.Report("Unused attribute for %ss %s='%s'",Set.name.c_str(),attr.name.c_str(),attr.value.c_str());
              } 
              ElementSet s = CreateSet(name).first;
              new_sets.push_back(s.GetHandle());
              set_comparators.push_back(comparator);
              Element::adj_type & hc = HighConn(s.GetHandle());
              Element::adj_type & lc = LowConn(s.GetHandle());
              //here we believe that all the connections are consistent
              //hc.resize(ElementSet::high_conn_reserved);
              if( parent != InvalidHandle() ) ElementSet::hParent(hc) = parent;
              if( child != InvalidHandle() ) ElementSet::hChild(hc) = child;
              if( sibling != InvalidHandle() ) ElementSet::hSibling(hc) = sibling;
              //s.SortSet(comparator);
              reader.ReadOpenContents();
              for(std::string val = reader.GetContentsWord(); !reader.isContentsEnded(); val = reader.GetContentsWord())
              {
                std::vector<std::pair<ElementType,int> > Vector;
                int Repeat;
                reader.ParseReference(val,Vector,Repeat,0);
                if( !Vector.empty() ) for(int l = 0; l < Repeat; ++l) 
                {
                  for(int q = 0; q < (int)Vector.size(); ++q)
                  {
                    if( ElementNum(Vector[q].first) < 4 )
                      lc.push_back(links[ElementNum(Vector[q].first)][Vector[q].second-offset]);
                    else if( ElementNum(Vector[q].first) == 4 ) //Set
                      lc.push_back(ComposeHandle(Vector[q].first,Vector[q].second-offset));
                    else //Mesh
                      lc.push_back(GetHandle());
					size++;
                  }
                    //s.PutElement(links[ElementNum(Vector[q].first)][Vector[q].second-offset]);
                  //s.PutElements(Vector[0],(enumerator)Vector.size());
                }
              }
              reader.ReadCloseContents();
				
				if( matchsize && size != expectsize ) reader.Report("Number %d of elements encountered do not match to the specified set size %d",size, expectsize);
            }
            
            reader.CloseTag(TagSetsData);
            if( !Set.Finalize() )
            {
              PassTag = Set;
              pass_tag = true;
            }
            
            if( matchsets && nsets != nsets_read ) reader.Report("Number %d of XML tags Set read do not match to the number %d specified",nsets_read,nsets);

            //correct links between sets
            for(int q = 0; q < (int)new_sets.size(); ++q)
            {
              Element::adj_type & lc = HighConn(new_sets[q]);
              for(Element::adj_type::iterator jt = lc.begin(); jt != lc.end(); ++jt)
              {
                if( GetHandleElementType(*jt) == ESET )
                  *jt = new_sets[GetHandleID(*jt)];
              }
              Element::adj_type & hc = HighConn(new_sets[q]);
              if( ElementSet::hParent(hc) != InvalidHandle() ) ElementSet::hParent(hc) = new_sets[GetHandleID(ElementSet::hParent(hc))-sets_offset];
              if( ElementSet::hChild(hc) != InvalidHandle() ) ElementSet::hChild(hc) = new_sets[GetHandleID(ElementSet::hChild(hc))-sets_offset];
              if( ElementSet::hSibling(hc) != InvalidHandle() ) ElementSet::hSibling(hc) = new_sets[GetHandleID(ElementSet::hSibling(hc))-sets_offset];
              if( set_comparators[q] != ElementSet::UNSORTED_COMPARATOR ) ElementSet(this,new_sets[q]).SortSet(set_comparators[q]);
            }
          }
          else if( TagSetsData.name == "Data" )
          {
            repeat = false;
            HandleType * links[5] =
            {
              new_nodes.empty() ? NULL : &new_nodes[0],
              new_edges.empty() ? NULL : &new_edges[0],
              new_faces.empty() ? NULL : &new_faces[0],
              new_cells.empty() ? NULL : &new_cells[0],
              new_sets.empty() ? NULL : &new_sets[0]
            };
            int ndata = 0, ndatapass = 0;
            bool matchndata = false;

            for(int q = 0; q < TagSetsData.NumAttrib(); ++q)
            {
              const XMLReader::XMLAttrib & attr = TagSetsData.GetAttrib(q);
              if( ToLower(attr.name) == "number" ) 
              {
                ndata = atoi(attr.value.c_str());
                matchndata = true;
              }
              else reader.Report("Unused attribute for %ss %s='%s'",TagSetsData.name.c_str(),attr.name.c_str(),attr.value.c_str());
            }

            XMLReader::XMLTag TagDataSet;
            for(TagDataSet = reader.OpenTag(); !TagDataSet.Finalize() && TagDataSet.name == "DataSet"; reader.CloseTag(TagDataSet), TagDataSet = reader.OpenTag())
            {
              ++ndatapass;
              int sparse_read = 2, offset = 0;
              std::string tagname = "", setname = "", meshname = "";
              Mesh * remote_mesh = NULL;
              ElementType etype = NONE;
              for(int q = 0; q < TagDataSet.NumAttrib(); ++q)
              {
                const XMLReader::XMLAttrib & attr = TagDataSet.GetAttrib(q);
                if( ToLower(attr.name) == "settype" ) 
                {
                  if( ToLower(attr.value) == "cells" ) etype = CELL;
                  else if( ToLower(attr.value) == "faces" ) etype = FACE;
                  else if( ToLower(attr.value) == "edges" ) etype = EDGE;
                  else if( ToLower(attr.value) == "nodes" ) etype = NODE;
                  else if( ToLower(attr.value) == "sets" ) etype = ESET;
                  else if( ToLower(attr.value) == "mesh" ) etype = MESH;
                  else if( ToLower(attr.value) == "setdata" ) etype = NONE;
				  else reader.Report("Cannot understand attribute value for %s %s='%s'",TagDataSet.name.c_str(),attr.name.c_str(),attr.value.c_str());
                }
                else if( ToLower(attr.name) == "tagname" ) tagname = attr.value;
                else if( ToLower(attr.name) == "setname" ) setname = attr.value;
                else if( ToLower(attr.name) == "meshname" ) meshname = attr.value;
                else if( ToLower(attr.name) == "sparse" ) sparse_read = reader.ParseBool(attr.value);
                else if( ToLower(attr.name) == "offset" ) offset = atoi(attr.value.c_str());
                else reader.Report("Unused attribute for %s %s='%s'",TagDataSet.name.c_str(),attr.name.c_str(),attr.value.c_str());
              }
              
              if( tagname == "" )
              {
                reader.Report("DataSet had no attribute TagName");
                throw BadFile;
              }
              if( !HaveTag(tagname) )
              {
                reader.Report("Tag %s do not exist",tagname.c_str());
                throw BadFile;
              }
              if( meshname != "" )
              {
                remote_mesh = Mesh::GetMesh(meshname);
                if( !remote_mesh )
                  reader.Report("Remote mesh %s do not exist, you should create it first inside of your application",meshname.c_str());
              }

              Tag t = GetTag(tagname);
              HandleType * set_elems = NULL, *it = NULL;
              enumerator set_size = 0;

              if( setname != "" )
              {
                if( etype != NONE ) reader.Report("Warning: SetType should be SetData for data specified for set.");
                ElementSet s = GetSet(setname);
                if( !s.isValid() )
                {
                  reader.Report("Set %s is not defined on the mesh",setname.c_str());
                  throw BadFile;
                }
                set_elems = s.getHandles();
                set_size = s.Size(); 
                if( sparse_read == 2 ) sparse_read = 0;
              }
              else
              {
                if( etype == NONE )
                {
                  for(ElementType test = NODE; test <= MESH; test = NextElementType(test) )
                    if( t.isDefined(test) ) etype |= test;
                  if( !OneType(etype) )
                  {
                    reader.Report("You must explicitly specify type of elements for which"
                                  "the data is written since tag %s is defined on multiple types of elements",tagname.c_str());
                    throw BadFile;
                  }
                }
                if( sparse_read == 2 )
                {
                  if( t.isSparse(etype) ) sparse_read = 1; //Expect data to be listed in sparse manner for sparse tag
                  else sparse_read = 0;
                }
                switch(etype)
                {
                case NODE: set_elems = &new_nodes[0]; set_size = (enumerator)new_nodes.size(); break;
                case EDGE: set_elems = &new_edges[0]; set_size = (enumerator)new_edges.size(); break;
                case FACE: set_elems = &new_faces[0]; set_size = (enumerator)new_faces.size(); break;
                case CELL: set_elems = &new_cells[0]; set_size = (enumerator)new_cells.size(); break;
                case ESET: set_elems = &new_sets[0]; set_size = (enumerator)new_sets.size(); break;
                case MESH: set_elems = &handle; set_size = 1; break;
                }
              }
              it = set_elems;
              reader.ReadOpenContents();
              for(std::string val = reader.GetContentsWord(); !reader.isContentsEnded(); val = reader.GetContentsWord())
              {
                if( sparse_read )
                {
				  int offset;
				  if( etype == ESET && ((val[0]=='\'' && val[val.size()-1] == '\'') || (val[0] == '"' && val[1] == '"')) )
				  {
					  //helps check that the set was not found
					  offset = -1;
					  //strip qoutes
					  val = val.substr(1,val.size()-2);
					  //find set with provided name
					  for(enumerator k = 0; k < set_size; ++k)
						  if( ElementSet(this,set_elems[k])->GetName() == val )
						  {
							  offset = k;
							  break;
						  }
					  if( offset == -1 )
					  {
						  reader.Report("Set with name %s was not found",val.c_str());
						  throw BadFile;
					  }
				  }
				  else offset = atoi(val.c_str());
                  it = set_elems + offset;
                  val = reader.GetContentsWord();
                }
                switch(t.GetDataType())
                {
                case DATA_REAL:
                  {
                    Storage::real_array data = RealArray(*it,t);
                    std::vector<Storage::real> Vector; int Repeat;
                    reader.ParseReal(val,Vector,Repeat,set_size);
					if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                    else if( t.GetSize() != Vector.size()*Repeat && sparse_read )
                    {
                      reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                      throw BadFile;
                    }
                      
                    for(int l = 0; l < Repeat; ++l)
                    {
                      for(int q = 0; q < (int)Vector.size(); ++q)
                      {
                        data[(q + l*((int)Vector.size()))%data.size()] = Vector[q];
                        if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                        {
                          ++it;
                          if( ((int)(it-set_elems)) < (int)set_size ) data = RealArray(*it,t);
                        }
                      }
                    }
					if( !sparse_read && t.GetSize() == ENUMUNDEF )
					{
						++it;
						if( ((int)(it-set_elems)) < (int)set_size ) data = RealArray(*it,t);
					}
                  }
                  break;
                case DATA_INTEGER:
                  {
                    Storage::integer_array data = IntegerArray(*it,t);
                    std::vector<Storage::integer> Vector; int Repeat;
                    reader.ParseInteger(val,Vector,Repeat,set_size);
                    if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                    else if( t.GetSize() != Vector.size()*Repeat && sparse_read)
                    {
                      reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                      throw BadFile;
                    }
                    for(int l = 0; l < Repeat; ++l)
                    {
                      for(int q = 0; q < (int)Vector.size(); ++q)
                      {
                        data[(q + l*((int)Vector.size()))%data.size()] = Vector[q];
                        if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                        {
                          ++it;
                          if( ((int)(it-set_elems)) < (int)set_size ) data = IntegerArray(*it,t);
                        }
                      }
                    }
					if( !sparse_read && t.GetSize() == ENUMUNDEF )
					{
						++it;
						if( ((int)(it-set_elems)) < (int)set_size ) data = IntegerArray(*it,t);
					}
                  }
                  break;
                case DATA_BULK:
                  {
                    Storage::bulk_array data = BulkArray(*it,t);
                    std::string Vector; int Repeat;
                    reader.ParseBulk(val,Vector,Repeat,set_size);
                    if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                    else if( t.GetSize() != Vector.size()*Repeat && sparse_read)
                    {
                      reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                      throw BadFile;
                    }
                    for(int l = 0; l < Repeat; ++l)
                    {
                      for(int q = 0; q < (int)Vector.size(); ++q)
                      {
                        data[(q + l*((int)Vector.size()))%data.size()] = Vector[q];
                        if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                        {
                          ++it;
                          if( ((int)(it-set_elems)) < (int)set_size ) data = BulkArray(*it,t);
                        }
                      }
                    }
					  if( !sparse_read && t.GetSize() == ENUMUNDEF )
					  {
						  ++it;
						  if( ((int)(it-set_elems)) < (int)set_size ) data = BulkArray(*it,t);
					  }
                  }
                  break;
                case DATA_REFERENCE:
                  {
                    Storage::reference_array data = ReferenceArray(*it,t);
                    std::vector<std::pair<ElementType,int> > Vector; int Repeat;
                    reader.ParseReference(val,Vector,Repeat,set_size);
                    if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                    else if( t.GetSize() != Vector.size()*Repeat && sparse_read)
                    {
                      reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                      throw BadFile;
                    }
                    for(int l = 0; l < Repeat; ++l)
                    {
                      for(int q = 0; q < (int)Vector.size(); ++q)
                      {
                        data[(q + l*((int)Vector.size()))%data.size()] = Element(this,links[ElementNum(Vector[q].first)][Vector[q].second-offset]);
                        if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                        {
                          ++it;
                          if( ((int)(it-set_elems)) < (int)set_size ) data = ReferenceArray(*it,t);
                        }
                      }
                    }
					  if( !sparse_read && t.GetSize() == ENUMUNDEF )
					  {
						  ++it;
						  if( ((int)(it-set_elems)) < (int)set_size ) data = ReferenceArray(*it,t);
					  }
                  }
                  break;
                case DATA_REMOTE_REFERENCE:
                  {
                    Storage::remote_reference_array data = RemoteReferenceArray(*it,t);
                    if( remote_mesh != NULL )
                    {
                      std::vector<std::pair<ElementType,int> > Vector; int Repeat;
                      reader.ParseReference(val,Vector,Repeat,set_size);
                      if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                      else if( t.GetSize() != Vector.size()*Repeat && sparse_read )
                      {
                        reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                        throw BadFile;
                      }
                      for(int l = 0; l < Repeat; ++l)
                      {
                        for(int q = 0; q < (int)Vector.size(); ++q)
                        {
                          data.at((q + l*((int)Vector.size()))%data.size()).first  = remote_mesh;
                          data.at((q + l*((int)Vector.size()))%data.size()).second = ComposeHandle(Vector[q].first,Vector[q].second-offset);
                          if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                          {
                            ++it;
                            if( ((int)(it-set_elems)) < (int)set_size ) data = RemoteReferenceArray(*it,t);
                          }
                        }
                      }
						if( !sparse_read && t.GetSize() == ENUMUNDEF )
						{
							++it;
							if( ((int)(it-set_elems)) < (int)set_size ) data = RemoteReferenceArray(*it,t);
						}
                    }
                    else
                    {
                      std::vector<std::pair<std::string,std::pair<ElementType,int> > > Vector; int Repeat;
                      reader.ParseRemoteReference(val,Vector,Repeat,set_size);
                      if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                      else if( t.GetSize() != Vector.size()*Repeat && sparse_read )
                      {
                        reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                        throw BadFile;
                      }
                      for(int l = 0; l < Repeat; ++l)
                      {
                        for(int q = 0; q < (int)Vector.size(); ++q)
                        {
                          Mesh * m = GetMesh(Vector[q].first);
                          if( m == NULL ) reader.Report("Cannot find remote mesh %s, you should create it first inside of your application",Vector[q].first.c_str());
                          data.at((q + l*((int)Vector.size()))%data.size()).first = m;
                          data.at((q + l*((int)Vector.size()))%data.size()).second = ComposeHandle(Vector[q].second.first,Vector[q].second.second-offset);
                          if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                          {
                            ++it;
                            if( ((int)(it-set_elems)) < (int)set_size ) data = RemoteReferenceArray(*it,t);
                          }
                        }
                      }
						if( !sparse_read && t.GetSize() == ENUMUNDEF )
						{
							++it;
							if( ((int)(it-set_elems)) < (int)set_size ) data = RemoteReferenceArray(*it,t);
						}
                    }
                  }
                  break;
#if defined(USE_AUTODIFF)
                case DATA_VARIABLE:
                  {
                    Storage::var_array data = VariableArray(*it,t);
                    std::vector<Storage::var> Vector; int Repeat;
                    reader.ParseVariable(val,Vector,Repeat,set_size);
                    if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                    else if( t.GetSize() != Vector.size()*Repeat && sparse_read )
                    {
                      reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                      throw BadFile;
                    }
                    for(int l = 0; l < Repeat; ++l)
                    {
                      for(int q = 0; q < (int)Vector.size(); ++q)
                      {
                        data[(q + l*((int)Vector.size()))%data.size()] = Vector[q];
                        if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                        {
                          ++it;
                          if( ((int)(it-set_elems)) < (int)set_size ) data = VariableArray(*it,t);
                        }
                      }
                    }
					  if( !sparse_read && t.GetSize() == ENUMUNDEF )
					  {
						  ++it;
						  if( ((int)(it-set_elems)) < (int)set_size ) data = VariableArray(*it,t);
					  }
                  }
                  break;
#endif
                }
              } 
              reader.ReadCloseContents();
            }

            if( matchndata && ndata != ndatapass ) reader.Report("warning: Number %d of DataSet tags encountered do not match to the specified %d",ndatapass,ndata);

            reader.CloseTag(TagSetsData);

            if( !TagDataSet.Finalize() )
            {
              PassTag = TagDataSet;
              pass_tag = true;
            }
          }
		  else if( TagSetsData.name == "" )
			  repeat = false;
          else
          {
            reader.Report("Unexpected tag %s, expected Tags or Sets or Data",TagSetsData.name.c_str());
            throw BadFile;
          }
        } while(repeat);
      }
      RepairGeometricTags();

      if( repair_orientation )
      {
        int numfixed = 0;
        for(int q = 0; q < (int)new_faces.size(); ++q)
          if( Face(this,new_faces[q]).FixNormalOrientation() ) numfixed++;
        if( numfixed ) std::cout << "Fixed orientation of " << numfixed << " faces" << std::endl;
      }
    }
    reader.CloseTag(TagParallelMesh);
    infile.close();
  }

  void Mesh::SaveXML(std::string File)
  {
    std::fstream fout(File.c_str(),std::ios::out);
    fout << "<ParallelMesh>\n";
    fout << "\t<Mesh>\n";
    Tag idx = CreateTag("TEMPORARY_XML_ENUMERATOR",DATA_INTEGER,MESH|ESET|CELL|FACE|EDGE|NODE,NONE,1);
    Integer(GetHandle(),idx) = 0;
    fout << "\t\t<Nodes Number=\"" << NumberOfNodes() << "\" Dimensions=\"" << GetDimensions() << "\">\n";
    fout << "\t\t\t<![CDATA[\n";
    dynarray<Storage::real,3> xyz(GetDimensions());
    int cnt = 0;
    for(iteratorNode it = BeginNode(); it != EndNode(); ++it)
    {
      fout << "\t\t\t";
      it->Centroid(xyz.data());
      for(int q = 0; q < GetDimensions(); ++q)
        fout << xyz[q] << " ";
      fout << "\n";
      it->Integer(idx) = cnt++;
    }
    fout << "\t\t\t]]>\n";
    fout << "\t\t</Nodes>\n";
    fout << "\t\t<Edges>\n";
    fout << "\t\t\t<Connections Number=\"" << NumberOfEdges() << "\" Type=\"Nodes\">\n";
    fout << "\t\t\t<![CDATA[\n";
    cnt = 0;
    for(iteratorEdge it = BeginEdge(); it != EndEdge(); ++it)
    {
      ElementArray<Node> nodes = it->getNodes();
      fout << "\t\t\t" << nodes.size();
      for(ElementArray<Node>::iterator jt = nodes.begin(); jt != nodes.end(); ++jt)
        fout << " " << jt->Integer(idx);
      fout << "\n";
      it->Integer(idx) = cnt++;
    }
    fout << "\t\t\t]]>\n";
    fout << "\t\t\t</Connections>\n";
    fout << "\t\t</Edges>\n";
    fout << "\t\t<Faces>\n";
    fout << "\t\t\t<Connections Number=\"" << NumberOfFaces() << "\" Type=\"Edges\">\n";
    fout << "\t\t\t<![CDATA[\n";
    cnt = 0;
    for(iteratorFace it = BeginFace(); it != EndFace(); ++it)
    {
      ElementArray<Edge> edges = it->getEdges();
      fout << "\t\t\t" << edges.size();
      for(ElementArray<Edge>::iterator jt = edges.begin(); jt != edges.end(); ++jt)
        fout << " " << jt->Integer(idx);
      fout << "\n";
      it->Integer(idx) = cnt++;
    }
    fout << "\t\t\t]]>\n";
    fout << "\t\t\t</Connections>\n";
    fout << "\t\t</Faces>\n";
    fout << "\t\t<Cells>\n";
    fout << "\t\t\t<Connections Number=\"" << NumberOfCells() << "\" Type=\"Faces\">\n";
    fout << "\t\t\t<![CDATA[\n";
    cnt = 0;
    for(iteratorCell it = BeginCell(); it != EndCell(); ++it)
    {
      ElementArray<Face> faces = it->getFaces();
      fout << "\t\t\t" << faces.size();
      for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
        fout << " " << jt->Integer(idx);
      fout << "\n";
      it->Integer(idx) = cnt++;
    }
    fout << "\t\t\t]]>\n";
    fout << "\t\t\t</Connections>\n";
    fout << "\t\t</Cells>\n";
    fout << "\t\t<Sets>\n";
    cnt = 0;
    for(iteratorSet it = BeginSet(); it != EndSet(); ++it)
      it->Integer(idx) = cnt++;
    for(iteratorSet it = BeginSet(); it != EndSet(); ++it)
    {
      fout << "\t\t\t<Set Name=\"" << it->GetName() << "\"\n";
      if( it->HaveParent() )
      fout << "\t\t\t     Parent=\"" << it->GetParent()->Integer(idx) << "\"\n";
      if( it->HaveChild() )
      fout << "\t\t\t     Child=\"" << it->GetChild()->Integer(idx) << "\"\n";
      if( it->HaveSibling() )
      fout << "\t\t\t     Sibling=\"" << it->GetSibling()->Integer(idx) << "\"\n";
      if( it->GetComparator() != ElementSet::UNSORTED_COMPARATOR )
      {
        fout << "\t\t\t     Comparator=\"";
        switch(it->GetComparator())
        {
        case ElementSet::GLOBALID_COMPARATOR: fout << "Identificator"; break;
        case ElementSet::CENTROID_COMPARATOR: fout << "Centroid"; break;
        case ElementSet::HIERARCHY_COMPARATOR: fout << "Hierarchy"; break;
        case ElementSet::HANDLE_COMPARATOR: fout << "Handle"; break;
        }
        fout << "\"\n";
      }
      fout << "\t\t\t     Size=\"" << it->Size() << "\">\n";
      fout << "\t\t\t<![CDATA[\n";
      int endl_count = 0;
      fout << "\t\t\t";
      for(ElementSet::iterator jt = it->Begin(); jt != it->End(); ++jt)
      {
        switch(jt->GetElementType())
        {
        case NODE: fout << "Node:"; break;
        case EDGE: fout << "Edge:"; break;
        case FACE: fout << "Face:"; break;
        case CELL: fout << "Cell:"; break;
        case ESET: fout << "Set:"; break;
        case MESH: fout << "Mesh:"; break;
        }
        fout << jt->Integer(idx) << " ";
        endl_count++;
        if( endl_count % 8 == 0 )
          fout << "\n\t\t\t";
      }
      fout << "\n\t\t\t]]>\n";
      fout << "\t\t\t</Set>\n";
    }
    fout << "\t\t</Sets>\n";
	fout << "\t\t<Tags>\n";
    for(int k = 0; k < (int)tags.size(); ++k)
    {
	  if( tags[k] == idx ) continue;
      if( tags[k].GetTagName().substr(0,9) == "PROTECTED" ) continue;
      std::string names[6] = {"Nodes","Edges","Faces","Cells","Sets","Mesh"};
      std::string definition = "", sparse = "", type = "";
      switch(tags[k].GetDataType())
      {
      case DATA_REAL: type = "Real"; break;
      case DATA_INTEGER: type = "Integer"; break;
      case DATA_BULK: type = "Bulk"; break;
      case DATA_REFERENCE: type = "Reference"; break;
      case DATA_REMOTE_REFERENCE: type = "RemoteReference"; break;
#if defined(USE_AUTODIFF)
      case DATA_VARIABLE: type = "Variable"; break;
#endif
      }
      for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype) )
      {
        if( tags[k].isDefined(etype) ) definition += names[ElementNum(etype)] + ",";
        if( tags[k].isSparse(etype) ) sparse += names[ElementNum(etype)] + ",";
      }
      if( !definition.empty() ) definition.resize(definition.size()-1); //remove trailing comma
      if( !sparse.empty() ) sparse.resize(sparse.size()-1); //remove trailing comma
      if( sparse == "" ) sparse = "None";
      fout << "\t\t\t<Tag Name  =\"" << tags[k].GetTagName() << "\"\n";
      if( tags[k].GetSize() !=ENUMUNDEF )
      fout << "\t\t\t     Size  =\"" << tags[k].GetSize() << "\"\n";
      fout << "\t\t\t     Type  =\"" << type << "\"\n";
      fout << "\t\t\t     Sparse=\"" << sparse << "\"\n";
      fout << "\t\t\t     Definition=\"" << definition << "\"/>\n";
    }
	fout << "\t\t</Tags>\n";
    fout << "\t\t<Data>\n";
    for(iteratorTag t = BeginTag(); t != EndTag(); ++t)
    {
      if( *t == idx ) continue;
      if( t->GetTagName().substr(0,9) == "PROTECTED" ) continue;
      for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype)) if( t->isDefined(etype) )
      {
        fout << "\t\t\t<DataSet TagName=\"" << t->GetTagName() << "\"\n";
        fout << "\t\t\t         SetType=\"";
        switch(etype)
        {
        case NODE: fout << "Nodes"; break;
        case EDGE: fout << "Edges"; break;
        case FACE: fout << "Faces"; break;
        case CELL: fout << "Cells"; break;
        case ESET: fout << "Sets"; break;
        case MESH: fout << "Mesh"; break;
        }
        fout << "\">\n";
        fout << "\t\t\t<![CDATA[\n\t\t\t";
        int endl_count = 0;
        switch(t->GetDataType())
        {
        case DATA_REAL:
          for(iteratorStorage jt = Begin(etype); jt != End(); ++jt) if( jt->HaveData(*t) )
          {
            Storage::real_array data = jt->RealArray(*t);
            if( !data.empty() )
            {
              if( t->isSparse(etype) ) 
              {
                fout << jt->Integer(idx) << " ";
                endl_count++;
              }
              if( data.size() == 1 ) 
              {
                fout << data[0];
              }
              else
              {
                fout << "{";
                for(int q = 0; q < (int)data.size()-1; ++q)
                  fout << data[q] << ",";
                fout << data.back();
                fout << "}";
              }
              endl_count+=data.size();
              if( endl_count > 8 )
              {
                fout << "\n\t\t\t";
                endl_count = 0;
              }
              else fout << " ";
            }
			else if( !t->isSparse(etype) ) fout << " { } ";

          }
          break;
        case DATA_INTEGER:
          for(iteratorStorage jt = Begin(etype); jt != End(); ++jt) if( jt->HaveData(*t) )
          {
            Storage::integer_array data = jt->IntegerArray(*t);
            if( !data.empty() )
            {
              if( t->isSparse(etype) ) 
              {
                fout << jt->Integer(idx) << " ";
                endl_count++;
              }
              if( data.size() == 1 ) fout << data[0];
              else
              {
                fout << "{";
                for(int q = 0; q < (int)data.size()-1; ++q)
                  fout << data[q] << ",";
                fout << data.back();
                fout << "}";
              }
              endl_count+=data.size();
              if( endl_count > 8 )
              {
                fout << "\n\t\t\t";
                endl_count = 0;
              }
              else fout << " ";
              //fout << std::endl;
            }
			else if( !t->isSparse(etype) ) fout << " { } ";
          }
          break;
        case DATA_BULK:
          for(iteratorStorage jt = Begin(etype); jt != End(); ++jt) if( jt->HaveData(*t) )
          {
            Storage::bulk_array data = jt->BulkArray(*t);
            if( !data.empty() )
            {
              if( t->isSparse(etype) ) 
              {
                fout << jt->Integer(idx) << " ";
                endl_count++;
              }
              if( data.size() == 1 ) fout << CharToHex(data[0]);
              else
              {
                fout << "{";
                for(int q = 0; q < (int)data.size()-1; ++q)
                  fout << CharToHex(data[q]) << ",";
                fout << CharToHex(data.back());
                fout << "}";
              }
              endl_count+=data.size();
              if( endl_count > 8 )
              {
                fout << "\n\t\t\t";
                endl_count = 0;
              }
              else fout << " ";
              //fout << std::endl;
            }
			else if( !t->isSparse(etype) ) fout << " { } ";
          }
          break;
        case DATA_REFERENCE:
          for(iteratorStorage jt = Begin(etype); jt != End(); ++jt) if( jt->HaveData(*t) )
          {
            Storage::reference_array data = jt->ReferenceArray(*t);
            if( !data.empty() )
            {
              if( t->isSparse(etype) ) 
              {
                fout << jt->Integer(idx) << " ";
                endl_count++;
              }
              if( data.size() == 1 ) 
              {
                if( data[0].isValid() )
                  fout << ReferenceToString(data.at(0),data[0].Integer(idx));
                else
                  fout << "None:0";
              }
              else
              {
                fout << "{";
                for(int q = 0; q < (int)data.size()-1; ++q)
                {
                  if( data[q].isValid() )
                    fout << ReferenceToString(data.at(q),data[q].Integer(idx)) << ",";
                  else
                    fout << "None:0,";
                }
                if( data[data.size()-1].isValid() )
                  fout << ReferenceToString(data.at(data.size()-1),data[data.size()-1].Integer(idx));
                else
                  fout << "None:0";
                fout << "}";
              }
              endl_count+=data.size();
              if( endl_count > 8 )
              {
                fout << "\n\t\t\t";
                endl_count = 0;
              }
              else fout << " ";
              //fout << std::endl;
            }
			else if( !t->isSparse(etype) ) fout << " { } ";
          }
          break;
        case DATA_REMOTE_REFERENCE:
          for(iteratorStorage jt = Begin(etype); jt != End(); ++jt) if( jt->HaveData(*t) )
          {
            Storage::remote_reference_array data = jt->RemoteReferenceArray(*t);
            if( !data.empty() )
            {
              if( t->isSparse(etype) ) 
              {
                fout << jt->Integer(idx) << " ";
                endl_count++;
              }
              if( data.size() == 1 ) 
              {
                if( data[0].isValid() )
                  fout << data.at(0).first->GetMeshName() << ":" << ReferenceToString(data.at(0).second,data[0].Integer(idx));
                else
                  fout << ":None:0";
              }
              else
              {
                fout << "{";
                for(int q = 0; q < (int)data.size()-1; ++q)
                {
                  if( data[q].isValid() )
                    fout << data.at(q).first->GetMeshName() << ":" << ReferenceToString(data.at(q).second,data[q].Integer(idx)) << ",";
                  else
                    fout << ":None:0,";
                }
                if( data[data.size()-1].isValid() )
                  fout << data.at(data.size()-1).first->GetMeshName() << ":" << ReferenceToString(data.at(data.size()-1).second,data[data.size()-1].Integer(idx));
                else
                  fout << ":None:0";
                fout << "}";
              }
              endl_count+=data.size();
              if( endl_count > 8 )
              {
                fout << "\n\t\t\t";
                endl_count = 0;
              }
              else fout << " ";
              //fout << std::endl;
            }
			else if( !t->isSparse(etype) ) fout << " { } ";
          }
          break;
#if defined(USE_AUTODIFF)
        case DATA_VARIABLE:
          for(iteratorStorage jt = Begin(etype); jt != End(); ++jt) if( jt->HaveData(*t) )
          {
            Storage::var_array data = jt->VariableArray(*t);
            if( !data.empty() )
            {
              if( t->isSparse(etype) ) 
              {
                fout << jt->Integer(idx) << " ";
                endl_count++;
              }
              if( data.size() == 1 ) fout << VariableToString(data.at(0));
              else
              {
                fout << "{";
                for(int q = 0; q < (int)data.size()-1; ++q)
                  fout << VariableToString(data.at(q)) << ",";
                fout << VariableToString(data.back());
                fout << "}";
              }
              endl_count+=data.size();
              if( endl_count > 8 )
              {
                fout << "\n\t\t\t";
                endl_count = 0;
              }
              else fout << " ";
              //fout << std::endl;
            }
			else if( !t->isSparse(etype) ) fout << " { } ";
          }
          break;
#endif
        }
        fout << "\n\t\t\t]]>\n";
        fout << "\t\t\t</DataSet>\n";
      }
    }
    DeleteTag(idx);
    fout << "\t\t</Data>\n";
    fout << "\t</Mesh>\n";
    fout << "</ParallelMesh>\n";
    fout.close();
  }
}

#endif
