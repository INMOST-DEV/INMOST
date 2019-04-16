#include "amesh.h"
#include <iomanip>
#include <set>


#if defined(USE_PARALLEL_WRITE_TIME)
#define REPORT_MPI(x) {m->WriteTab(m->out_time) << "<MPI><![CDATA[" << #x << "]]></MPI>" << std::endl; x;}
#define REPORT_STR(x) {m->WriteTab(m->out_time) << "<TEXT><![CDATA[" << x << "]]></TEXT>" << std::endl;}
#define REPORT_VAL(str,x) {m->WriteTab(m->out_time) << "<VALUE name=\"" << str << "\"> <CONTENT><![CDATA[" << x << "]]></CONTENT> <CODE><![CDATA[" << #x << "]]></CODE></VALUE>" << std::endl;}
#define ENTER_FUNC() double all_time = Timer(); {m->WriteTab(m->out_time) << "<FUNCTION name=\"" << __FUNCTION__ << "\" id=\"func" << m->GetFuncID()++ << "\">" << std::endl; m->Enter();}
#define ENTER_BLOCK() { double btime = Timer(); m->WriteTab(m->out_time) << "<FUNCTION name=\"" << __FUNCTION__ << ":" << NameSlash(__FILE__) << ":" << __LINE__ << "\" id=\"func" << m->GetFuncID()++ << "\">" << std::endl; m->Enter();
#define EXIT_BLOCK() m->WriteTab(m->out_time) << "<TIME>" << Timer() - btime << "</TIME>" << std::endl; m->Exit(); m->WriteTab(m->out_time) << "</FUNCTION>" << std::endl;}
#define EXIT_FUNC() {m->WriteTab(m->out_time) << "<TIME>" << Timer() - all_time << "</TIME>" << std::endl; m->Exit(); m->WriteTab(m->out_time) << "</FUNCTION>" << std::endl;}
#define EXIT_FUNC_DIE() {m->WriteTab(m->out_time) << "<TIME>" << -1 << "</TIME>" << std::endl; m->Exit(); m->WriteTab(m->out_time) << "</FUNCTION>" << std::endl;}
#else
#define REPORT_MPI(x) x
#define REPORT_STR(x) {}
#define REPORT_VAL(str,x) {}
#define ENTER_BLOCK()
#define EXIT_BLOCK()
#define ENTER_FUNC() {}
#define EXIT_FUNC() {}
#define EXIT_FUNC_DIE()  {}
#endif

static std::string NameSlash(std::string input)
{
	for(unsigned l = input.size(); l > 0; --l)
		if( input[l-1] == '/' || input[l-1] == '\\' )
			return std::string(input.c_str() + l);
	return input;
}


//#include "../../Source/Misc/base64.h"
//using namespace std;

//from inmost
//std::string base64_encode(unsigned char const* buf, unsigned int bufLen);

	std::string base64_encode_(unsigned char const* buf, unsigned int bufLen)
	{
		static const std::string base64_chars =
		"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		"abcdefghijklmnopqrstuvwxyz"
		"0123456789+/";
		std::string ret;
		int i = 0;
		int j = 0;
		unsigned char char_array_3[3];
		unsigned char char_array_4[4];
		
		while (bufLen--)
		{
			char_array_3[i++] = *(buf++);
			if (i == 3)
			{
				char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
				char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
				char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
				char_array_4[3] = char_array_3[2] & 0x3f;
				for(i = 0; (i <4) ; i++)
					ret += base64_chars[char_array_4[i]];
				i = 0;
			}
		}
		
		if (i)
		{
			for(j = i; j < 3; j++)
				char_array_3[j] = '\0';
	  
			char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
			char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
			char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
			char_array_4[3] = char_array_3[2] & 0x3f;
	  
			for (j = 0; (j < i + 1); j++)
				ret += base64_chars[char_array_4[j]];
	  
			while((i++ < 3))
				ret += '=';
		}
		
		return ret;
	}

/// todo:
/// 1. coarsment
/// 2. strategy for faces/edges with faults
/// 3. geom model supportbn
/// 4. make class abstract virtual for user implementation of refinement and coarsment indicators
/// see in code todo:
namespace INMOST
{
	void CleanupSets(ElementSet set)
	{
		ElementSet::iterator it = set.Begin();
		while(it != set.End())
		{
			if( it->isValid() ) ++it;
			else it = set.Erase(it);
		}
		for(ElementSet child = set.GetChild(); child.isValid(); child = child.GetSibling())
			CleanupSets(child);
	}
	
	void ReduceMax(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		(void) size;
		element->Integer(tag) = std::max(element->Integer(tag),*((const INMOST_DATA_INTEGER_TYPE *)data));
	}
	
	void ReduceMin(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		(void) size;
		element->Integer(tag) = std::min(element->Integer(tag),*((const INMOST_DATA_INTEGER_TYPE *)data));
	}
	
	void ReduceUnion(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		(void) size;
		const INMOST_DATA_INTEGER_TYPE * idata = (const INMOST_DATA_INTEGER_TYPE *)data;
		Storage::integer_array odata = element->IntegerArray(tag);
		std::vector<int> tmp(size+odata.size());
		tmp.resize(std::set_union(odata.begin(),odata.end(),idata,idata+size,tmp.begin())-tmp.begin());
		odata.replace(odata.begin(),odata.end(),tmp.begin(),tmp.end());
	}

	/*
    void AdaptiveMesh::PrintSetLocal(std::string offset, ElementSet it, std::stringstream& ss)
    {
        std::stringstream ss1;
        ss1 << offset << rank << ": Set : " << std::setw(5) << it.GetName() << " ";
        for (int i = ss1.str().length(); i < 23; i++) ss1 << " ";
        ss << ss1.str();
        ss << std::setw(6);
        if (it.GetStatus() == Element::Shared) ss << "shared";
        else if (it.GetStatus() == Element::Ghost) ss << "ghost";
        else if (it.GetStatus() == Element::Owned) ss << "owned";
        else ss << "none";

        ss << " tag_owner (" << it.IntegerDF(m->OwnerTag()) << ")";

        //ss << "   level (" << level[it.self()] << ")  ";
        ss << " tag_processors (";
        std::stringstream ss2;
        Storage::integer_array arr = it.IntegerArrayDV(m->ProcessorsTag());
        for (int i = 0; i < arr.size(); i++)
            ss2 << arr[i] << " ";
        ss << std::setw(5) << ss2.str() <<")";
        

        ElementSet::iterator p = it.Begin();
        ss << "     | Refs: ";
        int first = 0;
        while(p != it.End())
        {
        //    if (first++ == 0) ss << endl << offset << "   ";
            std::string type = "unknw";
            if (p->GetElementType() == CELL) type = "cell";
            if (p->GetElementType() == FACE) type = "face";
            if (p->GetElementType() == EDGE) type = "edge";
            if (p->GetElementType() == NODE) type = "node";
            ss << type << "-" << std::setw(2) << p->GlobalID() << " ";
            p++;
        }
        ss << std::endl;

        for(ElementSet child = it.GetChild(); child.isValid(); child = child.GetSibling())
        {
            PrintSetLocal(offset + "   ",child,ss);
        }
    }
	*/

	/*
    void AdaptiveMesh::PrintSet()
    {
        std::stringstream ss;
        for(Mesh::iteratorSet it = m->BeginSet(); it != m->EndSet(); ++it) 
        {
            //if (it->HaveParent()) continue;
            PrintSetLocal("",ElementSet(m,*it),ss);
        }
        std::cout << ss.str() << std::endl;
    }
	*/

	/*
    void PrintRefs(std::ostream& os, Storage::reference_array refs)
    {
        for(Storage::reference_array::size_type i = 0; i < refs.size(); ++i)
        {
            std::string type = "unknw";
            if (refs[i].GetElementType() == CELL) type = "cell";
            if (refs[i].GetElementType() == FACE) type = "face";
            if (refs[i].GetElementType() == EDGE) type = "edge";
            if (refs[i].GetElementType() == NODE) type = "node";
            os << "(" << type << "," << refs[i]->GlobalID() << ") ";
			Storage::real xyz[3] = {0,0,0};
            refs[i]->Centroid(xyz);
            os << "(" << xyz[0] << "," << xyz[1] << "," << xyz[2] <<")" <<  std::endl;
        }
    }
	*/

	/*
    void AdaptiveMesh::PrintMesh(std::ostream& os, int cell, int face, int edge, int node)
    {
        if (cell + face + edge + node == 0) return;
        std::stringstream ss;
        ss << "================= " << rank << " =====================" << std::endl;
        if (cell)
        {
            ss << "Cells: " << m->NumberOfCells() <<  std::endl;
            for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) 
            {
                ss << rank << ": " << it->GlobalID() << " - " << it->LocalID() << " - ";            
                if (it->GetStatus() == Element::Shared) ss << "shared";
                else if (it->GetStatus() == Element::Ghost) ss << "ghost";
                else ss << "none";
                ss << std::endl;
            }
        }

        if (face)
        {
            ss << "Faces: " << m->NumberOfFaces() <<  std::endl;
            for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
            {
                ss << rank << ": " << std::setw(2) << it->LocalID() << " " << std::setw(2) << it->GlobalID() << " - " ;            
                ss << std::setw(6);
                if (it->GetStatus() == Element::Shared) ss << "shared";
                else if (it->GetStatus() == Element::Ghost) ss << "ghost";
                else ss << "none";


                double xyz[3];
                it->Centroid(xyz);
                ss << "   (" << std::setw(5) << xyz[0] << " " << std::setw(5) << xyz[1] << " " << std::setw(5) << xyz[2] << ")";

                ss << "  " << m->GetMarker(*it,m->NewMarker());

                ss << " nc(" << it->getNodes().size() << ": "; 
                ElementArray<Node> nodes = it->getNodes();
                for (ElementArray<Node>::iterator node = nodes.begin(); node != nodes.end(); node++)
                    ss << std::setw(2) << node->GlobalID() << " ";
                ss << ")";
                ss << std::endl;
            }
        }

        if (edge)
        {
            ss << "Edges: " << m->NumberOfEdges() <<  std::endl;
            for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it) 
            {
                ss << rank << ": " << it->GlobalID() << " - " ;            
                if (it->GetStatus() == Element::Shared) ss << "shared";
                else if (it->GetStatus() == Element::Ghost) ss << "ghost";
                else ss << "none";
                ss << std::endl;
            }
        }

        if (node)
        {
            ss << "Nodes:" << std::endl;
            for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it) 
            {
                ss << rank << ": " << std::setw(2) << it->GlobalID() << " - " ;            
                ss << std::setw(6);
                if (it->GetStatus() == Element::Shared) ss << "shared";
                else if (it->GetStatus() == Element::Ghost) ss << "ghost";
                else ss << "none";

                {
                    ss << "  " << m->GetMarker(*it,m->NewMarker());

                    ss << "(" << 
                        std::setw(3) << it->RealArray(m->CoordsTag())[0] << " " << 
                        std::setw(3) << it->RealArray(m->CoordsTag())[1] << " " << 
                        std::setw(3) << it->RealArray(m->CoordsTag())[2] << ")";

                }
                ss << std::endl;
            }
        }

        ss << "=========================================" << std::endl;
        os << ss.str() << std::endl;
    }
	*/

	/*
    void AdaptiveMesh::UpdateStatus()
    {

        for(ElementType mask = CELL; mask >= NODE; mask = PrevElementType(mask))
        {
            for(Mesh::iteratorElement it = m->BeginElement(mask); it != m->EndElement(); it++)
            {
                int stat = 0;
                if (it->GetStatus() == Element::Shared) stat = 1;
                else if (it->GetStatus() == Element::Ghost)  stat = 2;

                tag_status[it->self()] = stat;
            }
        }
    }
	*/
    
	
	void AdaptiveMesh::PrintSet(std::ostream & fout, ElementSet set)
    {
		fout << "set " << set.GetName();
		fout << " size " << set.Size();
		if(set.HaveParent()) fout << " parent " << set.GetParent().GetName();
		fout << " elements ";
		for(ElementSet::iterator it = set.Begin(); it != set.End(); ++it)
			fout << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << ":" << it->GlobalID() << " ";
		fout << " children ";
		for(ElementSet child = set.GetChild(); child.isValid(); child = child.GetSibling())
            fout << child.GetName() << " ";
		fout << std::endl;
    }
	

	/*
    void AdaptiveMesh::SynchronizeSet(ElementSet set)
    {
    #ifdef USE_MPI
        int size = m->GetProcessorsNumber();
        int rank = m->GetProcessorRank();
        for (int i = 0; i < size; i++)
        {
            set.IntegerArray(m->SendtoTag()).push_back(i);
            m->ExchangeMarked();
        }
    #endif
    }
	*/

	/*
    void AdaptiveMesh::Test()
    {
        std::cout << rank << ": ================" << std::endl;
        PrintSet(root,"");
    }
	*/
	
	void AdaptiveMesh::ClearData()
	{
		set_id        = m->DeleteTag(set_id);
		level         = m->DeleteTag(level);
		hanging_nodes = m->DeleteTag(hanging_nodes);
		parent_set    = m->DeleteTag(parent_set);
		root.DeleteSetTree();
	}
	
	void AdaptiveMesh::PrepareSet()
	{
		//retrive set for coarsening, initialize set if is not present
		if( !root.isValid() )
		{
			root = m->GetSet("AM_ROOT_SET");
			if( root == InvalidElement() )
			{
				root = m->CreateSetUnique("AM_ROOT_SET").first;
				//root.SetExchange(ElementSet::SYNC_ELEMENTS_SHARED);
				level[root] = 0;
				for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				{
					root.PutElement(it->self());
					parent_set[it->self()] = root.GetHandle();
				}
				m->Enumerate(CELL,set_id);
			}
        }
		if( !m->HaveGlobalID(CELL) ) m->AssignGlobalID(CELL); //for unique set names
		m->ResolveSets();
	}
	
	AdaptiveMesh::AdaptiveMesh(Mesh & _m) : m(&_m)
	{
		model = NULL;
		//create a tag that stores maximal refinement level of each element
		level = m->CreateTag("REFINEMENT_LEVEL",DATA_INTEGER,CELL|FACE|EDGE|NODE|ESET,NONE,1);
		//tag_status = m->CreateTag("TAG_STATUS",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
		set_id = m->CreateTag("SET_ID",DATA_INTEGER,CELL,NONE,1);
		//tag_an = m->CreateTag("TAG_AN",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
		//ref_tag = m->CreateTag("REF",DATA_REFERENCE,CELL|FACE|EDGE|NODE,NONE);
		//create a tag that stores links to all the hanging nodes of the cell
		hanging_nodes = m->CreateTag("HANGING_NODES",DATA_REFERENCE,CELL|FACE,NONE);
		//create a tag that stores links to sets
		parent_set = m->CreateTag("PARENT_SET",DATA_REFERENCE,CELL,NONE,1);
	    size = m->GetProcessorsNumber();
    	rank = m->GetProcessorRank();
	}
	
	AdaptiveMesh::~AdaptiveMesh()
	{
		//do not delete tags, user may want to repetitively use this class
		//as extension of class mesh in limited code span
	}
	
	void AdaptiveMesh::CheckParentSet()
	{
		ENTER_FUNC();
		int err = 0;
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		{
			if( parent_set[*it] == InvalidHandle() )
			{
				REPORT_STR(m->GetProcessorRank() << " parent set not valid on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << " " << level[*it]);
				err++;
			}
			else if( GetHandleElementType(parent_set[*it]) != ESET )
			{
				REPORT_STR(m->GetProcessorRank() << " parent set is something else " << ElementTypeName(GetHandleElementType(parent_set[*it])) << ":" << GetHandleID(parent_set[*it]) << " on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << " " << level[*it]);
				err++;
			}
			else if( parent_set[*it] != root.GetHandle() )
			{
				ElementSet set(m,parent_set[*it]);
				if( !set.HaveParent() )
				{
					REPORT_STR(m->GetProcessorRank() << 
					" parent set " << ElementTypeName(GetHandleElementType(parent_set[*it])) << ":" << GetHandleID(parent_set[*it]) << 
					" name " << set.GetName() << " owner " << set.Integer(m->OwnerTag()) << " status " << Element::StatusName(set.GetStatus()) << 
					" does not have parent " <<
					" on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << " " << level[*it]);
					err++;
				}
				else
				{
					HandleType parent = set.GetParent().GetHandle();
					if( GetHandleElementType(parent) != ESET )
					{
						REPORT_STR(m->GetProcessorRank() << 
						" parent set " << ElementTypeName(GetHandleElementType(parent_set[*it])) << ":" << GetHandleID(parent_set[*it]) << 
						" name " << set.GetName() << " owner " << set.Integer(m->OwnerTag()) << " status " << Element::StatusName(set.GetStatus()) << 
						" has parent " << ElementTypeName(GetHandleElementType(parent)) << ":" << GetHandleID(parent) <<
						" on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << " " << level[*it]);
						err++;
					}
				}
			}
		}
		err = m->Integrate(err);
		EXIT_FUNC();
		if( err ) 
		{
			REPORT_STR(rank << " error in " << __FUNCTION__);
			std::cout << rank << " error in " << __FUNCTION__ << std::endl;
			exit(-1);
		}
		
	}
	
	bool AdaptiveMesh::Refine(TagInteger & indicator)
	{
		static int fi = 0;
        ENTER_FUNC();
		static int call_counter = 0;
		int ret = 0; //return number of refined cells
		//initialize tree structure
		//m->CheckCentroids(__FILE__,__LINE__);
		ENTER_BLOCK();
		PrepareSet();
		EXIT_BLOCK();

		//m->Save("before_refine"+std::to_string(fi)+".pvtk");
		//std::cout << "Save before_refine"+std::to_string(fi)+".pvtk" << std::endl;
		
		//ENTER_BLOCK();
		//for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it)
		//	if( it->getNodes().size() != 2 ) {REPORT_STR("edge " << it->LocalID() << " has " << it->getNodes().size() << " nodes ");}
		//EXIT_BLOCK();
		
		ENTER_BLOCK();
		m->ResolveSets();
		//m->CheckCentroids(__FILE__,__LINE__);
		m->ExchangeData(parent_set,CELL,0);
		//m->CheckCentroids(__FILE__,__LINE__);
		m->ResolveSets();
		//m->CheckCentroids(__FILE__,__LINE__);
		//m->ExchangeData(parent_set,CELL,0);
		//m->CheckCentroids(__FILE__,__LINE__);
		m->ExchangeData(hanging_nodes,CELL | FACE,0);
		//m->CheckCentroids(__FILE__,__LINE__);
		CheckParentSet();
		EXIT_BLOCK();
		
		//ENTER_BLOCK();
		//for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it)
		//	if( it->getNodes().size() != 2 ) {REPORT_STR("edge " << it->LocalID() << " has " << it->getNodes().size() << " nodes ");}
		//EXIT_BLOCK();

		//m->Save("before_refine_parent"+std::to_string(fi)+".pvtk");
		//std::cout << "Save before_refine_parent"+std::to_string(fi)+".pvtk" << std::endl;
		

		
		int schedule_counter = 1; //indicates order in which refinement will be scheduled
		int scheduled = 1; //indicates that at least one element was scheduled on current sweep
		ENTER_BLOCK();
		//0. Extend indicator for edges and faces
		indicator = m->CreateTag(indicator.GetTagName(),DATA_INTEGER,FACE|EDGE,NONE,1);
		while(scheduled)
		{
			REPORT_VAL("scheduled",scheduled);
			//1.Communicate indicator - it may be not synced
			m->ExchangeData(indicator,CELL,0);
			//2.Propogate indicator down to the faces,edges
			//  select schedule for them
			ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( indicator[c] == schedule_counter )
				{
					ElementArray<Element> adj = c.getAdjElements(FACE|EDGE);
					for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
					{
						if( level[adj[kt]] == level[c] ) //do not schedule finer or coarser elements
							indicator[adj[kt]] = schedule_counter; //refine together with current cell
					}
				}
			}
			EXIT_BLOCK();
			//3.Communicate indicator on faces and edges
			m->ExchangeData(indicator,FACE|EDGE,0);
			//4.Check for each cell if there is
			//  any hanging node with adjacent in a need to refine,
			//  schedule for refinement earlier.
			scheduled = 0;
			ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				//already scheduled cells may be required to be refined first
				//if( indicator[c] == 0 ) //some optimization
				{
					bool scheduled_c = false;
					//any finer level edge is scheduled to be refined first
					ElementArray<Edge> edges = c->getEdges();
					for(ElementArray<Edge>::size_type kt = 0; kt < edges.size() && !scheduled_c; ++kt)
					{
						//if a finer edge is scheduled
						//then this cell should be refined first
						if( indicator[edges[kt]] != 0 &&
							level[edges[kt]] > level[c] &&
							indicator[edges[kt]] >= indicator[c] )
						{
							indicator[c] = schedule_counter+1;
							scheduled++;
							scheduled_c = true;
						}
					}
				}
			}
			EXIT_BLOCK();
			//5.Go back to 1 until no new elements scheduled
			scheduled = m->Integrate(scheduled);
			if( scheduled ) schedule_counter++;
		}
		m->ExchangeData(indicator,CELL | FACE | EDGE,0);
		EXIT_BLOCK();
		//m->Save("indicator.pmf");
		
		//ENTER_BLOCK();
		//for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it)
		//	if( it->getNodes().size() != 2 ) {REPORT_STR("edge " << it->LocalID() << " has " << it->getNodes().size() << " nodes ");}
		//EXIT_BLOCK();
		
		//if( !Element::CheckConnectivity(m) ) std::cout << __FILE__ << ":" << __LINE__ << " broken connectivity" << std::endl;
		//m->CheckCentroids(__FILE__,__LINE__);
		//6.Refine
		ENTER_BLOCK();
		m->BeginModification();
		while(schedule_counter)
		{
			Storage::real xyz[3] = {0,0,0};
			//7.split all edges of the current schedule
			ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->EdgeLastLocalID(); ++it) if( m->isValidEdge(it) )
			{
				Edge e = m->EdgeByLocalID(it);
				if( !e.Hidden() && indicator[e] == schedule_counter )
				{
					//ENTER_BLOCK();
					//for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it)
					//	if( it->getNodes().size() != 2 ) {REPORT_STR("edge " << it->LocalID() << " has " << it->getNodes().size() << " nodes ");}
					//EXIT_BLOCK();
					//REPORT_STR("split edge " << e.LocalID() << " nodes " << e.getBeg().LocalID() << "," << e.getEnd().LocalID() << " level " << level[e] << " lc size " << m->LowConn(e.GetHandle()).size() );
					//ElementArray<Node> nodes = e.getNodes();
					//for(int q = 0; q < nodes.size(); ++q) REPORT_STR("node " << nodes[q].GetHandle() << " " << nodes[q].LocalID() << (nodes[q].Hidden()?" hidden " : " good ") );
					//remember adjacent faces that should get information about new hanging node
					ElementArray<Face> edge_faces = e.getFaces();
					//location on the center of the edge
					for(Storage::integer d = 0; d < m->GetDimensions(); ++d)
						xyz[d] = (e.getBeg().Coords()[d]+e.getEnd().Coords()[d])*0.5;
					//todo: request transformation of node location according to geometrical model
					//create middle node
					Node n = m->CreateNode(xyz);
					//set increased level for new node
					level[n] = level[e.getBeg()] = level[e.getEnd()] = level[e]+1;
					//for each face provide link to a new hanging node
					for(ElementArray<Face>::size_type kt = 0; kt < edge_faces.size(); ++kt)
						hanging_nodes[edge_faces[kt]].push_back(n);
					//split the edge by the middle node
					ElementArray<Edge> new_edges = Edge::SplitEdge(e,ElementArray<Node>(m,1,n.GetHandle()),0);
					//set increased level for new edges
					level[new_edges[0]] = level[new_edges[1]] = level[e]+1;
					
					//for(int q = 0; q < 2; ++q)
					//{
					//	REPORT_STR("new edges["<<q<<"]" << new_edges[q].LocalID() << " nodes " << new_edges[q].getBeg().LocalID() << "," << new_edges[q].getEnd().LocalID() << " level " << level[new_edges[q]]);
					//}
					
					//if( !Element::CheckConnectivity(m) ) std::cout << __FILE__ << ":" << __LINE__ << " broken connectivity" << std::endl;
				}
			}
			EXIT_BLOCK();
			//8.split all faces of the current schedule, using hanging nodes on edges
			ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->FaceLastLocalID(); ++it) if( m->isValidFace(it) )
			{
				Face f = m->FaceByLocalID(it);
				if( !f.Hidden() && indicator[f] == schedule_counter )
				{
					//connect face center to hanging nodes of the face
					Storage::reference_array face_hanging_nodes = hanging_nodes[f];
					//remember adjacent cells that should get information about new hanging node
					//and new hanging edges
					ElementArray<Cell> face_cells = f.getCells();
					//create node at face center
					//f->Centroid(xyz);
					for(int d = 0; d < 3; ++d) xyz[d] = 0.0;
					for(Storage::reference_array::size_type kt = 0; kt < face_hanging_nodes.size(); ++kt)
						for(int d = 0; d < 3; ++d) xyz[d] += face_hanging_nodes[kt].getAsNode().Coords()[d];
					for(int d = 0; d < 3; ++d) xyz[d] /= (Storage::real)face_hanging_nodes.size();
					//todo: request transformation of node location according to geometrical model
					//create middle node
					Node n = m->CreateNode(xyz);
					//set increased level for the new node
					level[n] = level[f]+1;
					//for each cell provide link to new hanging node
					for(ElementArray<Face>::size_type kt = 0; kt < face_cells.size(); ++kt)
						hanging_nodes[face_cells[kt]].push_back(n);
					ElementArray<Node> edge_nodes(m,2); //to create new edges
					ElementArray<Edge> hanging_edges(m,face_hanging_nodes.size());
					edge_nodes[0] = n;
					for(Storage::reference_array::size_type kt = 0; kt < face_hanging_nodes.size(); ++kt)
					{
						edge_nodes[1] = face_hanging_nodes[kt].getAsNode();
						hanging_edges[kt] = m->CreateEdge(edge_nodes).first;
						//set increased level for new edges
						level[hanging_edges[kt]] = level[f]+1;
					}
					//split the face by these edges
					ElementArray<Face> new_faces = Face::SplitFace(f,hanging_edges,0);
					//set increased level to new faces
					for(ElementArray<Face>::size_type kt = 0; kt < new_faces.size(); ++kt)
						level[new_faces[kt]] = level[f]+1;
				}
			}
			EXIT_BLOCK();
			//this tag helps recreate internal face
			TagReferenceArray internal_face_edges = m->CreateTag("INTERNAL_FACE_EDGES",DATA_REFERENCE,NODE,NODE,4);
			//this marker helps detect edges of current cell only
			MarkerType mark_cell_edges = m->CreateMarker();
			//this marker helps detect nodes hanging on edges of unrefined cell
			MarkerType mark_hanging_nodes = m->CreateMarker();
			//9.split all cells of the current schedule
			ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( !c.Hidden() && indicator[c] == schedule_counter )
				{
					Storage::reference_array cell_hanging_nodes = hanging_nodes[c]; //nodes to be connected
					//create node at cell center
					for(int d = 0; d < 3; ++d) xyz[d] = 0.0;
					for(Storage::reference_array::size_type kt = 0; kt < cell_hanging_nodes.size(); ++kt)
						for(int d = 0; d < 3; ++d) xyz[d] += cell_hanging_nodes[kt].getAsNode().Coords()[d];
					for(int d = 0; d < 3; ++d) xyz[d] /= (Storage::real)cell_hanging_nodes.size();
					//c->Centroid(xyz);
					//todo: request transformation of node location according to geometrical model
					//create middle node
					Node n = m->CreateNode(xyz);
					//set increased level for the new node
					level[n] = level[c]+1;
					//retrive all edges of current face to mark them
					ElementArray<Edge> cell_edges = c.getEdges();
					//mark all edges so that we can retive them later
					cell_edges.SetMarker(mark_cell_edges);
					//connect face center to centers of faces by edges
					ElementArray<Node> edge_nodes(m,2);
					ElementArray<Edge> edges_to_faces(m,cell_hanging_nodes.size());
					edge_nodes[0] = n;
					for(Storage::reference_array::size_type kt = 0; kt < cell_hanging_nodes.size(); ++kt)
					{
						assert(cell_hanging_nodes[kt].isValid());
						//todo: unmark hanging node on edge if no more cells depend on it
						edge_nodes[1] = cell_hanging_nodes[kt].getAsNode();
						edges_to_faces[kt] = m->CreateEdge(edge_nodes).first;
						//set increased level for new edges
						level[edges_to_faces[kt]] = level[c]+1;
						//for each node other then the hanging node of the face
						//(this is hanging node on the edge)
						//we record a pair of edges to reconstruct internal faces
						ElementArray<Edge> hanging_edges = cell_hanging_nodes[kt].getEdges(mark_cell_edges,0);
						for(ElementArray<Edge>::size_type lt = 0; lt < hanging_edges.size(); ++lt)
						{
							//get hanging node on the edge
							assert(hanging_edges[lt].getBeg() == cell_hanging_nodes[kt] || hanging_edges[lt].getEnd() == cell_hanging_nodes[kt]);
							Node v = hanging_edges[lt].getBeg() == cell_hanging_nodes[kt]? hanging_edges[lt].getEnd() : hanging_edges[lt].getBeg();
							//mark so that we can collect all of them
							v.SetMarker(mark_hanging_nodes);
							//fill the edges
							Storage::reference_array face_edges = internal_face_edges[v];
							//fill first two in forward order
							//this way we make a closed loop
							assert(face_edges[0] == InvalidElement() || face_edges[2] == InvalidElement());
							if( face_edges[0] == InvalidElement() )
							{
								face_edges[0] = edges_to_faces[kt];
								face_edges[1] = hanging_edges[lt];
							}
							else //last two in reverse
							{
								assert(face_edges[2] ==InvalidElement());
								face_edges[2] = hanging_edges[lt];
								face_edges[3] = edges_to_faces[kt];
							}
						}
					}
					//remove marker from cell edges
					cell_edges.RemMarker(mark_cell_edges);
					//now we have to create internal faces
					ElementArray<Node> edge_hanging_nodes = c.getNodes(mark_hanging_nodes,0);
					ElementArray<Face> internal_faces(m,edge_hanging_nodes.size());
					//unmark hanging nodes on edges
					edge_hanging_nodes.RemMarker(mark_hanging_nodes);
					for(ElementArray<Node>::size_type kt = 0; kt < edge_hanging_nodes.size(); ++kt)
					{
						//create a face based on collected edges
						Storage::reference_array face_edges = internal_face_edges[edge_hanging_nodes[kt]];
						assert(face_edges[0].isValid());
						assert(face_edges[1].isValid());
						assert(face_edges[2].isValid());
						assert(face_edges[3].isValid());
						internal_faces[kt] = m->CreateFace(ElementArray<Edge>(m,face_edges.begin(),face_edges.end())).first;
						//set increased level
						level[internal_faces[kt]] = level[c]+1;
						//clean up structure, so that other cells can use it
						edge_hanging_nodes[kt].DelData(internal_face_edges);
					}
					//if( c.GlobalID() == 228 )
					//{
					//	double cnt[3];
					//	c.Centroid(cnt);
					//	std::cout << "Split CELL:" << c.LocalID() << " " << c.GlobalID() << " " << Element::StatusName(c.GetStatus()) << " " << cnt[0] << " " << cnt[1] << " " << cnt[2] << std::endl;
					//	
					//}
					//split the cell
					ElementArray<Cell> new_cells = Cell::SplitCell(c,internal_faces,0);
					std::sort(new_cells.begin(),new_cells.end(),Mesh::CentroidComparator(m));
					//retrive parent set
					ElementSet parent(m,parent_set[c]);
					//create set corresponding to old coarse cell
					Storage::real cnt[3];
					c.Centroid(cnt);
					std::stringstream set_name;
					//set_name << parent.GetName() << "_C" << c.GlobalID(); //rand may be unsafe
					if( parent == root )
						set_name << "AM_R" << set_id[c];
					else
						set_name << parent.GetName() << "C" << set_id[c];
					//set_name << base64_encode_((unsigned char *)cnt,3*sizeof(double)/sizeof(unsigned char));
					
					ElementSet check_set = m->GetSet(set_name.str());
					if( check_set.isValid() )
					{
						std::cout << rank << " set " << set_name.str() << " for cell " << c.GlobalID() << " " << Element::StatusName(c.GetStatus()) << " already exists" << std::endl;
						if( check_set->HaveParent() )
							std::cout << rank << " parent is " << check_set->GetParent()->GetName() << " cell parent is " << parent.GetName() << std::endl;
						std::cout << rank << " Elements: ";
						for(ElementSet::iterator it = check_set.Begin(); it != check_set.End(); ++it)
							std::cout << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << "," << it->GlobalID() << "," << Element::StatusName(c.GetStatus()) << "," << level[*it] << " ";
						std::cout << std::endl;
						exit(-1);
					}
					
					ElementSet cell_set = m->CreateSetUnique(set_name.str()).first;
					//cell_set->SetExchange(ElementSet::SYNC_ELEMENTS_ALL);
					level[cell_set] = level[c]+1;
					//set up increased level for the new cells
					for(ElementArray<Cell>::size_type kt = 0; kt < new_cells.size(); ++kt)
					{
						set_id[new_cells[kt]] = kt;
						level[new_cells[kt]] = level[c]+1;
						cell_set.PutElement(new_cells[kt]);
						parent_set[new_cells[kt]] = cell_set.GetHandle();
					}
					/*
					if( check_set.isValid() )
					{
						std::cout << rank << " Elements: ";
						for(ElementSet::iterator it = check_set.Begin(); it != check_set.End(); ++it)
							std::cout << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << "," << it->GlobalID() << "," << Element::StatusName(c.GetStatus()) << "," << level[*it] << " ";
						std::cout << std::endl;
					}
					*/
					//if( !cell_set->HaveParent() )
					parent.AddChild(cell_set);
					//else assert(cell_set->GetParent() == parent);
					//increment number of refined cells
					ret++;
				}
			}
			EXIT_BLOCK();
			m->ReleaseMarker(mark_hanging_nodes);
			m->ReleaseMarker(mark_cell_edges);
			m->DeleteTag(internal_face_edges);
			//10.jump to later schedule, and go to 7.
			schedule_counter--;
		}

		//free created tag
		m->DeleteTag(indicator,FACE|EDGE);


		//11. Restore parallel connectivity, global ids
		m->ResolveModification();
		//m->SynchronizeMarker(m->NewMarker(),CELL|FACE|EDGE|NODE,SYNC_BIT_OR);
        //ExchangeGhost(3,NODE); // Construct Ghost cells in 2 layers connected via nodes
		//12. Let the user update their data
		//todo: call back function for New() cells
		if( model ) model->Adaptation(*m);
		//13. Delete old elements of the mesh
		m->ApplyModification();
		
		//m->ExchangeGhost(1,NODE,m->NewMarker());
		//14. Done
        //cout << rank << ": Before end " << std::endl;
		m->EndModification();
		EXIT_BLOCK();
		
		
		//m->Save("after_refine"+std::to_string(fi)+".pvtk");
		//std::cout << "Save after_refine"+std::to_string(fi)+".pvtk" << std::endl;
		fi++;
		
		//ExchangeData(hanging_nodes,CELL | FACE,0);
        //m->ResolveSets();

        //m->BeginModification();
        //    m->ExchangeGhost(1,NODE,marker_new); // Construct Ghost cells in 2 layers connected via nodes
        //    m->ReleaseMarker(marker_new,CELL|FACE|EDGE|NODE);
		//m->ApplyModification();
    	//m->EndModification();
    	//PrintSet();
    	//m->ResolveSets();
		
		//restore face orientation
		//BUG: bad orientation not fixed automatically
		
		/*
		ENTER_BLOCK();
		int nfixed = 0;
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
			if( !it->CheckNormalOrientation() )
			{
				it->FixNormalOrientation();
				nfixed++;
			}
			//std::cout << "Face " << it->LocalID() << " oriented incorrectly " << std::endl;
		if( nfixed ) REPORT_STR(rank << " fixed " << nfixed << " faces");
		EXIT_BLOCK();
		 */
		
		/*
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		{
			if( parent_set[*it] == InvalidHandle() )
				std::cout << m->GetProcessorRank() << " parent set not valid on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << std::endl;
			else if( GetHandleElementType(parent_set[*it]) != ESET )
				std::cout << m->GetProcessorRank() << " parent set is something else " << ElementTypeName(GetHandleElementType(parent_set[*it])) << ":" << GetHandleID(parent_set[*it]) << " on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << std::endl;
		}
		*/
		//std::cout << rank << " check parent_set at end" << std::endl;
		//CheckParentSet();

        //ExchangeData(hanging_nodes,CELL | FACE,0);
        //cout << rank << ": After end " << std::endl;
		//reorder element's data to free up space
		ENTER_BLOCK();
		m->ReorderEmpty(CELL|FACE|EDGE|NODE);
		EXIT_BLOCK();
		//return number of refined cells
		call_counter++;
		ret = m->Integrate(ret);
		REPORT_VAL("ret ",ret)
        EXIT_FUNC();
		return ret != 0;
	}

    void OperationMin(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
    {
        int value  =  (int)*((int*)data);


        if (value < element->Integer(tag))
        {
            element->Integer(tag) = value;
        }
        (void)size;
    }
	

	bool AdaptiveMesh::Coarse(TagInteger & indicator)
	{
		std::string file;
		ENTER_FUNC();
        //return false;
		static int call_counter = 0;
		//return number of coarsened cells
		int ret = 0;
		//initialize tree structure
		ENTER_BLOCK();
		PrepareSet();
		EXIT_BLOCK();

		//m->Save("before_coarse"+std::to_string(fi)+".pvtk");
		//std::cout << "Save before_coarse"+std::to_string(fi)+".pvtk" << std::endl;
		
		ENTER_BLOCK();
		//m->CheckCentroids(__FILE__,__LINE__);
		m->ResolveSets();
		m->ExchangeData(parent_set,CELL,0);
		//m->CheckCentroids(__FILE__,__LINE__);
		m->ResolveSets();
		//m->CheckCentroids(__FILE__,__LINE__);
		//m->ExchangeData(parent_set,CELL,0);
		//m->CheckCentroids(__FILE__,__LINE__);
		//m->ExchangeData(hanging_nodes,CELL | FACE,0);
		//m->CheckCentroids(__FILE__,__LINE__);
		//CheckParentSet();
		EXIT_BLOCK();

		//m->Save("before_coarse_parent"+std::to_string(fi)+".pvtk");
		//std::cout << "Save before_coarse_parent"+std::to_string(fi)+".pvtk" << std::endl;
		
		//ENTER_BLOCK();
		//SynchronizeIndicated(indicator);
		//EXIT_BLOCK();
		
		
		
		int schedule_counter = 1; //indicates order in which refinement will be scheduled
		int scheduled = 1, unscheduled = 0; //indicates that at least one element was scheduled on current sweep

		ENTER_BLOCK();
		//TagInteger coarsened = CreateTag("COARSENED",DATA_INTEGER,CELL,NONE,1);
		TagInteger coarse_indicator = m->CreateTag("COARSE_INDICATOR",DATA_INTEGER,ESET|EDGE,ESET,1); //used to find on fine cells indicator on coarse cells
		//0. Extend indicator for sets, edges and faces
		indicator = m->CreateTag(indicator.GetTagName(),DATA_INTEGER,ESET|FACE|EDGE,ESET,1);
		std::vector<Tag> indicators;
		indicators.push_back(coarse_indicator);
		indicators.push_back(indicator);
		//int mi = 0;
		//file = "w"+std::to_string(mi++)+".pvtk";
		//m->Save(file);
		//std::cout << "Save " << file << " " << __FILE__ << ":" << __LINE__ << std::endl;
		while(scheduled || unscheduled)
		{
			REPORT_VAL("scheduled",scheduled);
			REPORT_VAL("unscheduled",unscheduled);
			// rules
			// a) If there is adjacent finer edge that is not marked for coarsening
			// then this cell should not be coarsened
			// b) If there is adjacent coarser cell, then this cell should be coarsened
			// first
			//0.Communicate indicator - it may be not synced
			m->ExchangeData(indicator,CELL,0);
			//1. Mark each adjacent face/edge for coarsement schedule
			// problem: should mark so that if every adjacent cell is coarsened
			// then adjacent face/edge are also coarsened
			ENTER_BLOCK();
			for(ElementType etype = EDGE; etype <= FACE; etype = NextElementType(etype))
			{
				//for(Storage::integer it = 0; it < LastLocalID(etype); ++it) if( isValidElement(etype,it) )
				//	indicator[ElementByLocalID(etype,it)] = 0;
				for(Storage::integer it = 0; it < m->LastLocalID(etype); ++it) if( m->isValidElement(etype,it) )
				{
					Element e = m->ElementByLocalID(etype,it);
					ElementArray<Cell> adj = e.getCells();
					indicator[e] = INT_MAX;
					for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
						if( level[e] == level[adj[kt]]) indicator[e] = std::min(indicator[e],indicator[adj[kt]]);
					//if (indicator[e] == INT_MAX) cout << rank << ": " << ElementTypeName(e.GetElementType()) << e.GlobalID() << endl;
					//assert(indicator[e] != INT_MAX);
				}
			}
			EXIT_BLOCK();
			//2.Communicate indicator on faces and edges
			m->ReduceData(indicator,FACE|EDGE,0,ReduceMin);
			m->ExchangeData(indicator,FACE|EDGE,0);
			//3.If there is adjacent finer edge that are not marked for coarsening
			// then this cell should not be coarsened
			unscheduled = scheduled = 0;
			ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( indicator[c] )
				{
					ElementArray<Edge> edges = c.getEdges();
					for(ElementArray<Edge>::size_type kt = 0; kt < edges.size(); ++kt)
					{
						if( level[edges[kt]] > level[c] && indicator[edges[kt]] == 0 )
						{
							indicator[c] = 0;
							unscheduled++;
						}
					}
				}
			}
			EXIT_BLOCK();
			//file = "w"+std::to_string(mi++)+".pvtk";
			//m->Save(file);
			//std::cout << "Save " << file << " " << __FILE__ << ":" << __LINE__ << std::endl;
			//4. Propogate coarsement info over set tree to detect valid coarsenings.
			// go down over sets, if set does not have children and all of the cells
			// of the set are marked for coarsening, then mark the set for coarsement
			// otherwise unmark.
			// Unmark all cells that are not to be coarsened
			ENTER_BLOCK();
			//int v1 = 0, v2 = 0, v3 = 0;
			//std::fstream fout("set"+std::to_string(m->GetProcessorRank())+".txt",std::ios::out);
			for(Storage::integer it = 0; it < m->EsetLastLocalID(); ++it) if( m->isValidElementSet(it) )
			{
				ElementSet set = m->EsetByLocalID(it);
				if( set.GetName().substr(0,3) != "AM_" ) continue;
				if( !set.HaveChild() )
				{
					indicator[set] = INT_MAX;
					coarse_indicator[set] = 0;
					for(ElementSet::iterator jt = set.Begin(); jt != set.End(); ++jt)
					{
						assert(parent_set[*jt] == set.GetHandle());
						indicator[set] = std::min(indicator[set],indicator[*jt]);
						coarse_indicator[set] = std::max(coarse_indicator[set],indicator[*jt]);
					}
				}
			}
			EXIT_BLOCK();
			m->ReduceData(indicator,ESET,0,ReduceMin);
			m->ReduceData(coarse_indicator,ESET,0,ReduceMax);
			m->ExchangeData(indicators,ESET,0);
			ENTER_BLOCK()
			for(Storage::integer it = 0; it < m->EsetLastLocalID(); ++it) if( m->isValidElementSet(it) )
			{
				ElementSet set = m->EsetByLocalID(it);
				if( set.GetName().substr(0,3) != "AM_" ) continue;
				if( set.HaveChild() )
				{
					//we have finer elements to be coarsened first,
					//so unschedule all the cells
					for(ElementSet::iterator jt = set.Begin(); jt != set.End(); ++jt) if( indicator[*jt] != 0 )
					{
						indicator[*jt] = 0;
						unscheduled++;
					}
				}
				else
				{
					//check max and min of indicator over set matches
					if( coarse_indicator[set] != indicator[set] )
					{
						//we have to change schedule for some elements,
						//since elements of the set are scheduled in different order
						if( indicator[set] == 0 )
						{
							//some elements were denied from coarsening, so we have to deny the rest
							for(ElementSet::iterator jt = set.Begin(); jt != set.End(); ++jt) if( indicator[*jt] != 0 )
							{
								indicator[*jt] = 0;
								unscheduled++;
							}
						}
						else
						{
							for(ElementSet::iterator jt = set.Begin(); jt != set.End(); ++jt)
							{
								//some elements were scheduled later, so we have to reschedule the rest
								for(ElementSet::iterator jt = set.Begin(); jt != set.End(); ++jt) if( indicator[*jt] != coarse_indicator[set] )
								{
									indicator[*jt] = coarse_indicator[set];
									unscheduled++;
								}
							}
						}
					}
				}
			}
			EXIT_BLOCK();
			
			/*
			 ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
                //if (!isValidElement(c.GetHandle())) continue;
				if( indicator[c] )
				{
					ElementSet parent(m,parent_set[c]);
					//intermediate cell may not be coarsened
					//root set may not have coarsening cells
					if( parent.HaveChild() )//|| !parent.HaveParent() )
					{
						//PrintSet(fout,parent);
						indicator[parent] = 0;
						indicator[c] = 0;
						unscheduled++;
						//v1++;
					}
					else
					{
						Storage::integer schedule_first = 0;
						bool check = true;
						//check that all elements of the set are to be coarsened
						for(ElementSet::iterator it = parent.Begin(); it != parent.End(); ++it)
						{
							check &= (indicator[it->self()] != 0);
							schedule_first = std::max(schedule_first,indicator[it->self()]);
						}
						if(!check)
						{
							indicator[c] = 0;
							unscheduled++;
							//v2++;
						}
						else if( indicator[c] != schedule_first )
						{
							indicator[c] = schedule_first;
							unscheduled++;
							//v3++;
						}
					}
				}
			}
			 EXIT_BLOCK();
			 */
			//fout.close();
			//std::cout << "v1 " << v1 << " v2 " << v2 << " v3 " << v3 << std::endl;
			
			//file = "w"+std::to_string(mi++)+".pvtk";
			//m->Save(file);
			//std::cout << "Save " << file << " " << __FILE__ << ":" << __LINE__ << std::endl;
			//if( fi == 5 && mi == 3 ) exit(-1);
			//5.If there is an adjacent coarser element to be refined, then
			//   this one should be scheduled to be refined first
			//a) clean up coarse indicator tag
			ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->EdgeLastLocalID(); ++it) if( m->isValidEdge(it) )
				coarse_indicator[m->EdgeByLocalID(it)] = 0;
			EXIT_BLOCK();
			//b) each cell mark it's finer edges with cell's schedule
			ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( indicator[c] )
				{
					ElementArray<Element> adj = c.getAdjElements(EDGE);
					for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
					{
						if( level[adj[kt]] > level[c] ) //only finer edges
							coarse_indicator[adj[kt]] = std::max(coarse_indicator[adj[kt]],indicator[c]);
					}
				}
			}
			EXIT_BLOCK();
			//file = "w"+std::to_string(mi++)+".pvtk";
			//m->Save(file);
			//std::cout << "Save " << file << " " << __FILE__ << ":" << __LINE__ << std::endl;
			//c) data reduction to get maximum over mesh partition
			m->ReduceData(coarse_indicator,EDGE,0,ReduceMax);
			m->ExchangeData(coarse_indicator,EDGE,0);
			//d) look from cells if any edge is coarsened earlier
			ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( indicator[c] )
				{
					ElementArray<Element> adj = c.getAdjElements(EDGE);
					for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
					{
						if( level[c] == level[adj[kt]] && //do not look from coarse cell onto finer edge
						    indicator[c] <= coarse_indicator[adj[kt]])
						{
							indicator[c] = coarse_indicator[adj[kt]]+1;
							scheduled++;
						}
					}
				}
			}
			EXIT_BLOCK();
			m->ExchangeData(indicator,CELL|FACE|EDGE,0);
			//file = "w"+std::to_string(mi++)+".pvtk";
			//m->Save(file);
			//std::cout << "Save " << file << " " << __FILE__ << ":" << __LINE__ << std::endl;
			//5.Go back to 1 until no new elements scheduled
			scheduled = m->Integrate(scheduled);
			unscheduled = m->Integrate(unscheduled);
			if( scheduled ) schedule_counter++;
		}
		//cleanup
		coarse_indicator = m->DeleteTag(coarse_indicator);
		EXIT_BLOCK();
		
		
		//Order exchange of elemnets of the sets that are to be coarsened
		ENTER_BLOCK();
		for(Storage::integer it = 0; it < m->EsetLastLocalID(); ++it) if( m->isValidElementSet(it) )
		{
			ElementSet set = m->EsetByLocalID(it);
			if( indicator[set] != 0 )
				set.SynchronizeSetElements();
		}
		EXIT_BLOCK();
		m->ExchangeMarked();

		//CheckParentSet();
		
		ENTER_BLOCK();
		//m->CheckCentroids(__FILE__,__LINE__);
		m->ExchangeData(parent_set,CELL,0);
		//m->CheckCentroids(__FILE__,__LINE__);
		m->ResolveSets();
		//m->CheckCentroids(__FILE__,__LINE__);
		//m->ExchangeData(parent_set,CELL,0);
		//m->CheckCentroids(__FILE__,__LINE__);
		m->ExchangeData(hanging_nodes,CELL | FACE,0);
		//m->CheckCentroids(__FILE__,__LINE__);
		
		EXIT_BLOCK();


		ENTER_BLOCK();
		for(Storage::integer it = 0; it < m->EsetLastLocalID(); ++it) if( m->isValidElementSet(it) )
		{
			ElementSet set = m->EsetByLocalID(it);
			if( indicator[set] != 0 )
				set.SynchronizeSetParents();
		}
		EXIT_BLOCK();
		m->ExchangeMarked();
		
		CheckParentSet();
		
		//std::fstream fout("sets"+std::to_string(m->GetProcessorRank())+".txt",std::ios::out);
		//for(Mesh::iteratorSet it = m->BeginSet(); it != m->EndSet(); ++it)
		//	PrintSet(fout,it->self());
		
		
		
		//m->Save("unschdind"+std::to_string(fi)+".pvtk");
		//std::cout << "Save unschdind"+std::to_string(fi)+".pvtk" << std::endl;
		//Make schedule which elements should be refined earlier.
		ENTER_BLOCK();
		m->BeginModification();
		while(schedule_counter)
		{
			CheckParentSet();
			//CheckParentSet();
			//fout << "schedule_counter " << schedule_counter << std::endl;
			//unite cells
			//should find and set hanging nodes on faces
			//find single node at the center, all other nodes,
			//adjacent over edge are hanging nodes
			ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( !c.Hidden() && indicator[c] == schedule_counter )
				{
					double x[3];
					c.Centroid(x);
					//this set contains all the cells to be united
					ElementSet parent(m,parent_set[c]);
					ElementArray<Cell> unite_cells(m,parent.Size());
					//unmark indicator to prevent coarsement with next element
					Storage::integer kt = 0;
					for(ElementSet::iterator jt = parent.Begin(); jt != parent.End(); ++jt)
					{
						unite_cells[kt++] = jt->getAsCell();
						indicator[jt->self()] = 0; //prevent algorithm from visiting again
					}
					//fout << "parent set " << parent.GetName() << " size " << parent.Size() << " cell " << c.GlobalID() << " " << x[0] << " " << x[1] << " " << x[2] << std::endl;
					//fout << "cells(" << unite_cells.size() << "):" << std::endl;
					//for(kt = 0; kt < unite_cells.size(); ++kt)
					//{
					//	unite_cells[kt].Centroid(x);
					//	fout << unite_cells[kt].GlobalID();
					//	fout << " parent " << ElementSet(m,parent_set[unite_cells[kt]]).GetName();
					//	fout << " " << x[0] << " " << x[1] << " " << x[2];
					//	fout << " lvl " << level[unite_cells[kt]];
					//	fout << std::endl;
					//}
					//find a node common to all the cells
					ElementArray<Node> center_node = unite_cells[0].getNodes();
					for(kt = 1; kt < unite_cells.size(); ++kt)
						center_node.Intersect(unite_cells[kt].getNodes());
					
					//fout << "nodes(" << center_node.size() << "):" << std::endl;
					//for(kt = 0; kt < center_node.size(); ++kt)
					//	fout << center_node[kt].Coords()[0] << " " << center_node[kt].Coords()[1] << " " << center_node[kt].Coords()[2] << std::endl;
					//fout << "child sets: ";
					//for(ElementSet chld = parent.GetChild(); chld.isValid(); chld = chld.GetSibling())
					//	fout << " " << chld.GetName();
					//fout << std::endl;
					//only one should be found
					if( center_node.size() != 1 )
					{
						double x[3];
						std::cout << "call_counter " << call_counter << " schedule_counter " << schedule_counter << std::endl;
						c.Centroid(x);
						std::cout << "parent set " << parent.GetName() << " size " << parent.Size() << " cell " << c.GlobalID() << " " << x[0] << " " << x[1] << " " << x[2] << std::endl;
						std::cout << "cells(" << unite_cells.size() << "):" << std::endl;
						for(kt = 0; kt < unite_cells.size(); ++kt)
						{
							unite_cells[kt].Centroid(x);
							std::cout << unite_cells[kt].GlobalID() << " lid " << unite_cells[kt].LocalID();
							std::cout << " parent " << ElementSet(m,parent_set[unite_cells[kt]]).GetName();
							std::cout << " " << x[0] << " " << x[1] << " " << x[2];
							std::cout << std::endl;
						}
						std::cout << "nodes(" << center_node.size() << "):" << std::endl;
						for(kt = 0; kt < center_node.size(); ++kt)
							std::cout << center_node[kt].Coords()[0] << " " << center_node[kt].Coords()[1] << " " << center_node[kt].Coords()[2] << std::endl;
						std::cout << "child sets: ";
						for(ElementSet chld = parent.GetChild(); chld.isValid(); chld = chld.GetSibling())
							std::cout << " " << chld.GetName();
						std::cout << std::endl;
					}
					assert(center_node.size() == 1);
					ElementArray<Node> hanging = center_node[0].BridgeAdjacencies2Node(EDGE);
					Cell v = Cell::UniteCells(unite_cells,0);
					//connect hanging nodes to the cell
					assert(hanging_nodes[v].size() == 0);
					for(ElementArray<Node>::size_type kt = 0; kt < hanging.size(); ++kt)
						hanging_nodes[v].push_back(hanging[kt]);
					//set new parent
					parent_set[v] = parent.GetParent().GetHandle();
					//add cell to parent set
					ElementSet(m,parent_set[v]).PutElement(v);
					//set level for new cell
					level[v] = level[c]-1;
					
					v.Centroid(x);
					//fout << v.GlobalID() << " lid " << v.LocalID();
					//fout << " parent " << ElementSet(m,parent_set[v]).GetName();
					//fout << " " << x[0] << " " << x[1] << " " << x[2];
					//fout << " lvl " << level[v];
					//fout << std::endl;
					//delete set that contained cells
					//tree structure should be resolved on ApplyModification
					//fout << "delete set " << parent.GetName() << std::endl;
					parent.DeleteSet();
					//increment number of coarsened cells
					ret++;
				}
			}
			EXIT_BLOCK();
			//unite faces
			//should find and set hanging nodes on edges
			//find single node at the center, all other nodes,
			//adjacent over edge of the face are hanging nodes
			int numcoarsened = 0;
			ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->FaceLastLocalID(); ++it) if( m->isValidFace(it) )
			{
				Face f = m->FaceByLocalID(it);
				if( !f.Hidden() && indicator[f] == schedule_counter )
				{
					//one (or both) of the adjacent cells were coarsened and has lower level
					bool visited = false;
					ElementArray<Cell> cells = f.getCells();
					for(ElementArray<Cell>::size_type kt = 0; kt < cells.size(); ++kt)
					{
						assert(level[cells[kt]] < level[f]);
					}
					for(ElementArray<Cell>::size_type kt = 0; kt < cells.size(); ++kt)
					{
						if( level[cells[kt]] < level[f] )
						{
							//cell has one hanging node in common with current face
							ElementArray<Node> nodes = f.getNodes();
							Storage::reference_array search_hanging = hanging_nodes[cells[kt]];
							nodes.Intersect(search_hanging.data(),search_hanging.size());
							assert(nodes.size() == 1);
							//faces that hanging node shares with the cell are
							//those to be united
							ElementArray<Face> unite_faces = cells[kt].getFaces();
							unite_faces.Intersect(nodes[0].getFaces());
							//unmark faces to prevent visit
							for(ElementArray<Face>::size_type lt = 0; lt < unite_faces.size(); ++lt)
								indicator[unite_faces[lt]] = 0;
							//nodes connected by edges to hanging node and
							//common to the cell are hanging nodes on edges
							ElementArray<Node> hanging = cells[kt].getNodes();
							hanging.Intersect(nodes[0].BridgeAdjacencies(EDGE,NODE));
							//unite faces
							Face v = Face::UniteFaces(unite_faces,0);
							//set level for new face
							level[v] = level[f]-1;
							//connect new face to hanging nodes
							for(ElementArray<Node>::size_type lt = 0; lt < hanging.size(); ++lt)
								hanging_nodes[v].push_back(hanging[lt]);
							visited = true;
							numcoarsened++;
							break; //no need to visit the other cell
						}
					}
					assert(visited);
				}
			}
			EXIT_BLOCK();
			//unite edges
			ENTER_BLOCK();
			for(Storage::integer it = 0; it < m->EdgeLastLocalID(); ++it) if( m->isValidEdge(it) )
			{
				Edge e = m->EdgeByLocalID(it);
				if( !e.Hidden() && indicator[e] == schedule_counter )
				{
					//at least one face must have lower level
					bool visited = false;
					ElementArray<Face> faces = e.getFaces();
					for(ElementArray<Face>::size_type kt = 0; kt < faces.size(); ++kt)
					{
						if( level[faces[kt]] < level[e] )
						{
							//face has one hanging node in common with current edge
							ElementArray<Node> nodes = e.getNodes();
							Storage::reference_array search_hanging = hanging_nodes[faces[kt]];
							nodes.Intersect(search_hanging.data(),search_hanging.size());
							assert(nodes.size() == 1);
							//edges that hanging node shares with the face are those to
							//be united
							ElementArray<Edge> unite_edges = faces[kt].getEdges();
							unite_edges.Intersect(nodes[0].getEdges());
							//unmark edges to prevent visit
							for(ElementArray<Edge>::size_type lt = 0; lt < unite_edges.size(); ++lt)
								indicator[unite_edges[lt]] = 0;
							//unite edges
							Edge v = Edge::UniteEdges(unite_edges,0);
							//set level for new edge
							level[v] = level[e]-1;
							visited = true;
							break; //no need to visit any other face
						}
					}
					assert(visited);
				}
			}
			EXIT_BLOCK();
			/*
			for(Storage::integer it = 0; it < NodeLastLocalID(); ++it) if( isValidNode(it) )
			{
				Node e = NodeByLocalID(it);
				if( !e.Hidden() )
				{
					int my_level = -1;
					ElementArray<Edge> edges = e.getEdges();
					for(ElementArray<Edge>::iterator kt = edges.begin(); kt != edges.end(); ++kt)
						my_level = std::max(my_level,level[*kt]);
					level[e] = my_level;
				}
			}
			*/
			//jump to later schedule
			schedule_counter--;
		}
		//free created tag
		m->DeleteTag(indicator,FACE|EDGE);
		//todo:
		m->ResolveModification();
		//todo:
		//let the user update their data
		if( model ) model->Adaptation(*m);
		m->ApplyModification();
		//done
		m->EndModification();
		EXIT_BLOCK();
		//fout.close();
		
		//m->Save("after_coarse"+std::to_string(fi)+".pvtk");
		//std::cout << "Save after_coarse"+std::to_string(fi)+".pvtk" << std::endl;
		//exit(-1);
		
		ENTER_BLOCK();
		//m->CheckCentroids(__FILE__,__LINE__);
		//CheckCentroids();
		//cleanup null links to hanging nodes
		for(ElementType etype = FACE; etype <= CELL; etype = NextElementType(etype))
		{
			for(Storage::integer it = 0; it < m->LastLocalID(etype); ++it) if( m->isValidElement(etype,it) )
			{
				Storage::reference_array arr = hanging_nodes[m->ElementByLocalID(etype,it)];
				Storage::reference_array::size_type jt = 0;
				for(Storage::reference_array::size_type kt = 0; kt < arr.size(); ++kt)
					if( arr[kt] != InvalidElement() ) arr[jt++] = arr[kt];
				arr.resize(jt);
			}
		}
		EXIT_BLOCK();
		//m->ResolveSets();
		
		//cleanup null links in sets
		CleanupSets(root);
		
		//CheckParentSet();
		
		
		//restore face orientation
		//BUG: bad orientation not fixed automatically
		
		/*
		ENTER_BLOCK();
		int nfixed = 0;
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
			if( !it->CheckNormalOrientation() )
			{
				it->FixNormalOrientation();
				nfixed++;
			}
			//std::cout << "Face " << it->LocalID() << " oriented incorrectly " << std::endl;
		if( nfixed ) REPORT_STR("fixed " << nfixed << " faces");
		EXIT_BLOCK();
		 */
		
		//reorder element's data to free up space
		ENTER_BLOCK();
		m->ReorderEmpty(CELL|FACE|EDGE|NODE|ESET);
		EXIT_BLOCK();
		
		call_counter++;
		
		ret = m->Integrate(ret);
		REPORT_VAL("ret ",ret)
        EXIT_FUNC();
		return ret != 0;
	}
	
#if defined(USE_PARTITIONER)
	void AdaptiveMesh::SetNewOwner(ElementSet set, TagInteger owner)
	{
		if( set.HaveChild() )
		{
			if( set.Empty() ) owner[set] = INT_MAX;
			for(ElementSet chld = set.GetChild(); chld.isValid(); chld = chld.GetSibling())
			{
				SetNewOwner(chld,owner);
				owner[set] = std::min(owner[set],owner[chld]);
			}
		}
		//for(ElementSet::iterator it = set.Begin(); it != set.End(); ++it)
		//	owner[set] = std::min(owner[set],owner[*it]);
	}
	void AdaptiveMesh::SetNewProcs(ElementSet set, TagIntegerArray procs)
	{
		std::vector<int> tmp;
		for(ElementSet chld = set.GetChild(); chld.isValid(); chld = chld.GetSibling())
		{
			SetNewProcs(chld,procs);
			Storage::integer_array pm = procs[set];
			Storage::integer_array pc = procs[chld];
			tmp.resize(pm.size()+pc.size());
			tmp.resize(std::set_union(pm.begin(),pm.end(),pc.begin(),pc.end(),tmp.begin())-tmp.begin());
			pm.replace(pm.begin(),pm.end(),tmp.begin(),tmp.end());
		}
	}
	void AdaptiveMesh::RestoreParent(ElementSet set)
	{
		for(ElementSet chld = set.GetChild(); chld.isValid(); chld = chld.GetSibling())
			RestoreParent(chld);
		for(ElementSet::iterator it = set.Begin(); it != set.End(); ++it)
			parent_set[*it] = set.GetHandle();
	}
	void AdaptiveMesh::Repartition()
	{
		if( !root.isValid() ) return;
		TagInteger redist = m->RedistributeTag();
		TagInteger new_owner = m->CreateTag("TEMPORARY_NEW_OWNER",DATA_INTEGER,ESET, NONE,1); //workaround for sets
		TagIntegerArray new_procs =  m->CreateTag("TEMPORARY_NEW_PROCESSORS",DATA_INTEGER,ESET, NONE);
		for(Storage::integer it = 0; it < m->EsetLastLocalID(); ++it) if( m->isValidElementSet(it) )
		{
			ElementSet set = m->EsetByLocalID(it);
			new_owner[set] = INT_MAX;
			new_procs[set].clear();
		}
		std::cout << __FUNCTION__ << ":" << __FILE__ << ":" << __LINE__ << std::endl;
		for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
		{
			Cell c = m->CellByLocalID(it);
			ElementSet parent(m,parent_set[c]);
			std::cout << "cell:" << it << " gid " << c.GlobalID() << " has parent " << parent.LocalID() << " " << parent.GetName() << " redist " << redist[c] << std::endl;
			new_owner[parent] = std::min(new_owner[parent],redist[c]);
			Storage::integer_array procs = new_procs[parent];
			Storage::integer_array::iterator insp = std::lower_bound(procs.begin(),procs.end(),redist[c]);
			if( insp == procs.end() || *insp != redist[c] ) procs.insert(insp,redist[c]);
		}
		SetNewOwner(root,new_owner);
		SetNewProcs(root,new_procs);
		m->ReduceData(new_owner,ESET,0,ReduceMin);
		m->ExchangeData(new_owner,ESET,0);
		m->ReduceData(new_procs,ESET,0,ReduceUnion);
		m->ExchangeData(new_procs,ESET,0);
		for(Storage::integer it = 0; it < m->EsetLastLocalID(); ++it) if( m->isValidElementSet(it) )
		{
			ElementSet set = m->EsetByLocalID(it);
			Storage::integer_array old_p = set->IntegerArray(m->ProcessorsTag());
			Storage::integer_array new_p = new_procs[set];
			Storage::integer_array sendto = set->IntegerArray(m->SendtoTag());
			sendto.resize(std::max(old_p.size(),new_p.size()));
			sendto.resize(std::set_difference(new_p.begin(),new_p.end(),old_p.begin(),old_p.end(),sendto.begin())-sendto.begin());
		}
		
		std::cout << __FUNCTION__ << ":" << __FILE__ << ":" << __LINE__ << std::endl;
		for(Storage::integer it = 0; it < m->EsetLastLocalID(); ++it) if( m->isValidElementSet(it) )
		{
			ElementSet set = m->EsetByLocalID(it);
			if( new_owner[set] == INT_MAX ) std::cout << "no new owner for set " << it << " " << set.GetName() << std::endl;
		}
		for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
		{
			Cell c = m->CellByLocalID(it);
			if( parent_set[c] == InvalidHandle() ) std::cout << "No parent set for cell:" << it << " gid " << c.GlobalID() << std::endl;
		}
		
		m->BeginModification();
		m->Redistribute();
		m->ExchangeData(hanging_nodes,CELL|FACE,0); //maybe will work
		m->ApplyModification(); //this will correctly remove links
		m->EndModification();
		m->ReorderEmpty(CELL|FACE|EDGE|NODE);
		RestoreParent(root);
		std::cout << __FUNCTION__ << ":" << __FILE__ << ":" << __LINE__ << std::endl;
		for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
		{
			Cell c = m->CellByLocalID(it);
			if( parent_set[c] == InvalidHandle() ) std::cout << "No parent set for cell:" << it << " gid " << c.GlobalID() << std::endl;
		}
	}
#endif
}
