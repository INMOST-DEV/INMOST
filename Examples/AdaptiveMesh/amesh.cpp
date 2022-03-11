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

__INLINE std::string NameSlash(std::string input)
{
	for(unsigned l = static_cast<unsigned>(input.size()); l > 0; --l)
		if( input[l-1] == '/' || input[l-1] == '\\' )
			return std::string(input.c_str() + l);
	return input;
}

const bool check_orientation = false;
const bool check_convexity = false;


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
	void ReduceSum(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		(void) size;
		element.Real(tag) += *((const INMOST_DATA_REAL_TYPE *)data);
	}
	
	void ReduceUnion(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		(void) size;
		const INMOST_DATA_INTEGER_TYPE * idata = (const INMOST_DATA_INTEGER_TYPE *)data;
		Storage::integer_array odata = element->IntegerArray(tag);
		std::vector<int> tmp(static_cast<size_t>(size)+static_cast<size_t>(odata.size()));
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
		if(tri_hanging_edges.isValid())
			tri_hanging_edges = m->DeleteTag(tri_hanging_edges);
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
		//m->ResolveSets();
	}
	
	AdaptiveMesh::AdaptiveMesh(Mesh & _m, bool skip_tri) : m(&_m), skip_tri(skip_tri)
	{
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//		model = NULL;
//#endif
		//create a tag that stores maximal refinement level of each element
		level = m->CreateTag("REFINEMENT_LEVEL",DATA_INTEGER,CELL|FACE|EDGE|NODE|ESET,NONE,1);
		//tag_status = m->CreateTag("TAG_STATUS",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
		set_id = m->CreateTag("SET_ID",DATA_INTEGER,CELL|ESET,ESET,1);
		//tag_an = m->CreateTag("TAG_AN",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
		//ref_tag = m->CreateTag("REF",DATA_REFERENCE,CELL|FACE|EDGE|NODE,NONE);
		//create a tag that stores links to all the hanging nodes of the cell
		hanging_nodes = m->CreateTag("HANGING_NODES",DATA_REFERENCE,CELL|FACE,NONE);
		if (skip_tri)
			tri_hanging_edges = m->CreateTag("TRI_HANGING_EDGES", DATA_REFERENCE, CELL, CELL);
		//create a tag that stores links to sets
		parent_set = m->CreateTag("PARENT_SET",DATA_REFERENCE,CELL,NONE,1);
	    size = m->GetProcessorsNumber();
    	rank = m->GetProcessorRank();
	}

	void ReportSub(ElementSet root, int tab, std::fstream & fout)
	{
		if( root.HaveChild() )
		{
			for(ElementSet jt = root.GetChild(); jt.isValid(); jt = jt.GetSibling())
				ReportSub(jt,tab+1,fout);
		}
		Storage::integer_array parr = root.IntegerArray(root.GetMeshLink()->ProcessorsTag());
		for(int q = 0; q < tab; ++q) fout << "  ";
		fout << "set " << root.GetName() << " procs ";
		for(unsigned q = 0; q < parr.size(); ++q) fout << parr[q] << " ";
		fout << " owner " << root.Integer(root.GetMeshLink()->OwnerTag());
		fout << " kids " << root.Size();
		fout << " status " << Element::StatusName(root.GetStatus()) << std::endl;
	}

	void AdaptiveMesh::ReportSets(std::fstream & fout)
	{
		for(Mesh::iteratorSet it = m->BeginSet(); it != m->EndSet(); ++it)
			if( !it->HaveParent() ) ReportSub(it->self(),0,fout);
	}
	
	AdaptiveMesh::~AdaptiveMesh()
	{
		//do not delete tags, user may want to repetitively use this class
		//as extension of class mesh in limited code span
	}
	
	void AdaptiveMesh::CheckParentSet(std::string file, int line)//, TagInteger indicator)
	{
		(void)file,(void)line;
		ENTER_FUNC();
#if !defined(NDEBUG)
		Storage::integer err = 0;
		Storage::integer_array procs;
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		{
			if( parent_set[*it] == InvalidHandle() )
			{
				REPORT_STR(m->GetProcessorRank() << " parent set not valid on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << " " << level[*it]);
				std::cout << m->GetProcessorRank() << " parent set not valid on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << " " << level[*it] << std::endl;
				err++;
			}
			else if( GetHandleElementType(parent_set[*it]) != ESET )
			{
				REPORT_STR(m->GetProcessorRank() << " parent set is something else " << ElementTypeName(GetHandleElementType(parent_set[*it])) << ":" << GetHandleID(parent_set[*it]) << " on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << " " << level[*it]);
				std::cout << m->GetProcessorRank() << " parent set is something else " << ElementTypeName(GetHandleElementType(parent_set[*it])) << ":" << GetHandleID(parent_set[*it]) << " on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << " " << level[*it] << std::endl;
				err++;
			}
			else if( parent_set[*it] != root.GetHandle() )//&& (!indicator.isValid() || indicator[*it]) )
			{
				ElementSet set(m,parent_set[*it]);
				if( !set.HaveParent() )
				{
					REPORT_STR(m->GetProcessorRank() << 
					" parent set " << ElementTypeName(GetHandleElementType(parent_set[*it])) << ":" << GetHandleID(parent_set[*it]) << 
					" name " << set.GetName() << " owner " << set.Integer(m->OwnerTag()) << " status " << Element::StatusName(set.GetStatus()) << " kids " << set.Size() <<
					" does not have parent " <<
					" on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << " owner " << it->Integer(m->OwnerTag()) << " lvl " << level[*it]);
					std::cout << m->GetProcessorRank() << 
					" parent set " << ElementTypeName(GetHandleElementType(parent_set[*it])) << ":" << GetHandleID(parent_set[*it]) << 
					" name " << set.GetName() << " owner " << set.Integer(m->OwnerTag()) << " status " << Element::StatusName(set.GetStatus()) << " kids " << set.Size() <<
					" does not have parent " <<
					" on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << " owner " << it->Integer(m->OwnerTag()) << " lvl " << level[*it] << std::endl;
					procs = it->IntegerArray(m->ProcessorsTag());
					std::cout << "cell procs:"; for(unsigned q = 0; q < procs.size(); ++q) std::cout << " " << procs[q]; std::cout << std::endl;
					procs = set.IntegerArray(m->ProcessorsTag());
					std::cout << "eset procs:"; for(unsigned q = 0; q < procs.size(); ++q) std::cout << " " << procs[q]; std::cout << std::endl;
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
						std::cout << m->GetProcessorRank() << 
						" parent set " << ElementTypeName(GetHandleElementType(parent_set[*it])) << ":" << GetHandleID(parent_set[*it]) << 
						" name " << set.GetName() << " owner " << set.Integer(m->OwnerTag()) << " status " << Element::StatusName(set.GetStatus()) << 
						" has parent " << ElementTypeName(GetHandleElementType(parent)) << ":" << GetHandleID(parent) <<
						" on CELL:" << it->LocalID() << " " << Element::StatusName(it->GetStatus()) << " " << level[*it] << std::endl;
						err++;
					}
				}
			}
		}
		err = m->Integrate(err);
		if( err ) 
		{
			REPORT_STR(rank << " error in " << __FUNCTION__ << " " << file << ":" << line);
			std::cout << rank << " error in " << __FUNCTION__ << " " << file << ":" << line << std::endl;
			exit(-1);
		}
#endif //NDEBUG
		EXIT_FUNC();
	}
	
	bool AdaptiveMesh::Refine(TagInteger indicator)
	{
		static int fi = 0;
        ENTER_FUNC();
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//		if (model) model->PrepareAdaptation(*m);
//#endif
		if (!indicator.isValid() || !indicator.isDefined(CELL))
		{
			indicator = m->CreateTag(indicator.isValid() ? indicator.GetTagName() : "indicator", DATA_INTEGER, CELL, NONE, 1);
			for (Storage::integer i = 0; i < m->CellLastLocalID(); ++i) if (m->isValidCell(i))
				indicator[m->CellByLocalID(i)] = 0;
		}
		//TagInteger indicator = m->CreateTag(indicator.GetTagName(), DATA_INTEGER, CELL, NONE, 1);
		for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
			(*it)->RefineIndicator(*this,indicator);
		int refine = 0, tot = m->TotalNumberOf(CELL);
		for (Storage::integer i = 0; i < m->CellLastLocalID(); ++i) if (m->isValidCell(i))
		{
			Cell c = m->CellByLocalID(i);
			if (indicator[c]) refine++;
		}
		refine = m->Integrate(refine);
		if (m->GetProcessorRank() == 0)
			std::cout << __FUNCTION__ << " indicator marked " << refine << "/" << tot << std::endl;
		if (!refine) return false;
		m->ExchangeData(indicator, CELL);
		for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
			(*it)->BeginRefinement();
		static int call_counter = 0;
		Storage::integer ret = 0; //return number of refined cells
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
		m->CheckSetLinks(__FILE__,__LINE__);
		CheckParentSet(__FILE__,__LINE__);
		EXIT_BLOCK();
		/*
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
		//CheckParentSet(Tag());
		EXIT_BLOCK();
		*/
		//ENTER_BLOCK();
		//for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it)
		//	if( it->getNodes().size() != 2 ) {REPORT_STR("edge " << it->LocalID() << " has " << it->getNodes().size() << " nodes ");}
		//EXIT_BLOCK();

		//m->Save("before_refine_parent"+std::to_string(fi)+".pvtk");
		//std::cout << "Save before_refine_parent"+std::to_string(fi)+".pvtk" << std::endl;
		

		
		int schedule_counter = 1; //indicates order in which refinement will be scheduled
		Storage::integer scheduled = 1; //indicates that at least one element was scheduled on current sweep
		ENTER_BLOCK();
		//0. Extend indicator for edges and faces
		indicator = m->CreateTag(indicator.GetTagName(),DATA_INTEGER,FACE|EDGE,NONE,1);
		//1.Communicate indicator - it may be not synced
		m->ExchangeData(indicator,CELL,0);
		while(scheduled)
		{
			REPORT_VAL("scheduled",scheduled);
			//2.Propogate indicator down to the faces,edges
			//  select schedule for them
			ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( indicator[c] == schedule_counter )
				{
					ElementArray<Element> adj = c.getAdjElements(FACE|EDGE);
					for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
					{
						if (level[adj[kt]] == level[c]) //do not schedule finer or coarser elements
						{
#if defined(USE_OMP)
#pragma omp critical
#endif
							indicator[adj[kt]] = schedule_counter; //refine together with current cell
						}
					}
				}
			}
			EXIT_BLOCK();
			//3.Communicate indicator on faces and edges
			m->ReduceData(indicator,FACE|EDGE,0,ReduceMax);
			m->ExchangeData(indicator,FACE|EDGE,0);
			//4.Check for each cell if there is
			//  any hanging node with adjacent in a need to refine,
			//  schedule for refinement earlier.
			scheduled = 0;
			ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:scheduled)
#endif
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
			//5.Exchange indicator on cells
			m->ReduceData(indicator,CELL,0,ReduceMax);
			m->ExchangeData(indicator,CELL,0);
			//6.Go back to 1 until no new elements scheduled
			scheduled = m->Integrate(scheduled);
			if( scheduled ) schedule_counter++;
		}
		//m->ExchangeData(indicator,CELL | FACE | EDGE,0);
		EXIT_BLOCK();
		
		//m->ExchangeData(hanging_nodes,CELL|FACE,NONE);
		//m->Save("indicator.pmf");
		
		//ENTER_BLOCK();
		//for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it)
		//	if( it->getNodes().size() != 2 ) {REPORT_STR("edge " << it->LocalID() << " has " << it->getNodes().size() << " nodes ");}
		//EXIT_BLOCK();
		
		m->ExchangeData(hanging_nodes,CELL | FACE,0);
		if (tri_hanging_edges.isValid())
			m->ExchangeData(tri_hanging_edges, CELL, 0);
		
		//if( !Element::CheckConnectivity(m) ) std::cout << __FILE__ << ":" << __LINE__ << " broken connectivity" << std::endl;
		//m->CheckCentroids(__FILE__,__LINE__);
		//6.Refine

		assert(Element::CheckConnectivity(m));
		CheckClosure(__FILE__,__LINE__);
		ENTER_BLOCK();
		m->BeginModification();
		while(schedule_counter)
		{
			assert(Element::CheckConnectivity(m));
			CheckClosure(__FILE__,__LINE__);
			Storage::real xyz[3] = {0,0,0};
			//7.split all edges of the current schedule
			ENTER_BLOCK();
			{
				int new_nodes = 0, splits = 0;
				double t1, t2, tadj = 0, tcreate = 0, thanging = 0, tdata = 0, tsplit = 0;
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
						t1 = Timer();
						ElementArray<Face> edge_faces = e.getFaces();
						t2 = Timer(), tadj += t2 - t1, t1 = t2;
						//location on the center of the edge
						for(Storage::integer d = 0; d < m->GetDimensions(); ++d)
							xyz[d] = (e.getBeg().Coords()[d]+e.getEnd().Coords()[d])*0.5;
						//todo: request transformation of node location according to geometrical model
						//create middle node
						Node n = m->CreateNode(xyz);
						new_nodes++;
						//set increased level for new node
						level[n] = level[e.getBeg()] = level[e.getEnd()] = level[e]+1;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//						if (model) model->NewNode(e, n);
//#endif
						for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
							(*it)->NewNode(e,n);
						t2 = Timer(), tcreate += t2 - t1, t1 = t2;
						//for each face provide link to a new hanging node
						for(ElementArray<Face>::size_type kt = 0; kt < edge_faces.size(); ++kt)
							hanging_nodes[edge_faces[kt]].push_back(n);
						t2 = Timer(), thanging += t2 - t1, t1 = t2;
						//CheckClosure(__FILE__,__LINE__);
						//split the edge by the middle node
						ElementArray<Edge> new_edges = Edge::SplitEdge(e,ElementArray<Node>(m,1,n.GetHandle()),0);
						splits++;
						for(ElementArray<Face>::size_type kt = 0; kt < edge_faces.size(); ++kt) assert(edge_faces[kt].Closure());
						t2 = Timer(), tsplit += t2 - t1, t1 = t2;
						//set increased level for new edges
						level[new_edges[0]] = level[new_edges[1]] = level[e]+1;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//						if( model ) model->EdgeRefinement(e,new_edges);
//#endif
						for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
							(*it)->EdgeRefinement(e, new_edges);
						t2 = Timer(), tdata += t2 - t1, t1 = t2;
						//for(int q = 0; q < 2; ++q)
						//{
						//	REPORT_STR("new edges["<<q<<"]" << new_edges[q].LocalID() << " nodes " << new_edges[q].getBeg().LocalID() << "," << new_edges[q].getEnd().LocalID() << " level " << level[new_edges[q]]);
						//}
						//CheckClosure(__FILE__,__LINE__);
						//if( !Element::CheckConnectivity(m) ) std::cout << __FILE__ << ":" << __LINE__ << " broken connectivity" << std::endl;
					}
				}
				REPORT_VAL("adjacencies", tadj);
				REPORT_VAL("create (new nodes)", tcreate);
				REPORT_VAL("hanging", thanging);
				REPORT_VAL("data", tdata);
				REPORT_VAL("split", tsplit);
				REPORT_VAL("new nodes", new_nodes);
				REPORT_VAL("splits", splits);
			}
			EXIT_BLOCK();
			
			ENTER_BLOCK();
			assert(Element::CheckConnectivity(m));
			CheckClosure(__FILE__,__LINE__);
			EXIT_BLOCK();
			//8.split all faces of the current schedule, using hanging nodes on edges
			ENTER_BLOCK();
			{
				int new_nodes = 0, new_edges = 0, splits = 0;
				double t1, t2, tadj = 0, tcreate = 0, thanging = 0, tdata = 0, tsplit = 0;
				for(Storage::integer it = 0; it < m->FaceLastLocalID(); ++it) if( m->isValidFace(it) )
				{
					Face f = m->FaceByLocalID(it);
					if( !f.Hidden() && indicator[f] == schedule_counter )
					{
#if !defined(NDEBUG)
						ElementArray<Edge> face_edges = f.getEdges();
						for(ElementArray<Edge>::iterator jt = face_edges.begin(); jt != face_edges.end(); ++jt)
						{
							if( level[*jt] != level[f]+1 )
							{
								std::cout << m->GetProcessorRank() << " face " << f.LocalID();
								std::cout << " " << Element::StatusName(f.GetStatus()) << " owner " << f.Integer(m->OwnerTag());
								std::cout << " lvl " << level[f] << " ind " << indicator[f];
								std::cout << " edge " << jt->LocalID();
								std::cout << " " << Element::StatusName(jt->GetStatus()) << " owner " << jt->Integer(m->OwnerTag());
								std::cout << " lvl " << level[*jt] << " ind " << indicator[*jt];
								std::cout << std::endl;
								ElementArray<Cell> adj_cells;
								adj_cells = f.getCells();
								std::cout << "face cells ";
								for(ElementArray<Cell>::iterator kt = adj_cells.begin(); kt != adj_cells.end(); ++kt)
								{
									std::cout << m->GetProcessorRank() << " face " << kt->LocalID();
									std::cout << " " << Element::StatusName(kt->GetStatus()) << " owner " << kt->Integer(m->OwnerTag());
									std::cout << " lvl " << level[*kt] << " ind " << indicator[*kt];
									std::cout << std::endl;
								}
								std::cout << "edge cells ";
								adj_cells = jt->getCells();
								for(ElementArray<Cell>::iterator kt = adj_cells.begin(); kt != adj_cells.end(); ++kt)
								{
									std::cout << m->GetProcessorRank() << " face " << kt->LocalID();
									std::cout << " " << Element::StatusName(kt->GetStatus()) << " owner " << kt->Integer(m->OwnerTag());
									std::cout << " lvl " << level[*kt] << " ind " << indicator[*kt];
									std::cout << std::endl;
								}
							}
							assert(level[*jt] == level[f]+1);
						}
#endif //NDEBUG
						t1 = Timer();
						//remember adjacent cells that should get information about new hanging node
						//and new hanging edges
						ElementArray<Cell> face_cells = f.getCells();
						t2 = Timer(), tadj += t2 - t1, t1 = t2;
						//connect face center to hanging nodes of the face
						Storage::reference_array face_hanging_nodes = hanging_nodes[f];
						ElementArray<Node> edge_nodes(m, 2); //to create new edges
						ElementArray<Edge> hanging_edges(m);
						//skip triangle
						if (skip_tri && f.GetGeometricType() == Element::Tri)
						{
							hanging_edges.resize(3);
							for (Storage::reference_array::size_type kt = 0; kt < face_hanging_nodes.size(); ++kt)
							{
								edge_nodes[0] = face_hanging_nodes[kt].getAsNode();
								edge_nodes[1] = face_hanging_nodes[(kt+1)% face_hanging_nodes.size()].getAsNode();
								hanging_edges[kt] = m->CreateEdge(edge_nodes).first;
								level[hanging_edges[kt]] = level[f] + 1;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//								if (model) model->NewEdge(f, hanging_edges[kt]);
//#endif
								for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
									(*it)->NewEdge(f, hanging_edges[kt]);
								for (ElementArray<Face>::size_type kt = 0; kt < face_cells.size(); ++kt)
									tri_hanging_edges[face_cells[kt]].push_back(hanging_edges[kt]);
							}
							t2 = Timer(), thanging += t2 - t1, t1 = t2;
						}
						else
						{
							//create node at face center
							//f->Centroid(xyz);
							for (int d = 0; d < 3; ++d) xyz[d] = 0.0;
							for (Storage::reference_array::size_type kt = 0; kt < face_hanging_nodes.size(); ++kt)
								for (int d = 0; d < 3; ++d) xyz[d] += face_hanging_nodes[kt].getAsNode().Coords()[d];
							for (int d = 0; d < 3; ++d) xyz[d] /= (Storage::real)face_hanging_nodes.size();
							//todo: request transformation of node location according to geometrical model
							//create middle node
							Node n = m->CreateNode(xyz);
							new_nodes++;
							//set increased level for the new node
							level[n] = level[f] + 1;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//							if (model) model->NewNode(f, n, face_hanging_nodes);
//#endif
							for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
								(*it)->NewNode(f, n, face_hanging_nodes);
							t2 = Timer(), tcreate += t2 - t1, t1 = t2;
							//for each cell provide link to new hanging node
							for (ElementArray<Face>::size_type kt = 0; kt < face_cells.size(); ++kt)
								hanging_nodes[face_cells[kt]].push_back(n);
							edge_nodes[0] = n;
							hanging_edges.resize(face_hanging_nodes.size());
							for (Storage::reference_array::size_type kt = 0; kt < face_hanging_nodes.size(); ++kt)
							{
								edge_nodes[1] = face_hanging_nodes[kt].getAsNode();
								hanging_edges[kt] = m->CreateEdge(edge_nodes).first;
								new_edges++;
								//set increased level for new edges
								level[hanging_edges[kt]] = level[f] + 1;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//								if (model) model->NewEdge(f, hanging_edges[kt]);
//#endif
								for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
									(*it)->NewEdge(f, hanging_edges[kt]);
							}
							t2 = Timer(), thanging += t2 - t1, t1 = t2;
						}
						//split the face by these edges
						ElementArray<Face> new_faces = Face::SplitFace(f,hanging_edges,0);
						splits++;
						t2 = Timer(), tsplit += t2 - t1, t1 = t2;
						//set increased level to new faces
						for(ElementArray<Face>::size_type kt = 0; kt < new_faces.size(); ++kt)
							level[new_faces[kt]] = level[f]+1;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//						if( model ) model->FaceRefinement(f,new_faces);
//#endif
						for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
							(*it)->FaceRefinement(f, new_faces);
						t2 = Timer(), tdata += t2 - t1, t1 = t2;
					}
				}
				REPORT_VAL("adjacencies", tadj);
				REPORT_VAL("create (new nodes)", tcreate);
				REPORT_VAL("hanging (new edges)", thanging);
				REPORT_VAL("data", tdata);
				REPORT_VAL("split", tsplit);
				REPORT_VAL("new nodes", new_nodes);
				REPORT_VAL("new edges", new_edges);
				REPORT_VAL("splits", splits);
			}
			EXIT_BLOCK();
			
			ENTER_BLOCK();
			assert(Element::CheckConnectivity(m));
			CheckClosure(__FILE__,__LINE__);
			EXIT_BLOCK();
			
			TagReferenceArray internal_face_edges;
			MarkerType mark_cell_edges, mark_hanging_nodes;
			
			ENTER_BLOCK();
			//this tag helps recreate internal face
			internal_face_edges = m->CreateTag("INTERNAL_FACE_EDGES",DATA_REFERENCE,NODE,NODE,4);
			//this marker helps detect edges of current cell only
			mark_cell_edges = m->CreateMarker();
			//this marker helps detect nodes hanging on edges of unrefined cell
			mark_hanging_nodes = m->CreateMarker();
			EXIT_BLOCK();
			//9.split all cells of the current schedule
			ENTER_BLOCK();
			{
				int new_nodes = 0, new_edges = 0, new_faces = 0, new_sets = 0, splits = 0;
				double t1, t2, tadj = 0, tcreate = 0, thanging1 = 0, thanging2 = 0, tset = 0, tdata = 0, tsplit = 0, tsort = 0, tpconn = 0;
				for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
				{
					Cell c = m->CellByLocalID(it);
					if( !c.Hidden() && indicator[c] == schedule_counter )
					{
						t1 = Timer();
						Node n;
						Storage::reference_array cell_hanging_nodes = hanging_nodes[c]; //nodes to be connected
						if (!skip_tri || cell_hanging_nodes.size() > 3)
						{
							//create node at cell center
							for (int d = 0; d < 3; ++d) xyz[d] = 0.0;
							for (Storage::reference_array::size_type kt = 0; kt < cell_hanging_nodes.size(); ++kt)
								for (int d = 0; d < 3; ++d) xyz[d] += cell_hanging_nodes[kt].getAsNode().Coords()[d];
							for (int d = 0; d < 3; ++d) xyz[d] /= (Storage::real)cell_hanging_nodes.size();
							//c->Centroid(xyz);
							//todo: request transformation of node location according to geometrical model
							//create middle node
							n = m->CreateNode(xyz);
							new_nodes++;
							//set increased level for the new node
							level[n] = level[c] + 1;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//							if (model) model->NewNode(c, n, cell_hanging_nodes);
//#endif
							for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
								(*it)->NewNode(c, n, cell_hanging_nodes);
						}
						t2 = Timer(), tcreate += t2 - t1, t1 = t2;
						ElementArray<Face> internal_faces(m);
						if (n.isValid())
						{
							//retrive all edges of current face to mark them
							ElementArray<Edge> cell_edges = c.getEdges();
#if !defined(NDEBUG)
							for (ElementArray<Edge>::iterator jt = cell_edges.begin(); jt != cell_edges.end(); ++jt) assert(level[*jt] == level[c] + 1);
							ElementArray<Face> cell_faces = c.getFaces();
							for (ElementArray<Face>::iterator jt = cell_faces.begin(); jt != cell_faces.end(); ++jt) assert(level[*jt] == level[c] + 1);
#endif //NDEBUG
							t2 = Timer(), tadj += t2 - t1, t1 = t2;
							//mark all edges so that we can retive them later
							cell_edges.SetMarker(mark_cell_edges);
							//connect face center to centers of faces by edges
							ElementArray<Node> edge_nodes(m, 2);
							ElementArray<Edge> edges_to_faces(m, cell_hanging_nodes.size());
							edge_nodes[0] = n;
							for (Storage::reference_array::size_type kt = 0; kt < cell_hanging_nodes.size(); ++kt)
							{
								assert(cell_hanging_nodes[kt].isValid());
								//todo: unmark hanging node on edge if no more cells depend on it
								edge_nodes[1] = cell_hanging_nodes[kt].getAsNode();
								edges_to_faces[kt] = m->CreateEdge(edge_nodes).first;
								new_edges++;
								//set increased level for new edges
								level[edges_to_faces[kt]] = level[c] + 1;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//								if (model) model->NewEdge(c, edges_to_faces[kt]);
//#endif
								for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
									(*it)->NewEdge(c, edges_to_faces[kt]);
								//for each node other then the hanging node of the face
								//(this is hanging node on the edge)
								//we record a pair of edges to reconstruct internal faces
								ElementArray<Edge> hanging_edges = cell_hanging_nodes[kt].getEdges(mark_cell_edges, 0);
								for (ElementArray<Edge>::size_type lt = 0; lt < hanging_edges.size(); ++lt)
								{
									//get hanging node on the edge
									assert(hanging_edges[lt].getBeg() == cell_hanging_nodes[kt] || hanging_edges[lt].getEnd() == cell_hanging_nodes[kt]);
									Node v = hanging_edges[lt].getBeg() == cell_hanging_nodes[kt] ? hanging_edges[lt].getEnd() : hanging_edges[lt].getBeg();
									//mark so that we can collect all of them
									v.SetMarker(mark_hanging_nodes);
									//fill the edges
									Storage::reference_array face_edges = internal_face_edges[v];
									//fill first two in forward order
									//this way we make a closed loop
									assert(face_edges[0] == InvalidElement() || face_edges[2] == InvalidElement());
									if (face_edges[0] == InvalidElement())
									{
										face_edges[0] = edges_to_faces[kt];
										face_edges[1] = hanging_edges[lt];
									}
									else //last two in reverse
									{
										assert(face_edges[2] == InvalidElement());
										face_edges[2] = hanging_edges[lt];
										face_edges[3] = edges_to_faces[kt];
									}
								}
							}
							//remove marker from cell edges
							cell_edges.RemMarker(mark_cell_edges);
							t2 = Timer(), thanging1 += t2 - t1, t1 = t2;
							//now we have to create internal faces
							ElementArray<Node> edge_hanging_nodes = c.getNodes(mark_hanging_nodes, 0);
							internal_faces.resize(edge_hanging_nodes.size());
							//unmark hanging nodes on edges
							edge_hanging_nodes.RemMarker(mark_hanging_nodes);
							for (ElementArray<Node>::size_type kt = 0; kt < edge_hanging_nodes.size(); ++kt)
							{
								//create a face based on collected edges
								Storage::reference_array face_edges = internal_face_edges[edge_hanging_nodes[kt]];
								int valid = 0;
								for (int q = 0; q < 4; ++q)
									valid += face_edges[q].isValid() ? 1 : 0;
								assert(valid == 2 || valid == 4);
								if (valid == 4) //create a quadrilateral
									internal_faces[kt] = m->CreateFace(ElementArray<Edge>(m, face_edges.begin(), face_edges.end())).first;
								else if (valid == 2) //create a triangle (the neighbouring face was a triangle)
								{
									edge_nodes[1] = edge_hanging_nodes[kt].getAsNode();
									ElementArray<Edge> tri_edges(m, 3);
									tri_edges[0] = face_edges[0].getAsEdge();
									tri_edges[1] = face_edges[1].getAsEdge();
									tri_edges[2] = m->CreateEdge(edge_nodes).first;
									//double check, is this really a new edge???
									level[tri_edges[2]] = level[c] + 1;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//									if (model) model->NewEdge(c, tri_edges[2]);
//#endif
									for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
										(*it)->NewEdge(c, tri_edges[2]);
									internal_faces[kt] = m->CreateFace(tri_edges).first;
								}
								new_faces++;
								//set increased level
								level[internal_faces[kt]] = level[c] + 1;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//								if (model) model->NewFace(c, internal_faces[kt]);
//#endif
								for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
									(*it)->NewFace(c, internal_faces[kt]);
								//clean up structure, so that other cells can use it
								edge_hanging_nodes[kt].DelData(internal_face_edges);
							}
							if (skip_tri)
							{
								ElementArray<Node> tri_nodes(m, 3);
								Storage::reference_array cell_tri_hanging_edges = tri_hanging_edges[c]; //edges to be connected
								tri_nodes[0] = n;
								for (Storage::reference_array::size_type kt = 0; kt < cell_tri_hanging_edges.size(); ++kt)
								{
									tri_nodes[1] = cell_tri_hanging_edges[kt].getAsEdge().getBeg();
									tri_nodes[2] = cell_tri_hanging_edges[kt].getAsEdge().getEnd();
									internal_faces.push_back(m->CreateFace(tri_nodes).first);
									level[internal_faces.back()] = level[c] + 1;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//									if (model) model->NewFace(c, internal_faces.back());
//#endif
									for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
										(*it)->NewFace(c, internal_faces.back());
								}
							}
							t2 = Timer(), thanging2 += t2 - t1, t1 = t2;
						}
						//if( c.GlobalID() == 228 )
						//{
						//	double cnt[3];
						//	c.Centroid(cnt);
						//	std::cout << "Split CELL:" << c.LocalID() << " " << c.GlobalID() << " " << Element::StatusName(c.GetStatus()) << " " << cnt[0] << " " << cnt[1] << " " << cnt[2] << std::endl;
						//	
						//}
						//split the cell
						//retrive parent set
						ElementSet parent(m,parent_set[c]);
						//create set corresponding to old coarse cell
						//Storage::real cnt[3];
						//c.Centroid(cnt);
						std::stringstream set_name;
						//set_name << parent.GetName() << "_C" << c.GlobalID(); //rand may be unsafe
						if( parent == root )
							set_name << "AM_R" << set_id[c];
						else
							set_name << parent.GetName() << "C" << set_id[c];
						//set_name << base64_encode_((unsigned char *)cnt,3*sizeof(double)/sizeof(unsigned char));
#if !defined(NDEBUG)
						ElementSet check_set = m->GetSet(set_name.str());
						if( check_set.isValid() )
						{
							std::cout << rank << " set " << set_name.str() << " for cell " << c.GlobalID() << " " << Element::StatusName(c.GetStatus()) << " already exists" << std::endl;
							if( check_set->HaveParent() )
								std::cout << rank << " parent is " << check_set->GetParent()->GetName() << " cell parent is " << parent.GetName() << std::endl;
							std::cout << rank << " Elements of " << check_set.GetName() << ": ";
							for(ElementSet::iterator it = check_set.Begin(); it != check_set.End(); ++it)
								std::cout << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << "," << it->GlobalID() << "," << Element::StatusName(c.GetStatus()) << "," << level[*it] << "," << indicator[*it] << " ";
							std::cout << std::endl;
							std::cout << rank << " Elements of " << parent.GetName() << ": ";
							for(ElementSet::iterator it = parent.Begin(); it != parent.End(); ++it)
								std::cout << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << "," << it->GlobalID() << "," << Element::StatusName(c.GetStatus()) << "," << level[*it] << "," << indicator[*it] << " ";
							std::cout << std::endl;
							if( parent.HaveChild() )
							{
								std::cout << rank << " Children of " << parent.GetName() << ": ";
								for(ElementSet jt = parent.GetChild(); jt.isValid(); jt = jt.GetSibling() )
									std::cout << jt.GetName() << " size " << jt.Size() << " ";
								std::cout << std::endl;
							}
							exit(-1);
						}
#endif
						ElementSet cell_set = m->CreateSetUnique(set_name.str()).first;
						new_sets++;
						//cell_set->SetExchange(ElementSet::SYNC_ELEMENTS_ALL);
						level[cell_set] = level[c]+1;
						set_id[cell_set] = set_id[c];
						
						t2 = Timer(), tset += t2 - t1, t1 = t2;

						ElementArray<Cell> new_cells = Cell::SplitCell(c,internal_faces,0);
						splits++;
						t2 = Timer(), tsplit += t2 - t1, t1 = t2;
						
						std::sort(new_cells.begin(),new_cells.end(),Mesh::CentroidComparator(m));
						
						t2 = Timer(), tsort += t2 - t1, t1 = t2;
						//set up increased level for the new cells
						for(ElementArray<Cell>::size_type kt = 0; kt < new_cells.size(); ++kt)
						{
							set_id[new_cells[kt]] = (int)kt;
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
						t2 = Timer(), tpconn += t2 - t1, t1 = t2;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//						if (model) model->CellRefinement(c, new_cells, cell_set);
//#endif
						for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
							(*it)->CellRefinement(c, new_cells, cell_set);
						t2 = Timer(), tdata += t2 - t1, t1 = t2;
						//else assert(cell_set->GetParent() == parent);
						//increment number of refined cells
						ret++;
					}
				}
				REPORT_VAL("adjacencies", tadj);
				REPORT_VAL("create (new nodes)", tcreate);
				REPORT_VAL("hanging1 (new edges)", thanging1);
				REPORT_VAL("hanging2 (new faces)", thanging2);
				REPORT_VAL("create set", tset);
				REPORT_VAL("connect set", tpconn);
				REPORT_VAL("data", tdata);
				REPORT_VAL("split", tsplit);
				REPORT_VAL("sort cells", tsort);
				REPORT_VAL("new nodes", new_nodes);
				REPORT_VAL("new edges", new_edges);
				REPORT_VAL("new faces", new_faces);
				REPORT_VAL("new sets", new_sets);
				REPORT_VAL("splits", splits);
			}
			EXIT_BLOCK();
			ENTER_BLOCK();
			m->ReleaseMarker(mark_hanging_nodes);
			m->ReleaseMarker(mark_cell_edges);
			m->DeleteTag(internal_face_edges);
			EXIT_BLOCK();
			ENTER_BLOCK();
			assert(Element::CheckConnectivity(m));
			CheckClosure(__FILE__,__LINE__);
			EXIT_BLOCK();
			//10.jump to later schedule, and go to 7.
			REPORT_VAL("schedule counter",schedule_counter);
			schedule_counter--;
		}
		m->CheckSetLinks(__FILE__,__LINE__);
		//free created tag
		m->DeleteTag(indicator,FACE|EDGE);

		if( check_orientation )
		{
			ENTER_BLOCK();
			int nfixed = 0, nfixednew = 0, nfixedbnd = 0, nghost = 0;
			for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
				if( !it->CheckNormalOrientation() )
				{
					//it->FixNormalOrientation();
					nfixed++;
					if( it->New() )
						nfixednew++;
					if( it->Boundary() )
						nfixedbnd++;
					if( it->GetStatus() == Element::Ghost )
						nghost++;
					std::cout << "rank " << m->GetProcessorRank() << " face " << it->GetHandle() << std::endl;
				}
				//std::cout << "Face " << it->LocalID() << " oriented incorrectly " << std::endl;
			if( nfixed ) 
			{
				std::cout << __FILE__ << ":" << __LINE__ << " rank " << rank << " bad " << nfixed << " (new " << nfixednew << " bnd " << nfixedbnd << " ghost " << nghost << ") faces " << std::endl;
				REPORT_STR(rank << " bad " << nfixed << " faces");
			}
			EXIT_BLOCK();
		}
		
		if( check_convexity )
		{
			ENTER_BLOCK();
			int nbad = 0;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				if( !it->CheckConvexity() ) nbad++;
			if( nbad ) std::cout << __FILE__ << ":" << __LINE__ << " rank " << rank << " nonconvex cells: " << nbad << std::endl;
			EXIT_BLOCK();
		}
		ENTER_BLOCK();
		m->CheckSetLinks(__FILE__,__LINE__);
		EXIT_BLOCK();
		//ENTER_BLOCK();
		//m->Barrier();
		//EXIT_BLOCK();
		//11. Restore parallel connectivity, global ids
		m->ResolveModification();
		
		if( check_orientation )
		{
			ENTER_BLOCK();
			int nfixed = 0, nfixednew = 0, nfixedbnd = 0, nghost = 0;
			for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
				if( !it->CheckNormalOrientation() )
				{
					//it->FixNormalOrientation();
					nfixed++;
					if( it->New() )
						nfixednew++;
					if( it->Boundary() )
						nfixedbnd++;
					if( it->GetStatus() == Element::Ghost )
						nghost++;
					std::cout << "rank " << m->GetProcessorRank() << " face " << it->GetHandle() << std::endl;
				}
				//std::cout << "Face " << it->LocalID() << " oriented incorrectly " << std::endl;
			if( nfixed ) 
			{
				std::cout << __FILE__ << ":" << __LINE__ << " rank " << rank << " bad " << nfixed << " (new " << nfixednew << " bnd " << nfixedbnd << " ghost " << nghost << ") faces " << std::endl;
				REPORT_STR(rank << " bad " << nfixed << " faces");
			}
			EXIT_BLOCK();
		}
		
		if( check_convexity )
		{
			ENTER_BLOCK();
			int nbad = 0;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				if( !it->CheckConvexity() ) nbad++;
			if( nbad ) std::cout << __FILE__ << ":" << __LINE__ << " rank " << rank << " nonconvex cells: " << nbad << std::endl;
			EXIT_BLOCK();
		}
		//m->SynchronizeMarker(m->NewMarker(),CELL|FACE|EDGE|NODE,SYNC_BIT_OR);
        //ExchangeGhost(3,NODE); // Construct Ghost cells in 2 layers connected via nodes
		//12. Let the user update their data
		//todo: call back function for New( cells
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//		if( model ) model->Adaptation(*m);
//#endif
		ENTER_BLOCK();
		for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
			(*it)->Refinement();
		EXIT_BLOCK();
		m->CheckSetLinks(__FILE__,__LINE__);
		//13. Delete old elements of the mesh
		m->ApplyModification();


		if( check_orientation )
		{
			ENTER_BLOCK();
			int nfixed = 0, nfixednew = 0, nfixedbnd = 0, nghost = 0;
			for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
				if( !it->CheckNormalOrientation() )
				{
					//it->FixNormalOrientation();
					nfixed++;
					if( it->New() )
						nfixednew++;
					if( it->Boundary() )
						nfixedbnd++;
					if( it->GetStatus() == Element::Ghost )
						nghost++;
					std::cout << "rank " << m->GetProcessorRank() << " face " << it->GetHandle() << " cells " << it->nbAdjElements(CELL) << std::endl;
				}
				//std::cout << "Face " << it->LocalID() << " oriented incorrectly " << std::endl;
			if( nfixed ) 
			{
				std::cout << __FILE__ << ":" << __LINE__ << " rank " << rank << " bad " << nfixed << " (new " << nfixednew << " bnd " << nfixedbnd << " ghost " << nghost << ") faces " << std::endl;
				REPORT_STR(rank << " bad " << nfixed << " faces");
			}
			EXIT_BLOCK();
		}
		
		if( check_convexity )
		{
			ENTER_BLOCK();
			int nbad = 0;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				if( !it->CheckConvexity() ) nbad++;
			if( nbad ) std::cout << __FILE__ << ":" << __LINE__ << " rank " << rank << " nonconvex cells: " << nbad << std::endl;
			EXIT_BLOCK();
		}
		
		//m->ExchangeGhost(1,NODE,m->NewMarker());
		//14. Done
        //cout << rank << ": Before end " << std::endl;
		m->EndModification();
		EXIT_BLOCK();
		assert(Element::CheckConnectivity(m));
		CheckClosure(__FILE__,__LINE__);
		
		//keep links to prevent loss during balancing
		m->ExchangeData(parent_set,CELL,0);
		m->ExchangeData(hanging_nodes,CELL | FACE,0);
		if( tri_hanging_edges.isValid())
			m->ExchangeData(tri_hanging_edges, CELL, 0);
		/*
		ENTER_BLOCK();
		for(Storage::integer it = 0; it < m->EsetLastLocalID(); ++it) if( m->isValidElementSet(it) )
		{
			ElementSet set = m->EsetByLocalID(it);
			if( set.GetName().substr(0,3) == "AM_" )
			{
				//int imax = -1, imin = INT_MAX;
				//for(ElementSet::iterator jt = set.Begin(); jt != set.End(); ++jt)
				//{
				//	imax = std::max(imax,indicator[*jt]);
				//	imin = std::min(imin,indicator[*jt]);
				//}
				//std::cout << "on proc " << m->GetProcessorRank() << " set " << set.GetName() << " size " << set.Size() << " set indicator " << indicator[set] << " elements indicator " << imin << ":" << imax;
				//if( set.HaveParent() ) std::cout << " parent " << set.GetParent().GetName();
				//std::cout << std::endl;
				//if( !set.HaveChild() )
					set.SynchronizeSetParents();
			}
		}
		m->ExchangeMarked();
		EXIT_BLOCK();
		*/
		CheckParentSet(__FILE__,__LINE__);
		
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
		
		if( check_orientation )
		{
			ENTER_BLOCK();
			int nfixed = 0, nbnd = 0, nghost = 0;
			for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
				if( !it->CheckNormalOrientation() )
				{
					//it->FixNormalOrientation();
					nfixed++;
					if( it->Boundary() ) nbnd++;
					if( it->GetStatus() == Element::Ghost ) nghost++;
				}
				//std::cout << "Face " << it->LocalID() << " oriented incorrectly " << std::endl;
			if( nfixed ) 
			{
				std::cout << __FILE__ << ":" << __LINE__ << " rank " << rank << " bad " << nfixed << "(bnd " << nbnd << " ghost " << nghost << ") faces " << std::endl;
				REPORT_STR(rank << " bad " << nfixed << " faces");
			}
			EXIT_BLOCK();
		}
		
		if( check_convexity )
		{
			ENTER_BLOCK();
			int nbad = 0;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				if( !it->CheckConvexity() ) nbad++;
			if( nbad ) std::cout << __FILE__ << ":" << __LINE__ << " rank " << rank << " nonconvex cells: " << nbad << std::endl;
			EXIT_BLOCK();
		}
		
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
		ENTER_BLOCK();
		for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
			(*it)->EndRefinement();
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
	

	bool AdaptiveMesh::Coarse(TagInteger indicator)
	{
		std::string file;
		ENTER_FUNC();
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//		if (model) model->PrepareAdaptation(*m);
//#endif
		if (!indicator.isValid() || !indicator.isDefined(CELL))
		{
			indicator = m->CreateTag(indicator.isValid() ? indicator.GetTagName() : "indicator", DATA_INTEGER, CELL, NONE, 1);
			for (Storage::integer i = 0; i < m->CellLastLocalID(); ++i) if (m->isValidCell(i))
				indicator[m->CellByLocalID(i)] = 1;
		}
		for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
			(*it)->CoarseIndicator(*this,indicator);
		int coarse = 0, tot = m->TotalNumberOf(CELL);
		for (Storage::integer i = 0; i < m->CellLastLocalID(); ++i) if (m->isValidCell(i))
		{
			Cell c = m->CellByLocalID(i);
			if (GetLevel(c) && indicator[c] == 1)
				coarse++;
		}
		coarse = m->Integrate(coarse);
		if (m->GetProcessorRank() == 0)
			std::cout << __FUNCTION__ << " indicator marked " << coarse << "/" << tot << std::endl;
		if (!coarse) return false;
		m->ExchangeData(indicator, CELL);
		for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
			(*it)->BeginCoarsening();
        //return false;
		static int call_counter = 0;
		//return number of coarsened cells
		Storage::integer ret = 0;
		//initialize tree structure
		ENTER_BLOCK();
		PrepareSet();
		EXIT_BLOCK();

		//m->Save("before_coarse"+std::to_string(fi)+".pvtk");
		//std::cout << "Save before_coarse"+std::to_string(fi)+".pvtk" << std::endl;
		/*
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
		*/
		//m->Save("before_coarse_parent"+std::to_string(fi)+".pvtk");
		//std::cout << "Save before_coarse_parent"+std::to_string(fi)+".pvtk" << std::endl;
		
		//ENTER_BLOCK();
		//SynchronizeIndicated(indicator);
		//EXIT_BLOCK();
		
		CheckParentSet(__FILE__,__LINE__);
		
		int schedule_counter = 1; //indicates order in which refinement will be scheduled
		Storage::integer scheduled = 1, unscheduled = 0; //indicates that at least one element was scheduled on current sweep

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
#if defined(USE_OMP)
#pragma omp parallel for
#endif
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
#if defined(USE_OMP)
#pragma omp parallel for
#endif
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
#if defined(USE_OMP)
#pragma omp parallel for
#endif
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
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:unscheduled)
#endif
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
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(Storage::integer it = 0; it < m->EdgeLastLocalID(); ++it) if( m->isValidEdge(it) )
				coarse_indicator[m->EdgeByLocalID(it)] = 0;
			EXIT_BLOCK();
			//b) each cell mark it's finer edges with cell's schedule
			ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( indicator[c] )
				{
					ElementArray<Element> adj = c.getAdjElements(EDGE);
					for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
					{
						if (level[adj[kt]] > level[c]) //only finer edges
						{
#if defined(USE_OMP)
#pragma omp critical
#endif
							coarse_indicator[adj[kt]] = std::max(coarse_indicator[adj[kt]], indicator[c]);
						}
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
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:scheduled)
#endif
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
		
		
		//Order exchange of elements of the sets that are to be coarsened
		ENTER_BLOCK();
#if defined(USE_OMP)
#pragma omp parallel for //each set has unique set of elements, should be safe
#endif
		for(Storage::integer it = 0; it < m->EsetLastLocalID(); ++it) if( m->isValidElementSet(it) )
		{
			ElementSet set = m->EsetByLocalID(it);
			if( set.GetName().substr(0,3) == "AM_" )
			{
				if( indicator[set] != 0 && !set.Empty() )
					set.SynchronizeSetElements();
			}
		}
		EXIT_BLOCK();
		m->ExchangeMarked();


		m->ExchangeData(parent_set,CELL,0);
		m->ExchangeData(hanging_nodes,CELL | FACE,0);
		if(tri_hanging_edges.isValid())
			m->ExchangeData(tri_hanging_edges, CELL, 0);
		//m->ResolveSets();

		/*

		ENTER_BLOCK();
		for(Storage::integer it = 0; it < m->EsetLastLocalID(); ++it) if( m->isValidElementSet(it) )
		{
			ElementSet set = m->EsetByLocalID(it);
			if( set.GetName().substr(0,3) == "AM_" )
			{
				//int imax = -1, imin = INT_MAX;
				//for(ElementSet::iterator jt = set.Begin(); jt != set.End(); ++jt)
				//{
				//	imax = std::max(imax,indicator[*jt]);
				//	imin = std::min(imin,indicator[*jt]);
				//}
				//std::cout << "on proc " << m->GetProcessorRank() << " set " << set.GetName() << " size " << set.Size() << " set indicator " << indicator[set] << " elements indicator " << imin << ":" << imax;
				//if( set.HaveParent() ) std::cout << " parent " << set.GetParent().GetName();
				//std::cout << std::endl;
				//if( indicator[set] != 0 )
				//if( !set.HaveChild() )
				set.SynchronizeSetParents();
			}
		}
		EXIT_BLOCK();
		*/
		//m->Barrier();
		//std::cout << m->GetProcessorRank() << " call exchange marked" << std::endl;
		//m->ExchangeMarked();
		//std::cout << m->GetProcessorRank() << " finish exchange marked" << std::endl;
		//m->Barrier();
		ENTER_BLOCK();
		CheckParentSet(__FILE__,__LINE__);//,indicator);
		EXIT_BLOCK();
		//std::fstream fout("sets"+std::to_string(m->GetProcessorRank())+".txt",std::ios::out);
		//for(Mesh::iteratorSet it = m->BeginSet(); it != m->EndSet(); ++it)
		//	PrintSet(fout,it->self());
		ENTER_BLOCK();
		assert(Element::CheckConnectivity(m));
		CheckClosure(__FILE__,__LINE__);
		EXIT_BLOCK();
		//m->Save("unschdind"+std::to_string(fi)+".pvtk");
		//std::cout << "Save unschdind"+std::to_string(fi)+".pvtk" << std::endl;
		//Make schedule which elements should be refined earlier.
		ENTER_BLOCK();
		m->BeginModification();
		while(schedule_counter)
		{
			assert(Element::CheckConnectivity(m));
			CheckClosure(__FILE__,__LINE__);
			//CheckParentSet(__FILE__,__LINE__);//,indicator);
			//CheckParentSet();
			//fout << "schedule_counter " << schedule_counter << std::endl;
			//unite cells
			//should find and set hanging nodes on faces
			//find single node at the center, all other nodes,
			//adjacent over edge are hanging nodes
			ENTER_BLOCK();
			{
				int unites = 0;
				double t1, t2, tcollect = 0, tcenter = 0, thanging = 0, tunite = 0, thangcon = 0, tparent = 0, tdata = 0, tdelset = 0;
				for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
				{
					Cell c = m->CellByLocalID(it);
					if( !c.Hidden() && indicator[c] == schedule_counter )
					{
						t1 = Timer();
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
						t2 = Timer(), tcollect += t2 - t1, t1 = t2;
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
						for(kt = 1; kt < static_cast<Storage::integer>(unite_cells.size()); ++kt)
							center_node.Intersect(unite_cells[kt].getNodes());
						t2 = Timer(), tcenter += t2 - t1, t1 = t2;
						//fout << "nodes(" << center_node.size() << "):" << std::endl;
						//for(kt = 0; kt < center_node.size(); ++kt)
						//	fout << center_node[kt].Coords()[0] << " " << center_node[kt].Coords()[1] << " " << center_node[kt].Coords()[2] << std::endl;
						//fout << "child sets: ";
						//for(ElementSet chld = parent.GetChild(); chld.isValid(); chld = chld.GetSibling())
						//	fout << " " << chld.GetName();
						//fout << std::endl;
						//only one should be found
						/*
						if( center_node.size() != 1 )// || it == 127 )
						{
							INMOST_DATA_REAL_TYPE x[3];
							std::cout << "call_counter " << call_counter << " schedule_counter " << schedule_counter << std::endl;
							c.Centroid(x);
							std::cout << "parent set " << parent.GetName() << " size " << parent.Size() << " cell " << c.GlobalID() << " " << x[0] << " " << x[1] << " " << x[2] << std::endl;
							std::cout << "cells(" << unite_cells.size() << "):" << std::endl;
							for(kt = 0; kt < static_cast<Storage::integer>(unite_cells.size()); ++kt)
							{
								unite_cells[kt].Centroid(x);
								std::cout << unite_cells[kt].GlobalID() << " lid " << unite_cells[kt].LocalID();
								std::cout << " parent " << ElementSet(m,parent_set[unite_cells[kt]]).GetName();
								std::cout << " " << x[0] << " " << x[1] << " " << x[2];
								std::cout << std::endl;
							}
							std::cout << "nodes(" << center_node.size() << "):" << std::endl;
							for(kt = 0; kt < static_cast<Storage::integer>(center_node.size()); ++kt)
								std::cout << center_node[kt].Coords()[0] << " " << center_node[kt].Coords()[1] << " " << center_node[kt].Coords()[2] << std::endl;
							std::cout << "child sets: ";
							for(ElementSet chld = parent.GetChild(); chld.isValid(); chld = chld.GetSibling())
								std::cout << " " << chld.GetName();
							std::cout << std::endl;
							for(kt = 0; kt < static_cast<Storage::integer>(unite_cells.size()); ++kt)
							{
								ElementArray<Face> faces = unite_cells[kt].getFaces();
								std::cout << "cell " << unite_cells[kt].LocalID() << " faces " << faces.size() << " lc " << m->LowConn(unite_cells.at(kt)).size() << std::endl;
								std::cout << (unite_cells[kt].Closure() ? "closed" : "open") << " ";
								for(ElementArray<Face>::iterator qt = faces.begin(); qt != faces.end(); ++qt) std::cout << "f" << qt->LocalID() << " " << (qt->Closure() ? "closed" : "open") << " ";
								std::cout << std::endl;
								for(ElementArray<Face>::iterator qt = faces.begin(); qt != faces.end(); ++qt)
								{
									ElementArray<Edge> edges = qt->getEdges();
									for(ElementArray<Edge>::iterator mt = edges.begin(); mt != edges.end(); ++mt)
									{
										std::cout << "(" << mt->getBeg()->Coords()[0] << "," << mt->getBeg()->Coords()[1] << "," << mt->getBeg()->Coords()[2] << ")";
										std::cout << "<->";
										std::cout << "(" << mt->getEnd()->Coords()[0] << "," << mt->getEnd()->Coords()[1] << "," << mt->getEnd()->Coords()[2] << ")";
										std::cout << std::endl;
									}
								}
							}
						}
						*/
						assert(center_node.size() == 1);
						ElementArray<Node> hanging = center_node[0].BridgeAdjacencies2Node(EDGE);
						t2 = Timer(), thanging += t2 - t1, t1 = t2;
						Cell v = Cell::UniteCells(unite_cells,0);
						unites++;
						set_id[v] = set_id[parent];
						t2 = Timer(), tunite += t2 - t1, t1 = t2;
						//connect hanging nodes to the cell
						assert(hanging_nodes[v].size() == 0);
						for(ElementArray<Node>::size_type kt = 0; kt < hanging.size(); ++kt)
							hanging_nodes[v].push_back(hanging[kt]);
						t2 = Timer(), thangcon += t2 - t1, t1 = t2;
						//set new parent
						parent_set[v] = parent.GetParent().GetHandle();
						//add cell to parent set
						ElementSet(m,parent_set[v]).PutElement(v);
						t2 = Timer(), tparent += t2 - t1, t1 = t2;
						//set level for new cell
						level[v] = level[c]-1;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//						if( model ) model->CellCoarsening(unite_cells,v,parent);
//#endif
						for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
							(*it)->CellCoarsening(unite_cells,v,parent);
						t2 = Timer(), tdata += t2 - t1, t1 = t2;
						//~ v.Centroid(x);
						//fout << v.GlobalID() << " lid " << v.LocalID();
						//fout << " parent " << ElementSet(m,parent_set[v]).GetName();
						//fout << " " << x[0] << " " << x[1] << " " << x[2];
						//fout << " lvl " << level[v];
						//fout << std::endl;
						//delete set that contained cells
						//tree structure should be resolved on ApplyModification
						//fout << "delete set " << parent.GetName() << std::endl;
						parent.DeleteSet();
						t2 = Timer(), tdelset += t2 - t1, t1 = t2;
						//increment number of coarsened cells
						ret++;
					}
				}
				REPORT_VAL("collect", tcollect);
				REPORT_VAL("center node", tcenter);
				REPORT_VAL("collect faces", tcollect);
				REPORT_VAL("hanging",thanging);
				REPORT_VAL("unite", tunite);
				REPORT_VAL("connect hanging", thangcon);
				REPORT_VAL("reconnect parent", tparent);
				REPORT_VAL("data", tdata);
				REPORT_VAL("delete parent set", tdelset);
				REPORT_VAL("unites", unites);
			}
			EXIT_BLOCK();
			ENTER_BLOCK();
			assert(Element::CheckConnectivity(m));
			CheckClosure(__FILE__,__LINE__);
			EXIT_BLOCK();
			//unite faces
			//should find and set hanging nodes on edges
			//find single node at the center, all other nodes,
			//adjacent over edge of the face are hanging nodes
			int numcoarsened = 0;
			ENTER_BLOCK();
			{
				int unites = 0;
				double t1, t2, tadjcell = 0, tadjnode = 0, tcenter = 0, thanging = 0, tcollect = 0, tunite = 0, tdata = 0, thangcon = 0;
				for(Storage::integer it = 0; it < m->FaceLastLocalID(); ++it) if( m->isValidFace(it) )
				{
					Face f = m->FaceByLocalID(it);
					if( !f.Hidden() && indicator[f] == schedule_counter )
					{
						t1 = Timer();
						//one (or both) of the adjacent cells were coarsened and has lower level
						bool visited = false;
						(void)visited;
						ElementArray<Cell> cells = f.getCells();
						for(ElementArray<Cell>::size_type kt = 0; kt < cells.size(); ++kt)
						{
							assert(level[cells[kt]] < level[f]);
						}
						t2 = Timer(), tadjcell += t2 - t1;
						for(ElementArray<Cell>::size_type kt = 0; kt < cells.size(); ++kt)
						{
							if( level[cells[kt]] < level[f] )
							{
								t1 = Timer();
								//cell has one hanging node in common with current face
								ElementArray<Node> nodes = f.getNodes();
								t2 = Timer(), tadjnode += t2 - t1, t1 = t2;
								Storage::reference_array search_hanging = hanging_nodes[cells[kt]];
								nodes.Intersect(search_hanging.data(),search_hanging.size());
								assert(nodes.size() == 1);
								t2 = Timer(), tcenter += t2 - t1, t1 = t2;
								//faces that hanging node shares with the cell are
								//those to be united
								ElementArray<Face> unite_faces = cells[kt].getFaces();
								unite_faces.Intersect(nodes[0].getFaces());
								//unmark faces to prevent visit
								for(ElementArray<Face>::size_type lt = 0; lt < unite_faces.size(); ++lt)
									indicator[unite_faces[lt]] = 0;
								t2 = Timer(), tcollect += t2 - t1, t1 = t2;
								//nodes connected by edges to hanging node and
								//common to the cell are hanging nodes on edges
								ElementArray<Node> hanging = cells[kt].getNodes();
								hanging.Intersect(nodes[0].BridgeAdjacencies(EDGE,NODE));
								t2 = Timer(), thanging += t2 - t1, t1 = t2;
								//unite faces
								Face v = Face::UniteFaces(unite_faces,0);
								unites++;
								t2 = Timer(), tunite += t2 - t1, t1 = t2;
								//connect new face to hanging nodes
								for(ElementArray<Node>::size_type lt = 0; lt < hanging.size(); ++lt)
									hanging_nodes[v].push_back(hanging[lt]);
								t2 = Timer(), thangcon += t2 - t1, t1 = t2;
								//set level for new face
								level[v] = level[f]-1;
								visited = true;
								numcoarsened++;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//								if( model ) model->FaceCoarsening(unite_faces,v);
//#endif
								for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
									(*it)->FaceCoarsening(unite_faces,v);
								t2 = Timer(), tdata += t2 - t1, t1 = t2;
								break; //no need to visit the other cell
							}
						}
						assert(visited);
					}
				}
				REPORT_VAL("adj nodes", tadjnode);
				REPORT_VAL("adj cells", tadjcell);
				REPORT_VAL("center node", tcenter);
				REPORT_VAL("collect faces", tcollect);
				REPORT_VAL("hanging",thanging);
				REPORT_VAL("unite", tunite);
				REPORT_VAL("connect hanging", thangcon);
				REPORT_VAL("data", tdata);
				REPORT_VAL("unites", unites);
			}
			EXIT_BLOCK();
			ENTER_BLOCK();
			assert(Element::CheckConnectivity(m));
			CheckClosure(__FILE__,__LINE__);
			EXIT_BLOCK();
			//unite edges
			ENTER_BLOCK();
			{
				int unites = 0;
				double t1, t2, tadjface = 0, tadjnode = 0, tcenter = 0, tcollect = 0, tunite = 0, tdata = 0;
				for(Storage::integer it = 0; it < m->EdgeLastLocalID(); ++it) if( m->isValidEdge(it) )
				{
					Edge e = m->EdgeByLocalID(it);
					if( !e.Hidden() && indicator[e] == schedule_counter )
					{
						t1 = Timer();
						//at least one face must have lower level
						bool visited = false;
						(void)visited;
						ElementArray<Face> faces = e.getFaces();
						t2 = Timer(), tadjface += t2 - t1;
						for(ElementArray<Face>::size_type kt = 0; kt < faces.size(); ++kt)
						{
							if( level[faces[kt]] < level[e] )
							{
								t1 = Timer();
								//face has one hanging node in common with current edge
								ElementArray<Node> nodes = e.getNodes();
								t2 = Timer(), tadjnode += t2 - t1;
								Storage::reference_array search_hanging = hanging_nodes[faces[kt]];
								nodes.Intersect(search_hanging.data(),search_hanging.size());
								assert(nodes.size() == 1);
								t2 = Timer(), tcenter += t2 - t1;
								//edges that hanging node shares with the face are those to
								//be united
								ElementArray<Edge> unite_edges = faces[kt].getEdges();
								unite_edges.Intersect(nodes[0].getEdges());
								//unmark edges to prevent visit
								for(ElementArray<Edge>::size_type lt = 0; lt < unite_edges.size(); ++lt)
									indicator[unite_edges[lt]] = 0;
								t2 = Timer(), tcollect += t2 - t1;
								//unite edges
								Edge v = Edge::UniteEdges(unite_edges,0);
								unites++;
								t2 = Timer(), tunite += t2 - t1;
								//set level for new edge
								level[v] = level[e]-1;
								visited = true;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//								if( model ) model->EdgeCoarsening(unite_edges,v);
//#endif
								for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
									(*it)->EdgeCoarsening(unite_edges,v);
								t2 = Timer(), tdata += t2 - t1;
								break; //no need to visit any other face
							}
						}
						assert(visited);
					}
				}
				REPORT_VAL("adj faces", tadjface);
				REPORT_VAL("adj nodes", tadjnode);
				REPORT_VAL("center node", tcenter);
				REPORT_VAL("collect faces", tcollect);
				REPORT_VAL("unite", tunite);
				REPORT_VAL("data", tdata);
				REPORT_VAL("unites", unites);
			}
			EXIT_BLOCK();
			ENTER_BLOCK();
			assert(Element::CheckConnectivity(m));
			CheckClosure(__FILE__,__LINE__);
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
			REPORT_VAL("schedule counter", schedule_counter);
			//jump to later schedule
			schedule_counter--;
		}

		if( check_orientation )
		{
			ENTER_BLOCK();
			int nfixed = 0, nfixednew = 0, nfixedbnd = 0, nghost = 0;
			for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
				if( !it->CheckNormalOrientation() )
				{
					//it->FixNormalOrientation();
					nfixed++;
					if( it->New() )
						nfixednew++;
					if( it->Boundary() )
						nfixedbnd++;
					if( it->GetStatus() == Element::Ghost )
						nghost++;

				}
				//std::cout << "Face " << it->LocalID() << " oriented incorrectly " << std::endl;
			if( nfixed ) 
			{
				std::cout << __FILE__ << ":" << __LINE__ << " rank " << rank << " bad " << nfixed << " (new " << nfixednew << " bnd " << nfixedbnd << " ghost " << nghost << ") faces " << std::endl;
				REPORT_STR(rank << " bad " << nfixed << " faces");
			}
			EXIT_BLOCK();
		}
		if( check_convexity )
		{
			ENTER_BLOCK();
			int nbad = 0;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				if( !it->CheckConvexity() ) nbad++;
			if( nbad ) std::cout << __FILE__ << ":" << __LINE__ << " rank " << rank << " nonconvex cells: " << nbad << std::endl;
			EXIT_BLOCK();
		}

		//free created tag
		m->DeleteTag(indicator,FACE|EDGE);
		ENTER_BLOCK();
		m->Barrier();
		EXIT_BLOCK();
		//todo:
		m->ResolveModification();
		//todo:
		//let the user update their data
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//		if( model ) model->Adaptation(*m);
//#endif
		for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
			(*it)->Coarsening();
		m->ApplyModification();

		
		//done
		m->EndModification();
		EXIT_BLOCK();
		//fout.close();
		ENTER_BLOCK();
		assert(Element::CheckConnectivity(m));
		CheckClosure(__FILE__,__LINE__);
		EXIT_BLOCK();
		//restore links to prevent loss during balancing
		m->ExchangeData(parent_set,CELL,0);
		m->ExchangeData(hanging_nodes,CELL | FACE,0);
		if (tri_hanging_edges.isValid())
			m->ExchangeData(tri_hanging_edges, CELL, 0);
		
		/*
		ENTER_BLOCK();
		for(Storage::integer it = 0; it < m->EsetLastLocalID(); ++it) if( m->isValidElementSet(it) )
		{
			ElementSet set = m->EsetByLocalID(it);
			if( set.GetName().substr(0,3) == "AM_" )
			{
				//int imax = -1, imin = INT_MAX;
				//for(ElementSet::iterator jt = set.Begin(); jt != set.End(); ++jt)
				//{
				//	imax = std::max(imax,indicator[*jt]);
				//	imin = std::min(imin,indicator[*jt]);
				//}
				//std::cout << "on proc " << m->GetProcessorRank() << " set " << set.GetName() << " size " << set.Size() << " set indicator " << indicator[set] << " elements indicator " << imin << ":" << imax;
				//if( set.HaveParent() ) std::cout << " parent " << set.GetParent().GetName();
				//std::cout << std::endl;
				//if( !set.HaveChild() )
					set.SynchronizeSetParents();
			}
		}
		m->ExchangeMarked();
		EXIT_BLOCK();
		*/
		ENTER_BLOCK();
		CheckParentSet(__FILE__,__LINE__);
		EXIT_BLOCK();
		//m->Save("after_coarse"+std::to_string(fi)+".pvtk");
		//std::cout << "Save after_coarse"+std::to_string(fi)+".pvtk" << std::endl;
		//exit(-1);
		
		ENTER_BLOCK();
		//m->CheckCentroids(__FILE__,__LINE__);
		//CheckCentroids();
		//cleanup null links to hanging nodes
		for(ElementType etype = FACE; etype <= CELL; etype = NextElementType(etype))
		{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
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
		ENTER_BLOCK();
		CleanupSets(root);
		EXIT_BLOCK();
		
		//CheckParentSet();
		
		
		//restore face orientation
		//BUG: bad orientation not fixed automatically
		
		if( check_orientation )
		{
			ENTER_BLOCK();
			int nfixed = 0, nbnd = 0, nghost = 0;
			for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
				if( !it->CheckNormalOrientation() )
				{
					//it->FixNormalOrientation();
					nfixed++;
					if( it->Boundary() ) nbnd++;
					if( it->GetStatus() == Element::Ghost ) nghost++;
				}
				//std::cout << "Face " << it->LocalID() << " oriented incorrectly " << std::endl;
			if( nfixed ) 
			{
				std::cout << __FILE__ << ":" << __LINE__ << " rank " << rank << " bad " << nfixed << "(bnd " << nbnd << " ghost " << nghost << ") faces " << std::endl;
				REPORT_STR(rank << " bad " << nfixed << " faces");
			}
			EXIT_BLOCK();
		}
		
		if( check_convexity )
		{
			ENTER_BLOCK();
			int nbad = 0;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				if( !it->CheckConvexity() ) nbad++;
			if( nbad ) std::cout << __FILE__ << ":" << __LINE__ << " rank " << rank << " nonconvex cells: " << nbad << std::endl;
			EXIT_BLOCK();
		}
		 
		
		//reorder element's data to free up space
		ENTER_BLOCK();
		m->ReorderEmpty(CELL|FACE|EDGE|NODE|ESET);
		EXIT_BLOCK();
		ENTER_BLOCK();
		for (std::vector<AdaptiveMeshCallback*>::iterator it = callbacks.begin(); it != callbacks.end(); ++it)
			(*it)->EndCoarsening();
		EXIT_BLOCK();
		call_counter++;
		
		ret = m->Integrate(ret);
		REPORT_VAL("ret ",ret)
        EXIT_FUNC();
		return ret != 0;
	}
	
	void AdaptiveMesh::ComputeWeightCoarse(TagInteger indicator, TagReal wgt)
	{
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) wgt[*it] = 0;
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) if( it->GetStatus() != Element::Ghost )
		{
			if( indicator[*it] )
			{
				ElementSet parent(m,parent_set[*it]);
				for(ElementSet::iterator jt = parent.Begin(); jt != parent.End(); ++jt)
					if( jt->GetStatus() != Element::Ghost )
						wgt[*it]+=1;
			}
			else wgt[*it] = 0;
		}
		m->ReduceData(wgt,CELL,0,ReduceSum);
		m->ExchangeData(wgt,CELL,0);
	}
	
	void AdaptiveMesh::ComputeWeightRefine(TagInteger indicator, TagReal wgt)
	{
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		{
			if( indicator[*it] )
			{
				Storage::reference_array hanging;
				ElementArray<Node> nodes = it->getNodes();
				ElementArray<Face> faces = it->getFaces();
				hanging = hanging_nodes[*it];
				nodes.Subtract(hanging.data(),hanging.size());
				for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
				{
					Storage::reference_array hanging = hanging_nodes[*jt];
					nodes.Subtract(hanging.data(),hanging.size());
				}
				wgt[*it] = (INMOST_DATA_REAL_TYPE)nodes.size();
			}
			else wgt[*it] = 0;
		}
	}

	void AdaptiveMesh::CheckClosure(std::string file, int line)
	{
		ENTER_FUNC();
#if !defined(NDEBUG)
		for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
		{
			Cell c = m->CellByLocalID(it);
			assert(c.Hidden() || c.Closure());
		}
		for(Storage::integer it = 0; it < m->FaceLastLocalID(); ++it) if( m->isValidFace(it) )
		{
			Face c = m->FaceByLocalID(it);
			if( !c.Hidden() && !c.Closure() )
			{
				std::cout << "no closure face " << it << " at " << file << ":" << line << std::endl;
				ElementArray<Edge> edges = c.getEdges();
				std::cout << "edges: " << edges.size() << " lc: " << m->LowConn(c.GetHandle()).size() << std::endl;
				for(ElementArray<Edge>::iterator mt = edges.begin(); mt != edges.end(); ++mt)
				{
					std::cout << "(" << mt->getBeg()->Coords()[0] << "," << mt->getBeg()->Coords()[1] << "," << mt->getBeg()->Coords()[2] << ")";
					std::cout << "<->";
					std::cout << "(" << mt->getEnd()->Coords()[0] << "," << mt->getEnd()->Coords()[1] << "," << mt->getEnd()->Coords()[2] << ")";
					std::cout << std::endl;
				}
			}
			assert(c.Hidden() || c.Closure());
		}
#endif
		EXIT_FUNC();
	}
	
}
