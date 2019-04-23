#include "utils.h"

void RepartitionStatistics(INMOST::Mesh *m, std::ostream &os){
    int i;
    int NonGhostNum[4];
    int OwnedNum[4];
    int AllNum[4];
    NonGhostNum[0]=0;
    AllNum[0] = m->NumberOfNodes();
    for(INMOST::Mesh::iteratorNode it = m->BeginNode();it != m->EndNode();it++)
        if(it->GetStatus() != INMOST::Element::Ghost)
            NonGhostNum[0]++;
    NonGhostNum[1]=0;
    AllNum[1] = m->NumberOfEdges();
    for(INMOST::Mesh::iteratorEdge it = m->BeginEdge();it != m->EndEdge();it++)
        if(it->GetStatus() != INMOST::Element::Ghost)
            NonGhostNum[1]++;
    NonGhostNum[2]=0;
    AllNum[2] = m->NumberOfFaces();
    for(INMOST::Mesh::iteratorFace it = m->BeginFace();it != m->EndFace();it++)
        if(it->GetStatus() != INMOST::Element::Ghost)
            NonGhostNum[2]++;
    NonGhostNum[3]=0;
    AllNum[3] = m->NumberOfCells();
    for(INMOST::Mesh::iteratorCell it = m->BeginCell();it != m->EndCell();it++)
        if(it->GetStatus() != INMOST::Element::Ghost)
            NonGhostNum[3]++;


    OwnedNum[0]=0;
    for(INMOST::Mesh::iteratorNode it = m->BeginNode();it != m->EndNode();it++)
        if(it->GetStatus() == INMOST::Element::Owned)
            OwnedNum[0]++;
    OwnedNum[1]=0;

    for(INMOST::Mesh::iteratorEdge it = m->BeginEdge();it != m->EndEdge();it++)
        if(it->GetStatus() == INMOST::Element::Owned)
            OwnedNum[1]++;
    OwnedNum[2]=0;
    for(INMOST::Mesh::iteratorFace it = m->BeginFace();it != m->EndFace();it++)
        if(it->GetStatus() == INMOST::Element::Owned)
            OwnedNum[2]++;
    OwnedNum[3]=0;
    for(INMOST::Mesh::iteratorCell it = m->BeginCell();it != m->EndCell();it++)
        if(it->GetStatus() == INMOST::Element::Owned)
            OwnedNum[3]++;
    int MinOwn[4],MaxOwn[4],AvOwn[4];
    for(i=0;i<4;i++){
        MinOwn[i] = m->AggregateMin(OwnedNum[i]);
        MaxOwn[i] = m->AggregateMax(OwnedNum[i]);
        AvOwn[i] = m->Integrate(OwnedNum[i]);
        AvOwn[i] /= m->GetProcessorsNumber();
    }
    if(m->GetProcessorRank() == 0) {
        os << "Owned Nodes statistics : ";
        os << "Min " << MinOwn[0] << " Max " << MaxOwn[0] << " Average " << AvOwn[0] << std::endl;
        os << "Owned Edges statistics : " ;
        os << "Min " << MinOwn[1] << " Max " << MaxOwn[1] << " Average " << AvOwn[1] << std::endl;
        os << "Owned Faces statistics : " ;
        os << "Min " << MinOwn[2] << " Max " << MaxOwn[2] << " Average " << AvOwn[2] << std::endl;
        os << "Owned Cells statistics : ";
        os << "Min " << MinOwn[3] << " Max " << MaxOwn[3] << " Average " << AvOwn[3] << std::endl;
    }


    int Min[4],Max[4],Av[4];
    for(i=0;i<4;i++){
        Min[i] = m->AggregateMin(NonGhostNum[i]);
        Max[i] = m->AggregateMax(NonGhostNum[i]);
        Av[i] = m->Integrate(NonGhostNum[i]);
        Av[i] /= m->GetProcessorsNumber();
    }
    if(m->GetProcessorRank() == 0) {
        os << "Nodes statistics " << std::endl;
        os << "Min " << Min[0] << " Max " << Max[0] << " Average " << Av[0] << std::endl;
        os << "Edges statistics " << std::endl;
        os << "Min " << Min[1] << " Max " << Max[1] << " Average " << Av[1] << std::endl;
        os << "Faces statistics " << std::endl;
        os << "Min " << Min[2] << " Max " << Max[2] << " Average " << Av[2] << std::endl;
        os << "Cells statistics " << std::endl;
        os << "Min " << Min[3] << " Max " << Max[3] << " Average " << Av[3] << std::endl;
    }
}

void SaveMatrix(INMOST::Sparse::Matrix &mat, int MaxInd,MPI_Comm comm) {
    int processRank = 0, processorsCount = 1;
    int i;
    MPI_Comm_rank(comm, &processRank);  // Get the rank of the current process
    MPI_Comm_size(comm, &processorsCount); // Get the total number of processors used
    std::string name = mat.GetName();
    std::string filename = name + ".mtx";
    std::string ordname = name + ".ord";
    int *buf;
    if (processRank == 0){
        buf = (int *) malloc(processorsCount * sizeof(int));
    }
    std::cout<<"tst"<<std::endl;
    if(processorsCount>1)
        MPI_Gather(&MaxInd, 1, MPI_INT, buf, 1, MPI_INT, 0, comm);
    else
        buf[0] = MaxInd;
    std::cout<<"save"<<std::endl;
    mat.Save(filename.c_str());

    if(processRank == 0) {
        std::ofstream f_ord;
        f_ord.open(ordname.c_str());
        f_ord<<buf[processorsCount-1]+1<<std::endl;
        for ( i = 1; i <= buf[processorsCount-1]+1; i++) {
            f_ord<<i<<std::endl;
        }

        f_ord << processorsCount << std::endl;
        f_ord << 0 << std::endl;
        for ( i = 0; i < processorsCount; i++) {
            f_ord << buf[i] + 1 << std::endl;
        }
        f_ord << processorsCount << std::endl;
        f_ord.close();

        free(buf);
    }
}




void PrintCurrentMemory(std::ostream &os, INMOST::Mesh *m){
/*
    unsigned long int memory_test = getCurrentRSS();
    std::cout<<"Processor rank "<< m->GetProcessorRank()<<" current memory usage "<<memory_test << " B, or "<<memory_test/(1024L)
             <<" KB, or "<<memory_test/(1024L*1024L)<<" MB, or "<<memory_test/(1024L*1024L* 1024L) <<" GB "<<std::endl;
    unsigned long int total_memory;
    if (m->GetProcessorsNumber() > 1) {
        MPI_Barrier(INMOST_MPI_COMM_WORLD);
        MPI_Reduce(&memory_test, &total_memory, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, INMOST_MPI_COMM_WORLD);
        MPI_Barrier(INMOST_MPI_COMM_WORLD);
    }
    else
        total_memory = memory_test;
    if(m->GetProcessorRank()==0){
        std::cout<<"total current memory usage "<<total_memory << " B, or "<<total_memory/(1024L)
                 <<" KB, or "<<total_memory/(1024L*1024L)<<" MB, or "<<total_memory/(1024L*1024L* 1024L) <<" GB "<<std::endl;
    }

*/
}
void PrintPeakMemory(std::ostream &os, INMOST::Mesh *m) {
    /*
    unsigned long int memory_test = getPeakRSS();
    std::cout << "Processor rank " << m->GetProcessorRank() << " current memory usage " << memory_test << " B, or "
              << memory_test / (1024L)
              << " KB, or " << memory_test / (1024L * 1024L) << " MB, or " << memory_test / (1024L * 1024L * 1024L)
              << " GB " << std::endl;
    unsigned long int total_memory;// = m->Integrate(memory_test);
    if (m->GetProcessorsNumber() > 1) {
        MPI_Barrier(INMOST_MPI_COMM_WORLD);
        MPI_Reduce(&memory_test, &total_memory, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, INMOST_MPI_COMM_WORLD);
        MPI_Barrier(INMOST_MPI_COMM_WORLD);
    }
    else
        total_memory = memory_test;
    if(m->GetProcessorRank()==0){
        std::cout<<"total peak memory usage "<<total_memory << " B, or "<<total_memory/(1024L)
                 <<" KB, or "<<total_memory/(1024L*1024L)<<" MB, or "<<total_memory/(1024L*1024L* 1024L) <<" GB "<<std::endl;
    }*/
}




void WriteTags_element(std::string filename, INMOST::Mesh *m, std::vector<INMOST::Tag> tags, unsigned int Mask) {
    int i,j;
    INMOST::Tag set_id = m->CreateTag("TEMPORARY_ELEMENT_ID",INMOST::DATA_INTEGER, Mask | INMOST::NODE,INMOST::NONE,1);
    INMOST::Storage::integer cur_num = 0;
    for(INMOST::Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it) it->IntegerDF(set_id) = cur_num++;
    cur_num = 0;
    for(INMOST::Mesh::iteratorElement it = m->BeginElement(Mask); it != m->EndElement(); ++it) it->IntegerDF(set_id) = cur_num++;


    std::fstream fs;
    fs.open(filename.c_str(), std::fstream::out);
    if (!fs.is_open())
        std::cout << "DIDNT OPEN file "<<filename << "\n";
    fs<<std::setprecision(10);
    fs << "# vtk DataFile Version 2.0 "<<std::endl << "file is written by INMOST"<<std::endl << "ASCII"<<std::endl << "DATASET UNSTRUCTURED_GRID"<<std::endl;
    fs << "POINTS " << m->NumberOfNodes() << " double"<<std::endl;
    cur_num=0;
    for(INMOST::Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it) {
        fs << it->Coords()[0] << " " << it->Coords()[1] << " " << it->Coords()[2] << std::endl;
    }
    fs<<std::endl;



 //   int size_of_el = m->BeginElement(Mask)->getNodes().size();
    int num_elem = m->NumberOf(Mask);
    //std::cout<<"in in "<<std::endl;
/*
    fs << "CELLS " << num_elem << " " << (size_of_el+1)*num_elem << std::endl;
    for(INMOST::Mesh::iteratorElement it = m->BeginElement(Mask); it != m->EndElement(); it++)
    {
        fs << size_of_el<<" ";
        INMOST::ElementArray<INMOST::Node> nodes = it->getNodes();
        for(j=0;j<size_of_el;j++)
            fs<< nodes[j].Integer(set_id)<< " ";// << m->CellByLocalID(i).getNodes()[1].LocalID() << " " << m->CellByLocalID(i).getNodes()[2].LocalID()<< std::endl;;// tris[i].get_vert(0)->get_num() << " " << tris[i].get_vert(1)->get_num() << " " << tris[i].get_vert(2)->get_num() << std::endl;
        fs<<std::endl;
    }
    fs<<std::endl;*/

    long int count_cells=0;
    for(INMOST::Mesh::iteratorElement it = m->BeginElement(Mask); it != m->EndElement(); it++)
    {

        switch (it->GetGeometricType()) {
            case INMOST::Element::Vertex:
            case INMOST::Element::Line:
            case INMOST::Element::Tri:
            case INMOST::Element::Quad:
            case INMOST::Element::Polygon:
            case INMOST::Element::MultiLine:
            case INMOST::Element::Tet:
            case INMOST::Element::Hex :
            case INMOST::Element::Prism:
            case INMOST::Element::Pyramid : {
                count_cells += (it->getNodes().size() + 1);
                break;
            }
            case INMOST::Element::Polyhedron: {
                INMOST::ElementArray<INMOST::Face> faces = it->getFaces();
                count_cells += faces.size();
                for (i = 0; i < faces.size(); i++) {
                    count_cells += 1 + faces[i].getNodes().size();
                }
                break;
            }
            case INMOST::Element::Unset :
            default:
            std::cout<<m->GetProcessorRank()<<" proc have fail element num "<<it->Integer(set_id)<<std::endl;
                abort();
        }
    }
    fs << "CELLS " << num_elem << " " << count_cells << std::endl;
    for(INMOST::Mesh::iteratorElement it = m->BeginElement(Mask); it != m->EndElement(); it++)
    {
        switch (it->GetGeometricType()) {
            case INMOST::Element::Vertex:
            case INMOST::Element::Line:
            case INMOST::Element::Tri:
            case INMOST::Element::Quad:
            case INMOST::Element::Polygon:
            case INMOST::Element::MultiLine:
            case INMOST::Element::Tet:
            case INMOST::Element::Hex :
            case INMOST::Element::Prism:
            case INMOST::Element::Pyramid : {
                INMOST::ElementArray<INMOST::Node> nodes = it->getNodes();
                fs << nodes.size() << " ";
                for (j = 0; j < nodes.size(); j++)
                    fs << nodes[j].Integer(set_id)<< " ";
                       // << m->CellByLocalID(i).getNodes()[1].LocalID() << " " << m->CellByLocalID(i).getNodes()[2].LocalID()<< std::endl;;// tris[i].get_vert(0)->get_num() << " " << tris[i].get_vert(1)->get_num() << " " << tris[i].get_vert(2)->get_num() << std::endl;
                fs << std::endl;
                break;
            }
            case INMOST::Element::Polyhedron: {
                INMOST::ElementArray<INMOST::Face> faces = it->getFaces();
                fs << faces.size() << " ";
                for (i = 0; i < faces.size(); i++)
                    fs << faces[i].getNodes().size() << " ";
                for (i = 0; i < faces.size(); i++) {
                    INMOST::ElementArray<INMOST::Node> nodes = faces[i].getNodes();
                    for (j = 0; j < nodes.size(); j++)
                        fs << nodes[j].Integer(set_id)<< " ";
                           // << m->CellByLocalID(i).getNodes()[1].LocalID() << " " << m->CellByLocalID(i).getNodes()[2].LocalID()<< std::endl;;// tris[i].get_vert(0)->get_num() << " " << tris[i].get_vert(1)->get_num() << " " << tris[i].get_vert(2)->get_num() << std::endl;
                }
                fs << std::endl;
                break;
            }
            case INMOST::Element::Unset :
            default:
                std::cout<<m->GetProcessorRank()<<" proc have fail element num "<<it->Integer(set_id)<<std::endl;
                abort();

        }
    }
    fs<<std::endl;

   // std::cout<<"in in "<<std::endl;
    fs << "CELL_TYPES " << num_elem << std::endl;





   // MPI_Barrier(INMOST_MPI_COMM_WORLD);
    /*
    int vtk_label;
    switch (Mask){
        case INMOST::CELL :
            vtk_label = 10;
            break;
        case INMOST::FACE:
            vtk_label = 5;
            break;
        case INMOST::EDGE:
            vtk_label = 3;
            break;
        default:
        std::cout<<"wrong mask "<<Mask<<std::endl;
            abort();
    }
    for (i = 0; i < num_elem; i++){
        fs << vtk_label<<std::endl;
    }
*/

//    fs << "CELL_TYPES " << num_elem << std::endl;
    MPI_Barrier(INMOST_MPI_COMM_WORLD);
    for(INMOST::Mesh::iteratorElement it = m->BeginElement(Mask); it != m->EndElement(); it++) {

        switch (it->GetGeometricType()) {
            case INMOST::Element::Vertex:
                fs<< 1<<" ";
                break;
            case INMOST::Element::Line:
                fs<< 3<<" ";
                break;
            case INMOST::Element::Tri:
                fs<< 5<<" ";
                break;
            case INMOST::Element::Quad:
                fs<< 9<<" ";
                break;
            case INMOST::Element::Polygon:
                fs<< 7<<" ";
                break;
            case INMOST::Element::MultiLine:
                fs<< 4<<" ";
                break;
            case INMOST::Element::Tet:
                fs<< 10<<" ";
                break;
            case INMOST::Element::Hex :
                fs<< 12<<" ";
                break;
            case INMOST::Element::Prism:
                fs<< 13<<" ";
                break;
            case INMOST::Element::Pyramid :
                fs<< 14<<" ";
                break;
            case INMOST::Element::Polyhedron:
                fs<< 42<<" ";
                break;
        }
        fs<<std::endl;

    }
    MPI_Barrier(INMOST_MPI_COMM_WORLD);
    fs<<std::endl;
    for(i=0;i<tags.size();i++)
        if(tags[i].isDefined(Mask)) {
            fs << "CELL_DATA  " << num_elem << std::endl;
            break;
        }
        else{
            //std::cout<<"none tags are defined with this mask "<<Mask<<std::endl;
        }

    MPI_Barrier(INMOST_MPI_COMM_WORLD);

    for(i=0;i<tags.size();i++)
        if(tags[i].isDefined(Mask) ){
           // if(tags[i].GetSize()>1){
              //  std::cout<<"mesh proc "<<m->GetProcessorRank()<<" tag "<<tags[i].GetTagName()<< " has size more than one, "<<tags[i].GetSize()<<std::endl;
              //  continue;
           // }
            if (tags[i].GetSize()== 1) {
                fs << "SCALARS " << tags[i].GetTagName();
                if (tags[i].GetDataType() == INMOST::DATA_INTEGER) {
                    fs << " int " << std::endl << "LOOKUP_TABLE default " << std::endl;

                    for (INMOST::Mesh::iteratorElement it = m->BeginElement(Mask); it != m->EndElement(); it++) {
                        fs << it->Integer(tags[i]) << std::endl;
                    }
                } else {
                    if (tags[i].GetDataType() == INMOST::DATA_REAL) {
                        fs << " double  " << std::endl << "LOOKUP_TABLE default " << std::endl;
                        for (INMOST::Mesh::iteratorElement it = m->BeginElement(Mask); it != m->EndElement(); it++) {
                            fs << it->Real(tags[i]) << std::endl;
                        }
                    } else {
                        std::cerr << "Wrong datatype for writing " << tags[i].GetDataType() << std::endl;
                        exit(0);
                    }
                }
            }
            else{
                for(j=0;j<tags[i].GetSize();j++){
                    std::stringstream ss;
                    std::string name_tag;
                    ss << j;
                    name_tag=tags[i].GetTagName() + "_"+ss.str();
                    fs << "SCALARS " << name_tag;
                    if (tags[i].GetDataType() == INMOST::DATA_INTEGER) {
                        fs << " int " << std::endl << "LOOKUP_TABLE default " << std::endl;

                        for (INMOST::Mesh::iteratorElement it = m->BeginElement(Mask); it != m->EndElement(); it++) {
                            fs << it->IntegerArray(tags[i])[j] << std::endl;
                        }
                    } else {
                        if (tags[i].GetDataType() == INMOST::DATA_REAL) {
                            fs << " double  " << std::endl << "LOOKUP_TABLE default " << std::endl;
                            for (INMOST::Mesh::iteratorElement it = m->BeginElement(Mask); it != m->EndElement(); it++) {
                                fs << it->RealArray(tags[i])[j] << std::endl;
                            }
                        } else {
                            std::cerr << "Wrong datatype for writing " << tags[i].GetDataType() << std::endl;
                            exit(0);
                        }
                    }

                }

            }

        }

    MPI_Barrier(INMOST_MPI_COMM_WORLD);

    for(i=0;i<tags.size();i++)
        if(tags[i].isDefined(INMOST::NODE)) {
            fs << "POINT_DATA  " << m->NumberOfNodes() << std::endl;
            break;
        }
    MPI_Barrier(INMOST_MPI_COMM_WORLD);

    for(i=0;i<tags.size();i++)
        if(tags[i].isDefined(INMOST::NODE)){

            //if(tags[i].GetSize()>1){
             //   std::cout<<"mesh proc "<<m->GetProcessorRank()<<" tag "<<tags[i].GetTagName()<< " has size more than one, "<<tags[i].GetSize()<<std::endl;
             //   continue;
            //}
            if(tags[i].GetSize()==1) {
                fs << "SCALARS " << tags[i].GetTagName();
                if (tags[i].GetDataType() == INMOST::DATA_INTEGER) {
                    fs << " int  " << std::endl << "LOOKUP_TABLE default " << std::endl;
                    for (INMOST::Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); it++) {
                        fs << it->Integer(tags[i]) << std::endl;
                    }

                } else {
                    if (tags[i].GetDataType() == INMOST::DATA_REAL) {
                        fs << " double  " << std::endl << "LOOKUP_TABLE default " << std::endl;
                        for (INMOST::Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); it++) {
                            fs << it->Real(tags[i]) << std::endl;
                        }
                    } else {
                        std::cerr << "Wrong datatype for writing " << tags[i].GetDataType() << std::endl;
                        exit(0);
                    }
                }
            }
            else{
                for(j=0;j<tags[i].GetSize();j++){
                    std::stringstream ss;
                    std::string name_tag;
                    ss << j;
                    name_tag=tags[i].GetTagName() + "_"+ss.str();
                    fs << "SCALARS " << name_tag;
                    if (tags[i].GetDataType() == INMOST::DATA_INTEGER) {
                        fs << " int  " << std::endl << "LOOKUP_TABLE default " << std::endl;
                        for (INMOST::Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); it++) {
                            fs << it->IntegerArray(tags[i])[j] << std::endl;
                        }

                    } else {
                        if (tags[i].GetDataType() == INMOST::DATA_REAL) {
                            fs << " double  " << std::endl << "LOOKUP_TABLE default " << std::endl;
                            for (INMOST::Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); it++) {
                                fs << it->RealArray(tags[i])[j] << std::endl;
                            }
                        } else {
                            std::cerr << "Wrong datatype for writing " << tags[i].GetDataType() << std::endl;
                            exit(0);
                        }
                    }

                }
            }
        }

    MPI_Barrier(INMOST_MPI_COMM_WORLD);

    m->DeleteTag(set_id,Mask | INMOST::NODE);
}



void WriteTags_vtk(std::string filename, INMOST::Mesh *m, std::vector<INMOST::Tag> tags) {
    bool cell = false, face = false, edge = false;

    int i;
    for (i = 0; i < tags.size(); i++) {
        //std::cout<<"test bools inner, cell "<<tags[i].isDefined(INMOST::CELL)<<" face "<<tags[i].isDefined(INMOST::FACE)<<" edge "<<tags[i].isDefined(INMOST::EDGE)<<std::endl;
        //std::cout<<tags[i].GetTagName()<<" "<<tags[i].GetDataType()<<" "<<std::endl;
        if (tags[i].isDefined(INMOST::CELL))
            cell = true;
        if (tags[i].isDefined(INMOST::FACE))
            face = true;
        if (tags[i].isDefined(INMOST::EDGE))
            edge = true;
    }
    //std::cout<<"test bools, cell "<<cell<<" face "<<face<<" edge "<<edge<<std::endl;
   /* if (filename.find(".vtk") == std::string::npos) {
        std::cout << "fail to save vtk file" << std::endl;
        abort();
    }*/
    std::string name = filename;
   // std::string::size_type pos = name.rfind(".vtk");
 //   name.erase(pos);
    std::string::size_type l = name.find_last_of("/\\");
    std::string fname = name.substr(l + 1, name.length());
   // std::cout<<"in write tags"<<std::endl;
    //if(cell)
    {
        std::string end=fname+"_cells.vtk";
        WriteTags_element(end,m,tags,INMOST::CELL);
    }
    if(face){
        std::string end=fname+"_faces.vtk";
        WriteTags_element(end,m,tags,INMOST::FACE);
    }
    if(edge){
        std::string end=fname+"_edges.vtk";
        WriteTags_element(end,m,tags,INMOST::EDGE);
    }
}


void WriteTags_pvtk(std::string filename, INMOST::Mesh *m, std::vector<INMOST::Tag> tags) {
    bool cell = false, face = false, edge = false;

    int i;
    for (i = 0; i < tags.size(); i++) {
        if (tags[i].isDefined(INMOST::CELL))
            cell = true;
        if (tags[i].isDefined(INMOST::FACE))
            face = true;
        if (tags[i].isDefined(INMOST::EDGE))
            edge = true;
    }
   /* if (filename.find(".pvtk") == std::string::npos) {
        std::cout << "fail to save pvtk file" << std::endl;
        abort();
    }*/
    std::string name = filename;
  //  std::string::size_type pos = name.rfind(".pvtk");
  //  name.erase(pos);
    std::string::size_type l = name.find_last_of("/\\");
    std::string fname_init = name.substr(l + 1, name.length());

//    if(cell)
    {
        std::string fname = fname_init+"_cells";
        std::string f_pvtk = fname+ ".pvtk";
        if(m->GetProcessorRank()==0)
        {//create pvtk file
            std::stringstream ss;
            std::string numproc;
            ss << m->GetProcessorsNumber();
            numproc=ss.str();
            FILE  *f=fopen(f_pvtk.c_str(), "w");
            fprintf(f,"%s%s%s", "<File version=\"pvtk-1.0\" dataType=\"vtkUnstructuredGrid\" numberOfPieces=\"", numproc.c_str(),"\">\n");
            for (int i=0; i<m->GetProcessorsNumber(); i++)
            {
                std::stringstream ss;
                ss << i;//task: leading zeroes
                std::string end="_"+ss.str()+".vtk";

                std::string temp=fname;
                temp.append(end);
                fprintf(f, "%s%s%s", "<Piece fileName=\"" ,temp.c_str() ,"\"/>\n");
            }
            fprintf(f, "%s", "</File>");
            fclose(f);
        }
        std::stringstream ss;
        ss << m->GetProcessorRank();//task: leading zeroes
        std::string end="_"+ss.str()+".vtk";

        fname.append(end);
        WriteTags_element(fname,m,tags, INMOST::CELL);
    }
    if(edge){
        std::string fname = fname_init+"_edges";
        std::string f_pvtk = fname+ ".pvtk";
        if(m->GetProcessorRank()==0)
        {//create pvtk file
            std::stringstream ss;
            std::string numproc;
            ss << m->GetProcessorsNumber();
            numproc=ss.str();
            FILE  *f=fopen(f_pvtk.c_str(), "w");
            fprintf(f,"%s%s%s", "<File version=\"pvtk-1.0\" dataType=\"vtkUnstructuredGrid\" numberOfPieces=\"", numproc.c_str(),"\">\n");
            for (int i=0; i<m->GetProcessorsNumber(); i++)
            {
                std::stringstream ss;
                ss << i;//task: leading zeroes
                std::string end="_"+ss.str()+".vtk";

                std::string temp=fname;
                temp.append(end);
                fprintf(f, "%s%s%s", "<Piece fileName=\"" ,temp.c_str() ,"\"/>\n");
            }
            fprintf(f, "%s", "</File>");
            fclose(f);
        }
        std::stringstream ss;
        ss << m->GetProcessorRank();//task: leading zeroes
        std::string end="_"+ss.str()+".vtk";

        fname.append(end);
        WriteTags_element(fname,m,tags, INMOST::EDGE);
    }
    if(face){
        std::string fname = fname_init+"_faces";
        std::string f_pvtk = fname+ ".pvtk";
        if(m->GetProcessorRank()==0)
        {//create pvtk file
            std::stringstream ss;
            std::string numproc;
            ss << m->GetProcessorsNumber();
            numproc=ss.str();
            FILE  *f=fopen(f_pvtk.c_str(), "w");
            fprintf(f,"%s%s%s", "<File version=\"pvtk-1.0\" dataType=\"vtkUnstructuredGrid\" numberOfPieces=\"", numproc.c_str(),"\">\n");
            for (int i=0; i<m->GetProcessorsNumber(); i++)
            {
                std::stringstream ss;
                ss << i;//task: leading zeroes
                std::string end="_"+ss.str()+".vtk";

                std::string temp=fname;
                temp.append(end);
                fprintf(f, "%s%s%s", "<Piece fileName=\"" ,temp.c_str() ,"\"/>\n");
            }
            fprintf(f, "%s", "</File>");
            fclose(f);
        }
        std::stringstream ss;
        ss << m->GetProcessorRank();//task: leading zeroes
        std::string end="_"+ss.str()+".vtk";

        fname.append(end);
        WriteTags_element(fname,m,tags, INMOST::FACE);
    }
}

void WriteTags(std::string filename, INMOST::Mesh *m, std::vector<INMOST::Tag> tags){

    if(m->GetProcessorsNumber()==1){
        //std::string filename = "solution.vtk";
        WriteTags_vtk(filename,m,tags);
       // WriteTags_element(filename,m,tags,INMOST::CELL);
    }
    else{
        WriteTags_pvtk(filename,m,tags);
    }

}


//==============================================================================================================================================

void Write_global(std::string filename, INMOST::Mesh *m ) {
       bool cell = false, face = false, edge = false;

    int i;
    std::vector<INMOST::Tag> tags;

    std::vector<std::string> names;
    m->ListTagNames(names);
    for(i=0;i<names.size();i++){
        //std::cout<<names[i]<<std::endl;
        if(m->HaveTag(names[i])){
            INMOST::Tag tg = m->GetTag(names[i]);
           // std::cout<<tg.isValid()<<" "<<(tg.isSparse(INMOST::NODE))<<" "<<(tg.isSparse( INMOST::EDGE))<<" "<<(tg.isSparse(INMOST::FACE))<< " "<<(tg.isSparse(INMOST::CELL))<<std::endl;
           // std::cout<<( (tg.GetDataType() ==INMOST::DATA_REAL ) || (tg.GetDataType() ==INMOST::DATA_INTEGER ))<<" "<<(names[i].find("PROTECTED_") != std::string::npos)<<std::endl;
            if(tg.isValid()  &&
                    ( (tg.GetDataType() ==INMOST::DATA_REAL ) || (tg.GetDataType() ==INMOST::DATA_INTEGER )) &&
                    (tg != m->CoordsTag()) &&
                    (tg != m->SharedTag()) &&
                    (tg != m->SendtoTag()) &&
                    (tg != m->ProcessorsTag()) && (tg.GetTagName() != "TEMPORARY_ON_SKIN") )
            {
                std::cout<<tg.GetTagName()<<std::endl;
                    tags.push_back(tg);
            }
        }
    }

    for (i = 0; i < tags.size(); i++) {
        //std::cout<<"test bools inner, cell "<<tags[i].isDefined(INMOST::CELL)<<" face "<<tags[i].isDefined(INMOST::FACE)<<" edge "<<tags[i].isDefined(INMOST::EDGE)<<std::endl;
        //std::cout<<tags[i].GetTagName()<<" "<<tags[i].GetDataType()<<" "<<std::endl;
        if (tags[i].isDefined(INMOST::CELL))
            cell = true;
        if (tags[i].isDefined(INMOST::FACE))
            face = true;
        if (tags[i].isDefined(INMOST::EDGE))
            edge = true;
    }
    //std::cout<<"test bools, cell "<<cell<<" face "<<face<<" edge "<<edge<<std::endl;
    if (filename.find(".vtk") == std::string::npos) {
        std::cout << "fail to save vtk file" << std::endl;
        abort();
    }
    std::string name = filename;
    std::string::size_type pos = name.rfind(".vtk");
    name.erase(pos);
    std::string::size_type l = name.find_last_of("/\\");
    std::string fname = name.substr(l + 1, name.length());
   // std::cout<<"in write tags gl"<<std::endl;
    //if(cell)
    {
        std::string end=fname+"_cells.vtk";
        WriteTags_element(end,m,tags,INMOST::CELL);
    }
    if(face){
        std::string end=fname+"_faces.vtk";
        WriteTags_element(end,m,tags,INMOST::FACE);
    }
    if(edge){
        std::string end=fname+"_edges.vtk";
        WriteTags_element(end,m,tags,INMOST::EDGE);
    }
}
