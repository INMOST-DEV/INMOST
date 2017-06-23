#include "inmost_ani_mesh.h"



extern "C" {
    void refine_(int *number_of_refines ,int *nv ,int *nvmax ,int *nt ,int *ntmax ,int *nb ,int *nbmax,
                 double *vrt, int *tet, int *bnd, int *material, int *labelF,
                    double *MapMtr, int* Ref2MapMtr, int *iW, int *MaxWi,int* flag_write);

    void loadmani_(
         int * nvmax,int * nbmax,int * ntmax,
         int * nv,int * nb,int * nt,
         double * vrt,int * bnd,int * tet,int * labelF,int * material,
         int *nvfix,int * nbfix,int * ntfix,int * ivfix,int * ibfix,int * itfix,
         int *iW,int * iW1, const char * name);

}



Ani_Mesh::Ani_Mesh(){
    nv=0;nt=0;nb=0; n_curv=0; nvfix=0; nbfix=0; ntfix   =0;
    vrt=NULL; bnd = NULL; tet = NULL; labelF=NULL; labelT=NULL;
    ivfix = NULL; ibfix = NULL; itfix = NULL;
    nvfix = 0; nbfix = 0; ntfix = 0;
}
Ani_Mesh::~Ani_Mesh()
{
#ifdef STATISTICS
    //std::cout<<"Animesh destructor"<<std::endl;
#endif
	if(vrt !=NULL) {
#ifdef STATISTICS
//        std::cout<<"free vrt"<<std::endl;
#endif
        free(vrt);
    }
	if(tet !=NULL) {
#ifdef STATISTICS
//        std::cout<<"free tet"<<std::endl;
#endif
        free(tet);
    }
	if(labelT !=NULL) {
#ifdef STATISTICS
//        std::cout<<"free labelt"<<std::endl;
#endif
        free(labelT);
    }
	if(bnd !=NULL) {
#ifdef STATISTICS
//        std::cout<<"free bnd"<<std::endl;
#endif
        free(bnd);
    }
	if(labelF !=NULL) {
#ifdef STATISTICS
//        std::cout<<"free labelF"<<std::endl;
#endif
        free(labelF);
    }
	if(ivfix !=NULL) {
        free(ivfix);
    }
	if(ibfix !=NULL) {
        free(ibfix);
    }
	if(itfix !=NULL) {
        free(itfix);
    }

}

unsigned long int  Ani_Mesh::getsize(){
    return nvmax*3*sizeof(double) + (4+1)*ntmax*sizeof(int) + (3+1)*nbmax*sizeof(int) +
            (nvfix*10)*sizeof(int) + (nbfix*10)*sizeof(int) + (ntfix*10)*sizeof(int);
}
static inline int read_num(std::ifstream &is,std::string s )
{
    std::string str;
    int res;
    getline(is,str);
    std::size_t found= str.find(s);
//    std::cout<<str<<std::endl;
    std::string str_ = "CANNOT FIND " + s +  " IN THIS STRING !! \n Fail read ani file ";
	M_Assert((std::string::npos != found), str_.c_str() );
    sscanf(str.c_str() + found + s.size(), "%d",&res);

    return res;
}
void ReadAniFromFile(std::string filename,Ani_Mesh &m ){
	const char * name = filename.c_str();
    int *iW,*iF;
    std::ifstream file;
    file.open(filename.c_str());
    M_Assert((file.is_open()),"cannto open file in read ani");

    m.nv = read_num(file,"T points:");
    m.nb = read_num(file,"T faces:");
    m.nt = read_num(file,"T elements:");
    m.n_curv = read_num(file,"T curved faces:");
    m.nvfix = read_num(file,"T fixed points:");
    m.nbfix = read_num(file,"T fixed faces:");
    m.ntfix = read_num(file,"T fixed elements:");

    std::cout<<"Mesh parameters "<<m.nv<<" "<<m.nb<<" "<<m.nt<<std::endl;
    //std::cout<<m.n_curv<<std::endl;
	M_Assert((m.n_curv == 0), "Cannot create mesh yet. Curved Faces are in this mesh");
    file.close();

    m.nvmax = m.nv*2;
    m.nbmax = m.nb*2;
    m.ntmax = m.nt*2;
    m.vrt = (double *) malloc( 3*m.nvmax*sizeof(double));
    m.tet = (int *) malloc( 4*m.ntmax*sizeof(int));
    m.bnd = (int *) malloc( 3*m.nbmax*sizeof(int));

    m.labelT = (int *) malloc( m.ntmax*sizeof(int));
    m.labelF = (int *) malloc( m.nbmax*sizeof(int));


    iW = (int *) malloc( 100*m.ntmax*sizeof(int));
    iF = (int *) malloc( 100*m.ntmax*sizeof(int));
	
	if(m.nvfix>0)
	    m.ivfix = (int *) malloc((m.nvfix*10)*sizeof(int));
	if(m.nbfix >0)
	    m.ibfix = (int *) malloc((m.nbfix*10)*sizeof(int));
	if(m.ntfix>0)
	    m.itfix = (int *) malloc((m.ntfix*10)*sizeof(int));

    loadmani_(
         &m.nvmax,&m.nbmax,&m.ntmax,
         &m.nv,&m.nb,&m.nt,
         m.vrt,m.bnd,m.tet,m.labelF, m.labelT,
         &m.nvfix,&m.nbfix,&m.ntfix,m.ivfix,m.ibfix, m.itfix,
         iW,iF, filename.c_str());

	free(iW);free(iF);
    std::cout<<"Reading ended"<<std::endl;
}	

void LoadAFT(std::string filename,INMOST::Mesh *m)
{
    int i,j;
    const char * name = filename.c_str();
    std::ifstream file;
    file.open(filename.c_str());
    M_Assert((file.is_open()),"cannto open file in read aft");
    Tag Label_tag = m->CreateTag(LabelTag, DATA_INTEGER, CELL | FACE | EDGE | NODE, CELL | FACE | EDGE | NODE, 1);
    Tag Fortran_Cell_Tag = m->CreateTag(FortranCellTag, DATA_INTEGER, CELL , CELL , 1);
    Tag Fortran_Node_Tag = m->CreateTag(FortranNodeTag, DATA_INTEGER, NODE , NODE , 1);

    int nv,nt,nf;
    file >> nv;
    ElementArray<Node> newverts(m);
    newverts.reserve(nv);
     // std::cout<<"mesh begin "<< nv<<std::endl;
    for (i = 0; i < nv; i++) {
        Storage::real xyz[3];
        //    std::cout<<i<<std::endl;
        file>>xyz[0];
        file>>xyz[1];
        file>>xyz[2];
//        xyz[0] = animesh.vrt[3 * i + 0];
//        xyz[1] = animesh.vrt[3 * i + 1];
//        xyz[2] = animesh.vrt[3 * i + 2];
        //   std::cout<<i<<std::endl;
        Node nod =m->CreateNode(xyz);
        nod->Integer(Fortran_Node_Tag) = i;
        newverts.push_back(nod);
    }

    file >> nt;
  //   std::cout<<"mesh tet "<< nt<<std::endl;

 //   ElementArray<Cell> newtets(m);
   // newtets.reserve(nt);
    int num_node_tmp;
    for (i = 0; i < nt; i++) {
        ElementArray<Node> verts(m);
        for (j = 0; j < 4; j++) {
            file >> num_node_tmp;
            verts.push_back(newverts[num_node_tmp - 1]);
        }

        const INMOST_DATA_INTEGER_TYPE ne_face_nodes[12] = {0, 1, 2, 2, 1, 3, 1, 0, 3, 3, 0, 2};
        const INMOST_DATA_INTEGER_TYPE ne_num_nodes[4] = {3, 3, 3, 3};

        std::pair<Cell, bool> pair = m->CreateCell(verts, ne_face_nodes, ne_num_nodes,
                                                     4);
        file>>pair.first.Integer(Label_tag) ;
        pair.first.Integer(Fortran_Cell_Tag) = i;
        //std::cout<<i<<std::endl;
        //newtets.push_back(pair.first);
    }
    for (Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); it++) {
        // it->Integer(BC_tag) = No_BC_mark;
        it->Integer(Label_tag) = No_BC_mark;
    }

    file >> nf;
   // std::cout<<"mesh bnd "<< nf<<std::endl;
    int num_face[3];
    int labelF;
    int num_nodes[3];
    bool flag_find = true;
    for(i=0;i<nf;i++){
        for(j=0;j<3;j++) file>>num_face[j];
        file>>labelF;
        for (Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); it++) {
            for(j=0;j<3;j++) num_nodes[j] = it->getNodes()[j].Integer(Fortran_Node_Tag);
            flag_find = true;
            for(j=0;j<3;j++)
                if( (num_face[j] != num_nodes[0])&& (num_face[j] != num_nodes[1]) && (num_face[j] != num_nodes[2]) ){
                    flag_find = false;
                    break;
                }
            if(flag_find){
                it->Integer(Label_tag) = labelF;
                break;
            }
        }
    }
}

void CreateInmostMeshFromAni(INMOST::Mesh *m,Ani_Mesh &animesh) {
    int i, j, k;
    double center[3];


    ElementArray<Node> newverts(m);
    newverts.reserve(animesh.nv);
  //  std::cout<<"mesh begin "<< animesh.nv<<" max nv "<<animesh.nvmax<<std::endl;
    for (i = 0; i < animesh.nv; i++) {
        Storage::real xyz[3];
    //    std::cout<<i<<std::endl;
        xyz[0] = animesh.vrt[3 * i + 0];
        xyz[1] = animesh.vrt[3 * i + 1];
        xyz[2] = animesh.vrt[3 * i + 2];
     //   std::cout<<i<<std::endl;
        newverts.push_back(m->CreateNode(xyz));
    }

  //  std::cout<<"mesh tet "<< animesh.nt<<" max nv "<<animesh.ntmax<<std::endl;
    ElementArray<Cell> newtets(m);
    newtets.reserve(animesh.nt);
    for (i = 0; i < animesh.nt; i++) {
        ElementArray<Node> verts(m);
        for (j = 0; j < 4; j++) {
            verts.push_back(newverts[animesh.tet[4 * i + j] - 1]);
        }

        const INMOST_DATA_INTEGER_TYPE ne_face_nodes[12] = {0, 1, 2, 2, 1, 3, 1, 0, 3, 3, 0, 2};
        const INMOST_DATA_INTEGER_TYPE ne_num_nodes[4] = {3, 3, 3, 3};

        std::pair<Cell, bool> pair = m->CreateCell(verts, ne_face_nodes, ne_num_nodes,
                                                   4); // Create north-east prismatic cell
        newtets.push_back(pair.first);
    }


        //  std::cout<<"mesh created"<<std::endl;
    //Tag Material_tag = m->CreateTag(name_material_tag.c_str(), DATA_INTEGER, CELL, CELL, 1);
    Tag Label_tag = m->CreateTag(LabelTag, DATA_INTEGER, CELL | FACE | EDGE | NODE, CELL | FACE | EDGE | NODE, 1);
    Tag Fortran_Cell_Tag = m->CreateTag(FortranCellTag, DATA_INTEGER, CELL , CELL , 1);

    Tag Fortran_Node_Tag = m->CreateTag(FortranNodeTag, DATA_INTEGER, NODE , NODE , 1);
    for (i = 0; i < animesh.nt; i++) {
      //  newtets[i].Integer(Material_tag) = animesh.labelT[i];
        newtets[i].Integer(Label_tag) = animesh.labelT[i];
        newtets[i].Integer(Fortran_Cell_Tag) = i;
    }
    for (i = 0; i < animesh.nv; i++) {
        newverts[i].Integer(Fortran_Node_Tag) = i;
    }
  //  Tag BC_tag = m->CreateTag(name_bc_tag.c_str(), DATA_INTEGER, FACE, FACE, 1);
//    for(i=0;i<animesh.nb;i++){
  //     M_Assert((animesh.labelF[i] != No_BC_mark),"We have labelF woth NO_BC_MARK");
  //  }
    double **centroids;
    double dist;
    centroids = (double **) malloc(animesh.nb * sizeof(double *));
    for (i = 0; i < animesh.nb; i++) {
        centroids[i] = (double *) malloc(3 * sizeof(double));
        for (j = 0; j < 3; j++)
            centroids[i][j] = 0.0;
        for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
                centroids[i][k] += animesh.vrt[3 * (animesh.bnd[3 * i + j] - 1) + k];
        for (j = 0; j < 3; j++)
            centroids[i][j] /= 3.0;
    }
    int count_bc = 0;
    //std::cout<<"find bound faces "<<std::endl;
#ifdef STATISTICS
  //  std::cout<<"Sizeof centroids "<< animesh.nb *3 * sizeof(double)/(1048576)<<" MB , or "<<animesh.nb * 3 * sizeof(double)/(1024)<<" KB "<<std::endl;
#endif

    for (Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); it++) {
       // it->Integer(BC_tag) = No_BC_mark;
        it->Integer(Label_tag) = No_BC_mark;
    }

    for(Mesh::iteratorFace it = m->BeginFace();it != m->EndFace();it++)
        if(it->GetStatus() != Element::Ghost)
    {
        if( (it->Boundary()))
        {
            // std::cout<<it->LocalID()<<std::endl;
            bool flag_found = false;
            it->Centroid(center);
            for (i = 0; i < animesh.nb; i++) {
                dist = 0;
                for (j = 0; j < 3; j++)
                    dist += (center[j] - centroids[i][j]) * (center[j] - centroids[i][j]);
                dist = sqrt(dist);
                if (dist < EPS) {
                    // it->Integer(BC_tag) = animesh.labelF[i];
                    it->Integer(Label_tag) = animesh.labelF[i];
                    flag_found = true;
                    if (animesh.labelF[i] != No_BC_mark)
                        count_bc++;
                    break;
                }
            }
              if(!flag_found){
                 std::cout<<"DIDNT FIND BOUNDARY TAG!!! "<<it->DataLocalID()<<std::endl;
            // it->Integer(BC_tag) = No_BC_mark;
                  it->Integer(Label_tag) = No_BC_mark;

             }
        }
        else
        {
            // it->Integer(BC_tag) = No_BC_mark;
            it->Integer(Label_tag) = No_BC_mark;

        }
    }


    int control_bc=0;
    for(i=0;i<animesh.nb;i++)
        if(animesh.labelF[i] != No_BC_mark)
           control_bc++;
    //M_Warinig((control_bc == animesh.nb), " There are bnd with No_BC_Tag");
    if(control_bc != count_bc)
        std::cout<<"some faces are unmarked!! "<<control_bc<<" "<< count_bc<<std::endl;

//    std::cout<<"Processor rank "<< m->GetProcessorRank()<<" creating ended"<<std::endl;


    for(i=0;i<animesh.nb;i++)
        free(centroids[i]);
    free(centroids);
   // return m;


}


void CreateInmostMeshFromAni_with_interial_marks(INMOST::Mesh *m,Ani_Mesh &animesh) {
    int i, j, k;
    double center[3];


    ElementArray<Node> newverts(m);
    newverts.reserve(animesh.nv);
    //  std::cout<<"mesh begin "<< animesh.nv<<" max nv "<<animesh.nvmax<<std::endl;
    for (i = 0; i < animesh.nv; i++) {
        Storage::real xyz[3];
        //    std::cout<<i<<std::endl;
        xyz[0] = animesh.vrt[3 * i + 0];
        xyz[1] = animesh.vrt[3 * i + 1];
        xyz[2] = animesh.vrt[3 * i + 2];
        //   std::cout<<i<<std::endl;
        newverts.push_back(m->CreateNode(xyz));
    }

    //  std::cout<<"mesh tet "<< animesh.nt<<" max nv "<<animesh.ntmax<<std::endl;
    ElementArray<Cell> newtets(m);
    newtets.reserve(animesh.nt);
    for (i = 0; i < animesh.nt; i++) {
        ElementArray<Node> verts(m);
        for (j = 0; j < 4; j++) {
            verts.push_back(newverts[animesh.tet[4 * i + j] - 1]);
        }

        const INMOST_DATA_INTEGER_TYPE ne_face_nodes[12] = {0, 1, 2, 2, 1, 3, 1, 0, 3, 3, 0, 2};
        const INMOST_DATA_INTEGER_TYPE ne_num_nodes[4] = {3, 3, 3, 3};

        std::pair<Cell, bool> pair = m->CreateCell(verts, ne_face_nodes, ne_num_nodes,
                                                   4); // Create north-east prismatic cell
        newtets.push_back(pair.first);
    }


    //  std::cout<<"mesh created"<<std::endl;
    //Tag Material_tag = m->CreateTag(name_material_tag.c_str(), DATA_INTEGER, CELL, CELL, 1);
    Tag Label_tag = m->CreateTag(LabelTag, DATA_INTEGER, CELL | FACE | EDGE | NODE, CELL | FACE | EDGE | NODE, 1);
    Tag Fortran_Cell_Tag = m->CreateTag(FortranCellTag, DATA_INTEGER, CELL , CELL , 1);

    Tag Fortran_Node_Tag = m->CreateTag(FortranNodeTag, DATA_INTEGER, NODE , NODE , 1);
    for (i = 0; i < animesh.nt; i++) {
        //  newtets[i].Integer(Material_tag) = animesh.labelT[i];
        newtets[i].Integer(Label_tag) = animesh.labelT[i];
        newtets[i].Integer(Fortran_Cell_Tag) = i;
    }
    for (i = 0; i < animesh.nv; i++) {
        newverts[i].Integer(Fortran_Node_Tag) = i;
    }

//    Tag BC_tag = m->CreateTag(name_bc_tag.c_str(), DATA_INTEGER, FACE, FACE, 1);
    for(i=0;i<animesh.nb;i++){
        M_Assert((animesh.labelF[i] != No_BC_mark),"We have labelF woth NO_BC_MARK");
    }
    double **centroids;
    double dist;
    centroids = (double **) malloc(animesh.nb * sizeof(double *));
    for (i = 0; i < animesh.nb; i++) {
        centroids[i] = (double *) malloc(3 * sizeof(double));
        for (j = 0; j < 3; j++)
            centroids[i][j] = 0.0;
        for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
                centroids[i][k] += animesh.vrt[3 * (animesh.bnd[3 * i + j] - 1) + k];
        for (j = 0; j < 3; j++)
            centroids[i][j] /= 3.0;
    }
    int count_bc = 0;
    //std::cout<<"find bound faces "<<std::endl;
#ifdef STATISTICS
    //  std::cout<<"Sizeof centroids "<< animesh.nb *3 * sizeof(double)/(1048576)<<" MB , or "<<animesh.nb * 3 * sizeof(double)/(1024)<<" KB "<<std::endl;
#endif

    for (Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); it++) {
        // it->Integer(BC_tag) = No_BC_mark;
        it->Integer(Label_tag) = No_BC_mark;
    }

    for(Mesh::iteratorFace it = m->BeginFace();it != m->EndFace();it++) {
        //if( (it->Boundary()))
        if(it->GetStatus() != Element::Ghost)
        {
            // std::cout<<it->LocalID()<<std::endl;
            bool flag_found = false;
            it->Centroid(center);
            for (i = 0; i < animesh.nb; i++) {
                dist = 0;
                for (j = 0; j < 3; j++)
                    dist += (center[j] - centroids[i][j]) * (center[j] - centroids[i][j]);
                dist = sqrt(dist);
                if (dist < EPS) {
                    // it->Integer(BC_tag) = animesh.labelF[i];
                    it->Integer(Label_tag) = animesh.labelF[i];
                    flag_found = true;
                    if (animesh.labelF[i] != No_BC_mark)
                        count_bc++;
                    break;
                }
            }
        }

    }

/*
    for(i=0;i<animesh.nb;i++) {
        bool flag_found=false;
        for(Mesh::iteratorFace it = m->BeginFace();it != m->EndFace();it++) {
            it->Centroid(center);
            dist = 0;
            for(j=0;j<3;j++)
                dist += (center[j]- centroids[i][j])*(center[j]- centroids[i][j]);
            dist = sqrt(dist);
            if(dist < EPS){
                it->Integer(BC_tag) = animesh.labelF[i];
                it->Integer(Label_tag) = animesh.labelF[i];
                flag_found = true;
                if( animesh.labelF[i] != No_BC_mark)
                    count_bc++;
                break;
            }
        }
        M_Assert(flag_found,"DIDNT FIND MARKED FACE");
    }
*/
    int control_bc=0;
    for(i=0;i<animesh.nb;i++)
        if(animesh.labelF[i] != No_BC_mark)
            control_bc++;
    //M_Warinig((control_bc == animesh.nb), " There are bnd with No_BC_Tag");
    if(control_bc != count_bc)
        std::cout<<"some faces are unmarked!! "<<control_bc<<" "<< count_bc<<std::endl;

//    std::cout<<"Processor rank "<< m->GetProcessorRank()<<" creating ended"<<std::endl;


    for(i=0;i<animesh.nb;i++)
        free(centroids[i]);
    free(centroids);
    // return m;


}


void CreateAniMeshFromInmost(Mesh *m,  Ani_Mesh &animesh){
    int i,j,k;

//    Tag BC_tag_init = m->GetTag(name_bc_tag);
//    Tag Material_tag_init = m->GetTag(name_material_tag);
//    M_Assert(BC_tag_init.isValid(),"We have no valid BC tag");
//    M_Assert(Material_tag_init.isValid(),"We have no valid Material tag");

    Tag Label_tag = m->GetTag(LabelTag);
    M_Assert(Label_tag.isValid(),"We have no valid label tag");
    //std::cout<<animesh.nvmax<<" "<<m->NumberOfNodes()<<std::endl;
    M_Assert(animesh.nvmax>= m->NumberOfNodes(),"nvmax is small");
    M_Assert(animesh.nbmax>= m->NumberOfFaces(),"nbmax is small");
    M_Assert(animesh.ntmax>= m->NumberOfCells(),"ntmax is small");


    if(animesh.vrt != NULL) free(animesh.vrt);
    animesh.vrt = (double *) malloc( 3*animesh.nvmax*sizeof(double));

    if(animesh.tet != NULL) free(animesh.tet);
    animesh.tet = (int *) malloc( 4*animesh.ntmax*sizeof(int));

    if(animesh.bnd != NULL) free(animesh.bnd);
    animesh.bnd = (int *) malloc( 3*animesh.nbmax*sizeof(int));

    if(animesh.labelT != NULL) free(animesh.labelT);
    animesh.labelT = (int *) malloc( animesh.ntmax*sizeof(int));

    if(animesh.labelF != NULL) free(animesh.labelF);
    animesh.labelF = (int *) malloc( animesh.nbmax*sizeof(int));

    std::map<int,int> map_nodes;
    std::map<int,int>::iterator map_numbers_iterator;


    animesh.nt=0;
    animesh.nb=0;
    for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++)
        if((it->GetStatus() != Element::Ghost)){

            //find nodes. If not - mark
            ElementArray<Node> nodes = it->getNodes();
            for(i=0;i<nodes.size();i++){

                int id = nodes[i].DataLocalID();
                map_numbers_iterator = map_nodes.find(id);

                if(map_numbers_iterator != map_nodes.end()){
                    animesh.tet[animesh.nt*4+i] =map_numbers_iterator->second+1;
                }
                else{
                    k=map_nodes.size();
                    map_nodes[id] = k;
                    for(int i1= 0;i1<3;i1++)
                        animesh.vrt[3*k+i1] = nodes[i].Coords()[i1];
                    animesh.tet[animesh.nt*4+i] = k+1;
                }
            }
            animesh.labelT[animesh.nt] = it->Integer(Label_tag);
            ElementArray<Face> faces = it->getFaces();
            for(i=0;i<faces.size();i++){
                if((faces[i].Boundary()) ||
                        ((faces[i].BackCell().GetStatus() == Element::Ghost)
                         | (faces[i].FrontCell().GetStatus() == Element::Ghost))){

                        ElementArray<Node> nodes_f = faces[i].getNodes();
                        for(j=0;j<nodes_f.size();j++){
                        map_numbers_iterator = map_nodes.find(nodes_f[j].DataLocalID());
                        M_Assert((map_numbers_iterator != map_nodes.end()),"DIDNT FIND NODE!!! ");
                        animesh.bnd[animesh.nb*3 +j] = map_numbers_iterator->second+1;
                    }
                    if(faces[i].Boundary())
                        animesh.labelF[animesh.nb] =faces[i].Integer(Label_tag);
                    animesh.nb++;
                }

            }

            animesh.nt++;
        }
    //std::cout<<"Proc Rank "<<m->GetProcessorRank()<<" number of tets "<<animesh.nt<<std::endl;
    animesh.nv = map_nodes.size();

}

void RefineLocalMesh(Mesh *m,Mesh *m_result,   int NumberOfRefines){
    int i,j,k;
    int *iW ;
    int FlagWrite = 0;
    Ani_Mesh animesh;

//    Tag BC_tag_init = m->GetTag(FaceTag);
 //   M_Assert(BC_tag_init.isValid(),"We have no valid BC tag");

  //  Tag Material_tag_init = m->GetTag(CellTag);
  //  M_Assert(Material_tag_init.isValid(),"We have no valid Material tag");

    animesh.nvmax = m->NumberOfNodes()+10;
    animesh.nbmax = m->NumberOfFaces()+10;
    animesh.ntmax = m->NumberOfCells()+10;

    for(i=0;i<NumberOfRefines;i++){
        animesh.nvmax *=8;
        animesh.ntmax *=8;
        animesh.nbmax *=8;
    }
    CreateAniMeshFromInmost(m,  animesh);
//    std::cout<<"create animesh ended"<<std::endl;
    int MaxWi= 100*animesh.ntmax;
    double *MapMtr;
    int *Ref2MapMtr;

    MapMtr = (double *) malloc( 3*3*animesh.ntmax*sizeof(double));
    Ref2MapMtr = (int *) malloc( 10*animesh.ntmax*sizeof(int));
    iW = (int *) malloc( MaxWi*sizeof(int));
//    std::cout<<"MaxWI "<<MaxWi<<" "<<animesh.ntmax<<std::endl;


    refine_(&NumberOfRefines,&animesh.nv,&animesh.nvmax,&animesh.nt,&animesh.ntmax,&animesh.nb,&animesh.nbmax,
            animesh.vrt,animesh.tet,animesh.bnd,animesh.labelT,animesh.labelF,MapMtr,Ref2MapMtr,iW,&MaxWi,&FlagWrite);
//    exit(1);//db!!
//    std::cout<<"refine ended"<<std::endl;
  //  std::cout<<animesh.vrt[0]<<" "<<animesh.vrt[1]<<" "<<animesh.vrt[2]<<std::endl;
    CreateInmostMeshFromAni(m_result,animesh);
  //  m_result->ResolveShared();
  //  m_result->ExchangeGhost(1,NODE); // Construct Ghost cells in 1 layers connected via nodes
#ifdef STATISTICS
  //  std::cout<<"Sizeof animesh  after refine "<< animesh.getsize()/(1048576)<<" MB , or "<< animesh.getsize()/(1024)<<" KB "<<std::endl;
 //   std::cout<<"number of nodes "<<animesh.nv<<" "<<animesh.nvmax<<std::endl;
#endif
    free(MapMtr);
    free(Ref2MapMtr);
}

void ReadInmostFromAniFile(INMOST::Mesh *m,std::string filename ){
    Ani_Mesh animesh ;
    ReadAniFromFile(filename,animesh );
    CreateInmostMeshFromAni(m,animesh);
    //m->ResolveShared();
#ifdef STATISTICS
  //  std::cout<<"Sizeof animesh "<< animesh.getsize()/(1048576)<<" MB , or "<< animesh.getsize()/(1024)<<" KB "<<std::endl;
 //   std::cout<<"number of nodes "<<animesh.nv<<" "<<animesh.nvmax<<std::endl;
#endif

  //  m->ExchangeGhost(1,NODE); // Construct Ghost cells in 1 layers connected via nodes
}

