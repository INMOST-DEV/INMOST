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

void liste2r_(int *nP, int *nR, int *nE, int *IPE, int *IRE, int *nEP, int *IEP);

void liste2f_(int *nP, int *nF, int *nE, int *IPE, int *IFE, int *nEP, int *IEP);
void order2d_(int *n_t, int *nE, int *IPE);
//order2D_(&n_t,&animesh.nt,animesh.tet);
//      Subroutine listE2F(nP, nF, nE, IPE, IFE, nEP, IEP)

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

 //   std::cout<<"bef"<<std::endl;
    int n_t=4;
    order2d_(&n_t,&m.nt,m.tet);
   // std::cout<<"aft"<<std::endl;
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
    Tag Label_tag = m->CreateTag(LabelTag, DATA_INTEGER, CELL | FACE | EDGE | NODE, NONE, 1);
    Tag Fortran_Num_Tag = m->CreateTag(FortranNumTag, DATA_INTEGER, CELL | FACE | EDGE | NODE, NONE, 1);

 //   Tag Fortran_Cell_Tag = m->CreateTag(DebugCellTag, DATA_INTEGER, CELL , CELL , 1);
 //   Tag Fortran_Node_Tag = m->CreateTag(DebugNodeTag, DATA_INTEGER, NODE , NODE , 1);
 //   Tag Fortran_Edge_Tag = m->CreateTag(DebugEdgeTag, DATA_INTEGER, EDGE , EDGE , 1);
  //  Tag Fortran_Face_Tag = m->CreateTag(DebugFaceTag, DATA_INTEGER, EDGE , EDGE , 1);

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
        nod->Integer(Fortran_Num_Tag) = i;
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

        std::pair<Cell, bool> pair = m->CreateCell(verts, ne_face_nodes, ne_num_nodes,4);
        file>>pair.first.Integer(Label_tag) ;
        pair.first.Integer(Fortran_Num_Tag) = i;
        //std::cout<<i<<std::endl;
        //newtets.push_back(pair.first);
    }
    for (Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); it++) {
        // it->Integer(BC_tag) = No_BC_mark;
        it->Integer(Label_tag) = No_BC_mark;
        it->Integer(Fortran_Num_Tag) = -1;
    }

    for (Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); it++) {
        // it->Integer(BC_tag) = No_BC_mark;
        //it->Integer(Label_tag) = No_BC_mark;
        it->Integer(Fortran_Num_Tag) = -1;
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
            for(j=0;j<3;j++) num_nodes[j] = it->getNodes()[j].Integer(Fortran_Num_Tag);
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
    Tag Label_tag = m->CreateTag(LabelTag, DATA_INTEGER, CELL | FACE | EDGE | NODE,NONE, 1);
    Tag Fortran_Num_Tag = m->CreateTag(FortranNumTag, DATA_INTEGER, CELL | FACE | EDGE | NODE, NONE, 1);

   // std::cout<<"mesh begin "<< animesh.nv<<" max nv "<<animesh.nvmax<<std::endl;
    for (i = 0; i < animesh.nv; i++) {
        Storage::real xyz[3];
    //    std::cout<<i<<std::endl;
        xyz[0] = animesh.vrt[3 * i + 0];
        xyz[1] = animesh.vrt[3 * i + 1];
        xyz[2] = animesh.vrt[3 * i + 2];
     //   std::cout<<i<<std::endl;
        newverts.push_back(m->CreateNode(xyz));
        newverts[newverts.size()-1].Integer(Fortran_Num_Tag) = i;
    }

    //std::cout<<"mesh tet "<< animesh.nt<<" max nv "<<animesh.ntmax<<std::endl;
    ElementArray<Cell> newtets(m);
    newtets.reserve(animesh.nt);
    for (i = 0; i < animesh.nt; i++) {
//        std::cout<<i<<std::endl;
        ElementArray<Node> verts(m);
        for (j = 0; j < 4; j++) {
            verts.push_back(newverts[animesh.tet[4 * i + j] - 1]);
        }
        const INMOST_DATA_INTEGER_TYPE ne_face_nodes[12] = {0, 1, 2, 2, 1, 3, 1, 0, 3, 3, 0, 2};
        const INMOST_DATA_INTEGER_TYPE ne_num_nodes[4] = {3, 3, 3, 3};

        std::pair<Cell, bool> pair = m->CreateCell(verts, ne_face_nodes, ne_num_nodes,4);
        newtets.push_back(pair.first);
        newtets[newtets.size() - 1].Integer(Label_tag) = animesh.labelT[i];
        newtets[newtets.size() - 1].Integer(Fortran_Num_Tag) = i;
    }

    std::cout<<"numeration begin"<<std::endl;

    // std::cout<<animesh.nt<<std::endl;
    //std::cout<<m->NumberOfEdges()<<std::endl;
    int nEds;
    int nEP[animesh.nv];
   // int IRE[6*animesh.nt];
    int *IRE;
    IRE = new int[6*animesh.nt];
  //  std::cout<<"before"<<std::endl;
   // std::cout<<"bef"<<std::endl;


    for (Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); it++) {
        // it->Integer(BC_tag) = No_BC_mark;
        it->Integer(Fortran_Num_Tag)=-1;
    }
    for (Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); it++) {
        // it->Integer(BC_tag) = No_BC_mark;
        it->Integer(Fortran_Num_Tag)=-1;
    }


  //  int IEP[4*m->NumberOfEdges()];
    int *IEP;
    IEP = new int[4*m->NumberOfEdges()];

    //std::cout<<m->NumberOfEdges()<<std::endl;
    liste2r_(&animesh.nv,&nEds,&animesh.nt,animesh.tet,IRE,nEP,IEP);


    int nFs;
   // int nEP[animesh.nv];
    //int IFE[4*animesh.nt];
    int *IFE;
    IFE = new int[4*animesh.nt];
    //int IEP[4*m->NumberOfEdges()];
  //  Call order2D(4, nt,tet)
    int n_t= 4;
    int i1,i2,i3;
    liste2f_(&animesh.nv,&nFs,&animesh.nt,animesh.tet,IFE,nEP,IEP);

    std::cout<<"numeration"<<std::endl;
/*for(i=0;i<16;i++)
    std::cout<<IFE[i]-1<<" ";
std::cout<<std::endl;
for(i=0;i<16;i++)
    std::cout<<animesh.tet[i]-1<<" ";
std::cout<<std::endl;*/

    for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++){
        int num_tet = it->Integer(Fortran_Num_Tag);
        //edges
        ElementArray<Edge> edges = it->getEdges();
        int count_local =0;
        bool flag_find = false;




        for(i=0;i<3;i++)
            for(j=i+1;j<4;j++){
                for(k=0;k<edges.size();k++){
                    ElementArray<Node> nodes_f = edges[k].getNodes();//it->getEdges();
                    //std::cout<<nodes_f[0].Integer(Fortran_Num_Tag)<<" "<<animesh.tet[4*num_tet + i]-1<<std::endl;
                   // std::cout<<nodes_f[1].Integer(Fortran_Num_Tag)<<" "<<animesh.tet[4*num_tet + j]-1<<std::endl;
                   // std::cout<<std::endl;
                    if(((nodes_f[0].Integer(Fortran_Num_Tag) == animesh.tet[4*num_tet + i]-1)  && (nodes_f[1].Integer(Fortran_Num_Tag) == animesh.tet[4*num_tet + j]-1) ) ||
                            ((nodes_f[1].Integer(Fortran_Num_Tag) == animesh.tet[4*num_tet + i]-1)  && (nodes_f[0].Integer(Fortran_Num_Tag) == animesh.tet[4*num_tet + j]-1) )){
                        flag_find=true;

                        if(edges[k].Integer(Fortran_Num_Tag) != -1) {

                          if(edges[k].Integer(Fortran_Num_Tag)!= IRE[6*num_tet + count_local]-1){
                              std::cout<<"Already marked "<<edges[k].Integer(Fortran_Num_Tag)<<" "<<IRE[6*num_tet + count_local]-1<<  std::endl;  
                              exit(-1);                  
                           }  
                        }
                        else{
                          edges[k].Integer(Fortran_Num_Tag)= IRE[6*num_tet + count_local]-1;

                        }
                        count_local++;
                        break;
                    }
                }
                M_Assert(flag_find,"Didn't find edge num !!!");
                flag_find=false;

            }

        //faces
        ElementArray<Face> faces = it->getFaces();
        //count_local =0;
        flag_find = false;
        
       // i1=0;i2=1;i3=2;
        for(i=0;i<4;i++){
            i1= (0+i)%4;
            i2= (1+i)%4;
            i3= (2+i)%4;
           // std::cout<<animesh.tet[4*num_tet + i1]-1<<" "<<animesh.tet[4*num_tet + i2]-1<<" "<<animesh.tet[4*num_tet + i3]-1<<std::endl;

            for(int k1=0;k1<faces.size();k1++){
                ElementArray<Node> nodes_f = faces[k1].getNodes();//it->getEdges();
                int num_fort_loc[3]; for(int j1=0;j1<3;j1++) num_fort_loc[j1] = nodes_f[j1].Integer(Fortran_Num_Tag);
                if(comp_numb(animesh.tet[4*num_tet + i1]-1,num_fort_loc)&& comp_numb(animesh.tet[4*num_tet + i2]-1,num_fort_loc)&& comp_numb(animesh.tet[4*num_tet + i3]-1,num_fort_loc)){
              //      std::cout<<faces[k1].Integer(Fortran_Num_Tag)<<" "<<faces[k1].DataLocalID()<<std::endl;
                    flag_find=true;
                    if(faces[k1].Integer(Fortran_Num_Tag) !=-1){
                        if(faces[k1].Integer(Fortran_Num_Tag) != IFE[4*num_tet + i]-1){
                            std::cout<<"already marked "<<std::endl; 
                            std::cout<<"tet "<<num_tet<<std::endl;
                            std::cout<<i1<<" "<<i2<<" "<<i3<<std::endl;
                            std::cout<<i<<std::endl;
                            std::cout<<faces[k1].Integer(Fortran_Num_Tag) <<" "<<IFE[4*num_tet + i]-1<<std::endl;
                            std::cout<<animesh.tet[4*num_tet + i1]-1<<" "<<animesh.tet[4*num_tet + i2]-1<<" "<<animesh.tet[4*num_tet + i3]-1<<std::endl;
                            std::cout<<num_fort_loc[0]<<" "<<num_fort_loc[1]<<" "<<num_fort_loc[2]<<std::endl;
                            std::cout<<IFE[4*num_tet + 0]-1<<" "<<IFE[4*num_tet + 1]-1<<" "<<IFE[4*num_tet + 2]-1<<" "<<IFE[4*num_tet + 3]-1<<std::endl;
                            exit(-1);
                        }
                   }
                   else{
                       faces[k1].Integer(Fortran_Num_Tag) = IFE[4*num_tet + i]-1;
                     //  std::cout<<faces[k1].Integer(Fortran_Num_Tag)<<" "<<faces[k1].DataLocalID()<<std::endl;
                     //  std::cout<<faces[k1].Integer(Fortran_Num_Tag)<<std::endl;
                   }
                   break;

               }
           }
            if(! flag_find){
                std::cout<<i1<<" "<<i2<<" "<<i3<<std::endl;
            }
           M_Assert(flag_find,"Didn't find face num !!!");
           flag_find=false;
      }
     // std::cout<<std::endl;
    }
          std::cout<<"mesh created"<<std::endl;
    //Tag Material_tag = m->CreateTag(name_material_tag.c_str(), DATA_INTEGER, CELL, CELL, 1);
    //Tag Fortran_Cell_Tag = m->CreateTag(FortranCellTag, DATA_INTEGER, CELL , CELL , 1);

   // Tag Fortran_Node_Tag = m->CreateTag(FortranNodeTag, DATA_INTEGER, NODE , NODE , 1);

  //  Tag BC_tag = m->CreateTag(name_bc_tag.c_str(), DATA_INTEGER, FACE, FACE, 1);
   // for(i=0;i<animesh.nb;i++){
   //    M_Assert((animesh.labelF[i] != No_BC_mark),"We have labelF woth NO_BC_MARK ");
    //}


    for (Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); it++) {
        // it->Integer(BC_tag) = No_BC_mark;
        it->Integer(Label_tag) = No_BC_mark;
    }


    for(i=0;i<animesh.nb;i++){
        ElementArray<Node> verts(m);
        for (j = 0; j < 3; j++) {
            verts.push_back(newverts[animesh.bnd[3 * i + j] - 1]);
        }
        std::pair<Face, bool> pair = m->CreateFace(verts);
        M_Assert((!pair.second), "we didn't find face!!!");
        pair.first.Integer(Label_tag) = animesh.labelF[i];

    }


    /*
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
    std::cout<<"find bound faces "<<std::endl;
#ifdef STATISTICS
    std::cout<<"Sizeof centroids "<< animesh.nb *3 * sizeof(double)/(1048576)<<" MB , or "<<animesh.nb * 3 * sizeof(double)/(1024)<<" KB "<<std::endl;
#endif


    for(Mesh::iteratorFace it = m->BeginFace();it != m->EndFace();it++)
    if(it->GetStatus() != Element::Ghost)
    {
#ifdef NO_INTERNAL_MARKS
        if( (it->Boundary()))
#endif
        {
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
#ifdef NO_INTERNAL_MARKS
            if(!flag_found){
                 std::cout<<"DIDNT FIND BOUNDARY TAG!!! "<<it->DataLocalID()<<std::endl;
            // it->Integer(BC_tag) = No_BC_mark;
                  it->Integer(Label_tag) = No_BC_mark;

             }
#endif
        }
#ifdef NO_INTERNAL_MARKS
        else
        {
            // it->Integer(BC_tag) = No_BC_mark;
            it->Integer(Label_tag) = No_BC_mark;

        }
#endif
    }


*/
//           control_bc++;


    //
//    if(control_bc != count_bc)
//        std::cout<<"some faces are unmarked!! "<<control_bc<<" "<< count_bc<<std::endl;

//    std::cout<<"Processor rank "<< m->GetProcessorRank()<<" creating ended"<<std::endl;

//////////////////// Debug edge and face tags/////////
    delete[] IFE;
    delete[] IEP;
    delete [] IRE;
 //   for(i=0;i<animesh.nb;i++)
 //       free(centroids[i]);
  //  free(centroids);
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
    std::cout<<"refine ended"<<std::endl;
    std::cout<<animesh.vrt[0]<<" "<<animesh.vrt[1]<<" "<<animesh.vrt[2]<<std::endl;
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

