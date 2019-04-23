#include "inmost_ani_fem.h"
//
// Created by vasily on 3/17/17.
//


//#define DEBUG_FORTRAN_NUMERATION
//#define DEBUG_NUMBERS
/*struct Ani_discretization{
    Ani_discretization();
    int NumDof[4];
    int NumElem[4];
    int MinElem[4];
    int MaxElem[4];
    int InitIndex[4];
    Tag NumTags[4];
    Mesh *m;
    PrepareProblemAni(Mesh *m);
};*/
#define MIN_EPS 1e-16
extern "C" {

void fem3dext_(
        double * XY1,double *  XY2,double *  XY3,double *  XY4,
        int * lbE,int *  lbF,int *  lbR,int *  lbP,double *  DDATA,int *  IDATA,int *  iSYS,
        int * LDA,double *  A,double *  F,int *  nRow,int *  nCol,
        int * templateR,int *  templateC);


}


Ani_discretization::Ani_discretization(Mesh *m_){
    m= m_;
    MinInd = 0;
    MaxInd = 0;
    MatrSize=0;
    for(int i=0;i<4;i++){
        NumDof[i]=0;
        NumElem[i]=0;
        MinElem[i]=0;
        MaxElem[i]=0;
        InitIndex[i]=0;
    }
}


void Ani_discretization::PrepareProblemAni(int numNodeDofs, int numEdgeDofs, int numFaceDofs, int numCellDofs, ASSEMBLING_TYPE type_a_, ORDER_TYPE type_o_, std::string prefix)
{
    int i,j;
    int tmp_var;
    type_assembling = type_a_;
    type_order = type_o_;

    NumDof[0] = numNodeDofs;
    NumDof[1] = numEdgeDofs;
    NumDof[2] = numFaceDofs;
    NumDof[3] = numCellDofs;

    ////creating tag names
    for(i=0;i<4;i++){
        NumTagNames[i] = prefix+ NumTagNames_consts[i];
        FakeNumTagNames[i] = prefix+ FakeNumTagNames_consts[i];
    }
    ///Some tests
    tmp_var = 0;
    for(i=0;i<4;i++)
        tmp_var+= NumDof[i];
    M_Assert((tmp_var != 0), "Number DOFs is 0");
    for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++)
        M_Assert((it->GetGeometricType() == Element::Tet),"Mesh is not tetrahedral");

    ///creating mask for numerating all nodes, cells,faces, edges, if need
    int MASK = 0;
    for(i=0;i<4;i++)
        if(NumDof[i] >0)
            MASK = MASK | (NODE << i);
    MASK = MASK | CELL;
#ifndef DEBUG_FORTRAN_NUMERATION
    m->AssignGlobalID(MASK);
#else
    M_Assert((this->m->GetProcessorsNumber() ==1 ),"this test can be only on one proc");
    Tag Fortran_num_tag = m->GetTag(FortranNumTag);
    M_Assert(Fortran_num_tag.isValid(), "We have no valid Fortran_num_tag tag in test");

#endif
    tmp_var = m->TotalNumberOf(NODE|EDGE|FACE|CELL)+1;
    for(i=0;i<4;i++){
        MinElem[i] = tmp_var;
        MaxElem[i] = -1;    
    }
    for(Mesh::iteratorElement it = m->BeginElement(MASK); it != m->EndElement();it++){
        j=FastInnerIndex[it->GetElementType()/2];
#ifndef DEBUG_FORTRAN_NUMERATION
        tmp_var = it->GlobalID();
#else
        tmp_var = it->Integer(Fortran_num_tag);
#endif
        if(it->GetStatus() != Element::Ghost){
            NumElem[j]++;
            if(tmp_var<MinElem[j])
                MinElem[j] = tmp_var;
            if(tmp_var>MaxElem[j])
                MaxElem[j] = tmp_var;
        }
    }
   // for(i=0;i<4;i++)
  //      std::cout<<"Min Elem "<<MinElem[i]<<" max elem "<<MaxElem[i]<<std::endl;
    for(i=0;i<4;i++)
        if(NumDof[i] > 0)
            M_Assert((NumElem[i] != 0),"Zero number of element");
    for(i=0;i<4;i++)
        MinInd += NumDof[i]*MinElem[i];
    MaxInd = MinInd;
    for(i=0;i<4;i++)
        MaxInd += NumDof[i]*NumElem[i];
    MaxInd--;

    for(i=0;i<4;i++){
        if(NumDof[i]>0) {
            M_Assert((!m->HaveTag(NumTagNames[i])), "mesh already have a tag");
//            std::cout << i <<" "<<(ElementType )(NODE<<i)<<" "<<NODE<< std::endl;
            NumTags[i] = m->CreateTag(NumTagNames[i], DATA_INTEGER, (ElementType )(NODE << i), (ElementType )(NODE << i), NumDof[i]);
        }
    }


    for(i=0;i<4;i++)
        if(NumDof[i]>0) {
            M_Assert(NumTags[i].isValid(), "NumTag is not valid!!");
        }


    MatrSize = m->TotalNumberOf(NODE)*NumDof[0] + m->TotalNumberOf(FACE)*NumDof[2] + m->TotalNumberOf(EDGE)*NumDof[1] + m->TotalNumberOf(CELL)*NumDof[3];

    switch (type_order){
        case STRAIGHT:
            InitIndex[0] = MinInd;
            for(i=1;i<4;i++)
                InitIndex[i] = InitIndex[i-1] + NumDof[i-1]*NumElem[i-1];
            break;
        case REVERSE:
            InitIndex[3] = MinInd;
            for(i=2;i>=0;i--)
                InitIndex[i] = InitIndex[i+1] + NumDof[i+1]*NumElem[i+1];
            break;
        default:
            std::cout<<"Wrong type order"<<std::endl;
            abort();
    }


    for(Mesh::iteratorElement it = m->BeginElement(MASK); it != m->EndElement();it++) {
        j=FastInnerIndex[it->GetElementType()/2];
        if (it->GetStatus() != Element::Ghost) {
            for (i = 0; i < NumDof[j]; i++) {
                int num_inner= -1;
#ifndef DEBUG_FORTRAN_NUMERATION
                num_inner = it->GlobalID();
#else
                num_inner = it->Integer(Fortran_num_tag);
                M_Assert(!(num_inner == -1), "Failed to set num");
#endif


                switch (type_assembling) {
                    case MINIBLOCKS:
                        it->IntegerArray(NumTags[j])[i] =
                                (num_inner - MinElem[j]) * NumDof[j] + InitIndex[j] + i;
                        break;
                    case ANITYPE:
                        it->IntegerArray(NumTags[j])[i] =
                                (num_inner - MinElem[j]) + InitIndex[j] + i * NumElem[j];
                        break;
                    case OSEEN_FORTRA_STOKES:
                        if(j==0){
                            it->IntegerArray(NumTags[j])[i] =
                                    (num_inner - MinElem[j])  + i * (NumElem[0]+NumElem[1]);

                        }
                        else if(j==1){
                            it->IntegerArray(NumTags[j])[i] =
                                    (num_inner - MinElem[j])  + i * (NumElem[0]+NumElem[1]) +NumElem[0];
                        }
                        else{
                            std::cout<<"Diff j "<<j<<std::endl;
                            exit(0);
                        }
                        break;
                    default:
                        std::cerr << "Wrong type of assembling" << std::endl;
                        abort();
                }
                //int tst = it->IntegerArray(NumTags[j])[i];
                if((it->IntegerArray(NumTags[j])[i] > MaxInd))
                    std::cout<<MaxInd<<" "<<it->IntegerArray(NumTags[j])[i]<<std::endl;
                M_Assert((it->IntegerArray(NumTags[j])[i] <= MaxInd), "Wrong index min");

                M_Assert((it->IntegerArray(NumTags[j])[i] >= MinInd), "Wrong index max");

            }
        }
    }
    for(i=0;i<4;i++)
        if(NumDof[i]>0) m->ExchangeData(NumTags[i],(ElementType )(NODE<<i),0);


    //see template.f and FDofOrient of this tag
    FaceBackCells = m->CreateTag("FaceBackCellsForOrient", DATA_INTEGER,FACE,NONE,1);
    std::cout<<"orient faces"<<std::endl;
    for (Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); it++)
    if(it->GetStatus() != Element::Ghost){
#ifndef  DEBUG_FORTRAN_NUMERATION
        if(!it->Boundary()) {

            M_Assert((it->FrontCell().isValid()), "frontcell is not valid");
            M_Assert((it->BackCell().isValid()), "Backcell is not valid");
            if (it->BackCell().GlobalID() > it->FrontCell().GlobalID())
                it->Integer(FaceBackCells) = it->BackCell().GlobalID();
            else
                it->Integer(FaceBackCells) = it->FrontCell().GlobalID();
        }
        else{
            M_Assert((it->BackCell().isValid()), "Backcell bnd is not valid");
            it->Integer(FaceBackCells) = it->BackCell().GlobalID();
        }

#else
        if(!it->Boundary()) {
            if (it->BackCell().Integer(Fortran_num_tag) > it->FrontCell().Integer(Fortran_num_tag))
                it->Integer(FaceBackCells) = it->BackCell().Integer(Fortran_num_tag);
            else
                it->Integer(FaceBackCells) = it->FrontCell().Integer(Fortran_num_tag);
        }
        else{
            it->Integer(FaceBackCells) = it->BackCell().Integer(Fortran_num_tag);
        }
#endif
    }
    m->ExchangeData(FaceBackCells,FACE,0);
    std::cout<<"end prepare"<<std::endl;
}




void Ani_discretization::Assemble (Tag Label_Tag,Sparse::Matrix &matrix, Sparse::Vector &rhs, double *DDATA, int DDATA_size, int *IDATA, int IDATA_size, void (*fem3dext_new_)(
        double * XY1,double *  XY2,double *  XY3,double *  XY4,
        int * lbE,int *  lbF,int *  lbR,int *  lbP,double *  DDATA,int *  IDATA,int *  iSYS,
        int * LDA,double *  A,double *  F,int *  nRow,int *  nCol,
        int * templateR,int *  templateC)) {
//    void Ani_discretization::Assemble (Tag Label_Tag,Tag Cell_Tag, Tag Face_Tag,Tag Edge_Tag, Tag Node_Tag,Sparse::Matrix &matrix, Sparse::Vector &rhs, double *DDATA, int DDATA_size, int *IDATA, int IDATA_size) {

    int i, j, k;
    double XY[4][3];
    int lbE, lbF[4], lbR[6], lbP[4];
    int nRows;
    nRows = 4 * NumDof[0] + 6 * NumDof[1] + 4 * NumDof[2] + NumDof[3];

    //std::cout<<"Nrows "<<nRows<<std::endl;
    double Max_A, Min_A;
    Max_A= -1;
    Min_A = 1000000000000;

    double A[nRows * nRows];
    double F[nRows];
    int templateR[nRows], templateC[nRows];
    int templateNumR[nRows], templateNumC[nRows];
    int templateNumRDof[nRows], templateNumCDof[nRows];
    int indexesR[nRows], indexesC[nRows];
    bool index_ghostC[nRows];
    int iSYS[21]; ///system array for fortran
    int local_edge_index[6];
    int local_face_index[4];
    int count_DOFs[4];
    int count_elems[4];
    bool flag_symmetric;

    matrix.SetInterval(MinInd,MaxInd+1);
    rhs.SetInterval(MinInd,MaxInd+1);
   // std::cout<<"Nrows "<<nRows<<std::endl;

#ifndef DEBUG_FORTRAN_NUMERATION

#else

    M_Assert((this->m->GetProcessorsNumber() ==1 ),"this test can be only on one proc");
    Tag Fortran_num_tag = m->GetTag(FortranNumTag);
    M_Assert(Fortran_num_tag.isValid(), "We have no valid Fortran_num_tag tag in test");
#endif


   // std::cout<<"Nrows "<<nRows<<std::endl;

    for(i=MinInd;i<MaxInd+1;i++)
        matrix[i][i]=0.0;
   // std::cout<<"Nrows "<<nRows<<std::endl;
    for (Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++) {
    //    std::cout<<it->DataLocalID()<<std::endl;
        ElementArray<Node> nodes = it->getNodes();
        ElementArray<Edge> edges = it->getEdges();
        ElementArray<Face> faces = it->getFaces();
        int cnt = 0;
        //find local order fo edges
        for (i = 0; i < 3; i++) {
            for (j = i + 1; j < 4; j++) {
                bool flag_control = false;
#ifndef DEBUG_FORTRAN_NUMERATION
                int ind_loc[2] = {nodes[i].DataLocalID(), nodes[j].DataLocalID()};
#else
                int ind_loc[2] = {nodes[i].Integer(Fortran_num_tag), nodes[j].Integer(Fortran_num_tag)};

#endif
                for (k = 0; k < 6; k++) {
#ifndef DEBUG_FORTRAN_NUMERATION
                    int ind_loc_ed[2] = {edges[k].getNodes()[0].DataLocalID(), edges[k].getNodes()[1].DataLocalID()};
#else
                    int ind_loc_ed[2] = {edges[k].getNodes()[0].Integer(Fortran_num_tag), edges[k].getNodes()[1].Integer(Fortran_num_tag)};

#endif
                    if (((ind_loc[0] == ind_loc_ed[0]) && (ind_loc[1] == ind_loc_ed[1])) ||
                        ((ind_loc[1] == ind_loc_ed[0]) && (ind_loc[0] == ind_loc_ed[1]))) {
                        local_edge_index[cnt] = k;
                        cnt++;
                        flag_control = true;
                        break;
                    }
                }
                M_Assert(flag_control, "didnt find edge ");
            }
        }
        ///find local order of faces
        for (i = 0; i < 4; i++) {
            bool flag_control = false;
#ifndef DEBUG_FORTRAN_NUMERATION
            int ind_loc[3] = {nodes[i].DataLocalID(),nodes[(i + 1) % 4].DataLocalID(),nodes[(i + 2) % 4].DataLocalID() };
#else
            int ind_loc[3] = {nodes[i].Integer(Fortran_num_tag),nodes[(i + 1) % 4].Integer(Fortran_num_tag),nodes[(i + 2) % 4].Integer(Fortran_num_tag) };
#endif
            for (j = 0; j < 4; j++) {
#ifndef DEBUG_FORTRAN_NUMERATION
                int ind_loc2[3] = {faces[j].getNodes()[0].DataLocalID(),faces[j].getNodes()[1].DataLocalID(),faces[j].getNodes()[2].DataLocalID()};
#else
                int ind_loc2[3] = {faces[j].getNodes()[0].Integer(Fortran_num_tag),faces[j].getNodes()[1].Integer(Fortran_num_tag),faces[j].getNodes()[2].Integer(Fortran_num_tag)};
#endif
                if( comp_numb(ind_loc[0],ind_loc2) && comp_numb(ind_loc[1],ind_loc2)&& comp_numb(ind_loc[2],ind_loc2)){
                    local_face_index[i] = j;
                    flag_control = true;
                    break;
                }
            }
            M_Assert(flag_control, "didnt find face ");
        }

        for(i=0;i<4;i++)
            for(j=0;j<3;j++){
                XY[i][j] = nodes[i].Coords()[j];
            }
#ifndef DEBUG_FORTRAN_NUMERATION
        iSYS[2] = it->DataLocalID()+1;
        for(i=0;i<4;i++) iSYS[3+i] = nodes[i].DataLocalID()+1;
        for(i=0;i<6;i++) iSYS[7+i] = edges[local_edge_index[i]].DataLocalID()+1;
        for(i=0;i<4;i++) iSYS[13+i] = faces[local_face_index[i]].DataLocalID()+1;
#else
        iSYS[2] = it->Integer(Fortran_num_tag)+1;
        for(i=0;i<4;i++) iSYS[3+i] = nodes[i].Integer(Fortran_num_tag)+1;
        for(i=0;i<6;i++) iSYS[7+i] = edges[local_edge_index[i]].Integer(Fortran_num_tag)+1;
        for(i=0;i<4;i++) iSYS[13+i] = faces[local_face_index[i]].Integer(Fortran_num_tag)+1;
#endif




        iSYS[17] = m->NumberOfNodes();
        iSYS[18] = m->NumberOfEdges();
        iSYS[19] = m->NumberOfFaces();
        iSYS[20] = m->NumberOfCells();

        M_Assert(Label_Tag.isValid(),"Cell_tag is not valid");
        lbE = it->Integer(Label_Tag);
        for(i=0;i<4;i++) lbF[i] = faces[local_face_index[i]].Integer(Label_Tag);
        for(i=0;i<6;i++) lbR[i] = edges[local_edge_index[i]].Integer(Label_Tag);
        for(i=0;i<4;i++) lbP[i] = nodes[i].Integer(Label_Tag);

        //std::cout<<it->DataLocalID()<<std::endl;

        if(fem3dext_new_== nullptr) {
            fem3dext_(
                    XY[0], XY[1], XY[2], XY[3],
                    &lbE, lbF, lbR, lbP, DDATA, IDATA, iSYS,
                    &nRows, A, F, &nRows, &nRows,
                    templateR, templateC);
        }
        else{
            fem3dext_new_(
                    XY[0], XY[1], XY[2], XY[3],
                    &lbE, lbF, lbR, lbP, DDATA, IDATA, iSYS,
                    &nRows, A, F, &nRows, &nRows,
                    templateR, templateC);
        }
        //std::cout<<it->DataLocalID()<<std::endl;

        flag_symmetric = true;


        for(i=0;i<nRows;i++){
            if(( templateC[i] != templateR[i])  || (templateC[i] == 5) ||(templateC[i] == 6) )
            {
                flag_symmetric = false;
                break;
            }
        }
    //    for(i=0;i<nRows;i++){
     //       if((templateC[i]== 5 ) && (nodes[count_elems[Ndof]].GetStatus() != Element::Ghost))
     //   }
        for(i=0;i<4;i++){
            count_DOFs[i]=0;
            count_elems[i] = 0;
        }
        for(i=0;i<nRows;i++){
         //   std::cout<<templateC[i]<<std::endl;
            switch (templateC[i] - 1){
                case 0: //Nodes
                    M_Assert((count_DOFs[Ndof]<=NumDof[Ndof] ), "Wrong number of Node DOF!!!!!");
                    indexesC[i] = nodes[count_elems[Ndof]].IntegerArray(NumTags[Ndof])[count_DOFs[Ndof]];
                    index_ghostC[i] =( nodes[count_elems[Ndof]].GetStatus() == Element::Ghost);
                    count_elems[Ndof]++;
                    if(count_elems[Ndof] == 4){
                        count_elems[Ndof]=0;
                        count_DOFs[Ndof]++;
                    }
                    break;
                case 1:
                    M_Assert((count_DOFs[Fdof]<=NumDof[Fdof] ), "Wrong number of Face DOF!!!!!");
                    indexesC[i] =faces[local_face_index[count_elems[Fdof] ] ].IntegerArray(NumTags[Fdof])[ count_DOFs[Fdof]];
                    index_ghostC[i] =(faces[local_face_index[count_elems[Fdof] ] ].GetStatus() == Element::Ghost) ;
                    count_elems[Fdof]++;
                    if(count_elems[Fdof] == 4){
                        count_elems[Fdof]=0;
                        count_DOFs[Fdof]++;
                    }
                    break;
                case 2:
                    M_Assert((count_DOFs[Edof]<=NumDof[Edof] ), "Wrong number of Edge DOF!!!!!");
                    indexesC[i] =edges[local_edge_index[count_elems[Edof] ]].IntegerArray(NumTags[Edof])[count_DOFs[Edof]];
                    index_ghostC[i] =(edges[local_edge_index[count_elems[Edof] ] ].GetStatus() == Element::Ghost) ;
                    count_elems[Edof]++;
                    if(count_elems[Edof] == 6){
                        count_elems[Edof]=0;
                        count_DOFs[Edof]++;
                    }

                    break;
                case 3:
                    M_Assert((count_DOFs[Cdof]<=NumDof[Cdof] ), "Wrong number of Cell DOF!!!!!");
                    index_ghostC[i] = (it->GetStatus() == Element::Ghost);
                    indexesC[i] = it->IntegerArray(NumTags[Cdof])[count_DOFs[Cdof]];
                        count_DOFs[Cdof]++;
                    break;
                case 5:
                    M_Assert((count_DOFs[Fdof]<=NumDof[Fdof] ), "Wrong number of Face DOF!!!!!");
                    indexesC[i] =faces[local_face_index[count_elems[Fdof] ] ].IntegerArray(NumTags[Fdof])[ count_DOFs[Fdof]];
                    index_ghostC[i] =(faces[local_face_index[count_elems[Fdof] ] ].GetStatus() == Element::Ghost) ;
                    count_elems[Fdof]++;
                    if(count_elems[Fdof] == 4){
                        count_elems[Fdof]=0;
                        count_DOFs[Fdof]++;
                    }
                    //if(!index_ghostC[i] )
                    {
#ifndef DEBUG_FORTRAN_NUMERATION
                        if(it->DataLocalID() != faces[local_face_index[count_elems[Fdof] ] ].Integer(FaceBackCells)){
#else
                        if(it->Integer(Fortran_num_tag) != faces[local_face_index[count_elems[Fdof] ] ].Integer(FaceBackCells)){
#endif
                            for(j=0;j<nRows;j++)
                                A[i+j*nRows] = -A[i+j*nRows];
                            F[i] = -F[i];
                        }
                    }
                    break;
                case 6:
                    M_Assert((count_DOFs[Edof]<=NumDof[Edof] ), "Wrong number of Edge DOF!!!!!");
                    indexesC[i] =edges[local_edge_index[count_elems[Edof] ]].IntegerArray(NumTags[Edof])[count_DOFs[Edof]];
                    index_ghostC[i] =(edges[local_edge_index[count_elems[Edof] ] ].GetStatus() == Element::Ghost) ;
                    //if(!index_ghostC[i] )
                    { //if number of nodes[i].GlobalID()< nodes[j].GlobalID()
                        int cnt_tmp=0;
                        int i1,j1;
                        for(i1=0;i1<3;i1++)
                            for(j1=i1+1;j1<4;j1++){
                                if(cnt_tmp == count_elems[Edof]) {
                                    break;
                                }
                                cnt_tmp++;
                            }
                        M_Assert((cnt_tmp<5)," Fail to find edge in RdofOrient");
#ifndef DEBUG_FORTRAN_NUMERATION
                        if(nodes[i1].GlobalID()<nodes[j1].GlobalID()){
#else
                        if(nodes[i1].Integer(Fortran_num_tag)<nodes[j1].Integer(Fortran_num_tag)){
#endif
                            for(j=0;j<nRows;j++)
                                A[i+j*nRows] = -A[i+j*nRows];
                            F[i] = -F[i];
                        }
                    }


                    count_elems[Edof]++;
                    if(count_elems[Edof] == 6){
                        count_elems[Edof]=0;
                        count_DOFs[Edof]++;
                    }
                    break;
                default:
                    std::cout<<"wrong template "<<templateC[i] - 1<<std::endl;
                    abort();
            }
        }


        for(i=0;i<4;i++){
            count_DOFs[i]=0;
            count_elems[i] = 0;
        }
        if(flag_symmetric){
            for(i=0;i<nRows;i++)
                indexesR[i]= indexesC[i];
        }
        else {
            for(i=0;i<nRows;i++){
                switch (templateR[i] - 1){
                    case 0: //Nodes
                        M_Assert((count_DOFs[Ndof]<=NumDof[Ndof] ), "Wrong number of Node DOF!!!!!");

                        indexesR[i] = nodes[count_elems[Ndof]].IntegerArray(NumTags[Ndof])[count_DOFs[Ndof]];//nodes[templateNumR[i]-1].IntegerArray(NumTags[Ndof_fortran])[templateNumRDof[i]-1];
                        count_elems[Ndof]++;
                        if(count_elems[Ndof] == 4){
                            count_elems[Ndof]=0;
                            count_DOFs[Ndof]++;
                        }
                        break;
                    case 1:
                        M_Assert((count_DOFs[Fdof]<=NumDof[Fdof] ), "Wrong number of Face DOF!!!!!");

                        indexesR[i] = faces[local_face_index[count_elems[Fdof] ] ].IntegerArray(NumTags[Fdof])[ count_DOFs[Fdof]];;//faces[local_face_index[templateNumR[i]-1]].IntegerArray(NumTags[Fdof_fortran])[templateNumRDof[i]-1];
                        count_elems[Fdof]++;
                        if(count_elems[Fdof] == 4){
                            count_elems[Fdof]=0;
                            count_DOFs[Fdof]++;
                        }
                        break;
                    case 2:
                        M_Assert((count_DOFs[Edof]<=NumDof[Edof] ), "Wrong number of Edge DOF!!!!!");
                        indexesR[i] =  edges[local_edge_index[count_elems[Edof] ]].IntegerArray(NumTags[Edof])[count_DOFs[Edof]];
                        count_elems[Edof]++;
                        if(count_elems[Edof] == 6){
                            count_elems[Edof]=0;
                            count_DOFs[Edof]++;
                        }
                        break;
                    case 3:
                        M_Assert((count_DOFs[Cdof]<=NumDof[Cdof] ), "Wrong number of Cell DOF!!!!!");
                        indexesR[i] =it->IntegerArray(NumTags[Cdof])[count_DOFs[Cdof]];
                        count_DOFs[Cdof]++;


                    case 5:
                        M_Assert((count_DOFs[Fdof]<=NumDof[Fdof] ), "Wrong number of Face DOF!!!!!");
                        indexesR[i] =faces[local_face_index[count_elems[Fdof] ] ].IntegerArray(NumTags[Fdof])[ count_DOFs[Fdof]];
                        count_elems[Fdof]++;
                        if(count_elems[Fdof] == 4){
                            count_elems[Fdof]=0;
                            count_DOFs[Fdof]++;
                        }
                        //if(!index_ghostC[i] )
                        {
#ifndef DEBUG_FORTRAN_NUMERATION
                            if(it->DataLocalID() != faces[local_face_index[count_elems[Fdof] ] ].Integer(FaceBackCells)){
#else
                            if(it->Integer(Fortran_num_tag) != faces[local_face_index[count_elems[Fdof] ] ].Integer(FaceBackCells)){
#endif

                                for(j=0;j<nRows;j++)
                                    A[j+i*nRows] = -A[j+i*nRows];

                            }
                        }
                        break;
                    case 6:
                        M_Assert((count_DOFs[Edof]<=NumDof[Edof] ), "Wrong number of Edge DOF!!!!!");
                        indexesR[i] =edges[local_edge_index[count_elems[Edof] ]].IntegerArray(NumTags[Edof])[count_DOFs[Edof]];
                        //if(!index_ghostC[i] )
                        { //if number of nodes[i].GlobalID()< nodes[j].GlobalID()
                            int cnt_tmp=0;
                            int i1,j1;
                            for(i1=0;i1<3;i1++)
                                for(j1=i1+1;j1<4;j1++){
                                    if(cnt_tmp == count_elems[Edof]) {
                                        break;
                                    }
                                    cnt_tmp++;
                                }

                            M_Assert((cnt_tmp<5)," Fail to find edge in RdofOrient");
#ifndef DEBUG_FORTRAN_NUMERATION
                            if(nodes[i1].GlobalID()<nodes[j1].GlobalID()){
#else
                            if(nodes[i1].Integer(Fortran_num_tag)<nodes[j1].Integer(Fortran_num_tag)){
#endif
                                for(j=0;j<nRows;j++)
                                    A[j+i*nRows] = -A[j+i*nRows];
                            }
                        }


                        count_elems[Edof]++;
                        if(count_elems[Edof] == 6){
                            count_elems[Edof]=0;
                            count_DOFs[Edof]++;
                        }
                        break;
                    default:
                        std::cout<<"wrong template "<<templateR[i] - 1<<std::endl;
                        abort();
                }
            }
        }

        for(i=0;i<nRows;i++){
            if(!index_ghostC[i]){//cols
                if((indexesC[i]<MinInd) ||(indexesC[i]>MaxInd)){
                    std::cout<<"wrong index C "<<indexesC[i]<<" "<<MinInd<<" "<<MaxInd<<" "<<i<<std::endl;
                    abort();
                }
                for(j=0;j<nRows;j++){ //rows
                    {

                        if((indexesR[i]<0) ||(indexesR[i]>MatrSize)){
                            std::cout<<"Proc rank "<<m->GetProcessorRank()<<std::endl;
                            std::cout<<"wrong index R "<<indexesR[i]<<" "<<0<<" "<<MatrSize<<std::endl;
                            std::cout<<"flag_symmetric "<<flag_symmetric<<std::endl;
                            std::cout<<"Num elems"<<std::endl;
                            for(int i1=0;i1<4;i1++)
                                std::cout<<NumElem[i1]<<" ";
                            std::cout<<std::endl;
                            std::cout<<"Min elems"<<std::endl;
                            for(int i1=0;i1<4;i1++)
                                std::cout<<MinElem[i1]<<" ";
                            std::cout<<"Max elems"<<std::endl;
                            for(int i1=0;i1<4;i1++)
                                std::cout<<MaxElem[i1]<<" ";
                            std::cout<<std::endl;
                            std::cout<<"Min ind "<< MinInd<<" Max ind "<<MaxInd<< std::endl;
                            std::cout<<"indexes"<<std::endl;
                            for(int i1=0;i1<nRows;i1++)
                                std::cout<<indexesR[i1]<<" ";
                            std::cout<<std::endl;
                            abort();
                        }
                        if(A[j*nRows +i] != 0.0)
                            matrix[indexesC[i]][indexesR[j]] +=A[j*nRows +i];

                        if(std::isnan(A[j*nRows +i])){
                            std::cout<<"not a number in matrix"<<std::endl;
                            std::cout<<it->DataLocalID()<< " "<< A[j*nRows +i]<<std::endl;
                            for(int i1=0;i1<nRows;i1++) {
                                for (int j1 = 0; j1 < nRows; j1++) {
                                    std::cout << A[j1 * nRows + i1] << " ";
                                }
                                std::cout<<std::endl;
                            }
                            std::cout<<"num tet " <<iSYS[2]-1<<std::endl;
                            abort();
                        }
                        if(fabs(A[j*nRows +i])>Max_A)
                            Max_A = fabs(A[j*nRows +i]);
                        if((fabs(A[j*nRows +i])<Min_A) && (A[i*nRows +j] != 0.0))
                            Min_A = fabs(A[j*nRows +i]);
                    }
                }
                rhs[indexesC[i]] +=F[i];
            }

        }

        /////////////////////////////
    }
    //std::cout<<"Min elem "<<Min_A<<" Max elem "<<Max_A<<std::endl;
}






//====================================================================================leading dimension test====================================
/*

void Ani_discretization::PrepareProblemAni_LD( int numNodeDofs, int numEdgeDofs, int numFaceDofs, int numCellDofs, ASSEMBLING_TYPE type_a_, ORDER_TYPE type_o_)
{
    int i,j;
    int tmp_var;
    type_assembling = type_a_;
    type_order = type_o_;

    NumDof[0] = numNodeDofs;
    NumDof[1] = numEdgeDofs;
    NumDof[2] = numFaceDofs;
    NumDof[3] = numCellDofs;
    ///Some tests
    tmp_var = 0;
    for(i=0;i<4;i++)
        tmp_var+= NumDof[i];
    M_Assert((tmp_var != 0), "Number DOFs is 0");
    for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++)
        M_Assert((it->GetGeometricType() == Element::Tet),"Mesh is not tetrahedral");

    ///find leadig dimension
    Leading_DOF = NumDof[0];
    for(i=1;i<4;i++)
        if(NumDof[i]>Leading_DOF)
            Leading_DOF=NumDof[i];

    ///creating mask for numerating all nodes, cells,faces, edges, if need
    int MASK = 0;
    for(i=0;i<4;i++)
        if(NumDof[i] >0)
            MASK = MASK | (NODE << i);
    m->AssignGlobalID(MASK);

    tmp_var = m->TotalNumberOf(NODE|EDGE|FACE|CELL)+1;
    for(i=0;i<4;i++){
        MinElem[i] = tmp_var;
        MaxElem[i] = -1;
    }
    for(Mesh::iteratorElement it = m->BeginElement(MASK); it != m->EndElement();it++){
        j=FastInnerIndex[it->GetElementType()/2];
        tmp_var = it->GlobalID();
        if(it->GetStatus() != Element::Ghost){
            NumElem[j]++;
            if(tmp_var<MinElem[j])
                MinElem[j] = tmp_var;
            if(tmp_var>MaxElem[j])
                MaxElem[j] = tmp_var;
        }
    }

    for(i=0;i<4;i++)
        if(NumDof[i] > 0)
            M_Assert((NumElem[i] != 0),"Zero number of element");
    for(i=0;i<4;i++)
        if(NumDof[i]>0)
            MinInd += Leading_DOF*MinElem[i];
    MaxInd = MinInd;
    for(i=0;i<4;i++)
        if(NumDof[i]>0)
            MaxInd += Leading_DOF*NumElem[i];
    MaxInd--;

    for(i=0;i<4;i++){
        if(NumDof[i]>0) {
//            std::cout << i <<" "<<(ElementType )(NODE<<i)<<" "<<NODE<< std::endl;
            NumTags[i] = m->CreateTag(NumTagNames[i], DATA_INTEGER, (ElementType )(NODE << i), (ElementType )(NODE << i), NumDof[i]);
        }
    }
    for(i=0;i<4;i++){
        if(NumDof[i]>0) {
            if(NumDof[i]<Leading_DOF){
            std::cout << i <<" "<<Leading_DOF - NumDof[i]<< std::endl;
            FakeNumTags[i] = m->CreateTag(FakeNumTagNames[i], DATA_INTEGER, (ElementType )(NODE << i), (ElementType )(NODE << i), Leading_DOF - NumDof[i]);

            }
        }
    }


    for(i=0;i<4;i++)
        if(NumDof[i]>0) {
            M_Assert(NumTags[i].isValid(), "NumTag is not valid!!");
        }
    for(i=0;i<4;i++){
        if(NumDof[i]>0) {
            if(NumDof[i]<Leading_DOF)
//            std::cout << i <<" "<<(ElementType )(NODE<<i)<<" "<<NODE<< std::endl;
                M_Assert(FakeNumTags[i].isValid(), "NumTag is not valid!!");
//                FakeNumTags[i] = m->CreateTag(NumTagNames[i], DATA_INTEGER, (ElementType )(NODE << i), (ElementType )(NODE << i), Leading_DOF - NumDof[i]);
        }
    }
        MatrSize = 0;
    for(i=0;i<4;i++)
        if(NumDof[i]>0)
            MatrSize+= m->TotalNumberOf((ElementType )(NODE << i))* Leading_DOF;
//    MatrSize = m->TotalNumberOf(NODE)*NumDof[0] + m->TotalNumberOf(FACE)*NumDof[2] + m->TotalNumberOf(EDGE)*NumDof[1] + m->TotalNumberOf(CELL)*NumDof[3];


    switch (type_order){
        case STRAIGHT:
            InitIndex[0] = MinInd;
            for(i=1;i<4;i++)
                InitIndex[i] = InitIndex[i-1] + Leading_DOF*NumElem[i-1];
            break;
        case REVERSE:
            InitIndex[3] = MinInd;
            for(i=2;i>=0;i--)
                InitIndex[i] = InitIndex[i+1] + Leading_DOF*NumElem[i+1];
            break;
        default:
            std::cout<<"Wrong type order"<<std::endl;
            abort();
    }

    M_Warinig( (type_assembling == MINIBLOCKS),"Strange type assembling in leading dimension alg");
    for(Mesh::iteratorElement it = m->BeginElement(MASK); it != m->EndElement();it++) {
        j=FastInnerIndex[it->GetElementType()/2];
        if (it->GetStatus() != Element::Ghost) {
            //std::cout<<"!!!!!"<<std::endl;
            for (i = 0; i < NumDof[j]; i++) {
                switch (type_assembling) {
                    case MINIBLOCKS:
                        it->IntegerArray(NumTags[j])[i] =
                                (it->GlobalID() - MinElem[j]) * Leading_DOF + InitIndex[j] + i;
                        break;
                    case ANITYPE:
                        it->IntegerArray(NumTags[j])[i] =
                                (it->GlobalID() - MinElem[j]) + InitIndex[j] + i * NumElem[j];
                        break;
                    default:
                        std::cerr << "Wrong type of assembling" << std::endl;
                        abort();
                }
                int tst = it->IntegerArray(NumTags[j])[i];
                M_Assert((it->IntegerArray(NumTags[j])[i] <= MaxInd), "Wrong index");
                M_Assert((it->IntegerArray(NumTags[j])[i] >= MinInd), "Wrong index");
              //  std::cout<<it->IntegerArray(NumTags[j])[i]<<std::endl;

            }
            if(NumDof[j]<Leading_DOF){
               // std::cout<<"test "<<j<<" "<<FakeNumTags[j].isValid()<<" "<< FakeNumTags[j].GetSize()<< std::endl;
                for(i=0;i<Leading_DOF-NumDof[j];i++){
//                    std::cout<<i<<" "<<m->NumberOfNodes()<<std::endl;
                    //std::cout<<"!!!!!!! " <<(it->GlobalID() - MinElem[j])<<" "<<InitIndex[j]<< " "<<i+NumDof[j]<<" "<< (it->GlobalID() - MinElem[j]) * Leading_DOF + InitIndex[j] + i+NumDof[j]<< std::endl;
                    it->IntegerArray(FakeNumTags[j])[i] = (it->GlobalID() - MinElem[j]) * Leading_DOF + InitIndex[j] + i+NumDof[j];
                    M_Assert((it->IntegerArray(FakeNumTags[j])[i] <= MaxInd), "Wrong index");
                    M_Assert((it->IntegerArray(FakeNumTags[j])[i] >= MinInd), "Wrong index");

                }
            }
            //std::cout<<"test of one block of numbers "<<std::endl;
            if((j==1)&&(it->GlobalID()==0)){
                std::cout<<"test of one block of numbers "<<std::endl;
                for (i = 0; i < NumDof[j]; i++)
                         std::cout<<   it->IntegerArray(NumTags[j])[i] <<" ";

                if(NumDof[j]<Leading_DOF){
                    for(i=0;i<Leading_DOF-NumDof[j];i++)
                        std::cout<<it->IntegerArray(FakeNumTags[j])[i] <<" ";
                }
                std::cout<<" !!!!"<<std::endl;
            }

        }
    }

    for(i=0;i<4;i++)
        if(NumDof[i]>0) m->ExchangeData(NumTags[i],(ElementType )(NODE<<i),0);
    for(i=0;i<4;i++)
        if(NumDof[i]>0)
            if(NumDof[i]<Leading_DOF) m->ExchangeData(FakeNumTags[i],(ElementType )(NODE<<i),0);


    for(i=0;i<4;i++)
        if(this->FakeNumTags[i].isValid())
            std::cout<<" face tags "<<i<<" "<<Leading_DOF<<" "<<NumDof[i]<<std::endl;


    //see template.f and FDofOrient fof this tag
    FaceBackCells = m->CreateTag("FaceBackCellsForOrient", DATA_INTEGER,FACE,FACE,1);
    for (Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); it++)
        if(it->GetStatus() != Element::Ghost){
            if(!it->Boundary()) {
                if (it->BackCell().GlobalID() > it->FrontCell().GlobalID())
                    it->Integer(FaceBackCells) = it->BackCell().GlobalID();
                else
                    it->Integer(FaceBackCells) = it->FrontCell().GlobalID();
            }
            else{
                it->Integer(FaceBackCells) = it->BackCell().GlobalID();
            }
        }
    m->ExchangeData(FaceBackCells,FACE,0);
}

*/
/*
void Ani_discretization::Assemble_LD (Tag Label_Tag,Sparse::Matrix &matrix, Sparse::Vector &rhs, double *DDATA, int DDATA_size, int *IDATA, int IDATA_size) {
//    void Ani_discretization::Assemble (Tag Label_Tag,Tag Cell_Tag, Tag Face_Tag,Tag Edge_Tag, Tag Node_Tag,Sparse::Matrix &matrix, Sparse::Vector &rhs, double *DDATA, int DDATA_size, int *IDATA, int IDATA_size) {

    int i, j, k;
    double XY[4][3];
    int lbE, lbF[4], lbR[6], lbP[4];
    int nRows;
    nRows = 4 * NumDof[0] + 6 * NumDof[1] + 4 * NumDof[2] + NumDof[3];

    double Max_A, Min_A;
    Max_A= -1;
    Min_A = 1000000000000;

    double A[nRows * nRows];
    double F[nRows];
    int templateR[nRows], templateC[nRows];
    int templateNumR[nRows], templateNumC[nRows];
    int templateNumRDof[nRows], templateNumCDof[nRows];
    int indexesR[nRows], indexesC[nRows];
    bool index_ghostC[nRows];
    int iSYS[21]; ///system array for fortran
    int local_edge_index[6];
    int local_face_index[4];
    int count_DOFs[4];
    int count_elems[4];
    bool flag_symmetric;

    for(i=0;i<4;i++)
        if(this->FakeNumTags[i].isValid())
            std::cout<<" face tags "<<i<<std::endl;

            matrix.SetInterval(MinInd,MaxInd+1);
    rhs.SetInterval(MinInd,MaxInd+1);

    for(i=MinInd;i<MaxInd+1;i++)
        matrix[i][i]=0.0;
    for (Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++) {
        //std::cout<<it->DataLocalID()<<std::endl;
        ElementArray<Node> nodes = it->getNodes();
        ElementArray<Edge> edges = it->getEdges();
        ElementArray<Face> faces = it->getFaces();
        int cnt = 0;
        //find local order fo edges
        for (i = 0; i < 3; i++) {
            for (j = i + 1; j < 4; j++) {
                bool flag_control = false;
                int ind_loc[2] = {nodes[i].DataLocalID(), nodes[j].DataLocalID()};
                for (k = 0; k < 6; k++) {
                    int ind_loc_ed[2] = {edges[k].getNodes()[0].DataLocalID(), edges[k].getNodes()[1].DataLocalID()};
                    if (((ind_loc[0] == ind_loc_ed[0]) && (ind_loc[1] == ind_loc_ed[1])) ||
                        ((ind_loc[1] == ind_loc_ed[0]) && (ind_loc[0] == ind_loc_ed[1]))) {
                        local_edge_index[cnt] = k;
                        cnt++;
                        flag_control = true;
                        break;
                    }
                }
                M_Assert(flag_control, "didnt find edge ");
            }
        }
        ///find local order of faces
        for (i = 0; i < 4; i++) {
            bool flag_control = false;
            int ind_loc[3] = {nodes[i].DataLocalID(),nodes[(i + 1) % 4].DataLocalID(),nodes[(i + 2) % 4].DataLocalID() };
            for (j = 0; j < 4; j++) {
                int ind_loc2[3] = {faces[j].getNodes()[0].DataLocalID(),faces[j].getNodes()[1].DataLocalID(),faces[j].getNodes()[2].DataLocalID()};
                if( comp_numb(ind_loc[0],ind_loc2) && comp_numb(ind_loc[1],ind_loc2)&& comp_numb(ind_loc[2],ind_loc2)){
                    local_face_index[i] = j;
                    flag_control = true;
                    break;
                }
            }
            M_Assert(flag_control, "didnt find face ");
        }

        for(i=0;i<4;i++)
            for(j=0;j<3;j++){
                XY[i][j] = nodes[i].Coords()[j];
            }

        iSYS[2] = it->DataLocalID()+1;
        for(i=0;i<4;i++) iSYS[3+i] = nodes[i].DataLocalID()+1;
        for(i=0;i<6;i++) iSYS[7+i] = edges[local_edge_index[i]].DataLocalID()+1;
        for(i=0;i<4;i++) iSYS[13+i] = faces[local_face_index[i]].DataLocalID()+1;

        iSYS[17] = m->NumberOfNodes();
        iSYS[18] = m->NumberOfEdges();
        iSYS[19] = m->NumberOfFaces();
        iSYS[20] = m->NumberOfCells();

        M_Assert(Label_Tag.isValid(),"Cell_tag is not valid");
        lbE = it->Integer(Label_Tag);
        for(i=0;i<4;i++) lbF[i] = faces[local_face_index[i]].Integer(Label_Tag);
        for(i=0;i<6;i++) lbR[i] = edges[local_edge_index[i]].Integer(Label_Tag);
        for(i=0;i<4;i++) lbP[i] = nodes[i].Integer(Label_Tag);

        //std::cout<<it->DataLocalID()<<std::endl;

        fem3dext_(
                XY[0],XY[1],XY[2],XY[3],
                &lbE,lbF,lbR,lbP,DDATA, IDATA, iSYS,
                &nRows,A,F,&nRows,&nRows,
                templateR,templateC);
        //std::cout<<it->DataLocalID()<<std::endl;

        flag_symmetric = true;


        for(i=0;i<nRows;i++){
            if(( templateC[i] != templateR[i])  || (templateC[i] == 5) ||(templateC[i] == 6) )
            {
                flag_symmetric = false;
                break;
            }
        }
        //    for(i=0;i<nRows;i++){
        //       if((templateC[i]== 5 ) && (nodes[count_elems[Ndof]].GetStatus() != Element::Ghost))
        //   }
        for(i=0;i<4;i++){
            count_DOFs[i]=0;
            count_elems[i] = 0;
        }
        for(i=0;i<nRows;i++){
            //   std::cout<<templateC[i]<<std::endl;
            switch (templateC[i] - 1){
                case 0: //Nodes
                    M_Assert((count_DOFs[Ndof]<=NumDof[Ndof] ), "Wrong number of Node DOF!!!!!");
                    indexesC[i] = nodes[count_elems[Ndof]].IntegerArray(NumTags[Ndof])[count_DOFs[Ndof]];
                    index_ghostC[i] =( nodes[count_elems[Ndof]].GetStatus() == Element::Ghost);
                    count_elems[Ndof]++;
                    if(count_elems[Ndof] == 4){
                        count_elems[Ndof]=0;
                        count_DOFs[Ndof]++;
                    }
                    break;
                case 1:
                    M_Assert((count_DOFs[Fdof]<=NumDof[Fdof] ), "Wrong number of Face DOF!!!!!");
                    indexesC[i] =faces[local_face_index[count_elems[Fdof] ] ].IntegerArray(NumTags[Fdof])[ count_DOFs[Fdof]];
                    index_ghostC[i] =(faces[local_face_index[count_elems[Fdof] ] ].GetStatus() == Element::Ghost) ;
                    count_elems[Fdof]++;
                    if(count_elems[Fdof] == 4){
                        count_elems[Fdof]=0;
                        count_DOFs[Fdof]++;
                    }
                    break;
                case 2:
                    M_Assert((count_DOFs[Edof]<=NumDof[Edof] ), "Wrong number of Edge DOF!!!!!");
                    indexesC[i] =edges[local_edge_index[count_elems[Edof] ]].IntegerArray(NumTags[Edof])[count_DOFs[Edof]];
                    index_ghostC[i] =(edges[local_edge_index[count_elems[Edof] ] ].GetStatus() == Element::Ghost) ;
                    count_elems[Edof]++;
                    if(count_elems[Edof] == 6){
                        count_elems[Edof]=0;
                        count_DOFs[Edof]++;
                    }

                    break;
                case 3:
                    M_Assert((count_DOFs[Cdof]<=NumDof[Cdof] ), "Wrong number of Cell DOF!!!!!");
                    index_ghostC[i] = (it->GetStatus() == Element::Ghost);
                    indexesC[i] = it->IntegerArray(NumTags[Cdof])[count_DOFs[Cdof]];
                    count_DOFs[Cdof]++;
                    break;
                case 5:
                    M_Assert((count_DOFs[Fdof]<=NumDof[Fdof] ), "Wrong number of Face DOF!!!!!");
                    indexesC[i] =faces[local_face_index[count_elems[Fdof] ] ].IntegerArray(NumTags[Fdof])[ count_DOFs[Fdof]];
                    index_ghostC[i] =(faces[local_face_index[count_elems[Fdof] ] ].GetStatus() == Element::Ghost) ;
                    count_elems[Fdof]++;
                    if(count_elems[Fdof] == 4){
                        count_elems[Fdof]=0;
                        count_DOFs[Fdof]++;
                    }
                    //if(!index_ghostC[i] )
                    {
                        if(it->DataLocalID() != faces[local_face_index[count_elems[Fdof] ] ].Integer(FaceBackCells)){
                            for(j=0;j<nRows;j++)
                                A[i+j*nRows] = -A[i+j*nRows];
                            F[i] = -F[i];
                        }
                    }
                    break;
                case 6:
                    M_Assert((count_DOFs[Edof]<=NumDof[Edof] ), "Wrong number of Edge DOF!!!!!");
                    indexesC[i] =edges[local_edge_index[count_elems[Edof] ]].IntegerArray(NumTags[Edof])[count_DOFs[Edof]];
                    index_ghostC[i] =(edges[local_edge_index[count_elems[Edof] ] ].GetStatus() == Element::Ghost) ;
                    //if(!index_ghostC[i] )
                    { //if number of nodes[i].GlobalID()< nodes[j].GlobalID()
                        int cnt_tmp=0;
                        int i1,j1;
                        for(i1=0;i1<3;i1++)
                            for(j1=i1+1;j1<4;j1++){
                                if(cnt_tmp == count_elems[Edof]) {
                                    break;
                                }
                                cnt_tmp++;
                            }
                        M_Assert((cnt_tmp<5)," Fail to find edge in RdofOrient");
                        if(nodes[i1].GlobalID()<nodes[j1].GlobalID()){
                            for(j=0;j<nRows;j++)
                                A[i+j*nRows] = -A[i+j*nRows];
                            F[i] = -F[i];
                        }
                    }


                    count_elems[Edof]++;
                    if(count_elems[Edof] == 6){
                        count_elems[Edof]=0;
                        count_DOFs[Edof]++;
                    }
                    break;
                default:
                    std::cout<<"wrong template "<<templateC[i] - 1<<std::endl;
                    abort();
            }
        }


        for(i=0;i<4;i++){
            count_DOFs[i]=0;
            count_elems[i] = 0;
        }
        if(flag_symmetric){
            for(i=0;i<nRows;i++)
                indexesR[i]= indexesC[i];
        }
        else {
            for(i=0;i<nRows;i++){
                switch (templateR[i] - 1){
                    case 0: //Nodes
                        M_Assert((count_DOFs[Ndof]<=NumDof[Ndof] ), "Wrong number of Node DOF!!!!!");

                        indexesR[i] = nodes[count_elems[Ndof]].IntegerArray(NumTags[Ndof])[count_DOFs[Ndof]];//nodes[templateNumR[i]-1].IntegerArray(NumTags[Ndof_fortran])[templateNumRDof[i]-1];
                        count_elems[Ndof]++;
                        if(count_elems[Ndof] == 4){
                            count_elems[Ndof]=0;
                            count_DOFs[Ndof]++;
                        }
                        break;
                    case 1:
                        M_Assert((count_DOFs[Fdof]<=NumDof[Fdof] ), "Wrong number of Face DOF!!!!!");

                        indexesR[i] = faces[local_face_index[count_elems[Fdof] ] ].IntegerArray(NumTags[Fdof])[ count_DOFs[Fdof]];;//faces[local_face_index[templateNumR[i]-1]].IntegerArray(NumTags[Fdof_fortran])[templateNumRDof[i]-1];
                        count_elems[Fdof]++;
                        if(count_elems[Fdof] == 4){
                            count_elems[Fdof]=0;
                            count_DOFs[Fdof]++;
                        }
                        break;
                    case 2:
                        M_Assert((count_DOFs[Edof]<=NumDof[Edof] ), "Wrong number of Edge DOF!!!!!");
                        indexesR[i] =  edges[local_edge_index[count_elems[Edof] ]].IntegerArray(NumTags[Edof])[count_DOFs[Edof]];
                        count_elems[Edof]++;
                        if(count_elems[Edof] == 6){
                            count_elems[Edof]=0;
                            count_DOFs[Edof]++;
                        }
                        break;
                    case 3:
                        M_Assert((count_DOFs[Cdof]<=NumDof[Cdof] ), "Wrong number of Cell DOF!!!!!");
                        indexesR[i] =it->IntegerArray(NumTags[Cdof])[count_DOFs[Cdof]];
                        count_DOFs[Cdof]++;


                    case 5:
                        M_Assert((count_DOFs[Fdof]<=NumDof[Fdof] ), "Wrong number of Face DOF!!!!!");
                        indexesR[i] =faces[local_face_index[count_elems[Fdof] ] ].IntegerArray(NumTags[Fdof])[ count_DOFs[Fdof]];
                        count_elems[Fdof]++;
                        if(count_elems[Fdof] == 4){
                            count_elems[Fdof]=0;
                            count_DOFs[Fdof]++;
                        }
                        //if(!index_ghostC[i] )
                        {
                            if(it->DataLocalID() != faces[local_face_index[count_elems[Fdof] ] ].Integer(FaceBackCells)){
                                for(j=0;j<nRows;j++)
                                    A[j+i*nRows] = -A[j+i*nRows];

                            }
                        }
                        break;
                    case 6:
                        M_Assert((count_DOFs[Edof]<=NumDof[Edof] ), "Wrong number of Edge DOF!!!!!");
                        indexesR[i] =edges[local_edge_index[count_elems[Edof] ]].IntegerArray(NumTags[Edof])[count_DOFs[Edof]];
                        //if(!index_ghostC[i] )
                        { //if number of nodes[i].GlobalID()< nodes[j].GlobalID()
                            int cnt_tmp=0;
                            int i1,j1;
                            for(i1=0;i1<3;i1++)
                                for(j1=i1+1;j1<4;j1++){
                                    if(cnt_tmp == count_elems[Edof]) {
                                        break;
                                    }
                                    cnt_tmp++;
                                }
                            M_Assert((cnt_tmp<5)," Fail to find edge in RdofOrient");
                            if(nodes[i1].GlobalID()<nodes[j1].GlobalID()){
                                for(j=0;j<nRows;j++)
                                    A[j+i*nRows] = -A[j+i*nRows];
                            }
                        }


                        count_elems[Edof]++;
                        if(count_elems[Edof] == 6){
                            count_elems[Edof]=0;
                            count_DOFs[Edof]++;
                        }
                        break;
                    default:
                        std::cout<<"wrong template "<<templateR[i] - 1<<std::endl;
                        abort();
                }
            }
        }

        for(i=0;i<nRows;i++){
            if(!index_ghostC[i]){//cols
                if((indexesC[i]<MinInd) ||(indexesC[i]>MaxInd)){
                    std::cout<<"wrong index C "<<indexesC[i]<<" "<<MinInd<<" "<<MaxInd<<" "<<i<<std::endl;
                    abort();
                }
                for(j=0;j<nRows;j++){ //rows
                    {

                        if((indexesR[i]<0) ||(indexesR[i]>MatrSize)){
                            std::cout<<"Proc rank "<<m->GetProcessorRank()<<std::endl;
                            std::cout<<"wrong index R "<<indexesR[i]<<" "<<0<<" "<<MatrSize<<std::endl;
                            std::cout<<"flag_symmetric "<<flag_symmetric<<std::endl;
                            std::cout<<"Num elems"<<std::endl;
                            for(int i1=0;i1<4;i1++)
                                std::cout<<NumElem[i1]<<" ";
                            std::cout<<std::endl;
                            std::cout<<"Min elems"<<std::endl;
                            for(int i1=0;i1<4;i1++)
                                std::cout<<MinElem[i1]<<" ";
                            std::cout<<"Max elems"<<std::endl;
                            for(int i1=0;i1<4;i1++)
                                std::cout<<MaxElem[i1]<<" ";
                            std::cout<<std::endl;
                            std::cout<<"Min ind "<< MinInd<<" Max ind "<<MaxInd<< std::endl;
                            std::cout<<"indexes"<<std::endl;
                            for(int i1=0;i1<nRows;i1++)
                                std::cout<<indexesR[i1]<<" ";
                            std::cout<<std::endl;
                            abort();
                        }
                        if(A[j*nRows +i] != 0.0)
                            matrix[indexesC[i]][indexesR[j]] +=A[j*nRows +i];
                        if(fabs(A[j*nRows +i])>Max_A)
                            Max_A = fabs(A[j*nRows +i]);
                        if((fabs(A[j*nRows +i])<Min_A) && (A[j*nRows +i] != 0.0))
                            Min_A = fabs(A[j*nRows +i]);
                    }
                }
                rhs[indexesC[i]] +=F[i];
            }

        }
    }

  //  std::cout<<"@@@@@"<<std::endl;
        //std::cout<<"Min elem "<<Min_A<<" Max elem "<<Max_A<<std::endl;
    for(i=0;i<4;i++)
        if(NumDof[i]>0){

       // std::cout<<"tst "<<i<<std::endl;
            if(NumDof[i]<Leading_DOF){
                for(Mesh::iteratorElement it = m->BeginElement((ElementType)(NODE<<i)); it!= m->EndElement();it++)
                    if(it->GetStatus() != Element::Ghost)
                        for(j=0;j<Leading_DOF- NumDof[i];j++){
                            int index = it->IntegerArray(FakeNumTags[i])[j];
                           // std::cout<<"!! "<<index<<std::endl;
                            matrix[index][index] = 1.0;
                            rhs[index] = 0;
                        }
            }
        }
}
 */
//============================================================================================================================================
