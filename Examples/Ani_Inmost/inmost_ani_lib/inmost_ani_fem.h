//
// Created by vasily on 3/17/17.
//

#ifndef INMOST_ANI_FEM_H
#define INMOST_ANI_FEM_H

enum ASSEMBLING_TYPE {
    MINIBLOCKS,
    ANITYPE,
    OSEEN_FORTRA_STOKES
};
enum ORDER_TYPE {
    STRAIGHT,
    REVERSE
};
#include "inmost_ani_mesh.h"
static const int Ndof_fortran = 0; //NodeDOF;
static const int Fdof_fortran = 1; //FaceDOF
static const int Edof_fortran = 2; //EdgeDOF
static const int Cdof_fortran = 3; //CellDOF
static const int Ndof = 0; //NodeDOF;
static const int Edof = 1; //EdgeDOF
static const int Fdof = 2; //FaceDOF
static const int Cdof = 3; //CellDOF
static const char * NumTagNames_fortran[4] = {"NumNodeTag_fortran", "NumFaceTag_fortran","NumEdgeTag_fortran","NumCellTag_fortran"};
static const char * NumTagNames_consts[4] = {"NumNodeTag","NumEdgeTag","NumFaceTag","NumCellTag"};
static const char * FakeNumTagNames_consts[4] = {"FakeNumNodeTag","FakeNumEdgeTag","FakeNumFaceTag","FakeNumCellTag"};
static const int FastInnerIndex[5] = {0,1,2,-1,3};



struct Ani_discretization{
    Ani_discretization(Mesh *m_);
    int NumDof[4];
    int NumElem[4];
    int MinElem[4];
    int MaxElem[4];
    int InitIndex[4];
    Tag FaceBackCells;//for FDoforient
    Tag NumTags[4];
    Mesh *m;
    int MinInd,MaxInd,MatrSize;
    std::string NumTagNames[4];
    std::string FakeNumTagNames[4];
    void PrepareProblemAni ( int numNodeDofs, int numEdgeDOFS, int numFaceDofs, int numCellDofs,ASSEMBLING_TYPE type_a_,ORDER_TYPE type_o_, std::string prefix="");
   // void PrepareProblemAni_test_for_fsi ( int numNodeDofs, int numEdgeDOFS, int numFaceDofs, int numCellDofs,ASSEMBLING_TYPE type_a_,ORDER_TYPE type_o_);
   // void Assemble_test_for_fsi (Tag Label_Tag,Sparse::Matrix &matrix, Sparse::Vector &rhs, double *DDATA, int DDATA_size, int *IDATA, int IDATA_size) ;

        ASSEMBLING_TYPE  type_assembling;
    ORDER_TYPE  type_order;
    void Assemble (Tag Label_Tag,Sparse::Matrix &matrix, Sparse::Vector &rhs, double *DDATA, int DDATA_size, int *IDATA, int IDATA_size, void (*fem3dext_new_)(
            double * XY1,double *  XY2,double *  XY3,double *  XY4,
            int * lbE,int *  lbF,int *  lbR,int *  lbP,double *  DDATA,int *  IDATA,int *  iSYS,
            int * LDA,double *  A,double *  F,int *  nRow,int *  nCol,
            int * templateR,int *  templateC) =nullptr);
//    void Assemble (Tag Label_Tag,Tag Cell_Tag, Tag Face_Tag,Tag Edge_Tag, Tag Node_Tag,Sparse::Matrix &matrix, Sparse::Vector &rhs, double *DDATA, int DDATA_size, int *IDATA, int IDATA_size);
    ~Ani_discretization(){}
    int getMinInd(){ return MinInd;}
    int getMaxInd(){ return MaxInd;}
    //==================================================================================
   // void PrepareProblemAni_LD( int numNodeDofs, int numEdgeDofs, int numFaceDofs, int numCellDofs, ASSEMBLING_TYPE type_a_, ORDER_TYPE type_o_);
  //  void Assemble_LD (Tag Label_Tag,Sparse::Matrix &matrix, Sparse::Vector &rhs, double *DDATA, int DDATA_size, int *IDATA, int IDATA_size);
    Tag FakeNumTags[4];
    int Leading_DOF;


    //==========================================================================================

};



/*
struct Ani_discretization{
    Ani_discretization(Mesh *m_);
    int NumDof[4];
    int NumElem[4];
    int MinElem[4];
    int MaxElem[4];
    int InitIndex[4];
    Tag NumTags[4];
    Mesh *m;
    int MinInd,MaxInd,MatrSize;
    void PrepareProblemAni ( int numNodeDofs, int numEdgeDOFS, int numFaceDofs, int numCellDofs,ASSEMBLING_TYPE type_a_,ORDER_TYPE type_o_);
    void PrepareProblemAni_new ( int numNodeDofs, int numEdgeDOFS, int numFaceDofs, int numCellDofs,ASSEMBLING_TYPE type_a_,ORDER_TYPE type_o_);
    ASSEMBLING_TYPE  type_assembling;
    ORDER_TYPE  type_order;
    void Assemble (Tag Cell_Tag, Tag Face_Tag,Tag Edge_Tag, Tag Node_Tag,Sparse::Matrix &matrix, Sparse::Vector &rhs, double *DDATA, int DDATA_size, int *IDATA, int IDATA_size);
    ~Ani_discretization(){};
    int getMinInd(){ return MinInd;}
    int getMaxInd(){ return MaxInd;}
};*/
#endif //INMOST_ANI_FEM_H
