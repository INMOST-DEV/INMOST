//
// Created by vasily on 3/17/17.
//

#ifndef INMOST_ANI_FEM_H
#define INMOST_ANI_FEM_H

enum ASSEMBLING_TYPE {
    MINIBLOCKS,
    ANITYPE
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
static const char * NumTagNames[4] = {"NumNodeTag","NumEdgeTag","NumFaceTag","NumCellTag"};
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
    void PrepareProblemAni ( int numNodeDofs, int numEdgeDOFS, int numFaceDofs, int numCellDofs,ASSEMBLING_TYPE type_a_,ORDER_TYPE type_o_);
    void PrepareProblemAni_test_for_fsi ( int numNodeDofs, int numEdgeDOFS, int numFaceDofs, int numCellDofs,ASSEMBLING_TYPE type_a_,ORDER_TYPE type_o_);
    ASSEMBLING_TYPE  type_assembling;
    ORDER_TYPE  type_order;
    void Assemble (Tag Label_Tag,Sparse::Matrix &matrix, Sparse::Vector &rhs, double *DDATA, int DDATA_size, int *IDATA, int IDATA_size);
//    void Assemble (Tag Label_Tag,Tag Cell_Tag, Tag Face_Tag,Tag Edge_Tag, Tag Node_Tag,Sparse::Matrix &matrix, Sparse::Vector &rhs, double *DDATA, int DDATA_size, int *IDATA, int IDATA_size);
    ~Ani_discretization(){};
    int getMinInd(){ return MinInd;}
    int getMaxInd(){ return MaxInd;}
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
