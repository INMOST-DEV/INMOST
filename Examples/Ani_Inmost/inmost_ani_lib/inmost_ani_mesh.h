#ifndef INMOST_ANI_MESH_H
#define INMOST_ANI_MESH_H

#include <utility>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <map>
#include "inmost.h"
#include <fstream>
#include "assert.h"
using namespace INMOST;

//#define STATISTICS
#define NO_INTERNAL_MARKS

#define DebugCellTag "DebugCellTag"
#define DebugFaceTag "DebugFaceTag"
#define DebugEdgeTag "DebugEdgeTag"
#define DebugNodeTag "DebugNodeTag"


#define FortranCellTag "FrotranCellTag"
#define FortranNodeTag "FrotranNodeTag"




#define FortranNumTag "Fortran_number"
#define LabelTag "Label"
const int No_BC_mark = 0;
#define EPS 1e-10

struct Ani_Mesh {
	int nv,nt,nb;
	int nvmax,ntmax,nbmax;

	int n_curv;
	int nvfix,nbfix,ntfix;

	double *vrt;
	int *bnd,*tet;
	int *labelF,*labelT;
	int *ivfix,*ibfix,*itfix;
    unsigned long int getsize();
	Ani_Mesh();
	~Ani_Mesh();
};
inline bool comp_numb(int n, int num[3])
{
    if((n==num[0])||(n==num[1])||(n==num[2]))
        return true;
    else
        return false;
};


void CreateAniMeshFromInmost(Mesh *m,  Ani_Mesh &animesh);
void ReadInmostFromAniFile(INMOST::Mesh *m,std::string filename );
void RefineLocalMesh(Mesh *m,Mesh *m_result,   int NumberOfRefines);
void ReadAniFromFile(std::string filename,Ani_Mesh &m );
void LoadAFT(std::string filename,INMOST::Mesh *m);
void CreateInmostMeshFromAni(INMOST::Mesh *m,Ani_Mesh &animesh) ;
void CreateInmostMeshFromAni_with_interial_marks(INMOST::Mesh *m,Ani_Mesh &animesh) ;
#define M_Assert(Expr, Msg) \
    if (!Expr)  \
    { \
    std::cerr << "Assert failed:\t" << Msg << "\n" \
        << "Expected:\t" << #Expr << "\n" \
        << "Source:\t\t" << __FILE__ << ", line " << __LINE__ << "\n"; \
        abort(); \
    }

#define M_Warinig(Expr, Msg) \
    if (!Expr)  \
    { \
    std::cerr << "WARNING:\t" << Msg << "\n" \
        << "Expected:\t" << #Expr << "\n" \
        << "Source:\t\t" << __FILE__ << ", line " << __LINE__ << "\n"; \
    }

#endif
