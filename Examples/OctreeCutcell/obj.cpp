#include "obj.h"
#include "triang.h"
#include "tree.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "proj.hpp"
#include <vector>
#define DEEP_START 3
#define EPSCLOSE 10e-17
#define Eps 1e-8

int nowarn_obj = 0;
char * nowarn2 = NULL;

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define SQR(a)   ((a)*(a))



//~ extern Projection make_proj;

bool swapyz = true;


struct ObjVertex
{
	double v[3];
};

struct ObjNormal
{
	double n[3];
};

struct ObjTex
{
    double uv[2];
};

struct ObjMat
{
    int istex;
    float r,g,b,a;
    char texture[1024];
    char name[1024];
};

struct ObjMatuse
{
    char name[1024];
};

struct ObjPolygon
{
	int * vert, * norm, * tex;
	int nvert;
	int nnorm;
    int ntex;
    struct ObjMatuse material;
	double minx, maxx, miny,maxy, minz,maxz;
};



#define CELLBUF 1024


struct ObjInfo
{
	bool used;

	double scale;
	double coord[3];
//	double rot[3];
	struct quat rot;

	double inverse[16];


	double minx, miny, minz, maxx, maxy,maxz;
	struct ObjVertex * verts;
	struct ObjNormal * norms;
    struct ObjTex * texs;
    struct ObjMat * mats;
	struct point * p;
	struct tree * t;
	int nverts;
	int nnorms;
    int ntexs;
    int nmats;
	struct ObjPolygon * pols;
	int npols;

	char filename[1024];
//	char filedir[1024];

    int nmaterials;
} objs[OBJ_MAX];




int nobj = 0;


int NumberofObj()
{
	return nobj;
}

inline void vecmul(double * vecin, double * vecout, double * mat)
{
	int i,j;
	double temp[4] = {0.0,0.0,0.0,1.0};
	double tempout[4] = {0,0,0,0};
	temp[0] = vecin[0];
	temp[1] = vecin[1];
	temp[2] = vecin[2];
	//memcpy(temp,vecin,sizeof(double)*3);
	for(i = 0; i < 4; i++)
		for(j = 0; j < 4; j++)
			tempout[i] += temp[j]*mat[j+i*4];
	for(i = 0; i < 3; i++)
		tempout[i] /= tempout[3];
	//memcpy(vecout,tempout,sizeof(double)*3);
	vecout[0] = tempout[0];
	vecout[1] = tempout[1];
	vecout[2] = tempout[2];
}



inline double dotproduct(double * vecin1, double * vecin2)
{
	return vecin1[0]*vecin2[0]+vecin1[1]*vecin2[1]+vecin1[2]*vecin2[2];
}

inline void crossproduct(double * vecin1, double * vecin2, double * vecout)
{
	vecout[0] = vecin1[1]*vecin2[2] - vecin1[2]*vecin2[1];
	vecout[1] = vecin1[2]*vecin2[0] - vecin1[0]*vecin2[2];
	vecout[2] = vecin1[0]*vecin2[1] - vecin1[1]*vecin2[0];
}


inline void vecdiff(double * vecin1, double * vecin2, double * vecout)
{
	vecout[0] = vecin2[0] - vecin1[0];
	vecout[1] = vecin2[1] - vecin1[1];
	vecout[2] = vecin2[2] - vecin1[2];
}

inline void normalize(double * vec)
{
	double l = 1.0 / sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
	if( fabs(l) > 1e-15 )
	{
		vec[0] *= l;
		vec[1] *= l;
		vec[2] *= l;
	}
}

void quatprint(quat a)
{
	printf("quat: %f %f %f %f\n",a.vec[0],
		   a.vec[1],a.vec[2],a.w);
}
quat quatmul(quat a, quat b)
{
	quat ret;
	crossproduct(a.vec,b.vec,ret.vec);
	ret.vec[0] += a.w*b.vec[0];
	ret.vec[1] += a.w*b.vec[1];
	ret.vec[2] += a.w*b.vec[2];
	ret.vec[0] += b.w*a.vec[0];
	ret.vec[1] += b.w*a.vec[1];
	ret.vec[2] += b.w*a.vec[2];
	ret.w = a.w*b.w - dotproduct(a.vec,b.vec);
	printf("\n");
	return ret;
}

quat quatfromaxis(double x, double y, double z, double ang)
{
	quat ret;
	double hang = ang*0.5;
	double sina = sin(hang);
	ret.vec[0] = x*sina;
	ret.vec[1] = y*sina;
	ret.vec[2] = z*sina;
	ret.w =cos(hang);
	return ret;
}

quat quatdiv(quat a, double b)
{
	quat ret = a;
	ret.vec[0] /= b;
	ret.vec[1] /= b;
	ret.vec[2] /= b;
	ret.w /= b;
	return ret;
}

quat quatconj(quat a)
{
	quat ret = a;
	ret.vec[0] = -ret.vec[0];
	ret.vec[1] = -ret.vec[1];
	ret.vec[2] = -ret.vec[2];
	return ret;
}

void quatident(quat * a)
{
	a->vec[0] = 0.0;
	a->vec[1] = 0.0;
	a->vec[2] = 0.0;
	a->w = 1.0;
}
void quatzero(quat * a)
{
	a->vec[0] = 0.0;
	a->vec[1] = 0.0;
	a->vec[2] = 0.0;
	a->w = 0.0;
}
quat quatrot(quat a, quat b)
{
	return quatmul(quatmul(a,b),quatconj(a));
}

double quatnorm(quat a)
{
	return a.vec[0]*a.vec[0] + a.vec[1]*a.vec[1]
		 + a.vec[2]*a.vec[2] + a.w*a.w;
}

void quatmat(quat q, double * rot)
{
	rot[ 0] = (q.w*q.w + q.vec[0]*q.vec[0] - q.vec[1]*q.vec[1] - q.vec[2]*q.vec[2]);
	rot[ 1] = 2.*(q.vec[0]*q.vec[1] - q.w*q.vec[2]);
	rot[ 2] = 2.*(q.vec[0]*q.vec[2] + q.w*q.vec[1]);
	rot[ 3] = 0.0;
	rot[ 4] = 2.*(q.vec[0]*q.vec[1] + q.w*q.vec[2]);
	rot[ 5] = (q.w*q.w - q.vec[0]*q.vec[0] + q.vec[1]*q.vec[1] - q.vec[2]*q.vec[2]);
	rot[ 6] = 2.*(q.vec[1]*q.vec[2] - q.w*q.vec[0]);
	rot[ 7] = 0.0;
	rot[ 8] = 2.*(q.vec[0]*q.vec[2] - q.w*q.vec[1]);
	rot[ 9] = 2.*(q.vec[1]*q.vec[2] + q.w*q.vec[0]);
	rot[10] = (q.w*q.w - q.vec[0]*q.vec[0] - q.vec[1]*q.vec[1] + q.vec[2]*q.vec[2]);
	rot[11] = 0.0;
	rot[12] = 0.0;
	rot[13] = 0.0;
	rot[14] = 0.0;
	rot[15] = (q.w*q.w + q.vec[0]*q.vec[0] + q.vec[1]*q.vec[1] + q.vec[2]*q.vec[2]);
}

quat quatinv(quat a)
{
	return quatdiv(quatconj(a),quatnorm(a));
}

void matmul(double * a, double * b)
{
	int i,j,k;
	double ret[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	for(i = 0; i < 4; i++)
		for(j = 0; j < 4; j++)
			for(k = 0; k < 4; k++)
			{
				ret[i*4+j] += a[i*4+k]*b[k*4+j];
			}
	for(i = 0; i < 16; i++)
		a[i] = ret[i];
}


inline float cofactor_ij_v1(double * mat, int col,int row)
{
    static const int sel0[] = { 1,0,0,0 };
    static const int sel1[] = { 2,2,1,1 };
    static const int sel2[] = { 3,3,3,2 };
    // Let's first define the 3x3 matrix:
    const int col0 = sel0[col];
    const int col1 = sel1[col];
    const int col2 = sel2[col];
    const int row0 = sel0[row];
    const int row1 = sel1[row];
    const int row2 = sel2[row];
    // Compute the float sign-mask:
    const unsigned	int signpos  = (col + row) & 1;
    const unsigned	int signmask = signpos << 31;
    // Computer the det of the 3x3 matrix:
    //
    //   [ a b c ]
    //   [ d e f ] = aei - ahf + dhc - dbi + gbf - gec = (aei + dhc + gbf) - (ahf + dbi + gec)
    //   [ g h i ]
    //
    const double pos_part1 = mat[col0+row0*4] * mat[col1+row1*4] * mat[col2+row2*4]; // aei
    const double pos_part2 = mat[col0+row1*4] * mat[col1+row2*4] * mat[col2+row0*4]; // dhc
    const double pos_part3 = mat[col0+row2*4] * mat[col1+row0*4] * mat[col2+row1*4]; // gbf
    const double neg_part1 = mat[col0+row0*4] * mat[col1+row2*4] * mat[col2+row1*4]; // ahf
    const double neg_part2 = mat[col0+row1*4] * mat[col1+row0*4] * mat[col2+row2*4]; // dbi
    const double neg_part3 = mat[col0+row2*4] * mat[col1+row1*4] * mat[col2+row0*4]; // gec
    const double pos_part  = pos_part1 + pos_part2 + pos_part3;
    const double neg_part  = neg_part1 + neg_part2 + neg_part3;
    const double det3x3	  = pos_part - neg_part;
    // Now let's adjust the sign of the cofactor:
    union
    {
        float f;
        int	  i;
    } cofactor;
    cofactor.f  = det3x3;
    cofactor.i ^= signmask;
    return cofactor.f;
}


inline void cofactor_v1(double * output, double * source)
{
    int col,row;
    for ( col = 0 ; col < 4 ; col++ )
        for ( row = 0; row < 4; row++ )
            output[col+row*4] =  cofactor_ij_v1(source, col, row);
}
inline float determinant_v1(double * source, double * cofactor)
{
    int col;
    float det = 0.0f;
    for ( col = 0; col < 4; col++ )
        det += source[col] * cofactor[col];
    return det;
}

inline void transpose_v1(double * output, double * source)
{
    int index1,index2;
    for ( index1 = 0; index1 < 4; index1++ )
        for ( index2 = 0; index2 < 4; index2++ )
            output[index2+index1*4] = source[index1+index2*4];
}

inline void mul_v1(double * output, double * source, const float factor)
{
    int col,row;
    for ( col = 0 ; col < 4 ; col++ )
        for ( row = 0; row < 4; row++ )
            output[col+row*4] =  source[col+row*4] * factor;
}

inline void transpone(double * mat)
{
    int i,j;
    for(i = 0; i < 4; i++)
        for(j = i+1; j < 4; j++)
        {
            double temp = mat[i*4+j];
            mat[i*4+j] = mat[i+j*4];
            mat[i+j*4] = temp;
        }
}

inline void inverse_v1(double * output, double * source)
{
    double cof[16];
    double adj[16];
    float    oodet;
    // InvM = (1/det(M)) * Transpose(Cofactor(M))
    cofactor_v1(cof, source);
    oodet = 1.0f / determinant_v1(source, cof);
    transpose_v1(adj, cof);
    mul_v1(output, adj, oodet);
}

void PrepareMatricesObj(double * matrix, double * inverse, double * rotateonly, int n)
{
    if( n > nobj || n < 0 ) return;
	if( !objs[n].used ) return;
    double rotate[16];
	quatmat(objs[n].rot,rotate);
	int i;
    if( rotateonly != NULL )
	{
		for(i = 0; i < 16; i++)
			rotateonly[i] = rotate[i];
		//memcpy(rotateonly,rotate,sizeof(double)*16);
	}
		/*
    double prerot[16] =
    {
        1.0,0.0,0.0,0.0,
		0.0,1.0,0.0,0.0,
		0.0,0.0,1.0,0.0,
		-objs[n].center_mass[0]*objs[n].scale,-objs[n].center_mass[1]*objs[n].scale,-objs[n].center_mass[2]*objs[n].scale,1.0
    };
    transpone(prerot);
    double postrot[16] =
    {
        1.0,0.0,0.0,0.0,
		0.0,1.0,0.0,0.0,
		0.0,0.0,1.0,0.0,
		objs[n].center_mass[0]*objs[n].scale,objs[n].center_mass[1]*objs[n].scale,objs[n].center_mass[2]*objs[n].scale,1.0
    };
    transpone(postrot);*/
	double scale[16] =
	{
		objs[n].scale,0.0,0.0,0.0,
		0.0,objs[n].scale,0.0,0.0,
		0.0,0.0,objs[n].scale,0.0,
		0.0,0.0,0.0,1.0
	};
	double translate[16] =
	{
		1.0,0.0,0.0,0.0,
		0.0,1.0,0.0,0.0,
		0.0,0.0,1.0,0.0,
		objs[n].coord[0],objs[n].coord[1],objs[n].coord[2],1.0
	};
    transpone(translate);
	double identity[16] =
	{
		1.0,0.0,0.0,0.0,
		0.0,1.0,0.0,0.0,
		0.0,0.0,1.0,0.0,
		0.0,0.0,0.0,1.0
	};
	matmul(identity,translate);
//    matmul(identity,postrot);
	matmul(identity,rotate);
//    matmul(identity,prerot);
	matmul(identity,scale);
    if( matrix != NULL )
	{
		for(i = 0; i < 16; i++)
			matrix[i] = identity[i];
//        memcpy(matrix,identity,sizeof(double)*16);
	}
    if( inverse != NULL )
        inverse_v1(inverse,identity);
}


void SetRotateObj(double x, double y, double z, int n)
{
	quat a;
	if( n > nobj || n < 0 ) return;
	if( !objs[n].used ) return;
	quatident(&objs[n].rot);
	a = quatfromaxis(1.0,0.0,0.0,x);
	objs[n].rot = quatmul(objs[n].rot,a);
	a = quatfromaxis(0.0,1.0,0.0,y);
	objs[n].rot = quatmul(objs[n].rot,a);
	a = quatfromaxis(0.0,0.0,1.0,z);
	objs[n].rot = quatmul(objs[n].rot,a);
}
void SetScaleObj(double s,int n)
{
	if( n > nobj || n < 0 ) return;
	if( !objs[n].used ) return;
	objs[n].scale = s;
}
void SetTranslateObj(double x, double y, double z,int n)
{
	if( n > nobj || n < 0 ) return;
	if( !objs[n].used ) return;
	objs[n].coord[0] = x;
	objs[n].coord[1] = y;
	objs[n].coord[2] = z;
}

void RotateObj(double rx, double ry, double rz, int n)
{
	quat a;
	if( n > nobj || n < 0 ) return;
	if( !objs[n].used ) return;
	//printf("x\n");
	a = quatfromaxis(1.0,0.0,0.0,rx);
	objs[n].rot = quatmul(objs[n].rot,a);
	//printf("y\n");
	a = quatfromaxis(0.0,1.0,0.0,ry);
	objs[n].rot = quatmul(objs[n].rot,a);
	//printf("z\n");
	a = quatfromaxis(0.0,0.0,1.0,rz);
	objs[n].rot = quatmul(objs[n].rot,a);
	//PrintCoordObj(n);
}

void ScaleObj(double s, int n)
{
	if( n > nobj || n < 0 ) return;
	if( !objs[n].used ) return;
	objs[n].scale += s;
}
void TranslateObj(double dx, double dy, double dz,int n)
{
	if( n > nobj || n < 0 ) return;
	if( !objs[n].used ) return;
	objs[n].coord[0] += dx;
	objs[n].coord[1] += dy;
	objs[n].coord[2] += dz;
}

#define BIG_SIZE 1000
#define SMALL_SIZE 50

int findspace(char * str)
{
	int i;
	for(i = 0; str[i] != '\0'; i++)
		if( str[i] == ' ' )
			return i;
	return -1;
}


int ReadObj(char * file)
{
    char mtlname[1024];
	char str[4096];
	int i,j,k,m,cur;
	struct ObjInfo * n;
    strcpy(mtlname,"\0");
	cur = 0;
	for(i = 0; i < nobj; i++)
			if( objs[i].used == false )
			{
				cur = i;
				break;
			}
	if( cur == OBJ_MAX)
	{
		printf("No space for new objects, enlarge OBJ_MAX in obj.h\n");
		return -1;
	}
	else
	{
		cur = nobj;
	}
	n = objs + cur;
	strcpy(n->filename,file);
	printf("Reading %s\n",file);
//	nowarn2 = getcwd(n->filedir,1024);
	FILE * f = fopen(file,"r");
	if( f == NULL )
	{
		printf("ReadObj: cannot open file %s\n",file);
		return -1;
	}
	if( cur == nobj ) nobj++;
	for(i = 0; i < 3; i++)
	{
		n->coord[i] = 0;
//		n->rot[i] = 0;
		n->rot.vec[i] = 0;
		n->scale = 1.f;
	}
	n->rot.w = 1.0;
//	n->interactable = definteract;
	n->verts = (struct ObjVertex *) malloc(sizeof(struct ObjVertex)*BIG_SIZE);
	n->norms = (struct ObjNormal *) malloc(sizeof(struct ObjNormal)*BIG_SIZE);
	n->texs = (struct ObjTex *) malloc(sizeof(struct ObjTex)*BIG_SIZE);
    n->mats = (struct ObjMat *) malloc(sizeof(struct ObjMat)*SMALL_SIZE);
	n->nverts = 0;
	n->nnorms = 0;
	n->ntexs = 0;
    n->nmats = 0;
	n->pols = (struct ObjPolygon *) malloc(sizeof(struct ObjPolygon)*SMALL_SIZE);
    if( n->verts == NULL || n->norms == NULL || n->texs == NULL || n->pols == NULL || n->mats == NULL)
    {
        printf("Out of memory in ReadObj line %d\n Stop\n",__LINE__);
        exit(-1);
    }
	n->npols = 0;
	n->used = true;
	while( fgets(str,4096,f) != NULL)
	{
		//printf("read string: %s\n",str);
		switch(str[0])
		{
            case 'm':
                //printf("%s\n",str);
                if( !strncmp(str,"mtllib",6) )
                {
                    char str2[4096];
                    char dir[1024];
                    char tmp[1024];
                    char libname[1024];
                    FILE * fff;
                    i = strlen(file)-1;
                    while( file[i] != '/' && file[i] != '\\' && i >= 0 ) i--;
                    strncpy(dir,file,i+1);
                    dir[i+1] = '\0';
                    sscanf(str+7,"%s\n",libname);
                    strcpy(tmp,dir);
                    strcat(tmp,libname);
                    printf("Loading material library file %s\n",tmp);
                    fff = fopen(tmp,"r");
                    if( fff == NULL )
                    {
                        printf("Faild\n");
                        break;
                    }
					//printf("Found\n");
                    while( fgets(str2,4096,fff) != NULL )
                    {
                        switch(str2[0])
                        {
                            case '#': //comment
                                break;
                            case 'i':// skip illnum
                                break;
                            case 'n':
                                if( !strncmp(str2,"newmtl",6)  )
                                {

                                    sscanf(str2+7,"%s\n",n->mats[n->nmats].name);
//                                    printf("Found new material %s\n",n->mats[n->nmats].name);
                                    n->mats[n->nmats].istex = 0;
                                    n->mats[n->nmats].a = 1.0;
                                    n->nmats++;
                                    if( n->nmats % SMALL_SIZE == 0 )
                                    {
                                        n->mats = (struct ObjMat *) realloc(n->mats,
                                                                            sizeof(struct ObjMat)*
                                                                            (n->nmats/SMALL_SIZE+1)*SMALL_SIZE);
                                        if( n->mats == NULL )
                                        {
                                            printf("Out of memory in ReadObj line %d\n Stop\n",__LINE__);
                                            exit(-1);
                                        }
                                    }
                                }
                                break;
                            case 'T':
                                if( str2[1] == 'r' )
                                {
                                    sscanf(str2,"Tr %f\n",
                                           &n->mats[n->nmats-1].a);
                                }
                                break;
                            case 'K':
                                if( str2[1] == 'd' )
                                {
                                    sscanf(str2,"Kd %f %f %f\n",
                                           &n->mats[n->nmats-1].r,
                                           &n->mats[n->nmats-1].g,
                                           &n->mats[n->nmats-1].b);
                                }
                                break;
                            case 'm':
                                if( !strncmp(str2,"map_Kd",6)  )
                                {
                                    n->mats[n->nmats-1].istex = 1;
                                    i = strlen(str)-1;
                                    while( str2[i] != ' ' && i > 0 ) i--;
                                    sscanf(str2+i+1,"%s",n->mats[n->nmats-1].texture);
                                }
                                break;
                            case 'd':
                            	if( str2[1] == ' ' )
                                {
                                	sscanf(str2,"d %f\n",&n->mats[n->nmats-1].a);
                                }

                        }
                    }
                    fclose(fff);
                }
                break;
            case 'u':
                if( !strncmp(str,"usemtl",6) )
                {
                    sscanf(str+7,"%s\n",mtlname);
//                    printf("Found material usage %s\n",mtlname);
                }
                break;
			case '#':
				//printf("Comment found\n");
				//comments
				break;
			case 'f':
				//printf("Face found\n");
				n->pols[n->npols].nvert = 0;
				n->pols[n->npols].nnorm = 0;
				n->pols[n->npols].ntex = 0;
				n->pols[n->npols].vert = (int *)malloc(sizeof(int)*SMALL_SIZE);
				n->pols[n->npols].norm = (int *)malloc(sizeof(int)*SMALL_SIZE);
				n->pols[n->npols].tex = (int *)malloc(sizeof(int)*SMALL_SIZE);
                if(n->pols[n->npols].vert == NULL || n->pols[n->npols].norm == NULL ||
                    n->pols[n->npols].tex == NULL)
                {
                    printf("Out of memory in ReadObj line %d\n Stop\n",__LINE__);
                    exit(-1);
                }
				k = 0;
				m = 0;
				while((k = findspace(str+k)) != -1)
				{
					k = k+m+1;
				//	printf("Parsing %s, space after %d: %s\n",str,k,str+k);
					m = sscanf(str+k,"%d/%d/%d",
							   &(n->pols[n->npols].vert[n->pols[n->npols].nvert]),
							   &(n->pols[n->npols].tex[n->pols[n->npols].ntex]),
							   &(n->pols[n->npols].norm[n->pols[n->npols].nnorm]));
					if( m > 0 )
						n->pols[n->npols].nvert++;
                    if( m > 1 )
                        n->pols[n->npols].ntex++;
					if( m > 2 )
						n->pols[n->npols].nnorm++;
					if( n->pols[n->npols].nvert%SMALL_SIZE == 0 && n->pols[n->npols].nvert != 0)
                    {
						n->pols[n->npols].vert = (int *)realloc(n->pols[n->npols].vert,
																sizeof(int)*(n->pols[n->npols].nvert/SMALL_SIZE+1)*SMALL_SIZE);
                        if(n->pols[n->npols].vert == NULL)
                        {
                            printf("Out of memory in ReadObj line %d\n Stop\n",__LINE__);
                            exit(-1);
                        }
                    }
                    if( n->pols[n->npols].ntex%SMALL_SIZE == 0 && n->pols[n->npols].ntex != 0)
                    {
						n->pols[n->npols].tex = (int *)realloc(n->pols[n->npols].tex,
																sizeof(int)*(n->pols[n->npols].ntex/SMALL_SIZE+1)*SMALL_SIZE);
                        if(n->pols[n->npols].tex == NULL)
                        {
                            printf("Out of memory in ReadObj line %d\n Stop\n",__LINE__);
                            exit(-1);
                        }
                    }
					if( n->pols[n->npols].nnorm%SMALL_SIZE == 0 && n->pols[n->npols].nnorm != 0)
                    {
						n->pols[n->npols].norm = (int *)realloc(n->pols[n->npols].norm,
																sizeof(int)*(n->pols[n->npols].nnorm/SMALL_SIZE+1)*SMALL_SIZE);
                        if(n->pols[n->npols].norm == NULL)
                        {
                            printf("Out of memory in ReadObj line %d\n Stop\n",__LINE__);
                            exit(-1);
                        }
                    }
					m = k;
				}
                strcpy(n->pols[n->npols].material.name,mtlname);
				n->npols++;
				if( n->npols%SMALL_SIZE == 0)
				{
					n->pols = (struct ObjPolygon *) realloc(n->pols,sizeof(struct ObjPolygon)
															*(n->npols/SMALL_SIZE+1)*SMALL_SIZE);
                    if(n->pols == NULL)
                    {
                        printf("Out of memory in ReadObj line %d\n Stop\n",__LINE__);
                        exit(-1);
                    }
				}
				break;
			case 'v':
				//printf("Some vertex information found\n");
				switch(str[1])
				{
					case 't':
				//		printf("Texture coordinate\n");
						// texture coordinate, skip
                        sscanf(str+3,"%lf %lf",
							   &(n->texs[n->ntexs].uv[0]),
							   &(n->texs[n->ntexs].uv[1]));
						n->ntexs++;
						if(n->ntexs % BIG_SIZE == 0)
						{
							n->texs = (struct ObjTex *) realloc(n->texs,sizeof(struct ObjTex)
                                                                *(n->ntexs/BIG_SIZE+1)*BIG_SIZE);
                            if(n->texs == NULL)
                            {
                                printf("Out of memory in ReadObj line %d\n Stop\n",__LINE__);
                                exit(-1);
                            }
						}
						break;
					case 'n':
				//		printf("Normale\n");
						// normale
						sscanf(str+3,"%lf %lf %lf",
							   &(n->norms[n->nnorms].n[0]),
							   &(n->norms[n->nnorms].n[1]),
							   &(n->norms[n->nnorms].n[2]));
						n->nnorms++;
						if(n->nnorms % BIG_SIZE == 0)
						{
							n->norms = (struct ObjNormal *) realloc(n->norms,sizeof(struct ObjNormal)
																	*(n->nnorms/BIG_SIZE+1)*BIG_SIZE);
                            if(n->norms == NULL)
                            {
                                printf("Out of memory in ReadObj line %d\n Stop\n",__LINE__);
                                exit(-1);
                            }
						}
						break;
					case ' ':
				//		printf("Vertex\n");
						sscanf(str+2,"%lf %lf %lf",
							   &n->verts[n->nverts].v[0],
							   &n->verts[n->nverts].v[1],
							   &n->verts[n->nverts].v[2]);

						n->verts[n->nverts].v[0] /= 99.5;
						n->verts[n->nverts].v[1] /= 32;
						n->verts[n->nverts].v[2] /= 99.5;
						n->verts[n->nverts].v[0] -= 0.001;
						n->verts[n->nverts].v[1] -= 0.001;
						n->verts[n->nverts].v[2] += 0.001;
						n->verts[n->nverts].v[2] *= -1;
						if( swapyz )
						{
							double temp = n->verts[n->nverts].v[1];
							n->verts[n->nverts].v[1] = n->verts[n->nverts].v[2];
							n->verts[n->nverts].v[2] = temp;
						}
						double alpha = n->verts[n->nverts].v[2];
						double fl = 0;//make_proj.GetFL()(n->verts[n->nverts].v[0],n->verts[n->nverts].v[1]);
						double fh = 1;//make_proj.GetFH()(n->verts[n->nverts].v[0],n->verts[n->nverts].v[1]);
						n->verts[n->nverts].v[2] = fl+(fh-fl)*alpha;
						n->nverts++;
						if(n->nverts % BIG_SIZE == 0)
						{
							n->verts = (struct ObjVertex *) realloc(n->verts,sizeof(struct ObjVertex)
																	*(n->nverts/BIG_SIZE+1)*BIG_SIZE);
                            if(n->verts == NULL)
                            {
                                printf("Out of memory in ReadObj line %d\n Stop\n",__LINE__);
                                exit(-1);
                            }
						}

						// just vertex
						break;
				}
				break;
			case 'g':
				//printf("Group found\n");
				//groupping tag, skipping
				break;
			case 's':
				//printf("Group 2 found\n");
				//another grouping tag, skipping
				break;
//			case 'u':
				//printf("Material found\n");
				//material tag, skipping
				break;
			default:
				//printf("Something strange found\n");
				break;
		}
	}
	fclose(f);
	n->maxx = -10e20;
	n->minx = 10e20;
	n->maxy = -10e20;
	n->miny = 10e20;
	n->maxz = -10e20;
	n->minz = 10e20;
	n->p = (struct point *)malloc(sizeof(struct point)*n->npols);
	n->t = (struct tree *)malloc(sizeof(struct tree));
    if(n->p == NULL || n->t == NULL )
    {
        printf("Out of memory in ReadObj line %d\n Stop\n",__LINE__);
        exit(-1);
    }
	kdtree_mem_alloc(n->t,n->npols);
	for(i = 0; i < n->npols; i++)
	{
		n->pols[i].maxx = -10e20;
		n->pols[i].minx = 10e20;
		n->pols[i].maxy = -10e20;
		n->pols[i].miny = 10e20;
		n->pols[i].maxz = -10e20;
		n->pols[i].minz = 10e20;
		if( n->pols[i].nnorm == 0 && n->pols[i].nvert > 2)
		{
			double a[3],b[3],nv[3],l;
			a[0] = n->verts[n->pols[i].vert[0]-1].v[0] - n->verts[n->pols[i].vert[2]-1].v[0];
			a[1] = n->verts[n->pols[i].vert[0]-1].v[1] - n->verts[n->pols[i].vert[2]-1].v[1];
			a[2] = n->verts[n->pols[i].vert[0]-1].v[2] - n->verts[n->pols[i].vert[2]-1].v[2];
			b[0] = n->verts[n->pols[i].vert[1]-1].v[0] - n->verts[n->pols[i].vert[2]-1].v[0];
			b[1] = n->verts[n->pols[i].vert[1]-1].v[1] - n->verts[n->pols[i].vert[2]-1].v[1];
			b[2] = n->verts[n->pols[i].vert[1]-1].v[2] - n->verts[n->pols[i].vert[2]-1].v[2];
			nv[0] = a[1]*b[2] - a[2]*b[1];
			nv[1] = a[2]*b[0] - a[0]*b[2];
			nv[2] = a[0]*b[1] - a[1]*b[0];
			l = sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);
			nv[0] /= l;
			nv[1] /= l;
			nv[2] /= l;
//			printf("making normale %f %f %f\n",nv[0],nv[1],nv[2]);
			n->norms[n->nnorms].n[0] = nv[0];
			n->norms[n->nnorms].n[1] = nv[1];
			n->norms[n->nnorms].n[2] = nv[2];
			n->pols[i].norm[0] = n->nnorms+1;
			n->pols[i].nnorm++;
			n->nnorms++;
			if(n->nnorms % BIG_SIZE == 0)
			{
				n->norms = (struct ObjNormal *) realloc(n->norms,sizeof(struct ObjNormal)
														*(n->nnorms/BIG_SIZE+1)*BIG_SIZE);
                if(n->norms == NULL)
                {
                    printf("Out of memory in ReadObj line %d\n Stop\n",__LINE__);
                    exit(-1);
                }
			}
		}
		//p[i].coord[0] = 0.0;
		//p[i].coord[1] = 0.0;
		//p[i].coord[2] = 0.0;
		for(j = 0; j < n->pols[i].nvert; j++)
		{
			//p[i].coord[0] += n->verts[n->pols[i].vert[j]-1].v[0];
			//p[i].coord[1] += n->verts[n->pols[i].vert[j]-1].v[1];
			//p[i].coord[2] += n->verts[n->pols[i].vert[j]-1].v[2];

			if( n->verts[n->pols[i].vert[j]-1].v[0] > n->pols[i].maxx )
				n->pols[i].maxx = n->verts[n->pols[i].vert[j]-1].v[0];
			if( n->verts[n->pols[i].vert[j]-1].v[0] < n->pols[i].minx )
				n->pols[i].minx = n->verts[n->pols[i].vert[j]-1].v[0];

			if( n->verts[n->pols[i].vert[j]-1].v[1] > n->pols[i].maxy )
				n->pols[i].maxy = n->verts[n->pols[i].vert[j]-1].v[1];
			if( n->verts[n->pols[i].vert[j]-1].v[1] < n->pols[i].miny )
				n->pols[i].miny = n->verts[n->pols[i].vert[j]-1].v[1];

			if( n->verts[n->pols[i].vert[j]-1].v[2] > n->pols[i].maxz )
				n->pols[i].maxz = n->verts[n->pols[i].vert[j]-1].v[2];
			if( n->verts[n->pols[i].vert[j]-1].v[2] < n->pols[i].minz )
				n->pols[i].minz = n->verts[n->pols[i].vert[j]-1].v[2];
		}
		//p[i].coord[0] /= (float) n->pols[i].nvert;
		//p[i].coord[1] /= (float) n->pols[i].nvert;
		//p[i].coord[2] /= (float) n->pols[i].nvert;
		n->p[i].coord[0] = (n->pols[i].maxx + n->pols[i].minx)/2.0;
		n->p[i].coord[1] = (n->pols[i].maxy + n->pols[i].miny)/2.0;
		n->p[i].coord[2] = (n->pols[i].maxz + n->pols[i].minz)/2.0;
		n->p[i].radius[0] = (n->pols[i].maxx - n->pols[i].minx)/2.0;
		n->p[i].radius[1] = (n->pols[i].maxy - n->pols[i].miny)/2.0;
		n->p[i].radius[2] = (n->pols[i].maxz - n->pols[i].minz)/2.0;
		n->p[i].polynum = i;

		if( n->maxx < n->pols[i].maxx )
			n->maxx = n->pols[i].maxx;
		if( n->minx > n->pols[i].minx )
			n->minx = n->pols[i].minx;

		if( n->maxy < n->pols[i].maxy )
			n->maxy = n->pols[i].maxy;
		if( n->miny > n->pols[i].miny )
			n->miny = n->pols[i].miny;

		if( n->maxz < n->pols[i].maxz )
			n->maxz = n->pols[i].maxz;
		if( n->minz > n->pols[i].minz )
			n->minz = n->pols[i].minz;
	}
	printf("bounds: %g:%g %g:%g %g:%g\n",n->minx,n->maxx,n->miny,n->maxy,n->minz,n->maxz);
	kdtree_build(n->t,n->p,n->npols,0);
	/*
	double  center[3];
	center[0] = (n->maxx+n->minx)/2.0;
	center[1] = (n->maxy+n->miny)/2.0;
//	center[2] = (n->maxz+n->minz)/2.0;
	center[2] = (n->minz);
	double sx = fabs(n->pols[0].maxx-center[0]);
	for(i = 0; i < n->npols; i++)
	{
	if( sx < fabs(n->pols[i].maxx-center[0]))
	sx = fabs(n->pols[i].maxx-center[0]);
	if( sx < fabs(n->pols[i].minx-center[0]))
	sx = fabs(n->pols[i].minx-center[0]);


	if( sx < fabs(n->pols[i].maxy-center[1]))
	sx = fabs(n->pols[i].maxy-center[1]);
	if( sx < fabs(n->pols[i].miny-center[1]))
	sx = fabs(n->pols[i].miny-center[1]);


	if( sx < fabs(n->pols[i].maxz-center[2]))
	sx = fabs(n->pols[i].maxz-center[2]);
	if( sx < fabs(n->pols[i].minz-center[2]))
	sx = fabs(n->pols[i].minz-center[2]);
	}
	*/
	n->scale = 1.0;// /sx;
	n->coord[0] = 0;//-center[0]*n->scale+0.5;
	n->coord[1] = 0;//-center[1]*n->scale+0.5;
	n->coord[2] = 0;//-center[2]*n->scale;

/*
	n->scale = 1./ sx;
	n->coord[0] = -center[0]*n->scale+0.5;
	n->coord[1] = -center[1]*n->scale+0.5;
	n->coord[2] = -center[2]*n->scale;*/

	PrepareMatricesObj(NULL,objs[cur].inverse,NULL,cur);
	return cur;
}




double DistToTreeNode(double x, double y, double z, struct tree * node)
{
	double ret[3];
	if( node == NULL )
		return 10e20;
	if( node->size == 0 )
		return 10e20;
	if( x > node->center[0]+node->side[0])
		ret[0] = x - (node->center[0]+node->side[0]);
	else if( x < node->center[0]-node->side[0] )
		ret[0] = (node->center[0]-node->side[0]) - x;
	else
		ret[0] = 0;
	if( y > node->center[1]+node->side[1])
		ret[1] = y - (node->center[1]+node->side[1]);
	else if( y < node->center[1]-node->side[1] )
		ret[1] = (node->center[1]-node->side[1]) - y;
	else
		ret[1] = 0;

	if( z > node->center[2]+node->side[2])
		ret[2] = z - (node->center[2]+node->side[2]);
	else if( z < node->center[2]-node->side[2] )
		ret[2] = (node->center[2]-node->side[2]) - z;
	else
		ret[2] = 0;
    /*
	ret[0] = x - node->center[0];
	ret[1] = y - node->center[1];
	ret[2] = z - node->center[2];
    */
	return sqrt(dotproduct(ret,ret));
}

double DistToTreePolygons(double x, double y, double z, struct tree * node, int objid)
{
	TRIANGLE tri;
	XYZ p;
	int i;
//    int k1,k2,hits2;
	SCALE_TYPE k,j;
	struct ObjInfo * t;
//	double a[3],b[3],nv[3], D, l, realdist;
    double temp,v[3],distmin =10e20,distmin2 = 10e28;
	double * dists;
	if( objid > nobj || objid < 0 ) return distmin;
	if( !objs[objid].used ) return distmin;
	if( node == NULL )
		return distmin;
	if( node->size == 0 )
		return distmin;
	dists = (double *)malloc(sizeof(double)*node->size);
    if( dists == NULL )
    {
        printf("Out of memory in DistToTreePolygons line %d\n Stop\n",__LINE__);
        exit(-1);
    }
	t = objs + objid;

//	printf("0\n");
	for(j = 0; j < node->size; j++)
	{
		if( x > node->set[j].coord[0]+node->set[j].radius[0])
			v[0] = x - (node->set[j].coord[0]+node->set[j].radius[0]);
		else if( x < node->set[j].coord[0]-node->set[j].radius[0] )
			v[0] = (node->set[j].coord[0]-node->set[j].radius[0]) - x;
		else
			v[0] = 0;
		if( y > node->set[j].coord[1]+node->set[j].radius[1])
			v[1] = y - (node->set[j].coord[1]+node->set[j].radius[1]);
		else if( y < node->set[j].coord[1]-node->set[j].radius[1] )
			v[1] = (node->set[j].coord[1]-node->set[j].radius[1]) - y;
		else
			v[1] = 0;

		if( z > node->set[j].coord[2]+node->set[j].radius[2])
			v[2] = z - (node->set[j].coord[2]+node->set[j].radius[2]);
		else if( z < node->set[j].coord[2]-node->set[j].radius[2] )
			v[2] = (node->set[j].coord[2]-node->set[j].radius[2]) - z;
		else
			v[2] = 0;
		dists[j] = sqrt(dotproduct(v,v));
		if( dists[j] < distmin2 )
			distmin2 = dists[j];
	}
    p.x= x;
    p.y= y;
    p.z =z;
	for(k = 0; k < node->size; k++)
	{

//		if( fabs(dists[k] - distmin2) > 0.000001 )
//			continue;
		i = node->set[k].polynum;
		if( t->pols[i].nvert < 3 )
			continue;

        for(j = 0; j < 3; j++)
        {
        	tri.p[j].x = ( t->verts[t->pols[i].vert[j]-1].v[0]);
	        tri.p[j].y = ( t->verts[t->pols[i].vert[j]-1].v[1]);
    	   	tri.p[j].z = ( t->verts[t->pols[i].vert[j]-1].v[2]);
        }
        temp = DistTriPoint(tri,p);//return!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if( temp < distmin ){
        	distmin = temp;

        }

	}
	free(dists);
	return distmin;
}

double DistToObjSub(double x, double y, double z, struct tree * node, int objid,int deep)
{
	int i;
	double dists[WIDTH],distmin = 10e20, distmin2=10e20,temp;
	if( objid > nobj || objid < 0 ) return 10e20;
	if( !objs[objid].used ) return 10e20;
	if( node->children == NULL)
	{
//		printf("dist to tree polygons\n");
		return DistToTreePolygons(x,y,z,node,objid);
	}
	for( i = 0; i < WIDTH; i++)
	{
//		printf("dist to tree node %d\n",i);
		dists[i] = DistToTreeNode(x,y,z,node->children+i);
		if( dists[i] < distmin)
			distmin = dists[i];
	}
	for( i = 0; i < WIDTH; i++)
	{
		if( fabs(distmin-dists[i]) > 0.000001 && fabs(distmin) > 0.000001  && deep > DEEP_START)
			continue;
//		printf("dist to tree polygons %d\n",i);
		temp = DistToObjSub(x,y,z,node->children+i,objid,deep+1);
		if( temp < distmin2 )
			distmin2 = temp;
	}
	return distmin2;
}

double DistToObj2(double x, double y, double z, int n)
{
	double a[3];
	if( n > nobj || n < 0 ) return 10e20;
	if( !objs[n].used ) return 10e20;
	double inverse[16];
    double matrix[16];
    double dist;
	PrepareMatricesObj(matrix,inverse,NULL,n);
	a[0] = x;
	a[1] = y;
	a[2] = z;
	vecmul(a,a,inverse);

	x = a[0];
	y = a[1];
	z = a[2];
    dist = DistToObjSub(x,y,z,objs[n].t,n,0)*objs[n].scale;
	return dist;
}

double DistToObj(double x, double y, double z, int n)
{
	TRIANGLE tri;
    XYZ p;
	int i,j;
	double distmin = 10e20;
	struct ObjInfo * t;
	double temp,a[3];


//    int hits2,k1,k2;
//	double b[3],v[3],nv[3],l,D,realdist;

	if( n > nobj || n < 0 ) return distmin;
//	if( !objs[n].used || !objs[n].active ) return distmin;
	double inverse[16];
	PrepareMatricesObj(NULL,inverse,NULL,n);
	a[0] = x;
	a[1] = y;
	a[2] = z;

	vecmul(a,a,inverse);

	x = a[0];
	y = a[1];
	z = a[2];


  	p.x = ( a[0] );
	p.y = ( a[1] );
	p.z = ( a[2] );

	t = objs + n;

	for( i = 0; i < t->npols; i++)
	{
		if( t->pols[i].nvert < 3 )
			continue;

        for(j = 0; j < 3; j++)
        {
        	tri.p[j].x = (t->verts[t->pols[i].vert[j]-1].v[0]);
	        tri.p[j].y = (t->verts[t->pols[i].vert[j]-1].v[1]);
    	    tri.p[j].z = (t->verts[t->pols[i].vert[j]-1].v[2]);
        }
        temp = DistTriPoint(tri,p);
        if( temp < distmin )
        	distmin = temp;
	}
	return distmin;
}
double DistToAnyObj(double x, double y, double z)
{
	double distmin = 10e20,t;
	int i;
	for(i = 0; i < nobj; i++){
		t = DistToObj(x,y,z,i);
		if( t < distmin )
			distmin = t;
	}
	return distmin;
}

double DistToAnyObj2(double x, double y, double z)
{
	double distmin = 10e20,t;
	int i;
	for(i = 0; i < nobj; i++)
		if( (t = DistToObj2(x,y,z,i)) < distmin )
			distmin = t;
	return distmin;
}


bool PointInAnyObj2(double x, double y, double z)
{
	int i;
	for(i = 0; i < nobj; i++)
            if( PointInObj2(x,y,z,i) )
                return true;
	return false;
}


void InitObj()
{
	int i;
	for(i = 0; i < OBJ_MAX; i++)
	{
		objs[i].verts = NULL;
		objs[i].norms = NULL;
		objs[i].pols = NULL;
		objs[i].used = false;
	}
}



void DeleteObj(int i)
{
	int j;
	if( i >= nobj || i < 0 ) return;
	if( !objs[i].used ) return;
    if( objs[i].verts != NULL )
    {
        free(objs[i].verts);
        objs[i].verts = NULL;
    }
    if( objs[i].norms != NULL )
    {
        free(objs[i].norms);
        objs[i].norms = NULL;
    }
    if( objs[i].texs != NULL )
    {
        free(objs[i].texs);
        objs[i].texs = NULL;
    }
    if( objs[i].mats != NULL )
    {
        free(objs[i].mats);
        objs[i].mats = NULL;
    }
    if( objs[i].pols != NULL )
    {
        for( j = 0; j < objs[i].npols; j++)
        {
            if(objs[i].pols[j].vert != NULL )
            {
                free(objs[i].pols[j].vert);
                objs[i].pols[j].vert = NULL;
            }
            if(objs[i].pols[j].norm != NULL )
            {
                free(objs[i].pols[j].norm);
                objs[i].pols[j].norm = NULL;
            }
            if(objs[i].pols[j].tex != NULL )
            {
                free(objs[i].pols[j].tex);
                objs[i].pols[j].tex = NULL;
            }
        }
        free(objs[i].pols);
        objs[i].pols = NULL;
    }
    objs[i].used = false;
//	printf("Free points mem\n");
//	printf("Free tree mem\n");
	kdtree_mem_destroy(objs[i].t);
	free(objs[i].p);
	free(objs[i].t);
//	free(objs[i].memcells);
//	printf("Ok\n");
//??????????????????????????????????????????????????????????????????????
	for( j = i; j < nobj; j++)
		objs[j] = objs[j+1];
	objs[nobj-1].used = false;
	nobj--;
//??????????????????????????????????????????????????????????????????????
}

void DestroyObj()
{
	while(nobj != 0)
		DeleteObj(0);
}

///////////////////////////////////////////////////////////////
void rayrnd(double ray[3])
{
	ray[0] = (rand()/(double)RAND_MAX-0.5)*2.0;
	ray[1] = (rand()/(double)RAND_MAX-0.5)*2.0;
	ray[2] = (rand()/(double)RAND_MAX-0.5)*2.0;
}

int raytri(double tri[3][3], double pos[3], double ray[3], double * res, double cutout)
{
    double a[3],b[3],c[3],n[3], d, k, m;
    double dot00,dot01,dot02,dot11,dot12,invdenom, uq,vq;
    a[0] = tri[0][0] - tri[2][0];
    a[1] = tri[0][1] - tri[2][1];
    a[2] = tri[0][2] - tri[2][2];
    b[0] = tri[1][0] - tri[2][0];
    b[1] = tri[1][1] - tri[2][1];
    b[2] = tri[1][2] - tri[2][2];
    crossproduct(a,b,n);
    d = -dotproduct(n,tri[0]);
    m =  dotproduct(n,ray);
    if( fabs(m) < 1.0e-25 )
        return 0;
    k = -(d + dotproduct(n,pos))/m;
	if( k > cutout ) return 0; 
	if( res != NULL ) *res = k;
    if( k < 0 )
        return 0;
    c[0] = pos[0] + k*ray[0] - tri[2][0];
    c[1] = pos[1] + k*ray[1] - tri[2][1];
    c[2] = pos[2] + k*ray[2] - tri[2][2];
    dot00 = dotproduct(a,a);
    dot01 = dotproduct(a,b);
    dot02 = dotproduct(a,c);
    dot11 = dotproduct(b,b);
    dot12 = dotproduct(b,c);
    invdenom = (dot00*dot11 - dot01*dot01);
	uq = (dot11*dot02-dot01*dot12);
	vq = (dot00*dot12-dot01*dot02);
	if( fabs(invdenom) < 1.0e-25  && fabs(uq) >= 0.0 && fabs(vq) >= 0.0 )
		return 0;
	uq = uq/invdenom;
	vq = vq/invdenom;
    if( uq >= -1.0e-12 && vq >= -1.0e-12 && 1.0-(uq+vq) >= -1.0e-12 )
        return 1;
    return 0;
}

#include <iostream>

bool raybox(double box[6], double pos[3], double ray[3], double cutout)
{
	double tnear = -1.0e20, tfar = 1.0e20;
	double t1,t2,c;
	double cc[3] = {0,0,0};
	for(int i = 0; i < 3; i++)
	{
		if( fabs(ray[i]) < 1.0e-7 )
		{
			if( pos[i] < box[i*2]-Eps || pos[i] > box[i*2+1]+Eps )
				return false;
			//cc[i] = 0;
		}
		else
		{
			t1 = (box[i*2+0]-Eps - pos[i])/ray[i];
			t2 = (box[i*2+1]+Eps - pos[i])/ray[i];
			if( t1 > t2 )
			{
				c = t1;
				t1 = t2;
				t2 = c;
			}
			if( t1 > tnear ) tnear = t1;
			if( t2 < tfar ) tfar = t2;
			if( tnear > tfar ) return false;
			if( tfar < 0 ) return false;
			//cc[i] = tnear*fabs(ray[i]);
//			if( cc[i] < 0 ) cc[i] = 0;
//			else if( cc[i] > cutout ) return false;
			if( tnear*fabs(ray[i]) > cutout ) return false;
		}
	}
	//if( sqrt(cc[0]*cc[0] + cc[1]*cc[1] + cc[2]*cc[2]) > 1.0 ) return false;//underflow arithmetic
	return true;
}



int IntersectTreeNodePolygons(double pos[3],double ray[3], struct tree * node, double * intersect,int & nintersect, int n, double cutout)
{
	int i,j;
	SCALE_TYPE k;
	struct ObjInfo * t;
	double box[6], tri[3][3];
	double q;
	if( node == NULL ) return 0;
	if( node->size == 0 ) return 0;
	if( n > nobj || n < 0 ) return 0;
	if( !objs[n].used ) return 0;


	t = objs + n;
	for(k = 0; k < node->size; k++)
	{
		for(i = 0; i < 3; i++)
		{
			box[i*2+0] = node->set[k].coord[i]-node->set[k].radius[i];
			box[i*2+1] = node->set[k].coord[i]+node->set[k].radius[i];
		}
		if( !raybox(box,pos,ray,cutout) ) continue;
		j = node->set[k].polynum;
		for(i = 0; i < 3; i++)
		{
			tri[i][0] = t->verts[t->pols[j].vert[i]-1].v[0];
			tri[i][1] = t->verts[t->pols[j].vert[i]-1].v[1];
			tri[i][2] = t->verts[t->pols[j].vert[i]-1].v[2];
		}
		if(raytri(tri,pos,ray,&q,cutout))
		{
			int flag = 0;
			//if( fabs(q) < 1.e-8 )
			//	return 1;//пересечение с границей
			for(i = 0; i < nintersect; i++)
				if( fabs(intersect[i]-q) < 1.e-9 )
				{
					flag = 1;
					break;
				}
			if( !flag )
				intersect[nintersect++] = q;
		}
	}
	return 0;
}


int IntersectTreeNodePolygonsNormal(double pos[3],double ray[3], struct tree * node, double * intersect, double * normal, int & nintersect, int n, double cutout)
{
	int i,j;
	SCALE_TYPE k;
	struct ObjInfo * t;
	double box[6], tri[3][3], v1[3],v2[3];
	double q;
	if( node == NULL ) return 0;
	if( node->size == 0 ) return 0;
	if( n > nobj || n < 0 ) return 0;
	if( !objs[n].used ) return 0;


	t = objs + n;
	for(k = 0; k < node->size; k++)
	{
		for(i = 0; i < 3; i++)
		{
			box[i*2+0] = node->set[k].coord[i]-node->set[k].radius[i];
			box[i*2+1] = node->set[k].coord[i]+node->set[k].radius[i];
		}
		if( !raybox(box,pos,ray,cutout) ) continue;
		j = node->set[k].polynum;
		for(i = 0; i < 3; i++)
		{
			tri[i][0] = t->verts[t->pols[j].vert[i]-1].v[0];
			tri[i][1] = t->verts[t->pols[j].vert[i]-1].v[1];
			tri[i][2] = t->verts[t->pols[j].vert[i]-1].v[2];
		}
		if(raytri(tri,pos,ray,&q,cutout))
		{
			int flag = 0;
			//if( fabs(q) < 1.e-8 )
			//	return 1;//пересечение с границей
			for(i = 0; i < nintersect; i++)
				if( fabs(intersect[i]-q) < 1.e-9 )
				{
					flag = 1;
					break;
				}
			if( !flag )
			{
				intersect[nintersect] = q;
				vecdiff(tri[0],tri[2],v1);
				vecdiff(tri[1],tri[2],v2);
				crossproduct(v1,v2,normal+nintersect*3);
				normalize(normal+nintersect*3);
				nintersect++;
			}
		}
	}
	return 0;
}


bool IntersectTreeNode(double pos[3], double ray[3], struct tree * node, double cutout)
{
	double box[6];
	if( node == NULL ) return false;
	if( node->size == 0 ) return false;


	for(int i = 0; i < 3; i++)
	{
		box[i*2+0] = node->center[i]-node->side[i];
		box[i*2+1] = node->center[i]+node->side[i];
	}
	return raybox(box,pos,ray,cutout);
}

int PointInObjSub(double pos[3], double ray[3], struct tree * node, double * intersect, int & nintersect, int n, double cutout)
{
	int i;
	//int sum = 0;
	if( node == NULL ) return 0;
	if( n > nobj || n < 0 ) return 0;
	if( !objs[n].used ) return 0;
	if( node->size == 0 ) return 0;
	if( node->children == NULL )
		return IntersectTreeNodePolygons(pos,ray,node,intersect,nintersect, n,cutout);


	for(i = 0; i < WIDTH; i++)
		if( IntersectTreeNode(pos,ray,node->children+i,cutout) )
			if( PointInObjSub(pos,ray,node->children+i,intersect,nintersect,n,cutout) )
				return 1;


	return 0;
}


int PointInObjSubNormal(double pos[3], double ray[3], struct tree * node, double * intersect, double * normal, int & nintersect, int n, double cutout)
{
	int i;
	//int sum = 0;
	if( node == NULL ) return 0;
	if( n > nobj || n < 0 ) return 0;
	if( !objs[n].used ) return 0;
	if( node->size == 0 ) return 0;
	if( node->children == NULL )
		return IntersectTreeNodePolygonsNormal(pos,ray,node,intersect,normal,nintersect, n, cutout);


	for(i = 0; i < WIDTH; i++)
		if( IntersectTreeNode(pos,ray,node->children+i,cutout) )
			if( PointInObjSubNormal(pos,ray,node->children+i,intersect,normal,nintersect,n, cutout) )
				return 1;


	return 0;
}


bool PointInObj2(double x, double y, double z, int n)
{
	int test1,test2;
	double vec[3];
	if( n > nobj || n < 0 ) return false;
	if( !objs[n].used ) return false;


	test1 = test2 = 0;


//	if( objs[n].npols < 1000 ) return PointInObj(x,y,z,n);


	vec[0] = x;
	vec[1] = y;
	vec[2] = z;


	vecmul(vec,vec,objs[n].inverse);


	x = vec[0];
	y = vec[1];
	z = vec[2];
	vec[0] = objs[n].maxx;
	vec[1] = objs[n].maxy;
	vec[2] = objs[n].maxz;
	if( vec[0] < x || vec[1] < y || vec[2] < z)
		return false;
	vec[0] = objs[n].minx;
	vec[1] = objs[n].miny;
	vec[2] = objs[n].minz;
	if( vec[0] > x || vec[1] > y || vec[2] > z)
		return false;
	vec[0] = x;
	vec[1] = y;
	vec[2] = z;
	{
		double ray[3];// = {1.,1.,1.};
		double intersect[1000];
		int nintersect = 0;//, nintersect2 = 0;
		srand(0);
		rayrnd(ray);
		if( PointInObjSub(vec,ray,objs[n].t,intersect,nintersect,n,1.0e+25) )
			return true;
		test1 = nintersect%2;
		/*
		rayrnd(ray);
		if( PointInObjSub(vec,ray,objs[n].t,intersect,nintersect2,n) )
			return true;
		test1 = (nintersect%2+nintersect2%2)/2;
		*/
		return test1;
	}

	//return (test1+test2)/2;
}



double RayObjIntersection(double pos[3], double ray[3], int obj, double cutout)
{
	double intersect[1024], min = 1e20;
	int nintersect = 0;
	PointInObjSub(pos,ray,objs[obj].t,intersect, nintersect,obj,cutout);
	if( nintersect )
	{
		for(int j = 0; j < nintersect; j++)
			if( intersect[j] < min )
				min = intersect[j];
	}
	//printf("intersection %g\n",min);
	return min;
}


double RayObjNormal(double pos[3], double ray[3], int obj, double * nrm, double cutout)
{
	double intersect[1024], min = 1e20;
	double normal[3072];
	int nintersect = 0;
	PointInObjSubNormal(pos,ray,objs[obj].t,intersect,normal, nintersect,obj, cutout);
	if( nintersect )
	{
		for(int j = 0; j < nintersect; j++)
			if( intersect[j] < min )
			{
				min = intersect[j];
				nrm[0] = normal[j*3+0];
				nrm[1] = normal[j*3+1];
				nrm[2] = normal[j*3+2];
			}
	}
	//printf("intersection %g\n",min);
	return min;
}


double MiddleZ(int n)
{
	double z = 0;
	for(int i = 0; i < objs[n].nverts; i++)
		z += objs[n].verts[i].v[2];
	return z / (double) objs[n].nverts;
	
}
#if defined(__GRAPHICS__)
#include "my_glut.h"
inline void copy(double * veca, double * vecb)
{
	veca[0] = vecb[0];
	veca[1] = vecb[1];
	veca[2] = vecb[2];
}
void DrawObj(int n, int state, int highlight)
{
	int i,j,k,m;
	if( n > nobj || n < 0 ) return;
	if( !objs[n].used ) return;

    if( state > 1)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		/*
		if( highlight )
		{
			glColor3f(0.3,0.1,0.1);
		}
		else
		{
			glColor3f(0.1,0.1,0.1);
		}
		*/
        for(i = 0; i < objs[n].npols; i++)
        {
            glBegin(GL_POLYGON);
            for(j = 0; j < objs[n].pols[i].nvert; j++)
            {
                double v[3],nn[3];
                k = objs[n].pols[i].vert[j]-1;
                if ( j < objs[n].pols[i].nnorm )
                    m = objs[n].pols[i].norm[j]-1;
                else
                    m = -1;
                
                if( m >= 0 && m < objs[n].nnorms)
                {
                    copy(nn,objs[n].norms[m].n);
                    //vecmul(nn,nn,objs[n].rotonly);
                    glNormal3dv(nn);
                }
                copy(v,objs[n].verts[k].v);
                //vecmul(v,v,objs[n].matrix);
                glVertex3d(v[0],v[1],v[2]);
            }
            glEnd();
        }
    }

    if( state > 0  && state != 3)
    {
		if( highlight )
			glColor3f(0.0f,1.0f,1.0f);
		else 
			glColor3f(1.0f,1.0f,1.0f);
        glLineWidth(2.0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        for(i = 0; i < objs[n].npols; i++)
        {
            glBegin(GL_POLYGON);
            for(j = 0; j < objs[n].pols[i].nvert; j++)
            {
                double v[3],nn[3];
                k = objs[n].pols[i].vert[j]-1;
                if ( j < objs[n].pols[i].nnorm )
                    m = objs[n].pols[i].norm[j]-1;
                else
                    m = -1;
                
                if( m >= 0 && m < objs[n].nnorms)
                {
                    copy(nn,objs[n].norms[m].n);
                    //vecmul(nn,nn,objs[n].rotonly);
                    glNormal3dv(nn);
                }
                copy(v,objs[n].verts[k].v);
                //vecmul(v,v,objs[n].matrix);
                glVertex3d(v[0],v[1],v[2]);
            }
            glEnd();
        }
        glLineWidth(1.0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
	//DrawTreeObj(objs[n].t,0,identity);
}
#endif
