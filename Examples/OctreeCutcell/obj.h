
#ifndef OBJ_H_
#define OBJ_H_

#ifndef _OBJ_H
#define _OBJ_H

// Maximum for objects array
#define OBJ_MAX 100

struct quat
{
	double vec[3];
	double w;
};
void quatprint(quat a);
quat quatmul(quat a, quat b);
quat quatfromaxis(double x, double y, double z, double ang);
quat quatdiv(quat a, double b);
quat quatconj(quat a);
void quatident(quat * a);
void quatzero(quat * a);
quat quatrot(quat a, quat b);
double quatnorm(quat a);
void quatmat(quat q, double * rot);
quat quatinv(quat a);
void matmul(double * a, double * b);

int ReadObj(char * file);
void InitObj();

void DeleteObj(int i);
void DestroyObj();

void TranslateObj(double dx, double dy, double dz,int n);
void RotateObj(double rx, double ry, double rz,int n);
void ScaleObj(double s,int n);

void SetRotateObj(double x, double y, double z, int n);
void SetScaleObj(double s,int n);
void SetTranslateObj(double x, double y, double z,int n);

int NumberofObj();


double DistToAnyObj(double x, double y, double z);
double DistToObj(double x, double y, double z, int n);
double DistToAnyObj2(double x, double y, double z);
double DistToObj2(double x, double y, double z, int n);

bool PointInObj2(double x, double y, double z, int n);
bool PointInAnyObj2(double x, double y, double z);



double RayObjNormal(double pos[3], double ray[3], int obj, double * nrm,double cutout);
double RayObjIntersection(double pos[3], double ray[3], int obj,double cutout);
double MiddleZ(int n);


void DrawObj(int n, int state, int highlight);

#endif


#endif /* OBJ_H_ */
