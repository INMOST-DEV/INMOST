#pragma once
#ifndef _DISCR_H

#include "../../inmost.h"
using namespace INMOST;

//TODO:
//move orphan function inside class
//introduce discr_triplet_basic class 

void tensor_prod1(Storage::real * K, Storage::real v[3], Storage::real out[3]);
void tensor_prod2(Storage::real * K, Storage::real v[3], Storage::real out[3]);
void tensor_prod3(Storage::real * K, Storage::real v[3], Storage::real out[3]);
int tensor_prod(int Ktype, Storage::real_array & K, 
				Storage::real v[3], Storage::real out[3]);
void tensor_add1(Storage::real * K, Storage::real Kout[9]);
void tensor_add2(Storage::real * K, Storage::real Kout[9]);
void tensor_add3(Storage::real * K, Storage::real Kout[9]);
int tensor_add(int Ktype, Storage::real_array & K, Storage::real Kout[9]);
int average_tensor(int Ktype, Storage::real_array & K1, 
					Storage::real_array & K2, Storage::real Kout[9]);
void vec_diff(Storage::real x[3], Storage::real y[3], Storage::real out[3]);
Storage::real dot_prod(Storage::real x[3], Storage::real y[3]);


//type of the averaging point
#define AVG_NONE        0//material is homogeneous across cells - no point required
#define AVG_CENTROID    1//use centroid
#define AVG_HARMONIC    2//use harmonic averaging point
#define AVG_NONLINEAR   4//use nonlinear flux continuety to define value at the point
#define AVG_REVDIST     8//use reverse distance
#define AVG_NEUMANN     16 //use zero neumann boundary condition

#define CRITICAL_DIST_RATIO 10.0

class discr_basic
{
protected:
	Automatizator * aut;
	Mesh * m;
	Tag K;
	Tag bnd_conds;
	//Tag AvgCoord, AvgCoefs, AvgStats;
	int Ktype;
	int error;
	ElementType bnd_types;
	MIDType bnd_markers, add_markers;
	int rem_marker;
	int have_bnd;
	int allowed_types;
	// returns
	// 0 if point inside
	// 1 if distance to point by distance to farthest node is less then CRITICAL_DIST_RATIO
	// 2 otherwise
	int isPointInside(Face * f, Storage::real x[3]);
	//void InitializeSupportPoints(int allowed_types = AVG_NONE | AVG_CENTROID | AVG_HARMONIC | AVG_NONLINEAR | AVG_REVDIST);
	int FindHarmonicPoint(Face * fKL, Cell * cK, Cell * cL, 
						   Tag tensor, int tensor_type, 
						   Storage::real xK[3], Storage::real xL[3], 
						   Storage::real y[3], Storage::real & coef);
	int FindBoundaryPoint(Face * fK, Cell * cK, 
		                   Tag tensor, int tensor_type, 
						   Storage::real xK[3], Storage::real y[3]);
	virtual void Build(Face * f) = 0;
public:
	enum DiscrType {TPFA, MPFA, NTPFA, NMPFA};
	enum SchemeType{ Linear, Nonlinear };
public:
	discr_basic(Automatizator * aut, Mesh * m, std::string tensor_name);
	virtual bool needBuild(Face * f) const;
	virtual bool needUpdate(Face * f) const;
	virtual void Init() = 0;
	virtual DiscrType GetDiscrType() const = 0;
	virtual SchemeType GetSchemeType() const = 0;
	virtual expr Grad(const expr & param) const = 0;
	virtual expr BoundaryGrad(const expr & param) const = 0;
	virtual expr Interp(ElementType etype,const expr & param) const = 0;
	virtual void Update() = 0;
	virtual void Export(std::ostream & out, Storage::real trans_scale, Storage::real vol_scale) const; //export to ADGPRS
	// by default no boundary information is provided, schemes avoid using boundary
	// first marker is the marker that is used inside Automatizator to define dynamic tag domain
	// if some boundary element is not used inside scheme it will not be marked by usage_marker
	// second marker defines boundary faces/edges, if not set, it is computed automatically
	// third parameter defines what kinds of elements may be used
	void Boundary(MIDType bnd_marker = 0, ElementType bnd_types = FACE | EDGE);
	// mark additional collocation points that should be introduced into the system
	MIDType UnknownMarker() const { return add_markers; }
	MIDType BoundaryMarker() const { return bnd_markers; }
	virtual ~discr_basic();

	void AllowTypes(int types) {allowed_types = types;}
};

class discr_triplet_basic : public discr_basic
{
protected:
	int find_triplet(Element * r, Element * f, Face * harmonic_face, Storage::real v[3],
					 Storage::real dir, Storage::real area, 
					 Element * t[3], Storage::real coefs[3], int ptypes[3], 
					 int allowed = AVG_HARMONIC | AVG_NONLINEAR | AVG_REVDIST, bool f_necessery = false);
public:
	discr_triplet_basic(Automatizator * aut, Mesh * m, std::string tensor_name) : discr_basic(aut,m,tensor_name) {}
	virtual ~discr_triplet_basic();
};


class tpfa : public discr_basic
{
	INMOST_DATA_ENUM_TYPE stncl;
	Tag trans;
	Tag elems;
	void Build(Face * f);
public:
	tpfa(Automatizator * aut, Mesh *m, std::string tensor);
	DiscrType GetDiscrType() const { return TPFA; }
	SchemeType GetSchemeType() const { return Linear; }
	expr Grad(const expr & param) const;
	expr BoundaryGrad(const expr & param) const {return Grad(param);}
	expr Interp(ElementType etype,const expr & param) const;
	void Export(std::ostream & out, Storage::real trans_scale, Storage::real vol_scale) const;
	void Init();
	void Update();
	bool needBuild(Face * f) const { return discr_basic::needBuild(f); }
	bool needUpdate(Face * f) const { return discr_basic::needUpdate(f); }
	~tpfa();
};

class tpfa2 : public discr_basic
{
	INMOST_DATA_ENUM_TYPE stncl;
	Tag trans;
	Tag elems;
	void Build(Face * f);
public:
	tpfa2(Automatizator * aut, Mesh *m, std::string tensor);
	DiscrType GetDiscrType() const { return TPFA; }
	SchemeType GetSchemeType() const { return Linear; }
	expr Grad(const expr & param) const;
	expr BoundaryGrad(const expr & param) const {return Grad(param);}
	expr Interp(ElementType etype,const expr & param) const;
	void Export(std::ostream & out, Storage::real trans_scale, Storage::real vol_scale) const;
	void Init();
	void Update();
	bool needBuild(Face * f) const { return discr_basic::needBuild(f); }
	bool needUpdate(Face * f) const { return discr_basic::needUpdate(f); }
	~tpfa2();
};



class mpfa_o : public discr_basic
{
	INMOST_DATA_ENUM_TYPE stncl;
	Tag trans;
	Tag elems;
	void Build(Face * f);
public:
	mpfa_o(Automatizator * out, Mesh *m, std::string tensor);
	DiscrType GetDiscrType() const { return MPFA; }
	SchemeType GetSchemeType() const { return Linear; }
	expr Grad(const expr & param) const;
	expr BoundaryGrad(const expr & param) const {return Grad(param);}
	expr Interp(ElementType etype,const expr & param) const;
	void Export(std::ostream & out, Storage::real trans_scale, Storage::real vol_scale) const;
	void Init();
	void Update();
	bool needBuild(Face * f) const { return discr_basic::needBuild(f); } //TODO
	bool needUpdate(Face * f) const { return discr_basic::needBuild(f); } //TODO
	~mpfa_o();
};

class basic_a : public discr_triplet_basic
{
	MIDType mark_next;
	INMOST_DATA_ENUM_TYPE dwn, upw, intrp, add;
	INMOST_DATA_ENUM_TYPE trpr;
	INMOST_DATA_ENUM_TYPE trpl;
	Tag trpl_elems, trpr_elems;
	Tag trpl_coefs, trpr_coefs;
	Tag trans;
	Tag intrp_stncl, intrp_coefs;
	Tag pos;
	void MarkRecursive(Element * f);
	void Build(Face * f);
	static void get_r(Storage * e, Automatizator::stencil_pairs & out, void * user_data);
	static void get_l(Storage * e, Automatizator::stencil_pairs & out, void * user_data);
public:
	basic_a(Automatizator * out, Mesh *m, std::string tensor);
	void Init();
	void Update();
	bool needBuild(Face * f) const { return discr_basic::needBuild(f); } //TODO
	bool needUpdate(Face * f) const { return discr_basic::needBuild(f); } //TODO
	virtual ~basic_a();
};

class ntpfa_a : public discr_triplet_basic
{
	MIDType mark_next;
	INMOST_DATA_ENUM_TYPE dwn, upw, intrp;
	INMOST_DATA_ENUM_TYPE trpr, addr;
	INMOST_DATA_ENUM_TYPE trpl, addl;
	Tag trpl_elems, trpr_elems;
	Tag trpl_coefs, trpr_coefs; 
	Tag trpl_add, trpr_add;
	Tag intrp_stncl, intrp_coefs;
	Tag pos;
	void MarkRecursive(Element * f);
	void MarkRecursiveEdge(Element * f);
	void Build(Face * f);
	static void get_r(Storage * e, Automatizator::stencil_pairs & out, void * user_data);
	static void get_l(Storage * e, Automatizator::stencil_pairs & out, void * user_data);
public:
	ntpfa_a(Automatizator * out, Mesh *m, std::string tensor);
	DiscrType GetDiscrType() const { return NTPFA; }
	SchemeType GetSchemeType() const { return Nonlinear; }
	expr Grad(const expr & param) const;
	expr BoundaryGrad(const expr & param) const {return Grad(param);}
	expr Interp(ElementType etype,const expr & param) const;
	void Export(std::ostream & out, Storage::real trans_scale, Storage::real vol_scale) const;
	void Init();
	void Update();
	bool needBuild(Face * f) const { return discr_basic::needBuild(f); } //TODO
	bool needUpdate(Face * f) const { return discr_basic::needBuild(f); } //TODO
	~ntpfa_a();
};

class nmpfa_a : public basic_a
{
public:
	nmpfa_a(Automatizator * out, Mesh *m, std::string tensor);
	DiscrType GetDiscrType() const { return NMPFA; }
	SchemeType GetSchemeType() const { return Nonlinear; }
	expr Grad(const expr & param) const;
	expr BoundaryGrad(const expr & param) const {return Grad(param);}
	expr Interp(ElementType etype,const expr & param) const;
	void Export(std::ostream & out, Storage::real trans_scale, Storage::real vol_scale) const;
	~nmpfa_a();
};


class ntpfa_c : public discr_triplet_basic
{
	MIDType mark_next;
	INMOST_DATA_ENUM_TYPE dwn, upw, intrp;
	INMOST_DATA_ENUM_TYPE trpr1, trpl1, addl1, addr1;
	INMOST_DATA_ENUM_TYPE trpr2, trpl2, addl2, addr2;
	Tag trpl1_elems, trpl2_elems, trpr1_elems, trpr2_elems;
	Tag trpl1_coefs, trpl2_coefs, trpr1_coefs, trpr2_coefs; 
	Tag trpl1_add, trpl2_add, trpr1_add, trpr2_add;
	Tag intrp_stncl, intrp_coefs;
	void MarkRecursive(Element * f);
	void MarkRecursiveEdge(Element * f);
	void Build(Face * f);
	static void get_r(Storage * e, Automatizator::stencil_pairs & out, void * user_data);
	static void get_l(Storage * e, Automatizator::stencil_pairs & out, void * user_data);
public:
	ntpfa_c(Automatizator * out, Mesh *m, std::string tensor);
	DiscrType GetDiscrType() const { return NTPFA; }
	SchemeType GetSchemeType() const { return Nonlinear; }
	expr Grad(const expr & param) const;
	expr BoundaryGrad(const expr & param) const;
	expr Interp(ElementType etype,const expr & param) const;
	void Export(std::ostream & out, Storage::real trans_scale, Storage::real vol_scale) const;
	void Init();
	void Update();
	bool needBuild(Face * f) const { return discr_basic::needBuild(f); } //TODO
	bool needUpdate(Face * f) const { return discr_basic::needBuild(f); } //TODO
	~ntpfa_c();
};

class ntpfa_b : public discr_triplet_basic
{
	MIDType mark_next;
	INMOST_DATA_ENUM_TYPE dwn, upw, intrp;
	INMOST_DATA_ENUM_TYPE trpr, addr;
	INMOST_DATA_ENUM_TYPE trpl, addl;
	Tag trpl_elems, trpr_elems;
	Tag trpl_coefs, trpr_coefs; 
	Tag trpl_add, trpr_add;
	Tag intrp_stncl, intrp_coefs;
	Tag pos;
	void MarkRecursive(Element * f);
	void MarkRecursiveEdge(Element * f);
	void Build(Face * f);
	static void get_r(Storage * e, Automatizator::stencil_pairs & out, void * user_data);
	static void get_l(Storage * e, Automatizator::stencil_pairs & out, void * user_data);
public:
	ntpfa_b(Automatizator * out, Mesh *m, std::string tensor);
	DiscrType GetDiscrType() const { return NTPFA; }
	SchemeType GetSchemeType() const { return Nonlinear; }
	expr Grad(const expr & param) const;
	expr BoundaryGrad(const expr & param) const {return Grad(param);}
	expr Interp(ElementType etype,const expr & param) const;
	void Export(std::ostream & out, Storage::real trans_scale, Storage::real vol_scale) const;
	void Init();
	void Update();
	bool needBuild(Face * f) const { return discr_basic::needBuild(f); } //TODO
	bool needUpdate(Face * f) const { return discr_basic::needBuild(f); } //TODO
	~ntpfa_b();
};



class nmpfa : public discr_triplet_basic
{
	MIDType mark_next;
	INMOST_DATA_ENUM_TYPE dwn, upw, intrp, add;
	INMOST_DATA_ENUM_TYPE trpr;
	INMOST_DATA_ENUM_TYPE trpl;
	Tag trpl_elems, trpr_elems;
	Tag trpl_coefs, trpr_coefs;
	Tag addval;
	Tag intrp_stncl, intrp_coefs;
	Tag pos;
	void MarkRecursive(Element * f);
	void Build(Face * f);
	static void get_r(Storage * e, Automatizator::stencil_pairs & out, void * user_data);
	static void get_l(Storage * e, Automatizator::stencil_pairs & out, void * user_data);
public:
	nmpfa(Automatizator * out, Mesh *m, std::string tensor);
	DiscrType GetDiscrType() const { return NMPFA; }
	SchemeType GetSchemeType() const { return Nonlinear; }
	expr Grad(const expr & param) const;
	expr BoundaryGrad(const expr & param) const {return Interp(FACE,param);}
	expr Interp(ElementType etype,const expr & param) const;
	void Export(std::ostream & out, Storage::real trans_scale, Storage::real vol_scale) const;
	void Init();
	void Update();
	bool needBuild(Face * f) const { return discr_basic::needBuild(f); } //TODO
	bool needUpdate(Face * f) const { return discr_basic::needBuild(f); } //TODO
	~nmpfa();
};



class mpfa_t : public discr_triplet_basic
{
	MIDType mark_next;
	INMOST_DATA_ENUM_TYPE scheme, interp;
	Tag stncl, coefs, interp_stncl, interp_coefs;
	void MarkRecursive(Element * f);
	void Build(Face * f);
public:
	mpfa_t(Automatizator * out, Mesh *m, std::string tensor);
	DiscrType GetDiscrType() const { return NMPFA; }
	SchemeType GetSchemeType() const { return Linear; }
	expr Grad(const expr & param) const;
	expr BoundaryGrad(const expr & param) const {return Grad(param);}
	expr Interp(ElementType etype,const expr & param) const;
	void Export(std::ostream & out, Storage::real trans_scale, Storage::real vol_scale) const;
	void Init();
	void Update();
	bool needBuild(Face * f) const { return discr_basic::needBuild(f); } //TODO
	bool needUpdate(Face * f) const { return discr_basic::needBuild(f); } //TODO
	~mpfa_t();
};

class mpfa_y : public discr_triplet_basic
{
	MIDType mark_next;
	INMOST_DATA_ENUM_TYPE scheme, interp;
	Tag stncl, coefs, interp_stncl, interp_coefs;
	void MarkRecursive(Element * f);
	void Build(Face * f);
public:
	mpfa_y(Automatizator * out, Mesh *m, std::string tensor);
	DiscrType GetDiscrType() const { return NMPFA; }
	SchemeType GetSchemeType() const { return Linear; }
	expr Grad(const expr & param) const;
	expr BoundaryGrad(const expr & param) const {return Grad(param);}
	expr Interp(ElementType etype,const expr & param) const;
	void Export(std::ostream & out, Storage::real trans_scale, Storage::real vol_scale) const;
	void Init();
	void Update();
	bool needBuild(Face * f) const { return discr_basic::needBuild(f); } //TODO
	bool needUpdate(Face * f) const { return discr_basic::needBuild(f); } //TODO
	~mpfa_y();
};

#endif