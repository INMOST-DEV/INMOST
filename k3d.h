//------------------------------------------------------------------------------------------------
// File: k3d.h
//------------------------------------------------------------------------------------------------

#ifndef __K3D_H
#define __K3D_H

#include <iostream>
#include <fstream>
//#include <cstring>
//#include <string>
#include <cmath>
#include <vector>

using namespace std;

namespace k3d {

struct SSolverParams // Class that contains solver parameters
{
   int prec_float; ///< The type of preconditioner floating data to be used: prec_float == 1 (default) - use float type, prec_float == 2 - use double type
   int ncycle;     ///< The number of sparsity extension cycles when determining backward overlap for blocks, 0 < ncycle, default value ncycle = 1
   int ordtype;    ///< The ordering type when performing extended factorization: ordtype == 0 - no ordering, ordtype = 1 (default) - RCM ordering for initial, main and separator parts
   int collap;     ///< The parameter that controls print sparsity condensing: collap == -1 (defalut) - no print, 0 < collap - print sparsity as .ps files for each block with condensing into collap times
   int sctype;     ///< The type of scaling: sctype == -1 - no scaling, sctype == 0 (default) - usual scaling via diagonal value, sctype == 1 - balancing scaling is performed via rows/colums norms computations
   int nitersc;    ///< The number of balancing scaling iterations, 0 < nitersc, default value nitersc = 3
   int fcttype;    ///< The structural factorization type parameter, fcttype == -1 - no structural restrictions onto sparsity (default), 0 <= fcttype - elements of order fcttype at most are remaind in fct
   double pivmin;  ///< The minimal pivot value, -1 < pivmin < 1.0, for negative value pivoting is ignored, default value pivmin = -1.0e0
   double tau1;    ///< First  threshold value, 0 <= \tau1 < 1.0, default value \tau1 = 0.01e0
   double tau2;    ///< Second threshold value, 0 <= \tau2 <= \tau1, default value \tau2 = 0.0001e0
   double theta;   ///< Diagonal correction parameter, 0 <= \theta <= 1.0, default value \theta = 0.1e0
// Constructors and destructor
      /// Epmty constructor
   SSolverParams () {
      prec_float = 1; ncycle = 1; ordtype = 1; collap = -1; sctype = 0; nitersc = 3; 
      fcttype = -1; pivmin = -1.0e0; tau1 = 1.0e-2; tau2 = 1.0e-4; theta = 0.1e0;
   }; 
      /// Set defaluts
   void SetDefaults () {
      prec_float = 1; ncycle = 1; ordtype = 1; collap = -1; sctype = 0; nitersc = 3; 
      fcttype = -1; pivmin = -1.0e0; tau1 = 1.0e-2; tau2 = 1.0e-4; theta = 0.1e0;
   };
};

   /// Perform binary search
int BinarySearch (long long _isup, int _nblks, long long *_blks, int _iblkprev);

   /// Print array
void PrintArray (ostream &_fout, const char *_name, int _isize, int *_iarr);
void PrintArray (ostream &_fout, const char *_name, int _isize, long long *_iarr);
void PrintArray (ostream &_fout, const char *_name, int _isize, float *_farr);
void PrintArray (ostream &_fout, const char *_name, int _isize, double *_farr);

ostream &SetPw (ostream &stream); ///< Set ouput manip
void Round (double &dx, double &dy); ///< Round values

class CExchange
{
// Data
public:
// Control functions
   static void Synchronize (void *_comm); /// Barrier synchronization of the CPUs
   static int GetNproc (void *_comm);     /// Get the number of CPUs
   static int GetMyid (void *_comm);      /// Get the ID number of the CPU
// Exchange functions
      /// Collective array operation
   static void ExchangeArray (void *_comm, char _Datatype, char _Operation,
                              int _Length, void *_Array);
      /// Collective array operation (implementation)
   static void ExchangeArray_impl (void *_comm, char _Datatype, char _Operation,
                              int _Length, void *_Array);
      /// Major data exchange function
   static void DataExchange (void *_comm, vector<int> &_CpuIDSend, vector<vector<char> > &_ObjSend,
                              vector<int> &_CpuIDRecv, vector<vector<char> > &_ObjRecv);
      /// Major data exchange function (implementation)
   static void DataExchange_impl (void *_comm, vector<int> &_CpuIDSend, vector<vector<char> > &_ObjSend,
                                    vector<int> &_CpuIDRecv, vector<vector<char> > &_ObjRecv);
      /// Data exchange function (implementation)
   static void DataExchange_impl2 (void *_comm,
                                    int _NObjSend, int *_CpuIDSend, int *_ObjSizeSend, char **_ObjSend,
                                    int &_NObjRecv, int *&_CpuIDRecv, int *&_ObjSizeRecv, char **&_ObjRecv);
      /// Exchange data requests
   static void DataExchangeRequest (void *_comm,
                                    int _NObjSend, int *_CpuIDSend, int *_ObjSizeSend,
                                    int &_NObjRecv, int *&_CpuIDRecv, int *&_ObjSizeRecv);
// Send/recv/test/wait/probe functions
      /// Asynchronous send data to the other CPU
   static void ISend (void *_comm, int _rank, int _msgtag, int _length, char *_buffer, int _indrecv, void *_recv_arr);
      /// Asynchronous receive data from the other CPU
   static void IRecv (void *_comm, int _rank, int _msgtag, int _length, char *_buffer, int _indrecv, void *_recv_arr);
      /// Wait for completion of exchanges
   static void WaitAll  (int _count, void *_recvarr, void *_statarr);
      /// Allocate send/receive requests
   static void AllocateRecvs (int _nrecv, void *&_recvarr);
      /// Allocate statuses
   static void AllocateStats (int _nstat, void *&_statarr);
      /// Delete send/receive requests
   static void DeleteRecvs (void *_recvarr);
      /// Delete statuses
   static void DeleteStats (void *_statarr);
};

class CSortInt /// Support sorting of int data
{
public:
// Data
   int ival;  // ivalue contains integer value to be sorted
   int i2val; // i2val contains secondary integer value to be sorted
// Functions
// Constructors and destructor
   CSortInt () {}; // Default constructor
   CSortInt (int _ival, int _i2val) { // Memory allocation constructor with zero data
      ival = _ival;
      i2val = _i2val;
   };
   ~CSortInt () {}; // Destructor
// Operator functions
   CSortInt &operator= (const CSortInt &_ii2) { // Equality operator
      ival = _ii2.ival;
      i2val = _ii2.i2val;
      return *this;
   };
// Comparison function
   inline bool operator< (const CSortInt& _obj) const {
      return ival < _obj.ival;
   };
};

class CSortInt2 /// Support sorting of double int data
{
public:
// Data
   int ixval; // ixval contains x integer value to be sorted
   int iyval; // iyval contains y integer value to be sorted
   int itail; // itail contains secondary integer value to be sorted
// Functions
// Constructors and destructor
   CSortInt2 () {}; // Default constructor
   CSortInt2 (int _ixval, int _iyval, int _itail) { // Memory allocation constructor with zero data
      ixval = _ixval;
      iyval = _iyval;
      itail = _itail;
   };
   ~CSortInt2 () {}; // Destructor
// Operator functions
   CSortInt2 &operator= (const CSortInt2 &_ii2) { // Equality operator
      ixval = _ii2.ixval;
      iyval = _ii2.iyval;
      itail = _ii2.itail;
      return *this;
   };
// Comparison function
   inline bool operator< (const CSortInt2& _obj) const {
      if(ixval < _obj.ixval)
         return true;
      else if(ixval > _obj.ixval)
         return false;
      return iyval < _obj.iyval;
   };
};

template <typename _FltVect>
class CVect /// Class that supports vector operations
{
public:
   static void SetByZeroes (int _n, _FltVect *_x); ///< Set vector data by zeroes
   static _FltVect ScProd (int _n, _FltVect *_x, _FltVect *_y); ///< Compute scalar product
   static void CopyVector (int _n, _FltVect *_x, _FltVect *_y); ///< Copy array
   static void AddReplaceVector (int _n, _FltVect *_x1, _FltVect *_x1plus2); ///< Add vector data
   static void SubtractReplaceVector (int _n, _FltVect *_x1, _FltVect *_x2minus1); ///< Subtract vector data
   static void UpdateVector (int _n, _FltVect *_value, _FltVect *_arr_x, _FltVect *_arr_y); ///< Update array (axpy): y=a*x+y
   static void UpdateVectorMinus (int _n, _FltVect *_value, _FltVect *_arr_x, _FltVect *_arr_y); ///< Update array for minus alpha (axpy)
   static void UpdateVectorReversed (int _n, _FltVect *_value, _FltVect *_arr_x, _FltVect *_arr_y); ///< Update array reversed (aypx)
   static void OrderVector (int _n, const vector<int> &_order, 
                              const _FltVect *_x, _FltVect *_x_ord); ///< Order vector data
   static void InvOrderVector (int _n, const vector<int> &_order, 
                              const _FltVect *_x, _FltVect *_x_ord); ///< Inverse order vector data
};

template <typename _Int, typename _Flt>
class CIlu2_impl // Class that supports factorization computations
{
public:
      /// Compute optimal ordering for the matrix
   static void ComputeOptimalOrder (int _ordtype,
                                    const int _n, const vector<_Int> &_ia_alu, const vector<_Int> &_ja_alu,
                                    vector<int> &_order);
      /// Compute optimal ordering for the matrix via splitting into 3 blocks: 1-st order Schur complement for third block, which is the separator between first block and remaining data
   static void ComputeOptimalOrderSchur (int _ordtype, 
                                          const int _n, int _n1, const vector<_Int> &_ia_alu, const vector<_Int> &_ja_alu,
                                          int &_n2, vector<int> &_order);
      /// Compute ordered matrix (sparsity only)
   static void ReorderMatrixSp (int _n, const vector<int> &_order,
                              vector<_Int> &_ia_alu, vector<_Int> &_ja_alu,
                              vector<_Int> &_ia_alu_ord, vector<_Int> &_ja_alu_ord);
      /// Compute ordered matrix
   static void ReorderMatrix (int _n, const vector<int> &_order,
                              vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                              vector<_Int> &_ia_alu_ord, vector<_Int> &_ja_alu_ord, vector<_Flt> &_a_alu_ord);
      /// Perform ILU2 point factorization of the block with future diagonal modification
   static void Ilu2BlockTransform (int _sctype, int _nitersc, int _fcttype, double _pivmin, double _tau1, double _tau2, double _theta,
                                    int _n, 
                                    vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                                    vector<_Int> &_ia_lu, vector<_Int> &_ida_lu, vector<_Int> &_ja_lu, vector<_Flt> &_a_lu,
                                    double &_eigmin_att, double &_eigmax_att);
      /// Perform ILU2 point factorization of the block with future diagonal modification
   static void Ilu2Block (int _sctype, int _nitersc, int _fcttype, double _pivmin, double _tau1, double _tau2, double _theta,
                           int _n, 
                           vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                           vector<_Int> &_ia_l, vector<_Int> &_ja_l, vector<_Flt> &_a_l,
                           vector<_Int> &_ia_u, vector<_Int> &_ja_u, vector<_Flt> &_a_u,
                           double &_eigmin_att, double &_eigmax_att);
      /// Compute scaling
   static void ComputeScaling (int _sctype, int _nitersc,
                                 int _n, 
                                 vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                                 vector<_Flt> &_sclL, vector<_Flt> &_sclU);
      /// Perform explicit scaling
   static void MatrixScale (int _n, vector<_Flt> &_sclL, vector<_Flt> &_sclU,
                              vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu);
      /// Symmetrize sparsity
   static void SymmetrizeSparsity (int _n, const vector<_Int> &_ia_alu, const vector<_Int> &_ja_alu,
                                    vector<_Int> &_ia_symm, vector<_Int> &_ja_symm);
      /// Split matrix data into L and U parts for sparsity only
   static void SplitLUSp (int _n, 
                           const vector<_Int> &_ia_alu, const vector<_Int> &_ja_alu,
                           vector<_Int> &_ia_l, vector<_Int> &_ja_l,
                           vector<_Int> &_ia_u, vector<_Int> &_ja_u);
      /// Split matrix data into L and U parts
   static void SplitLU (int _n, 
                        vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                        vector<_Int> &_ia_l, vector<_Int> &_ja_l, vector<_Flt> &_a_l,
                        vector<_Int> &_ia_u, vector<_Int> &_ja_u, vector<_Flt> &_a_u);
      /// Transpose square matrix for sparsity only
   static void TransposeSp (int _n, 
                           vector<_Int> &_ia_a, vector<_Int> &_ja_a,
                           vector<_Int> &_ia_at, vector<_Int> &_ja_at);
      /// Transpose square matrix
   static void Transpose (int _n, 
                           vector<_Int> &_ia_a, vector<_Int> &_ja_a, vector<_Flt> &_a_a,
                           vector<_Int> &_ia_at, vector<_Int> &_ja_at, vector<_Flt> &_a_at);
      /// Combine sparsity of L and U
   static void CombineLUSp (int _n, 
                              vector<_Int> &_ia_al, vector<_Int> &_ja_al,
                              vector<_Int> &_ia_au, vector<_Int> &_ja_au,
                              vector<_Int> &_ia_alu, vector<_Int> &_ja_alu);
      /// Combine L and U data into extended pairs
   static void CombinePairs (int _n, 
                              vector<_Int> &_ia_al, vector<_Int> &_ja_al, vector<_Flt> &_a_al,
                              vector<_Int> &_ia_au, vector<_Int> &_ja_au, vector<_Flt> &_a_au,
                              vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu);
      /// Split pairs fct data into L and U parts with post filtering
   static void SplitPairsFilter (double _tau1, int _n, 
                                 vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                                 vector<_Int> &_ia_l, vector<_Int> &_ja_l, vector<_Flt> &_a_l,
                                 vector<_Int> &_ia_u, vector<_Int> &_ja_u, vector<_Flt> &_a_u);
      /// Combine L and U data into extended pairs
   static void CombineRowsLU (int _n, 
                              vector<_Int> &_ia_al, vector<_Int> &_ja_al, vector<_Flt> &_a_al,
                              vector<_Int> &_ia_au, vector<_Int> &_ja_au, vector<_Flt> &_a_au,
                              vector<_Int> &_ia_alu, vector<_Int> &_ida_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu);
      /// Perform ILU2 point factorization of the block with future diagonal modification (no structural control)
   static void Ilu2BlockIlu2 (double _pivmin, double _tau1, double _tau2, double _theta,
                              int _n, int _n_ini, 
                              vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                              vector<_Int> &_ia_lu, vector<_Int> &_ja_lu, vector<_Flt> &_a_lu,
                              double &_eigmin_att, double &_eigmax_att);
      /// Perform ILU2 point factorization of the block with future diagonal modification (with structural control)
   static void Ilu2BlockIlu2 (int _fcttype, double _pivmin, double _tau1, double _tau2, double _theta,
                              int _n, int _n_ini, 
                              vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<char> &_ja_char_alu, vector<_Flt> &_a_alu,
                              vector<_Int> &_ia_lu, vector<_Int> &_ja_lu, vector<char> &_ja_char_lu, vector<_Flt> &_a_lu,
                              double &_eigmin_att, double &_eigmax_att);
      /// Compute inverse scaling
   static void InverseDiag (int _n, vector<_Flt> &_sclU, vector<_Flt> &_invsclU);
      /// Rescale factor back
   static void RescaleU (int _n, vector<_Flt> &_sclU, vector<_Flt> &_invsclU,
                          vector<_Int> &_ia_u, vector<_Int> &_ja_u, vector<_Flt> &_a_u);
      /// Balance diagonal
   static void BalanceDiag (int _n, vector<_Int> &_ia_a, vector<_Int> &_ja_a, vector<_Flt> &_a_a,
                              vector<_Int> &_ia_l, vector<_Int> &_ja_l, vector<_Flt> &_a_l,
                              vector<_Int> &_ia_u, vector<_Int> &_ja_u, vector<_Flt> &_a_u,
                              double &_diacorr_min, double &_diacorr_max);
};

template <typename _Int, typename _Flt>
class CMatrix // Class that supports matrix data
{
   int n_list;               ///< The number of elements in the List and Ia arrays
   int n_list2;              ///< The number of elements in the List2 array
   int nz_ja;                ///< The number of elements in Ja array
   int nz_ja2;               ///< The number of elements in Ja2 array
   int nz_a;                 ///< The number of elements in A array
   vector<_Int> list_matr;   ///< List part of a matrix
   vector<_Int> list2_matr;  ///< List2 part of a matrix
   vector<_Int> ia_matr;     ///< Ia part of a matrix
   vector<_Int> ja_matr;     ///< Ja part of a matrix
   vector<_Int> ja2_matr;    ///< Ja2 part of a matrix
   vector<_Flt> a_matr;      ///< A  part of a matrix
// Get template names
   typedef _Int   TInt_type; ///< Reference to int type
   typedef _Flt TFloat_type; ///< Reference to float type
// Constructors and destructor
public:
      /// Epmty constructor
   CMatrix () { n_list = 0; list_matr.resize(1); list_matr[0] = 0; list2_matr.resize(1); list2_matr[0] = 0; 
                  ia_matr.resize(1); ia_matr[0] = 0; ja_matr.resize(1); ja_matr[0] = 0; ja2_matr.resize(1); ja2_matr[0] = 0; 
                  a_matr.resize(1); a_matr[0] = (_Flt)0.0e0; n_list2 = 0; nz_ja = 0; nz_ja2 = 0; nz_a = 0; }; 
      /// Copy constructor
   CMatrix (const CMatrix<_Int,_Flt> &_aa); 
      /// Equality operator
   CMatrix<_Int,_Flt> &operator= (const CMatrix<_Int,_Flt> &_aa); 
      /// Init by data
   CMatrix (int _n, _Int *_ia, _Int *_ja) {
      n_list = _n; n_list2 = 0; int nzja = (int)_ia[_n]; nz_ja = nzja; nz_ja2 = 0; nz_a = 0; 
      list_matr.resize(1); list_matr[0] = 0; list2_matr.resize(1); list2_matr[0] = 0; 
      ia_matr.resize(_n+1); ja_matr.resize(nzja+1); _Int *pia = &ia_matr[0]; _Int *pja = &ja_matr[0];
      int i; for (i=0;i<=_n;i++) pia[i] = _ia[i]; for (i=0;i<nzja;i++) pja[i] = _ja[i];
      ja2_matr.resize(1); ja2_matr[0] = 0; this->SetIdentityList();
   };
      /// Init by data
   CMatrix (int _nlist, _Int *_list, _Int *_ia, _Int *_ja) {
      n_list = _nlist; n_list2 = 0; int nzja = (int)_ia[_nlist]; nz_ja = nzja; nz_ja2 = 0; nz_a = 0; 
      list_matr.resize(_nlist+1); list2_matr.resize(1); list2_matr[0] = 0; 
      ia_matr.resize(_nlist+1); ja_matr.resize(nzja+1);
      _Int *plist = &list_matr[0]; _Int *pia = &ia_matr[0]; _Int *pja = &ja_matr[0];
      int i; for (i=0;i<_nlist;i++) plist[i] = _list[i]; for (i=0;i<=_nlist;i++) pia[i] = _ia[i]; for (i=0;i<nzja;i++) pja[i] = _ja[i];
      ja2_matr.resize(1); ja2_matr[0] = 0; a_matr.resize(1); a_matr[0] = (_Flt)0.0e0; 
   };
      /// Init by data
   CMatrix (int _nlist, _Int *_list, _Int *_ia, _Int *_ja, _Flt *_a) {
      n_list = _nlist; n_list2 = 0; int nzja = (int)_ia[_nlist]; nz_ja = nzja; nz_ja2 = 0; nz_a = 0; 
      list_matr.resize(_nlist+1); list2_matr.resize(1); list2_matr[0] = 0;
      ia_matr.resize(_nlist+1); ja_matr.resize(nzja+1); ja2_matr.resize(1); ja2_matr[0] = 0; a_matr.resize(nzja+1);
      _Int *plist = &list_matr[0]; _Int *pia = &ia_matr[0]; _Int *pja = &ja_matr[0]; _Flt *pa = &a_matr[0];
      int i; for (i=0;i<_nlist;i++) plist[i] = _list[i]; for (i=0;i<=_nlist;i++) pia[i] = _ia[i]; 
      for (i=0;i<nzja;i++) pja[i] = _ja[i]; for (i=0;i<nzja;i++) pa[i] = _a[i];
   };
      /// Init by data
   CMatrix (int _n, _Int *_ia, _Int *_ja, _Flt *_a) {
      n_list = _n; n_list2 = 0; int nzja = (int)_ia[_n]; nz_ja = nzja; nz_ja2 = 0; nz_a = nzja; 
      list_matr.resize(1); list_matr[0] = 0; list2_matr.resize(1); list2_matr[0] = 0; 
      ia_matr.resize(_n+1); ja_matr.resize(nzja+1); a_matr.resize(nzja+1);
      _Int *pia = &ia_matr[0]; _Int *pja = &ja_matr[0]; _Flt *pa = &a_matr[0];
      int i; for (i=0;i<=_n;i++) pia[i] = _ia[i]; for (i=0;i<nzja;i++) pja[i] = _ja[i]; for (i=0;i<nzja;i++) pa[i] = _a[i];
      this->SetIdentityList(); ja2_matr.resize(1); ja2_matr[0] = 0; 
   };
      /// Destructor
   ~CMatrix () {};
public:
// External functions
// Get/set functions
      /// Get N
   int GetN () {return n_list;}; 
   const int GetN () const {return n_list;}; 
      /// Get Nlist
   int GetNlist () {return n_list;}; 
   const int GetNlist () const {return n_list;}; 
      /// Get Nlist2
   int GetNlist2 () {return n_list2;}; 
   const int GetNlist2 () const {return n_list2;}; 
      /// Get Nzja
   int GetNzja () {return nz_ja;}; 
   const int GetNzja () const {return nz_ja;}; 
      /// Get Nzja2
   int GetNzja2 () {return nz_ja2;}; 
   const int GetNzja2 () const {return nz_ja2;}; 
      /// Get Nza
   int GetNza () {return nz_a;}; 
   const int GetNza () const {return nz_a;}; 
      /// Get List
   vector<_Int> * GetList () {return &list_matr;}; 
   const vector<_Int> * GetList () const {return &list_matr;}; 
      /// Get List array
   _Int * GetListArr () {return &list_matr[0];}; 
   const _Int * GetListArr () const {return &list_matr[0];}; 
      /// Get List2
   vector<_Int> * GetList2 () {return &list2_matr;}; 
   const vector<_Int> * GetList2 () const {return &list2_matr;}; 
      /// Get List2 array
   _Int * GetList2Arr () {return &list2_matr[0];}; 
   const _Int * GetList2Arr () const {return &list2_matr[0];}; 
      /// Get Ia
   vector<_Int> * GetIa () {return &ia_matr;}; 
   const vector<_Int> * GetIa () const {return &ia_matr;}; 
      /// Get Ia array
   _Int * GetIaArr () {return &ia_matr[0];}; 
   const _Int * GetIaArr () const {return &ia_matr[0];}; 
      /// Get Ja
   vector<_Int> * GetJa () {return &ja_matr;}; 
   const vector<_Int> * GetJa () const {return &ja_matr;}; 
      /// Get Ja array
   _Int * GetJaArr () {return &ja_matr[0];}; 
   const _Int * GetJaArr () const {return &ja_matr[0];}; 
      /// Get Ja2
   vector<_Int> * GetJa2 () {return &ja2_matr;}; 
   const vector<_Int> * GetJa2 () const {return &ja2_matr;}; 
      /// Get Ja2 array
   _Int * GetJa2Arr () {return &ja2_matr[0];}; 
   const _Int * GetJa2Arr () const {return &ja2_matr[0];}; 
      /// Get A
   vector<_Flt> * GetA () {return &a_matr;}; 
   const vector<_Flt> * GetA () const {return &a_matr;}; 
      /// Get A array
   _Flt * GetAArr () {return &a_matr[0];}; 
   const _Flt * GetAArr () const {return &a_matr[0];}; 
      /// Set N
   void SetN (int _n) {n_list = _n;}; 
      /// Set Nlist
   void SetNlist (int _n) {n_list = _n;}; 
      /// Set Nzja
   void SetNzja (int _nzja) {nz_ja = _nzja;}; 
      /// Set Nza
   void SetNza (int _nza) {nz_a = _nza;}; 
// Other
      /// Resize sp function
   void ResizeAndSetAllSp (int _nlist, int _nlist2, int _nzja, int _nzja2) {
      this->n_list = _nlist; this->n_list2 = _nlist2; this->nz_ja = _nzja; this->nz_ja2 = _nzja2; this->nz_a = 0; 
      this->list_matr.resize (_nlist+1); this->list2_matr.resize (_nlist2+1); this->ia_matr.resize (_nlist+1); this->ia_matr[0] = 0;
      this->ja_matr.resize (_nzja+1); this->ja2_matr.resize (_nzja2+1); this->a_matr.resize (1);
   };
      /// Resize function
   void ResizeAndSetAll (int _nlist, int _nlist2, int _nzja, int _nzja2, int _nza) {
      this->n_list = _nlist; this->n_list2 = _nlist2; this->nz_ja = _nzja; this->nz_ja2 = _nzja2; this->nz_a = _nza; 
      this->list_matr.resize (_nlist+1); this->list2_matr.resize (_nlist2+1); this->ia_matr.resize (_nlist+1); this->ia_matr[0] = 0;
      this->ja_matr.resize (_nzja+1); this->ja2_matr.resize (_nzja2+1); this->a_matr.resize (_nza+1);
   };
      /// Replace function
   void ReplaceFree (CMatrix<_Int,_Flt> &_aa) {
      this->n_list = _aa.n_list; this->n_list2 = _aa.n_list2; this->nz_ja = _aa.nz_ja; this->nz_ja2 = _aa.nz_ja2; this->nz_a = _aa.nz_a; 
      this->list_matr.swap (_aa.list_matr); this->list2_matr.swap (_aa.list2_matr); this->ia_matr.swap (_aa.ia_matr); 
      this->ja_matr.swap (_aa.ja_matr); this->ja2_matr.swap (_aa.ja2_matr); this->a_matr.swap (_aa.a_matr); _aa.Clean();
   };
      /// Clean matrix data
   void Clean () {
      vector<_Int> list_dummy(1); vector<_Int> list2_dummy(1); vector<_Int> ia_dummy(1); 
      vector<_Int> ja_dummy(1); vector<_Int> ja2_dummy(1); vector<_Flt> a_dummy(1);
      list_matr.swap (list_dummy); list2_matr.swap (list2_dummy); ia_matr.swap (ia_dummy); 
      ja_matr.swap (ja_dummy); ja2_matr.swap (ja2_dummy); a_matr.swap (a_dummy); 
      n_list = 0; n_list2 = 0; nz_ja = 0; nz_ja2 = 0; nz_a = 0; ia_matr[0] = 0;
   };
      /// Set identity list
   void SetIdentityList () {list_matr.resize(n_list+1); int i;  for (i=0;i<n_list;i++) list_matr[i] = (_Int)i;}; 
      /// Get sparsity
   void GetSparsity (const CMatrix<_Int,_Flt> &_a_sp);
      /// Compute transposed sparsity for incomplete list
   void TransposedSparsityListSp (int &_icycle, int *_imask, int *_indarr, int *_iptr, int *_listloc, int *_ialoc,
                                    CMatrix<_Int,_Flt> &_at) const;
      /// Compute transposed matrix for incomplete list
   void TransposedSparsityList (int &_icycle, int *_imask, int *_indarr, int *_iptr, int *_listloc, int *_ialoc,
                                 CMatrix<_Int,_Flt> &_at) const;
      /// Extend sparsity by zeroes
   void ExtendSparsity (const CMatrix<_Int,_Flt> &_a_sp,
                         CMatrix<_Int,_Flt> &_a_ext) const;
      /// Perform multiplication of matrices stored by rows, resulting matrix is also stored by rows
   static void MatrixByMatrixMultiply (int &_icycle, int *_imask, int *_imask1, int *_indarr, int *_listloc, _Flt *_fmask,
                                          CMatrix<_Int,_Flt> &_a, CMatrix<_Int,_Flt> &_b,
                                          CMatrix<_Int,_Flt> &_a_times_b);
      /// Compute the packed size
   int GetPackedSize () const;
      /// Fill char array by the packed data
   void FillPacked (int _length, char *_obj) const;
      /// Unpack data
   void UnPack (int _length, char *_obj);
};

template <typename _Int, typename _Flt, typename _FltVect>
class CMvmSlv_impl // Class that supports multiplication computations
{
public:
      /// Multiply by super sparse matrix by rows and add result into prescribed positions
   static void MvmA (int _nlist, const vector<_Int> &_list_alu, 
                        const vector<_Int> &_ia_alu, const vector<_Int> &_ja_alu, const vector<_Flt> &_a_alu,
                        const _FltVect *_x, _FltVect *_ax);
      /// Multiply by super sparse matrix by columns add result into prescribed positions
   static void MvmAT (int _nlist, const vector<_Int> &_list_alu, 
                        const vector<_Int> &_ia_alu, const vector<_Int> &_ja_alu, const vector<_Flt> &_a_alu,
                        const _FltVect *_x, _FltVect *_ax);
      /// Solve with L, L is stored by columns (diag is inverted)
   static void SolveL (int _n, const vector<_Int> &_ia_l, const vector<_Int> &_ja_l, const vector<_Flt> &_a_l,
                        const _FltVect *_x, _FltVect *_lx);
      /// Solve with U, U is stored by rows (diag is inverted)
   static void SolveU (int _n, const vector<_Int> &_ia_u, const vector<_Int> &_ja_u, const vector<_Flt> &_a_u,
                        const _FltVect *_x, _FltVect *_ux);
};

template <typename _Int, typename _Flt, typename _FltVect>
class CMvmSlv // Class that supports multiplication computations
{
public:
      /// Multiply by super sparse matrix by rows and add result into prescribed positions
   static void MvmA (const CMatrix<_Int,_Flt> &_a_matr,
                     const _FltVect *_x, _FltVect *_ax) {
      CMvmSlv_impl<_Int, _Flt, _FltVect>::MvmA (_a_matr.GetNlist(), *(_a_matr.GetList()), *(_a_matr.GetIa()), *(_a_matr.GetJa()), *(_a_matr.GetA()), 
                                                            _x, _ax);
   };
      /// Multiply by super sparse matrix by rows and add result into prescribed positions
   static void MvmAT (const CMatrix<_Int,_Flt> &_a_matr,
                     const _FltVect *_x, _FltVect *_ax) {
      CMvmSlv_impl<_Int, _Flt, _FltVect>::MvmAT (_a_matr.GetNlist(), *(_a_matr.GetList()), *(_a_matr.GetIa()), *(_a_matr.GetJa()), *(_a_matr.GetA()), 
                                                            _x, _ax);
   };
      /// Solve with L, L is stored by columns (diag is inverted)
   static void SolveL (const CMatrix<_Int,_Flt> &_l_matr,
                        const _FltVect *_x, _FltVect *_lx) {
      CMvmSlv_impl<_Int, _Flt, _FltVect>::SolveL (_l_matr.GetN(), *(_l_matr.GetIa()), *(_l_matr.GetJa()), *(_l_matr.GetA()), 
                                                            _x, _lx);
   };
      /// Solve with U, U is stored by rows (diag is inverted)
   static void SolveU (const CMatrix<_Int,_Flt> &_u_matr,
                        const _FltVect *_x, _FltVect *_ux) {
      CMvmSlv_impl<_Int, _Flt, _FltVect>::SolveU (_u_matr.GetN(), *(_u_matr.GetIa()), *(_u_matr.GetJa()), *(_u_matr.GetA()), 
                                                            _x, _ux);
   };
};

template <typename _Int, typename _Flt>
class CIlu2 // Class that supports factorization computations
{
public:
// Functions that support matrix structure as a whole
      /// Compute optimal ordering for the matrix
   static void ComputeOptimalOrder (const CMatrix<_Int,_Flt> &_a_matr, int _ordtype, vector<int> &_order) {
      CIlu2_impl<_Int,_Flt>::ComputeOptimalOrder (_ordtype, _a_matr.GetNlist(), *(_a_matr.GetIa()), *(_a_matr.GetJa()), _order);
   };
      /// Compute ordered matrix
   static void ReorderMatrix (CMatrix<_Int,_Flt> &_a_matr, const vector<int> &_order, CMatrix<_Int,_Flt> &_a_ord) {
      _a_ord.SetN (_a_matr.GetNlist()); _a_ord.SetNzja (_a_matr.GetNzja()); _a_ord.SetNza (_a_matr.GetNza());
      CIlu2_impl<_Int,_Flt>::ReorderMatrix (_a_matr.GetNlist(), _order, *(_a_matr.GetIa()), *(_a_matr.GetJa()), *(_a_matr.GetA()),
                                                *(_a_ord.GetIa()), *(_a_ord.GetJa()), *(_a_ord.GetA()));
      _a_ord.SetIdentityList ();
   };
      /// Perform ILU2 point factorization of the block with future diagonal modification
   static void Ilu2Matrix (CMatrix<_Int,_Flt> &_a_matr, SSolverParams &_params,
                           CMatrix<_Int,_Flt> &_l_matr, CMatrix<_Int,_Flt> &_u_matr,
                           double &_eigmin_att, double &_eigmax_att) {
      _l_matr.SetN (_a_matr.GetNlist()); _u_matr.SetN (_a_matr.GetNlist());
      CIlu2_impl<_Int,_Flt>::Ilu2Block (_params.sctype, _params.nitersc, _params.fcttype, _params.pivmin, _params.tau1, _params.tau2, _params.theta,
                                             _a_matr.GetNlist(), *(_a_matr.GetIa()), *(_a_matr.GetJa()), *(_a_matr.GetA()),
                                             *(_l_matr.GetIa()), *(_l_matr.GetJa()), *(_l_matr.GetA()), *(_u_matr.GetIa()), *(_u_matr.GetJa()), *(_u_matr.GetA()),
                                             _eigmin_att, _eigmax_att);
      int nlist_a = _a_matr.GetNlist(); _Int *pia_lmatr = _l_matr.GetIaArr(); _Int *pia_umatr = _u_matr.GetIaArr(); 
      int nzja_l = (int)pia_lmatr[nlist_a]; int nzja_u = (int)pia_umatr[nlist_a];
      _l_matr.SetNzja (nzja_l); _l_matr.SetNza (nzja_l); _u_matr.SetNzja (nzja_u); _u_matr.SetNza (nzja_u); 
      _l_matr.SetIdentityList (); _u_matr.SetIdentityList ();
   };
      /// Balance diagonal
   static void BalanceDiag (CMatrix<_Int,_Flt> &_a_matr, CMatrix<_Int,_Flt> &_l_matr, CMatrix<_Int,_Flt> &_u_matr,
                     double &_diacorr_min, double &_diacorr_max) {
      CIlu2_impl<_Int,_Flt>::BalanceDiag (
                                                _a_matr.GetNlist(), *(_a_matr.GetIa()), *(_a_matr.GetJa()), *(_a_matr.GetA()),
                                                *(_l_matr.GetIa()), *(_l_matr.GetJa()), *(_l_matr.GetA()),
                                                *(_u_matr.GetIa()), *(_u_matr.GetJa()), *(_u_matr.GetA()),
                                                _diacorr_min, _diacorr_max);
   };
      /// Print matrix data
   static void PrintMatrix (ofstream &_fout, CMatrix<_Int,_Flt> &_a_matr);
};

template <typename _Int, typename _Flt, typename _Flt2>
class CMatrixConv // Class that supports conversion of data
{
public:
      /// Convert matrix data
   static void InitAndConv (int _n, _Int *_ia, _Int *_ja, _Flt *_a, CMatrix<_Int,_Flt2> &_amatr);
};

template <typename _Int, typename _Flt>
class CHMatrix /// Class that supports parallel submatrix data
{
   int nzblk;                            ///< The number of blocks
   CMatrix<int,float> hmatr_str;         ///< Block sparsity of a submatrix
   vector<CMatrix<_Int,_Flt> > asub_arr; ///< Data of each block
// Get template names
   typedef _Int THInt_type;   ///< Reference to int type
   typedef _Flt THFloat_type; ///< Reference to float type
// Constructors and destructor
public:
      /// Epmty constructor
   CHMatrix () {nzblk=0; asub_arr.resize (1); }; 
      /// Copy constructor
   CHMatrix (const CHMatrix<_Int,_Flt> &_aa); 
      /// Оператор копирования
   CHMatrix<_Int,_Flt> &operator= (const CHMatrix<_Int,_Flt> &_aa); 
      /// Init by data
   CHMatrix (int _iblk, int _nlist, _Int *_list, _Int *_ia, _Int *_ja, _Flt *_a,
               int _nblks, long long *_blks, 
               int &_icycle, int *_imaskblk);
      /// Destructor
   ~CHMatrix () {};
// External functions
// Get/set functions
      /// Get nzblk
   int GetNzblk () {return nzblk;}; 
   const int GetNzblk () const {return nzblk;}; 
      /// Get hmatr_str
   CMatrix<int,float> * GetHMatrStr () { return &hmatr_str; };
   const CMatrix<int,float> * GetHMatrStr () const { return &hmatr_str; };
      /// Get asub_arr
   vector<CMatrix<_Int,_Flt> > * GetASub () {return &asub_arr;}; 
   const vector<CMatrix<_Int,_Flt> > * GetASub () const {return &asub_arr;}; 
      /// Get asub_arr
   CMatrix<_Int,_Flt> * GetASubArr () {return &asub_arr[0];}; 
   const CMatrix<_Int,_Flt> * GetASubArr () const {return &asub_arr[0];}; 
      /// Set nzblk
   void SetNzblk (int _nzblk) {nzblk = _nzblk;}; 
      /// Resize ASub
   void ResizeASub (int _nzblk) {asub_arr.resize(_nzblk+1);}; 
      /// Replace function
   void ReplaceFree (CHMatrix<_Int,_Flt> &_aa) {
      this->nzblk = _aa.nzblk; this->hmatr_str.ReplaceFree (_aa.hmatr_str); this->asub_arr.swap (_aa.asub_arr); _aa.Clean();
   };
      /// Clean hmatrix data
   void Clean () {
      CMatrix<int,float> HMatr_dummy; vector<CMatrix<_Int,_Flt> > ASub_dummy (1);
      hmatr_str.ReplaceFree(HMatr_dummy); asub_arr.swap (ASub_dummy); nzblk = 0;
   };
// Other functions
      /// Split structural matrix into the set of matrices
   static void SplitMatrSpIntoHMatrSp (int _nblks, long long *_blks, const CMatrix<_Int,_Flt> &_amatr, CHMatrix<_Int,_Flt> &_ahmatr); 
      /// Compute the symmetrized submatrices
   static void SymmetrizeSubmatrices (void *_comm,
                                       int _nblks, long long *_blks, int *_blk2cpu, 
                                       CHMatrix<_Int,_Flt> *_hmatr_arr, CHMatrix<_Int,_Flt> *_hmatr_symm_arr);
      /// Compute the extended lists
   static void ExtendedLists (void *_comm, int _ncycle,
                              int _nblks, long long *_blks, int *_blk2cpu, 
                              CHMatrix<_Int,_Flt> *_hmatr_arr,
                              int *_nlist_ext_arr, vector<int> *_list_ext_arr);
      /// Get the extended submatrices
   static void GetExtendedSubmatrices (void *_comm,
                                       int _nblks, long long *_blks, int *_blk2cpu, 
                                       CHMatrix<_Int,_Flt> *_hmatr_arr,
                                       int *_nlist_ext_arr, vector<int> *_list_ext_arr,
                                       CMatrix<_Int,_Flt> *_matr_ext_arr);
      /// Compute the packed size
   long long GetPackedSize () const;
      /// Fill char array by the packed data
   void FillPacked (long long _length, char *_obj) const;
      /// Unpack data
   void UnPack (long long _length, char *_obj);
// Low level support functions
      /// Filter the extended lists for backward only extention
   static void FilterListsBack (int _myid, int _nblks, int *_blk2cpu, 
                                  int *_nlist_ext_arr, vector<int> *_list_ext_arr);
      /// Compute Ja2 array by binary search
   static void ComputeJa2 (int _nblks, long long *_blks,
                           int _nlist, const _Int *_ia, const _Int *_ja, 
                           _Int *_ja2);
      /// Print sparsity with boxes that show nonzero blocks
   static void Str2PsBox (const CMatrix<_Int,_Flt> &_amatr, char *_fname, int _nblks, int *_blks);
      /// Print sparsity with boxes that show nonzero blocks
   static void Str2PsBox (int _collap, const CMatrix<_Int,_Flt> &_amatr, char *_fname, int _nblks, int *_blks);
      /// Print hmatrix data
   void PrintHMatrix (ofstream &_fout);
};

template <typename _Int, typename _Flt, typename _FltVect>
class CMvmPar // Class that supports parallel multiplication computations
{
   void *pcomm;                     ///< Pointer to MPI communicator
   int nblks;                       ///< Number of blocks
   long long *pblks;                ///< Pointer to blocks partitioning
   int *pblk2cpu;                   ///< Pointer to block to cpu distribution
   int ni_cpu;                      ///< The size of the local part of vector data
   vector<int> ibsblk;              ///< Base adresses of local parts of the vector data
   int nlistblk_own;                ///< The number of local blocks
   vector<int> listblk_own;         ///< The list of own blocks
   CHMatrix<_Int,_Flt> *phmatr;     ///< Pointer to hyperblocks array
   int nsends;                      ///< The number of sends for computations
   vector<int> snd2cpu;             ///< The list of cpus to be sended to for computations
   vector<int> ia_sends;            ///< The sizes of sends data for computations
   vector<int> ind_sends;           ///< The indices of local data for sends
   vector<_FltVect> x_send;         ///< The work memory for sends
   int nrecvs;                      ///< The number of recvs for computations
   vector<int> rcv2cpu;             ///< The list of cpus to be received data from for computations
   vector<int> ia_recvs;            ///< The sizes of recvs for computations
   vector<int> ja_recvs;            ///< The pairs of local data that describe recvs
   vector<_FltVect> x_recv;         ///< The work memory for recvs
   int nblks_recvs;                 ///< The number of received block parts
   vector<int> listblk_recvs;       ///< The list of block numbers of received blocks
   vector<int> ialist_recvs;        ///< The numbers of block rows contatining the block
   vector<int> japairs_recvs;       ///< The complete list and index numbers inside of block rows containing each of the received blocks
   vector<int> iablk_recvs;         ///< The sizes of copy data for each received block
   vector<int> ind_recvs;           ///< The list of indices inside block for received blocks
   vector<_FltVect> x_temp;         ///< The work memory for external blocks
// Constructors and destructor
private:
      /// Copy constructor
   CMvmPar (const CMvmPar<_Int,_Flt,_FltVect> &_aa) {}; 
      /// Copy operator
   CMvmPar<_Int,_Flt,_FltVect> &operator= (const CMvmPar<_Int,_Flt,_FltVect> &_aa) { return *this; }; 
public:
      /// Epmty constructor
   CMvmPar () {pcomm = NULL; nblks = 0; pblks = NULL; pblk2cpu = NULL; phmatr = NULL; 
                     ni_cpu = 0; nlistblk_own = 0; nsends = 0; nrecvs = 0; nblks_recvs = 0;}; 
      /// Destructor
   ~CMvmPar () {};
// External functions
      /// Init control data
   void InitControl (void *_pcomm, int _nblks, long long *_blks, int *_blk2cpu) {
      pcomm = _pcomm; int myid = CExchange::GetMyid (pcomm); nblks = _nblks; pblks = _blks; pblk2cpu = _blk2cpu; 
      ibsblk.resize (nblks+1); int *pibsblk = &ibsblk[0]; ni_cpu = 0; nlistblk_own = 0; int i; 
      for (i=0;i<nblks;i++) { if (pblk2cpu[i] == myid) { nlistblk_own++; pibsblk[i] = ni_cpu; ni_cpu += (int)(pblks[i+1]-pblks[i]); } else { pibsblk[i] = -1; }; };
      listblk_own.resize (nlistblk_own+1); int *plistblk_own = &listblk_own[0]; nlistblk_own = 0; 
      for (i=0;i<nblks;i++) { if (pblk2cpu[i] == myid) { plistblk_own[nlistblk_own] = i; nlistblk_own++; } };
   };
      /// Init MvmA data
   void InitMvmA (CHMatrix<_Int,_Flt> *_hmatr_arr);
      /// Perform MvmA computations
   void MvmA (const _FltVect *_x, _FltVect *_ax);
      /// Clean MvmA structure
   void Clean ();
};

template <typename _Int, typename _Flt, typename _FltVect>
class CSlvPar // Class that supports parallel solve computations (including L, U and ordering)
{
   void *pcomm;                     ///< Pointer to MPI communicator
   int nblks;                       ///< Number of blocks
   long long *pblks;                ///< Pointer to blocks partitioning
   long long *pblks_ext;            ///< Pointer to the extended blocks partitioning
   int *pblk2cpu;                   ///< Pointer to block to cpu distribution
   int ni_cpu;                      ///< The size of the local part of vector data
   vector<int> ibsblk;              ///< Base adresses of local parts of the vector data
   int nlistblk_own;                ///< The number of local blocks
   vector<int> listblk_own;         ///< The list of own blocks
   vector<int> *plistpairs_ext;     ///< Pointer to the extended data lists
   CMatrix<_Int,_Flt> *pmatrL;      ///< Pointer to L blocks array
   CMatrix<_Int,_Flt> *pmatrU;      ///< Pointer to U blocks array
   vector<int> *porderLU;           ///< Pointer to the extended ordering data if any
   int nsends;                      ///< The number of sends for computations
   vector<int> snd2cpu;             ///< The list of cpus to be sended to for computations
   vector<int> ia_sends;            ///< The sizes of sends data for computations
   vector<int> ind_sends;           ///< The indices of local data for sends
   vector<_FltVect> x_send;         ///< The work memory for sends
   int nrecvs;                      ///< The number of recvs for computations
   vector<int> rcv2cpu;             ///< The list of cpus to be received data from for computations
   vector<int> ia_recvs;            ///< The sizes of recvs for computations
   vector<_FltVect> x_recv;         ///< The work memory for recvs
   vector<int> ialist_recvs;        ///< For each local block the number of updates from received data
   vector<int> jalist_recvs;        ///< The lists of index numbers for received data for updates in current local block
   vector<_FltVect> x_temp;         ///< The work memory for extended blocks
// Constructors and destructor
private:
      /// Copy constructor
   CSlvPar (const CSlvPar<_Int,_Flt,_FltVect> &_aa) {}; 
      /// Copy operator
   CSlvPar<_Int,_Flt,_FltVect> &operator= (const CSlvPar<_Int,_Flt,_FltVect> &_aa) { return *this; }; 
public:
      /// Epmty constructor
   CSlvPar () {pcomm = NULL; nblks = 0; pblks = NULL; pblk2cpu = NULL; plistpairs_ext = NULL; pmatrL = NULL; pmatrU = NULL; porderLU = NULL;
                     ni_cpu = 0; nlistblk_own = 0; nsends = 0; nrecvs = 0; }; 
      /// Destructor
   ~CSlvPar () {};
// External functions
      /// Init control data
   void InitControl (void *_pcomm, int _nblks, long long *_blks, long long *_blks_ext, int *_blk2cpu) {
      pcomm = _pcomm; int myid = CExchange::GetMyid (pcomm); nblks = _nblks; pblks = _blks; pblks_ext = _blks_ext; pblk2cpu = _blk2cpu; 
      ibsblk.resize (nblks+1); int *pibsblk = &ibsblk[0]; ni_cpu = 0; nlistblk_own = 0; int i; 
      for (i=0;i<nblks;i++) { if (pblk2cpu[i] == myid) { nlistblk_own++; pibsblk[i] = ni_cpu; ni_cpu += (int)(pblks[i+1]-pblks[i]); } else { pibsblk[i] = -1; }; };
      listblk_own.resize (nlistblk_own+1); int *plistblk_own = &listblk_own[0]; nlistblk_own = 0; 
      for (i=0;i<nblks;i++) { if (pblk2cpu[i] == myid) { plistblk_own[nlistblk_own] = i; nlistblk_own++; } };
   };
      /// Init SolveLU data
   void InitSolveLU (vector<int> *_listpairs_ext, CMatrix<_Int,_Flt> *_matrL, CMatrix<_Int,_Flt> *_matrU, vector<int> *_orderLU);
      /// Perform SolveLU computations
   void SolveLU (const _FltVect *_x, _FltVect *_px);
      /// Clean Slv structure
   void Clean ();
};

template <typename _Int, typename _Flt, typename _FltVect>
class CSolver // Class that supports parallel solver computations
{
   void *pcomm;                                   ///< Pointer to MPI communicator
   SSolverParams params;                          ///< Solver params
   int nblks;                                     ///< Number of blocks
   vector<long long> blks;                        ///< Pointer to blocks partitioning
   vector<long long> blks_ext;                    ///< Pointer to the extended blocks partitioning
   vector<int> blk2cpu;                           ///< Pointer to block to cpu distribution
   vector<CHMatrix<_Int,_Flt> > hmatr_arr;        ///< Matrix data
   vector<int> nlist_ext_arr;                     ///< The number of elems in extension
   vector<vector<int> > list_ext_arr;             ///< The pairs of extension indices
   vector<vector<int> > order_LU;                 ///< Ordering data
   vector<CMatrix<_Int,float>  > matrL_float;     ///< Float  L part
   vector<CMatrix<_Int,double> > matrL_double;    ///< Double L part
   vector<CMatrix<_Int,float>  > matrU_float;     ///< Float  U part
   vector<CMatrix<_Int,double> > matrU_double;    ///< Double U part
   CMvmPar<_Int,_Flt,_FltVect> mvm;               ///< Support MvmA
   CSlvPar<_Int,float,  _FltVect> slv_float;      ///< Support float  Slv
   CSlvPar<_Int,double, _FltVect> slv_double;     ///< Support double Slv
// Constructors and destructor
private:
      /// Copy constructor
   CSolver (const CSolver<_Int,_Flt,_FltVect> &_aa) {}; 
      /// Copy operator
   CSolver<_Int,_Flt,_FltVect> &operator= (const CSolver<_Int,_Flt,_FltVect> &_aa) { return *this; }; 
public:
      /// Epmty constructor
   CSolver () {}; 
      /// Destructor
   ~CSolver () {};
// Get/set functions
      /// Get comm
   void * GetComm () { return pcomm; };
      /// Get nblks
   int GetNblks () {return nblks;}; 
      /// Get blks
   long long * GetBlks () { return &blks[0]; };
      /// Get blk2cpu
   int * GetBlk2cpu () { return &blk2cpu[0]; };
// External functions
      /// Prepare solver structures including performing parallel fct
   void PrepareSolver (void *_pcomm, SSolverParams *_params,
                        int _nblks, long long *_blks, int *_blk2cpu,
                        bool _b_store_matrix, _Int *_ia, _Int *_ja, _Flt *_a);
      /// Perform MvmA parallel computations
   void MvmA (const _FltVect *_x, _FltVect *_ax) {mvm.MvmA(_x,_ax);};
      /// Perform SlvLU parallel computations
   void SlvLU (const _FltVect *_x, _FltVect *_px) {if (params.prec_float == 1) { slv_float.SolveLU(_x,_px); } else { slv_double.SolveLU(_x,_px); }; };
      /// Clean MvmA structures
   void CleanMvmA () {
      mvm.Clean(); vector<CHMatrix<_Int,_Flt> > hmatr_arr_dummy; hmatr_arr.swap (hmatr_arr_dummy);
   };
      /// Clean Slv structures
   void CleanSlv () {
      vector<int> nlist_ext_arr_dummy; vector<vector<int> > list_ext_arr_dummy; vector<vector<int> > order_LU_dummy;
      vector<CMatrix<_Int,float>  > matrL_float_dummy; vector<CMatrix<_Int,double> > matrL_double_dummy;
      vector<CMatrix<_Int,float>  > matrU_float_dummy; vector<CMatrix<_Int,double> > matrU_double_dummy;
      nlist_ext_arr.swap (nlist_ext_arr_dummy); list_ext_arr.swap (list_ext_arr_dummy); order_LU.swap (order_LU_dummy);
      matrL_float.swap (matrL_float_dummy); matrL_double.swap (matrL_double_dummy); matrU_float.swap (matrU_float_dummy); matrU_double.swap (matrU_double_dummy); 
      slv_float.Clean(); slv_double.Clean();
   };
      /// Clean all structures
   void Clean () {
      pcomm = NULL; nblks = 0; vector<long long> blks_dummy; vector<long long> blks_ext_dummy; vector<int> blk2cpu_dummy;
      blks.swap (blks_dummy); blks_ext.swap (blks_ext_dummy); blk2cpu.swap (blk2cpu_dummy); CleanMvmA(); CleanSlv(); params.SetDefaults(); 
   };
      /// Perform iterations of the iterative scheme (BiCGStab)
   void BiCGStab (int _niter_max, double _eps, int _ichk, int _msglev, ofstream *_fout, 
                  _FltVect *_rhs, _FltVect *_sol,
                  double &_rhs_norm, double &_res_ini, int &_niter, double &_res_fin);
};

} // namespace k3d

#endif

