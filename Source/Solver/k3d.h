//------------------------------------------------------------------------------------------------
// File: k3d.h
//------------------------------------------------------------------------------------------------

#ifndef __K3D_H
#define __K3D_H

#include <iostream>
#include <fstream>
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
      /// Synchronous send data to the other CPU
   static void Send (void *_comm, int _rank, int _msgtag, int _length, char *_buffer);
      /// Asynchronous send data to the other CPU
   static void ISend (void *_comm, int _rank, int _msgtag, int _length, char *_buffer, int _indrecv, void *_recv_arr);
      /// Synchronous receive data from the other CPU
   static void Recv (void *_comm, int _rank, int _msgtag, int _length, char *_buffer);
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
// Get MPI wall time
   static double GetWallTimeMPI(); 
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
   static void SetByOnes (int _n, _FltVect *_x); ///< Set vector data by ones
   static _FltVect ScProd (int _n, const _FltVect *_x, const _FltVect *_y); ///< Compute scalar product
   static void CopyVector (int _n, const _FltVect *_x, _FltVect *_y); ///< Copy array
   static void AddReplaceVector (int _n, const _FltVect *_x1, _FltVect *_x1plus2); ///< Add vector data
   static void SubtractReplaceVector (int _n, const _FltVect *_x1, _FltVect *_x2minus1); ///< Subtract vector data
   static void UpdateVector (int _n, const _FltVect *_value, const _FltVect *_arr_x, _FltVect *_arr_y); ///< Update array (axpy): y=a*x+y
   static void UpdateVectorMinus (int _n, const _FltVect *_value, const _FltVect *_arr_x, _FltVect *_arr_y); ///< Update array for minus alpha (axpy)
   static void UpdateVectorReversed (int _n, const _FltVect *_value, const _FltVect *_arr_x, _FltVect *_arr_y); ///< Update array reversed (aypx)
   static void OrderVector (int _n, const vector<int> &_order, 
                              const _FltVect *_x, _FltVect *_x_ord); ///< Order vector data
   static void InvOrderVector (int _n, const vector<int> &_order, 
                              const _FltVect *_x, _FltVect *_x_ord); ///< Inverse order vector data
      /// Compute Hausholder transformation
   static void Housholder (int _m, _FltVect &_alpha, _FltVect *_x, _FltVect &_tau); 
      /// Compute QR decomposition
   static void QrdBlock (int _ncol, int _nrow,
                           _FltVect *_qblk, int _ldq, _FltVect *_tau);
      /// Multiply by Q
   static void MvmQ_Housholder (int _nrhs, int _nrows, int _ncols,
                                 _FltVect *_qblk, int _ldq, _FltVect *_tau, 
                                 _FltVect *_qx, int _ldqx);
      /// Multiply by QH
   static void MvmQH_Housholder (int _nrhs, int _nrows, int _ncols,
                                 _FltVect *_qblk, int _ldq, _FltVect *_tau,
                                 _FltVect *_qx, int _ldqx);
      /// Solve triangular system
   static void SolveR (char _slvtype, int _n, _FltVect *_rmatr, int _ldr, _FltVect *_rhs_sol);
      /// Compute polinomial via square upper Hessenberg matrix
   static void Polynomial (int _n, int _ncoef, _FltVect *_hmatrix, int _ldh, vector<double> &_coef);
      /// Compute Cholessky decomposition for matrices stored by columns
   static void CholesskyColumns (double &_diamod, int _n, _FltVect *_amatr, int _lda,
                                 _FltVect *_uarr, int _ldu,
                                 double *_dia_arr, double &_eigmin_att, double &_eigmax_att);
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
                           double &_sclmin_att, double &_sclmax_att, 
                           int &_nmodif, double &_eigmin_att, double &_eigmax_att);
      /// Compute scaling
   static void ComputeScaling (int _sctype, int _nitersc,
                                 int _n, 
                                 vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                                 vector<_Flt> &_sclL, vector<_Flt> &_sclU, 
                                 double &_sclmin_att, double &_sclmax_att);
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
                              int &_nmodif, double &_eigmin_att, double &_eigmax_att);
      /// Perform ILU2 point factorization of the block with future diagonal modification (with structural control)
   static void Ilu2BlockIlu2 (int _fcttype, double _pivmin, double _tau1, double _tau2, double _theta,
                              int _n, int _n_ini, 
                              vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<char> &_ja_char_alu, vector<_Flt> &_a_alu,
                              vector<_Int> &_ia_lu, vector<_Int> &_ja_lu, vector<char> &_ja_char_lu, vector<_Flt> &_a_lu,
                              int &_nmodif, double &_eigmin_att, double &_eigmax_att);
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
   int GetN () const {return n_list;}; 
      /// Get Nlist
   int GetNlist () {return n_list;}; 
   int GetNlist () const {return n_list;}; 
      /// Get Nlist2
   int GetNlist2 () {return n_list2;}; 
   int GetNlist2 () const {return n_list2;}; 
      /// Get Nzja
   int GetNzja () {return nz_ja;}; 
   int GetNzja () const {return nz_ja;}; 
      /// Get Nzja2
   int GetNzja2 () {return nz_ja2;}; 
   int GetNzja2 () const {return nz_ja2;}; 
      /// Get Nza
   int GetNza () {return nz_a;}; 
   int GetNza () const {return nz_a;}; 
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
                           double &_sclmin_att, double &_sclmax_att, 
                           int &_nmodif, double &_eigmin_att, double &_eigmax_att) {
      _l_matr.SetN (_a_matr.GetNlist()); _u_matr.SetN (_a_matr.GetNlist());
      CIlu2_impl<_Int,_Flt>::Ilu2Block (_params.sctype, _params.nitersc, _params.fcttype, _params.pivmin, _params.tau1, _params.tau2, _params.theta,
                                             _a_matr.GetNlist(), *(_a_matr.GetIa()), *(_a_matr.GetJa()), *(_a_matr.GetA()),
                                             *(_l_matr.GetIa()), *(_l_matr.GetJa()), *(_l_matr.GetA()), *(_u_matr.GetIa()), *(_u_matr.GetJa()), *(_u_matr.GetA()),
                                             _sclmin_att, _sclmax_att, _nmodif, _eigmin_att, _eigmax_att);
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
      /// Equality operator
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
   int GetNzblk () const {return nzblk;}; 
      /// Get hmatr_str
   CMatrix<int,float> * GetHMatrStr () { return &hmatr_str; };
   const CMatrix<int,float> * GetHMatrStr () const { return &hmatr_str; };
      /// Get asub_arr
   vector<CMatrix<_Int,_Flt> > * GetASub () {return &asub_arr;}; 
   const vector<CMatrix<_Int,_Flt> > * GetASub () const {return &asub_arr;}; 
      /// Get asub_arr
   CMatrix<_Int,_Flt> * GetASubArr () {return &asub_arr[0];}; 
   const CMatrix<_Int,_Flt> * GetASubArr () const {return &asub_arr[0];}; 
      /// Get nzatot
   int GetNzatot () {
      int nzatot = 0; CMatrix<_Int,_Flt> *pasub = &asub_arr[0]; int i; for (i=0;i<nzblk;i++) nzatot += pasub[i].GetNza();
      return nzatot;
   }; 
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

/// Class that supports tree data type
class CTree
{
   int root_id;                          ///< tree root
   int nlev;                             ///< number of levels
   int nnodes;                           ///< number of nodes
   vector<int> father;                   ///< father node for each node
   vector<int> nchilds;                  ///< number of childs for each node
   vector<vector<int> > childs_list;     ///< lists of childs for nodes
   vector<int> nodes_lev_id;             ///< level number for each node
   vector<int> nnodes_lev;               ///< number of nodes on each level
   vector<vector<int> > nodes_lev_list;  ///< list of nodes on each level
   vector<int> subtree_beg;              ///< subtree start for each node
public:
// Functions
// Constructors and destructor
      /// Default constructor
   CTree () {
      root_id = 0; nlev = 0; nnodes = 0; father.resize (1); nchilds.resize (1); childs_list.resize (1); 
      nodes_lev_id.resize (1); nnodes_lev.resize (1); nodes_lev_list.resize (1); subtree_beg.resize (1); 
   }; 
      /// Main constructor
   CTree (int _nnodes_ini, int _nchilds); 
      /// Destructor
   ~CTree () {}; 
// Operator functions
   CTree &operator= (CTree &_tree); ///< Copy operator
// Get/set interface functions
      /// Get root_id
   int GetRootId () const {return root_id;}; 
      /// Get nlev
   int GetNlev () const {return nlev;}; 
      /// Get nnodes
   int GetNnodes () const {return nnodes;}; 
      /// Get father
   int * GetFather () {return &father[0];}; 
      /// Get nchilds
   int * GetNchilds () {return &nchilds[0];}; 
      /// Get childs_list
   vector<int> * GetChildsList () {return &childs_list[0];}; 
      /// Get nodes_lev_id
   int * GetNodesLevId () {return &nodes_lev_id[0];}; 
      /// Get nnodes_lev
   int * GetNNodesLev () {return &nnodes_lev[0];}; 
      /// Get nodes_lev_list
   vector<int> * GetNodesLevList () {return &nodes_lev_list[0];}; 
      /// Get subtree_beg
   int * GetSubtreeBeg () {return &subtree_beg[0];}; 
// Compute functions
      /// Allocate memory
   void AllocateTree (int _nnodesmax, int _nlevmax) {
      father.resize (_nnodesmax+1); nchilds.resize (_nnodesmax+1); childs_list.resize (_nnodesmax+1);
      nodes_lev_id.resize (_nnodesmax+1); nnodes_lev.resize (_nlevmax+1); nodes_lev_list.resize (_nlevmax+1);
   }
   void FindOrderingOfNodesAccordingSubtrees (int *_ordernd); ///< Find ordering of nodes as subtrees
   void ApplyOrderingOfTreeNodes (int *_ordernd); ///< Apply ordering of nodes as subtrees
   void ComputeSubtreeStart (int *_subtree_start); ///< For each tree node find its subtree start
      /// Find subtree beg
   void InitSubtreeBeg () {
      this->subtree_beg.resize (this->nnodes+1); int *psubtree_beg = &this->subtree_beg[0]; this->ComputeSubtreeStart (psubtree_beg);
   }
   void PackTree (vector<char> &_obj); ///< Pack tree
   void UnPackTree (int _length, char *_obj); ///< Unpack tree
   int GetPackedTreeSize (); ///< Get packed tree size
   void FillPackedTree (int _length, char *_obj); ///< Fill packed tree
   void OutputTree (ostream &_stream); ///< Print tree
};

/// Class that supports base qrd
template <typename _FltVect>
class CQrdBase
{
protected:
   int nqblk;                        ///< Current number of q blocks
   int nqblk_alloc;                  ///< Allocated number of q blocks
   vector<int> qblksc;               ///< IA-type array that describes the number of columns in each q block
   vector<int> ncolarr_alloc;        ///< The list of allocate q column sizes 
   vector<int> nrowarr;              ///< The number of rows in each q block
   vector<int> nrowarr_alloc;        ///< The number of allocated rows in each q block
   vector<vector<_FltVect> > qarr;   ///< List of q blocks
   vector<vector<_FltVect> > tauarr; ///< List of tau blocks
public:
// Functions
// Constructors and destructor
   /// Default constructor
   CQrdBase () {
      nqblk = 0; nqblk_alloc = 0; qblksc.resize (2); qblksc[0] = 0; ncolarr_alloc.resize (1);
      nrowarr.resize (1); nrowarr_alloc.resize (1); qarr.resize (1); tauarr.resize (1);
   }; 
      /// Copy constructor
private:
   CQrdBase (const CQrdBase &_aa) {}; 
      /// Copy operator
   CQrdBase &operator= (const CQrdBase &_aa) {return *this; }; 
public:
      /// Destructor
   virtual ~CQrdBase () {}; 
// Operator functions
// Get/set functions
      /// Get nqblk
   int GetNqblk() const {return nqblk;}; 
      /// Get qblksc
   const int * GetQblksc() const {return &qblksc[0];}; 
      /// Get qblksc
   int * GetQblksc() {return &qblksc[0];}; 
      /// Get ncolarr_alloc
   const int * GetNcolarr_alloc() const {return &ncolarr_alloc[0];}; 
      /// Get ncolarr_alloc
   int * GetNcolarr_alloc() {return &ncolarr_alloc[0];}; 
      /// Get nrowarr
   const int * GetNrowarr() const {return &nrowarr[0];}; 
      /// Get nrowarr
   int * GetNrowarr() {return &nrowarr[0];}; 
      /// Get qarr
   const vector<_FltVect> * GetQarr() const {return &qarr[0];}; 
      /// Get qarr
   vector<_FltVect> * GetQarr() {return &qarr[0];}; 
      /// Get tauarr
   const vector<_FltVect> * GetTauarr() const {return &tauarr[0];}; 
      /// Get tauarr
   vector<_FltVect> * GetTauarr() {return &tauarr[0];}; 
      /// Set nqblk
   void SetNqblk (int _nqblk) { nqblk = _nqblk;}; 
// Qrd up level functions
   int GetNRows (); ///< Get number of rows
   int GetNCols (); ///< Get number of columns
   int GetAllocatedMemory (); ///< Get memory
      /// Update QR decomposition for current set of columns
   void UpdateQrdBlk (int _ncol, int _nrow, _FltVect *_ablk, int _lda);
      /// Multiply Q by given block (self version)
   void MvmQ (int _nrhs, _FltVect *_qx, int _ldqx);
      /// Multiply Q by given block
   void MvmQ (int _nrhs, int _nrows, _FltVect *_x, int _ldx, _FltVect *_qx, int _ldqx);
      /// Multiply Qh by given block (self version)
   void MvmQH (int _nrhs, _FltVect *_qx, int _ldqx);
      /// Multiply Qh by given block (self version)
   void MvmQHPart (int _nrhs, int _ibegQ, int _iendQ, _FltVect *_qx, int _ldqx);
      /// Multiply Qh by given block
   void MvmQH (int _nrhs, _FltVect *_x, int _ldx, _FltVect *_qx, int _ldqx);
      /// Get R part Qrd
   void GetRQrd (int _ibegc, int _iendc, int _ibegr, int _iendr,
                        _FltVect *_r, int _ldr); 
      /// Get R part Qrd (implementation)
   void GetRQrd_impl (int _ibegc, int _iendc, int _ibegr, int _iendr,
                        _FltVect *_r, int _ldr);
      /// Store R part Qrd
   void StoreR (int _ibegc, int _iendc, int _ibegr, int _iendr, _FltVect *_r, int _ldr); 
};

/// Class that supports QR decomposition with rows splitting
template <typename _FltVect>
class CQrdSet
{
protected:
   int nqrdsets;                      ///< number of subsets
   CQrdBase<_FltVect> *pqrd_head;     ///< pointer to head qrd structure
   CQrdBase<_FltVect> **ppqrd_childs; ///< set of pointers to childs of QR
public:
// Functions
// Constructors and destructor
      /// Default
   CQrdSet () {nqrdsets = 0;};
      /// Destructor
   virtual ~CQrdSet () {};
private:
      /// Copy
   CQrdSet (const CQrdSet &_aa) {};
      /// Copy operator
   CQrdSet &operator= (const CQrdSet &_aa) {return *this;}; 
public:
// Get/set interface functions
      /// Get nqrdsets
   int GetNqrdsets() const {return nqrdsets;}; 
      /// Get pqrd_head
   CQrdBase<_FltVect> * GetPQrdHead() const {return pqrd_head;}; 
      /// Get ppqrd_childs
   CQrdBase<_FltVect> ** GetPPQrdChilds() const {return ppqrd_childs;}; 
      /// Set nqrdsets
   void SetNqrdsets (int _nqrdsets) { nqrdsets = _nqrdsets;}; 
      /// Set pqrd_head
   void SetPQrdHead (CQrdBase<_FltVect> * _pqrd_head) { pqrd_head = _pqrd_head;}; 
      /// Set ppqrd_childs
   void SetPPQrdChilds (CQrdBase<_FltVect> ** _ppqrd_childs) { ppqrd_childs = _ppqrd_childs;}; 
// QrdSet up level functions
   void UpdateQrdHead (int _ncol); ///< Update head QR
      /// Multiply by head Q and prepare childs
   void MvmQ (int _nrhs, _FltVect *_x, int _ldx, _FltVect *_qx, int _ldqx);
};

/// Class that supports MPI QR decomposition
template <typename _FltVect>
class CQrdMPI 
{
protected:
   void *pComm;                       ///< Pointer to MPI comm
   int ni_myid;                       ///< Vector size of QRD for current processor
   CTree treeMPI;                     ///< MPI tree
   int node_myid_up;                  ///< Top level node for myid
   vector<int> nd2cpu;                ///< Number of cpu for each node
   vector<int> cpu2nd;                ///< Number of kernel node for each cpu
   CQrdSet<_FltVect> *pqrd_sets;      ///< Qrd sets
   CQrdBase<_FltVect> *pqrd_childs;   ///< Qrd childs
   CQrdBase<_FltVect> **ppqrd_childs; ///< Pointers to Qrd childs
public:
// Functions
// Constructors and destructor
      /// Default
   CQrdMPI () {
      pComm = NULL; ni_myid = 0; node_myid_up = -1; nd2cpu.resize (1); cpu2nd.resize (1);
      pqrd_sets = NULL; pqrd_childs = NULL; ppqrd_childs = NULL;
   }; 
      /// Destructor
   virtual ~CQrdMPI () {
      if (pqrd_sets != NULL) delete [] pqrd_sets;
      if (pqrd_childs != NULL) delete [] pqrd_childs;
      if (ppqrd_childs != NULL) delete [] ppqrd_childs;
   }; 
private:
      /// Copy
   CQrdMPI (const CQrdMPI &_aa) {};
      /// Copy operator
   CQrdMPI &operator= (const CQrdMPI &_aa) {return *this;}; 
public:
// Get/set interface functions
      /// Get number of rows
   int GetNRows () {return ni_myid;}; 
      /// Get tree
   CTree * GetTree () {return &treeMPI;}; 
      /// Get qrd_sets
   CQrdSet<_FltVect> * GetQrdSets () {return pqrd_sets;}; 
      /// Get qrd_childs
   CQrdBase<_FltVect> * GetQrdChilds () {return pqrd_childs;}; 
// Qrd up level functions
   void Init (void *_pComm, int _ni_myid); ///< Init the structure
   int GetNCols (); ///< Get the current number of columns
      /// Update QR decomposition
   void UpdateQrdMPI (int _ncol, int _nrows, _FltVect *_ablk, int _lda); 
      /// Multiply by Q
   void MvmQMPI (int _nrhs, _FltVect *_x, int _ldx, _FltVect *_qx, int _ldqx); 
      /// Get R part
   void GetRQrdMPI (int _ibegc, int _iendc, int _ibegr, int _iendr,
                        _FltVect *_r, int _ldr);
      /// Free qrd
   void FreeQrdMPI ();
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
      /// Get nblks
   int GetNiCpu () {return ni_cpu;}; 
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
   int ncoef_slv;                                 ///< Number of coefs in poly prec
   vector<double> coef_slv;                       ///< Coefs in poly prec
   vector<double> xwork;                          ///< Work memory for polynomial preconditioning
// Constructors and destructor
private:
      /// Copy constructor
   CSolver (const CSolver<_Int,_Flt,_FltVect> &_aa) {}; 
      /// Copy operator
   CSolver<_Int,_Flt,_FltVect> &operator= (const CSolver<_Int,_Flt,_FltVect> &_aa) { return *this; }; 
public:
      /// Epmty constructor
   CSolver () {nblks = 0; ncoef_slv = 0; /*cout<<"CSolver: empty constructor\n";*/}; 
      /// Destructor
   ~CSolver () {/*cout<<"CSolver: empty destructor\n";*/};
// Get/set functions
      /// Get comm
   void * GetComm () { return pcomm; };
      /// Get nblks
   int GetNblks () {return nblks;}; 
      /// Get blks
   long long * GetBlks () { return &blks[0]; };
      /// Get blk2cpu
   int * GetBlk2cpu () { return &blk2cpu[0]; };
      /// Set ncoef_slv
   void SetNcoef (int _ncoef) { ncoef_slv = _ncoef; };
      /// Get coef
   vector<double> * GetCoef () { return &coef_slv; };
// External functions
      /// Prepare solver structures including performing parallel fct
   void PrepareSolver (void *_pcomm, SSolverParams *_params,
                        int _nblks, long long *_blks, int *_blk2cpu,
                        bool _b_store_matrix, _Int *_ia, _Int *_ja, _Flt *_a,
                        double &_prec_extend, double &_density, double &_scpiv_min, double &_scpiv_max, 
                        int &_nmodif, double &_piv_min, double &_piv_max, double &_dtime_fct) 
   {
      PrepareMatrix (_pcomm, _nblks, _blks, _blk2cpu, _ia, _ja, _a);
      ComputeBILU2 (_params, _prec_extend, _density, _scpiv_min, _scpiv_max, 
                     _nmodif, _piv_min, _piv_max, _dtime_fct);
      if (!_b_store_matrix) CleanMvmA ();
   };
      /// Prepare matrix data only
   void PrepareMatrix (void *_pcomm, int _nblks, long long *_blks, int *_blk2cpu,
                        _Int *_ia, _Int *_ja, _Flt *_a);
      /// Peform BILU2 computations
   void ComputeBILU2 (SSolverParams *_params,
                        double &_prec_extend, double &_density, double &_scpiv_min, double &_scpiv_max, 
                        int &_nmodif, double &_piv_min, double &_piv_max, double &_dtime_fct);
      /// Perform iterations of the iterative scheme
   void SolveIter (int _ittype, int _niter_max, int _niter_cycle, int _ncoef, int _niter_cycle2, 
                     double _eps, int _ichk, int _msglev, ofstream *_fout, 
                     _FltVect *_rhs, _FltVect *_sol,
                     double &_rhs_norm, double &_res_ini, 
                     int &_niter, int &_nmvm, double &_res_fin, double &_dtime_iter);
      /// Perform iterations of the BiCGStab iterative scheme
   void BiCGStab (bool _b_use_poly,
                  int _niter_max, double _eps, int _ichk, int _msglev, ofstream *_fout, 
                  _FltVect *_rhs, _FltVect *_sol,
                  double &_rhs_norm, double &_res_ini, 
                  int &_niter, int &_nmvm, double &_res_fin, double &_dtime_iter);
      /// Perform iterations of the GMRES iterative scheme
   void Gmres (bool _b_use_poly, int _niter_max, int _niter_cycle, int _ncoef, double _eps, int _ichk, int _msglev, ofstream *_fout, 
                  _FltVect *_rhs, _FltVect *_sol,
                  double &_rhs_norm, double &_res_ini, 
                  int &_niter, int &_nmvm, double &_res_fin, double &_dtime_iter);
      /// Clean MvmA structures
   void CleanMvmA () {
      mvm.Clean(); vector<CHMatrix<_Int,_Flt> > hmatr_arr_dummy; hmatr_arr.swap (hmatr_arr_dummy);
   };
      /// Clean Slv structures
   void CleanSlv () {
      vector<int> nlist_ext_arr_dummy; vector<vector<int> > list_ext_arr_dummy; vector<vector<int> > order_LU_dummy;
      vector<CMatrix<_Int,float>  > matrL_float_dummy; vector<CMatrix<_Int,double> > matrL_double_dummy;
      vector<CMatrix<_Int,float>  > matrU_float_dummy; vector<CMatrix<_Int,double> > matrU_double_dummy;
      vector<double> coef_temp; vector<double> xwork_temp;
      nlist_ext_arr.swap (nlist_ext_arr_dummy); list_ext_arr.swap (list_ext_arr_dummy); order_LU.swap (order_LU_dummy);
      matrL_float.swap (matrL_float_dummy); matrL_double.swap (matrL_double_dummy); matrU_float.swap (matrU_float_dummy); matrU_double.swap (matrU_double_dummy); 
      slv_float.Clean(); slv_double.Clean(); ncoef_slv = 0; coef_slv.swap (coef_temp); xwork.swap (xwork_temp);
   };
      /// Clean all structures
   void Clean () {
      pcomm = NULL; nblks = 0; vector<long long> blks_dummy; vector<long long> blks_ext_dummy; vector<int> blk2cpu_dummy;
      blks.swap (blks_dummy); blks_ext.swap (blks_ext_dummy); blk2cpu.swap (blk2cpu_dummy); CleanMvmA(); CleanSlv(); params.SetDefaults(); 
   };
// Internal functions
      /// Perform MvmA parallel computations
   void MvmA (const _FltVect *_x, _FltVect *_ax) {mvm.MvmA(_x,_ax);};
      /// Perform polynomial preconditioning
   void PolySlvLU (bool _b_use_poly, const _FltVect *_x, _FltVect *_px);
      /// Perform SlvLU parallel computations
   void SlvLU (const _FltVect *_x, _FltVect *_px) {if (params.prec_float == 1) { slv_float.SolveLU(_x,_px); } else { slv_double.SolveLU(_x,_px); }; };
};

} // namespace k3d

#endif
