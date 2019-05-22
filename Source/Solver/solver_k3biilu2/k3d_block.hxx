//------------------------------------------------------------------------------------------------
// File: k3d_block.hxx
//------------------------------------------------------------------------------------------------

#ifndef __K3D_BLOCK_HXX
#define __K3D_BLOCK_HXX

#include "k3d_slv_thr.hxx"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

namespace k3d
{

   template < typename _Flt > class CBlock_BxB_traits   // Class that supports bxb dense blocks operations
   {
    public:
      /// Matrix by matrix square multiplication functions:
      static void MM_BxB (int _b_size, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AxB_blk);      ///< Multiply C  = A*B
      static void MM_Add_BxB (int _b_size, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AxB_blk);  ///< Multiply C += A*B
      static void MM_Sub_BxB (int _b_size, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AxB_blk);  ///< Multiply C -= A*B
      static void MtM_BxB (int _b_size, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AtxB_blk);    ///< Multiply C  = At*B
      static void MtM_Add_BxB (int _b_size, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AtxB_blk);        ///< Multiply C += At*B
      static void MtM_Sub_BxB (int _b_size, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AtxB_blk);        ///< Multiply C -= At*B
      static void MMt_BxB (int _b_size, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AxBt_blk);    ///< Multiply C  = A*Bt
      static void MMt_Add_BxB (int _b_size, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AxBt_blk);        ///< Multiply C += A*Bt
      static void MMt_Sub_BxB (int _b_size, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AxBt_blk);        ///< Multiply C -= A*Bt
      /// Square matrix by rect matrix multiplication functions:
      static void MM_BxB (int _b_size, int _ncol, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AxB_blk);   ///< Multiply C  = A*B
      static void MM_Add_BxB (int _b_size, int _ncol, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AxB_blk);       ///< Multiply C += A*B
      static void MM_Sub_BxB (int _b_size, int _ncol, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AxB_blk);       ///< Multiply C -= A*B
      static void MtM_BxB (int _b_size, int _ncol, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AtxB_blk); ///< Multiply C  = At*B
      static void MtM_Add_BxB (int _b_size, int _ncol, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AtxB_blk);     ///< Multiply C += At*B
      static void MtM_Sub_BxB (int _b_size, int _ncol, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AtxB_blk);     ///< Multiply C -= At*B
      /// Rect matrix by rect matrix multiplication functions:
      static void MM_BxB (int _b_size, int _b_size2, int _ncol, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AxB_blk);     ///< Multiply C  = A*B
      static void MM_Add_BxB (int _b_size, int _b_size2, int _ncol, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AxB_blk); ///< Multiply C += A*B
      static void MM_Sub_BxB (int _b_size, int _b_size2, int _ncol, const _Flt * _A_blk, const _Flt * _B_blk, _Flt * _AxB_blk); ///< Multiply C -= A*B
      /// Transpose block
      static void Transp_BxB (int _b_size, const _Flt * _A_blk, _Flt * _At_blk);
      /// Compute symmetric triangular fct
      static int FctSymm_BxB (double _piv_min, int _b_size, const _Flt * _AU_blk,
                              _Flt * _U_blk, double &_dia_min, double &_dia_max);
      /// Compute triangular fct
      static int Fct_BxB (double _piv_min, int _b_size, const _Flt * _A_blk,
                          _Flt * _Lt_blk, _Flt * _U_blk, double &_dia_min,
                          double &_dia_max);
      /// Inverse upper triangular
      static int InvU_BxB (int _dia_type, double _piv_min, int _b_size, _Flt * _U_blk,
                           _Flt * _U_inv_blk, _Flt * _Work, double *_dWork,
                           double &_dia_min, double &_dia_max)
      {
         int n_modif = 0;
         if (_dia_type == 0) {
            n_modif =
               CBlock_BxB_traits < _Flt >::InvU_BxB (_piv_min, _b_size, _U_blk,
                                                     _U_inv_blk, _dia_min, _dia_max);
         } else
         {
            n_modif =
               CBlock_BxB_traits < _Flt >::InvA_Svd_BxB (_piv_min, _b_size, _U_blk,
                                                         _U_inv_blk, _Work, _dWork,
                                                         _dia_min, _dia_max);
         }
         return n_modif;
      };
      /// Inverse upper triangular
      static void InvU_BxB (int _b_size, const _Flt * _U_blk, _Flt * _U_inv_blk);
      /// Inverse upper triangular
      static int InvU_BxB (double _piv_min, int _b_size, _Flt * _U_blk, _Flt * _U_inv_blk,
                           double &_dia_min, double &_dia_max);
      /// Inverse dense square matrix
      static int InvA_Svd_BxB (double _piv_min, int _b_size, _Flt * _U_blk,
                               _Flt * _U_inv_blk, _Flt * _Work, double *_dWork,
                               double &_dia_min, double &_dia_max);
      /// Compute symmetric triangular fct and its inverse
      static int FctSymmInv_BxB (double _piv_min, int _b_size, const _Flt * _AU_blk,
                                 _Flt * _U_blk, _Flt * _U_inv_blk, double &_dia_min,
                                 double &_dia_max)
      {
         int n_modif;
         n_modif =
            CBlock_BxB_traits < _Flt >::FctSymm_BxB (_piv_min, _b_size, _AU_blk, _U_blk,
                                                     _dia_min, _dia_max);
         CBlock_BxB_traits < _Flt >::InvU_BxB (_b_size, _U_blk, _U_inv_blk);
         return n_modif;
      }
      /// Compute triangular fct and its inverse
      static int FctInv_BxB (double _piv_min, int _b_size, const _Flt * _A_blk,
                             _Flt * _Lt_blk, _Flt * _U_blk, _Flt * _Lt_inv_blk,
                             _Flt * _U_inv_blk, double &_dia_min, double &_dia_max)
      {
         int n_modif;
         n_modif =
            CBlock_BxB_traits < _Flt >::Fct_BxB (_piv_min, _b_size, _A_blk, _Lt_blk,
                                                 _U_blk, _dia_min, _dia_max);
         CBlock_BxB_traits < _Flt >::InvU_BxB (_b_size, _Lt_blk, _Lt_inv_blk);
         CBlock_BxB_traits < _Flt >::InvU_BxB (_b_size, _U_blk, _U_inv_blk);
         return n_modif;
      }
      /// Compute fct and its inverse via two-sided QR
      static int FctInv_BiQrd_BxB (double _piv_min, int _b_size, const _Flt * _A_blk,
                                   _Flt * _Lt_blk, _Flt * _U_blk, _Flt * _Lt_inv_blk,
                                   _Flt * _U_inv_blk, _Flt * _Work, double &_dia_min,
                                   double &_dia_max);
      /// Compute fct and its inverse via SVD
      static int FctInv_Svd_BxB (double _piv_min, int _b_size, const _Flt * _A_blk,
                                 _Flt * _Lt_blk, _Flt * _U_blk, _Flt * _Lt_inv_blk,
                                 _Flt * _U_inv_blk, _Flt * _Work, double *_dWork,
                                 double &_dia_min, double &_dia_max);
      /// Modify diagonal when filtering
      static void ModifDia_BxB (int _b_size, double _theta, const _Flt * _L_JI_offd,
                                const _Flt * _U_IJ_offd, _Flt * _dia_I, _Flt * _dia_J);
      /// Compute norm of BxB block
      static double FNormValue_BxB (int _b_size, const _Flt * _A);
      /// Compute norm of BxB block
      static double FNormValue_BxB (int _b_size, const _Flt * _L_offd,
                                    const _Flt * _U_offd);
      /// Add diagonal modification elements
      static void AddModifDia_BxB (int _b_size, _Flt * _dia_modif, _Flt * _dia_U);
      /// Compute Svd
      static void ComputeSvd_BxB (int _b_size, const _Flt * _a_matr, _Flt * _Sv_arr,
                                  _Flt * _U_matr, _Flt * _V_matr, double *_dWork);
      /// Compute Qrd in column format
      static void ComputeQrd_BxB (int _m, int _b_size, const _Flt * _A_matr,
                                  _Flt * _Q_matr, _Flt * _R_matr, _Flt * _Work);
      /// Two sided bidiagonalization
      static void BiDiagonalize_BxB (int _n, const _Flt * _amatr, _Flt * _bidiag,
                                     _Flt * _uarr, _Flt * _vharr, _Flt * _Work);
   };

   template < typename _Int, typename _Flt > class CFct_bxb_impl        // Class that supports bxb blocks factorization computations
   {
    public:
      /// Find ordering as splitting diagonal blocks with very small pivot/singular values
      static void OrderSplitDiaPairs (SParams * _params, int _b_size, int _n,
                                      const vector < _Int > &_ia_au,
                                      const vector < _Int > &_ja_au,
                                      const vector < char >&_ja_char_au,
                                      const vector < _Flt > &_a_au, int *_order,
                                      int &_n1);
      /// Split matrix into submatrices by rows (pairs or not)
      static void SubmatricesByRows (bool _b_is_char, bool _b_is_pair, int _b_size,
                                     int _n, int _n1, const vector < _Int > &_ia,
                                     const vector < _Int > &_ja,
                                     const vector < char >&_ja_char,
                                     const vector < _Flt > &_a, vector < _Int > &_ia_ini,
                                     vector < _Int > &_ja_ini,
                                     vector < char >&_ja_char_ini,
                                     vector < _Flt > &_a_ini, vector < _Int > &_ia_last,
                                     vector < _Int > &_ja_last,
                                     vector < char >&_ja_char_last,
                                     vector < _Flt > &_a_last);
      /// Add sparse matrices with elements (pairs or not)
      static void AddMatrices (bool _b_is_char, bool _b_is_pair, int _b_size, int _n,
                               const vector < _Int > &_ia_1, const vector < _Int > &_ja_1,
                               const vector < char >&_ja_char_1,
                               const vector < _Flt > &_a_1, const vector < _Int > &_ia_2,
                               const vector < _Int > &_ja_2,
                               const vector < char >&_ja_char_2,
                               const vector < _Flt > &_a_2, vector < _Int > &_ia_sum,
                               vector < _Int > &_ja_sum, vector < char >&_ja_char_sum,
                               vector < _Flt > &_a_sum);
      /// Compute condensed block matrix
      static void CondenseMatrix (int _n_pt, const vector < _Int > &_ia_pt,
                                  const vector < _Int > &_ja_pt,
                                  const vector < _Flt > &_a_pt, int _b_size,
                                  vector < _Int > &_ia_cnd, vector < _Int > &_ja_cnd,
                                  vector < _Flt > &_a_cnd);
      /// Compute condensed block rectangular matrix
      static void CondenseRectMatrix (int _nlist_pt, const vector < _Int > &_list_pt,
                                      const vector < _Int > &_ia_pt,
                                      const vector < _Int > &_ja_pt,
                                      const vector < _Flt > &_a_pt, int _b_size,
                                      int &_nlist_cnd, vector < _Int > &_list_cnd,
                                      vector < _Int > &_ia_cnd, vector < _Int > &_ja_cnd,
                                      vector < _Flt > &_a_cnd, int _nimax, int &_icycle,
                                      int *_imask);
      /// Compute ordered block matrix
      static void ReorderMatrix (bool _b_is_char, bool _b_is_pair, int _b_size, int _n,
                                 const vector < int >&_order,
                                 const vector < _Int > &_ia_alu,
                                 const vector < _Int > &_ja_alu,
                                 const vector < char >&_ja_char_alu,
                                 const vector < _Flt > &_a_alu,
                                 vector < _Int > &_ia_alu_ord,
                                 vector < _Int > &_ja_alu_ord,
                                 vector < char >&_ja_char_alu_ord,
                                 vector < _Flt > &_a_alu_ord);
      /// Compute block scaling
      static void ComputeScaling (int _sctype, int _nitersc, double _scl_min, int _b_size,
                                  int _n, const vector < _Int > &_ia_alu,
                                  const vector < _Int > &_ja_alu,
                                  const vector < _Flt > &_a_alu, vector < _Flt > &_sclL,
                                  vector < _Flt > &_sclL_inv, vector < _Flt > &_sclU,
                                  vector < _Flt > &_sclU_inv, int &_nmodif,
                                  double &_sclmin_att, double &_sclmax_att);
      /// Perform explicit block scaling
      static void MatrixScale (int _b_size, int _n, const vector < _Flt > &_sclL,
                               const vector < _Flt > &_sclU,
                               const vector < _Int > &_ia_alu,
                               const vector < _Int > &_ja_alu, vector < _Flt > &_a_alu);
      /// Perform explicit block scaling
      static void MatrixScale (int _b_size, int _nlist, const _Flt * _sclL,
                               const _Flt * _sclU, const vector < _Int > &_list_alu,
                               const vector < _Int > &_ia_alu,
                               const vector < _Int > &_ja_alu, vector < _Flt > &_a_alu);
      /// Split block matrix data into L and U parts
      static void SplitLU (bool _b_is_char, int _b_size, int _n,
                           const vector < _Int > &_ia_alu, const vector < _Int > &_ja_alu,
                           const vector < char >&_jachar_alu,
                           const vector < _Flt > &_a_alu, vector < _Int > &_ia_l,
                           vector < _Int > &_ja_l, vector < char >&_jachar_l,
                           vector < _Flt > &_a_l, vector < _Int > &_ia_u,
                           vector < _Int > &_ja_u, vector < char >&_jachar_u,
                           vector < _Flt > &_a_u);
      /// Transpose square block matrix
      static void Transpose (bool _b_is_char, int _b_size, int _n,
                             const vector < _Int > &_ia_a, const vector < _Int > &_ja_a,
                             const vector < char >&_jachar_a, const vector < _Flt > &_a_a,
                             vector < _Int > &_ia_at, vector < _Int > &_ja_at,
                             vector < char >&_jachar_at, vector < _Flt > &_a_at);
      /// Combine block L and U
      static void CombineLU (bool _b_is_char, int _b_size, int _n,
                             const vector < _Int > &_ia_al, const vector < _Int > &_ja_al,
                             const vector < char >&_jachar_al,
                             const vector < _Flt > &_a_al, const vector < _Int > &_ia_au,
                             const vector < _Int > &_ja_au,
                             const vector < char >&_jachar_au,
                             const vector < _Flt > &_a_au, vector < _Int > &_ia_alu,
                             vector < _Int > &_ja_alu, vector < char >&_jachar_alu,
                             vector < _Flt > &_a_alu);
      /// Combine block L and U data into extended block pairs
      static void CombineLUPairs (bool _b_is_char, int _b_size, int _n,
                                  const vector < _Int > &_ia_l,
                                  const vector < _Int > &_ja_l,
                                  const vector < char >&_jachar_l,
                                  const vector < _Flt > &_a_l,
                                  const vector < _Int > &_ia_u,
                                  const vector < _Int > &_ja_u,
                                  const vector < char >&_jachar_u,
                                  const vector < _Flt > &_a_u, vector < _Int > &_ia_alu,
                                  vector < _Int > &_ja_alu, vector < char >&_jachar_alu,
                                  vector < _Flt > &_a_alu);
      /// Split block pairs fct data into block L and U parts
      static void SplitLUPairs (bool _b_is_char, int _b_size, int _n,
                                const vector < _Int > &_ia_alu,
                                const vector < _Int > &_ja_alu,
                                const vector < char >&_jachar_alu,
                                const vector < _Flt > &_a_alu, vector < _Int > &_ia_l,
                                vector < _Int > &_ja_l, vector < char >&_jachar_l,
                                vector < _Flt > &_a_l, vector < _Int > &_ia_u,
                                vector < _Int > &_ja_u, vector < char >&_jachar_u,
                                vector < _Flt > &_a_u);
      /// Rescale block factor back
      static void RescaleU (int _b_size, int _n, const vector < _Flt > &_sclU,
                            const vector < _Flt > &_invsclU, const vector < _Int > &_ia_u,
                            const vector < _Int > &_ja_u, vector < _Flt > &_a_u);
      /// Rescale block factor back
      static void RescaleU (char _type, int _b_size, int _nlist, const _Flt * _sclU,
                            const _Flt * _invsclU, const vector < _Int > &_list_u,
                            const vector < _Int > &_ia_u, const vector < _Int > &_ja_u,
                            vector < _Flt > &_a_u);
      /// Perform ICH2 small dense blocks sparse factorization of the block with future block diagonal modifications (with matrix degree structural control)
      static void FctBlockIch2Degree (SParams * _params, int _b_size, int _n, int _n_ini,
                                      const vector < _Int > &_ia_alu,
                                      const vector < _Int > &_ja_alu,
                                      const vector < char >&_ja_char_alu,
                                      const vector < _Flt > &_a_alu,
                                      vector < _Int > &_ia_lu, vector < _Int > &_ja_lu,
                                      vector < char >&_ja_char_lu, vector < _Flt > &_a_lu,
                                      int &_nmodif, double &_eigmin_att,
                                      double &_eigmax_att);
      /// Perform ILU2 small dense blocks sparse factorization of the block with future block diagonal modifications (with matrix degree structural control)
      static void FctBlockIlu2Degree (SParams * _params, int _b_size, int _n, int _n_ini,
                                      const vector < _Int > &_ia_alu,
                                      const vector < _Int > &_ja_alu,
                                      const vector < char >&_ja_char_alu,
                                      const vector < _Flt > &_a_alu,
                                      vector < _Int > &_ia_lu, vector < _Int > &_ja_lu,
                                      vector < char >&_ja_char_lu, vector < _Flt > &_a_lu,
                                      int &_nmodif, double &_eigmin_att,
                                      double &_eigmax_att);
      /// Perform one step fct for dinamically splitted and ordered triangular matrix
      static void DiaSplitFctSchur (SParams * _params, double _dia_split, int _b_size,
                                    int _n, const vector < _Int > &_ia_au,
                                    const vector < _Int > &_ja_au,
                                    const vector < char >&_ja_char_au,
                                    const vector < _Flt > &_a_au, int *_order,
                                    int &_n_ini, vector < _Int > &_ia_u,
                                    vector < _Int > &_ja_u, vector < char >&_ja_char_u,
                                    vector < _Flt > &_a_u, vector < _Int > &_ia_u_schur,
                                    vector < _Int > &_ja_u_schur,
                                    vector < char >&_ja_char_u_schur,
                                    vector < _Flt > &_a_u_schur, int &_nmodif,
                                    double &_eigmin_att, double &_eigmax_att);
      /// Perform one step fct for dinamically splitted and ordered triangular matrix
      static void DiaSplitFctSchur (SParams * _params, double _dia_split, int _b_size,
                                    int _n_row, int _n, int _n_ini_max,
                                    const vector < _Int > &_ia_au,
                                    const vector < _Int > &_ja_au,
                                    const vector < char >&_ja_char_au,
                                    const vector < _Flt > &_a_au, int *_order,
                                    int &_n_ini, int &_n_bad_diag, vector < _Int > &_ia_u,
                                    vector < _Int > &_ja_u, vector < char >&_ja_char_u,
                                    vector < _Flt > &_a_u, vector < _Int > &_ia_u_schur,
                                    vector < _Int > &_ja_u_schur,
                                    vector < char >&_ja_char_u_schur,
                                    vector < _Flt > &_a_u_schur, int &_nmodif,
                                    double &_eigmin_att, double &_eigmax_att);
      /// Perform multi step fct for dinamically splitted and ordered triangular matrix
      static void FctMLevDiaSplit (SParams * _params, int _nlev, double *_dia_lev,
                                   int _b_size, int _n, const vector < _Int > &_ia_au,
                                   const vector < _Int > &_ja_au,
                                   const vector < char >&_ja_char_au,
                                   const vector < _Flt > &_a_au, int *_order,
                                   int &_nblks_fct, int *_blks_fct,
                                   vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                                   vector < char >&_ja_char_u, vector < _Flt > &_a_u,
                                   int &_nmodif, double &_eigmin_att,
                                   double &_eigmax_att);
      /// Perform multi step fct for dinamically splitted and ordered extended triangular matrix
      static void FctMLevDiaSplit (SParams * _params, int _nlev, double *_dia_lev,
                                   int _b_size, int _n, int _n_ini,
                                   const vector < _Int > &_ia_au,
                                   const vector < _Int > &_ja_au,
                                   const vector < char >&_ja_char_au,
                                   const vector < _Flt > &_a_au, int *_order,
                                   int &_n_moved_diag, int &_nblks_fct, int *_blks_fct,
                                   vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                                   vector < char >&_ja_char_u, vector < _Flt > &_a_u,
                                   int &_nmodif, double &_eigmin_att,
                                   double &_eigmax_att);
   };

   template < typename _Int, typename _Flt > class CFct_bxb     // Class that supports factorization computations
   {
    public:
// Functions that support matrix structure as a whole
      /// Compute ordered matrix
      static void ReorderMatrix_BxB (bool _b_is_char, bool _b_is_pair, int _blksize,
                                     const CMatrix < _Int, _Flt > &_a_matr,
                                     vector < int >&_order, CMatrix < _Int,
                                     _Flt > &_a_ord)
      {
         _a_ord.SetN (_a_matr.GetNlist ());
         _a_ord.SetNzja (_a_matr.GetNzja ());
         _a_ord.SetNza (_a_matr.GetNza ());
         CFct_bxb_impl < _Int, _Flt >::ReorderMatrix (_b_is_char, _b_is_pair, _blksize,
                                                      _a_matr.GetNlist (), _order,
                                                      *(_a_matr.GetIa ()),
                                                      *(_a_matr.GetJa ()),
                                                      *(_a_matr.GetJaChar ()),
                                                      *(_a_matr.GetA ()),
                                                      *(_a_ord.GetIa ()),
                                                      *(_a_ord.GetJa ()),
                                                      *(_a_ord.GetJaChar ()),
                                                      *(_a_ord.GetA ()));
         _a_ord.SetIdentityList ();
      }
      /// Compute block scaling
      static void ComputeScaling (int _sctype, int _nitersc, double _scl_min, int _b_size,
                                  int _n, CMatrix < _Int, _Flt > &_a_matr,
                                  vector < _Flt > &_sclL, vector < _Flt > &_sclL_inv,
                                  vector < _Flt > &_sclU, vector < _Flt > &_sclU_inv,
                                  int &_nmodif, double &_sclmin_att, double &_sclmax_att)
      {
         CFct_bxb_impl < _Int, _Flt >::ComputeScaling (_sctype, _nitersc, _scl_min,
                                                       _b_size, _n, *(_a_matr.GetIa ()),
                                                       *(_a_matr.GetJa ()),
                                                       *(_a_matr.GetA ()), _sclL,
                                                       _sclL_inv, _sclU, _sclU_inv,
                                                       _nmodif, _sclmin_att, _sclmax_att);
      };
      /// Perform explicit block scaling
      static void MatrixScale (int _b_size, int _n, const vector < _Flt > &_sclL,
                               const vector < _Flt > &_sclU, CMatrix < _Int,
                               _Flt > &_a_matr)
      {
         CFct_bxb_impl < _Int, _Flt >::MatrixScale (_b_size, _n, _sclL, _sclU,
                                                    *(_a_matr.GetIa ()),
                                                    *(_a_matr.GetJa ()),
                                                    *(_a_matr.GetA ()));
      }
      /// Split block matrix data into L and U parts
      static void SplitLU (bool _b_is_char, int _b_size, int _n, const CMatrix < _Int,
                           _Flt > &_a_matr, CMatrix < _Int, _Flt > &_al_matr,
                           CMatrix < _Int, _Flt > &_au_matr)
      {
         CFct_bxb_impl < _Int, _Flt >::SplitLU (_b_is_char, _b_size, _n,
                                                *(_a_matr.GetIa ()), *(_a_matr.GetJa ()),
                                                *(_a_matr.GetJaChar ()),
                                                *(_a_matr.GetA ()), *(_al_matr.GetIa ()),
                                                *(_al_matr.GetJa ()),
                                                *(_al_matr.GetJaChar ()),
                                                *(_al_matr.GetA ()), *(_au_matr.GetIa ()),
                                                *(_au_matr.GetJa ()),
                                                *(_au_matr.GetJaChar ()),
                                                *(_au_matr.GetA ()));
         vector < _Int > *pia_al = _al_matr.GetIa ();
         int nzja_al = (int) (*pia_al)[_n];
         _al_matr.SetNlist (_n);
         _al_matr.SetNzja (nzja_al);
         if (_b_is_char)
            _al_matr.SetNzjaChar (nzja_al);
         _al_matr.SetNza (nzja_al * _b_size * _b_size);
         _al_matr.SetIdentityList ();
         vector < _Int > *pia_au = _au_matr.GetIa ();
         int nzja_au = (int) (*pia_au)[_n];
         _au_matr.SetNlist (_n);
         _au_matr.SetNzja (nzja_au);
         if (_b_is_char)
            _au_matr.SetNzjaChar (nzja_au);
         _au_matr.SetNza (nzja_au * _b_size * _b_size);
         _au_matr.SetIdentityList ();
      }
      /// Transpose square block matrix
      static void Transpose (bool _b_is_char, int _b_size, int _n, const CMatrix < _Int,
                             _Flt > &_a_matr, CMatrix < _Int, _Flt > &_at_matr)
      {
         CFct_bxb_impl < _Int, _Flt >::Transpose (_b_is_char, _b_size, _n,
                                                  *(_a_matr.GetIa ()),
                                                  *(_a_matr.GetJa ()),
                                                  *(_a_matr.GetJaChar ()),
                                                  *(_a_matr.GetA ()),
                                                  *(_at_matr.GetIa ()),
                                                  *(_at_matr.GetJa ()),
                                                  *(_at_matr.GetJaChar ()),
                                                  *(_at_matr.GetA ()));
         _at_matr.SetNlist (_a_matr.GetNlist ());
         _at_matr.SetNzja (_a_matr.GetNzja ());
         if (_b_is_char)
            _at_matr.SetNzjaChar (_a_matr.GetNzjaChar ());
         _at_matr.SetNza (_a_matr.GetNza ());
         _at_matr.SetIdentityList ();
      };
      /// Combine block L and U data into extended block pairs
      static void CombineLUPairs (bool _b_is_char, int _b_size, int _n,
                                  const CMatrix < _Int, _Flt > &_al_matr,
                                  const CMatrix < _Int, _Flt > &_au_matr, CMatrix < _Int,
                                  _Flt > &_alu_matr)
      {
         CFct_bxb_impl < _Int, _Flt >::CombineLUPairs (_b_is_char, _b_size, _n,
                                                       *(_al_matr.GetIa ()),
                                                       *(_al_matr.GetJa ()),
                                                       *(_al_matr.GetJaChar ()),
                                                       *(_al_matr.GetA ()),
                                                       *(_au_matr.GetIa ()),
                                                       *(_au_matr.GetJa ()),
                                                       *(_au_matr.GetJaChar ()),
                                                       *(_au_matr.GetA ()),
                                                       *(_alu_matr.GetIa ()),
                                                       *(_alu_matr.GetJa ()),
                                                       *(_alu_matr.GetJaChar ()),
                                                       *(_alu_matr.GetA ()));
         vector < _Int > *pia_alu = _alu_matr.GetIa ();
         int nzja_alu = (int) (*pia_alu)[_n];
         _alu_matr.SetNlist (_n);
         _alu_matr.SetNzja (nzja_alu);
         if (_b_is_char)
            _alu_matr.SetNzjaChar (nzja_alu);
         _alu_matr.SetNza (nzja_alu * _b_size * _b_size * 2);
         _alu_matr.SetIdentityList ();
      };
      /// Split block pairs fct data into block L and U parts
      static void SplitLUPairs (bool _b_is_char, int _b_size, int _n,
                                const CMatrix < _Int, _Flt > &_a_matr, CMatrix < _Int,
                                _Flt > &_al_matr, CMatrix < _Int, _Flt > &_au_matr)
      {
         CFct_bxb_impl < _Int, _Flt >::SplitLUPairs (_b_is_char, _b_size, _n,
                                                     *(_a_matr.GetIa ()),
                                                     *(_a_matr.GetJa ()),
                                                     *(_a_matr.GetJaChar ()),
                                                     *(_a_matr.GetA ()),
                                                     *(_al_matr.GetIa ()),
                                                     *(_al_matr.GetJa ()),
                                                     *(_al_matr.GetJaChar ()),
                                                     *(_al_matr.GetA ()),
                                                     *(_au_matr.GetIa ()),
                                                     *(_au_matr.GetJa ()),
                                                     *(_au_matr.GetJaChar ()),
                                                     *(_au_matr.GetA ()));
         vector < _Int > *pia_al = _al_matr.GetIa ();
         int nzja_al = (int) (*pia_al)[_n];
         _al_matr.SetNlist (_n);
         _al_matr.SetNzja (nzja_al);
         if (_b_is_char)
            _al_matr.SetNzjaChar (nzja_al);
         _al_matr.SetNza (nzja_al * _b_size * _b_size);
         _al_matr.SetIdentityList ();
         vector < _Int > *pia_au = _au_matr.GetIa ();
         int nzja_au = (int) (*pia_au)[_n];
         _au_matr.SetNlist (_n);
         _au_matr.SetNzja (nzja_au);
         if (_b_is_char)
            _au_matr.SetNzjaChar (nzja_au);
         _au_matr.SetNza (nzja_au * _b_size * _b_size);
         _au_matr.SetIdentityList ();
      };
      /// Perform ILU2 small dense blocks sparse factorization of the block with future block diagonal modifications (with matrix degree structural control)
      static void FctBlockIlu2Degree (SParams * _params, int _b_size, int _n, int _n_ini,
                                      const CMatrix < _Int, _Flt > &_alu_matr,
                                      CMatrix < _Int, _Flt > &_lu_matr, int &_nmodif,
                                      double &_eigmin_att, double &_eigmax_att)
      {
         CFct_bxb_impl < _Int, _Flt >::FctBlockIlu2Degree (_params, _b_size, _n, _n_ini,
                                                           *(_alu_matr.GetIa ()),
                                                           *(_alu_matr.GetJa ()),
                                                           *(_alu_matr.GetJaChar ()),
                                                           *(_alu_matr.GetA ()),
                                                           *(_lu_matr.GetIa ()),
                                                           *(_lu_matr.GetJa ()),
                                                           *(_lu_matr.GetJaChar ()),
                                                           *(_lu_matr.GetA ()), _nmodif,
                                                           _eigmin_att, _eigmax_att);
         vector < _Int > *pia_lu = _lu_matr.GetIa ();
         int nzja_lu = (int) (*pia_lu)[_n];
         _lu_matr.SetNlist (_n);
         _lu_matr.SetNzja (nzja_lu);
         _lu_matr.SetNzjaChar (nzja_lu);
         _lu_matr.SetNza (nzja_lu * _b_size * _b_size * 2);
         _lu_matr.SetIdentityList ();
      };
      /// Perform multiplication of matrices stored by rows, resulting matrix is also stored by rows
      static void MatrixByMatrixMultiply_BxB (int _blksize, int &_icycle, int *_imask,
                                              int *_imask1, int *_indarr, int *_listloc,
                                              _Flt * _fmask, const CMatrix < _Int,
                                              _Flt > &_a, const CMatrix < _Int,
                                              _Flt > &_b, CMatrix < _Int,
                                              _Flt > &_a_times_b);
      /// Perform hmatrix by hmatrix multiplication ( C = AxB )
      static void HMatrixByHMatrixMultiply_BxB (int _blksize, int _nblks, int *_blks,
                                                const CBMatrix < _Int,
                                                _Flt > &_hmatr_a_bxb,
                                                const CBMatrix < _Int,
                                                _Flt > &_hmatr_b_bxb, CBMatrix < _Int,
                                                _Flt > &_hmatr_c_bxb);
   };

   template < typename _Int, typename _Flt, typename _FltVect > class CMvmSlv_BxB_impl  // Class that supports bxb multiplication computations
   {
    public:
      /// Multiply by super sparse small blocks based matrix by rows and add result into prescribed positions
      static void MvmA (char _oper, int _b_size, int _nlist,
                        const vector < _Int > &_list_alu, const vector < _Int > &_ia_alu,
                        const vector < _Int > &_ja_alu, const vector < _Flt > &_a_alu,
                        const _FltVect * _x, _FltVect * _ax);
      /// Multiply by super sparse small blocks based matrix by columns add result into prescribed positions
      static void MvmAT (char _oper, int _b_size, int _nlist,
                         const vector < _Int > &_list_alu, const vector < _Int > &_ia_alu,
                         const vector < _Int > &_ja_alu, const vector < _Flt > &_a_alu,
                         const _FltVect * _x, _FltVect * _ax);
      /// Multiply by Q for super sparse small blocks based matrix Q stored as set of block Housholder transformations with external zero block
      static void MvmQ (int _b_size, int _nlist, int _nrows, const vector < _Int > &_ia_q,
                        const vector < _Int > &_ja_q, const vector < _Flt > &_a_q,
                        const vector < _Flt > &_q_diag, const vector < _Flt > &_tau_arr,
                        const _FltVect * _x_diag, const _FltVect * _x,
                        _FltVect * _qx_diag, _FltVect * _qx, _FltVect * _work);
      /// Multiply by QT for super sparse small blocks based matrix Q stored as set of block Housholder transformations with external zero block
      static void MvmQT (int _b_size, int _nlist, int _nrows,
                         const vector < _Int > &_ia_q, const vector < _Int > &_ja_q,
                         const vector < _Flt > &_a_q, const vector < _Flt > &_q_diag,
                         const vector < _Flt > &_tau_arr, const _FltVect * _x_diag,
                         const _FltVect * _x, _FltVect * _qtx_diag, _FltVect * _qtx,
                         _FltVect * _work);
      /// Solve with L, L is stored by small blocks columns (small blocks diag is inverted)
      static void SolveL (int _n, int _b_size, const vector < _Int > &_ia_l,
                          const vector < _Int > &_ja_l, const vector < _Flt > &_a_l,
                          const _FltVect * _x, _FltVect * _lx, _FltVect * _work);
      /// Solve with L, L is stored by small blocks rows (small blocks diag is inverted)
      static void SolveLRow (int _n, int _b_size, const vector < _Int > &_ia_l,
                             const vector < _Int > &_ja_l, const vector < _Flt > &_a_l,
                             const _FltVect * _x, _FltVect * _lx, _FltVect * _work);
      /// Solve with U, U is stored by small blocks rows (small blocks diag is inverted)
      static void SolveU (int _n, int _b_size, const vector < _Int > &_ia_u,
                          const vector < _Int > &_ja_u, const vector < _Flt > &_a_u,
                          const _FltVect * _x, _FltVect * _ux, _FltVect * _work);
      /// Solve with U, U is stored by small blocks columns (small blocks diag is inverted)
      static void SolveUCol (int _n, int _b_size, const vector < _Int > &_ia_u,
                             const vector < _Int > &_ja_u, const vector < _Flt > &_a_u,
                             const _FltVect * _x, _FltVect * _ux, _FltVect * _work);
   };

   template < typename _Int, typename _Flt, typename _FltVect > class CMvmSlv_BxB       // Class that supports multiplication computations
   {
    public:
      /// Multiply by super sparse matrix by rows and add result into prescribed positions
      static void MvmA (char _oper, int _blksize, const CMatrix < _Int, _Flt > &_a_matr,
                        const _FltVect * _x, _FltVect * _ax)
      {
         CMvmSlv_BxB_impl < _Int, _Flt, _FltVect >::MvmA (_oper, _blksize,
                                                          _a_matr.GetNlist (),
                                                          *(_a_matr.GetList ()),
                                                          *(_a_matr.GetIa ()),
                                                          *(_a_matr.GetJa ()),
                                                          *(_a_matr.GetA ()), _x, _ax);
      };
      /// Multiply by super sparse matrix by rows and add result into prescribed positions
      static void MvmAT (char _oper, int _blksize, const CMatrix < _Int, _Flt > &_a_matr,
                         const _FltVect * _x, _FltVect * _ax)
      {
         CMvmSlv_BxB_impl < _Int, _Flt, _FltVect >::MvmAT (_oper, _blksize,
                                                           _a_matr.GetNlist (),
                                                           *(_a_matr.GetList ()),
                                                           *(_a_matr.GetIa ()),
                                                           *(_a_matr.GetJa ()),
                                                           *(_a_matr.GetA ()), _x, _ax);
      };
      /// Multiply by Q for super sparse small blocks based matrix Q stored as set of block Housholder transformations with external zero block
      static void MvmQ (int _b_size, int _nrows, const CMatrix < _Int, _Flt > &_q_matr,
                        const vector < _Flt > &_q_diag, const vector < _Flt > &_tau_arr,
                        const _FltVect * _x_diag, const _FltVect * _x,
                        _FltVect * _qx_diag, _FltVect * _qx, _FltVect * _work)
      {
         CMvmSlv_BxB_impl < _Int, _Flt, _FltVect >::MvmQ (_b_size, _q_matr.GetNlist (),
                                                          _nrows, *(_q_matr.GetIa ()),
                                                          *(_q_matr.GetJa ()),
                                                          *(_q_matr.GetA ()), _q_diag,
                                                          _tau_arr, _x_diag, _x, _qx_diag,
                                                          _qx, _work);
      };
      /// Multiply by QT for super sparse small blocks based matrix Q stored as set of block Housholder transformations with external zero block
      static void MvmQT (int _b_size, int _nrows, const CMatrix < _Int, _Flt > &_q_matr,
                         const vector < _Flt > &_q_diag, const vector < _Flt > &_tau_arr,
                         const _FltVect * _x_diag, const _FltVect * _x,
                         _FltVect * _qtx_diag, _FltVect * _qtx, _FltVect * _work)
      {
         CMvmSlv_BxB_impl < _Int, _Flt, _FltVect >::MvmQT (_b_size, _q_matr.GetNlist (),
                                                           _nrows, *(_q_matr.GetIa ()),
                                                           *(_q_matr.GetJa ()),
                                                           *(_q_matr.GetA ()), _q_diag,
                                                           _tau_arr, _x_diag, _x,
                                                           _qtx_diag, _qtx, _work);
      };
      /// Solve with L, L is stored by columns (diag is inverted)
      static void SolveL (int _blksize, const CMatrix < _Int, _Flt > &_l_matr,
                          const _FltVect * _x, _FltVect * _lx, _FltVect * _work)
      {
         CMvmSlv_BxB_impl < _Int, _Flt, _FltVect >::SolveL (_blksize, _l_matr.GetN (),
                                                            *(_l_matr.GetIa ()),
                                                            *(_l_matr.GetJa ()),
                                                            *(_l_matr.GetA ()), _x, _lx,
                                                            _work);
      };
      /// Solve with U, U is stored by rows (diag is inverted)
      static void SolveUCol (int _blksize, const CMatrix < _Int, _Flt > &_u_matr,
                             const _FltVect * _x, _FltVect * _ux, _FltVect * _work)
      {
         CMvmSlv_BxB_impl < _Int, _Flt, _FltVect >::SolveUCol (_blksize, _u_matr.GetN (),
                                                               *(_u_matr.GetIa ()),
                                                               *(_u_matr.GetJa ()),
                                                               *(_u_matr.GetA ()), _x,
                                                               _ux, _work);
      };
      /// Solve with U, U is stored by rows (diag is inverted)
      static void SolveU (int _blksize, const CMatrix < _Int, _Flt > &_u_matr,
                          const _FltVect * _x, _FltVect * _ux, _FltVect * _work)
      {
         CMvmSlv_BxB_impl < _Int, _Flt, _FltVect >::SolveU (_blksize, _u_matr.GetN (),
                                                            *(_u_matr.GetIa ()),
                                                            *(_u_matr.GetJa ()),
                                                            *(_u_matr.GetA ()), _x, _ux,
                                                            _work);
      };
   };

   template < typename _Int, typename _Flt > class CFctBxBThreads       // Class that supports factorization computations with threads
   {
    public:
// Functions that support hmatrix structure as a whole
      /// Perform ILU2 block bxb factorization of the hmatrix with dinamic ordering and future diagonal modifications
      static void Ilu2HMatrix_BxB (bool _b_blk_wells, SParams & _params, int _blksize,
                                   int _nblks, int _nblks1, int _nblks2, long long *_blks,
                                   long long *_nzord_blks, CTree & _tree0, CTree & _tree1,
                                   CBMatrix < _Int, _Flt > &_a_matr, CBMatrix < _Int,
                                   _Flt > &_l_matr, CBMatrix < _Int, _Flt > &_u_matr,
                                   CTree & _tree2_new, int &_nblks_new,
                                   vector < long long >&_blks_new,
                                   CVectorData < int >&_ordernew, int &_nmodif_scl,
                                   double &_sclmin_att, double &_sclmax_att,
                                   int &_ndiasplit, int &_nmodif, double &_eigmin_att,
                                   double &_eigmax_att);
// Implementation functions
      /// Compute explicit scaling
      static void ComputeExplicitScaling_BxB (int _sctype, int _nitersc, double _sclmin,
                                              int _blksize, int _nblks, long long *_blks,
                                              CBMatrix < _Int, _Flt > &_a_matr,
                                              CVectorData < _Flt > &_scl_L,
                                              CVectorData < _Flt > &_scl_U,
                                              CVectorData < _Flt > &_scl_L_inv,
                                              CVectorData < _Flt > &_scl_U_inv,
                                              int &_n_modif, double &_sclmin_att,
                                              double &_sclmax_att);
      /// Perform explicit scaling
      static void HMatrixScale_BxB (int _blksize, int _nblks, long long *_blks,
                                    CBMatrix < _Int, _Flt > &_a_matr, _Flt * _sclL,
                                    _Flt * _sclU);
      /// Split LU and combine pairs
      static void SplitLUPair_BxB (bool _b_is_char, int _blksize, int _nblks,
                                   long long *_blks, CBMatrix < _Int, _Flt > &_a_matr,
                                   CBMatrix < _Int, _Flt > &_a_pair);
      /// Split pair into L and U
      static void SplitPair_BxB (bool _b_is_char, int _blksize, int _nblks,
                                 long long *_blks, CBMatrix < _Int, _Flt > &_a_pair,
                                 CBMatrix < _Int, _Flt > &_a_L, CBMatrix < _Int,
                                 _Flt > &_a_U);
      /// Perform ILU2 point factorization of the hmatrix with dynamic ordering and future diagonal modifications
      static void Ilu2BlockIlu2_BxB (bool _b_blk_wells, ofstream * _pfout, int _index_out,
                                     SParams * _pparams, int _blksize, int _nblks,
                                     int _nblks1, int _nblks2, long long *_blks,
                                     long long *_nzord_blks, CTree & _tree0,
                                     CTree & _tree1, CBMatrix < _Int, _Flt > &_alu_pair,
                                     CTree & _tree2_new, int &_nblks_new,
                                     vector < long long >&_blks_new,
                                     CVectorData < int >&_ordernew, CBMatrix < _Int,
                                     _Flt > &_lu_pair, int &_ndiasplit, int &_nmodif,
                                     double &_eigmin_att, double &_eigmax_att);
      /// Perform ICH2 block bxb factorization of the hmatrix with dinamic ordering and future diagonal modifications
      static void Ich2HMatrix_BxB (SParams & _params, int _blksize, int _nhblks,
                                   int *_hblk2blks, int _nblks, long long *_blks,
                                   CTree & _tree, CTree * _tree_arr, CBMatrix < _Int,
                                   _Flt > &_au_matr, CBMatrix < _Int, _Flt > &_u_matr,
                                   int &_nmodif, double &_eigmin_att,
                                   double &_eigmax_att);
      /// Rescale U
      static void RescaleU_BxB (int _blksize, int _nblks, long long *_blks,
                                CBMatrix < _Int, _Flt > &_U_matr, _Flt * _sclU,
                                _Flt * _inv_sclU);
      /// Set diagonal by zeroes
      static void SetZeroDiag_BxB (int _blksize, int _nblks, long long *_blks,
                                   CBMatrix < _Int, _Flt > &_A_matr);
      /// Perform filtering of the pairs block sparsity according to the tree with diagonal modifications
      static void PairsFilterTreeModif_BxB (CTree & _tree, double _theta, int _blksize,
                                            int _nblks, long long *_blks, CBMatrix < _Int,
                                            _Flt > &_AU_matr);
      /// Compute reordered by columns rectangular hmatrix
      static void ReorderHMatrixColumnsPairs_BxB (int _blksize, int _nblksR,
                                                  long long *_blksR, int _nblksC_ini,
                                                  long long *_blksC_ini, CBMatrix < _Int,
                                                  _Flt > &_hmatr_ini, int *_order,
                                                  int _nblksC_fin, long long *_blksC_fin,
                                                  CBMatrix < _Int, _Flt > &_hmatr_fin);
   };

 template < typename _Int, typename _Flt, typename _FltVect > class CMvmParBxBThreads:public CMvmParThreads < _Int, _Flt,
      _FltVect
      >
      // Class that supports combined parallel multiplication computations (MPI+thr version)
   {
      int blksize;
// Constructors and destructor
    private:
      /// Copy constructor
      CMvmParBxBThreads (const CMvmParBxBThreads < _Int, _Flt, _FltVect > &_aa)
      {
      };
      /// Copy operator
      CMvmParBxBThreads < _Int, _Flt,
         _FltVect > &operator= (const CMvmParBxBThreads < _Int, _Flt, _FltVect > &_aa)
      {
         return *this;
      };
    public:
      /// Epmty constructor
      CMvmParBxBThreads () {
         blksize = 1;
      };
      /// Destructor
      ~CMvmParBxBThreads () {
      };
// External functions
      /// Init control data
      void InitControl_BxB (void *_pcomm, int _blksize, int _nhblks, int *_hblk2cpu,
                            int *_hblk2blks, int *_blk2hblks, long long *_hblks,
                            int _nblks, long long *_blks)
      {
         this->InitControl (_pcomm, _nhblks, _hblk2cpu, _hblk2blks, _blk2hblks, _hblks,
                            _nblks, _blks);
         blksize = _blksize;
      };
      /// Init MvmA data
      void InitMvmA_BxB (CBMatrix < _Int, _Flt > *_hmatr_arr)
      {
         this->InitMvmA (_hmatr_arr);
         int isize = this->x_send.GetLength ();
         this->x_send.resize (isize * this->blksize);
         isize = this->x_recv.GetLength ();
         this->x_recv.resize (isize * this->blksize);
         isize = this->x_temp.GetLength ();
         this->x_temp.resize (isize * this->blksize);
      };
      /// Perform MvmA computations
      void MvmA_BxB (const _FltVect * _x, _FltVect * _ax);
      /// Clean MvmA structure
      void Clean_BxB ()
      {
         this->Clean ();
         this->blksize = 1;
      };
   };

 template < typename _Int, typename _Flt, typename _FltVect > class CSlvParBxBThreads:public CSlvParThreads < _Int, _Flt,
      _FltVect
      >
      // Class that supports parallel solve computations (including L, U and ordering) (MPI+thr version)
   {
      int blksize;
// Constructors and destructor
    private:
      /// Copy constructor
      CSlvParBxBThreads (const CSlvParBxBThreads < _Int, _Flt, _FltVect > &_aa)
      {
      };
      /// Copy operator
      CSlvParBxBThreads < _Int, _Flt,
         _FltVect > &operator= (const CSlvParBxBThreads < _Int, _Flt, _FltVect > &_aa)
      {
         return *this;
      };
    public:
      /// Epmty constructor
      CSlvParBxBThreads () {
         blksize = 1;
      };
      /// Destructor
      ~CSlvParBxBThreads () {
      };
      /// Get ni_cpu
      virtual int GetNiCpu ()
      {
         int ni_temp = this->GetNiCpuBase ();
         return ni_temp * blksize;
      };
// External functions
      /// Init control data
      void InitControl_BxB (void *_pcomm, int _blksize, int _nhblks, long long *_hblks,
                            int *_hblk2cpu, int *_hblk2blks, int *_blk2hblks,
                            long long *_blks)
      {
         this->InitControl (_pcomm, _nhblks, _hblks, _hblk2cpu, _hblk2blks, _blk2hblks,
                            _blks);
         blksize = _blksize;
      };
      /// Init SolveLU data
      void InitSolveLU_BxB (int _blksize, long long *_hblks_ext,
                            vector < int >*_listpairs_ext, int *_nblks_ext_arr,
                            vector < int >*_blksnum_ext_arr,
                            vector < long long >*_blks_ext_arr, CTree * _tree_arr,
                            int *_nblks_ilu2_arr, int *_nblks1_ilu2_arr,
                            int *_nblks2_ilu2_arr, vector < long long >*_blks_ilu2_arr,
                            CVectorData < int >*_orderLU, CBMatrix < _Int, _Flt > *_matrL,
                            CBMatrix < _Int, _Flt > *_matrU)
      {
         this->InitSolveLU (_hblks_ext, _listpairs_ext, _nblks_ext_arr, _blksnum_ext_arr,
                            _blks_ext_arr, _tree_arr, _nblks_ilu2_arr, _nblks1_ilu2_arr,
                            _nblks2_ilu2_arr, _blks_ilu2_arr, _orderLU, _matrL, _matrU);
         int isize = this->x_send.GetLength ();
         this->x_send.resize (isize * this->blksize);
         isize = this->x_recv.GetLength ();
         this->x_recv.resize (isize * this->blksize);
         isize = this->x_temp.GetLength ();
         this->x_temp.resize (isize * this->blksize);
      };
      /// Perform SolveLU computations
      void SolveLU_BxB (const _FltVect * _x, _FltVect * _px);
      /// In-place perform SolveU computations
      void SolveU_BxB (int _ihblk, _FltVect * _px);
      /// In-place perform SolveL computations
      void SolveL_BxB (int _ihblk, _FltVect * _px);
      /// Clean Slv structure
      void Clean_BxB ()
      {
         this->Clean ();
      };
   };

 template < typename _Int, typename _Flt, typename _FltVect > class CK3D_SolverBxBThreads:public CK3D_SolverThreads < _Int, _Flt,
      _FltVect
      >
      // Class that supports combined parallel small block solver computations (MPI+thr)
   {
      int blksize;
      CMvmParBxBThreads < _Int, _Flt, _FltVect > mvm_bxb;       ///< Support MvmA
      CSlvParBxBThreads < _Int, float, _FltVect > slv_float_bxb;        ///< Support float  Slv
      CSlvParBxBThreads < _Int, double, _FltVect > slv_double_bxb;      ///< Support double Slv
// Constructors and destructor
    private:
      /// Copy CK3D_SolverBxBThreads
      CK3D_SolverBxBThreads < _Int, _Flt, _FltVect > (const CK3D_SolverBxBThreads < _Int,
                                                      _Flt, _FltVect > &_aa) {
      }
      /// Copy operator
      CK3D_SolverBxBThreads < _Int, _Flt,
         _FltVect > &operator= (const CK3D_SolverBxBThreads < _Int, _Flt, _FltVect > &_aa)
      {
         return *this;
      }
    public:
      /// Epmty constructor
      CK3D_SolverBxBThreads () {
         blksize = 1;
      }
      /// Destructor
      ~CK3D_SolverBxBThreads () {
      }
// Get/set functions
      /// Set blksize
      void SetBlksize (int _blksize)
      {
         blksize = _blksize;
      }
      /// Get blksize
      int GetBlksize ()
      {
         return blksize;
      }
      /// Get ni
      virtual int GetNi ()
      {
         int ni_local = 0;
         if (this->params.prec_float == 1) {
            ni_local = slv_float_bxb.GetNiCpu ();
         } else {
            ni_local = slv_double_bxb.GetNiCpu ();
         }
         return ni_local;
      };
// Aggregated interfaces
      /// Prepare matrix structures
      void PrepareMatrix_BxB (void *_pcomm, SParams * _params, SStatData * _stats,
                              int _blksize, int _nhblks, int *_hblk2cpu, int *_hblk2blks,
                              long long *_hblks, int _nblks, long long *_blks,
                              int _format_type, _Int * _ia, _Int * _ja, _Flt * _a)
      {
         int myid = CMPIDataExchange::GetMyid (_pcomm);
         this->blksize = _blksize;
         if (_params->i_decomp_type == 0) {
            this->PrepareMatrixBase_BxB (_pcomm, _blksize, _nhblks, _hblk2cpu, _hblk2blks,
                                         _hblks, _nblks, _blks, _format_type, _ia, _ja,
                                         _a);
         } else if (_params->i_decomp_type == 1) {
            this->PrepareMatrixThreads_BxB (_pcomm, _params, _stats, _blksize, _nhblks,
                                            _hblk2cpu, _hblk2blks, _hblks, _nblks, _blks,
                                            _format_type, _ia, _ja, _a);
         } else if (_params->i_decomp_type == 2) {
            this->PrepareMatrixWells_BxB (_pcomm, _params, _stats, _blksize, _nhblks,
                                          _hblk2cpu, _hblk2blks, _hblks, _nblks, _blks,
                                          NULL, NULL, _format_type, _ia, _ja, _a);
         }
         if (_params->msglev >= 2 && myid == 0) {
            cout << "   Nhblks_ini = " << _nhblks << " Nblks_ini = " << _nblks;
            cout << " Nhblks = " << this->nhblks << " Nblks = " << this->
               nblks << " Nsup_well = " << _stats->n_well << endl;
         }
         if (_params->msglev >= 1 && myid == 0 && _params->pfout != NULL) {
            *(_params->
              pfout) << "   Nhblks_ini = " << _nhblks << " Nblks_ini = " << _nblks;
            *(_params->pfout) << " Nhblks = " << this->nhblks << " Nblks = " << this->
               nblks << " Nsup_well = " << _stats->n_well << endl;
         }
      }
      /// Perform fct and free matrix if necessary
      void ComputePreconditioner_BxB (SParams * _params, SStatData * _stats,
                                      bool _b_store_matrix)
      {
         this->ComputeBILU2_BxB (_params, _stats);
         if (!_b_store_matrix)
            this->CleanMvmA_BxB ();
      }
      /// Prepare solver structures including performing parallel fct
      void PrepareSolver_BxB (void *_pcomm, SParams * _params, SStatData * _stats,
                              int _blksize, int _nhblks, int *_hblk2cpu, int *_hblk2blks,
                              long long *_hblks, int _nblks, long long *_blks,
                              bool _b_store_matrix, int _format_type, _Int * _ia,
                              _Int * _ja, _Flt * _a)
      {
         this->PrepareMatrix_BxB (_pcomm, _params, _stats, _blksize, _nhblks, _hblk2cpu,
                                  _hblk2blks, _hblks, _nblks, _blks, _format_type, _ia,
                                  _ja, _a);
         this->ComputePreconditioner_BxB (_params, _stats, _b_store_matrix);
      }
      /// Perform iterations of the iterative scheme
      void SolveIter_BxB (SParams * _params, SStatData * _stats, _FltVect * _rhs,
                          _FltVect * _sol);
      /// Perform iterations of the iterative scheme (implementation)
      void SolveIter_BxB (SParams * _params, SStatData * _stats, void *_str_mvmA,
                          typename TMvmSlvFunc < _FltVect >::mvmA_f _mvm_f,
                          void *_str_slvLU,
                          typename TMvmSlvFunc < _FltVect >::slvLU_f _slv_f,
                          _FltVect * _rhs, _FltVect * _sol);
      /// Perform iterations of the iterative scheme (implementation) (no vector data transformation)
      void SolveIter_BxB_NoTransform (SParams * _params, SStatData * _stats,
                                      void *_str_mvmA,
                                      typename TMvmSlvFunc < _FltVect >::mvmA_f _mvm_f,
                                      void *_str_slvLU,
                                      typename TMvmSlvFunc < _FltVect >::slvLU_f _slv_f,
                                      _FltVect * _rhs, _FltVect * _sol);
      /// Transform vectors forward
      void TransformVectorsForward_BxB (_FltVect * _rhs_ini, _FltVect * _sol_ini,
                                        vector < _FltVect > &_rhs_new,
                                        vector < _FltVect > &_sol_new);
      /// Transform vectors backward
      void TransformVectorsBackward_BxB (vector < _FltVect > &_rhs_new,
                                         vector < _FltVect > &_sol_new,
                                         _FltVect * _rhs_fin, _FltVect * _sol_fin);
// External functions
      /// Store control data only
      void StoreControlData_BxB (void *_pcomm, int _blksize, int _nhblks, int *_hblk2cpu,
                                 int *_hblk2blks, long long *_hblks, int _nblks,
                                 long long *_blks)
      {
         this->StoreControlData (_pcomm, _nhblks, _hblk2cpu, _hblk2blks, _hblks, _nblks,
                                 _blks);
         int *phblk2cpu = this->GetHBlk2cpu ();
         int *phblk2blks = this->GetHBlk2blks ();
         int *pblk2hblks = this->GetBlk2hblks ();
         long long *phblks = this->GetHBlks ();
         long long *pblks = this->GetBlks ();
         this->mvm_bxb.InitControl_BxB (_pcomm, _blksize, _nhblks, phblk2cpu, phblk2blks,
                                        pblk2hblks, phblks, _nblks, pblks);
         this->slv_float_bxb.InitControl_BxB (_pcomm, _blksize, _nhblks, phblks,
                                              phblk2cpu, phblk2blks, pblk2hblks, pblks);
         this->slv_double_bxb.InitControl_BxB (_pcomm, _blksize, _nhblks, phblks,
                                               phblk2cpu, phblk2blks, pblk2hblks, pblks);
         this->blksize = _blksize;
      };
      /// Prepare matrix data only
      void PrepareMatrixBase_BxB (void *_pcomm, int _blksize, int _nhblks, int *_hblk2cpu,
                                  int *_hblk2blks, long long *_hblks, int _nblks,
                                  long long *_blks, int _format_type, _Int * _ia,
                                  _Int * _ja, _Flt * _a);
      /// Prepare matrix data only including preliminary repartitioning of the matrix
      void PrepareMatrixThreads_BxB (void *_pcomm, SParams * _params, SStatData * _stats,
                                     int _blksize, int _nhblks, int *_hblk2cpu,
                                     int *_hblk2blks_ini, long long *_hblks,
                                     int _nblks_ini, long long *_blks_ini,
                                     int _format_type, _Int * _ia, _Int * _ja, _Flt * _a);
      /// Prepare matrix data only including preliminary repartitioning of the matrix taking into account wells data
      void PrepareMatrixWells_BxB (void *_pcomm, SParams * _params, SStatData * _stats,
                                   int _blksize, int _nhblks_ini, int *_hblk2cpu_ini,
                                   int *_hblk2blks_ini, long long *_hblks_ini,
                                   int _nblks_ini, long long *_blks_ini,
                                   vector < int >*_p_ia_wells_ext,
                                   vector < int >*_p_ja_wells_ext, int _format_type,
                                   _Int * _ia, _Int * _ja, _Flt * _a);
      /// Clean partitioning data
      void CleanPartitioning_BxB ()
      {
         this->CleanPartitioning ();
         this->blksize = 1;
      }
      /// Clean MvmA structures
      void CleanMvmA_BxB ()
      {
         this->CleanMvmA ();
         this->mvm_bxb.Clean_BxB ();
      }
      /// Clean Slv structures
      void CleanSlv_BxB ()
      {
         this->CleanSlv ();
         this->slv_float_bxb.Clean_BxB ();
         this->slv_double_bxb.Clean_BxB ();
      }
      /// Clean all structures
      void Clean_BxB ()
      {
         this->Clean ();
         this->blksize = 1;
         this->CleanMvmA_BxB ();
         this->CleanSlv_BxB ();
      }
// Internal functions
    private:
      /// Store matrix data
      static void StoreMatrix_BxB (void *_pcomm, int _blksize, int _nhblks,
                                   int *_hblk2cpu, int *_hblk2blks, long long *_hblks,
                                   int _nblks, long long *_blks, int _format_type,
                                   _Int * _ia, _Int * _ja, _Flt * _a, CBMatrix < _Int,
                                   _Flt > *_hmatr_arr)
      {
         if (_format_type == 0) {
            CK3D_SolverBxBThreads < _Int, _Flt, _FltVect >::StoreMatrix_BxB_Point (_pcomm,
                                                                                   _blksize,
                                                                                   _nhblks,
                                                                                   _hblk2cpu,
                                                                                   _hblk2blks,
                                                                                   _hblks,
                                                                                   _nblks,
                                                                                   _blks,
                                                                                   _ia,
                                                                                   _ja,
                                                                                   _a,
                                                                                   _hmatr_arr);
         } else {
            CK3D_SolverBxBThreads < _Int, _Flt, _FltVect >::StoreMatrix_BxB_Block (_pcomm,
                                                                                   _blksize,
                                                                                   _nhblks,
                                                                                   _hblk2cpu,
                                                                                   _hblk2blks,
                                                                                   _hblks,
                                                                                   _nblks,
                                                                                   _blks,
                                                                                   _ia,
                                                                                   _ja,
                                                                                   _a,
                                                                                   _hmatr_arr);
         }
      }
      /// Store matrix data
      static void StoreMatrix_BxB_Point (void *_pcomm, int _blksize, int _nhblks,
                                         int *_hblk2cpu, int *_hblk2blks,
                                         long long *_hblks, int _nblks, long long *_blks,
                                         _Int * _ia, _Int * _ja, _Flt * _a,
                                         CBMatrix < _Int, _Flt > *_hmatr_arr);
      /// Store matrix data
      static void StoreMatrix_BxB_Block (void *_pcomm, int _blksize, int _nhblks,
                                         int *_hblk2cpu, int *_hblk2blks,
                                         long long *_hblks, int _nblks, long long *_blks,
                                         _Int * _ia, _Int * _ja, _Flt * _a,
                                         CBMatrix < _Int, _Flt > *_hmatr_arr);
      /// Prepare fast transform data
      void PrepareFastTransform_BxB (CBMatrix < _Int, _Flt > *_hmatr_ini_arr)
      {
      };
      /// Use fast transform data
      void UseFastTransform_BxB (CBMatrix < _Int, _Flt > *_hmatr_ini_arr)
      {
      };
      /// Peform BILU2 computations
      void ComputeBILU2_BxB (SParams * _params, SStatData * _stats);
      /// Peform BILU2 computations
      void ComputeBILU2_BxB_impl (SParams * _params, SStatData * _stats);
      /// Reorder rhs and solution from/to initial ordering
      void ReorderVectorDataIni_BxB (char _dir, _FltVect * _rhs, _FltVect * _sol);
      /// Reorder rhs and solution from/to wells ordering
      void ReorderVectorDataWells_BxB (char _dir, _FltVect * _rhs, _FltVect * _sol,
                                       CVectorData < _FltVect > &_rhs_ord,
                                       CVectorData < _FltVect > &_sol_ord);
      /// Perform MvmA parallel computations
      void MvmA_BxB (const _FltVect * _x, _FltVect * _ax)
      {
         mvm_bxb.MvmA_BxB (_x, _ax);
      };
      /// Perform SlvLU parallel computations
      void SlvLU_BxB (const _FltVect * _x, _FltVect * _px);
      /// Perform mvm by part of block scaling
      void MvmBScl (char _type, char _transp, const _FltVect * _x, _FltVect * _px);
    public:
      /// Multiply by A
      static void MvmA_BxB_static (void *_mvmA, const _FltVect * _x, _FltVect * _ax)
      {
         CK3D_SolverBxBThreads < _Int, _Flt, _FltVect > *p_mvmA =
            (CK3D_SolverBxBThreads < _Int, _Flt, _FltVect > *)_mvmA;
         p_mvmA->MvmA_BxB (_x, _ax);
      }
      /// Solve with LU for the block
      static void SlvLU_BxB_static (void *_slvLU, const _FltVect * _x, _FltVect * _px)
      {
         CK3D_SolverBxBThreads < _Int, _Flt, _FltVect > *p_slvLU =
            (CK3D_SolverBxBThreads < _Int, _Flt, _FltVect > *)_slvLU;
         p_slvLU->SlvLU_BxB (_x, _px);
      }
   };

}                               // namespace k3d

#endif
