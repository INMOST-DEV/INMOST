//------------------------------------------------------------------------------------------------
// File: k3d_slv.hxx
//------------------------------------------------------------------------------------------------

#ifndef __K3D_SLV_HXX
#define __K3D_SLV_HXX

#include "k3d_base.hxx"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

namespace k3d
{

   template < typename _Int, typename _Flt, typename _FltVect > class CMvmSlv_impl      // Class that supports multiplication computations
   {
    public:
      /// Multiply by super sparse matrix by rows and add result into prescribed positions
      static void MvmA (char _oper, int _nlist, const vector < _Int > &_list_alu,
                        const vector < _Int > &_ia_alu, const vector < _Int > &_ja_alu,
                        const vector < _Flt > &_a_alu, const _FltVect * _x,
                        _FltVect * _ax);
      /// Multiply by block diagonal by block rows, each diagonal block is blksize x blksize dense matrix
      static void MvmDiagBlksize (int _blksize, int _nlist, const vector < _Flt > &_a_alu,
                                  const _FltVect * _x, _FltVect * _ax);
      /// Multiply by block diagonal by block rows, each diagonal block is blksize x blksize dense matrix
      static void MvmDiagTBlksize (int _blksize, int _nlist,
                                   const vector < _Flt > &_a_alu, const _FltVect * _x,
                                   _FltVect * _ax);
      /// Multiply by super sparse matrix by columns add result into prescribed positions
      static void MvmAT (char _oper, int _nlist, const vector < _Int > &_list_alu,
                         const vector < _Int > &_ia_alu, const vector < _Int > &_ja_alu,
                         const vector < _Flt > &_a_alu, const _FltVect * _x,
                         _FltVect * _ax);
      /// Solve with L, L is stored by columns (diag is inverted)
      static void SolveL (int _n, const vector < _Int > &_ia_l,
                          const vector < _Int > &_ja_l, const vector < _Flt > &_a_l,
                          const _FltVect * _x, _FltVect * _lx);
      /// Solve with U, U is stored by rows (diag is inverted)
      static void SolveU (int _n, const vector < _Int > &_ia_u,
                          const vector < _Int > &_ja_u, const vector < _Flt > &_a_u,
                          const _FltVect * _x, _FltVect * _ux);
   };

   template < typename _Int, typename _Flt, typename _FltVect > class CMvmSlv   // Class that supports multiplication computations
   {
    public:
      /// Multiply by super sparse matrix by rows and add result into prescribed positions
      static void MvmA (char _oper, const CMatrix < _Int, _Flt > &_a_matr,
                        const _FltVect * _x, _FltVect * _ax)
      {
         CMvmSlv_impl < _Int, _Flt, _FltVect >::MvmA (_oper, _a_matr.GetNlist (),
                                                      *(_a_matr.GetList ()),
                                                      *(_a_matr.GetIa ()),
                                                      *(_a_matr.GetJa ()),
                                                      *(_a_matr.GetA ()), _x, _ax);
      };
      /// Multiply by super sparse matrix by rows and add result into prescribed positions
      static void MvmAT (char _oper, const CMatrix < _Int, _Flt > &_a_matr,
                         const _FltVect * _x, _FltVect * _ax)
      {
         CMvmSlv_impl < _Int, _Flt, _FltVect >::MvmAT (_oper, _a_matr.GetNlist (),
                                                       *(_a_matr.GetList ()),
                                                       *(_a_matr.GetIa ()),
                                                       *(_a_matr.GetJa ()),
                                                       *(_a_matr.GetA ()), _x, _ax);
      };
      /// Solve with L, L is stored by columns (diag is inverted)
      static void SolveL (const CMatrix < _Int, _Flt > &_l_matr, const _FltVect * _x,
                          _FltVect * _lx)
      {
         CMvmSlv_impl < _Int, _Flt, _FltVect >::SolveL (_l_matr.GetN (),
                                                        *(_l_matr.GetIa ()),
                                                        *(_l_matr.GetJa ()),
                                                        *(_l_matr.GetA ()), _x, _lx);
      };
      /// Solve with U, U is stored by rows (diag is inverted)
      static void SolveU (const CMatrix < _Int, _Flt > &_u_matr, const _FltVect * _x,
                          _FltVect * _ux)
      {
         CMvmSlv_impl < _Int, _Flt, _FltVect >::SolveU (_u_matr.GetN (),
                                                        *(_u_matr.GetIa ()),
                                                        *(_u_matr.GetJa ()),
                                                        *(_u_matr.GetA ()), _x, _ux);
      };
   };

   template < typename _Int, typename _Flt > class CFct // Class that supports factorization computations
   {
    public:
// Functions that support matrix structure as a whole
      /// Compute optimal ordering for the matrix
      static void ComputeOptimalOrder (const CMatrix < _Int, _Flt > &_a_matr,
                                       int _ordtype, vector < int >&_order)
      {
         CFct_impl < _Int, _Flt >::ComputeOptimalOrder (_ordtype, _a_matr.GetNlist (),
                                                        *(_a_matr.GetIa ()),
                                                        *(_a_matr.GetJa ()), _order);
      };
      /// Compute ordered matrix
      static void ReorderMatrix (CMatrix < _Int, _Flt > &_a_matr,
                                 const vector < int >&_order, CMatrix < _Int,
                                 _Flt > &_a_ord)
      {
         _a_ord.SetN (_a_matr.GetNlist ());
         _a_ord.SetNzja (_a_matr.GetNzja ());
         _a_ord.SetNza (_a_matr.GetNza ());
         CFct_impl < _Int, _Flt >::ReorderMatrix (_a_matr.GetNlist (), _order,
                                                  *(_a_matr.GetIa ()),
                                                  *(_a_matr.GetJa ()), *(_a_matr.GetA ()),
                                                  *(_a_ord.GetIa ()), *(_a_ord.GetJa ()),
                                                  *(_a_ord.GetA ()));
         _a_ord.SetIdentityList ();
      };
      /// Perform ILU2 point factorization of the block with future diagonal modification
      static void Ich2Matrix (CMatrix < _Int, _Flt > &_a_matr, SParams & _params,
                              CMatrix < _Int, _Flt > &_u_matr, double &_sclmin_att,
                              double &_sclmax_att, int &_nmodif, double &_eigmin_att,
                              double &_eigmax_att)
      {
         _u_matr.SetN (_a_matr.GetNlist ());
         CFct_impl < _Int, _Flt >::Ich2Block (_params.sclmin, _params.fcttype,
                                              _params.pivmin, _params.tau1, _params.tau2,
                                              _params.theta, _a_matr.GetNlist (),
                                              *(_a_matr.GetIa ()), *(_a_matr.GetJa ()),
                                              *(_a_matr.GetA ()), *(_u_matr.GetIa ()),
                                              *(_u_matr.GetJa ()), *(_u_matr.GetA ()),
                                              _sclmin_att, _sclmax_att, _nmodif,
                                              _eigmin_att, _eigmax_att);
         int nlist_a = _a_matr.GetNlist ();
         _Int *pia_umatr = _u_matr.GetIaArr ();
         int nzja_u = (int) pia_umatr[nlist_a];
         _u_matr.SetNzja (nzja_u);
         _u_matr.SetNza (nzja_u);
         _u_matr.SetIdentityList ();
      };
      /// Perform ILU2 point factorization of the block with future diagonal modification
      static void Ilu2Matrix (CMatrix < _Int, _Flt > &_a_matr, SParams & _params,
                              CMatrix < _Int, _Flt > &_l_matr, CMatrix < _Int,
                              _Flt > &_u_matr, double &_sclmin_att, double &_sclmax_att,
                              int &_nmodif, double &_eigmin_att, double &_eigmax_att)
      {
         _l_matr.SetN (_a_matr.GetNlist ());
         _u_matr.SetN (_a_matr.GetNlist ());
         CFct_impl < _Int, _Flt >::Ilu2Block (_params.sctype, _params.nitersc,
                                              _params.fcttype, _params.pivmin,
                                              _params.tau1, _params.tau2, _params.theta,
                                              _a_matr.GetNlist (), *(_a_matr.GetIa ()),
                                              *(_a_matr.GetJa ()), *(_a_matr.GetA ()),
                                              *(_l_matr.GetIa ()), *(_l_matr.GetJa ()),
                                              *(_l_matr.GetA ()), *(_u_matr.GetIa ()),
                                              *(_u_matr.GetJa ()), *(_u_matr.GetA ()),
                                              _sclmin_att, _sclmax_att, _nmodif,
                                              _eigmin_att, _eigmax_att);
         int nlist_a = _a_matr.GetNlist ();
         _Int *pia_lmatr = _l_matr.GetIaArr ();
         _Int *pia_umatr = _u_matr.GetIaArr ();
         int nzja_l = (int) pia_lmatr[nlist_a];
         int nzja_u = (int) pia_umatr[nlist_a];
         _l_matr.SetNzja (nzja_l);
         _l_matr.SetNza (nzja_l);
         _u_matr.SetNzja (nzja_u);
         _u_matr.SetNza (nzja_u);
         _l_matr.SetIdentityList ();
         _u_matr.SetIdentityList ();
      };
      /// Balance diagonal
      static void BalanceDiag (CMatrix < _Int, _Flt > &_a_matr, CMatrix < _Int,
                               _Flt > &_l_matr, CMatrix < _Int, _Flt > &_u_matr,
                               double &_diacorr_min, double &_diacorr_max)
      {
         CFct_impl < _Int, _Flt >::BalanceDiag (_a_matr.GetNlist (), *(_a_matr.GetIa ()),
                                                *(_a_matr.GetJa ()), *(_a_matr.GetA ()),
                                                *(_l_matr.GetIa ()), *(_l_matr.GetJa ()),
                                                *(_l_matr.GetA ()), *(_u_matr.GetIa ()),
                                                *(_u_matr.GetJa ()), *(_u_matr.GetA ()),
                                                _diacorr_min, _diacorr_max);
      };
      /// Print matrix data
      static void PrintMatrix (ofstream & _fout, CMatrix < _Int, _Flt > &_a_matr);
      /// Print matrix data over thresh
      static void PrintMatrixThresh (ofstream & _fout, CMatrix < _Int, _Flt > &_a_matr,
                                     double _thresh);
      /// Print diagonal part of matrix data
      static void PrintMatrixDiag (ofstream & _fout, CMatrix < _Int, _Flt > &_a_matr);
      /// Print diagonal part of matrix data in pairs format
      static void PrintMatrixDiagPair (ofstream & _fout, CMatrix < _Int, _Flt > &_a_matr);
   };

   template < typename _Int, typename _Flt > class CBMatrix     /// Class that supports parallel submatrix data
   {
      int nzblk;                ///< The number of blocks
      CMatrix < int, float >hmatr_str;  ///< Block sparsity of a submatrix
      vector < CMatrix < _Int, _Flt > >asub_arr;        ///< Data of each block
// Get template names
      typedef _Int THInt_type;  ///< Reference to int type
      typedef _Flt THFloat_type;        ///< Reference to float type
// Constructors and destructor
    public:
      /// Epmty constructor
      CBMatrix < _Int, _Flt > () {
         nzblk = 0;
         asub_arr.resize (1);
      };
      /// Copy constructor
      CBMatrix < _Int, _Flt > (const CBMatrix < _Int, _Flt > &_aa);
      /// Equality operator
      CBMatrix < _Int, _Flt > &operator= (const CBMatrix < _Int, _Flt > &_aa);
      /// Add replace operator
      CBMatrix < _Int, _Flt > &operator+= (const CBMatrix < _Int, _Flt > &_aa);
      /// Add replace operator for pairs
      CBMatrix < _Int, _Flt > &operator%= (const CBMatrix < _Int, _Flt > &_aa);
      /// Init by data
      CBMatrix < _Int, _Flt > (int _iblk, int _nlist, _Int * _list, _Int * _ia,
                               _Int * _ja, _Flt * _a, int _nblks, long long *_blks,
                               int &_icycle, int *_imaskblk);
      /// Init by data
      CBMatrix < _Int, _Flt > (int _blksize, int _iblk, int _nlist, _Int * _list,
                               _Int * _ia, _Int * _ja, _Flt * _a, int _nblks,
                               long long *_blks, int &_icycle, int *_imaskblk);
      /// Destructor
      ~CBMatrix () {
      };
// External functions
// Get/set functions
      /// Get nzblk
      int GetNzblk ()
      {
         return nzblk;
      };
      int GetNzblk () const
      {
         return nzblk;
      };
      /// Get hmatr_str
      CMatrix < int, float >*GetHMatrStr ()
      {
         return &hmatr_str;
      };
      const CMatrix < int, float >*GetHMatrStr () const
      {
         return &hmatr_str;
      };
      /// Get asub_arr
      vector < CMatrix < _Int, _Flt > >*GetASub () {
         return &asub_arr;
      };
      const vector < CMatrix < _Int, _Flt > >*GetASub () const
      {
         return &asub_arr;
      };
      /// Get asub_arr
      CMatrix < _Int, _Flt > *GetASubArr () {
         return &asub_arr[0];
      };
      const CMatrix < _Int, _Flt > *GetASubArr () const
      {
         return &asub_arr[0];
      };
      /// Get nzatot
      int GetNzatot ()
      {
         int nzatot = 0;
         CMatrix < _Int, _Flt > *pasub = &asub_arr[0];
         int i;
         for (i = 0; i < nzblk; i++)
            nzatot += pasub[i].GetNza ();
         return nzatot;
      };
      /// Set nzblk
      void SetNzblk (int _nzblk)
      {
         nzblk = _nzblk;
      };
      /// Resize ASub
      void ResizeASub (int _nzblk)
      {
         asub_arr.resize (_nzblk + 1);
      };
      /// Replace function
      void ReplaceFree (CBMatrix < _Int, _Flt > &_aa)
      {
         this->nzblk = _aa.nzblk;
         this->hmatr_str.ReplaceFree (_aa.hmatr_str);
         this->asub_arr.swap (_aa.asub_arr);
         _aa.Clean ();
      };
      /// Clean hmatrix data
      void Clean ()
      {
         CMatrix < int, float >HMatr_dummy;
         vector < CMatrix < _Int, _Flt > >ASub_dummy (1);
         hmatr_str.ReplaceFree (HMatr_dummy);
         asub_arr.swap (ASub_dummy);
         nzblk = 0;
      };
// Other functions
      /// Add replace operator for sparsity only (threads version)
      void AddReplaceSpThr (bool _b_use_thr, const CBMatrix < _Int, _Flt > &_aa);
      /// Add replace operator for pairs (threads version)
      void AddReplacePairsThr (bool _b_use_thr, const CBMatrix < _Int, _Flt > &_aa);
      /// Add replace operator (threads version)
      void AddReplaceBxBThr (bool _b_use_thr, int _blksize, const CBMatrix < _Int,
                             _Flt > &_aa);
      /// Add replace operator for pairs (threads version)
      void AddReplacePairsBxBThr (bool _b_use_thr, int _blksize, const CBMatrix < _Int,
                                  _Flt > &_aa);
      /// Split structural matrix into the set of matrices
      static void SplitMatrSpIntoHMatrSp (int _nblks, long long *_blks,
                                          const CMatrix < _Int, _Flt > &_amatr,
                                          CBMatrix < _Int, _Flt > &_ahmatr);
      /// Split matrix into the set of matrices
      static void SplitMatrIntoHMatr (int _nblks, long long *_blks, const CMatrix < _Int,
                                      _Flt > &_amatr, CBMatrix < _Int, _Flt > &_ahmatr);
      /// Split matrix into the set of matrices
      static void SplitMatrIntoHMatr_BxB (int _blksize, int _nblks, long long *_blks,
                                          const CMatrix < _Int, _Flt > &_amatr,
                                          CBMatrix < _Int, _Flt > &_ahmatr);
      /// Combine hmatrix into matrix
      void CombineHMatrIntoMatrSp (bool _b_is_char, int _nblks, long long *_blks,
                                   CMatrix < _Int, _Flt > &_amatr) const;
      /// Split structural rectangular matrix into the set of matrices
      static void SplitRectMatrSpIntoHMatrSp (int _nblksR, int *_blksR, int _nblksC,
                                              int *_blksC, const CMatrix < _Int,
                                              _Flt > &_amatr, CBMatrix < _Int,
                                              _Flt > &_ahmatr);
      /// Compute the symmetrized submatrices
      static void SymmetrizeSubmatrices (void *_comm, int _nblks, long long *_blks,
                                         int *_blk2cpu, CBMatrix < _Int,
                                         _Flt > *_hmatr_arr, CBMatrix < _Int,
                                         _Flt > *_hmatr_symm_arr);
      /// Compute the symmetrized submatrices
      static void SymmetrizeSubmatrices (void *_comm, int _nhblks, int *_hblk2cpu,
                                         int *_hblk2blks, int *_blk2hblks,
                                         long long *_blks, CBMatrix < _Int,
                                         _Flt > *_hmatr_arr, CBMatrix < _Int,
                                         _Flt > *_hmatr_symm_arr);
      /// Compute the symmetrized submatrices
      static void SymmetrizeSubmatrices_BxB (void *_comm, int _blksize, int _nhblks,
                                             int *_hblk2cpu, int *_hblk2blks,
                                             int *_blk2hblks, long long *_blks,
                                             CBMatrix < _Int, _Flt > *_hmatr_arr,
                                             CBMatrix < _Int, _Flt > *_hmatr_symm_arr);
      /// Compute the extended lists
      static void ExtendedLists (void *_comm, int _ncycle, int _nblks, long long *_blks,
                                 int *_blk2cpu, CBMatrix < _Int, _Flt > *_hmatr_arr,
                                 int *_nlist_ext_arr, vector < int >*_list_ext_arr);
      /// Compute the extended lists
      static void ExtendedLists (void *_comm, int _ncycle, int _nhblks, int *_hblk2cpu,
                                 int *_hblk2blk, int *_blk2hblk, int _nblks,
                                 long long *_blks, CBMatrix < _Int, _Flt > *_hmatr_arr,
                                 int *_nlist_ext_arr, vector < int >*_list_ext_arr);
      /// Get the extended submatrices
      static void GetExtendedSubmatrices (void *_comm, int _nblks, long long *_blks,
                                          int *_blk2cpu, CBMatrix < _Int,
                                          _Flt > *_hmatr_arr, int *_nlist_ext_arr,
                                          vector < int >*_list_ext_arr, CMatrix < _Int,
                                          _Flt > *_matr_ext_arr);
      /// Get the extended submatrices
      static void GetExtendedSubmatrices (void *_comm, int _nhblks, int *_hblk2cpu,
                                          int *_hblk2blk, int *_blk2hblk, int _nblks,
                                          long long *_blks, CBMatrix < _Int,
                                          _Flt > *_hmatr_arr, int *_nblks_ext_arr,
                                          vector < int >*_blksnum_ext_arr,
                                          vector < long long >*_pblks_ext_arr,
                                          int *_nlist_ext_arr,
                                          vector < int >*_list_ext_arr, CBMatrix < _Int,
                                          _Flt > *_hmatr_ext_arr);
      /// Get the extended submatrices
      static void GetExtendedSubmatrices_BxB (void *_comm, int _blksize, int _nhblks,
                                              int *_hblk2cpu, int *_hblk2blk,
                                              int *_blk2hblk, int _nblks,
                                              long long *_blks, CBMatrix < _Int,
                                              _Flt > *_hmatr_arr, int *_nblks_ext_arr,
                                              vector < int >*_blksnum_ext_arr,
                                              vector < long long >*_pblks_ext_arr,
                                              int *_nlist_ext_arr,
                                              vector < int >*_list_ext_arr,
                                              CBMatrix < _Int, _Flt > *_hmatr_ext_arr);
      /// Compute matrix decomposition with wells taken into account (MPI+thr parallelization)
      static void DecompWells (void *_comm, int _nneib_well, double _eps_dia_well,
                               double _thresh_max_well, int _ncycle_well, int _degree,
                               int _isize_max, int _isize_max2, int _nparts_bdecomp,
                               int _nparts, int _nparts_W, int _nhblks, int *_hblk2cpu,
                               int *_hblk2blk, long long *_blks,
                               vector < int >*_p_ia_wells_ext,
                               vector < int >*_p_ja_wells_ext, CBMatrix < _Int,
                               _Flt > *_hmatr_arr, int &_n_wells_new, int &_nhblks_new,
                               vector < int >&_hblk2cpu_new, vector < int >&_hblk2blk_new,
                               vector < int >&_blk2type_new,
                               vector < long long >&_blks_new, long long *_order);
      /// Find wells
      static void FindWells (void *_comm, int _nneib_well, double _eps_dia_well,
                             double _thresh_max_well, int _nhblks, int *_hblk2cpu,
                             int *_hblk2blk, long long *_blks, CBMatrix < _Int,
                             _Flt > *_hmatr_arr, vector < int >&_ia_wells_hblks,
                             vector < int >&_ja_3index_wells_hblks);
      /// Compute matrix decomposition (threads parallelization)
      static void DecompWeights_thr (bool _split_unconnected, int _degree, int _isize_max,
                                     int _nparts_decomp, int _nblks, long long *_blks,
                                     CBMatrix < _Int, _Flt > &_amatr, int &_nparts,
                                     vector < long long >&_parts, int *_order);
      /// Compute partitioning and binary tree according to sparsity and weights matrix
      static void DecompWeightsTree (CMatrix < _Int, _Flt > &_amatr_strW, int _npwr2_ext,
                                     CTree & _tree, int *_list2nd);
      /// Compute optimal ordering for the hmatrix via splitting into 3 sets of blocks
      static void ComputeOptimalOrderSchurWells (bool b_blk_wells, int _nparts_pwr2,
                                                 int _nblks, int _nblks1,
                                                 long long *_blks, CBMatrix < _Int,
                                                 _Flt > *_hmatr, CTree & _tree_0,
                                                 CTree & _tree_1, int &_nblks_new,
                                                 int &_nblks1_new, int &_nblks2_new,
                                                 vector < long long >&_blks_new,
                                                 vector < long long >&_nzord_new,
                                                 CVectorData < int >&_order);
      /// Compute optimal ordering for the hmatrix via splitting into 3 sets of blocks
      static void ComputeOptimalOrderSchur (int _nparts_pwr2, int _nblks, int _nblks1,
                                            long long *_blks, CBMatrix < _Int,
                                            _Flt > *_hmatr, CTree & _tree_0,
                                            CTree & _tree_1, int &_nblks_new,
                                            int &_nblks1_new, int &_nblks2_new,
                                            vector < long long >&_blks_new,
                                            vector < long long >&_nzord_new,
                                            CVectorData < int >&_order);
      /// Compute transposed hmatrix
      static void TransposeHMatrix (int _nblks, long long *_blks, CBMatrix < _Int,
                                    _Flt > &_hmatr_ini, CBMatrix < _Int,
                                    _Flt > &_hmatr_fin);
      /// Compute transposed hmatrix
      static void TransposeHMatrix_BxB (int _blksize, int _nblks, long long *_blks,
                                        CBMatrix < _Int, _Flt > &_hmatr_ini,
                                        CBMatrix < _Int, _Flt > &_hmatr_fin);
      /// Compute reordered hmatrix
      static void ReorderHMatrix (int _nblks_ini, long long *_blks_ini, CBMatrix < _Int,
                                  _Flt > &_hmatr_ini, int *_order, int _nblks_fin,
                                  long long *_blks_fin, CBMatrix < _Int,
                                  _Flt > &_hmatr_fin);
      /// Compute reordered hmatrix
      static void ReorderHMatrix_BxB (int _blksize, int _nblks_ini, long long *_blks_ini,
                                      const CBMatrix < _Int, _Flt > &_hmatr_ini,
                                      int *_order, int _nblks_fin, long long *_blks_fin,
                                      CBMatrix < _Int, _Flt > &_hmatr_fin);
      /// Compute inverse order
      static void InverseOrder (void *_comm, int _nhblks_ini, int *_hblk2cpu_ini,
                                int *_hblk2blk_ini, long long *_blks_ini,
                                long long *_order, int _nhblks_fin, int *_hblk2cpu_fin,
                                int *_hblk2blk_fin, long long *_blks_fin,
                                CVectorData < long long >&_iorder);
      /// Compute inverse order (local hblock diagonal ordering)
      static void InverseOrderDiag (void *_comm, int _nhblks_ini, int *_hblk2cpu_ini,
                                    int *_hblk2blk_ini, long long *_blks_ini,
                                    CVectorData < int >*_order_diag, int _nhblks_fin,
                                    int *_hblk2cpu_fin, int *_hblk2blk_fin,
                                    long long *_blks_fin,
                                    CVectorData < long long >&_iorder);
      /// Compute reordered sparsity hmatrix
      static void ReorderHMatrixSp (void *_comm, int _nhblks_ini, int *_hblk2cpu_ini,
                                    int *_hblk2blk_ini, long long *_blks_ini,
                                    CBMatrix < _Int, _Flt > *_hmatr_arr_ini,
                                    long long *_order, int _nhblks_fin,
                                    int *_hblk2cpu_fin, int *_hblk2blk_fin,
                                    long long *_blks_fin, CBMatrix < _Int,
                                    _Flt > *_hmatr_arr_fin);
      /// Compute reordered hmatrix
      static void ReorderHMatrix (void *_comm, int _nhblks_ini, int *_hblk2cpu_ini,
                                  int *_hblk2blk_ini, long long *_blks_ini,
                                  CBMatrix < _Int, _Flt > *_hmatr_arr_ini,
                                  long long *_order, int _nhblks_fin, int *_hblk2cpu_fin,
                                  int *_hblk2blk_fin, long long *_blks_fin,
                                  CBMatrix < _Int, _Flt > *_hmatr_arr_fin);
      /// Compute reordered hmatrix
      static void ReorderHMatrix_BxB (void *_comm, int _blksize, int _nhblks_ini,
                                      int *_hblk2cpu_ini, int *_hblk2blk_ini,
                                      long long *_blks_ini, CBMatrix < _Int,
                                      _Flt > *_hmatr_arr_ini, long long *_order,
                                      int _nhblks_fin, int *_hblk2cpu_fin,
                                      int *_hblk2blk_fin, long long *_blks_fin,
                                      CBMatrix < _Int, _Flt > *_hmatr_arr_fin);
      /// Compute reordered hmatrix rows
      static void ReorderHMatrixRows (int _nblksC, long long *_blksC, int _nblksR_ini,
                                      long long *_blksR_ini, CBMatrix < _Int,
                                      _Flt > &_hmatr_ini, int *_orderR, int _nblksR_fin,
                                      long long *_blksR_fin, CBMatrix < _Int,
                                      _Flt > &_hmatr_fin);
      /// Compute reordered hmatrix rows
      static void ReorderHMatrixRows_BxB (int _blksize, int _nblksC, long long *_blksC,
                                          int _nblksR_ini, long long *_blksR_ini,
                                          CBMatrix < _Int, _Flt > &_hmatr_ini,
                                          int *_orderR, int _nblksR_fin,
                                          long long *_blksR_fin, CBMatrix < _Int,
                                          _Flt > &_hmatr_fin);
      /// Compute reordered hmatrix cols
      static void ReorderHMatrixCols (int _nblksR, long long *_blksR, int _nblksC_ini,
                                      long long *_blksC_ini, CBMatrix < _Int,
                                      _Flt > &_hmatr_ini, int *_orderC, int _nblksC_fin,
                                      long long *_blksC_fin, CBMatrix < _Int,
                                      _Flt > &_hmatr_fin);
      /// Compute reordered hmatrix cols
      static void ReorderHMatrixCols_BxB (int _blksize, int _nblksR, long long *_blksR,
                                          int _nblksC_ini, long long *_blksC_ini,
                                          CBMatrix < _Int, _Flt > &_hmatr_ini,
                                          int *_orderC, int _nblksC_fin,
                                          long long *_blksC_fin, CBMatrix < _Int,
                                          _Flt > &_hmatr_fin);
      /// Compute reordered hmatrix (orderings are local inside hblocks)
      static void ReorderHMatrixDiag (void *_comm, int _nhblks, int *_hblk2cpu,
                                      int *_hblk2blk_ini, int _nblks_ini,
                                      long long *_blks_ini, CBMatrix < _Int,
                                      _Flt > *_hmatr_arr_ini,
                                      CVectorData < int >*_order_diag, int *_hblk2blk_fin,
                                      int _nblks_fin, long long *_blks_fin,
                                      CBMatrix < _Int, _Flt > *_hmatr_arr_fin);
      /// Compute reordered hmatrix (orderings are local inside hblocks)
      static void ReorderHMatrixDiag_BxB (void *_comm, int _blksize, int _nhblks,
                                          int *_hblk2cpu, int *_hblk2blk_ini,
                                          int _nblks_ini, long long *_blks_ini,
                                          CBMatrix < _Int, _Flt > *_hmatr_arr_ini,
                                          CVectorData < int >*_order_diag,
                                          int *_hblk2blk_fin, int _nblks_fin,
                                          long long *_blks_fin, CBMatrix < _Int,
                                          _Flt > *_hmatr_arr_fin);
      /// Compute comparison hmatrix for block hmatrix
      static void ComparisonHMatrix_BxB (int _blksize, CBMatrix < _Int,
                                         _Flt > &_hmatr_bxb, CBMatrix < _Int,
                                         _Flt > &_hmatr_pt);
      /// Compute block diagonal matrix
      static void BlockDiagonalHMatrix_BxB (int _blksize, int _nblks, long long *_blks,
                                            _Flt * _diag_data, CBMatrix < _Int,
                                            _Flt > &_hmatr_diag);
      /// Compute the packed size
      long long GetPackedSize () const;
      /// Fill char array by the packed data
      void FillPacked (long long _length, char *_obj) const;
      /// Fill char array by the packed data
      void FillPacked_thr (long long _length, char *_obj) const;
      /// Unpack data
      void UnPack (long long _length, char *_obj);
      /// Unpack data
      void UnPack_thr (long long _length, char *_obj);
// Low level support functions
      /// Compute reordered hmatrix (implementation)
      static void ReorderHMatrix_BxB_impl (bool _b_is_char, int _blksize, int _nblks_ini,
                                           long long *_blks_ini, const CBMatrix < _Int,
                                           _Flt > &_hmatr_ini, int *_order,
                                           int _nblks_fin, long long *_blks_fin,
                                           CBMatrix < _Int, _Flt > &_hmatr_fin);
      /// Compute matrix decomposition with wells taken into account (MPI+thr parallelization)
      static void DecompWells_impl (void *_comm, int _ncycle_well, int _degree,
                                    int _isize_max, int _isize_max2, int _nparts_bdecomp,
                                    int _nparts, int _nparts_W, int _nhblks,
                                    int *_hblk2cpu, int *_hblk2blk, long long *_blks,
                                    CBMatrix < _Int, _Flt > *_hmatr_arr,
                                    int *_ia_wells_hblks, int *_ja_wells_3index_hblks,
                                    int &_nhblks_new, vector < int >&_hblk2cpu_new,
                                    vector < int >&_hblk2blk_new,
                                    vector < int >&_blk2type_new,
                                    vector < long long >&_blks_new, long long *_order);
      /// Compute matrix decomposition (threads parallelization)
      static void DecompWeights_impl (bool _split_unconnected, int _degree,
                                      int _isize_max, int _nblks, long long *_blks,
                                      CBMatrix < _Int, _Flt > &_amatr, int _nparts,
                                      int *_partition);
      /// Prepare data for mvm functions
      static void PrepareMvmData (int _myid, int _nhblks, int *_hblk2cpu, int *_hblk2blk,
                                  long long *_blks, CBMatrix < _Int, _Flt > *_hmatr_arr,
                                  CMatrix < int, float >*_pAHBlkStr_arr,
                                  CVectorData < void *>*_ptr_arr,
                                  CVectorData < CVectorData < int > >*_listRC_offd);
      /// Perform in-place explicit block scaling
      static void ExplicitBlockScaling (void *_comm, int _blksize, int _nhblks,
                                        int *_hblk2cpu, int *_hblk2blk, long long *_blks,
                                        CBMatrix < _Int, _Flt > *_hmatr_arr,
                                        CBMatrix < _Int, _Flt > *_sclL_arr,
                                        CBMatrix < _Int, _Flt > *_sclU_arr);
      /// Filter the extended lists for backward only extention
      static void FilterListsBack (int _myid, int _nblks, int *_blk2cpu,
                                   int *_nlist_ext_arr, vector < int >*_list_ext_arr);
      /// Filter the extended lists for backward only extention
      static void FilterListsBack (int _myid, int _nhblks, int *_hblk2cpu,
                                   int *_blk2hblks, int *_nlist_ext_arr,
                                   std::vector < int >*_list_ext_arr);
      /// Compute Ja2 array by binary search
      static void ComputeJa2 (int _nblks, long long *_blks, int _nlist, const _Int * _ia,
                              const _Int * _ja, _Int * _ja2);
      /// Print sparsity with boxes that show nonzero blocks
      static void Str2PsBox (const CMatrix < _Int, _Flt > &_amatr, const char *_fname,
                             int _nblks, int *_blks);
      /// Print sparsity with boxes that show nonzero blocks
      static void Str2PsBox (int _collap, const CMatrix < _Int, _Flt > &_amatr,
                             const char *_fname, int _nblks, int *_blks);
      /// Print sparsity with boxes that show nonzero blocks
      static void StrRect2PsBox (const CMatrix < _Int, _Flt > &_amatr, const char *_fname,
                                 int _nblksR, int *_blksR, int _nblksC, int *_blksC);
      /// Print sparsity with boxes that show nonzero blocks
      static void StrRect2PsBox (int _collap, const CMatrix < _Int, _Flt > &_amatr,
                                 const char *_fname, int _nblksR, int *_blksR,
                                 int _nblksC, int *_blksC);
      /// Print sparsity with boxes that show nonzero blocks
      static void Str2PsBox2 (const CMatrix < _Int, _Flt > &_amatr, const char *_fname,
                              int _nhblks, int *_hblk2blks, int *_blks);
      /// Print sparsity with boxes that show nonzero blocks
      static void Str2PsBox2 (int _collap, const CMatrix < _Int, _Flt > &_amatr,
                              const char *_fname, int _nhblks, int *_hblk2blks,
                              int *_blks);
      /// Print hmatrix data
      void PrintHMatrix (ofstream & _fout);
      /// Print hmatrix data over thresh
      void PrintHMatrixThresh (ofstream & _fout, double _thresh);
      /// Print diag of hmatrix data
      void PrintHMatrixDiag (ofstream & _fout);
      /// Print diag of hmatrix data in pairs format
      void PrintHMatrixDiagPair (ofstream & _fout);
   };

   template < typename _Int, typename _Flt, typename _Flt2 > class CMatrixConv  // Class that supports conversion of data
   {
    public:
      /// Convert matrix data
      static void InitAndConv (int _n, _Int * _ia, _Int * _ja, _Flt * _a, CMatrix < _Int,
                               _Flt2 > &_amatr);
      static void InitAndConv (CMatrix < _Int, _Flt > &_amatr_ini, CMatrix < _Int,
                               _Flt2 > &_amatr_fin);
      static void InitAndConv (CBMatrix < _Int, _Flt > &_ahmatr_ini, CBMatrix < _Int,
                               _Flt2 > &_ahmatr_fin);
   };

/// Class that supports base qrd
   template < typename _FltVect > class CQrdBase {
    protected:
      int nqblk;                ///< Current number of q blocks
      int nqblk_alloc;          ///< Allocated number of q blocks
      vector < int >qblksc;     ///< IA-type array that describes the number of columns in each q block
      vector < int >ncolarr_alloc;      ///< The list of allocate q column sizes 
      vector < int >nrowarr;    ///< The number of rows in each q block
      vector < int >nrowarr_alloc;      ///< The number of allocated rows in each q block
      vector < CVectorData < _FltVect > >qarr;  ///< List of q blocks
      vector < CVectorData < _FltVect > >rdiag_arr;     ///< List of rdiag blocks
      vector < vector < _FltVect > >tauarr;     ///< List of tau blocks
    public:
// Functions
// Constructors and destructor
      /// Default constructor
      CQrdBase () {
         nqblk = 0;
         nqblk_alloc = 0;
         qblksc.resize (2);
         qblksc[0] = 0;
         ncolarr_alloc.resize (1);
         nrowarr.resize (1);
         nrowarr_alloc.resize (1);
         qarr.resize (1);
         rdiag_arr.resize (1);
         tauarr.resize (1);
      };
      /// Copy constructor
    private:
      CQrdBase (const CQrdBase & _aa)
      {
      };
      /// Copy operator
      CQrdBase & operator= (const CQrdBase & _aa)
      {
         return *this;
      };
    public:
      /// Destructor
      virtual ~ CQrdBase () {
      };
// Operator functions
// Get/set functions
      /// Get nqblk
      int GetNqblk () const
      {
         return nqblk;
      };
      /// Get qblksc
      const int *GetQblksc () const
      {
         return &qblksc[0];
      };
      /// Get qblksc
      int *GetQblksc ()
      {
         return &qblksc[0];
      };
      /// Get ncolarr_alloc
      const int *GetNcolarr_alloc () const
      {
         return &ncolarr_alloc[0];
      };
      /// Get ncolarr_alloc
      int *GetNcolarr_alloc ()
      {
         return &ncolarr_alloc[0];
      };
      /// Get nrowarr
      const int *GetNrowarr () const
      {
         return &nrowarr[0];
      };
      /// Get nrowarr
      int *GetNrowarr ()
      {
         return &nrowarr[0];
      };
      /// Get qarr
      const CVectorData < _FltVect > *GetQarr () const
      {
         return &qarr[0];
      };
      /// Get qarr
      CVectorData < _FltVect > *GetQarr () {
         return &qarr[0];
      };
      /// Get rdiag_arr
      const CVectorData < _FltVect > *GetRDiagarr () const
      {
         return &rdiag_arr[0];
      };
      /// Get rdiag_arr
      CVectorData < _FltVect > *GetRDiagarr () {
         return &rdiag_arr[0];
      };
      /// Get tauarr
      const vector < _FltVect > *GetTauarr () const
      {
         return &tauarr[0];
      };
      /// Get tauarr
      vector < _FltVect > *GetTauarr () {
         return &tauarr[0];
      };
      /// Set nqblk
      virtual void SetNqblk (int _nqblk)
      {
         nqblk = _nqblk;
      };
// Qrd up level functions
      int GetNRows ();          ///< Get number of rows
      int GetNCols ();          ///< Get number of columns
      virtual int GetAllocatedMemory ();        ///< Get memory
      /// Update QR decomposition for current set of columns
      void UpdateQrdBlk (int _ncol, int _nrow, _FltVect * _ablk, int _lda);
      /// Update QR decomposition for current set of columns (threads version)
      void UpdateQrdBlk_thr (int _ncol, int _nrow, _FltVect * _ablk, int _lda);
      /// Multiply Q by given block (self version)
      void MvmQ (int _nrhs, _FltVect * _qx, int _ldqx);
      /// Multiply Q by given block (self version) (threads version)
      void MvmQ_thr (int _nrhs, _FltVect * _qx, int _ldqx);
      /// Multiply Q by given block
      void MvmQ (int _nrhs, int _nrows, _FltVect * _x, int _ldx, _FltVect * _qx,
                 int _ldqx);
      /// Multiply Q by given block (threads version)
      void MvmQ_thr (int _nrhs, int _nrows, _FltVect * _x, int _ldx, _FltVect * _qx,
                     int _ldqx);
      /// Multiply Qh by given block (self version)
      void MvmQH (int _nrhs, _FltVect * _qx, int _ldqx);
      /// Multiply Qh by given block (self version) (threads version)
      void MvmQH_thr (int _nrhs, _FltVect * _qx, int _ldqx);
      /// Multiply Qh by given block (self version)
      void MvmQHPart (int _nrhs, int _ibegQ, int _iendQ, _FltVect * _qx, int _ldqx);
      /// Multiply Qh by given block
      void MvmQH (int _nrhs, _FltVect * _x, int _ldx, _FltVect * _qx, int _ldqx);
      /// Multiply Qh by given block (threads version)
      void MvmQH_thr (int _nrhs, _FltVect * _x, int _ldx, _FltVect * _qx, int _ldqx);
      /// Get R part Qrd
      void GetRQrd (int _ibegc, int _iendc, int _ibegr, int _iendr, _FltVect * _r,
                    int _ldr);
      /// Get R part Qrd (implementation)
      void GetRQrd_impl (int _ibegc, int _iendc, int _ibegr, int _iendr, _FltVect * _r,
                         int _ldr);
      /// Store R part Qrd
      void StoreR (int _ibegc, int _iendc, int _ibegr, int _iendr, _FltVect * _r,
                   int _ldr);
// Internal functions
    private:
      /// Update QR decomposition for current set of columns (implementation)
      void UpdateQrdBlk_impl (bool _b_use_threads, int _ncol, int _nrow, _FltVect * _ablk,
                              int _lda);
      /// Multiply Q by given block (self version) (implementation)
      void MvmQ_impl (bool _b_use_threads, int _nrhs, _FltVect * _qx, int _ldqx);
      /// Multiply Q by given block
      void MvmQ_impl (bool _b_use_threads, int _nrhs, int _nrows, _FltVect * _x, int _ldx,
                      _FltVect * _qx, int _ldqx);
      /// Multiply Qh by given block (self version) (implementation)
      void MvmQH_impl (bool _b_use_threads, int _nrhs, _FltVect * _qx, int _ldqx);
      /// Multiply Qh by given block (implementation)
      void MvmQH_impl (bool _b_use_threads, int _nrhs, _FltVect * _x, int _ldx,
                       _FltVect * _qx, int _ldqx);
   };

/// Class that supports QR decomposition with rows splitting
   template < typename _FltVect > class CQrdSet {
    protected:
      int nqrdsets;             ///< number of subsets
      CQrdBase < _FltVect > *pqrd_head; ///< pointer to head qrd structure
      CQrdBase < _FltVect > **ppqrd_childs;     ///< set of pointers to childs of QR
    public:
// Functions
// Constructors and destructor
      /// Default
      CQrdSet () {
         nqrdsets = 0;
      };
      /// Destructor
      virtual ~ CQrdSet () {
      };
    private:
      /// Copy
      CQrdSet (const CQrdSet & _aa)
      {
      };
      /// Copy operator
      CQrdSet & operator= (const CQrdSet & _aa)
      {
         return *this;
      };
    public:
// Get/set interface functions
      /// Get nqrdsets
      int GetNqrdsets () const
      {
         return nqrdsets;
      };
      /// Get pqrd_head
      CQrdBase < _FltVect > *GetPQrdHead ()const
      {
         return pqrd_head;
      };
      /// Get ppqrd_childs
      CQrdBase < _FltVect > **GetPPQrdChilds ()const
      {
         return ppqrd_childs;
      };
      /// Set nqrdsets
      void SetNqrdsets (int _nqrdsets)
      {
         nqrdsets = _nqrdsets;
      };
      /// Set pqrd_head
      void SetPQrdHead (CQrdBase < _FltVect > *_pqrd_head)
      {
         pqrd_head = _pqrd_head;
      };
      /// Set ppqrd_childs
      void SetPPQrdChilds (CQrdBase < _FltVect > **_ppqrd_childs)
      {
         ppqrd_childs = _ppqrd_childs;
      };
// QrdSet up level functions
      void UpdateQrdHead (int _ncol);   ///< Update head QR
      /// Multiply by head Q and prepare childs
      void MvmQ (int _nrhs, _FltVect * _x, int _ldx, _FltVect * _qx, int _ldqx);
   };

/// Class that supports MPI QR decomposition
   template < typename _FltVect > class CQrdMPI {
    protected:
      void *pComm;              ///< Pointer to MPI comm
      int ni_myid;              ///< Vector size of QRD for current processor
      CTree treeMPI;            ///< MPI tree
      int node_myid_up;         ///< Top level node for myid
      vector < int >nd2cpu;     ///< Number of cpu for each node
      vector < int >cpu2nd;     ///< Number of kernel node for each cpu
      CQrdSet < _FltVect > *pqrd_sets;  ///< Qrd sets
      CQrdBase < _FltVect > *pqrd_childs;       ///< Qrd childs
      CQrdBase < _FltVect > **ppqrd_childs;     ///< Pointers to Qrd childs
    public:
// Functions
// Constructors and destructor
      /// Default
      CQrdMPI () {
         pComm = NULL;
         ni_myid = 0;
         node_myid_up = -1;
         nd2cpu.resize (1);
         cpu2nd.resize (1);
         pqrd_sets = NULL;
         pqrd_childs = NULL;
         ppqrd_childs = NULL;
      };
      /// Destructor
      virtual ~ CQrdMPI () {
         if (pqrd_sets != NULL)
            delete[]pqrd_sets;
         if (pqrd_childs != NULL)
            delete[]pqrd_childs;
         if (ppqrd_childs != NULL)
            delete[]ppqrd_childs;
      };
    private:
      /// Copy
      CQrdMPI (const CQrdMPI & _aa)
      {
      };
      /// Copy operator
      CQrdMPI & operator= (const CQrdMPI & _aa)
      {
         return *this;
      };
    public:
// Get/set interface functions
      /// Get number of rows
      int GetNRows ()
      {
         return ni_myid;
      };
      /// Get tree
      CTree *GetTree ()
      {
         return &treeMPI;
      };
      /// Get qrd_sets
      CQrdSet < _FltVect > *GetQrdSets () {
         return pqrd_sets;
      };
      /// Get qrd_childs
      CQrdBase < _FltVect > *GetQrdChilds () {
         return pqrd_childs;
      };
// Qrd up level functions
      void Init (void *_pComm, int _ni_myid);   ///< Init the structure
      int GetNCols ();          ///< Get the current number of columns
      /// Update QR decomposition
      void UpdateQrdMPI (int _ncol, int _nrows, _FltVect * _ablk, int _lda);
      /// Update QR decomposition (threads version)
      void UpdateQrdMPI_thr (int _ncol, int _nrows, _FltVect * _ablk, int _lda);
      /// Multiply by Q
      void MvmQMPI (int _nrhs, _FltVect * _x, int _ldx, _FltVect * _qx, int _ldqx);
      /// Multiply by Q (threads version)
      void MvmQMPI_thr (int _nrhs, _FltVect * _x, int _ldx, _FltVect * _qx, int _ldqx);
      /// Get R part
      void GetRQrdMPI (int _ibegc, int _iendc, int _ibegr, int _iendr, _FltVect * _r,
                       int _ldr);
      /// Free qrd
      void FreeQrdMPI ();
// Internal functions
    private:
      /// Update QR decomposition (implementation)
      void UpdateQrdMPI_impl (bool _b_use_threads, int _ncol, int _nrows,
                              _FltVect * _ablk, int _lda);
      /// Multiply by Q (implementation)
      void MvmQMPI_impl (bool _b_use_threads, int _nrhs, _FltVect * _x, int _ldx,
                         _FltVect * _qx, int _ldqx);
   };

   template < typename _Int, typename _Flt, typename _FltVect > class CMvmPar   // Class that supports parallel multiplication computations
   {
      void *pcomm;              ///< Pointer to MPI communicator
      int nblks;                ///< Number of blocks
      long long *pblks;         ///< Pointer to blocks partitioning
      int *pblk2cpu;            ///< Pointer to block to cpu distribution
      int ni_cpu;               ///< The size of the local part of vector data
      vector < int >ibsblk;     ///< Base adresses of local parts of the vector data
      int nlistblk_own;         ///< The number of local blocks
      vector < int >listblk_own;        ///< The list of own blocks
      CBMatrix < _Int, _Flt > *phmatr;  ///< Pointer to hyperblocks array
      int nsends;               ///< The number of sends for computations
      vector < int >snd2cpu;    ///< The list of cpus to be sended to for computations
      vector < int >ia_sends;   ///< The sizes of sends data for computations
      vector < int >ind_sends;  ///< The indices of local data for sends
      vector < _FltVect > x_send;       ///< The work memory for sends
      int nrecvs;               ///< The number of recvs for computations
      vector < int >rcv2cpu;    ///< The list of cpus to be received data from for computations
      vector < int >ia_recvs;   ///< The sizes of recvs for computations
      vector < int >ja_recvs;   ///< The pairs of local data that describe recvs
      vector < _FltVect > x_recv;       ///< The work memory for recvs
      int nblks_recvs;          ///< The number of received block parts
      vector < int >listblk_recvs;      ///< The list of block numbers of received blocks
      vector < int >ialist_recvs;       ///< The numbers of block rows contatining the block
      vector < int >japairs_recvs;      ///< The complete list and index numbers inside of block rows containing each of the received blocks
      vector < int >iablk_recvs;        ///< The sizes of copy data for each received block
      vector < int >ind_recvs;  ///< The list of indices inside block for received blocks
      vector < _FltVect > x_temp;       ///< The work memory for external blocks
// Constructors and destructor
    private:
      /// Copy constructor
      CMvmPar (const CMvmPar < _Int, _Flt, _FltVect > &_aa)
      {
      };
      /// Copy operator
      CMvmPar < _Int, _Flt, _FltVect > &operator= (const CMvmPar < _Int, _Flt,
                                                   _FltVect > &_aa)
      {
         return *this;
      };
    public:
      /// Epmty constructor
      CMvmPar () {
         pcomm = NULL;
         nblks = 0;
         pblks = NULL;
         pblk2cpu = NULL;
         phmatr = NULL;
         ni_cpu = 0;
         nlistblk_own = 0;
         nsends = 0;
         nrecvs = 0;
         nblks_recvs = 0;
      };
      /// Destructor
      ~CMvmPar () {
      };
// External functions
      /// Init control data
      void InitControl (void *_pcomm, int _nblks, long long *_blks, int *_blk2cpu)
      {
         pcomm = _pcomm;
         int myid = CMPIDataExchange::GetMyid (pcomm);
         nblks = _nblks;
         pblks = _blks;
         pblk2cpu = _blk2cpu;
         ibsblk.resize (nblks + 1);
         int *pibsblk = &ibsblk[0];
         ni_cpu = 0;
         nlistblk_own = 0;
         int i;
         for (i = 0; i < nblks; i++) {
            if (pblk2cpu[i] == myid) {
               nlistblk_own++;
               pibsblk[i] = ni_cpu;
               ni_cpu += (int) (pblks[i + 1] - pblks[i]);
            } else {
               pibsblk[i] = -1;
            };
         };
         listblk_own.resize (nlistblk_own + 1);
         int *plistblk_own = &listblk_own[0];
         nlistblk_own = 0;
         for (i = 0; i < nblks; i++) {
            if (pblk2cpu[i] == myid) {
               plistblk_own[nlistblk_own] = i;
               nlistblk_own++;
            }
         };
      };
      /// Init MvmA data
      void InitMvmA (CBMatrix < _Int, _Flt > *_hmatr_arr);
      /// Perform MvmA computations
      void MvmA (const _FltVect * _x, _FltVect * _ax);
      /// Clean MvmA structure
      void Clean ();
   };

   template < typename _Int, typename _Flt, typename _FltVect > class CSlvPar   // Class that supports parallel solve computations (including L, U and ordering)
   {
      void *pcomm;              ///< Pointer to MPI communicator
      int nblks;                ///< Number of blocks
      long long *pblks;         ///< Pointer to blocks partitioning
      long long *pblks_ext;     ///< Pointer to the extended blocks partitioning
      int *pblk2cpu;            ///< Pointer to block to cpu distribution
      int ni_cpu;               ///< The size of the local part of vector data
      vector < int >ibsblk;     ///< Base adresses of local parts of the vector data
      int nlistblk_own;         ///< The number of local blocks
      vector < int >listblk_own;        ///< The list of own blocks
      vector < int >*plistpairs_ext;    ///< Pointer to the extended data lists
      CMatrix < _Int, _Flt > *pmatrL;   ///< Pointer to L blocks array
      CMatrix < _Int, _Flt > *pmatrU;   ///< Pointer to U blocks array
      vector < int >*porderLU;  ///< Pointer to the extended ordering data if any
      int nsends;               ///< The number of sends for computations
      vector < int >snd2cpu;    ///< The list of cpus to be sended to for computations
      vector < int >ia_sends;   ///< The sizes of sends data for computations
      vector < int >ind_sends;  ///< The indices of local data for sends
      vector < _FltVect > x_send;       ///< The work memory for sends
      int nrecvs;               ///< The number of recvs for computations
      vector < int >rcv2cpu;    ///< The list of cpus to be received data from for computations
      vector < int >ia_recvs;   ///< The sizes of recvs for computations
      vector < _FltVect > x_recv;       ///< The work memory for recvs
      vector < int >ialist_recvs;       ///< For each local block the number of updates from received data
      vector < int >jalist_recvs;       ///< The lists of index numbers for received data for updates in current local block
      vector < _FltVect > x_temp;       ///< The work memory for extended blocks
// Constructors and destructor
    private:
      /// Copy constructor
      CSlvPar (const CSlvPar < _Int, _Flt, _FltVect > &_aa)
      {
      };
      /// Copy operator
      CSlvPar < _Int, _Flt, _FltVect > &operator= (const CSlvPar < _Int, _Flt,
                                                   _FltVect > &_aa)
      {
         return *this;
      };
    public:
      /// Epmty constructor
      CSlvPar () {
         pcomm = NULL;
         nblks = 0;
         pblks = NULL;
         pblk2cpu = NULL;
         plistpairs_ext = NULL;
         pmatrL = NULL;
         pmatrU = NULL;
         porderLU = NULL;
         ni_cpu = 0;
         nlistblk_own = 0;
         nsends = 0;
         nrecvs = 0;
      };
      /// Destructor
      ~CSlvPar () {
      };
// External functions
      /// Get ni_cpu
      int GetNiCpu ()
      {
         return ni_cpu;
      };
      /// Init control data
      void InitControl (void *_pcomm, int _nblks, long long *_blks, long long *_blks_ext,
                        int *_blk2cpu)
      {
         pcomm = _pcomm;
         int myid = CMPIDataExchange::GetMyid (pcomm);
         nblks = _nblks;
         pblks = _blks;
         pblks_ext = _blks_ext;
         pblk2cpu = _blk2cpu;
         ibsblk.resize (nblks + 1);
         int *pibsblk = &ibsblk[0];
         ni_cpu = 0;
         nlistblk_own = 0;
         int i;
         for (i = 0; i < nblks; i++) {
            if (pblk2cpu[i] == myid) {
               nlistblk_own++;
               pibsblk[i] = ni_cpu;
               ni_cpu += (int) (pblks[i + 1] - pblks[i]);
            } else {
               pibsblk[i] = -1;
            };
         };
         listblk_own.resize (nlistblk_own + 1);
         int *plistblk_own = &listblk_own[0];
         nlistblk_own = 0;
         for (i = 0; i < nblks; i++) {
            if (pblk2cpu[i] == myid) {
               plistblk_own[nlistblk_own] = i;
               nlistblk_own++;
            }
         };
      };
      /// Init SolveLU data
      void InitSolveLU (vector < int >*_listpairs_ext, CMatrix < _Int, _Flt > *_matrL,
                        CMatrix < _Int, _Flt > *_matrU, vector < int >*_orderLU);
      /// Perform SolveLU computations
      void SolveLU (const _FltVect * _x, _FltVect * _px);
      /// Clean Slv structure
      void Clean ();
   };

   template < typename _Int, typename _Flt, typename _FltVect > class CK3D_Solver       // Class that supports parallel solver computations
   {
      void *pcomm;              ///< Pointer to MPI communicator
      SParams params;           ///< Solver params
      int nblks;                ///< Number of blocks
      vector < long long >blks; ///< Pointer to blocks partitioning
      vector < long long >blks_ext;     ///< Pointer to the extended blocks partitioning
      vector < int >blk2cpu;    ///< Pointer to block to cpu distribution
      vector < CBMatrix < _Int, _Flt > >hmatr_arr;      ///< Matrix data
      vector < int >nlist_ext_arr;      ///< The number of elems in extension
      vector < vector < int > >list_ext_arr;     ///< The pairs of extension indices
      vector < vector < int > >order_LU; ///< Ordering data
      vector < CMatrix < _Int, float > >matrL_float;     ///< Float  L part
      vector < CMatrix < _Int, double > >matrL_double;   ///< Double L part
      vector < CMatrix < _Int, float > >matrU_float;     ///< Float  U part
      vector < CMatrix < _Int, double > >matrU_double;   ///< Double U part
      CMvmPar < _Int, _Flt, _FltVect > mvm;     ///< Support MvmA
      CSlvPar < _Int, float, _FltVect > slv_float;      ///< Support float  Slv
      CSlvPar < _Int, double, _FltVect > slv_double;    ///< Support double Slv
      int ncoef_slv;            ///< Number of coefs in poly prec
      vector < double >coef_slv;        ///< Coefs in poly prec
      vector < double >xwork;   ///< Work memory for polynomial preconditioning
// Constructors and destructor
    private:
      /// Copy constructor
      CK3D_Solver (const CK3D_Solver < _Int, _Flt, _FltVect > &_aa)
      {
      };
      /// Copy operator
      CK3D_Solver < _Int, _Flt, _FltVect > &operator= (const CK3D_Solver < _Int, _Flt,
                                                       _FltVect > &_aa)
      {
         return *this;
      };
    public:
      /// Epmty constructor
      CK3D_Solver () {
         nblks = 0;
         ncoef_slv = 0;         /*cout<<"CK3D_Solver: empty constructor\n"; */
      };
      /// Destructor
      ~CK3D_Solver () {         /*cout<<"CK3D_Solver: empty destructor\n"; */
      };
// Get/set functions
      /// Get comm
      void *GetComm ()
      {
         return pcomm;
      };
      /// Get nblks
      int GetNblks ()
      {
         return nblks;
      };
      /// Get blks
      long long *GetBlks ()
      {
         return &blks[0];
      };
      /// Get blk2cpu
      int *GetBlk2cpu ()
      {
         return &blk2cpu[0];
      };
      /// Set ncoef_slv
      void SetNcoef (int _ncoef)
      {
         ncoef_slv = _ncoef;
      };
      /// Get coef
      vector < double >*GetCoef ()
      {
         return &coef_slv;
      };
// External functions
      /// Prepare solver structures including performing parallel fct
      void PrepareSolver (void *_pcomm, SParams * _params, int _nblks, long long *_blks,
                          int *_blk2cpu, bool _b_store_matrix, _Int * _ia, _Int * _ja,
                          _Flt * _a, double &_prec_extend, double &_density,
                          double &_scpiv_min, double &_scpiv_max, int &_nmodif,
                          double &_piv_min, double &_piv_max, double &_dtime_fct)
      {
         PrepareMatrix (_pcomm, _nblks, _blks, _blk2cpu, _ia, _ja, _a);
         ComputeBILU2 (_params, _prec_extend, _density, _scpiv_min, _scpiv_max, _nmodif,
                       _piv_min, _piv_max, _dtime_fct);
         if (!_b_store_matrix)
            CleanMvmA ();
      };
      /// Prepare matrix data only
      void PrepareMatrix (void *_pcomm, int _nblks, long long *_blks, int *_blk2cpu,
                          _Int * _ia, _Int * _ja, _Flt * _a);
      /// Peform BILU2 computations
      void ComputeBILU2 (SParams * _params, double &_prec_extend, double &_density,
                         double &_scpiv_min, double &_scpiv_max, int &_nmodif,
                         double &_piv_min, double &_piv_max, double &_dtime_fct);
      /// Perform iterations of the iterative scheme
      void SolveIter (int _ittype, int _niter_max, int _niter_cycle, int _ncoef,
                      int _niter_cycle2, double _eps, int _ichk, int _msglev,
                      ofstream * _fout, _FltVect * _rhs, _FltVect * _sol,
                      double &_rhs_norm, double &_res_ini, int &_niter, int &_nmvm,
                      double &_res_fin, double &_dtime_iter);
      /// Perform iterations of the BiCGStab iterative scheme
      void BiCGStab (bool _b_use_poly, int _niter_max, double _eps, int _ichk,
                     int _msglev, ofstream * _fout, _FltVect * _rhs, _FltVect * _sol,
                     double &_rhs_norm, double &_res_ini, int &_niter, int &_nmvm,
                     double &_res_fin, double &_dtime_iter);
      /// Perform iterations of the GMRES iterative scheme
      void Gmres (bool _b_use_poly, int _niter_max, int _niter_cycle, int _ncoef,
                  double _eps, int _ichk, int _msglev, ofstream * _fout, _FltVect * _rhs,
                  _FltVect * _sol, double &_rhs_norm, double &_res_ini, int &_niter,
                  int &_nmvm, double &_res_fin, double &_dtime_iter);
      /// Clean MvmA structures
      void CleanMvmA ()
      {
         mvm.Clean ();
         vector < CBMatrix < _Int, _Flt > >hmatr_arr_dummy;
         hmatr_arr.swap (hmatr_arr_dummy);
      };
      /// Clean Slv structures
      void CleanSlv ()
      {
         vector < int >nlist_ext_arr_dummy;
         vector < vector < int > >list_ext_arr_dummy;
         vector < vector < int > >order_LU_dummy;
         vector < CMatrix < _Int, float > >matrL_float_dummy;
         vector < CMatrix < _Int, double > >matrL_double_dummy;
         vector < CMatrix < _Int, float > >matrU_float_dummy;
         vector < CMatrix < _Int, double > >matrU_double_dummy;
         vector < double >coef_temp;
         vector < double >xwork_temp;
         nlist_ext_arr.swap (nlist_ext_arr_dummy);
         list_ext_arr.swap (list_ext_arr_dummy);
         order_LU.swap (order_LU_dummy);
         matrL_float.swap (matrL_float_dummy);
         matrL_double.swap (matrL_double_dummy);
         matrU_float.swap (matrU_float_dummy);
         matrU_double.swap (matrU_double_dummy);
         slv_float.Clean ();
         slv_double.Clean ();
         ncoef_slv = 0;
         coef_slv.swap (coef_temp);
         xwork.swap (xwork_temp);
      };
      /// Clean all structures
      void Clean ()
      {
         pcomm = NULL;
         nblks = 0;
         vector < long long >blks_dummy;
         vector < long long >blks_ext_dummy;
         vector < int >blk2cpu_dummy;
         blks.swap (blks_dummy);
         blks_ext.swap (blks_ext_dummy);
         blk2cpu.swap (blk2cpu_dummy);
         CleanMvmA ();
         CleanSlv ();
         params.SetDefaults ();
      };
// Internal functions
      /// Perform MvmA parallel computations
      void MvmA (const _FltVect * _x, _FltVect * _ax)
      {
         mvm.MvmA (_x, _ax);
      };
      /// Perform polynomial preconditioning
      void PolySlvLU (bool _b_use_poly, const _FltVect * _x, _FltVect * _px);
      /// Perform SlvLU parallel computations
      void SlvLU (const _FltVect * _x, _FltVect * _px)
      {
         if (params.prec_float == 1) {
            slv_float.SolveLU (_x, _px);
         } else {
            slv_double.SolveLU (_x, _px);
         };
      };
   };

}                               // namespace k3d

#endif
