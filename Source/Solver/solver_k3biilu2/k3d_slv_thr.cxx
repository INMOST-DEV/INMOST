//------------------------------------------------------------------------------------------------
// File: k3d_slv_thr.cxx
//------------------------------------------------------------------------------------------------
#ifndef K3D_SLV_THR_CXX
#define K3D_SLV_THR_CXX

#include "k3d_slv_thr.hxx"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iosfwd>
#include <cstring>
#include <iomanip>
#include <algorithm>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <map>

namespace k3d
{

//
// Perform ILU2 point factorization of the hmatrix with dinamic ordering and future diagonal modifications
//========================================================================================
   template < typename _Int, typename _Flt > void CFctThreads < _Int,
      _Flt >::Ilu2HMatrix (bool _b_blk_wells, SParams & _params, int _nblks, int _nblks1,
                           int _nblks2, long long *_blks, long long *_nzord_blks,
                           CTree & _tree0, CTree & _tree1, CBMatrix < _Int,
                           _Flt > &_a_matr, CBMatrix < _Int, _Flt > &_l_matr,
                           CBMatrix < _Int, _Flt > &_u_matr, CTree & _tree2_new,
                           int &_nblks_new, vector < long long >&_blks_new,
                           CVectorData < int >&_ordernew, double &_sclmin_att,
                           double &_sclmax_att, int &_nmodif, double &_eigmin_att,
                           double &_eigmax_att)
   {

      int ntot = (int) _blks[_nblks];

      int iparam1_out = _params.iparam1;

//      int iparam2_out = _params.iparam2;

//      char strbuff[256];
//      sprintf (strbuff,"ChkIlu2HMatrix_1x1_%i.dat",iparam1_out);
//      ofstream ffout (strbuff);

//      ffout << " nblks = " << _nblks << " nblks1 = " << _nblks1 << " nblks2 = " << _nblks2 << endl;

//      PrintArray (ffout, " Blks = ",_nblks+1,_blks);

//      ffout << " Tree 0 " << endl;
//      _tree0.OutputTree (ffout);
//      ffout << " Tree 1 " << endl;
//      _tree1.OutputTree (ffout);

// Get params

      int ordlevel = _params.ordlevel;
      int ordtype = _params.ordtype;
      int sctype = _params.sctype;
      double sclmin = _params.sclmin;
      int nitersc = _params.nitersc;
      int strtype = _params.strtype;
      int fcttype = _params.fcttype;
      double pivmin = _params.pivmin;
      double tau1 = _params.tau1;
      double tau2 = _params.tau2;
      double tau2_sch = _params.tau2_sch;
      double theta = _params.theta;

// Compute explicit scaling

//   ffout << " Ilu2HMatrix: Point 1" << endl;

        CVectorData < _Flt > scl_L;
        CVectorData < _Flt > scl_U;

        CFctThreads < _Int, _Flt >::ComputeExplicitScaling (sctype, nitersc, sclmin,
                                                            _nblks, _blks, _a_matr, scl_L,
                                                            scl_U, _sclmin_att,
                                                            _sclmax_att);

      _Flt *pscl_L = scl_L.Ptr ();
      _Flt *pscl_U = scl_U.Ptr ();

//      PrintArray (ffout, " Computed SclL ", ntot, pscl_L);
//      PrintArray (ffout, " Computed SclU ", ntot, pscl_U);

// Perform explicit scaling

//      ffout << " Ilu2HMatrix: Point 2" << endl;

//      ffout << " Diag bef scl" << endl;
//      if (ntot < 5000) _a_matr.PrintHMatrix (ffout);

        CFctThreads < _Int, _Flt >::HMatrixScale (_nblks, _blks, _a_matr, pscl_L, pscl_U);

//      ffout << " Diag aft scl" << endl;
//      if (ntot < 5000) _a_matr.PrintHMatrix (ffout);

// Split and combine by pairs

        CBMatrix < _Int, _Flt > alu_pair;

        CFctThreads < _Int, _Flt >::SplitLUPair (_nblks, _blks, _a_matr, alu_pair);

// Perform fct

        CBMatrix < _Int, _Flt > lu_pair;

//   ffout << " Ilu2HMatrix: Point 3" << endl;

//   CFctThreads<_Int,_Flt>::Ilu2BlockIlu2 (&ffout, iparam1_out, ordlevel, ordtype, 
        CFctThreads < _Int, _Flt >::Ilu2BlockIlu2 (_b_blk_wells, NULL, iparam1_out,
                                                   ordlevel, ordtype, strtype, fcttype,
                                                   fcttype, pivmin, tau1, tau2, tau2_sch,
                                                   theta, _nblks, _nblks1, _nblks2, _blks,
                                                   _nzord_blks, _tree0, _tree1, alu_pair,
                                                   _tree2_new, _nblks_new, _blks_new,
                                                   _ordernew, lu_pair, _nmodif,
                                                   _eigmin_att, _eigmax_att);

      long long *p_blks_new = &_blks_new[0];
      int *p_ordernew = _ordernew.Ptr ();

//      PrintArray (ffout, " Order ", ntot, p_ordernew);

// Split pair into L and U parts

//   ffout << " Ilu2HMatrix: Point 4" << endl;

        CFctThreads < _Int, _Flt >::SplitPair (_nblks_new, p_blks_new, lu_pair, _l_matr,
                                               _u_matr);

//      ffout << " Ilu2HMatrix: LU after split: " << endl;
//      ffout << "    L = " << endl;
//      if (ntot < 5000) _l_matr.PrintHMatrix (ffout);
//      ffout << "    U = " << endl;
//      if (ntot < 5000) _u_matr.PrintHMatrix (ffout);

// Reorder scaling and compute its inverse

        CVectorData < _Flt > inv_scl_L (ntot);
        CVectorData < _Flt > inv_scl_U (ntot);

      _Flt *pinv_scl_L = inv_scl_L.Ptr ();
      _Flt *pinv_scl_U = inv_scl_U.Ptr ();

        CVector < _Flt >::OrderVector_thr (ntot, p_ordernew, pscl_L, pinv_scl_L);
        CVector < _Flt >::OrderVector_thr (ntot, p_ordernew, pscl_U, pinv_scl_U);

        CVector < _Flt >::CopyVector_thr (ntot, pinv_scl_L, pscl_L);
        CVector < _Flt >::CopyVector_thr (ntot, pinv_scl_U, pscl_U);

        CVector < _Flt >::InverseVector_thr (ntot, pinv_scl_L);
        CVector < _Flt >::InverseVector_thr (ntot, pinv_scl_U);

//      PrintArray (ffout, " Order ", ntot, p_ordernew);

//      PrintArray (ffout, " Ordered SclL ", ntot, pscl_L);
//      PrintArray (ffout, " Ordered SclU ", ntot, pscl_U);

//      PrintArray (ffout, " Ordered InvSclL ", ntot, pinv_scl_L);
//      PrintArray (ffout, " Ordered InvSclU ", ntot, pinv_scl_U);

// Rescale factors

        CFctThreads < _Int, _Flt >::RescaleU (_nblks_new, p_blks_new, _l_matr, pscl_L,
                                              pinv_scl_L);
        CFctThreads < _Int, _Flt >::RescaleU (_nblks_new, p_blks_new, _u_matr, pscl_U,
                                              pinv_scl_U);

//      ffout << " Ilu2HMatrix: End Point " << endl;
//      ffout << "    L = " << endl;
//      if (ntot < 5000) _l_matr.PrintHMatrix (ffout);
//      ffout << "    U = " << endl;
//      if (ntot < 5000) _u_matr.PrintHMatrix (ffout);

   }

//
// Compute explicit scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CFctThreads < _Int,
      _Flt >::ComputeExplicitScaling (int _sctype, int _nitersc, double _sclmin,
                                      int _nblks, long long *_blks, CBMatrix < _Int,
                                      _Flt > &_a_matr, CVectorData < _Flt > &_scl_L,
                                      CVectorData < _Flt > &_scl_U, double &_sclmin_att,
                                      double &_sclmax_att)
   {

// Allocate scaling

//      ofstream ffout ("ChkScl_1x1.dat");

//      _a_matr.PrintHMatrix (ffout);

      int ntot = (int) _blks[_nblks];

      _scl_L.resize (ntot);
      _scl_U.resize (ntot);

      _Flt *p_sclL = _scl_L.Ptr ();
      _Flt *p_sclU = _scl_U.Ptr ();

// Open hmatrix

      CMatrix < int, float >*phmatr = _a_matr.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub = _a_matr.GetASubArr ();

      int nzja_hmatr = phmatr->GetNzja ();
      int *pia_hmatr = phmatr->GetIaArr ();
      int *pja_hmatr = phmatr->GetJaArr ();

// Simple diagonal based scaling

      _sclmin_att = 1.0e100;
      _sclmax_att = -1.0e100;

      if (_sctype == -1) {

         CVector < _Flt >::SetByOnes_thr (ntot, p_sclL);
         CVector < _Flt >::SetByOnes_thr (ntot, p_sclU);

         _Flt diag;

         diag = (_Flt) 1.;

         _sclmin_att = diag;
         _sclmax_att = diag;

      } else if (_sctype == 0) {

         int n_thr = 1;

#ifdef USE_THREADS
         n_thr = omp_get_max_threads ();
#endif

         vector < double >sclmin_thr (n_thr + 1);
         vector < double >sclmax_thr (n_thr + 1);

         double *psclmin_thr = &sclmin_thr[0];
         double *psclmax_thr = &sclmax_thr[0];

         int jjj;

         for (jjj = 0; jjj < n_thr; jjj++)
            psclmin_thr[jjj] = 1.0e100;
         for (jjj = 0; jjj < n_thr; jjj++)
            psclmax_thr[jjj] = -1.0e100;

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < _nblks; ipar++) {
            int my_thr = 0;
#ifdef USE_THREADS
            my_thr = omp_get_thread_num ();
#endif
            int niloc = (int) (_blks[ipar + 1] - _blks[ipar]);
            int ibeg = (int) _blks[ipar];

            CVector < _Flt >::SetByOnes (niloc, p_sclL + ibeg);
            CVector < _Flt >::SetByOnes (niloc, p_sclU + ibeg);

            int j, jj, ki, irow, kj, kjj;
            _Flt diag;

            for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
               jj = pja_hmatr[j];
               if (jj == ipar) {
                  int nlist_temp = pA_sub[j].GetNlist ();
                  _Int *plist_temp = pA_sub[j].GetListArr ();
                  _Int *pia_temp = pA_sub[j].GetIaArr ();
                  _Int *pja_temp = pA_sub[j].GetJaArr ();
                  _Flt *pa_temp = pA_sub[j].GetAArr ();
                  for (ki = 0; ki < nlist_temp; ki++) {
                     irow = (int) plist_temp[ki];
                     diag = 1.;
                     for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                        kjj = (int) pja_temp[kj];
                        if (kjj == irow) {
                           diag = pa_temp[kj];
                        }
                     }
                     if (diag >= 0.) {
                        if (diag < psclmin_thr[my_thr])
                           psclmin_thr[my_thr] = diag;
                        if (diag > psclmax_thr[my_thr])
                           psclmax_thr[my_thr] = diag;
                        if (diag < _sclmin)
                           diag = (_Flt) _sclmin;
                        diag = (_Flt) (1. / sqrt (diag));
                        p_sclL[ibeg + irow] = diag;
                        p_sclU[ibeg + irow] = diag;
                     } else {
                        diag = -diag;
                        if (diag < psclmin_thr[my_thr])
                           psclmin_thr[my_thr] = diag;
                        if (diag > psclmax_thr[my_thr])
                           psclmax_thr[my_thr] = diag;
                        if (diag < _sclmin)
                           diag = (_Flt) _sclmin;
                        diag = (_Flt) (1. / sqrt (diag));
                        p_sclL[ibeg + irow] = -diag;
                        p_sclU[ibeg + irow] = diag;
                     }

                  }
               }
            }
         }

         for (jjj = 0; jjj < n_thr; jjj++) {
            if (_sclmin_att > psclmin_thr[jjj])
               _sclmin_att = psclmin_thr[jjj];
            if (_sclmax_att < psclmax_thr[jjj])
               _sclmax_att = psclmax_thr[jjj];
         }

      } else {

// Compute transposed block sparsity

         CMatrix < int, float >hmatrT_str;
         vector < int >indHt2H (nzja_hmatr + 1);

         int *pindHt2H = &indHt2H[0];

         {

            CVectorData < int >imask (5 * _nblks + 1);
            int *pimask = imask.Ptr ();

            int i;

            for (i = 0; i < _nblks; i++)
               pimask[i] = -1;

            int icycle = -1;

            phmatr->TransposedSparsityListSp (icycle, pimask, pimask + _nblks,
                                              pimask + 2 * _nblks, pimask + 3 * _nblks,
                                              pimask + 4 * _nblks, hmatrT_str);

            int *pia_hmatrT = hmatrT_str.GetIaArr ();
            int *pja_hmatrT = hmatrT_str.GetJaArr ();

            for (i = 0; i < _nblks; i++)
               pimask[i] = pia_hmatr[i];

            {

               int j, jj, k;

               for (i = 0; i < _nblks; i++) {
                  for (j = pia_hmatrT[i]; j < pia_hmatrT[i + 1]; j++) {
                     jj = pja_hmatrT[j];
                     k = pimask[jj];
                     pindHt2H[j] = k;
                     pimask[jj]++;
                  }
               }

            }

         }

         int *pia_hmatrT = hmatrT_str.GetIaArr ();
         int *pja_hmatrT = hmatrT_str.GetJaArr ();

// Iterative rows/columns balancing scaling

         double sclRmin, sclRmax;
         double sclCmin, sclCmax;

         sclRmin = 1.0e100;
         sclRmax = 0.0e0;
         sclCmin = 1.0e100;
         sclCmax = 0.0e0;

         int n_thr = 1;

#ifdef USE_THREADS
         n_thr = omp_get_max_threads ();
#endif

         vector < double >sclRmin_thr (n_thr + 1);
         vector < double >sclRmax_thr (n_thr + 1);
         vector < double >sclCmin_thr (n_thr + 1);
         vector < double >sclCmax_thr (n_thr + 1);

         double *psclRmin_thr = &sclRmin_thr[0];
         double *psclRmax_thr = &sclRmax_thr[0];
         double *psclCmin_thr = &sclCmin_thr[0];
         double *psclCmax_thr = &sclCmax_thr[0];

         int jjj;

         for (jjj = 0; jjj < n_thr; jjj++)
            psclRmin_thr[jjj] = 1.0e100;
         for (jjj = 0; jjj < n_thr; jjj++)
            psclRmax_thr[jjj] = 0.0e0;
         for (jjj = 0; jjj < n_thr; jjj++)
            psclCmin_thr[jjj] = 1.0e100;
         for (jjj = 0; jjj < n_thr; jjj++)
            psclCmax_thr[jjj] = 0.0e0;

         {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < _nblks; ipar++) {

               int j, ki, irow, kj;
               _Flt aux1;

               int niloc = (int) (_blks[ipar + 1] - _blks[ipar]);
               int ibeg = (int) _blks[ipar];

               CVector < _Flt >::SetByZeroes (niloc, p_sclL + ibeg);

               for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
                  int nlist_temp = pA_sub[j].GetNlist ();
                  _Int *plist_temp = pA_sub[j].GetListArr ();
                  _Int *pia_temp = pA_sub[j].GetIaArr ();
                  _Flt *pa_temp = pA_sub[j].GetAArr ();
                  for (ki = 0; ki < nlist_temp; ki++) {
                     irow = (int) plist_temp[ki];
                     for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                        aux1 = pa_temp[kj];
                        p_sclL[ibeg + irow] += aux1 * aux1;
                     }
                  }
               }
            }
         }

         {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < _nblks; ipar++) {

               int my_thr = 0;
#ifdef USE_THREADS
               my_thr = omp_get_thread_num ();
#endif

               int j;
               _Flt aux;

               int niloc = (int) (_blks[ipar + 1] - _blks[ipar]);
               int ibeg = (int) _blks[ipar];

               for (j = 0; j < niloc; j++) {
                  aux = p_sclL[ibeg + j];
                  if (aux < psclRmin_thr[my_thr])
                     psclRmin_thr[my_thr] = aux;
                  if (aux > psclRmax_thr[my_thr])
                     psclRmax_thr[my_thr] = aux;
                  aux = (_Flt) (1.0 / aux);
                  aux = sqrt (aux);
                  p_sclL[ibeg + j] = aux;
               }
            }
         }

         {

            int j;

            for (j = 0; j < n_thr; j++) {
               if (sclRmin > psclRmin_thr[j])
                  sclRmin = psclRmin_thr[j];
               if (sclRmax < psclRmax_thr[j])
                  sclRmax = psclRmax_thr[j];
            }

         }

         sclRmin = sqrt (sclRmin);
         sclRmax = sqrt (sclRmax);

         CVector < _Flt >::SetByZeroes_thr (ntot, p_sclU);

         {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < _nblks; ipar++) {

               int ind, kk;
               int j, ki, kj;
               _Flt aux1;

               int ibeg = (int) _blks[ipar];

               for (j = pia_hmatrT[ipar]; j < pia_hmatrT[ipar + 1]; j++) {
                  ind = pindHt2H[j];
                  int nlist_temp = pA_sub[ind].GetNlist ();
                  _Int *pia_temp = pA_sub[ind].GetIaArr ();
                  _Int *pja_temp = pA_sub[ind].GetJaArr ();
                  _Flt *pa_temp = pA_sub[ind].GetAArr ();
                  for (ki = 0; ki < nlist_temp; ki++) {
                     for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                        kk = (int) pja_temp[kj];
                        aux1 = pa_temp[kj];
                        p_sclU[ibeg + kk] += aux1 * aux1;
                     }
                  }
               }
            }
         }

         {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < _nblks; ipar++) {

               int my_thr = 0;
#ifdef USE_THREADS
               my_thr = omp_get_thread_num ();
#endif
               int j;
               _Flt aux;

               int niloc = (int) (_blks[ipar + 1] - _blks[ipar]);
               int ibeg = (int) _blks[ipar];

               for (j = 0; j < niloc; j++) {
                  aux = p_sclU[ibeg + j];
                  if (aux < psclCmin_thr[my_thr])
                     psclCmin_thr[my_thr] = aux;
                  if (aux > psclCmax_thr[my_thr])
                     psclCmax_thr[my_thr] = aux;
               }
            }
         }

         {

            int j;

            for (j = 0; j < n_thr; j++) {
               if (sclCmin > psclCmin_thr[j])
                  sclCmin = psclCmin_thr[j];
               if (sclCmax < psclCmax_thr[j])
                  sclCmax = psclCmax_thr[j];
            }

         }

         sclCmin = sqrt (sclCmin);
         sclCmax = sqrt (sclCmax);

         _sclmin_att = sclRmin;
         if (sclCmin < _sclmin_att)
            _sclmin_att = sclCmin;
         _sclmax_att = sclRmax;
         if (sclCmax > _sclmax_att)
            _sclmax_att = sclCmax;

         int iter;

         for (iter = 0; iter < _nitersc; iter++) {

            CVector < _Flt >::SetByZeroes_thr (ntot, p_sclU);

            {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
               for (int ipar = 0; ipar < _nblks; ipar++) {

                  int ind, kk;
                  int j, jj, ki, irow, kj;
                  _Flt aux1;

                  int ibeg = (int) _blks[ipar];
                  for (j = pia_hmatrT[ipar]; j < pia_hmatrT[ipar + 1]; j++) {
                     jj = pja_hmatrT[j];
                     ind = pindHt2H[j];
                     int nlist_temp = pA_sub[ind].GetNlist ();
                     _Int *plist_temp = pA_sub[ind].GetListArr ();
                     _Int *pia_temp = pA_sub[ind].GetIaArr ();
                     _Int *pja_temp = pA_sub[ind].GetJaArr ();
                     _Flt *pa_temp = pA_sub[ind].GetAArr ();
                     int jbeg = (int) _blks[jj];
                     for (ki = 0; ki < nlist_temp; ki++) {
                        irow = (int) plist_temp[ki];
                        for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                           kk = (int) pja_temp[kj];
                           aux1 = pa_temp[kj];
                           p_sclU[ibeg + kk] += aux1 * aux1 * p_sclL[jbeg + irow];
                        }
                     }
                  }
               }
            }

            {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
               for (int ipar = 0; ipar < _nblks; ipar++) {

                  int j;
                  _Flt aux;

                  int niloc = (int) (_blks[ipar + 1] - _blks[ipar]);
                  int ibeg = (int) _blks[ipar];

                  for (j = 0; j < niloc; j++) {
                     aux = p_sclU[ibeg + j];
                     aux = (_Flt) (1.0 / aux);
                     p_sclU[ibeg + j] = aux;
                  }
               }

            }

            CVector < _Flt >::SetByZeroes_thr (ntot, p_sclL);

            {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
               for (int ipar = 0; ipar < _nblks; ipar++) {

                  int kk;
                  int j, jj, ki, irow, kj;
                  _Flt aux1;

                  int niloc = (int) (_blks[ipar + 1] - _blks[ipar]);
                  int ibeg = (int) _blks[ipar];
                  CVector < _Flt >::SetByZeroes (niloc, p_sclL + ibeg);
                  for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
                     jj = pja_hmatr[j];
                     int jbeg = (int) _blks[jj];
                     int nlist_temp = pA_sub[j].GetNlist ();
                     _Int *plist_temp = pA_sub[j].GetListArr ();
                     _Int *pia_temp = pA_sub[j].GetIaArr ();
                     _Int *pja_temp = pA_sub[j].GetJaArr ();
                     _Flt *pa_temp = pA_sub[j].GetAArr ();
                     for (ki = 0; ki < nlist_temp; ki++) {
                        irow = (int) plist_temp[ki];
                        for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                           kk = (int) pja_temp[kj];
                           aux1 = pa_temp[kj];
                           p_sclL[ibeg + irow] += aux1 * aux1 * p_sclU[jbeg + kk];
                        }
                     }
                  }
               }

            }

            {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
               for (int ipar = 0; ipar < _nblks; ipar++) {

                  int j;
                  _Flt aux;

                  int niloc = (int) (_blks[ipar + 1] - _blks[ipar]);
                  int ibeg = (int) _blks[ipar];

                  for (j = 0; j < niloc; j++) {
                     aux = p_sclL[ibeg + j];
                     aux = (_Flt) (1.0 / aux);
                     p_sclL[ibeg + j] = aux;
                  }
               }
            }

         }

         {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < _nblks; ipar++) {

               int j;

               int niloc = (int) (_blks[ipar + 1] - _blks[ipar]);
               int ibeg = (int) _blks[ipar];

               for (j = 0; j < niloc; j++) {
                  p_sclL[ibeg + j] = (_Flt) (sqrt (p_sclL[ibeg + j]));;
                  p_sclU[ibeg + j] = (_Flt) (sqrt (p_sclU[ibeg + j]));;
               }
            }
         }

      }

//      PrintArray (ffout, " Scl_L",ntot,p_sclL);
//      PrintArray (ffout, " Scl_U",ntot,p_sclU);

   }

//
// Perform explicit scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CFctThreads < _Int,
      _Flt >::HMatrixScale (int _nblks, long long *_blks, CBMatrix < _Int,
                            _Flt > &_a_matr, _Flt * _sclL, _Flt * _sclU)
   {

// Open hmatrix

//      ofstream ffout ("ChkScale_1x1.dat");

      CMatrix < int, float >*phmatr = _a_matr.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub = _a_matr.GetASubArr ();

      int *pia_hmatr = phmatr->GetIaArr ();
      int *pja_hmatr = phmatr->GetJaArr ();

#ifdef USE_THREADS
#pragma omp parallel for
#endif
      for (int ipar = 0; ipar < _nblks; ipar++) {

         int j, jj, ki, irow, kj, kk;
         _Flt aux;

         int ibeg = (int) _blks[ipar];

         for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
            jj = pja_hmatr[j];
            int jbeg = (int) _blks[jj];
            int nlist_temp = pA_sub[j].GetNlist ();
            _Int *plist_temp = pA_sub[j].GetListArr ();
            _Int *pia_temp = pA_sub[j].GetIaArr ();
            _Int *pja_temp = pA_sub[j].GetJaArr ();
            _Flt *pa_temp = pA_sub[j].GetAArr ();
            for (ki = 0; ki < nlist_temp; ki++) {
               irow = (int) plist_temp[ki];
               for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                  kk = (int) pja_temp[kj];
                  aux = pa_temp[kj];
                  pa_temp[kj] = aux * _sclL[ibeg + irow] * _sclU[jbeg + kk];
               }
            }
         }
      }

//      _a_matr.PrintHMatrix (ffout);

   }

//
// Split LU and combine pairs
//========================================================================================
   template < typename _Int, typename _Flt > void CFctThreads < _Int,
      _Flt >::SplitLUPair (int _nblks, long long *_blks, CBMatrix < _Int, _Flt > &_a_matr,
                           CBMatrix < _Int, _Flt > &_a_pair)
   {

// Open hmatrix

      CMatrix < int, float >*phmatr = _a_matr.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub = _a_matr.GetASubArr ();

      int nzja_hmatr = phmatr->GetNzja ();
      int *pia_hmatr = phmatr->GetIaArr ();
      int *pja_hmatr = phmatr->GetJaArr ();

// Compute transposed block sparsity

      CMatrix < int, float >hmatrT_str;
      vector < int >indH2Ht (nzja_hmatr + 1);

      int *pindH2Ht = &indH2Ht[0];

      {

         CVectorData < int >imask (5 * _nblks + 1);
         int *pimask = imask.Ptr ();

         int i;

         for (i = 0; i < _nblks; i++)
            pimask[i] = -1;

         int icycle = -1;

         phmatr->TransposedSparsityListSp (icycle, pimask, pimask + _nblks,
                                           pimask + 2 * _nblks, pimask + 3 * _nblks,
                                           pimask + 4 * _nblks, hmatrT_str);

         int *pia_hmatrT = hmatrT_str.GetIaArr ();
         int *pja_hmatrT = hmatrT_str.GetJaArr ();

         for (i = 0; i < _nblks; i++)
            pimask[i] = pia_hmatr[i];

         int j, jj, k;

         for (i = 0; i < _nblks; i++) {
            for (j = pia_hmatrT[i]; j < pia_hmatrT[i + 1]; j++) {
               jj = pja_hmatrT[j];
               k = pimask[jj];
               pindH2Ht[k] = j;
               pimask[jj]++;
            }
         }
      }

// Compute block sparsity of the triangular part of LU

      int nzja_pair = 0;

      int i, j, jj;

      for (i = 0; i < _nblks; i++) {
         for (j = pia_hmatr[i]; j < pia_hmatr[i + 1]; j++) {
            jj = pja_hmatr[j];
            if (jj >= i)
               nzja_pair++;
         }
      }

      vector < int >indpair2matr (nzja_pair + 1);
      int *pindpair2matr = &indpair2matr[0];

      _a_pair.SetNzblk (nzja_pair);
      _a_pair.ResizeASub (nzja_pair);

      CMatrix < int, float >*phmatr_pair = _a_pair.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_pair = _a_pair.GetASubArr ();

      phmatr_pair->ResizeAndSetAllSp (_nblks, _nblks, nzja_pair, nzja_pair);

      int *plist_pair = phmatr_pair->GetListArr ();
      int *plist2_pair = phmatr_pair->GetList2Arr ();
      int *pia_pair = phmatr_pair->GetIaArr ();
      int *pja_pair = phmatr_pair->GetJaArr ();
      int *pja2_pair = phmatr_pair->GetJa2Arr ();

      for (i = 0; i < _nblks; i++)
         plist_pair[i] = i;
      for (i = 0; i < _nblks; i++)
         plist2_pair[i] = 0;
      for (i = 0; i < nzja_pair; i++)
         pja2_pair[i] = 0;

      nzja_pair = 0;
      pia_pair[0] = 0;

      for (i = 0; i < _nblks; i++) {
         for (j = pia_hmatr[i]; j < pia_hmatr[i + 1]; j++) {
            jj = pja_hmatr[j];
            if (jj >= i) {
               pja_pair[nzja_pair] = jj;
               pindpair2matr[nzja_pair] = j;
               nzja_pair++;
            }
         }
         pia_pair[i + 1] = nzja_pair;
      }

// Compute maximal block size

      int nimax = 0;

      int niloc;

      for (i = 0; i < _nblks; i++) {
         niloc = (int) (_blks[i + 1] - _blks[i]);
         if (niloc > nimax)
            nimax = niloc;
      }

// Combine blocks

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      vector < int >icycle_thr (n_thr + 1);
      vector < CVectorData < int > >imaskblk_thr (n_thr + 1);

      int *picycle_thr = &icycle_thr[0];
      CVectorData < int >*pimaskblk_thr = &imaskblk_thr[0];

      for (i = 0; i < n_thr; i++)
         picycle_thr[i] = -1;

#ifdef USE_THREADS
#pragma omp parallel for
#endif
      for (int ipar = 0; ipar < _nblks; ipar++) {

         int my_thr = 0;
#ifdef USE_THREADS
         my_thr = omp_get_thread_num ();
#endif

         int j, jj;
         int ind1, ind2, k;
         vector < char >char_dummy;

         if (picycle_thr[my_thr] == -1) {
            pimaskblk_thr[my_thr].resize (nimax * 5 + 1);
            int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
            for (j = 0; j < nimax; j++)
               pimaskblk[j] = -1;
         }
         int icycleblk = picycle_thr[my_thr];
         int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
         for (j = pia_pair[ipar]; j < pia_pair[ipar + 1]; j++) {
            jj = pja_pair[j];
            if (jj == ipar) {

               ind2 = pindpair2matr[j];

               int nlist1 = pA_sub[ind2].GetNlist ();

               vector < _Int > *pia_lu = pA_sub[ind2].GetIa ();
               vector < _Int > *pja_lu = pA_sub[ind2].GetJa ();
               vector < _Flt > *pa_lu = pA_sub[ind2].GetA ();

               vector < _Int > ia_l;
               vector < _Int > ja_l;
               vector < _Flt > a_l;
               vector < _Int > ia_u;
               vector < _Int > ja_u;
               vector < _Flt > a_u;

               CFct_impl < _Int, _Flt >::SplitLU (false, nlist1, *pia_lu, *pja_lu,
                                                  char_dummy, *pa_lu, ia_l, ja_l,
                                                  char_dummy, a_l, ia_u, ja_u, char_dummy,
                                                  a_u);

               int nzja_u = (int) ia_u[nlist1];

               vector < _Int > ia_lt;
               vector < _Int > ja_lt;
               vector < _Flt > a_lt;

               CFct_impl < _Int, _Flt >::Transpose (false, nlist1, ia_l, ja_l, char_dummy,
                                                    a_l, ia_lt, ja_lt, char_dummy, a_lt);

               {
                  vector < _Int > ia_temp;
                  vector < _Int > ja_temp;
                  vector < _Flt > a_temp;
                  ia_l.swap (ia_temp);
                  ja_l.swap (ja_temp);
                  a_l.swap (a_temp);
               }

               CMatrix < _Int, _Flt > alu_pair;

               alu_pair.ResizeList (nlist1);

               _Int *plist_pair = alu_pair.GetListArr ();
               vector < _Int > *pia_pair = alu_pair.GetIa ();
               vector < _Int > *pja_pair = alu_pair.GetJa ();
               vector < _Flt > *pa_pair = alu_pair.GetA ();

               for (k = 0; k < nlist1; k++)
                  plist_pair[k] = k;

               CFct_impl < _Int, _Flt >::CombinePairs (false, nlist1, ia_lt, ja_lt,
                                                       char_dummy, a_lt, ia_u, ja_u,
                                                       char_dummy, a_u, *pia_pair,
                                                       *pja_pair, char_dummy, *pa_pair);

               alu_pair.SetNlist (nlist1);
               alu_pair.SetNzja (nzja_u);
               alu_pair.SetNza (2 * nzja_u);

               pA_sub_pair[j].ReplaceFree (alu_pair);

            } else {

               ind2 = pindpair2matr[j];
               ind1 = pindH2Ht[ind2];

               CMatrix < _Int, _Flt > al_blk;

               pA_sub[ind1].TransposedSparsityList (icycleblk, pimaskblk,
                                                    pimaskblk + nimax,
                                                    pimaskblk + 2 * nimax,
                                                    pimaskblk + 3 * nimax,
                                                    pimaskblk + 4 * nimax, al_blk);

               int nlist1 = al_blk.GetNlist ();
               int nzja1 = al_blk.GetNzja ();
               _Int *plist1 = al_blk.GetListArr ();
               vector < _Int > *pia1 = al_blk.GetIa ();
               vector < _Int > *pja1 = al_blk.GetJa ();
               vector < _Flt > *pa1 = al_blk.GetA ();

               int nlist2 = pA_sub[ind2].GetNlist ();
               int nzja2 = pA_sub[ind2].GetNzja ();
               vector < _Int > *pia2 = pA_sub[ind2].GetIa ();
               vector < _Int > *pja2 = pA_sub[ind2].GetJa ();
               vector < _Flt > *pa2 = pA_sub[ind2].GetA ();

               if (nlist1 != nlist2 || nzja1 != nzja2) {
                  throw " CFctThreads<_Int,_Flt>::SplitLUPair: error in nlist or nzja !";
               }

               CMatrix < _Int, _Flt > alu_blk;

               alu_blk.ResizeList (nlist1);

               _Int *plist_lu = alu_blk.GetListArr ();
               vector < _Int > *pia_lu = alu_blk.GetIa ();
               vector < _Int > *pja_lu = alu_blk.GetJa ();
               vector < _Flt > *pa_lu = alu_blk.GetA ();

               for (k = 0; k < nlist1; k++)
                  plist_lu[k] = plist1[k];

               CFct_impl < _Int, _Flt >::CombinePairs (false, nlist1, *pia1, *pja1,
                                                       char_dummy, *pa1, *pia2, *pja2,
                                                       char_dummy, *pa2, *pia_lu, *pja_lu,
                                                       char_dummy, *pa_lu);

               alu_blk.SetNlist (nlist1);
               alu_blk.SetNzja (nzja1);
               alu_blk.SetNza (2 * nzja1);

               pA_sub_pair[j].ReplaceFree (alu_blk);

            }
         }

         picycle_thr[my_thr] = icycleblk;
      }

   }

//
// Split pair into L and U
//========================================================================================
   template < typename _Int, typename _Flt > void CFctThreads < _Int,
      _Flt >::SplitPair (int _nblks, long long *_blks, CBMatrix < _Int, _Flt > &_a_pair,
                         CBMatrix < _Int, _Flt > &_a_L, CBMatrix < _Int, _Flt > &_a_U)
   {

// Open hmatrix

      CMatrix < int, float >*phmatr = _a_pair.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub = _a_pair.GetASubArr ();

      int nzblk = _a_pair.GetNzblk ();

      int *pia_hmatr = phmatr->GetIaArr ();

// Create block sparsity of L and U

      _a_L.ResizeASub (nzblk);
      _a_L.SetNzblk (nzblk);
      _a_U.ResizeASub (nzblk);
      _a_U.SetNzblk (nzblk);

      CMatrix < int, float >*phmatr_L = _a_L.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_L = _a_L.GetASubArr ();

      CMatrix < int, float >*phmatr_U = _a_U.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_U = _a_U.GetASubArr ();

      *phmatr_L = *phmatr;
      *phmatr_U = *phmatr;

// Split all blocks

#ifdef USE_THREADS
#pragma omp parallel for
#endif
      for (int ipar = 0; ipar < _nblks; ipar++) {

         vector < char >char_dummy;
         double tau1 = -1.0e0;

         int j, k;

         for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {

            int nlist_lu = pA_sub[j].GetNlist ();
            _Int *plist_lu = pA_sub[j].GetListArr ();
            vector < _Int > *pia_lu = pA_sub[j].GetIa ();
            vector < _Int > *pja_lu = pA_sub[j].GetJa ();
            vector < _Flt > *pa_lu = pA_sub[j].GetA ();

            pA_sub_L[j].ResizeList (nlist_lu);
            pA_sub_U[j].ResizeList (nlist_lu);

            _Int *plist_l = pA_sub_L[j].GetListArr ();
            vector < _Int > *pia_l = pA_sub_L[j].GetIa ();
            vector < _Int > *pja_l = pA_sub_L[j].GetJa ();
            vector < _Flt > *pa_l = pA_sub_L[j].GetA ();

            _Int *plist_u = pA_sub_U[j].GetListArr ();
            vector < _Int > *pia_u = pA_sub_U[j].GetIa ();
            vector < _Int > *pja_u = pA_sub_U[j].GetJa ();
            vector < _Flt > *pa_u = pA_sub_U[j].GetA ();

            for (k = 0; k < nlist_lu; k++)
               plist_l[k] = plist_lu[k];
            for (k = 0; k < nlist_lu; k++)
               plist_u[k] = plist_lu[k];

            CFct_impl < _Int, _Flt >::SplitPairsFilter (false, tau1, nlist_lu, *pia_lu,
                                                        *pja_lu, char_dummy, *pa_lu,
                                                        *pia_l, *pja_l, char_dummy, *pa_l,
                                                        *pia_u, *pja_u, char_dummy,
                                                        *pa_u);

            int nzja_L = (int) ((*pia_l)[nlist_lu]);
            int nzja_U = (int) ((*pia_u)[nlist_lu]);

            pA_sub_L[j].SetNlist (nlist_lu);
            pA_sub_L[j].SetNzja (nzja_L);
            pA_sub_L[j].SetNza (nzja_L);

            pA_sub_U[j].SetNlist (nlist_lu);
            pA_sub_U[j].SetNzja (nzja_U);
            pA_sub_U[j].SetNza (nzja_U);

         }
      }

   }

//========================================================================================
   template < typename _Int, typename _Flt > struct CBlockFctTree       /// Template structure that computes fct
   {
// Data
      int myid_cpu;             ///< Processor number
      int index2_out;           ///< Secondary index number
      int ordlevel;             ///< Ordering filtration level
      int ordtype;              ///< Ordering type
      int strtype;              ///< Structural control parameter
      int fcttype;              ///< Structural fct type
      double pivmin;            ///< Minimal pivot
      double tau1;              ///< First order filtering parameter
      double tau2;              ///< Second order filtering parameter
      double tau2_sch;          ///< Second order filtering parameter for Schur data
      double theta;             ///< Diagonal correction scaling
      int nblks;                ///< Total number of blocks
      int nblks1;               ///< Number of blocks in first  tree
      int nblks2;               ///< Number of blocks in second tree
      long long *pblks;         ///< Block partitioning
      int iblk;                 ///< Global block number
      CTree *ptree_1;           ///< First tree
      CTree *ptree_2;           ///< Second tree
      int *pia_ALU;             ///< Ia sparsity of ALU
      int *pja_ALU;             ///< Ja sparsity of ALU
        CMatrix < _Int, _Flt > *pASub_ALU;      ///< Blocks of ALU
        vector < int >*porder_arr;      ///< Set of ordering arrays for each block
        CBMatrix < _Int, _Flt > *pblockrowsLU_arr;      ///< Computed LU
        CBMatrix < _Int, _Flt > *pschur_arr;    ///< Schur data
      int nimax;                ///< Maximal block size
      int *picycle_thr;         ///< Array of cycle variables for threads
        CVectorData < int >*pimask_thr; ///< Mask arrays for threads
      int *pnmodif_thr;         ///< The number of modifs for threads
      double *peigmin_thr;      ///< Minimal pivots for threads
      double *peigmax_thr;      ///< Maximal pivots for threads
// Functions
      void Execute ()
      {

//         char strbuff[256];
//         if (index2_out == 0) {
//            sprintf (strbuff,"ChkFctK_1x1_%i_%i.dat",myid_cpu,iblk);
//         } else {
//            sprintf (strbuff,"ChkFctK_1x1_%i_Sch_%i.dat",myid_cpu,iblk);
//         }
//         ofstream ffout (strbuff);

//         cout << " Ihblk = " << myid_cpu << " Iblk = " << iblk << endl;
//         ffout << " Ihblk = " << myid_cpu << " Iblk = " << iblk << endl;

// Init mask arrays

         int my_thr = 0;

#ifdef USE_THREADS
           my_thr = omp_get_thread_num ();
#endif

         if (picycle_thr[my_thr] == -1) {
            pimask_thr[my_thr].resize (3 * nblks + 3 * nimax);
            int *pimaskblk = pimask_thr[my_thr].Ptr ();
            int *pimask = pimaskblk + nblks;
            int i;
            for (i = 0; i < nblks; i++)
                 pimaskblk[i] = -1;
            for (i = 0; i < nimax; i++)
                 pimask[i] = -1;
         }

         int icycleblk = picycle_thr[my_thr];
         int *pimaskblk = pimask_thr[my_thr].Ptr ();
         int *pimask = pimaskblk + nblks;
         int *pindblk = pimask + nimax;
         int *plistblk = pindblk + nblks;
         int *plist = plistblk + nblks;

// Determine list of childs of a node if any

         int ichild1 = -1;
         int ichild2 = -1;

         if (iblk < nblks1) {
            int *pnchilds = ptree_1->GetNchilds ();
            vector < int >*pchilds_list = ptree_1->GetChildsList ();
            int *ppchilds = &pchilds_list[iblk][0];
            if (pnchilds[iblk] == 2) {
               ichild1 = ppchilds[0];
               ichild2 = ppchilds[1];
            } else if (pnchilds[iblk] == 1) {
               ichild1 = ppchilds[0];
               if (ichild1 == iblk)
                  ichild1 = -1;
            } else {
               throw " CBlockFctTree<>::Execute: wrong tree !!! ";
            }
         } else {
            int iblk1 = iblk - nblks1;
            int *pnchilds = ptree_2->GetNchilds ();
            vector < int >*pchilds_list = ptree_2->GetChildsList ();
            int *ppchilds = &pchilds_list[iblk1][0];
            if (pnchilds[iblk1] == 2) {
               ichild1 = ppchilds[0] + nblks1;
               ichild2 = ppchilds[1] + nblks1;
            } else if (pnchilds[iblk1] == 1) {
               ichild1 = ppchilds[0];
               if (ichild1 == iblk1) {
                  ichild1 = -1;
               } else {
                  ichild1 += nblks1;
               }
            } else {
               throw " CBlockFctTree<>::Execute: wrong tree !!! ";
            }
         }

// Add Schur hblocks

         CBMatrix < _Int, _Flt > hblk_sum;

         if (ichild1 >= 0 && ichild2 >= 0) {
//            ffout << " Hblk Ichild 1: " << endl;
//            pschur_arr[ichild1].PrintHMatrix (ffout);
//            ffout << " Hblk Ichild 2: " << endl;
//            pschur_arr[ichild2].PrintHMatrix (ffout);
            hblk_sum.ReplaceFree (pschur_arr[ichild1]);
            hblk_sum %= pschur_arr[ichild2];
            pschur_arr[ichild1].Clean ();
            pschur_arr[ichild2].Clean ();
         } else if (ichild1 >= 0) {
            hblk_sum.ReplaceFree (pschur_arr[ichild1]);
            pschur_arr[ichild1].Clean ();
         } else if (ichild2 >= 0) {
            hblk_sum.ReplaceFree (pschur_arr[ichild2]);
            pschur_arr[ichild2].Clean ();
         }
//         ffout << " =========== Hblk_Schur_sum: " << endl;
//         hblk_sum.PrintHMatrix (ffout);

// Get set of blocks that correspond to current block row and split hblock for future add

         int nzja_curr = 0;
         vector < int >ja_curr;
         vector < CMatrix < _Int, _Flt > >ASub_curr;

         {

            int nzja_curr0 = pia_ALU[iblk + 1] - pia_ALU[iblk];
            int *pja_curr0 = pja_ALU + pia_ALU[iblk];
            CMatrix < _Int, _Flt > *pA_sub_curr0 = pASub_ALU + pia_ALU[iblk];

//            PrintArray (ffout," ja_curr0 = ",nzja_curr0,pja_curr0);

// Add zero jachar data if not available

            int i, j;

            for (i = 0; i < nzja_curr0; i++) {
               int nzja_temp = pA_sub_curr0[i].GetNzja ();
               int nzjachar_temp = pA_sub_curr0[i].GetNzjaChar ();
               if (nzja_temp != nzjachar_temp) {
                  pA_sub_curr0[i].ResizeJaChar (nzja_temp);
                  pA_sub_curr0[i].SetNzjaChar (nzja_temp);
                  char *pjachar_temp = pA_sub_curr0[i].GetJaCharArr ();
                  for (j = 0; j < nzja_temp; j++)
                     pjachar_temp[j] = 0;
               }
            }

            CMatrix < int, float >*phmatr_sum = hblk_sum.GetHMatrStr ();
            CMatrix < _Int, _Flt > *pA_sub_sum = hblk_sum.GetASubArr ();

            int nlist_hmatr_sum = phmatr_sum->GetNlist ();
            int nzja_hmatr_sum = phmatr_sum->GetNzja ();
            int *plist_hmatr_sum = phmatr_sum->GetListArr ();
            int *pia_hmatr_sum = phmatr_sum->GetIaArr ();
            int *pja_hmatr_sum = phmatr_sum->GetJaArr ();

            int ilist_curr = -1;

            for (i = 0; i < nlist_hmatr_sum; i++) {
               if (plist_hmatr_sum[i] == iblk)
                  ilist_curr = i;
            }

            int nzja_curr1 = 0;
            int *pja_curr1 = NULL;
            CMatrix < _Int, _Flt > *pA_sub_curr1 = NULL;

            if (ilist_curr >= 0) {

               nzja_curr1 = pia_hmatr_sum[ilist_curr + 1] - pia_hmatr_sum[ilist_curr];
               pja_curr1 = pja_hmatr_sum + pia_hmatr_sum[ilist_curr];
               pA_sub_curr1 = pA_sub_sum + pia_hmatr_sum[ilist_curr];

            }
//         PrintArray (ffout," ja_curr1 = ",nzja_curr1,pja_curr1);

            ja_curr.resize (nzja_curr0 + nzja_curr1 + 1);
            ASub_curr.resize (nzja_curr0 + nzja_curr1 + 1);

            int *pja_curr = &ja_curr[0];
            CMatrix < _Int, _Flt > *pASub_curr = &ASub_curr[0];

            nzja_curr = 0;

            int ip0 = 0;
            int ip1 = 0;

            int jj0, jj1;

            while (ip0 < nzja_curr0 || ip1 < nzja_curr1) {
               if (ip0 < nzja_curr0 && ip1 < nzja_curr1) {
                  jj0 = pja_curr0[ip0];
                  jj1 = pja_curr1[ip1];
                  if (jj0 == jj1) {
                     pja_curr[nzja_curr] = jj0;
                     CMatrix < _Int, _Flt > asum;
                     asum.ReplaceFree (pA_sub_curr0[ip0]);
                     asum %= pA_sub_curr1[ip1];
                     pASub_curr[nzja_curr].ReplaceFree (asum);
                     pA_sub_curr1[ip1].Clean ();
                     ip0++;
                     ip1++;
                     nzja_curr++;
                  } else if (jj0 < jj1) {
                     pja_curr[nzja_curr] = jj0;
                     pASub_curr[nzja_curr].ReplaceFree (pA_sub_curr0[ip0]);
                     ip0++;
                     nzja_curr++;
                  } else if (jj0 > jj1) {
                     pja_curr[nzja_curr] = jj1;
                     pASub_curr[nzja_curr].ReplaceFree (pA_sub_curr1[ip1]);
                     ip1++;
                     nzja_curr++;
                  }
               } else if (ip0 < nzja_curr0) {
                  pja_curr[nzja_curr] = pja_curr0[ip0];
                  pASub_curr[nzja_curr].ReplaceFree (pA_sub_curr0[ip0]);
                  ip0++;
                  nzja_curr++;
               } else if (ip1 < nzja_curr1) {
                  pja_curr[nzja_curr] = pja_curr1[ip1];
                  pASub_curr[nzja_curr].ReplaceFree (pA_sub_curr1[ip1]);
                  ip1++;
                  nzja_curr++;
               }
            }

            if (ilist_curr >= 0) {

               int nzja_flt = nzja_hmatr_sum - nzja_curr1;

               CBMatrix < _Int, _Flt > hblk_flt;

               hblk_flt.SetNzblk (nzja_flt);
               hblk_flt.ResizeASub (nzja_flt);

               CMatrix < int, float >*phmatr_flt = hblk_flt.GetHMatrStr ();
               CMatrix < _Int, _Flt > *pA_sub_flt = hblk_flt.GetASubArr ();

               int nlist_hmatr_flt = nlist_hmatr_sum - 1;

               phmatr_flt->ResizeAndSetAllSp (nlist_hmatr_flt, nlist_hmatr_flt, nzja_flt,
                                              nzja_flt);

               int *plist_flt = phmatr_flt->GetListArr ();
               int *plist2_flt = phmatr_flt->GetList2Arr ();
               int *pia_flt = phmatr_flt->GetIaArr ();
               int *pja_flt = phmatr_flt->GetJaArr ();
               int *pja2_flt = phmatr_flt->GetJa2Arr ();

               int i;

               for (i = 0; i < nlist_hmatr_flt; i++)
                  plist2_flt[i] = 0;
               for (i = 0; i < nzja_flt; i++)
                  pja2_flt[i] = 0;

               nlist_hmatr_flt = 0;
               nzja_flt = 0;

               int j;

               for (i = 0; i < nlist_hmatr_sum; i++) {
                  if (plist_hmatr_sum[i] != iblk) {
                     plist_flt[nlist_hmatr_flt] = plist_hmatr_sum[i];
                     for (j = pia_hmatr_sum[i]; j < pia_hmatr_sum[i + 1]; j++) {
                        pja_flt[nzja_flt] = pja_hmatr_sum[j];
                        pA_sub_flt[nzja_flt].ReplaceFree (pA_sub_sum[j]);
                        nzja_flt++;
                     }
                     pia_flt[nlist_hmatr_flt + 1] = nzja_flt;
                     nlist_hmatr_flt++;
                  }
               }

               hblk_sum.ReplaceFree (hblk_flt);

            }

         }

         int *pja_curr = &ja_curr[0];
         CMatrix < _Int, _Flt > *pASub_curr = &ASub_curr[0];

//         ffout << " Iblk = " << iblk << endl;
//         PrintArray (ffout," ja_curr = ",nzja_curr,pja_curr);

//         ffout << " =========== Hblk_Schur_sum flt: " << endl;
//         hblk_sum.PrintHMatrix (ffout);

// Compute optimal ordering via diagonal block

         int ni_diag = (int) (pblks[iblk + 1] - pblks[iblk]);

         int jblk0 = pja_curr[0];

//      ffout << " Point 1 check " << endl;
         if (jblk0 != iblk) {
            cout << " CBlockFctTree<>::Execute(): error in block sparsity !" << " Iblk = "
               << iblk << " jblk0 = " << jblk0 << " Ni_diag = " << ni_diag << endl;
            throw " CBlockFctTree<>::Execute(): error in block sparsity ! ";
         }
//      ffout << " Point 2 check " << endl;

         if (ordlevel < 0) {
            CFct < _Int, _Flt >::ComputeOptimalOrder (pASub_curr[0], ordtype,
                                                      porder_arr[iblk]);
//         porder_arr[iblk].resize (ni_diag+1);
//         int *pporder_temp = &porder_arr[iblk][0];
//         for (int iii=0;iii<ni_diag;iii++) pporder_temp[iii] = iii;
         } else {
//         ffout << " Point 3 check " << endl;
            int nzja_dia = pASub_curr[0].GetNzja ();
            _Int *pia_dia = pASub_curr[0].GetIaArr ();
            _Int *pja_dia = pASub_curr[0].GetJaArr ();
            char *pjachar_dia = pASub_curr[0].GetJaCharArr ();
            int nzja_flt = 0;
            int i;
            for (i = 0; i < nzja_dia; i++) {
               if (pjachar_dia[i] <= ordlevel)
                  nzja_flt++;
            }
            CMatrix < _Int, _Flt > dia_flt;
            dia_flt.ResizeAndSetAllSp (ni_diag, 0, nzja_flt, 0);
            _Int *plist_flt = dia_flt.GetListArr ();
            _Int *pia_flt = dia_flt.GetIaArr ();
            _Int *pja_flt = dia_flt.GetJaArr ();
            for (i = 0; i < ni_diag; i++)
               plist_flt[i] = (_Int) i;
            nzja_flt = 0;
            pia_flt[0] = 0;
            int j;
            for (i = 0; i < ni_diag; i++) {
               for (j = (int) pia_dia[i]; j < pia_dia[i + 1]; j++) {
                  if (pjachar_dia[j] <= ordlevel) {
                     pja_flt[nzja_flt] = (_Int) pja_dia[j];
                     nzja_flt++;
                  }
               }
               pia_flt[i + 1] = (_Int) nzja_flt;
            }
            CFct < _Int, _Flt >::ComputeOptimalOrder (dia_flt, ordtype, porder_arr[iblk]);
         }

//      ffout << " Point 4 check " << endl;
         int *pporder_curr = &porder_arr[iblk][0];
//         PrintArray (ffout," Order = ",ni_diag,pporder_curr);

// Reorder diagonal block

         {

// Split pair

            vector < _Int > *pia_alu_dia = pASub_curr[0].GetIa ();
            vector < _Int > *pja_alu_dia = pASub_curr[0].GetJa ();
            vector < char >*pjachar_alu_dia = pASub_curr[0].GetJaChar ();
            vector < _Flt > *pa_alu_dia = pASub_curr[0].GetA ();

            vector < _Int > ia_l_dia;
            vector < _Int > ja_l_dia;
            vector < char >jachar_l_dia;
            vector < _Flt > a_l_dia;

            vector < _Int > ia_u_dia;
            vector < _Int > ja_u_dia;
            vector < char >jachar_u_dia;
            vector < _Flt > a_u_dia;

            double tau1 = -1.0e0;

//      ffout << " Point 5 check " << endl;

            CFct_impl < _Int, _Flt >::SplitPairsFilter (true, tau1, ni_diag, *pia_alu_dia,
                                                        *pja_alu_dia, *pjachar_alu_dia,
                                                        *pa_alu_dia, ia_l_dia, ja_l_dia,
                                                        jachar_l_dia, a_l_dia, ia_u_dia,
                                                        ja_u_dia, jachar_u_dia, a_u_dia);

// Transpose L

            vector < _Int > ia_lt_dia;
            vector < _Int > ja_lt_dia;
            vector < char >jachar_lt_dia;
            vector < _Flt > a_lt_dia;

//      ffout << " Point 6 check " << endl;

            CFct_impl < _Int, _Flt >::Transpose (true, ni_diag, ia_l_dia, ja_l_dia,
                                                 jachar_l_dia, a_l_dia, ia_lt_dia,
                                                 ja_lt_dia, jachar_lt_dia, a_lt_dia);

            {
               vector < _Int > ia_dummy;
               vector < _Int > ja_dummy;
               vector < char >jachar_dummy;
               vector < _Flt > a_dummy;
               ia_l_dia.swap (ia_dummy);
               ja_l_dia.swap (ja_dummy);
               jachar_l_dia.swap (jachar_dummy);
               a_l_dia.swap (a_dummy);
            }

// Combine point L+U

            CMatrix < _Int, _Flt > amatr_lu;

            vector < _Int > *pia_lu = amatr_lu.GetIa ();
            vector < _Int > *pja_lu = amatr_lu.GetJa ();
            vector < char >*pjachar_lu = amatr_lu.GetJaChar ();
            vector < _Flt > *pa_lu = amatr_lu.GetA ();

//      ffout << " Point 6 check " << endl;

            CFct_impl < _Int, _Flt >::CombineLU (true, ni_diag, ia_lt_dia, ja_lt_dia,
                                                 jachar_lt_dia, a_lt_dia, ia_u_dia,
                                                 ja_u_dia, jachar_u_dia, a_u_dia, *pia_lu,
                                                 *pja_lu, *pjachar_lu, *pa_lu);
            {
               vector < _Int > ia_dummy;
               vector < _Int > ja_dummy;
               vector < char >jachar_dummy;
               vector < _Flt > a_dummy;
               ia_u_dia.swap (ia_dummy);
               ja_u_dia.swap (ja_dummy);
               jachar_u_dia.swap (jachar_dummy);
               a_u_dia.swap (a_dummy);
            }
            {
               vector < _Int > ia_dummy;
               vector < _Int > ja_dummy;
               vector < char >jachar_dummy;
               vector < _Flt > a_dummy;
               ia_lt_dia.swap (ia_dummy);
               ja_lt_dia.swap (ja_dummy);
               jachar_lt_dia.swap (jachar_dummy);
               a_lt_dia.swap (a_dummy);
            }

            int nzja_lu = (int) ((*pia_lu)[ni_diag]);

            amatr_lu.ResizeList (ni_diag);

            _Int *plist_lu = amatr_lu.GetListArr ();

            int i;

            for (i = 0; i < ni_diag; i++)
               plist_lu[i] = i;

            amatr_lu.SetNlist (ni_diag);
            amatr_lu.SetNzja (nzja_lu);
            amatr_lu.SetNzjaChar (nzja_lu);
            amatr_lu.SetNza (nzja_lu);

// Order

//      ffout << " Point 7 check " << endl;

            CMatrix < _Int, _Flt > amatr_ord;

            amatr_lu.OrderMtr (pporder_curr, amatr_ord);

            amatr_lu.Clean ();

            vector < _Int > *pia_ord = amatr_ord.GetIa ();
            vector < _Int > *pja_ord = amatr_ord.GetJa ();
            vector < char >*pjachar_ord = amatr_ord.GetJaChar ();
            vector < _Flt > *pa_ord = amatr_ord.GetA ();

// Split

//      ffout << " Point 8 check " << endl;

            CFct_impl < _Int, _Flt >::SplitLU (true, ni_diag, *pia_ord, *pja_ord,
                                               *pjachar_ord, *pa_ord, ia_l_dia, ja_l_dia,
                                               jachar_l_dia, a_l_dia, ia_u_dia, ja_u_dia,
                                               jachar_u_dia, a_u_dia);

            amatr_ord.Clean ();

// Transpose L back

            CFct_impl < _Int, _Flt >::Transpose (true, ni_diag, ia_l_dia, ja_l_dia,
                                                 jachar_l_dia, a_l_dia, ia_lt_dia,
                                                 ja_lt_dia, jachar_lt_dia, a_lt_dia);

            {
               vector < _Int > ia_dummy;
               vector < _Int > ja_dummy;
               vector < char >jachar_dummy;
               vector < _Flt > a_dummy;
               ia_l_dia.swap (ia_dummy);
               ja_l_dia.swap (ja_dummy);
               jachar_l_dia.swap (jachar_dummy);
               a_l_dia.swap (a_dummy);
            }

// Combine pair

//      ffout << " Point 9 check " << endl;

            CFct_impl < _Int, _Flt >::CombinePairs (true, ni_diag, ia_lt_dia, ja_lt_dia,
                                                    jachar_lt_dia, a_lt_dia, ia_u_dia,
                                                    ja_u_dia, jachar_u_dia, a_u_dia,
                                                    *pia_alu_dia, *pja_alu_dia,
                                                    *pjachar_alu_dia, *pa_alu_dia);

         }

// Reorder rows for all off-diagonal blocks

//      ffout << " Point 10 check " << endl;

         {
            int i;
            for (i = 1; i < nzja_curr; i++) {
               CMatrix < _Int, _Flt > aord;
               pASub_curr[i].OrderMtrRowsPairs (pporder_curr, aord);
               pASub_curr[i].ReplaceFree (aord);
            }
         }

//      ffout << " Point 11 check " << endl;

// Perform main fct computations

         {

// Create condense mask array

            int nj_blk = 0;

            int i, jblk;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               nj_blk += (int) (pblks[jblk + 1] - pblks[jblk]);
            }

            CVectorData < int >imask_col (nj_blk);
            CVectorData < int >list_col (nj_blk);
            CVectorData < int >ibsblk (nblks);

            int *pimask_col = imask_col.Ptr ();
            int *plist_col = list_col.Ptr ();
            int *pibsblk = ibsblk.Ptr ();

            for (i = 0; i < nj_blk; i++)
               pimask_col[i] = -1;

            nj_blk = 0;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               pibsblk[jblk] = nj_blk;
               nj_blk += (int) (pblks[jblk + 1] - pblks[jblk]);
            }

            CVectorData < int >ia_col (nzja_curr + 1);
            int *pia_col = ia_col.Ptr ();

            pia_col[0] = 0;

            int nlist_temp;
            _Int *plist_temp;
            _Int *pia_temp;
            _Int *pja_temp;
            char *pjachar_temp;

            int nlist_col_tot = 0;
            int nzja_tot = 0;

            int ibs, nlist_col, j, jj, k, kk;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               //nj_curr = (int) (pblks[jblk + 1] - pblks[jblk]);
               ibs = pibsblk[jblk];
               nlist_col = 0;
               nlist_temp = pASub_curr[i].GetNlist ();
               pia_temp = pASub_curr[i].GetIaArr ();
               pja_temp = pASub_curr[i].GetJaArr ();
               for (j = 0; j < nlist_temp; j++) {
                  for (k = (int) pia_temp[j]; k < pia_temp[j + 1]; k++) {
                     kk = (int) pja_temp[k];
                     if (pimask_col[ibs + kk] == -1) {
                        pimask_col[ibs + kk] = 1;
                        plist_col[ibs + nlist_col] = kk;
                        nlist_col++;
                     }
                  }
               }
               sort (plist_col + ibs, plist_col + ibs + nlist_col);
               for (j = 0; j < nlist_col; j++) {
                  jj = plist_col[ibs + j];
                  pimask_col[ibs + jj] = nlist_col_tot + j;
               }
               nlist_col_tot += nlist_col;
               pia_col[i + 1] = nlist_col_tot;
               nzja_tot += (int) pia_temp[nlist_temp];
//            ffout << " IndBlk = " << i << " jblk = " << jblk << " nj_curr = " << nj_curr << endl;
//            PrintArray (ffout," ListCol ",nlist_col,plist_col+ibs);
            }

//         PrintArray (ffout," IMaskCol ",nj_blk,pimask_col);

//      ffout << " Point 12 check " << endl;

// Compute condensed block as a sum of blocks

            CMatrix < _Int, _Flt > AU_blk;

            AU_blk.ResizeAndSetAll (nlist_col_tot, 0, nzja_tot, 0, 2 * nzja_tot);

            AU_blk.ResizeJaChar (nzja_tot);
            AU_blk.SetNzjaChar (nzja_tot);

            _Int *plist_AU_blk = AU_blk.GetListArr ();
            _Int *pia_AU_blk = AU_blk.GetIaArr ();
            _Int *pja_AU_blk = AU_blk.GetJaArr ();
            char *pjachar_AU_blk = AU_blk.GetJaCharArr ();
            _Flt *pa_AU_blk = AU_blk.GetAArr ();

            for (i = 0; i < nlist_col_tot; i++)
               plist_AU_blk[i] = (_Int) i;
            for (i = 0; i <= nlist_col_tot; i++)
               pia_AU_blk[i] = 0;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               ibs = pibsblk[jblk];
               nlist_temp = pASub_curr[i].GetNlist ();
               plist_temp = pASub_curr[i].GetListArr ();
               pia_temp = pASub_curr[i].GetIaArr ();
               for (j = 0; j < nlist_temp; j++) {
                  jj = (int) plist_temp[j];
                  pia_AU_blk[jj + 1] += (_Int) (pia_temp[j + 1] - pia_temp[j]);
               }
            }

            for (i = 0; i < nlist_col_tot; i++)
               pia_AU_blk[i + 1] = pia_AU_blk[i] + pia_AU_blk[i + 1];

            CVectorData < int >iptr (nlist_col_tot);
            int *piptr = iptr.Ptr ();

            for (i = 0; i < nlist_col_tot; i++)
               piptr[i] = (int) pia_AU_blk[i];

            _Flt *pa_temp;

            int kind, kknew;

//      ffout << " Point 13 check " << endl;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               ibs = pibsblk[jblk];
               nlist_temp = pASub_curr[i].GetNlist ();
               plist_temp = pASub_curr[i].GetListArr ();
               pia_temp = pASub_curr[i].GetIaArr ();
               pja_temp = pASub_curr[i].GetJaArr ();
               pjachar_temp = pASub_curr[i].GetJaCharArr ();
               pa_temp = pASub_curr[i].GetAArr ();
               for (j = 0; j < nlist_temp; j++) {
                  jj = (int) plist_temp[j];
                  kind = piptr[jj];
                  for (k = (int) pia_temp[j]; k < pia_temp[j + 1]; k++) {
                     kk = (int) pja_temp[k];
                     kknew = pimask_col[ibs + kk];
                     pja_AU_blk[kind] = (_Int) kknew;
                     pjachar_AU_blk[kind] = pjachar_temp[k];
                     pa_AU_blk[kind * 2] = pa_temp[k * 2];
                     pa_AU_blk[kind * 2 + 1] = pa_temp[k * 2 + 1];
                     kind++;
                  }
                  piptr[jj] = kind;
               }
            }

// Create backward transformation arrays

            CVectorData < int >listcolblk (nlist_col_tot);
            CVectorData < int >indexcol (nlist_col_tot);

            int *plistcolblk = listcolblk.Ptr ();
            int *pindexcol = indexcol.Ptr ();

            nlist_col_tot = 0;

            int njloc;

//      ffout << " Point 14 check " << endl;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               //nj_curr = (int) (pblks[jblk + 1] - pblks[jblk]);
               ibs = pibsblk[jblk];
               njloc = pia_col[i + 1] - pia_col[i];
               for (j = 0; j < njloc; j++) {
                  plistcolblk[nlist_col_tot] = jblk;
                  pindexcol[nlist_col_tot] = plist_col[ibs + j];
                  nlist_col_tot++;
               }
            }

// Perform factorization of the condensed block row

            int nmodif = 0;
            double eigmin = 1.0e100;
            double eigmax = -1.0e100;

            vector < _Int > *pia_AU_vect = AU_blk.GetIa ();
            vector < _Int > *pja_AU_vect = AU_blk.GetJa ();
            vector < char >*pjachar_AU_vect = AU_blk.GetJaChar ();
            vector < _Flt > *pa_AU_vect = AU_blk.GetA ();

            CMatrix < _Int, _Flt > U_blk;

            vector < _Int > *pia_U_vect = U_blk.GetIa ();
            vector < _Int > *pja_U_vect = U_blk.GetJa ();
            vector < char >*pjachar_U_vect = U_blk.GetJaChar ();
            vector < _Flt > *pa_U_vect = U_blk.GetA ();

//            ffout << " Before fct: AU pair: " << endl;
//            ffout << "   nlist_col_tot = " << nlist_col_tot << " ni_diag = " << ni_diag << endl;
//            ffout << "   tau2_sch = " << tau2_sch << endl;

//            CFct<_Int,_Flt>::PrintMatrix (ffout, AU_blk);

            if (strtype == 0) {

//               ofstream *pfout_debug = NULL;

//               if (index2_out == 0) {
//                  sprintf (strbuff,"ChkIluDeg_1x1_%i_%i.dat",myid_cpu,iblk);
//               } else {
//                  sprintf (strbuff,"ChkIluDeg_1x1_%i_Sch_%i.dat",myid_cpu,iblk);
//               }

//               ofstream fout_debug (strbuff);

//               pfout_debug = &fout_debug;

//               CFct_impl < _Int, _Flt >::Ilu2BlockIlu2Degree (pfout_debug, fcttype, fcttype_sch,
               CFct_impl < _Int, _Flt >::Ilu2BlockIlu2Degree (NULL, fcttype, fcttype,
                                                              pivmin, tau1, tau2,
                                                              tau2_sch, theta,
                                                              nlist_col_tot, ni_diag,
                                                              *pia_AU_vect, *pja_AU_vect,
                                                              *pjachar_AU_vect,
                                                              *pa_AU_vect, *pia_U_vect,
                                                              *pja_U_vect,
                                                              *pjachar_U_vect, *pa_U_vect,
                                                              nmodif, eigmin, eigmax);

            } else {

               CFct_impl < _Int, _Flt >::Ilu2BlockIlu2 (fcttype, fcttype, pivmin, tau1,
                                                        tau2, tau2_sch, theta,
                                                        nlist_col_tot, ni_diag,
                                                        *pia_AU_vect, *pja_AU_vect,
                                                        *pjachar_AU_vect, *pa_AU_vect,
                                                        *pia_U_vect, *pja_U_vect,
                                                        *pjachar_U_vect, *pa_U_vect,
                                                        nmodif, eigmin, eigmax);

            }

            pnmodif_thr[my_thr] += nmodif;
            if (eigmin < peigmin_thr[my_thr])
               peigmin_thr[my_thr] = eigmin;
            if (eigmax > peigmax_thr[my_thr])
               peigmax_thr[my_thr] = eigmax;

            U_blk.ResizeList (nlist_col_tot);
            _Int *plist_U = U_blk.GetListArr ();

            for (i = 0; i < nlist_col_tot; i++)
               plist_U[i] = i;

            int nzja_U_temp = (int) (*pia_U_vect)[nlist_col_tot];

            U_blk.SetNlist (nlist_col_tot);
            U_blk.SetNzja (nzja_U_temp);
            U_blk.SetNzjaChar (nzja_U_temp);
            U_blk.SetNza (nzja_U_temp * 2);

//            ffout << " After Fct: U pair: " << endl;
//            CFct<_Int,_Flt>::PrintMatrix (ffout, U_blk);

// Split U data

            int nlist_U = nlist_col_tot;
            _Int *pia_U = U_blk.GetIaArr ();
            _Int *pja_U = U_blk.GetJaArr ();
            char *pjachar_U = U_blk.GetJaCharArr ();
            _Flt *pa_U = U_blk.GetAArr ();

            if (nlist_U > 0) {

               icycleblk++;

               int nzjablk_U = 0;

               int kblk;

               for (j = 0; j < ni_diag; j++) {
                  for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                     kk = (int) pja_U[k];
                     kblk = plistcolblk[kk];
                     if (pimaskblk[kblk] != icycleblk) {
                        pimaskblk[kblk] = icycleblk;
                        plistblk[nzjablk_U] = kblk;
                        nzjablk_U++;
                     }
                  }
               }

               sort (plistblk, plistblk + nzjablk_U);

               pblockrowsLU_arr[iblk].ResizeASub (nzjablk_U);
               pblockrowsLU_arr[iblk].SetNzblk (nzjablk_U);

               CMatrix < _Int, _Flt > *pASub_U = pblockrowsLU_arr[iblk].GetASubArr ();

               CMatrix < int, float >*pHMatrU = pblockrowsLU_arr[iblk].GetHMatrStr ();

               pHMatrU->ResizeAndSetAllSp (1, 1, nzjablk_U, nzjablk_U);

               int *plist_UHblk = pHMatrU->GetListArr ();
               int *plist2_UHblk = pHMatrU->GetList2Arr ();
               int *pia_UHblk = pHMatrU->GetIaArr ();
               int *pja_UHblk = pHMatrU->GetJaArr ();
               int *pja2_UHblk = pHMatrU->GetJa2Arr ();

               plist_UHblk[0] = iblk;
               plist2_UHblk[0] = 0;
               pia_UHblk[0] = 0;
               pia_UHblk[1] = nzjablk_U;
               for (j = 0; j < nzjablk_U; j++)
                  pja_UHblk[j] = plistblk[j];
               for (j = 0; j < nzjablk_U; j++)
                  pja2_UHblk[j] = 0;

               for (j = 0; j < nzjablk_U; j++) {
                  jblk = plistblk[j];
                  pindblk[jblk] = j;
               }

               CVectorData < int >nzja_blkarr (nzjablk_U);
               int *pnzja_blkarr = nzja_blkarr.Ptr ();

               for (j = 0; j < nzjablk_U; j++)
                  pnzja_blkarr[j] = 0;

               int ind;

               for (j = 0; j < ni_diag; j++) {
                  for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                     kk = (int) pja_U[k];
                     kblk = plistcolblk[kk];
                     ind = pindblk[kblk];
                     pnzja_blkarr[ind]++;
                  }
               }

               CVectorData < CVectorData < _Int > >irows (nzjablk_U);
               CVectorData < CVectorData < _Int > >icols (nzjablk_U);
               CVectorData < CVectorData < _Flt > >elems (nzjablk_U);

               CVectorData < _Int > *pirows = irows.Ptr ();
               CVectorData < _Int > *picols = icols.Ptr ();
               CVectorData < _Flt > *pelems = elems.Ptr ();

               for (i = 0; i < nzjablk_U; i++) {
                  pirows[i].resize (pnzja_blkarr[i]);
                  picols[i].resize (pnzja_blkarr[i]);
                  pelems[i].resize (pnzja_blkarr[i] * 2);
               }

               for (j = 0; j < nzjablk_U; j++)
                  pnzja_blkarr[j] = 0;

               int kkk;

               for (j = 0; j < ni_diag; j++) {
                  for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                     kk = (int) pja_U[k];
                     kblk = plistcolblk[kk];
                     ind = pindblk[kblk];
                     _Int *ppirows = pirows[ind].Ptr ();
                     _Int *ppicols = picols[ind].Ptr ();
                     _Flt *ppelems = pelems[ind].Ptr ();
                     kkk = pnzja_blkarr[ind];
                     ppirows[kkk] = (_Int) j;
                     ppicols[kkk] = (_Int) pindexcol[kk];
                     ppelems[kkk * 2] = pa_U[k * 2];
                     ppelems[kkk * 2 + 1] = pa_U[k * 2 + 1];
                     pnzja_blkarr[ind]++;
                  }
               }

               int njrow, ibeg;

               for (j = 0; j < nzjablk_U; j++) {
                  icycleblk++;
                  nlist_temp = 0;
                  _Int *ppirows = pirows[j].Ptr ();
                  _Int *ppicols = picols[j].Ptr ();
                  _Flt *ppelems = pelems[j].Ptr ();
                  for (k = 0; k < pnzja_blkarr[j]; k++) {
                     kk = (int) ppirows[k];
                     if (pimask[kk] != icycleblk) {
                        pimask[kk] = icycleblk;
                        plist[nlist_temp] = kk;
                        nlist_temp++;
                     }
                  }
                  sort (plist, plist + nlist_temp);
                  pASub_U[j].ResizeAndSetAll (nlist_temp, 0, pnzja_blkarr[j], 0,
                                              2 * pnzja_blkarr[j]);
                  _Int *plist_temp = pASub_U[j].GetListArr ();
                  _Int *pia_temp = pASub_U[j].GetIaArr ();
                  _Int *pja_temp = pASub_U[j].GetJaArr ();
                  _Flt *pa_temp = pASub_U[j].GetAArr ();
                  for (k = 0; k < nlist_temp; k++)
                     plist_temp[k] = (_Int) plist[k];
                  for (k = 0; k < nlist_temp; k++) {
                     kk = (int) plist_temp[k];
                     plist[kk] = k;
                  }
                  for (k = 0; k <= nlist_temp; k++)
                     pia_temp[k] = 0;
                  for (k = 0; k < pnzja_blkarr[j]; k++) {
                     kk = (int) ppirows[k];
                     ind = plist[kk];
                     pia_temp[ind + 1]++;
                  }
                  for (k = 0; k < nlist_temp; k++)
                     pia_temp[k + 1] = pia_temp[k] + pia_temp[k + 1];
                  CVectorData < int >iptr (nlist_temp);
                  int *piptr = iptr.Ptr ();
                  for (k = 0; k < nlist_temp; k++)
                     piptr[k] = (int) pia_temp[k];
                  for (k = 0; k < pnzja_blkarr[j]; k++) {
                     kk = (int) ppirows[k];
                     ind = plist[kk];
                     kkk = piptr[ind];
                     pja_temp[kkk] = ppicols[k];
                     pa_temp[kkk * 2] = ppelems[k * 2];
                     pa_temp[kkk * 2 + 1] = ppelems[k * 2 + 1];
                     piptr[ind]++;
                  }
                  int njrowmax = 0;
                  for (k = 0; k < nlist_temp; k++) {
                     njrow = (int) (pia_temp[k + 1] - pia_temp[k]);
                     if (njrow > njrowmax)
                        njrowmax = njrow;
                  }
                  CVectorData < CSortInt > iiarr (njrowmax);
                  CVectorData < _Flt > elemsort (njrowmax * 2);
                  CSortInt *piiarr = iiarr.Ptr ();
                  _Flt *pelemsort = elemsort.Ptr ();
                  for (k = 0; k < nlist_temp; k++) {
                     ibeg = (int) pia_temp[k];
                     njrow = (int) (pia_temp[k + 1] - pia_temp[k]);
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        piiarr[kk - ibeg].ival = (int) pja_temp[kk];
                        piiarr[kk - ibeg].i2val = kk;
                     }
                     sort (piiarr, piiarr + njrow);
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        kkk = piiarr[kk - ibeg].i2val;
                        pelemsort[(kk - ibeg) * 2] = pa_temp[kkk * 2];
                        pelemsort[(kk - ibeg) * 2 + 1] = pa_temp[kkk * 2 + 1];
                     }
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        kkk = piiarr[kk - ibeg].ival;
                        pja_temp[kk] = (_Int) kkk;
                        pa_temp[kk * 2] = pelemsort[(kk - ibeg) * 2];
                        pa_temp[kk * 2 + 1] = pelemsort[(kk - ibeg) * 2 + 1];
                     }
                  }
               }

            } else {

               pblockrowsLU_arr[iblk].ResizeASub (1);
               pblockrowsLU_arr[iblk].SetNzblk (1);

               CMatrix < int, float >*pHMatrU = pblockrowsLU_arr[iblk].GetHMatrStr ();

               pHMatrU->ResizeAndSetAllSp (1, 1, 1, 1);

               int *plist_UHblk = pHMatrU->GetListArr ();
               int *plist2_UHblk = pHMatrU->GetList2Arr ();
               int *pia_UHblk = pHMatrU->GetIaArr ();
               int *pja_UHblk = pHMatrU->GetJaArr ();
               int *pja2_UHblk = pHMatrU->GetJa2Arr ();

               plist_UHblk[0] = iblk;
               plist2_UHblk[0] = 0;
               pia_UHblk[0] = 0;
               pia_UHblk[1] = 1;
               pja_UHblk[0] = iblk;
               pja2_UHblk[0] = 0;

            }

//            ffout << " U hblk fct: " << endl;
//            pblockrowsLU_arr[iblk].PrintHMatrix (ffout);

// Split Schur complement data

            {

               int nlistblk_schur = nzja_curr - 1;
               int nzjablk_schur = 0;

               int kblk, nlistblk_temp;

               for (i = 1; i < nzja_curr; i++) {
                  jblk = pja_curr[i];
                  icycleblk++;
                  njloc = pia_col[i + 1] - pia_col[i];
                  nlistblk_temp = 0;
                  for (j = pia_col[i]; j < pia_col[i + 1]; j++) {
                     for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                        kk = (int) pja_U[k];
                        kblk = plistcolblk[kk];
                        if (pimaskblk[kblk] != icycleblk) {
                           pimaskblk[kblk] = icycleblk;
                           nlistblk_temp++;
                        }
                     }
                  }
                  nzjablk_schur += nlistblk_temp;
               }

               pschur_arr[iblk].ResizeASub (nzjablk_schur);
               pschur_arr[iblk].SetNzblk (nzjablk_schur);

               CMatrix < _Int, _Flt > *pASub_Schur = pschur_arr[iblk].GetASubArr ();

               CMatrix < int, float >*pHMatrSchur = pschur_arr[iblk].GetHMatrStr ();

               pHMatrSchur->ResizeAndSetAllSp (nlistblk_schur, nlistblk_schur,
                                               nzjablk_schur, nzjablk_schur);

               int *plist_SchurHblk = pHMatrSchur->GetListArr ();
               int *plist2_SchurHblk = pHMatrSchur->GetList2Arr ();
               int *pia_SchurHblk = pHMatrSchur->GetIaArr ();
               int *pja_SchurHblk = pHMatrSchur->GetJaArr ();
               int *pja2_SchurHblk = pHMatrSchur->GetJa2Arr ();

               int iloc;

               for (i = 1; i < nzja_curr; i++) {
                  iloc = i - 1;
                  jblk = pja_curr[i];
                  plist_SchurHblk[iloc] = jblk;
                  plist2_SchurHblk[iloc] = 0;
               }

               pia_SchurHblk[0] = 0;
               nzjablk_schur = 0;

               for (i = 1; i < nzja_curr; i++) {
                  icycleblk++;
                  iloc = i;
                  jblk = pja_curr[i];
                  njloc = pia_col[iloc + 1] - pia_col[iloc];
                  nlistblk_temp = 0;
                  for (j = pia_col[iloc]; j < pia_col[iloc + 1]; j++) {
                     for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                        kk = (int) pja_U[k];
                        kblk = plistcolblk[kk];
                        if (pimaskblk[kblk] != icycleblk) {
                           pimaskblk[kblk] = icycleblk;
                           plistblk[nlistblk_temp] = kblk;
                           nlistblk_temp++;
                        }
                     }
                  }
                  sort (plistblk, plistblk + nlistblk_temp);
                  for (j = 0; j < nlistblk_temp; j++)
                     pja_SchurHblk[nzjablk_schur + j] = plistblk[j];
                  for (j = 0; j < nlistblk_temp; j++)
                     pja2_SchurHblk[nzjablk_schur + j] = 0;
                  nzjablk_schur += nlistblk_temp;
                  pia_SchurHblk[iloc] = nzjablk_schur;
               }

//            ffout << " Schur str: " << endl;
//            PrintArray (ffout," plist_SchurHblk = ",nlistblk_schur,plist_SchurHblk);
//            PrintArray (ffout," plist2_SchurHblk = ",nlistblk_schur,plist2_SchurHblk);
//            PrintArray (ffout," pia_SchurHblk = ",nlistblk_schur+1,pia_SchurHblk);
//            PrintArray (ffout," pja_SchurHblk = ",nzjablk_schur,pja_SchurHblk);
//            PrintArray (ffout," pja2_SchurHblk = ",nzjablk_schur,pja2_SchurHblk);

               CVectorData < int >nzja_blkarr (nzjablk_schur);
               int *pnzja_blkarr = nzja_blkarr.Ptr ();

               for (j = 0; j < nzjablk_schur; j++)
                  pnzja_blkarr[j] = 0;

               int ind;

               for (i = 1; i < nzja_curr; i++) {
                  iloc = i;
                  for (j = pia_SchurHblk[iloc - 1]; j < pia_SchurHblk[iloc]; j++) {
                     kblk = pja_SchurHblk[j];
                     pindblk[kblk] = j;
                  }
                  for (j = pia_col[iloc]; j < pia_col[iloc + 1]; j++) {
                     for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                        kk = (int) pja_U[k];
                        kblk = plistcolblk[kk];
                        ind = pindblk[kblk];
                        pnzja_blkarr[ind]++;
                     }
                  }
               }

               CVectorData < CVectorData < _Int > >irows (nzjablk_schur);
               CVectorData < CVectorData < _Int > >icols (nzjablk_schur);
               CVectorData < CVectorData < char > >chars (nzjablk_schur);
               CVectorData < CVectorData < _Flt > >elems (nzjablk_schur);

               CVectorData < _Int > *pirows = irows.Ptr ();
               CVectorData < _Int > *picols = icols.Ptr ();
               CVectorData < char >*pchars = chars.Ptr ();
               CVectorData < _Flt > *pelems = elems.Ptr ();

               for (i = 0; i < nzjablk_schur; i++) {
                  pirows[i].resize (pnzja_blkarr[i]);
                  picols[i].resize (pnzja_blkarr[i]);
                  pchars[i].resize (pnzja_blkarr[i]);
                  pelems[i].resize (pnzja_blkarr[i] * 2);
               }

               for (j = 0; j < nzjablk_schur; j++)
                  pnzja_blkarr[j] = 0;

               int kkk;

               for (i = 1; i < nzja_curr; i++) {
                  iloc = i;
                  for (j = pia_SchurHblk[iloc - 1]; j < pia_SchurHblk[iloc]; j++) {
                     kblk = pja_SchurHblk[j];
                     pindblk[kblk] = j;
                  }
                  for (j = pia_col[iloc]; j < pia_col[iloc + 1]; j++) {
                     for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                        kk = (int) pja_U[k];
                        kblk = plistcolblk[kk];
                        ind = pindblk[kblk];
                        _Int *ppirows = pirows[ind].Ptr ();
                        _Int *ppicols = picols[ind].Ptr ();
                        char *ppchars = pchars[ind].Ptr ();
                        _Flt *ppelems = pelems[ind].Ptr ();
                        kkk = pnzja_blkarr[ind];
                        ppirows[kkk] = (_Int) pindexcol[j];
                        ppicols[kkk] = (_Int) pindexcol[kk];
                        ppchars[kkk] = pjachar_U[k];
                        ppelems[kkk * 2] = pa_U[k * 2];
                        ppelems[kkk * 2 + 1] = pa_U[k * 2 + 1];
                        pnzja_blkarr[ind]++;
                     }
                  }
               }

               int njrow, ibeg;

               for (j = 0; j < nzjablk_schur; j++) {
                  icycleblk++;
                  nlist_temp = 0;
                  _Int *ppirows = pirows[j].Ptr ();
                  _Int *ppicols = picols[j].Ptr ();
                  char *ppchars = pchars[j].Ptr ();
                  _Flt *ppelems = pelems[j].Ptr ();
                  for (k = 0; k < pnzja_blkarr[j]; k++) {
                     kk = (int) ppirows[k];
                     if (pimask[kk] != icycleblk) {
                        pimask[kk] = icycleblk;
                        plist[nlist_temp] = kk;
                        nlist_temp++;
                     }
                  }
                  sort (plist, plist + nlist_temp);
                  pASub_Schur[j].ResizeAndSetAll (nlist_temp, 0, pnzja_blkarr[j], 0,
                                                  pnzja_blkarr[j] * 2);
                  pASub_Schur[j].ResizeJaChar (pnzja_blkarr[j]);
                  pASub_Schur[j].SetNzjaChar (pnzja_blkarr[j]);
                  _Int *plist_temp = pASub_Schur[j].GetListArr ();
                  _Int *pia_temp = pASub_Schur[j].GetIaArr ();
                  _Int *pja_temp = pASub_Schur[j].GetJaArr ();
                  char *pjachar_temp = pASub_Schur[j].GetJaCharArr ();
                  _Flt *pa_temp = pASub_Schur[j].GetAArr ();
                  for (k = 0; k < nlist_temp; k++)
                     plist_temp[k] = (_Int) plist[k];
                  for (k = 0; k < nlist_temp; k++) {
                     kk = (int) plist_temp[k];
                     plist[kk] = k;
                  }
                  for (k = 0; k <= nlist_temp; k++)
                     pia_temp[k] = 0;
                  for (k = 0; k < pnzja_blkarr[j]; k++) {
                     kk = (int) ppirows[k];
                     ind = plist[kk];
                     pia_temp[ind + 1]++;
                  }
                  for (k = 0; k < nlist_temp; k++)
                     pia_temp[k + 1] = pia_temp[k] + pia_temp[k + 1];
                  CVectorData < int >iptr (nlist_temp);
                  int *piptr = iptr.Ptr ();
                  for (k = 0; k < nlist_temp; k++)
                     piptr[k] = (int) pia_temp[k];
                  for (k = 0; k < pnzja_blkarr[j]; k++) {
                     kk = (int) ppirows[k];
                     ind = plist[kk];
                     kkk = piptr[ind];
                     pja_temp[kkk] = ppicols[k];
                     pjachar_temp[kkk] = ppchars[k];
                     pa_temp[kkk * 2] = ppelems[k * 2];
                     pa_temp[kkk * 2 + 1] = ppelems[k * 2 + 1];
                     piptr[ind]++;
                  }
                  int njrowmax = 0;
                  for (k = 0; k < nlist_temp; k++) {
                     njrow = (int) (pia_temp[k + 1] - pia_temp[k]);
                     if (njrow > njrowmax)
                        njrowmax = njrow;
                  }
                  CVectorData < CSortInt > iiarr (njrowmax);
                  CVectorData < char >charsort (njrowmax);
                  CVectorData < _Flt > elemsort (njrowmax * 2);
                  CSortInt *piiarr = iiarr.Ptr ();
                  char *pcharsort = charsort.Ptr ();
                  _Flt *pelemsort = elemsort.Ptr ();
                  for (k = 0; k < nlist_temp; k++) {
                     ibeg = (int) pia_temp[k];
                     njrow = (int) (pia_temp[k + 1] - pia_temp[k]);
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        piiarr[kk - ibeg].ival = (int) pja_temp[kk];
                        piiarr[kk - ibeg].i2val = kk;
                     }
                     sort (piiarr, piiarr + njrow);
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        kkk = piiarr[kk - ibeg].i2val;
                        pcharsort[kk - ibeg] = pjachar_temp[kkk];
                        pelemsort[(kk - ibeg) * 2] = pa_temp[kkk * 2];
                        pelemsort[(kk - ibeg) * 2 + 1] = pa_temp[kkk * 2 + 1];
                     }
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        kkk = piiarr[kk - ibeg].ival;
                        pja_temp[kk] = (_Int) kkk;
                        pjachar_temp[kk] = pcharsort[kk - ibeg];
                        pa_temp[kk * 2] = pelemsort[(kk - ibeg) * 2];
                        pa_temp[kk * 2 + 1] = pelemsort[(kk - ibeg) * 2 + 1];
                     }
                  }

               }

            }

//            ffout << " =========== Schur hblk fct: " << endl;
//            pschur_arr[iblk].PrintHMatrix (ffout);

         }

// Finally modify Schur data

         pschur_arr[iblk] %= hblk_sum;

//         ffout << " =========== Modified Schur hblk fct: " << endl;
//         pschur_arr[iblk].PrintHMatrix (ffout);

         picycle_thr[my_thr] = icycleblk;

      }
   };

//
// Perform ILU2 point factorization of the hmatrix with dynamic ordering and future diagonal modifications
//========================================================================================
   template < typename _Int, typename _Flt > void CFctThreads < _Int,
      _Flt >::Ilu2BlockIlu2 (bool _b_blk_wells, ofstream * _pfout, int _index_out,
                             int _ordlevel, int _ordtype, int _strtype, int _fcttype,
                             int _fcttype_sch, double _pivmin, double _tau1, double _tau2,
                             double _tau2_sch, double _theta, int _nblks, int _nblks1,
                             int _nblks2, long long *_blks, long long *_nzord_blks,
                             CTree & _tree0, CTree & _tree1, CBMatrix < _Int,
                             _Flt > &_alu_pair, CTree & _tree2_new, int &_nblks_new,
                             vector < long long >&_blks_new,
                             CVectorData < int >&_ordernew, CBMatrix < _Int,
                             _Flt > &_lu_pair, int &_nmodif, double &_eigmin_att,
                             double &_eigmax_att)
   {

//      char strbuff[256];
//      sprintf (strbuff,"ChkFct_1x1_%i.dat",_index_out);
//      ofstream ffout (strbuff);
//      _alu_pair.PrintHMatrix (ffout);

// Open hmatrix

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 1 " << endl;

      CMatrix < int, float >*phmatr_ALU = _alu_pair.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_ALU = _alu_pair.GetASubArr ();

      int *pia_hmatr_ALU = phmatr_ALU->GetIaArr ();
      int *pja_hmatr_ALU = phmatr_ALU->GetJaArr ();

// Allocate all necessary work arrays

      vector < vector < int > >order_arr (_nblks);
      vector < CBMatrix < _Int, _Flt > >blockrowsLU_arr (_nblks + 1);
      vector < CBMatrix < _Int, _Flt > >schur_arr (_nblks + 1);

      vector < int >*porder_arr = &order_arr[0];
      CBMatrix < _Int, _Flt > *pblockrowsLU_arr = &blockrowsLU_arr[0];
      CBMatrix < _Int, _Flt > *pschur_arr = &schur_arr[0];

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      vector < int >nmodif_thr (n_thr + 1);
      vector < double >eigmin_thr (n_thr + 1);
      vector < double >eigmax_thr (n_thr + 1);

      int *pnmodif_thr = &nmodif_thr[0];
      double *peigmin_thr = &eigmin_thr[0];
      double *peigmax_thr = &eigmax_thr[0];

      int i;

      for (i = 0; i < n_thr; i++)
         pnmodif_thr[i] = 0;
      for (i = 0; i < n_thr; i++)
         peigmin_thr[i] = 1.0e100;
      for (i = 0; i < n_thr; i++)
         peigmax_thr[i] = -1.0e100;

      int nimax = 0;
      int niloc;

      for (i = 0; i < _nblks; i++) {
         niloc = (int) (_blks[i + 1] - _blks[i]);
         if (niloc > nimax)
            nimax = niloc;
      }

      vector < int >icycle_thr (n_thr + 1);
      vector < CVectorData < int > >imask_thr (n_thr + 1);

      int *picycle_thr = &icycle_thr[0];
      CVectorData < int >*pimask_thr = &imask_thr[0];

      for (i = 0; i < n_thr; i++)
         picycle_thr[i] = -1;

// Perform factorization of first two subtrees level by level

      vector < CBlockFctTree < _Int, _Flt > >fct_12 (_nblks + 1);
      CBlockFctTree < _Int, _Flt > *pfct_12 = &fct_12[0];

      for (i = 0; i < _nblks1 + _nblks2; i++) {

         CBlockFctTree < _Int, _Flt > *task = pfct_12 + i;

         task->myid_cpu = _index_out;
         task->index2_out = 0;
         task->ordlevel = _ordlevel;
         task->ordtype = _ordtype;
         task->strtype = _strtype;
         task->fcttype = _fcttype;
         task->pivmin = _pivmin;
         task->tau1 = _tau1;
         task->tau2 = _tau2;
         task->tau2_sch = _tau2_sch;
         task->theta = _theta;
         task->nblks = _nblks;
         task->nblks1 = _nblks1;
         task->nblks2 = _nblks2;
         task->pblks = _blks;
         task->iblk = i;
         task->ptree_1 = &_tree0;
         task->ptree_2 = &_tree1;
         task->pia_ALU = pia_hmatr_ALU;
         task->pja_ALU = pja_hmatr_ALU;
         task->pASub_ALU = pA_sub_ALU;
         task->porder_arr = porder_arr;
         task->pblockrowsLU_arr = pblockrowsLU_arr;
         task->pschur_arr = pschur_arr;
         task->nimax = nimax;
         task->picycle_thr = picycle_thr;
         task->pimask_thr = pimask_thr;
         task->pnmodif_thr = pnmodif_thr;
         task->peigmin_thr = peigmin_thr;
         task->peigmax_thr = peigmax_thr;

      }

// Perform computations according to two tree's

      int nlev_1 = _tree0.GetNlev ();
      int *pnnodes_lev_1 = _tree0.GetNNodesLev ();
      vector < int >*pnodeslevlist_1 = _tree0.GetNodesLevList ();

      int nlev_2 = _tree1.GetNlev ();
      int *pnnodes_lev_2 = _tree1.GetNNodesLev ();
      vector < int >*pnodeslevlist_2 = _tree1.GetNodesLevList ();

      int nlev_12_min = nlev_1;
      if (nlev_2 < nlev_12_min)
         nlev_12_min = nlev_2;

      int ilev, nnodes_curr_1, nnodes_curr_2, iblk;
      int *ppnodeslevlist_1;
      int *ppnodeslevlist_2;

// Modify params

      ilev = nlev_1 - 1;

      nnodes_curr_1 = 0;
      ppnodeslevlist_1 = NULL;

      if (nlev_1 > 0) {
         nnodes_curr_1 = pnnodes_lev_1[ilev];
         ppnodeslevlist_1 = &pnodeslevlist_1[ilev][0];
      }

      for (i = 0; i < nnodes_curr_1; i++) {
         iblk = ppnodeslevlist_1[i];
         pfct_12[iblk].ordlevel = -1;
      }

      ilev = nlev_2 - 1;

      nnodes_curr_2 = 0;
      ppnodeslevlist_2 = NULL;

      if (nlev_2 > 0) {
         nnodes_curr_2 = pnnodes_lev_2[ilev];
         ppnodeslevlist_2 = &pnodeslevlist_2[ilev][0];
      }

      for (i = 0; i < nnodes_curr_2; i++) {
         iblk = _nblks1 + ppnodeslevlist_2[i];
         pfct_12[iblk].ordlevel = -1;
      }

// Run

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 2 " << endl;

      for (ilev = nlev_1 - 1; ilev >= nlev_12_min; ilev--) {

         nnodes_curr_1 = pnnodes_lev_1[ilev];
         ppnodeslevlist_1 = &pnodeslevlist_1[ilev][0];

//      if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 2.1 Ilev = " << ilev << " nnodes_curr_1 = " << nnodes_curr_1 << endl;

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_1; ipar++) {
            int iblk_temp = ppnodeslevlist_1[ipar];
            pfct_12[iblk_temp].Execute ();
         }

//      if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 2.2 After " << endl;

      }

      for (ilev = nlev_2 - 1; ilev >= nlev_12_min; ilev--) {

         nnodes_curr_2 = pnnodes_lev_2[ilev];
         ppnodeslevlist_2 = &pnodeslevlist_2[ilev][0];

//      if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 2.3 Ilev = " << ilev << " nnodes_curr_1 = " << nnodes_curr_1 << endl;

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_2; ipar++) {
            int iblk_temp = _nblks1 + ppnodeslevlist_2[ipar];
            pfct_12[iblk_temp].Execute ();
         }

//      if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 2.4 After " << endl;

      }

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 3 " << endl;

      for (ilev = nlev_12_min - 1; ilev >= 0; ilev--) {

         nnodes_curr_1 = pnnodes_lev_1[ilev];
         ppnodeslevlist_1 = &pnodeslevlist_1[ilev][0];

         nnodes_curr_2 = pnnodes_lev_2[ilev];
         ppnodeslevlist_2 = &pnodeslevlist_2[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_1 + nnodes_curr_2; ipar++) {
            int iblk_temp;
            if (ipar < nnodes_curr_1) {
               iblk_temp = ppnodeslevlist_1[ipar];
            } else {
               iblk_temp = _nblks1 + ppnodeslevlist_2[ipar - nnodes_curr_1];
            }
            pfct_12[iblk_temp].Execute ();
         }

      }

// Create last part of initial matrix data as separate hmatrix

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 4 " << endl;

      int nblks12 = _nblks1 + _nblks2;
      int nblks3 = _nblks - nblks12;

      CBMatrix < _Int, _Flt > alu_pair_last;

      {

         int nzjablk3 = 0;

         for (i = nblks12; i < _nblks; i++) {
            nzjablk3 += (pia_hmatr_ALU[i + 1] - pia_hmatr_ALU[i]);
         }

         alu_pair_last.SetNzblk (nzjablk3);
         alu_pair_last.ResizeASub (nzjablk3);

         CMatrix < int, float >*phmatr_ALU_3 = alu_pair_last.GetHMatrStr ();
         CMatrix < _Int, _Flt > *pA_sub_ALU_3 = alu_pair_last.GetASubArr ();

         phmatr_ALU_3->ResizeAndSetAllSp (nblks3, nblks3, nzjablk3, nzjablk3);

         int *plist_hmatr_ALU_3 = phmatr_ALU_3->GetListArr ();
         int *plist2_hmatr_ALU_3 = phmatr_ALU_3->GetList2Arr ();
         int *pia_hmatr_ALU_3 = phmatr_ALU_3->GetIaArr ();
         int *pja_hmatr_ALU_3 = phmatr_ALU_3->GetJaArr ();
         int *pja2_hmatr_ALU_3 = phmatr_ALU_3->GetJa2Arr ();

         for (i = 0; i < nblks3; i++)
            plist_hmatr_ALU_3[i] = i + nblks12;
         for (i = 0; i < nblks3; i++)
            plist2_hmatr_ALU_3[i] = 0;
         for (i = 0; i < nzjablk3; i++)
            pja2_hmatr_ALU_3[i] = 0;

         nzjablk3 = 0;
         pia_hmatr_ALU_3[0] = 0;

         int j;

         for (i = nblks12; i < _nblks; i++) {
            for (j = pia_hmatr_ALU[i]; j < pia_hmatr_ALU[i + 1]; j++) {
               pja_hmatr_ALU_3[nzjablk3] = pja_hmatr_ALU[j];
               pA_sub_ALU_3[nzjablk3].ReplaceFree (pA_sub_ALU[j]);
               nzjablk3++;
            }
            pia_hmatr_ALU_3[i + 1 - nblks12] = nzjablk3;
         }

      }

// Add jachar data

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 5 " << endl;

      {

         CMatrix < int, float >*phmatr_ALU_3 = alu_pair_last.GetHMatrStr ();
         CMatrix < _Int, _Flt > *pA_sub_ALU_3 = alu_pair_last.GetASubArr ();

         int nzja_hmatr_ALU_3 = phmatr_ALU_3->GetNzja ();

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nzja_hmatr_ALU_3; ipar++) {
            int nzja_temp = pA_sub_ALU_3[ipar].GetNzja ();
            pA_sub_ALU_3[ipar].ResizeJaChar (nzja_temp);
            pA_sub_ALU_3[ipar].SetNzjaChar (nzja_temp);
            char *pjachar_temp = pA_sub_ALU_3[ipar].GetJaCharArr ();
            int i;
            for (i = 0; i < nzja_temp; i++)
               pjachar_temp[i] = 0;
         }

      }

// Add Schur complements in parallel
/*
      {
//      *_pfout << " nblks = " << _nblks << endl;
//     *_pfout << " nblks1 = " << _nblks1 << endl;
//      *_pfout << " nblks2 = " << _nblks2 << endl;
         if (_nblks1 != 0) {
//         *_pfout << " Schur 1 " << endl;
            CMatrix < int, float >*phmatr_shur_1 =
               pschur_arr[_nblks1 - 1].GetHMatrStr ();
            int nzblk = phmatr_shur_1->GetNzja ();
            CMatrix < _Int, _Flt > *pA_sub_shur_1 =
               pschur_arr[_nblks1 - 1].GetASubArr ();
//         if (_pfout != NULL) CFct<int,float>::PrintMatrix (*_pfout, *phmatr_shur_1);
            if (_pfout != NULL) {
               int i;
               for (i = 0; i < nzblk; i++) {
//               *_pfout << " Block = " << i << endl;
//               CFct<_Int,_Flt>::PrintMatrix (*_pfout, pA_sub_shur_1[i]);
               }
            }
         }
         if (_nblks2 != 0) {
//         *_pfout << " Schur 2 " << endl;
            CMatrix < int, float >*phmatr_shur_2 =
               pschur_arr[nblks12 - 1].GetHMatrStr ();
            int nzblk = phmatr_shur_2->GetNzja ();
            CMatrix < _Int, _Flt > *pA_sub_shur_2 =
               pschur_arr[nblks12 - 1].GetASubArr ();
//         if (_pfout != NULL) CFct<int,float>::PrintMatrix (*_pfout, *phmatr_shur_2);
            if (_pfout != NULL) {
               int i;
               for (i = 0; i < nzblk; i++) {
//               *_pfout << " Block = " << i << endl;
                  CFct < _Int, _Flt >::PrintMatrix (*_pfout,
                                                    pA_sub_shur_2[i]);
               }
            }
         }
//      *_pfout << " Pair " << endl;
         CMatrix < int, float >*phmatr_pair = alu_pair_last.GetHMatrStr ();
         int nzblk = phmatr_pair->GetNzja ();
         CMatrix < _Int, _Flt > *pA_sub_pair = alu_pair_last.GetASubArr ();
//      if (_pfout != NULL) CFct<int,float>::PrintMatrix (*_pfout, *phmatr_pair);
         if (_pfout != NULL) {
            int i;
            for (i = 0; i < nzblk; i++) {
//            *_pfout << " Block = " << i << endl;
//            CFct<_Int,_Flt>::PrintMatrix (*_pfout, pA_sub_pair[i]);
            }
         }
      }
*/
      if (_nblks1 != 0) {
         alu_pair_last.AddReplacePairsThr (true, pschur_arr[_nblks1 - 1]);
      }

      if (_nblks2 != 0) {
         alu_pair_last.AddReplacePairsThr (true, pschur_arr[nblks12 - 1]);
      }
/*
      if (false) {
//      *_pfout << " Pair after add " << endl;
         CMatrix < int, float >*phmatr_pair = alu_pair_last.GetHMatrStr ();
         int nzblk = phmatr_pair->GetNzja ();
         CMatrix < _Int, _Flt > *pA_sub_pair = alu_pair_last.GetASubArr ();
//      if (_pfout != NULL) CFct<int,float>::PrintMatrix (*_pfout, *phmatr_pair);
         if (_pfout != NULL) {
            int i;
            for (i = 0; i < nzblk; i++) {
//            *_pfout << " Block = " << i << endl;
//            CFct<_Int,_Flt>::PrintMatrix (*_pfout, pA_sub_pair[i]);
            }
         }
      }
*/
      CMatrix < int, float >*phmatr_ALU_3 = alu_pair_last.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_ALU_3 = alu_pair_last.GetASubArr ();

      int nzja_hmatr_ALU_3 = phmatr_ALU_3->GetNzja ();
      int *plist_hmatr_ALU_3 = phmatr_ALU_3->GetListArr ();
      int *pia_hmatr_ALU_3 = phmatr_ALU_3->GetIaArr ();
      int *pja_hmatr_ALU_3 = phmatr_ALU_3->GetJaArr ();

// Renumber blocks and create reduced blocks partitioning

      vector < long long >blks3 (nblks3 + 1);
      long long *pblks3 = &blks3[0];

      for (i = 0; i <= nblks3; i++) {
         pblks3[i] = _blks[i + nblks12] - _blks[nblks12];
      }

      {

         for (i = 0; i < nblks3; i++)
            plist_hmatr_ALU_3[i] -= nblks12;
         for (i = 0; i < nzja_hmatr_ALU_3; i++)
            pja_hmatr_ALU_3[i] -= nblks12;

      }

      int ntot_3 = (int) pblks3[nblks3];

// Compute sparsity of filtered Schur complement matrix

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 6 " << endl;

      CMatrix < _Int, _Flt > astr_schur;

      {

         vector < int >ibsblkrow_arr (nblks3 + 1);
         int *pibsblkrow_arr = &ibsblkrow_arr[0];

         for (i = 0; i <= nblks3; i++)
            pibsblkrow_arr[i] = 0;

         {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < nblks3; ipar++) {

               int nz_curr, j, jj, ki, kii, kj, kjj, jj_char;

               nz_curr = 0;
               for (j = pia_hmatr_ALU_3[ipar]; j < pia_hmatr_ALU_3[ipar + 1]; j++) {
                  jj = pja_hmatr_ALU_3[j];
                  int nlist_temp = pA_sub_ALU_3[j].GetNlist ();
                  _Int *plist_temp = pA_sub_ALU_3[j].GetListArr ();
                  _Int *pia_temp = pA_sub_ALU_3[j].GetIaArr ();
                  _Int *pja_temp = pA_sub_ALU_3[j].GetJaArr ();
                  char *pjachar_temp = pA_sub_ALU_3[j].GetJaCharArr ();
                  for (ki = 0; ki < nlist_temp; ki++) {
                     kii = (int) plist_temp[ki];
                     for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                        kjj = (int) pja_temp[kj];
                        jj_char = pjachar_temp[kj];
                        if (((ipar == jj) && (kii == kjj))
                            || ((_ordlevel >= 0 && jj_char <= _ordlevel)
                                || _ordlevel < 0)) {
                           nz_curr++;
                        }
                     }
                  }
               }
               pibsblkrow_arr[ipar + 1] = nz_curr;
            }
         }

         for (i = 0; i < nblks3; i++)
            pibsblkrow_arr[i + 1] = pibsblkrow_arr[i] + pibsblkrow_arr[i + 1];

         int nzjatot_3 = pibsblkrow_arr[nblks3];

         astr_schur.ResizeAndSetAllSp (ntot_3, 0, nzjatot_3, 0);

         _Int *plist_schur = astr_schur.GetListArr ();
         _Int *pia_schur = astr_schur.GetIaArr ();
         _Int *pja_schur = astr_schur.GetJaArr ();

         for (i = 0; i < ntot_3; i++)
            plist_schur[i] = i;

         pia_schur[0] = 0;

         int iend;

         for (i = 0; i < nblks3; i++) {
            iend = (int) pblks3[i + 1] - 1;
            pia_schur[iend + 1] = pibsblkrow_arr[i + 1];
         }

         {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < nblks3; ipar++) {

               int i, j, jj, ki, kii, kj, kjj, jj_char;
               int ibeg, ibs, niloc, k;

               ibeg = (int) pblks3[ipar];
               niloc = (int) (pblks3[ipar + 1] - pblks3[ipar]);
               ibs = pibsblkrow_arr[ipar];
               vector < int >ialoc (niloc + 1);
               int *pialoc = &ialoc[0];
               for (j = 0; j <= niloc; j++)
                  pialoc[j] = 0;
               for (j = pia_hmatr_ALU_3[ipar]; j < pia_hmatr_ALU_3[ipar + 1]; j++) {
                  jj = pja_hmatr_ALU_3[j];
                  int nlist_temp = pA_sub_ALU_3[j].GetNlist ();
                  _Int *plist_temp = pA_sub_ALU_3[j].GetListArr ();
                  _Int *pia_temp = pA_sub_ALU_3[j].GetIaArr ();
                  _Int *pja_temp = pA_sub_ALU_3[j].GetJaArr ();
                  char *pjachar_temp = pA_sub_ALU_3[j].GetJaCharArr ();
                  for (ki = 0; ki < nlist_temp; ki++) {
                     kii = (int) plist_temp[ki];
                     for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                        kjj = (int) pja_temp[kj];
                        jj_char = pjachar_temp[kj];
                        if (((ipar == jj) && (kii == kjj))
                            || ((_ordlevel >= 0 && jj_char <= _ordlevel)
                                || _ordlevel < 0)) {
                           pialoc[kii + 1]++;
                        }
                     }
                  }
               }
               for (i = 0; i < niloc; i++)
                  pialoc[i + 1] = pialoc[i] + pialoc[i + 1];
               for (i = 1; i < niloc; i++)
                  pia_schur[ibeg + i] = ibs + pialoc[i];
               for (j = pia_hmatr_ALU_3[ipar]; j < pia_hmatr_ALU_3[ipar + 1]; j++) {
                  jj = pja_hmatr_ALU_3[j];
                  int njloc = (int) (pblks3[jj + 1] - pblks3[jj]);
                  if (njloc < 0) {
//                  if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Error !!! Strange block Iblk = " << ipar << " JJBlk = " << jj << " Nblks3 = " << nblks3 << " njloc = " << njloc << endl;
                  }
                  int nlist_temp = pA_sub_ALU_3[j].GetNlist ();
                  _Int *plist_temp = pA_sub_ALU_3[j].GetListArr ();
                  _Int *pia_temp = pA_sub_ALU_3[j].GetIaArr ();
                  _Int *pja_temp = pA_sub_ALU_3[j].GetJaArr ();
                  char *pjachar_temp = pA_sub_ALU_3[j].GetJaCharArr ();
                  for (ki = 0; ki < nlist_temp; ki++) {
                     kii = (int) plist_temp[ki];
                     for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                        kjj = (int) pja_temp[kj];
                        if (kjj < 0 || kjj >= njloc) {
//                        if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Error !!! Strange index kjj = " << kjj << " Njloc = " << njloc << " Ibs = " << ibs << " K = " << pialoc[kii] << endl;
                        }
                        jj_char = pjachar_temp[kj];
                        if (((ipar == jj) && (kii == kjj))
                            || ((_ordlevel >= 0 && jj_char <= _ordlevel)
                                || _ordlevel < 0)) {
                           k = pialoc[kii];
                           pja_schur[ibs + k] = kjj + (int) pblks3[jj];
                           pialoc[kii]++;
                        }
                     }
                  }
               }
            }
         }

      }

// Compute new ordering and partitioning

      int ni_wells = 0;
      if (_b_blk_wells)
         ni_wells = (int) (_blks[_nblks] - _blks[_nblks - 1]);

      int n_thr_part = 1;

#ifdef USE_THREADS
      n_thr_part = omp_get_max_threads ();
#endif

      int nparts = n_thr_part * 2;

      int nparts_pwr2 = 1;

      while (nparts_pwr2 < nparts)
         nparts_pwr2 *= 2;

      vector < int >order3 (1);
      int nblks_tree3 = 0;
      vector < long long >blks_tree3 (1);

      CTree tree_dummy;

      _tree2_new = tree_dummy;

      int ntot_3_wells = ntot_3 - ni_wells;

      if (ntot_3 > 0) {

         vector < _Int > *pia_schur_vect = astr_schur.GetIa ();
         vector < _Int > *pja_schur_vect = astr_schur.GetJa ();

// Symmetrize sparsity

         vector < _Int > ia_symm;
         vector < _Int > ja_symm;

         CFct_impl < _Int, _Flt >::SymmetrizeSparsity (ntot_3, *pia_schur_vect,
                                                       *pja_schur_vect, ia_symm, ja_symm);

         int nzja_symm = (int) ia_symm[ntot_3];

// Compute local ND ordering

         vector < _Int > ia_symm_subm;
         vector < _Int > ja_symm_subm;

         if (!_b_blk_wells) {
            CFct_impl < _Int, _Flt >::ComputeOptimalOrder (-1, ntot_3, ia_symm, ja_symm,
                                                           order3);
         } else {
            ia_symm_subm.resize (ntot_3_wells + 1);
            ja_symm_subm.resize (nzja_symm + 1);
            ia_symm_subm[0] = 0;
            int nzja_symm_subm = 0;
            int j, jj;
            for (i = 0; i < ntot_3_wells; i++) {
               for (j = (int) ia_symm[i]; j < ia_symm[i + 1]; j++) {
                  jj = (int) ja_symm[j];
                  if (jj < ntot_3_wells) {
                     ja_symm_subm[nzja_symm_subm] = jj;
                     nzja_symm_subm++;
                  }
               }
               ia_symm_subm[i + 1] = nzja_symm_subm;
            }
            CFct_impl < _Int, _Flt >::ComputeOptimalOrder (-1, ntot_3_wells, ia_symm_subm,
                                                           ja_symm_subm, order3);
            vector < int >order_ext (ntot_3 + 1);
            for (i = 0; i < ntot_3_wells; i++)
               order_ext[i] = order3[i];
            for (i = ntot_3_wells; i < ntot_3; i++)
               order_ext[i] = i;
            order3.swap (order_ext);
         }

// Compute reordered sparsity

         vector < _Int > ia_ord;
         vector < _Int > ja_ord;

         vector < _Int > ia_ord_subm;
         vector < _Int > ja_ord_subm;

         if (!_b_blk_wells) {
            CFct_impl < _Int, _Flt >::ReorderMatrixSp (ntot_3, order3, ia_symm, ja_symm,
                                                       ia_ord, ja_ord);
         } else {
            CFct_impl < _Int, _Flt >::ReorderMatrixSp (ntot_3_wells, order3, ia_symm_subm,
                                                       ja_symm_subm, ia_ord_subm,
                                                       ja_ord_subm);
         }

// Create extended tree

         CTree tree3 (nparts_pwr2, 2);

         int nblks_tree3_temp;
         vector < int >blks_tree3_temp;

         if (!_b_blk_wells) {
            CFct_impl < _Int, _Flt >::FindSeparatorsForTree (ntot_3, ia_ord, ja_ord,
                                                             tree3, nblks_tree3_temp,
                                                             blks_tree3_temp);
         } else {
            CFct_impl < _Int, _Flt >::FindSeparatorsForTree (ntot_3_wells, ia_ord_subm,
                                                             ja_ord_subm, tree3,
                                                             nblks_tree3_temp,
                                                             blks_tree3_temp);
            CTree tree3_ext (nparts_pwr2 * 2, 2);
            int nnodes_ini = tree3.GetNnodes ();
            int nnodes_ext = tree3_ext.GetNnodes ();
            int nblks_tree3_temp_ext = nnodes_ext;
            vector < int >blks_tree3_temp_ext (nnodes_ext + 1);
            for (i = 0; i <= nnodes_ext; i++)
               blks_tree3_temp_ext[i] = 0;
            for (i = 0; i < nnodes_ini; i++)
               blks_tree3_temp_ext[i + 1] = blks_tree3_temp[i + 1] - blks_tree3_temp[i];
            blks_tree3_temp_ext[nnodes_ext] = ni_wells;
            for (i = 0; i < nnodes_ext; i++)
               blks_tree3_temp_ext[i + 1] += blks_tree3_temp_ext[i];
            nblks_tree3_temp = nblks_tree3_temp_ext;
            blks_tree3_temp.swap (blks_tree3_temp_ext);
            tree3 = tree3_ext;
         }

         int *pblks_tree3_temp = &blks_tree3_temp[0];

// Perform filtering of zero data in a tree

         int nnodes_ini = tree3.GetNnodes ();
         int *psubtree_beg = tree3.GetSubtreeBeg ();

         vector < int >imasknd (nnodes_ini + 1);
         int *pimasknd = &imasknd[0];

         for (i = 0; i < nnodes_ini; i++)
            pimasknd[i] = 1;

         nblks_tree3 = 0;

         for (i = 0; i < nnodes_ini; i++) {
            int ibeg = psubtree_beg[i];
            if (pblks_tree3_temp[i + 1] - pblks_tree3_temp[ibeg] == 0) {
               pimasknd[i] = -1;
            } else {
               nblks_tree3++;
            }
         }

         tree3.FilterTree (pimasknd, _tree2_new);

         blks_tree3.resize (nblks_tree3 + 1);
         long long *pblks_tree3 = &blks_tree3[0];

         pblks_tree3[0] = 0;
         nblks_tree3 = 0;

         for (i = 0; i < nnodes_ini; i++) {
            if (pimasknd[i] > 0) {
               pblks_tree3[nblks_tree3 + 1] = pblks_tree3_temp[i + 1];
               nblks_tree3++;
            }
         }

      }

      int *porder3 = &order3[0];
      long long *pblks_tree3 = &blks_tree3[0];

      vector < int >iorder3 (ntot_3 + 1);
      int *piorder3 = &iorder3[0];

      for (i = 0; i < ntot_3; i++)
         piorder3[porder3[i]] = i;

// Reorder and repartition submatrix

// Split pair hblock

      CBMatrix < _Int, _Flt > al_last;
      CBMatrix < _Int, _Flt > au_last;

      CFctThreads < _Int, _Flt >::SplitPair (nblks3, pblks3, alu_pair_last, al_last,
                                             au_last);

      alu_pair_last.Clean ();

// Transpose hblock

      CBMatrix < _Int, _Flt > alt_last;

      CBMatrix < _Int, _Flt >::TransposeHMatrix (nblks3, pblks3, al_last, alt_last);

      al_last.Clean ();

// Set diagonal by zero

      CFctThreads < _Int, _Flt >::SetZeroDiag (nblks3, pblks3, alt_last);

// Combine into hmatrix

      CBMatrix < _Int, _Flt > alu_last;

      alu_last.ReplaceFree (alt_last);

      alu_last += au_last;

      alt_last.Clean ();
      au_last.Clean ();

// Compute reordered hmatrix

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 9 " << endl;

      CBMatrix < _Int, _Flt > alu_last_ord;

      CBMatrix < _Int, _Flt >::ReorderHMatrix (nblks3, pblks3, alu_last, porder3,
                                               nblks_tree3, pblks_tree3, alu_last_ord);

      alu_last.Clean ();

// Split into pairs

      CBMatrix < _Int, _Flt > alu_pair_last_ord;

      CFctThreads < _Int, _Flt >::SplitLUPair (nblks_tree3, pblks_tree3, alu_last_ord,
                                               alu_pair_last_ord);

      alu_last_ord.Clean ();

// Perform filtering of the block sparsity according to the tree with diagonal modifications

      CFctThreads < _Int, _Flt >::PairsFilterTreeModif (_tree2_new, _theta, nblks_tree3,
                                                        pblks_tree3, alu_pair_last_ord);

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 10 " << endl;

// Perform final fct for the last set of blocks

      CMatrix < int, float >*phmatr_3 = alu_pair_last_ord.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_3 = alu_pair_last_ord.GetASubArr ();

      int *pia_hmatr_3 = phmatr_3->GetIaArr ();
      int *pja_hmatr_3 = phmatr_3->GetJaArr ();

      for (i = 0; i < nblks_tree3; i++) {
         niloc = (int) (pblks_tree3[i + 1] - pblks_tree3[i]);
         if (niloc > nimax)
            nimax = niloc;
      }

      for (i = 0; i < n_thr; i++)
         picycle_thr[i] = -1;

      vector < vector < int > >order3_arr (nblks_tree3 + 1);
      vector < CBMatrix < _Int, _Flt > >blockrowsLU3_arr (nblks_tree3 + 1);
      vector < CBMatrix < _Int, _Flt > >schur3_arr (nblks_tree3 + 1);

      vector < int >*porder3_arr = &order3_arr[0];
      CBMatrix < _Int, _Flt > *pblockrowsLU3_arr = &blockrowsLU3_arr[0];
      CBMatrix < _Int, _Flt > *pschur3_arr = &schur3_arr[0];

      vector < CBlockFctTree < _Int, _Flt > >fct_3 (nblks_tree3 + 1);
      CBlockFctTree < _Int, _Flt > *pfct_3 = &fct_3[0];

      for (i = 0; i < nblks_tree3; i++) {

         CBlockFctTree < _Int, _Flt > *task = pfct_3 + i;

         task->myid_cpu = _index_out;
         task->index2_out = 1;
         task->ordlevel = _ordlevel;
         task->ordtype = _ordtype;
         task->strtype = _strtype;
         task->fcttype = _fcttype;
         task->pivmin = _pivmin;
         task->tau1 = _tau1;
         task->tau2 = _tau2;
         task->tau2_sch = _tau2_sch;
         task->theta = _theta;
         task->nblks = nblks_tree3;
         task->nblks1 = nblks_tree3;
         task->nblks2 = 0;
         task->pblks = pblks_tree3;
         task->iblk = i;
         task->ptree_1 = &_tree2_new;
         task->ptree_2 = &_tree2_new;
         task->pia_ALU = pia_hmatr_3;
         task->pja_ALU = pja_hmatr_3;
         task->pASub_ALU = pA_sub_3;
         task->porder_arr = porder3_arr;
         task->pblockrowsLU_arr = pblockrowsLU3_arr;
         task->pschur_arr = pschur3_arr;
         task->nimax = nimax;
         task->picycle_thr = picycle_thr;
         task->pimask_thr = pimask_thr;
         task->pnmodif_thr = pnmodif_thr;
         task->peigmin_thr = peigmin_thr;
         task->peigmax_thr = peigmax_thr;

      }

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 11 " << endl;

      int nlev_3 = _tree2_new.GetNlev ();
      int *pnnodes_lev_3 = _tree2_new.GetNNodesLev ();
      vector < int >*pnodeslevlist_3 = _tree2_new.GetNodesLevList ();

// Run

      int nnodes_curr_3;
      int *ppnodeslevlist_3;

      for (ilev = nlev_3 - 1; ilev >= 0; ilev--) {

         nnodes_curr_3 = pnnodes_lev_3[ilev];
         ppnodeslevlist_3 = &pnodeslevlist_3[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_3; ipar++) {
            int iblk_temp = ppnodeslevlist_3[ipar];
            pfct_3[iblk_temp].Execute ();
         }

      }

// Get subblock {12,3} into separate submatrix and repartition it

      CBMatrix < _Int, _Flt > hblock_12_3_ord;

      {
         int nzja_off = 0;

         int j, jj;

         for (i = 0; i < nblks12; i++) {
            CMatrix < int, float >*phmatr_temp = pblockrowsLU_arr[i].GetHMatrStr ();
            int nzja_temp = phmatr_temp->GetNzja ();
            int *pja_temp = phmatr_temp->GetJaArr ();
            for (j = 0; j < nzja_temp; j++) {
               jj = pja_temp[j];
               if (jj >= nblks12)
                  nzja_off++;
            }
         }

         CBMatrix < _Int, _Flt > hblock_12_3;

         hblock_12_3.SetNzblk (nzja_off);
         hblock_12_3.ResizeASub (nzja_off);

         CMatrix < int, float >*phmatr_12_3 = hblock_12_3.GetHMatrStr ();
         CMatrix < _Int, _Flt > *pASub_12_3 = hblock_12_3.GetASubArr ();

         phmatr_12_3->ResizeAndSetAllSp (nblks12, nblks12, nzja_off, nzja_off);

         int *plist_12_3 = phmatr_12_3->GetListArr ();
         int *plist2_12_3 = phmatr_12_3->GetList2Arr ();
         int *pia_12_3 = phmatr_12_3->GetIaArr ();
         int *pja_12_3 = phmatr_12_3->GetJaArr ();
         int *pja2_12_3 = phmatr_12_3->GetJa2Arr ();

         for (i = 0; i < nblks12; i++)
            plist_12_3[i] = i;
         for (i = 0; i < nblks12; i++)
            plist2_12_3[i] = 0;
         for (i = 0; i < nzja_off; i++)
            pja2_12_3[i] = 0;

         nzja_off = 0;
         pia_12_3[0] = 0;

         for (i = 0; i < nblks12; i++) {
            CMatrix < int, float >*phmatr_temp = pblockrowsLU_arr[i].GetHMatrStr ();
            CMatrix < _Int, _Flt > *pASub_temp = pblockrowsLU_arr[i].GetASubArr ();
            int nzja_temp = phmatr_temp->GetNzja ();
            int *pja_temp = phmatr_temp->GetJaArr ();
            for (j = 0; j < nzja_temp; j++) {
               jj = pja_temp[j];
               if (jj >= nblks12) {
                  pja_12_3[nzja_off] = jj - nblks12;
                  pASub_12_3[nzja_off].ReplaceFree (pASub_temp[j]);
                  nzja_off++;
               }
            }
            pia_12_3[i + 1] = nzja_off;
         }

//         ffout << " Before ReorderHMatrixColumnsPairs " << endl;
//         hblock_12_3.PrintHMatrix (ffout);

         CFctThreads < _Int, _Flt >::ReorderHMatrixColumnsPairs (nblks12, _blks, nblks3,
                                                                 pblks3, hblock_12_3,
                                                                 porder3, nblks_tree3,
                                                                 pblks_tree3,
                                                                 hblock_12_3_ord);

//         ffout << " After ReorderHMatrixColumnsPairs " << endl;
//         hblock_12_3_ord.PrintHMatrix (ffout);

      }

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 12 " << endl;

// Combine data into one hblock

      _nblks_new = nblks12 + nblks_tree3;
      _blks_new.resize (_nblks_new + 1);

      long long *p_blks_new = &_blks_new[0];

      for (i = 0; i <= nblks12; i++)
         p_blks_new[i] = _blks[i];
      for (i = 0; i < nblks_tree3; i++)
         p_blks_new[nblks12 + i + 1] = p_blks_new[nblks12] + pblks_tree3[i + 1];

      int nzja_fin = 0;

      int j, jj;

      for (i = 0; i < nblks12; i++) {
         CMatrix < int, float >*phmatr_temp = pblockrowsLU_arr[i].GetHMatrStr ();
         int nzja_temp = phmatr_temp->GetNzja ();
         int *pja_temp = phmatr_temp->GetJaArr ();
         for (j = 0; j < nzja_temp; j++) {
            jj = pja_temp[j];
            if (jj < nblks12)
               nzja_fin++;
         }
      }

      nzja_fin += hblock_12_3_ord.GetNzblk ();

      for (i = 0; i < nblks_tree3; i++) {
         nzja_fin += pblockrowsLU3_arr[i].GetNzblk ();
      }

      _lu_pair.SetNzblk (nzja_fin);
      _lu_pair.ResizeASub (nzja_fin);

      CMatrix < int, float >*phmatr_fin = _lu_pair.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pASub_fin = _lu_pair.GetASubArr ();

      phmatr_fin->ResizeAndSetAllSp (_nblks_new, _nblks_new, nzja_fin, nzja_fin);

      int *plist_fin = phmatr_fin->GetListArr ();
      int *plist2_fin = phmatr_fin->GetList2Arr ();
      int *pia_fin = phmatr_fin->GetIaArr ();
      int *pja_fin = phmatr_fin->GetJaArr ();
      int *pja2_fin = phmatr_fin->GetJa2Arr ();

      for (i = 0; i < _nblks_new; i++)
         plist_fin[i] = i;
      for (i = 0; i < _nblks_new; i++)
         plist2_fin[i] = 0;
      for (i = 0; i < nzja_fin; i++)
         pja2_fin[i] = 0;

      CMatrix < int, float >*phmatr_12_3_ord = hblock_12_3_ord.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pASub_12_3_ord = hblock_12_3_ord.GetASubArr ();

      int *pia_12_3_ord = phmatr_12_3_ord->GetIaArr ();
      int *pja_12_3_ord = phmatr_12_3_ord->GetJaArr ();

      nzja_fin = 0;
      pia_fin[0] = 0;

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 13 " << endl;

      for (i = 0; i < nblks12; i++) {
         CMatrix < int, float >*phmatr_temp = pblockrowsLU_arr[i].GetHMatrStr ();
         CMatrix < _Int, _Flt > *pASub_temp = pblockrowsLU_arr[i].GetASubArr ();
         int nzja_temp = phmatr_temp->GetNzja ();
         int *pja_temp = phmatr_temp->GetJaArr ();
         for (j = 0; j < nzja_temp; j++) {
            jj = pja_temp[j];
            if (jj < nblks12) {
               pja_fin[nzja_fin] = jj;
               pASub_fin[nzja_fin].ReplaceFree (pASub_temp[j]);
               nzja_fin++;
            }
         }
         for (j = pia_12_3_ord[i]; j < pia_12_3_ord[i + 1]; j++) {
            jj = pja_12_3_ord[j];
            pja_fin[nzja_fin] = jj + nblks12;
            pASub_fin[nzja_fin].ReplaceFree (pASub_12_3_ord[j]);
            nzja_fin++;
         }
         pia_fin[i + 1] = nzja_fin;
      }

      for (i = 0; i < nblks_tree3; i++) {
         CMatrix < int, float >*phmatr_temp = pblockrowsLU3_arr[i].GetHMatrStr ();
         CMatrix < _Int, _Flt > *pASub_temp = pblockrowsLU3_arr[i].GetASubArr ();
         int nzja_temp = phmatr_temp->GetNzja ();
         int *pja_temp = phmatr_temp->GetJaArr ();
         for (j = 0; j < nzja_temp; j++) {
            jj = pja_temp[j];
            pja_fin[nzja_fin] = jj + nblks12;
            pASub_fin[nzja_fin].ReplaceFree (pASub_temp[j]);
            nzja_fin++;
         }
         pia_fin[nblks12 + i + 1] = nzja_fin;
      }

// Reorder off diagonal data

      {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < _nblks_new; ipar++) {

            int j, jj;
            for (j = pia_fin[ipar]; j < pia_fin[ipar + 1]; j++) {
               jj = pja_fin[j];
               if (jj != ipar) {
                  int *porder_temp = NULL;
                  if (jj < nblks12) {
                     porder_temp = &porder_arr[jj][0];
                  } else {
                     porder_temp = &porder3_arr[jj - nblks12][0];
                  }

                  CMatrix < _Int, _Flt > aord;

                  pASub_fin[j].OrderMtrColsPairs (porder_temp, aord);

                  pASub_fin[j].ReplaceFree (aord);
               }
            }
         }
      }

// Form final result

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 14 " << endl;

      int ntot = (int) _blks[_nblks];
      int ntot_12 = (int) _blks[nblks12];

      _ordernew.resize (ntot);
      int *p_ordernew = _ordernew.Ptr ();

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nblks12; ipar++) {
            int ibeg, niloc, j;
            niloc = (int) (_blks[ipar + 1] - _blks[ipar]);
            ibeg = (int) _blks[ipar];
            int *pporder_temp = &porder_arr[ipar][0];
            for (j = 0; j < niloc; j++)
               p_ordernew[ibeg + j] = ibeg + pporder_temp[j];
         }
      }
      int iold;
      for (i = 0; i < nblks_tree3; i++) {
         niloc = (int) (pblks_tree3[i + 1] - pblks_tree3[i]);
         int ibeg = (int) pblks_tree3[i];
         int *pporder_temp = &porder3_arr[i][0];
         for (j = 0; j < niloc; j++) {
            iold = piorder3[ibeg + j];
            p_ordernew[ntot_12 + iold] = ntot_12 + ibeg + pporder_temp[j];
         }

      }

      _nmodif = 0;
      _eigmin_att = 1.0e100;
      _eigmax_att = -1.0e100;

      for (i = 0; i < n_thr; i++)
         _nmodif += pnmodif_thr[i];
      for (i = 0; i < n_thr; i++) {
         if (peigmin_thr[i] < _eigmin_att)
            _eigmin_att = peigmin_thr[i];
         if (peigmax_thr[i] > _eigmax_att)
            _eigmax_att = peigmax_thr[i];
      }

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: End point " << endl;
//      ffout << " Fct Result : " << endl;
//      _lu_pair.PrintHMatrix (ffout);

   }

//========================================================================================
   template < typename _Int, typename _Flt > struct CBlockSymbFctTree   /// Template structure that computes symbolic fct
   {
// Data
      int nblks;                ///< Total number of blocks
      long long *pblks;         ///< Block partitioning
      int inode;                ///< Global node number
      CTree *ptree;             ///< Tree
      SParams *pparams;         ///< Parameters
      int *pia_ALU_sp;          ///< Ia sparsity of ALU sp
      int *pja_ALU_sp;          ///< Ja sparsity of ALU sp
        CMatrix < _Int, _Flt > *pASub_ALU_sp;   ///< Blocks of ALU sp
      CTree *p_subtree_arr;     ///< Subtree for active nodes
      int *p_nblks_arr;         ///< Number of subblocks for active nodes
        vector < int >*p_blks_arr;      ///< Subblocks partitioning for active nodes
        vector < int >*porder_arr;      ///< Set of ordering arrays for each block
        CBMatrix < _Int, _Flt > *pblockrowsLU_arr;      ///< Computed LU
        CBMatrix < _Int, _Flt > *pschur_arr;    ///< Schur data
      int nimax;                ///< Maximal block size
      int *picycle_thr;         ///< Array of cycle variables for threads
        CVectorData < int >*pimask_thr; ///< Mask arrays for threads
// Functions
      void Execute ()
      {

//         char strbuff[256];
//         sprintf (strbuff,"ChkSymbFct_%i.dat",inode);
//         ofstream ffout (strbuff);

//         cout << " Inode Symb = " << inode << endl;
//         ffout << " Inode Symb = " << inode << endl;

         int ordtype = pparams->ordtype;

// Open tree

         int nnodes_tree = ptree->GetNnodes ();
         int *pnode2lev = ptree->GetNodesLevId ();
         int *pnchilds = ptree->GetNchilds ();
           vector < int >*pchilds_list = ptree->GetChildsList ();
         int *pnode2ind = ptree->GetNode2Ind ();

         int ilev = pnode2lev[inode];
         int iblk = pnode2ind[inode];

// Init mask arrays

         int my_thr = 0;

#ifdef USE_THREADS
           my_thr = omp_get_thread_num ();
#endif

         if (picycle_thr[my_thr] == -1) {
            pimask_thr[my_thr].resize (3 * nnodes_tree + 3 * nimax);
            int *pimaskblk = pimask_thr[my_thr].Ptr ();
            int *pimask = pimaskblk + nblks;
            int i;
            for (i = 0; i < nnodes_tree; i++)
                 pimaskblk[i] = -1;
            for (i = 0; i < nimax; i++)
                 pimask[i] = -1;
         }

         int icycleblk = picycle_thr[my_thr];

         int *pimaskblk = pimask_thr[my_thr].Ptr ();
         int *pimask = pimaskblk + nblks;
         int *pindblk = pimask + nimax;
         int *plistblk = pindblk + nnodes_tree;
         int *plist = plistblk + nnodes_tree;

// Determine list of childs of a node if any

         int ichild1 = -1;
         int ichild2 = -1;

         int *ppchilds = &pchilds_list[inode][0];
         if (pnchilds[inode] == 2) {
            ichild1 = ppchilds[0];
            ichild2 = ppchilds[1];
         } else if (pnchilds[inode] == 1) {
            ichild1 = ppchilds[0];
            if (ichild1 == inode)
               ichild1 = -1;
         } else {
            throw " CBlockSymbFctTree<>::Execute: wrong tree !!! ";
         }

// Add Schur hblocks

         CBMatrix < _Int, _Flt > hblk_sum;

         if (ichild1 >= 0 && ichild2 >= 0) {
//            ffout << " Hblk Ichild 1: " << endl;
//            pschur_arr[ichild1].PrintHMatrix (ffout);
//            ffout << " Hblk Ichild 2: " << endl;
//            pschur_arr[ichild2].PrintHMatrix (ffout);
            hblk_sum.ReplaceFree (pschur_arr[ichild1]);
            hblk_sum.AddReplaceSpThr (false, pschur_arr[ichild2]);
            pschur_arr[ichild1].Clean ();
            pschur_arr[ichild2].Clean ();
         } else if (ichild1 >= 0) {
            hblk_sum.ReplaceFree (pschur_arr[ichild1]);
            pschur_arr[ichild1].Clean ();
         } else if (ichild2 >= 0) {
            hblk_sum.ReplaceFree (pschur_arr[ichild2]);
            pschur_arr[ichild2].Clean ();
         }
//
//         ffout << " =========== Hblk_Schur_sum: " << endl;
//         hblk_sum.PrintHMatrix (ffout);

         if (iblk < 0) {
            pschur_arr[inode].ReplaceFree (hblk_sum);
            return;
         }
// Get set of blocks that correspond to current block row and split hblock for future add

         int nzja_curr = 0;
         vector < int >ja_curr;
         vector < CMatrix < _Int, _Flt > >ASub_curr;

         {

            int nzja_curr0 = pia_ALU_sp[iblk + 1] - pia_ALU_sp[iblk];
            int *pja_curr0 = pja_ALU_sp + pia_ALU_sp[iblk];
            CMatrix < _Int, _Flt > *pA_sub_curr0 = pASub_ALU_sp + pia_ALU_sp[iblk];

//            PrintArray (ffout," ja_curr0 = ",nzja_curr0,pja_curr0);

// Add zero jachar data if not available

            int i, j;

            for (i = 0; i < nzja_curr0; i++) {
               int nzja_temp = pA_sub_curr0[i].GetNzja ();
               int nzjachar_temp = pA_sub_curr0[i].GetNzjaChar ();
               if (nzja_temp != nzjachar_temp) {
                  pA_sub_curr0[i].ResizeJaChar (nzja_temp);
                  pA_sub_curr0[i].SetNzjaChar (nzja_temp);
                  char *pjachar_temp = pA_sub_curr0[i].GetJaCharArr ();
                  for (j = 0; j < nzja_temp; j++)
                     pjachar_temp[j] = 0;
               }
            }

            CMatrix < int, float >*phmatr_sum = hblk_sum.GetHMatrStr ();
            CMatrix < _Int, _Flt > *pA_sub_sum = hblk_sum.GetASubArr ();

            int nlist_hmatr_sum = phmatr_sum->GetNlist ();
            int nzja_hmatr_sum = phmatr_sum->GetNzja ();
            int *plist_hmatr_sum = phmatr_sum->GetListArr ();
            int *pia_hmatr_sum = phmatr_sum->GetIaArr ();
            int *pja_hmatr_sum = phmatr_sum->GetJaArr ();

            int ilist_curr = -1;

            for (i = 0; i < nlist_hmatr_sum; i++) {
               if (plist_hmatr_sum[i] == iblk)
                  ilist_curr = i;
            }

            int nzja_curr1 = 0;
            int *pja_curr1 = NULL;
            CMatrix < _Int, _Flt > *pA_sub_curr1 = NULL;

            if (ilist_curr >= 0) {

               nzja_curr1 = pia_hmatr_sum[ilist_curr + 1] - pia_hmatr_sum[ilist_curr];
               pja_curr1 = pja_hmatr_sum + pia_hmatr_sum[ilist_curr];
               pA_sub_curr1 = pA_sub_sum + pia_hmatr_sum[ilist_curr];

            }
//         PrintArray (ffout," ja_curr1 = ",nzja_curr1,pja_curr1);

            ja_curr.resize (nzja_curr0 + nzja_curr1 + 1);
            ASub_curr.resize (nzja_curr0 + nzja_curr1 + 1);

            int *pja_curr = &ja_curr[0];
            CMatrix < _Int, _Flt > *pASub_curr = &ASub_curr[0];

            nzja_curr = 0;

            int ip0 = 0;
            int ip1 = 0;

            int jj0, jj1;

            while (ip0 < nzja_curr0 || ip1 < nzja_curr1) {
               if (ip0 < nzja_curr0 && ip1 < nzja_curr1) {
                  jj0 = pja_curr0[ip0];
                  jj1 = pja_curr1[ip1];
                  if (jj0 == jj1) {
                     pja_curr[nzja_curr] = jj0;
                     pASub_curr[nzja_curr].AddBlocksSp (pA_sub_curr0[ip0],
                                                        pA_sub_curr1[ip1]);
                     pA_sub_curr0[ip0].Clean ();
                     pA_sub_curr1[ip1].Clean ();
                     ip0++;
                     ip1++;
                     nzja_curr++;
                  } else if (jj0 < jj1) {
                     pja_curr[nzja_curr] = jj0;
                     pASub_curr[nzja_curr].ReplaceFree (pA_sub_curr0[ip0]);
                     ip0++;
                     nzja_curr++;
                  } else if (jj0 > jj1) {
                     pja_curr[nzja_curr] = jj1;
                     pASub_curr[nzja_curr].ReplaceFree (pA_sub_curr1[ip1]);
                     ip1++;
                     nzja_curr++;
                  }
               } else if (ip0 < nzja_curr0) {
                  pja_curr[nzja_curr] = pja_curr0[ip0];
                  pASub_curr[nzja_curr].ReplaceFree (pA_sub_curr0[ip0]);
                  ip0++;
                  nzja_curr++;
               } else if (ip1 < nzja_curr1) {
                  pja_curr[nzja_curr] = pja_curr1[ip1];
                  pASub_curr[nzja_curr].ReplaceFree (pA_sub_curr1[ip1]);
                  ip1++;
                  nzja_curr++;
               }
            }

            if (ilist_curr >= 0) {

               int nzja_flt = nzja_hmatr_sum - nzja_curr1;

               CBMatrix < _Int, _Flt > hblk_flt;

               hblk_flt.SetNzblk (nzja_flt);
               hblk_flt.ResizeASub (nzja_flt);

               CMatrix < int, float >*phmatr_flt = hblk_flt.GetHMatrStr ();
               CMatrix < _Int, _Flt > *pA_sub_flt = hblk_flt.GetASubArr ();

               int nlist_hmatr_flt = nlist_hmatr_sum - 1;

               phmatr_flt->ResizeAndSetAllSp (nlist_hmatr_flt, nlist_hmatr_flt, nzja_flt,
                                              nzja_flt);

               int *plist_flt = phmatr_flt->GetListArr ();
               int *plist2_flt = phmatr_flt->GetList2Arr ();
               int *pia_flt = phmatr_flt->GetIaArr ();
               int *pja_flt = phmatr_flt->GetJaArr ();
               int *pja2_flt = phmatr_flt->GetJa2Arr ();

               int i;

               for (i = 0; i < nlist_hmatr_flt; i++)
                  plist2_flt[i] = 0;
               for (i = 0; i < nzja_flt; i++)
                  pja2_flt[i] = 0;

               nlist_hmatr_flt = 0;
               nzja_flt = 0;

               int j;

               for (i = 0; i < nlist_hmatr_sum; i++) {
                  if (plist_hmatr_sum[i] != iblk) {
                     plist_flt[nlist_hmatr_flt] = plist_hmatr_sum[i];
                     for (j = pia_hmatr_sum[i]; j < pia_hmatr_sum[i + 1]; j++) {
                        pja_flt[nzja_flt] = pja_hmatr_sum[j];
                        pA_sub_flt[nzja_flt].ReplaceFree (pA_sub_sum[j]);
                        nzja_flt++;
                     }
                     pia_flt[nlist_hmatr_flt + 1] = nzja_flt;
                     nlist_hmatr_flt++;
                  }
               }

               hblk_sum.ReplaceFree (hblk_flt);

            }

         }

         int *pja_curr = &ja_curr[0];
         CMatrix < _Int, _Flt > *pASub_curr = &ASub_curr[0];

//         ffout << " Iblk = " << iblk << endl;
//         PrintArray (ffout," ja_curr = ",nzja_curr,pja_curr);

//         ffout << " =========== Hblk_Schur_sum flt: " << endl;
//         hblk_sum.PrintHMatrix (ffout);

// Compute optimal ordering via diagonal block

         int ni_diag = (int) (pblks[iblk + 1] - pblks[iblk]);

         int jblk0 = pja_curr[0];

         if (jblk0 != iblk) {
            cout << " CBlockSymbFctTree<>::Execute(): error in block sparsity !" <<
               " Iblk = " << iblk << " jblk0 = " << jblk0 << " Ni_diag = " << ni_diag <<
               endl;
            throw " CBlockSymbFctTree<>::Execute(): error in block sparsity ! ";
         }

         int nlev_split = pparams->nlev_split;

         if ((ichild1 < 0 && ichild2 < 0) || ilev >= nlev_split) {
            CTree subtree (1, 2);
            p_subtree_arr[iblk] = subtree;
            p_nblks_arr[iblk] = 1;
            p_blks_arr[iblk].resize (2);
            p_blks_arr[iblk][0] = 0;
            p_blks_arr[iblk][1] = ni_diag;
            CFct < _Int, _Flt >::ComputeOptimalOrder (pASub_curr[0], ordtype,
                                                      porder_arr[iblk]);
         } else {
            pASub_curr[0].OrderNDSeparatorsForTree (nlev_split - ilev + 1, ordtype,
                                                    porder_arr[iblk], p_subtree_arr[iblk],
                                                    p_nblks_arr[iblk], p_blks_arr[iblk]);
         }

         int *pporder_curr = &porder_arr[iblk][0];

//         PrintArray (ffout," Order = ",ni_diag,pporder_curr);

// Reorder diagonal block

         {

// Split

            vector < _Int > *pia_alu_dia = pASub_curr[0].GetIa ();
            vector < _Int > *pja_alu_dia = pASub_curr[0].GetJa ();
            vector < char >*pjachar_alu_dia = pASub_curr[0].GetJaChar ();

            vector < _Int > ia_l_dia;
            vector < _Int > ja_l_dia;
            vector < char >jachar_l_dia;

            vector < _Int > ia_u_dia;
            vector < _Int > ja_u_dia;
            vector < char >jachar_u_dia;

//      ffout << " Point 5 check " << endl;

            CFct_impl < _Int, _Flt >::SplitLUSp (ni_diag, *pia_alu_dia, *pja_alu_dia,
                                                 *pjachar_alu_dia, ia_l_dia, ja_l_dia,
                                                 jachar_l_dia, ia_u_dia, ja_u_dia,
                                                 jachar_u_dia);

// Transpose L

            vector < _Int > ia_ut_dia;
            vector < _Int > ja_ut_dia;
            vector < char >jachar_ut_dia;

//      ffout << " Point 6 check " << endl;

            CFct_impl < _Int, _Flt >::TransposeSp (ni_diag, ia_u_dia, ja_u_dia,
                                                   jachar_u_dia, ia_ut_dia, ja_ut_dia,
                                                   jachar_ut_dia);

            {
               vector < _Int > ia_dummy;
               vector < _Int > ja_dummy;
               vector < char >jachar_dummy;
               ia_l_dia.swap (ia_dummy);
               ja_l_dia.swap (ja_dummy);
               jachar_l_dia.swap (jachar_dummy);
            }

// Combine sp L+U

            CMatrix < _Int, _Flt > amatr_lu;

            vector < _Int > *pia_lu = amatr_lu.GetIa ();
            vector < _Int > *pja_lu = amatr_lu.GetJa ();
            vector < char >*pjachar_lu = amatr_lu.GetJaChar ();

//      ffout << " Point 6 check " << endl;

            CFct_impl < _Int, _Flt >::CombineLUSp (ni_diag, ia_ut_dia, ja_ut_dia,
                                                   jachar_ut_dia, ia_u_dia, ja_u_dia,
                                                   jachar_u_dia, *pia_lu, *pja_lu,
                                                   *pjachar_lu);
            {
               vector < _Int > ia_dummy;
               vector < _Int > ja_dummy;
               vector < char >jachar_dummy;
               ia_u_dia.swap (ia_dummy);
               ja_u_dia.swap (ja_dummy);
               jachar_u_dia.swap (jachar_dummy);
            }
            {
               vector < _Int > ia_dummy;
               vector < _Int > ja_dummy;
               vector < char >jachar_dummy;
               ia_ut_dia.swap (ia_dummy);
               ja_ut_dia.swap (ja_dummy);
               jachar_ut_dia.swap (jachar_dummy);
            }

            int nzja_lu = (int) ((*pia_lu)[ni_diag]);

            amatr_lu.ResizeList (ni_diag);

            _Int *plist_lu = amatr_lu.GetListArr ();

            int i;

            for (i = 0; i < ni_diag; i++)
               plist_lu[i] = i;

            amatr_lu.SetNlist (ni_diag);
            amatr_lu.SetNzja (nzja_lu);
            amatr_lu.SetNzjaChar (nzja_lu);

// Order

//      ffout << " Point 7 check " << endl;

            CMatrix < _Int, _Flt > amatr_ord;

            amatr_lu.OrderMtrSp (pporder_curr, amatr_ord);

            amatr_lu.Clean ();

            vector < _Int > *pia_ord = amatr_ord.GetIa ();
            vector < _Int > *pja_ord = amatr_ord.GetJa ();
            vector < char >*pjachar_ord = amatr_ord.GetJaChar ();

// Split

//      ffout << " Point 8 check " << endl;

            CFct_impl < _Int, _Flt >::SplitLUSp (ni_diag, *pia_ord, *pja_ord,
                                                 *pjachar_ord, ia_l_dia, ja_l_dia,
                                                 jachar_l_dia, *pia_alu_dia, *pja_alu_dia,
                                                 *pjachar_alu_dia);

// Reorder rows for all off-diagonal blocks

//      ffout << " Point 10 check " << endl;

            {
               int i;
               for (i = 1; i < nzja_curr; i++) {
                  CMatrix < _Int, _Flt > aord;
                  pASub_curr[i].OrderMtrRowsSp (pporder_curr, aord);
                  pASub_curr[i].ReplaceFree (aord);
               }
            }

//      ffout << " Point 11 check " << endl;

// Perform main symbolic fct computations

// Create condense mask array

            int nj_blk = 0;

            int jblk;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               nj_blk += (int) (pblks[jblk + 1] - pblks[jblk]);
            }

            CVectorData < int >imask_col (nj_blk);
            CVectorData < int >list_col (nj_blk);
            CVectorData < int >ibsblk (nblks);

            int *pimask_col = imask_col.Ptr ();
            int *plist_col = list_col.Ptr ();
            int *pibsblk = ibsblk.Ptr ();

            for (i = 0; i < nj_blk; i++)
               pimask_col[i] = -1;

            nj_blk = 0;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               pibsblk[jblk] = nj_blk;
               nj_blk += (int) (pblks[jblk + 1] - pblks[jblk]);
            }

            CVectorData < int >ia_col (nzja_curr + 1);
            int *pia_col = ia_col.Ptr ();

            pia_col[0] = 0;

            int nlist_temp;
            _Int *plist_temp;
            _Int *pia_temp;
            _Int *pja_temp;
            char *pjachar_temp;

            int nlist_col_tot = 0;
            int nzja_tot = 0;

            int ibs, nlist_col, j, jj, k, kk;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               //nj_curr = (int) (pblks[jblk + 1] - pblks[jblk]);
               ibs = pibsblk[jblk];
               nlist_col = 0;
               nlist_temp = pASub_curr[i].GetNlist ();
               pia_temp = pASub_curr[i].GetIaArr ();
               pja_temp = pASub_curr[i].GetJaArr ();
               for (j = 0; j < nlist_temp; j++) {
                  for (k = (int) pia_temp[j]; k < pia_temp[j + 1]; k++) {
                     kk = (int) pja_temp[k];
                     if (pimask_col[ibs + kk] == -1) {
                        pimask_col[ibs + kk] = 1;
                        plist_col[ibs + nlist_col] = kk;
                        nlist_col++;
                     }
                  }
               }
               sort (plist_col + ibs, plist_col + ibs + nlist_col);
               for (j = 0; j < nlist_col; j++) {
                  jj = plist_col[ibs + j];
                  pimask_col[ibs + jj] = nlist_col_tot + j;
               }
               nlist_col_tot += nlist_col;
               pia_col[i + 1] = nlist_col_tot;
               nzja_tot += (int) pia_temp[nlist_temp];
//            ffout << " IndBlk = " << i << " jblk = " << jblk << " nj_curr = " << nj_curr << endl;
//            PrintArray (ffout," ListCol ",nlist_col,plist_col+ibs);
            }

//         PrintArray (ffout," IMaskCol ",nj_blk,pimask_col);

//      ffout << " Point 12 check " << endl;

// Compute condensed block as a sum of blocks

            CMatrix < _Int, _Flt > AU_blk;

            AU_blk.ResizeAndSetAllSp (nlist_col_tot, 0, nzja_tot, 0);

            AU_blk.ResizeJaChar (nzja_tot);
            AU_blk.SetNzjaChar (nzja_tot);

            _Int *plist_AU_blk = AU_blk.GetListArr ();
            _Int *pia_AU_blk = AU_blk.GetIaArr ();
            _Int *pja_AU_blk = AU_blk.GetJaArr ();
            char *pjachar_AU_blk = AU_blk.GetJaCharArr ();

            for (i = 0; i < nlist_col_tot; i++)
               plist_AU_blk[i] = (_Int) i;
            for (i = 0; i <= nlist_col_tot; i++)
               pia_AU_blk[i] = 0;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               ibs = pibsblk[jblk];
               nlist_temp = pASub_curr[i].GetNlist ();
               plist_temp = pASub_curr[i].GetListArr ();
               pia_temp = pASub_curr[i].GetIaArr ();
               for (j = 0; j < nlist_temp; j++) {
                  jj = (int) plist_temp[j];
                  pia_AU_blk[jj + 1] += (_Int) (pia_temp[j + 1] - pia_temp[j]);
               }
            }

            for (i = 0; i < nlist_col_tot; i++)
               pia_AU_blk[i + 1] = pia_AU_blk[i] + pia_AU_blk[i + 1];

            CVectorData < int >iptr (nlist_col_tot);
            int *piptr = iptr.Ptr ();

            for (i = 0; i < nlist_col_tot; i++)
               piptr[i] = (int) pia_AU_blk[i];

            int kind, kknew;

//      ffout << " Point 13 check " << endl;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               ibs = pibsblk[jblk];
               nlist_temp = pASub_curr[i].GetNlist ();
               plist_temp = pASub_curr[i].GetListArr ();
               pia_temp = pASub_curr[i].GetIaArr ();
               pja_temp = pASub_curr[i].GetJaArr ();
               pjachar_temp = pASub_curr[i].GetJaCharArr ();
               for (j = 0; j < nlist_temp; j++) {
                  jj = (int) plist_temp[j];
                  kind = piptr[jj];
                  for (k = (int) pia_temp[j]; k < pia_temp[j + 1]; k++) {
                     kk = (int) pja_temp[k];
                     kknew = pimask_col[ibs + kk];
                     pja_AU_blk[kind] = (_Int) kknew;
                     pjachar_AU_blk[kind] = pjachar_temp[k];
                     kind++;
                  }
                  piptr[jj] = kind;
               }
            }

// Create backward transformation arrays

            CVectorData < int >listcolblk (nlist_col_tot);
            CVectorData < int >indexcol (nlist_col_tot);

            int *plistcolblk = listcolblk.Ptr ();
            int *pindexcol = indexcol.Ptr ();

            nlist_col_tot = 0;

            int njloc;

//      ffout << " Point 14 check " << endl;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               //nj_curr = (int) (pblks[jblk + 1] - pblks[jblk]);
               ibs = pibsblk[jblk];
               njloc = pia_col[i + 1] - pia_col[i];
               for (j = 0; j < njloc; j++) {
                  plistcolblk[nlist_col_tot] = jblk;
                  pindexcol[nlist_col_tot] = plist_col[ibs + j];
                  nlist_col_tot++;
               }
            }

// Perform factorization of the condensed block row

            vector < _Int > *pia_AU_vect = AU_blk.GetIa ();
            vector < _Int > *pja_AU_vect = AU_blk.GetJa ();
            vector < char >*pjachar_AU_vect = AU_blk.GetJaChar ();

            CMatrix < _Int, _Flt > U_blk;

            vector < _Int > *pia_U_vect = U_blk.GetIa ();
            vector < _Int > *pja_U_vect = U_blk.GetJa ();
            vector < char >*pjachar_U_vect = U_blk.GetJaChar ();

//            ffout << " Before fct: AU pair: " << endl;
//            ffout << "   nlist_col_tot = " << nlist_col_tot << " ni_diag = " << ni_diag << endl;
//            ffout << "   tau2_sch = " << tau2_sch << endl;

//            CFct<_Int,_Flt>::PrintMatrix (ffout, AU_blk);

//            char strbuff[256];

//            ofstream *pfout_debug = NULL;

//            sprintf (strbuff,"ChkIluDeg_%i.dat",iblk);

//            ofstream fout_debug (strbuff);

//            pfout_debug = &fout_debug;

//            CFct_impl < _Int, _Flt >::Ilu2BlockIlu2DegreeSp (pfout_debug, pparams->fcttype, pparams->fcttype,
            CFct_impl < _Int, _Flt >::Ilu2BlockIlu2DegreeSp (NULL, pparams->fcttype,
                                                             pparams->fcttype,
                                                             nlist_col_tot, ni_diag,
                                                             *pia_AU_vect, *pja_AU_vect,
                                                             *pjachar_AU_vect,
                                                             *pia_U_vect, *pja_U_vect,
                                                             *pjachar_U_vect);

            U_blk.ResizeList (nlist_col_tot);
            _Int *plist_U = U_blk.GetListArr ();

            for (i = 0; i < nlist_col_tot; i++)
               plist_U[i] = i;

            int nzja_U_temp = (int) (*pia_U_vect)[nlist_col_tot];

            U_blk.SetNlist (nlist_col_tot);
            U_blk.SetNzja (nzja_U_temp);
            U_blk.SetNzjaChar (nzja_U_temp);

//            ffout << " After Fct: U pair: " << endl;
//            CFct<_Int,_Flt>::PrintMatrix (ffout, U_blk);

// Split U data

            int nlist_U = nlist_col_tot;
            _Int *pia_U = U_blk.GetIaArr ();
            _Int *pja_U = U_blk.GetJaArr ();
            char *pjachar_U = U_blk.GetJaCharArr ();

            if (nlist_U > 0) {

               icycleblk++;

               int nzjablk_U = 0;

               int kblk;

               for (j = 0; j < ni_diag; j++) {
                  for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                     kk = (int) pja_U[k];
                     kblk = plistcolblk[kk];
                     if (pimaskblk[kblk] != icycleblk) {
                        pimaskblk[kblk] = icycleblk;
                        plistblk[nzjablk_U] = kblk;
                        nzjablk_U++;
                     }
                  }
               }

               sort (plistblk, plistblk + nzjablk_U);

               pblockrowsLU_arr[iblk].ResizeASub (nzjablk_U);
               pblockrowsLU_arr[iblk].SetNzblk (nzjablk_U);

               CMatrix < _Int, _Flt > *pASub_U = pblockrowsLU_arr[iblk].GetASubArr ();

               CMatrix < int, float >*pHMatrU = pblockrowsLU_arr[iblk].GetHMatrStr ();

               pHMatrU->ResizeAndSetAllSp (1, 1, nzjablk_U, nzjablk_U);

               int *plist_UHblk = pHMatrU->GetListArr ();
               int *plist2_UHblk = pHMatrU->GetList2Arr ();
               int *pia_UHblk = pHMatrU->GetIaArr ();
               int *pja_UHblk = pHMatrU->GetJaArr ();
               int *pja2_UHblk = pHMatrU->GetJa2Arr ();

               plist_UHblk[0] = iblk;
               plist2_UHblk[0] = 0;
               pia_UHblk[0] = 0;
               pia_UHblk[1] = nzjablk_U;
               for (j = 0; j < nzjablk_U; j++)
                  pja_UHblk[j] = plistblk[j];
               for (j = 0; j < nzjablk_U; j++)
                  pja2_UHblk[j] = 0;

               for (j = 0; j < nzjablk_U; j++) {
                  jblk = plistblk[j];
                  pindblk[jblk] = j;
               }

               CVectorData < int >nzja_blkarr (nzjablk_U);
               int *pnzja_blkarr = nzja_blkarr.Ptr ();

               for (j = 0; j < nzjablk_U; j++)
                  pnzja_blkarr[j] = 0;

               int ind;

               for (j = 0; j < ni_diag; j++) {
                  for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                     kk = (int) pja_U[k];
                     kblk = plistcolblk[kk];
                     ind = pindblk[kblk];
                     pnzja_blkarr[ind]++;
                  }
               }

               CVectorData < CVectorData < _Int > >irows (nzjablk_U);
               CVectorData < CVectorData < _Int > >icols (nzjablk_U);
               CVectorData < CVectorData < char > >chars (nzjablk_U);

               CVectorData < _Int > *pirows = irows.Ptr ();
               CVectorData < _Int > *picols = icols.Ptr ();
               CVectorData < char >*pchars = chars.Ptr ();

               for (i = 0; i < nzjablk_U; i++) {
                  pirows[i].resize (pnzja_blkarr[i]);
                  picols[i].resize (pnzja_blkarr[i]);
                  pchars[i].resize (pnzja_blkarr[i]);
               }

               for (j = 0; j < nzjablk_U; j++)
                  pnzja_blkarr[j] = 0;

               int kkk;

               for (j = 0; j < ni_diag; j++) {
                  for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                     kk = (int) pja_U[k];
                     kblk = plistcolblk[kk];
                     ind = pindblk[kblk];
                     _Int *ppirows = pirows[ind].Ptr ();
                     _Int *ppicols = picols[ind].Ptr ();
                     char *ppchars = pchars[ind].Ptr ();
                     kkk = pnzja_blkarr[ind];
                     ppirows[kkk] = (_Int) j;
                     ppicols[kkk] = (_Int) pindexcol[kk];
                     ppchars[kkk] = pjachar_U[k];
                     pnzja_blkarr[ind]++;
                  }
               }

               int njrow, ibeg;

               for (j = 0; j < nzjablk_U; j++) {
                  icycleblk++;
                  nlist_temp = 0;
                  _Int *ppirows = pirows[j].Ptr ();
                  _Int *ppicols = picols[j].Ptr ();
                  char *ppchars = pchars[j].Ptr ();
                  for (k = 0; k < pnzja_blkarr[j]; k++) {
                     kk = (int) ppirows[k];
                     if (pimask[kk] != icycleblk) {
                        pimask[kk] = icycleblk;
                        plist[nlist_temp] = kk;
                        nlist_temp++;
                     }
                  }
                  sort (plist, plist + nlist_temp);
                  pASub_U[j].ResizeAndSetAllSp (nlist_temp, 0, pnzja_blkarr[j], 0);
                  pASub_U[j].ResizeJaChar (pnzja_blkarr[j]);
                  pASub_U[j].SetNzjaChar (pnzja_blkarr[j]);
                  _Int *plist_temp = pASub_U[j].GetListArr ();
                  _Int *pia_temp = pASub_U[j].GetIaArr ();
                  _Int *pja_temp = pASub_U[j].GetJaArr ();
                  char *pjachar_temp = pASub_U[j].GetJaCharArr ();
                  for (k = 0; k < nlist_temp; k++)
                     plist_temp[k] = (_Int) plist[k];
                  for (k = 0; k < nlist_temp; k++) {
                     kk = (int) plist_temp[k];
                     plist[kk] = k;
                  }
                  for (k = 0; k <= nlist_temp; k++)
                     pia_temp[k] = 0;
                  for (k = 0; k < pnzja_blkarr[j]; k++) {
                     kk = (int) ppirows[k];
                     ind = plist[kk];
                     pia_temp[ind + 1]++;
                  }
                  for (k = 0; k < nlist_temp; k++)
                     pia_temp[k + 1] = pia_temp[k] + pia_temp[k + 1];
                  CVectorData < int >iptr (nlist_temp);
                  int *piptr = iptr.Ptr ();
                  for (k = 0; k < nlist_temp; k++)
                     piptr[k] = (int) pia_temp[k];
                  for (k = 0; k < pnzja_blkarr[j]; k++) {
                     kk = (int) ppirows[k];
                     ind = plist[kk];
                     kkk = piptr[ind];
                     pja_temp[kkk] = ppicols[k];
                     pjachar_temp[kkk] = ppchars[k];
                     piptr[ind]++;
                  }
                  int njrowmax = 0;
                  for (k = 0; k < nlist_temp; k++) {
                     njrow = (int) (pia_temp[k + 1] - pia_temp[k]);
                     if (njrow > njrowmax)
                        njrowmax = njrow;
                  }
                  CVectorData < CSortInt > iiarr (njrowmax);
                  CVectorData < char >chararr (njrowmax);
                  CSortInt *piiarr = iiarr.Ptr ();
                  char *pchararr = chararr.Ptr ();
                  for (k = 0; k < nlist_temp; k++) {
                     ibeg = (int) pia_temp[k];
                     njrow = (int) (pia_temp[k + 1] - pia_temp[k]);
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        piiarr[kk - ibeg].ival = (int) pja_temp[kk];
                        piiarr[kk - ibeg].i2val = kk;
                        pchararr[kk - ibeg] = pjachar_temp[kk];
                     }
                     sort (piiarr, piiarr + njrow);
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        kkk = piiarr[kk - ibeg].ival;
                        pja_temp[kk] = (_Int) kkk;
                        kkk = piiarr[kk - ibeg].i2val;
                        pjachar_temp[kk] = pchararr[kkk - ibeg];
                     }
                  }
               }

            } else {

               pblockrowsLU_arr[iblk].ResizeASub (1);
               pblockrowsLU_arr[iblk].SetNzblk (1);

               CMatrix < int, float >*pHMatrU = pblockrowsLU_arr[iblk].GetHMatrStr ();

               pHMatrU->ResizeAndSetAllSp (1, 1, 1, 1);

               int *plist_UHblk = pHMatrU->GetListArr ();
               int *plist2_UHblk = pHMatrU->GetList2Arr ();
               int *pia_UHblk = pHMatrU->GetIaArr ();
               int *pja_UHblk = pHMatrU->GetJaArr ();
               int *pja2_UHblk = pHMatrU->GetJa2Arr ();

               plist_UHblk[0] = iblk;
               plist2_UHblk[0] = 0;
               pia_UHblk[0] = 0;
               pia_UHblk[1] = 1;
               pja_UHblk[0] = iblk;
               pja2_UHblk[0] = 0;

            }

//            ffout << " U hblk fct: " << endl;
//            pblockrowsLU_arr[iblk].PrintHMatrix (ffout);

// Split Schur complement data

            {

               int nlistblk_schur = nzja_curr - 1;
               int nzjablk_schur = 0;

               int kblk, nlistblk_temp;

               for (i = 1; i < nzja_curr; i++) {
                  jblk = pja_curr[i];
                  icycleblk++;
                  njloc = pia_col[i + 1] - pia_col[i];
                  nlistblk_temp = 0;
                  for (j = pia_col[i]; j < pia_col[i + 1]; j++) {
                     for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                        kk = (int) pja_U[k];
                        kblk = plistcolblk[kk];
                        if (pimaskblk[kblk] != icycleblk) {
                           pimaskblk[kblk] = icycleblk;
                           nlistblk_temp++;
                        }
                     }
                  }
                  nzjablk_schur += nlistblk_temp;
               }

               pschur_arr[inode].ResizeASub (nzjablk_schur);
               pschur_arr[inode].SetNzblk (nzjablk_schur);

               CMatrix < _Int, _Flt > *pASub_Schur = pschur_arr[inode].GetASubArr ();

               CMatrix < int, float >*pHMatrSchur = pschur_arr[inode].GetHMatrStr ();

               pHMatrSchur->ResizeAndSetAllSp (nlistblk_schur, nlistblk_schur,
                                               nzjablk_schur, nzjablk_schur);

               int *plist_SchurHblk = pHMatrSchur->GetListArr ();
               int *plist2_SchurHblk = pHMatrSchur->GetList2Arr ();
               int *pia_SchurHblk = pHMatrSchur->GetIaArr ();
               int *pja_SchurHblk = pHMatrSchur->GetJaArr ();
               int *pja2_SchurHblk = pHMatrSchur->GetJa2Arr ();

               int iloc;

               for (i = 1; i < nzja_curr; i++) {
                  iloc = i - 1;
                  jblk = pja_curr[i];
                  plist_SchurHblk[iloc] = jblk;
                  plist2_SchurHblk[iloc] = 0;
               }

               pia_SchurHblk[0] = 0;
               nzjablk_schur = 0;

               for (i = 1; i < nzja_curr; i++) {
                  icycleblk++;
                  iloc = i;
                  jblk = pja_curr[i];
                  njloc = pia_col[iloc + 1] - pia_col[iloc];
                  nlistblk_temp = 0;
                  for (j = pia_col[iloc]; j < pia_col[iloc + 1]; j++) {
                     for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                        kk = (int) pja_U[k];
                        kblk = plistcolblk[kk];
                        if (pimaskblk[kblk] != icycleblk) {
                           pimaskblk[kblk] = icycleblk;
                           plistblk[nlistblk_temp] = kblk;
                           nlistblk_temp++;
                        }
                     }
                  }
                  sort (plistblk, plistblk + nlistblk_temp);
                  for (j = 0; j < nlistblk_temp; j++)
                     pja_SchurHblk[nzjablk_schur + j] = plistblk[j];
                  for (j = 0; j < nlistblk_temp; j++)
                     pja2_SchurHblk[nzjablk_schur + j] = 0;
                  nzjablk_schur += nlistblk_temp;
                  pia_SchurHblk[iloc] = nzjablk_schur;
               }

//            ffout << " Schur str: " << endl;
//            PrintArray (ffout," plist_SchurHblk = ",nlistblk_schur,plist_SchurHblk);
//            PrintArray (ffout," plist2_SchurHblk = ",nlistblk_schur,plist2_SchurHblk);
//            PrintArray (ffout," pia_SchurHblk = ",nlistblk_schur+1,pia_SchurHblk);
//            PrintArray (ffout," pja_SchurHblk = ",nzjablk_schur,pja_SchurHblk);
//            PrintArray (ffout," pja2_SchurHblk = ",nzjablk_schur,pja2_SchurHblk);

               CVectorData < int >nzja_blkarr (nzjablk_schur);
               int *pnzja_blkarr = nzja_blkarr.Ptr ();

               for (j = 0; j < nzjablk_schur; j++)
                  pnzja_blkarr[j] = 0;

               int ind;

               for (i = 1; i < nzja_curr; i++) {
                  iloc = i;
                  for (j = pia_SchurHblk[iloc - 1]; j < pia_SchurHblk[iloc]; j++) {
                     kblk = pja_SchurHblk[j];
                     pindblk[kblk] = j;
                  }
                  for (j = pia_col[iloc]; j < pia_col[iloc + 1]; j++) {
                     for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                        kk = (int) pja_U[k];
                        kblk = plistcolblk[kk];
                        ind = pindblk[kblk];
                        pnzja_blkarr[ind]++;
                     }
                  }
               }

               CVectorData < CVectorData < _Int > >irows (nzjablk_schur);
               CVectorData < CVectorData < _Int > >icols (nzjablk_schur);
               CVectorData < CVectorData < char > >chars (nzjablk_schur);

               CVectorData < _Int > *pirows = irows.Ptr ();
               CVectorData < _Int > *picols = icols.Ptr ();
               CVectorData < char >*pchars = chars.Ptr ();

               for (i = 0; i < nzjablk_schur; i++) {
                  pirows[i].resize (pnzja_blkarr[i]);
                  picols[i].resize (pnzja_blkarr[i]);
                  pchars[i].resize (pnzja_blkarr[i]);
               }

               for (j = 0; j < nzjablk_schur; j++)
                  pnzja_blkarr[j] = 0;

               int kkk;

               for (i = 1; i < nzja_curr; i++) {
                  iloc = i;
                  for (j = pia_SchurHblk[iloc - 1]; j < pia_SchurHblk[iloc]; j++) {
                     kblk = pja_SchurHblk[j];
                     pindblk[kblk] = j;
                  }
                  for (j = pia_col[iloc]; j < pia_col[iloc + 1]; j++) {
                     for (k = (int) pia_U[j]; k < pia_U[j + 1]; k++) {
                        kk = (int) pja_U[k];
                        kblk = plistcolblk[kk];
                        ind = pindblk[kblk];
                        _Int *ppirows = pirows[ind].Ptr ();
                        _Int *ppicols = picols[ind].Ptr ();
                        char *ppchars = pchars[ind].Ptr ();
                        kkk = pnzja_blkarr[ind];
                        ppirows[kkk] = (_Int) pindexcol[j];
                        ppicols[kkk] = (_Int) pindexcol[kk];
                        ppchars[kkk] = pjachar_U[k];
                        pnzja_blkarr[ind]++;
                     }
                  }
               }

               int njrow, ibeg;

               for (j = 0; j < nzjablk_schur; j++) {
                  icycleblk++;
                  nlist_temp = 0;
                  _Int *ppirows = pirows[j].Ptr ();
                  _Int *ppicols = picols[j].Ptr ();
                  char *ppchars = pchars[j].Ptr ();
                  for (k = 0; k < pnzja_blkarr[j]; k++) {
                     kk = (int) ppirows[k];
                     if (pimask[kk] != icycleblk) {
                        pimask[kk] = icycleblk;
                        plist[nlist_temp] = kk;
                        nlist_temp++;
                     }
                  }
                  sort (plist, plist + nlist_temp);
                  pASub_Schur[j].ResizeAndSetAllSp (nlist_temp, 0, pnzja_blkarr[j], 0);
                  pASub_Schur[j].ResizeJaChar (pnzja_blkarr[j]);
                  pASub_Schur[j].SetNzjaChar (pnzja_blkarr[j]);
                  _Int *plist_temp = pASub_Schur[j].GetListArr ();
                  _Int *pia_temp = pASub_Schur[j].GetIaArr ();
                  _Int *pja_temp = pASub_Schur[j].GetJaArr ();
                  char *pjachar_temp = pASub_Schur[j].GetJaCharArr ();
                  for (k = 0; k < nlist_temp; k++)
                     plist_temp[k] = (_Int) plist[k];
                  for (k = 0; k < nlist_temp; k++) {
                     kk = (int) plist_temp[k];
                     plist[kk] = k;
                  }
                  for (k = 0; k <= nlist_temp; k++)
                     pia_temp[k] = 0;
                  for (k = 0; k < pnzja_blkarr[j]; k++) {
                     kk = (int) ppirows[k];
                     ind = plist[kk];
                     pia_temp[ind + 1]++;
                  }
                  for (k = 0; k < nlist_temp; k++)
                     pia_temp[k + 1] = pia_temp[k] + pia_temp[k + 1];
                  CVectorData < int >iptr (nlist_temp);
                  int *piptr = iptr.Ptr ();
                  for (k = 0; k < nlist_temp; k++)
                     piptr[k] = (int) pia_temp[k];
                  for (k = 0; k < pnzja_blkarr[j]; k++) {
                     kk = (int) ppirows[k];
                     ind = plist[kk];
                     kkk = piptr[ind];
                     pja_temp[kkk] = ppicols[k];
                     pjachar_temp[kkk] = ppchars[k];
                     piptr[ind]++;
                  }
                  int njrowmax = 0;
                  for (k = 0; k < nlist_temp; k++) {
                     njrow = (int) (pia_temp[k + 1] - pia_temp[k]);
                     if (njrow > njrowmax)
                        njrowmax = njrow;
                  }
                  CVectorData < CSortInt > iiarr (njrowmax);
                  CVectorData < char >charsort (njrowmax);
                  CSortInt *piiarr = iiarr.Ptr ();
                  char *pcharsort = charsort.Ptr ();

                  for (k = 0; k < nlist_temp; k++) {
                     ibeg = (int) pia_temp[k];
                     njrow = (int) (pia_temp[k + 1] - pia_temp[k]);
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        piiarr[kk - ibeg].ival = (int) pja_temp[kk];
                        piiarr[kk - ibeg].i2val = kk;
                     }
                     sort (piiarr, piiarr + njrow);
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        kkk = piiarr[kk - ibeg].i2val;
                        pcharsort[kk - ibeg] = pjachar_temp[kkk];
                     }
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        kkk = piiarr[kk - ibeg].ival;
                        pja_temp[kk] = (_Int) kkk;
                        pjachar_temp[kk] = pcharsort[kk - ibeg];
                     }

                  }
               }

            }

//            ffout << " =========== Schur hblk fct: " << endl;
//            pschur_arr[inode].PrintHMatrix (ffout);
         }

// Finally modify Schur data

         pschur_arr[inode].AddReplaceSpThr (false, hblk_sum);

//         ffout << " =========== Modified Schur hblk fct: " << endl;
//         pschur_arr[iblk].PrintHMatrix (ffout);

         picycle_thr[my_thr] = icycleblk;

      }
   };

//
// Perform symbolic factorization of the hmatrix with dynamic ordering and subtrees for schur blocks
//========================================================================================
   template < typename _Int, typename _Flt > void CFctThreads < _Int,
      _Flt >::SymbolicFct (ofstream * _pfout, SParams * _params, CTree & _tree,
                           int _nblks, long long *_blks, CBMatrix < _Int, _Flt > &_alu_sp,
                           CVectorData < CTree > &_subtree_arr,
                           CVectorData < int >&_nblks_arr,
                           CVectorData < vector < int > >&_blks_arr,
                           CVectorData < int >&_ordernew, CBMatrix < _Int, _Flt > &_lu_sp)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

//      char strbuff[256];
//      sprintf (strbuff,"ChkSymbFct_%i.dat",_index_out);
//      ofstream ffout (strbuff);
//      _alu_sp.PrintHMatrix (ffout);

// Open tree and check

      int nnodes_tree = _tree.GetNnodes ();
      int *pnode2ind = _tree.GetNode2Ind ();
      int *pind2node = _tree.GetInd2Node ();

      {
         int nblks_tree = 0;
         int i, ind;
         for (i = 0; i < nnodes_tree; i++) {
            ind = pnode2ind[i];
            if (ind >= 0) {
               if (ind != nblks_tree) {
                  cout <<
                     " CFctThreads<>::SymbolicFct: error: incorrect index number in a tree "
                     << endl;
                  throw
                     " CFctThreads<>::SymbolicFct: error: incorrect index number in a tree ";
               }
               if (pind2node[ind] != i) {
                  cout <<
                     " CFctThreads<>::SymbolicFct: error: incorrect node number for index in a tree "
                     << endl;
                  throw
                     " CFctThreads<>::SymbolicFct: error: incorrect node number for index in a tree ";
               }
               nblks_tree++;
            }
         }
         if (nblks_tree != _nblks) {
            cout <<
               " CFctThreads<>::SymbolicFct: error: incorrect number of blocks in a tree "
               << endl;
            throw
               " CFctThreads<>::SymbolicFct: error: incorrect number of blocks in a tree ";
         }
      }

// Open hmatrix

      CMatrix < int, float >*phmatr_sp = _alu_sp.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_sp = _alu_sp.GetASubArr ();

      int *pia_hmatr_sp = phmatr_sp->GetIaArr ();
      int *pja_hmatr_sp = phmatr_sp->GetJaArr ();

// Allocate output arrays

      _subtree_arr.resize (nnodes_tree);
      _nblks_arr.resize (nnodes_tree);
      _blks_arr.resize (nnodes_tree);

      CTree *p_subtree_arr = _subtree_arr.Ptr ();
      int *p_nblks_arr = _nblks_arr.Ptr ();
      vector < int >*p_blks_arr = _blks_arr.Ptr ();

// Allocate all necessary work arrays

      vector < vector < int > >order_arr (nnodes_tree);
      vector < CBMatrix < _Int, _Flt > >blockrowsLU_arr (nnodes_tree + 1);
      vector < CBMatrix < _Int, _Flt > >schur_arr (nnodes_tree + 1);

      vector < int >*porder_arr = order_arr.data ();
      CBMatrix < _Int, _Flt > *pblockrowsLU_arr = blockrowsLU_arr.data ();
      CBMatrix < _Int, _Flt > *pschur_arr = schur_arr.data ();

      int nimax = 0;
      int niloc;

      int i;

      for (i = 0; i < _nblks; i++) {
         niloc = (int) (_blks[i + 1] - _blks[i]);
         if (niloc > nimax)
            nimax = niloc;
      }

      vector < int >icycle_thr (n_thr + 1);
      vector < CVectorData < int > >imask_thr (n_thr + 1);

      int *picycle_thr = &icycle_thr[0];
      CVectorData < int >*pimask_thr = &imask_thr[0];

      for (i = 0; i < n_thr; i++)
         picycle_thr[i] = -1;

// Perform symbolic factorization of the tree level by level

      vector < CBlockSymbFctTree < _Int, _Flt > >fct (nnodes_tree + 1);
      CBlockSymbFctTree < _Int, _Flt > *pfct = &fct[0];

      for (i = 0; i < nnodes_tree; i++) {

         CBlockSymbFctTree < _Int, _Flt > *task = pfct + i;

         task->nblks = _nblks;
         task->pblks = _blks;
         task->inode = i;
         task->ptree = &_tree;
         task->pparams = _params;
         task->pia_ALU_sp = pia_hmatr_sp;
         task->pja_ALU_sp = pja_hmatr_sp;
         task->pASub_ALU_sp = pA_sub_sp;
         task->p_subtree_arr = p_subtree_arr;
         task->p_nblks_arr = p_nblks_arr;
         task->p_blks_arr = p_blks_arr;
         task->porder_arr = porder_arr;
         task->pblockrowsLU_arr = pblockrowsLU_arr;
         task->pschur_arr = pschur_arr;
         task->nimax = nimax;
         task->picycle_thr = picycle_thr;
         task->pimask_thr = pimask_thr;

      }

// Perform computations according to two tree's

      int nlev = _tree.GetNlev ();
      int *pnnodes_lev = _tree.GetNNodesLev ();
      vector < int >*pnodeslevlist = _tree.GetNodesLevList ();

// Run

      int ilev, nnodes_curr;
      int *ppnodeslevlist;

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 2 " << endl;

      for (ilev = nlev - 1; ilev >= 0; ilev--) {

         nnodes_curr = pnnodes_lev[ilev];
         ppnodeslevlist = &pnodeslevlist[ilev][0];

//      if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 2.1 Ilev = " << ilev << " nnodes_curr_1 = " << nnodes_curr_1 << endl;

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr; ipar++) {
            int inode = ppnodeslevlist[ipar];
            pfct[inode].Execute ();
         }

//      if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 2.2 After " << endl;

      }

// Combine data into one hblock

      int nzja_fin = 0;

      int j, jj;

      for (i = 0; i < _nblks; i++) {
         CMatrix < int, float >*phmatr_temp = pblockrowsLU_arr[i].GetHMatrStr ();
         int nzja_temp = phmatr_temp->GetNzja ();
         nzja_fin += nzja_temp;
      }

      _lu_sp.SetNzblk (nzja_fin);
      _lu_sp.ResizeASub (nzja_fin);

      CMatrix < int, float >*phmatr_fin = _lu_sp.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pASub_fin = _lu_sp.GetASubArr ();

      phmatr_fin->ResizeAndSetAllSp (_nblks, _nblks, nzja_fin, nzja_fin);

      int *plist_fin = phmatr_fin->GetListArr ();
      int *plist2_fin = phmatr_fin->GetList2Arr ();
      int *pia_fin = phmatr_fin->GetIaArr ();
      int *pja_fin = phmatr_fin->GetJaArr ();
      int *pja2_fin = phmatr_fin->GetJa2Arr ();

      for (i = 0; i < _nblks; i++)
         plist_fin[i] = i;
      for (i = 0; i < _nblks; i++)
         plist2_fin[i] = 0;
      for (i = 0; i < nzja_fin; i++)
         pja2_fin[i] = 0;

      nzja_fin = 0;
      pia_fin[0] = 0;

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 13 " << endl;

      for (i = 0; i < _nblks; i++) {

         CMatrix < int, float >*phmatr_temp = pblockrowsLU_arr[i].GetHMatrStr ();
         CMatrix < _Int, _Flt > *pASub_temp = pblockrowsLU_arr[i].GetASubArr ();

         int *pia_temp = phmatr_temp->GetIaArr ();
         int *pja_temp = phmatr_temp->GetJaArr ();

         int nzja_temp = pia_temp[1];

         for (j = 0; j < nzja_temp; j++) {
            jj = pja_temp[j];
            pja_fin[nzja_fin] = jj;
            pASub_fin[nzja_fin].ReplaceFree (pASub_temp[j]);
            nzja_fin++;
         }

         pia_fin[i + 1] = nzja_fin;

      }

// Reorder off diagonal data

      {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < _nblks; ipar++) {

            int j, jj;
            for (j = pia_fin[ipar]; j < pia_fin[ipar + 1]; j++) {
               jj = pja_fin[j];
               if (jj != ipar) {
                  int *porder_temp = &porder_arr[jj][0];

                  CMatrix < _Int, _Flt > aord;

                  pASub_fin[j].OrderMtrColsSp (porder_temp, aord);

                  pASub_fin[j].ReplaceFree (aord);
               }
            }
         }
      }

// Form global ordering array

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: Point 14 " << endl;

      int ntot = (int) _blks[_nblks];

      _ordernew.resize (ntot);
      int *p_ordernew = _ordernew.Ptr ();

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < _nblks; ipar++) {
            int ibeg, niloc, j;
            niloc = (int) (_blks[ipar + 1] - _blks[ipar]);
            ibeg = (int) _blks[ipar];
            int *pporder_temp = &porder_arr[ipar][0];
            for (j = 0; j < niloc; j++)
               p_ordernew[ibeg + j] = ibeg + pporder_temp[j];
         }
      }

//   if (_pfout != NULL) *_pfout << " Ilu2BlockIlu2: End point " << endl;
//      ffout << " Fct Result : " << endl;
//      _lu_pair.PrintHMatrix (ffout);

   }

// Compute ND order, create binary tree, find separators and condense tree, and compute subtree's for nodes
//========================================================================================
   template < typename _Int, typename _Flt > void CFctThreads < _Int,
      _Flt >::OrderNDSepSubTree (int _nlev, int _nlev_split, int _ordtype, int _fcttype,
                                 CMatrix < _Int, _Flt > &_amatr_sp, vector < int >&_order,
                                 CTree & _tree, int &_nhblks, vector < int >&_hblks,
                                 CVectorData < CTree > &_subtree_arr,
                                 vector < int >&_hblk2blks, vector < int >&_blks)
   {

// Compute major tree and ND ordering

      int ntot = _amatr_sp.GetNlist ();

      vector < int >order;

      _amatr_sp.OrderNDSeparatorsForTree (_nlev, -1, order, _tree, _nhblks, _hblks);

      int *porder = &order[0];

      int nnodes_tree = _tree.GetNnodes ();
      int *pblks = &_hblks[0];

      CMatrix < _Int, _Flt > amatr_ord_sp;

      _amatr_sp.OrderMtrSp (porder, amatr_ord_sp);

      if (true) {
         CBMatrix < _Int, _Flt >::Str2PsBox (3, amatr_ord_sp, "StrOrd_ini.ps", _nhblks,
                                             pblks);
      }
// Split sparsity into L and U

      CMatrix < _Int, _Flt > al_sp;
      CMatrix < _Int, _Flt > au_sp;

      amatr_ord_sp.SplitLUSp (al_sp, au_sp);

      amatr_ord_sp.Clean ();
      al_sp.Clean ();

// Split matrix into submatrices

      CVectorData < long long >blks_64 (_nhblks + 1);
      long long *pblks_64 = blks_64.Ptr ();

      int i;

      for (i = 0; i <= _nhblks; i++)
         pblks_64[i] = pblks[i];

      CBMatrix < _Int, _Flt > hmatr_ALU_sp;

      CBMatrix < _Int, _Flt >::SplitMatrSpIntoHMatrSp (_nhblks, pblks_64, au_sp,
                                                       hmatr_ALU_sp);

// Perform parallel symbolic fct

      CVectorData < int >nblks_arr;
      CVectorData < vector < int > >blks_arr;

      CVectorData < int >ordernew;

      CBMatrix < _Int, _Flt > hmatr_U_sp;

      SParams params;

      params.ordtype = _ordtype;
      params.nlev_split = _nlev_split;
      params.fcttype = _fcttype;

      CFctThreads < _Int, _Flt >::SymbolicFct (NULL, &params, _tree, _nhblks, pblks_64,
                                               hmatr_ALU_sp, _subtree_arr, nblks_arr,
                                               blks_arr, ordernew, hmatr_U_sp);

      int *pnblks_arr = nblks_arr.Ptr ();
      vector < int >*pblks_arr = blks_arr.Ptr ();
      int *pordernew = ordernew.Ptr ();

// Compute extended block partitioning

      _hblk2blks.resize (_nhblks + 1);
      int *phblk2blks = &_hblk2blks[0];

      int nblks_ext = 0;

      phblk2blks[0] = 0;

      for (i = 0; i < _nhblks; i++) {
         nblks_ext += pnblks_arr[i];
         phblk2blks[i + 1] = nblks_ext;
      }

      _blks.resize (nblks_ext + 1);
      int *pblks_ext = &_blks[0];

      nblks_ext = 0;
      pblks_ext[0] = 0;

      int j;

      for (i = 0; i < _nhblks; i++) {
         int *ppblks_arr = &pblks_arr[i][0];
         for (j = 0; j < pnblks_arr[i]; j++) {
            pblks_ext[nblks_ext + 1] =
               pblks_ext[nblks_ext] + (ppblks_arr[j + 1] - ppblks_arr[j]);
            nblks_ext++;
         }
      }

// Compute final order

      CVectorData < int >iorder (ntot);
      int *piorder = iorder.Ptr ();

      _order.resize (ntot + 1);
      int *p_order = &_order[0];

      for (i = 0; i < ntot; i++)
         piorder[porder[i]] = i;

      int iold;

      for (i = 0; i < ntot; i++) {
         iold = piorder[i];
         p_order[iold] = pordernew[i];
      }

   }

// Rescale U
//========================================================================================
   template < typename _Int, typename _Flt > void CFctThreads < _Int,
      _Flt >::RescaleU (int _nblks, long long *_blks, CBMatrix < _Int, _Flt > &_U_matr,
                        _Flt * _sclU, _Flt * _inv_sclU)
   {

// Open hmatrix

      CMatrix < int, float >*phmatr_U = _U_matr.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_U = _U_matr.GetASubArr ();

      int *pia_hmatr = phmatr_U->GetIaArr ();
      int *pja_hmatr = phmatr_U->GetJaArr ();

// Rescale data

#ifdef USE_THREADS
#pragma omp parallel for
#endif
      for (int ipar = 0; ipar < _nblks; ipar++) {

         int j, jj, ki, kj, kk;
         int ibeg = (int) _blks[ipar];

         for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
            jj = pja_hmatr[j];
            int jbeg = (int) _blks[jj];
            if (jj == ipar) {
               int nlist_temp = pA_sub_U[j].GetNlist ();
               _Int *pia_temp = pA_sub_U[j].GetIaArr ();
               _Int *pja_temp = pA_sub_U[j].GetJaArr ();
               _Flt *pa_temp = pA_sub_U[j].GetAArr ();
               for (ki = 0; ki < nlist_temp; ki++) {
                  for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                     kk = (int) pja_temp[kj];
                     if (kk != ki) {
                        pa_temp[kj] *= _inv_sclU[jbeg + kk];
                     } else {
                        pa_temp[kj] *= _sclU[ibeg + ki];
                     }
                  }
               }
            } else {
               int nlist_temp = pA_sub_U[j].GetNlist ();
               _Int *pia_temp = pA_sub_U[j].GetIaArr ();
               _Int *pja_temp = pA_sub_U[j].GetJaArr ();
               _Flt *pa_temp = pA_sub_U[j].GetAArr ();
               for (ki = 0; ki < nlist_temp; ki++) {
                  for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                     kk = (int) pja_temp[kj];
                     pa_temp[kj] *= _inv_sclU[jbeg + kk];
                  }
               }
            }
         }
      }

   }

// Set diagonal by zeroes
//========================================================================================
   template < typename _Int, typename _Flt > void CFctThreads < _Int,
      _Flt >::SetZeroDiag (int _nblks, long long *_blks, CBMatrix < _Int, _Flt > &_A_matr)
   {

// Open hmatrix

      CMatrix < int, float >*phmatr_A = _A_matr.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_A = _A_matr.GetASubArr ();

      int *pia_hmatr = phmatr_A->GetIaArr ();
      int *pja_hmatr = phmatr_A->GetJaArr ();

// Set diag values by zeroes

#ifdef USE_THREADS
#pragma omp parallel for
#endif
      for (int ipar = 0; ipar < _nblks; ipar++) {

         int j, jj, ki, irow, kj, kk;
         _Flt fzero = (_Flt) 0.0e0;

         for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
            jj = pja_hmatr[j];
            if (jj == ipar) {
               int nlist_temp = pA_sub_A[j].GetNlist ();
               _Int *plist_temp = pA_sub_A[j].GetListArr ();
               _Int *pia_temp = pA_sub_A[j].GetIaArr ();
               _Int *pja_temp = pA_sub_A[j].GetJaArr ();
               _Flt *pa_temp = pA_sub_A[j].GetAArr ();
               for (ki = 0; ki < nlist_temp; ki++) {
                  irow = (int) plist_temp[ki];
                  for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                     kk = (int) pja_temp[kj];
                     if (kk == irow) {
                        pa_temp[kj] = fzero;
                     }
                  }
               }
            }
         }
      }

   }

// Perform filtering of the pairs block sparsity according to the tree with diagonal modifications
//========================================================================================
   template < typename _Int, typename _Flt > void CFctThreads < _Int,
      _Flt >::PairsFilterTreeModif (CTree & _tree, double _theta, int _nblks,
                                    long long *_blks, CBMatrix < _Int, _Flt > &_AU_matr)
   {

// Open hmatrix

      CMatrix < int, float >*phmatr_AU = _AU_matr.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_AU = _AU_matr.GetASubArr ();

      int nzja_hmatr = phmatr_AU->GetNzja ();
      int *pia_hmatr = phmatr_AU->GetIaArr ();
      int *pja_hmatr = phmatr_AU->GetJaArr ();

// Compute transposed block sparsity

      CMatrix < int, float >hmatrT_str;
      vector < int >indHt2H (nzja_hmatr + 1);

      int *pindHt2H = &indHt2H[0];

      {

         CVectorData < int >imask (5 * _nblks + 1);
         int *pimask = imask.Ptr ();

         int i;

         for (i = 0; i < _nblks; i++)
            pimask[i] = -1;

         int icycle = -1;

         phmatr_AU->TransposedSparsityListSp (icycle, pimask, pimask + _nblks,
                                              pimask + 2 * _nblks, pimask + 3 * _nblks,
                                              pimask + 4 * _nblks, hmatrT_str);

         int *pia_hmatrT = hmatrT_str.GetIaArr ();
         int *pja_hmatrT = hmatrT_str.GetJaArr ();

         for (i = 0; i < _nblks; i++)
            pimask[i] = pia_hmatr[i];

         int j, jj, k;

         for (i = 0; i < _nblks; i++) {
            for (j = pia_hmatrT[i]; j < pia_hmatrT[i + 1]; j++) {
               jj = pja_hmatrT[j];
               k = pimask[jj];
               pindHt2H[j] = k;
               pimask[jj]++;
            }
         }

      }

      int *pia_hmatrT = hmatrT_str.GetIaArr ();
      int *pja_hmatrT = hmatrT_str.GetJaArr ();

// Mask initial data via transposed one

      vector < int >imask_flt (nzja_hmatr + 1);
      int *pimask_flt = &imask_flt[0];

      int i;

      for (i = 0; i < nzja_hmatr; i++)
         pimask_flt[i] = 1;

      int *psubtree_beg = _tree.GetSubtreeBeg ();

      int nzja_flt = 0;

      int ibeg, j, jj, ind;

      for (i = 0; i < _nblks; i++) {
         ibeg = psubtree_beg[i];
         for (j = pia_hmatrT[i]; j < pia_hmatrT[i + 1]; j++) {
            jj = pja_hmatrT[j];
            if (jj < ibeg) {
               ind = pindHt2H[j];
               pimask_flt[ind] = -1;
               nzja_flt++;
            }
         }
      }

// Move filtered block into separate data

      vector < int >ia_flt (_nblks + 1);
      vector < int >ja_flt (nzja_flt + 1);
      vector < CMatrix < _Int, _Flt > >blksflt_arr (nzja_flt + 1);

      int *pia_flt = &ia_flt[0];
      int *pja_flt = &ja_flt[0];
      CMatrix < _Int, _Flt > *pblksflt_arr = &blksflt_arr[0];

      pia_flt[0] = 0;
      nzja_flt = 0;

      for (i = 0; i < _nblks; i++) {
         for (j = pia_hmatr[i]; j < pia_hmatr[i + 1]; j++) {
            jj = pja_hmatr[j];
            if (pimask_flt[j] < 0) {
               pja_flt[nzja_flt] = jj;
               pblksflt_arr[nzja_flt].ReplaceFree (pA_sub_AU[j]);
               nzja_flt++;
            }
         }
         pia_flt[i + 1] = nzja_flt;
      }

// Reformat initial hblock data

      vector < int >ia_new (_nblks + 1);
      int *pia_new = &ia_new[0];

      int nzja_new = 0;
      pia_new[0] = 0;

      for (i = 0; i < _nblks; i++) {
         for (j = pia_hmatr[i]; j < pia_hmatr[i + 1]; j++) {
            jj = pja_hmatr[j];
            if (pimask_flt[j] > 0) {
               pja_hmatr[nzja_new] = jj;
               if (nzja_new != j) {
                  pA_sub_AU[nzja_new].ReplaceFree (pA_sub_AU[j]);
               }
               nzja_new++;
            }
         }
         pia_new[i + 1] = nzja_new;
      }

      for (i = 0; i <= _nblks; i++)
         pia_hmatr[i] = pia_new[i];

      phmatr_AU->SetNzja (nzja_new);
      phmatr_AU->SetNzja2 (nzja_new);

      _AU_matr.SetNzblk (nzja_new);

// Perform diagonal modifications computations

      int nimax = 0;

      int niloc;

      for (i = 0; i < _nblks; i++) {
         niloc = (int) (_blks[i + 1] - _blks[i]);
         if (niloc > nimax)
            nimax = niloc;
      }

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      vector < int >icycle_thr (n_thr + 1);
      vector < CVectorData < int > >imaskblk_thr (n_thr + 1);

      int *picycle_thr = &icycle_thr[0];
      CVectorData < int >*pimaskblk_thr = &imaskblk_thr[0];

      for (i = 0; i < n_thr; i++)
         picycle_thr[i] = -1;

      vector < CMatrix < _Int, _Flt > >diamodL_arr (nzja_flt + 1);
      vector < CMatrix < _Int, _Flt > >diamodU_arr (nzja_flt + 1);

      CMatrix < _Int, _Flt > *pdiamodL_arr = &diamodL_arr[0];
      CMatrix < _Int, _Flt > *pdiamodU_arr = &diamodU_arr[0];

      _Flt theta_flt = (_Flt) _theta;

#ifdef USE_THREADS
#pragma omp parallel for
#endif
      for (int ipar = 0; ipar < nzja_flt; ipar++) {

         int my_thr = 0;
#ifdef USE_THREADS
         my_thr = omp_get_thread_num ();
#endif

         int ki, kj, kk, j, jj;
         _Flt aux, aux1, auxL, auxU;

         if (picycle_thr[my_thr] == -1) {
            pimaskblk_thr[my_thr].resize (nimax * 5 + 1);
            int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
            for (j = 0; j < nimax; j++)
               pimaskblk[j] = -1;
         }
         int icycleblk = picycle_thr[my_thr];
         int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
         int *plistblk = pimaskblk + nimax;
         int *pindblk = plistblk + nimax;
         int nlist_temp = pblksflt_arr[ipar].GetNlist ();
         int nzja_temp = pblksflt_arr[ipar].GetNzja ();
         _Int *plist_temp = pblksflt_arr[ipar].GetListArr ();
         _Int *pia_temp = pblksflt_arr[ipar].GetIaArr ();
         _Int *pja_temp = pblksflt_arr[ipar].GetJaArr ();
         _Flt *pa_temp = pblksflt_arr[ipar].GetAArr ();
         pdiamodL_arr[ipar].ResizeAndSetAll (0, nlist_temp, 0, 0, nlist_temp);
         _Int *plistL_temp = pdiamodL_arr[ipar].GetList2Arr ();
         _Flt *pdiaL_temp = pdiamodL_arr[ipar].GetAArr ();
         for (j = 0; j < nlist_temp; j++)
            plistL_temp[j] = plist_temp[j];
         CVector < _Flt >::SetByZeroes (nlist_temp, pdiaL_temp);
         int nlist_diaU = 0;
         icycleblk++;
         for (j = 0; j < nzja_temp; j++) {
            jj = (int) pja_temp[j];
            if (pimaskblk[jj] != icycleblk) {
               plistblk[nlist_diaU] = jj;
               nlist_diaU++;
               pimaskblk[jj] = icycleblk;
            }
         }
         sort (plistblk, plistblk + nlist_diaU);
         for (j = 0; j < nlist_diaU; j++) {
            jj = plistblk[j];
            pindblk[jj] = j;
         }
         pdiamodU_arr[ipar].ResizeAndSetAll (0, nlist_diaU, 0, 0, nlist_diaU);
         _Int *plistU_temp = pdiamodU_arr[ipar].GetList2Arr ();
         _Flt *pdiaU_temp = pdiamodU_arr[ipar].GetAArr ();
         for (j = 0; j < nlist_diaU; j++)
            plistU_temp[j] = plistblk[j];
         CVector < _Flt >::SetByZeroes (nlist_diaU, pdiaU_temp);
         for (ki = 0; ki < nlist_temp; ki++) {
            for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
               kk = (int) pja_temp[kj];
               auxL = pa_temp[kj * 2];
               auxU = pa_temp[kj * 2 + 1];
               aux = auxL;
               if (aux < 0.0e0)
                  aux = -aux;
               aux1 = auxU;
               if (aux1 < 0.0e0)
                  aux1 = -aux1;
               if (aux1 > aux)
                  aux = aux1;
               aux *= theta_flt;
               pdiaL_temp[ki] += aux;
               ind = pindblk[kk];
               pdiaU_temp[ind] += aux;
            }
         }
         picycle_thr[my_thr] = icycleblk;
      }

// Compute transposed sparsity of filtered data

      vector < int >ia_fltT (_nblks + 1);
      vector < int >ja_fltT (nzja_flt + 1);
      vector < int >ind_fltT2flt (nzja_flt + 1);

      int *pia_fltT = &ia_fltT[0];
      int *pja_fltT = &ja_fltT[0];
      int *pind_fltT2flt = &ind_fltT2flt[0];

      for (i = 0; i <= _nblks; i++)
         pia_fltT[i] = 0;

      for (i = 0; i < _nblks; i++) {
         for (j = pia_flt[i]; j < pia_flt[i + 1]; j++) {
            jj = pja_flt[j];
            pia_fltT[jj + 1]++;
         }
      }

      for (i = 0; i < _nblks; i++)
         pia_fltT[i + 1] = pia_fltT[i] + pia_fltT[i + 1];
      for (i = 0; i < _nblks; i++)
         pia_new[i] = pia_fltT[i];

      int k;

      for (i = 0; i < _nblks; i++) {
         for (j = pia_flt[i]; j < pia_flt[i + 1]; j++) {
            jj = pja_flt[j];
            k = pia_new[jj];
            pja_fltT[k] = i;
            pind_fltT2flt[k] = j;
            pia_new[jj]++;
         }
      }

// Perform diagonal modifications

      {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < _nblks; ipar++) {
            int kk, j;
            int indT, ind_dia, ind, k;
            ind = pia_hmatr[ipar];
            if (pia_hmatr[ipar] < pia_hmatr[ipar + 1]) {
               _Int *pia_temp = pA_sub_AU[ind].GetIaArr ();
               _Flt *pa_temp = pA_sub_AU[ind].GetAArr ();
               for (j = pia_flt[ipar]; j < pia_flt[ipar + 1]; j++) {
                  int nlistL_temp = pdiamodL_arr[j].GetNlist ();
                  _Int *plistL_temp = pdiamodL_arr[j].GetList2Arr ();
                  _Flt *pdiaL_temp = pdiamodL_arr[j].GetAArr ();
                  for (k = 0; k < nlistL_temp; k++) {
                     kk = (int) plistL_temp[k];
                     ind_dia = (int) pia_temp[kk];
                     pa_temp[ind_dia * 2] += pdiaL_temp[k];
                     pa_temp[ind_dia * 2 + 1] += pdiaL_temp[k];
                  }
               }
               for (j = pia_fltT[ipar]; j < pia_fltT[ipar + 1]; j++) {
                  indT = pind_fltT2flt[j];
                  int nlistU_temp = pdiamodU_arr[indT].GetNlist ();
                  _Int *plistU_temp = pdiamodU_arr[indT].GetList2Arr ();
                  _Flt *pdiaU_temp = pdiamodU_arr[indT].GetAArr ();
                  for (k = 0; k < nlistU_temp; k++) {
                     kk = (int) plistU_temp[k];
                     ind_dia = (int) pia_temp[kk];
                     pa_temp[ind_dia * 2] += pdiaU_temp[k];
                     pa_temp[ind_dia * 2 + 1] += pdiaU_temp[k];
                  }
               }
            }
         }
      }
   }

// Compute reordered by columns rectangular hmatrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFctThreads < _Int,
      _Flt >::ReorderHMatrixColumnsPairs (int _nblksR, long long *_blksR, int _nblksC_ini,
                                          long long *_blksC_ini, CBMatrix < _Int,
                                          _Flt > &_hmatr_ini, int *_order,
                                          int _nblksC_fin, long long *_blksC_fin,
                                          CBMatrix < _Int, _Flt > &_hmatr_fin)
   {

// Open hmatr

      CMatrix < int, float >*phmatr = _hmatr_ini.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub = _hmatr_ini.GetASubArr ();

      int *pia_hmatr = phmatr->GetIaArr ();
      int *pja_hmatr = phmatr->GetJaArr ();

// Create arrays that describe block numbers for each index

      int ntotC = (int) _blksC_ini[_nblksC_ini];

      CVectorData < int >ind2blkC_ini (ntotC + 1);
      CVectorData < int >ind2blkC_fin (ntotC + 1);

      int *pind2blkC_ini = ind2blkC_ini.Ptr ();
      int *pind2blkC_fin = ind2blkC_fin.Ptr ();

      int i;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < _nblksC_ini; ipar++) {
            int i;
            for (i = (int) _blksC_ini[ipar]; i < _blksC_ini[ipar + 1]; i++) {
               pind2blkC_ini[i] = ipar;
            }
         }
      }
      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < _nblksC_fin; ipar++) {
            int i;
            for (i = (int) _blksC_fin[ipar]; i < _blksC_fin[ipar + 1]; i++) {
               pind2blkC_fin[i] = ipar;
            }
         }
      }

// Create set of combined initial block rows

      vector < CMatrix < _Int, _Flt > >blkrows (_nblksR + 1);
      CMatrix < _Int, _Flt > *pblkrows = &blkrows[0];

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < _nblksR; ipar++) {

            int i, j, jj, k, kk, kkk, kj, jold, jnew, jblknew;

            int nlist_curr = (int) (_blksR[ipar + 1] - _blksR[ipar]);
            int nzja_curr = 0;
            for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
               jj = pja_hmatr[j];
               nzja_curr += pA_sub[j].GetNzja ();
            }
            CMatrix < _Int, _Flt > ASub;
            ASub.ResizeAndSetAll (nlist_curr, 0, nzja_curr, nzja_curr, 2 * nzja_curr);
            _Int *pia_curr = ASub.GetIaArr ();
            _Int *pja_curr = ASub.GetJaArr ();
            _Int *pja2_curr = ASub.GetJa2Arr ();
            _Flt *pa_curr = ASub.GetAArr ();
            for (i = 0; i <= nlist_curr; i++)
               pia_curr[i] = 0;
            for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
               jj = pja_hmatr[j];
               int nlist_temp = pA_sub[j].GetNlist ();
               _Int *plist_temp = pA_sub[j].GetListArr ();
               _Int *pia_temp = pA_sub[j].GetIaArr ();
               for (k = 0; k < nlist_temp; k++) {
                  kk = (int) plist_temp[k];
                  pia_curr[kk + 1] += (pia_temp[k + 1] - pia_temp[k]);
               }
            }
            for (i = 0; i < nlist_curr; i++)
               pia_curr[i + 1] = pia_curr[i] + pia_curr[i + 1];
            vector < int >iptr (nlist_curr + 1);
            int *piptr = &iptr[0];
            for (i = 0; i < nlist_curr; i++)
               piptr[i] = (int) pia_curr[i];
            for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
               jj = pja_hmatr[j];
               int nlist_temp = pA_sub[j].GetNlist ();
               _Int *plist_temp = pA_sub[j].GetListArr ();
               _Int *pia_temp = pA_sub[j].GetIaArr ();
               _Int *pja_temp = pA_sub[j].GetJaArr ();
               _Flt *pa_temp = pA_sub[j].GetAArr ();
               for (k = 0; k < nlist_temp; k++) {
                  kk = (int) plist_temp[k];
                  kkk = piptr[kk];
                  for (kj = (int) pia_temp[k]; kj < pia_temp[k + 1]; kj++) {
                     jold = (int) (pja_temp[kj] + _blksC_ini[jj]);
                     jnew = _order[jold];
                     jblknew = pind2blkC_fin[jnew];
                     pja_curr[kkk] = (_Int) (jnew - _blksC_fin[jblknew]);
                     pja2_curr[kkk] = (_Int) jblknew;
                     pa_curr[kkk * 2] = pa_temp[kj * 2];
                     pa_curr[kkk * 2 + 1] = pa_temp[kj * 2 + 1];
                     kkk++;
                  }
                  piptr[kk] = kkk;
               }
            }
            pblkrows[ipar].ReplaceFree (ASub);
         }
      }

// Compute maximal size

      int nimax = 0;

      int niloc;

      for (i = 0; i < _nblksC_fin; i++) {
         niloc = (int) (_blksC_fin[i + 1] - _blksC_fin[i]);
         if (niloc > nimax)
            nimax = niloc;
      }

// Create set of new block rows

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      vector < int >icycle_thr (n_thr + 1);
      vector < CVectorData < int > >imaskblk_thr (n_thr + 1);

      int *picycle_thr = &icycle_thr[0];
      CVectorData < int >*pimaskblk_thr = &imaskblk_thr[0];

      for (i = 0; i < n_thr; i++)
         picycle_thr[i] = -1;

      vector < CBMatrix < _Int, _Flt > >hblkrows_new (_nblksR + 1);
      CBMatrix < _Int, _Flt > *phblkrows_new = &hblkrows_new[0];

      {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < _nblksR; ipar++) {

            int nlist_curr = (int) (_blksR[ipar + 1] - _blksR[ipar]);

            int my_thr = 0;
#ifdef USE_THREADS
            my_thr = omp_get_thread_num ();
#endif
            int i, j, jj, k, niloc;

            if (picycle_thr[my_thr] == -1) {
               pimaskblk_thr[my_thr].resize (3 * _nblksC_fin + 1);
               int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
               for (j = 0; j < _nblksC_fin; j++)
                  pimaskblk[j] = -1;
            }
            int icycleblk = picycle_thr[my_thr];
            int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
            int *plistblk = pimaskblk + _nblksC_fin;
            int *pindblk = plistblk + _nblksC_fin;

// Get block row

            CMatrix < _Int, _Flt > *pASub = pblkrows + ipar;

            int nzja_curr = pASub->GetNzja ();
            _Int *pia_curr = pASub->GetIaArr ();
            _Int *pja_curr = pASub->GetJaArr ();
            _Int *pja2_curr = pASub->GetJa2Arr ();
            _Flt *pa_curr = pASub->GetAArr ();

// Split into the set of blocks

            icycleblk++;

            int nlistblk = 0;

            int jblk;

            for (i = 0; i < nzja_curr; i++) {
               jblk = (int) pja2_curr[i];
               if (pimaskblk[jblk] != icycleblk) {
                  plistblk[nlistblk] = jblk;
                  nlistblk++;
                  pimaskblk[jblk] = icycleblk;
               }
            }

            sort (plistblk, plistblk + nlistblk);

            for (i = 0; i < nlistblk; i++) {
               jblk = plistblk[i];
               pindblk[jblk] = i;
            }

            vector < int >nzblk_arr (nlistblk + 1);
            int *pnzblk_arr = &nzblk_arr[0];

            for (i = 0; i < nlistblk; i++)
               pnzblk_arr[i] = 0;

            int ind;

            for (i = 0; i < nzja_curr; i++) {
               jblk = (int) pja2_curr[i];
               ind = pindblk[jblk];
               pnzblk_arr[ind]++;
            }

            vector < vector < _Int > >rows_arr (nlistblk + 1);
            vector < vector < _Int > >cols_arr (nlistblk + 1);
            vector < vector < _Flt > >elems_arr (nlistblk + 1);

            vector < _Int > *prows_arr = &rows_arr[0];
            vector < _Int > *pcols_arr = &cols_arr[0];
            vector < _Flt > *pelems_arr = &elems_arr[0];

            for (i = 0; i < nlistblk; i++) {
               prows_arr[i].resize (pnzblk_arr[i] + 1);
               pcols_arr[i].resize (pnzblk_arr[i] + 1);
               pelems_arr[i].resize (2 * pnzblk_arr[i] + 1);
            }

            for (i = 0; i < nlistblk; i++)
               pnzblk_arr[i] = 0;

            int jj2;

            for (i = 0; i < nlist_curr; i++) {
               for (j = (int) pia_curr[i]; j < pia_curr[i + 1]; j++) {
                  jj = (int) pja_curr[j];
                  jj2 = (int) pja2_curr[j];
                  ind = pindblk[jj2];
                  k = pnzblk_arr[ind];
                  _Int *pprows_arr = &prows_arr[ind][0];
                  _Int *ppcols_arr = &pcols_arr[ind][0];
                  _Flt *ppelems_arr = &pelems_arr[ind][0];
                  pprows_arr[k] = i;
                  ppcols_arr[k] = jj;
                  ppelems_arr[k * 2] = pa_curr[j * 2];
                  ppelems_arr[k * 2 + 1] = pa_curr[j * 2 + 1];
                  pnzblk_arr[ind]++;
               }
            }

// Compute hblock

            phblkrows_new[ipar].SetNzblk (nlistblk);
            phblkrows_new[ipar].ResizeASub (nlistblk);

            CMatrix < int, float >*phmatr_new = phblkrows_new[ipar].GetHMatrStr ();
            CMatrix < _Int, _Flt > *pA_sub_new = phblkrows_new[ipar].GetASubArr ();

            phmatr_new->SetNzja (nlistblk);
            phmatr_new->ResizeJa (nlistblk);

            int *pja_hmatr_new = phmatr_new->GetJaArr ();

            for (i = 0; i < nlistblk; i++)
               pja_hmatr_new[i] = plistblk[i];

            CVectorData < CSortInt > iiarr (nimax + 1);
            CVectorData < _Flt > elemsarr (nimax * 2 + 1);

            vector < int >imaskarr (nlist_curr + 1);
            vector < int >listarr (nlist_curr + 1);
            vector < int >indarr (nlist_curr + 1);

            CSortInt *piiarr = iiarr.Ptr ();
            _Flt *pelemsarr = elemsarr.Ptr ();

            int *pimaskarr = &imaskarr[0];
            int *plistarr = &listarr[0];
            int *pindarr = &indarr[0];

            for (i = 0; i < nlist_curr; i++)
               pimaskarr[i] = -1;

            int icycle = -1;
            int ibeg;

            for (i = 0; i < nlistblk; i++) {
               icycle++;
               int nzja_temp = pnzblk_arr[i];
               _Int *pprows_arr = &prows_arr[i][0];
               _Int *ppcols_arr = &pcols_arr[i][0];
               _Flt *ppelems_arr = &pelems_arr[i][0];
               int nlist_temp = 0;
               for (j = 0; j < nzja_temp; j++) {
                  jj = (int) pprows_arr[j];
                  if (pimaskarr[jj] != icycle) {
                     plistarr[nlist_temp] = jj;
                     nlist_temp++;
                     pimaskarr[jj] = icycle;
                  }
               }
               sort (plistarr, plistarr + nlist_temp);
               for (j = 0; j < nlist_temp; j++) {
                  jj = plistarr[j];
                  pindarr[jj] = j;
               }
               pA_sub_new[i].ResizeAndSetAll (nlist_temp, 0, nzja_temp, 0, nzja_temp * 2);
               _Int *plist_temp = pA_sub_new[i].GetListArr ();
               _Int *pia_temp = pA_sub_new[i].GetIaArr ();
               _Int *pja_temp = pA_sub_new[i].GetJaArr ();
               _Flt *pa_temp = pA_sub_new[i].GetAArr ();
               for (j = 0; j < nlist_temp; j++)
                  plist_temp[j] = plistarr[j];
               for (j = 0; j <= nlist_temp; j++)
                  pia_temp[j] = 0;
               for (j = 0; j < nzja_temp; j++) {
                  jj = (int) pprows_arr[j];
                  ind = pindarr[jj];
                  pia_temp[ind + 1]++;
               }
               for (j = 0; j < nlist_temp; j++)
                  pia_temp[j + 1] = pia_temp[j] + pia_temp[j + 1];
               for (j = 0; j < nlist_temp; j++)
                  plistarr[j] = (int) pia_temp[j];
               for (j = 0; j < nzja_temp; j++) {
                  jj = (int) pprows_arr[j];
                  ind = pindarr[jj];
                  k = plistarr[ind];
                  pja_temp[k] = ppcols_arr[j];
                  pa_temp[k * 2] = ppelems_arr[j * 2];
                  pa_temp[k * 2 + 1] = ppelems_arr[j * 2 + 1];
                  plistarr[ind]++;
               }
               for (j = 0; j < nlist_temp; j++) {
                  ibeg = (int) pia_temp[j];
                  niloc = (int) (pia_temp[j + 1] - pia_temp[j]);
                  for (k = (int) pia_temp[j]; k < pia_temp[j + 1]; k++) {
                     piiarr[k - ibeg].ival = (int) pja_temp[k];
                     piiarr[k - ibeg].i2val = k;
                     pelemsarr[(k - ibeg) * 2] = pa_temp[k * 2];
                     pelemsarr[(k - ibeg) * 2 + 1] = pa_temp[k * 2 + 1];
                  }
                  sort (piiarr, piiarr + niloc);
                  for (k = (int) pia_temp[j]; k < pia_temp[j + 1]; k++) {
                     pja_temp[k] = (_Int) piiarr[k - ibeg].ival;
                     ind = piiarr[k - ibeg].i2val;
                     pa_temp[k * 2] = pelemsarr[(ind - ibeg) * 2];
                     pa_temp[k * 2 + 1] = pelemsarr[(ind - ibeg) * 2 + 1];
                  }
               }
            }

            picycle_thr[my_thr] = icycleblk;

         }

      }

// Combine computed data into one hmatrix

      int nzja_tot = 0;

      for (i = 0; i < _nblksR; i++) {
         nzja_tot += phblkrows_new[i].GetNzblk ();
      }

      _hmatr_fin.SetNzblk (nzja_tot);
      _hmatr_fin.ResizeASub (nzja_tot);

      CMatrix < int, float >*phmatr_fin = _hmatr_fin.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_fin = _hmatr_fin.GetASubArr ();

      phmatr_fin->ResizeAndSetAllSp (_nblksR, _nblksR, nzja_tot, nzja_tot);

      int *plist_hmatr_fin = phmatr_fin->GetListArr ();
      int *plist2_hmatr_fin = phmatr_fin->GetList2Arr ();
      int *pia_hmatr_fin = phmatr_fin->GetIaArr ();
      int *pja_hmatr_fin = phmatr_fin->GetJaArr ();
      int *pja2_hmatr_fin = phmatr_fin->GetJa2Arr ();

      for (i = 0; i < _nblksR; i++)
         plist_hmatr_fin[i] = i;
      for (i = 0; i < _nblksR; i++)
         plist2_hmatr_fin[i] = 0;
      for (i = 0; i < nzja_tot; i++)
         pja2_hmatr_fin[i] = 0;

      nzja_tot = 0;
      pia_hmatr_fin[0] = 0;

      int j;

      for (i = 0; i < _nblksR; i++) {
         int nzblk_temp = phblkrows_new[i].GetNzblk ();
         CMatrix < int, float >*phmatr_temp = phblkrows_new[i].GetHMatrStr ();
         CMatrix < _Int, _Flt > *pA_sub_temp = phblkrows_new[i].GetASubArr ();
         int *pja_temp = phmatr_temp->GetJaArr ();
         for (j = 0; j < nzblk_temp; j++) {
            pja_hmatr_fin[nzja_tot] = pja_temp[j];
            pA_sub_fin[nzja_tot].ReplaceFree (pA_sub_temp[j]);
            nzja_tot++;
         }
         pia_hmatr_fin[i + 1] = nzja_tot;
      }

   }

// Init MvmA data
//========================================================================================
   template < typename _Int, typename _Flt > void CMvmParThreads_base < _Int,
      _Flt >::InitMvmA_base (CBMatrix < _Int, _Flt > *_hmatr_arr)
   {

// Get control data

      void *pcomm_loc = this->pcomm;

      int *phblk2cpu_loc = this->phblk2cpu;
      int *phblk2blks_loc = this->phblk2blks;
      int *pblk2hblks_loc = this->pblk2hblks;
      int nblks_loc = this->nblks;
      long long *pblks_loc = this->pblks;
      int *pibsblk_loc = &this->ibsblk[0];
      int nlisthblk_own_loc = this->nlisthblk_own;
      int *plisthblk_own_loc = &this->listhblk_own[0];

      int nproc = CMPIDataExchange::GetNproc (pcomm_loc);
      int myid = CMPIDataExchange::GetMyid (pcomm_loc);

// Compute maximal block size

      int nimax = 0;

      int i, niloc;

      for (i = 0; i < nblks_loc; i++) {
         niloc = (int) (pblks_loc[i + 1] - pblks_loc[i]);
         if (niloc > nimax)
            nimax = niloc;
      }

// Store pointer to matrix data

      this->phmatr = _hmatr_arr;

// Compute list of external blocks

      vector < int >imaskblk (nblks_loc + 1);
      vector < int >listblk (nblks_loc + 1);
      vector < int >indblk (nblks_loc + 1);

      int *pimaskblk = &imaskblk[0];
      int *plistblk = &listblk[0];
      int *pindblk = &indblk[0];

      for (i = 0; i < nblks_loc; i++)
         pimaskblk[i] = -1;

      int icycleblk = -1;

      int j, ihblk, jblk, jhblk, jblkgl;

      icycleblk++;

      int nlistblk_ext = 0;

      for (i = 0; i < nlisthblk_own_loc; i++) {
         ihblk = plisthblk_own_loc[i];
         CMatrix < int, float >*phmatr_str_loc = _hmatr_arr[ihblk].GetHMatrStr ();
         int nzja_hblk = phmatr_str_loc->GetNzja ();
         int *pja_hblk = phmatr_str_loc->GetJaArr ();
         int *pja2_hblk = phmatr_str_loc->GetJa2Arr ();
         for (j = 0; j < nzja_hblk; j++) {
            jblk = pja_hblk[j];
            jhblk = pja2_hblk[j];
            jblkgl = phblk2blks_loc[jhblk] + jblk;
            if (phblk2cpu_loc[jhblk] != myid && pimaskblk[jblkgl] != icycleblk) {
               plistblk[nlistblk_ext] = jblkgl;
               nlistblk_ext++;
               pimaskblk[jblkgl] = icycleblk;
            }
         }
      }

      sort (plistblk, plistblk + nlistblk_ext);

      int iblk;

      for (i = 0; i < nlistblk_ext; i++) {
         iblk = plistblk[i];
         pindblk[iblk] = i;
      }

      this->nblks_recvs = nlistblk_ext;
      this->listblk_recvs.resize (nlistblk_ext + 1);

      int *plistblk_recvs = &this->listblk_recvs[0];

      for (i = 0; i < nlistblk_ext; i++)
         plistblk_recvs[i] = plistblk[i];

// Compute the set of indices to be received and multiplications lists

      vector < int >ia_list_rcv (nlistblk_ext + 1);
      vector < vector < int > >ja_triple_rcv (nlistblk_ext + 1);

      int *pia_list_rcv = &ia_list_rcv[0];
      vector < int >*pja_triple_rcv = &ja_triple_rcv[0];

      for (i = 0; i <= nlistblk_ext; i++)
         pia_list_rcv[i] = 0;

      int ind, ilist;

      for (ilist = 0; ilist < nlisthblk_own_loc; ilist++) {
         ihblk = plisthblk_own_loc[ilist];
         CMatrix < int, float >*phmatr_str_loc = _hmatr_arr[ihblk].GetHMatrStr ();
         int nlist_hblk = phmatr_str_loc->GetNlist ();
         int *pia_hblk = phmatr_str_loc->GetIaArr ();
         int *pja_hblk = phmatr_str_loc->GetJaArr ();
         int *pja2_hblk = phmatr_str_loc->GetJa2Arr ();
         for (i = 0; i < nlist_hblk; i++) {
            for (j = pia_hblk[i]; j < pia_hblk[i + 1]; j++) {
               jblk = pja_hblk[j];
               jhblk = pja2_hblk[j];
               jblkgl = phblk2blks_loc[jhblk] + jblk;
               if (phblk2cpu_loc[jhblk] != myid) {
                  ind = pindblk[jblkgl];
                  pja_triple_rcv[ind].push_back (i);
                  pja_triple_rcv[ind].push_back (ihblk);
                  pja_triple_rcv[ind].push_back (j);
                  pia_list_rcv[ind + 1]++;
               }
            }
         }
      }

      for (i = 0; i < nlistblk_ext; i++)
         pia_list_rcv[i + 1] = pia_list_rcv[i] + pia_list_rcv[i + 1];

      this->ialist_recvs.resize (nlistblk_ext + 1);

      int *pialist_recvs = &this->ialist_recvs[0];

      for (i = 0; i <= nlistblk_ext; i++)
         pialist_recvs[i] = pia_list_rcv[i];

      int nzja_pairs = pialist_recvs[nlistblk_ext];

      this->jatriples_recvs.resize (3 * nzja_pairs + 1);

      int *pjatriples_recvs = &this->jatriples_recvs[0];

      int ibeg, jloc;

      for (i = 0; i < nlistblk_ext; i++) {
         int *ppja_triple_rcv = &(pja_triple_rcv[i][0]);
         ibeg = pia_list_rcv[i];
         for (j = pia_list_rcv[i]; j < pia_list_rcv[i + 1]; j++) {
            jloc = j - ibeg;
            pjatriples_recvs[j * 3] = ppja_triple_rcv[jloc * 3];
            pjatriples_recvs[j * 3 + 1] = ppja_triple_rcv[jloc * 3 + 1];
            pjatriples_recvs[j * 3 + 2] = ppja_triple_rcv[jloc * 3 + 2];
         }
      }

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      vector < int >icycle_thr (n_thr + 1);
      vector < CVectorData < int > >imask_thr (n_thr + 1);

      int *picycle_thr = &icycle_thr[0];
      CVectorData < int >*pimask_thr = &imask_thr[0];

      for (i = 0; i < n_thr; i++)
         picycle_thr[i] = -1;

      vector < int >ia_blk_rcv (nlistblk_ext + 1);
      vector < vector < int > >ja_blk_rcv (nlistblk_ext + 1);

      int *pia_blk_rcv = &ia_blk_rcv[0];
      vector < int >*pja_blk_rcv = &ja_blk_rcv[0];

      for (i = 0; i <= nlistblk_ext; i++)
         pia_blk_rcv[i] = 0;

      {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nlistblk_ext; ipar++) {

            int my_thr = 0;
#ifdef USE_THREADS
            my_thr = omp_get_thread_num ();
#endif

            int j, ihblk, ind, nlistloc, k, jj;

            if (picycle_thr[my_thr] == -1) {
               pimask_thr[my_thr].resize (2 * nimax + 1);
               int *pimask = pimask_thr[my_thr].Ptr ();
               for (j = 0; j < nimax; j++)
                  pimask[j] = -1;
            }
            int icycle_work = picycle_thr[my_thr];
            int *pimask_work = pimask_thr[my_thr].Ptr ();
            int *plist_work = pimask_work + nimax;
            icycle_work++;
            nlistloc = 0;
            for (j = pialist_recvs[ipar]; j < pialist_recvs[ipar + 1]; j++) {
               ihblk = pjatriples_recvs[j * 3 + 1];
               ind = pjatriples_recvs[j * 3 + 2];
               CMatrix < _Int, _Flt > *pasub_loc = _hmatr_arr[ihblk].GetASubArr ();
               int nzjaloc = pasub_loc[ind].GetNzja ();
               _Int *pjaloc = pasub_loc[ind].GetJaArr ();
               for (k = 0; k < nzjaloc; k++) {
                  jj = (int) pjaloc[k];
                  if (pimask_work[jj] != icycle_work) {
                     plist_work[nlistloc] = jj;
                     nlistloc++;
                     pimask_work[jj] = icycle_work;
                  }
               }
            }
            sort (plist_work, plist_work + nlistloc);
            pia_blk_rcv[ipar + 1] = nlistloc;
            pja_blk_rcv[ipar].resize (nlistloc + 1);
            int *ppja_blk_rcv = &(pja_blk_rcv[ipar][0]);
            for (j = 0; j < nlistloc; j++)
               ppja_blk_rcv[j] = plist_work[j];
            picycle_thr[my_thr] = icycle_work;
         }
      }

      for (i = 0; i < nlistblk_ext; i++)
         pia_blk_rcv[i + 1] = pia_blk_rcv[i] + pia_blk_rcv[i + 1];

      this->iablk_recvs.resize (nlistblk_ext + 1);

      int *piablk_recvs = &this->iablk_recvs[0];

      for (i = 0; i <= nlistblk_ext; i++)
         piablk_recvs[i] = pia_blk_rcv[i];

      int nzja_ind = pia_blk_rcv[nlistblk_ext];

      this->ind_recvs.resize (nzja_ind + 1);

      int *pind_recvs = &this->ind_recvs[0];

      for (i = 0; i < nlistblk_ext; i++) {
         int *ppja_pair_rcv = &(pja_blk_rcv[i][0]);
         ibeg = pia_blk_rcv[i];
         for (j = pia_blk_rcv[i]; j < pia_blk_rcv[i + 1]; j++) {
            jloc = j - ibeg;
            pind_recvs[j] = ppja_pair_rcv[jloc];
         }
      }

//      this->x_recv.resize (nzja_ind + 1);

      this->i_size_x_recv = nzja_ind + 1;

// Finally create cpu recv data

      vector < int >imaskcpu (nproc);
      vector < int >listcpu (nproc);
      vector < int >indcpu (nproc);

      int *pimaskcpu = &imaskcpu[0];
      int *plistcpu = &listcpu[0];
      int *pindcpu = &indcpu[0];

      for (i = 0; i < nproc; i++)
         pimaskcpu[i] = -1;

      int icyclecpu = -1;

      icyclecpu++;

      int nlistcpu = 0;

      int jcpu;

      for (i = 0; i < nlistblk_ext; i++) {
         jblkgl = plistblk_recvs[i];
         jhblk = pblk2hblks_loc[jblkgl];
         jcpu = phblk2cpu_loc[jhblk];
         if (pimaskcpu[jcpu] != icyclecpu) {
            plistcpu[nlistcpu] = jcpu;
            nlistcpu++;
            pimaskcpu[jcpu] = icyclecpu;
         }
      }

      sort (plistcpu, plistcpu + nlistcpu);

      for (i = 0; i < nlistcpu; i++) {
         jcpu = plistcpu[i];
         pindcpu[jcpu] = i;
      }

      this->nrecvs = nlistcpu;

      this->rcv2cpu.resize (nlistcpu + 1);
      this->ia_recvs.resize (nlistcpu + 1);
      this->ja_recvs.resize (2 * nzja_ind + 1);

      int *prcv2cpu = &this->rcv2cpu[0];
      int *pia_recvs = &this->ia_recvs[0];
      int *pja_recvs = &this->ja_recvs[0];

      for (i = 0; i < nlistcpu; i++)
         prcv2cpu[i] = plistcpu[i];
      for (i = 0; i <= nlistcpu; i++)
         pia_recvs[i] = 0;

      for (i = 0; i < nlistblk_ext; i++) {
         jblkgl = plistblk_recvs[i];
         jhblk = pblk2hblks_loc[jblkgl];
         jcpu = phblk2cpu_loc[jhblk];
         ind = pindcpu[jcpu];
         pia_recvs[ind + 1] += (pia_blk_rcv[i + 1] - pia_blk_rcv[i]);
      }

      for (i = 0; i < nlistcpu; i++)
         pia_recvs[i + 1] = pia_recvs[i] + pia_recvs[i + 1];

      for (i = 0; i < nlistcpu; i++)
         plistcpu[i] = pia_recvs[i];

      int k, jj;

      for (i = 0; i < nlistblk_ext; i++) {
         jblkgl = plistblk_recvs[i];
         jhblk = pblk2hblks_loc[jblkgl];
         jcpu = phblk2cpu_loc[jhblk];
         ind = pindcpu[jcpu];
         k = plistcpu[ind];
         for (j = pia_blk_rcv[i]; j < pia_blk_rcv[i + 1]; j++) {
            jj = pind_recvs[j];
            pja_recvs[k * 2] = jj;
            pja_recvs[k * 2 + 1] = jblkgl;
            pind_recvs[j] = k;
            k++;
         }
         plistcpu[ind] = k;
      }

// Prepare exchange data

      vector < CBMatrix < int, _Flt > >hblk_send (nlistcpu + 1);

      CBMatrix < int, _Flt > *phblk_send = &hblk_send[0];

      for (i = 0; i < nlistcpu; i++) {
         phblk_send[i].SetNzblk (1);
         phblk_send[i].ResizeASub (1);
         CMatrix < int, _Flt > *pA_sub = phblk_send[i].GetASubArr ();
         CMatrix < int, _Flt > ablk_temp;
         int niloc = pia_recvs[i + 1] - pia_recvs[i];
         ablk_temp.ResizeAndSetAllSp (0, 0, niloc * 2, 0);
         int *pjaloc = ablk_temp.GetJaArr ();
         ibeg = pia_recvs[i];
         for (j = pia_recvs[i]; j < pia_recvs[i + 1]; j++) {
            jloc = j - ibeg;
            pjaloc[jloc * 2] = pja_recvs[j * 2];
            pjaloc[jloc * 2 + 1] = pja_recvs[j * 2 + 1];
         }
         pA_sub->ReplaceFree (ablk_temp);
      }

// Pack send data

      vector < int >CpuIDSend (nlistcpu);
      vector < vector < char > >ObjSend (nlistcpu);

      int *pCpuIDSend = NULL;
      vector < char >*pObjSend = NULL;

      if (nlistcpu > 0) {

         pCpuIDSend = &CpuIDSend[0];
         pObjSend = &ObjSend[0];

      }

      long long isize;
      char *pobj;

      for (i = 0; i < nlistcpu; i++) {
         pCpuIDSend[i] = prcv2cpu[i];
         isize = phblk_send[i].GetPackedSize ();
         pObjSend[i].resize ((size_t) isize);
         pobj = &(pObjSend[i][0]);
         phblk_send[i].FillPacked (isize, pobj);
         phblk_send[i].Clean ();
      }

// Exchange

      vector < int >CpuIDRecv;
      vector < vector < char > >ObjRecv;

      CMPIDataExchange::DataExchange (pcomm_loc, CpuIDSend, ObjSend, CpuIDRecv, ObjRecv);

// Free send data

      {
         vector < int >CpuIDSend_temp;
         vector < vector < char > >ObjSend_temp;
         CpuIDSend.swap (CpuIDSend_temp);
         ObjSend.swap (ObjSend_temp);
      }

// Unpack receive data

      int nrecv_loc = (int) CpuIDRecv.size ();

      vector < char >*pObjRecv = NULL;
      int *pCpuIDRecv = NULL;

      if (nrecv_loc > 0) {
         pObjRecv = &ObjRecv[0];
         pCpuIDRecv = &CpuIDRecv[0];
      }

      vector < CBMatrix < int, _Flt > >hblk_recv (nrecv_loc + 1);

      CBMatrix < int, _Flt > *phblk_recv = &hblk_recv[0];

      for (i = 0; i < nrecv_loc; i++) {
         isize = (long long) pObjRecv[i].size ();
         pobj = &(pObjRecv[i][0]);
         phblk_recv[i].UnPack (isize, pobj);
      }

// Free recv data

      {
         vector < vector < char > >ObjRecv_temp;
         ObjRecv.swap (ObjRecv_temp);
      }

// Compute correct ordering of cpu data

      vector < CSortInt > iiarr (nrecv_loc + 1);
      CSortInt *piiarr = &iiarr[0];

      for (i = 0; i < nrecv_loc; i++) {
         piiarr[i].ival = pCpuIDRecv[i];
         piiarr[i].i2val = i;
      }

      sort (piiarr, piiarr + nrecv_loc);

// Store received data

      this->nsends = nrecv_loc;

      this->snd2cpu.resize (nrecv_loc + 1);
      this->ia_sends.resize (nrecv_loc + 1);

      int *psnd2cpu = &this->snd2cpu[0];
      int *pia_sends = &this->ia_sends[0];

      int nz_sends = 0;

      pia_sends[0] = 0;

      for (i = 0; i < nrecv_loc; i++) {
         ind = piiarr[i].i2val;
         psnd2cpu[i] = piiarr[i].ival;
         CMatrix < int, _Flt > *pA_sub = phblk_recv[ind].GetASubArr ();
         int nzjaloc = pA_sub->GetNzja () / 2;
         nz_sends += nzjaloc;
         pia_sends[i + 1] = nz_sends;
      }

      this->ind_sends.resize (nz_sends + 1);

      this->i_size_x_send = nz_sends + 1;

//      this->x_send.resize (nz_sends + 1);

      int *pind_sends = &this->ind_sends[0];

      nz_sends = 0;

      int jj2, ibs;

      for (i = 0; i < nrecv_loc; i++) {
         ind = piiarr[i].i2val;
         CMatrix < int, _Flt > *pA_sub = phblk_recv[ind].GetASubArr ();
         int nzjaloc = pA_sub->GetNzja () / 2;
         int *pjaloc = pA_sub->GetJaArr ();
         for (j = 0; j < nzjaloc; j++) {
            jj = pjaloc[j * 2];
            jj2 = pjaloc[j * 2 + 1];
            ibs = pibsblk_loc[jj2];
            pind_sends[nz_sends] = ibs + jj;
            nz_sends++;
         }
      }

      this->i_size_x_temp = nimax + 1;

//      this->x_temp.resize (nimax + 1);

   }

// Init MvmA data
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CMvmParThreads < _Int, _Flt,
      _FltVect >::InitMvmA (CBMatrix < _Int, _Flt > *_hmatr_arr)
   {

      this->InitMvmA_base (_hmatr_arr);

      this->x_recv.resize (this->i_size_x_recv);
      this->x_send.resize (this->i_size_x_send);
      this->x_temp.resize (this->i_size_x_temp);

   }

// Perform multiplication by A
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CMvmParThreads < _Int, _Flt,
      _FltVect >::MvmA (const _FltVect * _x, _FltVect * _ax)
   {

// Open mvm structure

      void *pcomm_loc = this->pcomm;

      int *phblk2cpu_loc = this->phblk2cpu;
      int *phblk2blks_loc = this->phblk2blks;
      int ni_cpu_loc = this->ni_cpu;
      int *pibsblk_loc = &this->ibsblk[0];
      int nlisthblk_own_loc = this->nlisthblk_own;
      int *plisthblk_own_loc = &this->listhblk_own[0];
      CBMatrix < _Int, _Flt > *phmatr_loc = this->phmatr;
      int nsends_loc = this->nsends;
      int *psnd2cpu_loc = &this->snd2cpu[0];
      int *pia_sends_loc = &this->ia_sends[0];
      int *pind_sends_loc = &this->ind_sends[0];
      _FltVect *px_send_loc = this->x_send.Ptr ();
      int nrecvs_loc = this->nrecvs;
      int *prcv2cpu_loc = &this->rcv2cpu[0];
      int *pia_recvs_loc = &this->ia_recvs[0];
      int *pja_recvs_loc = &this->ja_recvs[0];
      _FltVect *px_recv_loc = this->x_recv.Ptr ();
      int nblks_recvs_loc = this->nblks_recvs;
      int *plistblk_recvs_loc = &this->listblk_recvs[0];
      int *pialist_recvs_loc = &this->ialist_recvs[0];
      int *pjatriples_recvs_loc = &this->jatriples_recvs[0];
      int *piablk_recvs_loc = &this->iablk_recvs[0];
      int *pind_recvs_loc = &this->ind_recvs[0];
      _FltVect *px_temp_loc = this->x_temp.Ptr ();

      int myid = CMPIDataExchange::GetMyid (pcomm_loc);

// Init array ax by zeroes

      CVector < _FltVect >::SetByZeroes_thr (ni_cpu_loc, _ax);

// Prepare send

      int ni_send_loc = pia_sends_loc[nsends_loc];

      int i, ind;

      for (i = 0; i < ni_send_loc; i++) {
         ind = pind_sends_loc[i];
         px_send_loc[i] = _x[ind];
      }

      int ni_recv_loc = pia_recvs_loc[nrecvs_loc];

      CVector < _FltVect >::SetByZeroes_thr (ni_recv_loc, px_recv_loc);

// Init async recvs and sends

      void *psndrcv_recvs_loc;
      void *psndrcv_stats_loc;

      CMPIDataExchange::AllocateRecvs (nrecvs_loc + nsends_loc, psndrcv_recvs_loc);
      CMPIDataExchange::AllocateStats (nrecvs_loc + nsends_loc, psndrcv_stats_loc);

      int icpu, isize, ibs;

      for (i = 0; i < nrecvs_loc; i++) {
         icpu = prcv2cpu_loc[i];
         isize = (pia_recvs_loc[i + 1] - pia_recvs_loc[i]) * sizeof (_FltVect);
         ibs = pia_recvs_loc[i];
         CMPIDataExchange::IRecv (pcomm_loc, icpu, icpu, isize,
                                  (char *) (px_recv_loc + ibs), i, psndrcv_recvs_loc);
      }

      for (i = 0; i < nsends_loc; i++) {
         icpu = psnd2cpu_loc[i];
         isize = (pia_sends_loc[i + 1] - pia_sends_loc[i]) * sizeof (_FltVect);
         ibs = pia_sends_loc[i];
         CMPIDataExchange::ISend (pcomm_loc, icpu, myid, isize,
                                  (char *) (px_send_loc + ibs), i + nrecvs_loc,
                                  psndrcv_recvs_loc);
      }

// Perform local multiplicaions

      int ihblk, ibs_i, j, jblkgl, ilist, iblkgl;

      for (ilist = 0; ilist < nlisthblk_own_loc; ilist++) {
         ihblk = plisthblk_own_loc[ilist];
         CMatrix < _Int, _Flt > *pasub_loc = phmatr_loc[ihblk].GetASubArr ();
         CMatrix < int, float >*phmatr_str_loc = phmatr_loc[ihblk].GetHMatrStr ();
         int nlist_hblk = phmatr_str_loc->GetNlist ();
         int *pia_hblk = phmatr_str_loc->GetIaArr ();
         int *pja_hblk = phmatr_str_loc->GetJaArr ();
         int *pja2_hblk = phmatr_str_loc->GetJa2Arr ();
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nlist_hblk; ipar++) {
            int iblkgl = phblk2blks_loc[ihblk] + ipar;
            int ibs_i = pibsblk_loc[iblkgl];
            for (int j = pia_hblk[ipar]; j < pia_hblk[ipar + 1]; j++) {
               int jblk = pja_hblk[j];
               int jhblk = pja2_hblk[j];
               int jblkgl = phblk2blks_loc[jhblk] + jblk;
               if (phblk2cpu_loc[jhblk] == myid) {
                  int ibs_j = pibsblk_loc[jblkgl];
                  CMvmSlv < _Int, _Flt, _FltVect >::MvmA ('+', pasub_loc[j], _x + ibs_j,
                                                          _ax + ibs_i);
               }
            }
         }
      }

// Wait for completetion of sends/recvs

      CMPIDataExchange::WaitAll (nrecvs_loc + nsends_loc, psndrcv_recvs_loc,
                                 psndrcv_stats_loc);

// Perform remaining multiplications

      int jj, jj2, iblk;

      for (i = 0; i < nblks_recvs_loc; i++) {
         jblkgl = plistblk_recvs_loc[i];
         for (j = piablk_recvs_loc[i]; j < piablk_recvs_loc[i + 1]; j++) {
            ind = pind_recvs_loc[j];
            jj = pja_recvs_loc[ind * 2];
            jj2 = pja_recvs_loc[ind * 2 + 1];
            if (jj2 != jblkgl) {
               throw
                  " CTMvmSlvPar<_Int,_Flt,_FltVect>::MvmA: error: incorrect block number !";
            }
            px_temp_loc[jj] = px_recv_loc[ind];
         }
         for (j = pialist_recvs_loc[i]; j < pialist_recvs_loc[i + 1]; j++) {
            iblk = pjatriples_recvs_loc[j * 3];
            ihblk = pjatriples_recvs_loc[j * 3 + 1];
            ind = pjatriples_recvs_loc[j * 3 + 2];
            iblkgl = phblk2blks_loc[ihblk] + iblk;
            ibs_i = pibsblk_loc[iblkgl];
            CMatrix < _Int, _Flt > *pasub_loc = phmatr_loc[ihblk].GetASubArr ();
            CMatrix < int, float >*phmatr_str_loc = phmatr_loc[ihblk].GetHMatrStr ();
            int *pja_hblk = phmatr_str_loc->GetJaArr ();
            int *pja2_hblk = phmatr_str_loc->GetJa2Arr ();
            jj = pja_hblk[ind];
            jj2 = pja2_hblk[ind];
            if (phblk2blks_loc[jj2] + jj != jblkgl) {
               throw
                  " CTMvmSlvPar<_Int,_Flt,_FltVect>::MvmA: error 2: incorrect block number !";
            }
            CMvmSlv < _Int, _Flt, _FltVect >::MvmA ('+', pasub_loc[ind], px_temp_loc,
                                                    _ax + ibs_i);
         }
      }

      CMPIDataExchange::DeleteRecvs (psndrcv_recvs_loc);
      CMPIDataExchange::DeleteStats (psndrcv_stats_loc);

   }

// Clean MvmA structure
//========================================================================================
   template < typename _Int, typename _Flt > void CMvmParThreads_base < _Int,
      _Flt >::Clean_base ()
   {

      pcomm = NULL;
      nhblks = 0;
      phblk2cpu = NULL;
      phblk2blks = NULL;
      pblk2hblks = NULL;
      phblks = NULL;
      nblks = 0;
      pblks = NULL;
      phmatr = NULL;
      ni_cpu = 0;
      vector < int >ibsblk_dummy;
      ibsblk.swap (ibsblk_dummy);
      nlisthblk_own = 0;
      vector < int >listhblk_own_dummy;
      listhblk_own.swap (listhblk_own_dummy);
      nlistblk_own = 0;
      vector < int >listblk_own_dummy;
      listblk_own.swap (listblk_own_dummy);
      nsends = 0;
      vector < int >snd2cpu_dummy;
      vector < int >ia_sends_dummy;
      vector < int >ind_sends_dummy;
      snd2cpu.swap (snd2cpu_dummy);
      ia_sends.swap (ia_sends_dummy);
      ind_sends.swap (ind_sends_dummy);
      nrecvs = 0;
      vector < int >rcv2cpu_dummy;
      vector < int >ia_recvs_dummy;
      vector < int >ja_recvs_dummy;
      rcv2cpu.swap (rcv2cpu_dummy);
      ia_recvs.swap (ia_recvs_dummy);
      ja_recvs.swap (ja_recvs_dummy);
      nblks_recvs = 0;
      vector < int >listblk_recvs_dummy;
      vector < int >ialist_recvs_dummy;
      vector < int >jatriples_recvs_dummy;
      vector < int >iablk_recvs_dummy;
      vector < int >ind_recvs_dummy;
      listblk_recvs.swap (listblk_recvs_dummy);
      ialist_recvs.swap (ialist_recvs_dummy);
      jatriples_recvs.swap (jatriples_recvs_dummy);
      iablk_recvs.swap (iablk_recvs_dummy);
      ind_recvs.swap (ind_recvs_dummy);

      i_size_x_send = 0;
      i_size_x_recv = 0;
      i_size_x_temp = 0;

   }

// Clean MvmA structure
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CMvmParThreads < _Int, _Flt, _FltVect >::Clean ()
   {

      this->Clean_base ();

      CVectorData < _FltVect > x_send_dummy;
      CVectorData < _FltVect > x_recv_dummy;
      CVectorData < _FltVect > x_temp_dummy;

      x_send.swap (x_send_dummy);
      x_recv.swap (x_recv_dummy);
      x_temp.swap (x_temp_dummy);

   }

// Init SolveLU data
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CSlvParThreads < _Int, _Flt,
      _FltVect >::InitSolveLU (long long *_hblks_ext, vector < int >*_listpairs_ext,
                               int *_nblks_ext_arr, vector < int >*_blksnum_ext_arr,
                               vector < long long >*_blks_ext_arr, CTree * _tree_arr,
                               int *_nblks_ilu2_arr, int *_nblks1_ilu2_arr,
                               int *_nblks2_ilu2_arr, vector < long long >*_blks_ilu2_arr,
                               CVectorData < int >*_orderLU, CBMatrix < _Int,
                               _Flt > *_matrL, CBMatrix < _Int, _Flt > *_matrU)
   {

// Get control data

      void *pcomm_loc = this->pcomm;

      int nhblks_loc = this->nhblks;
      long long *phblks_loc = this->phblks;
      int *phblk2cpu_loc = this->phblk2cpu;
      int *pblk2hblks_loc = this->pblk2hblks;
      int *pibsblk_loc = &this->ibsblk[0];
      int nlisthblk_own_loc = this->nlisthblk_own;
      int *plisthblk_own_loc = &this->listhblk_own[0];

      int nproc = CMPIDataExchange::GetNproc (pcomm_loc);
      int myid = CMPIDataExchange::GetMyid (pcomm_loc);

// Store pointers

      this->phblks_ext = _hblks_ext;
      this->pnblks_ext_arr = _nblks_ext_arr;
      this->pblksnum_ext_arr = _blksnum_ext_arr;
      this->pblks_ext_arr = _blks_ext_arr;
      this->ptree_arr = _tree_arr;
      this->pnblks_ilu2_arr = _nblks_ilu2_arr;
      this->pnblks1_ilu2_arr = _nblks1_ilu2_arr;
      this->pnblks2_ilu2_arr = _nblks2_ilu2_arr;
      this->pblks_ilu2_arr = _blks_ilu2_arr;

      long long *phblks_ext_loc = this->phblks_ext;

// Compute maximal block size

      int nimax_ext = 0;

      int i, niloc;

      for (i = 0; i < nhblks_loc; i++) {
         niloc = (int) (phblks_ext_loc[i + 1] - phblks_ext_loc[i]);
         if (niloc > nimax_ext)
            nimax_ext = niloc;
      }

// Store pointers to other computed data

      this->porderLU = _orderLU;
      this->plistpairs_ext = _listpairs_ext;
      this->pmatrL = _matrL;
      this->pmatrU = _matrU;

// Create strLT structures

      this->strL_T.resize (nhblks_loc + 1);

      CMatrix < int, float >*pstrL_T = &this->strL_T[0];

      for (i = 0; i < nhblks_loc; i++) {
         if (phblk2cpu_loc[i] == myid) {
            CMatrix < int, float >*phmatr_temp = _matrL[i].GetHMatrStr ();
            int nlist_temp = phmatr_temp->GetNlist ();
            int nzja_temp = phmatr_temp->GetNzja ();
            int *pia_temp = phmatr_temp->GetIaArr ();

            vector < int >indHt2H (nzja_temp + 1);

            int *pindHt2H = &indHt2H[0];

            int j;

            {

               CVectorData < int >imask (5 * nlist_temp + 1);
               int *pimask = imask.Ptr ();

               for (j = 0; j < nlist_temp; j++)
                  pimask[j] = -1;

               int icycle = -1;

               phmatr_temp->TransposedSparsityListSp (icycle, pimask, pimask + nlist_temp,
                                                      pimask + 2 * nlist_temp,
                                                      pimask + 3 * nlist_temp,
                                                      pimask + 4 * nlist_temp,
                                                      pstrL_T[i]);

               int *pia_hmatrT = pstrL_T[i].GetIaArr ();
               int *pja_hmatrT = pstrL_T[i].GetJaArr ();

               for (j = 0; j < nlist_temp; j++)
                  pimask[j] = pia_temp[j];

               {

                  int ii, jj, k;

                  for (ii = 0; ii < nlist_temp; ii++) {
                     for (j = pia_hmatrT[ii]; j < pia_hmatrT[ii + 1]; j++) {
                        jj = pja_hmatrT[j];
                        k = pimask[jj];
                        pindHt2H[j] = k;
                        pimask[jj]++;
                     }
                  }

               }

            }

            pstrL_T[i].SetNzja2 (nzja_temp);
            pstrL_T[i].ResizeJa2 (nzja_temp);

            int *pja_hmatrT = pstrL_T[i].GetJa2Arr ();

            for (j = 0; j < nzja_temp; j++)
               pja_hmatrT[j] = pindHt2H[j];

         }
      }

// Create cpu recv data

      vector < int >imaskcpu (nproc);
      vector < int >listcpu (nproc);
      vector < int >indcpu (nproc);

      int *pimaskcpu = &imaskcpu[0];
      int *plistcpu = &listcpu[0];
      int *pindcpu = &indcpu[0];

      for (i = 0; i < nproc; i++)
         pimaskcpu[i] = -1;

      int icyclecpu = -1;

      icyclecpu++;

      int nlistcpu = 0;

      int jcpu, ihblk, jj2, jhblk, ni_ext, j;

      for (i = 0; i < nlisthblk_own_loc; i++) {
         ihblk = plisthblk_own_loc[i];
         ni_ext =
            (int) ((phblks_ext_loc[ihblk + 1] - phblks_ext_loc[ihblk]) -
                   (phblks_loc[ihblk + 1] - phblks_loc[ihblk]));
         int *pplistpairs_ext = &(_listpairs_ext[ihblk][0]);
         for (j = 0; j < ni_ext; j++) {
            jj2 = pplistpairs_ext[2 * j + 1];
            jhblk = pblk2hblks_loc[jj2];
            jcpu = phblk2cpu_loc[jhblk];
            if (jcpu != myid) {
               if (pimaskcpu[jcpu] != icyclecpu) {
                  plistcpu[nlistcpu] = jcpu;
                  nlistcpu++;
                  pimaskcpu[jcpu] = icyclecpu;
               }
            }
         }
      }

      sort (plistcpu, plistcpu + nlistcpu);

      for (i = 0; i < nlistcpu; i++) {
         jcpu = plistcpu[i];
         pindcpu[jcpu] = i;
      }

      this->nrecvs = nlistcpu;

      this->rcv2cpu.resize (nlistcpu + 1);
      this->ia_recvs.resize (nlistcpu + 1);

      int *prcv2cpu = &this->rcv2cpu[0];
      int *pia_recvs = &this->ia_recvs[0];

      for (i = 0; i < nlistcpu; i++)
         prcv2cpu[i] = plistcpu[i];
      for (i = 0; i <= nlistcpu; i++)
         pia_recvs[i] = 0;

      int ind;

      for (i = 0; i < nlisthblk_own_loc; i++) {
         ihblk = plisthblk_own_loc[i];
         ni_ext =
            (int) ((phblks_ext_loc[ihblk + 1] - phblks_ext_loc[ihblk]) -
                   (phblks_loc[ihblk + 1] - phblks_loc[ihblk]));
         int *pplistpairs_ext = &(_listpairs_ext[ihblk][0]);
         for (j = 0; j < ni_ext; j++) {
            jj2 = pplistpairs_ext[2 * j + 1];
            jhblk = pblk2hblks_loc[jj2];
            jcpu = phblk2cpu_loc[jhblk];
            if (jcpu != myid) {
               ind = pindcpu[jcpu];
               pia_recvs[ind + 1]++;
            }
         }
      }

      for (i = 0; i < nlistcpu; i++)
         pia_recvs[i + 1] = pia_recvs[i] + pia_recvs[i + 1];

      int nz_recvs = pia_recvs[nlistcpu];

      this->x_recv.resize (nz_recvs + 1);

      this->ialist_recvs.resize (nlisthblk_own_loc + 1);
      this->jalist_recvs.resize (2 * nz_recvs + 1);

      int *pialist_recvs = &this->ialist_recvs[0];
      int *pjalist_recvs = &this->jalist_recvs[0];

      nz_recvs = 0;

      for (i = 0; i < nlistcpu; i++)
         plistcpu[i] = pia_recvs[i];

      pialist_recvs[0] = 0;

      int k;

      for (i = 0; i < nlisthblk_own_loc; i++) {
         ihblk = plisthblk_own_loc[i];
         ni_ext =
            (int) ((phblks_ext_loc[ihblk + 1] - phblks_ext_loc[ihblk]) -
                   (phblks_loc[ihblk + 1] - phblks_loc[ihblk]));
         int *pplistpairs_ext = &(_listpairs_ext[ihblk][0]);
         for (j = 0; j < ni_ext; j++) {
            jj2 = pplistpairs_ext[2 * j + 1];
            jhblk = pblk2hblks_loc[jj2];
            jcpu = phblk2cpu_loc[jhblk];
            if (jcpu != myid) {
               ind = pindcpu[jcpu];
               k = plistcpu[ind];
               pjalist_recvs[nz_recvs * 2] = k;
               pjalist_recvs[nz_recvs * 2 + 1] = j;
               nz_recvs++;
               plistcpu[ind]++;
            }
         }
         pialist_recvs[i + 1] = nz_recvs;
      }

// Prepare exchange data

      vector < CBMatrix < int, _Flt > >hblk_send (nlistcpu + 1);

      CBMatrix < int, _Flt > *phblk_send = &hblk_send[0];

      for (i = 0; i < nlistcpu; i++) {
         phblk_send[i].SetNzblk (1);
         phblk_send[i].ResizeASub (1);
         CMatrix < int, _Flt > *pA_sub = phblk_send[i].GetASubArr ();
         CMatrix < int, _Flt > ablk_temp;
         int niloc = pia_recvs[i + 1] - pia_recvs[i];
         ablk_temp.ResizeAndSetAllSp (0, 0, niloc * 2, 0);
         pA_sub->ReplaceFree (ablk_temp);
      }

      nz_recvs = 0;

      for (i = 0; i < nlistcpu; i++)
         plistcpu[i] = pia_recvs[i];

      int jj;

      for (i = 0; i < nlisthblk_own_loc; i++) {
         ihblk = plisthblk_own_loc[i];
         ni_ext =
            (int) ((phblks_ext_loc[ihblk + 1] - phblks_ext_loc[ihblk]) -
                   (phblks_loc[ihblk + 1] - phblks_loc[ihblk]));
         int *pplistpairs_ext = &(_listpairs_ext[ihblk][0]);
         for (j = 0; j < ni_ext; j++) {
            jj = pplistpairs_ext[2 * j];
            jj2 = pplistpairs_ext[2 * j + 1];
            jhblk = pblk2hblks_loc[jj2];
            jcpu = phblk2cpu_loc[jhblk];
            if (jcpu != myid) {
               ind = pindcpu[jcpu];
               CMatrix < int, _Flt > *pA_sub = phblk_send[ind].GetASubArr ();
               int *pjaloc = pA_sub->GetJaArr ();
               k = plistcpu[ind] - pia_recvs[ind];
               pjaloc[k * 2] = jj;
               pjaloc[k * 2 + 1] = jj2;
               plistcpu[ind]++;
            }
         }
      }

// Pack send data

      vector < int >CpuIDSend (nlistcpu);
      vector < vector < char > >ObjSend (nlistcpu);

      int *pCpuIDSend = NULL;
      vector < char >*pObjSend = NULL;

      if (nlistcpu > 0) {
         pCpuIDSend = &CpuIDSend[0];
         pObjSend = &ObjSend[0];
      }

      long long isize;
      char *pobj;

      for (i = 0; i < nlistcpu; i++) {
         pCpuIDSend[i] = prcv2cpu[i];
         isize = phblk_send[i].GetPackedSize ();
         pObjSend[i].resize ((size_t) isize);
         pobj = &(pObjSend[i][0]);
         phblk_send[i].FillPacked (isize, pobj);
         phblk_send[i].Clean ();
      }

// Exchange

      vector < int >CpuIDRecv;
      vector < vector < char > >ObjRecv;

      CMPIDataExchange::DataExchange (pcomm_loc, CpuIDSend, ObjSend, CpuIDRecv, ObjRecv);

// Free send data

      {
         vector < int >CpuIDSend_temp;
         vector < vector < char > >ObjSend_temp;
         CpuIDSend.swap (CpuIDSend_temp);
         ObjSend.swap (ObjSend_temp);
      }

// Unpack receive data

      int nrecv_loc = (int) CpuIDRecv.size ();

      vector < char >*pObjRecv = NULL;
      int *pCpuIDRecv = NULL;

      if (nrecv_loc > 0) {
         pObjRecv = &ObjRecv[0];
         pCpuIDRecv = &CpuIDRecv[0];
      }

      vector < CBMatrix < int, _Flt > >hblk_recv (nrecv_loc + 1);

      CBMatrix < int, _Flt > *phblk_recv = &hblk_recv[0];

      for (i = 0; i < nrecv_loc; i++) {
         isize = (long long) pObjRecv[i].size ();
         pobj = &(pObjRecv[i][0]);
         phblk_recv[i].UnPack (isize, pobj);
      }

// Free recv data

      {
         vector < vector < char > >ObjRecv_temp;
         ObjRecv.swap (ObjRecv_temp);
      }

// Compute correct ordering of cpu data

      vector < CSortInt > iiarr (nrecv_loc + 1);
      CSortInt *piiarr = &iiarr[0];

      for (i = 0; i < nrecv_loc; i++) {
         piiarr[i].ival = pCpuIDRecv[i];
         piiarr[i].i2val = i;
      }

      sort (piiarr, piiarr + nrecv_loc);

// Store received data

      this->nsends = nrecv_loc;

      this->snd2cpu.resize (nrecv_loc + 1);
      this->ia_sends.resize (nrecv_loc + 1);

      int *psnd2cpu = &this->snd2cpu[0];
      int *pia_sends = &this->ia_sends[0];

      int nz_sends = 0;

      pia_sends[0] = 0;

      for (i = 0; i < nrecv_loc; i++) {
         ind = piiarr[i].i2val;
         psnd2cpu[i] = piiarr[i].ival;
         CMatrix < int, _Flt > *pA_sub = phblk_recv[ind].GetASubArr ();
         int nzjaloc = pA_sub->GetNzja () / 2;
         nz_sends += nzjaloc;
         pia_sends[i + 1] = nz_sends;
      }

      this->ind_sends.resize (nz_sends + 1);
      this->x_send.resize (nz_sends + 1);

      int *pind_sends = &this->ind_sends[0];

      nz_sends = 0;

      int ibs;

      for (i = 0; i < nrecv_loc; i++) {
         ind = piiarr[i].i2val;
         CMatrix < int, _Flt > *pA_sub = phblk_recv[ind].GetASubArr ();
         int nzjaloc = pA_sub->GetNzja () / 2;
         int *pjaloc = pA_sub->GetJaArr ();
         for (j = 0; j < nzjaloc; j++) {
            jj = pjaloc[j * 2];
            jj2 = pjaloc[j * 2 + 1];
            ibs = pibsblk_loc[jj2];
            pind_sends[nz_sends] = ibs + jj;
            nz_sends++;
         }
      }

      this->x_temp.resize (nimax_ext * 2 + 1);

   }

// Perform SolveLU computations
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CSlvParThreads < _Int, _Flt,
      _FltVect >::SolveLU (const _FltVect * _x, _FltVect * _px)
   {

// Open mvm structure

      void *pcomm_loc = this->pcomm;

      long long *phblks_loc = this->phblks;
      long long *phblks_ext_loc = this->phblks_ext;
      int *phblk2cpu_loc = this->phblk2cpu;
      int *phblk2blks_loc = this->phblk2blks;
      int *pblk2hblks_loc = this->pblk2hblks;
      int ni_cpu_loc = this->ni_cpu;
      int *pibsblk_loc = &this->ibsblk[0];
      int nlisthblk_own_loc = this->nlisthblk_own;
      int *plisthblk_own_loc = &this->listhblk_own[0];
      vector < int >*plistpairs_ext_loc = this->plistpairs_ext;
      int nsends_loc = this->nsends;
      int *psnd2cpu_loc = &this->snd2cpu[0];
      int *pia_sends_loc = &this->ia_sends[0];
      int *pind_sends_loc = &this->ind_sends[0];
      _FltVect *px_send_loc = this->x_send.Ptr ();
      int nrecvs_loc = this->nrecvs;
      int *prcv2cpu_loc = &this->rcv2cpu[0];
      int *pia_recvs_loc = &this->ia_recvs[0];
      _FltVect *px_recv_loc = this->x_recv.Ptr ();
      int *pialist_recvs_loc = &this->ialist_recvs[0];
      int *pjalist_recvs_loc = &this->jalist_recvs[0];
      _FltVect *px_temp_loc = this->x_temp.Ptr ();

      int myid = CMPIDataExchange::GetMyid (pcomm_loc);

// Init array ax by zeroes

      CVector < _FltVect >::SetByZeroes_thr (ni_cpu_loc, _px);

// Prepare send

      int ni_send_loc = pia_sends_loc[nsends_loc];

      int i;

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      {
         int ni_part = ni_send_loc / n_thr;
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {
            int ibeg = ipar * ni_part;
            int iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = ni_send_loc - 1;
            for (int i = ibeg; i <= iend; i++) {
               int ind = pind_sends_loc[i];
               px_send_loc[i] = _x[ind];
            }
         }
      }

      int ni_recv_loc = pia_recvs_loc[nrecvs_loc];

      CVector < _FltVect >::SetByZeroes_thr (ni_recv_loc, px_recv_loc);

// Init async recvs and sends

      void *psndrcv_recvs_loc;
      void *psndrcv_stats_loc;

      CMPIDataExchange::AllocateRecvs (nrecvs_loc + nsends_loc, psndrcv_recvs_loc);
      CMPIDataExchange::AllocateStats (nrecvs_loc + nsends_loc, psndrcv_stats_loc);

      int icpu, isize, ibs;

      for (i = 0; i < nrecvs_loc; i++) {
         icpu = prcv2cpu_loc[i];
         isize = (pia_recvs_loc[i + 1] - pia_recvs_loc[i]) * sizeof (_FltVect);
         ibs = pia_recvs_loc[i];
         CMPIDataExchange::IRecv (pcomm_loc, icpu, myid, isize,
                                  (char *) (px_recv_loc + ibs), i, psndrcv_recvs_loc);
      }

      for (i = 0; i < nsends_loc; i++) {
         icpu = psnd2cpu_loc[i];
         isize = (pia_sends_loc[i + 1] - pia_sends_loc[i]) * sizeof (_FltVect);
         ibs = pia_sends_loc[i];
         CMPIDataExchange::ISend (pcomm_loc, icpu, icpu, isize,
                                  (char *) (px_send_loc + ibs), i + nrecvs_loc,
                                  psndrcv_recvs_loc);
      }

// Wait for completetion of sends/recvs

      CMPIDataExchange::WaitAll (nrecvs_loc + nsends_loc, psndrcv_recvs_loc,
                                 psndrcv_stats_loc);

      CVector < _FltVect >::SetByZeroes_thr (ni_send_loc, px_send_loc);

// Perform local computations

      int iblk0, ihblk, ibs_i, niloc, niextloc, ni_ini;

      _FltVect *px1_temp;
      _FltVect *px2_temp;

      for (i = 0; i < nlisthblk_own_loc; i++) {
         ihblk = plisthblk_own_loc[i];
         iblk0 = phblk2blks_loc[ihblk];
         ibs_i = pibsblk_loc[iblk0];
         niloc = (int) (phblks_loc[ihblk + 1] - phblks_loc[ihblk]);
         niextloc = (int) (phblks_ext_loc[ihblk + 1] - phblks_ext_loc[ihblk]);
         ni_ini = niextloc - niloc;
         px1_temp = px_temp_loc;
         px2_temp = px1_temp + niextloc;
         CVector < _FltVect >::SetByZeroes_thr (ni_ini, px1_temp);
         CVector < _FltVect >::CopyVector_thr (niloc, _x + ibs_i, px1_temp + ni_ini);
         int *pplistpairs_loc = &(plistpairs_ext_loc[ihblk][0]);
         {
            int ni_part = ni_ini / n_thr;
#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < n_thr; ipar++) {
               int ibeg = ipar * ni_part;
               int iend = (ipar + 1) * ni_part - 1;
               if (ipar == n_thr - 1)
                  iend = ni_ini - 1;
               for (int j = ibeg; j <= iend; j++) {
                  int jj = pplistpairs_loc[j * 2];
                  int jj2 = pplistpairs_loc[j * 2 + 1];
                  int jhblk = pblk2hblks_loc[jj2];
                  if (phblk2cpu_loc[jhblk] == myid) {
                     int ibs_j = pibsblk_loc[jj2];
                     px1_temp[j] = _x[ibs_j + jj];
                  }
               }
            }
         }
         {
            int niloc = pialist_recvs_loc[i + 1] - pialist_recvs_loc[i];
            int ishift = pialist_recvs_loc[i];
            int ni_part = niloc / n_thr;
#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < n_thr; ipar++) {
               int ibeg = ishift + ipar * ni_part;
               int iend = ishift + (ipar + 1) * ni_part - 1;
               if (ipar == n_thr - 1)
                  iend = ishift + niloc - 1;
               for (int j = ibeg; j <= iend; j++) {
                  int ind = pjalist_recvs_loc[2 * j];
                  int ind1 = pjalist_recvs_loc[2 * j + 1];
                  px1_temp[ind1] = px_recv_loc[ind];
               }
            }
         }
         {
            int *pporderLU = this->porderLU[ihblk].Ptr ();
            CVector < _FltVect >::OrderVector_thr (niextloc, pporderLU, px1_temp,
                                                   px2_temp);
            this->SolveL (ihblk, px2_temp);
            CVector < _FltVect >::SetByZeroes_thr (ni_ini, px2_temp);
            this->SolveU (ihblk, px2_temp);
            CVector < _FltVect >::InvOrderVector_thr (niextloc, pporderLU, px2_temp,
                                                      px1_temp);
         }
         {
            int niloc = pialist_recvs_loc[i + 1] - pialist_recvs_loc[i];
            int ishift = pialist_recvs_loc[i];
            int ni_part = niloc / n_thr;
#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < n_thr; ipar++) {
               int ibeg = ishift + ipar * ni_part;
               int iend = ishift + (ipar + 1) * ni_part - 1;
               if (ipar == n_thr - 1)
                  iend = ishift + niloc - 1;
               for (int j = ibeg; j <= iend; j++) {
                  int ind = pjalist_recvs_loc[2 * j];
                  int ind1 = pjalist_recvs_loc[2 * j + 1];
                  px_recv_loc[ind] = px1_temp[ind1];
               }
            }
         }
         {
            int ni_part = ni_ini / n_thr;
#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < n_thr; ipar++) {
               int ibeg = ipar * ni_part;
               int iend = (ipar + 1) * ni_part - 1;
               if (ipar == n_thr - 1)
                  iend = ni_ini - 1;
               for (int j = ibeg; j <= iend; j++) {
                  int jj = pplistpairs_loc[j * 2];
                  int jj2 = pplistpairs_loc[j * 2 + 1];
                  int jhblk = pblk2hblks_loc[jj2];
                  if (phblk2cpu_loc[jhblk] == myid) {
                     int ibs_j = pibsblk_loc[jj2];
                     _px[ibs_j + jj] += px1_temp[j];
                  }
               }
            }
         }
         CVector < _FltVect >::AddReplaceVector_thr (niloc, px1_temp + ni_ini,
                                                     _px + ibs_i);
      }

// Backward exchanges

      for (i = 0; i < nsends_loc; i++) {
         icpu = psnd2cpu_loc[i];
         isize = (pia_sends_loc[i + 1] - pia_sends_loc[i]) * sizeof (_FltVect);
         ibs = pia_sends_loc[i];
         CMPIDataExchange::IRecv (pcomm_loc, icpu, myid, isize,
                                  (char *) (px_send_loc + ibs), i + nrecvs_loc,
                                  psndrcv_recvs_loc);
      }

      for (i = 0; i < nrecvs_loc; i++) {
         icpu = prcv2cpu_loc[i];
         isize = (pia_recvs_loc[i + 1] - pia_recvs_loc[i]) * sizeof (_FltVect);
         ibs = pia_recvs_loc[i];
         CMPIDataExchange::ISend (pcomm_loc, icpu, icpu, isize,
                                  (char *) (px_recv_loc + ibs), i, psndrcv_recvs_loc);
      }

// Wait for completetion of sends/recvs

      CMPIDataExchange::WaitAll (nrecvs_loc + nsends_loc, psndrcv_recvs_loc,
                                 psndrcv_stats_loc);

// Finally update received data

      {
         int ni_part = ni_send_loc / n_thr;
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {
            int ibeg = ipar * ni_part;
            int iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = ni_send_loc - 1;
            for (int i = ibeg; i <= iend; i++) {
               int ind = pind_sends_loc[i];
               _px[ind] += px_send_loc[i];
            }
         }
      }

      CMPIDataExchange::DeleteRecvs (psndrcv_recvs_loc);
      CMPIDataExchange::DeleteStats (psndrcv_stats_loc);

   }

// In-place perform SolveU computations
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CSlvParThreads < _Int, _Flt,
      _FltVect >::SolveU (int _ihblk, _FltVect * _px)
   {

// Open mvm structure

      CTree *ptree_arr_loc = this->ptree_arr;
      int *pnblks_ilu2_arr_loc = this->pnblks_ilu2_arr;
      int *pnblks1_ilu2_arr_loc = this->pnblks1_ilu2_arr;
      int *pnblks2_ilu2_arr_loc = this->pnblks2_ilu2_arr;
      vector < long long >*pblks_ilu2_arr_loc = this->pblks_ilu2_arr;
      CBMatrix < _Int, _Flt > *pmatrU_loc = this->pmatrU;

// Perform U solve according to the three tree's

      int nblks_curr = pnblks_ilu2_arr_loc[_ihblk];
      int nblks1_curr = pnblks1_ilu2_arr_loc[_ihblk];
      int nblks2_curr = pnblks2_ilu2_arr_loc[_ihblk];
      long long *ppblks_curr = &pblks_ilu2_arr_loc[_ihblk][0];

      CMatrix < int, float >*phmatr = pmatrU_loc[_ihblk].GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub = pmatrU_loc[_ihblk].GetASubArr ();

      int *pia_hmatr = phmatr->GetIaArr ();
      int *pja_hmatr = phmatr->GetJaArr ();

      int nblks12_curr = nblks1_curr + nblks2_curr;

      int nblks3_curr = nblks_curr - nblks12_curr;

      int ilev;

// Last tree

      if (nblks3_curr > 0) {

         int nlev_3 = ptree_arr_loc[_ihblk * 3 + 2].GetNlev ();
         int *pnnodes_lev_3 = ptree_arr_loc[_ihblk * 3 + 2].GetNNodesLev ();
         vector < int >*pnodeslevlist_3 =
            ptree_arr_loc[_ihblk * 3 + 2].GetNodesLevList ();

// Run

         int nnodes_curr_3;
         int *ppnodeslevlist_3;

         for (ilev = 0; ilev < nlev_3; ilev++) {

            nnodes_curr_3 = pnnodes_lev_3[ilev];
            ppnodeslevlist_3 = &pnodeslevlist_3[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < nnodes_curr_3; ipar++) {
               int iblk = ppnodeslevlist_3[ipar];
               int iblkgl = iblk + nblks12_curr;
               int ibs_i = (int) ppblks_curr[iblkgl];
               if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
                  for (int j = pia_hmatr[iblkgl + 1] - 1; j >= pia_hmatr[iblkgl] + 1; j--) {
                     int jblkgl = pja_hmatr[j];
                     int ibs_j = (int) ppblks_curr[jblkgl];
                     CMvmSlv < _Int, _Flt, _FltVect >::MvmA ('-', pA_sub[j], _px + ibs_j,
                                                             _px + ibs_i);
                  }
                  int j = pia_hmatr[iblkgl];
                  CMvmSlv < _Int, _Flt, _FltVect >::SolveU (pA_sub[j], _px + ibs_i,
                                                            _px + ibs_i);
               }
            }

         }

      }
// First and second tree

      int nlev_1 = ptree_arr_loc[_ihblk * 3].GetNlev ();
      int *pnnodes_lev_1 = ptree_arr_loc[_ihblk * 3].GetNNodesLev ();
      vector < int >*pnodeslevlist_1 = ptree_arr_loc[_ihblk * 3].GetNodesLevList ();

      int nlev_2 = ptree_arr_loc[_ihblk * 3 + 1].GetNlev ();
      int *pnnodes_lev_2 = ptree_arr_loc[_ihblk * 3 + 1].GetNNodesLev ();
      vector < int >*pnodeslevlist_2 = ptree_arr_loc[_ihblk * 3 + 1].GetNodesLevList ();

      int nlev_12_min = nlev_1;
      if (nlev_2 < nlev_12_min)
         nlev_12_min = nlev_2;

      int nnodes_curr_1, nnodes_curr_2;
      int *ppnodeslevlist_1;
      int *ppnodeslevlist_2;

// Run

      for (ilev = 0; ilev < nlev_12_min; ilev++) {

         nnodes_curr_1 = pnnodes_lev_1[ilev];
         ppnodeslevlist_1 = &pnodeslevlist_1[ilev][0];

         nnodes_curr_2 = pnnodes_lev_2[ilev];
         ppnodeslevlist_2 = &pnodeslevlist_2[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_1 + nnodes_curr_2; ipar++) {
            int iblkgl;
            if (ipar < nnodes_curr_1) {
               iblkgl = ppnodeslevlist_1[ipar];
            } else {
               iblkgl = nblks1_curr + ppnodeslevlist_2[ipar - nnodes_curr_1];
            }
            int ibs_i = (int) ppblks_curr[iblkgl];
            if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
               for (int j = pia_hmatr[iblkgl + 1] - 1; j >= pia_hmatr[iblkgl] + 1; j--) {
                  int jblkgl = pja_hmatr[j];
                  int ibs_j = (int) ppblks_curr[jblkgl];
                  CMvmSlv < _Int, _Flt, _FltVect >::MvmA ('-', pA_sub[j], _px + ibs_j,
                                                          _px + ibs_i);
               }
               int j = pia_hmatr[iblkgl];
               CMvmSlv < _Int, _Flt, _FltVect >::SolveU (pA_sub[j], _px + ibs_i,
                                                         _px + ibs_i);
            }
         }

      }

      for (ilev = nlev_12_min; ilev < nlev_1; ilev++) {

         nnodes_curr_1 = pnnodes_lev_1[ilev];
         ppnodeslevlist_1 = &pnodeslevlist_1[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_1; ipar++) {
            int iblkgl = ppnodeslevlist_1[ipar];
            int ibs_i = (int) ppblks_curr[iblkgl];
            if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
               for (int j = pia_hmatr[iblkgl + 1] - 1; j >= pia_hmatr[iblkgl] + 1; j--) {
                  int jblkgl = pja_hmatr[j];
                  int ibs_j = (int) ppblks_curr[jblkgl];
                  CMvmSlv < _Int, _Flt, _FltVect >::MvmA ('-', pA_sub[j], _px + ibs_j,
                                                          _px + ibs_i);
               }
               int j = pia_hmatr[iblkgl];
               CMvmSlv < _Int, _Flt, _FltVect >::SolveU (pA_sub[j], _px + ibs_i,
                                                         _px + ibs_i);
            }
         }

      }

      for (ilev = nlev_12_min; ilev < nlev_2; ilev++) {

         nnodes_curr_2 = pnnodes_lev_2[ilev];
         ppnodeslevlist_2 = &pnodeslevlist_2[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_2; ipar++) {
            int iblkgl = nblks1_curr + ppnodeslevlist_2[ipar];
            int ibs_i = (int) ppblks_curr[iblkgl];
            if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
               for (int j = pia_hmatr[iblkgl + 1] - 1; j >= pia_hmatr[iblkgl] + 1; j--) {
                  int jblkgl = pja_hmatr[j];
                  int ibs_j = (int) ppblks_curr[jblkgl];
                  CMvmSlv < _Int, _Flt, _FltVect >::MvmA ('-', pA_sub[j], _px + ibs_j,
                                                          _px + ibs_i);
               }
               int j = pia_hmatr[iblkgl];
               CMvmSlv < _Int, _Flt, _FltVect >::SolveU (pA_sub[j], _px + ibs_i,
                                                         _px + ibs_i);
            }
         }
      }
   }

// In-place perform SolveL computations
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CSlvParThreads < _Int, _Flt,
      _FltVect >::SolveL (int _ihblk, _FltVect * _px)
   {

// Open mvm structure

      CTree *ptree_arr_loc = this->ptree_arr;
      int *pnblks_ilu2_arr_loc = this->pnblks_ilu2_arr;
      int *pnblks1_ilu2_arr_loc = this->pnblks1_ilu2_arr;
      int *pnblks2_ilu2_arr_loc = this->pnblks2_ilu2_arr;
      vector < long long >*pblks_ilu2_arr_loc = this->pblks_ilu2_arr;
      CBMatrix < _Int, _Flt > *pmatrL_loc = this->pmatrL;
      CMatrix < int, float >*pstrL_T = &this->strL_T[0];

// Perform U solve according to the three tree's

      int nblks_curr = pnblks_ilu2_arr_loc[_ihblk];
      int nblks1_curr = pnblks1_ilu2_arr_loc[_ihblk];
      int nblks2_curr = pnblks2_ilu2_arr_loc[_ihblk];
      long long *ppblks_curr = &pblks_ilu2_arr_loc[_ihblk][0];

      CMatrix < _Int, _Flt > *pA_sub = pmatrL_loc[_ihblk].GetASubArr ();

      int *pia_hmatr = pstrL_T[_ihblk].GetIaArr ();
      int *pja_hmatr = pstrL_T[_ihblk].GetJaArr ();
      int *pind_lt2l = pstrL_T[_ihblk].GetJa2Arr ();

      int nblks12_curr = nblks1_curr + nblks2_curr;

      int nblks3_curr = nblks_curr - nblks12_curr;

      int ilev;

// First and second tree

      int nlev_1 = ptree_arr_loc[_ihblk * 3].GetNlev ();
      int *pnnodes_lev_1 = ptree_arr_loc[_ihblk * 3].GetNNodesLev ();
      vector < int >*pnodeslevlist_1 = ptree_arr_loc[_ihblk * 3].GetNodesLevList ();

      int nlev_2 = ptree_arr_loc[_ihblk * 3 + 1].GetNlev ();
      int *pnnodes_lev_2 = ptree_arr_loc[_ihblk * 3 + 1].GetNNodesLev ();
      vector < int >*pnodeslevlist_2 = ptree_arr_loc[_ihblk * 3 + 1].GetNodesLevList ();

      int nlev_12_min = nlev_1;
      if (nlev_2 < nlev_12_min)
         nlev_12_min = nlev_2;

      int nnodes_curr_1, nnodes_curr_2;
      int *ppnodeslevlist_1;
      int *ppnodeslevlist_2;

// Run

      for (ilev = nlev_1 - 1; ilev >= nlev_12_min; ilev--) {

         nnodes_curr_1 = pnnodes_lev_1[ilev];
         ppnodeslevlist_1 = &pnodeslevlist_1[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_1; ipar++) {
            int iblkgl = ppnodeslevlist_1[ipar];
            int ibs_i = (int) ppblks_curr[iblkgl];
            if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
               for (int j = pia_hmatr[iblkgl]; j < pia_hmatr[iblkgl + 1] - 1; j++) {
                  int jblkgl = pja_hmatr[j];
                  int ind = pind_lt2l[j];
                  int ibs_j = (int) ppblks_curr[jblkgl];
                  CMvmSlv < _Int, _Flt, _FltVect >::MvmAT ('-', pA_sub[ind], _px + ibs_j,
                                                           _px + ibs_i);
               }
               int j = pia_hmatr[iblkgl + 1] - 1;
               int ind = pind_lt2l[j];
               CMvmSlv < _Int, _Flt, _FltVect >::SolveL (pA_sub[ind], _px + ibs_i,
                                                         _px + ibs_i);
            }
         }

      }

      for (ilev = nlev_2 - 1; ilev >= nlev_12_min; ilev--) {

         nnodes_curr_2 = pnnodes_lev_2[ilev];
         ppnodeslevlist_2 = &pnodeslevlist_2[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_2; ipar++) {
            int iblkgl = nblks1_curr + ppnodeslevlist_2[ipar];
            int ibs_i = (int) ppblks_curr[iblkgl];
            if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
               for (int j = pia_hmatr[iblkgl]; j < pia_hmatr[iblkgl + 1] - 1; j++) {
                  int jblkgl = pja_hmatr[j];
                  int ind = pind_lt2l[j];
                  int ibs_j = (int) ppblks_curr[jblkgl];
                  CMvmSlv < _Int, _Flt, _FltVect >::MvmAT ('-', pA_sub[ind], _px + ibs_j,
                                                           _px + ibs_i);
               }
               int j = pia_hmatr[iblkgl + 1] - 1;
               int ind = pind_lt2l[j];
               CMvmSlv < _Int, _Flt, _FltVect >::SolveL (pA_sub[ind], _px + ibs_i,
                                                         _px + ibs_i);
            }
         }

      }

      for (ilev = nlev_12_min - 1; ilev >= 0; ilev--) {

         nnodes_curr_1 = pnnodes_lev_1[ilev];
         ppnodeslevlist_1 = &pnodeslevlist_1[ilev][0];

         nnodes_curr_2 = pnnodes_lev_2[ilev];
         ppnodeslevlist_2 = &pnodeslevlist_2[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_1 + nnodes_curr_2; ipar++) {
            int iblkgl;
            if (ipar < nnodes_curr_1) {
               iblkgl = ppnodeslevlist_1[ipar];
            } else {
               iblkgl = nblks1_curr + ppnodeslevlist_2[ipar - nnodes_curr_1];
            }
            int ibs_i = (int) ppblks_curr[iblkgl];
            if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
               for (int j = pia_hmatr[iblkgl]; j < pia_hmatr[iblkgl + 1] - 1; j++) {
                  int jblkgl = pja_hmatr[j];
                  int ind = pind_lt2l[j];
                  int ibs_j = (int) ppblks_curr[jblkgl];
                  CMvmSlv < _Int, _Flt, _FltVect >::MvmAT ('-', pA_sub[ind], _px + ibs_j,
                                                           _px + ibs_i);
               }
               int j = pia_hmatr[iblkgl + 1] - 1;
               int ind = pind_lt2l[j];
               CMvmSlv < _Int, _Flt, _FltVect >::SolveL (pA_sub[ind], _px + ibs_i,
                                                         _px + ibs_i);
            }
         }

      }

// Last tree

      if (nblks3_curr > 0) {

         int nlev_3 = ptree_arr_loc[_ihblk * 3 + 2].GetNlev ();
         int *pnnodes_lev_3 = ptree_arr_loc[_ihblk * 3 + 2].GetNNodesLev ();
         vector < int >*pnodeslevlist_3 =
            ptree_arr_loc[_ihblk * 3 + 2].GetNodesLevList ();

// Run

         int nnodes_curr_3;
         int *ppnodeslevlist_3;

         for (ilev = nlev_3 - 1; ilev >= 0; ilev--) {

            nnodes_curr_3 = pnnodes_lev_3[ilev];
            ppnodeslevlist_3 = &pnodeslevlist_3[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < nnodes_curr_3; ipar++) {
               int iblk = ppnodeslevlist_3[ipar];
               int iblkgl = iblk + nblks12_curr;
               int ibs_i = (int) ppblks_curr[iblkgl];
               if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
                  for (int j = pia_hmatr[iblkgl]; j < pia_hmatr[iblkgl + 1] - 1; j++) {
                     int jblkgl = pja_hmatr[j];
                     int ind = pind_lt2l[j];
                     int ibs_j = (int) ppblks_curr[jblkgl];
                     CMvmSlv < _Int, _Flt, _FltVect >::MvmAT ('-', pA_sub[ind],
                                                              _px + ibs_j, _px + ibs_i);
                  }
                  int j = pia_hmatr[iblkgl + 1] - 1;
                  int ind = pind_lt2l[j];
                  CMvmSlv < _Int, _Flt, _FltVect >::SolveL (pA_sub[ind], _px + ibs_i,
                                                            _px + ibs_i);
               }
            }

         }

      }

   }

// Clean Slv structure
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CSlvParThreads < _Int, _Flt, _FltVect >::Clean ()
   {

      pcomm = NULL;
      nhblks = 0;
      phblks = NULL;
      phblks_ext = NULL;
      phblk2cpu = NULL;
      phblk2blks = NULL;
      pblk2hblks = NULL;
      pblks = NULL;
      ni_cpu = 0;
      vector < int >ibsblk_dummy;
      ibsblk.swap (ibsblk_dummy);
      nlisthblk_own = 0;
      vector < int >listhblk_own_dummy;
      listhblk_own.swap (listhblk_own_dummy);
      plistpairs_ext = NULL;
      pnblks_ext_arr = NULL;
      pblksnum_ext_arr = NULL;
      pblks_ext_arr = NULL;
      ptree_arr = NULL;
      pnblks_ilu2_arr = NULL;
      pnblks1_ilu2_arr = NULL;
      pnblks2_ilu2_arr = NULL;
      pblks_ilu2_arr = NULL;
      porderLU = NULL;
      pmatrL = NULL;
      pmatrU = NULL;
      nsends = 0;
      vector < int >snd2cpu_dummy;
      vector < int >ia_sends_dummy;
      vector < int >ind_sends_dummy;
      CVectorData < _FltVect > x_send_dummy;
      snd2cpu.swap (snd2cpu_dummy);
      ia_sends.swap (ia_sends_dummy);
      ind_sends.swap (ind_sends_dummy);
      x_send.swap (x_send_dummy);
      nrecvs = 0;
      vector < int >rcv2cpu_dummy;
      vector < int >ia_recvs_dummy;
      CVectorData < _FltVect > x_recv_dummy;
      vector < int >ialist_recvs_dummy;
      vector < int >jalist_recvs_dummy;
      CVectorData < _FltVect > x_temp_dummy;
      rcv2cpu.swap (rcv2cpu_dummy);
      ia_recvs.swap (ia_recvs_dummy);
      x_recv.swap (x_recv_dummy);
      ialist_recvs.swap (ialist_recvs_dummy);
      jalist_recvs.swap (jalist_recvs_dummy);
      x_temp.swap (x_temp_dummy);
   }

// Prepare matrix data only including preliminary repartitioning of the matrix
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::PrepareMatrixThreads (void *_pcomm, SParams * _params,
                                        SStatData * _stats, int _nhblks, int *_hblk2cpu,
                                        int *_hblk2blks_ini, long long *_hblks,
                                        int _nblks_ini, long long *_blks_ini, _Int * _ia,
                                        _Int * _ja, _Flt * _a)
   {

      _stats->dtime_part = 0.0e0;
      _stats->dtime_ord = 0.0e0;

//      ofstream ffout ("ChkDecomp_1x1.dat");

// Store control data if necessary

      void *pcommloc = _pcomm;

      if (_params->b_new_part) {

         this->pcomm = _pcomm;
         this->nhblks_ini = _nhblks;
         this->nhblks = _nhblks;
         this->nblks_ini = _nblks_ini;

         this->hblks_ini.resize (_nhblks + 1);
         this->hblks.resize (_nhblks + 1);
         this->hblk2cpu_ini.resize (_nhblks + 1);
         this->hblk2cpu.resize (_nhblks + 1);
         this->hblk2blks_ini.resize (_nhblks + 1);
         this->blks_ini.resize (_nblks_ini + 1);

         long long *phblks_ini = &this->hblks_ini[0];
         long long *phblks = &this->hblks[0];
         int *phblk2cpu_ini = &this->hblk2cpu_ini[0];
         int *phblk2cpu = &this->hblk2cpu[0];
         int *phblk2blks_ini = &this->hblk2blks_ini[0];
         long long *pblks_ini = &this->blks_ini[0];

         int i;

         for (i = 0; i <= _nhblks; i++) {
            phblks_ini[i] = _hblks[i];
            phblks[i] = _hblks[i];
         }

         for (i = 0; i < _nhblks; i++) {
            phblk2cpu_ini[i] = _hblk2cpu[i];
            phblk2cpu[i] = _hblk2cpu[i];
         }

         for (i = 0; i <= _nhblks; i++)
            phblk2blks_ini[i] = _hblk2blks_ini[i];

         for (i = 0; i <= _nblks_ini; i++)
            pblks_ini[i] = _blks_ini[i];

      }

      vector < int >blk2hblks_ini (_nblks_ini);
      int *pblk2hblks_ini = &blk2hblks_ini[0];

      {

         int i, j;

         for (i = 0; i < _nhblks; i++) {
            for (j = _hblk2blks_ini[i]; j < _hblk2blks_ini[i + 1]; j++) {
               pblk2hblks_ini[j] = i;
            }
         }

      }

      int myid = CMPIDataExchange::GetMyid (_pcomm);

      vector < CBMatrix < _Int, _Flt > >hmatr_ini_arr (_nhblks + 1);
      CBMatrix < _Int, _Flt > *phmatr_ini_arr = &hmatr_ini_arr[0];

// Create matrix as set of hblocks

      int nimax = 0;

      int i, niloc;

      for (i = 0; i < _nblks_ini; i++) {
         niloc = (int) (_blks_ini[i + 1] - _blks_ini[i]);
         if (niloc > nimax)
            nimax = niloc;
      }

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      vector < int >icycle_thr (n_thr + 1);
      vector < CVectorData < _Int > >listloc_thr (n_thr + 1);
      vector < CVectorData < int > >imaskblk_thr (n_thr + 1);

      int *picycle_thr = &icycle_thr[0];
      CVectorData < _Int > *plistloc_thr = &listloc_thr[0];
      CVectorData < int >*pimaskblk_thr = &imaskblk_thr[0];

      for (i = 0; i < n_thr; i++)
         picycle_thr[i] = -1;

      int ibeg = 0;

      for (i = 0; i < _nhblks; i++) {
         if (_hblk2cpu[i] == myid) {

            int niblk = _hblk2blks_ini[i + 1] - _hblk2blks_ini[i];
            vector < CBMatrix < _Int, _Flt > >hblk_arr (niblk + 1);
            CBMatrix < _Int, _Flt > *phblk_arr = &hblk_arr[0];
            int iblk0 = _hblk2blks_ini[i];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = _hblk2blks_ini[i]; ipar < _hblk2blks_ini[i + 1]; ipar++) {

               int my_thr = 0;
#ifdef USE_THREADS
               my_thr = omp_get_thread_num ();
#endif

               int j;

               if (picycle_thr[my_thr] == -1) {
                  plistloc_thr[my_thr].resize (2 * nimax + 1);
                  pimaskblk_thr[my_thr].resize (3 * _nblks_ini + 1);
                  int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
                  for (j = 0; j < _nblks_ini; j++)
                     pimaskblk[j] = -1;
               }
               int icycleblk = picycle_thr[my_thr];
               _Int *plistloc = plistloc_thr[my_thr].Ptr ();
               _Int *pialoc = plistloc + nimax;
               int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
               icycleblk++;
               int niloc = (int) (_blks_ini[ipar + 1] - _blks_ini[ipar]);
               int ibs = (int) (_blks_ini[ipar] - _blks_ini[iblk0]);
               int ishift = (int) _ia[ibeg + ibs];
               for (j = 0; j < niloc; j++) {
                  plistloc[j] = (_Int) (j + _blks_ini[ipar]);
               }
               for (j = 0; j <= niloc; j++) {
                  pialoc[j] = _ia[ibeg + ibs + j] - _ia[ibeg + ibs];
               }
               CBMatrix < _Int, _Flt > hblk_loc (ipar, niloc, plistloc, pialoc,
                                                 _ja + ishift, _a + ishift, _nblks_ini,
                                                 _blks_ini, icycleblk, pimaskblk);
               picycle_thr[my_thr] = icycleblk;
               phblk_arr[ipar - iblk0].ReplaceFree (hblk_loc);
            }

// Combine array blocks into one hblock

            int nzja_tot = 0;

            int j;

            for (j = 0; j < niblk; j++)
               nzja_tot += phblk_arr[j].GetNzblk ();

            CBMatrix < _Int, _Flt > hblk;

            hblk.ResizeASub (nzja_tot);
            hblk.SetNzblk (nzja_tot);

            CMatrix < _Int, _Flt > *pASub = hblk.GetASubArr ();
            CMatrix < int, float >*phmatr_str = hblk.GetHMatrStr ();

            phmatr_str->ResizeAndSetAllSp (niblk, niblk, nzja_tot, nzja_tot);

            int *plist_str = phmatr_str->GetListArr ();
            int *plist2_str = phmatr_str->GetList2Arr ();
            int *pia_str = phmatr_str->GetIaArr ();
            int *pja_str = phmatr_str->GetJaArr ();
            int *pja2_str = phmatr_str->GetJa2Arr ();

            for (j = 0; j < niblk; j++)
               plist_str[j] = j;
            for (j = 0; j < niblk; j++)
               plist2_str[j] = i;

            nzja_tot = 0;

            pia_str[0] = 0;

            int k, kk, kk2;

            for (j = 0; j < niblk; j++) {
               CMatrix < _Int, _Flt > *pASub_temp = phblk_arr[j].GetASubArr ();
               CMatrix < int, float >*phmatr_str_temp = phblk_arr[j].GetHMatrStr ();
               int nzja_temp = phmatr_str_temp->GetNzja ();
               int *pja_temp = phmatr_str_temp->GetJaArr ();
               for (k = 0; k < nzja_temp; k++) {
                  kk = pja_temp[k];
                  kk2 = pblk2hblks_ini[kk];
                  pja_str[nzja_tot] = kk - _hblk2blks_ini[kk2];
                  pja2_str[nzja_tot] = kk2;
                  pASub[nzja_tot].ReplaceFree (pASub_temp[k]);
                  nzja_tot++;
               }
               pia_str[j + 1] = nzja_tot;
            }
            phmatr_ini_arr[i].ReplaceFree (hblk);
            ibeg += (int) (_hblks[i + 1] - _hblks[i]);
         }
      }

// Compute new local partitioning and ordering if necessary

      if (_params->b_new_part) {

         double time0;

         CMPIDataExchange::Synchronize (pcommloc);

         time0 = CMPIDataExchange::GetWallTimeMPI ();

         this->order_ini.resize (_nhblks + 1);

         vector < int >nblks_temp (_nhblks + 1);
         vector < vector < int > >blks_temp (_nhblks + 1);

         CVectorData < int >*porder_ini = &this->order_ini[0];

         int *pnblks_temp = &nblks_temp[0];
         vector < int >*pblks_temp = &blks_temp[0];

         for (i = 0; i < _nhblks; i++) {
            if (_hblk2cpu[i] == myid) {

               int niblk = _hblk2blks_ini[i + 1] - _hblk2blks_ini[i];
               int iblkbeg = _hblk2blks_ini[i];
               int iblkend = _hblk2blks_ini[i + 1] - 1;

               CMatrix < _Int, _Flt > *pASub_curr = phmatr_ini_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str_curr = phmatr_ini_arr[i].GetHMatrStr ();

               int *pia_str = phmatr_str_curr->GetIaArr ();
               int *pja_str = phmatr_str_curr->GetJaArr ();
               int *pja2_str = phmatr_str_curr->GetJa2Arr ();

               if (true) {

// Create block diagonal submatrix

                  CBMatrix < _Int, _Flt > hmatr_dia;

                  int nzblk = 0;

                  {

                     int ilist, j, jhblk;

                     for (ilist = 0; ilist < niblk; ilist++) {
                        for (j = pia_str[ilist]; j < pia_str[ilist + 1]; j++) {
                           jhblk = pja2_str[j];
                           if (jhblk == i)
                              nzblk++;
                        }
                     }

                  }

                  hmatr_dia.SetNzblk (nzblk);
                  hmatr_dia.ResizeASub (nzblk);

                  CMatrix < _Int, _Flt > *pASub_dia = hmatr_dia.GetASubArr ();
                  CMatrix < int, float >*phmatr_str_dia = hmatr_dia.GetHMatrStr ();

                  phmatr_str_dia->ResizeAndSetAllSp (niblk, niblk, nzblk, nzblk);

                  int *plist_dia = phmatr_str_dia->GetListArr ();
                  int *plist2_dia = phmatr_str_dia->GetList2Arr ();
                  int *pia_dia = phmatr_str_dia->GetIaArr ();
                  int *pja_dia = phmatr_str_dia->GetJaArr ();
                  int *pja2_dia = phmatr_str_dia->GetJa2Arr ();

                  {

                     int j;

                     for (j = 0; j < niblk; j++)
                        plist_dia[j] = j;
                     for (j = 0; j < niblk; j++)
                        plist2_dia[j] = 0;
                     for (j = 0; j < nzblk; j++)
                        pja2_dia[j] = 0;

                     nzblk = 0;

                     int ilist, jhblk;

                     pia_dia[0] = 0;

                     for (ilist = 0; ilist < niblk; ilist++) {
                        for (j = pia_str[ilist]; j < pia_str[ilist + 1]; j++) {
                           jhblk = pja2_str[j];
                           if (jhblk == i) {
                              pja_dia[nzblk] = pja_str[j];
                              nzblk++;
                           }
                        }
                        pia_dia[ilist + 1] = nzblk;
                     }

                  }

                  {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                     for (int ipar = 0; ipar < niblk; ipar++) {
                        int ibs = pia_dia[ipar];
                        int j, jhblk;
                        for (j = pia_str[ipar]; j < pia_str[ipar + 1]; j++) {
                           jhblk = pja2_str[j];
                           if (jhblk == i) {
                              pASub_dia[ibs] = pASub_curr[j];
                              ibs++;
                           }
                        }
                     }

                  }

// Symmetrize block sparsity

                  int nhblks_1cpu = 1;
                  int hblk2cpu_1cpu = 0;
                  int hblk2blks_1cpu[2];
                  hblk2blks_1cpu[0] = 0;
                  hblk2blks_1cpu[1] = niblk;

                  vector < int >blk2hblks_1cpu (niblk + 1);
                  vector < long long >blks_1cpu (niblk + 1);

                  int *pblk2hblks_1cpu = &blk2hblks_1cpu[0];
                  long long *pblks_1cpu = &blks_1cpu[0];

                  {
                     int j;
                     for (j = 0; j < niblk; j++)
                        pblk2hblks_1cpu[j] = 0;
                     for (j = 0; j <= niblk; j++)
                        pblks_1cpu[j] = _blks_ini[iblkbeg + j] - _blks_ini[iblkbeg];
                  }

                  CBMatrix < _Int, _Flt > hmatr_dia_symm;

                  CBMatrix < _Int, _Flt >::SymmetrizeSubmatrices (NULL, nhblks_1cpu,
                                                                  &hblk2cpu_1cpu,
                                                                  hblk2blks_1cpu,
                                                                  pblk2hblks_1cpu,
                                                                  pblks_1cpu, &hmatr_dia,
                                                                  &hmatr_dia_symm);

// Perform decomposition

                  int degree_decomp = 20;
                  int isize_max = 1000;

                  int nparts_new;

                  int nloc_dia = (int) pblks_1cpu[niblk];

                  porder_ini[i].resize (nloc_dia);
                  int *pporder_ini = porder_ini[i].Ptr ();

                  vector < long long >parts_new;

//                  ffout << " Ihblk = " << i << endl;
//                  ffout << " Degree = " << degree_decomp << " isize_max = " << isize_max << " nparts = " << _params->nparts << endl;
//                  hmatr_dia_symm.PrintHMatrix (ffout);

                  CBMatrix < _Int, _Flt >::DecompWeights_thr (false, degree_decomp,
                                                              isize_max, _params->nparts,
                                                              niblk, pblks_1cpu,
                                                              hmatr_dia_symm, nparts_new,
                                                              parts_new, pporder_ini);

                  pnblks_temp[i] = nparts_new;

                  pblks_temp[i].resize (nparts_new + 1);

                  int *ppblks_temp = &pblks_temp[i][0];

                  {
                     int j;
                     long long *pparts_new = &parts_new[0];
                     for (j = 0; j <= nparts_new; j++)
                        ppblks_temp[j] = (int) pparts_new[j];
                  }

//                  ffout << " nparts_new = " << nparts_new << endl;
//                  PrintArray (ffout, " ppblks_temp", nparts_new+1, ppblks_temp);
//                  PrintArray (ffout, " pporder_ini", nloc_dia, pporder_ini);

               } else {

// For each hblock combine its diagonal part into one matrix

                  vector < int >ibs_blks (niblk + 1);
                  int *pibs_blks = &ibs_blks[0];

                  pibs_blks[0] = 0;

                  {

                     int j, k, jj2;

                     for (j = 0; j < niblk; j++)
                        pibs_blks[j + 1] = 0;

                     for (j = 0; j < niblk; j++) {
                        for (k = pia_str[j]; k < pia_str[j + 1]; k++) {
                           jj2 = pja2_str[k];
                           if (jj2 == i) {
                              pibs_blks[j + 1] += pASub_curr[k].GetNzja ();
                           }
                        }
                     }

                     for (j = 0; j < niblk; j++)
                        pibs_blks[j + 1] = pibs_blks[j] + pibs_blks[j + 1];

                  }

                  int nloc_dia = (int) (_blks_ini[iblkend + 1] - _blks_ini[iblkbeg]);
                  int nzja_dia = pibs_blks[niblk];

                  CMatrix < _Int, _Flt > A_sub;

                  A_sub.ResizeAndSetAllSp (nloc_dia, 0, nzja_dia, 0);

                  _Int *plist_sub = A_sub.GetListArr ();
                  _Int *pia_sub = A_sub.GetIaArr ();
                  _Int *pja_sub = A_sub.GetJaArr ();

                  pia_sub[0] = 0;

                  {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                     for (int ipar = 0; ipar < niblk; ipar++) {
                        int nlist_curr =
                           (int) (_blks_ini[iblkbeg + ipar + 1] -
                                  _blks_ini[iblkbeg + ipar]);
                        int ibeg = (int) (_blks_ini[iblkbeg + ipar] - _blks_ini[iblkbeg]);
                        int j;
                        for (j = 0; j < nlist_curr; j++)
                           plist_sub[ibeg + j] = ibeg + j;
                        CVectorData < int >ia_temp (nlist_curr + 1);
                        int *pia_temp = ia_temp.Ptr ();
                        for (j = 0; j <= nlist_curr; j++)
                           pia_temp[j] = 0;
                        int jj, jj2, k, kk;
                        for (j = pia_str[ipar]; j < pia_str[ipar + 1]; j++) {
                           jj = pja_str[j];
                           jj2 = pja2_str[j];
                           if (jj2 == i) {
                              int nlist_temp = pASub_curr[j].GetNlist ();
                              _Int *plist_temp = pASub_curr[j].GetListArr ();
                              _Int *pia_temp1 = pASub_curr[j].GetIaArr ();
                              for (k = 0; k < nlist_temp; k++) {
                                 kk = (int) plist_temp[k];
                                 pia_temp[kk + 1] +=
                                    (int) (pia_temp1[k + 1] - pia_temp1[k]);
                              }
                           }
                        }
                        for (j = 0; j < nlist_curr; j++)
                           pia_temp[j + 1] = pia_temp[j] + pia_temp[j + 1];
                        CVectorData < int >iptr (nlist_curr + 1);
                        int *piptr = iptr.Ptr ();
                        for (j = 0; j < nlist_curr; j++)
                           piptr[j] = pia_temp[j];
                        int ibs0 = pibs_blks[ipar];
                        int jshift, ind, kj, kkj;
                        for (j = pia_str[ipar]; j < pia_str[ipar + 1]; j++) {
                           jj = pja_str[j];
                           jj2 = pja2_str[j];
                           if (jj2 == i) {
                              jshift =
                                 (int) (_blks_ini[iblkbeg + jj] - _blks_ini[iblkbeg]);
                              int nlist_temp = pASub_curr[j].GetNlist ();
                              _Int *plist_temp = pASub_curr[j].GetListArr ();
                              _Int *pia_temp1 = pASub_curr[j].GetIaArr ();
                              _Int *pja_temp = pASub_curr[j].GetJaArr ();
                              for (k = 0; k < nlist_temp; k++) {
                                 kk = (int) plist_temp[k];
                                 ind = piptr[kk];
                                 for (kj = (int) pia_temp1[k]; kj < pia_temp1[k + 1];
                                      kj++) {
                                    kkj = (int) pja_temp[kj];
                                    pja_sub[ibs0 + ind] = jshift + kkj;
                                    ind++;
                                 }
                                 piptr[kk] = ind;
                              }
                           }
                        }
                        for (j = 0; j < nlist_curr; j++)
                           pia_sub[ibeg + j + 1] = ibs0 + pia_temp[j + 1];
                     }
                  }

// Symmetrize matrix

                  vector < _Int > *pia_dia = A_sub.GetIa ();
                  vector < _Int > *pja_dia = A_sub.GetJa ();

                  CMatrix < _Int, _Flt > A_sub_symm;

                  vector < _Int > *pia_dia_symm = A_sub_symm.GetIa ();
                  vector < _Int > *pja_dia_symm = A_sub_symm.GetJa ();

                  CFct_impl < _Int, _Flt >::SymmetrizeSparsity (nloc_dia, *pia_dia,
                                                                *pja_dia, *pia_dia_symm,
                                                                *pja_dia_symm);

                  int nzja_symm = (int) ((*pia_dia_symm)[nloc_dia]);

                  A_sub_symm.SetNlist (nloc_dia);
                  A_sub_symm.SetNzja (nzja_symm);

                  A_sub_symm.ResizeList (nloc_dia);

                  _Int *plist_sub_symm = A_sub_symm.GetListArr ();
                  _Int *pia_sub_symm_arr = A_sub_symm.GetIaArr ();
                  _Int *pja_sub_symm_arr = A_sub_symm.GetJaArr ();

                  int j;

                  for (j = 0; j < nloc_dia; j++)
                     plist_sub_symm[j] = (_Int) j;

// Compute new ordering and partitioning

                  CMatrix < _Int, _Flt > amatr_strW;

                  amatr_strW.ResizeAndSetAllSp (nloc_dia, nloc_dia, nzja_symm, nzja_symm);

                  _Int *plist_W = amatr_strW.GetListArr ();
                  _Int *plist2_W = amatr_strW.GetList2Arr ();
                  _Int *pia_W = amatr_strW.GetIaArr ();
                  _Int *pja_W = amatr_strW.GetJaArr ();
                  _Int *pja2_W = amatr_strW.GetJa2Arr ();

                  for (j = 0; j < nloc_dia; j++)
                     plist_W[j] = (_Int) j;
                  for (j = 0; j < nloc_dia; j++)
                     plist2_W[j] = 1;
                  for (j = 0; j <= nloc_dia; j++)
                     pia_W[j] = pia_sub_symm_arr[j];
                  for (j = 0; j < nzja_symm; j++)
                     pja_W[j] = pja_sub_symm_arr[j];
                  for (j = 0; j < nzja_symm; j++)
                     pja2_W[j] = 1;

                  CVectorData < int >partition (nloc_dia);
                  int *ppartition = partition.Ptr ();

                  CMatrix < _Int, _Flt >::DecompWeights (false, amatr_strW,
                                                         _params->nparts, ppartition);

// Compute ordering and partitioning

                  porder_ini[i].resize (nloc_dia);
                  int *pporder_ini = porder_ini[i].Ptr ();

                  pblks_temp[i].resize (_params->nparts + 1);

                  int *ppblks_temp = &pblks_temp[i][0];

                  for (j = 0; j <= _params->nparts; j++)
                     ppblks_temp[j] = 0;

                  int jj;

                  for (j = 0; j < nloc_dia; j++) {
                     jj = ppartition[j];
                     ppblks_temp[jj + 1]++;
                  }

                  for (j = 0; j < _params->nparts; j++)
                     ppblks_temp[j + 1] = ppblks_temp[j] + ppblks_temp[j + 1];

                  vector < int >iptr (_params->nparts + 1);
                  int *piptr = &iptr[0];

                  for (j = 0; j < _params->nparts; j++)
                     piptr[j] = ppblks_temp[j];

                  int k;

                  for (j = 0; j < nloc_dia; j++) {
                     jj = ppartition[j];
                     k = piptr[jj];
                     pporder_ini[j] = k;
                     piptr[jj]++;
                  }

                  int nparts_new = 0;

                  for (j = 0; j < _params->nparts; j++) {
                     if (ppblks_temp[j + 1] > ppblks_temp[j]) {
                        ppblks_temp[nparts_new + 1] = ppblks_temp[j + 1];
                        nparts_new++;
                     }
                  }

                  pnblks_temp[i] = nparts_new;

               }
            }
         }

// Collect new global blocks partitioning

         CMPIDataExchange::ExchangeArray (_pcomm, 'I', '+', _nhblks, pnblks_temp);

         this->hblk2blks.resize (_nhblks + 1);

         int *phblk2blks = &this->hblk2blks[0];

         phblk2blks[0] = 0;

         for (i = 0; i < _nhblks; i++)
            phblk2blks[i + 1] = phblk2blks[i] + pnblks_temp[i];

         int nblks_new = phblk2blks[_nhblks];

         this->nblks = nblks_new;

         this->blks.resize (nblks_new + 1);

         long long *pblks = &this->blks[0];

         for (i = 0; i <= nblks_new; i++)
            pblks[i] = 0;

         int j;

         for (i = 0; i < _nhblks; i++) {
            if (_hblk2cpu[i] == myid) {
               int ibeg = phblk2blks[i];
               int *ppblks_temp = &pblks_temp[i][0];
               for (j = 0; j < pnblks_temp[i]; j++) {
                  pblks[ibeg + j + 1] = ppblks_temp[j + 1] - ppblks_temp[j];
               }
            }
         }

         CMPIDataExchange::ExchangeArray (_pcomm, 'L', '+', nblks_new + 1, pblks);

         for (i = 0; i < nblks_new; i++)
            pblks[i + 1] = pblks[i] + pblks[i + 1];

         this->blk2hblks.resize (nblks_new + 1);

         int *pblk2hblks = &this->blk2hblks[0];

         for (i = 0; i < _nhblks; i++) {
            for (j = phblk2blks[i]; j < phblk2blks[i + 1]; j++) {
               pblk2hblks[j] = i;
            }
         }

         if (_params->b_fast_transform) {
            this->PrepareFastTransform (phmatr_ini_arr);
         }

         double time1;

         CMPIDataExchange::Synchronize (pcommloc);

         time1 = CMPIDataExchange::GetWallTimeMPI ();

         _stats->dtime_part = time1 - time0;

      }
// Compute reordered matrix

      int nblks_loc = this->nblks;

      int *phblk2blks = &this->hblk2blks[0];
      int *pblk2hblks = &this->blk2hblks[0];
      long long *pblks = &this->blks[0];

      CVectorData < int >*porder_ini = &this->order_ini[0];

      this->hmatr_arr.resize (_nhblks + 1);

      CBMatrix < _Int, _Flt > *phmatr_arr = &this->hmatr_arr[0];

      CMPIDataExchange::Synchronize (pcommloc);

      double time0;

      time0 = CMPIDataExchange::GetWallTimeMPI ();

      if (this->b_fast_transform) {

         this->UseFastTransform (phmatr_ini_arr);

      } else {

         CBMatrix < _Int, _Flt >::ReorderHMatrixDiag (_pcomm, _nhblks, _hblk2cpu,
                                                      _hblk2blks_ini, _nblks_ini,
                                                      _blks_ini, phmatr_ini_arr,
                                                      porder_ini, phblk2blks, nblks_loc,
                                                      pblks, phmatr_arr);

      }

      phmatr_arr = &this->hmatr_arr[0];

      CMPIDataExchange::Synchronize (pcommloc);

      double time1;

      time1 = CMPIDataExchange::GetWallTimeMPI ();

      _stats->dtime_ord = time1 - time0;

// Symmetrize hmatr

      if (!this->b_fast_transform) {

         vector < CBMatrix < _Int, _Flt > >hmatr_symm_arr (_nhblks + 1);
         CBMatrix < _Int, _Flt > *phmatr_symm_arr = &hmatr_symm_arr[0];

         CBMatrix < _Int, _Flt >::SymmetrizeSubmatrices (_pcomm, _nhblks, _hblk2cpu,
                                                         phblk2blks, pblk2hblks, pblks,
                                                         phmatr_arr, phmatr_symm_arr);

         for (i = 0; i < _nhblks; i++) {
            if (_hblk2cpu[i] == myid) {
               phmatr_arr[i].ReplaceFree (phmatr_symm_arr[i]);
            }
         }

      }
// Init MvmA structures

      this->mvm.InitControl (_pcomm, _nhblks, _hblk2cpu, phblk2blks, pblk2hblks, _hblks,
                             nblks_loc, pblks);
      this->mvm.InitMvmA (phmatr_arr);

      this->b_use_ini = true;

   }

// Prepare matrix data only including preliminary repartitioning of the matrix taking into account wells data
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::PrepareMatrixWells (void *_pcomm, SParams * _params, SStatData * _stats,
                                      int _nhblks_ini, int *_hblk2cpu_ini,
                                      int *_hblk2blks_ini, long long *_hblks_ini,
                                      int _nblks_ini, long long *_blks_ini,
                                      vector < int >*_p_ia_wells_ext,
                                      vector < int >*_p_ja_wells_ext, _Int * _ia,
                                      _Int * _ja, _Flt * _a)
   {

      _stats->n_well = 0;

      _stats->dtime_part = 0.0e0;
      _stats->dtime_ord = 0.0e0;

// Store initial wells control data if necessary

      void *pcommloc = _pcomm;

      if (_params->b_new_part) {

         this->pcomm = _pcomm;
         this->nhblks_ini = _nhblks_ini;
         this->nblks_ini = _nblks_ini;

         this->hblks_ini.resize (_nhblks_ini + 1);
         this->hblk2cpu_ini.resize (_nhblks_ini + 1);
         this->hblk2blks_ini.resize (_nhblks_ini + 1);
         this->blks_ini.resize (_nblks_ini + 1);

         long long *phblks_ini = &this->hblks_ini[0];
         int *phblk2cpu_ini = &this->hblk2cpu_ini[0];
         int *phblk2blks_ini = &this->hblk2blks_ini[0];
         long long *pblks_ini = &this->blks_ini[0];

         int i;

         for (i = 0; i <= _nhblks_ini; i++)
            phblks_ini[i] = _hblks_ini[i];
         for (i = 0; i < _nhblks_ini; i++)
            phblk2cpu_ini[i] = _hblk2cpu_ini[i];
         for (i = 0; i <= _nhblks_ini; i++)
            phblk2blks_ini[i] = _hblk2blks_ini[i];

         for (i = 0; i <= _nblks_ini; i++)
            pblks_ini[i] = _blks_ini[i];

      }
// Create matrix as set of hblocks

      int myid = CMPIDataExchange::GetMyid (_pcomm);

      vector < CBMatrix < _Int, _Flt > >hmatr_ini_arr (_nhblks_ini + 1);
      CBMatrix < _Int, _Flt > *phmatr_ini_arr = &hmatr_ini_arr[0];

      vector < int >blk2hblks_ini (_nblks_ini);
      int *pblk2hblks_ini = &blk2hblks_ini[0];

      {

         int i, j;

         for (i = 0; i < _nhblks_ini; i++) {
            for (j = _hblk2blks_ini[i]; j < _hblk2blks_ini[i + 1]; j++) {
               pblk2hblks_ini[j] = i;
            }
         }

      }

      {

         int nimax = 0;

         int i, niloc;

         for (i = 0; i < _nblks_ini; i++) {
            niloc = (int) (_blks_ini[i + 1] - _blks_ini[i]);
            if (niloc > nimax)
               nimax = niloc;
         }

         int n_thr = 1;

#ifdef USE_THREADS
         n_thr = omp_get_max_threads ();
#endif

         vector < int >icycle_thr (n_thr + 1);
         vector < CVectorData < _Int > >listloc_thr (n_thr + 1);
         vector < CVectorData < int > >imaskblk_thr (n_thr + 1);

         int *picycle_thr = &icycle_thr[0];
         CVectorData < _Int > *plistloc_thr = &listloc_thr[0];
         CVectorData < int >*pimaskblk_thr = &imaskblk_thr[0];

         for (i = 0; i < n_thr; i++)
            picycle_thr[i] = -1;

         int ibeg = 0;

         for (i = 0; i < _nhblks_ini; i++) {
            if (_hblk2cpu_ini[i] == myid) {

               int niblk = _hblk2blks_ini[i + 1] - _hblk2blks_ini[i];
               vector < CBMatrix < _Int, _Flt > >hblk_arr (niblk + 1);
               CBMatrix < _Int, _Flt > *phblk_arr = &hblk_arr[0];
               int iblk0 = _hblk2blks_ini[i];

               {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = _hblk2blks_ini[i]; ipar < _hblk2blks_ini[i + 1]; ipar++) {

                     int my_thr = 0;
#ifdef USE_THREADS
                     my_thr = omp_get_thread_num ();
#endif

                     int j;

                     if (picycle_thr[my_thr] == -1) {
                        plistloc_thr[my_thr].resize (2 * nimax + 1);
                        pimaskblk_thr[my_thr].resize (3 * _nblks_ini + 1);
                        int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
                        for (j = 0; j < _nblks_ini; j++)
                           pimaskblk[j] = -1;
                     }
                     int icycleblk = picycle_thr[my_thr];
                     _Int *plistloc = plistloc_thr[my_thr].Ptr ();
                     _Int *pialoc = plistloc + nimax;
                     int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
                     icycleblk++;
                     int niloc = (int) (_blks_ini[ipar + 1] - _blks_ini[ipar]);
                     int ibs = (int) (_blks_ini[ipar] - _blks_ini[iblk0]);
                     int ishift = (int) _ia[ibeg + ibs];
                     for (j = 0; j < niloc; j++) {
                        plistloc[j] = (_Int) (j + _blks_ini[ipar]);
                     }
                     for (j = 0; j <= niloc; j++) {
                        pialoc[j] = _ia[ibeg + ibs + j] - _ia[ibeg + ibs];
                     }
                     CBMatrix < _Int, _Flt > hblk_loc (ipar, niloc, plistloc, pialoc,
                                                       _ja + ishift, _a + ishift,
                                                       _nblks_ini, _blks_ini, icycleblk,
                                                       pimaskblk);
                     picycle_thr[my_thr] = icycleblk;
                     phblk_arr[ipar - iblk0].ReplaceFree (hblk_loc);
                  }
               }

// Combine array blocks into one hblock

               int nzja_tot = 0;

               int j;

               for (j = 0; j < niblk; j++)
                  nzja_tot += phblk_arr[j].GetNzblk ();

               CBMatrix < _Int, _Flt > hblk;

               hblk.ResizeASub (nzja_tot);
               hblk.SetNzblk (nzja_tot);

               CMatrix < _Int, _Flt > *pASub = hblk.GetASubArr ();
               CMatrix < int, float >*phmatr_str = hblk.GetHMatrStr ();

               phmatr_str->ResizeAndSetAllSp (niblk, niblk, nzja_tot, nzja_tot);

               int *plist_str = phmatr_str->GetListArr ();
               int *plist2_str = phmatr_str->GetList2Arr ();
               int *pia_str = phmatr_str->GetIaArr ();
               int *pja_str = phmatr_str->GetJaArr ();
               int *pja2_str = phmatr_str->GetJa2Arr ();

               for (j = 0; j < niblk; j++)
                  plist_str[j] = j;
               for (j = 0; j < niblk; j++)
                  plist2_str[j] = i;

               nzja_tot = 0;

               pia_str[0] = 0;

               int k, kk, kk2;

               for (j = 0; j < niblk; j++) {
                  CMatrix < _Int, _Flt > *pASub_temp = phblk_arr[j].GetASubArr ();
                  CMatrix < int, float >*phmatr_str_temp = phblk_arr[j].GetHMatrStr ();
                  int nzja_temp = phmatr_str_temp->GetNzja ();
                  int *pja_temp = phmatr_str_temp->GetJaArr ();
                  for (k = 0; k < nzja_temp; k++) {
                     kk = pja_temp[k];
                     kk2 = pblk2hblks_ini[kk];
                     pja_str[nzja_tot] = kk - _hblk2blks_ini[kk2];
                     pja2_str[nzja_tot] = kk2;
                     pASub[nzja_tot].ReplaceFree (pASub_temp[k]);
                     nzja_tot++;
                  }
                  pia_str[j + 1] = nzja_tot;
               }
               phmatr_ini_arr[i].ReplaceFree (hblk);
               ibeg += (int) (_hblks_ini[i + 1] - _hblks_ini[i]);

            }
         }
      }

// Compute new local partitioning and ordering if necessary

      if (_params->b_new_part) {

         double time0;

         CMPIDataExchange::Synchronize (pcommloc);

         time0 = CMPIDataExchange::GetWallTimeMPI ();

// Symmetrize matrix data

         {

            vector < CBMatrix < _Int, _Flt > >hmatr_symm_arr (_nhblks_ini + 1);
            CBMatrix < _Int, _Flt > *phmatr_symm_arr = &hmatr_symm_arr[0];

            CBMatrix < _Int, _Flt >::SymmetrizeSubmatrices (_pcomm, _nhblks_ini,
                                                            _hblk2cpu_ini, _hblk2blks_ini,
                                                            pblk2hblks_ini, _blks_ini,
                                                            phmatr_ini_arr,
                                                            phmatr_symm_arr);

            for (int i = 0; i < _nhblks_ini; i++) {
               if (_hblk2cpu_ini[i] == myid) {
                  phmatr_ini_arr[i].ReplaceFree (phmatr_symm_arr[i]);
               }
            }

         }

// Perform decomposition and store decomp data

         {

            int nz_ord = 0;

            int ihblk;

            vector < long long >ibs_ord_wells (_nhblks_ini + 1);
            long long *pibs_ord_wells = &ibs_ord_wells[0];

            for (ihblk = 0; ihblk < _nhblks_ini; ihblk++)
               pibs_ord_wells[ihblk] = -1;

            for (ihblk = 0; ihblk < _nhblks_ini; ihblk++) {
               if (_hblk2cpu_ini[ihblk] == myid) {
                  int ibegblk = _hblk2blks_ini[ihblk];
                  int iendblk = _hblk2blks_ini[ihblk + 1] - 1;
                  int ni_hblk = (int) (_blks_ini[iendblk + 1] - _blks_ini[ibegblk]);
                  pibs_ord_wells[ihblk] = nz_ord;
                  nz_ord += ni_hblk;
               }
            }

            this->order_wells.resize (nz_ord);
            this->order2ind_wells.resize (2 * nz_ord);

            long long *porder_wells = this->order_wells.Ptr ();
            int *porder2ind_wells = this->order2ind_wells.Ptr ();

            int n_wells_new;
            int nhblks_new;
            vector < int >hblk2cpu_new;
            vector < int >hblk2blk_new;
            vector < int >blk2type_new;
            vector < long long >blks_new;

            CBMatrix < _Int, _Flt >::DecompWells (_pcomm, _params->nneib_well,
                                                  _params->eps_dia_well,
                                                  _params->thresh_max_well,
                                                  _params->ncycle_well, _params->degree,
                                                  _params->isize_max, _params->isize_max2,
                                                  _params->nparts_bdecomp,
                                                  _params->nparts, _params->nparts_W,
                                                  _nhblks_ini, _hblk2cpu_ini,
                                                  _hblk2blks_ini, _blks_ini,
                                                  _p_ia_wells_ext, _p_ja_wells_ext,
                                                  phmatr_ini_arr, n_wells_new, nhblks_new,
                                                  hblk2cpu_new, hblk2blk_new,
                                                  blk2type_new, blks_new, porder_wells);

            _stats->n_well = n_wells_new;

            this->nhblks = nhblks_new;
            this->hblk2cpu.swap (hblk2cpu_new);
            this->hblk2blks.swap (hblk2blk_new);
            this->blks.swap (blks_new);

            int *phblk2blks = &this->hblk2blks[0];
            int *pblk2type_new = &blk2type_new[0];
            long long *pblks = &this->blks[0];

            int nblks_new = phblk2blks[nhblks_new];

            this->nblks = nblks_new;

            this->hblks.resize (nhblks_new + 1);
            this->blk2hblks.resize (nblks_new + 1);

            long long *phblks = &this->hblks[0];
            int *pblk2hblks = &this->blk2hblks[0];

            int j;

            phblks[0] = 0;
            for (ihblk = 0; ihblk < nhblks_new; ihblk++) {
               int ibegblk = phblk2blks[ihblk];
               int iendblk = phblk2blks[ihblk + 1] - 1;
               int ni_hblk = (int) (pblks[iendblk + 1] - pblks[ibegblk]);
               phblks[ihblk + 1] = phblks[ihblk] + ni_hblk;
               for (j = ibegblk; j <= iendblk; j++)
                  pblk2hblks[j] = ihblk;
            }

// Compute 2index data assuming that order is in porder_wells

            {

               int n_thr = 1;

#ifdef USE_THREADS
               n_thr = omp_get_max_threads ();
#endif

               for (ihblk = 0; ihblk < _nhblks_ini; ihblk++) {
                  if (_hblk2cpu_ini[ihblk] == myid) {
                     int ibegblk = _hblk2blks_ini[ihblk];
                     int iendblk = _hblk2blks_ini[ihblk + 1] - 1;
                     int niblk = iendblk + 1 - ibegblk;
                     {
                        int nimax = 0;
                        {
                           int i, j, niloc;
                           for (i = 0; i < niblk; i++) {
                              j = ibegblk + i;
                              niloc = (int) (_blks_ini[j + 1] - _blks_ini[j]);
                              if (niloc > nimax)
                                 nimax = niloc;
                           }
                        }
                        CVectorData < CSortInt64 > iiarr_threads (nimax * n_thr);
                        CSortInt64 *piiarr_threads = iiarr_threads.Ptr ();
                        CVectorData < long long >listarr_threads (2 * nimax * n_thr);
                        long long *plistarr_threads = listarr_threads.Ptr ();

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                        for (int ipar = 0; ipar < niblk; ipar++) {
                           int my_thr = 0;
#ifdef USE_THREADS
                           my_thr = omp_get_thread_num ();
#endif

                           int jblk = ibegblk + ipar;
                           int niloc = (int) (_blks_ini[jblk + 1] - _blks_ini[jblk]);
                           int ibs =
                              (int) pibs_ord_wells[ihblk] + (int) (_blks_ini[jblk] -
                                                                   _blks_ini[ibegblk]);

                           CSortInt64 *piiarr_th = piiarr_threads + nimax * my_thr;
                           long long *plistarr_th = plistarr_threads + 2 * nimax * my_thr;
                           long long *plist2arr_th = plistarr_th + nimax;

                           int i;

                           for (i = 0; i < niloc; i++) {
                              piiarr_th[i].ival = porder_wells[ibs + i];
                              piiarr_th[i].i2val = ibs + i;
                           }
                           sort (piiarr_th, piiarr_th + niloc);

                           for (i = 0; i < niloc; i++) {
                              plistarr_th[i] = piiarr_th[i].ival;
                           }

                           long long ia_temp[2];

                           ia_temp[0] = 0;
                           ia_temp[1] = niloc;

                           CBMatrix < long long, float >::ComputeJa2 (nhblks_new, phblks,
                                                                      1, ia_temp,
                                                                      plistarr_th,
                                                                      plist2arr_th);

                           int jhblk_temp, jj_temp;
                           long long ind;

                           for (i = 0; i < niloc; i++) {
                              ind = piiarr_th[i].i2val;
                              jhblk_temp = (int) plist2arr_th[i];
                              jj_temp = (int) (plistarr_th[i] - phblks[jhblk_temp]);
                              porder2ind_wells[ind * 2] = jj_temp;
                              porder2ind_wells[ind * 2 + 1] = jhblk_temp;
                           }
                        }
                     }

                  }
               }

            }

            this->b_fast_transform = false;
            this->b_use_ini = false;
            this->b_use_wells = true;
            this->b_use_blksize = false;
            this->b_blk_wells = false;
            if (pblk2type_new[nblks_new - 1] == 0) {
               this->b_blk_wells = true;
            }

         }

// Prepare fast transform data if necessary

         if (_params->b_fast_transform) {
            this->PrepareFastTransform (phmatr_ini_arr);
         }

         double time1;

         CMPIDataExchange::Synchronize (pcommloc);

         time1 = CMPIDataExchange::GetWallTimeMPI ();

         _stats->dtime_part = time1 - time0;

      }
// Compute reordered matrix

      int nhblks_loc = this->nhblks;
      int nblks_loc = this->nblks;

      int *phblk2cpu = &this->hblk2cpu[0];
      long long *phblks = &this->hblks[0];
      int *phblk2blks = &this->hblk2blks[0];
      int *pblk2hblks = &this->blk2hblks[0];
      long long *pblks = &this->blks[0];

      long long *porder_wells = this->order_wells.Ptr ();

      this->hmatr_arr.resize (nhblks_loc + 1);

      CBMatrix < _Int, _Flt > *phmatr_arr = &this->hmatr_arr[0];

      CMPIDataExchange::Synchronize (pcommloc);

      double time0;

      time0 = CMPIDataExchange::GetWallTimeMPI ();

      if (this->b_fast_transform) {

         this->UseFastTransform (phmatr_ini_arr);

      } else {

         CBMatrix < _Int, _Flt >::ReorderHMatrix (_pcomm, _nhblks_ini, _hblk2cpu_ini,
                                                  _hblk2blks_ini, _blks_ini,
                                                  phmatr_ini_arr, porder_wells,
                                                  nhblks_loc, phblk2cpu, phblk2blks,
                                                  pblks, phmatr_arr);

      }

      phmatr_arr = &this->hmatr_arr[0];

      CMPIDataExchange::Synchronize (pcommloc);

      double time1;

      time1 = CMPIDataExchange::GetWallTimeMPI ();

      _stats->dtime_ord = time1 - time0;

// Init MvmA structures

      this->mvm.InitControl (_pcomm, nhblks_loc, phblk2cpu, phblk2blks, pblk2hblks,
                             phblks, nblks_loc, pblks);

      this->mvm.InitMvmA (phmatr_arr);

   }

// Prepare matrix data with blksize condensing
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::PrepareMatrixBlksize (void *_pcomm, SParams * _params,
                                        SStatData * _stats, int _nhblks_ini,
                                        int *_hblk2cpu_ini, int *_hblk2blks_ini,
                                        long long *_hblks_ini, int _nblks_ini,
                                        long long *_blks_ini, _Int * _ia, _Int * _ja,
                                        _Flt * _a)
   {

      if (_params->b_fast_transform) {
         cout <<
            " PrepareMatrixBlksize: Error: fast transform mode is not supported for blksize computations !!! "
            << endl;
         throw
            " PrepareMatrixBlksize: Error: fast transform mode is not supported for blksize computations !!! ";
      }
// Check block partitioning on entry

      int blksizeloc = _params->blksize_decomp;

      {
         int i, niloc;
         for (i = 0; i < _nblks_ini; i++) {
            niloc = (int) (_blks_ini[i + 1] - _blks_ini[i]);
            if (niloc % blksizeloc != 0) {
               cout <<
                  " PrepareMatrixBlksize: Error: incorrect block partitioning on entry !!! "
                  << endl;
               throw
                  " PrepareMatrixBlksize: Error: incorrect block partitioning on entry !!! ";
            }
         }
      }

// Fast return if necessary

      if (_params->b_new_part && _params->i_decomp_type == 3) {

         this->PrepareMatrixBase (_pcomm, _nhblks_ini, _hblk2cpu_ini, _hblk2blks_ini,
                                  _hblks_ini, _nblks_ini, _blks_ini, _ia, _ja, _a);

         this->b_use_blksize = true;
         this->blksize_bscl = blksizeloc;

         return;

      }
// Store initial wells control data if necessary

      _stats->dtime_part = 0.0e0;
      _stats->dtime_ord = 0.0e0;

      if (_params->b_new_part) {

         this->pcomm = _pcomm;
         this->nhblks_ini = _nhblks_ini;
         this->nblks_ini = _nblks_ini;

         this->hblks_ini.resize (_nhblks_ini + 1);
         this->hblk2cpu_ini.resize (_nhblks_ini + 1);
         this->hblk2blks_ini.resize (_nhblks_ini + 1);
         this->blks_ini.resize (_nblks_ini + 1);

         long long *phblks_ini = &this->hblks_ini[0];
         int *phblk2cpu_ini = &this->hblk2cpu_ini[0];
         int *phblk2blks_ini = &this->hblk2blks_ini[0];
         long long *pblks_ini = &this->blks_ini[0];

         int i;

         for (i = 0; i <= _nhblks_ini; i++)
            phblks_ini[i] = _hblks_ini[i];
         for (i = 0; i < _nhblks_ini; i++)
            phblk2cpu_ini[i] = _hblk2cpu_ini[i];
         for (i = 0; i <= _nhblks_ini; i++)
            phblk2blks_ini[i] = _hblk2blks_ini[i];

         for (i = 0; i <= _nblks_ini; i++)
            pblks_ini[i] = _blks_ini[i];

      }
// Create matrix as set of hblocks

      int myid = CMPIDataExchange::GetMyid (_pcomm);

      vector < CBMatrix < _Int, _Flt > >hmatr_ini_arr (_nhblks_ini + 1);
      CBMatrix < _Int, _Flt > *phmatr_ini_arr = &hmatr_ini_arr[0];

      vector < int >blk2hblks_ini (_nblks_ini);
      int *pblk2hblks_ini = &blk2hblks_ini[0];

      {

         int i, j;

         for (i = 0; i < _nhblks_ini; i++) {
            for (j = _hblk2blks_ini[i]; j < _hblk2blks_ini[i + 1]; j++) {
               pblk2hblks_ini[j] = i;
            }
         }

      }

      {

         int nimax = 0;

         int i, niloc;

         for (i = 0; i < _nblks_ini; i++) {
            niloc = (int) (_blks_ini[i + 1] - _blks_ini[i]);
            if (niloc > nimax)
               nimax = niloc;
         }

         int n_thr = 1;

#ifdef USE_THREADS
         n_thr = omp_get_max_threads ();
#endif

         vector < int >icycle_thr (n_thr + 1);
         vector < CVectorData < _Int > >listloc_thr (n_thr + 1);
         vector < CVectorData < int > >imaskblk_thr (n_thr + 1);

         int *picycle_thr = &icycle_thr[0];
         CVectorData < _Int > *plistloc_thr = &listloc_thr[0];
         CVectorData < int >*pimaskblk_thr = &imaskblk_thr[0];

         for (i = 0; i < n_thr; i++)
            picycle_thr[i] = -1;

         int ibeg = 0;

         for (i = 0; i < _nhblks_ini; i++) {
            if (_hblk2cpu_ini[i] == myid) {

               int niblk = _hblk2blks_ini[i + 1] - _hblk2blks_ini[i];
               vector < CBMatrix < _Int, _Flt > >hblk_arr (niblk + 1);
               CBMatrix < _Int, _Flt > *phblk_arr = &hblk_arr[0];
               int iblk0 = _hblk2blks_ini[i];

               {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = _hblk2blks_ini[i]; ipar < _hblk2blks_ini[i + 1]; ipar++) {

                     int my_thr = 0;
#ifdef USE_THREADS
                     my_thr = omp_get_thread_num ();
#endif

                     int j;

                     if (picycle_thr[my_thr] == -1) {
                        plistloc_thr[my_thr].resize (2 * nimax + 1);
                        pimaskblk_thr[my_thr].resize (3 * _nblks_ini + 1);
                        int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
                        for (j = 0; j < _nblks_ini; j++)
                           pimaskblk[j] = -1;
                     }
                     int icycleblk = picycle_thr[my_thr];
                     _Int *plistloc = plistloc_thr[my_thr].Ptr ();
                     _Int *pialoc = plistloc + nimax;
                     int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
                     icycleblk++;
                     int niloc = (int) (_blks_ini[ipar + 1] - _blks_ini[ipar]);
                     int ibs = (int) (_blks_ini[ipar] - _blks_ini[iblk0]);
                     int ishift = (int) _ia[ibeg + ibs];
                     for (j = 0; j < niloc; j++) {
                        plistloc[j] = (_Int) (j + _blks_ini[ipar]);
                     }
                     for (j = 0; j <= niloc; j++) {
                        pialoc[j] = _ia[ibeg + ibs + j] - _ia[ibeg + ibs];
                     }
                     CBMatrix < _Int, _Flt > hblk_loc (ipar, niloc, plistloc, pialoc,
                                                       _ja + ishift, _a + ishift,
                                                       _nblks_ini, _blks_ini, icycleblk,
                                                       pimaskblk);
                     picycle_thr[my_thr] = icycleblk;
                     phblk_arr[ipar - iblk0].ReplaceFree (hblk_loc);
                  }
               }

// Combine array blocks into one hblock

               int nzja_tot = 0;

               int j;

               for (j = 0; j < niblk; j++)
                  nzja_tot += phblk_arr[j].GetNzblk ();

               CBMatrix < _Int, _Flt > hblk;

               hblk.ResizeASub (nzja_tot);
               hblk.SetNzblk (nzja_tot);

               CMatrix < _Int, _Flt > *pASub = hblk.GetASubArr ();
               CMatrix < int, float >*phmatr_str = hblk.GetHMatrStr ();

               phmatr_str->ResizeAndSetAllSp (niblk, niblk, nzja_tot, nzja_tot);

               int *plist_str = phmatr_str->GetListArr ();
               int *plist2_str = phmatr_str->GetList2Arr ();
               int *pia_str = phmatr_str->GetIaArr ();
               int *pja_str = phmatr_str->GetJaArr ();
               int *pja2_str = phmatr_str->GetJa2Arr ();

               for (j = 0; j < niblk; j++)
                  plist_str[j] = j;
               for (j = 0; j < niblk; j++)
                  plist2_str[j] = i;

               nzja_tot = 0;

               pia_str[0] = 0;

               int k, kk, kk2;

               for (j = 0; j < niblk; j++) {
                  CMatrix < _Int, _Flt > *pASub_temp = phblk_arr[j].GetASubArr ();
                  CMatrix < int, float >*phmatr_str_temp = phblk_arr[j].GetHMatrStr ();
                  int nzja_temp = phmatr_str_temp->GetNzja ();
                  int *pja_temp = phmatr_str_temp->GetJaArr ();
                  for (k = 0; k < nzja_temp; k++) {
                     kk = pja_temp[k];
                     kk2 = pblk2hblks_ini[kk];
                     pja_str[nzja_tot] = kk - _hblk2blks_ini[kk2];
                     pja2_str[nzja_tot] = kk2;
                     pASub[nzja_tot].ReplaceFree (pASub_temp[k]);
                     nzja_tot++;
                  }
                  pia_str[j + 1] = nzja_tot;
               }
               phmatr_ini_arr[i].ReplaceFree (hblk);
               ibeg += (int) (_hblks_ini[i + 1] - _hblks_ini[i]);

            }
         }
      }

// Symmetrize matrix data

      {

         vector < CBMatrix < _Int, _Flt > >hmatr_symm_arr (_nhblks_ini + 1);
         CBMatrix < _Int, _Flt > *phmatr_symm_arr = &hmatr_symm_arr[0];

         CBMatrix < _Int, _Flt >::SymmetrizeSubmatrices (_pcomm, _nhblks_ini,
                                                         _hblk2cpu_ini, _hblk2blks_ini,
                                                         pblk2hblks_ini, _blks_ini,
                                                         phmatr_ini_arr, phmatr_symm_arr);

         for (int i = 0; i < _nhblks_ini; i++) {
            if (_hblk2cpu_ini[i] == myid) {
               phmatr_ini_arr[i].ReplaceFree (phmatr_symm_arr[i]);
            }
         }

      }

// Find set of well nodes

      vector < int >ia_wells_hblks (_nhblks_ini + 1);
      vector < int >ja_3index_wells_hblks;

      {
         for (int i = 0; i <= _nhblks_ini; i++)
            ia_wells_hblks[i] = 0;
      }

      if (_params->b_new_part && _params->i_decomp_type == 5) {

         CBMatrix < _Int, _Flt >::FindWells (_pcomm, _params->nneib_well,
                                             _params->eps_dia_well,
                                             _params->thresh_max_well, _nhblks_ini,
                                             _hblk2cpu_ini, _hblk2blks_ini, _blks_ini,
                                             phmatr_ini_arr, ia_wells_hblks,
                                             ja_3index_wells_hblks);

         int *pia_wells_hblks = &ia_wells_hblks[0];

         int n_wells = pia_wells_hblks[_nhblks_ini];

         if (myid == 0 && _params->msglev >= 2) {
            cout << " Nwells ini = " << n_wells << endl;
         }

         int *pja_3index_wells_hblks = NULL;
         if (n_wells > 0)
            pja_3index_wells_hblks = &ja_3index_wells_hblks[0];

// Transform wells as blksize blocks and filter

         vector < int >ia_wells_new (_nhblks_ini + 1);
         int *pia_wells_new = &ia_wells_new[0];

         int n_wells_flt = 0;
         pia_wells_new[0] = 0;

         int i, j, ihblk, iblk, inode, inode_new;

         int ihblk_prev = -1;
         int iblk_prev = -1;
         int inode_prev = -1;

         for (i = 0; i < _nhblks_ini; i++) {
            for (j = pia_wells_hblks[i]; j < pia_wells_hblks[i + 1]; j++) {
               inode = pja_3index_wells_hblks[j * 3];
               iblk = pja_3index_wells_hblks[j * 3 + 1];
               ihblk = pja_3index_wells_hblks[j * 3 + 2];
               inode_new = inode / blksizeloc;
               if (inode_new != inode_prev || iblk != iblk_prev || ihblk != ihblk_prev) {
                  pja_3index_wells_hblks[n_wells_flt * 3] = inode_new;
                  pja_3index_wells_hblks[n_wells_flt * 3 + 1] = iblk;
                  pja_3index_wells_hblks[n_wells_flt * 3 + 2] = ihblk;
                  n_wells_flt++;
                  inode_prev = inode_new;
                  iblk_prev = iblk;
                  ihblk_prev = ihblk;
               }
            }
            pia_wells_new[i + 1] = n_wells_flt;
         }

         ia_wells_hblks.swap (ia_wells_new);

      }
// Condensed blocks partitioning

      CVectorData < long long >hblks_ini_cnd (_nhblks_ini + 1);
      CVectorData < long long >blks_ini_cnd (_nblks_ini + 1);

      long long *phblks_ini_cnd = hblks_ini_cnd.Ptr ();
      long long *pblks_ini_cnd = blks_ini_cnd.Ptr ();

      {
         int i;
         for (i = 0; i <= _nhblks_ini; i++)
            phblks_ini_cnd[i] = _hblks_ini[i] / blksizeloc;
         for (i = 0; i <= _nblks_ini; i++)
            pblks_ini_cnd[i] = _blks_ini[i] / blksizeloc;
      }

// Compute condensed sparsity and values

      vector < CBMatrix < _Int, _Flt > >hmatr_ini_cnd_arr (_nhblks_ini + 1);
      CBMatrix < _Int, _Flt > *phmatr_ini_cnd_arr = &hmatr_ini_cnd_arr[0];

      {

         int nimax = 0;

         int i;

         {
            int niloc;
            for (i = 0; i < _nblks_ini; i++) {
               niloc = (int) (_blks_ini[i + 1] - _blks_ini[i]);
               if (niloc > nimax)
                  nimax = niloc;
            }
         }

         CVectorData < int >sp2blk (nimax + 1);
         int *psp2blk = sp2blk.Ptr ();

         {
            int nimax_cnd = nimax / blksizeloc;
            int j;
            for (i = 0; i < nimax_cnd; i++) {
               for (j = 0; j < blksizeloc; j++) {
                  psp2blk[i * blksizeloc + j] = i;
               }
            }
         }

         int n_thr = 1;

#ifdef USE_THREADS
         n_thr = omp_get_max_threads ();
#endif

         vector < int >icycle_thr (n_thr + 1);
         vector < CVectorData < int > >imask_thr (n_thr + 1);

         int *picycle_thr = &icycle_thr[0];
         CVectorData < int >*pimask_thr = &imask_thr[0];

         for (i = 0; i < n_thr; i++)
            picycle_thr[i] = -1;

         _Flt fone = (_Flt) 1.0e0;

         for (i = 0; i < _nhblks_ini; i++) {
            if (_hblk2cpu_ini[i] == myid) {

               int niblk = _hblk2blks_ini[i + 1] - _hblk2blks_ini[i];

               int nzblk = phmatr_ini_arr[i].GetNzblk ();
               CMatrix < _Int, _Flt > *pASub = phmatr_ini_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str = phmatr_ini_arr[i].GetHMatrStr ();

               int *pia_hblk = phmatr_str->GetIaArr ();
               int *pja_hblk = phmatr_str->GetJaArr ();
               int *pja2_hblk = phmatr_str->GetJa2Arr ();

               phmatr_ini_cnd_arr[i].SetNzblk (nzblk);
               phmatr_ini_cnd_arr[i].ResizeASub (nzblk);

               CMatrix < _Int, _Flt > *pASub_cnd = phmatr_ini_cnd_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str_cnd =
                  phmatr_ini_cnd_arr[i].GetHMatrStr ();

               *phmatr_str_cnd = *phmatr_str;

               {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < niblk; ipar++) {

                     int my_thr = 0;
#ifdef USE_THREADS
                     my_thr = omp_get_thread_num ();
#endif

                     int icycle_th = picycle_thr[my_thr];

                     if (icycle_th == -1) {
                        pimask_thr[my_thr].resize (nimax * 6 + 1);
                        int *ppimask_thr = pimask_thr[my_thr].Ptr ();
                        int j;
                        for (j = 0; j < nimax; j++)
                           ppimask_thr[j] = -1;
                     }

                     int *ppimask_thr = pimask_thr[my_thr].Ptr ();

                     int jind;

                     for (jind = pia_hblk[ipar]; jind < pia_hblk[ipar + 1]; jind++) {

                        int jblk = pja_hblk[jind];
                        int jhblk = pja2_hblk[jind];

                        pASub[jind].CondenseSparsitySp (psp2blk, icycle_th, nimax,
                                                        ppimask_thr, pASub_cnd[jind]);

                        int nlist_loc = pASub_cnd[jind].GetNlist ();
                        int nzja_loc = pASub_cnd[jind].GetNzja ();
                        _Int *pia_loc = pASub_cnd[jind].GetIaArr ();
                        _Int *pja_loc = pASub_cnd[jind].GetJaArr ();

                        pASub_cnd[jind].ResizeA (nzja_loc);
                        pASub_cnd[jind].SetNza (nzja_loc);

                        _Flt *pa_loc = pASub_cnd[jind].GetAArr ();

                        CVector < _Flt >::SetByZeroes (nzja_loc, pa_loc);

                        if (jhblk == i && jblk == ipar) {
                           int k, kj, kkk;
                           for (k = 0; k < nlist_loc; k++) {
                              for (kj = (int) pia_loc[k]; kj < pia_loc[k + 1]; kj++) {
                                 kkk = (int) pja_loc[kj];
                                 if (kkk == k)
                                    pa_loc[kj] = fone;
                              }
                           }
                        }

                     }

                     picycle_thr[my_thr] = icycle_th;

                  }
               }
            }

         }

      }

// Combine condensed matrix data into global CSR format

      CVectorData < _Int > ia_cnd;
      CVectorData < _Int > ja_cnd;
      CVectorData < _Flt > a_cnd;

      {

         int nblk_loc_cnd = 0;
         int nloc_cnd = 0;
         int nzja_cnd = 0;

         int i;

         for (i = 0; i < _nhblks_ini; i++) {
            if (_hblk2cpu_ini[i] == myid) {

               int niblk = _hblk2blks_ini[i + 1] - _hblk2blks_ini[i];

               int nzblk = phmatr_ini_cnd_arr[i].GetNzblk ();
               CMatrix < _Int, _Flt > *pASub = phmatr_ini_cnd_arr[i].GetASubArr ();

               nblk_loc_cnd += niblk;
               nloc_cnd += (int) (phblks_ini_cnd[i + 1] - phblks_ini_cnd[i]);

               int j;

               for (j = 0; j < nzblk; j++)
                  nzja_cnd += pASub[j].GetNzja ();

            }

         }

         ia_cnd.resize (nloc_cnd + 1);
         ja_cnd.resize (nzja_cnd + 1);
         a_cnd.resize (nzja_cnd + 1);

         _Int *pia_cnd = ia_cnd.Ptr ();
         _Int *pja_cnd = ja_cnd.Ptr ();
         _Flt *pa_cnd = a_cnd.Ptr ();

         CVectorData < int >ibs_blk (_nhblks_ini + 1);
         CVectorData < int >ibs_ia (nblk_loc_cnd + 1);
         CVectorData < int >ibs_ja (nblk_loc_cnd + 1);

         int *pibs_blk = ibs_blk.Ptr ();
         int *pibs_ia = ibs_ia.Ptr ();
         int *pibs_ja = ibs_ja.Ptr ();

         nblk_loc_cnd = 0;
         nzja_cnd = 0;

         pibs_blk[0] = 0;
         pibs_ia[0] = 0;
         pibs_ja[0] = 0;

         pia_cnd[0] = 0;

         for (i = 0; i < _nhblks_ini; i++) {
            if (_hblk2cpu_ini[i] == myid) {

               int niblk = _hblk2blks_ini[i + 1] - _hblk2blks_ini[i];
               int ibegblk = _hblk2blks_ini[i];

               CMatrix < _Int, _Flt > *pASub = phmatr_ini_cnd_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str = phmatr_ini_cnd_arr[i].GetHMatrStr ();

               int *pia_hblk = phmatr_str->GetIaArr ();

               int j, k;

               for (j = 0; j < niblk; j++) {
                  for (k = pia_hblk[j]; k < pia_hblk[j + 1]; k++) {
                     nzja_cnd += pASub[k].GetNzja ();
                  }
                  pibs_ia[nblk_loc_cnd + 1] =
                     pibs_ia[nblk_loc_cnd] + (int) (pblks_ini_cnd[ibegblk + j + 1] -
                                                    pblks_ini_cnd[ibegblk + j]);
                  pibs_ja[nblk_loc_cnd + 1] = nzja_cnd;
                  nblk_loc_cnd++;
               }

            }

            pibs_blk[i + 1] = nblk_loc_cnd;

         }

         for (i = 0; i < _nhblks_ini; i++) {
            if (_hblk2cpu_ini[i] == myid) {

               int ibegblk_cnd_loc = pibs_blk[i];

               int niblk = _hblk2blks_ini[i + 1] - _hblk2blks_ini[i];
               int ibegblk = _hblk2blks_ini[i];

               CMatrix < _Int, _Flt > *pASub = phmatr_ini_cnd_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str = phmatr_ini_cnd_arr[i].GetHMatrStr ();

               int *pia_hblk = phmatr_str->GetIaArr ();
               int *pja_hblk = phmatr_str->GetJaArr ();
               int *pja2_hblk = phmatr_str->GetJa2Arr ();

               {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < niblk; ipar++) {

                     int iblkgl_cnd = ibegblk + ipar;

                     int niloc =
                        (int) (pblks_ini_cnd[iblkgl_cnd + 1] - pblks_ini_cnd[iblkgl_cnd]);

                     int ibs_i = pibs_ia[ibegblk_cnd_loc + ipar];
                     int ibs_j = pibs_ja[ibegblk_cnd_loc + ipar];

                     CVectorData < int >ia_loc (niloc + 1);
                     CVectorData < int >iptr (niloc);

                     int *pia_loc = ia_loc.Ptr ();
                     int *piptr = iptr.Ptr ();

                     int ik, jk, kk;

                     for (ik = 0; ik <= niloc; ik++)
                        pia_loc[ik] = 0;

                     for (ik = pia_hblk[ipar]; ik < pia_hblk[ipar + 1]; ik++) {
                        int nlist_temp = pASub[ik].GetNlist ();
                        _Int *plist_temp = pASub[ik].GetListArr ();
                        _Int *pia_temp = pASub[ik].GetIaArr ();
                        for (jk = 0; jk < nlist_temp; jk++) {
                           kk = (int) plist_temp[jk];
                           pia_loc[kk + 1] += (int) (pia_temp[jk + 1] - pia_temp[jk]);
                        }
                     }

                     for (ik = 0; ik < niloc; ik++)
                        pia_loc[ik + 1] += pia_loc[ik];
                     for (ik = 0; ik < niloc; ik++)
                        piptr[ik] = ibs_j + pia_loc[ik];

                     int jblk_cnd, jhblk_cnd, jblkgl_cnd, ind, jkk;
                     long long jshift;

                     for (ik = pia_hblk[ipar]; ik < pia_hblk[ipar + 1]; ik++) {
                        jblk_cnd = pja_hblk[ik];
                        jhblk_cnd = pja2_hblk[ik];
                        jblkgl_cnd = _hblk2blks_ini[jhblk_cnd] + jblk_cnd;
                        jshift = pblks_ini_cnd[jblkgl_cnd];
                        int nlist_temp = pASub[ik].GetNlist ();
                        _Int *plist_temp = pASub[ik].GetListArr ();
                        _Int *pia_temp = pASub[ik].GetIaArr ();
                        _Int *pja_temp = pASub[ik].GetJaArr ();
                        _Flt *pa_temp = pASub[ik].GetAArr ();
                        for (jk = 0; jk < nlist_temp; jk++) {
                           kk = (int) plist_temp[jk];
                           ind = piptr[kk];
                           for (jkk = (int) pia_temp[jk]; jkk < pia_temp[jk + 1]; jkk++) {
                              pja_cnd[ind] = (_Int) (jshift + pja_temp[jkk]);
                              pa_cnd[ind] = pa_temp[jkk];
                              ind++;
                           }
                           piptr[kk] = ind;
                        }
                     }

                     for (ik = 0; ik < niloc; ik++)
                        pia_cnd[ibs_i + ik + 1] = ibs_j + pia_loc[ik + 1];

                  }

               }

            }
         }

      }

      _Int *pia_cnd = ia_cnd.Ptr ();
      _Int *pja_cnd = ja_cnd.Ptr ();
      _Flt *pa_cnd = a_cnd.Ptr ();

// Create condensed solver and compute partitioning for condensed matrix

      CK3D_SolverThreads < _Int, _Flt, _FltVect > solver_cnd;

      if (_params->b_new_part) {

         if (_params->i_decomp_type == 4) {

            solver_cnd.PrepareMatrixThreads (_pcomm, _params, _stats, _nhblks_ini,
                                             _hblk2cpu_ini, _hblk2blks_ini,
                                             phblks_ini_cnd, _nblks_ini, pblks_ini_cnd,
                                             pia_cnd, pja_cnd, pa_cnd);

         } else {

            solver_cnd.PrepareMatrixWells (_pcomm, _params, _stats, _nhblks_ini,
                                           _hblk2cpu_ini, _hblk2blks_ini, phblks_ini_cnd,
                                           _nblks_ini, pblks_ini_cnd, &ia_wells_hblks,
                                           &ja_3index_wells_hblks, pia_cnd, pja_cnd,
                                           pa_cnd);

         }

      }
// Get computed data

      if (_params->b_new_part) {

         int nhblks_cnd = solver_cnd.GetNhblks ();
         int nblks_cnd = solver_cnd.GetNblks ();
         long long *phblks_cnd = solver_cnd.GetHBlks ();
         int *phblk2cpu_cnd = solver_cnd.GetHBlk2cpu ();
         int *phblk2blks_cnd = solver_cnd.GetHBlk2blks ();
         long long *pblks_cnd = solver_cnd.GetBlks ();

// Transform computed data into the extended form

         this->nhblks = nhblks_cnd;
         this->nblks = nblks_cnd;

         this->hblks.resize (nhblks_cnd + 1);
         this->hblk2cpu.resize (nhblks_cnd + 1);
         this->hblk2blks.resize (nhblks_cnd + 1);
         this->blk2hblks.resize (nblks_cnd + 1);
         this->blks.resize (nblks_cnd + 1);

         long long *phblks = &this->hblks[0];
         int *phblk2cpu = &this->hblk2cpu[0];
         int *phblk2blks = &this->hblk2blks[0];
         int *pblk2hblks = &this->blk2hblks[0];
         long long *pblks = &this->blks[0];

         {
            int i;
            for (i = 0; i <= nhblks_cnd; i++)
               phblks[i] = phblks_cnd[i] * blksizeloc;
            for (i = 0; i < nhblks_cnd; i++)
               phblk2cpu[i] = phblk2cpu_cnd[i];
            for (i = 0; i <= nhblks_cnd; i++)
               phblk2blks[i] = phblk2blks_cnd[i];
            for (i = 0; i <= nblks_cnd; i++)
               pblks[i] = pblks_cnd[i] * blksizeloc;

            int j;

            for (i = 0; i < nhblks_cnd; i++) {
               int ibegblk = phblk2blks[i];
               int iendblk = phblk2blks[i + 1] - 1;
               for (j = ibegblk; j <= iendblk; j++)
                  pblk2hblks[j] = i;
            }

         }

      }

      int nhblks_loc = this->nhblks;
      int nblks_loc = this->nblks;
      long long *phblks = &this->hblks[0];
      int *phblk2cpu = &this->hblk2cpu[0];
      int *phblk2blks = &this->hblk2blks[0];
      int *pblk2hblks = &this->blk2hblks[0];
      long long *pblks = &this->blks[0];

// Transform ordering

      if (_params->b_new_part) {

         if (_params->i_decomp_type == 4) {

            CVectorData < int >*porder_ini_cnd = &solver_cnd.order_ini[0];

            this->order_ini.resize (nhblks_loc + 1);

            CVectorData < int >*porder_ini = &this->order_ini[0];

            {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
               for (int ipar = 0; ipar < nhblks_loc; ipar++) {
                  if (phblk2cpu[ipar] == myid) {
                     int *pporder_ini_cnd = porder_ini_cnd[ipar].Ptr ();
                     int ni_loc = (int) (phblks[ipar + 1] - phblks[ipar]);
                     int ni_cnd = ni_loc / blksizeloc;
                     porder_ini[ipar].resize (ni_loc);
                     int *pporder_ini = porder_ini[ipar].Ptr ();
                     int i, j, jj;
                     for (i = 0; i < ni_cnd; i++) {
                        jj = pporder_ini_cnd[i];
                        for (j = 0; j < blksizeloc; j++) {
                           pporder_ini[i * blksizeloc + j] = jj * blksizeloc + j;
                        }
                     }
                  }
               }

            }

         } else {

            int nz_ord = 0;

            int i;

            for (i = 0; i < _nhblks_ini; i++) {
               if (_hblk2cpu_ini[i] == myid) {
                  nz_ord += (int) (_hblks_ini[i + 1] - _hblks_ini[i]);
//            for (i=0;i<nhblks_loc;i++) {
//               if (phblk2cpu[i] == myid) {
//                  nz_ord += (int)(phblks[i+1]-phblks[i]);
               }
            }

            this->order_wells.resize (nz_ord + 1);
            this->order2ind_wells.resize (2 * nz_ord + 1);

            long long *porder_wells = this->order_wells.Ptr ();
            int *porder2ind_wells = this->order2ind_wells.Ptr ();

            int nz_ord_cnd = nz_ord / blksizeloc;

            long long *porder_wells_cnd = solver_cnd.order_wells.Ptr ();
            int *porder2ind_wells_cnd = solver_cnd.order2ind_wells.Ptr ();

            int n_thr = 1;

#ifdef USE_THREADS
            n_thr = omp_get_max_threads ();
#endif

            int ni_part_cnd = nz_ord_cnd / n_thr;

            {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
               for (int ipar = 0; ipar < n_thr; ipar++) {

                  int ibeg_cnd = ipar * ni_part_cnd;
                  int iend_cnd = (ipar + 1) * ni_part_cnd - 1;
                  if (ipar == n_thr - 1)
                     iend_cnd = nz_ord_cnd - 1;

                  int j, k, jj, jhblk;
                  long long jj_long;

                  for (j = ibeg_cnd; j <= iend_cnd; j++) {
                     jj_long = porder_wells_cnd[j];
                     for (k = 0; k < blksizeloc; k++) {
                        porder_wells[j * blksizeloc + k] = jj_long * blksizeloc + k;
                     }
                     jj = porder2ind_wells_cnd[j * 2];
                     jhblk = porder2ind_wells_cnd[j * 2 + 1];
                     for (k = 0; k < blksizeloc; k++) {
                        porder2ind_wells[j * blksizeloc * 2 + k * 2] =
                           jj * blksizeloc + k;
                        porder2ind_wells[j * blksizeloc * 2 + k * 2 + 1] = jhblk;
                     }
                  }


               }

            }

         }

      }
// Set control params

      if (_params->b_new_part) {
         if (_params->i_decomp_type == 4) {
            this->b_use_ini = true;
         } else if (_params->i_decomp_type == 5) {
            this->b_use_ini = false;
            this->b_use_wells = true;
            this->b_blk_wells = false;
            if (solver_cnd.b_blk_wells) {
               this->b_blk_wells = true;
            }
         }
         if (_params->b_fast_transform) {
            this->PrepareFastTransform (phmatr_ini_arr);
         }
      }
// Apply block diagonal ordering

      this->hmatr_arr.resize (nhblks_loc + 1);

      CBMatrix < _Int, _Flt > *phmatr_arr = &this->hmatr_arr[0];

      if (this->b_fast_transform) {

         this->UseFastTransform (phmatr_ini_arr);

      } else {

         if (_params->i_decomp_type == 4) {

            CVectorData < int >*porder_ini = &this->order_ini[0];

            CBMatrix < _Int, _Flt >::ReorderHMatrixDiag (_pcomm, _nhblks_ini,
                                                         _hblk2cpu_ini, _hblk2blks_ini,
                                                         _nblks_ini, _blks_ini,
                                                         phmatr_ini_arr, porder_ini,
                                                         phblk2blks, nblks_loc, pblks,
                                                         phmatr_arr);

         } else {

            long long *porder_wells = this->order_wells.Ptr ();

            CBMatrix < _Int, _Flt >::ReorderHMatrix (_pcomm, _nhblks_ini, _hblk2cpu_ini,
                                                     _hblk2blks_ini, _blks_ini,
                                                     phmatr_ini_arr, porder_wells,
                                                     nhblks_loc, phblk2cpu, phblk2blks,
                                                     pblks, phmatr_arr);

         }

      }

      phmatr_arr = &this->hmatr_arr[0];

// Init MvmA structures

      this->mvm.InitControl (_pcomm, nhblks_loc, phblk2cpu, phblk2blks, pblk2hblks,
                             phblks, nblks_loc, pblks);

      this->mvm.InitMvmA (phmatr_arr);

      this->b_use_blksize = false;
      this->blksize_bscl = blksizeloc;

      if (_params->b_new_part) {
         if (_params->i_decomp_type == 4) {
            this->b_use_ini = true;
         } else if (_params->i_decomp_type == 5) {
            this->b_use_ini = false;
            this->b_use_wells = true;
            this->b_blk_wells = false;
            if (solver_cnd.b_blk_wells) {
               this->b_blk_wells = true;
            }
         }
      }

   }

// Get local cpu vector size, current (with vtype == 1) or initial (with vtype == 0)
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > long long CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::GetVSize (int _vtype)
   {

      int myid = CMPIDataExchange::GetMyid (this->pcomm);

      long long isize = 0;

      if (_vtype == 1) {
         int nhblks_loc = this->GetNhblks ();
         long long *phblks_loc = this->GetHBlks ();
         int *phblk2cpu_loc = this->GetHBlk2cpu ();
         int i;
         for (i = 0; i < nhblks_loc; i++) {
            if (phblk2cpu_loc[i] == myid) {
               isize += (phblks_loc[i + 1] - phblks_loc[i]);
            }
         }
      } else {
         int nhblks_loc;
         long long *phblks_loc;
         int *phblk2cpu_loc;
         if (!this->b_use_wells) {
            nhblks_loc = this->GetNhblks ();
            phblks_loc = this->GetHBlks ();
            phblk2cpu_loc = this->GetHBlk2cpu ();
         } else {
            nhblks_loc = this->GetNhblksIni ();
            phblks_loc = this->GetHBlksIni ();
            phblk2cpu_loc = this->GetHBlk2cpuIni ();
         }
         int i;
         for (i = 0; i < nhblks_loc; i++) {
            if (phblk2cpu_loc[i] == myid) {
               isize += (phblks_loc[i + 1] - phblks_loc[i]);
            }
         }
      }

      return isize;

   }

// Store control data only
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::StoreControlData (void *_pcomm, int _nhblks, int *_hblk2cpu,
                                    int *_hblk2blks, long long *_hblks, int _nblks,
                                    long long *_blks)
   {

// Store control data

      this->pcomm = _pcomm;
      this->nhblks = _nhblks;
      this->nblks = _nblks;

      this->hblks.resize (_nhblks + 1);
      this->hblk2cpu.resize (_nhblks + 1);
      this->hblk2blks.resize (_nhblks + 1);
      this->blk2hblks.resize (_nblks + 1);
      this->blks.resize (_nblks + 1);

      long long *phblks = &this->hblks[0];
      int *phblk2cpu = &this->hblk2cpu[0];
      int *phblk2blks = &this->hblk2blks[0];
      int *pblk2hblks = &this->blk2hblks[0];
      long long *pblks = &this->blks[0];

      int i, j;

      for (i = 0; i <= _nhblks; i++)
         phblks[i] = _hblks[i];
      for (i = 0; i < _nhblks; i++)
         phblk2cpu[i] = _hblk2cpu[i];
      for (i = 0; i <= _nhblks; i++)
         phblk2blks[i] = _hblk2blks[i];
      for (i = 0; i < _nhblks; i++) {
         for (j = phblk2blks[i]; j < phblk2blks[i + 1]; j++) {
            pblk2hblks[j] = i;
         }
      }

      for (i = 0; i <= _nblks; i++)
         pblks[i] = _blks[i];

      this->mvm.InitControl (_pcomm, _nhblks, phblk2cpu, phblk2blks, pblk2hblks, phblks,
                             _nblks, pblks);
      this->slv_float.InitControl (_pcomm, _nhblks, phblks, phblk2cpu, phblk2blks,
                                   pblk2hblks, pblks);
      this->slv_double.InitControl (_pcomm, _nhblks, phblks, phblk2cpu, phblk2blks,
                                    pblk2hblks, pblks);

   }

// Prepare matrix data only
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::PrepareMatrixBase (void *_pcomm, int _nhblks, int *_hblk2cpu,
                                     int *_hblk2blks, long long *_hblks, int _nblks,
                                     long long *_blks, _Int * _ia, _Int * _ja, _Flt * _a)
   {

// Store control data

      this->pcomm = _pcomm;
      this->nhblks = _nhblks;
      this->nblks = _nblks;

      this->hblks.resize (_nhblks + 1);
      this->hblk2cpu.resize (_nhblks + 1);
      this->hblk2blks.resize (_nhblks + 1);
      this->blk2hblks.resize (_nblks + 1);
      this->blks.resize (_nblks + 1);

      long long *phblks = &this->hblks[0];
      int *phblk2cpu = &this->hblk2cpu[0];
      int *phblk2blks = &this->hblk2blks[0];
      int *pblk2hblks = &this->blk2hblks[0];
      long long *pblks = &this->blks[0];

      int i, j;

      for (i = 0; i <= _nhblks; i++)
         phblks[i] = _hblks[i];
      for (i = 0; i < _nhblks; i++)
         phblk2cpu[i] = _hblk2cpu[i];
      for (i = 0; i <= _nhblks; i++)
         phblk2blks[i] = _hblk2blks[i];
      for (i = 0; i < _nhblks; i++) {
         for (j = phblk2blks[i]; j < phblk2blks[i + 1]; j++) {
            pblk2hblks[j] = i;
         }
      }

      for (i = 0; i <= _nblks; i++)
         pblks[i] = _blks[i];

      int myid = CMPIDataExchange::GetMyid (_pcomm);

// Allocate arrays for matrix data

      this->hmatr_arr.resize (_nhblks + 1);

      CBMatrix < _Int, _Flt > *phmatr_arr = &this->hmatr_arr[0];

// Create matrix as set of hblocks

      int nimax = 0;

      int niloc;

      for (i = 0; i < _nblks; i++) {
         niloc = (int) (pblks[i + 1] - pblks[i]);
         if (niloc > nimax)
            nimax = niloc;
      }

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      vector < int >icycle_thr (n_thr + 1);
      vector < CVectorData < _Int > >listloc_thr (n_thr + 1);
      vector < CVectorData < int > >imaskblk_thr (n_thr + 1);

      int *picycle_thr = &icycle_thr[0];
      CVectorData < _Int > *plistloc_thr = &listloc_thr[0];
      CVectorData < int >*pimaskblk_thr = &imaskblk_thr[0];

      for (i = 0; i < n_thr; i++)
         picycle_thr[i] = -1;

      int ibeg = 0;

      for (i = 0; i < _nhblks; i++) {
         if (phblk2cpu[i] == myid) {

            int niblk = phblk2blks[i + 1] - phblk2blks[i];
            vector < CBMatrix < _Int, _Flt > >hblk_arr (niblk + 1);
            CBMatrix < _Int, _Flt > *phblk_arr = &hblk_arr[0];
            int iblk0 = phblk2blks[i];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = phblk2blks[i]; ipar < phblk2blks[i + 1]; ipar++) {

               int my_thr = 0;
#ifdef USE_THREADS
               my_thr = omp_get_thread_num ();
#endif

               int j;

               if (picycle_thr[my_thr] == -1) {
                  plistloc_thr[my_thr].resize (2 * nimax + 1);
                  pimaskblk_thr[my_thr].resize (3 * _nblks + 1);
                  int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
                  for (j = 0; j < _nblks; j++)
                     pimaskblk[j] = -1;
               }
               int icycleblk = picycle_thr[my_thr];
               _Int *plistloc = plistloc_thr[my_thr].Ptr ();
               _Int *pialoc = plistloc + nimax;
               int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
               icycleblk++;
               int niloc = (int) (pblks[ipar + 1] - pblks[ipar]);
               int ibs = (int) (pblks[ipar] - pblks[iblk0]);
               int ishift = (int) _ia[ibeg + ibs];
               for (j = 0; j < niloc; j++) {
                  plistloc[j] = (_Int) (j + pblks[ipar]);
               }
               for (j = 0; j <= niloc; j++) {
                  pialoc[j] = _ia[ibeg + ibs + j] - _ia[ibeg + ibs];
               }
               CBMatrix < _Int, _Flt > hblk_loc (ipar, niloc, plistloc, pialoc,
                                                 _ja + ishift, _a + ishift, _nblks, pblks,
                                                 icycleblk, pimaskblk);
               picycle_thr[my_thr] = icycleblk;
               phblk_arr[ipar - iblk0].ReplaceFree (hblk_loc);
            }

// Combine array blocks into one hblock

            int nzja_tot = 0;

            for (j = 0; j < niblk; j++)
               nzja_tot += phblk_arr[j].GetNzblk ();

            CBMatrix < _Int, _Flt > hblk;

            hblk.ResizeASub (nzja_tot);
            hblk.SetNzblk (nzja_tot);

            CMatrix < _Int, _Flt > *pASub = hblk.GetASubArr ();

            CMatrix < int, float >*phmatr_str = hblk.GetHMatrStr ();

            phmatr_str->ResizeAndSetAllSp (niblk, niblk, nzja_tot, nzja_tot);

            int *plist_str = phmatr_str->GetListArr ();
            int *plist2_str = phmatr_str->GetList2Arr ();
            int *pia_str = phmatr_str->GetIaArr ();
            int *pja_str = phmatr_str->GetJaArr ();
            int *pja2_str = phmatr_str->GetJa2Arr ();

            for (j = 0; j < niblk; j++)
               plist_str[j] = j;
            for (j = 0; j < niblk; j++)
               plist2_str[j] = i;

            nzja_tot = 0;

            pia_str[0] = 0;

            int k, kk, kk2;

            for (j = 0; j < niblk; j++) {
               CMatrix < _Int, _Flt > *pASub_temp = phblk_arr[j].GetASubArr ();
               CMatrix < int, float >*phmatr_str_temp = phblk_arr[j].GetHMatrStr ();
               int nzja_temp = phmatr_str_temp->GetNzja ();
               int *pja_temp = phmatr_str_temp->GetJaArr ();
               for (k = 0; k < nzja_temp; k++) {
                  kk = pja_temp[k];
                  kk2 = pblk2hblks[kk];
                  pja_str[nzja_tot] = kk - phblk2blks[kk2];
                  pja2_str[nzja_tot] = kk2;
                  pASub[nzja_tot].ReplaceFree (pASub_temp[k]);
                  nzja_tot++;
               }
               pia_str[j + 1] = nzja_tot;
            }
            phmatr_arr[i].ReplaceFree (hblk);
            ibeg += (int) (phblks[i + 1] - phblks[i]);
         }
      }

// Symmetrize hmatr

      vector < CBMatrix < _Int, _Flt > >hmatr_symm_arr (_nhblks + 1);
      CBMatrix < _Int, _Flt > *phmatr_symm_arr = &hmatr_symm_arr[0];

      CBMatrix < _Int, _Flt >::SymmetrizeSubmatrices (_pcomm, _nhblks, phblk2cpu,
                                                      phblk2blks, pblk2hblks, pblks,
                                                      phmatr_arr, phmatr_symm_arr);

      for (i = 0; i < _nhblks; i++) {
         if (phblk2cpu[i] == myid) {
            phmatr_arr[i].ReplaceFree (phmatr_symm_arr[i]);
         }
      }

// Init MvmA structures

      this->mvm.InitControl (_pcomm, _nhblks, phblk2cpu, phblk2blks, pblk2hblks, phblks,
                             _nblks, pblks);
      this->mvm.InitMvmA (phmatr_arr);

      this->b_fast_transform = false;
      this->b_use_ini = false;
      this->b_use_wells = false;
      this->b_use_blksize = false;
      this->b_blk_wells = false;

   }

// Perform iterations of the iterative scheme
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::SolveIter (SParams * _params, SStatData * _stats, _FltVect * _rhs,
                             _FltVect * _sol)
   {

      vector < _FltVect > rhs_new;
      vector < _FltVect > sol_new;

      this->TransformVectorsForward (_rhs, _sol, rhs_new, sol_new);

      _FltVect *prhs_new = NULL;
      _FltVect *psol_new = NULL;

      if (rhs_new.size () > 0)
         prhs_new = &rhs_new[0];
      if (sol_new.size () > 0)
         psol_new = &sol_new[0];

      _params->blksize_iter = 1;

      this->SolveIter (_params, _stats, this, CK3D_SolverThreads < _Int, _Flt,
                       _FltVect >::MvmA_static, this, CK3D_SolverThreads < _Int, _Flt,
                       _FltVect >::SlvLU_static, prhs_new, psol_new);

      this->TransformVectorsBackward (rhs_new, sol_new, _rhs, _sol);

   }

// Transform vectors forward
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::TransformVectorsForward (_FltVect * _rhs_ini, _FltVect * _sol_ini,
                                           vector < _FltVect > &_rhs_new,
                                           vector < _FltVect > &_sol_new)
   {

      long long isize_curr = this->GetVSize (1);

      _rhs_new.resize ((int) isize_curr + 1);
      _sol_new.resize ((int) isize_curr + 1);

      _FltVect *p_rhs_new = &_rhs_new[0];
      _FltVect *p_sol_new = &_sol_new[0];

      if (!this->b_use_wells) {

         if (this->b_use_ini) {
            CVector < _FltVect >::CopyVector_thr ((int) isize_curr, _rhs_ini, p_rhs_new);
            CVector < _FltVect >::CopyVector_thr ((int) isize_curr, _sol_ini, p_sol_new);
            this->ReorderVectorDataIni ('D', p_rhs_new, p_sol_new);
         } else {
            CVector < _FltVect >::CopyVector_thr ((int) isize_curr, _rhs_ini, p_rhs_new);
            CVector < _FltVect >::CopyVector_thr ((int) isize_curr, _sol_ini, p_sol_new);
         }

      } else {

         CVectorData < _FltVect > rhs_ord;
         CVectorData < _FltVect > sol_ord;

         this->ReorderVectorDataWells ('D', _rhs_ini, _sol_ini, rhs_ord, sol_ord);

         _FltVect *prhs_ord = rhs_ord.Ptr ();
         _FltVect *psol_ord = sol_ord.Ptr ();

         CVector < _FltVect >::CopyVector_thr ((int) isize_curr, prhs_ord, p_rhs_new);
         CVector < _FltVect >::CopyVector_thr ((int) isize_curr, psol_ord, p_sol_new);

      }

   }

// Transform vectors backward
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::TransformVectorsBackward (vector < _FltVect > &_rhs_new,
                                            vector < _FltVect > &_sol_new,
                                            _FltVect * _rhs_fin, _FltVect * _sol_fin)
   {

      long long isize_curr = this->GetVSize (1);

      _FltVect *p_rhs_new = &_rhs_new[0];
      _FltVect *p_sol_new = &_sol_new[0];

      if (!this->b_use_wells) {

         CVector < _FltVect >::CopyVector_thr ((int) isize_curr, p_rhs_new, _rhs_fin);
         CVector < _FltVect >::CopyVector_thr ((int) isize_curr, p_sol_new, _sol_fin);

         if (this->b_use_ini) {
            this->ReorderVectorDataIni ('I', _rhs_fin, _sol_fin);
         }

      } else {

         CVectorData < _FltVect > rhs_ord ((int) isize_curr);
         CVectorData < _FltVect > sol_ord ((int) isize_curr);

         _FltVect *prhs_ord = rhs_ord.Ptr ();
         _FltVect *psol_ord = sol_ord.Ptr ();

         CVector < _FltVect >::CopyVector_thr ((int) isize_curr, p_rhs_new, prhs_ord);
         CVector < _FltVect >::CopyVector_thr ((int) isize_curr, p_sol_new, psol_ord);

         this->ReorderVectorDataWells ('I', _rhs_fin, _sol_fin, rhs_ord, sol_ord);

      }

   }

// Prepare fast transform data
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::PrepareFastTransform (CBMatrix < _Int, _Flt > *_hmatr_ini_arr)
   {

// Get control data

      void *pcommloc = this->pcomm;

      int nhblks_ini_loc = this->nhblks_ini;
      int nhblks_loc = this->nhblks;

      int *phblk2cpu_ini = &this->hblk2cpu_ini[0];
      int *phblk2blks_ini = &this->hblk2blks_ini[0];
      long long *pblks_ini = &this->blks_ini[0];

      int *phblk2cpu = &this->hblk2cpu[0];
      int *phblk2blks = &this->hblk2blks[0];
      long long *pblks = &this->blks[0];

      long long *porder_wells = this->order_wells.Ptr ();

      int myid = CMPIDataExchange::GetMyid (pcommloc);

      bool b_order_dia = false;
      if (!this->b_use_ini && !this->b_use_wells) {
         b_order_dia = true;
      }

      CVectorData < int >*porder_dia = NULL;
      if (b_order_dia) {
         porder_dia = &this->order_ini[0];
      }
// Compute inverse order

      CVectorData < long long >iorder;

      if (b_order_dia) {
         CBMatrix < _Int, _Flt >::InverseOrderDiag (pcommloc, nhblks_ini_loc,
                                                    phblk2cpu_ini, phblk2blks_ini,
                                                    pblks_ini, porder_dia, nhblks_loc,
                                                    phblk2cpu, phblk2blks, pblks, iorder);
      } else {
         CBMatrix < _Int, _Flt >::InverseOrder (pcommloc, nhblks_ini_loc, phblk2cpu_ini,
                                                phblk2blks_ini, pblks_ini, porder_wells,
                                                nhblks_loc, phblk2cpu, phblk2blks, pblks,
                                                iorder);
      }

      long long *piorder = iorder.Ptr ();

// Symmetrize hmatr

      vector < int >blk2hblks_ini (nblks_ini + 1);
      int *pblk2hblks_ini = &blk2hblks_ini[0];

      {

         int i, j;

         for (i = 0; i < nhblks_ini_loc; i++) {
            for (j = phblk2blks_ini[i]; j < phblk2blks_ini[i + 1]; j++) {
               pblk2hblks_ini[j] = i;
            }
         }

      }

      vector < CBMatrix < _Int, _Flt > >hmatr_symm_arr (nhblks_ini_loc + 1);
      CBMatrix < _Int, _Flt > *phmatr_symm_arr = &hmatr_symm_arr[0];

      CBMatrix < _Int, _Flt >::SymmetrizeSubmatrices (pcommloc, nhblks_ini_loc,
                                                      phblk2cpu_ini, phblk2blks_ini,
                                                      pblk2hblks_ini, pblks_ini,
                                                      _hmatr_ini_arr, phmatr_symm_arr);

// Prepare sparsity to reorder

      this->hmatr_ini_FT_arr.resize (nhblks_ini_loc + 1);
      CBMatrix < _Int, double >*phmatr_ini_FT_arr = &this->hmatr_ini_FT_arr[0];

      {
         int i;
         for (i = 0; i < nhblks_ini_loc; i++) {
            if (phblk2cpu_ini[i] == myid) {
               int nblks_loc = phblk2blks_ini[i + 1] - phblk2blks_ini[i];

               int nzblk = phmatr_symm_arr[i].GetNzblk ();
               CMatrix < _Int, _Flt > *pASub = phmatr_symm_arr[i].GetASubArr ();
               CMatrix < int, float >*pHMatrStr = phmatr_symm_arr[i].GetHMatrStr ();

               phmatr_ini_FT_arr[i].ResizeASub (nzblk);
               phmatr_ini_FT_arr[i].SetNzblk (nzblk);
               CMatrix < _Int, double >*pASub_dbl = phmatr_ini_FT_arr[i].GetASubArr ();
               CMatrix < int, float >*pHMatrStr_dbl = phmatr_ini_FT_arr[i].GetHMatrStr ();

               *pHMatrStr_dbl = *pHMatrStr;

               int *pia_hmatr = pHMatrStr->GetIaArr ();

               {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < nblks_loc; ipar++) {
                     int j, k;
                     for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
                        int nlist_temp = pASub[j].GetNlist ();
                        int nzja_temp = pASub[j].GetNzja ();
                        _Int *plist_temp = pASub[j].GetListArr ();
                        _Int *pia_temp = pASub[j].GetIaArr ();
                        _Int *pja_temp = pASub[j].GetJaArr ();
                        pASub_dbl[j].ResizeAndSetAll (nlist_temp, 0, nzja_temp, 0,
                                                      nzja_temp);
                        _Int *plist_dbl = pASub_dbl[j].GetListArr ();
                        _Int *pia_dbl = pASub_dbl[j].GetIaArr ();
                        _Int *pja_dbl = pASub_dbl[j].GetJaArr ();
                        double *pa_dbl = pASub_dbl[j].GetAArr ();
                        for (k = 0; k < nlist_temp; k++)
                           plist_dbl[k] = plist_temp[k];
                        for (k = 0; k <= nlist_temp; k++)
                           pia_dbl[k] = pia_temp[k];
                        for (k = 0; k < nzja_temp; k++)
                           pja_dbl[k] = pja_temp[k];
                        for (k = 0; k < nzja_temp; k++)
                           pa_dbl[k] = 0.0e0;
                     }
                  }
               }
            }
         }
      }

// Reorder forward sparsity and zeroes of initial matrix

      this->hmatr_fin_FT_arr.resize (nhblks_loc + 1);
      CBMatrix < _Int, double >*phmatr_fin_FT_arr = &this->hmatr_fin_FT_arr[0];

      if (b_order_dia) {
         int nblks_ini = phblk2blks_ini[nhblks_ini_loc];
         int nblks_fin = phblk2blks[nhblks_loc];
         CBMatrix < _Int, double >::ReorderHMatrixDiag (pcommloc, nhblks_ini_loc,
                                                        phblk2cpu_ini, phblk2blks_ini,
                                                        nblks_ini, pblks_ini,
                                                        phmatr_ini_FT_arr, porder_dia,
                                                        phblk2blks, nblks_fin, pblks,
                                                        phmatr_fin_FT_arr);
      } else {
         CBMatrix < _Int, double >::ReorderHMatrix (pcommloc, nhblks_ini_loc,
                                                    phblk2cpu_ini, phblk2blks_ini,
                                                    pblks_ini, phmatr_ini_FT_arr,
                                                    porder_wells, nhblks_loc, phblk2cpu,
                                                    phblk2blks, pblks, phmatr_fin_FT_arr);
      }

      {
         vector < CBMatrix < _Int, double > >hmatr_temp;
         hmatr_ini_FT_arr.swap (hmatr_temp);
      }

// Fill reordered double data

      {
         int i;
         for (i = 0; i < nhblks_loc; i++) {
            if (phblk2cpu[i] == myid) {
               int nblks_loc = phblk2blks[i + 1] - phblk2blks[i];

               CMatrix < _Int, double >*pASub = phmatr_fin_FT_arr[i].GetASubArr ();
               CMatrix < int, float >*pHMatrStr = phmatr_fin_FT_arr[i].GetHMatrStr ();

               int *pia_hmatr = pHMatrStr->GetIaArr ();

               {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < nblks_loc; ipar++) {
                     int j, k;
                     SIndFT *pind = NULL;
                     for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
                        int nzja_temp = pASub[j].GetNzja ();
                        double *pa_temp = pASub[j].GetAArr ();
                        for (k = 0; k < nzja_temp; k++) {
                           pind = (SIndFT *) (pa_temp + k);
                           pind->ihblk_FT = (short) i;
                           pind->ind_blk_FT = (short) j;
                           pind->ind_elem_FT = (int) k;
                        }
                     }
                  }

               }
            }
         }
      }

// Reorder backward reordered matrix with filled double data

      this->hmatr_ini_FT_arr.resize (nhblks_ini_loc + 1);
      phmatr_ini_FT_arr = &this->hmatr_ini_FT_arr[0];

      CBMatrix < _Int, double >::ReorderHMatrix (pcommloc, nhblks_loc, phblk2cpu,
                                                 phblk2blks, pblks, phmatr_fin_FT_arr,
                                                 piorder, nhblks_ini_loc, phblk2cpu_ini,
                                                 phblk2blks_ini, pblks_ini,
                                                 phmatr_ini_FT_arr);

// Free double data in final FT blocks

      {
         int i;
         for (i = 0; i < nhblks_loc; i++) {
            if (phblk2cpu[i] == myid) {
               int nblks_loc = phblk2blks[i + 1] - phblk2blks[i];

               CMatrix < _Int, double >*pASub = phmatr_fin_FT_arr[i].GetASubArr ();
               CMatrix < int, float >*pHMatrStr = phmatr_fin_FT_arr[i].GetHMatrStr ();

               int *pia_hmatr = pHMatrStr->GetIaArr ();

               {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < nblks_loc; ipar++) {
                     int j;
                     for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
                        vector < double >*pa_temp = pASub[j].GetA ();
                        vector < double >a_dummy;
                        pa_temp->swap (a_dummy);
                     }
                  }

               }
            }
         }
      }

      this->b_fast_transform = true;

   }

// Use fast transform data
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::UseFastTransform (CBMatrix < _Int, _Flt > *_hmatr_ini_arr)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

// Get control data

      void *pcommloc = this->pcomm;

      int myid = CMPIDataExchange::GetMyid (pcommloc);

      int nhblks_ini_loc = this->nhblks_ini;
      int nhblks_loc = this->nhblks;

      int *phblk2cpu_ini = &this->hblk2cpu_ini[0];
      int *phblk2blks_ini = &this->hblk2blks_ini[0];

      int *phblk2cpu = &this->hblk2cpu[0];
      int *phblk2blks = &this->hblk2blks[0];

      CBMatrix < _Int, double >*phmatr_ini_FT_arr = &this->hmatr_ini_FT_arr[0];
      CBMatrix < _Int, double >*phmatr_fin_FT_arr = &this->hmatr_fin_FT_arr[0];

// Init resulting matrix by zero values

      _Flt fzero = (_Flt) 0.0e0;

      vector < CBMatrix < _Int, _Flt > >hmatr_ord (nhblks_loc + 1);
      CBMatrix < _Int, _Flt > *phmatr_ord = &hmatr_ord[0];

      {
         int i;
         for (i = 0; i < nhblks_loc; i++) {
            if (phblk2cpu[i] == myid) {
               int nblks_loc = phblk2blks[i + 1] - phblk2blks[i];

               int nzblk = phmatr_fin_FT_arr[i].GetNzblk ();
               CMatrix < _Int, double >*pASub_dbl = phmatr_fin_FT_arr[i].GetASubArr ();
               CMatrix < int, float >*pHMatrStr_dbl = phmatr_fin_FT_arr[i].GetHMatrStr ();

               phmatr_ord[i].ResizeASub (nzblk);
               phmatr_ord[i].SetNzblk (nzblk);
               CMatrix < _Int, _Flt > *pASub = phmatr_ord[i].GetASubArr ();
               CMatrix < int, float >*pHMatrStr = phmatr_ord[i].GetHMatrStr ();

               *pHMatrStr = *pHMatrStr_dbl;

               int *pia_hmatr = pHMatrStr->GetIaArr ();

               {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < nblks_loc; ipar++) {
                     int j, k;
                     for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
                        int nlist_temp = pASub_dbl[j].GetNlist ();
                        int nzja_temp = pASub_dbl[j].GetNzja ();
                        _Int *plist_temp = pASub_dbl[j].GetListArr ();
                        _Int *pia_temp = pASub_dbl[j].GetIaArr ();
                        _Int *pja_temp = pASub_dbl[j].GetJaArr ();
                        pASub[j].ResizeAndSetAll (nlist_temp, 0, nzja_temp, 0, nzja_temp);
                        _Int *plist_flt = pASub[j].GetListArr ();
                        _Int *pia_flt = pASub[j].GetIaArr ();
                        _Int *pja_flt = pASub[j].GetJaArr ();
                        _Flt *pa_flt = pASub[j].GetAArr ();
                        for (k = 0; k < nlist_temp; k++)
                           plist_flt[k] = plist_temp[k];
                        for (k = 0; k <= nlist_temp; k++)
                           pia_flt[k] = pia_temp[k];
                        for (k = 0; k < nzja_temp; k++)
                           pja_flt[k] = pja_temp[k];
                        for (k = 0; k < nzja_temp; k++)
                           pa_flt[k] = fzero;
                     }
                  }
               }
            }
         }
      }

// Scan pair of hblock matrices, store local data and prepare initial send data for nonlocal ones

      CVectorData < std::map < int, vector < pair < int, pair < int,
         _Flt > > > > >map_hblk2st_threads (n_thr + 1);
      std::map < int, vector < pair < int, pair < int, _Flt > > > >*pmap_hblk2st_threads =
         map_hblk2st_threads.Ptr ();

      {

         int i;
         for (i = 0; i < nhblks_ini_loc; i++) {
            if (phblk2cpu_ini[i] == myid) {

               int nblks_loc = phblk2blks_ini[i + 1] - phblk2blks_ini[i];

               CMatrix < int, float >*pHMatrStr = _hmatr_ini_arr[i].GetHMatrStr ();
               CMatrix < _Int, _Flt > *pASub = _hmatr_ini_arr[i].GetASubArr ();

               int *pia_hmatr = pHMatrStr->GetIaArr ();
               int *pja_hmatr = pHMatrStr->GetJaArr ();
               int *pja2_hmatr = pHMatrStr->GetJa2Arr ();

               CMatrix < int, float >*pHMatrStr_dbl = phmatr_ini_FT_arr[i].GetHMatrStr ();
               CMatrix < _Int, double >*pASub_dbl = phmatr_ini_FT_arr[i].GetASubArr ();

               int *pia_hmatr_dbl = pHMatrStr_dbl->GetIaArr ();
               int *pja_hmatr_dbl = pHMatrStr_dbl->GetJaArr ();
               int *pja2_hmatr_dbl = pHMatrStr_dbl->GetJa2Arr ();

               {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < nblks_loc; ipar++) {
                     int my_thr = 0;
#ifdef USE_THREADS
                     my_thr = omp_get_thread_num ();
#endif
                     if (pia_hmatr[ipar + 1] - pia_hmatr[ipar] !=
                         pia_hmatr_dbl[ipar + 1] - pia_hmatr_dbl[ipar]) {
                        cout <<
                           " CK3D_SolverThreads <>::UseFastTransform: error 1: incorrect block sparsity !!! "
                           << endl;
                        throw
                           " CK3D_SolverThreads <>::UseFastTransform: error 1: incorrect block sparsity !!! ";
                     }
                     int j, jblk1, jhblk1, jblk2, jhblk2;
                     for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
                        jblk1 = pja_hmatr[j];
                        jhblk1 = pja2_hmatr[j];
                        jblk2 = pja_hmatr_dbl[j];
                        jhblk2 = pja2_hmatr_dbl[j];
                        if (jblk1 != jblk2 || jhblk1 != jhblk2) {
                           cout <<
                              " CK3D_SolverThreads <>::UseFastTransform: error 2: incorrect block sparsity !!! "
                              << endl;
                           throw
                              " CK3D_SolverThreads <>::UseFastTransform: error 2: incorrect block sparsity !!! ";
                        }
                     }
                     int k, iend, iend_dbl, ip, ip_dbl, jj, jj_dbl;
                     int ihblk_ord, indblk_ord, ind_elem_ord;
                     SIndFT *pind = NULL;
                     CMatrix < _Int, _Flt > *pASub_ord = NULL;
                     _Flt *pa_ord = NULL;
                     typename std::map < int, vector < pair < int, pair < int,
                        _Flt > > > >::iterator it;
                     for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
                        int nlist_loc = pASub[j].GetNlist ();
                        _Int *plist_loc = pASub[j].GetListArr ();
                        _Int *pia_loc = pASub[j].GetIaArr ();
                        _Int *pja_loc = pASub[j].GetJaArr ();
                        _Flt *pa_loc = pASub[j].GetAArr ();
                        int nlist_dbl = pASub_dbl[j].GetNlist ();
                        _Int *plist_dbl = pASub_dbl[j].GetListArr ();
                        _Int *pia_dbl = pASub_dbl[j].GetIaArr ();
                        _Int *pja_dbl = pASub_dbl[j].GetJaArr ();
                        double *pa_dbl = pASub_dbl[j].GetAArr ();
                        if (nlist_loc != nlist_dbl) {
                           cout <<
                              " CK3D_SolverThreads <>::UseFastTransform: error 3: incorrect sparsity !!! "
                              << endl;
                           throw
                              " CK3D_SolverThreads <>::UseFastTransform: error 3: incorrect sparsity !!! ";
                        }
                        for (k = 0; k < nlist_loc; k++) {
                           if (plist_loc[k] != plist_dbl[k]) {
                              cout <<
                                 " CK3D_SolverThreads <>::UseFastTransform: error 4: incorrect sparsity !!! "
                                 << endl;
                              throw
                                 " CK3D_SolverThreads <>::UseFastTransform: error 4: incorrect sparsity !!! ";
                           }
                        }
                        for (k = 0; k < nlist_loc; k++) {
                           iend = (int) pia_loc[k + 1] - 1;
                           iend_dbl = (int) pia_dbl[k + 1] - 1;
                           ip = (int) pia_loc[k];
                           ip_dbl = (int) pia_dbl[k];
                           while (ip <= iend || ip_dbl <= iend_dbl) {
                              if (ip <= iend && ip_dbl <= iend_dbl) {
                                 jj = (int) pja_loc[ip];
                                 jj_dbl = (int) pja_dbl[ip_dbl];
                                 if (jj == jj_dbl) {
                                    pind = (SIndFT *) (pa_dbl + ip_dbl);
                                    ihblk_ord = (int) pind->ihblk_FT;
                                    indblk_ord = (int) pind->ind_blk_FT;
                                    ind_elem_ord = (int) pind->ind_elem_FT;
                                    if (phblk2cpu[ihblk_ord] == myid) {
                                       pASub_ord = phmatr_ord[ihblk_ord].GetASubArr ();
                                       pa_ord = pASub_ord[indblk_ord].GetAArr ();
                                       pa_ord[ind_elem_ord] = pa_loc[ip];
                                    } else {
                                       it = pmap_hblk2st_threads[my_thr].find (ihblk_ord);
                                       if (it == pmap_hblk2st_threads[my_thr].end ()) {
                                          vector < pair < int, pair < int,
                                             _Flt > > >vect_pair_new;
                                          pmap_hblk2st_threads[my_thr][ihblk_ord] =
                                             vect_pair_new;
                                          it =
                                             pmap_hblk2st_threads[my_thr].
                                             find (ihblk_ord);
                                       }
                                       pair < int, _Flt > pair_flt;
                                       pair_flt.first = ind_elem_ord;
                                       pair_flt.second = pa_loc[ip];
                                       pair < int, pair < int, _Flt > >pair_pair_flt;
                                       pair_pair_flt.first = indblk_ord;
                                       pair_pair_flt.second = pair_flt;
                                       it->second.push_back (pair_pair_flt);
                                    }
                                    ip++;
                                    ip_dbl++;
                                 } else if (jj > jj_dbl) {
                                    ip_dbl++;
                                 } else {
                                    cout <<
                                       " CK3D_SolverThreads <>::UseFastTransform: error 5: incorrect sparsity !!! "
                                       << endl;
                                    throw
                                       " CK3D_SolverThreads <>::UseFastTransform: error 5: incorrect sparsity !!! ";
                                 }
                              } else if (ip_dbl <= iend_dbl) {
                                 ip_dbl++;
                              } else if (ip <= iend) {
                                 cout <<
                                    " CK3D_SolverThreads <>::UseFastTransform: error 6: incorrect sparsity !!! "
                                    << endl;
                                 throw
                                    " CK3D_SolverThreads <>::UseFastTransform: error 6: incorrect sparsity !!! ";
                              }
                           }
                        }
                     }
                  }
               }

            }
         }

      }

// Move each separate hblock into usual block to be sent

      int nlist_hblk = 0;

      {
         typename std::map < int, vector < pair < int, pair < int,
            _Flt > > > >::iterator it;
         typename std::map < int, vector < pair < int, pair < int,
            _Flt > > > >::iterator it_end;
         int i;
         for (i = 0; i < n_thr; i++) {
            it = pmap_hblk2st_threads[i].begin ();
            it_end = pmap_hblk2st_threads[i].end ();
            while (it != it_end) {
               nlist_hblk++;
               it++;
            }
         }
      }

      CVectorData < int >list_hblk (nlist_hblk);
      CVectorData < vector < pair < int, pair < int, _Flt > > >*>ptr_vect (nlist_hblk);

      int *plist_hblk = list_hblk.Ptr ();
      vector < pair < int, pair < int, _Flt > > >**pptr_vect = ptr_vect.Ptr ();

      nlist_hblk = 0;

      {
         typename std::map < int, vector < pair < int, pair < int,
            _Flt > > > >::iterator it;
         typename std::map < int, vector < pair < int, pair < int,
            _Flt > > > >::iterator it_end;
         int i;
         for (i = 0; i < n_thr; i++) {
            it = pmap_hblk2st_threads[i].begin ();
            it_end = pmap_hblk2st_threads[i].end ();
            while (it != it_end) {
               plist_hblk[nlist_hblk] = it->first;
               pptr_vect[nlist_hblk] = &(it->second);
               nlist_hblk++;
               it++;
            }
         }
      }

      CVectorData < CMatrix < int, _Flt > >blks_send (nlist_hblk);
      CMatrix < int, _Flt > *pblks_send = blks_send.Ptr ();

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nlist_hblk; ipar++) {
            int nz = (int) pptr_vect[ipar]->size ();
            pblks_send[ipar].ResizeAndSetAll (0, nz, 0, nz, nz);
            int *plist_send = pblks_send[ipar].GetList2Arr ();
            int *pja_send = pblks_send[ipar].GetJa2Arr ();
            _Flt *pa_send = pblks_send[ipar].GetAArr ();
            int j;
            pair < int, pair < int, _Flt > >pair_pair_flt;
            pair < int, _Flt > pair_flt;
            for (j = 0; j < nz; j++) {
               pair_pair_flt = (*pptr_vect[ipar])[j];
               plist_send[j] = pair_pair_flt.first;
               pair_flt = pair_pair_flt.second;
               pja_send[j] = pair_flt.first;
               pa_send[j] = pair_flt.second;
            }
            {
               vector < pair < int, pair < int, _Flt > > >vect_dummy;
               pptr_vect[ipar]->swap (vect_dummy);
            }
         }
      }

// Prepare sends

      vector < int >CpuIDSend;
      vector < vector < char > >ObjSend;

      {

         std::map < int, pair < int, int > >map_cpu2nz;
         std::map < int, pair < int, int > >::iterator it;

         int i, ihblk, icpu;
         pair < int, int >*ppair_int;

         for (i = 0; i < nlist_hblk; i++) {
            ihblk = plist_hblk[i];
            icpu = phblk2cpu[ihblk];
            it = map_cpu2nz.find (icpu);
            if (it == map_cpu2nz.end ()) {
               map_cpu2nz.insert (std::pair < int, pair < int,
                                  int > >(icpu, pair < int, int >(-1, 0)));
               it = map_cpu2nz.find (icpu);
            }
            ppair_int = &(it->second);
            ppair_int->second++;
         }

         int nproc_send = 0;

         it = map_cpu2nz.begin ();
         while (it != map_cpu2nz.end ()) {
            ppair_int = &(it->second);
            ppair_int->first = nproc_send;
            nproc_send++;
            it++;
         }

         CpuIDSend.resize (nproc_send);
         ObjSend.resize (nproc_send);

         int *pCpuIDSend = NULL;
         vector < char >*pObjSend = NULL;

         if (nproc_send > 0) {
            pCpuIDSend = &CpuIDSend[0];
            pObjSend = &ObjSend[0];
         }

         CVectorData < CBMatrix < int, _Flt > >hblks_send (nproc_send);
         CBMatrix < int, _Flt > *phblks_send = hblks_send.Ptr ();

         nproc_send = 0;

         it = map_cpu2nz.begin ();
         while (it != map_cpu2nz.end ()) {
            pCpuIDSend[nproc_send] = it->first;
            ppair_int = &(it->second);
            int nzblk = ppair_int->second;
            phblks_send[nproc_send].SetNzblk (nzblk);
            phblks_send[nproc_send].ResizeASub (nzblk);
            CMatrix < int, float >*pHMatrStr = phblks_send[nproc_send].GetHMatrStr ();
            pHMatrStr->ResizeAndSetAllSp (0, nzblk, 0, 0);
            ppair_int->second = 0;
            nproc_send++;
            it++;
         }

         int ind, k;

         for (i = 0; i < nlist_hblk; i++) {
            ihblk = plist_hblk[i];
            icpu = phblk2cpu[ihblk];
            it = map_cpu2nz.find (icpu);
            ppair_int = &(it->second);
            ind = ppair_int->first;
            k = ppair_int->second;
            CMatrix < int, _Flt > *pASub = phblks_send[ind].GetASubArr ();
            CMatrix < int, float >*pHMatrStr = phblks_send[ind].GetHMatrStr ();
            int *plist2_hmatr = pHMatrStr->GetList2Arr ();
            plist2_hmatr[k] = ihblk;
            pASub[k].ReplaceFree (pblks_send[i]);
            ppair_int->second++;
         }

         long long isize;
         char *pobj;

         for (i = 0; i < nproc_send; i++) {
            isize = phblks_send[i].GetPackedSize ();
            pObjSend[i].resize ((size_t) isize);
            pobj = &(pObjSend[i][0]);
            phblks_send[i].FillPacked_thr (isize, pobj);
            phblks_send[i].Clean ();
         }

      }

      vector < int >CpuIDRecv;
      vector < vector < char > >ObjRecv;

      CMPIDataExchange::DataExchange (pcommloc, CpuIDSend, ObjSend, CpuIDRecv, ObjRecv);

// Free send data

      {
         vector < int >CpuIDSend_temp;
         vector < vector < char > >ObjSend_temp;
         CpuIDSend.swap (CpuIDSend_temp);
         ObjSend.swap (ObjSend_temp);
      }

// Unpack receive data

      int nrecv = (int) CpuIDRecv.size ();

      vector < char >*pObjRecv = NULL;

      if (nrecv > 0) {
         pObjRecv = &ObjRecv[0];
      }

      vector < CBMatrix < int, _Flt > >hblk_recv (nrecv + 1);

      CBMatrix < int, _Flt > *phblk_recv = &hblk_recv[0];

      {
         int i;
         long long isize;
         char *pobj;
         for (i = 0; i < nrecv; i++) {
            isize = (long long) pObjRecv[i].size ();
            pobj = &(pObjRecv[i][0]);
            phblk_recv[i].UnPack (isize, pobj);
         }
      }

// Free recv data

      {
         vector < vector < char > >ObjRecv_temp;
         ObjRecv.swap (ObjRecv_temp);
      }

// Use received data

      {

         int i;

         int nlist_hblk_recv = 0;

         for (i = 0; i < nrecv; i++) {
            nlist_hblk_recv += phblk_recv[i].GetNzblk ();
         }

         CVectorData < int >list_hblk_recv (nlist_hblk_recv);
         CVectorData < CMatrix < int, _Flt > *>ptr_blks_recv (nlist_hblk_recv);

         int *plist_hblk_recv = list_hblk_recv.Ptr ();
         CMatrix < int, _Flt > **pptr_blks_recv = ptr_blks_recv.Ptr ();

         nlist_hblk_recv = 0;

         int j;

         for (i = 0; i < nrecv; i++) {
            int nzblk = phblk_recv[i].GetNzblk ();
            CMatrix < int, _Flt > *pASub = phblk_recv[i].GetASubArr ();
            CMatrix < int, float >*pHMatrStr = phblk_recv[i].GetHMatrStr ();
            int *plist2_hmatr = pHMatrStr->GetList2Arr ();
            for (j = 0; j < nzblk; j++) {
               plist_hblk_recv[nlist_hblk_recv] = plist2_hmatr[j];
               pptr_blks_recv[nlist_hblk_recv] = pASub + j;
               nlist_hblk_recv++;
            }
         }

         {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < nlist_hblk_recv; ipar++) {
               int ihblk = plist_hblk_recv[ipar];
               if (ihblk < 0 || ihblk >= nhblks_loc) {
                  cout <<
                     " CK3D_SolverThreads <>::UseFastTransform: Error: ihblk is out of range !!! "
                     << endl;
                  throw
                     " CK3D_SolverThreads <>::UseFastTransform: Error: ihblk is out of range !!! ";
               }
               if (phblk2cpu[ihblk] != myid) {
                  cout <<
                     " CK3D_SolverThreads <>::UseFastTransform: Error: ihblk is not own !!! "
                     << endl;
                  throw
                     " CK3D_SolverThreads <>::UseFastTransform: Error: ihblk is not own !!! ";
               }
               int nzja = pptr_blks_recv[ipar]->GetNlist2 ();
               int *plist_recv = pptr_blks_recv[ipar]->GetList2Arr ();
               int *pja_recv = pptr_blks_recv[ipar]->GetJa2Arr ();
               _Flt *pa_recv = pptr_blks_recv[ipar]->GetAArr ();
               int nzblk_loc = phmatr_ord[ihblk].GetNzblk ();
               CMatrix < _Int, _Flt > *pASub_ord = phmatr_ord[ihblk].GetASubArr ();
               int ind_blk, ind_elem;
               _Flt *pa_ord = NULL;
               int j, nza_loc;
               for (j = 0; j < nzja; j++) {
                  ind_blk = plist_recv[j];
                  ind_elem = pja_recv[j];
                  if (ind_blk < 0 || ind_blk >= nzblk_loc) {
                     cout <<
                        " CK3D_SolverThreads <>::UseFastTransform Error: ind_blk is out of range !!! "
                        << endl;
                     throw
                        " CK3D_SolverThreads <>::UseFastTransform Error: ind_blk is out of range !!! ";
                  }
                  nza_loc = pASub_ord[ind_blk].GetNza ();
                  if (ind_elem < 0 || ind_elem >= nza_loc) {
                     cout <<
                        " CK3D_SolverThreads <>::UseFastTransform Error: ind_elem is out of range !!! "
                        << endl;
                     throw
                        " CK3D_SolverThreads <>::UseFastTransform Error: ind_elem is out of range !!! ";
                  }
                  pa_ord = pASub_ord[ind_blk].GetAArr ();
                  pa_ord[ind_elem] = pa_recv[j];
               }
            }
         }

      }

// Replace computed data

      this->hmatr_arr.swap (hmatr_ord);

   }

// Prepare solver structures including performing parallel fct
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::ComputeBILU2 (SParams * _params, SStatData * _stats)
   {

// Init output data

      _stats->prec_extend = 1.0;
      _stats->density = 0.0e0;
      _stats->scpiv_min = 1.0e100;
      _stats->scpiv_max = -1.0e100;
      _stats->nmodif = 0;
      _stats->piv_min = 1.0e100;
      _stats->piv_max = -1.0e100;
      _stats->dtime_fct = 0.0e0;

// Init timer

      void *pcommloc = this->pcomm;

      CMPIDataExchange::Synchronize (pcommloc);

      double time0;

      time0 = CMPIDataExchange::GetWallTimeMPI ();

// Store control data

      int nhblksloc = this->nhblks;
      int nblksloc = this->nblks;

      this->params = *_params;

      long long *phblks = &this->hblks[0];
      int *phblk2cpu = &this->hblk2cpu[0];
      int *phblk2blks = &this->hblk2blks[0];
      int *pblk2hblks = &this->blk2hblks[0];
      long long *pblks = &this->blks[0];

      int myid = CMPIDataExchange::GetMyid (pcommloc);

      _params->iparam1 = myid;

// Perform temporary explicit block scaling

      int sctype_save = _params->sctype;
      double scpiv_min_save = 1.0e100;
      double scpiv_max_save = -1.0e100;

      if (_params->sctype == 2) {
         this->b_use_blksize = true;
         this->ExplicitBlockScaling (_params, _stats);
         _params->sctype = -1;
         scpiv_min_save = _stats->scpiv_min;
         scpiv_max_save = _stats->scpiv_max;
      }
//      char strbuff_debug[256];
//      sprintf (strbuff_debug,"ChkCompBIlu_1x1_%i.dat",myid);
//      ofstream ffout_debug (strbuff_debug);

//   ffout_debug << " Begin check !!! " << endl;

// Allocate arrays for data

//   ffout_debug << " Point 1  " << endl;

      CBMatrix < _Int, _Flt > *phmatr_arr = &this->hmatr_arr[0];

      this->hblks_ext.resize (nhblksloc + 1);
      this->nlist_ext_arr.resize (nhblksloc + 1);
      this->list_ext_arr.resize (nhblksloc + 1);
      this->nblks_ext_arr.resize (nhblksloc + 1);
      this->blksnum_ext_arr.resize (nhblksloc + 1);
      this->blks_ext_arr.resize (nhblksloc + 1);
      this->tree_arr.resize (3 * nhblksloc + 1);
      this->nblks_ilu2_arr.resize (nhblksloc + 1);
      this->nblks1_ilu2_arr.resize (nhblksloc + 1);
      this->nblks2_ilu2_arr.resize (nhblksloc + 1);
      this->blks_ilu2_arr.resize (nhblksloc + 1);
      this->nzord_ilu2_arr.resize (nhblksloc + 1);
      this->order_LU.resize (nhblksloc + 1);
      this->matrL_float.resize (nhblksloc + 1);
      this->matrL_double.resize (nhblksloc + 1);
      this->matrU_float.resize (nhblksloc + 1);
      this->matrU_double.resize (nhblksloc + 1);

      long long *phblks_ext = &this->hblks_ext[0];
      int *pnlist_ext_arr = &this->nlist_ext_arr[0];
      vector < int >*plist_ext_arr = &this->list_ext_arr[0];
      int *pnblks_ext_arr = &this->nblks_ext_arr[0];
      vector < int >*pblksnum_ext_arr = &this->blksnum_ext_arr[0];
      vector < long long >*pblks_ext_arr = &this->blks_ext_arr[0];
      CTree *ptree_arr = &this->tree_arr[0];
      int *pnblks_ilu2_arr = &this->nblks_ilu2_arr[0];
      int *pnblks1_ilu2_arr = &this->nblks1_ilu2_arr[0];
      int *pnblks2_ilu2_arr = &this->nblks2_ilu2_arr[0];
      vector < long long >*pblks_ilu2_arr = &this->blks_ilu2_arr[0];
      vector < long long >*pnzord_ilu2_arr = &this->nzord_ilu2_arr[0];
      CVectorData < int >*porder_LU = &this->order_LU[0];
      CBMatrix < _Int, float >*pmatrL_float = &this->matrL_float[0];
      CBMatrix < _Int, double >*pmatrL_double = &this->matrL_double[0];
      CBMatrix < _Int, float >*pmatrU_float = &this->matrU_float[0];
      CBMatrix < _Int, double >*pmatrU_double = &this->matrU_double[0];

// Get nzatot

      long long nza_tot = 0;

      int i;

      for (i = 0; i < nhblksloc; i++) {
         if (phblk2cpu[i] == myid) {
            nza_tot += phmatr_arr[i].GetNzatot ();
         }
      }

      CMPIDataExchange::ExchangeArray (pcommloc, 'L', '+', 1, (void *) (&nza_tot));

// Compute extended lists

//   ffout_debug << " Point 2  " << endl;

      int ncycle_loc = _params->ncycle;

      CBMatrix < _Int, _Flt >::ExtendedLists (pcommloc, ncycle_loc, nhblksloc, phblk2cpu,
                                              phblk2blks, pblk2hblks, nblksloc, pblks,
                                              phmatr_arr, pnlist_ext_arr, plist_ext_arr);

// Perform filtering of the lists

      CBMatrix < _Int, _Flt >::FilterListsBack (myid, nhblksloc, phblk2cpu, pblk2hblks,
                                                pnlist_ext_arr, plist_ext_arr);

//   ffout_debug << " Point 3  " << endl;

// Compute extended hblocks array

      for (i = 0; i <= nhblksloc; i++)
         phblks_ext[i] = 0;

      for (i = 0; i < nhblksloc; i++) {
         if (phblk2cpu[i] == myid)
            phblks_ext[i + 1] = pnlist_ext_arr[i] + (phblks[i + 1] - phblks[i]);
      }

      CMPIDataExchange::ExchangeArray (pcommloc, 'L', '+', nhblksloc + 1,
                                       (void *) phblks_ext);

      for (i = 0; i < nhblksloc; i++)
         phblks_ext[i + 1] = phblks_ext[i] + phblks_ext[i + 1];

// Create separated initial extended blocks partitionings

      for (i = 0; i < nhblksloc; i++)
         pnblks_ext_arr[i] = 0;

      int j, i2, i2_prev;

//   ffout_debug << " Point 4  " << endl;

      for (i = 0; i < nhblksloc; i++) {
         if (phblk2cpu[i] == myid) {
            int *plist_ext_temp = NULL;
            if (pnlist_ext_arr[i] > 0)
               plist_ext_temp = &plist_ext_arr[i][0];
            int nblks_temp = 0;
            i2_prev = -1;
            for (j = 0; j < pnlist_ext_arr[i]; j++) {
               i2 = plist_ext_temp[j * 2 + 1];
               if (i2 != i2_prev) {
                  nblks_temp++;
                  i2_prev = i2;
               }
            }
            nblks_temp += (phblk2blks[i + 1] - phblk2blks[i]);
            pnblks_ext_arr[i] = nblks_temp;
            pblks_ext_arr[i].resize (nblks_temp + 2);
            pblksnum_ext_arr[i].resize (nblks_temp + 2);
            long long *pblks_ext_arr_temp = &pblks_ext_arr[i][0];
            int *pblksnum_ext_arr_temp = &pblksnum_ext_arr[i][0];
            long long nz_temp = 0;
            i2_prev = -1;
            nblks_temp = 0;
            pblks_ext_arr_temp[0] = 0;
            for (j = 0; j < pnlist_ext_arr[i]; j++) {
               i2 = plist_ext_temp[j * 2 + 1];
               if (i2 != i2_prev) {
                  pblksnum_ext_arr_temp[nblks_temp] = i2;
                  nblks_temp++;
                  i2_prev = i2;
               }
               nz_temp++;
               pblks_ext_arr_temp[nblks_temp] = nz_temp;
            }
            for (j = phblk2blks[i]; j < phblk2blks[i + 1]; j++) {
               pblksnum_ext_arr_temp[nblks_temp] = j;
               pblks_ext_arr_temp[nblks_temp + 1] =
                  pblks_ext_arr_temp[nblks_temp] + (pblks[j + 1] - pblks[j]);
               nblks_temp++;
            }
         }
      }

// Get extended submatrices

//   ffout_debug << " Point 5  " << endl;

      vector < CBMatrix < _Int, _Flt > >hmatr_ext_arr (nhblksloc);
      CBMatrix < _Int, _Flt > *phmatr_ext_arr = &hmatr_ext_arr[0];

      CBMatrix < _Int, _Flt >::GetExtendedSubmatrices (pcommloc, nhblksloc, phblk2cpu,
                                                       phblk2blks, pblk2hblks, nblksloc,
                                                       pblks, phmatr_arr, pnblks_ext_arr,
                                                       pblksnum_ext_arr, pblks_ext_arr,
                                                       pnlist_ext_arr, plist_ext_arr,
                                                       phmatr_ext_arr);

      long long ntot = pblks[nblksloc];
      long long ntot_ext = 0;

      for (i = 0; i < nhblksloc; i++) {
         if (phblk2cpu[i] == myid) {
            ntot_ext += (phblks[i + 1] - phblks[i]);
            ntot_ext += pnlist_ext_arr[i];
         }
      }

      CMPIDataExchange::ExchangeArray (pcommloc, 'L', '+', 1, (void *) (&ntot_ext));

      _stats->prec_extend = ((double) (ntot_ext) / (double) (ntot));

      int collap_loc = _params->collap;

      char strbuff[256];

      int blks_2[3];

      if (collap_loc > 0) {
         for (i = 0; i < nhblksloc; i++) {
            if (phblk2cpu[i] == myid) {
               CMatrix < int, float >*phmatr = phmatr_ext_arr[i].GetHMatrStr ();
               int nblks_temp = (int) (phblk2blks[i + 1] - phblk2blks[i]);
               sprintf (strbuff, "HBlkStrExt_%i.ps", i);
               blks_2[0] = 0;
               blks_2[1] = pnblks_ext_arr[i] - nblks_temp;
               blks_2[2] = pnblks_ext_arr[i];
               CBMatrix < int, float >::Str2PsBox (*phmatr, strbuff, 2, blks_2);
            }
         }
      }
// Compute new ordering of submatrices and reorder them

      vector < CVectorData < int > >order_arr (nhblksloc + 1);

      CVectorData < int >*porder_arr = &order_arr[0];

//   ffout_debug << " Point 6  " << endl;

#ifdef USE_THREADS
      //n_thr = omp_get_max_threads ();
#endif

      int n_thr_part = 1;

#ifdef USE_THREADS
      n_thr_part = omp_get_max_threads ();
#endif

      int nparts_pwr2 = 1;

      while (nparts_pwr2 < n_thr_part)
         nparts_pwr2 *= 2;

      bool b_blk_wells_temp = false;

      for (i = 0; i < nhblksloc; i++) {
         if (phblk2cpu[i] == myid) {

            int nblks_temp = phblk2blks[i + 1] - phblk2blks[i];
            long long *ppblks_ext_arr = &pblks_ext_arr[i][0];

            b_blk_wells_temp = false;
            if (i == nhblksloc - 1 && this->b_blk_wells)
               b_blk_wells_temp = true;

//            if (i == nhblksloc-1 && b_blk_wells_temp) {
//               ofstream fffout ("ChkHblkExt.dat");
//               phmatr_ext_arr[i].PrintHMatrix (fffout);
//            }

            int nparts_pwr2_temp = nparts_pwr2;
            if (b_blk_wells_temp)
               nparts_pwr2_temp = 1;

//            if (i == 0) {
//               ffout_debug << " Ihblk = " << i << " Before order " << endl;
//               phmatr_ext_arr[i].PrintHMatrix (ffout_debug);
//            }

            CBMatrix < _Int, _Flt >::ComputeOptimalOrderSchurWells (b_blk_wells_temp,
                                                                    nparts_pwr2_temp,
                                                                    pnblks_ext_arr[i],
                                                                    pnblks_ext_arr[i] -
                                                                    nblks_temp,
                                                                    ppblks_ext_arr,
                                                                    phmatr_ext_arr + i,
                                                                    ptree_arr[i * 3],
                                                                    ptree_arr[i * 3 + 1],
                                                                    pnblks_ilu2_arr[i],
                                                                    pnblks1_ilu2_arr[i],
                                                                    pnblks2_ilu2_arr[i],
                                                                    pblks_ilu2_arr[i],
                                                                    pnzord_ilu2_arr[i],
                                                                    porder_arr[i]);

            long long *ppblks_new_arr = &pblks_ilu2_arr[i][0];
            int *pporder_arr = porder_arr[i].Ptr ();

            CBMatrix < _Int, _Flt > hmatr_ord;

            CBMatrix < _Int, _Flt >::ReorderHMatrix (pnblks_ext_arr[i], ppblks_ext_arr,
                                                     phmatr_ext_arr[i], pporder_arr,
                                                     pnblks_ilu2_arr[i], ppblks_new_arr,
                                                     hmatr_ord);

//            if (i == nhblksloc-1 && b_blk_wells_temp) {
//               ofstream fffout ("ChkHblkExtOrd.dat");
//               hmatr_ord.PrintHMatrix (fffout);
//            }

            phmatr_ext_arr[i].ReplaceFree (hmatr_ord);

//            if (i == 0) {
//               ffout_debug << " Reordered hblock " << endl;
//               phmatr_ext_arr[i].PrintHMatrix (ffout_debug);
//            }

            if (collap_loc > 0) {
               sprintf (strbuff, "HBlkStrExtOrd_%i.ps", i);
               int blks_3[4];
               blks_3[0] = 0;
               blks_3[1] = pnblks1_ilu2_arr[i];
               blks_3[2] = pnblks1_ilu2_arr[i] + pnblks2_ilu2_arr[i];
               blks_3[3] = pnblks_ilu2_arr[i];
               CMatrix < int, float >*phmatr = phmatr_ext_arr[i].GetHMatrStr ();
               CBMatrix < int, float >::Str2PsBox (*phmatr, strbuff, 3, blks_3);
            }

         }
      }

// Perform computation of the ILU decomposition for all local ordered blocks

//   ffout_debug << " Point 7  " << endl;

      int prec_float_loc = _params->prec_float;

      vector < double >sclmin_arr (nhblksloc + 1);
      vector < double >sclmax_arr (nhblksloc + 1);
      vector < double >eigmin_arr (nhblksloc + 1);
      vector < double >eigmax_arr (nhblksloc + 1);

      double *psclmin_arr = &sclmin_arr[0];
      double *psclmax_arr = &sclmax_arr[0];
      double *peigmin_arr = &eigmin_arr[0];
      double *peigmax_arr = &eigmax_arr[0];

      for (i = 0; i < nhblksloc; i++)
         psclmin_arr[i] = 0;
      for (i = 0; i < nhblksloc; i++)
         psclmax_arr[i] = 0;
      for (i = 0; i < nhblksloc; i++)
         peigmin_arr[i] = 0;
      for (i = 0; i < nhblksloc; i++)
         peigmax_arr[i] = 0;

      long long nzlu_tot = 0;

      _stats->nmodif = 0;

      int nmodif_loc = 0;

      for (i = 0; i < nhblksloc; i++) {
         if (phblk2cpu[i] == myid) {

            int nblks_ilu2 = 0;
            vector < long long >blks_ilu2 (1);
            CVectorData < int >order_ilu2 (1);

            if (prec_float_loc == 1) {

               CBMatrix < _Int, float >*ptr_matr = NULL;
               CBMatrix < _Int, float >hmatr_conv;
               if (sizeof (_Flt) == sizeof (float)) {
                  ptr_matr = (CBMatrix < _Int, float >*) (phmatr_ext_arr + i);
               } else {
                  CMatrixConv < _Int, _Flt, float >::InitAndConv (phmatr_ext_arr[i],
                                                                  hmatr_conv);
                  phmatr_ext_arr[i].Clean ();
                  ptr_matr = &hmatr_conv;
               }

               long long *ppblks_new_arr = &pblks_ilu2_arr[i][0];
               long long *ppnzord_ilu2_arr = NULL;
//            long long *ppnzord_ilu2_arr = &pnzord_ilu2_arr[i][0];

//               _params->iparam1 = myid;
               _params->iparam1 = i;
               _params->iparam2 = 0;

               b_blk_wells_temp = false;
               if (i == nhblksloc - 1 && this->b_blk_wells)
                  b_blk_wells_temp = true;

               if (_params->msglev >= 2) {
                  cout << " Fct Ihblk: " << i << " Nblks = " << pnblks_ilu2_arr[i] -
                     pnblks1_ilu2_arr[i] << " Nblks_ext = " << pnblks_ilu2_arr[i] <<
                     " Ni_Loc = " << ppblks_new_arr[pnblks_ilu2_arr[i]] -
                     ppblks_new_arr[pnblks1_ilu2_arr[i]] << " Ni_Ext = " <<
                     ppblks_new_arr[pnblks1_ilu2_arr[i]] << endl;
               }

               CFctThreads < _Int, float >::Ilu2HMatrix (b_blk_wells_temp, *_params,
                                                         pnblks_ilu2_arr[i],
                                                         pnblks1_ilu2_arr[i],
                                                         pnblks2_ilu2_arr[i],
                                                         ppblks_new_arr, ppnzord_ilu2_arr,
                                                         ptree_arr[i * 3],
                                                         ptree_arr[i * 3 + 1], *ptr_matr,
                                                         pmatrL_float[i], pmatrU_float[i],
                                                         ptree_arr[i * 3 + 2], nblks_ilu2,
                                                         blks_ilu2, order_ilu2,
                                                         psclmin_arr[i], psclmax_arr[i],
                                                         nmodif_loc, peigmin_arr[i],
                                                         peigmax_arr[i]);

               _stats->nmodif += nmodif_loc;
               nzlu_tot += pmatrL_float[i].GetNzatot ();
               nzlu_tot += pmatrU_float[i].GetNzatot ();

            } else {

               CBMatrix < _Int, double >*ptr_matr = NULL;
               CBMatrix < _Int, double >hmatr_conv;
               if (sizeof (_Flt) == sizeof (double)) {
                  ptr_matr = (CBMatrix < _Int, double >*) (phmatr_ext_arr + i);
               } else {
                  CMatrixConv < _Int, _Flt, double >::InitAndConv (phmatr_ext_arr[i],
                                                                   hmatr_conv);
                  phmatr_ext_arr[i].Clean ();
                  ptr_matr = &hmatr_conv;
               }

               long long *ppblks_new_arr = &pblks_ilu2_arr[i][0];
               long long *ppnzord_ilu2_arr = NULL;
//            long long *ppnzord_ilu2_arr = &pnzord_ilu2_arr[i][0];

//               _params->iparam1 = myid;
               _params->iparam1 = i;
               _params->iparam2 = 0;

               b_blk_wells_temp = false;
               if (i == nhblksloc - 1 && this->b_blk_wells)
                  b_blk_wells_temp = true;

               if (_params->msglev >= 2) {
                  cout << " Fct Ihblk: " << i << " Nblks = " << pnblks_ilu2_arr[i] -
                     pnblks1_ilu2_arr[i] << " Nblks_ext = " << pnblks_ilu2_arr[i] <<
                     " Ni_Loc = " << ppblks_new_arr[pnblks_ilu2_arr[i]] -
                     ppblks_new_arr[pnblks1_ilu2_arr[i]] << " Ni_Ext = " <<
                     ppblks_new_arr[pnblks1_ilu2_arr[i]] << endl;
               }
//               sprintf (strbuff_debug,"ChkBIlu_1x1_%i.dat",i);
//               ofstream ffout_debug1 (strbuff_debug);

//               if (false) {
//                  ffout_debug1 << "  Ihblk = " << i << endl;
//                  ffout_debug1 << " nblks = " << pnblks_ilu2_arr[i] << " nblks1 = " <<  pnblks1_ilu2_arr[i] << " nblks2 = " <<  pnblks2_ilu2_arr[i] << endl;
//                  PrintArray (ffout_debug1, " blks ", pnblks_ilu2_arr[i]+1, ppblks_new_arr);
//                  ffout_debug1 << " Tree 1" << endl;
//                  ptree_arr[i * 3].OutputTree (ffout_debug1);
//                  ffout_debug1 << " Tree 2" << endl;
//                  ptree_arr[i * 3+1].OutputTree (ffout_debug1);
//                  ffout_debug1 << " HMatr " << endl;
//                  ptr_matr->PrintHMatrix (ffout_debug1);
//               }

               CFctThreads < _Int, double >::Ilu2HMatrix (b_blk_wells_temp, *_params,
                                                          pnblks_ilu2_arr[i],
                                                          pnblks1_ilu2_arr[i],
                                                          pnblks2_ilu2_arr[i],
                                                          ppblks_new_arr,
                                                          ppnzord_ilu2_arr,
                                                          ptree_arr[i * 3],
                                                          ptree_arr[i * 3 + 1], *ptr_matr,
                                                          pmatrL_double[i],
                                                          pmatrU_double[i],
                                                          ptree_arr[i * 3 + 2],
                                                          nblks_ilu2, blks_ilu2,
                                                          order_ilu2, psclmin_arr[i],
                                                          psclmax_arr[i], nmodif_loc,
                                                          peigmin_arr[i], peigmax_arr[i]);

//               if (false) {
//                  ffout_debug1 << " Tree 3" << endl;
//                  ptree_arr[i * 3+2].OutputTree (ffout_debug1);
//                  ffout_debug1 << " L = " << endl;
//                  pmatrL_double[i].PrintHMatrix (ffout_debug1);
//                  ffout_debug1 << " U = " << endl;
//                  pmatrU_double[i].PrintHMatrix (ffout_debug1);
//               }

               _stats->nmodif += nmodif_loc;
               nzlu_tot += pmatrL_double[i].GetNzatot ();
               nzlu_tot += pmatrU_double[i].GetNzatot ();

//               if (i == nhblksloc-1) {
//                  char strbuff[256];
//                  sprintf (strbuff,"PrecL_%i.dat",i);
//                  ofstream fffout (strbuff);
//                  fffout << " Nmodif = " << nmodif_loc << " Eigmin = " << peigmin_arr[i] << " Eigmax = " << peigmax_arr[i] << endl;
//                  pmatrL_double[i].PrintHMatrix (fffout);
//                  sprintf (strbuff,"PrecU_%i.dat",i);
//                  ofstream ffffout (strbuff);
//                  pmatrU_double[i].PrintHMatrix (ffffout);
//               }

            }

            if (_params->msglev >= 2) {
               cout << " Fct End Ihblk =  " << i << " PivMin = " << peigmin_arr[i] <<
                  " PivMax = " << peigmax_arr[i] << " NModif = " << nmodif_loc << endl;
            }
// Replace partitioning and ordering

            pnblks_ilu2_arr[i] = nblks_ilu2;
            pblks_ilu2_arr[i].swap (blks_ilu2);

            {

// Compute inverse initial order

               long long *ppblks_ilu2_arr = &pblks_ilu2_arr[i][0];

               int ntot_curr = (int) ppblks_ilu2_arr[nblks_ilu2];

               porder_LU[i].resize (ntot_curr);
               int *pporder_LU = porder_LU[i].Ptr ();

               int *pporder_arr = porder_arr[i].Ptr ();
               int *porder_ilu2 = order_ilu2.Ptr ();

               CVectorData < int >iorder_ini (ntot_curr);
               int *piorder_ini = iorder_ini.Ptr ();

               CVectorInt < int >::InvOrder_thr (ntot_curr, pporder_arr, piorder_ini);

               {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < nblks_ilu2; ipar++) {
                     for (int j = (int)ppblks_ilu2_arr[ipar];
                          j < ppblks_ilu2_arr[ipar + 1]; j++) {
                        int jold = piorder_ini[j];
                        pporder_LU[jold] = porder_ilu2[j];
                     }
                  }
               }

            }

// Print block sparsity if necessary

            if (collap_loc > 0) {
               sprintf (strbuff, "BlkStrU_%i.ps", i);
               int blks_3[4];
               blks_3[0] = 0;
               blks_3[1] = pnblks1_ilu2_arr[i];
               blks_3[2] = pnblks1_ilu2_arr[i] + pnblks2_ilu2_arr[i];
               blks_3[3] = pnblks_ilu2_arr[i];

               if (prec_float_loc == 1) {
                  CMatrix < int, float >*phmatr = pmatrU_float[i].GetHMatrStr ();
                  CBMatrix < int, float >::Str2PsBox (*phmatr, strbuff, 3, blks_3);
               } else {
                  CMatrix < int, float >*phmatr = pmatrU_double[i].GetHMatrStr ();
                  CBMatrix < int, float >::Str2PsBox (*phmatr, strbuff, 3, blks_3);
               }
            }

            phmatr_ext_arr[i].Clean ();

            if (_params->msglev >= 2) {
               cout << " End Fct Ihblk: " << i << endl;
            }

         }
      }

//   ffout_debug << " Point 8  " << endl;

      CMPIDataExchange::ExchangeArray (pcommloc, 'I', '+', 1, (void *) (&_stats->nmodif));

      CMPIDataExchange::ExchangeArray (pcommloc, 'D', '+', nhblksloc,
                                       (void *) (psclmin_arr));
      CMPIDataExchange::ExchangeArray (pcommloc, 'D', '+', nhblksloc,
                                       (void *) (psclmax_arr));
      CMPIDataExchange::ExchangeArray (pcommloc, 'D', '+', nhblksloc,
                                       (void *) (peigmin_arr));
      CMPIDataExchange::ExchangeArray (pcommloc, 'D', '+', nhblksloc,
                                       (void *) (peigmax_arr));

      _stats->scpiv_min = psclmin_arr[0];
      _stats->scpiv_max = psclmax_arr[0];
      _stats->piv_min = peigmin_arr[0];
      _stats->piv_max = peigmax_arr[0];

      for (i = 1; i < nhblksloc; i++) {
         if (psclmin_arr[i] < _stats->scpiv_min)
            _stats->scpiv_min = psclmin_arr[i];
         if (psclmax_arr[i] > _stats->scpiv_max)
            _stats->scpiv_max = psclmax_arr[i];
         if (peigmin_arr[i] < _stats->piv_min)
            _stats->piv_min = peigmin_arr[i];
         if (peigmax_arr[i] > _stats->piv_max)
            _stats->piv_max = peigmax_arr[i];
      }

      CMPIDataExchange::ExchangeArray (pcommloc, 'L', '+', 1, (void *) (&nzlu_tot));

      nzlu_tot -= ntot_ext;

      _stats->density = ((double) (nzlu_tot) / (double) (nza_tot));

// Init solve control structures

//   ffout_debug << " Point 9  " << endl;

      if (prec_float_loc == 1) {

         this->slv_float.InitControl (pcommloc, nhblksloc, phblks, phblk2cpu, phblk2blks,
                                      pblk2hblks, pblks);
         this->slv_float.InitSolveLU (phblks_ext, plist_ext_arr, pnblks_ext_arr,
                                      pblksnum_ext_arr, pblks_ext_arr, ptree_arr,
                                      pnblks_ilu2_arr, pnblks1_ilu2_arr, pnblks2_ilu2_arr,
                                      pblks_ilu2_arr, porder_LU, pmatrL_float,
                                      pmatrU_float);

      } else {

         this->slv_double.InitControl (pcommloc, nhblksloc, phblks, phblk2cpu, phblk2blks,
                                       pblk2hblks, pblks);
         this->slv_double.InitSolveLU (phblks_ext, plist_ext_arr, pnblks_ext_arr,
                                       pblksnum_ext_arr, pblks_ext_arr, ptree_arr,
                                       pnblks_ilu2_arr, pnblks1_ilu2_arr,
                                       pnblks2_ilu2_arr, pblks_ilu2_arr, porder_LU,
                                       pmatrL_double, pmatrU_double);

      }

// Return back unscaled matrix if necessary

      _params->sctype = sctype_save;

      if (this->b_use_blksize) {
         this->hmatr_arr.swap (this->hmatr_arr_temp);
         this->hmatr_arr_temp.resize (0);
         _stats->scpiv_min = scpiv_min_save;
         _stats->scpiv_max = scpiv_max_save;
      }
// Finalize timer

      CMPIDataExchange::Synchronize (pcommloc);

      double time1;

      time1 = CMPIDataExchange::GetWallTimeMPI ();

      _stats->dtime_fct = time1 - time0;

//   ffout_debug << " Return point !!! " << endl;

   }

// Peform explitic block scaling
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::ExplicitBlockScaling (SParams * _params, SStatData * _stats)
   {

// Open control data

      void *pcommloc = this->pcomm;

      int nhblksloc = this->nhblks;
      int nblksloc = this->nblks;

      this->params = *_params;

      int *phblk2cpu = &this->hblk2cpu[0];
      int *phblk2blks = &this->hblk2blks[0];
      long long *pblks = &this->blks[0];

      int myid = CMPIDataExchange::GetMyid (pcommloc);

      _params->iparam1 = myid;

      CBMatrix < _Int, _Flt > *phmatr_arr = &this->hmatr_arr[0];

// Allocate scaling data

      this->sclL.resize (nhblksloc + 1);
      this->sclU.resize (nhblksloc + 1);
      this->sclLInv.resize (nhblksloc + 1);
      this->sclUInv.resize (nhblksloc + 1);

      CBMatrix < _Int, _Flt > *pSclL_arr = &this->sclL[0];
      CBMatrix < _Int, _Flt > *pSclU_arr = &this->sclU[0];
      CBMatrix < _Int, _Flt > *pSclLInv_arr = &this->sclLInv[0];
      CBMatrix < _Int, _Flt > *pSclUInv_arr = &this->sclUInv[0];

// Transform matrix data into condensed block format

      vector < CBMatrix < _Int, _Flt > >hmatr_cnd_arr (nhblksloc + 1);
      CBMatrix < _Int, _Flt > *phmatr_cnd_arr = &hmatr_cnd_arr[0];

      int blksizeloc = _params->blksize_scale;

      this->blksize_bscl = blksizeloc;

      {
         int nimax = 0;
         {
            int i, ni;
            for (i = 0; i < nblksloc; i++) {
               ni = (int) (pblks[i + 1] - pblks[i]);
               if (ni > nimax)
                  nimax = ni;
            }
         }

         int nimax_cnd = nimax / blksizeloc;

         CVectorData < int >sp2blk (nimax);
         int *psp2blk = sp2blk.Ptr ();

         {
            int i, j;
            for (i = 0; i < nimax_cnd; i++) {
               for (j = 0; j < blksizeloc; j++) {
                  psp2blk[i * blksizeloc + j] = i;
               }
            }
         }

         int n_thr = 1;

#ifdef USE_THREADS
         n_thr = omp_get_max_threads ();
#endif

         vector < int >icycle_thr (n_thr + 1);
         vector < CVectorData < int > >imask_thr (n_thr + 1);

         int *picycle_thr = &icycle_thr[0];
         CVectorData < int >*pimask_thr = &imask_thr[0];

         int i;

         for (i = 0; i < n_thr; i++)
            picycle_thr[i] = -1;

         for (i = 0; i < nhblksloc; i++) {
            if (phblk2cpu[i] == myid) {

               int niblk = phblk2blks[i + 1] - phblk2blks[i];

               int nzblk = phmatr_arr[i].GetNzblk ();
               CMatrix < _Int, _Flt > *pASub = phmatr_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str = phmatr_arr[i].GetHMatrStr ();

               int *pia_hblk = phmatr_str->GetIaArr ();

               phmatr_cnd_arr[i].SetNzblk (nzblk);
               phmatr_cnd_arr[i].ResizeASub (nzblk);

               CMatrix < _Int, _Flt > *pASub_cnd = phmatr_cnd_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str_cnd = phmatr_cnd_arr[i].GetHMatrStr ();

               *phmatr_str_cnd = *phmatr_str;

               {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < niblk; ipar++) {

                     int my_thr = 0;
#ifdef USE_THREADS
                     my_thr = omp_get_thread_num ();
#endif

                     int icycle_th = picycle_thr[my_thr];

                     if (icycle_th == -1) {
                        pimask_thr[my_thr].resize (nimax * 6 + 1);
                        int *ppimask_thr = pimask_thr[my_thr].Ptr ();
                        int j;
                        for (j = 0; j < nimax; j++)
                           ppimask_thr[j] = -1;
                     }

                     int *ppimask_thr = pimask_thr[my_thr].Ptr ();

                     int jind;

                     for (jind = pia_hblk[ipar]; jind < pia_hblk[ipar + 1]; jind++) {

                        pASub[jind].CondenseSparsityBlksize (blksizeloc, psp2blk,
                                                             icycle_th, nimax,
                                                             ppimask_thr,
                                                             pASub_cnd[jind]);

                     }

                     picycle_thr[my_thr] = icycle_th;

                  }
               }
            }

         }

      }

// Compute explicit block scaling

      double sclmin_loc = _params->sclmin;
      double sclmin_att = 1.0e100;
      double sclmax_att = -1.0e100;
      int nmodif_att = 0;

      {

         int n_thr = 1;

#ifdef USE_THREADS
         n_thr = omp_get_max_threads ();
#endif

         CVectorData < double >sclmin_threads (n_thr);
         CVectorData < double >sclmax_threads (n_thr);
         CVectorData < int >nmodif_scl_threads (n_thr);

         double *psclmin_threads = sclmin_threads.Ptr ();
         double *psclmax_threads = sclmax_threads.Ptr ();
         int *pnmodif_scl_threads = nmodif_scl_threads.Ptr ();

         CVectorData < double >sclmin_hblk_threads (n_thr);
         CVectorData < double >sclmax_hblk_threads (n_thr);
         CVectorData < int >nmodif_scl_hblk_threads (n_thr);

         double *psclmin_hblk_threads = sclmin_hblk_threads.Ptr ();
         double *psclmax_hblk_threads = sclmax_hblk_threads.Ptr ();
         int *pnmodif_scl_hblk_threads = nmodif_scl_hblk_threads.Ptr ();

         int i;

         for (i = 0; i < n_thr; i++)
            psclmin_threads[i] = 1.0e100;
         for (i = 0; i < n_thr; i++)
            psclmax_threads[i] = -1.0e100;
         for (i = 0; i < n_thr; i++)
            pnmodif_scl_threads[i] = 0;

         for (i = 0; i < nhblksloc; i++) {
            if (phblk2cpu[i] == myid) {

               int niblk = phblk2blks[i + 1] - phblk2blks[i];
               int ibegblk = phblk2blks[i];

               pSclL_arr[i].SetNzblk (niblk);
               pSclL_arr[i].ResizeASub (niblk);

               CMatrix < _Int, _Flt > *pASub_sclL = pSclL_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str_sclL = pSclL_arr[i].GetHMatrStr ();

               phmatr_str_sclL->ResizeAndSetAllSp (niblk, niblk, niblk, niblk);

               int *plist_sclL = phmatr_str_sclL->GetListArr ();
               int *plist2_sclL = phmatr_str_sclL->GetList2Arr ();
               int *pia_sclL = phmatr_str_sclL->GetIaArr ();
               int *pja_sclL = phmatr_str_sclL->GetJaArr ();
               int *pja2_sclL = phmatr_str_sclL->GetJa2Arr ();

               {

                  int j;

                  for (j = 0; j < niblk; j++)
                     plist_sclL[j] = j;
                  for (j = 0; j < niblk; j++)
                     plist2_sclL[j] = i;
                  for (j = 0; j <= niblk; j++)
                     pia_sclL[j] = j;
                  for (j = 0; j < niblk; j++)
                     pja_sclL[j] = j;
                  for (j = 0; j < niblk; j++)
                     pja2_sclL[j] = i;

               }

               pSclU_arr[i].SetNzblk (niblk);
               pSclU_arr[i].ResizeASub (niblk);

               CMatrix < _Int, _Flt > *pASub_sclU = pSclU_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str_sclU = pSclU_arr[i].GetHMatrStr ();

               *phmatr_str_sclU = *phmatr_str_sclL;

               pSclLInv_arr[i].SetNzblk (niblk);
               pSclLInv_arr[i].ResizeASub (niblk);

               CMatrix < _Int, _Flt > *pASub_sclLInv = pSclLInv_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str_sclLInv = pSclLInv_arr[i].GetHMatrStr ();

               *phmatr_str_sclLInv = *phmatr_str_sclL;

               pSclUInv_arr[i].SetNzblk (niblk);
               pSclUInv_arr[i].ResizeASub (niblk);

               CMatrix < _Int, _Flt > *pASub_sclUInv = pSclUInv_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str_sclUInv = pSclUInv_arr[i].GetHMatrStr ();

               *phmatr_str_sclUInv = *phmatr_str_sclL;

               CMatrix < _Int, _Flt > *pASub = phmatr_cnd_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str = phmatr_cnd_arr[i].GetHMatrStr ();

               int *pia_hblk = phmatr_str->GetIaArr ();
               int *pja_hblk = phmatr_str->GetJaArr ();
               int *pja2_hblk = phmatr_str->GetJa2Arr ();

               {
                  int jjk;
                  for (jjk = 0; jjk < n_thr; jjk++)
                     psclmin_hblk_threads[jjk] = 1.0e100;
                  for (jjk = 0; jjk < n_thr; jjk++)
                     psclmax_hblk_threads[jjk] = -1.0e100;
                  for (jjk = 0; jjk < n_thr; jjk++)
                     pnmodif_scl_hblk_threads[jjk] = 0;
               }

               {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < niblk; ipar++) {

                     int my_thr = 0;
#ifdef USE_THREADS
                     my_thr = omp_get_thread_num ();
#endif

                     int iblkgl = ibegblk + ipar;
                     int nlistloc = (int) (pblks[iblkgl + 1] - pblks[iblkgl]);
                     nlistloc = nlistloc / blksizeloc;

                     int j, jblk, jhblk;
                     bool b_found = false;
                     int ind = -1;

                     for (j = pia_hblk[ipar]; j < pia_hblk[ipar + 1]; j++) {
                        jblk = pja_hblk[j];
                        jhblk = pja2_hblk[j];
                        if (jblk == ipar && jhblk == i) {
                           b_found = true;
                           ind = j;
                        }
                     }

                     if (!b_found) {
                        cout <<
                           " ExplicitBlockScaling: error: diagonal block is not found !!! "
                           << endl;
                        throw
                           " ExplicitBlockScaling: error: diagonal block is not found !!! ";
                     }

                     double sclmin_att_temp;
                     double sclmax_att_temp;
                     int nmodif_scl_temp;

                     pASub[ind].BlockScaling (blksizeloc, sclmin_loc, nlistloc,
                                              pASub_sclL[ipar], pASub_sclU[ipar],
                                              pASub_sclLInv[ipar], pASub_sclUInv[ipar],
                                              sclmin_att_temp, sclmax_att_temp,
                                              nmodif_scl_temp);

                     if (sclmin_att_temp < psclmin_threads[my_thr])
                        psclmin_threads[my_thr] = sclmin_att_temp;
                     if (sclmax_att_temp > psclmax_threads[my_thr])
                        psclmax_threads[my_thr] = sclmax_att_temp;
                     pnmodif_scl_threads[my_thr] += nmodif_scl_temp;

                     if (sclmin_att_temp < psclmin_hblk_threads[my_thr])
                        psclmin_hblk_threads[my_thr] = sclmin_att_temp;
                     if (sclmax_att_temp > psclmax_hblk_threads[my_thr])
                        psclmax_hblk_threads[my_thr] = sclmax_att_temp;
                     pnmodif_scl_hblk_threads[my_thr] += nmodif_scl_temp;

                  }
               }

               double sclmin_att_hblk = 1.0e100;
               double sclmax_att_hblk = -1.0e100;
               int nmodif_att_hblk = 0;

               {
                  int jjk;
                  for (jjk = 0; jjk < n_thr; jjk++) {
                     if (psclmin_hblk_threads[jjk] < sclmin_att_hblk)
                        sclmin_att_hblk = psclmin_hblk_threads[jjk];
                     if (psclmax_hblk_threads[jjk] > sclmax_att_hblk)
                        sclmax_att_hblk = psclmax_hblk_threads[jjk];
                     nmodif_att_hblk += pnmodif_scl_hblk_threads[jjk];
                  }
               }

               if (_params->msglev >= 2) {
                  cout << "  Ihblk = " << i << " NModif = " << nmodif_att_hblk <<
                     " SclMin = " << sclmin_att_hblk << " SclMax = " << sclmax_att_hblk <<
                     endl;
               }

            }
         }

         for (i = 0; i < n_thr; i++) {
            if (psclmin_threads[i] < sclmin_att)
               sclmin_att = psclmin_threads[i];
            if (psclmax_threads[i] > sclmax_att)
               sclmax_att = psclmax_threads[i];
            nmodif_att += pnmodif_scl_threads[i];
         }

         CMPIDataExchange::ExchangeArray (pcommloc, 'D', 'm', 1, &sclmin_att);
         CMPIDataExchange::ExchangeArray (pcommloc, 'D', 'M', 1, &sclmax_att);
         CMPIDataExchange::ExchangeArray (pcommloc, 'I', '+', 1, &nmodif_att);

      }

      _stats->scpiv_min = sclmin_att;
      _stats->scpiv_max = sclmax_att;
      _stats->nmodif_scl = nmodif_att;

// Apply explicit block scaling

      CBMatrix < _Int, _Flt >::ExplicitBlockScaling (pcommloc, blksizeloc, nhblksloc,
                                                     phblk2cpu, phblk2blks, pblks,
                                                     phmatr_cnd_arr, pSclL_arr,
                                                     pSclU_arr);

// Transform matrix back into point format

      this->hmatr_arr_temp.resize (nhblksloc + 1);
      CBMatrix < _Int, _Flt > *phmatr_arr_temp = &hmatr_arr_temp[0];

      {
         int nimax = 0;
         {
            int i, ni;
            for (i = 0; i < nblksloc; i++) {
               ni = (int) (pblks[i + 1] - pblks[i]);
               if (ni > nimax)
                  nimax = ni;
            }
         }

         int nimax_cnd = nimax / blksizeloc;

         CVectorData < int >sp2blk (nimax);
         int *psp2blk = sp2blk.Ptr ();

         {
            int i, j;
            for (i = 0; i < nimax_cnd; i++) {
               for (j = 0; j < blksizeloc; j++) {
                  psp2blk[i * blksizeloc + j] = i;
               }
            }
         }

         int n_thr = 1;

#ifdef USE_THREADS
         n_thr = omp_get_max_threads ();
#endif

         vector < int >icycle_thr (n_thr + 1);
         vector < CVectorData < int > >imask_thr (n_thr + 1);

         int *picycle_thr = &icycle_thr[0];
         CVectorData < int >*pimask_thr = &imask_thr[0];

         int i;

         for (i = 0; i < n_thr; i++)
            picycle_thr[i] = -1;

         for (i = 0; i < nhblksloc; i++) {
            if (phblk2cpu[i] == myid) {

               int niblk = phblk2blks[i + 1] - phblk2blks[i];
               int ibegblk = phblk2blks[i];

               CVectorData < int >nblk_colmax_arr (niblk);
               int *pnblk_colmax_arr = nblk_colmax_arr.Ptr ();

               int nzblk = phmatr_cnd_arr[i].GetNzblk ();
               CMatrix < _Int, _Flt > *pASub = phmatr_cnd_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str = phmatr_cnd_arr[i].GetHMatrStr ();

               int *pia_hblk = phmatr_str->GetIaArr ();

               phmatr_arr_temp[i].SetNzblk (nzblk);
               phmatr_arr_temp[i].ResizeASub (nzblk);

               CMatrix < _Int, _Flt > *pASub_ucnd = phmatr_arr_temp[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str_ucnd = phmatr_arr_temp[i].GetHMatrStr ();

               *phmatr_str_ucnd = *phmatr_str;

               {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < niblk; ipar++) {

                     int my_thr = 0;
#ifdef USE_THREADS
                     my_thr = omp_get_thread_num ();
#endif

                     int icycle_th = picycle_thr[my_thr];

                     if (icycle_th == -1) {
                        pimask_thr[my_thr].resize (nimax * 6 + 1);
                        int *ppimask_thr = pimask_thr[my_thr].Ptr ();
                        int j;
                        for (j = 0; j < nimax; j++)
                           ppimask_thr[j] = -1;
                     }

                     int *ppimask_thr = pimask_thr[my_thr].Ptr ();

                     int ni_cnd =
                        (int) ((pblks[ibegblk + ipar + 1] -
                                pblks[ibegblk + ipar]) / blksizeloc);

                     CVectorData < int >isize_arr (ni_cnd);
                     int *pisize_arr = isize_arr.Ptr ();

                     int k;

                     for (k = 0; k < ni_cnd; k++)
                        pisize_arr[k] = 0;

                     int jind, kk;

                     for (jind = pia_hblk[ipar]; jind < pia_hblk[ipar + 1]; jind++) {

                        int nlist_temp = pASub[jind].GetNlist ();
                        _Int *plist_temp = pASub[jind].GetListArr ();
                        _Int *pia_temp = pASub[jind].GetIaArr ();

                        for (k = 0; k < nlist_temp; k++) {
                           kk = (int) plist_temp[k];
                           pisize_arr[kk] += (int) (pia_temp[k + 1] - pia_temp[k]);
                        }

                        pASub[jind].UnCondenseSparsityBlksize (blksizeloc, psp2blk,
                                                               icycle_th, nimax,
                                                               ppimask_thr,
                                                               pASub_ucnd[jind]);

                     }

                     int isize_max = 0;

                     for (k = 0; k < ni_cnd; k++) {
                        if (pisize_arr[k] > isize_max)
                           isize_max = pisize_arr[k];
                     }

                     pnblk_colmax_arr[ipar] = isize_max;

                     picycle_thr[my_thr] = icycle_th;

                  }
               }

               {

                  int isize_hblk_max = 0;
                  int j;
                  for (j = 0; j < niblk; j++) {
                     if (pnblk_colmax_arr[j] > isize_hblk_max)
                        isize_hblk_max = pnblk_colmax_arr[j];
                  }

                  if (_params->msglev >= 2) {
                     cout << " Ihblk = " << i << " Max_Blk in Block Scale = " <<
                        isize_hblk_max << endl;
                  }

               }

            }

         }

      }

// Swap work and temp matrix data

      this->hmatr_arr.swap (this->hmatr_arr_temp);

   }

// Perform iterations of the BiCGStab iterative scheme
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > bool CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::BiCGStab (bool _b_use_poly, SParams * _params, SStatData * _stats,
                            void *_str_mvmA,
                            typename TMvmSlvFunc < _FltVect >::mvmA_f _mvm_f,
                            void *_str_slvLU,
                            typename TMvmSlvFunc < _FltVect >::slvLU_f _slv_f,
                            _FltVect * _rhs, _FltVect * _sol)
   {

      bool conv = false;

// Init output data

      int blksize_iter = _params->blksize_iter;

      _stats->rhs_norm = 0.0e0;
      _stats->res_norm = 0.0e0;
      _stats->niter = 0;
      _stats->nmvm = 0;
      _stats->res_fin = 0.0e0;
      _stats->dtime_iter = 0.0e0;

// Get control data

      void *pcomm_loc = this->GetComm ();
      int nhblks_loc = this->GetNhblks ();
      long long *phblks_loc = this->GetHBlks ();
      int *phblks2cpu_loc = this->GetHBlk2cpu ();

      int myid_loc = CMPIDataExchange::GetMyid (pcomm_loc);

// Init timer

      CMPIDataExchange::Synchronize (pcomm_loc);

      double time0;

      time0 = CMPIDataExchange::GetWallTimeMPI ();

// Get the local size of vector data

      int i;

      int n_local = 0;

      for (i = 0; i < nhblks_loc; i++) {
         if (phblks2cpu_loc[i] == myid_loc) {
            n_local += (int) (phblks_loc[i + 1] - phblks_loc[i]);
         }
      }

      n_local *= blksize_iter;

// Allocate work vector arrays

      CVectorData < _FltVect > rhs (n_local + 1);
      CVectorData < _FltVect > r (n_local + 1);
      CVectorData < _FltVect > p (n_local + 1);
      CVectorData < _FltVect > x (n_local + 1);
      CVectorData < _FltVect > u (n_local + 1);
      CVectorData < _FltVect > z (n_local + 1);
      CVectorData < _FltVect > v (n_local + 1);
      CVectorData < _FltVect > q (n_local + 1);

      _FltVect *pr = r.Ptr ();
      _FltVect *pp = p.Ptr ();
      _FltVect *px = x.Ptr ();
      _FltVect *pu = u.Ptr ();
      _FltVect *pz = z.Ptr ();
      _FltVect *pv = v.Ptr ();
      _FltVect *pq = q.Ptr ();

// Compute initial residual vector and its norm

      CVector < _FltVect >::CopyVector_thr (n_local, _sol, px);
      CVector < _FltVect >::CopyVector_thr (n_local, _rhs, pr);

      (_mvm_f) (_str_mvmA, px, pz);

      CVector < _FltVect >::SubtractReplaceVector_thr (n_local, pz, pr);

      _FltVect rhs_norm = 0.0e0;
      _FltVect resi0_norm = 0.0e0;

      rhs_norm = CVector < _FltVect >::ScProd_thr (n_local, _rhs, _rhs);
      resi0_norm = CVector < _FltVect >::ScProd_thr (n_local, pr, pr);

      double sum2_arr[2];

      sum2_arr[0] = (double) rhs_norm;
      sum2_arr[1] = (double) resi0_norm;

      CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', 2, sum2_arr);

      double d_rhs_norm = sum2_arr[0];
      double d_resi0_norm = sum2_arr[1];

      d_rhs_norm = sqrt (d_rhs_norm);
      d_resi0_norm = sqrt (d_resi0_norm);

      _stats->rhs_norm = d_rhs_norm;
      _stats->res_norm = d_resi0_norm;

      if (_params->msglev >= 2 && myid_loc == 0)
         std::cout << " Log10 || Rhs || = " << log10 (d_rhs_norm) << std::endl;
      if (_params->msglev >= 1 && _params->pfout != NULL && myid_loc == 0)
         *_params->pfout << " Log10 || Rhs || = " << log10 (d_rhs_norm) << std::endl;

      if (_params->msglev >= 2 && myid_loc == 0)
         std::cout << " Initial Log10 || Resi || = " << log10 (d_resi0_norm) << std::endl;
      if (_params->msglev >= 1 && _params->pfout != NULL && myid_loc == 0)
         *_params->pfout << " Initial Log10 || Resi || = " << log10 (d_resi0_norm) <<
            std::endl;

      if ((_stats->res_norm < _params->eps_rel * _stats->rhs_norm)
          || (_stats->res_norm < _params->eps_abs)) {
         _stats->res_fin = _stats->res_norm;
         _stats->niter = 0;
         conv = true;
         return conv;
      }
// Perform iterations starting from residual data

      int niter_perf = 0;
      double resi_norm = d_resi0_norm;
      double d_res_min = resi_norm;

// Choose initial direction vector

      CVector < _FltVect >::CopyVector_thr (n_local, pr, pz);   // z=r(0)

// p(1) = r(0)

      CVector < _FltVect >::CopyVector_thr (n_local, pr, pp);   // p=r(0)

// Main iterative cycle

      _FltVect alpha = 0.0e0, beta, omega = 0.0e0, rho, rhoold; // method parameters
      _FltVect uu, ur, zv;      // inner products

      _FltVect fone = (_FltVect) 1.0e0;

      rho = fone;

      int k, it;

      _FltVect faux;
      double daux;

      double d_resi;

      for (k = 1; k <= _params->niter_max; k++) {

         it = k - 1;            // the number of MVM is equal to 2*it

// rho(i-1) = z^T * r(i-1)

         rhoold = rho;

         rho = CVector < _FltVect >::ScProd_thr (n_local, pz, pr);

         daux = (double) rho;

         CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &daux);

         rho = (_FltVect) daux;

         if (k == 1) {
            alpha = omega = fone;
         } else {

// beta(i-1) = (rho(i-1) / rho(i-2)) * (alpha(i-1) / omega(i-1))

            beta = (rho / rhoold) * (alpha / omega);

// p(i) = r(i-1) + beta(i-1) * (p(i-1) - omega(i-1) * v(i-1))

            CVector < _FltVect >::UpdateVectorMinus_thr (n_local, &omega, pv, pp);
            CVector < _FltVect >::UpdateVectorReversed_thr (n_local, &beta, pr, pp);

         }

// q = M_solve * p(i)

         this->PolySlvLU (_b_use_poly, _str_mvmA, _mvm_f, _str_slvLU, _slv_f, pp, pq);

// v(i) = A * q

         (_mvm_f) (_str_mvmA, pq, pv);

// alpha(i) = rho(i-1) / z^T * v(i)

         zv = CVector < _FltVect >::ScProd_thr (n_local, pz, pv);

         daux = (double) zv;

         CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &daux);

         zv = (_FltVect) daux;

         alpha = rho / zv;

// x(i-1/2) = x(i-1) + alpha(i) * q

         CVector < _FltVect >::UpdateVector_thr (n_local, &alpha, pq, px);

// r(i-1/2) = r(i-1) - alpha(i) * v(i)

         CVector < _FltVect >::UpdateVectorMinus_thr (n_local, &alpha, pv, pr);

// The intermediate check of convergence can be added here

         niter_perf = (int) (it + 1);

         faux = CVector < _FltVect >::ScProd_thr (n_local, pr, pr);

         daux = (double) faux;

         CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &daux);

         d_resi = sqrt (daux);

         if (d_resi == 0.0e0)
            d_resi = d_rhs_norm * 1.0e-16;

         if ((d_resi < _params->eps_rel * d_rhs_norm) || (d_resi < _params->eps_abs)) {

            conv = true;

// Save the best attained solution

            if (d_resi < d_res_min) {
               d_res_min = d_resi;
               CVector < _FltVect >::CopyVector_thr (n_local, px, _sol);
            }

            break;
         }
// q = M_solve * r(i-1/2)

         this->PolySlvLU (_b_use_poly, _str_mvmA, _mvm_f, _str_slvLU, _slv_f, pr, pq);

// u = A * q

         (_mvm_f) (_str_mvmA, pq, pu);

// omega(i) = u^T * r(i-1/2) / u^T * u

         ur = CVector < _FltVect >::ScProd_thr (n_local, pu, pr);
         uu = CVector < _FltVect >::ScProd_thr (n_local, pu, pu);

         sum2_arr[0] = (double) ur;
         sum2_arr[1] = (double) uu;

         CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', 2, sum2_arr);

         ur = (_FltVect) sum2_arr[0];
         uu = (_FltVect) sum2_arr[1];

         omega = ur / uu;

// x(i) = x(i-1/2) + omega(i) * q

         CVector < _FltVect >::UpdateVector_thr (n_local, &omega, pq, px);

// r(i) = r(i-1/2) - omega(i) * u

         CVector < _FltVect >::UpdateVectorMinus_thr (n_local, &omega, pu, pr);

// Check the convergence

         niter_perf = (int) (it + 1);

         faux = CVector < _FltVect >::ScProd_thr (n_local, pr, pr);

         daux = (double) faux;

         CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &daux);

         d_resi = sqrt (daux);

         if ((it + _params->ichk) % _params->ichk == 0 && _params->msglev >= 2
             && myid_loc == 0)
            std::cout << " It = " << it << " Log10 || Resi || = " << log10 (d_resi) <<
               std::endl;
         if ((it + _params->ichk) % _params->ichk == 0 && _params->msglev >= 1
             && _params->pfout != NULL && myid_loc == 0)
            *_params->pfout << " It = " << it << " Log10 || Resi || = " << log10 (d_resi)
               << std::endl;

         if ((d_resi < _params->eps_rel * d_rhs_norm) || (d_resi < _params->eps_abs))
            conv = true;

// Save the best attained solution

         if (d_resi < d_res_min) {
            d_res_min = d_resi;
            CVector < _FltVect >::CopyVector_thr (n_local, px, _sol);
         }
// Break from iterations if converged

         if (conv) {
            break;
         }

      }                         // end of iterations

// Compute the final residual

      (_mvm_f) (_str_mvmA, _sol, pr);

      CVector < _FltVect >::SubtractReplaceVector_thr (n_local, _rhs, pr);      // r=b-A*x

// Compute the final residual norm

      resi_norm = CVector < _FltVect >::ScProd_thr (n_local, pr, pr);

      daux = (double) resi_norm;

      CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &daux);

      double d_resi_fin = sqrt (daux);

      if (d_resi_fin == 0.0e0)
         d_resi_fin = d_rhs_norm * 1.0e-16;

      if (_params->msglev >= 2 && myid_loc == 0)
         std::cout << " Final Log10 || Resi || = " << log10 (d_resi_fin) << std::endl;
      if (_params->msglev >= 1 && _params->pfout != NULL && myid_loc == 0)
         *_params->pfout << " Final Log10 || Resi || = " << log10 (d_resi_fin) << std::
            endl;

// Finalize timer

      CMPIDataExchange::Synchronize (pcomm_loc);

      double time1;

      time1 = CMPIDataExchange::GetWallTimeMPI ();

      _stats->dtime_iter = time1 - time0;

// Return statistics data

      _stats->niter = niter_perf;
      _stats->res_fin = d_resi_fin;
      _stats->nmvm = _stats->niter * 2;
      if (_b_use_poly) {
         int ncoef_loc = this->ncoef_slv;
         _stats->nmvm = _stats->niter * 2 * ncoef_loc;
      }

      return conv;

   }

// Perform iterations of the Gmres iterative scheme
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > bool CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::Gmres (bool _b_use_poly, SParams * _params, SStatData * _stats,
                         void *_str_mvmA,
                         typename TMvmSlvFunc < _FltVect >::mvmA_f _mvm_f,
                         void *_str_slvLU,
                         typename TMvmSlvFunc < _FltVect >::slvLU_f _slv_f,
                         _FltVect * _rhs, _FltVect * _sol) {

      bool conv = false;

// Init output data

      int blksize_loc = _params->blksize_iter;

      _stats->rhs_norm = 0.0e0;
      _stats->res_norm = 0.0e0;
      _stats->niter = 0;
      _stats->res_fin = 0.0e0;
      _stats->dtime_iter = 0.0e0;

// Get control data

      void *pcomm_loc = this->GetComm ();
      int nhblks_loc = this->GetNhblks ();
      long long *phblks_loc = this->GetHBlks ();
      int *phblks2cpu_loc = this->GetHBlk2cpu ();

      int myid_loc = CMPIDataExchange::GetMyid (pcomm_loc);

// Init timer

      CMPIDataExchange::Synchronize (pcomm_loc);

      double time0;

      time0 = CMPIDataExchange::GetWallTimeMPI ();

// Get the local size of vector data

      int i;

      int n_local = 0;

      for (i = 0; i < nhblks_loc; i++) {
         if (phblks2cpu_loc[i] == myid_loc) {
            n_local += (int) (phblks_loc[i + 1] - phblks_loc[i]);
         }
      }

      n_local *= blksize_loc;

// Allocate work vector arrays

      CVectorData < _FltVect > r (n_local + 1);
      CVectorData < _FltVect > x (n_local + 1);
      CVectorData < _FltVect > u (n_local + 1);
      CVectorData < _FltVect > z (n_local + 1);
      CVectorData < _FltVect > v (n_local + 1);

      _FltVect *pr = r.Ptr ();
      _FltVect *px = x.Ptr ();
      _FltVect *pu = u.Ptr ();
      _FltVect *pz = z.Ptr ();
      _FltVect *pv = v.Ptr ();

// Compute initial residual vector and its norm

      CVector < _FltVect >::CopyVector_thr (n_local, _sol, px);
      CVector < _FltVect >::CopyVector_thr (n_local, _rhs, pr);

      (_mvm_f) (_str_mvmA, px, pz);

      CVector < _FltVect >::SubtractReplaceVector_thr (n_local, pz, pr);

      CVector < _FltVect >::CopyVector_thr (n_local, pr, pu);

      _FltVect rhs_norm = 0.0e0;
      _FltVect resi0_norm = 0.0e0;

      rhs_norm = CVector < _FltVect >::ScProd_thr (n_local, _rhs, _rhs);
      resi0_norm = CVector < _FltVect >::ScProd_thr (n_local, pr, pr);

      double sum2_arr[2];

      sum2_arr[0] = (double) rhs_norm;
      sum2_arr[1] = (double) resi0_norm;

      CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', 2, sum2_arr);

      double d_rhs_norm = sum2_arr[0];
      double d_resi0_norm = sum2_arr[1];

      d_rhs_norm = sqrt (d_rhs_norm);
      d_resi0_norm = sqrt (d_resi0_norm);

      _stats->rhs_norm = d_rhs_norm;
      _stats->res_norm = d_resi0_norm;

      if (myid_loc == 0 && _params->msglev >= 2)
         std::cout << " Log10 || Rhs || = " << log10 (d_rhs_norm) << std::endl;
      if (myid_loc == 0 && _params->msglev >= 1 && _params->pfout != NULL)
         *_params->pfout << " Log10 || Rhs || = " << log10 (d_rhs_norm) << std::endl;

      if (myid_loc == 0 && _params->msglev >= 2)
         std::cout << " Initial Log10 || Resi || = " << log10 (d_resi0_norm) << std::endl;
      if (myid_loc == 0 && _params->msglev >= 1 && _params->pfout != NULL)
         *_params->pfout << " Initial Log10 || Resi || = " << log10 (d_resi0_norm) <<
            std::endl;

      if ((_stats->res_norm < _params->eps_rel * _stats->rhs_norm)
          || (_stats->res_norm < _params->eps_abs)) {
         _stats->res_fin = _stats->res_norm;
         _stats->niter = 0;
         conv = true;
         return conv;
      }
// Perform iterations starting from residual data

// Create qrd data

      CQrdMPI < _FltVect > qrdMPI_P;

      qrdMPI_P.Init (pcomm_loc, n_local);

      CQrdBase < _FltVect > qrdGIVENS;

// Get parameters

      int niter = _params->niter_max;
      int nitcycle = _params->niter_cycle;

      int nitmax = niter;

      int nitcycle2 = nitcycle + 2;

      int nitmax_alloc = nitmax;
      if (nitcycle > nitmax_alloc)
         nitmax_alloc = nitcycle;

// double eps_achieved = 1; //done_static;

// Estimate the maximal number of columns in R

      int ncolsRmax = nitcycle + 2;

// Determine the sizes arrays

      int ncolsmax = ncolsRmax;

// Allocate work arrays

      vector < _FltVect > alpha (ncolsmax + 1);
      vector < _FltVect > hmatr (ncolsmax * nitcycle2 + 1),
         rmatr (ncolsmax * nitcycle2 + 1);
      vector < double >diagR (ncolsmax + 1);

      _FltVect *palpha, *phmatr, *prmatr;
      double *pdiagR;

      palpha = &alpha[0];
      phmatr = &hmatr[0];
      prmatr = &rmatr[0];
      pdiagR = &diagR[0];

      CVector < _FltVect >::SetByZeroes (ncolsmax * nitcycle2, phmatr);
      CVector < _FltVect >::SetByZeroes (ncolsmax * nitcycle2, prmatr);

// Perform Gmres iterations

      int iterloc, ncolsP;
      double d_resi_norm;

      d_resi_norm = d_resi0_norm;

      int itgl = 0;
      int iter = 0;

      while (iter < nitmax && !conv) {

         if (myid_loc == 0 && _params->msglev >= 2)
            std::cout << " ItGl = " << itgl << std::endl;
         if (myid_loc == 0 && _params->msglev >= 1 && _params->pfout != NULL)
            *_params->pfout << " ItGl = " << itgl << std::endl;

// Check convergence for current guess

         {

            _FltVect resi0_2;

            resi0_2 = CVector < _FltVect >::ScProd_thr (n_local, pu, pu);

            double d_resi0_2 = (double) resi0_2;

            CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &d_resi0_2);

            d_resi_norm = sqrt (d_resi0_2);

//         if (myid_loc == 0 && _msglev >= 2) std::cout << "   It = " << iter << " Log10 || Resi || = " << log10(d_resi_norm) << std::endl;
//         if (myid_loc == 0 && _msglev >= 1 && _fout != NULL) *_fout   << "   It = " << iter << " Log10 || Resi || = " << log10(d_resi_norm) << std::endl;

            if ((d_resi_norm < _params->eps_rel * d_rhs_norm)
                || (d_resi_norm < _params->eps_abs))
               conv = true;

         }

// Fast return if necessary

         if (d_resi_norm == 0.0e0)
            conv = true;

         if (conv)
            break;

// Add residual data into the W qrd

         ncolsP = 0;

         qrdMPI_P.UpdateQrdMPI_thr (1, n_local, pu, n_local);
         ncolsP++;

// Get local rhs coefs

         CVector < _FltVect >::SetByZeroes (ncolsmax, palpha);

         d_resi_norm = 0.0e0;

         if (myid_loc == 0) {
            qrdMPI_P.GetRQrdMPI (ncolsP - 1, ncolsP - 1, 0, ncolsP - 1, palpha, ncolsmax);
            d_resi_norm = palpha[0] * palpha[0];
         }

         CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &d_resi_norm);

         d_resi_norm = sqrt (d_resi_norm);

// Perform local iterations

         itgl++;

         iterloc = 0;

         while (iterloc < nitcycle && iter < nitmax) {

// Check or estimate the convergence

            if (iterloc % _params->ichk == 0 && iterloc > 0) {
               if (myid_loc == 0 && _params->msglev >= 2)
                  std::cout << "   It = " << iter << " Log10 || Resi || = " <<
                     log10 (d_resi_norm) << std::endl;
               if (myid_loc == 0 && _params->msglev >= 1 && _params->pfout != NULL)
                  *_params->pfout << "   It = " << iter << " Log10 || Resi || = " <<
                     log10 (d_resi_norm) << std::endl;
            }

            if ((d_resi_norm < _params->eps_rel * d_rhs_norm)
                || (d_resi_norm < _params->eps_abs))
               conv = true;

            if (conv)
               break;

// New iteration

            iter++;
            iterloc++;

// Compute new direction vector

            if (myid_loc == 0) {
               CVector < _FltVect >::SetByZeroes (ncolsP, pz);
               CVector < _FltVect >::SetByOnes (1, pz + ncolsP - 1);
            }

            qrdMPI_P.MvmQMPI_thr (1, pz, ncolsP, pv, n_local);

// Multiply by the preconditioned matrix

            this->PolySlvLU (_b_use_poly, _str_mvmA, _mvm_f, _str_slvLU, _slv_f, pv, pz);
            (_mvm_f) (_str_mvmA, pz, pu);

// Update QR decomposition

            qrdMPI_P.UpdateQrdMPI_thr (1, n_local, pu, n_local);
            ncolsP++;

// Get current R part

            if (myid_loc == 0) {

               CVector < _FltVect >::SetByZeroes (ncolsmax,
                                                  phmatr + (iterloc - 1) * ncolsmax);
               qrdMPI_P.GetRQrdMPI (ncolsP - 1, ncolsP - 1, 0, ncolsP - 1,
                                    phmatr + (iterloc - 1) * ncolsmax, ncolsmax);

// Compute new Givens rotation and apply it to the current column

               qrdGIVENS.UpdateQrdBlk (1, ncolsP, phmatr + (iterloc - 1) * ncolsmax,
                                       ncolsmax);

               qrdGIVENS.GetRQrd (ncolsP - 2, ncolsP - 2, 0, ncolsP - 2,
                                  prmatr + (iterloc - 1) * ncolsmax, ncolsmax);

// Compute new reduced residual vector

               qrdGIVENS.MvmQHPart (1, iterloc - 1, iterloc - 1, palpha, ncolsmax);

            }
// Estimate residual norm

            if (myid_loc == 0) {
               _FltVect resi;
               resi =
                  CVector < _FltVect >::ScProd (1, palpha + iterloc, palpha + iterloc);
               d_resi_norm = (double) resi;
               d_resi_norm = sqrt (d_resi_norm);
            } else {
               d_resi_norm = 0.0e0;
            }

            CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &d_resi_norm);

         }

// Compute new local guess to the solution and coefs if necessary

         int ncoef = _params->ncoef;
         vector < double >coefs;

         if (myid_loc == 0) {

// Compute diagonal values of R

            for (i = 0; i < iterloc; i++) {
               pdiagR[i] = (double) prmatr[i * ncolsmax + i];
               if (pdiagR[i] < 0.0e0)
                  pdiagR[i] = -pdiagR[i];
            }

            double diagR_min = pdiagR[0];
            double diagR_max = pdiagR[0];

            for (i = 0; i < iterloc; i++) {
               if (pdiagR[i] > diagR_max)
                  diagR_max = pdiagR[i];
               if (pdiagR[i] < diagR_min)
                  diagR_min = pdiagR[i];
            }

            if (_params->msglev >= 2)
               cout << "   DiagR_min = " << diagR_min << " DiagR_max = " << diagR_max <<
                  endl;
            if (_params->msglev >= 1 && _params->pfout != NULL)
               *_params->pfout << "   DiagR_min = " << diagR_min << " DiagR_max = " <<
                  diagR_max << endl;

            if (_params->msglev >= 4)
               PrintArrayLow (cout, "DiagR", iterloc, pdiagR);
            if (_params->msglev >= 3 && _params->pfout != NULL)
               PrintArrayLow (*_params->pfout, "DiagR", iterloc, pdiagR);

// Compute new local guess

#ifdef USE_LAPACK

            {

               int n_temp = iterloc;

               CVectorData < _FltVect > RMatr (n_temp * n_temp);
               _FltVect *pRMatr = RMatr.Ptr ();

               CVector < _FltVect >::SetByZeroes (n_temp * n_temp, pRMatr);

               int i;

               for (i = 0; i < n_temp; i++) {
                  CVector < _FltVect >::CopyVector (i + 1, prmatr + i * ncolsmax,
                                                    pRMatr + i * n_temp);
               }

               CVectorData < _FltVect > work_svd (3 * n_temp * n_temp + 3 * n_temp);

               _FltVect *pSv_temp = work_svd.Ptr ();

               _FltVect *pU_temp = pSv_temp + n_temp;
               _FltVect *pV_temp = pU_temp + n_temp * n_temp;
               _FltVect *pX_temp = pV_temp + n_temp * n_temp;

               CVectorData < double >dwork_svd (4 * n_temp * n_temp + 20 * n_temp);

               double *pdwork_svd = dwork_svd.Ptr ();

               CVector < _FltVect >::ComputeSvd (n_temp, pRMatr, pSv_temp, pU_temp,
                                                 pV_temp, pdwork_svd);

               double svR_min = (double) pSv_temp[0];
               double svR_max = (double) pSv_temp[0];

               for (i = 0; i < n_temp; i++) {
                  if (pSv_temp[i] > svR_max)
                     svR_max = (double) pSv_temp[i];
                  if (pSv_temp[i] < svR_min)
                     svR_min = (double) pSv_temp[i];
               }

               if (_params->msglev >= 2)
                  cout << "   SvR_min = " << svR_min << " SvR_max = " << svR_max << endl;
               if (_params->msglev >= 1 && _params->pfout != NULL)
                  *_params->pfout << "   SvR_min = " << svR_min << " SvR_max = " <<
                     svR_max << endl;

               _FltVect fzero;

               CVector < _FltVect >::SetByZeroes (1, &fzero);

               int j;
               _FltVect aux;

               for (i = 0; i < n_temp; i++) {
                  aux = fzero;
                  for (j = 0; j < n_temp; j++) {
                     aux += pU_temp[i * n_temp + j] * palpha[j];
                  }
                  pX_temp[i] = aux / pSv_temp[i];
               }

               for (i = 0; i < n_temp; i++) {
                  aux = fzero;
                  for (j = 0; j < n_temp; j++) {
                     aux += pV_temp[j * n_temp + i] * pX_temp[j];
                  }
                  palpha[i] = aux;
               }

            }

#else
            CVector < _FltVect >::SolveR ('N', iterloc, prmatr, ncolsmax, palpha);
#endif

// Compute least squares polinomial if necessary

            if (ncoef > 1 && itgl == 1) {
               if (ncoef > iterloc)
                  ncoef = iterloc / 2;
               if (ncoef < 1) {
                  ncoef = -1;
               } else {
                  CVector < _FltVect >::Polynomial (iterloc, ncoef, phmatr, ncolsmax,
                                                    coefs);
               }
            }

         }
// Store coefs in preconditioner data if necessary

         if (itgl == 1) {
            if (myid_loc != 0)
               ncoef = 0;
            CMPIDataExchange::ExchangeArray (pcomm_loc, 'I', '+', 1, &ncoef);
            if (ncoef > 0) {
               if (myid_loc > 0) {
                  coefs.resize (ncoef + 1);
               }
               double *pcoefs = &coefs[0];
               if (myid_loc > 0) {
                  CVector < double >::SetByZeroes (ncoef, pcoefs);
               }
               CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', ncoef, pcoefs);
               this->SetNcoef (ncoef);
               vector < double >*ptr_coef_slv = this->GetCoef ();
               ptr_coef_slv->resize (ncoef + 1);
               double *pcoef_slv = &((*ptr_coef_slv)[0]);
               for (i = 0; i < ncoef; i++)
                  pcoef_slv[i] = pcoefs[i];
            }
         }
// Compute new guess to the solution

         qrdMPI_P.MvmQMPI_thr (1, palpha, ncolsP, pv, n_local);

         this->PolySlvLU (_b_use_poly, _str_mvmA, _mvm_f, _str_slvLU, _slv_f, pv, pu);

         CVector < _FltVect >::AddReplaceVector_thr (n_local, pu, px);

         (_mvm_f) (_str_mvmA, px, pv);

         CVector < _FltVect >::CopyVector_thr (n_local, _rhs, pr);

         CVector < _FltVect >::SubtractReplaceVector_thr (n_local, pv, pr);

         CVector < _FltVect >::CopyVector_thr (n_local, pr, pu);

         _FltVect resi;

         resi = CVector < _FltVect >::ScProd_thr (n_local, pu, pu);

         d_resi_norm = (double) resi;

         CMPIDataExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &d_resi_norm);

         d_resi_norm = sqrt (d_resi_norm);

         if (myid_loc == 0 && _params->msglev >= 2)
            std::cout << " Computed Log10 || Resi || = " << log10 (d_resi_norm)
               << std::endl;
         if (myid_loc == 0 && _params->msglev >= 1 && _params->pfout != NULL)
            *_params->pfout << " Computed Log10 || Resi || = " << log10 (d_resi_norm)
               << std::endl;

// Free QRD structures

         qrdGIVENS.SetNqblk (0);

         qrdMPI_P.FreeQrdMPI ();

      }

// Store solution

      CVector < _FltVect >::CopyVector_thr (n_local, px, _sol);

// Finalize timer

      CMPIDataExchange::Synchronize (pcomm_loc);

      double time1;

      time1 = CMPIDataExchange::GetWallTimeMPI ();

      _stats->dtime_iter = time1 - time0;

// Return statistics data

      _stats->niter = iter;
      _stats->res_fin = d_resi_norm;
      _stats->nmvm = _stats->niter;
      if (_b_use_poly) {
         int ncoef_loc = this->ncoef_slv;
         _stats->nmvm = _stats->niter * ncoef_loc;
      }

      return conv;

   }

// Perform iterations of the iterative scheme
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::SolveIter (SParams * _params, SStatData * _stats, void *_str_mvmA,
                             typename TMvmSlvFunc < _FltVect >::mvmA_f _mvm_f,
                             void *_str_slvLU,
                             typename TMvmSlvFunc < _FltVect >::slvLU_f _slv_f,
                             _FltVect * _rhs, _FltVect * _sol)
   {

      void *pcomm_loc = this->GetComm ();

      int myid_loc = CMPIDataExchange::GetMyid (pcomm_loc);

      bool conv;

      if (_params->ittype == 0 || _params->niter_cycle == 1) {

         conv =
            this->BiCGStab (false, _params, _stats, _str_mvmA, _mvm_f, _str_slvLU, _slv_f,
                            _rhs, _sol);

      } else if (_params->ittype == 1) {

         conv =
            this->Gmres (false, _params, _stats, _str_mvmA, _mvm_f, _str_slvLU, _slv_f,
                         _rhs, _sol);

      } else if (_params->ittype == 2) {

         int niter_max_save = _params->niter_max;

         _params->niter_max = _params->niter_cycle;

         double res_ini_temp, dtime_iter_temp;
         int niter_temp, nmvm_temp;

         conv =
            this->Gmres (false, _params, _stats, _str_mvmA, _mvm_f, _str_slvLU, _slv_f,
                         _rhs, _sol);

         _params->niter_max = niter_max_save;

         niter_temp = _stats->niter;
         nmvm_temp = _stats->nmvm;
         res_ini_temp = _stats->res_norm;
         dtime_iter_temp = _stats->dtime_iter;

         if (_stats->res_fin > _params->eps_rel * _stats->rhs_norm
             && _stats->res_fin > _params->eps_abs) {

            int niter_max_temp = niter_max_save - niter_temp;

            _params->niter_max = niter_max_temp;

            conv =
               this->BiCGStab (true, _params, _stats, _str_mvmA, _mvm_f, _str_slvLU,
                               _slv_f, _rhs, _sol);

            niter_temp += _stats->niter;
            nmvm_temp += _stats->nmvm;
            dtime_iter_temp += _stats->dtime_iter;

            _params->niter_max = niter_max_save;

            _stats->niter = niter_temp;
            _stats->nmvm = nmvm_temp;
            _stats->res_norm = res_ini_temp;
            _stats->dtime_iter = dtime_iter_temp;

         }

      } else if (_params->ittype == 3) {

         int niter_max_save = _params->niter_max;
         int niter_cycle_save = _params->niter_cycle;
         int ncoef_save = _params->ncoef;

         _params->niter_max = _params->niter_cycle;

         double res_ini_temp, dtime_iter_temp;
         int niter_temp, nmvm_temp;

         conv =
            this->Gmres (false, _params, _stats, _str_mvmA, _mvm_f, _str_slvLU, _slv_f,
                         _rhs, _sol);

         _params->niter_max = niter_max_save;

         niter_temp = _stats->niter;
         nmvm_temp = _stats->nmvm;
         res_ini_temp = _stats->res_norm;
         dtime_iter_temp = _stats->dtime_iter;

         if (_stats->res_fin > _params->eps_rel * _stats->rhs_norm
             && _stats->res_fin > _params->eps_abs) {

            int niter_max_temp = niter_max_save - niter_temp;

            _params->niter_max = niter_max_temp;
            _params->niter_cycle = _params->niter_cycle2;
            _params->ncoef = -1;

            conv =
               this->Gmres (true, _params, _stats, _str_mvmA, _mvm_f, _str_slvLU, _slv_f,
                            _rhs, _sol);

            niter_temp += _stats->niter;
            nmvm_temp += _stats->nmvm;
            dtime_iter_temp += _stats->dtime_iter;

            _params->niter_max = niter_max_save;
            _params->niter_cycle = niter_cycle_save;
            _params->ncoef = ncoef_save;

            _stats->niter = niter_temp;
            _stats->nmvm = nmvm_temp;
            _stats->res_norm = res_ini_temp;
            _stats->dtime_iter = dtime_iter_temp;

         }

      }
// Few GMRES post iterations if round off errors jumps in residual sligtly destroyed solution quality

      if (conv) {

         if (_stats->res_fin > _params->eps_rel * _stats->rhs_norm
             && _stats->res_fin > _params->eps_abs) {

            if (_params->msglev >= 2 && myid_loc == 0)
               std::cout << " Perform Post Iterations: " << std::endl;
            if (_params->msglev >= 1 && _params->pfout != NULL && myid_loc == 0)
               *_params->pfout << " Perform Post Iterations: " << std::endl;

            int niter_max_save = _params->niter_max;
            int niter_cycle_save = _params->niter_cycle;
            int ncoef_save = _params->ncoef;

            _params->niter_max = 3;
            _params->niter_cycle = 3;
            _params->ncoef = -1;

            int niter_temp = _stats->niter;
            int nmvm_temp = _stats->nmvm;
            double res_ini_temp = _stats->res_norm;
            double dtime_iter_temp = _stats->dtime_iter;

            conv =
               this->Gmres (false, _params, _stats, _str_mvmA, _mvm_f, _str_slvLU, _slv_f,
                            _rhs, _sol);

            niter_temp += _stats->niter;
            nmvm_temp += _stats->nmvm;
            dtime_iter_temp += _stats->dtime_iter;

            if (niter_temp <= niter_max_save) {
               _stats->niter = niter_temp;
               _stats->nmvm = nmvm_temp;
            }

            _stats->res_norm = res_ini_temp;
            _stats->dtime_iter = dtime_iter_temp;

            _params->niter_max = niter_max_save;
            _params->niter_cycle = niter_cycle_save;
            _params->ncoef = ncoef_save;

         }

      }

   }

// Perform polynomial preconditioning
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::PolySlvLU (bool _b_use_poly, void *_str_mvmA,
                             typename TMvmSlvFunc < _FltVect >::mvmA_f _mvm_f,
                             void *_str_slvLU,
                             typename TMvmSlvFunc < _FltVect >::slvLU_f _slv_f,
                             const _FltVect * _x, _FltVect * _px)
   {

// Fast return

      if (!_b_use_poly) {
         (_slv_f) (_str_slvLU, _x, _px);
         return;
      }
// Check and allocate work memory if necessary

      int ni_local = this->GetNi ();

      int isize = (int) (this->xwork.GetLength ());
      int isize_work = 3 * ni_local + 1;

      if (isize != isize_work) {
         this->xwork.resize (isize_work);
      }

      _FltVect *pxwork_f = this->xwork.Ptr ();

      _FltVect *pz1 = pxwork_f;
      _FltVect *pz2 = pz1 + ni_local;
      _FltVect *pz3 = pz2 + ni_local;

// Perform polynomially preconditioned computations

      int ncoef_loc = this->ncoef_slv;
      double *pcoef_loc = &(this->coef_slv[0]);

      CVector < _FltVect >::SetByZeroes_thr (ni_local, pz3);
      CVector < _FltVect >::CopyVector_thr (ni_local, _x, pz1);

      int i;
      _FltVect aux;

      for (i = 0; i < ncoef_loc; i++) {

         aux = (_FltVect) pcoef_loc[i];

         CVector < _FltVect >::UpdateVector_thr (ni_local, &aux, pz1, pz3);

         if (i < ncoef_loc - 1) {

            (_slv_f) (_str_slvLU, pz1, _px);
            (_mvm_f) (_str_mvmA, _px, pz2);

            CVector < _FltVect >::SubtractReplaceVector_thr (ni_local, pz2, pz1);

         }

      }

      (_slv_f) (_str_slvLU, pz3, _px);

   }

// Perform SlvLU parallel computations
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::SlvLU (const _FltVect * _x, _FltVect * _px)
   {
      if (!this->b_use_blksize) {
         if (this->params.prec_float == 1) {
            this->slv_float.SolveLU (_x, _px);
         } else {
            this->slv_double.SolveLU (_x, _px);
         }
      } else {

         int ni_local = 0;

         if (this->params.prec_float == 1) {
            ni_local = this->slv_float.GetNiCpu ();
         } else {
            ni_local = this->slv_double.GetNiCpu ();
         }

         int isize = (int) (this->xwork_bscl.GetLength ());
         int isize_work = 2 * ni_local + 1;

         if (isize != isize_work) {
            this->xwork_bscl.resize (isize_work);
         }

         _FltVect *pxwork_f = this->xwork_bscl.Ptr ();
         _FltVect *pxwork1_f = pxwork_f + ni_local;

         this->MvmBScl ('L', _x, pxwork_f);

         if (this->params.prec_float == 1) {
            this->slv_float.SolveLU (pxwork_f, pxwork1_f);
         } else {
            this->slv_double.SolveLU (pxwork_f, pxwork1_f);
         }

         this->MvmBScl ('U', pxwork1_f, _px);

      }
   }

// Perform mvm by part of block scaling
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::MvmBScl (char _type, const _FltVect * _x, _FltVect * _px)
   {

// Open data

      int nhblks_loc = this->nhblks;
      int *phblk2cpu = &this->hblk2cpu[0];
      int *phblk2blks = &this->hblk2blks[0];
      long long *phblks = &this->hblks[0];
      long long *pblks = &this->blks[0];

      int myid_loc = CMPIDataExchange::GetMyid (this->pcomm);

      CBMatrix < _Int, _Flt > *pmatr_arr = NULL;

      if (_type == 'L') {
         pmatr_arr = &this->sclL[0];
      } else if (_type == 'U') {
         pmatr_arr = &this->sclU[0];
      } else if (_type == 'l') {
         pmatr_arr = &this->sclLInv[0];
      } else if (_type == 'u') {
         pmatr_arr = &this->sclUInv[0];
      }

      int ihblk;

      int ibs_vect = 0;

      int blksizeloc = this->blksize_bscl;

      for (ihblk = 0; ihblk < nhblks_loc; ihblk++) {
         if (phblk2cpu[ihblk] == myid_loc) {

            int niblk = phblk2blks[ihblk + 1] - phblk2blks[ihblk];
            int ibegblk = phblk2blks[ihblk];

            CMatrix < _Int, _Flt > *pASub = pmatr_arr[ihblk].GetASubArr ();

            {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
               for (int ipar = 0; ipar < niblk; ipar++) {
                  int iblkgl = ibegblk + ipar;
                  int nlistloc = (int) (pblks[iblkgl + 1] - pblks[iblkgl]);
                  int nlistloc_cnd = nlistloc / blksizeloc;
                  int ishift = (int) (pblks[ibegblk + ipar] - pblks[ibegblk]);

                  vector < _Flt > *pa_matr = pASub[ipar].GetA ();

                  CMvmSlv_impl < _Int, _Flt, _FltVect >::MvmDiagBlksize (blksizeloc,
                                                                         nlistloc_cnd,
                                                                         *pa_matr,
                                                                         _x + ibs_vect +
                                                                         ishift,
                                                                         _px + ibs_vect +
                                                                         ishift);

               }

            }

            ibs_vect += (int) (phblks[ihblk + 1] - phblks[ihblk]);

         }
      }

   }

// Reorder rhs and solution from/to initial ordering
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::ReorderVectorDataIni (char _dir, _FltVect * _rhs, _FltVect * _sol)
   {

// Open data

      int nhblks_loc = this->nhblks;
      int *phblk2cpu = &this->hblk2cpu[0];
      long long *phblks = &this->hblks[0];
      CVectorData < int >*porder_ini = &this->order_ini[0];

      int myid_loc = CMPIDataExchange::GetMyid (this->pcomm);

// Reorder data

      int nimax = 0;

      int niloc, i;

      for (i = 0; i < nhblks_loc; i++) {
         niloc = (int) (phblks[i + 1] - phblks[i]);
         if (niloc > nimax)
            nimax = niloc;
      }

      CVectorData < _FltVect > ztemp (nimax);
      _FltVect *pztemp = ztemp.Ptr ();

      int ibs = 0;

      for (i = 0; i < nhblks_loc; i++) {
         niloc = (int) (phblks[i + 1] - phblks[i]);
         if (phblk2cpu[i] == myid_loc) {
            int *pord_ini = porder_ini[i].Ptr ();
            if (_dir == 'D') {
               CVector < _FltVect >::OrderVector_thr (niloc, pord_ini, _rhs + ibs,
                                                      pztemp);
               CVector < _FltVect >::CopyVector_thr (niloc, pztemp, _rhs + ibs);
               CVector < _FltVect >::OrderVector_thr (niloc, pord_ini, _sol + ibs,
                                                      pztemp);
               CVector < _FltVect >::CopyVector_thr (niloc, pztemp, _sol + ibs);
            } else {
               CVector < _FltVect >::InvOrderVector_thr (niloc, pord_ini, _rhs + ibs,
                                                         pztemp);
               CVector < _FltVect >::CopyVector_thr (niloc, pztemp, _rhs + ibs);
               CVector < _FltVect >::InvOrderVector_thr (niloc, pord_ini, _sol + ibs,
                                                         pztemp);
               CVector < _FltVect >::CopyVector_thr (niloc, pztemp, _sol + ibs);
            }
            ibs += niloc;
         }
      }

   }

// Reorder rhs and solution from/to wells ordering
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverThreads < _Int, _Flt,
      _FltVect >::ReorderVectorDataWells (char _dir, _FltVect * _rhs, _FltVect * _sol,
                                          CVectorData < _FltVect > &_rhs_ord,
                                          CVectorData < _FltVect > &_sol_ord)
   {

      int nproc_loc = CMPIDataExchange::GetNproc (this->pcomm);
      int myid_loc = CMPIDataExchange::GetMyid (this->pcomm);

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

// Open data

      int nhblks_ini_loc = this->nhblks_ini;
      int *phblk2cpu_ini = &this->hblk2cpu_ini[0];
      long long *phblks_ini = &this->hblks_ini[0];
      int *porder2ind_wells = this->order2ind_wells.Ptr ();

      int nhblks_loc = this->nhblks;
      int *phblk2cpu = &this->hblk2cpu[0];
      long long *phblks = &this->hblks[0];

// Create base reference data

      CVectorData < long long >ibs_ord_ini (nhblks_ini_loc + 1);
      CVectorData < long long >ibs_ord (nhblks_loc + 1);

      long long *pibs_ord_ini = ibs_ord_ini.Ptr ();
      long long *pibs_ord = ibs_ord.Ptr ();

      long long nz_ord_ini = 0;
      long long nz_ord = 0;

      {
         int ihblk;
         for (ihblk = 0; ihblk < nhblks_ini_loc; ihblk++) {
            if (phblk2cpu_ini[ihblk] == myid_loc) {
               pibs_ord_ini[ihblk] = nz_ord_ini;
               nz_ord_ini += (phblks_ini[ihblk + 1] - phblks_ini[ihblk]);
            } else {
               pibs_ord_ini[ihblk] = -1;
            }
         }
         for (ihblk = 0; ihblk < nhblks_loc; ihblk++) {
            if (phblk2cpu[ihblk] == myid_loc) {
               pibs_ord[ihblk] = nz_ord;
               nz_ord += (phblks[ihblk + 1] - phblks[ihblk]);
            } else {
               pibs_ord[ihblk] = -1;
            }
         }

      }

// Allocate vector data

      if (_dir == 'd' || _dir == 'D') {
         _rhs_ord.resize ((int) nz_ord);
         _sol_ord.resize ((int) nz_ord);
      }

      _FltVect *p_rhs_ord = _rhs_ord.Ptr ();
      _FltVect *p_sol_ord = _sol_ord.Ptr ();

// Split initial data into large parts for threads

      int nparts = n_thr;
      int isize_part = (int) (nz_ord_ini / nparts);
      CVectorData < int >parts (nparts + 1);
      int *pparts = parts.Ptr ();

      {
         int i;
         for (i = 0; i < nparts; i++)
            pparts[i] = i * isize_part;
         pparts[nparts] = (int) nz_ord_ini;
      }

// Create lists of sends/recvs to/from each cpu separately for parts

      CVectorData < int >nblks_cpu_arr (nparts);
      CVectorData < CVectorData < int > >list_cpu_arr (nparts);
      CVectorData < CVectorData < int > >blks_cpu_arr (nparts);
      CVectorData < int >index3_arr ((int) nz_ord_ini * 3);
      CVectorData < CVectorData < CMatrix < int, _FltVect > > >send_blks_arr (nparts);

      int *pnblks_cpu_arr = nblks_cpu_arr.Ptr ();
      CVectorData < int >*plist_cpu_arr = list_cpu_arr.Ptr ();
      CVectorData < int >*pblks_cpu_arr = blks_cpu_arr.Ptr ();
      int *pindex3_arr = index3_arr.Ptr ();
      CVectorData < CMatrix < int, _FltVect > >*psend_blks_arr = send_blks_arr.Ptr ();

      {

         CVectorData < int >imaskcpu_threads (3 * nproc_loc * n_thr);
         CVectorData < int >icycle_threads (n_thr);
         int *pimaskcpu_threads = imaskcpu_threads.Ptr ();
         int *picycle_threads = icycle_threads.Ptr ();
         {
            for (int i = 0; i < n_thr; i++)
               picycle_threads[i] = -1;
         }

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nparts; ipar++) {

            int my_thr = 0;
#ifdef USE_THREADS
            my_thr = omp_get_thread_num ();
#endif

            int *pimaskcpu_th = pimaskcpu_threads + 3 * nproc_loc * my_thr;
            int *plistcpu_th = pimaskcpu_th + nproc_loc;
            int *pindcpu_th = plistcpu_th + nproc_loc;
            int icycle_th = picycle_threads[my_thr];

            int i;

            if (icycle_th == -1) {
               for (i = 0; i < nproc_loc; i++)
                  pimaskcpu_th[i] = -1;
            }

            icycle_th++;

            int ibeg = pparts[ipar];
            int iend = pparts[ipar + 1] - 1;

            int j, jhblk, jcpu;

            int nlistcpu = 0;

            for (j = ibeg; j <= iend; j++) {
               jhblk = porder2ind_wells[j * 2 + 1];
               jcpu = phblk2cpu[jhblk];
               if (pimaskcpu_th[jcpu] != icycle_th) {
                  plistcpu_th[nlistcpu] = jcpu;
                  nlistcpu++;
                  pimaskcpu_th[jcpu] = icycle_th;
               }
            }

            sort (plistcpu_th, plistcpu_th + nlistcpu);

            for (j = 0; j < nlistcpu; j++) {
               jcpu = plistcpu_th[j];
               pindcpu_th[jcpu] = j;
            }

            pnblks_cpu_arr[ipar] = nlistcpu;
            plist_cpu_arr[ipar].resize (nlistcpu);
            pblks_cpu_arr[ipar].resize (nlistcpu + 1);
            psend_blks_arr[ipar].resize (nlistcpu);

            int *pplist_cpu_arr = plist_cpu_arr[ipar].Ptr ();
            int *ppblks_cpu_arr = pblks_cpu_arr[ipar].Ptr ();
            CMatrix < int, _FltVect > *ppsend_blks_arr = psend_blks_arr[ipar].Ptr ();

            for (j = 0; j < nlistcpu; j++)
               pplist_cpu_arr[j] = plistcpu_th[j];
            for (j = 0; j <= nlistcpu; j++)
               ppblks_cpu_arr[j] = 0;
            int ind, jj;
            for (j = ibeg; j <= iend; j++) {
               jj = porder2ind_wells[j * 2];
               jhblk = porder2ind_wells[j * 2 + 1];
               jcpu = phblk2cpu[jhblk];
               ind = pindcpu_th[jcpu];
               ppblks_cpu_arr[ind + 1]++;
            }
            for (j = 0; j < nlistcpu; j++)
               ppblks_cpu_arr[j + 1] += ppblks_cpu_arr[j];
            CVectorData < int >iptr (nlistcpu);
            int *piptr = iptr.Ptr ();
            for (j = 0; j < nlistcpu; j++)
               piptr[j] = ppblks_cpu_arr[j];
            int k;
            for (j = ibeg; j <= iend; j++) {
               jj = porder2ind_wells[j * 2];
               jhblk = porder2ind_wells[j * 2 + 1];
               jcpu = phblk2cpu[jhblk];
               ind = pindcpu_th[jcpu];
               k = ibeg + piptr[ind];
               pindex3_arr[k * 3] = j;
               pindex3_arr[k * 3 + 1] = jj;
               pindex3_arr[k * 3 + 2] = jhblk;
               piptr[ind]++;
            }

            int iibeg, iiend, nii_curr;

            for (i = 0; i < nlistcpu; i++) {
               iibeg = ibeg + ppblks_cpu_arr[i];
               iiend = ibeg + ppblks_cpu_arr[i + 1] - 1;
               nii_curr = iiend - iibeg + 1;
               if (_dir == 'd' || _dir == 'D') {
                  ppsend_blks_arr[i].ResizeAndSetAll (0, nii_curr * 2, 0, 0,
                                                      nii_curr * 2);
               } else {
                  ppsend_blks_arr[i].ResizeAndSetAll (0, nii_curr * 3, 0, 0, 0);
               }
               int *plist2_temp = ppsend_blks_arr[i].GetList2Arr ();
               _FltVect *pa_temp = ppsend_blks_arr[i].GetAArr ();
               if (_dir == 'd' || _dir == 'D') {
                  for (j = 0; j < nii_curr; j++) {
                     k = iibeg + j;
                     ind = pindex3_arr[k * 3];
                     plist2_temp[j * 2] = pindex3_arr[k * 3 + 1];
                     plist2_temp[j * 2 + 1] = pindex3_arr[k * 3 + 2];
                     pa_temp[j * 2] = _rhs[ind];
                     pa_temp[j * 2 + 1] = _sol[ind];
                  }
               } else {
                  for (j = 0; j < nii_curr; j++) {
                     k = iibeg + j;
                     plist2_temp[j * 3] = pindex3_arr[k * 3];
                     plist2_temp[j * 3 + 1] = pindex3_arr[k * 3 + 1];
                     plist2_temp[j * 3 + 2] = pindex3_arr[k * 3 + 2];
                  }
               }
            }

            picycle_threads[my_thr] = icycle_th;
         }

      }

// Prepare exchanges

      CVectorData < int >imaskcpu (3 * nproc_loc);

      int *pimaskcpu = imaskcpu.Ptr ();
      int *plistcpu = pimaskcpu + nproc_loc;
      int *pindcpu = plistcpu + nproc_loc;

      int nlistcpu = 0;

      {
         int i;
         for (i = 0; i < nproc_loc; i++)
            pimaskcpu[i] = -1;
         int j, jcpu;
         for (i = 0; i < nparts; i++) {
            int *pplist_cpu_arr = plist_cpu_arr[i].Ptr ();
            for (j = 0; j < pnblks_cpu_arr[i]; j++) {
               jcpu = pplist_cpu_arr[j];
               if (pimaskcpu[jcpu] < 0) {
                  plistcpu[nlistcpu] = jcpu;
                  nlistcpu++;
                  pimaskcpu[jcpu] = 1;
               }
            }
         }
         sort (plistcpu, plistcpu + nlistcpu);
         for (i = 0; i < nlistcpu; i++) {
            jcpu = plistcpu[i];
            pindcpu[jcpu] = i;
            pimaskcpu[jcpu] = 0;
         }
         for (i = 0; i < nparts; i++) {
            int *pplist_cpu_arr = plist_cpu_arr[i].Ptr ();
            for (j = 0; j < pnblks_cpu_arr[i]; j++) {
               jcpu = pplist_cpu_arr[j];
               pimaskcpu[jcpu]++;
            }
         }
      }

      CVectorData < CBMatrix < int, _FltVect > >hblk_send (nlistcpu);
      CBMatrix < int, _FltVect > *phblk_send = hblk_send.Ptr ();

      {
         int i, jcpu;
         for (i = 0; i < nlistcpu; i++) {
            jcpu = plistcpu[i];
            phblk_send[i].ResizeASub (pimaskcpu[jcpu]);
            phblk_send[i].SetNzblk (pimaskcpu[jcpu]);
            pimaskcpu[jcpu] = 0;
         }
         int j, k, ind;
         for (i = 0; i < nparts; i++) {
            int *pplist_cpu_arr = plist_cpu_arr[i].Ptr ();
            CMatrix < int, _FltVect > *ppsend_blks_arr = psend_blks_arr[i].Ptr ();
            for (j = 0; j < pnblks_cpu_arr[i]; j++) {
               jcpu = pplist_cpu_arr[j];
               ind = pindcpu[jcpu];
               CMatrix < int, _FltVect > *pASub = phblk_send[ind].GetASubArr ();
               k = pimaskcpu[jcpu];
               pASub[k].ReplaceFree (ppsend_blks_arr[j]);
               pimaskcpu[jcpu]++;
            }
         }
      }

// Pack send data

      vector < int >CpuIDSend (nlistcpu);
      vector < vector < char > >ObjSend (nlistcpu);

      int *pCpuIDSend = NULL;
      vector < char >*pObjSend = NULL;

      if (nlistcpu > 0) {
         pCpuIDSend = &CpuIDSend[0];
         pObjSend = &ObjSend[0];
      }

      {

         long long isize;
         char *pobj;

         int i;

         for (i = 0; i < nlistcpu; i++) {
            pCpuIDSend[i] = plistcpu[i];
            isize = phblk_send[i].GetPackedSize ();
            pObjSend[i].resize ((size_t) isize);
            pobj = &(pObjSend[i][0]);
            phblk_send[i].FillPacked_thr (isize, pobj);
            phblk_send[i].Clean ();
         }

      }

// Exchange

      vector < int >CpuIDRecv;
      vector < vector < char > >ObjRecv;

      CMPIDataExchange::DataExchange (this->pcomm, CpuIDSend, ObjSend, CpuIDRecv,
                                      ObjRecv);

// Free send data

      {
         vector < int >CpuIDSend_temp;
         vector < vector < char > >ObjSend_temp;
         CpuIDSend.swap (CpuIDSend_temp);
         ObjSend.swap (ObjSend_temp);
      }

// Unpack receive data

      int nrecv = (int) CpuIDRecv.size ();

      vector < char >*pObjRecv = NULL;

      if (nrecv > 0) {
         pObjRecv = &ObjRecv[0];
      }

      vector < CBMatrix < int, _FltVect > >hblk_recv (nrecv + 1);

      CBMatrix < int, _FltVect > *phblk_recv = &hblk_recv[0];

      {

         long long isize;
         char *pobj;

         int i;

         for (i = 0; i < nrecv; i++) {
            isize = (long long) pObjRecv[i].size ();
            pobj = &(pObjRecv[i][0]);
            phblk_recv[i].UnPack_thr (isize, pobj);
         }

      }

// Free recv data

      {
         vector < vector < char > >ObjRecv_temp;
         ObjRecv.swap (ObjRecv_temp);
      }

// Complete forward computations

      if (_dir == 'd' || _dir == 'D') {
         int irecv;
         for (irecv = 0; irecv < nrecv; irecv++) {
            int nzblk_loc = phblk_recv[irecv].GetNzblk ();
            CMatrix < int, _FltVect > *pASub = phblk_recv[irecv].GetASubArr ();
#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < nzblk_loc; ipar++) {
               int i, jj, jhblk, ibs;
               int nlist_temp = pASub[ipar].GetNlist2 () / 2;
               int *plist2_temp = pASub[ipar].GetList2Arr ();
               _FltVect *pa_temp = pASub[ipar].GetAArr ();
               for (i = 0; i < nlist_temp; i++) {
                  jj = plist2_temp[i * 2];
                  jhblk = plist2_temp[i * 2 + 1];
                  ibs = (int) pibs_ord[jhblk];
                  p_rhs_ord[ibs + jj] = pa_temp[i * 2];
                  p_sol_ord[ibs + jj] = pa_temp[i * 2 + 1];
               }
            }

         }

         return;

      } else {

// Store vector data

         int irecv;
         for (irecv = 0; irecv < nrecv; irecv++) {
            int nzblk_loc = phblk_recv[irecv].GetNzblk ();
            CMatrix < int, _FltVect > *pASub = phblk_recv[irecv].GetASubArr ();
#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < nzblk_loc; ipar++) {
               int i, jj, jhblk, ibs;
               int nlist_temp = pASub[ipar].GetNlist2 () / 3;
               int *plist2_temp = pASub[ipar].GetList2Arr ();
               pASub[ipar].ResizeA (nlist_temp);
               pASub[ipar].SetNza (nlist_temp);
               _FltVect *pa_temp = pASub[ipar].GetAArr ();
               for (i = 0; i < nlist_temp; i++) {
                  jj = plist2_temp[i * 3 + 1];
                  jhblk = plist2_temp[i * 3 + 2];
                  ibs = (int) pibs_ord[jhblk];
                  pa_temp[i] = p_sol_ord[ibs + jj];
               }
            }

         }

// Pack data back

         ObjRecv.resize (nrecv);

         vector < char >*pObjRecv = NULL;

         if (nrecv > 0) {
            pObjRecv = &ObjRecv[0];
         }

         {

            long long isize;
            char *pobj;

            int i;

            for (i = 0; i < nrecv; i++) {
               isize = phblk_recv[i].GetPackedSize ();
               pObjRecv[i].resize ((size_t) isize);
               pobj = &(pObjRecv[i][0]);
               phblk_recv[i].FillPacked_thr (isize, pobj);
               phblk_recv[i].Clean ();
            }

         }

// Exchange data back

         CMPIDataExchange::DataExchange (this->pcomm, CpuIDRecv, ObjRecv, CpuIDSend,
                                         ObjSend);

// Free send data

         {
            vector < int >CpuIDSend_temp;
            vector < vector < char > >ObjSend_temp;
            CpuIDRecv.swap (CpuIDSend_temp);
            ObjRecv.swap (ObjSend_temp);
         }

// Unpack receive data

         nrecv = (int) CpuIDSend.size ();

         pObjRecv = NULL;

         if (nrecv > 0) {
            pObjRecv = &ObjSend[0];
         }

         hblk_recv.resize (nrecv + 1);

         phblk_recv = &hblk_recv[0];

         {

            long long isize;
            char *pobj;

            int i;

            for (i = 0; i < nrecv; i++) {
               isize = (long long) pObjRecv[i].size ();
               pobj = &(pObjRecv[i][0]);
               phblk_recv[i].UnPack_thr (isize, pobj);
            }

         }

// Store vector data

         for (irecv = 0; irecv < nrecv; irecv++) {
            int nzblk_loc = phblk_recv[irecv].GetNzblk ();
            CMatrix < int, _FltVect > *pASub = phblk_recv[irecv].GetASubArr ();
#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < nzblk_loc; ipar++) {
               int i, ind;
               int nlist_temp = pASub[ipar].GetNlist2 () / 3;
               int *plist2_temp = pASub[ipar].GetList2Arr ();
               _FltVect *pa_temp = pASub[ipar].GetAArr ();
               for (i = 0; i < nlist_temp; i++) {
                  ind = plist2_temp[i * 3];
                  _sol[ind] = pa_temp[i];
               }
            }

         }

      }

   }

   template class CFctThreads < long long, double >;
   template class CFctThreads < int, double >;
   template class CFctThreads < long long, float >;
   template class CFctThreads < int, float >;
   template class CMvmParThreads < long long, double, double >;
   template class CMvmParThreads < int, double, double >;
   template class CMvmParThreads < long long, float, double >;
   template class CMvmParThreads < int, float, double >;
   template class CSlvParThreads < long long, double, double >;
   template class CSlvParThreads < int, double, double >;
   template class CSlvParThreads < long long, float, double >;
   template class CSlvParThreads < int, float, double >;
   template class CK3D_SolverThreads < long long, double, double >;
   template class CK3D_SolverThreads < int, double, double >;
   template class CK3D_SolverThreads < long long, float, double >;
   template class CK3D_SolverThreads < int, float, double >;

}                               // namespace k3d

#endif
