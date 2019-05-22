//------------------------------------------------------------------------------------------------
// File: k3d_block.cxx
//------------------------------------------------------------------------------------------------

#ifndef __K3D_BLOCK_CXX
#define __K3D_BLOCK_CXX

#include "k3d_block.hxx"

//#include "k3d_spbasis.h"

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
// Multiply C = A*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MM_BxB (int _b_size,
                                                                       const _Flt *
                                                                       _A_blk,
                                                                       const _Flt *
                                                                       _B_blk,
                                                                       _Flt * _AxB_blk)
   {
      static const _Flt fzero = (_Flt) 0.0e0;
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++)
      {
         for (j = 0; j < _b_size; j++)
         {
            aux = fzero;
            for (k = 0; k < _b_size; k++)
            {
               aux += _A_blk[k * _b_size + i] * _B_blk[j * _b_size + k];
            }
            _AxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C += A*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MM_Add_BxB (int _b_size,
                                                                           const _Flt *
                                                                           _A_blk,
                                                                           const _Flt *
                                                                           _B_blk,
                                                                           _Flt *
                                                                           _AxB_blk)
   {
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            aux = _AxB_blk[j * _b_size + i];
            for (k = 0; k < _b_size; k++) {
               aux += _A_blk[k * _b_size + i] * _B_blk[j * _b_size + k];
            }
            _AxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C -= A*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MM_Sub_BxB (int _b_size,
                                                                           const _Flt *
                                                                           _A_blk,
                                                                           const _Flt *
                                                                           _B_blk,
                                                                           _Flt *
                                                                           _AxB_blk)
   {
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            aux = _AxB_blk[j * _b_size + i];
            for (k = 0; k < _b_size; k++) {
               aux -= _A_blk[k * _b_size + i] * _B_blk[j * _b_size + k];
            }
            _AxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C = At*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MtM_BxB (int _b_size,
                                                                        const _Flt *
                                                                        _A_blk,
                                                                        const _Flt *
                                                                        _B_blk,
                                                                        _Flt * _AtxB_blk)
   {
      static const _Flt fzero = (_Flt) 0.0e0;
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            aux = fzero;
            for (k = 0; k < _b_size; k++) {
               aux += _A_blk[i * _b_size + k] * _B_blk[j * _b_size + k];
            }
            _AtxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C += At*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MtM_Add_BxB (int _b_size,
                                                                            const _Flt *
                                                                            _A_blk,
                                                                            const _Flt *
                                                                            _B_blk,
                                                                            _Flt *
                                                                            _AtxB_blk)
   {
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            aux = _AtxB_blk[j * _b_size + i];
            for (k = 0; k < _b_size; k++) {
               aux += _A_blk[i * _b_size + k] * _B_blk[j * _b_size + k];
            }
            _AtxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C -= At*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MtM_Sub_BxB (int _b_size,
                                                                            const _Flt *
                                                                            _A_blk,
                                                                            const _Flt *
                                                                            _B_blk,
                                                                            _Flt *
                                                                            _AtxB_blk)
   {
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            aux = _AtxB_blk[j * _b_size + i];
            for (k = 0; k < _b_size; k++) {
               aux -= _A_blk[i * _b_size + k] * _B_blk[j * _b_size + k];
            }
            _AtxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C = A*Bt
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MMt_BxB (int _b_size,
                                                                        const _Flt *
                                                                        _A_blk,
                                                                        const _Flt *
                                                                        _B_blk,
                                                                        _Flt * _AxBt_blk)
   {
      static const _Flt fzero = (_Flt) 0.0e0;
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            aux = fzero;
            for (k = 0; k < _b_size; k++) {
               aux += _A_blk[k * _b_size + i] * _B_blk[k * _b_size + j];
            }
            _AxBt_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C += A*Bt
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MMt_Add_BxB (int _b_size,
                                                                            const _Flt *
                                                                            _A_blk,
                                                                            const _Flt *
                                                                            _B_blk,
                                                                            _Flt *
                                                                            _AxBt_blk)
   {
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            aux = _AxBt_blk[j * _b_size + i];
            for (k = 0; k < _b_size; k++) {
               aux += _A_blk[k * _b_size + i] * _B_blk[k * _b_size + j];
            }
            _AxBt_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C -= A*Bt
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MMt_Sub_BxB (int _b_size,
                                                                            const _Flt *
                                                                            _A_blk,
                                                                            const _Flt *
                                                                            _B_blk,
                                                                            _Flt *
                                                                            _AxBt_blk)
   {
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            aux = _AxBt_blk[j * _b_size + i];
            for (k = 0; k < _b_size; k++) {
               aux -= _A_blk[k * _b_size + i] * _B_blk[k * _b_size + j];
            }
            _AxBt_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C = A*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MM_BxB (int _b_size,
                                                                       int _ncol,
                                                                       const _Flt *
                                                                       _A_blk,
                                                                       const _Flt *
                                                                       _B_blk,
                                                                       _Flt * _AxB_blk)
   {
      static const _Flt fzero = (_Flt) 0.0e0;
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _ncol; j++) {
            aux = fzero;
            for (k = 0; k < _b_size; k++) {
               aux += _A_blk[k * _b_size + i] * _B_blk[j * _b_size + k];
            }
            _AxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C += A*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MM_Add_BxB (int _b_size,
                                                                           int _ncol,
                                                                           const _Flt *
                                                                           _A_blk,
                                                                           const _Flt *
                                                                           _B_blk,
                                                                           _Flt *
                                                                           _AxB_blk)
   {
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _ncol; j++) {
            aux = _AxB_blk[j * _b_size + i];
            for (k = 0; k < _b_size; k++) {
               aux += _A_blk[k * _b_size + i] * _B_blk[j * _b_size + k];
            }
            _AxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C -= A*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MM_Sub_BxB (int _b_size,
                                                                           int _ncol,
                                                                           const _Flt *
                                                                           _A_blk,
                                                                           const _Flt *
                                                                           _B_blk,
                                                                           _Flt *
                                                                           _AxB_blk)
   {
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _ncol; j++) {
            aux = _AxB_blk[j * _b_size + i];
            for (k = 0; k < _b_size; k++) {
               aux -= _A_blk[k * _b_size + i] * _B_blk[j * _b_size + k];
            }
            _AxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C = At*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MtM_BxB (int _b_size,
                                                                        int _ncol,
                                                                        const _Flt *
                                                                        _A_blk,
                                                                        const _Flt *
                                                                        _B_blk,
                                                                        _Flt * _AtxB_blk)
   {
      static const _Flt fzero = (_Flt) 0.0e0;
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _ncol; j++) {
            aux = fzero;
            for (k = 0; k < _b_size; k++) {
               aux += _A_blk[i * _b_size + k] * _B_blk[j * _b_size + k];
            }
            _AtxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C += At*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MtM_Add_BxB (int _b_size,
                                                                            int _ncol,
                                                                            const _Flt *
                                                                            _A_blk,
                                                                            const _Flt *
                                                                            _B_blk,
                                                                            _Flt *
                                                                            _AtxB_blk)
   {
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _ncol; j++) {
            aux = _AtxB_blk[j * _b_size + i];
            for (k = 0; k < _b_size; k++) {
               aux += _A_blk[i * _b_size + k] * _B_blk[j * _b_size + k];
            }
            _AtxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C -= At*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MtM_Sub_BxB (int _b_size,
                                                                            int _ncol,
                                                                            const _Flt *
                                                                            _A_blk,
                                                                            const _Flt *
                                                                            _B_blk,
                                                                            _Flt *
                                                                            _AtxB_blk)
   {
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _ncol; j++) {
            aux = _AtxB_blk[j * _b_size + i];
            for (k = 0; k < _b_size; k++) {
               aux -= _A_blk[i * _b_size + k] * _B_blk[j * _b_size + k];
            }
            _AtxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C = A*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MM_BxB (int _b_size,
                                                                       int _b_size2,
                                                                       int _ncol,
                                                                       const _Flt *
                                                                       _A_blk,
                                                                       const _Flt *
                                                                       _B_blk,
                                                                       _Flt * _AxB_blk)
   {
      static const _Flt fzero = (_Flt) 0.0e0;
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _ncol; j++) {
            aux = fzero;
            for (k = 0; k < _b_size2; k++) {
               aux += _A_blk[k * _b_size + i] * _B_blk[j * _b_size2 + k];
            }
            _AxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C += A*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MM_Add_BxB (int _b_size,
                                                                           int _b_size2,
                                                                           int _ncol,
                                                                           const _Flt *
                                                                           _A_blk,
                                                                           const _Flt *
                                                                           _B_blk,
                                                                           _Flt *
                                                                           _AxB_blk)
   {
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _ncol; j++) {
            aux = _AxB_blk[j * _b_size + i];
            for (k = 0; k < _b_size2; k++) {
               aux += _A_blk[k * _b_size + i] * _B_blk[j * _b_size2 + k];
            }
            _AxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Multiply C -= A*B
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::MM_Sub_BxB (int _b_size,
                                                                           int _b_size2,
                                                                           int _ncol,
                                                                           const _Flt *
                                                                           _A_blk,
                                                                           const _Flt *
                                                                           _B_blk,
                                                                           _Flt *
                                                                           _AxB_blk)
   {
      int i, j, k;
      _Flt aux;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _ncol; j++) {
            aux = _AxB_blk[j * _b_size + i];
            for (k = 0; k < _b_size2; k++) {
               aux -= _A_blk[k * _b_size + i] * _B_blk[j * _b_size2 + k];
            }
            _AxB_blk[j * _b_size + i] = aux;
         }
      }
   }

//
// Transpose block
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::Transp_BxB (int _b_size,
                                                                           const _Flt *
                                                                           _A_blk,
                                                                           _Flt * _At_blk)
   {
      int i, j;
      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            _At_blk[i * _b_size + j] = _A_blk[j * _b_size + i];
         }
      }
   }

//
// Compute triangular fct
//========================================================================================
   template < typename _Flt > int CBlock_BxB_traits <
      _Flt >::FctSymm_BxB (double _piv_min, int _b_size, const _Flt * _AU_blk,
                           _Flt * _U_blk, double &_dia_min, double &_dia_max)
   {

      int n_modif = 0;

      static const _Flt fzero = (_Flt) 0.0e0;
      static const _Flt fone = (_Flt) 1.0e0;

      _Flt auxU;
      int i, j, k;
      double daux;

      for (j = 0; j < _b_size; j++) {
         for (i = 0; i <= j; i++) {
            _U_blk[i + j * _b_size] = _AU_blk[i + j * _b_size];
         }
      }
      for (j = 0; j < _b_size; j++) {
         for (i = j + 1; i < _b_size; i++) {
            _U_blk[i + j * _b_size] = fzero;
         }
      }
      for (i = 0; i < _b_size; i++) {
         for (k = 0; k < i; k++) {
            for (j = i; j < _b_size; j++) {
               _U_blk[i + j * _b_size] -=
                  _U_blk[k + i * _b_size] * _U_blk[k + j * _b_size];
            }
         }
         auxU = _U_blk[i + i * _b_size];
         daux = (double) auxU;
         if (i == 0 || daux < _dia_min)
            _dia_min = daux;
         if (i == 0 || daux > _dia_max)
            _dia_max = daux;
         if (auxU < _piv_min) {
            auxU = (_Flt) _piv_min;
            n_modif++;
         }
         if (auxU < 0.0e0) {
            cout << " CBlock_BxB_traits < _Flt >::FctSymm_BxB: error: pivot is negative !"
               << endl;
            throw " CBlock_BxB_traits < _Flt >::FctSymm_BxB: error: pivot is negative !";
         }
         auxU = sqrt (auxU);
         _U_blk[i + i * _b_size] = auxU;
         auxU = fone / auxU;
         for (j = i + 1; j < _b_size; j++) {
            _U_blk[i + j * _b_size] *= auxU;
         }
      }

      return n_modif;

   }

//
// Compute triangular fct
//========================================================================================
   template < typename _Flt > int CBlock_BxB_traits < _Flt >::Fct_BxB (double _piv_min,
                                                                       int _b_size,
                                                                       const _Flt *
                                                                       _A_blk,
                                                                       _Flt * _Lt_blk,
                                                                       _Flt * _U_blk,
                                                                       double &_dia_min,
                                                                       double &_dia_max)
   {

      int n_modif = 0;

      static const _Flt fzero = (_Flt) 0.0e0;
      static const _Flt fone = (_Flt) 1.0e0;

      _Flt auxL, auxU;
      bool is_neg;
      int i, j, k;
      double daux;

      for (j = 0; j < _b_size; j++) {
         for (i = 0; i <= j; i++) {
            _Lt_blk[i + j * _b_size] = _A_blk[j + i * _b_size];
            _U_blk[i + j * _b_size] = _A_blk[i + j * _b_size];
         }
      }
      for (j = 0; j < _b_size; j++) {
         for (i = j + 1; i < _b_size; i++) {
            _Lt_blk[i + j * _b_size] = fzero;
            _U_blk[i + j * _b_size] = fzero;
         }
      }
      for (i = 0; i < _b_size; i++) {
         for (k = 0; k < i; k++) {
            for (j = i; j < _b_size; j++) {
               _Lt_blk[i + j * _b_size] -=
                  _U_blk[k + i * _b_size] * _Lt_blk[k + j * _b_size];
               _U_blk[i + j * _b_size] -=
                  _Lt_blk[k + i * _b_size] * _U_blk[k + j * _b_size];
            }
         }
         auxU = _U_blk[i + i * _b_size];
         daux = (double) auxU;
         if (daux < 0)
            daux = -daux;
         if (i == 0 || daux < _dia_min)
            _dia_min = daux;
         if (i == 0 || daux > _dia_max)
            _dia_max = daux;
         is_neg = false;
         if (auxU < fzero) {
            auxU = -auxU;
            is_neg = true;
         }
         if (auxU < _piv_min) {
            auxU = (_Flt) _piv_min;
            n_modif++;
         }
         auxL = sqrt (auxU);
         auxU = auxL;
         if (is_neg)
            auxL = -auxL;
         _Lt_blk[i + i * _b_size] = auxL;
         _U_blk[i + i * _b_size] = auxU;
         auxU = fone / auxU;
         auxL = auxU;
         if (is_neg)
            auxL = -auxL;
         for (j = i + 1; j < _b_size; j++) {
            _Lt_blk[i + j * _b_size] *= auxU;
            _U_blk[i + j * _b_size] *= auxL;
         }
      }

      return n_modif;

   }

//
// Inverse dense square matrix
//========================================================================================
   template < typename _Flt > int CBlock_BxB_traits <
      _Flt >::InvA_Svd_BxB (double _piv_min, int _b_size, _Flt * _A_blk,
                            _Flt * _A_inv_blk, _Flt * _Work, double *_dWork,
                            double &_dia_min, double &_dia_max)
   {

      int b_2 = _b_size * _b_size;

      _Flt fone = (_Flt) 1.0e0;

// Svd

      _Flt *pSv = _Work;
      _Flt *pSvInv = pSv + _b_size;
      _Flt *pUMatr = pSvInv + _b_size;
      _Flt *pVMatr = pUMatr + b_2;
      _Flt *pUSMatr = pVMatr + b_2;
      _Flt *pVSMatr = pUSMatr + b_2;

      CBlock_BxB_traits < _Flt >::ComputeSvd_BxB (_b_size, _A_blk, pSv, pUMatr, pVMatr,
                                                  _dWork);

// Modif Sv

      int n_modif = 0;

      int i, j;

      _dia_min = pSv[0];
      _dia_max = pSv[0];

      for (i = 0; i < _b_size; i++) {
         if (pSv[i] < _dia_min)
            _dia_min = (double) pSv[i];
         if (pSv[i] > _dia_max)
            _dia_max = (double) pSv[i];
         if (pSv[i] < _piv_min) {
            n_modif++;
            pSv[i] = (_Flt) _piv_min;
         }
         pSv[i] = sqrt (pSv[i]);
         pSvInv[i] = fone / pSv[i];
      }

// Recompute A if modif

      if (n_modif > 0) {

         for (i = 0; i < _b_size; i++) {
            for (j = 0; j < _b_size; j++) {
               pUSMatr[i * _b_size + j] = pUMatr[i * _b_size + j] * pSv[i];
               pVSMatr[i * _b_size + j] = pVMatr[i * _b_size + j] * pSv[i];
            }
         }

         CBlock_BxB_traits < _Flt >::MMt_BxB (_b_size, pUSMatr, pVSMatr, _A_blk);

      }
// Compute Ainv

      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            pUSMatr[i * _b_size + j] = pVMatr[i * _b_size + j] * pSvInv[i];
            pVSMatr[i * _b_size + j] = pUMatr[i * _b_size + j] * pSvInv[i];
         }
      }

      CBlock_BxB_traits < _Flt >::MMt_BxB (_b_size, pUSMatr, pVSMatr, _A_inv_blk);

      return n_modif;

   }

//
// Compute fct and its inverse via two-sided QR
//========================================================================================
   template < typename _Flt > int CBlock_BxB_traits <
      _Flt >::FctInv_BiQrd_BxB (double _piv_min, int _b_size, const _Flt * _A_blk,
                                _Flt * _Lt_blk, _Flt * _U_blk, _Flt * _Lt_inv_blk,
                                _Flt * _U_inv_blk, _Flt * _Work, double &_dia_min,
                                double &_dia_max)
   {

      if (_b_size == 1) {
         int n_modif = 0;
         n_modif =
            CBlock_BxB_traits < _Flt >::FctInv_BxB (_piv_min, _b_size, _A_blk, _Lt_blk,
                                                    _U_blk, _Lt_inv_blk, _U_inv_blk,
                                                    _dia_min, _dia_max);
         return n_modif;
      }

      static const _Flt fone = (_Flt) 1.0e0;

      int b_2 = _b_size * _b_size;

      _Flt *pBiDiag = _Work;
      _Flt *pDiag = pBiDiag + _b_size * 2;
      _Flt *pDiagInv = pDiag + _b_size;
      _Flt *pUMatr = pDiagInv + _b_size;
      _Flt *pVMatr = pUMatr + b_2;
      _Flt *pBMatr = pVMatr + b_2;
      _Flt *pBMatrInv = pBMatr + b_2;
      _Flt *pWork = pBMatrInv + b_2;

      CBlock_BxB_traits < _Flt >::BiDiagonalize_BxB (_b_size, _A_blk, pBiDiag, pUMatr,
                                                     pVMatr, pWork);

      int n_modif = 0;

      int i, j;

      _Flt aux = pBiDiag[0];

      if (aux < 0.0e0)
         aux = -aux;

      _dia_min = aux;
      _dia_max = aux;

      for (i = 0; i < _b_size; i++) {
         aux = pBiDiag[i];
         if (aux < 0.0e0)
            aux = -aux;
         if (aux < _dia_min)
            _dia_min = aux;
         if (aux > _dia_max)
            _dia_max = aux;
         if (aux < _piv_min) {
            n_modif++;
            aux = (_Flt) _piv_min;
         }
         pDiag[i] = sqrt (aux);
         if (pBiDiag[i] > 0.0e0) {
            pBiDiag[i] = pDiag[i];
         } else {
            pBiDiag[i] = -pDiag[i];
         }
         pDiagInv[i] = fone / pDiag[i];
         pBiDiag[_b_size + i] *= pDiagInv[i];
      }

      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            _Lt_blk[j * _b_size + i] = pUMatr[i * _b_size + j] * pDiag[i];
            _Lt_inv_blk[i * _b_size + j] = pUMatr[i * _b_size + j] * pDiagInv[i];
         }
      }

      CVector < _Flt >::SetByZeroes (b_2, pBMatr);

      for (i = 0; i < _b_size; i++) {
         pBMatr[i * _b_size + i] = pBiDiag[i];
         if (i < _b_size - 1) {
            pBMatr[(i + 1) * _b_size + i] = pBiDiag[_b_size + i];
         }
      }

      CBlock_BxB_traits < _Flt >::MMt_BxB (_b_size, pBMatr, pVMatr, _U_blk);

      CBlock_BxB_traits < _Flt >::InvU_BxB (_b_size, pBMatr, pBMatrInv);

      CBlock_BxB_traits < _Flt >::MM_BxB (_b_size, pVMatr, pBMatrInv, _U_inv_blk);

      return n_modif;

   }

//
// Compute fct and its inverse via SVD
//========================================================================================
   template < typename _Flt > int CBlock_BxB_traits <
      _Flt >::FctInv_Svd_BxB (double _piv_min, int _b_size, const _Flt * _A_blk,
                              _Flt * _Lt_blk, _Flt * _U_blk, _Flt * _Lt_inv_blk,
                              _Flt * _U_inv_blk, _Flt * _Work, double *_dWork,
                              double &_dia_min, double &_dia_max)
   {

      static const _Flt fone = (_Flt) 1.0e0;

      int b_2 = _b_size * _b_size;

      _Flt *pSv = _Work;
      _Flt *pSvInv = pSv + _b_size;
      _Flt *pUMatr = pSvInv + _b_size;
      _Flt *pVMatr = pUMatr + b_2;

      CBlock_BxB_traits < _Flt >::ComputeSvd_BxB (_b_size, _A_blk, pSv, pUMatr, pVMatr,
                                                  _dWork);

      int n_modif = 0;

      int i, j;

      _dia_min = pSv[0];
      _dia_max = pSv[0];

      for (i = 0; i < _b_size; i++) {
         if (pSv[i] < _dia_min)
            _dia_min = (double) pSv[i];
         if (pSv[i] > _dia_max)
            _dia_max = (double) pSv[i];
         if (pSv[i] < _piv_min) {
            n_modif++;
            pSv[i] = (_Flt) _piv_min;
         }
         pSv[i] = sqrt (pSv[i]);
         pSvInv[i] = fone / pSv[i];
      }

      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            _Lt_blk[j * _b_size + i] = pUMatr[i * _b_size + j] * pSv[i];
            _U_blk[j * _b_size + i] = pVMatr[i * _b_size + j] * pSv[i];
            _Lt_inv_blk[i * _b_size + j] = pUMatr[i * _b_size + j] * pSvInv[i];
            _U_inv_blk[i * _b_size + j] = pVMatr[i * _b_size + j] * pSvInv[i];
         }
      }

      return n_modif;

   }

//
// Inverse upper triangular
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::InvU_BxB (int _b_size,
                                                                         const _Flt *
                                                                         _U_blk,
                                                                         _Flt *
                                                                         _U_inv_blk)
   {

      static const _Flt fzero = (_Flt) 0.0e0;
      static const _Flt fone = (_Flt) 1.0e0;

      int i, j, k;
      for (j = 0; j < _b_size; j++)
         for (i = 0; i < _b_size; i++)
            _U_inv_blk[i + j * _b_size] = fzero;
      for (i = 0; i < _b_size; i++)
         _U_inv_blk[i + i * _b_size] = fone;
      for (j = 0; j < _b_size; j++) {
         _U_inv_blk[j + j * _b_size] = fone / _U_blk[j + j * _b_size];
         for (i = j - 1; i >= 0; i--) {
            for (k = 0; k < j - i; k++)
               _U_inv_blk[i + j * _b_size] -=
                  _U_blk[i + (j - k) * _b_size] * _U_inv_blk[j - k + j * _b_size];
            _U_inv_blk[i + j * _b_size] *= _U_inv_blk[i + i * _b_size];
         }
      }
   }

//
// Inverse upper triangular
//========================================================================================
   template < typename _Flt > int CBlock_BxB_traits < _Flt >::InvU_BxB (double _piv_min,
                                                                        int _b_size,
                                                                        _Flt * _U_blk,
                                                                        _Flt * _U_inv_blk,
                                                                        double &_dia_min,
                                                                        double &_dia_max)
   {

      double aux = (double) _U_blk[0];

      if (aux < 0.0e0)
         aux = -aux;

      int n_modif = 0;
      _dia_min = aux;
      _dia_max = aux;

      int i;
      double aux1;

      for (i = 0; i < _b_size; i++) {
         aux = (double) _U_blk[i * _b_size + i];
         aux1 = aux;
         if (aux1 < 0.0e0)
            aux1 = -aux1;
         if (aux1 < _dia_min)
            _dia_min = aux1;
         if (aux1 > _dia_max)
            _dia_max = aux1;
         if (aux1 < _piv_min) {
            n_modif++;
            if (aux < 0.0e0) {
               _U_blk[i * _b_size + i] = -(_Flt) _piv_min;
            } else {
               _U_blk[i * _b_size + i] = (_Flt) _piv_min;
            }
         }
      }

      CBlock_BxB_traits < _Flt >::InvU_BxB (_b_size, _U_blk, _U_inv_blk);

      return n_modif;

   }

//
// Modify diagonal when filtering
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::ModifDia_BxB (int _b_size,
                                                                             double
                                                                             _theta,
                                                                             const _Flt *
                                                                             _L_JI_offd,
                                                                             const _Flt *
                                                                             _U_IJ_offd,
                                                                             _Flt *
                                                                             _dia_I,
                                                                             _Flt *
                                                                             _dia_J)
   {

      static const _Flt fzero = (_Flt) 0.0e0;

      _Flt ftheta = (_Flt) _theta;

      _Flt fnorm;
      _Flt faux;

      int i, j, k;

      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            k = j * _b_size + i;
            fnorm = 0.0e0;
            faux = _U_IJ_offd[k];
            if (faux < fzero)
               faux = -faux;
            if (faux > fnorm)
               fnorm = faux;
            faux = _L_JI_offd[k];
            if (faux < fzero)
               faux = -faux;
            if (faux > fnorm)
               fnorm = faux;
            fnorm *= ftheta;
            _dia_I[i * _b_size + i] += fnorm;
            _dia_J[j * _b_size + j] += fnorm;
         }
      }

   }

//
// Compute norm of BxB block
//========================================================================================
   template < typename _Flt > double CBlock_BxB_traits <
      _Flt >::FNormValue_BxB (int _b_size, const _Flt * _A)
   {

      static const _Flt fzero = (_Flt) 0.0e0;

      int b_2 = _b_size * _b_size;
      _Flt fnorm = fzero;

      _Flt faux;

      int i;

      for (i = 0; i < b_2; i++) {
         faux = _A[i];
         if (faux < fzero)
            faux = -faux;
         if (faux > fnorm)
            fnorm = faux;
      }

      double dnorm = (double) fnorm;

      return dnorm;

   }

//
// Compute norm of BxB block
//========================================================================================
   template < typename _Flt > double CBlock_BxB_traits <
      _Flt >::FNormValue_BxB (int _b_size, const _Flt * _L_offd, const _Flt * _U_offd)
   {

      static const _Flt fzero = (_Flt) 0.0e0;

      int b_2 = _b_size * _b_size;
      _Flt fnorm = fzero;

      _Flt faux;

      int i;

      for (i = 0; i < b_2; i++) {
         faux = _L_offd[i];
         if (faux < fzero)
            faux = -faux;
         if (faux > fnorm)
            fnorm = faux;
         faux = _U_offd[i];
         if (faux < fzero)
            faux = -faux;
         if (faux > fnorm)
            fnorm = faux;
      }

      double dnorm = (double) fnorm;

      return dnorm;

   }

//
// Add diagonal modification elements
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits <
      _Flt >::AddModifDia_BxB (int _b_size, _Flt * _dia_modif, _Flt * _dia_U)
   {

      static const _Flt fzero = (_Flt) 0.0e0;

      int i;

      for (i = 0; i < _b_size; i++) {
//         if (_dia_U[i*_b_size+i] >= fzero) {
         _dia_U[i * _b_size + i] += _dia_modif[i * _b_size + i];
//         } else {
//            _dia_U[i*_b_size+i] -= _dia_modif[i*_b_size+i];
//         }
      }

   }

//
// Compute Svd
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits <
      _Flt >::ComputeSvd_BxB (int _b_size, const _Flt * _a_matr, _Flt * _Sv_arr,
                              _Flt * _U_matr, _Flt * _V_matr, double *_dWork)
   {

#if defined(USE_LAPACK)

      int b_2 = _b_size * _b_size;

      double *psv_ini = _dWork;
      double *pamatr = psv_ini + _b_size;
      double *psv = pamatr + b_2;
      double *pu = psv + _b_size;
      double *pvt = pu + b_2;
      double *pwork = pvt + b_2;
      int lwork = 10 * _b_size;

      int i;

      for (i = 0; i < b_2; i++)
         pamatr[i] = (double) _a_matr[i];

      int info;

      dgesvd_ ("A", "A", &_b_size, &_b_size, pamatr, &_b_size, psv, pu, &_b_size, pvt,
               &_b_size, pwork, &lwork, &info);

      for (i = 0; i < _b_size; i++)
         psv_ini[i] = psv[i];
      for (i = 0; i < _b_size; i++)
         _Sv_arr[i] = (_Flt) psv[i];
      for (i = 0; i < b_2; i++)
         _U_matr[i] = (_Flt) pu[i];

      int j;

      for (i = 0; i < _b_size; i++) {
         for (j = 0; j < _b_size; j++) {
            _V_matr[i * _b_size + j] = (_Flt) pvt[j * _b_size + i];
         }
      }

#else
      cout << " CBlock_BxB_traits<>::ComputeSvd_BxB: error: lapack Svd is not found !!! "
         << endl;
      throw " CBlock_BxB_traits<>::ComputeSvd_BxB: error: lapack Svd is not found !!! ";
#endif

   }

//
// Compute Qrd in column format
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::ComputeQrd_BxB (int _m,
                                                                               int
                                                                               _b_size,
                                                                               const _Flt
                                                                               * _A_matr,
                                                                               _Flt *
                                                                               _Q_matr,
                                                                               _Flt *
                                                                               _R_matr,
                                                                               _Flt *
                                                                               _Work)
   {

      _Flt *pWork = _Work;

      _Flt *pqblk = pWork;
      pWork += _m * _b_size;

      _Flt *ptau = pWork;
      pWork += _b_size;

      CVector < _Flt >::CopyVector (_m * _b_size, _A_matr, pqblk);

      CVector < _Flt >::QrdBlock (_b_size, _m, pqblk, _m, ptau);

      CVector < _Flt >::SetByZeroes (_b_size * _b_size, _R_matr);

      int i;

      for (i = 0; i < _b_size; i++)
         CVector < _Flt >::CopyVector (i + 1, pqblk + i * _m, _R_matr + i * _b_size);

      CVector < _Flt >::SetByZeroes (_m * _b_size, _Q_matr);

      for (i = 0; i < _b_size; i++)
         CVector < _Flt >::SetByOnes (1, _Q_matr + i * _m + i);

      CVector < _Flt >::MvmQ_Housholder (_b_size, _m, _b_size, pqblk, _m, ptau, _Q_matr,
                                         _m);

   }

///
/// @brief Two sided bidiagonalization
//========================================================================================
   template < typename _Flt > void CBlock_BxB_traits < _Flt >::BiDiagonalize_BxB (int _n,
                                                                                  const
                                                                                  _Flt *
                                                                                  _amatr,
                                                                                  _Flt *
                                                                                  _bidiag,
                                                                                  _Flt *
                                                                                  _uarr,
                                                                                  _Flt *
                                                                                  _vharr,
                                                                                  _Flt *
                                                                                  _Work)
   {

      if (_n == 1) {
         _Flt fone = (_Flt) 1.0e0;
         _bidiag[0] = _amatr[0];
         _uarr[0] = fone;
         _vharr[0] = fone;
         return;
      }

      int n_2 = _n * _n;

      _Flt *pwork = _Work;

      _Flt *paloc = pwork;
      pwork += n_2;
      _Flt *ptau_u = pwork;
      pwork += _n;
      _Flt *ptau_v = pwork;
      pwork += _n;
      _Flt *pu_rot = pwork;
      pwork += n_2;
      _Flt *pv_rot = pwork;
      pwork += n_2;

      _Flt fzero;
      _Flt fone;

      CVector < _Flt >::SetByZeroes (1, &fzero);
      CVector < _Flt >::SetByOnes (1, &fone);

      CVector < _Flt >::CopyVector (n_2, _amatr, paloc);

      CVector < _Flt >::SetByZeroes (n_2, pu_rot);
      CVector < _Flt >::SetByZeroes (n_2, pv_rot);

      CVector < _Flt >::SetByZeroes (_n * 2, _bidiag);

      int i, j, k, kloc;
      _Flt scprod, scprod_minus, aux, alpha;
      _Flt *pvect;
      _Flt *pa;
      _Flt *ph;

      for (i = 0; i < _n; i++) {

// Compute current left transformation

         alpha = paloc[i * _n + i];

         CVector < _Flt >::CopyVector (_n - i - 1, paloc + i * _n + i + 1,
                                       pu_rot + i * _n + i + 1);

         pvect = NULL;

         if (i < _n - 1)
            pvect = pu_rot + i * _n + i + 1;

         CVector < _Flt >::Housholder (_n - i, alpha, pvect, ptau_u[i]);

// Apply transformation to the rest of the matrix

         _bidiag[i] = alpha;

         for (j = i; j < _n; j++) {

            pa = paloc + j * _n + i + 1;
            ph = pu_rot + i * _n + i + 1;

            scprod = CVector < _Flt >::ScProd (_n - i - 1, pa, ph);

            scprod += paloc[j * _n + i];

            scprod *= ptau_u[i];

            paloc[j * _n + i] -= scprod;

            scprod_minus = -scprod;

            CVector < _Flt >::UpdateVector (_n - i - 1, &scprod_minus, ph, pa);

         }

// Compute current right transformation

         if (i < _n - 1) {

            alpha = paloc[(i + 1) * _n + i];

            for (j = i + 2; j < _n; j++) {
               pv_rot[(i + 1) * _n + j] = paloc[j * _n + i];
            }

            pvect = NULL;

            if (i < _n - 2)
               pvect = pv_rot + (i + 1) * _n + i + 2;

            CVector < _Flt >::Housholder (_n - i - 1, alpha, pvect, ptau_v[i + 1]);

            _bidiag[_n + i] = alpha;

// Apply it to the rest of the matrix

            for (j = i; j < _n; j++) {
               aux = paloc[(i + 1) * _n + j];
               for (k = i + 2; k < _n; k++) {
                  kloc = k - i - 2;
                  aux += paloc[k * _n + j] * pvect[kloc];
               }
               aux *= (-ptau_v[i + 1]);
               paloc[(i + 1) * _n + j] += aux;
               for (k = i + 2; k < _n; k++) {
                  kloc = k - i - 2;
                  paloc[k * _n + j] += aux * pvect[kloc];
               }
            }

         }

      }

// Compute the rotation matrices explicitely

      for (j = 0; j < _n; j++) {
         for (i = 0; i < _n; i++) {
            _uarr[j * _n + i] = fzero;
         }
      }
      for (i = 0; i < _n; i++) {
         _uarr[i * _n + i] = fone;
         for (j = i; j >= 0; j--) {
            aux = _uarr[i * _n + j];
            for (k = j + 1; k < _n; k++) {
               aux += _uarr[i * _n + k] * pu_rot[j * _n + k];
            };
            aux *= -ptau_u[j];
            _uarr[i * _n + j] += aux;
            for (k = j + 1; k < _n; k++) {
               _uarr[i * _n + k] += aux * pu_rot[j * _n + k];
            }
         }
      }

      for (j = 0; j < _n; j++) {
         for (i = 0; i < _n; i++) {
            _vharr[j * _n + i] = fzero;
         }
      }
      for (i = 0; i < _n; i++) {
         _vharr[i * _n + i] = fone;
         for (j = i; j >= 1; j--) {
            aux = _vharr[i * _n + j];
            for (k = j + 1; k < _n; k++) {
               aux += _vharr[i * _n + k] * pv_rot[j * _n + k];
            };
            aux *= -ptau_v[j];
            _vharr[i * _n + j] += aux;
            for (k = j + 1; k < _n; k++) {
               _vharr[i * _n + k] += aux * pv_rot[j * _n + k];
            }
         }
      }

   }

//
// Find ordering as splitting diagonal blocks with very small pivot/singular values
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::OrderSplitDiaPairs (SParams * _params, int _b_size, int _n,
                                  const vector < _Int > &_ia_au,
                                  const vector < _Int > &_ja_au,
                                  const vector < char >&_ja_char_au,
                                  const vector < _Flt > &_a_au, int *_order, int &_n1)
   {

// Get params

      int fcttype = _params->fcttype;
      int fcttype_dia = _params->fcttype_dia;
      int ordtype = _params->ordtype;
      double pivmin = _params->pivmin;
      double dia_split = _params->dia_split;

// Fast return

      if (dia_split < 0.0e0) {

         vector < int >order;

         CFct_impl < _Int, _Flt >::ComputeOptimalOrder (ordtype, _n, _ia_au, _ja_au,
                                                        order);

         int *porder = &order[0];

         int i;
         for (i = 0; i < _n; i++)
            _order[i] = porder[i];

         _n1 = _n;

         return;

      }
// Open matrices data structures

      int b_2 = _b_size * _b_size;
      int b_2_2 = 2 * b_2;

      const _Int *p_ia = NULL;
      const _Int *p_ja = NULL;
      const char *p_ja_char = NULL;
      const _Flt *p_a = NULL;

      if (_ia_au.size () > 0)
         p_ia = &_ia_au[0];
      if (_ja_au.size () > 0)
         p_ja = &_ja_au[0];
      if (_ja_char_au.size () > 0)
         p_ja_char = &_ja_char_au[0];
      if (_a_au.size () > 0)
         p_a = &_a_au[0];

// Scan diagonal blocks

      CVectorData < _Flt > lu_diag (b_2_2);

      _Flt *pl_diag = lu_diag.Ptr ();
      _Flt *pu_diag = pl_diag + b_2;

      CVectorData < _Flt > lu_offd (b_2_2);

      _Flt *pl_offd = lu_offd.Ptr ();
      _Flt *pu_offd = pl_offd + b_2;

      CVectorData < _Flt > Work (10 * b_2);
      CVectorData < double >dWork (20 * b_2);

      _Flt *pWork = Work.Ptr ();
      double *pdWork = dWork.Ptr ();

      CVectorData < int >imask (_n);
      int *pimask = imask.Ptr ();

      int i;

      for (i = 0; i < _n; i++)
         pimask[i] = 1;

      int j;
      double eigmin_temp, eigmax_temp;

      for (i = 0; i < _n; i++) {
         j = (int) p_ia[i];
         if (fcttype_dia == 0) {

            CBlock_BxB_traits < _Flt >::FctInv_BxB (pivmin, _b_size,
                                                    p_a + j * b_2_2 + b_2, pl_offd,
                                                    pu_offd, pl_diag, pu_diag,
                                                    eigmin_temp, eigmax_temp);

         } else {

            CBlock_BxB_traits < _Flt >::FctInv_Svd_BxB (pivmin, _b_size,
                                                        p_a + j * b_2_2 + b_2, pl_offd,
                                                        pu_offd, pl_diag, pu_diag, pWork,
                                                        pdWork, eigmin_temp, eigmax_temp);

         }

         if (eigmin_temp < dia_split)
            pimask[i] = 0;

      }

// Compute initial order

      _n1 = 0;

      for (i = 0; i < _n; i++) {
         if (pimask[i] > 0)
            _n1++;
      }

      if (_n1 == 0 || _n1 == _n) {

// Again fast return

         vector < int >order;

         CFct_impl < _Int, _Flt >::ComputeOptimalOrder (ordtype, _n, _ia_au, _ja_au,
                                                        order);

         int *porder = &order[0];

         int i;
         for (i = 0; i < _n; i++)
            _order[i] = porder[i];

         return;

      }

      int i1 = 0;
      int i2 = _n1;

      vector < int >order (_n + 1);
      int *porder = &order[0];

      for (i = 0; i < _n; i++) {
         if (pimask[i] > 0) {
            porder[i] = i1;
            i1++;
         } else {
            porder[i] = i2;
            i2++;
         }
      }

// Symmetrize sparsity

      vector < _Int > ia_symm;
      vector < _Int > ja_symm;
      vector < char >jachar_symm;

      CFct_impl < _Int, _Flt >::SymmetrizeSparsity (_n, _ia_au, _ja_au, _ja_char_au,
                                                    ia_symm, ja_symm, jachar_symm);

// Order matrix sparsity

      vector < _Int > ia_ord;
      vector < _Int > ja_ord;

      CFct_impl < _Int, _Flt >::ReorderMatrixSp (_n, order, ia_symm, ja_symm, ia_ord,
                                                 ja_ord);

// Compute optimal suborder for the first submatrix

      vector < _Int > ia_sub;
      vector < _Int > ja_sub;

      CFct_impl < _Int, _Flt >::GetSubmatrixSp (_n, ia_ord, ja_ord, _n1, ia_sub, ja_sub);

      vector < int >order1;

      CFct_impl < _Int, _Flt >::ComputeOptimalOrder (ordtype, _n1, ia_sub, ja_sub,
                                                     order1);

      int *porder1 = &order1[0];

// Extend order

      vector < int >order1_ext (_n + 1);
      int *porder1_ext = &order1_ext[0];

      for (i = 0; i < _n1; i++)
         porder1_ext[i] = porder1[i];
      for (i = _n1; i < _n; i++)
         porder1_ext[i] = i;

      vector < int >iorder (_n + 1);
      int *piorder = &iorder[0];

      for (i = 0; i < _n; i++)
         piorder[porder[i]] = i;

      int iold;

      for (i = 0; i < _n; i++) {
         iold = piorder[i];
         porder[iold] = porder1_ext[i];
      }

// Fast return

      const int n_small = 10;

      if (_n - _n1 < n_small) {

         for (i = 0; i < _n; i++)
            _order[i] = porder[i];

         return;
      }
// Reorder matrix again

      vector < char >jachar_ord;

      CFct_impl < _Int, _Flt >::ReorderMatrixSp (_n, order, ia_symm, ja_symm, jachar_symm,
                                                 ia_ord, ja_ord, jachar_ord);

// Split sparsity into L and U parts

      vector < _Int > ia_AL;
      vector < _Int > ja_AL;
      vector < char >jachar_AL;
      vector < _Int > ia_AU;
      vector < _Int > ja_AU;
      vector < char >jachar_AU;

      CFct_impl < _Int, _Flt >::SplitLUSp (_n, ia_ord, ja_ord, jachar_ord, ia_AL, ja_AL,
                                           jachar_AL, ia_AU, ja_AU, jachar_AU);

// Get main and last submatrices with char values

      vector < char >jachar_sub;

      CFct_impl < _Int, _Flt >::GetSubmatrixSp (_n, ia_AU, ja_AU, jachar_AU, _n1, ia_sub,
                                                ja_sub, jachar_sub);

      vector < _Int > ia_sub_last;
      vector < _Int > ja_sub_last;

      CFct_impl < _Int, _Flt >::GetLastSubmatrixSp (_n, ia_AU, ja_AU, _n1, ia_sub_last,
                                                    ja_sub_last);

// Perform symbolic fct for the first submatrix 

      vector < _Int > ia_U;
      vector < _Int > ja_U;
      vector < char >jachar_U;

      CFct_impl < _Int, _Flt >::Ilu2BlockIlu2DegreeSp (NULL, fcttype, fcttype, _n, _n1,
                                                       ia_sub, ja_sub, jachar_sub, ia_U,
                                                       ja_U, jachar_U);

// Get last submatrix U

      vector < _Int > ia_U_last;
      vector < _Int > ja_U_last;

      CFct_impl < _Int, _Flt >::GetLastSubmatrixSp (_n, ia_U, ja_U, _n1, ia_U_last,
                                                    ja_U_last);

// Add sparsities

      int n2 = _n - _n1;

      int nzja_sub_last = (int) ia_sub_last[n2];

      CMatrix < _Int, _Flt > a_sub_last;

      a_sub_last.ResizeAndSetAllSp (n2, 0, nzja_sub_last, 0);

      _Int *plist_sub_last = a_sub_last.GetListArr ();
      vector < _Int > *pia_sub_set = a_sub_last.GetIa ();
      vector < _Int > *pja_sub_set = a_sub_last.GetJa ();

      for (i = 0; i < n2; i++)
         plist_sub_last[i] = i;

      pia_sub_set->swap (ia_sub_last);
      pja_sub_set->swap (ja_sub_last);

      int nzja_U_last = (int) ia_U_last[n2];

      CMatrix < _Int, _Flt > a_U_last;

      a_U_last.ResizeAndSetAllSp (n2, 0, nzja_U_last, 0);

      _Int *plist_U_last = a_U_last.GetListArr ();
      vector < _Int > *pia_U_set = a_sub_last.GetIa ();
      vector < _Int > *pja_U_set = a_sub_last.GetJa ();

      for (i = 0; i < n2; i++)
         plist_U_last[i] = i;

      pia_U_set->swap (ia_U_last);
      pja_U_set->swap (ja_U_last);

      CMatrix < _Int, _Flt > a_sum;

      a_sum.AddBlocksSp (a_sub_last, a_U_last);

      vector < _Int > *pia_sum_set = a_sum.GetIa ();
      vector < _Int > *pja_sum_set = a_sum.GetJa ();

// Final optimal order

      vector < int >order2;

      CFct_impl < _Int, _Flt >::ComputeOptimalOrder (ordtype, n2, *pia_sum_set,
                                                     *pja_sum_set, order2);

      int *porder2 = &order2[0];

      for (i = 0; i < _n1; i++)
         porder1_ext[i] = i;
      for (i = 0; i < n2; i++)
         porder1_ext[_n1 + i] = _n1 + porder2[i];

      for (i = 0; i < _n; i++)
         piorder[porder[i]] = i;

      for (i = 0; i < _n; i++) {
         iold = piorder[i];
         _order[iold] = porder1_ext[i];
      }

   }

//
// Split matrix into submatrices by rows (pairs or not)
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::SubmatricesByRows (bool _b_is_char, bool _b_is_pair, int _b_size, int _n,
                                 int _n1, const vector < _Int > &_ia,
                                 const vector < _Int > &_ja,
                                 const vector < char >&_ja_char,
                                 const vector < _Flt > &_a, vector < _Int > &_ia_ini,
                                 vector < _Int > &_ja_ini, vector < char >&_ja_char_ini,
                                 vector < _Flt > &_a_ini, vector < _Int > &_ia_last,
                                 vector < _Int > &_ja_last, vector < char >&_ja_char_last,
                                 vector < _Flt > &_a_last)
   {

// Open matrices data structures

      int b_2 = _b_size * _b_size;
      int b_2_2 = 2 * b_2;

      int b_2_curr = b_2;
      if (_b_is_pair)
         b_2_curr = b_2_2;

      const _Int *p_ia = _ia.data ();
      const _Int *p_ja = _ja.data ();
      const char *p_ja_char = _ja_char.data ();
      const _Flt *p_a = _a.data ();

// Count size of both data and fill ia

      int n_last = _n - _n1;

      int nzja_ini = (int) p_ia[_n1];
      int nzja_last = (int) (p_ia[_n] - nzja_ini);

      _ia_ini.resize (_n1 + 1);
      _ja_ini.resize (nzja_ini + 1);
      if (_b_is_char)
         _ja_char_ini.resize (nzja_ini + 1);
      _a_ini.resize (nzja_ini * b_2_curr + 1);

      _ia_last.resize (n_last + 1);
      _ja_last.resize (nzja_last + 1);
      if (_b_is_char)
         _ja_char_last.resize (nzja_last + 1);
      _a_last.resize (nzja_last * b_2_curr + 1);

      _Int *p_ia_ini = _ia_ini.data ();
      _Int *p_ja_ini = _ja_ini.data ();
      char *p_ja_char_ini = _ja_char_ini.data ();
      _Flt *p_a_ini = _a_ini.data ();

      _Int *p_ia_last = _ia_last.data ();
      _Int *p_ja_last = _ja_last.data ();
      char *p_ja_char_last = _ja_char_last.data ();
      _Flt *p_a_last = _a_last.data ();

      int i;

      for (i = 0; i <= _n1; i++)
         p_ia_ini[i] = p_ia[i];
      for (i = 0; i < nzja_ini; i++)
         p_ja_ini[i] = p_ja[i];
      if (_b_is_char) {
         for (i = 0; i < nzja_ini; i++)
            p_ja_char_ini[i] = p_ja_char[i];
      }
      CVector < _Flt >::CopyVector (b_2_curr * nzja_ini, p_a, p_a_ini);

      for (i = 0; i <= n_last; i++)
         p_ia_last[i] = p_ia[i + _n1] - nzja_ini;
      for (i = 0; i < nzja_last; i++)
         p_ja_last[i] = p_ja[i + nzja_ini];
      if (_b_is_char) {
         for (i = 0; i < nzja_last; i++)
            p_ja_char_last[i] = p_ja_char[i + nzja_ini];
      }
      CVector < _Flt >::CopyVector (b_2_curr * nzja_last, p_a + nzja_ini * b_2_curr,
                                    p_a_last);

   }

//
// Add sparse matrices with elements (pairs or not)
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::AddMatrices (bool _b_is_char, bool _b_is_pair, int _b_size, int _n,
                           const vector < _Int > &_ia_1, const vector < _Int > &_ja_1,
                           const vector < char >&_ja_char_1, const vector < _Flt > &_a_1,
                           const vector < _Int > &_ia_2, const vector < _Int > &_ja_2,
                           const vector < char >&_ja_char_2, const vector < _Flt > &_a_2,
                           vector < _Int > &_ia_sum, vector < _Int > &_ja_sum,
                           vector < char >&_ja_char_sum, vector < _Flt > &_a_sum)
   {

// Open matrices data structures

      int b_2 = _b_size * _b_size;
      int b_2_2 = 2 * b_2;

      int b_2_add = b_2;
      if (_b_is_pair)
         b_2_add = b_2_2;

      const _Int *p_ia_1 = NULL;
      const _Int *p_ja_1 = NULL;
      const char *p_ja_char_1 = NULL;
      const _Flt *p_a_1 = NULL;

      if (_ia_1.size () > 0)
         p_ia_1 = &_ia_1[0];
      if (_ja_1.size () > 0)
         p_ja_1 = &_ja_1[0];
      if (_b_is_char && _ja_char_1.size () > 0)
         p_ja_char_1 = &_ja_char_1[0];
      if (_a_1.size () > 0)
         p_a_1 = &_a_1[0];

      const _Int *p_ia_2 = NULL;
      const _Int *p_ja_2 = NULL;
      const char *p_ja_char_2 = NULL;
      const _Flt *p_a_2 = NULL;

      if (_ia_2.size () > 0)
         p_ia_2 = &_ia_2[0];
      if (_ja_2.size () > 0)
         p_ja_2 = &_ja_2[0];
      if (_b_is_char && _ja_char_2.size () > 0)
         p_ja_char_2 = &_ja_char_2[0];
      if (_a_2.size () > 0)
         p_a_2 = &_a_2[0];

// Count size of data and fill ia

      _ia_sum.resize (_n + 1);
      _Int *p_ia_sum = &_ia_sum[0];

      int nzja_sum = 0;
      p_ia_sum[0] = 0;

      int i, jp1, jp2, jend1, jend2, jj1, jj2;

      for (i = 0; i < _n; i++) {
         jend1 = (int) p_ia_1[i + 1] - 1;
         jend2 = (int) p_ia_2[i + 1] - 1;
         jp1 = (int) p_ia_1[i];
         jp2 = (int) p_ia_2[i];
         while (jp1 <= jend1 || jp2 <= jend2) {
            if (jp1 <= jend1 && jp2 <= jend2) {
               jj1 = (int) p_ja_1[jp1];
               jj2 = (int) p_ja_2[jp2];
               if (jj1 == jj2) {
                  nzja_sum++;
                  jp1++;
                  jp2++;
               } else if (jj1 < jj2) {
                  nzja_sum++;
                  jp1++;
               } else if (jj1 > jj2) {
                  nzja_sum++;
                  jp2++;
               }
            } else if (jp1 <= jend1) {
               nzja_sum++;
               jp1++;
            } else if (jp2 <= jend2) {
               nzja_sum++;
               jp2++;
            }
         }
         p_ia_sum[i + 1] = nzja_sum;
      }

// Allocate and init by zeroes

      _ja_sum.resize (nzja_sum + 1);
      if (_b_is_char) {
         _ja_char_sum.resize (nzja_sum + 1);
      }
      _a_sum.resize (nzja_sum * b_2_add + 1);

      _Int *p_ja_sum = &_ja_sum[0];
      char *p_ja_char_sum = NULL;
      if (_b_is_char) {
         p_ja_char_sum = &_ja_char_sum[0];
      }
      _Flt *p_a_sum = &_a_sum[0];

      nzja_sum = 0;
      p_ia_sum[0] = 0;

      int ibs1, ibs2, ibs;

      for (i = 0; i < _n; i++) {
         jend1 = (int) p_ia_1[i + 1] - 1;
         jend2 = (int) p_ia_2[i + 1] - 1;
         jp1 = (int) p_ia_1[i];
         jp2 = (int) p_ia_2[i];
         while (jp1 <= jend1 || jp2 <= jend2) {
            if (jp1 <= jend1 && jp2 <= jend2) {
               jj1 = (int) p_ja_1[jp1];
               jj2 = (int) p_ja_2[jp2];
               if (jj1 == jj2) {
                  p_ja_sum[nzja_sum] = jj1;
                  if (_b_is_char) {
                     p_ja_char_sum[nzja_sum] = p_ja_char_1[jp1];
                     if (p_ja_char_2[jp2] < p_ja_char_1[jp1]) {
                        p_ja_char_sum[nzja_sum] = p_ja_char_2[jp2];
                     }
                  }
                  ibs1 = jp1 * b_2_add;
                  ibs2 = jp2 * b_2_add;
                  ibs = nzja_sum * b_2_add;
                  CVector < _Flt >::AddVector (b_2_add, p_a_1 + ibs1, p_a_2 + ibs2,
                                               p_a_sum + ibs);
                  nzja_sum++;
                  jp1++;
                  jp2++;
               } else if (jj1 < jj2) {
                  p_ja_sum[nzja_sum] = jj1;
                  if (_b_is_char) {
                     p_ja_char_sum[nzja_sum] = p_ja_char_1[jp1];
                  }
                  ibs1 = jp1 * b_2_add;
                  ibs = nzja_sum * b_2_add;
                  CVector < _Flt >::CopyVector (b_2_add, p_a_1 + ibs1, p_a_sum + ibs);
                  nzja_sum++;
                  jp1++;
               } else if (jj1 > jj2) {
                  p_ja_sum[nzja_sum] = jj2;
                  if (_b_is_char) {
                     p_ja_char_sum[nzja_sum] = p_ja_char_2[jp2];
                  }
                  ibs2 = jp2 * b_2_add;
                  ibs = nzja_sum * b_2_add;
                  CVector < _Flt >::CopyVector (b_2_add, p_a_2 + ibs2, p_a_sum + ibs);
                  nzja_sum++;
                  jp2++;
               }
            } else if (jp1 <= jend1) {
               p_ja_sum[nzja_sum] = p_ja_1[jp1];
               if (_b_is_char) {
                  p_ja_char_sum[nzja_sum] = p_ja_char_1[jp1];
               }
               ibs1 = jp1 * b_2_add;
               ibs = nzja_sum * b_2_add;
               CVector < _Flt >::CopyVector (b_2_add, p_a_1 + ibs1, p_a_sum + ibs);
               nzja_sum++;
               jp1++;
            } else if (jp2 <= jend2) {
               p_ja_sum[nzja_sum] = p_ja_2[jp2];
               if (_b_is_char) {
                  p_ja_char_sum[nzja_sum] = p_ja_char_2[jp2];
               }
               ibs2 = jp2 * b_2_add;
               ibs = nzja_sum * b_2_add;
               CVector < _Flt >::CopyVector (b_2_add, p_a_2 + ibs2, p_a_sum + ibs);
               nzja_sum++;
               jp2++;
            }
         }
      }

   }

//
// Compute condensed block matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::CondenseMatrix (int _n_pt, const vector < _Int > &_ia_pt,
                              const vector < _Int > &_ja_pt, const vector < _Flt > &_a_pt,
                              int _b_size, vector < _Int > &_ia_cnd,
                              vector < _Int > &_ja_cnd, vector < _Flt > &_a_cnd)
   {

// Check block size

      if (_n_pt % _b_size != 0) {
         cout << " void CFct_bxb_impl<>::CondenseMatrix: error in block size !!! " <<
            endl;
         throw " void CFct_bxb_impl<>::CondenseMatrix: error in block size !!! ";
      }

      int nsup = _n_pt / _b_size;

// Open matrix data

      const _Int *p_ia_pt = NULL;
      const _Int *p_ja_pt = NULL;
      const _Flt *p_a_pt = NULL;

      if (_ia_pt.size () != 0)
         p_ia_pt = &_ia_pt[0];
      if (_ja_pt.size () != 0)
         p_ja_pt = &_ja_pt[0];
      if (_a_pt.size () != 0)
         p_a_pt = &_a_pt[0];

// Compute point to block array

      CVectorData < int >blk2pt (nsup + 1);
      CVectorData < int >pt2blk (_n_pt);

      int *pblk2pt = blk2pt.Ptr ();
      int *ppt2blk = pt2blk.Ptr ();

      int i, j;

      for (i = 0; i <= nsup; i++)
         pblk2pt[i] = i * _b_size;

      for (i = 0; i < nsup; i++) {
         for (j = pblk2pt[i]; j < pblk2pt[i + 1]; j++) {
            ppt2blk[j] = i;
         }
      }

// Compute block sparsity

      CVectorData < int >imaskblk (nsup);
      CVectorData < int >indblk (nsup);

      int *pimaskblk = imaskblk.Ptr ();
      int *pindblk = indblk.Ptr ();

      for (i = 0; i < nsup; i++)
         pimaskblk[i] = -1;

      int icycle_blk = -1;

      _ia_cnd.resize (nsup + 1);
      _Int *p_ia_cnd = &_ia_cnd[0];

      p_ia_cnd[0] = 0;

      int nzja_blk = 0;

      int k, kk, kk2;

      for (i = 0; i < nsup; i++) {
         icycle_blk++;
         for (j = pblk2pt[i]; j < pblk2pt[i + 1]; j++) {
            for (k = p_ia_pt[j]; k < p_ia_pt[j + 1]; k++) {
               kk = p_ja_pt[k];
               kk2 = ppt2blk[kk];
               if (pimaskblk[kk2] != icycle_blk) {
                  pimaskblk[kk2] = icycle_blk;
                  nzja_blk++;
               }
            }
         }
         p_ia_cnd[i + 1] = nzja_blk;
      }

// Store elems

      int b_2 = _b_size * _b_size;

      _ja_cnd.resize (nzja_blk + 1);
      _a_cnd.resize (nzja_blk * b_2 + 1);

      _Int *p_ja_cnd = &_ja_cnd[0];
      _Flt *p_a_cnd = &_a_cnd[0];

      CVector < _Flt >::SetByZeroes (nzja_blk * b_2, p_a_cnd);

      nzja_blk = 0;

      int nzja_blk0, ind, ishift, jshift;

      for (i = 0; i < nsup; i++) {
         icycle_blk++;
         nzja_blk0 = nzja_blk;
         for (j = pblk2pt[i]; j < pblk2pt[i + 1]; j++) {
            for (k = p_ia_pt[j]; k < p_ia_pt[j + 1]; k++) {
               kk = p_ja_pt[k];
               kk2 = ppt2blk[kk];
               if (pimaskblk[kk2] != icycle_blk) {
                  pimaskblk[kk2] = icycle_blk;
                  p_ja_cnd[nzja_blk] = kk2;
                  nzja_blk++;
               }
            }
         }
         sort (p_ja_cnd + nzja_blk0, p_ja_cnd + nzja_blk);
         for (j = nzja_blk0; j < nzja_blk; j++) {
            kk2 = p_ja_cnd[j];
            pindblk[kk2] = j;
         }
         for (j = pblk2pt[i]; j < pblk2pt[i + 1]; j++) {
            ishift = j - pblk2pt[i];
            for (k = p_ia_pt[j]; k < p_ia_pt[j + 1]; k++) {
               kk = p_ja_pt[k];
               kk2 = ppt2blk[kk];
               ind = pindblk[kk2];
               jshift = kk - kk2 * _b_size;
               p_a_cnd[ind * b_2 + jshift * _b_size + ishift] = p_a_pt[k];
            }
         }
      }

   }

//
// Compute condensed block rectangular matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::CondenseRectMatrix (int _nlist_pt, const vector < _Int > &_list_pt,
                                  const vector < _Int > &_ia_pt,
                                  const vector < _Int > &_ja_pt,
                                  const vector < _Flt > &_a_pt, int _b_size,
                                  int &_nlist_cnd, vector < _Int > &_list_cnd,
                                  vector < _Int > &_ia_cnd, vector < _Int > &_ja_cnd,
                                  vector < _Flt > &_a_cnd, int _nimax, int &_icycle,
                                  int *_imask)
   {

// Open work mask data

      int *pimask_th = _imask;
      int *plist_th = pimask_th + _nimax;
      int *pind_th = plist_th + _nimax;
      int *pind1_th = pind_th + _nimax;

// Open matrix data

      const _Int *p_list_pt = NULL;
      const _Int *p_ia_pt = NULL;
      const _Int *p_ja_pt = NULL;
      const _Flt *p_a_pt = NULL;

      if (_list_pt.size () != 0)
         p_list_pt = &_list_pt[0];
      if (_ia_pt.size () != 0)
         p_ia_pt = &_ia_pt[0];
      if (_ja_pt.size () != 0)
         p_ja_pt = &_ja_pt[0];
      if (_a_pt.size () != 0)
         p_a_pt = &_a_pt[0];

// Compute condensed list

      int nlist_cnd = 0;

      _icycle++;

      int i, jj, jblk;

      for (i = 0; i < _nlist_pt; i++) {
         jj = (int) p_list_pt[i];
         jblk = jj / _b_size;
         if (pimask_th[jblk] != _icycle) {
            pimask_th[jblk] = _icycle;
            plist_th[nlist_cnd] = jblk;
            nlist_cnd++;
         }
      }

      sort (plist_th, plist_th + nlist_cnd);

      _nlist_cnd = nlist_cnd;
      _list_cnd.resize (nlist_cnd + 1);

      _Int *p_list_cnd = &_list_cnd[0];

      for (i = 0; i < nlist_cnd; i++)
         p_list_cnd[i] = (_Int) plist_th[i];

// Index number for each column

      for (i = 0; i < nlist_cnd; i++) {
         jblk = (int) p_list_cnd[i];
         plist_th[jblk] = i;
      }

      int k;

      for (i = 0; i < _nlist_pt; i++) {
         jj = (int) p_list_pt[i];
         jblk = jj / _b_size;
         k = plist_th[jblk];
         pind_th[i] = k;
      }

// Compute ia_cnd array

      _ia_cnd.resize (nlist_cnd + 1);
      _Int *p_ia_cnd = &_ia_cnd[0];

      int nzja_cnd = 0;

      p_ia_cnd[0] = 0;

      int ip_beg = 0;
      int ilist_cnd = 0;

      int ip_end, jblk1, j;

      while (ip_beg < _nlist_pt) {
         ip_end = ip_beg;
         jblk = pind_th[ip_beg];
         while (ip_end + 1 < _nlist_pt) {
            jblk1 = pind_th[ip_end + 1];
            if (jblk1 == jblk) {
               ip_end++;
            } else {
               break;
            }
         }
         _icycle++;
         for (i = ip_beg; i <= ip_end; i++) {
            for (j = (int) p_ia_pt[i]; j < p_ia_pt[i + 1]; j++) {
               jj = (int) p_ja_pt[j];
               jblk = jj / _b_size;
               if (pimask_th[jblk] != _icycle) {
                  pimask_th[jblk] = _icycle;
                  nzja_cnd++;
               }
            }
         }
         p_ia_cnd[ilist_cnd + 1] = nzja_cnd;
         ilist_cnd++;
         ip_beg = ip_end + 1;
      }

// Allocate and fill ja_cnd, a_cnd

      int b_2 = _b_size * _b_size;

      _ja_cnd.resize (nzja_cnd + 1);
      _a_cnd.resize (nzja_cnd * b_2 + 1);

      _Int *p_ja_cnd = &_ja_cnd[0];
      _Flt *p_a_cnd = &_a_cnd[0];

      CVector < _Flt >::SetByZeroes (nzja_cnd * b_2, p_a_cnd);

      nzja_cnd = 0;

      p_ia_cnd[0] = 0;

      ip_beg = 0;

      int i0, j0, nlist_loc, iblk, ind, ilist_blkrow;

      while (ip_beg < _nlist_pt) {
         ip_end = ip_beg;
         jblk = pind_th[ip_beg];
         while (ip_end + 1 < _nlist_pt) {
            jblk1 = pind_th[ip_end + 1];
            if (jblk1 == jblk) {
               ip_end++;
            } else {
               break;
            }
         }
         _icycle++;
         nlist_loc = 0;
         for (i = ip_beg; i <= ip_end; i++) {
            for (j = (int) p_ia_pt[i]; j < p_ia_pt[i + 1]; j++) {
               jj = (int) p_ja_pt[j];
               jblk = jj / _b_size;
               if (pimask_th[jblk] != _icycle) {
                  plist_th[nlist_loc] = jblk;
                  nlist_loc++;
                  pimask_th[jblk] = _icycle;
               }
            }
         }
         sort (plist_th, plist_th + nlist_loc);
         for (i = 0; i < nlist_loc; i++) {
            jblk = plist_th[i];
            p_ja_cnd[nzja_cnd + i] = (_Int) jblk;
            pind1_th[jblk] = nzja_cnd + i;
         }
         for (i = ip_beg; i <= ip_end; i++) {
            jj = (int) p_list_pt[i];
            ilist_blkrow = pind_th[i];
            iblk = (int) p_list_cnd[ilist_blkrow];
            i0 = jj - iblk * _b_size;
            for (j = (int) p_ia_pt[i]; j < p_ia_pt[i + 1]; j++) {
               jj = (int) p_ja_pt[j];
               jblk = jj / _b_size;
               ind = pind1_th[jblk];
               j0 = jj - jblk * _b_size;
               p_a_cnd[ind * b_2 + j0 * _b_size + i0] = p_a_pt[j];
            }
         }
         nzja_cnd += nlist_loc;
         ip_beg = ip_end + 1;
      }
/*
      ofstream ffout ("ChkWrite.dat");
      for (i=0;i<_nlist_cnd;i++) {
         int irow = (int)p_list_cnd[i];
         for (j=(int)p_ia_cnd[i];j<p_ia_cnd[i+1];j++) {
            jj = (int)p_ja_cnd[j];
            if (jj == irow) {
               ffout << " Irow = " << irow << endl;
               PrintArray (ffout," Diag",b_2,p_a_cnd+j*b_2);
            }
         }
      }
*/
   }

//
// Compute ordered block matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::ReorderMatrix (bool _b_is_char, bool _b_is_pair, int _b_size, int _n,
                             const vector < int >&_order, const vector < _Int > &_ia_alu,
                             const vector < _Int > &_ja_alu,
                             const vector < char >&_ja_char_alu,
                             const vector < _Flt > &_a_alu, vector < _Int > &_ia_alu_ord,
                             vector < _Int > &_ja_alu_ord,
                             vector < char >&_ja_char_alu_ord,
                             vector < _Flt > &_a_alu_ord)
   {

// Compute inverse order

      vector < int >iord (_n + 1);

      int *piord = &iord[0];

      int i;

      for (i = 0; i < _n; i++)
         piord[_order[i]] = i;

// Reorder the matrix

      int b_2 = _b_size * _b_size;
      int b_2_2 = 2 * b_2;

      int b_2_ord = b_2;
      if (_b_is_pair)
         b_2_ord = b_2_2;

      int nzja = (int) _ia_alu[_n];

      _ia_alu_ord.resize (_n + 1);
      _ja_alu_ord.resize (nzja + 1);
      if (_b_is_char)
         _ja_char_alu_ord.resize (nzja + 1);
      _a_alu_ord.resize (nzja * b_2_ord + 1);

      _Int *piaord = &_ia_alu_ord[0];
      _Int *pjaord = &_ja_alu_ord[0];
      char *pjacharord = NULL;
      if (_b_is_char)
         pjacharord = &_ja_char_alu_ord[0];
      _Flt *paord = &_a_alu_ord[0];

      for (i = 0; i <= _n; i++)
         piaord[i] = 0;

      int nimax = 0;

      int ni;

      for (i = 0; i < _n; i++) {
         ni = (int) (_ia_alu[i + 1] - _ia_alu[i]);
         if (ni > nimax)
            nimax = ni;
      }

      vector < CSortInt > iiarr (nimax + 1);
      vector < char >charelems (nimax + 1);
      vector < _Flt > elems (nimax * b_2_ord + 1);

      CSortInt *piiarr = &iiarr[0];
      char *pcharelems = &charelems[0];
      _Flt *pelems = &elems[0];

      int j;

      for (i = 0; i < _n; i++) {
         j = _order[i];
         piaord[j + 1] = _ia_alu[i + 1] - _ia_alu[i];
      }

      for (i = 0; i < _n; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, nzloc, ibs, jold, jnew, k, ibs1, ibs2, ind;

      for (inew = 0; inew < _n; inew++) {
         i = piord[inew];
         nzloc = (int) (_ia_alu[i + 1] - _ia_alu[i]);
         ibs = nz;
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jold = (int) _ja_alu[j];
            jnew = _order[jold];
            pjaord[nz] = (_Int) jnew;
            if (_b_is_char)
               pjacharord[nz] = _ja_char_alu[j];
            ibs1 = nz * b_2_ord;
            ibs2 = j * b_2_ord;
            for (k = 0; k < b_2_ord; k++) {
               paord[ibs1 + k] = _a_alu[ibs2 + k];
            }
            nz++;
         }

// Sort elements

         for (j = 0; j < nzloc; j++) {
            piiarr[j].ival = (int) pjaord[ibs + j];
            piiarr[j].i2val = j;
         }

         sort (piiarr, piiarr + nzloc);

         if (_b_is_char) {
            for (j = 0; j < nzloc; j++)
               pcharelems[j] = pjacharord[ibs + j];
         }

         ibs1 = ibs * b_2_ord;

         for (j = 0; j < nzloc * b_2_ord; j++)
            pelems[j] = paord[ibs1 + j];

         for (j = 0; j < nzloc; j++) {
            pjaord[ibs + j] = (_Int) piiarr[j].ival;
            if (_b_is_char) {
               ind = piiarr[j].i2val;
               pjacharord[ibs + j] = pcharelems[ind];
            }
            ibs1 = (ibs + j) * b_2_ord;
            ibs2 = piiarr[j].i2val * b_2_ord;
            for (k = 0; k < b_2_ord; k++) {
               paord[ibs1 + k] = pelems[ibs2 + k];
            }
         }

      }

   }

//
// Compute block scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::ComputeScaling (int _sctype, int _nitersc, double _scl_min, int _b_size,
                              int _n, const vector < _Int > &_ia_alu,
                              const vector < _Int > &_ja_alu,
                              const vector < _Flt > &_a_alu, vector < _Flt > &_sclL,
                              vector < _Flt > &_sclL_inv, vector < _Flt > &_sclU,
                              vector < _Flt > &_sclU_inv, int &_nmodif,
                              double &_sclmin_att, double &_sclmax_att)
   {

      static const _Flt fone = (_Flt) 1.0e0;

// Allocate scl

      int b_2 = _b_size * _b_size;

      _sclL.resize (_n * b_2 + 1);
      _sclL_inv.resize (_n * b_2 + 1);
      _sclU.resize (_n * b_2 + 1);
      _sclU_inv.resize (_n * b_2 + 1);

      _Flt *p_sclL = &_sclL[0];
      _Flt *p_sclL_inv = &_sclL_inv[0];
      _Flt *p_sclU = &_sclU[0];
      _Flt *p_sclU_inv = &_sclU_inv[0];

// Open structures

      const _Int *p_ia_alu = NULL;
      const _Int *p_ja_alu = NULL;
      const _Flt *p_a_alu = NULL;

      if (_ia_alu.size () > 0)
         p_ia_alu = &_ia_alu[0];
      if (_ja_alu.size () > 0)
         p_ja_alu = &_ja_alu[0];
      if (_a_alu.size () > 0)
         p_a_alu = &_a_alu[0];

// Simple diagonal based block scaling

      _sclmin_att = 1.0e100;
      _sclmax_att = -1.0e100;

      int i, j, jj, ibs, jbs;

      CVectorData < _Flt > Work (20 * b_2);
      CVectorData < double >dWork (20 * b_2);

      _Flt *pWork = Work.Ptr ();
      double *pdWork = dWork.Ptr ();

      _nmodif = 0;

      if (_sctype == -1) {

         CVector < _Flt >::SetByZeroes (_n * b_2, p_sclL);
         CVector < _Flt >::SetByZeroes (_n * b_2, p_sclL_inv);
         CVector < _Flt >::SetByZeroes (_n * b_2, p_sclU);
         CVector < _Flt >::SetByZeroes (_n * b_2, p_sclU_inv);

         for (i = 0; i < _n; i++) {
            ibs = i * b_2;
            for (j = 0; j < _b_size; j++) {
               p_sclL[ibs + j * _b_size + j] = fone;
               p_sclL_inv[ibs + j * _b_size + j] = fone;
               p_sclU[ibs + j * _b_size + j] = fone;
               p_sclU_inv[ibs + j * _b_size + j] = fone;
            }
         }

         _sclmin_att = fone;
         _sclmax_att = fone;

      } else if (_sctype == 0 || _sctype == 1) {

         bool b_found;
         double sclmin_temp, sclmax_temp;

         _sclmin_att = 1.0e100;
         _sclmax_att = -1.0e100;

         for (i = 0; i < _n; i++) {

            b_found = false;
            ibs = i * b_2;

            for (j = (int) p_ia_alu[i]; j < p_ia_alu[i + 1]; j++) {
               jj = (int) p_ja_alu[j];
               if (i == jj) {
                  b_found = true;
                  jbs = j * b_2;
                  if (_sctype == 0) {
                     _nmodif +=
                        CBlock_BxB_traits < _Flt >::FctInv_BxB (_scl_min, _b_size,
                                                                p_a_alu + jbs,
                                                                p_sclL_inv + ibs,
                                                                p_sclU_inv + ibs,
                                                                p_sclL + ibs,
                                                                p_sclU + ibs, sclmin_temp,
                                                                sclmax_temp);
                  } else {
                     _nmodif +=
                        CBlock_BxB_traits < _Flt >::FctInv_Svd_BxB (_scl_min, _b_size,
                                                                    p_a_alu + jbs,
                                                                    p_sclL_inv + ibs,
                                                                    p_sclU_inv + ibs,
                                                                    p_sclL + ibs,
                                                                    p_sclU + ibs, pWork,
                                                                    pdWork, sclmin_temp,
                                                                    sclmax_temp);
                  }
                  if (sclmin_temp < _sclmin_att)
                     _sclmin_att = sclmin_temp;
                  if (sclmax_temp > _sclmax_att)
                     _sclmax_att = sclmax_temp;
               }
            }
            if (!b_found) {
               CVector < _Flt >::SetByZeroes (b_2, p_sclL + ibs);
               CVector < _Flt >::SetByZeroes (b_2, p_sclL_inv + ibs);
               CVector < _Flt >::SetByZeroes (b_2, p_sclU + ibs);
               CVector < _Flt >::SetByZeroes (b_2, p_sclU_inv + ibs);
               for (j = 0; j < _b_size; j++) {
                  p_sclL[ibs + j * _b_size + j] = fone;
                  p_sclL_inv[ibs + j * _b_size + j] = fone;
                  p_sclU[ibs + j * _b_size + j] = fone;
                  p_sclU_inv[ibs + j * _b_size + j] = fone;
               }
            }
         }

      } else {

// Iterative rows/columns block balancing scaling

// Compute transposed sparsity and reference to initial

         vector < _Int > ia_alu_t;
         vector < _Int > ja_alu_t;

         CFct_impl < _Int, _Flt >::TransposeSp (_n, _ia_alu, _ja_alu, ia_alu_t, ja_alu_t);

         _Int *pia_alu_t = NULL;
         _Int *pja_alu_t = NULL;

         if (ia_alu_t.size () > 0)
            pia_alu_t = &ia_alu_t[0];
         if (ja_alu_t.size () > 0)
            pja_alu_t = &ja_alu_t[0];

         int nzja_loc = (int) ja_alu_t.size ();

         CVectorData < _Int > iptr (_n);
         CVectorData < _Int > indt2ini (nzja_loc);

         _Int *piptr = iptr.Ptr ();
         _Int *pindt2ini = indt2ini.Ptr ();

         int i, j, jj, k;

         for (i = 0; i < _n; i++)
            piptr[i] = p_ia_alu[i];

         for (i = 0; i < _n; i++) {
            for (j = pia_alu_t[i]; j < pia_alu_t[i + 1]; j++) {
               jj = pja_alu_t[j];
               k = (int) piptr[jj];
               pindt2ini[j] = k;
               piptr[jj]++;
            }
         }

// Compute preliminary columns block scaling based on columns norms

         int ntot = _n * _b_size;

         CVectorData < _Flt > fnrm2 (ntot);
         _Flt *pfnrm2 = fnrm2.Ptr ();

         CVector < _Flt >::SetByZeroes (ntot, pfnrm2);

         int kii, kjj;
         const _Flt *pelems;

         _Flt aux;

         _sclmin_att = 1.0e100;
         _sclmax_att = -1.0e100;

         for (i = 0; i < _n; i++) {
            for (j = (int) p_ia_alu[i]; j < p_ia_alu[i + 1]; j++) {
               jj = (int) p_ja_alu[j];
               ibs = jj * _b_size;
               pelems = p_a_alu + j * b_2;
               for (kjj = 0; kjj < _b_size; kjj++) {
                  for (kii = 0; kii < _b_size; kii++) {
                     aux = pelems[kjj * _b_size + kii];
                     pfnrm2[ibs + kjj] += aux * aux;
                  }
               }
            }
         }

         CVector < _Flt >::SetByZeroes (_n * b_2, p_sclU);
         CVector < _Flt >::SetByZeroes (_n * b_2, p_sclU_inv);

         for (i = 0; i < _n; i++) {
            for (j = 0; j < _b_size; j++) {
               k = i * _b_size + j;
               aux = pfnrm2[k];
               aux = sqrt (aux);
               if (aux < _sclmin_att)
                  _sclmin_att = (double) aux;
               if (aux > _sclmax_att)
                  _sclmax_att = (double) aux;
               if (aux < _scl_min)
                  aux = (_Flt) _scl_min;
               aux = sqrt (aux);
               p_sclU_inv[i * b_2 + j * _b_size + j] = aux;
               aux = fone / aux;
               p_sclU[i * b_2 + j * _b_size + j] = aux;
            }
         }

// Prepare scaling support structures

         int nimax = 0;

         int niloc;

         for (i = 0; i < _n; i++) {
            niloc = p_ia_alu[i + 1] - p_ia_alu[i];
            if (niloc > nimax)
               nimax = niloc;
            niloc = pia_alu_t[i + 1] - pia_alu_t[i];
            if (niloc > nimax)
               nimax = niloc;
         }

         nimax++;

         CVectorData < _Flt > y_blk (nimax * b_2);
         CVectorData < _Flt > z_blk (nimax * b_2);
         CVectorData < _Flt > r_matr (b_2);
         CVectorData < _Flt > tau_matr (b_2);

         _Flt *py_blk = y_blk.Ptr ();
         _Flt *pz_blk = z_blk.Ptr ();
         _Flt *pr_matr = r_matr.Ptr ();
         _Flt *ptau_matr = tau_matr.Ptr ();

// Main scaling cycle

         int icycle, nloc, ind;

         for (icycle = 0; icycle < _nitersc; icycle++) {

// Collect data of each block row and compute its Qrd

            for (i = 0; i < _n; i++) {

               niloc = 0;

               for (j = p_ia_alu[i]; j < p_ia_alu[i + 1]; j++) {
                  jj = p_ja_alu[j];
                  CBlock_BxB_traits < _Flt >::MM_BxB (_b_size, p_a_alu + j * b_2,
                                                      p_sclU + jj * b_2,
                                                      py_blk + niloc * b_2);
                  niloc++;
               }

               nloc = niloc * _b_size;

               CVector < _Flt >::TransposeBlock (_b_size, nloc, py_blk, _b_size, pz_blk,
                                                 nloc);

               CVector < _Flt >::QrdBlock (_b_size, nloc, pz_blk, nloc, ptau_matr);

// Store R factor and its inverse

               CVector < _Flt >::SetByZeroes (b_2, pr_matr);

               for (kjj = 0; kjj < _b_size; kjj++) {
                  for (kii = 0; kii <= kjj; kii++) {
                     pr_matr[kjj * _b_size + kii] = pz_blk[kjj * nloc + kii];
                  }
               }

               CVector < _Flt >::CopyVector (b_2, pr_matr, p_sclL_inv + i * b_2);

               CBlock_BxB_traits < _Flt >::InvU_BxB (_b_size, pr_matr, p_sclL + i * b_2);

// Check block row scaling 

               if (false) {

                  niloc = 0;

                  for (j = p_ia_alu[i]; j < p_ia_alu[i + 1]; j++) {
                     jj = p_ja_alu[j];
                     CBlock_BxB_traits < _Flt >::MM_BxB (_b_size, p_a_alu + j * b_2,
                                                         p_sclU + jj * b_2, pr_matr);
                     CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, p_sclL + i * b_2,
                                                          pr_matr, py_blk + niloc * b_2);
                     niloc++;
                  }

                  nloc = niloc * _b_size;

                  CVector < _Flt >::TransposeBlock (_b_size, nloc, py_blk, _b_size,
                                                    pz_blk, nloc);

                  CVector < _Flt >::QrdBlock (_b_size, nloc, pz_blk, nloc, ptau_matr);

// Store R factor and its inverse

                  CVector < _Flt >::SetByZeroes (b_2, pr_matr);

                  for (kjj = 0; kjj < _b_size; kjj++) {
                     for (kii = 0; kii <= kjj; kii++) {
                        pr_matr[kjj * _b_size + kii] = pz_blk[kjj * nloc + kii];
                     }
                  }

               }
// Collect data of each block column and compute its QR

               niloc = 0;

               for (j = pia_alu_t[i]; j < pia_alu_t[i + 1]; j++) {
                  jj = (int) pja_alu_t[j];
                  ind = pindt2ini[j];
                  CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, p_a_alu + ind * b_2,
                                                       p_sclL + jj * b_2,
                                                       py_blk + niloc * b_2);
                  niloc++;
               }

               nloc = niloc * _b_size;

               CVector < _Flt >::TransposeBlock (_b_size, nloc, py_blk, _b_size, pz_blk,
                                                 nloc);

               CVector < _Flt >::QrdBlock (_b_size, nloc, pz_blk, nloc, ptau_matr);

// Store R factor and its inverse

               CVector < _Flt >::SetByZeroes (b_2, pr_matr);

               for (kjj = 0; kjj < _b_size; kjj++) {
                  for (kii = 0; kii <= kjj; kii++) {
                     pr_matr[kjj * _b_size + kii] = pz_blk[kjj * nloc + kii];
                  }
               }

               CVector < _Flt >::CopyVector (b_2, pr_matr, p_sclU_inv + i * b_2);

               CBlock_BxB_traits < _Flt >::InvU_BxB (_b_size, pr_matr, p_sclU + i * b_2);

            }

         }

      }

   }

//
// Perform explicit block scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::MatrixScale (int _b_size, int _n, const vector < _Flt > &_sclL,
                           const vector < _Flt > &_sclU, const vector < _Int > &_ia_alu,
                           const vector < _Int > &_ja_alu, vector < _Flt > &_a_alu)
   {

// Open structures

      int b_2 = _b_size * _b_size;

      const _Flt *p_sclL = NULL;
      const _Flt *p_sclU = NULL;

      if (_sclL.size () > 0)
         p_sclL = &_sclL[0];
      if (_sclU.size () > 0)
         p_sclU = &_sclU[0];

      const _Int *p_ia_alu = NULL;
      const _Int *p_ja_alu = NULL;
      _Flt *p_a_alu = NULL;

      if (_ia_alu.size () > 0)
         p_ia_alu = &_ia_alu[0];
      if (_ja_alu.size () > 0)
         p_ja_alu = &_ja_alu[0];
      if (_a_alu.size () > 0)
         p_a_alu = &_a_alu[0];

      int i, j, jj, ibs, jbs, kbs;

      CVectorData < _Flt > lu_temp (b_2);
      _Flt *plu_temp = lu_temp.Ptr ();

      for (i = 0; i < _n; i++) {
         ibs = i * b_2;
         for (j = (int) p_ia_alu[i]; j < p_ia_alu[i + 1]; j++) {
            jj = (int) p_ja_alu[j];
            jbs = jj * b_2;
            kbs = j * b_2;

            CBlock_BxB_traits < _Flt >::MM_BxB (_b_size, p_a_alu + kbs, p_sclU + jbs,
                                                plu_temp);
            CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, p_sclL + ibs, plu_temp,
                                                 p_a_alu + kbs);

         }
      }

   }

//
// Perform explicit block scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::MatrixScale (int _b_size, int _nlist, const _Flt * _sclL,
                           const _Flt * _sclU, const vector < _Int > &_list_alu,
                           const vector < _Int > &_ia_alu, const vector < _Int > &_ja_alu,
                           vector < _Flt > &_a_alu)
   {

// Open structures

      int b_2 = _b_size * _b_size;

      const _Flt *p_sclL = _sclL;
      const _Flt *p_sclU = _sclU;

      const _Int *p_list_alu = NULL;
      const _Int *p_ia_alu = NULL;
      const _Int *p_ja_alu = NULL;
      _Flt *p_a_alu = NULL;

      if (_list_alu.size () > 0)
         p_list_alu = &_list_alu[0];
      if (_ia_alu.size () > 0)
         p_ia_alu = &_ia_alu[0];
      if (_ja_alu.size () > 0)
         p_ja_alu = &_ja_alu[0];
      if (_a_alu.size () > 0)
         p_a_alu = &_a_alu[0];

      int ilist, i, j, jj, ibs, jbs, kbs;

      CVectorData < _Flt > lu_temp (b_2);
      _Flt *plu_temp = lu_temp.Ptr ();

      for (ilist = 0; ilist < _nlist; ilist++) {
         i = (int) p_list_alu[ilist];
         ibs = i * b_2;
         for (j = (int) p_ia_alu[ilist]; j < p_ia_alu[ilist + 1]; j++) {
            jj = (int) p_ja_alu[j];
            jbs = jj * b_2;
            kbs = j * b_2;

            CBlock_BxB_traits < _Flt >::MM_BxB (_b_size, p_a_alu + kbs, p_sclU + jbs,
                                                plu_temp);
            CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, p_sclL + ibs, plu_temp,
                                                 p_a_alu + kbs);

         }
      }

   }

//
// Split matrix data into L and U parts
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::SplitLU (bool _b_is_char, int _b_size, int _n,
                       const vector < _Int > &_ia_alu, const vector < _Int > &_ja_alu,
                       const vector < char >&_jachar_alu, const vector < _Flt > &_a_alu,
                       vector < _Int > &_ia_l, vector < _Int > &_ja_l,
                       vector < char >&_jachar_l, vector < _Flt > &_a_l,
                       vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                       vector < char >&_jachar_u, vector < _Flt > &_a_u)
   {

// Open structures

      int b_2 = _b_size * _b_size;

      const _Int *p_ia_alu = NULL;
      const _Int *p_ja_alu = NULL;
      const char *p_jachar_alu = NULL;
      const _Flt *p_a_alu = NULL;

      if (_ia_alu.size () > 0)
         p_ia_alu = &_ia_alu[0];
      if (_ja_alu.size () > 0)
         p_ja_alu = &_ja_alu[0];
      if (_b_is_char && _jachar_alu.size () > 0)
         p_jachar_alu = &_jachar_alu[0];
      if (_a_alu.size () > 0)
         p_a_alu = &_a_alu[0];

// Count number of elems

      int nzja_l = 0;
      int nzja_u = 0;

      int i, j, jj;

      bool b_found = false;

      int n_zero_diag = 0;

      for (i = 0; i < _n; i++) {
         b_found = false;
         for (j = (int) p_ia_alu[i]; j < p_ia_alu[i + 1]; j++) {
            jj = (int) p_ja_alu[j];
            if (jj <= i)
               nzja_l++;
            if (jj >= i)
               nzja_u++;
            if (jj == i)
               b_found = true;
         }
         if (!b_found) {
            n_zero_diag++;
            nzja_l++;
            nzja_u++;
         }
      }

      if (n_zero_diag > 0)
         cout << " CFct_bxb_impl <>::SplitLU: N_zero_diag = " << n_zero_diag << endl;

// Allocate and fill

      _ia_l.resize (_n + 1);
      _ja_l.resize (nzja_l + 1);
      _a_l.resize (nzja_l * b_2 + 1);

      _ia_u.resize (_n + 1);
      _ja_u.resize (nzja_u + 1);
      _a_u.resize (nzja_u * b_2 + 1);

      if (_b_is_char) {
         _jachar_l.resize (nzja_l + 1);
         _jachar_u.resize (nzja_u + 1);
      }

      _Int *p_ia_l = NULL;
      _Int *p_ja_l = NULL;
      char *p_jachar_l = NULL;
      _Flt *p_a_l = NULL;

      if (_ia_l.size () > 0)
         p_ia_l = &_ia_l[0];
      if (_ja_l.size () > 0)
         p_ja_l = &_ja_l[0];
      if (_b_is_char && _jachar_l.size () > 0)
         p_jachar_l = &_jachar_l[0];
      if (_a_l.size () > 0)
         p_a_l = &_a_l[0];

      _Int *p_ia_u = NULL;
      _Int *p_ja_u = NULL;
      char *p_jachar_u = NULL;
      _Flt *p_a_u = NULL;

      if (_ia_u.size () > 0)
         p_ia_u = &_ia_u[0];
      if (_ja_u.size () > 0)
         p_ja_u = &_ja_u[0];
      if (_b_is_char && _jachar_u.size () > 0)
         p_jachar_u = &_jachar_u[0];
      if (_a_u.size () > 0)
         p_a_u = &_a_u[0];

      nzja_l = 0;
      nzja_u = 0;

      p_ia_l[0] = 0;
      p_ia_u[0] = 0;

      for (i = 0; i < _n; i++) {
         b_found = false;
         for (j = (int) p_ia_alu[i]; j < p_ia_alu[i + 1]; j++) {
            jj = (int) p_ja_alu[j];
            if (jj == i)
               b_found = true;
         }
         if (!b_found) {
            p_ja_u[nzja_u] = i;
            if (_b_is_char) {
               p_jachar_u[nzja_u] = 0;
            }
            CVector < _Flt >::SetByZeroes (b_2, p_a_u + nzja_u * b_2);
            nzja_u++;
         }
         for (j = (int) p_ia_alu[i]; j < p_ia_alu[i + 1]; j++) {
            jj = (int) p_ja_alu[j];
            if (jj <= i) {
               p_ja_l[nzja_l] = (_Int) jj;
               if (_b_is_char) {
                  p_jachar_l[nzja_l] = p_jachar_alu[j];
               }
               CVector < _Flt >::CopyVector (b_2, p_a_alu + j * b_2,
                                             p_a_l + nzja_l * b_2);
               nzja_l++;
            }
            if (jj >= i) {
               p_ja_u[nzja_u] = (_Int) jj;
               if (_b_is_char) {
                  p_jachar_u[nzja_u] = p_jachar_alu[j];
               }
               CVector < _Flt >::CopyVector (b_2, p_a_alu + j * b_2,
                                             p_a_u + nzja_u * b_2);
               nzja_u++;
            }
         }
         if (!b_found) {
            p_ja_l[nzja_l] = i;
            if (_b_is_char) {
               p_jachar_l[nzja_l] = 0;
            }
            CVector < _Flt >::SetByZeroes (b_2, p_a_l + nzja_l * b_2);
            nzja_l++;
         }
         p_ia_l[i + 1] = nzja_l;
         p_ia_u[i + 1] = nzja_u;
      }

   }

//
// Transpose square block matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::Transpose (bool _b_is_char, int _b_size, int _n,
                         const vector < _Int > &_ia_a, const vector < _Int > &_ja_a,
                         const vector < char >&_jachar_a, const vector < _Flt > &_a_a,
                         vector < _Int > &_ia_at, vector < _Int > &_ja_at,
                         vector < char >&_jachar_at, vector < _Flt > &_a_at)
   {

// Open structures

      int b_2 = _b_size * _b_size;

      const _Int *p_ia_a = NULL;
      const _Int *p_ja_a = NULL;
      const char *p_jachar_a = NULL;
      const _Flt *p_a_a = NULL;

      if (_ia_a.size () > 0)
         p_ia_a = &_ia_a[0];
      if (_ja_a.size () > 0)
         p_ja_a = &_ja_a[0];
      if (_b_is_char && _jachar_a.size () > 0)
         p_jachar_a = &_jachar_a[0];
      if (_a_a.size () > 0)
         p_a_a = &_a_a[0];

      int nzja_a = (int) p_ia_a[_n];

// Allocate work array

      CVectorData < _Int > iptr (_n + 1);
      _Int *piptr = iptr.Ptr ();

// Allocate tranposed data

      _ia_at.resize (_n + 1);
      _ja_at.resize (nzja_a + 1);
      if (_b_is_char) {
         _jachar_at.resize (nzja_a + 1);
      }
      _a_at.resize (nzja_a * b_2 + 1);

      _Int *p_ia_at = NULL;
      _Int *p_ja_at = NULL;
      char *p_jachar_at = NULL;
      _Flt *p_a_at = NULL;

      if (_ia_at.size () > 0)
         p_ia_at = &_ia_at[0];
      if (_ja_at.size () > 0)
         p_ja_at = &_ja_at[0];
      if (_b_is_char && _jachar_at.size () > 0)
         p_jachar_at = &_jachar_at[0];
      if (_a_at.size () > 0)
         p_a_at = &_a_at[0];

      int i, j, jj, k;

      for (i = 0; i <= _n; i++)
         p_ia_at[i] = 0;

      for (i = 0; i < nzja_a; i++) {
         jj = (int) p_ja_a[i];
         p_ia_at[jj + 1]++;
      }

      for (i = 0; i < _n; i++)
         p_ia_at[i + 1] = p_ia_at[i] + p_ia_at[i + 1];

      for (i = 0; i < _n; i++)
         piptr[i] = p_ia_at[i];

      for (i = 0; i < _n; i++) {
         for (j = (int) p_ia_a[i]; j < p_ia_a[i + 1]; j++) {
            jj = (int) p_ja_a[j];
            k = (int) piptr[jj];
            p_ja_at[k] = i;
            if (_b_is_char) {
               p_jachar_at[k] = p_jachar_a[j];
            }
            CBlock_BxB_traits < _Flt >::Transp_BxB (_b_size, p_a_a + j * b_2,
                                                    p_a_at + k * b_2);
            piptr[jj]++;
         }
      }

   }

//
// Combine block L and U
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::CombineLU (bool _b_is_char, int _b_size, int _n,
                         const vector < _Int > &_ia_al, const vector < _Int > &_ja_al,
                         const vector < char >&_jachar_al, const vector < _Flt > &_a_al,
                         const vector < _Int > &_ia_au, const vector < _Int > &_ja_au,
                         const vector < char >&_jachar_au, const vector < _Flt > &_a_au,
                         vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                         vector < char >&_jachar_alu, vector < _Flt > &_a_alu)
   {

// Open structures

      int b_2 = _b_size * _b_size;

      const _Int *p_ia_l = NULL;
      const _Int *p_ja_l = NULL;
      const char *p_jachar_l = NULL;
      const _Flt *p_a_l = NULL;

      if (_ia_al.size () > 0)
         p_ia_l = &_ia_al[0];
      if (_ja_al.size () > 0)
         p_ja_l = &_ja_al[0];
      if (_b_is_char && _jachar_al.size () > 0)
         p_jachar_l = &_jachar_al[0];
      if (_a_al.size () > 0)
         p_a_l = &_a_al[0];

      const _Int *p_ia_u = NULL;
      const _Int *p_ja_u = NULL;
      const char *p_jachar_u = NULL;
      const _Flt *p_a_u = NULL;

      if (_ia_au.size () > 0)
         p_ia_u = &_ia_au[0];
      if (_ja_au.size () > 0)
         p_ja_u = &_ja_au[0];
      if (_b_is_char && _jachar_au.size () > 0)
         p_jachar_u = &_jachar_au[0];
      if (_a_au.size () > 0)
         p_a_u = &_a_au[0];

// Compute number of extended elems

      int nzja_ext = 0;

      int i, ipl, ipu, iendl, iendu, jj_l, jj_u;

      for (i = 0; i < _n; i++) {
         ipl = (int) p_ia_l[i];
         ipu = (int) p_ia_u[i];
         iendl = (int) p_ia_l[i + 1] - 1;
         iendu = (int) p_ia_u[i + 1] - 1;
         while (ipl <= iendl || ipu <= iendu) {
            if (ipl <= iendl && ipu <= iendu) {
               jj_l = (int) p_ja_l[ipl];
               jj_u = (int) p_ja_u[ipu];
               if (jj_l == jj_u) {
                  nzja_ext++;
                  ipl++;
                  ipu++;
               } else if (jj_l < jj_u) {
                  nzja_ext++;
                  ipl++;
               } else {
                  nzja_ext++;
                  ipu++;
               }
            } else if (ipl <= iendl) {
               nzja_ext++;
               ipl++;
            } else {
               nzja_ext++;
               ipu++;
            }
         }
      }

// Count number of elems

      _ia_alu.resize (_n + 1);
      _ja_alu.resize (nzja_ext + 1);
      if (_b_is_char) {
         _jachar_alu.resize (nzja_ext + 1);
      }
      _a_alu.resize (nzja_ext * b_2 + 1);

      _Int *p_ia_alu = NULL;
      _Int *p_ja_alu = NULL;
      char *p_jachar_alu = NULL;
      _Flt *p_a_alu = NULL;

      if (_ia_alu.size () > 0)
         p_ia_alu = &_ia_alu[0];
      if (_ja_alu.size () > 0)
         p_ja_alu = &_ja_alu[0];
      if (_b_is_char && _jachar_alu.size () > 0)
         p_jachar_alu = &_jachar_alu[0];
      if (_a_alu.size () > 0)
         p_a_alu = &_a_alu[0];

      p_ia_alu[0] = 0;
      nzja_ext = 0;

      for (i = 0; i < _n; i++) {
         ipl = (int) p_ia_l[i];
         ipu = (int) p_ia_u[i];
         iendl = (int) p_ia_l[i + 1] - 1;
         iendu = (int) p_ia_u[i + 1] - 1;
         while (ipl <= iendl || ipu <= iendu) {
            if (ipl <= iendl && ipu <= iendu) {
               jj_l = (int) p_ja_l[ipl];
               jj_u = (int) p_ja_u[ipu];
               if (jj_l == jj_u) {
                  p_ja_alu[nzja_ext] = jj_l;
                  if (_b_is_char) {
                     p_jachar_alu[nzja_ext] = p_jachar_l[ipl];
                  }
                  CVector < _Flt >::CopyVector (b_2, p_a_l + ipl * b_2,
                                                p_a_alu + nzja_ext * b_2);
                  nzja_ext++;
                  ipl++;
                  ipu++;
               } else if (jj_l < jj_u) {
                  p_ja_alu[nzja_ext] = jj_l;
                  if (_b_is_char) {
                     p_jachar_alu[nzja_ext] = p_jachar_l[ipl];
                  }
                  CVector < _Flt >::CopyVector (b_2, p_a_l + ipl * b_2,
                                                p_a_alu + nzja_ext * b_2);
                  nzja_ext++;
                  ipl++;
               } else {
                  p_ja_alu[nzja_ext] = jj_u;
                  if (_b_is_char) {
                     p_jachar_alu[nzja_ext] = p_jachar_u[ipu];
                  }
                  CVector < _Flt >::CopyVector (b_2, p_a_u + ipu * b_2,
                                                p_a_alu + nzja_ext * b_2);
                  nzja_ext++;
                  ipu++;
               }
            } else if (ipl <= iendl) {
               p_ja_alu[nzja_ext] = p_ja_l[ipl];
               if (_b_is_char) {
                  p_jachar_alu[nzja_ext] = p_jachar_l[ipl];
               }
               CVector < _Flt >::CopyVector (b_2, p_a_l + ipl * b_2,
                                             p_a_alu + nzja_ext * b_2);
               nzja_ext++;
               ipl++;
            } else {
               p_ja_alu[nzja_ext] = p_ja_u[ipu];
               if (_b_is_char) {
                  p_jachar_alu[nzja_ext] = p_jachar_u[ipu];
               }
               CVector < _Flt >::CopyVector (b_2, p_a_u + ipu * b_2,
                                             p_a_alu + nzja_ext * b_2);
               nzja_ext++;
               ipu++;
            }
         }
         p_ia_alu[i + 1] = nzja_ext;
      }

   }

//
// Combine block L and U data into extended block pairs
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::CombineLUPairs (bool _b_is_char, int _b_size, int _n,
                              const vector < _Int > &_ia_l, const vector < _Int > &_ja_l,
                              const vector < char >&_jachar_l,
                              const vector < _Flt > &_a_l, const vector < _Int > &_ia_u,
                              const vector < _Int > &_ja_u,
                              const vector < char >&_jachar_u,
                              const vector < _Flt > &_a_u, vector < _Int > &_ia_alu,
                              vector < _Int > &_ja_alu, vector < char >&_jachar_alu,
                              vector < _Flt > &_a_alu)
   {

// Open structures

      int b_2 = _b_size * _b_size;

      const _Int *p_ia_l = NULL;
      const _Int *p_ja_l = NULL;
      const char *p_jachar_l = NULL;
      const _Flt *p_a_l = NULL;

      if (_ia_l.size () > 0)
         p_ia_l = &_ia_l[0];
      if (_ja_l.size () > 0)
         p_ja_l = &_ja_l[0];
      if (_b_is_char && _jachar_l.size () > 0)
         p_jachar_l = &_jachar_l[0];
      if (_a_l.size () > 0)
         p_a_l = &_a_l[0];

      const _Int *p_ia_u = NULL;
      const _Int *p_ja_u = NULL;
      const char *p_jachar_u = NULL;
      const _Flt *p_a_u = NULL;

      if (_ia_u.size () > 0)
         p_ia_u = &_ia_u[0];
      if (_ja_u.size () > 0)
         p_ja_u = &_ja_u[0];
      if (_b_is_char && _jachar_u.size () > 0)
         p_jachar_u = &_jachar_u[0];
      if (_a_u.size () > 0)
         p_a_u = &_a_u[0];

// Compute number of extended elems

      int nzja_ext = 0;

      int i, ipl, ipu, iendl, iendu, jj_l, jj_u;

      for (i = 0; i < _n; i++) {
         ipl = (int) p_ia_l[i];
         ipu = (int) p_ia_u[i];
         iendl = (int) p_ia_l[i + 1] - 1;
         iendu = (int) p_ia_u[i + 1] - 1;
         while (ipl <= iendl || ipu <= iendu) {
            if (ipl <= iendl && ipu <= iendu) {
               jj_l = (int) p_ja_l[ipl];
               jj_u = (int) p_ja_u[ipu];
               if (jj_l == jj_u) {
                  nzja_ext++;
                  ipl++;
                  ipu++;
               } else if (jj_l < jj_u) {
                  nzja_ext++;
                  ipl++;
               } else {
                  nzja_ext++;
                  ipu++;
               }
            } else if (ipl <= iendl) {
               nzja_ext++;
               ipl++;
            } else {
               nzja_ext++;
               ipu++;
            }
         }
      }

// Count number of elems

      _ia_alu.resize (_n + 1);
      _ja_alu.resize (nzja_ext + 1);
      if (_b_is_char) {
         _jachar_alu.resize (nzja_ext + 1);
      }
      _a_alu.resize (nzja_ext * 2 * b_2 + 1);

      _Int *p_ia_alu = NULL;
      _Int *p_ja_alu = NULL;
      char *p_jachar_alu = NULL;
      _Flt *p_a_alu = NULL;

      if (_ia_alu.size () > 0)
         p_ia_alu = &_ia_alu[0];
      if (_ja_alu.size () > 0)
         p_ja_alu = &_ja_alu[0];
      if (_b_is_char && _jachar_alu.size () > 0)
         p_jachar_alu = &_jachar_alu[0];
      if (_a_alu.size () > 0)
         p_a_alu = &_a_alu[0];

      p_ia_alu[0] = 0;
      nzja_ext = 0;

      for (i = 0; i < _n; i++) {
         ipl = (int) p_ia_l[i];
         ipu = (int) p_ia_u[i];
         iendl = (int) p_ia_l[i + 1] - 1;
         iendu = (int) p_ia_u[i + 1] - 1;
         while (ipl <= iendl || ipu <= iendu) {
            if (ipl <= iendl && ipu <= iendu) {
               jj_l = (int) p_ja_l[ipl];
               jj_u = (int) p_ja_u[ipu];
               if (jj_l == jj_u) {
                  p_ja_alu[nzja_ext] = jj_l;
                  if (_b_is_char) {
                     p_jachar_alu[nzja_ext] = p_jachar_l[ipl];
                  }

                  CVector < _Flt >::CopyVector (b_2, p_a_l + ipl * b_2,
                                                p_a_alu + (nzja_ext * 2) * b_2);
                  CVector < _Flt >::CopyVector (b_2, p_a_u + ipu * b_2,
                                                p_a_alu + (nzja_ext * 2 + 1) * b_2);

                  nzja_ext++;
                  ipl++;
                  ipu++;
               } else if (jj_l < jj_u) {
                  p_ja_alu[nzja_ext] = jj_l;
                  if (_b_is_char) {
                     p_jachar_alu[nzja_ext] = p_jachar_l[ipl];
                  }

                  CVector < _Flt >::CopyVector (b_2, p_a_l + ipl * b_2,
                                                p_a_alu + (nzja_ext * 2) * b_2);
                  CVector < _Flt >::SetByZeroes (b_2, p_a_alu + (nzja_ext * 2 + 1) * b_2);

                  nzja_ext++;
                  ipl++;
               } else {
                  p_ja_alu[nzja_ext] = jj_u;
                  if (_b_is_char) {
                     p_jachar_alu[nzja_ext] = p_jachar_u[ipu];
                  }

                  CVector < _Flt >::SetByZeroes (b_2, p_a_alu + (nzja_ext * 2) * b_2);
                  CVector < _Flt >::CopyVector (b_2, p_a_u + ipu * b_2,
                                                p_a_alu + (nzja_ext * 2 + 1) * b_2);

                  nzja_ext++;
                  ipu++;
               }
            } else if (ipl <= iendl) {
               p_ja_alu[nzja_ext] = p_ja_l[ipl];
               if (_b_is_char) {
                  p_jachar_alu[nzja_ext] = p_jachar_l[ipl];
               }

               CVector < _Flt >::CopyVector (b_2, p_a_l + ipl * b_2,
                                             p_a_alu + (nzja_ext * 2) * b_2);
               CVector < _Flt >::SetByZeroes (b_2, p_a_alu + (nzja_ext * 2 + 1) * b_2);

               nzja_ext++;
               ipl++;
            } else {
               p_ja_alu[nzja_ext] = p_ja_u[ipu];
               if (_b_is_char) {
                  p_jachar_alu[nzja_ext] = p_jachar_u[ipu];
               }

               CVector < _Flt >::SetByZeroes (b_2, p_a_alu + (nzja_ext * 2) * b_2);
               CVector < _Flt >::CopyVector (b_2, p_a_u + ipu * b_2,
                                             p_a_alu + (nzja_ext * 2 + 1) * b_2);

               nzja_ext++;
               ipu++;
            }
         }
         p_ia_alu[i + 1] = nzja_ext;
      }

   }

//
// Split block pairs fct data into block L and U parts
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::SplitLUPairs (bool _b_is_char, int _b_size, int _n,
                            const vector < _Int > &_ia_alu,
                            const vector < _Int > &_ja_alu,
                            const vector < char >&_jachar_alu,
                            const vector < _Flt > &_a_alu, vector < _Int > &_ia_l,
                            vector < _Int > &_ja_l, vector < char >&_jachar_l,
                            vector < _Flt > &_a_l, vector < _Int > &_ia_u,
                            vector < _Int > &_ja_u, vector < char >&_jachar_u,
                            vector < _Flt > &_a_u)
   {

// Open structures

      int b_2 = _b_size * _b_size;

      const _Int *p_ia_alu = NULL;
      const _Int *p_ja_alu = NULL;
      const char *p_jachar_alu = NULL;
      const _Flt *p_a_alu = NULL;

      if (_ia_alu.size () > 0)
         p_ia_alu = &_ia_alu[0];
      if (_ja_alu.size () > 0)
         p_ja_alu = &_ja_alu[0];
      if (_b_is_char && _jachar_alu.size () > 0)
         p_jachar_alu = &_jachar_alu[0];
      if (_a_alu.size () > 0)
         p_a_alu = &_a_alu[0];

      int nzja_alu = (int) p_ia_alu[_n];

// Allocate and fill

      _ia_l.resize (_n + 1);
      _ja_l.resize (nzja_alu + 1);
      _a_l.resize (nzja_alu * b_2 + 1);
      _ia_u.resize (_n + 1);
      _ja_u.resize (nzja_alu + 1);
      _a_u.resize (nzja_alu * b_2 + 1);

      if (_b_is_char) {
         _jachar_l.resize (nzja_alu + 1);
         _jachar_u.resize (nzja_alu + 1);
      }

      _Int *p_ia_l = NULL;
      _Int *p_ja_l = NULL;
      char *p_jachar_l = NULL;
      _Flt *p_a_l = NULL;

      if (_ia_l.size () > 0)
         p_ia_l = &_ia_l[0];
      if (_ja_l.size () > 0)
         p_ja_l = &_ja_l[0];
      if (_b_is_char && _jachar_l.size () > 0)
         p_jachar_l = &_jachar_l[0];
      if (_a_l.size () > 0)
         p_a_l = &_a_l[0];

      _Int *p_ia_u = NULL;
      _Int *p_ja_u = NULL;
      char *p_jachar_u = NULL;
      _Flt *p_a_u = NULL;

      if (_ia_u.size () > 0)
         p_ia_u = &_ia_u[0];
      if (_ja_u.size () > 0)
         p_ja_u = &_ja_u[0];
      if (_b_is_char && _jachar_u.size () > 0)
         p_jachar_u = &_jachar_u[0];
      if (_a_u.size () > 0)
         p_a_u = &_a_u[0];

      int nzja_l = 0;
      int nzja_u = 0;

      int i, j, jj;

      p_ia_l[0] = 0;
      p_ia_u[0] = 0;

      for (i = 0; i < _n; i++) {
         for (j = (int) p_ia_alu[i]; j < p_ia_alu[i + 1]; j++) {
            jj = (int) p_ja_alu[j];
            p_ja_l[nzja_l] = jj;
            p_ja_u[nzja_u] = jj;
            if (_b_is_char) {
               p_jachar_l[nzja_l] = p_jachar_alu[j];
               p_jachar_u[nzja_u] = p_jachar_alu[j];
            }
            CVector < _Flt >::CopyVector (b_2, p_a_alu + (j * 2) * b_2,
                                          p_a_l + nzja_l * b_2);
            CVector < _Flt >::CopyVector (b_2, p_a_alu + (j * 2 + 1) * b_2,
                                          p_a_u + nzja_u * b_2);
            nzja_l++;
            nzja_u++;
         }
         p_ia_l[i + 1] = nzja_l;
         p_ia_u[i + 1] = nzja_u;
      }

   }

//
// Rescale block factor back
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::RescaleU (int _b_size, int _n, const vector < _Flt > &_sclU,
                        const vector < _Flt > &_invsclU, const vector < _Int > &_ia_u,
                        const vector < _Int > &_ja_u, vector < _Flt > &_a_u)
   {

// Open structures

      int b_2 = _b_size * _b_size;

      const _Flt *p_sclU = NULL;
      const _Flt *p_invsclU = NULL;

      if (_sclU.size () > 0)
         p_sclU = &_sclU[0];
      if (_invsclU.size () > 0)
         p_invsclU = &_invsclU[0];

      const _Int *p_ia_u = NULL;
      const _Int *p_ja_u = NULL;
      _Flt *p_a_u = NULL;

      if (_ia_u.size () > 0)
         p_ia_u = &_ia_u[0];
      if (_ja_u.size () > 0)
         p_ja_u = &_ja_u[0];
      if (_a_u.size () > 0)
         p_a_u = &_a_u[0];

// Rescaling cycle

      CVectorData < _Flt > xtemp (b_2);
      _Flt *pxtemp = xtemp.Ptr ();

      int i, j, jj;

      for (i = 0; i < _n; i++) {
         for (j = p_ia_u[i]; j < p_ia_u[i + 1]; j++) {
            jj = p_ja_u[j];
            if (i != jj) {

               CBlock_BxB_traits < _Flt >::MM_BxB (_b_size, p_a_u + j * b_2,
                                                   p_invsclU + jj * b_2, pxtemp);

               CVector < _Flt >::CopyVector (b_2, pxtemp, p_a_u + j * b_2);

            } else {

               CBlock_BxB_traits < _Flt >::MM_BxB (_b_size, p_sclU + i * b_2,
                                                   p_a_u + j * b_2, pxtemp);

               CVector < _Flt >::CopyVector (b_2, pxtemp, p_a_u + j * b_2);

            }
         }
      }

   }

//
// Rescale block factor back
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::RescaleU (char _type, int _b_size, int _nlist, const _Flt * _sclU,
                        const _Flt * _invsclU, const vector < _Int > &_list_u,
                        const vector < _Int > &_ia_u, const vector < _Int > &_ja_u,
                        vector < _Flt > &_a_u)
   {

// Open structures

      int b_2 = _b_size * _b_size;

      const _Flt *p_sclU = _sclU;
      const _Flt *p_invsclU = _invsclU;

      const _Int *p_list_u = NULL;
      const _Int *p_ia_u = NULL;
      const _Int *p_ja_u = NULL;
      _Flt *p_a_u = NULL;

      if (_list_u.size () > 0)
         p_list_u = &_list_u[0];
      if (_ia_u.size () > 0)
         p_ia_u = &_ia_u[0];
      if (_ja_u.size () > 0)
         p_ja_u = &_ja_u[0];
      if (_a_u.size () > 0)
         p_a_u = &_a_u[0];

// Rescaling cycle

      CVectorData < _Flt > xtemp (b_2);
      _Flt *pxtemp = xtemp.Ptr ();

      int ilist, i, j, jj;

      for (ilist = 0; ilist < _nlist; ilist++) {
         i = (int) p_list_u[ilist];
         for (j = (int) p_ia_u[ilist]; j < p_ia_u[ilist + 1]; j++) {
            jj = (int) p_ja_u[j];
            if (i != jj) {

               CBlock_BxB_traits < _Flt >::MM_BxB (_b_size, p_a_u + j * b_2,
                                                   p_invsclU + jj * b_2, pxtemp);

               CVector < _Flt >::CopyVector (b_2, pxtemp, p_a_u + j * b_2);

            } else {

               if (_type == 'D') {

                  CBlock_BxB_traits < _Flt >::MM_BxB (_b_size, p_sclU + i * b_2,
                                                      p_a_u + j * b_2, pxtemp);

                  CVector < _Flt >::CopyVector (b_2, pxtemp, p_a_u + j * b_2);

               } else {

                  CBlock_BxB_traits < _Flt >::MM_BxB (_b_size, p_a_u + j * b_2,
                                                      p_invsclU + jj * b_2, pxtemp);

                  CVector < _Flt >::CopyVector (b_2, pxtemp, p_a_u + j * b_2);

               }

            }
         }
      }

   }

//
// Perform ICH2 small dense blocks sparse factorization of the block with future block diagonal modifications (with matrix degree structural control)
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::FctBlockIch2Degree (SParams * _params, int _b_size, int _n, int _n_ini,
                                  const vector < _Int > &_ia_alu,
                                  const vector < _Int > &_ja_alu,
                                  const vector < char >&_ja_char_alu,
                                  const vector < _Flt > &_a_alu, vector < _Int > &_ia_lu,
                                  vector < _Int > &_ja_lu, vector < char >&_ja_char_lu,
                                  vector < _Flt > &_a_lu, int &_nmodif,
                                  double &_eigmin_att, double &_eigmax_att)
   {

//      ofstream ffout ("ChkFctDeg.dat");

//      int irow_chk = 5;

//      int index1 = 0;
//      int index2 = 189;
//      int index3 = 202;

// Open matrix data structures

      int b_2 = _b_size * _b_size;

      const _Int *p_ia_alu = _ia_alu.data ();
      const _Int *p_ja_alu = _ja_alu.data ();
      const char *p_ja_char_alu = _ja_char_alu.data ();
      const _Flt *p_a_alu = _a_alu.data ();

// Prepare mask arrays

      vector < _Int > ibegm (_n + 1);
      vector < _Int > madj (_n + 1);
      vector < _Int > iv (_n + 1);
      vector < _Int > imaskc (_n + 1);
      vector < char >imaskchar (_n + 1);
      vector < _Int > listc (_n + 1);

      vector < _Flt > fmaskc (b_2 * _n);
      vector < _Flt > a_dia (b_2 * _n);

      _Flt *pfmaskc = fmaskc.data ();
      _Flt *pa_dia = a_dia.data ();

      int i;

      for (i = 0; i <= _n; i++)
         ibegm[i] = -1;
      for (i = 0; i <= _n; i++)
         madj[i] = -1;
      for (i = 0; i <= _n; i++)
         imaskc[i] = -1;

      CVector < _Flt >::SetByZeroes (b_2 * _n, pa_dia);

      double eigmin_att = FLT_MAX;
      double eigmax_att = -FLT_MAX;

      int j, jj, k, kcolmn, irwprv, nlistcnew, ind, ind1, irow, irwpr1, jcolmn;

      _Flt *pauxL;
      _Flt *pauxU_row;

      CVectorData < _Flt > u_diag (b_2);

      _Flt *pu_diag = u_diag.Ptr ();

      CVectorData < _Flt > u_offd (b_2);

      _Flt *pu_offd = u_offd.Ptr ();

      CVectorData < _Flt > Work (10 * b_2);
      CVectorData < double >dWork (20 * b_2);

      _Flt *pWork = Work.Ptr ();
      double *pdWork = dWork.Ptr ();

      int fcttype_loc = _params->fcttype;
      int fcttype_dia_loc = _params->fcttype_dia;
      double pivmin_loc = _params->pivmin;
      double tau1_loc = _params->tau1;
      double tau2_loc = _params->tau2;
      double theta_loc = _params->theta;

//      ofstream *pfout_debug = _params->pfout;
//      *pfout_debug << " N = " << _n << " N_ini = " << _n_ini << endl;

      int icycle_int = -1;

      _ia_lu.resize (_n + 1);

      _ja_lu.resize (0);
      _ja_char_lu.resize (0);
      _a_lu.resize (0);

      _ja_lu.reserve (_n + 1);
      _ja_char_lu.reserve (_n + 1);
      _a_lu.reserve (_n * b_2 + 1);

      _ia_lu[0] = 0;

      int nzja_lu = 0;

      _nmodif = 0;

      char jjchar1, jjchar2, jjchar;
      double dfnorm;

      double eigmin_temp, eigmax_temp;

      for (i = 0; i < _n; i++) {

         irow = i;

//         if (irow == index1 || irow == index2 || irow == index3) {
//            if (pfout_debug != NULL) {
//               *pfout_debug << " <<< Irow = " << i << endl;
//            }
//         }

// Init current row

         int nlistcloc = 0;

         icycle_int++;

         if (i < _n_ini) {
            for (k = (int) p_ia_alu[i]; k < p_ia_alu[i + 1]; k++) {
               kcolmn = (int) p_ja_alu[k];
               listc[nlistcloc] = kcolmn;
               nlistcloc++;
               CVector < _Flt >::CopyVector (b_2, p_a_alu + k * b_2,
                                             pfmaskc + kcolmn * b_2);
               imaskc[kcolmn] = icycle_int;
               imaskchar[kcolmn] = p_ja_char_alu[k];
            }
            if (imaskc[i] != icycle_int) {
               listc[nlistcloc] = i;
               nlistcloc++;
               CVector < _Flt >::SetByZeroes (b_2, pfmaskc + i * b_2);
               imaskc[i] = icycle_int;
               imaskchar[i] = 0;
            }
         } else {
            listc[nlistcloc] = i;
            nlistcloc++;
            CVector < _Flt >::SetByZeroes (b_2, pfmaskc + i * b_2);
            imaskc[i] = icycle_int;
            imaskchar[i] = 0;
         }

//         if (irow == index1 || irow == index2 || irow == index3) {
//            if (pfout_debug != NULL) {
//               if (i < _n_ini) {
//                  int niloc = (int)(p_ia_alu[i + 1]-p_ia_alu[i]);
//                  int ibeg = (int)p_ia_alu[i];
//                  PrintArray (*pfout_debug, " JaRow_ini ", niloc, p_ja_alu+ibeg);
//                  PrintArray (*pfout_debug, " JaCharRow_ini ", niloc, p_ja_char_alu+ibeg);
//                  PrintArray (*pfout_debug, " ARow_ini ", niloc*b_2_2, p_a_alu+ibeg*b_2_2);
//               }
//            }
//         }

// Update current row

         irwprv = (int) ibegm[irow];

         while (irwprv != -1) {

//            if (i >= _n_ini && i < _n_ini+irow_chk) {
//               if (pfout_debug != NULL) {
//                  *pfout_debug << "       ++ Update by Irow prev = " << irwprv << endl;
//               }
//            }

            irwpr1 = (int) madj[irwprv];

            if (iv[irwprv] < _ia_lu[irwprv + 1]) {

               j = (int) iv[irwprv];

               jj = (int) _ja_lu[j];
               jjchar1 = _ja_char_lu[j] + 1;

               pauxL = _a_lu.data () + j * b_2;

               jcolmn = jj;
               if (jcolmn < 0)
                  jcolmn = -jcolmn - 1;

//               if (i >= _n_ini && i < _n_ini+irow_chk) {
//                  if (pfout_debug != NULL) {
//                     *pfout_debug << "       ++ Base Values jj = " << jj << " jjchar1 = " << (int)jjchar1 << " auxL = " << *pauxL << " auxU = " << *pauxU << endl;
//                  }
//               }

               if (jj >= 0) {

                  for (k = (int) iv[irwprv]; k < _ia_lu[irwprv + 1]; k++) {
                     kcolmn = (int) _ja_lu[k];
                     if (kcolmn < 0)
                        kcolmn = -kcolmn - 1;
                     jjchar2 = _ja_char_lu[k];
                     jjchar = jjchar1 + jjchar2;
//                     if (i >= _n_ini && i < _n_ini+irow_chk) {
//                        if (pfout_debug != NULL) {
//                           *pfout_debug << "          -- Row update kcolmn = " << kcolmn << " jjchar = " << (int)jjchar << endl;
//                        }
//                     }
                     if (kcolmn == irow || fcttype_loc < 0
                         || (fcttype_loc >= 0 && jjchar <= fcttype_loc)) {

                        if (imaskc[kcolmn] != icycle_int) {
                           listc[nlistcloc] = kcolmn;
                           nlistcloc++;
                           CVector < _Flt >::SetByZeroes (b_2, pfmaskc + kcolmn * b_2);
                           imaskc[kcolmn] = icycle_int;
                           imaskchar[kcolmn] = jjchar;
                        } else {
                           if (imaskchar[kcolmn] > jjchar)
                              imaskchar[kcolmn] = jjchar;
                        }

                        pauxU_row = _a_lu.data () + k * b_2;

                        CBlock_BxB_traits < _Flt >::MtM_Sub_BxB (_b_size, pauxL,
                                                                 pauxU_row,
                                                                 pfmaskc + kcolmn * b_2);

//                        if (i >= _n_ini && i < _n_ini+irow_chk) {
//                           if (pfout_debug != NULL) {
//                              *pfout_debug << "          == New values: L " << pfmaskc[kcolmn*b_2_2] << " U = " << pfmaskc[kcolmn * b_2_2 + b_2] << endl;
//                           }
//                        }

                     } else {

                        pauxU_row = _a_lu.data () + k * b_2;

                        CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, pauxL, pauxU_row,
                                                             pu_offd);

                        CBlock_BxB_traits < _Flt >::ModifDia_BxB (_b_size, theta_loc,
                                                                  pu_offd, pu_offd,
                                                                  pa_dia + irow * b_2,
                                                                  pa_dia + kcolmn * b_2);

//                        if (i >= _n_ini && i < _n_ini+irow_chk) {
//                           if (pfout_debug != NULL) {
//                              *pfout_debug << "          :: Modif: irow = " << pa_dia[irow*b_2]  << " kcolmn = " << pa_dia[kcolmn*b_2] << endl;
//                           }
//                        }

                     }
                  }
               } else {
                  for (k = (int) iv[irwprv]; k < _ia_lu[irwprv + 1]; k++) {
                     kcolmn = (int) _ja_lu[k];
                     if (kcolmn >= 0) {
                        jjchar2 = _ja_char_lu[k];
                        jjchar = jjchar1 + jjchar2;
//                        if (i >= _n_ini && i < _n_ini+irow_chk) {
//                           if (pfout_debug != NULL) {
//                              *pfout_debug << "          -- Row 2 update kcolmn = " << kcolmn << " jjchar = " << (int)jjchar << endl;
//                           }
//                        }
                        if (kcolmn == irow || fcttype_loc < 0
                            || (fcttype_loc >= 0 && jjchar <= fcttype_loc)) {
                           if (imaskc[kcolmn] != icycle_int) {
                              listc[nlistcloc] = kcolmn;
                              nlistcloc++;
                              CVector < _Flt >::SetByZeroes (b_2, pfmaskc + kcolmn * b_2);
                              imaskc[kcolmn] = icycle_int;
                              imaskchar[kcolmn] = jjchar;
                           } else {
                              if (imaskchar[kcolmn] > jjchar)
                                 imaskchar[kcolmn] = jjchar;
                           }

                           pauxU_row = _a_lu.data () + k * b_2;

                           CBlock_BxB_traits < _Flt >::MtM_Sub_BxB (_b_size, pauxL,
                                                                    pauxU_row,
                                                                    pfmaskc +
                                                                    kcolmn * b_2);

//                           if (i >= _n_ini && i < _n_ini+irow_chk) {
//                              if (pfout_debug != NULL) {
//                                 *pfout_debug << "          == New values 2: L " << pfmaskc[kcolmn * b_2_2] << " U = " << pfmaskc[kcolmn * b_2_2 + b_2] << endl;
//                              }
//                           }

                        } else {

                           pauxU_row = _a_lu.data () + k * b_2;

                           CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, pauxL, pauxU_row,
                                                                pu_offd);

                           CBlock_BxB_traits < _Flt >::ModifDia_BxB (_b_size, theta_loc,
                                                                     pu_offd, pu_offd,
                                                                     pa_dia + irow * b_2,
                                                                     pa_dia +
                                                                     kcolmn * b_2);

//                           if (i >= _n_ini && i < _n_ini+irow_chk) {
//                              if (pfout_debug != NULL) {
//                                 *pfout_debug << "          :: Modif 2: irow = " << pa_dia[irow*b_2]  << " kcolmn = " << pa_dia[kcolmn*b_2] << endl;
//                              }
//                           }

                        }
                     }
                  }
               }

            }

            irwprv = irwpr1;

         }

// Update transposed structure

         irwprv = (int) ibegm[irow];

         while (irwprv != -1) {

            irwpr1 = (int) madj[irwprv];

            if (iv[irwprv] >= _ia_lu[irwprv + 1] - 1) {
               madj[irwprv] = -1;
            } else {
               j = (int) iv[irwprv] + 1;
               jcolmn = (int) _ja_lu[j];
               if (jcolmn < 0)
                  jcolmn = -jcolmn - 1;
               madj[irwprv] = ibegm[jcolmn];
               ibegm[jcolmn] = irwprv;
            }

            iv[irwprv]++;

            irwprv = irwpr1;

         }

// Perform filtering of the data and modify diagonal

         if (i >= _n_ini) {
            fcttype_loc = _params->fcttype;
            tau2_loc = 0.0e0;
         }

         {

            nlistcnew = 0;

            for (j = 0; j < nlistcloc; j++) {
               ind = (int) listc[j];
               if (ind != irow) {

                  jjchar = imaskchar[ind];

                  dfnorm =
                     CBlock_BxB_traits < _Flt >::FNormValue_BxB (_b_size,
                                                                 pfmaskc + ind * b_2,
                                                                 pfmaskc + ind * b_2 +
                                                                 b_2);

                  if (jjchar > 0 && ((fcttype_loc < 0 && dfnorm < tau2_loc)
                                     || (fcttype_loc >= 0
                                         && (jjchar > fcttype_loc || dfnorm < tau2_loc)))
                     ) {

                     CBlock_BxB_traits < _Flt >::ModifDia_BxB (_b_size, theta_loc,
                                                               pfmaskc + ind * b_2,
                                                               pfmaskc + ind * b_2,
                                                               pa_dia + irow * b_2,
                                                               pa_dia + ind * b_2);

                  } else {
                     listc[nlistcnew] = ind;
                     nlistcnew++;
                  }
               } else {
                  listc[nlistcnew] = irow;
                  nlistcnew++;
               }
            }

            nlistcloc = nlistcnew;

         }

//         if (irow == index1 || irow == index2 || irow == index3) {
//            if (pfout_debug != NULL) {
//               *pfout_debug << "     Dia before modif add: L " << pfmaskc[irow * b_2_2] << " U = " << pfmaskc[irow * b_2_2 + b_2] << endl;
//            }
//         }

         CBlock_BxB_traits < _Flt >::AddModifDia_BxB (_b_size, pa_dia + irow * b_2,
                                                      pfmaskc + irow * b_2);

//         if (irow == index1 || irow == index2 || irow == index3) {
//            if (pfout_debug != NULL) {
//               *pfout_debug << "     Dia after modif add: L " << pfmaskc[irow * b_2_2] << " U = " << pfmaskc[irow * b_2_2 + b_2] << endl;
//            }
//         }

// Factorize current row and split into first and second order

         _Int *plistc = &listc[0];

         sort (plistc, plistc + nlistcloc);

//         if (irow == index1 || irow == index2 || irow == index3) {
//            if (pfout_debug != NULL) {
//               *pfout_debug << " Updated row before fct: " << endl;
//               for (j = 0; j < nlistcloc; j++) {
//                  ind = (int) listc[j];
//                  *pfout_debug << " Jcol = " << ind << endl;
//                  *pfout_debug << " Jchar = " << (int)imaskchar[ind] << endl;
//                  *pfout_debug << " Jcol = " << ind << endl;
//                  PrintArray (*pfout_debug, " Aelem ", b_2_2, pfmaskc+ind*b_2_2);
//               }
//            }
//         }

         if (i < _n_ini) {

//            if (fcttype_dia_loc == 0) {

            _nmodif +=
               CBlock_BxB_traits < _Flt >::FctSymmInv_BxB (pivmin_loc, _b_size,
                                                           pfmaskc + irow * b_2, pu_offd,
                                                           pu_diag, eigmin_temp,
                                                           eigmax_temp);

            CVector < _Flt >::CopyVector (b_2, pu_diag, pfmaskc + irow * b_2);

            if (eigmin_temp < eigmin_att)
               eigmin_att = eigmin_temp;
            if (eigmax_temp > eigmax_att)
               eigmax_att = eigmax_temp;

//            }
         }

         if (i < _n_ini) {

            for (j = 0; j < nlistcloc; j++) {
               ind = (int) listc[j];
               if (ind != irow) {

                  CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, pu_diag,
                                                       pfmaskc + ind * b_2, pu_offd);

                  CVector < _Flt >::CopyVector (b_2, pu_offd, pfmaskc + ind * b_2);

                  dfnorm =
                     CBlock_BxB_traits < _Flt >::FNormValue_BxB (_b_size,
                                                                 pfmaskc + ind * b_2,
                                                                 pfmaskc + ind * b_2);
                  jjchar = imaskchar[ind];

                  if (jjchar > 0 && dfnorm < tau1_loc) {
                     listc[j] = (_Int) (-ind - 1);
                  }
               }
            }

         }
// Store computed row elems

         _a_lu.resize ((nzja_lu + nlistcloc) * b_2);

         _Flt *p_a_lu_temp = _a_lu.data ();

         for (j = 0; j < nlistcloc; j++) {
            ind = (int) listc[j];
            _ja_lu.push_back (ind);
            ind1 = ind;
            if (ind1 < 0)
               ind1 = -ind1 - 1;
            jjchar = imaskchar[ind1];
            _ja_char_lu.push_back (jjchar);
            CVector < _Flt >::CopyVector (b_2, pfmaskc + ind1 * b_2,
                                          p_a_lu_temp + (nzja_lu + j) * b_2);
         }

//         if (irow == index1 || irow == index2 || irow == index3) {
//            if (pfout_debug != NULL) {
//               _Int *p_ja_lu_temp = &_ja_lu[0];
//               char *p_jachar_lu_temp = &_ja_char_lu[0];
//               _Flt *p_a_lu_temp = &_a_lu[0];
//               PrintArray (*pfout_debug, " Ja row ",nlistcloc,p_ja_lu_temp+nzja_lu);
//               PrintArray (*pfout_debug, " JaChar row ",nlistcloc,p_jachar_lu_temp+nzja_lu);
//               PrintArray (*pfout_debug, " A row ",nlistcloc*b_2_2,p_a_lu_temp+nzja_lu*b_2_2);
//            }
//         }

         nzja_lu += nlistcloc;

         _ia_lu[i + 1] = nzja_lu;

// Add current row into the transposed structures

         if (nlistcloc > 1) {
            if (i < _n_ini) {
               iv[irow] = _ia_lu[i] + 1;
               ind = (int) listc[1];
               if (ind < 0)
                  ind = -ind - 1;
               madj[irow] = ibegm[ind];
               ibegm[ind] = (_Int) irow;
            } else {
               iv[irow] = _ia_lu[i + 1] + 1;
            }
         }

      }

      _eigmin_att = eigmin_att;
      _eigmax_att = eigmax_att;

// Condense the result (filter second order elems)

      int nzja_new = 0;

      for (i = 0; i < nzja_lu; i++) {
         if (_ja_lu[i] >= 0)
            nzja_new++;
      }

      vector < _Int > jlu_cnd (nzja_new + 1);
      vector < char >jlu_char_cnd (nzja_new + 1);
      vector < _Flt > lu_cnd (b_2 * nzja_new + 1);

      _Int *pjlu_cnd = &jlu_cnd[0];
      char *pjlu_char_cnd = &jlu_char_cnd[0];
      _Flt *plu_cnd = &lu_cnd[0];

      nzja_new = 0;

      listc[0] = 0;

      _Flt *p_a_lu_temp = _a_lu.data ();

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_lu[i]; j < _ia_lu[i + 1]; j++) {
            if (_ja_lu[j] >= 0) {
               pjlu_cnd[nzja_new] = _ja_lu[j];
               pjlu_char_cnd[nzja_new] = _ja_char_lu[j];
               CVector < _Flt >::CopyVector (b_2, p_a_lu_temp + j * b_2,
                                             plu_cnd + nzja_new * b_2);
               nzja_new++;
            }
         }
         listc[i + 1] = nzja_new;
      }

      for (i = 0; i <= _n; i++)
         _ia_lu[i] = listc[i];

      _ja_lu.swap (jlu_cnd);
      _ja_char_lu.swap (jlu_char_cnd);
      _a_lu.swap (lu_cnd);

   }

//
// Perform ILU2 small dense blocks sparse factorization of the block with future block diagonal modifications (with matrix degree structural control)
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::FctBlockIlu2Degree (SParams * _params, int _b_size, int _n, int _n_ini,
                                  const vector < _Int > &_ia_alu,
                                  const vector < _Int > &_ja_alu,
                                  const vector < char >&_ja_char_alu,
                                  const vector < _Flt > &_a_alu, vector < _Int > &_ia_lu,
                                  vector < _Int > &_ja_lu, vector < char >&_ja_char_lu,
                                  vector < _Flt > &_a_lu, int &_nmodif,
                                  double &_eigmin_att, double &_eigmax_att)
   {

      ofstream ffout_debug;

      if (_params->b_write_file) {
         ffout_debug.open (_params->name_file);
      }

      int index1 = _params->iparam1;
      int index2 = _params->iparam2;
      int index3 = _params->iparam3;

// Open matrix data structures

      int b_2 = _b_size * _b_size;
      int b_2_2 = 2 * b_2;

      const _Int *p_ia_alu = NULL;
      const _Int *p_ja_alu = NULL;
      const char *p_ja_char_alu = NULL;
      const _Flt *p_a_alu = NULL;

      if (_ia_alu.size () > 0)
         p_ia_alu = &_ia_alu[0];
      if (_ja_alu.size () > 0)
         p_ja_alu = &_ja_alu[0];
      if (_ja_char_alu.size () > 0)
         p_ja_char_alu = &_ja_char_alu[0];
      if (_a_alu.size () > 0)
         p_a_alu = &_a_alu[0];

// Prepare mask arrays

      vector < _Int > ibegm (_n + 1);
      vector < _Int > madj (_n + 1);
      vector < _Int > iv (_n + 1);
      vector < _Int > imaskc (_n + 1);
      vector < char >imaskchar (_n + 1);
      vector < _Int > listc (_n + 1);

      vector < _Flt > fmaskc (b_2_2 * _n + 1);
      vector < _Flt > a_dia (b_2 * _n + 1);

      _Flt *pfmaskc = &fmaskc[0];
      _Flt *pa_dia = &a_dia[0];

      int i;

      for (i = 0; i <= _n; i++)
         ibegm[i] = -1;
      for (i = 0; i <= _n; i++)
         madj[i] = -1;
      for (i = 0; i <= _n; i++)
         imaskc[i] = -1;

      CVector < _Flt >::SetByZeroes (b_2 * _n, pa_dia);

      double eigmin_att = FLT_MAX;
      double eigmax_att = -FLT_MAX;

      int j, jj, k, kcolmn, irwprv, nlistcnew, ind, ind1, irow, irwpr1, jcolmn;

      _Flt *pauxL;
      _Flt *pauxU;

      _Flt *pauxL_col;
      _Flt *pauxU_row;

      CVectorData < _Flt > lu_diag (b_2_2);

      _Flt *pl_diag = lu_diag.Ptr ();
      _Flt *pu_diag = pl_diag + b_2;

      CVectorData < _Flt > lu_offd (b_2_2);

      _Flt *pl_offd = lu_offd.Ptr ();
      _Flt *pu_offd = pl_offd + b_2;

      CVectorData < _Flt > Work (20 * b_2);
      CVectorData < double >dWork (20 * b_2);

      _Flt *pWork = Work.Ptr ();
      double *pdWork = dWork.Ptr ();

      int fcttype_loc = _params->fcttype;
      int fcttype_dia_loc = _params->fcttype_dia;
      double pivmin_loc = _params->pivmin;
      double tau1_loc = _params->tau1;
      double tau2_loc = _params->tau2;
      double theta_loc = _params->theta;

      if (_params->b_write_file) {
         ffout_debug << " N = " << _n << " N_ini = " << _n_ini << endl;
//         ffout_debug << " Matrix on entry: " << endl;
//         for (int i=0;i<_n_ini;i++) {
//            ffout_debug << " Irow = " << i << endl;
//            PrintArray (ffout_debug, " Ja ",(int)(_ia_alu[i+1]-_ia_alu[i]),_ja_alu.data()+_ia_alu[i]);
//            PrintArray (ffout_debug, " A ",(int)((_ia_alu[i+1]-_ia_alu[i])*b_2_2),_a_alu.data()+_ia_alu[i]*b_2_2);
//         }
      }

      int icycle_int = -1;

      _ia_lu.resize (_n + 1);

      _ja_lu.resize (0);
      _ja_char_lu.resize (0);
      _a_lu.resize (0);

      _ja_lu.reserve (_n + 1);
      _ja_char_lu.reserve (_n + 1);
      _a_lu.reserve (_n * b_2_2 + 1);

      _ia_lu[0] = 0;

      int nzja_lu = 0;

      _nmodif = 0;

      char jjchar1, jjchar2, jjchar;
      double dfnorm;

      double eigmin_temp, eigmax_temp;

      for (i = 0; i < _n; i++) {

         irow = i;

         if (_params->b_write_file) {
            if (irow == index1 || irow == index2 || irow == index3) {
               ffout_debug << " <<< Irow = " << i << endl;
            }
         }
// Init current row

         int nlistcloc = 0;

         icycle_int++;

         if (i < _n_ini) {
            for (k = (int) p_ia_alu[i]; k < p_ia_alu[i + 1]; k++) {
               kcolmn = (int) p_ja_alu[k];
               listc[nlistcloc] = kcolmn;
               nlistcloc++;
               CVector < _Flt >::CopyVector (b_2_2, p_a_alu + k * b_2_2,
                                             pfmaskc + kcolmn * b_2_2);
               imaskc[kcolmn] = icycle_int;
               imaskchar[kcolmn] = p_ja_char_alu[k];
            }
            if (imaskc[i] != icycle_int) {
               listc[nlistcloc] = i;
               nlistcloc++;
               CVector < _Flt >::SetByZeroes (b_2_2, pfmaskc + i * b_2_2);
               imaskc[i] = icycle_int;
               imaskchar[i] = 0;
            }
         } else {
            listc[nlistcloc] = i;
            nlistcloc++;
            CVector < _Flt >::SetByZeroes (b_2_2, pfmaskc + i * b_2_2);
            imaskc[i] = icycle_int;
            imaskchar[i] = 0;
         }

         if (_params->b_write_file) {
            if (irow == index1 || irow == index2 || irow == index3) {
               if (i < _n_ini) {
                  int niloc = (int) (p_ia_alu[i + 1] - p_ia_alu[i]);
                  int ibeg = (int) p_ia_alu[i];
                  PrintArray (ffout_debug, " JaRow_ini ", niloc, p_ja_alu + ibeg);
                  PrintArray (ffout_debug, " JaCharRow_ini ", niloc,
                              p_ja_char_alu + ibeg);
                  PrintArray (ffout_debug, " ARow_ini ", niloc * b_2_2,
                              p_a_alu + ibeg * b_2_2);
               }
            }
         }
// Update current row

         irwprv = (int) ibegm[irow];

         while (irwprv != -1) {

            if (_params->b_write_file) {
               if (irow == index1 || irow == index2 || irow == index3) {
                  ffout_debug << "       ++ Update by Irow prev = " << irwprv << endl;
               }
            }

            irwpr1 = (int) madj[irwprv];

            if (iv[irwprv] < _ia_lu[irwprv + 1]) {

               j = (int) iv[irwprv];

               jj = (int) _ja_lu[j];
               jjchar1 = _ja_char_lu[j] + 1;

               pauxL = &_a_lu[0] + j * b_2_2;
               pauxU = pauxL + b_2;

               jcolmn = jj;
               if (jcolmn < 0)
                  jcolmn = -jcolmn - 1;

               if (_params->b_write_file) {
                  if (irow == index1 || irow == index2 || irow == index3) {
                     if (jcolmn == irow) {
                        ffout_debug << "       ++ Base Values jj = " << jj <<
                           " jjchar1 = " << (int) jjchar1 << endl;
                        PrintArray (ffout_debug, " OffL = ", b_2, pauxL);
                        PrintArray (ffout_debug, " OffU = ", b_2, pauxU);
                     }
                  }
               }

               if (jj >= 0) {

                  for (k = (int) iv[irwprv]; k < _ia_lu[irwprv + 1]; k++) {
                     kcolmn = (int) _ja_lu[k];
                     if (kcolmn < 0)
                        kcolmn = -kcolmn - 1;
                     jjchar2 = _ja_char_lu[k];
                     jjchar = jjchar1 + jjchar2;
                     if (_params->b_write_file) {
                        if (irow == index1 || irow == index2 || irow == index3) {
                           if (kcolmn == irow) {
                              ffout_debug << "          -- Row update kcolmn = " << kcolmn
                                 << " jjchar = " << (int) jjchar << endl;
                           }
                        }
                     }
                     if (kcolmn == irow || fcttype_loc < 0
                         || (fcttype_loc >= 0 && jjchar <= fcttype_loc)) {

                        if (imaskc[kcolmn] != icycle_int) {
                           listc[nlistcloc] = kcolmn;
                           nlistcloc++;
                           CVector < _Flt >::SetByZeroes (b_2_2,
                                                          pfmaskc + kcolmn * b_2_2);
                           imaskc[kcolmn] = icycle_int;
                           imaskchar[kcolmn] = jjchar;
                        } else {
                           if (imaskchar[kcolmn] > jjchar)
                              imaskchar[kcolmn] = jjchar;
                        }

                        pauxL_col = &_a_lu[0] + k * b_2_2;
                        pauxU_row = pauxL_col + b_2;

                        CBlock_BxB_traits < _Flt >::MtM_Sub_BxB (_b_size, pauxU,
                                                                 pauxL_col,
                                                                 pfmaskc +
                                                                 kcolmn * b_2_2);
                        CBlock_BxB_traits < _Flt >::MtM_Sub_BxB (_b_size, pauxL,
                                                                 pauxU_row,
                                                                 pfmaskc +
                                                                 kcolmn * b_2_2 + b_2);

                        if (_params->b_write_file) {
                           if (irow == index1 || irow == index2 || irow == index3) {
                              if (kcolmn == irow) {
                                 PrintArray (ffout_debug, " Mult L: U = ", b_2, pauxU);
                                 PrintArray (ffout_debug, "       : L col = ", b_2,
                                             pauxL_col);
                                 PrintArray (ffout_debug, " Mult U: L = ", b_2, pauxL);
                                 PrintArray (ffout_debug, "       : U row = ", b_2,
                                             pauxU_row);
                                 ffout_debug << "          == New values: " << endl;
                                 PrintArray (ffout_debug, " L = ", b_2,
                                             pfmaskc + kcolmn * b_2_2);
                                 PrintArray (ffout_debug, " U = ", b_2,
                                             pfmaskc + kcolmn * b_2_2 + b_2);
                              }
                           }
                        }

                     } else {

                        pauxL_col = &_a_lu[0] + k * b_2_2;
                        pauxU_row = pauxL_col + b_2;

                        CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, pauxU, pauxL_col,
                                                             pl_offd);
                        CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, pauxL, pauxU_row,
                                                             pu_offd);

                        CBlock_BxB_traits < _Flt >::ModifDia_BxB (_b_size, theta_loc,
                                                                  pl_offd, pu_offd,
                                                                  pa_dia + irow * b_2,
                                                                  pa_dia + kcolmn * b_2);

                        if (_params->b_write_file) {
                           if (irow == index1 || irow == index2 || irow == index3) {
                              if (kcolmn == irow) {
                                 PrintArray (ffout_debug, " Mult L: U = ", b_2, pauxU);
                                 PrintArray (ffout_debug, "       : L col = ", b_2,
                                             pauxL_col);
                                 PrintArray (ffout_debug, " Mult U: L = ", b_2, pauxL);
                                 PrintArray (ffout_debug, "       : U row = ", b_2,
                                             pauxU_row);
                                 ffout_debug << "          :: Modif: " << endl;
                                 PrintArray (ffout_debug, " Irow   = ", b_2,
                                             pa_dia + irow * b_2);
                                 PrintArray (ffout_debug, " Kcolmn = ", b_2,
                                             pa_dia + kcolmn * b_2);
                              }
                           }
                        }

                     }
                  }
               } else {
                  for (k = (int) iv[irwprv]; k < _ia_lu[irwprv + 1]; k++) {
                     kcolmn = (int) _ja_lu[k];
                     if (kcolmn >= 0) {
                        jjchar2 = _ja_char_lu[k];
                        jjchar = jjchar1 + jjchar2;
                        if (_params->b_write_file) {
                           if (irow == index1 || irow == index2 || irow == index3) {
                              if (kcolmn == irow) {
                                 ffout_debug << "          -- Row 2 update kcolmn = " <<
                                    kcolmn << " jjchar = " << (int) jjchar << endl;
                              }
                           }
                        }
                        if (kcolmn == irow || fcttype_loc < 0
                            || (fcttype_loc >= 0 && jjchar <= fcttype_loc)) {
                           if (imaskc[kcolmn] != icycle_int) {
                              listc[nlistcloc] = kcolmn;
                              nlistcloc++;
                              CVector < _Flt >::SetByZeroes (b_2_2,
                                                             pfmaskc + kcolmn * b_2_2);
                              imaskc[kcolmn] = icycle_int;
                              imaskchar[kcolmn] = jjchar;
                           } else {
                              if (imaskchar[kcolmn] > jjchar)
                                 imaskchar[kcolmn] = jjchar;
                           }

                           pauxL_col = &_a_lu[0] + k * b_2_2;
                           pauxU_row = pauxL_col + b_2;

                           CBlock_BxB_traits < _Flt >::MtM_Sub_BxB (_b_size, pauxU,
                                                                    pauxL_col,
                                                                    pfmaskc +
                                                                    kcolmn * b_2_2);
                           CBlock_BxB_traits < _Flt >::MtM_Sub_BxB (_b_size, pauxL,
                                                                    pauxU_row,
                                                                    pfmaskc +
                                                                    kcolmn * b_2_2 + b_2);

                           if (_params->b_write_file) {
                              if (irow == index1 || irow == index2 || irow == index3) {
                                 if (kcolmn == irow) {
                                    PrintArray (ffout_debug, " Mult L 2: U = ", b_2,
                                                pauxU);
                                    PrintArray (ffout_debug, "        2: L col = ", b_2,
                                                pauxL_col);
                                    PrintArray (ffout_debug, " Mult U 2: L = ", b_2,
                                                pauxL);
                                    PrintArray (ffout_debug, "        2: U row = ", b_2,
                                                pauxU_row);
                                    ffout_debug << "          == New values 2: " << endl;
                                    PrintArray (ffout_debug, " L = ", b_2,
                                                pfmaskc + kcolmn * b_2_2);
                                    PrintArray (ffout_debug, " U = ", b_2,
                                                pfmaskc + kcolmn * b_2_2 + b_2);
                                 }
                              }
                           }

                        } else {

                           pauxL_col = &_a_lu[0] + k * b_2_2;
                           pauxU_row = pauxL_col + b_2;

                           CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, pauxU, pauxL_col,
                                                                pl_offd);
                           CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, pauxL, pauxU_row,
                                                                pu_offd);

                           CBlock_BxB_traits < _Flt >::ModifDia_BxB (_b_size, theta_loc,
                                                                     pl_offd, pu_offd,
                                                                     pa_dia + irow * b_2,
                                                                     pa_dia +
                                                                     kcolmn * b_2);

                           if (_params->b_write_file) {
                              if (irow == index1 || irow == index2 || irow == index3) {
                                 if (kcolmn == irow) {
                                    PrintArray (ffout_debug, " Mult L 2: U = ", b_2,
                                                pauxU);
                                    PrintArray (ffout_debug, "        2: L col = ", b_2,
                                                pauxL_col);
                                    PrintArray (ffout_debug, " Mult U 2: L = ", b_2,
                                                pauxL);
                                    PrintArray (ffout_debug, "        2: U row = ", b_2,
                                                pauxU_row);
                                    ffout_debug << "          :: Modif 2: " << endl;
                                    PrintArray (ffout_debug, " Irow   = ", b_2,
                                                pa_dia + irow * b_2);
                                    PrintArray (ffout_debug, " Kcolmn = ", b_2,
                                                pa_dia + kcolmn * b_2);
                                 }
                              }
                           }

                        }

                     }
                  }
               }

            }

            irwprv = irwpr1;

         }

// Update transposed structure

         irwprv = (int) ibegm[irow];

         while (irwprv != -1) {

            irwpr1 = (int) madj[irwprv];

            if (iv[irwprv] >= _ia_lu[irwprv + 1] - 1) {
               madj[irwprv] = -1;
            } else {
               j = (int) iv[irwprv] + 1;
               jcolmn = (int) _ja_lu[j];
               if (jcolmn < 0)
                  jcolmn = -jcolmn - 1;
               madj[irwprv] = ibegm[jcolmn];
               ibegm[jcolmn] = irwprv;
            }

            iv[irwprv]++;

            irwprv = irwpr1;

         }

// Perform filtering of the data and modify diagonal

         if (i >= _n_ini) {
            fcttype_loc = _params->fcttype;
            tau2_loc = 0.0e0;
         }

         {

            nlistcnew = 0;

            for (j = 0; j < nlistcloc; j++) {
               ind = (int) listc[j];
               if (ind != irow) {

                  jjchar = imaskchar[ind];

                  dfnorm =
                     CBlock_BxB_traits < _Flt >::FNormValue_BxB (_b_size,
                                                                 pfmaskc + ind * b_2_2,
                                                                 pfmaskc + ind * b_2_2 +
                                                                 b_2);

                  if (jjchar > 0 && ((fcttype_loc < 0 && dfnorm < tau2_loc)
                                     || (fcttype_loc >= 0
                                         && (jjchar > fcttype_loc || dfnorm < tau2_loc)))
                     ) {

                     CBlock_BxB_traits < _Flt >::ModifDia_BxB (_b_size, theta_loc,
                                                               pfmaskc + ind * b_2_2,
                                                               pfmaskc + ind * b_2_2 +
                                                               b_2, pa_dia + irow * b_2,
                                                               pa_dia + ind * b_2);

                  } else {
                     listc[nlistcnew] = ind;
                     nlistcnew++;
                  }
               } else {
                  listc[nlistcnew] = irow;
                  nlistcnew++;
               }
            }

            nlistcloc = nlistcnew;

         }

         if (_params->b_write_file) {
            if (irow == index1 || irow == index2 || irow == index3) {
               ffout_debug << "          Dia before modif add: " << endl;
               PrintArray (ffout_debug, " Diag L = ", b_2, pfmaskc + irow * b_2_2);
               PrintArray (ffout_debug, " Diag U = ", b_2, pfmaskc + irow * b_2_2 + b_2);
            }
         }

         CBlock_BxB_traits < _Flt >::AddModifDia_BxB (_b_size, pa_dia + irow * b_2,
                                                      pfmaskc + irow * b_2_2);
         CBlock_BxB_traits < _Flt >::AddModifDia_BxB (_b_size, pa_dia + irow * b_2,
                                                      pfmaskc + irow * b_2_2 + b_2);

         if (_params->b_write_file) {
            if (irow == index1 || irow == index2 || irow == index3) {
               ffout_debug << "          Dia after  modif add: " << endl;
               PrintArray (ffout_debug, " Diag L = ", b_2, pfmaskc + irow * b_2_2);
               PrintArray (ffout_debug, " Diag U = ", b_2, pfmaskc + irow * b_2_2 + b_2);
            }
         }
// Factorize current row and split into first and second order

         _Int *plistc = &listc[0];

         sort (plistc, plistc + nlistcloc);

         if (i < _n_ini) {

            if (fcttype_dia_loc == -1) {

               _nmodif +=
                  CBlock_BxB_traits < _Flt >::FctInv_BxB (pivmin_loc, _b_size,
                                                          pfmaskc + irow * b_2_2 + b_2,
                                                          pl_offd, pu_offd, pl_diag,
                                                          pu_diag, eigmin_temp,
                                                          eigmax_temp);

               CVector < _Flt >::CopyVector (b_2, pl_diag, pfmaskc + irow * b_2_2);
               CVector < _Flt >::CopyVector (b_2, pu_diag, pfmaskc + irow * b_2_2 + b_2);

               if (eigmin_temp < eigmin_att)
                  eigmin_att = eigmin_temp;
               if (eigmax_temp > eigmax_att)
                  eigmax_att = eigmax_temp;

            } else if (fcttype_dia_loc == 0) {

               int nmodif_loc = 0;

               nmodif_loc =
                  CBlock_BxB_traits < _Flt >::FctInv_BiQrd_BxB (pivmin_loc, _b_size,
                                                                pfmaskc + irow * b_2_2 +
                                                                b_2, pl_offd, pu_offd,
                                                                pl_diag, pu_diag, pWork,
                                                                eigmin_temp, eigmax_temp);

               _nmodif += nmodif_loc;

               CVector < _Flt >::CopyVector (b_2, pl_diag, pfmaskc + irow * b_2_2);
               CVector < _Flt >::CopyVector (b_2, pu_diag, pfmaskc + irow * b_2_2 + b_2);

               if (eigmin_temp < eigmin_att)
                  eigmin_att = eigmin_temp;
               if (eigmax_temp > eigmax_att)
                  eigmax_att = eigmax_temp;

            } else {

               int nmodif_loc = 0;

               nmodif_loc =
                  CBlock_BxB_traits < _Flt >::FctInv_Svd_BxB (pivmin_loc, _b_size,
                                                              pfmaskc + irow * b_2_2 +
                                                              b_2, pl_offd, pu_offd,
                                                              pl_diag, pu_diag, pWork,
                                                              pdWork, eigmin_temp,
                                                              eigmax_temp);

               _nmodif += nmodif_loc;

               CVector < _Flt >::CopyVector (b_2, pl_diag, pfmaskc + irow * b_2_2);
               CVector < _Flt >::CopyVector (b_2, pu_diag, pfmaskc + irow * b_2_2 + b_2);

               if (eigmin_temp < eigmin_att)
                  eigmin_att = eigmin_temp;
               if (eigmax_temp > eigmax_att)
                  eigmax_att = eigmax_temp;

            }

         }

         if (i < _n_ini) {

            for (j = 0; j < nlistcloc; j++) {
               ind = (int) listc[j];
               if (ind != irow) {

                  CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, pl_diag,
                                                       pfmaskc + ind * b_2_2 + b_2,
                                                       pu_offd);
                  CBlock_BxB_traits < _Flt >::MtM_BxB (_b_size, pu_diag,
                                                       pfmaskc + ind * b_2_2, pl_offd);

                  CVector < _Flt >::CopyVector (b_2, pl_offd, pfmaskc + ind * b_2_2);
                  CVector < _Flt >::CopyVector (b_2, pu_offd,
                                                pfmaskc + ind * b_2_2 + b_2);

                  dfnorm =
                     CBlock_BxB_traits < _Flt >::FNormValue_BxB (_b_size,
                                                                 pfmaskc + ind * b_2_2,
                                                                 pfmaskc + ind * b_2_2 +
                                                                 b_2);
                  jjchar = imaskchar[ind];

                  if (jjchar > 0 && dfnorm < tau1_loc) {
                     listc[j] = (_Int) (-ind - 1);
                  }
               }
            }

         }
// Store computed row elems

         _a_lu.resize ((nzja_lu + nlistcloc) * b_2_2);

         _Flt *p_a_lu_temp = &_a_lu[0];

         for (j = 0; j < nlistcloc; j++) {
            ind = (int) listc[j];
            _ja_lu.push_back (ind);
            ind1 = ind;
            if (ind1 < 0)
               ind1 = -ind1 - 1;
            jjchar = imaskchar[ind1];
            _ja_char_lu.push_back (jjchar);
            CVector < _Flt >::CopyVector (b_2_2, pfmaskc + ind1 * b_2_2,
                                          p_a_lu_temp + (nzja_lu + j) * b_2_2);
         }

         nzja_lu += nlistcloc;

         _ia_lu[i + 1] = nzja_lu;

// Add current row into the transposed structures

         if (nlistcloc > 1) {
            if (i < _n_ini) {
               iv[irow] = _ia_lu[i] + 1;
               ind = (int) listc[1];
               if (ind < 0)
                  ind = -ind - 1;
               madj[irow] = ibegm[ind];
               ibegm[ind] = (_Int) irow;
            } else {
               iv[irow] = _ia_lu[i + 1] + 1;
            }
         }

      }

      _eigmin_att = eigmin_att;
      _eigmax_att = eigmax_att;

// Condense the result (filter second order elems)

      int nzja_new = 0;

      for (i = 0; i < nzja_lu; i++) {
         if (_ja_lu[i] >= 0)
            nzja_new++;
      }

      vector < _Int > jlu_cnd (nzja_new + 1);
      vector < char >jlu_char_cnd (nzja_new + 1);
      vector < _Flt > lu_cnd (b_2_2 * nzja_new + 1);

      _Int *pjlu_cnd = &jlu_cnd[0];
      char *pjlu_char_cnd = &jlu_char_cnd[0];
      _Flt *plu_cnd = &lu_cnd[0];

      nzja_new = 0;

      listc[0] = 0;

      _Flt *p_a_lu_temp = NULL;
      if (_a_lu.size () > 0)
         p_a_lu_temp = &_a_lu[0];

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_lu[i]; j < _ia_lu[i + 1]; j++) {
            if (_ja_lu[j] >= 0) {
               pjlu_cnd[nzja_new] = _ja_lu[j];
               pjlu_char_cnd[nzja_new] = _ja_char_lu[j];
               CVector < _Flt >::CopyVector (b_2_2, p_a_lu_temp + j * b_2_2,
                                             plu_cnd + nzja_new * b_2_2);
               nzja_new++;
            }
         }
         listc[i + 1] = nzja_new;
//         if (_params->b_write_file) {
//            ffout_debug << " Irow = " << i << endl;
//            PrintArray (ffout_debug, " Ja ",(int)(listc[i+1]-listc[i]),_ja_lu.data()+listc[i]);
//            PrintArray (ffout_debug, " A ",(int)((listc[i+1]-listc[i])*b_2_2),plu_cnd+listc[i]*b_2_2);
//         }
      }

      for (i = 0; i <= _n; i++)
         _ia_lu[i] = listc[i];

      _ja_lu.swap (jlu_cnd);
      _ja_char_lu.swap (jlu_char_cnd);
      _a_lu.swap (lu_cnd);

   }

//
// Perform one step fct for dinamically splitted and ordered triangular matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::DiaSplitFctSchur (SParams * _params, double _dia_split, int _b_size, int _n,
                                const vector < _Int > &_ia_au,
                                const vector < _Int > &_ja_au,
                                const vector < char >&_ja_char_au,
                                const vector < _Flt > &_a_au, int *_order, int &_n_ini,
                                vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                                vector < char >&_ja_char_u, vector < _Flt > &_a_u,
                                vector < _Int > &_ia_u_schur,
                                vector < _Int > &_ja_u_schur,
                                vector < char >&_ja_char_u_schur,
                                vector < _Flt > &_a_u_schur, int &_nmodif,
                                double &_eigmin_att, double &_eigmax_att)
   {

// Open matrix data structures

      int b_2 = _b_size * _b_size;
      int b_2_2 = 2 * b_2;

      const _Int *p_ia_au = NULL;
      const _Int *p_ja_au = NULL;
      const char *p_ja_char_au = NULL;
      const _Flt *p_a_au = NULL;

      if (_ia_au.size () > 0)
         p_ia_au = &_ia_au[0];
      if (_ja_au.size () > 0)
         p_ja_au = &_ja_au[0];
      if (_ja_char_au.size () > 0)
         p_ja_char_au = &_ja_char_au[0];
      if (_a_au.size () > 0)
         p_a_au = &_a_au[0];

// Compute list of active rows

      CVectorData < _Flt > work (20 * b_2);
      CVectorData < double >dwork (20 * b_2);

      _Flt *pSv = work.Ptr ();
      _Flt *pU_matr = pSv + _b_size;
      _Flt *pV_matr = pU_matr + b_2;
      _Flt *pwork1 = pV_matr + b_2;

      double *pdwork = dwork.Ptr ();

      CVectorData < int >imask (_n);
      int *pimask = imask.Ptr ();

      int i, j, jj, ind_dia;
      int kii;
      _Flt aux;

      if (_dia_split < 0.0e0) {
         for (i = 0; i < _n; i++)
            pimask[i] = 1;
      } else {

         for (i = 0; i < _n; i++) {
            ind_dia = -1;
            for (j = p_ia_au[i]; j < p_ia_au[i + 1]; j++) {
               jj = p_ja_au[j];
               if (jj == i)
                  ind_dia = (int) j;
            }
            if (ind_dia >= 0) {
               CBlock_BxB_traits < _Flt >::ComputeQrd_BxB (_b_size, _b_size,
                                                           p_a_au + ind_dia * b_2_2 + b_2,
                                                           pU_matr, pV_matr, pwork1);

               pimask[i] = 1;
               for (kii = 0; kii < _b_size; kii++) {
                  aux = pV_matr[_b_size * kii + kii];
                  if (aux < 0.0e0)
                     aux = -aux;
                  if (aux < _dia_split) {
                     pimask[i] = 0;
                  }
               }
            } else {
               pimask[i] = 0;
            }
         }

      }

// Split sparsity

      vector < int >order (_n + 1);
      int *porder = &order[0];

      _n_ini = 0;

      for (i = 0; i < _n; i++) {
         if (pimask[i] > 0)
            _n_ini++;
      }

      int n0 = 0;
      int n1 = _n_ini;

      for (i = 0; i < _n; i++) {
         if (pimask[i] > 0) {
            porder[i] = n0;
            n0++;
         } else {
            porder[i] = n1;
            n1++;
         }
      }

      vector < _Int > ia_au_ord;
      vector < _Int > ja_au_ord;

      CFct_impl < _Int, _Flt >::ReorderMatrixSp (_n, order, _ia_au, _ja_au, ia_au_ord,
                                                 ja_au_ord);

      _Int *pia_au_ord = NULL;
      _Int *pja_au_ord = NULL;

      if (ia_au_ord.size () > 0)
         pia_au_ord = &ia_au_ord[0];
      if (ja_au_ord.size () > 0)
         pja_au_ord = &ja_au_ord[0];

      int nzja0 = 0;

      for (i = 0; i < _n_ini; i++) {
         for (j = pia_au_ord[i]; j < pia_au_ord[i + 1]; j++) {
            jj = pja_au_ord[j];
            if (jj < _n_ini)
               nzja0++;
         }
      }

      vector < _Int > ia_au_ord0 (_n_ini + 1);
      vector < _Int > ja_au_ord0 (nzja0 + 1);

      _Int *pia_au_ord0 = &ia_au_ord0[0];
      _Int *pja_au_ord0 = &ja_au_ord0[0];

      nzja0 = 0;

      pia_au_ord0[0] = 0;

      for (i = 0; i < _n_ini; i++) {
         for (j = pia_au_ord[i]; j < pia_au_ord[i + 1]; j++) {
            jj = pja_au_ord[j];
            if (jj < _n_ini) {
               pja_au_ord0[nzja0] = jj;
               nzja0++;
            }
         }
         pia_au_ord0[i + 1] = nzja0;
      }

// Find optimal ordering

      vector < int >order0;

      CFct_impl < _Int, _Flt >::ComputeOptimalOrder (1, _n_ini, ia_au_ord0, ja_au_ord0,
                                                     order0);

      order0.resize (_n + 1);

      int *porder0 = &order0[0];

      for (i = _n_ini; i < _n; i++)
         porder0[i] = i;

      vector < int >iorder (_n + 1);
      int *piorder = &iorder[0];

      for (i = 0; i < _n; i++)
         piorder[porder[i]] = i;

      int iold;

      for (i = 0; i < _n; i++) {
         iold = piorder[i];
         _order[iold] = porder0[i];
         porder[iold] = porder0[i];
      }

// Free work arrays

      {
         vector < _Int > ia_au_ord_temp;
         vector < _Int > ja_au_ord_temp;
         vector < _Int > ia_au_ord0_temp;
         vector < _Int > ja_au_ord0_temp;
         ia_au_ord.swap (ia_au_ord_temp);
         ja_au_ord.swap (ja_au_ord_temp);
         ia_au_ord0.swap (ia_au_ord0_temp);
         ja_au_ord0.swap (ja_au_ord0_temp);
      }

// Split into L and U parts and reorder and back combine

      vector < _Int > ia_alu_ord_pair;
      vector < _Int > ja_alu_ord_pair;
      vector < char >ja_char_alu_ord_pair;
      vector < _Flt > a_alu_ord_pair;

      {

         vector < _Int > ia_l;
         vector < _Int > ja_l;
         vector < char >ja_char_l;
         vector < _Flt > a_l;

         vector < _Int > ia_u;
         vector < _Int > ja_u;
         vector < char >ja_char_u;
         vector < _Flt > a_u;

         CFct_bxb_impl < _Int, _Flt >::SplitLUPairs (true, _b_size, _n, _ia_au, _ja_au,
                                                     _ja_char_au, _a_au, ia_l, ja_l,
                                                     ja_char_l, a_l, ia_u, ja_u,
                                                     ja_char_u, a_u);

         vector < _Int > ia_lt;
         vector < _Int > ja_lt;
         vector < char >ja_char_lt;
         vector < _Flt > a_lt;

         CFct_bxb_impl < _Int, _Flt >::Transpose (true, _b_size, _n, ia_l, ja_l,
                                                  ja_char_l, a_l, ia_lt, ja_lt,
                                                  ja_char_lt, a_lt);

         {
            vector < _Int > ia_l_temp;
            vector < _Int > ja_l_temp;
            vector < char >ja_char_l_temp;
            vector < _Flt > a_l_temp;
            ia_l.swap (ia_l_temp);
            ja_l.swap (ja_l_temp);
            ja_char_l.swap (ja_char_l_temp);
            a_l.swap (a_l_temp);
         }

         vector < _Int > ia_alu;
         vector < _Int > ja_alu;
         vector < char >ja_char_alu;
         vector < _Flt > a_alu;

         CFct_bxb_impl < _Int, _Flt >::CombineLU (true, _b_size, _n, ia_lt, ja_lt,
                                                  ja_char_lt, a_lt, ia_u, ja_u, ja_char_u,
                                                  a_u, ia_alu, ja_alu, ja_char_alu,
                                                  a_alu);

         {
            vector < _Int > ia_lt_temp;
            vector < _Int > ja_lt_temp;
            vector < char >ja_char_lt_temp;
            vector < _Flt > a_lt_temp;
            ia_lt.swap (ia_lt_temp);
            ja_lt.swap (ja_lt_temp);
            ja_char_lt.swap (ja_char_lt_temp);
            a_lt.swap (a_lt_temp);
            vector < _Int > ia_u_temp;
            vector < _Int > ja_u_temp;
            vector < char >ja_char_u_temp;
            vector < _Flt > a_u_temp;
            ia_u.swap (ia_u_temp);
            ja_u.swap (ja_u_temp);
            ja_char_u.swap (ja_char_u_temp);
            a_u.swap (a_u_temp);
         }

         vector < _Int > ia_alu_ord;
         vector < _Int > ja_alu_ord;
         vector < char >ja_char_alu_ord;
         vector < _Flt > a_alu_ord;

         CFct_bxb_impl < _Int, _Flt >::ReorderMatrix (true, false, _b_size, _n, order,
                                                      ia_alu, ja_alu, ja_char_alu, a_alu,
                                                      ia_alu_ord, ja_alu_ord,
                                                      ja_char_alu_ord, a_alu_ord);

         {
            vector < _Int > ia_alu_temp;
            vector < _Int > ja_alu_temp;
            vector < char >ja_char_alu_temp;
            vector < _Flt > a_alu_temp;
            ia_alu.swap (ia_alu_temp);
            ja_alu.swap (ja_alu_temp);
            ja_char_alu.swap (ja_char_alu_temp);
            a_alu.swap (a_alu_temp);
         }

         CFct_bxb_impl < _Int, _Flt >::SplitLU (true, _b_size, _n, ia_alu_ord, ja_alu_ord,
                                                ja_char_alu_ord, a_alu_ord, ia_lt, ja_lt,
                                                ja_char_lt, a_lt, ia_u, ja_u, ja_char_u,
                                                a_u);

         {
            vector < _Int > ia_alu_temp;
            vector < _Int > ja_alu_temp;
            vector < char >ja_char_alu_temp;
            vector < _Flt > a_alu_temp;
            ia_alu_ord.swap (ia_alu_temp);
            ja_alu_ord.swap (ja_alu_temp);
            ja_char_alu_ord.swap (ja_char_alu_temp);
            a_alu_ord.swap (a_alu_temp);
         }

         CFct_bxb_impl < _Int, _Flt >::Transpose (true, _b_size, _n, ia_lt, ja_lt,
                                                  ja_char_lt, a_lt, ia_l, ja_l, ja_char_l,
                                                  a_l);

         {
            vector < _Int > ia_l_temp;
            vector < _Int > ja_l_temp;
            vector < char >ja_char_l_temp;
            vector < _Flt > a_l_temp;
            ia_lt.swap (ia_l_temp);
            ja_lt.swap (ja_l_temp);
            ja_char_lt.swap (ja_char_l_temp);
            a_lt.swap (a_l_temp);
         }

         CFct_bxb_impl < _Int, _Flt >::CombineLUPairs (true, _b_size, _n, ia_l, ja_l,
                                                       ja_char_l, a_l, ia_u, ja_u,
                                                       ja_char_u, a_u, ia_alu_ord_pair,
                                                       ja_alu_ord_pair,
                                                       ja_char_alu_ord_pair,
                                                       a_alu_ord_pair);

      }

      _Int *pia_alu_ord_pair = NULL;
      _Int *pja_alu_ord_pair = NULL;
      char *pja_char_alu_ord_pair = NULL;
      _Flt *pa_alu_ord_pair = NULL;

      if (ia_alu_ord_pair.size () > 0)
         pia_alu_ord_pair = &ia_alu_ord_pair[0];
      if (ja_alu_ord_pair.size () > 0)
         pja_alu_ord_pair = &ja_alu_ord_pair[0];
      if (ja_char_alu_ord_pair.size () > 0)
         pja_char_alu_ord_pair = &ja_char_alu_ord_pair[0];
      if (a_alu_ord_pair.size () > 0)
         pa_alu_ord_pair = &a_alu_ord_pair[0];

// Split into main and schur parts

      int n_last = _n - _n_ini;

      int nzja_ini = pia_alu_ord_pair[_n_ini];
      int nzja_last = pia_alu_ord_pair[_n] - pia_alu_ord_pair[_n_ini];

      vector < _Int > ia_ini (_n_ini + 1);
      vector < _Int > ja_ini (nzja_ini + 1);
      vector < char >ja_char_ini (nzja_ini + 1);
      vector < _Flt > a_ini (nzja_ini * b_2_2 + 1);

      _Int *pia_ini = NULL;
      _Int *pja_ini = NULL;
      char *pja_char_ini = NULL;
      _Flt *pa_ini = NULL;

      if (ia_ini.size () > 0)
         pia_ini = &ia_ini[0];
      if (ja_ini.size () > 0)
         pja_ini = &ja_ini[0];
      if (ja_char_ini.size () > 0)
         pja_char_ini = &ja_char_ini[0];
      if (a_ini.size () > 0)
         pa_ini = &a_ini[0];

      for (i = 0; i <= _n_ini; i++)
         pia_ini[i] = pia_alu_ord_pair[i];
      for (i = 0; i < nzja_ini; i++)
         pja_ini[i] = pja_alu_ord_pair[i];
      for (i = 0; i < nzja_ini; i++)
         pja_char_ini[i] = pja_char_alu_ord_pair[i];
      for (i = 0; i < nzja_ini * b_2_2; i++)
         pa_ini[i] = pa_alu_ord_pair[i];

      vector < _Int > ia_last (n_last + 1);
      vector < _Int > ja_last (nzja_last + 1);
      vector < char >ja_char_last (nzja_last + 1);
      vector < _Flt > a_last (nzja_last * b_2_2 + 1);

      _Int *pia_last = NULL;
      _Int *pja_last = NULL;
      char *pja_char_last = NULL;
      _Flt *pa_last = NULL;

      if (ia_last.size () > 0)
         pia_last = &ia_last[0];
      if (ja_last.size () > 0)
         pja_last = &ja_last[0];
      if (ja_char_last.size () > 0)
         pja_char_last = &ja_char_last[0];
      if (a_last.size () > 0)
         pa_last = &a_last[0];

      for (i = 0; i <= n_last; i++)
         pia_last[i] = pia_alu_ord_pair[i + _n_ini] - pia_alu_ord_pair[_n_ini];
      for (i = 0; i < nzja_last; i++)
         pja_last[i] = pja_alu_ord_pair[nzja_ini + i] - _n_ini;
      for (i = 0; i < nzja_last; i++)
         pja_char_last[i] = pja_char_alu_ord_pair[nzja_ini + i];
      int ishift = nzja_ini * b_2_2;
      for (i = 0; i < nzja_last * b_2_2; i++)
         pa_last[i] = pa_alu_ord_pair[ishift + i];

// Fct

      vector < _Int > ia_lu_ini;
      vector < _Int > ja_lu_ini;
      vector < char >ja_lu_char_ini;
      vector < _Flt > a_lu_ini;

      CFct_bxb_impl < _Int, _Flt >::FctBlockIlu2Degree (_params, _b_size, _n, _n_ini,
                                                        ia_ini, ja_ini, ja_char_ini,
                                                        a_ini, ia_lu_ini, ja_lu_ini,
                                                        ja_lu_char_ini, a_lu_ini, _nmodif,
                                                        _eigmin_att, _eigmax_att);

      _Int *pia_lu_ini = NULL;
      _Int *pja_lu_ini = NULL;
      char *pja_lu_char_ini = NULL;
      _Flt *pa_lu_ini = NULL;

      if (ia_lu_ini.size () > 0)
         pia_lu_ini = &ia_lu_ini[0];
      if (ja_lu_ini.size () > 0)
         pja_lu_ini = &ja_lu_ini[0];
      if (ja_lu_char_ini.size () > 0)
         pja_lu_char_ini = &ja_lu_char_ini[0];
      if (a_lu_ini.size () > 0)
         pa_lu_ini = &a_lu_ini[0];

// Finalize U and Schur complement

      int nzja_u = ia_lu_ini[_n_ini];
      int nzja_u2 = ia_lu_ini[_n] - ia_lu_ini[_n_ini];

      _ia_u.resize (_n_ini + 1);
      _ja_u.resize (nzja_u + 1);
      _ja_char_u.resize (nzja_u + 1);
      _a_u.resize (nzja_u * b_2_2 + 1);

      _Int *p_ia_u = NULL;
      _Int *p_ja_u = NULL;
      char *p_ja_char_u = NULL;
      _Flt *p_a_u = NULL;

      if (_ia_u.size () > 0)
         p_ia_u = &_ia_u[0];
      if (_ja_u.size () > 0)
         p_ja_u = &_ja_u[0];
      if (_ja_char_u.size () > 0)
         p_ja_char_u = &_ja_char_u[0];
      if (_a_u.size () > 0)
         p_a_u = &_a_u[0];

      for (i = 0; i <= _n_ini; i++)
         p_ia_u[i] = pia_lu_ini[i];
      for (i = 0; i < nzja_u; i++)
         p_ja_u[i] = pja_lu_ini[i];
      for (i = 0; i < nzja_u; i++)
         p_ja_char_u[i] = pja_lu_char_ini[i];
      for (i = 0; i < nzja_u * b_2_2; i++)
         p_a_u[i] = pa_lu_ini[i];

      vector < _Int > ia_u2 (n_last + 1);
      vector < _Int > ja_u2 (nzja_u2 + 1);
      vector < char >ja_char_u2 (nzja_u2 + 1);
      vector < _Flt > a_u2 (nzja_u2 * b_2_2 + 1);

      _Int *pia_u2 = NULL;
      _Int *pja_u2 = NULL;
      char *pja_char_u2 = NULL;
      _Flt *pa_u2 = NULL;

      if (ia_u2.size () > 0)
         pia_u2 = &ia_u2[0];
      if (ja_u2.size () > 0)
         pja_u2 = &ja_u2[0];
      if (ja_char_u2.size () > 0)
         pja_char_u2 = &ja_char_u2[0];
      if (a_u2.size () > 0)
         pa_u2 = &a_u2[0];

      for (i = 0; i <= n_last; i++)
         pia_u2[i] = pia_lu_ini[i + _n_ini] - pia_lu_ini[_n_ini];
      for (i = 0; i < nzja_u2; i++)
         pja_u2[i] = pja_lu_ini[nzja_u + i] - _n_ini;
      for (i = 0; i < nzja_u2; i++)
         pja_char_u2[i] = pja_lu_char_ini[nzja_u + i];
      ishift = nzja_u * b_2_2;
      for (i = 0; i < nzja_u2 * b_2_2; i++)
         pa_u2[i] = pa_lu_ini[ishift + i];

      CFct_bxb_impl < _Int, _Flt >::AddMatrices (true, true, _b_size, n_last, ia_last,
                                                 ja_last, ja_char_last, a_last, ia_u2,
                                                 ja_u2, ja_char_u2, a_u2, _ia_u_schur,
                                                 _ja_u_schur, _ja_char_u_schur,
                                                 _a_u_schur);

   }

//
// Perform one step fct for dinamically splitted and ordered triangular matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::DiaSplitFctSchur (SParams * _params, double _dia_split, int _b_size,
                                int _n_rows, int _n, int _n_ini_max,
                                const vector < _Int > &_ia_au,
                                const vector < _Int > &_ja_au,
                                const vector < char >&_ja_char_au,
                                const vector < _Flt > &_a_au, int *_order, int &_n_ini,
                                int &_n_bad_diag, vector < _Int > &_ia_u,
                                vector < _Int > &_ja_u, vector < char >&_ja_char_u,
                                vector < _Flt > &_a_u, vector < _Int > &_ia_u_schur,
                                vector < _Int > &_ja_u_schur,
                                vector < char >&_ja_char_u_schur,
                                vector < _Flt > &_a_u_schur, int &_nmodif,
                                double &_eigmin_att, double &_eigmax_att)
   {

//      ofstream *pfout_debug = _params->pfout;

      if (_n_rows != _n_ini_max && _n_rows != _n) {
         cout <<
            " CFct_bxb_impl < _Int,_Flt >::DiaSplitFctSchur: error in number of rows ! "
            << endl;
         throw
            " CFct_bxb_impl < _Int,_Flt >::DiaSplitFctSchur: error in number of rows ! ";
      }
// Open matrix data structures

      int b_2 = _b_size * _b_size;
      int b_2_2 = 2 * b_2;

      const _Int *p_ia_au = _ia_au.data ();
      const _Int *p_ja_au = _ja_au.data ();
      const char *p_ja_char_au = _ja_char_au.data ();
      const _Flt *p_a_au = _a_au.data ();

// Compute list of active rows

      CVectorData < _Flt > work (30 * b_2);
      CVectorData < double >dwork (30 * b_2);

      _Flt *pSv = work.Ptr ();
      _Flt *pU_matr = pSv + _b_size * 2;
      _Flt *pV_matr = pU_matr + b_2;
      _Flt *pwork1 = pV_matr + b_2;
      double *pdwork = dwork.Ptr ();

      CVectorData < int >imask (_n);
      int *pimask = imask.Ptr ();

      int i, j, jj, ind_dia;
      int kii;
      _Flt aux;

      if (_dia_split < 0.0e0) {
         for (i = 0; i < _n_ini_max; i++)
            pimask[i] = 1;
         for (i = _n_ini_max; i < _n; i++)
            pimask[i] = 0;
      } else {

         for (i = 0; i < _n_ini_max; i++) {
            ind_dia = -1;
            for (j = (int) p_ia_au[i]; j < p_ia_au[i + 1]; j++) {
               jj = (int) p_ja_au[j];
               if (jj == i)
                  ind_dia = (int) j;
            }
            if (ind_dia >= 0) {
               CBlock_BxB_traits < _Flt >::BiDiagonalize_BxB (_b_size,
                                                              p_a_au + ind_dia * b_2_2 +
                                                              b_2, pSv, pU_matr, pV_matr,
                                                              pwork1);
               pimask[i] = 1;
               for (kii = 0; kii < _b_size; kii++) {
                  aux = pSv[kii];
                  if (aux < 0.0e0)
                     aux = -aux;
                  if (aux < _dia_split) {
                     pimask[i] = 0;
                  }
               }
            } else {
               pimask[i] = 0;
            }
         }
         for (i = _n_ini_max; i < _n; i++)
            pimask[i] = 0;

      }

// Split sparsity

      vector < int >order (_n + 1);
      int *porder = &order[0];

      _n_ini = 0;

      for (i = 0; i < _n_ini_max; i++) {
         if (pimask[i] > 0)
            _n_ini++;
      }

      _n_bad_diag = _n_ini_max - _n_ini;

// Fast return

      if (_n_ini == 0)
         _n_ini = _n_ini_max;

      if (_n_ini == _n_ini_max) {

// Split by rows

         vector < _Int > ia_ini;
         vector < _Int > ja_ini;
         vector < char >ja_char_ini;
         vector < _Flt > a_ini;

         vector < _Int > ia_last;
         vector < _Int > ja_last;
         vector < char >ja_char_last;
         vector < _Flt > a_last;

         int n_last_base = _n_rows - _n_ini_max;

         if (n_last_base > 0) {

            CFct_bxb_impl < _Int, _Flt >::SubmatricesByRows (true, true, _b_size, _n_rows,
                                                             _n_ini_max, _ia_au, _ja_au,
                                                             _ja_char_au, _a_au, ia_ini,
                                                             ja_ini, ja_char_ini, a_ini,
                                                             ia_last, ja_last,
                                                             ja_char_last, a_last);

            int nzja_last_base = (int) ia_last[n_last_base];

            for (i = 0; i < nzja_last_base; i++)
               ja_last[i] = ja_last[i] - _n_ini_max;

         }
// Fct

         vector < _Int > ia_u_fct;
         vector < _Int > ja_u_fct;
         vector < char >ja_char_u_fct;
         vector < _Flt > a_u_fct;

         int iparam1 = _params->iparam1;
         int iparam2 = _params->iparam2;
         int iparam3 = _params->iparam3;
/*
         if (_n_ini_max != _n) {

            _params->b_write_file = true;

            sprintf (_params->name_file,"ChkDiaSplitSchur_%i_%i.dat",_n,_n_ini_max);

            _params->iparam1 = _n_ini_max;
            _params->iparam2 = -1;
            _params->iparam3 = -1;

         }
*/
         if (n_last_base > 0) {

            CFct_bxb_impl < _Int, _Flt >::FctBlockIlu2Degree (_params, _b_size, _n,
                                                              _n_ini_max, ia_ini, ja_ini,
                                                              ja_char_ini, a_ini,
                                                              ia_u_fct, ja_u_fct,
                                                              ja_char_u_fct, a_u_fct,
                                                              _nmodif, _eigmin_att,
                                                              _eigmax_att);

         } else {

            CFct_bxb_impl < _Int, _Flt >::FctBlockIlu2Degree (_params, _b_size, _n,
                                                              _n_ini_max, _ia_au, _ja_au,
                                                              _ja_char_au, _a_au,
                                                              ia_u_fct, ja_u_fct,
                                                              ja_char_u_fct, a_u_fct,
                                                              _nmodif, _eigmin_att,
                                                              _eigmax_att);

         }
/*
         if (_n_ini_max != _n) {

            _params->b_write_file = false;
            _params->iparam1 = iparam1;
            _params->iparam2 = iparam2;
            _params->iparam3 = iparam3;

         }
*/

// Split and shift

         CFct_bxb_impl < _Int, _Flt >::SubmatricesByRows (true, true, _b_size, _n,
                                                          _n_ini_max, ia_u_fct, ja_u_fct,
                                                          ja_char_u_fct, a_u_fct, _ia_u,
                                                          _ja_u, _ja_char_u, _a_u,
                                                          _ia_u_schur, _ja_u_schur,
                                                          _ja_char_u_schur, _a_u_schur);

         int n_last = _n - _n_ini_max;
         int nzja_U_Schur = (int) _ia_u_schur[n_last];

         for (i = 0; i < nzja_U_Schur; i++)
            _ja_u_schur[i] = _ja_u_schur[i] - _n_ini_max;

// Add

         if (n_last_base > 0) {

            vector < _Int > ia_sum;
            vector < _Int > ja_sum;
            vector < char >ja_char_sum;
            vector < _Flt > a_sum;

            CFct_bxb_impl < _Int, _Flt >::AddMatrices (true, true, _b_size, n_last,
                                                       ia_last, ja_last, ja_char_last,
                                                       a_last, _ia_u_schur, _ja_u_schur,
                                                       _ja_char_u_schur, _a_u_schur,
                                                       ia_sum, ja_sum, ja_char_sum,
                                                       a_sum);

            _ia_u_schur.swap (ia_sum);
            _ja_u_schur.swap (ja_sum);
            _ja_char_u_schur.swap (ja_char_sum);
            _a_u_schur.swap (a_sum);

         }

         for (i = 0; i < _n; i++)
            _order[i] = i;

         return;

      }

      int n0 = 0;
      int n1 = _n_ini;

      for (i = 0; i < _n; i++) {
         if (pimask[i] > 0) {
            porder[i] = n0;
            n0++;
         } else {
            porder[i] = n1;
            n1++;
         }
      }

      vector < _Int > ia_au_ord;
      vector < _Int > ja_au_ord;

      if (_n == _n_rows) {
         CFct_impl < _Int, _Flt >::ReorderMatrixSp (_n, order, _ia_au, _ja_au, ia_au_ord,
                                                    ja_au_ord);
      } else {

         vector < _Int > ia_au_ext (_n + 1);
         for (i = 0; i <= _n_rows; i++)
            ia_au_ext[i] = _ia_au[i];
         for (i = _n_rows + 1; i <= _n; i++)
            ia_au_ext[i] = ia_au_ext[i - 1];

         CFct_impl < _Int, _Flt >::ReorderMatrixSp (_n, order, ia_au_ext, _ja_au,
                                                    ia_au_ord, ja_au_ord);
      }

      _Int *pia_au_ord = NULL;
      _Int *pja_au_ord = NULL;

      if (ia_au_ord.size () > 0)
         pia_au_ord = &ia_au_ord[0];
      if (ja_au_ord.size () > 0)
         pja_au_ord = &ja_au_ord[0];

      int nzja0 = 0;

      for (i = 0; i < _n_ini; i++) {
         for (j = (int) pia_au_ord[i]; j < pia_au_ord[i + 1]; j++) {
            jj = (int) pja_au_ord[j];
            if (jj < _n_ini)
               nzja0++;
         }
      }

      vector < _Int > ia_au_ord0 (_n_ini + 1);
      vector < _Int > ja_au_ord0 (nzja0 + 1);

      _Int *pia_au_ord0 = &ia_au_ord0[0];
      _Int *pja_au_ord0 = &ja_au_ord0[0];

      nzja0 = 0;

      pia_au_ord0[0] = 0;

      for (i = 0; i < _n_ini; i++) {
         for (j = (int) pia_au_ord[i]; j < pia_au_ord[i + 1]; j++) {
            jj = (int) pja_au_ord[j];
            if (jj < _n_ini) {
               pja_au_ord0[nzja0] = jj;
               nzja0++;
            }
         }
         pia_au_ord0[i + 1] = nzja0;
      }

// Find optimal ordering

      vector < int >order0;

      CFct_impl < _Int, _Flt >::ComputeOptimalOrder (1, _n_ini, ia_au_ord0, ja_au_ord0,
                                                     order0);

      order0.resize (_n + 1);

      int *porder0 = &order0[0];

      for (i = _n_ini; i < _n; i++)
         porder0[i] = i;

      vector < int >iorder (_n + 1);
      int *piorder = &iorder[0];

      for (i = 0; i < _n; i++)
         piorder[porder[i]] = i;

      int iold;

      for (i = 0; i < _n; i++) {
         iold = piorder[i];
         _order[iold] = porder0[i];
         porder[iold] = porder0[i];
      }

// Free work arrays

      {
         vector < _Int > ia_au_ord_temp;
         vector < _Int > ja_au_ord_temp;
         vector < _Int > ia_au_ord0_temp;
         vector < _Int > ja_au_ord0_temp;
         ia_au_ord.swap (ia_au_ord_temp);
         ja_au_ord.swap (ja_au_ord_temp);
         ia_au_ord0.swap (ia_au_ord0_temp);
         ja_au_ord0.swap (ja_au_ord0_temp);
      }

// Split into L and U parts and reorder and back combine

      vector < _Int > ia_alu_ord_pair;
      vector < _Int > ja_alu_ord_pair;
      vector < char >ja_char_alu_ord_pair;
      vector < _Flt > a_alu_ord_pair;

      {

         vector < _Int > ia_l;
         vector < _Int > ja_l;
         vector < char >ja_char_l;
         vector < _Flt > a_l;

         vector < _Int > ia_u;
         vector < _Int > ja_u;
         vector < char >ja_char_u;
         vector < _Flt > a_u;

         {

            if (_n == _n_rows) {
               CFct_bxb_impl < _Int, _Flt >::SplitLUPairs (true, _b_size, _n, _ia_au,
                                                           _ja_au, _ja_char_au, _a_au,
                                                           ia_l, ja_l, ja_char_l, a_l,
                                                           ia_u, ja_u, ja_char_u, a_u);
            } else {

               int nzja = (int) _ia_au[_n_rows];

               vector < _Int > ia_au_ext (_n + 1);

               for (i = 0; i <= _n_rows; i++)
                  ia_au_ext[i] = _ia_au[i];
               for (i = _n_rows + 1; i <= _n; i++)
                  ia_au_ext[i] = ia_au_ext[i - 1] + 1;

               int nzja_ext = (int) ia_au_ext[_n];

               vector < _Int > ja_au_ext (nzja_ext + 1);
               vector < char >jachar_au_ext (nzja_ext + 1);
               vector < _Flt > a_au_ext (nzja_ext * b_2_2 + 1);

               const _Flt *p_a_au = &_a_au[0];
               _Flt *pa_au_ext = &a_au_ext[0];

               for (i = 0; i < nzja; i++)
                  ja_au_ext[i] = _ja_au[i];
               for (i = 0; i < nzja; i++)
                  jachar_au_ext[i] = _ja_char_au[i];

               CVector < _Flt >::CopyVector (nzja * b_2_2, p_a_au, pa_au_ext);

               int iloc;

               for (i = _n_rows; i < _n; i++) {
                  iloc = i - _n_rows;
                  ja_au_ext[nzja + iloc] = i;
                  jachar_au_ext[nzja + iloc] = 0;
               }

               CVector < _Flt >::SetByZeroes ((_n - _n_rows) * b_2_2,
                                              pa_au_ext + nzja * b_2_2);

               CFct_bxb_impl < _Int, _Flt >::SplitLUPairs (true, _b_size, _n, ia_au_ext,
                                                           ja_au_ext, jachar_au_ext,
                                                           a_au_ext, ia_l, ja_l,
                                                           ja_char_l, a_l, ia_u, ja_u,
                                                           ja_char_u, a_u);

            }

         }

         vector < _Int > ia_lt;
         vector < _Int > ja_lt;
         vector < char >ja_char_lt;
         vector < _Flt > a_lt;

         CFct_bxb_impl < _Int, _Flt >::Transpose (true, _b_size, _n, ia_l, ja_l,
                                                  ja_char_l, a_l, ia_lt, ja_lt,
                                                  ja_char_lt, a_lt);

         {
            vector < _Int > ia_l_temp;
            vector < _Int > ja_l_temp;
            vector < char >ja_char_l_temp;
            vector < _Flt > a_l_temp;
            ia_l.swap (ia_l_temp);
            ja_l.swap (ja_l_temp);
            ja_char_l.swap (ja_char_l_temp);
            a_l.swap (a_l_temp);
         }

         vector < _Int > ia_alu;
         vector < _Int > ja_alu;
         vector < char >ja_char_alu;
         vector < _Flt > a_alu;

         CFct_bxb_impl < _Int, _Flt >::CombineLU (true, _b_size, _n, ia_lt, ja_lt,
                                                  ja_char_lt, a_lt, ia_u, ja_u, ja_char_u,
                                                  a_u, ia_alu, ja_alu, ja_char_alu,
                                                  a_alu);

         {
            vector < _Int > ia_lt_temp;
            vector < _Int > ja_lt_temp;
            vector < char >ja_char_lt_temp;
            vector < _Flt > a_lt_temp;
            ia_lt.swap (ia_lt_temp);
            ja_lt.swap (ja_lt_temp);
            ja_char_lt.swap (ja_char_lt_temp);
            a_lt.swap (a_lt_temp);
            vector < _Int > ia_u_temp;
            vector < _Int > ja_u_temp;
            vector < char >ja_char_u_temp;
            vector < _Flt > a_u_temp;
            ia_u.swap (ia_u_temp);
            ja_u.swap (ja_u_temp);
            ja_char_u.swap (ja_char_u_temp);
            a_u.swap (a_u_temp);
         }

         vector < _Int > ia_alu_ord;
         vector < _Int > ja_alu_ord;
         vector < char >ja_char_alu_ord;
         vector < _Flt > a_alu_ord;

         CFct_bxb_impl < _Int, _Flt >::ReorderMatrix (true, false, _b_size, _n, order,
                                                      ia_alu, ja_alu, ja_char_alu, a_alu,
                                                      ia_alu_ord, ja_alu_ord,
                                                      ja_char_alu_ord, a_alu_ord);

         {
            vector < _Int > ia_alu_temp;
            vector < _Int > ja_alu_temp;
            vector < char >ja_char_alu_temp;
            vector < _Flt > a_alu_temp;
            ia_alu.swap (ia_alu_temp);
            ja_alu.swap (ja_alu_temp);
            ja_char_alu.swap (ja_char_alu_temp);
            a_alu.swap (a_alu_temp);
         }

         CFct_bxb_impl < _Int, _Flt >::SplitLU (true, _b_size, _n, ia_alu_ord, ja_alu_ord,
                                                ja_char_alu_ord, a_alu_ord, ia_lt, ja_lt,
                                                ja_char_lt, a_lt, ia_u, ja_u, ja_char_u,
                                                a_u);

         {
            vector < _Int > ia_alu_temp;
            vector < _Int > ja_alu_temp;
            vector < char >ja_char_alu_temp;
            vector < _Flt > a_alu_temp;
            ia_alu_ord.swap (ia_alu_temp);
            ja_alu_ord.swap (ja_alu_temp);
            ja_char_alu_ord.swap (ja_char_alu_temp);
            a_alu_ord.swap (a_alu_temp);
         }

         CFct_bxb_impl < _Int, _Flt >::Transpose (true, _b_size, _n, ia_lt, ja_lt,
                                                  ja_char_lt, a_lt, ia_l, ja_l, ja_char_l,
                                                  a_l);

         {
            vector < _Int > ia_l_temp;
            vector < _Int > ja_l_temp;
            vector < char >ja_char_l_temp;
            vector < _Flt > a_l_temp;
            ia_lt.swap (ia_l_temp);
            ja_lt.swap (ja_l_temp);
            ja_char_lt.swap (ja_char_l_temp);
            a_lt.swap (a_l_temp);
         }

         CFct_bxb_impl < _Int, _Flt >::CombineLUPairs (true, _b_size, _n, ia_l, ja_l,
                                                       ja_char_l, a_l, ia_u, ja_u,
                                                       ja_char_u, a_u, ia_alu_ord_pair,
                                                       ja_alu_ord_pair,
                                                       ja_char_alu_ord_pair,
                                                       a_alu_ord_pair);

      }

      _Int *pia_alu_ord_pair = NULL;
      _Int *pja_alu_ord_pair = NULL;
      char *pja_char_alu_ord_pair = NULL;
      _Flt *pa_alu_ord_pair = NULL;

      if (ia_alu_ord_pair.size () > 0)
         pia_alu_ord_pair = &ia_alu_ord_pair[0];
      if (ja_alu_ord_pair.size () > 0)
         pja_alu_ord_pair = &ja_alu_ord_pair[0];
      if (ja_char_alu_ord_pair.size () > 0)
         pja_char_alu_ord_pair = &ja_char_alu_ord_pair[0];
      if (a_alu_ord_pair.size () > 0)
         pa_alu_ord_pair = &a_alu_ord_pair[0];

// Split into main and schur parts

      int n_last = _n - _n_ini;

      int nzja_ini = (int) pia_alu_ord_pair[_n_ini];
      int nzja_last = (int) (pia_alu_ord_pair[_n] - pia_alu_ord_pair[_n_ini]);

      vector < _Int > ia_ini (_n_ini + 1);
      vector < _Int > ja_ini (nzja_ini + 1);
      vector < char >ja_char_ini (nzja_ini + 1);
      vector < _Flt > a_ini (nzja_ini * b_2_2 + 1);

      _Int *pia_ini = NULL;
      _Int *pja_ini = NULL;
      char *pja_char_ini = NULL;
      _Flt *pa_ini = NULL;

      if (ia_ini.size () > 0)
         pia_ini = &ia_ini[0];
      if (ja_ini.size () > 0)
         pja_ini = &ja_ini[0];
      if (ja_char_ini.size () > 0)
         pja_char_ini = &ja_char_ini[0];
      if (a_ini.size () > 0)
         pa_ini = &a_ini[0];

      for (i = 0; i <= _n_ini; i++)
         pia_ini[i] = pia_alu_ord_pair[i];
      for (i = 0; i < nzja_ini; i++)
         pja_ini[i] = pja_alu_ord_pair[i];
      for (i = 0; i < nzja_ini; i++)
         pja_char_ini[i] = pja_char_alu_ord_pair[i];
      for (i = 0; i < nzja_ini * b_2_2; i++)
         pa_ini[i] = pa_alu_ord_pair[i];

      vector < _Int > ia_last (n_last + 1);
      vector < _Int > ja_last (nzja_last + 1);
      vector < char >ja_char_last (nzja_last + 1);
      vector < _Flt > a_last (nzja_last * b_2_2 + 1);

      _Int *pia_last = NULL;
      _Int *pja_last = NULL;
      char *pja_char_last = NULL;
      _Flt *pa_last = NULL;

      if (ia_last.size () > 0)
         pia_last = &ia_last[0];
      if (ja_last.size () > 0)
         pja_last = &ja_last[0];
      if (ja_char_last.size () > 0)
         pja_char_last = &ja_char_last[0];
      if (a_last.size () > 0)
         pa_last = &a_last[0];

      for (i = 0; i <= n_last; i++)
         pia_last[i] = pia_alu_ord_pair[i + _n_ini] - pia_alu_ord_pair[_n_ini];
      for (i = 0; i < nzja_last; i++)
         pja_last[i] = pja_alu_ord_pair[nzja_ini + i] - _n_ini;
      for (i = 0; i < nzja_last; i++)
         pja_char_last[i] = pja_char_alu_ord_pair[nzja_ini + i];
      int ishift = nzja_ini * b_2_2;
      for (i = 0; i < nzja_last * b_2_2; i++)
         pa_last[i] = pa_alu_ord_pair[ishift + i];

// Fct

      vector < _Int > ia_lu_ini;
      vector < _Int > ja_lu_ini;
      vector < char >ja_lu_char_ini;
      vector < _Flt > a_lu_ini;

      CFct_bxb_impl < _Int, _Flt >::FctBlockIlu2Degree (_params, _b_size, _n, _n_ini,
                                                        ia_ini, ja_ini, ja_char_ini,
                                                        a_ini, ia_lu_ini, ja_lu_ini,
                                                        ja_lu_char_ini, a_lu_ini, _nmodif,
                                                        _eigmin_att, _eigmax_att);

      _Int *pia_lu_ini = NULL;
      _Int *pja_lu_ini = NULL;
      char *pja_lu_char_ini = NULL;
      _Flt *pa_lu_ini = NULL;

      if (ia_lu_ini.size () > 0)
         pia_lu_ini = &ia_lu_ini[0];
      if (ja_lu_ini.size () > 0)
         pja_lu_ini = &ja_lu_ini[0];
      if (ja_lu_char_ini.size () > 0)
         pja_lu_char_ini = &ja_lu_char_ini[0];
      if (a_lu_ini.size () > 0)
         pa_lu_ini = &a_lu_ini[0];

// Finalize U and Schur complement

      int nzja_u = (int) ia_lu_ini[_n_ini];
      int nzja_u2 = (int) (ia_lu_ini[_n] - ia_lu_ini[_n_ini]);

      _ia_u.resize (_n_ini + 1);
      _ja_u.resize (nzja_u + 1);
      _ja_char_u.resize (nzja_u + 1);
      _a_u.resize (nzja_u * b_2_2 + 1);

      _Int *p_ia_u = NULL;
      _Int *p_ja_u = NULL;
      char *p_ja_char_u = NULL;
      _Flt *p_a_u = NULL;

      if (_ia_u.size () > 0)
         p_ia_u = &_ia_u[0];
      if (_ja_u.size () > 0)
         p_ja_u = &_ja_u[0];
      if (_ja_char_u.size () > 0)
         p_ja_char_u = &_ja_char_u[0];
      if (_a_u.size () > 0)
         p_a_u = &_a_u[0];

      for (i = 0; i <= _n_ini; i++)
         p_ia_u[i] = pia_lu_ini[i];
      for (i = 0; i < nzja_u; i++)
         p_ja_u[i] = pja_lu_ini[i];
      for (i = 0; i < nzja_u; i++)
         p_ja_char_u[i] = pja_lu_char_ini[i];
      for (i = 0; i < nzja_u * b_2_2; i++)
         p_a_u[i] = pa_lu_ini[i];

      vector < _Int > ia_u2 (n_last + 1);
      vector < _Int > ja_u2 (nzja_u2 + 1);
      vector < char >ja_char_u2 (nzja_u2 + 1);
      vector < _Flt > a_u2 (nzja_u2 * b_2_2 + 1);

      _Int *pia_u2 = NULL;
      _Int *pja_u2 = NULL;
      char *pja_char_u2 = NULL;
      _Flt *pa_u2 = NULL;

      if (ia_u2.size () > 0)
         pia_u2 = &ia_u2[0];
      if (ja_u2.size () > 0)
         pja_u2 = &ja_u2[0];
      if (ja_char_u2.size () > 0)
         pja_char_u2 = &ja_char_u2[0];
      if (a_u2.size () > 0)
         pa_u2 = &a_u2[0];

      for (i = 0; i <= n_last; i++)
         pia_u2[i] = pia_lu_ini[i + _n_ini] - pia_lu_ini[_n_ini];
      for (i = 0; i < nzja_u2; i++)
         pja_u2[i] = pja_lu_ini[nzja_u + i] - _n_ini;
      for (i = 0; i < nzja_u2; i++)
         pja_char_u2[i] = pja_lu_char_ini[nzja_u + i];
      ishift = nzja_u * b_2_2;
      for (i = 0; i < nzja_u2 * b_2_2; i++)
         pa_u2[i] = pa_lu_ini[ishift + i];

      CFct_bxb_impl < _Int, _Flt >::AddMatrices (true, true, _b_size, n_last, ia_last,
                                                 ja_last, ja_char_last, a_last, ia_u2,
                                                 ja_u2, ja_char_u2, a_u2, _ia_u_schur,
                                                 _ja_u_schur, _ja_char_u_schur,
                                                 _a_u_schur);

   }

//
// Perform one step fct for dinamically splitted and ordered triangular matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::FctMLevDiaSplit (SParams * _params, int _nlev, double *_dia_lev,
                               int _b_size, int _n, const vector < _Int > &_ia_au,
                               const vector < _Int > &_ja_au,
                               const vector < char >&_ja_char_au,
                               const vector < _Flt > &_a_au, int *_order, int &_nblks_fct,
                               int *_blks_fct, vector < _Int > &_ia_u,
                               vector < _Int > &_ja_u, vector < char >&_ja_char_u,
                               vector < _Flt > &_a_u, int &_nmodif, double &_eigmin_att,
                               double &_eigmax_att)
   {

// Open matrix data structures

      int b_2 = _b_size * _b_size;
      int b_2_2 = 2 * b_2;

      const _Int *p_ia_au = NULL;
      const _Int *p_ja_au = NULL;
      const char *p_ja_char_au = NULL;
      const _Flt *p_a_au = NULL;

      if (_ia_au.size () > 0)
         p_ia_au = &_ia_au[0];
      if (_ja_au.size () > 0)
         p_ja_au = &_ja_au[0];
      if (_ja_char_au.size () > 0)
         p_ja_char_au = &_ja_char_au[0];
      if (_a_au.size () > 0)
         p_a_au = &_a_au[0];

// Allocate memory to store level data

      CVectorData < vector < _Int > >ia_U_lev (_nlev + 1);
      CVectorData < vector < _Int > >ja_U_lev (_nlev + 1);
      CVectorData < vector < char > >ja_char_U_lev (_nlev + 1);
      CVectorData < vector < _Flt > >a_U_lev (_nlev + 1);

      vector < _Int > *pia_U_lev = ia_U_lev.Ptr ();
      vector < _Int > *pja_U_lev = ja_U_lev.Ptr ();
      vector < char >*pja_char_U_lev = ja_char_U_lev.Ptr ();
      vector < _Flt > *pa_U_lev = a_U_lev.Ptr ();

      CVectorData < vector < _Int > >ia_USchur_lev (_nlev + 1);
      CVectorData < vector < _Int > >ja_USchur_lev (_nlev + 1);
      CVectorData < vector < char > >ja_char_USchur_lev (_nlev + 1);
      CVectorData < vector < _Flt > >a_USchur_lev (_nlev + 1);

      vector < _Int > *pia_USchur_lev = ia_USchur_lev.Ptr ();
      vector < _Int > *pja_USchur_lev = ja_USchur_lev.Ptr ();
      vector < char >*pja_char_USchur_lev = ja_char_USchur_lev.Ptr ();
      vector < _Flt > *pa_USchur_lev = a_USchur_lev.Ptr ();

      CVectorData < int >n_lev (_nlev + 2);
      CVectorData < int >n_ini_lev (_nlev + 1);
      CVectorData < vector < int > >order_lev (_nlev + 1);
      CVectorData < int >nmodif_lev (_nlev + 1);
      CVectorData < double >eigmin_lev (_nlev + 1);
      CVectorData < double >eigmax_lev (_nlev + 1);

      int *pn_lev = n_lev.Ptr ();
      int *pn_ini_lev = n_ini_lev.Ptr ();
      vector < int >*porder_lev = order_lev.Ptr ();
      int *pnmodif_lev = nmodif_lev.Ptr ();
      double *peigmin_lev = eigmin_lev.Ptr ();
      double *peigmax_lev = eigmax_lev.Ptr ();

// Perform level computations

      int ilev, n_loc, n_ini_loc;
      double dia_loc;

      pn_lev[0] = _n;

      for (ilev = 0; ilev < _nlev; ilev++) {

         n_loc = pn_lev[ilev];
         n_ini_loc = 0;

         dia_loc = _dia_lev[ilev];
         if (ilev == _nlev - 1)
            dia_loc = -1.0e0;

         pn_ini_lev[ilev] = 0;
         pnmodif_lev[ilev] = 0;
         peigmin_lev[ilev] = 1.0e100;
         peigmax_lev[ilev] = -1.0e100;

         porder_lev[ilev].resize (n_loc + 1);

         int *pporder_lev = &porder_lev[ilev][0];

         if (n_loc > 0) {

            if (ilev == 0) {
               CFct_bxb_impl < _Int, _Flt >::DiaSplitFctSchur (_params, dia_loc, _b_size,
                                                               n_loc, _ia_au, _ja_au,
                                                               _ja_char_au, _a_au,
                                                               pporder_lev,
                                                               pn_ini_lev[ilev],
                                                               pia_U_lev[ilev],
                                                               pja_U_lev[ilev],
                                                               pja_char_U_lev[ilev],
                                                               pa_U_lev[ilev],
                                                               pia_USchur_lev[ilev],
                                                               pja_USchur_lev[ilev],
                                                               pja_char_USchur_lev[ilev],
                                                               pa_USchur_lev[ilev],
                                                               pnmodif_lev[ilev],
                                                               peigmin_lev[ilev],
                                                               peigmax_lev[ilev]);
            } else {
               CFct_bxb_impl < _Int, _Flt >::DiaSplitFctSchur (_params, dia_loc, _b_size,
                                                               n_loc,
                                                               pia_USchur_lev[ilev - 1],
                                                               pja_USchur_lev[ilev - 1],
                                                               pja_char_USchur_lev[ilev -
                                                                                   1],
                                                               pa_USchur_lev[ilev - 1],
                                                               pporder_lev,
                                                               pn_ini_lev[ilev],
                                                               pia_U_lev[ilev],
                                                               pja_U_lev[ilev],
                                                               pja_char_U_lev[ilev],
                                                               pa_U_lev[ilev],
                                                               pia_USchur_lev[ilev],
                                                               pja_USchur_lev[ilev],
                                                               pja_char_USchur_lev[ilev],
                                                               pa_USchur_lev[ilev],
                                                               pnmodif_lev[ilev],
                                                               peigmin_lev[ilev],
                                                               peigmax_lev[ilev]);
               {
                  vector < _Int > ia_temp;
                  vector < _Int > ja_temp;
                  vector < char >ja_char_temp;
                  vector < _Flt > a_temp;
                  pia_USchur_lev[ilev - 1].swap (ia_temp);
                  pja_USchur_lev[ilev - 1].swap (ja_temp);
                  pja_char_USchur_lev[ilev - 1].swap (ja_char_temp);
                  pa_USchur_lev[ilev - 1].swap (a_temp);
               }
            }

            n_ini_loc = pn_ini_lev[ilev];

         }

         pn_lev[ilev + 1] = n_loc - n_ini_loc;

      }

// Backward cycle for ordering collection and reordering of fct bordering

      int n_last, i, iold;

      for (ilev = _nlev - 2; ilev >= 0; ilev--) {

// Compute extended ordering of bordering

         n_loc = pn_lev[ilev];
         n_ini_loc = pn_ini_lev[ilev];
         n_last = n_loc - n_ini_loc;

         vector < int >order_brd (n_loc + 1);
         int *porder_brd = &order_brd[0];

         for (i = 0; i < n_ini_loc; i++)
            porder_brd[i] = i;

         int *porder_last = &porder_lev[ilev + 1][0];

         for (i = 0; i < n_last; i++)
            porder_brd[n_ini_loc + i] = n_ini_loc + porder_last[i];

// Reorder and move

         vector < _Int > ia_ord;
         vector < _Int > ja_ord;
         vector < char >ja_char_ord;
         vector < _Flt > a_char_ord;

         CFct_bxb_impl < _Int, _Flt >::ReorderMatrix (true, true, _b_size, n_ini_loc,
                                                      order_brd, pia_U_lev[ilev],
                                                      pja_U_lev[ilev],
                                                      pja_char_U_lev[ilev],
                                                      pa_U_lev[ilev], ia_ord, ja_ord,
                                                      ja_char_ord, a_char_ord);

         pia_U_lev[ilev].swap (ia_ord);
         pja_U_lev[ilev].swap (ja_ord);
         pja_char_U_lev[ilev].swap (ja_char_ord);
         pa_U_lev[ilev].swap (a_char_ord);

// Update ordering

         int *pporder_lev = &porder_lev[ilev][0];

         vector < int >iorder (n_loc + 1);
         int *piorder = &iorder[0];

         for (i = 0; i < n_loc; i++)
            piorder[pporder_lev[i]] = i;

         for (i = 0; i < n_loc; i++) {
            iold = piorder[i];
            pporder_lev[iold] = porder_brd[i];
         }

      }

// Store final ordering

      int *pporder_lev = &porder_lev[0][0];

      for (i = 0; i < _n; i++)
         _order[i] = pporder_lev[i];

// Combine all computed fct data into one

      _nblks_fct = 0;
      _blks_fct[0] = 0;

      for (ilev = 0; ilev < _nlev; ilev++) {
         _blks_fct[_nblks_fct + 1] = _blks_fct[_nblks_fct] + pn_ini_lev[ilev];
         _nblks_fct++;
      }

      int nzja_tot = 0;

      for (ilev = 0; ilev < _nlev; ilev++) {
         n_ini_loc = pn_ini_lev[ilev];
         if (n_ini_loc > 0) {
            nzja_tot += pia_U_lev[ilev][n_ini_loc];
         }
      }

      _ia_u.resize (_n + 1);
      _ja_u.resize (nzja_tot + 1);
      _ja_char_u.resize (nzja_tot + 1);
      _a_u.resize (nzja_tot * b_2_2 + 1);

      _Int *p_ia_u = NULL;
      _Int *p_ja_u = NULL;
      char *p_ja_char_u = NULL;
      _Flt *p_a_u = NULL;

      if (_ia_u.size () > 0)
         p_ia_u = &_ia_u[0];
      if (_ja_u.size () > 0)
         p_ja_u = &_ja_u[0];
      if (_ja_char_u.size () > 0)
         p_ja_char_u = &_ja_char_u[0];
      if (_a_u.size () > 0)
         p_a_u = &_a_u[0];

      nzja_tot = 0;
      p_ia_u[0] = 0;

      int n_start, nzja_loc;

      for (ilev = 0; ilev < _nlev; ilev++) {

         n_start = _blks_fct[ilev];
         n_ini_loc = pn_ini_lev[ilev];

         nzja_loc = 0;

         if (n_ini_loc > 0) {
            nzja_loc = pia_U_lev[ilev][n_ini_loc];
         }

         _Int *pia_u_temp = NULL;
         _Int *pja_u_temp = NULL;
         char *pja_char_u_temp = NULL;
         _Flt *pa_u_temp = NULL;

         if (pia_U_lev[ilev].size () > 0)
            pia_u_temp = &pia_U_lev[ilev][0];
         if (pja_U_lev[ilev].size () > 0)
            pja_u_temp = &pja_U_lev[ilev][0];
         if (pja_char_U_lev[ilev].size () > 0)
            pja_char_u_temp = &pja_char_U_lev[ilev][0];
         if (pa_U_lev[ilev].size () > 0)
            pa_u_temp = &pa_U_lev[ilev][0];

         for (i = 0; i < n_ini_loc; i++)
            p_ia_u[n_start + i + 1] =
               p_ia_u[n_start + i] + (pia_u_temp[i + 1] - pia_u_temp[i]);
         for (i = 0; i < nzja_loc; i++)
            p_ja_u[nzja_tot + i] = n_start + pja_u_temp[i];
         for (i = 0; i < nzja_loc; i++)
            p_ja_char_u[nzja_tot + i] = pja_char_u_temp[i];
         CVector < _Flt >::CopyVector (nzja_loc * b_2_2, pa_u_temp,
                                       p_a_u + nzja_tot * b_2_2);

         nzja_tot += nzja_loc;

      }

// Filter blocks partitioning

      _nblks_fct = 0;

      for (ilev = 0; ilev < _nlev; ilev++) {
         if (pn_ini_lev[ilev] > 0) {
            _blks_fct[_nblks_fct + 1] = _blks_fct[_nblks_fct] + pn_ini_lev[ilev];
            _nblks_fct++;
         }
      }

// Get fct statistics

      _nmodif = 0;
      _eigmin_att = 1.0e100;
      _eigmax_att = -1.0e100;

      for (ilev = 0; ilev < _nlev; ilev++) {
         _nmodif += pnmodif_lev[ilev];
         if (peigmin_lev[ilev] < _eigmin_att)
            _eigmin_att = peigmin_lev[ilev];
         if (peigmax_lev[ilev] > _eigmax_att)
            _eigmax_att = peigmax_lev[ilev];
      }

   }

//
// Perform one step fct for dinamically splitted and ordered triangular matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb_impl < _Int,
      _Flt >::FctMLevDiaSplit (SParams * _params, int _nlev, double *_dia_lev,
                               int _b_size, int _n, int _n_ini,
                               const vector < _Int > &_ia_au,
                               const vector < _Int > &_ja_au,
                               const vector < char >&_ja_char_au,
                               const vector < _Flt > &_a_au, int *_order,
                               int &_n_moved_diag, int &_nblks_fct, int *_blks_fct,
                               vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                               vector < char >&_ja_char_u, vector < _Flt > &_a_u,
                               int &_nmodif, double &_eigmin_att, double &_eigmax_att)
   {

// Open matrix data structures

//      ofstream *pfout_debug = _params->pfout;

      int b_2 = _b_size * _b_size;
      int b_2_2 = 2 * b_2;

      const _Int *p_ia_au = NULL;
      const _Int *p_ja_au = NULL;
      const char *p_ja_char_au = NULL;
      const _Flt *p_a_au = NULL;

      if (_ia_au.size () > 0)
         p_ia_au = &_ia_au[0];
      if (_ja_au.size () > 0)
         p_ja_au = &_ja_au[0];
      if (_ja_char_au.size () > 0)
         p_ja_char_au = &_ja_char_au[0];
      if (_a_au.size () > 0)
         p_a_au = &_a_au[0];

// Allocate memory to store level data

      CVectorData < vector < _Int > >ia_U_lev (_nlev + 1);
      CVectorData < vector < _Int > >ja_U_lev (_nlev + 1);
      CVectorData < vector < char > >ja_char_U_lev (_nlev + 1);
      CVectorData < vector < _Flt > >a_U_lev (_nlev + 1);

      vector < _Int > *pia_U_lev = ia_U_lev.Ptr ();
      vector < _Int > *pja_U_lev = ja_U_lev.Ptr ();
      vector < char >*pja_char_U_lev = ja_char_U_lev.Ptr ();
      vector < _Flt > *pa_U_lev = a_U_lev.Ptr ();

      CVectorData < vector < _Int > >ia_USchur_lev (_nlev + 1);
      CVectorData < vector < _Int > >ja_USchur_lev (_nlev + 1);
      CVectorData < vector < char > >ja_char_USchur_lev (_nlev + 1);
      CVectorData < vector < _Flt > >a_USchur_lev (_nlev + 1);

      vector < _Int > *pia_USchur_lev = ia_USchur_lev.Ptr ();
      vector < _Int > *pja_USchur_lev = ja_USchur_lev.Ptr ();
      vector < char >*pja_char_USchur_lev = ja_char_USchur_lev.Ptr ();
      vector < _Flt > *pa_USchur_lev = a_USchur_lev.Ptr ();

      CVectorData < int >n_lev (_nlev + 2);
      CVectorData < int >n_ini_max_lev (_nlev + 2);
      CVectorData < int >n_ini_lev (_nlev + 1);
      CVectorData < int >n_bad_diag_lev (_nlev + 1);
      CVectorData < vector < int > >order_lev (_nlev + 1);
      CVectorData < int >nmodif_lev (_nlev + 1);
      CVectorData < double >eigmin_lev (_nlev + 1);
      CVectorData < double >eigmax_lev (_nlev + 1);

      int *pn_lev = n_lev.Ptr ();
      int *pn_ini_max_lev = n_ini_max_lev.Ptr ();
      int *pn_ini_lev = n_ini_lev.Ptr ();
      int *pn_bad_diag_lev = n_bad_diag_lev.Ptr ();
      vector < int >*porder_lev = order_lev.Ptr ();
      int *pnmodif_lev = nmodif_lev.Ptr ();
      double *peigmin_lev = eigmin_lev.Ptr ();
      double *peigmax_lev = eigmax_lev.Ptr ();

// Perform level computations

      int ilev, n_loc, n_ini_loc, n_ini_max_loc;
      double dia_loc;

      pn_lev[0] = _n;
      pn_ini_max_lev[0] = _n_ini;

      int nlev_curr = 0;

      for (ilev = 0; ilev < _nlev; ilev++) {

         n_loc = pn_lev[ilev];
         n_ini_max_loc = pn_ini_max_lev[ilev];
         n_ini_loc = 0;

         dia_loc = _dia_lev[ilev];
         if (ilev == _nlev - 1)
            dia_loc = -1.0e0;

         pn_ini_lev[ilev] = 0;
         pnmodif_lev[ilev] = 0;
         peigmin_lev[ilev] = 1.0e100;
         peigmax_lev[ilev] = -1.0e100;

         porder_lev[ilev].resize (n_loc + 1);

         int *pporder_lev = &porder_lev[ilev][0];

         if (n_ini_max_loc > 0) {

            if (ilev == 0) {
               CFct_bxb_impl < _Int, _Flt >::DiaSplitFctSchur (_params, dia_loc, _b_size,
                                                               n_ini_max_loc, n_loc,
                                                               n_ini_max_loc, _ia_au,
                                                               _ja_au, _ja_char_au, _a_au,
                                                               pporder_lev,
                                                               pn_ini_lev[ilev],
                                                               pn_bad_diag_lev[ilev],
                                                               pia_U_lev[ilev],
                                                               pja_U_lev[ilev],
                                                               pja_char_U_lev[ilev],
                                                               pa_U_lev[ilev],
                                                               pia_USchur_lev[ilev],
                                                               pja_USchur_lev[ilev],
                                                               pja_char_USchur_lev[ilev],
                                                               pa_USchur_lev[ilev],
                                                               pnmodif_lev[ilev],
                                                               peigmin_lev[ilev],
                                                               peigmax_lev[ilev]);
            } else {
               CFct_bxb_impl < _Int, _Flt >::DiaSplitFctSchur (_params, dia_loc, _b_size,
                                                               n_loc, n_loc,
                                                               n_ini_max_loc,
                                                               pia_USchur_lev[ilev - 1],
                                                               pja_USchur_lev[ilev - 1],
                                                               pja_char_USchur_lev[ilev -
                                                                                   1],
                                                               pa_USchur_lev[ilev - 1],
                                                               pporder_lev,
                                                               pn_ini_lev[ilev],
                                                               pn_bad_diag_lev[ilev],
                                                               pia_U_lev[ilev],
                                                               pja_U_lev[ilev],
                                                               pja_char_U_lev[ilev],
                                                               pa_U_lev[ilev],
                                                               pia_USchur_lev[ilev],
                                                               pja_USchur_lev[ilev],
                                                               pja_char_USchur_lev[ilev],
                                                               pa_USchur_lev[ilev],
                                                               pnmodif_lev[ilev],
                                                               peigmin_lev[ilev],
                                                               peigmax_lev[ilev]);
               {
                  vector < _Int > ia_temp;
                  vector < _Int > ja_temp;
                  vector < char >ja_char_temp;
                  vector < _Flt > a_temp;
                  pia_USchur_lev[ilev - 1].swap (ia_temp);
                  pja_USchur_lev[ilev - 1].swap (ja_temp);
                  pja_char_USchur_lev[ilev - 1].swap (ja_char_temp);
                  pa_USchur_lev[ilev - 1].swap (a_temp);
               }
            }

         } else {

            if (n_loc != 0 || n_ini_loc != 0) {
               cout <<
                  " CFct_bxb_impl < _Int,_Flt >::FctMLevDiaSplit: error in empty level step ! "
                  << endl;
               cout << " N = " << _n << " Nini = " << _n_ini << " Ilev = " << ilev <<
                  " N_loc " << n_loc << " N_ini_loc = " << n_ini_loc << endl;
               throw
                  " CFct_bxb_impl < _Int,_Flt >::FctMLevDiaSplit: error in empty level step ! ";
            }

            pn_ini_lev[ilev] = 0;
            pn_bad_diag_lev[ilev] = 0;
            pnmodif_lev[ilev] = 0;
            peigmin_lev[ilev] = 1.0e100;
            peigmax_lev[ilev] = -1.0e100;

            pia_U_lev[ilev].resize (2);
            pja_U_lev[ilev].resize (1);
            pja_char_U_lev[ilev].resize (1);
            pa_U_lev[ilev].resize (1);

            pia_U_lev[ilev][0] = 0;

            pia_USchur_lev[ilev].resize (2);
            pja_USchur_lev[ilev].resize (1);
            pja_char_USchur_lev[ilev].resize (1);
            pa_USchur_lev[ilev].resize (1);

            pia_USchur_lev[ilev][0] = 0;

         }

         n_ini_loc = pn_ini_lev[ilev];

         pn_lev[ilev + 1] = n_loc - n_ini_loc;
         pn_ini_max_lev[ilev + 1] = n_ini_max_loc - n_ini_loc;

         nlev_curr++;

         if (n_ini_loc == n_ini_max_loc)
            break;

      }

// Backward cycle for ordering collection and reordering of fct bordering

      int n_last, i, iold;

      for (ilev = nlev_curr - 2; ilev >= 0; ilev--) {

// Compute extended ordering of bordering

         n_loc = pn_lev[ilev];
         n_ini_loc = pn_ini_lev[ilev];
         n_last = n_loc - n_ini_loc;

         vector < int >order_brd (n_loc + 1);
         int *porder_brd = &order_brd[0];

         for (i = 0; i < n_ini_loc; i++)
            porder_brd[i] = i;

         int *porder_last = &porder_lev[ilev + 1][0];

         for (i = 0; i < n_last; i++)
            porder_brd[n_ini_loc + i] = n_ini_loc + porder_last[i];

// Reorder and move

         vector < _Int > ia_ord;
         vector < _Int > ja_ord;
         vector < char >ja_char_ord;
         vector < _Flt > a_ord;

         CFct_bxb_impl < _Int, _Flt >::ReorderMatrix (true, true, _b_size, n_ini_loc,
                                                      order_brd, pia_U_lev[ilev],
                                                      pja_U_lev[ilev],
                                                      pja_char_U_lev[ilev],
                                                      pa_U_lev[ilev], ia_ord, ja_ord,
                                                      ja_char_ord, a_ord);

         pia_U_lev[ilev].swap (ia_ord);
         pja_U_lev[ilev].swap (ja_ord);
         pja_char_U_lev[ilev].swap (ja_char_ord);
         pa_U_lev[ilev].swap (a_ord);

// Update ordering

         int *pporder_lev = &porder_lev[ilev][0];

         vector < int >iorder (n_loc + 1);
         int *piorder = &iorder[0];

         for (i = 0; i < n_loc; i++)
            piorder[pporder_lev[i]] = i;

         for (i = 0; i < n_loc; i++) {
            iold = piorder[i];
            pporder_lev[iold] = porder_brd[i];
         }

      }

// Store final ordering

      int *pporder_lev = &porder_lev[0][0];

      for (i = 0; i < _n; i++)
         _order[i] = pporder_lev[i];

// Combine all computed fct data and last Schur into one

      _nblks_fct = 0;
      _blks_fct[0] = 0;

      for (ilev = 0; ilev < nlev_curr; ilev++) {
         _blks_fct[_nblks_fct + 1] = _blks_fct[_nblks_fct] + pn_ini_lev[ilev];
         _nblks_fct++;
      }

      if (_blks_fct[_nblks_fct] != _n) {
         _blks_fct[_nblks_fct + 1] = _n;
         _nblks_fct++;
      }

      int nzja_tot = 0;

      for (ilev = 0; ilev < nlev_curr; ilev++) {
         n_ini_loc = pn_ini_lev[ilev];
         if (n_ini_loc > 0) {
            nzja_tot += (int) pia_U_lev[ilev][n_ini_loc];
         }
      }

      {
         ilev = nlev_curr - 1;
         n_ini_loc = pn_lev[ilev] - pn_ini_lev[ilev];
         if (n_ini_loc > 0) {
            nzja_tot += (int) pia_USchur_lev[ilev][n_ini_loc];
         }
      }

      _ia_u.resize (_n + 1);
      _ja_u.resize (nzja_tot + 1);
      _ja_char_u.resize (nzja_tot + 1);
      _a_u.resize (nzja_tot * b_2_2 + 1);

      _Int *p_ia_u = NULL;
      _Int *p_ja_u = NULL;
      char *p_ja_char_u = NULL;
      _Flt *p_a_u = NULL;

      if (_ia_u.size () > 0)
         p_ia_u = &_ia_u[0];
      if (_ja_u.size () > 0)
         p_ja_u = &_ja_u[0];
      if (_ja_char_u.size () > 0)
         p_ja_char_u = &_ja_char_u[0];
      if (_a_u.size () > 0)
         p_a_u = &_a_u[0];

      nzja_tot = 0;
      p_ia_u[0] = 0;

      int n_start, nzja_loc;

      for (ilev = 0; ilev < nlev_curr; ilev++) {

         n_start = _blks_fct[ilev];
         n_ini_loc = pn_ini_lev[ilev];

         nzja_loc = 0;

         if (n_ini_loc > 0) {
            nzja_loc = (int) pia_U_lev[ilev][n_ini_loc];
         }

         _Int *pia_u_temp = NULL;
         _Int *pja_u_temp = NULL;
         char *pja_char_u_temp = NULL;
         _Flt *pa_u_temp = NULL;

         if (pia_U_lev[ilev].size () > 0)
            pia_u_temp = &pia_U_lev[ilev][0];
         if (pja_U_lev[ilev].size () > 0)
            pja_u_temp = &pja_U_lev[ilev][0];
         if (pja_char_U_lev[ilev].size () > 0)
            pja_char_u_temp = &pja_char_U_lev[ilev][0];
         if (pa_U_lev[ilev].size () > 0)
            pa_u_temp = &pa_U_lev[ilev][0];

         for (i = 0; i < n_ini_loc; i++)
            p_ia_u[n_start + i + 1] =
               p_ia_u[n_start + i] + (pia_u_temp[i + 1] - pia_u_temp[i]);
         for (i = 0; i < nzja_loc; i++)
            p_ja_u[nzja_tot + i] = n_start + pja_u_temp[i];
         for (i = 0; i < nzja_loc; i++)
            p_ja_char_u[nzja_tot + i] = pja_char_u_temp[i];
         CVector < _Flt >::CopyVector (nzja_loc * b_2_2, pa_u_temp,
                                       p_a_u + nzja_tot * b_2_2);

         nzja_tot += nzja_loc;

      }

      {
         ilev = nlev_curr - 1;
         n_start = _blks_fct[ilev + 1];
         n_ini_loc = pn_lev[ilev] - pn_ini_lev[ilev];
         if (n_start != _n) {

            _Int *pia_u_temp = pia_USchur_lev[ilev].data ();
            _Int *pja_u_temp = pja_USchur_lev[ilev].data ();
            char *pja_char_u_temp = pja_char_USchur_lev[ilev].data ();
            _Flt *pa_u_temp = pa_USchur_lev[ilev].data ();

            nzja_loc = 0;
            if (n_ini_loc > 0) {
               nzja_loc = (int) pia_USchur_lev[ilev][n_ini_loc];
            }

            for (i = 0; i < n_ini_loc; i++)
               p_ia_u[n_start + i + 1] =
                  p_ia_u[n_start + i] + (pia_u_temp[i + 1] - pia_u_temp[i]);
            for (i = 0; i < nzja_loc; i++)
               p_ja_u[nzja_tot + i] = n_start + pja_u_temp[i];
            for (i = 0; i < nzja_loc; i++)
               p_ja_char_u[nzja_tot + i] = pja_char_u_temp[i];
            CVector < _Flt >::CopyVector (nzja_loc * b_2_2, pa_u_temp,
                                          p_a_u + nzja_tot * b_2_2);
         }
      }

// Get fct statistics

      _n_moved_diag = 0;
      _nmodif = 0;
      _eigmin_att = 1.0e100;
      _eigmax_att = -1.0e100;

      for (ilev = 0; ilev < nlev_curr; ilev++) {
         _n_moved_diag += pn_bad_diag_lev[ilev];
         _nmodif += pnmodif_lev[ilev];
         if (peigmin_lev[ilev] < _eigmin_att)
            _eigmin_att = peigmin_lev[ilev];
         if (peigmax_lev[ilev] > _eigmax_att)
            _eigmax_att = peigmax_lev[ilev];
      }

   }

// Perform multiplication of matrices stored by rows, resulting matrix is also stored by rows
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb < _Int,
      _Flt >::MatrixByMatrixMultiply_BxB (int _blksize, int &_icycle, int *_imask,
                                          int *_imask1, int *_indarr, int *_listloc,
                                          _Flt * _fmask, const CMatrix < _Int, _Flt > &_a,
                                          const CMatrix < _Int, _Flt > &_b,
                                          CMatrix < _Int, _Flt > &_a_times_b)
   {

      int b_2 = _blksize * _blksize;

// Open sparsities of A and B

      int nlist_a = _a.GetNlist ();
      const _Int *plist_a = _a.GetListArr ();
      const _Int *pia_a = _a.GetIaArr ();
      const _Int *pja_a = _a.GetJaArr ();
      const _Flt *pa_a = _a.GetAArr ();

      int nlist_b = _b.GetNlist ();
      const _Int *plist_b = _b.GetListArr ();
      const _Int *pia_b = _b.GetIaArr ();
      const _Int *pja_b = _b.GetJaArr ();
      const _Flt *pa_b = _b.GetAArr ();

// Perform multiplication

      vector < _Int > list_ab (nlist_a + 1);
      vector < _Int > ia_ab (nlist_a + 1);

      _Int *plist_ab = &list_ab[0];
      _Int *pia_ab = &ia_ab[0];

      vector < _Int > ja_ab (1);
      vector < _Flt > a_ab (1);

      int i;

      for (i = 0; i < nlist_a; i++)
         plist_ab[i] = plist_a[i];

      pia_ab[0] = 0;

      _icycle++;

      int icycle1 = _icycle;

      int jj;

      for (i = 0; i < nlist_b; i++) {
         jj = (int) plist_b[i];
         _imask1[jj] = icycle1;
         _indarr[jj] = i;
      }

      _Flt *p_a_ab_temp;

      int nzja_ab = 0;

      int j, k, kk, nlistloc, ind;

      for (i = 0; i < nlist_a; i++) {
         _icycle++;
         nlistloc = 0;
         for (j = (int) pia_a[i]; j < pia_a[i + 1]; j++) {
            jj = (int) pja_a[j];
            if (_imask1[jj] == icycle1) {
               ind = _indarr[jj];
               for (k = (int) pia_b[ind]; k < pia_b[ind + 1]; k++) {
                  kk = (int) pja_b[k];
                  if (_imask[kk] != _icycle) {
                     _listloc[nlistloc] = (int) kk;
                     nlistloc++;
                     _imask[kk] = _icycle;
                  }
               }
            }
         }
         sort (_listloc, _listloc + nlistloc);
         for (j = 0; j < nlistloc; j++) {
            jj = _listloc[j];
            CVector < _Flt >::SetByZeroes (b_2, _fmask + jj * b_2);
         }
         for (j = (int) pia_a[i]; j < pia_a[i + 1]; j++) {
            jj = (int) pja_a[j];
            if (_imask1[jj] == icycle1) {
               ind = _indarr[jj];
               for (k = (int) pia_b[ind]; k < pia_b[ind + 1]; k++) {
                  kk = (int) pja_b[k];
                  CBlock_BxB_traits < _Flt >::MM_Add_BxB (_blksize, pa_a + j * b_2,
                                                          pa_b + k * b_2,
                                                          _fmask + kk * b_2);
               }
            }
         }
         ja_ab.resize (nzja_ab + nlistloc + 1);
         a_ab.resize ((nzja_ab + nlistloc) * b_2 + 1);
         p_a_ab_temp = a_ab.data ();
         for (j = 0; j < nlistloc; j++)
            ja_ab[nzja_ab + j] = (_Int) _listloc[j];
         for (j = 0; j < nlistloc; j++) {
            jj = _listloc[j];
            CVector < _Flt >::CopyVector (b_2, _fmask + jj * b_2,
                                          p_a_ab_temp + (nzja_ab + j) * b_2);
         }
         nzja_ab += nlistloc;
         pia_ab[i + 1] = nzja_ab;
      }

// Store result

      vector < _Int > *p_list_ab = _a_times_b.GetList ();
      vector < _Int > *p_ia_ab = _a_times_b.GetIa ();
      vector < _Int > *p_ja_ab = _a_times_b.GetJa ();
      vector < _Flt > *p_a_ab = _a_times_b.GetA ();

      p_list_ab->swap (list_ab);
      p_ia_ab->swap (ia_ab);
      p_ja_ab->swap (ja_ab);
      p_a_ab->swap (a_ab);

      _a_times_b.SetNlist (nlist_a);
      _a_times_b.SetNzja (nzja_ab);
      _a_times_b.SetNza (nzja_ab * b_2);

   }

// Perform hmatrix by hmatrix multiplication ( C = AxB )
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_bxb < _Int,
      _Flt >::HMatrixByHMatrixMultiply_BxB (int _blksize, int _nblks, int *_blks,
                                            const CBMatrix < _Int, _Flt > &_hmatr_a_bxb,
                                            const CBMatrix < _Int, _Flt > &_hmatr_b_bxb,
                                            CBMatrix < _Int, _Flt > &_hmatr_c_bxb)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int b_2 = _blksize * _blksize;

// Open both submatrices

      const CMatrix < int, float >*pHMatr_a = _hmatr_a_bxb.GetHMatrStr ();
      const CMatrix < _Int, _Flt > *pASub_a = _hmatr_a_bxb.GetASubArr ();

      const int *pia_hmatr_a = pHMatr_a->GetIaArr ();
      const int *pja_hmatr_a = pHMatr_a->GetJaArr ();

      const CMatrix < int, float >*pHMatr_b = _hmatr_b_bxb.GetHMatrStr ();
      const CMatrix < _Int, _Flt > *pASub_b = _hmatr_b_bxb.GetASubArr ();

      int nzja_b = pHMatr_b->GetNzja ();
      const int *pia_hmatr_b = pHMatr_b->GetIaArr ();
      const int *pja_hmatr_b = pHMatr_b->GetJaArr ();

// Compute block sparsity of a product

      CMatrix < int, float >*pHMatr_c = _hmatr_c_bxb.GetHMatrStr ();

      {

         CVectorData < int >imask (4 * _nblks);
         int *pimask = imask.Ptr ();

         int i;

         for (i = 0; i < 2 * _nblks; i++)
            pimask[i] = -1;

         int icycle = -1;

         CMatrix < int, float >::MatrixByMatrixMultiplySp (icycle, pimask,
                                                           pimask + _nblks,
                                                           pimask + 2 * _nblks,
                                                           pimask + 3 * _nblks, *pHMatr_a,
                                                           *pHMatr_b, *pHMatr_c);

      }

      int *pia_hmatr_c = pHMatr_c->GetIaArr ();
      int *pja_hmatr_c = pHMatr_c->GetJaArr ();

// Compute transposed sparsity of b and index array to initial sparsity

      CMatrix < int, float >hmatr_b_T;

      {

         CVectorData < int >imask (5 * _nblks + 1);
         int *pimask = imask.Ptr ();

         int i;

         for (i = 0; i < _nblks; i++)
            pimask[i] = -1;

         int icycle = -1;

         pHMatr_b->TransposedSparsityListSp (icycle, pimask, pimask + _nblks,
                                             pimask + 2 * _nblks, pimask + 3 * _nblks,
                                             pimask + 4 * _nblks, hmatr_b_T);

      }

      int *pia_hmatr_b_T = hmatr_b_T.GetIaArr ();
      int *pja_hmatr_b_T = hmatr_b_T.GetJaArr ();

      CVectorData < int >indt2ini_b (nzja_b);

      int *pindt2ini_b = indt2ini_b.Ptr ();

      {

         CVectorData < int >iptr (_nblks);
         int *piptr = iptr.Ptr ();

         int i, j, jj, k;

         for (i = 0; i < _nblks; i++)
            piptr[i] = pia_hmatr_b[i];

         for (i = 0; i < _nblks; i++) {
            for (j = pia_hmatr_b_T[i]; j < pia_hmatr_b_T[i + 1]; j++) {
               jj = pja_hmatr_b_T[j];
               k = piptr[jj];
               pindt2ini_b[j] = k;
               piptr[jj]++;
            }
         }

      }

// Allocate and register blocks of C

      int nzja_C = pia_hmatr_c[_nblks];

      _hmatr_c_bxb.ResizeASub (nzja_C);
      _hmatr_c_bxb.SetNzblk (nzja_C);

      CMatrix < _Int, _Flt > *pASub_c = _hmatr_c_bxb.GetASubArr ();

// For each block of C prepare list of A and B multiplication pairs

      CVectorData < int >ia_mvm_pairs (nzja_C + 1);
      int *pia_mvm_pairs = ia_mvm_pairs.Ptr ();

      vector < const CMatrix < _Int, _Flt > *>mvm_pairs;

      {

         int nzja_pairs = 0;
         pia_mvm_pairs[0] = 0;

         int i, j, jj, ip_a, ip_bt, ipend_a, ipend_bt, jj_a, jj_bt, ind_b;

         for (i = 0; i < _nblks; i++) {
            for (j = pia_hmatr_c[i]; j < pia_hmatr_c[i + 1]; j++) {
               jj = pja_hmatr_c[j];
               ip_a = pia_hmatr_a[i];
               ip_bt = pia_hmatr_b_T[jj];
               ipend_a = pia_hmatr_a[i + 1] - 1;
               ipend_bt = pia_hmatr_b_T[jj + 1] - 1;
               while (ip_a <= ipend_a && ip_bt <= ipend_bt) {
                  jj_a = pja_hmatr_a[ip_a];
                  jj_bt = pja_hmatr_b_T[ip_bt];
                  if (jj_a == jj_bt) {
                     mvm_pairs.push_back (pASub_a + ip_a);
                     ind_b = pindt2ini_b[ip_bt];
                     mvm_pairs.push_back (pASub_b + ind_b);
                     nzja_pairs++;
                     ip_a++;
                     ip_bt++;
                  } else if (jj_a < jj_bt) {
                     ip_a++;
                  } else if (jj_a > jj_bt) {
                     ip_bt++;
                  }
               }
               pia_mvm_pairs[j + 1] = nzja_pairs;
            }
         }
      }

// Perform parallel multiplications by block rows

      {

         int nsup_max = 0;

         {
            int i, ni;
            for (i = 0; i < _nblks; i++) {
               ni = _blks[i + 1] - _blks[i];
               if (ni > nsup_max)
                  nsup_max = ni;
            }
         }

         CVectorData < int >icycle_thr (n_thr);
         CVectorData < CVectorData < int > >imask_thr (n_thr);
         CVectorData < CVectorData < _Flt > >fmask_thr (n_thr);

         int *picycle_thr = icycle_thr.Ptr ();
         CVectorData < int >*pimask_thr = imask_thr.Ptr ();
         CVectorData < _Flt > *pfmask_thr = fmask_thr.Ptr ();

         {
            int i;
            for (i = 0; i < n_thr; i++)
               picycle_thr[i] = -1;
         }

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nzja_C; ipar++) {

            int my_thr = 0;
#ifdef USE_THREADS
            my_thr = omp_get_thread_num ();
#endif

            int icycle_th = picycle_thr[my_thr];

            if (icycle_th == -1) {
               pimask_thr[my_thr].resize (nsup_max * 4);
               pfmask_thr[my_thr].resize (nsup_max * b_2);
               int *ppimask_th = pimask_thr[my_thr].Ptr ();
               int i;
               for (i = 0; i < nsup_max * 2; i++)
                  ppimask_th[i] = -1;
            }

            int *ppimask_th = pimask_thr[my_thr].Ptr ();
            _Flt *ppfmask_th = pfmask_thr[my_thr].Ptr ();

            int n_mult = pia_mvm_pairs[ipar + 1] - pia_mvm_pairs[ipar];

            CVectorData < CMatrix < _Int, _Flt > >mmm_arr (n_mult * 2);

            CMatrix < _Int, _Flt > *pmmm_arr = mmm_arr.Ptr ();

            int j;

// Multiply

            for (j = pia_mvm_pairs[ipar]; j < pia_mvm_pairs[ipar + 1]; j++) {

               int jloc = j - pia_mvm_pairs[ipar];

               const CMatrix < _Int, _Flt > *pa_temp = mvm_pairs[j * 2];
               const CMatrix < _Int, _Flt > *pb_temp = mvm_pairs[j * 2 + 1];

               CFct_bxb < _Int, _Flt >::MatrixByMatrixMultiply_BxB (_blksize, icycle_th,
                                                                    ppimask_th,
                                                                    ppimask_th + nsup_max,
                                                                    ppimask_th +
                                                                    2 * nsup_max,
                                                                    ppimask_th +
                                                                    3 * nsup_max,
                                                                    ppfmask_th, *pa_temp,
                                                                    *pb_temp,
                                                                    pmmm_arr[jloc]);

            }

// Add by binary adding algorithm

            CMatrix < _Int, _Flt > *pmmm_arr_1 = pmmm_arr;
            CMatrix < _Int, _Flt > *pmmm_arr_2 = pmmm_arr + n_mult;

            CMatrix < _Int, _Flt > *pmmm_arr_curr = pmmm_arr_1;
            CMatrix < _Int, _Flt > *pmmm_arr_next = pmmm_arr_2;

            int n_mult_curr = n_mult;

            while (n_mult_curr > 1) {
               int n_mult_next = 0;
               int n_add = n_mult_curr / 2;
               int k;
               for (k = 0; k < n_add; k++) {
                  pmmm_arr_next[k].AddBlocksBxB ('+', _blksize, pmmm_arr_curr[k * 2],
                                                 pmmm_arr_curr[k * 2 + 1]);
               }
               n_mult_next = n_add;
               if (n_mult_next * 2 != n_mult_curr) {
                  pmmm_arr_next[n_mult_next].ReplaceFree (pmmm_arr_curr[n_mult_curr - 1]);
                  n_mult_next++;
               }
               n_mult_curr = n_mult_next;
               if (pmmm_arr_curr == pmmm_arr_1) {
                  pmmm_arr_curr = pmmm_arr_1;
                  pmmm_arr_next = pmmm_arr_2;
               } else {
                  pmmm_arr_curr = pmmm_arr_2;
                  pmmm_arr_next = pmmm_arr_1;
               }
            }

            pASub_c[j].ReplaceFree (pmmm_arr_curr[0]);

            picycle_thr[my_thr] = icycle_th;

         }
      }


   }

//
// Multiply by super sparse small blocks based matrix by rows and add result into prescribed positions
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CMvmSlv_BxB_impl < _Int, _Flt,
      _FltVect >::MvmA (char _oper, int _b_size, int _nlist,
                        const vector < _Int > &_list_alu, const vector < _Int > &_ia_alu,
                        const vector < _Int > &_ja_alu, const vector < _Flt > &_a_alu,
                        const _FltVect * _x, _FltVect * _ax)
   {

      int b_2 = _b_size * _b_size;

      const _Int *p_list_a = &_list_alu[0];
      const _Int *p_ia_a = &_ia_alu[0];
      const _Int *p_ja_a = &_ja_alu[0];
      const _Flt *p_a_a = &_a_alu[0];

      int i, irow, j, jj, kii, kjj, ibs_i, ibs_j, ibs_a;

      if (_oper == '+') {

         for (i = 0; i < _nlist; i++) {
            irow = (int) p_list_a[i];
            ibs_i = irow * _b_size;
            for (j = (int) p_ia_a[i]; j < p_ia_a[i + 1]; j++) {
               jj = (int) p_ja_a[j];
               ibs_j = jj * _b_size;
               ibs_a = j * b_2;
               for (kii = 0; kii < _b_size; kii++) {
                  for (kjj = 0; kjj < _b_size; kjj++) {
                     _ax[ibs_i + kii] +=
                        p_a_a[ibs_a + kjj * _b_size + kii] * _x[ibs_j + kjj];
                  }
               }
            }
         }

      } else if (_oper == '-') {

         for (i = 0; i < _nlist; i++) {
            irow = (int) p_list_a[i];
            ibs_i = irow * _b_size;
            for (j = (int) p_ia_a[i]; j < p_ia_a[i + 1]; j++) {
               jj = (int) p_ja_a[j];
               ibs_j = jj * _b_size;
               ibs_a = j * b_2;
               for (kii = 0; kii < _b_size; kii++) {
                  for (kjj = 0; kjj < _b_size; kjj++) {
                     _ax[ibs_i + kii] -=
                        p_a_a[ibs_a + kjj * _b_size + kii] * _x[ibs_j + kjj];
                  }
               }
            }
         }

      }

   }

//
// Multiply by super sparse small blocks based matrix by columns add result into prescribed positions
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CMvmSlv_BxB_impl < _Int, _Flt,
      _FltVect >::MvmAT (char _oper, int _b_size, int _nlist,
                         const vector < _Int > &_list_alu, const vector < _Int > &_ia_alu,
                         const vector < _Int > &_ja_alu, const vector < _Flt > &_a_alu,
                         const _FltVect * _x, _FltVect * _ax)
   {

      int b_2 = _b_size * _b_size;

      const _Int *p_list_a = &_list_alu[0];
      const _Int *p_ia_a = &_ia_alu[0];
      const _Int *p_ja_a = &_ja_alu[0];
      const _Flt *p_a_a = &_a_alu[0];

      int i, irow, j, jj, kii, kjj, ibs_i, ibs_j, ibs_a;

      if (_oper == '+') {

         for (i = 0; i < _nlist; i++) {
            irow = (int) p_list_a[i];
            ibs_i = irow * _b_size;
            for (j = (int) p_ia_a[i]; j < p_ia_a[i + 1]; j++) {
               jj = (int) p_ja_a[j];
               ibs_j = jj * _b_size;
               ibs_a = j * b_2;
               for (kii = 0; kii < _b_size; kii++) {
                  for (kjj = 0; kjj < _b_size; kjj++) {
                     _ax[ibs_j + kjj] +=
                        p_a_a[ibs_a + kjj * _b_size + kii] * _x[ibs_i + kii];
                  }
               }
            }
         }

      } else if (_oper == '-') {

         for (i = 0; i < _nlist; i++) {
            irow = (int) p_list_a[i];
            ibs_i = irow * _b_size;
            for (j = (int) p_ia_a[i]; j < p_ia_a[i + 1]; j++) {
               jj = (int) p_ja_a[j];
               ibs_j = jj * _b_size;
               ibs_a = j * b_2;
               for (kii = 0; kii < _b_size; kii++) {
                  for (kjj = 0; kjj < _b_size; kjj++) {
                     _ax[ibs_j + kjj] -=
                        p_a_a[ibs_a + kjj * _b_size + kii] * _x[ibs_i + kii];
                  }
               }
            }
         }

      }

   }

//
// Multiply by Q for super sparse small blocks based matrix Q stored as set of block Housholder transformations with external zero block
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CMvmSlv_BxB_impl < _Int, _Flt,
      _FltVect >::MvmQ (int _b_size, int _nlist, int _nrows, const vector < _Int > &_ia_q,
                        const vector < _Int > &_ja_q, const vector < _Flt > &_a_q,
                        const vector < _Flt > &_q_diag, const vector < _Flt > &_tau_arr,
                        const _FltVect * _x_diag, const _FltVect * _x,
                        _FltVect * _qx_diag, _FltVect * _qx, _FltVect * _work)
   {

      int b_2 = _b_size * _b_size;

      const _Int *p_ia_q = &_ia_q[0];
      const _Int *p_ja_q = &_ja_q[0];
      const _Flt *p_a_q = &_a_q[0];

      const _Flt *p_q_diag = &_q_diag[0];
      const _Flt *p_tau_arr = &_tau_arr[0];

      CVector < _FltVect >::CopyVector (_nlist * _b_size, _x_diag, _qx_diag);
      CVector < _FltVect >::CopyVector (_nrows * _b_size, _x, _qx);

      _FltVect *pwork = _work;
      _FltVect *pwork1 = _work + _b_size;

      int i, j, jj, kii, kjj, ibs_j, ibs_q;

      for (i = _nlist - 1; i >= 0; i--) {

// BDot

         CVector < _FltVect >::SetByZeroes (_b_size, pwork);

         for (j = (int) p_ia_q[i]; j < p_ia_q[i + 1]; j++) {
            jj = (int) p_ja_q[j];
            ibs_j = jj * _b_size;
            ibs_q = j * b_2;
            for (kii = 0; kii < _b_size; kii++) {
               for (kjj = 0; kjj < _b_size; kjj++) {
                  pwork[kii] += p_a_q[ibs_q + kjj * _b_size + kii] * _qx[ibs_j + kjj];
               }
            }
         }
         ibs_j = i * _b_size;
         ibs_q = i * b_2;
         for (kii = 0; kii < _b_size; kii++) {
            for (kjj = 0; kjj < _b_size; kjj++) {
               pwork[kii] +=
                  p_q_diag[ibs_q + kjj * _b_size + kii] * _qx_diag[ibs_j + kjj];
            }
         }

// Apply tau

         CVector < _FltVect >::SetByZeroes (_b_size, pwork1);

         for (kii = 0; kii < _b_size; kii++) {
            for (kjj = 0; kjj < _b_size; kjj++) {
               pwork1[kii] += p_tau_arr[ibs_q + kii * _b_size + kjj] * pwork[kjj];
            }
         }

// BDaxpy

         for (j = (int) p_ia_q[i]; j < p_ia_q[i + 1]; j++) {
            jj = (int) p_ja_q[j];
            ibs_j = jj * _b_size;
            ibs_q = j * b_2;
            for (kii = 0; kii < _b_size; kii++) {
               for (kjj = 0; kjj < _b_size; kjj++) {
                  _qx[ibs_j + kii] += p_a_q[ibs_q + kii * _b_size + kjj] * pwork1[kjj];
               }
            }
         }
         ibs_j = i * _b_size;
         ibs_q = i * b_2;
         for (kii = 0; kii < _b_size; kii++) {
            for (kjj = 0; kjj < _b_size; kjj++) {
               _qx_diag[ibs_j + kii] +=
                  p_q_diag[ibs_q + kii * _b_size + kjj] * pwork1[kjj];
            }
         }

      }

   }

//
// Multiply by QT for super sparse small blocks based matrix Q stored as set of block Housholder transformations with external zero block
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CMvmSlv_BxB_impl < _Int, _Flt,
      _FltVect >::MvmQT (int _b_size, int _nlist, int _nrows,
                         const vector < _Int > &_ia_q, const vector < _Int > &_ja_q,
                         const vector < _Flt > &_a_q, const vector < _Flt > &_q_diag,
                         const vector < _Flt > &_tau_arr, const _FltVect * _x_diag,
                         const _FltVect * _x, _FltVect * _qx_diag, _FltVect * _qx,
                         _FltVect * _work)
   {

      int b_2 = _b_size * _b_size;

      const _Int *p_ia_q = &_ia_q[0];
      const _Int *p_ja_q = &_ja_q[0];
      const _Flt *p_a_q = &_a_q[0];

      const _Flt *p_q_diag = &_q_diag[0];
      const _Flt *p_tau_arr = &_tau_arr[0];

      CVector < _FltVect >::CopyVector (_nlist * _b_size, _x_diag, _qx_diag);
      CVector < _FltVect >::CopyVector (_nrows * _b_size, _x, _qx);

      _FltVect *pwork = _work;
      _FltVect *pwork1 = _work + _b_size;

      int i, j, jj, kii, kjj, ibs_j, ibs_q;

      for (i = 0; i < _nlist; i++) {

// BDot

         CVector < _FltVect >::SetByZeroes (_b_size, pwork);

         for (j = (int) p_ia_q[i]; j < p_ia_q[i + 1]; j++) {
            jj = (int) p_ja_q[j];
            ibs_j = jj * _b_size;
            ibs_q = j * b_2;
            for (kii = 0; kii < _b_size; kii++) {
               for (kjj = 0; kjj < _b_size; kjj++) {
                  pwork[kii] += p_a_q[ibs_q + kjj * _b_size + kii] * _qx[ibs_j + kjj];
               }
            }
         }
         ibs_j = i * _b_size;
         ibs_q = i * b_2;
         for (kii = 0; kii < _b_size; kii++) {
            for (kjj = 0; kjj < _b_size; kjj++) {
               pwork[kii] +=
                  p_q_diag[ibs_q + kjj * _b_size + kii] * _qx_diag[ibs_j + kjj];
            }
         }

// Apply tau

         CVector < _FltVect >::SetByZeroes (_b_size, pwork1);

         for (kii = 0; kii < _b_size; kii++) {
            for (kjj = 0; kjj < _b_size; kjj++) {
               pwork1[kii] += p_tau_arr[ibs_q + kjj * _b_size + kii] * pwork[kjj];
            }
         }

// BDaxpy

         for (j = (int) p_ia_q[i]; j < p_ia_q[i + 1]; j++) {
            jj = (int) p_ja_q[j];
            ibs_j = jj * _b_size;
            ibs_q = j * b_2;
            for (kii = 0; kii < _b_size; kii++) {
               for (kjj = 0; kjj < _b_size; kjj++) {
                  _qx[ibs_j + kii] += p_a_q[ibs_q + kii * _b_size + kjj] * pwork1[kjj];
               }
            }
         }
         ibs_j = i * _b_size;
         ibs_q = i * b_2;
         for (kii = 0; kii < _b_size; kii++) {
            for (kjj = 0; kjj < _b_size; kjj++) {
               _qx_diag[ibs_j + kii] +=
                  p_q_diag[ibs_q + kii * _b_size + kjj] * pwork1[kjj];
            }
         }

      }

   }

//
// Solve with L, L is stored by columns (diag is inverted)
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CMvmSlv_BxB_impl < _Int, _Flt,
      _FltVect >::SolveL (int _b_size, int _n, const vector < _Int > &_ia_l,
                          const vector < _Int > &_ja_l, const vector < _Flt > &_a_l,
                          const _FltVect * _x, _FltVect * _lx, _FltVect * _work)
   {

      static const _Flt fzero = (_FltVect) 0.0e0;

      int b_2 = _b_size * _b_size;

      int i;

      if (_lx != _x) {
         for (i = 0; i < _n * _b_size; i++)
            _lx[i] = _x[i];
      }

      const _Int *p_ia_l = &_ia_l[0];
      const _Int *p_ja_l = &_ja_l[0];
      const _Flt *p_a_l = &_a_l[0];

      int j, jj, ibeg, iend, kii, kjj, ibs_i, ibs_j, ibs_l;

      for (i = 0; i < _n; i++) {
         ibs_i = i * _b_size;
         ibeg = (int) p_ia_l[i];
         iend = (int) p_ia_l[i + 1] - 1;
         ibs_l = ibeg * b_2;
         for (kii = 0; kii < _b_size; kii++)
            _work[kii] = fzero;
         for (kii = 0; kii < _b_size; kii++) {
            for (kjj = 0; kjj < _b_size; kjj++) {
               _work[kjj] += p_a_l[ibs_l + kjj * _b_size + kii] * _lx[ibs_i + kii];
            }
         }
         for (kii = 0; kii < _b_size; kii++)
            _lx[ibs_i + kii] = _work[kii];
         for (j = ibeg + 1; j <= iend; j++) {
            jj = (int) p_ja_l[j];
            ibs_j = jj * _b_size;
            ibs_l = j * b_2;
            for (kii = 0; kii < _b_size; kii++) {
               for (kjj = 0; kjj < _b_size; kjj++) {
                  _lx[ibs_j + kjj] -=
                     p_a_l[ibs_l + kjj * _b_size + kii] * _lx[ibs_i + kii];
               }
            }
         }
      }

   }

//
// Solve with L, L is stored by small blocks rows (small blocks diag is inverted)
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CMvmSlv_BxB_impl < _Int, _Flt,
      _FltVect >::SolveLRow (int _b_size, int _n, const vector < _Int > &_ia_l,
                             const vector < _Int > &_ja_l, const vector < _Flt > &_a_l,
                             const _FltVect * _x, _FltVect * _lx, _FltVect * _work)
   {

      static const _Flt fzero = (_FltVect) 0.0e0;

      int b_2 = _b_size * _b_size;

      int i;

      const _Int *p_ia_l = &_ia_l[0];
      const _Int *p_ja_l = &_ja_l[0];
      const _Flt *p_a_l = &_a_l[0];

      if (_lx != _x) {
         for (i = 0; i < _n * _b_size; i++)
            _lx[i] = _x[i];
      }

      int j, jj, ibeg, iend, kii, kjj, ibs_i, ibs_j, ibs_l;

      for (i = 0; i < _n; i++) {
         ibeg = (int) p_ia_l[i];
         iend = (int) p_ia_l[i + 1] - 1;
         ibs_i = i * _b_size;
         for (j = ibeg; j < iend; j++) {
            jj = (int) p_ja_l[j];
            ibs_j = jj * _b_size;
            ibs_l = j * b_2;
            for (kii = 0; kii < _b_size; kii++) {
               for (kjj = 0; kjj < _b_size; kjj++) {
                  _lx[ibs_i + kii] -=
                     p_a_l[ibs_l + kjj * _b_size + kii] * _lx[ibs_j + kjj];
               }
            }
         }
         ibs_l = iend * b_2;
         for (kii = 0; kii < _b_size; kii++)
            _work[kii] = fzero;
         for (kii = 0; kii < _b_size; kii++) {
            for (kjj = 0; kjj < _b_size; kjj++) {
               _work[kii] += p_a_l[ibs_l + kjj * _b_size + kii] * _lx[ibs_i + kjj];
            }
         }
         for (kii = 0; kii < _b_size; kii++)
            _lx[ibs_i + kii] = _work[kii];
      }

   }

//
// Solve with U, U is stored by rows (diag is inverted)
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CMvmSlv_BxB_impl < _Int, _Flt,
      _FltVect >::SolveU (int _b_size, int _n, const vector < _Int > &_ia_u,
                          const vector < _Int > &_ja_u, const vector < _Flt > &_a_u,
                          const _FltVect * _x, _FltVect * _ux, _FltVect * _work)
   {

      static const _Flt fzero = (_FltVect) 0.0e0;

      int b_2 = _b_size * _b_size;

      int i;

      const _Int *p_ia_u = &_ia_u[0];
      const _Int *p_ja_u = &_ja_u[0];
      const _Flt *p_a_u = &_a_u[0];

      if (_ux != _x) {
         for (i = 0; i < _n * _b_size; i++)
            _ux[i] = _x[i];
      }

      int j, jj, ibeg, iend, kii, kjj, ibs_i, ibs_j, ibs_u;

      for (i = _n - 1; i >= 0; i--) {
         ibeg = (int) p_ia_u[i];
         iend = (int) p_ia_u[i + 1] - 1;
         ibs_i = i * _b_size;
         for (j = ibeg + 1; j <= iend; j++) {
            jj = (int) p_ja_u[j];
            ibs_j = jj * _b_size;
            ibs_u = j * b_2;
            for (kii = 0; kii < _b_size; kii++) {
               for (kjj = 0; kjj < _b_size; kjj++) {
                  _ux[ibs_i + kii] -=
                     p_a_u[ibs_u + kjj * _b_size + kii] * _ux[ibs_j + kjj];
               }
            }
         }
         ibs_u = ibeg * b_2;
         for (kii = 0; kii < _b_size; kii++)
            _work[kii] = fzero;
         for (kii = 0; kii < _b_size; kii++) {
            for (kjj = 0; kjj < _b_size; kjj++) {
               _work[kii] += p_a_u[ibs_u + kjj * _b_size + kii] * _ux[ibs_i + kjj];
            }
         }
         for (kii = 0; kii < _b_size; kii++)
            _ux[ibs_i + kii] = _work[kii];
      }

   }

//
// Solve with U, U is stored by small blocks columns (small blocks diag is inverted)
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CMvmSlv_BxB_impl < _Int, _Flt,
      _FltVect >::SolveUCol (int _b_size, int _n, const vector < _Int > &_ia_u,
                             const vector < _Int > &_ja_u, const vector < _Flt > &_a_u,
                             const _FltVect * _x, _FltVect * _ux, _FltVect * _work)
   {

      static const _Flt fzero = (_FltVect) 0.0e0;

      int b_2 = _b_size * _b_size;

      int i;

      if (_ux != _x) {
         for (i = 0; i < _n * _b_size; i++)
            _ux[i] = _x[i];
      }

      const _Int *p_ia_u = &_ia_u[0];
      const _Int *p_ja_u = &_ja_u[0];
      const _Flt *p_a_u = &_a_u[0];

      int j, jj, ibeg, iend, kii, kjj, ibs_i, ibs_j, ibs_u;

      for (i = _n - 1; i >= 0; i--) {
         ibs_i = i * _b_size;
         ibeg = (int) p_ia_u[i];
         iend = (int) p_ia_u[i + 1] - 1;
         ibs_u = iend * b_2;
         for (kii = 0; kii < _b_size; kii++)
            _work[kii] = fzero;
         for (kii = 0; kii < _b_size; kii++) {
            for (kjj = 0; kjj < _b_size; kjj++) {
               _work[kjj] += p_a_u[ibs_u + kjj * _b_size + kii] * _ux[ibs_i + kii];
            }
         }
         for (kii = 0; kii < _b_size; kii++)
            _ux[ibs_i + kii] = _work[kii];
         for (j = ibeg; j <= iend - 1; j++) {
            jj = (int) p_ja_u[j];
            ibs_j = jj * _b_size;
            ibs_u = j * b_2;
            for (kii = 0; kii < _b_size; kii++) {
               for (kjj = 0; kjj < _b_size; kjj++) {
                  _ux[ibs_j + kjj] -=
                     p_a_u[ibs_u + kjj * _b_size + kii] * _ux[ibs_i + kii];
               }
            }
         }
      }

   }

// Prepare matrix data only
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::StoreMatrix_BxB_Point (void *_pcomm, int _blksize, int _nhblks,
                                         int *_hblk2cpu, int *_hblk2blks,
                                         long long *_hblks, int _nblks, long long *_blks,
                                         _Int * _ia, _Int * _ja, _Flt * _a,
                                         CBMatrix < _Int, _Flt > *_hmatr_arr)
   {

      int myid = CMPIDataExchange::GetMyid (_pcomm);

// Compute extended blocks

      int b_2 = _blksize * _blksize;

      CVectorData < long long >hblks_pt (_nhblks + 1);
      CVectorData < long long >blks_pt (_nblks + 1);

      long long *phblks_pt = hblks_pt.Ptr ();
      long long *pblks_pt = blks_pt.Ptr ();

      int i;

      for (i = 0; i <= _nhblks; i++)
         phblks_pt[i] = _blksize * _hblks[i];
      for (i = 0; i <= _nblks; i++)
         pblks_pt[i] = _blksize * _blks[i];

      CVectorData < int >blk2hblks (_nblks);
      int *pblk2hblks = blk2hblks.Ptr ();

      {
         int j;
         for (i = 0; i < _nhblks; i++) {
            for (j = _hblk2blks[i]; j < _hblk2blks[i + 1]; j++) {
               pblk2hblks[j] = i;
            }
         }
      }

// Create matrix as set of hblocks

      int nimax = 0;

      int niloc;

      for (i = 0; i < _nblks; i++) {
         niloc = (int) (pblks_pt[i + 1] - pblks_pt[i]);
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
      CVectorData < int >imaskpt_thr (n_thr * nimax * 4);

      int *picycle_thr = &icycle_thr[0];
      CVectorData < _Int > *plistloc_thr = &listloc_thr[0];
      CVectorData < int >*pimaskblk_thr = &imaskblk_thr[0];
      int *pimaskpt_thr = imaskpt_thr.Ptr ();

      for (i = 0; i < n_thr; i++)
         picycle_thr[i] = -1;

      int ibeg = 0;

      for (i = 0; i < _nhblks; i++) {
         if (_hblk2cpu[i] == myid) {

            int niblk = _hblk2blks[i + 1] - _hblk2blks[i];
            vector < CBMatrix < _Int, _Flt > >hblk_arr (niblk + 1);
            CBMatrix < _Int, _Flt > *phblk_arr = &hblk_arr[0];
            int iblk0 = _hblk2blks[i];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = _hblk2blks[i]; ipar < _hblk2blks[i + 1]; ipar++) {

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
                  int *pimask = pimaskpt_thr + 4 * nimax * my_thr;
                  for (j = 0; j < nimax; j++)
                     pimask[j] = -1;
               }
               int icycleblk = picycle_thr[my_thr];
               _Int *plistloc = plistloc_thr[my_thr].Ptr ();
               _Int *pialoc = plistloc + nimax;
               int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
               int *pimask = pimaskpt_thr + 4 * nimax * my_thr;
               icycleblk++;
               int niloc = (int) (pblks_pt[ipar + 1] - pblks_pt[ipar]);
               int ibs = (int) (pblks_pt[ipar] - pblks_pt[iblk0]);
               int ishift = (int) _ia[ibeg + ibs];
               for (j = 0; j < niloc; j++) {
                  plistloc[j] = (_Int) (j + pblks_pt[ipar]);
               }
               for (j = 0; j <= niloc; j++) {
                  pialoc[j] = _ia[ibeg + ibs + j] - _ia[ibeg + ibs];
               }
               CBMatrix < _Int, _Flt > hblk_loc (ipar, niloc, plistloc, pialoc,
                                                 _ja + ishift, _a + ishift, _nblks,
                                                 pblks_pt, icycleblk, pimaskblk);
               int nzblk = hblk_loc.GetNzblk ();
               CMatrix < _Int, _Flt > *pASub = hblk_loc.GetASubArr ();
               for (j = 0; j < nzblk; j++) {

                  int nlist_temp = pASub[j].GetNlist ();
                  vector < _Int > *plist_temp = pASub[j].GetList ();
                  vector < _Int > *pia_temp = pASub[j].GetIa ();
                  vector < _Int > *pja_temp = pASub[j].GetJa ();
                  vector < _Flt > *pa_temp = pASub[j].GetA ();

                  int nlist_cnd;
                  vector < _Int > list_cnd;
                  vector < _Int > ia_cnd;
                  vector < _Int > ja_cnd;
                  vector < _Flt > a_cnd;

                  CFct_bxb_impl < _Int, _Flt >::CondenseRectMatrix (nlist_temp,
                                                                    *plist_temp,
                                                                    *pia_temp, *pja_temp,
                                                                    *pa_temp, _blksize,
                                                                    nlist_cnd, list_cnd,
                                                                    ia_cnd, ja_cnd, a_cnd,
                                                                    nimax, icycleblk,
                                                                    pimask);

                  int nzja_cnd = (int) ia_cnd[nlist_cnd];

                  pASub[j].SetBSize (_blksize);
                  pASub[j].SetNlist (nlist_cnd);
                  pASub[j].SetNzja (nzja_cnd);
                  pASub[j].SetNza (nzja_cnd * b_2);

                  plist_temp->swap (list_cnd);
                  pia_temp->swap (ia_cnd);
                  pja_temp->swap (ja_cnd);
                  pa_temp->swap (a_cnd);

               }
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
                  kk2 = pblk2hblks[kk];
                  pja_str[nzja_tot] = kk - _hblk2blks[kk2];
                  pja2_str[nzja_tot] = kk2;
                  pASub[nzja_tot].ReplaceFree (pASub_temp[k]);
                  nzja_tot++;
               }
               pia_str[j + 1] = nzja_tot;
            }
            _hmatr_arr[i].ReplaceFree (hblk);
            ibeg += (int) (phblks_pt[i + 1] - phblks_pt[i]);
         }
      }

   }

// Prepare matrix data only
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::StoreMatrix_BxB_Block (void *_pcomm, int _blksize, int _nhblks,
                                         int *_hblk2cpu, int *_hblk2blks,
                                         long long *_hblks, int _nblks, long long *_blks,
                                         _Int * _ia, _Int * _ja, _Flt * _a,
                                         CBMatrix < _Int, _Flt > *_hmatr_arr)
   {

      int myid = CMPIDataExchange::GetMyid (_pcomm);

// Compute extended blocks

      int b_2 = _blksize * _blksize;

      CVectorData < int >blk2hblks (_nblks);
      int *pblk2hblks = blk2hblks.Ptr ();

      int i;

      {
         int j;
         for (i = 0; i < _nhblks; i++) {
            for (j = _hblk2blks[i]; j < _hblk2blks[i + 1]; j++) {
               pblk2hblks[j] = i;
            }
         }
      }

// Create matrix as set of hblocks

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
      vector < CVectorData < _Int > >listloc_thr (n_thr + 1);
      vector < CVectorData < int > >imaskblk_thr (n_thr + 1);
      CVectorData < int >imaskpt_thr (n_thr * nimax * 4);

      int *picycle_thr = &icycle_thr[0];
      CVectorData < _Int > *plistloc_thr = &listloc_thr[0];
      CVectorData < int >*pimaskblk_thr = &imaskblk_thr[0];
      int *pimaskpt_thr = imaskpt_thr.Ptr ();

      for (i = 0; i < n_thr; i++)
         picycle_thr[i] = -1;

      int ibeg = 0;

      for (i = 0; i < _nhblks; i++) {
         if (_hblk2cpu[i] == myid) {

            int niblk = _hblk2blks[i + 1] - _hblk2blks[i];
            vector < CBMatrix < _Int, _Flt > >hblk_arr (niblk + 1);
            CBMatrix < _Int, _Flt > *phblk_arr = &hblk_arr[0];
            int iblk0 = _hblk2blks[i];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = _hblk2blks[i]; ipar < _hblk2blks[i + 1]; ipar++) {

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
                  int *pimask = pimaskpt_thr + 4 * nimax * my_thr;
                  for (j = 0; j < nimax; j++)
                     pimask[j] = -1;
               }
               int icycleblk = picycle_thr[my_thr];
               _Int *plistloc = plistloc_thr[my_thr].Ptr ();
               _Int *pialoc = plistloc + nimax;
               int *pimaskblk = pimaskblk_thr[my_thr].Ptr ();
               int *pimask = pimaskpt_thr + 4 * nimax * my_thr;
               icycleblk++;
               int niloc = (int) (_blks[ipar + 1] - _blks[ipar]);
               int ibs = (int) (_blks[ipar] - _blks[iblk0]);
               int ishift = (int) _ia[ibeg + ibs];
               for (j = 0; j < niloc; j++) {
                  plistloc[j] = (_Int) (j + _blks[ipar]);
               }
               for (j = 0; j <= niloc; j++) {
                  pialoc[j] = _ia[ibeg + ibs + j] - _ia[ibeg + ibs];
               }

               CBMatrix < _Int, _Flt > hblk_loc (_blksize, ipar, niloc, plistloc, pialoc,
                                                 _ja + ishift, _a + ishift * b_2, _nblks,
                                                 _blks, icycleblk, pimaskblk);

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
                  kk2 = pblk2hblks[kk];
                  pja_str[nzja_tot] = kk - _hblk2blks[kk2];
                  pja2_str[nzja_tot] = kk2;
                  pASub[nzja_tot].ReplaceFree (pASub_temp[k]);
                  nzja_tot++;
               }
               pia_str[j + 1] = nzja_tot;
            }
            _hmatr_arr[i].ReplaceFree (hblk);
            ibeg += (int) (_hblks[i + 1] - _hblks[i]);
         }
      }

   }

// Prepare matrix data only
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::PrepareMatrixBase_BxB (void *_pcomm, int _blksize, int _nhblks,
                                         int *_hblk2cpu, int *_hblk2blks,
                                         long long *_hblks, int _nblks, long long *_blks,
                                         int _format_type, _Int * _ia, _Int * _ja,
                                         _Flt * _a)
   {

// Store control data

      this->pcomm = _pcomm;
      this->blksize = _blksize;
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

      CK3D_SolverBxBThreads < _Int, _Flt, _FltVect >::StoreMatrix_BxB (_pcomm, _blksize,
                                                                       _nhblks, _hblk2cpu,
                                                                       _hblk2blks, _hblks,
                                                                       _nblks, _blks,
                                                                       _format_type, _ia,
                                                                       _ja, _a,
                                                                       phmatr_arr);

// Symmetrize hmatr

      vector < CBMatrix < _Int, _Flt > >hmatr_symm_arr (_nhblks + 1);
      CBMatrix < _Int, _Flt > *phmatr_symm_arr = &hmatr_symm_arr[0];

      CBMatrix < _Int, _Flt >::SymmetrizeSubmatrices_BxB (_pcomm, _blksize, _nhblks,
                                                          phblk2cpu, phblk2blks,
                                                          pblk2hblks, pblks, phmatr_arr,
                                                          phmatr_symm_arr);

      for (i = 0; i < _nhblks; i++) {
         if (phblk2cpu[i] == myid) {
            phmatr_arr[i].ReplaceFree (phmatr_symm_arr[i]);
         }
      }

// Init MvmA structures

      this->mvm_bxb.InitControl_BxB (_pcomm, _blksize, _nhblks, phblk2cpu, phblk2blks,
                                     pblk2hblks, phblks, _nblks, pblks);
      this->mvm_bxb.InitMvmA_BxB (phmatr_arr);

      this->b_fast_transform = false;
      this->b_use_ini = false;
      this->b_use_wells = false;
      this->b_use_blksize = false;
      this->b_blk_wells = false;

   }

// Prepare matrix data only including preliminary repartitioning of the matrix
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::PrepareMatrixThreads_BxB (void *_pcomm, SParams * _params,
                                            SStatData * _stats, int _blksize, int _nhblks,
                                            int *_hblk2cpu, int *_hblk2blks_ini,
                                            long long *_hblks, int _nblks_ini,
                                            long long *_blks_ini, int _format_type,
                                            _Int * _ia, _Int * _ja, _Flt * _a)
   {

      _stats->dtime_part = 0.0e0;
      _stats->dtime_ord = 0.0e0;

//      ofstream ffout ("ChkDecomp.dat");

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

      CK3D_SolverBxBThreads < _Int, _Flt, _FltVect >::StoreMatrix_BxB (_pcomm, _blksize,
                                                                       _nhblks, _hblk2cpu,
                                                                       _hblk2blks_ini,
                                                                       _hblks, _nblks_ini,
                                                                       _blks_ini,
                                                                       _format_type, _ia,
                                                                       _ja, _a,
                                                                       phmatr_ini_arr);

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

         int i;

         for (i = 0; i < _nhblks; i++) {
            if (_hblk2cpu[i] == myid) {

               int niblk = _hblk2blks_ini[i + 1] - _hblk2blks_ini[i];
               int iblkbeg = _hblk2blks_ini[i];
//               int iblkend = _hblk2blks_ini[i + 1] - 1;

               CMatrix < _Int, _Flt > *pASub_curr = phmatr_ini_arr[i].GetASubArr ();
               CMatrix < int, float >*phmatr_str_curr = phmatr_ini_arr[i].GetHMatrStr ();

               int *pia_str = phmatr_str_curr->GetIaArr ();
               int *pja_str = phmatr_str_curr->GetJaArr ();
               int *pja2_str = phmatr_str_curr->GetJa2Arr ();

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
//                           pASub_dia[ibs] = pASub_curr[j];
                           pASub_dia[ibs].GetSparsity (pASub_curr[j]);
                           pASub_dia[ibs].SetBSize (1);
                           int nzja_temp = pASub_dia[ibs].GetNzja ();
                           pASub_dia[ibs].ResizeA (nzja_temp);
                           pASub_dia[ibs].SetNza (nzja_temp);
                           _Flt *pa_temp = pASub_dia[ibs].GetAArr ();
                           CVector < _Flt >::SetByOnes (nzja_temp, pa_temp);
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

//               ffout << " Ihblk = " << i << endl;
//               ffout << " Degree = " << degree_decomp << " isize_max = " << isize_max << " nparts = " << _params->nparts << endl;
//               hmatr_dia_symm.PrintHMatrix (ffout);

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

//               ffout << " nparts_new = " << nparts_new << endl;
//               PrintArray (ffout, " ppblks_temp", nparts_new+1, ppblks_temp);
//               PrintArray (ffout, " pporder_ini", nloc_dia, pporder_ini);

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
            this->PrepareFastTransform_BxB (phmatr_ini_arr);
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

         this->UseFastTransform_BxB (phmatr_ini_arr);

      } else {

         CBMatrix < _Int, _Flt >::ReorderHMatrixDiag_BxB (_pcomm, _blksize, _nhblks,
                                                          _hblk2cpu, _hblk2blks_ini,
                                                          _nblks_ini, _blks_ini,
                                                          phmatr_ini_arr, porder_ini,
                                                          phblk2blks, nblks_loc, pblks,
                                                          phmatr_arr);

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

         CBMatrix < _Int, _Flt >::SymmetrizeSubmatrices_BxB (_pcomm, _blksize, _nhblks,
                                                             _hblk2cpu, phblk2blks,
                                                             pblk2hblks, pblks,
                                                             phmatr_arr, phmatr_symm_arr);

         int i;

         for (i = 0; i < _nhblks; i++) {
            if (_hblk2cpu[i] == myid) {
               phmatr_arr[i].ReplaceFree (phmatr_symm_arr[i]);
            }
         }

      }
// Init MvmA structures

//      ofstream ffout ("ChkMvm.dat");
//      {
//         int i;
//         for (i=0;i<_nhblks;i++) {
//            ffout << " Ihblk = " << endl;
//            phmatr_arr[i].PrintHMatrix (ffout);
//         }
//      }

      this->mvm_bxb.InitControl_BxB (_pcomm, _blksize, _nhblks, _hblk2cpu, phblk2blks,
                                     pblk2hblks, _hblks, nblks_loc, pblks);
      this->mvm_bxb.InitMvmA_BxB (phmatr_arr);

      this->b_use_ini = true;

   }

// Prepare matrix data only including preliminary repartitioning of the matrix taking into account wells data
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::PrepareMatrixWells_BxB (void *_pcomm, SParams * _params,
                                          SStatData * _stats, int _blksize,
                                          int _nhblks_ini, int *_hblk2cpu_ini,
                                          int *_hblk2blks_ini, long long *_hblks_ini,
                                          int _nblks_ini, long long *_blks_ini,
                                          vector < int >*_p_ia_wells_ext,
                                          vector < int >*_p_ja_wells_ext,
                                          int _format_type, _Int * _ia, _Int * _ja,
                                          _Flt * _a)
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

// Create matrix as set of hblocks

      vector < CBMatrix < _Int, _Flt > >hmatr_ini_arr (_nhblks_ini + 1);
      CBMatrix < _Int, _Flt > *phmatr_ini_arr = &hmatr_ini_arr[0];

      CK3D_SolverBxBThreads < _Int, _Flt, _FltVect >::StoreMatrix_BxB (_pcomm, _blksize,
                                                                       _nhblks_ini,
                                                                       _hblk2cpu_ini,
                                                                       _hblk2blks_ini,
                                                                       _hblks_ini,
                                                                       _nblks_ini,
                                                                       _blks_ini,
                                                                       _format_type, _ia,
                                                                       _ja, _a,
                                                                       phmatr_ini_arr);

// Compute new local partitioning and ordering if necessary

      if (_params->b_new_part) {

         double time0;

         CMPIDataExchange::Synchronize (pcommloc);

         time0 = CMPIDataExchange::GetWallTimeMPI ();

// Symmetrize matrix data

         {

            vector < CBMatrix < _Int, _Flt > >hmatr_symm_arr (_nhblks_ini + 1);
            CBMatrix < _Int, _Flt > *phmatr_symm_arr = &hmatr_symm_arr[0];

            CBMatrix < _Int, _Flt >::SymmetrizeSubmatrices_BxB (_pcomm, _blksize,
                                                                _nhblks_ini,
                                                                _hblk2cpu_ini,
                                                                _hblk2blks_ini,
                                                                pblk2hblks_ini, _blks_ini,
                                                                phmatr_ini_arr,
                                                                phmatr_symm_arr);

            for (int i = 0; i < _nhblks_ini; i++) {
               if (_hblk2cpu_ini[i] == myid) {
                  phmatr_ini_arr[i].ReplaceFree (phmatr_symm_arr[i]);
               }
            }

         }

// Compute point comparison matrix data

         vector < CBMatrix < _Int, _Flt > >hmatr_pt_arr (_nhblks_ini + 1);
         CBMatrix < _Int, _Flt > *phmatr_pt_arr = &hmatr_pt_arr[0];

         {
            int i;
            for (i = 0; i < _nhblks_ini; i++) {
               if (_hblk2cpu_ini[i] == myid) {
                  CBMatrix < _Int, _Flt >::ComparisonHMatrix_BxB (_blksize,
                                                                  phmatr_ini_arr[i],
                                                                  phmatr_pt_arr[i]);
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
                                                  phmatr_pt_arr, n_wells_new, nhblks_new,
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
            this->PrepareFastTransform_BxB (phmatr_ini_arr);
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

         this->UseFastTransform_BxB (phmatr_ini_arr);

      } else {

         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (_pcomm, _blksize, _nhblks_ini,
                                                      _hblk2cpu_ini, _hblk2blks_ini,
                                                      _blks_ini, phmatr_ini_arr,
                                                      porder_wells, nhblks_loc, phblk2cpu,
                                                      phblk2blks, pblks, phmatr_arr);

      }

      phmatr_arr = &this->hmatr_arr[0];

      CMPIDataExchange::Synchronize (pcommloc);

      double time1;

      time1 = CMPIDataExchange::GetWallTimeMPI ();

      _stats->dtime_ord = time1 - time0;

// Init MvmA structures

      this->mvm_bxb.InitControl_BxB (_pcomm, _blksize, nhblks_loc, phblk2cpu, phblk2blks,
                                     pblk2hblks, phblks, nblks_loc, pblks);
      this->mvm_bxb.InitMvmA_BxB (phmatr_arr);

   }

// Perform multiplication by A
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CMvmParBxBThreads < _Int, _Flt,
      _FltVect >::MvmA_BxB (const _FltVect * _x, _FltVect * _ax)
   {

// Open mvm structure

      void *pcomm_loc = this->pcomm;

      int blksize_loc = this->blksize;
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

      CVector < _FltVect >::SetByZeroes_thr (ni_cpu_loc * blksize, _ax);

// Prepare send

      int ni_send_loc = pia_sends_loc[nsends_loc];

      int i, ind;

      for (i = 0; i < ni_send_loc; i++) {
         ind = pind_sends_loc[i];
         CVector < _FltVect >::CopyVector (blksize_loc, _x + ind * blksize_loc,
                                           px_send_loc + i * blksize_loc);
      }

      int ni_recv_loc = pia_recvs_loc[nrecvs_loc];

      CVector < _FltVect >::SetByZeroes_thr (ni_recv_loc * blksize_loc, px_recv_loc);

// Init async recvs and sends

      void *psndrcv_recvs_loc;
      void *psndrcv_stats_loc;

      CMPIDataExchange::AllocateRecvs (nrecvs_loc + nsends_loc, psndrcv_recvs_loc);
      CMPIDataExchange::AllocateStats (nrecvs_loc + nsends_loc, psndrcv_stats_loc);

      int icpu, isize, ibs;

      for (i = 0; i < nrecvs_loc; i++) {
         icpu = prcv2cpu_loc[i];
         isize =
            blksize_loc * (pia_recvs_loc[i + 1] - pia_recvs_loc[i]) * sizeof (_FltVect);
         ibs = pia_recvs_loc[i] * blksize_loc;
         CMPIDataExchange::IRecv (pcomm_loc, icpu, icpu, isize,
                                  (char *) (px_recv_loc + ibs), i, psndrcv_recvs_loc);
      }

      for (i = 0; i < nsends_loc; i++) {
         icpu = psnd2cpu_loc[i];
         isize =
            blksize_loc * (pia_sends_loc[i + 1] - pia_sends_loc[i]) * sizeof (_FltVect);
         ibs = pia_sends_loc[i] * blksize_loc;
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
            int ibs_i = pibsblk_loc[iblkgl] * blksize_loc;
            for (int j = pia_hblk[ipar]; j < pia_hblk[ipar + 1]; j++) {
               int jblk = pja_hblk[j];
               int jhblk = pja2_hblk[j];
               int jblkgl = phblk2blks_loc[jhblk] + jblk;
               if (phblk2cpu_loc[jhblk] == myid) {
                  int ibs_j = pibsblk_loc[jblkgl] * blksize_loc;
                  CMvmSlv_BxB < _Int, _Flt, _FltVect >::MvmA ('+', blksize_loc,
                                                              pasub_loc[j], _x + ibs_j,
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
                  " CMvmParBxBThreads<_Int,_Flt,_FltVect>::MvmA: error: incorrect block number !";
            }
            CVector < _FltVect >::CopyVector (blksize_loc,
                                              px_recv_loc + ind * blksize_loc,
                                              px_temp_loc + jj * blksize_loc);
         }
         for (j = pialist_recvs_loc[i]; j < pialist_recvs_loc[i + 1]; j++) {
            iblk = pjatriples_recvs_loc[j * 3];
            ihblk = pjatriples_recvs_loc[j * 3 + 1];
            ind = pjatriples_recvs_loc[j * 3 + 2];
            iblkgl = phblk2blks_loc[ihblk] + iblk;
            ibs_i = pibsblk_loc[iblkgl] * blksize_loc;
            CMatrix < _Int, _Flt > *pasub_loc = phmatr_loc[ihblk].GetASubArr ();
            CMatrix < int, float >*phmatr_str_loc = phmatr_loc[ihblk].GetHMatrStr ();
            int *pja_hblk = phmatr_str_loc->GetJaArr ();
            int *pja2_hblk = phmatr_str_loc->GetJa2Arr ();
            jj = pja_hblk[ind];
            jj2 = pja2_hblk[ind];
            if (phblk2blks_loc[jj2] + jj != jblkgl) {
               throw
                  " CMvmParBxBThreads<_Int,_Flt,_FltVect>::MvmA: error 2: incorrect block number !";
            }
            CMvmSlv_BxB < _Int, _Flt, _FltVect >::MvmA ('+', blksize_loc, pasub_loc[ind],
                                                        px_temp_loc, _ax + ibs_i);
         }
      }

      CMPIDataExchange::DeleteRecvs (psndrcv_recvs_loc);
      CMPIDataExchange::DeleteStats (psndrcv_stats_loc);

   }

// In-place perform SolveU computations
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CSlvParBxBThreads < _Int, _Flt,
      _FltVect >::SolveU_BxB (int _ihblk, _FltVect * _px)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

// Open mvm structure

      int blksize_loc = this->blksize;
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

      int b_2 = blksize_loc * blksize_loc;

      CVectorData < _FltVect > work_thr (n_thr * b_2);
      _FltVect *pwork_thr = work_thr.Ptr ();

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
               int my_thr = 0;
#ifdef USE_THREADS
               my_thr = omp_get_thread_num ();
#endif
               int iblk = ppnodeslevlist_3[ipar];
               int iblkgl = iblk + nblks12_curr;
               int ibs_i = (int) ppblks_curr[iblkgl] * blksize_loc;
               if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
                  for (int j = pia_hmatr[iblkgl + 1] - 1; j >= pia_hmatr[iblkgl] + 1; j--) {
                     int jblkgl = pja_hmatr[j];
                     int ibs_j = (int) ppblks_curr[jblkgl] * blksize_loc;
                     CMvmSlv_BxB < _Int, _Flt, _FltVect >::MvmA ('-', blksize_loc,
                                                                 pA_sub[j], _px + ibs_j,
                                                                 _px + ibs_i);
                  }
                  int j = pia_hmatr[iblkgl];
                  CMvmSlv_BxB < _Int, _Flt, _FltVect >::SolveU (blksize_loc, pA_sub[j],
                                                                _px + ibs_i, _px + ibs_i,
                                                                pwork_thr + my_thr * b_2);
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
            int my_thr = 0;
#ifdef USE_THREADS
            my_thr = omp_get_thread_num ();
#endif
            int iblkgl;
            if (ipar < nnodes_curr_1) {
               iblkgl = ppnodeslevlist_1[ipar];
            } else {
               iblkgl = nblks1_curr + ppnodeslevlist_2[ipar - nnodes_curr_1];
            }
            int ibs_i = (int) ppblks_curr[iblkgl] * blksize_loc;
            if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
               for (int j = pia_hmatr[iblkgl + 1] - 1; j >= pia_hmatr[iblkgl] + 1; j--) {
                  int jblkgl = pja_hmatr[j];
                  int ibs_j = (int) ppblks_curr[jblkgl] * blksize_loc;
                  CMvmSlv_BxB < _Int, _Flt, _FltVect >::MvmA ('-', blksize_loc, pA_sub[j],
                                                              _px + ibs_j, _px + ibs_i);
               }
               int j = pia_hmatr[iblkgl];
               CMvmSlv_BxB < _Int, _Flt, _FltVect >::SolveU (blksize_loc, pA_sub[j],
                                                             _px + ibs_i, _px + ibs_i,
                                                             pwork_thr + my_thr * b_2);
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
            int my_thr = 0;
#ifdef USE_THREADS
            my_thr = omp_get_thread_num ();
#endif
            int iblkgl = ppnodeslevlist_1[ipar];
            int ibs_i = (int) ppblks_curr[iblkgl] * blksize_loc;
            if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
               for (int j = pia_hmatr[iblkgl + 1] - 1; j >= pia_hmatr[iblkgl] + 1; j--) {
                  int jblkgl = pja_hmatr[j];
                  int ibs_j = (int) ppblks_curr[jblkgl] * blksize_loc;
                  CMvmSlv_BxB < _Int, _Flt, _FltVect >::MvmA ('-', blksize_loc, pA_sub[j],
                                                              _px + ibs_j, _px + ibs_i);
               }
               int j = pia_hmatr[iblkgl];
               CMvmSlv_BxB < _Int, _Flt, _FltVect >::SolveU (blksize_loc, pA_sub[j],
                                                             _px + ibs_i, _px + ibs_i,
                                                             pwork_thr + my_thr * b_2);
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
            int my_thr = 0;
#ifdef USE_THREADS
            my_thr = omp_get_thread_num ();
#endif
            int iblkgl = nblks1_curr + ppnodeslevlist_2[ipar];
            int ibs_i = (int) ppblks_curr[iblkgl] * blksize_loc;
            if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
               for (int j = pia_hmatr[iblkgl + 1] - 1; j >= pia_hmatr[iblkgl] + 1; j--) {
                  int jblkgl = pja_hmatr[j];
                  int ibs_j = (int) ppblks_curr[jblkgl] * blksize_loc;
                  CMvmSlv_BxB < _Int, _Flt, _FltVect >::MvmA ('-', blksize_loc, pA_sub[j],
                                                              _px + ibs_j, _px + ibs_i);
               }
               int j = pia_hmatr[iblkgl];
               CMvmSlv_BxB < _Int, _Flt, _FltVect >::SolveU (blksize_loc, pA_sub[j],
                                                             _px + ibs_i, _px + ibs_i,
                                                             pwork_thr + my_thr * b_2);
            }
         }
      }
   }

// In-place perform SolveL computations
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CSlvParBxBThreads < _Int, _Flt,
      _FltVect >::SolveL_BxB (int _ihblk, _FltVect * _px)
   {

      int n_thr = 1;
#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

// Open mvm structure

      int blksize_loc = this->blksize;
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

      int b_2 = blksize_loc * blksize_loc;

      CVectorData < _FltVect > work_thr (n_thr * b_2);
      _FltVect *pwork_thr = work_thr.Ptr ();

      for (ilev = nlev_1 - 1; ilev >= nlev_12_min; ilev--) {

         nnodes_curr_1 = pnnodes_lev_1[ilev];
         ppnodeslevlist_1 = &pnodeslevlist_1[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_1; ipar++) {
            int my_thr = 0;
#ifdef USE_THREADS
            my_thr = omp_get_thread_num ();
#endif
            int iblkgl = ppnodeslevlist_1[ipar];
            int ibs_i = (int) ppblks_curr[iblkgl] * blksize_loc;
            if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
               for (int j = pia_hmatr[iblkgl]; j < pia_hmatr[iblkgl + 1] - 1; j++) {
                  int jblkgl = pja_hmatr[j];
                  int ind = pind_lt2l[j];
                  int ibs_j = (int) ppblks_curr[jblkgl] * blksize_loc;
                  CMvmSlv_BxB < _Int, _Flt, _FltVect >::MvmAT ('-', blksize_loc,
                                                               pA_sub[ind], _px + ibs_j,
                                                               _px + ibs_i);
               }
               int j = pia_hmatr[iblkgl + 1] - 1;
               int ind = pind_lt2l[j];
               CMvmSlv_BxB < _Int, _Flt, _FltVect >::SolveL (blksize_loc, pA_sub[ind],
                                                             _px + ibs_i, _px + ibs_i,
                                                             pwork_thr + my_thr * b_2);
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
            int my_thr = 0;
#ifdef USE_THREADS
            my_thr = omp_get_thread_num ();
#endif
            int iblkgl = nblks1_curr + ppnodeslevlist_2[ipar];
            int ibs_i = (int) ppblks_curr[iblkgl] * blksize_loc;
            if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
               for (int j = pia_hmatr[iblkgl]; j < pia_hmatr[iblkgl + 1] - 1; j++) {
                  int jblkgl = pja_hmatr[j];
                  int ind = pind_lt2l[j];
                  int ibs_j = (int) ppblks_curr[jblkgl] * blksize_loc;
                  CMvmSlv_BxB < _Int, _Flt, _FltVect >::MvmAT ('-', blksize_loc,
                                                               pA_sub[ind], _px + ibs_j,
                                                               _px + ibs_i);
               }
               int j = pia_hmatr[iblkgl + 1] - 1;
               int ind = pind_lt2l[j];
               CMvmSlv_BxB < _Int, _Flt, _FltVect >::SolveL (blksize_loc, pA_sub[ind],
                                                             _px + ibs_i, _px + ibs_i,
                                                             pwork_thr + my_thr * b_2);
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
            int my_thr = 0;
#ifdef USE_THREADS
            my_thr = omp_get_thread_num ();
#endif
            int iblkgl;
            if (ipar < nnodes_curr_1) {
               iblkgl = ppnodeslevlist_1[ipar];
            } else {
               iblkgl = nblks1_curr + ppnodeslevlist_2[ipar - nnodes_curr_1];
            }
            int ibs_i = (int) ppblks_curr[iblkgl] * blksize_loc;
            if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
               for (int j = pia_hmatr[iblkgl]; j < pia_hmatr[iblkgl + 1] - 1; j++) {
                  int jblkgl = pja_hmatr[j];
                  int ind = pind_lt2l[j];
                  int ibs_j = (int) ppblks_curr[jblkgl] * blksize_loc;
                  CMvmSlv_BxB < _Int, _Flt, _FltVect >::MvmAT ('-', blksize_loc,
                                                               pA_sub[ind], _px + ibs_j,
                                                               _px + ibs_i);
               }
               int j = pia_hmatr[iblkgl + 1] - 1;
               int ind = pind_lt2l[j];
               CMvmSlv_BxB < _Int, _Flt, _FltVect >::SolveL (blksize_loc, pA_sub[ind],
                                                             _px + ibs_i, _px + ibs_i,
                                                             pwork_thr + my_thr * b_2);
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
               int my_thr = 0;
#ifdef USE_THREADS
               my_thr = omp_get_thread_num ();
#endif
               int iblk = ppnodeslevlist_3[ipar];
               int iblkgl = iblk + nblks12_curr;
               int ibs_i = (int) ppblks_curr[iblkgl] * blksize_loc;
               if (pia_hmatr[iblkgl + 1] > pia_hmatr[iblkgl]) {
                  for (int j = pia_hmatr[iblkgl]; j < pia_hmatr[iblkgl + 1] - 1; j++) {
                     int jblkgl = pja_hmatr[j];
                     int ind = pind_lt2l[j];
                     int ibs_j = (int) ppblks_curr[jblkgl] * blksize_loc;
                     CMvmSlv_BxB < _Int, _Flt, _FltVect >::MvmAT ('-', blksize_loc,
                                                                  pA_sub[ind],
                                                                  _px + ibs_j,
                                                                  _px + ibs_i);
                  }
                  int j = pia_hmatr[iblkgl + 1] - 1;
                  int ind = pind_lt2l[j];
                  CMvmSlv_BxB < _Int, _Flt, _FltVect >::SolveL (blksize_loc, pA_sub[ind],
                                                                _px + ibs_i, _px + ibs_i,
                                                                pwork_thr + my_thr * b_2);
               }
            }

         }

      }

   }

// Perform SolveLU computations
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CSlvParBxBThreads < _Int, _Flt,
      _FltVect >::SolveLU_BxB (const _FltVect * _x, _FltVect * _px)
   {

// Open mvm structure

      void *pcomm_loc = this->pcomm;

      int blksize_loc = this->blksize;
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

      CVector < _FltVect >::SetByZeroes_thr (ni_cpu_loc * blksize_loc, _px);

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
               CVector < _FltVect >::CopyVector (blksize_loc, _x + ind * blksize_loc,
                                                 px_send_loc + i * blksize_loc);
            }
         }
      }

      int ni_recv_loc = pia_recvs_loc[nrecvs_loc];

      CVector < _FltVect >::SetByZeroes_thr (ni_recv_loc * blksize_loc, px_recv_loc);

// Init async recvs and sends

      void *psndrcv_recvs_loc;
      void *psndrcv_stats_loc;

      CMPIDataExchange::AllocateRecvs (nrecvs_loc + nsends_loc, psndrcv_recvs_loc);
      CMPIDataExchange::AllocateStats (nrecvs_loc + nsends_loc, psndrcv_stats_loc);

      int icpu, isize, ibs;

      for (i = 0; i < nrecvs_loc; i++) {
         icpu = prcv2cpu_loc[i];
         isize =
            blksize_loc * (pia_recvs_loc[i + 1] - pia_recvs_loc[i]) * sizeof (_FltVect);
         ibs = pia_recvs_loc[i] * blksize_loc;
         CMPIDataExchange::IRecv (pcomm_loc, icpu, myid, isize,
                                  (char *) (px_recv_loc + ibs), i, psndrcv_recvs_loc);
      }

      for (i = 0; i < nsends_loc; i++) {
         icpu = psnd2cpu_loc[i];
         isize =
            blksize_loc * (pia_sends_loc[i + 1] - pia_sends_loc[i]) * sizeof (_FltVect);
         ibs = pia_sends_loc[i] * blksize_loc;
         CMPIDataExchange::ISend (pcomm_loc, icpu, icpu, isize,
                                  (char *) (px_send_loc + ibs), i + nrecvs_loc,
                                  psndrcv_recvs_loc);
      }

// Wait for completetion of sends/recvs

      CMPIDataExchange::WaitAll (nrecvs_loc + nsends_loc, psndrcv_recvs_loc,
                                 psndrcv_stats_loc);

      CVector < _FltVect >::SetByZeroes_thr (ni_send_loc * blksize_loc, px_send_loc);

// Perform local computations

      int iblk0, ihblk, ibs_i, niloc, niextloc, ni_ini;

      _FltVect *px1_temp;
      _FltVect *px2_temp;

      for (i = 0; i < nlisthblk_own_loc; i++) {
         ihblk = plisthblk_own_loc[i];
         iblk0 = phblk2blks_loc[ihblk];
         ibs_i = pibsblk_loc[iblk0] * blksize_loc;
         niloc = (int) (phblks_loc[ihblk + 1] - phblks_loc[ihblk]);
         niextloc = (int) (phblks_ext_loc[ihblk + 1] - phblks_ext_loc[ihblk]);
         ni_ini = niextloc - niloc;
         px1_temp = px_temp_loc;
         px2_temp = px1_temp + niextloc * blksize_loc;
         CVector < _FltVect >::SetByZeroes_thr (ni_ini * blksize_loc, px1_temp);
         CVector < _FltVect >::CopyVector_thr (niloc * blksize_loc, _x + ibs_i,
                                               px1_temp + ni_ini * blksize_loc);
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
                     int ibs_j = pibsblk_loc[jj2] * blksize_loc;
                     CVector < _FltVect >::CopyVector (blksize_loc,
                                                       _x + ibs_j + jj * blksize_loc,
                                                       px1_temp + j * blksize_loc);
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
                  CVector < _FltVect >::CopyVector (blksize_loc,
                                                    px_recv_loc + ind * blksize_loc,
                                                    px1_temp + ind1 * blksize_loc);
               }
            }
         }
         {
            int *pporderLU = this->porderLU[ihblk].Ptr ();
            CVector < _FltVect >::OrderVector_thr (blksize_loc, niextloc, pporderLU,
                                                   px1_temp, px2_temp);
            this->SolveL_BxB (ihblk, px2_temp);
            CVector < _FltVect >::SetByZeroes_thr (ni_ini * blksize_loc, px2_temp);
            this->SolveU_BxB (ihblk, px2_temp);
            CVector < _FltVect >::InvOrderVector_thr (blksize_loc, niextloc, pporderLU,
                                                      px2_temp, px1_temp);
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
                  CVector < _FltVect >::CopyVector (blksize_loc,
                                                    px1_temp + ind1 * blksize_loc,
                                                    px_recv_loc + ind * blksize_loc);
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
                     CVector < _FltVect >::AddReplaceVector (blksize_loc,
                                                             px1_temp + j * blksize_loc,
                                                             _px + (ibs_j +
                                                                    jj) * blksize_loc);
                  }
               }
            }
         }
         CVector < _FltVect >::AddReplaceVector_thr (niloc * blksize_loc,
                                                     px1_temp + ni_ini * blksize_loc,
                                                     _px + ibs_i);
      }

// Backward exchanges

      for (i = 0; i < nsends_loc; i++) {
         icpu = psnd2cpu_loc[i];
         isize =
            blksize_loc * (pia_sends_loc[i + 1] - pia_sends_loc[i]) * sizeof (_FltVect);
         ibs = pia_sends_loc[i] * blksize_loc;
         CMPIDataExchange::IRecv (pcomm_loc, icpu, myid, isize,
                                  (char *) (px_send_loc + ibs), i + nrecvs_loc,
                                  psndrcv_recvs_loc);
      }

      for (i = 0; i < nrecvs_loc; i++) {
         icpu = prcv2cpu_loc[i];
         isize =
            blksize_loc * (pia_recvs_loc[i + 1] - pia_recvs_loc[i]) * sizeof (_FltVect);
         ibs = pia_recvs_loc[i] * blksize_loc;
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
               CVector < _FltVect >::AddReplaceVector (blksize_loc,
                                                       px_send_loc + i * blksize_loc,
                                                       _px + ind * blksize_loc);
            }
         }
      }

      CMPIDataExchange::DeleteRecvs (psndrcv_recvs_loc);
      CMPIDataExchange::DeleteStats (psndrcv_stats_loc);

   }

// Reorder rhs and solution from/to initial ordering
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::ReorderVectorDataIni_BxB (char _dir, _FltVect * _rhs, _FltVect * _sol)
   {

// Open data

      int blksize_loc = this->blksize;
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

      CVectorData < _FltVect > ztemp (nimax * blksize_loc);
      _FltVect *pztemp = ztemp.Ptr ();

      int ibs = 0;

      for (i = 0; i < nhblks_loc; i++) {
         niloc = (int) (phblks[i + 1] - phblks[i]);
         if (phblk2cpu[i] == myid_loc) {
            int *pord_ini = porder_ini[i].Ptr ();
            if (_dir == 'D') {
               CVector < _FltVect >::OrderVector_thr (blksize_loc, niloc, pord_ini,
                                                      _rhs + ibs, pztemp);
               CVector < _FltVect >::CopyVector_thr (niloc * blksize_loc, pztemp,
                                                     _rhs + ibs);
               CVector < _FltVect >::OrderVector_thr (blksize_loc, niloc, pord_ini,
                                                      _sol + ibs, pztemp);
               CVector < _FltVect >::CopyVector_thr (niloc * blksize_loc, pztemp,
                                                     _sol + ibs);
            } else {
               CVector < _FltVect >::InvOrderVector_thr (blksize_loc, niloc, pord_ini,
                                                         _rhs + ibs, pztemp);
               CVector < _FltVect >::CopyVector_thr (niloc * blksize_loc, pztemp,
                                                     _rhs + ibs);
               CVector < _FltVect >::InvOrderVector_thr (blksize_loc, niloc, pord_ini,
                                                         _sol + ibs, pztemp);
               CVector < _FltVect >::CopyVector_thr (niloc * blksize_loc, pztemp,
                                                     _sol + ibs);
            }
            ibs += niloc * blksize_loc;
         }
      }

   }

// Reorder rhs and solution from/to wells ordering
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::ReorderVectorDataWells_BxB (char _dir, _FltVect * _rhs, _FltVect * _sol,
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

      int blksize_loc = this->blksize;
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
         _rhs_ord.resize ((int) (nz_ord * blksize));
         _sol_ord.resize ((int) (nz_ord * blksize));
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
                                                      nii_curr * 2 * blksize_loc);
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
                     CVector < _FltVect >::CopyVector (blksize, _rhs + ind * blksize_loc,
                                                       pa_temp + j * 2 * blksize_loc);
                     CVector < _FltVect >::CopyVector (blksize, _sol + ind * blksize_loc,
                                                       pa_temp + (j * 2 +
                                                                  1) * blksize_loc);
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
                  CVector < _FltVect >::CopyVector (blksize,
                                                    pa_temp + i * 2 * blksize_loc,
                                                    p_rhs_ord + (ibs + jj) * blksize_loc);
                  CVector < _FltVect >::CopyVector (blksize,
                                                    pa_temp + (i * 2 + 1) * blksize_loc,
                                                    p_sol_ord + (ibs + jj) * blksize_loc);
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
               pASub[ipar].ResizeA (nlist_temp * blksize_loc);
               pASub[ipar].SetNza (nlist_temp * blksize_loc);
               _FltVect *pa_temp = pASub[ipar].GetAArr ();
               for (i = 0; i < nlist_temp; i++) {
                  jj = plist2_temp[i * 3 + 1];
                  jhblk = plist2_temp[i * 3 + 2];
                  ibs = (int) pibs_ord[jhblk];
                  CVector < _FltVect >::CopyVector (blksize,
                                                    p_sol_ord + (ibs + jj) * blksize_loc,
                                                    pa_temp + i * blksize_loc);
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
                  CVector < _FltVect >::CopyVector (blksize, pa_temp + i * blksize_loc,
                                                    _sol + ind * blksize_loc);
               }
            }

         }

      }

   }

// Perform SlvLU parallel computations
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::SlvLU_BxB (const _FltVect * _x, _FltVect * _px)
   {
      if (!this->b_use_blksize) {
         if (this->params.prec_float == 1) {
            this->slv_float_bxb.SolveLU_BxB (_x, _px);
         } else {
            this->slv_double_bxb.SolveLU_BxB (_x, _px);
         }
      } else {

         int ni_local = 0;

         if (this->params.prec_float == 1) {
            ni_local = this->slv_float_bxb.GetNiCpu ();
         } else {
            ni_local = this->slv_double_bxb.GetNiCpu ();
         }

         int isize = (int) (this->xwork_bscl.GetLength ());
         int isize_work = 2 * ni_local + 1;

         if (isize != isize_work) {
            this->xwork_bscl.resize (isize_work);
         }

         _FltVect *pxwork_f = this->xwork_bscl.Ptr ();
         _FltVect *pxwork1_f = pxwork_f + ni_local;

         this->MvmBScl ('L', 'T', _x, pxwork_f);

         if (this->params.prec_float == 1) {
            this->slv_float_bxb.SolveLU_BxB (pxwork_f, pxwork1_f);
         } else {
            this->slv_double_bxb.SolveLU_BxB (pxwork_f, pxwork1_f);
         }

         this->MvmBScl ('U', 'N', pxwork1_f, _px);

      }
   }

// Perform mvm by part of block scaling
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::MvmBScl (char _type, char _transp, const _FltVect * _x, _FltVect * _px)
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

      int blksize_loc = this->blksize;

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
                  int ishift =
                     (int) (pblks[ibegblk + ipar] - pblks[ibegblk]) * blksize_loc;

                  vector < _Flt > *pa_matr = pASub[ipar].GetA ();

                  if (_transp == 'N') {

                     CMvmSlv_impl < _Int, _Flt, _FltVect >::MvmDiagBlksize (blksize_loc,
                                                                            nlistloc,
                                                                            *pa_matr,
                                                                            _x +
                                                                            ibs_vect +
                                                                            ishift,
                                                                            _px +
                                                                            ibs_vect +
                                                                            ishift);

                  } else {

                     CMvmSlv_impl < _Int, _Flt, _FltVect >::MvmDiagTBlksize (blksize_loc,
                                                                             nlistloc,
                                                                             *pa_matr,
                                                                             _x +
                                                                             ibs_vect +
                                                                             ishift,
                                                                             _px +
                                                                             ibs_vect +
                                                                             ishift);

                  }

               }

            }

            ibs_vect += (int) (phblks[ihblk + 1] - phblks[ihblk]) * blksize_loc;

         }
      }

   }

// Transform vectors forward
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::TransformVectorsForward_BxB (_FltVect * _rhs_ini, _FltVect * _sol_ini,
                                               vector < _FltVect > &_rhs_new,
                                               vector < _FltVect > &_sol_new)
   {

      int blksize_loc = this->blksize;

      long long isize_curr = this->GetVSize (1);

      _rhs_new.resize ((int) isize_curr * blksize_loc + 1);
      _sol_new.resize ((int) isize_curr * blksize_loc + 1);

      _FltVect *p_rhs_new = &_rhs_new[0];
      _FltVect *p_sol_new = &_sol_new[0];

      if (!this->b_use_wells) {

         if (this->b_use_ini) {
            CVector < _FltVect >::CopyVector_thr ((int) isize_curr * blksize_loc,
                                                  _rhs_ini, p_rhs_new);
            CVector < _FltVect >::CopyVector_thr ((int) isize_curr * blksize_loc,
                                                  _sol_ini, p_sol_new);
            this->ReorderVectorDataIni_BxB ('D', p_rhs_new, p_sol_new);
         } else {
            CVector < _FltVect >::CopyVector_thr ((int) isize_curr * blksize_loc,
                                                  _rhs_ini, p_rhs_new);
            CVector < _FltVect >::CopyVector_thr ((int) isize_curr * blksize_loc,
                                                  _sol_ini, p_sol_new);
         }

      } else {

         CVectorData < _FltVect > rhs_ord;
         CVectorData < _FltVect > sol_ord;

         this->ReorderVectorDataWells_BxB ('D', _rhs_ini, _sol_ini, rhs_ord, sol_ord);

         _FltVect *prhs_ord = rhs_ord.Ptr ();
         _FltVect *psol_ord = sol_ord.Ptr ();

         CVector < _FltVect >::CopyVector_thr ((int) isize_curr * blksize_loc, prhs_ord,
                                               p_rhs_new);
         CVector < _FltVect >::CopyVector_thr ((int) isize_curr * blksize_loc, psol_ord,
                                               p_sol_new);

      }

   }

// Transform vectors backward
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::TransformVectorsBackward_BxB (vector < _FltVect > &_rhs_new,
                                                vector < _FltVect > &_sol_new,
                                                _FltVect * _rhs_fin, _FltVect * _sol_fin)
   {

      int blksize_loc = this->blksize;
      long long isize_curr = this->GetVSize (1);

      _FltVect *p_rhs_new = &_rhs_new[0];
      _FltVect *p_sol_new = &_sol_new[0];

      if (!this->b_use_wells) {

         CVector < _FltVect >::CopyVector_thr ((int) isize_curr * blksize_loc, p_rhs_new,
                                               _rhs_fin);
         CVector < _FltVect >::CopyVector_thr ((int) isize_curr * blksize_loc, p_sol_new,
                                               _sol_fin);

         if (this->b_use_ini) {
            this->ReorderVectorDataIni_BxB ('I', _rhs_fin, _sol_fin);
         }

      } else {

         CVectorData < _FltVect > rhs_ord ((int) isize_curr * blksize_loc);
         CVectorData < _FltVect > sol_ord ((int) isize_curr * blksize_loc);

         _FltVect *prhs_ord = rhs_ord.Ptr ();
         _FltVect *psol_ord = sol_ord.Ptr ();

         CVector < _FltVect >::CopyVector_thr ((int) isize_curr * blksize_loc, p_rhs_new,
                                               prhs_ord);
         CVector < _FltVect >::CopyVector_thr ((int) isize_curr * blksize_loc, p_sol_new,
                                               psol_ord);

         this->ReorderVectorDataWells_BxB ('I', _rhs_fin, _sol_fin, rhs_ord, sol_ord);

      }

   }

// Perform iterations of the iterative scheme
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::SolveIter_BxB (SParams * _params, SStatData * _stats, _FltVect * _rhs,
                                 _FltVect * _sol)
   {

      this->SolveIter_BxB (_params, _stats, this, CK3D_SolverBxBThreads < _Int, _Flt,
                           _FltVect >::MvmA_BxB_static, this,
                           CK3D_SolverBxBThreads < _Int, _Flt,
                           _FltVect >::SlvLU_BxB_static, _rhs, _sol);

   }

// Perform iterations of the iterative scheme
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::SolveIter_BxB (SParams * _params, SStatData * _stats, void *_str_mvmA,
                                 typename TMvmSlvFunc < _FltVect >::mvmA_f _mvm_f,
                                 void *_str_slvLU,
                                 typename TMvmSlvFunc < _FltVect >::slvLU_f _slv_f,
                                 _FltVect * _rhs, _FltVect * _sol)
   {

      vector < _FltVect > rhs_new;
      vector < _FltVect > sol_new;

      this->TransformVectorsForward_BxB (_rhs, _sol, rhs_new, sol_new);

      _FltVect *prhs_new = rhs_new.data ();
      _FltVect *psol_new = sol_new.data ();

      this->SolveIter_BxB_NoTransform (_params, _stats, this,
                                       CK3D_SolverBxBThreads < _Int, _Flt,
                                       _FltVect >::MvmA_BxB_static, this,
                                       CK3D_SolverBxBThreads < _Int, _Flt,
                                       _FltVect >::SlvLU_BxB_static, prhs_new, psol_new);

      this->TransformVectorsBackward_BxB (rhs_new, sol_new, _rhs, _sol);

   }

// Perform iterations of the iterative scheme (implementation) (no vector data transformation)
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::SolveIter_BxB_NoTransform (SParams * _params, SStatData * _stats,
                                             void *_str_mvmA,
                                             typename TMvmSlvFunc <
                                             _FltVect >::mvmA_f _mvm_f, void *_str_slvLU,
                                             typename TMvmSlvFunc <
                                             _FltVect >::slvLU_f _slv_f, _FltVect * _rhs,
                                             _FltVect * _sol)
   {

      int blksize_loc = this->blksize;

      _params->blksize_iter = blksize_loc;

      this->SolveIter (_params, _stats, _str_mvmA, _mvm_f, _str_slvLU, _slv_f, _rhs,
                       _sol);

   }

// Prepare solver structures including performing parallel fct
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::ComputeBILU2_BxB (SParams * _params, SStatData * _stats)
   {

      void *pcommloc = this->pcomm;

      int myid = CMPIDataExchange::GetMyid (pcommloc);
      int nproc = CMPIDataExchange::GetNproc (pcommloc);

      int sctype_expl_loc = _params->sctype_expl;

      if (sctype_expl_loc == -1) {

         this->ComputeBILU2_BxB_impl (_params, _stats);
         this->b_use_blksize = false;

      } else {

// Get control data

         int blksizeloc = this->blksize;
         int nhblksloc = this->nhblks;
         int nblksloc = this->nblks;

         long long *phblks = this->hblks.data ();
         int *phblk2cpu = this->hblk2cpu.data ();
         int *phblk2blks = this->hblk2blks.data ();
         int *pblk2hblks = this->blk2hblks.data ();
         long long *pblks = this->blks.data ();

// Save copy of original submatrix

         CBMatrix < _Int, _Flt > *phmatr_arr = this->hmatr_arr.data ();

         this->hmatr_arr_save.resize (nhblksloc);

         CBMatrix < _Int, _Flt > *phmatr_arr_save = this->hmatr_arr_save.data ();

         int i;

         for (i = 0; i < nhblksloc; i++) {
            if (phblk2cpu[i] == myid) {

               CBMatrix < _Int, _Flt > hblk_temp;

               hblk_temp.AddReplaceBxBThr (true, blksizeloc, phmatr_arr[i]);

               phmatr_arr_save[i].ReplaceFree (hblk_temp);

            }
         }

// Repartition the whole matrix onto the first cpu, 1 hblock

         int nhblks_1cpu = 1;
         int hblk2cpu_1cpu = 0;
         int hblk2blks_1cpu[2];

         hblk2blks_1cpu[0] = 0;
         hblk2blks_1cpu[1] = nblksloc;

         vector < CBMatrix < _Int, _Flt > >hmatr_arr_1cpu (nhblks_1cpu);
         CBMatrix < _Int, _Flt > *phmatr_arr_1cpu = hmatr_arr_1cpu.data ();

         int ni_ord = 0;

         for (i = 0; i < nhblksloc; i++) {
            if (phblk2cpu[i] == myid) {
               int ibeg_blk = phblk2blks[i];
               int iend_blk = phblk2blks[i + 1] - 1;
               int niloc = (int) (pblks[iend_blk + 1] - pblks[ibeg_blk]);
               ni_ord += niloc;
            }
         }

         vector < long long >order_1cpu (ni_ord);
         long long *porder_1cpu = order_1cpu.data ();

         ni_ord = 0;

         for (i = 0; i < nhblksloc; i++) {
            if (phblk2cpu[i] == myid) {
               int ibeg_blk = phblk2blks[i];
               int iend_blk = phblk2blks[i + 1] - 1;
               int niloc = (int) (pblks[iend_blk + 1] - pblks[ibeg_blk]);
               long long ibeg = pblks[ibeg_blk];
               CVectorInt < long long >::SetIdentity (niloc, ibeg, porder_1cpu + ni_ord);
               ni_ord += niloc;
            }
         }

         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (pcommloc, blksizeloc, nhblksloc,
                                                      phblk2cpu, phblk2blks, pblks,
                                                      phmatr_arr, porder_1cpu,
                                                      nhblks_1cpu, &hblk2cpu_1cpu,
                                                      hblk2blks_1cpu, pblks,
                                                      phmatr_arr_1cpu);

// Compute explicit block scaling on first cpu

         int sctype_expl = _params->sctype_expl;
         double sclmin = _params->sclmin;
         int nitersc_expl = _params->nitersc_expl;

         CVectorData < _Flt > scl_L_vect;
         CVectorData < _Flt > scl_U_vect;
         CVectorData < _Flt > scl_L_inv_vect;
         CVectorData < _Flt > scl_U_inv_vect;

         if (myid == 0) {

            int nmodif_scl;
            double sclmin_att, sclmax_att;

            CFctBxBThreads < _Int, _Flt >::ComputeExplicitScaling_BxB (sctype_expl,
                                                                       nitersc_expl,
                                                                       sclmin, blksizeloc,
                                                                       nblksloc, pblks,
                                                                       *phmatr_arr_1cpu,
                                                                       scl_L_vect,
                                                                       scl_U_vect,
                                                                       scl_L_inv_vect,
                                                                       scl_U_inv_vect,
                                                                       nmodif_scl,
                                                                       sclmin_att,
                                                                       sclmax_att);

         }

         _Flt *pscl_L_vect = scl_L_vect.Ptr ();
         _Flt *pscl_U_vect = scl_U_vect.Ptr ();
         _Flt *pscl_L_inv_vect = scl_L_inv_vect.Ptr ();
         _Flt *pscl_U_inv_vect = scl_U_inv_vect.Ptr ();

// Scale explicitely

         if (myid == 0) {

            CFctBxBThreads < _Int, _Flt >::HMatrixScale_BxB (blksizeloc, nblksloc, pblks,
                                                             *phmatr_arr_1cpu,
                                                             pscl_L_vect, pscl_U_vect);

         }
// Create block diagonal scaling matrices

         CBMatrix < _Int, _Flt > DL_hmatr;
         CBMatrix < _Int, _Flt > DU_hmatr;
         CBMatrix < _Int, _Flt > DLInv_hmatr;
         CBMatrix < _Int, _Flt > DUInv_hmatr;

         if (myid == 0) {

            CBMatrix < _Int, _Flt >::BlockDiagonalHMatrix_BxB (blksizeloc, nblksloc,
                                                               pblks, pscl_L_vect,
                                                               DL_hmatr);
            CBMatrix < _Int, _Flt >::BlockDiagonalHMatrix_BxB (blksizeloc, nblksloc,
                                                               pblks, pscl_U_vect,
                                                               DU_hmatr);
            CBMatrix < _Int, _Flt >::BlockDiagonalHMatrix_BxB (blksizeloc, nblksloc,
                                                               pblks, pscl_L_inv_vect,
                                                               DLInv_hmatr);
            CBMatrix < _Int, _Flt >::BlockDiagonalHMatrix_BxB (blksizeloc, nblksloc,
                                                               pblks, pscl_U_inv_vect,
                                                               DUInv_hmatr);

         }
// Repartition matrix and explicit scaling back to original partitioning

         long long nsup_tot = pblks[nblksloc];

         vector < long long >order_back (nsup_tot);
         long long *porder_back = order_back.data ();

         CVectorInt < long long >::SetIdentity ((int) nsup_tot, 0, porder_back);

         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (pcommloc, blksizeloc, nhblks_1cpu,
                                                      &hblk2cpu_1cpu, hblk2blks_1cpu,
                                                      pblks, phmatr_arr_1cpu, porder_back,
                                                      nhblksloc, phblk2cpu, phblk2blks,
                                                      pblks, phmatr_arr);

         vector < CBMatrix < _Int, _Flt > >sclL_temp (nhblksloc);
         vector < CBMatrix < _Int, _Flt > >sclU_temp (nhblksloc);
         vector < CBMatrix < _Int, _Flt > >sclLInv_temp (nhblksloc);
         vector < CBMatrix < _Int, _Flt > >sclUInv_temp (nhblksloc);

         CBMatrix < _Int, _Flt > *psclL_arr = sclL_temp.data ();
         CBMatrix < _Int, _Flt > *psclU_arr = sclU_temp.data ();
         CBMatrix < _Int, _Flt > *psclLInv_arr = sclLInv_temp.data ();
         CBMatrix < _Int, _Flt > *psclUInv_arr = sclUInv_temp.data ();

         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (pcommloc, blksizeloc, nhblks_1cpu,
                                                      &hblk2cpu_1cpu, hblk2blks_1cpu,
                                                      pblks, &DL_hmatr, porder_back,
                                                      nhblksloc, phblk2cpu, phblk2blks,
                                                      pblks, psclL_arr);

         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (pcommloc, blksizeloc, nhblks_1cpu,
                                                      &hblk2cpu_1cpu, hblk2blks_1cpu,
                                                      pblks, &DU_hmatr, porder_back,
                                                      nhblksloc, phblk2cpu, phblk2blks,
                                                      pblks, psclU_arr);

         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (pcommloc, blksizeloc, nhblks_1cpu,
                                                      &hblk2cpu_1cpu, hblk2blks_1cpu,
                                                      pblks, &DLInv_hmatr, porder_back,
                                                      nhblksloc, phblk2cpu, phblk2blks,
                                                      pblks, psclLInv_arr);

         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (pcommloc, blksizeloc, nhblks_1cpu,
                                                      &hblk2cpu_1cpu, hblk2blks_1cpu,
                                                      pblks, &DUInv_hmatr, porder_back,
                                                      nhblksloc, phblk2cpu, phblk2blks,
                                                      pblks, psclUInv_arr);

// Compute BIlu2

         this->ComputeBILU2_BxB_impl (_params, _stats);

// Restore scaling

         this->sclL.swap (sclL_temp);
         this->sclU.swap (sclU_temp);
         this->sclLInv.swap (sclLInv_temp);
         this->sclUInv.swap (sclUInv_temp);

// Restore original submatrix if not two level iterations

         if (_params->it_external == -1) {

            this->hmatr_arr.resize (nhblksloc);

            phmatr_arr = this->hmatr_arr.data ();

            for (i = 0; i < nhblksloc; i++) {
               if (phblk2cpu[i] == myid) {

                  phmatr_arr[i].ReplaceFree (phmatr_arr_save[i]);

                  CBMatrix < _Int, _Flt > hblk_temp;

                  phmatr_arr_save[i].ReplaceFree (hblk_temp);

               }
            }

            this->b_use_blksize = true;

         }

      }

   }

// Prepare solver structures including performing parallel fct
//========================================================================================
   template < typename _Int, typename _Flt,
      typename _FltVect > void CK3D_SolverBxBThreads < _Int, _Flt,
      _FltVect >::ComputeBILU2_BxB_impl (SParams * _params, SStatData * _stats)
   {

// Init output data

      _stats->prec_extend = 1.0;
      _stats->density = 0.0e0;
      _stats->scpiv_min = 1.0e100;
      _stats->scpiv_max = -1.0e100;
      _stats->ndiasplit = 0;
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

      int blksizeloc = this->blksize;
      int nhblksloc = this->nhblks;
      int nblksloc = this->nblks;

      this->params = *_params;

      long long *phblks = &this->hblks[0];
      int *phblk2cpu = &this->hblk2cpu[0];
      int *phblk2blks = &this->hblk2blks[0];
      int *pblk2hblks = &this->blk2hblks[0];
      long long *pblks = &this->blks[0];

      int myid = CMPIDataExchange::GetMyid (pcommloc);
      int nproc = CMPIDataExchange::GetNproc (pcommloc);

      _params->iparam1 = myid;

// Perform temporary explicit block scaling

      int sctype_save = _params->sctype;
//      double scpiv_min_save = 1.0e100;
//      double scpiv_max_save = -1.0e100;

//      char strbuff_debug[256];
//      sprintf (strbuff_debug,"ChkCompBIlu_%i.dat",myid);
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

      CBMatrix < _Int, _Flt >::GetExtendedSubmatrices_BxB (pcommloc, blksizeloc,
                                                           nhblksloc, phblk2cpu,
                                                           phblk2blks, pblk2hblks,
                                                           nblksloc, pblks, phmatr_arr,
                                                           pnblks_ext_arr,
                                                           pblksnum_ext_arr,
                                                           pblks_ext_arr, pnlist_ext_arr,
                                                           plist_ext_arr, phmatr_ext_arr);

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

            CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (blksizeloc, pnblks_ext_arr[i],
                                                         ppblks_ext_arr,
                                                         phmatr_ext_arr[i], pporder_arr,
                                                         pnblks_ilu2_arr[i],
                                                         ppblks_new_arr, hmatr_ord);

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

// Prepare statistics data for each cpu

      vector < int >ni2cpu (nproc);
      vector < int >niext2cpu (nproc);
      vector < int >nsplit2cpu (nproc);
      vector < int >nmodif2cpu (nproc);
      vector < long long >nza2cpu (nproc);
      vector < long long >nzlu2cpu (nproc);
      vector < double >dns2cpu (nproc);
      vector < double >dtime2cpu (nproc);

      int *pni2cpu = ni2cpu.data ();
      int *pniext2cpu = niext2cpu.data ();
      int *pnsplit2cpu = nsplit2cpu.data ();
      int *pnmodif2cpu = nmodif2cpu.data ();
      long long *pnza2cpu = nza2cpu.data ();
      long long *pnzlu2cpu = nzlu2cpu.data ();
      double *pdns2cpu = dns2cpu.data ();
      double *pdtime2cpu = dtime2cpu.data ();

      for (i = 0; i < nproc; i++)
         pni2cpu[i] = 0;
      for (i = 0; i < nproc; i++)
         pniext2cpu[i] = 0;
      for (i = 0; i < nproc; i++)
         pnsplit2cpu[i] = 0;
      for (i = 0; i < nproc; i++)
         pnmodif2cpu[i] = 0;
      for (i = 0; i < nproc; i++)
         pnza2cpu[i] = 0;
      for (i = 0; i < nproc; i++)
         pnzlu2cpu[i] = 0;
      for (i = 0; i < nproc; i++)
         pdns2cpu[i] = 0.0e0;
      for (i = 0; i < nproc; i++)
         pdtime2cpu[i] = 0.0e0;

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

      _stats->nmodif_scl = 0;
      _stats->ndiasplit = 0;
      _stats->nmodif = 0;

      int nmodif_scl_loc = 0;
      int ndiasplit_loc = 0;
      int nmodif_loc = 0;

      for (i = 0; i < nhblksloc; i++) {
         if (phblk2cpu[i] == myid) {

            int nblks_ilu2 = 0;
            vector < long long >blks_ilu2 (1);
            CVectorData < int >order_ilu2 (1);

            {

               long long *ppblks_new_arr = &pblks_ilu2_arr[i][0];

               pni2cpu[myid] +=
                  (int) ((ppblks_new_arr[pnblks_ilu2_arr[i]] -
                          ppblks_new_arr[pnblks1_ilu2_arr[i]]));
               pniext2cpu[myid] += (int) (ppblks_new_arr[pnblks1_ilu2_arr[i]]);

            }

            double time2;

            time2 = CMPIDataExchange::GetWallTimeMPI ();

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
                     " Nsup_Loc = " << ppblks_new_arr[pnblks_ilu2_arr[i]] -
                     ppblks_new_arr[pnblks1_ilu2_arr[i]] << " Nsup_Ext = " <<
                     ppblks_new_arr[pnblks1_ilu2_arr[i]] << endl;
               }
               if (_params->pfout != NULL && _params->msglev >= 1) {
                  *_params->
                     pfout << " Fct Ihblk: " << i << " Nblks = " << pnblks_ilu2_arr[i] -
                     pnblks1_ilu2_arr[i] << " Nblks_ext = " << pnblks_ilu2_arr[i] <<
                     " Ni_Loc = " << ppblks_new_arr[pnblks_ilu2_arr[i]] -
                     ppblks_new_arr[pnblks1_ilu2_arr[i]] << " Ni_Ext = " <<
                     ppblks_new_arr[pnblks1_ilu2_arr[i]] << endl;
               }

               CFctBxBThreads < _Int, float >::Ilu2HMatrix_BxB (b_blk_wells_temp,
                                                                *_params, blksizeloc,
                                                                pnblks_ilu2_arr[i],
                                                                pnblks1_ilu2_arr[i],
                                                                pnblks2_ilu2_arr[i],
                                                                ppblks_new_arr,
                                                                ppnzord_ilu2_arr,
                                                                ptree_arr[i * 3],
                                                                ptree_arr[i * 3 + 1],
                                                                *ptr_matr,
                                                                pmatrL_float[i],
                                                                pmatrU_float[i],
                                                                ptree_arr[i * 3 + 2],
                                                                nblks_ilu2, blks_ilu2,
                                                                order_ilu2,
                                                                nmodif_scl_loc,
                                                                psclmin_arr[i],
                                                                psclmax_arr[i],
                                                                ndiasplit_loc, nmodif_loc,
                                                                peigmin_arr[i],
                                                                peigmax_arr[i]);

               _stats->nmodif_scl += nmodif_scl_loc;
               _stats->ndiasplit += ndiasplit_loc;
               _stats->nmodif += nmodif_loc;

               nzlu_tot += pmatrL_float[i].GetNzatot ();
               nzlu_tot += pmatrU_float[i].GetNzatot ();

               pnza2cpu[myid] += ptr_matr->GetNzatot ();
               pnzlu2cpu[myid] += 2 * (pmatrL_float[i].GetNzatot ());

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

               if (_params->pfout != NULL && _params->msglev >= 1) {
                  *_params->
                     pfout << " Fct Ihblk: " << i << " Nblks = " << pnblks_ilu2_arr[i] -
                     pnblks1_ilu2_arr[i] << " Nblks_ext = " << pnblks_ilu2_arr[i] <<
                     " Ni_Loc = " << ppblks_new_arr[pnblks_ilu2_arr[i]] -
                     ppblks_new_arr[pnblks1_ilu2_arr[i]] << " Ni_Ext = " <<
                     ppblks_new_arr[pnblks1_ilu2_arr[i]] << endl;
               }
//               sprintf (strbuff_debug,"ChkBIlu_%i.dat",i);
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

               CFctBxBThreads < _Int, double >::Ilu2HMatrix_BxB (b_blk_wells_temp,
                                                                 *_params, blksizeloc,
                                                                 pnblks_ilu2_arr[i],
                                                                 pnblks1_ilu2_arr[i],
                                                                 pnblks2_ilu2_arr[i],
                                                                 ppblks_new_arr,
                                                                 ppnzord_ilu2_arr,
                                                                 ptree_arr[i * 3],
                                                                 ptree_arr[i * 3 + 1],
                                                                 *ptr_matr,
                                                                 pmatrL_double[i],
                                                                 pmatrU_double[i],
                                                                 ptree_arr[i * 3 + 2],
                                                                 nblks_ilu2, blks_ilu2,
                                                                 order_ilu2,
                                                                 nmodif_scl_loc,
                                                                 psclmin_arr[i],
                                                                 psclmax_arr[i],
                                                                 ndiasplit_loc,
                                                                 nmodif_loc,
                                                                 peigmin_arr[i],
                                                                 peigmax_arr[i]);

//               if (false) {
//                  ffout_debug1 << " Tree 3" << endl;
//                  ptree_arr[i * 3+2].OutputTree (ffout_debug1);
//                  ffout_debug1 << " L = " << endl;
//                  pmatrL_double[i].PrintHMatrix (ffout_debug1);
//                  ffout_debug1 << " U = " << endl;
//                  pmatrU_double[i].PrintHMatrix (ffout_debug1);
//               }

               _stats->nmodif_scl += nmodif_scl_loc;
               _stats->ndiasplit += ndiasplit_loc;
               _stats->nmodif += nmodif_loc;

               nzlu_tot += pmatrL_double[i].GetNzatot ();
               nzlu_tot += pmatrU_double[i].GetNzatot ();

               pnza2cpu[myid] += ptr_matr->GetNzatot ();
               pnzlu2cpu[myid] += 2 * (pmatrL_double[i].GetNzatot ());

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
               cout << " Fct End Ihblk =  " << i << " SclMin = " << psclmin_arr[i] <<
                  " SclMax = " << psclmax_arr[i] << " NModif_Scl = " << nmodif_scl_loc <<
                  endl;
               cout << "     NDiaSplit = " << ndiasplit_loc << " PivMin = " <<
                  peigmin_arr[i] << " PivMax = " << peigmax_arr[i] << " NModif = " <<
                  nmodif_loc << endl;
            }

            if (_params->pfout != NULL && _params->msglev >= 1) {
               *_params->
                  pfout << " Fct End Ihblk =  " << i << " SclMin = " << psclmin_arr[i] <<
                  " SclMax = " << psclmax_arr[i] << " NModif_Scl = " << nmodif_scl_loc <<
                  endl;
               *_params->
                  pfout << "     NDiaSplit = " << ndiasplit_loc << " PivMin = " <<
                  peigmin_arr[i] << " PivMax = " << peigmax_arr[i] << " NModif = " <<
                  nmodif_loc << endl;
            }

            pnsplit2cpu[myid] += ndiasplit_loc;
            pnmodif2cpu[myid] += nmodif_loc;

            double time3;

            time3 = CMPIDataExchange::GetWallTimeMPI ();

            pdtime2cpu[myid] += (time3 - time2);

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

         }
      }

//   ffout_debug << " Point 8  " << endl;

      CMPIDataExchange::ExchangeArray (pcommloc, 'I', '+', 1,
                                       (void *) (&_stats->nmodif_scl));
      CMPIDataExchange::ExchangeArray (pcommloc, 'I', '+', 1,
                                       (void *) (&_stats->ndiasplit));
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

         this->slv_float_bxb.InitControl_BxB (pcommloc, blksizeloc, nhblksloc, phblks,
                                              phblk2cpu, phblk2blks, pblk2hblks, pblks);
         this->slv_float_bxb.InitSolveLU_BxB (blksizeloc, phblks_ext, plist_ext_arr,
                                              pnblks_ext_arr, pblksnum_ext_arr,
                                              pblks_ext_arr, ptree_arr, pnblks_ilu2_arr,
                                              pnblks1_ilu2_arr, pnblks2_ilu2_arr,
                                              pblks_ilu2_arr, porder_LU, pmatrL_float,
                                              pmatrU_float);

      } else {

         this->slv_double_bxb.InitControl_BxB (pcommloc, blksizeloc, nhblksloc, phblks,
                                               phblk2cpu, phblk2blks, pblk2hblks, pblks);
         this->slv_double_bxb.InitSolveLU_BxB (blksizeloc, phblks_ext, plist_ext_arr,
                                               pnblks_ext_arr, pblksnum_ext_arr,
                                               pblks_ext_arr, ptree_arr, pnblks_ilu2_arr,
                                               pnblks1_ilu2_arr, pnblks2_ilu2_arr,
                                               pblks_ilu2_arr, porder_LU, pmatrL_double,
                                               pmatrU_double);

      }

// Return back unscaled matrix if necessary

      _params->sctype = sctype_save;

// Compute overall statistics for cpu's

      CMPIDataExchange::ExchangeArray (pcommloc, 'I', '+', nproc, pni2cpu);
      CMPIDataExchange::ExchangeArray (pcommloc, 'I', '+', nproc, pniext2cpu);
      CMPIDataExchange::ExchangeArray (pcommloc, 'I', '+', nproc, pnsplit2cpu);
      CMPIDataExchange::ExchangeArray (pcommloc, 'I', '+', nproc, pnmodif2cpu);
      CMPIDataExchange::ExchangeArray (pcommloc, 'L', '+', nproc, pnza2cpu);
      CMPIDataExchange::ExchangeArray (pcommloc, 'L', '+', nproc, pnzlu2cpu);
      CMPIDataExchange::ExchangeArray (pcommloc, 'D', '+', nproc, pdtime2cpu);

      for (i = 0; i < nproc; i++) {
         if (pni2cpu[i] > 0) {
            pdns2cpu[i] = (double) (pnzlu2cpu[i]) / (double) (pnza2cpu[i]);
         }
      }

      if (_params->msglev >= 4 && myid == 0) {
         PrintArray (cout, " Ni2Cpu ", nproc, pni2cpu);
         PrintArray (cout, " NiExt2Cpu ", nproc, pniext2cpu);
         PrintArray (cout, " NSplit2Cpu ", nproc, pnsplit2cpu);
         PrintArray (cout, " NModif2Cpu ", nproc, pnmodif2cpu);
         PrintArrayLow (cout, " Dnsty2Cpu ", nproc, pdns2cpu);
         PrintArrayLow (cout, " DTime2Cpu ", nproc, pdtime2cpu);
      }

      if (_params->msglev >= 3 && myid == 0 && _params->pfout != NULL) {
         PrintArray (*_params->pfout, " Ni2Cpu ", nproc, pni2cpu);
         PrintArray (*_params->pfout, " NiExt2Cpu ", nproc, pniext2cpu);
         PrintArray (*_params->pfout, " NSplit2Cpu ", nproc, pnsplit2cpu);
         PrintArray (*_params->pfout, " NModif2Cpu ", nproc, pnmodif2cpu);
         PrintArrayLow (*_params->pfout, " Dnsty2Cpu ", nproc, pdns2cpu);
         PrintArrayLow (*_params->pfout, " DTime2Cpu ", nproc, pdtime2cpu);
      }
// Finalize timer

      CMPIDataExchange::Synchronize (pcommloc);

      double time1;

      time1 = CMPIDataExchange::GetWallTimeMPI ();

      _stats->dtime_fct = time1 - time0;

//   ffout_debug << " Return point !!! " << endl;

   }

//
// Perform ILU2 point factorization of the hmatrix with dinamic ordering and future diagonal modifications
//========================================================================================
   template < typename _Int, typename _Flt > void CFctBxBThreads < _Int,
      _Flt >::Ilu2HMatrix_BxB (bool _b_blk_wells, SParams & _params, int _blksize,
                               int _nblks, int _nblks1, int _nblks2, long long *_blks,
                               long long *_nzord_blks, CTree & _tree0, CTree & _tree1,
                               CBMatrix < _Int, _Flt > &_a_matr, CBMatrix < _Int,
                               _Flt > &_l_matr, CBMatrix < _Int, _Flt > &_u_matr,
                               CTree & _tree2_new, int &_nblks_new,
                               vector < long long >&_blks_new,
                               CVectorData < int >&_ordernew, int &_nmodif_scl,
                               double &_sclmin_att, double &_sclmax_att, int &_ndiasplit,
                               int &_nmodif, double &_eigmin_att, double &_eigmax_att)
   {

      int b_2 = _blksize * _blksize;

      int nsuptot = (int) _blks[_nblks];

      int iparam1_out = _params.iparam1;

// Get params

      int sctype = _params.sctype;
      double sclmin = _params.sclmin;
      int nitersc = _params.nitersc;

// Compute explicit scaling

//   ffout << " Ilu2HMatrix: Point 1" << endl;

      CVectorData < _Flt > scl_L;
      CVectorData < _Flt > scl_U;
      CVectorData < _Flt > scl_L_inv;
      CVectorData < _Flt > scl_U_inv;

      CFctBxBThreads < _Int, _Flt >::ComputeExplicitScaling_BxB (sctype, nitersc, sclmin,
                                                                 _blksize, _nblks, _blks,
                                                                 _a_matr, scl_L, scl_U,
                                                                 scl_L_inv, scl_U_inv,
                                                                 _nmodif_scl, _sclmin_att,
                                                                 _sclmax_att);

      _Flt *pscl_L = scl_L.Ptr ();
      _Flt *pscl_U = scl_U.Ptr ();
      _Flt *pscl_L_inv = scl_L_inv.Ptr ();
      _Flt *pscl_U_inv = scl_U_inv.Ptr ();

// Create a_matr copy

      CBMatrix < _Int, _Flt > amatr_copy;

      amatr_copy.AddReplaceBxBThr (true, _blksize, _a_matr);

// Perform explicit scaling

      CFctBxBThreads < _Int, _Flt >::HMatrixScale_BxB (_blksize, _nblks, _blks,
                                                       amatr_copy, pscl_L, pscl_U);

// Split and combine by pairs

      CBMatrix < _Int, _Flt > alu_pair;

      CFctBxBThreads < _Int, _Flt >::SplitLUPair_BxB (false, _blksize, _nblks, _blks,
                                                      amatr_copy, alu_pair);

//      amatr_copy.Clean();

// Perform fct

      CBMatrix < _Int, _Flt > lu_pair;

      CFctBxBThreads < _Int, _Flt >::Ilu2BlockIlu2_BxB (_b_blk_wells, NULL, iparam1_out,
                                                        &_params, _blksize, _nblks,
                                                        _nblks1, _nblks2, _blks,
                                                        _nzord_blks, _tree0, _tree1,
                                                        alu_pair, _tree2_new, _nblks_new,
                                                        _blks_new, _ordernew, lu_pair,
                                                        _ndiasplit, _nmodif, _eigmin_att,
                                                        _eigmax_att);

      long long *p_blks_new = &_blks_new[0];
      int *p_ordernew = _ordernew.Ptr ();

// Split pair into L and U parts

      CFctBxBThreads < _Int, _Flt >::SplitPair_BxB (false, _blksize, _nblks_new,
                                                    p_blks_new, lu_pair, _l_matr,
                                                    _u_matr);
/*
      if (false) {

         CBMatrix < _Int, _Flt > a_hmatr_ord;

         int nblks_1blk = 1;

         long long blks_1blk[2];

         blks_1blk[0] = 0;
         blks_1blk[1] = nsuptot;

         vector<int> order_new_int (nsuptot);
         int *porder_new_int = order_new_int.data();

         for (int i=0;i<nsuptot;i++) porder_new_int[i] = (int)p_ordernew[i];

         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (_blksize, _nblks, _blks,
                                                        amatr_copy, p_ordernew,
                                                        nblks_1blk, blks_1blk,
                                                        a_hmatr_ord);

         CMatrix < _Int, _Flt > *pa_matr_ord = a_hmatr_ord.GetASubArr();

         vector<int> order_id (nsuptot);
         int *porder_id = order_id.data();

         CVectorInt<int>::SetIdentity (nsuptot, 0, porder_id);

         CBMatrix < _Int, _Flt > l_hmatr_1blk;
         CBMatrix < _Int, _Flt > u_hmatr_1blk;

         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (_blksize, _nblks_new, p_blks_new,
                                                        _l_matr, porder_id,
                                                        nblks_1blk, blks_1blk,
                                                        l_hmatr_1blk);
         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (_blksize, _nblks_new, p_blks_new,
                                                        _u_matr, porder_id,
                                                        nblks_1blk, blks_1blk,
                                                        u_hmatr_1blk);

         CMatrix < _Int, _Flt > *pl_matr_1blk = l_hmatr_1blk.GetASubArr();
         CMatrix < _Int, _Flt > *pu_matr_1blk = u_hmatr_1blk.GetASubArr();

         int ntot = nsuptot*_blksize;
         vector<double> x (ntot);
         vector<double> ax (ntot);
         vector<double> y (ntot);
         vector<double> z (ntot);
         vector<double> work (ntot);
         double *px = x.data();
         double *pax = ax.data();
         double *py = y.data();
         double *pz = z.data();
         double *pwork = work.data();

         int i;

         for (i=0;i<ntot;i++) px[i] = (double) (i%20) - 10.0e0+0.5e0;
         for (i=0;i<ntot;i++) pax[i] = 0.0e0;

         CMvmSlv_BxB<_Int, _Flt, double>::MvmA ('+', _blksize, *pa_matr_ord,
                                          px, pax);
         CMvmSlv_BxB<_Int, _Flt, double>::SolveL (_blksize, *pl_matr_1blk,
                                          pax, py, pwork);
         CMvmSlv_BxB<_Int, _Flt, double>::SolveU (_blksize, *pu_matr_1blk,
                                          py, pz, pwork);

         double diff = 0.0e0;
         double aux;

         for (i=0;i<ntot;i++) {
            aux = pz[i] - px[i];
            diff += aux*aux;
         }

         diff = sqrt(diff);

         cout << " Ntot = " << ntot << " Nblks_new = " << _nblks_new << " Scaled Diff = " << diff << endl;

// Split

         CMatrix<_Int,_Flt> al_matr;
         CMatrix<_Int,_Flt> au_matr;

         CFct_bxb<_Int,_Flt>::SplitLU (false, _blksize, nsuptot,
                                       *pa_matr_ord,
                                       al_matr, au_matr);

// Transpose

         CMatrix<_Int,_Flt> alt_matr;

         CFct_bxb<_Int,_Flt>::Transpose (false, _blksize, nsuptot,
                                          al_matr,
                                          alt_matr);

// Combine

         CMatrix<_Int,_Flt> alu_pair_matr;

         CFct_bxb<_Int,_Flt>::CombineLUPairs (false, _blksize, nsuptot,
                                                alt_matr, au_matr,
                                                alu_pair_matr);

         int nzja_alu = alu_pair_matr.GetNzja();

         alu_pair_matr.ResizeJaChar (nzja_alu);
         alu_pair_matr.SetNzjaChar (nzja_alu);

         {
            char *pjachar = alu_pair_matr.GetJaCharArr();
            int i;
            for (i=0;i<nzja_alu;i++) pjachar[i] = 0;
         }

// Compute block fct

         CMatrix<_Int,_Flt> lu_pair_matr;

         CMatrix<_Int,_Flt> alu_schur_matr;
         CMatrix<_Int,_Flt> lu_schur_matr;

         CMatrix<_Int,_Flt> lu_schur_sum;
         CMatrix<_Int,_Flt> lu_schur_fct;

         int nmodif;
         double eigmin_att, eigmax_att;

         if (_nblks_new > 1) {

            int nsup_ini = (int)p_blks_new[_nblks_new-1];

            int iparam1 = _params.iparam1;
            int iparam2 = _params.iparam2;
            int iparam3 = _params.iparam3;

            _params.b_write_file = true;

            sprintf (_params.name_file,"ChkSchur_%i_%i.dat",nsuptot,nsup_ini);

            _params.iparam1 = nsup_ini;
            _params.iparam2 = -1;
            _params.iparam3 = -1;

            CFct_bxb<_Int,_Flt>::FctBlockIlu2Degree (&_params, _blksize, nsuptot, nsup_ini,
                                                     alu_pair_matr,
                                                     lu_pair_matr,
                                                     nmodif, eigmin_att, eigmax_att);

            _params.b_write_file = false;
            _params.iparam1 = iparam1;
            _params.iparam2 = iparam2;
            _params.iparam3 = iparam3;

            int size_bxb = _blksize*_blksize*2;

            alu_pair_matr.GetLastSubmatrixBxB (true, size_bxb, nsup_ini, alu_schur_matr);
            lu_pair_matr.GetLastSubmatrixBxB (true, size_bxb, nsup_ini, lu_schur_matr);

            lu_schur_sum.AddBlocksPairsBxB ('+', _blksize, alu_schur_matr, lu_schur_matr);

            CFct_bxb<_Int,_Flt>::FctBlockIlu2Degree (&_params, _blksize, nsuptot-nsup_ini, nsuptot-nsup_ini,
                                                     lu_schur_sum,
                                                     lu_schur_fct,
                                                     nmodif, eigmin_att, eigmax_att);

         }

         CFct_bxb<_Int,_Flt>::FctBlockIlu2Degree (&_params, _blksize, nsuptot, nsuptot,
                                                  alu_pair_matr,
                                                  lu_pair_matr,
                                                  nmodif, eigmin_att, eigmax_att);

// Split result

         CMatrix<_Int,_Flt> l_matr;
         CMatrix<_Int,_Flt> u_matr;

         CFct_bxb<_Int,_Flt>::SplitLUPairs (false, _blksize, nsuptot,
                                             lu_pair_matr,
                                             l_matr, u_matr);

// Create 1 block hblocks

         CBMatrix<_Int,_Flt> l_hmatr_temp;
         CBMatrix<_Int,_Flt> u_hmatr_temp;

         l_hmatr_temp.ResizeASub (1);
         l_hmatr_temp.SetNzblk (1);
         u_hmatr_temp.ResizeASub (1);
         u_hmatr_temp.SetNzblk (1);

         CMatrix<int,float> *phblk_L = l_hmatr_temp.GetHMatrStr();

         phblk_L->ResizeAndSetAllSp (1, 1, 1, 1);

         int *plist_L  = phblk_L->GetListArr();
         int *plist2_L = phblk_L->GetList2Arr();
         int *pia_L    = phblk_L->GetIaArr();
         int *pja_L    = phblk_L->GetJaArr();
         int *pja2_L   = phblk_L->GetJa2Arr();

         plist_L[0]  = 0;
         plist2_L[0] = 0;
         pia_L[0]    = 0;
         pia_L[1]    = 1;
         pja_L[0]    = 0;
         pja2_L[0]   = 0;

         CMatrix<_Int,_Flt> *pasub_L = l_hmatr_temp.GetASubArr();

         pasub_L[0] = l_matr;

         CMatrix<int,float> *phblk_U = u_hmatr_temp.GetHMatrStr();

         phblk_U->ResizeAndSetAllSp (1, 1, 1, 1);

         int *plist_U  = phblk_U->GetListArr();
         int *plist2_U = phblk_U->GetList2Arr();
         int *pia_U    = phblk_U->GetIaArr();
         int *pja_U    = phblk_U->GetJaArr();
         int *pja2_U   = phblk_U->GetJa2Arr();

         plist_U[0]  = 0;
         plist2_U[0] = 0;
         pia_U[0]    = 0;
         pia_U[1]    = 1;
         pja_U[0]    = 0;
         pja2_U[0]   = 0;

         CMatrix<_Int,_Flt> *pasub_U = u_hmatr_temp.GetASubArr();

         pasub_U[0] = u_matr;

// Convert to hmatrix

         CBMatrix<_Int,_Flt> l_hmatr;
         CBMatrix<_Int,_Flt> u_hmatr;

         CBMatrix<_Int,_Flt>::ReorderHMatrix_BxB (_blksize, nblks_1blk, blks_1blk,
                                                  l_hmatr_temp, porder_id,
                                                  _nblks_new, p_blks_new,
                                                  l_hmatr);

         CBMatrix<_Int,_Flt>::ReorderHMatrix_BxB (_blksize, nblks_1blk, blks_1blk,
                                                  u_hmatr_temp, porder_id,
                                                  _nblks_new, p_blks_new,
                                                  u_hmatr);

// Print results

         char strbuff[256];

         sprintf (strbuff,"Fct_%i.dat",nsuptot);

         ofstream ffout (strbuff);

         ffout << " nblks1 = " << _nblks1 << " nblks2 = " << _nblks2 << endl;
         PrintArray (ffout," Blks ini ",_nblks+1,_blks);
         PrintArray (ffout," Blks ",_nblks_new+1,p_blks_new);

         ffout << " Ini matr for last block " << endl;
         alu_schur_matr.PrintMatrix (ffout);
         ffout << " Schur for last block " << endl;
         lu_schur_matr.PrintMatrix (ffout);

         ffout << " Schur sum for last block " << endl;
         lu_schur_sum.PrintMatrix (ffout);

         ffout << " Schur fct for last block " << endl;
         lu_schur_fct.PrintMatrix (ffout);

         sprintf (strbuff,"Fct_1blk_%i.dat",nsuptot);

         ofstream ffout1 (strbuff);

         ffout1 << " nblks1 = " << _nblks1 << " nblks2 = " << _nblks2 << endl;
         PrintArray (ffout1," Blks ini ",_nblks+1,_blks);
         PrintArray (ffout1," Blks ",_nblks_new+1,p_blks_new);

         ffout1 << " L comp" << endl;
         l_matr.PrintMatrixRows (ffout1, _blksize);
         ffout1 << " U comp" << endl;
         u_matr.PrintMatrixRows (ffout1, _blksize);

      }
*/
// Reorder scaling and compute its inverse

      CVectorData < _Flt > temp_LU (nsuptot * b_2);
      _Flt *ptemp_LU = temp_LU.Ptr ();

      CVector < _Flt >::OrderVector_thr (b_2, nsuptot, p_ordernew, pscl_L, ptemp_LU);
      CVector < _Flt >::CopyVector_thr (nsuptot * b_2, ptemp_LU, pscl_L);

      CVector < _Flt >::OrderVector_thr (b_2, nsuptot, p_ordernew, pscl_U, ptemp_LU);
      CVector < _Flt >::CopyVector_thr (nsuptot * b_2, ptemp_LU, pscl_U);

      CVector < _Flt >::OrderVector_thr (b_2, nsuptot, p_ordernew, pscl_L_inv, ptemp_LU);
      CVector < _Flt >::CopyVector_thr (nsuptot * b_2, ptemp_LU, pscl_L_inv);

      CVector < _Flt >::OrderVector_thr (b_2, nsuptot, p_ordernew, pscl_U_inv, ptemp_LU);
      CVector < _Flt >::CopyVector_thr (nsuptot * b_2, ptemp_LU, pscl_U_inv);

// Rescale factors

      CFctBxBThreads < _Int, _Flt >::RescaleU_BxB (_blksize, _nblks_new, p_blks_new,
                                                   _l_matr, pscl_L, pscl_L_inv);
      CFctBxBThreads < _Int, _Flt >::RescaleU_BxB (_blksize, _nblks_new, p_blks_new,
                                                   _u_matr, pscl_U, pscl_U_inv);

// Check fct
/*
      if (false) {

         CBMatrix < _Int, _Flt > a_hmatr_ord;

         int nblks_1blk = 1;

         long long blks_1blk[2];

         blks_1blk[0] = 0;
         blks_1blk[1] = nsuptot;

         vector<int> order_new_int (nsuptot);
         int *porder_new_int = order_new_int.data();

         for (int i=0;i<nsuptot;i++) porder_new_int[i] = (int)p_ordernew[i];

         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (_blksize, _nblks, _blks,
                                                        _a_matr, p_ordernew,
                                                        nblks_1blk, blks_1blk,
                                                        a_hmatr_ord);

         CMatrix < _Int, _Flt > *pa_matr_ord = a_hmatr_ord.GetASubArr();

         vector<int> order_id (nsuptot);
         int *porder_id = order_id.data();

         CVectorInt<int>::SetIdentity (nsuptot, 0, porder_id);

         CBMatrix < _Int, _Flt > l_hmatr_1blk;
         CBMatrix < _Int, _Flt > u_hmatr_1blk;

         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (_blksize, _nblks_new, p_blks_new,
                                                        _l_matr, porder_id,
                                                        nblks_1blk, blks_1blk,
                                                        l_hmatr_1blk);
         CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (_blksize, _nblks_new, p_blks_new,
                                                        _u_matr, porder_id,
                                                        nblks_1blk, blks_1blk,
                                                        u_hmatr_1blk);

         CMatrix < _Int, _Flt > *pl_matr_1blk = l_hmatr_1blk.GetASubArr();
         CMatrix < _Int, _Flt > *pu_matr_1blk = u_hmatr_1blk.GetASubArr();

         int ntot = nsuptot*_blksize;
         vector<double> x (ntot);
         vector<double> ax (ntot);
         vector<double> y (ntot);
         vector<double> z (ntot);
         vector<double> work (ntot);
         double *px = x.data();
         double *pax = ax.data();
         double *py = y.data();
         double *pz = z.data();
         double *pwork = work.data();

         int i;

         for (i=0;i<ntot;i++) px[i] = (double) (i%20) - 10.0e0+0.5e0;
         for (i=0;i<ntot;i++) pax[i] = 0.0e0;

         CMvmSlv_BxB<_Int, _Flt, double>::MvmA ('+', _blksize, *pa_matr_ord,
                                          px, pax);
         CMvmSlv_BxB<_Int, _Flt, double>::SolveL (_blksize, *pl_matr_1blk,
                                          pax, py, pwork);
         CMvmSlv_BxB<_Int, _Flt, double>::SolveU (_blksize, *pu_matr_1blk,
                                          py, pz, pwork);

         double diff = 0.0e0;
         double aux;

         for (i=0;i<ntot;i++) {
            aux = pz[i] - px[i];
            diff += aux*aux;
         }

         diff = sqrt(diff);

         cout << " Ntot = " << ntot << " Nblks_new = " << _nblks_new << " Unscaled Diff = " << diff << endl;

      }
*/
   }

//
// Compute explicit scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CFctBxBThreads < _Int,
      _Flt >::ComputeExplicitScaling_BxB (int _sctype, int _nitersc, double _sclmin,
                                          int _blksize, int _nblks, long long *_blks,
                                          CBMatrix < _Int, _Flt > &_a_matr,
                                          CVectorData < _Flt > &_scl_L,
                                          CVectorData < _Flt > &_scl_U,
                                          CVectorData < _Flt > &_scl_L_inv,
                                          CVectorData < _Flt > &_scl_U_inv, int &_n_modif,
                                          double &_sclmin_att, double &_sclmax_att)
   {

//      ofstream ffout ("ChkScl.dat");

//      _a_matr.PrintHMatrix (ffout);

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

// Allocate scaling

      int b_2 = _blksize * _blksize;

      int nsuptot = (int) _blks[_nblks];

      _scl_L.resize (nsuptot * b_2);
      _scl_U.resize (nsuptot * b_2);
      _scl_L_inv.resize (nsuptot * b_2);
      _scl_U_inv.resize (nsuptot * b_2);

      _Flt *p_sclL = _scl_L.Ptr ();
      _Flt *p_sclU = _scl_U.Ptr ();
      _Flt *p_sclL_inv = _scl_L_inv.Ptr ();
      _Flt *p_sclU_inv = _scl_U_inv.Ptr ();

// Open hmatrix

      CMatrix < int, float >*phmatr = _a_matr.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub = _a_matr.GetASubArr ();

      int nzja_hmatr = phmatr->GetNzja ();
      int *pia_hmatr = phmatr->GetIaArr ();
      int *pja_hmatr = phmatr->GetJaArr ();

// Simple diagonal based scaling

      _Flt fone;

      CVector < _Flt >::SetByOnes (1, &fone);

      _n_modif = 0;
      _sclmin_att = 1.0e100;
      _sclmax_att = -1.0e100;

      if (_sctype == -1) {

         CVector < _Flt >::SetByZeroes_thr (nsuptot * b_2, p_sclL);
         CVector < _Flt >::SetByZeroes_thr (nsuptot * b_2, p_sclU);
         CVector < _Flt >::SetByZeroes_thr (nsuptot * b_2, p_sclL_inv);
         CVector < _Flt >::SetByZeroes_thr (nsuptot * b_2, p_sclU_inv);

         int ni_part = nsuptot / n_thr;

         {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < n_thr; ipar++) {

               int ibeg, iend, ni_loc;

               ibeg = ipar * ni_part;
               iend = (ipar + 1) * ni_part - 1;
               if (ipar == n_thr - 1)
                  iend = nsuptot - 1;

               ni_loc = iend + 1 - ibeg;

               int k, kk;

               for (k = ibeg; k <= iend; k++) {
                  for (kk = 0; kk < _blksize; kk++) {
                     p_sclL[k * b_2 + kk * _blksize + kk] = fone;
                     p_sclU[k * b_2 + kk * _blksize + kk] = fone;
                     p_sclL_inv[k * b_2 + kk * _blksize + kk] = fone;
                     p_sclU_inv[k * b_2 + kk * _blksize + kk] = fone;
                  }
               }

            }

         }

         _Flt diag;

         diag = (_Flt) 1.;

         _sclmin_att = diag;
         _sclmax_att = diag;

      } else if (_sctype == 0) {

         vector < int >nmodif_thr (n_thr + 1);
         vector < double >sclmin_thr (n_thr + 1);
         vector < double >sclmax_thr (n_thr + 1);

         CVectorData < _Flt > Work_thr (20 * b_2 * n_thr);
         CVectorData < double >dWork_thr (20 * b_2 * n_thr);

         _Flt *pWork_thr = Work_thr.Ptr ();
         double *pdWork_thr = dWork_thr.Ptr ();

         int *pnmodif_thr = &nmodif_thr[0];
         double *psclmin_thr = &sclmin_thr[0];
         double *psclmax_thr = &sclmax_thr[0];

         int jjj;

         for (jjj = 0; jjj < n_thr; jjj++)
            pnmodif_thr[jjj] = 0;
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
            int iend = (int) _blks[ipar + 1] - 1;

            CVector < _Flt >::SetByZeroes_thr (niloc * b_2, p_sclL + ibeg * b_2);
            CVector < _Flt >::SetByZeroes_thr (niloc * b_2, p_sclU + ibeg * b_2);

            int k, kk;

            for (k = ibeg; k <= iend; k++) {
               for (kk = 0; kk < _blksize; kk++) {
                  p_sclL[k * b_2 + kk * _blksize + kk] = fone;
                  p_sclU[k * b_2 + kk * _blksize + kk] = fone;
               }
            }

            int j, jj, ki, irow, kj, kjj, jbs, ibs;
            double sclmin_temp, sclmax_temp;

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
                     ibs = (ibeg + irow) * b_2;
                     for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                        kjj = (int) pja_temp[kj];
                        if (kjj == irow) {
                           jbs = kj * b_2;
//                           ffout << "   Irow = " << irow << " kjj = " << kjj << endl;
//                           PrintArray (ffout," Diag ",b_2,pa_temp+jbs);
                           if (_sctype == 0) {
                              pnmodif_thr[my_thr] +=
                                 CBlock_BxB_traits < _Flt >::FctInv_BxB (_sclmin,
                                                                         _blksize,
                                                                         pa_temp + jbs,
                                                                         p_sclL_inv + ibs,
                                                                         p_sclU_inv + ibs,
                                                                         p_sclL + ibs,
                                                                         p_sclU + ibs,
                                                                         sclmin_temp,
                                                                         sclmax_temp);
                           } else {
                              pnmodif_thr[my_thr] +=
                                 CBlock_BxB_traits < _Flt >::FctInv_Svd_BxB (_sclmin,
                                                                             _blksize,
                                                                             pa_temp +
                                                                             jbs,
                                                                             p_sclL_inv +
                                                                             ibs,
                                                                             p_sclU_inv +
                                                                             ibs,
                                                                             p_sclL + ibs,
                                                                             p_sclU + ibs,
                                                                             pWork_thr +
                                                                             my_thr * 20 *
                                                                             b_2,
                                                                             pdWork_thr +
                                                                             my_thr * 20 *
                                                                             b_2,
                                                                             sclmin_temp,
                                                                             sclmax_temp);
                           }
                           if (sclmin_temp < psclmin_thr[my_thr])
                              psclmin_thr[my_thr] = sclmin_temp;
                           if (sclmax_temp > psclmax_thr[my_thr])
                              psclmax_thr[my_thr] = sclmax_temp;
                        }
                     }
                  }
               }
            }
         }

         for (jjj = 0; jjj < n_thr; jjj++) {
            _n_modif += pnmodif_thr[jjj];
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

// For all blocks compute its columns/rows front sizes

         CVectorData < CVectorData < int > >ia_cols (_nblks);
         CVectorData < CVectorData < int > >ia_rows (_nblks);

         CVectorData < int >*pia_cols = ia_cols.Ptr ();
         CVectorData < int >*pia_rows = ia_rows.Ptr ();

         {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < _nblks; ipar++) {

               int ni = (int) (_blks[ipar + 1] - _blks[ipar]);

               pia_cols[ipar].resize (ni + 1);
               pia_rows[ipar].resize (ni + 1);

               int *ppia_cols = pia_cols[ipar].Ptr ();
               int *ppia_rows = pia_rows[ipar].Ptr ();

               int j;

               for (j = 0; j <= ni; j++)
                  ppia_cols[j] = 0;
               for (j = 0; j <= ni; j++)
                  ppia_rows[j] = 0;

               int k, irow;

               for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
                  int nlist_temp = pA_sub[j].GetNlist ();
                  _Int *plist_temp = pA_sub[j].GetListArr ();
                  _Int *pia_temp = pA_sub[j].GetIaArr ();
                  for (k = 0; k < nlist_temp; k++) {
                     irow = (int) plist_temp[k];
                     ppia_rows[irow + 1] += (int) (pia_temp[k + 1] - pia_temp[k]);
                  }
               }

               int jj, ind, kj, kk;

               for (j = pia_hmatrT[ipar]; j < pia_hmatrT[ipar + 1]; j++) {
                  jj = pja_hmatrT[j];
                  ind = pindHt2H[j];
                  int nlist_temp = pA_sub[ind].GetNlist ();
                  _Int *pia_temp = pA_sub[ind].GetIaArr ();
                  _Int *pja_temp = pA_sub[ind].GetJaArr ();
                  for (k = 0; k < nlist_temp; k++) {
                     for (kj = (int) pia_temp[k]; kj < pia_temp[k + 1]; kj++) {
                        kk = (int) pja_temp[kj];
                        ppia_cols[kk + 1]++;
                     }
                  }
               }

               for (j = 0; j < ni; j++)
                  ppia_cols[j + 1] += ppia_cols[j];
               for (j = 0; j < ni; j++)
                  ppia_rows[j + 1] += ppia_rows[j];

//               ffout << " Iblk = " << ipar << endl;
//               PrintArray (ffout," ia_cols = ",ni+1,ppia_cols);
//               PrintArray (ffout," ia_rows = ",ni+1,ppia_rows);

            }
         }

// Compute maximal used profile memory for block row/column

         int niblk_max = 0;
         int nzja_max = 0;
         int ni_max = 0;

         {
            int i, j, njloc;
            for (i = 0; i < _nblks; i++) {
               int niloc = (int) (_blks[i + 1] - _blks[i]);
               if (niloc > niblk_max)
                  niblk_max = niloc;
               int *ppia_rows = pia_rows[i].Ptr ();
               for (j = 0; j < niloc; j++) {
                  njloc = ppia_rows[j + 1] - ppia_rows[j];
                  if (njloc > ni_max)
                     ni_max = njloc;
               }
               int *ppia_cols = pia_cols[i].Ptr ();
               for (j = 0; j < niloc; j++) {
                  njloc = ppia_cols[j + 1] - ppia_cols[j];
                  if (njloc > ni_max)
                     ni_max = njloc;
               }
               int nzja_loc = ppia_rows[niloc];
               if (nzja_max < nzja_loc)
                  nzja_max = nzja_loc;
               nzja_loc = ppia_cols[niloc];
               if (nzja_max < nzja_loc)
                  nzja_max = nzja_loc;
            }
         }

// Iterative rows/columns balancing block scaling via QR decomposition

         double sclRmin, sclRmax;
         double sclCmin, sclCmax;

         sclRmin = 1.0e100;
         sclRmax = -1.0e100;
         sclCmin = 1.0e100;
         sclCmax = -1.0e100;

         vector < double >sclRmin_thr (n_thr + 1);
         vector < double >sclRmax_thr (n_thr + 1);
         vector < double >sclCmin_thr (n_thr + 1);
         vector < double >sclCmax_thr (n_thr + 1);

         double *psclRmin_thr = &sclRmin_thr[0];
         double *psclRmax_thr = &sclRmax_thr[0];
         double *psclCmin_thr = &sclCmin_thr[0];
         double *psclCmax_thr = &sclCmax_thr[0];

         vector < int >nmodif_thr (n_thr + 1);
         int *pnmodif_thr = &nmodif_thr[0];

         {

            int jjj;

            for (jjj = 0; jjj < n_thr; jjj++)
               psclRmin_thr[jjj] = 1.0e100;
            for (jjj = 0; jjj < n_thr; jjj++)
               psclRmax_thr[jjj] = -1.0e100;
            for (jjj = 0; jjj < n_thr; jjj++)
               psclCmin_thr[jjj] = 1.0e100;
            for (jjj = 0; jjj < n_thr; jjj++)
               psclCmax_thr[jjj] = -1.0e100;
            for (jjj = 0; jjj < n_thr; jjj++)
               pnmodif_thr[jjj] = 0;

         }

         {

// Compute preliminary rows block scaling based on rows norms

            int ntot = nsuptot * _blksize;

            CVectorData < _Flt > fnrm2 (ntot);
            _Flt *pfnrm2 = fnrm2.Ptr ();

            CVector < _Flt >::SetByZeroes_thr (ntot, pfnrm2);

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < _nblks; ipar++) {

               int j, ki, kj, jcol, irow;
               int kii, kjj;
               _Flt *pelems;
               _Flt aux;

               int ibeg = (int) _blks[ipar];

               for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
//                  int jblk = pja_hmatr[j];
//                  int jbeg = (int) _blks[jblk];
                  int nlist_temp = pA_sub[j].GetNlist ();
                  _Int *plist_temp = pA_sub[j].GetListArr ();
                  _Int *pia_temp = pA_sub[j].GetIaArr ();
                  _Int *pja_temp = pA_sub[j].GetJaArr ();
                  _Flt *pa_temp = pA_sub[j].GetAArr ();
                  for (ki = 0; ki < nlist_temp; ki++) {
                     irow = (int) plist_temp[ki];
                     for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                        jcol = (int) pja_temp[kj];
                        pelems = pa_temp + kj * b_2;
                        for (kjj = 0; kjj < _blksize; kjj++) {
                           for (kii = 0; kii < _blksize; kii++) {
                              aux = pelems[kjj * _blksize + kii];
                              pfnrm2[(ibeg + irow) * _blksize + kii] += aux * aux;
                           }
                        }
                     }
                  }
               }

            }

//            PrintArray (ffout," pfnrm2 = ",ntot,pfnrm2);

            CVector < _Flt >::SetByZeroes_thr (nsuptot * b_2, p_sclL);
            CVector < _Flt >::SetByZeroes_thr (nsuptot * b_2, p_sclL_inv);

#ifdef USE_THREADS
#pragma omp parallel for
#endif
            for (int ipar = 0; ipar < _nblks; ipar++) {

//               int ibeg = (int) _blks[ipar];
//               int niloc = (int) (_blks[ipar+1]-_blks[ipar]);

               int j, kjj, k;
               _Flt aux;

               for (j = (int) _blks[ipar]; j < _blks[ipar + 1]; j++) {
                  for (kjj = 0; kjj < _blksize; kjj++) {
                     k = j * _blksize + kjj;
                     aux = pfnrm2[k];
                     aux = sqrt (aux);
                     if (aux < _sclmin)
                        aux = (_Flt) _sclmin;
                     aux = sqrt (aux);
                     p_sclL_inv[j * b_2 + kjj * _blksize + kjj] = aux;
                     aux = fone / aux;
                     p_sclL[j * b_2 + kjj * _blksize + kjj] = aux;
                  }
               }
            }

//            PrintArray (ffout," p_sclL ini = ",nsuptot*b_2,p_sclL);
//            PrintArray (ffout," p_sclL inv ini = ",nsuptot*b_2,p_sclL_inv);

         }

         {

// Work memory for each thread

            CVectorData < CVectorData < int > >iptr_thr (n_thr);
            CVectorData < CVectorData < _Flt > >a_front_thr (n_thr);
            CVectorData < CVectorData < _Flt > >work_front_thr (n_thr);

            CVectorData < int >*piptr_thr = iptr_thr.Ptr ();
            CVectorData < _Flt > *pa_front_thr = a_front_thr.Ptr ();
            CVectorData < _Flt > *pwork_front_thr = work_front_thr.Ptr ();

// Main scaling iterative cycle

            int iter;

            for (iter = 0; iter < _nitersc; iter++) {

//               ffout << " Iter scl = " << iter << endl;

// Columns block scaling

               {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < _nblks; ipar++) {

                     int niloc = (int) (_blks[ipar + 1] - _blks[ipar]);

                     int ibeg = (int) _blks[ipar];

                     int my_thr = 0;
#ifdef USE_THREADS
                     my_thr = omp_get_thread_num ();
#endif

// Memory support

                     if (piptr_thr[my_thr].GetLength () < niblk_max) {
                        piptr_thr[my_thr].resize (niblk_max);
                     }

                     int *ppiptr_th = piptr_thr[my_thr].Ptr ();


                     if (pa_front_thr[my_thr].GetLength () < nzja_max * b_2) {
                        pa_front_thr[my_thr].resize (nzja_max * b_2);
                     }

                     _Flt *ppa_front_th = pa_front_thr[my_thr].Ptr ();

                     int nz_work = b_2 * (ni_max + 1) * 2;

                     if (pwork_front_thr[my_thr].GetLength () < nz_work) {
                        pwork_front_thr[my_thr].resize (nz_work);
                     }

                     _Flt *ppwork_front_th = pwork_front_thr[my_thr].Ptr ();

                     _Flt *py_blk = ppwork_front_th;
                     _Flt *pz_blk = py_blk + ni_max * b_2;
                     _Flt *pr_matr = pz_blk + ni_max * b_2;
                     _Flt *ptau_matr = pr_matr + b_2;

                     int *ppia_cols = pia_cols[ipar].Ptr ();

                     int j, jblk, ind, ki, kj, kk, ind1, irow;

                     for (j = 0; j < niloc; j++)
                        ppiptr_th[j] = ppia_cols[j];

                     for (j = pia_hmatrT[ipar]; j < pia_hmatrT[ipar + 1]; j++) {
                        jblk = pja_hmatrT[j];
                        ind = pindHt2H[j];
                        int jbeg = (int) _blks[jblk];
                        int nlist_temp = pA_sub[ind].GetNlist ();
                        _Int *plist_temp = pA_sub[ind].GetListArr ();
                        _Int *pia_temp = pA_sub[ind].GetIaArr ();
                        _Int *pja_temp = pA_sub[ind].GetJaArr ();
                        _Flt *pa_temp = pA_sub[ind].GetAArr ();

                        for (ki = 0; ki < nlist_temp; ki++) {
                           irow = (int) plist_temp[ki];
                           for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                              kk = (int) pja_temp[kj];
                              ind1 = ppiptr_th[kk];
                              if (jbeg + irow >= nsuptot) {
                                 cout << " Error: too large index !!! " << endl;
                                 throw " Error: too large index !!! ";
                              }
                              CBlock_BxB_traits < _Flt >::MtM_BxB (_blksize,
                                                                   pa_temp + kj * b_2,
                                                                   p_sclL + (jbeg +
                                                                             irow) * b_2,
                                                                   ppa_front_th +
                                                                   ind1 * b_2);
                              ppiptr_th[kk]++;
                           }
                        }

                     }

// Compute QR and store new rows scaling

                     int i, njloc, nloc, ibs, kjj, kii;
                     _Flt aux, aux1;

                     for (i = 0; i < niloc; i++) {

                        njloc = ppia_cols[i + 1] - ppia_cols[i];
                        nloc = njloc * _blksize;

                        ibs = ppia_cols[i];

                        CVector < _Flt >::CopyVector (nloc * _blksize,
                                                      ppa_front_th + ibs * b_2, py_blk);

                        CVector < _Flt >::TransposeBlock (_blksize, nloc, py_blk,
                                                          _blksize, pz_blk, nloc);

                        CVector < _Flt >::QrdBlock (_blksize, nloc, pz_blk, nloc,
                                                    ptau_matr);

// Store R factor and its inverse

                        CVector < _Flt >::SetByZeroes (b_2, pr_matr);

                        for (kjj = 0; kjj < _blksize; kjj++) {
                           for (kii = 0; kii <= kjj; kii++) {
                              pr_matr[kjj * _blksize + kii] = pz_blk[kjj * nloc + kii];
                           }
                        }

                        if (iter == _nitersc - 1) {
                           for (kjj = 0; kjj < _blksize; kjj++) {
                              aux = pr_matr[kjj * _blksize + kjj];
                              if (aux < 0.0e0)
                                 aux = -aux;
                              if (psclCmin_thr[my_thr] > aux)
                                 psclCmin_thr[my_thr] = aux;
                              if (psclCmax_thr[my_thr] < aux)
                                 psclCmax_thr[my_thr] = aux;
                           }
                        }

                        for (kjj = 0; kjj < _blksize; kjj++) {
                           aux = pr_matr[kjj * _blksize + kjj];
                           aux1 = aux;
                           if (aux1 < 0.0)
                              aux1 = -aux1;
                           if (aux1 < _sclmin) {
                              if (iter == _nitersc - 1)
                                 pnmodif_thr[my_thr]++;
                              if (aux > 0.0e0) {
                                 pr_matr[kjj * _blksize + kjj] = (_Flt) _sclmin;
                              } else {
                                 pr_matr[kjj * _blksize + kjj] = -(_Flt) _sclmin;
                              }
                           }
                        }

                        CVector < _Flt >::CopyVector (b_2, pr_matr,
                                                      p_sclU_inv + (ibeg + i) * b_2);

                        CBlock_BxB_traits < _Flt >::InvU_BxB (_blksize, pr_matr,
                                                              p_sclU + (ibeg + i) * b_2);

                     }
                  }

               }

               if (iter == _nitersc - 1) {
//                  PrintArray (ffout," p_sclU curr = ",nsuptot*b_2,p_sclU);
//                  PrintArray (ffout," p_sclU inv curr = ",nsuptot*b_2,p_sclU_inv);
               }
// Rows block scaling

               {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < _nblks; ipar++) {

                     int niloc = (int) (_blks[ipar + 1] - _blks[ipar]);

                     int my_thr = 0;
#ifdef USE_THREADS
                     my_thr = omp_get_thread_num ();
#endif

// Memory support

                     if (piptr_thr[my_thr].GetLength () < niblk_max) {
                        piptr_thr[my_thr].resize (niblk_max);
                     }

                     int *ppiptr_th = piptr_thr[my_thr].Ptr ();


                     if (pa_front_thr[my_thr].GetLength () < nzja_max * b_2) {
                        pa_front_thr[my_thr].resize (nzja_max * b_2);
                     }

                     _Flt *ppa_front_th = pa_front_thr[my_thr].Ptr ();

                     int nz_work = b_2 * (ni_max + 1) * 2;

                     if (pwork_front_thr[my_thr].GetLength () < nz_work) {
                        pwork_front_thr[my_thr].resize (nz_work);
                     }

                     _Flt *ppwork_front_th = pwork_front_thr[my_thr].Ptr ();

                     _Flt *py_blk = ppwork_front_th;
                     _Flt *pz_blk = py_blk + ni_max * b_2;
                     _Flt *pr_matr = pz_blk + ni_max * b_2;
                     _Flt *ptau_matr = pr_matr + b_2;

// Compute and store A*U in condensed format

                     int ibeg = (int) _blks[ipar];

                     int *ppia_rows = pia_rows[ipar].Ptr ();

                     int j, ki, kj, jcol, ind, jj, irow;

                     for (j = 0; j < niloc; j++)
                        ppiptr_th[j] = ppia_rows[j];

                     for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
                        int jblk = pja_hmatr[j];
                        int jbeg = (int) _blks[jblk];
                        int nlist_temp = pA_sub[j].GetNlist ();
                        _Int *plist_temp = pA_sub[j].GetListArr ();
                        _Int *pia_temp = pA_sub[j].GetIaArr ();
                        _Int *pja_temp = pA_sub[j].GetJaArr ();
                        _Flt *pa_temp = pA_sub[j].GetAArr ();
                        for (ki = 0; ki < nlist_temp; ki++) {
                           irow = (int) plist_temp[ki];
                           for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
                              jcol = (int) pja_temp[kj];
                              ind = ppiptr_th[irow];
                              jj = jbeg + jcol;
                              CBlock_BxB_traits < _Flt >::MM_BxB (_blksize,
                                                                  pa_temp + kj * b_2,
                                                                  p_sclU + jj * b_2,
                                                                  ppa_front_th +
                                                                  ind * b_2);
                              ppiptr_th[irow]++;
                           }
                        }
                     }

// Compute QR and store new rows scaling

                     int i, njloc, nloc, ibs, kii, kjj;
                     _Flt aux, aux1;

                     for (i = 0; i < niloc; i++) {

                        njloc = ppia_rows[i + 1] - ppia_rows[i];
                        nloc = njloc * _blksize;

                        ibs = ppia_rows[i];

                        CVector < _Flt >::CopyVector (nloc * _blksize,
                                                      ppa_front_th + ibs * b_2, py_blk);

                        CVector < _Flt >::TransposeBlock (_blksize, nloc, py_blk,
                                                          _blksize, pz_blk, nloc);

                        CVector < _Flt >::QrdBlock (_blksize, nloc, pz_blk, nloc,
                                                    ptau_matr);

// Store R factor and its inverse

                        CVector < _Flt >::SetByZeroes (b_2, pr_matr);

                        for (kjj = 0; kjj < _blksize; kjj++) {
                           for (kii = 0; kii <= kjj; kii++) {
                              pr_matr[kjj * _blksize + kii] = pz_blk[kjj * nloc + kii];
                           }
                        }

                        if (iter == _nitersc - 1) {
                           for (kjj = 0; kjj < _blksize; kjj++) {
                              aux = pr_matr[kjj * _blksize + kjj];
                              if (aux < 0.0e0)
                                 aux = -aux;
                              if (psclRmin_thr[my_thr] > aux)
                                 psclRmin_thr[my_thr] = aux;
                              if (psclRmax_thr[my_thr] < aux)
                                 psclRmax_thr[my_thr] = aux;
                           }
                        }

                        for (kjj = 0; kjj < _blksize; kjj++) {
                           aux = pr_matr[kjj * _blksize + kjj];
                           aux1 = aux;
                           if (aux1 < 0.0)
                              aux1 = -aux1;
                           if (aux1 < _sclmin) {
                              if (iter == _nitersc - 1)
                                 pnmodif_thr[my_thr]++;
                              if (aux > 0.0e0) {
                                 pr_matr[kjj * _blksize + kjj] = (_Flt) _sclmin;
                              } else {
                                 pr_matr[kjj * _blksize + kjj] = -(_Flt) _sclmin;
                              }
                           }
                        }

                        CVector < _Flt >::CopyVector (b_2, pr_matr,
                                                      p_sclL_inv + (ibeg + i) * b_2);

                        CBlock_BxB_traits < _Flt >::InvU_BxB (_blksize, pr_matr,
                                                              p_sclL + (ibeg + i) * b_2);

                     }

                  }

               }

               if (iter == _nitersc - 1) {
//                  PrintArray (ffout," p_sclL curr = ",nsuptot*b_2,p_sclL);
//                  PrintArray (ffout," p_sclL inv curr = ",nsuptot*b_2,p_sclL_inv);
               }

            }
         }

         _n_modif = 0;
         _sclmin_att = 1.0e100;
         _sclmax_att = -1.0e100;

         {
            int i;
            for (i = 0; i < n_thr; i++) {
               _n_modif += pnmodif_thr[i];
               if (psclCmin_thr[i] < _sclmin_att)
                  _sclmin_att = psclCmin_thr[i];
               if (psclRmin_thr[i] < _sclmin_att)
                  _sclmin_att = psclRmin_thr[i];
               if (psclCmax_thr[i] > _sclmax_att)
                  _sclmax_att = psclCmax_thr[i];
               if (psclRmax_thr[i] > _sclmax_att)
                  _sclmax_att = psclRmax_thr[i];
            }
         }

         _sclmin_att *= _sclmin_att;
         _sclmax_att *= _sclmax_att;

      }

   }

//
// Perform explicit scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CFctBxBThreads < _Int,
      _Flt >::HMatrixScale_BxB (int _blksize, int _nblks, long long *_blks,
                                CBMatrix < _Int, _Flt > &_a_matr, _Flt * _sclL,
                                _Flt * _sclU)
   {

//      ofstream ffout ("ChkScale.dat");

// Open hmatrix

      int b_2 = _blksize * _blksize;

      CMatrix < int, float >*phmatr = _a_matr.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub = _a_matr.GetASubArr ();

      int *pia_hmatr = phmatr->GetIaArr ();
      int *pja_hmatr = phmatr->GetJaArr ();

#ifdef USE_THREADS
#pragma omp parallel for
#endif
      for (int ipar = 0; ipar < _nblks; ipar++) {

         int ibeg = (int) _blks[ipar];

         int j, jj;

         for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {

            jj = pja_hmatr[j];
            int jbeg = (int) _blks[jj];

            int nlist_temp = pA_sub[j].GetNlist ();

            CFct_bxb_impl < _Int, _Flt >::MatrixScale (_blksize, nlist_temp,
                                                       _sclL + ibeg * b_2,
                                                       _sclU + jbeg * b_2,
                                                       *pA_sub[j].GetList (),
                                                       *pA_sub[j].GetIa (),
                                                       *pA_sub[j].GetJa (),
                                                       *pA_sub[j].GetA ());

         }
      }

//      _a_matr.PrintHMatrix (ffout);

   }

//
// Split LU and combine pairs
//========================================================================================
   template < typename _Int, typename _Flt > void CFctBxBThreads < _Int,
      _Flt >::SplitLUPair_BxB (bool _b_is_char, int _blksize, int _nblks,
                               long long *_blks, CBMatrix < _Int, _Flt > &_a_matr,
                               CBMatrix < _Int, _Flt > &_a_pair)
   {

// Open hmatrix

      int b_2 = _blksize * _blksize;

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

// Check for error

         for (i = 0; i <= _nblks; i++) {
            if (pia_hmatr[i] != pia_hmatrT[i]) {
               cout << " CFctBxBThreads<>::SplitLUPair_BxB: error in block sparsity !" <<
                  endl;
               throw " CFctBxBThreads<>::SplitLUPair_BxB: error in block sparsity !";
            }
         }

         int nzja_temp = pia_hmatr[_nblks];

         for (i = 0; i < nzja_temp; i++) {
            if (pja_hmatr[i] != pja_hmatrT[i]) {
               cout << " CFctBxBThreads<>::SplitLUPair_BxB: error in block sparsity !" <<
                  endl;
               throw " CFctBxBThreads<>::SplitLUPair_BxB: error in block sparsity !";
            }
         }

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

               if (nlist1 != (_blks[ipar + 1] - _blks[ipar])) {
                  cout <<
                     " CFctBxBThreads<>::SplitLUPair_BxB: error in diagonal block sparsity !"
                     << endl;
                  throw
                     " CFctBxBThreads<>::SplitLUPair_BxB: error in diagonal block sparsity !";
               }

               vector < _Int > *pia_lu = pA_sub[ind2].GetIa ();
               vector < _Int > *pja_lu = pA_sub[ind2].GetJa ();
               vector < char >*pjachar_lu = pA_sub[ind2].GetJaChar ();
               vector < _Flt > *pa_lu = pA_sub[ind2].GetA ();

               vector < _Int > ia_l;
               vector < _Int > ja_l;
               vector < char >jachar_l;
               vector < _Flt > a_l;
               vector < _Int > ia_u;
               vector < _Int > ja_u;
               vector < char >jachar_u;
               vector < _Flt > a_u;

               CFct_bxb_impl < _Int, _Flt >::SplitLU (_b_is_char, _blksize, nlist1,
                                                      *pia_lu, *pja_lu, *pjachar_lu,
                                                      *pa_lu, ia_l, ja_l, jachar_l, a_l,
                                                      ia_u, ja_u, jachar_u, a_u);

               int nzja_u = (int) ia_u[nlist1];

               vector < _Int > ia_lt;
               vector < _Int > ja_lt;
               vector < char >jachar_lt;
               vector < _Flt > a_lt;

               CFct_bxb_impl < _Int, _Flt >::Transpose (_b_is_char, _blksize, nlist1,
                                                        ia_l, ja_l, jachar_l, a_l, ia_lt,
                                                        ja_lt, jachar_lt, a_lt);

               {
                  vector < _Int > ia_temp;
                  vector < _Int > ja_temp;
                  vector < char >jachar_temp;
                  vector < _Flt > a_temp;
                  ia_l.swap (ia_temp);
                  ja_l.swap (ja_temp);
                  jachar_l.swap (jachar_temp);
                  a_l.swap (a_temp);
               }

               CMatrix < _Int, _Flt > alu_pair;

               alu_pair.ResizeList (nlist1);

               _Int *plist_pair = alu_pair.GetListArr ();
               vector < _Int > *pia_pair = alu_pair.GetIa ();
               vector < _Int > *pja_pair = alu_pair.GetJa ();
               vector < char >*pjachar_pair = alu_pair.GetJaChar ();
               vector < _Flt > *pa_pair = alu_pair.GetA ();

               for (k = 0; k < nlist1; k++)
                  plist_pair[k] = k;

               CFct_bxb_impl < _Int, _Flt >::CombineLUPairs (_b_is_char, _blksize, nlist1,
                                                             ia_lt, ja_lt, jachar_lt,
                                                             a_lt, ia_u, ja_u, jachar_u,
                                                             a_u, *pia_pair, *pja_pair,
                                                             *pjachar_pair, *pa_pair);

               alu_pair.SetNlist (nlist1);
               alu_pair.SetNzja (nzja_u);
               if (_b_is_char)
                  alu_pair.SetNzjaChar (nzja_u);
               alu_pair.SetNza (2 * nzja_u * b_2);

               pA_sub_pair[j].ReplaceFree (alu_pair);

            } else {

               ind2 = pindpair2matr[j];
               ind1 = pindH2Ht[ind2];

               CMatrix < _Int, _Flt > al_blk;

               pA_sub[ind1].TransposedSparsityList_BxB (_blksize, icycleblk, pimaskblk,
                                                        pimaskblk + nimax,
                                                        pimaskblk + 2 * nimax,
                                                        pimaskblk + 3 * nimax,
                                                        pimaskblk + 4 * nimax, al_blk);

               int nlist1 = al_blk.GetNlist ();
               int nzja1 = al_blk.GetNzja ();
               _Int *plist1 = al_blk.GetListArr ();
               vector < _Int > *pia1 = al_blk.GetIa ();
               vector < _Int > *pja1 = al_blk.GetJa ();
               vector < char >*pjachar1 = al_blk.GetJaChar ();
               vector < _Flt > *pa1 = al_blk.GetA ();

               int nlist2 = pA_sub[ind2].GetNlist ();
               int nzja2 = pA_sub[ind2].GetNzja ();
               vector < _Int > *pia2 = pA_sub[ind2].GetIa ();
               vector < _Int > *pja2 = pA_sub[ind2].GetJa ();
               vector < char >*pjachar2 = pA_sub[ind2].GetJaChar ();
               vector < _Flt > *pa2 = pA_sub[ind2].GetA ();

               if (nlist1 != nlist2 || nzja1 != nzja2) {
                  throw
                     " CFctBxBThreads<_Int,_Flt>::SplitLUPair: error in nlist or nzja !";
               }

               CMatrix < _Int, _Flt > alu_blk;

               alu_blk.ResizeList (nlist1);

               _Int *plist_lu = alu_blk.GetListArr ();
               vector < _Int > *pia_lu = alu_blk.GetIa ();
               vector < _Int > *pja_lu = alu_blk.GetJa ();
               vector < char >*pjachar_lu = alu_blk.GetJaChar ();
               vector < _Flt > *pa_lu = alu_blk.GetA ();

               for (k = 0; k < nlist1; k++)
                  plist_lu[k] = plist1[k];

               CFct_bxb_impl < _Int, _Flt >::CombineLUPairs (_b_is_char, _blksize, nlist1,
                                                             *pia1, *pja1, *pjachar1,
                                                             *pa1, *pia2, *pja2,
                                                             *pjachar2, *pa2, *pia_lu,
                                                             *pja_lu, *pjachar_lu,
                                                             *pa_lu);

               alu_blk.SetNlist (nlist1);
               alu_blk.SetNzja (nzja1);
               if (_b_is_char)
                  alu_blk.SetNzjaChar (nzja1);
               alu_blk.SetNza (2 * nzja1 * b_2);

               pA_sub_pair[j].ReplaceFree (alu_blk);

            }
         }

         picycle_thr[my_thr] = icycleblk;
      }

   }

//
// Split pair into L and U
//========================================================================================
   template < typename _Int, typename _Flt > void CFctBxBThreads < _Int,
      _Flt >::SplitPair_BxB (bool _b_is_char, int _blksize, int _nblks, long long *_blks,
                             CBMatrix < _Int, _Flt > &_a_pair, CBMatrix < _Int,
                             _Flt > &_a_L, CBMatrix < _Int, _Flt > &_a_U)
   {

// Open hmatrix

      int b_2 = _blksize * _blksize;

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
//         double tau1 = -1.0e0;

         int j, k;

         for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {

            int nlist_lu = pA_sub[j].GetNlist ();
            _Int *plist_lu = pA_sub[j].GetListArr ();
            vector < _Int > *pia_lu = pA_sub[j].GetIa ();
            vector < _Int > *pja_lu = pA_sub[j].GetJa ();
            vector < char >*pjachar_lu = pA_sub[j].GetJaChar ();
            vector < _Flt > *pa_lu = pA_sub[j].GetA ();

            pA_sub_L[j].ResizeList (nlist_lu);
            pA_sub_U[j].ResizeList (nlist_lu);

            _Int *plist_l = pA_sub_L[j].GetListArr ();
            vector < _Int > *pia_l = pA_sub_L[j].GetIa ();
            vector < _Int > *pja_l = pA_sub_L[j].GetJa ();
            vector < char >*pjachar_l = pA_sub_L[j].GetJaChar ();
            vector < _Flt > *pa_l = pA_sub_L[j].GetA ();

            _Int *plist_u = pA_sub_U[j].GetListArr ();
            vector < _Int > *pia_u = pA_sub_U[j].GetIa ();
            vector < _Int > *pja_u = pA_sub_U[j].GetJa ();
            vector < char >*pjachar_u = pA_sub_U[j].GetJaChar ();
            vector < _Flt > *pa_u = pA_sub_U[j].GetA ();

            for (k = 0; k < nlist_lu; k++)
               plist_l[k] = plist_lu[k];
            for (k = 0; k < nlist_lu; k++)
               plist_u[k] = plist_lu[k];

            CFct_bxb_impl < _Int, _Flt >::SplitLUPairs (_b_is_char, _blksize, nlist_lu,
                                                        *pia_lu, *pja_lu, *pjachar_lu,
                                                        *pa_lu, *pia_l, *pja_l,
                                                        *pjachar_l, *pa_l, *pia_u, *pja_u,
                                                        *pjachar_u, *pa_u);

            int nzja_L = (int) ((*pia_l)[nlist_lu]);
            int nzja_U = (int) ((*pia_u)[nlist_lu]);

            pA_sub_L[j].SetNlist (nlist_lu);
            pA_sub_L[j].SetNzja (nzja_L);
            if (_b_is_char)
               pA_sub_L[j].SetNzjaChar (nzja_L);
            pA_sub_L[j].SetNza (nzja_L * b_2);

            pA_sub_U[j].SetNlist (nlist_lu);
            pA_sub_U[j].SetNzja (nzja_U);
            if (_b_is_char)
               pA_sub_U[j].SetNzjaChar (nzja_U);
            pA_sub_U[j].SetNza (nzja_U * b_2);

         }
      }

   }

// Rescale U
//========================================================================================
   template < typename _Int, typename _Flt > void CFctBxBThreads < _Int,
      _Flt >::RescaleU_BxB (int _blksize, int _nblks, long long *_blks, CBMatrix < _Int,
                            _Flt > &_U_matr, _Flt * _sclU, _Flt * _inv_sclU)
   {

// Open hmatrix

      int b_2 = _blksize * _blksize;

      CMatrix < int, float >*phmatr_U = _U_matr.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_U = _U_matr.GetASubArr ();

      int *pia_hmatr = phmatr_U->GetIaArr ();
      int *pja_hmatr = phmatr_U->GetJaArr ();

// Rescale data

#ifdef USE_THREADS
#pragma omp parallel for
#endif
      for (int ipar = 0; ipar < _nblks; ipar++) {

         int j, jj;
         int ibeg = (int) _blks[ipar];

         for (j = pia_hmatr[ipar]; j < pia_hmatr[ipar + 1]; j++) {
            jj = pja_hmatr[j];
            int jbeg = (int) _blks[jj];
            if (jj == ipar) {
               int nlist_temp = pA_sub_U[j].GetNlist ();
               CFct_bxb_impl < _Int, _Flt >::RescaleU ('D', _blksize, nlist_temp,
                                                       _sclU + ibeg * b_2,
                                                       _inv_sclU + ibeg * b_2,
                                                       *pA_sub_U[j].GetList (),
                                                       *pA_sub_U[j].GetIa (),
                                                       *pA_sub_U[j].GetJa (),
                                                       *pA_sub_U[j].GetA ());
            } else {
               int nlist_temp = pA_sub_U[j].GetNlist ();
               CFct_bxb_impl < _Int, _Flt >::RescaleU ('O', _blksize, nlist_temp,
                                                       _sclU + ibeg * b_2,
                                                       _inv_sclU + jbeg * b_2,
                                                       *pA_sub_U[j].GetList (),
                                                       *pA_sub_U[j].GetIa (),
                                                       *pA_sub_U[j].GetJa (),
                                                       *pA_sub_U[j].GetA ());
            }
         }
      }

   }

// Set diagonal by zeroes
//========================================================================================
   template < typename _Int, typename _Flt > void CFctBxBThreads < _Int,
      _Flt >::SetZeroDiag_BxB (int _blksize, int _nblks, long long *_blks,
                               CBMatrix < _Int, _Flt > &_A_matr)
   {

// Open hmatrix

      int b_2 = _blksize * _blksize;

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
                        CVector < _Flt >::SetByZeroes (b_2, pa_temp + kj * b_2);
                     }
                  }
               }
            }
         }
      }

   }

// Perform filtering of the pairs block sparsity according to the tree with diagonal modifications
//========================================================================================
   template < typename _Int, typename _Flt > void CFctBxBThreads < _Int,
      _Flt >::PairsFilterTreeModif_BxB (CTree & _tree, double _theta, int _blksize,
                                        int _nblks, long long *_blks, CBMatrix < _Int,
                                        _Flt > &_AU_matr)
   {

// Open hmatrix

      int b_2 = _blksize * _blksize;

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
         pdiamodL_arr[ipar].ResizeAndSetAll (0, nlist_temp, 0, 0, nlist_temp * b_2);
         _Int *plistL_temp = pdiamodL_arr[ipar].GetList2Arr ();
         _Flt *pdiaL_temp = pdiamodL_arr[ipar].GetAArr ();
         for (j = 0; j < nlist_temp; j++)
            plistL_temp[j] = plist_temp[j];
         CVector < _Flt >::SetByZeroes (nlist_temp * b_2, pdiaL_temp);
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
         pdiamodU_arr[ipar].ResizeAndSetAll (0, nlist_diaU, 0, 0, nlist_diaU * b_2);
         _Int *plistU_temp = pdiamodU_arr[ipar].GetList2Arr ();
         _Flt *pdiaU_temp = pdiamodU_arr[ipar].GetAArr ();
         for (j = 0; j < nlist_diaU; j++)
            plistU_temp[j] = plistblk[j];
         CVector < _Flt >::SetByZeroes (nlist_diaU * b_2, pdiaU_temp);
         for (ki = 0; ki < nlist_temp; ki++) {
            for (kj = (int) pia_temp[ki]; kj < pia_temp[ki + 1]; kj++) {
               kk = (int) pja_temp[kj];
               ind = pindblk[kk];
               CBlock_BxB_traits < _Flt >::ModifDia_BxB (_blksize, theta_flt,
                                                         pa_temp + kj * 2 * b_2,
                                                         pa_temp + (kj * 2 + 1) * b_2,
                                                         pdiaL_temp + ki * b_2,
                                                         pdiaU_temp + ind * b_2);
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
                     CBlock_BxB_traits < _Flt >::AddModifDia_BxB (_blksize,
                                                                  pdiaL_temp + k * b_2,
                                                                  pa_temp +
                                                                  ind_dia * 2 * b_2);
                     CBlock_BxB_traits < _Flt >::AddModifDia_BxB (_blksize,
                                                                  pdiaL_temp + k * b_2,
                                                                  pa_temp + (ind_dia * 2 +
                                                                             1) * b_2);
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
                     CBlock_BxB_traits < _Flt >::AddModifDia_BxB (_blksize,
                                                                  pdiaU_temp + k * b_2,
                                                                  pa_temp +
                                                                  ind_dia * 2 * b_2);
                     CBlock_BxB_traits < _Flt >::AddModifDia_BxB (_blksize,
                                                                  pdiaU_temp + k * b_2,
                                                                  pa_temp + (ind_dia * 2 +
                                                                             1) * b_2);
                  }
               }
            }
         }
      }
   }

// Compute reordered by columns rectangular hmatrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFctBxBThreads < _Int,
      _Flt >::ReorderHMatrixColumnsPairs_BxB (int _blksize, int _nblksR,
                                              long long *_blksR, int _nblksC_ini,
                                              long long *_blksC_ini, CBMatrix < _Int,
                                              _Flt > &_hmatr_ini, int *_order,
                                              int _nblksC_fin, long long *_blksC_fin,
                                              CBMatrix < _Int, _Flt > &_hmatr_fin)
   {

// Open hmatr

      int b_2 = _blksize * _blksize;
      int b_2_2 = b_2 * 2;

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
            ASub.ResizeAndSetAll (nlist_curr, 0, nzja_curr, nzja_curr, nzja_curr * b_2_2);
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
                     CVector < _Flt >::CopyVector (b_2_2, pa_temp + kj * b_2_2,
                                                   pa_curr + kkk * b_2_2);
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
               pelems_arr[i].resize (pnzblk_arr[i] * b_2_2 + 1);
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
                  CVector < _Flt >::CopyVector (b_2_2, pa_curr + j * b_2_2,
                                                ppelems_arr + k * b_2_2);
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
            CVectorData < _Flt > elemsarr (nimax * b_2_2 + 1);

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
               pA_sub_new[i].ResizeAndSetAll (nlist_temp, 0, nzja_temp, 0,
                                              nzja_temp * b_2_2);
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
                  CVector < _Flt >::CopyVector (b_2_2, ppelems_arr + j * b_2_2,
                                                pa_temp + k * b_2_2);
                  plistarr[ind]++;
               }
               for (j = 0; j < nlist_temp; j++) {
                  ibeg = (int) pia_temp[j];
                  niloc = (int) (pia_temp[j + 1] - pia_temp[j]);
                  for (k = (int) pia_temp[j]; k < pia_temp[j + 1]; k++) {
                     piiarr[k - ibeg].ival = (int) pja_temp[k];
                     piiarr[k - ibeg].i2val = k;
                     CVector < _Flt >::CopyVector (b_2_2, pa_temp + k * b_2_2,
                                                   pelemsarr + (k - ibeg) * b_2_2);
                  }
                  sort (piiarr, piiarr + niloc);
                  for (k = (int) pia_temp[j]; k < pia_temp[j + 1]; k++) {
                     pja_temp[k] = (_Int) piiarr[k - ibeg].ival;
                     ind = piiarr[k - ibeg].i2val;
                     CVector < _Flt >::CopyVector (b_2_2,
                                                   pelemsarr + (ind - ibeg) * b_2_2,
                                                   pa_temp + k * b_2_2);
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

//========================================================================================
   template < typename _Int, typename _Flt > struct CBlockFctTree_BxB   /// Template structure that computes fct
   {
// Data
      int index2_out;           ///< Secondary index number
      int ordlevel;             ///< Ordering filtration level
      SParams *pparams;         ///< Fct parameters
      int blksize;              ///< Small block size
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
      int *pndiasplit_thr;      ///< The number of diagonal splits
      int *pnmodif_thr;         ///< The number of modifs for threads
      double *peigmin_thr;      ///< Minimal pivots for threads
      double *peigmax_thr;      ///< Maximal pivots for threads
// Functions
      void Execute ()
      {

         int my_thr = 0;

#ifdef USE_THREADS
           my_thr = omp_get_thread_num ();
#endif

         int b_2 = this->blksize * this->blksize;
         int b_2_2 = b_2 * 2;

//         char strbuff[256];
//         if (this->index2_out == 0) {
//            sprintf (strbuff,"ChkFctK_%i_%i.dat",pparams->iparam1,iblk);
//         } else {
//            sprintf (strbuff,"ChkFctK_%i_Sch_%i.dat",pparams->iparam1,iblk);
//         }
//         ofstream ffout (strbuff);

//         cout << " MyThr = " << my_thr << " Ihblk = " << pparams->iparam1 << " Iblk = " << iblk << " Ni = " << pblks[iblk+1]-pblks[iblk] << endl;
//         ffout << " Ihblk = " << pparams->iparam1 << " Iblk = " << iblk << endl;

// Init mask arrays

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
               cout << " CBlockFctTree_BxB<>::Execute: wrong tree !!! " << endl;
               throw " CBlockFctTree_BxB<>::Execute: wrong tree !!! ";
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
               cout << " CBlockFctTree_BxB<>::Execute: wrong tree !!! " << endl;
               throw " CBlockFctTree_BxB<>::Execute: wrong tree !!! ";
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
            hblk_sum.AddReplacePairsBxBThr (false, this->blksize, pschur_arr[ichild2]);
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
//            PrintArray (ffout," ja_curr1 = ",nzja_curr1,pja_curr1);

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
                     asum.AddBlocksPairsBxB ('+', this->blksize, pA_sub_curr0[ip0],
                                             pA_sub_curr1[ip1]);
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

         if (jblk0 != iblk) {
            cout << " CBlockFctTree<>::Execute(): error in block sparsity !" << " Iblk = "
               << iblk << " jblk0 = " << jblk0 << " Ni_diag = " << ni_diag << endl;
            throw " CBlockFctTree<>::Execute(): error in block sparsity ! ";
         }

         vector < _Int > *pia_set = pASub_curr[0].GetIa ();
         vector < _Int > *pja_set = pASub_curr[0].GetJa ();
         vector < char >*pjachar_set = pASub_curr[0].GetJaChar ();
         vector < _Flt > *pa_set = pASub_curr[0].GetA ();

         porder_arr[iblk].resize (ni_diag + 1);
         int *pporder_temp = &porder_arr[iblk][0];

         CFct < _Int, _Flt >::ComputeOptimalOrder (pASub_curr[0], this->pparams->ordtype,
                                                   porder_arr[iblk]);

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

            CFct_bxb_impl < _Int, _Flt >::SplitLUPairs (true, this->blksize, ni_diag,
                                                        *pia_alu_dia, *pja_alu_dia,
                                                        *pjachar_alu_dia, *pa_alu_dia,
                                                        ia_l_dia, ja_l_dia, jachar_l_dia,
                                                        a_l_dia, ia_u_dia, ja_u_dia,
                                                        jachar_u_dia, a_u_dia);

// Transpose L

            vector < _Int > ia_lt_dia;
            vector < _Int > ja_lt_dia;
            vector < char >jachar_lt_dia;
            vector < _Flt > a_lt_dia;

            CFct_bxb_impl < _Int, _Flt >::Transpose (true, this->blksize, ni_diag,
                                                     ia_l_dia, ja_l_dia, jachar_l_dia,
                                                     a_l_dia, ia_lt_dia, ja_lt_dia,
                                                     jachar_lt_dia, a_lt_dia);

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

// Combine small block L+U

            CMatrix < _Int, _Flt > amatr_lu;

            vector < _Int > *pia_lu = amatr_lu.GetIa ();
            vector < _Int > *pja_lu = amatr_lu.GetJa ();
            vector < char >*pjachar_lu = amatr_lu.GetJaChar ();
            vector < _Flt > *pa_lu = amatr_lu.GetA ();

            CFct_bxb_impl < _Int, _Flt >::CombineLU (true, this->blksize, ni_diag,
                                                     ia_lt_dia, ja_lt_dia, jachar_lt_dia,
                                                     a_lt_dia, ia_u_dia, ja_u_dia,
                                                     jachar_u_dia, a_u_dia, *pia_lu,
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

            CMatrix < _Int, _Flt > amatr_ord;

            CFct_bxb < _Int, _Flt >::ReorderMatrix_BxB (true, false, this->blksize,
                                                        amatr_lu, porder_arr[iblk],
                                                        amatr_ord);

            amatr_lu.Clean ();

            vector < _Int > *pia_ord = amatr_ord.GetIa ();
            vector < _Int > *pja_ord = amatr_ord.GetJa ();
            vector < char >*pjachar_ord = amatr_ord.GetJaChar ();
            vector < _Flt > *pa_ord = amatr_ord.GetA ();

// Split

            CFct_bxb_impl < _Int, _Flt >::SplitLU (true, this->blksize, ni_diag, *pia_ord,
                                                   *pja_ord, *pjachar_ord, *pa_ord,
                                                   ia_l_dia, ja_l_dia, jachar_l_dia,
                                                   a_l_dia, ia_u_dia, ja_u_dia,
                                                   jachar_u_dia, a_u_dia);

            amatr_ord.Clean ();

// Transpose L back

            CFct_bxb_impl < _Int, _Flt >::Transpose (true, this->blksize, ni_diag,
                                                     ia_l_dia, ja_l_dia, jachar_l_dia,
                                                     a_l_dia, ia_lt_dia, ja_lt_dia,
                                                     jachar_lt_dia, a_lt_dia);

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

            CFct_bxb_impl < _Int, _Flt >::CombineLUPairs (true, this->blksize, ni_diag,
                                                          ia_lt_dia, ja_lt_dia,
                                                          jachar_lt_dia, a_lt_dia,
                                                          ia_u_dia, ja_u_dia,
                                                          jachar_u_dia, a_u_dia,
                                                          *pia_alu_dia, *pja_alu_dia,
                                                          *pjachar_alu_dia, *pa_alu_dia);

         }

// Reorder rows for all off-diagonal blocks

         {
            int i;
            for (i = 1; i < nzja_curr; i++) {
               CMatrix < _Int, _Flt > aord;
               pASub_curr[i].OrderMtrRowsPairs_BxB (true, this->blksize, pporder_curr,
                                                    aord);
               pASub_curr[i].ReplaceFree (aord);
            }
         }

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
            }

// Compute condensed block as a sum of blocks

            CMatrix < _Int, _Flt > AU_blk;

            AU_blk.ResizeAndSetAll (nlist_col_tot, 0, nzja_tot, 0, 2 * nzja_tot * b_2);

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
                     CVector < _Flt >::CopyVector (b_2_2, pa_temp + k * b_2_2,
                                                   pa_AU_blk + kind * b_2_2);
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

            vector < int >order1 (nlist_col_tot + 1);
            int *porder1 = &order1[0];

            for (i = 0; i < nlist_col_tot; i++)
               porder1[i] = i;

            int n_moved_diag = 0;

            int nlev = 10;
            vector < double >dia_lev (nlev + 1);

            double *pdia_lev = &dia_lev[0];

            for (i = 0; i < nlev; i++)
               pdia_lev[i] = this->pparams->dia_split;

            int nblks_fct = 2;
            vector < int >blks_fct (nlev + 2);

            int *pblks_fct = &blks_fct[0];

            pblks_fct[0] = 0;
            pblks_fct[1] = ni_diag;
            pblks_fct[2] = nlist_col_tot;

            pparams->iparam3 = iblk;

            double tau2_sch_save = this->pparams->tau2_sch;

            this->pparams->tau2_sch = 0.0e0;

            CFct_bxb_impl < _Int, _Flt >::FctMLevDiaSplit (this->pparams, nlev, pdia_lev,
                                                           this->blksize, nlist_col_tot,
                                                           ni_diag, *pia_AU_vect,
                                                           *pja_AU_vect, *pjachar_AU_vect,
                                                           *pa_AU_vect, porder1,
                                                           n_moved_diag, nblks_fct,
                                                           pblks_fct, *pia_U_vect,
                                                           *pja_U_vect, *pjachar_U_vect,
                                                           *pa_U_vect, nmodif, eigmin,
                                                           eigmax);

            this->pparams->tau2_sch = tau2_sch_save;

            pndiasplit_thr[my_thr] += n_moved_diag;

            {
               vector < int >iorder (ni_diag + 1);
               int *piorder = &iorder[0];

               for (i = 0; i < ni_diag; i++)
                  piorder[pporder_curr[i]] = i;

               int iold;

               for (i = 0; i < ni_diag; i++) {
                  iold = piorder[i];
                  pporder_curr[iold] = porder1[i];
               }
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
            U_blk.SetNza (nzja_U_temp * 2 * b_2);

//            ffout << " After FctBlockIlu2Degree: " << endl;
//            U_blk.OutputHead (ffout);

//            ffout << " After Fct: U pair: " << endl;
//            CFct<_Int,_Flt>::PrintMatrix (ffout, U_blk);

//            if (true) {
//               CBMatrix<_Int, _Flt>::Str2PsBox (10, U_blk,
//                                       "OrdU.ps", nblks_fct, pblks_fct);
//            }

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
                  pelems[i].resize (pnzja_blkarr[i] * 2 * b_2);
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
                     CVector < _Flt >::CopyVector (b_2_2, pa_U + k * b_2_2,
                                                   ppelems + kkk * b_2_2);
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
                                              2 * pnzja_blkarr[j] * b_2);
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
                     CVector < _Flt >::CopyVector (b_2_2, ppelems + k * b_2_2,
                                                   pa_temp + kkk * b_2_2);
                     piptr[ind]++;
                  }
                  int njrowmax = 0;
                  for (k = 0; k < nlist_temp; k++) {
                     njrow = (int) (pia_temp[k + 1] - pia_temp[k]);
                     if (njrow > njrowmax)
                        njrowmax = njrow;
                  }
                  CVectorData < CSortInt > iiarr (njrowmax);
                  CVectorData < _Flt > elemsort (njrowmax * 2 * b_2);
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
                        CVector < _Flt >::CopyVector (b_2_2, pa_temp + kkk * b_2_2,
                                                      pelemsort + (kk - ibeg) * b_2_2);
                     }
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        kkk = piiarr[kk - ibeg].ival;
                        pja_temp[kk] = (_Int) kkk;
                        CVector < _Flt >::CopyVector (b_2_2,
                                                      pelemsort + (kk - ibeg) * b_2_2,
                                                      pa_temp + kk * b_2_2);
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

//               ffout << " Schur str: " << endl;
//               PrintArray (ffout," plist_SchurHblk = ",nlistblk_schur,plist_SchurHblk);
//               PrintArray (ffout," plist2_SchurHblk = ",nlistblk_schur,plist2_SchurHblk);
//               PrintArray (ffout," pia_SchurHblk = ",nlistblk_schur+1,pia_SchurHblk);
//               PrintArray (ffout," pja_SchurHblk = ",nzjablk_schur,pja_SchurHblk);
//               PrintArray (ffout," pja2_SchurHblk = ",nzjablk_schur,pja2_SchurHblk);

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
                  pelems[i].resize (pnzja_blkarr[i] * 2 * b_2);
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
                        CVector < _Flt >::CopyVector (b_2_2, pa_U + k * b_2_2,
                                                      ppelems + kkk * b_2_2);
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
                                                  pnzja_blkarr[j] * 2 * b_2);
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
                     CVector < _Flt >::CopyVector (b_2_2, ppelems + k * b_2_2,
                                                   pa_temp + kkk * b_2_2);
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
                  CVectorData < _Flt > elemsort (njrowmax * 2 * b_2);
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
                        CVector < _Flt >::CopyVector (b_2_2, pa_temp + kkk * b_2_2,
                                                      pelemsort + (kk - ibeg) * b_2_2);
                     }
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        kkk = piiarr[kk - ibeg].ival;
                        pja_temp[kk] = (_Int) kkk;
                        pjachar_temp[kk] = pcharsort[kk - ibeg];
                        CVector < _Flt >::CopyVector (b_2_2,
                                                      pelemsort + (kk - ibeg) * b_2_2,
                                                      pa_temp + kk * b_2_2);
                     }
                  }

               }

            }

//            ffout << " =========== Schur hblk fct: " << endl;
//            pschur_arr[iblk].PrintHMatrix (ffout);

         }

// Finally modify Schur data

         pschur_arr[iblk].AddReplacePairsBxBThr (false, this->blksize, hblk_sum);

//         ffout << " =========== Modified Schur hblk fct: " << endl;
//         pschur_arr[iblk].PrintHMatrix (ffout);

         picycle_thr[my_thr] = icycleblk;

      }
   };

//
// Perform ILU2 point factorization of the hmatrix with dynamic ordering and future diagonal modifications
//========================================================================================
   template < typename _Int, typename _Flt > void CFctBxBThreads < _Int,
      _Flt >::Ilu2BlockIlu2_BxB (bool _b_blk_wells, ofstream * _pfout, int _index_out,
                                 SParams * _pparams, int _blksize, int _nblks,
                                 int _nblks1, int _nblks2, long long *_blks,
                                 long long *_nzord_blks, CTree & _tree0, CTree & _tree1,
                                 CBMatrix < _Int, _Flt > &_alu_pair, CTree & _tree2_new,
                                 int &_nblks_new, vector < long long >&_blks_new,
                                 CVectorData < int >&_ordernew, CBMatrix < _Int,
                                 _Flt > &_lu_pair, int &_ndiasplit, int &_nmodif,
                                 double &_eigmin_att, double &_eigmax_att)
   {

      int ordlevel = _pparams->ordlevel;

//      char strbuff[256];
//      sprintf (strbuff,"ChkFct_%i.dat",_index_out);
//      ofstream ffout (strbuff);
//      _alu_pair.PrintHMatrix (ffout);

// Open hmatrix

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

      vector < int >ndiasplit_thr (n_thr + 1);
      vector < int >nmodif_thr (n_thr + 1);
      vector < double >eigmin_thr (n_thr + 1);
      vector < double >eigmax_thr (n_thr + 1);

      int *pndiasplit_thr = &ndiasplit_thr[0];
      int *pnmodif_thr = &nmodif_thr[0];
      double *peigmin_thr = &eigmin_thr[0];
      double *peigmax_thr = &eigmax_thr[0];

      int i;

      for (i = 0; i < n_thr; i++)
         pndiasplit_thr[i] = 0;
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

      _pparams->iparam2 = 0;

      vector < CBlockFctTree_BxB < _Int, _Flt > >fct_12 (_nblks + 1);
      CBlockFctTree_BxB < _Int, _Flt > *pfct_12 = &fct_12[0];

      for (i = 0; i < _nblks1 + _nblks2; i++) {

         CBlockFctTree_BxB < _Int, _Flt > *task = pfct_12 + i;

         task->index2_out = 0;
         task->ordlevel = _pparams->ordlevel;
         task->pparams = _pparams;
         task->blksize = _blksize;
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
         task->pndiasplit_thr = pndiasplit_thr;
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

      for (ilev = nlev_1 - 1; ilev >= nlev_12_min; ilev--) {

         nnodes_curr_1 = pnnodes_lev_1[ilev];
         ppnodeslevlist_1 = &pnodeslevlist_1[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_1; ipar++) {
            int iblk_temp = ppnodeslevlist_1[ipar];
            pfct_12[iblk_temp].Execute ();
         }

      }

      for (ilev = nlev_2 - 1; ilev >= nlev_12_min; ilev--) {

         nnodes_curr_2 = pnnodes_lev_2[ilev];
         ppnodeslevlist_2 = &pnodeslevlist_2[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_2; ipar++) {
            int iblk_temp = _nblks1 + ppnodeslevlist_2[ipar];
            pfct_12[iblk_temp].Execute ();
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

      {

         CBMatrix < _Int, _Flt > alu_sum;

         if (_nblks1 != 0) {
//            ffout << " Schur 1: " << endl;
//            pschur_arr[_nblks1 - 1].PrintHMatrix (ffout);
            alu_sum.AddReplacePairsBxBThr (true, _blksize, pschur_arr[_nblks1 - 1]);
         }

         if (_nblks2 != 0) {
//            ffout << " Schur 2: " << endl;
//            pschur_arr[nblks12 - 1].PrintHMatrix (ffout);
            alu_sum.AddReplacePairsBxBThr (true, _blksize, pschur_arr[nblks12 - 1]);
         }
//         ffout << " Schur sum: " << endl;
//         alu_sum.PrintHMatrix (ffout);

//         ffout << " Before Schur add: " << endl;
//         alu_pair_last.PrintHMatrix (ffout);

         alu_pair_last.AddReplacePairsBxBThr (true, _blksize, alu_sum);

//         ffout << " Schur add: " << endl;
//         alu_pair_last.PrintHMatrix (ffout);

      }

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
                            || ((ordlevel >= 0 && jj_char <= ordlevel)
                                || ordlevel < 0)) {
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
                            || ((ordlevel >= 0 && jj_char <= ordlevel)
                                || ordlevel < 0)) {
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
                        jj_char = pjachar_temp[kj];
                        if (((ipar == jj) && (kii == kjj))
                            || ((ordlevel >= 0 && jj_char <= ordlevel)
                                || ordlevel < 0)) {
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

//      for (i = 0; i < ntot_3; i++) porder3[i] = i;  /// Check !!!!!!!!!!!!!!

      for (i = 0; i < ntot_3; i++)
         piorder3[porder3[i]] = i;

// Reorder and repartition submatrix

// Split pair hblock

      CBMatrix < _Int, _Flt > al_last;
      CBMatrix < _Int, _Flt > au_last;

      CFctBxBThreads < _Int, _Flt >::SplitPair_BxB (true, _blksize, nblks3, pblks3,
                                                    alu_pair_last, al_last, au_last);

      alu_pair_last.Clean ();

// Transpose hblock

      CBMatrix < _Int, _Flt > alt_last;

      CBMatrix < _Int, _Flt >::TransposeHMatrix_BxB (_blksize, nblks3, pblks3, al_last,
                                                     alt_last);

      al_last.Clean ();

// Set diagonal by zero

      CFctBxBThreads < _Int, _Flt >::SetZeroDiag_BxB (_blksize, nblks3, pblks3, alt_last);

// Combine into hmatrix

      CBMatrix < _Int, _Flt > alu_last;

      alu_last.ReplaceFree (alt_last);

      alu_last.AddReplaceBxBThr (true, _blksize, au_last);

      alt_last.Clean ();
      au_last.Clean ();

// Compute reordered hmatrix

      CBMatrix < _Int, _Flt > alu_last_ord;

      CBMatrix < _Int, _Flt >::ReorderHMatrix_BxB (_blksize, nblks3, pblks3, alu_last,
                                                   porder3, nblks_tree3, pblks_tree3,
                                                   alu_last_ord);

      alu_last.Clean ();

// Split into pairs

      CBMatrix < _Int, _Flt > alu_pair_last_ord;

      CFctBxBThreads < _Int, _Flt >::SplitLUPair_BxB (true, _blksize, nblks_tree3,
                                                      pblks_tree3, alu_last_ord,
                                                      alu_pair_last_ord);

      alu_last_ord.Clean ();

// Perform filtering of the block sparsity according to the tree with diagonal modifications

      CFctBxBThreads < _Int, _Flt >::PairsFilterTreeModif_BxB (_tree2_new,
                                                               _pparams->theta, _blksize,
                                                               nblks_tree3, pblks_tree3,
                                                               alu_pair_last_ord);

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

      vector < CBlockFctTree_BxB < _Int, _Flt > >fct_3 (nblks_tree3 + 1);
      CBlockFctTree_BxB < _Int, _Flt > *pfct_3 = &fct_3[0];

      _pparams->iparam2 = 1;

      for (i = 0; i < nblks_tree3; i++) {

         CBlockFctTree_BxB < _Int, _Flt > *task = pfct_3 + i;

         task->index2_out = 1;
         task->ordlevel = _pparams->ordlevel;
         task->pparams = _pparams;
         task->blksize = _blksize;
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
         task->pndiasplit_thr = pndiasplit_thr;
         task->pnmodif_thr = pnmodif_thr;
         task->peigmin_thr = peigmin_thr;
         task->peigmax_thr = peigmax_thr;

      }

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

         CFctBxBThreads < _Int, _Flt >::ReorderHMatrixColumnsPairs_BxB (_blksize, nblks12,
                                                                        _blks, nblks3,
                                                                        pblks3,
                                                                        hblock_12_3,
                                                                        porder3,
                                                                        nblks_tree3,
                                                                        pblks_tree3,
                                                                        hblock_12_3_ord);

      }

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

                  pASub_fin[j].OrderMtrColsPairs_BxB (true, _blksize, porder_temp, aord);

                  pASub_fin[j].ReplaceFree (aord);
               }
            }
         }
      }

// Form final result

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

      _ndiasplit = 0;
      _nmodif = 0;
      _eigmin_att = 1.0e100;
      _eigmax_att = -1.0e100;

      for (i = 0; i < n_thr; i++)
         _ndiasplit += pndiasplit_thr[i];
      for (i = 0; i < n_thr; i++)
         _nmodif += pnmodif_thr[i];
      for (i = 0; i < n_thr; i++) {
         if (peigmin_thr[i] < _eigmin_att)
            _eigmin_att = peigmin_thr[i];
         if (peigmax_thr[i] > _eigmax_att)
            _eigmax_att = peigmax_thr[i];
      }

   }

//========================================================================================
   template < typename _Int, typename _Flt > struct CBlockFctSymmTree_BxB       /// Template structure that computes symmetric fct
   {
// Data
      SParams *pparams;         ///< Fct parameters
      int blksize;              ///< Small block size
      int nblks;                ///< Total number of blocks
      long long *pblks;         ///< Block partitioning
      int inode;                ///< Extended node number
      int *pnode2blk;           ///< Block numbers for nodes
      int *pnode2dep;           ///< Dependency of the nodes
      int *pia_AU;              ///< Ia sparsity of AU
      int *pja_AU;              ///< Ja sparsity of AU
        CMatrix < _Int, _Flt > *pASub_AU;       ///< Blocks of AU
        CBMatrix < _Int, _Flt > *pblockrowsU_arr;       ///< Computed U
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

         int my_thr = 0;

#ifdef USE_THREADS
           my_thr = omp_get_thread_num ();
#endif

         int b_2 = this->blksize * this->blksize;

//         char strbuff[256];
//         if (this->index2_out == 0) {
//            sprintf (strbuff,"ChkFctK_%i_%i.dat",pparams->iparam1,iblk);
//         } else {
//            sprintf (strbuff,"ChkFctK_%i_Sch_%i.dat",pparams->iparam1,iblk);
//         }
//         ofstream ffout (strbuff);

//         cout << " MyThr = " << my_thr << " Ihblk = " << pparams->iparam1 << " Iblk = " << iblk << " Ni = " << pblks[iblk+1]-pblks[iblk] << endl;
//         ffout << " Ihblk = " << pparams->iparam1 << " Iblk = " << iblk << endl;

// Init mask arrays

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

// Update by the childs

         if (pnode2dep[inode * 2] >= 0) {
            int jnode = pnode2dep[inode * 2];
            pschur_arr[inode].AddReplacePairsBxBThr (false, this->blksize,
                                                     pschur_arr[jnode]);
            pschur_arr[jnode].Clean ();
         }
         if (pnode2dep[inode * 2 + 1] >= 0) {
            int jnode = pnode2dep[inode * 2 + 1];
            pschur_arr[inode].AddReplacePairsBxBThr (false, this->blksize,
                                                     pschur_arr[jnode]);
            pschur_arr[jnode].Clean ();
         }
// Factorize current block

         int iblk = pnode2blk[inode];

         if (iblk < 0)
            return;

         int ni_diag = (int) (pblks[iblk + 1] - pblks[iblk]);

         CBMatrix < _Int, _Flt > hblk_sum;

         hblk_sum.ReplaceFree (pschur_arr[inode]);

//      ffout << " =========== Hblk_Schur_sum: " << endl;
//      hblk_sum.PrintHMatrix (ffout);

// Get set of blocks that correspond to current block row and split hblock for future add

         int nzja_curr = 0;
         vector < int >ja_curr;
         vector < CMatrix < _Int, _Flt > >ASub_curr;

         {

            int nzja_curr0 = pia_AU[iblk + 1] - pia_AU[iblk];
            int *pja_curr0 = pja_AU + pia_AU[iblk];
            CMatrix < _Int, _Flt > *pA_sub_curr0 = pASub_AU + pia_AU[iblk];

//         PrintArray (ffout," ja_curr0 = ",nzja_curr0,pja_curr0);

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

            ja_curr.resize (nzja_curr0 + nzja_curr1);
            ASub_curr.resize (nzja_curr0 + nzja_curr1);

            int *pja_curr = ja_curr.data ();
            CMatrix < _Int, _Flt > *pASub_curr = ASub_curr.data ();

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
                     asum.AddBlocksBxB ('+', this->blksize, pA_sub_curr0[ip0],
                                        pA_sub_curr1[ip1]);
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

         int *pja_curr = ja_curr.data ();
         CMatrix < _Int, _Flt > *pASub_curr = ASub_curr.data ();

//      ffout << " Iblk = " << iblk << endl;
//      PrintArray (ffout," ja_curr = ",nzja_curr,pja_curr);

//      ffout << " =========== Hblk_Schur_sum flt: " << endl;
//       hblk_sum.PrintHMatrix (ffout);

//      ffout << " Point 11 check " << endl;

// Perform main fct computations

         {

// Create condense mask array

            int ni_max = 0;

            int nj_blk = 0;

            int i, jblk, ni;

            for (i = 0; i < nzja_curr; i++) {
               jblk = pja_curr[i];
               ni = (int) (pblks[jblk + 1] - pblks[jblk]);
               nj_blk += ni;
               if (ni > ni_max)
                  ni_max = ni;
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

//         ffout << " Point 12 check " << endl;

// Compute condensed block as a sum of blocks

            CMatrix < _Int, _Flt > AU_blk;

            AU_blk.ResizeAndSetAll (nlist_col_tot, 0, nzja_tot, 0, nzja_tot * b_2);

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
                     CVector < _Flt >::CopyVector (b_2, pa_temp + k * b_2,
                                                   pa_AU_blk + kind * b_2);
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
//            ffout << "   fcttype_sch = " << this->pparams->fcttype_sch << " tau2_sch = " <<  this->pparams->tau2_sch << endl;
//            CFct<_Int,_Flt>::PrintMatrix (ffout, AU_blk);

//            ofstream *pfout_save = this->pparams->pfout;

//            char strbuff1[256];
//            ofstream *pfout_debug = NULL;
//            if (index2_out == 0) {
//               sprintf (strbuff1,"ChkIluDeg_%i_%i.dat",pparams->iparam1,iblk);
//            } else {
//               sprintf (strbuff1,"ChkIluDeg_%i_Sch_%i.dat",pparams->iparam1,iblk);
//            }
//            ofstream fout_debug (strbuff1);
//            pfout_debug = &fout_debug;
//            this->pparams->pfout = pfout_debug;

            pparams->iparam3 = iblk;

            double tau2_sch_save = this->pparams->tau2_sch;

            this->pparams->tau2_sch = 0.0e0;

            CFct_bxb_impl < _Int, _Flt >::FctBlockIch2Degree (this->pparams,
                                                              this->blksize,
                                                              nlist_col_tot, ni_diag,
                                                              *pia_AU_vect, *pja_AU_vect,
                                                              *pjachar_AU_vect,
                                                              *pa_AU_vect, *pia_U_vect,
                                                              *pja_U_vect,
                                                              *pjachar_U_vect, *pa_U_vect,
                                                              nmodif, eigmin, eigmax);

            this->pparams->tau2_sch = tau2_sch_save;

//         this->pparams->pfout = pfout_save;

//         this->pparams->pfout = NULL;

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
            U_blk.SetNza (nzja_U_temp * b_2);

//            ffout << " After FctBlockIlu2Degree: " << endl;
//            U_blk.OutputHead (ffout);

//            ffout << " After Fct: U pair: " << endl;
//            CFct<_Int,_Flt>::PrintMatrix (ffout, U_blk);

//            if (true) {
//               CBMatrix<_Int, _Flt>::Str2PsBox (10, U_blk,
//                                       "OrdU.ps", nblks_fct, pblks_fct);
//            }

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

               pblockrowsU_arr[iblk].ResizeASub (nzjablk_U);
               pblockrowsU_arr[iblk].SetNzblk (nzjablk_U);

               CMatrix < _Int, _Flt > *pASub_U = pblockrowsU_arr[iblk].GetASubArr ();

               CMatrix < int, float >*pHMatrU = pblockrowsU_arr[iblk].GetHMatrStr ();

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
                  pelems[i].resize (pnzja_blkarr[i] * b_2);
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
                     CVector < _Flt >::CopyVector (b_2, pa_U + k * b_2,
                                                   ppelems + kkk * b_2);
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
                                              pnzja_blkarr[j] * b_2);
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
                     CVector < _Flt >::CopyVector (b_2, ppelems + k * b_2,
                                                   pa_temp + kkk * b_2);
                     piptr[ind]++;
                  }
                  int njrowmax = 0;
                  for (k = 0; k < nlist_temp; k++) {
                     njrow = (int) (pia_temp[k + 1] - pia_temp[k]);
                     if (njrow > njrowmax)
                        njrowmax = njrow;
                  }
                  CVectorData < CSortInt > iiarr (njrowmax);
                  CVectorData < _Flt > elemsort (njrowmax * b_2);
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
                        CVector < _Flt >::CopyVector (b_2, pa_temp + kkk * b_2,
                                                      pelemsort + (kk - ibeg) * b_2);
                     }
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        kkk = piiarr[kk - ibeg].ival;
                        pja_temp[kk] = (_Int) kkk;
                        CVector < _Flt >::CopyVector (b_2, pelemsort + (kk - ibeg) * b_2,
                                                      pa_temp + kk * b_2);
                     }
                  }
               }

            } else {

               pblockrowsU_arr[iblk].ResizeASub (1);
               pblockrowsU_arr[iblk].SetNzblk (1);

               CMatrix < int, float >*pHMatrU = pblockrowsU_arr[iblk].GetHMatrStr ();

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
               CVectorData < CVectorData < _Flt > >elems (nzjablk_schur);

               CVectorData < _Int > *pirows = irows.Ptr ();
               CVectorData < _Int > *picols = icols.Ptr ();
               CVectorData < char >*pchars = chars.Ptr ();
               CVectorData < _Flt > *pelems = elems.Ptr ();

               for (i = 0; i < nzjablk_schur; i++) {
                  pirows[i].resize (pnzja_blkarr[i]);
                  picols[i].resize (pnzja_blkarr[i]);
                  pchars[i].resize (pnzja_blkarr[i]);
                  pelems[i].resize (pnzja_blkarr[i] * b_2);
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
                        CVector < _Flt >::CopyVector (b_2, pa_U + k * b_2,
                                                      ppelems + kkk * b_2);
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
                                                  pnzja_blkarr[j] * b_2);
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
                     CVector < _Flt >::CopyVector (b_2, ppelems + k * b_2,
                                                   pa_temp + kkk * b_2);
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
                  CVectorData < _Flt > elemsort (njrowmax * b_2);
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
                        CVector < _Flt >::CopyVector (b_2, pa_temp + kkk * b_2,
                                                      pelemsort + (kk - ibeg) * b_2);
                     }
                     for (kk = (int) pia_temp[k]; kk < pia_temp[k + 1]; kk++) {
                        kkk = piiarr[kk - ibeg].ival;
                        pja_temp[kk] = (_Int) kkk;
                        pjachar_temp[kk] = pcharsort[kk - ibeg];
                        CVector < _Flt >::CopyVector (b_2, pelemsort + (kk - ibeg) * b_2,
                                                      pa_temp + kk * b_2);
                     }
                  }

               }

            }

//            ffout << " =========== Schur hblk fct: " << endl;
//            pschur_arr[iblk].PrintHMatrix (ffout);

// Finally modify Schur data

            pschur_arr[inode].AddReplacePairsBxBThr (false, this->blksize, hblk_sum);

//         ffout << " =========== Modified Schur hblk fct: " << endl;
//         pschur_arr[iblk].PrintHMatrix (ffout);

         }

         picycle_thr[my_thr] = icycleblk;

      }
   };

//
// Perform ICH2 block bxb factorization of the hmatrix with dinamic ordering and future diagonal modifications
//========================================================================================
   template < typename _Int, typename _Flt > void CFctBxBThreads < _Int,
      _Flt >::Ich2HMatrix_BxB (SParams & _params, int _blksize, int _nhblks,
                               int *_hblk2blks, int _nblks, long long *_blks,
                               CTree & _tree, CTree * _tree_arr, CBMatrix < _Int,
                               _Flt > &_au_matr, CBMatrix < _Int, _Flt > &_u_matr,
                               int &_nmodif, double &_eigmin_att, double &_eigmax_att)
   {

//   char strbuff[256];
//   sprintf (strbuff,"ChkFct_%i.dat",_index_out);
//   ofstream ffout (strbuff);
//   _alu_pair.PrintHMatrix (ffout);

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

// Open hmatrix

      CMatrix < int, float >*phmatr_AU = _au_matr.GetHMatrStr ();
      CMatrix < _Int, _Flt > *pA_sub_AU = _au_matr.GetASubArr ();

      int *pia_hmatr_AU = phmatr_AU->GetIaArr ();
      int *pja_hmatr_AU = phmatr_AU->GetJaArr ();

// Create mask arrays

      vector < int >icycle_thr (n_thr + 1);
      vector < CVectorData < int > >imask_thr (n_thr + 1);

      int *picycle_thr = icycle_thr.data ();
      CVectorData < int >*pimask_thr = imask_thr.data ();

      int nimax = 0;

      {

         int niloc;

         int i;

         for (i = 0; i < _nblks; i++) {
            niloc = (int) (_blks[i + 1] - _blks[i]);
            if (niloc > nimax)
               nimax = niloc;
         }

         for (i = 0; i < n_thr; i++)
            picycle_thr[i] = -1;

      }

// Compute arrays that support trees level by level computations 

      int nnodes_tree = _tree.GetNnodes ();
      int nlev_tree = _tree.GetNlev ();
      int *pnode2lev = _tree.GetNodesLevId ();
      int *pnnodeslev = _tree.GetNNodesLev ();
      vector < int >*pnodes_lev_list = _tree.GetNodesLevList ();
      int *pnchilds = _tree.GetNchilds ();
      vector < int >*pchilds_list = _tree.GetChildsList ();

      if (nnodes_tree != _nhblks) {
         cout << " CFctBxBThreads<>::Ich2HMatrix_BxB: error in hblks tree ! " << endl;
         throw " CFctBxBThreads<>::Ich2HMatrix_BxB: error in hblks tree ! ";
      }

      CVectorData < int >hblk2blks_ext (_nhblks + 1);
      CVectorData < int >lev2sublev (nlev_tree + 1);

      int *phblk2blks_ext = hblk2blks_ext.Ptr ();
      int *plev2sublev = lev2sublev.Ptr ();

      int nblks_ext = 0;

      {
         int i, ilev, ilev1;
         phblk2blks_ext[0] = 0;
         for (i = 0; i < _nhblks; i++) {
            nblks_ext += _tree_arr[i].GetNnodes ();
            phblk2blks_ext[i + 1] = nblks_ext;
         }
         plev2sublev[0] = 0;
         for (i = 1; i <= nlev_tree; i++)
            plev2sublev[i] = 1;
         for (i = 0; i < _nhblks; i++) {
            ilev = pnode2lev[i];
            ilev1 = nlev_tree - ilev - 1;
            int nlev_temp = _tree_arr[i].GetNlev ();
            if (plev2sublev[ilev1 + 1] < nlev_temp)
               plev2sublev[ilev1 + 1] = nlev_temp;
         }
         for (i = 0; i < nlev_tree; i++)
            plev2sublev[i + 1] += plev2sublev[i];
      }

      int nlev_ext = plev2sublev[nlev_tree];

      CVectorData < int >blkini2nodes_ext (_nblks);
      CVectorData < int >blk2nodes_ext (nblks_ext * 3);
      CVectorData < int >nodeext2blk (nblks_ext * 3);
      CVectorData < int >blk2lev_ext (nblks_ext);
      CVectorData < int >blk2dep_ext (nblks_ext * 2);

      int *pblkini2nodes_ext = blkini2nodes_ext.Ptr ();
      int *pblk2nodes_ext = blk2nodes_ext.Ptr ();
      int *pnodeext2blk = nodeext2blk.Ptr ();
      int *pblk2lev_ext = blk2lev_ext.Ptr ();
      int *pblk2dep_ext = blk2dep_ext.Ptr ();

      {
         int i, jj, j, ilev, ilev1, jlev, jlev1, k, ichild;
         for (i = 0; i < nblks_ext * 2; i++)
            pblk2dep_ext[i] = -1;
         nblks_ext = 0;
         for (i = 0; i < _nhblks; i++) {
            ilev = pnode2lev[i];
            ilev1 = nlev_tree - ilev - 1;
            int nnodes_temp = _tree_arr[i].GetNnodes ();
            int nlev_temp = _tree_arr[i].GetNlev ();
            int *pind2node_temp = _tree_arr[i].GetInd2Node ();
            int *pnode2lev_temp = _tree_arr[i].GetNodesLevId ();
            int *pnchilds_temp = _tree_arr[i].GetNchilds ();
            vector < int >*pchilds_list_temp = _tree_arr[i].GetChildsList ();
            int nblks_temp = _hblk2blks[i + 1] - _hblk2blks[i];
            for (j = 0; j < nnodes_temp; j++) {
               pblk2nodes_ext[(nblks_ext + j) * 3] = j;
               pblk2nodes_ext[(nblks_ext + j) * 3 + 1] = i;
               pblk2nodes_ext[(nblks_ext + j) * 3 + 2] = -1;
               pnodeext2blk[nblks_ext + j] = -1;
            }
            int ibegblk = _hblk2blks[i];
            for (j = 0; j < nblks_temp; j++) {
               jj = pind2node_temp[j];
               pblk2nodes_ext[(nblks_ext + jj) * 3 + 2] = ibegblk + j;
               pnodeext2blk[nblks_ext + jj] = ibegblk + j;
               pblkini2nodes_ext[ibegblk + j] = nblks_ext + jj;
            }
            for (j = 0; j < nnodes_temp; j++) {
               jlev = pnode2lev_temp[j];
               jlev1 = nlev_temp - jlev - 1;
               pblk2lev_ext[nblks_ext + j] = plev2sublev[ilev1] + jlev1;
            }
            for (j = 0; j < nnodes_temp; j++) {
               int nchilds = pnchilds_temp[j];
               int *pchilds = pchilds_list_temp[j].data ();
               for (k = 0; k < nchilds; k++) {
                  ichild = pchilds[k];
                  if (ichild != j) {
                     pblk2dep_ext[(nblks_ext + j) * 2 + k] = nblks_ext + ichild;
                  }
               }
            }
            nblks_ext += nnodes_temp;
         }
      }

// Create sublists of extended nodes for each extended node

      vector < int >ia_sublev (nlev_ext + 1);
      vector < int >ja_sublev (nblks_ext);

      int *pia_sublev = ia_sublev.data ();
      int *pja_sublev = ja_sublev.data ();

      {
         int i, jlev;
         for (i = 0; i <= nlev_ext; i++)
            pia_sublev[i] = 0;
         for (i = 0; i < nblks_ext; i++) {
            jlev = pblk2lev_ext[i];
            pia_sublev[jlev + 1]++;
         }
         for (i = 0; i < nlev_ext; i++)
            pia_sublev[i + 1] += pia_sublev[i];
         vector < int >iptr (nlev_ext);
         int *piptr = iptr.data ();
         for (i = 0; i < nlev_ext; i++)
            piptr[i] = pia_sublev[i];
         int k;
         for (i = 0; i < nblks_ext; i++) {
            jlev = pblk2lev_ext[i];
            k = piptr[jlev];
            pja_sublev[k] = i;
            piptr[jlev]++;
         }
      }

// Main cycle over hblocks

      vector < int >n_modif_thr (n_thr);
      vector < double >eigmin_thr (n_thr);
      vector < double >eigmax_thr (n_thr);

      int *pn_modif_thr = n_modif_thr.data ();
      double *peigmin_thr = eigmin_thr.data ();
      double *peigmax_thr = eigmax_thr.data ();

      {
         int j;
         for (j = 0; j < n_thr; j++)
            pn_modif_thr[j] = 0;
         for (j = 0; j < n_thr; j++)
            peigmin_thr[j] = 1.0e100;
         for (j = 0; j < n_thr; j++)
            peigmax_thr[j] = -1.0e100;
      }

      vector < CBMatrix < _Int, _Flt > >blockrowsU_arr (_nblks + 1);
      vector < CBMatrix < _Int, _Flt > >schur_arr (nblks_ext + 1);

      CBMatrix < _Int, _Flt > *pblockrowsU_arr = blockrowsU_arr.data ();
      CBMatrix < _Int, _Flt > *pschur_arr = schur_arr.data ();

      {

// Cycle over tree levels

         int ilev_gl;

         for (ilev_gl = nlev_tree - 1; ilev_gl >= 0; ilev_gl--) {

// Add in parallel Schur complements for hblocks

            int nnodes_lev = pnnodeslev[ilev_gl];
            int *ppnodes_lev_list = pnodes_lev_list[ilev_gl].data ();

            int ilev_gl1 = nlev_tree - 1 - ilev_gl;

            {

               if (nnodes_lev >= n_thr) {

#ifdef USE_THREADS
#pragma omp parallel for
#endif
                  for (int ipar = 0; ipar < nnodes_lev; ipar++) {
                     int inode = ppnodes_lev_list[ipar];
                     int inode_ext = phblk2blks_ext[inode + 1] - 1;
                     int nchilds = pnchilds[inode];
                     int *pchilds = pchilds_list[inode].data ();
                     int k, ichild, jnode, jnode_ext;
                     for (k = 0; k < nchilds; k++) {
                        ichild = pchilds[k];
                        if (ichild != inode) {
                           jnode = ichild;
                           jnode_ext = phblk2blks_ext[jnode + 1] - 1;
                           pschur_arr[inode_ext].AddReplaceBxBThr (false, _blksize,
                                                                   pschur_arr[jnode_ext]);
                        }
                     }
                  }

               } else {
                  for (int ipar = 0; ipar < nnodes_lev; ipar++) {
                     int inode = ppnodes_lev_list[ipar];
                     int inode_ext = phblk2blks_ext[inode + 1] - 1;
                     int nchilds = pnchilds[inode];
                     int *pchilds = pchilds_list[inode].data ();
                     int k, ichild, jnode, jnode_ext;
                     for (k = 0; k < nchilds; k++) {
                        ichild = pchilds[k];
                        if (ichild != inode) {
                           jnode = ichild;
                           jnode_ext = phblk2blks_ext[jnode + 1] - 1;
                           pschur_arr[inode_ext].AddReplaceBxBThr (true, _blksize,
                                                                   pschur_arr[jnode_ext]);
                        }
                     }
                  }
               }

            }

// Move some block rows of the current hblocks into the corresponding Schur data

            {

               for (int ipar = 0; ipar < nnodes_lev; ipar++) {
                  int inode = ppnodes_lev_list[ipar];
                  int inode_ext = phblk2blks_ext[inode + 1] - 1;
                  int iblkbeg = _hblk2blks[inode];
                  int iblkend = _hblk2blks[inode + 1] - 1;
                  if (iblkbeg < iblkend) {
                     CMatrix < int, float >*phmatr_hblk =
                        pschur_arr[inode_ext].GetHMatrStr ();
                     CMatrix < _Int, _Flt > *pASub_hblk =
                        pschur_arr[inode_ext].GetASubArr ();
                     int nlist_hblk = phmatr_hblk->GetNlist ();
                     int *plist_hblk = phmatr_hblk->GetListArr ();
                     int *pia_hblk = phmatr_hblk->GetIaArr ();
                     int *pja_hblk = phmatr_hblk->GetJaArr ();
                     int nzja_hblk = pia_hblk[nlist_hblk];
                     int nlist_last = 0;
                     int nzja_last = 0;
                     int j, jblk;
                     for (j = 0; j < nlist_hblk; j++) {
                        jblk = plist_hblk[j];
                        if (jblk >= iblkend) {
                           nlist_last++;
                           nzja_last += (pia_hblk[j + 1] - pia_hblk[j]);
                        }
                     }
                     int nlist_beg = nlist_hblk - nlist_last;
                     int nzja_beg = nzja_hblk - nzja_last;
                     int j1, nzja_beg_temp, jnode_ext;
                     for (j1 = 0; j1 < nlist_beg; j1++) {
                        nzja_beg_temp = pia_hblk[j1];
                        jblk = plist_hblk[j1];
                        jnode_ext = pblkini2nodes_ext[jblk];
                        int nzja_temp = (pia_hblk[j1 + 1] - pia_hblk[j1]);
                        CBMatrix < _Int, _Flt > hmatr_temp;
                        hmatr_temp.ResizeASub (nzja_temp);
                        hmatr_temp.SetNzblk (nzja_temp);
                        CMatrix < int, float >*phmatr_temp = hmatr_temp.GetHMatrStr ();
                        CMatrix < _Int, _Flt > *pASub_temp = hmatr_temp.GetASubArr ();
                        phmatr_temp->ResizeAndSetAllSp (1, 1, nzja_temp, nzja_temp);
                        int *plist_temp = phmatr_temp->GetListArr ();
                        int *plist2_temp = phmatr_temp->GetList2Arr ();
                        int *pia_temp = phmatr_temp->GetIaArr ();
                        int *pja_temp = phmatr_temp->GetJaArr ();
                        int *pja2_temp = phmatr_temp->GetJa2Arr ();
                        plist_temp[0] = jblk;
                        plist2_temp[0] = 0;
                        pia_temp[0] = 0;
                        pia_temp[1] = nzja_temp;
                        for (j = 0; j < nzja_temp; j++)
                           pja_temp[j] = pja_hblk[nzja_beg_temp + j];
                        for (j = 0; j < nzja_temp; j++)
                           pja2_temp[j] = 0;
                        for (j = 0; j < nzja_temp; j++)
                           pASub_temp[j].ReplaceFree (pASub_hblk[nzja_beg_temp + j]);
                        pschur_arr[jnode_ext].ReplaceFree (hmatr_temp);
                     }
                     CBMatrix < _Int, _Flt > hmatr_last;
                     hmatr_last.ResizeASub (nzja_last);
                     hmatr_last.SetNzblk (nzja_last);
                     CMatrix < int, float >*phmatr_last = hmatr_last.GetHMatrStr ();
                     CMatrix < _Int, _Flt > *pASub_last = hmatr_last.GetASubArr ();
                     phmatr_last->ResizeAndSetAllSp (nlist_last, nlist_last, nzja_last,
                                                     nzja_last);
                     int *plist_last = phmatr_last->GetListArr ();
                     int *plist2_last = phmatr_last->GetList2Arr ();
                     int *pia_last = phmatr_last->GetIaArr ();
                     int *pja_last = phmatr_last->GetJaArr ();
                     int *pja2_last = phmatr_last->GetJa2Arr ();
                     for (j = 0; j < nlist_last; j++)
                        plist_last[j] = plist_hblk[nlist_beg + j];
                     for (j = 0; j < nlist_last; j++)
                        plist2_last[j] = 0;
                     for (j = 0; j <= nlist_last; j++)
                        pia_last[j] = pia_hblk[nlist_beg + j] - pia_hblk[nlist_beg];
                     for (j = 0; j < nzja_last; j++)
                        pja_last[j] = pja_hblk[nzja_beg + j];
                     for (j = 0; j < nzja_last; j++)
                        pja2_last[j] = 0;
                     for (j = 0; j < nzja_last; j++)
                        pASub_last[j].ReplaceFree (pASub_hblk[nzja_beg + j]);
                     pschur_arr[inode_ext].ReplaceFree (hmatr_last);
                  }
               }

            }

// Cycle over sublevels

            int isublev;

            for (isublev = plev2sublev[ilev_gl1]; isublev < plev2sublev[ilev_gl1 + 1];
                 isublev++) {

               int nnodes_ext_lev = pia_sublev[isublev + 1] - pia_sublev[isublev];

               int ibeg_sublev = pia_sublev[isublev];

               vector < CBlockFctSymmTree_BxB < _Int, _Flt > >fctsymm (nnodes_ext_lev);
               CBlockFctSymmTree_BxB < _Int, _Flt > *pfctsymm = fctsymm.data ();

               int j;

               for (j = 0; j < nnodes_ext_lev; j++) {
                  pfctsymm[j].pparams = &_params;
                  pfctsymm[j].blksize = _blksize;
                  pfctsymm[j].nblks = _nblks;
                  pfctsymm[j].pblks = _blks;
                  pfctsymm[j].inode = pja_sublev[ibeg_sublev + j];
                  pfctsymm[j].pnode2blk = pnodeext2blk;
                  pfctsymm[j].pnode2dep = pblk2dep_ext;
                  pfctsymm[j].pia_AU = pia_hmatr_AU;
                  pfctsymm[j].pja_AU = pja_hmatr_AU;
                  pfctsymm[j].pASub_AU = pA_sub_AU;
                  pfctsymm[j].pblockrowsU_arr = pblockrowsU_arr;
                  pfctsymm[j].pschur_arr = pschur_arr;
                  pfctsymm[j].nimax = nimax;
                  pfctsymm[j].picycle_thr = picycle_thr;
                  pfctsymm[j].pimask_thr = pimask_thr;
                  pfctsymm[j].pnmodif_thr = pn_modif_thr;
                  pfctsymm[j].peigmin_thr = peigmin_thr;
                  pfctsymm[j].peigmax_thr = peigmax_thr;
               }

#ifdef USE_THREADS
#pragma omp parallel for
#endif
               for (int ipar = 0; ipar < nnodes_ext_lev; ipar++) {
                  pfctsymm[j].Execute ();
               }

            }

         }

      }


// Allocate all necessary work arrays
/*
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

      vector < int >ndiasplit_thr (n_thr + 1);
      vector < int >nmodif_thr (n_thr + 1);
      vector < double >eigmin_thr (n_thr + 1);
      vector < double >eigmax_thr (n_thr + 1);

      int *pndiasplit_thr = &ndiasplit_thr[0];
      int *pnmodif_thr = &nmodif_thr[0];
      double *peigmin_thr = &eigmin_thr[0];
      double *peigmax_thr = &eigmax_thr[0];

      int i;

      for (i = 0; i < n_thr; i++)
         pndiasplit_thr[i] = 0;
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

      _pparams->iparam2 = 0;

      vector < CBlockFctTree_BxB < _Int, _Flt > >fct_12 (_nblks + 1);
      CBlockFctTree_BxB < _Int, _Flt > *pfct_12 = &fct_12[0];

      for (i = 0; i < _nblks1 + _nblks2; i++) {

         CBlockFctTree_BxB < _Int, _Flt > *task = pfct_12 + i;

         task->index2_out = 0;
         task->ordlevel = _pparams->ordlevel;
         task->pparams = _pparams;
         task->blksize = _blksize;
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
         task->pndiasplit_thr = pndiasplit_thr;
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

      for (ilev = nlev_1 - 1; ilev >= nlev_12_min; ilev--) {

         nnodes_curr_1 = pnnodes_lev_1[ilev];
         ppnodeslevlist_1 = &pnodeslevlist_1[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_1; ipar++) {
            int iblk_temp = ppnodeslevlist_1[ipar];
            pfct_12[iblk_temp].Execute ();
         }

      }

      for (ilev = nlev_2 - 1; ilev >= nlev_12_min; ilev--) {

         nnodes_curr_2 = pnnodes_lev_2[ilev];
         ppnodeslevlist_2 = &pnodeslevlist_2[ilev][0];

#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < nnodes_curr_2; ipar++) {
            int iblk_temp = _nblks1 + ppnodeslevlist_2[ipar];
            pfct_12[iblk_temp].Execute ();
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

      if (_nblks1 != 0) {
         alu_pair_last.AddReplacePairsBxBThr (true, _blksize, pschur_arr[_nblks1 - 1]);
      }

      if (_nblks2 != 0) {
         alu_pair_last.AddReplacePairsBxBThr (true, _blksize, pschur_arr[nblks12 - 1]);
      }

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
*/
   }

// Perform instatiation of classes with template types 

   template class CBlock_BxB_traits < float >;
   template class CBlock_BxB_traits < double >;

   template class CFct_bxb_impl < int, float >;
   template class CFct_bxb_impl < int, double >;

   template class CFct_bxb < int, float >;
   template class CFct_bxb < int, double >;

   template class CMvmSlv_BxB_impl < long long, double, double >;
   template class CMvmSlv_BxB_impl < int, double, double >;
   template class CMvmSlv_BxB_impl < long long, float, double >;
   template class CMvmSlv_BxB_impl < int, float, double >;
   template class CMvmSlv_BxB_impl < long long, float, float >;
   template class CMvmSlv_BxB_impl < int, float, float >;

   template class CMvmSlv_BxB < long long, double, double >;
   template class CMvmSlv_BxB < int, double, double >;
   template class CMvmSlv_BxB < long long, float, double >;
   template class CMvmSlv_BxB < int, float, double >;
   template class CMvmSlv_BxB < long long, float, float >;
   template class CMvmSlv_BxB < int, float, float >;

   template class CFctBxBThreads < long long, double >;
   template class CFctBxBThreads < int, double >;
   template class CFctBxBThreads < long long, float >;
   template class CFctBxBThreads < int, float >;

   template class CMvmParBxBThreads < long long, double, double >;
   template class CMvmParBxBThreads < int, double, double >;
   template class CMvmParBxBThreads < long long, float, double >;
   template class CMvmParBxBThreads < int, float, double >;

   template class CSlvParBxBThreads < long long, double, double >;
   template class CSlvParBxBThreads < int, double, double >;
   template class CSlvParBxBThreads < long long, float, double >;
   template class CSlvParBxBThreads < int, float, double >;

   template class CK3D_SolverBxBThreads < long long, double, double >;
   template class CK3D_SolverBxBThreads < int, double, double >;
   template class CK3D_SolverBxBThreads < long long, float, double >;
   template class CK3D_SolverBxBThreads < int, float, double >;

}                               // namespace k3d

#endif
