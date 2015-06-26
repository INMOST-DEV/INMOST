//------------------------------------------------------------------------------------------------
// File: k3d.cpp
//------------------------------------------------------------------------------------------------
#include "k3d.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iosfwd>
#include <cstring>
#include <iomanip>
#include <algorithm>
#include <float.h>

#include <mpi.h>

namespace k3d {

/// @brief Perform binary search
//========================================================================================
int BinarySearch (long long _isup, int _nblks, long long *_blks, int _iblkprev)
{
   int iblkprev = _iblkprev;
   if (iblkprev < 0) iblkprev = 0;
   if (iblkprev >= _nblks) {
      iblkprev = _nblks-1;
   }
   if (_isup >= _blks[_iblkprev] && _isup < _blks[_iblkprev+1]) {
      return iblkprev;
   }
   int ibegblk = 0;
   int iendblk = _nblks;
   int iblk;
   while (true) {
      if (ibegblk == iendblk-1) {
         iblkprev = ibegblk;
         break;
      } else {
         iblk = (ibegblk+iendblk-1) / 2;
         if (iblk <= ibegblk) iblk = ibegblk+1;
         if (_isup >= _blks[iblk] && _isup < _blks[iblk+1]) {
            iblkprev = iblk;
            break;
         } else if (_isup < _blks[iblk]) {
            iendblk = iblk;
         } else {
            ibegblk = iblk+1;
         }
      }
   }
   return iblkprev;
};

/// @brief Print array
//========================================================================================
void PrintArray (ostream &_stream, const char *_name, int _isize, int *_iarr) {
   _stream << _name << " (Size = " << _isize << ")" << endl;
   for (int i=0;i<_isize;i++) {
      if ((i%5 == 0) && (i > 0)) _stream << endl;
      _stream << setw(12) << _iarr[i];
   }
   _stream << endl;
};

/// @brief Print array
//========================================================================================
void PrintArray (ostream &_stream, const char *_name, int _isize, long long *_iarr) {
   _stream << _name << " (Size = " << _isize << ")" << endl;
   for (int i=0;i<_isize;i++) {
      if ((i%5 == 0) && (i > 0)) _stream << endl;
      _stream << setw(18) << _iarr[i];
   }
   _stream << endl;
};

/// @brief Print array
//========================================================================================
void PrintArray (ostream &_stream, const char *_name, int _isize, float *_farr) {
   _stream << _name << " (Size = " << _isize << ")" << endl;
   for (int i=0;i<_isize;i++) {
      if ((i%5 == 0) && (i > 0)) _stream << endl;
      _stream << setw(14) << setprecision(8) << _farr[i] << " ";
   }
   _stream << endl;
};

/// @brief Print array
//========================================================================================
void PrintArray (ostream &_stream, const char *_name, int _isize, double *_darr) {
   _stream << _name << " (Size = " << _isize << ")" << endl;
   for (int i=0;i<_isize;i++) {
      if ((i%3 == 0) && (i > 0)) _stream << endl;
      _stream << setw(23) << setprecision(16) << _darr[i] << " ";
   }
   _stream << endl;
};

/// @brief Set manip
//========================================================================================
ostream &SetPw (ostream &stream) {
   stream.unsetf(ios::scientific);
   stream << ' ' << setprecision(5) << setw (9);
   return stream;
};

/// @brief Round two values
//========================================================================================
void Round (double &_dx, double &_dy) {

   long long i, i0;

   double dx1 = 1.0e1 * _dx;

   i0 = (long long) dx1;
   i = (long long) ((dx1-(double)i0) * 10000);
   dx1 = (double) i;
   dx1 = dx1 * 1.0e-4+(double)i0;

   _dx = dx1 * 1.e-1;

   double dy1 = 1.0e1 * _dy;

   i0 = (long long) dy1;
   i = (long long) ((dy1-(double)i0) * 10000);
   dy1 = (double) i;
   dy1 = dy1 * 1.0e-4+(double)i0;

   _dy = dy1 * 1.e-1;

};

// Compute optimal ordering for the matrix
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::ComputeOptimalOrder (int _ordtype,
                                                      const int _n, const vector<_Int> &_ia_alu, const vector<_Int> &_ja_alu,
                                                      vector<int> &_order) 
{

// Compute symmetrized matrix

   vector<_Int> ia_symm;
   vector<_Int> ja_symm;

   CIlu2_impl<_Int,_Flt>::SymmetrizeSparsity (_n, _ia_alu, _ja_alu, 
                                                   ia_symm, ja_symm);

// Call METIS ordering routine

   _order.resize (_n+1);

   int nzjas = (int)ia_symm[_n];

   if (_ordtype == -1) {
/*
      vector<idx_t> ialoc (_n+1);
      vector<idx_t> jaloc(nzjas+1);

      idx_t *pialoc = &ialoc[0];
      idx_t *pjaloc = &jaloc[0];

      idx_t i, j, jj;

      ialoc[0] = 0;
      idx_t nz = 0;

      for (i=0;i<_n;i++) {
         for (j=ia_symm[i];j<ia_symm[i+1];j++) {
            jj = ja_symm[j];
            if (jj != i) {
               pjaloc[nz] = jj;
               nz++;
            };
         };
         pialoc[i+1] = nz;
      };

      vector<idx_t> order_int (_n*2+1);

      idx_t *porder_int = &order_int[0];

      idx_t options_ND[METIS_NOPTIONS];

      METIS_SetDefaultOptions (options_ND);

      options_ND[METIS_OPTION_NUMBERING] = 0;
      options_ND[METIS_OPTION_SEED] = 0;

      idx_t nloc = _n;

      METIS_NodeND ((idx_t *)&nloc, (idx_t *)pialoc, (idx_t *)pjaloc, NULL,
                     (idx_t *)options_ND, (idx_t *)porder_int+nloc, (idx_t *)porder_int);

      for (i=0;i<_n;i++) {
         _order[i] = (int)porder_int[i];
      }
*/
   } else if (_ordtype == 0) {

// Identity ordering

      int i;

      for (i=0;i<_n;i++) {
         _order[i] = i;
      }

   } else if (_ordtype == 1) {

// Reversed Cuthill-McCee ordering

      vector<int> ialoc(_n+1);
      vector<int> jaloc(nzjas+1);

      int *pialoc = &ialoc[0];
      int *pjaloc = &jaloc[0];

      int i;

      for (i=0;i<=_n;i++) pialoc[i] = (int)(ia_symm[i] + 1);
      for (i=0;i<nzjas;i++) pjaloc[i] = (int)(ja_symm[i] + 1);

      vector<int> xls(_n+1);
      vector<int> mask(_n+1);

      int *pxls = &xls[0];
      int *pmask = &mask[0];

      vector<int> order_int(_n+1);

      int *porder_int = &order_int[0];

      void genrcm(int n, int * xadj, int * iadj, int * perm, int * xls, int * mask);

      genrcm (_n, pialoc, pjaloc, porder_int, pxls, pmask);

      for (i=0;i<_n;i++) {
         porder_int[i] -= 1;
      }

      for (i=0;i<_n;i++) {
         pmask[porder_int[i]] = i;
      }

      for (i=0;i<_n;i++) {
         _order[i] = pmask[i];
      }

   }

};

//
// Compute optimal ordering for the matrix via splitting into 3 blocks: 1-st order Schur complement for third block, which is the separator between first block and remaining data
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::ComputeOptimalOrderSchur (int _ordtype,
                                                            const int _n, int _n1, const vector<_Int> &_ia_alu, const vector<_Int> &_ja_alu,
                                                            int &_n2, vector<int> &_order) 
{

   _order.resize (_n+1);

   int *porder = &_order[0];

// Compute symmetrized matrix

   vector<_Int> ia_symm;
   vector<_Int> ja_symm;

   CIlu2_impl<_Int,_Flt>::SymmetrizeSparsity (_n, _ia_alu, _ja_alu, 
                                                   ia_symm, ja_symm);

   _Int *pia_symm = &ia_symm[0];
   _Int *pja_symm = &ja_symm[0];

// Get submatrix corresponding to the first block and compute optimal ordering for it

   int nzja_1 = 0;

   int i, j, jj;

   for (i=0;i<_n1;i++) {
      for (j=(int)pia_symm[i];j<pia_symm[i+1];j++) {
         jj = (int)pja_symm[j];
         if (jj < _n1) nzja_1++;
      }
   }

   vector<_Int> ia1 (_n1+1);
   vector<_Int> ja1 (nzja_1+1);

   _Int *pia1 = &ia1[0];
   _Int *pja1 = &ja1[0];

   nzja_1 = 0;
   pia1[0] = 0;

   for (i=0;i<_n1;i++) {
      for (j=(int)pia_symm[i];j<pia_symm[i+1];j++) {
         jj = (int)pja_symm[j];
         if (jj < _n1) {
            pja1[nzja_1] = jj;
            nzja_1++;
         }
      }
      pia1[i+1] = nzja_1;
   }

   vector<int> order1 (1);

   CIlu2_impl<_Int,_Flt>::ComputeOptimalOrder (_ordtype,
                                                      _n1, ia1, ja1,
                                                      order1);

   int *porder1 = &order1[0];

   for (i=0;i<_n1;i++) porder[i] = porder1[i];

   {
      vector<_Int> ia1_temp(1);
      vector<_Int> ja1_temp(1);
      vector<int> order1_temp(1);
      ia1.swap(ia1_temp);
      ja1.swap(ja1_temp);
      order1.swap(order1);
   }

// Register columns in the second part which is to be moved to the bordering

   vector<int> imask (_n+1);
   vector<int> list (_n+1);
   vector<int> indarr (_n+1);

   int *pimask = &imask[0];
   int *plist = &list[0];
   int *pindarr = &indarr[0];

   for (i=_n1;i<_n;i++) pimask[i] = -1;

   int nlist_bord = 0;

   for (i=0;i<_n1;i++) {
      for (j=(int)pia_symm[i];j<pia_symm[i+1];j++) {
         jj = (int)pja_symm[j];
         if (jj >= _n1 && pimask[jj] == -1) {
            plist[nlist_bord] = jj;
            nlist_bord++;
            pimask[jj] = 1;
         }
      }
   }

   sort (plist,plist+nlist_bord);

   _n2 = _n-nlist_bord;

   int ip2 = _n1;
   int ip3 = _n2;

   for (i=_n1;i<_n;i++) {
      if (pimask[i] == -1) {
         porder[i] = ip2;
         ip2++;
      } else {
         porder[i] = ip3;
         ip3++;
      }
   }

   vector<int> iorder (_n+1);

   int *piorder = &iorder[0];

   for (i=0;i<_n;i++) piorder[porder[i]] = i;

// Reorder the matrix

   vector<_Int> ia_ord(1);
   vector<_Int> ja_ord(1);

   CIlu2_impl<_Int,_Flt>::ReorderMatrixSp (_n, _order,
                                                ia_symm, ja_symm,
                                                ia_ord, ja_ord);

   _Int *pia_ord = &ia_ord[0];
   _Int *pja_ord = &ja_ord[0];

// Get submatrix 2 and compute its ordering

   int nzja_2 = 0;

   for (i=_n1;i<_n2;i++) {
      for (j=(int)pia_ord[i];j<pia_ord[i+1];j++) {
         jj = (int)pja_ord[j];
         if (jj >= _n1 && jj < _n2) nzja_2++;
      }
   }

   int n12 = _n2-_n1;

   vector<_Int> ia2 (n12+1);
   vector<_Int> ja2 (nzja_2+1);

   _Int *pia2 = &ia2[0];
   _Int *pja2 = &ja2[0];

   nzja_2 = 0;
   pia2[0] = 0;

   for (i=_n1;i<_n2;i++) {
      for (j=(int)pia_ord[i];j<pia_ord[i+1];j++) {
         jj = (int)pja_ord[j];
         if (jj >= _n1 && jj < _n2) {
            pja2[nzja_2] = jj-_n1;
            nzja_2++;
         }
      }
      pia2[i-_n1+1] = nzja_2;
   }

   vector<int> order2 (1);

   CIlu2_impl<_Int,_Flt>::ComputeOptimalOrder (_ordtype,
                                                      n12, ia2, ja2,
                                                      order2);

   int *porder2 = &order2[0];

   int iold;

   for (i=0;i<n12;i++) {
      iold = piorder[i+_n1];
      porder[iold] = porder2[i]+_n1;
   }

   {
      vector<_Int> ia2_temp (1);
      vector<_Int> ja2_temp (1);
      vector<int> order2_temp (1);
      ia2.swap(ia2_temp);
      ja2.swap(ja2_temp);
      order2.swap(order2);
   }

// Compute ends of second block in all rows

   vector<int> ia_end2 (_n+1);

   int *pia_end2 = &ia_end2[0];

   for (i=0;i<_n;i++) {
      pia_end2[i] = (int)pia_ord[i];
      for (j=(int)pia_ord[i];j<pia_ord[i+1];j++) {
         jj = (int)pja_ord[j];
         if (jj < _n2) pia_end2[i] = j;
      }
   }

// Compute extended bordering data

   int n23 = _n-_n2;

   for (i=0;i<_n;i++) pimask[i] = -1;

   int icycle = -1;

   vector<int> ia_bord (n23+1);

   int *pia_bord = &ia_bord[0];

   pia_bord[0] = 0;
   int nzja_bord = 0;

   vector<int> ja_bord (1);

   int nlistloc, k, kk;

   for (i=_n2;i<_n;i++) {
      icycle++;
      nlistloc = 0;
      for (j=(int)pia_ord[i];j<=pia_end2[i];j++) {
         jj = (int)pja_ord[j];
         for (k=(int)pia_ord[jj];k<=pia_end2[jj];k++) {
            kk = (int)pja_ord[k];
            if (pimask[kk] != icycle) {
               plist[nlistloc] = kk;
               nlistloc++;
               pimask[kk] = icycle;
            }
         }
      }
      sort (plist,plist+nlistloc);
      ja_bord.resize (nzja_bord+nlistloc+1);
      for (j=0;j<nlistloc;j++) ja_bord[nzja_bord+j] = plist[j];
      nzja_bord += nlistloc;
      pia_bord[i-_n2+1] = nzja_bord;
   }

   int *pja_bord = &ja_bord[0];

// Condense, renumber and extend

   icycle++;
   nlistloc = 0;
   for (j=0;j<nzja_bord;j++) {
      jj = pja_bord[j];
      if (pimask[jj] != icycle) {
         plist[nlistloc] = jj;
         nlistloc++;
         pimask[jj] = icycle;
      }
   }
   sort (plist,plist+nlistloc);
   for (j=0;j<nlistloc;j++) {
      jj = plist[j];
      pindarr[jj] = j;
   }
   for (j=0;j<nzja_bord;j++) {
      jj = pja_bord[j];
      pja_bord[j] = pindarr[jj];
   }

   ia_bord.resize (nlistloc+1);

   pia_bord = &ia_bord[0];

   for (j=n23+1;j<=nlistloc;j++) pia_bord[j] = nzja_bord;

// Transpose

   vector<int> ia_bord_t (1);
   vector<int> ja_bord_t (1);

   CIlu2_impl<int,_Flt>::TransposeSp (nlistloc, 
                                          ia_bord, ja_bord,
                                          ia_bord_t, ja_bord_t);

   int *pia_bord_t = &ia_bord_t[0];
   int *pja_bord_t = &ja_bord_t[0];

// Compute Schur complement via extended block data

   vector<_Int> ia3_Schur (n23+1);
   vector<_Int> ja3_Schur (1);

   _Int *pia3_Schur = &ia3_Schur[0];

   int nzja_Schur = 0;

   for (i=0;i<n23;i++) {
      icycle++;
      nlistloc = 0;
      for (j=pia_bord[i];j<pia_bord[i+1];j++) {
         jj = pja_bord[j];
         for (k=pia_bord_t[jj];k<pia_bord_t[jj+1];k++) {
            kk = pja_bord_t[k];
            if (pimask[kk] != icycle) {
               plist[nlistloc] = kk;
               nlistloc++;
               pimask[kk] = icycle;
            }
         }
      }
      sort (plist,plist+nlistloc);
      ja3_Schur.resize (nzja_Schur+nlistloc+1);
      for (j=0;j<nlistloc;j++) ja3_Schur[nzja_Schur+j] = plist[j];
      nzja_Schur += nlistloc;
      pia3_Schur[i+1] = nzja_Schur;
   }

   _Int *pja3_Schur = &ja3_Schur[0];

// Add submatrices

   vector<_Int> ia3 (n23+1);
   vector<_Int> ja3 (1);

   int nzja_3 = 0;

   _Int *pia3 = &ia3[0];

   for (i=0;i<n23;i++) {
      icycle++;
      nlistloc = 0;
      for (j=(int)pia3_Schur[i];j<pia3_Schur[i+1];j++) {
         jj = (int)pja3_Schur[j];
         plist[nlistloc] = jj;
         nlistloc++;
         pimask[jj] = icycle;
      }
      for (j=pia_end2[i+_n2]+1;j<pia_ord[i+_n2+1];j++) {
         jj = (int)pja_ord[j]-_n2;
         if (pimask[jj] != icycle) {
            plist[nlistloc] = jj;
            nlistloc++;
            pimask[jj] = icycle;
         }
      }
      sort (plist,plist+nlistloc);
      ja3.resize (nzja_3+nlistloc+1);
      for (j=0;j<nlistloc;j++) ja3[nzja_3+j] = plist[j];
      nzja_3 += nlistloc;
      pia3[i+1] = nzja_3;
   }

// Optimal order

   vector<int> order3 (1);

   CIlu2_impl<_Int,_Flt>::ComputeOptimalOrder (_ordtype,
                                                      n23, ia3, ja3,
                                                      order3);

   int *porder3 = &order3[0];

   for (i=0;i<n23;i++) {
      iold = piorder[i+_n2];
      porder[iold] = porder3[i]+_n2;
   }

};

//
// Compute ordered matrix (sparsity only)
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::ReorderMatrixSp (int _n, const vector<int> &_order,
                                                vector<_Int> &_ia_alu, vector<_Int> &_ja_alu,
                                                vector<_Int> &_ia_alu_ord, vector<_Int> &_ja_alu_ord)
{

// Compute inverse order

   vector<int> iord (_n+1);

   int *piord = &iord[0];

   int i;

   for (i=0;i<_n;i++) piord[_order[i]] = i;

// Reorder the matrix

   int nzja = (int)_ia_alu[_n];

   _ia_alu_ord.resize (_n+1);
   _ja_alu_ord.resize (nzja+1);

   _Int *piaord = &_ia_alu_ord[0];
   _Int *pjaord = &_ja_alu_ord[0];

   for (i=0;i<=_n;i++) piaord[i] = 0;

   int nimax = 0;

   int ni;

   for (i=0;i<_n;i++) {
      ni = (int)(_ia_alu[i+1]-_ia_alu[i]);
      if (ni>nimax) nimax = ni;
   }

   vector<CSortInt> iiarr (nimax+1);

   CSortInt *piiarr = &iiarr[0];

   int j;

   for (i=0;i<_n;i++) {
      j = _order[i];
      piaord[j+1] = _ia_alu[i+1]-_ia_alu[i];
   }

   for (i=0;i<_n;i++) piaord[i+1] += piaord[i];

   int nz = 0;

   int inew, nzloc, ibs, jold, jnew;

   for (inew=0;inew<_n;inew++) {
      i = piord[inew];
      nzloc = (int)(_ia_alu[i+1]-_ia_alu[i]);
      ibs = nz;
      for (j=(int)_ia_alu[i];j<_ia_alu[i+1];j++) {
         jold = (int)_ja_alu[j];
         jnew = _order[jold];
         pjaord[nz] = (_Int)jnew;
         nz++;
      }

// Sort elements

      for (j=0;j<nzloc;j++) {
         piiarr[j].ival = (int)pjaord[ibs+j];
         piiarr[j].i2val = j;
      }

      sort(piiarr,piiarr+nzloc);

      for (j=0;j<nzloc;j++) {
         pjaord[ibs+j] = (_Int)piiarr[j].ival;
      }

   }

};

//
// Compute ordered matrix
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::ReorderMatrix (int _n, const vector<int> &_order,
                                                vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                                                vector<_Int> &_ia_alu_ord, vector<_Int> &_ja_alu_ord, vector<_Flt> &_a_alu_ord)
{

// Compute inverse order

   vector<int> iord (_n+1);

   int *piord = &iord[0];

   int i;

   for (i=0;i<_n;i++) piord[_order[i]] = i;

// Reorder the matrix

   int nzja = (int)_ia_alu[_n];

   _ia_alu_ord.resize (_n+1);
   _ja_alu_ord.resize (nzja+1);
   _a_alu_ord.resize (nzja+1);

   _Int *piaord = &_ia_alu_ord[0];
   _Int *pjaord = &_ja_alu_ord[0];
   _Flt *paord = &_a_alu_ord[0];

   for (i=0;i<=_n;i++) piaord[i] = 0;

   int nimax = 0;

   int ni;

   for (i=0;i<_n;i++) {
      ni = (int)(_ia_alu[i+1]-_ia_alu[i]);
      if (ni>nimax) nimax = ni;
   }

   vector<CSortInt> iiarr (nimax+1);
   vector<_Flt> elems (nimax+1);

   CSortInt *piiarr = &iiarr[0];
   _Flt *pelems = &elems[0];

   int j;

   for (i=0;i<_n;i++) {
      j = _order[i];
      piaord[j+1] = _ia_alu[i+1]-_ia_alu[i];
   }

   for (i=0;i<_n;i++) piaord[i+1] += piaord[i];

   int nz = 0;

   int inew, nzloc, ibs, jold, jnew;

   for (inew=0;inew<_n;inew++) {
      i = piord[inew];
      nzloc = (int)(_ia_alu[i+1]-_ia_alu[i]);
      ibs = nz;
      for (j=(int)_ia_alu[i];j<_ia_alu[i+1];j++) {
         jold = (int)_ja_alu[j];
         jnew = _order[jold];
         pjaord[nz] = (_Int)jnew;
         paord[nz] = _a_alu[j];
         nz++;
      }

// Sort elements

      for (j=0;j<nzloc;j++) {
         piiarr[j].ival = (int)pjaord[ibs+j];
         piiarr[j].i2val = j;
      }

      sort(piiarr,piiarr+nzloc);

      for (j=0;j<nzloc;j++) pelems[j] = paord[ibs+j];

      for (j=0;j<nzloc;j++) {
         pjaord[ibs+j] = (_Int)piiarr[j].ival;
         paord[ibs+j] = pelems[piiarr[j].i2val];
      }

   }

};

//
// Perform ILU2 point factorization of the block with future diagonal modification
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::Ilu2BlockTransform (int _sctype, int _nitersc, int _fcttype, double _pivmin, double _tau1, double _tau2, double _theta,
                                                   int _n, 
                                                   vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                                                   vector<_Int> &_ia_lu, vector<_Int> &_ida_lu, 
                                                   vector<_Int> &_ja_lu, vector<_Flt> &_a_lu,
                                                   double &_eigmin_att, double &_eigmax_att) 
{

// Perform fct

   vector<_Int> ia_lt;
   vector<_Int> ja_lt;
   vector<_Flt> a_lt;
   vector<_Int> ia_u;
   vector<_Int> ja_u;
   vector<_Flt> a_u;

   CIlu2_impl<_Int,_Flt>::Ilu2Block (_sctype, _nitersc, _fcttype, _pivmin, _tau1, _tau2, _theta,
                                          _n, 
                                          _ia_alu, _ja_alu, _a_alu,
                                          ia_lt, ja_lt, a_lt,
                                          ia_u, ja_u, a_u,
                                          _eigmin_att, _eigmax_att);

// Transpose L

   vector<_Int> ia_l;
   vector<_Int> ja_l;
   vector<_Flt> a_l;

   CIlu2_impl<_Int,_Flt>::Transpose (_n,
                                          ia_lt, ja_lt, a_lt,
                                          ia_l, ja_l, a_l);

   {
      vector<_Int> ia_dummy;
      vector<_Int> ja_dummy;
      vector<_Flt> a_dummy;
      ia_lt.swap (ia_dummy);
      ja_lt.swap (ja_dummy);
      a_lt.swap (a_dummy);
   }

// Combine rows LU

   CIlu2_impl<_Int,_Flt>::CombineRowsLU (_n,
                                             ia_l, ja_l, a_l,
                                             ia_u, ja_u, a_u,
                                             _ia_lu, _ida_lu, _ja_lu, _a_lu);

};

//
// Perform ILU2 point factorization of the block with future diagonal modification
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::Ilu2Block (int _sctype, int _nitersc, int _fcttype, double _pivmin, double _tau1, double _tau2, double _theta,
                                             int _n, 
                                             vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                                             vector<_Int> &_ia_l, vector<_Int> &_ja_l, vector<_Flt> &_a_l,
                                             vector<_Int> &_ia_u, vector<_Int> &_ja_u, vector<_Flt> &_a_u,
                                             double &_eigmin_att, double &_eigmax_att) 
{

// Compute explicit scaling

   vector<_Flt> scl_L;
   vector<_Flt> scl_U;

   CIlu2_impl<_Int,_Flt>::ComputeScaling (_sctype, _nitersc,
                                                _n, 
                                                _ia_alu, _ja_alu, _a_alu,
                                                scl_L, scl_U);

// Perform explicit scaling

   int nza_alu = (int)_ia_alu[_n];

   vector<_Flt> a_scl_expl (nza_alu+1);

   _Flt *p_a_alu = &_a_alu[0];
   _Flt *pa_scl_expl = &a_scl_expl[0];

   memcpy (pa_scl_expl, p_a_alu, (size_t)(nza_alu*sizeof(_Flt)));

   CIlu2_impl<_Int,_Flt>::MatrixScale (_n, scl_L, scl_U,
                                             _ia_alu, _ja_alu, a_scl_expl);

// Split

   vector<_Int> ia_l;
   vector<_Int> ja_l;
   vector<_Flt> a_l;

   vector<_Int> ia_u;
   vector<_Int> ja_u;
   vector<_Flt> a_u;

   CIlu2_impl<_Int,_Flt>::SplitLU (_n,
                                       _ia_alu, _ja_alu, a_scl_expl,
                                       ia_l, ja_l, a_l,
                                       ia_u, ja_u, a_u);

   {
      vector<_Flt> a_dummy;
      a_scl_expl.swap (a_dummy);
   }

// Transpose L

   vector<_Int> ia_lt;
   vector<_Int> ja_lt;
   vector<_Flt> a_lt;

   CIlu2_impl<_Int,_Flt>::Transpose (_n,
                                          ia_l, ja_l, a_l,
                                          ia_lt, ja_lt, a_lt);

   {
      vector<_Int> ia_dummy;
      vector<_Int> ja_dummy;
      vector<_Flt> a_dummy;
      ia_l.swap (ia_dummy);
      ja_l.swap (ja_dummy);
      a_l.swap (a_dummy);
   }

// Combine into extended pairs

   vector<_Int> ia_alu_cnd;
   vector<_Int> ja_alu_cnd;
   vector<_Flt> a_alu_cnd;

   CIlu2_impl<_Int,_Flt>::CombinePairs (_n,
                                             ia_lt, ja_lt, a_lt,
                                             ia_u, ja_u, a_u,
                                             ia_alu_cnd, ja_alu_cnd, a_alu_cnd);

   {
      vector<_Int> ia_dummy;
      vector<_Int> ja_dummy;
      vector<_Flt> a_dummy;
      ia_lt.swap (ia_dummy);
      ja_lt.swap (ja_dummy);
      a_lt.swap (a_dummy);
      vector<_Int> ia_dummy1;
      vector<_Int> ja_dummy1;
      vector<_Flt> a_dummy1;
      ia_u.swap (ia_dummy1);
      ja_u.swap (ja_dummy1);
      a_u.swap (a_dummy1);
   }

// Perform fct

   vector<char> ja_char_alu_cnd;

   if (_fcttype >= 0) {
      int nzja_lu_cnd = (int)ia_alu_cnd[_n];
      ja_char_alu_cnd.resize (nzja_lu_cnd+1);
      int i;
      for (i=0;i<nzja_lu_cnd;i++) ja_char_alu_cnd[i] = 0;
   }

   vector<_Int> ia_lu_cnd;
   vector<_Int> ja_lu_cnd;
   vector<char> ja_char_lu_cnd;
   vector<_Flt> a_lu_cnd;

   CIlu2_impl<_Int,_Flt>::Ilu2BlockIlu2 (_fcttype, _pivmin, _tau1, _tau2, _theta, 
                                             _n, _n,
                                             ia_alu_cnd, ja_alu_cnd, ja_char_alu_cnd, a_alu_cnd,
                                             ia_lu_cnd, ja_lu_cnd, ja_char_lu_cnd, a_lu_cnd,
                                             _eigmin_att, _eigmax_att);

   {
      vector<_Int> ia_dummy;
      vector<_Int> ja_dummy;
      vector<char> ja_char_dummy;
      vector<_Flt> a_dummy;
      vector<char> ja_char1_dummy;
      ia_alu_cnd.swap (ia_dummy);
      ja_alu_cnd.swap (ja_dummy);
      ja_char_alu_cnd.swap (ja_char_dummy);
      a_alu_cnd.swap (a_dummy);
      ja_char_lu_cnd.swap (ja_char1_dummy);
   }

// Split pairs with filtering

   CIlu2_impl<_Int,_Flt>::SplitPairsFilter (_tau1, _n,
                                                ia_lu_cnd, ja_lu_cnd, a_lu_cnd,
                                                _ia_l, _ja_l, _a_l,
                                                _ia_u, _ja_u, _a_u);

   {
      vector<_Int> ia_dummy;
      vector<_Int> ja_dummy;
      vector<_Flt> a_dummy;
      ia_lu_cnd.swap (ia_dummy);
      ja_lu_cnd.swap (ja_dummy);
      a_lu_cnd.swap (a_dummy);
   }

// Perform explicit rescaling of compute triangular factors

   vector<_Flt> inv_scl_L (_n+1);
   vector<_Flt> inv_scl_U (_n+1);

   CIlu2_impl<_Int,_Flt>::InverseDiag (_n, scl_L, inv_scl_L);
   CIlu2_impl<_Int,_Flt>::InverseDiag (_n, scl_U, inv_scl_U);

   CIlu2_impl<_Int,_Flt>::RescaleU (_n, scl_L, inv_scl_L,
                                          _ia_l, _ja_l, _a_l);
   CIlu2_impl<_Int,_Flt>::RescaleU (_n, scl_U, inv_scl_U,
                                          _ia_u, _ja_u, _a_u);

   {
      vector<_Flt> dl_dummy;
      vector<_Flt> du_dummy;
      vector<_Flt> idl_dummy;
      vector<_Flt> idu_dummy;
      scl_L.swap (dl_dummy);
      scl_U.swap (du_dummy);
      inv_scl_L.swap (idl_dummy);
      inv_scl_U.swap (idu_dummy);
   }

};

//
// Compute scaling
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::ComputeScaling (int _sctype, int _nitersc,
                                                int _n, 
                                                vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                                                vector<_Flt> &_sclL, vector<_Flt> &_sclU)
{

   _sclL.resize (_n+1);
   _sclU.resize (_n+1);

// Simple diagonal based scaling

   int i, j, jj;

   _Flt diag;

   if (_sctype == -1) 
   {
      diag = (_Flt)1.;
      for (i=0;i<_n;i++) 
      {
         _sclL[i] = diag;
         _sclU[i] = diag;
      }
   } else if (_sctype == 0) 
   {
      for (i=0;i<_n;i++) 
      {
         diag = 1.;
         for (j=(int)_ia_alu[i];j<_ia_alu[i+1];j++) 
         {
            jj = (int)_ja_alu[j];
            if (i == jj) 
            {
               diag = _a_alu[j];
            }
         }
         if (diag > 0.) 
         {
            diag = (_Flt)(1./sqrt(diag));
            _sclL[i] = diag;
            _sclU[i] = diag;
         } else {
            diag = -diag;
            diag = (_Flt)(1./sqrt(diag));
            _sclL[i] = -diag;
            _sclU[i] = diag;
         }
      }
   } else {

// Iterative rows/columns balancing scaling

      _Flt aux, aux1;
      int iter;

      for (i=0;i<_n;i++) 
      {
         aux = 0.;
         for (j=(int)_ia_alu[i];j<_ia_alu[i+1];j++) 
         {
            jj = (int)_ja_alu[j];
            aux1 = _a_alu[j];
            aux += aux1*aux1;
         }
         aux = (_Flt)(1.0 / sqrt(aux));
         _sclL[i] = aux;
      }
      for (iter=0;iter<_nitersc;iter++) 
      {
         for (i=0;i<_n;i++) _sclU[i] = 0.;
         for (i=0;i<_n;i++) 
         {
            for (j=(int)_ia_alu[i];j<_ia_alu[i+1];j++) 
            {
               jj = (int)_ja_alu[j];
               aux1 = _a_alu[j];
               _sclU[jj] += _sclL[i]*aux1*aux1;
            }
         }
         for (i=0;i<_n;i++) 
         {
            _sclU[i] = (_Flt)(1./_sclU[i]);
         }
         for (i=0;i<_n;i++) 
         {
            aux = 0.;
            for (j=(int)_ia_alu[i];j<_ia_alu[i+1];j++) 
            {
               jj = (int)_ja_alu[j];
               aux1 = _a_alu[j];
               aux += _sclU[jj]*aux1*aux1;
            }
            aux = (_Flt)(1.0 / aux);
            _sclL[i] = aux;
         }
      }
      for (i=0;i<_n;i++) _sclL[i] = (_Flt)(sqrt(_sclL[i]));
      for (i=0;i<_n;i++) _sclU[i] = (_Flt)(sqrt(_sclU[i]));
   }

};

//
// Perform explicit scaling
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::MatrixScale (int _n, vector<_Flt> &_sclL, vector<_Flt> &_sclU,
                                                vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu)
{

   int i, j, jj;

   for (i=0;i<_n;i++) 
   {
      for (j=(int)_ia_alu[i];j<_ia_alu[i+1];j++) 
      {
         jj = (int)_ja_alu[j];
         _a_alu[j] = _sclL[i]*_a_alu[j]*_sclU[jj];
      }
   }

}

//
// Symmetrize sparsity
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::SymmetrizeSparsity (int _n, 
                                                      const vector<_Int> &_ia_alu, const vector<_Int> &_ja_alu,
                                                      vector<_Int> &_ia_symm, vector<_Int> &_ja_symm)
{

// Split

   vector<_Int> ia_l;
   vector<_Int> ja_l;

   vector<_Int> ia_u;
   vector<_Int> ja_u;

   CIlu2_impl<_Int,_Flt>::SplitLUSp (_n,
                                       _ia_alu, _ja_alu,
                                       ia_l, ja_l,
                                       ia_u, ja_u);

// Transpose L

   vector<_Int> ia_lt;
   vector<_Int> ja_lt;

   CIlu2_impl<_Int,_Flt>::TransposeSp (_n,
                                             ia_l, ja_l,
                                             ia_lt, ja_lt);

   {
      vector<_Int> ia_dummy;
      vector<_Int> ja_dummy;
      ia_l.swap (ia_dummy);
      ja_l.swap (ja_dummy);
   }

// Combine into extended pairs

   vector<_Int> ia_alu_cnd;
   vector<_Int> ja_alu_cnd;

   CIlu2_impl<_Int,_Flt>::CombineLUSp (_n,
                                             ia_lt, ja_lt,
                                             ia_u, ja_u,
                                             ia_alu_cnd, ja_alu_cnd);

   {
      vector<_Int> ia_dummy;
      vector<_Int> ja_dummy;
      ia_lt.swap (ia_dummy);
      ja_lt.swap (ja_dummy);
      vector<_Int> ia_dummy1;
      vector<_Int> ja_dummy1;
      ia_u.swap (ia_dummy1);
      ja_u.swap (ja_dummy1);
   }

// Transpose combined sparsity

   vector<_Int> ia_alu_cnd_t;
   vector<_Int> ja_alu_cnd_t;

   CIlu2_impl<_Int,_Flt>::TransposeSp (_n,
                                             ia_alu_cnd, ja_alu_cnd,
                                             ia_alu_cnd_t, ja_alu_cnd_t);

// Compute symmetrized matrix

   CIlu2_impl<_Int,_Flt>::CombineLUSp (_n,
                                             ia_alu_cnd_t, ja_alu_cnd_t,
                                             ia_alu_cnd, ja_alu_cnd,
                                             _ia_symm, _ja_symm);

};

//
// Split matrix data into L and U parts for sparsity only
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::SplitLUSp (int _n, 
                                          const vector<_Int> &_ia_alu, const vector<_Int> &_ja_alu,
                                          vector<_Int> &_ia_l, vector<_Int> &_ja_l,
                                          vector<_Int> &_ia_u, vector<_Int> &_ja_u)
{

   int nzja_alu = (int)_ia_alu[_n];

// Count number of elems

   int nzja_l = 0;
   int nzja_u = 0;

   int i, j, jj;

   for (i=0;i<_n;i++) 
   {
      for (j=(int)_ia_alu[i];j<_ia_alu[i+1];j++) 
      {
         jj = (int)_ja_alu[j];
         if (jj <= i) nzja_l++;
         if (jj >= i) nzja_u++;
      }
   }

// Allocate and fill

   _ia_l.resize (_n+1);
   _ja_l.resize (nzja_l+1);
   _ia_u.resize (_n+1);
   _ja_u.resize (nzja_u+1);

   nzja_l = 0;
   nzja_u = 0;

   _ia_l[0] = 0;
   _ia_u[0] = 0;

   for (i=0;i<_n;i++) 
   {
      for (j=(int)_ia_alu[i];j<_ia_alu[i+1];j++) 
      {
         jj = (int)_ja_alu[j];
         if (jj <= i) 
         {
            _ja_l[nzja_l] = jj;
            nzja_l++;
         }
         if (jj >= i) 
         {
            _ja_u[nzja_u] = jj;
            nzja_u++;
         }
      }
      _ia_l[i+1] = nzja_l;
      _ia_u[i+1] = nzja_u;
   }

};

//
// Split matrix data into L and U parts
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::SplitLU (int _n, 
                                          vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                                          vector<_Int> &_ia_l, vector<_Int> &_ja_l, vector<_Flt> &_a_l,
                                          vector<_Int> &_ia_u, vector<_Int> &_ja_u, vector<_Flt> &_a_u)
{

   int nzja_alu = (int)_ia_alu[_n];

// Count number of elems

   int nzja_l = 0;
   int nzja_u = 0;

   int i, j, jj;

   for (i=0;i<_n;i++) 
   {
      for (j=(int)_ia_alu[i];j<_ia_alu[i+1];j++) 
      {
         jj = (int)_ja_alu[j];
         if (jj <= i) nzja_l++;
         if (jj >= i) nzja_u++;
      }
   }

// Allocate and fill

   _ia_l.resize (_n+1);
   _ja_l.resize (nzja_l+1);
   _a_l.resize (nzja_l+1);
   _ia_u.resize (_n+1);
   _ja_u.resize (nzja_u+1);
   _a_u.resize (nzja_u+1);

   nzja_l = 0;
   nzja_u = 0;

   _ia_l[0] = 0;
   _ia_u[0] = 0;

   for (i=0;i<_n;i++) 
   {
      for (j=(int)_ia_alu[i];j<_ia_alu[i+1];j++) 
      {
         jj = (int)_ja_alu[j];
         if (jj <= i) 
         {
            _ja_l[nzja_l] = jj;
            _a_l[nzja_l] = _a_alu[j];
            nzja_l++;
         }
         if (jj >= i) 
         {
            _ja_u[nzja_u] = jj;
            _a_u[nzja_u] = _a_alu[j];
            nzja_u++;
         }
      }
      _ia_l[i+1] = nzja_l;
      _ia_u[i+1] = nzja_u;
   }

};

//
// Transpose square matrix for sparsity only
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::TransposeSp (int _n, 
                                             vector<_Int> &_ia_a, vector<_Int> &_ja_a,
                                             vector<_Int> &_ia_at, vector<_Int> &_ja_at)
{

   int nzja_a = (int)_ia_a[_n];

// Allocate work array

   vector<_Int> iptr (_n+1);

// Allocate tranposed data

   _ia_at.resize (_n+1);
   _ja_at.resize (nzja_a+1);

   int i, j, jj, k;

   for (i=0;i<=_n;i++) _ia_at[i] = 0;

   for (i=0;i<nzja_a;i++) 
   {
      jj = (int)_ja_a[i];
      _ia_at[jj+1]++;
   }

   for (i=0;i<_n;i++) _ia_at[i+1] = _ia_at[i]+_ia_at[i+1];

   for (i=0;i<_n;i++) iptr[i] = _ia_at[i];

   for (i=0;i<_n;i++) 
   {
      for (j=(int)_ia_a[i];j<_ia_a[i+1];j++) 
      {
         jj = (int)_ja_a[j];
         k = (int)iptr[jj];
         _ja_at[k] = i;
         iptr[jj]++;
      }
   }

};

//
// Transpose square matrix
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::Transpose (int _n, 
                                             vector<_Int> &_ia_a, vector<_Int> &_ja_a, vector<_Flt> &_a_a,
                                             vector<_Int> &_ia_at, vector<_Int> &_ja_at, vector<_Flt> &_a_at)
{

   int nzja_a = (int)_ia_a[_n];

// Allocate work array

   vector<_Int> iptr (_n+1);

// Allocate tranposed data

   _ia_at.resize (_n+1);
   _ja_at.resize (nzja_a+1);
   _a_at.resize (nzja_a+1);

   int i, j, jj, k;

   for (i=0;i<=_n;i++) _ia_at[i] = 0;

   for (i=0;i<nzja_a;i++) 
   {
      jj = (int)_ja_a[i];
      _ia_at[jj+1]++;
   }

   for (i=0;i<_n;i++) _ia_at[i+1] = _ia_at[i]+_ia_at[i+1];

   for (i=0;i<_n;i++) iptr[i] = _ia_at[i];

   for (i=0;i<_n;i++) 
   {
      for (j=(int)_ia_a[i];j<_ia_a[i+1];j++) 
      {
         jj = (int)_ja_a[j];
         k = (int)iptr[jj];
         _ja_at[k] = i;
         _a_at[k] = _a_a[j];
         iptr[jj]++;
      }
   }

};

//
// Combine sprsity of L and U
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::CombineLUSp (int _n, 
                                                vector<_Int> &_ia_l, vector<_Int> &_ja_l,
                                                vector<_Int> &_ia_u, vector<_Int> &_ja_u,
                                                vector<_Int> &_ia_alu, vector<_Int> &_ja_alu)
{

// Compute number of extended elems

   int nzja_ext = 0;

   int i, ipl, ipu, iendl, iendu, jj_l, jj_u;

   for (i=0;i<_n;i++) 
   {
      ipl = (int)_ia_l[i];
      ipu = (int)_ia_u[i];
      iendl = (int)_ia_l[i+1]-1;
      iendu = (int)_ia_u[i+1]-1;
      while (ipl <= iendl || ipu <= iendu) 
      {
         if (ipl <= iendl && ipu <= iendu) 
         {
            jj_l = (int)_ja_l[ipl];
            jj_u = (int)_ja_u[ipu];
            if (jj_l == jj_u) {
               nzja_ext++;
               ipl++;
               ipu++;
            }
            else if (jj_l < jj_u)
            {
               nzja_ext++;
               ipl++;
            }
            else
            {
               nzja_ext++;
               ipu++;
            }
         }
         else if (ipl <= iendl) 
         {
            nzja_ext++;
            ipl++;
         } 
         else
         {
            nzja_ext++;
            ipu++;
         }
      }
   }

// Count number of elems

   _ia_alu.resize (_n+1);
   _ja_alu.resize (nzja_ext+1);

   _ia_alu[0] = 0;
   nzja_ext = 0;

   for (i=0;i<_n;i++) 
   {
      ipl = (int)_ia_l[i];
      ipu = (int)_ia_u[i];
      iendl = (int)_ia_l[i+1]-1;
      iendu = (int)_ia_u[i+1]-1;
      while (ipl <= iendl || ipu <= iendu) 
      {
         if (ipl <= iendl && ipu <= iendu) 
         {
            jj_l = (int)_ja_l[ipl];
            jj_u = (int)_ja_u[ipu];
            if (jj_l == jj_u) {
               _ja_alu[nzja_ext] = jj_l;
               nzja_ext++;
               ipl++;
               ipu++;
            }
            else if (jj_l < jj_u)
            {
               _ja_alu[nzja_ext] = jj_l;
               nzja_ext++;
               ipl++;
            }
            else
            {
               _ja_alu[nzja_ext] = jj_u;
               nzja_ext++;
               ipu++;
            }
         }
         else if (ipl <= iendl) 
         {
            _ja_alu[nzja_ext] = _ja_l[ipl];
            nzja_ext++;
            ipl++;
         } 
         else
         {
            _ja_alu[nzja_ext] = _ja_u[ipu];
            nzja_ext++;
            ipu++;
         }
      }
      _ia_alu[i+1] = nzja_ext;
   }

};

//
// Combine L and U data into extended pairs
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::CombinePairs (int _n, 
                                                vector<_Int> &_ia_l, vector<_Int> &_ja_l, vector<_Flt> &_a_l,
                                                vector<_Int> &_ia_u, vector<_Int> &_ja_u, vector<_Flt> &_a_u,
                                                vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu)
{

// Compute number of extended elems

   int nzja_ext = 0;

   int i, ipl, ipu, iendl, iendu, jj_l, jj_u;

   for (i=0;i<_n;i++) 
   {
      ipl = (int)_ia_l[i];
      ipu = (int)_ia_u[i];
      iendl = (int)_ia_l[i+1]-1;
      iendu = (int)_ia_u[i+1]-1;
      while (ipl <= iendl || ipu <= iendu) 
      {
         if (ipl <= iendl && ipu <= iendu) 
         {
            jj_l = (int)_ja_l[ipl];
            jj_u = (int)_ja_u[ipu];
            if (jj_l == jj_u) {
               nzja_ext++;
               ipl++;
               ipu++;
            }
            else if (jj_l < jj_u)
            {
               nzja_ext++;
               ipl++;
            }
            else
            {
               nzja_ext++;
               ipu++;
            }
         }
         else if (ipl <= iendl) 
         {
            nzja_ext++;
            ipl++;
         } 
         else
         {
            nzja_ext++;
            ipu++;
         }
      }
   }

// Count number of elems

   _ia_alu.resize (_n+1);
   _ja_alu.resize (nzja_ext+1);
   _a_alu.resize (nzja_ext*2+1);

   _ia_alu[0] = 0;
   nzja_ext = 0;

   for (i=0;i<_n;i++) 
   {
      ipl = (int)_ia_l[i];
      ipu = (int)_ia_u[i];
      iendl = (int)_ia_l[i+1]-1;
      iendu = (int)_ia_u[i+1]-1;
      while (ipl <= iendl || ipu <= iendu) 
      {
         if (ipl <= iendl && ipu <= iendu) 
         {
            jj_l = (int)_ja_l[ipl];
            jj_u = (int)_ja_u[ipu];
            if (jj_l == jj_u) {
               _ja_alu[nzja_ext] = jj_l;
               _a_alu[nzja_ext*2] = _a_l[ipl];
               _a_alu[nzja_ext*2+1] = _a_u[ipu];
               nzja_ext++;
               ipl++;
               ipu++;
            }
            else if (jj_l < jj_u)
            {
               _ja_alu[nzja_ext] = jj_l;
               _a_alu[nzja_ext*2] = _a_l[ipl];
               _a_alu[nzja_ext*2+1] = 0.;
               nzja_ext++;
               ipl++;
            }
            else
            {
               _ja_alu[nzja_ext] = jj_u;
               _a_alu[nzja_ext*2] = 0.;
               _a_alu[nzja_ext*2+1] = _a_u[ipu];
               nzja_ext++;
               ipu++;
            }
         }
         else if (ipl <= iendl) 
         {
            _ja_alu[nzja_ext] = _ja_l[ipl];
            _a_alu[nzja_ext*2] = _a_l[ipl];
            _a_alu[nzja_ext*2+1] = 0.;
            nzja_ext++;
            ipl++;
         } 
         else
         {
            _ja_alu[nzja_ext] = _ja_u[ipu];
            _a_alu[nzja_ext*2] = 0.;
            _a_alu[nzja_ext*2+1] = _a_u[ipu];
            nzja_ext++;
            ipu++;
         }
      }
      _ia_alu[i+1] = nzja_ext;
   }

};

//
// Split pairs fct data into L and U parts with post filtering
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::SplitPairsFilter (double _tau1, int _n, 
                                                   vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                                                   vector<_Int> &_ia_l, vector<_Int> &_ja_l, vector<_Flt> &_a_l,
                                                   vector<_Int> &_ia_u, vector<_Int> &_ja_u, vector<_Flt> &_a_u)
{

   int nzja_alu = (int)_ia_alu[_n];

// Allocate and fill

   _ia_l.resize (_n+1);
   _ja_l.resize (nzja_alu+1);
   _a_l.resize (nzja_alu+1);
   _ia_u.resize (_n+1);
   _ja_u.resize (nzja_alu+1);
   _a_u.resize (nzja_alu+1);

   int nzja_l = 0;
   int nzja_u = 0;

   int i, j, jj;
   double auxL, auxU;

   _ia_l[0] = 0;
   _ia_u[0] = 0;

   for (i=0;i<_n;i++) 
   {
      for (j=(int)_ia_alu[i];j<_ia_alu[i+1];j++) 
      {
         jj = (int)_ja_alu[j];
         auxL = _a_alu[j*2];
         auxU = _a_alu[j*2+1];
         if (auxL < 0.) auxL = -auxL;
         if (auxU < 0.) auxU = -auxU;
         if (jj == i || auxL >= _tau1) 
         {
            _ja_l[nzja_l] = jj;
            _a_l[nzja_l] = _a_alu[j*2];
            nzja_l++;
         }
         if (jj == i || auxU >= _tau1) 
         {
            _ja_u[nzja_u] = jj;
            _a_u[nzja_u] = _a_alu[j*2+1];
            nzja_u++;
         }
      }
      _ia_l[i+1] = nzja_l;
      _ia_u[i+1] = nzja_u;
   }

};

//
// Combine L and U data into extended pairs
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::CombineRowsLU (int _n, 
                                                vector<_Int> &_ia_l, vector<_Int> &_ja_l, vector<_Flt> &_a_l,
                                                vector<_Int> &_ia_u, vector<_Int> &_ja_u, vector<_Flt> &_a_u,
                                                vector<_Int> &_ia_alu, vector<_Int> &_ida_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu)
{

// Compute number of extended elems

   int nzja_ext = (int)(_ia_l[_n]+_ia_u[_n]);

// Count number of elems

   _ia_alu.resize (_n+1);
   _ida_alu.resize (_n+1);
   _ja_alu.resize (nzja_ext+1);
   _a_alu.resize (nzja_ext+1);

   _ia_alu[0] = 0;

   int i, j;

   nzja_ext = 0;

   for (i=0;i<_n;i++) 
   {
      for (j=(int)_ia_l[i];j<_ia_l[i+1];j++) 
      {
         _ja_alu[nzja_ext] = _ja_l[j];
         _a_alu[nzja_ext] = _a_l[j];
         nzja_ext++;
      }
      _ida_alu[i] = nzja_ext;
      for (j=(int)_ia_u[i];j<_ia_u[i+1];j++) 
      {
         _ja_alu[nzja_ext] = _ja_u[j];
         _a_alu[nzja_ext] = _a_u[j];
         nzja_ext++;
      }
      _ia_alu[i+1] = nzja_ext;
   }

};

//
// Compute inverse scaling
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::InverseDiag (int _n, vector<_Flt> &_sclU, vector<_Flt> &_invsclU)
{

   int i;

   for (i=0;i<_n;i++) 
   {
      _invsclU[i] = (_Flt)(1. / _sclU[i]);
   }

}

//
// Rescale factor back
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::RescaleU (int _n, vector<_Flt> &_sclU, vector<_Flt> &_invsclU,
                                            vector<_Int> &_ia_u, vector<_Int> &_ja_u, vector<_Flt> &_a_u)
{

   int i, j, jj;

   for (i=0;i<_n;i++) 
   {
      for (j=(int)_ia_u[i];j<_ia_u[i+1];j++) 
      {
         jj = (int)_ja_u[j];
         if (jj != i) {
            _a_u[j] *= _invsclU[jj];
         } else {
            _a_u[j] *= _sclU[i];
         }
      }
   }

};

//
// Perform ILU2 point factorization of the block with future diagonal modification (no structural control)
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::Ilu2BlockIlu2 (double _pivmin, double _tau1, double _tau2, double _theta,
                                                int _n, int _n_ini, 
                                                vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<_Flt> &_a_alu,
                                                vector<_Int> &_ia_lu, vector<_Int> &_ja_lu, vector<_Flt> &_a_lu,
                                                double &_eigmin_att, double &_eigmax_att) 
{

// Open structures

   int nzja_alu = (int)_ia_alu[_n];

// Prepare mask arrays

   vector<_Int> ibegm  (_n+1);
   vector<_Int> madj   (_n+1);
   vector<_Int> iv     (_n+1);
   vector<_Int> imaskc (_n+1);
   vector<_Int> listc  (_n+1);

   vector<_Flt> fmaskc (2*_n+1);
   vector<_Flt> a_dia (_n+1);

   _Int *plistc = &listc[0];

   int i;

   for (i=0;i<=_n;i++) ibegm[i] = -1;
   for (i=0;i<=_n;i++) madj[i] = -1;
   for (i=0;i<=_n;i++) imaskc[i] = -1;
   for (i=0;i<_n;i++) a_dia[i] = 0.;

   double eigmin_att = FLT_MAX;
   double eigmax_att = -FLT_MAX;

   _Flt fzero = 0.;
   _Flt fone = 1.;

   int j, jj, k, kcolmn, irwprv, nlistcnew, ind, ind1, irow, irwpr1, jcolmn;
   double auxL, auxU;
   double aux1, aux2, dfnorm;

   int icycle_int = -1;

   _ia_lu.resize(_n+1);

   _ja_lu.reserve(_n+1);
   _a_lu.reserve(_n*2+1);

   _ia_lu[0] = 0;

   int nzja_lu = 0;

   for (i=0;i<_n;i++) 
   {

      irow = i;

// Init current row

      int nlistcloc = 0;

      icycle_int++;

      if (i < _n_ini) {

         for (k=(int)_ia_alu[i];k<_ia_alu[i+1];k++) 
         {
            kcolmn = (int)_ja_alu[k];
            listc[nlistcloc] = kcolmn;
            nlistcloc++;
            fmaskc[kcolmn*2] = _a_alu[k*2];
            fmaskc[kcolmn*2+1] = _a_alu[k*2+1];
            imaskc[kcolmn] = icycle_int;
         }

      } else {
         listc[nlistcloc] = i;
         nlistcloc++;
         fmaskc[i*2] = fzero;
         fmaskc[i*2+1] = fzero;
         imaskc[i] = icycle_int;
      }

// Update current row

      irwprv = (int)ibegm[irow];

      while (irwprv != -1) 
      {
         irwpr1 = (int)madj[irwprv];
         if (iv[irwprv] < _ia_lu[irwprv+1]) 
         {

            j = (int)iv[irwprv];

            jj = (int)_ja_lu[j];

            auxL = _a_lu[j*2];
            auxU = _a_lu[j*2+1];

            jcolmn = jj;
            if (jcolmn < 0) jcolmn = -jcolmn-1;

            if (jj >= 0) 
            {
               for (k=(int)iv[irwprv];k<_ia_lu[irwprv+1];k++) 
               {
                  kcolmn = (int)_ja_lu[k];
                  if (kcolmn < 0) kcolmn = -kcolmn-1;
                  if (imaskc[kcolmn] != icycle_int) 
                  {
                     listc[nlistcloc] = kcolmn;
                     nlistcloc++;
                     fmaskc[kcolmn*2] = fzero;
                     fmaskc[kcolmn*2+1] = fzero;
                     imaskc[kcolmn] = icycle_int;
                  }
                  fmaskc[kcolmn*2] -= (_Flt)(auxU * _a_lu[k*2]);
                  fmaskc[kcolmn*2+1] -= (_Flt)(auxL * _a_lu[k*2+1]);
               }
            }
            else 
            {
               for (k=(int)iv[irwprv];k<_ia_lu[irwprv+1];k++) 
               {
                  kcolmn = (int)_ja_lu[k];
                  if (kcolmn >= 0) 
                  {
                     if (imaskc[kcolmn] != icycle_int) 
                     {
                        listc[nlistcloc] = kcolmn;
                        nlistcloc++;
                        fmaskc[kcolmn*2] = fzero;
                        fmaskc[kcolmn*2+1] = fzero;
                        imaskc[kcolmn] = icycle_int;
                     }
                     fmaskc[kcolmn*2] -= (_Flt)(auxU * _a_lu[k*2]);
                     fmaskc[kcolmn*2+1] -= (_Flt)(auxL * _a_lu[k*2+1]);
                  }
               }
            }

         }

         irwprv = irwpr1;

      }

// Update transposed structure

      irwprv = (int)ibegm[irow];

      while (irwprv != -1) 
      {

         irwpr1 = (int)madj[irwprv];

         if (iv[irwprv] >= _ia_lu[irwprv+1]-1) 
         {
            madj[irwprv] = -1;
         } 
         else 
         {
            j = (int)iv[irwprv]+1;
            jcolmn = (int)_ja_lu[j];
            if (jcolmn < 0) jcolmn = -jcolmn-1;
            madj[irwprv] = ibegm[jcolmn];
            ibegm[jcolmn] = irwprv;
         }

         iv[irwprv]++;

         irwprv = irwpr1;

      }

// Perform filtering of the data and modify diagonal

      if (i<_n_ini) 
      {

         nlistcnew = 0;

         for (j=0;j<nlistcloc;j++) 
         {
            ind = (int)listc[j];
            if (ind != irow) 
            {
               aux1 = fmaskc[ind*2];
               aux2 = fmaskc[ind*2+1];
               if (aux1 < fzero) aux1 = -aux1;
               if (aux2 < fzero) aux2 = -aux2;
               dfnorm = (aux1 > aux2) ? aux1 : aux2;
               if (dfnorm < _tau2) {
                  a_dia[irow] += (_Flt)(dfnorm);
                  a_dia[ind] += (_Flt)(dfnorm);
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

      if (fmaskc[irow*2] > fzero) 
      {
         fmaskc[irow*2] += (_Flt)(a_dia[irow]*_theta);
      } else {
         fmaskc[irow*2] -= (_Flt)(a_dia[irow]*_theta);
      }

// Factorize current row and split into first and second order

      sort (plistc, plistc+nlistcloc);

      if (i<_n_ini) 
      {

         aux1 = fmaskc[irow*2];

         if (aux1 > 0.) {
            if (aux1 < eigmin_att) eigmin_att = aux1;
            if (aux1 > eigmax_att) eigmax_att = aux1;
            if (aux1 < _pivmin) aux1 = _pivmin;
            aux1 = sqrt(aux1);
            aux1 = fone / aux1;
            auxL = aux1;
            auxU = aux1;
         } else {
            aux1 = -aux1;
            if (aux1 < eigmin_att) eigmin_att = aux1;
            if (aux1 > eigmax_att) eigmax_att = aux1;
            if (aux1 < _pivmin) aux1 = _pivmin;
            aux1 = sqrt(aux1);
            aux1 = fone / aux1;
            auxL = -aux1;
            auxU = aux1;
         }
         fmaskc[irow*2] = (_Flt)(auxL);
         fmaskc[irow*2+1] = (_Flt)(auxU);
      }

      if (i<_n_ini) 
      {
         for (j=0;j<nlistcloc;j++) {
            ind = (int)listc[j];
            if (ind != irow) {
               fmaskc[ind*2] *= (_Flt)(auxU);
               fmaskc[ind*2+1] *= (_Flt)(auxL);
               aux1 = fmaskc[ind*2];
               aux2 = fmaskc[ind*2+1];
               if (aux1 < fzero) aux1 = -aux1;
               if (aux2 < fzero) aux2 = -aux2;
               dfnorm = (aux1 > aux2) ? aux1 : aux2;
               if (dfnorm < _tau1) {
                  listc[j] = (_Int) (-ind-1);
               }
            }
         }
      }

// Store computed row elems

      for (j=0;j<nlistcloc;j++) {
         ind = (int)plistc[j];
         _ja_lu.push_back (ind);
         ind1 = ind;
         if (ind1 < 0) ind1 = -ind1-1;
         _a_lu.push_back (fmaskc[ind1*2]);
         _a_lu.push_back (fmaskc[ind1*2+1]);
      };

      nzja_lu += nlistcloc;

      _ia_lu[i+1] = nzja_lu;

// Add current row into the transposed structures

      if (nlistcloc > 1) 
      {
         if (i < _n_ini) {
            iv[irow] = _ia_lu[i]+1;
            ind = (int)listc[1];
            if (ind < 0) ind = -ind-1;
            madj[irow] = ibegm[ind];
            ibegm[ind] = (_Int) irow;
         } else {
            iv[irow] = _ia_lu[i+1]+1;
         }
      }

   }

   _eigmin_att = eigmin_att;
   _eigmax_att = eigmax_att;

// Condense the result (filter second order elems)

   int nzja_new = 0;

   for (i=0;i<nzja_lu;i++) 
   {
      if (_ja_lu[i] >= 0) nzja_new++;
   }

   vector<_Int> jlu_cnd (nzja_new+1);
   vector<_Flt> lu_cnd (2*nzja_new+1);

   nzja_new = 0;

   listc[0] = 0;

   for (i=0;i<_n;i++) 
   {
      for (j=(int)_ia_lu[i];j<_ia_lu[i+1];j++) 
      {
         if (_ja_lu[j] >= 0) 
         {
            jlu_cnd[nzja_new] = _ja_lu[j];
            lu_cnd[nzja_new*2] = _a_lu[j*2];
            lu_cnd[nzja_new*2+1] = _a_lu[j*2+1];
            nzja_new++;
         }
      }
      listc[i+1] = nzja_new;
   }

   for (i=0;i<=_n;i++) _ia_lu[i] = listc[i];

   _ja_lu.swap (jlu_cnd);
   _a_lu.swap (lu_cnd);

};

//
// Perform ILU2 point factorization of the block with future diagonal modification (with structural control)
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::Ilu2BlockIlu2 (int _fcttype, double _pivmin, double _tau1, double _tau2, double _theta,
                                                int _n, int _n_ini, 
                                                vector<_Int> &_ia_alu, vector<_Int> &_ja_alu, vector<char> &_ja_char_alu, vector<_Flt> &_a_alu,
                                                vector<_Int> &_ia_lu, vector<_Int> &_ja_lu, vector<char> &_ja_char_lu, vector<_Flt> &_a_lu,
                                                double &_eigmin_att, double &_eigmax_att) 
{

// Fast return

   if (_fcttype == -1) {
      CIlu2_impl<_Int,_Flt>::Ilu2BlockIlu2 (_pivmin, _tau1, _tau2, _theta,
                                                _n, _n_ini, 
                                                _ia_alu, _ja_alu, _a_alu,
                                                _ia_lu, _ja_lu, _a_lu,
                                                _eigmin_att, _eigmax_att);
      return;
   }

// Open structures

   int nzja_alu = (int)_ia_alu[_n];

// Prepare mask arrays

   vector<_Int> ibegm    (_n+1);
   vector<_Int> madj     (_n+1);
   vector<_Int> iv       (_n+1);
   vector<_Int> imaskc   (_n+1);
   vector<char> imaskchar (_n+1);
   vector<_Int> listc    (_n+1);

   vector<_Flt> fmaskc (2*_n+1);
   vector<_Flt> a_dia (_n+1);

   _Int *plistc = &listc[0];

   int i;

   for (i=0;i<=_n;i++) ibegm[i] = -1;
   for (i=0;i<=_n;i++) madj[i] = -1;
   for (i=0;i<=_n;i++) imaskc[i] = -1;
   for (i=0;i<_n;i++) a_dia[i] = 0.;

   double eigmin_att = FLT_MAX;
   double eigmax_att = -FLT_MAX;

   int j, jj, k, kcolmn, irwprv, nlistcnew, ind, ind1, irow, irwpr1, jcolmn;
   double auxL, auxU;
   double aux1, aux2, dfnorm;

   _Flt fzero = 0.;
   _Flt fone = 1.;

   int icycle_int = -1;

   _ia_lu.resize(_n+1);

   _ja_lu.reserve(_n+1);
   _ja_char_lu.reserve(_n+1);
   _a_lu.reserve(_n*2+1);

   _ia_lu[0] = 0;

   int nzja_lu = 0;

   char jjchar1, jjchar2, jjchar;

   for (i=0;i<_n;i++) 
   {

      irow = i;

// Init current row

      int nlistcloc = 0;

      icycle_int++;

      if (i<_n_ini) 
      {
         for (k=(int)_ia_alu[i];k<_ia_alu[i+1];k++) 
         {
            kcolmn = (int)_ja_alu[k];
            listc[nlistcloc] = kcolmn;
            nlistcloc++;
            fmaskc[kcolmn*2] = _a_alu[k*2];
            fmaskc[kcolmn*2+1] = _a_alu[k*2+1];
            imaskc[kcolmn] = icycle_int;
            imaskchar[kcolmn] = _ja_char_alu[k];
         }
      } else {
         listc[nlistcloc] = i;
         nlistcloc++;
         fmaskc[i*2] = fzero;
         fmaskc[i*2+1] = fzero;
         imaskc[i] = icycle_int;
         imaskchar[i] = 0;
      }

// Update current row

      irwprv = (int)ibegm[irow];

      while (irwprv != -1) 
      {
         irwpr1 = (int)madj[irwprv];
         if (iv[irwprv] < _ia_lu[irwprv+1]) 
         {

            j = (int)iv[irwprv];

            jj = (int)_ja_lu[j];
            jjchar1 = _ja_char_lu[j]+1;

            auxL = _a_lu[j*2];
            auxU = _a_lu[j*2+1];

            jcolmn = jj;
            if (jcolmn < 0) jcolmn = -jcolmn-1;

            if (jj >= 0) 
            {
               for (k=(int)iv[irwprv];k<_ia_lu[irwprv+1];k++) 
               {
                  kcolmn = (int)_ja_lu[k];
                  if (kcolmn < 0) kcolmn = -kcolmn-1;
                  jjchar2 = _ja_char_lu[k]+1;
                  jjchar = (jjchar1 < jjchar2) ? jjchar2 : jjchar1;
                  if (imaskc[kcolmn] != icycle_int) 
                  {
                     listc[nlistcloc] = kcolmn;
                     nlistcloc++;
                     fmaskc[kcolmn*2] = fzero;
                     fmaskc[kcolmn*2+1] = fzero;
                     imaskc[kcolmn] = icycle_int;
                     imaskchar[kcolmn] = jjchar;
                  } else {
                     if (imaskchar[kcolmn] > jjchar) imaskchar[kcolmn] = jjchar;
                  }
                  fmaskc[kcolmn*2] -= (_Flt)(auxU * _a_lu[k*2]);
                  fmaskc[kcolmn*2+1] -= (_Flt)(auxL * _a_lu[k*2+1]);
               }
            }
            else 
            {
               for (k=(int)iv[irwprv];k<_ia_lu[irwprv+1];k++) 
               {
                  kcolmn = (int)_ja_lu[k];
                  if (kcolmn >= 0) 
                  {
                     jjchar2 = _ja_char_lu[k]+1;
                     jjchar = (jjchar1 < jjchar2) ? jjchar2 : jjchar1;
                     if (imaskc[kcolmn] != icycle_int) 
                     {
                        listc[nlistcloc] = kcolmn;
                        nlistcloc++;
                        fmaskc[kcolmn*2] = fzero;
                        fmaskc[kcolmn*2+1] = fzero;
                        imaskc[kcolmn] = icycle_int;
                        imaskchar[kcolmn] = jjchar;
                     } else {
                        if (imaskchar[kcolmn] > jjchar) imaskchar[kcolmn] = jjchar;
                     }
                     fmaskc[kcolmn*2] -= (_Flt)(auxU * _a_lu[k*2]);
                     fmaskc[kcolmn*2+1] -= (_Flt)(auxL * _a_lu[k*2+1]);
                  }
               }
            }

         }

         irwprv = irwpr1;

      }

// Update transposed structure

      irwprv = (int)ibegm[irow];

      while (irwprv != -1) 
      {

         irwpr1 = (int)madj[irwprv];

         if (iv[irwprv] >= _ia_lu[irwprv+1]-1) 
         {
            madj[irwprv] = -1;
         } 
         else 
         {
            j = (int)iv[irwprv]+1;
            jcolmn = (int)_ja_lu[j];
            if (jcolmn < 0) jcolmn = -jcolmn-1;
            madj[irwprv] = ibegm[jcolmn];
            ibegm[jcolmn] = irwprv;
         }

         iv[irwprv]++;

         irwprv = irwpr1;

      }

// Perform filtering of the data and modify diagonal

      if (i<_n_ini) 
      {

         nlistcnew = 0;

         for (j=0;j<nlistcloc;j++) 
         {
            ind = (int)listc[j];
            if (ind != irow) 
            {
               jjchar = imaskchar[ind];
               aux1 = fmaskc[ind*2];
               aux2 = fmaskc[ind*2+1];
               if (aux1 < 0.0) aux1 = -aux1;
               if (aux2 < 0.0) aux2 = -aux2;
               dfnorm = (aux1 > aux2) ? aux1 : aux2;
               if (jjchar > 0 && (dfnorm < _tau2 || jjchar > _fcttype)) {
                  a_dia[irow] += (_Flt)dfnorm;
                  a_dia[ind] += (_Flt)dfnorm;
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

      if (fmaskc[irow*2] > fzero) 
      {
         fmaskc[irow*2] += (_Flt) (a_dia[irow]*_theta);
      } else {
         fmaskc[irow*2] -= (_Flt) (a_dia[irow]*_theta);
      }

// Factorize current row and split into first and second order

      sort (plistc, plistc+nlistcloc);

      if (i<_n_ini) 
      {
         aux1 = fmaskc[irow*2];

         if (aux1 > fzero) {
            if (aux1 < eigmin_att) eigmin_att = aux1;
            if (aux1 > eigmax_att) eigmax_att = aux1;
            if (aux1 < _pivmin) aux1 = _pivmin;
            aux1 = sqrt(aux1);
            aux1 = fone / aux1;
            auxL = aux1;
            auxU = aux1;
         } else {
            aux1 = -aux1;
            if (aux1 < eigmin_att) eigmin_att = aux1;
            if (aux1 > eigmax_att) eigmax_att = aux1;
            if (aux1 < _pivmin) aux1 = _pivmin;
            aux1 = sqrt(aux1);
            aux1 = fone / aux1;
            auxL = -aux1;
            auxU = aux1;
         }
         fmaskc[irow*2] = (_Flt)(auxL);
         fmaskc[irow*2+1] = (_Flt)(auxU);
      }

      if (i<_n_ini) 
      {
         for (j=0;j<nlistcloc;j++) {
            ind = (int)listc[j];
            if (ind != irow) {
               fmaskc[ind*2] *= (_Flt)(auxU);
               fmaskc[ind*2+1] *= (_Flt)(auxL);
               aux1 = fmaskc[ind*2];
               aux2 = fmaskc[ind*2+1];
               if (aux1 < 0.0) aux1 = -aux1;
               if (aux2 < 0.0) aux2 = -aux2;
               dfnorm = (aux1 > aux2) ? aux1 : aux2;
               jjchar = imaskchar[ind];
               if (jjchar > 0 && dfnorm < _tau1) {
                  listc[j] = (_Int) (-ind-1);
               }
            }
         }
      }

// Store computed row elems

      for (j=0;j<nlistcloc;j++) {
         ind = (int)plistc[j];
         _ja_lu.push_back (ind);
         ind1 = ind;
         if (ind1 < 0) ind1 = -ind1-1;
         jjchar = imaskchar[ind1];
         _ja_char_lu.push_back (jjchar);
         _a_lu.push_back (fmaskc[ind1*2]);
         _a_lu.push_back (fmaskc[ind1*2+1]);
      };

      nzja_lu += nlistcloc;

      _ia_lu[i+1] = nzja_lu;

// Add current row into the transposed structures

      if (nlistcloc > 1) 
      {
         if (i < _n_ini) {
            iv[irow] = _ia_lu[i]+1;
            ind = (int)listc[1];
            if (ind < 0) ind = -ind-1;
            madj[irow] = ibegm[ind];
            ibegm[ind] = (_Int) irow;
         } else {
            iv[irow] = _ia_lu[i+1]+1;
         }
      }

   }

   _eigmin_att = eigmin_att;
   _eigmax_att = eigmax_att;

// Condense the result (filter second order elems)

   int nzja_new = 0;

   for (i=0;i<nzja_lu;i++) 
   {
      if (_ja_lu[i] >= 0) nzja_new++;
   }

   vector<_Int> jlu_cnd (nzja_new+1);
   vector<char> jlu_char_cnd (nzja_new+1);
   vector<_Flt> lu_cnd (2*nzja_new+1);

   nzja_new = 0;

   listc[0] = 0;

   for (i=0;i<_n;i++) 
   {
      for (j=(int)_ia_lu[i];j<_ia_lu[i+1];j++) 
      {
         if (_ja_lu[j] >= 0) 
         {
            jlu_cnd[nzja_new] = _ja_lu[j];
            jlu_char_cnd[nzja_new] = _ja_char_lu[j];
            lu_cnd[nzja_new*2] = _a_lu[j*2];
            lu_cnd[nzja_new*2+1] = _a_lu[j*2+1];
            nzja_new++;
         }
      }
      listc[i+1] = nzja_new;
   }

   for (i=0;i<=_n;i++) _ia_lu[i] = listc[i];

   _ja_lu.swap (jlu_cnd);
   _ja_char_lu.swap (jlu_char_cnd);
   _a_lu.swap (lu_cnd);

};

//
// Balance diagonal
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2_impl<_Int,_Flt>::BalanceDiag (int _n, vector<_Int> &_ia_a, vector<_Int> &_ja_a, vector<_Flt> &_a_a,
                                             vector<_Int> &_ia_l, vector<_Int> &_ja_l, vector<_Flt> &_a_l,
                                             vector<_Int> &_ia_u, vector<_Int> &_ja_u, vector<_Flt> &_a_u,
                                             double &_diacorr_min, double &_diacorr_max) 
{

// Get arrays

   _Int *pia_a  = &_ia_a[0];
   _Int *pja_a  = &_ja_a[0];
   _Flt *pa_a = &_a_a[0];

   _Int *pia_l  = &_ia_l[0];
   _Int *pja_l  = &_ja_l[0];
   _Flt *pa_l = &_a_l[0];

   _Int *pia_u  = &_ia_u[0];
   _Int *pja_u  = &_ja_u[0];
   _Flt *pa_u = &_a_u[0];

// Compute sum array

   vector<_Flt> sum_arr (_n+1);
   _Flt *psum_arr = &sum_arr[0];

   _Flt fzero = (_Flt) 0.0e0;
   _Flt fone = (_Flt) 1.0e0;

   int i;

   for (i=0;i<_n;i++) psum_arr[i] = fzero;

   int ipL, ipU, jjL, jjU, iendL, iendU;

   for (i=0;i<_n;i++) {
      iendL = (int)pia_l[i+1]-1;
      iendU = (int)pia_u[i+1]-1;
      ipL = (int)pia_l[i]+1;
      ipU = (int)pia_u[i]+1;
      while (ipL <= iendL && ipU <= iendU) {
         jjL = (int)pja_l[ipL];
         jjU = (int)pja_u[ipU];
         if (jjL == jjU) {
            psum_arr[jjL] += (pa_l[ipL]*pa_u[ipU]);
            ipL++;
            ipU++;
         } else if (jjL < jjU) {
            ipL++;
         } else {
            ipU++;
         }
      }
   }

// Correct diagonal values and get stat data

   _diacorr_min = 1.0e20;
   _diacorr_max = -1.0e20;

   int j, jj;
   _Flt diagA, invDiagL, invDiagU, sum, aux, lu_inv_prod;
   double daux;
   bool is_found;

   for (i=0;i<_n;i++) {
      is_found = false;
      for (j=(int)pia_a[i];j<pia_a[i+1];j++) {
         jj = (int)pja_a[j];
         if (jj == i) {
            is_found = true;
            diagA = pa_a[j];
         }
      }
      if (!is_found) {
         throw " CIlu2_impl<_Int,_Flt>::BalanceDiag: error: diagonal value is not found! ";
      }
      ipL = (int)pia_l[i];
      invDiagL = pa_l[ipL];
      ipU = (int)pia_u[i];
      invDiagU = pa_u[ipU];
      sum = diagA-psum_arr[i];
      lu_inv_prod = invDiagL * invDiagU;
      aux = sum * lu_inv_prod;
      if (aux < _diacorr_min) _diacorr_min = aux;
      if (aux > _diacorr_max) _diacorr_max = aux;
      if (aux < 0.0e0) aux = 1.0e0;
      daux = (double)aux;
      daux = sqrt(daux);
      aux = (_Flt)daux;
      aux = fone / aux;
      pa_l[ipL] *= aux;
      pa_u[ipU] *= aux;
   }

};

//
// Multiply by super sparse matrix by rows and add result into prescribed positions
//========================================================================================
template <typename _Int, typename _Flt, typename _FltVect>
void CMvmSlv_impl<_Int,_Flt,_FltVect>::MvmA (int _nlist, const vector<_Int> &_list_alu, 
                                             const vector<_Int> &_ia_alu, const vector<_Int> &_ja_alu, const vector<_Flt> &_a_alu,
                                             const _FltVect *_x, _FltVect *_ax)
{

   const _Int *p_list_a = &_list_alu[0];
   const _Int *p_ia_a = &_ia_alu[0];
   const _Int *p_ja_a = &_ja_alu[0];
   const _Flt *p_a_a = &_a_alu[0];

   int i, irow, j, jj;

   for (i=0;i<_nlist;i++) {
      irow = (int)p_list_a[i];
      for (j=(int)p_ia_a[i];j<p_ia_a[i+1];j++) {
         jj = (int)p_ja_a[j];
         _ax[irow] += p_a_a[j]*_x[jj];
      }
   }

};

//
// Multiply by super sparse matrix by columns and add result into prescribed positions
//========================================================================================
template <typename _Int, typename _Flt, typename _FltVect>
void CMvmSlv_impl<_Int,_Flt,_FltVect>::MvmAT (int _nlist, const vector<_Int> &_list_alu, 
                                                const vector<_Int> &_ia_alu, const vector<_Int> &_ja_alu, const vector<_Flt> &_a_alu,
                                                const _FltVect *_x, _FltVect *_ax)
{

   const _Int *p_list_a = &_list_alu[0];
   const _Int *p_ia_a = &_ia_alu[0];
   const _Int *p_ja_a = &_ja_alu[0];
   const _Flt *p_a_a = &_a_alu[0];

   int i, irow, j, jj;

   for (i=0;i<_nlist;i++) {
      irow = (int)p_list_a[i];
      for (j=(int)p_ia_a[i];j<p_ia_a[i+1];j++) {
         jj = (int)p_ja_a[j];
         _ax[jj] += p_a_a[j]*_x[irow];
      }
   }

};

//
// Solve with L, L is stored by columns (diag is inverted)
//========================================================================================
template <typename _Int, typename _Flt, typename _FltVect>
void CMvmSlv_impl<_Int,_Flt,_FltVect>::SolveL (int _n, const vector<_Int> &_ia_l, const vector<_Int> &_ja_l, const vector<_Flt> &_a_l,
                                                const _FltVect *_x, _FltVect *_lx)
{

   int i;

   for (i=0;i<_n;i++) _lx[i] = _x[i];

   const _Int *p_ia_l = &_ia_l[0];
   const _Int *p_ja_l = &_ja_l[0];
   const _Flt *p_a_l = &_a_l[0];

   int j, jj, ibeg, iend;

   for (i=0;i<_n;i++) {
      ibeg = (int)p_ia_l[i];
      iend = (int)p_ia_l[i+1]-1;
      _lx[i] *= p_a_l[ibeg];
      for (j=ibeg+1;j<=iend;j++) {
         jj = (int)p_ja_l[j];
         _lx[jj] -= p_a_l[j]*_lx[i];
      }
   }

};

//
// Solve with U, U is stored by rows (diag is inverted)
//========================================================================================
template <typename _Int, typename _Flt, typename _FltVect>
void CMvmSlv_impl<_Int,_Flt,_FltVect>::SolveU (int _n, const vector<_Int> &_ia_u, const vector<_Int> &_ja_u, const vector<_Flt> &_a_u,
                                                const _FltVect *_x, _FltVect *_ux)
{

   int i;

   const _Int *p_ia_u = &_ia_u[0];
   const _Int *p_ja_u = &_ja_u[0];
   const _Flt *p_a_u = &_a_u[0];

   for (i=0;i<_n;i++) _ux[i] = _x[i];

   int j, jj, ibeg, iend;

   for (i=_n-1;i>=0;i--) {
      ibeg = (int)p_ia_u[i];
      iend = (int)p_ia_u[i+1]-1;
      for (j=ibeg+1;j<=iend;j++) {
         jj = (int)p_ja_u[j];
         _ux[i] -= p_a_u[j]*_ux[jj];
      }
      _ux[i] *= p_a_u[ibeg];
   }

};

// Copy constructor
//========================================================================================
template <typename _Int, typename _Flt>
CMatrix<_Int,_Flt>::CMatrix (const CMatrix<_Int,_Flt> &_aa) {

   int nlistloc = _aa.GetNlist ();
   int nlist2loc = _aa.GetNlist2 ();
   int nzjaloc = _aa.GetNzja ();
   int nzja2loc = _aa.GetNzja2 ();
   int nzaloc = _aa.GetNza ();

   const _Int *plist_aa = &(_aa.list_matr[0]);
   const _Int *plist2_aa = &(_aa.list2_matr[0]);
   const _Int *pia_aa = &(_aa.ia_matr[0]);
   const _Int *pja_aa = &(_aa.ja_matr[0]);
   const _Int *pja2_aa = &(_aa.ja2_matr[0]);
   const _Flt *pa_aa = &(_aa.a_matr[0]);

   this->n_list = nlistloc;
   this->n_list2 = nlist2loc;
   this->nz_ja = nzjaloc;
   this->nz_ja2 = nzja2loc;
   this->nz_a = nzaloc;

   this->list_matr.resize(nlistloc+1);
   this->list2_matr.resize(nlist2loc+1);
   this->ia_matr.resize(nlistloc+1);
   this->ja_matr.resize(nzjaloc+1);
   this->ja2_matr.resize(nzja2loc+1);
   this->a_matr.resize(nzaloc+1);

   _Int *plist = &(this->list_matr[0]);
   _Int *plist2 = &(this->list2_matr[0]);
   _Int *pia = &(this->ia_matr[0]);
   _Int *pja = &(this->ja_matr[0]);
   _Int *pja2 = &(this->ja2_matr[0]);
   _Flt *pa = &(this->a_matr[0]);

   int i;

   for (i=0;i<nlistloc;i++) plist[i] = plist_aa[i];
   for (i=0;i<=nlistloc;i++) pia[i] = pia_aa[i];
   for (i=0;i<nlist2loc;i++) plist2[i] = plist2_aa[i];
   for (i=0;i<nzjaloc;i++) pja[i] = pja_aa[i];
   for (i=0;i<nzja2loc;i++) pja2[i] = pja2_aa[i];

   for (i=0;i<nzaloc;i++) pa[i] = pa_aa[i];

};

// Equality operator
//========================================================================================
template <typename _Int, typename _Flt>
CMatrix<_Int,_Flt> &CMatrix<_Int,_Flt>::operator= (const CMatrix<_Int,_Flt> &_aa) {

   int nlistloc = _aa.GetNlist ();
   int nlist2loc = _aa.GetNlist2 ();
   int nzjaloc = _aa.GetNzja ();
   int nzja2loc = _aa.GetNzja2 ();
   int nzaloc = _aa.GetNza ();

   const _Int *plist_aa = &(_aa.list_matr[0]);
   const _Int *plist2_aa = &(_aa.list2_matr[0]);
   const _Int *pia_aa = &(_aa.ia_matr[0]);
   const _Int *pja_aa = &(_aa.ja_matr[0]);
   const _Int *pja2_aa = &(_aa.ja2_matr[0]);
   const _Flt *pa_aa = &(_aa.a_matr[0]);

   this->n_list = nlistloc;
   this->n_list2 = nlist2loc;
   this->nz_ja = nzjaloc;
   this->nz_ja2 = nzja2loc;
   this->nz_a = nzaloc;

   this->list_matr.resize(nlistloc+1);
   this->list2_matr.resize(nlist2loc+1);
   this->ia_matr.resize(nlistloc+1);
   this->ja_matr.resize(nzjaloc+1);
   this->ja2_matr.resize(nzja2loc+1);
   this->a_matr.resize(nzaloc+1);

   _Int *plist = &(this->list_matr[0]);
   _Int *plist2 = &(this->list2_matr[0]);
   _Int *pia = &(this->ia_matr[0]);
   _Int *pja = &(this->ja_matr[0]);
   _Int *pja2 = &(this->ja2_matr[0]);
   _Flt *pa = &(this->a_matr[0]);

   int i;

   for (i=0;i<nlistloc;i++) plist[i] = plist_aa[i];
   for (i=0;i<=nlistloc;i++) pia[i] = pia_aa[i];
   for (i=0;i<nlist2loc;i++) plist2[i] = plist2_aa[i];
   for (i=0;i<nzjaloc;i++) pja[i] = pja_aa[i];
   for (i=0;i<nzja2loc;i++) pja2[i] = pja2_aa[i];

   for (i=0;i<nzaloc;i++) pa[i] = pa_aa[i];

   return *this;

};

// Get sparsity
//========================================================================================
template <typename _Int, typename _Flt>
void CMatrix<_Int,_Flt>::GetSparsity (const CMatrix<_Int,_Flt> &_a_sp) {

   const int nlistloc = _a_sp.GetNlist ();
   const int nlist2loc = _a_sp.GetNlist2 ();
   const int nzjaloc = _a_sp.GetNzja ();
   const int nzja2loc = _a_sp.GetNzja2 ();

   const _Int *plist_aa = &(_a_sp.list_matr[0]);
   const _Int *plist2_aa = &(_a_sp.list2_matr[0]);
   const _Int *pia_aa = &(_a_sp.ia_matr[0]);
   const _Int *pja_aa = &(_a_sp.ja_matr[0]);
   const _Int *pja2_aa = &(_a_sp.ja2_matr[0]);

   this->n_list = nlistloc;
   this->n_list2 = nlist2loc;
   this->nz_ja = nzjaloc;
   this->nz_ja2 = nzja2loc;
   this->nz_a = 0;

   this->list_matr.resize(nlistloc+1);
   this->list2_matr.resize(nlist2loc+1);
   this->ia_matr.resize(nlistloc+1);
   this->ja_matr.resize(nzjaloc+1);
   this->ja2_matr.resize(nzja2loc+1);
   this->a_matr.resize(1);

   _Int *plist = &(this->list_matr[0]);
   _Int *plist2 = &(this->list2_matr[0]);
   _Int *pia = &(this->ia_matr[0]);
   _Int *pja = &(this->ja_matr[0]);
   _Int *pja2 = &(this->ja2_matr[0]);

   int i;

   for (i=0;i<nlistloc;i++) plist[i] = plist_aa[i];
   for (i=0;i<=nlistloc;i++) pia[i] = pia_aa[i];
   for (i=0;i<nlist2loc;i++) plist2[i] = plist2_aa[i];
   for (i=0;i<nzjaloc;i++) pja[i] = pja_aa[i];
   for (i=0;i<nzja2loc;i++) pja2[i] = pja2_aa[i];

};

// Compute transposed sparsity for incomplete list
//========================================================================================
template <typename _Int, typename _Flt>
void CMatrix<_Int,_Flt>::TransposedSparsityListSp (int &_icycle, int *_imask, 
                                                         int *_indarr, int *_iptr, int *_listloc, int *_ialoc,
                                                         CMatrix<_Int,_Flt> &_at) const
{

   const int nlist_ini = this->n_list;
   const int nzja_ini = this->nz_ja;
   const _Int *plist = &this->list_matr[0];
   const _Int *pia = &this->ia_matr[0];
   const _Int *pja = &this->ja_matr[0];

// Count the sparsity of at

   int nlistloc = 0;

   _icycle++;

   int i, j, jj;

   for (i=0;i<nlist_ini;i++) {
      for (j=(int)pia[i];j<pia[i+1];j++) {
         jj = (int)pja[j];
         if (_imask[jj] != _icycle) {
            _listloc[nlistloc] = jj;
            nlistloc++;
            _imask[jj] = _icycle;
            _indarr[jj] = 0;
         }
         _indarr[jj]++;
      }
   }

   sort (_listloc,_listloc+nlistloc);

   _ialoc[0] = 0;

   for (i=0;i<nlistloc;i++) {
      jj = _listloc[i];
      _ialoc[i+1] = _ialoc[i]+_indarr[jj];
      _iptr[i] = _ialoc[i];
      _indarr[jj] = i;
   }

// Init transposed matrix

   CMatrix<_Int,_Flt> at;

   at.ResizeAndSetAllSp (nlistloc, 0, nzja_ini, 0);

   _Int *plist_at = at.GetListArr ();
   _Int *pia_at = at.GetIaArr ();
   _Int *pja_at = at.GetJaArr ();

   for (i=0;i<nlistloc;i++) plist_at[i] = (_Int)_listloc[i];
   for (i=0;i<=nlistloc;i++) pia_at[i] = (_Int)_ialoc[i];

   int k, irow, ind;

   for (i=0;i<nlist_ini;i++) {
      irow = (int)plist[i];
      for (j=(int)pia[i];j<pia[i+1];j++) {
         jj = (int)pja[j];
         ind = _indarr[jj];
         _iptr[ind]++;
         k = _iptr[ind];
         pja_at[k-1] = (_Int)irow;
      }
   }

   _at.ReplaceFree (at);

}

// Compute transposed matrix for incomplete list
//========================================================================================
template <typename _Int, typename _Flt>
void CMatrix<_Int,_Flt>::TransposedSparsityList (int &_icycle, int *_imask, 
                                                         int *_indarr, int *_iptr, int *_listloc, int *_ialoc,
                                                         CMatrix<_Int,_Flt> &_at) const
{

   const int nlist_ini = this->n_list;
   const int nzja_ini = this->nz_ja;
   const _Int *plist = &this->list_matr[0];
   const _Int *pia = &this->ia_matr[0];
   const _Int *pja = &this->ja_matr[0];
   const _Flt *pa = &this->a_matr[0];

// Count the sparsity of at

   int nlistloc = 0;

   _icycle++;

   int i, j, jj;

   for (i=0;i<nlist_ini;i++) {
      for (j=(int)pia[i];j<pia[i+1];j++) {
         jj = (int)pja[j];
         if (_imask[jj] != _icycle) {
            _listloc[nlistloc] = jj;
            nlistloc++;
            _imask[jj] = _icycle;
            _indarr[jj] = 0;
         }
         _indarr[jj]++;
      }
   }

   sort (_listloc,_listloc+nlistloc);

   _ialoc[0] = 0;

   for (i=0;i<nlistloc;i++) {
      jj = _listloc[i];
      _ialoc[i+1] = _ialoc[i]+_indarr[jj];
      _iptr[i] = _ialoc[i];
      _indarr[jj] = i;
   }

// Init transposed matrix

   CMatrix<_Int,_Flt> at;

   at.ResizeAndSetAll (nlistloc, 0, nzja_ini, 0, nzja_ini);

   _Int *plist_at = at.GetListArr ();
   _Int *pia_at = at.GetIaArr ();
   _Int *pja_at = at.GetJaArr ();
   _Flt *pa_at = at.GetAArr ();

   for (i=0;i<nlistloc;i++) plist_at[i] = (_Int)_listloc[i];
   for (i=0;i<=nlistloc;i++) pia_at[i] = (_Int)_ialoc[i];

   int k, irow, ind;

   for (i=0;i<nlist_ini;i++) {
      irow = (int)plist[i];
      for (j=(int)pia[i];j<pia[i+1];j++) {
         jj = (int)pja[j];
         ind = _indarr[jj];
         _iptr[ind]++;
         k = _iptr[ind];
         pja_at[k-1] = (_Int)irow;
         pa_at[k-1] = pa[j];
      }
   }

   _at.ReplaceFree (at);

}

// Extend sparsity by zeroes
//========================================================================================
template <typename _Int, typename _Flt>
void CMatrix<_Int,_Flt>::ExtendSparsity (const CMatrix<_Int,_Flt> &_a_sp,
                                                CMatrix<_Int,_Flt> &_a_ext) const
{

// Open sparsities

   const int nlist_ini = this->n_list;
   const int nzja_ini = this->nz_ja;
   const _Int *plist = &this->list_matr[0];
   const _Int *pia = &this->ia_matr[0];
   const _Int *pja = &this->ja_matr[0];
   const _Flt *pa = &this->a_matr[0];

   const int nlist_sp = _a_sp.GetNlist();
   const int nzja_sp = _a_sp.GetNzja();
   const _Int *plist_sp = _a_sp.GetListArr();
   const _Int *pia_sp = _a_sp.GetIaArr();
   const _Int *pja_sp = _a_sp.GetJaArr();

// Compute extended by zeroes data

   int nlist_ext_max = nlist_ini+nlist_sp;
   int nzja_ext_max = nzja_ini+nzja_sp;

   vector<_Int> list_ext (nlist_ext_max+1);
   vector<_Int> ia_ext (nlist_ext_max+1);
   vector<_Int> ja_ext (nzja_ext_max+1);
   vector<_Flt> a_ext (nzja_ext_max+1);

   _Int *plist_ext = &list_ext[0];
   _Int *pia_ext = &ia_ext[0];
   _Int *pja_ext = &ja_ext[0];
   _Flt *pa_ext = &a_ext[0];

   _Flt fzero = (_Flt) 0.0e0;

   int nlist_ext = 0;
   int nzja_ext = 0;

   int ip_list, ip_list_sp, irow, irow_sp, j, jj, jj_sp, jp, jp_sp;

   ip_list = 0;
   ip_list_sp = 0;

   pia_ext[0] = 0;

   while (ip_list < nlist_ini || ip_list_sp < nlist_sp) {
      if (ip_list < nlist_ini && ip_list_sp < nlist_sp) {
         irow = (int)plist[ip_list];
         irow_sp = (int)plist_sp[ip_list_sp];
         if (irow == irow_sp) {
            plist_ext[nlist_ext] = irow;
            nlist_ext++;
            jp = (int)pia[ip_list];
            jp_sp = (int)pia_sp[ip_list_sp];
            while (jp < pia[ip_list+1] || jp_sp < pia_sp[ip_list_sp+1]) {
               if (jp < pia[ip_list+1] && jp_sp < pia_sp[ip_list_sp+1]) {
                  jj = (int)pja[jp];
                  jj_sp = (int)pja_sp[jp_sp];
                  if (jj == jj_sp) {
                     pja_ext[nzja_ext] = pja[jp];
                     pa_ext[nzja_ext] = pa[jp];
                     nzja_ext++;
                     jp++;
                     jp_sp++;
                  } else if (jj < jj_sp) {
                     pja_ext[nzja_ext] = pja[jp];
                     pa_ext[nzja_ext] = pa[jp];
                     nzja_ext++;
                     jp++;
                  } else {
                     pja_ext[nzja_ext] = pja_sp[jp_sp];
                     pa_ext[nzja_ext] = fzero;
                     nzja_ext++;
                     jp_sp++;
                  }
               } else if (jp < pia[ip_list+1]) {
                  pja_ext[nzja_ext] = pja[jp];
                  pa_ext[nzja_ext] = pa[jp];
                  nzja_ext++;
                  jp++;
               } else {
                  pja_ext[nzja_ext] = pja_sp[jp_sp];
                  pa_ext[nzja_ext] = fzero;
                  nzja_ext++;
                  jp_sp++;
               }
            }
            ip_list++;
            ip_list_sp++;
         } else if (irow < irow_sp) {
            plist_ext[nlist_ext] = irow;
            nlist_ext++;
            for (j=(int)pia[ip_list];j<pia[ip_list+1];j++) {
               pja_ext[nzja_ext] = pja[j];
               pa_ext[nzja_ext] = pa[j];
               nzja_ext++;
            }
            ip_list++;
         } else {
            plist_ext[nlist_ext] = irow_sp;
            nlist_ext++;
            for (j=(int)pia_sp[ip_list_sp];j<pia_sp[ip_list_sp+1];j++) {
               pja_ext[nzja_ext] = pja_sp[j];
               pa_ext[nzja_ext] = fzero;
               nzja_ext++;
            }
            ip_list_sp++;
         }
      } else if (ip_list < nlist_ini) {
         plist_ext[nlist_ext] = plist[ip_list];
         nlist_ext++;
         for (j=(int)pia[ip_list];j<pia[ip_list+1];j++) {
            pja_ext[nzja_ext] = pja[j];
            pa_ext[nzja_ext] = pa[j];
            nzja_ext++;
         }
         ip_list++;
      } else if (ip_list_sp < nlist_sp) {
         plist_ext[nlist_ext] = plist_sp[ip_list_sp];
         nlist_ext++;
         for (j=(int)pia_sp[ip_list_sp];j<pia_sp[ip_list_sp+1];j++) {
            pja_ext[nzja_ext] = pja_sp[j];
            pa_ext[nzja_ext] = fzero;
            nzja_ext++;
         }
         ip_list_sp++;
      }
      pia_ext[nlist_ext] = nzja_ext;
   }

// Store computed data

   CMatrix<_Int,_Flt> a_new;

   a_new.ResizeAndSetAll (nlist_ext, 0, nzja_ext, 0, nzja_ext);

   _Int *plist_new = a_new.GetListArr();
   _Int *pia_new = a_new.GetIaArr();
   _Int *pja_new = a_new.GetJaArr();
   _Flt *pa_new = a_new.GetAArr();

   int i;

   for (i=0;i<nlist_ext;i++) plist_new[i] = plist_ext[i];
   for (i=0;i<=nlist_ext;i++) pia_new[i] = pia_ext[i];
   for (i=0;i<nzja_ext;i++) pja_new[i] = pja_ext[i];
   for (i=0;i<nzja_ext;i++) pa_new[i] = pa_ext[i];

   _a_ext.ReplaceFree(a_new);

}

// Perform multiplication of matrices stored by rows, resulting matrix is also stored by rows
//========================================================================================
template <typename _Int, typename _Flt>
void CMatrix<_Int,_Flt>::MatrixByMatrixMultiply (int &_icycle, int *_imask, int *_imask1, int *_indarr, int *_listloc, _Flt *_fmask,
                                                      CMatrix<_Int,_Flt> &_a, CMatrix<_Int,_Flt> &_b,
                                                      CMatrix<_Int,_Flt> &_a_times_b)
{

// Open sparsities of A and B

   int nlist_a = _a.GetNlist();
   _Int *plist_a = _a.GetListArr();
   _Int *pia_a = _a.GetIaArr();
   _Int *pja_a = _a.GetJaArr();
   _Flt *pa_a = _a.GetAArr();

   int nlist_b = _b.GetNlist();
   _Int *plist_b = _b.GetListArr();
   _Int *pia_b = _b.GetIaArr();
   _Int *pja_b = _b.GetJaArr();
   _Flt *pa_b = _b.GetAArr();

// Perform multiplication

   vector<_Int> list_ab (nlist_a+1);
   vector<_Int> ia_ab (nlist_a+1);

   _Int *plist_ab = &list_ab[0];
   _Int *pia_ab = &ia_ab[0];

   vector<_Int> ja_ab (1);
   vector<_Flt> a_ab (1);

   int i;

   for (i=0;i<nlist_a;i++) plist_ab[i] = plist_a[i];

   pia_ab[0] = 0;

   _icycle++;

   int icycle1 = _icycle;

   int jj;

   for (i=0;i<nlist_b;i++) {
      jj = (int)plist_b[i];
      _imask1[jj] = icycle1;
      _indarr[jj] = i;
   }

   _Flt fzero = (_Flt) 0.0e0;

   int nzja_ab = 0;

   int j, k, kk, nlistloc, ind;
   _Flt aux;

   for (i=0;i<nlist_a;i++) {
      _icycle++;
      nlistloc = 0;
      for (j=(int)pia_a[i];j<pia_a[i+1];j++) {
         jj = (int)pja_a[j];
         if (_imask1[jj] == icycle1) {
            ind = _indarr[jj];
            for (k=(int)pia_b[ind];k<pia_b[ind+1];k++) {
               kk = (int)pja_b[k];
               if (_imask[kk] != _icycle) {
                  _listloc[nlistloc] = (int)kk;
                  nlistloc++;
                  _imask[kk] = _icycle;
               }
            }
         }
      }
      sort (_listloc,_listloc+nlistloc);
      for (j=0;j<nlistloc;j++) {
         jj = _listloc[j];
         _fmask[jj] = fzero;
      }
      for (j=(int)pia_a[i];j<pia_a[i+1];j++) {
         jj = (int)pja_a[j];
         if (_imask1[jj] == icycle1) {
            aux = pa_a[j];
            ind = _indarr[jj];
            for (k=(int)pia_b[ind];k<pia_b[ind+1];k++) {
               kk = (int)pja_b[k];
               _fmask[kk] += aux*pa_b[k];
            }
         }
      }
      ja_ab.resize (nzja_ab+nlistloc+1);
      a_ab.resize (nzja_ab+nlistloc+1);
      for (j=0;j<nlistloc;j++) ja_ab[nzja_ab+j] = (_Int)_listloc[j];
      for (j=0;j<nlistloc;j++) {
         jj = _listloc[j];
         a_ab[nzja_ab+j] = _fmask[jj];
      }
      nzja_ab += nlistloc;
      pia_ab[i+1] = nzja_ab;
   }

// Store result

   vector<_Int> *p_list_ab = _a_times_b.GetList();
   vector<_Int> *p_ia_ab = _a_times_b.GetIa();
   vector<_Int> *p_ja_ab = _a_times_b.GetJa();
   vector<_Flt> *p_a_ab = _a_times_b.GetA();

   p_list_ab->swap (list_ab);
   p_ia_ab->swap (ia_ab);
   p_ja_ab->swap (ja_ab);
   p_a_ab->swap (a_ab);

   _a_times_b.SetNlist (nlist_a);
   _a_times_b.SetNzja (nzja_ab);
   _a_times_b.SetNza (nzja_ab);

}

// Compute the packed size
//========================================================================================
template <typename _Int, typename _Flt>
int CMatrix<_Int,_Flt>::GetPackedSize () const {

   int nlistloc = this->n_list;
   int nlist2loc = this->n_list2;
   int nzjaloc = this->nz_ja;
   int nzja2loc = this->nz_ja2;
   int nzaloc = this->nz_a;

   int isize = 5*sizeof(int) + (2*nlistloc+1+nlist2loc+nzjaloc+nzja2loc)*sizeof(_Int) + nzaloc*sizeof(_Flt);

   return isize;

}

// Fill by packed data
//========================================================================================
template <typename _Int, typename _Flt>
void CMatrix<_Int,_Flt>::FillPacked (int _length, char *_obj) const {

   int nlistloc = this->n_list;
   int nlist2loc = this->n_list2;
   int nzjaloc = this->nz_ja;
   int nzja2loc = this->nz_ja2;
   int nzaloc = this->nz_a;

   char* pLoc;

   pLoc = _obj;

   int *pHead = NULL;
   _Int *plist_obj = NULL;
   _Int *plist2_obj = NULL;
   _Int *pia_obj = NULL;
   _Int *pja_obj = NULL;
   _Int *pja2_obj = NULL;
   _Flt *pa_obj = NULL;

   pHead = (int *) pLoc;
   pLoc += 5 * sizeof(int);

   plist_obj = (_Int *) pLoc;
   pLoc += nlistloc * sizeof(_Int);

   plist2_obj = (_Int *) pLoc;
   pLoc += nlist2loc * sizeof(_Int);

   pia_obj = (_Int *) pLoc;
   pLoc += (nlistloc+1) * sizeof(_Int);

   pja_obj = (_Int *) pLoc;
   pLoc += (nzjaloc) * sizeof(_Int);

   pja2_obj = (_Int *) pLoc;
   pLoc += (nzja2loc) * sizeof(_Int);

   pa_obj = (_Flt *) pLoc;
   pLoc += (nzaloc) * sizeof(_Flt);

   pHead[0]  = nlistloc;
   pHead[1]  = nlist2loc;
   pHead[2]  = nzjaloc;
   pHead[3]  = nzja2loc;
   pHead[4]  = nzaloc;

   if (pLoc-_obj != _length) {
      throw " CMatrix<_Int,_Flt>::FillPacked: incorrect length on entry ";
   }

// Fill arrays

   const _Int *plist  = &(this->list_matr[0]);
   const _Int *plist2 = &(this->list2_matr[0]);
   const _Int *pia    = &(this->ia_matr[0]);
   const _Int *pja    = &(this->ja_matr[0]);
   const _Int *pja2   = &(this->ja2_matr[0]);
   const _Flt *pa   = &(this->a_matr[0]);

   int j;

   for (j = 0; j < nlistloc;   j++) plist_obj[j]  = plist[j];
   for (j = 0; j < nlist2loc;  j++) plist2_obj[j] = plist2[j];
   for (j = 0; j < nlistloc+1; j++) pia_obj[j]    = pia[j];
   for (j = 0; j < nzjaloc;    j++) pja_obj[j]    = pja[j];
   for (j = 0; j < nzja2loc;   j++) pja2_obj[j]   = pja2[j];
   for (j = 0; j < nzaloc;     j++) pa_obj[j]     = pa[j];

}

// Fill by packed data
//========================================================================================
template <typename _Int, typename _Flt>
void CMatrix<_Int,_Flt>::UnPack (int _length, char *_obj) {

// Get head data

   char *pLoc;

   pLoc = _obj;

   int *pHead;

   pHead = (int *) pLoc;
   pLoc += 5 * sizeof(int);

   int nlistloc  = pHead[0];
   int nlist2loc = pHead[1];
   int nzjaloc   = pHead[2];
   int nzja2loc  = pHead[3];
   int nzaloc    = pHead[4];

   _Int *plist_obj = NULL;
   _Int *plist2_obj = NULL;
   _Int *pia_obj = NULL;
   _Int *pja_obj = NULL;
   _Int *pja2_obj = NULL;
   _Flt *pa_obj = NULL;

   plist_obj = (_Int *) pLoc;
   pLoc += nlistloc * sizeof(_Int);

   plist2_obj = (_Int *) pLoc;
   pLoc += nlist2loc * sizeof(_Int);

   pia_obj = (_Int *) pLoc;
   pLoc += (nlistloc+1) * sizeof(_Int);

   pja_obj = (_Int *) pLoc;
   pLoc += (nzjaloc) * sizeof(_Int);

   pja2_obj = (_Int *) pLoc;
   pLoc += (nzja2loc) * sizeof(_Int);

   pa_obj = (_Flt *) pLoc;
   pLoc += (nzaloc) * sizeof(_Flt);

// Store data

   this->n_list = nlistloc;
   this->n_list2 = nlist2loc;
   this->nz_ja = nzjaloc;
   this->nz_ja2 = nzja2loc;
   this->nz_a = nzaloc;

   this->list_matr.resize(nlistloc+1);
   this->list2_matr.resize(nlist2loc+1);
   this->ia_matr.resize(nlistloc+1);
   this->ja_matr.resize(nzjaloc+1);
   this->ja2_matr.resize(nzja2loc+1);
   this->a_matr.resize(nzaloc+1);

   int i;
   _Int *piarr;
   _Flt *paarr;

   piarr = &(this->list_matr[0]);
   for (i=0;i<nlistloc;i++) piarr[i] = plist_obj[i];

   piarr = &(this->list2_matr[0]);
   for (i=0;i<nlist2loc;i++) piarr[i] = plist2_obj[i];

   piarr = &(this->ia_matr[0]);
   for (i=0;i<=nlistloc;i++) piarr[i] = pia_obj[i];

   piarr = &(this->ja_matr[0]);
   for (i=0;i<nzjaloc;i++) piarr[i] = pja_obj[i];

   piarr = &(this->ja2_matr[0]);
   for (i=0;i<nzja2loc;i++) piarr[i] = pja2_obj[i];

   paarr = &(this->a_matr[0]);
   for (i=0;i<nzaloc;i++) paarr[i] = pa_obj[i];

};

//
// Set vector data by zeroes
//========================================================================================
template <typename _FltVect>
void CVect<_FltVect>::SetByZeroes (int _n, _FltVect *_x)
{

   int i;

   _FltVect fzero = (_FltVect)0.0e0;

   for (i=0;i<_n;i++) _x[i] = fzero;

};

//
// Compute scalar product
//========================================================================================
template <typename _FltVect>
void CVect<_FltVect>::CopyVector (int _n, _FltVect *_x, _FltVect *_y)
{

   int i;

   for (i=0;i<_n;i++) _y[i] = _x[i];

};

//
// Compute scalar product
//========================================================================================
template <typename _FltVect>
_FltVect CVect<_FltVect>::ScProd (int _n, _FltVect *_x, _FltVect *_y)
{

   _FltVect fsum = (_FltVect)0.0e0;

   int i;

   for (i=0;i<_n;i++) fsum += _x[i]*_y[i];

   return fsum;

};

//
// Add vector data
//========================================================================================
template <typename _FltVect>
void CVect<_FltVect>::AddReplaceVector (int _n, _FltVect *_x1, _FltVect *_x1plus2)
{

   int i;

   for (i=0;i<_n;i++) _x1plus2[i] += _x1[i];

};

//
// Subtract vector data
//========================================================================================
template <typename _FltVect>
void CVect<_FltVect>::SubtractReplaceVector (int _n, _FltVect *_x1, _FltVect *_x2minus1)
{

   int i;

   for (i=0;i<_n;i++) _x2minus1[i] -= _x1[i];

};

// Update array for minus alpha (axpy)
//========================================================================================
template <typename _FltVect>
void CVect<_FltVect>::UpdateVectorMinus (int _n, _FltVect *_value, _FltVect *_arr_x, _FltVect *_arr_y) {

   _FltVect value = *_value;

   int i;

   for (i=0;i<_n;i++) _arr_y[i] = _arr_y[i] - value * _arr_x[i];

};

// Update array reversed (aypx)
//========================================================================================
template <typename _FltVect>
void CVect<_FltVect>::UpdateVectorReversed (int _n, _FltVect *_value, _FltVect *_arr_x, _FltVect *_arr_y) {

   _FltVect value = *_value;

   for (int i=0;i<_n;i++) _arr_y[i] =  value * _arr_y[i] + _arr_x[i];

};

// Update array (axpy)
//========================================================================================
template <typename _FltVect>
void CVect<_FltVect>::UpdateVector (int _n, _FltVect *_value, _FltVect *_arr_x, _FltVect *_arr_y) { // Update array (axpy): y=a*x+y

   _FltVect value = *_value;

   for (int i=0;i<_n;i++) _arr_y[i] += value * _arr_x[i];

};

//
// Order vector data
//========================================================================================
template <typename _FltVect>
void CVect<_FltVect>::OrderVector (int _n, const vector<int> &_order,
                                       const _FltVect *_x, _FltVect *_x_ord)
{

   int i;

   for (i=0;i<_n;i++) _x_ord[_order[i]] = _x[i];

};

//
// Inverse order vector data
//========================================================================================
template <typename _FltVect>
void CVect<_FltVect>::InvOrderVector (int _n, const vector<int> &_order,
                                          const _FltVect *_x, _FltVect *_x_ord)
{

   int i;

   for (i=0;i<_n;i++) _x_ord[i] = _x[_order[i]];

};

//
// Print matrix data
//========================================================================================
template <typename _Int, typename _Flt>
void CIlu2<_Int,_Flt>::PrintMatrix (ofstream &_fout, CMatrix<_Int,_Flt> &_a_matr) 
{

   _fout << " CMatrix:" << endl;

   _fout << "    Nlist = " << _a_matr.GetNlist() << " Nlist2 = " << _a_matr.GetNlist2 () << endl;
   _fout << "    Nzja  = " << _a_matr.GetNzja () << " Nzja2  = " << _a_matr.GetNzja2 () << " Nza = " << _a_matr.GetNza () << endl;

   if (_a_matr.GetNlist() > 0)  PrintArray (_fout, "    List  = ",_a_matr.GetNlist(),  _a_matr.GetListArr());
   if (_a_matr.GetNlist2() > 0) PrintArray (_fout, "    List2 = ",_a_matr.GetNlist2(), _a_matr.GetList2Arr());
   if (_a_matr.GetNlist() > 0)  PrintArray (_fout, "    Ia    = ",_a_matr.GetNlist()+1,_a_matr.GetIaArr());
   if (_a_matr.GetNzja() > 0)   PrintArray (_fout, "    Ja    = ",_a_matr.GetNzja(),   _a_matr.GetJaArr());
   if (_a_matr.GetNzja2() > 0)  PrintArray (_fout, "    Ja2   = ",_a_matr.GetNzja2(),  _a_matr.GetJa2Arr());
   if (_a_matr.GetNza() > 0)    PrintArray (_fout, "    A     = ",_a_matr.GetNza(),    _a_matr.GetAArr());

};

//
// Init by data
//========================================================================================
template <typename _Int, typename _Flt, typename _Flt2>
void CMatrixConv<_Int,_Flt,_Flt2>::InitAndConv (int _n, _Int *_ia, _Int *_ja, _Flt *_a, CMatrix<_Int,_Flt2> &_amatr) 
{

   CMatrix<_Int,_Flt2> a_temp (_n, _ia, _ja);

   vector<_Flt2> *pa = a_temp.GetA();

   int nzja = (int)_ia[_n];

   pa->resize(nzja+1);

   _Flt2 *p_ptr_a = &((*pa)[0]);

   int i;

   for (i=0;i<nzja;i++) p_ptr_a[i] = (_Flt2)_a[i];

   _amatr.ReplaceFree (a_temp);

};

// Copy constructor
//========================================================================================
template <typename _Int, typename _Flt>
CHMatrix<_Int,_Flt>::CHMatrix (const CHMatrix<_Int,_Flt> &_aa) {

   this->nzblk = _aa.nzblk;
   this->hmatr_str = _aa.hmatr_str;
   this->asub_arr.resize (_aa.nzblk+1);
   const CMatrix<_Int,_Flt> *pasub_arr_aa = &(_aa.asub_arr[0]);
   CMatrix<_Int,_Flt> *pasub_arr = &(this->asub_arr[0]);
   for (int i=0;i<_aa.nzblk;i++) pasub_arr[i] = pasub_arr_aa[i];

};

// Equality operator
//========================================================================================
template <typename _Int, typename _Flt>
CHMatrix<_Int,_Flt> &CHMatrix<_Int,_Flt>::operator= (const CHMatrix<_Int,_Flt> &_aa) {

   this->nzblk = _aa.nzblk;
   this->hmatr_str = _aa.hmatr_str;
   this->asub_arr.resize (_aa.nzblk+1);
   const CMatrix<_Int,_Flt> *pasub_arr_aa = &(_aa.asub_arr[0]);
   CMatrix<_Int,_Flt> *pasub_arr = &(this->asub_arr[0]);
   for (int i=0;i<_aa.nzblk;i++) pasub_arr[i] = pasub_arr_aa[i];

   return *this;

};

//
// Init by data
//========================================================================================
template <typename _Int, typename _Flt>
CHMatrix<_Int,_Flt>::CHMatrix (int _iblk, int _nlist, _Int *_list, _Int *_ia, _Int *_ja, _Flt *_a,
                                       int _nblks, long long *_blks,
                                       int &_icycle, int *_imaskblk) 
{

// Check list array

   int i;

   long long jj;

   for (i=0;i<_nlist;i++) {
      jj = _list[i];
      if (jj < _blks[_iblk] || jj >= _blks[_iblk+1]) {
         throw " CHMatrix<>::CHMatrix: error: incorrect list array ";
      }
   }

// Compute ja2 array

   int nzja = (int)_ia[_nlist];

   vector<_Int> ja2 (nzja+1);
   _Int *pja2 = &ja2[0];

   CHMatrix<_Int,_Flt>::ComputeJa2 (_nblks, _blks,
                                          _nlist, _ia, _ja, 
                                          pja2);

// Create list of block numbers

   _icycle++;

   int nlistblk = 0;

   vector<int> listblk;

   int jblk;

   for (i=0;i<nzja;i++) {
      jblk = (int)pja2[i];
      if (_imaskblk[jblk] != _icycle) {
         _imaskblk[jblk] = _icycle;
         nlistblk++;
         listblk.push_back (jblk);
      }
   }

   int *plistblk = NULL;
   if (nlistblk > 0) plistblk = &listblk[0];

   sort (plistblk,plistblk+nlistblk);

   vector<int> indlistblk (_nblks);
   int *pindlistblk = &indlistblk[0];

   for (i=0;i<nlistblk;i++) {
      jblk = (int)plistblk[i];
      pindlistblk[jblk] = i;
   }

// Split block row into blocks

   vector<int> nzblk (_nblks);
   int *pnzblk = &nzblk[0];

   for (i=0;i<nlistblk;i++) {
      jblk = plistblk[i];
      pnzblk[jblk] = 0;
   }

   for (i=0;i<nzja;i++) {
      jblk = (int)pja2[i];
      pnzblk[jblk]++;
   }

   vector<vector<_Int> >   irowarr_blk (nlistblk+1);
   vector<vector<_Int> >   icolarr_blk (nlistblk+1);
   vector<vector<_Flt> > valarr_blk  (nlistblk+1);

   vector<_Int> *pirowarr_blk  = &irowarr_blk[0];
   vector<_Int> *picolarr_blk  = &icolarr_blk[0];
   vector<_Flt> *pvalarr_blk = &valarr_blk[0];

   for (i=0;i<nlistblk;i++) {
      jblk = (int)plistblk[i];
      pirowarr_blk[i].resize (pnzblk[jblk]+1);
      picolarr_blk[i].resize (pnzblk[jblk]+1);
      pvalarr_blk [i].resize (pnzblk[jblk]+1);
   }

   for (i=0;i<nlistblk;i++) {
      jblk = plistblk[i];
      pnzblk[jblk] = 0;
   }

   _Int *ppirowarr_blk;
   _Int *ppicolarr_blk;
   _Flt *ppvalarr_blk;

   int j, ind, k;

   for (i=0;i<_nlist;i++) {
      for (j=(int)_ia[i];j<_ia[i+1];j++) {
         jj = (int)_ja[j];
         jblk = (int)pja2[j];
         jj -= (int)_blks[jblk];
         ind = (int)pindlistblk[jblk];
         ppirowarr_blk = &(pirowarr_blk[ind][0]);
         ppicolarr_blk = &(picolarr_blk[ind][0]);
         ppvalarr_blk = &(pvalarr_blk[ind][0]);
         k = (int)pnzblk[jblk];
         ppirowarr_blk[k] = (_Int)(_list[i]-_blks[_iblk]);
         ppicolarr_blk[k] = (_Int)jj;
         ppvalarr_blk [k] = _a[j];
         pnzblk[jblk]++;
      }
   }

// Store block sparsity matrix

   this->nzblk = nlistblk;

   int ia_1blk[2];

   ia_1blk[0] = 0;
   ia_1blk[1] = nlistblk;

   CMatrix<int,float> hmatr_str_temp (1, &_iblk, ia_1blk, plistblk);

   this->hmatr_str.ReplaceFree (hmatr_str_temp);

// For each block compute the list of rows and ia arrays

   int nimax = 0;

   int ni;

   for (i=0;i<nlistblk;i++) {
      jblk = plistblk[i];
      ni = (int)(_blks[jblk+1]-_blks[jblk]);
      if (ni > nimax) nimax = ni;
   }

   vector<int> imask (nimax+1);
   vector<int> list (nimax+1);
   vector<int> indlist (nimax+1);

   int *pimask = &imask[0];
   int *plist = &list[0];
   int *pindlist = &indlist[0];

   int icycle = -1;

   for (i=0;i<nimax;i++) pimask[i] = -1;

   vector<int> nlistarr_blk (nlistblk+1);
   vector<vector<_Int> > listarr_blk (nlistblk+1);
   vector<vector<_Int> > iaarr_blk (nlistblk+1);

   int *pnlistarr_blk = &nlistarr_blk[0];
   vector<_Int> *plistarr_blk = &listarr_blk[0];
   vector<_Int> *piaarr_blk = &iaarr_blk[0];

   for (i=0;i<nlistblk;i++) pnlistarr_blk[i] = 0;

   int nlistloc;
   _Int *pplistarr_blk;
   _Int *ppiaarr_blk;

   for (i=0;i<nlistblk;i++) {
      jblk = plistblk[i];
      icycle++;
      ppirowarr_blk = &(pirowarr_blk[i][0]);
      nlistloc = 0;
      for (j=0;j<pnzblk[jblk];j++) {
         jj = ppirowarr_blk[j];
         if (jj < 0 || jj > nimax) {
            cout << " Wait !!! " << endl;
         }
         if (pimask[jj] != icycle) {
            plist[nlistloc] = (int)jj;
            nlistloc++;
            pimask[jj] = icycle;
         }
      }
      pnlistarr_blk[i] = nlistloc;
      sort (plist,plist+nlistloc);
      for (j=0;j<nlistloc;j++) {
         jj = plist[j];
         pindlist[jj] = j;
      }
      plistarr_blk[i].resize(nlistloc+1);
      piaarr_blk[i].resize(nlistloc+1);
      pplistarr_blk = &(plistarr_blk[i][0]);
      ppiaarr_blk = &(piaarr_blk[i][0]);
      for (j=0;j<nlistloc;j++) pplistarr_blk[j] = plist[j];
      for (j=0;j<=nlistloc;j++) ppiaarr_blk[j] = 0;
      for (j=0;j<pnzblk[jblk];j++) {
         jj = ppirowarr_blk[j];
         ind = pindlist[jj];
         ppiaarr_blk[ind+1]++;
      }
      for (j=0;j<nlistloc;j++) ppiaarr_blk[j+1] = ppiaarr_blk[j]+ppiaarr_blk[j+1];
   }

// Store results as hmatrix data

   vector<CMatrix<_Int,_Flt> > asub_arr_temp (nlistblk+1); 

   CMatrix<_Int,_Flt> *pasub_arr_temp = &asub_arr_temp[0];

   for (i=0;i<nlistblk;i++) {
      jblk = plistblk[i];
      vector<_Int> *plist_temp = pasub_arr_temp[i].GetList();
      vector<_Int> *pia_temp = pasub_arr_temp[i].GetIa();
      vector<_Int> *pja_temp = pasub_arr_temp[i].GetJa();
      vector<_Flt> *pa_temp = pasub_arr_temp[i].GetA();
      pasub_arr_temp[i].SetNlist(pnlistarr_blk[i]);
      pasub_arr_temp[i].SetNzja(pnzblk[jblk]);
      pasub_arr_temp[i].SetNza(pnzblk[jblk]);
      plist_temp->swap(plistarr_blk[i]);
      pia_temp->swap(piaarr_blk[i]);
      pja_temp->swap(picolarr_blk[i]);
      pa_temp->swap(pvalarr_blk[i]);
   }

   this->asub_arr.resize (nlistblk+1);

   CMatrix<_Int,_Flt> *pasub_arr = &(this->asub_arr[0]);

   for (i=0;i<nlistblk;i++) pasub_arr[i].ReplaceFree(pasub_arr_temp[i]);

};

// Compute the symmetrized submatrices
//========================================================================================
template <typename _Int, typename _Flt>
void CHMatrix<_Int,_Flt>::SymmetrizeSubmatrices (void *_comm,
                                                   int _nblks, long long *_blks, int *_blk2cpu, 
                                                   CHMatrix<_Int,_Flt> *_hmatr_arr, CHMatrix<_Int,_Flt> *_hmatr_symm_arr) 
{

   int myid = CExchange::GetMyid (_comm);
   int nproc = CExchange::GetNproc (_comm);

// Create cpu reg data

   int icyclecpu = -1;

   vector<int> imaskcpu (nproc);
   vector<int> listcpu (nproc);
   vector<int> indcpu (nproc);

   int *pimaskcpu = &imaskcpu[0];
   int *plistcpu = &listcpu[0];
   int *pindcpu = &indcpu[0];

   int i;

   for (i=0;i<nproc;i++) {
      pimaskcpu[i] = icyclecpu;
   }

// For own blocks compute the list of cpus

   icyclecpu++;

   int nlistcpu = 0;

   int iblk, j, jj, jcpu;

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         CMatrix<int,float> *pHMatr_sub = _hmatr_arr[iblk].GetHMatrStr();
         int nlist_hmatr = pHMatr_sub->GetNlist();
         int *pia_hmatr = pHMatr_sub->GetIaArr();
         int *pja_hmatr = pHMatr_sub->GetJaArr();
         for (i=0;i<nlist_hmatr;i++) {
            for (j=pia_hmatr[i];j<pia_hmatr[i+1];j++) {
               jj = pja_hmatr[j];
               if (jj != iblk) {
                  jcpu = _blk2cpu[jj];
                  if (pimaskcpu[jcpu] != icyclecpu) {
                     plistcpu[nlistcpu] = jcpu;
                     nlistcpu++;
                     pimaskcpu[jcpu] = icyclecpu;
                  }
               }
            }
         }
      }
   }

   sort (plistcpu, plistcpu+nlistcpu);

   for (i=0;i<nlistcpu;i++) {
      jcpu = plistcpu[i];
      pindcpu[jcpu] = i;
   }

// Compute number of blocks send to each cpu (except diagonal one)

   vector<int> nzblk_cpu (nlistcpu+1);
   int *pnzblk_cpu = &nzblk_cpu[0];

   for (i=0;i<nlistcpu;i++) pnzblk_cpu[i] = 0;

   int ind;

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         CMatrix<int,float> *pHMatr_sub = _hmatr_arr[iblk].GetHMatrStr();
         int nlist_hmatr = pHMatr_sub->GetNlist();
         int *pia_hmatr = pHMatr_sub->GetIaArr();
         int *pja_hmatr = pHMatr_sub->GetJaArr();
         for (i=0;i<nlist_hmatr;i++) {
            for (j=pia_hmatr[i];j<pia_hmatr[i+1];j++) {
               jj = pja_hmatr[j];
               if (jj != iblk) {
                  jcpu = _blk2cpu[jj];
                  ind = pindcpu[jcpu];
                  pnzblk_cpu[ind]++;
               }
            }
         }
      }
   }

// Prepare send hblock data

   vector<CHMatrix<_Int,_Flt> > hblk_send (nlistcpu+1);

   CHMatrix<_Int,_Flt> *phblk_send = &hblk_send[0];

   for (i=0;i<nlistcpu;i++) {
      CMatrix<int,float> astr_temp;
      astr_temp.ResizeAndSetAllSp (0, 0, pnzblk_cpu[i], pnzblk_cpu[i]);
      CMatrix<int,float> *pHMatrStr = phblk_send[i].GetHMatrStr();
      pHMatrStr->ReplaceFree (astr_temp);
      phblk_send[i].SetNzblk (pnzblk_cpu[i]);
      phblk_send[i].ResizeASub (pnzblk_cpu[i]+1);
   }

   int nimax = 0;

   int niloc;

   for (i=0;i<_nblks;i++) {
      niloc = (int)(_blks[i+1]-_blks[i]);
      if (niloc > nimax) nimax = niloc;
   }

   vector<int> imaskblk (nimax+1);
   vector<int> listblk (nimax+1);
   vector<int> indblk (nimax+1);
   vector<int> iablk (nimax+1);
   vector<int> iptrblk (nimax+1);

   int *pimaskblk = &imaskblk[0];
   int *plistblk  = &listblk[0];
   int *pindblk   = &indblk[0];
   int *piablk    = &iablk[0];
   int *piptrblk  = &iptrblk[0];

   int icycleblk = -1;

   for (i=0;i<nimax;i++) pimaskblk[i] = icycleblk;

   for (i=0;i<nlistcpu;i++) pnzblk_cpu[i] = 0;

   int k;

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         CMatrix<_Int,_Flt> *pA_sub = _hmatr_arr[iblk].GetASubArr();
         CMatrix<int,float> *pHMatr_sub = _hmatr_arr[iblk].GetHMatrStr();
         int nlist_hmatr = pHMatr_sub->GetNlist();
         int *pia_hmatr = pHMatr_sub->GetIaArr();
         int *pja_hmatr = pHMatr_sub->GetJaArr();
         for (i=0;i<nlist_hmatr;i++) {
            for (j=pia_hmatr[i];j<pia_hmatr[i+1];j++) {
               jj = pja_hmatr[j];
               if (jj != iblk) {
                  jcpu = _blk2cpu[jj];
                  ind = pindcpu[jcpu];
                  k = pnzblk_cpu[ind];
                  CMatrix<_Int,_Flt> *pA_cpu = phblk_send[ind].GetASubArr();
                  CMatrix<int,float> *pHMatr_cpu = phblk_send[ind].GetHMatrStr();
                  int *pja_cpu = pHMatr_cpu->GetJaArr();
                  int *pja2_cpu = pHMatr_cpu->GetJa2Arr();
                  pja_cpu[k] = jj;
                  pja2_cpu[k] = iblk;
                  CMatrix<_Int,_Flt> a_sp;
                  a_sp.GetSparsity (pA_sub[j]);
                  a_sp.TransposedSparsityListSp (icycleblk, pimaskblk, pindblk, 
                                                   piptrblk, plistblk, piablk, pA_cpu[k]);
                  pnzblk_cpu[ind]++;
               }
            }
         }
      }
   }

// Pack send data

   vector<int> CpuIDSend (nlistcpu);
   vector<vector<char> > ObjSend (nlistcpu);

   int *pCpuIDSend = NULL;
   vector<char> *pObjSend = NULL;

   if (nlistcpu > 0) {
      pCpuIDSend = &CpuIDSend[0];
      pObjSend = &ObjSend[0];
   }

   long long isize;
   char *pobj;

   for (i=0;i<nlistcpu;i++) {
      pCpuIDSend[i] = plistcpu[i];
      isize = phblk_send[i].GetPackedSize();
      pObjSend[i].resize ((size_t)isize);
      pobj = &(pObjSend[i][0]);
      phblk_send[i].FillPacked (isize, pobj);
      phblk_send[i].Clean ();
   }

// Exchange

   vector<int> CpuIDRecv;
   vector<vector<char> > ObjRecv;

   CExchange::DataExchange (_comm, CpuIDSend, ObjSend, CpuIDRecv, ObjRecv);

// Free send data

   {
      vector<int> CpuIDSend_temp;
      vector<vector<char> > ObjSend_temp;
      CpuIDSend.swap (CpuIDSend_temp);
      ObjSend.swap (ObjSend_temp);
   }

// Unpack receive data

   int nrecv = (int) CpuIDRecv.size();

   vector<char> *pObjRecv = NULL;

   if (nrecv > 0) {
      pObjRecv = &ObjRecv[0];
   }

   vector<CHMatrix<_Int,_Flt> > hblk_recv (nrecv+1);

   CHMatrix<_Int,_Flt> *phblk_recv = &hblk_recv[0];

   for (i=0;i<nrecv;i++) {
      isize = (long long) pObjRecv[i].size();
      pobj = &(pObjRecv[i][0]);
      phblk_recv[i].UnPack (isize, pobj);
   }

// Free recv data

   {
      vector<int> CpuIDRecv_temp;
      vector<vector<char> > ObjRecv_temp;
      CpuIDRecv.swap (CpuIDRecv_temp);
      ObjRecv.swap (ObjRecv_temp);
   }

// Store received block sparsity for each own block row

   vector<int> iablk_recv (_nblks+1);
   int *piablk_recv = &iablk_recv[0];

   for (i=0;i<=_nblks;i++) piablk_recv[i] = 0;

   for (i=0;i<nrecv;i++) {
      CMatrix<int,float> *pHMatr_sub = phblk_recv[i].GetHMatrStr();
      int nzja_hmatr = pHMatr_sub->GetNzja();
      int *pja_hmatr = pHMatr_sub->GetJaArr();
      for (j=0;j<nzja_hmatr;j++) {
         iblk = pja_hmatr[j];
         piablk_recv[iblk+1]++;
      }
   }

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         piablk_recv[iblk+1]++;
      }
   }

   for (i=0;i<_nblks;i++) piablk_recv[i+1] = piablk_recv[i]+piablk_recv[i+1];

   int nzblk_recv = piablk_recv[_nblks];

   vector<int> iptr_recv (_nblks);
   vector<int> jablk_recv (nzblk_recv+1);
   vector<CMatrix<_Int,_Flt> > ablk_recv (nzblk_recv+1);

   int *piptr_recv = &iptr_recv[0];
   int *pjablk_recv = &jablk_recv[0];
   CMatrix<_Int,_Flt> *pablk_recv = &ablk_recv[0];

   for (i=0;i<_nblks;i++) piptr_recv[i] = piablk_recv[i];

   int jblk;

   for (i=0;i<nrecv;i++) {
      CMatrix<_Int,_Flt> *pA_sub = phblk_recv[i].GetASubArr();
      CMatrix<int,float> *pHMatr_sub = phblk_recv[i].GetHMatrStr();
      int nzja_hmatr = pHMatr_sub->GetNzja();
      int *pja_hmatr = pHMatr_sub->GetJaArr();
      int *pja2_hmatr = pHMatr_sub->GetJa2Arr();
      for (j=0;j<nzja_hmatr;j++) {
         iblk = pja_hmatr[j];
         jblk = pja2_hmatr[j];
         k = piptr_recv[iblk];
         pjablk_recv[k] = jblk;
         pablk_recv[k].ReplaceFree (pA_sub[j]);
         piptr_recv[iblk]++;
      }
   }

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         k = piptr_recv[iblk];
         pjablk_recv[k] = iblk;
         CMatrix<_Int,_Flt> *pA_sub = _hmatr_arr[iblk].GetASubArr();
         CMatrix<int,float> *pHMatr_sub = _hmatr_arr[iblk].GetHMatrStr();
         int nzja_hmatr = pHMatr_sub->GetNzja();
         int *pja_hmatr = pHMatr_sub->GetJaArr();
         ind = -1;
         for (j=0;j<nzja_hmatr;j++) {
            jblk = pja_hmatr[j];
            if (jblk == iblk) ind = j;
         }
         if (ind < 0) {
            throw " CHMatrix<>::SymmetrizeSubmatrices: error: incorrect block sparsity ";
         }
         CMatrix<_Int,_Flt> a_sp;
         a_sp.GetSparsity (pA_sub[ind]);
         a_sp.TransposedSparsityListSp (icycleblk, pimaskblk, pindblk, 
                                          piptrblk, plistblk, piablk, pablk_recv[k]);
         piptr_recv[iblk]++;
      }
   }

// Free received hblock data

   for (i=0;i<nrecv;i++) {
      phblk_recv[i].Clean ();
   }

// Compute symmetrized blocks

   vector<CSortInt> iiarr (_nblks);
   CSortInt *piiarr = &iiarr[0];

   int nlist_t, ibeg, nlist_new, ip, ip_t, jj_t;

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         nlist_t = piablk_recv[iblk+1]-piablk_recv[iblk];
         ibeg = piablk_recv[iblk];
         for (j=0;j<nlist_t;j++) {
            piiarr[j].ival = pjablk_recv[j+ibeg];
            piiarr[j].i2val = j+ibeg;
         }
         sort (piiarr,piiarr+nlist_t);
         CMatrix<_Int,_Flt> *pA_sub = _hmatr_arr[iblk].GetASubArr();
         CMatrix<int,float> *pHMatr_sub = _hmatr_arr[iblk].GetHMatrStr();
         int nzja_hmatr = pHMatr_sub->GetNzja();
         int *pja_hmatr = pHMatr_sub->GetJaArr();
         nlist_new = 0;
         ip = 0;
         ip_t = 0;
         while (ip < nzja_hmatr || ip_t < nlist_t) {
            if (ip < nzja_hmatr && ip_t < nlist_t) {
               jj = pja_hmatr[ip];
               jj_t = piiarr[ip_t].ival;
               if (jj == jj_t) {
                  nlist_new++;
                  ip++;
                  ip_t++;
               } else if (jj < jj_t) {
                  nlist_new++;
                  ip++;
               } else {
                  nlist_new++;
                  ip_t++;
               }
            } else if (ip < nzja_hmatr) {
               nlist_new++;
               ip++;
            } else {
               nlist_new++;
               ip_t++;
            }
         }
         _hmatr_symm_arr[iblk].ResizeASub (nlist_new+1);
         _hmatr_symm_arr[iblk].SetNzblk (nlist_new);
         CMatrix<_Int,_Flt> *pASymm_sub = _hmatr_symm_arr[iblk].GetASubArr();
         CMatrix<int,float> *pHMatrSymm_sub = _hmatr_symm_arr[iblk].GetHMatrStr();
         pHMatrSymm_sub->ResizeAndSetAllSp (1, 0, nlist_new, 0);
         int *plist_hmatrsymm = pHMatrSymm_sub->GetListArr();
         int *pia_hmatrsymm = pHMatrSymm_sub->GetIaArr();
         int *pja_hmatrsymm = pHMatrSymm_sub->GetJaArr();
         plist_hmatrsymm[0] = iblk;
         pia_hmatrsymm[0] = 0;
         pia_hmatrsymm[1] = nlist_new;
         nlist_new = 0;
         ip = 0;
         ip_t = 0;
         while (ip < nzja_hmatr || ip_t < nlist_t) {
            if (ip < nzja_hmatr && ip_t < nlist_t) {
               jj = pja_hmatr[ip];
               jj_t = piiarr[ip_t].ival;
               if (jj == jj_t) {
                  pja_hmatrsymm[nlist_new] = jj;
                  ind = piiarr[ip_t].i2val;
                  pA_sub[ip].ExtendSparsity (pablk_recv[ind], pASymm_sub[nlist_new]);
                  nlist_new++;
                  ip++;
                  ip_t++;
               } else if (jj < jj_t) {
                  pja_hmatrsymm[nlist_new] = jj;
                  pASymm_sub[nlist_new] = pA_sub[ip];
                  nlist_new++;
                  ip++;
               } else {
                  pja_hmatrsymm[nlist_new] = jj_t;
                  ind = piiarr[ip_t].i2val;
                  CMatrix<_Int,_Flt> atemp;
                  atemp.ExtendSparsity (pablk_recv[ind], pASymm_sub[nlist_new]);
                  nlist_new++;
                  ip_t++;
               }
            } else if (ip < nzja_hmatr) {
               pja_hmatrsymm[nlist_new] = pja_hmatr[ip];
               pASymm_sub[nlist_new] = pA_sub[ip];
               nlist_new++;
               ip++;
            } else {
               pja_hmatrsymm[nlist_new] = piiarr[ip_t].ival;
               ind = piiarr[ip_t].i2val;
               CMatrix<_Int,_Flt> atemp;
               atemp.ExtendSparsity (pablk_recv[ind], pASymm_sub[nlist_new]);
               nlist_new++;
               ip_t++;
            }
         }
      }
   }

}

// Compute the extended lists
//========================================================================================
template <typename _Int, typename _Flt>
void CHMatrix<_Int,_Flt>::ExtendedLists (void *_comm, int _ncycle,
                                                int _nblks, long long *_blks, int *_blk2cpu, 
                                                CHMatrix<_Int,_Flt> *_hmatr_arr,
                                                int *_nlist_ext_arr, vector<int> *_list_ext_arr) {

   int myid = CExchange::GetMyid (_comm);
   int nproc = CExchange::GetNproc (_comm);

// Compute the maximal block size

   int nimax = 0;

   int i;
   int niloc;

   for (i=0;i<_nblks;i++) {
      niloc = (int)(_blks[i+1]-_blks[i]);
      if (niloc > nimax) nimax = niloc;
   }

   vector<int> imaskblk (nimax+1);
   vector<int> listblk (nimax+1);
   vector<int> iablk (nimax+1);

   int *pimaskblk = &imaskblk[0];
   int *plistblk = &listblk[0];
   int *piablk = &iablk[0];

   int icycleblk = -1;

   for (i=0;i<nimax;i++) pimaskblk[i] = icycleblk;

// Compute output data for simple cases

   for (i=0;i<_nblks;i++) _nlist_ext_arr[i] = 0;

   int iblk, j, jj, k, kk, nlistnew;

   if (_ncycle == 0) {
      return;
   } else if (_ncycle == 1) {

      for (iblk=0;iblk<_nblks;iblk++) {
         if (_blk2cpu[iblk] == myid) {
            CMatrix<_Int,_Flt> *pA_sub = _hmatr_arr[iblk].GetASubArr();
            CMatrix<int,float> *pHMatr_sub = _hmatr_arr[iblk].GetHMatrStr();
            int nzja_hmatr = pHMatr_sub->GetNzja();
            int *pja_hmatr = pHMatr_sub->GetJaArr();
            nlistnew = 0;
            for (j=0;j<nzja_hmatr;j++) {
               jj = pja_hmatr[j];
               if (jj != iblk) {
                  icycleblk++;
                  int nzja_sub = pA_sub[j].GetNzja();
                  _Int *pja_sub = pA_sub[j].GetJaArr();
                  for (k=0;k<nzja_sub;k++) {
                     kk = (int)pja_sub[k];
                     if (pimaskblk[kk] != icycleblk) {
                        nlistnew++;
                        pimaskblk[kk] = icycleblk;
                     }
                  }
               }
            }
            _nlist_ext_arr[iblk] = nlistnew;
            _list_ext_arr[iblk].resize (nlistnew*2+1);
            int *p_list_ext_arr = &(_list_ext_arr[iblk][0]);
            nlistnew = 0;
            for (j=0;j<nzja_hmatr;j++) {
               jj = pja_hmatr[j];
               if (jj != iblk) {
                  icycleblk++;
                  int nzja_sub = pA_sub[j].GetNzja();
                  _Int *pja_sub = pA_sub[j].GetJaArr();
                  int nlist_temp = 0;
                  for (k=0;k<nzja_sub;k++) {
                     kk = (int)pja_sub[k];
                     if (pimaskblk[kk] != icycleblk) {
                        plistblk[nlist_temp] = kk;
                        nlist_temp++;
                        pimaskblk[kk] = icycleblk;
                     }
                  }
                  sort (plistblk,plistblk+nlist_temp);
                  for (k=0;k<nlist_temp;k++) {
                     p_list_ext_arr[nlistnew*2] = plistblk[k];
                     p_list_ext_arr[nlistnew*2+1] = jj;
                     nlistnew++;
                  }
               }
            }
         }
      }
      return;
   }

// Compute data as two index rows block structures

   vector<CMatrix<_Int,_Flt> > blkrows (_nblks);
   CMatrix<_Int,_Flt> *pblkrows = &blkrows[0];

   int nzjanew, ind, kj;

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         CMatrix<_Int,_Flt> *pA_sub = _hmatr_arr[iblk].GetASubArr();
         CMatrix<int,float> *pHMatr_sub = _hmatr_arr[iblk].GetHMatrStr();
         int nzja_hmatr = pHMatr_sub->GetNzja();
         int *pja_hmatr = pHMatr_sub->GetJaArr();
         nzjanew = 0;
         for (j=0;j<nzja_hmatr;j++) {
            nzjanew += pA_sub[j].GetNzja();
         }
         nlistnew = (int)(_blks[iblk+1]-_blks[iblk]);
         CMatrix<_Int,_Flt> anew;
         anew.ResizeAndSetAllSp (nlistnew,nlistnew,nzjanew,nzjanew);
         _Int *plist_new = anew.GetListArr();
         _Int *plist2_new = anew.GetList2Arr();
         _Int *pia_new = anew.GetIaArr();
         _Int *pja_new = anew.GetJaArr();
         _Int *pja2_new = anew.GetJa2Arr();
         for (j=0;j<nlistnew;j++) plist_new[j] = (_Int)j;
         for (j=0;j<nlistnew;j++) plist2_new[j] = (_Int)iblk;
         for (j=0;j<=nlistnew;j++) pia_new[j] = 0;
         for (j=0;j<nzja_hmatr;j++) {
            int nlistloc = pA_sub[j].GetNlist();
            _Int *plistloc = pA_sub[j].GetListArr();
            _Int *pialoc = pA_sub[j].GetIaArr();
            for (k=0;k<nlistloc;k++) {
               kk = (int)plistloc[k];
               pia_new[kk+1] += (pialoc[k+1]-pialoc[k]);
            }
         }
         for (j=0;j<nlistnew;j++) pia_new[j+1] = pia_new[j]+pia_new[j+1];
         for (j=0;j<nlistnew;j++) piablk[j] = (int)pia_new[j];
         for (j=0;j<nzja_hmatr;j++) {
            jj = (int)pja_hmatr[j];
            int nlistloc = pA_sub[j].GetNlist();
            _Int *plistloc = pA_sub[j].GetListArr();
            _Int *pialoc = pA_sub[j].GetIaArr();
            _Int *pjaloc = pA_sub[j].GetJaArr();
            for (k=0;k<nlistloc;k++) {
               kk = (int)plistloc[k];
               ind = (int)piablk[kk];
               for (kj=(int)pialoc[k];kj<pialoc[k+1];kj++) {
                  pja_new[ind] = pjaloc[kj];
                  pja2_new[ind] = jj;
                  ind++;
               }
               piablk[kk] = ind;
            }
         }
         pblkrows[iblk].ReplaceFree (anew);
      }
   }

// Prepare initial lists data

   vector<int> nlist_arr (_nblks);
   vector<vector<int> > list_arr (_nblks);
   vector<int> nlist_arr_prev (_nblks);
   vector<vector<int> > list_arr_prev (_nblks);

   int *pnlist_arr_prev = &nlist_arr_prev[0];
   int *pnlist_arr = &nlist_arr[0];
   vector<int> *plist_arr = &list_arr[0];
   vector<int> *plist_arr_prev = &list_arr_prev[0];

   for (iblk=0;iblk<_nblks;iblk++) pnlist_arr_prev[iblk] = 0;
   for (iblk=0;iblk<_nblks;iblk++) pnlist_arr[iblk] = 0;

   int jj2;

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         nlistnew = pblkrows[iblk].GetNlist();
         nzjanew = pblkrows[iblk].GetNzja();
         _Int *pia_row = pblkrows[iblk].GetIaArr();
         _Int *pja_row = pblkrows[iblk].GetJaArr();
         _Int *pja2_row = pblkrows[iblk].GetJa2Arr();
         int nzjaflt = 0;
         for (j=0;j<nzjanew;j++) {
            jj2 = (int)pja2_row[j];
            if (jj2 != iblk) nzjaflt++;
         }
         vector<CSortInt2> ii2arr (nzjaflt+1);
         CSortInt2 *pii2arr = &ii2arr[0];
         nzjaflt = 0;
         for (j=0;j<nzjanew;j++) {
            jj = (int)pja_row[j];
            jj2 = (int)pja2_row[j];
            if (jj2 != iblk) {
               pii2arr[nzjaflt].ixval = jj2;
               pii2arr[nzjaflt].iyval = jj;
               pii2arr[nzjaflt].itail = j;
               nzjaflt++;
            }
         }
         sort (pii2arr,pii2arr+nzjaflt);
         int nlist_temp = 0;
         int ix_prev = -1;
         int iy_prev = -1;
         int ix_curr, iy_curr, ival_curr;
         for (j=0;j<nzjaflt;j++) {
            ix_curr = pii2arr[j].ixval;
            iy_curr = pii2arr[j].iyval;
            ival_curr = pii2arr[j].itail;
            if (ix_curr != ix_prev || iy_curr != iy_prev) {
               pii2arr[nlist_temp].ixval = ix_curr;
               pii2arr[nlist_temp].iyval = iy_curr;
               pii2arr[nlist_temp].itail = ival_curr;
               nlist_temp++;
               ix_prev = ix_curr;
               iy_prev = iy_curr;
            }
         }
         plist_arr[iblk].resize(2*nlist_temp+1);
         plist_arr_prev[iblk].resize(2*nlist_temp+1);
         int *pplist_arr = &(plist_arr[iblk][0]);
         int *pplist_arr_prev = &(plist_arr_prev[iblk][0]);
         for (j=0;j<nlist_temp;j++) {
            pplist_arr[j*2] = pii2arr[j].iyval;
            pplist_arr[j*2+1] = pii2arr[j].ixval;
            pplist_arr_prev[j*2] = pii2arr[j].iyval;
            pplist_arr_prev[j*2+1] = pii2arr[j].ixval;
         }
         pnlist_arr[iblk] = nlist_temp;
         pnlist_arr_prev[iblk] = nlist_temp;
      }
   }

// Create mask data

   int icyclecpu = -1;

   vector<int> imaskcpu (nproc);
   vector<int> listcpu (nproc);
   vector<int> indcpu (nproc);

   int *pimaskcpu = &imaskcpu[0];
   int *plistcpu = &listcpu[0];
   int *pindcpu = &indcpu[0];

   for (i=0;i<nproc;i++) {
      pimaskcpu[i] = icyclecpu;
   }

   int icycletot = -1;

   vector<int> imasktot (_nblks);
   vector<int> listtot (_nblks);
   vector<int> indtot (_nblks);

   int *pimasktot = &imasktot[0];
   int *plisttot = &listtot[0];
   int *pindtot = &indtot[0];

   for (i=0;i<_nblks;i++) {
      pimasktot[i] = icycletot;
   }

// Main extention cycle

   int icycle_ext, jcpu;

   for (icycle_ext=0;icycle_ext<_ncycle-1;icycle_ext++) {

// Create list of cpus and count sizes

      icyclecpu++;

      int nlistcpu = 0;

      for (iblk=0;iblk<_nblks;iblk++) {
         if (_blk2cpu[iblk] == myid) {
            int *pplist_arr_prev = &(plist_arr_prev[iblk][0]);
            for (j=0;j<pnlist_arr_prev[iblk];j++) {
               jj2 = pplist_arr_prev[j*2+1];
               jcpu = _blk2cpu[jj2];
               if (pimaskcpu[jcpu] != icyclecpu) {
                  plistcpu[nlistcpu] = jcpu;
                  nlistcpu++;
                  pimaskcpu[jcpu] = icyclecpu;
               }
            }
         }
      }

      sort (plistcpu,plistcpu+nlistcpu);

      for (i=0;i<nlistcpu;i++) {
         jcpu = plistcpu[i];
         pindcpu[jcpu] = i;
      }

      vector<int> nz_send (nlistcpu+1);

      int *pnz_send = &nz_send[0];

      for (i=0;i<nlistcpu;i++) pnz_send[i] = 0;

      for (iblk=0;iblk<_nblks;iblk++) {
         if (_blk2cpu[iblk] == myid) {
            int *pplist_arr_prev = &(plist_arr_prev[iblk][0]);
            for (j=0;j<pnlist_arr_prev[iblk];j++) {
               jj2 = pplist_arr_prev[j*2+1];
               jcpu = _blk2cpu[jj2];
               ind = pindcpu[jcpu];
               pnz_send[ind]++;
            }
         }
      }

// Prepare send data

      vector<CHMatrix<_Int,_Flt> > hblk_send (nlistcpu+1);

      CHMatrix<_Int,_Flt> *phblk_send = &hblk_send[0];

      for (i=0;i<nlistcpu;i++) {
         phblk_send[i].SetNzblk (1);
         phblk_send[i].ResizeASub (1);
         CMatrix<_Int,_Flt> *pA_sub = phblk_send[i].GetASubArr();
         CMatrix<_Int,_Flt> ablk_temp;
         ablk_temp.ResizeAndSetAllSp (0, 0, 3*pnz_send[i], 0);
         pA_sub->ReplaceFree (ablk_temp);
      }

      for (i=0;i<nlistcpu;i++) pnz_send[i] = 0;

      for (iblk=0;iblk<_nblks;iblk++) {
         if (_blk2cpu[iblk] == myid) {
            int *pplist_arr_prev = &(plist_arr_prev[iblk][0]);
            for (j=0;j<pnlist_arr_prev[iblk];j++) {
               jj = pplist_arr_prev[j*2];
               jj2 = pplist_arr_prev[j*2+1];
               jcpu = _blk2cpu[jj2];
               ind = pindcpu[jcpu];
               CMatrix<_Int,_Flt> *pA_sub = phblk_send[ind].GetASubArr();
               _Int *pja_asub = pA_sub->GetJaArr();
               k = pnz_send[ind];
               pja_asub[3*k] = iblk;
               pja_asub[3*k+1] = jj;
               pja_asub[3*k+2] = jj2;
               pnz_send[ind]++;
            }
         }
      }

// Pack send data

      vector<int> CpuIDSend (nlistcpu);
      vector<vector<char> > ObjSend (nlistcpu);

      int *pCpuIDSend = NULL;
      vector<char> *pObjSend = NULL;

      if (nlistcpu > 0) {
         pCpuIDSend = &CpuIDSend[0];
         pObjSend = &ObjSend[0];
      }

      long long isize;
      char *pobj;

      for (i=0;i<nlistcpu;i++) {
         pCpuIDSend[i] = plistcpu[i];
         isize = phblk_send[i].GetPackedSize();
         pObjSend[i].resize ((size_t)isize);
         pobj = &(pObjSend[i][0]);
         phblk_send[i].FillPacked (isize, pobj);
         phblk_send[i].Clean ();
      }

// Exchange

      vector<int> CpuIDRecv;
      vector<vector<char> > ObjRecv;

      CExchange::DataExchange (_comm, CpuIDSend, ObjSend, CpuIDRecv, ObjRecv);

// Free send data

      {
         vector<int> CpuIDSend_temp;
         vector<vector<char> > ObjSend_temp;
         CpuIDSend.swap (CpuIDSend_temp);
         ObjSend.swap (ObjSend_temp);
      }

// Unpack receive data

      int nrecv = (int) CpuIDRecv.size();

      vector<char> *pObjRecv = NULL;

      if (nrecv > 0) {
         pObjRecv = &ObjRecv[0];
      }

      vector<CHMatrix<_Int,_Flt> > hblk_recv (nrecv+1);

      CHMatrix<_Int,_Flt> *phblk_recv = &hblk_recv[0];

      for (i=0;i<nrecv;i++) {
         isize = (long long) pObjRecv[i].size();
         pobj = &(pObjRecv[i][0]);
         phblk_recv[i].UnPack (isize, pobj);
      }

// Free recv data

      {
         vector<vector<char> > ObjRecv_temp;
         ObjRecv.swap (ObjRecv_temp);
      }

// Prepare the answer for all received data

      for (i=0;i<nrecv;i++) {
         CMatrix<_Int,_Flt> *pA_sub = phblk_recv[i].GetASubArr();
         int nzjaloc = pA_sub->GetNzja();
         _Int *pja_sub = pA_sub->GetJaArr();
         int nlistloc = nzjaloc / 3;
         int nzja_new = 0;
         for (j=0;j<nlistloc;j++) {
            jj = (int)pja_sub[j*3+1];
            jj2 = (int)pja_sub[j*3+2];
            _Int *pia_rows = pblkrows[jj2].GetIaArr();
            nzja_new += (int)(pia_rows[jj+1]-pia_rows[jj]);
         }
         CMatrix<_Int,_Flt> a_new;
         a_new.ResizeAndSetAllSp (nlistloc, 0, nzja_new, nzja_new);
         _Int *plist_new = a_new.GetListArr();
         _Int *pia_new = a_new.GetIaArr();
         _Int *pja_new = a_new.GetJaArr();
         _Int *pja2_new = a_new.GetJa2Arr();
         nzja_new = 0;
         pia_new[0] = 0;
         for (j=0;j<nlistloc;j++) {
            plist_new[j] = pja_sub[j*3];
            jj = (int)pja_sub[j*3+1];
            jj2 = (int)pja_sub[j*3+2];
            _Int *pia_rows = pblkrows[jj2].GetIaArr();
            _Int *pja_rows = pblkrows[jj2].GetJaArr();
            _Int *pja2_rows = pblkrows[jj2].GetJa2Arr();
            for (k=(int)pia_rows[jj];k<pia_rows[jj+1];k++) {
               pja_new[nzja_new] = pja_rows[k];
               pja2_new[nzja_new] = pja2_rows[k];
               nzja_new++;
            }
            pia_new[j+1] = nzja_new;
         }
         pA_sub->ReplaceFree (a_new);
      }

// Pack send data

      ObjRecv.resize (nrecv);

      pObjRecv = NULL;

      if (nrecv > 0) {
         pObjRecv = &ObjRecv[0];
      }

      for (i=0;i<nrecv;i++) {
         isize = phblk_recv[i].GetPackedSize();
         pObjRecv[i].resize ((size_t)isize);
         pobj = &(pObjRecv[i][0]);
         phblk_recv[i].FillPacked (isize, pobj);
         phblk_recv[i].Clean ();
      }

// Exchange data back

      CExchange::DataExchange (_comm, CpuIDRecv, ObjRecv, CpuIDSend, ObjSend);

// Free send data

      {
         vector<int> CpuIDSend_temp;
         vector<vector<char> > ObjSend_temp;
         CpuIDRecv.swap (CpuIDSend_temp);
         ObjRecv.swap (ObjSend_temp);
      }

// Unpack receive data

      nrecv = (int) CpuIDSend.size();

      pObjSend = NULL;
      if (nrecv > 0) {
         pObjSend = &ObjSend[0];
      }

      hblk_send.resize (nrecv+1);

      phblk_send = &hblk_send[0];

      for (i=0;i<nrecv;i++) {
         isize = (long long) pObjSend[i].size();
         pobj = &(pObjSend[i][0]);
         phblk_send[i].UnPack (isize, pobj);
      }

// Free recv data

      {
         vector<int> CpuIDRecv_temp;
         vector<vector<char> > ObjRecv_temp;
         CpuIDSend.swap (CpuIDRecv_temp);
         ObjSend.swap (ObjRecv_temp);
      }

// Combine received data into the new lists

      icycletot++;

      int nlisttot = 0;

      for (i=0;i<nrecv;i++) {
         CMatrix<_Int,_Flt> *pA_sub = phblk_send[i].GetASubArr();
         int nlist_sub = pA_sub->GetNlist();
         _Int *plist_sub = pA_sub->GetListArr();
         for (j=0;j<nlist_sub;j++) {
            jj2 = (int)plist_sub[j];
            if (pimasktot[jj2] != icycletot) {
               plisttot[nlisttot] = jj2;
               nlisttot++;
               pimasktot[jj2] = icycletot;
            }
         }
      }

      sort (plisttot,plisttot+nlisttot);

      for (i=0;i<nlisttot;i++) {
         jj2 = plisttot[i];
         pindtot[jj2] = i;
      }

      vector<int> nz_list_new (nlisttot+1);
      int *pnz_list_new = &nz_list_new[0];

      for (i=0;i<nlisttot;i++) pnz_list_new[i] = 0;

      for (i=0;i<nrecv;i++) {
         CMatrix<_Int,_Flt> *pA_sub = phblk_send[i].GetASubArr();
         int nlist_sub = pA_sub->GetNlist();
         _Int *plist_sub = pA_sub->GetListArr();
         _Int *pia_sub = pA_sub->GetIaArr();
         for (j=0;j<nlist_sub;j++) {
            jj2 = (int)plist_sub[j];
            ind = (int)pindtot[jj2];
            pnz_list_new[ind] += (int)(pia_sub[j+1]-pia_sub[j]);
         }
      }

      vector<vector<int> > list_new (nlisttot+1);
      vector<int> *plist_new = &list_new[0];

      for (i=0;i<nlisttot;i++) {
         plist_new[i].resize(2*pnz_list_new[i]+1);
      }

      for (i=0;i<nlisttot;i++) pnz_list_new[i] = 0;

      int kk2;

      for (i=0;i<nrecv;i++) {
         CMatrix<_Int,_Flt> *pA_sub = phblk_send[i].GetASubArr();
         int nlist_sub = pA_sub->GetNlist();
         _Int *plist_sub = pA_sub->GetListArr();
         _Int *pia_sub = pA_sub->GetIaArr();
         _Int *pja_sub = pA_sub->GetJaArr();
         _Int *pja2_sub = pA_sub->GetJa2Arr();
         for (j=0;j<nlist_sub;j++) {
            jj2 = (int)plist_sub[j];
            ind = (int)pindtot[jj2];
            k = (int)pnz_list_new[ind];
            int *pplist_new = &(plist_new[ind][0]);
            for (kj=(int)pia_sub[j];kj<pia_sub[j+1];kj++) {
               kk = (int)pja_sub[kj];
               kk2 = (int)pja2_sub[kj];
               pplist_new[k*2] = kk;
               pplist_new[k*2+1] = kk2;
               k++;
            }
            pnz_list_new[ind] = k;
         }
      }

// Sort and filter

      for (i=0;i<nlisttot;i++) {

// Sort and filter coincident and own block numbers

         iblk = plisttot[i];
         int *pplist_new = &(plist_new[i][0]);
         vector<CSortInt2> ii2arr (pnz_list_new[i]+1);
         CSortInt2 *pii2arr = &ii2arr[0];
         for (j=0;j<pnz_list_new[i];j++) {
            pii2arr[j].ixval = pplist_new[j*2+1];
            pii2arr[j].iyval = pplist_new[j*2];
            pii2arr[j].itail = j;
         }
         sort (pii2arr,pii2arr+pnz_list_new[i]);
         int nznew = 0;
         int ix_prev = -1;
         int iy_prev = -1;
         int ix_curr, iy_curr;
         for (j=0;j<pnz_list_new[i];j++) {
            ix_curr = pii2arr[j].ixval;
            iy_curr = pii2arr[j].iyval;
            if (ix_curr != iblk && (ix_curr != ix_prev || iy_curr != iy_prev)) {
               pii2arr[nznew].ixval = ix_curr;
               pii2arr[nznew].iyval = iy_curr;
               pii2arr[nznew].itail = nznew;
               nznew++;
               ix_prev = ix_curr;
               iy_prev = iy_curr;
            }
         }

// Filter by previous data and extend previous sorted list data

         int *pplist_arr = &(plist_arr[iblk][0]);

         vector<int> list_add (2*(pnlist_arr[iblk]+nznew)+1);
         vector<int> list_flt (2*nznew+1);

         int *plist_add = &list_add[0];
         int *plist_flt = &list_flt[0];

         int nlistadd = 0;
         int nlistflt = 0;

         int jj_new, jj2_new;

         int ip = 0;
         int ip_new = 0;

         while (ip < pnlist_arr[iblk] || ip_new < nznew) {
            if (ip < pnlist_arr[iblk] && ip_new < nznew) {
               jj = pplist_arr[2*ip];
               jj2 = pplist_arr[2*ip+1];
               jj_new = pii2arr[ip_new].iyval;
               jj2_new = pii2arr[ip_new].ixval;
               if (jj == jj_new && jj2 == jj2_new) {
                  plist_add[2*nlistadd] = jj;
                  plist_add[2*nlistadd+1] = jj2;
                  nlistadd++;
                  ip++;
                  ip_new++;
               } else if (jj2 < jj2_new || (jj2 == jj2_new && jj < jj_new)) {
                  plist_add[2*nlistadd] = jj;
                  plist_add[2*nlistadd+1] = jj2;
                  nlistadd++;
                  ip++;
               } else {
                  plist_add[2*nlistadd] = jj_new;
                  plist_add[2*nlistadd+1] = jj2_new;
                  nlistadd++;
                  plist_flt[2*nlistflt] = jj_new;
                  plist_flt[2*nlistflt+1] = jj2_new;
                  nlistflt++;
                  ip_new++;
               }
            } else if (ip < pnlist_arr[iblk]) {
               plist_add[2*nlistadd] = pplist_arr[2*ip];
               plist_add[2*nlistadd+1] = pplist_arr[2*ip+1];
               nlistadd++;
               ip++;
            } else {
               plist_add[2*nlistadd] = pii2arr[ip_new].iyval;
               plist_add[2*nlistadd+1] = pii2arr[ip_new].ixval;
               nlistadd++;
               plist_flt[2*nlistflt] = pii2arr[ip_new].iyval;
               plist_flt[2*nlistflt+1] = pii2arr[ip_new].ixval;
               nlistflt++;
               ip_new++;
            }
         }

// Replace

         pnlist_arr[iblk] = nlistadd;
         pnlist_arr_prev[iblk] = nlistflt;

         plist_arr[iblk].swap(list_add);
         plist_arr_prev[iblk].swap(list_flt);

      }

   }

// Prepare final data

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         int *pplist_arr = &(plist_arr[iblk][0]);
         _nlist_ext_arr[iblk] = pnlist_arr[iblk];
         _list_ext_arr[iblk].resize(2*pnlist_arr[iblk]+1);
         int *pplist_ext = &(_list_ext_arr[iblk][0]);
         for (j=0;j<2*(pnlist_arr[iblk]);j++) pplist_ext[j] = pplist_arr[j];
      }
   }

}

// Get the extended submatrices
//========================================================================================
template <typename _Int, typename _Flt>
void CHMatrix<_Int,_Flt>::GetExtendedSubmatrices (void *_comm,
                                                         int _nblks, long long *_blks, int *_blk2cpu, 
                                                         CHMatrix<_Int,_Flt> *_hmatr_arr,
                                                         int *_nlist_ext_arr, vector<int> *_list_ext_arr,
                                                         CMatrix<_Int,_Flt> *_matr_ext_arr)
{

   int myid = CExchange::GetMyid (_comm);
   int nproc = CExchange::GetNproc (_comm);

// Compute the maximal block size

   int nimax = 0;

   int i;
   int niloc;

   for (i=0;i<_nblks;i++) {
      niloc = (int)(_blks[i+1]-_blks[i]);
      if (niloc > nimax) nimax = niloc;
   }

   vector<int> imaskblk (nimax+1);
   vector<int> listblk (nimax+1);
   vector<int> iablk (nimax+1);

   int *pimaskblk = &imaskblk[0];
   int *plistblk = &listblk[0];
   int *piablk = &iablk[0];

   int icycleblk = -1;

   for (i=0;i<nimax;i++) pimaskblk[i] = icycleblk;

// Compute data as two index rows block structures

   vector<CMatrix<_Int,_Flt> > blkrows (_nblks);
   CMatrix<_Int,_Flt> *pblkrows = &blkrows[0];

   int nzjanew, ind, kj;
   int iblk, j, jj, k, kk, nlistnew;

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         CMatrix<_Int,_Flt> *pA_sub = _hmatr_arr[iblk].GetASubArr();
         CMatrix<int,float> *pHMatr_sub = _hmatr_arr[iblk].GetHMatrStr();
         int nzja_hmatr = pHMatr_sub->GetNzja();
         int *pja_hmatr = pHMatr_sub->GetJaArr();
         nzjanew = 0;
         for (j=0;j<nzja_hmatr;j++) {
            nzjanew += pA_sub[j].GetNzja();
         }
         nlistnew = (int)(_blks[iblk+1]-_blks[iblk]);
         CMatrix<_Int,_Flt> anew;
         anew.ResizeAndSetAll (nlistnew,nlistnew,nzjanew,nzjanew,nzjanew);
         _Int *plist_new = anew.GetListArr();
         _Int *plist2_new = anew.GetList2Arr();
         _Int *pia_new = anew.GetIaArr();
         _Int *pja_new = anew.GetJaArr();
         _Int *pja2_new = anew.GetJa2Arr();
         _Flt *pa_new = anew.GetAArr();
         for (j=0;j<nlistnew;j++) plist_new[j] = (_Int)j;
         for (j=0;j<nlistnew;j++) plist2_new[j] = (_Int)iblk;
         for (j=0;j<=nlistnew;j++) pia_new[j] = 0;
         for (j=0;j<nzja_hmatr;j++) {
            int nlistloc = pA_sub[j].GetNlist();
            _Int *plistloc = pA_sub[j].GetListArr();
            _Int *pialoc = pA_sub[j].GetIaArr();
            for (k=0;k<nlistloc;k++) {
               kk = (int)plistloc[k];
               pia_new[kk+1] += (pialoc[k+1]-pialoc[k]);
            }
         }
         for (j=0;j<nlistnew;j++) pia_new[j+1] = pia_new[j]+pia_new[j+1];
         for (j=0;j<nlistnew;j++) piablk[j] = (int)pia_new[j];
         for (j=0;j<nzja_hmatr;j++) {
            jj = pja_hmatr[j];
            int nlistloc = pA_sub[j].GetNlist();
            _Int *plistloc = pA_sub[j].GetListArr();
            _Int *pialoc = pA_sub[j].GetIaArr();
            _Int *pjaloc = pA_sub[j].GetJaArr();
            _Flt *paloc = pA_sub[j].GetAArr();
            for (k=0;k<nlistloc;k++) {
               kk = (int)plistloc[k];
               ind = (int)piablk[kk];
               for (kj=(int)pialoc[k];kj<pialoc[k+1];kj++) {
                  pja_new[ind] = pjaloc[kj];
                  pja2_new[ind] = jj;
                  pa_new[ind] = paloc[kj];
                  ind++;
               }
               piablk[kk] = ind;
            }
         }
         pblkrows[iblk].ReplaceFree (anew);
      }
   }

// Create mask data

   int icyclecpu = -1;

   vector<int> imaskcpu (nproc);
   vector<int> listcpu (nproc);
   vector<int> indcpu (nproc);

   int *pimaskcpu = &imaskcpu[0];
   int *plistcpu = &listcpu[0];
   int *pindcpu = &indcpu[0];

   for (i=0;i<nproc;i++) {
      pimaskcpu[i] = icyclecpu;
   }

   int icycletot = -1;

   vector<int> imasktot (_nblks);
   vector<int> listtot (_nblks);
   vector<int> indtot (_nblks);
   vector<int> indtot2 (_nblks);

   int *pimasktot = &imasktot[0];
   int *plisttot = &listtot[0];
   int *pindtot = &indtot[0];
   int *pindtot2 = &indtot2[0];

   for (i=0;i<_nblks;i++) {
      pimasktot[i] = icycletot;
   }

// Create list of cpus and count sizes

   icyclecpu++;

   int nlistcpu = 0;

   int jj2, jcpu;

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         int *pplist_ext_arr = &(_list_ext_arr[iblk][0]);
         for (j=0;j<_nlist_ext_arr[iblk];j++) {
            jj2 = (int)pplist_ext_arr[j*2+1];
            jcpu = _blk2cpu[jj2];
            if (pimaskcpu[jcpu] != icyclecpu) {
               plistcpu[nlistcpu] = jcpu;
               nlistcpu++;
               pimaskcpu[jcpu] = icyclecpu;
            }
         }
      }
   }

   sort (plistcpu,plistcpu+nlistcpu);

   for (i=0;i<nlistcpu;i++) {
      jcpu = plistcpu[i];
      pindcpu[jcpu] = i;
   }

   vector<int> nz_send (nlistcpu+1);

   int *pnz_send = &nz_send[0];

   for (i=0;i<nlistcpu;i++) pnz_send[i] = 0;

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         int *pplist_ext_arr = &(_list_ext_arr[iblk][0]);
         for (j=0;j<_nlist_ext_arr[iblk];j++) {
            jj2 = pplist_ext_arr[j*2+1];
            jcpu = _blk2cpu[jj2];
            ind = pindcpu[jcpu];
            pnz_send[ind]++;
         }
      }
   }

// Prepare send data

   vector<CHMatrix<_Int,_Flt> > hblk_send (nlistcpu+1);

   CHMatrix<_Int,_Flt> *phblk_send = &hblk_send[0];

   for (i=0;i<nlistcpu;i++) {
      phblk_send[i].SetNzblk (1);
      phblk_send[i].ResizeASub (1);
      CMatrix<_Int,_Flt> *pA_sub = phblk_send[i].GetASubArr();
      CMatrix<_Int,_Flt> ablk_temp;
      ablk_temp.ResizeAndSetAllSp (0, 0, 3*pnz_send[i], 0);
      pA_sub->ReplaceFree (ablk_temp);
   }

   for (i=0;i<nlistcpu;i++) pnz_send[i] = 0;

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         int *pplist_ext_arr = &(_list_ext_arr[iblk][0]);
         for (j=0;j<_nlist_ext_arr[iblk];j++) {
            jj = (int)pplist_ext_arr[j*2];
            jj2 = (int)pplist_ext_arr[j*2+1];
            jcpu = _blk2cpu[jj2];
            ind = pindcpu[jcpu];
            CMatrix<_Int,_Flt> *pA_sub = phblk_send[ind].GetASubArr();
            _Int *pja_asub = pA_sub->GetJaArr();
            k = pnz_send[ind];
            pja_asub[3*k] = iblk;
            pja_asub[3*k+1] = jj;
            pja_asub[3*k+2] = jj2;
            pnz_send[ind]++;
         }
      }
   }

// Pack send data

   vector<int> CpuIDSend (nlistcpu);
   vector<vector<char> > ObjSend (nlistcpu);

   int *pCpuIDSend = NULL;
   vector<char> *pObjSend = NULL;

   if (nlistcpu > 0) {
      pCpuIDSend = &CpuIDSend[0];
      pObjSend = &ObjSend[0];
   }

   long long isize;
   char *pobj;

   for (i=0;i<nlistcpu;i++) {
      pCpuIDSend[i] = plistcpu[i];
      isize = phblk_send[i].GetPackedSize();
      pObjSend[i].resize ((size_t)isize);
      pobj = &(pObjSend[i][0]);
      phblk_send[i].FillPacked (isize, pobj);
      phblk_send[i].Clean ();
   }

// Exchange

   vector<int> CpuIDRecv;
   vector<vector<char> > ObjRecv;

   CExchange::DataExchange (_comm, CpuIDSend, ObjSend, CpuIDRecv, ObjRecv);

// Free send data

   {
      vector<int> CpuIDSend_temp;
      vector<vector<char> > ObjSend_temp;
      CpuIDSend.swap (CpuIDSend_temp);
      ObjSend.swap (ObjSend_temp);
   }

// Unpack receive data

   int nrecv = (int) CpuIDRecv.size();

   vector<char> *pObjRecv = NULL;

   if (nrecv > 0) {
      pObjRecv = &ObjRecv[0];
   }

   vector<CHMatrix<_Int,_Flt> > hblk_recv (nrecv+1);

   CHMatrix<_Int,_Flt> *phblk_recv = &hblk_recv[0];

   for (i=0;i<nrecv;i++) {
      isize = (long long) pObjRecv[i].size();
      pobj = &(pObjRecv[i][0]);
      phblk_recv[i].UnPack (isize, pobj);
   }

// Free recv data

   {
      vector<vector<char> > ObjRecv_temp;
      ObjRecv.swap (ObjRecv_temp);
   }

// Prepare the answer for all received data

   for (i=0;i<nrecv;i++) {
      CMatrix<_Int,_Flt> *pA_sub = phblk_recv[i].GetASubArr();
      int nzjaloc = pA_sub->GetNzja();
      _Int *pja_sub = pA_sub->GetJaArr();
      int nlistloc = nzjaloc / 3;
      int nzja_new = 0;
      for (j=0;j<nlistloc;j++) {
         jj = (int)pja_sub[j*3+1];
         jj2 = (int)pja_sub[j*3+2];
         _Int *pia_rows = pblkrows[jj2].GetIaArr();
         nzja_new += (int)(pia_rows[jj+1]-pia_rows[jj]);
      }
      CMatrix<_Int,_Flt> a_new;
      a_new.ResizeAndSetAll (nlistloc, nlistloc*2, nzja_new, nzja_new, nzja_new);
      _Int *plist_new = a_new.GetListArr();
      _Int *plist2_new = a_new.GetList2Arr();
      _Int *pia_new = a_new.GetIaArr();
      _Int *pja_new = a_new.GetJaArr();
      _Int *pja2_new = a_new.GetJa2Arr();
      _Flt *pa_new = a_new.GetAArr();
      nzja_new = 0;
      pia_new[0] = 0;
      for (j=0;j<nlistloc;j++) {
         plist_new[j] = pja_sub[j*3];
         jj = (int)pja_sub[j*3+1];
         jj2 = (int)pja_sub[j*3+2];
         plist2_new[j*2] = jj;
         plist2_new[j*2+1] = jj2;
         _Int *pia_rows = pblkrows[jj2].GetIaArr();
         _Int *pja_rows = pblkrows[jj2].GetJaArr();
         _Int *pja2_rows = pblkrows[jj2].GetJa2Arr();
         _Flt *pa_rows = pblkrows[jj2].GetAArr();
         for (k=(int)pia_rows[jj];k<pia_rows[jj+1];k++) {
            pja_new[nzja_new] = pja_rows[k];
            pja2_new[nzja_new] = pja2_rows[k];
            pa_new[nzja_new] = pa_rows[k];
            nzja_new++;
         }
         pia_new[j+1] = nzja_new;
      }
      pA_sub->ReplaceFree (a_new);
   }

// Pack send data

   ObjRecv.resize (nrecv);

   pObjRecv = NULL;

   if (nrecv > 0) {
      pObjRecv = &ObjRecv[0];
   }

   for (i=0;i<nrecv;i++) {
      isize = phblk_recv[i].GetPackedSize();
      pObjRecv[i].resize ((size_t)isize);
      pobj = &(pObjRecv[i][0]);
      phblk_recv[i].FillPacked (isize, pobj);
      phblk_recv[i].Clean ();
   }

// Exchange data back

   CExchange::DataExchange (_comm, CpuIDRecv, ObjRecv, CpuIDSend, ObjSend);

// Free send data

   {
      vector<int> CpuIDSend_temp;
      vector<vector<char> > ObjSend_temp;
      CpuIDRecv.swap (CpuIDSend_temp);
      ObjRecv.swap (ObjSend_temp);
   }

// Unpack receive data

   nrecv = (int) CpuIDSend.size();

   pObjSend = NULL;

   if (nrecv > 0) {
      pObjSend = &ObjSend[0];
   }

   hblk_send.resize (nrecv+1);

   phblk_send = &hblk_send[0];

   for (i=0;i<nrecv;i++) {
      isize = (long long) pObjSend[i].size();
      pobj = &(pObjSend[i][0]);
      phblk_send[i].UnPack (isize, pobj);
   }

// Free recv data

   {
      vector<int> CpuIDRecv_temp;
      vector<vector<char> > ObjRecv_temp;
      CpuIDSend.swap (CpuIDRecv_temp);
      ObjSend.swap (ObjRecv_temp);
   }

// Combine received data into the new lists and add local data

   icycletot++;

   int nlisttot = 0;

   for (i=0;i<nrecv;i++) {
      CMatrix<_Int,_Flt> *pA_sub = phblk_send[i].GetASubArr();
      int nlist_sub = pA_sub->GetNlist();
      _Int *plist_sub = pA_sub->GetListArr();
      for (j=0;j<nlist_sub;j++) {
         jj2 = (int)plist_sub[j];
         if (pimasktot[jj2] != icycletot) {
            plisttot[nlisttot] = jj2;
            nlisttot++;
            pimasktot[jj2] = icycletot;
         }
      }
   }

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         if (pimasktot[iblk] != icycletot) {
            plisttot[nlisttot] = iblk;
            nlisttot++;
            pimasktot[iblk] = icycletot;
         }
      }
   }

   sort (plisttot,plisttot+nlisttot);

   for (i=0;i<nlisttot;i++) {
      jj2 = (int)plisttot[i];
      pindtot[jj2] = i;
   }

   vector<int> nlist_new (_nblks);
   vector<int> nz_list_new (_nblks);

   int *pnlist_new = &nlist_new[0];
   int *pnz_list_new = &nz_list_new[0];

   for (i=0;i<_nblks;i++) pnlist_new[i] = 0;
   for (i=0;i<_nblks;i++) pnz_list_new[i] = 0;

   for (i=0;i<nrecv;i++) {
      CMatrix<_Int,_Flt> *pA_sub = phblk_send[i].GetASubArr();
      int nlist_sub = pA_sub->GetNlist();
      _Int *plist_sub = pA_sub->GetListArr();
      _Int *pia_sub = pA_sub->GetIaArr();
      for (j=0;j<nlist_sub;j++) {
         jj2 = (int)plist_sub[j];
         pnlist_new[jj2]++;
         pnz_list_new[jj2] += (int)(pia_sub[j+1]-pia_sub[j]);
      }
   }

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         int nlistloc = pblkrows[iblk].GetNlist();
         int nzjaloc = pblkrows[iblk].GetNzja();
         pnlist_new[iblk] += nlistloc;
         pnz_list_new[iblk] += nzjaloc;
      }
   }

   vector<CMatrix<_Int,_Flt> > matr_arr_new (_nblks);
   CMatrix<_Int,_Flt> *pmatr_arr_new = &matr_arr_new[0];

   for (i=0;i<nlisttot;i++) {
      iblk = plisttot[i];
      pmatr_arr_new[iblk].ResizeAndSetAll (pnlist_new[iblk], pnlist_new[iblk], pnz_list_new[iblk], pnz_list_new[iblk], pnz_list_new[iblk]);
      _Int *pia_new = pmatr_arr_new[iblk].GetIaArr();
      pia_new[0] = 0;
      pnlist_new[iblk] = 0;
      pnz_list_new[iblk] = 0;
   }

// Fill received data

   int kk2;

   for (i=0;i<nrecv;i++) {
      CMatrix<_Int,_Flt> *pA_sub = phblk_send[i].GetASubArr();
      int nlist_sub = pA_sub->GetNlist();
      _Int *plist_sub = pA_sub->GetListArr();
      _Int *plist2_sub = pA_sub->GetList2Arr();
      _Int *pia_sub = pA_sub->GetIaArr();
      _Int *pja_sub = pA_sub->GetJaArr();
      _Int *pja2_sub = pA_sub->GetJa2Arr();
      _Flt *pa_sub = pA_sub->GetAArr();
      for (j=0;j<nlist_sub;j++) {
         jj2 = (int)plist_sub[j];
         _Int *plist_new = pmatr_arr_new[jj2].GetListArr();
         _Int *plist2_new = pmatr_arr_new[jj2].GetList2Arr();
         _Int *pia_new = pmatr_arr_new[jj2].GetIaArr();
         _Int *pja_new = pmatr_arr_new[jj2].GetJaArr();
         _Int *pja2_new = pmatr_arr_new[jj2].GetJa2Arr();
         _Flt *pa_new = pmatr_arr_new[jj2].GetAArr();
         k = (int)pnlist_new[jj2];
         plist_new[k] = plist2_sub[j*2];
         plist2_new[k] = plist2_sub[j*2+1];
         pia_new[k+1] = pia_new[k]+(pia_sub[j+1]-pia_sub[j]);
         pnlist_new[jj2]++;
         k = (int)pnz_list_new[jj2];
         for (kj=(int)pia_sub[j];kj<pia_sub[j+1];kj++) {
            kk = (int)pja_sub[kj];
            kk2 = (int)pja2_sub[kj];
            pja_new[k] = kk;
            pja2_new[k] = kk2;
            pa_new[k] = pa_sub[kj];
            k++;
         }
         pnz_list_new[jj2] = k;
      }
   }

// Add own data

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {
         int nlist_sub = pblkrows[iblk].GetNlist();
         _Int *plist_sub = pblkrows[iblk].GetListArr();
         _Int *plist2_sub = pblkrows[iblk].GetList2Arr();
         _Int *pia_sub = pblkrows[iblk].GetIaArr();
         _Int *pja_sub = pblkrows[iblk].GetJaArr();
         _Int *pja2_sub = pblkrows[iblk].GetJa2Arr();
         _Flt *pa_sub = pblkrows[iblk].GetAArr();
         _Int *plist_new = pmatr_arr_new[iblk].GetListArr();
         _Int *plist2_new = pmatr_arr_new[iblk].GetList2Arr();
         _Int *pia_new = pmatr_arr_new[iblk].GetIaArr();
         _Int *pja_new = pmatr_arr_new[iblk].GetJaArr();
         _Int *pja2_new = pmatr_arr_new[iblk].GetJa2Arr();
         _Flt *pa_new = pmatr_arr_new[iblk].GetAArr();
         for (j=0;j<nlist_sub;j++) {
            k = (int)pnlist_new[iblk];
            plist_new[k] = plist_sub[j];
            plist2_new[k] = plist2_sub[j];
            pia_new[k+1] = pia_new[k]+(pia_sub[j+1]-pia_sub[j]);
            pnlist_new[iblk]++;
            k = (int)pnz_list_new[iblk];
            for (kj=(int)pia_sub[j];kj<pia_sub[j+1];kj++) {
               kk = (int)pja_sub[kj];
               kk2 = (int)pja2_sub[kj];
               pja_new[k] = kk;
               pja2_new[k] = kk2;
               pa_new[k] = pa_sub[kj];
               k++;
            }
            pnz_list_new[iblk] = k;
         }
      }
   }

// Perform filtering and renumbering of the extended data

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == myid) {

         int nlist_new = pmatr_arr_new[iblk].GetNlist();
         int nzja_new = pmatr_arr_new[iblk].GetNzja();
         _Int *plist_new = pmatr_arr_new[iblk].GetListArr();
         _Int *plist2_new = pmatr_arr_new[iblk].GetList2Arr();
         _Int *pia_new = pmatr_arr_new[iblk].GetIaArr();
         _Int *pja_new = pmatr_arr_new[iblk].GetJaArr();
         _Int *pja2_new = pmatr_arr_new[iblk].GetJa2Arr();
         _Flt *pa_new = pmatr_arr_new[iblk].GetAArr();
         icycletot++;
         nlisttot = 0;
         for (j=0;j<nlist_new;j++) {
            jj2 = (int)plist2_new[j];
            if (pimasktot[jj2] != icycletot) {
               plisttot[nlisttot] = jj2;
               nlisttot++;
               pimasktot[jj2] = icycletot;
            }
         }

         sort (plisttot,plisttot+nlisttot);

         int jblk;

         for (j=0;j<nlisttot;j++) {
            jblk = (int)plisttot[j];
            pindtot2[jblk] = j;
         }

         vector<int> ia_list (nlisttot+1);
         vector<int> ja_list (nlist_new+1);
         vector<int> iptr (nlisttot+1);

         int *pia_list = &ia_list[0];
         int *pja_list = &ja_list[0];
         int *piptr = &iptr[0];

         for (j=0;j<=nlisttot;j++) pia_list[j] = 0;

         for (j=0;j<nlist_new;j++) {
            jj2 = (int)plist2_new[j];
            ind = pindtot2[jj2];
            pia_list[ind+1]++;
         }

         for (j=0;j<nlisttot;j++) pia_list[j+1] = pia_list[j]+pia_list[j+1];
         for (j=0;j<nlisttot;j++) piptr[j] = pia_list[j];

         for (j=0;j<nlist_new;j++) {
            jj = (int)plist_new[j];
            jj2 = (int)plist2_new[j];
            ind = (int)pindtot2[jj2];
            k = (int)piptr[ind];
            pja_list[k] = jj;
            piptr[ind]++;
         }

         for (j=0;j<nlisttot;j++) {
            sort (pja_list+pia_list[j],pja_list+pia_list[j+1]);
         }

         vector<int> imask_elems (nzja_new+1);
         vector<int> ia_elems (nlisttot+1);
         vector<int> ind_elems (nzja_new+1);

         int *pimask_elems = &imask_elems[0];
         int *pia_elems = &ia_elems[0];
         int *pind_elems = &ind_elems[0];

         for (j=0;j<nzja_new;j++) pimask_elems[j] = -1;
         for (j=0;j<=nlisttot;j++) pia_elems[j] = 0;

         for (j=0;j<nzja_new;j++) {
            jj2 = (int)pja2_new[j];
            if (pimasktot[jj2] == icycletot) {
               ind = pindtot2[jj2];
               pia_elems[ind+1]++;
            }
         }

         for (j=0;j<nlisttot;j++) pia_elems[j+1] = pia_elems[j]+pia_elems[j+1];

         for (j=0;j<nlisttot;j++) piptr[j] = pia_elems[j];

         for (j=0;j<nzja_new;j++) {
            jj2 = (int)pja2_new[j];
            if (pimasktot[jj2] == icycletot) {
               ind = (int)pindtot2[jj2];
               k = (int)piptr[ind];
               pind_elems[k] = j;
               piptr[ind]++;
            }
         }

         for (j=0;j<nlisttot;j++) {
            icycleblk++;
            for (k=pia_list[j];k<pia_list[j+1];k++) {
               kk = (int)pja_list[k];
               pimaskblk[kk] = icycleblk;
               piablk[kk] = k;
            }
            for (k=pia_elems[j];k<pia_elems[j+1];k++) {
               ind = (int)pind_elems[k];
               jj = (int)pja_new[ind];
               if (pimaskblk[jj] == icycleblk) {
                  pimask_elems[ind] = piablk[jj];
               }
            }
         }

         int nzja_flt = 0;

         for (j=0;j<nzja_new;j++) {
            if (pimask_elems[j] >= 0) nzja_flt++;
         }

         CMatrix<_Int,_Flt> a_flt;

         a_flt.ResizeAndSetAll (nlist_new, 0, nzja_flt, 0, nzja_flt);

         _Int *plist_flt = a_flt.GetListArr();
         _Int *pia_flt = a_flt.GetIaArr();
         _Int *pja_flt = a_flt.GetJaArr();
         _Flt *pa_flt = a_flt.GetAArr();

         for (j=0;j<nlist_new;j++) plist_flt[j] = (_Int)j;

         nzja_flt = 0;
         pia_flt[0] = 0;

         for (i=0;i<nlist_new;i++) {
            for (j=(int)pia_new[i];j<pia_new[i+1];j++) {
               if (pimask_elems[j] >= 0) {
                  pja_flt[nzja_flt] = pimask_elems[j];
                  pa_flt[nzja_flt] = pa_new[j];
                  nzja_flt++;
               }
            }
            pia_flt[i+1] = nzja_flt;
         }

         int njmax = 0;

         int njloc = 0;

         for (j=0;j<nlist_new;j++) {
            njloc = (int)(pia_flt[j+1]-pia_flt[j]);
            if (njloc > njmax) njmax = njloc;
         }

         vector<CSortInt> iiarr (njmax+1);
         vector<_Flt> elems (njmax+1);

         CSortInt *piiarr = &iiarr[0];
         _Flt *pelems = &elems[0];

         int ibeg, jloc;

         for (i=0;i<nlist_new;i++) {
            ibeg = (int)pia_flt[i];
            njloc = (int)(pia_flt[i+1]-pia_flt[i]);
            for (j=(int)pia_flt[i];j<pia_flt[i+1];j++) {
               jloc = j-ibeg;
               piiarr[jloc].ival = (int)pja_flt[j];
               piiarr[jloc].i2val = jloc;
               pelems[jloc] = pa_flt[j];
            }
            sort (piiarr,piiarr+njloc);
            for (j=(int)pia_flt[i];j<pia_flt[i+1];j++) {
               jloc = j-ibeg;
               pja_flt[j] = piiarr[jloc].ival;
               pa_flt[j] = pelems[piiarr[jloc].i2val];
            }
         }

         _matr_ext_arr[iblk].ReplaceFree (a_flt);

      }
   }

}

// Compute the packed size
//========================================================================================
template <typename _Int, typename _Flt>
long long CHMatrix<_Int,_Flt>::GetPackedSize () const 
{

   long long length = sizeof(int);

   length += (this->nzblk+2) * sizeof(long long);

   length += this->hmatr_str.GetPackedSize();

   if (this->nzblk > 0) {

      const CMatrix<_Int,_Flt> *pABlocks = &(this->asub_arr[0]);

      int i;

      for (i=0;i<this->nzblk;i++) {
         length += pABlocks[i].GetPackedSize();
      }

   }

   return length;

}

// Fill by packed data
//========================================================================================
template <typename _Int, typename _Flt>
void CHMatrix<_Int,_Flt>::FillPacked (long long _length, char *_obj) const 
{

   char* pLoc;

   pLoc = _obj;

   int *pHead;
   long long *pibsblk;

   pHead = (int *) pLoc;
   pLoc += sizeof(int);

   pibsblk = (long long *) pLoc;
   pLoc += (this->nzblk+2) * sizeof(long long);

   int sizeloc = this->hmatr_str.GetPackedSize();

   pibsblk[0] = pLoc-_obj;
   pLoc += sizeloc;

   pibsblk[1] = pLoc-_obj;

   int i;

   if (this->nzblk > 0) {

      const CMatrix<_Int,_Flt> *pABlocks = &(this->asub_arr[0]);

      for (i=0;i<this->nzblk;i++) {

         sizeloc = pABlocks[i].GetPackedSize ();

         pibsblk[i+1] = pLoc-_obj;
         pLoc += sizeloc;

      }

      pibsblk[this->nzblk+1] = pLoc-_obj;

   }

   pHead[0] = nzblk;

   if (pLoc-_obj != _length) {
      throw " CHMatrix<_Int,_Flt>::FillPacked: incorrect length on entry ";
   }

// Pack sparsity block

   sizeloc = (int)(pibsblk[1]-pibsblk[0]);

   this->hmatr_str.FillPacked (sizeloc, _obj+pibsblk[0]);

// Pack blocks

   if (this->nzblk > 0) {

      const CMatrix<_Int,_Flt> *pABlocks = &(this->asub_arr[0]);

      for (i=0;i<this->nzblk;i++) {
         sizeloc = (int)(pibsblk[i+2]-pibsblk[i+1]);
         pABlocks[i].FillPacked (sizeloc, _obj+pibsblk[i+1]);
      }

   }

}

// Fill by packed data
//========================================================================================
template <typename _Int, typename _Flt>
void CHMatrix<_Int,_Flt>::UnPack (long long _length, char *_obj) 
{

// Get head data

   char *pLoc;

   pLoc = _obj;

   int *pHead;

   pHead = (int *) pLoc;
   pLoc += sizeof(int);

   int nzblk_loc = pHead[0];

   long long *pibsblk;

   pibsblk = (long long *) pLoc;
   pLoc += (nzblk_loc+2)*sizeof(long long);

   this->nzblk = nzblk_loc;

// Unpack sparsity

   int sizeloc = (int)(pibsblk[1]-pibsblk[0]);

   this->hmatr_str.UnPack (sizeloc, pLoc);

   pLoc += sizeloc;

// Pack blocks

   this->asub_arr.resize (nzblk_loc+1);

   CMatrix<_Int,_Flt> *pABlocks = &(this->asub_arr[0]);

   int i;

   for (i=0;i<nzblk;i++) {

      sizeloc = (int)(pibsblk[i+2]-pibsblk[i+1]);
      pABlocks[i].UnPack (sizeloc, pLoc);
      pLoc += sizeloc;

   }

};

// Filter the extended lists for backward only extention
//========================================================================================
template <typename _Int, typename _Flt>
void CHMatrix<_Int,_Flt>::FilterListsBack (int _myid, int _nblks, int *_blk2cpu, 
                                                  int *_nlist_ext_arr, vector<int> *_list_ext_arr) 
{

   int iblk, i, jj, jj2;

   for (iblk=0;iblk<_nblks;iblk++) {
      if (_blk2cpu[iblk] == _myid) {
         int *plist_ext = &(_list_ext_arr[iblk][0]);
         vector<int> list_flt (2*_nlist_ext_arr[iblk]+1);
         int *plist_flt = &list_flt[0];
         int nlistflt = 0;
         for (i=0;i<_nlist_ext_arr[iblk];i++) {
            jj = plist_ext[2*i];
            jj2 = plist_ext[2*i+1];
            if (jj2 < iblk) {
               plist_flt[nlistflt*2] = jj;
               plist_flt[nlistflt*2+1] = jj2;
               nlistflt++;
            }
         }
         _nlist_ext_arr[iblk] = nlistflt;
         _list_ext_arr[iblk].swap(list_flt);
      }
   }

};

// Compute Ja2 array by binary search
//========================================================================================
template <typename _Int, typename _Flt>
void CHMatrix<_Int,_Flt>::ComputeJa2 (int _nblks, long long *_blks,
                                             int _nlist, const _Int *_ia, const _Int *_ja, 
                                             _Int *_ja2) 
{
   int nzja = (int)_ia[_nlist];

   int i;
   int iblk, iblkprev;
   long long jj;

   iblkprev = -1;

   for (i=0;i<nzja;i++) {
      jj = _ja[i];
      if (iblkprev >= 0 && iblkprev < _nblks) {
         if (jj >= _blks[iblkprev] && jj < _blks[iblkprev+1]) {
            iblk = iblkprev;
         } else {
            iblk = BinarySearch (jj, _nblks, _blks, iblkprev);
         } 
      } else {
         iblk = BinarySearch (jj, _nblks, _blks, iblkprev);
      }

      _ja2[i] = iblk;
      iblkprev = iblk;

   }
};

/// @brief Split structural matrix into the set of matrices
//========================================================================================
template <typename _Int, typename _Flt>
void CHMatrix<_Int,_Flt>::SplitMatrSpIntoHMatrSp (int _nblks, long long *_blks, const CMatrix<_Int,_Flt> &_amatr, CHMatrix<_Int,_Flt> &_ahmatr)
{

// Open sparsity

   const int nlistloc_ablk = _amatr.GetNlist();
   const int nzjaloc_ablk = _amatr.GetNzja();
   const _Int *pialoc_ablk = _amatr.GetIaArr();
   const _Int *pjaloc_ablk = _amatr.GetJaArr();

   if (nlistloc_ablk != _blks[_nblks]) throw " CHMatrix<_Int,_Flt>::SplitMatrSpIntoHMatrSp: incorrect matrix on entry ";

// Compute block indices

   vector<_Int> ja2loc(nzjaloc_ablk+1);

   _Int *pja2loc = &ja2loc[0];

   CHMatrix<_Int,_Flt>::ComputeJa2 (_nblks, _blks,
                                          nlistloc_ablk, pialoc_ablk, pjaloc_ablk,
                                          pja2loc);

// Compute block sparsity

   vector<int>  imaskblk(_nblks);
   vector<int>  listblk(_nblks);

   int *pimaskblk = &imaskblk[0];
   int *plistblk = &listblk[0];

   int i;

   for (i=0;i<_nblks;i++) pimaskblk[i] = -1;

   int icycle = -1;

   vector<int> ia_blk(_nblks+1);

   int *pia_blk = &ia_blk[0];

   for (i=0;i<=_nblks;i++) pia_blk[i] = 0;

   int ihblk, nlistloc, j, jhblk;

   for (ihblk=0;ihblk<_nblks;ihblk++) {
      icycle++;
      nlistloc = 0;
      for (i=(int)_blks[ihblk];i<_blks[ihblk+1];i++) {
         for (j=(int)pialoc_ablk[i];j<pialoc_ablk[i+1];j++) {
            jhblk = (int)pja2loc[j];
            if (pimaskblk[jhblk] != icycle) {
               nlistloc++;
               pimaskblk[jhblk] = icycle;
            }
         }
      }
      pia_blk[ihblk+1] = nlistloc;
   }

   for (i=0;i<_nblks;i++) pia_blk[i+1] = pia_blk[i]+pia_blk[i+1];

   int nzjablk = pia_blk[_nblks];

   vector<int> ja_blk(nzjablk+1);

   int *pja_blk = &ja_blk[0];

   int ibeg;

   for (ihblk=0;ihblk<_nblks;ihblk++) {
      icycle++;
      nlistloc = 0;
      for (i=(int)_blks[ihblk];i<_blks[ihblk+1];i++) {
         for (j=(int)pialoc_ablk[i];j<pialoc_ablk[i+1];j++) {
            jhblk = (int)pja2loc[j];
            if (pimaskblk[jhblk] != icycle) {
               plistblk[nlistloc] = jhblk;
               nlistloc++;
               pimaskblk[jhblk] = icycle;
            }
         }
      }
      sort(plistblk,plistblk+nlistloc);
      ibeg = pia_blk[ihblk];
      for (i=0;i<nlistloc;i++) pja_blk[ibeg+i] = plistblk[i];
   }

// Count numbers of elements in each block

   vector<int> nz_blk(nzjablk+1);

   int *pnz_blk = &nz_blk[0];

   for (i=0;i<nzjablk;i++) pnz_blk[i] = 0;

   int ind;

   for (ihblk=0;ihblk<_nblks;ihblk++) {
      for (i=pia_blk[ihblk];i<pia_blk[ihblk+1];i++) {
         jhblk = pja_blk[i];
         plistblk[jhblk] = i;
      }
      for (i=(int)_blks[ihblk];i<_blks[ihblk+1];i++) {
         for (j=(int)pialoc_ablk[i];j<pialoc_ablk[i+1];j++) {
            jhblk = (int)pja2loc[j];
            ind = plistblk[jhblk];
            pnz_blk[ind]++;
         }
      }
   }

// Rewrite data block by block

   vector<int> iptr(nzjablk+1);
   vector<vector<_Int> > iindarr(nzjablk+1);
   vector<vector<_Int> > jindarr(nzjablk+1);
   vector<vector<_Int> > j2indarr(nzjablk+1);

   int *piptr = &iptr[0];
   vector<_Int> *piindarr = &iindarr[0];
   vector<_Int> *pjindarr = &jindarr[0];

   for (i=0;i<nzjablk;i++) piptr[i] = 0;

   int nzjmax = 0;

   for (i=0;i<nzjablk;i++) {
      piindarr[i].resize(pnz_blk[i]+1);
      pjindarr[i].resize(pnz_blk[i]+1);
      if (pnz_blk[i] > nzjmax) nzjmax = pnz_blk[i];
   };

   vector<CSortInt> iiarr(nzjmax+1);

   CSortInt *piiarr = &iiarr[0];

   _ahmatr.ResizeASub (nzjablk);
   _ahmatr.SetNzblk(nzjablk);

   CMatrix<_Int,_Flt> *pABlocks = _ahmatr.GetASubArr();

   int jj, k, iloc, jloc, niloc;

   for (ihblk=0;ihblk<_nblks;ihblk++) {
      for (i=pia_blk[ihblk];i<pia_blk[ihblk+1];i++) {
         jhblk = pja_blk[i];
         plistblk[jhblk] = i;
      }
      for (i=(int)_blks[ihblk];i<_blks[ihblk+1];i++) {
         for (j=(int)pialoc_ablk[i];j<pialoc_ablk[i+1];j++) {
            jj = (int)pjaloc_ablk[j];
            jhblk = (int)pja2loc[j];
            ind = plistblk[jhblk];
            iloc = i-(int)_blks[ihblk];
            jloc = jj-(int)_blks[jhblk];
            k = piptr[ind];
            _Int *ppiindarr = &(piindarr[ind][0]);
            _Int *ppjindarr = &(pjindarr[ind][0]);
            ppiindarr[k] = (_Int)iloc;
            ppjindarr[k] = (_Int)jloc;
            piptr[ind]++;
         }
      }
      for (i=pia_blk[ihblk];i<pia_blk[ihblk+1];i++) {
         jhblk = pja_blk[i];
         _Int *ppiindarr =  &(piindarr[i][0]);
         _Int *ppjindarr = &(pjindarr[i][0]);
         for (k=0;k<pnz_blk[i];k++) {
            piiarr[k].ival = (int)ppiindarr[k];
            piiarr[k].i2val = k;
         }
         sort (piiarr,piiarr+pnz_blk[i]);
         int nlistloc = 0;
         int irow_prev = -1;
         for (k=0;k<pnz_blk[i];k++) {
            if (piiarr[k].ival != irow_prev) {
               irow_prev = piiarr[k].ival;
               nlistloc++;
            }
         }

         CMatrix<_Int,_Flt> ablk;
         ablk.ResizeAndSetAllSp (nlistloc,0,pnz_blk[i],0);

         _Int *plist_loc = ablk.GetListArr();
         _Int *pia_loc = ablk.GetIaArr();
         _Int *pja_loc = ablk.GetJaArr();

         irow_prev = -1;

         pia_loc[0] = 0;

         nlistloc = 0;

         for (k=0;k<pnz_blk[i];k++) {
            if (piiarr[k].ival != irow_prev) {
               irow_prev = piiarr[k].ival;
               plist_loc[nlistloc] = (_Int)irow_prev;
               nlistloc++;
            }
            pia_loc[nlistloc] = (_Int)(k+1);
         }
         for (k=0;k<pnz_blk[i];k++) {
            pja_loc[k] = ppjindarr[piiarr[k].i2val];
         }
         for (k=0;k<nlistloc;k++) {
            ibeg = (int)pia_loc[k];
            niloc = (int)(pia_loc[k+1]-pia_loc[k]);
            sort(pja_loc+ibeg,pja_loc+ibeg+niloc);
         }

         pABlocks[i].ReplaceFree (ablk);

      }
   }

   CMatrix<int,float> *pAHBlkStr = _ahmatr.GetHMatrStr();

   pAHBlkStr->ResizeAndSetAllSp (_nblks,0,nzjablk,0);

   vector<int> *pia_HBlk = pAHBlkStr->GetIa();
   vector<int> *pja_HBlk = pAHBlkStr->GetJa();

   pia_HBlk->swap (ia_blk);
   pja_HBlk->swap (ja_blk);

   int *plist_blk = pAHBlkStr->GetListArr();

   for (i=0;i<_nblks;i++) plist_blk[i] = i;

};

// Print hmatrix data
//========================================================================================
template <typename _Int, typename _Flt>
void CHMatrix<_Int,_Flt>::PrintHMatrix (ofstream &_fout) 
{

   _fout << " CHMatrix:" << endl;

   _fout << "    HMatrStr:" << endl;

   CIlu2<int,float>::PrintMatrix (_fout, this->hmatr_str);

   _fout << "    Nzblk = " << this->nzblk << endl;

   _fout << "    ABlocks: " << endl;

   if (this->nzblk > 0) {

      int i;

      CMatrix<_Int,_Flt> *pABlocks = this->GetASubArr ();

      for (i=0;i<this->nzblk;i++) {
         _fout << "       Iblk = " << i << endl;
         CIlu2<_Int,_Flt>::PrintMatrix (_fout, pABlocks[i]);
      }

   }

}

// Print sparsity with boxes that show nonzero blocks
//========================================================================================
template <typename _Int, typename _Flt>
void CHMatrix<_Int,_Flt>::Str2PsBox (int _collap, const CMatrix<_Int,_Flt> &_amatr, char *_fname, int _nblks, int *_blks)
{

// Compute new blocks partitioning

   vector<int> blksnw (_nblks+2);

   int *pblksnw = &blksnw[0];

   const int nsupmx = _amatr.GetNlist();

   if (_nblks == 0) {
      pblksnw[0] = 0;
      pblksnw[1] = (nsupmx+_collap-1) / _collap;
   } else {
      pblksnw[0] = 0;
      for (int iblk=0;iblk<_nblks;iblk++) {
         int ni = _blks[iblk+1]-_blks[iblk];
         int niloc = (ni+_collap-1) / _collap;
         pblksnw[iblk+1] = pblksnw[iblk]+niloc;
      }
   }

// Compute nd2sp and sp2nd arrays

   int nnew;

   if (_nblks == 0) {
      nnew = pblksnw[_nblks+1];
   } else {
      nnew = pblksnw[_nblks];
   }

   vector<int> nd2sp (_amatr.GetNlist()+1);
   vector<int> sp2nd (nnew+1);

   int *pnd2sp = &nd2sp[0];
   int *psp2nd = &sp2nd[0];

   int i, j;

   if (_nblks == 0) {
      psp2nd[0] = 0;
      for (i=0;i<nnew-1;i++) {
         psp2nd[i+1] = psp2nd[i] + _collap;
      }
      psp2nd[nnew] = nsupmx;
   } else {
      psp2nd[0] = 0;
      for (int iblk=0;iblk<_nblks;iblk++) {
         for (i=pblksnw[iblk];i<pblksnw[iblk+1]-1;i++) {
            psp2nd[i+1] = psp2nd[i] + _collap;
         }
         psp2nd[pblksnw[iblk+1]] = _blks[iblk+1];
      }
   }

   for (i=0;i<nnew;i++) {
      for (j=psp2nd[i];j<psp2nd[i+1];j++) {
         pnd2sp[j] = i;
      }
   }

// Count the total number of elements in the collapsed matrix

   int isup, jsup, jj;

   int icycle=0;

   vector<int> imask (nnew+1);
   vector<int> lstloc (nnew+1);

   int *pimask = &imask[0];
   int *plstloc = &lstloc[0];

   for (i=0;i<nnew;i++) pimask[i] = icycle;

   int nz = 0;
   int nlstloc;

   const _Int *pia = _amatr.GetIaArr();
   const _Int *pja = _amatr.GetJaArr();

   for (isup=0;isup<nnew;isup++) {
      icycle++;
      nlstloc = 0;
      for (i=psp2nd[isup];i<psp2nd[isup+1];i++) {
         for (j=(int)pia[i];j<pia[i+1];j++) {
            jj = (int)pja[j];
            jsup = pnd2sp[jj];
            if (pimask[jsup] != icycle) {
               plstloc[nlstloc] = jsup;
               pimask[jsup] = icycle;
               nlstloc++;
            }
         }
      }
      nz += nlstloc;
   }

// Compute collapsed matrix

   CMatrix<int,float> temp;

   temp.ResizeAndSetAllSp(nnew,0,nz,0);

   int *plistt = temp.GetListArr();
   int *piat = temp.GetIaArr();
   int *pjat = temp.GetJaArr();

   piat[0] = 0;
   for (isup=0;isup<nnew;isup++) {
      icycle++;
      nlstloc = 0;
      for (i=psp2nd[isup];i<psp2nd[isup+1];i++) {
         for (j=(int)pia[i];j<pia[i+1];j++) {
            jj = (int)pja[j];
            jsup = pnd2sp[jj];
            if (pimask[jsup] != icycle) {
               plstloc[nlstloc] = jsup;
               pimask[jsup] = icycle;
               nlstloc++;
            }
         }
      }
      if (nlstloc != 0) sort (plstloc, plstloc + nlstloc);
      for (j=0;j<nlstloc;j++) pjat[piat[isup]+j] = (int)plstloc[j];
      piat[isup+1] = piat[isup]+(int)nlstloc;
   }

   for (i=0;i<nnew;i++) plistt[i] = (int)i;

// Print collapsed sparsity structure

   CHMatrix<int,float>::Str2PsBox (temp, _fname, _nblks, pblksnw);

};


// Print sparsity with boxes that show nonzero blocks
//========================================================================================
template <typename _Int, typename _Flt>
void CHMatrix<_Int,_Flt>::Str2PsBox (const CMatrix<_Int,_Flt> &_amatr, char *_fname, int _nblks, int *_blks)
{

// Open output file

   ofstream fout (_fname);

   if (!fout.is_open()) {
      cout << " Error: File named " << _fname << " is not opened !!!" << endl;
      return;
   }

   static float default_bg[3] = {(float)0.9, (float)0.9, (float)1.0};
   static float default_fg[3] = {(float)0.7, (float)0.7, (float)1.0};

// Write the header information

   fout << "%!PS-Adobe-2.0" << endl;
   fout << "%%BoundingBox: 0 0 600 600" << endl;
   fout << "/m {moveto} def % x y" << endl;
   fout << "/l {lineto} def % x y" << endl;
   fout << "/s {stroke} def % x y" << endl;
   fout << "/n {newpath} def % x y" << endl;
   fout << "/c {closepath} def % x y" << endl;

// Parameters of the local window

   const int nsupmx = _amatr.GetNlist();

   double s, s1;

   s1 = 50.0e0;
   s = 500.0e0 / ((double) nsupmx);

// Print the bounding window

   double dx, dy;

   fout << SetPw << 0.03 << " setlinewidth" << endl;

   fout << " n" << endl;
   dx = s1+s*0.5; dy = s1+s*0.5;
   Round (dx,dy);
   fout << SetPw << dx << SetPw << dy << " m" << endl;
   dx = s1+s*(nsupmx+0.5); dy = s1+s*0.5;
   Round (dx,dy);
   fout << SetPw << dx << SetPw << dy << " l" << endl;
   dx = s1+s*(nsupmx+0.5); dy = s1+s*(nsupmx+0.5);
   Round (dx,dy);
   fout << SetPw << dx << SetPw << dy << " l" << endl;
   dx = s1+s*0.5; dy = s1+s*(nsupmx+0.5);
   Round (dx,dy);
   fout << SetPw << dx << SetPw << dy << " l" << endl;
   dx = s1+s*0.5; dy = s1+s*0.5;
   Round (dx,dy);
   fout << SetPw << dx << SetPw << dy << " l s c" << endl;

// Print blocks partitioning

   int i, j, k;
   double x, y;
   double x1, y1;

   if (_nblks > 0) {

      vector<long long> blks_new (_nblks+2);

      long long *pblks_new = &blks_new[0];

      int nblks_new = 0;

      pblks_new[0] = 0;

      for (i=0;i<_nblks;i++) {
         if (_blks[i+1] > _blks[i]) {
            pblks_new[nblks_new+1] = _blks[i+1];
            nblks_new++;
         }
      }

      if (pblks_new[nblks_new] < _amatr.GetNlist()) {
         pblks_new[nblks_new+1] = _amatr.GetNlist();
         nblks_new++;
      }

      CHMatrix<_Int,_Flt> ahblk;

      CHMatrix<_Int,_Flt>::SplitMatrSpIntoHMatrSp (nblks_new, pblks_new, _amatr, ahblk);

      CMatrix<int,float> *pahblkstr = ahblk.GetHMatrStr();

      int nlistloc = pahblkstr->GetNlist();
      int *plistloc = pahblkstr->GetListArr();
      int *pialoc = pahblkstr->GetIaArr();
      int *pjaloc = pahblkstr->GetJaArr();

      for(int ilist = 0;ilist < nlistloc;ilist++) {
         int i = plistloc[ilist];
         for(int jlist=pialoc[ilist];jlist<pialoc[ilist+1];jlist++) {
            j = pjaloc[jlist];
            x = pblks_new[j]+0.5;
            x1 = pblks_new[j + 1]+0.5;
            y = _amatr.GetNlist()-((int)pblks_new[i]+0.5)+1;
            y1 = _amatr.GetNlist()-((int)pblks_new[i + 1]+0.5)+1;
            fout << default_bg[0] << " " << default_bg[1] << " " << default_bg[2] << " setrgbcolor" << endl;
            fout << "[" << s1+s*x << " " << s1+s*y << " " << s*(x1 - x) << " " << s*(y1 - y) << "] rectfill" << endl;
            fout << default_fg[0] << " " << default_fg[1] << " " << default_fg[2] << " setrgbcolor" << endl;
            fout << "[" << s1+s*x << " " << s1+s*y << " " << s*(x1 - x) << " " << s*(y1 - y) << "] rectstroke" << endl;
         }
      }
      fout << "0 0 0 setrgbcolor" << endl;
   }

// Compute radius of the circle

   double r;

   r = (log10(1.0e0*nsupmx)+1.0) / 4.0 ;
   if (r > 1.0e0) r=1.0e0;

   r = s * r / 2.0;

   if (r < 0.01e0) r=0.01e0;

   fout << SetPw << 2*r << " setlinewidth" << endl;
   fout << " /d {moveto currentpoint " << SetPw << r;
   fout << " 0 360 arc fill} def % x y" << endl;

// Print the sparsity

   int i1, j1;

   const _Int *plist = _amatr.GetListArr();
   const _Int *pia = _amatr.GetIaArr();
   const _Int *pja = _amatr.GetJaArr();

   fout << " n" << endl;
   for (int ilist=0;ilist<_amatr.GetNlist();ilist++) {
      i = (int)plist[ilist];
      for (k=(int)pia[ilist];k<pia[ilist+1];k++) {
         j = (int)pja[k];
         i1 = i+1;
         j1 = j+1;
         x = (double) j1;
         y = (double) (_amatr.GetNlist()-i1+1);
         dx = s1+s*x-r; dy = s1+s*y;
         Round (dx,dy);
         fout << SetPw << dx << SetPw << dy << " m" << endl;
         dx = s1+s*x+r; dy = s1+s*y;
         Round (dx,dy);
         fout << SetPw << dx << SetPw << dy << " l" << endl;
      }
   }

// Write the footer

   fout << " s c" << endl;
   fout << " showpage" << endl;

// Close output file

   fout.close();

};

// Init MvmA data
//========================================================================================
template <typename _Int, typename _Flt, typename _FltVect>
void CMvmPar<_Int,_Flt,_FltVect>::InitMvmA (CHMatrix<_Int,_Flt> *_hmatr_arr) 
{

// Get control data

   void *pcomm_loc = this->pcomm;

   int nblks_loc            = this->nblks;
   long long *pblks_loc     = this->pblks;
   int *pblk2cpu_loc        = this->pblk2cpu;
   int ni_cpu_loc           = this->ni_cpu;
   int *pibsblk_loc         = &this->ibsblk[0];
   int nlistblk_own_loc     = this->nlistblk_own;
   int *plistblk_own_loc    = &this->listblk_own[0];

   int nproc = CExchange::GetNproc (pcomm_loc);
   int myid = CExchange::GetMyid (pcomm_loc);

// Compute maximal block size

   int nimax = 0;

   int i, niloc;

   for (i=0;i<nblks_loc;i++) {
      niloc = (int)(pblks_loc[i+1]-pblks_loc[i]);
      if (niloc > nimax) nimax = niloc;
   }

// Store pointer to matrix data

   this->phmatr = _hmatr_arr;

// Compute list of external blocks

   vector<int> imaskblk (nblks_loc+1);
   vector<int> listblk (nblks_loc+1);
   vector<int> indblk (nblks_loc+1);

   int *pimaskblk = &imaskblk[0];
   int *plistblk = &listblk[0];
   int *pindblk = &indblk[0];

   for (i=0;i<nblks_loc;i++) pimaskblk[i] = -1;

   int icycleblk = -1;

   int j, iblk, jblk;

   icycleblk++;

   int nlistblk_ext = 0;

   for (i=0;i<nlistblk_own_loc;i++) {
      iblk = plistblk_own_loc[i];
      CMatrix<int,float> *phmatr_str_loc = _hmatr_arr[iblk].GetHMatrStr();
      int nzja_hblk = phmatr_str_loc->GetNzja();
      int *pja_hblk = phmatr_str_loc->GetJaArr();
      for (j=0;j<nzja_hblk;j++) {
         jblk = pja_hblk[j];
         if (pblk2cpu_loc[jblk] != myid && pimaskblk[jblk] != icycleblk) {
            plistblk[nlistblk_ext] = jblk;
            nlistblk_ext++;
            pimaskblk[jblk] = icycleblk;
         }
      }
   }

   sort (plistblk,plistblk+nlistblk_ext);

   for (i=0;i<nlistblk_ext;i++) {
      iblk = plistblk[i];
      pindblk[iblk] = i;
   }

   this->nblks_recvs = nlistblk_ext;
   this->listblk_recvs.resize (nlistblk_ext+1);

   int *plistblk_recvs = &this->listblk_recvs[0];

   for (i=0;i<nlistblk_ext;i++) plistblk_recvs[i] = plistblk[i];

// Compute the set of indices to be received and multiplications lists

   vector<int> ia_list_rcv (nlistblk_ext+1);
   vector<vector<int> > ja_pair_rcv (nlistblk_ext+1);

   int *pia_list_rcv = &ia_list_rcv[0];
   vector<int> *pja_pair_rcv = &ja_pair_rcv[0];

   for (i=0;i<=nlistblk_ext;i++) pia_list_rcv[i] = 0;

   int ind;

   for (i=0;i<nlistblk_own_loc;i++) {
      iblk = plistblk_own_loc[i];
      CMatrix<int,float> *phmatr_str_loc = _hmatr_arr[iblk].GetHMatrStr();
      int nzja_hblk = phmatr_str_loc->GetNzja();
      int *pja_hblk = phmatr_str_loc->GetJaArr();
      for (j=0;j<nzja_hblk;j++) {
         jblk = pja_hblk[j];
         if (pblk2cpu_loc[jblk] != myid) {
            ind = pindblk[jblk];
            pja_pair_rcv[ind].push_back (iblk);
            pja_pair_rcv[ind].push_back (j);
            pia_list_rcv[ind+1]++;
         }
      }
   }

   for (i=0;i<nlistblk_ext;i++) pia_list_rcv[i+1] = pia_list_rcv[i]+pia_list_rcv[i+1];

   this->ialist_recvs.resize (nlistblk_ext+1);

   int *pialist_recvs = &this->ialist_recvs[0];

   for (i=0;i<=nlistblk_ext;i++) pialist_recvs[i] = pia_list_rcv[i];

   int nzja_pairs = pialist_recvs[nlistblk_ext];

   this->japairs_recvs.resize (2*nzja_pairs+1);

   int *pjapairs_recvs = &this->japairs_recvs[0];

   int ibeg, jloc;

   for (i=0;i<nlistblk_ext;i++) {
      int *ppja_pair_rcv = &(pja_pair_rcv[i][0]);
      ibeg = pia_list_rcv[i];
      for (j=pia_list_rcv[i];j<pia_list_rcv[i+1];j++) {
         jloc = j-ibeg;
         pjapairs_recvs[j*2] = ppja_pair_rcv[jloc*2];
         pjapairs_recvs[j*2+1] = ppja_pair_rcv[jloc*2+1];
      }
   }

   vector<int> imask_work (nimax+1);
   vector<int> list_work (nimax+1);

   int *pimask_work = &imask_work[0];
   int *plist_work = &list_work[0];

   for (i=0;i<nimax;i++) pimask_work[i] = -1;

   int icycle_work = -1;

   vector<int> ia_blk_rcv (nlistblk_ext+1);
   vector<vector<int> > ja_blk_rcv (nlistblk_ext+1);

   int *pia_blk_rcv = &ia_blk_rcv[0];
   vector<int> *pja_blk_rcv = &ja_blk_rcv[0];

   for (i=0;i<=nlistblk_ext;i++) pia_blk_rcv[i] = 0;

   int nlistloc, k, jj;

   for (i=0;i<nlistblk_ext;i++) {
      icycle_work++;
      nlistloc = 0;
      for (j=pialist_recvs[i];j<pialist_recvs[i+1];j++) {
         iblk = pjapairs_recvs[j*2];
         ind = pjapairs_recvs[j*2+1];
         CMatrix<_Int,_Flt> *pasub_loc = _hmatr_arr[iblk].GetASubArr();
         int nzjaloc = pasub_loc[ind].GetNzja ();
         _Int *pjaloc = pasub_loc[ind].GetJaArr ();
         for (k=0;k<nzjaloc;k++) {
            jj = (int)pjaloc[k];
            if (pimask_work[jj] != icycle_work) {
               plist_work[nlistloc] = jj;
               nlistloc++;
               pimask_work[jj] = icycle_work;
            }
         }
      }
      sort (plist_work,plist_work+nlistloc);
      pia_blk_rcv[i+1] = nlistloc;
      pja_blk_rcv[i].resize (nlistloc+1);
      int *ppja_blk_rcv = &(pja_blk_rcv[i][0]);
      for (j=0;j<nlistloc;j++) ppja_blk_rcv[j] = plist_work[j];
   }

   for (i=0;i<nlistblk_ext;i++) pia_blk_rcv[i+1] = pia_blk_rcv[i]+pia_blk_rcv[i+1];

   this->iablk_recvs.resize (nlistblk_ext+1);

   int *piablk_recvs = &this->iablk_recvs[0];

   for (i=0;i<=nlistblk_ext;i++) piablk_recvs[i] = pia_blk_rcv[i];

   int nzja_ind = pia_blk_rcv[nlistblk_ext];

   this->ind_recvs.resize (nzja_ind+1);

   int *pind_recvs = &this->ind_recvs[0];

   for (i=0;i<nlistblk_ext;i++) {
      int *ppja_pair_rcv = &(pja_blk_rcv[i][0]);
      ibeg = pia_blk_rcv[i];
      for (j=pia_blk_rcv[i];j<pia_blk_rcv[i+1];j++) {
         jloc = j-ibeg;
         pind_recvs[j] = ppja_pair_rcv[jloc];
      }
   }

   this->x_recv.resize (nzja_ind+1);

// Finally create cpu recv data

   vector<int> imaskcpu (nproc);
   vector<int> listcpu (nproc);
   vector<int> indcpu (nproc);

   int *pimaskcpu = &imaskcpu[0];
   int *plistcpu = &listcpu[0];
   int *pindcpu = &indcpu[0];

   for (i=0;i<nproc;i++) pimaskcpu[i] = -1;

   int icyclecpu = -1;

   icyclecpu++;

   int nlistcpu = 0;

   int jcpu;

   for (i=0;i<nlistblk_ext;i++) {
      jblk = plistblk_recvs[i];
      jcpu = pblk2cpu_loc[jblk];
      if (pimaskcpu[jcpu] != icyclecpu) {
         plistcpu[nlistcpu] = jcpu;
         nlistcpu++;
         pimaskcpu[jcpu] = icyclecpu;
      }
   }

   sort (plistcpu,plistcpu+nlistcpu);

   for (i=0;i<nlistcpu;i++) {
      jcpu = plistcpu[i];
      pindcpu[jcpu] = i;
   }

   this->nrecvs = nlistcpu;

   this->rcv2cpu.resize (nlistcpu+1);
   this->ia_recvs.resize (nlistcpu+1);
   this->ja_recvs.resize (2*nzja_ind+1);

   int *prcv2cpu = &this->rcv2cpu[0];
   int *pia_recvs = &this->ia_recvs[0];
   int *pja_recvs = &this->ja_recvs[0];

   for (i=0;i<nlistcpu;i++) prcv2cpu[i] = plistcpu[i];
   for (i=0;i<=nlistcpu;i++) pia_recvs[i] = 0;

   for (i=0;i<nlistblk_ext;i++) {
      jblk = plistblk_recvs[i];
      jcpu = pblk2cpu_loc[jblk];
      ind = pindcpu[jcpu];
      pia_recvs[ind+1] += (pia_blk_rcv[i+1]-pia_blk_rcv[i]);
   }

   for (i=0;i<nlistcpu;i++) pia_recvs[i+1] = pia_recvs[i]+pia_recvs[i+1];

   for (i=0;i<nlistcpu;i++) plistcpu[i] = pia_recvs[i];

   for (i=0;i<nlistblk_ext;i++) {
      jblk = plistblk_recvs[i];
      jcpu = pblk2cpu_loc[jblk];
      ind = pindcpu[jcpu];
      k = plistcpu[ind];
      for (j=pia_blk_rcv[i];j<pia_blk_rcv[i+1];j++) {
         jj = pind_recvs[j];
         pja_recvs[k*2] = jj;
         pja_recvs[k*2+1] = jblk;
         pind_recvs[j] = k;
         k++;
      }
      plistcpu[ind] = k;
   }

// Prepare exchange data

   vector<CHMatrix<int,_Flt> > hblk_send (nlistcpu+1);

   CHMatrix<int,_Flt> *phblk_send = &hblk_send[0];

   for (i=0;i<nlistcpu;i++) {
      phblk_send[i].SetNzblk (1);
      phblk_send[i].ResizeASub (1);
      CMatrix<int,_Flt> *pA_sub = phblk_send[i].GetASubArr();
      CMatrix<int,_Flt> ablk_temp;
      int niloc = pia_recvs[i+1]-pia_recvs[i];
      ablk_temp.ResizeAndSetAllSp (0, 0, niloc*2, 0);
      int *pjaloc = ablk_temp.GetJaArr ();
      ibeg = pia_recvs[i];
      for (j=pia_recvs[i];j<pia_recvs[i+1];j++) {
         jloc = j-ibeg;
         pjaloc[jloc*2] = pja_recvs[j*2];
         pjaloc[jloc*2+1] = pja_recvs[j*2+1];
      }
      pA_sub->ReplaceFree (ablk_temp);
   }

// Pack send data

   vector<int> CpuIDSend (nlistcpu);
   vector<vector<char> > ObjSend (nlistcpu);

   int *pCpuIDSend = NULL;
   vector<char> *pObjSend = NULL;

   if (nlistcpu > 0) {

      pCpuIDSend = &CpuIDSend[0];
      pObjSend = &ObjSend[0];

   }

   long long isize;
   char *pobj;

   for (i=0;i<nlistcpu;i++) {
      pCpuIDSend[i] = prcv2cpu[i];
      isize = phblk_send[i].GetPackedSize();
      pObjSend[i].resize ((size_t)isize);
      pobj = &(pObjSend[i][0]);
      phblk_send[i].FillPacked (isize, pobj);
      phblk_send[i].Clean ();
   }

// Exchange

   vector<int> CpuIDRecv;
   vector<vector<char> > ObjRecv;

   CExchange::DataExchange (pcomm_loc, CpuIDSend, ObjSend, CpuIDRecv, ObjRecv);

// Free send data

   {
      vector<int> CpuIDSend_temp;
      vector<vector<char> > ObjSend_temp;
      CpuIDSend.swap (CpuIDSend_temp);
      ObjSend.swap (ObjSend_temp);
   }

// Unpack receive data

   int nrecv_loc = (int) CpuIDRecv.size();

   vector<char> *pObjRecv = NULL;
   int *pCpuIDRecv = NULL;

   if (nrecv_loc > 0) {
      pObjRecv = &ObjRecv[0];
      pCpuIDRecv = &CpuIDRecv[0];
   }

   vector<CHMatrix<int,_Flt> > hblk_recv (nrecv_loc+1);

   CHMatrix<int,_Flt> *phblk_recv = &hblk_recv[0];

   for (i=0;i<nrecv_loc;i++) {
      isize = (long long) pObjRecv[i].size();
      pobj = &(pObjRecv[i][0]);
      phblk_recv[i].UnPack (isize, pobj);
   }

// Free recv data

   {
      vector<vector<char> > ObjRecv_temp;
      ObjRecv.swap (ObjRecv_temp);
   }

// Compute correct ordering of cpu data

   vector<CSortInt> iiarr (nrecv_loc+1);
   CSortInt *piiarr = &iiarr[0];

   for (i=0;i<nrecv_loc;i++) {
      piiarr[i].ival = pCpuIDRecv[i];
      piiarr[i].i2val = i;
   }

   sort (piiarr,piiarr+nrecv_loc);

// Store received data

   this->nsends = nrecv_loc;

   this->snd2cpu.resize (nrecv_loc+1);
   this->ia_sends.resize (nrecv_loc+1);

   int *psnd2cpu = &this->snd2cpu[0];
   int *pia_sends = &this->ia_sends[0];

   int nz_sends = 0;

   pia_sends[0] = 0;

   for (i=0;i<nrecv_loc;i++) {
      ind = piiarr[i].i2val;
      psnd2cpu[i] = piiarr[i].ival;
      CMatrix<int,_Flt> *pA_sub = phblk_recv[ind].GetASubArr();
      int nzjaloc = pA_sub->GetNzja() / 2;
      nz_sends += nzjaloc;
      pia_sends[i+1] = nz_sends;
   }

   this->ind_sends.resize (nz_sends+1);
   this->x_send.resize (nz_sends+1);

   int *pind_sends = &this->ind_sends[0];

   nz_sends = 0;

   int jj2, ibs;

   for (i=0;i<nrecv_loc;i++) {
      ind = piiarr[i].i2val;
      CMatrix<int,_Flt> *pA_sub = phblk_recv[ind].GetASubArr();
      int nzjaloc = pA_sub->GetNzja() / 2;
      int *pjaloc = pA_sub->GetJaArr();
      for (j=0;j<nzjaloc;j++) {
         jj = pjaloc[j*2];
         jj2 = pjaloc[j*2+1];
         ibs = pibsblk_loc[jj2];
         pind_sends[nz_sends] = ibs+jj;
         nz_sends++;
      }
   }

   this->x_temp.resize (nimax+1);

};

//
// Perform multiplication by A
//========================================================================================
template <typename _Int, typename _Flt, typename _FltVect>
void CMvmPar<_Int,_Flt,_FltVect>::MvmA (const _FltVect *_x, _FltVect *_ax) 
{

// Open mvm structure

   void *pcomm_loc = this->pcomm;

   int nblks_loc            = this->nblks;
   long long *pblks_loc     = this->pblks;
   int *pblk2cpu_loc        = this->pblk2cpu;
   int ni_cpu_loc           = this->ni_cpu;
   int *pibsblk_loc         = &this->ibsblk[0];
   int nlistblk_own_loc     = this->nlistblk_own;
   int *plistblk_own_loc    = &this->listblk_own[0];
   CHMatrix<_Int,_Flt> *phmatr_loc = this->phmatr;
   int nsends_loc           = this->nsends;
   int *psnd2cpu_loc        = &this->snd2cpu[0];
   int *pia_sends_loc       = &this->ia_sends[0];
   int *pind_sends_loc      = &this->ind_sends[0];
   _FltVect *px_send_loc    = &this->x_send[0];
   int nrecvs_loc           = this->nrecvs;
   int *prcv2cpu_loc        = &this->rcv2cpu[0];
   int *pia_recvs_loc       = &this->ia_recvs[0];
   int *pja_recvs_loc       = &this->ja_recvs[0];
   _FltVect *px_recv_loc    = &this->x_recv[0];
   int nblks_recvs_loc      = this->nblks_recvs;
   int *plistblk_recvs_loc  = &this->listblk_recvs[0];
   int *pialist_recvs_loc   = &this->ialist_recvs[0];
   int *pjapairs_recvs_loc  = &this->japairs_recvs[0];
   int *piablk_recvs_loc    = &this->iablk_recvs[0];
   int *pind_recvs_loc      = &this->ind_recvs[0];
   _FltVect *px_temp_loc    = &this->x_temp[0];

   int nproc = CExchange::GetNproc (pcomm_loc);
   int myid = CExchange::GetMyid (pcomm_loc);

// Init array ax by zeroes

   CVect<_FltVect>::SetByZeroes (ni_cpu_loc, _ax);

// Prepare send

   int ni_send_loc = pia_sends_loc[nsends_loc];

   int i, ind;

   for (i=0;i<ni_send_loc;i++) {
      ind = pind_sends_loc[i];
      px_send_loc[i] = _x[ind];
   }

   int ni_recv_loc = pia_recvs_loc[nrecvs_loc];

   CVect<_FltVect>::SetByZeroes (ni_recv_loc, px_recv_loc);

// Init async recvs and sends

   void *psndrcv_recvs_loc;
   void *psndrcv_stats_loc;

   CExchange::AllocateRecvs (nrecvs_loc+nsends_loc, psndrcv_recvs_loc);
   CExchange::AllocateStats (nrecvs_loc+nsends_loc, psndrcv_stats_loc);

   int icpu, isize, ibs;

   for (i=0;i<nrecvs_loc;i++) {
      icpu = prcv2cpu_loc[i];
      isize = (pia_recvs_loc[i+1]-pia_recvs_loc[i])*sizeof (_FltVect);
      ibs = pia_recvs_loc[i];
      CExchange::IRecv (pcomm_loc, icpu, icpu, isize, (char *)(px_recv_loc+ibs), i, psndrcv_recvs_loc);
   }

   for (i=0;i<nsends_loc;i++) {
      icpu = psnd2cpu_loc[i];
      isize = (pia_sends_loc[i+1]-pia_sends_loc[i])*sizeof (_FltVect);
      ibs = pia_sends_loc[i];
      CExchange::ISend (pcomm_loc, icpu, myid, isize, (char *)(px_send_loc+ibs), i+nrecvs_loc, psndrcv_recvs_loc);
   }

// Perform local multiplicaions

   int iblk, ibs_i, j, jblk, ibs_j;

   for (i=0;i<nlistblk_own_loc;i++) {
      iblk = plistblk_own_loc[i];
      ibs_i = pibsblk_loc[iblk];
      CMatrix<_Int,_Flt> *pasub_loc = phmatr_loc[iblk].GetASubArr();
      CMatrix<int,float> *phmatr_str_loc = phmatr_loc[iblk].GetHMatrStr();
      int nzja_hblk = phmatr_str_loc->GetNzja();
      int *pja_hblk = phmatr_str_loc->GetJaArr();
      for (j=0;j<nzja_hblk;j++) {
         jblk = pja_hblk[j];
         if (pblk2cpu_loc[jblk] == myid) {
            ibs_j = pibsblk_loc[jblk];
            CMvmSlv<_Int,_Flt,_FltVect>::MvmA (pasub_loc[j], _x+ibs_j, _ax+ibs_i);
         }
      }
   }

// Wait for completetion of sends/recvs

   CExchange::WaitAll (nrecvs_loc+nsends_loc, psndrcv_recvs_loc, psndrcv_stats_loc);

// Perform remaining multiplications

   int jj, jj2;

   for (i=0;i<nblks_recvs_loc;i++) {
      jblk = plistblk_recvs_loc[i];
      for (j=piablk_recvs_loc[i];j<piablk_recvs_loc[i+1];j++) {
         ind = pind_recvs_loc[j];
         jj = pja_recvs_loc[ind*2];
         jj2 = pja_recvs_loc[ind*2+1];
         if (jj2 != jblk) {
            throw " CTMvmSlvPar<_Int,_Flt,_FltVect>::MvmA: error: incorrect block number !";
         }
         px_temp_loc[jj] = px_recv_loc[ind];
      }
      for (j=pialist_recvs_loc[i];j<pialist_recvs_loc[i+1];j++) {
         iblk = pjapairs_recvs_loc[j*2];
         ind = pjapairs_recvs_loc[j*2+1];
         ibs_i = pibsblk_loc[iblk];
         CMatrix<_Int,_Flt> *pasub_loc = phmatr_loc[iblk].GetASubArr();
         CMatrix<int,float> *phmatr_str_loc = phmatr_loc[iblk].GetHMatrStr();
         int *pja_hblk = phmatr_str_loc->GetJaArr();
         if (pja_hblk[ind] != jblk) {
            throw " CTMvmSlvPar<_Int,_Flt,_FltVect>::MvmA: error 2: incorrect block number !";
         }
         CMvmSlv<_Int,_Flt,_FltVect>::MvmA (pasub_loc[ind], px_temp_loc, _ax+ibs_i);
      }
   }

   CExchange::DeleteRecvs (psndrcv_recvs_loc);
   CExchange::DeleteStats (psndrcv_stats_loc);

};

// Clean MvmA structure
//========================================================================================
template <typename _Int, typename _Flt, typename _FltVect>
void CMvmPar<_Int,_Flt,_FltVect>::Clean () 
{

   pcomm = NULL;
   nblks = 0;
   pblks = NULL;
   pblk2cpu = NULL;
   ni_cpu = 0;
   vector<int> ibsblk_dummy; 
   ibsblk.swap (ibsblk_dummy);
   nlistblk_own = 0;
   vector<int> listblk_own_dummy; 
   listblk_own.swap (listblk_own_dummy);
   phmatr = NULL;
   nsends = 0;
   vector<int> snd2cpu_dummy;
   vector<int> ia_sends_dummy;
   vector<int> ind_sends_dummy;
   vector<_FltVect> x_send_dummy;
   snd2cpu.swap (snd2cpu_dummy);
   ia_sends.swap (ia_sends_dummy);
   ind_sends.swap (ind_sends_dummy);
   x_send.swap (x_send_dummy);
   nrecvs = 0;
   vector<int> rcv2cpu_dummy;
   vector<int> ia_recvs_dummy;
   vector<int> ja_recvs_dummy;
   vector<_FltVect> x_recv_dummy;
   rcv2cpu.swap (rcv2cpu_dummy);
   ia_recvs.swap (ia_recvs_dummy);
   ja_recvs.swap (ja_recvs_dummy);
   x_recv.swap (x_recv_dummy);
   nblks_recvs = 0;
   vector<int> listblk_recvs_dummy;
   vector<int> ialist_recvs_dummy;
   vector<int> japairs_recvs_dummy;
   vector<int> iablk_recvs_dummy;
   vector<int> ind_recvs_dummy;
   vector<_FltVect> x_temp_dummy;
   listblk_recvs.swap (listblk_recvs_dummy);
   ialist_recvs.swap (ialist_recvs_dummy);
   japairs_recvs.swap (japairs_recvs_dummy);
   iablk_recvs.swap (iablk_recvs_dummy);
   ind_recvs.swap (ind_recvs_dummy);
   x_temp.swap (x_temp_dummy);

};

// Init SolveLU data
//========================================================================================
template <typename _Int, typename _Flt, typename _FltVect>
void CSlvPar<_Int,_Flt,_FltVect>::InitSolveLU (vector<int> *_listpairs_ext, 
                                                CMatrix<_Int,_Flt> *_matrL, CMatrix<_Int,_Flt> *_matrU, 
                                                vector<int> *_orderLU) 
{

// Get control data

   void *pcomm_loc = this->pcomm;

   int nblks_loc            = this->nblks;
   long long *pblks_loc     = this->pblks;
   long long *pblks_ext_loc = this->pblks_ext;
   int *pblk2cpu_loc        = this->pblk2cpu;
   int ni_cpu_loc           = this->ni_cpu;
   int *pibsblk_loc         = &this->ibsblk[0];
   int nlistblk_own_loc     = this->nlistblk_own;
   int *plistblk_own_loc    = &this->listblk_own[0];

   int nproc = CExchange::GetNproc (pcomm_loc);
   int myid = CExchange::GetMyid (pcomm_loc);

// Compute maximal block size

   int nimax_ext = 0;

   int i, niloc;

   for (i=0;i<nblks_loc;i++) {
      niloc = (int)(pblks_ext_loc[i+1]-pblks_ext_loc[i]);
      if (niloc > nimax_ext) nimax_ext = niloc;
   }

// Store pointers to data

   this->plistpairs_ext = _listpairs_ext;
   this->pmatrL = _matrL;
   this->pmatrU = _matrU;
   this->porderLU = _orderLU;

// Create cpu recv data

   vector<int> imaskcpu (nproc);
   vector<int> listcpu (nproc);
   vector<int> indcpu (nproc);

   int *pimaskcpu = &imaskcpu[0];
   int *plistcpu = &listcpu[0];
   int *pindcpu = &indcpu[0];

   for (i=0;i<nproc;i++) pimaskcpu[i] = -1;

   int icyclecpu = -1;

   icyclecpu++;

   int nlistcpu = 0;

   int jcpu, iblk, jj2, ni_ext, j;

   for (i=0;i<nlistblk_own_loc;i++) {
      iblk = plistblk_own_loc[i];
      ni_ext = (int)((pblks_ext_loc[iblk+1]-pblks_ext_loc[iblk])-(pblks_loc[iblk+1]-pblks_loc[iblk]));
      int *pplistpairs_ext = &(_listpairs_ext[iblk][0]);
      for (j=0;j<ni_ext;j++) {
         jj2 = pplistpairs_ext[2*j+1];
         jcpu = pblk2cpu_loc[jj2];
         if (jcpu != myid) {
            if (pimaskcpu[jcpu] != icyclecpu) {
               plistcpu[nlistcpu] = jcpu;
               nlistcpu++;
               pimaskcpu[jcpu] = icyclecpu;
            }
         }
      }
   }

   sort (plistcpu,plistcpu+nlistcpu);

   for (i=0;i<nlistcpu;i++) {
      jcpu = plistcpu[i];
      pindcpu[jcpu] = i;
   }

   this->nrecvs = nlistcpu;

   this->rcv2cpu.resize (nlistcpu+1);
   this->ia_recvs.resize (nlistcpu+1);

   int *prcv2cpu = &this->rcv2cpu[0];
   int *pia_recvs = &this->ia_recvs[0];

   for (i=0;i<nlistcpu;i++) prcv2cpu[i] = plistcpu[i];
   for (i=0;i<=nlistcpu;i++) pia_recvs[i] = 0;

   int ind;

   for (i=0;i<nlistblk_own_loc;i++) {
      iblk = plistblk_own_loc[i];
      ni_ext = (int)((pblks_ext_loc[iblk+1]-pblks_ext_loc[iblk])-(pblks_loc[iblk+1]-pblks_loc[iblk]));
      int *pplistpairs_ext = &(_listpairs_ext[iblk][0]);
      for (j=0;j<ni_ext;j++) {
         jj2 = pplistpairs_ext[2*j+1];
         jcpu = pblk2cpu_loc[jj2];
         if (jcpu != myid) {
            ind = pindcpu[jcpu];
            pia_recvs[ind+1]++;
         }
      }
   }

   for (i=0;i<nlistcpu;i++) pia_recvs[i+1] = pia_recvs[i]+pia_recvs[i+1];

   int nz_recvs = pia_recvs[nlistcpu];

   this->x_recv.resize (nz_recvs+1);

   this->ialist_recvs.resize (nlistblk_own_loc+1);
   this->jalist_recvs.resize (2*nz_recvs+1);

   int *pialist_recvs = &this->ialist_recvs[0];
   int *pjalist_recvs = &this->jalist_recvs[0];

   nz_recvs = 0;

   for (i=0;i<nlistcpu;i++) plistcpu[i] = pia_recvs[i];

   pialist_recvs[0] = 0;

   int k;

   for (i=0;i<nlistblk_own_loc;i++) {
      iblk = plistblk_own_loc[i];
      ni_ext = (int)((pblks_ext_loc[iblk+1]-pblks_ext_loc[iblk])-(pblks_loc[iblk+1]-pblks_loc[iblk]));
      int *pplistpairs_ext = &(_listpairs_ext[iblk][0]);
      for (j=0;j<ni_ext;j++) {
         jj2 = pplistpairs_ext[2*j+1];
         jcpu = pblk2cpu_loc[jj2];
         if (jcpu != myid) {
            ind = pindcpu[jcpu];
            k = plistcpu[ind];
            pjalist_recvs[nz_recvs*2] = k;
            pjalist_recvs[nz_recvs*2+1] = j;
            nz_recvs++;
            plistcpu[ind]++;
         }
      }
      pialist_recvs[i+1] = nz_recvs;
   }

// Prepare exchange data

   vector<CHMatrix<int,_Flt> > hblk_send (nlistcpu+1);

   CHMatrix<int,_Flt> *phblk_send = &hblk_send[0];

   for (i=0;i<nlistcpu;i++) {
      phblk_send[i].SetNzblk (1);
      phblk_send[i].ResizeASub (1);
      CMatrix<int,_Flt> *pA_sub = phblk_send[i].GetASubArr();
      CMatrix<int,_Flt> ablk_temp;
      int niloc = pia_recvs[i+1]-pia_recvs[i];
      ablk_temp.ResizeAndSetAllSp (0, 0, niloc*2, 0);
      pA_sub->ReplaceFree (ablk_temp);
   }

   nz_recvs = 0;

   for (i=0;i<nlistcpu;i++) plistcpu[i] = pia_recvs[i];

   int jj;

   for (i=0;i<nlistblk_own_loc;i++) {
      iblk = plistblk_own_loc[i];
      ni_ext = (int)((pblks_ext_loc[iblk+1]-pblks_ext_loc[iblk])-(pblks_loc[iblk+1]-pblks_loc[iblk]));
      int *pplistpairs_ext = &(_listpairs_ext[iblk][0]);
      for (j=0;j<ni_ext;j++) {
         jj = pplistpairs_ext[2*j];
         jj2 = pplistpairs_ext[2*j+1];
         jcpu = pblk2cpu_loc[jj2];
         if (jcpu != myid) {
            ind = pindcpu[jcpu];
            CMatrix<int,_Flt> *pA_sub = phblk_send[ind].GetASubArr();
            int *pjaloc = pA_sub->GetJaArr ();
            k = plistcpu[ind]-pia_recvs[ind];
            pjaloc[k*2] = jj;
            pjaloc[k*2+1] = jj2;
            plistcpu[ind]++;
         }
      }
   }

// Pack send data

   vector<int> CpuIDSend (nlistcpu);
   vector<vector<char> > ObjSend (nlistcpu);

   int *pCpuIDSend = NULL;
   vector<char> *pObjSend = NULL;

   if (nlistcpu > 0) {
      pCpuIDSend = &CpuIDSend[0];
      pObjSend = &ObjSend[0];
   }

   long long isize;
   char *pobj;

   for (i=0;i<nlistcpu;i++) {
      pCpuIDSend[i] = prcv2cpu[i];
      isize = phblk_send[i].GetPackedSize();
      pObjSend[i].resize ((size_t)isize);
      pobj = &(pObjSend[i][0]);
      phblk_send[i].FillPacked (isize, pobj);
      phblk_send[i].Clean ();
   }

// Exchange

   vector<int> CpuIDRecv;
   vector<vector<char> > ObjRecv;

   CExchange::DataExchange (pcomm_loc, CpuIDSend, ObjSend, CpuIDRecv, ObjRecv);

// Free send data

   {
      vector<int> CpuIDSend_temp;
      vector<vector<char> > ObjSend_temp;
      CpuIDSend.swap (CpuIDSend_temp);
      ObjSend.swap (ObjSend_temp);
   }

// Unpack receive data

   int nrecv_loc = (int) CpuIDRecv.size();

   vector<char> *pObjRecv = NULL;
   int *pCpuIDRecv = NULL;

   if (nrecv_loc > 0) {
      pObjRecv = &ObjRecv[0];
      pCpuIDRecv = &CpuIDRecv[0];
   }

   vector<CHMatrix<int,_Flt> > hblk_recv (nrecv_loc+1);

   CHMatrix<int,_Flt> *phblk_recv = &hblk_recv[0];

   for (i=0;i<nrecv_loc;i++) {
      isize = (long long) pObjRecv[i].size();
      pobj = &(pObjRecv[i][0]);
      phblk_recv[i].UnPack (isize, pobj);
   }

// Free recv data

   {
      vector<vector<char> > ObjRecv_temp;
      ObjRecv.swap (ObjRecv_temp);
   }

// Compute correct ordering of cpu data

   vector<CSortInt> iiarr (nrecv_loc+1);
   CSortInt *piiarr = &iiarr[0];

   for (i=0;i<nrecv_loc;i++) {
      piiarr[i].ival = pCpuIDRecv[i];
      piiarr[i].i2val = i;
   }

   sort (piiarr,piiarr+nrecv_loc);

// Store received data

   this->nsends = nrecv_loc;

   this->snd2cpu.resize (nrecv_loc+1);
   this->ia_sends.resize (nrecv_loc+1);

   int *psnd2cpu = &this->snd2cpu[0];
   int *pia_sends = &this->ia_sends[0];

   int nz_sends = 0;

   pia_sends[0] = 0;

   for (i=0;i<nrecv_loc;i++) {
      ind = piiarr[i].i2val;
      psnd2cpu[i] = piiarr[i].ival;
      CMatrix<int,_Flt> *pA_sub = phblk_recv[ind].GetASubArr();
      int nzjaloc = pA_sub->GetNzja() / 2;
      nz_sends += nzjaloc;
      pia_sends[i+1] = nz_sends;
   }

   this->ind_sends.resize (nz_sends+1);
   this->x_send.resize (nz_sends+1);

   int *pind_sends = &this->ind_sends[0];

   nz_sends = 0;

   int ibs;

   for (i=0;i<nrecv_loc;i++) {
      ind = piiarr[i].i2val;
      CMatrix<int,_Flt> *pA_sub = phblk_recv[ind].GetASubArr();
      int nzjaloc = pA_sub->GetNzja() / 2;
      int *pjaloc = pA_sub->GetJaArr();
      for (j=0;j<nzjaloc;j++) {
         jj = pjaloc[j*2];
         jj2 = pjaloc[j*2+1];
         ibs = pibsblk_loc[jj2];
         pind_sends[nz_sends] = ibs+jj;
         nz_sends++;
      }
   }

   this->x_temp.resize (nimax_ext*2+1);

};

// Perform SolveLU computations
//========================================================================================
template <typename _Int, typename _Flt, typename _FltVect>
void CSlvPar<_Int,_Flt,_FltVect>::SolveLU (const _FltVect *_x, _FltVect *_px) 
{

// Open mvm structure

   void *pcomm_loc = this->pcomm;

   int nblks_loc            = this->nblks;
   long long *pblks_loc     = this->pblks;
   long long *pblks_ext_loc = this->pblks_ext;
   int *pblk2cpu_loc        = this->pblk2cpu;
   int ni_cpu_loc           = this->ni_cpu;
   int *pibsblk_loc         = &this->ibsblk[0];
   int nlistblk_own_loc     = this->nlistblk_own;
   int *plistblk_own_loc    = &this->listblk_own[0];
   vector<int> *plistpairs_ext_loc = this->plistpairs_ext;
   CMatrix<_Int,_Flt> *pmatrL_loc  = this->pmatrL;
   CMatrix<_Int,_Flt> *pmatrU_loc  = this->pmatrU;
   vector<int> *porderLU_loc       = this->porderLU;
   int nsends_loc           = this->nsends;
   int *psnd2cpu_loc        = &this->snd2cpu[0];
   int *pia_sends_loc       = &this->ia_sends[0];
   int *pind_sends_loc      = &this->ind_sends[0];
   _FltVect *px_send_loc    = &this->x_send[0];
   int nrecvs_loc           = this->nrecvs;
   int *prcv2cpu_loc        = &this->rcv2cpu[0];
   int *pia_recvs_loc       = &this->ia_recvs[0];
   _FltVect *px_recv_loc    = &this->x_recv[0];
   int *pialist_recvs_loc   = &this->ialist_recvs[0];
   int *pjalist_recvs_loc   = &this->jalist_recvs[0];
   _FltVect *px_temp_loc    = &this->x_temp[0];

   int nproc = CExchange::GetNproc (pcomm_loc);
   int myid = CExchange::GetMyid (pcomm_loc);

// Init array ax by zeroes

   CVect<_FltVect>::SetByZeroes (ni_cpu_loc, _px);

// Prepare send

   int ni_send_loc = pia_sends_loc[nsends_loc];

   int i, ind;

   for (i=0;i<ni_send_loc;i++) {
      ind = pind_sends_loc[i];
      px_send_loc[i] = _x[ind];
   }

   int ni_recv_loc = pia_recvs_loc[nrecvs_loc];

   CVect<_FltVect>::SetByZeroes (ni_recv_loc, px_recv_loc);

// Init async recvs and sends

   void *psndrcv_recvs_loc;
   void *psndrcv_stats_loc;

   CExchange::AllocateRecvs (nrecvs_loc+nsends_loc, psndrcv_recvs_loc);
   CExchange::AllocateStats (nrecvs_loc+nsends_loc, psndrcv_stats_loc);

   int icpu, isize, ibs;

   for (i=0;i<nrecvs_loc;i++) {
      icpu = prcv2cpu_loc[i];
      isize = (pia_recvs_loc[i+1]-pia_recvs_loc[i])*sizeof (_FltVect);
      ibs = pia_recvs_loc[i];
      CExchange::IRecv (pcomm_loc, icpu, myid, isize, (char *)(px_recv_loc+ibs), i, psndrcv_recvs_loc);
   }

   for (i=0;i<nsends_loc;i++) {
      icpu = psnd2cpu_loc[i];
      isize = (pia_sends_loc[i+1]-pia_sends_loc[i])*sizeof (_FltVect);
      ibs = pia_sends_loc[i];
      CExchange::ISend (pcomm_loc, icpu, icpu, isize, (char *)(px_send_loc+ibs), i+nrecvs_loc, psndrcv_recvs_loc);
   }

// Wait for completetion of sends/recvs

   CExchange::WaitAll (nrecvs_loc+nsends_loc, psndrcv_recvs_loc, psndrcv_stats_loc);

   CVect<_FltVect>::SetByZeroes (ni_send_loc, px_send_loc);

// Perform local computations

   int iblk, ibs_i, j, jj, jj2, ibs_j, ind1, niloc, niextloc, ni_ini;

   _FltVect *px1_temp;
   _FltVect *px2_temp;

   for (i=0;i<nlistblk_own_loc;i++) {
      iblk = plistblk_own_loc[i];
      ibs_i = pibsblk_loc[iblk];
      niloc = (int)(pblks_loc[iblk+1]-pblks_loc[iblk]);
      niextloc = (int)(pblks_ext_loc[iblk+1]-pblks_ext_loc[iblk]);
      ni_ini = niextloc-niloc;
      px1_temp = px_temp_loc;
      px2_temp = px1_temp+niextloc;
      CVect<_FltVect>::SetByZeroes (ni_ini, px1_temp);
      memcpy (px1_temp+ni_ini,_x+ibs_i,(size_t)(niloc*sizeof(_FltVect)));
      int *pplistpairs_loc = &(plistpairs_ext_loc[iblk][0]);
      for (j=0;j<ni_ini;j++) {
         jj = pplistpairs_loc[j*2];
         jj2 = pplistpairs_loc[j*2+1];
         if (pblk2cpu_loc[jj2] == myid) {
            ibs_j = pibsblk_loc[jj2];
            px1_temp[j] = _x[ibs_j+jj];
         }
      }
      for (j=pialist_recvs_loc[i];j<pialist_recvs_loc[i+1];j++) {
         ind = pjalist_recvs_loc[2*j];
         ind1 = pjalist_recvs_loc[2*j+1];
         px1_temp[ind1] = px_recv_loc[ind];
      }
      if (porderLU != NULL) {
         CVect<_FltVect>::OrderVector    (niextloc, porderLU[iblk], px1_temp, px2_temp);
         CMvmSlv<_Int,_Flt,_FltVect>::SolveL (pmatrL_loc[iblk], px2_temp, px1_temp);
         CVect<_FltVect>::SetByZeroes (ni_ini, px1_temp);
         CMvmSlv<_Int,_Flt,_FltVect>::SolveU (pmatrU_loc[iblk], px1_temp, px2_temp);
         CVect<_FltVect>::InvOrderVector (niextloc, porderLU[iblk], px2_temp, px1_temp);
      } else {
         CMvmSlv<_Int,_Flt,_FltVect>::SolveL (pmatrL_loc[iblk], px1_temp, px2_temp);
         CVect<_FltVect>::SetByZeroes (ni_ini, px2_temp);
         CMvmSlv<_Int,_Flt,_FltVect>::SolveU (pmatrU_loc[iblk], px2_temp, px1_temp);
      }
      for (j=pialist_recvs_loc[i];j<pialist_recvs_loc[i+1];j++) {
         ind = pjalist_recvs_loc[2*j];
         ind1 = pjalist_recvs_loc[2*j+1];
         px_recv_loc[ind] = px1_temp[ind1];
      }
      for (j=0;j<ni_ini;j++) {
         jj = pplistpairs_loc[j*2];
         jj2 = pplistpairs_loc[j*2+1];
         if (pblk2cpu_loc[jj2] == myid) {
            ibs_j = pibsblk_loc[jj2];
            _px[ibs_j+jj] += px1_temp[j];
         }
      }
      CVect<_FltVect>::AddReplaceVector (niloc, px1_temp+ni_ini, _px+ibs_i);
   }

// Backward exchanges

   for (i=0;i<nsends_loc;i++) {
      icpu = psnd2cpu_loc[i];
      isize = (pia_sends_loc[i+1]-pia_sends_loc[i])*sizeof (_FltVect);
      ibs = pia_sends_loc[i];
      CExchange::IRecv (pcomm_loc, icpu, myid, isize, (char *)(px_send_loc+ibs), i+nrecvs_loc, psndrcv_recvs_loc);
   }

   for (i=0;i<nrecvs_loc;i++) {
      icpu = prcv2cpu_loc[i];
      isize = (pia_recvs_loc[i+1]-pia_recvs_loc[i])*sizeof (_FltVect);
      ibs = pia_recvs_loc[i];
      CExchange::ISend (pcomm_loc, icpu, icpu, isize, (char *)(px_recv_loc+ibs), i, psndrcv_recvs_loc);
   }

// Wait for completetion of sends/recvs

   CExchange::WaitAll (nrecvs_loc+nsends_loc, psndrcv_recvs_loc, psndrcv_stats_loc);

// Finally update received data

   for (i=0;i<ni_send_loc;i++) {
      ind = pind_sends_loc[i];
      _px[ind] += px_send_loc[i];
   }

   CExchange::DeleteRecvs (psndrcv_recvs_loc);
   CExchange::DeleteStats (psndrcv_stats_loc);

};

// Clean Slv structure
//========================================================================================
template <typename _Int, typename _Flt, typename _FltVect>
void CSlvPar<_Int,_Flt,_FltVect>::Clean () 
{

   pcomm = NULL;
   nblks = 0;
   pblks = NULL;
   pblks_ext = NULL;
   pblk2cpu = NULL;
   ni_cpu = 0;
   vector<int> ibsblk_dummy; 
   ibsblk.swap (ibsblk_dummy);
   nlistblk_own = 0;
   vector<int> listblk_own_dummy; 
   listblk_own.swap (listblk_own_dummy);
   pmatrL = NULL;
   pmatrU = NULL;
   porderLU = NULL;
   nsends = 0;
   vector<int> snd2cpu_dummy;
   vector<int> ia_sends_dummy;
   vector<int> ind_sends_dummy;
   vector<_FltVect> x_send_dummy;
   snd2cpu.swap (snd2cpu_dummy);
   ia_sends.swap (ia_sends_dummy);
   ind_sends.swap (ind_sends_dummy);
   x_send.swap (x_send_dummy);
   nrecvs = 0;
   vector<int> rcv2cpu_dummy;
   vector<int> ia_recvs_dummy;
   vector<_FltVect> x_recv_dummy;
   vector<int> ialist_recvs_dummy;
   vector<int> jalist_recvs_dummy;
   vector<_FltVect> x_temp_dummy;
   rcv2cpu.swap (rcv2cpu_dummy);
   ia_recvs.swap (ia_recvs_dummy);
   x_recv.swap (x_recv_dummy);
   ialist_recvs.swap (ialist_recvs_dummy);
   jalist_recvs.swap (jalist_recvs_dummy);
   x_temp.swap (x_temp_dummy);

}

// Prepare solver structures including performing parallel fct
//========================================================================================
template <typename _Int, typename _Flt, typename _FltVect>
void CSolver<_Int,_Flt,_FltVect>::PrepareSolver (void *_pcomm, SSolverParams *_params,
                                                   int _nblks, long long *_blks, int *_blk2cpu,
                                                   bool _b_store_matrix, _Int *_ia, _Int *_ja, _Flt *_a) 
{

// Store control data

   this->pcomm  = _pcomm;
   this->params = *_params;
   this->nblks  = _nblks;

   this->blks.resize     (_nblks+1);
   this->blks_ext.resize (_nblks+1);
   this->blk2cpu.resize  (_nblks+1);

   long long *pblks     = &this->blks[0];
   long long *pblks_ext = &this->blks_ext[0];
   int *pblk2cpu        = &this->blk2cpu[0];

   int i;

   for (i=0;i<=_nblks;i++) pblks[i] = _blks[i];
   for (i=0;i<=_nblks;i++) pblks_ext[i] = _blks[i];
   for (i=0;i<_nblks;i++) pblk2cpu[i] = _blk2cpu[i];

   int myid = CExchange::GetMyid(_pcomm);

// Allocate arrays for data

   this->hmatr_arr.resize (_nblks);
   this->nlist_ext_arr.resize (_nblks);
   this->list_ext_arr.resize (_nblks);
   this->order_LU.resize (_nblks);
   this->matrL_float.resize (_nblks);
   this->matrL_double.resize (_nblks);
   this->matrU_float.resize (_nblks);
   this->matrU_double.resize (_nblks);

   CHMatrix<_Int,_Flt> *phmatr_arr     = &this->hmatr_arr[0];
   int *pnlist_ext_arr                 = &this->nlist_ext_arr[0];
   vector<int> *plist_ext_arr          = &this->list_ext_arr[0];
   vector<int> *porder_LU              = &this->order_LU[0];
   CMatrix<_Int,float>  *pmatrL_float  = &this->matrL_float[0];
   CMatrix<_Int,double> *pmatrL_double = &this->matrL_double[0];
   CMatrix<_Int,float>  *pmatrU_float  = &this->matrU_float[0];
   CMatrix<_Int,double> *pmatrU_double = &this->matrU_double[0];

// Create matrix as set of hblocks

   int nimax = 0;

   int niloc;

   for (i=0;i<_nblks;i++) {
      niloc = (int)(pblks[i+1]-pblks[i]);
      if (niloc > nimax) nimax = niloc;
   }

   vector<_Int> listloc (nimax+1);
   vector<_Int> ialoc (nimax+1);

   _Int *plistloc = &listloc[0];
   _Int *pialoc = &ialoc[0];

   vector<int> imaskblk (_nblks);
   int *pimaskblk = &imaskblk[0];

   for (i=0;i<_nblks;i++) pimaskblk[i] = -1;

   int icycleblk = -1;

   int ibeg = 0;

   int ishift, j;

   for (i=0;i<_nblks;i++) {
      if (pblk2cpu[i] == myid) {
         niloc = (int)(pblks[i+1]-pblks[i]);
         ishift = (int)_ia[ibeg];
         for (j=0;j<niloc;j++) {
            plistloc[j] = (_Int)(j+pblks[i]);
         }
         for (j=0;j<=niloc;j++) {
            pialoc[j] = _ia[ibeg+j]-_ia[ibeg];
         }
         CHMatrix<_Int,_Flt> hblk (i, niloc, plistloc, pialoc, _ja+ishift, _a+ishift,
                                          _nblks, pblks, 
                                          icycleblk, pimaskblk);
         ibeg += niloc;
         phmatr_arr[i].ReplaceFree (hblk);
      }
   }

// Symmetrize hmatr

   vector<CHMatrix<_Int,_Flt> > hmatr_symm_arr (_nblks);
   CHMatrix<_Int,_Flt> *phmatr_symm_arr = &hmatr_symm_arr[0];

   CHMatrix<_Int,_Flt>::SymmetrizeSubmatrices (_pcomm,
                                                      _nblks, pblks, pblk2cpu, 
                                                      phmatr_arr, phmatr_symm_arr);

   for (i=0;i<_nblks;i++) {
      if (pblk2cpu[i] == myid) {
         phmatr_arr[i].ReplaceFree (phmatr_symm_arr[i]);
      }
   }

// Compute extended lists

   int ncycle_loc = _params->ncycle;

   CHMatrix<_Int,_Flt>::ExtendedLists (_pcomm, ncycle_loc,
                                             _nblks, pblks, pblk2cpu, 
                                             phmatr_arr,
                                             pnlist_ext_arr, plist_ext_arr);

// Perform filtering of the lists

   CHMatrix<_Int,_Flt>::FilterListsBack (myid, _nblks, pblk2cpu, 
                                                   pnlist_ext_arr, plist_ext_arr);

// Create extended blocks partitioning

   for (i=0;i<=_nblks;i++) pblks_ext[i] = 0;

   for (i=0;i<_nblks;i++) {
      if (pblk2cpu[i] == myid) {
         pblks_ext[i+1] = (pblks[i+1]-pblks[i])+pnlist_ext_arr[i];
      }
   }

   CExchange::ExchangeArray (_pcomm, 'L', '+',
                              _nblks+1, (void *)pblks_ext);

   for (i=0;i<_nblks;i++) pblks_ext[i+1] = pblks_ext[i]+pblks_ext[i+1];

// Get extended submatrices

   vector<CMatrix<_Int,_Flt> > matr_ext_arr (_nblks);
   CMatrix<_Int,_Flt> *pmatr_ext_arr = &matr_ext_arr[0];

   CHMatrix<_Int,_Flt>::GetExtendedSubmatrices (_pcomm,
                                                         _nblks, pblks, pblk2cpu, 
                                                         phmatr_arr,
                                                         pnlist_ext_arr, plist_ext_arr,
                                                         pmatr_ext_arr);

   int collap_loc = _params->collap;

   char strbuff[256];

   int blks_2[3];

   if (collap_loc > 0) {
      for (i=0;i<_nblks;i++) {
         if (pblk2cpu[i] == myid) {
            sprintf (strbuff,"BlkStrExt_%i.ps",i);
            blks_2[0] = 0;
            blks_2[1] = pnlist_ext_arr[i];
            blks_2[2] = pmatr_ext_arr[i].GetNlist();
            CHMatrix<_Int,_Flt>::Str2PsBox (collap_loc, pmatr_ext_arr[i], strbuff, 2, blks_2);
         }
      }
   }

   if (!_b_store_matrix) {
      this->CleanMvmA ();
   }

// Compute new ordering

   vector<int> n2_blocks (_nblks+1);

   int *pn2_blocks = &n2_blocks[0];

   for (i=0;i<_nblks;i++) pn2_blocks[i] = pmatr_ext_arr[i].GetNlist();

   int blks_3[4];

   int ordtype_loc = _params->ordtype;

   if (ordtype_loc > 0) {

      vector<CMatrix<_Int,_Flt> > matr_ord_arr (_nblks);
      CMatrix<_Int,_Flt> *pmatr_ord_arr = &matr_ord_arr[0];

      for (i=0;i<_nblks;i++) {
         if (pblk2cpu[i] == myid) {
            CIlu2_impl<_Int,_Flt>::ComputeOptimalOrderSchur (ordtype_loc, 
                                                                  pmatr_ext_arr[i].GetNlist(), pnlist_ext_arr[i], 
                                                                  *(pmatr_ext_arr[i].GetIa()), *(pmatr_ext_arr[i].GetJa()),
                                                                  pn2_blocks[i], porder_LU[i]);
            CIlu2<_Int,_Flt>::ReorderMatrix (pmatr_ext_arr[i], porder_LU[i], pmatr_ord_arr[i]);
            if (collap_loc > 0) {
               sprintf (strbuff,"BlkStrExtOrd_%i.ps",i);
               blks_3[0] = 0;
               blks_3[1] = pnlist_ext_arr[i];
               blks_3[2] = pn2_blocks[i];
               blks_3[3] = pmatr_ext_arr[i].GetNlist();
               CHMatrix<_Int,_Flt>::Str2PsBox (collap_loc, pmatr_ord_arr[i], strbuff, 3, blks_3);
            }
            pmatr_ext_arr[i].ReplaceFree (pmatr_ord_arr[i]);
         }
      }

   }

// Perform computation of the ILU decomposition for all local ordered blocks

   int prec_float_loc = _params->prec_float;

   vector<double> eigmin_arr (nblks);
   vector<double> eigmax_arr (nblks);

   double *peigmin_arr = &eigmin_arr[0];
   double *peigmax_arr = &eigmax_arr[0];

   for (i=0;i<_nblks;i++) {
      if (pblk2cpu[i] == myid) {
         if (prec_float_loc == 1) {
            CMatrix<_Int,float> *ptr_matr = NULL;
            CMatrix<_Int,float> matr_conv;
            if (sizeof(_Flt) == sizeof(float)) {
               ptr_matr = (CMatrix<_Int,float> *)(pmatr_ext_arr+i);
            } else {
               int nlist_temp = pmatr_ext_arr[i].GetNlist();
               _Int *pia_temp = pmatr_ext_arr[i].GetIaArr();
               _Int *pja_temp = pmatr_ext_arr[i].GetJaArr();
               _Flt *pa_temp = pmatr_ext_arr[i].GetAArr();
               CMatrixConv<_Int,_Flt,float>::InitAndConv (nlist_temp, pia_temp, pja_temp, pa_temp, matr_conv);
               ptr_matr = &matr_conv;
            }
            CIlu2<_Int,float>::Ilu2Matrix (*ptr_matr, *_params, 
                                                      pmatrL_float[i], pmatrU_float[i], 
                                                      peigmin_arr[i], peigmax_arr[i]);
         } else {
            CMatrix<_Int,double> *ptr_matr = NULL;
            CMatrix<_Int,double> matr_conv;
            if (sizeof(_Flt) == sizeof(double)) {
               ptr_matr = (CMatrix<_Int,double> *)(pmatr_ext_arr+i);
            } else {
               int nlist_temp = pmatr_ext_arr[i].GetNlist();
               _Int *pia_temp = pmatr_ext_arr[i].GetIaArr();
               _Int *pja_temp = pmatr_ext_arr[i].GetJaArr();
               _Flt *pa_temp = pmatr_ext_arr[i].GetAArr();
               CMatrixConv<_Int,_Flt,double>::InitAndConv (nlist_temp, pia_temp, pja_temp, pa_temp, matr_conv);
               ptr_matr = &matr_conv;
            }
            CIlu2<_Int,double>::Ilu2Matrix (*ptr_matr, *_params, 
                                                      pmatrL_double[i], pmatrU_double[i], 
                                                      peigmin_arr[i], peigmax_arr[i]);
         }
         if (collap_loc > 0) {
            sprintf (strbuff,"BlkStrU_%i.ps",i);
            blks_3[0] = 0;
            blks_3[1] = pnlist_ext_arr[i];
            blks_3[2] = pn2_blocks[i];
            blks_3[3] = pmatr_ext_arr[i].GetNlist();
            if (prec_float_loc == 1) {
               CHMatrix<_Int,float>::Str2PsBox (collap_loc, pmatrU_float[i], strbuff, 3, blks_3);
            } else {
               CHMatrix<_Int,double>::Str2PsBox (collap_loc, pmatrU_double[i], strbuff, 3, blks_3);
            }
         }
         pmatr_ext_arr[i].Clean ();
      }
   }

// Init multiplication and solve control structures

   if (_b_store_matrix) {

      this->mvm.InitControl (_pcomm, _nblks, pblks, pblk2cpu);
      this->mvm.InitMvmA (phmatr_arr);

   }

   if (prec_float_loc == 1) {

      this->slv_float.InitControl (_pcomm, _nblks, pblks, pblks_ext, pblk2cpu);
      if (ordtype_loc > 0) {
         this->slv_float.InitSolveLU (plist_ext_arr, pmatrL_float, pmatrU_float, porder_LU);
      } else {
         this->slv_float.InitSolveLU (plist_ext_arr, pmatrL_float, pmatrU_float, NULL);
      }

   } else {

      this->slv_double.InitControl (_pcomm, _nblks, pblks, pblks_ext, pblk2cpu);
      if (ordtype_loc > 0) {
         this->slv_double.InitSolveLU (plist_ext_arr, pmatrL_double, pmatrU_double, porder_LU);
      } else {
         this->slv_double.InitSolveLU (plist_ext_arr, pmatrL_double, pmatrU_double, NULL);
      }

   }

}

// Perform iterations of the iterative scheme (BiCGStab)
//========================================================================================
template <typename _Int, typename _Flt, typename _FltVect>
void CSolver<_Int,_Flt,_FltVect>::BiCGStab (int _niter_max, double _eps, int _ichk, int _msglev, ofstream *_fout, 
                                             _FltVect *_rhs, _FltVect *_sol,
                                             double &_rhs_norm, double &_res_ini, int &_niter, double &_res_fin)
{

// Init output data

   _rhs_norm = 0.0e0;
   _res_ini = 0.0e0;
   _niter = 0;
   _res_fin = 0.0e0;

// Get control data

   void *pcomm_loc      = this->GetComm();
   int nblks_loc        = this->GetNblks();
   long long *pblks_loc = this->GetBlks();
   int *pblks2cpu_loc   = this->GetBlk2cpu();

   int nproc_loc = CExchange::GetNproc (pcomm_loc);
   int myid_loc  = CExchange::GetMyid  (pcomm_loc);

// Get the local size of vector data

   int i;

   int n_local = 0;

   for (i=0;i<nblks_loc;i++) {
      if (pblks2cpu_loc[i] == myid_loc) {
         n_local += (int)(pblks_loc[i+1]-pblks_loc[i]);
      }
   }

// Allocate work vector arrays

   vector<_FltVect> rhs (n_local+1);
   vector<_FltVect> r (n_local+1);
   vector<_FltVect> p (n_local+1);
   vector<_FltVect> x (n_local+1);
   vector<_FltVect> u (n_local+1);
   vector<_FltVect> z (n_local+1);
   vector<_FltVect> v (n_local+1);
   vector<_FltVect> q (n_local+1);

   _FltVect *prhs = &rhs[0];
   _FltVect *pr = &r[0];
   _FltVect *pp = &p[0];
   _FltVect *px = &x[0];
   _FltVect *pu = &u[0];
   _FltVect *pz = &z[0];
   _FltVect *pv = &v[0];
   _FltVect *pq = &q[0];

// Compute initial residual vector and its norm

   CVect<_FltVect>::CopyVector (n_local, _sol, px);
   CVect<_FltVect>::CopyVector (n_local, _rhs, pr);

   this->MvmA (px, pz);

   CVect<_FltVect>::SubtractReplaceVector (n_local, pz, pr);

   _FltVect rhs_norm = 0.0e0;
   _FltVect resi0_norm = 0.0e0;

   rhs_norm  = CVect<_FltVect>::ScProd (n_local, _rhs,   _rhs);
   resi0_norm = CVect<_FltVect>::ScProd (n_local, pr, pr);

   double sum2_arr[2];

   sum2_arr[0] = (double)rhs_norm;
   sum2_arr[1] = (double)resi0_norm;

   CExchange::ExchangeArray (pcomm_loc, 'D', '+', 2, sum2_arr);

   double d_rhs_norm = sum2_arr[0];
   double d_resi0_norm = sum2_arr[1];

   d_rhs_norm = sqrt(d_rhs_norm);
   d_resi0_norm = sqrt(d_resi0_norm);

   _rhs_norm = d_rhs_norm;
   _res_ini = d_resi0_norm;

   if (_res_ini < _eps * _rhs_norm) {
      _res_fin = _res_ini;
      _niter = 0;
      return;
   }

   if (_msglev >= 2) std::cout << " Log10 || Rhs || = " << log10(_rhs_norm) << std::endl;
   if (_msglev >= 1 && _fout != NULL) *_fout   << " Log10 || Rhs || = " << log10(_rhs_norm) << std::endl;

   if (_msglev >= 2) std::cout << " Initial Log10 || Resi || = " << log10(d_resi0_norm) << std::endl;
   if (_msglev >= 1 && _fout != NULL) *_fout   << " Initial Log10 || Resi || = " << log10(d_resi0_norm) << std::endl;

// Perform iterations starting from residual data

   int niter_perf = 0;
   double resi_norm = d_resi0_norm;
   double d_res_min = resi_norm;

// Choose initial direction vector

   CVect<_FltVect>::CopyVector (n_local, pr, pz); // z=r(0)

// p(1) = r(0)

   CVect<_FltVect>::CopyVector (n_local, pr, pp); // p=r(0)

// Main iterative cycle

   _FltVect alpha = 0.0e0, beta, omega = 0.0e0, rho, rhoold; // method parameters
   _FltVect uu, ur, zv;                                      // inner products

   _FltVect fzero = (_FltVect) 0.0e0;
   _FltVect fone = (_FltVect) 1.0e0;

   rho = fone;

   int k, it;

   _FltVect faux;
   double daux;

   double d_resi;

   bool conv = false;

   for (k=1;k<=_niter_max;k++) {

      it = k-1; // the number of MVM is equal to 2*it

// rho(i-1) = z^T * r(i-1)

      rhoold = rho;

      rho = CVect<_FltVect>::ScProd (n_local, pz, pr);

      daux = (double) rho;

      CExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &daux);

      rho = (_FltVect) daux;

      if (k == 1) {
         alpha = omega = fone;
      } else {

// beta(i-1) = (rho(i-1) / rho(i-2)) * (alpha(i-1) / omega(i-1))

         beta = (rho / rhoold) * (alpha / omega);

// p(i) = r(i-1) + beta(i-1) * (p(i-1) - omega(i-1) * v(i-1))

         CVect<_FltVect>::UpdateVectorMinus (n_local, &omega, pv, pp);
         CVect<_FltVect>::UpdateVectorReversed (n_local, &beta, pr, pp);

      }

// q = M_solve * p(i)

      this->SlvLU (pp, pq);

// v(i) = A * q

      this->MvmA (pq, pv);

// alpha(i) = rho(i-1) / z^T * v(i)

      zv = CVect<_FltVect>::ScProd (n_local, pz, pv);

      daux = (double) zv;

      CExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &daux);

      zv = (_FltVect) daux;

      alpha = rho / zv;

// x(i-1/2) = x(i-1) + alpha(i) * q

      CVect<_FltVect>::UpdateVector (n_local, &alpha, pq, px);

// r(i-1/2) = r(i-1) - alpha(i) * v(i)

      CVect<_FltVect>::UpdateVectorMinus (n_local, &alpha, pv, pr);

// The intermediate check of convergence can be added here

// q = M_solve * r(i-1/2)

      this->SlvLU (pr, pq);

// u = A * q

      this->MvmA (pq, pu);

// omega(i) = u^T * r(i-1/2) / u^T * u

      ur = CVect<_FltVect>::ScProd (n_local, pu, pr);
      uu = CVect<_FltVect>::ScProd (n_local, pu, pu);

      sum2_arr[0] = (double)ur;
      sum2_arr[1] = (double)uu;

      CExchange::ExchangeArray (pcomm_loc, 'D', '+', 2, sum2_arr);

      ur = (_FltVect)sum2_arr[0];
      uu = (_FltVect)sum2_arr[1];

      omega = ur / uu;

// x(i) = x(i-1/2) + omega(i) * q

      CVect<_FltVect>::UpdateVector (n_local, &omega, pq, px);

// r(i) = r(i-1/2) - omega(i) * u

      CVect<_FltVect>::UpdateVectorMinus (n_local, &omega, pu, pr);

// Check the convergence

      niter_perf = (int)(it+1);

      faux = CVect<_FltVect>::ScProd (n_local, pr, pr);

      daux = (double) faux;

      CExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &daux);

      d_resi = sqrt(daux);

      if ((it+_ichk)%_ichk == 0 && _msglev >= 2) std::cout << " It = " << it << " Log10 || Resi || = " << log10(d_resi) << std::endl;
      if ((it+_ichk)%_ichk == 0 && _msglev >= 1 && _fout != NULL) *_fout   << " It = " << it << " Log10 || Resi || = " << log10(d_resi) << std::endl;

      if (d_resi < _eps * d_rhs_norm) conv = true;

// Save the best attained solution

      if (d_resi < d_res_min) {
         d_res_min = d_resi;
         CVect<_FltVect>::CopyVector (n_local, px, _sol);
      }

// Break from iterations if converged

      if (conv) {
         break;
      }

   } // end of iterations

// Compute the final residual

   this->MvmA (_sol, pr);

   CVect<_FltVect>::SubtractReplaceVector (n_local, _rhs, pr); // r=b-A*x

// Compute the final residual norm

   resi_norm = CVect<_FltVect>::ScProd (n_local, pr, pr);

   daux = (double) resi_norm;

   CExchange::ExchangeArray (pcomm_loc, 'D', '+', 1, &daux);

   double d_resi_fin = sqrt (daux);

   if (_msglev >= 2) std::cout << " Final Log10 || Resi || = " << log10(d_resi_fin) << std::endl;
   if (_msglev >= 1 && _fout != NULL) *_fout   << " Final Log10 || Resi || = " << log10(d_resi_fin) << std::endl;

// Return statistics data

   _niter = niter_perf;
   _res_fin = d_resi_fin;

}

template class CVect       <double>;
template class CVect       <float>;
template class CIlu2_impl  <long long,double>;
template class CIlu2_impl  <int,      double>;
template class CIlu2_impl  <long long,float>;
template class CIlu2_impl  <int,      float>;
template class CMatrix     <long long,double>;
template class CMatrix     <int,      double>;
template class CMatrix     <long long,float>;
template class CMatrix     <int,      float>;
template class CMvmSlv_impl<long long,double,double>;
template class CMvmSlv_impl<int,      double,double>;
template class CMvmSlv_impl<long long,float, double>;
template class CMvmSlv_impl<int,      float, double>;
template class CMvmSlv     <long long,double,double>;
template class CMvmSlv     <int,      double,double>;
template class CMvmSlv     <long long,float, double>;
template class CMvmSlv     <int,      float, double>;
template class CIlu2       <long long,double>;
template class CIlu2       <int,      double>;
template class CIlu2       <long long,float>;
template class CIlu2       <int,      float>;
template class CHMatrix    <long long,double>;
template class CHMatrix    <int,      double>;
template class CHMatrix    <long long,float>;
template class CHMatrix    <int,      float>;
template class CMatrixConv <long long,double,double>;
template class CMatrixConv <int,      double,double>;
template class CMatrixConv <long long,double,float>;
template class CMatrixConv <int,      double,float>;
template class CMvmPar     <long long,double,double>;
template class CMvmPar     <int,      double,double>;
template class CMvmPar     <long long,float, double>;
template class CMvmPar     <int,      float, double>;
template class CSlvPar     <long long,double,double>;
template class CSlvPar     <int,      double,double>;
template class CSlvPar     <long long,float, double>;
template class CSlvPar     <int,      float, double>;
template class CSolver     <long long,double,double>;
template class CSolver     <int,      double,double>;
template class CSolver     <long long,float, double>;
template class CSolver     <int,      float, double>;

////////////////////////////////////////////////////////////////////////////////////////////////
//...reordering algorithms (reverse Cuthill-Mckee ordering);
void genrcm(int n, int * xadj, int * iadj, int * perm, int * xls, int * mask);
void subrcm(int * xadj, int * iadj, int * mask, int nsubg, int * subg, int * perm, int * xls, int n);
/////////////////////////
//...parameters:
//   n          - the dimension of the matrix;
//   xadj, iadj - the matrix structure: xadj[n+1], iadj[*]; information about row i is stored 
//                in xadj[i-1] -- xadj[i]-1 of the the adjacency structure iadj[*];
//                for each row, it contains the column indices of the nonzero entries;
//   perm[n]    - contains the rcm ordering;
//   mask[n]    - marks variables that have been numbered (working array);
//   xls[n+1]   - the index vector for a level structure; the level structure is stored 
//                in the currently unused spaces in the permutation vector perm;
//   nsubg      - the size of the subgraph;
//   subg[n]    - contains the nodes in subgraph (which may be disconnected);
/////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
//...generates the connected level structure rooted at a given node;
void rootls(int root, int * xadj, int * iadj, int * mask, int & nlvl, int * xls, int * ls, int)
{
   int i, iccsze, j, jstop, jstrt, lbegin, lvlend, lvsize, nbr, node;

   mask[root-1] = 0;
   ls[0]        = root;
   nlvl         = 0;
   lvlend       = 0;
   iccsze       = 1;

   do {
       lbegin      = lvlend + 1;
       lvlend      = iccsze;
       xls[nlvl++] = lbegin;

/////////////////////////////////////////////////////////////////////////////////////////////
//...Generate the next level by finding all the masked neighbors of nodes in the current level;
       for (i = lbegin; i <= lvlend;  i++) {
            jstrt = xadj[(node = ls[i-1])-1];
            jstop = xadj[ node]-1;

            for (j = jstrt; j <= jstop; j++)
            if  (mask[(nbr = iadj[j-1])-1] != 0) {
                 ls[iccsze++] = nbr;
                 mask[nbr-1]  = 0;
            }
       }
   }
   while ((lvsize = iccsze-lvlend) > 0);

////////////////////////////////////////////////////////////
//...Reset MASK to one for the nodes in the level structure;
   for (xls[nlvl] = lvlend+1, i = 1; i <= iccsze; i++)
       mask[(node = ls[i-1])-1] = 1;
}

///////////////////////////////////
//...finds pseudo-peripheral nodes;
void fnroot(int & root, int * xadj, int * iadj, int * mask, int & nlvl, int * xls, int * ls, int n)
{
   int iccsze, j, jstrt, k, kstop, kstrt, mindeg, nabor, ndeg, node, nunlvl;

   rootls(root, xadj, iadj, mask, nlvl, xls, ls, n);
   if (nlvl == 1 || nlvl == (iccsze = xls[nlvl]-1)) return;

   do {
       mindeg = iccsze;
       root   = ls[(jstrt = xls[nlvl-1])-1];

       if (iccsze > jstrt) {
           for (j = jstrt; j <= iccsze; j++) {
                ndeg  = 0;
                kstrt = xadj[(node = ls[j-1])-1];
                kstop = xadj[ node]-1;

                for (k = kstrt; k <= kstop; k++) 
                if (mask[(nabor = iadj[k-1])-1] > 0 ) ndeg++;

                if (ndeg < mindeg) {
                    root   = node;
                    mindeg = ndeg;
                }
           }
       }

       rootls (root, xadj, iadj, mask, nunlvl, xls, ls, n);
       if (nunlvl <= nlvl) return;
   }
   while ((nlvl = nunlvl) < iccsze);
}

//////////////////////////////////////////////////////////////////
//...computes the degrees of the nodes in the connected component;
void degree (int root, int * xadj, int * iadj, int * mask, int * deg, int & iccsze, int * ls, int)
{
   int i, ideg, j, jstop, jstrt, lbegin, lvlend, lvsize, nbr, node;

   ls[0]        = root;
   lvlend       = 0;
   iccsze       = 1;
   xadj[root-1] = -xadj[root-1];

   do {
       lbegin = lvlend+1;
       lvlend = iccsze;

       for (i = lbegin; i <= lvlend; i++) {
            jstrt = -xadj[(node = ls[i-1])-1];
//          jstop =  abs(xadj[node])-1;
            jstop = xadj[node];
            if (jstop < 0) jstop = -jstop;
            jstop--;
            ideg  = 0;

            for (j = jstrt; j <= jstop; j++)
            if  (mask[(nbr = iadj[j-1])-1] != 0) {
                 ideg = ideg+1;
                 if (xadj[nbr-1] >= 0) {
                     xadj[nbr-1]  = -xadj[nbr-1];
                     ls[iccsze++] = nbr;
                 }
            }
            deg[node-1] = ideg;
       }
   }
   while ((lvsize = iccsze - lvlend) > 0);

///////////////////////////////////////////////
//...Reset XADJ to its correct sign and return;
   for (i = 1; i <= iccsze; i++) {
        node         = ls[i-1];
        xadj[node-1] = -xadj[node-1];
   }
}


////////////////////////////////////////////////
//...reverses the elements of an integer vector;
#define iSWAP(A, B) { int iSWAP_temp = (A); (A) = (B); (B) = iSWAP_temp; }

void ivec_reverse (int n, int * a)
{
  int  m, i;
  for (m = n/2, i = 1; i <= m; i++)
       iSWAP(a[i-1], a[n-i]);
}

#undef iSWAP

/////////////////////////////////////////////////////////////////////////////
//...numbers a connected component using the reverse Cuthill McKee algorithm;
void rcm(int root, int * xadj, int * iadj, int * mask, int * perm, int & iccsze, int * deg, int n)
{
   int fnbr, i, j, jstop, jstrt, k, l, lbegin, lnbr, lperm, lvlend, nbr, node;

   degree (root, xadj, iadj, mask, deg, iccsze, perm, n);
   mask[root-1] = 0;
   if ( iccsze <= 1) return;

   lvlend = 0;
   lnbr   = 1;

   do {
       lbegin = lvlend+1;
       lvlend = lnbr;

       for (i = lbegin; i <= lvlend; i++) {
            jstrt = xadj[(node = perm[i-1])-1];
            jstop = xadj[ node]-1;
            fnbr  = lnbr+1;

            for (j = jstrt; j <= jstop; j++)
            if  (mask[(nbr = iadj[j-1])-1] != 0) {
                 mask[ nbr-1] = 0;
                 perm[lnbr++] = nbr;
            }
/////////////////////////////////////////////////////////////
//...Sort the neighbors of node in increasing order by degree;
            if (fnbr < lnbr) {
                k = fnbr;
                do {
                    l   = k;
                    nbr = perm[k++];
label40:
                    if (l > fnbr && deg[(lperm = perm[l-1])-1] > deg[nbr-1]) {
                        perm[l--] = lperm;
                        goto label40;
                    }
                    perm[l] = nbr;
                }
                while (k < lnbr);
            }
       }
   }
   while (lnbr > lvlend);

////////////////////////////////////////
//...Reverse the Cuthill-McKee ordering;
   ivec_reverse(iccsze, perm);
}

//////////////////////////////////////////////////////////////////
//...finds the reverse Cuthill-Mckee ordering for a general graph;
void genrcm(int n, int * xadj, int * iadj, int * perm, int * xls, int * mask)
{
   int i, iccsze, nlvl, num, root;

   for (         i = 1; i <= n; i++) mask[i-1] = 1;
   for (num = 1, i = 1; i <= n; i++)
   if  (mask[i-1] != 0) {
       fnroot(root = i, xadj, iadj, mask, nlvl, xls,  perm+num-1,  n);
       rcm   (root,     xadj, iadj, mask, perm+num-1, iccsze, xls, n);

       if ((num += iccsze) > n) return;
   }
}

///////////////////////////////////////////////////////////////////
//...finds the reverse Cuthill-Mckee ordering for a given subgraph;
void subrcm(int * xadj, int * iadj, int * mask, int nsubg, int * subg, int * perm, int * xls, int n)
{
   int  i, iccsze, nlvl, node, num;
   for (i = 1; i <= nsubg; i++)
        mask[(node = subg[i-1])-1] = 1;

   for (num = 0, i = 1; i <= nsubg; i++)
   if  (mask[(node = subg[i-1])-1] > 0 ) {
       fnroot(node, xadj, iadj, mask, nlvl, xls,  perm+num,  n);
       rcm   (node, xadj, iadj, mask, perm+num, iccsze, xls, n);

       if ((num += iccsze) >= nsubg ) return;
   }
}

// Barrier synchronization of the CPUs
//========================================================================================
void CExchange::Synchronize (void *_comm) 
{

   MPI_Comm *pcomm = (MPI_Comm *)_comm;
   MPI_Barrier (*pcomm);

}

// Get the number of CPUs
//========================================================================================
int CExchange::GetNproc (void *_comm) 
{

   MPI_Comm *pcomm = (MPI_Comm *)_comm;

   int nproc;
   MPI_Comm_size (*pcomm,&nproc);

   return nproc;

}

// Get the ID number of the CPU
//========================================================================================
int CExchange::GetMyid (void *_comm) 
{

   MPI_Comm *pcomm = (MPI_Comm *)_comm;

   int myid;
   MPI_Comm_rank (*pcomm,&myid);

   return myid;

}

// Description: Asynchronous send data to the other CPU (implementation)
//========================================================================================
void CExchange::ISend (void *_comm, int _rank, int _msgtag, int _length, char *_buffer, int _indrecv, void *_recv_arr) 
{

   MPI_Comm *pcomm = (MPI_Comm *)_comm;

   MPI_Request *p_recv_arr = (MPI_Request *)_recv_arr;

   int ierr = MPI_Isend (_buffer, _length, MPI_CHAR, _rank, _msgtag, *(pcomm), p_recv_arr+_indrecv);

}

// Description: Asynchronous receive data from the other CPU (implementation)
//========================================================================================
void CExchange::IRecv (void *_comm, int _rank, int _msgtag, int _length, char *_buffer, int _indrecv, void *_recv_arr) 
{

   MPI_Comm *pcomm = (MPI_Comm *)_comm;

   MPI_Request *p_recv_arr = (MPI_Request *)_recv_arr;

   int ierr = MPI_Irecv (_buffer, _length, MPI_CHAR, _rank, _msgtag, *(pcomm), p_recv_arr+_indrecv);

}

// Description: Allocate send/receive requests (implementation)
//========================================================================================
void CExchange::AllocateRecvs (int _nrecv, void *&_recvarr) 
{

   MPI_Request *p_recv_arr = new MPI_Request [_nrecv];

   _recvarr = (void *)p_recv_arr;

}

// Description: Delete send/receive requests (implementation)
//========================================================================================
void CExchange::DeleteRecvs (void *_recvarr) 
{

   MPI_Request *p_recv_arr = (MPI_Request *)_recvarr;

   delete [] p_recv_arr;

}

// Description: Allocate statuses (implementation)
//========================================================================================
void CExchange::AllocateStats (int _nstat, void *&_statarr) 
{

   MPI_Status *p_stat_arr = new MPI_Status [_nstat];

   _statarr = (void *)p_stat_arr;

}

// Description: Delete statuses (implementation)
//========================================================================================
void CExchange::DeleteStats (void *_statarr) 
{

   MPI_Status *p_stat_arr = (MPI_Status *)_statarr;

   delete [] p_stat_arr;

}

// Description: Wait for completion of exchanges (implementation)
//========================================================================================
void CExchange::WaitAll (int _count, void *_recvarr, void *_statarr) 
{

   MPI_Request *p_recv_arr = (MPI_Request *)_recvarr;
   MPI_Status *p_stat_arr = (MPI_Status *)_statarr;

   int ierr = MPI_Waitall (_count, p_recv_arr, p_stat_arr);

}

//========================================================================================
void CExchange::ExchangeArray (void *_comm,
                                    char _Datatype, char _Operation,
                                    int _Length, void *_Array) 
{

   static int maxbufsize = 2408*1024;

   MPI_Comm *pcomm = (MPI_Comm *)_comm;

   if (_Length < maxbufsize) {

      CExchange::ExchangeArray_impl(_comm, 
                                          _Datatype, _Operation,
                                          _Length, _Array);

   } else {

      int myidloc;
      MPI_Comm_rank (*pcomm,&myidloc);

      char *parr = (char *)_Array;

      int iscale = 0;

      if (_Datatype == 'C') iscale = sizeof(char);
      if (_Datatype == 'I') iscale = sizeof(int);
      if (_Datatype == 'L') iscale = sizeof(long long);
      if (_Datatype == 'F') iscale = sizeof(float);
      if (_Datatype == 'D') iscale = sizeof(double);

      int ibeg, iend, ni;

      iend = -1;
      while (iend < _Length-1) {
         ibeg = iend+1;
         iend = ibeg+maxbufsize-1;
         if (iend > _Length-1) iend = _Length-1;
         ni = iend-ibeg+1;
         CExchange::ExchangeArray_impl (_comm,
                                             _Datatype, _Operation,
                                             ni, (void *)(parr+ibeg*iscale));
      }

   }

}

//========================================================================================
void OpCharAdd (void *in, void *inout, int *len, MPI_Datatype * )
{
   char *inL = (char*)in;
   char *inoutL = (char*)inout;
   int i; 
   char c; 

   for (i=0; i< *len; ++i) { 
      c = *inoutL+*inL;
      *inoutL = c; 
      inL++; inoutL++; 
   }
};

//========================================================================================
void OpCharMaximum (void *in, void *inout, int *len, MPI_Datatype *dptr )
{
   char *inL = (char*)in;
   char *inoutL = (char*)inout;
   int i; 
   char c; 

   for (i=0; i< *len; ++i) { 
      c = *inoutL>*inL ? *inoutL : *inL;
      *inoutL = c; 
      inL++; inoutL++; 
   }
};

//========================================================================================
void OpCharMinimum (void *in, void *inout, int *len, MPI_Datatype *dptr )
{
   char *inL = (char*)in;
   char *inoutL = (char*)inout;
   int i; 
   char c; 

   for (i=0; i< *len; ++i) { 
      c = *inoutL<*inL ? *inoutL : *inL;
      *inoutL = c; 
      inL++; inoutL++; 
   }
};

//========================================================================================
void OpLongLongAdd (void *in, void *inout, int *len, MPI_Datatype *dptr )
{
   long long *inL = (long long*)in;
   long long *inoutL = (long long*)inout;
   int i; 
   long long c; 

   for (i=0; i< *len; ++i) { 
      c = *inoutL+*inL;
      *inoutL = c; 
      inL++; inoutL++; 
   }
};

// Collective array operation
//========================================================================================
void CExchange::ExchangeArray_impl (void *_comm,
                                    char _Datatype, char _Operation,
                                    int _Length, void *_Array) 
{

   MPI_Comm *pcomm = (MPI_Comm *)_comm;

   int nprocloc;
   MPI_Comm_size (*pcomm,&nprocloc);

   MPI_Datatype datatype;
   MPI_Op op;

   int size_of_data = 0;

   bool op_create_flag = false;

   if (_Operation == '+') {
      op = MPI_SUM;
   } else if (_Operation == 'M') {
      op = MPI_MAX;
   } else if (_Operation == 'm') {
      op = MPI_MIN;
   }

   if (_Datatype == 'C') {

      datatype = MPI_CHAR;
      size_of_data = sizeof(char);

      if (_Operation == '+') {
         MPI_Op_create (OpCharAdd, true, &op);
         op_create_flag = true;
      } else if (_Operation == 'M') {
         MPI_Op_create (OpCharMaximum, true, &op);
         op_create_flag = true;
      } else if (_Operation == 'm') {
         MPI_Op_create (OpCharMinimum, true, &op);
         op_create_flag = true;
      }

   } else if (_Datatype == 'I') {

      datatype = MPI_INT;
      size_of_data = sizeof(int);

   } else if (_Datatype == 'L') {

      datatype = MPI_LONG_LONG_INT;
      size_of_data = sizeof(long long);
      if (_Operation == '+')
      {
         MPI_Op_create (OpLongLongAdd, true, &op);
         op_create_flag = true;
      }

   } else if (_Datatype == 'F') {

      datatype = MPI_FLOAT;
      size_of_data = sizeof(float);

   } else if (_Datatype == 'D') {

      datatype = MPI_DOUBLE;
      size_of_data = sizeof(double);

   }

   if (nprocloc != 1) {
      void* buffer = new char[_Length*size_of_data];
      MPI_Allreduce(_Array, buffer, _Length, datatype, op, *(pcomm));
      memcpy(_Array, buffer, _Length*size_of_data);
      delete [] (char *) buffer;
   }

   if (op_create_flag) {
      MPI_Op_free (&op);
   }

}

// Major data exchange function
//========================================================================================
void CExchange::DataExchange (void *_comm,
                              vector<int> &_CpuIDSend, vector<vector<char> > &_ObjSend,
                              vector<int> &_CpuIDRecv, vector<vector<char> > &_ObjRecv) 
{

   int nproc = CExchange::GetNproc(_comm);

   if (nproc == 1) {
      _CpuIDRecv.swap (_CpuIDSend);
      _ObjRecv.swap (_ObjSend);
   } else {
      CExchange::DataExchange_impl (_comm,
                                    _CpuIDSend, _ObjSend,
                                    _CpuIDRecv, _ObjRecv);
   }

}

// Description: Packed data exchange (implementation)
//========================================================================================
void CExchange::DataExchange_impl (void *_comm,
                                          vector<int> &_CpuIDSend, vector<vector<char> > &_ObjSend,
                                          vector<int> &_CpuIDRecv, vector<vector<char> > &_ObjRecv) {

// Exchange the data

   int NObjSend = (int)_CpuIDSend.size();

   vector<int> CpuIDSend   (NObjSend);
   vector<int> ObjSizeSend (NObjSend);
   vector<char *> ObjSend  (NObjSend);

   int* pCpuIDSend   = NULL;
   int* pObjSizeSend = NULL;
   char** pObjSend   = NULL;

   if (NObjSend > 0) {
      pCpuIDSend   = &CpuIDSend[0];
      pObjSizeSend = &ObjSizeSend[0];
      pObjSend     = &ObjSend[0];
   }

   int *p_data_send_cpuid = NULL;
   vector<char> *p_data_send_ch = NULL;

   if (NObjSend > 0) {
      p_data_send_cpuid = &_CpuIDSend[0];
      p_data_send_ch    = &_ObjSend[0];
   }

   int i;

   for (i=0;i<NObjSend;i++) {
      pCpuIDSend[i] = p_data_send_cpuid[i];
      pObjSizeSend[i] = (int)p_data_send_ch[i].size();
      pObjSend[i] = &(p_data_send_ch[i][0]);
   }

   int NObjRecv;
   int* CpuIDRecv;
   int* ObjSizeRecv;
   char** ObjRecv;

   CExchange::DataExchange_impl2 (_comm,
                                          (int) NObjSend, pCpuIDSend, pObjSizeSend, pObjSend,
                                          NObjRecv, CpuIDRecv, ObjSizeRecv, ObjRecv);

// Prepare reply

   _CpuIDRecv.resize (NObjRecv);
   _ObjRecv.resize (NObjRecv);

   int *p_data_recv_cpuid = NULL;
   vector<char> *p_data_recv_ch = NULL;

   if (NObjRecv > 0) {
      p_data_recv_cpuid = &_CpuIDRecv[0];
      p_data_recv_ch = &_ObjRecv[0];
   }

   char *pwork;

   for (i = 0; i < NObjRecv; i++) {
      p_data_recv_cpuid[i] = CpuIDRecv[i];
      p_data_recv_ch[i].resize (ObjSizeRecv[i]);
      pwork = NULL;
      if (ObjSizeRecv[i] != 0) {
         pwork = &(p_data_recv_ch[i][0]);
         memcpy (pwork, ObjRecv[i], ObjSizeRecv[i]);
      }
      delete [] ObjRecv[i];
   }

// Free receive structures

   delete [] CpuIDRecv;
   delete [] ObjSizeRecv;
   delete [] ObjRecv;

}

// Data exchange function
//========================================================================================
void CExchange::DataExchange_impl2 (void *_comm,
                                       int _NObjSend, int *_CpuIDSend, int *_ObjSizeSend, char **_ObjSend,
                                       int &_NObjRecv, int *&_CpuIDRecv, int *&_ObjSizeRecv, char **&_ObjRecv) 
{

// Take the number of cpus and cpu id

   MPI_Comm *pcomm = (MPI_Comm *)_comm;

   int nproc;
   int myid;

   MPI_Comm_size (*pcomm,&nproc);
   MPI_Comm_rank (*pcomm,&myid);

// Exchange send control data and transform it into receive control data

   CExchange::DataExchangeRequest (_comm,
                                          _NObjSend, _CpuIDSend, _ObjSizeSend,
                                          _NObjRecv, _CpuIDRecv, _ObjSizeRecv);

// Create structures that support sparse search

   int *iasend, *jasend;
   int *iarecv, *jarecv;
   int *iptr;

   iasend = new int [nproc+1];
   jasend = new int [_NObjSend];
   iarecv = new int [nproc+1];
   jarecv = new int [_NObjRecv];
   iptr = new int [nproc];

   int i, iproc;

   for (i=0;i<=nproc;i++) iasend[i] = 0;

   for (i=0;i<_NObjSend;i++) {
      iproc = _CpuIDSend[i];
      iasend[iproc+1]++;
   }

   for (i=0;i<nproc;i++) iasend[i+1] = iasend[i]+iasend[i+1];

   for (i=0;i<nproc;i++) iptr[i] = iasend[i];

   int k;

   for (i=0;i<_NObjSend;i++) {
      iproc = _CpuIDSend[i];
      k = iptr[iproc];
      jasend[k] = i;
      iptr[iproc]++;
   }

   for (i=0;i<=nproc;i++) iarecv[i] = 0;

   for (i=0;i<_NObjRecv;i++) {
      iproc = _CpuIDRecv[i];
      iarecv[iproc+1]++;
   }

   for (i=0;i<nproc;i++) iarecv[i+1] = iarecv[i]+iarecv[i+1];

   for (i=0;i<nproc;i++) iptr[i] = iarecv[i];

   for (i=0;i<_NObjRecv;i++) {
      iproc = _CpuIDRecv[i];
      k = iptr[iproc];
      jarecv[k] = i;
      iptr[iproc]++;
   }

// Create the set of send/receive MPI datatypes

   int *nobjsend, *nobjrecv;
   size_t *sizesendtot, *sizerecvtot;

   nobjsend = new int [nproc];
   nobjrecv = new int [nproc];
   sizesendtot = new size_t [nproc];
   sizerecvtot = new size_t [nproc];

   int j;

   for (iproc=0;iproc<nproc;iproc++) {
      nobjsend[iproc] = 0;
      sizesendtot[iproc] = 0;
      for (j=iasend[iproc];j<iasend[iproc+1];j++) {
         i = jasend[j];
         sizesendtot[iproc] += _ObjSizeSend[i];
         nobjsend[iproc]++;
      }
      nobjrecv[iproc] = 0;
      sizerecvtot[iproc] = 0;
      for (j=iarecv[iproc];j<iarecv[iproc+1];j++) {
         i = jarecv[j];
         sizerecvtot[iproc] += _ObjSizeRecv[i];
         nobjrecv[iproc]++;
      }
   }

// Allocate receive data

   _ObjRecv = new char* [_NObjRecv];

   for (i=0;i<_NObjRecv;i++) {
      _ObjRecv[i] = new char [_ObjSizeRecv[i]];
   }

// Perform local copy

   int ipsend, iprecv, nloc;

   if (nobjsend[myid] != 0) {
      ipsend = 0;
      iprecv = 0;
      nloc = 0;
      while (nloc != nobjsend[myid]) {
         while (_CpuIDSend[ipsend] != myid) ipsend++;
         while (_CpuIDRecv[iprecv] != myid) iprecv++;
         memcpy (_ObjRecv[iprecv], _ObjSend[ipsend], _ObjSizeSend[ipsend] * sizeof(char));
         ipsend++;
         iprecv++;
         nloc++;
      }
   }

// Create the set of send/receive MPI datatypes

   size_t isizesend = 0;
   size_t isizerecv = 0;

   for (iproc=0;iproc<nproc;iproc++) {
      if (iproc != myid) {
         if (sizesendtot[iproc] > isizesend) isizesend = sizesendtot[iproc];
         if (sizerecvtot[iproc] > isizerecv) isizerecv = sizerecvtot[iproc];
      }
   }

   size_t maxsendrecvbuf = 128*1024*1024;

   size_t isizemax = isizesend;
   if (isizerecv > isizemax) isizemax = isizerecv;

   size_t ncycle_send = 1;
   if (isizemax > maxsendrecvbuf) ncycle_send = (isizemax+maxsendrecvbuf-1) / maxsendrecvbuf;

   char *objsendbuf;
   char *objrecvbuf;

   objsendbuf = new char [isizesend];
   objrecvbuf = new char [isizerecv];

// Main exchange cycles

   int iproc_send, iproc_recv;

   MPI_Status stat;

   for (k=1;k<nproc;k++) {

      iproc_send = (myid+k) % nproc;
      iproc_recv = (myid-k+nproc) % nproc;

      isizesend = 0;
      for (j=iasend[iproc_send];j<iasend[iproc_send+1];j++) {
         i = jasend[j];
         memcpy (objsendbuf+isizesend, _ObjSend[i], _ObjSizeSend[i] * sizeof(char));
         isizesend += _ObjSizeSend[i];
      }
      isizerecv = 0;
      for (j=iarecv[iproc_recv];j<iarecv[iproc_recv+1];j++) {
         i = jarecv[j];
         isizerecv += _ObjSizeRecv[i];
      }

//      if (isizesend > 0 || isizerecv > 0) {

         size_t iend_send = -1;
         size_t iend_recv = -1;

         for (int icycle=0;icycle<ncycle_send;icycle++) {

            size_t ibeg_send = iend_send+1;
            size_t ibeg_recv = iend_recv+1;

            iend_send = ibeg_send+maxsendrecvbuf;
            iend_recv = ibeg_recv+maxsendrecvbuf;

            if (iend_send >= isizesend) iend_send = isizesend-1;
            if (iend_recv >= isizerecv) iend_recv = isizerecv-1;

            size_t ni_send = iend_send-ibeg_send+1;
            size_t ni_recv = iend_recv-ibeg_recv+1;

//            if (ni_send > 0 || ni_recv > 0) {

               MPI_Sendrecv (objsendbuf+ibeg_send, (int)(ni_send*sizeof(char)), MPI_CHAR, iproc_send, iproc_send, 
                                    objrecvbuf+ibeg_recv, (int)(ni_recv*sizeof(char)), MPI_CHAR, iproc_recv, myid, 
                                    *(pcomm), &stat);

//            }

         }

//      }

      isizerecv = 0;
      for (j=iarecv[iproc_recv];j<iarecv[iproc_recv+1];j++) {
         i = jarecv[j];
         memcpy (_ObjRecv[i], objrecvbuf+isizerecv, _ObjSizeRecv[i] * sizeof(char));
         isizerecv += _ObjSizeRecv[i];
      }

   }

// Free work arrays

   delete [] iasend;
   delete [] jasend;
   delete [] iarecv;
   delete [] jarecv;
   delete [] iptr;
   delete [] nobjsend;
   delete [] nobjrecv;
   delete [] sizesendtot;
   delete [] sizerecvtot;
   delete [] objsendbuf;
   delete [] objrecvbuf;

};

// Description: Exchange sizes data requests
//========================================================================================
void CExchange::DataExchangeRequest (void *_comm,
                                          int _NObjSend, int *_CpuIDSend, int *_ObjSizeSend,
                                          int &_NObjRecv, int *&_CpuIDRecv, int *&_ObjSizeRecv) {

// Take the number of cpus and cpu id

   MPI_Comm *pcomm = (MPI_Comm *)_comm;

   int nproc;
   int myid;

   MPI_Comm_size (*pcomm,&nproc);
   MPI_Comm_rank (*pcomm,&myid);

// Collect the numbers of objects from all cpu's

   int *cpu2obj;

   cpu2obj = new int [nproc+1];

   int i, iproc;

   for (i=0;i<=nproc;i++) cpu2obj[i] = 0;

   cpu2obj[myid+1] = _NObjSend;

   if (nproc > 1) {

      CExchange::ExchangeArray (_comm, 'I', '+', nproc+1, cpu2obj);

   }

   for (i=0;i<nproc;i++) cpu2obj[i+1] = cpu2obj[i] + cpu2obj[i+1];

// Allocate and combine all the data into one collected set of arrays

   int nobjtot = cpu2obj[nproc];

   int *objsizetot;
   int *cpuidtot;

   objsizetot = new int [nobjtot];
   cpuidtot = new int [nobjtot];

   for (i=0;i<nobjtot;i++) {

      cpuidtot  [i] = 0;
      objsizetot[i] = 0;

   }

   int ibeg = cpu2obj[myid];

   for (i=cpu2obj[myid];i<cpu2obj[myid+1];i++) {

      cpuidtot  [i] = _CpuIDSend  [i-ibeg];
      objsizetot[i] = _ObjSizeSend[i-ibeg];

   }

   if (nproc > 1) {

      CExchange::ExchangeArray (_comm, 'I', '+', nobjtot, cpuidtot);
      CExchange::ExchangeArray (_comm, 'I', '+', nobjtot, objsizetot);

   }

// Prepare the final data

   _NObjRecv = 0;

   for (i=0;i<nobjtot;i++) {
      if (cpuidtot[i] == myid) _NObjRecv++;
   }

   _ObjSizeRecv = new int [_NObjRecv];
   _CpuIDRecv = new int [_NObjRecv];

   _NObjRecv = 0;

   for (iproc=0;iproc<nproc;iproc++) {
      for (i=cpu2obj[iproc];i<cpu2obj[iproc+1];i++) {
         if (cpuidtot[i] == myid) {
            _ObjSizeRecv[_NObjRecv] = objsizetot[i];
            _CpuIDRecv[_NObjRecv] = iproc;
            _NObjRecv++;
         }
      }
   }

// Free work arrays

   delete [] cpu2obj;
   delete [] objsizetot;
   delete [] cpuidtot;

};

} // namespace k3d

