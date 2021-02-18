//------------------------------------------------------------------------------------------------
// File: k3d_base.cxx
//------------------------------------------------------------------------------------------------
#ifndef K3D_BASE_CXX
#define K3D_BASE_CXX

#include "k3d_base.hxx"

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

/// @brief Perform binary search
//========================================================================================
   int BinarySearch (long long _isup, int _nblks, long long *_blks, int _iblkprev)
   {
      int iblkprev = _iblkprev;
      if (iblkprev < 0)
           iblkprev = 0;
      if (iblkprev >= _nblks)
      {
         iblkprev = _nblks - 1;
      }
      if (_isup >= _blks[iblkprev] && _isup < _blks[iblkprev + 1])
      {
         return iblkprev;
      }
      int ibegblk = 0;
      int iendblk = _nblks;
      int iblk;
      while (true) {
         if (ibegblk == iendblk - 1) {
            iblkprev = ibegblk;
            break;
         } else {
            iblk = (ibegblk + iendblk - 1) / 2;
            if (iblk <= ibegblk)
               iblk = ibegblk + 1;
            if (_isup >= _blks[iblk] && _isup < _blks[iblk + 1]) {
               iblkprev = iblk;
               break;
            } else if (_isup < _blks[iblk]) {
               iendblk = iblk;
            } else {
               ibegblk = iblk + 1;
            }
         }
      }
      return iblkprev;
   }

/// @brief Print array
//========================================================================================
   void PrintArray (ostream & _stream, const char *_name, int _isize, const char *_charr)
   {
      _stream << _name << " (Size = " << _isize << ")" << endl;
      for (int i = 0; i < _isize; i++) {
         if ((i % 5 == 0) && (i > 0))
            _stream << endl;
         _stream << setw (12) << (int) _charr[i];
      }
      _stream << endl;
   }

/// @brief Print array
//========================================================================================
   void PrintArray (ostream & _stream, const char *_name, int _isize, const int *_iarr)
   {
      _stream << _name << " (Size = " << _isize << ")" << endl;
      for (int i = 0; i < _isize; i++) {
         if ((i % 5 == 0) && (i > 0))
            _stream << endl;
         _stream << setw (12) << _iarr[i];
      }
      _stream << endl;
   }

/// @brief Print array
//========================================================================================
   void PrintArray (ostream & _stream, const char *_name, int _isize,
                    const long long *_iarr)
   {
      _stream << _name << " (Size = " << _isize << ")" << endl;
      for (int i = 0; i < _isize; i++) {
         if ((i % 5 == 0) && (i > 0))
            _stream << endl;
         _stream << setw (18) << _iarr[i];
      }
      _stream << endl;
   }

/// @brief Print array
//========================================================================================
   void PrintArray (ostream & _stream, const char *_name, int _isize, const float *_farr)
   {
      _stream << _name << " (Size = " << _isize << ")" << endl;
      for (int i = 0; i < _isize; i++) {
         if ((i % 5 == 0) && (i > 0))
            _stream << endl;
         _stream << setw (14) << setprecision (8) << _farr[i] << " ";
      }
      _stream << endl;
   }

/// @brief Print array
//========================================================================================
   void PrintArray (ostream & _stream, const char *_name, int _isize, const double *_darr)
   {
      _stream << _name << " (Size = " << _isize << ")" << endl;
      for (int i = 0; i < _isize; i++) {
         if ((i % 3 == 0) && (i > 0))
            _stream << endl;
         _stream << setw (23) << setprecision (16) << _darr[i] << " ";
      }
      _stream << endl;
   }

/// @brief Print double array in low accuracy
//========================================================================================
   void PrintArrayLow (ostream & _stream, const char *_name, int _isize,
                       const double *_darr)
   {
      _stream << _name << " (Size = " << _isize << ")" << endl;
      for (int i = 0; i < _isize; i++) {
         if ((i % 5 == 0) && (i > 0))
            _stream << endl;
         _stream << setw (14) << setprecision (8) << _darr[i] << " ";
      }
      _stream << endl;
   }

/// @brief Set manip
//========================================================================================
   ostream & SetPw (ostream & stream) {
      stream.unsetf (ios::scientific);
      stream << ' ' << setprecision (6) << setw (10);
      return stream;
   }

/// @brief Round two values
//========================================================================================
   void Round (double &_dx, double &_dy)
   {
/*
      long long i, i0;

      double dx1 = 1.0e1 * _dx;

      i0 = (long long) dx1;
      i = (long long) ((dx1 - (double) i0) * 10000);
      dx1 = (double) i;
      dx1 = dx1 * 1.0e-4 + (double) i0;

      _dx = dx1 * 1.e-1;

      double dy1 = 1.0e1 * _dy;

      i0 = (long long) dy1;
      i = (long long) ((dy1 - (double) i0) * 10000);
      dy1 = (double) i;
      dy1 = dy1 * 1.0e-4 + (double) i0;

      _dy = dy1 * 1.e-1;
*/
      long long i, i0;

      double dx1 = 1.0e2 * _dx;

      i0 = (long long) dx1;
      i = (long long) ((dx1 - (double) i0) * 100000);
      dx1 = (double) i;
      dx1 = dx1 * 1.0e-5 + (double) i0;

      _dx = dx1 * 1.e-2;

      double dy1 = 1.0e2 * _dy;

      i0 = (long long) dy1;
      i = (long long) ((dy1 - (double) i0) * 100000);
      dy1 = (double) i;
      dy1 = dy1 * 1.0e-5 + (double) i0;

      _dy = dy1 * 1.e-2;
   }

/// @brief Read char array (direct access)
//========================================================================================
   void FGet (FILE * _file, size_t _size, char *_charr, size_t _offset)
   {

      size_t i;

      size_t j = sizeof (char) * _offset;

#ifdef _WINDOWS
      i = _fseeki64 (_file, j, SEEK_SET);
#else
      i = fseek (_file, j, SEEK_SET);
#endif
      if (i != 0) {
         cout << " Error when moving cursor for direct IO file " << endl;
         throw "Error when moving cursor for direct IO file";
      }

      i = fread (_charr, sizeof (char), _size, _file);
      if (i != _size) {
         cout << " Error when reading from direct IO file " << endl;
         throw "Error when reading from direct IO file";
      }

   }

/// @brief Determine the length of ascii file
//========================================================================================
   size_t LengthOfFile (const char *_filename)
   {

      size_t length = 0;
      FILE *BLN;
      if ((BLN = fopen (_filename, "rb")) != NULL) {
#ifdef _WINDOWS
         _fseeki64 (BLN, 0L, SEEK_END);
         length = (size_t) _ftelli64 (BLN);
#else
         fseek (BLN, 0L, SEEK_END);
         length = ftell (BLN);
#endif
         fclose (BLN);
      }
      return (length);

   }

/// @brief Read current file into char variable
//========================================================================================
   char *ReadFile (const char *_filename, size_t & _length, bool _bprint_out)
   {

      _length = 0;
      _length = LengthOfFile (_filename);

      if (_bprint_out) {
         cout << " Read from Ascii file {" << _filename << "}" << " Length = " << _length
            << endl;
      }

      char *filechar = 0;
      if (_length) {
         filechar = new char[_length + 1];
         if (filechar != NULL) {
            strcpy (filechar, "");
            FILE *file = fopen (_filename, "rb");
            if (ferror (file)) {
               char str[256] = "";
               perror (str);
               cout << " Error in open = " << str << endl;
               throw " ReadAsciiFile: error in open";
            }
            FGet (file, _length, filechar, 0);
            strcpy (filechar + _length, "");
            fclose (file);
         }
      }
      if (_bprint_out) {
         cout << " Reading completed successfully" << endl;
      }
      return (filechar);
   }

/// @brief Read current char from file
//========================================================================================
   char ReadChar (char *_file, size_t & _k, size_t _fend)
   {
      if (_k >= _fend)
         return ('\xFF');
      char buf = toupper (_file[_k++]);
      if (buf == ';') {
         while (buf != '\n' && _k < _fend)
            buf = _file[_k++];
         return ('\xFF');
      }
      if (isalpha (buf) || isdigit (buf) || buf == '+' || (buf == ':' && (_k == 1 || isalpha (_file[_k - 2])))  //-Werror=parentheses //suggest parentheses around ‘&&’ within ‘||’
          || buf == '-' || buf == '.' || buf == '#' || buf == '_' || buf == '\\'
          || buf == '!' || buf == '*' || buf == '%') {
         return (buf);
      } else
         return ('\xFF');
   }

/// @brief Read current word from file
//========================================================================================
   size_t ReadWord (char *_buf, char *_file, size_t & _k, size_t _fend)
   {

      unsigned int CHAR_BUF_SIZE = 256;

      size_t l = 0;
      while (l < CHAR_BUF_SIZE && _k < _fend) {
         char buf0 = ReadChar (_file, _k, _fend);
         if (buf0 != '\xFF') {
            _buf[l++] = buf0;
         } else {
            if (l != 0) {
               _buf[l] = '\x0';
               return (1);
            }
         }
      }
      _buf[l] = '\x0';
      return (l != 0);
   }

// Description: Read current word in braces from file
//========================================================================================
   size_t ReadWordInBraces (char *_buf, char *_file, size_t & _k, size_t _fend)
   {

      unsigned int CHAR_BUF_SIZE = 256;

      size_t fendloc = _fend;
      size_t i = ReadSymbol (_file, _k, fendloc, '"');
      if (i == 0)
         return i;
      _k++;
      size_t kloc = _k;
      fendloc = _fend;
      i = ReadSymbol (_file, _k, fendloc, '"');
      if (i == 0)
         return i;
      _k++;
      size_t l = 0;
      while (l < CHAR_BUF_SIZE && kloc < _k - 1) {
         _buf[l++] = _file[kloc++];
      }
      _buf[l] = '\x0';
      return (l != 0);

   }

// Description: Read prescribed symbol from file
//========================================================================================
   size_t ReadSymbol (char *_file, size_t & _k, size_t _fend, char _symbol)
   {

      char buf = _file[_k];

      while (buf != _symbol && buf != '\x0') {
         if (_k >= _fend - 1) {
            _k++;
            break;
         }
         buf = _file[++_k];
      }

      if (buf == '\x0')
         return (0);
      else
         return (1);

   }

/// @brief Replace one symbol by another
//========================================================================================
   void ReplaceSymbol (size_t _length, char *_arr, char _symb_ini, char _symb_fin)
   {
      int i;
      for (i = 0; i < (int) _length; i++) {
         if (_arr[i] == _symb_ini) {
            _arr[i] = _symb_fin;
         }
         if (_arr[i] == 0x0D) {
            _arr[i] = ' ';
         }
      }
   }

/// @brief Output hystogram for set of values
//========================================================================================
   void Hyst (ostream & _ffout, int _nhyst, double *_hyst, int _n_val, double *_values)
   {
      vector < int >num_arr (_nhyst + 1);
      int *pnum_arr = &num_arr[0];

      int i, j;

      for (i = 0; i < _nhyst; i++) {
         pnum_arr[i] = 0;
      }

      for (i = 0; i < _nhyst; i++) {
         for (j = 0; j < _n_val; j++) {
            if (_values[j] > _hyst[i])
               pnum_arr[i]++;
         }
      }

      _ffout << "  Hystogram of values: " << endl;

      for (i = 0; i < _nhyst; i++) {
         _ffout << "     Ihyst = " << i << " Thr = " << _hyst[i] << " Num = " <<
            pnum_arr[i] << endl;
      }

   }

// Compute optimal ordering for the matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::ComputeOptimalOrder (int _ordtype, const int _n,
                                   const vector < _Int > &_ia_alu,
                                   const vector < _Int > &_ja_alu, vector < int >&_order)
   {

// Compute symmetrized matrix

      vector < _Int > ia_symm;
      vector < _Int > ja_symm;

      CFct_impl < _Int, _Flt >::SymmetrizeSparsity (_n, _ia_alu, _ja_alu, ia_symm,
                                                    ja_symm);

// Call METIS ordering routine

#ifdef USE_METIS
      typedef idx_t TINT_TYPE;
#else
      typedef int TINT_TYPE;
#endif

      _order.resize (_n + 1);

      int nzjas = (int) ia_symm[_n];

      if (_ordtype == -1) {

         vector < TINT_TYPE > ialoc (_n + 1);
         vector < TINT_TYPE > jaloc (nzjas + 1);

         TINT_TYPE *pialoc = &ialoc[0];
         TINT_TYPE *pjaloc = &jaloc[0];

         int i, j, jj;

         ialoc[0] = 0;
         TINT_TYPE nz = 0;

         for (i = 0; i < _n; i++) {
            for (j = (int) ia_symm[i]; j < ia_symm[i + 1]; j++) {
               jj = (int) ja_symm[j];
               if (jj != i) {
                  pjaloc[nz] = jj;
                  nz++;
               }
            }
            pialoc[i + 1] = nz;
         }

         vector < TINT_TYPE > order_int (_n * 2 + 1);

         TINT_TYPE *porder_int = &order_int[0];

         TINT_TYPE nloc = _n;

         if (_n > 0) {

#ifdef USE_METIS

            idx_t options_ND[METIS_NOPTIONS];

            METIS_SetDefaultOptions (options_ND);

            options_ND[METIS_OPTION_NUMBERING] = 0;
            options_ND[METIS_OPTION_SEED] = 0;

            METIS_NodeND ((idx_t *) & nloc, (idx_t *) pialoc, (idx_t *) pjaloc, NULL,
                          (idx_t *) options_ND, (idx_t *) porder_int + nloc,
                          (idx_t *) porder_int);
#else
            cout <<
               " CFct_impl<>::ComputeOptimalOrder: error: Metis NodeND is not found !!! "
               << endl;
            throw
               " CFct_impl<>::ComputeOptimalOrder: error: Metis NodeND is not found !!! ";
#endif

         }

         for (i = 0; i < _n; i++) {
            _order[i] = (int) porder_int[i];
         }

      } else if (_ordtype == 0) {

// Identity ordering

         int i;

         for (i = 0; i < _n; i++) {
            _order[i] = i;
         }

      } else if (_ordtype == 1) {

// Reversed Cuthill-McCee ordering

         vector < int >ialoc (_n + 1);
         vector < int >jaloc (nzjas + 1);

         int *pialoc = &ialoc[0];
         int *pjaloc = &jaloc[0];

         int i;

         for (i = 0; i <= _n; i++)
            pialoc[i] = (int) (ia_symm[i] + 1);
         for (i = 0; i < nzjas; i++)
            pjaloc[i] = (int) (ja_symm[i] + 1);

         vector < int >xls (_n + 1);
         vector < int >mask (_n + 1);

         int *pxls = &xls[0];
         int *pmask = &mask[0];

         vector < int >order_int (_n + 1);

         int *porder_int = &order_int[0];

         void genrcm (int n, int *xadj, int *iadj, int *perm, int *xls, int *mask);

         genrcm (_n, pialoc, pjaloc, porder_int, pxls, pmask);

         for (i = 0; i < _n; i++) {
            porder_int[i] -= 1;
         }

         for (i = 0; i < _n; i++) {
            pmask[porder_int[i]] = i;
         }

         for (i = 0; i < _n; i++) {
            _order[i] = pmask[i];
         }

      }

   }

/// @brief Get main submatrix sparsity
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::GetSubmatrixSp (int _n, const vector < _Int > &_ia,
                              const vector < _Int > &_ja, int _n1,
                              vector < _Int > &_ia_asub, vector < _Int > &_ja_asub)
   {

      int nzja_loc = (int) _ia[_n1];

      _ia_asub.resize (_n1 + 1);
      _ja_asub.resize (nzja_loc + 1);

      _Int *p_ia_asub = &_ia_asub[0];
      _Int *p_ja_asub = &_ja_asub[0];

      int i, j, jj;

      p_ia_asub[0] = 0;

      nzja_loc = 0;

      for (i = 0; i < _n1; i++) {
         for (j = (int) _ia[i]; j < _ia[i + 1]; j++) {
            jj = (int) _ja[j];
            if (jj < _n1) {
               p_ja_asub[nzja_loc] = jj;
               nzja_loc++;
            }
         }
         p_ia_asub[i + 1] = nzja_loc;
      }

   }

/// @brief Get main submatrix sparsity
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::GetSubmatrixSp (int _n, const vector < _Int > &_ia,
                              const vector < _Int > &_ja, const vector < char >&_jachar,
                              int _n1, vector < _Int > &_ia_asub,
                              vector < _Int > &_ja_asub, vector < char >&_jachar_asub)
   {

      int nzja_loc = (int) _ia[_n1];

      _ia_asub.resize (_n1 + 1);
      _ja_asub.resize (nzja_loc + 1);
      _jachar_asub.resize (nzja_loc + 1);

      _Int *p_ia_asub = &_ia_asub[0];
      _Int *p_ja_asub = &_ja_asub[0];
      char *p_jachar_asub = &_jachar_asub[0];

      int i, j, jj;

      p_ia_asub[0] = 0;

      nzja_loc = 0;

      for (i = 0; i < _n1; i++) {
         for (j = (int) _ia[i]; j < _ia[i + 1]; j++) {
            jj = (int) _ja[j];
            if (jj < _n1) {
               p_ja_asub[nzja_loc] = jj;
               p_jachar_asub[nzja_loc] = _jachar[j];
               nzja_loc++;
            }
         }
         p_ia_asub[i + 1] = nzja_loc;
      }

   }

/// @brief Get last submatrix sparsity
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::GetLastSubmatrixSp (int _n, const vector < _Int > &_ia,
                                  const vector < _Int > &_ja, int _n1,
                                  vector < _Int > &_ia_asub_last,
                                  vector < _Int > &_ja_asub_last)
   {

      int n2 = _n - _n1;

      int nzja_loc = (int) (_ia[_n] - _ia[_n1]);

      _ia_asub_last.resize (n2 + 1);
      _ja_asub_last.resize (nzja_loc + 1);

      _Int *p_ia_asub_last = &_ia_asub_last[0];
      _Int *p_ja_asub_last = &_ja_asub_last[0];

      int i, j, jj;

      p_ia_asub_last[0] = 0;

      nzja_loc = 0;

      for (i = 0; i < n2; i++) {
         for (j = (int) _ia[_n1 + i]; j < _ia[_n1 + i + 1]; j++) {
            jj = (int) _ja[j];
            if (jj >= _n1) {
               p_ja_asub_last[nzja_loc] = jj - _n1;
               nzja_loc++;
            }
         }
         p_ia_asub_last[i + 1] = nzja_loc;
      }

   }

/// @brief Get last submatrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::GetLastSubmatrixBxB (bool _b_is_char, int _n, int _size_bxb,
                                   const vector < _Int > &_ia, const vector < _Int > &_ja,
                                   const vector < char >&_ja_char,
                                   const vector < _Flt > &_a, int _n1,
                                   vector < _Int > &_ia_asub_last,
                                   vector < _Int > &_ja_asub_last,
                                   vector < char >&_jachar_asub_last,
                                   vector < _Flt > &_a_asub_last)
   {

      int n2 = _n - _n1;

      const _Flt *p_a = _a.data ();

      int nzja_loc = (int) (_ia[_n] - _ia[_n1]);

      _ia_asub_last.resize (n2 + 1);
      _ja_asub_last.resize (nzja_loc + 1);
      if (_b_is_char) {
         _jachar_asub_last.resize (nzja_loc + 1);
      }
      _a_asub_last.resize (nzja_loc * _size_bxb + 1);

      _Int *p_ia_asub_last = _ia_asub_last.data ();
      _Int *p_ja_asub_last = _ja_asub_last.data ();
      char *p_jachar_asub_last = _jachar_asub_last.data ();
      _Flt *p_a_asub_last = _a_asub_last.data ();

      int i, j, jj;

      p_ia_asub_last[0] = 0;

      nzja_loc = 0;

      for (i = 0; i < n2; i++) {
         for (j = (int) _ia[_n1 + i]; j < _ia[_n1 + i + 1]; j++) {
            jj = (int) _ja[j];
            if (jj >= _n1) {
               p_ja_asub_last[nzja_loc] = jj - _n1;
               if (_b_is_char)
                  p_jachar_asub_last[nzja_loc] = _ja_char[j];
               CVector < _Flt >::CopyVector (_size_bxb, p_a + j * _size_bxb,
                                             p_a_asub_last + nzja_loc * _size_bxb);
               nzja_loc++;
            }
         }
         p_ia_asub_last[i + 1] = nzja_loc;
      }

   }

/// @brief Find all separators for ND ordered matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::FindAllSeparators (const int _n, const vector < _Int > &_ia_alu,
                                 const vector < _Int > &_ja_alu, int &_nblks,
                                 vector < int >&_blks)
   {

      int nloc = _n;

      const _Int *pia = &_ia_alu[0];
      const _Int *pja = &_ja_alu[0];

// Allocate work arrays

      vector < int >profile (nloc + 1);
      vector < int >profile_monotone (nloc + 1);
      vector < int >list_sep (nloc + 1);
      vector < int >list_bnd (nloc + 1);
      vector < int >list_ord (2 * nloc + 3);

      int *pprofile = &profile[0];
      int *pprofile_monotone = &profile_monotone[0];
      int *plist_sep = &list_sep[0];
      int *plist_bnd = &list_bnd[0];
      int *plist_ord = &list_ord[0];

// Compute initial profile

      int i, ibeg, jj;

      for (i = 0; i < nloc; i++) {
         pprofile[i] = i;
         ibeg = (int) pia[i];
         if (ibeg < pia[i + 1]) {
            jj = (int) pja[ibeg];
            if (jj < i)
               pprofile[i] = jj;
         }
         pprofile_monotone[i] = pprofile[i];
      }

// Compute initial separators columns (up diagonal part is empty)

      int nlist_sep = 0;

      for (i = 0; i < nloc; i++) {
         if (pprofile[i] == i) {
            plist_sep[nlist_sep] = i;
            nlist_sep++;
         }
      }

// For each separator find its secondary boundary

      int i_sep, i_sep_next, ival_curr, jval, j;

      for (i = nlist_sep - 1; i >= 0; i--) {

         i_sep = plist_sep[i];

// Initial monotonization

         i_sep_next = nloc;
         if (i < nlist_sep - 1)
            i_sep_next = plist_sep[i + 1];

         ival_curr = i_sep;

         for (j = i_sep; j < nloc; j++) {
            jval = pprofile_monotone[j];
            if (j < i_sep_next) {
               if (jval < ival_curr) {
                  ival_curr = jval;
               }
            } else {
               if (jval < ival_curr) {
                  break;
               }
            }
            pprofile_monotone[j] = ival_curr;
         }

// Find secondary boundary

         plist_bnd[i] = nloc;

         for (j = i_sep; j < nloc; j++) {
            if (pprofile_monotone[j] < i_sep) {
               plist_bnd[i] = j;
               break;
            }
         }

      }

// Create the final blocks partitioning

      for (i = 0; i < nlist_sep; i++) {
         plist_ord[i * 2] = plist_sep[i];
         plist_ord[i * 2 + 1] = plist_bnd[i];
      }
      plist_ord[2 * nlist_sep] = 0;
      plist_ord[2 * nlist_sep + 1] = nloc;

      sort (plist_ord, plist_ord + 2 * nlist_sep + 2);

      _nblks = 0;

      int ind_curr = -1;

      for (i = 0; i < 2 * nlist_sep + 2; i++) {
         j = plist_ord[i];
         if (j != ind_curr) {
            ind_curr = j;
            plist_ord[_nblks] = ind_curr;
            _nblks++;
         }
      }

      _nblks--;

      _blks.resize (_nblks + 2);

      int *p_blks = &_blks[0];

      for (i = 0; i <= _nblks; i++)
         p_blks[i] = plist_ord[i];

   }

/// @brief Find set of separators for given tree
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::FindSeparatorsForTree (const int _n, const vector < _Int > &_ia_alu,
                                     const vector < _Int > &_ja_alu, CTree & _tree,
                                     int &_nblks, vector < int >&_blks)
   {

      int nloc = _n;

      const _Int *pia = &_ia_alu[0];
      const _Int *pja = &_ja_alu[0];

// Find all initial separators

      int nblksloc_sep;
      vector < int >blks_sep;

      CFct_impl < _Int, _Flt >::FindAllSeparators (_n, _ia_alu, _ja_alu, nblksloc_sep,
                                                   blks_sep);

      int *pblks_sep = &blks_sep[0];

// Find separators and their boundaries

      vector < int >list_sep (nblksloc_sep + 1);
      vector < int >list_bnd (nblksloc_sep + 1);

      int *plist_sep = &list_sep[0];
      int *plist_bnd = &list_bnd[0];

      int n_sep = 0;

      int i, j, i_bnd, i_sep, ibeg, jj, j_sep, jbeg, kk;

      for (i = 1; i < nblksloc_sep; i++) {
         i_sep = pblks_sep[i];
         ibeg = (int) pia[i_sep];
         if (ibeg < pia[i_sep + 1]) {
            jj = (int) pja[ibeg];
            if (jj >= i_sep) {
               i_bnd = nloc;
               for (j = i + 1; j < nblksloc_sep; j++) {
                  j_sep = pblks_sep[j];
                  jbeg = (int) pia[j_sep];
                  if (jbeg < pia[j_sep + 1]) {
                     kk = (int) pja[jbeg];
                     if (kk < i_sep) {
                        i_bnd = j_sep;
                        break;
                     }
                  }
               }
               plist_sep[n_sep] = i_sep;
               plist_bnd[n_sep] = i_bnd;
               n_sep++;
            }
         }
      }

// Create initial partitioning

      int nnodes_loc = _tree.GetNnodes ();
      int root_loc = _tree.GetRootId ();
      int nlev_loc = _tree.GetNlev ();

      vector < int >list_nd_1 (nnodes_loc + 1);
      vector < int >list_nd_2 (nnodes_loc + 1);
      vector < int >ibegblk (nnodes_loc + 1);
      vector < int >iendblk (nnodes_loc + 1);
      vector < int >ibegblk_sep (nnodes_loc + 1);
      vector < int >iendblk_sep (nnodes_loc + 1);

      int *plist_nd_1 = &list_nd_1[0];
      int *plist_nd_2 = &list_nd_2[0];
      int *pibegblk = &ibegblk[0];
      int *piendblk = &iendblk[0];
      int *pibegblk_sep = &ibegblk_sep[0];
      int *piendblk_sep = &iendblk_sep[0];

      int *plist_nd_curr = plist_nd_1;
      int *plist_nd_next = plist_nd_2;

      int nlist_nd_curr = 1;
      plist_nd_curr[0] = root_loc;

      for (i = 0; i < nnodes_loc; i++)
         pibegblk[i] = -1;
      for (i = 0; i < nnodes_loc; i++)
         piendblk[i] = -1;

      pibegblk[root_loc] = 0;
      piendblk[root_loc] = nloc - 1;

      pibegblk_sep[root_loc] = 0;
      piendblk_sep[root_loc] = n_sep - 1;

// Cycle over the tree levels

      int *pnchilds = _tree.GetNchilds ();
      vector < int >*pchilds_list = _tree.GetChildsList ();

      int nlist_nd_next, ilist, inode_curr, nchilds_curr, ichild1, ichild2, ibegloc,
         iendloc, ind1, ind2;

      int *pchilds_curr;

      int ilevloc = 1;

      int ibeg_sep, iend_sep, ibeg_sep_loc, iend_sep_loc, isep, ind_opt, jj_bnd, niloc,
         njloc;
      double square, square_opt;

      while (ilevloc < nlev_loc) {

         ilevloc++;

// Cycle over the nodes on current level

         nlist_nd_next = 0;

         for (ilist = 0; ilist < nlist_nd_curr; ilist++) {

            inode_curr = plist_nd_curr[ilist];

            nchilds_curr = pnchilds[inode_curr];
            pchilds_curr = &pchilds_list[inode_curr][0];

            if (nchilds_curr != 2)
               throw
                  " CTBlockSp<_TInt>::FindSeparatorsForBinaryTree: incorrect binary tree ";

            ichild1 = pchilds_curr[0];
            ichild2 = pchilds_curr[1];

            ibegloc = pibegblk[inode_curr];
            iendloc = piendblk[inode_curr];

            ibeg_sep = pibegblk_sep[inode_curr];
            iend_sep = piendblk_sep[inode_curr];

// Find optimal separator

            if (iendloc >= ibegloc) {
               if (ibeg_sep <= iend_sep) {
                  ibeg_sep_loc = -1;
                  iend_sep_loc = -1;
                  for (i = ibeg_sep; i <= iend_sep; i++) {
                     jj = plist_sep[i];
                     if (jj >= ibegloc && jj <= iendloc) {
                        if (ibeg_sep_loc == -1)
                           ibeg_sep_loc = i;
                        iend_sep_loc = i;
                     }
                  }
                  if (ibeg_sep_loc == -1) {
                     ind1 = -1;
                     ind2 = -1;
                     isep = -1;
                  } else {
                     square_opt = -1;
                     ind_opt = -1;
                     for (i = ibeg_sep_loc; i <= iend_sep_loc; i++) {
                        jj = plist_sep[i];
                        jj_bnd = plist_bnd[i];
                        if (jj_bnd > iendloc)
                           jj_bnd = iendloc + 1;
                        niloc = jj - ibegloc;
                        njloc = jj_bnd - jj;
                        square = ((double) niloc) * ((double) njloc);
                        if (square > square_opt) {
                           square_opt = square;
                           ind_opt = i;
                        }
                     }
                     ind1 = plist_sep[ind_opt];
                     ind2 = plist_bnd[ind_opt];
                     if (ind2 > iendloc)
                        ind2 = iendloc + 1;
                     isep = ind_opt;
                  }
               } else {
                  ind1 = -1;
                  ind2 = -1;
                  isep = -1;
               }
            } else {
               ind1 = -1;
               ind2 = -1;
               isep = -1;
            }

            if (ind1 < 0 || ind2 < 0) {
               pibegblk[ichild1] = ibegloc;
               piendblk[ichild1] = ibegloc - 1;
               pibegblk[ichild2] = ibegloc;
               piendblk[ichild2] = ibegloc - 1;
               pibegblk_sep[ichild1] = ibeg_sep;
               piendblk_sep[ichild1] = ibeg_sep - 1;
               pibegblk_sep[ichild2] = ibeg_sep;
               piendblk_sep[ichild2] = ibeg_sep - 1;
            } else {
               pibegblk[ichild1] = ibegloc;
               piendblk[ichild1] = ind1 - 1;
               pibegblk[ichild2] = ind1;
               piendblk[ichild2] = ind2 - 1;

               pibegblk_sep[ichild1] = ibeg_sep;
               piendblk_sep[ichild1] = isep - 1;
               pibegblk_sep[ichild2] = isep + 1;
               piendblk_sep[ichild2] = iend_sep;

               pibegblk[inode_curr] = ind2;
            }

            plist_nd_next[nlist_nd_next] = ichild1;
            plist_nd_next[nlist_nd_next + 1] = ichild2;
            nlist_nd_next += 2;

         }

// Switch data

         if (plist_nd_curr == plist_nd_1) {
            plist_nd_curr = plist_nd_2;
            plist_nd_next = plist_nd_1;
         } else {
            plist_nd_curr = plist_nd_1;
            plist_nd_next = plist_nd_2;
         }

         nlist_nd_curr = nlist_nd_next;

      }

// Prepare return data

      _nblks = nnodes_loc;

      _blks.resize (_nblks + 2);

      int *p_blks = &_blks[0];

      for (i = 0; i < nnodes_loc; i++)
         p_blks[i] = pibegblk[i];

      p_blks[nnodes_loc] = nloc;

      sort (p_blks, p_blks + nnodes_loc + 1);

   }

//
// Compute optimal ordering for the matrix via splitting into 3 blocks: 1-st order Schur complement for third block, which is the separator between first block and remaining data
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::ComputeOptimalOrderSchur (int _ordtype, const int _n, int _n1,
                                        const vector < _Int > &_ia_alu,
                                        const vector < _Int > &_ja_alu, int &_n2,
                                        vector < int >&_order)
   {

      _order.resize (_n + 1);

      int *porder = &_order[0];

// Compute symmetrized matrix

      vector < _Int > ia_symm;
      vector < _Int > ja_symm;

      CFct_impl < _Int, _Flt >::SymmetrizeSparsity (_n, _ia_alu, _ja_alu, ia_symm,
                                                    ja_symm);

      _Int *pia_symm = &ia_symm[0];
      _Int *pja_symm = &ja_symm[0];

// Get submatrix corresponding to the first block and compute optimal ordering for it

      int nzja_1 = 0;

      int i, j, jj;

      for (i = 0; i < _n1; i++) {
         for (j = (int) pia_symm[i]; j < pia_symm[i + 1]; j++) {
            jj = (int) pja_symm[j];
            if (jj < _n1)
               nzja_1++;
         }
      }

      vector < _Int > ia1 (_n1 + 1);
      vector < _Int > ja1 (nzja_1 + 1);

      _Int *pia1 = &ia1[0];
      _Int *pja1 = &ja1[0];

      nzja_1 = 0;
      pia1[0] = 0;

      for (i = 0; i < _n1; i++) {
         for (j = (int) pia_symm[i]; j < pia_symm[i + 1]; j++) {
            jj = (int) pja_symm[j];
            if (jj < _n1) {
               pja1[nzja_1] = jj;
               nzja_1++;
            }
         }
         pia1[i + 1] = nzja_1;
      }

      vector < int >order1 (1);

      CFct_impl < _Int, _Flt >::ComputeOptimalOrder (_ordtype, _n1, ia1, ja1, order1);

      int *porder1 = &order1[0];

      for (i = 0; i < _n1; i++)
         porder[i] = porder1[i];

      {
         vector < _Int > ia1_temp (1);
         vector < _Int > ja1_temp (1);
         vector < int >order1_temp (1);
         ia1.swap (ia1_temp);
         ja1.swap (ja1_temp);
         order1.swap (order1);
      }

// Register columns in the second part which is to be moved to the bordering

      vector < int >imask (_n + 1);
      vector < int >list (_n + 1);
      vector < int >indarr (_n + 1);

      int *pimask = &imask[0];
      int *plist = &list[0];
      int *pindarr = &indarr[0];

      for (i = _n1; i < _n; i++)
         pimask[i] = -1;

      int nlist_bord = 0;

      for (i = 0; i < _n1; i++) {
         for (j = (int) pia_symm[i]; j < pia_symm[i + 1]; j++) {
            jj = (int) pja_symm[j];
            if (jj >= _n1 && pimask[jj] == -1) {
               plist[nlist_bord] = jj;
               nlist_bord++;
               pimask[jj] = 1;
            }
         }
      }

      sort (plist, plist + nlist_bord);

      _n2 = _n - nlist_bord;

      int ip2 = _n1;
      int ip3 = _n2;

      for (i = _n1; i < _n; i++) {
         if (pimask[i] == -1) {
            porder[i] = ip2;
            ip2++;
         } else {
            porder[i] = ip3;
            ip3++;
         }
      }

      vector < int >iorder (_n + 1);

      int *piorder = &iorder[0];

      for (i = 0; i < _n; i++)
         piorder[porder[i]] = i;

// Reorder the matrix

      vector < _Int > ia_ord (1);
      vector < _Int > ja_ord (1);

      CFct_impl < _Int, _Flt >::ReorderMatrixSp (_n, _order, ia_symm, ja_symm, ia_ord,
                                                 ja_ord);

      _Int *pia_ord = &ia_ord[0];
      _Int *pja_ord = &ja_ord[0];

// Get submatrix 2 and compute its ordering

      int nzja_2 = 0;

      for (i = _n1; i < _n2; i++) {
         for (j = (int) pia_ord[i]; j < pia_ord[i + 1]; j++) {
            jj = (int) pja_ord[j];
            if (jj >= _n1 && jj < _n2)
               nzja_2++;
         }
      }

      int n12 = _n2 - _n1;

      vector < _Int > ia2 (n12 + 1);
      vector < _Int > ja2 (nzja_2 + 1);

      _Int *pia2 = &ia2[0];
      _Int *pja2 = &ja2[0];

      nzja_2 = 0;
      pia2[0] = 0;

      for (i = _n1; i < _n2; i++) {
         for (j = (int) pia_ord[i]; j < pia_ord[i + 1]; j++) {
            jj = (int) pja_ord[j];
            if (jj >= _n1 && jj < _n2) {
               pja2[nzja_2] = jj - _n1;
               nzja_2++;
            }
         }
         pia2[i - _n1 + 1] = nzja_2;
      }

      vector < int >order2 (1);

      CFct_impl < _Int, _Flt >::ComputeOptimalOrder (_ordtype, n12, ia2, ja2, order2);

      int *porder2 = &order2[0];

      int iold;

      for (i = 0; i < n12; i++) {
         iold = piorder[i + _n1];
         porder[iold] = porder2[i] + _n1;
      }

      {
         vector < _Int > ia2_temp (1);
         vector < _Int > ja2_temp (1);
         vector < int >order2_temp (1);
         ia2.swap (ia2_temp);
         ja2.swap (ja2_temp);
         order2.swap (order2);
      }

// Compute ends of second block in all rows

      vector < int >ia_end2 (_n + 1);

      int *pia_end2 = &ia_end2[0];

      for (i = 0; i < _n; i++) {
         pia_end2[i] = (int) pia_ord[i];
         for (j = (int) pia_ord[i]; j < pia_ord[i + 1]; j++) {
            jj = (int) pja_ord[j];
            if (jj < _n2)
               pia_end2[i] = j;
         }
      }

// Compute extended bordering data

      int n23 = _n - _n2;

      for (i = 0; i < _n; i++)
         pimask[i] = -1;

      int icycle = -1;

      vector < int >ia_bord (n23 + 1);

      int *pia_bord = &ia_bord[0];

      pia_bord[0] = 0;
      int nzja_bord = 0;

      vector < int >ja_bord (1);

      int nlistloc, k, kk;

      for (i = _n2; i < _n; i++) {
         icycle++;
         nlistloc = 0;
         for (j = (int) pia_ord[i]; j <= pia_end2[i]; j++) {
            jj = (int) pja_ord[j];
            for (k = (int) pia_ord[jj]; k <= pia_end2[jj]; k++) {
               kk = (int) pja_ord[k];
               if (pimask[kk] != icycle) {
                  plist[nlistloc] = kk;
                  nlistloc++;
                  pimask[kk] = icycle;
               }
            }
         }
         sort (plist, plist + nlistloc);
         ja_bord.resize (nzja_bord + nlistloc + 1);
         for (j = 0; j < nlistloc; j++)
            ja_bord[nzja_bord + j] = plist[j];
         nzja_bord += nlistloc;
         pia_bord[i - _n2 + 1] = nzja_bord;
      }

      int *pja_bord = &ja_bord[0];

// Condense, renumber and extend

      icycle++;
      nlistloc = 0;
      for (j = 0; j < nzja_bord; j++) {
         jj = pja_bord[j];
         if (pimask[jj] != icycle) {
            plist[nlistloc] = jj;
            nlistloc++;
            pimask[jj] = icycle;
         }
      }
      sort (plist, plist + nlistloc);
      for (j = 0; j < nlistloc; j++) {
         jj = plist[j];
         pindarr[jj] = j;
      }
      for (j = 0; j < nzja_bord; j++) {
         jj = pja_bord[j];
         pja_bord[j] = pindarr[jj];
      }

      ia_bord.resize (nlistloc + 1);

      pia_bord = &ia_bord[0];

      for (j = n23 + 1; j <= nlistloc; j++)
         pia_bord[j] = nzja_bord;

// Transpose

      vector < int >ia_bord_t (1);
      vector < int >ja_bord_t (1);

      CFct_impl < int, _Flt >::TransposeSp (nlistloc, ia_bord, ja_bord, ia_bord_t,
                                            ja_bord_t);

      int *pia_bord_t = &ia_bord_t[0];
      int *pja_bord_t = &ja_bord_t[0];

// Compute Schur complement via extended block data

      vector < _Int > ia3_Schur (n23 + 1);
      vector < _Int > ja3_Schur (1);

      _Int *pia3_Schur = &ia3_Schur[0];

      int nzja_Schur = 0;

      for (i = 0; i < n23; i++) {
         icycle++;
         nlistloc = 0;
         for (j = pia_bord[i]; j < pia_bord[i + 1]; j++) {
            jj = pja_bord[j];
            for (k = pia_bord_t[jj]; k < pia_bord_t[jj + 1]; k++) {
               kk = pja_bord_t[k];
               if (pimask[kk] != icycle) {
                  plist[nlistloc] = kk;
                  nlistloc++;
                  pimask[kk] = icycle;
               }
            }
         }
         sort (plist, plist + nlistloc);
         ja3_Schur.resize (nzja_Schur + nlistloc + 1);
         for (j = 0; j < nlistloc; j++)
            ja3_Schur[nzja_Schur + j] = plist[j];
         nzja_Schur += nlistloc;
         pia3_Schur[i + 1] = nzja_Schur;
      }

      _Int *pja3_Schur = &ja3_Schur[0];

// Add submatrices

      vector < _Int > ia3 (n23 + 1);
      vector < _Int > ja3 (1);

      int nzja_3 = 0;

      _Int *pia3 = &ia3[0];

      for (i = 0; i < n23; i++) {
         icycle++;
         nlistloc = 0;
         for (j = (int) pia3_Schur[i]; j < pia3_Schur[i + 1]; j++) {
            jj = (int) pja3_Schur[j];
            plist[nlistloc] = jj;
            nlistloc++;
            pimask[jj] = icycle;
         }
         for (j = pia_end2[i + _n2] + 1; j < pia_ord[i + _n2 + 1]; j++) {
            jj = (int) pja_ord[j] - _n2;
            if (pimask[jj] != icycle) {
               plist[nlistloc] = jj;
               nlistloc++;
               pimask[jj] = icycle;
            }
         }
         sort (plist, plist + nlistloc);
         ja3.resize (nzja_3 + nlistloc + 1);
         for (j = 0; j < nlistloc; j++)
            ja3[nzja_3 + j] = plist[j];
         nzja_3 += nlistloc;
         pia3[i + 1] = nzja_3;
      }

// Optimal order

      vector < int >order3 (1);

      CFct_impl < _Int, _Flt >::ComputeOptimalOrder (_ordtype, n23, ia3, ja3, order3);

      int *porder3 = &order3[0];

      for (i = 0; i < n23; i++) {
         iold = piorder[i + _n2];
         porder[iold] = porder3[i] + _n2;
      }

   }

//
// Compute ordered matrix (sparsity only)
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::ReorderMatrixSp (int _n, const vector < int >&_order,
                               const vector < _Int > &_ia_alu,
                               const vector < _Int > &_ja_alu,
                               vector < _Int > &_ia_alu_ord, vector < _Int > &_ja_alu_ord)
   {

// Compute inverse order

      vector < int >iord (_n + 1);

      int *piord = &iord[0];

      int i;

      for (i = 0; i < _n; i++)
         piord[_order[i]] = i;

// Reorder the matrix

      int nzja = (int) _ia_alu[_n];

      _ia_alu_ord.resize (_n + 1);
      _ja_alu_ord.resize (nzja + 1);

      _Int *piaord = &_ia_alu_ord[0];
      _Int *pjaord = &_ja_alu_ord[0];

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

      CSortInt *piiarr = &iiarr[0];

      int j;

      for (i = 0; i < _n; i++) {
         j = _order[i];
         piaord[j + 1] = _ia_alu[i + 1] - _ia_alu[i];
      }

      for (i = 0; i < _n; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, nzloc, ibs, jold, jnew;

      for (inew = 0; inew < _n; inew++) {
         i = piord[inew];
         nzloc = (int) (_ia_alu[i + 1] - _ia_alu[i]);
         ibs = nz;
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jold = (int) _ja_alu[j];
            jnew = _order[jold];
            pjaord[nz] = (_Int) jnew;
            nz++;
         }

// Sort elements

         for (j = 0; j < nzloc; j++) {
            piiarr[j].ival = (int) pjaord[ibs + j];
            piiarr[j].i2val = j;
         }

         sort (piiarr, piiarr + nzloc);

         for (j = 0; j < nzloc; j++) {
            pjaord[ibs + j] = (_Int) piiarr[j].ival;
         }

      }

   }

//
// Compute ordered matrix (sparsity only)
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::ReorderMatrixSp (int _n, const vector < int >&_order,
                               const vector < _Int > &_ia_alu,
                               const vector < _Int > &_ja_alu,
                               const vector < char >&_jachar_alu,
                               vector < _Int > &_ia_alu_ord, vector < _Int > &_ja_alu_ord,
                               vector < char >&_jachar_alu_ord)
   {

// Compute inverse order

      vector < int >iord (_n + 1);

      int *piord = &iord[0];

      int i;

      for (i = 0; i < _n; i++)
         piord[_order[i]] = i;

// Reorder the matrix

      int nzja = (int) _ia_alu[_n];

      _ia_alu_ord.resize (_n + 1);
      _ja_alu_ord.resize (nzja + 1);
      _jachar_alu_ord.resize (nzja + 1);

      _Int *piaord = &_ia_alu_ord[0];
      _Int *pjaord = &_ja_alu_ord[0];
      char *pjacharord = &_jachar_alu_ord[0];

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
      vector < char >chars (nimax + 1);

      CSortInt *piiarr = &iiarr[0];
      char *pchars = &chars[0];

      int j;

      for (i = 0; i < _n; i++) {
         j = _order[i];
         piaord[j + 1] = _ia_alu[i + 1] - _ia_alu[i];
      }

      for (i = 0; i < _n; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, nzloc, ibs, jold, jnew, ind;

      for (inew = 0; inew < _n; inew++) {
         i = piord[inew];
         nzloc = (int) (_ia_alu[i + 1] - _ia_alu[i]);
         ibs = nz;
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jold = (int) _ja_alu[j];
            jnew = _order[jold];
            pjaord[nz] = (_Int) jnew;
            pjacharord[nz] = _jachar_alu[j];
            nz++;
         }

// Sort elements

         for (j = 0; j < nzloc; j++) {
            piiarr[j].ival = (int) pjaord[ibs + j];
            piiarr[j].i2val = j;
            pchars[j] = pjacharord[ibs + j];
         }

         sort (piiarr, piiarr + nzloc);

         for (j = 0; j < nzloc; j++) {
            ind = piiarr[j].i2val;
            pjaord[ibs + j] = (_Int) piiarr[j].ival;
            pjacharord[ibs + j] = pchars[ind];
         }

      }

   }

//
// Compute ordered matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::ReorderMatrix (int _n, const vector < int >&_order,
                             vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                             vector < _Flt > &_a_alu, vector < _Int > &_ia_alu_ord,
                             vector < _Int > &_ja_alu_ord, vector < _Flt > &_a_alu_ord)
   {

// Compute inverse order

      vector < int >iord (_n + 1);

      int *piord = &iord[0];

      int i;

      for (i = 0; i < _n; i++)
         piord[_order[i]] = i;

// Reorder the matrix

      int nzja = (int) _ia_alu[_n];

      _ia_alu_ord.resize (_n + 1);
      _ja_alu_ord.resize (nzja + 1);
      _a_alu_ord.resize (nzja + 1);

      _Int *piaord = &_ia_alu_ord[0];
      _Int *pjaord = &_ja_alu_ord[0];
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
      vector < _Flt > elems (nimax + 1);

      CSortInt *piiarr = &iiarr[0];
      _Flt *pelems = &elems[0];

      int j;

      for (i = 0; i < _n; i++) {
         j = _order[i];
         piaord[j + 1] = _ia_alu[i + 1] - _ia_alu[i];
      }

      for (i = 0; i < _n; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, nzloc, ibs, jold, jnew;

      for (inew = 0; inew < _n; inew++) {
         i = piord[inew];
         nzloc = (int) (_ia_alu[i + 1] - _ia_alu[i]);
         ibs = nz;
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jold = (int) _ja_alu[j];
            jnew = _order[jold];
            pjaord[nz] = (_Int) jnew;
            paord[nz] = _a_alu[j];
            nz++;
         }

// Sort elements

         for (j = 0; j < nzloc; j++) {
            piiarr[j].ival = (int) pjaord[ibs + j];
            piiarr[j].i2val = j;
         }

         sort (piiarr, piiarr + nzloc);

         for (j = 0; j < nzloc; j++)
            pelems[j] = paord[ibs + j];

         for (j = 0; j < nzloc; j++) {
            pjaord[ibs + j] = (_Int) piiarr[j].ival;
            paord[ibs + j] = pelems[piiarr[j].i2val];
         }

      }

   }

//
// Perform ILU2 point factorization of the block with future diagonal modification
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::Ilu2BlockTransform (int _sctype, int _nitersc, int _fcttype, double _pivmin,
                                  double _tau1, double _tau2, double _theta, int _n,
                                  vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                                  vector < _Flt > &_a_alu, vector < _Int > &_ia_lu,
                                  vector < _Int > &_ida_lu, vector < _Int > &_ja_lu,
                                  vector < _Flt > &_a_lu, double &_eigmin_att,
                                  double &_eigmax_att)
   {

// Perform fct

      vector < _Int > ia_lt;
      vector < _Int > ja_lt;
      vector < _Flt > a_lt;
      vector < _Int > ia_u;
      vector < _Int > ja_u;
      vector < _Flt > a_u;

      double sclmin_att, sclmax_att;

      int nmodif;

      CFct_impl < _Int, _Flt >::Ilu2Block (_sctype, _nitersc, _fcttype, _pivmin, _tau1,
                                           _tau2, _theta, _n, _ia_alu, _ja_alu, _a_alu,
                                           ia_lt, ja_lt, a_lt, ia_u, ja_u, a_u,
                                           sclmin_att, sclmax_att, nmodif, _eigmin_att,
                                           _eigmax_att);

// Transpose L

      vector < char >char_dummy;

      vector < _Int > ia_l;
      vector < _Int > ja_l;
      vector < _Flt > a_l;

      CFct_impl < _Int, _Flt >::Transpose (false, _n, ia_lt, ja_lt, char_dummy, a_lt,
                                           ia_l, ja_l, char_dummy, a_l);

      {
         vector < _Int > ia_dummy;
         vector < _Int > ja_dummy;
         vector < _Flt > a_dummy;
         ia_lt.swap (ia_dummy);
         ja_lt.swap (ja_dummy);
         a_lt.swap (a_dummy);
      }

// Combine rows LU

      CFct_impl < _Int, _Flt >::CombineRowsLU (_n, ia_l, ja_l, a_l, ia_u, ja_u, a_u,
                                               _ia_lu, _ida_lu, _ja_lu, _a_lu);

   }

//
// Perform ICH2 point factorization of the block with future diagonal modification
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::Ich2Block (double _sclmin, int _fcttype, double _pivmin, double _tau1,
                         double _tau2, double _theta, int _n, vector < _Int > &_ia_alu,
                         vector < _Int > &_ja_alu, vector < _Flt > &_a_alu,
                         vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                         vector < _Flt > &_a_u, double &_sclmin_att, double &_sclmax_att,
                         int &_nmodif, double &_eigmin_att, double &_eigmax_att)
   {

// Compute explicit symm scaling

      vector < _Flt > scl_U;

      CFct_impl < _Int, _Flt >::ComputeScalingSymm (_sclmin, _n, _ia_alu, _ja_alu, _a_alu,
                                                    scl_U, _sclmin_att, _sclmax_att);

// Perform explicit scaling

      int nza_alu = (int) _ia_alu[_n];

      vector < _Flt > a_scl_expl (nza_alu + 1);

      _Flt *p_a_alu = &_a_alu[0];
      _Flt *pa_scl_expl = &a_scl_expl[0];

      memcpy (pa_scl_expl, p_a_alu, (size_t) (nza_alu * sizeof (_Flt)));

      CFct_impl < _Int, _Flt >::MatrixScale (_n, scl_U, scl_U, _ia_alu, _ja_alu,
                                             a_scl_expl);

// Split

      vector < char >char_dummy;

      vector < _Int > ia_l;
      vector < _Int > ja_l;
      vector < _Flt > a_l;

      vector < _Int > ia_u;
      vector < _Int > ja_u;
      vector < _Flt > a_u;

      CFct_impl < _Int, _Flt >::SplitLU (false, _n, _ia_alu, _ja_alu, char_dummy,
                                         a_scl_expl, ia_l, ja_l, char_dummy, a_l, ia_u,
                                         ja_u, char_dummy, a_u);

      {
         vector < _Int > ia_l_dummy;
         vector < _Int > ja_l_dummy;
         vector < _Flt > a_l_dummy;
         vector < _Flt > a_dummy;
         ia_l.swap (ia_l_dummy);
         ja_l.swap (ja_l_dummy);
         a_l.swap (a_l_dummy);
         a_scl_expl.swap (a_dummy);
      }

// Perform fct

      CFct_impl < _Int, _Flt >::Ich2BlockIlu2 (_pivmin, _tau1, _tau2, _theta, _n, _n,
                                               ia_u, ja_u, a_u, _ia_u, _ja_u, _a_u,
                                               _nmodif, _eigmin_att, _eigmax_att);

// Perform explicit rescaling of compute triangular factors

      vector < _Flt > inv_scl_U (_n + 1);

      CFct_impl < _Int, _Flt >::InverseDiag (_n, scl_U, inv_scl_U);

      CFct_impl < _Int, _Flt >::RescaleU (_n, scl_U, inv_scl_U, _ia_u, _ja_u, _a_u);

   }

//
// Perform ILU2 point factorization of the block with future diagonal modification
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::Ilu2Block (int _sctype, int _nitersc, int _fcttype, double _pivmin,
                         double _tau1, double _tau2, double _theta, int _n,
                         vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                         vector < _Flt > &_a_alu, vector < _Int > &_ia_l,
                         vector < _Int > &_ja_l, vector < _Flt > &_a_l,
                         vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                         vector < _Flt > &_a_u, double &_sclmin_att, double &_sclmax_att,
                         int &_nmodif, double &_eigmin_att, double &_eigmax_att)
   {

      vector < char >char_dummy;

// Compute explicit scaling

      vector < _Flt > scl_L;
      vector < _Flt > scl_U;

      CFct_impl < _Int, _Flt >::ComputeScaling (_sctype, _nitersc, _n, _ia_alu, _ja_alu,
                                                _a_alu, scl_L, scl_U, _sclmin_att,
                                                _sclmax_att);

// Perform explicit scaling

      int nza_alu = (int) _ia_alu[_n];

      vector < _Flt > a_scl_expl (nza_alu + 1);

      _Flt *p_a_alu = &_a_alu[0];
      _Flt *pa_scl_expl = &a_scl_expl[0];

      memcpy (pa_scl_expl, p_a_alu, (size_t) (nza_alu * sizeof (_Flt)));

      CFct_impl < _Int, _Flt >::MatrixScale (_n, scl_L, scl_U, _ia_alu, _ja_alu,
                                             a_scl_expl);

// Split

      vector < _Int > ia_l;
      vector < _Int > ja_l;
      vector < _Flt > a_l;

      vector < _Int > ia_u;
      vector < _Int > ja_u;
      vector < _Flt > a_u;

      CFct_impl < _Int, _Flt >::SplitLU (false, _n, _ia_alu, _ja_alu, char_dummy,
                                         a_scl_expl, ia_l, ja_l, char_dummy, a_l, ia_u,
                                         ja_u, char_dummy, a_u);

      {
         vector < _Flt > a_dummy;
         a_scl_expl.swap (a_dummy);
      }

// Transpose L

      vector < _Int > ia_lt;
      vector < _Int > ja_lt;
      vector < _Flt > a_lt;

      CFct_impl < _Int, _Flt >::Transpose (false, _n, ia_l, ja_l, char_dummy, a_l, ia_lt,
                                           ja_lt, char_dummy, a_lt);

      {
         vector < _Int > ia_dummy;
         vector < _Int > ja_dummy;
         vector < _Flt > a_dummy;
         ia_l.swap (ia_dummy);
         ja_l.swap (ja_dummy);
         a_l.swap (a_dummy);
      }

// Combine into extended pairs

      vector < _Int > ia_alu_cnd;
      vector < _Int > ja_alu_cnd;
      vector < _Flt > a_alu_cnd;

      CFct_impl < _Int, _Flt >::CombinePairs (false, _n, ia_lt, ja_lt, char_dummy, a_lt,
                                              ia_u, ja_u, char_dummy, a_u, ia_alu_cnd,
                                              ja_alu_cnd, char_dummy, a_alu_cnd);

      {
         vector < _Int > ia_dummy;
         vector < _Int > ja_dummy;
         vector < _Flt > a_dummy;
         ia_lt.swap (ia_dummy);
         ja_lt.swap (ja_dummy);
         a_lt.swap (a_dummy);
         vector < _Int > ia_dummy1;
         vector < _Int > ja_dummy1;
         vector < _Flt > a_dummy1;
         ia_u.swap (ia_dummy1);
         ja_u.swap (ja_dummy1);
         a_u.swap (a_dummy1);
      }

// Perform fct

      vector < char >ja_char_alu_cnd;

      {
         int nzja_lu_cnd = (int) ia_alu_cnd[_n];
         ja_char_alu_cnd.resize (nzja_lu_cnd + 1);
         int i;
         for (i = 0; i < nzja_lu_cnd; i++)
            ja_char_alu_cnd[i] = 0;
      }

      vector < _Int > ia_lu_cnd;
      vector < _Int > ja_lu_cnd;
      vector < char >ja_char_lu_cnd;
      vector < _Flt > a_lu_cnd;

      int fcttype_sch = _fcttype;
      double tau2_sch = 0.0e0;

      CFct_impl < _Int, _Flt >::Ilu2BlockIlu2 (_fcttype, fcttype_sch, _pivmin, _tau1,
                                               _tau2, tau2_sch, _theta, _n, _n,
                                               ia_alu_cnd, ja_alu_cnd, ja_char_alu_cnd,
                                               a_alu_cnd, ia_lu_cnd, ja_lu_cnd,
                                               ja_char_lu_cnd, a_lu_cnd, _nmodif,
                                               _eigmin_att, _eigmax_att);

      {
         vector < _Int > ia_dummy;
         vector < _Int > ja_dummy;
         vector < char >ja_char_dummy;
         vector < _Flt > a_dummy;
         vector < char >ja_char1_dummy;
         ia_alu_cnd.swap (ia_dummy);
         ja_alu_cnd.swap (ja_dummy);
         ja_char_alu_cnd.swap (ja_char_dummy);
         a_alu_cnd.swap (a_dummy);
         ja_char_lu_cnd.swap (ja_char1_dummy);
      }

// Split pairs with filtering

      CFct_impl < _Int, _Flt >::SplitPairsFilter (false, _tau1, _n, ia_lu_cnd, ja_lu_cnd,
                                                  char_dummy, a_lu_cnd, _ia_l, _ja_l,
                                                  char_dummy, _a_l, _ia_u, _ja_u,
                                                  char_dummy, _a_u);

      {
         vector < _Int > ia_dummy;
         vector < _Int > ja_dummy;
         vector < _Flt > a_dummy;
         ia_lu_cnd.swap (ia_dummy);
         ja_lu_cnd.swap (ja_dummy);
         a_lu_cnd.swap (a_dummy);
      }

// Perform explicit rescaling of compute triangular factors

      vector < _Flt > inv_scl_L (_n + 1);
      vector < _Flt > inv_scl_U (_n + 1);

      CFct_impl < _Int, _Flt >::InverseDiag (_n, scl_L, inv_scl_L);
      CFct_impl < _Int, _Flt >::InverseDiag (_n, scl_U, inv_scl_U);

      CFct_impl < _Int, _Flt >::RescaleU (_n, scl_L, inv_scl_L, _ia_l, _ja_l, _a_l);
      CFct_impl < _Int, _Flt >::RescaleU (_n, scl_U, inv_scl_U, _ia_u, _ja_u, _a_u);

      {
         vector < _Flt > dl_dummy;
         vector < _Flt > du_dummy;
         vector < _Flt > idl_dummy;
         vector < _Flt > idu_dummy;
         scl_L.swap (dl_dummy);
         scl_U.swap (du_dummy);
         inv_scl_L.swap (idl_dummy);
         inv_scl_U.swap (idu_dummy);
      }

   }

//
// Compute symmectric scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::ComputeScalingSymm (double _sclmin, int _n, vector < _Int > &_ia_alu,
                                  vector < _Int > &_ja_alu, vector < _Flt > &_a_alu,
                                  vector < _Flt > &_sclU, double &_sclmin_att,
                                  double &_sclmax_att)
   {

      _sclU.resize (_n + 1);

// Simple diagonal based scaling

      double thresh_zero = 1.0e-13;

      _sclmin_att = 1.0e100;
      _sclmax_att = -1.0e100;

      int i, j, jj;

      _Flt diag;

      for (i = 0; i < _n; i++) {
         diag = 1.;
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jj = (int) _ja_alu[j];
            if (i == jj)
               diag = _a_alu[j];
         }
         if (diag < _sclmin_att)
            _sclmin_att = diag;
         if (diag > _sclmax_att)
            _sclmax_att = diag;
         if (diag < _sclmin)
            diag = (_Flt) _sclmin;
         if (diag < 0.0e0)
            diag = -diag;
         if (diag < thresh_zero)
            diag = (_Flt) thresh_zero;
         diag = (_Flt) (1. / sqrt (diag));
         _sclU[i] = diag;
      }

   }


//
// Compute scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::ComputeScaling (int _sctype, int _nitersc, int _n, vector < _Int > &_ia_alu,
                              vector < _Int > &_ja_alu, vector < _Flt > &_a_alu,
                              vector < _Flt > &_sclL, vector < _Flt > &_sclU,
                              double &_sclmin_att, double &_sclmax_att)
   {

      _sclL.resize (_n + 1);
      _sclU.resize (_n + 1);

// Simple diagonal based scaling

      _sclmin_att = 1.0e100;
      _sclmax_att = -1.0e100;

      int i, j, jj;

      _Flt diag;

      if (_sctype == -1) {
         diag = (_Flt) 1.;
         for (i = 0; i < _n; i++) {
            _sclL[i] = diag;
            _sclU[i] = diag;
         }
         _sclmin_att = diag;
         _sclmax_att = diag;
      } else if (_sctype == 0) {
         for (i = 0; i < _n; i++) {
            diag = 1.;
            for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
               jj = (int) _ja_alu[j];
               if (i == jj) {
                  diag = _a_alu[j];
               }
            }
            if (diag > 0.) {
               if (diag < _sclmin_att)
                  _sclmin_att = diag;
               if (diag > _sclmax_att)
                  _sclmax_att = diag;
               diag = (_Flt) (1. / sqrt (diag));
               _sclL[i] = diag;
               _sclU[i] = diag;
            } else {
               diag = -diag;
               if (diag < _sclmin_att)
                  _sclmin_att = diag;
               if (diag > _sclmax_att)
                  _sclmax_att = diag;
               diag = (_Flt) (1. / sqrt (diag));
               _sclL[i] = -diag;
               _sclU[i] = diag;
            }
         }
      } else {

// Iterative rows/columns balancing scaling

         _Flt aux, aux1;
         int iter;
         double sclRmin, sclRmax;
         double sclCmin, sclCmax;

         int nsmall_R = 0;
         int nsmall_C = 0;
         double thresh_small = 1.0e-8;

         sclRmin = 1.0e100;
         sclRmax = 0.0e0;
         sclCmin = 1.0e100;
         sclCmax = 0.0e0;

         for (i = 0; i < _n; i++) {
            aux = 0.;
            for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
               jj = (int) _ja_alu[j];
               aux1 = _a_alu[j];
               aux += aux1 * aux1;
            }
            if (aux < sclRmin)
               sclRmin = aux;
            if (aux > sclRmax)
               sclRmax = aux;
            if (aux < thresh_small) {
               nsmall_R++;
            }
//         aux = sqrt(aux);
//         aux = (_Flt)(1.0 / sqrt(aux));
            aux = (_Flt) (1.0 / aux);
            _sclL[i] = aux;
         }
         sclRmin = sqrt (sclRmin);
         sclRmax = sqrt (sclRmax);
         for (i = 0; i < _n; i++)
            _sclU[i] = 0.;
         for (i = 0; i < _n; i++) {
            for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
               jj = (int) _ja_alu[j];
               aux1 = _a_alu[j];
               _sclU[jj] += aux1 * aux1;
            }
         }
         for (i = 0; i < _n; i++) {
            aux = _sclU[i];
            if (aux < thresh_small) {
               nsmall_C++;
            }
            if (aux < sclCmin)
               sclCmin = aux;
            if (aux > sclCmax)
               sclCmax = aux;
         }
         sclCmin = sqrt (sclCmin);
         sclCmax = sqrt (sclCmax);
//      cout << "  thresh_small " << thresh_small << " nsmall_R = " << nsmall_R << "  nsmall_C = " << nsmall_C;
//      cout << " sclRmin = " << sclRmin << "  sclRmax = " << sclRmax << "  sclCmin = " << sclCmin << "  sclCmax = " << sclCmax << endl;
         _sclmin_att = sclRmin;
         if (sclCmin < _sclmin_att)
            _sclmin_att = sclCmin;
         _sclmax_att = sclRmax;
         if (sclCmax > _sclmax_att)
            _sclmax_att = sclCmax;
         for (iter = 0; iter < _nitersc; iter++) {
            for (i = 0; i < _n; i++)
               _sclU[i] = 0.;
            for (i = 0; i < _n; i++) {
               for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
                  jj = (int) _ja_alu[j];
                  aux1 = _a_alu[j];
                  _sclU[jj] += _sclL[i] * aux1 * aux1;
               }
            }
            for (i = 0; i < _n; i++) {
               _sclU[i] = (_Flt) (1. / _sclU[i]);
            }
            for (i = 0; i < _n; i++) {
               aux = 0.;
               for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
                  jj = (int) _ja_alu[j];
                  aux1 = _a_alu[j];
                  aux += aux1 * aux1 * _sclU[jj];
               }
               aux = (_Flt) (1.0 / aux);
               _sclL[i] = aux;
            }
         }
         for (i = 0; i < _n; i++)
            _sclL[i] = (_Flt) (sqrt (_sclL[i]));
         for (i = 0; i < _n; i++)
            _sclU[i] = (_Flt) (sqrt (_sclU[i]));
      }

   }

//
// Perform explicit scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::MatrixScale (int _n, vector < _Flt > &_sclL, vector < _Flt > &_sclU,
                           vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                           vector < _Flt > &_a_alu)
   {

      int i, j, jj;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jj = (int) _ja_alu[j];
            _a_alu[j] = _sclL[i] * _a_alu[j] * _sclU[jj];
         }
      }

   }

//
// Symmetrize sparsity
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::SymmetrizeSparsity (int _n, const vector < _Int > &_ia_alu,
                                  const vector < _Int > &_ja_alu,
                                  vector < _Int > &_ia_symm, vector < _Int > &_ja_symm)
   {

// Split

      vector < _Int > ia_l;
      vector < _Int > ja_l;

      vector < _Int > ia_u;
      vector < _Int > ja_u;

      CFct_impl < _Int, _Flt >::SplitLUSp (_n, _ia_alu, _ja_alu, ia_l, ja_l, ia_u, ja_u);

// Transpose L

      vector < _Int > ia_lt;
      vector < _Int > ja_lt;

      CFct_impl < _Int, _Flt >::TransposeSp (_n, ia_l, ja_l, ia_lt, ja_lt);

      {
         vector < _Int > ia_dummy;
         vector < _Int > ja_dummy;
         ia_l.swap (ia_dummy);
         ja_l.swap (ja_dummy);
      }

// Combine into extended pairs

      vector < _Int > ia_alu_cnd;
      vector < _Int > ja_alu_cnd;

      CFct_impl < _Int, _Flt >::CombineLUSp (_n, ia_lt, ja_lt, ia_u, ja_u, ia_alu_cnd,
                                             ja_alu_cnd);

      {
         vector < _Int > ia_dummy;
         vector < _Int > ja_dummy;
         ia_lt.swap (ia_dummy);
         ja_lt.swap (ja_dummy);
         vector < _Int > ia_dummy1;
         vector < _Int > ja_dummy1;
         ia_u.swap (ia_dummy1);
         ja_u.swap (ja_dummy1);
      }

// Transpose combined sparsity

      vector < _Int > ia_alu_cnd_t;
      vector < _Int > ja_alu_cnd_t;

      CFct_impl < _Int, _Flt >::TransposeSp (_n, ia_alu_cnd, ja_alu_cnd, ia_alu_cnd_t,
                                             ja_alu_cnd_t);

// Compute symmetrized matrix

      CFct_impl < _Int, _Flt >::CombineLUSp (_n, ia_alu_cnd_t, ja_alu_cnd_t, ia_alu_cnd,
                                             ja_alu_cnd, _ia_symm, _ja_symm);

   }

//
// Symmetrize sparsity
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::SymmetrizeSparsity (int _n, const vector < _Int > &_ia_alu,
                                  const vector < _Int > &_ja_alu,
                                  const vector < char >&_jachar_alu,
                                  vector < _Int > &_ia_symm, vector < _Int > &_ja_symm,
                                  vector < char >&_jachar_symm)
   {

// Split

      vector < _Int > ia_l;
      vector < _Int > ja_l;
      vector < char >jachar_l;

      vector < _Int > ia_u;
      vector < _Int > ja_u;
      vector < char >jachar_u;

      CFct_impl < _Int, _Flt >::SplitLUSp (_n, _ia_alu, _ja_alu, _jachar_alu, ia_l, ja_l,
                                           jachar_l, ia_u, ja_u, jachar_u);

// Transpose L

      vector < _Int > ia_lt;
      vector < _Int > ja_lt;
      vector < char >jachar_lt;

      CFct_impl < _Int, _Flt >::TransposeSp (_n, ia_l, ja_l, jachar_l, ia_lt, ja_lt,
                                             jachar_lt);

      {
         vector < _Int > ia_dummy;
         vector < _Int > ja_dummy;
         vector < char >jachar_dummy;
         ia_l.swap (ia_dummy);
         ja_l.swap (ja_dummy);
         jachar_l.swap (jachar_dummy);
      }

// Combine into extended pairs

      vector < _Int > ia_alu_cnd;
      vector < _Int > ja_alu_cnd;
      vector < char >jachar_alu_cnd;

      CFct_impl < _Int, _Flt >::CombineLUSp (_n, ia_lt, ja_lt, jachar_lt, ia_u, ja_u,
                                             jachar_u, ia_alu_cnd, ja_alu_cnd,
                                             jachar_alu_cnd);

      {
         vector < _Int > ia_dummy;
         vector < _Int > ja_dummy;
         vector < char >jachar_dummy;
         ia_lt.swap (ia_dummy);
         ja_lt.swap (ja_dummy);
         jachar_lt.swap (jachar_dummy);
         vector < _Int > ia_dummy1;
         vector < _Int > ja_dummy1;
         vector < char >jachar_dummy1;
         ia_u.swap (ia_dummy1);
         ja_u.swap (ja_dummy1);
         jachar_u.swap (jachar_dummy1);
      }

// Transpose combined sparsity

      vector < _Int > ia_alu_cnd_t;
      vector < _Int > ja_alu_cnd_t;
      vector < char >jachar_alu_cnd_t;

      CFct_impl < _Int, _Flt >::TransposeSp (_n, ia_alu_cnd, ja_alu_cnd, jachar_alu_cnd,
                                             ia_alu_cnd_t, ja_alu_cnd_t,
                                             jachar_alu_cnd_t);

// Compute symmetrized matrix

      CFct_impl < _Int, _Flt >::CombineLUSp (_n, ia_alu_cnd_t, ja_alu_cnd_t,
                                             jachar_alu_cnd_t, ia_alu_cnd, ja_alu_cnd,
                                             jachar_alu_cnd, _ia_symm, _ja_symm,
                                             _jachar_symm);

   }

//
// Split matrix data into L and U parts for sparsity only
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::SplitLUSp (int _n, const vector < _Int > &_ia_alu,
                         const vector < _Int > &_ja_alu, vector < _Int > &_ia_l,
                         vector < _Int > &_ja_l, vector < _Int > &_ia_u,
                         vector < _Int > &_ja_u)
   {

// Count number of elems

      int nzja_l = 0;
      int nzja_u = 0;

      int i, j, jj;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jj = (int) _ja_alu[j];
            if (jj <= i)
               nzja_l++;
            if (jj >= i)
               nzja_u++;
         }
      }

// Allocate and fill

      _ia_l.resize (_n + 1);
      _ja_l.resize (nzja_l + 1);
      _ia_u.resize (_n + 1);
      _ja_u.resize (nzja_u + 1);

      nzja_l = 0;
      nzja_u = 0;

      _ia_l[0] = 0;
      _ia_u[0] = 0;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jj = (int) _ja_alu[j];
            if (jj <= i) {
               _ja_l[nzja_l] = jj;
               nzja_l++;
            }
            if (jj >= i) {
               _ja_u[nzja_u] = jj;
               nzja_u++;
            }
         }
         _ia_l[i + 1] = nzja_l;
         _ia_u[i + 1] = nzja_u;
      }

   }

//
// Split matrix data into L and U parts for sparsity only
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::SplitLUSp (int _n, const vector < _Int > &_ia_alu,
                         const vector < _Int > &_ja_alu,
                         const vector < char >&_jachar_alu, vector < _Int > &_ia_l,
                         vector < _Int > &_ja_l, vector < char >&_jachar_l,
                         vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                         vector < char >&_jachar_u)
   {

// Count number of elems

      int nzja_l = 0;
      int nzja_u = 0;

      int i, j, jj;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jj = (int) _ja_alu[j];
            if (jj <= i)
               nzja_l++;
            if (jj >= i)
               nzja_u++;
         }
      }

// Allocate and fill

      _ia_l.resize (_n + 1);
      _ja_l.resize (nzja_l + 1);
      _jachar_l.resize (nzja_l + 1);
      _ia_u.resize (_n + 1);
      _ja_u.resize (nzja_u + 1);
      _jachar_u.resize (nzja_u + 1);

      nzja_l = 0;
      nzja_u = 0;

      _ia_l[0] = 0;
      _ia_u[0] = 0;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jj = (int) _ja_alu[j];
            if (jj <= i) {
               _ja_l[nzja_l] = jj;
               _jachar_l[nzja_l] = _jachar_alu[j];
               nzja_l++;
            }
            if (jj >= i) {
               _ja_u[nzja_u] = jj;
               _jachar_u[nzja_u] = _jachar_alu[j];
               nzja_u++;
            }
         }
         _ia_l[i + 1] = nzja_l;
         _ia_u[i + 1] = nzja_u;
      }

   }

//
// Split matrix data into L and U parts
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::SplitLU (bool _b_is_char, int _n, const vector < _Int > &_ia_alu,
                       const vector < _Int > &_ja_alu, const vector < char >&_jachar_alu,
                       const vector < _Flt > &_a_alu, vector < _Int > &_ia_l,
                       vector < _Int > &_ja_l, vector < char >&_jachar_l,
                       vector < _Flt > &_a_l, vector < _Int > &_ia_u,
                       vector < _Int > &_ja_u, vector < char >&_jachar_u,
                       vector < _Flt > &_a_u)
   {

// Count number of elems

      int nzja_l = 0;
      int nzja_u = 0;

      _Flt fzero = (_Flt) 0.0e0;

      int i, j, jj;

      bool b_found = false;

      int n_zero_diag = 0;

      for (i = 0; i < _n; i++) {
         b_found = false;
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jj = (int) _ja_alu[j];
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
         cout << " N_zero_diag = " << n_zero_diag << endl;

// Allocate and fill

      _ia_l.resize (_n + 1);
      _ja_l.resize (nzja_l + 1);
      _a_l.resize (nzja_l + 1);
      _ia_u.resize (_n + 1);
      _ja_u.resize (nzja_u + 1);
      _a_u.resize (nzja_u + 1);

      if (_b_is_char) {
         _jachar_l.resize (nzja_l + 1);
         _jachar_u.resize (nzja_u + 1);
      }

      nzja_l = 0;
      nzja_u = 0;

      _ia_l[0] = 0;
      _ia_u[0] = 0;

      for (i = 0; i < _n; i++) {
         b_found = false;
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jj = (int) _ja_alu[j];
            if (jj == i)
               b_found = true;
         }
         if (!b_found) {
            _ja_u[nzja_u] = i;
            if (_b_is_char) {
               _jachar_u[nzja_u] = 0;
            }
            _a_u[nzja_u] = fzero;
            nzja_u++;
         }
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jj = (int) _ja_alu[j];
            if (jj <= i) {
               _ja_l[nzja_l] = jj;
               if (_b_is_char) {
                  _jachar_l[nzja_l] = _jachar_alu[j];
               }
               _a_l[nzja_l] = _a_alu[j];
               nzja_l++;
            }
            if (jj >= i) {
               _ja_u[nzja_u] = jj;
               if (_b_is_char) {
                  _jachar_u[nzja_u] = _jachar_alu[j];
               }
               _a_u[nzja_u] = _a_alu[j];
               nzja_u++;
            }
         }
         if (!b_found) {
            _ja_l[nzja_l] = i;
            if (_b_is_char) {
               _jachar_l[nzja_l] = 0;
            }
            _a_l[nzja_l] = fzero;
            nzja_l++;
         }
         _ia_l[i + 1] = nzja_l;
         _ia_u[i + 1] = nzja_u;
      }

   }

//
// Transpose square matrix for sparsity only
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::TransposeSp (int _n, const vector < _Int > &_ia_a,
                           const vector < _Int > &_ja_a, vector < _Int > &_ia_at,
                           vector < _Int > &_ja_at)
   {

      int nzja_a = (int) _ia_a[_n];

// Allocate work array

      vector < _Int > iptr (_n + 1);

// Allocate tranposed data

      _ia_at.resize (_n + 1);
      _ja_at.resize (nzja_a + 1);

      int i, j, jj, k;

      for (i = 0; i <= _n; i++)
         _ia_at[i] = 0;

      for (i = 0; i < nzja_a; i++) {
         jj = (int) _ja_a[i];
         _ia_at[jj + 1]++;
      }

      for (i = 0; i < _n; i++)
         _ia_at[i + 1] = _ia_at[i] + _ia_at[i + 1];

      for (i = 0; i < _n; i++)
         iptr[i] = _ia_at[i];

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_a[i]; j < _ia_a[i + 1]; j++) {
            jj = (int) _ja_a[j];
            k = (int) iptr[jj];
            _ja_at[k] = i;
            iptr[jj]++;
         }
      }

   }

//
// Transpose square matrix for sparsity only
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::TransposeSp (int _n, const vector < _Int > &_ia_a,
                           const vector < _Int > &_ja_a, const vector < char >&_jachar_a,
                           vector < _Int > &_ia_at, vector < _Int > &_ja_at,
                           vector < char >&_jachar_at)
   {

      int nzja_a = (int) _ia_a[_n];

// Allocate work array

      vector < _Int > iptr (_n + 1);

// Allocate tranposed data

      _ia_at.resize (_n + 1);
      _ja_at.resize (nzja_a + 1);
      _jachar_at.resize (nzja_a + 1);

      int i, j, jj, k;

      for (i = 0; i <= _n; i++)
         _ia_at[i] = 0;

      for (i = 0; i < nzja_a; i++) {
         jj = (int) _ja_a[i];
         _ia_at[jj + 1]++;
      }

      for (i = 0; i < _n; i++)
         _ia_at[i + 1] = _ia_at[i] + _ia_at[i + 1];

      for (i = 0; i < _n; i++)
         iptr[i] = _ia_at[i];

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_a[i]; j < _ia_a[i + 1]; j++) {
            jj = (int) _ja_a[j];
            k = (int) iptr[jj];
            _ja_at[k] = i;
            _jachar_at[k] = _jachar_a[j];
            iptr[jj]++;
         }
      }

   }

//
// Transpose square matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::Transpose (bool _b_is_char, int _n, vector < _Int > &_ia_a,
                         vector < _Int > &_ja_a, vector < char >&_jachar_a,
                         vector < _Flt > &_a_a, vector < _Int > &_ia_at,
                         vector < _Int > &_ja_at, vector < char >&_jachar_at,
                         vector < _Flt > &_a_at)
   {

      int nzja_a = (int) _ia_a[_n];

// Allocate work array

      vector < _Int > iptr (_n + 1);

// Allocate tranposed data

      _ia_at.resize (_n + 1);
      _ja_at.resize (nzja_a + 1);
      if (_b_is_char) {
         _jachar_at.resize (nzja_a + 1);
      }
      _a_at.resize (nzja_a + 1);

      int i, j, jj, k;

      for (i = 0; i <= _n; i++)
         _ia_at[i] = 0;

      for (i = 0; i < nzja_a; i++) {
         jj = (int) _ja_a[i];
         _ia_at[jj + 1]++;
      }

      for (i = 0; i < _n; i++)
         _ia_at[i + 1] = _ia_at[i] + _ia_at[i + 1];

      for (i = 0; i < _n; i++)
         iptr[i] = _ia_at[i];

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_a[i]; j < _ia_a[i + 1]; j++) {
            jj = (int) _ja_a[j];
            k = (int) iptr[jj];
            _ja_at[k] = i;
            if (_b_is_char) {
               _jachar_at[k] = _jachar_a[j];
            }
            _a_at[k] = _a_a[j];
            iptr[jj]++;
         }
      }

   }

//
// Combine sparsity of L and U
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::CombineLUSp (int _n, vector < _Int > &_ia_l, vector < _Int > &_ja_l,
                           vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                           vector < _Int > &_ia_alu, vector < _Int > &_ja_alu)
   {

// Compute number of extended elems

      int nzja_ext = 0;

      int i, ipl, ipu, iendl, iendu, jj_l, jj_u;

      for (i = 0; i < _n; i++) {
         ipl = (int) _ia_l[i];
         ipu = (int) _ia_u[i];
         iendl = (int) _ia_l[i + 1] - 1;
         iendu = (int) _ia_u[i + 1] - 1;
         while (ipl <= iendl || ipu <= iendu) {
            if (ipl <= iendl && ipu <= iendu) {
               jj_l = (int) _ja_l[ipl];
               jj_u = (int) _ja_u[ipu];
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

      _ia_alu[0] = 0;
      nzja_ext = 0;

      for (i = 0; i < _n; i++) {
         ipl = (int) _ia_l[i];
         ipu = (int) _ia_u[i];
         iendl = (int) _ia_l[i + 1] - 1;
         iendu = (int) _ia_u[i + 1] - 1;
         while (ipl <= iendl || ipu <= iendu) {
            if (ipl <= iendl && ipu <= iendu) {
               jj_l = (int) _ja_l[ipl];
               jj_u = (int) _ja_u[ipu];
               if (jj_l == jj_u) {
                  _ja_alu[nzja_ext] = jj_l;
                  nzja_ext++;
                  ipl++;
                  ipu++;
               } else if (jj_l < jj_u) {
                  _ja_alu[nzja_ext] = jj_l;
                  nzja_ext++;
                  ipl++;
               } else {
                  _ja_alu[nzja_ext] = jj_u;
                  nzja_ext++;
                  ipu++;
               }
            } else if (ipl <= iendl) {
               _ja_alu[nzja_ext] = _ja_l[ipl];
               nzja_ext++;
               ipl++;
            } else {
               _ja_alu[nzja_ext] = _ja_u[ipu];
               nzja_ext++;
               ipu++;
            }
         }
         _ia_alu[i + 1] = nzja_ext;
      }

   }

//
// Combine sparsity of L and U
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::CombineLUSp (int _n, vector < _Int > &_ia_l, vector < _Int > &_ja_l,
                           vector < char >&_jachar_l, vector < _Int > &_ia_u,
                           vector < _Int > &_ja_u, vector < char >&_jachar_u,
                           vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                           vector < char >&_jachar_alu)
   {

// Compute number of extended elems

      int nzja_ext = 0;

      int i, ipl, ipu, iendl, iendu, jj_l, jj_u;

      for (i = 0; i < _n; i++) {
         ipl = (int) _ia_l[i];
         ipu = (int) _ia_u[i];
         iendl = (int) _ia_l[i + 1] - 1;
         iendu = (int) _ia_u[i + 1] - 1;
         while (ipl <= iendl || ipu <= iendu) {
            if (ipl <= iendl && ipu <= iendu) {
               jj_l = (int) _ja_l[ipl];
               jj_u = (int) _ja_u[ipu];
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
      _jachar_alu.resize (nzja_ext + 1);

      _ia_alu[0] = 0;
      nzja_ext = 0;

      for (i = 0; i < _n; i++) {
         ipl = (int) _ia_l[i];
         ipu = (int) _ia_u[i];
         iendl = (int) _ia_l[i + 1] - 1;
         iendu = (int) _ia_u[i + 1] - 1;
         while (ipl <= iendl || ipu <= iendu) {
            if (ipl <= iendl && ipu <= iendu) {
               jj_l = (int) _ja_l[ipl];
               jj_u = (int) _ja_u[ipu];
               if (jj_l == jj_u) {
                  _ja_alu[nzja_ext] = jj_l;
                  _jachar_alu[nzja_ext] = _jachar_l[ipl];
                  nzja_ext++;
                  ipl++;
                  ipu++;
               } else if (jj_l < jj_u) {
                  _ja_alu[nzja_ext] = jj_l;
                  _jachar_alu[nzja_ext] = _jachar_l[ipl];
                  nzja_ext++;
                  ipl++;
               } else {
                  _ja_alu[nzja_ext] = jj_u;
                  _jachar_alu[nzja_ext] = _jachar_u[ipu];
                  nzja_ext++;
                  ipu++;
               }
            } else if (ipl <= iendl) {
               _ja_alu[nzja_ext] = _ja_l[ipl];
               _jachar_alu[nzja_ext] = _jachar_l[ipl];
               nzja_ext++;
               ipl++;
            } else {
               _ja_alu[nzja_ext] = _ja_u[ipu];
               _jachar_alu[nzja_ext] = _jachar_u[ipu];
               nzja_ext++;
               ipu++;
            }
         }
         _ia_alu[i + 1] = nzja_ext;
      }

   }

//
// Combine L and U
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::CombineLU (bool _b_is_char, int _n, const vector < _Int > &_ia_l,
                         const vector < _Int > &_ja_l, const vector < char >&_jachar_l,
                         const vector < _Flt > &_a_l, const vector < _Int > &_ia_u,
                         const vector < _Int > &_ja_u, const vector < char >&_jachar_u,
                         const vector < _Flt > &_a_u, vector < _Int > &_ia_alu,
                         vector < _Int > &_ja_alu, vector < char >&_jachar_alu,
                         vector < _Flt > &_a_alu)
   {

// Compute number of extended elems

      int nzja_ext = 0;

      int i, ipl, ipu, iendl, iendu, jj_l, jj_u;

      for (i = 0; i < _n; i++) {
         ipl = (int) _ia_l[i];
         ipu = (int) _ia_u[i];
         iendl = (int) _ia_l[i + 1] - 1;
         iendu = (int) _ia_u[i + 1] - 1;
         while (ipl <= iendl || ipu <= iendu) {
            if (ipl <= iendl && ipu <= iendu) {
               jj_l = (int) _ja_l[ipl];
               jj_u = (int) _ja_u[ipu];
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
      _a_alu.resize (nzja_ext + 1);

      _ia_alu[0] = 0;
      nzja_ext = 0;

      for (i = 0; i < _n; i++) {
         ipl = (int) _ia_l[i];
         ipu = (int) _ia_u[i];
         iendl = (int) _ia_l[i + 1] - 1;
         iendu = (int) _ia_u[i + 1] - 1;
         while (ipl <= iendl || ipu <= iendu) {
            if (ipl <= iendl && ipu <= iendu) {
               jj_l = (int) _ja_l[ipl];
               jj_u = (int) _ja_u[ipu];
               if (jj_l == jj_u) {
                  _ja_alu[nzja_ext] = jj_l;
                  if (_b_is_char) {
                     _jachar_alu[nzja_ext] = _jachar_l[ipl];
                  }
                  _a_alu[nzja_ext] = _a_l[ipl];
                  nzja_ext++;
                  ipl++;
                  ipu++;
               } else if (jj_l < jj_u) {
                  _ja_alu[nzja_ext] = jj_l;
                  if (_b_is_char) {
                     _jachar_alu[nzja_ext] = _jachar_l[ipl];
                  }
                  _a_alu[nzja_ext] = _a_l[ipl];
                  nzja_ext++;
                  ipl++;
               } else {
                  _ja_alu[nzja_ext] = jj_u;
                  if (_b_is_char) {
                     _jachar_alu[nzja_ext] = _jachar_u[ipu];
                  }
                  _a_alu[nzja_ext] = _a_u[ipu];
                  nzja_ext++;
                  ipu++;
               }
            } else if (ipl <= iendl) {
               _ja_alu[nzja_ext] = _ja_l[ipl];
               if (_b_is_char) {
                  _jachar_alu[nzja_ext] = _jachar_l[ipl];
               }
               _a_alu[nzja_ext] = _a_l[ipl];
               nzja_ext++;
               ipl++;
            } else {
               _ja_alu[nzja_ext] = _ja_u[ipu];
               if (_b_is_char) {
                  _jachar_alu[nzja_ext] = _jachar_u[ipu];
               }
               _a_alu[nzja_ext] = _a_u[ipu];
               nzja_ext++;
               ipu++;
            }
         }
         _ia_alu[i + 1] = nzja_ext;
      }

   }

//
// Combine L and U data into extended pairs
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::CombinePairs (bool _b_is_char, int _n, vector < _Int > &_ia_l,
                            vector < _Int > &_ja_l, vector < char >&_jachar_l,
                            vector < _Flt > &_a_l, vector < _Int > &_ia_u,
                            vector < _Int > &_ja_u, vector < char >&_jachar_u,
                            vector < _Flt > &_a_u, vector < _Int > &_ia_alu,
                            vector < _Int > &_ja_alu, vector < char >&_jachar_alu,
                            vector < _Flt > &_a_alu)
   {

// Compute number of extended elems

      int nzja_ext = 0;

      int i, ipl, ipu, iendl, iendu, jj_l, jj_u;

      for (i = 0; i < _n; i++) {
         ipl = (int) _ia_l[i];
         ipu = (int) _ia_u[i];
         iendl = (int) _ia_l[i + 1] - 1;
         iendu = (int) _ia_u[i + 1] - 1;
         while (ipl <= iendl || ipu <= iendu) {
            if (ipl <= iendl && ipu <= iendu) {
               jj_l = (int) _ja_l[ipl];
               jj_u = (int) _ja_u[ipu];
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
      _a_alu.resize (nzja_ext * 2 + 1);

      _ia_alu[0] = 0;
      nzja_ext = 0;

      for (i = 0; i < _n; i++) {
         ipl = (int) _ia_l[i];
         ipu = (int) _ia_u[i];
         iendl = (int) _ia_l[i + 1] - 1;
         iendu = (int) _ia_u[i + 1] - 1;
         while (ipl <= iendl || ipu <= iendu) {
            if (ipl <= iendl && ipu <= iendu) {
               jj_l = (int) _ja_l[ipl];
               jj_u = (int) _ja_u[ipu];
               if (jj_l == jj_u) {
                  _ja_alu[nzja_ext] = jj_l;
                  if (_b_is_char) {
                     _jachar_alu[nzja_ext] = _jachar_l[ipl];
                  }
                  _a_alu[nzja_ext * 2] = _a_l[ipl];
                  _a_alu[nzja_ext * 2 + 1] = _a_u[ipu];
                  nzja_ext++;
                  ipl++;
                  ipu++;
               } else if (jj_l < jj_u) {
                  _ja_alu[nzja_ext] = jj_l;
                  if (_b_is_char) {
                     _jachar_alu[nzja_ext] = _jachar_l[ipl];
                  }
                  _a_alu[nzja_ext * 2] = _a_l[ipl];
                  _a_alu[nzja_ext * 2 + 1] = 0.;
                  nzja_ext++;
                  ipl++;
               } else {
                  _ja_alu[nzja_ext] = jj_u;
                  if (_b_is_char) {
                     _jachar_alu[nzja_ext] = _jachar_u[ipu];
                  }
                  _a_alu[nzja_ext * 2] = 0.;
                  _a_alu[nzja_ext * 2 + 1] = _a_u[ipu];
                  nzja_ext++;
                  ipu++;
               }
            } else if (ipl <= iendl) {
               _ja_alu[nzja_ext] = _ja_l[ipl];
               if (_b_is_char) {
                  _jachar_alu[nzja_ext] = _jachar_l[ipl];
               }
               _a_alu[nzja_ext * 2] = _a_l[ipl];
               _a_alu[nzja_ext * 2 + 1] = 0.;
               nzja_ext++;
               ipl++;
            } else {
               _ja_alu[nzja_ext] = _ja_u[ipu];
               if (_b_is_char) {
                  _jachar_alu[nzja_ext] = _jachar_u[ipu];
               }
               _a_alu[nzja_ext * 2] = 0.;
               _a_alu[nzja_ext * 2 + 1] = _a_u[ipu];
               nzja_ext++;
               ipu++;
            }
         }
         _ia_alu[i + 1] = nzja_ext;
      }

   }

//
// Split pairs fct data into L and U parts with post filtering
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::SplitPairsFilter (bool _b_is_char, double _tau1, int _n,
                                vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                                vector < char >&_jachar_alu, vector < _Flt > &_a_alu,
                                vector < _Int > &_ia_l, vector < _Int > &_ja_l,
                                vector < char >&_jachar_l, vector < _Flt > &_a_l,
                                vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                                vector < char >&_jachar_u, vector < _Flt > &_a_u)
   {

      int nzja_alu = (int) _ia_alu[_n];

// Allocate and fill

      _ia_l.resize (_n + 1);
      _ja_l.resize (nzja_alu + 1);
      _a_l.resize (nzja_alu + 1);
      _ia_u.resize (_n + 1);
      _ja_u.resize (nzja_alu + 1);
      _a_u.resize (nzja_alu + 1);

      if (_b_is_char) {
         _jachar_l.resize (nzja_alu + 1);
         _jachar_u.resize (nzja_alu + 1);
      }

      int nzja_l = 0;
      int nzja_u = 0;

      int i, j, jj;
      double auxL, auxU;

      _ia_l[0] = 0;
      _ia_u[0] = 0;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_alu[i]; j < _ia_alu[i + 1]; j++) {
            jj = (int) _ja_alu[j];
            auxL = _a_alu[j * 2];
            auxU = _a_alu[j * 2 + 1];
            if (auxL < 0.)
               auxL = -auxL;
            if (auxU < 0.)
               auxU = -auxU;
            if (jj == i || auxL >= _tau1) {
               _ja_l[nzja_l] = jj;
               if (_b_is_char) {
                  _jachar_l[nzja_l] = _jachar_alu[j];
               }
               _a_l[nzja_l] = _a_alu[j * 2];
               nzja_l++;
            }
            if (jj == i || auxU >= _tau1) {
               _ja_u[nzja_u] = jj;
               if (_b_is_char) {
                  _jachar_u[nzja_u] = _jachar_alu[j];
               }
               _a_u[nzja_u] = _a_alu[j * 2 + 1];
               nzja_u++;
            }
         }
         _ia_l[i + 1] = nzja_l;
         _ia_u[i + 1] = nzja_u;
      }

   }

//
// Combine L and U data into extended pairs
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::CombineRowsLU (int _n, vector < _Int > &_ia_l, vector < _Int > &_ja_l,
                             vector < _Flt > &_a_l, vector < _Int > &_ia_u,
                             vector < _Int > &_ja_u, vector < _Flt > &_a_u,
                             vector < _Int > &_ia_alu, vector < _Int > &_ida_alu,
                             vector < _Int > &_ja_alu, vector < _Flt > &_a_alu)
   {

// Compute number of extended elems

      int nzja_ext = (int) (_ia_l[_n] + _ia_u[_n]);

// Count number of elems

      _ia_alu.resize (_n + 1);
      _ida_alu.resize (_n + 1);
      _ja_alu.resize (nzja_ext + 1);
      _a_alu.resize (nzja_ext + 1);

      _ia_alu[0] = 0;

      int i, j;

      nzja_ext = 0;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_l[i]; j < _ia_l[i + 1]; j++) {
            _ja_alu[nzja_ext] = _ja_l[j];
            _a_alu[nzja_ext] = _a_l[j];
            nzja_ext++;
         }
         _ida_alu[i] = nzja_ext;
         for (j = (int) _ia_u[i]; j < _ia_u[i + 1]; j++) {
            _ja_alu[nzja_ext] = _ja_u[j];
            _a_alu[nzja_ext] = _a_u[j];
            nzja_ext++;
         }
         _ia_alu[i + 1] = nzja_ext;
      }

   }

//
// Compute inverse scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::InverseDiag (int _n, vector < _Flt > &_sclU, vector < _Flt > &_invsclU)
   {

      int i;

      for (i = 0; i < _n; i++) {
         _invsclU[i] = (_Flt) (1. / _sclU[i]);
      }

   }

//
// Rescale factor back
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::RescaleU (int _n, vector < _Flt > &_sclU, vector < _Flt > &_invsclU,
                        vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                        vector < _Flt > &_a_u)
   {

      int i, j, jj;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_u[i]; j < _ia_u[i + 1]; j++) {
            jj = (int) _ja_u[j];
            if (jj != i) {
               _a_u[j] *= _invsclU[jj];
            } else {
               _a_u[j] *= _sclU[i];
            }
         }
      }

   }

//
// Perform ICH2 point factorization of the block with future diagonal modification (no structural control)
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::Ich2BlockIlu2 (double _pivmin, double _tau1, double _tau2, double _theta,
                             int _n, int _n_ini, vector < _Int > &_ia_au,
                             vector < _Int > &_ja_au, vector < _Flt > &_a_au,
                             vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                             vector < _Flt > &_a_u, int &_nmodif, double &_eigmin_att,
                             double &_eigmax_att)
   {

// Prepare mask arrays

      vector < _Int > ibegm (_n + 1);
      vector < _Int > madj (_n + 1);
      vector < _Int > iv (_n + 1);
      vector < _Int > imaskc (_n + 1);
      vector < _Int > listc (_n + 1);

      vector < _Flt > fmaskc (_n + 1);
      vector < _Flt > a_dia (_n + 1);

      _Int *plistc = &listc[0];

      int i;

      for (i = 0; i <= _n; i++)
         ibegm[i] = -1;
      for (i = 0; i <= _n; i++)
         madj[i] = -1;
      for (i = 0; i <= _n; i++)
         imaskc[i] = -1;
      for (i = 0; i < _n; i++)
         a_dia[i] = 0.;

      double eigmin_att = FLT_MAX;
      double eigmax_att = -FLT_MAX;

      double thresh_zero = 1.0e-13;

      _Flt fzero = 0.;
      _Flt fone = 1.;

      int j, jj, k, kcolmn, irwprv, nlistcnew, ind, ind1, irow, irwpr1, jcolmn;
      double auxL;
      double aux1, dfnorm;

      int icycle_int = -1;

      _ia_u.resize (0);
      _ja_u.resize (0);
      _a_u.resize (0);

      _ia_u.resize (_n + 1);

      _ja_u.reserve (_n + 1);
      _a_u.reserve (_n * 2 + 1);

      _ia_u[0] = 0;

      int nzja_u = 0;

      _nmodif = 0;

      for (i = 0; i < _n; i++) {

         irow = i;

// Init current row

         int nlistcloc = 0;

         icycle_int++;

         if (i < _n_ini) {

            for (k = (int) _ia_au[i]; k < _ia_au[i + 1]; k++) {
               kcolmn = (int) _ja_au[k];
               listc[nlistcloc] = kcolmn;
               nlistcloc++;
               fmaskc[kcolmn] = _a_au[k];
               imaskc[kcolmn] = icycle_int;
            }

            if (imaskc[i] != icycle_int) {
               listc[nlistcloc] = i;
               nlistcloc++;
               fmaskc[i] = fzero;
               imaskc[i] = icycle_int;
            }

         } else {
            listc[nlistcloc] = i;
            nlistcloc++;
            fmaskc[i] = fzero;
            imaskc[i] = icycle_int;
         }

// Update current row

         irwprv = (int) ibegm[irow];

         while (irwprv != -1) {
            irwpr1 = (int) madj[irwprv];
            if (iv[irwprv] < _ia_u[irwprv + 1]) {

               j = (int) iv[irwprv];

               jj = (int) _ja_u[j];

               auxL = _a_u[j];

               jcolmn = jj;
               if (jcolmn < 0)
                  jcolmn = -jcolmn - 1;

               if (jj >= 0) {
                  for (k = (int) iv[irwprv]; k < _ia_u[irwprv + 1]; k++) {
                     kcolmn = (int) _ja_u[k];
                     if (kcolmn < 0)
                        kcolmn = -kcolmn - 1;
                     if (imaskc[kcolmn] != icycle_int) {
                        listc[nlistcloc] = kcolmn;
                        nlistcloc++;
                        fmaskc[kcolmn] = fzero;
                        imaskc[kcolmn] = icycle_int;
                     }
                     fmaskc[kcolmn] -= (_Flt) (auxL * _a_u[k]);
                  }
               } else {
                  for (k = (int) iv[irwprv]; k < _ia_u[irwprv + 1]; k++) {
                     kcolmn = (int) _ja_u[k];
                     if (kcolmn >= 0) {
                        if (imaskc[kcolmn] != icycle_int) {
                           listc[nlistcloc] = kcolmn;
                           nlistcloc++;
                           fmaskc[kcolmn] = fzero;
                           imaskc[kcolmn] = icycle_int;
                        }
                        fmaskc[kcolmn] -= (_Flt) (auxL * _a_u[k]);
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

            if (iv[irwprv] >= _ia_u[irwprv + 1] - 1) {
               madj[irwprv] = -1;
            } else {
               j = (int) iv[irwprv] + 1;
               jcolmn = (int) _ja_u[j];
               if (jcolmn < 0)
                  jcolmn = -jcolmn - 1;
               madj[irwprv] = ibegm[jcolmn];
               ibegm[jcolmn] = irwprv;
            }

            iv[irwprv]++;

            irwprv = irwpr1;

         }

// Perform filtering of the data and modify diagonal

         if (i < _n_ini) {

            nlistcnew = 0;

            for (j = 0; j < nlistcloc; j++) {
               ind = (int) listc[j];
               if (ind != irow) {
                  aux1 = fmaskc[ind];
                  if (aux1 < fzero)
                     aux1 = -aux1;
                  dfnorm = aux1;
                  if (dfnorm < _tau2) {
                     a_dia[irow] += (_Flt) (dfnorm);
                     a_dia[ind] += (_Flt) (dfnorm);
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

         if (fmaskc[irow] > fzero) {
            fmaskc[irow] += (_Flt) (a_dia[irow] * _theta);
         } else {
            fmaskc[irow] -= (_Flt) (a_dia[irow] * _theta);
         }

// Factorize current row and split into first and second order

         sort (plistc, plistc + nlistcloc);

         if (i < _n_ini) {

            aux1 = fmaskc[irow];

            if (aux1 < eigmin_att)
               eigmin_att = aux1;
            if (aux1 > eigmax_att)
               eigmax_att = aux1;
            if (aux1 < _pivmin) {
               aux1 = _pivmin;
               _nmodif++;
            } else if (aux1 < fzero) {
               aux1 = -aux1;
               _nmodif++;
            }
            if (aux1 < thresh_zero)
               aux1 = (_Flt) thresh_zero;
            aux1 = sqrt (aux1);
            aux1 = fone / aux1;
            auxL = aux1;
            fmaskc[irow] = (_Flt) (auxL);
         }

         if (i < _n_ini) {
            for (j = 0; j < nlistcloc; j++) {
               ind = (int) listc[j];
               if (ind != irow) {
                  fmaskc[ind] *= (_Flt) (auxL);
                  aux1 = fmaskc[ind];
                  if (aux1 < fzero)
                     aux1 = -aux1;
                  dfnorm = aux1;
                  if (dfnorm < _tau1) {
                     listc[j] = (_Int) (-ind - 1);
                  }
               }
            }
         }
// Store computed row elems

         for (j = 0; j < nlistcloc; j++) {
            ind = (int) plistc[j];
            _ja_u.push_back (ind);
            ind1 = ind;
            if (ind1 < 0)
               ind1 = -ind1 - 1;
            _a_u.push_back (fmaskc[ind1]);
         }

         nzja_u += nlistcloc;

         _ia_u[i + 1] = nzja_u;

// Add current row into the transposed structures

         if (nlistcloc > 1) {
            if (i < _n_ini) {
               iv[irow] = _ia_u[i] + 1;
               ind = (int) listc[1];
               if (ind < 0)
                  ind = -ind - 1;
               madj[irow] = ibegm[ind];
               ibegm[ind] = (_Int) irow;
            } else {
               iv[irow] = _ia_u[i + 1] + 1;
            }
         }

      }

      _eigmin_att = eigmin_att;
      _eigmax_att = eigmax_att;

// Condense the result (filter second order elems)

      int nzja_new = 0;

      for (i = 0; i < nzja_u; i++) {
         if (_ja_u[i] >= 0)
            nzja_new++;
      }

      vector < _Int > ju_cnd (nzja_new + 1);
      vector < _Flt > u_cnd (nzja_new + 1);

      nzja_new = 0;

      listc[0] = 0;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_u[i]; j < _ia_u[i + 1]; j++) {
            if (_ja_u[j] >= 0) {
               ju_cnd[nzja_new] = _ja_u[j];
               u_cnd[nzja_new] = _a_u[j];
               nzja_new++;
            }
         }
         listc[i + 1] = nzja_new;
      }

      for (i = 0; i <= _n; i++)
         _ia_u[i] = listc[i];

      _ja_u.swap (ju_cnd);
      _a_u.swap (u_cnd);

   }

//
// Perform ILU2 point factorization of the block with future diagonal modification (no structural control)
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::Ilu2BlockIlu2 (double _pivmin, double _tau1, double _tau2, double _theta,
                             int _n, int _n_ini, vector < _Int > &_ia_alu,
                             vector < _Int > &_ja_alu, vector < _Flt > &_a_alu,
                             vector < _Int > &_ia_lu, vector < _Int > &_ja_lu,
                             vector < _Flt > &_a_lu, int &_nmodif, double &_eigmin_att,
                             double &_eigmax_att)
   {

// Prepare mask arrays

      vector < _Int > ibegm (_n + 1);
      vector < _Int > madj (_n + 1);
      vector < _Int > iv (_n + 1);
      vector < _Int > imaskc (_n + 1);
      vector < _Int > listc (_n + 1);

      vector < _Flt > fmaskc (2 * _n + 1);
      vector < _Flt > a_dia (_n + 1);

      _Int *plistc = &listc[0];

      int i;

      for (i = 0; i <= _n; i++)
         ibegm[i] = -1;
      for (i = 0; i <= _n; i++)
         madj[i] = -1;
      for (i = 0; i <= _n; i++)
         imaskc[i] = -1;
      for (i = 0; i < _n; i++)
         a_dia[i] = 0.;

      double eigmin_att = FLT_MAX;
      double eigmax_att = -FLT_MAX;

      _Flt fzero = 0.;
      _Flt fone = 1.;

      int j, jj, k, kcolmn, irwprv, nlistcnew, ind, ind1, irow, irwpr1, jcolmn;
      double auxL, auxU;
      double aux1, aux2, dfnorm;

      int icycle_int = -1;

      _ia_lu.resize (_n + 1);

      _ja_lu.reserve (_n + 1);
      _a_lu.reserve (_n * 2 + 1);

      _ia_lu[0] = 0;

      int nzja_lu = 0;

      _nmodif = 0;

      for (i = 0; i < _n; i++) {

         irow = i;

// Init current row

         int nlistcloc = 0;

         icycle_int++;

         if (i < _n_ini) {

            for (k = (int) _ia_alu[i]; k < _ia_alu[i + 1]; k++) {
               kcolmn = (int) _ja_alu[k];
               listc[nlistcloc] = kcolmn;
               nlistcloc++;
               fmaskc[kcolmn * 2] = _a_alu[k * 2];
               fmaskc[kcolmn * 2 + 1] = _a_alu[k * 2 + 1];
               imaskc[kcolmn] = icycle_int;
            }

            if (imaskc[i] != icycle_int) {
               listc[nlistcloc] = i;
               nlistcloc++;
               fmaskc[2 * i] = fzero;
               fmaskc[2 * i + 1] = fzero;
               imaskc[i] = icycle_int;
            }

         } else {
            listc[nlistcloc] = i;
            nlistcloc++;
            fmaskc[i * 2] = fzero;
            fmaskc[i * 2 + 1] = fzero;
            imaskc[i] = icycle_int;
         }

// Update current row

         irwprv = (int) ibegm[irow];

         while (irwprv != -1) {
            irwpr1 = (int) madj[irwprv];
            if (iv[irwprv] < _ia_lu[irwprv + 1]) {

               j = (int) iv[irwprv];

               jj = (int) _ja_lu[j];

               auxL = _a_lu[j * 2];
               auxU = _a_lu[j * 2 + 1];

               jcolmn = jj;
               if (jcolmn < 0)
                  jcolmn = -jcolmn - 1;

               if (jj >= 0) {
                  for (k = (int) iv[irwprv]; k < _ia_lu[irwprv + 1]; k++) {
                     kcolmn = (int) _ja_lu[k];
                     if (kcolmn < 0)
                        kcolmn = -kcolmn - 1;
                     if (imaskc[kcolmn] != icycle_int) {
                        listc[nlistcloc] = kcolmn;
                        nlistcloc++;
                        fmaskc[kcolmn * 2] = fzero;
                        fmaskc[kcolmn * 2 + 1] = fzero;
                        imaskc[kcolmn] = icycle_int;
                     }
                     fmaskc[kcolmn * 2] -= (_Flt) (auxU * _a_lu[k * 2]);
                     fmaskc[kcolmn * 2 + 1] -= (_Flt) (auxL * _a_lu[k * 2 + 1]);
                  }
               } else {
                  for (k = (int) iv[irwprv]; k < _ia_lu[irwprv + 1]; k++) {
                     kcolmn = (int) _ja_lu[k];
                     if (kcolmn >= 0) {
                        if (imaskc[kcolmn] != icycle_int) {
                           listc[nlistcloc] = kcolmn;
                           nlistcloc++;
                           fmaskc[kcolmn * 2] = fzero;
                           fmaskc[kcolmn * 2 + 1] = fzero;
                           imaskc[kcolmn] = icycle_int;
                        }
                        fmaskc[kcolmn * 2] -= (_Flt) (auxU * _a_lu[k * 2]);
                        fmaskc[kcolmn * 2 + 1] -= (_Flt) (auxL * _a_lu[k * 2 + 1]);
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

         if (i < _n_ini) {

            nlistcnew = 0;

            for (j = 0; j < nlistcloc; j++) {
               ind = (int) listc[j];
               if (ind != irow) {
                  aux1 = fmaskc[ind * 2];
                  aux2 = fmaskc[ind * 2 + 1];
                  if (aux1 < fzero)
                     aux1 = -aux1;
                  if (aux2 < fzero)
                     aux2 = -aux2;
                  dfnorm = (aux1 > aux2) ? aux1 : aux2;
                  if (dfnorm < _tau2) {
                     a_dia[irow] += (_Flt) (dfnorm);
                     a_dia[ind] += (_Flt) (dfnorm);
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

         if (fmaskc[irow * 2] > fzero) {
            fmaskc[irow * 2] += (_Flt) (a_dia[irow] * _theta);
         } else {
            fmaskc[irow * 2] -= (_Flt) (a_dia[irow] * _theta);
         }

// Factorize current row and split into first and second order

         sort (plistc, plistc + nlistcloc);

         if (i < _n_ini) {

            aux1 = fmaskc[irow * 2];

            if (aux1 > 0.) {
               if (aux1 < eigmin_att)
                  eigmin_att = aux1;
               if (aux1 > eigmax_att)
                  eigmax_att = aux1;
               if (aux1 < _pivmin) {
                  aux1 = _pivmin;
                  _nmodif++;
               }
               aux1 = sqrt (aux1);
               aux1 = fone / aux1;
               auxL = aux1;
               auxU = aux1;
            } else {
               aux1 = -aux1;
               if (aux1 < eigmin_att)
                  eigmin_att = aux1;
               if (aux1 > eigmax_att)
                  eigmax_att = aux1;
               if (aux1 < _pivmin) {
                  aux1 = _pivmin;
                  _nmodif++;
               }
               aux1 = sqrt (aux1);
               aux1 = fone / aux1;
               auxL = -aux1;
               auxU = aux1;
            }
            fmaskc[irow * 2] = (_Flt) (auxL);
            fmaskc[irow * 2 + 1] = (_Flt) (auxU);
         }

         if (i < _n_ini) {
            for (j = 0; j < nlistcloc; j++) {
               ind = (int) listc[j];
               if (ind != irow) {
                  fmaskc[ind * 2] *= (_Flt) (auxU);
                  fmaskc[ind * 2 + 1] *= (_Flt) (auxL);
                  aux1 = fmaskc[ind * 2];
                  aux2 = fmaskc[ind * 2 + 1];
                  if (aux1 < fzero)
                     aux1 = -aux1;
                  if (aux2 < fzero)
                     aux2 = -aux2;
                  dfnorm = (aux1 > aux2) ? aux1 : aux2;
                  if (dfnorm < _tau1) {
                     listc[j] = (_Int) (-ind - 1);
                  }
               }
            }
         }
// Store computed row elems

         for (j = 0; j < nlistcloc; j++) {
            ind = (int) plistc[j];
            _ja_lu.push_back (ind);
            ind1 = ind;
            if (ind1 < 0)
               ind1 = -ind1 - 1;
            _a_lu.push_back (fmaskc[ind1 * 2]);
            _a_lu.push_back (fmaskc[ind1 * 2 + 1]);
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
      vector < _Flt > lu_cnd (2 * nzja_new + 1);

      nzja_new = 0;

      listc[0] = 0;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_lu[i]; j < _ia_lu[i + 1]; j++) {
            if (_ja_lu[j] >= 0) {
               jlu_cnd[nzja_new] = _ja_lu[j];
               lu_cnd[nzja_new * 2] = _a_lu[j * 2];
               lu_cnd[nzja_new * 2 + 1] = _a_lu[j * 2 + 1];
               nzja_new++;
            }
         }
         listc[i + 1] = nzja_new;
      }

      for (i = 0; i <= _n; i++)
         _ia_lu[i] = listc[i];

      _ja_lu.swap (jlu_cnd);
      _a_lu.swap (lu_cnd);

   }

//
// Perform ILU2 point factorization of the block with future diagonal modification (with structural control)
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::Ilu2BlockIlu2 (int _fcttype, int _fcttype_sch, double _pivmin, double _tau1,
                             double _tau2, double _tau2_sch, double _theta, int _n,
                             int _n_ini, vector < _Int > &_ia_alu,
                             vector < _Int > &_ja_alu, vector < char >&_ja_char_alu,
                             vector < _Flt > &_a_alu, vector < _Int > &_ia_lu,
                             vector < _Int > &_ja_lu, vector < char >&_ja_char_lu,
                             vector < _Flt > &_a_lu, int &_nmodif, double &_eigmin_att,
                             double &_eigmax_att)
   {

// Prepare mask arrays

      vector < _Int > ibegm (_n + 1);
      vector < _Int > madj (_n + 1);
      vector < _Int > iv (_n + 1);
      vector < _Int > imaskc (_n + 1);
      vector < char >imaskchar (_n + 1);
      vector < _Int > listc (_n + 1);

      vector < _Flt > fmaskc (2 * _n + 1);
      vector < _Flt > a_dia (_n + 1);

      _Int *plistc = &listc[0];

      int i;

      for (i = 0; i <= _n; i++)
         ibegm[i] = -1;
      for (i = 0; i <= _n; i++)
         madj[i] = -1;
      for (i = 0; i <= _n; i++)
         imaskc[i] = -1;
      for (i = 0; i < _n; i++)
         a_dia[i] = 0.;

      double eigmin_att = FLT_MAX;
      double eigmax_att = -FLT_MAX;

      int j, jj, k, kcolmn, irwprv, nlistcnew, ind, ind1, irow, irwpr1, jcolmn;
      double auxL, auxU;
      double aux1, aux2, dfnorm;

      int fcttype_loc = _fcttype;
      double tau2_loc = _tau2;

      _Flt fzero = 0.;
      _Flt fone = 1.;

      int icycle_int = -1;

      _ia_lu.resize (_n + 1);

      _ja_lu.resize (0);
      _ja_char_lu.resize (0);
      _a_lu.resize (0);

      _ja_lu.reserve (_n + 1);
      _ja_char_lu.reserve (_n + 1);
      _a_lu.reserve (_n * 2 + 1);

      _ia_lu[0] = 0;

      int nzja_lu = 0;

      _nmodif = 0;

      char jjchar1, jjchar2, jjchar;

      for (i = 0; i < _n; i++) {

         irow = i;

// Init current row

         int nlistcloc = 0;

         icycle_int++;

         if (i < _n_ini) {
            for (k = (int) _ia_alu[i]; k < _ia_alu[i + 1]; k++) {
               kcolmn = (int) _ja_alu[k];
               listc[nlistcloc] = kcolmn;
               nlistcloc++;
               fmaskc[kcolmn * 2] = _a_alu[k * 2];
               fmaskc[kcolmn * 2 + 1] = _a_alu[k * 2 + 1];
               imaskc[kcolmn] = icycle_int;
               imaskchar[kcolmn] = _ja_char_alu[k];
            }
            if (imaskc[i] != icycle_int) {
               listc[nlistcloc] = i;
               nlistcloc++;
               fmaskc[2 * i] = fzero;
               fmaskc[2 * i + 1] = fzero;
               imaskc[i] = icycle_int;
               imaskchar[i] = 0;
            }
         } else {
            listc[nlistcloc] = i;
            nlistcloc++;
            fmaskc[i * 2] = fzero;
            fmaskc[i * 2 + 1] = fzero;
            imaskc[i] = icycle_int;
            imaskchar[i] = 0;
         }

// Update current row

         irwprv = (int) ibegm[irow];

         while (irwprv != -1) {
            irwpr1 = (int) madj[irwprv];
            if (iv[irwprv] < _ia_lu[irwprv + 1]) {

               j = (int) iv[irwprv];

               jj = (int) _ja_lu[j];
               jjchar1 = _ja_char_lu[j] + 1;

               auxL = _a_lu[j * 2];
               auxU = _a_lu[j * 2 + 1];

               jcolmn = jj;
               if (jcolmn < 0)
                  jcolmn = -jcolmn - 1;

               if (jj >= 0) {
                  for (k = (int) iv[irwprv]; k < _ia_lu[irwprv + 1]; k++) {
                     kcolmn = (int) _ja_lu[k];
                     if (kcolmn < 0)
                        kcolmn = -kcolmn - 1;
                     jjchar2 = _ja_char_lu[k] + 1;
                     jjchar = (jjchar1 < jjchar2) ? jjchar2 : jjchar1;
                     if (imaskc[kcolmn] != icycle_int) {
                        listc[nlistcloc] = kcolmn;
                        nlistcloc++;
                        fmaskc[kcolmn * 2] = fzero;
                        fmaskc[kcolmn * 2 + 1] = fzero;
                        imaskc[kcolmn] = icycle_int;
                        imaskchar[kcolmn] = jjchar;
                     } else {
                        if (imaskchar[kcolmn] > jjchar)
                           imaskchar[kcolmn] = jjchar;
                     }
                     fmaskc[kcolmn * 2] -= (_Flt) (auxU * _a_lu[k * 2]);
                     fmaskc[kcolmn * 2 + 1] -= (_Flt) (auxL * _a_lu[k * 2 + 1]);
                  }
               } else {
                  for (k = (int) iv[irwprv]; k < _ia_lu[irwprv + 1]; k++) {
                     kcolmn = (int) _ja_lu[k];
                     if (kcolmn >= 0) {
                        jjchar2 = _ja_char_lu[k] + 1;
                        jjchar = (jjchar1 < jjchar2) ? jjchar2 : jjchar1;
                        if (imaskc[kcolmn] != icycle_int) {
                           listc[nlistcloc] = kcolmn;
                           nlistcloc++;
                           fmaskc[kcolmn * 2] = fzero;
                           fmaskc[kcolmn * 2 + 1] = fzero;
                           imaskc[kcolmn] = icycle_int;
                           imaskchar[kcolmn] = jjchar;
                        } else {
                           if (imaskchar[kcolmn] > jjchar)
                              imaskchar[kcolmn] = jjchar;
                        }
                        fmaskc[kcolmn * 2] -= (_Flt) (auxU * _a_lu[k * 2]);
                        fmaskc[kcolmn * 2 + 1] -= (_Flt) (auxL * _a_lu[k * 2 + 1]);
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
            fcttype_loc = _fcttype;
            tau2_loc = 0.0e0;
         }

         {

            nlistcnew = 0;

            for (j = 0; j < nlistcloc; j++) {
               ind = (int) listc[j];
               if (ind != irow) {
                  jjchar = imaskchar[ind];
                  aux1 = fmaskc[ind * 2];
                  aux2 = fmaskc[ind * 2 + 1];
                  if (aux1 < 0.0)
                     aux1 = -aux1;
                  if (aux2 < 0.0)
                     aux2 = -aux2;
                  dfnorm = (aux1 > aux2) ? aux1 : aux2;
                  if (jjchar > 0 && ((fcttype_loc < 0 && dfnorm < tau2_loc)
                                     || (fcttype_loc >= 0
                                         && (jjchar > fcttype_loc
                                             || dfnorm < tau2_loc)))) {
                     a_dia[irow] += (_Flt) dfnorm;
                     a_dia[ind] += (_Flt) dfnorm;
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

         if (fmaskc[irow * 2] > fzero) {
            fmaskc[irow * 2] += (_Flt) (a_dia[irow] * _theta);
            fmaskc[irow * 2 + 1] += (_Flt) (a_dia[irow] * _theta);
         } else {
            fmaskc[irow * 2] -= (_Flt) (a_dia[irow] * _theta);
            fmaskc[irow * 2 + 1] -= (_Flt) (a_dia[irow] * _theta);
         }

// Factorize current row and split into first and second order

         sort (plistc, plistc + nlistcloc);

         if (i < _n_ini) {
            aux1 = fmaskc[irow * 2];

            if (aux1 > fzero) {
               if (aux1 < eigmin_att)
                  eigmin_att = aux1;
               if (aux1 > eigmax_att)
                  eigmax_att = aux1;
               if (aux1 < _pivmin) {
                  aux1 = _pivmin;
                  _nmodif++;
               }
               aux1 = sqrt (aux1);
               aux1 = fone / aux1;
               auxL = aux1;
               auxU = aux1;
            } else {
               aux1 = -aux1;
               if (aux1 < eigmin_att)
                  eigmin_att = aux1;
               if (aux1 > eigmax_att)
                  eigmax_att = aux1;
               if (aux1 < _pivmin) {
                  aux1 = _pivmin;
                  _nmodif++;
               }
               aux1 = sqrt (aux1);
               aux1 = fone / aux1;
               auxL = -aux1;
               auxU = aux1;
            }
            fmaskc[irow * 2] = (_Flt) (auxL);
            fmaskc[irow * 2 + 1] = (_Flt) (auxU);
         }

         if (i < _n_ini) {
            for (j = 0; j < nlistcloc; j++) {
               ind = (int) listc[j];
               if (ind != irow) {
                  fmaskc[ind * 2] *= (_Flt) (auxU);
                  fmaskc[ind * 2 + 1] *= (_Flt) (auxL);
                  aux1 = fmaskc[ind * 2];
                  aux2 = fmaskc[ind * 2 + 1];
                  if (aux1 < 0.0)
                     aux1 = -aux1;
                  if (aux2 < 0.0)
                     aux2 = -aux2;
                  dfnorm = (aux1 > aux2) ? aux1 : aux2;
                  jjchar = imaskchar[ind];
                  if (jjchar > 0 && dfnorm < _tau1) {
                     listc[j] = (_Int) (-ind - 1);
                  }
               }
            }
         }
// Store computed row elems

         for (j = 0; j < nlistcloc; j++) {
            ind = (int) plistc[j];
            _ja_lu.push_back (ind);
            ind1 = ind;
            if (ind1 < 0)
               ind1 = -ind1 - 1;
            jjchar = imaskchar[ind1];
            _ja_char_lu.push_back (jjchar);
            _a_lu.push_back (fmaskc[ind1 * 2]);
            _a_lu.push_back (fmaskc[ind1 * 2 + 1]);
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
      vector < _Flt > lu_cnd (2 * nzja_new + 1);

      nzja_new = 0;

      listc[0] = 0;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_lu[i]; j < _ia_lu[i + 1]; j++) {
            if (_ja_lu[j] >= 0) {
               jlu_cnd[nzja_new] = _ja_lu[j];
               jlu_char_cnd[nzja_new] = _ja_char_lu[j];
               lu_cnd[nzja_new * 2] = _a_lu[j * 2];
               lu_cnd[nzja_new * 2 + 1] = _a_lu[j * 2 + 1];
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
// Perform symbolic factorization of the block with with matrix degree structural control
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::Ilu2BlockIlu2DegreeSp (ofstream * _pfout_debug, int _fcttype,
                                     int _fcttype_sch, int _n, int _n_ini,
                                     vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                                     vector < char >&_ja_char_alu,
                                     vector < _Int > &_ia_lu, vector < _Int > &_ja_lu,
                                     vector < char >&_ja_char_lu)
   {

// Prepare mask arrays

      vector < _Int > ibegm (_n + 1);
      vector < _Int > madj (_n + 1);
      vector < _Int > iv (_n + 1);
      vector < _Int > imaskc (_n + 1);
      vector < char >imaskchar (_n + 1);
      vector < _Int > listc (_n + 1);

      _Int *plistc = &listc[0];

      int i;

      for (i = 0; i <= _n; i++)
         ibegm[i] = -1;
      for (i = 0; i <= _n; i++)
         madj[i] = -1;
      for (i = 0; i <= _n; i++)
         imaskc[i] = -1;

      int j, jj, k, kcolmn, irwprv, ind, irow, irwpr1, jcolmn;

      int fcttype_loc = _fcttype;

      int icycle_int = -1;

      _ia_lu.resize (_n + 1);

      _ja_lu.resize (0);
      _ja_char_lu.resize (0);

      _ja_lu.reserve (_n + 1);
      _ja_char_lu.reserve (_n + 1);

      _ia_lu[0] = 0;

      int nzja_lu = 0;

      char jjchar1, jjchar2, jjchar;

      for (i = 0; i < _n; i++) {

         irow = i;

//         if (_pfout_debug != NULL) {
//            *_pfout_debug << " Irow = " << irow << endl;
//         }

// Init current row

         int nlistcloc = 0;

         icycle_int++;

         if (i < _n_ini) {
            for (k = (int) _ia_alu[i]; k < _ia_alu[i + 1]; k++) {
               kcolmn = (int) _ja_alu[k];
               listc[nlistcloc] = kcolmn;
               nlistcloc++;
               imaskc[kcolmn] = icycle_int;
               imaskchar[kcolmn] = _ja_char_alu[k];
            }
            if (imaskc[i] != icycle_int) {
               listc[nlistcloc] = i;
               nlistcloc++;
               imaskc[i] = icycle_int;
               imaskchar[i] = 0;
            }
         } else {
            listc[nlistcloc] = i;
            nlistcloc++;
            imaskc[i] = icycle_int;
            imaskchar[i] = 0;
         }

//         if (_pfout_debug != NULL) {
//            *_pfout_debug << " Initial data: " << endl;
//            for (k=0;k<nlistcloc;k++) {
//               kcolmn = (int)listc[k];
//               *_pfout_debug << "    k = " << k << " Icol = " << kcolmn << " Ichar = " << (int)imaskchar[kcolmn] << endl;
//            }
//         }

// Update current row

         irwprv = (int) ibegm[irow];

         while (irwprv != -1) {

            irwpr1 = (int) madj[irwprv];

            if (iv[irwprv] < _ia_lu[irwprv + 1]) {

               j = (int) iv[irwprv];

               jj = (int) _ja_lu[j];
               jjchar1 = _ja_char_lu[j] + 1;

               jcolmn = jj;

//               if (_pfout_debug != NULL) {
//                  *_pfout_debug << "    Update by irwprv = " << irwprv << " jj = " << jj << " jjchar1 = " << (int) jjchar1 << endl;
//               }

               for (k = (int) iv[irwprv]; k < _ia_lu[irwprv + 1]; k++) {
                  kcolmn = (int) _ja_lu[k];
                  jjchar2 = _ja_char_lu[k];
                  jjchar = jjchar1 + jjchar2;
                  if (kcolmn == irow || fcttype_loc < 0
                      || (fcttype_loc >= 0 && jjchar <= fcttype_loc)) {
                     if (imaskc[kcolmn] != icycle_int) {
                        listc[nlistcloc] = kcolmn;
                        nlistcloc++;
                        imaskc[kcolmn] = icycle_int;
                        imaskchar[kcolmn] = jjchar;
//                        if (_pfout_debug != NULL) {
//                           *_pfout_debug << "    Add elem kcolmn = " << kcolmn << " jjchar = " << (int)jjchar << endl;
//                        }
                     } else {
                        if (imaskchar[kcolmn] > jjchar) {
//                           if (_pfout_debug != NULL) {
//                              *_pfout_debug << "    For kcolmn = " << kcolmn << " Reduce jchar from = " << (int)imaskchar[kcolmn] << " to jjchar = " << (int)jjchar << endl;
//                           }
                           imaskchar[kcolmn] = jjchar;
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
               madj[irwprv] = ibegm[jcolmn];
               ibegm[jcolmn] = irwprv;
            }

            iv[irwprv]++;

            irwprv = irwpr1;

         }

// Perform filtering of the data and modify diagonal

         if (i >= _n_ini) {
            fcttype_loc = _fcttype;
         }
// Factorize current row and split into first and second order

         sort (plistc, plistc + nlistcloc);

// Store computed row elems

         int nlistc_flt = 0;

//         if (_pfout_debug != NULL) {
//            *_pfout_debug << "    Final row: " << endl;
//         }

         for (j = 0; j < nlistcloc; j++) {
            ind = (int) plistc[j];
            jjchar = imaskchar[ind];
            if (ind == irow || fcttype_loc < 0
                || (fcttype_loc >= 0 && jjchar <= fcttype_loc)) {
               _ja_lu.push_back (ind);
               _ja_char_lu.push_back (jjchar);
               nzja_lu++;
               plistc[nlistc_flt] = ind;
               nlistc_flt++;
//               *_pfout_debug << "    k = " << nlistc_flt << " Icol = " << ind << " Ichar = " << (int)jjchar << endl;
            }
         }

         nlistcloc = nlistc_flt;

         _ia_lu[i + 1] = nzja_lu;

// Add current row into the transposed structures

         if (nlistcloc > 1) {
            if (i < _n_ini) {
               iv[irow] = _ia_lu[i] + 1;
               ind = (int) listc[1];
               madj[irow] = ibegm[ind];
               ibegm[ind] = (_Int) irow;
            } else {
               iv[irow] = _ia_lu[i + 1] + 1;
            }
         }

      }

// Condense the result (filter second order elems)

      int nzja_new = 0;

      for (i = 0; i < nzja_lu; i++) {
         if (_ja_lu[i] >= 0)
            nzja_new++;
      }

      vector < _Int > jlu_cnd (nzja_new + 1);
      vector < char >jlu_char_cnd (nzja_new + 1);

      nzja_new = 0;

      listc[0] = 0;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_lu[i]; j < _ia_lu[i + 1]; j++) {
            if (_ja_lu[j] >= 0) {
               jlu_cnd[nzja_new] = _ja_lu[j];
               jlu_char_cnd[nzja_new] = _ja_char_lu[j];
               nzja_new++;
            }
         }
         listc[i + 1] = nzja_new;
      }

      for (i = 0; i <= _n; i++)
         _ia_lu[i] = listc[i];

      _ja_lu.swap (jlu_cnd);
      _ja_char_lu.swap (jlu_char_cnd);

//      if (_pfout_debug != NULL) {
//         _Int *p_ia_lu = &_ia_lu[0];
//         _Int *p_ja_lu = &_ja_lu[0];
//         char *p_jachar_lu = &_ja_char_lu[0];
//         *_pfout_debug << "   <<<<<< Final output > > > > > > " << endl;
//         for (i=0;i<_n;i++) {
//            *_pfout_debug << " Irow = " << i << endl;
//            PrintArray (*_pfout_debug, " Cols  ",(int)(p_ia_lu[i+1]-p_ia_lu[i]), p_ja_lu+p_ia_lu[i]);
//            PrintArray (*_pfout_debug, " Chars ",(int)(p_ia_lu[i+1]-p_ia_lu[i]), p_jachar_lu+p_ia_lu[i]);
//         }
//      }

   }

//
// Perform ILU2 point factorization of the block with future diagonal modification (with matrix degree structural control)
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::Ilu2BlockIlu2Degree (ofstream * _pfout_debug, int _fcttype,
                                   int _fcttype_sch, double _pivmin, double _tau1,
                                   double _tau2, double _tau2_sch, double _theta, int _n,
                                   int _n_ini, vector < _Int > &_ia_alu,
                                   vector < _Int > &_ja_alu, vector < char >&_ja_char_alu,
                                   vector < _Flt > &_a_alu, vector < _Int > &_ia_lu,
                                   vector < _Int > &_ja_lu, vector < char >&_ja_char_lu,
                                   vector < _Flt > &_a_lu, int &_nmodif,
                                   double &_eigmin_att, double &_eigmax_att)
   {

// Prepare mask arrays

      int irow_chk = 5;

      vector < _Int > ibegm (_n + 1);
      vector < _Int > madj (_n + 1);
      vector < _Int > iv (_n + 1);
      vector < _Int > imaskc (_n + 1);
      vector < char >imaskchar (_n + 1);
      vector < _Int > listc (_n + 1);

      vector < _Flt > fmaskc (2 * _n + 1);
      vector < _Flt > a_dia (_n + 1);

      _Int *plistc = &listc[0];

      int i;

      for (i = 0; i <= _n; i++)
         ibegm[i] = -1;
      for (i = 0; i <= _n; i++)
         madj[i] = -1;
      for (i = 0; i <= _n; i++)
         imaskc[i] = -1;
      for (i = 0; i < _n; i++)
         a_dia[i] = 0.;

      double eigmin_att = FLT_MAX;
      double eigmax_att = -FLT_MAX;

      int j, jj, k, kcolmn, irwprv, nlistcnew, ind, ind1, irow, irwpr1, jcolmn;
      double auxL, auxU;
      double aux1, aux2, dfnorm;

      int fcttype_loc = _fcttype;
      double tau2_loc = _tau2;

      _Flt fzero = 0.;
      _Flt fone = 1.;

      int icycle_int = -1;

      _ia_lu.resize (_n + 1);

      _ja_lu.resize (0);
      _ja_char_lu.resize (0);
      _a_lu.resize (0);

      _ja_lu.reserve (_n + 1);
      _ja_char_lu.reserve (_n + 1);
      _a_lu.reserve (_n * 2 + 1);

      _ia_lu[0] = 0;

      int nzja_lu = 0;

      _nmodif = 0;

      char jjchar1, jjchar2, jjchar;

      for (i = 0; i < _n; i++) {

//         if (i % irow_chk == 0 || (i >= _n_ini && i < _n_ini+irow_chk)) {
//            if (_pfout_debug != NULL) {
//               *_pfout_debug << " <<< Irow = " << i << endl;
//            }
//         }

         irow = i;

// Init current row

         int nlistcloc = 0;

         icycle_int++;

         if (i < _n_ini) {
            for (k = (int) _ia_alu[i]; k < _ia_alu[i + 1]; k++) {
               kcolmn = (int) _ja_alu[k];
               listc[nlistcloc] = kcolmn;
               nlistcloc++;
               fmaskc[kcolmn * 2] = _a_alu[k * 2];
               fmaskc[kcolmn * 2 + 1] = _a_alu[k * 2 + 1];
               imaskc[kcolmn] = icycle_int;
               imaskchar[kcolmn] = _ja_char_alu[k];
            }
            if (imaskc[i] != icycle_int) {
               listc[nlistcloc] = i;
               nlistcloc++;
               fmaskc[2 * i] = fzero;
               fmaskc[2 * i + 1] = fzero;
               imaskc[i] = icycle_int;
               imaskchar[i] = 0;
            }
         } else {
            listc[nlistcloc] = i;
            nlistcloc++;
            fmaskc[i * 2] = fzero;
            fmaskc[i * 2 + 1] = fzero;
            imaskc[i] = icycle_int;
            imaskchar[i] = 0;
         }

// Update current row

         irwprv = (int) ibegm[irow];

         while (irwprv != -1) {

//            if (i >= _n_ini && i < _n_ini+irow_chk) {
//               if (_pfout_debug != NULL) {
//                  *_pfout_debug << "       ++ Update by Irow prev = " << irwprv << endl;
//               }
//            }

            irwpr1 = (int) madj[irwprv];

            if (iv[irwprv] < _ia_lu[irwprv + 1]) {

               j = (int) iv[irwprv];

               jj = (int) _ja_lu[j];
               jjchar1 = _ja_char_lu[j] + 1;

               auxL = _a_lu[j * 2];
               auxU = _a_lu[j * 2 + 1];

//               if (i >= _n_ini && i < _n_ini+irow_chk) {
//                  if (_pfout_debug != NULL) {
//                     *_pfout_debug << "       ++ Base Values jj = " << jj << " jjchar1 = " << (int)jjchar1 << " auxL = " << auxL << " auxU = " << auxU << endl;
//                  }
//               }

               jcolmn = jj;
               if (jcolmn < 0)
                  jcolmn = -jcolmn - 1;

               if (jj >= 0) {
                  for (k = (int) iv[irwprv]; k < _ia_lu[irwprv + 1]; k++) {
                     kcolmn = (int) _ja_lu[k];
                     if (kcolmn < 0)
                        kcolmn = -kcolmn - 1;
                     jjchar2 = _ja_char_lu[k];
                     jjchar = jjchar1 + jjchar2;
//                     if (i >= _n_ini && i < _n_ini+irow_chk) {
//                        if (_pfout_debug != NULL) {
//                           *_pfout_debug << "          -- Row update kcolmn = " << kcolmn << " jjchar = " << (int)jjchar << endl;
//                        }
//                     }
                     if (kcolmn == irow || fcttype_loc < 0
                         || (fcttype_loc >= 0 && jjchar <= fcttype_loc)) {
                        if (imaskc[kcolmn] != icycle_int) {
                           listc[nlistcloc] = kcolmn;
                           nlistcloc++;
                           fmaskc[kcolmn * 2] = fzero;
                           fmaskc[kcolmn * 2 + 1] = fzero;
                           imaskc[kcolmn] = icycle_int;
                           imaskchar[kcolmn] = jjchar;
                        } else {
                           if (imaskchar[kcolmn] > jjchar)
                              imaskchar[kcolmn] = jjchar;
                        }
                        fmaskc[kcolmn * 2] -= (_Flt) (auxU * _a_lu[k * 2]);
                        fmaskc[kcolmn * 2 + 1] -= (_Flt) (auxL * _a_lu[k * 2 + 1]);
//                        if (i >= _n_ini && i < _n_ini+irow_chk) {
//                           if (_pfout_debug != NULL) {
//                              *_pfout_debug << "          == New values: L " << fmaskc[kcolmn * 2] << " U = " << fmaskc[kcolmn * 2 + 1] << endl;
//                           }
//                        }
                     } else {
                        aux1 = auxU * _a_lu[k * 2];
                        aux2 = auxL * _a_lu[k * 2 + 1];
                        if (aux1 < 0.0)
                           aux1 = -aux1;
                        if (aux2 < 0.0)
                           aux2 = -aux2;
                        dfnorm = (aux1 > aux2) ? aux1 : aux2;
                        a_dia[irow] += (_Flt) dfnorm;
                        a_dia[kcolmn] += (_Flt) dfnorm;
//                        if (i >= _n_ini && i < _n_ini+irow_chk) {
//                           if (_pfout_debug != NULL) {
//                              *_pfout_debug << "          :: Modif: irow = " << a_dia[irow]  << " kcolmn = " << a_dia[kcolmn] << endl;
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
                        if (i >= _n_ini && i < _n_ini + irow_chk) {
                           if (_pfout_debug != NULL) {
                              *_pfout_debug << "          -- Row 2 update kcolmn = " <<
                                 kcolmn << " jjchar = " << (int) jjchar << endl;
                           }
                        }
                        if (kcolmn == irow || fcttype_loc < 0
                            || (fcttype_loc >= 0 && jjchar <= fcttype_loc)) {
                           if (imaskc[kcolmn] != icycle_int) {
                              listc[nlistcloc] = kcolmn;
                              nlistcloc++;
                              fmaskc[kcolmn * 2] = fzero;
                              fmaskc[kcolmn * 2 + 1] = fzero;
                              imaskc[kcolmn] = icycle_int;
                              imaskchar[kcolmn] = jjchar;
                           } else {
                              if (imaskchar[kcolmn] > jjchar)
                                 imaskchar[kcolmn] = jjchar;
                           }
                           fmaskc[kcolmn * 2] -= (_Flt) (auxU * _a_lu[k * 2]);
                           fmaskc[kcolmn * 2 + 1] -= (_Flt) (auxL * _a_lu[k * 2 + 1]);
//                           if (i >= _n_ini && i < _n_ini+irow_chk) {
//                              if (_pfout_debug != NULL) {
//                                 *_pfout_debug << "          == New values 2: L " << fmaskc[kcolmn * 2] << " U = " << fmaskc[kcolmn * 2 + 1] << endl;
//                              }
//                           }
                        } else {
                           aux1 = auxU * _a_lu[k * 2];
                           aux2 = auxL * _a_lu[k * 2 + 1];
                           if (aux1 < 0.0)
                              aux1 = -aux1;
                           if (aux2 < 0.0)
                              aux2 = -aux2;
                           dfnorm = (aux1 > aux2) ? aux1 : aux2;
                           a_dia[irow] += (_Flt) dfnorm;
                           a_dia[kcolmn] += (_Flt) dfnorm;
//                           if (i >= _n_ini && i < _n_ini+irow_chk) {
//                              if (_pfout_debug != NULL) {
//                                 *_pfout_debug << "          :: Modif 2: irow = " << a_dia[irow]  << " kcolmn = " << a_dia[kcolmn] << endl;
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
            fcttype_loc = _fcttype;
            tau2_loc = 0.0e0;
         }

         {

            nlistcnew = 0;

            for (j = 0; j < nlistcloc; j++) {
               ind = (int) listc[j];
               if (ind != irow) {
                  jjchar = imaskchar[ind];
                  aux1 = fmaskc[ind * 2];
                  aux2 = fmaskc[ind * 2 + 1];
                  if (aux1 < 0.0)
                     aux1 = -aux1;
                  if (aux2 < 0.0)
                     aux2 = -aux2;
                  dfnorm = (aux1 > aux2) ? aux1 : aux2;
                  if (jjchar > 0 && ((fcttype_loc < 0 && dfnorm < tau2_loc)
                                     || (fcttype_loc >= 0
                                         && (jjchar > fcttype_loc || dfnorm < tau2_loc)))
                     ) {
                     a_dia[irow] += (_Flt) dfnorm;
                     a_dia[ind] += (_Flt) dfnorm;
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

//         if (i >= _n_ini && i < _n_ini+irow_chk) {
//            if (_pfout_debug != NULL) {
//               *_pfout_debug << "     Dia before modif add: L " << fmaskc[irow * 2] << " U = " << fmaskc[irow * 2 + 1] << endl;
//            }
//         }

         if (fmaskc[irow * 2] > fzero) {
            fmaskc[irow * 2] += (_Flt) (a_dia[irow] * _theta);
            fmaskc[irow * 2 + 1] += (_Flt) (a_dia[irow] * _theta);
         } else {
            fmaskc[irow * 2] -= (_Flt) (a_dia[irow] * _theta);
            fmaskc[irow * 2 + 1] -= (_Flt) (a_dia[irow] * _theta);
         }

//         if (i >= _n_ini && i < _n_ini+irow_chk) {
//            if (_pfout_debug != NULL) {
//               *_pfout_debug << "     Dia after modif add: L " << fmaskc[irow * 2] << " U = " << fmaskc[irow * 2 + 1] << endl;
//            }
//         }

// Factorize current row and split into first and second order

         sort (plistc, plistc + nlistcloc);

//         if (i >= _n_ini && i < _n_ini+irow_chk) {
//            if (_pfout_debug != NULL) {
//               PrintArray (*_pfout_debug," Filtered List",nlistcloc,plistc);
//            }
//         }

         if (i < _n_ini) {
            aux1 = fmaskc[irow * 2];

            if (aux1 > fzero) {
               if (aux1 < eigmin_att)
                  eigmin_att = aux1;
               if (aux1 > eigmax_att)
                  eigmax_att = aux1;
               if (aux1 < _pivmin) {
                  aux1 = _pivmin;
                  _nmodif++;
               }
               aux1 = sqrt (aux1);
               aux1 = fone / aux1;
               auxL = aux1;
               auxU = aux1;
            } else {
               aux1 = -aux1;
               if (aux1 < eigmin_att)
                  eigmin_att = aux1;
               if (aux1 > eigmax_att)
                  eigmax_att = aux1;
               if (aux1 < _pivmin) {
                  aux1 = _pivmin;
                  _nmodif++;
               }
               aux1 = sqrt (aux1);
               aux1 = fone / aux1;
               auxL = -aux1;
               auxU = aux1;
            }
            fmaskc[irow * 2] = (_Flt) (auxL);
            fmaskc[irow * 2 + 1] = (_Flt) (auxU);
         }

         if (i < _n_ini) {
            for (j = 0; j < nlistcloc; j++) {
               ind = (int) listc[j];
               if (ind != irow) {
                  fmaskc[ind * 2] *= (_Flt) (auxU);
                  fmaskc[ind * 2 + 1] *= (_Flt) (auxL);
                  aux1 = fmaskc[ind * 2];
                  aux2 = fmaskc[ind * 2 + 1];
                  if (aux1 < 0.0)
                     aux1 = -aux1;
                  if (aux2 < 0.0)
                     aux2 = -aux2;
                  dfnorm = (aux1 > aux2) ? aux1 : aux2;
                  jjchar = imaskchar[ind];
                  if (jjchar > 0 && dfnorm < _tau1) {
                     listc[j] = (_Int) (-ind - 1);
                  }
               }
            }
         }
// Store computed row elems

         for (j = 0; j < nlistcloc; j++) {
            ind = (int) plistc[j];
            _ja_lu.push_back (ind);
            ind1 = ind;
            if (ind1 < 0)
               ind1 = -ind1 - 1;
            jjchar = imaskchar[ind1];
            _ja_char_lu.push_back (jjchar);
            _a_lu.push_back (fmaskc[ind1 * 2]);
            _a_lu.push_back (fmaskc[ind1 * 2 + 1]);
         }

//         if (i % irow_chk == 0 || (i >= _n_ini && i < _n_ini+irow_chk)) {
//            if (_pfout_debug != NULL) {
//               _Int *p_ja_lu_temp = &_ja_lu[0];
//               char *p_jachar_lu_temp = &_ja_char_lu[0];
//               _Flt *p_a_lu_temp = &_a_lu[0];
//               PrintArray (*_pfout_debug, " Ja row ",nlistcloc,p_ja_lu_temp+nzja_lu);
//               PrintArray (*_pfout_debug, " JaChar row ",nlistcloc,p_jachar_lu_temp+nzja_lu);
//               PrintArray (*_pfout_debug, " A row ",nlistcloc*2,p_a_lu_temp+nzja_lu*2);
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
      vector < _Flt > lu_cnd (2 * nzja_new + 1);

      nzja_new = 0;

      listc[0] = 0;

      for (i = 0; i < _n; i++) {
         for (j = (int) _ia_lu[i]; j < _ia_lu[i + 1]; j++) {
            if (_ja_lu[j] >= 0) {
               jlu_cnd[nzja_new] = _ja_lu[j];
               jlu_char_cnd[nzja_new] = _ja_char_lu[j];
               lu_cnd[nzja_new * 2] = _a_lu[j * 2];
               lu_cnd[nzja_new * 2 + 1] = _a_lu[j * 2 + 1];
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
// Balance diagonal
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::BalanceDiag (int _n, vector < _Int > &_ia_a, vector < _Int > &_ja_a,
                           vector < _Flt > &_a_a, vector < _Int > &_ia_l,
                           vector < _Int > &_ja_l, vector < _Flt > &_a_l,
                           vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                           vector < _Flt > &_a_u, double &_diacorr_min,
                           double &_diacorr_max)
   {

// Get arrays

      _Int *pia_a = &_ia_a[0];
      _Int *pja_a = &_ja_a[0];
      _Flt *pa_a = &_a_a[0];

      _Int *pia_l = &_ia_l[0];
      _Int *pja_l = &_ja_l[0];
      _Flt *pa_l = &_a_l[0];

      _Int *pia_u = &_ia_u[0];
      _Int *pja_u = &_ja_u[0];
      _Flt *pa_u = &_a_u[0];

// Compute sum array

      vector < _Flt > sum_arr (_n + 1);
      _Flt *psum_arr = &sum_arr[0];

      _Flt fzero = (_Flt) 0.0e0;
      _Flt fone = (_Flt) 1.0e0;

      int i;

      for (i = 0; i < _n; i++)
         psum_arr[i] = fzero;

      int ipL, ipU, jjL, jjU, iendL, iendU;

      for (i = 0; i < _n; i++) {
         iendL = (int) pia_l[i + 1] - 1;
         iendU = (int) pia_u[i + 1] - 1;
         ipL = (int) pia_l[i] + 1;
         ipU = (int) pia_u[i] + 1;
         while (ipL <= iendL && ipU <= iendU) {
            jjL = (int) pja_l[ipL];
            jjU = (int) pja_u[ipU];
            if (jjL == jjU) {
               psum_arr[jjL] += (pa_l[ipL] * pa_u[ipU]);
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

      for (i = 0; i < _n; i++) {
         is_found = false;
         for (j = (int) pia_a[i]; j < pia_a[i + 1]; j++) {
            jj = (int) pja_a[j];
            if (jj == i) {
               is_found = true;
               diagA = pa_a[j];
            }
         }
         if (!is_found) {
            throw
               " CFct_impl<_Int,_Flt>::BalanceDiag: error: diagonal value is not found! ";
         }
         ipL = (int) pia_l[i];
         invDiagL = pa_l[ipL];
         ipU = (int) pia_u[i];
         invDiagU = pa_u[ipU];
         sum = diagA - psum_arr[i];
         lu_inv_prod = invDiagL * invDiagU;
         aux = sum * lu_inv_prod;
         if (aux < _diacorr_min)
            _diacorr_min = aux;
         if (aux > _diacorr_max)
            _diacorr_max = aux;
         if (aux < 0.0e0)
            aux = 1.0e0;
         daux = (double) aux;
         daux = sqrt (daux);
         aux = (_Flt) daux;
         aux = fone / aux;
         pa_l[ipL] *= aux;
         pa_u[ipU] *= aux;
      }

   }

//
// Dense scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CFct_impl < _Int,
      _Flt >::DenseScaling (double _sclmin, int _n, _Flt * _a, _Flt * _sclL, _Flt * _sclU,
                            _Flt * _sclLInv, _Flt * _sclUInv, _Flt * _work,
                            double *_dwork, double &_sclmin_att, double &_sclmax_att,
                            int &_nmodif)
   {

// Work memory

      int n_2 = _n * _n;

      _Flt *pu = _work;
      _Flt *pv = pu + n_2;
      _Flt *psv = pv + n_2;

// Compute Svd

      CVector < _Flt >::ComputeSvd (_n, _a, psv, pu, pv, _dwork);

// Modify singular values if necessary

      _sclmin_att = 1.0e100;
      _sclmax_att = -1.0e100;
      _nmodif = 0;

      int i;
      _Flt aux, aux1;

      _Flt fone = (_Flt) 1.0e0;

      for (i = 0; i < _n; i++) {
         if (psv[i] < _sclmin_att)
            _sclmin_att = psv[i];
         if (psv[i] > _sclmax_att)
            _sclmax_att = psv[i];
         if (psv[i] < _sclmin) {
            psv[i] = (_Flt) _sclmin;
            _nmodif++;
         }
         aux = psv[i];
         aux1 = sqrt (aux);
         psv[i] = aux1;
      }

      int j;

      for (i = 0; i < _n; i++) {
         for (j = 0; j < _n; j++) {
            _sclUInv[i * _n + j] = pv[j * _n + i] * psv[j];
         }
      }
      for (i = 0; i < _n; i++) {
         for (j = 0; j < _n; j++) {
            _sclLInv[i * _n + j] = pu[i * _n + j] * psv[i];
         }
      }

      for (i = 0; i < _n; i++)
         psv[i] = fone / psv[i];

      for (i = 0; i < _n; i++) {
         for (j = 0; j < _n; j++) {
            _sclU[i * _n + j] = pv[i * _n + j] * psv[i];
         }
      }
      for (i = 0; i < _n; i++) {
         for (j = 0; j < _n; j++) {
            _sclL[i * _n + j] = pu[j * _n + i] * psv[j];
         }
      }

   }

/// @brief Add and replace operator
//========================================================================================
   template < typename _Int, typename _Flt > CMatrix < _Int, _Flt > &CMatrix < _Int,
      _Flt >::operator+= (const CMatrix < _Int, _Flt > &_aa)
   {

      CMatrix < _Int, _Flt > sum;

      sum.AddBlocks ('+', *this, _aa);

      this->ReplaceFree (sum);

      return *this;

   }

/// @brief Add and replace operator for pairs
//========================================================================================
   template < typename _Int, typename _Flt > CMatrix < _Int, _Flt > &CMatrix < _Int,
      _Flt >::operator%= (const CMatrix < _Int, _Flt > &_aa) {

      CMatrix < _Int, _Flt > sum;

      sum.AddBlocksPairs ('+', *this, _aa);

      this->ReplaceFree (sum);

      return *this;

   }

/// @brief Add sparsities
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::AddBlocksSp (const CMatrix < _Int, _Flt > &_aa, const CMatrix < _Int,
                           _Flt > &_bb)
   {

      int nlistloc_aa = _aa.GetNlist ();
      int nlist2loc_aa = _aa.GetNlist2 ();
      int nzjaloc_aa = _aa.GetNzja ();
      int nzja2loc_aa = _aa.GetNzja2 ();
      int nzjacharloc_aa = _aa.GetNzjaChar ();

      const _Int *plist_aa = _aa.GetListArr ();
      const _Int *plist2_aa = _aa.GetList2Arr ();
      const _Int *pia_aa = _aa.GetIaArr ();
      const _Int *pja_aa = _aa.GetJaArr ();
      const _Int *pja2_aa = _aa.GetJa2Arr ();
      const char *pja_char_aa = _aa.GetJaCharArr ();

      int nlistloc_bb = _bb.GetNlist ();
      int nlist2loc_bb = _bb.GetNlist2 ();
      int nzjaloc_bb = _bb.GetNzja ();
      int nzja2loc_bb = _bb.GetNzja2 ();
      int nzjacharloc_bb = _bb.GetNzjaChar ();

      const _Int *plist_bb = _bb.GetListArr ();
      const _Int *plist2_bb = _bb.GetList2Arr ();
      const _Int *pia_bb = _bb.GetIaArr ();
      const _Int *pja_bb = _bb.GetJaArr ();
      const _Int *pja2_bb = _bb.GetJa2Arr ();
      const char *pja_char_bb = _bb.GetJaCharArr ();

// Separate and check the cases of 1 and 2 indices data

      int icase = 0;

      if (nlist2loc_aa == 0 && nlist2loc_bb == 0 && nzja2loc_aa == 0 && nzja2loc_bb == 0) {

         icase = 1;

      } else if (nlist2loc_aa == nlistloc_aa && nlist2loc_bb == nlistloc_bb
                 && nzja2loc_aa == nzjaloc_aa && nzja2loc_bb == nzjaloc_bb) {

         icase = 2;

      } else {
         throw
            " void CMatrix<_Int,_Flt>::AddBlocksSp: error 1: incompartible block sparsity formats on entry ";
      }

// Separate and check the cases of char indices data

      bool b_is_char = false;

      if (nzjacharloc_aa == 0 && nzjacharloc_bb == 0) {
         b_is_char = false;
      } else if (nzjaloc_aa == nzjacharloc_aa && nzjaloc_bb == nzjacharloc_bb) {
         b_is_char = true;
      } else {
         throw
            " void CMatrix<_Int,_Flt>::AddBlocksSp: error 2: incompartible block sparsity formats on entry ";
      }

// 1 indices data

      if (icase == 1) {

         int nlistmax = nlistloc_aa + nlistloc_bb;

         vector < _Int > list_add (nlistmax + 1);
         _Int *plist_add = &list_add[0];

         _Int isup1, isup2;

         int nlistloc = 0;

         int iend1 = nlistloc_aa - 1;
         int iend2 = nlistloc_bb - 1;

         int ip1 = 0;
         int ip2 = 0;

         while (ip1 <= iend1 || ip2 <= iend2) {
            if (ip1 <= iend1 && ip2 <= iend2) {
               isup1 = plist_aa[ip1];
               isup2 = plist_bb[ip2];
               if (isup1 == isup2) {
                  plist_add[nlistloc] = isup1;
                  ip1++;
                  ip2++;
               } else if (isup1 < isup2) {
                  plist_add[nlistloc] = isup1;
                  ip1++;
               } else {
                  plist_add[nlistloc] = isup2;
                  ip2++;
               }
            } else if (ip1 <= iend1) {
               isup1 = plist_aa[ip1];
               plist_add[nlistloc] = isup1;
               ip1++;
            } else {
               isup2 = plist_bb[ip2];
               plist_add[nlistloc] = isup2;
               ip2++;
            }
            nlistloc++;
         }

         vector < _Int > ia_add (nlistloc + 1);
         _Int *pia_add = &ia_add[0];

// Fill ja array

         int nzjamax = nzjaloc_aa + nzjaloc_bb;

         vector < _Int > ja_add (nzjamax + 1);
         vector < char >ja_char_add (1);

         if (b_is_char) {
            ja_char_add.resize (nzjamax + 1);
         }

         _Int *pja_add = &ja_add[0];
         char *pja_char_add = &ja_char_add[0];

         int ilist;
         int nz, iblk1, iblk2;
         _Int irow, irow1, irow2, jj1, jj2, jp1, jp2, jend1, jend2;
         char jj1_char, jj2_char, jj_char;

         jj1_char = 0;
         jj2_char = 0;
         jj_char = 0;

         pia_add[0] = 0;
         nz = 0;

         ip1 = 0;
         ip2 = 0;

         for (ilist = 0; ilist < nlistloc; ilist++) {
            irow = plist_add[ilist];
            iblk1 = 0;
            iblk2 = 0;
            if (ip1 <= iend1) {
               irow1 = plist_aa[ip1];
               if (irow1 != irow)
                  iblk1 = -1;
            } else {
               iblk1 = -1;
            }
            if (ip2 <= iend2) {
               irow2 = plist_bb[ip2];
               if (irow2 != irow)
                  iblk2 = -1;
            } else {
               iblk2 = -1;
            }
            if (iblk1 >= 0 && iblk2 >= 0) {
               jend1 = pia_aa[ip1 + 1] - 1;
               jend2 = pia_bb[ip2 + 1] - 1;
               jp1 = pia_aa[ip1];
               jp2 = pia_bb[ip2];
               while (jp1 <= jend1 || jp2 <= jend2) {
                  if (jp1 <= jend1 && jp2 <= jend2) {
                     jj1 = pja_aa[jp1];
                     jj2 = pja_bb[jp2];
                     if (b_is_char) {
                        jj1_char = pja_char_aa[jp1];
                        jj2_char = pja_char_bb[jp2];
                     }
                     if (jj1 == jj2) {
                        pja_add[nz] = jj1;
                        if (b_is_char) {
                           jj_char = jj1_char;
                           if (jj2_char < jj_char)
                              jj_char = jj2_char;
                           pja_char_add[nz] = jj_char;
                        }
                        nz++;
                        jp1++;
                        jp2++;
                     } else if (jj1 < jj2) {
                        pja_add[nz] = jj1;
                        if (b_is_char) {
                           pja_char_add[nz] = jj1_char;
                        }
                        nz++;
                        jp1++;
                     } else if (jj1 > jj2) {
                        pja_add[nz] = jj2;
                        if (b_is_char) {
                           pja_char_add[nz] = jj2_char;
                        }
                        nz++;
                        jp2++;
                     }
                  } else if (jp1 <= jend1) {
                     pja_add[nz] = pja_aa[jp1];
                     if (b_is_char) {
                        pja_char_add[nz] = pja_char_aa[jp1];
                     }
                     nz++;
                     jp1++;
                  } else if (jp2 <= jend2) {
                     pja_add[nz] = pja_bb[jp2];
                     if (b_is_char) {
                        pja_char_add[nz] = pja_char_bb[jp2];
                     }
                     nz++;
                     jp2++;
                  }
               }
               ip1++;
               ip2++;
            } else if (iblk1 >= 0) {
               for (jp1 = pia_aa[ip1]; jp1 < pia_aa[ip1 + 1]; jp1++) {
                  pja_add[nz] = pja_aa[jp1];
                  if (b_is_char) {
                     pja_char_add[nz] = pja_char_aa[jp1];
                  }
                  nz++;
               }
               ip1++;
            } else {
               for (jp2 = pia_bb[ip2]; jp2 < pia_bb[ip2 + 1]; jp2++) {
                  pja_add[nz] = pja_bb[jp2];
                  if (b_is_char) {
                     pja_char_add[nz] = pja_char_bb[jp2];
                  }
                  nz++;
               }
               ip2++;
            }
            pia_add[ilist + 1] = (_Int) nz;
         }

         this->ResizeAndSetAllSp (nlistloc, 0, nz, 0);

         _Int *plistnew = this->GetListArr ();
         _Int *pianew = this->GetIaArr ();
         _Int *pjanew = this->GetJaArr ();

         int i;

         for (i = 0; i < nlistloc; i++)
            plistnew[i] = plist_add[i];
         for (i = 0; i <= nlistloc; i++)
            pianew[i] = pia_add[i];
         for (i = 0; i < nz; i++)
            pjanew[i] = pja_add[i];

         if (b_is_char) {
            this->SetNzjaChar (nz);
            this->ResizeJaChar (nz);
            char *pjacharnew = this->GetJaCharArr ();
            for (i = 0; i < nz; i++)
               pjacharnew[i] = pja_char_add[i];
         }
// 2 indices data

      } else if (icase == 2) {

         int nlistmax = nlistloc_aa + nlistloc_bb;

         vector < _Int > list_add (nlistmax + 1);
         vector < _Int > list2_add (nlistmax + 1);

         _Int *plist_add = &list_add[0];
         _Int *plist2_add = &list2_add[0];

         _Int isup1, isup2, iisup1, iisup2;

         int nlistloc = 0;

         int iend1 = nlistloc_aa - 1;
         int iend2 = nlistloc_bb - 1;

         int ip1 = 0;
         int ip2 = 0;

         while (ip1 <= iend1 || ip2 <= iend2) {
            if (ip1 <= iend1 && ip2 <= iend2) {
               isup1 = plist_aa[ip1];
               iisup1 = plist2_aa[ip1];
               isup2 = plist_bb[ip2];
               iisup2 = plist2_bb[ip2];
               if (iisup1 == iisup2) {
                  if (isup1 == isup2) {
                     plist_add[nlistloc] = isup1;
                     plist2_add[nlistloc] = iisup1;
                     ip1++;
                     ip2++;
                  } else if (isup1 < isup2) {
                     plist_add[nlistloc] = isup1;
                     plist2_add[nlistloc] = iisup1;
                     ip1++;
                  } else {
                     plist_add[nlistloc] = isup2;
                     plist2_add[nlistloc] = iisup2;
                     ip2++;
                  }
               } else if (iisup1 < iisup2) {
                  plist_add[nlistloc] = isup1;
                  plist2_add[nlistloc] = iisup1;
                  ip1++;
               } else {
                  plist_add[nlistloc] = isup2;
                  plist2_add[nlistloc] = iisup2;
                  ip2++;
               }
            } else if (ip1 <= iend1) {
               isup1 = plist_aa[ip1];
               iisup1 = plist2_aa[ip1];
               plist_add[nlistloc] = isup1;
               plist2_add[nlistloc] = iisup1;
               ip1++;
            } else {
               isup2 = plist_bb[ip2];
               iisup2 = plist2_bb[ip2];
               plist_add[nlistloc] = isup2;
               plist2_add[nlistloc] = iisup2;
               ip2++;
            }
            nlistloc++;
         }

         vector < _Int > ia_add (nlistloc + 1);

         _Int *pia_add = &ia_add[0];

// Fill ja and ja2 arrays

         int nzjamax = nzjaloc_aa + nzjaloc_bb;

         vector < _Int > ja_add (nzjamax + 1);
         vector < _Int > ja2_add (nzjamax + 1);
         vector < char >ja_char_add (1);
         if (b_is_char) {
            ja_char_add.resize (nzjamax + 1);
         }

         _Int *pja_add = &ja_add[0];
         _Int *pja2_add = &ja2_add[0];
         char *pja_char_add = &ja_char_add[0];

         int ilist;
         _Int jp1, jp2, jend1, jend2;
         int nz, iblk1, iblk2;
         _Int irow, irow1, irow2, jj1, jj2;
         _Int iirow, iirow1, iirow2, jjj1, jjj2;
         char jj1_char, jj2_char, jj_char;

         jj1_char = 0;
         jj2_char = 0;
         jj_char = 0;

         pia_add[0] = 0;
         nz = 0;

         ip1 = 0;
         ip2 = 0;

         for (ilist = 0; ilist < nlistloc; ilist++) {
            irow = plist_add[ilist];
            iirow = plist2_add[ilist];
            iblk1 = 0;
            iblk2 = 0;
            if (ip1 <= iend1) {
               irow1 = plist_aa[ip1];
               iirow1 = plist2_aa[ip1];
               if (irow1 != irow || iirow1 != iirow)
                  iblk1 = -1;
            } else {
               iblk1 = -1;
            }
            if (ip2 <= iend2) {
               irow2 = plist_bb[ip2];
               iirow2 = plist2_bb[ip2];
               if (irow2 != irow || iirow2 != iirow)
                  iblk2 = -1;
            } else {
               iblk2 = -1;
            }
            if (iblk1 >= 0 && iblk2 >= 0) {
               jend1 = pia_aa[ip1 + 1] - 1;
               jend2 = pia_bb[ip2 + 1] - 1;
               jp1 = pia_aa[ip1];
               jp2 = pia_bb[ip2];
               while (jp1 <= jend1 || jp2 <= jend2) {
                  if (jp1 <= jend1 && jp2 <= jend2) {
                     jj1 = pja_aa[jp1];
                     jjj1 = pja2_aa[jp1];
                     jj2 = pja_bb[jp2];
                     jjj2 = pja2_bb[jp2];
                     if (b_is_char) {
                        jj1_char = pja_char_aa[jp1];
                        jj2_char = pja_char_bb[jp2];
                     }
                     if (jjj1 == jjj2) {
                        if (jj1 == jj2) {
                           pja_add[nz] = jj1;
                           pja2_add[nz] = jjj1;
                           if (b_is_char) {
                              jj_char = jj1_char;
                              if (jj2_char < jj_char)
                                 jj_char = jj2_char;
                              pja_char_add[nz] = jj_char;
                           }
                           nz++;
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           pja_add[nz] = jj1;
                           pja2_add[nz] = jjj1;
                           if (b_is_char) {
                              pja_char_add[nz] = jj1_char;
                           }
                           nz++;
                           jp1++;
                        } else if (jj1 > jj2) {
                           pja_add[nz] = jj2;
                           pja2_add[nz] = jjj2;
                           if (b_is_char) {
                              pja_char_add[nz] = jj2_char;
                           }
                           nz++;
                           jp2++;
                        }
                     } else if (jjj1 < jjj2) {
                        pja_add[nz] = jj1;
                        pja2_add[nz] = jjj1;
                        if (b_is_char) {
                           pja_char_add[nz] = jj1_char;
                        }
                        nz++;
                        jp1++;
                     } else {
                        pja_add[nz] = jj2;
                        pja2_add[nz] = jjj2;
                        if (b_is_char) {
                           pja_char_add[nz] = jj2_char;
                        }
                        nz++;
                        jp2++;
                     }
                  } else if (jp1 <= jend1) {
                     pja_add[nz] = pja_aa[jp1];
                     pja2_add[nz] = pja2_aa[jp1];
                     if (b_is_char) {
                        pja_char_add[nz] = pja_char_aa[jp1];
                     }
                     nz++;
                     jp1++;
                  } else if (jp2 <= jend2) {
                     pja_add[nz] = pja_bb[jp2];
                     pja2_add[nz] = pja2_bb[jp2];
                     if (b_is_char) {
                        pja_char_add[nz] = pja_char_bb[jp2];
                     }
                     nz++;
                     jp2++;
                  }
               }
               ip1++;
               ip2++;
            } else if (iblk1 >= 0) {
               for (jp1 = pia_aa[ip1]; jp1 < pia_aa[ip1 + 1]; jp1++) {
                  pja_add[nz] = pja_aa[jp1];
                  pja2_add[nz] = pja2_aa[jp1];
                  if (b_is_char) {
                     pja_char_add[nz] = pja_char_aa[jp1];
                  }
                  nz++;
               }
               ip1++;
            } else {
               for (jp2 = pia_bb[ip2]; jp2 < pia_bb[ip2 + 1]; jp2++) {
                  pja_add[nz] = pja_bb[jp2];
                  pja2_add[nz] = pja2_bb[jp2];
                  if (b_is_char) {
                     pja_char_add[nz] = pja_char_bb[jp2];
                  }
                  nz++;
               }
               ip2++;
            }
            pia_add[ilist + 1] = (_Int) nz;
         }

         this->ResizeAndSetAllSp (nlistloc, nlistloc, nz, nz);

         _Int *plistnew = this->GetListArr ();
         _Int *plist2new = this->GetList2Arr ();
         _Int *pianew = this->GetIaArr ();
         _Int *pjanew = this->GetJaArr ();
         _Int *pja2new = this->GetJa2Arr ();

         int i;

         for (i = 0; i < nlistloc; i++)
            plistnew[i] = plist_add[i];
         for (i = 0; i < nlistloc; i++)
            plist2new[i] = plist2_add[i];
         for (i = 0; i <= nlistloc; i++)
            pianew[i] = pia_add[i];
         for (i = 0; i < nz; i++)
            pjanew[i] = pja_add[i];
         for (i = 0; i < nz; i++)
            pja2new[i] = pja2_add[i];

         if (b_is_char) {
            this->SetNzjaChar (nz);
            this->ResizeJaChar (nz);
            char *pjacharnew = this->GetJaCharArr ();
            for (i = 0; i < nz; i++)
               pjacharnew[i] = pja_char_add[i];
         }

      }

   }

/// @brief Add, subtruct or replace blocks
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::AddBlocks (char _oper, const CMatrix < _Int, _Flt > &_aa,
                         const CMatrix < _Int, _Flt > &_bb)
   {

      int blksize_aa = _aa.GetBSize ();
      int nlistloc_aa = _aa.GetNlist ();
      int nlist2loc_aa = _aa.GetNlist2 ();
      int nzjaloc_aa = _aa.GetNzja ();
      int nzja2loc_aa = _aa.GetNzja2 ();
      int nzaloc_aa = _aa.GetNza ();

      int blksize_bb = _bb.GetBSize ();
      int nlistloc_bb = _bb.GetNlist ();
      int nlist2loc_bb = _bb.GetNlist2 ();
      int nzjaloc_bb = _bb.GetNzja ();
      int nzja2loc_bb = _bb.GetNzja2 ();
      int nzaloc_bb = _bb.GetNza ();

// Separate and check the cases of 1 and 2 indices data

      if (blksize_aa != blksize_bb) {
         throw " CMatrix<_Int,_Flt>::AddBlocks: incompartible small block sizes ";
      }

      if (nlist2loc_aa == 0 && nlist2loc_bb == 0 && nzja2loc_aa == 0 && nzja2loc_bb == 0) {
      } else if (nlist2loc_aa == nlistloc_aa && nlist2loc_bb == nlistloc_bb
                 && nzja2loc_aa == nzjaloc_aa && nzja2loc_bb == nzjaloc_bb) {
      } else {
         throw
            " CMatrix<_Int,_Flt>::AddBlocks: incompartible block sparsity formats on entry ";
      }

      if (nzaloc_aa != nzjaloc_aa || nzaloc_bb != nzjaloc_bb) {
         throw " CMatrix<_Int,_Flt>::AddBlocks: matrix is not in point format ";
      }
// Add sparsities first

      this->AddBlocksSp (_aa, _bb);

      int nzjaloc_add = this->GetNzja ();

// Allocate and init by zeroes float data

      this->ResizeA (nzjaloc_add);
      this->SetNza (nzjaloc_add);

      _Flt *pa_add = this->GetAArr ();

      int i;

      _Flt fzero = (_Flt) 0.0e0;

      for (i = 0; i < nzjaloc_add; i++)
         pa_add[i] = fzero;

// Add float data

      this->AddValues ('+', _aa);
      this->AddValues (_oper, _bb);

      this->b_size = blksize_aa;
   }

/// @brief Add, subtruct or replace for data
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::AddBlocksBxB (char _oper, int _blksize, const CMatrix < _Int, _Flt > &_aa,
                            const CMatrix < _Int, _Flt > &_bb)
   {

      int b_2 = _blksize * _blksize;

      int nlistloc_aa = _aa.GetNlist ();
      int nlist2loc_aa = _aa.GetNlist2 ();
      int nzjaloc_aa = _aa.GetNzja ();
      int nzja2loc_aa = _aa.GetNzja2 ();
      int nzaloc_aa = _aa.GetNza ();

      int nlistloc_bb = _bb.GetNlist ();
      int nlist2loc_bb = _bb.GetNlist2 ();
      int nzjaloc_bb = _bb.GetNzja ();
      int nzja2loc_bb = _bb.GetNzja2 ();
      int nzaloc_bb = _bb.GetNza ();

// Separate and check the cases of 1 and 2 indices data

      if (nlist2loc_aa == 0 && nlist2loc_bb == 0 && nzja2loc_aa == 0 && nzja2loc_bb == 0) {
      } else if (nlist2loc_aa == nlistloc_aa && nlist2loc_bb == nlistloc_bb
                 && nzja2loc_aa == nzjaloc_aa && nzja2loc_bb == nzjaloc_bb) {
      } else {
         throw
            " CMatrix<_Int,_Flt>::AddBlocksPairsBxB: incompartible block sparsity formats on entry ";
      }

      if (nzaloc_aa != nzjaloc_aa * b_2 || nzaloc_bb != nzjaloc_bb * b_2) {
         throw " CMatrix<_Int,_Flt>::AddBlocksPairsBxB: matrix is not in block format ";
      }
// Add sparsities first

      this->AddBlocksSp (_aa, _bb);

      int nzjaloc_add = this->GetNzja ();

// Allocate and init by zeroes float data

      this->ResizeA (nzjaloc_add * b_2);
      this->SetNza (nzjaloc_add * b_2);

      _Flt *pa_add = this->GetAArr ();

      int i;

      _Flt fzero = (_Flt) 0.0e0;

      for (i = 0; i < nzjaloc_add * b_2; i++)
         pa_add[i] = fzero;

// Add float data

      this->AddValuesBxB ('+', _blksize, _aa);
      this->AddValuesBxB (_oper, _blksize, _bb);

   }

/// @brief Add, subtruct or replace for pairs data
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::AddBlocksPairs (char _oper, const CMatrix < _Int, _Flt > &_aa,
                              const CMatrix < _Int, _Flt > &_bb)
   {

      int blksize_aa = _aa.GetBSize ();
      int nlistloc_aa = _aa.GetNlist ();
      int nlist2loc_aa = _aa.GetNlist2 ();
      int nzjaloc_aa = _aa.GetNzja ();
      int nzja2loc_aa = _aa.GetNzja2 ();
      int nzaloc_aa = _aa.GetNza ();

      int blksize_bb = _bb.GetBSize ();
      int nlistloc_bb = _bb.GetNlist ();
      int nlist2loc_bb = _bb.GetNlist2 ();
      int nzjaloc_bb = _bb.GetNzja ();
      int nzja2loc_bb = _bb.GetNzja2 ();
      int nzaloc_bb = _bb.GetNza ();

// Separate and check the cases of 1 and 2 indices data

      if (blksize_aa != blksize_bb) {
         throw " CMatrix<_Int,_Flt>::AddBlocksPairs: incompartible small block sizes ";
      }

      if (nlist2loc_aa == 0 && nlist2loc_bb == 0 && nzja2loc_aa == 0 && nzja2loc_bb == 0) {
      } else if (nlist2loc_aa == nlistloc_aa && nlist2loc_bb == nlistloc_bb
                 && nzja2loc_aa == nzjaloc_aa && nzja2loc_bb == nzjaloc_bb) {
      } else {
         throw
            " CMatrix<_Int,_Flt>::AddBlocksPairs: incompartible block sparsity formats on entry ";
      }

      if (nzaloc_aa != 2 * nzjaloc_aa || nzaloc_bb != 2 * nzjaloc_bb) {
         throw " CMatrix<_Int,_Flt>::AddBlocksPairs: matrix is not in pairs format ";
      }
// Add sparsities first

      this->AddBlocksSp (_aa, _bb);

      int nzjaloc_add = this->GetNzja ();

// Allocate and init by zeroes float data

      this->ResizeA (nzjaloc_add * 2);
      this->SetNza (nzjaloc_add * 2);

      _Flt *pa_add = this->GetAArr ();

      int i;

      _Flt fzero = (_Flt) 0.0e0;

      for (i = 0; i < 2 * nzjaloc_add; i++)
         pa_add[i] = fzero;

// Add float data

      this->AddValuesPairs ('+', _aa);
      this->AddValuesPairs (_oper, _bb);

      this->b_size = blksize_aa;

   }

/// @brief Add, subtruct or replace for pairs data
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::AddBlocksPairsBxB (char _oper, int _blksize, const CMatrix < _Int,
                                 _Flt > &_aa, const CMatrix < _Int, _Flt > &_bb)
   {

      int b_2 = _blksize * _blksize;

      int nlistloc_aa = _aa.GetNlist ();
      int nlist2loc_aa = _aa.GetNlist2 ();
      int nzjaloc_aa = _aa.GetNzja ();
      int nzja2loc_aa = _aa.GetNzja2 ();
      int nzaloc_aa = _aa.GetNza ();

      int nlistloc_bb = _bb.GetNlist ();
      int nlist2loc_bb = _bb.GetNlist2 ();
      int nzjaloc_bb = _bb.GetNzja ();
      int nzja2loc_bb = _bb.GetNzja2 ();
      int nzaloc_bb = _bb.GetNza ();

// Separate and check the cases of 1 and 2 indices data

      if (nlist2loc_aa == 0 && nlist2loc_bb == 0 && nzja2loc_aa == 0 && nzja2loc_bb == 0) {
      } else if (nlist2loc_aa == nlistloc_aa && nlist2loc_bb == nlistloc_bb
                 && nzja2loc_aa == nzjaloc_aa && nzja2loc_bb == nzjaloc_bb) {
      } else {
         throw
            " CMatrix<_Int,_Flt>::AddBlocksPairsBxB: incompartible block sparsity formats on entry ";
      }

      if (nzaloc_aa != 2 * nzjaloc_aa * b_2 || nzaloc_bb != 2 * nzjaloc_bb * b_2) {
         throw
            " CMatrix<_Int,_Flt>::AddBlocksPairsBxB: matrix is not in block pairs format ";
      }
// Add sparsities first

      this->AddBlocksSp (_aa, _bb);

      int nzjaloc_add = this->GetNzja ();

// Allocate and init by zeroes float data

      this->ResizeA (nzjaloc_add * 2 * b_2);
      this->SetNza (nzjaloc_add * 2 * b_2);

      _Flt *pa_add = this->GetAArr ();

      int i;

      _Flt fzero = (_Flt) 0.0e0;

      for (i = 0; i < 2 * nzjaloc_add * b_2; i++)
         pa_add[i] = fzero;

// Add float data

      this->AddValuesPairsBxB ('+', _blksize, _aa);
      this->AddValuesPairsBxB (_oper, _blksize, _bb);

   }

/// @brief Add aa values, sparsity of aa is assumed to be included into the sparsity of current block
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::AddValues (char _oper, const CMatrix < _Int, _Flt > &_aa)
   {

      int nlistloc_aa = _aa.GetNlist ();
      int nlist2loc_aa = _aa.GetNlist2 ();
      int nzjaloc_aa = _aa.GetNzja ();
      int nzja2loc_aa = _aa.GetNzja2 ();
      int nzaloc_aa = _aa.GetNza ();

      const _Int *plist_aa = _aa.GetListArr ();
      const _Int *plist2_aa = _aa.GetList2Arr ();
      const _Int *pia_aa = _aa.GetIaArr ();
      const _Int *pja_aa = _aa.GetJaArr ();
      const _Int *pja2_aa = _aa.GetJa2Arr ();
      const _Flt *pa_aa = _aa.GetAArr ();

      int nlistloc = this->GetNlist ();
      int nlist2loc = this->GetNlist2 ();
      int nzjaloc = this->GetNzja ();
      int nzja2loc = this->GetNzja2 ();
      int nzaloc = this->GetNza ();

      _Int *plist = this->GetListArr ();
      _Int *plist2 = this->GetList2Arr ();
      _Int *pia = this->GetIaArr ();
      _Int *pja = this->GetJaArr ();
      _Int *pja2 = this->GetJa2Arr ();
      _Flt *pa = this->GetAArr ();

// Separate and check the cases of 1 and 2 indices data

      int icase = 0;

      if (nlist2loc_aa == 0 && nlist2loc == 0 && nzja2loc_aa == 0 && nzja2loc == 0) {

         icase = 1;

      } else if (nlist2loc_aa == nlistloc_aa && nlist2loc == nlistloc
                 && nzja2loc_aa == nzjaloc_aa && nzja2loc == nzjaloc) {

         icase = 2;

      } else {
         throw
            " CMatrix<_Int,_Flt>::AddValues: incompartible block sparsity formats on entry ";
      }

      if (nzaloc_aa != nzjaloc_aa || nzaloc != nzjaloc) {
         throw " CMatrix<_Int,_Flt>::AddValues: matrix is not in point format ";
      }
// Scan matrices and add values

// 1 index case

      if (icase == 1) {

         if (_oper == '+') {

            int iend2 = nlistloc_aa - 1;

            int ilist;
            _Int ip2, iblk2, jend1, jend2, jp1, jp2;

            _Int irow, irow2, jj1, jj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  if (irow2 != irow) {
                     iblk2 = -1;
                     if (irow2 < irow)
                        throw
                           " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jj2 = pja_aa[jp2];
                        if (jj1 == jj2) {
                           pa[jp1] += pa_aa[jp2];
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           jp1++;
                        } else if (jj1 > jj2) {
                           throw
                              " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '-') {

            int iend2 = nlistloc_aa - 1;

            int ilist;
            _Int ip2, iblk2, jend1, jend2, jp1, jp2;

            _Int irow, irow2, jj1, jj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  if (irow2 != irow) {
                     iblk2 = -1;
                     if (irow2 < irow)
                        throw
                           " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jj2 = pja_aa[jp2];
                        if (jj1 == jj2) {
                           pa[jp1] -= pa_aa[jp2];
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           jp1++;
                        } else if (jj1 > jj2) {
                           throw
                              " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '=') {

            int iend2 = nlistloc_aa - 1;

            int ilist;
            _Int ip2, iblk2, jend1, jend2, jp1, jp2;

            _Int irow, irow2, jj1, jj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  if (irow2 != irow) {
                     iblk2 = -1;
                     if (irow2 < irow)
                        throw
                           " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jj2 = pja_aa[jp2];
                        if (jj1 == jj2) {
                           pa[jp1] = pa_aa[jp2];
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           jp1++;
                        } else if (jj1 > jj2) {
                           throw
                              " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         }
// 2 index case

      } else {

         if (_oper == '+') {

            int iend2 = nlistloc_aa - 1;

            int ilist;

            _Int ip2, iblk2, jend1, jend2, jp1, jp2;
            _Int irow, iirow, irow2, iirow2, jj1, jjj1, jj2, jjj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iirow = plist2[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  iirow2 = plist2_aa[ip2];
                  if (irow2 != irow || iirow2 != iirow) {
                     iblk2 = -1;
                     if (iirow2 < iirow || (iirow2 == iirow && irow2 < irow)) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                     }
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jjj1 = pja2[jp1];
                        jj2 = pja_aa[jp2];
                        jjj2 = pja2_aa[jp2];
                        if (jjj1 == jjj2) {
                           if (jj1 == jj2) {
                              pa[jp1] += pa_aa[jp2];
                              jp1++;
                              jp2++;
                           } else if (jj1 < jj2) {
                              jp1++;
                           } else if (jj1 > jj2) {
                              throw
                                 " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                           }
                        } else if (jjj1 < jjj2) {
                           jp1++;
                        } else {
                           throw
                              " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '-') {

            int iend2 = nlistloc_aa - 1;

            int ilist;

            _Int ip2, iblk2, jend1, jend2, jp1, jp2;
            _Int irow, iirow, irow2, iirow2, jj1, jjj1, jj2, jjj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iirow = plist2[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  iirow2 = plist2_aa[ip2];
                  if (irow2 != irow || iirow2 != iirow) {
                     iblk2 = -1;
                     if (iirow2 < iirow || (iirow2 == iirow && irow2 < irow)) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                     }
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jjj1 = pja2[jp1];
                        jj2 = pja_aa[jp2];
                        jjj2 = pja2_aa[jp2];
                        if (jjj1 == jjj2) {
                           if (jj1 == jj2) {
                              pa[jp1] -= pa_aa[jp2];
                              jp1++;
                              jp2++;
                           } else if (jj1 < jj2) {
                              jp1++;
                           } else if (jj1 > jj2) {
                              throw
                                 " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                           }
                        } else if (jjj1 < jjj2) {
                           jp1++;
                        } else {
                           throw
                              " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '=') {

            int iend2 = nlistloc_aa - 1;

            int ilist;

            _Int ip2, iblk2, jend1, jend2, jp1, jp2;
            _Int irow, iirow, irow2, iirow2, jj1, jjj1, jj2, jjj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iirow = plist2[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  iirow2 = plist2_aa[ip2];
                  if (irow2 != irow || iirow2 != iirow) {
                     iblk2 = -1;
                     if (iirow2 < iirow || (iirow2 == iirow && irow2 < irow)) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                     }
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jjj1 = pja2[jp1];
                        jj2 = pja_aa[jp2];
                        jjj2 = pja2_aa[jp2];
                        if (jjj1 == jjj2) {
                           if (jj1 == jj2) {
                              pa[jp1] = pa_aa[jp2];
                              jp1++;
                              jp2++;
                           } else if (jj1 < jj2) {
                              jp1++;
                           } else if (jj1 > jj2) {
                              throw
                                 " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                           }
                        } else if (jjj1 < jjj2) {
                           jp1++;
                        } else {
                           throw
                              " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValues: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         }

      }

   }

/// @brief Add aa pairs of values, sparsity of aa is assumed to be included into the sparsity of current block
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::AddValuesBxB (char _oper, int _blksize, const CMatrix < _Int, _Flt > &_aa)
   {

      int b_2 = _blksize * _blksize;

      int nlistloc_aa = _aa.GetNlist ();
      int nlist2loc_aa = _aa.GetNlist2 ();
      int nzjaloc_aa = _aa.GetNzja ();
      int nzja2loc_aa = _aa.GetNzja2 ();
      int nzaloc_aa = _aa.GetNza ();

      const _Int *plist_aa = _aa.GetListArr ();
      const _Int *plist2_aa = _aa.GetList2Arr ();
      const _Int *pia_aa = _aa.GetIaArr ();
      const _Int *pja_aa = _aa.GetJaArr ();
      const _Int *pja2_aa = _aa.GetJa2Arr ();
      const _Flt *pa_aa = _aa.GetAArr ();

      int nlistloc = this->GetNlist ();
      int nlist2loc = this->GetNlist2 ();
      int nzjaloc = this->GetNzja ();
      int nzja2loc = this->GetNzja2 ();
      int nzaloc = this->GetNza ();

      _Int *plist = this->GetListArr ();
      _Int *plist2 = this->GetList2Arr ();
      _Int *pia = this->GetIaArr ();
      _Int *pja = this->GetJaArr ();
      _Int *pja2 = this->GetJa2Arr ();
      _Flt *pa = this->GetAArr ();

// Separate and check the cases of 1 and 2 indices data

      int icase = 0;

      if (nlist2loc_aa == 0 && nlist2loc == 0 && nzja2loc_aa == 0 && nzja2loc == 0) {

         icase = 1;

      } else if (nlist2loc_aa == nlistloc_aa && nlist2loc == nlistloc
                 && nzja2loc_aa == nzjaloc_aa && nzja2loc == nzjaloc) {

         icase = 2;

      } else {
         throw
            " CMatrix<_Int,_Flt>::AddValuesPairsBxB: incompartible block sparsity formats on entry ";
      }

      if (nzaloc_aa != nzjaloc_aa * b_2 || nzaloc != nzjaloc * b_2) {
         throw " CMatrix<_Int,_Flt>::AddValuesPairsBxB: matrix is not in block format ";
      }
// Scan matrices and add values

// 1 index case

      if (icase == 1) {

         if (_oper == '+') {

            int iend2 = nlistloc_aa - 1;

            int ilist;
            _Int ip2, iblk2, jend1, jend2, jp1, jp2;

            _Int irow, irow2, jj1, jj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  if (irow2 != irow) {
                     iblk2 = -1;
                     if (irow2 < irow)
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jj2 = pja_aa[jp2];
                        if (jj1 == jj2) {
                           CVector < _Flt >::AddReplaceVector (b_2, pa_aa + jp2 * b_2,
                                                               pa + jp1 * b_2);
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           jp1++;
                        } else if (jj1 > jj2) {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '-') {

            int iend2 = nlistloc_aa - 1;

            int ilist;
            _Int ip2, iblk2, jend1, jend2, jp1, jp2;

            _Int irow, irow2, jj1, jj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  if (irow2 != irow) {
                     iblk2 = -1;
                     if (irow2 < irow)
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jj2 = pja_aa[jp2];
                        if (jj1 == jj2) {
                           CVector < _Flt >::SubtractReplaceVector (b_2,
                                                                    pa_aa + jp2 * b_2,
                                                                    pa + jp1 * b_2);
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           jp1++;
                        } else if (jj1 > jj2) {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '=') {

            int iend2 = nlistloc_aa - 1;

            int ilist;
            _Int ip2, iblk2, jend1, jend2, jp1, jp2;

            _Int irow, irow2, jj1, jj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  if (irow2 != irow) {
                     iblk2 = -1;
                     if (irow2 < irow)
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jj2 = pja_aa[jp2];
                        if (jj1 == jj2) {
                           CVector < _Flt >::CopyVector (b_2, pa_aa + jp2 * b_2,
                                                         pa + jp1 * b_2);
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           jp1++;
                        } else if (jj1 > jj2) {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         }
// 2 index case

      } else {

         if (_oper == '+') {

            int iend2 = nlistloc_aa - 1;

            int ilist;

            _Int ip2, iblk2, jend1, jend2, jp1, jp2;
            _Int irow, iirow, irow2, iirow2, jj1, jjj1, jj2, jjj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iirow = plist2[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  iirow2 = plist2_aa[ip2];
                  if (irow2 != irow || iirow2 != iirow) {
                     iblk2 = -1;
                     if (iirow2 < iirow || (iirow2 == iirow && irow2 < irow)) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                     }
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jjj1 = pja2[jp1];
                        jj2 = pja_aa[jp2];
                        jjj2 = pja2_aa[jp2];
                        if (jjj1 == jjj2) {
                           if (jj1 == jj2) {
                              CVector < _Flt >::AddReplaceVector (b_2, pa_aa + jp2 * b_2,
                                                                  pa + jp1 * b_2);
                              jp1++;
                              jp2++;
                           } else if (jj1 < jj2) {
                              jp1++;
                           } else if (jj1 > jj2) {
                              throw
                                 " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                           }
                        } else if (jjj1 < jjj2) {
                           jp1++;
                        } else {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '-') {

            int iend2 = nlistloc_aa - 1;

            int ilist;

            _Int ip2, iblk2, jend1, jend2, jp1, jp2;
            _Int irow, iirow, irow2, iirow2, jj1, jjj1, jj2, jjj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iirow = plist2[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  iirow2 = plist2_aa[ip2];
                  if (irow2 != irow || iirow2 != iirow) {
                     iblk2 = -1;
                     if (iirow2 < iirow || (iirow2 == iirow && irow2 < irow)) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                     }
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jjj1 = pja2[jp1];
                        jj2 = pja_aa[jp2];
                        jjj2 = pja2_aa[jp2];
                        if (jjj1 == jjj2) {
                           if (jj1 == jj2) {
                              CVector < _Flt >::SubtractReplaceVector (b_2,
                                                                       pa_aa + jp2 * b_2,
                                                                       pa + jp1 * b_2);
                              jp1++;
                              jp2++;
                           } else if (jj1 < jj2) {
                              jp1++;
                           } else if (jj1 > jj2) {
                              throw
                                 " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                           }
                        } else if (jjj1 < jjj2) {
                           jp1++;
                        } else {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '=') {

            int iend2 = nlistloc_aa - 1;

            int ilist;

            _Int ip2, iblk2, jend1, jend2, jp1, jp2;
            _Int irow, iirow, irow2, iirow2, jj1, jjj1, jj2, jjj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iirow = plist2[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  iirow2 = plist2_aa[ip2];
                  if (irow2 != irow || iirow2 != iirow) {
                     iblk2 = -1;
                     if (iirow2 < iirow || (iirow2 == iirow && irow2 < irow)) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                     }
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jjj1 = pja2[jp1];
                        jj2 = pja_aa[jp2];
                        jjj2 = pja2_aa[jp2];
                        if (jjj1 == jjj2) {
                           if (jj1 == jj2) {
                              CVector < _Flt >::CopyVector (b_2, pa_aa + jp2 * b_2,
                                                            pa + jp1 * b_2);
                              jp1++;
                              jp2++;
                           } else if (jj1 < jj2) {
                              jp1++;
                           } else if (jj1 > jj2) {
                              throw
                                 " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                           }
                        } else if (jjj1 < jjj2) {
                           jp1++;
                        } else {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesBxB: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         }

      }

   }

/// @brief Add aa pairs of values, sparsity of aa is assumed to be included into the sparsity of current block
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::AddValuesPairs (char _oper, const CMatrix < _Int, _Flt > &_aa)
   {

      int nlistloc_aa = _aa.GetNlist ();
      int nlist2loc_aa = _aa.GetNlist2 ();
      int nzjaloc_aa = _aa.GetNzja ();
      int nzja2loc_aa = _aa.GetNzja2 ();
      int nzaloc_aa = _aa.GetNza ();

      const _Int *plist_aa = _aa.GetListArr ();
      const _Int *plist2_aa = _aa.GetList2Arr ();
      const _Int *pia_aa = _aa.GetIaArr ();
      const _Int *pja_aa = _aa.GetJaArr ();
      const _Int *pja2_aa = _aa.GetJa2Arr ();
      const _Flt *pa_aa = _aa.GetAArr ();

      int nlistloc = this->GetNlist ();
      int nlist2loc = this->GetNlist2 ();
      int nzjaloc = this->GetNzja ();
      int nzja2loc = this->GetNzja2 ();
      int nzaloc = this->GetNza ();

      _Int *plist = this->GetListArr ();
      _Int *plist2 = this->GetList2Arr ();
      _Int *pia = this->GetIaArr ();
      _Int *pja = this->GetJaArr ();
      _Int *pja2 = this->GetJa2Arr ();
      _Flt *pa = this->GetAArr ();

// Separate and check the cases of 1 and 2 indices data

      int icase = 0;

      if (nlist2loc_aa == 0 && nlist2loc == 0 && nzja2loc_aa == 0 && nzja2loc == 0) {

         icase = 1;

      } else if (nlist2loc_aa == nlistloc_aa && nlist2loc == nlistloc
                 && nzja2loc_aa == nzjaloc_aa && nzja2loc == nzjaloc) {

         icase = 2;

      } else {
         throw
            " CMatrix<_Int,_Flt>::AddValuesPairs: incompartible block sparsity formats on entry ";
      }

      if (nzaloc_aa != 2 * nzjaloc_aa || nzaloc != 2 * nzjaloc) {
         throw " CMatrix<_Int,_Flt>::AddValuesPairs: matrix is not in pairs format ";
      }
// Scan matrices and add values

// 1 index case

      if (icase == 1) {

         if (_oper == '+') {

            int iend2 = nlistloc_aa - 1;

            int ilist;
            _Int ip2, iblk2, jend1, jend2, jp1, jp2;

            _Int irow, irow2, jj1, jj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  if (irow2 != irow) {
                     iblk2 = -1;
                     if (irow2 < irow)
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jj2 = pja_aa[jp2];
                        if (jj1 == jj2) {
                           pa[jp1 * 2] += pa_aa[jp2 * 2];
                           pa[jp1 * 2 + 1] += pa_aa[jp2 * 2 + 1];
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           jp1++;
                        } else if (jj1 > jj2) {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '-') {

            int iend2 = nlistloc_aa - 1;

            int ilist;
            _Int ip2, iblk2, jend1, jend2, jp1, jp2;

            _Int irow, irow2, jj1, jj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  if (irow2 != irow) {
                     iblk2 = -1;
                     if (irow2 < irow)
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jj2 = pja_aa[jp2];
                        if (jj1 == jj2) {
                           pa[jp1 * 2] -= pa_aa[jp2 * 2];
                           pa[jp1 * 2 + 1] -= pa_aa[jp2 * 2 + 1];
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           jp1++;
                        } else if (jj1 > jj2) {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '=') {

            int iend2 = nlistloc_aa - 1;

            int ilist;
            _Int ip2, iblk2, jend1, jend2, jp1, jp2;

            _Int irow, irow2, jj1, jj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  if (irow2 != irow) {
                     iblk2 = -1;
                     if (irow2 < irow)
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jj2 = pja_aa[jp2];
                        if (jj1 == jj2) {
                           pa[jp1 * 2] = pa_aa[jp2 * 2];
                           pa[jp1 * 2 + 1] = pa_aa[jp2 * 2 + 1];
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           jp1++;
                        } else if (jj1 > jj2) {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         }
// 2 index case

      } else {

         if (_oper == '+') {

            int iend2 = nlistloc_aa - 1;

            int ilist;

            _Int ip2, iblk2, jend1, jend2, jp1, jp2;
            _Int irow, iirow, irow2, iirow2, jj1, jjj1, jj2, jjj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iirow = plist2[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  iirow2 = plist2_aa[ip2];
                  if (irow2 != irow || iirow2 != iirow) {
                     iblk2 = -1;
                     if (iirow2 < iirow || (iirow2 == iirow && irow2 < irow)) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                     }
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jjj1 = pja2[jp1];
                        jj2 = pja_aa[jp2];
                        jjj2 = pja2_aa[jp2];
                        if (jjj1 == jjj2) {
                           if (jj1 == jj2) {
                              pa[jp1 * 2] += pa_aa[jp2 * 2];
                              pa[jp1 * 2 + 1] += pa_aa[jp2 * 2 + 1];
                              jp1++;
                              jp2++;
                           } else if (jj1 < jj2) {
                              jp1++;
                           } else if (jj1 > jj2) {
                              throw
                                 " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                           }
                        } else if (jjj1 < jjj2) {
                           jp1++;
                        } else {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '-') {

            int iend2 = nlistloc_aa - 1;

            int ilist;

            _Int ip2, iblk2, jend1, jend2, jp1, jp2;
            _Int irow, iirow, irow2, iirow2, jj1, jjj1, jj2, jjj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iirow = plist2[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  iirow2 = plist2_aa[ip2];
                  if (irow2 != irow || iirow2 != iirow) {
                     iblk2 = -1;
                     if (iirow2 < iirow || (iirow2 == iirow && irow2 < irow)) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                     }
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jjj1 = pja2[jp1];
                        jj2 = pja_aa[jp2];
                        jjj2 = pja2_aa[jp2];
                        if (jjj1 == jjj2) {
                           if (jj1 == jj2) {
                              pa[jp1 * 2] -= pa_aa[jp2 * 2];
                              pa[jp1 * 2 + 1] -= pa_aa[jp2 * 2 + 1];
                              jp1++;
                              jp2++;
                           } else if (jj1 < jj2) {
                              jp1++;
                           } else if (jj1 > jj2) {
                              throw
                                 " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                           }
                        } else if (jjj1 < jjj2) {
                           jp1++;
                        } else {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '=') {

            int iend2 = nlistloc_aa - 1;

            int ilist;

            _Int ip2, iblk2, jend1, jend2, jp1, jp2;
            _Int irow, iirow, irow2, iirow2, jj1, jjj1, jj2, jjj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iirow = plist2[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  iirow2 = plist2_aa[ip2];
                  if (irow2 != irow || iirow2 != iirow) {
                     iblk2 = -1;
                     if (iirow2 < iirow || (iirow2 == iirow && irow2 < irow)) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                     }
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jjj1 = pja2[jp1];
                        jj2 = pja_aa[jp2];
                        jjj2 = pja2_aa[jp2];
                        if (jjj1 == jjj2) {
                           if (jj1 == jj2) {
                              pa[jp1 * 2] = pa_aa[jp2 * 2];
                              pa[jp1 * 2 + 1] = pa_aa[jp2 * 2 + 1];
                              jp1++;
                              jp2++;
                           } else if (jj1 < jj2) {
                              jp1++;
                           } else if (jj1 > jj2) {
                              throw
                                 " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                           }
                        } else if (jjj1 < jjj2) {
                           jp1++;
                        } else {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairs: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         }

      }

   }

/// @brief Add aa pairs of values, sparsity of aa is assumed to be included into the sparsity of current block
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::AddValuesPairsBxB (char _oper, int _blksize, const CMatrix < _Int,
                                 _Flt > &_aa)
   {

      int b_2 = _blksize * _blksize;
      int b_2_2 = b_2 * 2;

      int nlistloc_aa = _aa.GetNlist ();
      int nlist2loc_aa = _aa.GetNlist2 ();
      int nzjaloc_aa = _aa.GetNzja ();
      int nzja2loc_aa = _aa.GetNzja2 ();
      int nzaloc_aa = _aa.GetNza ();

      const _Int *plist_aa = _aa.GetListArr ();
      const _Int *plist2_aa = _aa.GetList2Arr ();
      const _Int *pia_aa = _aa.GetIaArr ();
      const _Int *pja_aa = _aa.GetJaArr ();
      const _Int *pja2_aa = _aa.GetJa2Arr ();
      const _Flt *pa_aa = _aa.GetAArr ();

      int nlistloc = this->GetNlist ();
      int nlist2loc = this->GetNlist2 ();
      int nzjaloc = this->GetNzja ();
      int nzja2loc = this->GetNzja2 ();
      int nzaloc = this->GetNza ();

      _Int *plist = this->GetListArr ();
      _Int *plist2 = this->GetList2Arr ();
      _Int *pia = this->GetIaArr ();
      _Int *pja = this->GetJaArr ();
      _Int *pja2 = this->GetJa2Arr ();
      _Flt *pa = this->GetAArr ();

// Separate and check the cases of 1 and 2 indices data

      int icase = 0;

      if (nlist2loc_aa == 0 && nlist2loc == 0 && nzja2loc_aa == 0 && nzja2loc == 0) {

         icase = 1;

      } else if (nlist2loc_aa == nlistloc_aa && nlist2loc == nlistloc
                 && nzja2loc_aa == nzjaloc_aa && nzja2loc == nzjaloc) {

         icase = 2;

      } else {
         throw
            " CMatrix<_Int,_Flt>::AddValuesPairsBxB: incompartible block sparsity formats on entry ";
      }

      if (nzaloc_aa != 2 * nzjaloc_aa * b_2 || nzaloc != 2 * nzjaloc * b_2) {
         throw " CMatrix<_Int,_Flt>::AddValuesPairsBxB: matrix is not in pairs format ";
      }
// Scan matrices and add values

// 1 index case

      if (icase == 1) {

         if (_oper == '+') {

            int iend2 = nlistloc_aa - 1;

            int ilist;
            _Int ip2, iblk2, jend1, jend2, jp1, jp2;

            _Int irow, irow2, jj1, jj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  if (irow2 != irow) {
                     iblk2 = -1;
                     if (irow2 < irow)
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jj2 = pja_aa[jp2];
                        if (jj1 == jj2) {
                           CVector < _Flt >::AddReplaceVector (b_2_2, pa_aa + jp2 * b_2_2,
                                                               pa + jp1 * b_2_2);
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           jp1++;
                        } else if (jj1 > jj2) {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '-') {

            int iend2 = nlistloc_aa - 1;

            int ilist;
            _Int ip2, iblk2, jend1, jend2, jp1, jp2;

            _Int irow, irow2, jj1, jj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  if (irow2 != irow) {
                     iblk2 = -1;
                     if (irow2 < irow)
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jj2 = pja_aa[jp2];
                        if (jj1 == jj2) {
                           CVector < _Flt >::SubtractReplaceVector (b_2_2,
                                                                    pa_aa + jp2 * b_2_2,
                                                                    pa + jp1 * b_2_2);
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           jp1++;
                        } else if (jj1 > jj2) {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '=') {

            int iend2 = nlistloc_aa - 1;

            int ilist;
            _Int ip2, iblk2, jend1, jend2, jp1, jp2;

            _Int irow, irow2, jj1, jj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  if (irow2 != irow) {
                     iblk2 = -1;
                     if (irow2 < irow)
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jj2 = pja_aa[jp2];
                        if (jj1 == jj2) {
                           CVector < _Flt >::CopyVector (b_2_2, pa_aa + jp2 * b_2_2,
                                                         pa + jp1 * b_2_2);
                           jp1++;
                           jp2++;
                        } else if (jj1 < jj2) {
                           jp1++;
                        } else if (jj1 > jj2) {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         }
// 2 index case

      } else {

         if (_oper == '+') {

            int iend2 = nlistloc_aa - 1;

            int ilist;

            _Int ip2, iblk2, jend1, jend2, jp1, jp2;
            _Int irow, iirow, irow2, iirow2, jj1, jjj1, jj2, jjj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iirow = plist2[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  iirow2 = plist2_aa[ip2];
                  if (irow2 != irow || iirow2 != iirow) {
                     iblk2 = -1;
                     if (iirow2 < iirow || (iirow2 == iirow && irow2 < irow)) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                     }
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jjj1 = pja2[jp1];
                        jj2 = pja_aa[jp2];
                        jjj2 = pja2_aa[jp2];
                        if (jjj1 == jjj2) {
                           if (jj1 == jj2) {
                              CVector < _Flt >::AddReplaceVector (b_2_2,
                                                                  pa_aa + jp2 * b_2_2,
                                                                  pa + jp1 * b_2_2);
                              jp1++;
                              jp2++;
                           } else if (jj1 < jj2) {
                              jp1++;
                           } else if (jj1 > jj2) {
                              throw
                                 " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                           }
                        } else if (jjj1 < jjj2) {
                           jp1++;
                        } else {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '-') {

            int iend2 = nlistloc_aa - 1;

            int ilist;

            _Int ip2, iblk2, jend1, jend2, jp1, jp2;
            _Int irow, iirow, irow2, iirow2, jj1, jjj1, jj2, jjj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iirow = plist2[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  iirow2 = plist2_aa[ip2];
                  if (irow2 != irow || iirow2 != iirow) {
                     iblk2 = -1;
                     if (iirow2 < iirow || (iirow2 == iirow && irow2 < irow)) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                     }
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jjj1 = pja2[jp1];
                        jj2 = pja_aa[jp2];
                        jjj2 = pja2_aa[jp2];
                        if (jjj1 == jjj2) {
                           if (jj1 == jj2) {
                              CVector < _Flt >::SubtractReplaceVector (b_2_2,
                                                                       pa_aa +
                                                                       jp2 * b_2_2,
                                                                       pa + jp1 * b_2_2);
                              jp1++;
                              jp2++;
                           } else if (jj1 < jj2) {
                              jp1++;
                           } else if (jj1 > jj2) {
                              throw
                                 " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                           }
                        } else if (jjj1 < jjj2) {
                           jp1++;
                        } else {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         } else if (_oper == '=') {

            int iend2 = nlistloc_aa - 1;

            int ilist;

            _Int ip2, iblk2, jend1, jend2, jp1, jp2;
            _Int irow, iirow, irow2, iirow2, jj1, jjj1, jj2, jjj2;

            ip2 = 0;

            for (ilist = 0; ilist < nlistloc; ilist++) {
               irow = plist[ilist];
               iirow = plist2[ilist];
               iblk2 = 0;
               if (ip2 <= iend2) {
                  irow2 = plist_aa[ip2];
                  iirow2 = plist2_aa[ip2];
                  if (irow2 != irow || iirow2 != iirow) {
                     iblk2 = -1;
                     if (iirow2 < iirow || (iirow2 == iirow && irow2 < irow)) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                     }
                  }
               } else {
                  iblk2 = -1;
               }
               if (iblk2 >= 0) {
                  jend1 = pia[ilist + 1] - 1;
                  jend2 = pia_aa[ip2 + 1] - 1;
                  jp1 = pia[ilist];
                  jp2 = pia_aa[ip2];
                  while (jp1 <= jend1 || jp2 <= jend2) {
                     if (jp1 <= jend1 && jp2 <= jend2) {
                        jj1 = pja[jp1];
                        jjj1 = pja2[jp1];
                        jj2 = pja_aa[jp2];
                        jjj2 = pja2_aa[jp2];
                        if (jjj1 == jjj2) {
                           if (jj1 == jj2) {
                              CVector < _Flt >::CopyVector (b_2_2, pa_aa + jp2 * b_2_2,
                                                            pa + jp1 * b_2_2);
                              jp1++;
                              jp2++;
                           } else if (jj1 < jj2) {
                              jp1++;
                           } else if (jj1 > jj2) {
                              throw
                                 " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                           }
                        } else if (jjj1 < jjj2) {
                           jp1++;
                        } else {
                           throw
                              " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                        }
                     } else if (jp1 <= jend1) {
                        jp1++;
                     } else if (jp2 <= jend2) {
                        throw
                           " CMatrix<_Int,_Flt>::AddValuesPairsBxB: sparsity structure is not inclusive ";
                     }
                  }
                  ip2++;
               }
            }

         }

      }

   }

/// @brief Compute column list for transposed matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::ComputeList2 (int &_icycle, int *_imask, int *_listloc)
   {

      int nzjaloc = this->GetNzja ();
      _Int *pja = this->GetJaArr ();

      _icycle++;

      int nlistloc = 0;

      int i, jj;

      for (i = 0; i < nzjaloc; i++) {
         jj = (int) pja[i];
         if (_imask[jj] != _icycle) {
            _listloc[nlistloc] = jj;
            nlistloc++;
            _imask[jj] = _icycle;
         }
      }

      sort (_listloc, _listloc + nlistloc);

      this->ResizeList2 (nlistloc);
      this->SetNlist2 (nlistloc);

      _Int *plist2 = this->GetList2Arr ();

      for (i = 0; i < nlistloc; i++)
         plist2[i] = (_Int) _listloc[i];

   }

/// @brief Get main submatrix
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::GetMainSubmatrixSp (int _ibeg, int _iend, CMatrix < _Int,
                                  _Flt > &_asub) const
   {

      int nlistloc = this->GetNlist ();
      const _Int *pia = this->GetIaArr ();
      const _Int *pja = this->GetJaArr ();

      if (_ibeg >= nlistloc || _iend >= nlistloc)
      {
         cout << " CMatrix<>::GetMainSubmatrixSp: error in index interval ! " << endl;
         throw " CMatrix<>::GetMainSubmatrixSp: error in index interval ! ";
      }

      int niloc = _iend - _ibeg + 1;
      int nzja_sub = (int) (pia[_iend + 1] - pia[_ibeg]);

      _asub.ResizeAndSetAllSp (niloc, 0, nzja_sub, 0);

      _Int *plist_sub = _asub.GetListArr ();
      _Int *pia_sub = _asub.GetIaArr ();
      _Int *pja_sub = _asub.GetJaArr ();

      int i;

      for (i = 0; i < niloc; i++)
         plist_sub[i] = (_Int) (_ibeg + i);

      int nzja_new = 0;
      pia_sub[0] = 0;

      int i1, jj, j;

      for (i = 0; i < niloc; i++) {
         i1 = i + _ibeg;
         for (j = (int) pia[i1]; j < pia[i1 + 1]; j++) {
            jj = (int) pja[j];
            if (jj >= _ibeg && jj <= _iend) {
               pja_sub[nzja_new] = (_Int) jj;
               nzja_new++;
            }
         }
         pia_sub[i + 1] = (_Int) nzja_new;
      }

      _asub.SetNzja (nzja_new);

   }

/// @brief Get last submatrix
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::GetLastSubmatrixBxB (bool _b_is_char, int _size_bxb, int _n_ini,
                                   CMatrix < _Int, _Flt > &_asub_last) const
   {

      int nlistloc = this->GetNlist ();

      const vector < _Int > *pia = this->GetIa ();
      const vector < _Int > *pja = this->GetJa ();
      const vector < char >*pjachar = this->GetJaChar ();
      const vector < _Flt > *pa = this->GetA ();

        vector < _Int > *pia_last = _asub_last.GetIa ();
        vector < _Int > *pja_last = _asub_last.GetJa ();
        vector < char >*pjachar_last = _asub_last.GetJaChar ();
        vector < _Flt > *pa_last = _asub_last.GetA ();

        CFct_impl < _Int, _Flt >::GetLastSubmatrixBxB (_b_is_char, nlistloc, _size_bxb,
                                                       *pia, *pja, *pjachar, *pa, _n_ini,
                                                       *pia_last, *pja_last,
                                                       *pjachar_last, *pa_last);

      int nlist_last = nlistloc - _n_ini;

      int nzja_last = (int) (*pia_last)[nlist_last];

        _asub_last.SetNlist (nlist_last);
        _asub_last.SetNzja (nzja_last);
      if (_b_is_char)
           _asub_last.SetNzjaChar (nzja_last);
        _asub_last.SetNza (nzja_last * _size_bxb);

        _asub_last.ResizeList (nlist_last);

      _Int *plist_last = _asub_last.GetListArr ();

      int i;

      for (i = 0; i < nlist_last; i++)
           plist_last[i] = (_Int) i;

   }

// Extend sparsity by zeroes
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::ExtendSparsity (const CMatrix < _Int, _Flt > &_a_sp, CMatrix < _Int,
                              _Flt > &_a_ext) const
   {

// Open sparsities

      const int nlist_ini = this->n_list;
      const int nzja_ini = this->nz_ja;
      const int nzjachar_ini = this->nz_jachar;
      const _Int *plist = &this->list_matr[0];
      const _Int *pia = &this->ia_matr[0];
      const _Int *pja = &this->ja_matr[0];
      const char *pjachar = &this->jachar_matr[0];
      const _Flt *pa = &this->a_matr[0];

      const int nlist_sp = _a_sp.GetNlist ();
      const int nzja_sp = _a_sp.GetNzja ();
      const int nzjachar_sp = _a_sp.GetNzjaChar ();
      const _Int *plist_sp = _a_sp.GetListArr ();
      const _Int *pia_sp = _a_sp.GetIaArr ();
      const _Int *pja_sp = _a_sp.GetJaArr ();
      const char *pjachar_sp = _a_sp.GetJaCharArr ();

      bool b_is_char = false;
      if (nzjachar_ini + nzjachar_sp > 0)
      {
         if (nzjachar_ini != nzja_ini) {
            throw " CMatrix<_Int,_Flt>::ExtendSparsity: incorrect sparsity! ";
         }
         if (nzjachar_sp != nzja_sp)
         {
            throw " CMatrix<_Int,_Flt>::ExtendSparsity: incorrect sparsity! ";
         }
         b_is_char = true;
      }
// Compute extended by zeroes data

      int nlist_ext_max = nlist_ini + nlist_sp;
      int nzja_ext_max = nzja_ini + nzja_sp;
      int nzjachar_ext_max = nzjachar_ini + nzjachar_sp;

      vector < _Int > list_ext (nlist_ext_max + 1);
      vector < _Int > ia_ext (nlist_ext_max + 1);
      vector < _Int > ja_ext (nzja_ext_max + 1);
      vector < char >jachar_ext (nzjachar_ext_max + 1);
      vector < _Flt > a_ext (nzja_ext_max + 1);

      _Int *plist_ext = &list_ext[0];
      _Int *pia_ext = &ia_ext[0];
      _Int *pja_ext = &ja_ext[0];
      char *pjachar_ext = &jachar_ext[0];
      _Flt *pa_ext = &a_ext[0];

      _Flt fzero = (_Flt) 0.0e0;

      int nlist_ext = 0;
      int nzja_ext = 0;

      int ip_list, ip_list_sp, irow, irow_sp, j, jj, jj_sp, jp, jp_sp;
      char ichar1, ichar2, ichar;

      ip_list = 0;
      ip_list_sp = 0;

      pia_ext[0] = 0;

      while (ip_list < nlist_ini || ip_list_sp < nlist_sp) {
         if (ip_list < nlist_ini && ip_list_sp < nlist_sp) {
            irow = (int) plist[ip_list];
            irow_sp = (int) plist_sp[ip_list_sp];
            if (irow == irow_sp) {
               plist_ext[nlist_ext] = irow;
               nlist_ext++;
               jp = (int) pia[ip_list];
               jp_sp = (int) pia_sp[ip_list_sp];
               while (jp < pia[ip_list + 1] || jp_sp < pia_sp[ip_list_sp + 1]) {
                  if (jp < pia[ip_list + 1] && jp_sp < pia_sp[ip_list_sp + 1]) {
                     jj = (int) pja[jp];
                     jj_sp = (int) pja_sp[jp_sp];
                     if (jj == jj_sp) {
                        pja_ext[nzja_ext] = pja[jp];
                        if (b_is_char) {
                           ichar1 = pjachar[jp];
                           ichar2 = pjachar_sp[jp_sp];
                           ichar = ichar1;
                           if (ichar > ichar2)
                              ichar = ichar2;
                           pjachar_ext[nzja_ext] = ichar;
                        }
                        pa_ext[nzja_ext] = pa[jp];
                        nzja_ext++;
                        jp++;
                        jp_sp++;
                     } else if (jj < jj_sp) {
                        pja_ext[nzja_ext] = pja[jp];
                        if (b_is_char) {
                           pjachar_ext[nzja_ext] = pjachar[jp];
                        }
                        pa_ext[nzja_ext] = pa[jp];
                        nzja_ext++;
                        jp++;
                     } else {
                        pja_ext[nzja_ext] = pja_sp[jp_sp];
                        if (b_is_char) {
                           pjachar_ext[nzja_ext] = pjachar_sp[jp_sp];
                        }
                        pa_ext[nzja_ext] = fzero;
                        nzja_ext++;
                        jp_sp++;
                     }
                  } else if (jp < pia[ip_list + 1]) {
                     pja_ext[nzja_ext] = pja[jp];
                     if (b_is_char) {
                        pjachar_ext[nzja_ext] = pjachar[jp];
                     }
                     pa_ext[nzja_ext] = pa[jp];
                     nzja_ext++;
                     jp++;
                  } else {
                     pja_ext[nzja_ext] = pja_sp[jp_sp];
                     if (b_is_char) {
                        pjachar_ext[nzja_ext] = pjachar_sp[jp_sp];
                     }
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
               for (j = (int) pia[ip_list]; j < pia[ip_list + 1]; j++) {
                  pja_ext[nzja_ext] = pja[j];
                  if (b_is_char) {
                     pjachar_ext[nzja_ext] = pjachar[j];
                  }
                  pa_ext[nzja_ext] = pa[j];
                  nzja_ext++;
               }
               ip_list++;
            } else {
               plist_ext[nlist_ext] = irow_sp;
               nlist_ext++;
               for (j = (int) pia_sp[ip_list_sp]; j < pia_sp[ip_list_sp + 1]; j++) {
                  pja_ext[nzja_ext] = pja_sp[j];
                  if (b_is_char) {
                     pjachar_ext[nzja_ext] = pjachar_sp[j];
                  }
                  pa_ext[nzja_ext] = fzero;
                  nzja_ext++;
               }
               ip_list_sp++;
            }
         } else if (ip_list < nlist_ini) {
            plist_ext[nlist_ext] = plist[ip_list];
            nlist_ext++;
            for (j = (int) pia[ip_list]; j < pia[ip_list + 1]; j++) {
               pja_ext[nzja_ext] = pja[j];
               if (b_is_char) {
                  pjachar_ext[nzja_ext] = pjachar[j];
               }
               pa_ext[nzja_ext] = pa[j];
               nzja_ext++;
            }
            ip_list++;
         } else if (ip_list_sp < nlist_sp) {
            plist_ext[nlist_ext] = plist_sp[ip_list_sp];
            nlist_ext++;
            for (j = (int) pia_sp[ip_list_sp]; j < pia_sp[ip_list_sp + 1]; j++) {
               pja_ext[nzja_ext] = pja_sp[j];
               if (b_is_char) {
                  pjachar_ext[nzja_ext] = pjachar_sp[j];
               }
               pa_ext[nzja_ext] = fzero;
               nzja_ext++;
            }
            ip_list_sp++;
         }
         pia_ext[nlist_ext] = nzja_ext;
      }

// Store computed data

      CMatrix < _Int, _Flt > a_new;

      a_new.ResizeAndSetAll (nlist_ext, 0, nzja_ext, 0, nzja_ext);

      if (b_is_char) {
         a_new.ResizeJaChar (nzja_ext);
         a_new.SetNzjaChar (nzja_ext);
      }

      _Int *plist_new = a_new.GetListArr ();
      _Int *pia_new = a_new.GetIaArr ();
      _Int *pja_new = a_new.GetJaArr ();
      char *pjachar_new = a_new.GetJaCharArr ();
      _Flt *pa_new = a_new.GetAArr ();

      int i;

      for (i = 0; i < nlist_ext; i++)
         plist_new[i] = plist_ext[i];
      for (i = 0; i <= nlist_ext; i++)
         pia_new[i] = pia_ext[i];
      for (i = 0; i < nzja_ext; i++)
         pja_new[i] = pja_ext[i];
      if (b_is_char) {
         for (i = 0; i < nzja_ext; i++)
            pjachar_new[i] = pjachar_ext[i];
      }
      for (i = 0; i < nzja_ext; i++)
         pa_new[i] = pa_ext[i];

      _a_ext.ReplaceFree (a_new);

   }

// Extend sparsity by zeroes
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::ExtendSparsity_BxB (int _blksize, const CMatrix < _Int, _Flt > &_a_sp,
                                  CMatrix < _Int, _Flt > &_a_ext) const
   {

// Open sparsities

      int b_2 = _blksize * _blksize;

      const int nlist_ini = this->n_list;
      const int nzja_ini = this->nz_ja;
      const int nzjachar_ini = this->nz_jachar;
      const _Int *plist = &this->list_matr[0];
      const _Int *pia = &this->ia_matr[0];
      const _Int *pja = &this->ja_matr[0];
      const char *pjachar = &this->jachar_matr[0];
      const _Flt *pa = &this->a_matr[0];

      const int nlist_sp = _a_sp.GetNlist ();
      const int nzja_sp = _a_sp.GetNzja ();
      const int nzjachar_sp = _a_sp.GetNzjaChar ();
      const _Int *plist_sp = _a_sp.GetListArr ();
      const _Int *pia_sp = _a_sp.GetIaArr ();
      const _Int *pja_sp = _a_sp.GetJaArr ();
      const char *pjachar_sp = _a_sp.GetJaCharArr ();

      bool b_is_char = false;
      if (nzjachar_ini + nzjachar_sp > 0)
      {
         if (nzjachar_ini != nzja_ini) {
            throw " CMatrix<_Int,_Flt>::ExtendSparsity_BxB: incorrect sparsity! ";
         }
         if (nzjachar_sp != nzja_sp)
         {
            throw " CMatrix<_Int,_Flt>::ExtendSparsity_BxB: incorrect sparsity! ";
         }
         b_is_char = true;
      }
// Compute extended by zeroes data

      int nlist_ext_max = nlist_ini + nlist_sp;
      int nzja_ext_max = nzja_ini + nzja_sp;
      int nzjachar_ext_max = nzjachar_ini + nzjachar_sp;

      vector < _Int > list_ext (nlist_ext_max + 1);
      vector < _Int > ia_ext (nlist_ext_max + 1);
      vector < _Int > ja_ext (nzja_ext_max + 1);
      vector < char >jachar_ext (nzjachar_ext_max + 1);
      vector < _Flt > a_ext (nzja_ext_max * b_2 + 1);

      _Int *plist_ext = &list_ext[0];
      _Int *pia_ext = &ia_ext[0];
      _Int *pja_ext = &ja_ext[0];
      char *pjachar_ext = &jachar_ext[0];
      _Flt *pa_ext = &a_ext[0];

      int nlist_ext = 0;
      int nzja_ext = 0;

      int ip_list, ip_list_sp, irow, irow_sp, j, jj, jj_sp, jp, jp_sp;
      char ichar1, ichar2, ichar;

      ip_list = 0;
      ip_list_sp = 0;

      pia_ext[0] = 0;

      while (ip_list < nlist_ini || ip_list_sp < nlist_sp) {
         if (ip_list < nlist_ini && ip_list_sp < nlist_sp) {
            irow = (int) plist[ip_list];
            irow_sp = (int) plist_sp[ip_list_sp];
            if (irow == irow_sp) {
               plist_ext[nlist_ext] = irow;
               nlist_ext++;
               jp = (int) pia[ip_list];
               jp_sp = (int) pia_sp[ip_list_sp];
               while (jp < pia[ip_list + 1] || jp_sp < pia_sp[ip_list_sp + 1]) {
                  if (jp < pia[ip_list + 1] && jp_sp < pia_sp[ip_list_sp + 1]) {
                     jj = (int) pja[jp];
                     jj_sp = (int) pja_sp[jp_sp];
                     if (jj == jj_sp) {
                        pja_ext[nzja_ext] = pja[jp];
                        if (b_is_char) {
                           ichar1 = pjachar[jp];
                           ichar2 = pjachar_sp[jp_sp];
                           ichar = ichar1;
                           if (ichar > ichar2)
                              ichar = ichar2;
                           pjachar_ext[nzja_ext] = ichar;
                        }
                        CVector < _Flt >::CopyVector (b_2, pa + jp * b_2,
                                                      pa_ext + nzja_ext * b_2);
                        nzja_ext++;
                        jp++;
                        jp_sp++;
                     } else if (jj < jj_sp) {
                        pja_ext[nzja_ext] = pja[jp];
                        if (b_is_char) {
                           pjachar_ext[nzja_ext] = pjachar[jp];
                        }
                        CVector < _Flt >::CopyVector (b_2, pa + jp * b_2,
                                                      pa_ext + nzja_ext * b_2);
                        nzja_ext++;
                        jp++;
                     } else {
                        pja_ext[nzja_ext] = pja_sp[jp_sp];
                        if (b_is_char) {
                           pjachar_ext[nzja_ext] = pjachar_sp[jp_sp];
                        }
                        CVector < _Flt >::SetByZeroes (b_2, pa_ext + nzja_ext * b_2);
                        nzja_ext++;
                        jp_sp++;
                     }
                  } else if (jp < pia[ip_list + 1]) {
                     pja_ext[nzja_ext] = pja[jp];
                     if (b_is_char) {
                        pjachar_ext[nzja_ext] = pjachar[jp];
                     }
                     CVector < _Flt >::CopyVector (b_2, pa + jp * b_2,
                                                   pa_ext + nzja_ext * b_2);
                     nzja_ext++;
                     jp++;
                  } else {
                     pja_ext[nzja_ext] = pja_sp[jp_sp];
                     if (b_is_char) {
                        pjachar_ext[nzja_ext] = pjachar_sp[jp_sp];
                     }
                     CVector < _Flt >::SetByZeroes (b_2, pa_ext + nzja_ext * b_2);
                     nzja_ext++;
                     jp_sp++;
                  }
               }
               ip_list++;
               ip_list_sp++;
            } else if (irow < irow_sp) {
               plist_ext[nlist_ext] = irow;
               nlist_ext++;
               for (j = (int) pia[ip_list]; j < pia[ip_list + 1]; j++) {
                  pja_ext[nzja_ext] = pja[j];
                  if (b_is_char) {
                     pjachar_ext[nzja_ext] = pjachar[j];
                  }
                  CVector < _Flt >::CopyVector (b_2, pa + j * b_2,
                                                pa_ext + nzja_ext * b_2);
                  nzja_ext++;
               }
               ip_list++;
            } else {
               plist_ext[nlist_ext] = irow_sp;
               nlist_ext++;
               for (j = (int) pia_sp[ip_list_sp]; j < pia_sp[ip_list_sp + 1]; j++) {
                  pja_ext[nzja_ext] = pja_sp[j];
                  if (b_is_char) {
                     pjachar_ext[nzja_ext] = pjachar_sp[j];
                  }
                  CVector < _Flt >::SetByZeroes (b_2, pa_ext + nzja_ext * b_2);
                  nzja_ext++;
               }
               ip_list_sp++;
            }
         } else if (ip_list < nlist_ini) {
            plist_ext[nlist_ext] = plist[ip_list];
            nlist_ext++;
            for (j = (int) pia[ip_list]; j < pia[ip_list + 1]; j++) {
               pja_ext[nzja_ext] = pja[j];
               if (b_is_char) {
                  pjachar_ext[nzja_ext] = pjachar[j];
               }
               CVector < _Flt >::CopyVector (b_2, pa + j * b_2, pa_ext + nzja_ext * b_2);
               nzja_ext++;
            }
            ip_list++;
         } else if (ip_list_sp < nlist_sp) {
            plist_ext[nlist_ext] = plist_sp[ip_list_sp];
            nlist_ext++;
            for (j = (int) pia_sp[ip_list_sp]; j < pia_sp[ip_list_sp + 1]; j++) {
               pja_ext[nzja_ext] = pja_sp[j];
               if (b_is_char) {
                  pjachar_ext[nzja_ext] = pjachar_sp[j];
               }
               CVector < _Flt >::SetByZeroes (b_2, pa_ext + nzja_ext * b_2);
               nzja_ext++;
            }
            ip_list_sp++;
         }
         pia_ext[nlist_ext] = nzja_ext;
      }

// Store computed data

      CMatrix < _Int, _Flt > a_new;

      a_new.ResizeAndSetAll (nlist_ext, 0, nzja_ext, 0, nzja_ext * b_2);

      a_new.SetBSize (_blksize);

      if (b_is_char) {
         a_new.ResizeJaChar (nzja_ext);
         a_new.SetNzjaChar (nzja_ext);
      }

      _Int *plist_new = a_new.GetListArr ();
      _Int *pia_new = a_new.GetIaArr ();
      _Int *pja_new = a_new.GetJaArr ();
      char *pjachar_new = a_new.GetJaCharArr ();
      _Flt *pa_new = a_new.GetAArr ();

      int i;

      for (i = 0; i < nlist_ext; i++)
         plist_new[i] = plist_ext[i];
      for (i = 0; i <= nlist_ext; i++)
         pia_new[i] = pia_ext[i];
      for (i = 0; i < nzja_ext; i++)
         pja_new[i] = pja_ext[i];
      if (b_is_char) {
         for (i = 0; i < nzja_ext; i++)
            pjachar_new[i] = pjachar_ext[i];
      }
      for (i = 0; i < nzja_ext * b_2; i++)
         pa_new[i] = pa_ext[i];

      _a_ext.ReplaceFree (a_new);

   }

// Condense sparsity only
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::CondenseSparsitySp (int *_sp2blk, int &_icycle, int _nimax, int *_imask,
                                  CMatrix < _Int, _Flt > &_a_sp_cnd)
   {

// Open work memory

      int *pimask_w = _imask;
      int *plist_w = _imask + _nimax;
      int *pia_w = _imask + 4 * _nimax;

// Open sparsity

      int nlistloc = this->GetNlist ();
      _Int *plist = this->GetListArr ();
      _Int *pia = this->GetIaArr ();
      _Int *pja = this->GetJaArr ();

// Create the lists of nodes to be scanned per each block

      int i, j, jj, iblk, jblk;

      _icycle++;

      int nlistnew = 0;

      for (i = 0; i < nlistloc; i++) {
         j = (int) plist[i];
         jj = _sp2blk[j];
         if (pimask_w[jj] != _icycle) {
            plist_w[nlistnew] = jj;
            nlistnew++;
            pimask_w[jj] = _icycle;
         }
      }

      sort (plist_w, plist_w + nlistnew);

      for (i = 0; i < nlistnew; i++) {
         jj = plist_w[i];
         pia_w[jj] = 0;
      }

      for (i = 0; i < nlistloc; i++) {
         j = (int) plist[i];
         jj = _sp2blk[j];
         pia_w[jj]++;
      }

// Count the total number of condensed blocks

      int nz = 0;

      int nlistblk = 0;
      int nzblk = 0;

      int ilistblk;

      for (ilistblk = 0; ilistblk < nlistnew; ilistblk++) {
         iblk = plist_w[ilistblk];
         nlistblk++;
         _icycle++;
         for (i = 0; i < pia_w[iblk]; i++) {
            for (j = (int) pia[nz]; j < pia[nz + 1]; j++) {
               jj = (int) pja[j];
               jblk = _sp2blk[jj];
               if (pimask_w[jblk] != _icycle) {
                  nzblk++;
                  pimask_w[jblk] = _icycle;
               }
            }
            nz++;
         }
      }

// Allocate and init condensed matrix

      _a_sp_cnd.ResizeAndSetAll (nlistblk, 0, nzblk, 0, 0);

      _Int *plistnew = _a_sp_cnd.GetListArr ();
      _Int *pianew = _a_sp_cnd.GetIaArr ();
      _Int *pjanew = _a_sp_cnd.GetJaArr ();

      for (i = 0; i < nlistblk; i++)
         plistnew[i] = (_Int) plist_w[i];

      nzblk = 0;
      nz = 0;

      pianew[0] = 0;

      for (ilistblk = 0; ilistblk < nlistblk; ilistblk++) {
         iblk = (int) plistnew[ilistblk];
         _icycle++;
         int nlistloc = 0;
         for (i = 0; i < pia_w[iblk]; i++) {
            for (j = (int) pia[nz]; j < pia[nz + 1]; j++) {
               jj = (int) pja[j];
               jblk = _sp2blk[jj];
               if (pimask_w[jblk] != _icycle) {
                  plist_w[nlistloc] = jblk;
                  nlistloc++;
                  pimask_w[jblk] = _icycle;
               }
            }
            nz++;
         }
         sort (plist_w, plist_w + nlistloc);
         for (j = 0; j < nlistloc; j++)
            pjanew[nzblk + j] = (_Int) plist_w[j];
         nzblk += nlistloc;
         pianew[ilistblk + 1] = nzblk;
      }

   }

// Uncondense sparsity only
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::UnCondenseSparsityBlksizeSp (int _blksize, CMatrix < _Int,
                                           _Flt > &_a_sp_point)
   {

// Open sparsity

      int nlist_cnd = this->GetNlist ();
      int nzja_cnd = this->GetNzja ();
      _Int *plist_cnd = this->GetListArr ();
      _Int *pia_cnd = this->GetIaArr ();
      _Int *pja_cnd = this->GetJaArr ();

// Allocate memory

      int blksize_2 = _blksize * _blksize;

      int nlistloc = nlist_cnd * _blksize;
      int nzjaloc = nzja_cnd * blksize_2;

      _a_sp_point.ResizeAndSetAll (nlistloc, 0, nzjaloc, 0, 0);

      _Int *plistnew = _a_sp_point.GetListArr ();
      _Int *pianew = _a_sp_point.GetIaArr ();
      _Int *pjanew = _a_sp_point.GetJaArr ();

// Store arrays

      int i, j, iblk, ibeg;
      int ilist = 0;

      for (i = 0; i < nlist_cnd; i++) {
         iblk = (int) plist_cnd[i];
         ibeg = iblk * _blksize;
         for (j = 0; j < _blksize; j++) {
            plistnew[ilist] = (_Int) (ibeg + j);
            ilist++;
         }
      }
      ilist = 0;
      pianew[0] = 0;
      int nz = 0;
      int jblk, jbeg, ki, kj;
      for (i = 0; i < nlist_cnd; i++) {
         for (j = 0; j < _blksize; j++) {
            for (ki = (int) pia_cnd[i]; ki < pia_cnd[i + 1]; ki++) {
               jblk = (int) pja_cnd[ki];
               jbeg = jblk * _blksize;
               for (kj = 0; kj < _blksize; kj++) {
                  pjanew[nz] = (_Int) (jbeg + kj);
                  nz++;
               }
            }
            pianew[ilist + 1] = (_Int) nz;
            ilist++;
         }
      }

   }

// Condense sparsity and elems
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::CondenseSparsityBlksize (int _blksize, int *_sp2blk, int &_icycle,
                                       int _nimax, int *_imask, CMatrix < _Int,
                                       _Flt > &_a_cnd)
   {

// Open work memory

      int *pibs_w = _imask + 3 * _nimax;

// Open data

      int nlistloc = this->GetNlist ();
      _Int *plist = this->GetListArr ();
      _Int *pia = this->GetIaArr ();
      _Int *pja = this->GetJaArr ();
      _Flt *pa = this->GetAArr ();

// Compute condensed sparsity

      this->CondenseSparsitySp (_sp2blk, _icycle, _nimax, _imask, _a_cnd);

// Open condensed data

      int nlist_cnd = _a_cnd.GetNlist ();
      int nzja_cnd = _a_cnd.GetNzja ();
      _Int *plist_cnd = _a_cnd.GetListArr ();
      _Int *pia_cnd = _a_cnd.GetIaArr ();
      _Int *pja_cnd = _a_cnd.GetJaArr ();

      int blksize_2 = _blksize * _blksize;

      int nza_cnd = nzja_cnd * blksize_2;

      _a_cnd.ResizeA (nza_cnd);
      _a_cnd.SetNza (nza_cnd);

      _Flt *pa_cnd = _a_cnd.GetAArr ();

      CVector < _Flt >::SetByZeroes (nza_cnd, pa_cnd);

// Fill matrix data

      int ilist_cnd, j, kii, kjj, jblk, iblk, iblk1, irow, jcol, ibs;

      int ilist = 0;

      for (ilist_cnd = 0; ilist_cnd < nlist_cnd; ilist_cnd++) {
         iblk = (int) plist_cnd[ilist_cnd];
         for (j = (int) pia_cnd[ilist_cnd]; j < pia_cnd[ilist_cnd + 1]; j++) {
            jblk = (int) pja_cnd[j];
            pibs_w[jblk] = j * blksize_2;
         }
         while (ilist < nlistloc) {
            irow = (int) plist[ilist];
            iblk1 = _sp2blk[irow];
            if (iblk1 == iblk) {
               kii = irow - iblk1 * _blksize;
               for (j = (int) pia[ilist]; j < pia[ilist + 1]; j++) {
                  jcol = (int) pja[j];
                  jblk = _sp2blk[jcol];
                  ibs = pibs_w[jblk];
                  kjj = jcol - jblk * _blksize;
                  pa_cnd[ibs + kjj * _blksize + kii] = pa[j];
               }
            } else {
               break;
            }
            ilist++;
         }
      }

   }

// Uncondense sparsity and elems
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::UnCondenseSparsityBlksize (int _blksize, int *_sp2blk, int &_icycle,
                                         int _nimax, int *_imask, CMatrix < _Int,
                                         _Flt > &_a_point)
   {

// Open work memory

      int *pibs_w = _imask + 3 * _nimax;

// Open condensed data

      int nlist_cnd = this->GetNlist ();
      _Int *plist_cnd = this->GetListArr ();
      _Int *pia_cnd = this->GetIaArr ();
      _Int *pja_cnd = this->GetJaArr ();
      _Flt *pa_cnd = this->GetAArr ();

      int blksize_2 = _blksize * _blksize;

// Compute uncondensed sparsity

      this->UnCondenseSparsityBlksizeSp (_blksize, _a_point);

// Open data

      int nlistloc = _a_point.GetNlist ();
      int nzjaloc = _a_point.GetNzja ();
      _Int *plist = _a_point.GetListArr ();
      _Int *pia = _a_point.GetIaArr ();
      _Int *pja = _a_point.GetJaArr ();

      _a_point.ResizeA (nzjaloc);
      _a_point.SetNza (nzjaloc);

      _Flt *pa = _a_point.GetAArr ();

      CVector < _Flt >::SetByZeroes (nzjaloc, pa);

// Fill matrix data

      int ilist_cnd, j, kii, kjj, jblk, iblk, iblk1, irow, jcol, ibs;

      int ilist = 0;

      for (ilist_cnd = 0; ilist_cnd < nlist_cnd; ilist_cnd++) {
         iblk = (int) plist_cnd[ilist_cnd];
         for (j = (int) pia_cnd[ilist_cnd]; j < pia_cnd[ilist_cnd + 1]; j++) {
            jblk = (int) pja_cnd[j];
            pibs_w[jblk] = j * blksize_2;
         }
         while (ilist < nlistloc) {
            irow = (int) plist[ilist];
            iblk1 = _sp2blk[irow];
            if (iblk1 == iblk) {
               kii = irow - iblk1 * _blksize;
               for (j = (int) pia[ilist]; j < pia[ilist + 1]; j++) {
                  jcol = (int) pja[j];
                  jblk = _sp2blk[jcol];
                  ibs = pibs_w[jblk];
                  kjj = jcol - jblk * _blksize;
                  pa[j] = pa_cnd[ibs + kjj * _blksize + kii];
               }
            } else {
               break;
            }
            ilist++;
         }
      }

   }

// Copy constructor
//========================================================================================
   template < typename _Int, typename _Flt > CMatrix < _Int,
      _Flt >::CMatrix (const CMatrix < _Int, _Flt > &_aa)
   {

      int bsizeloc = _aa.GetBSize ();
      int nlistloc = _aa.GetNlist ();
      int nlist2loc = _aa.GetNlist2 ();
      int nzjaloc = _aa.GetNzja ();
      int nzja2loc = _aa.GetNzja2 ();
      int nzjacharloc = _aa.GetNzjaChar ();
      int nzaloc = _aa.GetNza ();

      const _Int *plist_aa = &(_aa.list_matr[0]);
      const _Int *plist2_aa = &(_aa.list2_matr[0]);
      const _Int *pia_aa = &(_aa.ia_matr[0]);
      const _Int *pja_aa = &(_aa.ja_matr[0]);
      const _Int *pja2_aa = &(_aa.ja2_matr[0]);
      const char *pjachar_aa = &(_aa.jachar_matr[0]);
      const _Flt *pa_aa = &(_aa.a_matr[0]);

      this->b_size = bsizeloc;
      this->n_list = nlistloc;
      this->n_list2 = nlist2loc;
      this->nz_ja = nzjaloc;
      this->nz_ja2 = nzja2loc;
      this->nz_jachar = nzjacharloc;
      this->nz_a = nzaloc;

      this->list_matr.resize (nlistloc + 1);
      this->list2_matr.resize (nlist2loc + 1);
      this->ia_matr.resize (nlistloc + 1);
      this->ja_matr.resize (nzjaloc + 1);
      this->ja2_matr.resize (nzja2loc + 1);
      this->jachar_matr.resize (nzjacharloc + 1);
      this->a_matr.resize (nzaloc + 1);

      _Int *plist = &(this->list_matr[0]);
      _Int *plist2 = &(this->list2_matr[0]);
      _Int *pia = &(this->ia_matr[0]);
      _Int *pja = &(this->ja_matr[0]);
      _Int *pja2 = &(this->ja2_matr[0]);
      char *pjachar = &(this->jachar_matr[0]);
      _Flt *pa = &(this->a_matr[0]);

      int i;

      for (i = 0; i < nlistloc; i++)
         plist[i] = plist_aa[i];
      for (i = 0; i <= nlistloc; i++)
         pia[i] = pia_aa[i];
      for (i = 0; i < nlist2loc; i++)
         plist2[i] = plist2_aa[i];
      for (i = 0; i < nzjaloc; i++)
         pja[i] = pja_aa[i];
      for (i = 0; i < nzja2loc; i++)
         pja2[i] = pja2_aa[i];
      for (i = 0; i < nzjacharloc; i++)
         pjachar[i] = pjachar_aa[i];

      for (i = 0; i < nzaloc; i++)
         pa[i] = pa_aa[i];

   }

// Equality operator
//========================================================================================
   template < typename _Int, typename _Flt > CMatrix < _Int, _Flt > &CMatrix < _Int,
      _Flt >::operator= (const CMatrix < _Int, _Flt > &_aa)
   {

      int bsizeloc = _aa.GetBSize ();
      int nlistloc = _aa.GetNlist ();
      int nlist2loc = _aa.GetNlist2 ();
      int nzjaloc = _aa.GetNzja ();
      int nzja2loc = _aa.GetNzja2 ();
      int nzjacharloc = _aa.GetNzjaChar ();
      int nzaloc = _aa.GetNza ();

      const _Int *plist_aa = &(_aa.list_matr[0]);
      const _Int *plist2_aa = &(_aa.list2_matr[0]);
      const _Int *pia_aa = &(_aa.ia_matr[0]);
      const _Int *pja_aa = &(_aa.ja_matr[0]);
      const _Int *pja2_aa = &(_aa.ja2_matr[0]);
      const char *pjachar_aa = &(_aa.jachar_matr[0]);
      const _Flt *pa_aa = &(_aa.a_matr[0]);

      this->b_size = bsizeloc;
      this->n_list = nlistloc;
      this->n_list2 = nlist2loc;
      this->nz_ja = nzjaloc;
      this->nz_ja2 = nzja2loc;
      this->nz_jachar = nzjacharloc;
      this->nz_a = nzaloc;

      this->list_matr.resize (nlistloc + 1);
      this->list2_matr.resize (nlist2loc + 1);
      this->ia_matr.resize (nlistloc + 1);
      this->ja_matr.resize (nzjaloc + 1);
      this->ja2_matr.resize (nzja2loc + 1);
      this->jachar_matr.resize (nzjacharloc + 1);
      this->a_matr.resize (nzaloc + 1);

      _Int *plist = &(this->list_matr[0]);
      _Int *plist2 = &(this->list2_matr[0]);
      _Int *pia = &(this->ia_matr[0]);
      _Int *pja = &(this->ja_matr[0]);
      _Int *pja2 = &(this->ja2_matr[0]);
      char *pjachar = &(this->jachar_matr[0]);
      _Flt *pa = &(this->a_matr[0]);

      int i;

      for (i = 0; i < nlistloc; i++)
         plist[i] = plist_aa[i];
      for (i = 0; i <= nlistloc; i++)
         pia[i] = pia_aa[i];
      for (i = 0; i < nlist2loc; i++)
         plist2[i] = plist2_aa[i];
      for (i = 0; i < nzjaloc; i++)
         pja[i] = pja_aa[i];
      for (i = 0; i < nzja2loc; i++)
         pja2[i] = pja2_aa[i];
      for (i = 0; i < nzjacharloc; i++)
         pjachar[i] = pjachar_aa[i];

      for (i = 0; i < nzaloc; i++)
         pa[i] = pa_aa[i];

      return *this;

   }

// Get sparsity
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::GetSparsity (const CMatrix < _Int, _Flt > &_a_sp)
   {

      const int nlistloc = _a_sp.GetNlist ();
      const int nlist2loc = _a_sp.GetNlist2 ();
      const int nzjaloc = _a_sp.GetNzja ();
      const int nzja2loc = _a_sp.GetNzja2 ();
      const int nzjacharloc = _a_sp.GetNzjaChar ();

      const _Int *plist_aa = &(_a_sp.list_matr[0]);
      const _Int *plist2_aa = &(_a_sp.list2_matr[0]);
      const _Int *pia_aa = &(_a_sp.ia_matr[0]);
      const _Int *pja_aa = &(_a_sp.ja_matr[0]);
      const _Int *pja2_aa = &(_a_sp.ja2_matr[0]);
      const char *pjachar_aa = &(_a_sp.jachar_matr[0]);

      this->n_list = nlistloc;
      this->n_list2 = nlist2loc;
      this->nz_ja = nzjaloc;
      this->nz_ja2 = nzja2loc;
      this->nz_jachar = nzjacharloc;
      this->nz_a = 0;

      this->list_matr.resize (nlistloc + 1);
      this->list2_matr.resize (nlist2loc + 1);
      this->ia_matr.resize (nlistloc + 1);
      this->ja_matr.resize (nzjaloc + 1);
      this->ja2_matr.resize (nzja2loc + 1);
      this->jachar_matr.resize (nzjacharloc + 1);
      this->a_matr.resize (1);

      _Int *plist = &(this->list_matr[0]);
      _Int *plist2 = &(this->list2_matr[0]);
      _Int *pia = &(this->ia_matr[0]);
      _Int *pja = &(this->ja_matr[0]);
      _Int *pja2 = &(this->ja2_matr[0]);
      char *pjachar = &(this->jachar_matr[0]);

      int i;

      for (i = 0; i < nlistloc; i++)
         plist[i] = plist_aa[i];
      for (i = 0; i <= nlistloc; i++)
         pia[i] = pia_aa[i];
      for (i = 0; i < nlist2loc; i++)
         plist2[i] = plist2_aa[i];
      for (i = 0; i < nzjaloc; i++)
         pja[i] = pja_aa[i];
      for (i = 0; i < nzja2loc; i++)
         pja2[i] = pja2_aa[i];
      for (i = 0; i < nzjacharloc; i++)
         pjachar[i] = pjachar_aa[i];

   }

/// @brief Create simple matrix test (diag and dense)
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::CreateTestMatrixSimple (int _itype, int _nx)
   {

      if (_itype == 0) {

         this->ResizeAndSetAll (_nx, 0, _nx, 0, _nx);

         _Int *plist = &(this->list_matr[0]);
         _Int *pia = &(this->ia_matr[0]);
         _Int *pja = &(this->ja_matr[0]);
         _Flt *pa = this->GetAArr ();

         int i;

         pia[0] = 0;

         for (i = 0; i < _nx; i++) {
            plist[i] = i;
            pia[i + 1] = i + 1;
            pja[i] = i;
            pa[i] = (_Flt) ((i % 4) + 3.3e0);
         }

      } else {

         this->ResizeAndSetAll (_nx, 0, _nx * _nx, 0, _nx * _nx);

         _Int *plist = &(this->list_matr[0]);
         _Int *pia = &(this->ia_matr[0]);
         _Int *pja = &(this->ja_matr[0]);
         _Flt *pa = this->GetAArr ();

         int i, j;
         double aux;

         pia[0] = 0;

         int nz = 0;

         for (i = 0; i < _nx; i++) {
            plist[i] = i;
            for (j = 0; j < _nx; j++) {
               pja[nz] = j;
               aux = (double) (10 * i + j);
               if (i == j)
                  aux += 1.0e2;
//               aux = 2.0e0 / ((double)(i-j)+0.5e0);
               pa[nz] = (_Flt) aux;
               nz++;
            }
            pia[i + 1] = nz;
         }

      }

   }

/// @brief Create sparsity of a part of test matrix (7 point)
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::Sparsity3D (int _nx, int _ny, int _nz, int _ibegz, int _iendz)
   {

// Allocate the memory

      int nzloc = _iendz - _ibegz + 1;

      int nloc = _nx * _ny * nzloc;
      int nxy = _nx * _ny;

      this->ResizeAndSetAllSp (nloc, 0, nloc * 7, 0);

      _Int *plist = this->GetListArr ();
      _Int *pia = this->GetIaArr ();
      _Int *pja = this->GetJaArr ();

      pia[0] = 0;

      int i;

      for (i = 0; i < nloc * 7; i++)
         pja[i] = -1;

      for (i = 1; i <= nloc; i++)
         pia[i] = 1;

      long long ix, iy, iz;
      long long ind;

      for (iz = _ibegz; iz <= _iendz; iz++) {
         for (iy = 0; iy < _ny; iy++) {
            for (ix = 0; ix < _nx; ix++) {
               ind = (iz - _ibegz) * nxy + iy * _nx + ix;
               if (ix > 0)
                  pia[ind + 1]++;
               if (ix < _nx - 1)
                  pia[ind + 1]++;
               if (iy > 0)
                  pia[ind + 1]++;
               if (iy < _ny - 1)
                  pia[ind + 1]++;
               if (iz > 0)
                  pia[ind + 1]++;
               if (iz < _nz - 1)
                  pia[ind + 1]++;
            }
         }
      }

      for (i = 0; i < nloc; i++)
         pia[i + 1] = pia[i] + pia[i + 1];

      CVectorData < int >iptr (nloc + 1);

      int *piptr = iptr.Ptr ();

      for (i = 0; i < nloc; i++)
         piptr[i] = (int) pia[i];

      long long k, jj;

      for (iz = _ibegz; iz <= _iendz; iz++) {
         for (iy = 0; iy < _ny; iy++) {
            for (ix = 0; ix < _nx; ix++) {
               ind = (iz - _ibegz) * nxy + iy * _nx + ix;
               jj = iz * nxy + iy * _nx + ix;
               k = piptr[ind];
               pja[k] = (_Int) jj;
               piptr[ind]++;
            }
         }
      }

      long long ind1;

      for (iz = _ibegz; iz <= _iendz; iz++) {
         for (iy = 0; iy < _ny; iy++) {
            for (ix = 0; ix < _nx; ix++) {
               ind = (iz - _ibegz) * nxy + iy * _nx + ix;
               jj = iz * nxy + iy * _nx + ix;
               if (ix > 0) {
                  k = piptr[ind];
                  ind1 = iz * nxy + iy * _nx + (ix - 1);
                  pja[k] = (_Int) ind1;
                  piptr[ind]++;
               }
               if (ix < _nx - 1) {
                  k = piptr[ind];
                  ind1 = iz * nxy + iy * _nx + (ix + 1);
                  pja[k] = (_Int) ind1;
                  piptr[ind]++;
               }
               if (iy > 0) {
                  k = piptr[ind];
                  ind1 = iz * nxy + (iy - 1) * _nx + ix;
                  pja[k] = (_Int) ind1;
                  piptr[ind]++;
               }
               if (iy < _ny - 1) {
                  k = piptr[ind];
                  ind1 = iz * nxy + (iy + 1) * _nx + ix;
                  pja[k] = (_Int) ind1;
                  piptr[ind]++;
               }
               if (iz > 0) {
                  k = piptr[ind];
                  ind1 = (iz - 1) * nxy + iy * _nx + ix;
                  pja[k] = (_Int) ind1;
                  piptr[ind]++;
               }
               if (iz < _nz - 1) {
                  k = piptr[ind];
                  ind1 = (iz + 1) * nxy + iy * _nx + ix;
                  pja[k] = (_Int) ind1;
                  piptr[ind]++;
               }
               sort (pja + pia[ind], pja + pia[ind + 1]);
            }
         }
      }

      for (iz = _ibegz; iz <= _iendz; iz++) {
         for (iy = 0; iy < _ny; iy++) {
            for (ix = 0; ix < _nx; ix++) {
               ind = (iz - _ibegz) * nxy + iy * _nx + ix;
               jj = iz * nxy + iy * _nx + ix;
               plist[ind] = (_Int) jj;
            }
         }
      }

      this->SetNzja ((int) pia[nloc]);

   }

/// @brief Create a part of test matrix (7 point)
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::CreateTestMatrix3D (int _nx, int _ny, int _nz, int _ibegz, int _iendz,
                                  double _eps, double _dunsy)
   {

      this->Sparsity3D (_nx, _ny, _nz, _ibegz, _iendz);

      this->SetNza (this->nz_ja);
      this->ResizeA (this->nz_ja);

      _Flt *pa = this->GetAArr ();

// Init a

      double eps2 = _eps * _eps;

      double s[] = { -1.0e0, -_eps, -eps2 - _dunsy / 2.0e0, 2.0e0 * (1.0e0 + _eps + eps2),
         -eps2 + _dunsy / 2.0e0, -_eps, -1.0e0
      };

      int is[] = { 0, 0, -1, 0, 1, 0, 0 };
      int js[] = { 0, -1, 0, 0, 0, 1, 0 };
      int ks[] = { -1, 0, 0, 0, 0, 0, 1 };

      double d;

      int ii, jj, kk;
      int i, j, k, t;
      long long iii, jjj, indd = 0;
      int nxy = _nx * _ny;
      int nz = 0;

      for (kk = _ibegz; kk <= _iendz; kk++) {
         for (jj = 0; jj < _ny; jj++) {
            for (ii = 0; ii < _nx; ii++) {
               iii = kk * nxy + jj * _nx + ii;
               d = 0.0e0;
               for (t = 0; t < 7; t++) {
                  i = ii + is[t];
                  j = jj + js[t];
                  k = kk + ks[t];
                  jjj = k * nxy + j * _nx + i;
                  if (iii == jjj)
                     indd = nz;
                  if (k < 0)
                     d += s[0];
                  if (j < 0)
                     d += s[1];
                  if (j >= _ny)
                     d += s[5];
                  if (k >= _nz)
                     d += s[6];
                  if (i >= 0 && i < _nx && j >= 0 && j < _ny && k >= 0 && k < _nz) {
                     pa[nz] = (_Flt) s[t];
                     nz++;
                  }
               }
               pa[indd] += (_Flt) d;
            }
         }
      }

   }

/// @brief Create sparsity of a part of simple 2D test matrix (5-point Laplace)
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::Sparsity2D (int _nx, int _ny, int _ibegy, int _iendy)
   {

// Allocate the memory

      int nyloc = _iendy - _ibegy + 1;

      int nloc = _nx * nyloc;
      int nx = _nx;

      this->ResizeList (nloc);
      this->ResizeIa (nloc + 1);
      this->SetNlist (nloc);

      _Int *plist = this->GetListArr ();
      _Int *pia = this->GetIaArr ();

      pia[0] = 0;

      int i;

      for (i = 1; i <= nloc; i++)
         pia[i] = 1;

      int ix, iy, ind;

      for (iy = _ibegy; iy <= _iendy; iy++) {
         for (ix = 0; ix < _nx; ix++) {
            ind = (iy - _ibegy) * nx + ix;
            if (ix > 0)
               pia[ind + 1]++;
            if (ix < _nx - 1)
               pia[ind + 1]++;
            if (iy > 0)
               pia[ind + 1]++;
            if (iy < _ny - 1)
               pia[ind + 1]++;
         }
      }

      for (i = 0; i < nloc; i++)
         pia[i + 1] = pia[i] + pia[i + 1];

      int nzjaloc = (int) pia[nloc];

      this->ResizeJa (nzjaloc);
      this->SetNzja (nzjaloc);

      _Int *pja = this->GetJaArr ();

      for (i = 0; i < nzjaloc; i++)
         pja[i] = -1;

      CVectorData < int >iptr (nloc);

      int *piptr = iptr.Ptr ();

      for (i = 0; i < nloc; i++)
         piptr[i] = (int) pia[i];

      int k, jj;

      for (iy = _ibegy; iy <= _iendy; iy++) {
         for (ix = 0; ix < _nx; ix++) {
            ind = (iy - _ibegy) * nx + ix;
            jj = iy * _nx + ix;
            k = piptr[ind];
            pja[k] = (_Int) jj;
            piptr[ind]++;
         }
      }

      int ind1;

      for (iy = _ibegy; iy <= _iendy; iy++) {
         for (ix = 0; ix < _nx; ix++) {
            ind = (iy - _ibegy) * nx + ix;
            jj = iy * _nx + ix;
            if (ix > 0) {
               k = piptr[ind];
               ind1 = iy * _nx + (ix - 1);
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            if (ix < _nx - 1) {
               k = piptr[ind];
               ind1 = iy * _nx + (ix + 1);
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            if (iy > 0) {
               k = piptr[ind];
               ind1 = (iy - 1) * _nx + ix;
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            if (iy < _ny - 1) {
               k = piptr[ind];
               ind1 = (iy + 1) * _nx + ix;
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            sort (pja + pia[ind], pja + pia[ind + 1]);
         }
      }

      for (iy = _ibegy; iy <= _iendy; iy++) {
         for (ix = 0; ix < _nx; ix++) {
            ind = (iy - _ibegy) * nx + ix;
            jj = iy * _nx + ix;
            plist[ind] = (_Int) jj;
         }
      }

   }

/// @brief Create sparsity of a part of simple 2D test matrix (5-point Laplace)
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::CreateTestMatrix2DLaplace (int _nx, int _ny, int _ibegy, int _iendy)
   {

      this->Sparsity2D (_nx, _ny, _ibegy, _iendy);

      this->SetNza (this->nz_ja);
      this->ResizeA (this->nz_ja);

      _Flt *pa = this->GetAArr ();

// Init a

      int is[] = { 0, -1, 0, 1, 0 };
      int js[] = { -1, 0, 0, 0, 1 };
      double s[] = { -1, -1, 4, -1, -1 };

      int ii, jj;
      int i, j, t;
      int nz = 0;

      for (jj = _ibegy; jj <= _iendy; jj++) {
         for (ii = 0; ii < _nx; ii++) {
            for (t = 0; t < 5; t++) {
               i = ii + is[t];
               j = jj + js[t];
               if (i >= 0 && i < _nx && j >= 0 && j < _ny) {
                  pa[nz] = (_Flt) s[t];
                  nz++;
               }
            }
         }
      }

   }

/// @brief Create sparsity of a part of complicated 2D test matrix (13-point Laplace^2)
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::SparsityExtended2D (int _nx, int _ny, int _ibegy, int _iendy)
   {

// Allocate the memory

      int nyloc = _iendy - _ibegy + 1;

      int nloc = _nx * nyloc;
      int nx = _nx;

      this->ResizeList (nloc);
      this->ResizeIa (nloc + 1);
      this->SetNlist (nloc);

      _Int *plist = this->GetListArr ();
      _Int *pia = this->GetIaArr ();

      pia[0] = 0;

      int i;

      for (i = 1; i <= nloc; i++)
         pia[i] = 1;

      int ix, iy, ind;

      for (iy = _ibegy; iy <= _iendy; iy++) {
         for (ix = 0; ix < _nx; ix++) {
            ind = (iy - _ibegy) * nx + ix;
            if (ix > 0)
               pia[ind + 1]++;
            if (ix > 1)
               pia[ind + 1]++;
            if (ix < _nx - 1)
               pia[ind + 1]++;
            if (ix < _nx - 2)
               pia[ind + 1]++;
            if (iy > 0)
               pia[ind + 1]++;
            if (iy > 1)
               pia[ind + 1]++;
            if (iy < _ny - 1)
               pia[ind + 1]++;
            if (iy < _ny - 2)
               pia[ind + 1]++;
            if (ix > 0 && iy > 0)
               pia[ind + 1]++;
            if (ix > 0 && iy < _ny - 1)
               pia[ind + 1]++;
            if (ix < _nx - 1 && iy > 0)
               pia[ind + 1]++;
            if (ix < _nx - 1 && iy < _ny - 1)
               pia[ind + 1]++;
         }
      }

      for (i = 0; i < nloc; i++)
         pia[i + 1] = pia[i] + pia[i + 1];

      int nzjaloc = (int) pia[nloc];

      this->ResizeJa (nzjaloc);
      this->SetNzja (nzjaloc);

      _Int *pja = this->GetJaArr ();

      for (i = 0; i < nzjaloc; i++)
         pja[i] = -1;

      CVectorData < int >iptr (nloc);

      int *piptr = iptr.Ptr ();

      for (i = 0; i < nloc; i++)
         piptr[i] = (int) pia[i];

      int k, jj;

      for (iy = _ibegy; iy <= _iendy; iy++) {
         for (ix = 0; ix < _nx; ix++) {
            ind = (iy - _ibegy) * nx + ix;
            jj = iy * _nx + ix;
            k = piptr[ind];
            pja[k] = (_Int) jj;
            piptr[ind]++;
         }
      }

      int ind1;

      for (iy = _ibegy; iy <= _iendy; iy++) {
         for (ix = 0; ix < _nx; ix++) {
            ind = (iy - _ibegy) * nx + ix;
            jj = iy * _nx + ix;
            if (ix > 0) {
               k = piptr[ind];
               ind1 = iy * _nx + (ix - 1);
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            if (ix > 1) {
               k = piptr[ind];
               ind1 = iy * _nx + (ix - 2);
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            if (ix < _nx - 1) {
               k = piptr[ind];
               ind1 = iy * _nx + (ix + 1);
               pja[k] = ind1;
               piptr[ind]++;
            }
            if (ix < _nx - 2) {
               k = piptr[ind];
               ind1 = iy * _nx + (ix + 2);
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            if (iy > 0) {
               k = piptr[ind];
               ind1 = (iy - 1) * _nx + ix;
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            if (iy > 1) {
               k = piptr[ind];
               ind1 = (iy - 2) * _nx + ix;
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            if (iy < _ny - 1) {
               k = piptr[ind];
               ind1 = (iy + 1) * _nx + ix;
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            if (iy < _ny - 2) {
               k = piptr[ind];
               ind1 = (iy + 2) * _nx + ix;
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            if (ix > 0 && iy > 0) {
               k = piptr[ind];
               ind1 = (iy - 1) * _nx + ix - 1;
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            if (ix > 0 && iy < _ny - 1) {
               k = piptr[ind];
               ind1 = (iy + 1) * _nx + ix - 1;
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            if (ix < _nx - 1 && iy > 0) {
               k = piptr[ind];
               ind1 = (iy - 1) * _nx + ix + 1;
               pja[k] = ind1;
               piptr[ind]++;
            }
            if (ix < _nx - 1 && iy < _ny - 1) {
               k = piptr[ind];
               ind1 = (iy + 1) * _nx + ix + 1;
               pja[k] = (_Int) ind1;
               piptr[ind]++;
            }
            sort (pja + pia[ind], pja + pia[ind + 1]);
         }
      }

      for (iy = _ibegy; iy <= _iendy; iy++) {
         for (ix = 0; ix < _nx; ix++) {
            ind = (iy - _ibegy) * nx + ix;
            jj = iy * _nx + ix;
            plist[ind] = (_Int) jj;
         }
      }

   }

/// @brief Create a part of complicated 2D test matrix (13-point Laplace^2)
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::CreateTestMatrix2D (int _nx, int _ny, int _ibegy, int _iendy)
   {

      this->SparsityExtended2D (_nx, _ny, _ibegy, _iendy);

      this->SetNza (nz_ja);
      this->ResizeA (nz_ja);

      _Flt *pa = this->GetAArr ();

// Init a

      int is[] = { 0, -1, 0, 1, -2, -1, 0, 1, 2, -1, 0, 1, 0 };
      int js[] = { -2, -1, -1, -1, 0, 0, 0, 0, 0, 1, 1, 1, 2 };
      double s[] = { 1, 2, -8, 2, 1, -8, 20, -8, 1, 2, -8, 2, 1 };

      int ii, jj;
      int i, j;
      int t;
      int nz = 0;

      for (jj = _ibegy; jj <= _iendy; jj++) {
         for (ii = 0; ii < _nx; ii++) {
            for (t = 0; t < 13; t++) {
               i = ii + is[t];
               j = jj + js[t];
               if (i >= 0 && i < _nx && j >= 0 && j < _ny) {
                  pa[nz] = (_Flt) s[t];
                  nz++;
               }
            }
         }
      }

   }

// Compute transposed sparsity for incomplete list
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::TransposedSparsityListSp (int &_icycle, int *_imask, int *_indarr,
                                        int *_iptr, int *_listloc, int *_ialoc,
                                        CMatrix < _Int, _Flt > &_at) const
   {

      const int nlist_ini = this->n_list;
      const int nzja_ini = this->nz_ja;
      const int nzjachar_ini = this->nz_jachar;
      const _Int *plist = &this->list_matr[0];
      const _Int *pia = &this->ia_matr[0];
      const _Int *pja = &this->ja_matr[0];
      const char *pjachar = &this->jachar_matr[0];

      bool b_is_char = false;
      if (nzjachar_ini > 0)
      {
         if (nzjachar_ini != nzja_ini) {
            throw " CMatrix<_Int,_Flt>::TransposedSparsityListSp: incorrect sparsity! ";
         }
         b_is_char = true;
      }
// Count the sparsity of at

      int nlistloc = 0;

      _icycle++;

      int i, j, jj;

      for (i = 0; i < nlist_ini; i++) {
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jj = (int) pja[j];
            if (_imask[jj] != _icycle) {
               _listloc[nlistloc] = jj;
               nlistloc++;
               _imask[jj] = _icycle;
               _indarr[jj] = 0;
            }
            _indarr[jj]++;
         }
      }

      sort (_listloc, _listloc + nlistloc);

      _ialoc[0] = 0;

      for (i = 0; i < nlistloc; i++) {
         jj = _listloc[i];
         _ialoc[i + 1] = _ialoc[i] + _indarr[jj];
         _iptr[i] = _ialoc[i];
         _indarr[jj] = i;
      }

// Init transposed matrix

      CMatrix < _Int, _Flt > at;

      at.ResizeAndSetAllSp (nlistloc, 0, nzja_ini, 0);

      if (b_is_char) {
         at.ResizeJaChar (nzja_ini);
         at.SetNzjaChar (nzja_ini);
      }

      _Int *plist_at = at.GetListArr ();
      _Int *pia_at = at.GetIaArr ();
      _Int *pja_at = at.GetJaArr ();
      char *pjachar_at = at.GetJaCharArr ();

      for (i = 0; i < nlistloc; i++)
         plist_at[i] = (_Int) _listloc[i];
      for (i = 0; i <= nlistloc; i++)
         pia_at[i] = (_Int) _ialoc[i];

      int k, irow, ind;

      for (i = 0; i < nlist_ini; i++) {
         irow = (int) plist[i];
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jj = (int) pja[j];
            ind = _indarr[jj];
            _iptr[ind]++;
            k = _iptr[ind];
            pja_at[k - 1] = (_Int) irow;
            if (b_is_char) {
               pjachar_at[k - 1] = pjachar[j];
            }
         }
      }

      _at.ReplaceFree (at);

   }

// Compute transposed matrix for incomplete list
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::TransposedSparsityList (int &_icycle, int *_imask, int *_indarr, int *_iptr,
                                      int *_listloc, int *_ialoc, CMatrix < _Int,
                                      _Flt > &_at) const
   {

      const int nlist_ini = this->n_list;
      const int nzja_ini = this->nz_ja;
      const int nzjachar_ini = this->nz_jachar;
      const _Int *plist = &this->list_matr[0];
      const _Int *pia = &this->ia_matr[0];
      const _Int *pja = &this->ja_matr[0];
      const char *pjachar = &this->jachar_matr[0];
      const _Flt *pa = &this->a_matr[0];

      bool b_is_char = false;
      if (nzjachar_ini > 0)
      {
         if (nzjachar_ini != nzja_ini) {
            throw " CMatrix<_Int,_Flt>::TransposedSparsityList: incorrect sparsity! ";
         }
         b_is_char = true;
      }
// Count the sparsity of at

      int nlistloc = 0;

      _icycle++;

      int i, j, jj;

      for (i = 0; i < nlist_ini; i++) {
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jj = (int) pja[j];
            if (_imask[jj] != _icycle) {
               _listloc[nlistloc] = jj;
               nlistloc++;
               _imask[jj] = _icycle;
               _indarr[jj] = 0;
            }
            _indarr[jj]++;
         }
      }

      sort (_listloc, _listloc + nlistloc);

      _ialoc[0] = 0;

      for (i = 0; i < nlistloc; i++) {
         jj = _listloc[i];
         _ialoc[i + 1] = _ialoc[i] + _indarr[jj];
         _iptr[i] = _ialoc[i];
         _indarr[jj] = i;
      }

// Init transposed matrix

      CMatrix < _Int, _Flt > at;

      at.ResizeAndSetAll (nlistloc, 0, nzja_ini, 0, nzja_ini);

      if (b_is_char) {
         at.ResizeJaChar (nzja_ini);
         at.SetNzjaChar (nzja_ini);
      }

      _Int *plist_at = at.GetListArr ();
      _Int *pia_at = at.GetIaArr ();
      _Int *pja_at = at.GetJaArr ();
      char *pjachar_at = at.GetJaCharArr ();
      _Flt *pa_at = at.GetAArr ();

      for (i = 0; i < nlistloc; i++)
         plist_at[i] = (_Int) _listloc[i];
      for (i = 0; i <= nlistloc; i++)
         pia_at[i] = (_Int) _ialoc[i];

      int k, irow, ind;

      for (i = 0; i < nlist_ini; i++) {
         irow = (int) plist[i];
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jj = (int) pja[j];
            ind = _indarr[jj];
            _iptr[ind]++;
            k = _iptr[ind];
            pja_at[k - 1] = (_Int) irow;
            if (b_is_char) {
               pjachar_at[k - 1] = pjachar[j];
            }
            pa_at[k - 1] = pa[j];
         }
      }

      _at.ReplaceFree (at);

   }

// Compute transposed matrix for incomplete list
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::TransposedSparsityList_BxB (int _blksize, int &_icycle, int *_imask,
                                          int *_indarr, int *_iptr, int *_listloc,
                                          int *_ialoc, CMatrix < _Int, _Flt > &_at) const
   {

      int b_2 = _blksize * _blksize;

      const int nlist_ini = this->n_list;
      const int nzja_ini = this->nz_ja;
      const int nzjachar_ini = this->nz_jachar;
      const _Int *plist = &this->list_matr[0];
      const _Int *pia = &this->ia_matr[0];
      const _Int *pja = &this->ja_matr[0];
      const char *pjachar = &this->jachar_matr[0];
      const _Flt *pa = &this->a_matr[0];

      bool b_is_char = false;
      if (nzjachar_ini > 0)
      {
         if (nzjachar_ini != nzja_ini) {
            throw " CMatrix<_Int,_Flt>::TransposedSparsityList: incorrect sparsity! ";
         }
         b_is_char = true;
      }
// Count the sparsity of at

      int nlistloc = 0;

      _icycle++;

      int i, j, jj;

      for (i = 0; i < nlist_ini; i++) {
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jj = (int) pja[j];
            if (_imask[jj] != _icycle) {
               _listloc[nlistloc] = jj;
               nlistloc++;
               _imask[jj] = _icycle;
               _indarr[jj] = 0;
            }
            _indarr[jj]++;
         }
      }

      sort (_listloc, _listloc + nlistloc);

      _ialoc[0] = 0;

      for (i = 0; i < nlistloc; i++) {
         jj = _listloc[i];
         _ialoc[i + 1] = _ialoc[i] + _indarr[jj];
         _iptr[i] = _ialoc[i];
         _indarr[jj] = i;
      }

// Init transposed matrix

      CMatrix < _Int, _Flt > at;

      at.ResizeAndSetAll (nlistloc, 0, nzja_ini, 0, nzja_ini * b_2);

      if (b_is_char) {
         at.ResizeJaChar (nzja_ini);
         at.SetNzjaChar (nzja_ini);
      }

      _Int *plist_at = at.GetListArr ();
      _Int *pia_at = at.GetIaArr ();
      _Int *pja_at = at.GetJaArr ();
      char *pjachar_at = at.GetJaCharArr ();
      _Flt *pa_at = at.GetAArr ();

      for (i = 0; i < nlistloc; i++)
         plist_at[i] = (_Int) _listloc[i];
      for (i = 0; i <= nlistloc; i++)
         pia_at[i] = (_Int) _ialoc[i];

      int k, irow, ind;

      for (i = 0; i < nlist_ini; i++) {
         irow = (int) plist[i];
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jj = (int) pja[j];
            ind = _indarr[jj];
            _iptr[ind]++;
            k = _iptr[ind];
            pja_at[k - 1] = (_Int) irow;
            if (b_is_char) {
               pjachar_at[k - 1] = pjachar[j];
            }
            CVector < _Flt >::TransposeBlock (_blksize, _blksize, pa + j * b_2, _blksize,
                                              pa_at + (k - 1) * b_2, _blksize);
         }
      }

      _at.ReplaceFree (at);

   }

// Compute symmetrized matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::SymmetrizeAndAddZeroes (CMatrix < _Int, _Flt > &_asymm) const
   {

      const int nlist_ini = this->n_list;

// Split L and U

        vector < _Int > ia_L;
        vector < _Int > ja_L;
        vector < char >jachar_L;
        vector < _Flt > a_L;

        vector < _Int > ia_U;
        vector < _Int > ja_U;
        vector < char >jachar_U;
        vector < _Flt > a_U;

        CFct_impl < _Int, _Flt >::SplitLU (false, nlist_ini, *(this->GetIa ()),
                                           *(this->GetJa ()), *(this->GetJaChar ()),
                                           *(this->GetA ()), ia_L, ja_L, jachar_L, a_L,
                                           ia_U, ja_U, jachar_U, a_U);

// Tranpose

        vector < _Int > ia_Lt;
        vector < _Int > ja_Lt;
        vector < char >jachar_Lt;
        vector < _Flt > a_Lt;

        CFct_impl < _Int, _Flt >::Transpose (false, nlist_ini, ia_L, ja_L, jachar_L, a_L,
                                             ia_Lt, ja_Lt, jachar_Lt, a_Lt);

// Combine pairs

        vector < _Int > ia_LU;
        vector < _Int > ja_LU;
        vector < char >jachar_LU;
        vector < _Flt > a_LU;

        CFct_impl < _Int, _Flt >::CombinePairs (false, nlist_ini, ia_Lt, ja_Lt, jachar_Lt,
                                                a_Lt, ia_U, ja_U, jachar_U, a_U, ia_LU,
                                                ja_LU, jachar_LU, a_LU);

// Split pairs

        CFct_impl < _Int, _Flt >::SplitPairsFilter (false, -1.0e0, nlist_ini, ia_LU,
                                                    ja_LU, jachar_LU, a_LU, ia_Lt, ja_Lt,
                                                    jachar_Lt, a_Lt, ia_U, ja_U, jachar_U,
                                                    a_U);

// Extend sparsity by zeroes

        CFct_impl < _Int, _Flt >::Transpose (false, nlist_ini, ia_Lt, ja_Lt, jachar_Lt,
                                             a_Lt, ia_L, ja_L, jachar_L, a_L);

// Combine

        CFct_impl < _Int, _Flt >::CombineLU (false, nlist_ini, ia_L, ja_L, jachar_L, a_L,
                                             ia_U, ja_U, jachar_U, a_U,
                                             *(_asymm.GetIa ()), *(_asymm.GetJa ()),
                                             *(_asymm.GetJaChar ()), *(_asymm.GetA ()));

        vector < _Int > *pia_LU = _asymm.GetIa ();

      int nzja_lu = (int) ((*pia_LU)[nlist_ini]);

        _asymm.SetNlist (nlist_ini);
        _asymm.SetNzja (nzja_lu);
        _asymm.SetNza (nzja_lu);

        _asymm.SetIdentityList ();

   }

/// @brief Split matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::SplitMatrix (int _degree_max, int _isize_max, int &_nparts,
                           vector < long long >&_parts, int *_order)
   {

// Open matrix structure

      int nlistloc = this->GetNlist ();
      _Int *pia_loc = this->GetIaArr ();
      _Int *pja_loc = this->GetJaArr ();

// Mask data

      CVectorData < int >imask (nlistloc);
      CVectorData < int >listloc (nlistloc);

      int *pimask = imask.Ptr ();
      int *plistloc = listloc.Ptr ();

      int i;

      for (i = 0; i < nlistloc; i++) {
         pimask[i] = -1;
      }

// Main cycle

      int nlist_ini = 0;
      int ipart = 0;
      int ibeg_ini = 0;

      int i_curr, nlist_ext, nlist_curr, nlist_ini0, idegree, j, jj, irow, ibeg_ini0;

      while (true) {

// Find new initial point

         i_curr = -1;
         for (i = ibeg_ini; i < nlistloc; i++) {
            if (pimask[i] < 0) {
               i_curr = i;
               break;
            }
         }

         if (i_curr < 0)
            break;

// Self extend list

         nlist_ini0 = nlist_ini;

         nlist_curr = nlist_ini + 1;
         pimask[i_curr] = ipart;
         plistloc[nlist_ini] = i_curr;

         idegree = 0;

         while (nlist_curr > nlist_ini) {
            idegree++;
            nlist_ext = nlist_curr;
            for (i = nlist_ini; i < nlist_curr; i++) {
               irow = plistloc[i];
               for (j = (int) pia_loc[irow]; j < pia_loc[irow + 1]; j++) {
                  jj = (int) pja_loc[j];
                  if (pimask[jj] < 0) {
                     pimask[jj] = ipart;
                     plistloc[nlist_ext] = jj;
                     nlist_ext++;
                  }
               }
            }
            if (idegree >= _degree_max || nlist_ext - nlist_ini0 >= _isize_max) {
               nlist_ini = nlist_ext;
               nlist_curr = nlist_ext;
            } else {
               nlist_ini = nlist_curr;
               nlist_curr = nlist_ext;
            }
         }

// Modify initial search point

         ibeg_ini0 = ibeg_ini;

         for (i = ibeg_ini0; i < nlistloc; i++) {
            if (pimask[i] >= 0) {
               ibeg_ini++;
            } else {
               break;
            }
         }

         ipart++;

      }

// Finally compute order and partitioning

      _nparts = ipart;

      _parts.resize (_nparts + 1);

      long long *p_parts = &_parts[0];

      for (i = 0; i <= _nparts; i++)
         p_parts[i] = 0;

      for (i = 0; i < nlistloc; i++) {
         jj = pimask[i];
         p_parts[jj + 1]++;
      }

      for (i = 0; i < _nparts; i++)
         p_parts[i + 1] += p_parts[i];

      for (i = 0; i < _nparts; i++)
         plistloc[i] = (int) p_parts[i];

      int k;

      for (i = 0; i < nlistloc; i++) {
         jj = pimask[i];
         k = plistloc[jj];
         _order[i] = k;
         plistloc[jj]++;
      }

   }

/// @brief Compute block scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::BlockScaling (int _blksize, double _sclmin, int _nlist, CMatrix < _Int,
                            _Flt > &_sclL, CMatrix < _Int, _Flt > &_sclU, CMatrix < _Int,
                            _Flt > &_sclLInv, CMatrix < _Int, _Flt > &_sclUInv,
                            double &_sclmin_att, double &_sclmax_att, int &_nmodif_scl)
   {

      int blksize_2 = _blksize * _blksize;

// Open matrix structure

      int nlistloc = this->GetNlist ();
      _Int *pia_loc = this->GetIaArr ();
      _Int *pja_loc = this->GetJaArr ();
      _Flt *pa_loc = this->GetAArr ();

      if (nlistloc != _nlist) {
         cout << " BlockScaling: Error: incorrect list !" << endl;
         throw " BlockScaling: Error: incorrect list !";
      }
// Allocate and init block sparsities of scl

      _sclL.ResizeAndSetAll (nlistloc, 0, nlistloc, 0, nlistloc * blksize_2);
      _sclU.ResizeAndSetAll (nlistloc, 0, nlistloc, 0, nlistloc * blksize_2);
      _sclLInv.ResizeAndSetAll (nlistloc, 0, nlistloc, 0, nlistloc * blksize_2);
      _sclUInv.ResizeAndSetAll (nlistloc, 0, nlistloc, 0, nlistloc * blksize_2);

      int i;

      _Int *pitemp;

      pitemp = _sclL.GetListArr ();
      for (i = 0; i < nlistloc; i++)
         pitemp[i] = (_Int) i;
      pitemp = _sclU.GetListArr ();
      for (i = 0; i < nlistloc; i++)
         pitemp[i] = (_Int) i;
      pitemp = _sclLInv.GetListArr ();
      for (i = 0; i < nlistloc; i++)
         pitemp[i] = (_Int) i;
      pitemp = _sclUInv.GetListArr ();
      for (i = 0; i < nlistloc; i++)
         pitemp[i] = (_Int) i;

      pitemp = _sclL.GetIaArr ();
      for (i = 0; i <= nlistloc; i++)
         pitemp[i] = (_Int) i;
      pitemp = _sclU.GetIaArr ();
      for (i = 0; i <= nlistloc; i++)
         pitemp[i] = (_Int) i;
      pitemp = _sclLInv.GetIaArr ();
      for (i = 0; i <= nlistloc; i++)
         pitemp[i] = (_Int) i;
      pitemp = _sclUInv.GetIaArr ();
      for (i = 0; i <= nlistloc; i++)
         pitemp[i] = (_Int) i;

      pitemp = _sclL.GetJaArr ();
      for (i = 0; i < nlistloc; i++)
         pitemp[i] = (_Int) i;
      pitemp = _sclU.GetJaArr ();
      for (i = 0; i < nlistloc; i++)
         pitemp[i] = (_Int) i;
      pitemp = _sclLInv.GetJaArr ();
      for (i = 0; i < nlistloc; i++)
         pitemp[i] = (_Int) i;
      pitemp = _sclUInv.GetJaArr ();
      for (i = 0; i < nlistloc; i++)
         pitemp[i] = (_Int) i;

      _Flt *psclL = _sclL.GetAArr ();
      _Flt *psclU = _sclU.GetAArr ();
      _Flt *psclLInv = _sclLInv.GetAArr ();
      _Flt *psclUInv = _sclUInv.GetAArr ();

      CVectorData < _Flt > work (blksize_2 * 3);
      CVectorData < double >dwork (blksize_2 * 5 + 20 * _blksize);

      _Flt *pwork = work.Ptr ();
      double *pdwork = dwork.Ptr ();

      _sclmin_att = 1.0e100;
      _sclmax_att = -1.0e100;

      double sclmin_temp, sclmax_temp;

      bool b_found;

      int j, jj, nmodif_loc;
      int nmodif_tot = 0;

      for (i = 0; i < nlistloc; i++) {
         b_found = false;
         for (j = (int) pia_loc[i]; j < pia_loc[i + 1]; j++) {
            jj = (int) pja_loc[j];
            if (jj == i) {
               b_found = true;
               CFct_impl < _Int, _Flt >::DenseScaling (_sclmin, _blksize,
                                                       pa_loc + j * blksize_2,
                                                       psclL + i * blksize_2,
                                                       psclU + i * blksize_2,
                                                       psclLInv + i * blksize_2,
                                                       psclUInv + i * blksize_2, pwork,
                                                       pdwork, sclmin_temp, sclmax_temp,
                                                       nmodif_loc);
               if (sclmin_temp < _sclmin_att)
                  _sclmin_att = sclmin_temp;
               if (sclmax_temp > _sclmax_att)
                  _sclmax_att = sclmax_temp;
               nmodif_tot += nmodif_loc;
            }
         }
         if (!b_found) {
            cout << " BlockScaling: Error: diagonal block is not found !" << endl;
            throw " BlockScaling: Error: diagonal block is not found !";
         }
      }

      _nmodif_scl = nmodif_tot;

   }

// Apply ordering to sparse matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OrderMtr (int *_order, CMatrix < _Int, _Flt > &_aord) const
   {

// Reorder the matrix

      int nlistloc = this->GetNlist ();
      int nzjaloc = this->GetNzja ();
      int nzjacharloc = this->GetNzjaChar ();
      const _Int *pia = this->GetIaArr ();
      const _Int *pja = this->GetJaArr ();
      const char *pjachar = this->GetJaCharArr ();
      const _Flt *pa = this->GetAArr ();

      bool is_char = false;
      if (nzjacharloc > 0)
           is_char = true;

        _aord.ResizeAndSetAll (nlistloc, 0, nzjaloc, 0, nzjaloc);

      if (is_char)
      {
         _aord.SetNzjaChar (nzjaloc);
         _aord.ResizeJaChar (nzjaloc);
      }

      _Int *plistord = _aord.GetListArr ();
      _Int *piaord = _aord.GetIaArr ();
      _Int *pjaord = _aord.GetJaArr ();
      char *pjacharord = _aord.GetJaCharArr ();
      _Flt *paord = _aord.GetAArr ();

// Compute inverse order

      vector < int >iord (nlistloc + 1);

      int *piord = &iord[0];

      int i;

      for (i = 0; i < nlistloc; i++)
         piord[_order[i]] = i;

// Reorder the matrix

      for (i = 0; i < nlistloc; i++)
         plistord[i] = (_Int) i;

      for (i = 0; i <= nlistloc; i++)
         piaord[i] = 0;

      int nimax = 0;

      int ni;

      for (i = 0; i < nlistloc; i++) {
         ni = (int) (pia[i + 1] - pia[i]);
         if (ni > nimax)
            nimax = ni;
      }

      vector < CSortInt > iiarr (nimax + 1);
      vector < char >chararr (nimax + 1);
      vector < _Flt > elems (nimax + 1);

      CSortInt *piiarr = &iiarr[0];
      char *pchararr = &chararr[0];
      _Flt *pelems = &elems[0];

      int j;

      for (i = 0; i < nlistloc; i++) {
         j = _order[i];
         piaord[j + 1] = pia[i + 1] - pia[i];
      }

      for (i = 0; i < nlistloc; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, nzloc, ibs, jold, jnew, ind;

      for (inew = 0; inew < nlistloc; inew++) {
         i = piord[inew];
         nzloc = (int) (pia[i + 1] - pia[i]);
         ibs = nz;
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jold = (int) pja[j];
            jnew = _order[jold];
            pjaord[nz] = (_Int) jnew;
            if (is_char)
               pjacharord[nz] = pjachar[j];
            paord[nz] = pa[j];
            nz++;
         }

// Sort elements

         for (j = 0; j < nzloc; j++) {
            piiarr[j].ival = (int) pjaord[ibs + j];
            piiarr[j].i2val = j;
         }

         sort (piiarr, piiarr + nzloc);

         for (j = 0; j < nzloc; j++)
            pelems[j] = paord[ibs + j];
         if (is_char)
            for (j = 0; j < nzloc; j++)
               pchararr[j] = pjacharord[ibs + j];

         for (j = 0; j < nzloc; j++) {
            pjaord[ibs + j] = (_Int) piiarr[j].ival;
            ind = piiarr[j].i2val;
            if (is_char)
               pjacharord[ibs + j] = pchararr[ind];
            paord[ibs + j] = pelems[ind];
         }

      }

   }

// Apply ordering to sparse matrix (sparsity only)
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OrderMtrSp (int *_order, CMatrix < _Int, _Flt > &_aord) const
   {

// Reorder the matrix

      int nlistloc = this->GetNlist ();
      int nzjaloc = this->GetNzja ();
      int nzjacharloc = this->GetNzjaChar ();
      const _Int *pia = this->GetIaArr ();
      const _Int *pja = this->GetJaArr ();
      const char *pjachar = this->GetJaCharArr ();

      bool is_char = false;
      if (nzjacharloc > 0)
           is_char = true;

        _aord.ResizeAndSetAllSp (nlistloc, 0, nzjaloc, 0);

      if (is_char)
      {
         _aord.SetNzjaChar (nzjaloc);
         _aord.ResizeJaChar (nzjaloc);
      }

      _Int *plistord = _aord.GetListArr ();
      _Int *piaord = _aord.GetIaArr ();
      _Int *pjaord = _aord.GetJaArr ();
      char *pjacharord = _aord.GetJaCharArr ();

// Compute inverse order

      vector < int >iord (nlistloc + 1);

      int *piord = &iord[0];

      int i;

      for (i = 0; i < nlistloc; i++)
         piord[_order[i]] = i;

// Reorder the matrix

      for (i = 0; i < nlistloc; i++)
         plistord[i] = (_Int) i;

      for (i = 0; i <= nlistloc; i++)
         piaord[i] = 0;

      int nimax = 0;

      int ni;

      for (i = 0; i < nlistloc; i++) {
         ni = (int) (pia[i + 1] - pia[i]);
         if (ni > nimax)
            nimax = ni;
      }

      vector < CSortInt > iiarr (nimax + 1);
      vector < char >chararr (nimax + 1);

      CSortInt *piiarr = &iiarr[0];
      char *pchararr = &chararr[0];

      int j;

      for (i = 0; i < nlistloc; i++) {
         j = _order[i];
         piaord[j + 1] = pia[i + 1] - pia[i];
      }

      for (i = 0; i < nlistloc; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, nzloc, ibs, jold, jnew, ind;

      for (inew = 0; inew < nlistloc; inew++) {
         i = piord[inew];
         nzloc = (int) (pia[i + 1] - pia[i]);
         ibs = nz;
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jold = (int) pja[j];
            jnew = _order[jold];
            pjaord[nz] = (_Int) jnew;
            if (is_char)
               pjacharord[nz] = pjachar[j];
            nz++;
         }

// Sort elements

         for (j = 0; j < nzloc; j++) {
            piiarr[j].ival = (int) pjaord[ibs + j];
            piiarr[j].i2val = j;
         }

         sort (piiarr, piiarr + nzloc);

         if (is_char)
            for (j = 0; j < nzloc; j++)
               pchararr[j] = pjacharord[ibs + j];

         for (j = 0; j < nzloc; j++) {
            pjaord[ibs + j] = (_Int) piiarr[j].ival;
            ind = piiarr[j].i2val;
            if (is_char)
               pjacharord[ibs + j] = pchararr[ind];
         }

      }

   }

// Apply ordering to sparse matrix and weights data
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OrderMtrWeights (int *_order, CMatrix < _Int, _Flt > &_aordW) const
   {

// Reorder the matrix

      int nlistloc = this->GetNlist ();
      int nzjaloc = this->GetNzja ();
      int nzjacharloc = this->GetNzjaChar ();
      const _Int *plist2 = this->GetList2Arr ();
      const _Int *pia = this->GetIaArr ();
      const _Int *pja = this->GetJaArr ();
      const _Int *pja2 = this->GetJa2Arr ();
      const char *pjachar = this->GetJaCharArr ();

      bool is_char = false;
      if (nzjacharloc > 0)
           is_char = true;

        _aordW.ResizeAndSetAll (nlistloc, nlistloc, nzjaloc, nzjaloc, 0);

      if (is_char)
      {
         _aordW.SetNzjaChar (nzjaloc);
         _aordW.ResizeJaChar (nzjaloc);
      }

      _Int *plistord = _aordW.GetListArr ();
      _Int *plist2ord = _aordW.GetList2Arr ();
      _Int *piaord = _aordW.GetIaArr ();
      _Int *pjaord = _aordW.GetJaArr ();
      _Int *pja2ord = _aordW.GetJa2Arr ();
      char *pjacharord = _aordW.GetJaCharArr ();

// Compute inverse order

      vector < int >iord (nlistloc + 1);

      int *piord = &iord[0];

      int i;

      for (i = 0; i < nlistloc; i++)
         piord[_order[i]] = i;

// Reorder the matrix

      for (i = 0; i < nlistloc; i++)
         plistord[i] = (_Int) i;

      for (i = 0; i <= nlistloc; i++)
         piaord[i] = 0;

      int nimax = 0;

      int ni;

      for (i = 0; i < nlistloc; i++) {
         ni = (int) (pia[i + 1] - pia[i]);
         if (ni > nimax)
            nimax = ni;
      }

      vector < CSortInt > iiarr (nimax + 1);
      vector < _Int > ja2arr (nimax + 1);
      vector < char >chararr (nimax + 1);

      CSortInt *piiarr = &iiarr[0];
      _Int *pja2arr = &ja2arr[0];
      char *pchararr = &chararr[0];

      int j;

      for (i = 0; i < nlistloc; i++) {
         j = _order[i];
         piaord[j + 1] = pia[i + 1] - pia[i];
         plist2ord[j] = plist2[i];
      }

      for (i = 0; i < nlistloc; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, nzloc, ibs, jold, jnew, ind;

      for (inew = 0; inew < nlistloc; inew++) {
         i = piord[inew];
         nzloc = (int) (pia[i + 1] - pia[i]);
         ibs = nz;
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jold = (int) pja[j];
            jnew = _order[jold];
            pjaord[nz] = (_Int) jnew;
            pja2ord[nz] = pja2[j];
            if (is_char)
               pjacharord[nz] = pjachar[j];
            nz++;
         }

// Sort elements

         for (j = 0; j < nzloc; j++) {
            piiarr[j].ival = (int) pjaord[ibs + j];
            piiarr[j].i2val = j;
         }

         sort (piiarr, piiarr + nzloc);

         for (j = 0; j < nzloc; j++) {
            pja2arr[j] = pja2ord[ibs + j];
         }

         if (is_char)
            for (j = 0; j < nzloc; j++)
               pchararr[j] = pjacharord[ibs + j];

         for (j = 0; j < nzloc; j++) {
            pjaord[ibs + j] = (_Int) piiarr[j].ival;
            ind = piiarr[j].i2val;
            pja2ord[ibs + j] = pja2arr[ind];
            if (is_char)
               pjacharord[ibs + j] = pchararr[ind];
         }

      }

   }

// Apply ordering to matrix rows
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OrderMtrRows (int *_order, CMatrix < _Int, _Flt > &_aord) const
   {

// Reorder the matrix

      int nlistloc = this->GetNlist ();
      int nzjaloc = this->GetNzja ();
      int nzjacharloc = this->GetNzjaChar ();
      const _Int *plist = this->GetListArr ();
      const _Int *pia = this->GetIaArr ();
      const _Int *pja = this->GetJaArr ();
      const char *pjachar = this->GetJaCharArr ();
      const _Flt *pa = this->GetAArr ();

      bool is_char = false;
      if (nzjacharloc > 0)
           is_char = true;

        _aord.ResizeAndSetAll (nlistloc, 0, nzjaloc, 0, nzjaloc);

      if (is_char)
      {
         _aord.SetNzjaChar (nzjaloc);
         _aord.ResizeJaChar (nzjaloc);
      }

      _Int *plistord = _aord.GetListArr ();
      _Int *piaord = _aord.GetIaArr ();
      _Int *pjaord = _aord.GetJaArr ();
      char *pjacharord = _aord.GetJaCharArr ();
      _Flt *paord = _aord.GetAArr ();

      vector < CSortInt > iiarr (nlistloc + 1);

      CSortInt *piiarr = &iiarr[0];

      int i;

      for (i = 0; i < nlistloc; i++) {
         piiarr[i].ival = _order[plist[i]];
         piiarr[i].i2val = i;
      }

      sort (piiarr, piiarr + nlistloc);

      vector < int >orderloc (nlistloc + 1);
      vector < int >iorderloc (nlistloc + 1);

      int *porderloc = &orderloc[0];
      int *piorderloc = &iorderloc[0];

      for (i = 0; i < nlistloc; i++)
         plistord[i] = (_Int) piiarr[i].ival;
      for (i = 0; i < nlistloc; i++)
         piorderloc[i] = piiarr[i].i2val;
      for (i = 0; i < nlistloc; i++)
         porderloc[piorderloc[i]] = i;

      for (i = 0; i <= nlistloc; i++)
         piaord[i] = 0;

      int j;

      for (i = 0; i < nlistloc; i++) {
         j = porderloc[i];
         piaord[j + 1] = pia[i + 1] - pia[i];
      }

      for (i = 0; i < nlistloc; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, jold, jnew;

      for (inew = 0; inew < nlistloc; inew++) {
         i = piorderloc[inew];
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jold = (int) pja[j];
            jnew = jold;
            pjaord[nz] = (_Int) jnew;
            if (is_char)
               pjacharord[nz] = pjachar[j];
            paord[nz] = pa[j];
            nz++;
         }
      }

   }

// Apply ordering to matrix rows (sparsity only)
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OrderMtrRowsSp (int *_order, CMatrix < _Int, _Flt > &_aord) const
   {

// Reorder the matrix

      int nlistloc = this->GetNlist ();
      int nzjaloc = this->GetNzja ();
      int nzjacharloc = this->GetNzjaChar ();
      const _Int *plist = this->GetListArr ();
      const _Int *pia = this->GetIaArr ();
      const _Int *pja = this->GetJaArr ();
      const char *pjachar = this->GetJaCharArr ();

      bool is_char = false;
      if (nzjacharloc > 0)
           is_char = true;

        _aord.ResizeAndSetAllSp (nlistloc, 0, nzjaloc, 0);

      if (is_char)
      {
         _aord.SetNzjaChar (nzjaloc);
         _aord.ResizeJaChar (nzjaloc);
      }

      _Int *plistord = _aord.GetListArr ();
      _Int *piaord = _aord.GetIaArr ();
      _Int *pjaord = _aord.GetJaArr ();
      char *pjacharord = _aord.GetJaCharArr ();

      vector < CSortInt > iiarr (nlistloc + 1);

      CSortInt *piiarr = &iiarr[0];

      int i;

      for (i = 0; i < nlistloc; i++) {
         piiarr[i].ival = _order[plist[i]];
         piiarr[i].i2val = i;
      }

      sort (piiarr, piiarr + nlistloc);

      vector < int >orderloc (nlistloc + 1);
      vector < int >iorderloc (nlistloc + 1);

      int *porderloc = &orderloc[0];
      int *piorderloc = &iorderloc[0];

      for (i = 0; i < nlistloc; i++)
         plistord[i] = (_Int) piiarr[i].ival;
      for (i = 0; i < nlistloc; i++)
         piorderloc[i] = piiarr[i].i2val;
      for (i = 0; i < nlistloc; i++)
         porderloc[piorderloc[i]] = i;

      for (i = 0; i <= nlistloc; i++)
         piaord[i] = 0;

      int j;

      for (i = 0; i < nlistloc; i++) {
         j = porderloc[i];
         piaord[j + 1] = pia[i + 1] - pia[i];
      }

      for (i = 0; i < nlistloc; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, jold, jnew;

      for (inew = 0; inew < nlistloc; inew++) {
         i = piorderloc[inew];
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jold = (int) pja[j];
            jnew = jold;
            pjaord[nz] = (_Int) jnew;
            if (is_char)
               pjacharord[nz] = pjachar[j];
            nz++;
         }
      }

   }

// Apply ordering to matrix rows with pairs data
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OrderMtrRowsPairs (int *_order, CMatrix < _Int, _Flt > &_aord) const
   {

// Reorder the matrix

      int nlistloc = this->GetNlist ();
      int nzjaloc = this->GetNzja ();
      int nzjacharloc = this->GetNzjaChar ();
      const _Int *plist = this->GetListArr ();
      const _Int *pia = this->GetIaArr ();
      const _Int *pja = this->GetJaArr ();
      const char *pjachar = this->GetJaCharArr ();
      const _Flt *pa = this->GetAArr ();

      bool is_char = false;
      if (nzjacharloc > 0)
           is_char = true;

        _aord.ResizeAndSetAll (nlistloc, 0, nzjaloc, 0, 2 * nzjaloc);

      if (is_char)
      {
         _aord.SetNzjaChar (nzjaloc);
         _aord.ResizeJaChar (nzjaloc);
      }

      _Int *plistord = _aord.GetListArr ();
      _Int *piaord = _aord.GetIaArr ();
      _Int *pjaord = _aord.GetJaArr ();
      char *pjacharord = _aord.GetJaCharArr ();
      _Flt *paord = _aord.GetAArr ();

      vector < CSortInt > iiarr (nlistloc + 1);

      CSortInt *piiarr = &iiarr[0];

      int i;

      for (i = 0; i < nlistloc; i++) {
         piiarr[i].ival = _order[plist[i]];
         piiarr[i].i2val = i;
      }

      sort (piiarr, piiarr + nlistloc);

      vector < int >orderloc (nlistloc + 1);
      vector < int >iorderloc (nlistloc + 1);

      int *porderloc = &orderloc[0];
      int *piorderloc = &iorderloc[0];

      for (i = 0; i < nlistloc; i++)
         plistord[i] = (_Int) piiarr[i].ival;
      for (i = 0; i < nlistloc; i++)
         piorderloc[i] = piiarr[i].i2val;
      for (i = 0; i < nlistloc; i++)
         porderloc[piorderloc[i]] = i;

      for (i = 0; i <= nlistloc; i++)
         piaord[i] = 0;

      int j;

      for (i = 0; i < nlistloc; i++) {
         j = porderloc[i];
         piaord[j + 1] = pia[i + 1] - pia[i];
      }

      for (i = 0; i < nlistloc; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, jold, jnew;

      for (inew = 0; inew < nlistloc; inew++) {
         i = piorderloc[inew];
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jold = (int) pja[j];
            jnew = jold;
            pjaord[nz] = (_Int) jnew;
            if (is_char)
               pjacharord[nz] = pjachar[j];
            paord[nz * 2] = pa[j * 2];
            paord[nz * 2 + 1] = pa[j * 2 + 1];
            nz++;
         }
      }

   }

// Apply ordering to matrix rows with pairs data
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OrderMtrRowsPairs_BxB (bool _do_pairs, int _blksize, int *_order,
                                     CMatrix < _Int, _Flt > &_aord) const
   {

// Reorder the matrix

      int b_2 = _blksize * _blksize;
      int b_2_2 = b_2;

      if (_do_pairs)
           b_2_2 *= 2;

      int nlistloc = this->GetNlist ();
      int nzjaloc = this->GetNzja ();
      int nzjacharloc = this->GetNzjaChar ();
      const _Int *plist = this->GetListArr ();
      const _Int *pia = this->GetIaArr ();
      const _Int *pja = this->GetJaArr ();
      const char *pjachar = this->GetJaCharArr ();
      const _Flt *pa = this->GetAArr ();

      bool is_char = false;
      if (nzjacharloc > 0)
           is_char = true;

        _aord.ResizeAndSetAll (nlistloc, 0, nzjaloc, 0, 2 * nzjaloc * b_2);

      if (is_char)
      {
         _aord.SetNzjaChar (nzjaloc);
         _aord.ResizeJaChar (nzjaloc);
      }

      _Int *plistord = _aord.GetListArr ();
      _Int *piaord = _aord.GetIaArr ();
      _Int *pjaord = _aord.GetJaArr ();
      char *pjacharord = _aord.GetJaCharArr ();
      _Flt *paord = _aord.GetAArr ();

      vector < CSortInt > iiarr (nlistloc + 1);

      CSortInt *piiarr = &iiarr[0];

      int i;

      for (i = 0; i < nlistloc; i++) {
         piiarr[i].ival = _order[plist[i]];
         piiarr[i].i2val = i;
      }

      sort (piiarr, piiarr + nlistloc);

      vector < int >orderloc (nlistloc + 1);
      vector < int >iorderloc (nlistloc + 1);

      int *porderloc = &orderloc[0];
      int *piorderloc = &iorderloc[0];

      for (i = 0; i < nlistloc; i++)
         plistord[i] = (_Int) piiarr[i].ival;
      for (i = 0; i < nlistloc; i++)
         piorderloc[i] = piiarr[i].i2val;
      for (i = 0; i < nlistloc; i++)
         porderloc[piorderloc[i]] = i;

      for (i = 0; i <= nlistloc; i++)
         piaord[i] = 0;

      int j;

      for (i = 0; i < nlistloc; i++) {
         j = porderloc[i];
         piaord[j + 1] = pia[i + 1] - pia[i];
      }

      for (i = 0; i < nlistloc; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, jold, jnew;

      for (inew = 0; inew < nlistloc; inew++) {
         i = piorderloc[inew];
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jold = (int) pja[j];
            jnew = jold;
            pjaord[nz] = (_Int) jnew;
            if (is_char)
               pjacharord[nz] = pjachar[j];
            CVector < _Flt >::CopyVector (b_2_2, pa + j * b_2_2, paord + nz * b_2_2);
            nz++;
         }
      }

   }

// Apply ordering to matrix cols
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OrderMtrCols (int *_order, CMatrix < _Int, _Flt > &_aord) const
   {

// Reorder the matrix

      int nlistloc = this->GetNlist ();
      int nzjaloc = this->GetNzja ();
      int nzjacharloc = this->GetNzjaChar ();
      const _Int *plist = this->GetListArr ();
      const _Int *pia = this->GetIaArr ();
      const _Int *pja = this->GetJaArr ();
      const char *pjachar = this->GetJaCharArr ();
      const _Flt *pa = this->GetAArr ();

      bool is_char = false;
      if (nzjacharloc > 0)
           is_char = true;

        _aord.ResizeAndSetAll (nlistloc, 0, nzjaloc, 0, nzjaloc);

      if (is_char)
      {
         _aord.SetNzjaChar (nzjaloc);
         _aord.ResizeJaChar (nzjaloc);
      }

      _Int *plistord = _aord.GetListArr ();
      _Int *piaord = _aord.GetIaArr ();
      _Int *pjaord = _aord.GetJaArr ();
      char *pjacharord = _aord.GetJaCharArr ();
      _Flt *paord = _aord.GetAArr ();

      int i;

      for (i = 0; i < nlistloc; i++)
         plistord[i] = plist[i];

      for (i = 0; i <= nlistloc; i++)
         piaord[i] = 0;

      int nimax = 0;

      int ni;

      for (i = 0; i < nlistloc; i++) {
         ni = (int) (pia[i + 1] - pia[i]);
         if (ni > nimax)
            nimax = ni;
      }

      vector < CSortInt > iiarr (nimax + 1);
      vector < _Flt > elems (nimax + 1);
      vector < char >chararr (nimax + 1);

      char *pchararr = &chararr[0];
      _Flt *pelems = &elems[0];
      CSortInt *piiarr = &iiarr[0];

      int j;

      for (i = 0; i < nlistloc; i++) {
         j = i;
         piaord[j + 1] = pia[i + 1] - pia[i];
      }

      for (i = 0; i < nlistloc; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, nzloc, ibs, jold, jnew, ind;

      for (inew = 0; inew < nlistloc; inew++) {
         i = inew;
         nzloc = (int) (pia[i + 1] - pia[i]);
         ibs = nz;
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jold = (int) pja[j];
            jnew = _order[jold];
            pjaord[nz] = (_Int) jnew;
            if (is_char)
               pjacharord[nz] = pjachar[j];
            paord[nz] = pa[j];
            nz++;
         }

// Sort elements

         for (j = 0; j < nzloc; j++) {
            piiarr[j].ival = (int) pjaord[ibs + j];
            piiarr[j].i2val = j;
         }

         sort (piiarr, piiarr + nzloc);

         for (j = 0; j < nzloc; j++)
            pelems[j] = paord[ibs + j];
         if (is_char)
            for (j = 0; j < nzloc; j++)
               pchararr[j] = pjacharord[ibs + j];

         for (j = 0; j < nzloc; j++) {
            pjaord[ibs + j] = (_Int) piiarr[j].ival;
            ind = piiarr[j].i2val;
            if (is_char)
               pjacharord[ibs + j] = pchararr[ind];
            paord[ibs + j] = pelems[ind];
         }

      }

   }

// Apply ordering to matrix cols (sparsity only)
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OrderMtrColsSp (int *_order, CMatrix < _Int, _Flt > &_aord) const
   {

// Reorder the matrix

      int nlistloc = this->GetNlist ();
      int nzjaloc = this->GetNzja ();
      int nzjacharloc = this->GetNzjaChar ();
      const _Int *plist = this->GetListArr ();
      const _Int *pia = this->GetIaArr ();
      const _Int *pja = this->GetJaArr ();
      const char *pjachar = this->GetJaCharArr ();

      bool is_char = false;
      if (nzjacharloc > 0)
           is_char = true;

        _aord.ResizeAndSetAllSp (nlistloc, 0, nzjaloc, 0);

      if (is_char)
      {
         _aord.SetNzjaChar (nzjaloc);
         _aord.ResizeJaChar (nzjaloc);
      }

      _Int *plistord = _aord.GetListArr ();
      _Int *piaord = _aord.GetIaArr ();
      _Int *pjaord = _aord.GetJaArr ();
      char *pjacharord = _aord.GetJaCharArr ();

      int i;

      for (i = 0; i < nlistloc; i++)
         plistord[i] = plist[i];

      for (i = 0; i <= nlistloc; i++)
         piaord[i] = 0;

      int nimax = 0;

      int ni;

      for (i = 0; i < nlistloc; i++) {
         ni = (int) (pia[i + 1] - pia[i]);
         if (ni > nimax)
            nimax = ni;
      }

      vector < CSortInt > iiarr (nimax + 1);
      vector < char >chararr (nimax + 1);

      char *pchararr = &chararr[0];
      CSortInt *piiarr = &iiarr[0];

      int j;

      for (i = 0; i < nlistloc; i++) {
         j = i;
         piaord[j + 1] = pia[i + 1] - pia[i];
      }

      for (i = 0; i < nlistloc; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, nzloc, ibs, jold, jnew, ind;

      for (inew = 0; inew < nlistloc; inew++) {
         i = inew;
         nzloc = (int) (pia[i + 1] - pia[i]);
         ibs = nz;
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jold = (int) pja[j];
            jnew = _order[jold];
            pjaord[nz] = (_Int) jnew;
            if (is_char)
               pjacharord[nz] = pjachar[j];
            nz++;
         }

// Sort elements

         for (j = 0; j < nzloc; j++) {
            piiarr[j].ival = (int) pjaord[ibs + j];
            piiarr[j].i2val = j;
         }

         sort (piiarr, piiarr + nzloc);

         if (is_char)
            for (j = 0; j < nzloc; j++)
               pchararr[j] = pjacharord[ibs + j];

         for (j = 0; j < nzloc; j++) {
            pjaord[ibs + j] = (_Int) piiarr[j].ival;
            ind = piiarr[j].i2val;
            if (is_char)
               pjacharord[ibs + j] = pchararr[ind];
         }

      }

   }

// Apply ordering to matrix cols with pairs data
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OrderMtrColsPairs (int *_order, CMatrix < _Int, _Flt > &_aord) const
   {


// Reorder the matrix

      int nlistloc = this->GetNlist ();
      int nzjaloc = this->GetNzja ();
      int nzjacharloc = this->GetNzjaChar ();
      const _Int *plist = this->GetListArr ();
      const _Int *pia = this->GetIaArr ();
      const _Int *pja = this->GetJaArr ();
      const char *pjachar = this->GetJaCharArr ();
      const _Flt *pa = this->GetAArr ();

      bool is_char = false;
      if (nzjacharloc > 0)
           is_char = true;

        _aord.ResizeAndSetAll (nlistloc, 0, nzjaloc, 0, 2 * nzjaloc);

      if (is_char)
      {
         _aord.SetNzjaChar (nzjaloc);
         _aord.ResizeJaChar (nzjaloc);
      }

      _Int *plistord = _aord.GetListArr ();
      _Int *piaord = _aord.GetIaArr ();
      _Int *pjaord = _aord.GetJaArr ();
      char *pjacharord = _aord.GetJaCharArr ();
      _Flt *paord = _aord.GetAArr ();

      int i;

      for (i = 0; i < nlistloc; i++)
         plistord[i] = plist[i];

      for (i = 0; i <= nlistloc; i++)
         piaord[i] = 0;

      int nimax = 0;

      int ni;

      for (i = 0; i < nlistloc; i++) {
         ni = (int) (pia[i + 1] - pia[i]);
         if (ni > nimax)
            nimax = ni;
      }

      vector < CSortInt > iiarr (nimax + 1);
      vector < _Flt > elems (nimax * 2 + 1);
      vector < char >chararr (nimax + 1);

      char *pchararr = &chararr[0];
      _Flt *pelems = &elems[0];
      CSortInt *piiarr = &iiarr[0];

      int j;

      for (i = 0; i < nlistloc; i++) {
         j = i;
         piaord[j + 1] = pia[i + 1] - pia[i];
      }

      for (i = 0; i < nlistloc; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, nzloc, ibs, jold, jnew, ind;

      for (inew = 0; inew < nlistloc; inew++) {
         i = inew;
         nzloc = (int) (pia[i + 1] - pia[i]);
         ibs = nz;
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jold = (int) pja[j];
            jnew = _order[jold];
            pjaord[nz] = (_Int) jnew;
            if (is_char)
               pjacharord[nz] = pjachar[j];
            paord[nz * 2] = pa[j * 2];
            paord[nz * 2 + 1] = pa[j * 2 + 1];
            nz++;
         }

// Sort elements

         for (j = 0; j < nzloc; j++) {
            piiarr[j].ival = (int) pjaord[ibs + j];
            piiarr[j].i2val = j;
         }

         sort (piiarr, piiarr + nzloc);

         for (j = 0; j < 2 * nzloc; j++)
            pelems[j] = paord[ibs * 2 + j];
         if (is_char)
            for (j = 0; j < nzloc; j++)
               pchararr[j] = pjacharord[ibs + j];

         for (j = 0; j < nzloc; j++) {
            pjaord[ibs + j] = (_Int) piiarr[j].ival;
            ind = piiarr[j].i2val;
            if (is_char)
               pjacharord[ibs + j] = pchararr[ind];
            paord[ibs * 2 + j * 2] = pelems[ind * 2];
            paord[ibs * 2 + j * 2 + 1] = pelems[ind * 2 + 1];
         }

      }

   }

// Apply ordering to matrix cols with pairs data
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OrderMtrColsPairs_BxB (bool _do_pairs, int _blksize, int *_order,
                                     CMatrix < _Int, _Flt > &_aord) const
   {

// Reorder the matrix

      int b_2 = _blksize * _blksize;
      int b_2_2 = b_2;

      if (_do_pairs)
           b_2_2 *= 2;

      int nlistloc = this->GetNlist ();
      int nzjaloc = this->GetNzja ();
      int nzjacharloc = this->GetNzjaChar ();
      const _Int *plist = this->GetListArr ();
      const _Int *pia = this->GetIaArr ();
      const _Int *pja = this->GetJaArr ();
      const char *pjachar = this->GetJaCharArr ();
      const _Flt *pa = this->GetAArr ();

      bool is_char = false;
      if (nzjacharloc > 0)
           is_char = true;

        _aord.ResizeAndSetAll (nlistloc, 0, nzjaloc, 0, 2 * nzjaloc * b_2);

      if (is_char)
      {
         _aord.SetNzjaChar (nzjaloc);
         _aord.ResizeJaChar (nzjaloc);
      }

      _Int *plistord = _aord.GetListArr ();
      _Int *piaord = _aord.GetIaArr ();
      _Int *pjaord = _aord.GetJaArr ();
      char *pjacharord = _aord.GetJaCharArr ();
      _Flt *paord = _aord.GetAArr ();

      int i;

      for (i = 0; i < nlistloc; i++)
         plistord[i] = plist[i];

      for (i = 0; i <= nlistloc; i++)
         piaord[i] = 0;

      int nimax = 0;

      int ni;

      for (i = 0; i < nlistloc; i++) {
         ni = (int) (pia[i + 1] - pia[i]);
         if (ni > nimax)
            nimax = ni;
      }

      vector < CSortInt > iiarr (nimax + 1);
      vector < _Flt > elems (nimax * 2 * b_2 + 1);
      vector < char >chararr (nimax + 1);

      char *pchararr = &chararr[0];
      _Flt *pelems = &elems[0];
      CSortInt *piiarr = &iiarr[0];

      int j;

      for (i = 0; i < nlistloc; i++) {
         j = i;
         piaord[j + 1] = pia[i + 1] - pia[i];
      }

      for (i = 0; i < nlistloc; i++)
         piaord[i + 1] += piaord[i];

      int nz = 0;

      int inew, nzloc, ibs, jold, jnew, ind;

      for (inew = 0; inew < nlistloc; inew++) {
         i = inew;
         nzloc = (int) (pia[i + 1] - pia[i]);
         ibs = nz;
         for (j = (int) pia[i]; j < pia[i + 1]; j++) {
            jold = (int) pja[j];
            jnew = _order[jold];
            pjaord[nz] = (_Int) jnew;
            if (is_char)
               pjacharord[nz] = pjachar[j];
            CVector < _Flt >::CopyVector (b_2_2, pa + j * b_2_2, paord + nz * b_2_2);
            nz++;
         }

// Sort elements

         for (j = 0; j < nzloc; j++) {
            piiarr[j].ival = (int) pjaord[ibs + j];
            piiarr[j].i2val = j;
         }

         sort (piiarr, piiarr + nzloc);

         CVector < _Flt >::CopyVector (b_2_2 * nzloc, paord + ibs * b_2_2, pelems);

         if (is_char)
            for (j = 0; j < nzloc; j++)
               pchararr[j] = pjacharord[ibs + j];

         for (j = 0; j < nzloc; j++) {
            pjaord[ibs + j] = (_Int) piiarr[j].ival;
            ind = piiarr[j].i2val;
            if (is_char)
               pjacharord[ibs + j] = pchararr[ind];
            CVector < _Flt >::CopyVector (b_2_2, pelems + ind * b_2_2,
                                          paord + (ibs + j) * b_2_2);
         }

      }

   }

/// @brief Pack zero rows
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int, _Flt >::PackZeroRows ()
   {

// Check data on entry

      bool b_index2 = false;

      if (this->n_list2 != 0 || this->nz_ja2 != 0) {
         b_index2 = true;
//      throw " CTBlockSp<_TInt>::PackZeroRows: incorrect sparsity on entry ";
      }
// Count the new number of rows

      _Int *plist = this->GetListArr ();
      _Int *plist2 = this->GetList2Arr ();
      _Int *pia = this->GetIaArr ();

      int nlistnew = 0;

      int i;

      for (i = 0; i < this->n_list; i++) {
         if (pia[i + 1] > pia[i])
            nlistnew++;
      }

// Create new arrays

      vector < _Int > listnew (nlistnew + 1);
      vector < _Int > list2new (nlistnew + 1);
      vector < _Int > ianew (nlistnew + 1);

      _Int *plistnew = &listnew[0];
      _Int *plist2new = &list2new[0];
      _Int *pianew = &ianew[0];

      pianew[0] = 0;

      nlistnew = 0;

      for (i = 0; i < this->n_list; i++) {
         if (pia[i + 1] > pia[i]) {
            plistnew[nlistnew] = plist[i];
            if (b_index2) {
               plist2new[nlistnew] = plist2[i];
            }
            pianew[nlistnew + 1] = pia[i + 1];
            nlistnew++;
         }
      }

// Replace and register new arrays

      this->n_list = nlistnew;

      this->list_matr.swap (listnew);
      if (b_index2) {
         this->n_list2 = nlistnew;
         this->list2_matr.swap (list2new);
      }
      this->ia_matr.swap (ianew);

   }

// Compute ND order, create binary tree, find separators and condense tree
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OrderNDSeparatorsForTree (int _nlev, int _ordtype, vector < int >&_order,
                                        CTree & _tree, int &_nblks, vector < int >&_blks)
   {

// Compute ND order

      this->ComputeOptimalOrder (-1, _order);

      int *porder = &_order[0];

// Apply order to sparsity

      CMatrix < _Int, _Flt > a_ord_sp;

      this->OrderMtrSp (porder, a_ord_sp);

// Symmetrize 

      CMatrix < _Int, _Flt > a_ord_symm;

      a_ord_sp.SymmetrizeSparsitySp (a_ord_symm);

      a_ord_sp.Clean ();

// Create initial binary tree

      int nnodes_base = 1;

      int i;

      for (i = 0; i < _nlev - 1; i++)
         nnodes_base *= 2;

      CTree tree (nnodes_base, 2);

      int nnodes_tree = tree.GetNnodes ();
      int *pnchilds = tree.GetNchilds ();
      vector < int >*pchilds_list = tree.GetChildsList ();
      int *pnode2ind = tree.GetNode2Ind ();
      int *pind2node = tree.GetInd2Node ();

      for (i = 0; i < nnodes_tree; i++) {
         pnode2ind[i] = -1;
         pind2node[i] = -1;
      }

// Find separators and filter blocks partitionings

      a_ord_symm.FindSeparatorsForTree (tree, _nblks, _blks);

      if (_nblks != nnodes_tree) {
         cout <<
            " CMatrix<>::OrderNDSeparatorsForTree: error: incorrect number of blocks ! "
            << endl;
         throw
            " CMatrix<>::OrderNDSeparatorsForTree: error: incorrect number of blocks ! ";
      }

      int *pblks = &_blks[0];

      vector < int >blks_flt (_nblks + 1);
      int *pblks_flt = &blks_flt[0];

      int nblks_flt = 0;

      pblks_flt[0] = 0;

      for (i = 0; i < _nblks; i++) {
         if (pblks[i + 1] > pblks[i]) {
            pnode2ind[i] = nblks_flt;
            pind2node[nblks_flt] = i;
            pblks_flt[nblks_flt + 1] = pblks[i + 1];
            nblks_flt++;
         } else {
            pnode2ind[i] = -1;
         }
      }

      _nblks = nblks_flt;
      _blks.swap (blks_flt);

      pblks = &_blks[0];

// Filter the tree

      CVectorData < int >imasknd (nnodes_tree);
      int *pimasknd = imasknd.Ptr ();

      for (i = 0; i < nnodes_tree; i++)
         pimasknd[i] = -1;

      int nnodes_flt = 0;

      for (i = 0; i < nnodes_tree; i++) {
         if (pnode2ind[i] != -1) {
            pimasknd[i] = nnodes_flt;
            nnodes_flt++;
         } else {
            int j, ichild;
            int nchilds_loc = pnchilds[i];
            int *ppchilds_list = &pchilds_list[i][0];
            if (nchilds_loc > 1 || ppchilds_list[0] != i) {
               bool b_found = false;
               for (j = 0; j < nchilds_loc; j++) {
                  ichild = ppchilds_list[j];
                  if (pimasknd[ichild] >= 0) {
                     b_found = true;
                  }
               }
               if (b_found) {
                  pimasknd[i] = nnodes_flt;
                  nnodes_flt++;
               }
            }
         }
      }

      if (nnodes_flt != nnodes_tree) {
         CTree tree_flt;
         tree.FilterTree (pimasknd, tree_flt);
         tree = tree_flt;
      }

      _tree = tree;

// Optimize profile of diagonal subblocks

      int nsupA = pblks[_nblks];

      {
         int ibeg, iend, niloc, j;
         CVectorData < int >order_new (nsupA);
         int *porder_new = order_new.Ptr ();
         for (i = 0; i < _nblks; i++) {
            ibeg = pblks[i];
            iend = pblks[i + 1] - 1;
            niloc = iend - ibeg + 1;
            CMatrix < _Int, _Flt > a_sub;
            a_ord_symm.GetMainSubmatrixSp (ibeg, iend, a_sub);
            int nzja_sub = a_sub.GetNzja ();
            _Int *plist_sub = a_sub.GetListArr ();
            _Int *pja_sub = a_sub.GetJaArr ();
            for (j = 0; j < niloc; j++)
               plist_sub[j] -= ibeg;
            for (j = 0; j < nzja_sub; j++)
               pja_sub[j] -= ibeg;
            vector < int >ord_sub;
            a_sub.ComputeOptimalOrder (_ordtype, ord_sub);
            int *pord_sub = &ord_sub[0];
            for (j = 0; j < niloc; j++)
               porder_new[ibeg + j] = ibeg + pord_sub[j];
         }
         CVectorData < int >iorder (nsupA);
         int *piorder = iorder.Ptr ();
         for (i = 0; i < nsupA; i++)
            piorder[porder[i]] = i;
         int iold;
         for (i = 0; i < nsupA; i++) {
            iold = piorder[i];
            porder[iold] = porder_new[i];
         }
      }

   }

// Perform multiplication of matrices stored by rows, resulting matrix is also stored by rows
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::MatrixByMatrixMultiplySp (int &_icycle, int *_imask, int *_imask1,
                                        int *_indarr, int *_listloc, const CMatrix < _Int,
                                        _Flt > &_a, const CMatrix < _Int, _Flt > &_b,
                                        CMatrix < _Int, _Flt > &_a_times_b)
   {

// Open sparsities of A and B

      int nlist_a = _a.GetNlist ();
      const _Int *plist_a = _a.GetListArr ();
      const _Int *pia_a = _a.GetIaArr ();
      const _Int *pja_a = _a.GetJaArr ();

      int nlist_b = _b.GetNlist ();
      const _Int *plist_b = _b.GetListArr ();
      const _Int *pia_b = _b.GetIaArr ();
      const _Int *pja_b = _b.GetJaArr ();

// Perform multiplication

      vector < _Int > list_ab (nlist_a + 1);
      vector < _Int > ia_ab (nlist_a + 1);

      _Int *plist_ab = &list_ab[0];
      _Int *pia_ab = &ia_ab[0];

      vector < _Int > ja_ab (1);

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
         ja_ab.resize (nzja_ab + nlistloc + 1);
         for (j = 0; j < nlistloc; j++)
            ja_ab[nzja_ab + j] = (_Int) _listloc[j];
         nzja_ab += nlistloc;
         pia_ab[i + 1] = nzja_ab;
      }

// Store result

      vector < _Int > *p_list_ab = _a_times_b.GetList ();
      vector < _Int > *p_ia_ab = _a_times_b.GetIa ();
      vector < _Int > *p_ja_ab = _a_times_b.GetJa ();

      p_list_ab->swap (list_ab);
      p_ia_ab->swap (ia_ab);
      p_ja_ab->swap (ja_ab);

      _a_times_b.SetNlist (nlist_a);
      _a_times_b.SetNzja (nzja_ab);
      _a_times_b.SetNza (0);

   }

// Perform multiplication of matrices stored by rows, resulting matrix is also stored by rows
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::MatrixByMatrixMultiply (int &_icycle, int *_imask, int *_imask1,
                                      int *_indarr, int *_listloc, _Flt * _fmask,
                                      const CMatrix < _Int, _Flt > &_a,
                                      const CMatrix < _Int, _Flt > &_b, CMatrix < _Int,
                                      _Flt > &_a_times_b)
   {

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

      _Flt fzero = (_Flt) 0.0e0;

      int nzja_ab = 0;

      int j, k, kk, nlistloc, ind;
      _Flt aux;

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
            _fmask[jj] = fzero;
         }
         for (j = (int) pia_a[i]; j < pia_a[i + 1]; j++) {
            jj = (int) pja_a[j];
            if (_imask1[jj] == icycle1) {
               aux = pa_a[j];
               ind = _indarr[jj];
               for (k = (int) pia_b[ind]; k < pia_b[ind + 1]; k++) {
                  kk = (int) pja_b[k];
                  _fmask[kk] += aux * pa_b[k];
               }
            }
         }
         ja_ab.resize (nzja_ab + nlistloc + 1);
         a_ab.resize (nzja_ab + nlistloc + 1);
         for (j = 0; j < nlistloc; j++)
            ja_ab[nzja_ab + j] = (_Int) _listloc[j];
         for (j = 0; j < nlistloc; j++) {
            jj = _listloc[j];
            a_ab[nzja_ab + j] = _fmask[jj];
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
      _a_times_b.SetNza (nzja_ab);

   }

/// @brief Compute matrix decomposition
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::DecompWeights (bool _split_unconnected, CMatrix < _Int, _Flt > &_amatr_strW,
                             int _nparts, int *_partition)
   {

// Compute transposed matrix

      int ntot = _amatr_strW.GetNlist ();
      int nztot = _amatr_strW.GetNzja ();

      _Int *plisttot_w = _amatr_strW.GetList2Arr ();
      _Int *piatot_sp = _amatr_strW.GetIaArr ();
      _Int *pjatot_sp = _amatr_strW.GetJaArr ();
      _Int *pjatot_w = _amatr_strW.GetJa2Arr ();

      CMatrix < _Int, _Flt > at_sp;

      at_sp.ResizeAndSetAllSp (ntot, ntot, nztot, nztot);

      _Int *plist_t_sp = at_sp.GetListArr ();
      _Int *pia_t_sp = at_sp.GetIaArr ();
      _Int *pja_t_sp = at_sp.GetJaArr ();
      _Int *pja_t_w = at_sp.GetJa2Arr ();

      int i, j, jj;

      for (i = 0; i <= ntot; i++)
         pia_t_sp[i] = 0;

      for (i = 0; i < ntot; i++) {
         for (j = (int) piatot_sp[i]; j < piatot_sp[i + 1]; j++) {
            jj = (int) pjatot_sp[j];
            pia_t_sp[jj + 1]++;
         }
      }

      for (i = 0; i < ntot; i++)
         pia_t_sp[i + 1] = pia_t_sp[i] + pia_t_sp[i + 1];
      for (i = 0; i < ntot; i++)
         plist_t_sp[i] = pia_t_sp[i];

      int k;

      for (i = 0; i < ntot; i++) {
         for (j = (int) piatot_sp[i]; j < piatot_sp[i + 1]; j++) {
            jj = (int) pjatot_sp[j];
            k = (int) plist_t_sp[jj];
            pja_t_sp[k] = (int) i;
            pja_t_w[k] = pjatot_w[j];
            plist_t_sp[jj]++;
         }
      }

// Create filtered symmetrized matrix (sparsity and weights)

#ifdef USE_METIS
      typedef idx_t TINT_TYPE;
      typedef real_t TFLT_TYPE;
#else
      typedef int TINT_TYPE;
      typedef double TFLT_TYPE;
#endif

      vector < TINT_TYPE > listsymm_w (ntot + 1);
      vector < TINT_TYPE > iasymm (ntot + 1);
      vector < TINT_TYPE > jasymm (2 * nztot + 1);
      vector < TINT_TYPE > jasymm_w (2 * nztot + 1);

      TINT_TYPE *plistsymm_w = &listsymm_w[0];
      TINT_TYPE *piasymm = &iasymm[0];
      TINT_TYPE *pjasymm = &jasymm[0];
      TINT_TYPE *pjasymm_w = &jasymm_w[0];

      for (i = 0; i < ntot; i++)
         plistsymm_w[i] = (TINT_TYPE) plisttot_w[i];

      piasymm[0] = 0;

      int nznew = 0;

      int ibeg, iend, ibegt, iendt, ip, ipt, jjt;

      for (i = 0; i < ntot; i++) {
         ibeg = (int) piatot_sp[i];
         iend = (int) piatot_sp[i + 1] - 1;
         ibegt = (int) pia_t_sp[i];
         iendt = (int) pia_t_sp[i + 1] - 1;
         ip = ibeg;
         ipt = ibegt;
         while (ip <= iend || ipt <= iendt) {
            if (ip <= iend && ipt <= iendt) {
               jj = (int) pjatot_sp[ip];
               jjt = (int) pja_t_sp[ipt];
               if (jj == jjt) {
                  if (jj != i) {
                     pjasymm[nznew] = (TINT_TYPE) jj;
                     pjasymm_w[nznew] = (TINT_TYPE) ((pjatot_w[ip] + pja_t_w[ipt]) / 2);
                     nznew++;
                  }
                  ip++;
                  ipt++;
               } else if (jj < jjt) {
                  if (jj != i) {
                     pjasymm[nznew] = (TINT_TYPE) jj;
                     pjasymm_w[nznew] = (TINT_TYPE) pjatot_w[ip];
                     nznew++;
                  }
                  ip++;
               } else if (jj > jjt) {
                  if (jjt != i) {
                     pjasymm[nznew] = (TINT_TYPE) jjt;
                     pjasymm_w[nznew] = (TINT_TYPE) pja_t_w[ipt];
                     nznew++;
                  }
                  ipt++;
               }
            } else if (ip <= iend) {
               jj = (int) pjatot_sp[ip];
               if (jj != i) {
                  pjasymm[nznew] = (TINT_TYPE) jj;
                  pjasymm_w[nznew] = (TINT_TYPE) pjatot_w[ip];
                  nznew++;
               }
               ip++;
            } else if (ipt <= iendt) {
               jjt = (int) pja_t_sp[ipt];
               if (jjt != i) {
                  pjasymm[nznew] = (TINT_TYPE) jjt;
                  pjasymm_w[nznew] = (TINT_TYPE) pja_t_w[ipt];
                  nznew++;
               }
               ipt++;
            }
         }
         piasymm[i + 1] = (TINT_TYPE) nznew;
      }

      at_sp.Clean ();

// Call Metis

      TINT_TYPE ncon = 1;
      TINT_TYPE nparts = _nparts;

      CVectorData < TFLT_TYPE > tpwgts ((int) (ncon * nparts));
      CVectorData < TFLT_TYPE > ubvec ((int) ncon);

      TFLT_TYPE *ptpwgts = tpwgts.Ptr ();
      TFLT_TYPE *pubvec = ubvec.Ptr ();

      double aux = 1.0e0 / ((double) nparts);
      TFLT_TYPE faux = (TFLT_TYPE) aux;

      for (i = 0; i < nparts; i++)
         ptpwgts[i] = faux;

      pubvec[0] = (TFLT_TYPE) 1.05e0;

      if (nparts == 1) {

         for (i = 0; i < ntot; i++)
            _partition[i] = 0;

      } else {

         TINT_TYPE ntot_int = ntot;

         if (ntot_int > 0) {

            CVectorData < TINT_TYPE > partition (ntot);

            TINT_TYPE *ppartition = partition.Ptr ();

//         if (sizeof(idx_t) != sizeof(long long)) throw "  Error: incompartible data types ! ";

#ifdef USE_METIS

            idx_t options[METIS_NOPTIONS];

            METIS_SetDefaultOptions (options);

            options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
            options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
            options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;

            options[METIS_OPTION_NUMBERING] = 0;
            options[METIS_OPTION_SEED] = 0;

            idx_t edgecut;

            METIS_PartGraphRecursive ((idx_t *) & ntot_int, (idx_t *) & ncon,
                                      (idx_t *) piasymm, (idx_t *) pjasymm,
                                      (idx_t *) plistsymm_w, NULL, (idx_t *) pjasymm_w,
                                      (idx_t *) & nparts, ptpwgts, pubvec,
                                      (idx_t *) options, (idx_t *) & edgecut,
                                      (idx_t *) ppartition);

#else
            cout <<
               " CMatrix<>::DecompWeights: error: Metis PartGraphRecursive is not found !!! "
               << endl;
            throw
               " CMatrix<>::DecompWeights: error: Metis PartGraphRecursive is not found !!! ";
#endif

            for (i = 0; i < ntot; i++)
               _partition[i] = (int) ppartition[i];

         }

      }

// Count the number of parts

      if (_split_unconnected) {

         int nparts_tot = 0;

         for (i = 0; i < ntot; i++) {
            if (_partition[i] > nparts_tot)
               nparts_tot = _partition[i];
         }

         nparts_tot++;

// Create lists of parts

         CVectorData < int >ibs_part (nparts_tot + 1);
         CVectorData < int >iptr_part (nparts_tot);
         CVectorData < int >list_part (ntot);

         int *pibs_part = ibs_part.Ptr ();
         int *piptr_part = iptr_part.Ptr ();
         int *plist_part = list_part.Ptr ();

         for (i = 0; i <= nparts_tot; i++)
            pibs_part[i] = 0;

         for (i = 0; i < ntot; i++) {
            j = _partition[i];
            pibs_part[j + 1]++;
         }

         for (i = 0; i < nparts_tot; i++)
            pibs_part[i + 1] = pibs_part[i] + pibs_part[i + 1];

         for (i = 0; i < nparts_tot; i++)
            piptr_part[i] = pibs_part[i];

         for (i = 0; i < ntot; i++) {
            j = _partition[i];
            k = piptr_part[j];
            plist_part[k] = (int) i;
            piptr_part[j]++;
         }

// Compute splitting of parts into connected subparts

         CVectorData < int >part_new (ntot);
         CVectorData < int >list_new (ntot);

         int *ppart_new = part_new.Ptr ();
         int *plist_new = list_new.Ptr ();

         for (i = 0; i < ntot; i++)
            ppart_new[i] = -1;

         int nparts_new = 0;
         int nlist_new = 0;

         int ipart;

         int ind, nlist_new_ini, kk;

         for (ipart = 0; ipart < nparts_tot; ipart++) {
            while (true) {
               ind = -1;
               for (i = pibs_part[ipart]; i < pibs_part[ipart + 1]; i++) {
                  j = plist_part[i];
                  if (ppart_new[j] == -1) {
                     ind = j;
                     break;
                  }
               }
               if (ind < 0)
                  break;
               nlist_new_ini = nlist_new;
               plist_new[nlist_new] = (int) ind;
               ppart_new[ind] = nparts_new;
               nlist_new++;
               for (i = nlist_new_ini; i < nlist_new; i++) {
                  j = plist_new[i];
                  for (k = (int) piasymm[j]; k < piasymm[j + 1]; k++) {
                     kk = (int) pjasymm[k];
                     if (_partition[kk] == ipart && ppart_new[kk] == -1) {
                        plist_new[nlist_new] = (int) kk;
                        ppart_new[kk] = nparts_new;
                        nlist_new++;
                     }
                  }
               }
               nparts_new++;
            }
         }

         for (i = 0; i < ntot; i++) {
            _partition[i] = ppart_new[i];
         }

      }

   }

/// @brief Compute 2 level matrix decomposition
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::DecompWeights2Level (int _nblks_max, int _nparts2blk_max, CMatrix < _Int,
                                   _Flt > &_amatr_strW, int *_order, int &_nblks,
                                   vector < int >&_blk2parts, vector < int >&_parts)
   {

// Perform initial global partitioning

      int ntot = _amatr_strW.GetNlist ();

      CVectorData < int >partition (ntot);

      int *ppartition = partition.Ptr ();

      int nblks_decomp = _nblks_max;
      if (nblks_decomp > ntot)
         nblks_decomp = ntot;

      if (nblks_decomp > 1) {
         CMatrix < _Int, _Flt >::DecompWeights (false, _amatr_strW, nblks_decomp,
                                                ppartition);
      } else {
         int i;
         for (i = 0; i < ntot; i++)
            ppartition[i] = 0;
      }

// Initial order and partitioning

      int nblks_curr = 0;

      int i, ind;

      for (i = 0; i < ntot; i++) {
         ind = ppartition[i];
         if (ind > nblks_curr)
            nblks_curr = ind;
      }

      nblks_curr++;

      CVectorData < int >blks_curr (nblks_curr + 1);
      int *pblks_curr = blks_curr.Ptr ();

      for (i = 0; i <= nblks_curr; i++)
         pblks_curr[i] = 0;

      for (i = 0; i < ntot; i++) {
         ind = ppartition[i];
         pblks_curr[ind + 1]++;
      }

      for (i = 0; i < nblks_curr; i++)
         pblks_curr[i + 1] += pblks_curr[i];

      CVectorData < int >iptr (ntot + 1);
      int *piptr = iptr.Ptr ();

      for (i = 0; i < nblks_curr; i++)
         piptr[i] = pblks_curr[i];

      int k;

      for (i = 0; i < ntot; i++) {
         ind = ppartition[i];
         k = piptr[ind];
         _order[i] = k;
         piptr[ind]++;
      }

      int nblks_flt = 0;

      for (i = 0; i < nblks_curr; i++) {
         if (pblks_curr[i + 1] > pblks_curr[i]) {
            pblks_curr[nblks_flt + 1] = pblks_curr[i + 1];
            nblks_flt++;
         }
      }

      nblks_curr = nblks_flt;

      _nblks = nblks_flt;

// Reorder explicitely matrix and weights data

      CMatrix < _Int, _Flt > aordW;

      _amatr_strW.OrderMtrWeights (_order, aordW);

// Perform diagonal blocks computations

      _Int *plist2_ord = aordW.GetList2Arr ();
      _Int *pia_ord = aordW.GetIaArr ();
      _Int *pja_ord = aordW.GetJaArr ();
      _Int *pja2_ord = aordW.GetJa2Arr ();

      CVectorData < int >nparts_arr (nblks_curr + 1);
      vector < CVectorData < int > >parts_arr (nblks_curr + 1);
      vector < CVectorData < int > >order_arr (nblks_curr + 1);

      int *pnparts_arr = nparts_arr.Ptr ();
      CVectorData < int >*pparts_arr = &parts_arr[0];
      CVectorData < int >*porder_arr = &order_arr[0];

      int ibegblk, iendblk, niloc, j, jrow, kcol, nzloc;

      for (i = 0; i < nblks_curr; i++) {
         ibegblk = pblks_curr[i];
         iendblk = pblks_curr[i + 1] - 1;
         niloc = iendblk + 1 - ibegblk;
         nzloc = 0;
         for (j = 0; j < niloc; j++) {
            jrow = ibegblk + j;
            for (k = (int) pia_ord[jrow]; k < pia_ord[jrow + 1]; k++) {
               kcol = (int) pja_ord[k];
               if (kcol >= ibegblk && kcol <= iendblk)
                  nzloc++;
            }
         }
         CMatrix < _Int, _Flt > asubW;
         asubW.ResizeAndSetAll (niloc, niloc, nzloc, nzloc, 0);
         _Int *plist_sub = asubW.GetListArr ();
         _Int *plist2_sub = asubW.GetList2Arr ();
         _Int *pia_sub = asubW.GetIaArr ();
         _Int *pja_sub = asubW.GetJaArr ();
         _Int *pja2_sub = asubW.GetJa2Arr ();
         for (j = 0; j < niloc; j++) {
            jrow = ibegblk + j;
            plist_sub[j] = j;
            plist2_sub[j] = plist2_ord[jrow];
         }
         pia_sub[0] = 0;
         nzloc = 0;
         for (j = 0; j < niloc; j++) {
            jrow = ibegblk + j;
            for (k = (int) pia_ord[jrow]; k < pia_ord[jrow + 1]; k++) {
               kcol = (int) pja_ord[k];
               if (kcol >= ibegblk && kcol <= iendblk) {
                  pja_sub[nzloc] = kcol - ibegblk;
                  pja2_sub[nzloc] = pja2_ord[k];
                  nzloc++;
               }
            }
            pia_sub[j + 1] = nzloc;
         }
         porder_arr[i].resize (niloc);
         int *porder_loc = porder_arr[i].Ptr ();
         CVectorData < int >partition_sub (niloc);
         int *ppartition_sub = partition_sub.Ptr ();
         int nparts_decomp = _nparts2blk_max;
         if (nparts_decomp > niloc)
            nparts_decomp = niloc;
         if (nparts_decomp > 1) {
            CMatrix < _Int, _Flt >::DecompWeights (false, asubW, nparts_decomp,
                                                   ppartition_sub);
         } else {
            for (j = 0; j < niloc; j++)
               ppartition_sub[j] = 0;
         }
         int nparts_curr = 0;
         for (j = 0; j < niloc; j++) {
            ind = ppartition_sub[j];
            if (ind > nparts_curr)
               nparts_curr = ind;
         }
         nparts_curr++;
         pparts_arr[i].resize (nparts_curr + 1);
         int *pparts_curr = pparts_arr[i].Ptr ();
         for (j = 0; j <= nparts_curr; j++)
            pparts_curr[j] = 0;
         for (j = 0; j < niloc; j++) {
            ind = ppartition_sub[j];
            pparts_curr[ind + 1]++;
         }
         for (j = 0; j < nparts_curr; j++)
            pparts_curr[j + 1] += pparts_curr[j];
         for (j = 0; j < nparts_curr; j++)
            piptr[j] = pparts_curr[j];
         for (j = 0; j < niloc; j++) {
            ind = ppartition_sub[j];
            k = piptr[ind];
            porder_loc[j] = k;
            piptr[ind]++;
         }
         int nparts_flt = 0;

         for (j = 0; j < nparts_curr; j++) {
            if (pparts_curr[j + 1] > pparts_curr[j]) {
               pparts_curr[nparts_flt + 1] = pparts_curr[j + 1];
               nparts_flt++;
            }
         }
         nparts_curr = nparts_flt;
         pnparts_arr[i] = nparts_curr;
      }

// Prepare final data

      _blk2parts.resize (nblks_curr + 1);
      int *p_blk2parts = &_blk2parts[0];

      CVectorData < int >iorder (ntot);
      CVectorData < int >order2 (ntot);

      int *piorder = iorder.Ptr ();
      int *porder2 = order2.Ptr ();

      p_blk2parts[0] = 0;
      int nparts_tot = 0;

      for (i = 0; i < nblks_curr; i++) {
         ibegblk = pblks_curr[i];
         iendblk = pblks_curr[i + 1] - 1;
         niloc = iendblk + 1 - ibegblk;
         int *porder_loc = porder_arr[i].Ptr ();
         for (j = 0; j < niloc; j++) {
            porder2[ibegblk + j] = ibegblk + porder_loc[j];
         }
         nparts_tot += pnparts_arr[i];
         p_blk2parts[i + 1] = nparts_tot;
      }

      _parts.resize (nparts_tot + 1);
      int *p_parts = &_parts[0];

      nparts_tot = 0;

      p_parts[0] = 0;

      for (i = 0; i < nblks_curr; i++) {
         int *pparts_loc = pparts_arr[i].Ptr ();
         for (j = 0; j < pnparts_arr[i]; j++) {
            p_parts[nparts_tot + 1] =
               p_parts[nparts_tot] + (pparts_loc[j + 1] - pparts_loc[j]);
            nparts_tot++;
         }
      }

      for (i = 0; i < ntot; i++)
         piorder[_order[i]] = i;

      for (i = 0; i < ntot; i++) {
         ind = piorder[i];
         _order[ind] = porder2[i];
      }
/*
   if (true) {
      _amatr_strW.OrderMtrWeights (_order, aordW);
      CBMatrix<_Int,_Flt>::Str2PsBox2 (aordW, "StrW_ord.ps", _nblks, p_blk2parts, p_parts);
   }
*/
   }

/// @brief Perform in-place explicit block scaling
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::ExplicitBlockScale (int _blksize, CMatrix < _Int, _Flt > &_sclL,
                                  int *_ref_jcol, CMatrix < _Int, _Flt > &_sclU,
                                  _Flt * _work)
   {

      int blksize_2 = _blksize * _blksize;

// Open matrix structure

      int nlistloc = this->GetNlist ();
      _Int *plist_loc = this->GetListArr ();
      _Int *pia_loc = this->GetIaArr ();
      _Int *pja_loc = this->GetJaArr ();
      _Flt *pa_loc = this->GetAArr ();

      _Flt *psclL = _sclL.GetAArr ();
      _Flt *psclU = _sclU.GetAArr ();

// Scale

      int i, j, jj, irow, indU;

      for (i = 0; i < nlistloc; i++) {
         irow = (int) plist_loc[i];
         for (j = (int) pia_loc[i]; j < pia_loc[i + 1]; j++) {
            jj = (int) pja_loc[j];
            if (_ref_jcol == NULL) {
               indU = jj;
            } else {
               indU = _ref_jcol[jj];
            }
            CVector < _Flt >::MMM (_blksize, psclL + irow * blksize_2,
                                   pa_loc + j * blksize_2, _work);
            CVector < _Flt >::MMM (_blksize, _work, psclU + indU * blksize_2,
                                   pa_loc + j * blksize_2);
         }
      }

   }

/// @brief Compute comparison matrix
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::ComparisonMatrix_BxB (char _diatype, int _blksize, const CMatrix < _Int,
                                    _Flt > &_amatr_bxb, CMatrix < _Int, _Flt > &_amatr_pt)
   {

      int blksize_2 = _blksize * _blksize;

// Open matrix structure

      int nlistloc = _amatr_bxb.GetNlist ();
      int nzjaloc = _amatr_bxb.GetNzja ();
      int nzaloc = _amatr_bxb.GetNza ();
      const _Int *plist_loc = _amatr_bxb.GetListArr ();
      const _Int *pia_loc = _amatr_bxb.GetIaArr ();
      const _Int *pja_loc = _amatr_bxb.GetJaArr ();
      const _Flt *pa_loc = _amatr_bxb.GetAArr ();

      if (nzaloc != nzjaloc * blksize_2) {
         cout << " CMatrix<>::ComparisonMatrix_BxB: incorrect matrix on entry ! " << endl;
         throw " CMatrix<>::ComparisonMatrix_BxB: incorrect matrix on entry ! ";
      }

      _amatr_pt.GetSparsity (_amatr_bxb);

      _amatr_pt.ResizeA (nzjaloc);
      _amatr_pt.SetNza (nzjaloc);

      _Flt *pa_pt_loc = _amatr_pt.GetAArr ();

      int i, j, irow, icol, ibs, kk;
      _Flt aux, diamin, offdmax;

      for (i = 0; i < nlistloc; i++) {
         irow = (int) plist_loc[i];
         for (j = (int) pia_loc[i]; j < pia_loc[i + 1]; j++) {
            icol = (int) pja_loc[j];
            ibs = j * blksize_2;
            if (irow == icol && _diatype == 'D') {
               aux = pa_loc[ibs];
               if (aux < 0.0e0)
                  aux = -aux;
               diamin = aux;
               for (kk = 1; kk < _blksize; kk++) {
                  aux = pa_loc[ibs + kk * _blksize + kk];
                  if (aux < 0.0e0)
                     aux = -aux;
                  if (diamin > aux)
                     diamin = aux;
               }
               pa_pt_loc[j] = diamin;
            } else {
               aux = pa_loc[ibs];
               if (aux < 0.0e0)
                  aux = -aux;
               offdmax = aux;
               for (kk = 1; kk < blksize_2; kk++) {
                  aux = pa_loc[ibs + kk];
                  if (aux < 0.0e0)
                     aux = -aux;
                  if (offdmax < aux)
                     offdmax = aux;
               }
               pa_pt_loc[j] = offdmax;
            }
         }
      }

   }

//
// Print matrix data
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::PrintMatrix (ofstream & _fout)
   {

      _fout << " CMatrix:" << endl;

      _fout << "    Nlist = " << this->GetNlist () << " Nlist2 = " << this->
         GetNlist2 () << endl;
      _fout << "    Nzja  = " << this->GetNzja () << " Nzja2  = " << this->
         GetNzja2 () << " Nzjachar  = " << this->GetNzjaChar () << " Nza = " << this->
         GetNza () << endl;

      if (this->GetNlist () > 0)
         PrintArray (_fout, "    List  = ", this->GetNlist (), this->GetListArr ());
      if (this->GetNlist2 () > 0)
         PrintArray (_fout, "    List2 = ", this->GetNlist2 (), this->GetList2Arr ());
      if (this->GetNlist () > 0)
         PrintArray (_fout, "    Ia    = ", this->GetNlist () + 1, this->GetIaArr ());
      if (this->GetNzja () > 0)
         PrintArray (_fout, "    Ja    = ", this->GetNzja (), this->GetJaArr ());
      if (this->GetNzja2 () > 0)
         PrintArray (_fout, "    Ja2   = ", this->GetNzja2 (), this->GetJa2Arr ());
      if (this->GetNzjaChar () > 0)
         PrintArray (_fout, "    JaChar = ", this->GetNzjaChar (), this->GetJaCharArr ());
      if (this->GetNza () > 0)
         PrintArray (_fout, "    A     = ", this->GetNza (), this->GetAArr ());

   }

//
// Print matrix data by rows
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::PrintMatrixRows (ofstream & _fout, int _blksize)
   {

      _fout << " CMatrix:" << endl;

      _fout << "    Nlist = " << this->GetNlist () << " Nlist2 = " << this->
         GetNlist2 () << endl;
      _fout << "    Nzja  = " << this->GetNzja () << " Nzja2  = " << this->
         GetNzja2 () << " Nzjachar  = " << this->GetNzjaChar () << " Nza = " << this->
         GetNza () << endl;

      int nlistloc = this->GetNlist ();
      _Int *plist = this->GetListArr ();
      _Int *pia = this->GetIaArr ();
      _Int *pja = this->GetJaArr ();
      _Flt *pa = this->GetAArr ();

      int i, ni, ibeg, irow;

      for (i = 0; i < nlistloc; i++) {
         irow = (int) plist[i];
         _fout << " Ilist = " << i << " Irow = " << irow << endl;
         ibeg = (int) pia[i];
         ni = (int) (pia[i + 1] - pia[i]);
         PrintArray (_fout, " Colmns = ", ni, pja + ibeg);
         PrintArray (_fout, " Elems = ", ni * _blksize * _blksize,
                     pa + ibeg * _blksize * _blksize);
      }

   }

#ifdef TURN_ON_EX
#pragma GCC push_options
#pragma GCC optimize("O0")
#endif

// Compute the packed size
//========================================================================================
   template < typename _Int, typename _Flt > int CMatrix < _Int, _Flt >::GetPackedSize () const
   {

      int nlistloc = this->n_list;
      int nlist2loc = this->n_list2;
      int nzjaloc = this->nz_ja;
      int nzja2loc = this->nz_ja2;
      int nzaloc = this->nz_a;

      int isize =
         5 * sizeof (int) + (2 * nlistloc + 1 + nlist2loc + nzjaloc +
                             nzja2loc) * sizeof (_Int) + nzaloc * sizeof (_Flt);

        return isize;

   }

// Fill by packed data
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::FillPacked (int _length, char *_obj) const
   {

      int nlistloc = this->n_list;
      int nlist2loc = this->n_list2;
      int nzjaloc = this->nz_ja;
      int nzja2loc = this->nz_ja2;
      int nzaloc = this->nz_a;

      char *pLoc;

        pLoc = _obj;

      int *pHead = NULL;
      _Int *plist_obj = NULL;
      _Int *plist2_obj = NULL;
      _Int *pia_obj = NULL;
      _Int *pja_obj = NULL;
      _Int *pja2_obj = NULL;
      _Flt *pa_obj = NULL;

        pHead = (int *) pLoc;
        pLoc += 5 * sizeof (int);

        plist_obj = (_Int *) pLoc;
        pLoc += nlistloc * sizeof (_Int);

        plist2_obj = (_Int *) pLoc;
        pLoc += nlist2loc * sizeof (_Int);

        pia_obj = (_Int *) pLoc;
        pLoc += (nlistloc + 1) * sizeof (_Int);

        pja_obj = (_Int *) pLoc;
        pLoc += (nzjaloc) * sizeof (_Int);

        pja2_obj = (_Int *) pLoc;
        pLoc += (nzja2loc) * sizeof (_Int);

        pa_obj = (_Flt *) pLoc;
        pLoc += (nzaloc) * sizeof (_Flt);

        pHead[0] = nlistloc;
        pHead[1] = nlist2loc;
        pHead[2] = nzjaloc;
        pHead[3] = nzja2loc;
        pHead[4] = nzaloc;

      if (pLoc - _obj != _length)
      {
         throw " CMatrix<_Int,_Flt>::FillPacked: incorrect length on entry ";
      }
// Fill arrays

      const _Int *plist = &(this->list_matr[0]);
      const _Int *plist2 = &(this->list2_matr[0]);
      const _Int *pia = &(this->ia_matr[0]);
      const _Int *pja = &(this->ja_matr[0]);
      const _Int *pja2 = &(this->ja2_matr[0]);
      const _Flt *pa = &(this->a_matr[0]);

      int j;

      for (j = 0; j < nlistloc; j++)
         plist_obj[j] = plist[j];
      for (j = 0; j < nlist2loc; j++)
         plist2_obj[j] = plist2[j];
      for (j = 0; j < nlistloc + 1; j++)
         pia_obj[j] = pia[j];
      for (j = 0; j < nzjaloc; j++)
         pja_obj[j] = pja[j];
      for (j = 0; j < nzja2loc; j++)
         pja2_obj[j] = pja2[j];
      for (j = 0; j < nzaloc; j++)
         pa_obj[j] = pa[j];

   }

// Fill by packed data
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::UnPack (int _length, char *_obj)
   {

// Get head data

      char *pLoc;

      pLoc = _obj;

      int *pHead;

      pHead = (int *) pLoc;
      pLoc += 5 * sizeof (int);

      int nlistloc = pHead[0];
      int nlist2loc = pHead[1];
      int nzjaloc = pHead[2];
      int nzja2loc = pHead[3];
      int nzaloc = pHead[4];

      _Int *plist_obj = NULL;
      _Int *plist2_obj = NULL;
      _Int *pia_obj = NULL;
      _Int *pja_obj = NULL;
      _Int *pja2_obj = NULL;
      _Flt *pa_obj = NULL;

      plist_obj = (_Int *) pLoc;
      pLoc += nlistloc * sizeof (_Int);

      plist2_obj = (_Int *) pLoc;
      pLoc += nlist2loc * sizeof (_Int);

      pia_obj = (_Int *) pLoc;
      pLoc += (nlistloc + 1) * sizeof (_Int);

      pja_obj = (_Int *) pLoc;
      pLoc += (nzjaloc) * sizeof (_Int);

      pja2_obj = (_Int *) pLoc;
      pLoc += (nzja2loc) * sizeof (_Int);

      pa_obj = (_Flt *) pLoc;
      pLoc += (nzaloc) * sizeof (_Flt);

// Store data

      this->n_list = nlistloc;
      this->n_list2 = nlist2loc;
      this->nz_ja = nzjaloc;
      this->nz_ja2 = nzja2loc;
      this->nz_a = nzaloc;

      this->list_matr.resize (nlistloc + 1);
      this->list2_matr.resize (nlist2loc + 1);
      this->ia_matr.resize (nlistloc + 1);
      this->ja_matr.resize (nzjaloc + 1);
      this->ja2_matr.resize (nzja2loc + 1);
      this->a_matr.resize (nzaloc + 1);

      int i;
      _Int *piarr;
      _Flt *paarr;

      piarr = &(this->list_matr[0]);
      for (i = 0; i < nlistloc; i++)
         piarr[i] = plist_obj[i];

      piarr = &(this->list2_matr[0]);
      for (i = 0; i < nlist2loc; i++)
         piarr[i] = plist2_obj[i];

      piarr = &(this->ia_matr[0]);
      for (i = 0; i <= nlistloc; i++)
         piarr[i] = pia_obj[i];

      piarr = &(this->ja_matr[0]);
      for (i = 0; i < nzjaloc; i++)
         piarr[i] = pja_obj[i];

      piarr = &(this->ja2_matr[0]);
      for (i = 0; i < nzja2loc; i++)
         piarr[i] = pja2_obj[i];

      paarr = &(this->a_matr[0]);
      for (i = 0; i < nzaloc; i++)
         paarr[i] = pa_obj[i];

   }

#ifdef TURN_ON_EX
#pragma GCC pop_options
#endif

// Output header data
//========================================================================================
   template < typename _Int, typename _Flt > void CMatrix < _Int,
      _Flt >::OutputHead (ostream & _fout)
   {
      _fout << " CMatrix head data:" << endl;
      _fout << "     Nlist = " << this->n_list << " Nlist2 = " << this->n_list2 << endl;
      _fout << "     Nzja = " << this->nz_ja << " Nzja2 = " << this->nz_ja2 << endl;
      _fout << "     Nzja_char = " << this->nz_jachar << " Nza = " << this->nz_a << endl;
      _fout << "   Arrays sizes:" << endl;
      _fout << "     List = " << this->list_matr.size () << " List2 = " << this->
         list2_matr.size () << endl;
      _fout << "     Ja = " << this->ja_matr.size () << " Ja2 = " << this->ja_matr.
         size () << endl;
      _fout << "     Ja_char = " << this->jachar_matr.size () << " A = " << this->a_matr.
         size () << endl;

   }

//
// Set vector data by zeroes
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::SetByZeroes (int _n,
                                                                          _FltVect * _x)
   {

      int i;

      _FltVect fzero = (_FltVect) 0.0e0;

      for (i = 0; i < _n; i++)
         _x[i] = fzero;

   }

//
// Set vector data by zeroes (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::SetByZeroes_thr (int _n,
                                                                              _FltVect *
                                                                              _x)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::SetByZeroes (ni_loc, _x + ibeg);

         }

      }

   }

//
// Set vector data by ones
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::SetByOnes (int _n,
                                                                        _FltVect * _x)
   {

      int i;

      _FltVect fone = (_FltVect) 1.0e0;

      for (i = 0; i < _n; i++)
         _x[i] = fone;

   }

//
// Set vector data by zeroes (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::SetByOnes_thr (int _n,
                                                                            _FltVect * _x)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::SetByOnes (ni_loc, _x + ibeg);

         }

      }

   }

//
// Copy vector
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::CopyVector (int _n,
                                                                         const _FltVect *
                                                                         _x,
                                                                         _FltVect * _y)
   {

      int i;

      for (i = 0; i < _n; i++)
         _y[i] = _x[i];

   }

//
// Copy vector (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::CopyVector_thr (int _n,
                                                                             const
                                                                             _FltVect *
                                                                             _x,
                                                                             _FltVect *
                                                                             _y)
   {

//      CVector<_FltVect>::CopyVector (_n, _x, _y);

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::CopyVector (ni_loc, _x + ibeg, _y + ibeg);

         }

      }

   }

//
// Compute scalar product
//========================================================================================
   template < typename _FltVect > _FltVect CVector < _FltVect >::ScProd (int _n,
                                                                         const _FltVect *
                                                                         _x,
                                                                         const _FltVect *
                                                                         _y)
   {

      _FltVect fsum = (_FltVect) 0.0e0;

      int i;

      for (i = 0; i < _n; i++)
         fsum += _x[i] * _y[i];

      return fsum;

   }

//
// Compute scalar product (threads version)
//========================================================================================
   template < typename _FltVect > _FltVect CVector < _FltVect >::ScProd_thr (int _n,
                                                                             const
                                                                             _FltVect *
                                                                             _x,
                                                                             const
                                                                             _FltVect *
                                                                             _y) {

      _FltVect fsum = (_FltVect) 0.0e0;

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      vector < _FltVect > fsum_arr (n_thr + 1);
      _FltVect *pfsum_arr = &fsum_arr[0];

      CVector < _FltVect >::SetByZeroes (n_thr, pfsum_arr);

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               pfsum_arr[ipar] =
                  CVector < _FltVect >::ScProd (ni_loc, _x + ibeg, _y + ibeg);

         }
      }

      int i;

      for (i = 0; i < n_thr; i++)
         fsum += pfsum_arr[i];

      return fsum;

   }

//
// Add vector data
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::AddVector (int _n,
                                                                        const _FltVect *
                                                                        _x1,
                                                                        const _FltVect *
                                                                        _x2,
                                                                        _FltVect *
                                                                        _x1plus2)
   {

      int i;

      for (i = 0; i < _n; i++)
         _x1plus2[i] = _x1[i] + _x2[i];

   }

//
// Add vector data (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::AddVector_thr (int _n,
                                                                            const _FltVect
                                                                            * _x1,
                                                                            const _FltVect
                                                                            * _x2,
                                                                            _FltVect *
                                                                            _x1plus2)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::AddVector (ni_loc, _x1 + ibeg, _x2 + ibeg,
                                                _x1plus2 + ibeg);

         }

      }

   }

//
// Add vector data
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::AddReplaceVector (int _n,
                                                                               const
                                                                               _FltVect *
                                                                               _x1,
                                                                               _FltVect *
                                                                               _x1plus2)
   {

      int i;

      for (i = 0; i < _n; i++)
         _x1plus2[i] += _x1[i];

   }

//
// Add vector data (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::AddReplaceVector_thr (int _n,
                                                                                   const
                                                                                   _FltVect
                                                                                   * _x1,
                                                                                   _FltVect
                                                                                   *
                                                                                   _x1plus2)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::AddReplaceVector (ni_loc, _x1 + ibeg,
                                                       _x1plus2 + ibeg);

         }

      }

   }

//
// Subtract vector data
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::SubtractReplaceVector (int _n, const _FltVect * _x1,
                                         _FltVect * _x2minus1)
   {

      int i;

      for (i = 0; i < _n; i++)
         _x2minus1[i] -= _x1[i];

   }

//
// Subtract vector data (threads version)
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::SubtractReplaceVector_thr (int _n, const _FltVect * _x1,
                                             _FltVect * _x2minus1)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::SubtractReplaceVector (ni_loc, _x1 + ibeg,
                                                            _x2minus1 + ibeg);

         }

      }

   }

// Update array (axpy)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::UpdateVector (int _n,
                                                                           const _FltVect
                                                                           * _value,
                                                                           const _FltVect
                                                                           * _arr_x,
                                                                           _FltVect *
                                                                           _arr_y)
   {

      _FltVect value = *_value;

      for (int i = 0; i < _n; i++)
         _arr_y[i] += value * _arr_x[i];

   }

//
// Update array (axpy) (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::UpdateVector_thr (int _n,
                                                                               const
                                                                               _FltVect *
                                                                               _value,
                                                                               const
                                                                               _FltVect *
                                                                               _arr_x,
                                                                               _FltVect *
                                                                               _arr_y)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::UpdateVector (ni_loc, _value, _arr_x + ibeg,
                                                   _arr_y + ibeg);

         }

      }

   }

// Update array for minus alpha (axpy)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::UpdateVectorMinus (int _n,
                                                                                const
                                                                                _FltVect *
                                                                                _value,
                                                                                const
                                                                                _FltVect *
                                                                                _arr_x,
                                                                                _FltVect *
                                                                                _arr_y)
   {

      _FltVect value = *_value;

      int i;

      for (i = 0; i < _n; i++)
         _arr_y[i] = _arr_y[i] - value * _arr_x[i];

   }

//
// Update array for minus alpha (axpy) (threads version)
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::UpdateVectorMinus_thr (int _n, const _FltVect * _value,
                                         const _FltVect * _arr_x, _FltVect * _arr_y)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::UpdateVectorMinus (ni_loc, _value, _arr_x + ibeg,
                                                        _arr_y + ibeg);

         }

      }

   }

// Update array reversed (aypx)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::UpdateVectorReversed (int _n,
                                                                                   const
                                                                                   _FltVect
                                                                                   *
                                                                                   _value,
                                                                                   const
                                                                                   _FltVect
                                                                                   *
                                                                                   _arr_x,
                                                                                   _FltVect
                                                                                   *
                                                                                   _arr_y)
   {

      _FltVect value = *_value;

      for (int i = 0; i < _n; i++)
         _arr_y[i] = value * _arr_y[i] + _arr_x[i];

   }

//
// Update array reversed (aypx) (threads version)
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::UpdateVectorReversed_thr (int _n, const _FltVect * _value,
                                            const _FltVect * _arr_x, _FltVect * _arr_y)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::UpdateVectorReversed (ni_loc, _value, _arr_x + ibeg,
                                                           _arr_y + ibeg);

         }

      }

   }


// Multiply array by value y=a*y
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::MultiplyVectorValue (int _n,
                                                                                  const
                                                                                  _FltVect
                                                                                  *
                                                                                  _value,
                                                                                  _FltVect
                                                                                  *
                                                                                  _arr_y)
   {

      _FltVect value = *_value;

      for (int i = 0; i < _n; i++)
         _arr_y[i] *= value;

   }

//
// Multiply array by value y=a*y (threads version)
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::MultiplyVectorValue_thr (int _n, const _FltVect * _value,
                                           _FltVect * _arr_y)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::MultiplyVectorValue (ni_loc, _value, _arr_y + ibeg);

         }

      }

   }

//
// Compute vector inverse
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::InverseVector (int _n,
                                                                            _FltVect * _x)
   {

      int i;

      _FltVect fone = (_FltVect) 1.0;

      for (i = 0; i < _n; i++)
         _x[i] = fone / _x[i];

   }

//
// Compute vector inverse (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::InverseVector_thr (int _n,
                                                                                _FltVect *
                                                                                _x)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::InverseVector (ni_loc, _x + ibeg);

         }

      }

   }

///
/// @brief Compute block dot (C = A^t*B)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::BlockDot (int _m, int _n,
                                                                       const _FltVect *
                                                                       _amatr,
                                                                       const _FltVect *
                                                                       _bmatr,
                                                                       _FltVect * _cmatr)
   {

      _FltVect fzero;

      CVector < _FltVect >::SetByZeroes (1, &fzero);

      int i, j, k;

      _FltVect aux;

      for (i = 0; i < _n; i++) {
         for (j = 0; j < _n; j++) {
            aux = fzero;
            for (k = 0; k < _m; k++) {
               aux += _amatr[k * _n + i] * _bmatr[k * _n + j];
            }
            _cmatr[i * _n + j] = aux;
         }
      }
   }

//
// Compute block dot (C = A^t*B) (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::BlockDot_thr (int _m, int _n,
                                                                           const _FltVect
                                                                           * _amatr,
                                                                           const _FltVect
                                                                           * _bmatr,
                                                                           _FltVect *
                                                                           _cmatr,
                                                                           _FltVect *
                                                                           _work)
   {

      int n_2 = _n * _n;

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _m / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _m - 1;

            ni_loc = iend + 1 - ibeg;

            CVector < _FltVect >::BlockDot (ni_loc, _n, _amatr + ibeg * _n,
                                            _bmatr + ibeg * _n, _work + ipar * n_2);

         }

         CVector < _FltVect >::SetByZeroes (n_2, _cmatr);

         int i;

         for (i = 0; i < n_thr; i++)
            CVector < _FltVect >::AddReplaceVector (n_2, _work + i * n_2, _cmatr);

      }

   }

///
/// @brief Compute block dot (C = A^t*B)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::BlockDot (int _m, int _na,
                                                                       int _nb,
                                                                       const _FltVect *
                                                                       _amatr,
                                                                       const _FltVect *
                                                                       _bmatr,
                                                                       _FltVect * _cmatr)
   {

      _FltVect fzero;

      CVector < _FltVect >::SetByZeroes (1, &fzero);

      int i, j, k;

      _FltVect aux;

      for (i = 0; i < _na; i++) {
         for (j = 0; j < _nb; j++) {
            aux = fzero;
            for (k = 0; k < _m; k++) {
               aux += _amatr[k * _na + i] * _bmatr[k * _nb + j];
            }
            _cmatr[i * _nb + j] = aux;
         }
      }
   }

//
// Compute block dot (C = A^t*B) (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::BlockDot_thr (int _m,
                                                                           int _na,
                                                                           int _nb,
                                                                           const _FltVect
                                                                           * _amatr,
                                                                           const _FltVect
                                                                           * _bmatr,
                                                                           _FltVect *
                                                                           _cmatr,
                                                                           _FltVect *
                                                                           _work)
   {

      int n_2 = _na * _nb;

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _m / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _m - 1;

            ni_loc = iend + 1 - ibeg;

            CVector < _FltVect >::BlockDot (ni_loc, _na, _nb, _amatr + ibeg * _na,
                                            _bmatr + ibeg * _nb, _work + ipar * n_2);

         }

         CVector < _FltVect >::SetByZeroes (n_2, _cmatr);

         int i;

         for (i = 0; i < n_thr; i++)
            CVector < _FltVect >::AddReplaceVector (n_2, _work + i * n_2, _cmatr);

      }

   }

///
/// @brief Compute block daxpy (C += A*B)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::BlockDaxpy (int _m, int _n,
                                                                         const _FltVect *
                                                                         _amatr,
                                                                         const _FltVect *
                                                                         _bmatr,
                                                                         _FltVect *
                                                                         _cmatr)
   {

      _FltVect fzero;

      CVector < _FltVect >::SetByZeroes (1, &fzero);

      int i, j, k;

      _FltVect aux;

      for (i = 0; i < _m; i++) {
         for (j = 0; j < _n; j++) {
            aux = fzero;
            for (k = 0; k < _n; k++) {
               aux += _amatr[i * _n + k] * _bmatr[k * _n + j];
            }
            _cmatr[i * _n + j] += aux;
         }
      }

   }

//
// Compute block daxpy (C += A*B) (C = A^t*B) (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::BlockDaxpy_thr (int _m,
                                                                             int _n,
                                                                             const
                                                                             _FltVect *
                                                                             _amatr,
                                                                             const
                                                                             _FltVect *
                                                                             _bmatr,
                                                                             _FltVect *
                                                                             _cmatr)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _m / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _m - 1;

            ni_loc = iend + 1 - ibeg;

            CVector < _FltVect >::BlockDaxpy (ni_loc, _n, _amatr + ibeg * _n, _bmatr,
                                              _cmatr + ibeg * _n);

         }

      }

   }

///
/// @brief Compute block daxpy (C += A*B)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::BlockDaxpy (int _m, int _na,
                                                                         int _nb,
                                                                         const _FltVect *
                                                                         _amatr,
                                                                         const _FltVect *
                                                                         _bmatr,
                                                                         _FltVect *
                                                                         _cmatr)
   {

      _FltVect fzero;

      CVector < _FltVect >::SetByZeroes (1, &fzero);

      int i, j, k;

      _FltVect aux;

      for (i = 0; i < _m; i++) {
         for (j = 0; j < _nb; j++) {
            aux = fzero;
            for (k = 0; k < _na; k++) {
               aux += _amatr[i * _na + k] * _bmatr[k * _nb + j];
            }
            _cmatr[i * _nb + j] += aux;
         }
      }

   }

//
// Compute block daxpy (C += A*B) (C = A^t*B) (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::BlockDaxpy_thr (int _m,
                                                                             int _na,
                                                                             int _nb,
                                                                             const
                                                                             _FltVect *
                                                                             _amatr,
                                                                             const
                                                                             _FltVect *
                                                                             _bmatr,
                                                                             _FltVect *
                                                                             _cmatr)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _m / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _m - 1;

            ni_loc = iend + 1 - ibeg;

            CVector < _FltVect >::BlockDaxpy (ni_loc, _na, _nb, _amatr + ibeg * _na,
                                              _bmatr, _cmatr + ibeg * _nb);

         }

      }

   }

///
/// @brief Compute block Housholder transformation
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::BlockHousholder (bool _b_use_thr, int _m, int _n, _FltVect * _q,
                                   _FltVect * _rdiag, _FltVect * _tau, _FltVect * _work)
   {

// Split work memory

      _FltVect *pwork = _work;

      _FltVect *pqtemp = pwork;
      pwork += _n * _m;

      _FltVect *ptau = pwork;
      pwork += _n;

      _FltVect *pscprod = pwork;
      pwork += _n;

      _FltVect *pvect = pwork;
      pwork += _n;

      _FltVect *pqtq_matr = pwork;
      pwork += _n * _n;

// Transpose block

      int i, j;

      for (i = 0; i < _m; i++) {
         for (j = 0; j < _n; j++) {
            pqtemp[j * _m + i] = _q[i * _n + j];
         }
      }

// Compute QR decomposition

      if (_b_use_thr) {
         CVector < _FltVect >::QrdBlock_thr (_n, _m, pqtemp, _m, ptau);
      } else {
         CVector < _FltVect >::QrdBlock (_n, _m, pqtemp, _m, ptau);
      }

// Store R part

      CVector < _FltVect >::SetByZeroes (_n * _n, _rdiag);

      for (i = 0; i < _n; i++) {
         for (j = i; j < _n; j++) {
            _rdiag[i * _n + j] = pqtemp[j * _m + i];
         }
      }

// Store Q part

      for (i = 0; i < _m; i++) {
         for (j = 0; j < _n; j++) {
            _q[i * _n + j] = pqtemp[j * _m + i];
         }
      }

      for (i = 0; i < _n; i++) {
         CVector < _FltVect >::SetByZeroes (_n - i, _q + i * _n + i);
         CVector < _FltVect >::SetByOnes (1, _q + i * _n + i);
      }

// Compute QtQ

      if (_b_use_thr) {
         CVector < _FltVect >::BlockDot_thr (_m, _n, _q, _q, pqtq_matr, pwork);
      } else {
         CVector < _FltVect >::BlockDot (_m, _n, _q, _q, pqtq_matr);
      }

// Compute Tau

      CVector < _FltVect >::SetByZeroes (_n * _n, _tau);

      for (i = 0; i < _n; i++) {

// Init diagonal

         _tau[i * _n + i] = ptau[i];

// Scalar products

         for (j = 0; j < i; j++) {
            pscprod[j] = pqtq_matr[j * _n + i];
         }

// Multiply

         CVector < _FltVect >::Mvm ('T', i, _tau, _n, pscprod, pvect);

// Store

         for (j = 0; j < i; j++) {
            _tau[j * _n + i] = -pvect[j] * ptau[i];
         }

      }

// Scale tau array

      _FltVect mone = (_FltVect) - 1.0e0;

      CVector < _FltVect >::MultiplyVectorValue (_n * _n, &mone, _tau);

   }

//
// Order vector data
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::OrderVector (int _n,
                                                                          const int
                                                                          *_order,
                                                                          const _FltVect *
                                                                          _x,
                                                                          _FltVect *
                                                                          _x_ord)
   {

      int i;

      for (i = 0; i < _n; i++)
         _x_ord[_order[i]] = _x[i];

   }

//
// Order vector data (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::OrderVector_thr (int _n,
                                                                              const int
                                                                              *_order,
                                                                              const
                                                                              _FltVect *
                                                                              _x,
                                                                              _FltVect *
                                                                              _x_ord)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::OrderVector (ni_loc, _order + ibeg, _x + ibeg,
                                                  _x_ord);

         }

      }

   }

//
// Order vector data
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::OrderVector (int _blksize,
                                                                          int _n,
                                                                          const int
                                                                          *_order,
                                                                          const _FltVect *
                                                                          _x,
                                                                          _FltVect *
                                                                          _x_ord)
   {

      int i;

      for (i = 0; i < _n; i++) {
         CVector < _FltVect >::CopyVector (_blksize, _x + i * _blksize,
                                           _x_ord + _order[i] * _blksize);
      }

   }

//
// Order vector data (threads version)
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::OrderVector_thr (int _blksize, int _n, const int *_order,
                                   const _FltVect * _x, _FltVect * _x_ord)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::OrderVector (_blksize, ni_loc, _order + ibeg,
                                                  _x + ibeg * _blksize, _x_ord);

         }

      }

   }

//
// Inverse order vector data
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::InvOrderVector (int _n,
                                                                             const int
                                                                             *_order,
                                                                             const
                                                                             _FltVect *
                                                                             _x,
                                                                             _FltVect *
                                                                             _x_ord)
   {

      int i;

      for (i = 0; i < _n; i++)
         _x_ord[i] = _x[_order[i]];

   }

//
// Inverse order vector data (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::InvOrderVector_thr (int _n,
                                                                                 const int
                                                                                 *_order,
                                                                                 const
                                                                                 _FltVect
                                                                                 * _x,
                                                                                 _FltVect
                                                                                 * _x_ord)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::InvOrderVector (ni_loc, _order + ibeg, _x,
                                                     _x_ord + ibeg);

         }

      }

   }

//
// Inverse order vector data
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::InvOrderVector (int _blksize,
                                                                             int _n,
                                                                             const int
                                                                             *_order,
                                                                             const
                                                                             _FltVect *
                                                                             _x,
                                                                             _FltVect *
                                                                             _x_ord)
   {

      int i;

      for (i = 0; i < _n; i++) {
         CVector < _FltVect >::CopyVector (_blksize, _x + _order[i] * _blksize,
                                           _x_ord + i * _blksize);
      }

   }

//
// Inverse order vector data (threads version)
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::InvOrderVector_thr (int _blksize, int _n, const int *_order,
                                      const _FltVect * _x, _FltVect * _x_ord)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVector < _FltVect >::InvOrderVector (_blksize, ni_loc, _order + ibeg, _x,
                                                     _x_ord + ibeg * _blksize);

         }

      }

   }

/// @brief Compute Housholder transformation
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::Housholder (int _m,
                                                                         _FltVect &
                                                                         _alpha,
                                                                         _FltVect * _x,
                                                                         _FltVect & _tau)
   {

// Compute X norm

      _FltVect fzero;
      _FltVect fone;

      CVector < _FltVect >::SetByZeroes (1, &fzero);
      CVector < _FltVect >::SetByOnes (1, &fone);

      _FltVect xnorm_2 = fzero;

      xnorm_2 = CVector < _FltVect >::ScProd (_m - 1, _x, _x);

      if (xnorm_2 == fzero) {

         _tau = fzero;

      } else {

         _FltVect rnorm = _alpha * _alpha + xnorm_2;

         rnorm = sqrt (rnorm);

         _FltVect beta;

         if (_alpha >= fzero) {
            beta = -rnorm;
         } else {
            beta = rnorm;
         }

         _FltVect diff = _alpha - beta;

         _tau = -diff / beta;

         _FltVect gamma = fone / diff;

         CVector < _FltVect >::MultiplyVectorValue (_m - 1, &gamma, _x);

         _alpha = beta;

      }

   }

/// @brief Compute Housholder transformation (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::Housholder_thr (int _m,
                                                                             _FltVect &
                                                                             _alpha,
                                                                             _FltVect *
                                                                             _x,
                                                                             _FltVect &
                                                                             _tau)
   {

// Compute X norm

      _FltVect fzero;
      _FltVect fone;

      CVector < _FltVect >::SetByZeroes (1, &fzero);
      CVector < _FltVect >::SetByOnes (1, &fone);

      _FltVect xnorm_2 = fzero;

      xnorm_2 = CVector < _FltVect >::ScProd_thr (_m - 1, _x, _x);

      if (xnorm_2 == fzero) {

         _tau = fzero;

      } else {

         _FltVect rnorm = _alpha * _alpha + xnorm_2;

         rnorm = sqrt (rnorm);

         _FltVect beta;

         if (_alpha >= fzero) {
            beta = -rnorm;
         } else {
            beta = rnorm;
         }

         _FltVect diff = _alpha - beta;

         _tau = -diff / beta;

         _FltVect gamma = fone / diff;

         CVector < _FltVect >::MultiplyVectorValue_thr (_m - 1, &gamma, _x);

         _alpha = beta;

      }

   }

/// @brief Compute QR decomposition for the current block
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::QrdBlock (int _ncol,
                                                                       int _nrow,
                                                                       _FltVect * _qblk,
                                                                       int _ldq,
                                                                       _FltVect * _tau)
   {

      int i, j, k;
      _FltVect scprod;

      _FltVect *pa, *ph;

// Main cycle over columns

      for (i = 0; i < _ncol; i++) {

// Apply previous columns to the current column

         for (j = 0; j < i; j++) {

            pa = _qblk + i * _ldq + j + 1;
            ph = _qblk + j * _ldq + j + 1;

            scprod = CVector < _FltVect >::ScProd (_nrow - j - 1, pa, ph);

            scprod += _qblk[i * _ldq + j];

            scprod *= _tau[j];

            _qblk[i * _ldq + j] -= scprod;

            pa = _qblk + i * _ldq + j + 1;
            ph = _qblk + j * _ldq + j + 1;

            _FltVect scprod_minus = -scprod;

            CVector < _FltVect >::UpdateVector (_nrow - j - 1, &scprod_minus, ph, pa);

         }

// Compute new transformation

         j = _nrow - 1;
         if (i + 1 < _nrow - 1)
            j = i + 1;

         k = _nrow - i;

         CVector < _FltVect >::Housholder (k, _qblk[i * _ldq + i], _qblk + i * _ldq + j,
                                           _tau[i]);

      }

   }

/// @brief Compute QR decomposition for the current block (threads version)
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::QrdBlock_thr (int _ncol,
                                                                           int _nrow,
                                                                           _FltVect *
                                                                           _qblk,
                                                                           int _ldq,
                                                                           _FltVect *
                                                                           _tau)
   {

      int i, j, k;
      _FltVect scprod;

      _FltVect *pa, *ph;

// Main cycle over columns

      for (i = 0; i < _ncol; i++) {

// Apply previous columns to the current column

         for (j = 0; j < i; j++) {

            pa = _qblk + i * _ldq + j + 1;
            ph = _qblk + j * _ldq + j + 1;

            scprod = CVector < _FltVect >::ScProd_thr (_nrow - j - 1, pa, ph);

            scprod += _qblk[i * _ldq + j];

            scprod *= _tau[j];

            _qblk[i * _ldq + j] -= scprod;

            pa = _qblk + i * _ldq + j + 1;
            ph = _qblk + j * _ldq + j + 1;

            _FltVect scprod_minus = -scprod;

            CVector < _FltVect >::UpdateVector_thr (_nrow - j - 1, &scprod_minus, ph, pa);

         }

// Compute new transformation

         j = _nrow - 1;
         if (i + 1 < _nrow - 1)
            j = i + 1;

         k = _nrow - i;

         CVector < _FltVect >::Housholder_thr (k, _qblk[i * _ldq + i],
                                               _qblk + i * _ldq + j, _tau[i]);

      }

   }

///
/// @brief Multiply Q factor by the current block
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::MvmQ_Housholder (int _nrhs,
                                                                              int _nrows,
                                                                              int _ncols,
                                                                              _FltVect *
                                                                              _qblk,
                                                                              int _ldq,
                                                                              _FltVect *
                                                                              _tau,
                                                                              _FltVect *
                                                                              _qx,
                                                                              int _ldqx)
   {

// Main cycle over columns

      _FltVect *qblk = (_FltVect *) _qblk;
      _FltVect *tau = (_FltVect *) _tau;
      _FltVect *qxtot = (_FltVect *) _qx;

// Apply Q to the current column

      int j, irhs;

      _FltVect scprod, scprod1;
      _FltVect *pq, *pqx;
      _FltVect *qx;

      for (irhs = 0; irhs < _nrhs; irhs++) {

         qx = qxtot + irhs * _ldqx;

         for (j = _ncols - 1; j >= 0; j--) {

            scprod = qx[j];

            pq = qblk + j * _ldq + j + 1;
            pqx = qx + j + 1;

            scprod1 = CVector < _FltVect >::ScProd (_nrows - j - 1, pq, pqx);
            scprod += scprod1;

            scprod = -scprod * tau[j];

            qx[j] += scprod;

            pq = qblk + j * _ldq + j + 1;
            pqx = qx + j + 1;

            CVector < _FltVect >::UpdateVector (_nrows - j - 1, &scprod, pq, pqx);

         }
      }

   }

///
/// @brief Multiply Q factor by the current block (threads version)
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::MvmQ_Housholder_thr (int _nrhs, int _nrows, int _ncols,
                                       _FltVect * _qblk, int _ldq, _FltVect * _tau,
                                       _FltVect * _qx, int _ldqx)
   {

// Main cycle over columns

      _FltVect *qblk = (_FltVect *) _qblk;
      _FltVect *tau = (_FltVect *) _tau;
      _FltVect *qxtot = (_FltVect *) _qx;

// Apply Q to the current column

      int j, irhs;

      _FltVect scprod, scprod1;
      _FltVect *pq, *pqx;
      _FltVect *qx;

      for (irhs = 0; irhs < _nrhs; irhs++) {

         qx = qxtot + irhs * _ldqx;

         for (j = _ncols - 1; j >= 0; j--) {

            scprod = qx[j];

            pq = qblk + j * _ldq + j + 1;
            pqx = qx + j + 1;

            scprod1 = CVector < _FltVect >::ScProd_thr (_nrows - j - 1, pq, pqx);
            scprod += scprod1;

            scprod = -scprod * tau[j];

            qx[j] += scprod;

            pq = qblk + j * _ldq + j + 1;
            pqx = qx + j + 1;

            CVector < _FltVect >::UpdateVector_thr (_nrows - j - 1, &scprod, pq, pqx);

         }
      }

   }

///
/// @brief Multiply QH factor by the current block
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::MvmQH_Housholder (int _nrhs,
                                                                               int _nrows,
                                                                               int _ncols,
                                                                               _FltVect *
                                                                               _qblk,
                                                                               int _ldq,
                                                                               _FltVect *
                                                                               _tau,
                                                                               _FltVect *
                                                                               _qx,
                                                                               int _ldqx)
   {

// Main cycle over columns

      _FltVect *qblk = (_FltVect *) _qblk;
      _FltVect *tau = (_FltVect *) _tau;
      _FltVect *qxtot = (_FltVect *) _qx;

// Apply Q to the current column

      int j, irhs;

      _FltVect scprod, scprod1;
      _FltVect *pq, *pqx;
      _FltVect *qx;

      for (irhs = 0; irhs < _nrhs; irhs++) {

         qx = qxtot + irhs * _ldqx;

         for (j = 0; j < _ncols; j++) {

            scprod = qx[j];

            pq = qblk + j * _ldq + j + 1;
            pqx = qx + j + 1;

            scprod1 = CVector < _FltVect >::ScProd (_nrows - j - 1, pq, pqx);
            scprod += scprod1;

            scprod = -scprod * tau[j];

            qx[j] += scprod;

            pq = qblk + j * _ldq + j + 1;
            pqx = qx + j + 1;

            CVector < _FltVect >::UpdateVector (_nrows - j - 1, &scprod, pq, pqx);

         }
      }

   }

///
/// @brief Multiply QH factor by the current block (threads version)
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::MvmQH_Housholder_thr (int _nrhs, int _nrows, int _ncols,
                                        _FltVect * _qblk, int _ldq, _FltVect * _tau,
                                        _FltVect * _qx, int _ldqx)
   {

// Main cycle over columns

      _FltVect *qblk = (_FltVect *) _qblk;
      _FltVect *tau = (_FltVect *) _tau;
      _FltVect *qxtot = (_FltVect *) _qx;

// Apply Q to the current column

      int j, irhs;

      _FltVect scprod, scprod1;
      _FltVect *pq, *pqx;
      _FltVect *qx;

      for (irhs = 0; irhs < _nrhs; irhs++) {

         qx = qxtot + irhs * _ldqx;

         for (j = 0; j < _ncols; j++) {

            scprod = qx[j];

            pq = qblk + j * _ldq + j + 1;
            pqx = qx + j + 1;

            scprod1 = CVector < _FltVect >::ScProd_thr (_nrows - j - 1, pq, pqx);
            scprod += scprod1;

            scprod = -scprod * tau[j];

            qx[j] += scprod;

            pq = qblk + j * _ldq + j + 1;
            pqx = qx + j + 1;

            CVector < _FltVect >::UpdateVector_thr (_nrows - j - 1, &scprod, pq, pqx);

         }
      }

   }

///
/// @brief Multiply by block Housholder transformation
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::MvmByBlockHousholder (bool _b_use_thr, char _transp, int _nrow,
                                        int _ncol, int _nrhs, _FltVect * _qx,
                                        const _FltVect * _q, const _FltVect * _tau,
                                        _FltVect * _work)
   {

// Split work memory

      _FltVect *pwork = _work;

      _FltVect *pmatr = pwork;
      pwork += _ncol * _nrhs;

      _FltVect *pmatr1 = pwork;
      pwork += _ncol * _nrhs;

// Multiply

      if (_b_use_thr) {
         CVector < _FltVect >::BlockDot_thr (_nrow, _ncol, _nrhs, _q, _qx, pmatr1, pwork);
      } else {
         CVector < _FltVect >::BlockDot (_nrow, _ncol, _nrhs, _q, _qx, pmatr1);
      }

      if (_transp == 'N' || _transp == 'n') {

         CVector < _FltVect >::MMM ('N', 'N', _nrhs, _ncol, _ncol, pmatr1, _nrhs, _tau,
                                    _ncol, pmatr, _nrhs);

      } else if (_transp == 'T' || _transp == 't') {

         CVector < _FltVect >::MMM ('N', 'T', _nrhs, _ncol, _ncol, pmatr1, _nrhs, _tau,
                                    _ncol, pmatr, _nrhs);

      }

      if (_b_use_thr) {
         CVector < _FltVect >::BlockDaxpy_thr (_nrow, _ncol, _nrhs, _q, pmatr, _qx);
      } else {
         CVector < _FltVect >::BlockDaxpy (_nrow, _ncol, _nrhs, _q, pmatr, _qx);
      }

   }

///
/// @brief Multiply by block Housholder transformation
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::MvmByBlockHousholder (bool _b_use_thr, char _transp, int _nrow,
                                        int _ncol, int _nrhs, _FltVect * _qxdiag,
                                        _FltVect * _qx, const _FltVect * _qdiag,
                                        const _FltVect * _q, _FltVect * _tau,
                                        _FltVect * _work)
   {

// Split work memory

      _FltVect *pwork = _work;

      _FltVect *pmatr = pwork;
      pwork += _ncol * _nrhs;

      _FltVect *pmatr1 = pwork;
      pwork += _ncol * _nrhs;

// Multiply

      CVector < _FltVect >::BlockDot (_ncol, _ncol, _nrhs, _qdiag, _qxdiag, pmatr);

      if (_b_use_thr) {
         CVector < _FltVect >::BlockDot_thr (_nrow, _ncol, _nrhs, _q, _qx, pmatr1, pwork);
      } else {
         CVector < _FltVect >::BlockDot (_nrow, _ncol, _nrhs, _q, _qx, pmatr1);
      }

      CVector < _FltVect >::AddReplaceVector (_ncol * _nrhs, pmatr, pmatr1);

      if (_transp == 'N' || _transp == 'n') {

         CVector < _FltVect >::MMM ('N', 'N', _nrhs, _ncol, _ncol, pmatr1, _nrhs, _tau,
                                    _ncol, pmatr, _nrhs);

      } else if (_transp == 'T' || _transp == 't') {

         CVector < _FltVect >::MMM ('N', 'T', _nrhs, _ncol, _ncol, pmatr1, _nrhs, _tau,
                                    _ncol, pmatr, _nrhs);

      }

      CVector < _FltVect >::BlockDaxpy (_ncol, _ncol, _nrhs, _qdiag, pmatr, _qxdiag);

      if (_b_use_thr) {
         CVector < _FltVect >::BlockDaxpy_thr (_nrow, _ncol, _nrhs, _q, pmatr, _qx);
      } else {
         CVector < _FltVect >::BlockDaxpy (_nrow, _ncol, _nrhs, _q, pmatr, _qx);
      }

   }

/// @brief Solve in-place triangular system
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::SolveR (char _slvtype,
                                                                     int _n,
                                                                     _FltVect * _rmatr,
                                                                     int _ldr,
                                                                     _FltVect * _rhs_sol)
   {

      int i, j;
      if (_slvtype == 'N' || _slvtype == 'n') {
         for (i = _n - 1; i >= 0; i--) {
            _rhs_sol[i] = _rhs_sol[i] / _rmatr[i * _ldr + i];
            for (j = 0; j < i; j++) {
               _rhs_sol[j] -= _rhs_sol[i] * _rmatr[i * _ldr + j];
            }
         }
      } else {
         for (i = 0; i < _n; i++) {
            for (j = 0; j < i; j++) {
               _rhs_sol[i] -= _rhs_sol[j] * _rmatr[i * _ldr + j];
            }
            _rhs_sol[i] = _rhs_sol[i] / _rmatr[i * _ldr + i];
         }
      }

   }

/// @brief Compute Cholessky for matrices stored by columns
//========================================================================================
   template < typename _FltVect > void CVector <
      _FltVect >::CholesskyColumns (double &_diamod, int _n, _FltVect * _amatr, int _lda,
                                    _FltVect * _uarr, int _ldu, double *_dia_arr,
                                    double &_eigmin_att, double &_eigmax_att)
   {

      _FltVect aux;
      double daux;
      int i, j, k;

      for (j = 0; j < _n; j++) {
         for (i = 0; i <= j; i++) {
            _uarr[i + j * _ldu] = _amatr[i + j * _lda];
         }
      }

      _FltVect fzero = (_FltVect) 0.0e0;
      _FltVect fone = (_FltVect) 1.0e0;

      for (j = 0; j < _n; j++) {
         for (i = j + 1; i < _n; i++) {
            _uarr[i + j * _ldu] = fzero;
         }
      }

      for (i = 0; i < _n; i++) {

         for (k = 0; k < i; k++) {
            for (j = i; j < _n; j++) {
               _uarr[i + j * _ldu] -= _uarr[k + i * _ldu] * _uarr[k + j * _ldu];
            }
         }

         aux = _uarr[i + i * _ldu];

         daux = (double) aux;

         _dia_arr[i] = daux;

         if (i == 0 || daux < _eigmin_att)
            _eigmin_att = daux;
         if (i == 0 || daux > _eigmax_att)
            _eigmax_att = daux;

         if (aux < _diamod)
            aux = (_FltVect) _diamod;

         aux = sqrt (aux);

         _uarr[i + i * _ldu] = aux;

         aux = fone / aux;

         for (j = i + 1; j < _n; j++) {
            _uarr[i + j * _ldu] *= aux;
         }

      }

   }

/// @brief Compute SVD
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::ComputeSvd (int _n,
                                                                         _FltVect *
                                                                         _amatr,
                                                                         _FltVect * _sv,
                                                                         _FltVect * _u,
                                                                         _FltVect * _v,
                                                                         double *_work)
   {

      int n_2 = _n * _n;

      double *pamatr = _work;
      double *psv = pamatr + n_2;
      double *pu = psv + _n;
      double *pvh = pu + n_2;
      double *pwork = pvh + n_2;

      int lwork = 10 * (int) _n;

      int i;

      for (i = 0; i < n_2; i++)
         pamatr[i] = (double) _amatr[i];

      if (_n > 0) {

         bool b_found = false;

#ifdef USE_LAPACK
         int n_int = (int) _n;
         int info;
         b_found = true;
         dgesvd_ ("A", "A", &n_int, &n_int, pamatr, &n_int, psv, pu, &n_int, pvh, &n_int,
                  pwork, &lwork, &info);
#else
         (void) lwork;
         (void) pwork;
#endif
         if (!b_found) {
            cout << " Error: Lapack routine DGESVD is not found ! " << endl;
            throw " Error: Lapack routine DGESVD is not found ! ";
         }

      }

      for (i = 0; i < _n; i++)
         _sv[i] = (_FltVect) psv[i];
      for (i = 0; i < n_2; i++)
         _u[i] = (_FltVect) pu[i];

      int j;

      for (i = 0; i < _n; i++) {
         for (j = 0; j < _n; j++) {
            _v[i * _n + j] = (_FltVect) pvh[j * _n + i];
         }
      }

   }

/// @brief Multiply matrix by matrix
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::MMM (int _n,
                                                                  const _FltVect *
                                                                  _x_matr,
                                                                  const _FltVect *
                                                                  _y_matr,
                                                                  _FltVect *
                                                                  _x_times_y_matr)
   {

      _FltVect fzero = (_FltVect) 0.0e0;

      int i, j, k;
      _FltVect aux;

      for (i = 0; i < _n; i++) {
         for (j = 0; j < _n; j++) {
            aux = fzero;
            for (k = 0; k < _n; k++) {
               aux += _x_matr[k * _n + i] * _y_matr[j * _n + k];
            }
            _x_times_y_matr[j * _n + i] = aux;
         }
      }

   }

///
/// @brief  Multiply matrix by matrix (C = op(A)*op(B))
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::MMM (char _atype,
                                                                  char _btype, int _m,
                                                                  int _n, int _k,
                                                                  const _FltVect * _amatr,
                                                                  int _lda,
                                                                  const _FltVect * _bmatr,
                                                                  int _ldb,
                                                                  _FltVect * _cmatr,
                                                                  int _ldc)
   {

      _FltVect fzero;

      CVector < _FltVect >::SetByZeroes (1, &fzero);

      _FltVect *pamatr = (_FltVect *) _amatr;
      _FltVect *pbmatr = (_FltVect *) _bmatr;
      _FltVect *pcmatr = (_FltVect *) _cmatr;

      int i, j, k;
      _FltVect aux;

      if (_atype == 'N' || _atype == 'n') {
         if (_btype == 'N' || _btype == 'n') {
            for (i = 0; i < _m; i++) {
               for (j = 0; j < _n; j++) {
                  aux = fzero;
                  for (k = 0; k < _k; k++) {
                     aux += pamatr[k * _lda + i] * pbmatr[j * _ldb + k];
                  }
                  pcmatr[j * _ldc + i] = aux;
               }
            }
         } else {
            for (i = 0; i < _m; i++) {
               for (j = 0; j < _n; j++) {
                  aux = fzero;
                  for (k = 0; k < _k; k++) {
                     aux += pamatr[k * _lda + i] * pbmatr[k * _ldb + j];
                  }
                  pcmatr[j * _ldc + i] = aux;
               }
            }
         }
      } else {
         if (_btype == 'N' || _btype == 'n') {
            for (i = 0; i < _m; i++) {
               for (j = 0; j < _n; j++) {
                  aux = fzero;
                  for (k = 0; k < _k; k++) {
                     aux += pamatr[i * _lda + k] * pbmatr[j * _ldb + k];
                  }
                  pcmatr[j * _ldc + i] = aux;
               }
            }
         } else {
            for (i = 0; i < _m; i++) {
               for (j = 0; j < _n; j++) {
                  aux = fzero;
                  for (k = 0; k < _k; k++) {
                     aux += pamatr[i * _lda + k] * pbmatr[k * _ldb + j];
                  }
                  pcmatr[j * _ldc + i] = aux;
               }
            }
         }
      }

   }

/// @brief Multiply matrix by vector
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::Mvm (char _transp, int _n,
                                                                  _FltVect * _x_matr,
                                                                  int _ldx,
                                                                  _FltVect * _y_vect,
                                                                  _FltVect *
                                                                  _x_times_y_vect)
   {

      _FltVect fzero = (_FltVect) 0.0e0;

      int i, j;
      _FltVect aux;

      if (_transp == 'N') {

         for (i = 0; i < _n; i++) {
            aux = fzero;
            for (j = 0; j < _n; j++) {
               aux += _x_matr[j * _ldx + i] * _y_vect[j];
            }
            _x_times_y_vect[i] = aux;
         }

      } else {

         for (i = 0; i < _n; i++) {
            aux = fzero;
            for (j = 0; j < _n; j++) {
               aux += _x_matr[i * _ldx + j] * _y_vect[j];
            }
            _x_times_y_vect[i] = aux;
         }

      }

   }

/// @brief Transpose rectangular block
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::TransposeBlock (int _m,
                                                                             int _n,
                                                                             const
                                                                             _FltVect *
                                                                             _a_matr,
                                                                             int _lda,
                                                                             _FltVect *
                                                                             _at_matr,
                                                                             int _ldat)
   {

      int i, j;

      for (i = 0; i < _m; i++) {
         for (j = 0; j < _n; j++) {
            _at_matr[i * _ldat + j] = _a_matr[j * _lda + i];
         }
      }

   }

/// @brief Compute polinomial via square upper Hessenberg matrix
//========================================================================================
   template < typename _FltVect > void CVector < _FltVect >::Polynomial (int _n,
                                                                         int _ncoef,
                                                                         _FltVect *
                                                                         _hmatrix,
                                                                         int _ldh,
                                                                         vector <
                                                                         double >&_coef)
   {

// Allocate work data

      int k = _n;
      int d = _ncoef;

      int k_d = k - d;

      int k_times_kd = k * k_d;

      vector < _FltVect > bmatr_arr ((d + 1) * k_times_kd + 1);
      _FltVect *pbmatr_arr = &bmatr_arr[0];

// Compute b matrices

      _FltVect *pb_0 = pbmatr_arr;
      _FltVect *pb_arr = pbmatr_arr + k_times_kd;

// B_0:

      CVector < _FltVect >::SetByZeroes (k_times_kd, pb_0);

      int i;

      for (i = 0; i < k_d; i++)
         CVector < _FltVect >::SetByOnes (1, pb_0 + i * k + i);

// B_1:

      for (i = 0; i < k_d; i++)
         CVector < _FltVect >::CopyVector (k, _hmatrix + i * _ldh, pb_arr + i * k);

// B_j, j>1:

      int ii, jj, kk;
      _FltVect aux;

      _FltVect *pb_i;
      _FltVect *pb_i_prev;

      for (i = 1; i < d; i++) {

         pb_i = pb_arr + i * k_times_kd;
         pb_i_prev = pb_arr + (i - 1) * k_times_kd;

         CVector < _FltVect >::SetByZeroes (k_times_kd, pb_i);

         for (ii = 0; ii < k; ii++) {
            for (jj = 0; jj < k_d; jj++) {
               aux = (_FltVect) 0.0e0;
               for (kk = 0; kk < k; kk++) {
                  aux += _hmatrix[kk * _ldh + ii] * pb_i_prev[jj * k + kk];
               }
               pb_i[jj * k + ii] = pb_i_prev[jj * k + ii] - aux;
            }
         }

      }

// Compute matrix and rhs

      vector < _FltVect > Ymatr (d * d + 1);
      vector < _FltVect > y (d + 1);

      _FltVect *pYmatr = &Ymatr[0];
      _FltVect *py = &y[0];

      CVector < _FltVect >::SetByZeroes (d * d, pYmatr);
      CVector < _FltVect >::SetByZeroes (d, py);

// Rhs:

      for (i = 0; i < d; i++) {
         pb_i = pb_arr + i * k_times_kd;
         for (ii = 0; ii < k_d; ii++) {
            aux = (_FltVect) 0.0e0;
            for (kk = 0; kk < k; kk++) {
               aux += pb_i[ii * k + kk] * pb_0[ii * k + kk];
            }
            py[i] += aux;
         }

      }

// Matrix:

      _FltVect *pb_j;

      int j;

      for (i = 0; i < d; i++) {
         pb_i = pb_arr + i * k_times_kd;
         for (j = 0; j < d; j++) {
            pb_j = pb_arr + j * k_times_kd;
            for (ii = 0; ii < k_d; ii++) {
               aux = (_FltVect) 0.0e0;
               for (kk = 0; kk < k; kk++) {
                  aux += pb_i[ii * k + kk] * pb_j[ii * k + kk];
               }
               pYmatr[i * d + j] += aux;
            }
         }
      }

// Solve SLE

// Compute fct

      vector < _FltVect > Umatr (d * d + 1);
      vector < _FltVect > x (d + 1);

      vector < double >dia_arr (d + 1);

      _FltVect *pUmatr = &Umatr[0];
      _FltVect *px = &x[0];
      double *pdia_arr = &dia_arr[0];

      double diamod = 1.0e-16;
      double eigmin_att, eigmax_att;

      CVector < _FltVect >::CholesskyColumns (diamod, d, pYmatr, d, pUmatr, d, pdia_arr,
                                              eigmin_att, eigmax_att);

// Solve

      CVector < _FltVect >::CopyVector (d, py, px);

      CVector < _FltVect >::SolveR ('T', d, pUmatr, d, px);
      CVector < _FltVect >::SolveR ('N', d, pUmatr, d, px);

// Store coefs

      _coef.resize (d + 1);
      double *pcoef = &_coef[0];

      for (i = 0; i < d; i++)
         pcoef[i] = (double) px[i];

   }

//
// Set vector data by zeroes
//========================================================================================
   template < typename _IntVect > void CVectorInt < _IntVect >::SetByZeroes (int _n,
                                                                             _IntVect *
                                                                             _ix)
   {

      int i;

      _IntVect izero = (_IntVect) 0.0e0;

      for (i = 0; i < _n; i++)
         _ix[i] = izero;

   }

//
// Set vector data by zeroes (threads version)
//========================================================================================
   template < typename _IntVect > void CVectorInt < _IntVect >::SetByZeroes_thr (int _n,
                                                                                 _IntVect
                                                                                 * _ix)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVectorInt < _IntVect >::SetByZeroes (ni_loc, _ix + ibeg);
         }

      }

   }

//
// Set vector data by ones
//========================================================================================
   template < typename _IntVect > void CVectorInt < _IntVect >::SetByOnes (int _n,
                                                                           _IntVect * _ix)
   {

      int i;

      _IntVect ione = (_IntVect) 1.0e0;

      for (i = 0; i < _n; i++)
         _ix[i] = ione;

   }

//
// Set vector data by zeroes (threads version)
//========================================================================================
   template < typename _IntVect > void CVectorInt < _IntVect >::SetByOnes_thr (int _n,
                                                                               _IntVect *
                                                                               _ix)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVectorInt < _IntVect >::SetByOnes (ni_loc, _ix + ibeg);

         }

      }

   }

//
// Copy vector
//========================================================================================
   template < typename _IntVect > void CVectorInt < _IntVect >::CopyVectorInt (int _n,
                                                                               const
                                                                               _IntVect *
                                                                               _ix,
                                                                               _IntVect *
                                                                               _iy)
   {

      int i;

      for (i = 0; i < _n; i++)
         _iy[i] = _ix[i];

   }

//
// Copy vector (threads version)
//========================================================================================
   template < typename _IntVect > void CVectorInt < _IntVect >::CopyVectorInt_thr (int _n,
                                                                                   const
                                                                                   _IntVect
                                                                                   * _ix,
                                                                                   _IntVect
                                                                                   * _iy)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVectorInt < _IntVect >::CopyVectorInt (ni_loc, _ix + ibeg, _iy + ibeg);

         }

      }

   }

//
// Shift array
//========================================================================================
   template < typename _IntVect > void CVectorInt < _IntVect >::ShiftVectorInt (int _n,
                                                                                const
                                                                                _IntVect
                                                                                _ishift,
                                                                                _IntVect *
                                                                                _iy)
   {

      int i;

      for (i = 0; i < _n; i++)
         _iy[i] += _ishift;

   }

//
// Shift array (threads version)
//========================================================================================
   template < typename _IntVect > void CVectorInt <
      _IntVect >::ShiftVectorInt_thr (int _n, const _IntVect _ishift, _IntVect * _iy)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            if (ni_loc > 0)
               CVectorInt < _IntVect >::ShiftVectorInt (ni_loc, _ishift, _iy + ibeg);

         }

      }

   }

//
// Set shifted identity
//========================================================================================
   template < typename _IntVect > void CVectorInt < _IntVect >::SetIdentity (int _n,
                                                                             const
                                                                             _IntVect
                                                                             _ishift,
                                                                             _IntVect *
                                                                             _iy)
   {

      int i;

      for (i = 0; i < _n; i++)
         _iy[i] = i + _ishift;

   }

//
// Set shifted identity (threads version)
//========================================================================================
   template < typename _IntVect > void CVectorInt < _IntVect >::SetIdentity_thr (int _n,
                                                                                 const
                                                                                 _IntVect
                                                                                 _ishift,
                                                                                 _IntVect
                                                                                 * _iy)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend, ni_loc;
            _IntVect ishift_curr;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            ni_loc = iend + 1 - ibeg;

            ishift_curr = _ishift + ibeg;

            if (ni_loc > 0)
               CVectorInt < _IntVect >::SetIdentity (ni_loc, ishift_curr, _iy + ibeg);

         }

      }

   }

//
// Inverse order ata
//========================================================================================
   template < typename _IntVect > void CVectorInt < _IntVect >::InvOrder (int _n,
                                                                          const _IntVect *
                                                                          _order,
                                                                          _IntVect *
                                                                          _iorder)
   {

      int i;

      for (i = 0; i < _n; i++)
         _iorder[_order[i]] = (_IntVect) i;

   }

//
// Set shifted identity (threads version)
//========================================================================================
   template < typename _IntVect > void CVectorInt < _IntVect >::InvOrder_thr (int _n,
                                                                              const
                                                                              _IntVect *
                                                                              _order,
                                                                              _IntVect *
                                                                              _iorder)
   {

      int n_thr = 1;

#ifdef USE_THREADS
      n_thr = omp_get_max_threads ();
#endif

      int ni_part = _n / n_thr;

      {
#ifdef USE_THREADS
#pragma omp parallel for
#endif
         for (int ipar = 0; ipar < n_thr; ipar++) {

            int ibeg, iend;

            ibeg = ipar * ni_part;
            iend = (ipar + 1) * ni_part - 1;
            if (ipar == n_thr - 1)
               iend = _n - 1;

            int i;

            for (i = ibeg; i <= iend; i++) {
               _iorder[_order[i]] = (_IntVect) i;
            }

         }

      }

   }

/// @brief Create tree
//========================================================================================
   CTree::CTree (int _nnodes_ini, int _nchilds)
   {

// Perform initial tree computation

// Allocate the memory

      int nnodesmax = 2 * _nnodes_ini + 5;

      father.resize (nnodesmax);
      nchilds.resize (nnodesmax);
      childs_list.resize (nnodesmax);
      nodes_lev_id.resize (nnodesmax);

      int *pfather = &father[0];
      int *pnchilds = &nchilds[0];
      vector < int >*pchilds_list = &childs_list[0];
      int *pnodes_lev_id = &nodes_lev_id[0];

// Register first childs for up tree level

      int nnodes_curr = _nnodes_ini;
      int nnodes_prev = 0;
      int ilev_curr = 0;

      int i;

      for (i = 0; i < nnodes_curr; i++) {
         pnchilds[i] = 1;
         pchilds_list[i].resize (1);
         int *ppchilds_list = &pchilds_list[i][0];
         ppchilds_list[0] = i;
         pnodes_lev_id[i] = ilev_curr;
      }

// Init all tree levels

      int nnodes_new;
      int nnodes_prev_curr;

      int j, nchildsloc;

      while (nnodes_curr - nnodes_prev != 1) {

         nnodes_new = nnodes_curr;
         nnodes_prev_curr = nnodes_prev;

         while (nnodes_prev_curr != nnodes_curr) {
            nchildsloc = _nchilds;
            if (nnodes_prev_curr + nchildsloc > nnodes_curr)
               nchildsloc = nnodes_curr - nnodes_prev_curr;
            if (nnodes_curr == nnodes_prev_curr + nchildsloc + 1)
               nchildsloc++;
            for (j = nnodes_prev_curr; j < nnodes_prev_curr + nchildsloc; j++) {
               pfather[j] = nnodes_new;
            }
            pnchilds[nnodes_new] = nchildsloc;
            pchilds_list[nnodes_new].resize (nchildsloc + 1);
            int *ppchilds_list = &pchilds_list[nnodes_new][0];
            for (j = 0; j < nchildsloc; j++)
               ppchilds_list[j] = nnodes_prev_curr + j;
            nnodes_prev_curr += nchildsloc;
            nnodes_new++;
         }

         ilev_curr++;

         for (i = nnodes_curr; i < nnodes_new; i++) {
            pnodes_lev_id[i] = ilev_curr;
         }

         nnodes_prev = nnodes_curr;
         nnodes_curr = nnodes_new;

      }

// Finalize computations with the root node

      root_id = nnodes_curr - 1;

      pfather[root_id] = root_id;

      nlev = ilev_curr + 1;

      nnodes = nnodes_curr;

// Invert level numbers

      int ilev;

      for (i = 0; i < nnodes; i++) {
         ilev = pnodes_lev_id[i];
         pnodes_lev_id[i] = nlev - ilev - 1;
      }

// Fill level arrays

      nnodes_lev.resize (nlev + 1);

      int *pnnodes_lev = &nnodes_lev[0];

      for (i = 0; i < nlev; i++)
         pnnodes_lev[i] = 0;

      for (i = 0; i < nnodes; i++) {
         ilev = pnodes_lev_id[i];
         pnnodes_lev[ilev]++;
      }

      nodes_lev_list.resize (nlev + 1);

      vector < int >*pnodes_lev_list = &nodes_lev_list[0];

      for (i = 0; i < nlev; i++) {
         pnodes_lev_list[i].resize (pnnodes_lev[i] + 1);
      }

      vector < int >iptr (nlev + 1);

      int *piptr = &iptr[0];

      for (i = 0; i < nlev; i++)
         piptr[i] = 0;

      int k;

      for (i = 0; i < nnodes; i++) {
         ilev = pnodes_lev_id[i];
         k = piptr[ilev];
         int *ppnodes_lev_list = &pnodes_lev_list[ilev][0];
         ppnodes_lev_list[k] = i;
         piptr[ilev]++;
      }

      subtree_beg.resize (nnodes + 1);

      int *psubtree_beg = &subtree_beg[0];

      for (i = 0; i < nnodes; i++)
         psubtree_beg[i] = -1;

      this->node2cpu.resize (nnodes + 1);
      this->node2ind.resize (nnodes + 1);
      this->ind2node.resize (nnodes + 1);

      int *pnode2cpu = &this->node2cpu[0];
      int *pnode2ind = &this->node2ind[0];
      int *pind2node = &this->ind2node[0];

      for (i = 0; i < nnodes; i++)
         pnode2cpu[i] = 0;
      for (i = 0; i < nnodes; i++)
         pnode2ind[i] = -1;
      for (i = 0; i < nnodes; i++)
         pind2node[i] = -1;

// Finalize tree computations

      vector < int >ordernd (nnodes + 1);
      int *pordernd = &ordernd[0];

      this->FindOrderingOfNodesAccordingSubtrees (pordernd);

      this->ApplyOrderingOfTreeNodes (pordernd);

      this->InitSubtreeBeg ();

   }

/// @brief Copy operator
//========================================================================================
   CTree & CTree::operator= (const CTree & _tree)
   {

      root_id = _tree.root_id;
      nlev = _tree.nlev;
      nnodes = _tree.nnodes;

      father.resize (nnodes + 1);
      nchilds.resize (nnodes + 1);
      childs_list.resize (nnodes + 1);
      nodes_lev_id.resize (nnodes + 1);
      nnodes_lev.resize (nlev + 1);
      nodes_lev_list.resize (nlev + 1);
      node2cpu.resize (nnodes + 1);
      node2ind.resize (nnodes + 1);
      ind2node.resize (nnodes + 1);
      subtree_beg.resize (nnodes + 1);

      int *pfather = &father[0];
      int *pnchilds = &nchilds[0];
      vector < int >*pchilds_list = &childs_list[0];
      int *pnodes_lev_id = &nodes_lev_id[0];
      int *pnnodes_lev = &nnodes_lev[0];
      vector < int >*pnodes_lev_list = &nodes_lev_list[0];
      int *pnode2cpu = &node2cpu[0];
      int *pnode2ind = &node2ind[0];
      int *pind2node = &ind2node[0];
      int *psubtree_beg = &subtree_beg[0];

      const int *pfather_tree = &_tree.father[0];
      const int *pnchilds_tree = &_tree.nchilds[0];
      const vector < int >*pchilds_list_tree = &_tree.childs_list[0];
      const int *pnodes_lev_id_tree = &_tree.nodes_lev_id[0];
      const int *pnnodes_lev_tree = &_tree.nnodes_lev[0];
      const vector < int >*pnodes_lev_list_tree = &_tree.nodes_lev_list[0];
      const int *pnode2cpu_tree = &_tree.node2cpu[0];
      const int *pnode2ind_tree = &_tree.node2ind[0];
      const int *pind2node_tree = &_tree.ind2node[0];
      const int *psubtree_beg_tree = &_tree.subtree_beg[0];

      int i, j;

      for (i = 0; i < nnodes; i++) {
         pfather[i] = pfather_tree[i];
         pnchilds[i] = pnchilds_tree[i];
         pnodes_lev_id[i] = pnodes_lev_id_tree[i];
         pnode2cpu[i] = pnode2cpu_tree[i];
         pnode2ind[i] = pnode2ind_tree[i];
         pind2node[i] = pind2node_tree[i];
         psubtree_beg[i] = psubtree_beg_tree[i];
      }
      for (i = 0; i < nnodes; i++) {
         pchilds_list[i].resize (pnchilds_tree[i] + 1);
         int *ppchilds_list = &pchilds_list[i][0];
         const int *ppchilds_list_tree = &pchilds_list_tree[i][0];
         for (j = 0; j < pnchilds_tree[i]; j++)
            ppchilds_list[j] = ppchilds_list_tree[j];
      }
      for (i = 0; i < nlev; i++) {
         pnnodes_lev[i] = pnnodes_lev_tree[i];
      }
      for (i = 0; i < nlev; i++) {
         pnodes_lev_list[i].resize (pnnodes_lev_tree[i] + 1);
         int *ppnodes_lev_list = &pnodes_lev_list[i][0];
         const int *ppnodes_lev_list_tree = &pnodes_lev_list_tree[i][0];
         for (j = 0; j < pnnodes_lev_tree[i]; j++)
            ppnodes_lev_list[j] = ppnodes_lev_list_tree[j];
      }

      return *this;

   }

/// @brief Find ordering of nodes as subtrees
//========================================================================================
   void CTree::FindOrderingOfNodesAccordingSubtrees (int *_ordernd)
   {

// Open tree structures

      int *pfather = &father[0];
      int *pnchilds = &nchilds[0];
      vector < int >*pchilds_list = &childs_list[0];
      int *pnodes_lev_id = &nodes_lev_id[0];

// Check monotonicity of node levels

      vector < int >nz_lev (nlev + 2);
      vector < int >ibs_lev (nlev + 2);

      int *pnz_lev = &nz_lev[0];
      int *pibs_lev = &ibs_lev[0];

      int i;

      for (i = 0; i <= nlev; i++)
         pnz_lev[i] = 0;
      for (i = 0; i <= nlev; i++)
         pibs_lev[i] = 0;

      int ilev = nlev - 1;

      for (i = 0; i < nnodes; i++) {
         if (pnodes_lev_id[i] < ilev) {
            ilev = pnodes_lev_id[i];
         } else if (pnodes_lev_id[i] > ilev) {
            throw
               " CTree::FindOrderingOfNodesAccordingSubtrees: node levels are not monotone ";
         }
         pnz_lev[ilev]++;
      }

// Create ordering array

      int inode = 0;
      int inodenew = 0;

      int ifather, inodecurr, nchildsloc;
      int *ppchilds;

      while (inode < pnz_lev[nlev - 1]) {
         _ordernd[inode] = inodenew;
         inodenew++;
         inodecurr = inode;
         while (true) {
            ifather = pfather[inodecurr];
            if (ifather != inodecurr) {
               ppchilds = &pchilds_list[ifather][0];
               nchildsloc = pnchilds[ifather];
               if (ppchilds[nchildsloc - 1] == inodecurr) {
                  _ordernd[ifather] = inodenew;
                  inodenew++;
                  inodecurr = ifather;
               } else {
                  break;
               }
            } else {
               break;
            }
         }
         inode++;
      }

   }

/// @brief Apply ordering of nodes as subtrees
//========================================================================================
   void CTree::ApplyOrderingOfTreeNodes (int *_ordernd)
   {

// Open tree structures

      int root_id_loc = root_id;

      int root_id_new = _ordernd[root_id_loc];

      root_id = root_id_new;

      int *pfather = &father[0];
      int *pnchilds = &nchilds[0];

      int i;

      vector < int >fathernew (nnodes + 1);
      vector < int >nchildsnew (nnodes + 1);

      int *pfathernew = &fathernew[0];
      int *pnchildsnew = &nchildsnew[0];

      int inew, ifather, ifathernew;

      for (i = 0; i < nnodes; i++) {
         inew = _ordernd[i];
         ifather = pfather[i];
         ifathernew = _ordernd[ifather];
         pfathernew[inew] = ifathernew;
      }

      father.swap (fathernew);

      int nchildsloc;

      for (i = 0; i < nnodes; i++) {
         inew = _ordernd[i];
         nchildsloc = pnchilds[i];
         pnchildsnew[inew] = nchildsloc;
      }

      nchilds.swap (nchildsnew);

      pnchilds = &nchilds[0];

      int *pnodes_lev_id = &nodes_lev_id[0];

      vector < int >nodes_lev_id_new (nnodes);

      int *pnodes_lev_id_new = &nodes_lev_id_new[0];

      int ilev;

      for (i = 0; i < nnodes; i++) {
         inew = _ordernd[i];
         ilev = pnodes_lev_id[i];
         pnodes_lev_id_new[inew] = ilev;
      }

      nodes_lev_id.swap (nodes_lev_id_new);

      vector < int >*pchilds_list = &childs_list[0];

      vector < vector < int > >childs_list_new (nnodes);

      vector < int >*pchilds_list_new = &childs_list_new[0];

      int *pchildsloc;

      int j, ichild, ichildnew;

      for (i = 0; i < nnodes; i++) {
         inew = _ordernd[i];
         nchildsloc = pnchilds[inew];
         pchilds_list_new[inew].swap (pchilds_list[i]);
         pchildsloc = &pchilds_list_new[inew][0];
         for (j = 0; j < nchildsloc; j++) {
            ichild = pchildsloc[j];
            ichildnew = _ordernd[ichild];
            pchildsloc[j] = ichildnew;
         }
      }

      childs_list.swap (childs_list_new);

      vector < int >iptrlev (nlev + 1);

      int *piptrlev = &iptrlev[0];

      for (i = 0; i <= nlev; i++)
         piptrlev[i] = 0;

      vector < int >*pnodes_lev_list = &nodes_lev_list[0];

      pnodes_lev_id = &nodes_lev_id[0];

      int *pnodeslevloc;

      int k;

      for (i = 0; i < nnodes; i++) {
         ilev = pnodes_lev_id[i];
         pnodeslevloc = &pnodes_lev_list[ilev][0];
         k = piptrlev[ilev];
         pnodeslevloc[k] = i;
         piptrlev[ilev]++;
      }

   }

/// @brief Compute common node in up part of a tree
//========================================================================================
   int CTree::FindCommonNode (int _nlistnd, int *_listnd, int &_icycle, int *_imasknd)
   {

      int inode_curr = -1;

      if (_nlistnd == 0)
         return inode_curr;
      if (_nlistnd == 1)
         return _listnd[0];

// Find lowest level number

      int nnodes_loc = this->GetNnodes ();
      int *pnodes_lev_id = &nodes_lev_id[0];
      int *pfather = &father[0];

      int inodeloc = _listnd[0];

      int ilev_min = pnodes_lev_id[inodeloc];

      int i, ilevloc;

      for (i = 1; i < _nlistnd; i++) {
         inodeloc = _listnd[i];
         ilevloc = pnodes_lev_id[inodeloc];
         if (ilev_min > ilevloc)
            ilev_min = ilevloc;
      }

// Create initial list of nodes on the same level

      int *pwork1 = _imasknd + nnodes_loc;
      int *pwork2 = pwork1 + nnodes_loc;

      _icycle++;

      int nlist_curr = 0;

      for (i = 0; i < _nlistnd; i++) {
         inodeloc = _listnd[i];
         ilevloc = pnodes_lev_id[inodeloc];
         while (ilevloc != ilev_min) {
            inodeloc = pfather[inodeloc];
            ilevloc = pnodes_lev_id[inodeloc];
         }
         if (_imasknd[inodeloc] != _icycle) {
            pwork1[nlist_curr] = inodeloc;
            nlist_curr++;
            _imasknd[inodeloc] = _icycle;
         }
      }

      if (nlist_curr > nnodes_loc) {
         cout << " CTree::FindCommonNode: error: incorrect number of nodes in a list ! "
            << endl;
      }
// Find up level common root

      int *pwork_curr = pwork1;
      int *pwork_next = pwork2;

      int nlist_next;

      while (nlist_curr != 1) {

         _icycle++;

         nlist_next = 0;

         for (i = 0; i < nlist_curr; i++) {
            inodeloc = pwork_curr[i];
            inodeloc = pfather[inodeloc];
            if (_imasknd[inodeloc] != _icycle) {
               pwork_next[nlist_next] = inodeloc;
               nlist_next++;
               _imasknd[inodeloc] = _icycle;
            }
         }

         if (pwork_curr == pwork1) {
            pwork_curr = pwork2;
            pwork_next = pwork1;
         } else {
            pwork_curr = pwork1;
            pwork_next = pwork2;
         }

         nlist_curr = nlist_next;

      }

      inode_curr = pwork_curr[0];

      return inode_curr;

   }

/// @brief For each tree node find its subtree start
//========================================================================================
   void CTree::ComputeSubtreeStart (int *_subtree_start)
   {

      int *pnchilds = &nchilds[0];
      vector < int >*pchilds_list = &childs_list[0];

      int i, ichild0;
      int *pchildsloc;

      for (i = 0; i < nnodes; i++)
         _subtree_start[i] = -1;

      for (i = 0; i < nnodes; i++) {
         if (pnchilds[i] == 1) {
            _subtree_start[i] = i;
         } else {
            pchildsloc = &pchilds_list[i][0];
            ichild0 = pchildsloc[0];
            _subtree_start[i] = _subtree_start[ichild0];
         }
      }

   }

/// @brief Perform tree filtering
//========================================================================================
   void CTree::FilterTree (int *_imasknd, CTree & _tree_flt)
   {

// Compute old2new ordering

      int *pnodes_lev_id_curr = this->GetNodesLevId ();

      vector < int >ordernd (this->nnodes + 1);
      vector < int >listnd (this->nnodes + 1);

      int *pordernd = &ordernd[0];
      int *plistnd = &listnd[0];

      int i;

      for (i = 0; i < this->nnodes; i++)
         pordernd[i] = -1;

      int nnodes_new = 0;
      int nlev_new = 0;
      int ilev;

      for (i = 0; i < this->nnodes; i++) {
         if (_imasknd[i] > 0) {
            pordernd[i] = nnodes_new;
            plistnd[nnodes_new] = i;
            nnodes_new++;
            ilev = pnodes_lev_id_curr[i];
            if (ilev >= nlev_new)
               nlev_new = ilev + 1;
         }
      }

      _tree_flt.AllocateTree (nnodes_new, nlev_new);

      _tree_flt.root_id = pordernd[this->root_id];
      _tree_flt.nlev = nlev_new;
      _tree_flt.nnodes = nnodes_new;

      int *pfather = this->GetFather ();
      int *pfather_new = _tree_flt.GetFather ();

      int inode_old, father_old;

      for (i = 0; i < nnodes_new; i++) {
         inode_old = plistnd[i];
         father_old = pfather[inode_old];
         pfather_new[i] = pordernd[father_old];
      };

      int *pnchilds = this->GetNchilds ();
      vector < int >*pchilds_list = this->GetChildsList ();

      int *pnchilds_new = _tree_flt.GetNchilds ();
      vector < int >*pchilds_list_new = _tree_flt.GetChildsList ();

      int nchilds_old, nchilds_new, j, ichild, ichild_new;

      for (i = 0; i < nnodes_new; i++) {
         inode_old = plistnd[i];
         nchilds_old = pnchilds[inode_old];
         int *pchilds_old = &pchilds_list[inode_old][0];
         nchilds_new = 0;
         for (j = 0; j < nchilds_old; j++) {
            ichild = pchilds_old[j];
            ichild_new = pordernd[ichild];
            if (ichild_new >= 0)
               nchilds_new++;
         }
         if (nchilds_new == 0)
            nchilds_new = 1;
         pnchilds_new[i] = nchilds_new;
      }

      for (i = 0; i < nnodes_new; i++)
         pchilds_list_new[i].resize (pnchilds_new[i] + 1);

      for (i = 0; i < nnodes_new; i++) {
         inode_old = plistnd[i];
         nchilds_old = pnchilds[inode_old];
         int *pchilds_old = &pchilds_list[inode_old][0];
         int *pchilds_new = &pchilds_list_new[i][0];
         nchilds_new = 0;
         for (j = 0; j < nchilds_old; j++) {
            ichild = pchilds_old[j];
            ichild_new = pordernd[ichild];
            if (ichild_new >= 0) {
               pchilds_new[nchilds_new] = ichild_new;
               nchilds_new++;
            }
         }
         if (nchilds_new == 0) {
            pchilds_new[0] = i;
         }
      }

      int *pnodes_lev_id = this->GetNodesLevId ();
      int *pnodes_lev_id_new = _tree_flt.GetNodesLevId ();

      for (i = 0; i < nnodes_new; i++) {
         inode_old = plistnd[i];
         pnodes_lev_id_new[i] = pnodes_lev_id[inode_old];
      }

      _tree_flt.nnodes_lev.resize (nlev_new + 1);

      int *pnnodes_lev_new = _tree_flt.GetNNodesLev ();

      for (i = 0; i < nlev_new; i++)
         pnnodes_lev_new[i] = 0;

      for (i = 0; i < nnodes_new; i++) {
         ilev = pnodes_lev_id_new[i];
         pnnodes_lev_new[ilev]++;
      }

      vector < int >*pnodes_lev_list = &_tree_flt.nodes_lev_list[0];

      for (i = 0; i < nlev_new; i++)
         pnodes_lev_list[i].resize (pnnodes_lev_new[i] + 1);

      for (i = 0; i < nlev_new; i++)
         pnnodes_lev_new[i] = 0;

      for (i = 0; i < nnodes_new; i++) {
         ilev = pnodes_lev_id_new[i];
         int *ppnodes_lev_list = &pnodes_lev_list[ilev][0];
         ppnodes_lev_list[pnnodes_lev_new[ilev]] = i;
         pnnodes_lev_new[ilev]++;
      }

      int *pnode2cpu = this->GetNode2Cpu ();
      int *pnode2cpu_new = _tree_flt.GetNode2Cpu ();

      for (i = 0; i < nnodes_new; i++) {
         inode_old = plistnd[i];
         pnode2cpu_new[i] = pnode2cpu[inode_old];
      }

      int *pnode2ind = this->GetNode2Ind ();
      int *pnode2ind_new = _tree_flt.GetNode2Ind ();

      for (i = 0; i < nnodes_new; i++) {
         inode_old = plistnd[i];
         pnode2ind_new[i] = pnode2ind[inode_old];
      }

      int *pind2node = this->GetInd2Node ();
      (void)pind2node;
      int *pind2node_new = _tree_flt.GetInd2Node ();

      for (i = 0; i < nnodes_new; i++)
         pind2node_new[i] = -1;

      int ind;

      for (i = 0; i < nnodes_new; i++) {
         ind = pnode2ind_new[i];
         if (ind >= nnodes_new) {
            cout << " CTree::FilterTree: error: incorrect index number !" << endl;
            throw " CTree::FilterTree: error: incorrect index number !";
         }
         if (ind >= 0) {
            pind2node_new[ind] = i;
         }
      }

      _tree_flt.InitSubtreeBeg ();

   }

/// @brief Pack tree
//========================================================================================
   void CTree::PackTree (vector < char >&_obj)
   {

      int length = this->GetPackedTreeSize ();

      _obj.resize (length + 1);

      char *pobj = &_obj[0];

      FillPackedTree (length, pobj);

   }

/// @brief Get packed tree size
//========================================================================================
   int CTree::GetPackedTreeSize ()
   {

      int *pnchilds = &nchilds[0];

      int nchildstot = 0;

      int i;

      for (i = 0; i < nnodes; i++)
         nchildstot += pnchilds[i];

      int size = 4 * sizeof (int) + (nnodes * 7 + nlev + nchildstot) * sizeof (int);

      return size;

   }

/// @brief Fill packed tree
//========================================================================================
   void CTree::FillPackedTree (int _length, char *_obj)
   {

      int *pfather = &father[0];
      int *pnchilds = &nchilds[0];
      vector < int >*pchilds_list = &childs_list[0];
      int *pnodes_lev_id = &nodes_lev_id[0];
      int *pnnodes_lev = &nnodes_lev[0];
      vector < int >*pnodes_lev_list = &nodes_lev_list[0];
      int *pnode2cpu = &node2cpu[0];
      int *pnode2ind = &node2ind[0];
      int *pind2node = &ind2node[0];

      int nchildstot = 0;

      int i;

      for (i = 0; i < nnodes; i++)
         nchildstot += pnchilds[i];

      char *pLoc;

      pLoc = _obj;

      int *pHead;
      int *pfather_obj;
      int *pnchilds_obj;
      int *pchilds_obj;
      int *pnodeslev_obj;
      int *pnnodeslev_obj;
      int *pnodeslevlist_obj;
      int *pnode2cpu_obj;
      int *pnode2ind_obj;
      int *pind2node_obj;

      pHead = (int *) pLoc;
      pLoc += 4 * sizeof (int);

      pHead[0] = root_id;
      pHead[1] = nlev;
      pHead[2] = nnodes;
      pHead[3] = nchildstot;

      pfather_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);
      pnchilds_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);
      pchilds_obj = (int *) pLoc;
      pLoc += nchildstot * sizeof (int);
      pnodeslev_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);
      pnnodeslev_obj = (int *) pLoc;
      pLoc += nlev * sizeof (int);
      pnodeslevlist_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);
      pnode2cpu_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);
      pnode2ind_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);
      pind2node_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);

      if (pLoc - _obj != _length) {
         throw " CTree::FillPackedTree: incorrect length on entry ";
      }
// Fill arrays

      for (i = 0; i < nnodes; i++)
         pfather_obj[i] = pfather[i];
      for (i = 0; i < nnodes; i++)
         pnchilds_obj[i] = pnchilds[i];
      for (i = 0; i < nnodes; i++)
         pnodeslev_obj[i] = pnodes_lev_id[i];
      for (i = 0; i < nlev; i++)
         pnnodeslev_obj[i] = pnnodes_lev[i];

      nchildstot = 0;

      int j;

      for (i = 0; i < nnodes; i++) {
         int *pchildsloc = &pchilds_list[i][0];
         for (j = 0; j < pnchilds[i]; j++) {
            pchilds_obj[nchildstot] = pchildsloc[j];
            nchildstot++;
         }
      }

      int nnodesloc = 0;

      for (i = 0; i < nlev; i++) {
         int *ppnodes_lev_list = &pnodes_lev_list[i][0];
         for (j = 0; j < pnnodes_lev[i]; j++) {
            pnodeslevlist_obj[nnodesloc] = ppnodes_lev_list[j];
            nnodesloc++;
         }
      }

      for (i = 0; i < nnodes; i++)
         pnode2cpu_obj[i] = pnode2cpu[i];

      for (i = 0; i < nnodes; i++)
         pnode2ind_obj[i] = pnode2ind[i];

      for (i = 0; i < nnodes; i++)
         pind2node_obj[i] = pind2node[i];

   }

/// @brief Unpack tree
//========================================================================================
   void CTree::UnPackTree (int _length, char *_obj)
   {

// Get head data

      char *pLoc;

      pLoc = _obj;

      int *pHead;

      pHead = (int *) pLoc;
      pLoc += 4 * sizeof (int);

      root_id = pHead[0];

      int nlev_curr = pHead[1];
      int nnodes_curr = pHead[2];
      int nchildstot = pHead[3];

      this->AllocateTree (nnodes_curr, nlev_curr);

      nnodes = nnodes_curr;
      nlev = nlev_curr;

      int *pfather_obj;
      int *pnchilds_obj;
      int *pchilds_obj;
      int *pnodeslev_obj;
      int *pnnodeslev_obj;
      int *pnodeslevlist_obj;
      int *pnode2cpu_obj;
      int *pnode2ind_obj;
      int *pind2node_obj;

      pfather_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);
      pnchilds_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);
      pchilds_obj = (int *) pLoc;
      pLoc += nchildstot * sizeof (int);
      pnodeslev_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);
      pnnodeslev_obj = (int *) pLoc;
      pLoc += nlev * sizeof (int);
      pnodeslevlist_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);
      pnode2cpu_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);
      pnode2ind_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);
      pind2node_obj = (int *) pLoc;
      pLoc += nnodes * sizeof (int);

      if (pLoc - _obj != _length) {
         throw " CTree::FillPackedTree: incorrect length on entry ";
      }

      int *pfather = &father[0];
      int *pnchilds = &nchilds[0];
      vector < int >*pchilds_list = &childs_list[0];
      int *pnodes_lev_id = &nodes_lev_id[0];
      int *pnnodes_lev = &nnodes_lev[0];
      vector < int >*pnodes_lev_list = &nodes_lev_list[0];
      int *pnode2cpu = &node2cpu[0];
      int *pnode2ind = &node2ind[0];
      int *pind2node = &ind2node[0];

      int i, j;

      for (i = 0; i < nnodes; i++)
         pfather[i] = pfather_obj[i];
      for (i = 0; i < nnodes; i++)
         pnchilds[i] = pnchilds_obj[i];
      for (i = 0; i < nnodes; i++)
         pnodes_lev_id[i] = pnodeslev_obj[i];
      for (i = 0; i < nlev; i++)
         pnnodes_lev[i] = pnnodeslev_obj[i];

      for (i = 0; i < nnodes; i++)
         pchilds_list[i].resize (pnchilds[i]);
      for (i = 0; i < nlev; i++)
         pnodes_lev_list[i].resize (pnnodes_lev[i]);

      nchildstot = 0;

      for (i = 0; i < nnodes; i++) {
         int *pchildsloc = &pchilds_list[i][0];
         for (j = 0; j < pnchilds[i]; j++) {
            pchildsloc[j] = pchilds_obj[nchildstot];
            nchildstot++;
         }
      }

      int nnodesloc = 0;

      for (i = 0; i < nlev; i++) {
         int *ppnodes_lev_list = &pnodes_lev_list[i][0];
         for (j = 0; j < pnnodes_lev[i]; j++) {
            ppnodes_lev_list[j] = pnodeslevlist_obj[nnodesloc];
            nnodesloc++;
         }
      }

      for (i = 0; i < nnodes; i++)
         pnode2cpu[i] = pnode2cpu_obj[i];

      for (i = 0; i < nnodes; i++)
         pnode2ind[i] = pnode2ind_obj[i];

      for (i = 0; i < nnodes; i++)
         pind2node[i] = pind2node_obj[i];

   }

/// @brief Output tree
//========================================================================================
   void CTree::OutputTree (ostream & _stream)
   {

      _stream << " CTTree: " << endl;
      _stream << "    Root_id = " << root_id << endl;
      _stream << "    Nlev = " << nlev << endl;
      _stream << "    Nnodes = " << nnodes << endl;

      int *pfather = &father[0];
      int *pnchilds = &nchilds[0];
      vector < int >*pchilds_list = &childs_list[0];
      int *pnodes_lev_id = &nodes_lev_id[0];
      int *pnnodes_lev = &nnodes_lev[0];
      vector < int >*pnodes_lev_list = &nodes_lev_list[0];
      int *pnode2cpu = &node2cpu[0];
      int *pnode2ind = &node2ind[0];
      int *pind2node = &ind2node[0];

      PrintArray (_stream, " Father ", (int) nnodes, pfather);
      PrintArray (_stream, " Nchilds ", (int) nnodes, pnchilds);
      int i;
      for (i = 0; i < nnodes; i++) {
         _stream << " Inode = " << i << endl;
         int *ppchilds_list = &pchilds_list[i][0];
         PrintArray (_stream, " Childs ", (int) pnchilds[i], ppchilds_list);
      }
      PrintArray (_stream, " Nodes_lev_id ", (int) nnodes, pnodes_lev_id);
      PrintArray (_stream, " Nnodes_lev ", (int) nlev, pnnodes_lev);
      for (i = 0; i < nlev; i++) {
         _stream << " Ilev = " << i << endl;
         int *ppnodes_lev_list = &pnodes_lev_list[i][0];
         PrintArray (_stream, " Nodes ", (int) pnnodes_lev[i], ppnodes_lev_list);
      }
      PrintArray (_stream, " Node2Cpu ", (int) nnodes, pnode2cpu);
      PrintArray (_stream, " Node2Ind ", (int) nnodes, pnode2ind);
      PrintArray (_stream, " Ind2Node ", (int) nnodes, pind2node);

   }

   template class CVector < double >;
   template class CVector < float >;
   template class CVectorInt < int >;
   template class CVectorInt < long long >;
   template class CFct_impl < long long, double >;
   template class CFct_impl < int, double >;
   template class CFct_impl < long long, float >;
   template class CFct_impl < int, float >;
   template class CMatrix < long long, double >;
   template class CMatrix < int, double >;
   template class CMatrix < long long, float >;
   template class CMatrix < int, float >;

////////////////////////////////////////////////////////////////////////////////////////////////
//...reordering algorithms (reverse Cuthill-Mckee ordering);
   void genrcm (int n, int *xadj, int *iadj, int *perm, int *xls, int *mask);
   void subrcm (int *xadj, int *iadj, int *mask, int nsubg, int *subg, int *perm,
                int *xls, int n);
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
   void rootls (int root, int *xadj, int *iadj, int *mask, int &nlvl, int *xls, int *ls,
                int)
   {
      int i, iccsze, j, jstop, jstrt, lbegin, lvlend, lvsize, nbr, node;

      mask[root - 1] = 0;
      ls[0] = root;
      nlvl = 0;
      lvlend = 0;
      iccsze = 1;

      do {
         lbegin = lvlend + 1;
         lvlend = iccsze;
         xls[nlvl++] = lbegin;

/////////////////////////////////////////////////////////////////////////////////////////////
//...Generate the next level by finding all the masked neighbors of nodes in the current level;
         for (i = lbegin; i <= lvlend; i++) {
            jstrt = xadj[(node = ls[i - 1]) - 1];
            jstop = xadj[node] - 1;

            for (j = jstrt; j <= jstop; j++)
               if (mask[(nbr = iadj[j - 1]) - 1] != 0) {
                  ls[iccsze++] = nbr;
                  mask[nbr - 1] = 0;
               }
         }
      }
      while ((lvsize = iccsze - lvlend) > 0);

////////////////////////////////////////////////////////////
//...Reset MASK to one for the nodes in the level structure;
      for (xls[nlvl] = lvlend + 1, i = 1; i <= iccsze; i++)
         mask[(node = ls[i - 1]) - 1] = 1;
   }

///////////////////////////////////
//...finds pseudo-peripheral nodes;
   void fnroot (int &root, int *xadj, int *iadj, int *mask, int &nlvl, int *xls, int *ls,
                int n)
   {
      int iccsze, j, jstrt, k, kstop, kstrt, mindeg, nabor, ndeg, node, nunlvl;

      rootls (root, xadj, iadj, mask, nlvl, xls, ls, n);
      if (nlvl == 1 || nlvl == (iccsze = xls[nlvl] - 1))
         return;

      do {
         mindeg = iccsze;
         root = ls[(jstrt = xls[nlvl - 1]) - 1];

         if (iccsze > jstrt) {
            for (j = jstrt; j <= iccsze; j++) {
               ndeg = 0;
               kstrt = xadj[(node = ls[j - 1]) - 1];
               kstop = xadj[node] - 1;

               for (k = kstrt; k <= kstop; k++)
                  if (mask[(nabor = iadj[k - 1]) - 1] > 0)
                     ndeg++;

               if (ndeg < mindeg) {
                  root = node;
                  mindeg = ndeg;
               }
            }
         }

         rootls (root, xadj, iadj, mask, nunlvl, xls, ls, n);
         if (nunlvl <= nlvl)
            return;
      }
      while ((nlvl = nunlvl) < iccsze);
   }

//////////////////////////////////////////////////////////////////
//...computes the degrees of the nodes in the connected component;
   void degree (int root, int *xadj, int *iadj, int *mask, int *deg, int &iccsze, int *ls,
                int)
   {
      int i, ideg, j, jstop, jstrt, lbegin, lvlend, lvsize, nbr, node;

      ls[0] = root;
      lvlend = 0;
      iccsze = 1;
      xadj[root - 1] = -xadj[root - 1];

      do {
         lbegin = lvlend + 1;
         lvlend = iccsze;

         for (i = lbegin; i <= lvlend; i++) {
            jstrt = -xadj[(node = ls[i - 1]) - 1];
//          jstop =  abs(xadj[node])-1;
            jstop = xadj[node];
            if (jstop < 0)
               jstop = -jstop;
            jstop--;
            ideg = 0;

            for (j = jstrt; j <= jstop; j++)
               if (mask[(nbr = iadj[j - 1]) - 1] != 0) {
                  ideg = ideg + 1;
                  if (xadj[nbr - 1] >= 0) {
                     xadj[nbr - 1] = -xadj[nbr - 1];
                     ls[iccsze++] = nbr;
                  }
               }
            deg[node - 1] = ideg;
         }
      }
      while ((lvsize = iccsze - lvlend) > 0);

///////////////////////////////////////////////
//...Reset XADJ to its correct sign and return;
      for (i = 1; i <= iccsze; i++) {
         node = ls[i - 1];
         xadj[node - 1] = -xadj[node - 1];
      }
   }


////////////////////////////////////////////////
//...reverses the elements of an integer vector;
#define iSWAP(A, B) { int iSWAP_temp = (A); (A) = (B); (B) = iSWAP_temp; }

   void ivec_reverse (int n, int *a)
   {
      int m, i;
      for (m = n / 2, i = 1; i <= m; i++)
         iSWAP (a[i - 1], a[n - i]);
   }

#undef iSWAP

/////////////////////////////////////////////////////////////////////////////
//...numbers a connected component using the reverse Cuthill McKee algorithm;
   void rcm (int root, int *xadj, int *iadj, int *mask, int *perm, int &iccsze, int *deg,
             int n)
   {
      int fnbr, i, j, jstop, jstrt, k, l, lbegin, lnbr, lperm, lvlend, nbr, node;

      degree (root, xadj, iadj, mask, deg, iccsze, perm, n);
      mask[root - 1] = 0;
      if (iccsze <= 1)
         return;

      lvlend = 0;
      lnbr = 1;

      do {
         lbegin = lvlend + 1;
         lvlend = lnbr;

         for (i = lbegin; i <= lvlend; i++) {
            jstrt = xadj[(node = perm[i - 1]) - 1];
            jstop = xadj[node] - 1;
            fnbr = lnbr + 1;

            for (j = jstrt; j <= jstop; j++)
               if (mask[(nbr = iadj[j - 1]) - 1] != 0) {
                  mask[nbr - 1] = 0;
                  perm[lnbr++] = nbr;
               }
/////////////////////////////////////////////////////////////
//...Sort the neighbors of node in increasing order by degree;
            if (fnbr < lnbr) {
               k = fnbr;
               do {
                  l = k;
                  nbr = perm[k++];
                label40:
                  if (l > fnbr && deg[(lperm = perm[l - 1]) - 1] > deg[nbr - 1]) {
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
      ivec_reverse (iccsze, perm);
   }

//////////////////////////////////////////////////////////////////
//...finds the reverse Cuthill-Mckee ordering for a general graph;
   void genrcm (int n, int *xadj, int *iadj, int *perm, int *xls, int *mask)
   {
      int i, iccsze, nlvl, num, root;

      for (i = 1; i <= n; i++)
         mask[i - 1] = 1;
      for (num = 1, i = 1; i <= n; i++)
         if (mask[i - 1] != 0) {
            fnroot (root = i, xadj, iadj, mask, nlvl, xls, perm + num - 1, n);
            rcm (root, xadj, iadj, mask, perm + num - 1, iccsze, xls, n);

            if ((num += iccsze) > n)
               return;
         }
   }

///////////////////////////////////////////////////////////////////
//...finds the reverse Cuthill-Mckee ordering for a given subgraph;
   void subrcm (int *xadj, int *iadj, int *mask, int nsubg, int *subg, int *perm,
                int *xls, int n)
   {
      int i, iccsze, nlvl, node, num;
      for (i = 1; i <= nsubg; i++)
         mask[(node = subg[i - 1]) - 1] = 1;

      for (num = 0, i = 1; i <= nsubg; i++)
         if (mask[(node = subg[i - 1]) - 1] > 0) {
            fnroot (node, xadj, iadj, mask, nlvl, xls, perm + num, n);
            rcm (node, xadj, iadj, mask, perm + num, iccsze, xls, n);

            if ((num += iccsze) >= nsubg)
               return;
         }
   }

// Barrier synchronization of the CPUs
//========================================================================================
   void CMPIDataExchange::Synchronize (void *_comm)
   {
#ifdef USE_MPI
      MPI_Comm *pcomm = (MPI_Comm *) _comm;
      MPI_Barrier (*pcomm);
#endif
   }

// Get the number of CPUs
//========================================================================================
   int CMPIDataExchange::GetNproc (void *_comm)
   {
#ifdef USE_MPI
      MPI_Comm *pcomm = (MPI_Comm *) _comm;

      int nproc;
      MPI_Comm_size (*pcomm, &nproc);

      return nproc;
#else
      return 1;
#endif
   }

// Get the ID number of the CPU
//========================================================================================
   int CMPIDataExchange::GetMyid (void *_comm)
   {
#ifdef USE_MPI
      MPI_Comm *pcomm = (MPI_Comm *) _comm;

      int myid;
      MPI_Comm_rank (*pcomm, &myid);

      return myid;
#else
      return 0;
#endif
   }

// Description: Synchronous send data to the other CPU (implementation)
//========================================================================================
   void CMPIDataExchange::Send (void *_comm, int _rank, int _msgtag, int _length,
                                char *_buffer)
   {
#ifdef USE_MPI
      MPI_Comm *pcomm = (MPI_Comm *) _comm;

      MPI_Send (_buffer, _length, MPI_CHAR, _rank, _msgtag, *(pcomm));
#endif
   }

// Description: Asynchronous send data to the other CPU (implementation)
//========================================================================================
   void CMPIDataExchange::ISend (void *_comm, int _rank, int _msgtag, int _length,
                                 char *_buffer, int _indrecv, void *_recv_arr)
   {
#ifdef USE_MPI
      MPI_Comm *pcomm = (MPI_Comm *) _comm;

      MPI_Request *p_recv_arr = (MPI_Request *) _recv_arr;

      MPI_Isend (_buffer, _length, MPI_CHAR, _rank, _msgtag, *(pcomm),
                 p_recv_arr + _indrecv);
#endif
   }

// Description: Synchronous receive data from the other CPU (implementation)
//========================================================================================
   void CMPIDataExchange::Recv (void *_comm, int _rank, int _msgtag, int _length,
                                char *_buffer)
   {
#ifdef USE_MPI
      MPI_Comm *pcomm = (MPI_Comm *) _comm;

      MPI_Status stat;

      MPI_Recv (_buffer, _length, MPI_CHAR, _rank, _msgtag, *(pcomm), &stat);
#endif
   }

// Description: Asynchronous receive data from the other CPU (implementation)
//========================================================================================
   void CMPIDataExchange::IRecv (void *_comm, int _rank, int _msgtag, int _length,
                                 char *_buffer, int _indrecv, void *_recv_arr)
   {
#ifdef USE_MPI
      MPI_Comm *pcomm = (MPI_Comm *) _comm;

      MPI_Request *p_recv_arr = (MPI_Request *) _recv_arr;

      MPI_Irecv (_buffer, _length, MPI_CHAR, _rank, _msgtag, *(pcomm),
                 p_recv_arr + _indrecv);
#endif
   }

// Description: Allocate send/receive requests (implementation)
//========================================================================================
   void CMPIDataExchange::AllocateRecvs (int _nrecv, void *&_recvarr)
   {
#ifdef USE_MPI
      MPI_Request *p_recv_arr = new MPI_Request[_nrecv];
      _recvarr = (void *) p_recv_arr;
#else
      char *p_recv_arr = new char[_nrecv];
      _recvarr = (void *) p_recv_arr;
#endif
   }

// Description: Delete send/receive requests (implementation)
//========================================================================================
   void CMPIDataExchange::DeleteRecvs (void *_recvarr)
   {
#ifdef USE_MPI
      MPI_Request *p_recv_arr = (MPI_Request *) _recvarr;
      delete[]p_recv_arr;
#else
      char *p_recv_arr = (char *) _recvarr;
      delete[]p_recv_arr;
#endif
   }

// Description: Allocate statuses (implementation)
//========================================================================================
   void CMPIDataExchange::AllocateStats (int _nstat, void *&_statarr)
   {
#ifdef USE_MPI
      MPI_Status *p_stat_arr = new MPI_Status[_nstat];
      _statarr = (void *) p_stat_arr;
#else
      char *p_stat_arr = new char[_nstat];
      _statarr = (void *) p_stat_arr;
#endif
   }

// Description: Delete statuses (implementation)
//========================================================================================
   void CMPIDataExchange::DeleteStats (void *_statarr)
   {
#ifdef USE_MPI
      MPI_Status *p_stat_arr = (MPI_Status *) _statarr;
      delete[]p_stat_arr;
#else
      char *p_stat_arr = (char *) _statarr;
      delete[]p_stat_arr;
#endif
   }

// Description: Wait for completion of exchanges (implementation)
//========================================================================================
   void CMPIDataExchange::WaitAll (int _count, void *_recvarr, void *_statarr)
   {
#ifdef USE_MPI
      MPI_Request *p_recv_arr = (MPI_Request *) _recvarr;
      MPI_Status *p_stat_arr = (MPI_Status *) _statarr;

      MPI_Waitall (_count, p_recv_arr, p_stat_arr);
#endif
   }

//========================================================================================
   void CMPIDataExchange::ExchangeArray (void *_comm, char _Datatype, char _Operation,
                                         int _Length, void *_Array)
   {
#ifdef USE_MPI
      static int maxbufsize = 2408 * 1024;

      MPI_Comm *pcomm = (MPI_Comm *) _comm;

      if (_Length < maxbufsize) {

         CMPIDataExchange::ExchangeArray_impl (_comm, _Datatype, _Operation, _Length,
                                               _Array);

      } else {

         int myidloc;
         MPI_Comm_rank (*pcomm, &myidloc);

         char *parr = (char *) _Array;

         int iscale = 0;

         if (_Datatype == 'C')
            iscale = sizeof (char);
         if (_Datatype == 'I')
            iscale = sizeof (int);
         if (_Datatype == 'L')
            iscale = sizeof (long long);
         if (_Datatype == 'F')
            iscale = sizeof (float);
         if (_Datatype == 'D')
            iscale = sizeof (double);

         int ibeg, iend, ni;

         iend = -1;
         while (iend < _Length - 1) {
            ibeg = iend + 1;
            iend = ibeg + maxbufsize - 1;
            if (iend > _Length - 1)
               iend = _Length - 1;
            ni = iend - ibeg + 1;
            CMPIDataExchange::ExchangeArray_impl (_comm, _Datatype, _Operation, ni,
                                                  (void *) (parr + ibeg * iscale));
         }

      }
#endif
   }

//========================================================================================
//   void MPIAPI OpCharAdd (void *in, void *inout, int *len, MPI_Datatype *)
   void OpCharAdd (void *in, void *inout, int *len, MPI_Datatype *)
   {
      char *inL = (char *) in;
      char *inoutL = (char *) inout;
      int i;
      char c;

      for (i = 0; i < *len; ++i) {
         c = *inoutL + *inL;
         *inoutL = c;
         inL++;
         inoutL++;
      }
   }

//========================================================================================
//   void MPIAPI OpCharMaximum (void *in, void *inout, int *len, MPI_Datatype * dptr)
   void OpCharMaximum (void *in, void *inout, int *len, MPI_Datatype * dptr)
   {
      char *inL = (char *) in;
      char *inoutL = (char *) inout;
      int i;
      char c;

      for (i = 0; i < *len; ++i) {
         c = *inoutL > *inL ? *inoutL : *inL;
         *inoutL = c;
         inL++;
         inoutL++;
      }
   }

//========================================================================================
//   void MPIAPI OpCharMinimum (void *in, void *inout, int *len, MPI_Datatype * dptr)
   void OpCharMinimum (void *in, void *inout, int *len, MPI_Datatype * dptr)
   {
      char *inL = (char *) in;
      char *inoutL = (char *) inout;
      int i;
      char c;

      for (i = 0; i < *len; ++i) {
         c = *inoutL < *inL ? *inoutL : *inL;
         *inoutL = c;
         inL++;
         inoutL++;
      }
   }

//========================================================================================
//   void MPIAPI OpLongLongAdd (void *in, void *inout, int *len, MPI_Datatype * dptr)
   void OpLongLongAdd (void *in, void *inout, int *len, MPI_Datatype * dptr)
   {
      long long *inL = (long long *) in;
      long long *inoutL = (long long *) inout;
      int i;
      long long c;

      for (i = 0; i < *len; ++i) {
         c = *inoutL + *inL;
         *inoutL = c;
         inL++;
         inoutL++;
      }
   }

// Collective array operation
//========================================================================================
   void CMPIDataExchange::ExchangeArray_impl (void *_comm, char _Datatype,
                                              char _Operation, int _Length, void *_Array)
   {
#ifdef USE_MPI
      MPI_Comm *pcomm = (MPI_Comm *) _comm;

      int nprocloc;
      MPI_Comm_size (*pcomm, &nprocloc);

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
         size_of_data = sizeof (char);

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
         size_of_data = sizeof (int);

      } else if (_Datatype == 'L') {

         datatype = MPI_LONG_LONG_INT;
         size_of_data = sizeof (long long);
         if (_Operation == '+') {
            MPI_Op_create (OpLongLongAdd, true, &op);
            op_create_flag = true;
         }

      } else if (_Datatype == 'F') {

         datatype = MPI_FLOAT;
         size_of_data = sizeof (float);

      } else if (_Datatype == 'D') {

         datatype = MPI_DOUBLE;
         size_of_data = sizeof (double);

      }

      if (nprocloc != 1) {
         void *buffer = new char[_Length * size_of_data];
         MPI_Allreduce (_Array, buffer, _Length, datatype, op, *(pcomm));
         memcpy (_Array, buffer, _Length * size_of_data);
         delete[](char *) buffer;
      }

      if (op_create_flag) {
         MPI_Op_free (&op);
      }
#endif
   }

// Major data exchange function
//========================================================================================
   void CMPIDataExchange::DataExchange (void *_comm, vector < int >&_CpuIDSend,
                                        vector < vector < char > >&_ObjSend,
                                        vector < int >&_CpuIDRecv,
                                        vector < vector < char > >&_ObjRecv)
   {
#ifdef USE_MPI

      int nproc = CMPIDataExchange::GetNproc (_comm);

      if (nproc == 1) {
         _CpuIDRecv.swap (_CpuIDSend);
         _ObjRecv.swap (_ObjSend);
      } else {
         CMPIDataExchange::DataExchange_impl (_comm, _CpuIDSend, _ObjSend, _CpuIDRecv,
                                              _ObjRecv);
      }
#else
      _CpuIDRecv.swap (_CpuIDSend);
      _ObjRecv.swap (_ObjSend);
#endif
   }

// Description: Packed data exchange (implementation)
//========================================================================================
   void CMPIDataExchange::DataExchange_impl (void *_comm, vector < int >&_CpuIDSend,
                                             vector < vector < char > >&_ObjSend,
                                             vector < int >&_CpuIDRecv,
                                             vector < vector < char > >&_ObjRecv)
   {
#ifdef USE_MPI

// Exchange the data

      int NObjSend = (int) _CpuIDSend.size ();

      vector < int >CpuIDSend (NObjSend);
      vector < int >ObjSizeSend (NObjSend);
      vector < char *>ObjSend (NObjSend);

      int *pCpuIDSend = NULL;
      int *pObjSizeSend = NULL;
      char **pObjSend = NULL;

      if (NObjSend > 0) {
         pCpuIDSend = &CpuIDSend[0];
         pObjSizeSend = &ObjSizeSend[0];
         pObjSend = &ObjSend[0];
      }

      int *p_data_send_cpuid = NULL;
      vector < char >*p_data_send_ch = NULL;

      if (NObjSend > 0) {
         p_data_send_cpuid = &_CpuIDSend[0];
         p_data_send_ch = &_ObjSend[0];
      }

      int i;

      for (i = 0; i < NObjSend; i++) {
         pCpuIDSend[i] = p_data_send_cpuid[i];
         pObjSizeSend[i] = (int) p_data_send_ch[i].size ();
         pObjSend[i] = &(p_data_send_ch[i][0]);
      }

      int NObjRecv;
      int *CpuIDRecv;
      int *ObjSizeRecv;
      char **ObjRecv;

      CMPIDataExchange::DataExchange_impl2 (_comm, (int) NObjSend, pCpuIDSend,
                                            pObjSizeSend, pObjSend, NObjRecv, CpuIDRecv,
                                            ObjSizeRecv, ObjRecv);

// Prepare reply

      _CpuIDRecv.resize (NObjRecv);
      _ObjRecv.resize (NObjRecv);

      int *p_data_recv_cpuid = NULL;
      vector < char >*p_data_recv_ch = NULL;

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
         delete[]ObjRecv[i];
      }

// Free receive structures

      delete[]CpuIDRecv;
      delete[]ObjSizeRecv;
      delete[]ObjRecv;
#endif
   }

#ifdef USE_MASSIVE_SENDRECV

// Data exchange function
//========================================================================================
   void CMPIDataExchange::DataExchange_impl2 (void *_comm, int _NObjSend, int *_CpuIDSend,
                                              int *_ObjSizeSend, char **_ObjSend,
                                              int &_NObjRecv, int *&_CpuIDRecv,
                                              int *&_ObjSizeRecv, char **&_ObjRecv)
   {
#ifdef USE_MPI

// Take the number of cpus and cpu id

      MPI_Comm *pcomm = (MPI_Comm *) _comm;

      int nproc;
      int myid;

      MPI_Comm_size (*pcomm, &nproc);
      MPI_Comm_rank (*pcomm, &myid);

// Exchange send control data and transform it into receive control data

      int *pind_send;
      int *pind_recv;

      CMPIDataExchange::DataExchangeRequest (_comm, _NObjSend, _CpuIDSend, _ObjSizeSend,
                                             _NObjRecv, _CpuIDRecv, _ObjSizeRecv,
                                             pind_send, pind_recv);

// Create structures that support sparse search

      int *iasend, *jasend;
      int *iarecv, *jarecv;
      int *iptr;

      iasend = new int[nproc + 1];
      jasend = new int[_NObjSend];
      iarecv = new int[nproc + 1];
      jarecv = new int[_NObjRecv];
      iptr = new int[nproc];

      int i, iproc;

      for (i = 0; i <= nproc; i++)
         iasend[i] = 0;

      for (i = 0; i < _NObjSend; i++) {
         iproc = _CpuIDSend[i];
         iasend[iproc + 1]++;
      }

      for (i = 0; i < nproc; i++)
         iasend[i + 1] = iasend[i] + iasend[i + 1];

      for (i = 0; i < nproc; i++)
         iptr[i] = iasend[i];

      int k;

      for (i = 0; i < _NObjSend; i++) {
         iproc = _CpuIDSend[i];
         k = iptr[iproc];
         jasend[k] = i;
         iptr[iproc]++;
      }

      for (i = 0; i <= nproc; i++)
         iarecv[i] = 0;

      for (i = 0; i < _NObjRecv; i++) {
         iproc = _CpuIDRecv[i];
         iarecv[iproc + 1]++;
      }

      for (i = 0; i < nproc; i++)
         iarecv[i + 1] = iarecv[i] + iarecv[i + 1];

      for (i = 0; i < nproc; i++)
         iptr[i] = iarecv[i];

      for (i = 0; i < _NObjRecv; i++) {
         iproc = _CpuIDRecv[i];
         k = iptr[iproc];
         jarecv[k] = i;
         iptr[iproc]++;
      }

// Create the set of send/receive MPI datatypes

      int *nobjsend, *nobjrecv;
      size_t *sizesendtot, *sizerecvtot;

      nobjsend = new int[nproc];
      nobjrecv = new int[nproc];
      sizesendtot = new size_t[nproc];
      sizerecvtot = new size_t[nproc];

      int j;

      for (iproc = 0; iproc < nproc; iproc++) {
         nobjsend[iproc] = 0;
         sizesendtot[iproc] = 0;
         for (j = iasend[iproc]; j < iasend[iproc + 1]; j++) {
            i = jasend[j];
            sizesendtot[iproc] += _ObjSizeSend[i];
            nobjsend[iproc]++;
         }
         nobjrecv[iproc] = 0;
         sizerecvtot[iproc] = 0;
         for (j = iarecv[iproc]; j < iarecv[iproc + 1]; j++) {
            i = jarecv[j];
            sizerecvtot[iproc] += _ObjSizeRecv[i];
            nobjrecv[iproc]++;
         }
      }

// Allocate receive data

      _ObjRecv = new char *[_NObjRecv];

      for (i = 0; i < _NObjRecv; i++) {
         _ObjRecv[i] = new char[_ObjSizeRecv[i]];
      }

// Perform local copy

      int ipsend, iprecv, nloc;

      if (nobjsend[myid] != 0) {
         ipsend = 0;
         iprecv = 0;
         nloc = 0;
         while (nloc != nobjsend[myid]) {
            while (_CpuIDSend[ipsend] != myid)
               ipsend++;
            while (_CpuIDRecv[iprecv] != myid)
               iprecv++;
            memcpy (_ObjRecv[iprecv], _ObjSend[ipsend],
                    _ObjSizeSend[ipsend] * sizeof (char));
            ipsend++;
            iprecv++;
            nloc++;
         }
      }
// Create the set of send/receive MPI datatypes

      size_t isizesend = 0;
      size_t isizerecv = 0;

      for (iproc = 0; iproc < nproc; iproc++) {
         if (iproc != myid) {
            if (sizesendtot[iproc] > isizesend)
               isizesend = sizesendtot[iproc];
            if (sizerecvtot[iproc] > isizerecv)
               isizerecv = sizerecvtot[iproc];
         }
      }

      size_t maxsendrecvbuf = 128 * 1024 * 1024;
//      size_t maxsendrecvbuf = 256 * 1024 * 1024;

      size_t isizemax = isizesend;
      if (isizerecv > isizemax)
         isizemax = isizerecv;

      size_t ncycle_send = 1;
      if (isizemax > maxsendrecvbuf)
         ncycle_send = (isizemax + maxsendrecvbuf - 1) / maxsendrecvbuf;

      int ncycle_send_int = (int) ncycle_send;

      CMPIDataExchange::ExchangeArray (_comm, 'I', 'M', 1, &ncycle_send_int);

      char *objsendbuf;
      char *objrecvbuf;

      objsendbuf = new char[isizesend];
      objrecvbuf = new char[isizerecv];

// Main exchange cycles

      int iproc_send, iproc_recv;

      MPI_Status stat;

      for (k = 1; k < nproc; k++) {

         iproc_send = (myid + k) % nproc;
         iproc_recv = (myid - k + nproc) % nproc;

         isizesend = 0;
         for (j = iasend[iproc_send]; j < iasend[iproc_send + 1]; j++) {
            i = jasend[j];
            memcpy (objsendbuf + isizesend, _ObjSend[i], _ObjSizeSend[i] * sizeof (char));
            isizesend += _ObjSizeSend[i];
         }
         isizerecv = 0;
         for (j = iarecv[iproc_recv]; j < iarecv[iproc_recv + 1]; j++) {
            i = jarecv[j];
            isizerecv += _ObjSizeRecv[i];
         }

//      if (isizesend > 0 || isizerecv > 0) {

         size_t iend_send = -1;
         size_t iend_recv = -1;

         for (int icycle = 0; icycle < ncycle_send_int; icycle++) {

            size_t ibeg_send = iend_send + 1;
            size_t ibeg_recv = iend_recv + 1;

            iend_send = ibeg_send + maxsendrecvbuf;
            iend_recv = ibeg_recv + maxsendrecvbuf;

            if (iend_send >= isizesend)
               iend_send = isizesend - 1;
            if (iend_recv >= isizerecv)
               iend_recv = isizerecv - 1;

            size_t ni_send = iend_send - ibeg_send + 1;
            size_t ni_recv = iend_recv - ibeg_recv + 1;

//            if (ni_send > 0 || ni_recv > 0) {

            MPI_Sendrecv (objsendbuf + ibeg_send, (int) (ni_send * sizeof (char)),
                          MPI_CHAR, iproc_send, iproc_send + icycle * nproc,
                          objrecvbuf + ibeg_recv, (int) (ni_recv * sizeof (char)),
                          MPI_CHAR, iproc_recv, myid + icycle * nproc, *(pcomm), &stat);

//            }

         }

//      }

         isizerecv = 0;
         for (j = iarecv[iproc_recv]; j < iarecv[iproc_recv + 1]; j++) {
            i = jarecv[j];
            memcpy (_ObjRecv[i], objrecvbuf + isizerecv, _ObjSizeRecv[i] * sizeof (char));
            isizerecv += _ObjSizeRecv[i];
         }

      }

// Free work arrays

      delete[]iasend;
      delete[]jasend;
      delete[]iarecv;
      delete[]jarecv;
      delete[]iptr;
      delete[]nobjsend;
      delete[]nobjrecv;
      delete[]sizesendtot;
      delete[]sizerecvtot;
      delete[]objsendbuf;
      delete[]objrecvbuf;
      delete[]indsend;
      delete[]indrecv;

#endif
   }

#else

// Data exchange function
//========================================================================================
   void CMPIDataExchange::DataExchange_impl2 (void *_comm, int _NObjSend, int *_CpuIDSend,
                                              int *_ObjSizeSend, char **_ObjSend,
                                              int &_NObjRecv, int *&_CpuIDRecv,
                                              int *&_ObjSizeRecv, char **&_ObjRecv)
   {
#ifdef USE_MPI

// Take the number of cpus and cpu id

      MPI_Comm *pcomm = (MPI_Comm *) _comm;

      int nproc;
      int myid;

      MPI_Comm_size (*pcomm, &nproc);
      MPI_Comm_rank (*pcomm, &myid);

// Exchange send control data and transform it into receive control data

      int *pind_send;
      int *pind_recv;

      CMPIDataExchange::DataExchangeRequest (_comm, _NObjSend, _CpuIDSend, _ObjSizeSend,
                                             _NObjRecv, _CpuIDRecv, _ObjSizeRecv,
                                             pind_send, pind_recv);

// Allocate receive data

      int i;

      _ObjRecv = new char *[_NObjRecv];

      for (i = 0; i < _NObjRecv; i++) {
         _ObjRecv[i] = new char[_ObjSizeRecv[i]];
      }

// Copy own data

      int ip_send = 0;
      int ip_recv = 0;

      int iproc, iproc1;

      while (ip_send < _NObjSend && ip_recv < _NObjRecv) {
         iproc = _CpuIDSend[ip_send];
         iproc1 = _CpuIDRecv[ip_recv];
         if (iproc != myid && iproc1 != myid) {
            ip_send++;
            ip_recv++;
         } else if (iproc != myid && iproc1 == myid) {
            ip_send++;
         } else if (iproc == myid && iproc1 != myid) {
            ip_recv++;
         } else {
            if (_ObjSizeSend[ip_send] != _ObjSizeRecv[ip_recv]) {
               cout <<
                  " CMPIDataExchange::DataExchange_impl2: error: incorrect size sended to itself !!! ";
               throw
                  " CMPIDataExchange::DataExchange_impl2: error: incorrect size sended to itself !!! ";
            } else {
               memcpy (_ObjRecv[ip_recv], _ObjSend[ip_send],
                       _ObjSizeSend[ip_send] * sizeof (char));
            }
            ip_send++;
            ip_recv++;
         }
      }

// Allocate reqs and stats

      void *precvarr;
      void *pstatarr;

      CMPIDataExchange::AllocateRecvs (_NObjSend + _NObjRecv, precvarr);
      CMPIDataExchange::AllocateStats (_NObjSend + _NObjRecv, pstatarr);

// Perform async sends and recvs

      int nsends = 0;

      int msgtag;

      for (i = 0; i < _NObjSend; i++) {
         iproc = _CpuIDSend[i];
         msgtag = pind_send[i];
         if (iproc != myid) {
            CMPIDataExchange::ISend (_comm, iproc, msgtag, _ObjSizeSend[i], _ObjSend[i],
                                     nsends, precvarr);
            nsends++;
         }
      }

      int nrecvs = 0;

      for (i = 0; i < _NObjRecv; i++) {
         iproc = _CpuIDRecv[i];
         msgtag = pind_recv[i];
         if (iproc != myid) {
            CMPIDataExchange::IRecv (_comm, iproc, msgtag, _ObjSizeRecv[i], _ObjRecv[i],
                                     nrecvs + nsends, precvarr);
            nrecvs++;
         }
      }

      CMPIDataExchange::WaitAll (nsends + nrecvs, precvarr, pstatarr);

      CMPIDataExchange::DeleteRecvs (precvarr);
      CMPIDataExchange::DeleteStats (pstatarr);

      delete[]pind_send;
      delete[]pind_recv;

#endif
   }

#endif

// Description: Exchange sizes data requests
//========================================================================================
   void CMPIDataExchange::DataExchangeRequest (void *_comm, int _NObjSend,
                                               int *_CpuIDSend, int *_ObjSizeSend,
                                               int &_NObjRecv, int *&_CpuIDRecv,
                                               int *&_ObjSizeRecv, int *&_IndSend,
                                               int *&_IndRecv)
   {
#ifdef USE_MPI

// Take the number of cpus and cpu id

      MPI_Comm *pcomm = (MPI_Comm *) _comm;

      int nproc;
      int myid;

      MPI_Comm_size (*pcomm, &nproc);
      MPI_Comm_rank (*pcomm, &myid);

// Collect the numbers of objects from all cpu's

      int *cpu2obj;

      cpu2obj = new int[nproc + 1];

      int i, iproc;

      for (i = 0; i <= nproc; i++)
         cpu2obj[i] = 0;

      cpu2obj[myid + 1] = _NObjSend;

      if (nproc > 1) {

         CMPIDataExchange::ExchangeArray (_comm, 'I', '+', nproc + 1, cpu2obj);

      }

      for (i = 0; i < nproc; i++)
         cpu2obj[i + 1] = cpu2obj[i] + cpu2obj[i + 1];

// Allocate and combine all the data into one collected set of arrays

      int nobjtot = cpu2obj[nproc];

      int *objsizetot;
      int *cpuidtot;

      objsizetot = new int[nobjtot];
      cpuidtot = new int[nobjtot];

      _IndSend = new int[_NObjSend];

      for (i = 0; i < nobjtot; i++) {

         cpuidtot[i] = 0;
         objsizetot[i] = 0;

      }

      int ibeg = cpu2obj[myid];

      for (i = cpu2obj[myid]; i < cpu2obj[myid + 1]; i++) {

         cpuidtot[i] = _CpuIDSend[i - ibeg];
         objsizetot[i] = _ObjSizeSend[i - ibeg];
         _IndSend[i - ibeg] = i;

      }

      if (nproc > 1) {

         CMPIDataExchange::ExchangeArray (_comm, 'I', '+', nobjtot, cpuidtot);
         CMPIDataExchange::ExchangeArray (_comm, 'I', '+', nobjtot, objsizetot);

      }
// Prepare the final data

      _NObjRecv = 0;

      for (i = 0; i < nobjtot; i++) {
         if (cpuidtot[i] == myid)
            _NObjRecv++;
      }

      _ObjSizeRecv = new int[_NObjRecv];
      _CpuIDRecv = new int[_NObjRecv];
      _IndRecv = new int[_NObjRecv];

      _NObjRecv = 0;

      for (iproc = 0; iproc < nproc; iproc++) {
         for (i = cpu2obj[iproc]; i < cpu2obj[iproc + 1]; i++) {
            if (cpuidtot[i] == myid) {
               _ObjSizeRecv[_NObjRecv] = objsizetot[i];
               _CpuIDRecv[_NObjRecv] = iproc;
               _IndRecv[_NObjRecv] = i;
               _NObjRecv++;
            }
         }
      }

// Free work arrays

      delete[]cpu2obj;
      delete[]objsizetot;
      delete[]cpuidtot;

#endif
   }

// Description: Get MPI wall time
//========================================================================================
   double CMPIDataExchange::GetWallTimeMPI ()
   {
#ifdef USE_MPI
      return MPI_Wtime ();
#else
      clock_t time1;
      double wtime1;

      time1 = clock ();
      wtime1 = (double) (time1) / (double) CLOCKS_PER_SEC;

      return wtime1;
#endif
   }

}                               // namespace k3d

#endif
