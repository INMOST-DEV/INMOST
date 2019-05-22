//------------------------------------------------------------------------------------------------
// File: k3d_slv_thr.hxx
//------------------------------------------------------------------------------------------------

#ifndef __K3D_SLV_THR_HXX
#define __K3D_SLV_THR_HXX

#include "k3d_slv.hxx"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

namespace k3d
{

   template < typename _Int, typename _Flt > class CFctThreads  // Class that supports factorization computations with threads
   {
    public:
// Functions that support hmatrix structure as a whole
      /// Perform ILU2 point factorization of the hmatrix with dinamic ordering and future diagonal modifications
      static void Ilu2HMatrix (bool _b_blk_wells, SParams & _params, int _nblks,
                               int _nblks1, int _nblks2, long long *_blks,
                               long long *_nzord_blks, CTree & _tree0, CTree & _tree1,
                               CBMatrix < _Int, _Flt > &_a_matr, CBMatrix < _Int,
                               _Flt > &_l_matr, CBMatrix < _Int, _Flt > &_u_matr,
                               CTree & _tree2_new, int &_nblks_new,
                               vector < long long >&_blks_new,
                               CVectorData < int >&_ordernew, double &_sclmin_att,
                               double &_sclmax_att, int &_nmodif, double &_eigmin_att,
                               double &_eigmax_att);
// Implementation functions
      /// Compute explicit scaling
      static void ComputeExplicitScaling (int _sctype, int _nitersc, double _sclmin,
                                          int _nblks, long long *_blks, CBMatrix < _Int,
                                          _Flt > &_a_matr, CVectorData < _Flt > &_scl_L,
                                          CVectorData < _Flt > &_scl_U,
                                          double &_sclmin_att, double &_sclmax_att);
      /// Perform explicit scaling
      static void HMatrixScale (int _nblks, long long *_blks, CBMatrix < _Int,
                                _Flt > &_a_matr, _Flt * _sclL, _Flt * _sclU);
      /// Split LU and combine pairs
      static void SplitLUPair (int _nblks, long long *_blks, CBMatrix < _Int,
                               _Flt > &_a_matr, CBMatrix < _Int, _Flt > &_a_pair);
      /// Split pair into L and U
      static void SplitPair (int _nblks, long long *_blks, CBMatrix < _Int,
                             _Flt > &_a_pair, CBMatrix < _Int, _Flt > &_a_L,
                             CBMatrix < _Int, _Flt > &_a_U);
      /// Perform ILU2 point factorization of the hmatrix with dynamic ordering and future diagonal modifications
      static void Ilu2BlockIlu2 (bool _b_blk_wells, ofstream * _pfout, int _index_out,
                                 int _ordlevel, int _ordtype, int _strtype, int _fcttype,
                                 int _fcttype_sch, double _pivmin, double _tau1,
                                 double _tau2, double _tau2_sch, double _theta,
                                 int _nblks, int _nblks1, int _nblks2, long long *_blks,
                                 long long *_nzord_blks, CTree & _tree0, CTree & _tree1,
                                 CBMatrix < _Int, _Flt > &_alu_pair, CTree & _tree2_new,
                                 int &_nblks_new, vector < long long >&_blks_new,
                                 CVectorData < int >&_ordernew, CBMatrix < _Int,
                                 _Flt > &_lu_pair, int &_nmodif, double &_eigmin_att,
                                 double &_eigmax_att);
      /// Perform symbolic factorization of the hmatrix with dynamic ordering and subtrees for schur blocks
      static void SymbolicFct (ofstream * _pfout, SParams * _params, CTree & _tree,
                               int _nblks, long long *_blks, CBMatrix < _Int,
                               _Flt > &_alu_sp, CVectorData < CTree > &_subtree_arr,
                               CVectorData < int >&_nblks_arr,
                               CVectorData < vector < int > >&_blks_arr,
                               CVectorData < int >&_ordernew, CBMatrix < _Int,
                               _Flt > &_lu_sp);
      /// Compute ND order, create binary tree, find separators and condense tree, and compute subtree's for nodes
      static void OrderNDSepSubTree (int _nlev, int _nlev_split, int _ordtype,
                                     int _fcttype, CMatrix < _Int, _Flt > &_amatr_sp,
                                     vector < int >&_order, CTree & _tree, int &_nhblks,
                                     vector < int >&_hblks,
                                     CVectorData < CTree > &_subtree_arr,
                                     vector < int >&_hblk2blks, vector < int >&_blks);
      /// Rescale U
      static void RescaleU (int _nblks, long long *_blks, CBMatrix < _Int,
                            _Flt > &_U_matr, _Flt * _sclU, _Flt * _inv_sclU);
      /// Set diagonal by zeroes
      static void SetZeroDiag (int _nblks, long long *_blks, CBMatrix < _Int,
                               _Flt > &_A_matr);
      /// Perform filtering of the pairs block sparsity according to the tree with diagonal modifications
      static void PairsFilterTreeModif (CTree & _tree, double _theta, int _nblks,
                                        long long *_blks, CBMatrix < _Int,
                                        _Flt > &_AU_matr);
      /// Compute reordered by columns rectangular hmatrix
      static void ReorderHMatrixColumnsPairs (int _nblksR, long long *_blksR,
                                              int _nblksC_ini, long long *_blksC_ini,
                                              CBMatrix < _Int, _Flt > &_hmatr_ini,
                                              int *_order, int _nblksC_fin,
                                              long long *_blksC_fin, CBMatrix < _Int,
                                              _Flt > &_hmatr_fin);
   };

   /// Base class that supports combined parallel multiplication computations (MPI+thr version)
   template < typename _Int, typename _Flt > class CMvmParThreads_base {
    protected:
      void *pcomm;              ///< Pointer to MPI communicator
      int nhblks;               ///< Number of hblocks
      int *phblk2cpu;           ///< Pointer to hblock to cpu distribution
      int *phblk2blks;          ///< Pointer to hblock to block distribution
      int *pblk2hblks;          ///< Pointer to block to hblocks distribution
      long long *phblks;        ///< Pointer to hblocks partitioning
      int nblks;                ///< Number of blocks
      long long *pblks;         ///< Pointer to blocks partitioning
      CBMatrix < _Int, _Flt > *phmatr;  ///< Pointer to hyperblocks array
      int ni_cpu;               ///< The size of the local part of vector data
      vector < int >ibsblk;     ///< Base adresses of local parts of the vector data
      int nlisthblk_own;        ///< The number of local hblocks
      vector < int >listhblk_own;       ///< The list of own hblocks
      int nlistblk_own;         ///< The number of local blocks
      vector < int >listblk_own;        ///< The list of own blocks
      int nsends;               ///< The number of sends for computations
      vector < int >snd2cpu;    ///< The list of cpus to be sended to for computations
      vector < int >ia_sends;   ///< The sizes of sends data for computations
      vector < int >ind_sends;  ///< The indices of local data for sends
      int nrecvs;               ///< The number of recvs for computations
      vector < int >rcv2cpu;    ///< The list of cpus to be received data from for computations
      vector < int >ia_recvs;   ///< The sizes of recvs for computations
      vector < int >ja_recvs;   ///< The pairs of local data that describe recvs
      int nblks_recvs;          ///< The number of received block parts
      vector < int >listblk_recvs;      ///< The list of block numbers of received blocks
      vector < int >ialist_recvs;       ///< The numbers of block rows contatining the block
      vector < int >jatriples_recvs;    ///< The complete list (pair) and index numbers inside of block rows containing each of the received blocks
      vector < int >iablk_recvs;        ///< The sizes of copy data for each received block
      vector < int >ind_recvs;  ///< The list of indices inside block for received blocks
      int i_size_x_send;        ///< Size of work send array
      int i_size_x_recv;        ///< Size of work recv array
      int i_size_x_temp;        ///< Size of work temp array
// Constructors and destructor
    private:
      /// Copy constructor
      CMvmParThreads_base (const CMvmParThreads_base < _Int, _Flt > &_aa)
      {
      };
      /// Copy operator
      CMvmParThreads_base < _Int, _Flt > &operator= (const CMvmParThreads_base < _Int,
                                                     _Flt > &_aa)
      {
         return *this;
      };
    public:
      /// Epmty constructor
      CMvmParThreads_base () {
         pcomm = NULL;
         nhblks = 0;
         phblk2cpu = NULL;
         phblk2blks = NULL;
         phblks = NULL;
         nblks = 0;
         pblks = NULL;
         phmatr = NULL;
         ni_cpu = 0;
         nlisthblk_own = 0;
         nlistblk_own = 0;
         nsends = 0;
         nrecvs = 0;
         nblks_recvs = 0;
         i_size_x_send = 0;
         i_size_x_recv = 0;
         i_size_x_temp = 0;
      };
      /// Destructor
      ~CMvmParThreads_base () {
      };
// External functions
      /// Init control data
      void InitControl_base (void *_pcomm, int _nhblks, int *_hblk2cpu, int *_hblk2blks,
                             int *_blk2hblks, long long *_hblks, int _nblks,
                             long long *_blks)
      {
         pcomm = _pcomm;
         int myid = CMPIDataExchange::GetMyid (pcomm);
         nhblks = _nhblks;
         phblk2cpu = _hblk2cpu;
         phblk2blks = _hblk2blks;
         pblk2hblks = _blk2hblks;
         phblks = _hblks;
         nblks = _nblks;
         pblks = _blks;
         ibsblk.resize (nblks + 1);
         int *pibsblk = &ibsblk[0];
         ni_cpu = 0;
         nlisthblk_own = 0;
         nlistblk_own = 0;
         int i, j;
         for (i = 0; i < nhblks; i++) {
            if (phblk2cpu[i] == myid) {
               nlisthblk_own++;
               for (j = phblk2blks[i]; j < phblk2blks[i + 1]; j++) {
                  nlistblk_own++;
                  pibsblk[j] = ni_cpu;
                  ni_cpu += (int) (pblks[j + 1] - pblks[j]);
               }
            } else {
               for (j = phblk2blks[i]; j < phblk2blks[i + 1]; j++)
                  pibsblk[j] = -1;
            }
         }
         listhblk_own.resize (nlisthblk_own + 1);
         int *plisthblk_own = &listhblk_own[0];
         nlisthblk_own = 0;
         listblk_own.resize (nlistblk_own + 1);
         int *plistblk_own = &listblk_own[0];
         nlistblk_own = 0;
         for (i = 0; i < nhblks; i++) {
            if (phblk2cpu[i] == myid) {
               plisthblk_own[nlisthblk_own] = i;
               nlisthblk_own++;
               for (j = phblk2blks[i]; j < phblk2blks[i + 1]; j++) {
                  plistblk_own[nlistblk_own] = j;
                  nlistblk_own++;
               }
            }
         }
      };
      /// Init MvmA data
      void InitMvmA_base (CBMatrix < _Int, _Flt > *_hmatr_arr);
      /// Clean MvmA structure
      void Clean_base ();
   };

   /// Class that supports combined parallel multiplication computations (MPI+thr version)
 template < typename _Int, typename _Flt, typename _FltVect > class CMvmParThreads:public CMvmParThreads_base < _Int,
      _Flt >
   {
    protected:
      CVectorData < _FltVect > x_send;  ///< The work memory for sends
      CVectorData < _FltVect > x_recv;  ///< The work memory for recvs
      CVectorData < _FltVect > x_temp;  ///< The work memory for external blocks
// Constructors and destructor
    private:
      /// Copy constructor
      CMvmParThreads (const CMvmParThreads < _Int, _Flt, _FltVect > &_aa)
      {
      };
      /// Copy operator
      CMvmParThreads < _Int, _Flt, _FltVect > &operator= (const CMvmParThreads < _Int,
                                                          _Flt, _FltVect > &_aa)
      {
         return *this;
      };
    public:
      /// Epmty constructor
      CMvmParThreads () {
      };
      /// Destructor
      ~CMvmParThreads () {
      };
// External functions
      /// Init control data
      void InitControl (void *_pcomm, int _nhblks, int *_hblk2cpu, int *_hblk2blks,
                        int *_blk2hblks, long long *_hblks, int _nblks, long long *_blks)
      {
         this->InitControl_base (_pcomm, _nhblks, _hblk2cpu, _hblk2blks, _blk2hblks,
                                 _hblks, _nblks, _blks);
      };
      /// Init MvmA data
      void InitMvmA (CBMatrix < _Int, _Flt > *_hmatr_arr);
      /// Perform MvmA computations
      void MvmA (const _FltVect * _x, _FltVect * _ax);
      /// Clean MvmA structure
      void Clean ();
   };

   template < typename _Int, typename _Flt, typename _FltVect > class CSlvParThreads    // Class that supports parallel solve computations (including L, U and ordering) (MPI+thr version)
   {
    protected:
      void *pcomm;              ///< Pointer to MPI communicator
      int nhblks;               ///< Number of hblocks
      int *phblk2cpu;           ///< Pointer to hblock to cpu distribution
      long long *phblks;        ///< Pointer to hblocks partitioning
      long long *phblks_ext;    ///< Pointer to extended hblocks partitioning
      int *phblk2blks;          ///< Pointer to hblock to blks partitioning
      int *pblk2hblks;          ///< Pointer to block to hblks partitioning
      long long *pblks;         ///< Pointer to blocks partitioning
      int nlisthblk_own;        ///< The number of local hblocks
      vector < int >listhblk_own;       ///< The list of own hblocks
      int ni_cpu;               ///< The size of the local part of vector data
      vector < int >ibsblk;     ///< Base adresses of local parts of the vector data
      vector < int >*plistpairs_ext;    ///< Pointer to the extended data lists
      int *pnblks_ext_arr;      ///< Pointer to the initial number of local blocks in extension
      vector < int >*pblksnum_ext_arr;  ///< Pointer to the initial block numbers of local hblocks in extension
      vector < long long >*pblks_ext_arr;       ///< Pointer to the initial blocks partitioning of local hblocks in extension
      CTree *ptree_arr;         ///< Pointer to the set of subtrees for fct
      int *pnblks_ilu2_arr;     ///< Pointer to the total number of local blocks in fct
      int *pnblks1_ilu2_arr;    ///< Pointer to the total number of local first  set of blocks in fct
      int *pnblks2_ilu2_arr;    ///< Pointer to the total number of local second set of blocks in fct
      vector < long long >*pblks_ilu2_arr;      ///< Pointer to the blocks partitioning of local hblocks in fct
      CVectorData < int >*porderLU;     ///< Pointer to the extended ordering data
      CBMatrix < _Int, _Flt > *pmatrL;  ///< Pointer to L blocks array
      CBMatrix < _Int, _Flt > *pmatrU;  ///< Pointer to U blocks array
      vector < CMatrix < int, float > >strL_T;   ///< Set of transposed block sparsities that support solve with L with parallelization
      int nsends;               ///< The number of sends for computations
      vector < int >snd2cpu;    ///< The list of cpus to be sended to for computations
      vector < int >ia_sends;   ///< The sizes of sends data for computations
      vector < int >ind_sends;  ///< The indices of local data for sends
      CVectorData < _FltVect > x_send;  ///< The work memory for sends
      int nrecvs;               ///< The number of recvs for computations
      vector < int >rcv2cpu;    ///< The list of cpus to be received data from for computations
      vector < int >ia_recvs;   ///< The sizes of recvs for computations
      CVectorData < _FltVect > x_recv;  ///< The work memory for recvs
      vector < int >ialist_recvs;       ///< For each local block the number of updates from received data
      vector < int >jalist_recvs;       ///< The lists of index numbers for received data for updates in current local block
      CVectorData < _FltVect > x_temp;  ///< The work memory for extended blocks
// Constructors and destructor
    private:
      /// Copy constructor
      CSlvParThreads (const CSlvParThreads < _Int, _Flt, _FltVect > &_aa)
      {
      };
      /// Copy operator
      CSlvParThreads < _Int, _Flt, _FltVect > &operator= (const CSlvParThreads < _Int,
                                                          _Flt, _FltVect > &_aa)
      {
         return *this;
      };
    public:
      /// Epmty constructor
      CSlvParThreads () {
         pcomm = NULL;
         nhblks = 0;
         phblks = NULL;
         phblks_ext = NULL;
         phblk2cpu = NULL;
         phblk2blks = NULL;
         pblk2hblks = NULL;
         pblks = NULL;
         ni_cpu = 0;
         nlisthblk_own = 0;
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
         nrecvs = 0;
      };
      /// Destructor
      ~CSlvParThreads () {
      };
// External functions
      /// Get ni_cpu
      virtual int GetNiCpu ()
      {
         return ni_cpu;
      };
      /// Get ni_cpu
      int GetNiCpuBase ()
      {
         return ni_cpu;
      };
      /// Init control data
      void InitControl (void *_pcomm, int _nhblks, long long *_hblks, int *_hblk2cpu,
                        int *_hblk2blks, int *_blk2hblks, long long *_blks)
      {
         pcomm = _pcomm;
         nhblks = _nhblks;
         phblks = _hblks;
         phblk2cpu = _hblk2cpu;
         phblk2blks = _hblk2blks;
         pblk2hblks = _blk2hblks;
         pblks = _blks;
         int myid = CMPIDataExchange::GetMyid (pcomm);
         int nblks = phblk2blks[nhblks];
         ibsblk.resize (nblks + 1);
         int *pibsblk = &ibsblk[0];
         ni_cpu = 0;
         nlisthblk_own = 0;
         int i, j;
         for (i = 0; i < nhblks; i++) {
            if (phblk2cpu[i] == myid) {
               nlisthblk_own++;
               for (j = _hblk2blks[i]; j < _hblk2blks[i + 1]; j++) {
                  pibsblk[j] = ni_cpu;
                  ni_cpu += (int) (pblks[j + 1] - pblks[j]);
               }
            } else {
               for (j = _hblk2blks[i]; j < _hblk2blks[i + 1]; j++)
                  pibsblk[j] = -1;
            }
         }
         listhblk_own.resize (nlisthblk_own + 1);
         int *plisthblk_own = &listhblk_own[0];
         nlisthblk_own = 0;
         for (i = 0; i < nhblks; i++) {
            if (phblk2cpu[i] == myid) {
               plisthblk_own[nlisthblk_own] = i;
               nlisthblk_own++;
            }
         };
      };
      /// Init SolveLU data
      void InitSolveLU (long long *_hblks_ext, vector < int >*_listpairs_ext,
                        int *_nblks_ext_arr, vector < int >*_blksnum_ext_arr,
                        vector < long long >*_blks_ext_arr, CTree * _tree_arr,
                        int *_nblks_ilu2_arr, int *_nblks1_ilu2_arr,
                        int *_nblks2_ilu2_arr, vector < long long >*_blks_ilu2_arr,
                        CVectorData < int >*_orderLU, CBMatrix < _Int, _Flt > *_matrL,
                        CBMatrix < _Int, _Flt > *_matrU);
      /// Perform SolveLU computations
      void SolveLU (const _FltVect * _x, _FltVect * _px);
      /// In-place perform SolveU computations
      void SolveU (int _ihblk, _FltVect * _px);
      /// In-place perform SolveL computations
      void SolveL (int _ihblk, _FltVect * _px);
      /// Clean Slv structure
      void Clean ();
   };

   /// Mvm and Slv interface functions
   template < typename _FltVect > struct TMvmSlvFunc
   {
      typedef void (*mvmA_f) (void *_obj_mvm, const _FltVect * _x, _FltVect * _ax);
      typedef void (*slvLU_f) (void *_obj_slv, const _FltVect * _x, _FltVect * _px);
   };

   /// Structure that contains indices for Fast Transformation (FT)
   struct SIndFT
   {
      short ihblk_FT;           ///< Hblock number
      short ind_blk_FT;         ///< Index of block in hblock
      int ind_elem_FT;          ///< Index of element in a block
   };

   template < typename _Int, typename _Flt, typename _FltVect > class CK3D_SolverThreads        // Class that supports combined parallel solver computations (MPI+thr)
   {
    protected:
      void *pcomm;              ///< Pointer to MPI communicator
      SParams params;           ///< Solver params
      bool b_fast_transform;    ///< True if fast transformation is set and prepared
      bool b_use_ini;           ///< Use or not initial data decomposition
      bool b_use_wells;         ///< Use or not well decomposition
      bool b_use_blksize;       ///< Use or not blksize based ordering and block scaling
      bool b_blk_wells;         ///< True if well partitioning have wells nodes itself
      int nhblks_ini;           ///< Initial number of hblocks
      int nhblks;               ///< Number of hblocks
      int nblks_ini;            ///< Initial number of blocks
      int nblks;                ///< Number of blocks
      vector < long long >hblks_ini;    ///< Initial hblocks partitioning
      vector < int >hblk2cpu_ini;       ///< Initial block to cpu distribution
      vector < int >hblk2blks_ini;      ///< Initial hblocks to blocks partitioning
      vector < long long >blks_ini;     ///< Initial blocks partitioning
      vector < CVectorData < int > >order_ini;   ///< Ordering of initial data
      CVectorData < long long >order_wells;     ///< Ordering of wells data
      CVectorData < int >order2ind_wells;       ///< 2index (hblk,ind) ordering of wells data
      vector < long long >hblks;        ///< Hblocks partitioning
      vector < int >hblk2cpu;   ///< Pointer to tblock to cpu distribution
      vector < int >hblk2blks;  ///< Hblocks to blocks partitioning
      vector < int >blk2hblks;  ///< Blocks to hblocks partitioning
      vector < long long >blks; ///< Blocks partitioning
      vector < CBMatrix < _Int, _Flt > >hmatr_arr;      ///< Matrix data
      vector < CBMatrix < _Int, _Flt > >hmatr_arr_save; ///< Temporary store of matrix data (temporary save before explicit scaling)
      vector < CBMatrix < _Int, double > >hmatr_ini_FT_arr;      ///< Fast transform initial matrix data
      vector < CBMatrix < _Int, double > >hmatr_fin_FT_arr;      ///< Fast transform final matrix data
      CMvmParThreads < _Int, _Flt, _FltVect > mvm;      ///< Support MvmA
      int blksize_bscl;         ///< Block scaling parameter
      vector < CBMatrix < _Int, _Flt > >sclL;   ///< Left  blksize scaling (explicit block scaling also)
      vector < CBMatrix < _Int, _Flt > >sclU;   ///< Right blksize scaling (explicit block scaling also)
      vector < CBMatrix < _Int, _Flt > >sclLInv;        ///< Inverse left  blksize scaling (explicit block scaling also)
      vector < CBMatrix < _Int, _Flt > >sclUInv;        ///< Inverse right blksize scaling (explicit block scaling also)
      vector < long long >hblks_ext;    ///< Hblocks partitioning
      vector < int >nlist_ext_arr;      ///< The numbers of elems in extension
      vector < vector < int > >list_ext_arr;     ///< The pairs of extension indices
      vector < int >nblks_ext_arr;      ///< The initial number of local blocks in extension
      vector < vector < int > >blksnum_ext_arr;  ///< The initial block numbers of local hblocks in extension
      vector < vector < long long > >blks_ext_arr;       ///< The initial blocks partitioning of local hblocks in extension
      vector < CTree > tree_arr;        ///< The set of subtrees for fct
      vector < int >nblks_ilu2_arr;     ///< The total number of local blocks in fct
      vector < int >nblks1_ilu2_arr;    ///< The total number of local first  set of blocks in fct
      vector < int >nblks2_ilu2_arr;    ///< The total number of local second set of blocks in fct
      vector < vector < long long > >blks_ilu2_arr;      ///< The blocks partitioning of local hblocks in fct
      vector < vector < long long > >nzord_ilu2_arr;     ///< Parts of diagonal blocks that could be ordered in fct
      vector < CVectorData < int > >order_LU;    ///< Ordering data
      vector < CBMatrix < _Int, float > >matrL_float;    ///< Float  L part
      vector < CBMatrix < _Int, double > >matrL_double;  ///< Double L part
      vector < CBMatrix < _Int, float > >matrU_float;    ///< Float  U part
      vector < CBMatrix < _Int, double > >matrU_double;  ///< Double U part
      CSlvParThreads < _Int, float, _FltVect > slv_float;       ///< Support float  Slv
      CSlvParThreads < _Int, double, _FltVect > slv_double;     ///< Support double Slv
      int ncoef_slv;            ///< Number of coefs in poly prec
      vector < double >coef_slv;        ///< Coefs in poly prec
      CVectorData < _FltVect > xwork;   ///< Work memory for polynomial preconditioning
      CVectorData < _FltVect > xwork_bscl;      ///< Work memory for block scaling computations
      vector < CBMatrix < _Int, _Flt > >hmatr_arr_temp; ///< Temporary storage for matrix data
// Constructors and destructor
    private:
      /// Copy constructor
      CK3D_SolverThreads < _Int, _Flt, _FltVect > (const CK3D_SolverThreads < _Int, _Flt,
                                                   _FltVect > &_aa) {
      }
      /// Copy operator
      CK3D_SolverThreads < _Int, _Flt,
         _FltVect > &operator= (const CK3D_SolverThreads < _Int, _Flt, _FltVect > &_aa)
      {
         return *this;
      }
    public:
      /// Epmty constructor
      CK3D_SolverThreads () {
         b_use_ini = false;
         b_use_wells = false;
         b_use_blksize = false;
         b_blk_wells = false;
         b_fast_transform = false;
         nhblks_ini = 0;
         nhblks = 0;
         nblks_ini = 0;
         nblks = 0;
         blksize_bscl = -1;
         ncoef_slv = 0;
      }
      /// Destructor
      ~CK3D_SolverThreads () {
      }
// Get/set functions
      /// Get comm
      void *GetComm ()
      {
         return pcomm;
      }
      /// Get nhblks
      int GetNhblks ()
      {
         return nhblks;
      }
      /// Get nhblks_ini
      int GetNhblksIni ()
      {
         return nhblks_ini;
      }
      /// Get nblks
      int GetNblks ()
      {
         return nblks;
      }
      /// Get tblks
      long long *GetHBlks ()
      {
         return &hblks[0];
      }
      /// Get tblks
      long long *GetHBlksIni ()
      {
         return &hblks_ini[0];
      }
      /// Get tblk2cpu
      int *GetHBlk2cpu ()
      {
         return &hblk2cpu[0];
      }
      /// Get tblk2cpu_ini
      int *GetHBlk2cpuIni ()
      {
         return &hblk2cpu_ini[0];
      }
      /// Get tblk2blks
      int *GetHBlk2blks ()
      {
         return &hblk2blks[0];
      }
      /// Get blk2hblks
      int *GetBlk2hblks ()
      {
         return &blk2hblks[0];
      }
      /// Get blks
      long long *GetBlks ()
      {
         return &blks[0];
      }
      /// Get order_wells
      long long *GetOrderWells ()
      {
         return order_wells.Ptr ();
      }
      /// Set ncoef_slv
      void SetNcoef (int _ncoef)
      {
         ncoef_slv = _ncoef;
      }
      /// Get coef
      vector < double >*GetCoef ()
      {
         return &coef_slv;
      }
      /// Get local cpu vector size, current (with vtype == 1) or initial (with vtype == 0)
      long long GetVSize (int _vtype);
      /// Get ni
      virtual int GetNi ()
      {
         int ni_local = 0;
         if (this->params.prec_float == 1) {
            ni_local = slv_float.GetNiCpu ();
         } else {
            ni_local = slv_double.GetNiCpu ();
         }
         return ni_local;
      };
      /// Get matrix or scaling data data 
      CBMatrix < _Int, _Flt > *GetHMatr (char _type) {
         CBMatrix < _Int, _Flt > *phmatr_arr = NULL;
         if (_type == 'A')
            phmatr_arr = this->hmatr_arr.data ();
         if (_type == 'a')
            phmatr_arr = this->hmatr_arr_save.data ();
         if (_type == 'L')
            phmatr_arr = this->sclL.data ();
         if (_type == 'l')
            phmatr_arr = this->sclLInv.data ();
         if (_type == 'U')
            phmatr_arr = this->sclU.data ();
         if (_type == 'u')
            phmatr_arr = this->sclUInv.data ();
         return phmatr_arr;
      }
// Aggregated interfaces
      /// Prepare matrix structures
      void PrepareMatrix (void *_pcomm, SParams * _params, SStatData * _stats,
                          int _nhblks, int *_hblk2cpu, int *_hblk2blks, long long *_hblks,
                          int _nblks, long long *_blks, _Int * _ia, _Int * _ja, _Flt * _a)
      {
         if (_params->i_decomp_type == 0) {
            this->PrepareMatrixBase (_pcomm, _nhblks, _hblk2cpu, _hblk2blks, _hblks,
                                     _nblks, _blks, _ia, _ja, _a);
         } else if (_params->i_decomp_type == 1) {
            this->PrepareMatrixThreads (_pcomm, _params, _stats, _nhblks, _hblk2cpu,
                                        _hblk2blks, _hblks, _nblks, _blks, _ia, _ja, _a);
         } else if (_params->i_decomp_type == 2) {
            this->PrepareMatrixWells (_pcomm, _params, _stats, _nhblks, _hblk2cpu,
                                      _hblk2blks, _hblks, _nblks, _blks, NULL, NULL, _ia,
                                      _ja, _a);
         } else if (_params->i_decomp_type > 2) {
            this->PrepareMatrixBlksize (_pcomm, _params, _stats, _nhblks, _hblk2cpu,
                                        _hblk2blks, _hblks, _nblks, _blks, _ia, _ja, _a);
         }
      }
      /// Perform fct and free matrix if necessary
      void ComputePreconditioner (SParams * _params, SStatData * _stats,
                                  bool _b_store_matrix)
      {
         this->ComputeBILU2 (_params, _stats);
         if (!_b_store_matrix)
            this->CleanMvmA ();
      }
      /// Prepare solver structures including performing parallel fct
      void PrepareSolver (void *_pcomm, SParams * _params, SStatData * _stats,
                          int _nhblks, int *_hblk2cpu, int *_hblk2blks, long long *_hblks,
                          int _nblks, long long *_blks, bool _b_store_matrix, _Int * _ia,
                          _Int * _ja, _Flt * _a)
      {
         this->PrepareMatrix (_pcomm, _params, _stats, _nhblks, _hblk2cpu, _hblk2blks,
                              _hblks, _nblks, _blks, _ia, _ja, _a);
         this->ComputePreconditioner (_params, _stats, _b_store_matrix);
      }
      /// Perform iterations of the iterative scheme
      void SolveIter (SParams * _params, SStatData * _stats, _FltVect * _rhs,
                      _FltVect * _sol);
      /// Perform iterations of the iterative scheme (implementation)
      void SolveIter (SParams * _params, SStatData * _stats, void *_str_mvmA,
                      typename TMvmSlvFunc < _FltVect >::mvmA_f _mvm_f, void *_str_slvLU,
                      typename TMvmSlvFunc < _FltVect >::slvLU_f _slv_f, _FltVect * _rhs,
                      _FltVect * _sol);
      /// Transform vectors forward
      void TransformVectorsForward (_FltVect * _rhs_ini, _FltVect * _sol_ini,
                                    vector < _FltVect > &_rhs_new,
                                    vector < _FltVect > &_sol_new);
      /// Transform vectors backward
      void TransformVectorsBackward (vector < _FltVect > &_rhs_new,
                                     vector < _FltVect > &_sol_new, _FltVect * _rhs_fin,
                                     _FltVect * _sol_fin);
// External functions
      /// Store control data only
      void StoreControlData (void *_pcomm, int _nhblks, int *_hblk2cpu, int *_hblk2blks,
                             long long *_hblks, int _nblks, long long *_blks);
      /// Prepare matrix data only
      void PrepareMatrixBase (void *_pcomm, int _nhblks, int *_hblk2cpu, int *_hblk2blks,
                              long long *_hblks, int _nblks, long long *_blks, _Int * _ia,
                              _Int * _ja, _Flt * _a);
      /// Prepare matrix data only including preliminary repartitioning of the matrix
      void PrepareMatrixThreads (void *_pcomm, SParams * _params, SStatData * _stats,
                                 int _nhblks, int *_hblk2cpu, int *_hblk2blks_ini,
                                 long long *_hblks, int _nblks_ini, long long *_blks_ini,
                                 _Int * _ia, _Int * _ja, _Flt * _a);
      /// Prepare matrix data only including preliminary repartitioning of the matrix taking into account wells data
      void PrepareMatrixWells (void *_pcomm, SParams * _params, SStatData * _stats,
                               int _nhblks_ini, int *_hblk2cpu_ini, int *_hblk2blks_ini,
                               long long *_hblks_ini, int _nblks_ini,
                               long long *_blks_ini, vector < int >*_p_ia_wells_ext,
                               vector < int >*_p_ja_wells_ext, _Int * _ia, _Int * _ja,
                               _Flt * _a);
      /// Prepare matrix data with blksize condensing
      void PrepareMatrixBlksize (void *_pcomm, SParams * _params, SStatData * _stats,
                                 int _nhblks_ini, int *_hblk2cpu_ini, int *_hblk2blks_ini,
                                 long long *_hblks_ini, int _nblks_ini,
                                 long long *_blks_ini, _Int * _ia, _Int * _ja, _Flt * _a);
      /// Clean partitioning data
      void CleanPartitioning ()
      {
         b_use_ini = false;
         b_use_wells = false;
         b_use_blksize = false;
         b_blk_wells = false;
         nhblks_ini = 0;
         nblks_ini = 0;
         vector < long long >hblks_ini_dummy;
         vector < int >hblk2cpu_ini_dummy;
         vector < int >hblk2blks_ini_dummy;
         vector < long long >blks_ini_dummy;
         vector < CVectorData < int > >order_ini_dummy;
         CVectorData < long long >order_wells_dummy;
         CVectorData < int >order2ind_wells_dummy;
         hblks_ini.swap (hblks_ini_dummy);
         hblk2blks_ini.swap (hblk2blks_ini_dummy);
         hblk2cpu_ini.swap (hblk2cpu_ini_dummy);
         blks_ini.swap (blks_ini_dummy);
         order_ini.swap (order_ini_dummy);
         order_wells.swap (order_wells_dummy);
         order2ind_wells.swap (order2ind_wells_dummy);
         CleanFastTransform ();
      }
      /// Clean MvmA structures
      void CleanMvmA ()
      {
         mvm.Clean ();
         vector < CBMatrix < _Int, _Flt > >hmatr_arr_dummy;
         hmatr_arr.swap (hmatr_arr_dummy);
      }
      /// Clean Slv structures
      void CleanSlv ()
      {
         vector < CBMatrix < _Int, _Flt > >sclL_dummy;
         vector < CBMatrix < _Int, _Flt > >sclU_dummy;
         vector < CBMatrix < _Int, _Flt > >sclLInv_dummy;
         vector < CBMatrix < _Int, _Flt > >sclUInv_dummy;
         vector < long long >hblks_ext_dummy;
         vector < int >nlist_ext_arr_dummy;
         vector < vector < int > >list_ext_arr_dummy;
         vector < CVectorData < int > >order_LU_dummy;
         vector < int >nblks_ext_arr_dummy;
         vector < vector < int > >blksnum_ext_arr_dummy;
         vector < vector < long long > >blks_ext_arr_dummy;
         vector < CTree > tree_arr_dummy;
         vector < int >nblks_ilu2_arr_dummy;
         vector < int >nblks1_ilu2_arr_dummy;
         vector < int >nblks2_ilu2_arr_dummy;
         vector < vector < long long > >blks_ilu2_arr_dummy;
         vector < CBMatrix < _Int, float > >matrL_float_dummy;
         vector < CBMatrix < _Int, double > >matrL_double_dummy;
         vector < CBMatrix < _Int, float > >matrU_float_dummy;
         vector < CBMatrix < _Int, double > >matrU_double_dummy;
         vector < double >coef_temp;
         CVectorData < _FltVect > xwork_temp;
         CVectorData < _FltVect > xwork_bscl_temp;
         sclL.swap (sclL_dummy);
         sclU.swap (sclU_dummy);
         sclLInv.swap (sclLInv_dummy);
         sclUInv.swap (sclUInv_dummy);
         hblks_ext.swap (hblks_ext_dummy);
         nlist_ext_arr.swap (nlist_ext_arr_dummy);
         list_ext_arr.swap (list_ext_arr_dummy);
         order_LU.swap (order_LU_dummy);
         nblks_ext_arr.swap (nblks_ext_arr_dummy);
         blksnum_ext_arr.swap (blksnum_ext_arr_dummy);
         blks_ext_arr.swap (blks_ext_arr_dummy);
         tree_arr.swap (tree_arr_dummy);
         nblks_ilu2_arr.swap (nblks_ilu2_arr_dummy);
         nblks1_ilu2_arr.swap (nblks1_ilu2_arr_dummy);
         nblks2_ilu2_arr.swap (nblks2_ilu2_arr_dummy);
         blks_ilu2_arr.swap (blks_ilu2_arr_dummy);
         matrL_float.swap (matrL_float_dummy);
         matrL_double.swap (matrL_double_dummy);
         matrU_float.swap (matrU_float_dummy);
         matrU_double.swap (matrU_double_dummy);
         slv_float.Clean ();
         slv_double.Clean ();
         ncoef_slv = 0;
         coef_slv.swap (coef_temp);
         xwork.swap (xwork_temp);
         xwork_bscl.swap (xwork_bscl_temp);
      }
      /// Clean all structures
      void Clean ()
      {
         pcomm = NULL;
         nhblks = 0;
         nblks = 0;
         vector < long long >hblks_dummy;
         vector < int >hblk2cpu_dummy;
         vector < int >hblk2blks_dummy;
         vector < int >blk2hblks_dummy;
         vector < long long >blks_dummy;
         hblks.swap (hblks_dummy);
         hblk2cpu.swap (hblk2cpu_dummy);
         hblk2blks.swap (hblk2blks_dummy);
         blk2hblks.swap (blk2hblks_dummy);
         blks.swap (blks_dummy);
         CleanPartitioning ();
         CleanMvmA ();
         CleanSlv ();
         CleanFastTransform ();
         params.SetDefaults ();
      }
// Internal functions
    private:
      /// Prepare fast transform data
      void PrepareFastTransform (CBMatrix < _Int, _Flt > *_hmatr_ini_arr);
      /// Use fast transform data
      void UseFastTransform (CBMatrix < _Int, _Flt > *_hmatr_ini_arr);
      /// Clean fast transform data
      void CleanFastTransform ()
      {
         b_fast_transform = false;
         vector < CBMatrix < _Int, double > >hmatr_ini_FT_arr_temp;
         vector < CBMatrix < _Int, double > >hmatr_fin_FT_arr_temp;
         this->hmatr_ini_FT_arr.swap (hmatr_ini_FT_arr_temp);
         this->hmatr_fin_FT_arr.swap (hmatr_fin_FT_arr_temp);
      }
      /// Peform BILU2 computations
      void ComputeBILU2 (SParams * _params, SStatData * _stats);
      /// Peform explitic block scaling
      void ExplicitBlockScaling (SParams * _params, SStatData * _stats);
      /// Perform iterations of the BiCGStab iterative scheme
      bool BiCGStab (bool _b_use_poly, SParams * _params, SStatData * _stats,
                     void *_str_mvmA, typename TMvmSlvFunc < _FltVect >::mvmA_f _mvm_f,
                     void *_str_slvLU, typename TMvmSlvFunc < _FltVect >::slvLU_f _slv_f,
                     _FltVect * _rhs, _FltVect * _sol);
      /// Perform iterations of the GMRES iterative scheme
      bool Gmres (bool _b_use_poly, SParams * _params, SStatData * _stats,
                  void *_str_mvmA, typename TMvmSlvFunc < _FltVect >::mvmA_f _mvm_f,
                  void *_str_slvLU, typename TMvmSlvFunc < _FltVect >::slvLU_f _slv_f,
                  _FltVect * _rhs, _FltVect * _sol);
      /// Reorder rhs and solution from/to initial ordering
      void ReorderVectorDataIni (char _dir, _FltVect * _rhs, _FltVect * _sol);
      /// Reorder rhs and solution from/to wells ordering
      void ReorderVectorDataWells (char _dir, _FltVect * _rhs, _FltVect * _sol,
                                   CVectorData < _FltVect > &_rhs_ord,
                                   CVectorData < _FltVect > &_sol_ord);
      /// Perform MvmA parallel computations
      void MvmA (const _FltVect * _x, _FltVect * _ax)
      {
         mvm.MvmA (_x, _ax);
      };
      /// Perform polynomial preconditioning
      void PolySlvLU (bool _b_use_poly, void *_str_mvmA,
                      typename TMvmSlvFunc < _FltVect >::mvmA_f _mvm_f, void *_str_slvLU,
                      typename TMvmSlvFunc < _FltVect >::slvLU_f _slv_f,
                      const _FltVect * _x, _FltVect * _px);
      /// Perform SlvLU parallel computations
      void SlvLU (const _FltVect * _x, _FltVect * _px);
      /// Perform mvm by part of block scaling
      void MvmBScl (char _type, const _FltVect * _x, _FltVect * _px);
    public:
      /// Multiply by A
      static void MvmA_static (void *_mvmA, const _FltVect * _x, _FltVect * _ax)
      {
         CK3D_SolverThreads < _Int, _Flt, _FltVect > *p_mvmA =
            (CK3D_SolverThreads < _Int, _Flt, _FltVect > *)_mvmA;
         p_mvmA->MvmA (_x, _ax);
      }
      /// Solve with LU for the block
      static void SlvLU_static (void *_slvLU, const _FltVect * _x, _FltVect * _px)
      {
         CK3D_SolverThreads < _Int, _Flt, _FltVect > *p_slvLU =
            (CK3D_SolverThreads < _Int, _Flt, _FltVect > *)_slvLU;
         p_slvLU->SlvLU (_x, _px);
      }

   };

}                               // namespace k3d

#endif
