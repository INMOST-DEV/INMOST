//------------------------------------------------------------------------------------------------
// File: k3d_base.hxx
//------------------------------------------------------------------------------------------------

#ifndef __K3D_BASE_HXX
#define __K3D_BASE_HXX

#define TURN_ON_EX

#if defined(TURN_ON_EX)

#ifndef USE_MPI
#define USE_MPI
#endif

#ifndef USE_SOLVER_METIS
#define USE_SOLVER_METIS
#endif

#ifndef USE_LAPACK
// #define USE_LAPACK
#endif

#else

#ifndef USE_MPI
#define USE_MPI
#endif

#ifndef USE_SOLVER_METIS
#define USE_SOLVER_METIS
#endif

#ifndef USE_LAPACK
#define USE_LAPACK
#endif

#endif

#if defined(USE_MPI)
#include <mpi.h>
#endif

#if defined(USE_SOLVER_METIS)
#define USE_METIS
#endif

#if defined(USE_METIS)
#define METIS_EXPORT
#include <metis.h>
#endif

#if defined(USE_OMP)
#define USE_THREADS
#include <omp.h>
#endif

#ifndef USE_MPI
#include <ctime>
#endif

#undef INT32_MIN
#undef INT32_MAX
#undef INT64_MIN
#undef INT64_MAX

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

namespace k3d
{

    /// Structure that contains solver parameters
    struct SParams
    {
// General params
        int msglev;               ///< Message level:
        ///<      msglev == 0 - no output,
        ///<      msglev >= 1 - output into file only,
        ///<      msglev >= 2 - output into file and screen
        ofstream *pfout;          ///< If nonzero it is the pointer to the file stream for stats output
// Decomp params
        bool b_new_part;          ///< Compute new partitioning or use previous ones
        bool b_fast_transform;    ///< Use or not fast transformation for stored partitioning (compute it on first partitioning)
        int i_decomp_type;        ///< Type of decomposition used in fct:
        ///<      i_decomp_type == 0 - user defined initial partitioning,
        ///<      i_decomp_type == 1 - threads partitioning,
        ///<      i_decomp_type == 2 - wells MPI partitioning,
        ///<      i_decomp_type == 3 - blksize based user defined initial partitioning,
        ///<      i_decomp_type == 4 - blksize based threads partitioning,
        ///<      i_decomp_type == 5 - blksize based wells MPI partitioning
        int blksize_decomp;       ///< Blksize parameter for decomp computations if applicable
        int nparts;               ///< The number of local subparts (threads) in decomposition (wells and threads interfaces only)
        int nparts_bdecomp;       ///< The number of global parts (MPI) in decomposition (wells interfaces only)
        int nneib_well;           ///< The minimal number of elems in a row to be treated as well (wells interfaces only)
        double eps_dia_well;      ///< The threshold diagonal value for row to be treated as well (wells interfaces only)
        double thresh_max_well;   ///< The threshold value for elems; rows/cols having elements with modulus larger than thresh_max_well are moved to wells; ignored if negative (default) (wells interfaces only)
        int ncycle_well;          ///< The number of wells extention cycles
        int nparts_W;             ///< The number of local subparts (threads) in decomposition of wells hblock (wells interface only)
        int degree;               ///< The maximal degree of matrix sparsity when performing decomposition (wells and threads interfaces only)
        int isize_max;            ///< The initial   (small) maximal size of set of rows treated as block in decomposition (wells and threads interfaces only)
        int isize_max2;           ///< The secondary (large) maximal size of set of rows treated as block in decomposition (wells and threads interfaces only)
// Fct params
        int nlev;                 ///< For symbolic fct the number of tree levels
        int nlev_split;           ///< For symbolic fct the first tree level for Schur complement splitting
        int prec_float;           ///< The type of preconditioner floating data to be used: prec_float == 1 (default) - use float type, prec_float == 2 - use double type
        int ncycle;               ///< The number of sparsity extension cycles when determining backward overlap for blocks, 0 < ncycle, default value ncycle = 1
        int ordlevel;             ///< The ordering filtering for Schur complement data, ordlevel == -1 - no filtering of Schur complement data
        int ordtype;              ///< The ordering type when performing extended factorization:
        ///<      ordtype ==-1 - ND ordering via METIS,
        ///<      ordtype == 0 - no ordering,
        ///<      ordtype == 1 (default) - RCM ordering for initial, main and separator parts
        int collap;               ///< The parameter that controls print sparsity condensing:
        ///<      collap == -1 (defalut) - no print,
        ///<      0 < collap - print sparsity as .ps files for each block with condensing into collap times
        int sctype_expl;          ///< The type of explicit scaling:
        ///<      sctype_expl == -1 (default) - no explicit scaling,
        ///<      sctype_expl == 0 - usual scaling via diagonal value,
        ///<      sctype_expl == 1 - balancing scaling is performed via rows/colums norms computations (block row/col QR)
        int sctype;               ///< The type of scaling:
        ///<      sctype == -1 - no scaling,
        ///<      sctype == 0 (default) - usual scaling via diagonal value/block
        ///<      sctype == 1 - balancing scaling is performed via rows/colums norms computations (block row/col QR)
        ///<      sctype == 2 - block scaling is performed (usually assuming block partitioning)
        int blksize_scale;        ///< Blksize parameter for scaling computations if block scaling is used
        int nitersc_expl;         ///< The number of explicit balancing scaling iterations, 0 < nitersc_expl, default value nitersc_expl = 3
        int nitersc;              ///< The number of balancing scaling iterations, 0 < nitersc, default value nitersc = 3
        double sclmin;            ///< The minimal scaling value, 0 < sclmin, for negative value scaling pivoting is ignored, default value sclmin = -1.0e0
        int strtype;              ///< The structural control type during fct,
        ///<      strtype == 0 - degree type structural control (default),
        ///<      strtype > 0 - fill-in type structural control
        int fcttype;              ///< The structural factorization type parameter,
        ///<      fcttype == -1 - no structural restrictions onto sparsity (default),
        ///<      0 <= fcttype - elements of order fcttype at most are remaind in fct
        int fcttype_dia;          ///< The type of fct used for dense diagonal blocks
        ///<      fcttype_dia == 0 - triangular LU fct (default),
        ///<      fcttype_dia == 1 - SVD
        double pivmin;            ///< The minimal pivot value, -1 < pivmin < 1.0, for negative value pivoting is ignored, default value pivmin = -1.0e0
        double tau1;              ///< First  threshold value, 0 <= \tau1 < 1.0, default value \tau1 = 0.01e0
        double tau2;              ///< Second threshold value, 0 <= \tau2 <= \tau1, default value \tau2 = 0.0001e0
        double tau2_sch;          ///< Second order threshold value for Schur complement, 0 <= \tau2_sch <= \tau1, default value \tau2_sch = 0.01e0
        double theta;             ///< Diagonal correction parameter, 0 <= \theta <= 1.0, default value \theta = 0.1e0
        double dia_split;         ///< Diagonal splitting parameter
// Iter params
        int it_external;          ///< it_external parameter specifies the accuracy of the external iterative scheme (-1 - no, 0 - double, 1 - quad, 2 - octal)
        int blksize_iter;         ///< Blksize parameter for iterative scheme computations if block scaling is used
        int ittype;               ///< The type of iterative scheme used:
        ///<      ittype == 0 - BiCgStab,
        ///<      ittype == 1 - GMRES,
        ///<      ittype == 2 - (1cycle GMRES, Poly)+BiCgStab(use Poly),
        ///<      ittype == 3 - (1cycle GMRES, Poly)+GMRES(use Poly)
        int niter_max;            ///< The maximal number of iterations
        int niter_cycle;          ///< The size of first GMRES cycle
        int ncoef;                ///< The number of coefs in polynomial (default value ncoef = 5)
        int niter_cycle2;         ///< The size of second GMRES cycle
        double eps_abs;           ///< The absolute stopping criterion
        double eps_rel;           ///< The relative stopping criterion
        int ichk;                 ///< The number of iterations per convergence check and per print output if any
// Dummy params
        bool b_write_file;        ///< Write data to file or not
        char name_file[256];      ///< Filename  for debug
        int iparam1;              ///< Parameter for debug
        int iparam2;              ///< Parameter for debug
        int iparam3;              ///< Parameter for debug
// Constructors and destructor
        /// Epmty constructor
        SParams ()
        {
            this->SetDefaults ();
        }
        /// Set defaults
        void SetDefaults ()
        {

            msglev = 0;
            pfout = NULL;

            b_new_part = true;
            b_fast_transform = false;
            i_decomp_type = 0;
            blksize_decomp = 1;
            nparts = 10;
            nparts_bdecomp = 5;
            nneib_well = 500;
            eps_dia_well = 1.0e-15;
            thresh_max_well = -1.0e0;
            ncycle_well = 2;
            nparts_W = 3;
            degree = 20;
            isize_max = 500;
            isize_max2 = 2000;

            nlev = 5;
            nlev_split = 3;
            prec_float = 1;
            ncycle = 1;
            ordlevel = -1;
            ordtype = 1;
            collap = -1;
            sctype_expl = -1;
//         sctype_expl = 1;
            sctype = 0;
            blksize_scale = 1;
            nitersc_expl = 5;
            nitersc = 3;
            sclmin = -1.0e0;
            strtype = 0;
            fcttype = -1;
            fcttype_dia = 0;
            pivmin = -1.0e0;
            tau1 = 1.0e-2;
            tau2 = 1.0e-4;
            tau2_sch = 1.0e-2;
            theta = 0.1e0;
            dia_split = -1.0e0;

            it_external = -1;
            blksize_iter = 1;
            ittype = 0;
            niter_max = 1000;
            niter_cycle = 30;
            ncoef = 5;
            niter_cycle2 = 5;
            eps_abs = 1.0e-10;
            eps_rel = 1.0e-10;
            ichk = 5;

            b_write_file = false;
            sprintf (name_file, "%s", "");
            iparam1 = -1;
            iparam2 = -1;
            iparam3 = -1;
        }
    };

    /// Structure that contains statistics data for output
    struct SStatData
    {
// Decomp stats
        int n_well;               ///< The total number of wells
        double dtime_part;        ///< The partitioning time
        double dtime_ord;         ///< The ordering time
// Fct stats
        double prec_extend;       ///< The scaling factor of matrix size extension due to backward overlapping onto ncycle times
        double density;           ///< The relative preconditioner size as compared to matrix size
        int nmodif_scl;           ///< The total number of pivot modifications during scaling
        double scpiv_min;         ///< The minimal scaling pivot
        double scpiv_max;         ///< The maximal scaling pivot
        int ndiasplit;            ///< The total number of diagonal blocks moved back before fct
        int nmodif;               ///< The total number of pivot modifications during fct
        double piv_min;           ///< The minimal fct pivot
        double piv_max;           ///< The maximal fct pivot
        double dtime_fct;         ///< The fct time
// Iter stats
        double rhs_norm;          ///< The rhs norm
        double res_norm;          ///< The norm of initial residual
        int niter;                ///< The number of iterations
        int nmvm;                 ///< The number of multiplications
        double res_fin;           ///< The norm of final residual
        double dtime_iter;        ///< The iterations time
// Constructors and destructor
        /// Epmty constructor
        SStatData ()
        {
            this->SetDefaults ();
        }
        /// Set defaults
        void SetDefaults ()
        {

            n_well = 0;
            dtime_part = 0.0e0;
            dtime_ord = 0.0e0;

            prec_extend = 0.0e0;
            density = 0.0e0;
            nmodif_scl = 0;
            scpiv_min = 1.0e100;
            scpiv_max = -1.0e100;
            ndiasplit = 0;
            nmodif = 0;
            piv_min = 1.0e100;
            piv_max = -1.0e100;
            dtime_fct = 0.0e0;

            rhs_norm = -1.0e0;
            res_norm = -1.0e0;
            niter = 0;
            nmvm = 0;
            res_fin = -1.0e0;
            dtime_iter = 0.0e0;

        }
    };

    /// Perform binary search
    int BinarySearch (long long _isup, int _nblks, long long *_blks, int _iblkprev);

    /// Print array
    void PrintArray (ostream & _fout, const char *_name, int _isize, const char *_iarr);
    void PrintArray (ostream & _fout, const char *_name, int _isize, const int *_iarr);
    void PrintArray (ostream & _fout, const char *_name, int _isize,
                     const long long *_iarr);
    void PrintArray (ostream & _fout, const char *_name, int _isize, const float *_farr);
    void PrintArray (ostream & _fout, const char *_name, int _isize, const double *_farr);
    void PrintArrayLow (ostream & _fout, const char *_name, int _isize,
                        const double *_farr);

    ostream & SetPw (ostream & stream);  ///< Set output manip
    void Round (double &dx, double &dy); ///< Round values

    void FGet (FILE * _file, size_t _size, char *_charr, size_t _offset);        ///< Read  char array (direct access)
    size_t LengthOfFile (const char *_filename); ///< Determine the length of ascii file
    char *ReadFile (const char *_filename, size_t & _length, bool _bprint_out);  ///< Read current file into char variable
    char ReadChar (char *_file, size_t & _k, size_t _fend);      ///< Read current char from file
    size_t ReadWord (char *_buf, char *_file, size_t & _k, size_t _fend);        ///< Read current word from file
    size_t ReadWordInBraces (char *_buf, char *_file, size_t & _k, size_t _fend);        // Read current word in braces from file
    size_t ReadSymbol (char *_file, size_t & _k, size_t _fend, char _symbol);    // Read prescribed symbol from file
    void ReplaceSymbol (size_t _length, char *_arr, char _symb_ini, char _symb_fin);     // Replace symbol in char
    void Hyst (ostream & _ffout, int _nhyst, double *_hyst, int _n_val, double *_values);        // Output hystogram for set of values

// Lapack interfaces

    extern "C"
    {

    int dlarfg_ (int *n, double *alpha, double *x, int *incx, double *tau);

    int dsyev_ (const char *jobz, const char *uplo, int *n, double *a, int *lda,
                double *w, double *work, int *lwork, int *info);

    int dgesvd_ (const char *jobu, const char *jobvt, int *m, int *n, double *a,
                 int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt,
                 double *work, int *lwork, int *info);
    }

    template < typename _T > class CVectorData   /// Support data container
    {
    protected:
        _T * pData;               ///< Data pointer
        int Length;               ///< Length of data
    private:
// Functions
        void Alloc (int _N)
        {
            if (pData != NULL) {
                delete[]pData;
            }
            Length = _N;
            pData = NULL;
            if (Length > 0) {
                pData = new _T[(size_t) Length];
            }
        }
    public:
// Constructors
        CVectorData () {
            Length = 0;
            pData = NULL;
        }
        CVectorData (int _N)
        {
            Length = 0;
            pData = NULL;
            resize (_N);
        }
// Destructor
        virtual ~ CVectorData () {
            if (pData != NULL) {
                delete[]pData;
            }
            pData = NULL;
            Length = 0;
        }
// Get/set functions
        int GetLength () const
        {
            return Length;
        };
        inline _T *Ptr ()
        {
            return pData;
        };
        const _T *Ptr () const
        {
            return pData;
        };
        void swap (CVectorData < _T > &_elem)
        {
            std::swap (pData, _elem.pData);
            std::swap (Length, _elem.Length);
        }
// Main interface function
        void resize (int _MaxCount)
        {
            Alloc (_MaxCount);
        }
        void ReplaceData (_T * _pData)
        {
            pData = _pData;
        }
        void ReplaceLength (int Length)
        {
            Length = Length;
        }
    };

    class CMPIDataExchange
    {
// Data
    public:
// Control functions
        static void Synchronize (void *_comm);    /// Barrier synchronization of the CPUs
        static int GetNproc (void *_comm);        /// Get the number of CPUs
        static int GetMyid (void *_comm); /// Get the ID number of the CPU
// Exchange functions
        /// Collective array operation
        static void ExchangeArray (void *_comm, char _Datatype, char _Operation,
                                   int _Length, void *_Array);
        /// Collective array operation (implementation)
        static void ExchangeArray_impl (void *_comm, char _Datatype, char _Operation,
                                        int _Length, void *_Array);
        /// Major data exchange function
        static void DataExchange (void *_comm, vector < int >&_CpuIDSend,
                                  vector < vector < char > >&_ObjSend,
                                  vector < int >&_CpuIDRecv,
                                  vector < vector < char > >&_ObjRecv);
        /// Major data exchange function (implementation)
        static void DataExchange_impl (void *_comm, vector < int >&_CpuIDSend,
                                       vector < vector < char > >&_ObjSend,
                                       vector < int >&_CpuIDRecv,
                                       vector < vector < char > >&_ObjRecv);
        /// Data exchange function (implementation)
        static void DataExchange_impl2 (void *_comm, int _NObjSend, int *_CpuIDSend,
                                        int *_ObjSizeSend, char **_ObjSend, int &_NObjRecv,
                                        int *&_CpuIDRecv, int *&_ObjSizeRecv,
                                        char **&_ObjRecv);
        /// Exchange data requests
        static void DataExchangeRequest (void *_comm, int _NObjSend, int *_CpuIDSend,
                                         int *_ObjSizeSend, int &_NObjRecv,
                                         int *&_CpuIDRecv, int *&_ObjSizeRecv,
                                         int *&_IndSend, int *&_IndRecv);
// Send/recv/test/wait/probe functions
        /// Synchronous send data to the other CPU
        static void Send (void *_comm, int _rank, int _msgtag, int _length, char *_buffer);
        /// Asynchronous send data to the other CPU
        static void ISend (void *_comm, int _rank, int _msgtag, int _length, char *_buffer,
                           int _indrecv, void *_recv_arr);
        /// Synchronous receive data from the other CPU
        static void Recv (void *_comm, int _rank, int _msgtag, int _length, char *_buffer);
        /// Asynchronous receive data from the other CPU
        static void IRecv (void *_comm, int _rank, int _msgtag, int _length, char *_buffer,
                           int _indrecv, void *_recv_arr);
        /// Wait for completion of exchanges
        static void WaitAll (int _count, void *_recvarr, void *_statarr);
        /// Allocate send/receive requests
        static void AllocateRecvs (int _nrecv, void *&_recvarr);
        /// Allocate statuses
        static void AllocateStats (int _nstat, void *&_statarr);
        /// Delete send/receive requests
        static void DeleteRecvs (void *_recvarr);
        /// Delete statuses
        static void DeleteStats (void *_statarr);
// Get MPI wall time
        static double GetWallTimeMPI ();
    };

    class CSortInt               /// Support sorting of int data
    {
    public:
// Data
        int ival;                 // ival contains integer value to be sorted
        int i2val;                // i2val contains secondary integer value to be sorted
// Functions
// Constructors and destructor
        CSortInt ()
        {
        };                        // Default constructor
        CSortInt (int _ival, int _i2val)
        {                         // Memory allocation constructor with zero data
            ival = _ival;
            i2val = _i2val;
        };
        ~CSortInt () {
        };                        // Destructor
// Operator functions
        CSortInt & operator= (const CSortInt & _ii2)
        {                         // Equality operator
            ival = _ii2.ival;
            i2val = _ii2.i2val;
            return *this;
        };
// Comparison function
        inline bool operator< (const CSortInt & _obj) const
        {
            return ival < _obj.ival;
        };
    };

    class CSortInt64             /// Support sorting of long int data
    {
    public:
// Data
        long long ival;           // ival contains integer value to be sorted
        long long i2val;          // i2val contains secondary integer value to be sorted
// Functions
// Constructors and destructor
        CSortInt64 ()
        {
        };                        // Default constructor
        CSortInt64 (long long _ival, long long _i2val)
        {                         // Memory allocation constructor with zero data
            ival = _ival;
            i2val = _i2val;
        };
        ~CSortInt64 () {
        };                        // Destructor
// Operator functions
        CSortInt64 & operator= (const CSortInt64 & _ii2)
        {                         // Equality operator
            ival = _ii2.ival;
            i2val = _ii2.i2val;
            return *this;
        };
// Comparison function
        inline bool operator< (const CSortInt64 & _obj) const
        {
            return ival < _obj.ival;
        };
    };

    class CSortInt2              /// Support sorting of double int data
    {
    public:
// Data
        int ixval;                // ixval contains x integer value to be sorted
        int iyval;                // iyval contains y integer value to be sorted
        int itail;                // itail contains secondary integer value to be sorted
// Functions
// Constructors and destructor
        CSortInt2 ()
        {
        };                        // Default constructor
        CSortInt2 (int _ixval, int _iyval, int _itail)
        {                         // Memory allocation constructor with zero data
            ixval = _ixval;
            iyval = _iyval;
            itail = _itail;
        };
        ~CSortInt2 () {
        };                        // Destructor
// Operator functions
        CSortInt2 & operator= (const CSortInt2 & _ii2)
        {                         // Equality operator
            ixval = _ii2.ixval;
            iyval = _ii2.iyval;
            itail = _ii2.itail;
            return *this;
        };
// Comparison function
        inline bool operator< (const CSortInt2 & _obj) const
        {
            if (ixval < _obj.ixval)
                return true;
            else if (ixval > _obj.ixval)
                return false;
            return iyval < _obj.iyval;
        };
    };

    template < typename _FltVect > class CVector /// Class that supports vector operations
    {
    public:
        static void SetByZeroes (int _n, _FltVect * _x);  ///< Set vector data by zeroes
        static void SetByZeroes_thr (int _n, _FltVect * _x);      ///< Set vector data by zeroes (threads version)
        static void SetByOnes (int _n, _FltVect * _x);    ///< Set vector data by ones
        static void SetByOnes_thr (int _n, _FltVect * _x);        ///< Set vector data by ones (threads version)
        static _FltVect ScProd (int _n, const _FltVect * _x, const _FltVect * _y);        ///< Compute scalar product
        static _FltVect ScProd_thr (int _n, const _FltVect * _x, const _FltVect * _y);    ///< Compute scalar product (threads version)
        static void CopyVector (int _n, const _FltVect * _x, _FltVect * _y);      ///< Copy array
        static void CopyVector_thr (int _n, const _FltVect * _x, _FltVect * _y);  ///< Copy array (threads version)
        static void AddVector (int _n, const _FltVect * _x1, const _FltVect * _x2, _FltVect * _x1plus2);  ///< Add vector data
        static void AddVector_thr (int _n, const _FltVect * _x1, const _FltVect * _x2, _FltVect * _x1plus2);      ///< Add vector data (threads version)
        static void AddReplaceVector (int _n, const _FltVect * _x1, _FltVect * _x1plus2); ///< Add vector data
        static void AddReplaceVector_thr (int _n, const _FltVect * _x1, _FltVect * _x1plus2);     ///< Add vector data (threads version)
        static void SubtractReplaceVector (int _n, const _FltVect * _x1, _FltVect * _x2minus1);   ///< Subtract vector data
        static void SubtractReplaceVector_thr (int _n, const _FltVect * _x1, _FltVect * _x2minus1);       ///< Subtract vector data (threads version)
        static void UpdateVector (int _n, const _FltVect * _value, const _FltVect * _arr_x, _FltVect * _arr_y);   ///< Update array (axpy): y=a*x+y
        static void UpdateVector_thr (int _n, const _FltVect * _value, const _FltVect * _arr_x, _FltVect * _arr_y);       ///< Update array (axpy): y=a*x+y (threads version)
        static void UpdateVectorMinus (int _n, const _FltVect * _value, const _FltVect * _arr_x, _FltVect * _arr_y);      ///< Update array for minus alpha (axpy)
        static void UpdateVectorMinus_thr (int _n, const _FltVect * _value, const _FltVect * _arr_x, _FltVect * _arr_y);  ///< Update array for minus alpha (axpy) (threads version)
        static void UpdateVectorReversed (int _n, const _FltVect * _value, const _FltVect * _arr_x, _FltVect * _arr_y);   ///< Update array reversed (aypx)
        static void UpdateVectorReversed_thr (int _n, const _FltVect * _value, const _FltVect * _arr_x, _FltVect * _arr_y);       ///< Update array reversed (aypx) (threads version)
        static void MultiplyVectorValue (int _n, const _FltVect * _value, _FltVect * _arr_y);     ///< Multiply array by value y=a*y
        static void MultiplyVectorValue_thr (int _n, const _FltVect * _value, _FltVect * _arr_y); ///< Multiply array by value y=a*y (threads version)
        static void InverseVector (int _n, _FltVect * _x);        ///< Compute vector inverse
        static void InverseVector_thr (int _n, _FltVect * _x);    ///< Compute vector inverse (threads version)
        static void BlockDot (int _m, int _n, const _FltVect * _amatr, const _FltVect * _bmatr, _FltVect * _cmatr);       ///< Compute block dot (C = A^t*B)
        static void BlockDot_thr (int _m, int _n, const _FltVect * _amatr, const _FltVect * _bmatr, _FltVect * _cmatr, _FltVect * _work); ///< Compute block dot (C = A^t*B) (threads version)
        static void BlockDot (int _m, int _na, int _nb, const _FltVect * _amatr, const _FltVect * _bmatr, _FltVect * _cmatr);     ///< Compute block dot (C = A^t*B)
        static void BlockDot_thr (int _m, int _na, int _nb, const _FltVect * _amatr, const _FltVect * _bmatr, _FltVect * _cmatr, _FltVect * _work);       ///< Compute block dot (C = A^t*B) (threads version)
        static void BlockDaxpy (int _m, int _n, const _FltVect * _amatr, const _FltVect * _bmatr, _FltVect * _cmatr);     ///< Compute block daxpy (C += A*B)
        static void BlockDaxpy_thr (int _m, int _n, const _FltVect * _amatr, const _FltVect * _bmatr, _FltVect * _cmatr); ///< Compute block  daxpy (C += A*B) (threads version)
        static void BlockDaxpy (int _m, int _na, int _nb, const _FltVect * _amatr, const _FltVect * _bmatr, _FltVect * _cmatr);   ///< Compute block  daxpy (C += A*B)
        static void BlockDaxpy_thr (int _m, int _na, int _nb, const _FltVect * _amatr, const _FltVect * _bmatr, _FltVect * _cmatr);       ///< Compute block  daxpy (C += A*B) (threads version)
        static void OrderVector (int _n, const int *_order, const _FltVect * _x, _FltVect * _x_ord);      ///< Order vector data
        static void OrderVector_thr (int _n, const int *_order, const _FltVect * _x, _FltVect * _x_ord);  ///< Order vector data (threads version)
        static void OrderVector (int _blksize, int _n, const int *_order, const _FltVect * _x, _FltVect * _x_ord);        ///< Order block vector data
        static void OrderVector_thr (int _blksize, int _n, const int *_order, const _FltVect * _x, _FltVect * _x_ord);    ///< Order block vector data (threads version)
        static void InvOrderVector (int _n, const int *_order, const _FltVect * _x, _FltVect * _x_ord);   ///< Inverse order vector data
        static void InvOrderVector_thr (int _n, const int *_order, const _FltVect * _x, _FltVect * _x_ord);       ///< Inverse order vector data (threads version)
        static void InvOrderVector (int _blksize, int _n, const int *_order, const _FltVect * _x, _FltVect * _x_ord);     ///< Inverse order vector data
        static void InvOrderVector_thr (int _blksize, int _n, const int *_order, const _FltVect * _x, _FltVect * _x_ord); ///< Inverse order vector data (threads version)
        /// Compute Housholder transformation
        static void Housholder (int _m, _FltVect & _alpha, _FltVect * _x, _FltVect & _tau);
        /// Compute Housholder transformation (threads version)
        static void Housholder_thr (int _m, _FltVect & _alpha, _FltVect * _x,
                                    _FltVect & _tau);
        /// Compute block Housholder transformation
        static void BlockHousholder (bool _b_use_thr, int _m, int _n, _FltVect * _q,
                                     _FltVect * _rdiag, _FltVect * _tau, _FltVect * _work);
        /// Compute QR decomposition
        static void QrdBlock (int _ncol, int _nrow, _FltVect * _qblk, int _ldq,
                              _FltVect * _tau);
        /// Compute QR decomposition (threads version)
        static void QrdBlock_thr (int _ncol, int _nrow, _FltVect * _qblk, int _ldq,
                                  _FltVect * _tau);
        /// Multiply by Q
        static void MvmQ_Housholder (int _nrhs, int _nrows, int _ncols, _FltVect * _qblk,
                                     int _ldq, _FltVect * _tau, _FltVect * _qx, int _ldqx);
        /// Multiply by Q (threads version)
        static void MvmQ_Housholder_thr (int _nrhs, int _nrows, int _ncols,
                                         _FltVect * _qblk, int _ldq, _FltVect * _tau,
                                         _FltVect * _qx, int _ldqx);
        /// Multiply by QH
        static void MvmQH_Housholder (int _nrhs, int _nrows, int _ncols, _FltVect * _qblk,
                                      int _ldq, _FltVect * _tau, _FltVect * _qx, int _ldqx);
        /// Multiply by QH (threads version)
        static void MvmQH_Housholder_thr (int _nrhs, int _nrows, int _ncols,
                                          _FltVect * _qblk, int _ldq, _FltVect * _tau,
                                          _FltVect * _qx, int _ldqx);
        /// Multiply by block Housholder transformation
        static void MvmByBlockHousholder (bool _b_use_thr, char _transp, int _nrow,
                                          int _ncol, int _nrhs, _FltVect * _qx,
                                          const _FltVect * _q, const _FltVect * _tau,
                                          _FltVect * _work);
        /// Multiply by block Housholder transformation
        static void MvmByBlockHousholder (bool _b_use_thr, char _transp, int _nrow,
                                          int _ncol, int _nrhs, _FltVect * _qxdiag,
                                          _FltVect * _qx, const _FltVect * _qdiag,
                                          const _FltVect * _q, _FltVect * _tau,
                                          _FltVect * _work);
        /// Solve triangular system
        static void SolveR (char _slvtype, int _n, _FltVect * _rmatr, int _ldr,
                            _FltVect * _rhs_sol);
        /// Compute polinomial via square upper Hessenberg matrix
        static void Polynomial (int _n, int _ncoef, _FltVect * _hmatrix, int _ldh,
                                vector < double >&_coef);
        /// Compute Cholessky decomposition for matrices stored by columns
        static void CholesskyColumns (double &_diamod, int _n, _FltVect * _amatr, int _lda,
                                      _FltVect * _uarr, int _ldu, double *_dia_arr,
                                      double &_eigmin_att, double &_eigmax_att);
        /// Compute SVD
        static void ComputeSvd (int _n, _FltVect * _amatr, _FltVect * _sv, _FltVect * _u,
                                _FltVect * _v, double *_work);
        /// Multiply matrix by matrix
        static void MMM (int _n, const _FltVect * _x_matr, const _FltVect * _y_matr,
                         _FltVect * _x_times_y_matr);
        /// Multiply matrix by matrix (C = op(A)*op(B))
        static void MMM (char _atype, char _btype, int _m, int _n, int _k,
                         const _FltVect * _amatr, int _lda, const _FltVect * _bmatr,
                         int _ldb, _FltVect * _cmatr, int _ldc);
        /// Multiply matrix by vector
        static void Mvm (char _transp, int _n, _FltVect * _x_matr, int _ldx,
                         _FltVect * _y_vect, _FltVect * _x_times_y_vect);
        /// Transpose rectangular block
        static void TransposeBlock (int _m, int _n, const _FltVect * _a_matr, int _lda,
                                    _FltVect * _at_matr, int _ldat);
    };

    template < typename _IntVect > class CVectorInt      /// Class that supports int vector operations
    {
    public:
        static void SetByZeroes (int _n, _IntVect * _ix); ///< Set vector data by zeroes
        static void SetByZeroes_thr (int _n, _IntVect * _ix);     ///< Set vector data by zeroes (threads version)
        static void SetByOnes (int _n, _IntVect * _ix);   ///< Set vector data by ones
        static void SetByOnes_thr (int _n, _IntVect * _ix);       ///< Set vector data by ones (threads version)
        static void CopyVectorInt (int _n, const _IntVect * _ix, _IntVect * _iy); ///< Copy array
        static void CopyVectorInt_thr (int _n, const _IntVect * _ix, _IntVect * _iy);     ///< Copy array (threads version)
        static void ShiftVectorInt (int _n, const _IntVect _ishift, _IntVect * _iy);      ///< Shift array
        static void ShiftVectorInt_thr (int _n, const _IntVect _ishift, _IntVect * _iy);  ///< Shift array (threads version)
        static void SetIdentity (int _n, const _IntVect _ishift, _IntVect * _iy); ///< Set shifted identity
        static void SetIdentity_thr (int _n, const _IntVect _ishift, _IntVect * _iy);     ///< Set shifted identity (threads version)
        static void InvOrder (int _n, const _IntVect * _order, _IntVect * _iorder);       ///< Compute inverse order
        static void InvOrder_thr (int _n, const _IntVect * _order, _IntVect * _iorder);   ///< Compute inverse order (threads version)
    };

/// Class that supports tree data type
    class CTree
    {
        int root_id;              ///< tree root
        int nlev;                 ///< number of levels
        int nnodes;               ///< number of nodes
        vector < int >father;   ///< father node for each node
        vector < int >nchilds;  ///< number of childs for each node
        vector < vector < int > >childs_list;    ///< lists of childs for nodes
        vector < int >nodes_lev_id;     ///< level number for each node
        vector < int >nnodes_lev;       ///< number of nodes on each level
        vector < vector < int > >nodes_lev_list; ///< list of nodes on each level
        vector < int >node2cpu; ///< cpu   number for each node
        vector < int >node2ind; ///< index number for each node
        vector < int >ind2node; ///< node number for each index
        vector < int >subtree_beg;      ///< subtree start for each node
    public:
// Functions
// Constructors and destructor
        /// Default constructor
        CTree ()
        {
            root_id = 0;
            nlev = 0;
            nnodes = 0;
            father.resize (1);
            nchilds.resize (1);
            childs_list.resize (1);
            nodes_lev_id.resize (1);
            nnodes_lev.resize (1);
            nodes_lev_list.resize (1);
            node2cpu.resize (1);
            node2ind.resize (1);
            ind2node.resize (1);
            subtree_beg.resize (1);
        };
        /// Main constructor
        CTree (int _nnodes_ini, int _nchilds);
        /// Destructor
        ~CTree () {
        };
// Operator functions
        CTree & operator= (const CTree & _tree);  ///< Copy operator
// Get/set interface functions
        /// Get root_id
        int GetRootId () const
        {
            return root_id;
        };
        /// Get nlev
        int GetNlev () const
        {
            return nlev;
        };
        /// Get nnodes
        int GetNnodes () const
        {
            return nnodes;
        };
        /// Get father
        int *GetFather ()
        {
            return father.data ();
        };
        /// Get nchilds
        int *GetNchilds ()
        {
            return &nchilds[0];
        };
        /// Get childs_list
        vector < int >*GetChildsList ()
        {
            return childs_list.data ();
        };
        /// Get nodes_lev_id
        int *GetNodesLevId ()
        {
            return nodes_lev_id.data ();
        };
        /// Get nnodes_lev
        int *GetNNodesLev ()
        {
            return nnodes_lev.data ();
        };
        /// Get nodes_lev_list
        vector < int >*GetNodesLevList ()
        {
            return nodes_lev_list.data ();
        };
        /// Get node2cpu
        int *GetNode2Cpu ()
        {
            return node2cpu.data ();
        };
        /// Get node2ind
        int *GetNode2Ind ()
        {
            return node2ind.data ();
        };
        /// Get ind2node
        int *GetInd2Node ()
        {
            return ind2node.data ();
        };
        /// Get subtree_beg
        int *GetSubtreeBeg ()
        {
            return subtree_beg.data ();
        };
// Compute functions
        /// Allocate memory
        void AllocateTree (int _nnodesmax, int _nlevmax)
        {
            father.resize (_nnodesmax + 1);
            nchilds.resize (_nnodesmax + 1);
            childs_list.resize (_nnodesmax + 1);
            nodes_lev_id.resize (_nnodesmax + 1);
            nnodes_lev.resize (_nlevmax + 1);
            nodes_lev_list.resize (_nlevmax + 1);
            node2cpu.resize (_nnodesmax + 1);
            node2ind.resize (_nnodesmax + 1);
            ind2node.resize (_nnodesmax + 1);
        }
        void FindOrderingOfNodesAccordingSubtrees (int *_ordernd);        ///< Find ordering of nodes as subtrees
        void ApplyOrderingOfTreeNodes (int *_ordernd);    ///< Apply ordering of nodes as subtrees
        /// Compute common node in up part of a tree
        int FindCommonNode (int _nlistnd, int *_listnd, int &_icycle, int *_imasknd);
        /// For each tree node find its subtree start
        void ComputeSubtreeStart (int *_subtree_start);
        /// Find subtree beg
        void InitSubtreeBeg ()
        {
            this->subtree_beg.resize (this->nnodes + 1);
            int *psubtree_beg = &this->subtree_beg[0];
            this->ComputeSubtreeStart (psubtree_beg);
        }
        /// Perform tree filtering
        void FilterTree (int *_imasknd, CTree & _tree_flt);
        void PackTree (vector < char >&_obj);     ///< Pack tree
        void UnPackTree (int _length, char *_obj);        ///< Unpack tree
        int GetPackedTreeSize (); ///< Get packed tree size
        void FillPackedTree (int _length, char *_obj);    ///< Fill packed tree
        void OutputTree (ostream & _stream);      ///< Print tree
    };

    template < typename _Int, typename _Flt > class CFct_impl    // Class that supports factorization computations
    {
    public:
        /// Compute optimal ordering for the matrix
        static void ComputeOptimalOrder (int _ordtype, const int _n,
                                         const vector < _Int > &_ia_alu,
                                         const vector < _Int > &_ja_alu,
                                         vector < int >&_order);
        /// Get main submatrix sparsity
        static void GetSubmatrixSp (int _n, const vector < _Int > &_ia,
                                    const vector < _Int > &_ja, int _n1,
                                    vector < _Int > &_ia_asub, vector < _Int > &_ja_asub);
        /// Get main submatrix sparsity
        static void GetSubmatrixSp (int _n, const vector < _Int > &_ia,
                                    const vector < _Int > &_ja,
                                    const vector < char >&_jachar, int _n1,
                                    vector < _Int > &_ia_asub, vector < _Int > &_ja_asub,
                                    vector < char >&_jachar_asub);
        /// Get last submatrix sparsity
        static void GetLastSubmatrixSp (int _n, const vector < _Int > &_ia,
                                        const vector < _Int > &_ja, int _n1,
                                        vector < _Int > &_ia_asub_last,
                                        vector < _Int > &_ja_asub_last);
        /// Get last submatrix
        static void GetLastSubmatrixBxB (bool _b_is_char, int _n, int _size_bxb,
                                         const vector < _Int > &_ia,
                                         const vector < _Int > &_ja,
                                         const vector < char >&_jachar,
                                         const vector < _Flt > &_a, int _n1,
                                         vector < _Int > &_ia_asub_last,
                                         vector < _Int > &_ja_asub_last,
                                         vector < char >&_jachar_asub_last,
                                         vector < _Flt > &_a_asub_last);
        /// Find all separators for ND ordered matrix
        static void FindAllSeparators (const int _n, const vector < _Int > &_ia_alu,
                                       const vector < _Int > &_ja_alu, int &_nblks,
                                       vector < int >&_blks);
        /// Find set of separators for given tree
        static void FindSeparatorsForTree (const int _n, const vector < _Int > &_ia_alu,
                                           const vector < _Int > &_ja_alu, CTree & _tree,
                                           int &_nblks, vector < int >&_blks);
        /// Compute optimal ordering for the matrix via splitting into 3 blocks: 1-st order Schur complement for third block, which is the separator between first block and remaining data
        static void ComputeOptimalOrderSchur (int _ordtype, const int _n, int _n1,
                                              const vector < _Int > &_ia_alu,
                                              const vector < _Int > &_ja_alu, int &_n2,
                                              vector < int >&_order);
        /// Compute ordered matrix (sparsity only)
        static void ReorderMatrixSp (int _n, const vector < int >&_order,
                                     const vector < _Int > &_ia_alu,
                                     const vector < _Int > &_ja_alu,
                                     vector < _Int > &_ia_alu_ord,
                                     vector < _Int > &_ja_alu_ord);
        /// Compute ordered matrix (sparsity only)
        static void ReorderMatrixSp (int _n, const vector < int >&_order,
                                     const vector < _Int > &_ia_alu,
                                     const vector < _Int > &_ja_alu,
                                     const vector < char >&_jachar_alu,
                                     vector < _Int > &_ia_alu_ord,
                                     vector < _Int > &_ja_alu_ord,
                                     vector < char >&_jachar_alu_ord);
        /// Compute ordered matrix
        static void ReorderMatrix (int _n, const vector < int >&_order,
                                   vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                                   vector < _Flt > &_a_alu, vector < _Int > &_ia_alu_ord,
                                   vector < _Int > &_ja_alu_ord,
                                   vector < _Flt > &_a_alu_ord);
        /// Perform ICH2 point factorization of the block with future diagonal modification
        static void Ich2Block (double _sclmin, int _fcttype, double _pivmin, double _tau1,
                               double _tau2, double _theta, int _n,
                               vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                               vector < _Flt > &_a_alu, vector < _Int > &_ia_u,
                               vector < _Int > &_ja_u, vector < _Flt > &_a_u,
                               double &_sclmin_att, double &_sclmax_att, int &_nmodif,
                               double &_eigmin_att, double &_eigmax_att);
        /// Perform ILU2 point factorization of the block with future diagonal modification
        static void Ilu2BlockTransform (int _sctype, int _nitersc, int _fcttype,
                                        double _pivmin, double _tau1, double _tau2,
                                        double _theta, int _n, vector < _Int > &_ia_alu,
                                        vector < _Int > &_ja_alu, vector < _Flt > &_a_alu,
                                        vector < _Int > &_ia_lu, vector < _Int > &_ida_lu,
                                        vector < _Int > &_ja_lu, vector < _Flt > &_a_lu,
                                        double &_eigmin_att, double &_eigmax_att);
        /// Perform ILU2 point factorization of the block with future diagonal modification
        static void Ilu2Block (int _sctype, int _nitersc, int _fcttype, double _pivmin,
                               double _tau1, double _tau2, double _theta, int _n,
                               vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                               vector < _Flt > &_a_alu, vector < _Int > &_ia_l,
                               vector < _Int > &_ja_l, vector < _Flt > &_a_l,
                               vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                               vector < _Flt > &_a_u, double &_sclmin_att,
                               double &_sclmax_att, int &_nmodif, double &_eigmin_att,
                               double &_eigmax_att);
        /// Compute symmectric scaling
        static void ComputeScalingSymm (double _sclmin, int _n, vector < _Int > &_ia_alu,
                                        vector < _Int > &_ja_alu, vector < _Flt > &_a_alu,
                                        vector < _Flt > &_sclU, double &_sclmin_att,
                                        double &_sclmax_att);
        /// Compute scaling
        static void ComputeScaling (int _sctype, int _nitersc, int _n,
                                    vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                                    vector < _Flt > &_a_alu, vector < _Flt > &_sclL,
                                    vector < _Flt > &_sclU, double &_sclmin_att,
                                    double &_sclmax_att);
        /// Perform explicit scaling
        static void MatrixScale (int _n, vector < _Flt > &_sclL, vector < _Flt > &_sclU,
                                 vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                                 vector < _Flt > &_a_alu);
        /// Symmetrize sparsity
        static void SymmetrizeSparsity (int _n, const vector < _Int > &_ia_alu,
                                        const vector < _Int > &_ja_alu,
                                        vector < _Int > &_ia_symm,
                                        vector < _Int > &_ja_symm);
        /// Symmetrize sparsity
        static void SymmetrizeSparsity (int _n, const vector < _Int > &_ia_alu,
                                        const vector < _Int > &_ja_alu,
                                        const vector < char >&_jachar_alu,
                                        vector < _Int > &_ia_symm,
                                        vector < _Int > &_ja_symm,
                                        vector < char >&_jachar_symm);
        /// Split matrix data into L and U parts for sparsity only
        static void SplitLUSp (int _n, const vector < _Int > &_ia_alu,
                               const vector < _Int > &_ja_alu, vector < _Int > &_ia_l,
                               vector < _Int > &_ja_l, vector < _Int > &_ia_u,
                               vector < _Int > &_ja_u);
        /// Split matrix data into L and U parts for sparsity only
        static void SplitLUSp (int _n, const vector < _Int > &_ia_alu,
                               const vector < _Int > &_ja_alu,
                               const vector < char >&_jachar_alu, vector < _Int > &_ia_l,
                               vector < _Int > &_ja_l, vector < char >&_jachar_l,
                               vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                               vector < char >&_jachar_u);
        /// Split matrix data into L and U parts
        static void SplitLU (bool _b_is_char, int _n, const vector < _Int > &_ia_alu,
                             const vector < _Int > &_ja_alu,
                             const vector < char >&_jachar_alu,
                             const vector < _Flt > &_a_alu, vector < _Int > &_ia_l,
                             vector < _Int > &_ja_l, vector < char >&_jachar_l,
                             vector < _Flt > &_a_l, vector < _Int > &_ia_u,
                             vector < _Int > &_ja_u, vector < char >&_jachar_u,
                             vector < _Flt > &_a_u);
        /// Transpose square matrix for sparsity only
        static void TransposeSp (int _n, const vector < _Int > &_ia_a,
                                 const vector < _Int > &_ja_a, vector < _Int > &_ia_at,
                                 vector < _Int > &_ja_at);
        /// Transpose square matrix for sparsity only
        static void TransposeSp (int _n, const vector < _Int > &_ia_a,
                                 const vector < _Int > &_ja_a,
                                 const vector < char >&_jachar_a, vector < _Int > &_ia_at,
                                 vector < _Int > &_ja_at, vector < char >&_jachar_at);
        /// Transpose square matrix
        static void Transpose (bool _b_is_char, int _n, vector < _Int > &_ia_a,
                               vector < _Int > &_ja_a, vector < char >&_jachar_a,
                               vector < _Flt > &_a_a, vector < _Int > &_ia_at,
                               vector < _Int > &_ja_at, vector < char >&_jachar_at,
                               vector < _Flt > &_a_at);
        /// Combine sparsity of L and U
        static void CombineLUSp (int _n, vector < _Int > &_ia_al, vector < _Int > &_ja_al,
                                 vector < _Int > &_ia_au, vector < _Int > &_ja_au,
                                 vector < _Int > &_ia_alu, vector < _Int > &_ja_alu);
        /// Combine sparsity of L and U
        static void CombineLUSp (int _n, vector < _Int > &_ia_al, vector < _Int > &_ja_al,
                                 vector < char >&_jachar_al, vector < _Int > &_ia_au,
                                 vector < _Int > &_ja_au, vector < char >&_jachar_au,
                                 vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                                 vector < char >&_jachar_alu);
        /// Combine L and U
        static void CombineLU (bool _b_is_char, int _n, const vector < _Int > &_ia_al,
                               const vector < _Int > &_ja_al,
                               const vector < char >&_jachar_al,
                               const vector < _Flt > &_a_al, const vector < _Int > &_ia_au,
                               const vector < _Int > &_ja_au,
                               const vector < char >&_jachar_au,
                               const vector < _Flt > &_a_au, vector < _Int > &_ia_alu,
                               vector < _Int > &_ja_alu, vector < char >&_jachar_alu,
                               vector < _Flt > &_a_alu);
        /// Combine L and U data into extended pairs
        static void CombinePairs (bool _b_is_char, int _n, vector < _Int > &_ia_al,
                                  vector < _Int > &_ja_al, vector < char >&_jachar_al,
                                  vector < _Flt > &_a_al, vector < _Int > &_ia_au,
                                  vector < _Int > &_ja_au, vector < char >&_jachar_au,
                                  vector < _Flt > &_a_au, vector < _Int > &_ia_alu,
                                  vector < _Int > &_ja_alu, vector < char >&_jachar_alu,
                                  vector < _Flt > &_a_alu);
        /// Split pairs fct data into L and U parts with post filtering
        static void SplitPairsFilter (bool _b_is_char, double _tau1, int _n,
                                      vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                                      vector < char >&_jachar_alu, vector < _Flt > &_a_alu,
                                      vector < _Int > &_ia_l, vector < _Int > &_ja_l,
                                      vector < char >&_jachar_l, vector < _Flt > &_a_l,
                                      vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                                      vector < char >&_jachar_u, vector < _Flt > &_a_u);
        /// Combine L and U data into extended pairs
        static void CombineRowsLU (int _n, vector < _Int > &_ia_al, vector < _Int > &_ja_al,
                                   vector < _Flt > &_a_al, vector < _Int > &_ia_au,
                                   vector < _Int > &_ja_au, vector < _Flt > &_a_au,
                                   vector < _Int > &_ia_alu, vector < _Int > &_ida_alu,
                                   vector < _Int > &_ja_alu, vector < _Flt > &_a_alu);
        /// Perform ILU2 point factorization of the block with future diagonal modification (no structural control)
        static void Ich2BlockIlu2 (double _pivmin, double _tau1, double _tau2,
                                   double _theta, int _n, int _n_ini,
                                   vector < _Int > &_ia_au, vector < _Int > &_ja_au,
                                   vector < _Flt > &_a_au, vector < _Int > &_ia_u,
                                   vector < _Int > &_ja_u, vector < _Flt > &_a_u,
                                   int &_nmodif, double &_eigmin_att, double &_eigmax_att);
        /// Perform ILU2 point factorization of the block with future diagonal modification (no structural control)
        static void Ilu2BlockIlu2 (double _pivmin, double _tau1, double _tau2,
                                   double _theta, int _n, int _n_ini,
                                   vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                                   vector < _Flt > &_a_alu, vector < _Int > &_ia_lu,
                                   vector < _Int > &_ja_lu, vector < _Flt > &_a_lu,
                                   int &_nmodif, double &_eigmin_att, double &_eigmax_att);
        /// Perform ILU2 point factorization of the block with future diagonal modification (with structural control)
        static void Ilu2BlockIlu2 (int _fcttype, int _fcttype_sch, double _pivmin,
                                   double _tau1, double _tau2, double _tau2_sch,
                                   double _theta, int _n, int _n_ini,
                                   vector < _Int > &_ia_alu, vector < _Int > &_ja_alu,
                                   vector < char >&_ja_char_alu, vector < _Flt > &_a_alu,
                                   vector < _Int > &_ia_lu, vector < _Int > &_ja_lu,
                                   vector < char >&_ja_char_lu, vector < _Flt > &_a_lu,
                                   int &_nmodif, double &_eigmin_att, double &_eigmax_att);
        /// Perform symbolic factorization of the block with with matrix degree structural control
        static void Ilu2BlockIlu2DegreeSp (ofstream * _pfout_debug, int _fcttype,
                                           int _fcttype_sch, int _n, int _n_ini,
                                           vector < _Int > &_ia_alu,
                                           vector < _Int > &_ja_alu,
                                           vector < char >&_ja_char_alu,
                                           vector < _Int > &_ia_lu, vector < _Int > &_ja_lu,
                                           vector < char >&_ja_char_lu);
        /// Perform ILU2 point factorization of the block with future diagonal modification (with matrix degree structural control)
        static void Ilu2BlockIlu2Degree (ofstream * _pfout_debug, int _fcttype,
                                         int _fcttype_sch, double _pivmin, double _tau1,
                                         double _tau2, double _tau2_sch, double _theta,
                                         int _n, int _n_ini, vector < _Int > &_ia_alu,
                                         vector < _Int > &_ja_alu,
                                         vector < char >&_ja_char_alu,
                                         vector < _Flt > &_a_alu, vector < _Int > &_ia_lu,
                                         vector < _Int > &_ja_lu,
                                         vector < char >&_ja_char_lu,
                                         vector < _Flt > &_a_lu, int &_nmodif,
                                         double &_eigmin_att, double &_eigmax_att);
        /// Compute inverse scaling
        static void InverseDiag (int _n, vector < _Flt > &_sclU, vector < _Flt > &_invsclU);
        /// Rescale factor back
        static void RescaleU (int _n, vector < _Flt > &_sclU, vector < _Flt > &_invsclU,
                              vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                              vector < _Flt > &_a_u);
        /// Balance diagonal
        static void BalanceDiag (int _n, vector < _Int > &_ia_a, vector < _Int > &_ja_a,
                                 vector < _Flt > &_a_a, vector < _Int > &_ia_l,
                                 vector < _Int > &_ja_l, vector < _Flt > &_a_l,
                                 vector < _Int > &_ia_u, vector < _Int > &_ja_u,
                                 vector < _Flt > &_a_u, double &_diacorr_min,
                                 double &_diacorr_max);
        /// Dense scaling
        static void DenseScaling (double _sclmin, int _n, _Flt * _a, _Flt * _sclL,
                                  _Flt * _sclU, _Flt * _sclLInv, _Flt * _sclUInv,
                                  _Flt * _work, double *_dwork, double &_sclmin_att,
                                  double &_sclmax_att, int &_nmodif);
    };

    template < typename _Int, typename _Flt > class CMatrix      // Class that supports matrix data
    {
        int b_size;               ///< The block size if any
        int n_list;               ///< The number of elements in the List and Ia arrays
        int n_list2;              ///< The number of elements in the List2 array
        int nz_ja;                ///< The number of elements in Ja array
        int nz_ja2;               ///< The number of elements in Ja2 array
        int nz_jachar;            ///< The number of elements in JaChar array
        int nz_a;                 ///< The number of elements in A array
        vector < _Int > list_matr;        ///< List part of a matrix
        vector < _Int > list2_matr;       ///< List2 part of a matrix
        vector < _Int > ia_matr;  ///< Ia part of a matrix
        vector < _Int > ja_matr;  ///< Ja part of a matrix
        vector < _Int > ja2_matr; ///< Ja2 part of a matrix
        vector < char >jachar_matr;       ///< JaChar part of a matrix
        vector < _Flt > a_matr;   ///< A  part of a matrix
// Get template names
        typedef _Int TInt_type;   ///< Reference to int type
        typedef _Flt TFlt_type;   ///< Reference to float type
// Constructors and destructor
    public:
        /// Epmty constructor
        CMatrix () {
            b_size = 1;
            n_list = 0;
            list_matr.resize (1);
            list_matr[0] = 0;
            list2_matr.resize (1);
            list2_matr[0] = 0;
            ia_matr.resize (1);
            ia_matr[0] = 0;
            ja_matr.resize (1);
            ja_matr[0] = 0;
            ja2_matr.resize (1);
            ja2_matr[0] = 0;
            jachar_matr.resize (1);
            a_matr.resize (1);
            a_matr[0] = (_Flt) 0.0e0;
            n_list2 = 0;
            nz_ja = 0;
            nz_ja2 = 0;
            nz_a = 0;
            nz_jachar = 0;
        };
        /// Copy constructor
        CMatrix (const CMatrix < _Int, _Flt > &_aa);
        /// Equality operator
        CMatrix < _Int, _Flt > &operator= (const CMatrix < _Int, _Flt > &_aa);
        /// Add replace operator
        CMatrix < _Int, _Flt > &operator+= (const CMatrix < _Int, _Flt > &_aa);
        /// Add replace operator for pairs
        CMatrix < _Int, _Flt > &operator%= (const CMatrix < _Int, _Flt > &_aa);
        /// Init by data
        CMatrix (int _n, _Int * _ia, _Int * _ja)
        {
            b_size = 1;
            n_list = _n;
            n_list2 = 0;
            int nzja = (int) _ia[_n];
            nz_ja = nzja;
            nz_ja2 = 0;
            nz_jachar = 0;
            nz_a = 0;
            list_matr.resize (1);
            list_matr[0] = 0;
            list2_matr.resize (1);
            list2_matr[0] = 0;
            ia_matr.resize (_n + 1);
            ja_matr.resize (nzja + 1);
            _Int *pia = &ia_matr[0];
            _Int *pja = &ja_matr[0];
            int i;
            for (i = 0; i <= _n; i++)
                pia[i] = _ia[i];
            for (i = 0; i < nzja; i++)
                pja[i] = _ja[i];
            ja2_matr.resize (1);
            ja2_matr[0] = 0;
            jachar_matr.resize (1);
            this->SetIdentityList ();
        };
        /// Init by data
        CMatrix (int _nlist, _Int * _list, _Int * _ia, _Int * _ja)
        {
            b_size = 1;
            n_list = _nlist;
            n_list2 = 0;
            int nzja = (int) _ia[_nlist];
            nz_ja = nzja;
            nz_ja2 = 0;
            nz_jachar = 0;
            nz_a = 0;
            list_matr.resize (_nlist + 1);
            list2_matr.resize (1);
            list2_matr[0] = 0;
            ia_matr.resize (_nlist + 1);
            ja_matr.resize (nzja + 1);
            _Int *plist = &list_matr[0];
            _Int *pia = &ia_matr[0];
            _Int *pja = &ja_matr[0];
            int i;
            for (i = 0; i < _nlist; i++)
                plist[i] = _list[i];
            for (i = 0; i <= _nlist; i++)
                pia[i] = _ia[i];
            for (i = 0; i < nzja; i++)
                pja[i] = _ja[i];
            ja2_matr.resize (1);
            ja2_matr[0] = 0;
            jachar_matr.resize (1);
            a_matr.resize (1);
            a_matr[0] = (_Flt) 0.0e0;
        };
        /// Init by data
        CMatrix (int _nlist, _Int * _list, _Int * _ia, _Int * _ja, _Flt * _a)
        {
            b_size = 1;
            n_list = _nlist;
            n_list2 = 0;
            int nzja = (int) _ia[_nlist];
            nz_ja = nzja;
            nz_ja2 = 0;
            nz_jachar = 0;
            nz_a = 0;
            list_matr.resize (_nlist + 1);
            list2_matr.resize (1);
            list2_matr[0] = 0;
            ia_matr.resize (_nlist + 1);
            ja_matr.resize (nzja + 1);
            ja2_matr.resize (1);
            ja2_matr[0] = 0;
            jachar_matr.resize (1);
            a_matr.resize (nzja + 1);
            _Int *plist = &list_matr[0];
            _Int *pia = &ia_matr[0];
            _Int *pja = &ja_matr[0];
            _Flt *pa = &a_matr[0];
            int i;
            for (i = 0; i < _nlist; i++)
                plist[i] = _list[i];
            for (i = 0; i <= _nlist; i++)
                pia[i] = _ia[i];
            for (i = 0; i < nzja; i++)
                pja[i] = _ja[i];
            for (i = 0; i < nzja; i++)
                pa[i] = _a[i];
        };
        /// Init by data
        CMatrix (int _n, _Int * _ia, _Int * _ja, _Flt * _a)
        {
            b_size = 1;
            n_list = _n;
            n_list2 = 0;
            int nzja = (int) _ia[_n];
            nz_ja = nzja;
            nz_ja2 = 0;
            nz_jachar = 0;
            nz_a = nzja;
            list_matr.resize (1);
            list_matr[0] = 0;
            list2_matr.resize (1);
            list2_matr[0] = 0;
            ia_matr.resize (_n + 1);
            ja_matr.resize (nzja + 1);
            a_matr.resize (nzja + 1);
            _Int *pia = &ia_matr[0];
            _Int *pja = &ja_matr[0];
            _Flt *pa = &a_matr[0];
            int i;
            for (i = 0; i <= _n; i++)
                pia[i] = _ia[i];
            for (i = 0; i < nzja; i++)
                pja[i] = _ja[i];
            for (i = 0; i < nzja; i++)
                pa[i] = _a[i];
            this->SetIdentityList ();
            ja2_matr.resize (1);
            ja2_matr[0] = 0;
            jachar_matr.resize (1);
        };
        /// Destructor
        ~CMatrix () {
        };
    public:
// External functions
// Get/set functions
        /// Get b_size
        int GetBSize ()
        {
            return b_size;
        };
        int GetBSize () const
        {
            return b_size;
        };
        /// Get N
        int GetN ()
        {
            return n_list;
        };
        int GetN () const
        {
            return n_list;
        };
        /// Get Nlist
        int GetNlist ()
        {
            return n_list;
        };
        int GetNlist () const
        {
            return n_list;
        };
        /// Get Nlist2
        int GetNlist2 ()
        {
            return n_list2;
        };
        int GetNlist2 () const
        {
            return n_list2;
        };
        /// Get Nzja
        int GetNzja ()
        {
            return nz_ja;
        };
        int GetNzja () const
        {
            return nz_ja;
        };
        /// Get Nzja2
        int GetNzja2 ()
        {
            return nz_ja2;
        };
        int GetNzja2 () const
        {
            return nz_ja2;
        };
        /// Get Nzjachar
        int GetNzjaChar ()
        {
            return nz_jachar;
        };
        int GetNzjaChar () const
        {
            return nz_jachar;
        };
        /// Get Nza
        int GetNza ()
        {
            return nz_a;
        };
        int GetNza () const
        {
            return nz_a;
        };
        /// Get List
        vector < _Int > *GetList () {
            return &list_matr;
        };
        const vector < _Int > *GetList () const
        {
            return &list_matr;
        };
        /// Get List array
        _Int *GetListArr ()
        {
            return list_matr.data ();
        };
        const _Int *GetListArr () const
        {
            return list_matr.data ();
        };
        /// Get List2
        vector < _Int > *GetList2 () {
            return &list2_matr;
        };
        const vector < _Int > *GetList2 () const
        {
            return &list2_matr;
        };
        /// Get List2 array
        _Int *GetList2Arr ()
        {
            return list2_matr.data ();
        };
        const _Int *GetList2Arr () const
        {
            return list2_matr.data ();
        };
        /// Get Ia
        vector < _Int > *GetIa () {
            return &ia_matr;
        };
        const vector < _Int > *GetIa () const
        {
            return &ia_matr;
        };
        /// Get Ia array
        _Int *GetIaArr ()
        {
            return ia_matr.data ();
        };
        const _Int *GetIaArr () const
        {
            return ia_matr.data ();
        };
        /// Get Ja
        vector < _Int > *GetJa () {
            return &ja_matr;
        };
        const vector < _Int > *GetJa () const
        {
            return &ja_matr;
        };
        /// Get Ja array
        _Int *GetJaArr ()
        {
            return ja_matr.data ();
        };
        const _Int *GetJaArr () const
        {
            return ja_matr.data ();
        };
        /// Get Ja2
        vector < _Int > *GetJa2 () {
            return &ja2_matr;
        };
        const vector < _Int > *GetJa2 () const
        {
            return &ja2_matr;
        };
        /// Get Ja2 array
        _Int *GetJa2Arr ()
        {
            return ja2_matr.data ();
        };
        const _Int *GetJa2Arr () const
        {
            return ja2_matr.data ();
        };
        /// Get JaChar
        vector < char >*GetJaChar ()
        {
            return &jachar_matr;
        };
        const vector < char >*GetJaChar () const
        {
            return &jachar_matr;
        };
        /// Get JaChar array
        char *GetJaCharArr ()
        {
            return jachar_matr.data ();
        };
        const char *GetJaCharArr () const
        {
            return jachar_matr.data ();
        };
        /// Get A
        vector < _Flt > *GetA () {
            return &a_matr;
        };
        const vector < _Flt > *GetA () const
        {
            return &a_matr;
        };
        /// Get A array
        _Flt *GetAArr ()
        {
            return a_matr.data ();
        };
        const _Flt *GetAArr () const
        {
            return a_matr.data ();
        };
        /// Set b_size
        void SetBSize (int _b_size)
        {
            b_size = _b_size;
        };
        /// Set n
        void SetN (int _n)
        {
            n_list = _n;
        };
        /// Set nlist
        void SetNlist (int _n)
        {
            n_list = _n;
        };
        /// Set nlist2
        void SetNlist2 (int _n)
        {
            n_list2 = _n;
        };
        /// Set nzja
        void SetNzja (int _nzja)
        {
            nz_ja = _nzja;
        };
        /// Set nzja2
        void SetNzja2 (int _nzja2)
        {
            nz_ja2 = _nzja2;
        };
        /// Set nzjachar
        void SetNzjaChar (int _nzjachar)
        {
            nz_jachar = _nzjachar;
        };
        /// Set nza
        void SetNza (int _nza)
        {
            nz_a = _nza;
        };
        /// Resize List function
        void ResizeList (int _nlist)
        {
            this->list_matr.resize (_nlist + 1);
        }
        /// Resize List2 function
        void ResizeList2 (int _nlist)
        {
            this->list2_matr.resize (_nlist + 1);
        }
        /// Resize Ia function
        void ResizeIa (int _nlist)
        {
            this->ia_matr.resize (_nlist + 1);
        }
        /// Resize Ja function
        void ResizeJa (int _nzaja)
        {
            this->ja_matr.resize (_nzaja + 1);
        }
        /// Resize Ja2 function
        void ResizeJa2 (int _nzaja)
        {
            this->ja2_matr.resize (_nzaja + 1);
        }
        /// Resize Jachar function
        void ResizeJaChar (int _nzajachar)
        {
            this->jachar_matr.resize (_nzajachar + 1);
        }
        /// Resize A function
        void ResizeA (int _nza)
        {
            this->a_matr.resize (_nza + 1);
        }
// Other
        /// Resize sp function
        void ResizeAndSetAllSp (int _nlist, int _nlist2, int _nzja, int _nzja2)
        {
            this->b_size = 1;
            this->n_list = _nlist;
            this->n_list2 = _nlist2;
            this->nz_ja = _nzja;
            this->nz_ja2 = _nzja2;
            this->nz_jachar = 0;
            this->nz_a = 0;
            this->list_matr.resize (_nlist + 1);
            this->list2_matr.resize (_nlist2 + 1);
            this->ia_matr.resize (_nlist + 1);
            this->ia_matr[0] = 0;
            this->ja_matr.resize (_nzja + 1);
            this->ja2_matr.resize (_nzja2 + 1);
            this->jachar_matr.resize (1);
            this->a_matr.resize (1);
        };
        /// Resize function
        void ResizeAndSetAll (int _nlist, int _nlist2, int _nzja, int _nzja2, int _nza)
        {
            this->b_size = 1;
            this->n_list = _nlist;
            this->n_list2 = _nlist2;
            this->nz_ja = _nzja;
            this->nz_ja2 = _nzja2;
            this->nz_a = _nza;
            this->list_matr.resize (_nlist + 1);
            this->list2_matr.resize (_nlist2 + 1);
            this->ia_matr.resize (_nlist + 1);
            this->ia_matr[0] = 0;
            this->ja_matr.resize (_nzja + 1);
            this->ja2_matr.resize (_nzja2 + 1);
            this->jachar_matr.resize (1);
            this->a_matr.resize (_nza + 1);
        };
        /// Replace function
        void ReplaceFree (CMatrix < _Int, _Flt > &_aa)
        {
            this->b_size = _aa.b_size;
            this->n_list = _aa.n_list;
            this->n_list2 = _aa.n_list2;
            this->nz_ja = _aa.nz_ja;
            this->nz_ja2 = _aa.nz_ja2;
            this->nz_jachar = _aa.nz_jachar;
            this->nz_a = _aa.nz_a;
            this->list_matr.swap (_aa.list_matr);
            this->list2_matr.swap (_aa.list2_matr);
            this->ia_matr.swap (_aa.ia_matr);
            this->ja_matr.swap (_aa.ja_matr);
            this->ja2_matr.swap (_aa.ja2_matr);
            this->jachar_matr.swap (_aa.jachar_matr);
            this->a_matr.swap (_aa.a_matr);
            _aa.Clean ();
        };
        /// Clean matrix data
        void Clean ()
        {
            vector < _Int > list_dummy (1);
            vector < _Int > list2_dummy (1);
            vector < _Int > ia_dummy (1);
            vector < _Int > ja_dummy (1);
            vector < _Int > ja2_dummy (1);
            vector < char >jachar_dummy (1);
            vector < _Flt > a_dummy (1);
            list_matr.swap (list_dummy);
            list2_matr.swap (list2_dummy);
            ia_matr.swap (ia_dummy);
            ja_matr.swap (ja_dummy);
            ja2_matr.swap (ja2_dummy);
            jachar_matr.swap (jachar_dummy);
            a_matr.swap (a_dummy);
            b_size = 1;
            n_list = 0;
            n_list2 = 0;
            nz_ja = 0;
            nz_ja2 = 0;
            nz_jachar = 0;
            nz_a = 0;
            ia_matr[0] = 0;
        };
        /// Create simple matrix test (diag and dense)
        void CreateTestMatrixSimple (int _itype, int _nx);
        /// Create sparsity of a part of test matrix (7 point)
        void Sparsity3D (int _nx, int _ny, int _nz, int _ibegz, int _iendz);
        /// Create a part of test matrix (7 point)
        void CreateTestMatrix3D (int _nx, int _ny, int _nz, int _ibegz, int _iendz,
                                 double _eps, double _dunsy);
        /// Create sparsity of a part of simple 2D test matrix (5-point Laplace)
        void Sparsity2D (int _nx, int _ny, int _ibegy, int _iendy);
        /// Create a part of simple 2D test matrix (5-point Laplace)
        void CreateTestMatrix2DLaplace (int _nx, int _ny, int _ibegy, int _iendy);
        /// Create sparsity of a part of complicated 2D test matrix (13-point Laplace^2)
        void SparsityExtended2D (int _nx, int _ny, int _ibegy, int _iendy);
        /// Create a part of complicated 2D test matrix (13-point Laplace^2)
        void CreateTestMatrix2D (int _nx, int _ny, int _ibegy, int _iendy);
        /// Set identity list
        void SetIdentityList ()
        {
            list_matr.resize (n_list + 1);
            int i;
            for (i = 0; i < n_list; i++)
                list_matr[i] = (_Int) i;
        };
        /// Get sparsity
        void GetSparsity (const CMatrix < _Int, _Flt > &_a_sp);
        /// Add sparsities
        void AddBlocksSp (const CMatrix < _Int, _Flt > &_aa, const CMatrix < _Int,
                _Flt > &_bb);
        /// Add, subtruct or replace
        void AddBlocks (char _oper, const CMatrix < _Int, _Flt > &_aa, const CMatrix < _Int,
                _Flt > &_bb);
        /// Add, subtruct or replace
        void AddBlocksBxB (char _oper, int _blksize, const CMatrix < _Int, _Flt > &_aa,
                           const CMatrix < _Int, _Flt > &_bb);
        /// Add, subtruct or replace for pairs data
        void AddBlocksPairs (char _oper, const CMatrix < _Int, _Flt > &_aa,
                             const CMatrix < _Int, _Flt > &_bb);
        /// Add, subtruct or replace for pairs data
        void AddBlocksPairsBxB (char _oper, int _blksize, const CMatrix < _Int, _Flt > &_aa,
                                const CMatrix < _Int, _Flt > &_bb);
        /// Add aa values, sparsity of aa is assumed to be included into the sparsity of current block
        void AddValues (char _oper, const CMatrix < _Int, _Flt > &_aa);
        /// Add aa values, sparsity of aa is assumed to be included into the sparsity of current block
        void AddValuesBxB (char _oper, int _blksize, const CMatrix < _Int, _Flt > &_aa);
        /// Add aa pairs of values, sparsity of aa is assumed to be included into the sparsity of current block
        void AddValuesPairs (char _oper, const CMatrix < _Int, _Flt > &_aa);
        /// Add aa pairs of values, sparsity of aa is assumed to be included into the sparsity of current block
        void AddValuesPairsBxB (char _oper, int _blksize, const CMatrix < _Int,
                _Flt > &_aa);
        /// Compute transposed sparsity for incomplete list
        void TransposedSparsityListSp (int &_icycle, int *_imask, int *_indarr, int *_iptr,
                                       int *_listloc, int *_ialoc, CMatrix < _Int,
                _Flt > &_at) const;
        /// Compute transposed matrix for incomplete list
        void TransposedSparsityList (int &_icycle, int *_imask, int *_indarr, int *_iptr,
                                     int *_listloc, int *_ialoc, CMatrix < _Int,
                _Flt > &_at) const;
        /// Compute transposed matrix for incomplete list
        void TransposedSparsityList_BxB (int _blksize, int &_icycle, int *_imask,
                                         int *_indarr, int *_iptr, int *_listloc,
                                         int *_ialoc, CMatrix < _Int, _Flt > &_at) const;
        /// Split matrix data into L and U parts for sparsity only
        void SplitLUSp (CMatrix < _Int, _Flt > &_AL_sp, CMatrix < _Int, _Flt > &_AU_sp) const
        {
            int nlist_loc = this->GetNlist ();
            int nzja_loc = this->GetNzja ();
            const vector < _Int > *pia_alu = this->GetIa ();
            const vector < _Int > *pja_alu = this->GetJa ();
            vector < _Int > *pia_l_sp = _AL_sp.GetIa ();
            vector < _Int > *pja_l_sp = _AL_sp.GetJa ();
            vector < _Int > *pia_u_sp = _AU_sp.GetIa ();
            vector < _Int > *pja_u_sp = _AU_sp.GetJa ();
            CFct_impl < _Int, _Flt >::SplitLUSp (nlist_loc, *pia_alu, *pja_alu, *pia_l_sp,
                                                 *pja_l_sp, *pia_u_sp, *pja_u_sp);
            int nzja_l = (int) (*pia_l_sp)[nlist_loc];
            _AL_sp.SetNlist (nlist_loc);
            _AL_sp.SetNzja (nzja_l);
            _AL_sp.ResizeList (nlist_loc);
            _Int *plist_l = _AL_sp.GetListArr ();
            int i;
            for (i = 0; i < nlist_loc; i++)
                plist_l[i] = (_Int) i;
            int nzja_u = (int) (*pia_u_sp)[nlist_loc];
            _AU_sp.SetNlist (nlist_loc);
            _AU_sp.SetNzja (nzja_u);
            _AU_sp.ResizeList (nlist_loc);
            _Int *plist_u = _AU_sp.GetListArr ();
            for (i = 0; i < nlist_loc; i++)
                plist_u[i] = (_Int) i;
        }
        /// Compute symmetrized matrix
        void SymmetrizeSparsitySp (CMatrix < _Int, _Flt > &_asymm) const
        {
            int nlist_loc = this->GetNlist ();
            int nzja_loc = this->GetNzja ();
            const vector < _Int > *pia_alu = this->GetIa ();
            const vector < _Int > *pja_alu = this->GetJa ();
            vector < _Int > *pia_alu_symm = _asymm.GetIa ();
            vector < _Int > *pja_alu_symm = _asymm.GetJa ();
            CFct_impl < _Int, _Flt >::SymmetrizeSparsity (nlist_loc, *pia_alu, *pja_alu,
                                                          *pia_alu_symm, *pja_alu_symm);
            int nzja_symm = (int) (*pia_alu_symm)[nlist_loc];
            _asymm.SetNlist (nlist_loc);
            _asymm.SetNzja (nzja_symm);
            _asymm.ResizeList (nlist_loc);
            _Int *plist_symm = _asymm.GetListArr ();
            int i;
            for (i = 0; i < nlist_loc; i++)
                plist_symm[i] = (_Int) i;
        }
        /// Compute symmetrized matrix
        void SymmetrizeAndAddZeroes (CMatrix < _Int, _Flt > &_asymm) const;
        /// Compute column list for transposed matrix
        void ComputeList2 (int &_icycle, int *_imask, int *_listloc);
        /// Get main submatrix
        void GetMainSubmatrixSp (int _ibeg, int _iend, CMatrix < _Int, _Flt > &_asub) const;
        /// Get last submatrix
        void GetLastSubmatrixBxB (bool _b_is_char, int _size_bxb, int _n_ini,
                                  CMatrix < _Int, _Flt > &_asub_last) const;
        /// Extend sparsity by zeroes
        void ExtendSparsity (const CMatrix < _Int, _Flt > &_a_sp, CMatrix < _Int,
                _Flt > &_a_ext) const;
        /// Extend sparsity by zeroes
        void ExtendSparsity_BxB (int _blksize, const CMatrix < _Int, _Flt > &_a_sp,
                                 CMatrix < _Int, _Flt > &_a_ext) const;
        /// Condense sparsity only
        void CondenseSparsitySp (int *_sp2blk, int &_icycle, int _nimax, int *_imask,
                                 CMatrix < _Int, _Flt > &_a_sp_cnd);
        /// Uncondense sparsity only
        void UnCondenseSparsityBlksizeSp (int _blksize, CMatrix < _Int,
                _Flt > &_a_sp_point);
        /// Condense sparsity and elems
        void CondenseSparsityBlksize (int _blksize, int *_sp2blk, int &_icycle, int _nimax,
                                      int *_imask, CMatrix < _Int, _Flt > &_a_cnd);
        /// Uncondense sparsity and elems
        void UnCondenseSparsityBlksize (int _blksize, int *_sp2blk, int &_icycle,
                                        int _nimax, int *_imask, CMatrix < _Int,
                _Flt > &_a_point);
        /// Compute optimal ordering for the matrix
        void ComputeOptimalOrder (int _ordtype, vector < int >&_order) const
        {
            int nlist = this->GetNlist ();
            const vector < _Int > *pia_alu = this->GetIa ();
            const vector < _Int > *pja_alu = this->GetJa ();
            CFct_impl < _Int, _Flt >::ComputeOptimalOrder (_ordtype, nlist, *pia_alu,
                                                           *pja_alu, _order);
        }
        /// Apply ordering to sparse matrix
        void OrderMtr (int *_order, CMatrix < _Int, _Flt > &_aord) const;
        /// Apply ordering to sparse matrix (sparsity only)
        void OrderMtrSp (int *_order, CMatrix < _Int, _Flt > &_aord) const;
        /// Apply ordering to sparse matrix and weights data
        void OrderMtrWeights (int *_order, CMatrix < _Int, _Flt > &_aordW) const;
        /// Apply ordering to matrix rows
        void OrderMtrRows (int *_order, CMatrix < _Int, _Flt > &_aord) const;
        /// Apply ordering to matrix rows (sparsity only)
        void OrderMtrRowsSp (int *_order, CMatrix < _Int, _Flt > &_aord) const;
        /// Apply ordering to matrix rows with pairs data
        void OrderMtrRowsPairs (int *_order, CMatrix < _Int, _Flt > &_aord) const;
        /// Apply ordering to matrix rows with pairs data
        void OrderMtrRowsPairs_BxB (bool _do_pairs, int _blksize, int *_order,
                                    CMatrix < _Int, _Flt > &_aord) const;
        /// Apply ordering to matrix cols
        void OrderMtrCols (int *_order, CMatrix < _Int, _Flt > &_aord) const;
        /// Apply ordering to matrix cols (sparsity only)
        void OrderMtrColsSp (int *_order, CMatrix < _Int, _Flt > &_aord) const;
        /// Apply ordering to matrix cols with pairs data
        void OrderMtrColsPairs (int *_order, CMatrix < _Int, _Flt > &_aord) const;
        /// Apply ordering to matrix cols with pairs data
        void OrderMtrColsPairs_BxB (bool _do_pairs, int _blksize, int *_order,
                                    CMatrix < _Int, _Flt > &_aord) const;
        /// Pack zero rows
        void PackZeroRows ();
        /// Compute ND order, create binary tree, find separators and condense tree
        void OrderNDSeparatorsForTree (int _nlev, int _ordtype, vector < int >&_order,
                                       CTree & _tree, int &_nblks, vector < int >&_blks);
        /// Find set of separators for given tree
        void FindSeparatorsForTree (CTree & _tree, int &_nblks, vector < int >&_blks)
        {
            int nlist = this->GetNlist ();
            const vector < _Int > *pia_alu = this->GetIa ();
            const vector < _Int > *pja_alu = this->GetJa ();
            CFct_impl < _Int, _Flt >::FindSeparatorsForTree (nlist, *pia_alu, *pja_alu,
                                                             _tree, _nblks, _blks);
        }
        /// Perform multiplication of matrices stored by rows, resulting matrix is also stored by rows
        static void MatrixByMatrixMultiplySp (int &_icycle, int *_imask, int *_imask1,
                                              int *_indarr, int *_listloc,
                                              const CMatrix < _Int, _Flt > &_a,
                                              const CMatrix < _Int, _Flt > &_b,
                                              CMatrix < _Int, _Flt > &_a_times_b);
        /// Perform multiplication of matrices stored by rows, resulting matrix is also stored by rows
        static void MatrixByMatrixMultiply (int &_icycle, int *_imask, int *_imask1,
                                            int *_indarr, int *_listloc, _Flt * _fmask,
                                            const CMatrix < _Int, _Flt > &_a,
                                            const CMatrix < _Int, _Flt > &_b,
                                            CMatrix < _Int, _Flt > &_a_times_b);
        /// Compute matrix decomposition
        static void DecompWeights (bool _split_unconnected, CMatrix < _Int,
                _Flt > &_amatr_strW, int _nparts, int *_partition);
        /// Compute 2 level matrix decomposition
        static void DecompWeights2Level (int _nblks_max, int _nparts2blk_max,
                                         CMatrix < _Int, _Flt > &_amatr_strW, int *_order,
                                         int &_nblks, vector < int >&_blk2parts,
                                         vector < int >&_parts);
        /// Split matrix
        void SplitMatrix (int _degree_max, int _isize_max, int &_nparts,
                          vector < long long >&_parts, int *_order);
        /// Compute block scaling
        void BlockScaling (int _blksize, double _sclmin, int _nlist, CMatrix < _Int,
                _Flt > &_sclL, CMatrix < _Int, _Flt > &_sclU, CMatrix < _Int,
                _Flt > &_sclLInv, CMatrix < _Int, _Flt > &_sclUInv,
                           double &_sclmin_att, double &_sclmax_att, int &_nmodif_scl);
        /// Perform in-place explicit block scaling
        void ExplicitBlockScale (int _blksize, CMatrix < _Int, _Flt > &_sclL,
                                 int *_ref_jcol, CMatrix < _Int, _Flt > &_sclU,
                                 _Flt * _work);
        /// Compute comparison matrix
        static void ComparisonMatrix_BxB (char _diatype, int _blksize, const CMatrix < _Int,
                _Flt > &_amatr_bxb, CMatrix < _Int,
                _Flt > &_amatr_pt);
        /// Print matrix data
        void PrintMatrix (ofstream & _fout);
        /// Print matrix data by rows
        void PrintMatrixRows (ofstream & _fout, int _blksize);
        /// Compute the packed size
        int GetPackedSize () const;
        /// Fill char array by the packed data
        void FillPacked (int _length, char *_obj) const;
        /// Unpack data
        void UnPack (int _length, char *_obj);
        /// Output header data
        void OutputHead (ostream & _fout);
    };

}                               // namespace k3d

#endif