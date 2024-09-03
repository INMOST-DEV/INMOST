// Collection of routines used in adaptive partitioner

#include <cfloat>
#include <cmath>
#include <cstring>

#if defined(_MSC_VER)
  #include <limits.h>
  #include <stdint.h>
#else
  #include <inttypes.h>
  #include <sys/types.h>
#endif
#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <mpi.h> 
#define _GKLIB_H_
#define METIS_EXPORT
#include <metis.h>

// Data structures for adaptive partitioning
//------------------------------------------------------------------------
// Information on vertex degrees for partitioning
typedef struct ckrinfo_t {
 int id;     // The internal degree of a vertex (sum of weights)
 int ed;     // The total external degree of a vertex
 int nnbrs;  // The number of neighboring subdomains
 int inbr;   // The index in the pool array of pairs where the nnbrs list 
             // of neighbors is stored
} ckrinfo_t;

//------------------------------------------------------------------------
// Sparse matrix in CSR format
// The diagonal entry is in the first position of each row
typedef struct matrix_t {
  int nrows, nnzs;
  int *rowptr;
  int *colind;
  double *values;
  double *transfer;
} matrix_t;

//------------------------------------------------------------------------
// Key-value pairs for real and int values
typedef struct {
  int key;
  int val;
} ikv_t;

typedef struct {
  double key;
  int    val;
} rkv_t;

//------------------------------------------------------------------------
// Class for storing refinement data
class NbrPool {
private:
  size_t size_;     // The number of entries that have been allocated
  size_t pos_;      // The position of the first free entry in the array
  size_t reallocs_; // The number of times the pool was resized

  std::vector<std::pair<int,int> > pool_; // Stores pairs (pid,ed)

public:
  NbrPool(size_t size = 0);
  ~NbrPool() {}

  size_t size()     const { return size_; }
  size_t pos()      const { return pos_; }
  size_t reallocs() const { return reallocs_; }

  void   reserve(size_t size);
  void   reset();
  int    get_next(size_t n);

  std::pair<int,int> *data() { return pool_.data(); }
};


//------------------------------------------------------------------------
struct Graph;
// Control parameters for adaptive partitioning
struct Control
{
  int   nconstraints;   // The number of balancing constraints
  int   coarse_size;    // The size of the coarsest graph
  int   nparts;         // The number of partitions

  int     local_seed;     // Local random number seed
  int     global_seed;    // Global random number seed

  std::vector<double> part_weights; // Target partition weights
  std::vector<double> inv_weights;  // Inverted constraint weights
  std::vector<double> disbalance;   // Disbalance factor (per constraint)

  bool    is_adaptive;    // flag to use adaptive partition
  bool    is_coupled;     // flag indicating if number of subdomains is matching the number of cores

  double redist_ratio;    // ratio of edge/node weights to edge/node sizes (redistribution ratio)
  double itr_ratio;       // ratio of inter-processor communication to data redistribution times
  double edge_node_ratio; // ratio of total number of edges to number of nodes

  int      myrank, ncpus; // Parallel system data
  int      nneighbors;    // The maximum number of communication neighbors this
  bool     free_gcomm;    // True if gcomm to be freed
  MPI_Comm gcomm;         // A duplicate communicator
  MPI_Comm comm;	  // The current communicator
                          // cpus needs to communicate with.
  std::vector<MPI_Request> sreqs;  // MPI send requests
  std::vector<MPI_Request> rreqs;  // MPI receive requests
  std::vector<MPI_Status>  wstats; // MPI wait statuses for send/receive requests
  MPI_Status               status;

  // These areMemory pool for use by the k-way refinement routines
  NbrPool pool;

  // Ctor/Dtor
  Control();
  ~Control();

  void set_inv_weights(const Graph&);
  void set_nnbrs(int nnbrs);
};

//------------------------------------------------------------------------
// Input graph
struct Graph {
  int gnvtxs, nvtxs, nedges, ncon;

  std::vector<int> xadj;    // Pointers to the locally stored vertices
  std::vector<int> adjncy;  // The adjacency edges of nvtxs
  std::vector<int> adjwgt;  // The weights of the adjacency edges
  std::vector<int> vtxdist; // Distribution of vertices
  std::vector<int> home;    // The initial partition of the vertex
  std::vector<int> vsize;   // Vertex size
  std::vector<int> vwgt;    // Vertex weights (int)

  std::vector<double> nvwgt;// Vertex weights (double)

  // Coarsening structures
  std::vector<int> match;
  std::vector<int> cmap;

  // Used during initial partitioning
  std::vector<int> label;

  // Communication/Setup parameters
  int nnbrs;                  // The number of neighboring processors
  int nrecv;                  // The total number of remote vertices that need to 
                              // be received. nrecv == recvptr[nnbrs]
  int nsend;                  // The total number of local vertices that need to 
                              // be sent. This corresponds to the communication 
                              // volume of each pe, in the sense that if a vertex 
                              // needs to be sent multiple times, it is accounted 
                              // in nsend. nsend == sendptr[nnbrs]
  std::vector<int> peind;     // Array of size nnbrs storing the neighboring cpus
  std::vector<int> sendptr;   // CSR format of the vertices that are sent to each
  std::vector<int> sendind;   // of the neighboring cpus
  std::vector<int> recvptr;   // CSR format of the vertices that are received from
  std::vector<int> recvind;   // each of the neighboring cpus
  std::vector<int> imap;      // The inverse map of local to global indices

  std::vector<int> pexadj;    // CSR format of this cpu each vertex is adjacent to
  std::vector<int> peadjncy;  // along with the location in the sendind of the
  std::vector<int> peadjloc;  // non-local adjacent vertices
  std::vector<int> lperm;     // lperm[0:nlocal] points to interior vertices,
                              // the rest are interface
  int nlocal;                 // Number of interior vertices

  // Communication parameters for projecting the partition. 
  // These are computed during CreateCoarseGraph and used during projection 
  std::vector<int> rlens;     // Arrays of size nnbrs of how many vertices
  std::vector<int> slens;     // this cpu is sending/receiving
  std::vector<ikv_t> rcand;


  // Partition parameters
  std::vector<int>    where;
  std::vector<double> lnpwgts;
  std::vector<double> gnpwgts;
  std::vector<ckrinfo_t> ckrinfo;

  int lmincut, mincut;

  int level;
  bool global_match;

  struct Graph *coarser, *finer;

  // Ctor/Dtor
  Graph();
  ~Graph();

  void set_weights(const Control&);
  void remap(const Control&);
  void renumber_conns(int *adjncy);
};

//-------------------------------------------------------------------------

static const int DEFAULT_SEED = 12345;

static const int OPTIONAL_SEED_IDX = 2;
static const int COUPLING_IDX      = 3;

static const int NGD_STEPS    = 10; // Maximal number of diffusion steps
static const int NGR_STEPS    =  4; // Number of greedy refinement passes

static const double COARSEN_RATIO    = 0.75; // Node reduction between succesive coarsening levels
static const double DISBALANCE_RATIO = 1.05; // Desired disbalance ratio

static const double AVG_THRESHOLD = 0.03;
static const double SMALL_FLOAT   = 1.e-6;

static const int  UNMATCHED     = -1;
static const int  MAYBE_MATCHED = -2;
static const int  TOO_HEAVY     = -3;

//-------------------------------------------------------------------------
// Blas routines

template <typename T>
T sum(size_t n, const T *arr, size_t incx)
{
  T sum = (T)0;
  for(size_t i=0; i < n; arr+=incx, ++i) sum += (*arr);
  return sum;
}

template int    sum<int>   (size_t n, const int    *arr, size_t incx);
template double sum<double>(size_t n, const double *arr, size_t incx);
int    isum(size_t n, const int    *arr, size_t incx) { return sum<int>   (n, arr, incx); }
double rsum(size_t n,       double *arr, size_t incx) { return sum<double>(n, arr, incx); }

template <typename T>
T *set(size_t n, T val, T *arr)
{
  for(size_t i=0; i < n; ++i) arr[i] = val;
  return arr;
}

template int    *set<int>   (size_t n, int    val, int    *arr);
template double *set<double>(size_t n, double val, double *arr);
int  *iset(size_t n, int  val, int  *arr) { return set<int>   (n, val, arr); }
double *rset(size_t n, double val, double *arr) { return set<double>(n, val, arr); }

template <typename T>
T *copy(size_t n, const T *src, T *dest)
{ return (T*)std::memmove(dest, src, sizeof(T)*n); }

template int    *copy<int>   (size_t n, const int    *src, int    *dest);
template double *copy<double>(size_t n, const double *src, double *dest);
template ikv_t  *copy<ikv_t> (size_t n, const ikv_t  *src, ikv_t  *dest);
int    *icopy  (size_t n, int         *src, int    *dest) { return copy<int>   (n, src, dest); }
double *rcopy  (size_t n, double      *src, double *dest) { return copy<double>(n, src, dest); }
ikv_t  *ikvcopy(size_t n, const ikv_t *src, ikv_t  *dest) { return copy<ikv_t>(n, src, dest); }

int *iincset(size_t n, int base, int *x)
{
  for(size_t i=0; i < n; ++i) x[i] = base + (int)i;
  return x;
}

double rmax(size_t n, double *x, size_t incx)
{
  size_t i;
  double max;
  if( n <= 0 ) return 0.;
  for(max=(*x), x+=incx, i=1; i < n; ++i, x+=incx)
    max = ((*x) > max ? (*x) : max);
  return max;
}

double rmin(size_t n, double *x, size_t incx)
{
  size_t i;
  double min;
  if( n <= 0 ) return 0.;
  for(min=(*x), x+=incx, i=1; i < n; ++i, x+=incx)
    min = ((*x) < min ? (*x) : min);
  return min;
}

double rnorm2(size_t n, double *x, size_t incx)
{
  double sum = 0;
  for(size_t i=0; i < n; ++i, x+=incx)  sum += (*x) * (*x);
  return (sum > 0. ? std::sqrt(sum) : 0.); 
}

double rdot(size_t n, double *x, size_t incx, double *y, size_t incy)
{
  double sum = 0;
  for(size_t i=0; i < n; ++i, x+=incx, y+=incy)  sum += (*x) * (*y);
  return sum;
}

double *raxpy(size_t n, double alpha, double *x, size_t incx, double *y, size_t incy)
{
  double *y_in = y;
  for(size_t i=0; i < n; ++i, x+=incx, y+=incy)  *y += alpha*(*x);
  return y_in;
}

//------------------------------------------------------------------------
template <typename T>
int max_index(size_t n, const T *x, size_t incx=1)
{
  size_t max_ind = 0;
  n *= incx;
  for(size_t i=incx; i < n; i+=incx)  max_ind = (x[i] > x[max_ind] ? i : max_ind);
  return (int)(max_ind / incx);
}

template int max_index<int>   (size_t n, const int    *arr, size_t incx);
template int max_index<double>(size_t n, const double *arr, size_t incx);

// Return the index of the almost maximum element in a vector
int max_index2(size_t n, double *x)
{
  size_t max1, max2;

  if( x[0] > x[1] ) {
    max1 = 0;
    max2 = 1;
  } else {
    max1 = 1;
    max2 = 0;
  }

  for(size_t i=2; i < n; ++i) {
    if( x[i] > x[max1] ) {
      max2 = max1;
      max1 = i;
    } else if( x[i] > x[max2] )
      max2 = i;
  }

  return (int)max2;
}

size_t iargmax(size_t n, int *x, size_t incx) // TODO: remove after testing
{
  size_t i, j, max=0;
  for(i=1, j=incx; i < n; ++i, j+=incx)  max = (x[j] > x[max] ? j : max);
  return (size_t)(max / incx);
}

size_t rargmax_strd(size_t n, double *x, size_t incx) // TODO: remove after testing
{
  size_t max = 0;
  n *= incx;
  for(size_t i=incx; i < n; i+=incx)  max = (x[i] > x[max] ? i : max);
  return max / incx;
}

double average(size_t n, double *x)
{
  double sum = 0.0;
  for(size_t i=0; i < n; ++i)  sum += x[i];
  return sum / n;
}


//-------------------------------------------------------------------------
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wregister"
#endif

// GKlib routines
#include <gk_macros.h>
#include <gk_mksort.h>
//extern "C" { void *gk_malloc(size_t, char*); }

// Sorting routines
void isorti(size_t n, int *base)
{
#define i_lt(a, b) ((*a) < (*b))
  GK_MKQSORT(int, base, n, i_lt);
#undef i_lt
}

void ikvsorti(size_t n, ikv_t *base)
{
#define ikey_lt(a, b) ((a)->key < (b)->key)
  GK_MKQSORT(ikv_t, base, n, ikey_lt);
#undef ikey_lt
}

// Sorts based both on key and val
void ikvsortii(size_t n, ikv_t *base)
{
#define ikeyval_lt(a, b) ((a)->key < (b)->key || ((a)->key == (b)->key && (a)->val < (b)->val))
  GK_MKQSORT(ikv_t, base, n, ikeyval_lt);
#undef ikeyval_lt
}

void ikvsortd(size_t n, ikv_t *base)
{
#define ikey_gt(a, b) ((a)->key > (b)->key)
  GK_MKQSORT(ikv_t, base, n, ikey_gt);
#undef ikey_gt
}

void rkvsorti(size_t n, rkv_t *base)
{
#define rkey_lt(a, b) ((a)->key < (b)->key)
  GK_MKQSORT(rkv_t, base, n, rkey_lt);
#undef rkey_lt
}

void rkvsortd(size_t n, rkv_t *base)
{
#define rkey_gt(a, b) ((a)->key > (b)->key)
  GK_MKQSORT(rkv_t, base, n, rkey_gt);
#undef rkey_gt
}

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
//------------------------------------------------------------------------
#if defined(_MSC_VER)
typedef ptrdiff_t ssize_t;
#endif

typedef struct {
  size_t   nnodes;
  size_t   maxnodes;
  rkv_t   *heap;
  ssize_t *locator;
} rpq_t;

// Create and initialize a priority queue
rpq_t *rpqCreate(size_t maxnodes)
{
  rpq_t *queue = new rpq_t;

  queue->nnodes = 0;
  queue->maxnodes = maxnodes;

  queue->heap    = new rkv_t[maxnodes];
  queue->locator = new ssize_t[maxnodes];
  for(size_t i=0; i < maxnodes; ++i) queue->locator[i] = -1;

  return queue;
}

// Reset the priority queue
void rpqReset(rpq_t *queue)
{
  for(ssize_t i=queue->nnodes-1; i >= 0; --i)
    queue->locator[queue->heap[i].val] = -1;
  queue->nnodes = 0;
}

// Free the internal datastructures of the priority queue and the queue itself
void rpqDestroy(rpq_t *queue)
{
  if( !queue ) return;

  delete queue->heap;
  delete queue->locator;
  queue->maxnodes = 0;

  delete queue;
}

// Return the length of the queue
size_t rpqLength(rpq_t *queue)
{
  return queue->nnodes;
}

// Check the consistency of the heap
int rpqCheckHeap(rpq_t *queue)
{
  if( queue->nnodes == 0 )  return 1;

  assert(queue->locator[queue->heap[0].val] == 0);
  for(size_t i=1; i < queue->nnodes; ++i) {
    assert(queue->locator[queue->heap[i].val] == i);
    assert(!(queue->heap[i].key > queue->heap[(i-1)/2].key));
  }
  for(size_t i=1; i < queue->nnodes; ++i)
    ASSERT(!(queue->heap[i].key > queue->heap[0].key));

  size_t j = 0;
  for(size_t i=0; i < queue->maxnodes; i++ ) {
    if( queue->locator[i] != -1 )
      ++j;
  }
  assert(j == queue->nnodes);
  return 1;
}

// Add an item in the queue
int rpqInsert(rpq_t *queue, int node, double key)
{
  ssize_t *locator = queue->locator;
  rkv_t   *heap    = queue->heap;

  assert(rpqCheckHeap(queue));
  assert(locator[node] == -1);

  ssize_t i = queue->nnodes++;
  while( i > 0 ) {
    ssize_t j = (i-1)>>1;
    if( key > heap[j].key ) {
      heap[i] = heap[j];
      locator[heap[i].val] = i;
      i = j;
    } else
      break;
  }
  assert( i >= 0 );
  heap[i].key   = key;
  heap[i].val   = node;
  locator[node] = i;

  assert(rpqCheckHeap(queue));
  return 0;
}

// Delete an item from the priority queue
int rpqDelete(rpq_t *queue, int node)
{
  ssize_t *locator = queue->locator;
  rkv_t   *heap    = queue->heap;

  assert(locator[node] != -1);
  assert(heap[locator[node]].val == node);
  assert(rpqCheckHeap(queue));

  ssize_t j, i = locator[node];
  locator[node] = -1;

  ssize_t nnodes;
  if( --queue->nnodes > 0  &&  heap[queue->nnodes].val != node ) {
    node = heap[queue->nnodes].val;
    double newkey = heap[queue->nnodes].key;
    double oldkey = heap[i].key;

    if( newkey > oldkey ) { // Filter-up 
      while( i > 0 ) {
        j = (i-1)>>1;
        if( newkey > heap[j].key ) {
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        } else
          break;
      }
    } else { // Filter down
      nnodes = queue->nnodes;
      while( (j=(i<<1)+1) < nnodes ) {
        if( heap[j].key > newkey ) {
          if( j+1 < nnodes  &&  (heap[j+1].key > heap[j].key) )
            ++j;
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        } else if( j+1 < nnodes  &&  (heap[j+1].key > newkey) ) {
          ++j;
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        } else
          break;
      }
    }

    heap[i].key   = newkey;
    heap[i].val   = node;
    locator[node] = i;
  }

  assert(rpqCheckHeap(queue));

  return 0;
}

// Update the key values associated for a particular item
void rpqUpdate(rpq_t *queue, int node, double newkey)
{
  ssize_t *locator = queue->locator;
  rkv_t   *heap    = queue->heap;

  double oldkey = heap[locator[node]].key;
  if( !(newkey > oldkey)  &&  !(oldkey > newkey) ) return;

  assert(locator[node] != -1);
  assert(heap[locator[node]].val == node);
  assert(rpqCheckHeap(queue));

  ssize_t j, i = locator[node];
  ssize_t nnodes;
  if( newkey > oldkey ) { // Filter-up
    while( i > 0 ) {
      j = (i-1)>>1;
      if( newkey > heap[j].key ) {
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      } else
        break;
    }
  } else { // Filter down
    nnodes = queue->nnodes;
    while( (j=(i<<1)+1) < nnodes ) {
      if( heap[j].key > newkey ) {
        if( j+1 < nnodes  &&  (heap[j+1].key > heap[j].key) )
          ++j;
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      } else if( j+1 < nnodes  &&  (heap[j+1].key > newkey) ) {
        ++j;
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      } else
        break;
    }
  }

  heap[i].key   = newkey;
  heap[i].val   = node;
  locator[node] = i;

  assert(rpqCheckHeap(queue));
  return;
}

// Return the item at the top of the queue and remove it from the queue
int rpqGetTop(rpq_t *queue)
{
  assert(rpqCheckHeap(queue));

  if( queue->nnodes == 0 ) return -1;

  queue->nnodes--;

  rkv_t   *heap    = queue->heap;
  ssize_t *locator = queue->locator;

  int vtx = heap[0].val;
  locator[vtx] = -1;

  ssize_t i, j;
  ssize_t nnodes = queue->nnodes;
  if( (i = nnodes) > 0 ) {
    double key  = heap[i].key;
    int  node = heap[i].val;
    i = 0;
    while( (j=2*i+1) < nnodes) {
      if( (heap[j].key > key) ) {
        if( j+1 < nnodes  &&  (heap[j+1].key > heap[j].key) )
          j = j+1;
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      } else if( j+1 < nnodes  &&  (heap[j+1].key > key) ) {
        j = j+1;
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      } else
        break;
    }

    heap[i].key   = key;
    heap[i].val   = node;
    locator[node] = i;
  }

  assert(rpqCheckHeap(queue));
  return vtx;
}

// Return the item at the top of the queue; it will not be deleted from the queue
int rpqSeeTopVal(rpq_t *queue)
{
  return (queue->nnodes == 0 ? -1 : queue->heap[0].val);
}

// Return the key of the top item; it will not be deleted from the queue
double rpqSeeTopKey(rpq_t *queue)
{
  return (queue->nnodes == 0 ? REAL_MAX : queue->heap[0].key);
}

// Return the key of a specific item
double rpqSeeKey(rpq_t *queue, int node)
{
  return queue->heap[queue->locator[node]].key;
}

//-------------------------------------------------------------------------

#ifndef ASSERT
// Debugging macros
#ifndef NDEBUG
# define ASSERT(expr)                                       \
  if( !(expr) ) {                                           \
    std::cout << "**ASSERTION failed on line " << __LINE__  \
              << " of file " << __FILE__ << ": " << #expr   \
	      << std::endl;                                 \
    throw std::runtime_error("error");                      \
  }
#else
# define ASSERT(expr) ;
#endif
#endif

//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Check the validity of the input arguments
int check_inputs(const int *vtxdist, const int *xadj,  const int *adjncy, 
                 const int *vwgt,    const int *vsize, const int *adjwgt,
                 int weight_flag, int ncons, int nparts,
                 const double *part_weights, const double *disbalance,
                 double itr_ratio, const int *options, const int *part,
                 MPI_Comm *comm)
{
  // Check the validity of supplied pointers
  if( !comm ) {
    std::runtime_error("AdaptivePartitioner: invalid communicator");
  }
  int myrank;
  MPI_Comm_rank(*comm, &myrank);

  ASSERT(vtxdist);
  ASSERT(xadj);
  ASSERT(adjncy);
  //ASSERT(vsize); // unnecessary parameter
  ASSERT(part_weights);
  ASSERT(disbalance);
  ASSERT(options);
  ASSERT(part);

  // Check the validity of supplied data
  if( weight_flag == 1 || weight_flag == 3 ) {
    ASSERT(adjwgt);
  }
  if( weight_flag == 2 || weight_flag == 3 ) {
    ASSERT(vwgt);
    for(int i=0; i < ncons; ++i) {
      int isum = sum(vtxdist[myrank+1]-vtxdist[myrank], vwgt+i, ncons);
      int imin;
      MPI_Allreduce(&isum, &imin, 1, MPI_INT, MPI_SUM, *comm);

      if( imin == 0 ) {
        if( myrank == 0 )
          std::cout << "AdaptivePartitioner: zero weight for constraint " << i << std::endl;
        throw std::runtime_error("error");
      }
    }
  }

  if( vtxdist[myrank+1]-vtxdist[myrank] < 1 ) {
    std::cout << "AdaptivePartitioner: bad vertex distribution: "
	      << "CPU(" << myrank << ") has no vertices assigned"
	      << std::endl;
    return 0;
  }

  for(int i=0; i < ncons; ++i) {
    double rsum = sum(nparts, part_weights+i, ncons);
    if( std::abs(rsum - 1.0) > 0.001 ) {
      std::cout << "AdaptivePartitioner: the sum of part_weights for constraint "
		<< i << " is not equal to 1"
		<< std::endl;
      return 0;
    }
  }

  for(int i=0; i < ncons; ++i) {
    for(int j=0; j < nparts; ++j) {
      if( part_weights[j*ncons+i] < 0.0 || part_weights[j*ncons+i] > 1.001) {
        std::cout << "AdaptivePartitioner: the part_weights for constraint "
          << i << " and partition " << j << " is out of bounds [0,1]"
          << std::endl;
        return 0;
      }
    }
  }

  for(int i=0; i < ncons; ++i) {
    if( disbalance[i] <= 1.0 ) {
      std::cout << "AdaptivePartitioner: the disbalance for constraint "
		<< i << " must be > 1"
		<< std::endl;
      return 0;
    }
  }

  if( itr_ratio < 0.001 || itr_ratio > 1000.0 ) {
    std::cout << "AdaptivePartitioner: the itr_ratio value should be "
        << "between [0.001, 1000]"
        << std::endl;
    return 0;
  }

  return 1;
}


//-------------------------------------------------------------------------
// Class NbrPool
NbrPool::NbrPool(size_t size)
  : size_(0), pos_(0), reallocs_(0), pool_()
{
  reserve(size);
}

void NbrPool::reserve(size_t size)
{
  size_ = size;
  pool_.resize(size_);
}

void NbrPool::reset()
{
  pos_ = 0;
}

// Get next free index from the pool
int NbrPool::get_next(size_t nnbrs)
{
  pos_ += nnbrs;

  if( pos_ > size_ ) {
    size_ += std::max(20*nnbrs, size_/2);
    pool_.resize(size_);
    ++reallocs_;
  }

  return (int)(pos_ - nnbrs);
}


//------------------------------------------------------------------------
Control::Control()
  : nconstraints(0),
    coarse_size(0),
    nparts(0),
    local_seed(0),
    global_seed(0),
    part_weights(),
    inv_weights(),
    disbalance(),
    is_adaptive(false),
    is_coupled(false),
    redist_ratio(0.),
    itr_ratio(0.),
    edge_node_ratio(0.),
    myrank(0),
    ncpus(0),
    nneighbors(0),
    free_gcomm(false)
{}

Control::~Control()
{
  if( free_gcomm )
    MPI_Comm_free(&gcomm);
}

// Compute the inv_weights of a graph and stores them in ctrl
void Control::set_inv_weights(const Graph& graph)
{
  inv_weights.resize(graph.ncon);

  for(int j=0; j < graph.ncon; ++j) {
    int wgt_sum_loc = sum(graph.nvtxs, &graph.vwgt[j], graph.ncon);
    int wgt_sum_glb;
    MPI_Allreduce(&wgt_sum_loc, &wgt_sum_glb, 1, MPI_INT, MPI_SUM, comm);

    inv_weights[j] = 1. / wgt_sum_glb;
  }
}

// Update the sresq/rreqs/wstats arrays to new number of cpus
void Control::set_nnbrs(int nnbrs)
{
  if( nneighbors >= nnbrs ) return;

  nneighbors = nnbrs;
  wstats.resize(nnbrs);
  rreqs .resize(nnbrs);
  sreqs .resize(nnbrs);
}

//------------------------------------------------------------------------
// Set the Control structure
Control *setup_controls(bool adaptive, int nconstraints, int nparts, const int *options, 
			double *part_weights_in, double *disbalance_in, MPI_Comm comm)
{
  Control *ctrl = new Control();

  // Communicator-related info
  MPI_Comm_dup(comm, &(ctrl->gcomm));
  ctrl->comm = ctrl->gcomm;
  ctrl->free_gcomm = true;
  MPI_Comm_rank(ctrl->gcomm, &ctrl->myrank);
  MPI_Comm_size(ctrl->gcomm, &ctrl->ncpus);

  // Options[]-related info
  bool defaults = (options == NULL ? true : options[0] == 0);
  if( adaptive ) {
    ctrl->is_adaptive = true;
    ctrl->is_coupled  = ((ctrl->ncpus == nparts)
			 ? (defaults ? true : options[COUPLING_IDX] == 1)
			 : false);
  } else {
    ctrl->is_adaptive = false;
    ctrl->is_coupled  = false;
  }

  ctrl->local_seed  = (defaults ? DEFAULT_SEED : options[OPTIONAL_SEED_IDX]);
  MPI_Allreduce(&ctrl->local_seed, &ctrl->global_seed, 1, MPI_INT, MPI_MAX, ctrl->comm);

  ctrl->local_seed  = (ctrl->local_seed == 0 ? ctrl->myrank : ctrl->local_seed*ctrl->myrank);

  // Common info
  ctrl->nconstraints = nconstraints;    
  ctrl->nparts       = nparts;    
  ctrl->redist_ratio = 1.;

  // Setup target partition weights
  ctrl->part_weights.resize(nparts*nconstraints);
  if( part_weights_in ) {
    ctrl->part_weights.assign(part_weights_in, part_weights_in + nparts*nconstraints);
  } else {
    for(int i=0; i < nparts; ++i) {
      for(int j=0; j < nconstraints; ++j)
        ctrl->part_weights[i*nconstraints+j] = 1. / nparts;
    }
  }

  // Setup disbalance
  ctrl->disbalance.assign(nconstraints, DISBALANCE_RATIO);
  if( disbalance_in )
    copy(nconstraints, disbalance_in, ctrl->disbalance.data());

  // Initialize the random number generator
  srand(ctrl->local_seed);

  return ctrl;
}


//------------------------------------------------------------------------
// De-allocate memory allocated for the control structures
void free_controls(Control **r_ctrl)
{
  Control *ctrl = *r_ctrl;

  delete ctrl;

  *r_ctrl = NULL;
}


//------------------------------------------------------------------------
// Perform the following operations:
//   - determine the cpus containing adjacent vertices and setup
//     the infrastructure for efficient communication
//   - localize the numbering of the adjancency lists
void setup_comm(Control *ctrl, Graph *graph)
{
  if( !graph->lperm.empty() )
    return; // The communication structure has already been setup

  int nvtxs  = graph->nvtxs;
  graph->lperm.resize(nvtxs);
  iincset(nvtxs, 0, graph->lperm.data());

  int first_vtx = graph->vtxdist[ctrl->myrank];
  int last_vtx  = graph->vtxdist[ctrl->myrank+1];

  //-------------------------------------------------------------
  // Determine data to be sent/received
  // Step 1: determine nadj and interior/interface vertices
  int nlocal = 0, nadj = 0;
  for(int i=0; i < nvtxs; ++i) {
    bool islocal = true;
    for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
      int k = graph->adjncy[j];
      if( k < first_vtx || k >= last_vtx ) { // remote vertex
        ++nadj;
        islocal = false;
      }
    }
    if( islocal ) {
      graph->lperm[i]        = graph->lperm[nlocal];
      graph->lperm[nlocal++] = i;
    }
  }
  graph->nlocal = nlocal;

  std::vector<ikv_t> adjpairs(nadj+1);

  // Step 2: rewrite locale entries and populate remote edges
  nadj = 0;
  for(int i=0; i < nvtxs; ++i) {
    for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
      int k = graph->adjncy[j];
      if( k >= first_vtx && k < last_vtx ) { // local vertex
        graph->adjncy[j] = k - first_vtx; 
      } else { // Remote vertex
        adjpairs[nadj].key   = k;
        adjpairs[nadj++].val = j;
      }
    }
  }

  // Sort pairs
  ikvsorti(nadj, adjpairs.data());
  adjpairs[nadj].key = graph->gnvtxs+1;

  // Determinenumber of vertices to be received
  int nrecv = 0;
  for(int i=0; i < nadj; ++i) {
    if( adjpairs[i].key != adjpairs[i+1].key )
      ++nrecv;
  }
  graph->nrecv = nrecv;

  // Allocate space for vertices to be received
  graph->recvind.resize(nrecv);
  int *recvind = graph->recvind.data();

  // Store distinct vertices in recvind array and re-write adjncy
  nrecv = 0;
  for(int i=0; i < nadj; ++i) {
    graph->adjncy[adjpairs[i].val] = nvtxs + nrecv;
    if( adjpairs[i].key != adjpairs[i+1].key )
      recvind[nrecv++] = adjpairs[i].key;
  }
  ASSERT(nrecv == graph->nrecv);

  // Determine the number of neighboring cpus
  int nnbrs = 0;
  for(int j=0, i=0, icpu=0; icpu < ctrl->ncpus; ++icpu) {
    for(j=i; j < nrecv; ++j) {
      if( recvind[j] >= graph->vtxdist[icpu+1] )
        break;
    }
    if( j > i ) {
      ++nnbrs;
      i = j;
    }
  }
  graph->nnbrs = nnbrs;

  // Update communication part of ctrl arrays
  ctrl->set_nnbrs(nnbrs);

  // Allocate space for peind/recvptr part of the recvinfo
  graph->peind.resize(nnbrs);
  graph->recvptr.resize(nnbrs+1);
  int *recvptr = graph->recvptr.data();

  // Populate receiving part of arrays
  nnbrs = 0;
  recvptr[0] = 0;
  for(int j=0, i=0, icpu=0; icpu < ctrl->ncpus; ++icpu) {
    for(j=i; j < nrecv; ++j) {
      if( recvind[j] >= graph->vtxdist[icpu+1] )
        break;
    }
    if( j > i ) {
      graph->peind[nnbrs++] = icpu;
      recvptr[nnbrs] = j;
      i = j;
    }
  }
  ASSERT(nnbrs == graph->nnbrs);

  //-------------------------------------------------------------
  // Determine data to send 
  std::vector<ikv_t> recv_reqs(ctrl->ncpus);
  std::vector<ikv_t> send_reqs(ctrl->ncpus);
  memset(recv_reqs.data(), 0, sizeof(ikv_t)*ctrl->ncpus);
  for(int i=0; i < nnbrs; ++i) {
    recv_reqs[graph->peind[i]].key = recvptr[i+1] - recvptr[i];
    recv_reqs[graph->peind[i]].val = nvtxs+recvptr[i];
  }
  MPI_Alltoall(recv_reqs.data(), 2, MPI_INT, send_reqs.data(), 2, MPI_INT, ctrl->comm);

  std::vector<int> startsind(nnbrs);
  graph->sendptr.resize(nnbrs+1);
  int *sendptr = graph->sendptr.data();

  int count = 0;
  for(int i=0; i < ctrl->ncpus; ++i) {
    if( send_reqs[i].key > 0 ) {
      sendptr[count]   = send_reqs[i].key;
      startsind[count] = send_reqs[i].val;
      ++count;
    }
  }
  ASSERT(count == nnbrs);

  for(int i=1; i < nnbrs; ++i) sendptr[i] += sendptr[i-1];
  for(int i=nnbrs; i > 0; --i) sendptr[i]  = sendptr[i-1];
  sendptr[0] = 0;

  graph->nsend   = sendptr[nnbrs];
  graph->sendind.resize(graph->nsend);
  int *sendind = graph->sendind.data();

  // Issue the receives first
  for(int i=0; i < nnbrs; ++i) {
    MPI_Irecv(&sendind[sendptr[i]], sendptr[i+1]-sendptr[i], MPI_INT,
              graph->peind[i], 1, ctrl->comm, &ctrl->rreqs[i]);
  }

  // Issue the sends next
  for(int i=0; i < nnbrs; ++i) {
    MPI_Isend(&recvind[recvptr[i]], recvptr[i+1]-recvptr[i], MPI_INT,
              graph->peind[i], 1, ctrl->comm, &ctrl->sreqs[i]);
  }

  MPI_Waitall(nnbrs, &ctrl->rreqs[0], &ctrl->wstats[0]);
  MPI_Waitall(nnbrs, &ctrl->sreqs[0], &ctrl->wstats[0]);

  // Create the peadjncy data structure for sparse boundary exchanges
  graph->peadjncy.resize(graph->nsend);
  graph->peadjloc.resize(graph->nsend);
  graph->pexadj.assign(nvtxs+1, 0);
  int *pexadj = graph->pexadj.data();

  for(int i=0; i< graph->nsend; ++i) {
    ++pexadj[sendind[i] - first_vtx];
  }

  for(int i=1; i < nvtxs; ++i) pexadj[i] += pexadj[i-1];
  for(int i=nvtxs; i > 0; --i) pexadj[i]  = pexadj[i-1];
  pexadj[0] = 0;

  for(int i=0; i < nnbrs; ++i) {
    for(int j=sendptr[i]; j < sendptr[i+1]; ++j) {
      int k = pexadj[sendind[j] - first_vtx]++;
      graph->peadjncy[k] = i;  // peind[i] is the actual cpu number
      graph->peadjloc[k] = startsind[i]++;
    }
  }
  ASSERT(pexadj[nvtxs] == graph->nsend);

  for(int i=nvtxs; i > 0; --i) pexadj[i] = pexadj[i-1];
  pexadj[0] = 0;

  // Create the inverse map from ladjncy to adjncy
  graph->imap.resize(nvtxs+nrecv);
  for(int i=0; i < nvtxs; ++i)
    graph->imap[i] = first_vtx + i;
  for(int i=0; i < nrecv; ++i)
    graph->imap[nvtxs+i] = recvind[i];
}


//------------------------------------------------------------------------
// Perform distribution for the boundary vertices
void distribute_bdry_data(Control *ctrl, Graph *graph, int *data, 
			  int *recvvector)
{
  int first_vtx = graph->vtxdist[ctrl->myrank];
  int *sendptr  = graph->sendptr.data();
  int *recvptr  = graph->recvptr.data();

  // Issue the receives first
  for(int i=0; i < graph->nnbrs; ++i) {
    MPI_Irecv(&recvvector[recvptr[i]], recvptr[i+1]-recvptr[i], MPI_INT,
              graph->peind[i], 1, ctrl->comm, &ctrl->rreqs[i]);
  }

  // Issue the sends next
  int k = sendptr[graph->nnbrs];
  std::vector<int> send_vec(k);
  for(int i=0; i < k; ++i) 
    send_vec[i] = data[graph->sendind[i] - first_vtx];

  for(int i=0; i < graph->nnbrs; ++i) {
    MPI_Isend(&send_vec[sendptr[i]], sendptr[i+1]-sendptr[i], MPI_INT,
              graph->peind[i], 1, ctrl->comm, &ctrl->sreqs[i]); 
  }

  // Wait for the operations to finish
  MPI_Waitall(graph->nnbrs, &ctrl->rreqs[0], &ctrl->wstats[0]);
  MPI_Waitall(graph->nnbrs, &ctrl->sreqs[0], &ctrl->wstats[0]);
}


//------------------------------------------------------------------------
// Perform the exchange of the boundary vertices
void exchange_bdry_data(Control *ctrl, Graph *graph, int nchanged, 
			int *changed, int *data,
			ikv_t *sendpairs, ikv_t *recvpairs)
{
  int first_vtx = graph->vtxdist[ctrl->myrank];
  int *sendptr  = graph->sendptr.data();
  int *recvptr  = graph->recvptr.data();

  // Issue the receives first
  for(int i=0; i < graph->nnbrs; ++i) {
    MPI_Irecv(&recvpairs[recvptr[i]], 2*(recvptr[i+1]-recvptr[i]), MPI_INT, 
              graph->peind[i], 1, ctrl->comm, &ctrl->rreqs[i]);
  }

  if( nchanged != 0 ) {
    std::vector<int> psendptr(sendptr, sendptr + graph->nnbrs);

    // Copy the changed values into the sendvector
    for(int i=0; i < nchanged; ++i) {
      int j = changed[i];
      for(int k=graph->pexadj[j]; k < graph->pexadj[j+1]; ++k) {
        int icpu = graph->peadjncy[k];
        sendpairs[psendptr[icpu]].key = graph->peadjloc[k];
        sendpairs[psendptr[icpu]].val = data[j];
        ++psendptr[icpu];
      }
    }

    for(int i=0; i < graph->nnbrs; ++i) {
      MPI_Isend((sendpairs+sendptr[i]), 2*(psendptr[i]-sendptr[i]), MPI_INT, 
                graph->peind[i], 1, ctrl->comm, &ctrl->sreqs[i]);
    }
  } else {
    for(int i=0; i < graph->nnbrs; ++i) 
      MPI_Isend(sendpairs, 0, MPI_INT, graph->peind[i], 1, ctrl->comm, &ctrl->sreqs[i]);
  }

  // Wait for the operations to finish
  for(int i=0; i < graph->nnbrs; ++i) {
    MPI_Wait(&ctrl->rreqs[i], &(ctrl->status));
    int nreceived;
    MPI_Get_count(&ctrl->status, MPI_INT, &nreceived);
    if( nreceived != 0 ) {
      nreceived /= 2;
      ikv_t *pairs = recvpairs + graph->recvptr[i];
      for(int k=0; k < nreceived; ++k) 
        data[pairs[k].key] = pairs[k].val;
    }
  }

  MPI_Waitall(graph->nnbrs, &ctrl->sreqs[0], &ctrl->wstats[0]);
}

//------------------------------------------------------------------------
// Randomly permute the contents of an array
//   flag == 0, don't initialize perm
//   flag == 1, set p[i] = i 
void random_permute(int n, int *p, int flag)
{
  if( flag == 1 ) {
    for(int i=0; i < n; i++)
      p[i] = i;
  }

  for(int i=0; i < n; ++i) {
    int v = RandomInRange(n);
    int u = RandomInRange(n);
    std::swap(p[v], p[u]);
  }
}

//------------------------------------------------------------------------
// Randomly permute the contents of an array
//   flag == 0, don't initialize perm
//   flag == 1, set p[i] = i
void random_permute_fast(int n, int *p, int flag)
{
  // Fo small arrays
  if( n < 25 ) {
    random_permute(n, p, flag);
    return;
  }

  if( flag == 1 ) {
    for(int i=0; i < n; ++i)
      p[i] = i;
  }

  for(int i=0; i < n; i+=8) {
    int v = RandomInRange(n-4);
    int u = RandomInRange(n-4);
    std::swap(p[v],   p[u]);
    std::swap(p[v+1], p[u+1]);
    std::swap(p[v+2], p[u+2]);
    std::swap(p[v+3], p[u+3]);
  }
}

//------------------------------------------------------------------------
// Check if the pairwise balance between two partitions will improve by
// moving the vertex v from pfrom to pto, accounting for partition weights
int check_better_balance_wt(int ncon, double *pfrom, double *pto,
                              double *nvwgt, double *disbalance)
{
  double blb1 = 0., alb1 = 0.;
  double blb2 = 0., alb2 = 0.;
  double sblb = 0., salb = 0.;

  for(int i=0; i < ncon; ++i) {
    double temp = std::max(pfrom[i], pto[i]) / disbalance[i];
    if( blb1 < temp ) {
      blb2 = blb1;
      blb1 = temp;
    } else if( blb2 < temp ) {
      blb2 = temp;
    }
    sblb += temp;

    temp = std::max(pfrom[i]-nvwgt[i], pto[i]+nvwgt[i]) / disbalance[i];
    if( alb1 < temp ) {
      alb2 = alb1;
      alb1 = temp;
    } else if( alb2 < temp ) {
      alb2 = temp;
    }
    salb += temp;
  }

  if( alb1 < blb1 ) return 1;
  if( blb1 < alb1 ) return 0;
  if( alb2 < blb2 ) return 1;
  if( blb2 < alb2 ) return 0;

  return salb < sblb;
}


//------------------------------------------------------------------------
// Check if it will be better to move a vertex to pt2 than to pt1 subject
// to their target weights of pt1 and pt2, respectively
int check_better_balance_pt(int ncon, double *pt1, double *pt2,
                              double *nvwgt, double *disbalance)
{
  double m11 = 0., m12 = 0.;
  double m21 = 0., m22 = 0.;
  double sm1 = 0., sm2 = 0.;

  for(int i=0; i < ncon; ++i) {
    double temp = (pt1[i] + nvwgt[i]) / disbalance[i];
    if( m11 < temp ) {
      m12 = m11;
      m11 = temp;
    } else if( m12 < temp ) {
      m12 = temp;
    }
    sm1 += temp;

    temp = (pt2[i] + nvwgt[i]) / disbalance[i];
    if( m21 < temp ) {
      m22 = m21;
      m21 = temp;
    } else if (m22 < temp) {
      m22 = temp;
    }
    sm2 += temp;
  }

  if( m21 < m11 ) return 1;
  if( m21 > m11 ) return 0;
  if( m22 < m12 ) return 1;
  if( m22 > m12 ) return 0;

  return sm2 < sm1;
}

//------------------------------------------------------------------------
// Compute the balance of the serial partitioning
void get_balance_local(Control *ctrl, Graph *graph, int *where, double *disbalance)
{
  int ncon = graph->ncon; // convinience shortcut

  std::vector<int> pwgts (ncon*ctrl->nparts, 0);
  std::vector<int> tvwgts(ncon, 0);

  for(int i=0; i < graph->nvtxs; ++i) {
    for(int j=0; j < ncon; ++j) {
      pwgts[where[i]*ncon+j] += graph->vwgt[i*ncon+j];
      tvwgts[j] += graph->vwgt[i*ncon+j];
    }
  }

  for(int j=0; j < ncon; ++j) {
    double maximb = 0.;
    for(int i=0; i < ctrl->nparts; ++i)
      maximb =std::max(maximb,
		       (1. + pwgts[i*ncon+j]) /
		       (1. + (ctrl->part_weights[i*ncon+j] * tvwgts[j])));
    disbalance[j] = maximb;
  }
}

//------------------------------------------------------------------------
// Compute the balance of the distributed partitioning
void get_balance_global(Control *ctrl, Graph *graph, int *where, double *disbalance)
{
  int ncon = graph->ncon;  // convinience shortcut

  std::vector<double> lminvwgts(ncon, 1.);
  std::vector<double> gminvwgts(ncon);
  std::vector<double> lnpwgts  (ncon*ctrl->nparts, 0.);
  std::vector<double> gnpwgts  (ncon*ctrl->nparts);

  for(int i=0; i < graph->nvtxs; ++i) {
    for(int j=0; j < ncon; ++j) {
      lnpwgts[where[i]*ncon+j] += graph->nvwgt[i*ncon+j];

      // Special care for part_weights[] == 0.
      if( graph->nvwgt[i*ncon+j] > 0. && lminvwgts[j] > graph->nvwgt[i*ncon+j] )
	lminvwgts[j] = graph->nvwgt[i*ncon+j];
    }
  }

  MPI_Allreduce(lnpwgts.data(),   gnpwgts.data(),   ncon*ctrl->nparts, MPI_DOUBLE, MPI_SUM, ctrl->comm);
  MPI_Allreduce(lminvwgts.data(), gminvwgts.data(), ncon,              MPI_DOUBLE, MPI_MIN, ctrl->comm);

  for(int j=0; j < ncon; ++j) {
    double maximb = 0.;
    for(int i=0; i < ctrl->nparts; ++i)
      maximb = std::max(maximb,
			(gminvwgts[j] + gnpwgts[i*ncon+j]) /
			(gminvwgts[j] + ctrl->part_weights[i*ncon+j]));
    disbalance[j] = maximb;
  }
}

//------------------------------------------------------------------------
// Perform refinement
void refine_adaptive_kw(Control *ctrl, Graph *graph, int npasses)
{
  int  nvtxs   = graph->nvtxs;   // convinience shortcut
  int  ncon    = graph->ncon;    // convinience shortcut
  int *xadj    = graph->xadj.data();    // convinience shortcut
  int *sendptr = graph->sendptr.data(); // convinience shortcut
  int  nparts  = ctrl->nparts;   // convinience shortcut

  int first_vtx = graph->vtxdist[ctrl->myrank];

  // Allocate working space
  std::vector<double> lbvec(ncon);
  std::vector<double> badmaxpwgt(ncon*nparts);
  std::vector<double> movewgts  (ncon*nparts);
  std::vector<double> ognpwgts  (ncon*nparts);
  std::vector<double> pgnpwgts  (ncon*nparts);
  std::vector<double> overfill  (ncon*nparts);

  std::vector<int> nupd_cpu (ctrl->ncpus);
  std::vector<int> pperm    (nparts);
  std::vector<int> old_ed   (nvtxs);
  std::vector<int> changed  (nvtxs);
  std::vector<int> perm     (nvtxs);
  std::vector<int> update   (nvtxs);
  std::vector<int> moved    (nvtxs);
  std::vector<int> htable   (nvtxs + graph->nrecv, 0);
  std::vector<int> where_tmp(graph->where);
  std::vector<int> supdate  (graph->nrecv);
  std::vector<int> rupdate  (graph->nsend);

  std::vector<ikv_t> rwchanges(graph->nrecv);
  std::vector<ikv_t> swchanges(graph->nsend);

  for(int i=0; i < nparts; ++i) {
    for(int j=0; j < ncon; ++j) 
      badmaxpwgt[i*ncon+j] = ctrl->disbalance[j] * ctrl->part_weights[i*ncon+j];
  }

  // Store external degrees of the vertices before refinement
  for(int i=0; i < nvtxs; ++i)
    old_ed[i] = graph->ckrinfo[i].ed;

  //-------------------------------------------------------------
  // Perform several walk ins through the vertices
  for(int pass=0; pass < npasses; ++pass) {
    int oldcut = graph->mincut;
    if (ctrl->myrank == 0)
      random_permute(nparts, pperm.data(), 1);

    MPI_Bcast(pperm.data(), nparts, MPI_INT, 0, ctrl->comm);

    // Move dirty vertices first 
    int ndirty = 0;
    for(int i=0; i < nvtxs; ++i)
      if( graph->where[i] != ctrl->myrank )
        ++ndirty;

    int dptr = 0;
    for(int i=0; i < nvtxs; ++i)
      if( graph->where[i] != ctrl->myrank)
        perm[dptr++] = i;
      else
        perm[ndirty++] = i;

    ASSERT(ndirty == nvtxs);
    ndirty = dptr;
    int nclean = nvtxs - dptr;
    random_permute_fast(ndirty, &perm[0], 0);
    if( nclean > 0 )
    random_permute_fast(nclean, &perm[ndirty], 0);

    // Compute unbalance of the current partition
    get_balance_global(ctrl, graph, graph->where.data(), lbvec.data());
    double ubavg = average(ncon, ctrl->disbalance.data());
    double lbavg = average(ncon, lbvec.data());
    bool imbalanced = (lbavg > ubavg);

    for(int step=0; step < 2; ++step) {
      rcopy(ncon*nparts, graph->gnpwgts.data(), ognpwgts.data());
      rset(ncon*nparts, 0.0, movewgts.data());
      int nmoved = 0;

      // Step 1 - record staistics for desired moves
      int to = -1;
      for(int ii=0; ii < nvtxs; ++ii) {
        int i     = perm[ii];
        int from  = where_tmp[i];
        double *nvwgt = &graph->nvwgt[i*ncon];

        double vsize = (double)(graph->vsize[i]);
	bool weight_coincide = false;
        for(int jj=0; jj < ncon; ++jj) {
          if( std::abs(nvwgt[jj] - graph->gnpwgts[from*ncon+jj]) < SMALL_FLOAT ) {
	    weight_coincide = true;
            break;
	  }
        }

        if( weight_coincide )
          continue;

        // Check for possible improvement
        ckrinfo_t *myrinfo = &graph->ckrinfo[i];
        if( myrinfo->ed <= 0 )
          continue;

        ASSERT(myrinfo->inbr != -1);
        std::pair<int,int> *mynbrs = ctrl->pool.data() + myrinfo->inbr;

	int k;
        for(k=myrinfo->nnbrs-1; k >= 0; --k) {
          to = mynbrs[k].first;
	  bool proper_side = ((step == 0) ? (pperm[from] - pperm[to] < 0)
                                          : (pperm[from] - pperm[to] > 0));
	  if( proper_side ) {
	    int jj;
            for(jj=0; jj < ncon; ++jj) {
              if( graph->gnpwgts[to*ncon+jj]+nvwgt[jj] > badmaxpwgt[to*ncon+jj] && nvwgt[jj] > 0. )
                break;
            }

            if( jj == ncon )
              break;
          }
        }

        // Skip the rest if no candidate found
        if( k < 0 )
          continue;

        int oldto = to;

	bool better, worse;
        if( ctrl->is_coupled ) {
          better = (oldto == ctrl->myrank);
          worse  = (from  == ctrl->myrank);
	} else {
          better = (oldto == graph->home[i]);
          worse  = (from  == graph->home[i]);
	}

        double oldgain = ctrl->itr_ratio * (mynbrs[k].second - myrinfo->id);
        if( better ) oldgain += ctrl->redist_ratio * vsize;
        if( worse )  oldgain -= ctrl->redist_ratio * vsize;

	double gain;
        for(int j=k-1; j >= 0; --j) {
          to = mynbrs[j].first;
	  bool proper_side = ((step == 0) ? (pperm[from] - pperm[to] < 0)
                                          : (pperm[from] - pperm[to] > 0));
	  if( proper_side ) {
            if( ctrl->is_coupled )
              better = (to == ctrl->myrank);
	    else
              better = (to == graph->home[i]);

            gain = ctrl->itr_ratio * (mynbrs[j].second - myrinfo->id);
            if( better ) gain += ctrl->redist_ratio * vsize;
            if( worse )  gain -= ctrl->redist_ratio * vsize;

	    int jj;
            for(jj=0; jj < ncon; ++jj) {
              if( graph->gnpwgts[to*ncon+jj]+nvwgt[jj] > badmaxpwgt[to*ncon+jj] && nvwgt[jj] > 0. )
                break;
            }

            if( jj == ncon ) {
              if( gain > oldgain || 
                  (std::abs(gain-oldgain) < SMALL_FLOAT &&
                   check_better_balance_pt(ncon, &graph->gnpwgts[oldto*ncon],
                                           &graph->gnpwgts[to*ncon],
                                           nvwgt, ctrl->disbalance.data())) ) {
                oldgain = gain;
                oldto   = to;
                k       = j;
              }
            }
          }
        }
        to   = oldto;
        gain = oldgain;

        if( gain > 0. ||
            (gain > -SMALL_FLOAT &&
             (imbalanced ||  graph->level > 3  || ii % 8 == 0) &&
             check_better_balance_wt(ncon, &graph->gnpwgts[from*ncon],
                                     &graph->gnpwgts[to*ncon],
                                     nvwgt, ctrl->disbalance.data())) ) {
          // Update temp arrays of the moved vertex 
          where_tmp[i] = to;
          moved[nmoved++] = i;
          for(int jj=0; jj < ncon; ++jj) {
	    graph->lnpwgts[to*ncon+jj]   += nvwgt[jj];
	    graph->lnpwgts[from*ncon+jj] -= nvwgt[jj];
	    graph->gnpwgts[to*ncon+jj]   += nvwgt[jj];
	    graph->gnpwgts[from*ncon+jj] -= nvwgt[jj];
	    movewgts[to*ncon+jj]  += nvwgt[jj];
	    movewgts[from*ncon+jj]-= nvwgt[jj];
          }

          myrinfo->ed += myrinfo->id - mynbrs[k].second;
          std::swap(myrinfo->id, mynbrs[k].second);
          if( mynbrs[k].second == 0 )
            mynbrs[k] = mynbrs[--myrinfo->nnbrs];
          else
            mynbrs[k].first = from;

          // Update the degrees of adjacent vertices
          for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
            // Skip verts from other cpus
            if( graph->adjncy[j] >= nvtxs )
              continue;

            int me = graph->adjncy[j];
            int this_part = where_tmp[me];

            myrinfo = &graph->ckrinfo[me];
            if( myrinfo->inbr == -1 ) {
              myrinfo->inbr  = ctrl->pool.get_next(graph->xadj[me+1] - graph->xadj[me]+1);
              myrinfo->nnbrs = 0;
            }
            mynbrs = ctrl->pool.data() + myrinfo->inbr;

            if( this_part == from ) {
	      myrinfo->ed += graph->adjwgt[j];
	      myrinfo->id -= graph->adjwgt[j];
            } else {
              if( this_part == to ) {
		myrinfo->id += graph->adjwgt[j];
		myrinfo->ed -= graph->adjwgt[j];
              }
            }

            // Remove edge weight in 'from'
            if( this_part != from ) {
              for(int k=0; k < myrinfo->nnbrs; ++k) {
                if( mynbrs[k].first == from ) {
                  if( mynbrs[k].second == graph->adjwgt[j] )
                    mynbrs[k] = mynbrs[--myrinfo->nnbrs];
                  else
                    mynbrs[k].second -= graph->adjwgt[j];
                  break;
                }
              }
            }

            // Add edge weight in 'to'
            if( this_part != to ) {
	      bool exist = false;
              for(int k=0; k < myrinfo->nnbrs; ++k) {
                if( mynbrs[k].first == to ) {
                  mynbrs[k].second += graph->adjwgt[j];
		  exist = true;
                  break;
                }
              }
	      if( !exist ) {
                mynbrs[myrinfo->nnbrs].first = to;
                mynbrs[myrinfo->nnbrs].second= graph->adjwgt[j];
                myrinfo->nnbrs++;
              }
            }
          }
        }
      }

      //-------------------------------------------------------------
      // After commiting all transfers propagate subdomain weights
      MPI_Allreduce(graph->lnpwgts.data(), pgnpwgts.data(), nparts*ncon,
                    MPI_DOUBLE, MPI_SUM, ctrl->comm);

      // Compute overfill array
      bool is_overweight = false;
      for(int j=0; j < nparts; ++j) {
        for(int jj=0; jj < ncon; ++jj) {
          if( pgnpwgts[j*ncon+jj] > ognpwgts[j*ncon+jj] ) {
            overfill[j*ncon+jj] = (pgnpwgts[j*ncon+jj] - badmaxpwgt[j*ncon+jj]) /
                                  (pgnpwgts[j*ncon+jj] - ognpwgts  [j*ncon+jj]);
          } else {
            overfill[j*ncon+jj] = 0.;
          }

          overfill[j*ncon+jj]  = std::max(overfill[j*ncon+jj], 0.);
          overfill[j*ncon+jj] *= movewgts[j*ncon+jj];

          if( overfill[j*ncon+jj] > 0. )
            is_overweight = true;
        }
      }

      // Based on overfill values undo some transfers
      if( is_overweight ) {
        for(int ii=0; ii < nmoved; ++ii) {
          int i     = moved[ii];
          int oldto = where_tmp[i];
          double *nvwgt = &graph->nvwgt[i*ncon];

          ckrinfo_t *myrinfo = &graph->ckrinfo[i];
          std::pair<int,int> *mynbrs = ctrl->pool.data() + myrinfo->inbr;
          ASSERT(myrinfo->nnbrs == 0 || myrinfo->inbr != -1);

	  int k;
          for(k=0; k < myrinfo->nnbrs; ++k) {
            if( mynbrs[k].first == graph->where[i] )
              break;
          }

	  int jj;
          for(jj=0; jj < ncon; ++jj) {
            if( nvwgt[jj] > 0. && overfill[oldto*ncon+jj] > nvwgt[jj]/4. )
              break;
          }

          // Cancel this transfer if necessary
          if( k != myrinfo->nnbrs && jj != ncon ) {
            moved[ii]  = -1;
            int from = oldto;
            to         = graph->where[i];

            for(jj=0; jj < ncon; ++jj) 
              overfill[oldto*ncon+jj] = std::max(overfill[oldto*ncon+jj] - nvwgt[jj], 0.);

            where_tmp[i] = to;
            myrinfo->ed += myrinfo->id - mynbrs[k].second;
            std::swap(myrinfo->id, mynbrs[k].second);
            if (mynbrs[k].second == 0)
              mynbrs[k] = mynbrs[--myrinfo->nnbrs];
            else
              mynbrs[k].first = from;

            for(jj=0; jj < ncon; ++jj) {
	      graph->lnpwgts[to*ncon+jj]   += nvwgt[jj];
	      graph->lnpwgts[from*ncon+jj] -= nvwgt[jj];
            }

            // Update the degrees of adjacent vertices
            for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
              // Skip vertices on other cpus
              if( graph->adjncy[j] >= nvtxs )
                continue;

              int me = graph->adjncy[j];
              int this_part = where_tmp[me];

              myrinfo = &graph->ckrinfo[me];
              if( myrinfo->inbr == -1 ) {
                myrinfo->inbr  = ctrl->pool.get_next(graph->xadj[me+1] - graph->xadj[me]+1);
                myrinfo->nnbrs = 0;
              }
              mynbrs = ctrl->pool.data() + myrinfo->inbr;

              if( this_part == from ) {
		myrinfo->ed += graph->adjwgt[j];
		myrinfo->id -= graph->adjwgt[j];
              } else {
                if( this_part == to ) {
		  myrinfo->id += graph->adjwgt[j];
		  myrinfo->ed -= graph->adjwgt[j];
                }
              }

	      // Remove edge weight in 'from'
              if( this_part != from ) {
                for(int k=0; k < myrinfo->nnbrs; ++k) {
                  if( mynbrs[k].first == from ) {
                    if( mynbrs[k].second == graph->adjwgt[j] )
                      mynbrs[k] = mynbrs[--myrinfo->nnbrs];
                    else 
                      mynbrs[k].second -= graph->adjwgt[j];
                    break;
                  }
                }
              }

	      // Add edge weight i 'to'
              if( this_part != to ) {
       		bool exist = false;
                for(int k=0; k < myrinfo->nnbrs; ++k) {
                  if( mynbrs[k].first == to ) {
                    mynbrs[k].second += graph->adjwgt[j];
		    exist = true;
                    break;
                  }
                }
                if( !exist ) {
                  mynbrs[myrinfo->nnbrs].first = to;
                  mynbrs[myrinfo->nnbrs].second= graph->adjwgt[j];
                  myrinfo->nnbrs++;
                }
              }
            }
          }
        }
      }

      // Step 2 - store the rest info for the moves
      int nlupd=0, nsupd=0, nmoves=0, nchanged=0;
      for(int ii=0; ii < nmoved; ++ii) {
        int i = moved[ii];
        if( i == -1 )
          continue;

        graph->where[i] = where_tmp[i]; 

        // Update the vertex information
        if( htable[i] == 0 ) {
          htable[i] = 1;
          update[nlupd++] = i;
        }

        // Put the vertices adjacent to i into the update array
        for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
          int k = graph->adjncy[j];
          if( htable[k] == 0 ) {
            htable[k] = 1;
            if( k < nvtxs )
              update[nlupd++] = k;
            else
              supdate[nsupd++] = k;
          }
        }
        ++nmoves;

        if( graph->pexadj[i+1] - graph->pexadj[i] > 0 )
          changed[nchanged++] = i;
      }

      // Propagate to neigboring cpus new where[] info for the interface vertices
      exchange_bdry_data(ctrl, graph, nchanged,
			 changed.data(), graph->where.data(),
			 swchanges.data(), rwchanges.data()); 

      //-------------------------------------------------------------
      // Send to other cpus the verts which degrees are to be updated
      // Issue the receives first
      for(int i=0; i < graph->nnbrs; ++i) {
        MPI_Irecv(&rupdate[sendptr[i]], sendptr[i+1]-sendptr[i], MPI_INT,
		  graph->peind[i], 1, ctrl->comm, &ctrl->rreqs[i]);
      }

      // Issue the sends next
      for(int i=0; i < nsupd; ++i) {
        htable[supdate[i]] = 0;
        supdate[i] = graph->imap[supdate[i]];
      }
      isorti(nsupd, supdate.data());

      for(int k=0, j=0, i=0; i < graph->nnbrs; ++i) {
        int mylast_vtx = graph->vtxdist[graph->peind[i]+1];
        for(k=j; (k < nsupd) && (supdate[k] < mylast_vtx); ++k); 
        MPI_Isend(&supdate[j], k-j, MPI_INT, graph->peind[i], 1, ctrl->comm, &ctrl->sreqs[i]);
        j = k;
      }

      // Wait for the operations to finish
      MPI_Waitall(graph->nnbrs, &ctrl->rreqs[0], &ctrl->wstats[0]);
      for(int i=0; i < graph->nnbrs; ++i) 
        MPI_Get_count(&ctrl->wstats[i], MPI_INT, &nupd_cpu[i]);
      MPI_Waitall(graph->nnbrs, &ctrl->sreqs[0], &ctrl->wstats[0]);

      //-------------------------------------------------------------
      // Place the recieved marked vertices into update[] array
      for(int i=0; i < graph->nnbrs; ++i) {
        int *pe_updates = &rupdate[sendptr[i]];
        for(int j=0; j < nupd_cpu[i]; ++j) {
          int k = pe_updates[j];
          if( htable[k - first_vtx] == 0 ) {
            htable[k - first_vtx] = 1;
            update[nlupd] = k - first_vtx;
	    ++nlupd;
          }
        }
      }

      //-------------------------------------------------------------
      // Update information of the vertices in the update[] array
      for(int ii=0; ii < nlupd; ++ii) {
        int i = update[ii];
        ASSERT(htable[i] == 1);

        htable[i] = 0;

        int this_part = graph->where[i];
        ckrinfo_t *myrinfo = &graph->ckrinfo[i];

        if( myrinfo->inbr == -1 )
          myrinfo->inbr  = ctrl->pool.get_next(graph->xadj[i+1] - graph->xadj[i]+1);
        std::pair<int,int> *mynbrs = ctrl->pool.data() + myrinfo->inbr;

        graph->lmincut -= old_ed[i];
        myrinfo->nnbrs  = 0;
        myrinfo->id     = 0;
        myrinfo->ed     = 0;

        for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
          int other_part = graph->where[graph->adjncy[j]];
          if( this_part != other_part ) {
            myrinfo->ed += graph->adjwgt[j];

	    bool exist = false;
            for(int k=0; k < myrinfo->nnbrs; ++k) {
              if( mynbrs[k].first == other_part ) {
                mynbrs[k].second += graph->adjwgt[j];
		exist = true;
                break;
              }
            }
            if( !exist ) {
              mynbrs[myrinfo->nnbrs].first = other_part;
              mynbrs[myrinfo->nnbrs].second= graph->adjwgt[j];
              myrinfo->nnbrs++;
            }
            ASSERT(myrinfo->nnbrs <= graph->xadj[i+1] - graph->xadj[i]);
          } else {
            myrinfo->id += graph->adjwgt[j];
          }
        }
        graph->lmincut += myrinfo->ed;
        old_ed[i]       = myrinfo->ed; // for the next iteration
      }

      // Get sump of the partition weights
      MPI_Allreduce(graph->lnpwgts.data(), graph->gnpwgts.data(), nparts*ncon,
                    MPI_DOUBLE, MPI_SUM, ctrl->comm);
    }

    MPI_Allreduce(&graph->lmincut, &graph->mincut, 1, MPI_INT, MPI_SUM, ctrl->comm);
    graph->mincut /= 2;

    if( graph->mincut == oldcut )
      break;
  }
}

//------------------------------------------------------------------------
// Return true if target weights are similar
bool are_tweights_similar(const double *target_weights, int ncon, int s1, int s2)
{
  static const double small_real = 1.e-6;

  for(int i=0; i < ncon; ++i)
    if( std::abs(target_weights[s1*ncon+i]-target_weights[s2*ncon+i]) > small_real )
      return false;

  return true;
}

//------------------------------------------------------------------------
// Compute the assignment minimizing of the total volume of moved data
void remap_total_volume(const Control& ctrl, const std::vector<int>& lpwgts,
                        std::vector<int>& map, int npasses, int ncon)
{
  int nparts = ctrl.nparts; // convinience shortcut

  std::vector<int> rowmap(nparts, -1);
  std::vector<int> mylpwgts(lpwgts);
  std::vector<ikv_t> recv(nparts);

  int nmapped = 0, done = 0;
  ikv_t send;
  for(int pass=0; pass < npasses; ++pass) {
    int maxipwgt = (int)iargmax(nparts, mylpwgts.data(), 1);

    if( mylpwgts[maxipwgt] > 0 && !done ) {
      send.key = -mylpwgts[maxipwgt];
      send.val = ctrl.myrank*nparts + maxipwgt;
    } else {
      send.key = 0;
      send.val = -1;
    }

    // each processor sends its selection
    MPI_Allgather(&send, 2, MPI_INT, recv.data(), 2, MPI_INT, ctrl.comm); 

    ikvsorti(nparts, recv.data());
    if( recv[0].key == 0 )
      break;

    // now make as many assignments as possible
    for(int ii=0; ii < nparts; ++ii) {
      int i = recv[ii].val;

      if( i == -1 )
        continue;

      int j = i % nparts;
      int k = i / nparts;
      if( map[j] == -1 && rowmap[k] == -1 &&
          are_tweights_similar(ctrl.part_weights.data(), ncon, j, k) ) {
        map[j] = k;
        rowmap[k] = j;
        nmapped++;
        mylpwgts[j] = 0;
        if( ctrl.myrank == k )
          done = 1;
      }

      if( nmapped == nparts )
        break;
    }

    if( nmapped == nparts )
      break;
  }

  // Map unmapped partitions
  if( nmapped < nparts ) {
    for(int i=0, j=0; j < nparts && nmapped < nparts; ++j) {
      if( map[j] == -1 ) {
        for(; i < nparts; ++i) {
          if( rowmap[i] == -1 &&
              are_tweights_similar(ctrl.part_weights.data(), ncon, i, j) ) {
            map[j] = i;
            rowmap[i] = j;
            nmapped++;
            break;
          }
        }
      }
    }
  }

  // Check to see if remapping fails (due to dis-similar part_weights)
  // If remapping fails, revert to original mapping
  if( nmapped < nparts ) {
    for(int i=0; i < nparts; ++i)
      map[i] = i; 
  } else {
    // check for a savings
    int oldwgt = lpwgts[ctrl.myrank];
    int newwgt = lpwgts[rowmap[ctrl.myrank]];
    int nsaved  = newwgt - oldwgt;

    int gnsaved;
    MPI_Allreduce(&nsaved, &gnsaved, 1, MPI_INT, MPI_SUM, ctrl.comm);

    // undo everything if we don't see a savings
    if( gnsaved <= 0 ) {
      for(int i=0; i < nparts; ++i)
        map[i] = i;
    }
  }
}

//------------------------------------------------------------------------
Graph::Graph()
  : gnvtxs (-1),
    nvtxs  (-1),
    nedges (-1),
    ncon   (0),
    nnbrs  (-1),
    nrecv  (-1),
    nsend  (-1),
    nlocal (-1),
    lmincut(0),
    mincut (0),
    level  (0),
    global_match(false),
    coarser(NULL),
    finer  (NULL),
    xadj    (),
    adjncy  (),
    adjwgt  (),
    vtxdist (),
    home    (),
    vsize   (),
    vwgt    (),
    nvwgt   (),
    match   (),
    cmap    (),
    label   (),
    peind   (),
    sendptr (),
    sendind (),
    recvptr (),
    recvind (),
    imap    (),
    pexadj  (),
    peadjncy(),
    peadjloc(),
    lperm   (),
    rlens   (),
    slens   (),
    rcand   (),
    lnpwgts (),
    gnpwgts (),
    where   (),
    ckrinfo ()
{}

Graph::~Graph()
{
  xadj.clear();
  adjncy.clear();
  adjwgt.clear();
  vtxdist.clear();
  home.clear();
  vsize.clear();
  vwgt.clear();
  nvwgt.clear();

  match.clear();
  cmap.clear();
  label.clear();
  peind.clear();
  sendptr.clear();
  sendind.clear();
  recvptr.clear();
  recvind.clear();
  imap.clear();
  pexadj.clear();
  peadjncy.clear();
  peadjloc.clear();
  lperm.clear();
  rlens.clear();
  slens.clear();
  rcand.clear();
  lnpwgts.clear();
  gnpwgts.clear();
  where.clear();
  ckrinfo.clear();
}

//------------------------------------------------------------------------
// Create the nvwgt of the graph
void Graph::set_weights(const Control& ctrl)
{
  nvwgt.resize(nvtxs * ncon);
  for(int i=0; i < nvtxs; ++i) {
    for(int j=0; j < ncon; ++j)
      nvwgt[i*ncon+j] = ctrl.inv_weights[j] * vwgt[i*ncon+j];
  }
}

//------------------------------------------------------------------------
// Create the nvwgt of the graph
void Graph::remap(const Control& ctrl)
{
  static const int NREMAP_PASSES = 4;

  if( ctrl.ncpus != ctrl.nparts ) return;

  std::vector<int> lpwgts(ctrl.nparts, 0);

  if( !vsize.empty() ) {
    for(int i=0; i < nvtxs; ++i)
      lpwgts[where[i]] += vsize[i];
  } else {
    for(int i=0; i < nvtxs; ++i)
      lpwgts[where[i]] += 1;
  }

  std::vector<int> map(ctrl.nparts, -1);
  remap_total_volume(ctrl, lpwgts, map, NREMAP_PASSES, ncon);

  for(int i=0; i < nvtxs; ++i)
    where[i] = map[where[i]];
}

//------------------------------------------------------------------------
// Renumber connections list from local to global
void Graph::renumber_conns(int *adjncy_out)
{
  if( imap.empty() ) return;

  // Apply local to global index mapping
  for(int i=0; i < nedges; ++i)
    adjncy_out[i] = imap[adjncy[i]];
}

//------------------------------------------------------------------------
// Create a Graph data structure
Graph *create_graph()
{
  Graph *graph = new Graph();
  return graph;
}

//------------------------------------------------------------------------
// Create the graph from the inputs 
Graph *setup_graph(Control *ctrl, int ncon, int *vtxdist, int *xadj, 
		   int *vwgt, int *vsize, int *adjncy, int *adjwgt, 
		   int wgtflag)
{
  Graph *graph = create_graph();

  graph->level   = 0;
  graph->gnvtxs  = vtxdist[ctrl->ncpus];
  graph->nvtxs   = vtxdist[ctrl->myrank+1]-vtxdist[ctrl->myrank];
  graph->ncon    = ncon;
  graph->nedges  = xadj[graph->nvtxs];

  graph->xadj   .assign(xadj,    xadj + graph->nvtxs + 1);;
  graph->adjncy .assign(adjncy,  adjncy + graph->nedges);
  graph->vtxdist.assign(vtxdist, vtxdist + ctrl->ncpus + 1);

  // Allocate memory for weight arrays if not provided
  if( (wgtflag & 2) == 0 ||
      !vwgt )     graph->vwgt.assign(graph->nvtxs*ncon, 1);
  else            graph->vwgt.assign(vwgt, vwgt + graph->nvtxs*ncon);

  if( (wgtflag & 1) == 0 ||
      !adjwgt )   graph->adjwgt.assign(graph->nedges, 1);
  else            graph->adjwgt.assign(adjwgt, adjwgt + graph->nedges);

  // Allocate memory for special arrays
  if( vsize )     graph->vsize.assign(vsize, vsize + graph->nvtxs);
  else            graph->vsize.assign(graph->nvtxs, 1);

  graph->home.assign(graph->nvtxs, 1);

  int edge_weight_loc = isum(graph->nedges, graph->adjwgt.data(), 1);
  int node_weight_loc = isum(graph->nvtxs,  graph->vsize.data(),  1);

  int edge_weight_glb, node_weight_glb;
  MPI_Allreduce(&edge_weight_loc, &edge_weight_glb, 1, MPI_INT, MPI_SUM, ctrl->comm);
  MPI_Allreduce(&node_weight_loc, &node_weight_glb, 1, MPI_INT, MPI_SUM, ctrl->comm);

  ctrl->edge_node_ratio = (edge_weight_glb + 0.01) / (node_weight_glb + 0.01);

  // Compute invtvwgts
  ctrl->set_inv_weights(*graph);

  // Compute nvwgts
  graph->set_weights(*ctrl);

  return graph;
}

//------------------------------------------------------------------------
// Deallocate any memory stored in a graph
void free_graph(Graph *graph) 
{
  delete graph;
}


//------------------------------------------------------------------------
// Keep graph parts
void keep_parts(Graph *graph, int *part, int mypart)
{
  std::vector<int> reindex(graph->nvtxs);
 
  for(int count=0, i=0; i < graph->nvtxs; ++i) {
    if( part[i] == mypart )
      reindex[i] = count++;
  }

  int nv = 0, ne = 0;
  for(int j=graph->xadj[0], i=0; i < graph->nvtxs; ++i) {
    if( part[i] == mypart ) {
      for(; j < graph->xadj[i+1]; ++j) {
        int k = graph->adjncy[j];
        if( part[k] == mypart ) {
          graph->adjncy[ne]   = reindex[k];
          graph->adjwgt[ne++] = graph->adjwgt[j];
        }
      }
      j = graph->xadj[i+1];  // Save xadj[i+1] for later use

      for(int k=0; k < graph->ncon; ++k)
        graph->vwgt[nv*graph->ncon+k] = graph->vwgt[i*graph->ncon+k];

      graph->label[nv]  = graph->label[i];
      graph->xadj[++nv] = ne;
    } else {
      j = graph->xadj[i+1];  // Save xadj[i+1] for later use
    }
  }

  graph->nvtxs  = nv;
  graph->nedges = ne;
}



//------------------------------------------------------------------------
// Sort the vwgts of a vertex and return a hashed value 
int get_rank_hash(int ncon, const int *vwgt)
{
  int multiplier = 1;
  int retval = 0;
  for(int i=0; i < ncon; ++i) {
    multiplier *= (i+1);
    retval += vwgt[ncon-1-i] * multiplier;
  }

  return retval;
}

//------------------------------------------------------------------------
// Queue selection is hardcoded for up to 4 constraints 
void select_queue_dynamic(int nqueues, int ncon, int subdomain1, 
                          int subdomain2, int *currentq, double *flows,
                          int *from, int *qnum, int minval,
                          double avgvwgt, double maxdiff)
{
  static const double BALANCE_THRESHOLD = 0.2;

  *qnum = -1;

  // Allocate memory
  std::vector<int> not_needed(ncon, 1);
  std::vector<rkv_t> array(ncon);

  double sign = 0.;
  int index = -1;
  if( *from == -1 ) {
    for(int i=0; i < ncon; ++i) {
      array[i].key = std::abs(flows[i]);
      array[i].val = i;
    }

    // May be need to check the correct direction of the sort
    rkvsorti(ncon, array.data());

    if( flows[array[ncon-1].val] > avgvwgt ) {
      *from = subdomain1;
      sign  = 1.0;
      index = 0;
    }

    if( flows[array[ncon-1].val] < -avgvwgt ) {
      *from = subdomain2;
      sign  = -1.0;
      index = nqueues;
    }

    if( *from == -1 )
      return;
  }

  assert(*from == subdomain1 || *from == subdomain2);

  if( *from == subdomain1 ) {
    sign  = 1.;
    index = 0;
  } else {
    sign  = -1.;
    index = nqueues;
  }

  for(int i=0; i < ncon; ++i) {
    array[i].key = flows[i] * sign;
    array[i].val = i;
  }

  // May be need to check the direction of the sort
  rkvsorti(ncon, array.data());

  for(int current=0, i=0; i < ncon-1; ++i) {
    if( (array[i+1].key - array[i].key < maxdiff * BALANCE_THRESHOLD) &&
	(not_needed[current] < ncon-1) ) {
      not_needed[current]++;
      not_needed[i+1] = 0;
    } else {
      current = i+1;
    }
  }

  int nperms, perm[24][5];
  switch (ncon) {
    // ----------------------
    case 2:
      nperms = 1;
      perm[0][0] = 0;   perm[0][1] = 1;
      break;

    // ----------------------
    case 3:
      // if the first and second flows are close
      if( not_needed[0] == 2 && not_needed[1] == 0 && not_needed[2] == 1 ) {
        nperms = 4;
        perm[0][0] = 0;   perm[0][1] = 1;   perm[0][2] = 2;
        perm[1][0] = 1;   perm[1][1] = 0;   perm[1][2] = 2;
        perm[2][0] = 0;   perm[2][1] = 2;   perm[2][2] = 1;
        perm[3][0] = 1;   perm[3][1] = 2;   perm[3][2] = 0;
        break;
      }

      // if the second and third flows are close
      if( not_needed[0] == 1 && not_needed[1] == 2 && not_needed[2] == 0 ) {
        nperms = 4;
        perm[0][0] = 0;   perm[0][1] = 1;   perm[0][2] = 2;
        perm[1][0] = 0;   perm[1][1] = 2;   perm[1][2] = 1;
        perm[2][0] = 1;   perm[2][1] = 0;   perm[2][2] = 2;
        perm[3][0] = 2;   perm[3][1] = 0;   perm[3][2] = 1;
        break;
      }

      // all or none of the flows are close
      nperms = 3;
      perm[0][0] = 0;   perm[0][1] = 1;   perm[0][2] = 2;
      perm[1][0] = 1;   perm[1][1] = 0;   perm[1][2] = 2;
      perm[2][0] = 0;   perm[2][1] = 2;   perm[2][2] = 1;

      break;

    // ----------------------
    case 4:
      if( not_needed[0] == 2 && not_needed[1] == 0 &&
          not_needed[2] == 1 && not_needed[3] == 1 ) {
        nperms = 14;
        perm[0][0] =  0;   perm[0][1] =  1;   perm[0][2] =  2;   perm[0][3] =  3;
        perm[1][0] =  1;   perm[1][1] =  0;   perm[1][2] =  2;   perm[1][3] =  3;
        perm[2][0] =  0;   perm[2][1] =  2;   perm[2][2] =  1;   perm[2][3] =  3;
        perm[3][0] =  1;   perm[3][1] =  2;   perm[3][2] =  0;   perm[3][3] =  3;
        perm[4][0] =  0;   perm[4][1] =  1;   perm[4][2] =  3;   perm[4][3] =  2;
        perm[5][0] =  1;   perm[5][1] =  0;   perm[5][2] =  3;   perm[5][3] =  2;
        
        perm[6][0] =  0;   perm[6][1] =  3;   perm[6][2] =  1;   perm[6][3] =  2;
        perm[7][0] =  1;   perm[7][1] =  3;   perm[7][2] =  0;   perm[7][3] =  2;

        perm[8][0] =  0;   perm[8][1] =  2;   perm[8][2] =  3;   perm[8][3] =  1;
        perm[9][0] =  1;   perm[9][1] =  2;   perm[9][2] =  3;   perm[9][3] =  0;

        perm[10][0] = 2;   perm[10][1] = 0;   perm[10][2] = 1;   perm[10][3] = 3;
        perm[11][0] = 2;   perm[11][1] = 1;   perm[11][2] = 0;   perm[11][3] = 3;
        
        perm[12][0] = 0;   perm[12][1] = 3;   perm[12][2] = 2;   perm[12][3] = 1;
        perm[13][0] = 1;   perm[13][1] = 3;   perm[13][2] = 2;   perm[13][3] = 0;
        break;
      }

      if( not_needed[0] == 1 && not_needed[1] == 1 &&
          not_needed[2] == 2 && not_needed[3] == 0 ) {
        nperms = 14;
        perm[0][0] =  0;   perm[0][1] =  1;   perm[0][2] =  2;   perm[0][3] =  3;
        perm[1][0] =  0;   perm[1][1] =  1;   perm[1][2] =  3;   perm[1][3] =  2;
        perm[2][0] =  0;   perm[2][1] =  2;   perm[2][2] =  1;   perm[2][3] =  3;
        perm[3][0] =  0;   perm[3][1] =  3;   perm[3][2] =  1;   perm[3][3] =  2;
        perm[4][0] =  1;   perm[4][1] =  0;   perm[4][2] =  2;   perm[4][3] =  3;
        perm[5][0] =  1;   perm[5][1] =  0;   perm[5][2] =  3;   perm[5][3] =  2;

        perm[6][0] =  1;   perm[6][1] =  2;   perm[6][2] =  0;   perm[6][3] =  3;
        perm[7][0] =  1;   perm[7][1] =  3;   perm[7][2] =  0;   perm[7][3] =  2;

        perm[8][0] =  2;   perm[8][1] =  0;   perm[8][2] =  1;   perm[8][3] =  3;
        perm[9][0] =  3;   perm[9][1] =  0;   perm[9][2] =  1;   perm[9][3] =  2;

        perm[10][0] = 0;   perm[10][1] = 2;   perm[10][2] = 3;   perm[10][3] = 1;
        perm[11][0] = 0;   perm[11][1] = 3;   perm[11][2] = 2;   perm[11][3] = 1;

        perm[12][0] = 2;   perm[12][1] = 1;   perm[12][2] = 0;   perm[12][3] = 3;
        perm[13][0] = 3;   perm[13][1] = 1;   perm[13][2] = 0;   perm[13][3] = 2;
        break;
      }

      if( not_needed[0] == 2 && not_needed[1] == 0 &&
          not_needed[2] == 2 && not_needed[3] == 0 ) {
        nperms = 14;
        perm[0][0] =  0;   perm[0][1] =  1;   perm[0][2] =  2;   perm[0][3] =  3;
        perm[1][0] =  1;   perm[1][1] =  0;   perm[1][2] =  2;   perm[1][3] =  3;
        perm[2][0] =  0;   perm[2][1] =  1;   perm[2][2] =  3;   perm[2][3] =  2;
        perm[3][0] =  1;   perm[3][1] =  0;   perm[3][2] =  3;   perm[3][3] =  2;

        perm[4][0] =  0;   perm[4][1] =  2;   perm[4][2] =  1;   perm[4][3] =  3;
        perm[5][0] =  1;   perm[5][1] =  2;   perm[5][2] =  0;   perm[5][3] =  3;
        perm[6][0] =  0;   perm[6][1] =  3;   perm[6][2] =  1;   perm[6][3] =  2;
        perm[7][0] =  1;   perm[7][1] =  3;   perm[7][2] =  0;   perm[7][3] =  2;

        perm[8][0] = 2;    perm[8][1] = 0;    perm[8][2] = 1;    perm[8][3] = 3;
        perm[9][0] = 0;    perm[9][1] = 2;    perm[9][2] = 3;    perm[9][3] = 1;
        perm[10][0] = 2;   perm[10][1] = 1;   perm[10][2] = 0;   perm[10][3] = 3;
        perm[11][0] = 0;   perm[11][1] = 3;   perm[11][2] = 2;   perm[11][3] = 1;
        perm[12][0] = 3;   perm[12][1] = 0;   perm[12][2] = 1;   perm[12][3] = 2;
        perm[13][0] = 1;   perm[13][1] = 2;   perm[13][2] = 3;   perm[13][3] = 0;
        break;
      }

      if( not_needed[0] == 3 && not_needed[1] == 0 &&
          not_needed[2] == 0 && not_needed[3] == 1 ) {
        nperms = 14;
        perm[0][0] =  0;   perm[0][1] =  1;   perm[0][2] =  2;   perm[0][3] =  3;
        perm[1][0] =  0;   perm[1][1] =  2;   perm[1][2] =  1;   perm[1][3] =  3;
        perm[2][0] =  1;   perm[2][1] =  0;   perm[2][2] =  2;   perm[2][3] =  3;
        perm[3][0] =  2;   perm[3][1] =  0;   perm[3][2] =  1;   perm[3][3] =  3;
        perm[4][0] =  1;   perm[4][1] =  2;   perm[4][2] =  0;   perm[4][3] =  3;
        perm[5][0] =  2;   perm[5][1] =  1;   perm[5][2] =  0;   perm[5][3] =  3;

        perm[6][0] =  0;   perm[6][1] =  1;   perm[6][2] =  3;   perm[6][3] =  2;
        perm[7][0] =  1;   perm[7][1] =  0;   perm[7][2] =  3;   perm[7][3] =  2;
        perm[8][0] =  0;   perm[8][1] =  2;   perm[8][2] =  3;   perm[8][3] =  1;
        perm[9][0] =  2;   perm[9][1] =  0;   perm[9][2] =  3;   perm[9][3] =  1;
        perm[10][0] = 1;   perm[10][1] = 2;   perm[10][2] = 3;   perm[10][3] = 0;
        perm[11][0] = 2;   perm[11][1] = 1;   perm[11][2] = 3;   perm[11][3] = 0;

        perm[12][0] = 0;   perm[12][1] = 3;   perm[12][2] = 1;   perm[12][3] = 2;
        perm[13][0] = 0;   perm[13][1] = 3;   perm[13][2] = 2;   perm[13][3] = 1;
        break;
      }

      if( not_needed[0] == 1 && not_needed[1] == 3 &&
          not_needed[2] == 0 && not_needed[3] == 0 ) {
        nperms = 14;
        perm[0][0] =  0;   perm[0][1] =  1;   perm[0][2] =  2;   perm[0][3] =  3;
        perm[1][0] =  0;   perm[1][1] =  2;   perm[1][2] =  1;   perm[1][3] =  3;
        perm[2][0] =  0;   perm[2][1] =  1;   perm[2][2] =  3;   perm[2][3] =  2;
        perm[3][0] =  0;   perm[3][1] =  2;   perm[3][2] =  3;   perm[3][3] =  1;
        perm[4][0] =  0;   perm[4][1] =  3;   perm[4][2] =  1;   perm[4][3] =  2;
        perm[5][0] =  0;   perm[5][1] =  3;   perm[5][2] =  2;   perm[5][3] =  1;

        perm[6][0] =  1;   perm[6][1] =  0;   perm[6][2] =  2;   perm[6][3] =  3;
        perm[7][0] =  1;   perm[7][1] =  0;   perm[7][2] =  3;   perm[7][3] =  2;
        perm[8][0] =  2;   perm[8][1] =  0;   perm[8][2] =  1;   perm[8][3] =  3;
        perm[9][0] =  2;   perm[9][1] =  0;   perm[9][2] =  3;   perm[9][3] =  1;
        perm[10][0] = 3;   perm[10][1] = 0;   perm[10][2] = 1;   perm[10][3] = 2;
        perm[11][0] = 3;   perm[11][1] = 0;   perm[11][2] = 2;   perm[11][3] = 1;

        perm[12][0] = 1;   perm[12][1] = 2;   perm[12][2] = 0;   perm[12][3] = 3;
        perm[13][0] = 2;   perm[13][1] = 1;   perm[13][2] = 0;   perm[13][3] = 3;

        break;
      }

      nperms = 14;
      perm[0][0] =  0;   perm[0][1] =  1;   perm[0][2] =  2;   perm[0][3] =  3;
      perm[1][0] =  1;   perm[1][1] =  0;   perm[1][2] =  2;   perm[1][3] =  3;
      perm[2][0] =  0;   perm[2][1] =  2;   perm[2][2] =  1;   perm[2][3] =  3;
      perm[3][0] =  0;   perm[3][1] =  1;   perm[3][2] =  3;   perm[3][3] =  2;
      perm[4][0] =  1;   perm[4][1] =  0;   perm[4][2] =  3;   perm[4][3] =  2;

      perm[5][0] =  2;   perm[5][1] =  0;   perm[5][2] =  1;   perm[5][3] =  3;
      perm[6][0] =  0;   perm[6][1] =  2;   perm[6][2] =  3;   perm[6][3] =  1;

      perm[7][0] =  1;   perm[7][1] =  2;   perm[7][2] =  0;   perm[7][3] =  3;
      perm[8][0] =  0;   perm[8][1] =  3;   perm[8][2] =  1;   perm[8][3] =  2;

      perm[9][0] =  2;   perm[9][1] =  1;   perm[9][2] =  0;   perm[9][3] =  3;
      perm[10][0] = 0;   perm[10][1] = 3;   perm[10][2] = 2;   perm[10][3] = 1;
      perm[11][0] = 2;   perm[11][1] = 0;   perm[11][2] = 3;   perm[11][3] = 1;

      perm[12][0] = 3;   perm[12][1] = 0;   perm[12][2] = 1;   perm[12][3] = 2;
      perm[13][0] = 1;   perm[13][1] = 2;   perm[13][2] = 3;   perm[13][3] = 0;
      break;

    // ----------------------
    default:
      return;
  }

  {
    std::vector<int> candidate(ncon);
    std::vector<int> rank     (ncon);
    for(int i=0; i < nperms; ++i) {
      for(int j=0; j < ncon; ++j)
	candidate[j] = array[perm[i][j]].val;

      for(int j=0; j < ncon; ++j)
	rank[candidate[j]] = j;

      int hash = get_rank_hash(ncon, rank.data()) - minval;
      if( currentq[hash+index] > 0 ) {
	*qnum = hash;
	return;
      }
    }
  }
}


//------------------------------------------------------------------------
// Extract a subgraph from a graph given an indicator array
Graph *extract_graph(Control *ctrl, Graph *graph, int *indicator,
		     int *map, int *rmap)
{
  int count = 0;
  for(int i=0; i < graph->nvtxs; ++i) {
    if( indicator[i] == 1 ) {
      map[count] = i;
      rmap[i] = count;
      ++count;
    }
  }

  if( count == 0 ) 
    return NULL;

  //------------------------------------------------
  // Allocate memory
  Graph *egraph = create_graph();
  egraph->ncon  = graph->ncon;
  egraph->nvtxs = count;
  int envtxs = egraph->nvtxs;

  egraph->xadj.resize(envtxs+1);
  int *exadj = egraph->xadj.data();

  egraph->where.resize(envtxs);
  egraph->vsize.resize(envtxs);

  egraph->nvwgt.resize(envtxs*graph->ncon);

  //------------------------------------------------
  // Compute xadj, where, nvwgt, and vsize arrays
  iset(envtxs+1, 0, exadj);
  for(int i=0; i < envtxs; ++i) {
    int im = map[i];

    egraph->where[i] = graph->where[im];
    for(int j=0; j < graph->ncon; ++j)
      egraph->nvwgt[i*graph->ncon+j] = graph->nvwgt[im*graph->ncon+j];

    if( ctrl->is_adaptive )
      egraph->vsize[i] = graph->vsize[im];

    for(int j=graph->xadj[im]; j < graph->xadj[im+1]; ++j)
      if( indicator[graph->adjncy[j]] == 1 )
        exadj[i]++;
  }

  for(int i=1; i < envtxs; i++) exadj[i] += exadj[i-1];
  for(int i=envtxs; i > 0; i--) exadj[i] = exadj[i-1];
  exadj[0] = 0;

  //------------------------------------------------
  // Compute adjncy and adjwgt arrays
  egraph->nedges = exadj[envtxs];
  egraph->adjncy.resize(egraph->nedges);
  egraph->adjwgt.resize(egraph->nedges);

  for(int i=0; i < envtxs; ++i) {
    int im = map[i];
    for(int j=graph->xadj[im]; j < graph->xadj[im+1]; ++j) {
      if( indicator[graph->adjncy[j]] == 1 ) {
        egraph->adjncy[exadj[i]]   = rmap[graph->adjncy[j]];
        egraph->adjwgt[exadj[i]++] = graph->adjwgt[j];
      }
    }
  }

  for(int i=envtxs; i > 0; --i) exadj[i] = exadj[i-1];
  exadj[0] = 0;

  return egraph;
}

//------------------------------------------------------------------------
// Find a matching using the HEM heuristic
void find_match(matrix_t *matrix, int *match, int *mlist,
		int *skip, int ncon)
{
  int     nrows    = matrix->nrows;    // convinience shortcut
  int    *rowptr   = matrix->rowptr;   // convinience shortcut
  int    *colind   = matrix->colind;   // convinience shortcut
  double *transfer = matrix->transfer; // convinience shortcut

  iset(nrows, UNMATCHED, match);

  std::vector<rkv_t> links(nrows);

  for(int i=0; i < nrows; ++i) {
    links[i].key = 0.0;
    links[i].val = i;
    for(int j=rowptr[i]; j < rowptr[i+1]; ++j) {
      for(int k=0;  k < ncon; ++k) {
	double avalue = std::abs(transfer[j*ncon+k]);
        if( links[i].key < avalue )
          links[i].key = avalue;
      }
    }
  }

  rkvsortd(nrows, links.data());

  for(int count=0, ii=0; ii < nrows; ++ii) {
    int i = (int)links[ii].val;

    if( match[i] == UNMATCHED ) {
      int  maxidx = i;
      double maxwgt = 0.0;

      // Find a heavy-edge matching
      for(int j=rowptr[i]; j < rowptr[i+1]; ++j) {
        int edge = colind[j];
        if( match[edge] == UNMATCHED && edge != i && skip[j] == 0 ) {
	  bool large_found = false;
	  double avalue;
          for(int k=0; k < ncon; ++k) {
	    avalue = std::abs(transfer[j*ncon+k]);
            if( maxwgt < avalue ) {
	      large_found = true;
              break;
	    }
	  }

          if( large_found ) {
            maxwgt = avalue;
            maxidx = edge;
          }
        }
      }

      if( maxidx != i ) {
        match[i] = maxidx;
        match[maxidx] = i;
        mlist[count++] = std::max(i, maxidx);
        mlist[count++] = std::min(i, maxidx);
      }
    }
  }
}


//------------------------------------------------------------------------
// Sort nvwgts of a vertex and return a hashed value
int get_weights_hash(int ncon, const double *nvwgt)
{
  std::vector<rkv_t> array(ncon);
  for(int i=0; i < ncon; ++i) {
    array[i].key = nvwgt[i];
    array[i].val = i;
  }

  rkvsorti(ncon, array.data());

  std::vector<int> rank(ncon);
  for(int i=0; i < ncon; ++i)
    rank[array[i].val] = i;

  int multiplier = 1;
  int retval = 0;
  for(int i=0; i < ncon; ++i) {
    multiplier *= (i+1);
    retval += rank[ncon-i-1] * multiplier;
  }

  return retval;
}

//------------------------------------------------------------------------
// Perform an edge-based FM refinement
int balance_links(Control *ctrl, Graph *graph, int *home, int me,
          int you, double *flows, double maxdiff, double *diff_cost, 
          double *diff_lbavg, double avgvwgt)
{
  static const int BALANCE_STEPS = 4;

  int     nvtxs = graph->nvtxs; // convinience shortcut
  int     ncon  = graph->ncon;  // convinience shortcut
  int    *where = graph->where.data(); // convinience shortcut
  double *nvwgt = graph->nvwgt.data(); // convinience shortcut

  std::vector<double> target_weights(2*ncon);
  for(int i=0; i < ncon; ++i) {
    target_weights[i]      = -flows[i];
    target_weights[ncon+i] =  flows[i];
  }

  std::vector<double> pwgts(2*ncon, 0.);
  for(int i=0; i < nvtxs; ++i) {
    if( where[i] == me ) {
      for(int j=0; j < ncon; ++j) {
        target_weights[j] += nvwgt[i*ncon+j];
        pwgts[j] += nvwgt[i*ncon+j];
      }
    } else {
      assert(where[i] == you);
      for(int j=0; j < ncon; ++j) {
        target_weights[ncon+j] += nvwgt[i*ncon+j];
        pwgts[ncon+j]          += nvwgt[i*ncon+j];
      }
    }
  }

  // All target_weights should be non-negative
  for(int i=0; i < ncon; ++i) {
    if( target_weights[i] < 0. ) {
      target_weights[ncon+i] += target_weights[i];
      target_weights[i] = 0.;
    }

    if( target_weights[ncon+i] < 0. ) {
      target_weights[i] += target_weights[ncon+i];
      target_weights[ncon+i] = 0.;
    }
  }

  //-------------------------------------------------------------
  // Insert vertices into queues
  int minval = 0, maxval = 0;
  int multiplier = 1;
  for(int i=0; i < ncon; ++i) {
    multiplier *= (i+1);
    maxval += i*multiplier;
    minval += (ncon-1-i)*multiplier;
  }

  std::vector<int> hval(nvtxs);
  for(int i=0; i < nvtxs; ++i)
    hval[i] = get_weights_hash(ncon, nvwgt+i*ncon) - minval;

  std::vector<int> nvpq(nvtxs, 0);
  for(int i=0; i < nvtxs; ++i)
    ++nvpq[hval[i]];

  int nqueues = maxval - minval + 1;
  std::vector<int> ptr(nqueues+1);
  ptr[0] = 0;
  for(int i=0; i < nqueues; ++i)
    ptr[i+1] = ptr[i] + nvpq[i];

  std::vector<int> map (nvtxs);
  std::vector<int> rmap(nvtxs);
  for(int i=0; i < nvtxs; ++i) {
    map[i] = ptr[hval[i]];
    rmap[ptr[hval[i]]++] = i;
  }

  for(int i=nqueues; i > 0; --i) ptr[i] = ptr[i-1];
  ptr[0] = 0;

  // Initialize queues
  std::vector<rpq_t *> queues(2*nqueues);
  for(int i=0; i < nqueues; ++i)
    if( nvpq[i] > 0 ) {
      queues[i]         = rpqCreate(nvpq[i]);
      queues[nqueues+i] = rpqCreate(nvpq[i]);
    }

  // Compute internal/external degrees
  std::vector<int> id(nvtxs, 0);
  std::vector<int> ed(nvtxs, 0);
  for(int i=0; i<nvtxs; ++i) {
    for(int j=graph->xadj[i]; j<graph->xadj[i+1]; ++j) {
      if( where[graph->adjncy[j]] == where[i] )
        id[i] += graph->adjwgt[j];
      else 
        ed[i] += graph->adjwgt[j];
    }
  }

  std::vector<int> myqueue(nvtxs);
  std::vector<int> changes(nvtxs);
  std::vector<int> inq(2*nqueues);
  int nswaps = 0;
  for(int pass=0; pass < BALANCE_STEPS; ++pass) {
    iset(nvtxs, -1, myqueue.data()); 
    iset(nqueues*2, 0, inq.data());

    // Insert vertices into correct queues
    for(int j=0; j < nvtxs; ++j) {
      int index = (where[j] == me) ? 0 : nqueues;

      double newgain = ctrl->itr_ratio*(double)(ed[j]-id[j]);
      if( home[j] == me || home[j] == you ) {
        if( where[j] == home[j] )
          newgain -= ctrl->redist_ratio * graph->vsize[j];
        else
          newgain += ctrl->redist_ratio * graph->vsize[j];
      }

      rpqInsert(queues[hval[j]+index], map[j] - ptr[hval[j]], newgain);
      myqueue[j] = (where[j] == me) ? 0 : 1;
      inq[hval[j]+index]++;
    }

    int max_indx = 0;
    for(int i=0; i < ncon; ++i) {
      if( std::abs(flows[i]) > std::abs(flows[max_indx]) )
        max_indx = i;
    }
    double bestflow = std::abs(flows[max_indx]);

    int qnum, nchanges = 0, nmoves = 0;
    for(int ii=0; ii < nvtxs/2; ++ii) {
      int from = -1;
      select_queue_dynamic(nqueues, ncon, me, you, inq.data(), flows, 
                           &from, &qnum, minval, avgvwgt, maxdiff);

      // Vertex is not found in one subdomain, try the other
      if( from != -1 && qnum == -1 ) {
        from = (from == me) ? you : me;

	bool found = false;
        if( from == me ) {
          for(int j=0; j < ncon; ++j) {
            if( flows[j] > avgvwgt ) {
	      found = true;
              break;
	    }
          }
        } else {
          for(int j=0; j < ncon; ++j) {
            if( flows[j] < -avgvwgt ) {
	      found = true;
              break;
	    }
          }
        }

        if( found )
          select_queue_dynamic(nqueues, ncon, me, you, inq.data(), flows, 
                               &from, &qnum, minval, avgvwgt, maxdiff);
      }

      if( qnum == -1 )
        break;

      int to = (from == me) ? you : me;
      int index = (from == me) ? 0 : nqueues;
      int higain = rpqGetTop(queues[qnum+index]);
      --inq[qnum+index];
      assert(higain != -1);

      //-------------------------------------------------------------
      // Make the swap
      int vtx = rmap[higain+ptr[qnum]];
      myqueue[vtx] = -1;
      where[vtx] = to;
      ++nswaps;
      ++nmoves;

      // Update the flows
      for(int j=0; j < ncon; ++j)
        flows[j] += (to == me) ? nvwgt[vtx*ncon+j] : -nvwgt[vtx*ncon+j];

      int max_indx = 0;
      for(int j=0; j < ncon; ++j) {
        if(std::abs(flows[j]) > std::abs(flows[max_indx]) ) 
          max_indx = j;
      }
      double ftmp = std::abs(flows[max_indx]);

      if( ftmp < bestflow ) {
        bestflow = ftmp;
        nchanges = 0;
      } else {
        changes[nchanges++] = vtx;
      }

      std::swap(id[vtx], ed[vtx]);

      for(int j=graph->xadj[vtx]; j < graph->xadj[vtx+1]; ++j) {
        int edge = graph->adjncy[j];

        int tmp = (to == where[edge] ? graph->adjwgt[j] : -graph->adjwgt[j]);
	id[edge] += tmp;
	ed[edge] -= tmp;

        if( myqueue[edge] != -1 ) {
          double newgain = ctrl->itr_ratio * (ed[edge]-id[edge]);
          if( home[edge] == me || home[edge] == you ) {
            if( where[edge] == home[edge] )
              newgain -= ctrl->redist_ratio * graph->vsize[edge];
            else
              newgain += ctrl->redist_ratio * graph->vsize[edge];
          }

          rpqUpdate(queues[hval[edge] + (nqueues*myqueue[edge])], 
		    map[edge] - ptr[hval[edge]], newgain);
        }
      }
    }

    //-------------------------------------------------------------
    // Get back to the best flow
    nswaps -= nchanges;
    nmoves -= nchanges;
    for(int i=0; i < nchanges; ++i) {
      int vtx = changes[i];
      int from = where[vtx];
      int to   = (from == me) ? you : me;
      where[vtx] = to;

      std::swap(id[vtx], ed[vtx]);
      for(int j=graph->xadj[vtx]; j < graph->xadj[vtx+1]; ++j) {
        int edge = graph->adjncy[j];
        int tmp = (to == where[edge] ? graph->adjwgt[j] : -graph->adjwgt[j]);
	id[edge] += tmp;
	ed[edge] -= tmp;
      }
    }

    for(int i=0; i < nqueues; ++i) {
      if( nvpq[i] > 0 ) {
        rpqReset(queues[i]);
        rpqReset(queues[i+nqueues]);
      }
    }

    if( nmoves == 0 )
      break;
  }

  //-------------------------------------------------------------
  // Compute imbalance
  std::vector<double> local_weights(ncon, 0.);
  for(int i=0; i < nvtxs; ++i) {
    if( where[i] == me ) {
      for(int j=0; j < ncon; ++j)
        local_weights[j] += nvwgt[i*ncon+j];
    }
  }

  std::vector<double> lbvec(ncon);
  for(int i=0; i < ncon; ++i) {
    double ftmp = (pwgts[i] + pwgts[ncon+i]) / 2.0;
    if( std::abs(ftmp) > 1.e-12 )
      lbvec[i] = std::abs(local_weights[i] - target_weights[i]) / ftmp;
    else
      lbvec[i] = 0.0;
  }
  double lbavg = average(ncon, lbvec.data());
  *diff_lbavg = lbavg;

  //-------------------------------------------------------------
  // Compute the cost
  int cut = 0, totalv = 0;
  for(int i=0; i < nvtxs; ++i) {
    if( where[i] != home[i] )
      totalv += graph->vsize[i];

      for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
        if( where[graph->adjncy[j]] != where[i] )
          cut += graph->adjwgt[j];
      }
  }
  cut /= 2;

  double mycost = cut*ctrl->itr_ratio + totalv*ctrl->redist_ratio;
  *diff_cost = mycost;

  // Free memory
  for(int i=0; i < nqueues; ++i) {
    if( nvpq[i] > 0 ) {
      rpqDestroy(queues[i]);
      rpqDestroy(queues[i+nqueues]);
    }
  }

  return nswaps;
}


//------------------------------------------------------------------------
// Create matrix form the graph
void setup_matrix(const Graph *graph, matrix_t *matrix, int *workspace)
{
  int     nvtxs  = graph->nvtxs;

  int     nrows  = matrix->nrows;
  int    *rowptr = matrix->rowptr;
  int    *colind = matrix->colind;
  double *values = matrix->values;

  int    *perm   = workspace;
  int    *marker = workspace + nvtxs;
  int    *pcounts= workspace + nvtxs + nrows;

  iset(nrows,  -1, marker);
  iset(nrows+1, 0, pcounts);

  for(int i=0; i < nvtxs; ++i)
    pcounts[graph->where[i]]++;

  for(int i=1; i < nrows; ++i) pcounts[i] += pcounts[i-1];
  for(int i=nrows; i > 0; --i) pcounts[i]  = pcounts[i-1];
  pcounts[0] = 0;

  for(int i=0; i < nvtxs; ++i)
    perm[pcounts[graph->where[i]]++] = i;

  for(int i=nrows; i > 0; --i)
    pcounts[i] = pcounts[i-1];
  pcounts[0] = 0;

  // Construct the matrix
  int k = 0;
  rowptr[0] = 0;
  for(int ii=0; ii < nrows; ++ii) {
    colind[k++] = ii;
    marker[ii]  = ii;

    for(int jj=pcounts[ii]; jj < pcounts[ii+1]; ++jj) {
      int i = perm[jj];
      for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
        int kk = graph->where[graph->adjncy[j]];
        if( marker[kk] != ii ) {
          colind[k]   = kk;
          values[k++] = -1.0;
          marker[kk]  = ii;
        }
      }
    }
    values[rowptr[ii]] = (double)(k - rowptr[ii] - 1);
    rowptr[ii+1] = k;
  }
  matrix->nnzs = rowptr[nrows];
}

//------------------------------------------------------------------------
// Compute the total weight of a serial graph
int compute_totalw(Graph *graph, int *home)
{
  int totalv = 0;

  for(int i=0; i < graph->nvtxs; ++i) {
    if( graph->where[i] != home[i] )
      totalv += (graph->vsize.empty() ? graph->vwgt[i*graph->ncon] : graph->vsize[i]);
  }

  return totalv;
}

//------------------------------------------------------------------------
// Compute the load for each subdomain
void compute_load(Graph *graph, int nparts, double *load,
		  double *target_weights, int index)
{
  int ncon = graph->ncon;

  rset(nparts, 0., load);

  for(int i=0; i < graph->nvtxs; ++i)
    load[graph->where[i]] += graph->nvwgt[i*ncon+index];

  assert(std::abs(rsum(nparts, load, 1)-1.0) < 0.001);

  for(int i=0; i < nparts; ++i) {
    load[i] -= target_weights[i*ncon+index];
  }

  return;
}

//------------------------------------------------------------------------
// Perform matrix-vector multiplication
void mv_mult(matrix_t *A, double *v, double *w)
{
  for(int i=0; i < A->nrows; ++i)
    w[i] = 0.;

  for(int i=0; i < A->nrows; ++i)
    for(int j = A->rowptr[i]; j < A->rowptr[i+1]; ++j)
      w[i] += A->values[j] * v[A->colind[j]];
}

//------------------------------------------------------------------------
// CG solver used in diffusion algorithms
void cg_solve(matrix_t *A, double *b, double *x, double tol, double *workspace)
{
  static const double RMACH = 1.e-16;

  int n = A->nrows;

  // Initial setup
  double *p = workspace;
  double *r = workspace + n;
  double *q = workspace + 2*n;
  double *z = workspace + 3*n;
  double *M = workspace + 4*n;

  for(int i=0; i < n; ++i) {
    x[i] = 0.;
    double dval = A->values[A->rowptr[i]];
    M[i] = ((std::abs(dval) > RMACH) ? 1. / dval : 0.);
  }

  // r = b - Ax
  mv_mult(A, x, r);
  for(int i=0; i < n; ++i)
    r[i] = b[i] - r[i];

  double bnrm2 = rnorm2(n, b, 1);
  if( bnrm2 > 0. ) {
    double error = rnorm2(n, r, 1) / bnrm2;

    double rho_1 = -1.;
    if( error > tol ) {
      // Start iterations
      for(int k=0; k < n; ++k) {
        for(int i=0; i < n; ++i)
          z[i] = r[i] * M[i];

        double rho = rdot(n, r, 1, z, 1);

        if( k == 0 ) {
          rcopy(n, z, p);
        } else {
	  double beta = 0.;
          if( std::abs(rho_1) > RMACH )
            beta = rho / rho_1;
          for(int i=0; i < n; ++i)
            p[i] = z[i] + beta*p[i];
        }

        mv_mult(A, p, q); // q = A*p

        double tmp = rdot(n, p, 1, q, 1);
	double alpha;
        if( tmp != 0.0 )
          alpha = rho / tmp;
        else
          alpha = 0.0;
        raxpy(n,  alpha, p, 1, x, 1);   // x = x + alpha*p
        raxpy(n, -alpha, q, 1, r, 1);   // r = r - alpha*q
        error = rnorm2(n, r, 1) / bnrm2;
        if( error < tol )
          break;

        rho_1 = rho;
      }
    }
  }
}

//------------------------------------------------------------------------
// Set up the transfer vectors
void compute_transfer(int ncon, matrix_t *matrix, double *solution,
		      double *transfer, int index)
{
  for(int j=0; j < matrix->nrows; ++j) {
    for(int k=matrix->rowptr[j]+1; k < matrix->rowptr[j+1]; ++k) {
      if( solution[j] > solution[matrix->colind[k]] )
        transfer[k*ncon+index] = solution[j] - solution[matrix->colind[k]];
      else
        transfer[k*ncon+index] = 0.;
    }
  }
}


//------------------------------------------------------------------------
// Compute the initial id/ed
void get_part_params_serial(Control *ctrl, Graph *graph, int nparts)
{
  int ncon = graph->ncon;  // convinience shortcut

  ASSERT(!graph->ckrinfo.empty());
  memset(graph->ckrinfo.data(), 0, sizeof(ckrinfo_t)*graph->nvtxs);
  ctrl->pool.reset();

  //------------------------------------------------------------
  // Compute id/ed degrees
  //------------------------------------------------------------
  int mincut = 0;
  for(int i=0; i < graph->nvtxs; ++i) {
    int me = graph->where[i];
    rcopy(ncon, &graph->nvwgt[i*ncon], &graph->gnpwgts[me*ncon]);

    ckrinfo_t *myrinfo = &graph->ckrinfo[i];

    for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
      if( me == graph->where[graph->adjncy[j]] ) {
        myrinfo->id += graph->adjwgt[j];
      } else {
        myrinfo->ed += graph->adjwgt[j];
      }
    }

    mincut += myrinfo->ed;

    // Compute the external degrees
    if( myrinfo->ed > 0 ) {
      myrinfo->inbr = ctrl->pool.get_next(graph->xadj[i+1] - graph->xadj[i]);
      std::pair<int,int> *mynbrs = ctrl->pool.data() + myrinfo->inbr;

      for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
        int other = graph->where[graph->adjncy[j]];
        if( me != other ) {
	  int k;
          for(k=0; k < myrinfo->nnbrs; ++k) {
            if( mynbrs[k].first == other ) {
              mynbrs[k].second += graph->adjwgt[j];
              break;
            }
          }
          if( k == myrinfo->nnbrs ) {
            mynbrs[k].first = other;
            mynbrs[k].second= graph->adjwgt[j];
            myrinfo->nnbrs++;
          }
        }
      }
    } else {
      myrinfo->inbr = -1;
    }
  }

  graph->mincut = mincut/2;
}

//------------------------------------------------------------------------
// Compute the initial id/ed degrees
void get_part_params_initial(Graph *graph)
{
  double *npwgts = graph->gnpwgts.data();
  for(int i=0; i < 2*graph->ncon; ++i) {
    npwgts[i] = 0.;
  }

  for(int i=0; i < graph->nvtxs; ++i) {
    graph->sendind[i] = 0;
    graph->recvind[i] = 0;
    graph->sendptr[i] =-1;
  }
  int *id     = graph->sendind.data();
  int *ed     = graph->recvind.data();
  int *bndptr = graph->sendptr.data();
  int *bndind = graph->recvptr.data();

  //------------------------------------------------------------
  int nbnd=0, mincut=0;
  for(int i=0; i < graph->nvtxs; ++i) {
    int me = graph->where[i];
    raxpy(graph->ncon, 1.0, &graph->nvwgt[i*graph->ncon], 1, npwgts+me*graph->ncon, 1);

    for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
      if( me == graph->where[graph->adjncy[j]] )
        id[i] += graph->adjwgt[j];
      else
        ed[i] += graph->adjwgt[j];
    }

    if( ed[i] > 0 || graph->xadj[i] == graph->xadj[i+1] ) {
      mincut += ed[i];
      bndind[nbnd] = i;
      bndptr[i] = nbnd++;
    }
  }

  graph->mincut = mincut/2;
  graph->gnvtxs = nbnd;
}



//------------------------------------------------------------------------
// Compute the load imbalance over all the constrains
void get_load_imbalance_kw(int ncon, int nparts,
                           const double *npwgts, double *lbvec)
{
  for(int i=0; i < ncon; ++i) {
    double max = 0.;
    for(int j=0; j < nparts; ++j) {
      if( npwgts[j*ncon+i] > max )
        max = npwgts[j*ncon+i];
    }

    lbvec[i] = max*nparts;
  }
}

//------------------------------------------------------------------------
// Compute the load imbalance over all the constrains
double get_load_imbalance_2w(int ncon, double *npwgts, double *target_weights)
{
  double max = 0.;

  for(int i=0; i < ncon; ++i) {
    double tmp = 0.;
    if( std::abs(target_weights[i]) > 1.e-10 )
      tmp = std::abs(target_weights[i]-npwgts[i]) / target_weights[i];

    if( max < tmp ) max = tmp;
  }
  
  return 1. + max;
}

//------------------------------------------------------------------------
// Check if vertex weights of two vertices are below a given limit
bool are_weights_below(int ncon,
		       double alpha, const double *vwgt1,
		       double beta,  const double *vwgt2, const double *limit)
{
  for(int i=0; i < ncon; ++i)
    if( alpha*vwgt1[i] + beta*vwgt2[i] > limit[i] )
      return false;

  return true;
}

//------------------------------------------------------------------------
// Check if the vertex weights of two vertices are below a given set of values
bool are_any_weights_below(int ncon, const double *vwgt1, const double *limit)
{
  for(int i=0; i < ncon; ++i)
    if( vwgt1[i] < limit[i] )
      return true;

  return false;
}

//------------------------------------------------------------------------
// Perform refinement
void refine_serial_kw(Control *ctrl, Graph *graph, int nparts, 
                      int *home, double *orig_disbalance, int npasses)
{
  int     ncon   = graph->ncon;           // convinience shortcut
  double *npwgts = graph->gnpwgts.data(); // convinience shortcut

  // Setup the weight intervals of the various partitions
  std::vector<double> disbalance(ncon);
  get_load_imbalance_kw(ncon, nparts, npwgts, disbalance.data());
  for(int i=0; i < ncon; ++i)
    disbalance[i] = std::max(disbalance[i], orig_disbalance[i]);

  std::vector<double> minwgt(nparts*ncon);
  std::vector<double> maxwgt(nparts*ncon);
  for(int i=0; i < nparts; ++i) {
    for(int j=0; j < ncon; ++j) {
      maxwgt[i*ncon+j] = disbalance[j] / nparts;
      minwgt[i*ncon+j] = disbalance[j] * nparts;
    }
  }

  std::vector<ikv_t> candidate(graph->nvtxs);
  for(int pass=0; pass < npasses; ++pass) {
    int oldcut = graph->mincut;

    for(int i=0; i < graph->nvtxs; ++i) {
      candidate[i].key = graph->ckrinfo[i].ed-graph->ckrinfo[i].id;
      candidate[i].val = i;
    }
    ikvsortd(graph->nvtxs, candidate.data());

    for(int iii=0; iii < graph->nvtxs; ++iii) {
      int i = candidate[iii].val;

      ckrinfo_t *myrinfo = &graph->ckrinfo[i];

      if( myrinfo->ed >= myrinfo->id ) {
        int from    = graph->where[i];
        int myhome  = home[i];
        double *nvwgt = &graph->nvwgt[i*ncon];

        if( myrinfo->id > 0 &&
            are_weights_below(ncon, 1.0, npwgts+from*ncon, -1.0, nvwgt, &minwgt[from*ncon]) )
          continue;

        std::pair<int,int> *mynbrs = ctrl->pool.data() + myrinfo->inbr;

	int to, k;
        for(k=myrinfo->nnbrs-1; k >= 0; --k) {
          to = mynbrs[k].first;
          if( (mynbrs[k].second - myrinfo->id) >= 0 && 
	      (are_weights_below(ncon, 1.0, npwgts+to*ncon, 1.0, nvwgt, &maxwgt[to*ncon]) ||
	       check_better_balance_wt(ncon, npwgts+from*ncon, npwgts+to*ncon,
                                       nvwgt, disbalance.data())) ) {
            break;
          }
        }

        // break out if you did not find a candidate
        if( k < 0 )
          continue;

        for(int j=k-1; j >= 0; --j) {
          to = mynbrs[j].first;
          bool fit_in_to = are_weights_below(ncon, 1., npwgts+to*ncon, 1., nvwgt, 
					     &maxwgt[to*ncon]);

          if( ((mynbrs[j].second >  mynbrs[k].second) &&
	       (fit_in_to ||
		check_better_balance_wt(ncon, npwgts+from*ncon,
                                        npwgts+to*ncon, nvwgt, disbalance.data())))
            ||
	      ((mynbrs[j].second == mynbrs[k].second) &&
	       ((fit_in_to && (myhome == to)) ||
		check_better_balance_pt(ncon, npwgts+mynbrs[k].first*ncon,
                                        npwgts+to*ncon, nvwgt, disbalance.data())))
	    ) {
            k = j;
          }
        }

        to = mynbrs[k].first;

        bool fit_in_from    = are_weights_below(ncon, 1.0, npwgts+from*ncon, 0.0, 
						npwgts+from*ncon, &maxwgt[from*ncon]);
        bool better_balance = check_better_balance_wt(ncon,
                                                      npwgts+from*ncon, npwgts+to*ncon,
                                                      nvwgt, disbalance.data());

        if( (mynbrs[k].second == myrinfo->id) && (myhome != to) &&
	    !better_balance && fit_in_from)
          continue;

        //---------------------------------------------------------------------
        // At that stage the vertex can be moved from 'from' to 'to' 
        graph->mincut -= mynbrs[k].second - myrinfo->id;

        // Update where, weight, and ID/ED information of the vertex you moved
        raxpy(ncon,  1.0, nvwgt, 1, npwgts+to*ncon,   1);
        raxpy(ncon, -1.0, nvwgt, 1, npwgts+from*ncon, 1);
        graph->where[i] = to;
        myrinfo->ed += myrinfo->id-mynbrs[k].second;
        std::swap(myrinfo->id, mynbrs[k].second);

        if( mynbrs[k].second == 0 )
          mynbrs[k] = mynbrs[--myrinfo->nnbrs];
        else
          mynbrs[k].first = from;

        // Update the degrees of adjacent vertices
        for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
          int ii = graph->adjncy[j];
          int me = graph->where[ii];

          myrinfo = &graph->ckrinfo[ii];
          if( myrinfo->inbr == -1 ) {
            myrinfo->inbr  = ctrl->pool.get_next(graph->xadj[ii+1]-graph->xadj[ii]);
            myrinfo->nnbrs = 0;
          }
          mynbrs = ctrl->pool.data() + myrinfo->inbr;

          if( me == from ) {
	    myrinfo->ed += graph->adjwgt[j];
	    myrinfo->id -= graph->adjwgt[j];
          } else if(me == to ) {
	    myrinfo->id += graph->adjwgt[j];
	    myrinfo->ed -= graph->adjwgt[j];
          }

          // Remove contribution of the ed from 'from'
          if( me != from ) {
            for(int k=0; k < myrinfo->nnbrs; ++k) {
              if( mynbrs[k].first == from ) {
                if( mynbrs[k].second == graph->adjwgt[j] )
                  mynbrs[k] = mynbrs[--myrinfo->nnbrs];
                else
                  mynbrs[k].second -= graph->adjwgt[j];
                break;
              }
            }
          }

          // Add contribution of the ed to 'to'
          if( me != to ) {
            for(int k=0; k < myrinfo->nnbrs; ++k) {
              if( mynbrs[k].first == to ) {
                mynbrs[k].second += graph->adjwgt[j];
                break;
              }
            }
            if( k == myrinfo->nnbrs ) {
              mynbrs[k].first = to;
              mynbrs[k].second= graph->adjwgt[j];
              myrinfo->nnbrs++;
            }
          }
        }
      }
    }

    if( graph->mincut == oldcut )
      break;
  }
}

//------------------------------------------------------------------------
// Select the partition number and the queue from which vertices will be moved
int select_queue_oneway(int ncon,
                          const double *npwgts, const double *target_weights, 
                          int from, rpq_t **queues[2])
{
  int cnum = -1;
  double max = 0.;

  for(int i=0; i < ncon; ++i) {
    if( npwgts[from*ncon+i]-target_weights[from*ncon+i] >= max &&
        rpqLength(queues[0][i]) + rpqLength(queues[1][i]) > 0 ) {
      max = npwgts[from*ncon+i] - target_weights[i];
      cnum = i;
    }
  }

  return cnum;
}


//------------------------------------------------------------------------
// Balance two partitions by moving nodes from the overweight part to
// the underweight part
void balance_fm_init(Graph *graph, double *target_weights)
{
  int     ncon   = graph->ncon;           // convinience shortcut
  int    *id     = graph->sendind.data(); // convinience shortcut
  int    *ed     = graph->recvind.data(); // convinience shortcut
  int    *bndptr = graph->sendptr.data(); // convinience shortcut
  int    *bndind = graph->recvptr.data(); // convinience shortcut
  double *nvwgt  = graph->nvwgt.data();   // convinience shortcut
  double *npwgts = graph->gnpwgts.data(); // convinience shortcut

  std::vector<int> qnum(graph->nvtxs);

  rpq_t **parts[2];
  std::vector<rpq_t *> rpq_space(2*ncon);
  parts[0] = rpq_space.data();
  parts[1] = rpq_space.data() + ncon;

  // This is called for initial partitioning so we know from where to pick nodes
  int from = 1;
  int to = (from+1)%2;

  for(int i=0; i < ncon; ++i) {
    parts[0][i] = rpqCreate(graph->nvtxs);
    parts[1][i] = rpqCreate(graph->nvtxs);
  }

  // Compute the queues in which each vertex will be assigned to
  for(int i=0; i < graph->nvtxs; ++i)
    qnum[i] = max_index(ncon, nvwgt+i*ncon);

  // Insert the nodes of the proper partition in the appropriate priority queue
  for(int i=0; i < graph->nvtxs; ++i) {
    if( graph->where[i] == from ) {
      if( ed[i] > 0 )
        rpqInsert(parts[0][qnum[i]], i, (double)(ed[i]-id[i]));
      else
        rpqInsert(parts[1][qnum[i]], i, (double)(ed[i]-id[i]));
    }
  }

  int mincut = graph->mincut;
  int nbnd   = graph->gnvtxs;
  for(int nswaps=0; nswaps < graph->nvtxs; ++nswaps) {
    if( are_any_weights_below(ncon, npwgts+from*ncon, target_weights+from*ncon) )
      break;

    int cnum = select_queue_oneway(ncon, npwgts, target_weights, from, parts);
    if( cnum == -1 )
      break;

    int higain = rpqGetTop(parts[0][cnum]);
    if( higain == -1 )
      higain = rpqGetTop(parts[1][cnum]);

    mincut -= (ed[higain]-id[higain]);
    raxpy(ncon,  1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon,   1);
    raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);

    graph->where[higain] = to;

    // Update the id[i]/ed[i] values of the affected nodes
    std::swap(id[higain], ed[higain]);
    if( ed[higain] == 0 && bndptr[higain] != -1 && graph->xadj[higain] < graph->xadj[higain+1] )
    {
      bndind[bndptr[higain]] = bndind[--nbnd];
      bndptr[bndind[nbnd]]   = bndptr[higain];
      bndptr[higain] = -1;
    }
    if( ed[higain] >  0 && bndptr[higain] == -1 )
    {
      bndind[nbnd]   = higain;
      bndptr[higain] = nbnd++;
    }

    for(int j=graph->xadj[higain]; j < graph->xadj[higain+1]; ++j) {
      int k = graph->adjncy[j];

      int kwgt = graph->adjwgt[j];
      if( to != graph->where[k] ) kwgt = -kwgt;
      id[k] += kwgt; ed[k] -= kwgt;

      // Update the queue position
      if( graph->where[k] == from ) {
        if( ed[k] > 0 && bndptr[k] == -1 ) {  // It moves in boundary
          rpqDelete(parts[1][qnum[k]], k);
          rpqInsert(parts[0][qnum[k]], k, (double)(ed[k]-id[k]));
        } else { // It must be in the boundary already
          rpqUpdate(parts[0][qnum[k]], k, (double)(ed[k]-id[k]));
        }
      }

      // Update its boundary information
      if( ed[k] == 0 && bndptr[k] != -1 )
      {
        bndind[bndptr[k]]    = bndind[--nbnd];
        bndptr[bndind[nbnd]] = bndptr[k];
        bndptr[k] = -1;
      }
      else if (ed[k] > 0 && bndptr[k] == -1)
      {
        bndind[nbnd] = k;
        bndptr[k]    = nbnd++;
      }
    }
  }

  graph->mincut = mincut;
  graph->gnvtxs = nbnd;

  for(int i=0; i < ncon; ++i) {
    rpqDestroy(parts[0][i]);
    rpqDestroy(parts[1][i]);
  }
}



//------------------------------------------------------------------------
// Select the partition number and the queue from which the vertices will be moved
void select_queue(int ncon, double *npwgts, double *target_weights, int *from, 
                  int *cnum, rpq_t **queues[2])
{
  *from = -1;
  *cnum = -1;

  // First determine the side and the queue, irrespective of the presence of nodes
  double maxdiff=0.0;
  for(int part=0; part < 2; ++part) {
    for(int i=0; i < ncon; ++i) {
      if( npwgts[part*ncon+i] - target_weights[part*ncon+i] >= maxdiff ) {
        maxdiff = npwgts[part*ncon+i] - target_weights[part*ncon+i];
        *from = part;
        *cnum = i;
      }
    }
  }

  if( *from != -1 && rpqLength(queues[*from][*cnum]) == 0 ) {
    // The desired queue is empty, select a node from that side anyway
    double max = -1.0;
    int i;
    for(i=0; i < ncon; ++i) {
      if( rpqLength(queues[*from][i]) > 0 ) {
        max = npwgts[(*from)*ncon + i];
        *cnum = i;
        break;
      }
    }

    for(++i; i<ncon; ++i) {
      if( npwgts[(*from)*ncon + i] > max && rpqLength(queues[*from][i]) > 0 ) {
        max = npwgts[(*from)*ncon + i];
        *cnum = i;
      }
    }
  }

  // Check to see if you can focus on the cut
  if( maxdiff <= 0.0 || *from == -1 ) {
    double maxgain = -100000.0;

    for(int part=0; part < 2; ++part) {
      for(int i=0; i < ncon; ++i) {
        if( rpqLength(queues[part][i]) > 0 &&
            rpqSeeTopKey(queues[part][i]) > maxgain ) {
          maxgain = rpqSeeTopKey(queues[part][i]);
          *from = part;
          *cnum = i;
        }
      }
    }
  }
}

//------------------------------------------------------------------------
// Check if the balance achieved is better than the diff
bool is_balance_better(int ncon, double *npwgts, double *target_weights, 
                       double *diff, double *tmpdiff)
{
  for(int i=0; i < ncon; ++i)
    tmpdiff[i] = std::abs(target_weights[i] - npwgts[i]);

  return (rnorm2(ncon, tmpdiff, 1) < rnorm2(ncon, diff, 1));
}


//------------------------------------------------------------------------
// Perform an edge-based FM refinement
void balance_fm(Graph *graph, double *target_weights, double lbfactor)
{
  int     ncon   = graph->ncon;           // convinience shortcut
  int    *where  = graph->where.data();   // convinience shortcut
  int    *id     = graph->sendind.data(); // convinience shortcut
  int    *ed     = graph->recvind.data(); // convinience shortcut
  int    *bndptr = graph->sendptr.data(); // convinience shortcut
  int    *bndind = graph->recvptr.data(); // convinience shortcut
  double *nvwgt  = graph->nvwgt.data();   // convinience shortcut
  double *npwgts = graph->gnpwgts.data(); // convinience shortcut

  std::vector<double> mindiff_vec(ncon);
  double *mindiff = mindiff_vec.data();
  std::vector<double> tmpdiff_vec(ncon);
  double *tmpdiff = tmpdiff_vec.data();

  rpq_t **parts[2];
  std::vector<rpq_t *> rpq_space(2*ncon);
  parts[0] = rpq_space.data();
  parts[1] = rpq_space.data() + ncon;

  int *qsizes[2];
  std::vector<int> qsize_space(2*ncon, 0);
  qsizes[0] = qsize_space.data();
  qsizes[1] = qsize_space.data() + ncon;

  int limit = (int)std::min(std::max((int)graph->nvtxs/100, 15), 100);

  // Initialize the queues
  for(int i=0; i < ncon; ++i) {
    parts[0][i] = rpqCreate(graph->nvtxs);
    parts[1][i] = rpqCreate(graph->nvtxs);
  }

  std::vector<int> qnum(graph->nvtxs);
  for(int i=0; i < graph->nvtxs; ++i) {
    qnum[i] = max_index(ncon, nvwgt+i*ncon);
    qsizes[qnum[i]][where[i]]++;
  }

  for(int from=0; from < 2; ++from) {
    for(int j=0; j < ncon; ++j) {
      if( qsizes[j][from] == 0 ) {
        for(int i=0; i < graph->nvtxs; ++i) {
          if( where[i] != from )
            continue;

          int k = max_index2(ncon, nvwgt+i*ncon);
          if( k == j &&
	      qsizes[qnum[i]][from] > qsizes[j][from] &&
	      nvwgt[i*ncon+qnum[i]] < 1.3*nvwgt[i*ncon+j] ) {
            qsizes[qnum[i]][from]--;
            qsizes[j][from]++;
            qnum[i] = j;
          }
        }
      }
    }
  }

  for(int i=0; i < ncon; ++i)
    mindiff[i] = std::abs(target_weights[i] - npwgts[i]);
  double origbal = get_load_imbalance_2w(ncon, npwgts, target_weights);
  double minbal  = origbal;
  int  mincut, newcut, mincut_order = -1;
  mincut = newcut = graph->mincut;

  std::vector<int> moved(graph->nvtxs, -1);
  std::vector<int> swaps(graph->nvtxs);

  // Insert all nodes in the priority queues
  int nbnd = graph->gnvtxs;
  for(int i=0; i < graph->nvtxs; ++i) 
    rpqInsert(parts[where[i]][qnum[i]], i, (double)(ed[i]-id[i]));

  int nswaps, from, to, cnum, higain = -1;
  for(nswaps=0; nswaps < graph->nvtxs; ++nswaps) {
    if( minbal < lbfactor )
      break;

    select_queue(ncon, npwgts, target_weights, &from, &cnum, parts);
    to = (from+1)%2;

    if( from == -1 || (higain = rpqGetTop(parts[from][cnum])) == -1 )
      break;

    raxpy(ncon , 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon,   1);
    raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
    newcut -= (ed[higain]-id[higain]);
    double newbal = get_load_imbalance_2w(ncon, npwgts, target_weights);

    if( newbal < minbal ||
	(newbal == minbal &&
	 (newcut < mincut ||
	  (newcut == mincut &&
	   is_balance_better(ncon, npwgts, target_weights, mindiff, tmpdiff)))) ) {
      mincut = newcut;
      minbal = newbal;
      mincut_order = nswaps;
      for(int i=0; i < ncon; ++i)
        mindiff[i] = std::abs(target_weights[i] - npwgts[i]);
    } else if( nswaps - mincut_order > limit ) {
      newcut += (ed[higain]-id[higain]);
      raxpy(ncon,  1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
      raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon,   1);
      break;
    }

    where[higain] = to;
    moved[higain] = nswaps;
    swaps[nswaps] = higain;

    // Update the id[i]/ed[i] values of the affected nodes
    std::swap(id[higain], ed[higain]);
    if( ed[higain] == 0 && bndptr[higain] != -1 && graph->xadj[higain] < graph->xadj[higain+1] ) {
      bndind[bndptr[higain]] = bndind[--nbnd];
      bndptr[bndind[nbnd]]   = bndptr[higain];
      bndptr[higain] = -1;
    }
    if( ed[higain] > 0 && bndptr[higain] == -1 ) {
      bndind[nbnd]   = higain;
      bndptr[higain] = nbnd++;
    }

    for(int j=graph->xadj[higain]; j < graph->xadj[higain+1]; ++j) {
      int k = graph->adjncy[j];

      int kwgt = graph->adjwgt[j];
      if( to != where[k] ) kwgt = -kwgt;
      id[k] += kwgt; ed[k] -= kwgt;

      // Update the queue position
      if( moved[k] == -1 )
        rpqUpdate(parts[where[k]][qnum[k]], k, (double)(ed[k]-id[k]));

      // Update its boundary information
      if( ed[k] == 0 && bndptr[k] != -1 ) {
	bndind[bndptr[k]]    = bndind[--nbnd];
	bndptr[bndind[nbnd]] = bndptr[k];
	bndptr[k] = -1;
      } else if (ed[k] > 0 && bndptr[k] == -1) {
	bndind[nbnd] = k;
	bndptr[k]    = nbnd++;
      }
    }
  }

  // Roll back computations
  for(nswaps--; nswaps > mincut_order; --nswaps) {
    higain = swaps[nswaps];

    to = where[higain] = (where[higain]+1)%2;
    std::swap(id[higain], ed[higain]);
    if( ed[higain] == 0 && bndptr[higain] != -1 && graph->xadj[higain] < graph->xadj[higain+1] ) {
      bndind[bndptr[higain]] = bndind[--nbnd];
      bndptr[bndind[nbnd]]   = bndptr[higain];
      bndptr[higain] = -1;
    } else if( ed[higain] > 0 && bndptr[higain] == -1 ) {
      bndind[nbnd]   = higain;
      bndptr[higain] = nbnd++;
    }

    raxpy(ncon,  1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon,         1);
    raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+((to+1)%2)*ncon, 1);
    for(int j=graph->xadj[higain]; j < graph->xadj[higain+1]; ++j) {
      int k = graph->adjncy[j];

      int kwgt = graph->adjwgt[j];
      if( to != where[k] ) kwgt = -kwgt;
      id[k] += kwgt; ed[k] -= kwgt;

      if( bndptr[k] != -1 && ed[k] == 0 ) {
	bndind[bndptr[k]]    = bndind[--nbnd];
	bndptr[bndind[nbnd]] = bndptr[k];
	bndptr[k] = -1;
      }
      if( bndptr[k] == -1 && ed[k] > 0) {
	bndind[nbnd] = k;
	bndptr[k]    = nbnd++;
      }
    }
  }

  graph->mincut = mincut;
  graph->gnvtxs = nbnd;


  for(int i=0; i < ncon; ++i) {
    rpqDestroy(parts[0][i]);
    rpqDestroy(parts[1][i]);
  }
}


//------------------------------------------------------------------------
// Perform an edge-based FM refinement
void refine_fm(Graph *graph, double *target_weights, int npasses)
{
  int     ncon   = graph->ncon;           // convinience shortcut
  int    *where  = graph->where.data();   // convinience shortcut
  int    *id     = graph->sendind.data(); // convinience shortcut
  int    *ed     = graph->recvind.data(); // convinience shortcut
  int    *bndptr = graph->sendptr.data(); // convinience shortcut
  int    *bndind = graph->recvptr.data(); // convinience shortcut
  double *nvwgt  = graph->nvwgt.data();   // convinience shortcut
  double *npwgts = graph->gnpwgts.data(); // convinience shortcut

  rpq_t **parts[2];
  std::vector<rpq_t *> rpq_space(2*ncon);
  parts[0] = rpq_space.data();
  parts[1] = rpq_space.data() + ncon;

  int limit = (int)std::min(std::max((int)graph->nvtxs/100, 25), 150);

  // Initialize the queues
  for(int i=0; i < ncon; ++i) {
    parts[0][i] = rpqCreate(graph->nvtxs);
    parts[1][i] = rpqCreate(graph->nvtxs);
  }

  std::vector<int> qnum(graph->nvtxs);
  for(int i=0; i < graph->nvtxs; ++i)
    qnum[i] = max_index(ncon, nvwgt+i*ncon);

  double origbal = get_load_imbalance_2w(ncon, npwgts, target_weights);

  std::vector<double> rtpwgts_vec(2*ncon);
  double *rtpwgts = rtpwgts_vec.data();

  for(int i=0; i < ncon; ++i) {
    rtpwgts[i]      = origbal*target_weights[i];
    rtpwgts[ncon+i] = origbal*target_weights[ncon+i];
  }

  std::vector<double> mindiff_vec(ncon);
  double *mindiff = mindiff_vec.data();
  std::vector<double> tmpdiff_vec(ncon);
  double *tmpdiff = tmpdiff_vec.data();

  std::vector<int> swaps( graph->nvtxs);
  std::vector<int> moved(graph->nvtxs, -1);
  for(int pass=0; pass < npasses; ++pass) {
    for(int i=0; i < ncon; ++i) {
      rpqReset(parts[0][i]);
      rpqReset(parts[1][i]);
    }

    int newcut, mincut, initcut, mincut_order = -1;
    newcut = mincut = initcut = graph->mincut;
    for(int i=0; i < ncon; ++i)
      mindiff[i] = std::abs(target_weights[i] - npwgts[i]);
    double minbal = get_load_imbalance_2w(ncon, npwgts, target_weights);

    // Insert boundary nodes in the priority queues
    int nbnd = graph->gnvtxs;
    for(int ii=0; ii < nbnd; ++ii) {
      int i = bndind[ii];
      rpqInsert(parts[where[i]][qnum[i]], i, (double)(ed[i]-id[i]));
    }

    int nswaps, from, to, cnum, higain = -1;
    for(nswaps=0; nswaps < graph->nvtxs; ++nswaps) {
      select_queue(ncon, npwgts, rtpwgts, &from, &cnum, parts);
      to = (from+1)%2;

      if( from == -1 || (higain = rpqGetTop(parts[from][cnum])) == -1 )
        break;

      raxpy(ncon,  1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon,   1);
      raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);

      newcut -= (ed[higain]-id[higain]);
      double newbal = get_load_imbalance_2w(ncon, npwgts, target_weights);

      if( (newcut < mincut && newbal-origbal <= .00001) ||
          (newcut == mincut &&
           (newbal < minbal ||
            (newbal == minbal && 
             is_balance_better(ncon, npwgts, target_weights, mindiff, tmpdiff))))) {
        mincut = newcut;
        minbal = newbal;
        mincut_order = nswaps;
        for(int i=0; i < ncon; ++i)
          mindiff[i] = std::abs(target_weights[i] - npwgts[i]);
      } else if( (nswaps - mincut_order) > limit ) {
        newcut += (ed[higain]-id[higain]);
        raxpy(ncon,  1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
        raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon,   1);
        break;
      }

      where[higain] = to;
      moved[higain] = nswaps;
      swaps[nswaps] = higain;

      // Update the id[i]/ed[i] values of the affected nodes
      std::swap(id[higain], ed[higain]);
      if( ed[higain] == 0 && graph->xadj[higain] < graph->xadj[higain+1] ) {
	bndind[bndptr[higain]] = bndind[--nbnd];
	bndptr[bndind[nbnd]]   = bndptr[higain];
	bndptr[higain] = -1;
      }

      for(int j=graph->xadj[higain]; j < graph->xadj[higain+1]; ++j) {
        int k = graph->adjncy[j];

	int kwgt = graph->adjwgt[j];
	if( to != where[k] ) kwgt = -kwgt;
	id[k] += kwgt; ed[k] -= kwgt;

        // Update its boundary information and queue position
        if( bndptr[k] != -1 ) { // If k was a boundary vertex
          if( ed[k] == 0 ) { // Not a boundary vertex any more
	    bndind[bndptr[k]]    = bndind[--nbnd];
	    bndptr[bndind[nbnd]] = bndptr[k];
	    bndptr[k] = -1;
            if( moved[k] == -1 )  // Remove it if in the queues
              rpqDelete(parts[where[k]][qnum[k]], k);
          } else { // If it has not been moved, update its position in the queue
            if( moved[k] == -1 )
              rpqUpdate(parts[where[k]][qnum[k]], k, (double)(ed[k]-id[k]));
          }
        } else {
          if( ed[k] > 0 ) {  // It will now become a boundary vertex
	    bndind[nbnd] = k;
	    bndptr[k]    = nbnd++;
            if( moved[k] == -1 )
              rpqInsert(parts[where[k]][qnum[k]], k, (double)(ed[k]-id[k]));
          }
        }
      }
    }

    // Roll back computations
    for(int i=0; i < nswaps; ++i)
      moved[swaps[i]] = -1;  // reset moved array
    for(nswaps--; nswaps > mincut_order; --nswaps) {
      higain = swaps[nswaps];

      to = where[higain] = (where[higain]+1)%2;
      std::swap(id[higain], ed[higain]);
      if( ed[higain] == 0 && bndptr[higain] != -1 && graph->xadj[higain] < graph->xadj[higain+1] ) {
	bndind[bndptr[higain]] = bndind[--nbnd];
	bndptr[bndind[nbnd]]   = bndptr[higain];
	bndptr[higain] = -1;
      } else if (ed[higain] > 0 && bndptr[higain] == -1) {
	bndind[nbnd]   = higain;
	bndptr[higain] = nbnd++;
      }

      raxpy(ncon,  1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon,         1);
      raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+((to+1)%2)*ncon, 1);
      for(int j=graph->xadj[higain]; j < graph->xadj[higain+1]; ++j) {
        int k = graph->adjncy[j];

 	int kwgt = graph->adjwgt[j];
	if( to != where[k] ) kwgt = -kwgt;
	id[k] += kwgt; ed[k] -= kwgt;

        if( bndptr[k] != -1 && ed[k] == 0 ) {
	  bndind[bndptr[k]]    = bndind[--nbnd];
	  bndptr[bndind[nbnd]] = bndptr[k];
	  bndptr[k] = -1;
	}
        if( bndptr[k] == -1 && ed[k] > 0 ) {
	  bndind[nbnd] = k;
	  bndptr[k]    = nbnd++;
	}
      }
    }

    graph->mincut = mincut;
    graph->gnvtxs = nbnd;

    if( mincut_order == -1 || mincut == initcut )
      break;
  }

  for(int i=0; i < ncon; ++i) {
    rpqDestroy(parts[0][i]);
    rpqDestroy(parts[1][i]);
  }
}


//------------------------------------------------------------------------
// Comparison function for serial remap

// Triplet structure for double keys
typedef struct i2kv_t {
  int key1, key2;
  int val;
} i2kv_t;

int Comp2keys(const void *fptr, const void *sptr)
{
  i2kv_t *first, *second;

  first  = (i2kv_t *)(fptr);
  second = (i2kv_t *)(sptr);

  if( first->key1 > second->key1 )  return  1;
  if( first->key1 < second->key1 )  return -1;
  if( first->key2 < second->key2 )  return  1;
  if( first->key2 > second->key2 )  return -1;

  return 0;
}

//------------------------------------------------------------------------
// Remap a partitioning on a single processor
void remap_serial(Graph *graph, int nparts, int *base, int *scratch,
                  int *remap, double *target_weights)
{
  std::vector<i2kv_t> sortvtx(graph->nvtxs);
  for(int i=0; i < graph->nvtxs; ++i) {
    sortvtx[i].key1 = base[i];
    sortvtx[i].key2 = graph->vsize[i];
    sortvtx[i].val = i;
  }

  qsort((void *)sortvtx.data(), (size_t)graph->nvtxs, (size_t)sizeof(i2kv_t), Comp2keys);

  int max_mult = std::min(20, nparts); // Uisng not less than 20 parts
  std::vector<ikv_t> bestflow(nparts*max_mult);
  std::vector<ikv_t> flowto(nparts);
  for(int i=0; i < nparts; ++i) {
    flowto[i].key = 0;
    flowto[i].val = i;
  }

  std::vector<int> htable(2*nparts, -1);
  int small_count = 0, big_count = 0, current_from = 0;
  for(int ii=0; ii < graph->nvtxs; ++ii) {
    int i = sortvtx[ii].val;
    int from = base[i];
    int to   = scratch[i];

    if( from > current_from ) {
      // Reset the hash table
      for(int j=0; j < small_count; ++j)
        htable[flowto[j].val] = -1;
      assert(isum(nparts, htable.data(), 1) == -nparts);

      ikvsorti(small_count, flowto.data());

      for(int j=0; j < std::min(small_count, max_mult); ++big_count, ++j) {
        bestflow[big_count].key = flowto[j].key;
        bestflow[big_count].val = flowto[j].val + current_from*nparts;
      }

      small_count = 0;
      current_from = from;
    }

    if( htable[to] == -1 ) {
      htable[to] = small_count;
      flowto[small_count].key = -graph->vsize[i];
      flowto[small_count].val = to;
      ++small_count;
    } else {
      flowto[htable[to]].key += -graph->vsize[i];
    }
  }

  // Reset the hash table
  for(int j=0; j < small_count; ++j)
    htable[flowto[j].val] = -1;
  assert(isum(nparts, htable.data(), 1) == -nparts);

  ikvsorti(small_count, flowto.data());

  for(int j=0; j < std::min(small_count, max_mult); ++big_count, ++j) {
    bestflow[big_count].key = flowto[j].key;
    bestflow[big_count].val = flowto[j].val + current_from*nparts;
  }
  ikvsorti(big_count, bestflow.data());

  int *map    = htable.data();
  int *rowmap = map + nparts;
  assert(isum(nparts, map,    1) == -nparts);
  assert(isum(nparts, rowmap, 1) == -nparts);
  int nmapped = 0;

  for(int ii=0; ii < big_count; ++ii) {
    int i = bestflow[ii].val;
    int j = i % nparts;  // to
    int k = i / nparts;  // from

    if( map[j] == -1 && rowmap[k] == -1 &&
        are_tweights_similar(target_weights, graph->ncon, j, k) ) {
      map[j] = k;
      rowmap[k] = j;
      nmapped++;
    }

    if( nmapped == nparts )
      break;
  }

  // Remap the leftowers
  if( nmapped < nparts ) {
    for(int j=0; j < nparts && nmapped<nparts; ++j) {
      if( map[j] == -1 ) {
        for(int ii=0; ii < nparts; ++ii) {
          int i = (j+ii) % nparts;
          if( rowmap[i] == -1 &&
              are_tweights_similar(target_weights, graph->ncon, i, j) ) {
            map[j] = i;
            rowmap[i] = j;
            nmapped++;
            break;
          }
        }
      }
    }
  }

  // Check if remapping fails; in that case revert to original mapping
  if( nmapped < nparts )
    for(int i=0; i < nparts; ++i)
      map[i] = i;

  for(int i=0; i < graph->nvtxs; ++i)
    remap[i] = map[remap[i]];
}


//------------------------------------------------------------------------
// Perform an edge-based FM refinement
void refine_links(Control *ctrl, Graph *graph, int *home, int me,
		  int you, double *flows, double *sr_cost, double *sr_lbavg)
{
  static const int REFINE_STEPS = 10;

  int  ncon  = graph->ncon;  // convinience shortcut
  int  nvtxs = graph->nvtxs; // convinience shortcut
  int *where = graph->where.data(); // convinience shortcut

  std::vector<int> costwhere(nvtxs);
  std::vector<int> lbwhere(nvtxs);
  std::vector<int> perm(nvtxs);

  std::vector<double> pweights(2*ncon, 0.);
  std::vector<double> target_weights(2*ncon);

  random_permute(nvtxs, perm.data(), 1);
  costwhere.assign(where, where + nvtxs);
  lbwhere  .assign(where, where + nvtxs);

  // Compute target weights
  for(int i=0; i < ncon; ++i) {
    target_weights[i]      = -flows[i];
    target_weights[ncon+i] =  flows[i];
  }

  for(int i=0; i < nvtxs; ++i) {
    if( where[i] == me ) {
      for(int j=0; j < ncon; ++j) {
        target_weights[j] += graph->nvwgt[i*ncon+j];
        pweights[j]       += graph->nvwgt[i*ncon+j];
      }
    } else {
      assert(where[i] == you);
      for(int j=0; j < ncon; ++j) {
        target_weights[ncon+j] += graph->nvwgt[i*ncon+j];
        pweights[ncon+j]       += graph->nvwgt[i*ncon+j];
      }
    }
  }

  // All weights should be positive
  for(int i=0; i < ncon; ++i) {
    if( target_weights[i] < 0. ) {
      target_weights[i+ncon] += target_weights[i];
      target_weights[i] = 0.;
    }

    if( target_weights[ncon+i] < 0. ) {
      target_weights[i] += target_weights[i+ncon];
      target_weights[i+ncon] = 0.;
    }
  } 

  // Compute new bisection
  double bestcost = (isum(graph->nedges, graph->adjwgt.data(), 1)*ctrl->itr_ratio + 
                     isum(nvtxs,         graph->vsize.data(),  1)*ctrl->redist_ratio);
  double best_lbavg = 10.0;

  graph->gnpwgts.resize(2*ncon);

  std::vector<double> lbvec(ncon);
  double other_lbavg = -1., other_cost = -1.;
  int lastseed = 0;
  for(int pass=REFINE_STEPS; pass > 0; --pass) {
    for(int i=0; i < nvtxs; ++i) where[i] = 1;

    // Find seed vertices
    int ii = perm[lastseed] % nvtxs;
    lastseed = (lastseed+1) % nvtxs;
    where[ii] = 0;

    get_part_params_initial(graph);
    balance_fm_init(graph, target_weights.data());
    refine_fm      (graph, target_weights.data(), 4);
    balance_fm     (graph, target_weights.data(), 1.02);
    refine_fm      (graph, target_weights.data(), 4);

    for(int i=0; i < nvtxs; ++i)
      where[i] = (where[i] == 0) ? me : you;

    for(int i=0; i < ncon; ++i) {
      double ftmp = (pweights[i] + pweights[ncon+i])/2.0;
      lbvec[i] = ((std::abs(ftmp) < 1.e-10) ? 0. :
		  std::abs(graph->gnpwgts[i] - target_weights[i]) / ftmp);
    }
    double lbavg = average(ncon, lbvec.data());

    int total_vsize = 0;
    for(int i=0; i < nvtxs; ++i)
      if( where[i] != home[i] )
        total_vsize += graph->vsize[i];

    double mycost = graph->mincut*ctrl->itr_ratio + total_vsize*ctrl->redist_ratio;

    if( bestcost >= mycost ) {
      bestcost = mycost;
      other_lbavg = lbavg;
      costwhere.assign(where, where + nvtxs);
    }

    if( best_lbavg >= lbavg ) {
      best_lbavg = lbavg;
      other_cost = mycost;
      lbwhere.assign(where, where + nvtxs);
    }
  }

  if( other_lbavg <= .05 ) {
    icopy(nvtxs, costwhere.data(), where);
    *sr_cost = bestcost;
    *sr_lbavg = other_lbavg;
  } else {
    icopy(nvtxs, lbwhere.data(), where);
    *sr_cost = other_cost;
    *sr_lbavg = best_lbavg;
  }
}


//------------------------------------------------------------------------
// Entry point of the initial partitioning algorithm which assembles
// the graph to all the cpus and proceeds serially
void process_diffusion(Control *ctrl, Graph *graph, int *vtxdist, int *where, 
		       int *home, int npasses)
{
  if( graph->ncon > 3 )
    return;

  int   nvtxs  = graph->nvtxs; // convinience shortcut
  int   ncon   = graph->ncon;  // convinience shortcut
  int   nparts = ctrl->nparts; // convinience shortcut
  double ubavg = average(ncon, ctrl->disbalance.data());

  // Initialize variables and allocate memory
  std::vector<double> lbvec(ncon);
  std::vector<double> flows_diff(ncon);
  std::vector<double> flows_sr(ncon);

  std::vector<double> load(nparts);
  std::vector<double> solution(nparts);
  graph->gnpwgts.resize(ncon*nparts);
  std::vector<double> values(graph->nedges);
  std::vector<double> transfer(ncon*graph->nedges);

  std::vector<int> proc2sub(std::max(nparts, 2*ctrl->ncpus));
  std::vector<int> sub2proc(nparts);
  std::vector<int> match(nparts);
  std::vector<int> rowptr(nparts+1);
  std::vector<int> colind(graph->nedges);
  
  std::vector<int> rcount(ctrl->ncpus);
  std::vector<int> rdispl(ctrl->ncpus+1);

  std::vector<int> pack(nvtxs);
  std::vector<int> unpack(nvtxs);
  std::vector<int> rbuffer(nvtxs);
  std::vector<int> sbuffer(nvtxs);
  std::vector<int> map(nvtxs);
  std::vector<int> rmap(nvtxs);
  std::vector<int> diff_where(nvtxs);
  std::vector<int> ehome(nvtxs);

  size_t wsize = std::max(sizeof(double)*nparts*6, sizeof(int)*(nvtxs+nparts*2+1));
  std::vector<double> workspace(wsize);

  graph->ckrinfo.resize(nvtxs);

  // Construct subdomain connectivity graph
  matrix_t matrix;
  matrix.nrows    = nparts;
  matrix.values   = values.data();
  matrix.transfer = transfer.data();
  matrix.rowptr   = rowptr.data();
  matrix.colind   = colind.data();
  setup_matrix(graph, &matrix, (int *)workspace.data());
  int nlinks = (matrix.nnzs - nparts) / 2;

  std::vector<int> lvisited(matrix.nnzs);
  std::vector<int> gvisited(matrix.nnzs);
  int nvisited;
  for(int pass=0; pass < npasses; ++pass) {
    rset(matrix.nnzs*ncon, 0.0, transfer.data());
    iset(matrix.nnzs, 0, gvisited.data());
    iset(matrix.nnzs, 0, lvisited.data());
    int iter = 0;
    nvisited = 0;

    // Compute ncon flow solutions
    for(int i=0; i < ncon; ++i) {
      rset(nparts, 0.0, solution.data());
      compute_load(graph, nparts, load.data(), ctrl->part_weights.data(), i);

      lbvec[i] = (rmax(nparts, load.data(), 1)+1.0/nparts) * nparts;

      cg_solve(&matrix, load.data(), solution.data(), 0.001, workspace.data());
      compute_transfer(ncon, &matrix, solution.data(), transfer.data(), i);
    }

    double lbavg_old = average(ncon, lbvec.data());
    int nswaps = -1, allnswaps = 0;
    double max_diff = 0.;
    for(int i=0; i < nparts; ++i) {
      for(int j=rowptr[i]; j < rowptr[i+1]; ++j) {
        double maxflow = rmax(ncon, &transfer[j*ncon], 1);
        double minflow = rmin(ncon, &transfer[j*ncon], 1);
	max_diff = std::max(max_diff, maxflow - minflow);
      }
    }

    while( nvisited < nlinks ) {
      // Compute independent sets of subdomains
      iset(std::max(nparts, 2*ctrl->ncpus), UNMATCHED, proc2sub.data());
      find_match(&matrix, match.data(), proc2sub.data(), gvisited.data(), ncon);

      // Set up the packing arrays
      iset(nparts, UNMATCHED, sub2proc.data());
      for(int i=0; i < 2*ctrl->ncpus; ++i) {
        if( proc2sub[i] == UNMATCHED )
          break;

        sub2proc[proc2sub[i]] = i/2;
      }

      iset(ctrl->ncpus, 0, rcount.data());
      for(int i=0; i < nvtxs; ++i) {
        int cpu = sub2proc[where[i]];
        if( cpu != UNMATCHED ) 
          rcount[cpu]++;
      }

      rdispl[0] = 0;
      for(int i=1; i < ctrl->ncpus+1; ++i)
        rdispl[i] = rdispl[i-1] + rcount[i-1];

      iset(nvtxs, UNMATCHED, unpack.data());
      for(int i=0; i < nvtxs; ++i) {
        int cpu = sub2proc[where[i]];
        if( cpu != UNMATCHED ) 
          unpack[rdispl[cpu]++] = i;
      }

      for(int i=ctrl->ncpus; i > 0; --i) rdispl[i] = rdispl[i-1];
      rdispl[0] = 0;

      iset(nvtxs, UNMATCHED, pack.data());
      for(int i=0; i < rdispl[ctrl->ncpus]; ++i) {
        assert(unpack[i] != UNMATCHED);
        int domain = where[unpack[i]];
        int cpu = sub2proc[domain];
        if( cpu != UNMATCHED ) 
          pack[unpack[i]] = i;
      }

      // Compute the flows
      if( proc2sub[2*ctrl->myrank] != UNMATCHED ) {
        int me  = proc2sub[2*ctrl->myrank];
        int you = proc2sub[2*ctrl->myrank+1];
        assert(me != you);

        for(int j=rowptr[me]; j < rowptr[me+1]; ++j) {
          if( colind[j] == you ) {
            lvisited[j] = 1;
            rcopy(ncon, &transfer[j*ncon], flows_diff.data());
            break;
          }
        }

        for(int j=rowptr[you]; j < rowptr[you+1]; ++j) {
          if( colind[j] == me ) {
            lvisited[j] = 1;
            for(int ii=0; ii < ncon; ++ii) {
              if( transfer[j*ncon+ii] > 0.0 )
                flows_diff[ii] = -transfer[j*ncon+ii];
            }
            break;
          }
        } 

        nswaps = 1;
        rcopy(ncon, flows_diff.data(), flows_sr.data());

        iset(nvtxs, 0, sbuffer.data());
        for(int i=0; i < nvtxs; ++i) {
          if( where[i] == me || where[i] == you )
            sbuffer[i] = 1;
        }

        Graph *egraph = extract_graph(ctrl, graph, sbuffer.data(), map.data(), rmap.data());

        if( egraph != NULL ) {
          icopy(egraph->nvtxs, egraph->where.data(), diff_where.data()); //<<<<<
          for(int j=0; j < egraph->nvtxs; ++j)
            ehome[j] = home[map[j]];

	  double lbavg_sr, cost_sr;
          refine_links(ctrl, egraph, ehome.data(), me, you, flows_sr.data(), &cost_sr, &lbavg_sr);

          if( ncon <= 4 ) {
            std::vector<int> sr_where(egraph->where);
            egraph->where = diff_where;

	    double lbavg_diff, cost_diff;
            nswaps = balance_links(ctrl, egraph, ehome.data(), me, you, flows_diff.data(), max_diff, 
				   &cost_diff, &lbavg_diff, 1.0/(double)nvtxs);
	    
            if( (lbavg_sr < lbavg_diff &&
		 (lbavg_diff >= ubavg-1.0 || cost_sr == cost_diff)) ||
                (lbavg_sr < ubavg-1.0 && cost_sr < cost_diff) ) {
              for(int i=0; i < egraph->nvtxs; ++i)
                where[map[i]] = sr_where[i];
            } else {
              for(int i=0; i < egraph->nvtxs; ++i)
                where[map[i]] = diff_where[i];
            }
          } else {
            for(int i=0; i < egraph->nvtxs; ++i)
              where[map[i]] = egraph->where[i];
          }

          delete egraph;
        }

        // Pack the flow data
        iset(nvtxs, UNMATCHED, sbuffer.data());
        for(int i=0; i < nvtxs; ++i) {
          int domain = where[i];
          if( domain == you || domain == me )
            sbuffer[pack[i]] = where[i];
        }
      }

      // Broadcast the flow data
      MPI_Allgatherv((void *)&sbuffer[rdispl[ctrl->myrank]], rcount[ctrl->myrank], MPI_INT, 
		     (void *)rbuffer.data(), rcount.data(), rdispl.data(), MPI_INT, ctrl->comm);

      // Unpack the flow data
      for(int i=0; i < rdispl[ctrl->ncpus]; ++i) {
        if( rbuffer[i] != UNMATCHED )
          where[unpack[i]] = rbuffer[i];
      }

      // Do other stuff
      MPI_Allreduce((void *)lvisited.data(), (void *)gvisited.data(), matrix.nnzs, MPI_INT,
		    MPI_MAX, ctrl->comm);
      nvisited = isum(matrix.nnzs, gvisited.data(), 1)/2;

      int swap_sum;
      MPI_Allreduce((void *)&nswaps, (void *)&swap_sum, 1, MPI_INT, MPI_SUM, ctrl->comm);

      allnswaps += swap_sum;

      if( iter++ >= NGD_STEPS )
        break;
    }

    // Perform serial refinement
    get_part_params_serial(ctrl, graph, nparts);
    refine_serial_kw(ctrl, graph, nparts, home, ctrl->disbalance.data(), 10);

    // Check for early breakout
    for(int i=0; i < ncon; ++i) {
      lbvec[i] = (double)(nparts) *
        graph->gnpwgts[rargmax_strd(nparts,&graph->gnpwgts[i],ncon)*ncon+i];
    }
    double lbavg = average(ncon, lbvec.data());

    int done = 0;
    if( allnswaps == 0 || lbavg >= lbavg_old || lbavg <= ubavg + AVG_THRESHOLD )
      done = 1;

    int alldone;
    MPI_Allreduce((void *)&done, (void *)&alldone, 1, MPI_INT, MPI_MAX, ctrl->comm);
    if( alldone == 1 )
      break;
  }

  // Clean up
  graph->ckrinfo.clear();
}


//------------------------------------------------------------------------
// Compute the edge-cut of a serial graph
int get_edge_cut(Graph *graph)
{
  int cut = 0;
  for(int i=0; i < graph->nvtxs; ++i) {
    for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j)
      if( graph->where[i] != graph->where[graph->adjncy[j]] )
        cut += graph->adjwgt[j];
  }
  graph->mincut = cut/2;

  return graph->mincut;
}


//------------------------------------------------------------------------
// Compute the indices of top three values of a double array
void get_3_max(int n, double *x, int *first, int *second, int *third)
{
  if( n <= 0 ) {
    *first = *second = *third = -1;
    return;
  }

  *second = *third = -1;
  *first = 0;

  for(int i=1; i < n; ++i) {
    if( x[i] > x[*first] ) {
      *third = *second;
      *second = *first;
      *first = i;
      continue;
    }

    if( *second == -1 || x[i] > x[*second] ) {
      *third = *second;
      *second = i;
      continue;
    }

    if( *third == -1 || x[i] > x[*third] )
      *third = i;
  }
}

//------------------------------------------------------------------------
// Perform directed diffusion
double frontal_diffusion(Control *ctrl, Graph *graph, int *home)
{
  int  nparts = ctrl->nparts;  // convinience shortcut
  int  nvtxs  = graph->nvtxs;  // convinience shortcut
  int  nedges = graph->nedges; // convinience shortcut
  int *where  = graph->where.data();  // convinience shortcut

  double flowFactor = 0.35;
  flowFactor = (ctrl->myrank == 2) ? 0.50 : flowFactor;
  flowFactor = (ctrl->myrank == 3) ? 0.75 : flowFactor;
  flowFactor = (ctrl->myrank == 4) ? 1.00 : flowFactor;

  matrix_t matrix;
  matrix.nrows = nparts;

  // Allocate memory
  std::vector<double> rspace(6*nparts+2*nedges);
  double *solution = rspace.data();
  double *tmpvec   = solution + nparts;           // nparts
  double *npwgts   = solution + 2*nparts;         // nparts
  double *load     = solution + 3*nparts;         // nparts
  matrix.values    = solution + 4*nparts;         // nparts+nedges
  matrix.transfer  = solution + 5*nparts + nedges;// nparts+nedges

  std::vector<int> ispace(2*nvtxs+3*nparts+nedges+1);
  int *perm     = ispace.data();
  int *ed       = perm + nvtxs;                   // nvtxs
  int *psize    = perm + 2*nvtxs;                 // nparts
  matrix.rowptr = perm + 2*nvtxs + nparts;        // nparts+1
  matrix.colind = perm + 2*nvtxs + 2*nparts + 1;  // nparts+nedges

  std::vector<double> rwspace(nparts*6);
  std::vector<int>    iwspace(nparts*2+nvtxs+1);
  std::vector<ikv_t>  candidate(nvtxs);

  // Populate empty subdomains
  iset(nparts, 0, psize);
  for(int i=0; i < nvtxs; ++i)
    psize[where[i]]++;

  bool empty_part_exist = false;
  for(int i=0; i < nparts; ++i) {
    if( psize[i] == 0 ) {
      empty_part_exist = true;
      break;
    }
  }

  if( empty_part_exist ) {
    random_permute_fast(nvtxs, perm, 1);
    for(int mind=0; mind < nparts; ++mind) {
      if( psize[mind] > 0 )
        continue;

      int maxd = (int)iargmax(nparts, psize, 1);
      if( psize[maxd] == 1 )
        break;  // Pass out if the heaviest subdomain contains only one vertex

      for(int i=0; i < nvtxs; ++i) {
        int k = perm[i];
        if( where[k] == maxd ) {
          where[k] = mind;
          psize[mind]++;
          psize[maxd]--;
          break;
        }
      }
    }
  }

  // Compute external degrees of the vertices
  iset(nvtxs, 0, ed);
  rset(nparts, 0.0, npwgts);
  for(int i=0; i < nvtxs; ++i) {
    npwgts[where[i]] += graph->nvwgt[i];
    for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j)
      ed[i] += (where[i] != where[graph->adjncy[j]] ? graph->adjwgt[j] : 0);
  }

  compute_load(graph, nparts, load, ctrl->part_weights.data(), 0);

  rset(nparts, 0., tmpvec);

  int npasses = std::min(nparts/2, NGD_STEPS);
  for(int done=0, pass=0; pass < npasses; ++pass) {
    // Set-up and solve the diffusion equation
    int nswaps = 0;

    // Solve flow equations
    setup_matrix(graph, &matrix, iwspace.data());

    // Check for disconnected subdomains
    bool is_connected = true;
    for(int i=0; i < matrix.nrows; ++i) {
      if( matrix.rowptr[i]+1 == matrix.rowptr[i+1] ) {
	is_connected = false;
        break;
      }
    }

    if( is_connected ) {
      const double tolerance = 0.001;
      cg_solve(&matrix, load, solution, tolerance, rwspace.data());
      compute_transfer(1, &matrix, solution,  matrix.transfer, 0);

      int first, second, third;
      get_3_max(nparts, load, &first, &second, &third);
  
      if( pass % 3 == 0 ) {
        random_permute_fast(nvtxs, perm, 1);
      } else {
        // Move dirty vertices first
        int ndirty = 0;
        for(int i=0; i < nvtxs; ++i) {
          if( where[i] != home[i] )
            ndirty++;
        }
  
        int dptr = 0;
        for(int i=0; i < nvtxs; ++i) {
          if( where[i] != home[i] ) perm[dptr++]   = i;
          else                      perm[ndirty++] = i;
        }

        ASSERT(ndirty == nvtxs);
        ndirty = dptr;
        int nclean = nvtxs - dptr;
        random_permute_fast(ndirty, perm, 0);
        random_permute_fast(nclean, perm+ndirty, 0);
      }
  
      if( ctrl->myrank == 0 ) {
	int k = 0;
        for(int j=nvtxs, ii=0; ii < nvtxs; ++ii) {
          int i = perm[ii];
          if( ed[i] != 0 ) {
            candidate[k].key = -ed[i];
            candidate[k++].val = i;
          } else {
            candidate[--j].key = 0;
            candidate[j].val = i;
          }
        }
        ikvsorti(k, candidate.data());
      }

      for(int ii=0; ii < nvtxs/3; ++ii) {
        int i = (ctrl->myrank == 0) ? candidate[ii].val : perm[ii];
        int from = where[i];
  
        // Subdomain should have at least one vertex left
        if( psize[from] == 1 )
          continue;
  
        bool clean = (from == home[i]);
  
        // Only move from top three or dirty vertices
        if( from != first && from != second && from != third && clean )
          continue;
  
        // Scatter the sparse transfer row into the dense tmpvec row
        for(int j=matrix.rowptr[from]+1; j < matrix.rowptr[from+1]; ++j)
          tmpvec[matrix.colind[j]] =  matrix.transfer[j];
  
        for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
          int to = where[graph->adjncy[j]];
          if( from != to ) {
            if( tmpvec[to] > (flowFactor * graph->nvwgt[i]) ) {
	      double weight = graph->nvwgt[i];
              tmpvec[to]   -= weight;
	      npwgts[to]   += weight;
	      npwgts[from] -= weight;
	      load[to]     += weight;
	      load[from]   -= weight;
	      psize[to]    += 1;
	      psize[from]  -= 1;
              where[i] = to;
              ++nswaps;
  
              // Update external degrees
              ed[i] = 0;
              for(int k=graph->xadj[i]; k < graph->xadj[i+1]; ++k) {
                int edge = graph->adjncy[k];
                ed[i] += (to != where[edge] ? graph->adjwgt[k] : 0);
  
                if( where[edge] == from ) ed[edge] += graph->adjwgt[k];
                if( where[edge] == to   ) ed[edge] -= graph->adjwgt[k];
              }
              break;
            }
          }
        }
  
        // Gather the dense tmpvec row into the sparse transfer row
        for(int j=matrix.rowptr[from]+1; j < matrix.rowptr[from+1]; ++j) {
	  matrix.transfer[j] = tmpvec[matrix.colind[j]];
          tmpvec[matrix.colind[j]] = 0.;
        }
      }
    }

    if( pass % 2 == 1 ) {
      double ubfactor = ctrl->disbalance[0];
      double balance  = rmax(nparts, npwgts, 1)*nparts;
      if( balance < ubfactor + AVG_THRESHOLD )
        done = 1;

      int done_any;
      MPI_Allreduce(&done, &done_any, 1, MPI_INT, MPI_MAX, ctrl->comm);
      if( done_any > 0 )
        break;

      int noswaps = (nswaps > 0 ? 0 : 1);
      int noswaps_total;
      MPI_Allreduce(&noswaps, &noswaps_total, 1, MPI_INT, MPI_SUM, ctrl->comm);
      if( noswaps_total > ctrl->ncpus/2 )
        break;
    }
  }

  graph->mincut = get_edge_cut(graph);
  int totalv  = compute_totalw(graph, home);
  double cost   = ctrl->itr_ratio * graph->mincut + ctrl->redist_ratio * totalv;

  return cost;
}


//------------------------------------------------------------------------
// Assemble graph on a single processor
Graph *assemble_adaptive_graph(Control *ctrl, Graph *graph)
{
  int gnvtxs  = graph->gnvtxs; // convinience shortcut
  int nvtxs   = graph->nvtxs;  // convinience shortcut
  int ncon    = graph->ncon;   // convinience shortcut
  int nedges  = graph->xadj[nvtxs];

  // Determine the number of itemst to receive from each cpu
  std::vector<int> rcounts(ctrl->ncpus);

  int mysize = (1+ncon)*nvtxs + 2*nedges;
  if( ctrl->is_adaptive )
    mysize = (2+ncon)*nvtxs + 2*nedges;

  MPI_Allgather(&mysize, 1, MPI_INT, rcounts.data(), 1, MPI_INT, ctrl->comm);

  std::vector<int> rdispls(ctrl->ncpus+1);
  rdispls[0] = 0;
  for(int i=1; i < ctrl->ncpus+1; ++i)
    rdispls[i] = rdispls[i-1] + rcounts[i-1];

  // Allocate memory for recv buffer of the assembled graph
  int gsize = rdispls[ctrl->ncpus];
  std::vector<int> ggraph(gsize);

  // Construct the array storage for assembled graph
  std::vector<int> mygraph(mysize);

  for(int k=0, i=0; i < nvtxs; ++i) {
    mygraph[k++] = graph->xadj[i+1] - graph->xadj[i];
    for(int j=0; j < ncon; ++j)
      mygraph[k++] = graph->vwgt[i*ncon+j];
    if( ctrl->is_adaptive )
      mygraph[k++] = graph->vsize[i];
    for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
      mygraph[k++] = graph->imap[graph->adjncy[j]];
      mygraph[k++] = graph->adjwgt[j];
    }
  }

  // Assemble the entire graph
  MPI_Allgatherv(mygraph.data(), mysize, MPI_INT, ggraph.data(), 
		 rcounts.data(), rdispls.data(), MPI_INT, ctrl->comm);

  Graph *agraph = create_graph();
  agraph->nvtxs = gnvtxs;
  agraph->ncon  = ncon;

  if( ctrl->is_adaptive ) agraph->nedges = (gsize-(2+ncon)*gnvtxs)/2;
  else                    agraph->nedges = (gsize-(1+ncon)*gnvtxs)/2;
  int gnedges = agraph->nedges;

  // Allocate memory for the assembled graph
  agraph->xadj  .resize(gnvtxs+1);
  agraph->adjncy.resize(gnedges);
  agraph->nvwgt .resize(gnvtxs*ncon);
  agraph->vwgt  .resize(gnvtxs*ncon);
  agraph->adjwgt.resize(gnedges);

  agraph->label.resize(gnvtxs);
  if( ctrl->is_adaptive )
    agraph->vsize.resize(gnvtxs);

  for(int k=0, j=0, i=0; i < gnvtxs; ++i) {
    agraph->xadj[i] = ggraph[k++];
    for(int jj=0; jj < ncon; ++jj)
      agraph->vwgt[i*ncon+jj] = ggraph[k++];
    if( ctrl->is_adaptive )
      agraph->vsize[i] = ggraph[k++];
    for(int jj=0; jj < agraph->xadj[i]; ++jj) {
      agraph->adjncy[j] = ggraph[k++];
      agraph->adjwgt[j] = ggraph[k++];
      ++j;
    }
  }

  // Reconcile the received graph
  for(int i=1; i < gnvtxs; ++i) agraph->xadj[i] += agraph->xadj[i-1];
  for(int i=gnvtxs; i > 0; --i) agraph->xadj[i]  = agraph->xadj[i-1];
  agraph->xadj[0] = 0;

  for(int i=0; i < gnvtxs; ++i) {
    for(int j=0; j < ncon; ++j)
      agraph->nvwgt[i*ncon+j] = ctrl->inv_weights[j] * agraph->vwgt[i*ncon+j];
  }

  iincset(gnvtxs, 0, &agraph->label[0]);

  return agraph;
}

//------------------------------------------------------------------------
// Entry point of the initial balancing algorithm.
// The algorithm assembles the graph on all cpus and proceeds with balancing.
void balance_partition(Control *ctrl, Graph *graph)
{
  struct {
    double cost;
    int rank;
  } lpecost, gpecost;

  static const int    SPLIT_FACTOR      = 8;
  static const int    DIFF_STEPS        = 6;
  static const int    MAX_NCON          = 2;  // For diffusion
  static const double REDISTRIBUTE_WGT  = 2.0;

  int *vtxdist = graph->vtxdist.data();
  Graph *agraph = assemble_adaptive_graph(ctrl, graph);

  Graph cgraph;
  int nvtxs   = cgraph.nvtxs  = agraph->nvtxs;
  int nedges  = cgraph.nedges = agraph->nedges;
  int ncon    = cgraph.ncon   = agraph->ncon;

  cgraph.xadj  .assign(agraph->xadj  .cbegin(), agraph->xadj  .cbegin() + nvtxs + 1);
  cgraph.vwgt  .assign(agraph->vwgt  .cbegin(), agraph->vwgt  .cbegin() + nvtxs*ncon);
  cgraph.vsize .assign(agraph->vsize .cbegin(), agraph->vsize .cbegin() + nvtxs);
  cgraph.adjncy.assign(agraph->adjncy.cbegin(), agraph->adjncy.cbegin() + nedges);
  cgraph.adjwgt.assign(agraph->adjwgt.cbegin(), agraph->adjwgt.cbegin() + nedges);

  agraph->where.resize(nvtxs);
  int *part = agraph->where.data();

  std::vector<int>  lwhere(nvtxs);
  std::vector<int>  home  (nvtxs);
  std::vector<double> lbvec (graph->ncon);

  //-------------------------------------------------------------
  if( ctrl->is_coupled ) {
    for(int i=0; i < ctrl->ncpus; ++i) {
      for(int j=vtxdist[i]; j < vtxdist[i+1]; ++j)
        part[j] = home[j] = i;
    }
  } else {
    std::vector<int> rcounts(ctrl->ncpus);
    std::vector<int> rdispls(ctrl->ncpus+1);

    for(int i=0; i < ctrl->ncpus; ++i) 
      rdispls[i] = rcounts[i] = vtxdist[i+1] - vtxdist[i];

    for(int i=1; i < ctrl->ncpus; ++i) rdispls[i] += rdispls[i-1];
    for(int i=ctrl->ncpus; i > 0; --i) rdispls[i]  = rdispls[i-1];
    rdispls[0] = 0;

    MPI_Allgatherv(graph->home.data(), graph->nvtxs, MPI_INT,
		   part, rcounts.data(), rdispls.data(), MPI_INT, ctrl->comm);

    for(int i=0; i < agraph->nvtxs; ++i)
      home[i] = part[i];
  }

  // Fix partitioning
  for(int i=0; i < agraph->nvtxs; ++i) {
    if( part[i] >= ctrl->nparts )
      part[i] = home[i] = part[i] % ctrl->nparts;
    if( part[i] < 0 )
      part[i] = home[i] = (-1*part[i]) % ctrl->nparts;
  }

  //-------------------------------------------------------------
  // Split cpus into 2 groups
  int sr = (ctrl->myrank % 2 == 0) ? 1 : 0;
  int gd = (ctrl->myrank % 2 == 1) ? 1 : 0;

  if( graph->ncon > MAX_NCON || ctrl->ncpus == 1 ) {
    sr = 1;
    gd = 0;
  }

  int sr_pe = 0;
  int gd_pe = 1;

  int myrank_sr, myrank_gd;
  int ncpus_sr, ncpus_gd;
  MPI_Comm ipcomm;
  MPI_Comm_split(ctrl->gcomm, sr, 0, &ipcomm);
  MPI_Comm_rank(ipcomm, &myrank_gd);
  MPI_Comm_size(ipcomm, &ncpus_gd);

  if( sr == 1 ) { // Half of cpus do scratch-remap
    int edgecut;
    int ngroups = std::max(std::min(SPLIT_FACTOR, ncpus_gd), 1);
    MPI_Comm srcomm;
    MPI_Comm_split(ipcomm, myrank_gd % ngroups, 0, &srcomm);
    MPI_Comm_rank(srcomm, &myrank_sr);
    MPI_Comm_size(srcomm, &ncpus_sr);

    int moptions[METIS_NOPTIONS];
    METIS_SetDefaultOptions(moptions);
    moptions[METIS_OPTION_SEED] = ctrl->global_seed + (myrank_gd % ngroups) + 1;

    iset(nvtxs, 0, lwhere.data());
    int lnparts = ctrl->nparts;
    int fpart = 0, fpe = 0;

    double *target_weights1 = ctrl->part_weights.data();
    std::vector<double> target_weights2(2*ncon);
    int ncpus_l = ncpus_sr;
    int twoparts = 2;
    while( ncpus_l > 1 && lnparts > 1 ) {
      ASSERT(agraph->nvtxs > 1);
      // Determine the weights
      for(int j=(lnparts>>1), i=0; i < ncon; ++i) {
        target_weights2[i]      = rsum(j, target_weights1+fpart*ncon+i, ncon);
        target_weights2[ncon+i] = rsum(lnparts-j, target_weights1+(fpart+j)*ncon+i, ncon);
        double wsum             = 1.0/(target_weights2[i] + target_weights2[ncon+i]);
        target_weights2[i]      *= wsum;
        target_weights2[ncon+i] *= wsum;
      }

      METIS_PartGraphRecursive(&agraph->nvtxs, &ncon, agraph->xadj.data(), 
                               agraph->adjncy.data(), agraph->vwgt.data(), NULL, agraph->adjwgt.data(),
                               &twoparts, target_weights2.data(),
                               NULL, moptions, &edgecut, part);

      // Pick one of the branches
      if( myrank_sr < fpe+ncpus_l/2 ) {
        keep_parts(agraph, part, 0);
        ncpus_l /= 2;
        lnparts = lnparts/2;
      }
      else {
        keep_parts(agraph, part, 1);
        fpart   = fpart + lnparts/2;
        fpe     = fpe + ncpus_l/2;
        ncpus_l-= ncpus_l/2;
        lnparts = lnparts - lnparts/2;
      }
    }

    if( lnparts == 1 ) { // Case in which ncpus_sr >= nparts
      // Only the first cpu will assign labels
      if( myrank_sr == fpe ) {
        for(int i=0; i < agraph->nvtxs; ++i)
          lwhere[agraph->label[i]] = fpart;
      }
    } else { // Case in which ncpus_sr < nparts
      // Create normalized target_weights for the lnparts
      std::vector<double> target_weights(lnparts*ncon);
      for(int j=0; j < ncon; ++j) {
	double wsum = 0.;
        for(int i=0; i < lnparts; ++i) {
          target_weights[i*ncon+j] = ctrl->part_weights[(fpart+i)*ncon+j];
          wsum += target_weights[i*ncon+j];
        }
        for(int i=0; i < lnparts; ++i)
          target_weights[i*ncon+j] /= wsum;
      }

      METIS_PartGraphKway(&agraph->nvtxs, &ncon, agraph->xadj.data(), agraph->adjncy.data(),
                          agraph->vwgt.data(), NULL, agraph->adjwgt.data(), &lnparts,
                          target_weights.data(), NULL, moptions, &edgecut, part);

      for(int i=0; i < agraph->nvtxs; ++i)
        lwhere[agraph->label[i]] = fpart + part[i];
    }

    MPI_Allreduce(lwhere.data(), part, nvtxs, MPI_INT, MPI_SUM, srcomm);

    cgraph.where = agraph->where;
    edgecut = get_edge_cut(&cgraph);
    get_balance_local(ctrl, &cgraph, part, lbvec.data());

    double lbsum = rsum(ncon, lbvec.data(), 1);
    int max_cut;
    double min_lbsum;
    MPI_Allreduce(&edgecut, &max_cut, 1, MPI_INT, MPI_MAX, ipcomm);
    MPI_Allreduce(&lbsum, &min_lbsum, 1, MPI_DOUBLE, MPI_MIN, ipcomm);
    lpecost.rank = ctrl->myrank;
    lpecost.cost = lbsum;
    if( min_lbsum < DISBALANCE_RATIO * (double)(ncon) ) {
      if( lbsum < DISBALANCE_RATIO * (double)(ncon) )
        lpecost.cost = (double)edgecut;
      else
        lpecost.cost = (double)max_cut + lbsum;
    }
    MPI_Allreduce(&lpecost, &gpecost, 1, MPI_DOUBLE_INT, MPI_MINLOC, ipcomm);

    if( ctrl->myrank == gpecost.rank && ctrl->myrank != sr_pe )
      MPI_Send(part, nvtxs, MPI_INT, sr_pe, 1, ctrl->comm);

    MPI_Status status;
    if( ctrl->myrank != gpecost.rank && ctrl->myrank == sr_pe )
      MPI_Recv(part, nvtxs, MPI_INT, gpecost.rank, 1, ctrl->comm, &status);

    if( ctrl->myrank == sr_pe ) {
      icopy(nvtxs, part, lwhere.data());
      remap_serial(&cgraph, ctrl->nparts, home.data(), lwhere.data(),
                   part, ctrl->part_weights.data());
    }

    MPI_Comm_free(&srcomm);
  } else { // The other half do global diffusion
    // Setup ctrl for the diffusion
    Control *myctrl = new Control();

    myctrl->myrank        = myrank_gd;
    myctrl->ncpus         = ncpus_gd;
    myctrl->free_gcomm    = false;
    myctrl->comm          = ipcomm;
    myctrl->global_seed   = ctrl->global_seed;
    myctrl->local_seed    = ctrl->local_seed;
    myctrl->nparts        = ctrl->nparts;
    myctrl->nconstraints  = ctrl->nconstraints;
    myctrl->itr_ratio     = ctrl->itr_ratio;
    myctrl->redist_ratio  = 1.;
    myctrl->is_adaptive   = true;
    myctrl->is_coupled    = false;

    myctrl->part_weights.resize(myctrl->nparts*myctrl->nconstraints);
    myctrl->inv_weights.resize (myctrl->nconstraints);
    myctrl->disbalance .resize (myctrl->nconstraints);

    myctrl->part_weights.assign(ctrl->part_weights.data(),
                                ctrl->part_weights.data() + myctrl->nparts*myctrl->nconstraints);
    myctrl->inv_weights .assign(ctrl->inv_weights.data(),
                                ctrl->inv_weights.data() + myctrl->nconstraints);
    myctrl->disbalance  .assign(ctrl->disbalance.data(),
                                ctrl->disbalance.data() + myctrl->nconstraints);

    myctrl->pool.reserve(agraph->nvtxs);

    // This is required to balance out the sr MPI_Comm_split
    MPI_Comm srcomm;
    MPI_Comm_split(ipcomm, MPI_UNDEFINED, 0, &srcomm);

    if( ncon == 1 ) {
      double rating = frontal_diffusion(myctrl, agraph, home.data());
      get_balance_local(ctrl, &cgraph, part, lbvec.data());
      double lbsum = rsum(ncon, lbvec.data(), 1);

      // Determine which PE computed the best partitioning
      double max_rating, min_lbsum;
      MPI_Allreduce(&rating, &max_rating, 1, MPI_DOUBLE, MPI_MAX, ipcomm);
      MPI_Allreduce(&lbsum, &min_lbsum, 1, MPI_DOUBLE, MPI_MIN, ipcomm);

      lpecost.rank = ctrl->myrank;
      lpecost.cost = lbsum;
      if( min_lbsum < DISBALANCE_RATIO * (double)(ncon) ) {
        if( lbsum < DISBALANCE_RATIO * (double)(ncon) )
          lpecost.cost = rating;
        else
          lpecost.cost = max_rating + lbsum;
      }

      MPI_Allreduce(&lpecost, &gpecost, 1, MPI_DOUBLE_INT, MPI_MINLOC, ipcomm);

      // Now send this to the coordinating processor
      if( ctrl->myrank == gpecost.rank && ctrl->myrank != gd_pe )
        MPI_Send(part, nvtxs, MPI_INT, gd_pe, 1, ctrl->comm);

      MPI_Status status;
      if( ctrl->myrank != gpecost.rank && ctrl->myrank == gd_pe )
        MPI_Recv(part, nvtxs, MPI_INT, gpecost.rank, 1, ctrl->comm, &status);

      if( ctrl->myrank == gd_pe ) {
        icopy(nvtxs, part, lwhere.data());
        remap_serial(&cgraph, ctrl->nparts, home.data(), lwhere.data(),
                     part, ctrl->part_weights.data());
      }
    } else {
      process_diffusion(myctrl, agraph, graph->vtxdist.data(),
                        agraph->where.data(), home.data(), DIFF_STEPS);
    }

    free_controls(&myctrl);
  }

  int who_wins;
  if( graph->ncon <= MAX_NCON ) {
    double my_balance = -1., my_cost = -1.;
    double other_balance = -1., other_cost = -1.;
    if( ctrl->myrank == sr_pe  || ctrl->myrank == gd_pe ) {
      //-------------------------------------------------------------
      // Use balancing to decide on the best partitioning
      cgraph.where = agraph->where;
      double my_cut = (double) get_edge_cut(&cgraph);
      double my_totalv = (double) compute_totalw(&cgraph, home.data());
      get_balance_local(ctrl, &cgraph, part, lbvec.data());

      my_balance = rsum(cgraph.ncon, lbvec.data(), 1);
      my_balance /= (double) cgraph.ncon;
      my_cost = ctrl->itr_ratio * my_cut + REDISTRIBUTE_WGT * my_totalv;

      double buffer[2];
      if( ctrl->myrank == gd_pe ) {
        buffer[0] = my_cost;
        buffer[1] = my_balance;
        MPI_Send(buffer, 2, MPI_DOUBLE, sr_pe, 1, ctrl->comm);
      } else {
	MPI_Status status;
        MPI_Recv(buffer, 2, MPI_DOUBLE, gd_pe, 1, ctrl->comm, &status);
        other_cost    = buffer[0];
        other_balance = buffer[1];
      }
    }

    if( ctrl->myrank == sr_pe ) {
      who_wins = gd_pe;
      if( (my_balance < 1.1 && other_balance > 1.1) ||
          (my_balance < 1.1 && other_balance < 1.1 && my_cost < other_cost) ||
          (my_balance > 1.1 && other_balance > 1.1 && my_balance < other_balance) ) {
        who_wins = sr_pe;
      }
    }

    MPI_Bcast(&who_wins, 1, MPI_INT, sr_pe, ctrl->comm);
  } else {
    who_wins = sr_pe;
  }

  MPI_Bcast(part, nvtxs, MPI_INT, who_wins, ctrl->comm);
  icopy(graph->nvtxs, part+vtxdist[ctrl->myrank], graph->where.data());

  MPI_Comm_free(&ipcomm);

  agraph->where.clear();
  free_graph(agraph);
}


//------------------------------------------------------------------------
// Do a binary search on an array for given key and return its index
int binary_search(int n, int *array, int key)
{
  int c, a = 0, b = n;

  while(b-a > 8) {
    c = (a+b)>>1;
    if( array[c] > key )
      b = c;
    else
      a = c;
  }

  for(c=a; c < b; ++c) {
    if( array[c] == key )
      return c;
  }

  std::runtime_error("Key not found!");
  return 0;
}

//------------------------------------------------------------------------
// Check if (v+u2) provides a better balance in the weight vector than (v+u1)
bool is_wbalance_better(int ncon, double *vwgt, double *u1wgt, double *u2wgt)
{
  if( ncon == 1 )
    return (u1wgt[0] - u1wgt[0]);

  double sum1 = 0., sum2 = 0.;
  for(int i=0; i < ncon; ++i) {
    sum1 += vwgt[i] + u1wgt[i];
    sum2 += vwgt[i] + u2wgt[i];
  }
  sum1 = sum1 / ncon;
  sum2 = sum2 / ncon;

  double diff1 = 0., diff2 = 0.;
  for(int i=0; i < ncon; ++i) {
    diff1 += std::abs(sum1 - (vwgt[i] + u1wgt[i]));
    diff2 += std::abs(sum2 - (vwgt[i] + u2wgt[i]));
  }

  return (diff1 - diff2) >= 0.;
}

//------------------------------------------------------------------------
#define HASH_SIZE  8192 // Should be a power of two
#define HASH_MASK  8191 // Should be equal to HASH_SIZE-1
#define KEEP_FLAG 0x40000000L
//------------------------------------------------------------------------
// Create the coarser graph after global matching
void create_cgraph_global(Control *ctrl, Graph *graph, int cnvtxs)
{
  int  nvtxs   = graph->nvtxs;   // convinience shortcut
  int  ncon    = graph->ncon;    // convinience shortcut
  int  nnbrs   = graph->nnbrs;   // convinience shortcut
  int *match   = graph->match.data(); // convinience shortcut
  int *peind   = graph->peind.data(); // convinience shortcut
  int *recvind = graph->recvind.data(); // convinience shortcut

  int first_vtx = graph->vtxdist[ctrl->myrank];
  int last_vtx  = graph->vtxdist[ctrl->myrank+1];

  graph->cmap.assign(graph->nvtxs+graph->nrecv, -1);
  int *cmap = graph->cmap.data();

  // Initialize the coarser graph
  Graph *cgraph = create_graph();
  cgraph->nvtxs  = cnvtxs;
  cgraph->ncon   = ncon;
  cgraph->level  = graph->level+1;
  cgraph->finer  = graph;
  graph->coarser = cgraph;

  // Obtain the vtxdist of the coarser graph
  cgraph->vtxdist.resize(ctrl->ncpus+1);
  int *cvtxdist = cgraph->vtxdist.data();

  cvtxdist[ctrl->ncpus] = cnvtxs;  // Use last position in the cvtxdist as a temp buffer

  MPI_Allgather(&cvtxdist[ctrl->ncpus], 1, MPI_INT, cvtxdist, 1, MPI_INT, ctrl->comm);

  for(int i=1; i < ctrl->ncpus; ++i) cvtxdist[i] += cvtxdist[i-1];
  for(int i=ctrl->ncpus; i > 0; --i) cvtxdist[i] = cvtxdist[i-1];
  cvtxdist[0] = 0;

  cgraph->gnvtxs = cvtxdist[ctrl->ncpus];

  //-------------------------------------------------------------
  // Construct the cmap vector
  int cfirst_vtx = cvtxdist[ctrl->myrank];
  int clast_vtx  = cvtxdist[ctrl->myrank+1];

  cnvtxs = 0;
  // Create the cmap of what you know so far locally.
  for(int i=0; i < nvtxs; ++i) {
    if( match[i] >= KEEP_FLAG ) {
      int k = match[i] - KEEP_FLAG;
      if( k >= first_vtx && k < first_vtx+i )
        continue;  // Both i, k are local and i has been matched .

      cmap[i] = cfirst_vtx + cnvtxs;
      ++cnvtxs;
      if( (k != first_vtx + i) && (k >= first_vtx && k < last_vtx)) { // Local match
        cmap [k-first_vtx]  = cmap[i];
        match[k-first_vtx] += KEEP_FLAG;  // Add the KEEP_FLAG to simplify marking
      }
    }
  }
  ASSERT(cnvtxs == clast_vtx-cfirst_vtx);

  distribute_bdry_data(ctrl, graph, cmap, cmap+nvtxs);

  // Update cmap for locally stored vertices that will go away. 
  // The remote cpu assigned cmap for them.
  for(int i=0; i < nvtxs; ++i) {
    if( match[i] < KEEP_FLAG ) {
      cmap[i] = cmap[nvtxs+binary_search(graph->nrecv, recvind, match[i])];
    }
  }

  distribute_bdry_data(ctrl, graph, cmap, cmap+nvtxs);

  //-------------------------------------------------------------
  // Determine how many adjcency lists are needed to send/receive.
  // First setp: determine sizes
  int nsend = 0, nrecv = 0;
  for(int i=0; i < nvtxs; ++i) {
    if( match[i] < KEEP_FLAG ) {
      ++nsend;
    } else {
      int k = match[i] - KEEP_FLAG;
      if( k < first_vtx || k >= last_vtx )
        ++nrecv;
    }
  }

  // Candidates for send/receive
  std::vector<ikv_t> cand_s(nsend);
  graph->rcand.resize(nrecv);
  ikv_t *cand_r = graph->rcand.data();

  // Second step: place sizes in the appropriate arrays
  int nkeepsize = 0;
  nsend = nrecv = 0;
  for(int i=0; i < nvtxs; ++i) {
    if( match[i] < KEEP_FLAG ) {
      cand_s[nsend].key = match[i];
      cand_s[nsend].val = i;
      ++nsend;
    } else {
      nkeepsize += (graph->xadj[i+1] - graph->xadj[i]);

      int k = match[i] - KEEP_FLAG;
      if( k < first_vtx || k >= last_vtx ) {
        cand_r[nrecv].key = k;
        cand_r[nrecv].val = cmap[i] - cfirst_vtx;
        ASSERT(cand_r[nrecv].val>=0 && cand_r[nrecv].val<cnvtxs);
        ++nrecv;
      }
    }
  }

  //-------------------------------------------------------------
  // Determine how many lists and their sizes need to send/receive
  // for each of the neighboring cpus
  graph->rlens.resize(nnbrs+1);
  int *rlens = graph->rlens.data();
  graph->slens.resize(nnbrs+1);
  int *slens = graph->slens.data();

  std::vector<int> rsizes(nnbrs, 0);
  std::vector<int> ssizes(nnbrs, 0);

  // Take care of sending data
  ikvsortii(nsend, cand_s.data());
  slens[0] = 0;
  for(int k=0, i=0; i < nnbrs; ++i) {
    int other_last_vtx = graph->vtxdist[peind[i]+1];
    for(; k < nsend && cand_s[k].key < other_last_vtx; ++k)
      ssizes[i] += (graph->xadj[cand_s[k].val+1] - graph->xadj[cand_s[k].val]);
    slens[i+1] = k;
  }

  // Then take care of receiving data
  ikvsortii(nrecv, cand_r);
  rlens[0] = 0;
  for(int k=0, i=0; i < nnbrs; ++i) {
    int other_last_vtx = graph->vtxdist[peind[i]+1];
    for(; k < nrecv && cand_r[k].key < other_last_vtx; ++k);
    rlens[i+1] = k;
  }

  //-------------------------------------------------------------
  // Exchange the sizes info
  // Always issue the receives first
  for(int i=0; i < nnbrs; ++i) {
    if( rlens[i+1]-rlens[i] > 0 )
      MPI_Irecv(&rsizes[i], 1, MPI_INT, peind[i], 1, ctrl->comm, &ctrl->rreqs[i]);
  }

  // Issue the sends next
  for(int i=0; i < nnbrs; ++i) {
    if( slens[i+1]-slens[i] > 0 )
      MPI_Isend(&ssizes[i], 1, MPI_INT, peind[i], 1, ctrl->comm, &ctrl->sreqs[i]);
  }

  // Wait for the operations to finish
  for(int i=0; i < nnbrs; ++i) {
    if( rlens[i+1]-rlens[i] > 0 )
      MPI_Wait(&ctrl->rreqs[i], &ctrl->status);
  }
  for(int i=0; i < nnbrs; ++i) {
    if( slens[i+1]-slens[i] > 0 )
      MPI_Wait(&ctrl->sreqs[i], &ctrl->status);
  }

  //-------------------------------------------------------------
  // Allocate memory for received/sent graphs and start sending 
  // and receiving data.
  // rgraph and sgraph are used to facilitate single message exchange
  int nrecvsize = isum(nnbrs, rsizes.data(), 1);
  int nsendsize = isum(nnbrs, ssizes.data(), 1);
  std::vector<int> rgraph((4+ncon)*nrecv + 2*nrecvsize);

  {
    std::vector<int> sgraph((4+ncon)*nsend + 2*nsendsize);

    // Deal with the received portion first
    for(int k=0, i=0; i < nnbrs; ++i) {
      // Issue a receive only if something is coming in
      if( rlens[i+1] - rlens[i] > 0 ) {
	MPI_Irecv(&rgraph[k], (4+ncon)*(rlens[i+1]-rlens[i]) + 2*rsizes[i], 
		  MPI_INT, peind[i], 1, ctrl->comm, &ctrl->rreqs[i]);
	k += (4+ncon)*(rlens[i+1]-rlens[i]) + 2*rsizes[i];
      }
    }

    // Deal with the sent portion now
    for(int kk=0, k=0, i=0; i < nnbrs; ++i) {
      if( slens[i+1]-slens[i] > 0 ) {
	for(int j=slens[i]; j < slens[i+1]; ++j) {
	  int jj = cand_s[j].val;
	  sgraph[kk++] = first_vtx + jj;
	  sgraph[kk++] = graph->xadj[jj+1] - graph->xadj[jj];
	  for(int ii=0; ii < ncon; ++ii)
	    sgraph[kk++] = graph->vwgt[jj*ncon+ii];

	  if( ctrl->is_adaptive ) {
	    sgraph[kk++] = graph->vsize[jj];
	    sgraph[kk++] = graph->home [jj];
	  } else {
	    sgraph[kk++] = -1;
	    sgraph[kk++] = -1;
	  }

	  for(int ii=graph->xadj[jj]; ii < graph->xadj[jj+1]; ++ii) {
	    sgraph[kk++] = cmap[graph->adjncy[ii]];
	    sgraph[kk++] = graph->adjwgt[ii];
	  }
	}

	ASSERT(kk-k == (4+ncon)*(slens[i+1]-slens[i])+2*ssizes[i]);

	MPI_Isend(&sgraph[k], kk-k, MPI_INT, peind[i], 1, ctrl->comm, &ctrl->sreqs[i]);
	k = kk;
      }
    }

    // Wait for the operations to finish
    for(int i=0; i < nnbrs; ++i) {
      if( rlens[i+1]-rlens[i] > 0 )
	MPI_Wait(&ctrl->rreqs[i], &ctrl->status);
    }

    for(int i=0; i < nnbrs; ++i) {
      if( slens[i+1]-slens[i] > 0 )
	MPI_Wait(&ctrl->sreqs[i], &ctrl->status);
    }
  }

  //-------------------------------------------------------------
  // Setup mapping btw indices returned by binary_search and the stored ones
  std::vector<int> perm(graph->nrecv, -1);
  for(int j=0, i=0; i < nrecv; ++i) {
    perm[binary_search(graph->nrecv, recvind, rgraph[j])] = j+1;
    j += (4+ncon)+2*rgraph[j+1];
  }

  //-------------------------------------------------------------
  // Create the coarser graph
  cgraph->xadj.resize(cnvtxs+1);
  int *cxadj = cgraph->xadj.data();

  cgraph->vwgt.resize(cnvtxs*ncon);
  int *cvwgt = cgraph->vwgt.data();

  cgraph->nvwgt.resize(cnvtxs*ncon);
  double *cnvwgt = cgraph->nvwgt.data();

  if( ctrl->is_adaptive ) {
    cgraph->vsize.resize(cnvtxs);
    cgraph->home .resize(cnvtxs);
  }

  if( !graph->where.empty() )
    cgraph->where.resize(cnvtxs);

  // Use upper bound estimate for these sizes
  std::vector<int> cadjncy(nkeepsize + nrecvsize);
  std::vector<int> cadjwgt(nkeepsize + nrecvsize);

  int htable[HASH_SIZE], htableidx[HASH_SIZE];
  iset(HASH_SIZE, -1, htable);

  cxadj[0] = cnvtxs = 0;
  int cnedges = 0;
  for(int i=0; i < nvtxs; ++i) {
    if( match[i] >= KEEP_FLAG ) {
      int v = first_vtx + i; 
      int u = match[i] - KEEP_FLAG;

      if( u >= first_vtx && u < last_vtx && v > u )
        continue;  // Pair (u,v) is already collapsed

      // Collapse the v vertex first, which you know is local
      for(int j=0; j < ncon; ++j)
        cvwgt[cnvtxs*ncon+j] = graph->vwgt[i*ncon+j];

      if( ctrl->is_adaptive ) {
        cgraph->vsize[cnvtxs] = graph->vsize[i];
        cgraph->home [cnvtxs] = graph->home [i];
      }

      if( !graph->where.empty() )
        cgraph->where[cnvtxs] = graph->where[i];

      int nedges = 0;
      // Collapse the vertices
      for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
        int k = cmap[graph->adjncy[j]];
        if( k < 0 ) {
	  throw std::runtime_error("Bad logic: negative index in CreateCoarseGraph!");
	}

        if( k != cfirst_vtx + cnvtxs ) {  // The edge is not internal
          int km = k & HASH_MASK;
          if( htable[km] == -1 ) { // Seeing this for first time
            htable[km] = k;
            htableidx[km] = cnedges + nedges;
            cadjncy[cnedges+nedges] = k; 
            cadjwgt[cnedges+nedges++] = graph->adjwgt[j];
          } else if( htable[km] == k ) {
            cadjwgt[htableidx[km]] += graph->adjwgt[j];
          } else { // Expensive case - need to search
	    int ie;
            for(ie=0; ie < nedges; ++ie) {
              if( cadjncy[cnedges+ie] == k )
                break;
            }
            if( ie < nedges ) {
              cadjwgt[cnedges+ie] += graph->adjwgt[j];
            } else {
              cadjncy[cnedges+nedges] = k; 
              cadjwgt[cnedges+nedges++] = graph->adjwgt[j];
            }
          }
        }
      }

      // Collapse the u vertex next
      if( v != u ) {
        if( u >= first_vtx && u < last_vtx) { // Local vertex
          u -= first_vtx;
          for(int j=0; j < ncon; ++j)
            cvwgt[cnvtxs*ncon+j] += graph->vwgt[u*ncon+j];

          if( ctrl->is_adaptive )
            cgraph->vsize[cnvtxs] += graph->vsize[u];

          for(int j=graph->xadj[u]; j < graph->xadj[u+1]; ++j) {
            int k = cmap[graph->adjncy[j]];
            if( k != cfirst_vtx + cnvtxs ) {  // If this is not an internal edge
              int km = k & HASH_MASK;
              if( htable[km] == -1 ) { // Seeing this for first time
                htable[km] = k;
                htableidx[km] = cnedges + nedges;
                cadjncy[cnedges+nedges] = k; 
                cadjwgt[cnedges+nedges] = graph->adjwgt[j];
                ++nedges;
              } else if( htable[km] == k ) {
                cadjwgt[htableidx[km]] += graph->adjwgt[j];
              } else { // Expensive case - need to search
		int ie;
                for(ie=0; ie < nedges; ++ie) {
                  if( cadjncy[cnedges+ie] == k )
                    break;
                }
                if( ie < nedges ) {
                  cadjwgt[cnedges+ie] += graph->adjwgt[j];
                } else {
                  cadjncy[cnedges+nedges] = k; 
                  cadjwgt[cnedges+nedges] = graph->adjwgt[j];
                  nedges++;
                }
              }
            }
          }
        } else { // Remote vertex
          u = perm[binary_search(graph->nrecv, recvind, u)];

          for(int j=0; j < ncon; ++j)
            // Remember that the +1 stores the vertex weight
            cvwgt[cnvtxs*ncon+j] += rgraph[(u+1)+j];

          if( ctrl->is_adaptive ) {
            cgraph->vsize[cnvtxs] += rgraph[u+1+ncon];
            cgraph->home [cnvtxs]  = rgraph[u+2+ncon];
          }

          for(int j=0; j < rgraph[u]; ++j) {
            int k = rgraph[u+3+ncon+2*j];
            if( k != cfirst_vtx+cnvtxs ) {  // If this is not an internal edge
              int km = k & HASH_MASK;
              if( htable[km] == -1 ) { // Seeing this for first time
                htable[km] = k;
                htableidx[km] = cnedges+nedges;
                cadjncy[cnedges+nedges] = k; 
                cadjwgt[cnedges+nedges] = rgraph[u+3+ncon+2*j+1];
                ++nedges;
              } else if( htable[km] == k ) {
                cadjwgt[htableidx[km]] += rgraph[u+3+ncon+2*j+1];
              } else { // Expensive case - need to search
		int ie;
                for(ie=0; ie < nedges; ++ie) {
                  if( cadjncy[cnedges+ie] == k )
                    break;
                }
                if( ie < nedges ) {
                  cadjwgt[cnedges+ie] += rgraph[u+3+ncon+2*j+1];
                } else {
                  cadjncy[cnedges+nedges] = k; 
                  cadjwgt[cnedges+nedges] = rgraph[u+3+ncon+2*j+1];
                  nedges++;
                }
              }
            }
          }
        }
      }

      cnedges += nedges;
      for(int j=cxadj[cnvtxs]; j < cnedges; ++j)
        htable[cadjncy[j] & HASH_MASK] = -1;  // reset the htable
      cxadj[++cnvtxs] = cnedges;
    }
  }

  cgraph->nedges = cnedges;

  for(int j=0; j < cnvtxs; ++j) {
    for(int i=0; i < ncon; ++i)
      cgraph->nvwgt[j*ncon+i] = ctrl->inv_weights[i] * cvwgt[j*ncon+i];
  }

  cgraph->adjncy.assign(cadjncy.data(), cadjncy.data() + cnedges);
  cgraph->adjwgt.assign(cadjwgt.data(), cadjwgt.data() + cnedges);

  //delete cgraph;

  graph->where.clear();
}


//------------------------------------------------------------------------
// Create the coarser graph after local matching
void create_cgraph_local(Control *ctrl, Graph *graph, int cnvtxs)
{
  int ncon = graph->ncon;  // convinience shortcut

  int first_vtx = graph->vtxdist[ctrl->myrank];

  graph->cmap.assign(graph->nvtxs+graph->nrecv, -1);
  int *cmap = graph->cmap.data();

  // Initialize the coarser graph
  Graph *cgraph = create_graph();
  cgraph->nvtxs  = cnvtxs;
  cgraph->ncon   = ncon;
  cgraph->level  = graph->level+1;
  cgraph->finer  = graph;
  graph->coarser = cgraph;

  // Obtain the vtxdist of the coarser graph
  cgraph->vtxdist.resize(ctrl->ncpus+1);
  int *cvtxdist = cgraph->vtxdist.data();

  cvtxdist[ctrl->ncpus] = cnvtxs;  // Use last position in the cvtxdist as a temp buffer

  MPI_Allgather(&cvtxdist[ctrl->ncpus], 1, MPI_INT, cvtxdist, 1, MPI_INT, ctrl->comm);

  for(int i=1; i < ctrl->ncpus; ++i) cvtxdist[i] += cvtxdist[i-1];
  for(int i=ctrl->ncpus; i > 0; --i) cvtxdist[i]  = cvtxdist[i-1];
  cvtxdist[0] = 0;

  cgraph->gnvtxs = cvtxdist[ctrl->ncpus];

  // Construct the cmap vector
  int cfirst_vtx = cvtxdist[ctrl->myrank];

  // Create the cmap of what you know so far locally
  cnvtxs = 0;
  for(int i=0; i < graph->nvtxs; ++i) {
    if( graph->match[i] >= KEEP_FLAG ) {
      int k = graph->match[i] - KEEP_FLAG;
      if( k < first_vtx + i )
        continue;  // i has been matched via the (k,i) side

      cmap[i] = cfirst_vtx + cnvtxs;
      ++cnvtxs;
      if( k != first_vtx + i ) {
        cmap[k - first_vtx] = cmap[i];
        graph->match[k - first_vtx] += KEEP_FLAG;
      }
    }
  }

  distribute_bdry_data(ctrl, graph, cmap, cmap + graph->nvtxs);

  // Create the coarser graph
  cgraph->xadj.resize(cnvtxs+1);
  int *cxadj = cgraph->xadj.data();

  cgraph->vwgt.resize(cnvtxs*ncon);
  int *cvwgt = cgraph->vwgt.data();

  cgraph->nvwgt.resize(cnvtxs*ncon);
  double *cnvwgt = cgraph->nvwgt.data();

  if( ctrl->is_adaptive )      cgraph->home.resize(cnvtxs);
  if( !graph->vsize.empty() )  cgraph->vsize.resize(cnvtxs);
  if( !graph->where.empty() )  cgraph->where.resize(cnvtxs);

  std::vector<int> cadjncy(graph->nedges);
  std::vector<int> cadjwgt(graph->nedges);

  int htable[HASH_SIZE], htableidx[HASH_SIZE];
  iset(HASH_SIZE, -1, htable);

  cxadj[0] = cnvtxs = 0;
  int cnedges = 0;

  for(int i=0; i < graph->nvtxs; ++i) {
    int v = first_vtx + i; 
    int u = graph->match[i] - KEEP_FLAG;

    if( v > u ) 
      continue;  // Pair (u,v) is already collapsed

    // Collapse the v vertex first (it should be local)
    for(int j=0; j < ncon; ++j)
      cvwgt[cnvtxs*ncon+j] = graph->vwgt[i*ncon+j];

    if( ctrl->is_adaptive )     cgraph->home [cnvtxs] = graph->home [i];
    if( !graph->vsize.empty() ) cgraph->vsize[cnvtxs] = graph->vsize[i];
    if( !graph->where.empty() ) cgraph->where[cnvtxs] = graph->where[i];

    int nedges = 0;
    for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
      int k = cmap[graph->adjncy[j]];
      if( k != cfirst_vtx + cnvtxs ) {  // If this is not an internal edge
        int km = k & HASH_MASK;
        if( htable[km] == -1 ) { // Seeing this for first time
          htable[km] = k;
          htableidx[km] = cnedges + nedges;
          cadjncy[cnedges+nedges] = k; 
          cadjwgt[cnedges+nedges] = graph->adjwgt[j];
	  ++nedges;
        } else if( htable[km] == k ) {
          cadjwgt[htableidx[km]] += graph->adjwgt[j];
        } else { // Expensive case - need to search
	  int ie;
          for(ie=0; ie<nedges; ++ie) {
            if( cadjncy[cnedges+ie] == k )
              break;
          }
          if( ie < nedges ) {
            cadjwgt[cnedges+ie] += graph->adjwgt[j];
          } else {
            cadjncy[cnedges+nedges] = k; 
            cadjwgt[cnedges+nedges] = graph->adjwgt[j];
	    ++nedges;
          }
        }
      }
    }

    // Collapse the u vertex next
    if( v != u ) { 
      u -= first_vtx;
      for(int j=0; j < ncon; ++j)
        cvwgt[cnvtxs*ncon+j] += graph->vwgt[u*ncon+j];

      if( !graph->vsize.empty() )
	cgraph->vsize[cnvtxs] += graph->vsize[u];

      if( !graph->where.empty() && (cgraph->where[cnvtxs] != graph->where[u]) ) {
	std::cout << "Something wrong with <where> local matching: "
		  << cgraph->where[cnvtxs] << ", " << graph->where[u]
		  << std::endl;
	throw std::runtime_error("Something wrong with <where> local matching");
      }

      for(int j=graph->xadj[u]; j < graph->xadj[u+1]; ++j) {
        int k = cmap[graph->adjncy[j]];
        if( k != cfirst_vtx + cnvtxs ) {  // If this is not an internal edge
          int km = k & HASH_MASK;
          if( htable[km] == -1 ) { // Seeing this for first time
            htable[km] = k;
            htableidx[km] = cnedges + nedges;
            cadjncy[cnedges+nedges] = k; 
            cadjwgt[cnedges+nedges] = graph->adjwgt[j];
	    ++nedges;
          } else if( htable[km] == k ) {
            cadjwgt[htableidx[km]] += graph->adjwgt[j];
          } else { // Expensive case - need to search
	    int ie;
            for(ie=0; ie < nedges; ++ie) {
              if( cadjncy[cnedges+ie] == k )
                break;
            }
            if( ie < nedges ) {
              cadjwgt[cnedges+ie] += graph->adjwgt[j];
            } else {
              cadjncy[cnedges+nedges] = k; 
              cadjwgt[cnedges+nedges] = graph->adjwgt[j];
	      ++nedges;
            }
          }
        }
      }
    }

    cnedges += nedges;
    for(int j=cxadj[cnvtxs]; j < cnedges; ++j)
      htable[cadjncy[j] & HASH_MASK] = -1;  // reset the htable
    cxadj[++cnvtxs] = cnedges;
  }

  cgraph->nedges = cnedges;

  for(int j=0; j < cnvtxs; ++j) {
    for(int i=0; i < ncon; ++i) 
      cgraph->nvwgt[j*ncon+i] = ctrl->inv_weights[i] * cvwgt[j*ncon+i];
  }

  cgraph->adjncy.assign(cadjncy.data(), cadjncy.data() + cnedges);
  cgraph->adjwgt.assign(cadjwgt.data(), cadjwgt.data() + cnedges);

  //delete cgraph;

  graph->where.clear();
}


//------------------------------------------------------------------------
// Find an HEM matching involving both local and remote vertices
void match_global(Control *ctrl, Graph *graph)
{
  static const int MATCH_STEPS = 4;

  double maxnvwgt = 0.75 / ctrl->coarse_size;

  graph->global_match = true;

  int     nvtxs   = graph->nvtxs;   // convinience shortcut
  int     ncon    = graph->ncon;    // convinience shortcut
  int     nnbrs   = graph->nnbrs;   // convinience shortcut
  int    *peind   = graph->peind.data();   // convinience shortcut
  int    *sendptr = graph->sendptr.data(); // convinience shortcut
  int    *recvptr = graph->recvptr.data(); // convinience shortcut
  double *nvwgt   = graph->nvwgt.data();   // convinience shortcut

  int first_vtx = graph->vtxdist[ctrl->myrank];

  graph->match.assign(nvtxs+graph->nrecv, UNMATCHED);
  int *match  = graph->match.data();

  // Local arrays
  std::vector<int> myhome(nvtxs+graph->nrecv, UNMATCHED);
  std::vector<int> nreqs_sr(nnbrs, 0);
  std::vector<int> perm(nvtxs);
  std::vector<int> perm_inv(nvtxs);
  std::vector<int> perm_new(nnbrs);
  std::vector<int> changed(nvtxs);

  std::vector<ikv_t> match_send(graph->nsend);
  std::vector<ikv_t> match_recv(graph->nrecv);

  // If coasening for adaptive partition, exchange home information
  if( ctrl->is_adaptive ) {
    ASSERT(!graph->home.empty());
    icopy(nvtxs, graph->home.data(), myhome.data());
    distribute_bdry_data(ctrl, graph, myhome.data(), &myhome[nvtxs]);
  }

  // Create traversal order
  random_permute_fast(nvtxs, perm.data(), 1);
  for(int i=0; i < nvtxs; ++i)
    perm_inv[perm[i]] = i;
  iincset(nnbrs, 0, perm_new.data());

  // Mark all heavy vertices so they will not be matched
  int nchanged = 0;
  for(int i=0; i < nvtxs; ++i) {
    for(int j=0; j < ncon; ++j) {
      if( nvwgt[i*ncon+j] > maxnvwgt ) {
        match[i] = TOO_HEAVY;
        ++nchanged;
        break;
      }
    }
  }

  // If no heavy vertices, mark one at random to have at least one vertex
  // left in each partition.
  if( nchanged == 0 )
    match[RandomInRange(nvtxs)] = TOO_HEAVY;

  distribute_bdry_data(ctrl, graph, match, match+nvtxs);

  // Set initial value of nkept
  int nkept = graph->gnvtxs / ctrl->ncpus - nvtxs;

  // Find a matching in several passes
  int nmatched = 0;
  for(int pass=0; pass < MATCH_STEPS; ++pass) {
    int nreqs = 0;
    int wside = (graph->level + pass)%2;
    nchanged = 0;
    for(int unmatched=nmatched, ii=nmatched; ii < nvtxs; ++ii) {
      int i = perm[ii];
      if( match[i] == UNMATCHED ) {
	int maxi = -1,  maxidx = i;

        if( graph->xadj[i] == graph->xadj[i+1] ) {
          unmatched = std::max(ii, unmatched) + 1;
          for(; unmatched < nvtxs; ++unmatched) {
            int k = perm[unmatched];
            if( match[k] == UNMATCHED && myhome[i] == myhome[k] ) {
              match[i] = first_vtx + k + (i <= k ? KEEP_FLAG : 0);
              match[k] = first_vtx + i + (i >  k ? KEEP_FLAG : 0);
              changed[nchanged++] = i;
              changed[nchanged++] = k;
              break;
            }
          }
          continue;
        }

        // Find heavy-edge matching
        for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
          int k = graph->adjncy[j];
          if( match[k] == UNMATCHED && myhome[k] == myhome[i] ) { 
            if( ncon == 1 ) {
              if( maxi == -1 || graph->adjwgt[maxi] < graph->adjwgt[j] ||
                  (graph->adjwgt[maxi] == graph->adjwgt[j] &&
		   RandomInRange(graph->xadj[i+1]-graph->xadj[i]) == 0) ) {
                maxi   = j;
                maxidx = k;
              }
            } else {
              if( maxi == -1 || graph->adjwgt[maxi] < graph->adjwgt[j] ||
                  (graph->adjwgt[maxi] == graph->adjwgt[j] && maxidx < nvtxs && k < nvtxs &&
                   is_wbalance_better(ncon, nvwgt+i*ncon, nvwgt+maxidx*ncon, nvwgt+k*ncon)) ) {
                maxi   = j;
                maxidx = k;
              }
            }
          }
        }

        if( maxi != -1 ) {
          int k = graph->adjncy[maxi];
          if( k < nvtxs ) { // Take care the local vertices first
            // Give preference the local matching
            match[i] = first_vtx + k + (i <= k ? KEEP_FLAG : 0);
            match[k] = first_vtx + i + (i >  k ? KEEP_FLAG : 0);
            changed[nchanged++] = i;
            changed[nchanged++] = k;
          } else { // Take care any remote boundary vertices
            match[k] = MAYBE_MATCHED;
            // Alternate among which vertices will issue the requests
            if( (wside == 0 && first_vtx+i < graph->imap[k]) || 
                (wside == 1 && first_vtx+i > graph->imap[k])) { 
              match[i] = MAYBE_MATCHED;
              match_send[nreqs].key = graph->imap[k];
              match_send[nreqs].val = first_vtx+i;
              ++nreqs;
            }
          }
        }
      }
    }

    // Requests for this cpu are stored in match_recv
    for(int i=0; i < nnbrs; ++i) {
      MPI_Irecv(&match_recv[recvptr[i]], 2*(recvptr[i+1]-recvptr[i]), MPI_INT,
                peind[i], 1, ctrl->comm, &ctrl->rreqs[i]);
    }

    // Issue the sends next
    ikvsorti(nreqs, match_send.data());
    for(int j=0, i=0; i < nnbrs; ++i) {
      int k, other_last_vtx = graph->vtxdist[peind[i]+1];
      for(k=j; k < nreqs && match_send[k].key < other_last_vtx; ++k);
      MPI_Isend(&match_send[j], 2*(k-j), MPI_INT, peind[i], 1, 
		ctrl->comm, &ctrl->sreqs[i]);
      j = k;
    }

    // Wait for the operations to finish
    MPI_Waitall(nnbrs, &ctrl->rreqs[0], &ctrl->wstats[0]);
    for(int i=0; i < nnbrs; ++i) {
      MPI_Get_count(&ctrl->wstats[i], MPI_INT, &nreqs_sr[i]);
      nreqs_sr[i] = nreqs_sr[i] / 2;  // Adjust for pairs
    }
    MPI_Waitall(nnbrs, &ctrl->sreqs[0], &ctrl->wstats[0]);

    // Process the requests that were received in match_recv
    random_permute(nnbrs, perm_new.data(), 0);
    for(int ii=0; ii < nnbrs; ++ii) {
      int i = perm_new[ii];
      ikv_t *recv_reqs = &match_recv[recvptr[i]];
      for(int j=0; j < nreqs_sr[i]; ++j) {
        int k = recv_reqs[j].key;

        if( match[k - first_vtx] == UNMATCHED ) { // Need to grant this request
          changed[nchanged++] = k - first_vtx;
          if( nkept >= 0 ) { // Deciding which to keep based on local balance
            match[k - first_vtx] = recv_reqs[j].val + KEEP_FLAG;
            --nkept;
          } else {
            match[k - first_vtx] = recv_reqs[j].val;
            recv_reqs[j].key += KEEP_FLAG;
            ++nkept;
          }
        } else {
          recv_reqs[j].key = UNMATCHED;
        }
      }
    }

    // Exchange match_recv information. It is stored in match_send
    for(int i=0; i < nnbrs; ++i) {
      MPI_Irecv(&match_send[sendptr[i]], 2*(sendptr[i+1]-sendptr[i]), MPI_INT,
                peind[i], 1, ctrl->comm, &ctrl->rreqs[i]);
    }

    // Issue the sends next
    for(int i=0; i < nnbrs; ++i) {
      MPI_Isend(&match_recv[recvptr[i]], 2*nreqs_sr[i], MPI_INT, 
                peind[i], 1, ctrl->comm, &ctrl->sreqs[i]);
    }

    // Wait for the operations to finish
    MPI_Waitall(nnbrs, &ctrl->rreqs[0], &ctrl->wstats[0]);
    for(int i=0; i < nnbrs; ++i) {
      MPI_Get_count(&ctrl->wstats[i], MPI_INT, &nreqs_sr[i]);
      nreqs_sr[i] = nreqs_sr[i]/2;  // Adjust for pairs
    }
    MPI_Waitall(nnbrs, &ctrl->sreqs[0], &ctrl->wstats[0]);

    // Go through all match_send and update local match information
    // for the matchings that were granted
    for(int i=0; i < nnbrs; ++i) {
      ikv_t *send_reqs = &match_send[sendptr[i]];
      for(int j=0; j < nreqs_sr[i]; ++j) {
        match[send_reqs[j].val - first_vtx] = send_reqs[j].key;
        if( send_reqs[j].key != UNMATCHED )
          changed[nchanged++] = send_reqs[j].val-first_vtx;
      }
    }

    for(int i=0; i < nchanged; ++i) {
      int ii = perm_inv[changed[i]];
      perm[ii] = perm[nmatched];
      perm_inv[perm[nmatched]] = ii;
      ++nmatched;
    }

    exchange_bdry_data(ctrl, graph, nchanged, changed.data(), match, 
		       match_send.data(), match_recv.data());
  }

  // Unmatched vertices will be matched with themselves
  int cnvtxs = 0;
  for(int i=0; i < nvtxs; ++i) {
    if( match[i] == UNMATCHED || match[i] == TOO_HEAVY ) {
      match[i] = (first_vtx + i) + KEEP_FLAG;
      ++cnvtxs;
    } else if (match[i] >= KEEP_FLAG) {  // A matched vertex which is kept
      ++cnvtxs;
    }
  }

   create_cgraph_global(ctrl, graph, cnvtxs);
}


//------------------------------------------------------------------------
// Find an HEM matching involving only local vertices
void match_local(Control *ctrl, Graph *graph)
{
  double maxnvwgt = 0.75 / ctrl->coarse_size;

  graph->global_match = false;

  int     nvtxs = graph->nvtxs; // convinience shortcut
  int     ncon  = graph->ncon;  // convinience shortcut
  double *nvwgt = graph->nvwgt.data(); // convinience shortcut

  int first_vtx = graph->vtxdist[ctrl->myrank];

  graph->match.resize(nvtxs+graph->nrecv);
  int *match  = graph->match.data();

  // Local arrays
  std::vector<int> myhome(nvtxs+graph->nrecv, UNMATCHED);
  std::vector<int> perm(nvtxs);

  // If coasening for adaptive partition, exchange home information
  if( ctrl->is_adaptive ) {
    icopy(nvtxs, graph->home.data(), myhome.data());
    distribute_bdry_data(ctrl, graph, myhome.data(), &myhome[nvtxs]);
  }

  //  Find a local matching
  iset(nvtxs, UNMATCHED, match);
  iset(graph->nrecv, 0, match+nvtxs);

  random_permute_fast(nvtxs, perm.data(), 1);
  int cnvtxs = 0;
  for(int ii=0; ii < nvtxs; ++ii) {
    int i = perm[ii];
    if( match[i] == UNMATCHED ) {
      int maxidx = -1, maxi = -1;

      // Find a heavy-edge matching, if the weight of the vertex is OK
      int j;
      for(j=0; j < ncon; ++j)
        if( nvwgt[i*ncon+j] > maxnvwgt )
          break;

      if( j == ncon && ii < nvtxs ) {
        for(int jj=graph->xadj[i]; jj < graph->xadj[i+1]; ++jj) {
          int edge = graph->adjncy[jj];

          // Match only with local vertices
          if( myhome[edge] != myhome[i] || edge >= nvtxs )
            continue;

          for(j=0; j < ncon; ++j)
            if( nvwgt[edge*ncon+j] > maxnvwgt ) 
              break;

          if( j == ncon ) {
            if( match[edge] == UNMATCHED &&
		(maxi == -1 ||
		 graph->adjwgt[maxi] < graph->adjwgt[jj] ||
		 (graph->adjwgt[maxi] == graph->adjwgt[jj] &&
		  is_wbalance_better(ncon, nvwgt+i*ncon, nvwgt+maxidx*ncon, nvwgt+edge*ncon))) ) {
              maxi = jj;
              maxidx = edge;
            }
          }
        }
      }

      if( maxi != -1 ) {
        int k = graph->adjncy[maxi];
        match[i] = first_vtx + k + (i <= k ? KEEP_FLAG : 0);
        match[k] = first_vtx + i + (i >  k ? KEEP_FLAG : 0);
      } else {
        match[i] = first_vtx + i + KEEP_FLAG;
      }
      ++cnvtxs;
    }
  }

  distribute_bdry_data(ctrl, graph, match, match+nvtxs);

  create_cgraph_local(ctrl, graph, cnvtxs);
}


//------------------------------------------------------------------------
// Project partition
void project_partition(Control *ctrl, Graph *graph)
{
  Graph *cgraph = graph->coarser; // convinience shortcut

  graph->where.resize(graph->nvtxs+graph->nrecv);

  int cfirst_vtx = cgraph->vtxdist[ctrl->myrank];
  int first_vtx  = graph->vtxdist[ctrl->myrank];

  std::vector<ikv_t> scand;
  if( graph->global_match ) {
    //---------------------------------------------------
    // Start exchanging the remote <where> information 
    ikv_t *rcand = graph->rcand.data(); // convinience shortcut

    scand.resize(graph->slens[graph->nnbrs]);

    // Always issue the receives first
    for(int i=0; i < graph->nnbrs; ++i) {
      int nitems = graph->slens[i+1] - graph->slens[i];
      if( nitems > 0 )
        MPI_Irecv(&scand[graph->slens[i]], 2*nitems, MPI_INT, 
		  graph->peind[i], 1, ctrl->comm, &ctrl->rreqs[i]);
    }

    // Put where[rcand[].key] into the val field of the candidate
    for(int i=0; i < graph->rlens[graph->nnbrs]; ++i) {
      ASSERT((rcand[i].val >= 0) && (rcand[i].val < cgraph->nvtxs));
      rcand[i].val = cgraph->where[rcand[i].val];
    }

    // Issue the sends after receives
    for(int i=0; i < graph->nnbrs; ++i) {
      int nitems = graph->rlens[i+1]-graph->rlens[i];
      if( nitems > 0 )
        MPI_Isend(&rcand[graph->rlens[i]], 2*nitems, MPI_INT, 
		  graph->peind[i], 1, ctrl->comm, &ctrl->sreqs[i]);
    }
  }

  //------------------------------------------------------------
  // Project local vertices
  for(int i=0; i < graph->nvtxs; ++i) {
    if( graph->match[i] >= KEEP_FLAG ) {
      ASSERT((graph->cmap[i] - cfirst_vtx >= 0) &&
             (graph->cmap[i] - cfirst_vtx < cgraph->nvtxs));
      graph->where[i] = cgraph->where[graph->cmap[i] - cfirst_vtx];
    }
  }

  if( graph->global_match ) {
    //---------------------------------------------------
    // Wait for nonblocking operations to finish
    for(int i=0; i < graph->nnbrs; ++i) {
      if( graph->rlens[i+1] - graph->rlens[i] > 0 )  
        MPI_Wait(&ctrl->sreqs[i], &ctrl->status);
    }
    for(int i=0; i < graph->nnbrs; ++i) {
      if( graph->slens[i+1] - graph->slens[i] > 0 )  
        MPI_Wait(&ctrl->rreqs[i], &ctrl->status);
    }

    //---------------------------------------------------
    // Project received vertices
    for(int i=0; i < graph->slens[graph->nnbrs]; ++i) {
      graph->where[scand[i].key - first_vtx] = scand[i].val;
    }
  }

  free_graph(graph->coarser);
  graph->coarser = NULL;
}


//------------------------------------------------------------------------
// Compute the initial id/ed
void get_part_params(Control *ctrl, Graph *graph)
{
  int nvtxs = graph->nvtxs; // convinience shortcut
  int ncon  = graph->ncon;  // convinience shortcut

  graph->ckrinfo.resize(nvtxs);
  std::memset(graph->ckrinfo.data(), 0, sizeof(ckrinfo_t)*nvtxs);

  graph->lnpwgts.assign(ctrl->nparts*ncon, 0.0);
  graph->gnpwgts.resize(ctrl->nparts*ncon);

  // Exchange information (where) of interface vertices
  distribute_bdry_data(ctrl, graph, &graph->where[0], &graph->where[nvtxs]); 

  // Compute now the id/ed degrees
  graph->lmincut = 0;
  for(int i=0; i < nvtxs; ++i) {
    int me = graph->where[i];
    ckrinfo_t *myrinfo = &graph->ckrinfo[i];

    for(int j=0; j < ncon; ++j)
      graph->lnpwgts[me*ncon+j] += graph->nvwgt[i*ncon+j];

    for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
      if( me == graph->where[graph->adjncy[j]] )
        myrinfo->id += graph->adjwgt[j];
      else
        myrinfo->ed += graph->adjwgt[j];
    }

    if( myrinfo->ed > 0 ) {
      graph->lmincut += myrinfo->ed;

      myrinfo->inbr = ctrl->pool.get_next(graph->xadj[i+1] - graph->xadj[i]);
      std::pair<int,int> *mynbrs = ctrl->pool.data() + myrinfo->inbr;

      for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
        int other = graph->where[graph->adjncy[j]];
        if( me != other ) {
	  bool exists = false;
          for(int k=0; k < myrinfo->nnbrs; ++k) {
            if( mynbrs[k].first == other ) {
              mynbrs[k].second += graph->adjwgt[j];
	      exists = true;
              break;
            }
          }
          if( !exists ) {
            mynbrs[myrinfo->nnbrs].first  = other;
            mynbrs[myrinfo->nnbrs].second = graph->adjwgt[j];
            myrinfo->nnbrs++;
          }
          ASSERT(myrinfo->nnbrs <= graph->xadj[i+1]-graph->xadj[i]);
        }
      }
    } else {
      myrinfo->inbr = -1;
    }
  }

  // Sum partition weights
  MPI_Allreduce(graph->lnpwgts.data(), graph->gnpwgts.data(), ctrl->nparts*ncon, 
		MPI_DOUBLE, MPI_SUM, ctrl->comm);

  MPI_Allreduce(&graph->lmincut, &graph->mincut, 1, MPI_INT, MPI_SUM, ctrl->comm);
  graph->mincut /= 2;
}


//------------------------------------------------------------------------
// Perform refinement
void balance_kw(Control *ctrl, Graph *graph, int npasses)
{
  int  nvtxs   = graph->nvtxs;   // convinience shortcut
  int  ncon    = graph->ncon;    // convinience shortcut
  int *peind   = graph->peind.data();  // convinience shortcut
  int *sendptr = graph->sendptr.data(); // convinience shortcut
  int  nparts  = ctrl->nparts;   // convinience shortcut

  int first_vtx = graph->vtxdist[ctrl->myrank];

  // Allocate working space
  std::vector<int> nupd_cpu (ctrl->ncpus);
  std::vector<int> pperm    (nparts);
  std::vector<int> old_ed   (nvtxs);
  std::vector<int> changed  (nvtxs);
  std::vector<int> perm     (nvtxs);
  std::vector<int> update   (nvtxs);
  std::vector<int> moved    (nvtxs);
  std::vector<int> htable   (nvtxs+graph->nrecv, 0);
  std::vector<int> where_tmp(nvtxs+graph->nrecv);
  std::vector<int> supdate  (graph->nrecv);
  std::vector<int> rupdate  (graph->nsend);

  std::vector<ikv_t> rwchanges(graph->nrecv);
  std::vector<ikv_t> swchanges(graph->nsend);

  icopy(nvtxs+graph->nrecv, graph->where.data(), where_tmp.data());

  // Store external degrees of the vertices before refinement
  for(int i=0; i < nvtxs; ++i)
    old_ed[i] = graph->ckrinfo[i].ed;

  //-------------------------------------------------------------
  // Perform several walk ins through the vertices
  for(int nswaps=0, pass=0; pass < npasses; ++pass) {
    int oldcut = graph->mincut;
    if( ctrl->myrank == 0 )
      random_permute(nparts, pperm.data(), 1);
    MPI_Bcast(pperm.data(), nparts, MPI_INT, 0, ctrl->comm);
    random_permute_fast(nvtxs, perm.data(), 1);

    for(int step=0; step < 2; ++step) {
      int nmoved = 0;

      // Step 1 - record staistics for desired moves
      int to = -1;
      for(int ii=0; ii < nvtxs; ++ii) {
        int i = perm[ii];
        int from = where_tmp[i];
        double *nvwgt = &graph->nvwgt[i*ncon];

	bool weight_coincide = false;
        for(int jj=0; jj < ncon; ++jj) {
          if( std::abs(nvwgt[jj] - graph->gnpwgts[from*ncon+jj]) < SMALL_FLOAT ) {
	    weight_coincide = true;
            break;
	  }
        }

        if( weight_coincide )
          continue;

        // Check for possible improvement
        ckrinfo_t *myrinfo = &graph->ckrinfo[i];
        if( myrinfo->ed == 0 || myrinfo->ed < myrinfo->id )
          continue;

        ASSERT(myrinfo->inbr != -1);
        std::pair<int,int> *mynbrs = ctrl->pool.data() + myrinfo->inbr;

	int candidate_indx = -1;
        for(int k=myrinfo->nnbrs-1; k >= 0; --k) {
          to = mynbrs[k].first;
	  bool proper_side = ((step == 0) ? (pperm[from] - pperm[to] < 0)
                                          : (pperm[from] - pperm[to] > 0));
	  if( proper_side &&
              check_better_balance_wt(ncon,
                                      &graph->gnpwgts[from*ncon],
                                      &graph->gnpwgts[to*ncon],
                                      nvwgt, ctrl->disbalance.data()) ) {
	    candidate_indx = k;
            break;
          }
        }

        // Skip the rest if no candidate found
        if( candidate_indx < 0 )
          continue;

        int oldto = to;
        for(int j=candidate_indx-1; j >= 0; --j) {
          to = mynbrs[j].first;
	  bool proper_side = ((step == 0) ? (pperm[from] - pperm[to] < 0)
                                          : (pperm[from] - pperm[to] > 0));
	  if( proper_side &&
              check_better_balance_pt(ncon,
                                      &graph->gnpwgts[oldto*ncon],
                                      &graph->gnpwgts[to*ncon],
                                      nvwgt, ctrl->disbalance.data()) ) {
            oldto = to;
            candidate_indx = j;
          }
        }
        to = oldto;

        if( ii % ctrl->ncpus == 0 ) {
          // Update temp arrays of the moved vertex 
          where_tmp[i] = to;
          moved[nmoved++] = i;
          for(int j=0; j < ncon; ++j) {
	    graph->lnpwgts[to*ncon+j]   += nvwgt[j];
	    graph->lnpwgts[from*ncon+j] -= nvwgt[j];
	    graph->gnpwgts[to*ncon+j]   += nvwgt[j];
	    graph->gnpwgts[from*ncon+j] -= nvwgt[j];
          }

          myrinfo->ed += myrinfo->id - mynbrs[candidate_indx].second;
          std::swap(myrinfo->id, mynbrs[candidate_indx].second);
          if( mynbrs[candidate_indx].second == 0 )
            mynbrs[candidate_indx] = mynbrs[--myrinfo->nnbrs];
          else
            mynbrs[candidate_indx].first = from;

          // Update the degrees of adjacent vertices
          for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
            // Skip verts from other cpus
            if( graph->adjncy[j] >= nvtxs )
              continue;

            int me = graph->adjncy[j];
            int this_part = where_tmp[me];

            myrinfo = &graph->ckrinfo[me];
            if( myrinfo->inbr == -1 ) {
              myrinfo->inbr  = ctrl->pool.get_next(graph->xadj[me+1] - graph->xadj[me]);
              myrinfo->nnbrs = 0;
            }
            mynbrs = ctrl->pool.data() + myrinfo->inbr;

            if( this_part == from ) {
	      myrinfo->ed += graph->adjwgt[j];
	      myrinfo->id -= graph->adjwgt[j];
            } else {
              if( this_part == to ) {
		myrinfo->id += graph->adjwgt[j];
		myrinfo->ed -= graph->adjwgt[j];
              }
            }

            // Remove contribution from .ed of donor
            if( this_part != from ) {
              for(int k=0; k < myrinfo->nnbrs; ++k) {
                if( mynbrs[k].first == from ) {
                  if( mynbrs[k].second == graph->adjwgt[j] )
                    mynbrs[k] = mynbrs[--myrinfo->nnbrs];
                  else
                    mynbrs[k].second -= graph->adjwgt[j];
                  break;
                }
              }
            }

            // Add contribution to <ed> of recepient
            if( this_part != to ) {
	      bool exist = false;
              for(int k=0; k < myrinfo->nnbrs; ++k) {
                if( mynbrs[k].first == to ) {
                  mynbrs[k].second += graph->adjwgt[j];
		  exist = true;
                  break;
                }
              }
	      if( !exist ) {
                mynbrs[myrinfo->nnbrs].first = to;
                mynbrs[myrinfo->nnbrs].second= graph->adjwgt[j];
                myrinfo->nnbrs++;
              }
            }
          }
        }
      }

      // Step 2 - store the rest info for the moves
      int nlupd=0, nsupd=0, nmoves=0, nchanged=0;
      for(int ii=0; ii < nmoved; ++ii) {
        int i = moved[ii];
        if( i == -1 )
          continue;

        graph->where[i] = where_tmp[i];

        // Update the vertex information
        if( htable[i] == 0 ) {
          htable[i] = 1;
          update[nlupd++] = i;
        }

        // Put the vertices adjacent to i into the update array
        for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
          int k = graph->adjncy[j];
          if( htable[k] == 0 ) {
            htable[k] = 1;
            if( k < nvtxs )
              update[nlupd++] = k;
            else
              supdate[nsupd++] = k;
          }
        }
        ++nmoves;
        ++nswaps;

        if( graph->pexadj[i+1] - graph->pexadj[i] > 0 )
          changed[nchanged++] = i;
      }

      // Propagate to neighboring cpus the new where[] info for the interface vertices
      exchange_bdry_data(ctrl, graph, nchanged,
			 changed.data(), graph->where.data(),
			 swchanges.data(), rwchanges.data()); 

      //-------------------------------------------------------------
      // Send to other cpus the verts which degrees are to be updated
      // Issue the receives first
      for(int i=0; i < graph->nnbrs; ++i) {
        MPI_Irecv(&rupdate[sendptr[i]], sendptr[i+1]-sendptr[i], MPI_INT,
		  peind[i], 1, ctrl->comm, &ctrl->rreqs[i]);
      }

      // Issue the sends next
      for(int i=0; i < nsupd; ++i) {
        htable[supdate[i]] = 0;
        supdate[i] = graph->imap[supdate[i]];
      }
      isorti(nsupd, supdate.data());

      for(int k=0, j=0, i=0; i < graph->nnbrs; ++i) {
        int mylast_vtx = graph->vtxdist[peind[i]+1];
        for(k=j; (k < nsupd) && (supdate[k] < mylast_vtx); ++k); 
        MPI_Isend(&supdate[j], k-j, MPI_INT, peind[i], 1, ctrl->comm, &ctrl->sreqs[i]);
        j = k;
      }

      // Wait for the operations to finish
      MPI_Waitall(graph->nnbrs, &ctrl->rreqs[0], &ctrl->wstats[0]);
      for(int i=0; i < graph->nnbrs; ++i) 
        MPI_Get_count(&ctrl->wstats[i], MPI_INT, &nupd_cpu[i]);
      MPI_Waitall(graph->nnbrs, &ctrl->sreqs[0], &ctrl->wstats[0]);

      //-------------------------------------------------------------
      // Place the recieved marked vertices into update[] array
      for(int i=0; i < graph->nnbrs; ++i) {
        int *pe_updates = &rupdate[sendptr[i]];
        for(int j=0; j < nupd_cpu[i]; ++j) {
          int k = pe_updates[j];
          if( htable[k - first_vtx] == 0 ) {
            htable[k - first_vtx] = 1;
            update[nlupd] = k - first_vtx;
	    ++nlupd;
          }
        }
      }

      //-------------------------------------------------------------
      // Update information of the vertices in the update[] array
      for(int ii=0; ii < nlupd; ++ii) {
        int i = update[ii];
        ASSERT(htable[i] == 1);

        htable[i] = 0;

        int this_part = graph->where[i];
        ckrinfo_t *myrinfo = &graph->ckrinfo[i];

        if( myrinfo->inbr == -1 )
          myrinfo->inbr  = ctrl->pool.get_next(graph->xadj[i+1] - graph->xadj[i]);
        std::pair<int,int> *mynbrs = ctrl->pool.data() + myrinfo->inbr;

        graph->lmincut -= old_ed[i];
        myrinfo->nnbrs  = 0;
        myrinfo->id     = 0;
        myrinfo->ed     = 0;

        for(int j=graph->xadj[i]; j < graph->xadj[i+1]; ++j) {
          int other_part = graph->where[graph->adjncy[j]];
          if( this_part != other_part ) {
            myrinfo->ed += graph->adjwgt[j];

	    bool exist = false;
            for(int k=0; k < myrinfo->nnbrs; ++k) {
              if( mynbrs[k].first == other_part ) {
                mynbrs[k].second += graph->adjwgt[j];
		exist = true;
                break;
              }
            }
            if( !exist ) {
              mynbrs[myrinfo->nnbrs].first = other_part;
              mynbrs[myrinfo->nnbrs].second= graph->adjwgt[j];
              myrinfo->nnbrs++;
            }
            ASSERT(myrinfo->nnbrs <= graph->xadj[i+1]-graph->xadj[i]);
          } else {
            myrinfo->id += graph->adjwgt[j];
          }
        }
        graph->lmincut += myrinfo->ed;
        old_ed[i]       = myrinfo->ed; // for the next iteration
      }

      // Get sump of the partition weights
      MPI_Allreduce(graph->lnpwgts.data(), graph->gnpwgts.data(), nparts*ncon,
		    MPI_DOUBLE, MPI_SUM, ctrl->comm);
    }

    MPI_Allreduce(&graph->lmincut, &graph->mincut, 1, MPI_INT, MPI_SUM, ctrl->comm);
    graph->mincut /= 2;

    if( graph->mincut == oldcut )
      break;
  }
}


//------------------------------------------------------------------------
// The main driver for the adaptive refinement algorithm
void adaptive_partition(Control *ctrl, Graph *graph)
{
  // Set up data structures
  setup_comm(ctrl, graph);

  std::vector<double> lbvec(graph->ncon);

  double ubavg  = average(graph->ncon,   ctrl->disbalance.data());
  double tewgt  = isum(graph->nedges, graph->adjwgt.data(), 1);
  double tvsize = isum(graph->nvtxs,  graph->vsize.data(),  1);

  double small_add = 1. / graph->gnvtxs; // to avoid exceptions
  double gtewgt;
  MPI_Allreduce(&tewgt, &gtewgt, 1, MPI_DOUBLE, MPI_SUM, ctrl->comm);
  gtewgt += small_add;

  double gtvsize;
  MPI_Allreduce(&tvsize, &gtvsize, 1, MPI_DOUBLE, MPI_SUM, ctrl->comm);
  gtvsize += small_add;

  ctrl->redist_ratio = (gtewgt/gtvsize) / ctrl->edge_node_ratio;

  if( (graph->gnvtxs < (int)(1.3*ctrl->coarse_size)) ||
      (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_RATIO) ) {
    ctrl->pool.reserve(2*graph->nedges);

    //---------------------------------------------------
    // Balance the partition on the coarsest graph
    graph->where.resize(graph->nvtxs + graph->nrecv);
    icopy(graph->nvtxs, graph->home.data(), graph->where.data());

    get_balance_global(ctrl, graph, graph->where.data(), lbvec.data());
    double lbavg = average(graph->ncon, lbvec.data());

    if( lbavg > ubavg + AVG_THRESHOLD )
      balance_partition(ctrl, graph);

    // Check if no coarsening took place
    if( graph->finer == NULL ) {
      get_part_params(ctrl, graph);
      balance_kw(ctrl, graph, graph->ncon);
      refine_adaptive_kw(ctrl, graph, NGR_STEPS);
    }
  } else {
    // Coarsen and partition
    if( ctrl->is_coupled ) match_local (ctrl, graph);
    else                   match_global(ctrl, graph);

    adaptive_partition(ctrl, graph->coarser);

    //---------------------------------------------------
    // Project partition and refine
    project_partition(ctrl, graph);
    get_part_params(ctrl, graph);

    if( graph->ncon > 1 && graph->level < 4 ) {
      get_balance_global(ctrl, graph, graph->where.data(), lbvec.data());
      double lbavg = average(graph->ncon, lbvec.data());

      if( lbavg > ubavg + 0.025 ) {
        balance_kw(ctrl, graph, graph->ncon);
      }
    }

    refine_adaptive_kw(ctrl, graph, NGR_STEPS);
  }
}


//------------------------------------------------------------------------
// The entry point of the parallel multilevel local diffusion algorithm.
int AdaptiveRepartitioner(int *vtxdist, int *xadj,  int *adjncy,
                          int *vwgt,    int *vsize, int *adjwgt,
                          double *part_weights, double *disbalance,
                          int wgtflag, int ncon, int nparts,
                          double itr_ratio, int *options,
                          int *part, MPI_Comm *comm)
{
  // Check the input parameters and return in case of an error
  int lstatus = check_inputs(vtxdist, xadj, adjncy, vwgt, vsize, adjwgt, wgtflag, ncon, nparts,
                             part_weights, disbalance, itr_ratio, options, part, comm);
  int gstatus;
  MPI_Allreduce(&lstatus, &gstatus, 1, MPI_INT, MPI_MIN, *comm);
  if( gstatus == 0 )
    return -1; // error

  // Set up the control
  Control *ctrl = setup_controls(true, ncon, nparts, options, part_weights, disbalance, *comm);

  // Take care the nparts == 1 case
  if( nparts == 1 ) {
    iset(vtxdist[ctrl->myrank+1]-vtxdist[ctrl->myrank], 0, part); 

  } else {

    Graph *graph = setup_graph(ctrl, ncon, vtxdist, xadj, vwgt, vsize, adjncy, adjwgt, wgtflag);
  
    if( ctrl->is_coupled ) {
      iset(graph->nvtxs, ctrl->myrank, graph->home.data());
    } else {
      // Downgrade the partition numbers if part[] has more partitions that nparts
      for(int i=0; i < graph->nvtxs; ++i)
        if( part[i] >= ctrl->nparts ) part[i] = 0;

      icopy(graph->nvtxs, part, graph->home.data());
    }

    // Partition and Remap
    ctrl->itr_ratio   = itr_ratio;
    ctrl->coarse_size = std::min(graph->gnvtxs+1,
          (std::max(ctrl->ncpus, nparts) > 256 ? 20 : 50)*ncon*std::max(ctrl->ncpus, nparts));

    adaptive_partition(ctrl, graph);
    graph->remap(*ctrl);

    icopy(graph->nvtxs, graph->where.data(), part);
    //int edgecut = graph->mincut; // May report number of cuts if necessary

    graph->renumber_conns(adjncy);
    free_graph(graph);
  }

  free_controls(&ctrl);

  return 1;
}
