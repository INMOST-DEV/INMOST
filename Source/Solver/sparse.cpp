#define _CRT_SECURE_NO_WARNINGS
#include "inmost_sparse.h"
#include <fstream>
#include <sstream>

namespace INMOST
{
  namespace Sparse
  {
    const std::string stubstring = "";

    bool _hasRowEntryType = false;
    INMOST_MPI_Type RowEntryType = INMOST_MPI_DATATYPE_NULL;


    INMOST_MPI_Type GetRowEntryType() {return RowEntryType;}

    bool HaveRowEntryType() {return _hasRowEntryType;}

    void CreateRowEntryType()
    {
#if defined(USE_MPI)
      if( !HaveRowEntryType() )
      {
        int ierr;
        MPI_Datatype type[3] = { INMOST_MPI_DATA_ENUM_TYPE, INMOST_MPI_DATA_REAL_TYPE, MPI_UB};
		    int blocklen[3] = { 1, 1, 1 };
		    MPI_Aint disp[3];
		    disp[0] = offsetof(Sparse::Row::entry,first);
		    disp[1] = offsetof(Sparse::Row::entry,second);
		    disp[2] = sizeof(Sparse::Row::entry);
		    ierr = MPI_Type_create_struct(3, blocklen, disp, type, &RowEntryType);
		    if( ierr != MPI_SUCCESS )
		    {
			    std::cout << __FILE__ << ":" << __LINE__ << "problem in MPI_Type_create_struct" << std::endl;
		    }
		    ierr = MPI_Type_commit(&RowEntryType);
		    if( ierr != MPI_SUCCESS )
		    {
			    std::cout << __FILE__ << ":" << __LINE__ << "problem in MPI_Type_commit" << std::endl;
		    }
      }
#endif
      _hasRowEntryType = true;
    }

    void DestroyRowEntryType()
    {
#if defined(USE_MPI)
      if( HaveRowEntryType() )
      {
		    MPI_Type_free(&RowEntryType);
        RowEntryType = INMOST_MPI_DATATYPE_NULL;
      }
#endif
      _hasRowEntryType = false;
    }

    RowMerger::RowMerger() : Sorted(true), Nonzeros(0) {}

    INMOST_DATA_REAL_TYPE & RowMerger::operator[] (INMOST_DATA_ENUM_TYPE pos)
    {
      if( LinkedList[pos+1].first != UNDEF ) return LinkedList[pos+1].second;
      else 
      {
        INMOST_DATA_ENUM_TYPE index = LinkedList.get_interval_beg(), next;
        if( Sorted )
        {
          next = index;
          while(next < pos+1)
          {
            index = next;
            next = LinkedList[index].first;
          }
          assert(index < pos+1);
          assert(pos+1 < next);
          ++Nonzeros;
          LinkedList[index].first = pos+1;
          LinkedList[pos+1].first = next;
          return LinkedList[pos+1].second;
        }
        else
        {
          INMOST_DATA_ENUM_TYPE index = LinkedList.get_interval_beg();
          ++Nonzeros;
          LinkedList[pos+1].first = LinkedList[index].first;
          LinkedList[index].first = pos+1;
          return LinkedList[pos+1].second;
        }
      }
    }

    INMOST_DATA_REAL_TYPE RowMerger::operator[] (INMOST_DATA_ENUM_TYPE pos) const
    {
      if( LinkedList[pos+1].first != UNDEF ) return LinkedList[pos+1].second;
      else throw -1;
    }

	  RowMerger::RowMerger(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end, bool Sorted) 
				  : Sorted(Sorted), Nonzeros(0), LinkedList(interval_begin,interval_end+1,Row::make_entry(UNDEF,0.0)) 
	  {
		  LinkedList.begin()->first = EOL;
	  }

	  void RowMerger::Resize(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end, bool _Sorted) 
	  {
		  LinkedList.set_interval_beg(interval_begin);
		  LinkedList.set_interval_end(interval_end+1);
		  std::fill(LinkedList.begin(),LinkedList.end(),Row::make_entry(UNDEF,0.0));
		  LinkedList.begin()->first = EOL;
		  Nonzeros = 0;
		  Sorted = _Sorted;
	  }

	  void RowMerger::Resize(Matrix & A, bool _Sorted)
	  {
		  INMOST_DATA_ENUM_TYPE mbeg, mend;
		  A.GetInterval(mbeg,mend);
		  Resize(mbeg,mend,_Sorted);
	  }

	  RowMerger::RowMerger(Matrix & A, bool Sorted) : Sorted(Sorted), Nonzeros(0)
	  {
		  INMOST_DATA_ENUM_TYPE mbeg, mend;
		  A.GetInterval(mbeg,mend);
		  LinkedList.set_interval_beg(mbeg);
		  LinkedList.set_interval_end(mend+1);
		  std::fill(LinkedList.begin(),LinkedList.end(),Row::make_entry(UNDEF,0.0));
		  LinkedList.begin()->first = EOL;
	  }

	  RowMerger::~RowMerger() {}

	  void RowMerger::Clear()
	  {
		  INMOST_DATA_ENUM_TYPE i = LinkedList.begin()->first, j;
		  LinkedList.begin()->first = EOL;
		  while( i != EOL )
		  {
			  j = LinkedList[i].first;
			  LinkedList[i].first = UNDEF;
        LinkedList[i].second = 0.0;
			  i = j;
		  }
		  Nonzeros = 0;
	  }

	  void RowMerger::PushRow(INMOST_DATA_REAL_TYPE coef, Row & r, bool PreSortRow)
	  {
		  if( Sorted && PreSortRow ) std::sort(r.Begin(),r.End());
      //if( Nonzeros != 0 ) printf("nnz %d link %p proc %d\n",Nonzeros,this,omp_get_thread_num());
		  assert(Nonzeros == 0); //Linked list should be empty
		  assert(LinkedList.begin()->first == EOL); //again check that list is empty
		  INMOST_DATA_ENUM_TYPE index = LinkedList.get_interval_beg();
		  Row::iterator it = r.Begin(), jt;
		  while( it != r.End() )
		  {
			  LinkedList[index].first = it->first+1;
			  LinkedList[it->first+1].first = EOL;
			  LinkedList[it->first+1].second = it->second*coef;
			  index = it->first+1;
			  ++Nonzeros;
			  jt = it;
			  ++it;
			  assert(!Sorted || it == r.End() || jt->first < it->first);
		  }
	  }

	  void RowMerger::AddRow(INMOST_DATA_REAL_TYPE coef, Row & r, bool PreSortRow)
	  {
		  if( Sorted && PreSortRow ) std::sort(r.Begin(),r.End());
		  INMOST_DATA_ENUM_TYPE index = LinkedList.get_interval_beg(), next;
		  Row::iterator it = r.Begin(), jt;
		  while( it != r.End() )
		  {
			  if( LinkedList[it->first+1].first != UNDEF )
				  LinkedList[it->first+1].second += coef*it->second;
			  else if( Sorted )
			  {
				  next = index;
				  while(next < it->first+1)
				  {
					  index = next;
					  next = LinkedList[index].first;
				  }
				  assert(index < it->first+1);
				  assert(it->first+1 < next);
				  LinkedList[index].first = it->first+1;
				  LinkedList[it->first+1].first = next;
				  LinkedList[it->first+1].second = coef*it->second;
				  ++Nonzeros;
			  }
			  else
			  {
				  LinkedList[it->first+1].first = LinkedList[index].first;
				  LinkedList[it->first+1].second = coef*it->second;
				  LinkedList[index].first = it->first+1;
				  ++Nonzeros;
			  }
			  jt = it;
			  ++it;
			  assert(!Sorted || it == r.End() || jt->first < it->first);
		  }
	  }

	  void RowMerger::RetriveRow(Row & r)
	  {
		  r.Resize(Nonzeros);
		  INMOST_DATA_ENUM_TYPE i = LinkedList.begin()->first, k = 0;
		  while( i != EOL )
		  {
			  if( LinkedList[i].second )
			  {
				  r.GetIndex(k) = i-1;
				  r.GetValue(k) = LinkedList[i].second;
				  ++k;
			  }
			  i = LinkedList[i].first;
		  }
		  r.Resize(k);
	  }









    Vector::Vector(std::string _name, INMOST_DATA_ENUM_TYPE start, INMOST_DATA_ENUM_TYPE end, INMOST_MPI_Comm _comm) :data(start,end) 
	  {
		  comm = _comm;
		  name = _name;
		  is_parallel = false;
	  }
	
	  Vector::Vector(const Vector & other) : data(other.data) 
	  {
		  comm = other.comm;
		  name = other.name;
		  is_parallel = other.is_parallel;
	  }
	
	  Vector & Vector::operator =(Vector const & other) 
	  {
		  comm = other.comm;
		  data = other.data; 
		  name = other.name; 
		  is_parallel = other.is_parallel;
		  return *this;
	  }
	
	  Vector::~Vector() 
	  {
	  }

    INMOST_DATA_REAL_TYPE   Row::RowVec(Vector & x) const
	  {
		  INMOST_DATA_REAL_TYPE ret = 0;
		  INMOST_DATA_ENUM_TYPE end = Size();
		  for(INMOST_DATA_ENUM_TYPE i = 0; i < end; i++) ret = ret + x[GetIndex(i)]*GetValue(i);
		  return ret;
	  }
	
	  void Matrix::MatVec(INMOST_DATA_REAL_TYPE alpha, Vector & x, INMOST_DATA_REAL_TYPE beta, Vector & out) const //y = alpha*A*x + beta * y
	  {
		  INMOST_DATA_ENUM_TYPE mbeg, mend;
		  INMOST_DATA_INTEGER_TYPE ind, imbeg, imend;
		  if( out.Empty() )
		  {
			  INMOST_DATA_ENUM_TYPE vbeg,vend;
			  GetInterval(vbeg,vend);
			  out.SetInterval(vbeg,vend);
		  }
		  //CHECK SOMEHOW FOR DEBUG THAT PROVIDED VECTORS ARE OK
		  //~ assert(GetFirstIndex() == out.GetFirstIndex());
		  //~ assert(Size() == out.Size());
		  GetInterval(mbeg,mend);
		  imbeg = mbeg;
		  imend = mend;
  #if defined(USE_OMP)
  #pragma omp for private(ind)
  #endif
		  for(ind = imbeg; ind < imend; ++ind) //iterate rows of matrix
			  out[ind] = beta * out[ind] + alpha * (*this)[ind].RowVec(x);
		  // outer procedure should update out vector, if needed
	  }

    
	  void Matrix::MatVecTranspose(INMOST_DATA_REAL_TYPE alpha, Vector & x, INMOST_DATA_REAL_TYPE beta, Vector & out) const //y = alpha*A*x + beta * y
	  {
		  INMOST_DATA_ENUM_TYPE mbeg, mend;
		  INMOST_DATA_INTEGER_TYPE ind, imbeg, imend;
		  if( out.Empty() )
		  {
			  INMOST_DATA_ENUM_TYPE vbeg,vend;
			  GetInterval(vbeg,vend);
			  out.SetInterval(vbeg,vend);
		  }
		  //CHECK SOMEHOW FOR DEBUG THAT PROVIDED VECTORS ARE OK
		  //~ assert(GetFirstIndex() == out.GetFirstIndex());
		  //~ assert(Size() == out.Size());
		  GetInterval(mbeg,mend);
		  imbeg = mbeg;
		  imend = mend;
		  if( beta ) for(Vector::iterator it = out.Begin(); it != out.End(); ++it) (*it) *= beta;
  #if defined(USE_OMP)
  #pragma omp for private(ind)
  #endif
		  for(ind = imbeg; ind < imend; ++ind)
		  {
			  for(Row::const_iterator it = (*this)[ind].Begin(); it != (*this)[ind].End(); ++it)
				  out[it->first] += alpha * x[ind] * it->second;
		  }
		  // outer procedure should update out vector, if needed
	  }
	

	  Matrix::Matrix(std::string _name, INMOST_DATA_ENUM_TYPE start, INMOST_DATA_ENUM_TYPE end, INMOST_MPI_Comm _comm) 
	  :data(start,end) 
	  {
		  is_parallel = false;
		  comm = _comm;
		  SetInterval(start,end);
		  name = _name;
	  }
	
	  Matrix::Matrix(const Matrix & other) :data(other.data) 
	  {
		  comm = other.comm;
		  name = other.name;
	  }
	
	  Matrix & Matrix::operator =(Matrix const & other) 
	  {
		  comm = other.comm;
		  data = other.data; 
		  name = other.name; return *this;
	  }
	
	  Matrix::~Matrix() 
	  {
	  }
	
    void      Matrix::MoveRows(INMOST_DATA_ENUM_TYPE from, INMOST_DATA_ENUM_TYPE to, INMOST_DATA_ENUM_TYPE size)
	  {
		  INMOST_DATA_ENUM_TYPE i = to + size, j = from + size;
		  if( size > 0 && to != from )
			  while( j != from ) data[--j].MoveRow(data[--i]);
      if( !text.empty() )
      {
        i = to + size;
        j = from + size;
        if( size > 0 && to != from )
			    while( j != from ) text[--j] = text[--i];
      }
	  }

    
	  void     Matrix::Load(std::string file, INMOST_DATA_ENUM_TYPE mbeg, INMOST_DATA_ENUM_TYPE mend)
	  {
		  char str[16384];
		  std::ifstream input(file.c_str());
		  if( input.fail() ) throw -1;
		  int state = 0, k;
		  INMOST_DATA_ENUM_TYPE mat_size, max_lines, row, col, mat_block;
		  INMOST_DATA_REAL_TYPE val;
		  int size = 1, rank = 0;
  #if defined(USE_MPI)
      int flag = 0;
      MPI_Initialized(&flag);
		  if( flag && mend == ENUMUNDEF && mbeg == ENUMUNDEF )
		  {
			  MPI_Comm_rank(GetCommunicator(),&rank);
			  MPI_Comm_size(GetCommunicator(),&size);
		  }
  #endif
		  int line = 0;
		  while( !input.getline(str,16384).eof() )
		  {
			  line++;
			  k = 0; while( isspace(str[k]) ) k++;
			  if( str[k] == '%' || str[k] == '\0' ) continue;
			  std::istringstream istr(str+k);
			  switch(state)
			  {
			  case 0: 
				  istr >> mat_size >> mat_size >> max_lines; state = 1; 
				  mat_block = mat_size/size;
				  if( mbeg == ENUMUNDEF ) mbeg = rank*mat_block;
				  if( mend == ENUMUNDEF ) 
				  {
					  if( rank == size-1 ) mend = mat_size;
					  else mend = mbeg+mat_block;
				  }
				  SetInterval(mbeg,mend);
				  //~ std::cout << rank << " my interval " << mbeg << ":" << mend << std::endl;
			  break;
			  case 1: 
				  istr >> row >> col >> val;
				  row--; col--;
				  if( row >= mbeg && row < mend ) data[row][col] = val; 
			  break;
			  }
		  }
		  int nonzero = 0;
		  for(iterator it = Begin(); it != End(); ++it) nonzero += it->Size();
		  //~ std::cout << rank << " total nonzero " << max_lines << " my nonzero " << nonzero << std::endl;
		  input.close();
	  }

    void     Vector::Load(std::string file, INMOST_DATA_ENUM_TYPE mbeg, INMOST_DATA_ENUM_TYPE mend)
	  {
		  char str[16384];
		  std::ifstream input(file.c_str());
		  if( input.fail() ) throw -1;
		  int state = 0, k;
		  INMOST_DATA_ENUM_TYPE vec_size, vec_block, ind = 0;
		  INMOST_DATA_REAL_TYPE val;
		  int size = 1, rank = 0;
  #if defined(USE_MPI)
      int flag = 0;
      MPI_Initialized(&flag);
		  if( flag && mend == ENUMUNDEF && mbeg == ENUMUNDEF )
		  {
			  MPI_Comm_rank(GetCommunicator(),&rank);
			  MPI_Comm_size(GetCommunicator(),&size);
		  }
  #endif
		  while( !input.getline(str,16384).eof() )
		  {
			  k = 0; while( isspace(str[k]) ) k++;
			  if( str[k] == '%' || str[k] == '\0' ) continue;
			  std::istringstream istr(str+k);
			  switch(state)
			  {
			  case 0: 
				  istr >> vec_size; state = 1; 
				  vec_block = vec_size/size;
				  if( mbeg == ENUMUNDEF ) mbeg = rank*vec_block;
				  if( mend == ENUMUNDEF ) 
				  {
					  if( rank == size-1 ) mend = vec_size;
					  else mend = mbeg+vec_block;
				  }
				  SetInterval(mbeg,mend);
			  break;
			  case 1: 
				  istr >> val;
				  if( ind >= mbeg && ind < mend ) data[ind] = val;
				  ind++;
			  break;
			  }
		  }
		  input.close();
	  }

    
	  void     Vector::Save(std::string file)
	  {
		  INMOST_DATA_ENUM_TYPE vecsize = Size();
		
#if defined(USE_MPI)
		  int rank = 0, size = 1;
		  {
			  MPI_Comm_rank(GetCommunicator(),&rank);
			  MPI_Comm_size(GetCommunicator(),&size);
			  INMOST_DATA_ENUM_TYPE temp = vecsize;
			  MPI_Allreduce(&temp,&vecsize,1,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,GetCommunicator());
		  }
#endif
		  std::stringstream rhs(std::ios::in | std::ios::out);
		  rhs << std::scientific;
		  rhs.precision(15);
		  for(iterator it = Begin(); it != End(); ++it) rhs << *it << std::endl;
#if defined(USE_MPI) && defined(USE_MPI_FILE) // Use mpi files
		  { 
			  int ierr;
			  MPI_File fh;
			  MPI_Status stat;
			  ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(file.c_str()), MPI_MODE_CREATE | MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			  ierr = MPI_File_close(&fh);
			  if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			  ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(file.c_str()),MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&fh);
			  if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			  if( rank == 0 )
			  {
				  std::stringstream header;
				  //header << "% vector " << name << std::endl;
				  //header << "% is written by INMOST" << std::endl;
				  //header << "% by MPI_File_* api" << std::endl;
				  header << vecsize << std::endl;
				  ierr = MPI_File_write_shared(fh,const_cast<char *>(header.str().c_str()),static_cast<int>(header.str().size()),MPI_CHAR,&stat);
				  if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			  }
			  ierr = MPI_File_write_ordered(fh,const_cast<char *>(rhs.str().c_str()),static_cast<int>(rhs.str().size()),MPI_CHAR,&stat);
			  if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			  ierr = MPI_File_close(&fh);
			  if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
		  }
#elif defined(USE_MPI) //USE_MPI alternative
		  std::string senddata = rhs.str(), recvdata;
		  int sendsize = static_cast<int>(senddata.size());
		  std::vector<int> recvsize(size), displ(size);
		  MPI_Gather(&sendsize,1,MPI_INT,&recvsize[0],1,MPI_INT,0,GetCommunicator());
		  if( rank == 0 )
		  {
			  int totsize = recvsize[0];
			
			  displ[0] = 0;
			  for(int i = 1; i < size; i++) 
			  {
				  totsize += recvsize[i];
				  displ[i] = displ[i-1]+recvsize[i-1];
			  }
			  recvdata.resize(totsize);
		  }
		  else recvdata.resize(1); //protect from dereferencing null
		  MPI_Gatherv(&senddata[0],sendsize,MPI_CHAR,&recvdata[0],&recvsize[0],&displ[0],MPI_CHAR,0,GetCommunicator());
		  if( rank == 0 )
		  {
			  std::fstream output(file.c_str(),std::ios::out);
			  output << vecsize << std::endl;
			  output << recvdata;
		  }
#else
		  std::fstream output(file.c_str(),std::ios::out);
		  //output << "% vector " << name << std::endl;
		  //output << "% is written by INMOST" << std::endl;
		  //output << "% by sequential write" << std::endl;
		  output << vecsize << std::endl;
		  output << rhs.rdbuf();
#endif
	  }


    void     Matrix::Save(std::string file)
	  {
		  INMOST_DATA_ENUM_TYPE matsize = Size(), nonzero = 0, row = GetFirstIndex()+1;
		
		  for(iterator it = Begin(); it != End(); ++it) nonzero += it->Size();
#if defined(USE_MPI)
		  int rank = 0, size = 1;
		  {
			  MPI_Comm_rank(GetCommunicator(),&rank);
			  MPI_Comm_size(GetCommunicator(),&size);
			  INMOST_DATA_ENUM_TYPE temp_two[2] = {matsize,nonzero}, two[2];
			  MPI_Allreduce(temp_two,two,2,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,GetCommunicator());
			  matsize = two[0];
			  nonzero = two[1];
		  }
#endif
		  std::stringstream mtx(std::ios::in | std::ios::out);
		  mtx << std::scientific;
		  mtx.precision(15);
		  for(iterator it = Begin(); it != End(); ++it)
		  {
        if( !text.empty() ) mtx << "% " << Annotation(it-Begin()).c_str() << "\n";
			  for(Row::iterator jt = it->Begin(); jt != it->End(); ++jt)
        {
				  mtx << row << " " << jt->first+1 << " " << jt->second << "\n";
        }
			  ++row;
		  }
#if defined(USE_MPI) && defined(USE_MPI_FILE) // USE_MPI2?
		  {
			  int ierr;
			  MPI_File fh;
			  MPI_Status stat;
			  ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(file.c_str()), MPI_MODE_CREATE | MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			  ierr = MPI_File_close(&fh);
			  if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			  ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(file.c_str()),MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&fh);
			  if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			  if( rank == 0 )
			  {
				  std::stringstream header;
				  header << "%%MatrixMarket matrix coordinate real general" << std::endl;
				  header << "% matrix " << name << std::endl;
				  header << "% is written by INMOST" << std::endl;
				  header << "% by MPI_File_* api" << std::endl;
				  header << matsize << " " << matsize << " " << nonzero << std::endl;
				  //std::string header_data(header.str());
				  ierr = MPI_File_write_shared(fh,const_cast<char *>(header.str().c_str()),static_cast<int>(header.str().size()),MPI_CHAR,&stat);
				  if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			  }
			  ierr = MPI_File_write_ordered(fh,const_cast<char *>(mtx.str().c_str()),static_cast<int>(mtx.str().size()),MPI_CHAR,&stat);
			  if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
			  ierr = MPI_File_close(&fh);
			  if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),-1);
		  }
#elif defined(USE_MPI)//USE_MPI alternative
		  std::string senddata = mtx.str(), recvdata;
		  int sendsize = static_cast<int>(senddata.size());
		  std::vector<int> recvsize(size), displ(size);
		  MPI_Gather(&sendsize,1,MPI_INT,&recvsize[0],1,MPI_INT,0,GetCommunicator());
		  if( rank == 0 )
		  {
			  int totsize = recvsize[0];
			
			  displ[0] = 0;
			  for(int i = 1; i < size; i++) 
			  {
				  totsize += recvsize[i];
				  displ[i] = displ[i-1]+recvsize[i-1];
			  }
			  recvdata.resize(totsize);
		  }
		  else recvdata.resize(1); //protect from dereferencing null
		  MPI_Gatherv(&senddata[0],sendsize,MPI_CHAR,&recvdata[0],&recvsize[0],&displ[0],MPI_CHAR,0,GetCommunicator());
		  if( rank == 0 )
		  {
			  std::fstream output(file.c_str(),std::ios::out);
			  output << "%%MatrixMarket matrix coordinate real general" << std::endl;
			  output << "% matrix " << name << std::endl;
			  output << "% is written by INMOST" << std::endl;
			  output << "% by MPI_Gather* api and sequential write" << std::endl;
			  output << matsize << " " << matsize << " " << nonzero << std::endl;
			  output << recvdata;
		  }
#else
		  std::fstream output(file.c_str(),std::ios::out);
		  output << "%%MatrixMarket matrix coordinate real general" << std::endl;
		  output << "% matrix " << name << std::endl;
		  output << "% is written by INMOST" << std::endl;
		  output << "% by sequential write " << std::endl;
		  output << matsize << " " << matsize << " " << nonzero << std::endl;
		  output << mtx.rdbuf();
#endif
	  }
    std::string & Matrix::Annotation(INMOST_DATA_ENUM_TYPE row) 
    {
      if( text.empty() ) 
      {
        text.set_interval_beg(GetFirstIndex());
        text.set_interval_end(GetLastIndex());
      }
      return text[row];
    }
    const std::string & Matrix::Annotation(INMOST_DATA_ENUM_TYPE row) const
    {
      if( text.empty() ) return stubstring;
      return text[row];
    }
    Row::Row() :data() 
    {
#if defined(USE_OMP)
      omp_init_lock(&lock);
#endif
      modified_pattern = marker = false;
    }
    Row::Row(const Row & other) :marker(other.marker),data(other.data) 
    { 
#if defined(USE_OMP)
      omp_init_lock(&lock);
#endif
      modified_pattern = other.modified_pattern; 
    }
    Row::Row(entry * pbegin, entry * pend) :data(pbegin, pend) 
    { 
#if defined(USE_OMP)
      omp_init_lock(&lock);
#endif
      modified_pattern = true; marker = false; 
    }
    Row & Row::operator = (Row const & other) 
    { 
      data = other.data; 
      marker = other.marker; 
      return *this; 
    }
    Row::~Row()
    {
#if defined(USE_OMP)
      omp_destroy_lock(&lock);
#endif
    }
    void  Row::Swap(Row & other) 
    { 
      data.swap(other.data); 
      bool tmp = marker; 
      marker = other.marker; 
      other.marker = tmp; 
#if defined(USE_OMP)
      //swap locks?
#endif
    }
    void Matrix::SetInterval(INMOST_DATA_ENUM_TYPE   start, INMOST_DATA_ENUM_TYPE   end)
    {
      data.set_interval_beg(start); 
      data.set_interval_end(end);
      if( !text.empty() )
      {
        text.set_interval_beg(start);
        text.set_interval_end(end);
      }
    }
    void Matrix::ShiftInterval(INMOST_DATA_ENUM_TYPE shift) 
    {
      data.shift_interval(shift);
      if( !text.empty() ) text.shift_interval(shift);
    }


  }
}