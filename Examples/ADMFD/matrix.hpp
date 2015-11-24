#ifndef __MATRIX_H
#define __MATRIX_H
#include "inmost.h"

using namespace INMOST;
// Matrix with n columns and m rows
//   __m__
//  |     |
// n|     |
//  |_____|
//
// todo: 
// 1. expression templates for operations
// 2. (ok) template matrix type for AD variables
template<typename Var>
class Matrix
{
public:
  typedef unsigned enumerator;
protected:
  array<Var> space;
  enumerator n, m;
public:
  static Matrix Unit(enumerator pn)
  {
    Matrix ret(pn,pn);
    for(enumerator i = 0; i < pn; ++i)
    {
      for(enumerator j = 0; j < pn; ++j)
      {
        if( i == j ) 
          ret(i,j) = 1;
        else
          ret(i,j) = 0;
      }
    }
    return ret;
  }
  Matrix(Var * pspace, enumerator pn, enumerator pm) : space(pspace,pspace+pn*pm), n(pn), m(pm) {}
  Matrix(enumerator pn, enumerator pm) : space(pn*pm), n(pn), m(pm) {}
  Matrix(const Matrix & other) : space(other.n*other.m), n(other.n), m(other.m) 
  {
    for(enumerator i = 0; i < n*m; ++i)
      space[i] = other.space[i];
  }
  ~Matrix() {remove();}
  Matrix & operator =(Matrix const & other)
  {
    remove();
    space.resize(other.n*other.m);
    for(enumerator i = 0; i < other.n*other.m; ++i)
      space[i] = other.space[i];
    n = other.n;
    m = other.m;
    return *this;
  }
  void remove() {}
  // i is in [0,n] - row index
  // j is in [0,m] - column index
  Var & operator()(enumerator i, enumerator j)
  {
    assert(i >= 0 && i < n);
    assert(j >= 0 && j < m);
    assert(i*m+j < n*m); //overflow check?
    return space[i*m+j];
  }
  Var operator()(enumerator i, enumerator j) const
  {
    assert(i >= 0 && i < n);
    assert(j >= 0 && j < m);
    assert(i*m+j < n*m); //overflow check?
    return space[i*m+j];
  }
  Matrix operator-(const Matrix & other) const
  {
    assert(n == other.n);
    assert(m == other.m);
    Matrix ret(n,m); //check RVO
    for(enumerator k = 0; k < n*m; ++k) ret.space[k] = space[k]-other.space[k];
    return ret;
  }
  Matrix & operator-=(const Matrix & other)
  {
    assert(n == other.n);
    assert(m == other.m);
    for(enumerator k = 0; k < n*m; ++k) space[k] -= other.space[k];
    return *this;
  }
  Matrix operator+(const Matrix & other) const
  {
    assert(n == other.n);
    assert(m == other.m);
    Matrix ret(n,m); //check RVO
    for(enumerator k = 0; k < n*m; ++k) ret.space[k] = space[k]+other.space[k];
    return ret;
  }
  Matrix & operator+=(const Matrix & other)
  {
    assert(n == other.n);
    assert(m == other.m);
    for(enumerator k = 0; k < n*m; ++k) space[k] += other.space[k];
    return *this;
  }
  Matrix operator*(Var coef) const
  {
    Matrix ret(n,m); //check RVO
    for(enumerator k = 0; k < n*m; ++k) ret.space[k] = space[k]*coef;
    return ret;
  }
  Matrix & operator*=(Var coef)
  {
    for(enumerator k = 0; k < n*m; ++k) space[k] *= coef;
    return *this;
  }
  Matrix operator/(Var coef) const
  {
    Matrix ret(n,m); //check RVO
    for(enumerator k = 0; k < n*m; ++k) ret.space[k] = space[k]/coef;
    return ret;
  }
  Matrix & operator/=(Var coef)
  {
    for(enumerator k = 0; k < n*m; ++k) space[k] /= coef;
    return *this;
  }
  Matrix operator*(const Matrix & other) const
  {
    assert(m == other.n);
    Matrix ret(n,other.m); //check RVO
    for(enumerator i = 0; i < n; ++i) //loop rows
    {
      for(enumerator j = 0; j < other.m; ++j) //loop columns
      {
        Var tmp = 0.0;
        for(enumerator k = 0; k < m; ++k)
          tmp += (*this)(i,k)*other(k,j);
        ret(i,j) = tmp;
      }
    }
    return ret;
  }
  /*
  template<typename InType>
  Matrix<InType> operator*(const Matrix<InType> & other)
  {
    assert(m == other.n);
    Matrix<InType> ret(n,other.m); //check RVO
    for(enumerator i = 0; i < n; ++i) //loop rows
    {
      for(enumerator j = 0; j < other.m; ++j) //loop columns
      {
        InType tmp = 0.0;
        for(enumerator k = 0; k < m; ++k)
          tmp += (*this)(i,k)*other(k,j);
        ret(i,j) = tmp;
      }
    }
    return ret;
  }
  template<typename InType>
  Matrix operator*(const Matrix<InType> & other)
  {
    assert(m == other.n);
    Matrix ret(n,other.m); //check RVO
    for(enumerator i = 0; i < n; ++i) //loop rows
    {
      for(enumerator j = 0; j < other.m; ++j) //loop columns
      {
        Var tmp = 0.0;
        for(enumerator k = 0; k < m; ++k)
          tmp += (*this)(i,k)*other(k,j);
        ret(i,j) = tmp;
      }
    }
    return ret;
  }
  */
  Matrix Transpose() const
  {
    Matrix ret(m,n);
    for(enumerator i = 0; i < n; ++i)
    {
      for(enumerator j = 0; j < m; ++j)
      {
        ret(j,i) = (*this)(i,j);
      }
    }
    return ret;
  }
  Matrix Invert() const
  {
    Matrix ret(m,n);
    Matrix At = Transpose(); //m by n matrix
    Matrix AtB = At*Unit(n); //m by n matrix
    Matrix AtA = At*(*this); //m by m matrix
    enumerator * order = new enumerator [m];
    for(enumerator i = 0; i < m; ++i) order[i] = i;
    for(enumerator i = 0; i < m; i++)
	  {
		  enumerator maxk = i, maxq = i, temp2;
      Var max, temp;
		  max = fabs(AtA(maxk,maxq));
		  //Find best pivot
      for(enumerator k = i; k < m; k++) // over rows
		  {
			  for(enumerator q = i; q < m; q++) // over columns
			  {
				  if( fabs(AtA(k,q)) > max )
				  {
					  max = fabs(AtA(k,q));
					  maxk = k;
					  maxq = q;
				  }
			  }
		  }
		  //Exchange rows
		  if( maxk != i ) 
		  {
			  for(enumerator q = 0; q < m; q++) // over columns of A
			  {
				  temp = AtA(maxk,q);
				  AtA(maxk,q) = AtA(i,q);
				  AtA(i,q) = temp;
			  }
			  //exchange rhs
			  for(enumerator q = 0; q < n; q++) // over columns of B
			  {
				  temp = AtB(maxk,q);
				  AtB(maxk,q) = AtB(i,q);
				  AtB(i,q) = temp;
			  }
		  }
		  //Exchange columns
		  if( maxq != i ) 
		  {
			  for(enumerator k = 0; k < m; k++) //over rows
			  {
				  temp = AtA(k,maxq);
				  AtA(k,maxq) = AtA(k,i);
				  AtA(k,i) = temp;
			  }
			  //remember order in sol
			  {
				  temp2 = order[maxq];
				  order[maxq] = order[i];
				  order[i] = temp2;
			  }
		  }

		  if( fabs(AtA(i,i)) < 1.0e-54 )
		  {
			  bool ok = true;
			  for(enumerator k = 0; k < n; k++) // over columns of B
			  {
				  if( fabs(AtB(i,k)/1.0e-54) > 1 )
				  {
					  ok = false;
					  break;
				  }
			  }
			  if( ok ) AtA(i,i) = AtA(i,i) < 0.0 ? - 1.0e-12 : 1.0e-12;
			  else assert(false && "Matrix is singular");
		  }
		  for(enumerator k = i+1; k < m; k++)
		  {
			  AtA(i,k) /= AtA(i,i);
			  AtA(k,i) /= AtA(i,i);
		  }
		  for(enumerator k = i+1; k < m; k++)
		  for(enumerator q = i+1; q < m; q++)
		  {
			  AtA(k,q) -= AtA(k,i) * AtA(i,i) * AtA(i,q);
		  }
		  for(enumerator k = 0; k < n; k++)
		  {
			  for(enumerator j = i+1; j < m; j++) //iterate over columns of L
			  {
				  AtB(j,k) -= AtB(i,k) * AtA(j,i);
			  }
			  AtB(i,k) /= AtA(i,i);
		  }
	  }

	  for(enumerator k = 0; k < n; k++)
	  {
		  for(enumerator i = m; i-- > 0; ) //iterate over rows of U
			  for(enumerator j = i+1; j < m; j++) 
			  {
				  AtB(i,k) -= AtB(j,k) * AtA(i,j);
			  }
		  for(enumerator i = 0; i < m; i++)
			  ret(order[i],k) = AtB(i,k);
	  }
    delete [] order;
    return ret;
  }
  void Zero()
  {
    for(enumerator i = 0; i < n*m; ++i) space[i] = 0.0;
  }
  Var Trace() const
  {
    assert(n == m);
    Var ret = 0.0;
    for(enumerator i = 0; i < n; ++i) ret += (*this)(i,i);
    return ret;
  }
  Var * data() {return space.data();}
  const Var * data() const {return space.data();}
  enumerator Rows() const {return n;}
  enumerator Cols() const {return m;}
  void Print(Storage::real threshold = 1.0e-10) const
  {
    for(enumerator k = 0; k < n; ++k)
    {
      for(enumerator l = 0; l < m; ++l) 
      {
        if( fabs((*this)(k,l)) > threshold )
          std::cout << (*this)(k,l);
        else 
          std::cout << 0;
        std::cout << " ";
      }
      std::cout << std::endl;
    }
  }
  bool isSymmetric() const
  {
    if( n != m ) return false;
    for(enumerator k = 0; k < n; ++k)
    {
      for(enumerator l = k+1; l < n; ++l)
        if( fabs((*this)(k,l)-(*this)(l,k)) > 1.0e-7 )
          return false;
    }
    return true;
  }
  Var DotProduct(const Matrix & other) const
  {
    assert(n == other.n);
    assert(m == other.m);
    Var ret = 0.0;
    for(enumerator i = 0; i < n*m; ++i) ret += space[i]*other.space[i];
    return ret;
  }
};


typedef Matrix<Storage::real> rMatrix; //shortcut for real matrix
typedef Matrix<Storage::var> vMatrix; //shortcut for matrix with variations




// Multiplication of matrices with promotion
__INLINE vMatrix operator *(const rMatrix & A, const vMatrix & B)
{
  assert(A.Cols() == B.Rows());
  vMatrix ret(A.Rows(),B.Cols()); //check RVO
  for(rMatrix::enumerator i = 0; i < A.Rows(); ++i) //loop rows
  {
    for(rMatrix::enumerator j = 0; j < B.Cols(); ++j) //loop columns
    {
      variable tmp = 0.0;
      for(rMatrix::enumerator k = 0; k < A.Cols(); ++k)
        tmp += A(i,k)*B(k,j);
      ret(i,j) = tmp;
    }
  }
  return ret;
}

// Multiplication of matrices with promotion
__INLINE vMatrix operator *(const vMatrix & A, const rMatrix & B)
{
  assert(A.Cols() == B.Rows());
  vMatrix ret(A.Rows(),B.Cols()); //check RVO
  for(rMatrix::enumerator i = 0; i < A.Rows(); ++i) //loop rows
  {
    for(rMatrix::enumerator j = 0; j < B.Cols(); ++j) //loop columns
    {
      variable tmp = 0.0;
      for(rMatrix::enumerator k = 0; k < A.Cols(); ++k)
        tmp += A(i,k)*B(k,j);
      ret(i,j) = tmp;
    }
  }
  return ret;
}

//retrive permeability matrix from stored data
__INLINE rMatrix FromTensor(Storage::real_array & K)
{
  rMatrix Kc(3,3);
  switch(K.size())
  {
  case 1: //scalar permeability tensor
    Kc.Zero();
    Kc(0,0) = Kc(1,1) = Kc(2,2) = K[0];
  case 3:
    Kc.Zero(); //diagonal permeability tensor 
    Kc(0,0) = K[0]; //KXX
    Kc(1,1) = K[1]; //KYY
    Kc(2,2) = K[2]; //KZZ
    break;
  case 6: //symmetric permeability tensor
    Kc(0,0) = K[0]; //KXX
    Kc(0,1) = Kc(1,0) = K[1]; //KXY
    Kc(0,2) = Kc(2,0) = K[2]; //KXZ
    Kc(1,1) = K[3]; //KYY
    Kc(1,2) = Kc(2,1) = K[4]; //KYZ
    Kc(2,2) = K[5]; //KZZ
    break;
  case 9: //full permeability tensor
    Kc(0,0) = K[0]; //KXX
    Kc(0,1) = K[1]; //KXY
    Kc(0,2) = K[2]; //KXZ
    Kc(1,0) = K[3]; //KYX
    Kc(1,1) = K[4]; //KYY
    Kc(1,2) = K[5]; //KYZ
    Kc(2,0) = K[6]; //KZX
    Kc(2,1) = K[7]; //KZY
    Kc(2,2) = K[8]; //KZZ
    break;
  }
  return Kc;
}

//retrive vector in matrix form from array
__INLINE rMatrix FromVector(Storage::real n[3])
{
  return rMatrix(n,3,1);
}

//create diagonal matrix from array
__INLINE rMatrix FromDiagonal(Storage::real * r, rMatrix::enumerator size)
{
  rMatrix ret(size,size);
  ret.Zero();
  for(rMatrix::enumerator k = 0; k < size; ++k) ret(k,k) = r[k];
  return ret;
}

//create diagonal matrix from array of inversed values
__INLINE rMatrix FromDiagonalInverse(Storage::real * r, rMatrix::enumerator size)
{
  rMatrix ret(size,size);
  ret.Zero();
  for(rMatrix::enumerator k = 0; k < size; ++k) ret(k,k) = 1.0/r[k];
  return ret;
}

#endif