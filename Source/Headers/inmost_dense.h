#ifndef INMOST_DENSE_INCLUDED
#define INMOST_DENSE_INCLUDED
#include "inmost_common.h"
#if defined(USE_AUTODIFF)
#include "inmost_expression.h"
#endif

// Matrix with n columns and m rows
//   __m__
//  |     |
// n|     |
//  |_____|
//
// todo: 
// 1. expression templates for operations
//    (???) how to for multiplication?
// 2. (ok) template matrix type for AD variables


namespace INMOST
{
  

  template<class A, class B> struct Promote;
  template<> struct Promote<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> {typedef INMOST_DATA_REAL_TYPE type;};
#if defined(USE_AUTODIFF)
  template<> struct Promote<INMOST_DATA_REAL_TYPE, variable>  {typedef variable type;};
  template<> struct Promote<variable, INMOST_DATA_REAL_TYPE>  {typedef variable type;};
  template<> struct Promote<variable, variable> {typedef variable type;};
#else
  INMOST_DATA_REAL_TYPE get_value(INMOST_DATA_REAL_TYPE x) {return x;}
#endif

  template<typename Var>
  class Matrix
  {
  public:
    typedef unsigned enumerator;
  protected:
    array<Var> space;
    enumerator n, m;
  

    static Var sign_func(const Var & a, const Var & b) {return (b >= 0.0 ? 1.0 : -1.0)*fabs(a);}
	  static Var max_func(const Var & x, const Var & y) { return x > y ? x : y; }
	  static Var pythag(const Var & a, const Var & b)
	  {
		  Var at = fabs(a), bt = fabs(b), ct, result;
		  if (at > bt)       { ct = bt / at; result = sqrt(at) * sqrt(at + ct * bt); }
		  else if (bt > 0.0) { ct = at / bt; result = sqrt(bt) * sqrt(bt + ct * at); }
		  else result = 0.0;
		  return result;
	  }
  public:
    void Swap(Matrix & b)
    {
      space.swap(b.space);
      std::swap(n,b.n);
      std::swap(m,b.m);
    }
    bool SVD(Matrix & U, Matrix & Sigma, Matrix & VT)
    {
      int flag, i, its, j, jj, k, l, nm;
      Var c, f, h, s, x, y, z;
      Var anorm = 0.0, g = 0.0, scale = 0.0;
      if (n < m) 
      {
        bool success = Transpose().SVD(U,Sigma,VT);
        if( success )
        {
          U.Swap(VT);
          U = U.Transpose();
          VT = VT.Transpose();
          return true;
        }
        else return false;
      } //m <= n
      array<Var> _rv1(m);
      shell<Var> rv1(_rv1);
      U = (*this);
      Sigma.Resize(m,m);
      Sigma.Zero();
      VT.Resize(m,m);
  
      std::swap(n,m); //this how original algorithm takes it
      // Householder reduction to bidiagonal form
      for (i = 0; i < n; i++) 
      {
        // left-hand reduction
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
          for (k = i; k < m; k++) scale += fabs(U[k][i]);
          if (get_value(scale)) 
          {
            for (k = i; k < m; k++) 
            {
              U[k][i] /= scale;
              s += U[k][i] * U[k][i];
            }
            f = U[i][i];
            g = -sign_func(sqrt(s), f);
            h = f * g - s;
            U[i][i] = f - g;
            if (i != n - 1) 
            {
              for (j = l; j < n; j++) 
              {
                for (s = 0.0, k = i; k < m; k++) s += U[k][i] * U[k][j];
                f = s / h;
                for (k = i; k < m; k++) U[k][j] += f * U[k][i];
              }
            }
            for (k = i; k < m; k++) U[k][i] *= scale;
          }
        }
        Sigma[i][i] = scale * g;
        // right-hand reduction
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
          for (k = l; k < n; k++) scale += fabs(U[i][k]);
          if (get_value(scale)) 
          {
            for (k = l; k < n; k++) 
            {
              U[i][k] = U[i][k]/scale;
              s += U[i][k] * U[i][k];
            }
            f = U[i][l];
            g = -sign_func(sqrt(s), f);
            h = f * g - s;
            U[i][l] = f - g;
            for (k = l; k < n; k++) rv1[k] = U[i][k] / h;
            if (i != m - 1) 
            {
              for (j = l; j < m; j++) 
              {
                for (s = 0.0, k = l; k < n; k++) s += U[j][k] * U[i][k];
                for (k = l; k < n; k++) U[j][k] += s * rv1[k];
              }
            }
            for (k = l; k < n; k++) U[i][k] *= scale;
          }
        }
        anorm = max_func(anorm,fabs(Sigma[i][i]) + fabs(rv1[i]));
      }

      // accumulate the right-hand transformation
      for (i = n - 1; i >= 0; i--) 
      {
        if (i < n - 1) 
        {
          if (get_value(g)) 
          {
            for (j = l; j < n; j++) VT[j][i] = ((U[i][j] / U[i][l]) / g);
            // double division to avoid underflow
            for (j = l; j < n; j++) 
            {
              for (s = 0.0, k = l; k < n; k++) s += U[i][k] * VT[k][j];
              for (k = l; k < n; k++) VT[k][j] += s * VT[k][i];
            }
          }
          for (j = l; j < n; j++) VT[i][j] = VT[j][i] = 0.0;
        }
        VT[i][i] = 1.0;
        g = rv1[i];
        l = i;
      }

      // accumulate the left-hand transformation
      for (i = n - 1; i >= 0; i--) 
      {
        l = i + 1;
        g = Sigma[i][i];
        if (i < n - 1) 
          for (j = l; j < n; j++) 
            U[i][j] = 0.0;
        if (get_value(g)) 
        {
          g = 1.0 / g;
          if (i != n - 1) 
          {
            for (j = l; j < n; j++) 
            {
              for (s = 0.0, k = l; k < m; k++) s += (U[k][i] * U[k][j]);
              f = (s / U[i][i]) * g;
              for (k = i; k < m; k++) U[k][j] += f * U[k][i];
            }
          }
          for (j = i; j < m; j++) U[j][i] = U[j][i]*g;
        }
        else for (j = i; j < m; j++) U[j][i] = 0.0;
        U[i][i] += 1;
      }

      // diagonalize the bidiagonal form
      for (k = n - 1; k >= 0; k--) 
      {// loop over singular values
        for (its = 0; its < 30; its++) 
        {// loop over allowed iterations
          flag = 1;
          for (l = k; l >= 0; l--) 
          {// test for splitting
            nm = l - 1;
            if (fabs(rv1[l]) + anorm == anorm) 
            {
              flag = 0;
              break;
            }
            if (fabs(Sigma[nm][nm]) + anorm == anorm) 
              break;
          }
          if (flag) 
          {
            c = 0.0;
            s = 1.0;
            for (i = l; i <= k; i++) 
            {
              f = s * rv1[i];
              if (fabs(f) + anorm != anorm) 
              {
                g = Sigma[i][i];
                h = pythag(f, g);
                Sigma[i][i] = h; 
                h = 1.0 / h;
                c = g * h;
                s = (- f * h);
                for (j = 0; j < m; j++) 
                {
                  y = U[j][nm];
                  z = U[j][i];
                  U[j][nm] = (y * c + z * s);
                  U[j][i] = (z * c - y * s);
                }
              }
            }
          }
          z = Sigma[k][k];
          if (l == k) 
          {// convergence
            if (z < 0.0) 
            {// make singular value nonnegative
              Sigma[k][k] = -z;
              for (j = 0; j < n; j++) VT[j][k] = -VT[j][k];
            }
            break;
          }
          if (its >= 30) 
          {
            std::cout << "No convergence after " << its << " iterations" << std::endl;
            std::swap(n,m);
            return false;
          }
          // shift from bottom 2 x 2 minor
          x = Sigma[l][l];
          nm = k - 1;
          y = Sigma[nm][nm];
          g = rv1[nm];
          h = rv1[k];
          f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
          g = pythag(f, 1.0);
          f = ((x - z) * (x + z) + h * ((y / (f + sign_func(g, f))) - h)) / x;
          // next QR transformation
          c = s = 1.0;
          for (j = l; j <= nm; j++) 
          {
            i = j + 1;
            g = rv1[i];
            y = Sigma[i][i];
            h = s * g;
            g = c * g;
            z = pythag(f, h);
            rv1[j] = z;
            c = f / z;
            s = h / z;
            f = x * c + g * s;
            g = g * c - x * s;
            h = y * s;
            y = y * c;
            for (jj = 0; jj < n; jj++) 
            {
              x = VT[jj][j];
              z = VT[jj][i];
              VT[jj][j] = (x * c + z * s);
              VT[jj][i] = (z * c - x * s);
            }
            z = pythag(f, h);
            Sigma[j][j] = z;
            if (z) 
            {
              z = 1.0 / z;
              c = f * z;
              s = h * z;
            }
            f = (c * g) + (s * y);
            x = (c * y) - (s * g);
            for (jj = 0; jj < m; jj++) 
            {
              y = U[jj][j];
              z = U[jj][i];
              U[jj][j] = (y * c + z * s);
              U[jj][i] = (z * c - y * s);
            }
          }
          rv1[l] = 0.0;
          rv1[k] = f;
          Sigma[k][k] = x;
        }
      }
      std::swap(n,m);
      return true;
    }
  
    Matrix(Var * pspace, enumerator pn, enumerator pm) : space(pspace,pspace+pn*pm), n(pn), m(pm) {}
    Matrix(enumerator pn, enumerator pm) : space(pn*pm), n(pn), m(pm) {}
    Matrix(const Matrix & other) : space(other.n*other.m), n(other.n), m(other.m) 
    {
      for(enumerator i = 0; i < n*m; ++i)
        space[i] = other.space[i];
    }
    ~Matrix() {}
    void Resize(enumerator nrows, enumerator mcols)
    {
      if( space.size() != mcols*nrows )
        space.resize(mcols*nrows);
      n = nrows;
      m = mcols;
    }
    Matrix & operator =(Matrix const & other)
    {
      space.resize(other.n*other.m);
      for(enumerator i = 0; i < other.n*other.m; ++i)
        space[i] = other.space[i];
      n = other.n;
      m = other.m;
      return *this;
    }
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
    Matrix operator-() const
    {
      Matrix ret(n,m);
      for(enumerator k = 0; k < n*m; ++k) ret.space[k] = -space[k];
      return ret;
    }
    template<typename typeB>
    Matrix<typename Promote<Var,typeB>::type> operator-(const Matrix<typeB> & other) const
    {
      assert(Rows() == other.Rows());
      assert(Cols() == other.Cols());
      Matrix<typename Promote<Var,typeB>::type> ret(n,m); //check RVO
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
    template<typename typeB>
    Matrix<typename Promote<Var,typeB>::type> operator+(const Matrix<typeB> & other) const
    {
      assert(Rows() == other.Rows());
      assert(Cols() == other.Cols());
      Matrix<typename Promote<Var,typeB>::type> ret(n,m); //check RVO
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
    template<typename typeB>
    Matrix<typename Promote<Var,typeB>::type> operator*(typeB coef) const
    {
      Matrix<typename Promote<Var,typeB>::type> ret(n,m); //check RVO
      for(enumerator k = 0; k < n*m; ++k) ret.space[k] = space[k]*coef;
      return ret;
    }
    Matrix & operator*=(Var coef)
    {
      for(enumerator k = 0; k < n*m; ++k) space[k] *= coef;
      return *this;
    }
    template<typename typeB>
    Matrix<typename Promote<Var,typeB>::type> operator/(typeB coef) const
    {
      Matrix<typename Promote<Var,typeB>::type> ret(n,m); //check RVO
      for(enumerator k = 0; k < n*m; ++k) ret.space[k] = space[k]/coef;
      return ret;
    }
    Matrix & operator/=(Var coef)
    {
      for(enumerator k = 0; k < n*m; ++k) space[k] /= coef;
      return *this;
    }
    template<typename typeB>
    Matrix<typename Promote<Var,typeB>::type> operator*(const Matrix<typeB> & other) const
    {
      assert(Cols() == other.Rows());
      Matrix<typename Promote<Var,typeB>::type> ret(Rows(),other.Cols()); //check RVO
      for(enumerator i = 0; i < Rows(); ++i) //loop rows
      {
        for(enumerator j = 0; j < other.Cols(); ++j) //loop columns
        {
          typename Promote<Var,typeB>::type tmp = 0.0;
          for(enumerator k = 0; k < Cols(); ++k)
            tmp += (*this)(i,k)*other(k,j);
          ret(i,j) = tmp;
        }
      }
      return ret;
    }
    /// performs A*B^{-1}
    /// checks existence of B^{-1} only in debug mode.
    template<typename typeB>
    Matrix<typename Promote<Var,typeB>::type> operator/(const Matrix<typeB> & other) const
    {
      std::pair<Matrix<typeB>,bool> other_inv = other.Invert();
      assert(other_inv.second);
      assert(Cols() == other_inv.Rows());
      Matrix<typename Promote<Var,typeB>::type> ret(n,other.m); //check RVO
      for(enumerator i = 0; i < Rows(); ++i) //loop rows
      {
        for(enumerator j = 0; j < other_inv.Cols(); ++j) //loop columns
        {
          typename Promote<Var,typeB>::type tmp = 0.0;
          for(enumerator k = 0; k < Cols(); ++k)
            tmp += (*this)(i,k)*other_inv.first(k,j);
          ret(i,j) = tmp;
        }
      }
      return ret;
    }
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
    std::pair<Matrix,bool> Invert() const
    {
      std::pair<Matrix,bool> ret = std::make_pair(Matrix(m,n),true);
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
			    else
          {
            ret.second = false;
            return ret;
          }
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
			    ret.first(order[i],k) = AtB(i,k);
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
    void Print(INMOST_DATA_REAL_TYPE threshold = 1.0e-10) const
    {
      for(enumerator k = 0; k < n; ++k)
      {
        for(enumerator l = 0; l < m; ++l) 
        {
          if( fabs(get_value((*this)(k,l))) > threshold )
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
    Var FrobeniusNorm()
    {
      Var ret = 0;
      for(enumerator i = 0; i < n*m; ++i) ret += space[i]*space[i];
      return sqrt(ret);
    }
    static Matrix<Var> FromTensor(Var * K, enumerator size)
    {
      Matrix<Var> Kc(3,3);
      switch(size)
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
    static Matrix<Var> FromVector(Var * n, enumerator size)
    {
      return Matrix(n,size,1);
    }
    //create diagonal matrix from array
    static Matrix<Var> FromDiagonal(Var * r, enumerator size)
    {
      Matrix ret(size,size);
      ret.Zero();
      for(enumerator k = 0; k < size; ++k) ret(k,k) = r[k];
      return ret;
    }
    //create diagonal matrix from array of inversed values
    static Matrix<Var> FromDiagonalInverse(Var * r, enumerator size)
    {
      Matrix ret(size,size);
      ret.Zero();
      for(enumerator k = 0; k < size; ++k) ret(k,k) = 1.0/r[k];
      return ret;
    }
    //Unit matrix
    static Matrix Unit(enumerator pn)
    {
      Matrix ret(pn,pn);
      ret.Zero();
      for(enumerator i = 0; i < pn; ++i) ret(i,i) = 1.0;
      return ret;
    }
  };


  typedef Matrix<INMOST_DATA_REAL_TYPE> rMatrix; //shortcut for real matrix
#if defined(USE_AUTODIFF)
  typedef Matrix<variable> vMatrix; //shortcut for matrix with variations
#endif
}


#endif //INMOST_DENSE_INCLUDED