#ifndef INMOST_DENSE_INCLUDED
#define INMOST_DENSE_INCLUDED
#include "inmost_common.h"
#if defined(USE_AUTODIFF)
#include "inmost_expression.h"
#endif
#include <iomanip>

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
// 3. template container type for data storage.
// 4. option for wrapper container around provided data storage. (to perform matrix operations with existing data)
// 5. class Subset for fortran-like access to matrix.

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
  

    static Var sign_func(const Var & a, const Var & b) {return (b >= 0.0 ? fabs(a) : -fabs(a));}
	  static INMOST_DATA_REAL_TYPE max_func(INMOST_DATA_REAL_TYPE x, INMOST_DATA_REAL_TYPE y) { return x > y ? x : y; }
	  static Var pythag(const Var & a, const Var & b)
	  {
		  Var at = fabs(a), bt = fabs(b), ct, result;
		  if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
		  else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
		  else result = 0.0;
		  return result;
	  }
  public:
    void RemoveRow(enumerator row)
    {
      for(enumerator k = row+1; k < n; ++k)
      {
        for(enumerator l = 0; l < m; ++l)
          (*this)(k-1,l) = (*this)(k,l);
      }
      space.resize((n-1)*m);
      --n;
    }
    void RemoveRows(enumerator first, enumerator last)
    {
      enumerator shift = last-first;
      for(enumerator k = last+1; k < n; ++k)
      {
        for(enumerator l = 0; l < m; ++l)
        (*this)(k-shift-1,l) = (*this)(k,l);
      }
      space.resize((n-shift)*m);
      n-=shift;
    }
    void RemoveColumn(enumerator col)
    {
      Matrix<Var> tmp(n,m-1);
      for(enumerator k = 0; k < n; ++k)
      {
        for(enumerator l = 0; l < col; ++l)
            tmp(k,l) = (*this)(k,l);
        for(enumerator l = col+1; l < m; ++l)
          tmp(k,l-1) = (*this)(k,l);
      }
      this->Swap(tmp);
    }
      void RemoveColumns(enumerator first, enumerator last)
      {
          enumerator shift = last-first;
          Matrix<Var> tmp(n,m-shift);
          for(enumerator k = 0; k < n; ++k)
          {
              for(enumerator l = 0; l < first; ++l)
                  tmp(k,l) = (*this)(k,l);
              for(enumerator l = last+1; l < m; ++l)
                  tmp(k,l-shift-1) = (*this)(k,l);
          }
          this->Swap(tmp);
      }
      void RemoveSubset(enumerator firstrow, enumerator lastrow, enumerator firstcol, enumerator lastcol)
      {
          enumerator shiftrow = lastrow-firstrow;
          enumerator shiftcol = lastcol-firstcol;
          Matrix<Var> tmp(n-shiftrow, m-shiftcol);
          for(enumerator k = 0; k < firstrow; ++k)
          {
              for(enumerator l = 0; l < firstcol; ++l)
                  tmp(k,l) = (*this)(k,l);
              for(enumerator l = lastcol+1; l < m; ++l)
                  tmp(k,l-shiftcol-1) = (*this)(k,l);
          }
          for(enumerator k = lastrow+1; k < n; ++k)
          {
              for(enumerator l = 0; l < firstcol; ++l)
                  tmp(k-shiftrow-1,l) = (*this)(k,l);
              for(enumerator l = lastcol+1; l < m; ++l)
                  tmp(k-shiftrow-1,l-shiftcol-1) = (*this)(k,l);
          }
          this->Swap(tmp);
      }
    void Swap(Matrix & b)
    {
      space.swap(b.space);
      std::swap(n,b.n);
      std::swap(m,b.m);
    }
    /// Singular value decomposition.
    /// Reconstruct matrix: A = U*Sigma*V.Transpose().
    /// @param U Left unitary matrix, U^T U = I.
    /// @param Sigma Diagonal matrix with singular values.
    /// @param V Right unitary matrix, not transposed.
    /// @param order_singular_values
    bool SVD(Matrix & U, Matrix & Sigma, Matrix & V, bool order_singular_values = true)
    {
      int flag, i, its, j, jj, k, l, nm;
      Var c, f, h, s, x, y, z;
      Var g = 0.0, scale = 0.0;
      INMOST_DATA_REAL_TYPE anorm = 0.0;
      if (n < m) 
      {
        bool success = Transpose().SVD(U,Sigma,V);
        if( success )
        {
          U.Swap(V);
          U = U.Transpose();
          V = V.Transpose();
          return true;
        }
        else return false;
      } //m <= n
      array<Var> _rv1(m);
      shell<Var> rv1(_rv1);
      U = (*this);
      Sigma.Resize(m,m);
      Sigma.Zero();
      V.Resize(m,m);
  
      std::swap(n,m); //this how original algorithm takes it
      // Householder reduction to bidiagonal form
      for (i = 0; i < (int)n; i++) 
      {
        // left-hand reduction
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < (int)m) 
        {
          for (k = i; k < (int)m; k++) scale += fabs(U(k,i));
          if (get_value(scale)) 
          {
            for (k = i; k < (int)m; k++) 
            {
              U(k,i) /= scale;
              s += U(k,i) * U(k,i);
            }
            f = U(i,i);
            g = -sign_func(sqrt(s), f);
            h = f * g - s;
            U(i,i) = f - g;
            if (i != n - 1) 
            {
              for (j = l; j < (int)n; j++) 
              {
                for (s = 0.0, k = i; k < (int)m; k++) s += U(k,i) * U(k,j);
                f = s / h;
                for (k = i; k < (int)m; k++) U(k,j) += f * U(k,i);
              }
            }
            for (k = i; k < (int)m; k++) U(k,i) *= scale;
          }
        }
        Sigma(i,i) = scale * g;
        // right-hand reduction
        g = s = scale = 0.0;
        if (i < (int)m && i != n - 1) 
        {
          for (k = l; k < (int)n; k++) scale += fabs(U(i,k));
          if (get_value(scale)) 
          {
            for (k = l; k < (int)n; k++) 
            {
              U(i,k) = U(i,k)/scale;
              s += U(i,k) * U(i,k);
            }
            f = U(i,l);
            g = -sign_func(sqrt(s), f);
            h = f * g - s;
            U(i,l) = f - g;
            for (k = l; k < (int)n; k++) rv1[k] = U(i,k) / h;
            if (i != m - 1) 
            {
              for (j = l; j < (int)m; j++) 
              {
                for (s = 0.0, k = l; k < (int)n; k++) s += U(j,k) * U(i,k);
                for (k = l; k < (int)n; k++) U(j,k) += s * rv1[k];
              }
            }
            for (k = l; k < (int)n; k++) U(i,k) *= scale;
          }
        }
        anorm = max_func(anorm,fabs(get_value(Sigma(i,i))) + fabs(get_value(rv1[i])));
      }

      // accumulate the right-hand transformation
      for (i = n - 1; i >= 0; i--) 
      {
        if (i < (int)(n - 1)) 
        {
          if (get_value(g)) 
          {
            for (j = l; j < (int)n; j++) V(j,i) = ((U(i,j) / U(i,l)) / g);
            // double division to avoid underflow
            for (j = l; j < (int)n; j++) 
            {
              for (s = 0.0, k = l; k < (int)n; k++) s += U(i,k) * V(k,j);
              for (k = l; k < (int)n; k++) V(k,j) += s * V(k,i);
            }
          }
          for (j = l; j < (int)n; j++) V(i,j) = V(j,i) = 0.0;
        }
        V(i,i) = 1.0;
        g = rv1[i];
        l = i;
      }

      // accumulate the left-hand transformation
      for (i = n - 1; i >= 0; i--) 
      {
        l = i + 1;
        g = Sigma(i,i);
        if (i < (int)(n - 1)) 
          for (j = l; j < (int)n; j++) 
            U(i,j) = 0.0;
        if (get_value(g)) 
        {
          g = 1.0 / g;
          if (i != n - 1) 
          {
            for (j = l; j < (int)n; j++) 
            {
              for (s = 0.0, k = l; k < (int)m; k++) s += (U(k,i) * U(k,j));
              f = (s / U(i,i)) * g;
              for (k = i; k < (int)m; k++) U(k,j) += f * U(k,i);
            }
          }
          for (j = i; j < (int)m; j++) U(j,i) = U(j,i)*g;
        }
        else for (j = i; j < (int)m; j++) U(j,i) = 0.0;
        U(i,i) += 1;
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
            if (fabs(get_value(rv1[l])) + anorm == anorm) 
            {
              flag = 0;
              break;
            }
            if (fabs(get_value(Sigma(nm,nm))) + anorm == anorm) 
              break;
          }
          if (flag) 
          {
            c = 0.0;
            s = 1.0;
            for (i = l; i <= k; i++) 
            {
              f = s * rv1[i];
              if (fabs(get_value(f)) + anorm != anorm) 
              {
                g = Sigma(i,i);
                h = pythag(f, g);
                Sigma(i,i) = h; 
                h = 1.0 / h;
                c = g * h;
                s = (- f * h);
                for (j = 0; j < (int)m; j++) 
                {
                  y = U(j,nm);
                  z = U(j,i);
                  U(j,nm) = (y * c + z * s);
                  U(j,i) = (z * c - y * s);
                }
              }
            }
          }
          z = Sigma(k,k);
          if (l == k) 
          {// convergence
            if (z < 0.0) 
            {// make singular value nonnegative
              Sigma(k,k) = -z;
              for (j = 0; j < (int)n; j++) V(j,k) = -V(j,k);
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
          x = Sigma(l,l);
          nm = k - 1;
          y = Sigma(nm,nm);
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
            y = Sigma(i,i);
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
            for (jj = 0; jj < (int)n; jj++) 
            {
              x = V(jj,j);
              z = V(jj,i);
              V(jj,j) = (x * c + z * s);
              V(jj,i) = (z * c - x * s);
            }
            z = pythag(f, h);
            Sigma(j,j) = z;
            if (z) 
            {
              z = 1.0 / z;
              c = f * z;
              s = h * z;
            }
            f = (c * g) + (s * y);
            x = (c * y) - (s * g);
            for (jj = 0; jj < (int)m; jj++) 
            {
              y = U(jj,j);
              z = U(jj,i);
              U(jj,j) = (y * c + z * s);
              U(jj,i) = (z * c - y * s);
            }
          }
          rv1[l] = 0.0;
          rv1[k] = f;
          Sigma(k,k) = x;
        }
      }
      //CHECK THIS!
      if( order_singular_values )
      {
        for(i = 0; i < (int)n; i++)
		    {
			    k = i;
			    for(j = i+1; j < (int)n; ++j)
				    if( Sigma(k,k) < Sigma(j,j) ) k = j;
			    Var temp;
			    if( Sigma(k,k) > Sigma(i,i) )
			    {
				    temp       = Sigma(k,k);
				    Sigma(k,k) = Sigma(i,i);
				    Sigma(i,i) = temp;
            // U is m by n
            for(int j = 0; j < (int)m; ++j)
				    {
					    temp   = U(j,k);
					    U(j,k) = U(j,i);
					    U(j,i) = temp;
            }
            // V is n by n
            for(int j = 0; j < (int)n; ++j)
            {
					    temp   = V(j,k);
					    V(j,k) = V(j,i);
					    V(j,i) = temp;
				    }
			    }
		    }
      }
      
      std::swap(n,m);
      return true;
    }
    Matrix() : space(0),n(0),m(0) {}
    Matrix(Var * pspace, enumerator pn, enumerator pm) : space(pspace,pspace+pn*pm), n(pn), m(pm) {}
    Matrix(enumerator pn, enumerator pm) : space(pn*pm), n(pn), m(pm) {}
    Matrix(const Matrix & other) : space(other.n*other.m), n(other.n), m(other.m) 
    {
      for(enumerator i = 0; i < n*m; ++i)
        space[i] = other.space[i];
    }
      template<typename typeB>
      Matrix(const Matrix<typeB> & other) : space(other.Cols()*other.Rows()), n(other.Rows()), m(other.Cols())
      {
          for(enumerator i = 0; i < n; ++i)
              for(enumerator j = 0; j < m; ++j)
                  (*this)(i,j) = get_value(other(i,j));
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
      if( n*m != other.n*other.m ) space.resize(other.n*other.m);
      for(enumerator i = 0; i < other.n*other.m; ++i)
        space[i] = other.space[i];
      n = other.n;
      m = other.m;
      return *this;
    }
      template<typename typeB>
      Matrix & operator =(Matrix<typeB> const & other)
      {
          if( n*m != other.n*other.m ) space.resize(other.n*other.m);
          for(enumerator i = 0; i < other.n*other.m; ++i)
              space[i] = get_value(other.space[i]);
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
    const Var & operator()(enumerator i, enumerator j) const
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
        for(enumerator i = 0; i < Rows(); ++i)
            for(enumerator j = 0; j < Cols(); ++j)
                ret(i,j) = (*this)(i,j)-other(i,j);
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
      for(enumerator i = 0; i < Rows(); ++i)
          for(enumerator j = 0; j < Cols(); ++j)
              ret(i,j) = (*this)(i,j)+other(i,j);
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
      for(enumerator i = 0; i < Rows(); ++i)
        for(enumerator j = 0; j < Cols(); ++j) ret(i,j) = (*this)(i,j)*coef;
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
        for(enumerator i = 0; i < Rows(); ++i)
            for(enumerator j = 0; j < Cols(); ++j) ret(i,j) = (*this)(i,j)/coef;
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
    std::pair<Matrix,bool> Invert(bool print_fail = false) const
    {
      std::pair<Matrix,bool> ret = std::make_pair(Matrix(m,n),true);
      Matrix At = Transpose(); //m by n matrix
      Matrix AtB = At; //m by n matrix
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
					if( print_fail ) std::cout << "Failed to invert matrix" << std::endl;
					ret.second = false;
					delete [] order;
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
#if defined(USE_AUTODIFF)
            std::cout << std::setw(10) << get_value((*this)(k,l));
#else
            std::cout << std::setw(10) << (*this)(k,l);
#endif
          else
            std::cout << std::setw(10) << 0;
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
    template<typename typeB>
    typename Promote<Var,typeB>::type DotProduct(const Matrix<typeB> & other) const
    {
      assert(n == other.n);
      assert(m == other.m);
      typename Promote<Var,typeB>::type ret = 0.0;
      for(enumerator i = 0; i < n; ++i)
        for(enumerator j = 0; j < m; ++j)
          ret += ((*this)(i,j))*other(i,j);
      return ret;
    }
      template<typename typeB>
      typename Promote<Var,typeB>::type operator ^(const Matrix<typeB> & other) const
      {
          return DotProduct(other);
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
    ///Retrive vector in matrix form from array
    static Matrix<Var> FromVector(Var * n, enumerator size)
    {
      return Matrix(n,size,1);
    }
    ///Create diagonal matrix from array
    static Matrix<Var> FromDiagonal(Var * r, enumerator size)
    {
      Matrix ret(size,size);
      ret.Zero();
      for(enumerator k = 0; k < size; ++k) ret(k,k) = r[k];
      return ret;
    }
    ///Create diagonal matrix from array of inversed values
    static Matrix<Var> FromDiagonalInverse(Var * r, enumerator size)
    {
      Matrix ret(size,size);
      ret.Zero();
      for(enumerator k = 0; k < size; ++k) ret(k,k) = 1.0/r[k];
      return ret;
    }
	static Matrix CrossProduct(Var vec[3])
	{
		// |  0  -z   y |
		// |  z   0  -x |
		// | -y   x   0 |
		Matrix ret(3,3);
		ret(0,0) = 0.0;
		ret(0,1) = -vec[2]; //-z
		ret(0,2) = vec[1]; //y
		ret(1,0) = vec[2]; //z
		ret(1,1) = 0;
		ret(1,2) = -vec[0]; //-x
		ret(2,0) = -vec[1]; //-y
		ret(2,1) = vec[0]; //x
		ret(2,2) = 0;
		return ret;
	}
    ///Unit matrix
    static Matrix Unit(enumerator pn)
    {
      Matrix ret(pn,pn);
      ret.Zero();
      for(enumerator i = 0; i < pn; ++i) ret(i,i) = 1.0;
      return ret;
    }
    /// Concatenate B matrix as columns of current matrix.
    /// Assumes that number of rows of current matrix is
    /// equal to number of rows of B matrix.
    Matrix ConcatCols(const Matrix & B)
    {
        assert(Rows() == B.Rows());
        Matrix ret(Rows(),Cols()+B.Cols());
        Matrix & A = *this;
        for(enumerator i = 0; i < Rows(); ++i)
        {
            for(enumerator j = 0; j < Cols(); ++j)
                ret(i,j) = A(i,j);
            for(enumerator j = 0; j < B.Cols(); ++j)
                ret(i,j+Cols()) = B(i,j);
        }
        return ret;
    }
    /// Concatenate B matrix as rows of current matrix.
    /// Assumes that number of colums of current matrix is
    /// equal to number of columns of B matrix.
    Matrix ConcatRows(const Matrix & B)
    {
        assert(Cols() == B.Cols());
        Matrix ret(Rows()+B.Rows(),Cols());
        Matrix & A = *this;
        for(enumerator i = 0; i < Rows(); ++i)
        {
            for(enumerator j = 0; j < Cols(); ++j)
                ret(i,j) = A(i,j);
        }
        for(enumerator i = 0; i < B.Rows(); ++i)
        {
            for(enumerator j = 0; j < Cols(); ++j)
                ret(i+Rows(),j) = B(i,j);
        }
        return ret;
    }

    /// Joint diagonalization algorithm by Cardoso
    /// http://perso.telecom-paristech.fr/~cardoso/Algo/Joint_Diag/joint_diag_r.m
    /// Current matrix should have size n by n*m
    /// And represent concatination of m n by n matrices
    Matrix JointDiagonalization(INMOST_DATA_REAL_TYPE threshold = 1.0e-7)
    {
        enumerator N = Rows();
        enumerator M = Cols() / Rows();
        Matrix V = Matrix::Unit(m);
        Matrix R(2,M);
        Matrix G(2,2);
        Matrix & A = *this;
        Var ton, toff, theta, c, s, Ap, Aq, Vp, Vq;
        bool repeat;
        do
        {
            repeat = false;
            for(enumerator p = 0; p < N-1; ++p)
            {
                for(enumerator q = p+1; q < N; ++q)
                {
                    for(enumerator k = 0; k < M; ++k)
                    {
                        R(0,k) = A(p,p + k*N) - A(q,q + k*N);
                        R(1,k) = A(p,q + k*N) + A(q,p + k*N);
                    }
                    G = R*R.Transpose();
                    Var ton  = G(0,0) - G(1,1);
                    Var toff = G(0,1) + G(1,0);
                    Var theta = 0.5 * atan2( toff, ton + sqrt(ton*ton + toff*toff) );
                    Var c = cos(theta);
                    Var s = sin(theta);
                    if( fabs(s) > threshold )
                    {
                        //std::cout << "p,q: " << p << "," << q << " c,s: " << c << "," << s << std::endl;
                        repeat = true;
                        for(enumerator k = 0; k < M; ++k)
                        {
                            for(enumerator i = 0; i < N; ++i)
                            {
                                Ap = A(i,p + k*N);
                                Aq = A(i,q + k*N);
                                A(i,p + k*N) = Ap*c + Aq*s;
                                A(i,q + k*N) = Aq*c - Ap*s;
                            }
                        }
                        for(enumerator k = 0; k < M; ++k)
                        {
                            for(enumerator j = 0; j < N; ++j)
                            {
                                Ap = A(p,j + k*N);
                                Aq = A(q,j + k*N);
                                A(p,j + k*N) = Ap*c + Aq*s;
                                A(q,j + k*N) = Aq*c - Ap*s;
                            }
                        }
                        for(enumerator i = 0; i < N; ++i)
                        {
                            Vp = V(i,p);
                            Vq = V(i,q);
                            V(i,p) = Vp*c + Vq*s;
                            V(i,q) = Vq*c - Vp*s;
                        }
                    }
                }
            }
            //Print();
        } while( repeat );
        return V;
    }

  };
    
	typedef Matrix<INMOST_DATA_REAL_TYPE> rMatrix; //shortcut for real matrix
#if defined(USE_AUTODIFF)
	typedef Matrix<variable> vMatrix; //shortcut for matrix with variations
#endif
    
}

template<typename typeB>
INMOST::Matrix<typename INMOST::Promote<INMOST_DATA_REAL_TYPE,typeB>::type> operator *(INMOST_DATA_REAL_TYPE coef, const INMOST::Matrix<typeB> & other)
{return other*coef;}
	
#if defined(USE_AUTODIFF)
template<typename typeB>
INMOST::Matrix<typename INMOST::Promote<INMOST::variable,typeB>::type> operator *(const INMOST::variable & coef, const INMOST::Matrix<typeB> & other)
{return other*coef;}
#endif


#endif //INMOST_DENSE_INCLUDED