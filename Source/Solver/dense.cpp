#define _CRT_SECURE_NO_WARNINGS
#include "inmost_dense.h"

namespace INMOST
{
	
	/*
	template<>
	template<>
	Matrix<typename Promote<INMOST_DATA_REAL_TYPE,variable>::type, pool_array_t<typename Promote<INMOST_DATA_REAL_TYPE,variable>::type> >
	AbstractMatrix<INMOST_DATA_REAL_TYPE>::operator*<variable>(const AbstractMatrix<variable> & other) const
	{
		std::cout << __FUNCTION__ << std::endl;
		assert(Cols() == other.Rows());
		Matrix<typename Promote<INMOST_DATA_REAL_TYPE,variable>::type, pool_array_t<typename Promote<INMOST_DATA_REAL_TYPE,variable>::type> > ret(Rows(),other.Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i) //loop rows
		{
			for(enumerator j = 0; j < other.Cols(); ++j) //loop columns
			{
				if( CheckCurrentAutomatizator() )
				{
					Sparse::RowMerger & merger = GetCurrentMerger();
					double value = 0.0;
					for(enumerator k = 0; k < Cols(); ++k)
					{
						value += (*this)(i,k)*other(k,j).GetValue();
						merger.AddRow((*this)(i,k),other(k,j).GetRow());
					}
					ret(i,j).SetValue(value);
					merger.RetriveRow(ret(i,j).GetRow());
					merger.Clear();
				}
				else
				{
					//~ typename Promote<INMOST_DATA_REAL_TYPE,variable>::type tmp = 0.0;
#pragma unroll
					for(enumerator k = 0; k < Cols(); ++k)
						ret(i,j) += (*this)(i,k)*other(k,j);
					//~ ret(i,j) = tmp;
				}
			}
		}
		return ret;
	}
	
	template<>
	template<>
	Matrix<typename Promote<variable,INMOST_DATA_REAL_TYPE>::type, pool_array_t<typename Promote<variable,INMOST_DATA_REAL_TYPE>::type> >
	AbstractMatrix<variable>::operator*<INMOST_DATA_REAL_TYPE>(const AbstractMatrix<INMOST_DATA_REAL_TYPE> & other) const
	{
		std::cout << __FUNCTION__ << std::endl;
		assert(Cols() == other.Rows());
		Matrix<typename Promote<variable,INMOST_DATA_REAL_TYPE>::type, pool_array_t<typename Promote<variable,INMOST_DATA_REAL_TYPE>::type> > ret(Rows(),other.Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i) //loop rows
		{
			for(enumerator j = 0; j < other.Cols(); ++j) //loop columns
			{
				if( CheckCurrentAutomatizator() )
				{
					Sparse::RowMerger & merger = GetCurrentMerger();
					double value = 0.0;
					for(enumerator k = 0; k < Cols(); ++k)
					{
						value += (*this)(i,k).GetValue()*other(k,j);
						merger.AddRow(other(k,j),(*this)(i,k).GetRow());
					}
					ret(i,j).SetValue(value);
					merger.RetriveRow(ret(i,j).GetRow());
					merger.Clear();
				}
				else
				{
					//~ typename Promote<INMOST_DATA_REAL_TYPE,variable>::type tmp = 0.0;
#pragma unroll
					for(enumerator k = 0; k < Cols(); ++k)
						ret(i,j) += (*this)(i,k)*other(k,j);
					//~ ret(i,j) = tmp;
				}
			}
		}
		return ret;
	}
	
	template<>
	template<>
	Matrix<typename Promote<variable,variable>::type, pool_array_t<typename Promote<variable,variable>::type> >
	AbstractMatrix<variable>::operator*<variable>(const AbstractMatrix<variable> & other) const
	{
		std::cout << __FUNCTION__ << std::endl;
		assert(Cols() == other.Rows());
		Matrix<typename Promote<variable,variable>::type, pool_array_t<typename Promote<variable,variable>::type> > ret(Rows(),other.Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i) //loop rows
		{
			for(enumerator j = 0; j < other.Cols(); ++j) //loop columns
			{
				if( CheckCurrentAutomatizator() )
				{
					Sparse::RowMerger & merger = GetCurrentMerger();
					double value = 0.0;
					for(enumerator k = 0; k < Cols(); ++k)
					{
						value += (*this)(i,k).GetValue()*other(k,j).GetValue();
						merger.AddRow(other(k,j).GetValue(),(*this)(i,k).GetRow());
						merger.AddRow((*this)(i,k).GetValue(),other(k,j).GetRow());
					}
					ret(i,j).SetValue(value);
					merger.RetriveRow(ret(i,j).GetRow());
					merger.Clear();
				}
				else
				{
					//~ typename Promote<INMOST_DATA_REAL_TYPE,variable>::type tmp = 0.0;
#pragma unroll
					for(enumerator k = 0; k < Cols(); ++k)
						ret(i,j) += (*this)(i,k)*other(k,j);
					//~ ret(i,j) = tmp;
				}
			}
		}
		return ret;
	}
	extern template class Matrix<typename Promote<INMOST_DATA_REAL_TYPE,variable>::type, pool_array_t<typename Promote<INMOST_DATA_REAL_TYPE,variable>::type> > AbstractMatrix<INMOST_DATA_REAL_TYPE>::operator*<variable>(const AbstractMatrix<variable> & other) const;
	extern template class Matrix<typename Promote<variable,INMOST_DATA_REAL_TYPE>::type, pool_array_t<typename Promote<variable,INMOST_DATA_REAL_TYPE>::type> > AbstractMatrix<variable>::operator*<INMOST_DATA_REAL_TYPE>(const AbstractMatrix<INMOST_DATA_REAL_TYPE> & other) const;
	extern template class Matrix<typename Promote<variable,variable>::type, pool_array_t<typename Promote<variable,variable>::type> > AbstractMatrix<variable>::operator*<variable>(const AbstractMatrix<variable> & other) const;
	*/
	//~ Matrix<typename Promote<INMOST_DATA_REAL_TYPE,variable>::type, pool_array_t<typename Promote<INMOST_DATA_REAL_TYPE,variable>::type> >
	//~ AbstractMatrix<INMOST_DATA_REAL_TYPE>::operator*<variable>(const AbstractMatrix<variable> & other) const;
	
	//~ Matrix<typename Promote<variable,INMOST_DATA_REAL_TYPE>::type, pool_array_t<typename Promote<variable,INMOST_DATA_REAL_TYPE>::type> >
	//~ AbstractMatrix<variable>::operator*<INMOST_DATA_REAL_TYPE>(const AbstractMatrix<INMOST_DATA_REAL_TYPE> & other) const
	
	//~ Matrix<typename Promote<variable,variable>::type, pool_array_t<typename Promote<variable,variable>::type> >
	//~ AbstractMatrix<variable>::operator*<variable>(const AbstractMatrix<variable> & other) const;
}
