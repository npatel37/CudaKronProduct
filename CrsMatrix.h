/*! \file CrsMatrix.h
 *
 *  A class to represent a sparse matrix in Compressed Row Storage
 *
 */

#ifndef CRSMATRIX_HEADER_H
#define CRSMATRIX_HEADER_H
#include <algorithm>
#include "Matrix.h"
#include <cassert>
#include <stdexcept>
#include <vector>
#include <complex>


//! A Sparse Matrix in Compressed Row Storage (CRS) format.
/**
	The CRS format puts the subsequent nonzero elements of the matrix rows
	in contiguous memory locations. We create 3 vectors: one for complex numbers
	containing the values of the
	matrix entries
	and the other two for integers ($colind$ and $rowptr$).
	The vector $values$ stores the values of the non-zero elements of the matrix,
	as they are traversed in a row-wise fashion.
	The $colind$ vector stores the column indices of the elements of the $values$
	vector. That is, if $values[k] = a[i][j]$ then $colind[k] = j$.
	The $rowptr$ vector stores the locations in the $values$ vector that start
	a row, that is $values[k] = a[i][j]$ if $rowptr[i] \le i < rowptr[i + 1]$.
	By convention, we define $rowptr[N_{dim}]$ to be equal to the number of non-zero elements,
	$n_z$, in the matrix. The storage savings of this approach are significant
	because instead of
	storing $N_{dim}^2$ elements, we need only $2n_z + N_{dim} + 1$ storage locations.\\
	To illustrate how the CRS format works, consider the non-symmetric matrix defined by
	\begin{equation}
		A=\left[\begin{tabular}{llllll}

		10 &  0 & 0 & 0  & -2 & 0 \\
		3 &  9 &  0 &  0 &  0 &  3 \\
		0 &  7 &  8 &  7 &  0 &  0 \\
		3 &  0 &  8 &  7  & 5 &  0 \\
		0 &   8 &  0 &  9 &  9 & 13 \\
		0 &  4 &  0 &  0 &  2&  -1 \\
	\end{tabular}\right]\end{equation}
	The CRS format for this matrix is then specified by the arrays:\\
	\begin{tt}
		values = [10 -2  3  9  3  7  8  7  3 ... 9 13  4  2 -1 ]\\
		colind = [ 0  4  0  1  5  1  2  3  0 ... 4  5  1  4  5 ]\\
		rowptr = [ 0  2  5  8 12 16 19 ]\\
	\end{tt}
	*/
template<class T>
class CrsMatrix {

public:

	typedef T MatrixElementType;
	typedef T value_type;

	CrsMatrix() 
		: nrow_(0),ncol_(0) 
		{ }

	~CrsMatrix() {  }

	CrsMatrix(size_t nrow,size_t ncol)
	    : nrow_(nrow),ncol_(ncol)
	{
		resize(nrow,ncol);
	}

	template<typename S>
	CrsMatrix(const CrsMatrix<S>& a)
	{
		colind_=a.colind_;
		rowptr_=a.rowptr_;
		values_=a.values_;
		nrow_ = a.nrow_;
		ncol_ = a.ncol_;
	}

	template<typename S>
	CrsMatrix(const CrsMatrix<std::complex<S> >& a)
	{
		colind_=a.colind_;
		rowptr_=a.rowptr_;
		values_=a.values_;
		nrow_ = a.nrow_;
		ncol_ = a.ncol_;
	}

	explicit CrsMatrix(const Matrix<T>& a)
	{
		int counter=0;
		double eps = 0;

		resize(a.n_row(),a.n_col());

		for (size_t i = 0; i < a.n_row(); i++) {
			setRow(i,counter);
			for (size_t j=0;j<a.n_col();j++) {
				if (std::norm(a(i,j))<=eps) continue;
				pushValue(a(i,j));
				pushCol(j);
				counter++;
			}

		}
		setRow(a.n_row(),counter);
	}

	void resize(size_t nrow,size_t ncol)
	{
		colind_.clear();
		values_.clear();
		rowptr_.clear();
		rowptr_.resize(nrow+1);
		nrow_ = nrow;
		ncol_ = ncol;
	}

	void clear()
	{
		colind_.clear();
		values_.clear();
		rowptr_.clear();
		nrow_=ncol_=0;
	}

	void resize(size_t nrow,size_t ncol,size_t nonzero)
	{
		resize(nrow,ncol);
		colind_.resize(nonzero);
		values_.resize(nonzero);
	}
        	
	void resizecv(size_t nonzero)
	{
		colind_.resize(nonzero);
		values_.resize(nonzero);
	}

	void setRow(size_t n,size_t v)
	{
		assert(n<rowptr_.size());
		rowptr_[n]=v;
	}

	void setCol(int n,int v) {
		colind_[n]=v;
	}

	void setValues(int n,const T &v) {
		values_[n]=v;
	}

	void operator*=(T x)
	{
		values_ *= x;
	}

	void operator+=(const CrsMatrix& m)
	{
		CrsMatrix c;
		add(c,m,1.0);
		*this = c;
	}

	T element(int i,int j) const
	{
		for (int k=rowptr_[i];k<rowptr_[i+1];k++) {
			if (colind_[k]==j) {
				return values_[k];
			}
		}
		return static_cast<T>(0.0);
	}

	int nonZero() const { 
		return colind_.size(); 
	}

	/** performs x = x + A * y
		 ** where x and y are vectors and A is a sparse matrix in
		 ** row-compressed format */
	template<typename VectorLikeType>
	void matrixVectorProduct(VectorLikeType& x, const VectorLikeType& y) const
	{
		assert(x.size()==y.size());
		for (size_t i = 0; i < y.size(); i++) {
			assert(i+1<rowptr_.size());
			for (int j = rowptr_[i]; j < rowptr_[i + 1]; j++) {
				assert(size_t(j)<values_.size());
				assert(size_t(j)<colind_.size());
				assert(size_t(colind_[j])<y.size());
				x[i] += values_[j] * y[colind_[j]];
			}
		}
	}

	size_t row() const { return nrow_; }

	size_t col() const { return ncol_; }

	size_t rank() const
	{
		if (nrow_!=ncol_)
			throw std::runtime_error("CrsMatrix: rank(): only for square matrices\n");
		return nrow_;
	}

	void pushCol(size_t i) { 
		colind_.push_back(i); 
	}

	void pushValue(T const &value) { 
		values_.push_back(value); 
	}

	//! Make a diagonal CRS matrix with value "value"
	void makeDiagonal(size_t row,T const &value=0)
	{
		nrow_=row;
		ncol_=row;
		rowptr_.resize(row+1);
		values_.resize(row);
		colind_.resize(row);

		for (size_t i=0;i<row;i++) {
			values_[i]=value;
			colind_[i]=i;
			rowptr_[i]=i;
		}
		rowptr_[row]=row;
	}

	const int& getRowPtr(size_t i) const { 
		assert(i<rowptr_.size()); 
		return rowptr_[i]; 
	}

	const int& getCol(size_t i) const { 
		assert(i<colind_.size()); 
		return colind_[i]; 
	}

	const T& getValue(size_t i) const { 
		assert(i<values_.size()); 
		return values_[i];
	}

	Matrix<T> toDense() const
	{
		Matrix<T> m;
		crsMatrixToFullMatrix(m,*this);
		return m;
	}

	void checkValidity() const
	{
		size_t n = nrow_;
		assert(n+1==rowptr_.size());
		assert(nrow_>0 && ncol_>0);
		for (size_t i=0;i<n;i++) {
			typename std::vector<size_t> p(ncol_,0);
			for (int k=rowptr_[i];k<rowptr_[i+1];k++) {
				size_t col = colind_[k];
				assert(p[col]==0);
				p[col] = 1;
			}
		}
	}

    void DensePrint() const
    {
      Matrix<T> m;
      crsMatrixToFullMatrix(m,*this);
      m.print();
    }

	//serializr start class CrsMatrix
	typename std::vector<int> rowptr_;	//serializr normal rowptr_
	typename std::vector<int> colind_;	//serializr normal colind_
	typename std::vector<T> values_;	//serializr normal values_
	size_t nrow_;	//serializr normal nrow_
	size_t ncol_;	//serializr normal ncol_

private:

	template<typename T1>
	void add(CrsMatrix<T>& c, const CrsMatrix<T>& m, const T1& t1) const
	{
		assert(m.row()==m.col());

		const T1 one = 1.0;
        if (nrow_>=m.row()) {
            //c=this*one+m*t1
            operatorPlus(c,*this,one,m,t1);
        } else {
            //c=m*t1+this*one
            operatorPlus(c,m,t1,*this,one);
        }
	}



}; // class CrsMatrix


//! Transforms a Compressed-Row-Storage (CRS) into a full Matrix (Fast version)
template<typename T>
void crsMatrixToFullMatrix(Matrix<T>& m,const CrsMatrix<T>& crsMatrix)
{
	m.reset(crsMatrix.row(),crsMatrix.col());
	for (size_t i = 0; i < crsMatrix.row() ; i++) {
		for (size_t k=0;k<crsMatrix.row();k++) m(i,k)=0;
		for (int k=crsMatrix.getRowPtr(i);k<crsMatrix.getRowPtr(i+1);k++)
			m(i,crsMatrix.getCol(k))=crsMatrix.getValue(k);
	}

}


//! Transforms a full matrix into a Compressed-Row-Storage (CRS) Matrix
// Use the constructor if possible
template<typename T>
void fullMatrixToCrsMatrix(CrsMatrix<T>& crsMatrix, const Matrix<T>& a)
{
	crsMatrix.resize(a.n_row(),a.n_col());

	size_t counter = 0;
	for (size_t i = 0; i < a.n_row(); i++) {
		crsMatrix.setRow(i,counter);
		for (size_t j=0;j<a.n_col();j++) {
			if (a(i,j)==static_cast<T>(0)) continue;
			crsMatrix.pushValue(a(i,j));
			crsMatrix.pushCol(j);
			counter++;
		}

	}
	crsMatrix.setRow(crsMatrix.row(),counter);
	crsMatrix.checkValidity();
}


//! Sets A=B*b1+C*c1, restriction: B.size has to be larger or equal than C.size
template<typename T, typename T1>
void operatorPlus(CrsMatrix<T>& A,
                  const CrsMatrix<T>& B,
                  T1& b1,
                  const CrsMatrix<T>& C,
                  T1& c1)
{
    size_t n = B.row();
    assert(B.row()==B.col());
    assert(C.row()==C.col());

    T tmp;

    assert(n>=C.row());

    typename std::vector<T>  valueTmp(n);
    typename std::vector<int> index;
    A.resize(n,B.col());

    size_t counter=0;
    for (size_t k2=0;k2<n;k2++) valueTmp[k2]= static_cast<T>(0.0);

    for (size_t i = 0; i < n; i++) {
        int k;
        A.setRow(i,counter);

        if (i<C.row()) {
            // inspect this
            index.clear();
            for (k=B.getRowPtr(i);k<B.getRowPtr(i+1);k++) {
                if (B.getCol(k)<0 || size_t(B.getCol(k))>=n)
                    throw std::runtime_error("operatorPlus (1)\n");
                valueTmp[B.getCol(k)]=B.getValue(k)*b1;
                index.push_back(B.getCol(k));
            }

            // inspect C
            for (k=C.getRowPtr(i);k<C.getRowPtr(i+1);k++) {
                tmp = C.getValue(k)*c1;
                if (C.getCol(k)>=int(valueTmp.size()) || C.getCol(k)<0)
                    throw std::runtime_error("operatorPlus (2)\n");

                valueTmp[C.getCol(k)] += tmp;
                index.push_back(C.getCol(k));
            }
            std::sort(index.begin(),index.end());
            k= -1;
            for (size_t kk=0;kk<index.size();kk++) {
                if (k==index[kk]) continue;
                k=index[kk];
                if (k<0 || size_t(k)>=n)
                    throw std::runtime_error("operatorPlus (3)\n");
                tmp = valueTmp[k];
                if (tmp!=static_cast<T>(0.0)) {
                    A.pushCol(k);
                    A.pushValue(tmp);
                    counter++;
                    valueTmp[k]=static_cast<T>(0.0);
                }
            }
        } else {
            for (k=B.getRowPtr(i);k<B.getRowPtr(i+1);k++) {
                tmp = B.getValue(k)*b1;
                A.pushCol(B.getCol(k));
                A.pushValue(tmp);
                counter++;
            }
        }
    }

    A.setRow(n,counter);
}



#endif

