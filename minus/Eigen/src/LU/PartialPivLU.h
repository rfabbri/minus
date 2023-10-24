// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2006-2009 Benoit Jacob <jacob.benoit.1@gmail.com>
// Copyright (C) 2009 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//------------------------------------------------------------------------------
// Largely simplified by Fabbri
//

#ifndef EIGEN_PARTIALLU_H
#define EIGEN_PARTIALLU_H

namespace Eigen {

namespace internal {
template<typename _MatrixType> struct traits<PartialPivLU<_MatrixType> >
 : traits<_MatrixType>
{
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef int StorageIndex;
  typedef traits<_MatrixType> BaseTraits;
  enum {
    Flags = BaseTraits::Flags & RowMajorBit,
    CoeffReadCost = Dynamic
  };
};

template<typename T,typename Derived>
struct enable_if_ref;
// {
//   typedef Derived type;
// };

template<typename T,typename Derived>
struct enable_if_ref<Ref<T>,Derived> {
  typedef Derived type;
};

} // end namespace internal

/** \ingroup LU_Module
  *
  * \class PartialPivLU
  *
  * \brief LU decomposition of a matrix with partial pivoting, and related features
  *
  * \tparam _MatrixType the type of the matrix of which we are computing the LU decomposition
  *
  * This class represents a LU decomposition of a \b square \b invertible matrix, with partial pivoting: the matrix A
  * is decomposed as A = PLU where L is unit-lower-triangular, U is upper-triangular, and P
  * is a permutation matrix.
  *
  * Typically, partial pivoting LU decomposition is only considered numerically stable for square invertible
  * matrices. Thus LAPACK's dgesv and dgesvx require the matrix to be square and invertible. The present class
  * does the same. It will assert that the matrix is square, but it won't (actually it can't) check that the
  * matrix is invertible: it is your task to check that you only use this decomposition on invertible matrices.
  *
  * The guaranteed safe alternative, working for all matrices, is the full pivoting LU decomposition, provided
  * by class FullPivLU.
  *
  * This is \b not a rank-revealing LU decomposition. Many features are intentionally absent from this class,
  * such as rank computation. If you need these features, use class FullPivLU.
  *
  * This LU decomposition is suitable to invert invertible matrices. It is what MatrixBase::inverse() uses
  * in the general case.
  * On the other hand, it is \b not suitable to determine whether a given matrix is invertible.
  *
  * The data of the LU decomposition can be directly accessed through the methods matrixLU(), permutationP().
  *
  * This class supports the \link InplaceDecomposition inplace decomposition \endlink mechanism.
  * 
  * \sa MatrixBase::partialPivLu(), MatrixBase::determinant(), MatrixBase::inverse(), MatrixBase::computeInverse(), class FullPivLU
  */
template<typename _MatrixType> class PartialPivLU
  : public SolverBase<PartialPivLU<_MatrixType> >
{
  public:

    typedef _MatrixType MatrixType;
    typedef SolverBase<PartialPivLU> Base;
    friend class SolverBase<PartialPivLU>;

    EIGEN_GENERIC_PUBLIC_INTERFACE(PartialPivLU)
    enum {
      MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };
    typedef PermutationMatrix<RowsAtCompileTime, MaxRowsAtCompileTime> PermutationType;
    typedef Transpositions<RowsAtCompileTime, MaxRowsAtCompileTime> TranspositionType;
    typedef typename MatrixType::PlainObject PlainObject;

    /**
      * \brief Default Constructor.
      *
      * The default constructor is useful in cases in which the user intends to
      * perform decompositions via PartialPivLU::compute(const MatrixType&).
      */
    PartialPivLU();

    template<typename InputType>
    explicit PartialPivLU(EigenBase<InputType>& matrix);
    
    typedef Map<Matrix<Scalar, Dynamic, Dynamic, 0> > MapLU;

    
    /** \internal performs the LU decomposition in-place of the matrix \a lu
      * using an unblocked algorithm.
      *
      * In addition, this function returns the row transpositions in the
      * vector \a row_transpositions which must have a size equal to the number
      * of columns of the matrix \a lu, and an integer \a nb_transpositions
      * which returns the actual number of transpositions.
      *
      * \returns The index of the first pivot which is exactly zero if any, or a negative number otherwise.
      */

    typedef Block<MapLU, Dynamic, Dynamic> MatrixType2;
    // XXX modified by Fabbri to suit Chicago problem
    static __attribute__((always_inline)) void unblocked_lu(
        MatrixType2& lu, 
        typename TranspositionType::StorageIndex* row_transpositions, 
        typename TranspositionType::StorageIndex& nb_transpositions)
    {
      typedef internal::scalar_score_coeff_op<Scalar> Scoring;
      typedef typename Scoring::result_type Score;
      static constexpr Index rows = 14;
      static constexpr Index cols = 14;
      nb_transpositions = 0;
      Index first_zero_pivot = -1;
      for(Index k = 0; k < 14; ++k)
      {
        Index rrows = rows-k-1;
        Index rcols = cols-k-1;

  //      Score biggest_in_corner
  //        = lu.col(k).tail(rows-k).unaryExpr(Scoring()).maxCoeff(&row_of_biggest_in_col);
        
        Index row_of_biggest_in_col(k);
        Score biggest_in_corner = std::norm(lu.coeff(k,k));
        
        for (unsigned j=rows-1; j != k; --j) {
          if (std::norm(lu.coeff(j,k)) > biggest_in_corner*1000) {
              row_of_biggest_in_col = j;
              biggest_in_corner = std::norm(lu.coeff(j,k));
              break;
          }
        }

        row_transpositions[k] = typename TranspositionType::StorageIndex(row_of_biggest_in_col);

        if(biggest_in_corner != Score(0)) {
          if(k != row_of_biggest_in_col) {
            lu.row(k).swap(lu.row(row_of_biggest_in_col));
            ++nb_transpositions;
          }

          // FIXME shall we introduce a safe quotient expression in cas 1/lu.coeff(k,k)
          // overflow but not the actual quotient?
          lu.col(k).tail(rrows) /= lu.coeff(k,k);
        } else if(first_zero_pivot==-1)
          // the pivot is exactly zero, we record the index of the first pivot which is exactly 0,
          // and continue the factorization such we still have A = PLU
          first_zero_pivot = k;

        if (k < rows-1)
          lu.bottomRightCorner(rrows,rcols).noalias() -= lu.col(k).tail(rrows) * lu.row(k).tail(rcols);
      }
    }

    template<typename InputType> inline __attribute__((always_inline)) 
    PartialPivLU& compute(const EigenBase<InputType>& matrix) {
      m_lu = matrix.derived();

      typename TranspositionType::StorageIndex nb_transpositions;
      TranspositionType m_rowsTranspositions;
      {
        MapLU lu1(&m_lu.coeffRef(0,0),m_lu.outerStride(),14);
        Block<MapLU, Dynamic, Dynamic> lu(lu1,0,0,14,14);

        // if the matrix is too small, no blocking:
        unblocked_lu(lu, &m_rowsTranspositions.coeffRef(0), nb_transpositions);
        // template<typename Scalar, int StorageOrder, typename PivIndex>
      }

      m_p = m_rowsTranspositions;
      return *this;
    }

    #ifdef EIGEN_PARSED_BY_DOXYGEN
    /** This method returns the solution x to the equation Ax=b, where A is the matrix of which
      * *this is the LU decomposition.
      *
      * \param b the right-hand-side of the equation to solve. Can be a vector or a matrix,
      *          the only requirement in order for the equation to make sense is that
      *          b.rows()==A.rows(), where A is the matrix of which *this is the LU decomposition.
      *
      * \returns the solution.
      *
      * Example: \include PartialPivLU_solve.cpp
      * Output: \verbinclude PartialPivLU_solve.out
      *
      * Since this PartialPivLU class assumes anyway that the matrix A is invertible, the solution
      * theoretically exists and is unique regardless of b.
      *
      * \sa TriangularView::solve(), inverse(), computeInverse()
      */
    template<typename Rhs>
    inline __attribute__((always_inline)) const Solve<PartialPivLU, Rhs>
    solve(const MatrixBase<Rhs>& b) const;
    #endif

    inline Index rows() const { return m_lu.rows(); }
    inline Index cols() const { return m_lu.cols(); }

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    template<typename RhsType, typename DstType>
    EIGEN_DEVICE_FUNC
    __attribute__((always_inline)) void _solve_impl(const RhsType &rhs, DstType &dst) const {
     /* The decomposition PA = LU can be rewritten as A = P^{-1} L U.
      * So we proceed as follows:
      * Step 1: compute c = Pb.
      * Step 2: replace c by the solution x to Lx = c.
      * Step 3: replace c by the solution x to Ux = c.
      */

      // Step 1
      dst = m_p * rhs;

      // Step 2
      m_lu.template triangularView<UnitLower>().solveInPlace(dst);

      // Step 3
      m_lu.template triangularView<Upper>().solveInPlace(dst);
    }
    #endif

  protected:

    MatrixType m_lu;
    PermutationType m_p;
};

template<typename MatrixType>
PartialPivLU<MatrixType>::PartialPivLU()
  : m_lu(),
    m_p()
{
}


} // end namespace Eigen

#endif // EIGEN_PARTIALLU_H
