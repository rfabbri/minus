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
    typedef PermutationMatrix<14, 14> PermutationType;
    typedef Transpositions<14, 14> TranspositionType;

    /**
      * \brief Default Constructor.
      *
      * The default constructor is useful in cases in which the user intends to
      * perform decompositions via PartialPivLU::compute(const MatrixType&).
      */
    PartialPivLU() : m(), m_p() { }

    template<typename InputType>
    explicit PartialPivLU(EigenBase<InputType>& matrix);

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

    // XXX modified by Fabbri to suit Chicago problem
    static __attribute__((always_inline)) void unblocked_lu(
        MatrixType &lu, 
        typename TranspositionType::StorageIndex* row_transpositions, 
        typename TranspositionType::StorageIndex& nb_transpositions)
    {
      typedef internal::scalar_score_coeff_op<Scalar> Scoring;
      typedef typename Scoring::result_type Score;
      static constexpr Index rows = 14;
      static constexpr Index cols = 14;
      nb_transpositions = 0;
      Index first_zero_pivot = -1;
      for(Index k = 0; k < 14; ++k) {
        Index rrows = rows-k-1;
        Index rcols = cols-k-1;

  //      Score biggest_in_corner
  //        = lu.col(k).tail(rows-k).unaryExpr(Scoring()).maxCoeff(&row_of_biggest_in_col);
        
        Index row_of_biggest_in_col(k);
        Score biggest_in_corner = std::norm(lu(k,k));// std::norm(lu.coeff(k,k));
        for (unsigned j=rows-1; j != k; --j) {
          Score tmp;
          if ((tmp = std::norm(lu(j,k))) > biggest_in_corner*1000) {
              biggest_in_corner = tmp;
              row_of_biggest_in_col = j;
              break;
          }
        }

        row_transpositions[k] = typename TranspositionType::StorageIndex(row_of_biggest_in_col);

        if (biggest_in_corner != Score(0)) {
          if (k != row_of_biggest_in_col) {
            lu.row(k).swap(lu.row(row_of_biggest_in_col));
            ++nb_transpositions;
          }

          lu.col(k).tail(rrows) /= lu(k,k);
        } else if (first_zero_pivot==-1)
          // the pivot is exactly zero, we record the index of the first pivot which is exactly 0,
          // and continue the factorization such we still have A = PLU
          first_zero_pivot = k;

        if (k < rows-1)
          lu.bottomRightCorner(rrows,rcols).noalias() -= lu.col(k).tail(rrows) * lu.row(k).tail(rcols);
      }
    }

    template<typename InputType> inline __attribute__((always_inline)) 
    PartialPivLU& compute(const EigenBase<InputType>& matrix) {
      m = matrix.derived();
      typename TranspositionType::StorageIndex nb_transpositions;
      TranspositionType m_rowsTranspositions;
      unblocked_lu(m, &m_rowsTranspositions.coeffRef(0), nb_transpositions);

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

    inline Index rows() const { return m.rows(); }
    inline Index cols() const { return m.cols(); }

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    template<typename RhsType, typename DstType>
    EIGEN_DEVICE_FUNC
    __attribute__((always_inline)) void _solve_impl(const RhsType &rhs, DstType &d) const {
     /* The decomposition PA = LU can be rewritten as A = P^{-1} L U.
      * So we proceed as follows:
      * Step 1: compute c = Pb.
      * Step 2: replace c by the solution x to Lx = c.
      * Step 3: replace c by the solution x to Ux = c.
      */

      // Step 1
      d = m_p * rhs;

      // TODO: use block indexing and std::vector-std::vector multiplication
      d(1)  -= m(1,0)*d(0);
      d(2)  -= m(2,0)*d(0)+ m(2,1)*d(1);
      d(3)  -= m(3,0)*d(0)+ m(3,1)*d(1)+ m(3,2)*d(2);
      d(4)  -= m(4,0)*d(0)+ m(4,1)*d(1)+ m(4,2)*d(2)+ m(4,3)*d(3);
      d(5)  -= m(5,0)*d(0)+ m(5,1)*d(1)+ m(5,2)*d(2)+ m(5,3)*d(3)+ m(5,4)*d(4);
      d(6)  -= m(6,0)*d(0)+ m(6,1)*d(1)+ m(6,2)*d(2)+ m(6,3)*d(3)+ m(6,4)*d(4)+ m(6,5)*d(5);
      d(7)  -= m(7,0)*d(0)+ m(7,1)*d(1)+ m(7,2)*d(2)+ m(7,3)*d(3)+ m(7,4)*d(4)+ m(7,5)*d(5)+ m(7,6)*d(6);
      d(8)  -= m(8,0)*d(0)+ m(8,1)*d(1)+ m(8,2)*d(2)+ m(8,3)*d(3)+ m(8,4)*d(4)+ m(8,5)*d(5)+ m(8,6)*d(6)+ m(8,7)*d(7);
      d(9)  -= m(9,0)*d(0)+ m(9,1)*d(1)+ m(9,2)*d(2)+ m(9,3)*d(3)+ m(9,4)*d(4)+ m(9,5)*d(5)+ m(9,6)*d(6)+ m(9,7)*d(7)+ m(9,8)*d(8);
      d(10) -= m(10,0)*d(0)+ m(10,1)*d(1)+ m(10,2)*d(2)+ m(10,3)*d(3)+ m(10,4)*d(4)+ m(10,5)*d(5)+ m(10,6)*d(6)+ m(10,7)*d(7)+ m(10,8)*d(8)+ m(10,9)*d(9);
      d(11) -= m(11,0)*d(0)+ m(11,1)*d(1)+ m(11,2)*d(2)+ m(11,3)*d(3)+ m(11,4)*d(4)+ m(11,5)*d(5)+ m(11,6)*d(6)+ m(11,7)*d(7)+ m(11,8)*d(8)+ m(11,9)*d(9)+ m(11,10)*d(10);
      d(12) -= m(12,0)*d(0)+ m(12,1)*d(1)+ m(12,2)*d(2)+ m(12,3)*d(3)+ m(12,4)*d(4)+ m(12,5)*d(5)+ m(12,6)*d(6)+ m(12,7)*d(7)+ m(12,8)*d(8)+ m(12,9)*d(9)+ m(12,10)*d(10)+ m(12,11)*d(11);
      d(13) -= m(13,0)*d(0)+ m(13,1)*d(1)+ m(13,2)*d(2)+ m(13,3)*d(3)+ m(13,4)*d(4)+ m(13,5)*d(5)+ m(13,6)*d(6)+ m(13,7)*d(7)+ m(13,8)*d(8)+ m(13,9)*d(9)+ m(13,10)*d(10)+ m(13,11)*d(11)+ m(13,12)*d(12);

      
      // Step 2
      //m.template triangularView<UnitLower>().solveInPlace(d);

      d(13) /= m(13,13);
      d(12) -= m(12,13)*d(13);
      d(12) /= m(12,12);
      d(11) -= (m(11,12)*d(12)+ m(11,13)*d(13));
      d(11) /= m(11,11);
      d(10) -= (m(10,11)*d(11)+ m(10,12)*d(12)+ m(10,13)*d(13));
      d(10) /= m(10,10);
      d(9)  -= (m(9,10)*d(10)+ m(9,11)*d(11)+ m(9,12)*d(12)+ m(9,13)*d(13));
      d(9) /= m(9,9);
      d(8)  -= (m(8,9)*d(9)+ m(8,10)*d(10)+ m(8,11)*d(11)+ m(8,12)*d(12)+ m(8,13)*d(13));
      d(8) /= m(8,8);
      d(7)  -= (m(7,8)*d(8)+ m(7,9)*d(9)+ m(7,10)*d(10)+ m(7,11)*d(11)+ m(7,12)*d(12)+ m(7,13)*d(13));
      d(7) /= m(7,7);
      d(6)  -= (m(6,7)*d(7)+ m(6,8)*d(8)+ m(6,9)*d(9)+ m(6,10)*d(10)+ m(6,11)*d(11)+ m(6,12)*d(12)+ m(6,13)*d(13));
      d(6) /= m(6,6);
      d(5)  -= (m(5,6)*d(6)+ m(5,7)*d(7)+ m(5,8)*d(8)+ m(5,9)*d(9)+ m(5,10)*d(10)+ m(5,11)*d(11)+ m(5,12)*d(12)+ m(5,13)*d(13));
      d(5) /= m(5,5);
      d(4)  -= (m(4,5)*d(5)+ m(4,6)*d(6)+ m(4,7)*d(7)+ m(4,8)*d(8)+ m(4,9)*d(9)+ m(4,10)*d(10)+ m(4,11)*d(11)+ m(4,12)*d(12)+ m(4,13)*d(13));
      d(4) /= m(4,4);
      d(3)  -= (m(3,4)*d(4)+ m(3,5)*d(5)+ m(3,6)*d(6)+ m(3,7)*d(7)+ m(3,8)*d(8)+ m(3,9)*d(9)+ m(3,10)*d(10)+ m(3,11)*d(11)+ m(3,12)*d(12)+ m(3,13)*d(13));
      d(3) /= m(3,3);
      d(2)  -= (m(2,3)*d(3)+ m(2,4)*d(4)+ m(2,5)*d(5)+ m(2,6)*d(6)+ m(2,7)*d(7)+ m(2,8)*d(8)+ m(2,9)*d(9)+ m(2,10)*d(10)+ m(2,11)*d(11)+ m(2,12)*d(12)+ m(2,13)*d(13));
      d(2) /= m(2,2);
      d(1)  -= (m(1,2)*d(2)+ m(1,3)*d(3)+ m(1,4)*d(4)+ m(1,5)*d(5)+ m(1,6)*d(6)+ m(1,7)*d(7)+ m(1,8)*d(8)+ m(1,9)*d(9)+ m(1,10)*d(10)+ m(1,11)*d(11)+ m(1,12)*d(12)+ m(1,13)*d(13));
      d(1) /= m(1,1);
      d(0)  -= (m(0,1)*d(1)+ m(0,2)*d(2)+ m(0,3)*d(3)+ m(0,4)*d(4)+ m(0,5)*d(5)+ m(0,6)*d(6)+ m(0,7)*d(7)+ m(0,8)*d(8)+ m(0,9)*d(9)+ m(0,10)*d(10)+ m(0,11)*d(11)+ m(0,12)*d(12)+ m(0,13)*d(13));
      d(0) /= m(0,0);
      
      // d(0)  = d(0) - (m.row(0).tail(12) * d.tail(12).transpose());

      // Step 3
      //m.template triangularView<Upper>().solveInPlace(d);
    }
    #endif

  protected:

    MatrixType m; // matrix holding LU together
    PermutationType m_p;
};


} // end namespace Eigen

#endif // EIGEN_PARTIALLU_H
