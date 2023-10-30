// Specific to Chicago
template <problem P, typename F>
__attribute__((always_inline)) inline void
lsolve(
    Map<Matrix<C<F>, minus_core<P,F>::f::nve, minus_core<P,F>::f::nve>,Aligned> & __restrict m, 
    Map<const Matrix<C<F>, minus_core<P,F>::f::nve, 1>, Aligned > & __restrict b,
    Map<Matrix<C<F>, minus_core<P,F>::f::nve, 1>,Aligned> & __restrict x)
{
  typedef minus_core<P, F> M;
  typedef Matrix<C<F>, M::f::nve, M::f::nve>  MatrixType;
  typedef PermutationMatrix<M::f::nve, M::f::nve> PermutationType;
  typedef Transpositions<M::f::nve, M::f::nve> TranspositionType;
  PermutationType m_p; TranspositionType m_rowsTranspositions;
  typename TranspositionType::StorageIndex* row_transpositions = &m_rowsTranspositions.coeffRef(0);
  static constexpr Index rows = M::f::nve;
  for(Index k = 0; k < M::f::nve; ++k) {
    Index rrows = rows-k-1, row_of_biggest_in_col(k);
    F biggest_in_corner = std::norm(m(k,k));
    for (unsigned j=rows-1; j != k; --j) {
      F tmp;
      if ((tmp = std::norm(m(j,k))) > biggest_in_corner*1000) {
          biggest_in_corner = tmp;
          row_of_biggest_in_col = j;
          break;
      }
    }
    row_transpositions[k] = typename TranspositionType::StorageIndex(row_of_biggest_in_col);
    if (k != row_of_biggest_in_col) m.row(k).swap(m.row(row_of_biggest_in_col));

    m.col(k).tail(rrows) /= m(k,k);

    if (k < rows-1)
      m.bottomRightCorner(rrows,rrows).noalias() -= m.col(k).tail(rrows) * m.row(k).tail(rrows);
  }
  m_p = m_rowsTranspositions;
  x = m_p * b;

  // TODO: use block indexing and std::vector-std::vector multiplication
  x(1)  -= m(1,0)*x(0);
  x(2)  -= m(2,0)*x(0)+m(2,1)*x(1);
  x(3)  -= m(3,0)*x(0)+m(3,1)*x(1)+m(3,2)*x(2);
  x(4)  -= m(4,0)*x(0)+m(4,1)*x(1)+m(4,2)*x(2)+m(4,3)*x(3);
  x(5)  -= m(5,0)*x(0)+m(5,1)*x(1)+m(5,2)*x(2)+m(5,3)*x(3)+m(5,4)*x(4);
  x(6)  -= m(6,0)*x(0)+m(6,1)*x(1)+m(6,2)*x(2)+m(6,3)*x(3)+m(6,4)*x(4)+m(6,5)*x(5);
  x(7)  -= m(7,0)*x(0)+m(7,1)*x(1)+m(7,2)*x(2)+m(7,3)*x(3)+m(7,4)*x(4)+m(7,5)*x(5)+m(7,6)*x(6);
  x(8)  -= m(8,0)*x(0)+m(8,1)*x(1)+m(8,2)*x(2)+m(8,3)*x(3)+m(8,4)*x(4)+m(8,5)*x(5)+m(8,6)*x(6)+m(8,7)*x(7);
  x(9)  -= m(9,0)*x(0)+m(9,1)*x(1)+m(9,2)*x(2)+m(9,3)*x(3)+m(9,4)*x(4)+m(9,5)*x(5)+m(9,6)*x(6)+m(9,7)*x(7)+m(9,8)*x(8);
  x(10) -= m(10,0)*x(0)+m(10,1)*x(1)+m(10,2)*x(2)+m(10,3)*x(3)+m(10,4)*x(4)+m(10,5)*x(5)+m(10,6)*x(6)+m(10,7)*x(7)+m(10,8)*x(8)+m(10,9)*x(9);
  x(11) -= m(11,0)*x(0)+m(11,1)*x(1)+m(11,2)*x(2)+m(11,3)*x(3)+m(11,4)*x(4)+m(11,5)*x(5)+m(11,6)*x(6)+m(11,7)*x(7)+m(11,8)*x(8)+m(11,9)*x(9)+m(11,10)*x(10);
  x(12) -= m(12,0)*x(0)+m(12,1)*x(1)+m(12,2)*x(2)+m(12,3)*x(3)+m(12,4)*x(4)+m(12,5)*x(5)+m(12,6)*x(6)+m(12,7)*x(7)+m(12,8)*x(8)+m(12,9)*x(9)+m(12,10)*x(10)+m(12,11)*x(11);
  x(13) -= (m(13,0)*x(0)+m(13,1)*x(1)+m(13,2)*x(2)+m(13,3)*x(3)+m(13,4)*x(4)+m(13,5)*x(5)+m(13,6)*x(6)+m(13,7)*x(7)+m(13,8)*x(8)+m(13,9)*x(9)+m(13,10)*x(10)+m(13,11)*x(11)+m(13,12)*x(12));

  // Step 2
  //m.template triangularView<UnitLower>().solveInPlace(x);

  x(13) /= m(13,13);
  x(12) -= m(12,13)*x(13); x(12) /= m(12,12);
  x(11) -= (m(11,12)*x(12)+m(11,13)*x(13)); x(11) /= m(11,11);
  x(10) -= (m(10,11)*x(11)+m(10,12)*x(12)+m(10,13)*x(13)); x(10) /= m(10,10);
  x(9)  -= (m(9,10)*x(10)+m(9,11)*x(11)+m(9,12)*x(12)+m(9,13)*x(13)); x(9) /= m(9,9);
  x(8)  -= (m(8,9)*x(9)+m(8,10)*x(10)+m(8,11)*x(11)+m(8,12)*x(12)+m(8,13)*x(13)); x(8) /= m(8,8);
  x(7)  -= (m(7,8)*x(8)+m(7,9)*x(9)+m(7,10)*x(10)+m(7,11)*x(11)+m(7,12)*x(12)+m(7,13)*x(13)); x(7) /= m(7,7);
  x(6)  -= (m(6,7)*x(7)+m(6,8)*x(8)+m(6,9)*x(9)+m(6,10)*x(10)+m(6,11)*x(11)+m(6,12)*x(12)+m(6,13)*x(13)); x(6) /= m(6,6);
  x(5)  -= (m(5,6)*x(6)+m(5,7)*x(7)+m(5,8)*x(8)+m(5,9)*x(9)+m(5,10)*x(10)+m(5,11)*x(11)+m(5,12)*x(12)+m(5,13)*x(13)); x(5) /= m(5,5);
  x(4)  -= (m(4,5)*x(5)+m(4,6)*x(6)+m(4,7)*x(7)+m(4,8)*x(8)+m(4,9)*x(9)+m(4,10)*x(10)+m(4,11)*x(11)+m(4,12)*x(12)+m(4,13)*x(13)); x(4) /= m(4,4);
  x(3)  -= (m(3,4)*x(4)+m(3,5)*x(5)+m(3,6)*x(6)+m(3,7)*x(7)+m(3,8)*x(8)+m(3,9)*x(9)+m(3,10)*x(10)+m(3,11)*x(11)+m(3,12)*x(12)+m(3,13)*x(13)); x(3) /= m(3,3);
  x(2)  -= (m(2,3)*x(3)+m(2,4)*x(4)+m(2,5)*x(5)+m(2,6)*x(6)+m(2,7)*x(7)+m(2,8)*x(8)+m(2,9)*x(9)+m(2,10)*x(10)+m(2,11)*x(11)+m(2,12)*x(12)+m(2,13)*x(13)); x(2) /= m(2,2);
  x(1)  -= (m(1,2)*x(2)+m(1,3)*x(3)+m(1,4)*x(4)+m(1,5)*x(5)+m(1,6)*x(6)+m(1,7)*x(7)+m(1,8)*x(8)+m(1,9)*x(9)+m(1,10)*x(10)+m(1,11)*x(11)+m(1,12)*x(12)+m(1,13)*x(13)); x(1) /= m(1,1);
  x(0)  -= (m(0,1)*x(1)+m(0,2)*x(2)+m(0,3)*x(3)+m(0,4)*x(4)+m(0,5)*x(5)+m(0,6)*x(6)+m(0,7)*x(7)+m(0,8)*x(8)+m(0,9)*x(9)+m(0,10)*x(10)+m(0,11)*x(11)+m(0,12)*x(12)+m(0,13)*x(13)); x(0) /= m(0,0);
}
