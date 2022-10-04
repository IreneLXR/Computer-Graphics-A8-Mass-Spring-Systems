#include "signed_incidence_matrix_dense.h"

void signed_incidence_matrix_dense(
  const int n,
  const Eigen::MatrixXi & E,
  Eigen::MatrixXd & A)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  A = Eigen::MatrixXd::Zero(E.rows(),n);
  // for each pair(i, j): A(i, j) += v
  for (int i = 0; i < E.rows(); i ++)
  {
    for (int j = 0; j < n; j ++)
    {
      if (E(i, 0) == j)
      {
        A(i, j) += 1;
      }
      else if (E(i, 1) == j)
      {
        A(i, j) -= 1;
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////////
}
