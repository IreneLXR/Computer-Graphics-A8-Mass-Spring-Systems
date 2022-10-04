#include "signed_incidence_matrix_sparse.h"
#include <vector>

void signed_incidence_matrix_sparse(
  const int n,
  const Eigen::MatrixXi & E,
  Eigen::SparseMatrix<double>  & A)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  std::vector<Eigen::Triplet<double> > ijv;
  for (int i = 0; i < E.rows(); i ++)
  {
    for (int j = 0; j < n; j ++)
    {
      if (E(i, 0) == j)
      {
        ijv.emplace_back(i, j, 1.0);
      }
      else if (E(i, 1) == j)
      {
        ijv.emplace_back(i, j, - 1.0);
      }
    }
  }
  A.resize(E.rows(),n);
  A.setFromTriplets(ijv.begin(),ijv.end());
  //////////////////////////////////////////////////////////////////////////////
}
