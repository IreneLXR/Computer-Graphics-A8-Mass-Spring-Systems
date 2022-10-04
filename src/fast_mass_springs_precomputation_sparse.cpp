#include "fast_mass_springs_precomputation_sparse.h"
#include "signed_incidence_matrix_sparse.h"
#include <vector>

bool fast_mass_springs_precomputation_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::SparseMatrix<double>  & M,
  Eigen::SparseMatrix<double>  & A,
  Eigen::SparseMatrix<double>  & C,
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization)
{
  /////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  std::vector<Eigen::Triplet<double> > ijv;
  const int n = V.rows();
  //for(int i = 0;i<n;i++) ijv.emplace_back(i,i,1);
  Eigen::SparseMatrix<double> Q(n,n);
  //Q.setFromTriplets(ijv.begin(),ijv.end());
  r = Eigen::VectorXd::Zero(E.rows());
  for (int i = 0; i < E.rows(); i ++)
  {
    r[i] = (V.row(E(i, 0)) - V.row(E(i, 1))).norm();
  }
  M.resize(V.rows(), V.rows());
  for(int i = 0; i < V.rows(); i ++)
  {
    ijv.emplace_back(i, i, m[i]);
  }
  M.setFromTriplets(ijv.begin(), ijv.end());
  ijv.clear();
  signed_incidence_matrix_sparse(V.rows(), E, A);
  C.resize(b.rows(), V.rows());
  for (int i = 0; i < b.rows(); i ++)
  {
    ijv.emplace_back(i, b[i], 1);
  }
  C.setFromTriplets(ijv.begin(), ijv.end());
  ijv.clear();
  Q = k * A.transpose() * A + 1.0 / (pow(delta_t, 2)) * M + C.transpose() * C * 1e10;
  /////////////////////////////////////////////////////////////////////////////
  prefactorization.compute(Q);
  return prefactorization.info() != Eigen::NumericalIssue;
}
