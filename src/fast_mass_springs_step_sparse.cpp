#include "fast_mass_springs_step_sparse.h"
#include <igl/matlab_format.h>

void fast_mass_springs_step_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXi & b,
  const double delta_t,
  const Eigen::MatrixXd & fext,
  const Eigen::VectorXd & r,
  const Eigen::SparseMatrix<double>  & M,
  const Eigen::SparseMatrix<double>  & A,
  const Eigen::SparseMatrix<double>  & C,
  const Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization,
  const Eigen::MatrixXd & Uprev,
  const Eigen::MatrixXd & Ucur,
  Eigen::MatrixXd & Unext)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  /*
  for(int iter = 0;iter < 50;iter++)
  {
    const Eigen::MatrixXd l = Ucur;
    Unext = prefactorization.solve(l);
  }
  */
  Unext = Ucur;
  Eigen::MatrixXd d = Eigen::MatrixXd::Zero(E.rows(), 3);
  Eigen::MatrixXd y = 1.0 / (pow(delta_t, 2)) * M * (2 * Ucur - Uprev) + fext + 1e10 * C.transpose() * C * V;
  for(int iter = 0;iter < 50;iter++)
  {
    // step 1 - local - determine dij
    // step 2 - global - find p
    for (int i = 0; i < E.rows(); i ++)
    {
      d.row(i) = (Unext.row(E(i, 0)) - Unext.row(E(i, 1))).normalized() * r[i];
    }
    Eigen::MatrixXd b_prime = k * A.transpose() * d + y;
    Unext = prefactorization.solve(b_prime);
  }
  //////////////////////////////////////////////////////////////////////////////
}
