
#ifndef UTILS_SPMM_HPP_
#define UTILS_SPMM_HPP_

template <typename sc_t>
struct Spmm
{
  // C = beta*C + alpha*A*X
  // C in R^(n,k)
  // A sparse size = n,m
  // X in R^(m,k)
  template <typename ordinal_t>
  static void spmm(const std::size_t nnz,
		   const std::size_t nrows,
		   const std::size_t k,
		   double & memCostMB,
		   double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    const auto ordsz = sizeof(ordinal_t);

    const double matsize = 1.0*(nnz*(scsz+ordsz) + nrows*ordsz);
    const double x_r = 1.*nnz*k*scsz;
    const double c_r = 1.*nrows*k*scsz;
    const double c_w = 1.*nrows*k*scsz;

    memCostMB = (matsize + x_r + c_r + c_w)/1024./1024.;
    flops = 2.*nnz*k /*Ax*/ + nrows*k /*c+..*/;
  }

};
#endif
