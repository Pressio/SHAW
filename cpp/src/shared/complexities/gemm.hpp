
#ifndef UTILS_GEMM_HPP_
#define UTILS_GEMM_HPP_

template <typename sc_t>
struct Gemm
{
  // C = beta*C + alpha*A*X
  // C in R^(m,n)
  // A in R^(m,k)
  // X in R^(k,n)

  // C = beta*C + alpha*A*X
  static void gemm(const std::size_t m,
		   const std::size_t k,
		   const std::size_t n,
		   double & memCostMB,
		   double & flops)
  {
    const auto scsz  = sizeof(sc_t);

    const double matread = 1.0*(m*n + m*k + k*n)*scsz;
    const double mat_w   = 1.*m*n*scsz;
    memCostMB = (matread + mat_w)/1024./1024.;
    flops = 2.*m*k*n + 2.*m*n;
  }

  // C = C + alpha*A*X
  static void gemm_beta_one(const std::size_t m,
			    const std::size_t k,
			    const std::size_t n,
			    double & memCostMB,
			    double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    const double matread = 1.0*(m*n + m*k + k*n)*scsz;
    const double mat_w   = 1.*m*n*scsz;
    memCostMB = (matread + mat_w)/1024./1024.;
    flops = 2.*m*k*n + 1.*m*n;
  }

};
#endif
