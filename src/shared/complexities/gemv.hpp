
#ifndef UTILS_GEMV_HPP_
#define UTILS_GEMV_HPP_

template <typename sc_t>
struct Gemv
{
  // y = beta*y + alpha*A*x
  // y in R^m
  // A in R^(m,n)
  // x in R^n

  // general case: y = beta*y + alpha*op(A)*x
  static void gemv(const std::size_t m,
		   const std::size_t n,
		   double & memCostMB,
		   double & flops)
  {
    const auto scsz  = sizeof(sc_t);

    const double matread = 1.0*m*n*scsz;
    const double vec_r   = 1.*(m+n)*scsz;
    const double vec_w   = 1.*m*scsz;
    memCostMB = (matread + vec_r + vec_w)/1024./1024.;
    flops = 2.*m*n + 2.*m;
  }

  // y = y + alpha*op(A)*x
  static void gemv_beta_one(const std::size_t m,
			    const std::size_t n,
			    double & memCostMB,
			    double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    const double matread = 1.0*m*n*scsz;
    const double vec_r   = 1.*(m+n)*scsz;
    const double vec_w   = 1.*m*scsz;
    memCostMB = (matread + vec_r + vec_w)/1024./1024.;
    flops = 2.*m*n + 1.*m;
  }

};
#endif
