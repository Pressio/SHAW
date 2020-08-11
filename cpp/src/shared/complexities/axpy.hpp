
#ifndef UTILS_AXPY_HPP_
#define UTILS_AXPY_HPP_

template <typename sc_t>
struct Axpy
{
  // y = y + alpha*x

  // y,x rank-2 R^(m,n) and alpha = rank-1 R^n
  static void axpy(const std::size_t m,
		   const std::size_t n,
		   double & memCostMB,
		   double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    const double reads = 1.0*(2.*m*n + 1.*n)*scsz;
    const double write = 1.0*(m*n)*scsz;
    memCostMB = (reads+write)/1024./1024.;
    flops = 2.*n;
  }

  // alpha = scalar
  static void axpy(const std::size_t n,
		   double & memCostMB,
		   double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    memCostMB = (3.0*n*scsz)/1024./1024.;
    flops = 2.*n;
  }

  // alpha = 0
  static void axpy_alpha_zero(const std::size_t n,
			      double & memCostMB,
			      double & flops)
  {
    memCostMB = 0.;
    flops = 0.;
  }

  // alpha=1
  static void axpy_alpha_one(const std::size_t n,
			     double & memCostMB,
			     double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    memCostMB = (3.0*n*scsz)/1024./1024.;
    flops = 1.*n;
  }

};
#endif
