
#ifndef UTILS_AXPBY_HPP_
#define UTILS_AXPBY_HPP_

#include "axpy.hpp"

template <typename sc_t>
struct Axpby
{
  // y = beta*y + alpha*x

  // y = x
  static void axpby_alpha_one_beta_zero(const std::size_t m,
					double & memCostMB,
					double & flops)
  {
    Axpy<sc_t>::axpy_alpha_one(m, memCostMB, flops);
  }

  // y = y
  static void axpby_alpha_zero_beta_one(const std::size_t m,
					double & memCostMB,
					double & flops)
  {
    Axpy<sc_t>::axpy_alpha_one(m, memCostMB, flops);
  }

  // y =  y + x
  static void axpby_alpha_one_beta_one(const std::size_t m,
				       double & memCostMB,
				       double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    const double reads = 1.0*(2.*m)*scsz;
    const double write = 1.0*m*scsz;
    memCostMB = (reads+write)/1024./1024.;
    flops = 1.*m;
  }

  // y = y + alpha*x
  static void axpby_beta_one(const std::size_t m,
			     double & memCostMB,
			     double & flops)
  {
    Axpy<sc_t>::axpy(m, memCostMB, flops);
  }

  // y = beta*y + x
  static void axpby_alpha_one(const std::size_t m,
			      double & memCostMB,
			      double & flops)
  {
    Axpy<sc_t>::axpy(m, memCostMB, flops);
  }

  // beta,alpha are scalars, x,y are rank-1
  static void axpby(const std::size_t m,
		    double & memCostMB,
		    double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    const double reads = 1.0*(2.*m)*scsz;
    const double write = 1.0*m*scsz;
    memCostMB = (reads+write)/1024./1024.;
    flops = 3.*m;
  }

  // y = beta*y + alpha*x where beta,alpha are rank-1
  // y.extent = m,n
  // x.extent = m,n
  static void axpby(const std::size_t m,
		    const std::size_t n,
		    double & memCostMB,
		    double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    const double reads = 1.0*(2.*m*n + 2.*n)*scsz;
    const double write = 1.0*(m*n)*scsz;
    memCostMB = (reads+write)/1024./1024.;
    flops = 3.*m*n;
  }

};
#endif
