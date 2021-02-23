
#ifndef UTILS_MULT_HPP_
#define UTILS_MULT_HPP_

template <typename sc_t>
struct Mult
{
  // y = beta*y + alpha*a*x

  // regular case
  static void mult(const std::size_t n,
		   double & memCostMB,
		   double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    memCostMB = (4.0*n*scsz)/1024./1024.;
    flops = 4.0*n;
  }

  // // y = beta*y
  // static void mult_alpha_zero(const std::size_t n,
  // 			      double & memCostMB,
  // 			      double & flops)
  // {
  //   const auto scsz  = sizeof(sc_t);
  //   memCostMB = (2.0*n*scsz)/1024./1024.;
  //   flops = 1.0*n;
  // }

  // // y = a*x
  // static void mult_alpha_one_beta_zero(const std::size_t n,
  // 				       double & memCostMB,
  // 				       double & flops)
  // {
  //   const auto scsz  = sizeof(sc_t);
  //   memCostMB = (3.0*n*scsz)/1024./1024.;
  //   flops = 1.0*n;
  // }

  // y = y + alpha*a*x
  static void mult_beta_one(const std::size_t n,
			    double & memCostMB,
			    double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    memCostMB = (4.0*n*scsz)/1024./1024.;
    flops = 3.0*n;
  }

  // // y = beta*y + a*x
  // static void mult_alpha_one(const std::size_t n,
  // 			     double & memCostMB,
  // 			     double & flops)
  // {
  //   const auto scsz  = sizeof(sc_t);
  //   memCostMB = (4.0*n*scsz)/1024./1024.;
  //   flops = 3.0*n;
  // }

  // y = y + a*x
  static void mult_alpha_one_beta_one(const std::size_t n,
				      double & memCostMB,
				      double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    memCostMB = (4.0*n*scsz)/1024./1024.;
    flops = 2.0*n;
  }

};
#endif
