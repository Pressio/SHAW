
#ifndef UTILS_SCAL_HPP_
#define UTILS_SCAL_HPP_

template <typename sc_t>
struct Scal
{

public:
  // ------------------------------
  // y = alpha*x
  // y,x = rank-1
  // ------------------------------

  // when alpha=0, nothing happens
  static void scal_alpha_zero(const std::size_t n,
			     double & memCostMB,
			     double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    memCostMB = (1.0*n*scsz)/1024./1024.;
    flops = 0.;
  }

  // when alpha=1, y=x
  static void scal_alpha_one(const std::size_t n,
			     double & memCostMB,
			     double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    memCostMB = (2.0*n*scsz)/1024./1024.;
    flops = 0.;
  }

  // regular case
  static void scal(const std::size_t n,
		   double & memCostMB,
		   double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    memCostMB = (2.0*n*scsz)/1024./1024.;
    flops = 1.0*n;
  }
};

#endif
