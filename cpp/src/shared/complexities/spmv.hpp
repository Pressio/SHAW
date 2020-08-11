
#ifndef UTILS_SPMV_HPP_
#define UTILS_SPMV_HPP_

template <typename sc_t>
struct Spmv
{
  //https://books.google.it/books?id=O22cDwAAQBAJ&pg=PA97&lpg=PA97&dq=gemm+flops+supercomputing+frontiers&source=bl&ots=tiCjQSEypI&sig=ACfU3U2D2zTZHLGPVXNzIFV8SBO_jhz-1w&hl=en&sa=X&ved=2ahUKEwjNieO9nsXpAhVuwqYKHaRJDZIQ6AEwAHoECAoQAQ#v=onepage&q=gemm%20flops%20supercomputing%20frontiers&f=false

  // y = beta*y + alpha*A*x
  template <typename ordinal_t>
  static void spmv(const std::size_t nnz,
		   const std::size_t nrows,
		   double & memCostMB,
		   double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    const auto ordsz = sizeof(ordinal_t);

    const double matsize = 1.0*(nnz*(scsz+ordsz) + nrows*ordsz);
    const double vec_x_r = 1.*nnz*scsz;
    const double vec_y_r = 1.*nrows*scsz;
    const double vec_y_w = 1.*nrows*scsz;

    memCostMB = (matsize + vec_x_r + vec_y_r + vec_y_w)/1024./1024.;
    flops = 2.0*nnz /*Ax*/ + nrows /*y+...*/;
  }

};
#endif
