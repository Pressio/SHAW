
#ifndef UTILS_COMPLEXITY_HPP_
#define UTILS_COMPLEXITY_HPP_

#include "./complexities/scal.hpp"
#include "./complexities/spmv.hpp"
#include "./complexities/gemm.hpp"
#include "./complexities/gemv.hpp"
#include "./complexities/axpy.hpp"
#include "./complexities/axpby.hpp"
#include "./complexities/mult.hpp"
#include "./complexities/spmm.hpp"

template <typename sc_t>
struct Complexity
  : Scal<sc_t>,  Spmv<sc_t>,
    Axpby<sc_t>, Axpy<sc_t>,
    Mult<sc_t>,  Spmm<sc_t>,
    Gemv<sc_t>,  Gemm<sc_t>
{};

#endif
