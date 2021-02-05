
#ifndef READ_REF_STATE_HPP_
#define READ_REF_STATE_HPP_

#ifdef SHW_HAVE_TPL_EIGEN
template <typename sc_t, typename int_t, typename v_t>
typename std::enable_if< is_vector_eigen<v_t>::value, v_t >::type
readRefState(const std::string fileName, const int_t useBinary)
{
  v_t a;
  throw std::runtime_error("readRefState not implemented yet");
  // if (useBinary == 1){
  //   readBinaryVectorWithSize(fileName, a);
  // }
  return a;
}
#endif

#endif
