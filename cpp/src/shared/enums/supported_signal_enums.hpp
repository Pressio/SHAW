
#ifndef UTILS_SIGNAL_ENUMS_HPP_
#define UTILS_SIGNAL_ENUMS_HPP_

enum class signalKind {unknown, ricker, sinusoid, gaussDer};

std::string signalKindToString(const signalKind e){
  switch (e){
  case signalKind::ricker:    return "ricker";
  case signalKind::sinusoid:  return "sinusoid";
  case signalKind::gaussDer:  return "gaussianDeriv";
  default:		      return "unknown";
  }
}


signalKind stringToSignalKind(const std::string s){
  if (s == "ricker" or s=="Ricker")
    return signalKind::ricker;
  else if (s == "sinusoid" or s=="Sinusoid")
    return signalKind::sinusoid;
  else if (s == "gaussDer" or s=="gaussianDeriv")
    return signalKind::gaussDer;
  else
    return signalKind::unknown;
}

#endif
