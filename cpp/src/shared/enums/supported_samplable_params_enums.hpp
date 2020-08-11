
#ifndef UTILS_SAMPLABLE_PARAMS_ENUMS_HPP_
#define UTILS_SAMPLABLE_PARAMS_ENUMS_HPP_

enum class samplable {unknown, signalPeriod};

std::string samplableToString(const samplable e){
  switch (e){
  case samplable::signalPeriod:  return "signalPeriod";
  default:		      return "unknown";
  }
}

samplable stringToSamplable(const std::string s){
  if (s == "signalPeriod")
    return samplable::signalPeriod;
  else
    return samplable::unknown;
}

#endif
