
#ifndef UTILS_SUPPORTED_MATERIAL_MODEL_ENUMS_HPP_
#define UTILS_SUPPORTED_MATERIAL_MODEL_ENUMS_HPP_


enum class materialModelKind {unknown, unilayer, bilayer, prem, ak135f};

std::string materialModelKindToString(const materialModelKind e){
  switch (e){
  case materialModelKind::unilayer: return "unilayer";
  case materialModelKind::bilayer:  return "bilayer";
  case materialModelKind::prem:	    return "prem";
  default:			    return "unknown";
  }
}

#endif
