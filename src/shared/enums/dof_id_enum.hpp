
#ifndef UTILS_DOF_ID_ENUM_HPP_
#define UTILS_DOF_ID_ENUM_HPP_

enum class dofId {unknown, vp, sp};

std::string dofIdToString(const dofId e){
  switch (e){
  default:	     return "unknown";
  case dofId::vp:    return "vp"; // velocity dof
  case dofId::sp:    return "sp"; // stress dof
  }
}

#endif
