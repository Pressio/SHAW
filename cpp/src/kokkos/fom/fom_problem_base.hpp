
#ifndef FOM_PROBLEM_BASE_HPP_
#define FOM_PROBLEM_BASE_HPP_

namespace kokkosapp{

class FomProblemBase{
public:
  virtual void execute() = 0;
};

}//end namespace
#endif
