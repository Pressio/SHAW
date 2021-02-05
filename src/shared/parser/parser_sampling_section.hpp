
#ifndef SHAXIPP_PARSER_MIXIN_SAMPLING_SECTION_HPP_
#define SHAXIPP_PARSER_MIXIN_SAMPLING_SECTION_HPP_

template <typename scalar_t, typename int_t>
struct ParserSamplingSection
{
private:
  bool samplingOn_  = false;
  int_t numSamples_ = {};
  int_t numParams_  = {};

  std::vector<samplable> paramNames_;
  std::vector<std::vector<scalar_t>> values_;

  bool enableRankTwoForcing_ = false;
  int_t forcingSize_ = 1;

public:
  auto enableSampling() const{ return samplingOn_; };
  auto getNumSamples() const{ return numSamples_; }
  auto getNumParams() const{ return numParams_; }

  auto getNameParamToSample(int_t i) const{
    if (i >= numParams_)
      throw std::runtime_error("Trying to index non-existing sampling parameter");

    return paramNames_[i];
  }
  auto getValues(int_t i) const{ return values_[i]; }
  auto getForcingSize() const{ return forcingSize_; };
  auto enableMultiForcing() const{ return enableRankTwoForcing_; };

public:
  void parseSampling(const std::string inputFile)
  {
    const YAML::Node node0 = YAML::LoadFile(inputFile);

    const auto node = node0["sampling"];
    if (node)
    {
      // if node found, then sampling is on
      samplingOn_ = true;

      {
	const auto n = node["params"];
	if (n){
	  const auto pars = n.as<std::vector<std::string>>();
	  paramNames_.resize(pars.size());
	  values_.resize(pars.size());
	  numParams_ = pars.size();
	  if (numParams_ > 1)
	    throw std::runtime_error("Currently supporting sampling for only 1 param= signalPeriod");
	  if (pars[0] != "signalPeriod")
	    throw std::runtime_error("Currently supporting only sampling for signalPeriod");

	  for (auto i=0; i<pars.size(); ++i)
	    paramNames_[i] = stringToSamplable(pars[i]);
	}
	else throw std::runtime_error("You must set names of params to target");
      }

      {
	const auto n = node["values"];
	if (n){
	  const auto vals = n.as<std::vector<scalar_t>>();
	  numSamples_ = vals.size();
	  values_[0] = vals;
	}
	else
	  throw std::runtime_error("You must provide the values for the sampling param as list");
      }

      {
	const auto n = node["forcingSize"];
	if (n){
      	  forcingSize_ = n.as<int_t>();
      	  // enable rank2 forcing only if size > 1
      	  enableRankTwoForcing_ = forcingSize_ > 1;

	  // make sure the number of samples if divisible by the size
	  if (numSamples_ % forcingSize_ != 0)
	    throw std::runtime_error("numSamples not divisible by forcing size, it must be");
	}
      }

      this->validate();
    }
    this->print();
  }

private:
  void validate() const
  {
    // if (vpSamplingSize_<=0) throw std::runtime_error("Cannot have vpSamplingSize <=0 ");
    // if (spSamplingSize_<=0) throw std::runtime_error("Cannot have spSamplingSize <=0 ");
  }

  void print() const
  {
    std::cout << std::endl;
    std::cout << "enableSampling = " << std::boolalpha << samplingOn_ << " \n";
    if (samplingOn_){
      std::cout << "numParams_ = "  << numParams_ << " \n"
		<< "paramName_ = "  << samplableToString(paramNames_[0]) << " \n";

      std::cout << "values = ";
      for (auto & it : values_[0])
	std::cout << it << " ";

      std::cout << "\nenableRankTwoF_ = "  << enableRankTwoForcing_  << " \n"
		<< "forcingSize_ = "  << forcingSize_ << " \n";
    }
  }
};

#endif
