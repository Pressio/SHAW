
#ifndef SHAXIPP_PARSER_MIXIN_FORCING_SECTION_HPP_
#define SHAXIPP_PARSER_MIXIN_FORCING_SECTION_HPP_

// forward declaration
template <typename sc_t>
class Signal;

template <typename scalar_t>
class ParserForcingSection
{

  bool multiForcing_	= false;
  bool enableRank2Mode_	= false;
  int forcingSize_ = 1;

  signalKind kind_ = {};
  std::vector<scalar_t> depths_  = {}; // km;
  std::vector<scalar_t> angles_  = {}; // degrees;
  std::vector<scalar_t> periods_ = {}; // seconds;
  std::vector<scalar_t> delays_  = {}; // seconds;

public:
  bool rank2Enabled() const {
    return enableRank2Mode_;
  }

  bool multiForcing() const {
    return multiForcing_;
  }

  auto getSourceSignalKind() const{
    return kind_;
  }

  const std::vector<scalar_t> & viewDepths() const { return depths_; }
  const std::vector<scalar_t> & viewAngles() const { return angles_; }
  const std::vector<scalar_t> & viewPeriods() const { return periods_; }
  const std::vector<scalar_t> & viewDelays() const { return delays_; }

public:
  void parseForcing(const std::string & inputFile)
  {
    YAML::Node node = YAML::LoadFile(inputFile);

    // check if section is present, otherwise throw
    const auto node1 = node["source"];
    if (node1)
    {
      const auto node2 = node1["signal"];
      if (node2)
      {
	// check if subnode for "kind" exists
	auto n = node2["kind"];
      	if (n){
	  kind_ = stringToSignalKind(n.as<std::string>());
	  std::cout << "1\n";
	}
      	else{
	  throw std::runtime_error("You must set the kind of the signal");
	}

	// check if subnode for "depth" exists
	auto n1 = node2["depth"];
       	if (n1){
	  // check if it is a scalar or an array
	  if (n1.IsScalar()){
	    depths_.push_back(n1.as<scalar_t>());
	  }
	  else{
	    depths_= n1.as<std::vector<scalar_t>>();
	  }
	  std::cout << "2\n";
	}
      	else{
	  throw std::runtime_error("You must set the depth of the signal");
	}

	auto n2 = node2["angle"];
       	if (n2){
	  // check if it is a scalar or an array
	  if (n2.IsScalar()){
	    angles_.push_back(n2.as<scalar_t>());
	  }
	  else{
	    angles_ = n2.as<std::vector<scalar_t>>();
	  }
	  std::cout << "3\n";
	}
      	else{
	  std::cout << "Angle of source not specified, default==0";
	}

	auto n3 = node2["period"];
       	if (n3){
	  // check if it is a scalar or an array
	  if (n3.IsScalar()){
	    periods_.push_back(n3.as<scalar_t>());
	  }
	  else{
	    periods_ = n3.as<std::vector<scalar_t>>();
	  }
	  std::cout << "4\n";
	}
      	else{
	  throw std::runtime_error("You must set the period of the signal");
	}

	auto n4 = node2["delay"];
       	if (n4){
	  // check if it is a scalar or an array
	  if (n4.IsScalar()){
	    delays_.push_back(n4.as<scalar_t>());
	  }
	  else{
	    delays_ = n4.as<std::vector<scalar_t>>();
	  }
	  std::cout << "5\n";
	}
      	else{
	  throw std::runtime_error("You must set the delay of the signal");
	}

	auto n5 = node2["forcingSize"];
       	if (n5){
	  // if forcingSize is present and >=2, then rank-2 is on
	  forcingSize_ = n5.as<int>();
	  enableRank2Mode_ = true;
	}

	// signalKind sKind = {};
	// scalar_t sDepth = {};
	// scalar_t sAngle = static_cast<scalar_t>(0);
	// scalar_t sPeriod = {};
	// scalar_t sDelay = {};

	// // check if subnode for kind exists
	// auto n = node2["kind"];
      	// if (n){
	//   // note exits, check if it is a vector or just one value
	//   sKind = stringToSignalKind(n.as<std::string>());
	// }
      	// else{
	//   throw std::runtime_error("You must set the kind of the signal");
	// }

	// n = node2["depth"];
       	// if (n){
	//   // const auto vals = n.as<std::vector<scalar_t>>();
	//   // if (vals.size() = 1) std::cout << "SINGLE\n";
	//   sDepth = n.as<scalar_t>();
	// }
      	// else{
	//   throw std::runtime_error("You must set the depth of the signal");
	// }

       	// if (node2["angle"]) sAngle = node2["angle"].as<scalar_t>();
      	// else{
	//   std::cout << "Angle of source not specified, default==0";
	// }

      	// n = node2["period"];
       	// if (n) sPeriod = n.as<scalar_t>();
      	// else throw std::runtime_error("You must set the period of the signal");

      	// n = node2["delay"];
       	// if (n) sDelay = n.as<scalar_t>();
      	// else throw std::runtime_error("You must set the delay of the signal");

	// Signal<scalar_t> signal(sKind, sDelay, sPeriod);
	// source_ = SourceInfo<scalar_t>(sDepth, sAngle, signal);
      }
      else{
      	throw std::runtime_error
	  ("Signal node in inputfile is missing, it is mandatory!");
      }
    }
    else{
      throw std::runtime_error
	("Cannot find source section in inputfile, must be set!");
    }

    setMultiForcingFlag();
    if (multiForcing_){
      validateMultiForcing();
      printMultiForcing();
    }
    else{
      validateSingleForcing();
      printSingleForcing();
    }
  }

private:
  // set multiForcing to true if input file demands so
  void setMultiForcingFlag()
  {
    // we have a multifFOrcing run if any of the source fields is an array
    if (depths_.size() >= 2 or
	angles_.size() >= 2 or
	periods_.size() >= 2 or
	delays_.size() >= 2)
      {
	multiForcing_ = true;
      }
  }

  void validateMultiForcing() const
  {
    if (kind_ == signalKind::unknown){
      const auto ssKstr = signalKindToString(kind_);
      const auto msg = "Invalid source kind: " + ssKstr + ", maybe you mispelled it in input file?\n";
      throw std::runtime_error(msg);
    }

    for (auto & it : depths_){
      if (it<=0.){
	throw std::runtime_error("The depth of the source must be a positive number");
      }
    }

    for (auto & it : angles_){
      if (it<0.){
	throw std::runtime_error("The angle of source must be zero or a positive value");
      }
    }

    for (auto & it : periods_){
      if (it<=0.){
	throw std::runtime_error("Period for signal must be positive. ");
      }
    }

    for (auto & it : delays_){
      if (it<0.){
	throw std::runtime_error("Delay time of the signal must be >=0");
      }
    }
  }

  void validateSingleForcing() const
  {
    if (kind_ == signalKind::unknown){
      const auto ssKstr = signalKindToString(kind_);
      const auto msg = "Invalid source kind: " + ssKstr + ", maybe you mispelled it in input file?\n";
      throw std::runtime_error(msg);
    }

    if (depths_[0]<=0.){
      throw std::runtime_error("The depth of the source must be a positive number");
    }

    if (angles_[0]<0.){
      throw std::runtime_error("The angle of source must be zero or a positive value");
    }

    if (periods_[0]<=0.){
      throw std::runtime_error("Period for signal must be positive. ");
    }

    if (delays_[0]<0.){
      throw std::runtime_error("Delay time of the signal must be >=0");
    }
  }

  void printMultiForcing() const
  {}


  void printSingleForcing() const
  {
    std::cout << std::endl;
    std::cout << "signal "	 << " "
	      << "type= "	 << signalKindToString(kind_) << " "
	      << "depth[km]= "   << depths_[0]  << " "
	      << "angle[deg]= "  << angles_[0]  << " "
	      << "period[sec]= " << periods_[0] << " "
	      << "delay[sec]= "  << delays_[0] << " \n";
  }
};

#endif
