
#ifndef SHAXIPP_PARSER_MIXIN_FORCING_SECTION_HPP_
#define SHAXIPP_PARSER_MIXIN_FORCING_SECTION_HPP_

template <typename sc_t>
class Signal;

template <typename scalar_t, typename int_t>
class ParserForcingSection
{

  template <typename sc_t>
  struct SourceInfo
  {
    sc_t depth_	   = {}; // km;
    sc_t angle_	   = {}; // degrees;
    Signal<sc_t> signal_;

    SourceInfo() = default;
    SourceInfo(sc_t dep, sc_t ang, const Signal<sc_t> & signalIn)
      : depth_(dep), angle_(ang), signal_(signalIn){}
  };

  SourceInfo<scalar_t> source_;

public:
  auto getSourceSignalKind() const{
    return source_.signal_.getKind();
  }

  auto getSignal() const{
    return source_.signal_;
  }

  scalar_t getSourceProperty(const std::string & propName) const
  {
    if (propName == "depth")
      return source_.depth_;
    else if (propName == "angle")
      return source_.angle_;
    else if (propName == "period")
      return source_.signal_.getPeriod();
    else if (propName == "delay")
      return source_.signal_.getDelay();
    else{
      const std::string message = "You queried an invalid property: " + propName;
      throw std::runtime_error(message);
    }
    return 0.;
  }

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

	signalKind sKind = {};
	scalar_t sDepth = {};
	scalar_t sAngle = static_cast<scalar_t>(0);
	scalar_t sPeriod = {};
	scalar_t sDelay = {};

	auto n = node2["kind"];
      	if (n) sKind = stringToSignalKind(n.as<std::string>());
      	else throw std::runtime_error("You must set the kind of the signal");

	n = node2["depth"];
       	if (n) sDepth = n.as<scalar_t>();
      	else throw std::runtime_error("You must set the depth of the signal");

       	if (node2["angle"]) sAngle = node2["angle"].as<scalar_t>();
      	else{
	  std::cout << "Angle of source not specified, default==0";
	}

      	n = node2["period"];
       	if (n) sPeriod = n.as<scalar_t>();
      	else throw std::runtime_error("You must set the period of the signal");

      	n = node2["delay"];
       	if (n) sDelay = n.as<scalar_t>();
      	else throw std::runtime_error("You must set the delay of the signal");

	Signal<scalar_t> signal(sKind, sDelay, sPeriod);
	source_ = SourceInfo<scalar_t>(sDepth, sAngle, signal);
      }
      else{
      	throw std::runtime_error
	  ("Cannot find signal entry within the source section in inputfile, it is mandatory!");
      }
    }
    else{
      throw std::runtime_error
	("Cannot find source section in inputfile, must be set!");
    }

    this->validate();
    this->print();
  }

private:
  void validate() const
  {
    if (source_.signal_.getKind() == signalKind::unknown){
      const auto ssKstr = signalKindToString(source_.signal_.getKind());
      const auto msg = "Invalid source kind: " + ssKstr + ", maybe you mispelled it in input file?\n";
      throw std::runtime_error(msg);
    }

    if (source_.depth_<=0.){
      throw std::runtime_error("The depth of the source must be a positive number");
    }

    if (source_.angle_<0.){
      throw std::runtime_error("The angle of source must be zero or a positive value");
    }

    if (source_.signal_.getPeriod()<=0.){
      throw std::runtime_error("Period for signal must be positive. ");
    }

    if (source_.signal_.getDelay()<0.){
      throw std::runtime_error("Delay time of the signal must be >=0");
    }
}

  void print() const
  {
    const auto & thisSrc = source_;
    std::cout << std::endl;
    std::cout << "signal "	 << " "
	      << "type= "	 << signalKindToString(thisSrc.signal_.getKind()) << " "
	      << "depth[km]= "   << thisSrc.depth_     << " "
	      << "angle[deg]= "  << thisSrc.angle_	 << " "
	      << "period[sec]= " << thisSrc.signal_.getPeriod()    << " "
	      << "delay[sec]= "  << thisSrc.signal_.getDelay() << " \n";
  }
};

#endif
