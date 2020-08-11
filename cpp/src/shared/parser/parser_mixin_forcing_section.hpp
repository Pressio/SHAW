
#ifndef SHAXIPP_PARSER_MIXIN_FORCING_SECTION_HPP_
#define SHAXIPP_PARSER_MIXIN_FORCING_SECTION_HPP_

template <typename sc_t>
class Signal;

template <typename scalar_t, typename int_t>
struct ParserForcingSectionMixin
{
private:

  template <typename sc_t>
  struct SourceInfo{
    sc_t depth_	   = {}; // km;
    sc_t angle_	   = {}; // degrees;
    Signal<sc_t> signal_;

    SourceInfo() = default;
    SourceInfo(sc_t dep, sc_t ang, const Signal<sc_t> & signalIn)
      : depth_(dep), angle_(ang), signal_(signalIn){}
  };

private:
  SourceInfo<scalar_t> source_;

public:

  ParserForcingSectionMixin() = default;

  auto getSourceSignalKind() const{
    return source_.signal_.getKind();
  }

  auto getSignal() const{
    return source_.signal_;
  }

  scalar_t getSourceProperty(std::string propName) const
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
  void parseForcing(const std::string inputFile)
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
	scalar_t sAngle = {};
	scalar_t sPeriod = {};
	scalar_t sDelay = {};

	auto n = node2["kind"];
      	if (n) sKind = stringToSignalKind(n.as<std::string>());
      	else throw std::runtime_error("You must set kind of the signal");

	n = node2["depth"];
       	if (n) sDepth = n.as<scalar_t>();
      	else throw std::runtime_error("You must set depth of the signal");

	n = node2["angle"];
       	if (n) sAngle = n.as<scalar_t>();
      	else throw std::runtime_error("You must set angle of the signal");

      	n = node2["period"];
       	if (n) sPeriod = n.as<scalar_t>();
      	else throw std::runtime_error("You must set period of the signal");

      	n = node2["delay"];
       	if (n) sDelay = n.as<scalar_t>();
      	else throw std::runtime_error("You must set delay of the signal");

	Signal<scalar_t> signal(sKind, sDelay, sPeriod);
	source_ = SourceInfo<scalar_t>(sDepth, sAngle, signal);
      }
      else
      	throw std::runtime_error("Cannot find signal entry in input");

    }
    else
      throw std::runtime_error("Cannot find source section in yaml");

    this->validate();
    this->print();
  }

private:
  void validate() const
  {
    const auto & thisSrc = source_;
    if (thisSrc.signal_.getKind() == signalKind::unknown){
      std::cerr << "Invalid source type\n";
    }

    if (thisSrc.depth_<=0.){
      std::cerr << "Cannot have negative or zero depth for source" << std::endl;
    }

    if (thisSrc.angle_<0.){
      std::cerr << "Cannot have negative angle for source" << std::endl;
    }

    if (thisSrc.signal_.getPeriod()<=0.){
      std::cerr << "Period for signal should be positive. " << std::endl;
    }

    if (thisSrc.signal_.getDelay()<0.){
      std::cerr << "Delay time of the signal should be >=0" << std::endl;
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
