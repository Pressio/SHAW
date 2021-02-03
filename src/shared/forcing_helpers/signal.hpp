
#ifndef SIGNAL_HPP_
#define SIGNAL_HPP_

template <typename sc_t>
class Signal
{
  static constexpr auto one  = constants<sc_t>::one();
  static constexpr auto two  = constants<sc_t>::two();
  static constexpr auto myPISq = M_PI*M_PI;

  signalKind  myKind_ = {};
  sc_t delayTime_     = {};
  sc_t period_        = {};
  sc_t frequency_     = {};
  sc_t frequencySq_   = frequency_*frequency_;

public:
  Signal() = default;

  Signal(const signalKind kind, const sc_t delayTimeIn, const sc_t periodIn)
    : myKind_{kind},
      delayTime_{delayTimeIn},
      period_{periodIn},
      frequency_{constants<sc_t>::one()/period_}
  {}

  auto getKind() const{ return myKind_; }
  const sc_t getFrequency() const{ return frequency_; }
  const sc_t getPeriod() const{ return period_; }
  const sc_t getDelay() const{ return delayTime_; }

  void resetPeriod(sc_t period){
    period_ = period;
    frequency_ = one/period_;
    frequencySq_ = frequency_*frequency_;
  }

  void operator()(const sc_t & t, sc_t & result) const
  {
    switch (myKind_)
    {
      case signalKind::ricker:
	{
	  const auto tDiffSq = (t-delayTime_)*(t-delayTime_);
	  const auto expTerm = std::exp( -myPISq*frequencySq_*tDiffSq );
	  result = (one - two*myPISq*tDiffSq*frequencySq_) * expTerm;
	  break;
	}

	// sinusoidal (as in original shaxi fortran)
      case signalKind::sinusoid:
	{
	  constexpr sc_t rn    = 0.001;
	  constexpr auto term2 = ( rn/(rn+two) );

	  if (t < period_){
	    const sc_t term1 = std::sin( rn*M_PI*t/period_ );
	    const sc_t term3 = std::sin( (rn+two)*M_PI*t/period_ );
	    result = term1 - term2 * term3;
	  }
	  else{
	    constexpr auto term1 = std::sin( rn*M_PI );
	    constexpr auto term3 = std::sin( (rn+two)*M_PI );
	    result = term1 - term2 * term3;
	  }
	  break;
	}

	// first derivative of Gaussian
      case signalKind::gaussDer:{
	const auto tDiff   = (t-delayTime_);
	const auto tDiffSq = tDiff*tDiff;
	const auto expTerm = std::exp( -frequencySq_*tDiffSq );
	result = -two*tDiff*frequencySq_*expTerm;
	break;
      }

    }
  }//end ()

};
#endif
