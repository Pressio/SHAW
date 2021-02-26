/*
//@HEADER
// ************************************************************************
//
// signal.hpp
//                     		Pressio/SHAW
//                         Copyright 2019
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
	  using namespace std;
    	  const auto tDiffSq = (t-delayTime_)*(t-delayTime_);
    	  const auto expTerm = exp( -myPISq*frequencySq_*tDiffSq );
    	  result = (one - two*myPISq*tDiffSq*frequencySq_) * expTerm;
    	  break;
    	}

    	// sinusoidal (as in original shaxi fortran)
      case signalKind::sinusoid:
    	{
    	  constexpr auto rn    = static_cast<sc_t>(0.001);
    	  constexpr auto term2 = rn/(rn+two);

    	  if (t < period_){
    	    const auto term1 = std::sin(rn*M_PI*t/period_ );
    	    const auto term3 = std::sin((rn+two)*M_PI*t/period_ );
    	    result = term1 - term2 * term3;
    	  }
    	  else{
    	    const auto term1 = std::sin( rn*M_PI );
    	    const auto term3 = std::sin((rn+two)*M_PI );
    	    result = term1 - term2 * term3;
    	  }
    	  break;
    	}

    	// first derivative of Gaussian
      case signalKind::gaussDer:
	{
	  const auto tDiff   = (t-delayTime_);
	  const auto tDiffSq = tDiff*tDiff;
	  const auto expTerm = std::exp( -frequencySq_*tDiffSq );
	  result = -two*tDiff*frequencySq_*expTerm;
	  break;
	}
    }
  }
};

#endif
