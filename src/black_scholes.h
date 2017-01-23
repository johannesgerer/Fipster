#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <cmath>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/distributions/normal.hpp>
#endif

#include "exceptions.h"

namespace fipster {
		
	namespace _black_scholes {

	using namespace boost::math;

//http://en.wikipedia.org/wiki/Black%E2%80%93Scholes#Black.E2.80.93Scholes_formula
struct black_scholes{
	const static bool n_eq_1;
	double S, K, T, volatility, dividend, interest;
	//helpers
	double dplus,dminus;
	normal s;
	unsigned int n;

	black_scholes(){};

	black_scholes(double S, double K, double T, double interest, double volatility, double dividend)
		:S(S),K(K),T(T),volatility(volatility),dividend(dividend),interest(interest),n(1)
	{}

	/** adjust the dividend and volatility to be able to model an option on the geometric average of
		n assets.
		*/
	black_scholes(double S, double K, double T, double interest, double volatility, double dividend,  double correlation, int n)
		:S(S),K(K),T(T),volatility(volatility),dividend(dividend),interest(interest),n(n)
	{
		if(n!=1){
			double f=1.0/n;
			this->dividend += (1-correlation)*(1-f)*volatility*volatility*0.5;
			this->volatility*=sqrt( (1-correlation)*f+correlation);
		}
	}
		
	void update_ds(double S){
		double nenner=volatility*sqrt(T);
		dplus=(log1p((S-K)/K)+(interest-dividend+volatility*volatility*0.5)*T)/nenner;
		dminus=dplus-nenner;
	}

	double call(double S){		
		update_ds(S);
		return cdf(s,dplus)*S*exp(-dividend*T)-cdf(s,dminus)*K*exp(-interest*T);
	}
	double put(double S){		
		update_ds(S);
		return -cdf(s,-dplus)*S*exp(-dividend*T)+cdf(s,-dminus)*K*exp(-interest*T);
	}
	double callDelta(double S){
		update_ds(S);
		return exp(-dividend*T)*cdf(s,dplus);
	}
	double putDelta(double S){
		update_ds(S);
		return exp(-dividend*T)*(cdf(s,dplus)-1);
	}
	double callGamma(double S){
		update_ds(S);
		return exp(-dividend*T)*pdf(s,dplus)/volatility/S/sqrt(T);
	}
	double putGamma(double S){
		return callGamma(S);
	}
	static void test()
	{{
		black_scholes bs(0,100,0.25,0.02,0.2,0.1);
		double S=10;
		const int N=6;
		double v[N]={1.6656018380511777e-120,
		89.748148798984914,
		3.8910939812958383e-119,
		-0.97530991202833262,
		9.0346688971466682e-118,
		9.0346688971466682e-118};
		double (black_scholes::*m[N])(double) ={&black_scholes::call,
												&black_scholes::put,
												&black_scholes::callDelta,
												&black_scholes::putDelta,
												&black_scholes::callGamma,
												&black_scholes::putGamma};
		for(int i=0;i<N;i++)
			FIPSTER_ASSERT((bs.*m[i])(S)==v[i]);
	}{
		black_scholes bs(0,100,0.25,0.02,0.2,0.1,0.5,10);
		double S=10;
		const int N=6;
		double v[N]={3.7942844439466054e-217,
89.77006860297864,
1.6089851348458302e-215,
-0.97311793162895988,
6.800006619321876e-214,
6.800006619321876e-214};
		double (black_scholes::*m[N])(double) ={&black_scholes::call,
												&black_scholes::put,
												&black_scholes::callDelta,
												&black_scholes::putDelta,
												&black_scholes::callGamma,
												&black_scholes::putGamma};
		for(int i=0;i<N;i++)
			FIPSTER_ASSERT((bs.*m[i])(S)==v[i]);
	}}
};

	}
	using _black_scholes::black_scholes;

}
