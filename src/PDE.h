#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/lexical_cast.hpp>
#endif

#include "finite_difference_weights.h"
#include "integer_types.h"
#include "configured_objects.h"
#include "configuration_parser.h"

namespace fipster { 
	
	/** The computation of expectation values is achieved via solution
		of a corresponding Feynman-Kac formula, which is a
		terminal-boundary value linear complementarity problem described
		by:

		\code
		min ( f_t + A f , f - g ) = 0
		\endcode
		with known terminal value (and boundary) values:
		\code
		f(T)
		\endcode

		The above choice of the relative sign of A and the (positive) time
		derivative will dictate a negative M-Matrix property for the
		discretization of A.

		This file contains the PDE coefficient functions of A.

	The (reduced) concept for a PDE class:
	\code
	struct PDEconcept: PDE_base {
		template<class T>
		FORCE_INLINE double operator()(T& state) const{
			//return the coefficient of the non-derivative term
		};

		template<class T,class A>
		FORCE_INLINE double operator()(A& pdIndex,T& state,uint D) const{
			//calculate the coefficient of the d/(dx_i)^n_i/(dx_j)^n_j... term,
			//	where n_i are the entries of pdIndex
		}
	};
	\endcode
	*/
	namespace finite_difference_weights {

		/** Black Scholes isometric

				Calculates the expectation of continuously discounted payoffs,
				where all quantities (payoff and result) are in undiscounted
				form. (Theorem 6.4.3 Shreve II)
				
				For actual black scholes (which calculates the expectation
				value under the risk free measure, i.e. where the stock with
				dividends reinvested grows with the risk free rate) set:

				  discounting_rate = risk_free_rate
   				drift = risk_free_rate - dividend
						
				For the conditional expectation with respect to a
				geom. Brownian motion set:
  
          discounting_rate = 0
 */
		struct BSiso: PDE_base {

			double discounting_rate, volatility, correlation,drift;
			// double interest, volatility, correlation, dividend;

			BSiso(string fdw,const ptree& pt,string type):PDE_base(type){
				discounting_rate = pt.get<double>("discountingRate");
				if(discounting_rate != 0)
					cout<<"careful: discounting_rate only applies to the expectation, thus it is not suitable to calcuate higher moments of discounted quantites by solving this pde with undiscounted final values"<<endl;
				volatility = pt.get<double>("volatility");
				correlation = pt.get("correlation",0.0);
				drift = betterGet<double>(pt,"drift",fdw);
			}
		
			string comparison_string()const{
				using namespace boost;
				return	lexical_cast<string>(discounting_rate)	 +";" +
						lexical_cast<string>(volatility) + ";" +
						lexical_cast<string>(drift)+ ";" +
						lexical_cast<string>(correlation);
			}

			template<class T>
			FORCE_INLINE double operator()(T& state) const
			{return -discounting_rate;};

			template<class T,class A>
			FORCE_INLINE double operator()(const A& pdIndex,T& state,uint D) const{

				#define state_(a) state[*a-1]
				#define multiindex_(a) pdIndex(*a-1)

						switch(pdIndex.sum()){
						case 0: return this->operator()(state);
						case 1:
							for (uint i = 1; i <= D; ++i)
								if (multiindex_(&i) == 1)
									return drift* state_(&i);
						case 2:
							for (uint i = 1; i <= D; ++i) {
								if (multiindex_(&i) == 2) {
									// Computing 2nd power 
									double t = state_(&i);
									return volatility * volatility / 2 * (t*t);
								} else if (multiindex_(&i) == 1) {
									double ret_val = volatility * volatility * correlation * state_(&i);
									for (uint j = i + 1; j <= D; ++j) {
										if (multiindex_(&j) == 1) {
											return ret_val * state_(&j);
										}
									}
								}
							}
						default:
							return 0.;
						}
		
				#undef state
				#undef multiindex
			}
		};

		struct BSiso_geom_average : BSiso {

			string finiteDifferenceWeights_name;

			BSiso_geom_average(string fdw, const ptree& pt, objs_t objs )
				: BSiso(dynamic_cast<const BSiso&>
								(*betterAt
								 (objs.fdweights_configs
									,betterGet<string>(pt,"finiteDifferenceWeights",fdw+"_BSIsoGeom")
									,fdw+"_BSIsoGeom")->PDE.get()))
				, finiteDifferenceWeights_name
					(betterGet<string>(pt,"finiteDifferenceWeights",fdw+"_BSIsoGeom"))
			{

				auto pt_D=pt.get_optional<uint>("dimension");
				
				uint D;
				if(pt_D) D=*pt_D;
				else{
					auto grid=pt.get<string>("dimension.<xmlattr>.fromGrid");
					D = betterAt(objs.grids,grid,fdw+"_BSIsoGeom")->D;
				}
				double f=1.0/D;

				assert(!"checked again");
				drift -= (1-correlation)*(1-f)*volatility*volatility*0.5;
				//was : dividend += (1-correlation)*(1-f)*volatility*volatility*0.5;
				volatility*=sqrt( (1-correlation)*f+correlation);
			}
		};

	}
}
