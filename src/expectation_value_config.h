#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/property_tree/ptree.hpp>
#include <boost/optional.hpp>
#include <memory>
#include <array>
#include <comphash.h>
#endif

#include "set_of_times.h"
#include "utility.h"
#include "forward.h"
#include "option.h"

namespace fipster { namespace expectation_values {


		using namespace std;

		// #################################################
		// ##############   PHASE 1:			############
		// ##############   CONFIGURATION       ############
		// #################################################
		using namespace boost::property_tree;

		struct details_t : comphashable<details_t>
									,stringifyable<details_t>
		{
			uint power;
			boost::optional<state_t> hedging;
			boost::optional<double> entropic_theta;
			double discounting_factor;

			
			template<class T> FORCE_INLINE void my_combine(T& t)const{
				t(&details_t::power
					,&details_t::hedging
					,&details_t::discounting_factor
					,&details_t::entropic_theta);
			}

			details_t():power(1),discounting_factor(1){};
		};

		struct expect_future_t : comphashable<expect_future_t>
													 ,stringifyable<expect_future_t>
		{
			btime field_time; // point of time of the field's value = x in E(x)
			plain_spacetime_field field;
			bool rannacher; //indicates, whether the future value is the the
											//bare x (with possibly discontinuous
											//derivatives)
			plain_spacetime_field exercise; //lcp constraint (optional, can be nullptr)
			details_t details;
			
			template<class T> FORCE_INLINE void my_combine(T& t)const{
				t(&expect_future_t::field_time
					,&expect_future_t::field
					,&expect_future_t::rannacher
					,&expect_future_t::exercise
					,&expect_future_t::details);
			}
			
			expect_future_t(btime field_time, btime stop
											,plain_spacetime_field exercise
											,boost::optional<double> entropic_theta
											,double discounting_factor);
		};
		
		typedef _fields::field_arg_t<expect_future_t> ev_field_arg;

		// information used to calculate the current step. (this
		// information will be determined from the node_factory's argument
		// consisting of the current time and expect_future)
		struct current_step{
			btime step_size;
			double theta;
			double timestep;
			bool strang_symmetrization;
			details_t now;
			details_t future;
			array<double,2> betas();
		};
		
		struct config_t{
		
			int rannacher_steps;
			string boundary_conditions_s,fd_weights_s,id;
			comb steps;
			btime rannacher_step_size;
		
			//actually used in body:
			double theta,rannacher_theta;
			bool strang_symmetrization,rannacher_strang_symmetrization;

			//ctor
			config_t(const ptree& pt);
		
			boost::optional<current_step>
			calculate_current_step( btime time, expect_future_t future);
			void check_arbitrage(shared_ptr<option_config_t> o_conf
													 ,shared_const_objs_t objs
													 , btime fut_time, btime arg_time);
		};
	
		typedef shared_ptr<config_t> config_ptr;

	}
	using expectation_values::ev_field_arg;
	using expectation_values::expect_future_t;
	
}
