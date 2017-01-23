#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <map>
#include <array>
#include <memory>
#endif

#include "spacetime_field.h"
#include "configured_objects.h"
#include "expectation_value_config.h"


namespace fipster { namespace expectation_values {

	
	using namespace std;

	// #################################################
	// ##############   PHASE 2:	      		############
	// ##############   Node Factory        ############
	// #################################################

		/** conditional expectation of future timeslices of
				spacetime_fields. if given it adds the value of the hedging
				possition at future time (see expect_future_t), except if
				time==future.time (see CodeRefH)
		 */
		struct node_factory :spacetime_field<J_SGS,expect_future_t>
	{
		config_ptr config;
		shared_const_objs_t objs;

		//Constructor #########################################################	
		node_factory(config_ptr config,shared_const_objs_t objs);

		//returns the actual step_size, taking into account, that
		//the stepsizes do not have to be aligned with the 
		//stopping time limits!
		// btime calculate_step_size_and_theta(btime time,double& theta, bool& strang_symmetrization);

		/** spacetime_field: setup_node ####################################
		*/
		sender_ptr inner_setup(const arg_t& arg);
		bv_sender_ptr boundary_setup(const arg_t& arg);
		
		current_step calculate_current_step(const arg_t& arg);
	};

	
	/** pair consisting of the field_arg_t and a pointer to the node_factory
		of the expectation value for which the 
		inhomogeneous part of the boundary conditions should be calculated
	*/
		typedef pair<ev_field_arg
								 ,shared_ptr<node_factory> >
		bc_inhomogentities_arg_t;

}

	void time_limitation(btime time, btime last, string c);

	template<typename A>
	A time_supported(A ot, btime time, string c)
	{
		if(!ot) BOOST_THROW_EXCEPTION(runtime_error(c+": "+toS(time)+
																							 " is no supported time"));
		return ot;
	}
}
