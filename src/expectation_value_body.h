#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <iostream>
#include <map>
#include <utility>
#include <array>
#include <boost/optional.hpp>
#endif

#include "spacetime_field.h"
#include "auto_connecting_body.h"
#include "expectation_value_config.h"
#include "timestepping_operator.h"
#include "grid_fields.h"

namespace fipster { namespace expectation_values {


	using namespace std;

	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################

	//body_t ######################################################################
	/** Body for the main timestepping work
	//\todo warum nicht J_SGS?
	*/
	struct body_t : auto_connecting_body
	     <	body_t,
					//result
					discretized_space_field<1>,
					//future 0:
					discretized_space_field<1>, //Exercise value
					//future 1:
					discretized_space_field<1>, //Future value
					//future 2:
					timestepping_op_t,			//vector of timestepping operators
					//future 3:
					boundary_field_t<J_SGS>,		//boundary inhomogenities
					//future 4:
					boundary_field_t<J_SGS>		//future boundary inhomogenities
		>
	{ 
		string id;
		current_step cur_step;
		ev_field_arg fieldargs;
		string meta_info; //stored to be put into result
		double vola; //to be used for guessing higher powers at the
								 //boundary
		btime time;

		body_t(const config_ptr& config,const current_step& cur_step
					 ,const ev_field_arg& fieldargs, const string& meta_info
					 ,double vola, btime);

		result_sptr_t operator()(const input_t& futures);

	};
		
	//######################################################################
	/** Body for the boundary value extraction
	*/
	struct boundary_values_body_t : 
		auto_connecting_body<	boundary_values_body_t,	//body
								boundary_field_t<J_SGS>,		//result
								discretized_space_field<J_SGS>,	//future 1 (interior values)
								B_t,					//future 2
								boundary_field_t<J_SGS>		//future 3
							>
	{
		result_sptr_t operator()(const input_t& futures);		
	};



}


}


