#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <memory>
#include <boost/property_tree/ptree.hpp>
#include <array>
#include <comphash.h>
#endif

#include "state_variables_config.h"
#include "boundary_enum.h"

namespace fipster { 
	namespace _boundary_conditions {

	using namespace std;
	using namespace boost::property_tree;

	// #################################################
	// ##############   PHASE 1:			############
	// ##############   CONFIGURATION       ############
	// #################################################

	enum inhomSource_t {  exercise_value, future_value};	
	enum type_t { zero_value, vonNeumann , Dirichlet };


	// #################### bc  ###########################
	//This object encodes the information needed to describe boundary 
	//conditions for one bounding plane of the hypercube. It will be 
	//used internally by the class "boudnary_conditions" which handles
	//the boundary conditions as a whole.
	struct bc : comphashable<bc> {
		
		inhomSource_t inhom_source;
		type_t type;
		type_t inhom_value;
		int time_steps;

		void init(const ptree& pt);

		template<class T> void my_combine( T& t ) const
		{
			t(	&bc::time_steps,
				&bc::inhom_value,
				&bc::type,
				&bc::inhom_source);
		}
	};


	// #################### wrapper type for one state variable's two bcs  ###########################
	struct state_variable_t : comphashable<state_variable_t> {


		array<bc,2> lowerupper;

		state_variable_t(const ptree& pt);

		template<class T> void my_combine(T& t) const
		{
			t(&state_variable_t::lowerupper);
		}
	};



	/** boundary_conditions::config_t 
	*/
	struct config_t : configuration::state_variables_config<state_variable_t> {
		config_t(const ptree& pt);
		double tolerance;
	};

	
}

}
