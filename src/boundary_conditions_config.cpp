#include "precompiled.h"
#include "boundary_conditions_config.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/property_tree/ptree.hpp>
#include <boost/throw_exception.hpp>
#include <exception>
#endif

#include "exceptions.h"

namespace fipster { namespace _boundary_conditions {
	
	// #################################################
	// ##############   PHASE 1:			############
	// ##############   CONFIGURATION       ############
	// #################################################

	void bc::init(const ptree& pt)
	{ 
		string stype = pt.get<string>("type");
		string inhomValue = pt.get("inhomValue",stype); //optional
		 
		if(stype=="vonNeumann")
			type=vonNeumann;
		else if(stype=="Dirichlet")
			type=Dirichlet;
		else
			BOOST_THROW_EXCEPTION(invalid_argument(stype+" is no valid Boundary Condition Type"));

		
		if(inhomValue=="vonNeumann")
			inhom_value=vonNeumann;
		else if(inhomValue=="Dirichlet")
			inhom_value=Dirichlet;
		else if(inhomValue=="zero")
			inhom_value=zero_value;
		else
			BOOST_THROW_EXCEPTION(invalid_argument(inhomValue+" is no valid Boundary Condition inhomValue"));

		if(inhom_value!=zero_value){
			
			string inhomSource = pt.get<string>("inhomSource");

			if(inhomSource=="exerciseValue")
				inhom_source=exercise_value;
			else if(inhomSource=="futureValue")
				inhom_source = future_value;
			else
				BOOST_THROW_EXCEPTION(invalid_argument(inhomSource+" is no valid Boundary Condition inhomSource"));
		}

		//CodeRefB
		if(inhom_source == future_value)
			time_steps = pt.get("timeSteps",0);
		else
			time_steps = 0;

	}



	state_variable_t::state_variable_t(const ptree& pt)
	{
		lowerupper[lower_ind].init(pt.get_child("lower.<xmlattr>"));
		lowerupper[upper_ind].init(pt.get_child("upper.<xmlattr>"));
	}



	config_t::config_t(const ptree& pt)
		:state_variables_config_t(pt)
	{
		tolerance = pt.get<double>("<xmlattr>.tolerance");
	}

}}
