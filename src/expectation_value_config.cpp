#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/throw_exception.hpp>
#include <boost/property_tree/ptree.hpp>
#include <string>
#include <pretty_printer.h>
#endif

#include "exceptions.h"
#include "configuration_parser.h"
#include "expectation_value_config.h"
#include "Math.h"
#include "global_config.h"

using namespace fipster;
using namespace fipster::expectation_values;

// #################################################
// ##############   PHASE 1:			############
// ##############   CONFIGURATION       ############
// #################################################
		
expect_future_t::expect_future_t(btime field_time, btime rannacher_time
																 , plain_spacetime_field exercise
																 ,boost::optional<double> entropic_theta
																 ,double discounting_factor)
	: field_time(field_time)
	, rannacher(field_time==rannacher_time)
		//apply rannacher to damp the oscillations caused by
		//discontinuities in the derivative of the payoff
	,exercise(exercise)
{
	details.entropic_theta=entropic_theta;
	details.discounting_factor = discounting_factor;
}

config_t::config_t( const ptree& pt )
	: id(betterGet<string>(pt,"<xmlattr>.id","expectation value"))
	, steps(betterChild(pt,"timeStepping",id),id+".timeStepping")
{

	boundary_conditions_s = pt.get<string>("boundaryConditions");
	fd_weights_s = pt.get<string>("finiteDifferenceWeights");

	strang_symmetrization = pt.get<bool>("strangSymmetrization");

	auto ran=betterChild(pt,"rannacherSteps",id);
	auto id2=id+".rannacherSteps";
	rannacher_step_size = getBtime(ran,"stepSize",id2);
	rannacher_steps = betterGet<int>(ran,"<xmlattr>.steps",id2);
	rannacher_theta = betterGet<double>(ran,"theta",id2);
	rannacher_strang_symmetrization = betterGet<bool>(ran,"strangSymmetrization",id2);

	theta = betterGet<double>(pt,"theta",id);

	if(theta<0 || theta>1)
		FIPSTER_THROW_EXCEPTION(runtime_error(toS(theta)+" is no valid Time Stepping theta value"));
		
}


