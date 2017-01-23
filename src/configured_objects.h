#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <unordered_map>
#include <memory>
#include <string>
#include <list>
#include <boost/filesystem.hpp>
#endif

#include "forward.h"

namespace fipster { 

	using namespace std;
	
	
namespace configuration {

	typedef boost::filesystem::path fpath;
	
	struct configured_objects{

		//typedefs
		typedef shared_ptr<_fields::spacetime_field<J_SGS>> field_ptr_t;
		typedef unordered_map<string,field_ptr_t> spacetime_fields_t;

		typedef shared_ptr<expectation_values::node_factory> exp_ptr_t;
		typedef unordered_map<string,exp_ptr_t> expectation_values_t;

		typedef shared_ptr<option_node_factory> opt_ptr_t;
		typedef unordered_map<string,opt_ptr_t> options_t;

		typedef unordered_map<string,
			shared_ptr<_boundary_conditions::config_t>> bcs_config_t;

		typedef list<_solution::solution> solutions_t;

		typedef unordered_map<string,
			shared_ptr<_grid::Grid>> grids_t;

		typedef unordered_map<string,
			shared_ptr<const finite_difference_weights::config_t>> fdweights_configs_t;

		//associative containers of parsed objects
		expectation_values_t expectations;
		options_t options;
		spacetime_fields_t spacetime_fields;
		bcs_config_t bcs_configs;
		grids_t grids;
		solutions_t solutions;
		fdweights_configs_t fdweights_configs;

		void print_all_map_sizes();
	};

}
using configuration::configured_objects;
typedef shared_ptr<configured_objects> shared_objs_t;
typedef const configured_objects& objs_t;
typedef shared_ptr<const configured_objects> shared_const_objs_t;
}
