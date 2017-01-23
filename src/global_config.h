#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <string>
#include <boost/property_tree/ptree.hpp>
#endif

#include "exceptions.h"
#include "configured_objects.h"

namespace fipster { namespace configuration {

using namespace std;
using boost::property_tree::ptree;

/** Contains global configuration (e.g. to configure singletons)
*/
struct global_config{ 
	bool test_splitting;
	double test_splitting_tolerance;
	bool test_tridiagonally_split_grid_operator;
	bool test_bc_inversion;
	bool test_A_o;
	bool test_grid_iterators;
	bool initialized;
	bool test_boundary_condition_iterator;
	double years_per_btick;
	bool flow_graph_serial;
	bool test_gsi_gather_scatter;
	int flow_graph_processors;
	string output_path;
	int plot_count;

	global_config():initialized(false),plot_count(1){};

	static global_config& get(bool require_initialized=true){
		static global_config singleton;
		FIPSTER_ASSERT(!require_initialized || singleton.initialized);
		return singleton;
	}

	void reinit(const fpath base, const ptree& pt);
};

}

using configuration::global_config;

}

