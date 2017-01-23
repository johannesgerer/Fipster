#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/property_tree/ptree.hpp>
#include <string>
#include <vector>
#include <set>
#endif

#include "set_of_times.h"
#include "spacetime_field.h"
#include "configured_objects.h"
#include "grid.h"
#include "option.h"

namespace fipster { namespace _solution {  

	using namespace std;
	using namespace boost::property_tree;

	struct solution{

		//configuration parameters
		grid_ptr grid;
		const string id, grid_name;
		ptree fields_tree;
		vector<referenced_field_t> fields;
		boost::optional<btime> fixed_plot_time;
		const combB times;
		bool active;
		boost::optional<string> output_as_geometric_average_grid_name,interval_from;
		vector<state_t> single_values;
		set<btime> not_done;
		bool show_hedging_pos;
		   
		//solution storage (unused)
		//spacetime_field<J_SGS>::result_sptr_t storage;

		//other configured objects
		shared_const_objs_t objs;

		//constructor from XML parsed property tree ################
		solution(const ptree& pt,shared_const_objs_t objs);

		void setup_nodes();
		
		void print();
	};
}

using _solution::solution;

ostream& operator<<(ostream& out,const solution& sol);

}
