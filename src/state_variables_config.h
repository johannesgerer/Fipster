#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/property_tree/ptree.hpp>
#include <vector>
#include <algorithm>
#include <map>
#include <type_traits>
#include <comphash.h>
#include <pretty_printer.h>
#endif

#include "utility.h"

namespace fipster { namespace configuration {


	using namespace std;
	using namespace boost::property_tree;

	// ### Concept of state_variable_t ######
	/*
	struct state_variable_t{
		state_variable_t(const ptree& pt){...}
	};
	*/

template<class state_variable_t>
struct state_variables_config : comphashable<state_variables_config<state_variable_t> > {

	// typedef state_variable_t state_variable_t;
	typedef state_variables_config<state_variable_t> state_variables_config_t;


	//state variable vector
	vector<state_variable_t> sv_vec;

	//ctor
	state_variables_config(const ptree& pt)
	{
		//create a temporary associate container that keeps elements sorted
		map<int,state_variable_t> tc;

		//Loop through all childes that are stateVariable
		for(auto a=begin(pt);a!=end(pt);a++)
			if(a->first == "stateVariable")
			{
				//get state variable identifier
				int nr = a->second.get<int>("<xmlattr>.nr");

				auto ins = tc.insert (	
					make_pair(nr,state_variable_t(a->second)) );

				//if entry already existed
				if(!ins.second)
					BOOST_THROW_EXCEPTION(invalid_argument("duplicate entry for stateVariable "+toS(nr)));
			}else if(!is_xml_special(a->first) && a->first != "maxSubgrid")
				BOOST_THROW_EXCEPTION(runtime_error(a->first +" should be tag <stateVariable>"));

			//check, if largest stateVariable's number equals the total number of entries
		  if(tc.size() != (uint)tc.rbegin()->first)
				BOOST_THROW_EXCEPTION(invalid_argument("boundary conditions not specified for all state variables"));

			//insert the sorted elements into the main storage
			for(auto a=begin(tc);a!=end(tc);a++)
				sv_vec.push_back(a->second);
	}

	template<class T> void my_combine(T& t) const{
		t(&state_variables_config_t::sv_vec);
	}

};


}}
