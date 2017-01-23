#include "precompiled.h"
#ifndef AUTOMATIC_PRECOMPILATION
#include <pretty_printer.h>
#endif
#include "configured_objects.h"
#include "boundary_conditions.h"
#include "finite_difference_weights.h"
#include "spacetime_field.h"
#include "timestepping_operator.h"
#include "option.h"
#include "expectation_value_node_factory.h"

namespace fipster { 
namespace configuration {


	template<class T>
	void print_map_size(string name,const T& m)
	{
		cout<<"Cache size of '"+name+"':" +toS(m.cache.size())<<endl;
	}
	
	template<class T>
	void print_map_sizes(string name,const T& m){
		for(auto it=begin(m);it!=end(m);it++)
			print_map_size(it->first,*it->second);
	}

	void configured_objects::print_all_map_sizes(){
		print_map_size("boundary_conditions",boundary_conditions::singleton());
		print_map_size("fdweights",fdweights::singleton());
		print_map_sizes("spacetime_field",spacetime_fields);
		print_map_sizes("options",options);
		print_map_sizes("expectations",expectations);
		
		// for(auto i:expectations)
		// 	for(auto j: i.second->cache)
		// 		cout << j.first.time << ": " << j.second.get() <<" " << j.first.carg << endl;

		{
			using namespace _timestepping_operator;
			print_map_size("splitter",splitter::singleton());
			print_map_size("timestepping_operator",timestepping_operator::singleton());
		}
		
	}
}}
