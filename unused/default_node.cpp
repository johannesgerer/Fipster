unused 

#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#endif

#include "tbb_flow_graph.h"

namespace fipster { 
	namespace _flow_graph {

	using namespace std;
	using namespace tbb::flow;

		template<class T>
		struct default_node_t : continue_receiver {
			runtime_join_node& my_parent;
			//ctor
			continue_port_t(runtime_join_node& my_parent)
				:my_parent(my_parent)	{}
			
			task * execute(){ return my_parent.complete(); }
		};
	}
	using _flow_graph::default_node;
}
