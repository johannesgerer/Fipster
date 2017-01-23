#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <forward_list>
#include <boost/exception/all.hpp>
#endif

#include "rooted_graph.h"
#include "tbb_flow_graph.h"
#include "logging.h"
#include "named_counter.h"
#include "global_config.h"

namespace fipster { namespace _flow_graph {

	using namespace tbb;

	void rooted_graph::start_wait_for_all()
	{
		thread_logger()<<"GRAPH  | Start and wait for all"<<endl;

		root_node.try_put(continue_msg());
		tbb_graph.wait_for_all();
	}

	void rooted_graph::add_node( string s,graph_node_ptr node )
	{
		named_counter::increment(s);
		rooted_graph::get().nodes[s].
			push_front(node);
	}

	
}


}
