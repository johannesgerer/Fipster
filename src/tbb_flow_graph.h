#pragma once

// works with arch version 4.2_20140122-1


#ifndef AUTOMATIC_PRECOMPILATION
#include <tbb/flow_graph.h>
#endif

#include "logging.h"

namespace fipster { namespace _flow_graph {

	using namespace tbb::flow;
	using namespace std;

	//! Makes an edge between a single predecessor and a single successor
//template< typename T >
//inline void make_edge2( sender<T> &p, receiver<T> &s ) {
//	cout<<&p<<" "<<&s<<endl;
//    p.register_successor( s );
//}

}}
