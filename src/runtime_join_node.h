#pragma once

//todo gcc: theoretisch könnte man die auto_commecting_body für
//argumente spezialiseren, die als 'vector' markiert sind und dann
//automatisch den runtime_join_node davor schalten.

#ifndef AUTOMATIC_PRECOMPILATION
#include <memory>
#include <map>
#include <exception>
#include <utility>
#include <algorithm>
#include <boost/throw_exception.hpp>
#include <tbb/spin_mutex.h>
#include <tbb/atomic.h>
#endif



#include "tbb_flow_graph.h"
#include "type_traits.h"
#include "rooted_graph.h"
#include "named_counter.h"
#include "exceptions.h"


namespace fipster { 
	namespace _flow_graph {

	using namespace std;
	using tbb::spin_mutex;
	using tbb::task;
	using namespace tbb::flow;

	/**

		 T must be default initializable and copyassignable

	this class provides a new type of join node, that joins at
	runtime.

	*/
	template<class T>
  struct runtime_join_node : graph_node, 
		function_node<nothing,shared_ptr<const vector<T> > >::fOutput_type
	{
		typedef vector<T> vector_t;
		typedef typename vector_t::size_type size_type;

		struct continue_port_t : continue_receiver {
			runtime_join_node& my_parent;
			//ctor
			continue_port_t(runtime_join_node& my_parent)
				:my_parent(my_parent)	{}

			task * execute(){ return my_parent.complete(); }
		};

    /** \brief input ports for the results 
		 */ 
		struct input_port_t : receiver<T>
		{
			using input_type = typename input_port_t::input_type;

			bool closed;	
			spin_mutex input_port_mutex;
			runtime_join_node& my_parent;
			size_type storage_location;

			//ctor
			input_port_t(runtime_join_node& my_parent, size_type i)
				:closed(false),my_parent(my_parent),storage_location(i)	{};

			task * try_put_task( const input_type &t );
			void reset_receiver(reset_flags f = rf_reset_protocol)
			{ FIPSTER_ASSERT(0);/*not implemented*/	}

			~input_port_t(){
				//BOOST_THROW_EXCEPTION(runtime_error("input_port destructor: is this node a successor of any other node?"));
			}
			
		};	

		continue_port_t continue_port;
		vector<unique_ptr<input_port_t> > input_ports;
		shared_ptr<vector_t> output;
		size_type number_closed;
		spin_mutex join_node_mutex;

		input_port_t& get_pushed();

		task * input_port_closed();

		task* complete(){
			auto r = make_shared<vector_t>(input_ports.size());
			swap(r,output);
			return this->successors().try_put_task(r);
		}

		//ctor
		runtime_join_node(graph& g) 
			: graph_node(g),continue_port(*this)
			, output(make_shared<vector_t>(0))
			, number_closed(0) {}

		void reset_node(reset_flags f=rf_reset_protocol)
		{ FIPSTER_ASSERT(0);/*not implemented*/	}
	};

		//derive runtime_join_node type from container of
		//shared_ptr<sender<input_t>>
		template<class T> struct from{ 
			typedef runtime_join_node<typename
						 		T::value_type::element_type::output_type> rjn; };

	/** \brief static factory
			\tparam sender_t should be a container of shared_ptr<sender<input_t>>-Objects
		*/
		template<class T>
		static shared_ptr<typename from<T>::rjn> 
				new_runtime_join_node(T& senders)
		{
			typedef typename from<T>::rjn rjn;

			//create new node
			auto& g = rooted_graph::get();
			auto node = shared_ptr<rjn>(new rjn(g.tbb_graph));
			auto size = senders.size();

			if(size>0)
				for(auto& i : senders)
					make_edge(*i,node->get_pushed());
			else //connect to root_node
				make_edge(g.root_node, node->continue_port );
			
			
			named_counter::increment("edges",min((decltype(size))1,size));
			g.add_node("runtime_join_nodes",node);
			return node;
		}
		
	template<class T>
  typename runtime_join_node<T>::input_port_t&  runtime_join_node<T>::get_pushed()
	{
		input_ports.push_back(unique_ptr<input_port_t>
													(new input_port_t(*this,output->size())));
		output->push_back(T());
		return *input_ports.back();
	}


	template<class T>
  task*  runtime_join_node<T>::input_port_closed()
	{
		spin_mutex::scoped_lock lock(join_node_mutex);

		number_closed++;

		if(number_closed==input_ports.size()){
			for(auto& i : input_ports)	i->closed=false;

			number_closed=0;
			auto r = make_shared<vector_t>(input_ports.size());
			swap(r,output);
			return this->successors().try_put_task(r);
		}else
			return 0;
	}
		
	//  fyi: no race possible. 'closed' is only changed to 'true' by
	//  one thread, only in this function, and only if it was
	//  'false'. and the 'input_port_closed' function does not read
	//  the 'closed' values.

	//! Overwriting (it seems dangurous to overwrite try_put, as
	// many times try_put_task is called directly!?
	// [[file:/usr/include/tbb/flow_graph.h::task%20*res%20%3D%20try_put_task(t)%3B][file:/usr/include/tbb/flow_graph.h::task
	// *res = try_put_task(t);]]
	template<class T>
  task * runtime_join_node<T>::input_port_t::try_put_task( const input_type &t )
	{
		{
			spin_mutex::scoped_lock lock(input_port_mutex);
			if(closed)
				return 0;
			closed=true;
		}
		(*my_parent.output)[storage_location]=t;
		task* res = my_parent.input_port_closed();
		if(!res) res = tbb::flow::interface9::SUCCESSFULLY_ENQUEUED;

		return res;
	}

	 
}

using _flow_graph::new_runtime_join_node;
}
