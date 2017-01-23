//OLD UNUSED

#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <memory>
#include <tuple>
#include <map>
#include <exception>
#include <utility>
#include <algorithm>
#include <boost/throw_exception.hpp>
#endif


#include "tbb_flow_graph.h"
#include "type_traits.h"
#include "rooted_graph.h"
#include "named_counter.h"

namespace fipster { namespace _flow_graph {

	
	template<class arg_t,class result_t>
	struct runtime_join_node : 
		graph_node,
		sender<shared_ptr<const vector<shared_ptr<const result_t>>>>
	{
		//Forward sender<T> interface to the first ouput_port
		bool register_successor( successor_type &r ){
			return output_port<0>(my_node).register_successor(r);
		}
		bool remove_successor( successor_type &r ){
			return output_port<0>(my_node).remove_successor(r);
		}


		typedef shared_ptr<const pair<shared_ptr<const result_t>,arg_t>> input_t;  
		typedef multioutput_function_node<input_t,tuple<output_type>> node_t;

		typedef vector<shared_ptr<const result_t>> result_vector_t;
		typedef shared_ptr<result_vector_t> editable_result;		
		typedef shared_ptr<runtime_join_node> node_ptr;

		typedef multimap<arg_t,int> arg_map_t;

		//fields
		node_t my_node;

		//ctor
		runtime_join_node(const arg_map_t& arg_map)
			: my_node(rooted_graph::get(),serial,body_t(arg_map))
		{}

		//body struct
		struct body_t{

			//fields
			editable_result my_result;///\todo shared_ptr<vector<>> ?
			const arg_map_t arg_map; //maps argument value to result slots

			//ctor
			body_t(const arg_map_t& arg_map)
				:arg_map(arg_map),
				my_result(new result_vector_t)
			{
				my_result->reserve(max_element(begin(arg_map),end(arg_map),
					[](const arg_map_t::value_type& x, const arg_map_t::value_type& y ){
						return x.second < y.second;
				})->second);
			}

			void operator()(const input_t &v, typename node_t::output_ports_type &p)
			{
				//##################  PERFORMANCE critical  ##########################
				///\todo alles
				/*my_result->push_back(v->first);

				if(my_result->size()==n){ ///\todo if alles da?
					output_type result=my_result;
					my_result.reset();
					get<0>(p).put(result);
				}*/

				//##################  PERFORMANCE critical  ##########################
			}
		};

		//static factory
		template<class F>
		static node_ptr create(const arg_map_t& arg_map,
								F get_sender){

			//create new node
			node_ptr node(new runtime_join_node(arg_map));

			//register in the graph
			rooted_graph::get().add_node("multioutput_function_node",node);


			for(auto it=begin(arg_map);it!=end(arg_map);it++)
				make_edge<input_t>(*get_sender(it->first),node->my_node);

			named_counter::increment("edges",arg_map.size());

			return node;
		}

		////conversion
		//operator (){
		//sender_r t(){
		//	return sender_ptr(&output_port<0>(*this));///\todo:shared_ptr muss durch refernecen oder so ersetzt werden
		//}


	};
}

using _flow_graph::runtime_join_node;

}

