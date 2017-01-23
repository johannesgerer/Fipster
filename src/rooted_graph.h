#pragma once


#ifndef AUTOMATIC_PRECOMPILATION
#include <forward_list>
#include <memory> 
#endif

#include "named_counter.h"
#include "tbb_flow_graph.h"
#include "logging.h"

namespace fipster { namespace _flow_graph {


	//########################
	// Rooted graph
	//########################
	// A graph that contains also a root_node.  The nodes who have open
	// receivers should be connected to it to be able to start them all
	// collectively with the command: root_node.try_put(continue_msg);
	struct rooted_graph {
		
		//typedefs
		
		typedef shared_ptr<graph_node> graph_node_ptr; 

		typedef forward_list<graph_node_ptr> graph_node_list;

		graph tbb_graph;
		broadcast_node<continue_msg> root_node;
		//assures the survival of nodes (that otherwise are only passed
		//around tbb's flow_graph as pointers)
		map<string,graph_node_list> nodes;
		

		//add_node (to to list of nodes, according to string s) ##########
		void add_node(string s,graph_node_ptr node );

		//conversion #################################################
		operator graph&(){
			return tbb_graph;
		}

		//constructor #################################################
		rooted_graph():root_node(tbb_graph){};

		//singleton #################################################
		static rooted_graph& get(){
			static rooted_graph singleton_instance;
			return singleton_instance;
		}

		//sink_receiver class definition #################################
		//A receiver that stores in the received value in the result_sink
		template<class T>
		struct sink_receiver : public receiver<T>,public graph_node {

			//reference to the result
			T& target;
			//ctor
			sink_receiver(T& target,graph& g) 
				: target(target), graph_node(g)
			{}

			bool try_put(const T& v){
				thread_logger()<<"SINK RECEIVED"<<endl;
				target = v;
				return true;
			}
		};

		//#################################################
		//creates a sink_receiver that will fill the target on receive
		template<class T,class sender_t>
		void attach_result_to_sink(shared_ptr<T>& target, sender_t sender){
			auto receiver = new sink_receiver<shared_ptr<const T>>(target,tbb_graph);

			add_node("sink_receivers",graph_node_ptr(receiver));

			//connect the receiver to the sender
			//cout<<"sol sink"<<endl;
			make_edge(*sender,*receiver);
			named_counter::increment("edges");
		}

		void start_wait_for_all();
	};

}

using _flow_graph::rooted_graph;

}
