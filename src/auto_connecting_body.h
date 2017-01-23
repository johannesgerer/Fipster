#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <tuple>
#include <type_traits>
#endif

#include "type_traits.h"
#include "rooted_graph.h"
#include "named_counter.h"
#include "tbb_flow_graph.h"
#include "logging.h"

// Documentation:
/*
The class template auto_connecting_body (automatically connecting body) serves two purposes:

1.  It is a factory for graph_nodes that creates readily 
	connected node according to the "Body". 
	Its static member ::create(...) will be return a shared pointer to a 
		* continue_node for nullary bodies
		* function_node for unary bodies
		* join_node + function_node for n>1-ary bodies

	The concepts for bodies as described in the TBB reference apply correspondingly, but
	the result and argument types (for n>0-ary bodies) will be transformed into shared_ptrs

2.	It provides and simplifies the use of types connected with the programming
	of new bodies. To make use of this, simply derive the new body_t from the
	corresponding auto_connecting_body:

	
	//Example
	struct body_t1 : struct auto_connecting_body<body_t,result_t,future_t>  {
		
		//unary: typedef shared_ptr<future_t> input_t;
		//nullary: no named(!) futures argument for "operator()";
		//n-ary (Tuple): typedef typename apply<shared_ptr,future_t>::type input_t; 

		result_sptr_t operator()(const input_t &arg){
			//The result can be created using:
			auto result = make_shared<result_t>(...);
			return result;
		}
	};
	
	
*/


namespace fipster { namespace _flow_graph {

//########################
// make::edges
//########################
// This is utility function, that connects the first n senders contained
// in sender_tuple to the input_ports of the joinnode
// via recursive template expansion
template<int n> struct make{
	template<class A,class B>
	static void edges(A& sender_tuple, B& joinnode){
		make_edge(*get<n-1>(sender_tuple),input_port<n-1>(joinnode));
		make<n-1>::edges(sender_tuple,joinnode);
	}
};

template<> struct make<0>{
	template<class A,class B>
	static void edges(A& sender_tuple,B& joinnode){}
};



//########################
// Worker Base class
//########################
template<class result_tpl>  //The result type of the worker
struct auto_connecting_body_base{

	//########################
	//Type definitions
	//########################
	
	//The result type of the worker
	typedef result_tpl result_t;
	//Shared pointer type to the result. This is actually returned by the operator()
	typedef shared_ptr<const result_t> result_sptr_t;

	//Shared pointer type to create the result:
	typedef shared_ptr<result_t> editable_result;

};

//(Not used anymore) this was needed for (older) runtime_join_node
////########################
//// Worker Base class Specializations
////########################
//template<class arg_t,class result_t>
//struct auto_connecting_body_base<pair<arg_t,shared_ptr<const result_t>>>{
//	typedef result_t result_t;
//	typedef shared_ptr<const result_t> editable_result;
//	typedef pair<arg_t,shared_ptr<const result_t>> result_sptr_t;	
//};

//this is used for partial template specialization
struct nothing {};



//########################
//auto_connecting_body
//########################
template<class body_t, class result_t,  OPTIONAL_TUPLE_PARAMS>
struct auto_connecting_body : auto_connecting_body_base<result_t>
{
	//########################
	//Type definitions
	//########################

	typedef typename auto_connecting_body::result_sptr_t result_sptr_t;

	//Tuple of Argument types
	typedef tuple<TUPLE_ARGS> future_tupel_t;	

	//Tuple of shared pointers to Argument types. This is actually passed to operator()
	typedef typename apply_meta_apply<shared_ptr,add_const, future_tupel_t>::type input_t; 

	//Compile time checks
	static const int n = tuple_size<future_tupel_t>::value;
	static_assert(n>1,"###########  Don't use the tuple version for null-/unary workers ###########");
#ifdef _MSC_FULL_VER
	static_assert(n<=10,"###########  Not implemented for tuples with size larger than 10 ###########");
#endif

	//Tuple of pointers to senders. From these senders the args are received
	//(raw pointer const type is used instead of reference to 
	//be able to compare the target's address)
	typedef typename apply_meta<add_const,
		typename apply_apply_meta<add_pointer,tbb::flow::sender,input_t>::type
	>::type sender_t;


	//########################
	// Custom node object
	//########################
	struct node_t : public function_node<input_t,result_sptr_t> {

		join_node<input_t> joinnode;

		void make_edges(const sender_t& sender_tuple){
			make<n>::edges(sender_tuple,joinnode);
			
		}

		//constructor
		node_t(rooted_graph& g,const body_t& body,const sender_t* sender_tuple)
			:	function_node<input_t,result_sptr_t>(g,serial,body), joinnode(g)
		{
			if(sender_tuple){
				//connect the tuples of senders to the receivers of the join_node
				make_edges(*sender_tuple);
				named_counter::increment("edges",n);
			}
			make_edge(joinnode,*this);
			named_counter::increment("edges");
		};

	};

	typedef shared_ptr<node_t> node_ptr;


	static node_ptr new_node(rooted_graph& g,const body_t& body,const sender_t* sender_tuple=nullptr)
	{
		auto node = node_ptr(new node_t(g,body,sender_tuple));
		g.add_node("joined_function_node",node);
		return node;
	}
};

//########################
//Unary specialization
//########################
template<class body_t, class result_t, class T0>
struct auto_connecting_body<body_t,result_t, NINE_DEFAULT_TUPLE_ARGS>
	: auto_connecting_body_base<result_t> {
	
	//########################
	//Type definitions
	//########################
	typedef typename auto_connecting_body::result_sptr_t result_sptr_t;

	//Shared pointer to future_t. This is actually passed to operator()
	typedef shared_ptr<const T0> input_t;

	//pointers to the sender. From this senders the arg is received
	//(raw pointer const type is used instead of reference to 
	//be able to compare the target's address)
	// \todo0: why?
	typedef sender<input_t>* const sender_t;
	typedef function_node<input_t,result_sptr_t> node_t;
	typedef shared_ptr<node_t> node_ptr;
	

	static node_ptr new_node(rooted_graph& g, const body_t& body,const sender_t* sender=nullptr)
	{
		//create new node
		auto node = node_ptr(new node_t(g,serial,body));

		if(sender)
			//connect the sender to the new node
			make_edge(**sender,*node);

		named_counter::increment("edges");
		
		g.add_node("function_nodes",node);
		return node;
	}
};




//########################
//Nullary Specialization
//########################
template<class body_t, class result_t>
struct auto_connecting_body<body_t,result_t TEN_DEFAULT_TUPLE_ARGS> 
			: auto_connecting_body_base<result_t>{
	
	//########################
	//Type definitions
	//########################
	typedef typename auto_connecting_body::result_sptr_t result_sptr_t;

	//continue_msg is actually passed to operator()
	typedef continue_msg input_t;

	typedef continue_node<result_sptr_t> node_t;
	typedef shared_ptr<node_t> node_ptr;


	static node_ptr new_node(rooted_graph& g, const body_t& body)
	{	
		auto node = node_ptr(new node_t(g.tbb_graph,body));

		//connect the root_node to the receiver of the continue_node
		make_edge(g.root_node,*node);
		g.add_node("continue_nodes",node);
		named_counter::increment("edges");

		return node;
	}

};




//########################//########################//########################
/** \brief Node factory for auto_connecting_bodies with futures
	\param sender either a sender<input_t>*const for unary bodies or a tuple of these.
*/
template<class body_t>
typename body_t::node_ptr 
	create_node(const body_t& body,const typename body_t::sender_t& sender)
{
	return body_t::new_node(rooted_graph::get(), body, &sender);
}

/** \brief Node factory for auto_connecting_bodies with no futures
*/
template<class body_t>
typename body_t::node_ptr 
create_node(const body_t& body)
{
	return body_t::new_node(rooted_graph::get(), body);
}


		/** body that merely returns an empty pointer
		 */
		template<class T>
		struct default_body_t : auto_connecting_body<
			default_body_t<T>, T>
		{
			typename default_body_t::result_sptr_t 
			operator()(const typename default_body_t::input_t&){
				return nullptr;
			}
		};
		
		/** node that merely sends an empty pointer
		 */
		template<class T>
		typename default_body_t<T>::node_ptr 
		default_node()
		{
			static auto n = create_node(default_body_t<T>());
			return n;
		}

}

using _flow_graph::auto_connecting_body;
using _flow_graph::create_node;
using _flow_graph::default_node;

}//end namespaces

