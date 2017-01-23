#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <map>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <memory>
#include <type_traits>
#endif

#include "tbb_flow_graph.h"
#include "named_counter.h"
#include "auto_connecting_body.h"

// Documentation:


namespace fipster { namespace _flow_graph {
		

/** \brief Base class for custom node factories classes.

	These factory classes are a convenient way to setup functors that 
	create graph_nodes and posses their own context. To add memoization to 
	a node factory class it has to derive from memoized_node_factory
	and comply to the following concept:
	\code
struct my_factory : memoized_node_factory<my_factory, arg_t, result_t {

	//obligatory:
	sender_ptr setup_node(const arg_t& arg){
		//Best practice: Use the free \ref create_node function:
		return create_node(auto_connecting_body::new_node(some_arg));
	};
	
	//optional: (overwriting)
	virtual sender_t get_node(const arg_t& arg){
		//Trivial (do nothing) example
		return node_factory::get_node(arg);
	}
};
	\endcode

	Then use the inherited "sender_t get_node(const arg_t& arg)" function to register the
	sender<T> with other nodes using through 
	auto_connecting_body, runtime_join_node,
	rooted_graph::attach_result_to_sink, tbb::flow::make_edges and similar

	\tparam derived_factory_t The derived_factory_t template argument is used for static polymorphism
	and, whats more important, as a means to be able to have different
	static instances of the class with distinct maps for memoization.
*/

template<class result_t>
		struct mtypes{
	typedef typename auto_connecting_body_base<result_t>::result_sptr_t
	result_sptr_t; 

	typedef shared_ptr<sender<result_sptr_t> > sender_ptr;
};

template<class derived_factory_tpl, class arg_tpl, class result_tpl>
struct memoized_node_factory{

	//Typedefs

	typedef derived_factory_tpl derived_factory_t; ///< The type of the derived factory
	typedef arg_tpl arg_t;						//Argument
	typedef result_tpl result_t;					//Result	
	typedef typename mtypes<result_t>::result_sptr_t result_sptr_t;	//Result pointer

	
	typedef typename mtypes<result_t>::sender_ptr sender_ptr;

	//this type
	typedef memoized_node_factory<derived_factory_t,arg_t,result_t>
		node_factory_t;
	typedef shared_ptr<node_factory_t> node_factory_ptr;

	//possible cache types
	typedef unordered_map<arg_t,sender_ptr,boost::hash<arg_t> > unordered_map_t;	
	typedef map<arg_t,sender_ptr> map_t;
	typedef unordered_map_t cache_t;	//cache type


	cache_t cache; ///< %A (hash) map of all the memorized nodes

	void decr(typename map_t::iterator& t){ --t; }
	template <class T> void decr(T& t){ }
	
	/**\brief "The" main memoization method.

		\returns either an existing node (\ref sender_t) corresponding to arg,
		or creates a new by calling the  derived_factory_t::setup_node(arg)	

		\note This method at first sight, does not have to be virtual, 
				as a static_cast could be employed. But it is needed, to
				allow runtime polymorphism for objects derived from derived_factory_t
				(like spacetime_field)
				
*/
	virtual sender_ptr get_node(const arg_t& arg){ 
		//Look for arg (key)
		auto range = cache.equal_range(arg);
		
		//if entry did not already exist, create it:
		if(range.first == range.second)
		{
			named_counter::increment("cache miss");

			if(cache.empty())
				decr(range.first);

			//cout<<boost::hash<decltype(arg)>()(arg)<<endl;

			//attention: recursion in setup_node means quasi-parallel access
			//to this section. 

			//\todo1 gcc: is this a problem due to "unordered_map::insert: If
			//rehashing occurs due to the insertion, all iterators are
			//invalidated." ?
			//http://en.cppreference.com/w/cpp/container/unordered_map/insert

			//\todo1 gcc: the hint does not do anything for the unordered_map (?)
			// todo: put hint back in: range.first or use second return
			// value of insert instead
			auto it = cache.insert(make_pair(arg, 
				static_cast<derived_factory_t*>(this)->setup_node(arg)
														));
			return it.first->second;
		}else{
			named_counter::increment("cache hit");
			return range.first->second;
		}
	}
};



/** \brief Singleton Memoized Node factory
	
	Derive from this class to add a static singleton
	version of \ref get_node (called \ref sget_node). This is needed/usefull for
	node factories, that need to setup (possibily interconnected)
	nodes, but the bodies ctors' arguments can be fully deduced
	from the memoized argument. I.e. there is no additional information
	required, which makes instances of such factories superflous, and thus
	a singleton is used.

	In addition to the <B>derived_factory_t concept</B> required at \ref memoized_node_factory,
	it needs <B>default constructor</B>.
*/
template<class derived_factory_t, class arg_t, class result_t>
struct memoized_node_factory_singleton : 
	memoized_node_factory<derived_factory_t,arg_t,result_t>
{
	/**\brief Type of the singleton */
	typedef memoized_node_factory_singleton<derived_factory_t,arg_t,result_t>
		node_factory_singleton_t;

	/**\brief Singleton accessor */
	static node_factory_singleton_t& singleton(){
		static node_factory_singleton_t my_singleton;
		return my_singleton;
	}

	/**\brief "The" main singleton memoization method.
		\returns either an existing node (\ref sender_t) corresponding to arg,
		or creates a new by calling the  derived_factory_t::setup_node(arg)		
	*/
	static typename memoized_node_factory_singleton::node_factory_t::sender_ptr
	sget_node(const arg_t& arg){
		return static_cast<derived_factory_t&>(singleton()).get_node(arg);
	}
};




//Singleton Derivation Helper Stuff ######################################
/* NOT needed so far
struct nothing2{};
//forward declaration
template<class body_t,class cond_derived_factory_t=nothing2>
struct memoized_node_factory_singleton_simple;

template<class body_t,class cond_derived_factory_t>
struct cond_derivation{
	typedef memoized_node_factory_singleton_simple<body_t,nothing2> singleton_factory_t;

	typedef typename conditional<
		is_same<cond_derived_factory_t,nothing2>::value,
		singleton_factory_t,cond_derived_factory_t>::type type;
};
*/

/** \brief Singleton Derivation with the most simple concept

	a singleton derivation where there is
	no differenc between the type of the body_t::body_t's argument
	and the memoization arg_t. Classes deriving from this class
	do not need to implement and members.
*/
template<class derived_factory_t, class arg_t, class body_t>
struct memoized_node_factory_singleton_simple 
	: memoized_node_factory_singleton<
							derived_factory_t,
							arg_t, 
							typename body_t::result_t>
{
public:
	/** \brief the concept required member
		As the arguments of body_t's ctor and
		the arguments for the sepcial memoized node factory are the same,
		it simply forwards them to create_node() with a new body.
	*/
	typename memoized_node_factory_singleton_simple::sender_ptr 
	setup_node(const arg_t& arg)
	{
		return create_node(body_t(arg));
	}
};


}

using _flow_graph::memoized_node_factory;
using _flow_graph::memoized_node_factory_singleton;
using _flow_graph::memoized_node_factory_singleton_simple;
using _flow_graph::mtypes;
}


