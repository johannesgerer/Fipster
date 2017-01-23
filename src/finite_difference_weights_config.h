#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <utility>
#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <comphash.h>
#endif

#include "configured_objects.h"

namespace fipster { namespace finite_difference_weights {

	using namespace std;
	using namespace boost::property_tree;

	// #################################################
	// ##############   PHASE 1:			############
	// ##############   CONFIGURATION       ############
	// #################################################
	enum stenciltype_t { StencilFull , StencilNoDiagonals, StencilOneDiagonal, StencilTwoDiagonals };
	enum PDE_type_t { BlackScholesIsometric };

	/** This base class only serves one purpose: Make the derived classes
	comparable in a (type-)safe way. The purely virtual method comparison_string has to
	return a string that uniquely identifies or describes the objects.
	The actual "double operator()(...)" members are accessed using static polymorphism at the
	body_t level in:

	fdweights::sender_ptr fdweights::setup_node( const arg_t& arg )
	*/
	struct PDE_base : comphashable<PDE_base>,stringifyable<PDE_base> {		
		string type;
		//ctor
		PDE_base(string type):type(type){};
		template<class T>
		void my_combine(T& t)const
		{
			t(&PDE_base::type,&PDE_base::comparison_string);
		}
		virtual string comparison_string() const=0;
	};

	//Not necessarily const:
	typedef shared_ptr<const PDE_base> PDE_base_ptr;

	struct config_t : comphashable<config_t>,stringifyable<config_t>{

		//fields
		string id;
		PDE_base_ptr PDE;
		stenciltype_t stencil_type;
		int partial_derivative_order;
		double higher_order_weight;
		
		//ctor
		config_t( const ptree& pt,objs_t objs );

		template<class T> void my_combine(T& t)const{
			t(	byValue(&config_t::PDE),
				&config_t::stencil_type,
				&config_t::partial_derivative_order,
				&config_t::higher_order_weight);
		}
	};

}
}
