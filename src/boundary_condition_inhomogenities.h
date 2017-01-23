#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <vector>
#endif


#include "auto_connecting_body.h"
#include "memoized_node_factory.h"
#include "boundary_conditions.h"
#include "spacetime_field.h"
#include "expectation_value_node_factory.h"

namespace fipster { namespace _boundary_conditions {

		using namespace Eigen;
		using expectation_values::bc_inhomogentities_arg_t;

	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################
	
	// #################################################
	/**  */
	struct inhom_body_t : 
		auto_connecting_body<	inhom_body_t	//body
								,boundary_field_t<J_SGS>	//result
								//future 0 (vector of interior values):
								,vector<spacetime_field<J_SGS>::result_sptr_t>
								//future 1 (vector of boundary values):
								,vector<spacetime_field<J_SGS>::bv_result_sptr_t>
							>
		, boundary_conditions_line<body_t::sg>
	{
		const bc_arg_t arg;
		const bc_inhomogentities_arg_t bcinhom_arg;
		
		inhom_body_t(const bc_arg_t& arg,const bc_inhomogentities_arg_t& bcinhom_arg);

		typedef vector<array<array<int,2>,2> > sources_t ;
		sources_t sources;		

		result_sptr_t operator()(const input_t& futures);	

		void add_to_interior(	boundary_index<sg>& row,grid_index<sg+1>& col,
			double v,uint variable,bstate_t state);

		void add_to_boundary(	boundary_index<sg>& row,boundary_index<sg>& col,
			double v,uint variable,bstate_t state);
		
		void set_zero(	boundary_index<sg>& row,uint variable,bstate_t state);

		const vector<spacetime_field<J_SGS>::result_sptr_t>* interiors;
		const vector<spacetime_field<J_SGS>::bv_result_sptr_t>* boundarys;
		editable_result result;
	};


	// #################################################
	// ##############   PHASE 2:			############
	// ##############   Node Factory        ############
	// #################################################



	/** this represents the calculation of the so called Boundary Condition Inhomogenities.
		They constitute the inhomogenious part of the boundary conditions of a expectation_value, 
		thus they are calculated FOR an expectation_value and not FROM it.
		(as opposed to the boundary values of an expectation_value).

		\code
		B_i.x_i + B_o.x_o = BCInhom
		\endcode
		They can either be zero or be defined using other spacetime_fields:
		\code
		BCInhom := B_i.y_o + B_o.y_o
		\endcode

		where y is a combination of future and underlying value.

	*/
	struct bc_inhomogentities
		: memoized_node_factory_singleton<
			bc_inhomogentities,			//body
			bc_inhomogentities_arg_t,	//arg
			boundary_field_t<J_SGS>		//result
		>
	{
		sender_ptr setup_node(const arg_t& arg);
	};
}

using _boundary_conditions::bc_inhomogentities;

}
