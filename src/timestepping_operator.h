#pragma once
#ifndef AUTOMATIC_PRECOMPILATION
#include <vector>
#include <utility>
#include <map>
#include <array>
#include <memory>
#include <vector>
#endif


#include "memoized_node_factory.h"
#include "auto_connecting_body.h"
#include "tridiagonally_split_grid_operator.h"
#include "boundary_conditions.h"
#include "finite_difference_weights.h"

namespace fipster { 
	namespace _timestepping_operator {

	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################

	/** time stepping operators: A_o, B and A_i split for every direction
		if strang_symmetrization is set, all direction except the last one 
		will have their "beta"-values split in half.
		if applicable there will be two version with "beta" and "-beta+1", which
		is needed for certain implicit timestepping schemes.
	*/

	struct timestepping_op_t {
		
		array<tridiagonally_split_grid_operator,2> op;
		shared_ptr<const B_t>	B;
		array<double,2> betas;
		
		timestepping_op_t (	shared_ptr<const B_t> B, 
							tridiagonally_split_grid_operator::gsic_ptr gsi,
							array<double,2> betas)
		:B(B),betas(betas)
		{
			for(uint i=0;i<2;i++)
				if(betas[i]!=0.0)
					op[i].init(gsi);
		}
	};

	// #################################################
	/**\brief Operator splitter (node body)

		Takes no arguments, but two futures:
			- Future 1: Operator A in raw form (fd_weights_result_t)
			- Future 2: (inverted) Boundary operators (B_t)

		And splits the raw A and combines it with B in a timestepping_op_t object
		which is everything that is needed for time stepping.
	*/
	struct splitter_body_t : 
		auto_connecting_body<	splitter_body_t,	//body
								tridiagonally_split_grid_operator,	//result
								B_t,				//future 0
								fd_weights_result_t	//future 1
							>
	{
		const finite_difference_weights::config_t& fdw_conf;
		
		splitter_body_t(const finite_difference_weights::config_t& fdw_conf)
			: fdw_conf(fdw_conf){};
		
		result_sptr_t operator()(const input_t& futures);
		
	};

	// #################################################
	

	/**\brief Body to perform Result = 1 + beta*op	*/
	struct multiply_adder_body_t 
		: auto_connecting_body<	multiply_adder_body_t,//body
								timestepping_op_t,//result
								//future 0
								tridiagonally_split_grid_operator,	//splitting result
								//future 1
								B_t,				
								//future 2
								fd_weights_result_t>
	{
		array<double,2> betas;
		bool strang_symmetrization;
		//Ctor
		multiply_adder_body_t(array<double,2> betas, bool strang_symmetrization);
		result_sptr_t operator()(const input_t& futures);
	};

	// #################################################
	// ##############   PHASE 2:			############
	// ##############   Node Factory        ############
	// #################################################

	

	/** \brief The argument-less splitter_body is memoized on its 
				futures' senders (BC and fd_weights)
	*/
		typedef pair<splitter_body_t::sender_t
								 ,const finite_difference_weights::config_t*>
		splitter_arg_t;
		
	struct splitter 
		: memoized_node_factory_singleton<splitter, //node_factory_t
		splitter_arg_t,	//arg
		tridiagonally_split_grid_operator> //result_t
	{
		sender_ptr setup_node(const arg_t& arg);
	};
	
	/** \brief argument type for the multiply_adder memoized node factory

			 - get<0>(tuple) := ("beta_A","beta_B") (might be the trivial 0)
			 - get<1>(tuple) := use strang symmetrization
			 - get<2>(tuple) := future sender expected by the body (splitter)
			 
			 uses sender_t and not shared_ptr (passed internally to
			 make_edge), as it does not store these objects and thus storage
			 must behandled somewhere else.
	*/
	typedef tuple<array<double,2>,bool,splitter_arg_t> 
		timestepping_operator_body_arg_t;

	
	struct timestepping_operator 
		: memoized_node_factory_singleton<timestepping_operator,
		timestepping_operator_body_arg_t,timestepping_op_t>
	{
		sender_ptr setup_node(const arg_t& arg);
	};
}
using _timestepping_operator::timestepping_operator;
using _timestepping_operator::timestepping_op_t;
}
