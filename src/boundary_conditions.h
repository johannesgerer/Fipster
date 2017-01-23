#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <memory>
#endif

#include "grid_fields.h"
#include "auto_connecting_body.h"
#include "memoized_node_factory.h"
#include "boundary_condition_iteration.h"

namespace fipster { namespace _boundary_conditions {

	using namespace std;

	struct sg_t{
		const static int sg = 0;//specify the sub-grid to work on
	};

	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################

	


	typedef sparse_field<boundary_index<J_SGS-1>,grid_index<J_SGS>,RowMajor> Bi_t;
	typedef sparse_field<boundary_index<J_SGS-1>,boundary_index<J_SGS-1>,ColMajor> Bo_t;

	struct B_t {
		Bi_t Bi;
		Bo_t Bo1;
		Bo_t Bo;
	};


	

	struct body_t : auto_connecting_body<body_t,B_t> {
		const static int sg = J_SGS-1;
		const bc_arg_t arg;
		body_t(const bc_arg_t& arg);
		result_sptr_t operator()(const input_t&);
	};

	

	// #################################################
	// ##############   PHASE 2:			############
	// ##############   Node Factory        ############
	// #################################################
	struct boundary_conditions 
		: memoized_node_factory_singleton_simple<
				boundary_conditions,
				bc_arg_t,
				body_t>
	{};

}

using _boundary_conditions::B_t;
typedef shared_ptr<const B_t> B_ptr_t;

using _boundary_conditions::boundary_conditions;

}
