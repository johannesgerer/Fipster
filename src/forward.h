#pragma once

namespace fipster { 

	//forward declarations:
	struct nothing;
	namespace _grid { class Grid; }
	namespace _fields { template<int n,class carg_t=nothing>struct spacetime_field;
		template<int n,class carg_t>struct spacetime_field_proxy;
		template<int n,class carg_t=nothing>struct spacetime_field_base;
	  template<class carg_t=nothing> struct field_arg_t; }
	namespace expectation_values { struct node_factory; }
	struct option_node_factory; 
	namespace _boundary_conditions { struct config_t; }
	namespace finite_difference_weights { struct config_t; }
	namespace _solution { struct solution; }


	typedef shared_ptr<_fields::spacetime_field_base<J_SGS> >
	plain_spacetime_field;
}
