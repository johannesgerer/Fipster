#pragma once
#include "spacetime_field.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/property_tree/ptree.hpp>
#endif


#include "configured_objects.h"


namespace fipster { namespace explicit_fields {

	using boost::property_tree::ptree;

	// #################################################
	// ##############   PHASE 1:			############
	// ##############   CONFIGURATION       ############
	// #################################################

		enum type_t {CallPayoff,Call,PutPayoff,Put,BinaryCallPayoff,BinaryCall};

	/*configuration is very heterogeneous among the different
	  explicit fields and is referred to the calculating class
	*/

	// #################################################
	// ##############   PHASE 2:			############
	// ##############   Node Factory        ############
	// #################################################


class node_factory : public spacetime_field<J_SGS>
{
public:
	const ptree pt;
	type_t type;
	shared_const_objs_t objs;

	sender_ptr inner_setup(const arg_t& arg);
	bv_sender_ptr boundary_setup(const arg_t& arg);

	template<int sg,class result2_t,class iterator2_t,class sender2_t>
	sender2_t selecting_setup( const node_factory::arg_t& arg);

	//Constructor #########################################################
	node_factory(const ptree& pt,shared_const_objs_t objs);
	// spacetime_field_ptr boundary_value_delegation(arg_t& arg){
	// 	cout<<"SHIT"<<endl;
	// 	return 0;
	// }

	//Overwriting to allow for constant payoffs
	boost::optional<arg_t> delegation(const arg_t& arg);
};



}}
