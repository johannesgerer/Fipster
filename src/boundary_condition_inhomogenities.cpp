#include "precompiled.h"
#include "boundary_condition_inhomogenities.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <vector>
#include <tuple>
#endif

#include "spacetime_field.h"
#include "configured_objects.h"
#include "utility.h"
#include "runtime_join_node.h"
#include "boundary_conditions.h"

namespace fipster { namespace _boundary_conditions {
	
	using namespace std;

	enum {boundaries_selector, interiors_selector};

	// #################################################
	// ##############   PHASE 2:			      ############
	// ##############   Node Factory        ############
	// #################################################
	inhom_body_t::inhom_body_t(const bc_arg_t& arg,
		const bc_inhomogentities_arg_t& bcinhom_arg)
		:arg(arg),bcinhom_arg(bcinhom_arg){};

	bc_inhomogentities::sender_ptr 
		bc_inhomogentities::setup_node( const arg_t& arg )
	{
		//alias
		auto& expectation_factory = arg.second;
		auto& spacefield_args=arg.first;
		auto& expectation_config = *expectation_factory->config;
		auto& objs=*expectation_factory->objs;
		auto& BCs=*betterAt(objs.bcs_configs,expectation_config.boundary_conditions_s
												,expectation_config.id);

		vector<spacetime_field<J_SGS>::sender_ptr> required_spacefields;
		vector<spacetime_field<J_SGS>::bv_sender_ptr> required_boundaryvalues;

		inhom_body_t body(make_pair(BCs,spacefield_args.grid),arg);

    // gcc: i think the following means:

    //it is possible to use a boundary condition object, that
    // describes more state variables than the the spacefield

    // \todo how is this handled in other cases (i.e. the homogenious
    // BCs,...) and why is it the spacefield whose dimension is
    // important. the receiver that uses the BCs should dictate what
    // dimensions it requires and these cannot surpass the BCs's
    // specified dimensions and the available spacefields' dimensions.
   
		FIPSTER_ASSERT(BCs.sv_vec.size() >= spacefield_args.grid->D);
    uint size=spacefield_args.grid->D;  //gcc: was  min(D,size())
		body.sources.resize(size);

	//1. create nodes for fields that are involved in the calculation 

		//loop over all boundary planes:
		auto source=begin(body.sources);
		for(uint it = 0;it<size;it++,++source){
			//loop over upper and lower plane
			for(int i=0;i<2;i++){

				auto& bc_conf=BCs.sv_vec[it].lowerupper[i];
		
				if(bc_conf.inhom_value==zero_value){
					(*source)[i][boundaries_selector]=-1;
					(*source)[i][interiors_selector]=-1;
				}else{
					//get pointer to source of the boundary value computation
					shared_ptr<spacetime_field_base<J_SGS> > spacetime_field;
					switch(bc_conf.inhom_source){
						case exercise_value:
							{ auto e=spacefield_args.carg.exercise;
							if(!e) BOOST_THROW_EXCEPTION(runtime_error
							 ("no exercise value to be used as BC inhoms' source. (nr "+toS(it)+")"));
							spacetime_field = e; }
							break;
						case future_value:
							spacetime_field = 
								expectation_factory->create_proxy(spacefield_args.carg);
							break;
					  default:
							FIPSTER_THROW_EXCEPTION(runtime_error("bad bc_conf.inhom_source"));
					}

					//calculate timeshift as specified in boundary condition config file
					//CodeRefB
					//create a new argument with shifted time
					auto arg2(spacefield_args.strip());
					if(bc_conf.time_steps==0)
						arg2.time=spacefield_args.carg.field_time;
					else{
						assert(bc_conf.time_steps>0);
						arg2.time+=min(spacefield_args.carg.field_time
													 ,expectation_config.steps.step_size*bc_conf.time_steps);
					}

					//get nodes for interior and boundary and register them in 
					//the body's 'sources' vector
					(*source)[i][interiors_selector]=get_index(required_spacefields,
									spacetime_field->get_node(arg2));
					(*source)[i][boundaries_selector]=get_index(required_boundaryvalues,
									spacetime_field->bv_get_node(arg2));
				}
			}
		}
		
	//2. create node
		return create_node(body,make_tuple(
			new_runtime_join_node(required_spacefields).get(),
			new_runtime_join_node(required_boundaryvalues).get()
			));
	}


	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################

	void inhom_body_t::add_to_interior(	boundary_index<sg>& row,grid_index<sg+1>& col,
			double v,uint variable,bstate_t state){
		result->at(row) += v * (*interiors)[sources[variable][state][interiors_selector]]->at(col);
	};

	void inhom_body_t::add_to_boundary(	boundary_index<sg>& row,boundary_index<sg>& col,
			double v,uint variable,bstate_t state){
		result->at(row) += v * (*boundarys)[sources[variable][state][boundaries_selector]]->at(col);
	};
	
	void inhom_body_t::set_zero(boundary_index<sg>& row,uint variable,bstate_t state){
		result->at(row) =0;
	};


	/** calculates BCinhom as specified in the configurations
	*/
	inhom_body_t::result_sptr_t inhom_body_t::operator()( const input_t& futures )
	{
		//logging
		if(0)thread_logger()<<"BCInhom("<<bcinhom_arg.first.time<<","<<bcinhom_arg.second->config->id<<")| START"<<endl;

		interiors=get<0>(futures).get();
		boundarys=get<1>(futures).get();

		result=make_shared<boundary_field_t<J_SGS>>();
		result->resize(*arg.second);
		
		(*result)->setZero();
		 
		for(boundary_condition_iterator<body_t::sg,&bc::inhom_value> bcit(arg);
				bcit.valid(); ++bcit)
			bcit.handle_bc_line(*this);

		if(0)thread_logger()<<"BCInhom| DONE"<<endl;
		
		return move(result);
	}

}}
