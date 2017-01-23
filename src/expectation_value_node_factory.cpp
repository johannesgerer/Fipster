#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <tuple>
#include <utility>
#include <cstdlib> //rand
#include <pretty_printer.h>
#endif

#include "expectation_value_body.h"
#include "expectation_value_node_factory.h"
#include "timestepping_operator.h"
#include "boundary_conditions.h"
#include "finite_difference_weights.h"
#include "grid.h"
#include "spacetime_field.h"
#include "task_queue.h"
#include "tbb_flow_graph.h"
#include "global_config.h"
#include "boundary_condition_inhomogenities.h"
#include "Math.h"
#include "PDE.h"

using namespace fipster;
using namespace fipster::expectation_values;


using namespace tbb::flow;

// #################################################
// ##############   PHASE 2:			############
// ##############   Node Factory        ############
// #################################################

void fipster::time_limitation(btime time, btime last, string c)
{
	if(time > last )
		//we have a problem, as the conditional expectation of a
		//function of a past markov state is not function of the current
		//markov state and thus no spacetime_field.
		BOOST_THROW_EXCEPTION
			(runtime_error("Conditional expectation (at "+toS(time)+
										 ") of past values (at "+toS(time)+" not implemented!"));
}



boost::optional<current_step>
config_t::calculate_current_step( btime time, expect_future_t future)
{
	current_step r;
	comb rann(future.field_time,rannacher_step_size);
	auto r_start = future.field_time - rannacher_steps*rannacher_step_size;
	comb* use;

	//find out if we are in rannacher zone:
	if(future.rannacher 
		 && rannacher_steps > 0
		 && time >= r_start)
		{
			if(!rann.contains(time))
				return boost::none;

			r.theta=rannacher_theta;
			r.strang_symmetrization = rannacher_strang_symmetrization;
			use = &rann;
		}else{
		if(!steps.contains(time))
			return boost::none;
		r.theta=theta;
		r.strang_symmetrization = strang_symmetrization;
		use = &steps;
	}

	//never skip the maximum future.timestepping time
	btime next_time = use->next(time,future.field_time);

	// //never skip the minimum time
	// if( start && time < *start )
	// 	next_time = min(next_time,*start);

	//never skip the rannacher starting time
	if(future.rannacher && time < r_start )
		next_time = min(next_time,r_start);

	//if it is the last step, then add the hedging stuff, CodeRefH
	if(future.field_time==next_time)
		r.now=future.details;

	r.future = future.details; //CodeRefBad
																				
	
	if((r.step_size = next_time-time)!=use->step_size){
		cout<<"WARNING: time stepping steps are not aligned, is this intended?\n"
				<<"ID: "<<id<<"\n time: "<<
			time<<", defined: "<<use->step_size<<", derived: "<<r.step_size<<endl;
	}
	
	r.timestep = global_config::get().years_per_btick*r.step_size;
	return r;
}

node_factory::bv_sender_ptr
node_factory::boundary_setup( const arg_t& arg )
{
	auto& time = arg.time;
	auto& future = arg.carg;
	if(time==future.field_time)
		return future.field->bv_get_node(arg.strip());

	calculate_current_step(arg);

	auto node = create_node(boundary_values_body_t());
	auto joinnode = &node->joinnode;

	//use task queue to break recursion:
	task_queue::enqueue([=](){
		make_edge(*get_node(arg),input_port<0>(*joinnode)); //interior values
		make_edge(*boundary_conditions::sget_node(make_pair
  							(*betterAt(this->objs->bcs_configs,config->boundary_conditions_s
													 ,config->id)
								 ,arg.grid))
							,input_port<1>(*joinnode));//boundary conditions
		make_edge(*bc_inhomogentities::sget_node
							(make_pair(arg,
												 static_pointer_cast<node_factory>(shared_from_this()))),
							input_port<2>(*joinnode)); //bc inhomogeneities
	});
		
	return node;
}

/** rannacher stepping convention: 
		1. rannacher_steps:
		determines the number of small (rannacher) time steps before 
		the maximum stopping time.
		2. stepSize:
		determines the step size of the small (rannacher) steps.
		3. Rannacher steps are always performed in requested number and size,
		except, if the minimum stopping time, splits on rannacher step in two.
		(Due to a mismatch in step size, offset, maximum stopping time 
		and rannacher step size, its possible that the last non-rannacher 
		step is smaller than the step_size, to ensure 3.)

		\returns boost::none is the step is not supported
*/

array<double,2> current_step::betas(){
	return { -theta*timestep, (1-theta)*timestep};
}


current_step node_factory::calculate_current_step(const arg_t& arg){
	auto& time = arg.time;
	auto& future = arg.carg;

  time_limitation(time,future.field_time,config->id); 

	//calculate the parameters of this time step
	return *time_supported(config->calculate_current_step(time,future)
												 ,time
												 ,config->id);
}

/**
	 creates the nodes that produce the futures required to perform
	 a expectation value time-stepping step.
	 Due the recursive nature of time-stepping for expectation values 
	 this method dispatches tasks using the singleton task_queue
	 object, to prevent stack overflow.
*/
node_factory::sender_ptr node_factory::inner_setup( const arg_t& arg )
{
	// cout<<"EX "<< arg.time << ": " <<arg.carg <<endl;
	auto& time = arg.time;
	auto& future = arg.carg;

	if(time==future.field_time)
		//do not not exponentiate, raise to raise to a power or add current
		//hedging position's value CodeRefH
		return future.field->get_node(arg.strip());

	auto c = calculate_current_step(arg);

	//Sender: finite difference weights
	auto fdw_conf = betterAt(objs->fdweights_configs
													 ,config->fd_weights_s
													 ,config->id);

	auto pde=dynamic_pointer_cast
		<const finite_difference_weights::BSiso>(fdw_conf->PDE);

	double vola=0;
	if(pde) vola=pde->volatility;

	body_t body=body_t(config
										 ,c
										 ,arg
										 ,"id: "+config->id+", time: "+toS(time)+
										 ", step_size: "+toS(c.step_size)+
										 ", strang_symmetrization: "+toS(c.strang_symmetrization)+
										 ", theta: "+boost::lexical_cast<string>(c.theta),
										 vola,time);

	auto fd_weights = fdweights::sget_node
		(make_pair(*fdw_conf,arg.grid));

	//Sender: boundary conditions
	auto bc = boundary_conditions::sget_node
		(make_pair(*betterAt(objs->bcs_configs
												 ,config->boundary_conditions_s
												 ,config->id)
							 ,arg.grid));

	//create a expectation value body using the create_node function for
	// unary bodies, because future have to be registered separately to break the recursion
	// induced by the dependency of underlying_field and future_expectation
	// (see below)
	auto node = create_node(body);
	auto joinnode = &node->joinnode;

	make_edge(*timestepping_operator::sget_node
						(make_tuple
						 (c.betas(),c.strang_symmetrization
							,make_pair(make_tuple(bc.get(),fd_weights.get())
												 ,fdw_conf.get())))
						,input_port<2>(*joinnode));
		
	//icc bug workaround: could not use "this" directly in lambda expression!
	auto this2 = static_pointer_cast<node_factory>(shared_from_this());
	auto lcp = future.exercise;

	//use task queue to break recursion:
	task_queue::enqueue([=](){
			
			if(lcp)
				make_edge(*lcp->get_node(arg.strip()),input_port<0>(*joinnode));
			else
				make_edge(*default_node<discretized_space_field<1> >()
									,input_port<0>(*joinnode));

			auto future_arg=arg.time_shifted(c.step_size);
			//Sender: future expectation value
			auto future_expectation = this2->get_node(future_arg);

			make_edge(*future_expectation,input_port<1>(*joinnode));

			//Sender: boundary inhomogeneities
			auto bc_inhom = bc_inhomogentities::sget_node(make_pair(arg,this2));
			//Sender: future boundary inhomogeneities
			auto bc_inhom2 = bc_inhomogentities::sget_node(make_pair(future_arg,this2));
			/*\remark:	the case, where the future time, is the latest stopping 
				time, (i.e. maturity: the future value is the financial terminal value)
				is not handled specially. 
				For more info see: boundary_value_delegation
			*/

			//theoretically bc_inhom[i] is only needed if the corresponding
			//beta is not zero, but as the node requires 5 futures and I
			//have no mechanism to disable a future, i leave it, because
			//this is no overhead, as at some point the inhomogenity is
			//needed anyhow.
			make_edge(*bc_inhom,input_port<3>(*joinnode));
			make_edge(*bc_inhom2,input_port<4>(*joinnode));
		});

	return node;
}

node_factory::node_factory( config_ptr config,shared_const_objs_t objs )
	: config(config),objs(objs)
{}
		
	
/** old: check bc inhom setup and this boundary_setup for current behavior

		overwrite boundary_value_delegation:

		This is due to this convention:
		The boundary values of the expectation at the latest stopping time,
		which are known, as they coincide with the boundary values of the 
		underlying field process, are not used, but instead also subject to
		the expectation's boundary conditions.

		The alternative would be:
		Use the underlying's boundary values. This however implies either
		of the following changes:
		1) The elimination of boundary values and subsequent correction of the
		time stepping matrix and inhomogenities have to to be deactivated and 
		instead the correct inhomogenities araising from known boundary values
		have to be used.
		2) The boundary conditions have to be such, that the resulting boundary values
		coincide with the boundary values of the underlying. This can be achieved
		by using "transparent" boundary values:
		inhomValue == type, inhomSource == exerciseValue and timeSteps=0
		(This would require some time dependece in the specifications 
		of boundary conditions.)
*/
// node_factory::bv_sender_ptr
// node_factory::boundary_delegation(arg_t& arg){
// 	// INFO: calling this virtual function generally works (as testet
// 	// by explicit_filds), but it seems, that the boundary values of
// 	// this example (commit 31139ed2a) are not called at all due to
// 	// inhomSource=exerciseValue and so this function is never called.

// 	//perform same argument transformations as the (interior spacetime_field)
// 	field_delegation(arg);
// 	//but never delegate to underlying:
// 	return 0;
// }


void config_t::check_arbitrage(shared_ptr<option_config_t> o_conf
															 ,shared_const_objs_t objs
															 ,btime fut_time
															 ,btime arg_time){
	if(o_conf->entropic_theta || !o_conf->hedging())
		return;
		
	double timestep  = global_config::get().years_per_btick
		*(fut_time-arg_time);

	auto&   k1 =o_conf->hedging()->k1;
	double  ks =(1+k1)/(1-k1);
	auto   pde=dynamic_pointer_cast
		<const finite_difference_weights::BSiso>(betterAt(objs->fdweights_configs
																											,fd_weights_s
																											,id)->PDE);
	FIPSTER_ASSERT(pde);
		
  double s=pde->volatility;
	double e=exp(((pde->discounting_rate + o_conf->discounting_rate)
								- pde->drift)*timestep);
			//check for arbitrage opportunity
	double A = o_conf->sigma_factor * sqrt(exp(s*s*timestep)-1);
	double B = 1-e*ks;
	double C = -1+e/ks;
	bool arb = A < max(B,C);
	if(arb){
		cout<<"\nsigma_factor                         "<< o_conf->sigma_factor
				<<"\nr=sqrt(exp(vol*vol*timestep)-1)     "<< sqrt(exp(s*s*timestep)-1)
				<<"\nA=sigma_factor*r                       "
				<< o_conf->sigma_factor*sqrt(exp(s*s*timestep)-1)
				<<"\ne=exp((discounting - drift)*timestep) "<<e
				<<"\nks=(1+k1)/(1-k1)                    "<<ks
				<<"\nB=1-e*ks                              "<<1-e*ks
				<<"\nC=-1+e/ks                             "<<-1+e/ks
				<<endl;
		if(A<B)
			cout << "A<B => LONG ARBITRAGE\n" << endl;
		if(A<C)
			cout << "A<C => SHORT ARBITRAGE" << endl;
		cout<< "option id: "<<o_conf->id << "\nepectation id: "
				<<o_conf->expectation_s
				<<"\n time: "<<arg_time<<"  future: "<<fut_time
				<<"  step: "<<timestep<<"\n"<<endl;
	}
	// FIPSTER_ASSERT(!arb);
	
}
