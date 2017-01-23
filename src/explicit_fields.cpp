#include "precompiled.h"
#include "explicit_fields.h"
#include "auto_connecting_body.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <string>
#include <boost/optional.hpp>
#include <pretty_printer.h>
#endif

#include "logging.h"
#include "grid_fields.h"
#include "black_scholes.h"
#include "global_config.h"
#include "Math.h"
#include "finite_difference_weights_config.h"

namespace fipster { namespace explicit_fields {

	using namespace std; 

	// #################################################
	// ##############   PHASE 1:			############
	// ##############   CONFIGURATION       ############
	// #################################################

	node_factory::node_factory(const ptree& pt,shared_const_objs_t objs)
		: pt(pt),objs(objs)
	{
		string stype(pt.get<string>("type"));
		if(stype == "CallPayoff") type = CallPayoff;
		else if(stype == "PutPayoff") type = PutPayoff;
		else if(stype == "BinaryCallPayoff") type = BinaryCallPayoff;
		else if(stype == "Put") type = Put;
		else if(stype == "Call") type = Call;
		else BOOST_THROW_EXCEPTION(runtime_error(stype+" is no valid Explicit Field Type"));
	}

	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################
	
	/** this class calculates the payoff of a Call/Put on the first
			assert
	 */
	template<type_t type>
	struct VanillaPayoff_t{
		double strike;
		
		/** constructor that reads in configuration */
		VanillaPayoff_t
		(const ptree& pt,const field_arg_t<>&,shared_const_objs_t)
			:strike(pt.get<double>("strike")){}

		/** core function, computing the value for a given state */
		double getValue(const state_t& state){
			/*if(state[0]>120){
				auto s2=state;
				s2[0]-=120;
				return getValue(s2);
			}
			if(state[0]>75) return 0;*/
			double a=30;
			if(type==BinaryCall) return state[0] >= strike ? 1 : 0;
			
			return max( (type==Call ? +1 : -1 )*(state[0]-strike-0*a),0.0);
// -
// 			max( (type==Call ? +1 : -1 )*(state[0]-strike-a),0.0)-
// 				max( (type==Call ? +1 : -1 )*(state[0]-strike-2*a),0.0)+
// 				max( (type==Call ? +1 : -1 )*(state[0]-strike-3*a),0.0);
		};
	};
	
	/** this class calculates the Black Scholes Values for a Call/Put on 
	the geometric average of grid.D assets (using equal correlation, dividend and volatility)
	*/
	template<type_t type>
	struct BS_t{
		double strike,interest,volatility,correlation;
		black_scholes bs;

		/** constructor that reads in configuration */
		BS_t(const ptree& pt,const field_arg_t<>& arg,shared_const_objs_t objs){
		
			auto same_as = pt.get_optional<string>("sameAs");

			if(!same_as)
				bs=black_scholes(0,
					pt.get<double>("strike"),
					global_config::get().years_per_btick*(pt.get<btime>("expiration")-arg.time),
					pt.get<double>("interest"),
					pt.get<double>("volatility"),
					pt.get("dividend",0.0),
					pt.get("correlation",0.0),
					arg.grid->D);
			else{
				FIPSTER_THROW_EXCEPTION(runtime_error("not implemented"));
				/*
				//get the expectation value to take the interval from (if available)
				auto nf=dynamic_cast<const expectation_values::node_factory*>
					(objs->spacetime_fields.at(*same_as).get());
				if(!nf) FIPSTER_THROW_EXCEPTION(runtime_error("sameAs: "+*same_as+" is no expcetation_value"));
				auto fdw=objs->fdweights_configs.at(nf->config->fd_weights_s);
				fdw->
				*/
			}
		}

		/** core function, computing the value for a given state */
		double getValue(const state_t& state);
	};

	template<>
	double BS_t<Call>::getValue(const state_t& state)
	{
		return bs.call(geometric_average(state,bs.n));
	};

	template<>
	double BS_t<Put>::getValue(const state_t& state)
	{
		return bs.put(geometric_average(state,bs.n));
	};

	//#######################  Node body ##############################
	template<int sg, class result_t, class state_it,class core_t>//specify the sub-grid to work on
	struct explicit_field_body : 
		auto_connecting_body<explicit_field_body<sg,result_t,state_it,core_t>,//body_t
							result_t>//result
	{
		core_t core;
		const btime time;
		grid_ptr grid;
		string meta_info;
		double factor;

		using result_sptr_t = typename explicit_field_body::result_sptr_t;
		using input_t = typename explicit_field_body::input_t;

		explicit_field_body(field_arg_t<> field_args,const ptree& pt,shared_const_objs_t objs)
			:core(core_t(pt,field_args,objs))
			,time(field_args.time)
			,grid(field_args.grid)
			,meta_info(pt.get<string>("type")
								 + " EXPL FL ("+toS(time)+") "+pt.get<string>("<xmlattr>.id"))
			,factor(pt.get("factor",1.0))
		{};


		//computation function #####################################################
		result_sptr_t operator()(const input_t&){
			// logging
			if(1)thread_logger()<<"EXPL FL ("<<time<<") | START"<<endl;
			auto result = make_shared<result_t>
				(meta_info);

			//allocate new space for the result on the specified Grid
			auto end = result->resize(*grid);

			for(state_it it(*grid); it!=end; ++it)
				result->at(it.index()) = factor*core.getValue(it.state());
				
	//tbb::this_tbb_thread::sleep(tbb::tick_count::interval_t(0.2));

			if(1)thread_logger()<<"EXPL FL| DONE"<<endl;
			return result;
		};

	};

	// #################################################
	// ##############   PHASE 2:			############
	// ##############   Node Factory        ############
	// #################################################

	boost::optional<node_factory::arg_t>
	node_factory::delegation(const arg_t& arg){
		if((type==CallPayoff || type==PutPayoff || type==BinaryCallPayoff) && arg.time != 0 ){
			arg_t a=arg; a.time = 0;
			return a;
		}
		return boost::none;
	}

	template<int sg,class result2_t,class iterator2_t,class sender2_t>
	sender2_t node_factory::selecting_setup( const node_factory::arg_t& arg)
	{
		switch(type){
		case CallPayoff:
			return create_node
				(explicit_field_body<sg,result2_t,iterator2_t,VanillaPayoff_t<Call> >(arg,pt,objs));
		case PutPayoff:
			return create_node
				(explicit_field_body<sg,result2_t,iterator2_t,VanillaPayoff_t<Put> >(arg,pt,objs));
		case BinaryCallPayoff:
			return create_node
				(explicit_field_body<sg,result2_t,iterator2_t,VanillaPayoff_t<BinaryCall> >(arg,pt,objs));
		case Call:
			return create_node
				(explicit_field_body<sg,result2_t,iterator2_t,BS_t<Call> >          (arg,pt,objs));
		case Put:
			return create_node
				(explicit_field_body<sg,result2_t,iterator2_t,BS_t<Put> >           (arg,pt,objs));
		default:
			BOOST_THROW_EXCEPTION(runtime_error(toS(type)+" is no valid Explicit Field Type"));
		}
	}
		
	// spacetime_field: create_node
	node_factory::sender_ptr node_factory::inner_setup(const arg_t& arg){
		auto a = delegation(arg); if(a) return get_node(*a);
		return selecting_setup<J_SGS,
								discretized_space_field<J_SGS>,
								grid_iterator<J_SGS,true>,
								sender_ptr
							>(arg);
	}

	node_factory::bv_sender_ptr 
		node_factory::boundary_setup( const arg_t& arg ){ 
		auto a = delegation(arg); if(a) return bv_get_node(*a);
		return selecting_setup<J_SGS,
								boundary_field_t<J_SGS>,
								boundary_iterator<J_SGS-1,false,true>,
								bv_sender_ptr
							>(arg);
	}

}}
