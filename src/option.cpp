#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/utility/in_place_factory.hpp>
#endif

#include "option.h"
#include "runtime_join_node.h"
#include "global_config.h"
#include "auto_connecting_body.h"
#include "expectation_value_node_factory.h"
#include "finite_difference_weights.h"
#include "PDE.h"


using namespace fipster;
using namespace options;

// #################################################
// ##############   PHASE 1:	      		############
// ##############   Config              ############
// #################################################
hedging_t::hedging_t(const ptree& pt,string c)
	// : combA(pt, c)
	: id(c+".hedging")
	, k0(betterGet<double>(pt,"transactionCosts.<xmlattr>.k0",id))
	, k1(betterGet<double>(pt,"transactionCosts.<xmlattr>.k1",id))
	, infinite(pt.get("transactionCosts.<xmlattr>.infinite",false))
	, initialized(false)
	, grid_s(betterGet<string>(pt,"grid",id))
{}
		
delta_hedging_t::delta_hedging_t(const ptree& pt, const string& context)
	:id(context+".delta_hedging")
	,pt(pt){};

void hedging_t::init(shared_const_objs_t objs){
	if(!initialized){
		grid = betterAt(objs->grids,grid_s,id);
		thread_logger()<<"Hedging grid ("<<grid_s<<"):\n"<<*grid<<endl;
		initialized=true;
	}
}

referenced_field_t delta_hedging_t::source(shared_const_objs_t objs){
	if(!my_source){
		vector<referenced_field_t> fields;
		read_referenced_fields(id,nullptr,pt,fields,objs);
		if(fields.size() == 1)
			my_source = fields.front();
		else
			FIPSTER_ASSERT(fields.empty());
		// \CodeRef15-11-05 not having a delta source should be a way to specify delta=0 zero?
	}
	return *my_source;
}
		
const boost::optional<hedging_t>& option_config_t::hedging() const{
	FIPSTER_ASSERT(!my_hedging || my_hedging->initialized);
	return my_hedging;
}

const boost::optional<hedging_t>&
option_config_t::hedging(shared_const_objs_t objs){
	if(my_hedging)
		my_hedging->init(objs);
	return my_hedging;
}

exercise_t::exercise_t(const ptree& pt,string c)
	// : combA(pt, c)
	: field_s(betterGet<string>(pt,"value",c  +".exercise"))
	, implicit(betterGet<bool>(pt,"implicit",c+".exercise"))
{
	if(implicit) FIPSTER_ASSERT(!"implemented"); //\todoI
}
		
option_config_t::option_config_t( const ptree& pt )
	: id(betterGet<string>(pt,"<xmlattr>.id","field: option"))
	, expectation_s(betterGet<string>(pt,"expectation",id))
	, exercise(betterChild(pt,"exercise",id),id)
	, times(betterChild(pt,"times",id),id)
	, entropic_theta(pt.get_optional<double>("entropicTheta"))
	, sigma_factor(pt.get<double>("sigmaFactor",0.0))
	, powers(sigma_factor==0?1:2)
	, discounting_rate(betterGet<double>(pt,"discountingRate",id))
{
	if(pt.get<uint>("iterations",1) != 1 && sigma_factor !=0)
		FIPSTER_THROW_EXCEPTION
			(runtime_error("iterations !=0 only makes sense for non-time "
										 "consistent pricing functions.\nContext: "+toS(id)));
		
	auto hed=pt.get_child_optional("hedging");
	if(hed) my_hedging = boost::in_place(*hed,id);
  if(sigma_factor != 0 && entropic_theta)
		FIPSTER_THROW_EXCEPTION
			(runtime_error("sigma_factor != 0 and entropic_theta cannot be combined\nContext: "+toS(id)));
}

// #################################################
// ##############   PHASE 2:		      	############
// ##############   Node Factory        ############
// #################################################


option_carg_t option_carg_t::to_marketarb() const{
	auto carg = *this;
	carg.exercise = never;
	if(carg.hedging){
		carg.hedging->delta = boost::none; //nullptr;//zero hedging CodeRef15-11-05
		carg.hedging->skip = 1;
	}
	return carg;
}

option_carg_t option_carg_t::change_hedging_index(grid_iterator<0,true> it) const
{
	auto carg = *this;
	FIPSTER_ASSERT(hedging);
	carg.hedging = hedging_arg_t(it,carg.hedging->skip,carg.hedging->delta);
	return carg;
}

option_node_factory::sender_ptr
option_node_factory::inner_setup(const arg_t& arg){

	auto default_n = default_node<discretized_space_field<1> >();
	auto& times = config->times;
	auto& time = arg.time;
	auto& id = config->id;
	auto& exercise_config = config->exercise;
	plain_spacetime_field exercise_field = 
		betterAt(objs->spacetime_fields,exercise_config.field_s,id);
	auto expectation = betterAt(objs->expectations,config->expectation_s,id);
	
	if(!config->market_arb_constant && config->entropic_theta){
		auto fdw_conf = betterAt(objs->fdweights_configs
														 ,expectation->config->fd_weights_s
														 ,config->id);
		
		auto pde=dynamic_pointer_cast
			<const finite_difference_weights::BSiso>(fdw_conf->PDE);

		config->market_arb_constant = (config->discounting_rate - pde->drift)
			/ pde->volatility / pde->volatility / *config->entropic_theta;
	}
		
	// cout<<"OP "<< arg.time << ": " <<arg.carg <<endl;
	
	time_limitation(time,times.stop,id);
	time_supported(times.contains(time),time,id);
	
	//with possible exercise?
	bool do_exercise = (arg.carg.exercise == always
											|| (arg.carg.exercise == expiration
													&& times.stop==time))
		&& (!times.start || arg.time >= *times.start);

	//\todoI implicit also needs marketarb somehow, thats why it is
	//deactivated right now (see todoI)

	auto exp_time = times.next(time,times.stop);
	auto rate_per_btick = config->discounting_rate
		*global_config::get().years_per_btick;
	
  auto discounting_factor = exp(-rate_per_btick*exp_time);
	
	// generate time_step
	expect_future_t fut
		(exp_time, times.stop
		 ,do_exercise && exercise_config.implicit ? 
		 exercise_field : 0
		 ,config->entropic_theta
		 ,discounting_factor);

	auto explicit_exercise_field = do_exercise && 
		(!exercise_config.implicit || time==times.stop) 
		? exercise_field->get_node(arg.strip()) :
		default_n;
											
	auto marketarb = do_exercise && time!=times.stop 
		//why not at time.stop? because it is only required for the
		//exercise-or-not comparison before the last step. at the last
		//step, both alternatives contain p(H) and it is easier to
		//calculate this last p(H) explicitely (by using the normal
		//transaction cost mechanism already used in the other
		//steps. furthermore, this approach allows to use identical
		//routines to calculate the marketarb and all other oprices.
		? get_node(arg.endow(arg.carg.to_marketarb())) : 
													 //CodeRefMA todoPerformance: instead of calculation the price of an
													 //actual zero claim, which could be cached for different options,
													 //this calculates the price of the specific claim with exercise
													 //disabled.
		default_n;

  spacetime_field<J_SGS>::sender_ptr delta_source=default_n;

	vector<spacetime_field<J_SGS>::sender_ptr> continuations;

	bool infinite = false;
	
	if(time!=times.stop){ //the last exercise/hedging time
		//generate vector of senders for every hedging grid point and power
		auto ps = config->powers;
		auto reg_conts = [&](){
			for(uint p=1;p<=ps;++p){
				fut.details.power = p;
				continuations.push_back(expectation->
																get_node(ev_field_arg(arg.strip(),fut)));
			}
		};
		if(!arg.carg.hedging){
			fut.details.hedging = boost::none;
			fut.field = create_proxy(arg.carg);
			reg_conts();
		}else{
			auto h = config->hedging();
			auto ha = *arg.carg.hedging;
			FIPSTER_ASSERT(h || !"option configured for hedging. (cf. "
											"option_carg_t's constructor)");

			infinite =  !config->times.contains(time,ha.skip) || h->infinite;
				// thread_logger()<<config->times<<" time: "<<time<<
				// 	" ha.skip: "<<ha.skip<<" h->infinite: "<<h->infinite<<endl;

			expectation->config->check_arbitrage
				(config,objs,fut.field_time,arg.time);
		
			//\todo Performance: it is inefficient to calculate future values
			//for all hedging ratios, given that only one will be
			//selected. but how to use the results of the delta_source, to
			//select a hedging ratio to use in fut.details.hedging????
			if(ha.delta && *ha.delta){
				FIPSTER_ASSERT(do_exercise);
				delta_source = (*arg.carg.hedging->delta)->
					             source(objs).get_node(time,arg.grid);
			}

			auto& hgrid=*h->grid;
			if(hgrid.D != arg.grid->D)
				BOOST_THROW_EXCEPTION
					(runtime_error("Hedging only implemented for identical grid sizes ("
												 + hgrid.id + ", " + arg.grid->id+")"));	
			for(grid_iterator<0,true> i(hgrid);i.valid();++i){
				fut.details.hedging = i.state();
				fut.field = create_proxy
					(arg.carg.change_hedging_index(i));
				reg_conts();
			}
		}

	}
			
	return create_node
		(body_t(arg,config,infinite)
		 ,make_tuple
		 (new_runtime_join_node(continuations).get()
			,explicit_exercise_field.get()
			,marketarb.get()
			,delta_source.get()
			));

	// //select european or american stepping
	// {
	// 	auto& e=*(option_config::exercise_t*)0; //war c.exercise; 
	// 	if( e.implicit )
	// 		//Sender: underlying field (argument is simply forwarded)
	// 		make_edge(*betterAt(objs->spacetime_fields
	// 												,e.field_s
	// 												,config->id)->get_node(arg)
	// 							,input_port<0>(*joinnode));
	// 	else
	// 		make_edge(*default_node<discretized_space_field<1> >()
	// 							,input_port<0>(*joinnode));
	// }
			
	// cout << "ASD: "<<arg.carg.size()<<" " <<arg.carg << endl << endl;
}

//forward to exercise value. this is not part of some consistent
//scheme. //\todo
option_node_factory::bv_sender_ptr
option_node_factory::boundary_setup(const arg_t& arg){
	auto& time = arg.time;
	auto& times = config->times;
	auto& id = config->id;
	
	time_limitation(time,times.stop,id);

	//redirect all requests to the stop time see \todo4
	if(time!=times.stop){
		auto arg2(arg);
		arg2.time=times.stop;
		return bv_get_node(arg2);
	}

	auto& exercise_config = config->exercise;
	plain_spacetime_field exercise_field = 
		betterAt(objs->spacetime_fields,exercise_config.field_s,id);
	auto explicit_exercise_field = arg.carg.exercise != never
		? exercise_field->bv_get_node(arg.strip()) :
		default_node<boundary_field_t<J_SGS> >();
	return create_node(boundary_body_t(arg,config)
										 ,explicit_exercise_field.get());
}

void options::read_referenced_fields(const string& id
																		 ,const grid_ptr& grid
																		 ,const ptree& fields_tree
																		 ,vector<referenced_field_t>& fields 
																		 ,const shared_const_objs_t& objs){

	for(auto it : fields_tree){
		if(is_xml_special(it.first)) continue;
		auto grid1d_name=it.second.get_optional<string>
			("<xmlattr>.calculateAsGeomAverageOn");
		auto grid1d = grid1d_name	? betterAt(objs->grids,*grid1d_name,id)
			: nullptr;
		if(grid1d) FIPSTER_ASSERT(grid && grid->D==1);

		auto fname =betterGet<string>(it.second,"<xmlattr>.id",id+".fields");
		if(it.first=="field"){
			fields.emplace_back
				(fname,betterAt(objs->spacetime_fields,fname,id),grid1d,boost::none);
		}else if(it.first=="option"){
			auto option=betterAt(objs->options,fname,id);

			Exercise exercise = never;
			auto exS = it.second.get<string>("<xmlattr>.exercise");
			if(exS=="always") exercise=always;
			else if(exS=="expiration") exercise=expiration;
			else FIPSTER_ASSERT(exS=="never");

			auto context = id+".fields.option."+fname;

			auto delta_pt = it.second.get_child_optional("deltaHedging");
			auto delta_hedging = delta_pt ? boost::make_optional
				(exercise==never ? nullptr : make_shared<delta_hedging_t>(*delta_pt,context))
				: boost::none;

			auto hedging_bands = it.second.get("<xmlattr>.hedgingBands",false);
			
			auto hedging = option->config->hedging(objs);
			if(hedging){
				auto hgrid=hedging->grid;
				auto poss = read_vecs<int>("at","n",it.second,context);
				list<grid_iterator<0,true> > hs;
				for(auto& p : poss){
					hs.emplace_back(*hgrid,p);
					hs.back().check_bounds();
					FIPSTER_ASSERT(p==hs.back().index().get_position(*hgrid));
				}
				if(hs.empty())  // use every possible hedging position if non are given
					for(grid_iterator<0,true> i(*hgrid);i.valid();++i)
						hs.push_back(i);

				int n_h = hs.size();

				if(hedging_bands)  // use every possible hedging position for calculation of hedging bands
					for(grid_iterator<0,true> i(*hgrid);i.valid();++i)
						hs.push_back(i);

				uint skip = it.second.get("<xmlattr>.skip",1);
				FIPSTER_ASSERT(skip>0);
					
				for(auto& h : hs){
					auto fname2 = fname+"("+exS+")";
					boost::optional<string> name;
					switch(n_h){
					case 0:
						name = "lowerOptimal_{" + fname2 + "}";
						break;
					case -1:
						name = "upperOptimal_{" + fname2 + "}";
						break;
					case -2:
						name = "lowerBound_{" + fname2 + "}";
						break;
					case -3:
						name = "upperBound_{" + fname2 + "}";
						break;
					default:
						if(n_h>0){
							auto a= fname+"("+toS(h.state())+","+exS+
								(delta_hedging?",deltaH(skip="+toS(skip)+")":"")+")";
							name = a+"\tH_{"+a+"}\tE_{"+a+"}";
						}
					}
					
					fields.emplace_back
						(name, option->
						 create_proxy(option_carg_t
													(exercise,hedging_arg_t
													 (h,skip,delta_hedging)))
						 ,grid1d,referenced_field_t::hedging_t(hgrid,h.index(),n_h <= 0)); //\todo CodeRefBadHack
					--n_h;
				}
			}else{
				auto a= fname;
				fields.emplace_back
					(a+"\tempty\tE_{"+a+"}", option-> create_proxy
					 (option_carg_t(exercise)),grid1d,boost::none);
			}
		}
	}
}

referenced_field_t::referenced_field_t(const boost::optional<string>& id,
																			 plain_spacetime_field field,
																			 grid_ptr grid,
																			 const boost::optional<hedging_t>& hedging)
	: id(id),field(field),grid(grid),hedging(hedging)
{}
