#pragma once
#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/property_tree/ptree.hpp>
#include <comphash.h>
#endif

#include "spacetime_field.h"
#include "grid.h"
#include "set_of_times.h"

namespace fipster {

	// #################################################
	// ##############   PHASE 1:	      		############
	// ##############   Config              ############
	// #################################################
	namespace options {
		
		using namespace std;
		using namespace boost::property_tree;
	
  	struct referenced_field_t{
  		boost::optional<string> id;
  		plain_spacetime_field field;
  		grid_ptr grid; //optional
  
			// tuple (hedging grid, hedging index, use this result for
			// hedging band calculation only (\todo CodeRefBadHack)
  		typedef tuple<grid_ptr,grid_index<0>,bool >	 hedging_t;
  
  		boost::optional<hedging_t> hedging;
  
  		referenced_field_t(const boost::optional<string>& id,
  											 plain_spacetime_field field,
  											 grid_ptr grid,
  											 const boost::optional<hedging_t>& hedging);
  
  		
  		spacetime_field<J_SGS>::sender_ptr get_node(btime time, grid_ptr grid2){
  			return field->get_node
  				(field_arg_t<>(time,grid ? grid : grid2,nothing()));
  		}
  	};
		
    struct delta_hedging_t {
			string id;
			referenced_field_t source(shared_const_objs_t objs);
			
			delta_hedging_t(const ptree& pt, const string& contetx);
			
		private:
			ptree pt;
			boost::optional<referenced_field_t> my_source;
		};
		
		struct hedging_t{
			string id;
			double k0,k1; // transaction costs
			bool infinite;
			grid_ptr grid;

			void init(shared_const_objs_t objs);
			
			hedging_t(const ptree& pt,string c);

			bool initialized;
		private:
			const string grid_s;
		};
		
		struct exercise_t{
			string field_s;
			bool implicit;
			
			exercise_t(const ptree& pt,string c);
		};
		
		struct option_config_t{
		
			const string id, expectation_s;
			exercise_t exercise;
			combA times;
			boost::optional<double> entropic_theta;
			double sigma_factor;
			uint powers;
			double discounting_rate;
			boost::optional<double> market_arb_constant;
			
			const boost::optional<hedging_t>& hedging(shared_const_objs_t objs);
			const boost::optional<hedging_t>& hedging() const;
			
			//ctor
			option_config_t(const ptree& pt);

		private:
			boost::optional<hedging_t> my_hedging;
		};

		//get field pointers and (optional) 1D grids from a ptree
		void read_referenced_fields(const string& id
																,const grid_ptr& grid
																,const ptree& fields_tree
																,vector<referenced_field_t>& fields 
																,const shared_const_objs_t& objs);
	}

	using options::hedging_t;
	using options::delta_hedging_t;
	using options::exercise_t;
	using options::option_config_t;
	using options::read_referenced_fields;
	using options::referenced_field_t;

	typedef shared_ptr<const exercise_t> exercise_ptr;
	// #################################################
	// ##############   PHASE 2:		      	############
	// ##############   Node Factory        ############
	// #################################################
	enum Exercise  { never,
									 expiration,
									 always };

    struct hedging_arg_t : comphashable<hedging_arg_t>
											 , stringifyable<hedging_arg_t>{
			state_t state;
			grid_index<0> ind;
			uint skip;

		
			// none: the optimal hedging strategy is selected
		  // nullptr: use 0 hedge	(for marketarb/zero claim)
			// pointer: delta hedging from given field
			boost::optional<shared_ptr<delta_hedging_t> > delta;

			template<class T> FORCE_INLINE void my_combine(T& t)const{
				t(&hedging_arg_t::ind
					,&hedging_arg_t::skip
					,&hedging_arg_t::delta);
			}

			hedging_arg_t(grid_iterator<0,true> it
										, uint skip
										, boost::optional<shared_ptr<delta_hedging_t> > delta)
				: state(it.state()), ind(it.index()), skip(skip), delta(delta) {}
		};

	struct option_carg_t : comphashable<option_carg_t>
											 , stringifyable<option_carg_t>
	{
		Exercise exercise;
		
		//none: if the result does not depend on the hedging position
		//      (i.e. is constant with respect to that parameter)
		boost::optional<hedging_arg_t> hedging;

		option_carg_t(Exercise ex
									,boost::optional<hedging_arg_t>&& hedging_arg = boost::none)
			: exercise(ex)
			, hedging(hedging_arg){}

		option_carg_t to_marketarb() const;

		option_carg_t change_hedging_index(grid_iterator<0,true> it) const;

		template<class T> FORCE_INLINE void my_combine(T& t)const{
			t(&option_carg_t::exercise
				,&option_carg_t::hedging
				);
		}
	};

	typedef shared_ptr<option_config_t> option_config_ptr;

	//pairs of (exercise or not, hedging position)
	typedef pair<bool, boost::optional<grid_index<0> > > decision_t;
	typedef discretized_space_field_plus<decision_t> option_field_t;
	typedef field_arg_t<option_carg_t> option_field_arg;

	/** \brief results are actually of type option_field_t
	 */
	struct option_node_factory : spacetime_field<J_SGS,option_carg_t>
	{
		option_config_ptr config;
		shared_const_objs_t objs;
			
		option_node_factory(option_config_ptr c,
												shared_const_objs_t o)
			: config(c),objs(o) {}
			
		sender_ptr inner_setup(const arg_t& arg);
		bv_sender_ptr boundary_setup(const arg_t& arg);
	};

// #################################################
// ##############   PHASE 3:		      	############
// ##############   Node Bodies         ############
// #################################################

	namespace options {
struct body_t : auto_connecting_body
<body_t
 //result
 ,discretized_space_field<J_SGS>
 //future 0:
 , vector<spacetime_field<J_SGS>::result_sptr_t> // expectations for different hedging ratios/powers
 //future 1:
 ,discretized_space_field<1> //exercise value (optional)
 //future 2:
 ,discretized_space_field<1> // market arbitrage price (optional) //codeRefMA
 //future 3:
 ,discretized_space_field<1> // field to give the delta for hedging (optional)
 >
{
	option_field_arg fieldargs;
	shared_ptr<const option_config_t> config;
	bool infinite;
	result_sptr_t operator()(const input_t& futures);
	
	body_t(option_field_arg f
				 ,option_config_ptr c
				 ,bool infinite);
};
struct boundary_body_t : auto_connecting_body
<boundary_body_t
 //result
 ,boundary_field_t<J_SGS> 
 //future 0:
 ,boundary_field_t<J_SGS>  //exercise value
 >
{
	option_field_arg fieldargs;
	shared_ptr<const option_config_t> config;
	result_sptr_t operator()(const input_t& futures);
	
	boundary_body_t(option_field_arg f
				 ,option_config_ptr c);
};
}
}
