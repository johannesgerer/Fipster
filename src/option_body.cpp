#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <limits>
#endif

#include "utility.h"
#include "option.h"
#include "expectation_value_config.h"
#include "global_config.h"
// #include "interpolation_hedging.h"

using namespace fipster;
using namespace options;

body_t::body_t(field_arg_t<option_carg_t> f
							 ,option_config_ptr c
							 ,bool infinite)
	: fieldargs(f),config(c),infinite(infinite)
{}

boundary_body_t::boundary_body_t(field_arg_t<option_carg_t> f
							 ,option_config_ptr c)
	: fieldargs(f),config(c) {}

typedef  pair<bool,boost::optional<grid_index<0> > > add_t;

template<class Iterator,class Body,class Res>
struct action{

	typedef const typename Body::result_sptr_t field_t;

	const Body& body;
	const boost::optional<hedging_arg_t>& hed0;
	const option_config_t& config;
	const boost::optional<hedging_t>& hedging;
	shared_ptr<Res> result;
	field_t& exercv;
	

	action(Body& body, field_t& exercv)
		:body(body)
		,hed0(body.fieldargs.carg.hedging)
		,config(*body.config)
		,hedging(config.hedging())
		,result(new Res("option body result with decision info bla"))
		,exercv(exercv)
	{}
	
	void	last_step(){
		auto& g=*body.fieldargs.grid;
		auto gend = result->resize(g);
		bool zero_tc = !hed0 || 
			hed0->state == state_t::Constant(hedging->grid->D,0);
			// see CodeRefTHR:
			//((hed0.state - state_t::Constant(config.hedging->grid()->D,0)).array().abs()).matrix().maxCoeff() < 1e-15;
		
		//calculate max{g,0} + K_n (where hed1 = 0)
		for(Iterator i(g);i!=gend;++i)
			//CodeRefBO
			exercise(i, 0, boost::none
							 ,zero_tc ? 0 : transaction_cost(i,0)
							 , field_t());
	}

	FORCE_INLINE
	double transaction_cost(const Iterator& i,double hed1){
		auto& hed_conf = *hedging;
		auto& k0 = hed_conf.k0;
		auto& k1 = hed_conf.k1;
		return k0 + ( k1>0 ? 
									//proportional to the change in total value of the position
									k1 * abs(hed1 - hed0->state.dot(i.state())) : 0);
	}


	FORCE_INLINE
	void exercise( const Iterator& i
								 , double c, const boost::optional<grid_index<0> >& hedg, double add
								, field_t& marketarb){
		double ev;
		bool ex = exercv; //check if there actually is an exercise value
		if(ex){ //CodeRefMA
			ev=exercv->at(i.index()) +  (marketarb ? marketarb->at(i.index()) : 0);
			ex = ev > c;
		}
		
		result->at(i.index()) = add + (ex ? ev : c);
		store_additional(i,make_pair(ex, hedg));
	}

  void store_additional(Iterator i, add_t&& a);
};

template<class Iterator,class Body,class Res>
void action<Iterator,Body,Res>::store_additional
(Iterator i, add_t&& a){
	result->additional.at(i.index())=std::move(a);
}

//template specialization for boundary (without additional information)
template<>
void action<boundary_iterator<0,false,true>
				 ,boundary_body_t
				 ,boundary_field_t<J_SGS> >::store_additional
(boundary_iterator<0,false,true> i, add_t&& a){}

boundary_body_t::result_sptr_t boundary_body_t::operator()(const input_t& exercv){
	// logging
	if(0)thread_logger()<<config->id<<"("<<fieldargs.time<<", "
											<<fieldargs.carg.hedging->state<<", "<<
				 fieldargs.carg.exercise<<")| START"<<endl;
	
	// auto& g=*fieldargs.grid;
	// auto& hed0 = fieldargs.carg;
	
	action<boundary_iterator<0,false,true>
				 ,boundary_body_t
				 ,boundary_field_t<J_SGS> > ac(*this,exercv);

	result_sptr_t const_result;
	
	if(1){ // || config->hedging()){ // todoH: was soll das?
		auto& times = config->times;
		
		if(fieldargs.time==times.stop){ //the last exercise/hedging time
			ac.last_step();
		}else{
			FIPSTER_ASSERT(!"implemented");
		}
		const_result = ac.result;
	}else{
		FIPSTER_ASSERT(!"implemented");
	}
	return const_result;
}

const int dimension=0;

//second order accurate Delta Calculation (see Mathematica/15-10-20\ finite\ differences.nb)
template<class T, class I, class I2, class H>
uint select_delta_hedge(const T& d_source, grid_ref g, grid_ref hg,
												const I& begin, const I2& endMinusOne, const I& cur,
												H& hed1,uint hedge_n, const discretized_space_field<1>& marketarb_h){

	double delta = marketarb_h.at(cur.index());
	

	if(d_source){
		auto& ds = *d_source;
	
		auto p = cur.pos()[dimension];
		double
			a = g.stepsizes[dimension][p],
			a2 = a*a,
			b = g.stepsizes[dimension][p-1],
			b2 = b*b;

			auto f = [&](uint i){return ds.at(cur.index()+i);};
	
		if(1){ //actually use the calculated delta?
			if(cur==begin)
				delta +=	(f(1) - f(0))/a;
			else if(cur==endMinusOne)
				delta +=	(f(0) - f(-1))/b;
			else
				delta +=	((a2 - b2)*f(0) + b2*f(1) - a2*f(-1))
					/(a*b*(a + b));
		}
	}
	
	double diff=numeric_limits<double>::max();
	grid_index<0> last_index;
		
	//find closest hedging position
	for(uint count = 0; count < hedge_n; ++count, ++*hed1){
		double h = hed1->state()[0];
		double temp = delta - h;
		if(temp>0 && count < hedge_n - 1){
			diff = temp;
			last_index = hed1->index();
		}else{
			if(diff < -temp){
				hed1=make_shared<grid_iterator<0,true> >(hg,last_index);
			  --count;
			}
	// thread_logger()<<scientific<<setprecision(17)
	// 							 <<cur.state()[0]<<"\t"
	// 							 <<delta<<"\t"<<hed1->state()[0]<<"\t"<<count<<endl;
								 // <<"\t"marketarb<<endl;
			// f(0)<<"\t"<<a<<"\t"<<a2;
			// <<"\t"<<b<<"\t"<<f(1)<<"\t"<<f(-1)<<"\t"<<delta<<endl;

			// thread_logger()<<"select_delta_hedge "<<d_source<<" "<<delta<<" "<<count<<endl;
			return count;
		}
	}
	FIPSTER_ASSERT(!"this should never be reached");
}

template<class I, class I2>
void generate_optiomal_zero_hedging(double factor, 
																		grid_ref g, 
																		const I& gbegin, const I2& gend,
																		discretized_space_field<1>& result){

	// thread_logger()<<"generate_optiomal_zero_hedging "<<t<<endl;

	result.resize(g);
	grid_iterator<1,true> i = gbegin;

	for(++i;i!=gend;++i){
		result.at(i.index()) = 0 ? 0 : factor/i.state()[0];
	}
	// thread_logger()<<t<<endl;
}


body_t::result_sptr_t body_t::operator()(const input_t& futures){
	// logging
	if(0)thread_logger()<<config->id<<"("<<fieldargs.time<<", "
											<<fieldargs.carg.hedging->state<<", "<<
				 fieldargs.carg.exercise<<")| START"<<endl;
	
	auto& g=*fieldargs.grid;
	auto& exps = *get<0>(futures); //expectations
	auto& exercv = get<1>(futures); //exercise value
	auto& marketarb = get<2>(futures);
	auto& delta_source = get<3>(futures);
	auto& hed0 = fieldargs.carg;

	double discounting_factor = exp
		(config->discounting_rate*fieldargs.time*global_config::get().years_per_btick);
	
	result_sptr_t const_result;
	action<grid_iterator<1,true>,body_t,option_field_t> ac(*this,exercv);
	
	auto& times = config->times;
	auto& e_theta = config->entropic_theta;
		
  FIPSTER_ASSERT(!delta_source || g.D==1);
	
	if(fieldargs.time==times.stop){ //the last exercise/hedging time
		ac.last_step();
	}else{
		auto gend = ac.result->resize(g);
		auto gendMinus1 = gend-1;
		grid_iterator<1,true> gbegin(g);
			
		uint hedge_n;
		shared_ptr<grid_iterator<0,true> > hed1_begin;
		grid_ptr hg;

		if(!hed0.hedging){
			hedge_n=1;
		}else{
			hg=config->hedging()->grid;
			hedge_n=hg->N[0];
			hed1_begin=make_shared<grid_iterator<0,true> >(*hg);
		}
		auto ps = config->powers;
		FIPSTER_ASSERT(exps.size() == hedge_n*ps);
		auto& sigma_factor = config->sigma_factor;

    // auto marketarb_o = dynamic_pointer_cast<const option_field_t>(marketarb);
		// FIPSTER_ASSERT(!delta_source || marketarb_o);
		discretized_space_field<1> marketarb_h("");
		// if(delta_source && marketarb_o)
		// 	interpolate_hedging(*marketarb_o,g,*hg,gbegin,gend,marketarb_h);
		uint max_count = hedge_n;
		if(hed0.hedging && hed0.hedging->delta && !this->infinite){
			generate_optiomal_zero_hedging(*config->market_arb_constant * discounting_factor
																		 ,g,gbegin,gend,marketarb_h);
			max_count=1;
		}
		
		//loop over all states and find optimal policy at each site
		// or keep current hedging position (for infinite transaction costs)
		// or change position according to externally given delta
		for(auto i = gbegin ;i!=gend;++i){

			shared_ptr<grid_iterator<0,true> > hed1;

			auto exp = exps.begin();

			if(hed0.hedging){
				hed1 =  make_shared<grid_iterator<0,true> >(*hed1_begin);

				if(hed0.hedging->delta && !this->infinite){
					exp += ps * select_delta_hedge
						(delta_source,g,*hg,gbegin,gendMinus1,i,hed1,hedge_n,marketarb_h);
				}
			}

			double min=numeric_limits<double>::max();
			boost::optional<grid_index<0> > min_i = boost::none;

			for(uint count = 0;count < max_count; ++count){

				if(count>0){
					hed1->operator++();
					exp+=ps;
				}
				
				// thread_logger()<<hed0.ind<<" vs "<<hed1->index()<<endl;
				
				double hv = 0;
				if(hed0.hedging){
					if(this->infinite && hed1->index() != hed0.hedging->ind)
						continue;
					hv = hed1->state().dot(i.state());
				}

				double con = (*exp)->at(i.index());
				double v= hv; //start with current position value

				if(sigma_factor!=0){
					double con2=(*(exp+1))->at(i.index());
					double temp=con2 - con*con;
					//add the expectation and standard deviation
					v+= discounting_factor
						*(con + sigma_factor*(temp>0 ? sqrt(temp) : 0));
					//\todo3 nan schutz, wegen fortpflanzung

					// cout<<v<<"\n"<<con<<"\n"<<con2<<"\n"<<sigma_factor<<"\n"<<endl;
				}else{
					v += discounting_factor
						*(e_theta ? -log1p(con)/ *e_theta : con); //CodeRefExp
				}


				if(hed0.hedging){
					if(hed1->index()!=hed0.hedging->ind) //if hedging position is changed
						v += ac.transaction_cost(i,hv);	
				}
				
					// cout<<v<<endl;
				// assert(v < 1e100);
				// assert(result->at(i.index()) < 1e100);
				if(v < min){
						min=v;
						if(hed0.hedging){
							min_i = boost::in_place(hed1->index());
						}
				}
				// if(0 && this->infinite)
				// 	thread_logger()<<"v: "<<v<<", min: "<<
				// 		min<<", hed1: "<<*hed1<<", min_i: "<<min_i<<endl;
				// if(hed1->index().ind==3 || hed1->index().ind==4 )
				// 	thread_logger()<<"v: "<<v<<", i: "<<hed1->index()<<" min: "<<min<<", min_i: "<<min_i<<endl;
			}
			ac.exercise(i, min, min_i,0,marketarb);
		}
	}
	const_result=ac.result;
	return const_result;
}

