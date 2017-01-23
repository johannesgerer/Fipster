#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <tuple>
#include <tbb/combinable.h>
#include <cmath>
#endif

#include "logging.h"
#include "expectation_value_body.h"
#include "cryer.h"
#include "Math.h"
#include "tbb_flow_graph.h"
#include "grid_fields.h"
#include "tridiagonally_split_grid_operator.h"
#include "expectation_body_stor.h"
#include "expectation_value_node_factory.h"

using namespace fipster;
using namespace expectation_values;
	
	
	
body_t::body_t(const config_ptr& config
							 ,const current_step& cur_step
							 ,const ev_field_arg& fieldargs
							 , const string& meta_info,
							 double vola,btime time)
	:id(config->id),cur_step(cur_step)
	,fieldargs(fieldargs),meta_info(meta_info),vola(vola)
	,time(time)
{}

// #################################################
// ##############   PHASE 3:			############
// ##############   Node Bodies         ############
// #################################################
	
template<bool g0_is_zero,non_zero_t non_zero>
void perform_step2(	const body_t::input_t& futures,
										body_t::result_t& result,
										const body_t& body);

template<bool g0_is_zero>
void perform_step(	const body_t::input_t& futures,
										body_t::result_t& result,
										const body_t& body){
	if(get<2>(futures)->betas[0]==0)
		perform_step2<g0_is_zero,only_second>(futures,result,body);
	else if(get<2>(futures)->betas[1]==0)
		perform_step2<g0_is_zero,only_first>(futures,result,body);
	else
		perform_step2<g0_is_zero,both>(futures,result,body);

}
tbb::combinable<storage> tls;
	
typedef discretized_space_field<J_SGS> dsf;
		
/** add the for each state the scalar product of a with the state and
		raise to power. 

    \note this is performed for EVERY time step. However, the
node_factory will set power=1 if it is not the last time step.

factor = -1 is used for the ask price
 */
template<int factor=1,class A,class B>
shared_ptr<const B>
final_transform(grid_ref g
								,const shared_ptr<const B >& payoff
								,const details_t& more
								,const body_t& body
								,const bool boundary
								,const double correction=1){

	auto& hedge_pos=more.hedging;
	auto& e_theta=more.entropic_theta;
	auto& power=more.power;

	// printf("%.17e\t%.17e\t%.17e\n",discounting_factor,body.fieldargs.carg.discounting_rate,cstep.timestep);

	if(!hedge_pos && power==1 && more.discounting_factor == 1
		 && correction == 1
		 && !e_theta) return payoff;

	double dividend_factor = 1 ? 1 : exp(0.1e7*1e-8*.05);
	
	auto result=make_shared<B>("");
	auto gend=result->resize(g);
	A gbegin(g);
	for(auto i = gbegin;i!=gend;++i){
		double v = more.discounting_factor
			* (payoff->at(i.index()) + 
				 ( hedge_pos ? factor*hedge_pos->dot(i.state())*dividend_factor : 0));
				
		//\todo: Fix. this was only a quick hack with static values sigma=.5, drift-r=.01
		double optimalZeroClaim =
			0 && (body.fieldargs.carg.rannacher || boundary) && e_theta && i.state()(0) > 0 ?
			-log(i.state()(0))/.25*.03 : 0;
			
		result->at(i.index()) = correction *
			( e_theta ? expm1(optimalZeroClaim - *e_theta*v)
				: pow2(v, power)); //CodeRefExp

		// 	thread_logger()<<boundary<<": "<<optimalZeroClaim<<endl;
		if(0 && body.time == 99890000 && body.fieldargs.carg.rannacher){
		// 			// if(v!=0)
		// 	thread_logger()<<"expect: "<<payoff->at(i.index())<<"v: "<<v<<endl;
		// 	// if(i == gbegin.index())
				thread_logger()//<<optimalZeroClaim<<" "<<hedge_pos<<" e_theta"<<e_theta<<endl
					<<boundary<<" "<<body.meta_info<<endl;
				thread_logger()<<i.state()(0)<<", hedge_pos="<<hedge_pos
											 <<", v="<<v<<", r="<< result->at(i.index()) <<endl;
		}
	}
	return result;
}



/** performs all fractional time steps for all systems
 */
template<bool g0_is_zero,non_zero_t non_zero>
void perform_step2(	const body_t::input_t& futures,
										body_t::result_t& result,
										const body_t& body)
{
	//deduce the index of of the beta slot that is not zero (there
	//might be two)
	const int first_non_zero_access = non_zero==only_second?1:0;
	// array<bool,2> is_nonzero={ non_zero!=only_second , non_zero!=only_first };

	//aliases
	auto& underlying_field	=  get<0>(futures);
	auto& future_val			  =  get<1>(futures);
	auto& ts_ops		       	= *get<2>(futures);
	auto bc_inhoms = make_tuple(get<3>(futures), get<4>(futures));
	// cout<<"FUT USE COUNT: "<<get<1>(futures).use_count()<<endl;
	// cout<<"UNDER USE COUNT: "<<get<0>(futures).use_count()<<endl;
	auto& gsic				= *ts_ops.op[first_non_zero_access].gsic;
	auto& g					= *gsic.grid;
	auto& cstep = body.cur_step;

	//\todo4 bad hack, which is required due to the delegation of
	//boundary values to the future time. see \todo4

  // add hedging position (optional) and raise to power or exponentiate CodeRefH 
	auto future_exp = final_transform<-1,grid_iterator<1,true>,dsf>
		(g,future_val,cstep.now,body,false);

	result.resize(g);

	//get a direction iterator
	auto dir_it	= ts_ops.op[first_non_zero_access].dir_it();
		
	//deduce number of passes (two passes for strang in more than on dimension)
	uint npass = cstep.strang_symmetrization && ts_ops.op[first_non_zero_access].nDirs()>1 
		? 2 : 1;

	//get storages:
	storage &stor = tls.local();
		
	double correction = cstep.future.power == 2 ? exp(body.vola*body.vola) : 1;
			//todoSigma: still good idea?
	
	// if(body.fieldargs.carg.rannacher){
	// 	thread_logger()<<cstep.now.hedging<<" "<<cstep.future.hedging<<endl;
	// 	thread_logger()<<cstep.now.entropic_theta<<" "<<cstep.future.entropic_theta<<endl;
	// }
	//\todo CodeRefBad: very inefficient, works only for dirichlet,...
	auto bc_inhoms2 = 0 ? bc_inhoms : make_tuple
		( final_transform<-1,boundary_iterator<0,false,true>,boundary_field_t<J_SGS> >
			(g,get<0>(bc_inhoms),cstep.future,body,true,correction)

		,final_transform<-1,boundary_iterator<0,false,true>,boundary_field_t<J_SGS> >
			(g,get<1>(bc_inhoms),cstep.future,body,true,correction)
			);

	if(0 && cstep.now.hedging){
	thread_logger()
		<<"hedging:\n"<< *cstep.now.hedging
		<<"\nboundary :\n"<<**get<0>(bc_inhoms)
		<<"\nboundary2:\n"<<**get<0>(bc_inhoms2)
		<<"\nfi       :\n"<<**future_exp <<endl <<endl;
	}
	//calculate Bo1 * ( beta[1]*BC_inhom[1] - beta[0]*BC_inhom[0] )
	//see (amoung others) Scan 12-05-07 Final LCP and how Bo1_BC_inhoms is used later in this procedure.
	//The use of first_non_zero_access in the following lines, prevents the need to initialize Bo1_BC_inhom_temp:
	*  stor.Bo1_BC_inhom_temp  = (first_non_zero_access==0 ? -1 : 1) * ts_ops.betas[  first_non_zero_access] * **get<  first_non_zero_access>(bc_inhoms2);
	if(ts_ops.betas[1-first_non_zero_access]!=0)
		*stor.Bo1_BC_inhom_temp -= (first_non_zero_access==0 ? -1 : 1) * ts_ops.betas[1-first_non_zero_access] * **get<1-first_non_zero_access>(bc_inhoms2);

	*stor.Bo1_BC_inhoms = *ts_ops.B->Bo1 * *stor.Bo1_BC_inhom_temp;
		
	//loop over passes (there are two for strang symmetrization)
	for(uint pass=0;pass< npass;pass++){
		//if the second pass starts, set iterator to the (d-1)-th,
		//where d is the number of directions.
		if(pass==1){
			--dir_it;--dir_it;
		}

		//loop over directions: (once)
		while(dir_it.valid()){

			double strang_factor = npass==2 && !dir_it.is_last() ? 0.5 : 1;

			//thread_logger()<<"pass "<<pass<<" dir "<<dir_it.at()<<endl;

			//this is to reserve enough space now, to prevent allocations later
			stor.resize(gsic.at(dir_it).max_length);

			array<tridiag_systems_collection::const_sys_iterator,2> sys;
			auto& sys_for_iteration = sys[first_non_zero_access];
				
			if(non_zero!=only_second)
				sys[0]=ts_ops.op[0].at(dir_it).begin();
			if(non_zero!=only_first)
				sys[1]=ts_ops.op[1].at(dir_it).begin();

			//loop over systems:
			for(;sys_for_iteration.valid();++sys_for_iteration){
	
				uint length=sys_for_iteration.M.length;
				stor.resize(length);

				//MAIN WORK BODY:
				//compare Scan 12-05-07 Final LCP!
						
				//Gather vectors

				//copy f1 to local storage:
				sys_for_iteration.gather(*future_exp,stor.f1);

				//copy g0 to local storage:
				if(!g0_is_zero) //this assures that underlying_field is not null
					sys_for_iteration.gather(*underlying_field,stor.g0);

				//set zero:
				stor.q.setZero();

				//Calculate beta[i]*A_o.Bo1.BC_inhom[i]
				//Compare Scan 11-10-21 BC fundamentals
				// CodeRefA
				//the choice of sys[i].M used here does not matter, as they
				//only differ in the factor beta that was used in their generation.
				//But as A_o is included in M without modification, just use
				//the first non zero M.
				//The actual selection of whose boundary values to include at all (f0's or f1's or both)
				//is performed above in the construction of Bo1_BC_inhoms.
				stor.q[0] = strang_factor*sys_for_iteration.M.subdiag(0)
					*stor.Bo1_BC_inhoms.at(sys_for_iteration.connected_set_it->begin_bi);
				stor.q[length-1] = strang_factor*sys_for_iteration.M.superdiag(length-1)
					*stor.Bo1_BC_inhoms.at(sys_for_iteration.connected_set_it->end_bi);
					
				//add M1.f1
				if(0 || non_zero==only_first) //is M1==1?
					stor.q += stor.f1;
				else
					sys[1].M.times_vec<VectorXd,1>(stor.q,stor.f1);
					
				if(0 || non_zero==only_second){ //if the first beta=zero (which is the case for explicit Euler)
					stor.f0= g0_is_zero ? stor.q : stor.g0.array().max(stor.q.array());
				}else{
					if(g0_is_zero){
						//solve M0.u0 -q = 0
						stor.tridiagonal(sys[0]);
					}else{
						//subtract M0.g0
						sys[0].M.times_vec<VectorXd,-1>(stor.q,stor.g0);
						//solve M0.u0 - q > 0, u0 > 0 and (M0.u0 - q).u0=0
						stor.cryer_lcp(sys[0]);
						//add g0
						stor.f0 += stor.g0;
					}
				}				
						
				sys_for_iteration.scatter(result,stor.f0);
						
				if(1 && non_zero==both)
					++sys[1-first_non_zero_access];
			}

			//depending on the current pass (forward or backward)
			//the dir iterator is progressed
			pass==0 ? ++dir_it : --dir_it;
		}
		/*else{ WHY?
			++dir_it;++dir_it;
			}*/
	}
}

body_t::result_sptr_t body_t::operator()(const input_t& futures)
{
	auto& a=cur_step.now.hedging;
	// logging
	if(0)thread_logger()<<id<<"("<<fieldargs.time<<
				 ","<< (a?toS(*a):"")<<","<< cur_step.now.power <<" )| START"<<endl;

	auto result = make_shared<result_t>(meta_info);
		
	//Siehe 12-02-29-7 inhomogeneous coordinates final.jpg
	
	//select european or american stepping (for european g0 is set to zero)
	if(!get<0>(futures))
		perform_step<true>(futures,*result,*this);
	else
		perform_step<false>(futures,*result,*this);

	if(0)thread_logger()<<id<<"| DONE"<<endl;
	return result;
}


boundary_values_body_t::result_sptr_t 
boundary_values_body_t::operator()(const input_t& futures){

	if(0)thread_logger()<<"BV     | START"<<endl;

	//aliases:
	auto& x=*get<0>(futures);
	auto& BC=*get<1>(futures);
	auto& BCinhom=*get<2>(futures);

	auto result = make_shared<boundary_field_t<J_SGS>>();
		
	**result = *BC.Bo1**BCinhom - *BC.Bi**x;

	if(0)thread_logger()<<"BV     | DONE"<<endl;
	return result;
}
