#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION

#endif

#include "timestepping_operator.h"
#include "runtime_join_node.h"
#include "logging.h"

namespace fipster { namespace _timestepping_operator {
	
	using namespace Eigen;

	// #################################################
	// ##############   PHASE 2:			############
	// ##############   Node Factory        ############
	// #################################################
	
	multiply_adder_body_t::
		multiply_adder_body_t( array<double,2> betas, bool strang_symmetrization )
		:betas(betas),strang_symmetrization(strang_symmetrization)
	{}


	timestepping_operator::sender_ptr 
	timestepping_operator::setup_node(const arg_t& arg){
		auto split = splitter::sget_node(get<2>(arg));
		
		auto sender= make_tuple(split.get(),					//future 0 (tridiagonally_split_grid_operator)
								get<0>(get<2>(arg).first),	//future 1 (B_t)
								get<1>(get<2>(arg).first)		//future 2 (fd_weights_result_t	)
								);

		return create_node(
			multiply_adder_body_t(	get<0>(arg),//betas
									get<1>(arg)	//strang
								),sender);
	}
	
	splitter::sender_ptr
	splitter::setup_node(const arg_t& arg){
		return create_node(splitter_body_t(*arg.second),arg.first);
	}

	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################
	
	multiply_adder_body_t::result_sptr_t 
		multiply_adder_body_t::operator()(const input_t& futures)
	{
		if(0)thread_logger()<<"MUTLADD| START"<<endl;

		//aliases
		const auto& split_op = *get<0>(futures);
		const B_ptr_t&				B_ptr	= get<1>(futures);
		const fd_weights_result_t&	weights	= *get<2>(futures);
		
		//create result:		
		auto result = make_shared<timestepping_op_t>
			(B_ptr,weights.grid_striding_infos,betas);
		
		for(uint beta=0;beta<2;beta++){
			if(betas[beta]!=0.0)
				
				//loop over directions:
				for(auto dir = result->op[beta].dir_it();dir.valid();++dir)
				{
					//TODO: is it correct to do this? (It was missing)
					//if strang is used, every but the last direction are split in half.
					double actual_beta = betas[beta]*
							(strang_symmetrization && !dir.is_last() ? 0.5 : 1 );
					
					//loop over systems:
					auto source_sys=split_op.at(dir).begin();
			
					for(auto sys=result->op[beta].at(dir).begin();
						sys.valid();++source_sys,++sys){
						//loop over system:
						int l=sys.connected_set_it->length;
						for(int i=0;i<l;i++){
								sys.M.diag_ref(i)		=1.0 + actual_beta*source_sys.M.diag(i);
								//dont use beta for i==0 and i==l-1:
								//\todo why?
								//Copy A_o without modification (see CodeRefA)
								//\todo performance: make this check at compiletime?
								sys.M.subdiag_ref(i)	=	   (i>0?actual_beta:1)*source_sys.M.subdiag(i);
								sys.M.superdiag_ref(i)	=	   (i<l-1?actual_beta:1)*source_sys.M.superdiag(i);
						}
					}

					;
				}
		}

		if(0)thread_logger()<<"MUTLADD| DONE"<<endl;
		return result;
	}

}}
