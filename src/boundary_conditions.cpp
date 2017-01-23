#include "precompiled.h"
#include "boundary_conditions.h"
#include "boundary_condition_iteration.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/throw_exception.hpp>
#include <tbb/tick_count.h>
#include <pretty_printer.h>
#endif

#include "logging.h"
#include "grid_iterators.h"
#include "logging.h"
#include "sparse_sherman_morisson.h"
#include "utility.h"
#include "global_config.h"

namespace fipster { namespace _boundary_conditions {
	
	using namespace Eigen;
	
	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################
	struct bc_line_handler : boundary_conditions_line<body_t::sg>{
		Bo_t& Bo;
		Bi_t& Bi;

		bc_line_handler(Bo_t& Bo, Bi_t& Bi)
			:Bo(Bo),Bi(Bi)
		{}

		void add_to_interior(	boundary_index<sg>& row,grid_index<sg+1>& col,
									double v, uint, bstate_t){
					Bi.insert(row,col)=v;
		}
		void add_to_boundary(	boundary_index<sg>& row,boundary_index<sg>& col,
								double v, uint, bstate_t){
			Bo.insert(row,col)=v;
		}

		void set_zero(	boundary_index<sg>&,uint variable,bstate_t state){
			//nothing to do for a zero row in a sparse matrix
		};
	};


	body_t::result_sptr_t body_t::operator()( const input_t& )
	{
		//logging
		if(1)thread_logger()<<"BCs    | START"<<endl;
		using tbb::tick_count;

		tick_count t0 = tick_count::now();

		auto& g= *arg.second;
		auto result = make_shared<result_t>();
		auto& Bi=result->Bi.resize(g);
		auto& Bo1 = result->Bo1.resize(g);
		auto& Bo = result->Bo.resize(g);
		
		bc_line_handler bcl(Bo,Bi);

		for(auto bcit=boundary_condition_iterator<sg,&bc::type>(arg);
				bcit.valid(); ++bcit)
			bcit.handle_bc_line(bcl);

		//cout<<_sparse_sherman_morisson::test_sparse_sherman_morisson(100)<<endl;

		if(0)
			invert_using_rowwise_sparse_sherman_morisson(*Bo1,*Bo,false);
		else
			invert_using_single_entry_sparse_sherman_morisson(*Bo1,*Bo,!FIPSTER_NDEBUG);
		
		//exit(0);
		if(global_config::get().test_bc_inversion){
			double BCinversion_norm = inversion_residuum(*Bo1,*Bo);
			if(BCinversion_norm>arg.first.tolerance)
				BOOST_THROW_EXCEPTION(runtime_error("BC residuum norm "+toS(BCinversion_norm)+" > BC config tolerance"));
		}

		//Eigen Bug: (A*B).pruned() did not only not prune the result, 
		// instead id added zero entries(?!!!!!)
		//*Bi = (*Bo1 * *Bi).eval(); ///\todo TEST : is the product really the product
		
		if(1)thread_logger()<<"BCs    | DONE, assembled and inverted: "<<scientific <<setprecision(3)<<
			(tick_count::now()-t0).seconds()<<"s"<<endl;
			
		return result;
	}

	body_t::body_t(const bc_arg_t& arg)
		: arg(arg)
	{}

}}
