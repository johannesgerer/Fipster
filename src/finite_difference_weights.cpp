#include "precompiled.h"
#include "finite_difference_weights.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <utility>
#include <tbb/tick_count.h>
#include <tbb/partitioner.h>
#include <boost/throw_exception.hpp>
#endif

#include "Eigen.h"
#include "utility.h"
#include "Math.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <pretty_printer.h>
#endif

#include "finite_difference_stencils.h"
#include "logging.h"
#include "finite_difference_weights_algo.h"
#include "grid_striding_info_inline.h"
#include "global_config.h"
#include "exceptions.h"
#include "PDE.h"

namespace fipster { namespace finite_difference_weights {

	// #################################################
	// ##############   PHASE 2:			############
	// ##############   Node Factories      ############
	// #################################################

	fdweights::sender_ptr fdweights::setup_node( const arg_t& arg )
	{
		auto PDE=arg.first.PDE;

		auto pde = dynamic_pointer_cast<const BSiso>(PDE);
		if(pde) return create_node(body_t<const BSiso>(arg,pde));

		// check for different PDEs
		// auto pde = dynamic_cast<const DifferenPDE* >(PDE);
		// if(pde) return create_node(body_t<const ...>(arg,pde));

		FIPSTER_THROW_EXCEPTION(runtime_error("PDE type "+PDE->type+" not implemented!"));
		return 0;
	}


	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################


	using namespace Eigen;

	template<class PDE_core>
	body_t<PDE_core>::body_t( const fd_weights_arg_t& args
														, shared_ptr<const PDE_core> pde_ptr ) 
		:	pde(pde_ptr),args(args),D(args.second->D)
	{}


	//############################################################################
	//Constructs the multi_indices corresponding to partial derivatives    
	int references::generatePDindices(uint order)
	{
		uint i__,i__1,i__2,i__3,k,j;
		bool more;
		/* Calculate number of partial derivatives to include */
		int numpds = 0;
		i__1 = order;
		for (i__ = 1; i__ <= i__1; ++i__) {
			i__2 = i__ + D - 1;
			i__3 = D - 1;
			numpds += (int) binom(i__2, i__3);
		}

		multi_indices.resize(numpds,D);
		position_t ind(D);
		more=0;
		j = 0;
		/* Calculate multi indices for partial derivatives */
		i__1 = numpds;
		int h,t;
		for (i__ = 0; i__ < i__1; ++i__) {
			if (! more) {
				++j;
				if (j > order)
					BOOST_THROW_EXCEPTION(invalid_argument(
					"Number of Partial derivatives is to high. Specified order exceeded"));
			}
			comp_next(j, D, &ind[0], &more, &h, &t);
			for (k = 0; k < D; ++k) {
				multi_indices(i__,k) = ind[k];
			}
		}
		
		return numpds;
	}

	//############################################################################
	// Compute line-wise constant factorials in denominator for weight optimization systems
	void references::generateFactorialDenoms(double HigherOrderWeight){

		factorialdenoms.resize(multi_indices.rows());

		uint j;	

		for (int i = 0; i < factorialdenoms.size(); ++i) {
			uint sum=0;
			factorialdenoms[i]=1.0;
			for (j=0; j < D; ++j) {
				factorialdenoms[i] /= factorial(multi_indices(i,j));
				if(multi_indices(i,j)!=0) sum++;
			}
			if (sum>2)
				factorialdenoms[i] *= HigherOrderWeight;
		}
	}

	template<class F,class G, class H, class I >
	void test_A_o(const F& A_o, const G& A_io, const H& gsi, const I& stencil_points){
		auto& grid=*gsi.grid;
		auto end_b = grid.template b_end<J_SGS>();			
		set<boundary_index<J_SGS-1>> obinds;
		for(auto bit = boundary_iterator<J_SGS>(grid); bit != end_b; ++bit)
		{
			obinds.clear();

			//loop over neighbors
			const auto s_end=end(gsi.sorted_neighs);
			for(auto s=begin(gsi.sorted_neighs);s!=s_end;++s){
				//get neighbor position
				auto pos=bit.pos()+stencil_points.at(s->second);
				//check, if points lies on outer boundary:
				boundary_index<J_SGS-1> obind;
				if(obind.make_ind_from_pos<true>(pos,grid)){
					FIPSTER_ASSERT(A_o.coeff(bit.b_ind,obind)
										==A_io.at(s->second).at(bit.g_ind));
					if(A_io.at(s->second).at(bit.g_ind)!=0)
						obinds.insert(obind);
				}
			}

			//compare number of nonzeros in this A_o row and the number of unique accesses
			FIPSTER_ASSERT(obinds.size()==
											(uint)A_o->innerVector(bit.b_ind.ind).nonZeros());
		}
	}

	//############################################################################
	/** Generate the finite difference weights for the selected stencil,
		to get a finite difference approximation of the differential operator.
		
		It's exactly done as described in  the paper 
		"Ito, K., & Toivanen, J. (2009). Lagrange Multiplier Approach with 
		Optimized Finite Difference Stencils for Pricing American Options 
		under Stochastic Volatility. SIAM Journal on Scientific Computing, 
		31(4), 2646. doi:10.1137/07070574X"
	
		With the following additions:
		1. The sign choices specified in PDE.h dicate that the weights 
			constitute a NEGATIVE M-MATRIX. 
		2. This is achieved by using non-negative least squares for the off-diagonal weights,
			as opposed to the proposed solution of a quadratic programming problem.
		3. Extension to higher dimensional PDEs
		4. Eq. 40 contains (dx_i)^delta_i. I use 0^0 = 1 (Reason: If there is no 
			derivative present in one direction (delta_i=0) then it is 
			ignored (i.e. factor 1). In all other cases a stepsize of 
			zero in one direction (dx_i=0) will yield a 0 in the matrix d.)
	
		 \todo verbesserung: wie vorgeschlagen, die partiellen Ableitungen
		 die man mit einem gröberen Gitter erhällt für die optimierung
		 benutzen...
	
	\return 1. A 2D array, where each row corresponds to one point in the multi-dim 
			grid and the i-th column represent the weights of the i-th stencil point,
			if sorted lexicically by its direction vector. The last weight 
			corresponds to the unfinished (see note) center weight.
			2. The so called "grid_striding_info"

	\note To get the final center weight all non-center weights have
			to be subtracted.

	*/
	template<class PDE_core>
	typename body_t<PDE_core>::result_sptr_t body_t<PDE_core>::operator()( const input_t& )
	{
		const grid_ptr& grid = args.second;
		
		//gcc: was editable_result(new result_t) (i.e. possibly without initizalitaion?)
		fd_weights_result_ptr_t result = make_shared<result_t>(); 

		//create object to hold stuff
		references r(D);

		//#### Create Stencil Points (excl. center point), sorted, such that:
		// [ points(i) , points(points.size()-i-1) ] is
		// the i-th direction
		
		try{

			r.stencil_points.construct(D,args.first.stencil_type);

		}catch(boost::exception& e){
			e<<boost::errinfo_file_name(("Stencil Type "+toS(args.first.stencil_type)).c_str());
			throw;
		}
	
		uint proper_neighs=r.stencil_points.proper_neighs();

		//#### Allocate result memory (+1 for center weights)
		auto end_ind = result->A_io.resize(*grid,r.stencil_points.nDirs());

		//#### Construct the multi_indices corresponding to partial derivatives
		uint numpds = r.generatePDindices(args.first.partial_derivative_order);
		
		thread_logger()<<"WEIGHTS| for "+toS(proper_neighs)+" neighbors and "+toS(numpds)+" of PDE terms"<<endl;
		
		if(0){cout<<"X"<<endl;
		cout<<"WEIGHTS| for "+toS(proper_neighs)+" neighbors and "+toS(numpds)+" of PDE terms"<<endl;
		cout<<args.first<<endl;
		}
		//#### Check, if weight (unconstrained) determination is not overdetermined: 
		if (numpds == proper_neighs) 
			thread_logger()<<"WEIGHTS| (unconstrained weights) are NOT overdetermined"<<endl;
		else if (numpds < proper_neighs)
			thread_logger()<<"WEIGHTS| (unconstrained weights) are UNDERDETERMINED and not unique!"<<endl;

		//####  Compute line wise constant factorials in denominator 
		r.generateFactorialDenoms(args.first.higher_order_weight);

		//###### Optimize weights in parallel ( minimze norm: | D*w - b |) 
		{	using namespace tbb;
			optimize_body<fd_weights_result_t::A_io_t,body_t<PDE_core>> body
				(r,result->A_io,pde,*grid);
			thread_logger()<<"WEIGHTS| nSites: "<<end_ind.ind<<endl;
			tick_count t0 = tick_count::now();
			parallel_reduce(
				blocked_range<grid_index<sg>>(grid_index<sg>(),end_ind,100),
				body/*,simple_partitioner()*/);
			tick_count t1 = tick_count::now();
			thread_logger()<<"WEIGHTS| DONE, maxNorm:"<<scientific <<setprecision(3)<<
				body.max_norm <<" in " <<(t1-t0).seconds()<<"s"<<endl;
		}

		//###### create information for decoupled tridiagonal systems in each splitting direction
		//###### AND A_o
		tick_count t0 = tick_count::now();

		auto b_end = grid->b_end<sg>();

		auto& gsi = *new grid_striding_info_collection<sg>(r.stencil_points.nDirs(),r.stencil_points,grid);
		thread_logger()<< "&gsi: " << &gsi << endl;
		result->grid_striding_infos.reset(&gsi);		

		auto& A_io=result->A_io;
		auto& A_o=*result->A_o;
		A_o.resize(*grid);

		//test for uniqueness of the systems:
		//map<pair<int,int>,int> m;

		A_o->reserve(VectorXi::Constant(grid->Nbound[J_SGS],proper_neighs));

		//iterate over boundary and check for system beginnings
		for(connected_set_setup_iterator<sg> it(*grid);it!=b_end;++it){

			//cout<<"at "<<it.inner_start_g_ind().ind<<endl;

			auto i = r.stencil_points.dir_it();
			for(;i.valid();++i){//iterate through all directions
				auto& dir = r.stencil_points[i];
				if(it.apply_dir(dir)){//if a system in this direction exists
					
					//cout<<"dir: "<<dir.transpose()<<endl;

					/*
					//test for uniqueness of the system:
					//The naive test, that identifies a system
					//via its order start and end point on the inner boundary
					//fails in D>2, as a corner (of the inner boundary) can 
					//belong to more than one system (in D>2)!
					//int a=it.inner_start_b_ind().ind;int b=it.inner_end_b_ind.ind;
					//use this instead:
					int a=it.outer_start_b_ind().ind;int b=it.outer_end_b_ind().ind;
					if(a==2 && b==2){
						cout<<"pos: "<<it.it.pos.transpose()<<endl;
						cout<<"dir: "<<dir.transpose()<<endl;
					}
					m[make_pair(min(a,b),max(a,b))]++;
					if(m[make_pair(min(a,b),max(a,b))]>1){
						cout<<m[make_pair(min(a,b),max(a,b))]<<" Systems found between "<<min(a,b)<<" and "<<max(a,b)<<endl;
					}
					*/

					//add subsystem
					gsi.at(i).pushback_connected_set(it);
					
					bool debug = false;///\todo: set false

					//add A_o corresponding to start of system 
					A_o.insert(it.inner_start_b_ind(),it.outer_start_b_ind,
								(debug?1:A_io.at_opposite(i)[it.inner_start_g_ind()]));
					
					//cout<<it.inner_start_g_ind().ind<<": "<<A_io.at_opposite(i)[it.inner_start_g_ind()]<<endl;
					//cout<<it.inner_end_g_ind.ind<<": "<<A_io.at(i)[it.inner_end_g_ind]<<endl;

					//add A_o corresponding to end of system 
					A_o.insert(it.inner_end_b_ind,it.outer_end_b_ind,
								(debug?1:A_io.at(i)[it.inner_end_g_ind]));

					///\todo: check for natural boundaries (i.e. all A_o=0 on bound)
				}
			}
		}

		A_o->makeCompressed();

		if(global_config::get().test_A_o){
			test_A_o(A_o,A_io,gsi,r.stencil_points);
		}

		//A_o.print();

		tick_count t1 = tick_count::now();
		thread_logger()<<"WEIGHTS| DONE systems information in "<<scientific <<setprecision(3)
			<<(t1-t0).seconds()<<"s"<<endl;

		return result;
	}



}}
