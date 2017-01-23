#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#endif

#include "logging.h"
#include "boundary_conditions.h"
#include "timestepping_operator.h"
#include "grid_fields.h"
#include "exceptions.h"
#include "global_config.h"
#include "Math.h"
#include "type_traits.h"

namespace fipster { namespace _timestepping_operator {
	

	using namespace Eigen;

	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################

	template<bool start, bool end>
	struct main_loop_body{
		template<class A, class B, class C, class D, class E,class F, class G,class H>
		///\todo, check if inlined
		FORCE_INLINE static void go(const A& nn_map,const H& bend,const B& i,
						const C& grid_ind, D& M,const E& A_io,
						const F& Ao_Bo1_Bi_dense,const G& dir, double oneOverDirs){
			/**\todo: this question only has to be asked on a system
			that lies on the (inner) boundary
			if there is a row in the A.Bo-1.Bi matrix*/
			auto b_ind=nn_map.at(grid_ind);
			if(b_ind!=bend)
				inner_body<true>::go(b_ind,i,grid_ind,M,A_io,Ao_Bo1_Bi_dense,dir,oneOverDirs);
			else
				inner_body<false>::go(b_ind,i,grid_ind,M,A_io,Ao_Bo1_Bi_dense,dir,oneOverDirs);
		}

		template<bool on_boundary>
		struct inner_body{
			template<class A, class B, class C, class D, class E,class F, class G>
		///\todo, check if inlined
		FORCE_INLINE static void go(const A& b_ind,const B& i,const C& grid_ind, 
						D& M,const E& A_io, const F& Ao_Bo1_Bi_dense,const G& dir,
						double oneOverDirs)
		{
			//copy super and sub diagonals. IMPORTANT: this is also done for
			//boundary points, thus yielding tridiagonal matrices plus
			//subdiag(0) and superdiag(n-1) this is used at CodeRefA

			double super = A_io.at(dir).at(grid_ind)
								-(on_boundary&&!end?
									Ao_Bo1_Bi_dense.at(dir).at(b_ind):0);

			double sub	= A_io.at_opposite(dir).at(grid_ind)
								-(on_boundary&&!start?
									Ao_Bo1_Bi_dense.at_opposite(dir).at(b_ind):0);
	
			
			M.subdiag_ref(i)	= sub; 
			M.superdiag_ref(i)	= super;

			M.diag_ref(i)=oneOverDirs*(A_io.at_center(dir).at(grid_ind)
								- (on_boundary?Ao_Bo1_Bi_dense.at_center(dir).at(b_ind):0)
								) - (!start?sub:0) - (!end?super:0);
			
			if(0 && !FIPSTER_NDEBUG){
#if 1
				cout<<"center: "<<M.diag(i)<<endl;
				cout<<A_io.at_center(dir).at(grid_ind)<<endl;
				cout<<(on_boundary?Ao_Bo1_Bi_dense.at_center(dir).at(b_ind):0)<<endl;
				cout<<"super"<<endl;
				cout<<A_io.at(dir).at(grid_ind)<<endl;
				cout<<(on_boundary&&!end?
									Ao_Bo1_Bi_dense.at(dir).at(b_ind):0)<<endl;
				cout<<"sub"<<endl;
				cout<<A_io.at_opposite(dir).at(grid_ind)<<endl;
				cout<<(on_boundary&&!start?
									Ao_Bo1_Bi_dense.at_opposite(dir).at(b_ind):0)<<endl;
				cout<<"offs"<<endl;
				cout<<super<<endl;
				cout<<sub<<endl;
#endif				
				M.template check_M_property<true,start,end>(i,false);
			}
		}
		};
	};

	/** This test checks, whether the result actually 
		contains what is specified in Scan (12-02-22 1).
	*/
	template<class A, class B, class C, class E,class F,class G>
	void perform_test(const string& name,	const A& result,const B& Ao_Bo1_Bi_sparse, const C& A_io, 
						const E& A_o,const F& grid,const G& gsi){

		thread_logger()<<"SPLITIN| TESTING Split"<<endl;

		int ndirs = result.nDirs();
		grid_field_dir_container<grid_index<J_SGS>,Matrix<int,Dynamic,Dynamic,ColMajor>> n_tests;
		n_tests.resize(grid,ndirs);
		n_tests.setConstant(0);
		grid_field<grid_index<J_SGS>,VectorXd> diag_weights;
		grid_field<grid_index<J_SGS>,VectorXd> offdiag_weights;
		offdiag_weights.resize(grid);
		offdiag_weights->setZero();
		diag_weights.resize(grid);
		diag_weights->setZero();

		

		//loop over directions:
		for(auto dir = result.dir_it();dir.valid();++dir){
			//loop over systems:
			for(auto sys=result.at(dir).begin();sys.valid();++sys){
				auto grid_ind = sys.connected_set_it->begin_gi;
				
				auto stride = gsi.at(dir).stride;
				auto length = sys.connected_set_it->length;
				int i=0;


				//system interior
				for(; i<length;	i++,grid_ind+=stride)
				{
					//check, if point lies on (inner) boundary:
					position_t pos=grid_ind.get_position(grid);
					boundary_index<J_SGS> b_ind;
					bool on_b=b_ind.make_ind_from_pos<true>(pos,grid,grid_ind);

					//accumulate all center entries of the result + all off-diagonal
					//entries of A

					// cout<<dir<<": "<<center_weights.at(grid_ind)<<"+"<<sys.M.diag(i)<<"+"
					// 		<<A_io.at_opposite(dir).at(grid_ind)<<"+"
					// 		<<A_io.at(dir).at(grid_ind);

					diag_weights.at(grid_ind)+=sys.M.diag(i);
					offdiag_weights.at(grid_ind)+=A_io.at_opposite(dir).at(grid_ind)
						+A_io.at(dir).at(grid_ind);

					// cout <<"="<<center_weights.at(grid_ind)<< "=="
					// 		 <<A_io.at_center(dir).at(grid_ind)<<endl;
					
					// //this only is meaningful if nDirs == 1
					// if(center_weights.at(grid_ind) A_io.at_center(dir).at(grid_ind)

					FIPSTER_ASSERT(sys.M.subdiag(i)==A_io.at_opposite(dir).at(grid_ind)
													- (on_b&&i!=0?
														Ao_Bo1_Bi_sparse.coeff(b_ind,grid_ind-stride):0));
					FIPSTER_ASSERT(sys.M.superdiag(i)==A_io.at(dir).at(grid_ind)
													- (on_b&& i != length-1?
														Ao_Bo1_Bi_sparse.coeff(b_ind,grid_ind+stride):0));

					n_tests.at(dir).at(grid_ind)+=ndirs;
					n_tests.at_opposite(dir).at(grid_ind)+=ndirs;
					n_tests.at_center(dir).at(grid_ind)++;
				}
			}
		}

		//check, if all grid points were processed:
		auto end=grid.template end<J_SGS>();
		auto dir = A_io.dir_it();
		double tolerance=global_config::get().test_splitting_tolerance;
		for(auto grid_ind=grid.template begin<J_SGS>();grid_ind!=end;++grid_ind){

			//check, if point lies on (inner) boundary:
			position_t pos=grid_ind.get_position(grid);
			boundary_index<J_SGS> b_ind;
			bool on_b=b_ind.make_ind_from_pos<true>(pos,grid,grid_ind);

			//check, if center sum equals the unsplit diagonal elements.
			//in analogy to the Scan this corresponds to:
			// Sum(split result diags) == A_diag - AoBo1Bi_diag - Sum(offdiags A) 

			double rd=reldiff(offdiag_weights.at(grid_ind),
					A_io.at_center(dir).at(grid_ind)
												- (on_b?Ao_Bo1_Bi_sparse.coeff(b_ind,grid_ind):0)
												- diag_weights.at(grid_ind)
												,true);
			if(rd>=tolerance){
				cout<<scientific<<setprecision(17);
				cout<<"grid_ind\n"<<grid_ind.ind<<endl;
				cout<<"offdiag_weights.at(grid_ind): "<<endl<<offdiag_weights.at(grid_ind)<<endl;
				cout<<"diag_weights.at(grid_ind): "<<endl<<diag_weights.at(grid_ind)<<endl;
				cout<<"A_io.at_center(dir).at(grid_ind): "<<endl<<A_io.at_center(dir).at(grid_ind)<<endl;
				cout<<"(on_b?Ao_Bo1_Bi_sparse.coeff(b_ind,grid_ind):0): "<<endl<<(on_b?Ao_Bo1_Bi_sparse.coeff(b_ind,grid_ind):0)<<endl;
				cout<<__FILE__<<":"<<__LINE__<<": error: "<<rd<<endl;
				cout<<name<<endl;
				// FIPSTER_ASSERT(rd<tolerance);
			}



			//check all neighbors
			for(auto j=n_tests.neigh_it();j.valid();++j)
				if((uint)n_tests.at(j).at(grid_ind)!=result.nDirs()){
					cout<<grid_ind.ind<<": "<<n_tests.at(j).at(grid_ind)<<endl;
					FIPSTER_THROW_EXCEPTION(runtime_error("Problem in splitting's perform_test (1)"));
				}
		}

 		thread_logger()<<"SPLITIN| DONE TESTING Split"<<endl;
	}

	/** The splitter has to perform three steps (not necessarily in this order!)
	
		1. Get the endomorphism corresponding to an elimination of the boundary variables:
			A_i - A_o B_o^(-1) B_i  (See scan 11-10-21 BCs fundamentals and specifics.jpg)
			
		2. The center weights have to be "finished"
			(NOTE: See "finite_difference_weights::body_t" for explanation)

		3. Split this operator:
			Every pair of opposite stencil points together with the center 
			is used to create a tridiagonal operator. The center weight is split 
			equally among the directions.

	*/
	splitter_body_t::result_sptr_t 
		splitter_body_t::operator()(const input_t& futures){
			thread_logger()<<"SPLITIN| START"<<endl;

			//aliases
			const B_ptr_t&				B_ptr	=  get<0>(futures);
			const fd_weights_result_t&	weights	= *get<1>(futures);

			//create result:
			auto result = make_shared<tridiagonally_split_grid_operator>();

			//aliases
			const auto& gsi = *weights.grid_striding_infos;
			grid_ref grid = *gsi.grid;
			const auto& A_io = weights.A_io;
			const auto& A_o = *weights.A_o;

			result->init(weights.grid_striding_infos);
			auto& op = *result;


		//1. calculate A_o*Bi
			
			{
				using namespace _boundary_conditions;
				static_assert( is_same<A_o_t::col_index_t,Bo_t::row_index_t>::value,
											"is_same<A_o_t::col_index_t,Bo_t::col_index_t>::value");
				static_assert( is_same<Bi_t::row_index_t,Bo_t::col_index_t>::value,
											"is_same<Bi_t::row_index_t,Bo_t::col_index_t>::value");
			}

			//store temporarily in sparse matrix
			typedef sparse_field<boundary_index<J_SGS>,grid_index<J_SGS>,RowMajor> 
				Ao_Bo1_Bi_sparse_t;

			Ao_Bo1_Bi_sparse_t Ao_Bo1_Bi_sparse;
			*Ao_Bo1_Bi_sparse= **weights.A_o**B_ptr->Bo1**B_ptr->Bi; ///\todo TEST : is the EIGEN product correct?

			///\todo Performance: Ao_Bo1_Bi_sparse is possibly very (row-)sparse. This could be exploited?

		//loop through the (inner) boundary and copy the values into a dense matrix:
			grid_field<grid_index<J_SGS>,Matrix<boundary_index<J_SGS>,Dynamic,1>> nn_map;
			nn_map.resize(grid);
			auto bend=grid.b_end<J_SGS>();
			nn_map->setConstant(bend);

			grid_field_dir_container<boundary_index<J_SGS>> Ao_Bo1_Bi_dense;
			auto end_b = Ao_Bo1_Bi_dense.resize(grid,A_io.nDirs());
			Ao_Bo1_Bi_dense.setZero();

			bool do_throw=true;

			auto dir_it=Ao_Bo1_Bi_dense.dir_it();

			//for every boundary_index (row) find the direction 
			//corresponding to the non-zeros in the sparse array.
			//Every non-zero has to correspond to a direction
			//(while the opposite is not true)
			for(auto bit = boundary_iterator<J_SGS>(grid); bit != end_b; ++bit)
			{
				double row_sum=0;

				//get non-zero iterator
				Ao_Bo1_Bi_sparse_t::sm_t::InnerIterator 
					nonzero_it(*Ao_Bo1_Bi_sparse,bit.b_ind.ind);

				bool no_entires=!nonzero_it;

				//loop over directions
				const auto s_end=end(gsi.sorted_neighs);
				for(auto s=begin(gsi.sorted_neighs);s!=s_end && nonzero_it;++s){
					grid_index<J_SGS> neighbor_ind = bit.g_ind+s->first;

					//if non-zero entry's index matches
					if(nonzero_it.col()==neighbor_ind.ind){
						Ao_Bo1_Bi_dense.at(s->second).at(bit.b_ind)=nonzero_it.value();
						row_sum+=nonzero_it.value();
						++nonzero_it;
						continue;

					}else if(nonzero_it.col()<neighbor_ind.ind){
						// int tt=nonzero_it.col();
						if(do_throw) BOOST_THROW_EXCEPTION(runtime_error(
						"Boundary Conditions not compatible with FD stencil (1)!"));
					}else{
						/*auto neighbor_ind2 = 
						position_t pos=neighbor_ind.get_position(grid);
						cout<<pos<<endl;
						boundary_index<J_SGS-1> b_ind;
						if(b_ind.make_ind_from_pos<true>(pos,grid))
							FIPSTER_ASSERT(A_o.coeff(bit.b_ind,b_ind)==
							A_io.at(strides[i].second).at(bit.ind));*/
					}
				}
				if(nonzero_it)
					if(do_throw) BOOST_THROW_EXCEPTION(runtime_error(
						"Boundary Conditions not compatible with FD stencil (2)!"));
				
				//As concluded in Scan 22.2.2012/1+2:
				//Add the sum of the entries of the current row of A_o to the center entry:
				//get non-zero iterator
				finite_difference_weights::A_o_t::sm_t::InnerIterator
					it(*A_o,bit.b_ind.ind);
				no_entires = no_entires && !it;
				//sum up values
				while(it){
					
					row_sum += it.value();					
					++it;}

				Ao_Bo1_Bi_dense.at_center(dir_it).at(bit.b_ind)=row_sum;

				//if there are entries, register this boundary index into the nn_map
				if(!no_entires)
					nn_map.at(bit.g_ind)=bit.b_ind;
			}

			//B_ptr->Bo1.print();cout<<endl;
			//B_ptr->Bi.print();cout<<endl;
			//Ao_Bo1_Bi_sparse.print();cout<<endl;
			
		//2. fill the tridiagonal matrices with A_i - A_o B_o^(-1) B_i
			op.init(weights.grid_striding_infos);
			double oneOverDirs = 1.0/op.nDirs();

			do_throw=true;

			//loop over directions:
			for(auto dir = op.dir_it();dir.valid();++dir)
			{
				//loop over systems:
				for(auto sys=op.at(dir).begin();sys.valid();++sys){

					auto grid_ind = sys.connected_set_it->begin_gi;
					auto stride = gsi.at(dir).stride;
					auto length = sys.connected_set_it->length;

					int i=0;

					if(length==1)
						main_loop_body<true,true>::go(nn_map,bend,i,grid_ind,sys.M,A_io,Ao_Bo1_Bi_dense,dir,oneOverDirs);
					else{
						//system start
						main_loop_body<true,false>::go(nn_map,bend,i,grid_ind,sys.M,A_io,Ao_Bo1_Bi_dense,dir,oneOverDirs);

						//system interior
						for(i++,grid_ind+=stride; i<length-1;i++,grid_ind+=stride)
							main_loop_body<false,false>::go(nn_map,bend,i,grid_ind,sys.M,A_io,Ao_Bo1_Bi_dense,dir,oneOverDirs);

						//system end:
						main_loop_body<false,true>::go(nn_map,bend,i,grid_ind,sys.M,A_io,Ao_Bo1_Bi_dense,dir,oneOverDirs);
					}
					
				}//loop over systems:
			}//loop over directions:

			if(global_config::get().test_splitting)
				perform_test(fdw_conf.id,op,Ao_Bo1_Bi_sparse,A_io,A_o,grid,gsi);

			thread_logger()<<"SPLITIN| DONE"<<endl;
			return result;
	}
	
}}
