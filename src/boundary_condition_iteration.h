#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <set>
#endif

#include "boundary_conditions_config.h"
#include "grid_iterators.h"
#include "grid.h"
#include "global_config.h"
#include "exceptions.h"

namespace fipster { namespace _boundary_conditions {


	typedef pair<const config_t, grid_ptr> bc_arg_t;

	/** this class can be derived to define actions that should be taken for 
		single every line (equation) of the boundary conditions.
	*/
	template<int sg_tpl>
	struct boundary_conditions_line{
		const static int sg=sg_tpl;
		virtual void add_to_interior(	boundary_index<sg>& row,
										grid_index<sg+1>& col,
										double v,
										uint variable,
										bstate_t state)=0;
		virtual void add_to_boundary(	boundary_index<sg>& row,
										boundary_index<sg>& col,
										double v,
										uint variable,
										bstate_t state)=0;
		virtual void set_zero(	boundary_index<sg>& row,
										uint variable,
										bstate_t state)=0;
	};

		template<int sg,type_t(bc::*bc_type_access)>
		struct bcl_test;

	/** an iterator to iterate through all boundary points and process the
		boundary condition accessed through bc_type_access and calling
		(virtual) members on a given boundary_conditions_line object
		*/
	template<int sg,type_t(bc::*bc_type_access)>
	struct boundary_condition_iterator{
		
		typedef boundary_iterator<sg,true> boundary_iterator_t;

		const config_t& config;
		boundary_iterator_t it;
		const boundary_index<sg> end;
		position_t temp;
		const Grid& g;

		boundary_condition_iterator(const bc_arg_t& arg, bool test=true)
			:	config(arg.first),
				it(*arg.second),
				end(arg.second->b_end<sg>()),
				g(*arg.second)
		{
			if(test && global_config::get().test_boundary_condition_iterator)
				bcl_test<sg,bc_type_access>().perform_test(arg);
		}

		void handle_bc_line(boundary_conditions_line<sg>& bcl){
			const auto ab = it.get_active_boundary();
			auto& bc=config.sv_vec[ab.first].lowerupper[ab.second];
			double inward_inverse_stepsize;
			switch(bc.*bc_type_access){
				case vonNeumann:
					//find the "other" point for the derivative
					temp = it.pos();

					if(ab.second==upper_ind){
						temp[ab.first]--;
						inward_inverse_stepsize = -1./g.stepsizes[ab.first][temp[ab.first]];
					}else{
						temp[ab.first]++;
						inward_inverse_stepsize = 1./g.stepsizes[ab.first][sg];
					}

					//add to diagonal
					bcl.add_to_boundary(it.b_ind,it.b_ind,-inward_inverse_stepsize,ab.first,ab.second);
					

					//if "other" point does not lie on the boundary
					if(it.n_active==1){
						grid_index<sg+1> op(temp,g);
						bcl.add_to_interior(it.b_ind,op,inward_inverse_stepsize,ab.first,ab.second);						
					}else{//if "other" point does lie on the boundary
						boundary_index<sg> op(temp,g);						
						bcl.add_to_boundary(it.b_ind,op,inward_inverse_stepsize,ab.first,ab.second);
					}

					break;

				case Dirichlet:
					//add to diagonal
					bcl.add_to_boundary(it.b_ind,it.b_ind,1.0,ab.first,ab.second);
					break;
				case zero_value:
					bcl.set_zero(it.b_ind,ab.first,ab.second);
			}
		}

		bool valid() const{
			return it!=end;
		}

		void operator++(){
			++it;
		}
	};

	/** a class to test the boundary_conditions_line and 
		boundary_condition_iterator functionality
		*/ 
	template<int sg,type_t(bc::*bc_type_access)>
	struct bcl_test : boundary_conditions_line<sg>{

		vector<boundary_index<sg>> rows;
		vector<tuple<boundary_index<sg>,grid_index<sg+1>,double>> interiors;
		vector<tuple<boundary_index<sg>,boundary_index<sg>,double>> boundaries;
		uint variable2;
		bstate_t state2;
		bool new_row;

		void update_state_variable(uint variable,bstate_t state){
			if(!new_row)
				FIPSTER_ASSERT(variable2==variable && state2==state);
			else{
				variable2=variable;
				state2=state;
				new_row=false;
			}
		}

		void perform_test(const bc_arg_t& arg){
			auto& g=*arg.second;
			uint count=0;
			for(auto bcit=boundary_condition_iterator<sg,bc_type_access>(arg,false);
				bcit.valid(); ++bcit){
				count++;
				new_row=true;
				interiors.clear();
				boundaries.clear();
				bcit.handle_bc_line(*this);
				switch(arg.first.sv_vec[variable2].lowerupper[state2].*bc_type_access){
				case vonNeumann:{
					FIPSTER_ASSERT(interiors.size()+boundaries.size()==2);
					double step;
					if(state2==upper_ind){							
						step=-g.stepsizes[variable2][g.sizes[variable2]-2];
					}else{
						step=g.stepsizes[variable2][sg];
					}
					if(interiors.size()==1)
					{
						FIPSTER_ASSERT(get<0>(boundaries[0])==get<1>(boundaries[0]));
						FIPSTER_ASSERT(get<2>(boundaries[0])==-1*get<2>(interiors[0]));
						FIPSTER_ASSERT(get<0>(boundaries[0])==get<0>(interiors[0]));
						auto pos=get<1>(interiors[0]).get_position(g);
						pos[variable2]+=(state2==upper_ind?1:-1);
						FIPSTER_ASSERT(get<0>(boundaries[0])==boundary_index<sg>(pos,g));
						FIPSTER_ASSERT(get<2>(interiors[0])==1/step);
					}else{//both points lie on boundary
						FIPSTER_ASSERT(get<2>(boundaries[0])==-1*get<2>(boundaries[1]));
						bool found=false;
						for(uint i=0;i<2;i++)
							if(get<0>(boundaries[i])==get<1>(boundaries[i]))
							{
								found=true;
								FIPSTER_ASSERT(get<1>(boundaries[i]).ind==get<1>(boundaries[1-i]).ind+
											(state2==upper_ind?1:-1)*
											g.strides.coeff(sg,variable2));
								FIPSTER_ASSERT(get<2>(boundaries[1-i])==1/step);
								break;
							}
						FIPSTER_ASSERT(found);

					}
				}
					break;
				case Dirichlet:
					FIPSTER_ASSERT(interiors.size()==0);
					FIPSTER_ASSERT(boundaries.size()==1);
					FIPSTER_ASSERT(get<0>(boundaries[0])==get<1>(boundaries[0]));
					FIPSTER_ASSERT(get<2>(boundaries[0])==1.0);
					break;
				case zero_value:
					FIPSTER_ASSERT(interiors.size()==0);
					FIPSTER_ASSERT(boundaries.size()==0);
					break;
				}
			}

			FIPSTER_ASSERT(arg.second->Nbound[sg]==count);
			FIPSTER_ASSERT(arg.second->Nbound[sg]==rearrange_and_get_number_of_uniques(rows));
		}
		
		void add_to_boundary(	boundary_index<sg>& row,boundary_index<sg>& col,
								double v,uint variable,bstate_t state){
			rows.push_back(row);
			boundaries.push_back(make_tuple(row,col,v));
			update_state_variable(variable,state);
		}
		void add_to_interior(	boundary_index<sg>& row,grid_index<sg+1>& col,double v,
								uint variable,bstate_t state){
			rows.push_back(row);
			interiors.push_back(make_tuple(row,col,v));
			update_state_variable(variable,state);
		}
		void set_zero(	boundary_index<sg>& row,uint variable,bstate_t state){
			rows.push_back(row);
			update_state_variable(variable,state);
		}
	};
}

	using _boundary_conditions::boundary_conditions_line;
	using _boundary_conditions::boundary_condition_iterator;
}
