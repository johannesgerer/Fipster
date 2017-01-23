#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <iterator>
#include <vector>
#include <stack>
#include <boost/throw_exception.hpp>
#include <type_traits>
#include <comphash.h>
#endif

#include "integer_types.h"
#include "Eigen.h"
#include "grid_indices.h"
#include "boundary_enum.h"

namespace fipster { 

	using namespace std;

	namespace _grid {

	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################

	struct without_state{
		without_state(uint D){};
	};

	/** \brief represents a position on a sg-subgrid with corresponding state
	 */
	template<int sg,class Grid_t,bool with_state=false>
	struct position_state
	{

		typedef position_state<sg,Grid_t,with_state> position_state_t;
		typedef typename conditional<with_state,state_t,without_state>::type cond_state_t;

	private:
		position_t position;
		cond_state_t internal_state;
	public:
		
		// static const int sg=sg;
		
		const Grid_t& g;
		const uint D;

		position_state(const Grid_t& g,const position_t& pos)
			: position(pos),internal_state(g.D),g(g),D(g.D)
		{
			if(with_state)
				for(uint i=0;i<D;i++)
					update_state(i,internal_state,g,position);
			static_assert(sg<Grid_t::implemented_inner_grids,"sg>=implemented_inner_grids !!");
		}
		void change_pos(uint j,int p){
			position[j]=p;
			update_state(j,internal_state,g,position);
		}
		const position_t& pos() const{
			return position;
		}
		const state_t& state() const{
			return internal_state;
		}

		void check_bounds() const{
			if(with_state)
				for(uint i=0;i<D;i++)
					(g.coords.at(i)).at(position[i]);
		}
	};

	template<class Grid_t>
	void update_state(uint j,state_t& state,const Grid_t& g,const position_t& pos){
			state[j]= g.coords[j][pos[j]];
	}

	template<class Grid_t>
	void update_state(uint j,without_state& state,const Grid_t& g,const position_t& pos){};

		/**##########################################################################
	//position and state iterator template
	//it's templated on Grid_t, to avoid cyclic dependencies of 
	//class definitions with "Grid"
	slow usage: 
		for(grid_iterator<sg>i(grid);i.valid();++i)
	*/
	template<int sg,bool with_state=false,class Grid_t=Grid>
	struct grid_iterator : position_state<sg,Grid_t,with_state>
	                     , comphashable<grid_iterator<sg,with_state,Grid_t> >
	                     , stringifyable<grid_iterator<sg,with_state,Grid_t> >
	{
		typedef grid_iterator<sg,with_state,Grid> my_t;
		typedef typename grid_iterator::position_state position_state_t;

		using position_state_t::pos;
		using position_state_t::change_pos;
		using position_state_t::g;

		grid_index<sg,Grid_t> ind;

		bool on_leading_boundary()const{
			return pos()[0]==sg;
		}

		const grid_index<sg,Grid_t>& index() const{
			return ind;
		}

		grid_index<sg,Grid_t>& index(){
			return ind;
		}

		//initialized the position to point to
		//the origin of the sg-Grid
		grid_iterator(const Grid_t& g)
			: position_state_t(g,position_t::Constant(g.D,sg)),ind()
		{}

		////non-origin ctors:
		grid_iterator(const Grid_t& g,grid_index<sg,Grid_t> ind2)
			: position_state_t(g,ind2.get_position(g)),ind(ind2)
		{}
		grid_iterator(const Grid_t& g,const position_t& p)
			: position_state_t(g,p),ind(p,g)
		{}

		//prefix increment operator
		my_t& operator++(){

			++ind;

			for(uint j=0; j < this->D; ++j) 
				if (pos()[j] == g.sizes[j]-sg-1) 
					change_pos(j,sg);
				else {
					change_pos(j,pos()[j]+1);
					break;
				}

			return *this;
		}

		//this is a bad idea, as this enables a '<' operator which produces
		//wrong results?
		// operator bool()const{
		// 	return ind!=ind.end(g);
		// }
		bool valid()const{
			return ind!=ind.end(g);
		}
		
		bool operator==(const my_t& y)const{
			return operator==(y.ind);
		}

		bool operator==(grid_index<sg> y)const{
			return ind==y;
		}

		bool operator!=(grid_index<sg> y)const{
			return !operator==(y);
		}


		template<class T> FORCE_INLINE void my_combine(T& t)const{
			t(&grid_iterator::ind); 
			//compare only the ind, which uniquely defines a grid_iterator
			//on a given grid
		}
	};
	

	//##########################################################################
	template<int sg,bool active_boundary_set=false,bool with_state=false,class Grid_t=Grid>
	struct boundary_iterator : position_state<sg,Grid_t,with_state>{

		// static const bool active_boundary_set = active_boundary_set;

		typedef boundary_iterator<sg,active_boundary_set,with_state,Grid_t> my_t;
		typedef boundary_index<sg,Grid_t> b_ind_t;
		typedef pair<uint,bstate_t> active_boundary_t;

		typedef typename boundary_iterator::position_state position_state_t;
		using position_state_t::change_pos;
		using position_state_t::g;
		using position_state_t::pos;
		using position_state_t::D;

		vector<bstate_t> active_set;

		grid_index<sg,Grid_t> g_ind;
		b_ind_t b_ind;

		const b_ind_t& index() const {
			return b_ind;
		}

		b_ind_t& index(){
			return b_ind;
		}

		uint n_active;/* (number of (intersecting) border layers, the current ind is located on,
							ignoring the first dimension) + 1*/

		//initialized the position to point to
		//the origin of the sg-Grid
		boundary_iterator(const Grid_t& g) 
			: position_state_t(g,position_t::Constant(g.D,sg))
			,active_set(active_boundary_set?g.D:0,lower_ind)
			,g_ind(),b_ind(),n_active(g.D){}

	private:

		void overrun(uint j)
		{
			change_pos(j,sg);
			g_ind.ind += (sg-(g.sizes[j]-1-sg))*g.strides(sg,j);			
		}
		void jump_to_right(uint j)
		{
			change_pos(j, g.sizes[j]-1-sg);
			g_ind.ind += (pos()(j)-sg)*g.strides(sg,j);
		}
		void increment(uint j){
			g_ind.ind += g.strides(sg,j);
			change_pos(j,pos()(j)+1);
		}
		

		void change_upper_to_lower(uint j){
			if(active_boundary_set)
				active_set[j] = lower_ind;
		}
		void change_lower_to_upper(uint j){
			if(active_boundary_set)
				active_set[j] = upper_ind;
		}
		void deactivate(uint j){
			if(active_boundary_set)
				active_set[j] = inactive_ind;
			n_active--;
		}
		void activate_upper(uint j){
			if(active_boundary_set)
				active_set[j] = upper_ind;
			n_active++;
		}

	public:
		//prefix increment operator
		my_t& operator++(){
			++b_ind;
			
			//increment higher dimensions
			for (uint j = 0; j < D; ++j){

				//if pos sits on the right boundary, overrun
				if( pos()(j) == -2+g.sizes[j]-sg+1 ){
					change_upper_to_lower(j);
					overrun(j);
					continue;
				}

				//if pos() sits on the left boundary and is about to leave it,
				if(pos()(j) == sg){

					//first dimension, extra check:
					//if pos() sits on the left boundary and is about to leave it, 
					// and no other boundary is active
					if(j==0 && n_active == 1){
						change_lower_to_upper(0);
						jump_to_right(0);
						return *this; //true
					}else
						deactivate(j);

				}//else if pos() is about to enter the boundary it,
				else if(pos()(j) == -2+g.sizes[j] - sg){
					activate_upper(j);
				}

				increment(j);
				return *this; //true
			}

			//if every dimension has overrun, we are done
			return *this; //false
		}

		active_boundary_t get_active_boundary(){
			//FEATURE: the selection could be prioritized
			//			by defining the iteration order through the active boundary set 
			for(uint i=0;i<D;i++)
				if(active_set[i]!=inactive_ind)
					return make_pair(i,active_set[i]);
						
			BOOST_THROW_EXCEPTION(runtime_error("no active boundary found!"));
		}

		bool operator==(b_ind_t y)const{
			return b_ind==y;
		}

		bool operator!=(b_ind_t y)const{
			return !operator==(y);
		}
	
	};

	
}

	using _grid::grid_iterator;
	using _grid::boundary_iterator;
	using _grid::position_state;
}
