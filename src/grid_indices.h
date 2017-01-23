#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <iterator>
#include <vector>
#include <stack>
#include <boost/throw_exception.hpp>
#include <comphash.h>
#endif


#include "Eigen.h"
#include "utility.h"
#include "exceptions.h"

namespace fipster { 
	
	using namespace std;

	//Standard grid vector types
	//Their (int n) c'tors do not need to default(zero) initialize!
	//Their () c'tors have to empty initialize! \todo ??? was heiﬂt das???
	//Much code relies on this
	typedef Eigen::VectorXi position_t;
	//typedef Eigen::VectorXi position_t;
	typedef Eigen::VectorXd state_t;
	
	namespace _grid {


	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################

	class Grid;

		/** class representing strides, i.e. offsets in a linearized
				representation of a grid corresponding to a vector
				(i.e. connection between two points) on a that grid */
	template<int sg,class Grid_t=Grid>
	struct grid_stride{

		// static const int sg=sg;

		int s;///\todo: make private

		//ctor
		grid_stride():s(0){}
		grid_stride(const position_t& dir,const Grid_t& g)
			:s(g.strides.row(sg)*dir){}

		//copy assignment
		grid_stride<sg,Grid_t>& operator=(const grid_stride<sg,Grid_t>& that)
		{
			s=that.s;
			return *this;
		}

		grid_stride<sg> operator*(int i) const{
			return grid_stride<sg>(s*i);
		}

	/*	conversion:
	operator int()const{
			return s;
		}*/

	private:
		//ctor
		grid_stride(int s):s(s){}
	};

	//##########################################################################
	//represents an index into the flattened D-dimensional sg-subgrid
	template<int sg,class Grid_t=Grid>
	struct grid_index :
	                      comphashable<grid_index<sg,Grid_t> >
		, stringifyable<grid_index<sg,Grid_t> >
{

		typedef grid_index<sg,Grid> my_t;

		// static const int sg=sg;
		int ind;///\todo: make private

		int& operator*(){
			return ind;
		}

		//ctor
		grid_index():ind(0){}


		//ctor
		//get index of a position
		grid_index(const position_t& pos,const Grid_t& g)
			:ind(-g.origin_offset(sg)+grid_stride<sg>(pos,g).s){
		}

		//ctor		
		friend Grid_t;
	private: //public would defy the whole idea of this class: 
		grid_index(int i):ind(i){}

	public:

		position_t get_position(const Grid_t& g) const{
			int ind2 = ind;
			uint i; position_t pos(g.D);
			for (i = 0; i < g.D-1; ++i) {
				int t = g.sizes[i]-sg*2;
				pos[i] = ind2 % t + sg;
				ind2 /= t;
			}
			pos[i] = ind2 + sg;
			return pos;
		}

		template<class T> FORCE_INLINE void my_combine(T& t) const{
			t(&my_t::ind);
		}

		FORCE_INLINE my_t operator+(const grid_stride<sg> k)const{
			return my_t(ind+k.s);
		}
		FORCE_INLINE my_t operator-(const grid_stride<sg> k)const{
			return my_t(ind-k.s);
		}
		FORCE_INLINE my_t& operator+=(const grid_stride<sg> k){
			ind+=k.s;
			return *this;
		}

		FORCE_INLINE my_t operator+(size_t k)const{
			return my_t(ind+k);
		}
	
		FORCE_INLINE my_t operator-(size_t k)const{
			return my_t(ind-k);
		}

		FORCE_INLINE int operator-(const my_t& y)const{
			return ind-y.ind;
		}

		//prefix increment operator
		FORCE_INLINE my_t& operator++(){
			ind++;
			return *this;
		}

		//prefix decrement operator
		FORCE_INLINE my_t& operator--(){
			ind--;
			return *this;
		}

		grid_index<sg> static end(const Grid_t& g){
			return g.template end<sg>();
		}

		int static total(const Grid_t& g){
			return g.N[sg];
		}

		my_t& operator=(const my_t& y){
			ind = y.ind;
			return *this;
		}

	};

	
	
	//##########################################################################
	//represents an index into the flattened boundary points of the D-dimensional sg-subgrid
	template<int sg,class Grid_t=Grid>
	struct boundary_index : comparable<boundary_index<sg,Grid_t> >{
		
		typedef boundary_index<sg,Grid> my_t;
		
		int ind;///\todo: make private
		// static const int sg=sg;

		int& operator*(){
			return ind;
		}

		int operator*() const{
			return ind;
		}

		//ctor
		boundary_index():ind(0){}

		///** create from pos */
		//template<bool check_for_boundary>
		//boundary_index<sg,Grid_t> create(const position_t& pos, const Grid_t& g)
		//{
		//	boundary_index<sg,Grid_t> ret;
		//	ret.make_ind_from_pos<check_for_boundary>(pos,g,grid_index<sg>(pos,g)))
		//	return ret;
		//}

		///** ctor from pos (and also export grid_index) */
		//template<bool check_for_boundary>
		//boundary_index<sg,Grid_t> create(const position_t& pos, const Grid_t& g,grid_index<sg>& g_ind)
		//{
		//	boundary_index<sg,Grid_t> ret;
		//	g_ind=grid_index<sg>(pos,g);
		//	ret.make_ind_from_pos<!FIPSTER_NDEBUG>(pos,g,g_ind);
		//		//BOOST_THROW_EXCEPTION(runtime_error("get_boundary_array_offset: pos lies outside of boundary!"));
		//}

		/** ctor from pos */
		boundary_index(const position_t& pos, const Grid_t& g)
		{
			if(!make_ind_from_pos<!FIPSTER_NDEBUG>(pos,g))
				FIPSTER_THROW_EXCEPTION(runtime_error("get_boundary_array_offset: pos lies outside of boundary!"));
		}

		/** ctor from pos with precalculated g_ind */
		boundary_index(const position_t& pos, const Grid_t& g,const grid_index<sg>& g_ind)
		{			
			if(!make_ind_from_pos<!FIPSTER_NDEBUG>(pos,g,g_ind))				
				FIPSTER_THROW_EXCEPTION(runtime_error("get_boundary_array_offset: pos lies outside of boundary!"));
		}

		

//		So oder so ‰hnlich (siehe Scan vom 14.2.12)
//
//		static position_t get_position(const Grid_t& g,uint ind) {
//			position_t pos(g.D);
//			get_position_from_ind2(g,g.D-1,ind,pos);
//			return pos;
//		}
//		
//		static void get_position_from_ind2(const Grid_t& g,uint j,uint ind2,position_t& pos){
//			if(j==0){
//				pos(0)=
//				cout<<pos.transpose()<<endl;
//			}else if(ind2<(uint)g.strides(sg,j)){ 
//				get_position_partial(g,j,ind2,pos);
//				cout<<pos.transpose()<<endl;
//			}else if(ind2>g.Nbound[sg]-g.strides(sg,j)){
//				ind2+=g.strides(sg+1,j)*(g.sizes[j]-sg-1-1-sg);
//				get_position_partial(g,j,ind2,pos);
//				cout<<pos.transpose()<<endl;
//			}else{
//				int diff = g.strides(sg,j)-g.strides(sg+1,j);
//				ind2-=g.strides(sg,j);
//				pos(j)=ind2/diff+sg+1;
//				ind2-=(pos(j)-sg-1)*diff;
//				get_position_from_ind2(g,j-1,ind2,pos);
//			}
//		}
///** fills pos up to (and including) index j */
//		static void get_position_partial(const Grid_t& g,uint j,uint ind2,position_t& pos) {
//			for (uint i = 0; i < j; ++i) {
//				int t = g.sizes[i]-sg*2;
//				pos[i] = ind2 % t + sg;
//				ind2 /= t;
//			}
//			pos[j] = ind2 + sg;
//		}


		template<bool check_for_boundary>
		bool make_ind_from_pos(const position_t& pos, const Grid_t& g)
		{
			return make_ind_from_pos<check_for_boundary>(pos,g,grid_index<sg>(pos,g));
		}
		
		/** extracts the boundary index from a given position (vector)
		\param pos position vector on the grid
		\param g Grid
		\param ind has to be set to grid_index<sg>(pos,g).ind;
		\tparam check_for_boundary indicates, whether to check if the points actually lies on boundary
		*/
		template<bool check_for_boundary>
		bool make_ind_from_pos(const position_t& pos, const Grid_t& g,const grid_index<sg>& g_ind)
		{
			if(!FIPSTER_NDEBUG || check_for_boundary)
				FIPSTER_ASSERT(grid_index<sg>(pos,g)==g_ind);

			ind=g_ind.ind;
			int i;

			//this checks, if pos lies either on boundary or interior and not outside
			if(check_for_boundary)
				for(uint j=0;j<g.D;j++)
					if(pos(j)<sg || pos(j)>g.sizes[j]-1-sg)
						return false;

			//loop through every dimension until "break"
			for(i=g.D-1;i>=0;--i){

				//if the point sits on the lower bound, the internal 
				//offset for all following dimensions is easily calculated (zero)
				if (pos(i)==sg){
					//if(with_active_bound) active_bound = 2*i;
					break;
				}else{
					//if not, reduce the offset by a value corresponding to a shift to the 
					//lower boundary
					ind -= g.strides(sg+1,i)*(pos(i)-1-sg);

					//if the point sits on the upper bound, 
					if(pos(i)==g.sizes[i]-sg-1){
						//if(with_active_bound) active_bound = 2*i+1;
						break;
					}
				}	
			}

			//this checks, that pos does not lie on interior
			if(check_for_boundary && i==-1)
				return false;
			
			return true;
		}

		//ctor		
		friend Grid_t;
	private: //public would defy the whole idea of this class: 
		boundary_index(int i):ind(i){}
	public:

		template<class T> FORCE_INLINE void my_combine(T& t) const{
			t(&my_t::ind);
		}	

		
		boundary_index<sg> static end(const Grid_t& g){
			return g.template b_end<sg>();
		}

		int static total(const Grid_t& g){
			return g.Nbound[sg];
		}

		//prefix increment operator
		FORCE_INLINE my_t& operator++(){
			ind++;
			return *this;
		}
	};
	
}

	using _grid::boundary_index;
	using _grid::grid_index;
	using _grid::grid_stride;
}
