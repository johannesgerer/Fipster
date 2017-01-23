#pragma once
#ifndef AUTOMATIC_PRECOMPILATION
#include <vector>
#endif

#include "dir_container.h"
#include "grid.h"

namespace fipster { 

	using namespace std;
	
	namespace _grid_striding_info  {
	

	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################


/** \brief is used to compile information about connected sets.

	This class contains a boundary iterator and generates information 
	about the connected set starting from the corresponding boundary 
	point in a certain direction.

*/
template<int sg_tpl>
struct connected_set_setup_iterator{

	const static int sg = sg_tpl;

	boundary_iterator<sg> it;
	boundary_index<sg>	inner_end_b_ind;
	boundary_index<sg-1> outer_start_b_ind,outer_end_b_ind;
	grid_index<sg>		inner_end_g_ind;
	position_t inner_end_pos,outer_end_pos,outer_start_pos;
	int length;
	grid_ref g;
	uint D;

	const boundary_index<sg>& inner_start_b_ind()const{
		return it.b_ind;
	}
	const grid_index<sg>& inner_start_g_ind()const{
		return it.g_ind;
	}

	bool operator!=(const boundary_index<sg>& y)const{ return it!=y; }
	//ctor
	connected_set_setup_iterator(const Grid& g)
		:it(g),inner_end_pos(g.D),outer_end_pos(g.D)
		,outer_start_pos(g.D),g(g),D(g.D)
	{}

	/** \brief Extracts information about the connected set in direction "dir"
		to use this: include "grid_striding_info_inline.h"
	*/
	bool apply_dir(const position_t& dir);

	/**
		\brief iterate one step on the boundary
	*/
	void operator++(){ ++it; }
};


template<int sg>
struct connected_set_t{

		grid_index<sg> begin_gi;
		boundary_index<sg-1> begin_bi,end_bi;
		uint cumulative_offset; ///< Equals the sum of the lengths of all connected_sets up to and excluding this
		int length;


		connected_set_t(const connected_set_setup_iterator<sg>& it,uint cumulative_offset)
		:begin_gi(it.inner_start_g_ind())
		,begin_bi(it.outer_start_b_ind)
		,end_bi(it.outer_end_b_ind)
		,cumulative_offset(cumulative_offset)
		,length(it.length)
		{}
	};

//#####################################################################
/** \brief contains information about a strided view or storage of a grid.
	
	It essentially contains a vector<connected_set_t> containing info about
	sets of grid points that are directly connected via the given grid_stride.
*/
template<int sg_tpl>
struct grid_striding_info_t{

	const static int sg = sg_tpl;


//#####################################################################
//types


	
	typedef vector<connected_set_t<sg>>  vec_t;
	typedef typename vec_t::const_iterator it_base_t;
	
	
//#####################################################################
//data
	vec_t connected_sets;
	grid_stride<sg> stride;
	int max_length;///< Equals the maximum of the lengths of all connected point sets
	int cumulative_length;///< Equals the sum of the lengths of all connected point sets


//#####################################################################
//members
	void pushback_connected_set(const connected_set_setup_iterator<sg>& it){
		connected_sets.emplace_back(connected_set_t<sg>(it,cumulative_length));
		max_length=max(max_length, it.length);
		cumulative_length += it.length;
	}

	//ctor
	void init(size_t reserve, grid_stride<sg> stride_)
	{ 
		cumulative_length=0;
		connected_sets.reserve(reserve);	
		max_length = 0;
		stride = stride_;
	}

	struct iterator : private it_base_t {
		int stride;

		using it_base_t::operator==;
		using it_base_t::operator++;

		iterator(int stride,it_base_t it)
			: stride(stride), it_base_t(it)
		{}
	};

	/*iterator begin(){
		return iterator(stride,std::begin(connected_sets));
	}

	it_base_t end(){
		return std::end(connected_sets);
	}*/


};



template<int sg>
struct grid_striding_info_collection 
	: dir_container<grid_striding_info_t<sg>,_dir::simple> {

	using dir_cont_t = typename grid_striding_info_collection::dir_cont_t;
	using dir_cont_t::at;

	typedef pair<grid_stride<sg>,neighbor_iterator> stride_map;
	
	struct strides_map_compare{
		bool operator()(const stride_map& a,const stride_map& b){
		return a.first.s<b.first.s;
	}};


	//to use this: #include "grid_striding_info_inline.h"
	template<class T>
	grid_striding_info_collection(const uint nDirs,const T& stencil_points,grid_ptr grid);

	grid_ptr grid;

	/** create a vector containing the strides and neighbor_iterators
		(including center) sorted by stride in ascending order
	*/
	vector<stride_map> sorted_neighs;
};


}

using _grid_striding_info::grid_striding_info_collection;
using _grid_striding_info::connected_set_setup_iterator;
using _grid_striding_info::grid_striding_info_t;
using _grid_striding_info::connected_set_t;


}
