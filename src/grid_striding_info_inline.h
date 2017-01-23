#pragma once
#ifndef AUTOMATIC_PRECOMPILATION
#include <vector>
#include <utility>
#endif

#include "grid_striding_info.h"

namespace fipster { 

	using namespace std;

	namespace _grid_striding_info  {


		// #################################################
		// ##############   PHASE 3:			############
		// ##############   Node Bodies         ############
		// #################################################

template<int sg>
bool connected_set_setup_iterator<sg>::apply_dir(const position_t& dir)
{
	outer_start_pos = it.pos() - dir;

	//check if points lies on outer boundary
	if(outer_start_b_ind.template make_ind_from_pos<true>(outer_start_pos,g))
	{
		//if it does, update all relevant information:			
		length = g.get_system_length<sg,false>(dir,it);
		inner_end_pos = it.pos()+(length-1)*dir;
		outer_end_pos = inner_end_pos +dir;
		inner_end_g_ind = grid_index<sg>(inner_end_pos,g);
		inner_end_b_ind = boundary_index<sg>(inner_end_pos,g,inner_end_g_ind );
		outer_end_b_ind = boundary_index<sg-1>(outer_end_pos,g);
		return true;
	}

	//return false, if point is no boundary point
	return false;
}



// #################################################
template<int sg>
template<class T>
	grid_striding_info_collection<sg>::grid_striding_info_collection
		(const uint nDirs,const T& stencil_points,grid_ptr grid)
	: dir_cont_t(nDirs),grid(grid)
{
	sorted_neighs.reserve(nDirs*2+1);

	//prepare directions
	auto i=this->dir_it();
	sorted_neighs.push_back(make_pair(grid_stride<J_SGS>(),i.at_center_it()));
	for(;i.valid();++i){//iterate through all directions
		
		int total=0;

		for(uint j=0;j<grid->D;j++){
			//estimate upper bound for number of subsystems
			if(stencil_points[i](j)!=0)//if the boundary layer normal to "j" is active
				total += grid->N(sg)/(grid->sizes[j]-sg*2);//estimate
		}

		//reserve enough space for the subsystems and set stride
		at(i).init(total,grid_stride<sg>(stencil_points[i],*grid));
		
		sorted_neighs.push_back(make_pair(at(i).stride,			i.at_it()));
		sorted_neighs.push_back(make_pair(at(i).stride*(-1),	i.at_opposite_it()));
	}
		
	FIPSTER_ASSERT(sorted_neighs.size()==nDirs*2+1);
	sort(begin(sorted_neighs),end(sorted_neighs),strides_map_compare());
}

}
}
