#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <boost/optional.hpp>
#include <type_traits>
#endif

#include "Eigen.h"
#include "state_variables_config.h"
#include "grid_iterators.h"

namespace fipster { namespace _grid {

	using namespace std;
	using namespace boost::property_tree;
	using namespace Eigen;

	//Represents one Grid direction. The data of double1 holds the coordinates
	struct state_variable_t
	{

		//(Total, including boundary) number of points in this direction.
		uint n;

		//Lower and upper bound of the corresponding state variable 
		// (the 0-th and the n-1-th point)
		double lower_bound,upper_bound;
		double tolerance;

		string spacing_type;
		const ptree& spacing_node;
		bool hits_zero;

		state_variable_t(const ptree& pt);

		template<class A>
		void fill_coords(A& coords, A& stepsizes);
	};

//Represents a regular rectangular D-dimensional grid.
//The outer grid has the in each directions i the dimensions:
//  0 ... sizes[i]-1
//Besides the outer grid, there are several inner grids. They are numbered, such that 
// the 0-grid is the outer grid, the 1-grid is the interior grid from 1 ... sizes[i]-2, 
// where one layer of boundary points is excluded.
//###############################################################################	
class Grid : configuration::state_variables_config<state_variable_t>  {
public:

	//this type is chosen for safety (bound checking for dimensions of unequal sizes)
	//performance did not change compared to continuous 2D array
	typedef vector<vector<double>> stepsizes_t;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	static const uint implemented_inner_grids = 3;

	//Dimension of the hypercube
	uint D;

	//Number of grid points (N[0]: total number, N[1]: interior points)
	Matrix<uint,implemented_inner_grids,1> N;

	//Number of grid points on each inner grid's inner boundary
	Matrix<uint,implemented_inner_grids-1,1> Nbound;

	//Dimensions of the hypercube.
	position_t sizes;
	int maxsize;
	
	string id;
	boost::optional<uint> max_subgrid;

	//coordinates:	
	//coords[i][j] is the coordinate of the j-th point in the i-th dimension
	//with j < sizes[i]
	stepsizes_t coords;
	//Step sizes: stepsize[i][j] = coords[i][j+1] - coords[i][j],
	//with j < sizes[i]-1
	stepsizes_t stepsizes;

	//strides(i,j) contains the stride for dimension j
	//on the i-th grid. 
	Matrix<int,implemented_inner_grids,Dynamic,RowMajor> strides;
	Array<int,implemented_inner_grids,1> origin_offset;

	//(const) Grid Functions #########################################################
	
	// Iterators #########################################################
	template<int sg>//specify the sub-grid to work on
	boundary_index<sg> b_end() const{ return boundary_index<sg>(Nbound[sg]); }
	template<int sg>//specify the sub-grid to work on
	grid_index<sg> end() const{ return grid_index<sg>(N[sg]); }

	template<int sg>//specify the sub-grid to work on
	boundary_index<sg> b_begin() const{ return boundary_index<sg>(0); }
	template<int sg>//specify the sub-grid to work on
	grid_index<sg> begin() const{ return grid_index<sg>(0); }

	//################################################################################
	//inline and template functions
	template<int sg, bool with_boundingdim>
	int  get_system_length( const position_t& dir,  const boundary_iterator<sg>& it, int& boundingdim=*(int*)nullptr)  const;

	// Constructor #########################################################
	Grid(const ptree& pt);

	//printing
	friend ostream& operator<<(ostream& out,const Grid& grid);

};


//###############################################################################
/* Returns the length of the 1D system starting from "ind" in direction "dir" 
   The dimension, that stops the system will be saved in "boundingdim".
*/
template<int sg, bool with_boundingdim>
int Grid::get_system_length( const position_t& dir, const  boundary_iterator<sg>& it, int& boundingdim)  const 
{
    int ret_val = maxsize + 1;
    for(uint j=0; j < D; ++j) {
		if (dir[j] > 0){
			int temp = sizes[j]-1-it.pos()(j);
			if (ret_val > temp ) {
				if(with_boundingdim) boundingdim = 2*j+1;
				ret_val = temp;
			}
		} else if (dir[j] < 0) {
			if (ret_val > it.pos()(j)) {
				if(with_boundingdim) boundingdim = 2*j;
				ret_val = it.pos()(j);
			}
		}
    }
    return ret_val - sg+1;
} 

}

using _grid::Grid;
typedef std::shared_ptr<const Grid> grid_ptr;
typedef const Grid& grid_ref;

}
