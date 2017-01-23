#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/property_tree/ptree.hpp>
#include <boost/throw_exception.hpp>
#include <exception>
#include <memory>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#endif

#include "grid.h"
#include "utility.h"
#include "grid_tests.h"
#include "utility.h"
#include "global_config.h"
#include "grid_iterators.h"
#include "configuration_parser.h"


void cut_off(double& a){ // or you need something like CodeRefTHR
	a = abs(a)<1e-15 ? 0 : a;
}


namespace fipster { namespace _grid {

	state_variable_t::state_variable_t( const ptree& pt ) 
		:	n(betterGet<int>(pt,"n","state variable"))		
		,	spacing_node(betterChild(pt,"spacing","state variable "+toS(n)))
		, hits_zero(spacing_node.get("<xmlattr>.hitsZero",false))
	{
		lower_bound = pt.get<double>("lowerBound");
		upper_bound = pt.get<double>("upperBound");
		FIPSTER_ASSERT(lower_bound<upper_bound);
		tolerance = pt.get("tolerance",1e-17);

		spacing_type = spacing_node.get<string>("<xmlattr>.type");
	}

	template<class A>
	void state_variable_t::fill_coords( A& coords, A& stepsizes )
	{
#define coords(i) coords[i]
#define stepsizes(i) stepsizes[i]
		// homogeneous #########################################
		if(spacing_type == "homogeneous"){
			double s;
			s = (upper_bound-lower_bound)/(n-1);
			for(uint i = 0; i<n ; ++i )
				coords(i) = i*s + lower_bound;					

			// exponential #########################################
		}else if(spacing_type == "exponential"){
			double a,b,h = spacing_node.get<double>("h");
			a = (upper_bound-lower_bound)/(exp(h)-1.0);
			b = h/(n-1);
			for(uint i=0; i < n; ++i)
				coords(i) = a*(exp(i*b)-1.0)+lower_bound;

			// not implemented #########################################
		}else
			BOOST_THROW_EXCEPTION(invalid_argument("Grid spacing '"+spacing_type+"' not implemented"));

    //shift to include 0
		double shift = 0;
		if(hits_zero){
			vector<double> temp=coords;
			for(auto& t: temp) t=abs(t);
			shift=coords(min_element(begin(temp),end(temp))-begin(temp));
			for(auto& t: coords) t-=shift;
		}
	
		//produce stepsizes
		for(uint i=0;i<n-1;i++){
			cut_off(coords(i));
			stepsizes(i) = coords(i+1)-coords(i);
			
		}

		// validate coordinates and stepsizes  #########################################
		double s=abs(coords(n-1)-upper_bound+shift);
		if(s >= tolerance*(upper_bound-lower_bound) 
			|| s!=s ) //This checks for NaN and Inf and stuff
			BOOST_THROW_EXCEPTION(underflow_error
														("coords did to not add up (rel.Err: "+toString(s)+
														 ", tolerance: "+toString(tolerance)+")"));

#undef coords
#undef stepsizes
	}

	


	Grid::Grid( const ptree& pt )
		: state_variables_config_t(pt)
		, id(betterGet<string>(pt,"<xmlattr>.id","grid"))
		, max_subgrid(pt.get_optional<double>("maxSubgrid"))
	{
		//infer the dimension of the grid
		D = sv_vec.size();
		sizes.resize(D);

		maxsize=0;
		for(uint i=0;i<D;i++){
			sizes[i] = sv_vec[i].n;
			maxsize = max(maxsize,sizes[i]);
		}

		//Strides and Number of points
		N.resize(implemented_inner_grids);
		strides.resize(implemented_inner_grids,D);
		Nbound.resize(implemented_inner_grids-1);

		for(uint j=0; j < implemented_inner_grids; ++j) 
		{
			//Strides
			strides(j,0) = 1;
			for(uint i=1; i < D; ++i)
				strides(j,i)=strides(j,i-1)*(sizes[i-1]-j*2);

			//Number of points
			N[j]=strides(j,D-1)*(sizes[D-1]-j*2);

			//check if sizes are sufficient to implement up to implemented_inner_grids!
			if(max_subgrid && *max_subgrid>=j && sizes[D-1]<=(int)j*2)
				BOOST_THROW_EXCEPTION(
					runtime_error("grid size in dimension "+toS(j)+" is too small for inner_grid "+toS(j)));
		}

		//Number of boundary points
		for(uint j=0; j < implemented_inner_grids-1; ++j) 
			Nbound[j] = N[j] - N[j+1];

		//origin offsets
		for(uint j=0; j < implemented_inner_grids; ++j) 
			origin_offset[j] = strides.row(j).sum()*j;


		//reserve space
		coords.resize(D);
		stepsizes.resize(D);

		//fill coordinates and stepsizes
		for(uint i=0;i<D;i++){
			coords[i].resize(sizes[i]);
			stepsizes[i].resize(sizes[i]-1);
			sv_vec[i].fill_coords(coords[i],stepsizes[i]);
		}

		//test iterator logic:
		if(global_config::get().test_grid_iterators)
			test_iterator_logic(*this);
	}


	ostream& operator<<(ostream& out,const Grid& grid)
	{
		auto state_it = grid_iterator<0,true>(grid);
		for(;state_it!=grid.end<0>();++state_it){
			out<<state_it.ind<<"\t";
			for(uint i=0;i<grid.D;i++)
				out<<boost::lexical_cast<string>(state_it.state()[i])<<"\t";
			out<<"\n";
		}
		out.flush();		
		return out;		
	}

		
}}
