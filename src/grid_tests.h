#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/throw_exception.hpp>
#include <exception>
#endif

#include "grid.h"
#include "grid_iterators.h"

namespace fipster { namespace _grid {

	


template<int sg,class Grid_t>
struct iterator_tester{


	static void test(const Grid_t& g){
		if(!g.max_subgrid || *g.max_subgrid >= sg){ 
//PUT TESTS HERE ##################################################################
		if(sg<Grid_t::implemented_inner_grids-1){ //Boundary tests

			auto b_end=g.template b_end<sg>();
			uint i=0;
			int last_ind=-100000,last_b_ind=-10000;
			typedef boundary_iterator<sg,true,true> bi_t;
			boundary_index<sg> bind;
			auto it2=boundary_iterator<sg>(g);
			for(auto it = bi_t(g);it!=b_end && it2!=b_end;++it2,++it,i++){


				for(uint j=0;j<g.D;j++)
					switch(betterAt(it.active_set,j,"")){
					case lower_ind:
						FIPSTER_ASSERT(it.pos()[j]==sg);
						break;
					case upper_ind:
						FIPSTER_ASSERT(it.pos()[j]==g.sizes[j]-sg-1);
						break;
					case inactive_ind:
						FIPSTER_ASSERT(it.pos()[j]<g.sizes[j]-sg-1 && it.pos()[j]>sg);
						break;
					}

				
				grid_index<sg> asd(it.pos(),g);
				bind.template make_ind_from_pos<true>(it.pos(),g,asd);

				
				auto stit=grid_iterator<sg,true,Grid_t>(g,asd);
				FIPSTER_ASSERT(it.state() == stit.state());

				FIPSTER_ASSERT(it.g_ind==asd);
				FIPSTER_ASSERT(it.b_ind==bind);
				FIPSTER_ASSERT(it.b_ind==it2.b_ind);
				FIPSTER_ASSERT(it.pos()==it2.pos());
				FIPSTER_ASSERT(it.g_ind==it2.g_ind);
				FIPSTER_ASSERT(i==0|| (last_b_ind<it2.b_ind.ind && last_ind<it2.g_ind.ind));
				last_b_ind=it2.b_ind.ind;
				last_ind=it2.g_ind.ind;
			}
			//this together with the monotonicity check assures,
			//that all boundary points 
			FIPSTER_ASSERT(i==g.Nbound[sg]);
		}
		{//interior tests
			grid_iterator<sg,true,Grid_t> t(g);
			auto e=g.template end<sg>();
			uint i=0;
			for(auto it=g.template begin<sg>();it!=e;i++,++t,++it){
				auto stit=grid_iterator<sg,true,Grid_t>(g,it);
				FIPSTER_ASSERT((t.pos() == stit.pos()));
				FIPSTER_ASSERT(t.state() == stit.state());
			}
			FIPSTER_ASSERT(i==g.N[sg]);
		}
//PUT TESTS HERE - END ###########################################################
		}
		iterator_tester<sg-1,Grid_t>::test(g);
	}
};

template<class Grid_t>
struct iterator_tester<-1,Grid_t>{
	static void test(const Grid_t& g){}
};

template<class Grid_t>
void test_iterator_logic(const Grid_t& g){
	iterator_tester<Grid_t::implemented_inner_grids-1,Grid_t>::test(g);
}

}
using _grid::test_iterator_logic;
}

/*
template<int n> 
void Grid::test_grid_interateboundary_canonical_order() const
{
	int k=0,i=0,te,ab;
	int1 ind(D),work(D);

	list<string> l;
	k=1;
	while(this->iterateonboundary_canonical_order<n,false,0>(i,ind)){
		te=get_boundary_array_offset<n,true>(ind,ab);
		l.push_back(toS(ind)+" "+toS(k)+" "+toS(te));

		int j;
		for( j=1 ; j < D+1; ++j){ 
			if( (ind(D-j)==n) != (ab==2*(D-j)) ||
				(ind(D-j)==sizes[D-j]-n) != (ab==2*(D-j)+1)
				)
				trace("wrong \"active boundary\"!",new std::exception(""));
			if(ind(D-j)==n || ind(D-j)==sizes[D-j]-n) break;
		}
		if(j==D+1) trace("not on boundary",new std::exception(""));
		else trace("boundary "+toS(D-j)+" "+toS(ab));
		if(k!=te) trace("get_boundary_array_offset error",new std::exception(""));
		k++;
	}

	if(l.size() != nboundarylayerpoints(n))
		trace("problmlebmlebmlembl",new std::exception(""));

}*/
