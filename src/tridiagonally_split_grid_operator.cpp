#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#endif

#include "tridiagonally_split_grid_operator.h"
#include "spacetime_field.h"
#include "exceptions.h"

namespace fipster {

void tridiag_systems_collection::test_gather_scatter(const Grid& g){
	
	vector<double> lc(striding_info->max_length);

	discretized_space_field<J_SGS> gl("test_gather_scatter.gl"),gl2("test_gather_scatter.gl2");
	gl.resize(g);
	gl2.resize(g);

	//fill in test values:
	double j=0;
	for(auto i=g.begin<J_SGS>();i!=g.end<J_SGS>();++i,j++){
		gl.at(i)=j;
		gl2.at(i)=-1;
	}

	//copy values via gather scatter
	for(auto sys=begin();sys.valid();++sys){
		sys.gather(gl,lc);
		sys.scatter(gl2,lc);
	}

	//check for test-values:
	j=0;
	for(auto i=g.begin<J_SGS>();i!=g.end<J_SGS>();++i,j++){
		FIPSTER_ASSERT(gl2.at(i)==j);
	}
	FIPSTER_ASSERT((int)j==striding_info->cumulative_length);
}


/**
tests, if the tridiagonally_split_grid_operator really references N distinct
memory locations, where N is the number of matrix elements in all tridiagonal matrices
(including boundary elements)
*/
void tridiagonally_split_grid_operator::test_locations() const{
	
	uint N = gsic->grid->N[J_SGS]*nDirs()*3;
	vector<const double*> addresses;
	addresses.reserve(N);

	//loop over directions:
	for(auto dir = dir_it();dir.valid();++dir)
	{
		//loop over systems:
		for(auto sys=at(dir).begin();sys.valid();++sys){

			//loop over system:
			int l=sys.connected_set_it->length;
			for(int i=0;i<l;i++){
				addresses.push_back(&sys.M.diag(i));
				addresses.push_back(&sys.M.superdiag(i));
				addresses.push_back(&sys.M.subdiag(i));
			}					
		}
	}

	FIPSTER_ASSERT(rearrange_and_get_number_of_uniques(addresses)==N);
}

}
