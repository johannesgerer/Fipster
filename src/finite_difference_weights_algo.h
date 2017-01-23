#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <vector>
#include <cstdlib>
#include <cmath>
#endif

#include "Eigen.h"
#include "logging.h"
#include "TLS.h"
#include "least_squares.h"
#include "finite_difference_stencils.h"

namespace fipster { namespace finite_difference_weights {

	using namespace Eigen;
	using namespace tbb;

	
	
	struct references{
		stencil_points_t			stencil_points;//nPoints
		Matrix<int,Dynamic,Dynamic,RowMajor> multi_indices;//numPDs,D
		VectorXd					factorialdenoms;//numPDs
		const uint D;

		references(	uint D)
			: D(D){}

		int generatePDindices(uint order);
		void generateFactorialDenoms(double HigherOrderWeight);
	};

	combinable<least_squares_stor> tls;

	//############################################################################
	/** \brief Body for parallel weights calculation
	*/
	template<typename result_data_t,class body_t>
	struct optimize_body{

		typedef Matrix<double,1,Dynamic> single_t;

		references				&refs;
		result_data_t			&result_data;
		grid_ref				grid;
		double					max_norm;
		const uint				D;
		const int				numPDs;

		typedef const typename body_t::PDE_core_ptr pde_ptr;
		pde_ptr					PDE;

		optimize_body(references& refs,result_data_t &result_data
									,pde_ptr PDE,grid_ref grid)
			:refs(refs)
			,result_data(result_data)
			,grid(grid)
			,max_norm(0)
			,D(refs.D)
			,numPDs(refs.factorialdenoms.size())
			,PDE(PDE)
		{}

		optimize_body(optimize_body& y,split)
			:refs(y.refs)
			,result_data(y.result_data)
			,grid(y.grid)
			,max_norm(0)
			,D(y.D)
			,numPDs(y.numPDs)
			,PDE(y.PDE)
		{}

		void join(optimize_body& y){max_norm=max(max_norm,y.max_norm);}

		//(parallel) WORK
		void operator()(const blocked_range<grid_index<body_t::sg>>& range){
			
			//thread_logger()<<"von "<<range.begin()<<" bis "<<range.end()<<endl;
			//tbb::this_tbb_thread::sleep(tbb::tick_count::interval_t(0.01));

			least_squares_stor &temps = tls.local();
			temps.resize(numPDs,refs.stencil_points.proper_neighs(),1);
			

			//Main loop
			auto state_it = grid_iterator<body_t::sg,true>(grid,range.begin());
			for(;state_it!=range.end();++state_it){

				
				//cout<<state_it.ind.ind<<": "<<toS(state_it.state)<<endl;

				/*########## construct the rhs-vector b */
				for(int k=0;k<numPDs;k++)
					temps.b[k] = (*PDE)(refs.multi_indices.row(k),state_it.state(),D);

					//cout<<temps.b.transpose()<<endl;

				/*######### construct the coefficient matrix D */
				{auto dir_it = refs.stencil_points.dir_it();
				for(int i=0;dir_it.valid();++dir_it,i+=2) {//iterate through all directions

					auto& point = refs.stencil_points.at(dir_it);
					auto coeffrow0 = temps.d.col(i);
					auto coeffrow1 = temps.d.col(i+1);
					//copy denominators
					coeffrow0 = refs.factorialdenoms;
					coeffrow1 = refs.factorialdenoms;


					for (int k = 0; k < numPDs; ++k) {
						auto& coeff0 = coeffrow0(k);
						auto& coeff1 = coeffrow1(k);
						auto multi_index = refs.multi_indices.row(k);
					
						for (uint l = 0; l < D ; ++l) {
							/* if the derivative in this direction is present */
							int& ind = multi_index(l);
							if (ind > 0) {
#define stepsizes(i,j) grid.stepsizes[i][j]	
								switch(point[l]){
								case 1:		/* positive component */
									coeff0 *= pow(stepsizes(l,state_it.pos()[l]),ind);
									coeff1 *= pow(-stepsizes(l,state_it.pos()[l]-1),ind);
									break;
								case -1:	/* negative component */
									coeff0 *= pow(-stepsizes(l,state_it.pos()[l]-1),ind);
									coeff1 *= pow(stepsizes(l,state_it.pos()[l]),ind);
									break;
								case 0:	/* no component */
									coeff0 = 0.0;
									coeff1 = 0.0;
									goto nextCoeff;
									break;
								}
#undef stepsizes
							}
						}	
nextCoeff:;
					}
				}}



				/*######### find optimal weights */

				#if 1 //using non-negative general least squares
					max_norm = max(max_norm,temps.calculate());
				#else //using general least squares
					max_norm = max(max_norm,temps.gels());	
				#endif

				auto dir_it=refs.stencil_points.dir_it();


				/* 3.2.4 calculate center coeff. REMARK: THIS IS NOT CENTER WEIGHT */
				result_data.at_center(dir_it).at(state_it.ind) = (*PDE)(state_it.state()); 
				/* TO GET CENTER WEIGHT ONE NEEDS: - SUM(StencilWeights(1:stencil.nPoints,j)) 
				this is done in splitter in a way that yields 1dimensional negative M-Matrices */
				
				//cout<<temps.result.transpose()<<endl;

				//######### copy the results (deliberately differing storage order)
				for(int i=0;dir_it.valid();i+=2,++dir_it){
					result_data.at(dir_it).at(state_it.ind) = temps.result(i); ///\todo TEST : do the weights really approximate the PDE as given by the norm and are non-negative? (außerhalb dieser function natürlich)
					result_data.at_opposite(dir_it).at(state_it.ind) = temps.result(i+1);
				}
			}

		//	thread_logger()<<"DONE, maxNorm in Thread:"<<scientific << max_norm <<endl;
		}
	};
}}
