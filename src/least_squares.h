#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/throw_exception.hpp>
#endif

#include "Eigen.h"
#include "nnls.h"
#include "logging.h"
#include "exceptions.h"

namespace fipster { namespace _least_squares {


	using namespace Eigen;

	struct least_squares_stor{
		Matrix<double,Dynamic,Dynamic,ColMajor>	d;
		VectorXd						b;
		VectorXd						result;
		VectorXd						w;
		VectorXi						index;
		int in,im;
		//ctor
		least_squares_stor():in(0),im(0){};
		void resize(int m,int n,int additional_result_size){
			if(im!=m || in!=n)
			{
				thread_logger()<<"WEIGHTS| Accessing Least Squares TLS"<<endl;
				im=m;in=n;
				d.resize(m,n);
				b.resize(m);
				w.resize(n+m);
				index.resize(n);
				result.resize(n+additional_result_size);
			}
		}
		//returns the residuum norm
		double calculate(){
			double norm;
			switch(nnls_c(d.data(),&im,&in,b.data(),result.data(),&norm,w.data(),index.data())){
			case 2:
				BOOST_THROW_EXCEPTION(invalid_argument("THE DIMENSIONS OF THE PROBLEM ARE BAD. EITHER M .LE. 0 OR N .LE. 0."));
				break;
			case 3:
				BOOST_THROW_EXCEPTION(invalid_argument("ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS."));
				break;
			}

			return norm;
		}

		// double gels(){ //optimize without constraint (using general least squares)
		// 	//#include <mkl_lapack.h>
		// 	//LLS routienen allgemein http://software.intel.com/sites/products/documentation/hpc/composerxe/en-us/mklxe/mkl_manual_win_mac/lse/lse_drllsp.htm#tbl4-8
		// 	//GELS http://software.intel.com/sites/products/documentation/hpc/composerxe/en-us/mklxe/mkl_manual_win_mac/lse/functn_gels.htm
		// 	// LAPACK and C++ http://software.intel.com/sites/products/documentation/hpc/mkl/mkl_userguide_lnx/index.htm#GUID-ABCC618B-43C4-4DCD-ADA2-6F061B5116CD.htm
		// 	double residuumnorm=numeric_limits<double>().quiet_NaN();
		// 	int info=-1;
		// 	int one=1;
		// 	int lwork=im+in;

		// 	if (im >= in) {
		// 		DGELS("N",&im,&in,&one,d.data(),&im,b.data(),&im,w.data(),&lwork,&info);
		// 		//copy b (which contains the LLS solution) into b
		// 		result.topRows(in)=b.topRows(in);
		// 		residuumnorm=b.bottomRows(im-in).squaredNorm();
		// 	} else {//im < in 
		// 		//Das problem bei diesem "underdetermined system" ist, dass 
		// 		//die "minimum norm" solution aus allen möglichen ausgewählt wird.
		// 		//diese hat i.A. garnichts, mit der von NNLS zu tun.
		// 		//Alternativ könnte man durch expliziete auswahl von Algorithmen
		// 		//oder einer Abwandlung des LLS Problems näher kannkommen:
		// 		//Mehr dazu siehe hier:
		// 		//http://eigen.tuxfamily.org/dox-devel/TutorialLinearAlgebra.html#TutorialLinAlgLeastsquares
		// 		//LLS allgemein http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/cpp/win/mkl/refman/lse/lse_intro.html
			
		// 		//copy b into result 
		// 		result.topRows(im)=b.topRows(im);
				
		// 		/* gels(D,StencilWeights(1:,j),'N',info) */
		// 		DGELS("N",&im,&in,&one,d.data(),&im,result.data(),&in,w.data(),&lwork,&info);
		// 		residuumnorm = 0.;
		// 	}
		// 	if(info != 0)
		// 		FIPSTER_THROW_EXCEPTION(runtime_error("error in GELS"+toS(info)+"  Res.norm:  "+toS(residuumnorm)));
		// 	return residuumnorm;
		// }
	};

}

using _least_squares::least_squares_stor;

}
