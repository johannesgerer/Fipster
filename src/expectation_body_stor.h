#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/throw_exception.hpp>
//gcc: #include <mkl_lapack.h>
#endif

#include "Eigen.h"
#include "integer_types.h"
#include "cryer.h"
#include "logging.h"
#include "exceptions.h"
#include "tridiagonally_split_grid_operator.h"
#include "grid_fields.h"

namespace fipster { namespace expectation_values {

	
	/** indicates which of the beta is not 0.0 (trivial): 0,1,2(both)
	*/
	enum non_zero_t{
		only_first,only_second,both
	};

	using namespace Eigen;

	struct storage{

		//per time_stepp stuff:

		boundary_field_t<J_SGS> Bo1_BC_inhom_temp,Bo1_BC_inhoms;
		cryer<tridiag_systems_collection::tridag_matrix_t,VectorXd,VectorXd> Cryer;



		//per system stuff:

		VectorXd			g0,f1,q,f0;	

		storage()
			:Cryer(&q,&f0,1e-10),length(0)
		{
			thread_logger()<<"EXPECT | Accessing Time Stepping TLS"<<endl;
		};

		uint length;

		void resize(uint n){
			length=n;
			g0.resize(n);
			f1.resize(n);
			q.resize(n);
			f0.resize(n);
		}
		template<class T>
		void cryer_lcp(T& sys,bool check=false){
			Cryer.compute<1,true,false>(&sys.M,length);
		}

		template<class T>
		void tridiagonal(T& sys,bool check=false){
			Cryer.compute<false,true,true>(&sys.M,length);
		}
	};

}

}
