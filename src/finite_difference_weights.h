#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <utility>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <algorithm>
#endif

#include "finite_difference_weights_config.h"
#include "auto_connecting_body.h"
#include "grid_fields.h"
#include "utility.h"
#include "memoized_node_factory.h"
#include "grid_striding_info.h"
#include "dir_container.h"
#include "Eigen.h"

namespace fipster { namespace finite_difference_weights {

	using namespace std;
	using namespace Eigen;

	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################

	
	/**\brief The internal storage for (raw) A_io.
			
		It represents a banded matrix that is stored in a dense
		Eigen::Matrix. It is only accessed via the dir_container interface, from
		which A_io_t derives.
		*/
	typedef Matrix<double,Dynamic,Dynamic,ColMajor> data_t;


	/** \brief Storage wrapper type that endows the raw grid_field
		with dir_container interface	

			The columns of the internal storage #data will be
			references as dir_collection making it possible to safely access
			the left, right and center weight for each grid point and direction via:
			\code
at(some_dir_index)[some_grid_index]; //left
at_opposite(some_dir_index)[some_grid_index]; //right
at_center(some_dir_index)[some_grid_index]; //middle
			\endcode
		*/
	template<typename index_t,typename data_t2=data_t>
	struct grid_field_dir_container : 
		dir_container<grid_field<index_t,typename data_t2::ColXpr>,_dir::center>{

			typedef grid_field<index_t,data_t2> grid_field_data_t;
	  	using vec_t = typename grid_field_dir_container::vec_t;
	  	using dir_cont_t = typename grid_field_dir_container::dir_cont_t;

	private:
			grid_field_data_t data;

	public:
			index_t resize(grid_ref grid,uint nDirs){
				uint nPoints = this->detemine_size(nDirs);
				//resize the 2D Array
				auto ret=data.resize(grid,nPoints);
				//create columns from that array
				vec_t temp; temp.reserve(nPoints);
				for(uint i=0;i<nPoints;i++)
					temp.push_back(data.col(i));
				//use columns in dir_container
				dir_cont_t::reset(move(temp));
				return ret;
			}

			/*grid_field_data_t& UNSAFE_data_access(){
				return data;
			}*/

			void setZero(){
				data->setZero();
			}
			template<class A>
			void setConstant(const A& t){
				data->setConstant(t);
			}
			/* TEST in "finite_difference_weights.cpp":
			assert(&A_io[i][it.inner_start_g_ind()]==&A_io.data.at2(it.inner_start_g_ind(),i.left()));
			assert(&A_io.at_opposite(i)[it.inner_end_g_ind]==&A_io.data.at2(it.inner_end_g_ind,i.right()));
			*/
		};

	typedef pair<const config_t, grid_ptr> fd_weights_arg_t;

	/** \brief contains the part of the operator, that involves the boundary.
			
			As the stencil is compact (only nearest neighbor), there is only
			a connection between the first inner boundary "j_sgs" and the first outer
			boundary "j_sgs-1".
		*/
	typedef sparse_field<boundary_index<J_SGS>,boundary_index<J_SGS-1>,RowMajor> A_o_t;

	/** \brief Result type of the finite difference weights body
	*/
	struct fd_weights_result_t{
		static const int sg = J_SGS;///< Specifies on what subgrid weights should be calculated
		
		typedef shared_ptr<grid_striding_info_collection<sg>> gsic_ptr; 
		gsic_ptr grid_striding_infos;///<Contains the stridigin information for the used stencil type
		
		
		typedef grid_field_dir_container<grid_index<sg>> A_io_t; ; ///<operator A_io in a custom storage format
			
		A_io_t A_io;

		
		shared_ptr<A_o_t> A_o; ///<sparse matrix A_o ("outer")

		//ctor
		fd_weights_result_t():A_o(new A_o_t){}
	};


	// body_t ######################################################################
	template<class PDE_core_tpl>
	struct body_t : auto_connecting_body<body_t<PDE_core_tpl>,fd_weights_result_t> {

		
		typedef PDE_core_tpl PDE_core;
		typedef shared_ptr<PDE_core> PDE_core_ptr;
		PDE_core_ptr pde;
		
		using result_sptr_t = typename body_t::result_sptr_t;
		using result_t = typename body_t::result_t;
		using input_t = typename body_t::input_t;

		static const int sg = fd_weights_result_t::sg;

		const fd_weights_arg_t args;
		uint D;
		
		body_t(const fd_weights_arg_t& args,shared_ptr<const PDE_core> pde_ptr);

		result_sptr_t operator()(const input_t&);
	};



	// #################################################
	// ##############   PHASE 2:			############
	// ##############   Node Factory        ############
	// #################################################
	struct fdweights
		: memoized_node_factory_singleton<fdweights,fd_weights_arg_t,fd_weights_result_t>
	{
		sender_ptr setup_node(const arg_t& arg);
	};

}

using finite_difference_weights::A_o_t;
using finite_difference_weights::fd_weights_result_t;
typedef shared_ptr<fd_weights_result_t> fd_weights_result_ptr_t;
using finite_difference_weights::fd_weights_arg_t;
using finite_difference_weights::fdweights;
using finite_difference_weights::grid_field_dir_container;

}
