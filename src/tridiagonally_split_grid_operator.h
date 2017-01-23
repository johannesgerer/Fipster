#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
//gcc #include <mkl_lapack.h>
#include <vector>
#endif

#include "integer_types.h"
#include "dir_container.h"
#include "grid_striding_info.h"
#include "exceptions.h"
#include "global_config.h"

namespace fipster { namespace _tridiagonally_split_grid_operator {

	
//typedef vector<double>::iterator matrix_data_it_t;

enum storage_type { rowwise, diagonalwise,systemwise_diagonals};

/** this is used to access the data of a tridiagonal matrix
	It does not possess own data. As this class handles the data access,
	it is the place to define different storage orders/schemes.
	*/
template<storage_type storage, class matrix_data_it_t>
class tridiag_matrix{

public:
	/** The routine performs a matrix-vector product of the form:
	y += alpha*this*x
	*/
	template<class T,int alpha>
	void times_vec(T& y,const T& x)const{
		

		if(length==1)
			y[0]+=alpha*x[0]*diag(0);
		else{
			//\todo performance: this can be sped up using Eigen auto vectorization, but requires aligned data to work!
			uint i=0;y[i]+=alpha*(x[i]*diag(i)+x[i+1]*superdiag(i));
			for(i++;i<length-1;i++)
				y[i]+=alpha*(x[i]*diag(i)+x[i+1]*superdiag(i)+x[i-1]*subdiag(i));
			y[i]+=alpha*(x[i]*diag(i)+y[i-1]*subdiag(i));
			
		}
	}


#if 0
	
	//Mathematical operations #########################################################
	
	//multiply with vector:
	//The routine performs a matrix-vector product of the form:
	// b := alpha*THIS*x + beta*b
	//(alpha ,beta may be 0.0, 1.0, or -1.0)
	void times_vec(const double1& x, double1& b,const double alpha=1.0,const double beta=0.0){
		if(bound_checking)
			if(d.length>x.length || d.length > b.length)
				trace("Error in times_vec, wrong dimensions"+toS(d.length)+" "+
				toS(x.length)+" "+toS(b.length),new invalid_argument(""));

		const char trans='T';
		const int one = 1;

		//http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/cpp/win/mkl/refman/lau/functn_lagtm.html
		DLAGTM(&trans,&d.length,&one,&alpha,dl.getData(),d.getData(),du.getData(),x.getData(),
			&d.length,&beta,b.getData(),&d.length);
		///\todo  offset and dl(0) and du(n-1)
	}

#endif

	matrix_data_it_t data_begin;
	array<matrix_data_it_t,3> its;
	uint cumulative_length;

	int stride() const{

		static_assert(storage==rowwise ||storage==systemwise_diagonals||storage==diagonalwise,
			"bad storage order");

		switch(storage){
		case rowwise:
			return 3;
		case diagonalwise:
		case systemwise_diagonals:
			return 1;
		default:
			return 4893484;
		}
	}

public:
	
	uint length;

	tridiag_matrix(){};

	tridiag_matrix(matrix_data_it_t data_begin, int cumulative_length)
		:data_begin(data_begin),cumulative_length(cumulative_length)
	{}

	string spos(int i){
		return "at "+toS(i)+"/"+toS(length);
	}

	/** checks for M-Matrix property 
	\param negative indicates, that the negative version of the matrix is checked
	*/
	template<bool negative, bool start, bool end>
	bool check_M_property(int i, bool do_throw){
		int factor=1;
		if(negative) factor=-1;
		
		if(factor*diag(i)<=0)
			if(do_throw)
				FIPSTER_THROW_EXCEPTION(runtime_error("Diagonal not positive "+spos(i)));
			else return false;
		else if((!start && factor*subdiag(i)>0) || (!end && factor*superdiag(i)>0) )
			if(do_throw)
				FIPSTER_THROW_EXCEPTION(runtime_error("Off-Diagonal positive "+spos(i)));
			else return false;
		else{
			double temp=diag(i)+(!end?superdiag(i):0)+(!start?subdiag(i):0);
			if(factor*temp<0)
				if(do_throw)
					FIPSTER_THROW_EXCEPTION(runtime_error("not diagonally dominant "+spos(i)));
				else return false;
			else if(temp==0){
				if(do_throw)
					FIPSTER_THROW_EXCEPTION(runtime_error("not strictly diagonally dominant "+spos(i)));
				else return false;
			}
		}
		return true;
	}


	const double& diag(int i) const{ return its[0][stride()*i]; }
	const double& subdiag(int i) const{ return its[1][stride()*i]; };
	const double& superdiag(int i) const{ return its[2][stride()*i]; };

	double& diag_ref(int i) { return its[0][stride()*i]; }
	double& subdiag_ref(int i) { return its[1][stride()*i]; };
	double& superdiag_ref(int i) { return its[2][stride()*i]; };

	void init(int cumulative_offset, int length)
	{
		this->length=length;

		switch(storage){
		case rowwise:
			its[0]=data_begin+3*cumulative_offset;
			its[1]=its[0]+1;
			its[2]=its[1]+1;
			break;
		case diagonalwise:
			its[0]=data_begin+cumulative_offset;
			its[1]=its[0]+cumulative_length;
			its[2]=its[1]+cumulative_length;
			break;
		case systemwise_diagonals:
			its[0]=data_begin+3*cumulative_offset;
			its[1]=its[0]+length;
			its[2]=its[1]+length;
		}
	}
};

template<storage_type storage, class matrix_data_it_t>
ostream& operator<<(ostream& out,const tridiag_matrix<storage,matrix_data_it_t>& M){
	cout<<"super: ";
	for(uint i=0;i<M.length;i++)
		cout<<M.superdiag(i)<<" ";
	cout<<endl<<"diag: ";
	for(uint i=0;i<M.length;i++)
		cout<<M.diag(i)<<" ";
	cout<<endl<<"sub: ";
	for(uint i=0;i<M.length;i++)
		cout<<M.subdiag(i)<<" ";
	return cout;
}

/**  \brief Iterator, that takes care of everything needed to be performed on/with a tridiagonal system:

	1. handles access to the diagonals of a tridiagonal matrix via the use of #tridiag_matrix.
	2. Performs gather/scatter between local and global storage schemes for the vectors, the tridiag matrices should act on.
	
	It does not possess own data. It instead manages the bulk data of the tridiag_systems_collection.

	M_ii = d_i, M_(i,i-1) = dl_i, M_(i,i+1) = du_i 
	NOTE: du(n-1) and dl(0) are not part of the square matrix , but are required! See CodeRefA
*/
template<class matrix_data_it_t>
class tridiag_systems_iterator{
	
public:
	
	//rowwise, systemwise_diagonals or diagonalwise,
	typedef tridiag_matrix<rowwise,matrix_data_it_t> tridag_matrix_t;
	tridag_matrix_t M;

	//grid striding info stuff:
	typedef grid_striding_info_t<J_SGS>::vec_t::const_iterator connected_set_it_t;
	connected_set_it_t connected_set_it,connected_set_end;
	grid_stride<J_SGS> stride;

	tridiag_systems_iterator(){};

	/** construct the iterator, by setting up the iterator to the first connected_set 
		and passing the data iterator to the tridiag_matrix member #M
		*/
	template<class A>
	tridiag_systems_iterator(A& container)
		:M(container.data.begin(),container.striding_info->cumulative_length)
		,connected_set_it(container.striding_info->connected_sets.begin())
		,connected_set_end(container.striding_info->connected_sets.end())
		,stride(container.striding_info->stride)
	{
		reinit_M();
	}


	/** comparison, to be used in loops
	*/
	bool operator!=(const connected_set_it_t& that)const
	{
		return connected_set_it!=that;
	}

	/** prefix increment, to be used in loops
	*/
	void operator++(){
		connected_set_it++;
		if(valid())
			reinit_M();
	}

	bool valid(){
		return connected_set_it!=connected_set_end;
	}

	void reinit_M(){
		M.init(connected_set_it->cumulative_offset,connected_set_it->length);
	}

	/*gather the distributed elements of a strided system within a global vector (spacetime_field)
	and put them in a local dense smaller vector
	*/
	template<class global_t, class local_t>
	void gather(const global_t& gl, local_t& lc){
		auto it=connected_set_it->begin_gi;
		const uint l=connected_set_it->length;
		for(uint i=0;i<l;i++,it+=stride){
			lc[i] = gl.at(it);
		}
	}
	/*scatter a dense local vector representing a strided system
	into the larger global vector (spacetime_field)
	*/
	template<class global_t, class local_t>
	void scatter(global_t& gl, const local_t& lc){
		auto it=connected_set_it->begin_gi;
		const uint l=connected_set_it->length;
		for(uint i=0;i<l;i++,it+=stride){
			gl.at(it)=lc[i];
		}
	}
};


/** Contains an interface to perform tridiagonal matrix operations
	on strided systems on the grid using the information stored
	in a grid_striding_into_t object and allocated memory for the
	coefficients.
	*/
struct tridiag_systems_collection
{
	
	typedef const grid_striding_info_t<J_SGS>* striding_info_ptr;	

	striding_info_ptr striding_info;
		
	typedef vector<double> data_t;
	typedef tridiag_systems_iterator<data_t::const_iterator> const_sys_iterator ;
	typedef tridiag_systems_iterator<data_t::iterator> sys_iterator ;
	data_t data;
	

	typedef tridiag_systems_iterator<data_t::const_iterator>::tridag_matrix_t 
		tridag_matrix_t;

	void test_gather_scatter(const Grid& g);

	sys_iterator begin(){
		return sys_iterator(*this);
	}

	const_sys_iterator begin() const {
		return const_sys_iterator(*this);
	}

	void init(striding_info_ptr striding_info_){
		striding_info=striding_info_;
		data.resize(striding_info->cumulative_length*3);
	}

};

/** \brief A operator on the grid, which is split in tridiagonal operators,
			one for each direction
*/
struct tridiagonally_split_grid_operator 
	: dir_container<tridiag_systems_collection,_dir::simple>{

		typedef dir_container<tridiag_systems_collection,_dir::simple> parent_t;

		typedef shared_ptr<grid_striding_info_collection<J_SGS>> gsic_ptr; 
		gsic_ptr gsic; ///< Grid striding info collection 


		//ctor
		tridiagonally_split_grid_operator(){};

		/** initialize the tridiagonal grid operator using
			grid_striding_information.
			A shared_ptr to it will be kept
		*/
		void init(gsic_ptr gsic_)
		{
			gsic = gsic_;
			parent_t::resize(gsic->nDirs());
			for(auto i = gsic->dir_it();i.valid();++i){
				at(i).init(&gsic->at(i));
				if(global_config::get().test_gsi_gather_scatter)
					at(i).test_gather_scatter(*gsic->grid);
			}

			if(global_config::get().test_tridiagonally_split_grid_operator)
				test_locations();
		}

		void test_locations() const;
};

}
using _tridiagonally_split_grid_operator::tridiagonally_split_grid_operator;
using _tridiagonally_split_grid_operator::tridiag_systems_collection;

}
