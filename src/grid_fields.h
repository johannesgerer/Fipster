#pragma once
#include "grid.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/lexical_cast.hpp>
#include <string>
#endif

namespace fipster { namespace _grid {

/** \brief This wraps an arbitrary type T with interface that takes index_t.
	
	As methods are only instantiated as needed, the methods in this class
	are optional. If a member function uses "data" in an unsupported way, just dont call it!
*/
template<class index_tpl,class T>
class grid_field{

	T data;

public:

	typedef T data_t;
	typedef index_tpl index_t;
	typedef grid_field<index_t,T> gf_t;

	//ctor
	grid_field(const T& t):data(t){}

	//ctor
	grid_field(){};

	//assignment
	gf_t& operator=(const T& t){
		data=t;
	}

	template<class scalar_t>
	void setConstant(const scalar_t& t){
		data.setConstant(data.rows(),data.cols(),t);
	}
	
	//returns the grid_intex of the end of the sub-grid
	index_t resize(grid_ref grid){
		data.resize(index_t::total(grid));
		return index_t::end(grid);
	}

	//returns the grid_intex of the end of the sub-grid
	index_t resize(grid_ref grid,uint d){
		data.resize(index_t::total(grid),d);		
		return index_t::end(grid);
	}

	T* operator->(){
		return &data;
	}

	T& operator*(){
		return data;
	}
	const T* operator->()const{
		return &data;
	}
	const T& operator*()const {
		return data;
	}

	double& operator[](const index_t i){
		return at(i);
	}
	auto at(const index_t i) -> decltype(T(data)(i.ind)) {
		return data(i.ind);
	}
	auto at(const index_t i) const
		-> decltype((*(new const T(data)))(i.ind)) {
		return data(i.ind);
	}
	auto at2(const index_t i,int j) -> decltype(T(data)(i.ind,j)){
		return data(i.ind,j);
	}

	typename T::RowXpr row(index_t i){
		return data.row(i.ind);
	}

	typename T::ColXpr col(int i){
		return data.col(i);
	}

	void print(ostream& out, grid_ref grid)const {
		for(index_t ind;ind!=ind.end(grid);++ind)
			out<<boost::lexical_cast<string>(at(ind))<<"\n";
		out.flush();
	}
};

//##########################################
/** Scalar and Vector fields on the (spacial) grid
*/
template<int sg,int dim=1,class T = double, int Majority=ColMajor>
struct discretized_space_field : 
		grid_field<grid_index<sg>,
							 Matrix<T,Dynamic,dim,Majority> >{
	string meta_info;
	
	//ctor
	discretized_space_field(){};

	discretized_space_field(string meta_info)
		:meta_info(meta_info)
		{
			// cout<<"CONSTRUCTOR: "<<meta_info<<endl;
		}
	//dtor
	~discretized_space_field()
		{
			// if(meta_info.length()>0)
			// 	cout<<"DESTRUCTOR: "<<meta_info<<endl;
		}
	virtual void i_want_to_be_polymorphic(){}
};

		/** \brief a discretized_space_field with another field for
				additional information
		 */
template<class T,int sg=J_SGS,int dim=1,class K=double, int Majority=ColMajor>
struct discretized_space_field_plus : 
		discretized_space_field<sg,dim,K,Majority> {
	
	using index_t = typename discretized_space_field_plus::index_t;

	//ctor
	discretized_space_field_plus(string meta_info)
		:discretized_space_field<sg,dim,K,Majority> (meta_info){}
	

	index_t resize(grid_ref grid){
		additional.resize(grid);
		return discretized_space_field<sg,dim,K,Majority>::resize(grid);
	}

	grid_field<grid_index<sg>,
						 Matrix<T,Dynamic,dim,Majority> > additional;
};

	 


/** Scalar field on the boundary of the (spacial) grid
*/
template<int sg>
struct boundary_field_t : grid_field<boundary_index<sg-1>,VectorXd>{
	boundary_field_t(string meta_info=""){}
};


template<class row_index_tpl, class col_index_tpl, int Majority=ColMajor>
class sparse_field{
public:
	typedef row_index_tpl row_index_t;
	typedef col_index_tpl col_index_t;

	typedef sparse_field<row_index_t,col_index_t,Majority> my_t;
	typedef SparseMatrix<double,Majority> sm_t;
private:
	sm_t sm;

public:

	//Ctor
	//sparse_field(){};

	my_t& resize(grid_ref g){
		sm.resize(row_index_t::total(g),col_index_t::total(g));
		return *this;
	}

	double& coeffRef(row_index_t row, col_index_t col)
	{
		return sm.coeffRef(*row,*col);
	}

	double coeff(row_index_t row, col_index_t col) const
	{
		return sm.coeff(*row,*col);
	}

	template<class T>
	double row_times_dense(row_index_t row,const T& v) const
	{
		static_assert(Majority==RowMajor,"Majority!=RowMajor");
		
		double r=0;
		for (typename sm_t::InnerIterator it(sm,*row); it; ++it)
			r+= it.value()*(*v)(it.index());
		return r;
	}

	typename sm_t::InnerIterator row_it(row_index_t row)const{
		static_assert(Majority==RowMajor,"Majority!=RowMajor");
		return sm_t::InnerIterator(sm,*row);
	}

	double& insert(row_index_t row, col_index_t col){
		return sm.insert(*row,*col);
	}

	bool insert(row_index_t row, col_index_t col,double w){
		if(w==0.) return false;
		insert(row,col)=w;
		return true; 
	}

	sm_t& operator*(){
		return sm;
	}
	const sm_t& operator*() const{
		return sm;
	}
	sm_t* operator->(){
		return &sm;
	}
	const sm_t* operator->() const{
		return &sm;
	}

	void print() const{
	for (int k=0; k<sm.outerSize(); ++k)
		for (typename sm_t::InnerIterator it(sm,k); it; ++it)
		{
			cout<<it.row()<<" "<<it.col()<<": "<<it.value()<<endl;
			//it.index(); 
		}
	}

	


};


}

	using _grid::sparse_field;
	using _grid::discretized_space_field;
	using _grid::discretized_space_field_plus;
	using _grid::grid_field;
	using _grid::boundary_field_t;
}
