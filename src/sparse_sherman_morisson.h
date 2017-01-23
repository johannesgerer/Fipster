#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/throw_exception.hpp>
#include <limits>
#include <type_traits>
#include <cstdlib>	 //rand
#include <memory>
#include <array>
#endif

#include "Eigen.h"
#include "utility.h"
#include "Math.h"
#include "exceptions.h"

namespace fipster { namespace _sparse_sherman_morisson {

using namespace Eigen;
using namespace std;


template<int Inverse_Majority>
class sparse_sherman_morisson;

/**

Sparse Entry Wise Sherman-Morrison Inverter 

This method is suitable for very sparse A, where the inverse
itself is also very sparse. (For example the inversion of boundary conditions).
The algorithm starts with the inverse of the diagonal part of A.
It then iterates through every off-diagonal entry of A and internally updates the inverse
to reflect this change.

For storage of the sparse representations the library Eigen is used.

The trick to yield good performance is, that the resulting matrix is stored twice: One copy using row major
and another one using column major storage.

For further notes and formulas see "Notizen Scan/11-10-25 Sparse BCs.jpg"
\returns the residuum norm if check was set to true

\param A contains the matrix to be inverted 
\param A1 on exit contains the inverse

*/
template<int A_Majority,int A1_Majority>
double invert_using_single_entry_sparse_sherman_morisson(
	SparseMatrix<double,A1_Majority>& A1,
	const SparseMatrix<double,A_Majority>& A,
	bool check=false){	

	return sparse_sherman_morisson<A1_Majority>::
		create_single_entry_inverse(A1,A,check);		
}

/**

Sparse Row Wise Sherman-Morrison Inverter 

This method is suitable for very sparse A, where the inverse
itself is also very sparse. (For example the inversion of boundary conditions).
The algorithm starts with the inverse of the diagonal part of A.
It then iterates through row of A (without its diagonal) and internally updates the inverse
to reflect this change.

For storage of the sparse representations the library Eigen is used.

The trick to yield good performance is, that the resulting matrix is stored twice: One copy using row major
and another one using column major storage.

For further notes and formulas see "Notizen Scan/11-10-25 Sparse BCs.jpg"
\returns the residuum norm if check was set to true

\param A contains the matrix to be inverted 
\param A1 on exit contains the inverse

*/
template<int A_Majority,int A1_Majority>
double invert_using_rowwise_sparse_sherman_morisson(
	SparseMatrix<double,A1_Majority>& A1,
	SparseMatrix<double,A_Majority>& A,
	bool check=false){	

	return sparse_sherman_morisson<A1_Majority>::
		create_inverse_and_finalize(A1,A,check);		
}

/**

Template specialization if the storage order of the original
and inverse matrix coincide. In this case, a copy of the original
matrix into different storage order has to be performed upfront.

*/
template<int A1_Majority>
double invert_using_rowwise_sparse_sherman_morisson(
	SparseMatrix<double,A1_Majority>& A1,
	SparseMatrix<double,A1_Majority>& A,
	bool check=false){	

	typedef sparse_sherman_morisson<A1_Majority> ssm_t;

	typename ssm_t::A_t Astor=A;//COPY
	return ssm_t::create_inverse_and_finalize(A1,Astor,check);		
}


template<int Inverse_Majority>
class sparse_sherman_morisson
{
public:
	const static int A_Majority = Inverse_Majority==RowMajor?ColMajor:RowMajor;
	const static int A1_Majority = Inverse_Majority;

	typedef SparseMatrix<double,Inverse_Majority> A1_t;
	typedef SparseMatrix<double,A_Majority> A_t;

	typedef SparseMatrix<double,ColMajor> cms_t;
	typedef SparseMatrix<double,RowMajor> rms_t;

	typedef typename A1_t::InnerIterator A1_inner_it_t;
	typedef typename A_t::InnerIterator A_inner_it_t;
	typedef SparseVector<double>::InnerIterator vec_it_t;

	typedef SparseVector<double> tmp_vec;
	typedef array<tmp_vec,2> tmp_vecs;
	
	template<class A_t_e>
	static double create_single_entry_inverse(A1_t& A1, const A_t_e& A,
										bool check=false){	
		
		
		int n=A.rows();
		FIPSTER_ASSERT(A.rows()==A.cols());

		//resize (also zero-initialize) all storages (including external)
		tmp_vecs tmpvecs;
		tmpvecs[0].resize(n);
		tmpvecs[1].resize(n);

		//resize and set zero 
		A1.resize(n,n);
		A1_t A1_T(n,n);

		//if(0){
		//	A.makeCompressed();
		//	A.prune(0.);//BUG 395 don't use prune on uncompressed matrices
		//}

		//create trivial inverse of diagonal matrix and remove it from A
		for(int i=0;i<n;++i){
			double dd=A.coeff(i,i);
			A1.insert(i,i) = 1./dd;
			A1_T.insert(i,i) = 1./dd;
		}

		for (int k=0; k<n; ++k){
			for (typename A_t_e::InnerIterator it(A,k); it; ++it)
				if(it.row()!=it.col())
					update_inverse(A1,A1_T,tmpvecs,it.row(),it.col(),it.value());
		}

		A1.makeCompressed();
		A1.prune(0.);//BUG 395 don't use prune on uncompressed matrices

		if(check) 
			return inversion_residuum(A1,A);
		else
			return numeric_limits<double>().infinity();
	}


	static double create_inverse_and_finalize(A1_t& A1,A_t& A,
										bool check=false){	
		
		
		int n=A.rows();
		FIPSTER_ASSERT(A.rows()==A.cols());

		//resize (also zero-initialize) all storages (including external)
		tmp_vecs tmpvecs;
		tmpvecs[0].resize(n);
		tmpvecs[1].resize(n);

		//resize and set zero 
		A1.resize(n,n);
		A1_t A1_T(n,n);

		
		VectorXd diag(n);

		if(0){
			A.makeCompressed();
			A.prune(0.);//BUG 395 don't use prune on uncompressed matrices
		}

		//create trivial inverse of diagonal matrix and remove it from A
		for(int i=0;i<n;++i){
			double dd=A.coeff(i,i);
			diag[i]=dd;
			A.coeffRef(i,i)=0;
			A1.insert(i,i) = 1./dd;
			A1_T.insert(i,i) = 1./dd;
		}

		//A.makeCompressed();
		//A.prune(0.);//BUG 395 don't use prune on uncompressed matrices
		
		VectorXd vd(VectorXd::Zero(n));
		VectorXd wd(VectorXd::Zero(n));


		for (int k=0; k<n; ++k){
			
			if(0)update_inverse_rowwise2(tmpvecs,k,A1,A,n);
			
			else if(1) update_inverse_rowwise(tmpvecs,wd,vd,k,A1,A,n);
			
			//entry wise:
			for (A_inner_it_t it(A,k); it; ++it)
				update_inverse(A1,A1_T,tmpvecs,it.row(),it.col(),it.value());
		}

		A1.makeCompressed();
		A1.prune(0.);//BUG 395 don't use prune on uncompressed matrices

		
		//restore diagonal of A
		for(int i=0;i<n;++i){
			A.coeffRef(i,i)=diag[i];
		}

		if(check) 
			return inversion_residuum(A1,A);
		else
			return numeric_limits<double>().infinity();
	}


private:

	//Helper sparse_dot_product ########################################################
	template<class T1, class T2>
	FORCE_INLINE static double sparse_dot_product(T1 a,T2 b){
		double result=0;

		while(a && b){
			if(a.index()==b.index()){
				result += a.value()*b.value();
				++a;++b;
			}else if(a.index()<b.index())
				++a;
			else
				++b;
		}
		return result;
	}


	//gcc: war: __declspec(noinline) keine ahnung warum
	template<bool one>
	FORCE_INLINE static void fill_temp2(tmp_vec& w, A1_t& A1,A_inner_it_t v,double la, int n){
		double value;
		for(int i=0;i<n;i++)
			if( (value = sparse_dot_product(A1_inner_it_t(A1,i),v)) !=0){
				if(one) w.coeffRef(i) = value;
				else w.coeffRef(i) = value*la;
			}
	}

	

	//Core function ########################################################
	//Performs an update of the inverse matrix, corresponding
	//to the change in A(row_or_col)+=row_or_col without diagonal	
	static void update_inverse_rowwise(tmp_vecs& tmpvecs,VectorXd& wd, VectorXd&vd, int row_or_col,A1_t& A1,A_t& A,int n){
		
		A_inner_it_t v(A,row_or_col);

		if(v){

			//copy into dense vector
			vd.setZero();
			for(auto v2=v;v2;++v2)
				vd[v2.index()]=v2.value();
			
			auto& z = tmpvecs[0];
			z = A1.innerVector(row_or_col);

			vec_it_t z_it(z);
		
			double la = 1.0/(1+sparse_dot_product(z_it,v));

			bool transpose= A1_Majority==ColMajor;

			// w = Transpose(Inverse(A))*v
			if(transpose)
				wd=A1.transpose()*vd;
			else
				wd=A1*vd;

			//A1 -= (z*w.transpose());
			for(int j = 0 ;j<n ; ++j)
				if(wd[j]!=0.0){					
					double w2=la*wd[j];
					for(auto i = z_it;i;++i)
						if(transpose)
							A1.coeffRef(i.index(),j) -= i.value()*w2;
						else
							A1.coeffRef(j,i.index()) -= i.value()*w2;
				}
		}
	}
	//Core function ########################################################
	//Performs an update of the inverse matrix, corresponding
	//to the change in A(row_or_col)+=row_or_col without diagonal	
	static void update_inverse_rowwise2(tmp_vecs& tmpvecs, int row_or_col,A1_t& A1,A_t& A,int n){
		
		A_inner_it_t v(A,row_or_col);

		if(v){

			auto& z = tmpvecs[0];
			auto& w = tmpvecs[1];
			w.setZero();
			z = A1.innerVector(row_or_col);

			vec_it_t z_it(z);
		
		
			double la = 1.0/(1+sparse_dot_product(z_it,v));

			// w = Transpose(Inverse(A))*v
			if(la==1.0)
				fill_temp2<true>(w,A1,v,la,n);
			else
				fill_temp2<false>(w,A1,v,la,n);
		
			vec_it_t w_it(w);

			//A1 -= (z*w.transpose());
			for(auto i = z_it;i;++i)
				for(auto j = w_it ;j ; ++j)
					if(Inverse_Majority==ColMajor)
						A1.coeffRef(i.index(),j.index()) -= i.value()*j.value();
					else
						A1.coeffRef(j.index(),i.index()) -= i.value()*j.value();
		}
	}
	

	//Performs an update of the inverse matrix, corresponding
	//to the change in A(row,col)+=val
	static void update_inverse(A1_t& A1,A1_t& A1_T,tmp_vecs& tmpvecs,int row, int col, double val)
	{
		if(val==0.0) return;

		auto& colVec = tmpvecs[0];
		auto& rowVec = tmpvecs[1];

		//Get copies of the old inverse' entries:

		if(A1_Majority==ColMajor){
			colVec=A1.innerVector(row);
			rowVec=A1_T.innerVector(col);
		}else{
			colVec=A1_T.innerVector(row);
			rowVec=A1.innerVector(col);
		}

		//The old inverse element at (row,col)
		double Aminus1nm = rowVec.coeff(row); //or equally: =colVec.coeff(col)
		
		//calculate factor:
		double c=  -val/(1.+val*Aminus1nm);
		
		for(vec_it_t i(colVec) ;i;++i)
		{
			double t=c*i.value();
			for(vec_it_t j(rowVec);j ; ++j){
				double t2=t*j.value();
				A1.coeffRef(i.index(),j.index()) += t2;
				A1_T.coeffRef(j.index(),i.index()) += t2;	
			}
		}

	}

	/*

	//Core function ########################################################
	//Performs an update of the inverse matrix, corresponding
	//to the change in A(row_or_col)+=row_or_col without diagonal	
	void update_inverse_rowwise2(int row_or_col,A1_t& A1)
	{
		auto& w = tempA;
		auto& z = tempA1;

		z = A1.innerVector(row_or_col);

		double la = 1.0/(1+sparse_dot_product(A1_inner_it_t(z,0),
												A_inner_it_t(A,row_or_col)));

		//This is slower than the other version:

		if(A_Majority==ColMajor){
			w.resize(n,1);
			w = (A1*A.innerVector(row_or_col)).pruned()*la;
			A1 -= A1_t(w*z); 
		}else{
			w.resize(1,n);
			w = (A.innerVector(row_or_col)*A1).pruned()*la;
			A1 -= A1_t(z*w);
		}
	
	}


	*/


};



//Testing function ########################################################
template<int A1,int A>
double  test_helper(int n){
	SparseMatrix<double,A1> r(n,n);
	SparseMatrix<double,A> a(n,n);

	srand( (unsigned)time( NULL ));
	
	//a= Matrix<double,Dynamic,Dynamic>::Random(n,n);

	//fill
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++){
			double aa=(double)rand()-RAND_MAX;
			a.insert(i,j)=aa;
		}


	cout<<"rowwise:"<<endl;
	double rnom2=invert_using_sparse_sherman_morisson(r,a,true);
	cout<<"internal residuum "<<A1<<" "<<A<<" : "<<rnom2<<endl;	

	double rnorm=inversion_residuum(r,a);
	cout<<"residuum "<<A1<<" "<<A<<" : "<<rnorm<<endl;

	cout<<"single entries:"<<endl;
	rnom2=invert_using_single_entry_sparse_sherman_morisson(r,a,true);
	cout<<"internal residuum "<<A1<<" "<<A<<" : "<<rnom2<<endl;	

	//rnorm=inversion_residuum(r,a);
	//cout<<"residuum "<<A1<<" "<<A<<" : "<<rnorm<<endl;

	return rnorm;
}

template<class T>
double test_sparse_sherman_morisson(T n=10)
{
	double rnorm=
	test_helper<ColMajor,RowMajor>(n);
	rnorm+=test_helper<RowMajor,RowMajor>(n);
	rnorm+=test_helper<RowMajor,ColMajor>(n);
	rnorm+=test_helper<ColMajor,ColMajor>(n);

	cout<<rnorm<<endl;
	return rnorm;
}

}
using _sparse_sherman_morisson::invert_using_rowwise_sparse_sherman_morisson;
using _sparse_sherman_morisson::invert_using_single_entry_sparse_sherman_morisson;
}
