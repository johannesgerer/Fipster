#pragma once


#ifndef AUTOMATIC_PRECOMPILATION
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Core>
#include <comphash.h>
#include <boost/functional/hash.hpp>
#endif

namespace Eigen {

	template<class T>
	const T* begin(const Matrix<T,Dynamic,1>& x){
		return x.data();
	}

	template<class T>
	const T* end(const Matrix<T,Dynamic,1>& x){
		return x.data()+x.size();
	}

	template<class A>
	auto v_r() -> decltype(range(begin<A>,end<A>)) {
		return range(begin<A>,end<A>);
	}
	
	//more expandable: t is variadic
	template<class A,class T>
	void combine(T& t){
		t(range(begin<A>,end<A>));
	}

	template<class T>
	std::size_t hash_value(const Matrix<T,Dynamic,1>& x){
		// return comphash::Hash(x)(combine<T>);
		return comphash::Hash(x)(v_r<T>());
	}
 
	template<class T>
	bool operator==(const Matrix<T,Dynamic,1>&& x,const  Matrix<T,Dynamic,1>&& y){
		return comphash::Equal(x,y)(combine<T>);
	}
 
}

namespace fipster {
	using Eigen::RowMajor;
	using Eigen::ColMajor;
}
