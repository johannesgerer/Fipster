#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <algorithm>
#include <cmath>
#include <boost/math/distributions/normal.hpp>
#include <exception>
#endif

#include "integer_types.h"

namespace fipster {


//integer power function
template<class A> A pow2(A a,uint b)
{
	A r=1;
	for(uint i = 0; i<b ; ++i )  r*=a;
	return r;
}


/* \brief Solves the tridiagonal system M.r=u

	Where M is tridiagonal and
	M_(i,i-1) = a_i, M_ii = b_i, M_(i,i+1) = c_i
	NOTE: c(n-1) and b(0) are not used 

	gam(0...n-1) is workspace

	Code taken from "numerical recipes":
	"Solves for a vector u[0..n-1] the tridiagonal linear set given by equation (2.4.1). a[0..n-1],
	b[0..n-1], c[0..n-1], and r[0..n-1] are input vectors and are not modified."
	(2.4.1):  a - sub diagonal, b - diagonal, c - superdiagonal
	"Notice that a0 and cN1 are undefined and are not referenced by the routine that
	follows."
template<class T>
void tridag(const T &a,const T &b,const T &c, const T &r, T &u,unsigned int n, T &gam)
{
typedef double Doub;
typedef int Int;
Int j;
Doub bet;
if (b[0] == 0.0) throw("Error 1 in tridag");
If this happens, then you should rewrite your equations as a set of order N  1, with u1
trivially eliminated.
u[0]=r[0]/(bet=b[0]);
for (j=1;j<n;j++) { Decomposition and forward substitution.
gam[j]=c[j-1]/bet;
bet=b[j]-a[j]*gam[j];
2.4 Tridiagonal and Band-Diagonal Systems of Equations 57
if (bet == 0.0) throw("Error 2 in tridag"); Algorithm fails; see below.
u[j]=(r[j]-a[j]*u[j-1])/bet;
}
for (j=(n-2);j>=0;j--)
u[j] -= gam[j+1]*u[j+1]; Backsubstitution.
}
	*/

template<class T>
double geometric_average(const T& v,uint n){
	double r=1;

	for(uint i=0;i<n;i++)
		r*= v[i];

	return n==1 ? r : pow(r,1.0/n);
}


/** a collection of interpolation methods
*/
struct interpolate{
	/** 1D interpolation: get the interpolated value f(c) of the function
	f given by points f(*(x+i))=y[offset+i] for x <= x+i < e
	x has to be sorted in ascending order
	*/
	template<class A,class B,class dB>
	static double oneD(const A& x,const A& e,const B& y,dB offset, double c){
		auto d=equal_range(x,e,c);
		//if exact match exists
		if(d.first < d.second)
			return y[offset+d.first-x];
		else{ //interpolate if no exact entry was found
			if(!(d.first >  x && d.first < e))
				BOOST_THROW_EXCEPTION(std::range_error("interpolation outside of x"));
			//http://de.wikipedia.org/wiki/Interpolation_(Mathematik)#Lineare_Interpolation
			auto h=d.first-1;
			auto f=h-x+offset;		
			return y[f]+(y[f+1]-y[f])/(*(h+1)-*h)*(c-*h);
		}
		return c; // return something
	}	
	/** this method tests the implemented logic */
	static void test();
};

//Returns Ln[Gamma[xx]]
double gammln(const double xx);

//Returns n! as double
double factorial(const unsigned int n);

//Returns Ln[n!]
double factln(const int n);

//Returns the binomial coefficient N over K as double
double binom(const int n, const int k);

//COMP_NEXT computes the compositions of the integer N into K parts.
void comp_next ( int n, int k, int a[], bool *more, int *h, int *t );

//round to nearest integer represented as double
double round(double r);

/** calculates the relative difference of two floating point numbers
 abs_if_zero: use abs difference, if b is zero
 */
double reldiff(double a, double b, bool abs_if_zero);

template<class A,class B>
double maxreldiff( A& a, B &b )
{ 
	if(a.size()!=b.size())
		return -1.0;

	double r=0;
	for(int i = 0; i<a.size() ; ++i ) 
		r=max(r,reldiff(a[i],b[i]));
	return r;
}

//double maxreldiff(double2& a, double2 &b);



/** calculate the residuum of a sparse matrix inversion: 
	\returns Norm(A*A1-1)+Norm(A1*A-1)
*/
template<class A1_t,class A_t>
double inversion_residuum(A1_t& A1,A_t& A){
	
	A1_t C = A1*A;
	A1_t D = A*A1;

	for(int i = 0; i<A.rows() ; ++i  ){
		C.coeffRef(i,i)-=1.0;
		D.coeffRef(i,i)-=1.0;
	}
	
	if(0){

		C.prune(0.);
		// cout<<A<<endl<<endl;
		// cout<<A1<<endl<<endl;
		// cout<<C<<endl;
		
	}

	return (C.norm()+D.norm());///min(A1.norm(),A.norm());
}

	

}
