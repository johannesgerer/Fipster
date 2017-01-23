#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <utility>
#include <string>
#include <iostream>
#include <sstream> 
#include <list>
#include <tbb/spin_mutex.h>
#include <algorithm>
#include <cstdint>
#include <exception>
#include <comphash.h>
#endif


#ifdef _MSC_FULL_VER
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

/** \todo make sure, this is used in the correct way: 
// all tests, that cannot fail in production, that should
// be tested automatically, should be conditional on FIPSTER_NDEBUG
// All tests, that have to performed every time (or depending on runtime)
 must not use this!*/
//see also precompiled.h
#ifndef DEBUG
#define FIPSTER_NDEBUG 1
#else
#define FIPSTER_NDEBUG 0
#endif


//freopen( "c:\\test.dat" ,"w",stdout);

namespace fipster {

	using namespace std;

	//check if all elements in container are unique
	template<class T>
	auto rearrange_and_get_number_of_uniques(T& t)
		-> decltype(begin(t)-end(t))
	{
		sort(begin(t),end(t));
		return unique(begin(t),end(t)) - begin(t);
	}


	//Free reserved storage
	template<class T>
	void free_reserved(T& t)
	{
		T temp;
		swap(t,temp);
	}

	//this is used for default template parameters
	struct nothing : comphashable<nothing> {
		int empty() const{ return 0; }
		template<class T> FORCE_INLINE void my_combine(T& t)const{
			t(&nothing::empty);
		}
	};


/** \brief Business time data type (INTEGER!!!!!!!)
	The length of one tick of btime is given by years_per_btick
	Needs to be signed to be able to safely calculate differences (e.g. with offsets or max/min times) 
*/
typedef int64_t btime;



class Field;

//key type for effective boundary values (that depend on the time and an underlying field)
typedef pair<btime,Field*> eff_bv_key_t;

//Different indices representing different array linearizations on the grid
typedef int arrayoffset;
typedef int interiorarrayoffset;
typedef int boundaryarrayoffset;

// a = max(a,b)
template<typename A,typename B>void applymax(A& a, B& b){ if(a<b) a=b; }
template<typename A,typename B>void applymin(A& a, B& b){ if(a>b) a=b; }

////STL Pair output
//template<class a,class b>ostream& operator<<(ostream& out,const pair<a,b>& A)
//{
//	out<<"("<<A.first<<", "<<A.second<<")";
//	return out;
//}
//
//template<typename type> string toString(const type& a)
//{
//	ostringstream s;
//	s<<a;
//	return s.str();
//}



/*This routine splits the the work, given 
!as a number of independent tasks, that are numbered from 
!1 to "work" amoung a number of threads "n".
!The result is an array workLoads, which contains the upper bounds 
!of the calcuation and should be used like this:
!The i-th thread execute the numbers 
!from workLoads(i)+1 to workLoads(i+1)*/
/*TODO? inline void SplitWorkEqually(int work, int1& workloads, int threads)
{
	workloads.reset(threads+1);
	workloads(0)=0;
	int mod=work%threads;

	for(int i=1; i <= threads; ++i)
		workloads(i) = workloads(i-1)+work/threads+(i<=mod?1:0);
		
}*/



//NOT WORKING
//template<typename T> class Clone{
//public:
//    virtual T* clone() { 
//		return new T(static_cast<const T&>(*this)); 
//	}
//	/*virtual T* clone() { 
//		return new T(*(static_cast<T*>(this))); 
//	}*/
//};
//

//This function performs a cast from typeof(source) to typeof(target)
//The cast is dynamic and checks for error in Debug Mode and static, if FIPSTER_NDEBUG==1
template< class newtype,class oldtype> void my_cast(oldtype *source,newtype*& target)
{	
	if(FIPSTER_NDEBUG)
		target = static_cast<newtype*>(source);
	else{
/* from  http://msdn.microsoft.com/en-us/library/c36yw7x9.aspx

	 In general you use static_cast when you want to convert numeric data types such as enums to ints or ints to floats, and you are certain of the data types involved in the conversion. static_cast conversions are not as safe as dynamic_cast conversions, because static_cast does no run-time type check, while dynamic_cast does. A dynamic_cast to an ambiguous pointer will fail, while a static_cast returns as if nothing were wrong; this can be dangerous. Although dynamic_cast conversions are safer, dynamic_cast only works on pointers or references, and the run-time type check is an overhead. For more information, see dynamic_cast Operator.*/
		try{
			target = dynamic_cast<newtype*>(source);
    }catch(std::bad_cast& e){
			cout<<"Error in dynamic_cast:  bad_cast: "<< e.what()  <<endl;
 			assert(false);  // __debugbreak(); GCC
			throw &e;
		}
#ifdef _MSC_FULL_VER
    catch(std::__non_rtti_object& e){
			cout<<"Error in dynamic_cast: __non_rtti_object "<<endl;
 			assert(false);  // __debugbreak(); GCC
			throw &e;
		}
#endif

		if( target==nullptr)
		{
			cout<<"Error in dynamic_cast"<<endl;
 			assert(false);  // __debugbreak();
			throw new bad_cast();
		}
	}
};



inline bool is_xml_special(string tag)
{
	return (tag == "<xmlattr>" || 
			tag == "<xmlcomment>" || 
			tag == "<xmltext>");
}

	/**
	if a does not contain b, b is appended to a (using push_back).
	the function returns the index of b in a.
	*/
	template<class A,class B>
	auto get_index(A& a,B b) -> decltype(end(a)-begin(a))
	{
		auto i=find(begin(a),end(a),b);
		if(i==end(a)){
			a.push_back(b);
			return (--end(a))-begin(a);
		}
		return i-begin(a);
	}

template<class A, class B,class C>
auto betterAt(A& a, const B& b, const C& c) -> decltype(a.at(b)) {
	try{ return a.at(b); }
	catch(out_of_range& e){
		BOOST_THROW_EXCEPTION(out_of_range
													(string(e.what())+"\nKey: "+toS(b)+"\nContext: "+toS(c)));
	}
}

}



/* bad practice:
#define fo(a,b,c) for( a = ( b ); a < ( c ); ++ a )
#define fr(a,b) fo( a, 0, ( b ) )
#define fi(a) fr( i, ( a ) )
#define fj(a) fr( j, ( a ) )
#define fk(a) fr( k, ( a ) )
#define fui(a) for(uint i = 0; i<( a ); ++ i )
#define fii1(a) for(int i = 1; i<=( a ); ++ i )
#define fii(a) for(int i = 0; i<( a ); ++ i )
#define fij(a) for(int j = 0; j<( a ); ++ j )
*/
