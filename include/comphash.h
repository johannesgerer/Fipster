/*
@@@@@@@@@ THIS FILE IS OBSOLETE @@@@@@@@@@

Please use 

https://github.com/johannesgerer/Enhance

instead!

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 */
#pragma once

#include <boost/functional/hash.hpp>
#include <array>
#include <boost/optional.hpp>

#include <type_traits>
#include "pretty_printer.h"

#ifndef FORCE_INLINE
#ifdef _MSC_FULL_VER
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif
#endif

/** \brief external base classes
		
"c++ enhance" zero-cost c++ class enhancer    

SumItUp - Zero over-head automatic derivation of hash, comparison,
serialization, pretty printing

		base classes, that can be used to inherit predefined (working) 
		hash_value or comparison_combiner operator interfaces
*/
namespace comphash {


			//plain data-member overload
	template<class A1,class D>
			FORCE_INLINE const A1& access(A1(D::*a1),const D& x)
				//-> decltype((x.*a1)) 
				//decltype(()): The inner parentheses cause the statement to be 
				//evaluated as an expression instead of a member access. 
				//And because a is declared as a const pointer, the type is a
				//reference to const double.
				//see http://msdn.microsoft.com/de-de/library/dd537655.aspx
			{
				return x.*a1;
			}
			
			//nullary  function-(const!!!)member overload
	template<class A1,class D>
			FORCE_INLINE A1 access(A1(D::*a1)()const,const D& x)
				//	-> decltype((x.*a1)())
			{
				return (x.*a1)();
			}

			//operator overload
			template<class A1,class D>
			FORCE_INLINE auto access(A1 a1,const D& x)
				-> decltype(a1(x)){
				return a1(x);
			}

		template<class Accessor>
		struct CombineByValue{
			Accessor m;
			CombineByValue(Accessor m):m(m){};
			template<class D>
			auto operator()(const D& d) const
				->decltype(*(d.*m)){ 
					//there are no pointers to references, so reference qualifiers are not important here
				return *(d.*m);
			}
		};
	/** \brief Utility wrapper for pointer members, that should by accessed by value
	 */
	template<class Accessor>
	comphash::CombineByValue<Accessor> byValue(Accessor a){
		return comphash::CombineByValue<Accessor>(a);
	}

	template<class A, class B>
		struct Range{
		A a;
		B b;
			Range(A a,B b):a(a),b(b){};
		};
	template<class A,class B>
	Range<A,B> range(A a,B b){
		return Range<A,B>(a,b);
	}
	

	// the following would make 'access' ambiguous'. it would be possible with another wrapper.
	// //values
	// template<class A1,class D>
	// FORCE_INLINE A1 access(A1 a1,const D& x)
	// 	{
	// 	return a1;
	// }

	enum { LessEnum, EqualEnum, GreaterEnum };

	namespace internal{
		template<class kernel_t,class D,class result_t>
		struct recursion{

			typedef recursion<kernel_t,D,result_t> recursion_base_t;
			const D& x;
			result_t result;

			//ctor
			recursion(const D& x,result_t result=result_t())
				:x(x),result(result){};
			
			//apply kernel
			template<class A1>
			FORCE_INLINE bool is_done_after_kernel(A1 a1)
			{ 
				/*FORCE_INLINE void finalize(){
					overwrite to implement a final transformation of the result
				}
				bool t=static_cast<kernel_t&>(*this).kernel(a1);
				if(t) static_cast<kernel_t&>(*this).finalize();*/
				return static_cast<kernel_t&>(*this).kernel(a1);
			}
						
#ifndef _MSC_FULL_VER
		private:
			FORCE_INLINE void operator()(){}			
		public:
			template<class A, class ... Rest>
			FORCE_INLINE result_t operator()(A a, Rest  ... rest)
			{	if(!is_done_after_kernel(a)) operator()(rest ...); return result; }
#else
			template<class A1>
			FORCE_INLINE result_t operator()(A1 a1)
			{	is_done_after_kernel(a1); return result; }			
			template<class A1,class A2>
			FORCE_INLINE result_t operator()(A1 a1, A2 a2)
			{	if(!is_done_after_kernel(a1)) operator()(a2); return result; }
			template<class A1,class A2,class A3>
			FORCE_INLINE result_t operator()(A1 a1, A2 a2,A3 a3)
			{	if(!is_done_after_kernel(a1)) operator()(a2,a3); return result; }
			template<class A1,class A2,class A3,class A4>
			FORCE_INLINE result_t operator()(A1 a1, A2 a2,A3 a3,A4 a4)
			{	if(!is_done_after_kernel(a1)) operator()(a2,a3,a4); return result; }
			template<class A1,class A2,class A3,class A4,class A5>
			FORCE_INLINE result_t operator()(A1 a1, A2 a2,A3 a3,A4 a4,A5 a5)
			{	if(!is_done_after_kernel(a1)) operator()(a2,a3,a4,a5); return result; }
			template<class A1,class A2,class A3,class A4,class A5,class A6>
			FORCE_INLINE result_t operator()(A1 a1, A2 a2,A3 a3,A4 a4,A5 a5,A6 a6)
			{	if(!is_done_after_kernel(a1)) operator()(a2,a3,a4,a5,a6); return result; }
			template<class A1,class A2,class A3,class A4,class A5,class A6,class A7>
			FORCE_INLINE result_t operator()(A1 a1, A2 a2,A3 a3,A4 a4,A5 a5,A6 a6,A7 a7)
			{	if(!is_done_after_kernel(a1)) operator()(a2,a3,a4,a5,a6,a7); return result; }
			template<class A1,class A2,class A3,class A4,class A5,class A6,class A7,class A8>
			FORCE_INLINE result_t operator()(A1 a1, A2 a2,A3 a3,A4 a4,A5 a5,A6 a6,A7 a7,A8 a8)
			{	if(!is_done_after_kernel(a1)) operator()(a2,a3,a4,a5,a6,a7,a8); return result; }
			template<class A1,class A2,class A3,class A4,class A5,class A6,class A7,class A8,class A9>
			FORCE_INLINE result_t operator()(A1 a1, A2 a2,A3 a3,A4 a4,A5 a5,A6 a6,A7 a7,A8 a8,A9 a9)
			{	if(!is_done_after_kernel(a1)) operator()(a2,a3,a4,a5,a6,a7,a8,a9); return result; }
			template<class A1,class A2,class A3,class A4,class A5,class A6,class A7,class A8,class A9,class A10>
			FORCE_INLINE result_t operator()(A1 a1, A2 a2,A3 a3,A4 a4,A5 a5,A6 a6,A7 a7,A8 a8,A9 a9,A10 a10)
			{	if(!is_done_after_kernel(a1)) operator()(a2,a3,a4,a5,a6,a7,a8,a9,a10); return result; }
#endif


			FORCE_INLINE result_t operator()(void (*p)(kernel_t&)){
				p(static_cast<kernel_t&>(*this));
				return result;
			}

			FORCE_INLINE result_t compute_result(){
				x.my_combine2(static_cast<kernel_t&>(*this));
				return result;
			}

			FORCE_INLINE operator result_t(){
				x.my_combine(static_cast<kernel_t&>(*this));
				return result;
			}
		};
		
		template<class D,char sep>
		struct string_combiner : recursion<string_combiner<D,sep>,D,string>{
						
			FORCE_INLINE string_combiner(const D& x)
				: string_combiner::recursion(x,""){};

			template<class Accessor>
			FORCE_INLINE bool kernel(Accessor a) {
				this->result+=toS(access(a,this->x));
				this->result.push_back(sep);
				return false;//false implies that the recursion does not stop until the last 
							 //member to combine has been handled
			}
		};

		template<class D>
		struct hash_combiner : recursion<hash_combiner<D>,D,size_t>{
						
			FORCE_INLINE hash_combiner(const D& x)
				: hash_combiner::recursion(x,0){};

			template<class Accessor>
			FORCE_INLINE bool kernel(Accessor a) {
        
				boost::hash_combine(this->result,access(a,this->x));
				return false;
				//false implies that the recursion does not stop until the
				//last member to combine has been handled
			}

			template<class A,class B>
			FORCE_INLINE bool kernel(Range<A,B> ac) {
				auto   b = access(ac.a,this->x);
				// r-values for functions and cast to reference for l-values
				auto&& e = access(ac.b,this->x);
				for(; b<e; ++b)
					boost::hash_combine(this->result,*b);
				return false;
			}
			
		};


		template<class comp_t,class D>
		struct comparison_combiner : recursion<comparison_combiner<comp_t,D>,D,bool> {

			const D& y;

			/*this constructor is used, to ensure the best possible
			  RVO, by not requiring to use a named temporary to be passed
			  to a my_combine all, which would require NRVO
			  */
			FORCE_INLINE comparison_combiner(const D& x,const D& y)
				:comparison_combiner::recursion(x),y(y){};

			template<class Accessor>
			FORCE_INLINE bool kernel(Accessor ac) {
				// r-values for functions and cast to reference for l-values
				auto&& a = access(ac,this->x);
				auto&& b = access(ac,this->y);
				return comp_t::go(this->result,a,b);
			}

			template<class A,class B>
			FORCE_INLINE bool kernel(Range<A,B> ac) {
				auto   b_x = access(ac.a,this->x);
				auto   b_y = access(ac.a,this->y);
				// r-values for functions and cast to reference for l-values
				auto&& e_x = access(ac.b,this->x);
				auto&& e_y = access(ac.b,this->y);

				assert(e_x-b_x==e_y-b_y); 
				//only compare ranges of equal length (alternative
				//implementations are possible, esp for ==
				for(; b_x<e_x; ++b_x, ++b_y)
					if(comp_t::go(this->result,*b_x,*b_y))
						return true;
				return false;
			}
		};

	}


	// Hash Interface ####################################################
	//template<class D>
	//std::size_t combined_hash(const D& x){
	//	return internal::hash_combiner<D>(x).seed;
	//}

	//base class
	template<class D>
	struct hashable{
		typedef internal::hash_combiner<D> Hash;
		FORCE_INLINE friend std::size_t hash_value(const D& x){
			return Hash(x);
		}
	};

	// Stringify Interface ####################################################
	
	template<class D,char sep=';'>
	struct stringifyable{
		typedef internal::string_combiner<D,sep> Stringify;
		FORCE_INLINE friend ostream& operator<<(ostream& out,const D& x){
			out<<(string)Stringify(x);
			return out;
		}
	};

	

	// Compare Interface ####################################################

	struct LessOp{
		template<class A>
		static bool go(bool& r, const A& a, const A& b){
			return (r = a < b) || (a != b);
		}
	};
	struct EqualOp{
		template<class A>
		static bool go(bool& r, const A& a, const A& b){
			return !(r = a == b);
		}
	};
	struct GreaterOp{
		template<class A>
		static bool go(bool& r, const A& a, const A& b){
			return (r = a > b) || (a != b);
		}
	};
	struct ComparisonOp{};

	template<class D>
	struct comparable_base{
		typedef internal::comparison_combiner<LessOp,D> Less;
		typedef internal::comparison_combiner<GreaterOp,D> Greater;
		typedef internal::comparison_combiner<EqualOp,D> Equal;
	};

	template<class D,class comp_t=ComparisonOp>
	struct comparable
		:comparable<D,GreaterOp>,comparable<D,EqualOp>,comparable<D,LessOp>
	{};

	//  LESS
	template<class D>
	struct comparable<D,LessOp> : comparable_base<D> {
		FORCE_INLINE bool operator<(const D& y)const {
					return Less(static_cast<const D&>(*this),y);
		}
		//Automatically deduced comparison_combiner operator
		FORCE_INLINE bool operator>=(const D& y)const {
			return !(static_cast<const D&>(*this) < y);
		}
	};

	//  GREATER
	template<class D>
	struct comparable<D,GreaterOp> : comparable_base<D> {
		FORCE_INLINE bool operator>(const D& y)const {
			return Greater(static_cast<const D&>(*this),y);
		}
		//Automatically deduced comparison_combiner operator
		FORCE_INLINE bool operator<=(const D& y)const {
			return !(static_cast<const D&>(*this)>y);
		}
	};

	//  EQUAL
	template<class D>
	struct comparable<D,EqualOp> : comparable_base<D> {
		FORCE_INLINE bool operator==(const D& y)const {
			return Equal(static_cast<const D&>(*this),y);
		}
		//Automatically deduced comparison_combiner operator
		FORCE_INLINE bool operator!=(const D& y)const {
			return !(static_cast<const D&>(*this)==y);
		}
	};

	
	// Compare & Hash Interface ####################################################
	template<class D>
	struct comphashable : comparable<D>,hashable<D>{};

	// Free Compare & Hash combiner factories ####################################################
	// These functions can be used in comparison hash_value functions in classes,
	//that cannot derive from a comphash base class.

	template<class D>
	internal::string_combiner<D,';'> Stringify(const D& x){
		return internal::string_combiner<D,';'>(x);
	}
	

	template<class D>
	internal::hash_combiner<D> Hash(const D& x){
		return internal::hash_combiner<D>(x);
	}
	
	template<class T>
	internal::comparison_combiner<LessOp,T> Less(const T& x,const T& y){
		return internal::comparison_combiner<LessOp,T>(x,y);
	}

	template<class T>
	internal::comparison_combiner<GreaterOp,T> Greater(const T& x,const T& y){
		return internal::comparison_combiner<GreaterOp,T>(x,y);
	}

	template<class T>
	internal::comparison_combiner<EqualOp,T> Equal(const T& x,const T& y){
		return internal::comparison_combiner<EqualOp,T>(x,y);
	}

}


using comphash::comparable;
using comphash::comphashable;
using comphash::stringifyable;
using comphash::hashable;
using comphash::byValue;
using comphash::range;









// not needed anymore //###############################################################
// //##################### Tuple and Array hash (std::tr1) #########

// //hash compile over static range (array or tuple) for use in array and tuple hash
// namespace comphash {

// template<int n>	struct static_hash{
// 	template<class T> static void range(const T& t,size_t& seed)
// 	{
// 		boost::hash_combine(seed,std::get<n-1>(t));
// 		static_hash<n-1>::range(t,seed);
// 	}
// };

// template<>	struct static_hash<0>{
// 	template<class T> static void range(const T& t,size_t& seed)
// 	{}
// };

// }

// //Hash_value for arrays, tuples and shared_ptr,etc ...!!!!!!!!!!!
// namespace std{
// namespace tr1{
	
// //Hash_value for arrays & tuples 
// template<class T>
// std::size_t hash_value(const T& x){
// 	std::size_t seed = 0;
// 	//todo replace with boost::hash_range 
// // http://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#idp53740664-bb
// 	comphash::static_hash<std::tuple_size<T>::value>::range(x,seed);
// 	return seed;
// }
	
// }}

// //Hash_value for shared_ptr
// template<class T>
// std::size_t hash_value(const shared_ptr<T>& x){
// 	return boost::hash<T*>()(x.get());
// }
namespace boost{

	//Hash_value for optional
	template<class T>
	std::size_t hash_value(const boost::optional<T>& x){
		std::size_t seed = 0;
		boost::hash_combine(seed,(bool)x);
		if(x)
			boost::hash_combine(seed,*x);
		return seed;
	}
}


//TODO:
/*
explain http://msdn.microsoft.com/query/dev10.query?appId=Dev10IDEF1&l=DE-DE&k=k(C2898);k(VS.ERRORLIST)&rd=true

https://www.google.com/search?q=c%2B%2B+automatic+derivation+operators&ie=utf-8&oe=utf-8&client=firefox-b#q=c%2B%2B+automatic+derivation+operators+-differentiation


todo?
boost.operators

boost.serialize

https://github.com/philsquared/Catch

*/
