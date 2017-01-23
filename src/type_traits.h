#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <tuple>
#include <vector>
#include <map>
#endif

namespace fipster { namespace _type_traits {

using namespace std;

#if defined(_MSC_VER) && _MSC_VER == 1600
#define MANUAL_TUPLE_ARGUMENT_EXPANSION
#endif

#ifdef MANUAL_TUPLE_ARGUMENT_EXPANSION
#define TTN2 std::tr1::_Nil
#define TTN =TTN2
#define TUPLE_PARAMS \
	typename T0, typename T1, typename T2, typename T3, typename T4, \
	typename T5, typename T6, typename T7, typename T8, typename T9
#define OPTIONAL_TUPLE_PARAMS \
	typename T0 TTN, typename T1 TTN, typename T2 TTN, typename T3 TTN, typename T4 TTN, \
	typename T5 TTN, typename T6 TTN, typename T7 TTN, typename T8 TTN, typename T9 TTN
#define TUPLE_ARGS T0, T1, T2, T3, T4, T5, T6, T7, T8, T9
#define TEN_DEFAULT_TUPLE_ARGS  TTN2,TTN2,TTN2,TTN2,TTN2,TTN2,TTN2,TTN2,TTN2,TTN2
#define NINE_DEFAULT_TUPLE_ARGS , T0,  TTN2,TTN2,TTN2,TTN2,TTN2,TTN2,TTN2,TTN2,TTN2

#else
#define TUPLE_PARAMS typename... Types
#define TUPLE_ARGS Types...
#define OPTIONAL_TUPLE_PARAMS TUPLE_PARAMS
#define NINE_DEFAULT_TUPLE_ARGS T0
#define TEN_DEFAULT_TUPLE_ARGS  
#endif


/** \brief meta function to generate a std::vector<T> with default arguments
*/
template<class T>
struct default_vector_apply{
	typedef vector<T> type;
};

/** \brief meta function to generate for a given K a meta function which generates std::map<K,T> with default arguments for a given T
*/
template<class K>
struct default_map{
	template<class T>
	struct apply{
		typedef map<K,T> type;
	};
};

//############################
// apply_meta
//############################
/** \brief X metafunction
// Applies X<a> to every element a of A. The resulting type can be accessed using
// apply<X,A>::type
*/
template<template <class> class X,class A>
struct apply_meta{
	typedef typename X<A>::type type;
};

/** \brief  Tuple specialization */
template<template <class> class X, TUPLE_PARAMS>
struct apply_meta<X,tuple<TUPLE_ARGS>>
{
	#ifdef MANUAL_TUPLE_ARGUMENT_EXPANSION

	typedef tuple<TUPLE_ARGS> input_tuple;
	static const int n = tuple_size<input_tuple>::value;
	typedef tuple< 
		typename conditional< (0>=n) , T0, typename X<T0>::type >::type,
		typename conditional< (1>=n) , T1, typename X<T1>::type >::type,
		typename conditional< (2>=n) , T2, typename X<T2>::type >::type,
		typename conditional< (3>=n) , T3, typename X<T3>::type >::type,
		typename conditional< (4>=n) , T4, typename X<T4>::type >::type,
		typename conditional< (5>=n) , T5, typename X<T5>::type >::type,
		typename conditional< (6>=n) , T6, typename X<T6>::type >::type,
		typename conditional< (7>=n) , T7, typename X<T7>::type >::type,
		typename conditional< (8>=n) , T8, typename X<T8>::type >::type,
		typename conditional< (9>=n) , T9, typename X<T9>::type >::type
	> type;

	#else
	typedef tuple<typename X<Types>::type...> type;
	#endif
};

/** \brief  Tuple<> specialization */
template<template <class> class X>
struct apply_meta<X,tuple<>>
{
	typedef tuple<> type;
};

//############################
// apply
//############################
/** \brief applies a template to a class: X<A>
	
	\tparam X template
*/
template<template <class> class X,class A>
struct apply{

	template<class B> struct temp{
		typedef X<B> type;
	};

	typedef typename apply_meta<temp,A>::type type;
};

//############################
// apply_meta_apply
//############################
// X template, Y metafunction
template<template <class> class X,template <class> class Y,class A>
struct apply_meta_apply{

	template<class B> struct temp{
		typedef X<typename Y<B>::type> type;
	};

	typedef typename apply_meta<temp,A>::type type;
};


//############################
// apply_apply_meta
//############################
// X metafunction, Y template
template<template <class> class X,template <class> class Y,class A>
struct apply_apply_meta{

	template<class B> struct temp{
		typedef typename X<Y<B>>::type type;
	};

	typedef typename apply_meta<temp,A>::type type;
};

//############################
// apply_apply
//############################
// X template, Y template
template<template <class> class X,template <class> class Y,class A>
struct apply_apply{

	template<class B> struct temp{
		typedef X<Y<B>> type;
	};

	typedef typename apply_meta<temp,A>::type type;
};

}

using _type_traits::apply;
using _type_traits::apply_meta;
using _type_traits::apply_meta_apply;
using _type_traits::apply_apply;
using _type_traits::apply_apply_meta;
using _type_traits::default_vector_apply;
using _type_traits::default_map;

}
