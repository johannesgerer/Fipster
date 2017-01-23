#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <tbb/combinable.h>
#endif

/** singleton Thread Local Storage. Use with
 someClass& c = TLS<someClass>::local();
 DO NOT USE!
	\remark: this does no work to well. Es the combinable object should 
	be created (somehow) in sequential code before the parallel stuff. 
	(e.g. everything broke when calling this in a parallel_reduce body)
	Also the signleton design does not make sense, because it basically
	is a global map of values to thread local storages. What I my application needs to 
	different int's per thread?

	*/
/*
namespace fipster { namespace _TLS {

	using namespace tbb;

	template<class T>
	struct TLS {
		static combinable<T>& get(){
			static combinable<T> cTLS;
			return cTLS;
		}
		static T& local(bool& exists){
			return get().local(exists);
		}
		static T& local(){
			return get().local();
		}
	};	

}

using _TLS::TLS;

}
*/

//################################
//Using BOOST:
//################################
//#ifndef AUTOMATIC_PRECOMPILATION
//#include <boost/thread.hpp>
//#endif
//
////############   Thread local storage  ################
//
//// This class wraps a "thread_specific_ptr" to a storage object
//// and adds to the "get" member function the possibility to 
//// specify the storage size.
//
//template<class Storage> class TLS : private boost::thread_specific_ptr<Storage>{
//public: 
//	Storage* get(int n)
//	{
//		Storage* storage = thread_specific_ptr<Storage>::get();
//
//		if(storage==NULL)
//		{
//			storage = new Storage;
//			reset(storage);
//		}
//
//		storage->set(n);
//		return storage;
//	};
//};
