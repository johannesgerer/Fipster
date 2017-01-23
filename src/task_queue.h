#pragma once


#ifndef AUTOMATIC_PRECOMPILATION
//#include <boost/asio.hpp>
#include <functional>
#include <queue>
#include <vector>
#include <list>
#include <algorithm>
#endif



namespace fipster { namespace _task_queue {

	using namespace std;
	 
/** \brief %A queue with tasks to be completed. (Recursion breaker)

This is not thread safe

A THREADED EXTENSION WOULD REQUIRE BOOST::ASIO with io_service and thread_group
OR TBB on a higher level OR hand-written locking

*/
class task_queue
{
public:

	typedef function<void()> functor_t;
	#if 0
	typedef queue<functor_t,list<functor_t>> queue_t;
#else
	typedef queue<function<void()>> queue_t;
#endif
	queue_t q;
	unsigned int m;
	function<void()> next;
	bool is_next;

	//template<class F>
	//static void set_next(F&& f)
	//	get().next = forward<F>(f);
	//	get().is_next = true;
	//}

	//static void do_next(){
	//	get().is_next = false;
	//	get().next();
	//}

	template<class F>
	static void enqueue(F&& f){
		get().q.push(forward<F>(f));
	}

	static void finish_tasks();

	//get singleton
	static task_queue& get(){
		static task_queue q;
		return q;
	}
};

}
using _task_queue::task_queue;

}
