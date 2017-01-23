#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <string>
#include <map>
#include <iostream>
#include <tbb/spin_mutex.h>
#include <algorithm>
#include <tbb/compat/thread>
#include <vector>
#include <sstream>
#endif

namespace fipster { namespace _thread_logger {

	using namespace std;
	using tbb::spin_mutex;

	/** use like this: thread_logger()<<"blabla"<<...<<endl;
	*/
struct thread_logger{

	
	static spin_mutex& get_mutex();

	//tbb::spin_mutex::scoped_lock lock;
	ostringstream buffer;
	
	/**ctor*/
	thread_logger();

	typedef decltype(tbb::this_tbb_thread::get_id()) thread_id_t;

	const static unsigned int min_size;
	
	~thread_logger();

	template<typename T>
	ostream& operator<<(const T &value)
	{
		buffer<<value;return buffer;	
	}

	static string get_current_name();
};

}

using _thread_logger::thread_logger;

}
