#include "precompiled.h"
#include "logging.h"

namespace fipster {
	namespace _thread_logger {

		const unsigned int thread_logger::min_size = 6;

	/**ctor*/
	thread_logger::thread_logger()
		//:lock(get_mutex())
	{
		operator<<(get_current_name())<<
			//" ("<<tbb::this_tbb_thread::get_id()<<
			": ";
	}

	string thread_logger::get_current_name(){
		static vector<thread_id_t> ids;
		static vector<string> names;
		static size_t s=0;

		auto id = tbb::this_tbb_thread::get_id();
		
		size_t i;
		for(i=0;i<s;i++)
			if(ids[i]==id)
				return names[i];
		
		{
			tbb::spin_mutex::scoped_lock lock(get_mutex());
			ids.push_back(id);
			string r(max((size_t)min_size,ids.size()),'.');
			r[i]='x';
			names.push_back(r);
			s++;
			return r;
		}
	}

	spin_mutex& thread_logger::get_mutex(){
		static spin_mutex my_mutex;
		return my_mutex; 
	}
	
	thread_logger::~thread_logger(){
		if(1&&buffer.str().length()>0)
		{
			//lock for sequential IO
			//tbb::spin_mutex::scoped_lock lock=get_mutex();
			cout<<buffer.str();
			cout.flush();
		}
			//cout.flush();
		//buffer.clear();
	}
}}
