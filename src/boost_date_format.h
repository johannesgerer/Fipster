#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/date_time/posix_time/posix_time.hpp>
#include <string>
#endif



namespace fipster {

	namespace date {
		
		using namespace std;
		using namespace boost::posix_time;

		string format_date(const string& format,const ptime& date_time=second_clock::local_time());
	}
	using date::format_date;
}
