#include "precompiled.h"
#include "boost_date_format.h"

namespace fipster{
	namespace date{

string format_date(const string& format,const ptime& date_time)
{
time_facet * facet = new time_facet(format.c_str());
std::ostringstream stream;
stream.imbue(std::locale(stream.getloc(), facet));
stream << date_time;
return stream.str();
}

	}}
