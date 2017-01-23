#include "precompiled.h"
#include "set_of_times.h"
#include "Math.h"

namespace fipster { 
	namespace configuration {

		comb::comb(const ptree &pt, string context)
			: offset    (getBtime(pt,"offset",context))
			, step_size (getBtime(pt,"stepSize",context)) {};

		comb::comb(btime offset, btime step_size)
			: offset(offset),step_size(step_size) {}

		bool comb::contains(btime time) const{
			return contains(time,1);
		}

		bool comb::contains(btime time,uint skip) const{
			return !((time-offset)%(step_size*skip));
		}

		stopped::stopped(const ptree &pt, string c)
			: comb(pt,c) 
			, stop(getBtime(pt,"stop",c)) {};

		combA::combA(const ptree &pt, string c)
			: stopped(pt,c)
			, start(pt.get_optional<btime>("start")) {};

		combB::combB(const ptree &pt, string c)
			: stopped(pt,c)
			, start(getBtime(pt,"start",c)) {};

		set<btime> combB::times() const{
			return times(start,stop);
		}

		// bool combB::contains(btime time) const{
		// 	return time <= stop && time >= start && comb::contains(time);


		stopped::stopped(btime stop, btime offset, btime step_size)
			: comb(offset,step_size),stop(stop) {}

		combB::combB(btime start, btime stop, btime offset, btime step_size)
				: stopped(stop,offset,step_size),start(start) {}

		
		/** \returns the smallest (largest) number "step_size*n + offset" with some integer n,
				that is larger (smaller) than "last", if stepsize is positive (negative).
				see tests
		*/
		template<class T>
		T next_offset_multiple(const T& last, const T& offset, const T& step_size){
			T d=(last-offset)%step_size;
			return last + ( (step_size>0 ? d<0 : d>0) ? 0 : step_size) - d;
		}
		
		btime comb::next(btime time,btime stop)const{
			return min(stop
								 ,next_offset_multiple(time,offset,step_size));
		}

		set<btime> comb::times(btime start, btime stop) const{
			set<btime> times;
			btime t=start;
			while(t<=stop)
				{
					times.insert(t);
					if(t==stop) break;
					t=next(t,stop);
				}
			return times;
		}
	}
}
