#pragma once
#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/property_tree/ptree.hpp>
#include <boost/optional.hpp>
#include <string>
#include <set>
#include <comphash.h>
#endif

#include "configuration_parser.h"

namespace fipster { 

	// struct neg_inf{};
	// typedef boost::variant<btime,neg_inf>; //extended btime
	
	namespace configuration {

		using namespace boost::property_tree;

		struct comb : stringifyable<comb> {
			btime offset, step_size;
			comb(const ptree &pt, string context);
			comb(btime offset, btime step_size);
		 	set<btime> times(btime start, btime stop) const;
			bool contains(btime time) const;
			bool contains(btime time,uint skip) const;
			btime next(btime time, btime stop) const;

			template<class T> FORCE_INLINE void my_combine(T& t)const{
				t(&comb::offset
					,&comb::step_size
					);
			}
		};
	
		struct stopped : comb {
			btime stop;
			stopped(const ptree &pt, string context);
			stopped(btime stop, btime offset, btime step_size);
		};
		
		struct combA : stopped {
			boost::optional<btime> start;
			combA(const ptree &pt, string context);
		};
		
		struct combB : stopped {
			using comb::times;
			btime start;
			combB(const ptree &pt, string context);
			combB(btime start, btime stop, btime offset, btime step_size);
		 	set<btime> times() const;
			// bool contains(btime time) const;
		};
	
	}
	using configuration::comb;
	using configuration::combA;
	using configuration::combB;
}
