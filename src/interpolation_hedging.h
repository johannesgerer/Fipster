#pragma once
#include "option.h"

namespace fipster { namespace options {

const int dimension=0;

		typedef const boost::optional<grid_index<0> >&  res_t;
		
		struct interpolate_option{
			const option_field_t& o;
			interpolate_option(const option_field_t& o) : o(o) {}
			res_t at(grid_index<1> i) const{
				return o.additional.at(i).second;
			}
		};
		
			typedef discretized_space_field<1,1,array<boost::optional<grid_index<0> >,4 > >
				bands_field;
		
		struct interpolate_band{
			const bands_field& o;
			uint j;
			interpolate_band(const bands_field& o, uint j) : o(o),j(j) {}
			res_t at(grid_index<1> i) const{
				auto& a = o.at(i)[j];
				// FIPSTER_ASSERT(a);
				return a;
			}
		};

		template<class I, class I2, class O>
void interpolate_hedging(O source,
												 grid_ref g, grid_ref hg,
												 const I& gbegin, const I2& gend,
												 discretized_space_field<1>& result){

	result.resize(g);
	
	grid_iterator<1,true> i = gbegin;
	grid_index<1> write_i = gbegin.index();

	auto get_hi = [&]()
		{ return source.at(i.index()); };

	while(!get_hi()){
		result.at(write_i)=numeric_limits<double>::quiet_NaN();
		++i; ++write_i;
		if(i==gend) return;
	}
		
	double cur_x,cur_y;

	grid_index<0> cur_hi = *get_hi();

	auto get_x = [&]()
		{ return i.state()[dimension]; };

	auto update = [&](){
		cur_y = hg.coords[0][cur_hi.ind];
		cur_x = get_x();
	};

	update();

	double cur_interpol = cur_y;
	
	for(++i;i!=gend;++i){
		res_t b_temp_hi = get_hi();
		if(!b_temp_hi) break;
		
		auto temp_hi = *b_temp_hi;
		if(0){ //do nothing?
			result.at(++write_i) = hg.coords[0][temp_hi.ind];
			continue;
		}
		if(cur_hi != temp_hi){
			bool jump = abs(cur_hi.ind - temp_hi.ind) != 1;
			cur_hi = temp_hi;
			double last_x = cur_x;
			double last_y = cur_y;
			double last_interpol = cur_interpol;
			update();
			cur_interpol = jump ? last_y : 0.5*(cur_y + last_y);

			double slope = (cur_interpol - last_interpol) /  (cur_x - last_x);

			result.at(write_i)=last_interpol;
			for(++write_i;i!=write_i;++write_i)
				result.at(write_i)= last_interpol
					+  slope * (result.at(write_i) - last_x);
			if(jump)
				cur_interpol = cur_y;
		}else
			result.at(i.index()) = get_x();
	}

	for(;i.index()!=write_i;++write_i)
		result.at(write_i)= cur_interpol;
}

}}
