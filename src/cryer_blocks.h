#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <array>
#include <cassert>
#endif
#include "exceptions.h"

namespace fipster { namespace _cryer {

	using namespace std;

	typedef array<bool,2> bool2;
	typedef array<int,2> int2;


	struct forward  {const static size_t sp=0;const static int step=1;};
	struct backward {const static size_t sp=1;const static int step=-1;};

//block_t class ####################
	class block_t{
	public:
		//block_t::iterator class ####################
		template<class pass_tpl>
		class iterator{
			typedef pass_tpl pass_t;
			typedef iterator<pass_t> it_t;
			int it;
			friend block_t;

		public:
			//ctor
			iterator(int it):it(it)
			{};
			it_t& operator++(){
				it+=pass_t::step;
				return *this;
			}
			it_t& operator--(){
				it-=pass_t::step;
				return *this;
			}
			bool operator!=(const it_t& y)const{ return it!=y.it; }
			bool operator==(const it_t& y)const{ return !(*this!=y); }
			/**	\todo to slow?	*/
			it_t next()const{ return ++it_t(it); }
			/**	\todo to slow?	*/
			it_t previous()const{ return --it_t(it); }
			int operator*()const{ check(); return it; }

			void check()const{
				//This does not make sense, as the crossing block bounds does not invalidate an interator
				//#if(!FIPSTER_NDEBUG)
				//				assert(it*forward::step < forward::step*(*upper_bounds)[forward::sp] &&
				//					it*backward::step < backward::step*(*upper_bounds)[backward::sp]);
				//#endif
			}
		};
	private:
		bool2 done;
		bool2 is_on_margin_;
		int2 upper_bounds;

		void init(){upper_bounds[0]=-10;upper_bounds[1]=-10;}

	public :

		block_t(){init();};

		//ctor
		block_t(int begin_,int end_,int2& system_bounds)
		{
			init();
			begin_+=backward::step;

			//check if block lies within bounds
			FIPSTER_ASSERT(end_*forward::step <= forward::step*system_bounds[forward::sp] &&
					begin_*backward::step <= backward::step*system_bounds[backward::sp]);

			//set block bounds to begin and end and check,
			//if they lie in system_bounds
			set_end(iterator<forward>(end_));
			set_end(iterator<backward>(begin_));

			is_on_margin<forward>() = end_ == system_bounds[forward::sp];
			is_on_margin<backward>() = begin_ == system_bounds[backward::sp];
		}
		template<class pass_t>
		bool set_end(const iterator<pass_t>& it){

			if(it.it != upper_bounds[pass_t::sp]){
				//set upper bound
				upper_bounds[pass_t::sp]= it.it; //cannot use indirection operator, 
				//set other not done 
				mark_other_as_todo<pass_t>();
				return true;
			}
			return false;
		}

		/** \brief merges this block into it two blocks, if they are adjacent
			\param [in,out] it the iterator pointing to the
			following block, that is tried to be merged with this.
			\returns if a merge has taken place
		*/
		template<class block_it_t,class pass_t>
		void merge(block_it_t& it,pass_t pass){
			set_end(it->end(pass));
			is_on_margin(pass) = it->is_on_margin(pass);
			it.erase();
			is_done(pass)=false;
		}



		template<class pass_t>
		iterator<pass_t> begin(pass_t pass=pass_t()){
			return ++iterator<pass_t>(upper_bounds[1-pass_t::sp]);
		}
		template<class pass_t>
		iterator<pass_t> end(pass_t pass=pass_t()){
			return iterator<pass_t>(upper_bounds[pass_t::sp]);
		}
		template<class pass_t>
		bool& is_done(pass_t=pass_t()){
			return done[pass_t::sp];
		}
		template<class pass_t>
		void mark_other_as_todo(pass_t=pass_t()){
				done[1-pass_t::sp]=false;
		}
		template<class pass_t>
		bool& is_on_margin(pass_t=pass_t()){
			return is_on_margin_[pass_t::sp];
		}
	};
}}
