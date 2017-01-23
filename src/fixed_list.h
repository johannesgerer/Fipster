#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <vector>
#include <limits>
#include <new>
#include <utility>
#include <boost/throw_exception.hpp>
#endif

#include "utility.h"


namespace fipster { namespace _fixed_list {

	using namespace std;


/** \brief List Class that offers a fast insertion phase and then
			remove-only fixed phase. 
			
			See #fix method for motivation.
*/
template<class T,class my_size_tpl=unsigned int>
class fixed_list{

	typedef my_size_tpl my_size_t;
	typedef T value_type;
	typedef fixed_list<value_type,my_size_t> flist_t;

	struct list_element;//forward dec
	typedef vector<list_element> v_t;
	typedef typename v_t::iterator v_it_t;

	struct list_element{

		list_element()
#if(!FIPSTER_NDEBUG)
		:i(-1)
#endif
		{}
		value_type value;
		v_it_t next;
		v_it_t previous;
#if(!FIPSTER_NDEBUG)
		int i;
#endif
	};

	v_t v;
	my_size_t size_;
	bool fixed;

	//################  Iterators  #######################
	template<bool is_forward>
	struct iterator{
		friend flist_t;
		typedef iterator<is_forward> it_t;
		v_it_t it;
		//ctor
		iterator(v_it_t it):it(it){};
		bool operator!=(const it_t& y)const{
			return it!=y.it;
		}
		bool operator==(const it_t& y)const{ return !(*this!=y); }
		void erase(){ check2();
			it->previous->next = it->next;
			it->next->previous = it->previous;
#if(!FIPSTER_NDEBUG)
			//set the increment to itself (required for check())
			(is_forward ? it->next : it->previous)= it;
#endif
		}
		it_t& operator++(){
			it= is_forward ? it->next : it->previous;
			return *this;
		}
		T* operator->(){ 
			check2(); //was macht überhaupt hier und dann in ~_vector_val
			return &it->value;
		}
		T& operator*(){
			check2(); //und hier ??
			return it->value;
		}
		void check2(){
#if(!FIPSTER_NDEBUG)
			auto cp = *this;
			if(++cp == *this)
				BOOST_THROW_EXCEPTION(range_error("iterator not dereferencable: out of range"));
#endif
		}
	};

public:

	int size()const{
		return size_;
	}
	
	/** \brief  0th-layer storage access to fixed endings (in both directions) using underlying vector iterators*/
	iterator<1> end(){ require_fixed(); 
		return iterator<1>(v.begin()+size_+1); }
	iterator<0> rend(){ require_fixed(); 
		return iterator<0>(v.begin()); }

	/** \brief 1st-layer storage access to begin iterators (in both directions) using end(),rend() and links
	*/
	iterator<1> begin(){ return iterator<1>(rend().it->next); }
	iterator<0> rbegin(){ return iterator<0>(end().it->previous); }

	/** \brief 2nd-layer to front and back element using begin, rbegin */
	T& front(){ return *begin(); }
	T& back(){ return *rbegin(); }

	void require_fixed(){
#if(!FIPSTER_NDEBUG)
		if(!fixed)
			BOOST_THROW_EXCEPTION(
			runtime_error("container is not fixed. call fix() before calling this method"));
#endif
	}

	/** \brief reserves space for n list elements */
	void clear(){
		fixed=false;
		size_=0; 
	}

	/** \brief reserves space for n list elements 
	*/
	void resize(my_size_t n){
		v.resize(n+2); //init enough space (for size_ + end and rend)
#if(!FIPSTER_NDEBUG)
		//inscribe internal positions in list objects
		for(my_size_t i=0;i<n+2;i++) v[i].i=i;
#endif
	}

	/** \brief extends the list by one element
		\returns a reference to it
		*/
	template<class H>
	void emplace_back(H&& value){
		if(fixed)
			BOOST_THROW_EXCEPTION(
			runtime_error("container is fixed, user clear() before calling this method"));

		size_++;

		
		if(size_+2>v.size()) //(size_ + end and rend)
			v.resize(size_+2);

		//copy construct new element at the old position:
		//new (v[size_].value) T(move(value));
		//copy assign new element:
		v[size_].value = move(value);
	}

	/** \brief fixed the list (to prevent further extension) and builds links

		This is done for speed and to allow for reallocation
		within fill phase, which would invalidate iterators, that are
		used as safe next and previous pointers.		
	*/
	void fix(){
		if(!fixed){
			//build doubly linked list
			auto it = v.begin();
			auto myend = it+size_+1;
			it->previous = it;
			for(;it!=myend;it++){
				it->next = it+1;
				(it+1)->previous = it;
			}
			it->next = it;	
			
			fixed=true;
		}
	}
};


}
using _fixed_list::fixed_list;
}
