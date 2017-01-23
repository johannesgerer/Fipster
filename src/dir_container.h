#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <type_traits>
#include <fstream>
#endif

#include "integer_types.h"

namespace fipster { namespace _dir {

	using namespace std;
	
	class dir_iterator;

	
	/** this class can be used to safely iterate over 
		all neighbors in dir_containers. 
		The (internal) access convention correspond to dir_iterator's
		internal conventions.
	*/ 
	class neighbor_iterator{
		
		uint nNeighs;
		bool _valid;
		uint neig;

		friend class dir_iterator; 

		neighbor_iterator(uint neig,uint nDirs)
			:nNeighs(2*nDirs+1)
			,_valid(neig<2*nDirs+1)
			,neig(neig)
		{};
		
	public:
		
		neighbor_iterator& operator=(const neighbor_iterator& that){
			neig=that.neig;
			nNeighs=that.nNeighs;
			_valid=that._valid;
			return *this;
		}

		neighbor_iterator(uint nDirs)
			:nNeighs(2*nDirs+1)
			,_valid(true),neig(0)
		{};

		uint operator()() const{
			return neig;///\todo check for validity?
		};

		bool valid()const{ return _valid;}
		neighbor_iterator& operator++(){
			neig++;
			_valid=neig<nNeighs;
			return *this;
		}
	};

	/** this class can be used to safely iterate over 
		all directions in dir_containers. 
		The (internal) access convention is only used explicitly once
		in stencil_points_type::construct(...)
		Every call to grid_field_dir_container::UNSAFE_data_access()
		has to line up with the internal convention.
	*/ 
	class dir_iterator{
		uint _left,_right,_center;
		bool _valid;
	public:

		const uint nDirs;

		dir_iterator(uint nDirs)
			:nDirs(nDirs)
		{ 
			if(nDirs==0){
				_valid=false;
				return;
			}
			_valid=true;
			_left=0;
			_center=nDirs*2;
			_right=_center-1;
		}
		///\todo check for validity?
		neighbor_iterator at_center_it()const { return neighbor_iterator(at_center(),nDirs); }
		neighbor_iterator at_it()const { return neighbor_iterator(at(),nDirs); }
		neighbor_iterator at_opposite_it()const { return neighbor_iterator(at_opposite(),nDirs); }
		uint at_center()const{ return _center; }///\todo make private
		uint at_opposite()const{ return _right; }///\todo make private
		uint at()const{ return _left; }///\todo make private
		bool valid()const{ return _valid; }
		dir_iterator& operator++(){ 
			_left++;
			_right--;
			_valid = _left<nDirs;
			return *this;
		}
		/** \note this uses unsigned integers, so there are no negative values,
			they overrun instead.
			*/
		dir_iterator& operator--(){  
			_left--;
			_right++;
			_valid = _left<nDirs; 
			return *this;
		}

		bool is_last() const{
			return _left + 1 == nDirs;
		}

		friend ostream& operator<<(ostream& out,const dir_iterator& dit){
			out<<"left: "<<dit._left<<" right: "<<dit._right<<" center: "<<dit._center<<" nDirs: "<<dit.nDirs;
			return out;
		}
	};

	

	enum {simple,center,opposite};
	
	class stencil_points_t;

	//############################################################
	/** \brief This class wraps a vector<T>, containing one object for each direction
		in a stencil

		\tparam flag specifies, if each direction contains only one object (simple), or two objects (opposite), of if there are two for each direction and one for the "center" (center)
	*/
	template<class T,int flag=simple>
	class dir_container{
		
	protected:

		uint detemine_size(){ return detemine_size(_nDirs); }
		uint detemine_size(uint nDirs){			
			switch(flag){
			case simple: return nDirs;
			case opposite: return nDirs*2;
			case center: return nDirs*2+1;
			}
			return -1;
		};
		uint detemine_nDirs(){
			uint size = vec.size();
			switch(flag){
			case simple: return size;
			case opposite: return size/2;
			case center: return (size-1)/2;
			}
			return -1;
		};


		typedef T value_t;
		typedef vector<T> vec_t;
		typedef dir_container<T,flag> dir_cont_t;

		uint _nDirs;
		
		vec_t vec;

	public:
		//ctor zero
		dir_container():_nDirs(0){}

		//ctor
		dir_container(uint nDirs)
		{
			resize(nDirs);
		}

		void resize(uint nDirs)
		{
			_nDirs=nDirs;
			vec.resize(detemine_size());
		}

		//reset using existing vector
		void reset(vec_t&& vec_)
		{
			vec = move(vec_);
			_nDirs = detemine_nDirs();
		}

		uint nDirs()const{return _nDirs;}

		//const
		const T& operator[](const dir_iterator& it)const { return at(it); }
		const T& at(const dir_iterator& it)const
		{ static_assert(flag == simple||flag == opposite||flag == center,"not available");
			return vec[it.at()]; }
		const T& at(const neighbor_iterator& it)const
		{ static_assert(flag == center,"not available");
			return vec[it()]; }
		const T& at_opposite(const dir_iterator& it)const
		{ static_assert(flag == opposite||flag == center,"not available");
			return vec[it.at_opposite()]; }
		const T& at_center(const dir_iterator& it)const
		{ static_assert(flag == center,"not available");
			return vec[it.at_center()]; }

		//non-const
		T& operator[](const dir_iterator& it){ return at(it); }
		T& at(const dir_iterator& it)
		{ static_assert(flag == simple||flag == opposite||flag == center,"not available");
		return vec[it.at()]; }
		T& at(const neighbor_iterator& it)
		{ static_assert(flag == center,"not available");
			return vec[it()]; }
		T& at_opposite(const dir_iterator& it)
		{ static_assert(flag == opposite||flag == center,"not available");
		return vec[it.at_opposite()]; }
		T& at_center(const dir_iterator& it)
		{ static_assert(flag == center,"not available");
		return vec[it.at_center()]; }

		/** Usage:
		\code
			for(auto i = dir_it();i.valid();++i);
		\endcode
		\return dir_iterator pointing at first dir
		*/
		dir_iterator dir_it() const{
			return dir_iterator(_nDirs);
		}

		neighbor_iterator neigh_it() const{
			return neighbor_iterator(_nDirs);
		}
	};



}

using _dir::dir_iterator;
using _dir::neighbor_iterator;
using _dir::dir_container;


}
