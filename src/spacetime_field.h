#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <memory>
#include <tuple>
#include <functional>
#include <comphash.h>
#include <boost/functional/hash.hpp>
#endif

#include "forward.h"
#include "grid_fields.h"
#include "memoized_node_factory.h"
#include "utility.h"

namespace fipster { namespace _fields {

		using namespace std;

		//arg
		template<class carg_t>
		struct field_arg_t : comphashable<field_arg_t<carg_t> > {
			btime time;
			grid_ptr grid;
			carg_t carg;

			field_arg_t(const field_arg_t<>& a, carg_t c)
				: time(a.time),grid(a.grid),carg(c) {}

			field_arg_t(btime time,const grid_ptr& grid, carg_t c):
				time(time),grid(grid),carg(c){};

			field_arg_t<> strip() const{
				field_arg_t<> a(time,grid,nothing());
				return a;
			}

			field_arg_t time_shifted(btime shift) const{
				auto a=(*this);
				a.time += shift;
				return a;
			}

			field_arg_t endow(carg_t&& carg) const{
				return field_arg_t(this->strip(),carg);
			}

			template<class T> FORCE_INLINE void my_combine(T& t)const{
				t(&field_arg_t::time,&field_arg_t::grid,&field_arg_t::carg);
			}
		};
		
             


		/**This pure class represents a scalar field in space time. With
			 two get_node methods for inner and boundary values creating
			 nodes for calculation of the fields' discretized values on a
			 certain grid, time-slice and possibly custom arg (carg_t).
		*/
		template<int sg         //specify the sub-grid to work on
						 ,class carg_t> // custom arg (default value declared in
		// configured_objects.h)
		struct spacetime_field_base
			: enable_shared_from_this<spacetime_field_base<sg,carg_t> >
		{
			typedef discretized_space_field<sg> v_t;
			typedef boundary_field_t<sg> bv_t;
			typedef typename mtypes<v_t>::sender_ptr sender_ptr;
			typedef typename mtypes<bv_t>::sender_ptr bv_sender_ptr;
			typedef field_arg_t<carg_t>              arg_t;
			
			typedef shared_ptr<spacetime_field_proxy<sg,carg_t> > proxy_t;
			
			virtual sender_ptr get_node(const arg_t& arg)=0;
			virtual bv_sender_ptr bv_get_node(const arg_t& arg)=0;

			proxy_t create_proxy(const carg_t& a) {
				auto this2 = this->shared_from_this();
				typedef shared_ptr<spacetime_field_base<sg,carg_t> > my_t;
				typedef shared_ptr<spacetime_field_proxy<sg,carg_t> > p_t;
				typedef pair<carg_t,my_t> k_t;
				static unordered_map<	k_t, p_t,boost::hash<k_t> > cache;
				auto p = make_pair(a,this2);
				auto c = cache.find(p);
				if(c!=cache.end()) return c->second;
				auto r = cache.insert(make_pair
															(p,
															 make_shared<spacetime_field_proxy<sg,carg_t> >
															 (a,this2)
															 ));
				return r.first->second;
			}
		};

		/** class with actuall node factories, and setup/delegation
				methods to be overloaded by derived classes.
		*/
		template<int sg, class carg_t> 
		struct spacetime_field : spacetime_field_base<sg,carg_t> {

			typedef spacetime_field<sg,carg_t> sf_t;
			typedef shared_ptr<sf_t> spacetime_field_ptr;
			

			using arg_t = typename spacetime_field::arg_t;

			/**************************************************** 
        Factory for inner values
			*/
			struct inner_f : memoized_node_factory<inner_f
																						 ,arg_t
																						 ,discretized_space_field<sg> >
			{
				// using node_factory_ptr  = typename spacetime_field::node_factory_ptr;
				using sender_ptr  = typename inner_f::sender_ptr;
				sf_t& parent;
		
				inner_f(sf_t& parent):
					parent(parent) {};

				sender_ptr setup_node(const arg_t& arg){
					return this->parent.inner_setup(arg);
				}
			};

			/**************************************************** 
        Factory for boundary values
			*/
			struct boundary_f : memoized_node_factory<boundary_f
																								,arg_t
																								,boundary_field_t<sg> >
			{
				using sender_ptr  = typename boundary_f::sender_ptr;

				sf_t& parent;
		
				boundary_f(sf_t& parent):
					parent(parent) {};

				sender_ptr setup_node(const arg_t& arg){
					return this->parent.boundary_setup(arg);
				}
			};

			// Pass down Types  ######################################
			typedef typename inner_f::sender_ptr          sender_ptr;
			typedef typename inner_f::result_sptr_t       result_sptr_t;
			typedef typename inner_f::result_t            result_t;
			typedef typename boundary_f::sender_ptr    bv_sender_ptr;
			typedef typename boundary_f::result_sptr_t bv_result_sptr_t;
			typedef typename boundary_f::result_t      bv_result_t;

			// Factories  ######################################
		private:
			shared_ptr<inner_f> inner;
			shared_ptr<boundary_f> boundary;
		public:

			sender_ptr get_node(const arg_t& arg) final 
			{return inner->get_node(arg);}
			bv_sender_ptr bv_get_node(const arg_t& arg) final
			{return boundary->get_node(arg);}
	
			//provide access to inner's cache (for debugging)
			const typename inner_f::cache_t& cache;

			// Pure Virtual Setup Functions ######################################
			virtual sender_ptr inner_setup(const arg_t& arg)=0;
			///the nodes can also return nullptr, indicating, that zero be
			///used as value (not implemented yet)
			virtual bv_sender_ptr	boundary_setup(const arg_t& arg)=0; 
	
	
			// GCC #pragma warning( disable: 4355 ) //	 'this' : used in base member initializer list
			/**ctor
				 \notes Warning 4355 is Ok, because this pointer is only stored for later use.
				 No virtual functions or data members are accessed in the constructor.
				 Later the virtual function setup_boundary_value_node will be called.
			*/
			// #pragma warning ( default: 4355 )
			//ctor
			spacetime_field()
				: inner(make_shared<inner_f>(*this))
				, boundary(make_shared<boundary_f>(*this))
				, cache(inner->cache)
			{};


		};
		
		/** proxy class, that provides a spacetime field with no carg and
				forwards requests to another spacetime field with cargs */
		template<int sg, class carg_t> 
		struct spacetime_field_proxy : spacetime_field_base<sg> {

			typedef shared_ptr<spacetime_field_base<sg,carg_t> > spacetime_field_ptr;
			typedef typename spacetime_field_proxy::sender_ptr sender_ptr;
			typedef typename spacetime_field_proxy::bv_sender_ptr bv_sender_ptr;
			typedef typename spacetime_field_proxy::arg_t arg_t;
			typedef field_arg_t<carg_t> arg2_t;
			
			carg_t carg;
			spacetime_field_ptr field;

			spacetime_field_proxy(carg_t carg,spacetime_field_ptr field)
				: carg(carg), field(field) {}
			
			sender_ptr get_node(const arg_t& arg) final override
			{return field->get_node(arg2_t(arg,carg));}
			bv_sender_ptr bv_get_node(const arg_t& arg) final override
			{return field->bv_get_node(arg2_t(arg,carg));}
			
		};
			
	}

	using _fields::spacetime_field;
	using _fields::spacetime_field_proxy;
	using _fields::spacetime_field_base;
	using _fields::field_arg_t;

}
