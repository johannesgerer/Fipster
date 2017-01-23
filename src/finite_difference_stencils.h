#pragma once
#include "finite_difference_weights_config.h"

#include "utility.h"
#include "Math.h"
#include "grid.h"

#include "Eigen.h"
#include "exceptions.h"
#include "global_config.h"

#include "dir_container.h"

namespace fipster { namespace finite_difference_weights {

	using namespace Eigen;

	typedef position_t stencil_point_t;
	typedef vector<stencil_point_t> stencil_point_vec_t;

	
	// #################################################
	// ##############   PHASE 3:			############
	// ##############   Node Bodies         ############
	// #################################################

	

	template<class T>
	inline void write_diagonal_stencils(int sign,T& it, int k, int j){
		(*++it)[k] = 1;
		(*  it)[j] = sign;
		(*++it)[k] = -1;
		(*  it)[j] = -sign;
	}

	/** \brief append Stencil Points to container points (without the center point!)
		\todo For StencilOneDiagonal: choose sign of write_diagonal_stencils according to problem
	*/
	void push_back_stencil_points(int D, stenciltype_t type, stencil_point_vec_t& points)
	{
		switch(type){

		case StencilFull: 	
			{
				uint nPoints = pow2(3,D)-1;
				stencil_point_t p(D);
				p.setZero();

				while(points.size() < nPoints) {
					for (int j = 0; j < D; ++j) {
						if(p[j]==1){	p[j]=-1;
						}else{ 			p[j]++;	 break;					
						}
					}
					points.push_back(p);
				}
			}
			break;

		case StencilNoDiagonals:
			{
				points.resize(points.size()+2*D,stencil_point_t::Zero(D));
				auto it = points.begin();
				for(int i=0;i<D;i++,++it){
					(*it)[i]=1;
					(*++it)[i]=-1;
				}
				break;
			}
			
		case StencilTwoDiagonals:
		case StencilOneDiagonal:
			{	
				int more=(D*D-D)*(type==StencilTwoDiagonals?2:1);
				points.reserve(points.size()+2*D+more);

				//first get points without diagonals
				push_back_stencil_points(D,StencilNoDiagonals,points);

				auto it = --points.end();

				points.resize(points.size()+more,stencil_point_t::Zero(D));
			
				//iterate through every 2D plain in the cube
				for (int k = 0; k < D; ++k) {
					for (int j = k + 1; j < D ; ++j) {
						int sign=1;
						write_diagonal_stencils(sign,it,k,j);
						if(type==StencilTwoDiagonals)
							write_diagonal_stencils(-sign,it,k,j);
					}
				}
			}
			break;
			
		}

	}

	//utility function
	template<class T>
	bool is_opposite(T& A, T& B)
	{
		int l=max(A.size(),B.size());
		for(int i=0;i<l;i++)
			if(A[i] != -B[i])
				return false;

		return true;
	}

	
	bool pointless(const stencil_point_t& a,const stencil_point_t& b){
		for(int i=a.size()-1;i>0;i--)
			if(a(i)!=b(i)) return a(i)<b(i);
		return a(0)<b(0);
	}
	bool pointequal(const stencil_point_t& a,const stencil_point_t& b){
		return (a.array()==b.array()).all();
	}

	struct stencil_points_t : dir_container<stencil_point_t,_dir::center>{
		/**
		\returns the number of points = (2*nDirs)
		*/
		uint proper_neighs()const{return 2*nDirs();}
		
		/** construct stencil points including center point */
		void construct(int D, stenciltype_t type){

			push_back_stencil_points(D,type,vec);
			
			uint nPoints=0;

			switch(type){	
			case StencilFull: 			
				nPoints = pow2(3,D)-1;	
				break;
			case StencilNoDiagonals:	
				nPoints = 2*D;
				break;
			case StencilTwoDiagonals:
				nPoints = 2*D*D;
				break;
			case  StencilOneDiagonal:  
				nPoints = D+D*D;
				break;
			}

			if(nPoints!=vec.size())
				FIPSTER_THROW_EXCEPTION(runtime_error("Wrong number of stencil points"));
	

			//Sort points to get order corresponding to directions
			sort(begin(vec),end(vec),pointless);			

			auto new_end=unique(begin(vec),end(vec),pointequal);		

			if(nPoints!=new_end-begin(vec))
				FIPSTER_THROW_EXCEPTION(runtime_error("Duplicate stencil points"));

			if(vec.size()%2!=0)
				FIPSTER_THROW_EXCEPTION(
				runtime_error("even number of stencil points required for splitting. nPoints: "+toS(vec.size())));

			
			_nDirs=vec.size()/2;

			auto i = dir_it();

			//Check, if the dir_iterators point to distinct points
			set<decltype(&at(i))> dirits;
			
			for(;i.valid();++i){
				if(!is_opposite(at(i),at_opposite(i)))
					FIPSTER_THROW_EXCEPTION(
					runtime_error("Stencil contains incomplete splitting direction."));
				dirits.insert(&at(i));
				dirits.insert(&at_opposite(i));
			}
			
			//Add and check origin
			vec.push_back(position_t::Zero(D));
			dirits.insert(&at_center(i));
			FIPSTER_ASSERT(at_center(i)==position_t::Zero(D));

			//Check, if the dir_iterators point to distinct points
			FIPSTER_ASSERT(dirits.size()==nPoints+1);
			
		}
	};

}}

