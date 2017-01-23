#pragma once


#include "cryer_blocks.h"
#include "fixed_list.h"
#include "utility.h"
#include "exceptions.h"

namespace fipster { namespace _cryer {
	
	using namespace std;

/**  Contains the routines and needed temporary storage to solve a LCP

	Solves LCP: M.z - q >= 0, z >= 0 , (M.z - q).z = 0

	Where M is tridiagonal and a negative M-Matrix!
	M_ii = a_i, M_(i,i-1) = b_i, M_(i,i+1) = c_i
	NOTE: c(n-1) and b(0) are not used 

	\tparam  M_t specifies the type of object representing the tridiagonal minkowsky matrix M, subject to the following concept:
	\code
	class M_t {
		double diag(int i) const;
		double subdiag(int i) const ;
		double superdiag(int i) const;
	}
	\endcode
	\tparam  q_t specifies the type of object representing the vector q, subject to the following concept:
	\code
	class q_t {
		double operator[](int i) const;
	}
	\endcode
	\tparam  z_t specifies the type of object representing the solution vector z, subject to the following concept:
	\code
	class z_t {
		double& operator[](int i);
	}
	\endcode
	
	\todo NOT IMPORTANT, DUE TO PRIMARY TESTS: write explicit tests for substitution functions and inplace_tridiagonal_solve, (as its hard to tell, if an error will have consequences at all in a complete Cryr run)
*/
template<typename M_t, typename q_t, typename z_t, bool calc_flop=false>
class cryer{

	//input and solution
	//vec_it_t q,z;
	M_t* M;
	q_t* q;
	z_t* z;

	double tol;

	static const int flop_per_div = 4, flop_per_multiplyadd = 1;

	vector<double> temp[2];
	int n,npass;
	int2 upper_bounds;
	fixed_list<block_t> blocks;
	int flop;
	const static bool check=false;

	//#########  public methods: ###############
public:
	cryer(q_t* q, z_t* z,double tol)
		:q(q),z(z),tol(tol),n(0){}

	int get_flop()const{ cerr<<"fkio"<<flop<<endl; return flop; };

	template<int validate_solution,bool validate_input, bool assume_q_positive>
		void compute(M_t* M, int n_);

private:
	
	template<class T>
	void set_z(const T& it,double v){
		if(check && (v<0 || v>1e70)){
			cout<<v<<endl;
		}
		z->operator[](*it) = v;
	}

	void fill_test(M_t* M);

	template<bool do_throw>
	void validate_single_line(int i, double tol);

	template<class block_it_t,class pass_t>
	bool perform_pass(block_it_t block,block_it_t block_end,pass_t pass);

	template<bool assume_q_positive>
	void reinit(int n_);

	void validate_sol(double tol);
	void validate_in();

	void check_initial_block();
	void build_initial_block_structure();

	//################# memory/storage access #################
	template<class pass_t>
	double sub_diag(const block_t::iterator<pass_t>& it){
		return pass_t::step == 1 ? M->subdiag(*it) : M->superdiag(*it);
	}
	template<class pass_t>
	double super_diag(const block_t::iterator<pass_t>& it){
		return pass_t::step == -1 ? M->subdiag(*it) : M->superdiag(*it);
	}
	template<class pass_t>
	double diag(const block_t::iterator<pass_t>& it){
		return M->diag(*it);
	}
	template<class pass_t>
	double z_temp(const block_t::iterator<pass_t>& it){
		return temp[pass_t::sp][2**it];
	}
	template<class pass_t>
	void set_z_temp(const block_t::iterator<pass_t>& it,double v){
		if(check && (v<0 || v>1e70)){
			cout<<"set_z_temp: "<<*it<<" : "<<v<<endl;
		}
		temp[pass_t::sp][2**it]=v;
	}

	template<class pass_t>
	double& a_temp(const block_t::iterator<pass_t>& it){
		return temp[pass_t::sp][1+2**it];
	}

	//#########  generic tridiagonal solvers: ###############
	template<class block_it>
	FORCE_INLINE void outofplace_forward_subsitution_step(const block_it& it);
	template<class block_it>
	FORCE_INLINE void outofplace_forward_subsitution(block_it& it,const block_it& end);
	template<class pass_t>
	FORCE_INLINE void outofplace_backward_subsitution(block_t& block);
	template<class it_t>
	FORCE_INLINE void inplace_tridiagonal_solve(it_t it, it_t end);


	//#########  utility: ###############
	/*template<class it_t,class block_it>
	FORCE_INLINE void finish_block_extension(const it_t& it,bool& changed,block_it& block);*/
	template<class pass_t,class block_it>
	FORCE_INLINE block_t::iterator<pass_t> init_block(block_it& block,pass_t);
};

//#################################################################
//#################################################################
//Definitions:
//#################################################################
//#################################################################

/** \brief performs one forward substitution step using two temporary arrays without altering the solution
	\tparam block_it the pass_t of the passed iterator specifies in which direction the forward pass should be performed
	\param it marks the location at which the step should be performed
*/
template<typename M_t, typename q_t, typename z_t, bool calc_flop>
template<class block_it>
FORCE_INLINE void cryer<M_t,q_t,z_t,calc_flop>
	::outofplace_forward_subsitution_step(const block_it& it)
{
	if(calc_flop) flop+=flop_per_div+2*flop_per_multiplyadd;
	double r = sub_diag(it.next()) / a_temp(it);
	a_temp(it.next()) = diag(it.next()) - r * super_diag(it);
	if(check) cout<<a_temp(it.next())<<endl;
	set_z_temp(it.next(),q->operator[](*it.next()) - r * z_temp(it));
}
	
/** \brief performs a forward substitution using two temporary arrays without altering the solution
	\tparam block_it the pass_t of the passed iterator specifies, in which direction the forward pass should be performed
	\param it marks the location at which the forward substitution should start
	\param end marks the location at which the forward substitution should stop
*/
template<typename M_t, typename q_t, typename z_t, bool calc_flop>
template<class block_it>
FORCE_INLINE void cryer<M_t,q_t,z_t,calc_flop>
	::outofplace_forward_subsitution(block_it& it,const block_it& end)
{
	for (;it.next()!=end;++it)
		outofplace_forward_subsitution_step(it);
}

/** \brief performs the corresponding backward substitution for an
		earlier outofplace_forward_subsitution().
			
	To do a backward substitution, a block::iterator<pass_t> is decremented in a loop!

	\tparam pass_t has to match the pass_t of the earlier outofplace_forward_subsitution().
*/
template<typename M_t, typename q_t, typename z_t, bool calc_flop>
template<class pass_t>
FORCE_INLINE void cryer<M_t,q_t,z_t,calc_flop>
	::outofplace_backward_subsitution(block_t& block)
{
	auto it=block.end(pass_t()).previous();
	auto begin = block.begin(pass_t()).previous();
	if(calc_flop) flop+=(flop_per_div+flop_per_multiplyadd)*(*it-*begin-1);
	for (--it;it!=begin;--it){
		set_z(it,(z_temp(it) - super_diag(it)*z->operator[](*it.next()))/a_temp(it));
	}
	///\todo the last is already valid, if the block is on no margin. (check that and use it)
}

/** \brief Solves the tridiagonal system in one block

		Code taken from "numerical recipes". Current pass direction does not matter
			
	*/

template<typename M_t, typename q_t, typename z_t, bool calc_flop>
template<class it_t>
FORCE_INLINE void cryer<M_t,q_t,z_t,calc_flop>
	::inplace_tridiagonal_solve(it_t begin, it_t end)
{
	auto it=begin;

	if(calc_flop){
		int n=*end-*it;
		flop+=flop_per_div*((n-1)*2+1)+flop_per_multiplyadd*3*(n-1);
	}
	//Solves for a vector z[0..n-1] the tridiagonal linear set given by equation (2.4.1). sub_diag(0..n-1),
	//b[0..n-1], super_diag(0..n-1), and r[0..n-1] are input vectors and are not modified.
	double bet;
	FIPSTER_ASSERT(diag(it) != 0.0);
	//If this happens, then you should rewrite your equations as a set of order N  1, with u1
	//trivially eliminated.
	set_z(it,q->operator[](*it)/(bet=diag(it)));
	for (++it;it!=end;++it) { //Decomposition and forward substitution.
		temp[0][*it]=super_diag(it.previous())/bet;
		bet=diag(it)-sub_diag(it)*temp[0][*it];
		//2.4 Tridiagonal and Band-Diagonal Systems of Equations 57
		FIPSTER_ASSERT(bet != 0.0);// Algorithm fails; see below.
		set_z(it,(q->operator[](*it)-sub_diag(it)*z->operator[](*it.previous()))/bet);
	}
	for (--it;it!=begin;--it){
		set_z(it.previous(),
			z->operator[](*it.previous()) - temp[0][*it]*z->operator[](*it)); //Back substitution.
	}
}
	
//template<class it_t,class block_it>
//FORCE_INLINE void cryer<M_t,q_t,z_t,calc_flop>
//	::finish_block_extension(const it_t& it,bool& changed,block_it& block)
//{
//	block->set_end(it);
//}

/** sets a_temp and z_temp and returns an iterator to the 
	beginning of the block (in accordance with the path)
*/
template<typename M_t, typename q_t, typename z_t, bool calc_flop>
template<class pass_t,class block_it>
/*FORCE_INLINE*/ block_t::iterator<pass_t> cryer<M_t,q_t,z_t,calc_flop>
	::init_block(block_it& block,pass_t)
{
	auto it=block->begin(pass_t());
	a_temp(it) = diag(it);
 	set_z_temp(it,q->operator[](*it));
	return it;
} //NRVO?

/** \brief calculate the LCP's solution and optionally validate the solution and block structure for given problem vectors
	\tparam validate_solution 1 - Validate LCP (ignoring the given block structure); 2 - DEBUG block structure and solution
	\param tol is a tolerance for rounding errors occurring in the numerical solution of a tridiagonal system within a "block". Outside of "blocks", the solution is simply set zero, thus no inaccuracies induced here.
	\tparam assume_q_positive if this is set true, effectively this equation is solved: M.u = q
*/
template<typename M_t, typename q_t, typename z_t, bool calc_flop>
template<int validate_solution,bool validate_input, bool assume_q_positive>
void cryer<M_t,q_t,z_t,calc_flop>
	::compute(M_t* M_,int n_)
{
	M=M_;
	reinit<assume_q_positive>(n_);

	if(validate_input)
		validate_in();

	if(assume_q_positive){
		inplace_tridiagonal_solve(block_t::iterator<forward>(0),
									block_t::iterator<forward>(n));
	}else{
		build_initial_block_structure();
	
		//set the solution to zero (this is important)
		for (int j = 0; j < n; ++j) z->operator[](j) = 0.;

		//check, if the LCP does not posses the trivial solution (z=0):
		if(blocks.size()>0){

			blocks.fix();

			check_initial_block();

			/* 4. Iterate, performing forward and backward passes: */
			//choose a direction to start with based on which direction's first block
			//is closer to the beginning:
			bool forwardpass = *blocks.front().begin<forward>() <= 
									n - *(--blocks.back().end<forward>());
			bool changed; npass=0;
			do{
				changed = forwardpass ? perform_pass(blocks.begin(),blocks.end(),forward())
					: perform_pass(blocks.rbegin(),blocks.rend(),backward());
				forwardpass=!forwardpass;
				npass++;
			}while(changed || npass<2); /* if(nothingChanged .and. nPass>=2) exit passes */

			//iterator through all final blocks and compute the final solution
			auto end = blocks.end();
			for(auto block = blocks.begin();block != end;++block)
			{
				//PERFORMANCE: there are cases, where  part of a 
				// out-of-place substitution has already happened for part of block,
				// it as not been completed, because the block has been extended (this also includes joins)
				// to the boundary. In such a case, the substitution is canceled,
				// because it's not needed in this direction,
				// but it is likely to not be the final pass, thus to be changed in the future
				// HOWEVER: The position could marked instead of a simple "done"

				//if forward substitution in forward direction has been performed, 
				//do the corresponding backward substitution
				if(block->is_done(forward())) 
					outofplace_backward_subsitution<forward>(*block);
				//if there was no forward substitution in the forward direction, check the backward direction
				else if(block->is_done(backward())) 
					outofplace_backward_subsitution<backward>(*block);
				//if there was no forward substitution performed at all, solve the block from scratch
				else
					inplace_tridiagonal_solve(block->begin(forward()),block->end(forward()));
			}
		}
	}

	switch(validate_solution){
	case 1:
		//Validate LCP (ignoring the given block structure)
		for(int i=0;i<n;i++) 
			validate_single_line<true>(i,tol);
		break;
	case 2:
		//DEBUG block structure and solution;
		validate_sol(tol);
		break;
	}

} 

/**
\returns a boolean value indicating whether this pass changed the block structure or not
*/
template<typename M_t, typename q_t, typename z_t, bool calc_flop>
template<class block_it_t,class pass_t>
bool cryer<M_t,q_t,z_t,calc_flop>
	::perform_pass(block_it_t block,block_it_t block_end,pass_t pass){

	bool changed=false;

 	auto it = init_block(block,pass);
	//There is always at least one block, because the case of no blocks
	//is excluded earlier!! (so no check is needed at the beginning) WHERE?

	while(true){

		//if the block lies on the margin, it cannot grow, 
		//thus we don't even need z and we are done 
		if (block->is_on_margin(pass)) 
			return changed;

		if (block->is_done(pass)){
			/* if z(end) for this block has been solved yet */
			//process next block
			if(++block==block_end) 
				return changed;
			it = init_block(block,pass);
			continue;
		}


		/* perform forward substitution */
		outofplace_forward_subsitution(it,block->end(pass));
 		set_z(it,z_temp(it) / a_temp(it));

		if(calc_flop) flop+=flop_per_div;

		// extend block until something happens
		while(true)
		{
			double d = z->operator[](*it) * sub_diag(it.next()) - q->operator[](*it.next());
			if(calc_flop) flop+=flop_per_multiplyadd;
			++it;
			//if we are only one step from the global boundary 
			if (*it.next()==upper_bounds[pass_t::sp]){
				//extended block if d<0
				if (d < 0.){ //\todo: what about floating points special cases NAN ...
					++it; 
					block->is_on_margin(pass)=true;
				}else
					block->is_done(pass) = true;

				//register new end and eventual changes from block extensions
				changed=block->set_end(it) || changed;
				//quit this pass (as we reached the global bound)
				return changed;
			}else{ 
				/** \todo  z is only (structural) non-zero, if the it.next
						is the start of a new block. only then a merge could happen.*/

				//make more complicated check, too see if block is finished
				if(calc_flop) flop+=flop_per_multiplyadd;
				if(d >= -z->operator[](*it.next())*super_diag(it)){ //\todo as the product is always positive, d must be positive. check that?
					block->is_done(pass) = true;

					//register new end and eventual changes from block extensions
					changed=block->set_end(it) || changed;
					//and process next block:
					if(++block==block_end) 
						return changed;
					it = init_block(block,pass);
					break;
				}else{//if the block has to be extended
					//try to merge it with
					//the next block (if they touch)
					auto next_block = block; ++next_block;
					if(next_block!=block_end && next_block->begin(pass)==it.next()){
						block->merge(next_block,pass);
						//process merged block in the main loop (e.g. outofplace_forward_subsitution)
						--it;
						changed=true;
						break;

					//if the block was not merged
					}else{//update the forward substitution and z
						/* p = p + 1 \todo: cryer wants to count block extensions. is this needed? */
						outofplace_forward_subsitution_step(it.previous());
						if(calc_flop) flop+=flop_per_div;
						set_z(it,z_temp(it) / a_temp(it));
					}//if there is next block touching && merged
				}//if not extend
			}//if on bound
		}
	}// while blocks

	return changed;
}

template<typename M_t, typename q_t, typename z_t, bool calc_flop>
template<bool assume_q_positive>
void cryer<M_t,q_t,z_t,calc_flop>
	::reinit(int n_){
	if(n_>n){
		blocks.resize((n_+1)/2);//there are at most (n+1)/2 blocks			
		temp[1].resize(2*n_);
		temp[0].resize(2*n_);
	}

	n=n_;
	if(calc_flop) 
		flop=0;
	if(!assume_q_positive){
		upper_bounds[forward::sp]=forward::step > 0 ? n : -1;
		upper_bounds[backward::sp]=backward::step > 0 ? n : -1;
	}
}
/** \brief validates a single line of the LCP for a given values
	\param i is the index of the LCP line in question
	\param tol is a tolerance for rounding errors occurring in the numerical solution of a tridiagonal system within a "block". Outside of "blocks", z is simply set zero, thus no inaccuracies induced here.
*/
template<typename M_t, typename q_t, typename z_t, bool calc_flop>
template<bool do_throw>
void cryer<M_t,q_t,z_t,calc_flop>
	::validate_single_line(int i, double tol){
		//cout<<"validate_single_line: "<<i<<endl;
	bool print = false;
	double t = M->diag(i)*z->operator[](i)-q->operator[](i);	
	
	if(i>0)
		t += M->subdiag(i)*z->operator[](i-1);	
	if(i<n-1)
		t += M->superdiag(i)*z->operator[](i+1);

	if(z->operator[](i)<0 || t<-tol){
		print = true;
		cout<<"inequality violation: t="<<t<<" z="<<z<<" z(i)="<<z->operator[](i)<<endl;
		if(do_throw)
			FIPSTER_THROW_EXCEPTION(runtime_error("LCP inequality violation in Cryer"));
	}else if(z->operator[](i)!=0.0 && t>tol){
		print = true;
		cout<<"complementarity violation: t="<<t<<" z="<<z->operator[](i)<<endl;
		if(do_throw)
			FIPSTER_THROW_EXCEPTION(runtime_error("LCP complementarity violation in Cryer"));
	}
	if(print)
		cout<<i<<": "<<z<<"  "<<t<<endl;
}

/** \brief validates a the solution and its block structure
	
	The block structure might be considered of only minor importance, as a correct 
	solution stays correct even, if the block structure used in its generation is flawed.

	\param tol is a tolerance for rounding errors occurring in the numerical solution of a tridiagonal system within a "block". Outside of "blocks", the solution is simply set zero, thus no inaccuracies induced here.
*/
template<typename M_t, typename q_t, typename z_t, bool calc_flop>
void cryer<M_t,q_t,z_t,calc_flop>
	::validate_sol(double tol){
	
	// validate block structure and 
	//Reminder: A block marks a region, where the trivial constraint (z=0) is not active

	int consistency_counter=0;
	int i=0;

	if(blocks.size()>0){
		auto block=blocks.begin();

		for(;block!=blocks.end();++block)
		{
			auto next_block_begin=block->begin<forward>();
			//check if there is a separation between the two blocks
			FIPSTER_ASSERT(i==0 || i<*next_block_begin);
			for(;i<*next_block_begin;++i){
				consistency_counter++;
				validate_single_line<false>(i,tol);
			}

			//cout<<"START"<<endl;

			auto next_block_end=block->end<forward>();
			for(;i<*next_block_end;++i){
				consistency_counter++;
				validate_single_line<false>(i,tol);
			}
			//cout<<"END"<<endl;
		}
	}
	for(;i<n;++i){
		consistency_counter++;
		validate_single_line<false>(i,tol);
	}
	FIPSTER_ASSERT(consistency_counter==n);
}

/** 2. determine set Q = { i | q_i > 0 } and use its 
   sign changes as the initial block structure: */
template<typename M_t, typename q_t, typename z_t, bool calc_flop>
void cryer<M_t,q_t,z_t,calc_flop>
	::build_initial_block_structure(){

	blocks.clear();

	/* search for sign changes */
	int block_begin=-1;
	for (int j = 0; j < n; ++j) {
		//If j-1 was in a block */
		if (block_begin>=0) {
			//if the sign changed,terminate block */
			if (q->operator[](j) <= 0.){
				blocks.emplace_back(block_t(block_begin,j,upper_bounds));
				block_begin = -1;
			}
		} else if (q->operator[](j) > 0.) //j-1 was not in a block and if a new block starts:
			block_begin=j;
	}
	//if j=n was in a block */
	if (block_begin>=0)
		blocks.emplace_back(block_t(block_begin,n,upper_bounds));

	//see if cryer algorithm is needed (and brennan schwartz(?) would not apply)
	//FIPSTER_ASSERT(blocks.size()<2);
}
template<typename M_t, typename q_t, typename z_t, bool calc_flop>
void cryer<M_t,q_t,z_t,calc_flop>
	::check_initial_block(){
	
	// validate block structure and 
	//Reminder: A block marks a region, where q>0
	int consistency_counter=0;
	int i=0;

	if(blocks.size()>0){
		auto block=blocks.begin();

		for(;block!=blocks.end();++block)
		{
			auto next_block_begin=block->begin<forward>();
			//check if there is a separation between the two blocks
			FIPSTER_ASSERT(i==0 || i<*next_block_begin);
			for(;i<*next_block_begin;++i){
				consistency_counter++;
				if(q->operator[](i) > 0)
					cout<<"q>0 outside of block: "<<q->operator[](i)<<endl;
			}

			//cout<<"START"<<endl;

			auto next_block_end=block->end<forward>();
			for(;i<*next_block_end;++i){
				consistency_counter++;
				if(q->operator[](i) <= 0)
					cout<<"q<=0 inside block: "<<q->operator[](i)<<endl;
			}
			//cout<<"END"<<endl;
		}
	}
	for(;i<n;++i){
		consistency_counter++;
		if(q->operator[](i) > 0)
			cout<<"q>0 outside of block: "<<q->operator[](i)<<endl;
	}
	FIPSTER_ASSERT(consistency_counter==n);
}

/** \brief check for negative M-Matrix property */
template<typename M_t, typename q_t, typename z_t, bool calc_flop>
void cryer<M_t,q_t,z_t,calc_flop>
	::validate_in(){

	for(int i=0;i<n;i++){
		double	sub		= i==0	? 0 : M->subdiag(i),
				super	= i==n-1? 0 : M->superdiag(i);

		if(sub>0 || super>0) FIPSTER_THROW_EXCEPTION(
			runtime_error("wrong input matrix: positive offdiagonals"));
		else if(M->diag(i)<-super-sub) FIPSTER_THROW_EXCEPTION(
			runtime_error("wrong input matrix: not diagonally dominant"));
		else if(M->diag(i)==-super-sub) FIPSTER_THROW_EXCEPTION(
			runtime_error("wrong input matrix: not strictly diagonally dominant"));
	}
}

/** \brief fills the Matrix with a slightly diagonally dominant Laplace matrix for (testing purposes)
*/
template<typename M_t, typename q_t, typename z_t, bool calc_flop>
void cryer<M_t,q_t,z_t,calc_flop>
	::fill_test(M_t* M){
	
	for(int i=0;i<n;i++)
	{
		M->diag_ref(i)=2.01;
		M->subdiag_ref(i)=-1;
		M->superdiag_ref(i)=-1;
	}
}

struct M_test{
	vector<double> sub,super,d;
	M_test(int n){
		sub.resize(n);
		super.resize(n);
		d.resize(n);
	}
	double& diag(int i) { return d[i];}
	double& subdiag(int i) { return sub[i];}
	double& superdiag(int i) { return super[i];}

};



}
using _cryer::cryer;

// #include <cstdlib>

// template<bool a>
// void cryer_tridiagonal_subroutine_test(int n){
// 	cryer<_cryer::M_test,vector<double>,vector<double>,true> C;
// 	_cryer::M_test M(n);
// 	vector<double> q(n),z(n);
// 	for(int i=0;i<n;i++) q[i]=rand()-RAND_MAX/2;
// 	q[0]=9470;
// 	q[1]=-15468;
// 	q[2]=-1249;
// 	q[3]=-13464;
// 	q[4]=8959;
// 	q[5]=-1200;
// 	q[6]=-5172;
// 	q[7]=10749;
// 	q[8]=3243;
// 	q[9]=-6063;
// 	C.compute<1,true,true>(&M,&q,&z,n,1e-10);
// }
}
