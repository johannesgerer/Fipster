//OLD UNUSED

#include "precompiled.h"

#include "financial_derivative.h"
#include "Math.h"

TLS<time_stepping_work_storage> financial_derivative::storage;

////############################  get_result  ######################################
///
//Task* financial_derivative::get_result( btime time,Task* caller/*=nullptr*/,int result_slot/*=-1*/ )
//{
//
//	//check, if this time is to be calculated at all
//	//the first condition is, that the time is before expiration 
//	//the second condition, is that the time is on of the time slices that are planned
//	//to be calculated. The timeslices are t_offset + i*t_stepsize
//	if(t_expiration<time || (time-t_offset)%t_stepsize)
//		trace("Requested time "+toS(time)+" is not to be calculated for this option",new bad_exception());
//
//	if(time == t_expiration)
//		return payoff->get_result(t_expiration,caller,result_slot);
//	else
//		return Field::get_result(time,caller,result_slot);
//}

//############################  setup_task  ######################################
void financial_derivative::setup_task(Task* t,btime time){
	//reserve space for payoff, future option price and (up to four) fd operators
	t->reset( 6 );

	btime t_timestep=time+t_stepsize>=t_expiration ? t_expiration - time : t_stepsize;
	double timestep = global_config::get().years_per_btick*t_timestep;

	//put exercise value in 0. slot
	payoff->get_result(time,t,0);

	//the future value in 1. slot
	this->get_result(time + t_timestep, t, 1);


	//the time stepping operator in 2. slot:


	//The slack variable  given by
	//slack := (1 + mu h A) (u(k+1)+g(k+1)) - (1 - (1-mu) h A) (u(k)+g(k))
	//	can be transformed in standard form:
	//
	//		slack :=  M.u - q
	//
	//		Where
	//		M = 1 + mu h A
	//		q = (1 - (1-mu) h A) * lastsol - (1 + mu h A) g(k+1)


	double beta;
	switch(time_stepping_scheme)
	{
		case ImplicitEuler:
			//Implicit Euler
            //The LCP becomes:
            //       M = 1 + hA
            //       q = Sol(k) - M.g(k+1)
            //       sol >= 0, M.sol - q >= 0, (M.sol - q)*sol = 0
        	require_fd_operator(timestep,t,2);
			break;		
		case CrankNicolson:
			//generalized Crank-Nicolson
			//The LCP becomes:
			//   M = 1 + 0.5 h A
			//   q = ( 1 - 0.5 h A)*Sol(k) - M*g(k+1)
			//       sol >= 0, M.sol - q >= 0, (M.sol - q)*sol = 0
			require_fd_operator(beta=0.5*timestep,t,2);
			require_fd_operator(-beta,t,4);
			break;
		case ThetaMethod:
			//generalized Crank-Nicolson
			//The LCP becomes:
			//   M = 1 + theta h A
			//   q = ( 1 - (1-theta) h A)*Sol(k) - M*g(k+1)
			//       sol >= 0, M.sol - q >= 0, (M.sol - q)*sol = 0
			require_fd_operator(beta = theta*timestep,t,2);
			require_fd_operator(beta - timestep,t,4); //corresponds to : -(1.0-theta)*timestep;
			break;
		case ExplicitEuler:
			//Explicit Euler
			//The LCP becomes trivial:
			//       q = (1-hA)*Sol(k) - g(k+1)
			//       sol >= 0, sol - q >= 0, (sol - q)*sol = 0
			//
			//with the solution:     sol = Max(0,q)
			require_fd_operator(-timestep,t,2);
			break;
		default:
			trace("time_stepping_scheme ("+toS(time_stepping_scheme)+") Not implemented", 
				new invalid_argument(""));
			break;
	}
	
	//register computation routine with custom  signature
	t->computation=boost::bind(&financial_derivative::perform_computation,this,_1,timestep,time);
}

	// (operator 1+beta*A)
void financial_derivative::require_fd_operator(double beta, Task* t,int slot)
{
	weights->require(bcs.get(),beta,t,slot);
	//and another fd operator for strang symmetrized splitting
	if(strang_symmetrization)
		weights->require(bcs.get(),beta*0.5,t,slot+1);
}

//############################  perform_computation  ######################################
void financial_derivative::perform_computation( Task* t, double timestep,btime time )
{	
	//Get hold of thread local storage for Cryer LCP algo and other stuff
		time_stepping_work_storage* stor = storage.get(weights->max_length);
		stor->init_solutions(grid->N[0]);
	//current solution vector
		double1* current_solution = stor->get_next_solution();

	//get task dependencies
		const double1* payoff_data, *last_solution;
		//time stepping operators (of the form 1+beta*A)
		const ts_op_t *M[2]={nullptr}, *M_minus[2]={nullptr};

		my_cast(*t->at(0),payoff_data);
		my_cast(*t->at(1),last_solution);
		my_cast(*t->at(2),M[0]);
		if(strang_symmetrization)
			my_cast(*t->at(3),M[1]);
		if(time_stepping_scheme!=ImplicitEuler)
		{
			my_cast(*t->at(4),M_minus[0]);
			if(strang_symmetrization)
				my_cast(*t->at(5),M_minus[1]);
		}

	int strang_selector;
	
	//trace("perform calc "+toS(time)+": "+toS(this_thread::get_id()));

	
	//The number of fractional time steps to perform
	int n_fractional_steps = strang_symmetrization ?  2*weights->n_dir-1 :  weights->n_dir;
	const split_dir* dir;
	const split_dir_op_t *split_M=nullptr,*split_M_minus=nullptr;
	int dir_i;

	for(int frac_step=1;frac_step<=n_fractional_steps;frac_step++)
	{
		//determine fractional step size and fractional direction
		dir = determine_fractonial_dir(frac_step,strang_selector,dir_i);

		//get pointer to the required split time step operators
		split_M=&M[strang_selector]->at(dir_i);
		if(time_stepping_scheme!=ImplicitEuler)
			split_M_minus=&M_minus[strang_selector]->at(dir_i);

		//ITERATE THROUGH EVERY SYSTEM:
		for(uint sys=0;sys < dir->systems.size(); sys++)
		{
			/* TODO: for  constant g, there are a lot of */
			/* operations that only have to be done once (and not at every timestep) */

			//CALCULATE THE Q VARIABLE
			switch(time_stepping_scheme)
			{
				case ImplicitEuler:
					//   q = Sol(k) - M.g(k+1)
					dir->copy_2_global_to_local(sys,last_solution,	stor->q,
													payoff_data,	stor->g);
					split_M->at(sys).times_vec(stor->g,stor->q,-1.0,1.0);

					stor->cryer.solve_lcp(split_M->at(sys),stor->q,stor->sol);

					break;		
				case CrankNicolson:
					//   q = ( 1 - 0.5 h A)*Sol(k) - M*g(k+1)
				case ThetaMethod:
					//   q = ( 1 - (1-theta) h A)*Sol(k) - M*g(k+1)

					dir->copy_2_global_to_local(sys,last_solution,	stor->sol,
													payoff_data,	stor->g);
					// first step: q := M*g
					split_M->at(sys).times_vec(stor->g,stor->q,1.0,0.0);
					// second step: q :=  M_minus*sol - q
					split_M_minus->at(sys).times_vec(stor->sol,stor->q,1.0,-1.0);

					//third step: solve
					stor->cryer.solve_lcp(split_M->at(sys),stor->q,stor->sol);

					break;
				
				case ExplicitEuler:///\todo: this does not need M!
					//       q = (1-hA)*Sol(k) - g(k+1)
					//with the solution:     sol = Max(0,q)
				
					dir->copy_2_global_to_local(sys,last_solution,	stor->sol,
													payoff_data,	stor->q);
					//first step: q := M_minus*sol
					split_M_minus->at(sys).times_vec(stor->sol,stor->q,1.0,0.0);

					for(int i = 0; i<stor->q.length ; ++i ) 
						if(stor->q(i)>stor->g(i))
							stor->q(i)-=stor->g(i);
						else
							stor->q(i)=0.0;

					break;
			}

			//Copy solution in global solution vector
			dir->copy_local_to_global(sys,stor->sol,current_solution);
		}

		last_solution = current_solution;
		current_solution = stor->get_next_solution();
	}

	//save the result in the current task;
	const_cast<double1*>(last_solution)->setName("finder solution",t);
	t->result = stor->save_solution();

}


split_dir* financial_derivative::determine_fractonial_dir(int frac_step,int& strang_selector,int& dir_i)
{
	if(strang_symmetrization)
		if(frac_step>weights->n_dir){
			strang_selector =1;
			dir_i = 2*weights->n_dir-frac_step-1;
			return weights->dirs[dir_i];
		}else{
			if(frac_step==weights->n_dir)
				strang_selector = 0;
			else
				strang_selector = 1;
			dir_i = frac_step-1;
			return weights->dirs[dir_i];
		}
	else{
		strang_selector = 0;
		dir_i = frac_step-1;
		return weights->dirs[dir_i];
	}
}
