#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <iostream>
#include <boost/exception/all.hpp>
#include <boost/program_options.hpp>
#endif

#include "DebugSink.h"
#include "configuration_parser.h"
#include "utility.h"

#include "grid_tests.h"
#include "global_config.h"

using namespace fipster;

void testmain();
vector<string> mainWithOptions(int, char*[]);

#define MAIN_CATCH_BLOCK 0

int main(int ac, char* av[])
{
#if MAIN_CATCH_BLOCK
try{
#endif
	DebugSink::activate_std();
	srand( (unsigned)time( NULL ));
	testmain();

			
	parser p(mainWithOptions(ac,av));
	DebugSink::activate_std(true,global_config::get().output_path);
	p.create_graph_and_run();
#if MAIN_CATCH_BLOCK
}
catch (...)
{
	cerr << boost::current_exception_diagnostic_information();
	// throw;
}
#endif
	return 0;
}

// old?
#if 0
#include "grid.h"
#include "fdstencils.h"
#include "fdweights.h"
#include "GridDirGenerators.h"
#include "payoff_node_factory.h"
#include "financial_derivative.h"
#include "Cryer.h"




class testC{};

void test1(string s,Task* t){
	
	cout<<s<<": asdasdasd!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
}

void start(Call* c, Task* t,financial_derivative* ff){
	((Field*)c)->get_result(5,t,0);
	c->get_result(10,t,1);
	c->get_result(10,t,0);
	c->get_result(-10,t,2);	
	ff->get_result(0,t,2);
	t->ready();
}


int main2()
{

	DebugSink::activate();

	propery_collection<boundary_conditions> bcs("boundaryConditions.xml");
	//propery_collection<boundary_conditions>::propCollectionBase* b = &bcs;

	int n=10,m=3;

	Grid g(2 ,	new HomogenousCoords(6,0.0,1.0),new HomogenousCoords(11,0.0,1.0),
					new ExpCoords(6,-0.2,0.1,1),new ExpCoords(5,-0.2,0.1,1));

	g.test_grid_interateboundary_canonical_order<1>();

	fdweights* t1 = g.get_fdweights(StencilOneDiagonal,new BlackScholesIsometric(0.05,0.2,0.1),6,0);
	
	g.checkiterdirstartpoints<2>(t1->stencil);
	/*fdweights* t2 = g.get_fdweights(StencilNoDiagonals,t1->pde,6);
	fdweights* t3 = g.get_fdweights(StencilOneDiagonal,t1->pde,6);
	fdweights* t4 = g.get_fdweights(StencilTwoDiagonals,t1->pde,6);*/
	
	Call c( -.15 , 1 );
	c.setGrid(&g);
	Task t(boost::bind(test1,"cool",_1),3);
	financial_derivative ff(&c,t1,bcs["testBC"],1000,10);
	financial_derivative ff2(&ff,t1,bcs["testBC"],1000,10,ImplicitEuler);
	ff2.setName("fin_der2");
	ff.setName("fin_der1");
	Task::dispatch_wrapped(boost::bind(start,&c,&t,&ff2));

	Task::finish_computations();
	
	//c.perform_computation(234.0);
	//financial_derivative f(&c,t1);
	cout<<"END END END END END END END END !"<<endl;
	cout<<Task::n_unfinished_max<<endl;

	return 0;
}


#endif
