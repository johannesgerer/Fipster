#include "precompiled.h"

#include "auto_connecting_body.h"
#include "memoized_node_factory.h"
#include "DebugSink.h"

using namespace fipster;

struct expI : auto_connecting_body<expI,double,int,int> {

	expI(bool b){};

	result_sptr_t operator()(const input_t& arg)
	{

		return result_sptr_t(new double(exp((double)1)+*get<0>(arg)-2));
	}

};

struct tt : auto_connecting_body<tt,int>{
	shared_ptr<int> i;
	
	tt(int i):i(new int(i)){};

	result_sptr_t operator()(input_t){
		return i;
	}
};


//
//struct myfactory : memoized_node_factory<myfactory,int, double>
//{
//
//};
//
using namespace tbb::flow;
//
//class asd{
//public:
//	asd(int){
//		auto a = this;
//	}
//};


struct mymem : public memoized_node_factory<mymem,int,double> {
	
	sender_ptr setup_node(arg_t i)
	{
		return create_node(expI(true),
			make_tuple(create_node(tt(i+100)).get(),create_node(tt(-i+100)).get()));
	}
};


void main2()
{
	DebugSink::activate_std();
	bool ass=true;
	
	rooted_graph& g=rooted_graph::get();

	write_once_node<expI::result_sptr_t>w(g.tbb_graph);

	int count=0;

	/*auto dd=expI::sender_ptr(new source_node<shared_ptr<int>>(g.tbb_graph,[&](shared_ptr<int>& i)->bool
	{ 
		if(count==0){
			count++;
			i.reset(new int(13));return true;
		}
		return false;
	}
	));

	auto exp1 = expI::create(g,dd);
	auto exp1 = expI::create(g,make_tuple(tt::create(g,14),false);
	
	make_edge(*exp1,w);
*/
	
	
	shared_ptr<double> r1,r2 ;
	//g.attach_result_to_sink(r1,&w);
	//g.attach_result_to_sink(r2,exp1);
	
	expI::result_sptr_t r ;

	mymem asdmem;

	vector<shared_ptr<const double>> v(20);
	for(int i=0;i<10;i++){
		g.attach_result_to_sink(v[i],asdmem.get_node(i));
		//g.attach_result_to_sink(v[10+i],asdmem(9-i));
	}
	
	g.root_node.try_put(tbb::flow::continue_msg());

	//w.try_get(r);

	g.tbb_graph.wait_for_all();
		
	
	
	if(w.try_get(r))
	{
		cout<<*r<<endl;
	}
	
	
	for(int i=0;i<10;i++)
		cout<<*(shared_ptr<const double>)v[i]<<endl;
	fipster::named_counter::print();
/*
	auto ac = myWorker::create(g,make_tuple((myWorker3::sender_t)nullptr,(myWorker3::sender_t)nullptr));
	
	myWorker::node_t(*new rooted_graph(), myWorker::sender_t());
	myWorker::sender_tuple_t a;
	auto& b = *get<0>(a);*/
}
