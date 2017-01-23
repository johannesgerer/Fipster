#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <string>
#include <unordered_set>
#include <memory>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <tbb/task_scheduler_init.h>
#endif

#include "configuration_parser.h"
#include "explicit_fields.h"
#include "spacetime_field.h"
#include "grid.h"
#include "solution.h"
#include "expectation_value_node_factory.h"
#include "boundary_conditions_config.h"
#include "finite_difference_weights_config.h"
#include "named_counter.h"
#include "exceptions.h"
#include "task_queue.h"
#include "global_config.h"
#include "option.h"


using namespace fipster;
using namespace fipster::configuration;
using namespace tbb;

// #################################################
// ##############   PHASE 1:			############
// ##############   CONFIGURATION       ############
// #################################################

string getId(const ptree& pt,string tag)
{
	return betterGet<string>(pt,"<xmlattr>.id",tag);
}

void parser::parseTree(const fpath base, const ptree& pt)
{
	//main loop
	for(auto& it : pt){

		auto& node = it.second;
		auto& tag = it.first;
		
		// file #######################################################
		if(tag == "file"){
			auto p = base; p/=node.get_value<string>();
			parseFile(p.string());

			// boundaryConditions #######################################################
		}else if(tag == "boundaryConditions"){

			objs->bcs_configs[getId(node,tag)].reset(new _boundary_conditions::config_t(node));

			// field #######################################################
		}else if(tag == "field"){

			string type = node.get<string>("<xmlattr>.type");

			auto& new_field = objs->spacetime_fields[getId(node,tag)];

			if(type == "explicit_field")
				new_field.reset(new explicit_fields::node_factory(node,objs));
			else
				BOOST_THROW_EXCEPTION
					(runtime_error("'"+ type +"' is not a valid spacetime_field"));

			// options #######################################################
		}else if(tag == "option"){

			objs->options[getId(node,tag)].reset
				(new option_node_factory(make_shared<option_config_t>(node),objs));

			// expectation #######################################################
		}else if(tag == "expectation"){

			objs->expectations[getId(node,tag)].reset
				(new expectation_values::node_factory
				 (make_shared<expectation_values::config_t>(node),objs));

			// grid #######################################################
		}else if(tag == "grid"){

			objs->grids[getId(node,tag)].reset(new Grid(node));

			// solution #######################################################
		}else if(tag == "solution"){

			objs->solutions.push_back(solution(node,objs));

			// finiteDifferenceWeights #######################################################
		}else if(tag == "finiteDifferenceWeights"){

			objs->fdweights_configs[getId(node,tag)].reset
				(new finite_difference_weights::config_t(node,*objs));

			// <xml...> #######################################################
		}else if(tag == "globalConfig"){
			global_config::get(false).reinit(base, node);
			// <xml...> #######################################################
		}else if(!is_xml_special(tag))
			BOOST_THROW_EXCEPTION(runtime_error("'"+ tag +"' is not an supported configuration tag"));
	}
}

			

void parser::parseFile(string path_to_xml){
#define CATCH_IT 1
#ifdef  CATCH_IT
	try{
#endif
		//try to insert new file
		auto ins = processed_files.insert(path_to_xml);
		//if file is already processed, return
		if(!ins.second) return;

		ptree propertytree;
		read_xml(path_to_xml,propertytree);
		parseTree(fpath(path_to_xml).parent_path()
							,propertytree.get_child("jML"));

#ifdef  CATCH_IT
	}catch(boost::exception& e){
		e << boost::errinfo_file_name(path_to_xml);
		cerr << boost::current_exception_diagnostic_information();
		abort();
		throw;
	}
#endif

}

//Ctor
parser::parser(vector<string> paths_to_xml):
	objs(new configured_objects)
{
	for(auto &f : paths_to_xml)
		parseFile(f);
}


// #################################################
// ##############   PHASE 2:			############
// ##############   Node Factory        ############
// #################################################

void parser::create_flow_graph_nodes(){
	// try{
	//start the global node setup using the solutions' setup_nodes members
	for(auto it=begin(objs->solutions);it!=end(objs->solutions);it++)
		it->setup_nodes();
	// }catch(...){//boost::exception& e){
	// 	// e << boost::errinfo_file_name(path_to_xml);
	// 	cerr << boost::current_exception_diagnostic_information();
	// 	throw;
	// }

	named_counter::print();	

	objs->print_all_map_sizes();
}

// #################################################
// ##############   PHASE 3:			############
// ##############   Node Factory        ############
// #################################################

	struct sleeping_body{	
		void operator()(){
			thread_logger()<<"GRAPH  | Going to eternal sleep"<<endl;
			this_tbb_thread::sleep(tick_count::interval_t(10000000000.0));
		}
	};

void parser::create_graph_and_run()
{
		int processors=global_config::get().flow_graph_processors;
		if(processors == 0)
			processors = task_scheduler_init::automatic;
		
		bool serial=global_config::get().flow_graph_serial;// || processors==1;

		thread_logger()<<"Parser  | create_graph_and_run ("<<(serial?"serial":"parallel")<<")"<<endl;

		if(serial){ // \todo fails with 
			FIPSTER_ASSERT(0);
			/*  Assertion t->prefix().ref_count==0 failed on line 500 of file ../../src/tbb/custom_scheduler.h
					Detailed description: Task still has children after it has been executed */
			//override default scheduler
			task_scheduler_init serial_scheduler(1);
			//workaround: TTB always starts one additional worker, that has
			//to be sent to sleep manually:
			rooted_graph::get().tbb_graph.run(sleeping_body());
			rooted_graph::get().tbb_graph.decrement_wait_count();
			//let this thread sleep to give the other worker thread a 
			//chance to run sleeping_body before "try_put"
			this_tbb_thread::sleep(tick_count::interval_t(0.1));
		}
			
		task_scheduler_init serial_scheduler(serial ? 1 : processors);
			
		create_flow_graph_nodes();
		start_and_check();
}


void parser::start_and_check()
{
	{
		using namespace boost;
		using namespace boost::filesystem;

		//copy config files
		if(false)
			for(auto it=processed_files.begin();
					it!=processed_files.end();it++){
				fpath p(global_config::get().output_path+*it);
				create_directories(p.parent_path());
				copy_file(*it, p, copy_option::overwrite_if_exists);
			}

		if(0)
			//print state grid
			for(auto it=objs->grids.begin();it!=objs->grids.end();it++){
				std::ofstream out(global_config::get().output_path+it->first+".dat");
				out<<*it->second;
				out.close();
			}
		
	}

	rooted_graph::get().start_wait_for_all();

	//check if every solution sink is filled
	for(auto& it : objs->solutions)
		{
			if(it.not_done.empty()) continue;
			thread_logger()<<"!! Graph not connected properly: The following solution sink was not filled:"<<endl;
			thread_logger()<<it<<", at times:"<<endl;
			for(auto t : it.not_done)
				thread_logger()<<t<<endl;
		}
}

btime fipster::configuration::getBtime
(const ptree& pt, const string& s, const string& c){
	double d = betterGet<double>(pt,s,c);
	btime i = round(d);
	if((double)i != d) BOOST_THROW_EXCEPTION
											 (boost::property_tree::ptree_bad_data
												("Unable to convert to integer (via double)\nData: "+
												 toS(d)+"\nContext: "+toS(c),d));
	return i;
}


