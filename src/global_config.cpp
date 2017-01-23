#include "precompiled.h"
#include "global_config.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <boost/filesystem/operations.hpp>
#endif

#include "boost_date_format.h"


namespace fipster { namespace configuration {

	void global_config::reinit(const fpath base, const ptree& pt){
		FIPSTER_ASSERT(!initialized);//("WARNING: OVERWRITING the global_config");
		initialized=true;
		{
			using namespace boost::posix_time;
			using namespace boost::filesystem;
			auto v = base; v/=pt.get<string>("outputPath");
			output_path = v.string()+"/";
			// if(is_directory(output_path+"newest"))
				// rename(output_path+"newest",output_path+format_date("%Y-%m-%d %H-%M-%S"));
			output_path += "newest/";
			if(!is_directory(output_path+"newest"))
				create_directories(output_path);
		}
		test_splitting = pt.get<bool>("tests.splitting");
		test_A_o = pt.get<bool>("tests.A_o");
		test_splitting_tolerance = pt.get<double>("tests.splitting.<xmlattr>.tolerance");
		test_tridiagonally_split_grid_operator = pt.get<bool>("tests.tridiagonally_split_grid_operator");
		test_bc_inversion = pt.get<bool>("tests.bc_inversion");
		flow_graph_serial = pt.get<bool>("flow_graph.serial");
		test_gsi_gather_scatter = pt.get<bool>("tests.gsi_gather_scatter");		
		flow_graph_processors = pt.get<int>("flow_graph.processors");
		years_per_btick = pt.get<double>("years_per_btick");
		test_grid_iterators = pt.get<bool>("tests.grid_iterators");
		test_boundary_condition_iterator = pt.get<bool>("tests.boundary_condition_iterator");
	}
}}
