
#undef _GLIBCXX_DEBUG
// http://stackoverflow.com/questions/19729036/boost-program-options-wont-work-with-glibcxx-debug

#include <boost/program_options.hpp>
#include <vector>
#include <debug/vector>
#include <iostream>
using namespace std;
using namespace boost;
using namespace boost::program_options;

template<class A>
void usage(const char* p,A desc){
	cerr << "Usage: " << p << " [options] [config files] " << endl  << desc << endl;
}

// typedef vector<string> rv; //using this typedef in place of vector<string> strangely breaks stuff below

#ifdef NDEBUG
typedef vector<string> dv;
#else	 
typedef __gnu_debug::vector<string> dv;
#endif

dv mainWithOptions(int ac, char* av[])
{
	variables_map vm;
	options_description desc("Options");
	//http://www.boost.org/doc/libs/1_58_0/doc/html/program_options/tutorial.html
	
	positional_options_description p;
		

	desc.add_options()
		("help,h", "produce help message")
	  ("config", value< vector<string> >()->default_value({"config.xml"},"config.xml"), "configuration files (XML)");
		;
	
	
	p.add("configuration", -1);
	
	store(command_line_parser(ac, av).
				options(desc).positional(p).run(), vm);
	notify(vm);

	if(vm.count("help")){
		usage(*av,desc);
		exit(0);
	}
	if(!vm.count("config")){
		cerr<<"No configuration file given"<<endl<<endl;
		usage(*av,desc);
		exit(1);
	}
	
	auto v=vm["config"].as<vector<string> >();
	return v;
}

