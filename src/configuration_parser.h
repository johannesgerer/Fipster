#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <string>
#include <unordered_set>
#include <boost/property_tree/ptree.hpp>
#include <pretty_printer.h>
#endif

#include "configured_objects.h"
#include "utility.h"
#include "grid_indices.h"


namespace fipster { namespace configuration {

	using namespace std;
	using boost::property_tree::ptree;

// ##############  main parser  ##################	
struct parser{


	shared_objs_t objs;
	unordered_set<string> processed_files;

	void parseTree(const fpath base, const ptree& pt);

	void parseFile(string path_to_xml);

	//Ctor
	parser(vector<string> path_to_xml);

	//Start flow graph nodes creation
	void create_flow_graph_nodes(); 

	void create_graph_and_run(); 

	//Start computation and wait for all nodes to finish
	void  start_and_check();
};

		


template<class A>
A betterGet(const ptree& pt, const string& s, const string& c){
	try{ 
		if(s=="") return pt.get_value<A>();
		return pt.get<A>(s); 
	}catch(boost::property_tree::ptree_bad_data& e){
    auto d = e.template data<string>();
		BOOST_THROW_EXCEPTION
			(boost::property_tree::ptree_bad_data
			 (string(e.what())+"\nData: "+
			 d+"\nContext: "+toS(c)+" "+toS(s),d));
	}catch(boost::property_tree::ptree_bad_path& e){
    auto p = e.template path<boost::property_tree::path>();
		BOOST_THROW_EXCEPTION
			(boost::property_tree::ptree_bad_path
			 (string(e.what())+"\nContext: "+c,p));
	}
}

template<class T>
T& betterChild(T& pt, const string& s, const string& c){
	try{ return pt.get_child(s); }
	catch(boost::property_tree::ptree_bad_path& e){
    auto p = e.template path<boost::property_tree::path>();
		BOOST_THROW_EXCEPTION
			(boost::property_tree::ptree_bad_path
			 (string(e.what())+"\nContext: "+c,p));
	}
}
		
btime getBtime(const ptree& pt, const string& s, const string& c);

//extract a state for the current grid (dimension) from a config ptree
template<class T>
Eigen::Matrix<T,Eigen::Dynamic,1> read_vec(string tag, const ptree& pt, string c){
	vector<T> t;
	for(auto j : pt){
		if(is_xml_special(j.first)) continue;
		if(j.first!=tag)
			BOOST_THROW_EXCEPTION
				(runtime_error("'"+ j.first +"' should be '"+tag+"' when parsing a state ("+c+")"));
		t.push_back(betterGet<T>(j.second,"",c+"."+tag));
	}
	return Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1> >(t.data(),t.size());
}

template<class T>
vector<Eigen::Matrix<T,Eigen::Dynamic,1> > 
read_vecs(string tag0, string tag, const ptree& pt, string c){
	vector<Eigen::Matrix<T,Eigen::Dynamic,1> > res;
	for(auto i : pt){
		if(is_xml_special(i.first) || i.first!=tag0) continue;
		// FIPSTER_ASSERT(i.first==tag0); why should'nt there be other
		// elements allowed, that are simply ignored
		res.push_back(read_vec<T>(tag,i.second,c+".at"));
	}
	return res;
}
	}

using configuration::parser;
using configuration::betterGet;
using configuration::betterChild;
using configuration::getBtime;
using configuration::read_vec;
using configuration::read_vecs;

}
