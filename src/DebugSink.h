#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <iostream>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/tee.hpp>
#include <string>
#include <fstream>
#endif



namespace fipster { namespace _DebugSink {


	using namespace std;
	namespace io = boost::iostreams;


	//buffer class 
	struct my_buffer
	{
		typedef char char_type;
		typedef io::sink_tag category;

		shared_ptr<ofstream> ofs;

		my_buffer(shared_ptr<ofstream> ofs);

		std::vector<char> buffer;

		// write #####################
		std::streamsize write(const char *s, std::streamsize n);
	};


/** its static "activate()" member replaces std::cout's buffer with a TeeBuffer, 
 that copies the stream into a file and a DebugSink, that will write everything to 
 Debug Window using OutputDebugStringA.
 \note there is no deactivate, so if the object is destroyed cout is broken
*/
struct DebugSink{

	typedef io::tee_device<my_buffer, std::streambuf> TeeDevice;
	typedef io::stream_buffer<TeeDevice>  TeeBuf;
	typedef shared_ptr<TeeDevice> tee_ptr;
	typedef shared_ptr<TeeBuf> buf_ptr;

	tee_ptr tee;
	buf_ptr buf;

	shared_ptr<ofstream> ofs;

	// activate #####################
	static void activate_std(bool save_output=false,string path="");
	
	void set_path(string path);

	// custom #####################
	DebugSink(ostream& os);

	~DebugSink();

};


}


using _DebugSink::DebugSink;

}

