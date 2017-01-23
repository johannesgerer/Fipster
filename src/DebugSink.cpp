#include "precompiled.h"

#ifndef AUTOMATIC_PRECOMPILATION
#include <iostream>
#include <memory>
#include <vector>
#endif

//gcc: #ifndef NOMINMAX
// #define NOMINMAX
// #endif
// #include <windows.h>
// #undef NOMINMAX

#include "DebugSink.h"

namespace fipster { namespace _DebugSink {

	

	std::streamsize my_buffer::write(const char *s, std::streamsize n)
	{
		unsigned int size=buffer.size();
		buffer.reserve(static_cast<unsigned int>(n)+size);
		buffer.resize(size-1);
		buffer.insert(buffer.end(),s, s+n);
		buffer.push_back(0); // we must null-terminate for WINAPI
		//gcc: \todo replace:  OutputDebugStringA(&buffer[size-1]);		
		if(ofs->is_open()){
			(*ofs)<<&buffer[0];
			ofs->flush();
			buffer.resize(1);
		}
		return n;
	}

	my_buffer::my_buffer( shared_ptr<ofstream> ofs)
		: ofs(ofs)
	{
		buffer.resize(1);
	}

	DebugSink::DebugSink(ostream& os)
		:ofs(new ofstream())
	{
#if defined(_WIN32) || defined(_WIN32__)  
		tee = tee_ptr(new TeeDevice(my_buffer(ofs), *os.rdbuf()));
		buf = buf_ptr(new TeeBuf(*tee));
		os.rdbuf(buf.get());
#endif // _WIN32
	}

	void DebugSink::activate_std(bool save_output,string path)
	{
		static DebugSink singleton2(cout);
		static DebugSink singleton1(cerr);
		if(save_output){
			singleton2.set_path(path+"output-cout");
			singleton1.set_path(path+"output-cerr");
		}
	}

	void DebugSink::set_path( string path )
	{
		ofs->open(path+".txt");
	}

	DebugSink::~DebugSink()
	{
		ofs->close();	}

}}
