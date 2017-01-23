#pragma once

#ifndef AUTOMATIC_PRECOMPILATION
#include <exception>
#include <boost/throw_exception.hpp>
#include <boost/exception/all.hpp>
#endif


/* Example */
#if 0

try{

	//this try block is needed, if additional
	//information should be added to boost::exceptions
	try{

		//this try block is needed, if plain std::exceptions should
		//be wrapped in boost::exception
		//(which is only needed if information should be added to them)
		try
		{
			//preferred (this will generate File(Line) info to
			//be used when visual studio tools do not help
			//(e.g. inside tbb tasking)
			BOOST_THROW_EXCEPTION(runtime_error("asd"))

			//legacy/third-party code
			throw runtime_error("asd");

		}
		//this macro is needed, if plain std::exceptions should
		//be wrapped in boost::exception
		//(which is only needed if information should be added to them)
		CATCH_BOOST_STD_EXCEPTION ;

	//this catch block is needed, if additional
	//information should be added to boost::exceptions
	}catch(boost::exception& e){
		
		e << boost::errinfo_file_name(path_to_xml)
			<< boost::errinfo_errno(no)
			<< boost::errinfo_api_function("asd")
			<< [custom made: boost::error_info ](...);
		throw;
	}

//this catch block will ouput information
//about any kind of exception. 
catch (...)
{
	cerr << boost::current_exception_diagnostic_information();
	throw;
}
#endif

#define CATCH_BOOST_STD_EXCEPTION catch (std::exception& e)\
		{									\
			if(boost::current_exception_cast<boost::exception>())	\
				throw;						\
			else							\
				boost::throw_exception(e);	\
		}

/** \todo make sure, this is used in the correct way: 
// do not use assert, as does not work in production mode.
// If you want a test only to run in (my custom) debug mode, make it
// optional on FIPSTER_NDEBUG
// for ALL TESTS use on of these two. They use assert/abort, if available
// and throw an exception otherwise. 
// the difference is only the condition in FIPSTER_ASSERT

GCC: diabled assert version in FIPSTER_THROW_EXCEPTION \todo: bring it back
*/
#ifdef NDEBUG 
#define FIPSTER_THROW_EXCEPTION(x) BOOST_THROW_EXCEPTION(x)
#else
#define FIPSTER_THROW_EXCEPTION(x) assert((x,0))
#endif

#ifdef NDEBUG
#define FIPSTER_ASSERT(x) (x?[](){}():(BOOST_THROW_EXCEPTION(std::runtime_error("Assertion failed: "+ string(#x) + "\n at " + __FILE__ + ":" +to_string(__LINE__)))))
#else
#define FIPSTER_ASSERT(x) assert(x)
#endif

