#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE  MyTests
#include <boost/test/unit_test.hpp>
#include "Math.h"

using namespace fipster;

BOOST_AUTO_TEST_CASE( offset_tests )
{
	BOOST_CHECK_EQUAL ( next_offset_multiple(-100,-100,1), -99 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(-100,-100,50), -50 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(-101,-100,50), -100 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(-150,-100,50), -100 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(150,-100,50), 200 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(150,1050,50), 200 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(149,-100,50), 150 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(151,-100,50), 200 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(152,152,50), 202 );

	BOOST_CHECK_EQUAL ( next_offset_multiple(-100,-100,-1), -101 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(-100,-100,-50), -150 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(-101,-100,-50), -150 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(-150,-100,-50), -200 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(150,-100,-50), 100 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(149,-100,-50), 100 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(151,-100,-50), 150 );
	BOOST_CHECK_EQUAL ( next_offset_multiple(152,152,-50), 102 );
}
