/*
 * IHWTest.cpp
 *
 *  Created on: 7 Jul 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include <string>
#define BOOST_TEST_MODULE IHW
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(IHWCheck)

BOOST_AUTO_TEST_CASE(first)
{
		std::string cmd = std::string ("Rscript --vanilla -e 'if (!require(IHW)) {quit(status = 11)} ' ");

		int sys = system(cmd.c_str());

		BOOST_TEST( sys == 0 );
}

BOOST_AUTO_TEST_SUITE_END()
