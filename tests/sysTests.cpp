/*
 * sysTests.cpp
 *
 *  Created on: 12 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#define BOOST_TEST_MODULE SystemCheck
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(SystemCheck)

BOOST_AUTO_TEST_CASE(systemCheck) {
	int count = 3;
	int query = 4;

	BOOST_CHECK_NE(count, query);
	BOOST_CHECK_EQUAL(count, 3);

}

BOOST_AUTO_TEST_SUITE_END()
