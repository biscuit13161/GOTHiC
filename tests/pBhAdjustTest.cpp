/*
 * pBhAdjustTest.cpp
 *
 *  Created on: 21 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "padjust.h"
#define BOOST_TEST_MODULE pBhAdjustTest
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(pBhAdjustTest)

BOOST_AUTO_TEST_CASE(first) {

	double P = 6.079281e-10;
	double n = 28679;
	double o = pBhAdjust(P, n);

	BOOST_CHECK_CLOSE(1.743477e-05, o, 5e-7 ); // tolerance is a PERCENTAGE
}

BOOST_AUTO_TEST_SUITE_END()


