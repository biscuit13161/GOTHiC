/*
 * binomTestTest.cpp
 *
 *  Created on: 26 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "binomTest.h"
#include "hicupData.h"
#include "BinomData.h"
#include "pbinom.h"
#include <vector>
#define BOOST_TEST_MODULE binomTestTest
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(binomTestTest)

BOOST_AUTO_TEST_CASE(first) {

	double P1 = binomTest(2, 28679, 6.079281e-10, "greater");
	double P2 = binomTest(2, 28679, 6.079281e-10, "less");

	BOOST_CHECK_CLOSE(1.519785e-10, P1, 5e-5);
	BOOST_CHECK_CLOSE(1, P2, 5e-5);
}

BOOST_AUTO_TEST_SUITE_END()
