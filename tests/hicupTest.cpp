/*
 * hicupTests.cpp
 *
 *  Created on: 11 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "hicupData.h"
#include "BinomData.h"
#include "Utils.h"
#include <vector>
#include <string>
#define BOOST_TEST_MODULE HiCUPTests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(HiCUPTests)

BOOST_AUTO_TEST_CASE(constructor) {
	Interaction A = Interaction("chr1","chr2",12553,15273,3);

	BOOST_CHECK_EQUAL(A.getFreq(), 3);
	BOOST_TEST(A.getInt1() == "chr1:12553");
	BOOST_TEST(A.getInt2() == "chr2:15273");
}

BOOST_AUTO_TEST_SUITE_END()
