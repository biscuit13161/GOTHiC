/*
 * importHicupTests.cpp
 *
 *  Created on: 12 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "../src/hicupData.h"
#include "../src/BinomData.h"
#include "../src/Utils.h"
#include <vector>
#include "tbb/concurrent_vector.h"
#include <fstream>
#define BOOST_TEST_MODULE importHicupTests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(importHicupTests)

BOOST_AUTO_TEST_CASE(constructor) {
	tbb::concurrent_vector<Interaction> interactions;

	importHicupTxt("../../example_data/example_Test_Data.txt", interactions);

	BOOST_CHECK_EQUAL(interactions.size(), 2);
	BOOST_CHECK_EQUAL(interactions[0].getInt1(), "chr17_gl000204_random:51677");
	BOOST_CHECK_EQUAL(interactions[0].getInt2(), "chr2:155678032");
	BOOST_CHECK_EQUAL(interactions[1].getInt1(), "chr3:101663949");
	BOOST_CHECK_EQUAL(interactions[1].getInt2(), "chr2:230303040");
}

BOOST_AUTO_TEST_SUITE_END()
