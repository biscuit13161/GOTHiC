/*
 * hicupCompInterTest.cpp
 *
 *  Created on: 15 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "hicupData.h"
#include "BinomData.h"
#include "Utils.h"
#include <vector>
#include <string>
#define BOOST_TEST_MODULE HiCUPCompTests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(HiCUPCompTests)

BOOST_AUTO_TEST_CASE(compInteractions) {

	Interaction Q1 = Interaction("chr1","chr2",12553,15273);
	Interaction Q2 = Interaction("chr1","chr2",12553,15273);

	BOOST_CHECK_EQUAL(Q1.getFreq(), 1);
	BOOST_CHECK_EQUAL(Q1, Q2);

}

BOOST_AUTO_TEST_SUITE_END()


