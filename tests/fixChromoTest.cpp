/*
 * fixChromoTest.cpp
 *
 *  Created on: 15 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "hicupData.h"
#include "BinomData.h"
#include "Utils.h"
#include <vector>
#define BOOST_TEST_MODULE SystemCheck
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(SystemCheck)

BOOST_AUTO_TEST_CASE(systemCheck) {

	BOOST_CHECK_EQUAL(fixChromosomeNames("1"), "chr1");
	BOOST_CHECK_EQUAL(fixChromosomeNames("X"), "chrX");
	BOOST_CHECK_EQUAL(fixChromosomeNames("chr2"), "chr2");
	BOOST_CHECK_EQUAL(fixChromosomeNames("CHR3"), "chr3");

}

BOOST_AUTO_TEST_SUITE_END()
