/*
 * calcFreqTest.cpp
 *
 *  Created on: 15 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "BinomData.h"
#include "hicupData.h"
#include <vector>
#include "tbb/concurrent_vector.h"
#define BOOST_TEST_MODULE calcFreqTest
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(calcFreqTest)

BOOST_AUTO_TEST_CASE(first) {

	tbb::concurrent_vector<Interaction> interactions;
	interactions.push_back(Interaction("chr2","chr1",12553,15273,6));
	interactions.push_back(Interaction("chr1","chr1",17753,15273,2));
	interactions.push_back(Interaction("chrX","chr7",1255,1020,1));

	std::map<std::string,int> cov;
    double coverage = 0;
    int max = 0;
    int numberOfReadPairs = 0;

    calcFreq(interactions, cov, numberOfReadPairs, coverage, max);

    //ASSERT_TRUE(true);
    BOOST_CHECK_EQUAL(coverage,18);
    BOOST_CHECK_EQUAL(max, 8);
    //EXPECT_TRUE(numberOfReadPairs == 7);
    BOOST_CHECK_EQUAL(cov["chr1:17753"], 2);
    BOOST_CHECK_EQUAL(cov["chr2:12553"], 6);
    BOOST_CHECK_EQUAL(cov["chr1:15273"], 8);
    BOOST_CHECK_EQUAL(cov["chrX:1255"], 1);
    BOOST_CHECK_EQUAL(cov["chr7:1020"], 1);

}

BOOST_AUTO_TEST_SUITE_END()
