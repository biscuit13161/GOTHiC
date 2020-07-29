/*
 * findTransTest.cpp
 *
 *  Created on: 15 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "BinomData.h"
#include "hicupData.h"
#include <vector>
#include "tbb/concurrent_vector.h"
#define BOOST_TEST_MODULE findTransTest
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(findTransTest)

BOOST_AUTO_TEST_CASE(test) {

	tbb::concurrent_vector<Interaction> interactions;
	interactions.push_back(Interaction("chr2","chr1",12553,15273));
	interactions.push_back(Interaction("chr1","chr1",17753,15273));
	interactions.push_back(Interaction("chrX","chr7",1255,1020));

	tbb::concurrent_vector<Interaction> binned_df_filtered;

	findTrans(binned_df_filtered, interactions);

	BOOST_CHECK_EQUAL(binned_df_filtered.size(), 2);
}


BOOST_AUTO_TEST_SUITE_END()
