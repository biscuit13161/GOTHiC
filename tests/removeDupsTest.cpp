/*
 *  removeDupsTest.cpp
 *
 *  Created on: 13 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */


#include "BinomData.h"
#include "hicupData.h"
#include <vector>
#include "tbb/concurrent_vector.h"
#define BOOST_TEST_MODULE SystemCheck
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(SystemCheck)

BOOST_AUTO_TEST_CASE(systemCheck) {

	tbb::concurrent_vector<Interaction> interactions;
	interactions.push_back(Interaction("chr2","chr1",12553,15273));
	interactions.push_back(Interaction("chr1","chr1",17753,17753));
	interactions.push_back(Interaction("chr1","chr1",17753,15273));
	interactions.push_back(Interaction("chrX","chr7",1255,1020));
	interactions.push_back(Interaction("chrX","chrX",1255,1255));

	tbb::concurrent_vector<Interaction> sources;
	
	removeDuplicates(interactions, sources);

	BOOST_CHECK_EQUAL(interactions.size(), 5);
	BOOST_CHECK_EQUAL(sources.size(), 3);
	BOOST_CHECK_NE(sources[0].getInt1(), sources[0].getInt2());
	BOOST_CHECK_NE(sources[1].getInt1(), sources[1].getInt2());
	BOOST_CHECK_NE(sources[2].getInt1(), sources[2].getInt2());
}

BOOST_AUTO_TEST_SUITE_END()
