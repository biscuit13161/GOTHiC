/*
 * sysTests.cpp
 *
 *  Created on: 12 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "BinomData.h"
#include "hicupData.h"
#include <vector>
#include "tbb/concurrent_vector.h"
#define BOOST_TEST_MODULE sortPositionsTest
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(sortPositionsTest)

BOOST_AUTO_TEST_CASE(first) {

	tbb::concurrent_vector<Interaction> interactions;
	interactions.push_back(Interaction("chr2","chr1",12553,15273));
	interactions.push_back(Interaction("chr1","chr1",17753,15273));
	interactions.push_back(Interaction("chrX","chr7",1255,1020));

	int iSize = interactions.size();
	std::vector<halfInteraction> sources;
	std::vector<halfInteraction> targets;
	for (int i = 0; i < iSize; i++)
	{
		sources.push_back(halfInteraction(interactions[i].getChr1(), interactions[i].getLocus1() ));
		targets.push_back(halfInteraction(interactions[i].getChr2(), interactions[i].getLocus2() ));
	}

	sortPositions(interactions, iSize, sources, targets);

	BOOST_CHECK_EQUAL(interactions[0].getInt2(), "chr2:12553");
	BOOST_CHECK_EQUAL(interactions[1].getInt1(), "chr1:15273");
	BOOST_CHECK_EQUAL(interactions[2].getInt2(), "chrX:1255");

}

BOOST_AUTO_TEST_SUITE_END()
