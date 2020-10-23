/*
 * randomTest.cpp
 *
 *  Created on: 22 Oct 2020
 *      Author: rich
 */

#include "random.h"
#define BOOST_TEST_MODULE randomTests
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
	sources.push_back(Interaction("chr2","chr1",12553,15273));
	sources.push_back(Interaction("chrX","chr7",1255,1020));
	sources.push_back(Interaction("chrX","chrX",1255,1255));

	tbb::concurrent_vector<Interaction> interactions2;
	interactions2.push_back(Interaction("chr2","chr1",12553,15273));
	interactions2.push_back(Interaction("chr1","chr1",17753,17753));
	interactions2.push_back(Interaction("chr1","chr1",17753,15273));
	interactions2.push_back(Interaction("chrX","chr7",1255,1020));
	interactions2.push_back(Interaction("chrX","chrX",1255,1255));

	tbb::concurrent_vector<Interaction> interactions3;
	interactions3.push_back(Interaction("chr2","chr1",12553,15273));
	interactions3.push_back(Interaction("chr1","chr1",17753,17753));
	interactions3.push_back(Interaction("chr1","chr1",17753,15273));
	interactions3.push_back(Interaction("chrX","chr7",1255,1020));
	interactions3.push_back(Interaction("chrX","chrX",1255,1255));

	BOOST_CHECK_EQUAL(interactions.size(), 5);
	BOOST_CHECK_EQUAL(sources.size(), 3);
	BOOST_CHECK_EQUAL(interactions2.size(), 5);
	BOOST_CHECK_EQUAL(interactions3.size(), 5);

	RandomSubset(interactions, sources);
	RandomChoose(sources, interactions2);
	RandomChoose(interactions3, sources);

	BOOST_CHECK_EQUAL(interactions.size(), 3);
	BOOST_CHECK_EQUAL(interactions2.size(), 3);
	BOOST_CHECK_EQUAL(interactions3.size(), 3);
}

BOOST_AUTO_TEST_SUITE_END()
