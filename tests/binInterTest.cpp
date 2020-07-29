/*
 * binInterTests.cpp
 *
 *  Created on: 13 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "BinomData.h"
#include "hicupData.h"
#include "SetupData.h"
#include "Utils.h"
#include <vector>
#include "tbb/concurrent_vector.h"
#include <map>
#define BOOST_TEST_MODULE binInterTests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(binInterTests)

BOOST_AUTO_TEST_CASE(constructor) {

	tbb::concurrent_vector<Interaction> interactions;
	SetupData setupValues;
	setupValues.setRes(1000);
	interactions.push_back(Interaction("chr2","chr1",12553,15273));
	interactions.push_back(Interaction("chr1","chr1",17753,15273));
	interactions.push_back(Interaction("chrX","chr7",1255,1020));

	std::map<std::string,int> list;
	list["chr1:17753"] = 1;
	list["chr2:12553"] = 1;
	list["chr1:15273"] = 1;
	list["chrX:1255"] = 1;
	list["chr7:1020"] = 1;

	std::map<std::string,int> list2;
	list2["chr1:17000"] = 1;
	list2["chr2:12000"] = 1;
	list2["chr1:15000"] = 1;
	list2["chrX:1000"] = 1;
	list2["chr7:1000"] = 1;

	binInteractions(interactions, setupValues);


	BOOST_CHECK_NE(list[interactions[0].getInt1()], 1);
	BOOST_CHECK_EQUAL(list2[interactions[0].getInt1()], 1);
	BOOST_CHECK_NE(list[interactions[2].getInt2()], 1);
	BOOST_CHECK_EQUAL(list2[interactions[1].getInt2()], 1);

}

BOOST_AUTO_TEST_SUITE_END()

