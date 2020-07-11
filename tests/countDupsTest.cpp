/*
 * countDupsTests.cpp
 *
 *  Created on: 11 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "hicupData.h"
#include "BinomData.h"
#include "Utils.h"
#include <vector>
#include <string>
#include <map>
#define BOOST_TEST_MODULE countDupsTests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(countDupsTests)

BOOST_AUTO_TEST_CASE(testCountDuplicates) {
	std::vector<Interaction> interactions;
	interactions.push_back(Interaction("chr1","chr2",12553,15273));
	interactions.push_back(Interaction("chr1","chr2",12553,15273));
	interactions.push_back(Interaction("chr1","chr2",12553,15273));
	interactions.push_back(Interaction("chr6","chrX",125585523,1063441));
	interactions.push_back(Interaction("chr6","chrX",125585523,1063441));
	interactions.push_back(Interaction("chr10","chr5",1064473,1505273));

	std::map<std::string,int> list;
	list["chr10:1064473"] = 1;
	list["chr6:125585523"] = 2;
	list["chr1:12553"] = 3;
	list["chr5:1505273"] = 1;
	list["chrX:1063441"] = 2;
	list["chr2:15273"] = 3;


	countDuplicates(interactions);

	BOOST_CHECK_EQUAL(interactions.size(), 3);

	for (auto i : interactions)
	{
		BOOST_CHECK( i.getFreq() == list[i.getInt1()]);
		BOOST_TEST( i.getFreq() == list[i.getInt2()]);
	}

}//*/

BOOST_AUTO_TEST_SUITE_END()

