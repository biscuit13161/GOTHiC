/*
 *  removeDupsTest.cpp
 *
 *  Created on: 13 May 2020
 *      Author: rich
 */


#include "BinomData.h"
#include "hicupData.h"
#include <gtest/gtest.h>
#include <vector>


TEST(removeDupsTest, first)
{
	std::vector<Interaction> interactions;
	interactions.push_back(Interaction("chr2","chr1",12553,15273));
	interactions.push_back(Interaction("chr1","chr1",17753,17753));
	interactions.push_back(Interaction("chr1","chr1",17753,15273));
	interactions.push_back(Interaction("chrX","chr7",1255,1020));
	interactions.push_back(Interaction("chrX","chrX",1255,1255));

	std::vector<Interaction> sources;
	
	removeDuplicates(interactions, sources);

	ASSERT_TRUE(interactions.size() == 5);
	EXPECT_TRUE(sources.size() == 3);
	EXPECT_FALSE(sources[0].getInt1() == sources[0].getInt2());
	EXPECT_FALSE(sources[1].getInt1() == sources[1].getInt2());
	EXPECT_FALSE(sources[2].getInt1() == sources[2].getInt2());
}

