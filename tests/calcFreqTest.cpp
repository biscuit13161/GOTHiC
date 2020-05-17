/*
 * calcFreqTest.cpp
 *
 *  Created on: 15 May 2020
 *      Author: rich
 */

#include "BinomData.h"
#include "hicupData.h"
#include <gtest/gtest.h>
#include <vector>

TEST(calcFreqTest, first)
{
	std::vector<Interaction> interactions;
	interactions.push_back(Interaction("chr2","chr1",12553,15273,6));
	interactions.push_back(Interaction("chr1","chr1",17753,15273,2));
	interactions.push_back(Interaction("chrX","chr7",1255,1020,1));

	std::map<std::string,int> cov;
    double coverage = 0;
    int max = 0;

    calcFreq(interactions, cov, coverage, max);

    ASSERT_TRUE(true);
    EXPECT_TRUE(coverage == 18);
    EXPECT_TRUE(max == 8);

}

