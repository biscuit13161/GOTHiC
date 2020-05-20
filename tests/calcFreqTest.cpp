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

    //ASSERT_TRUE(true);
    ASSERT_TRUE(coverage == 18);
    EXPECT_TRUE(max == 8);
    EXPECT_TRUE(cov["chr1:17753"] == 2);
    EXPECT_TRUE(cov["chr2:12553"] == 6);
    EXPECT_TRUE(cov["chr1:15273"] == 8);
    EXPECT_TRUE(cov["chrX:1255"] == 1);
	EXPECT_TRUE(cov["chr7:1020"] == 1);

}

