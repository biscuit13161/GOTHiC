/*
 * binInterTests.cpp
 *
 *  Created on: 13 May 2020
 *      Author: rich
 */

#include "BinomData.h"
#include "hicupData.h"
#include "Utils.h"
#include <gtest/gtest.h>
#include <vector>


TEST(binInterTests, constructor)
{
	std::vector<Interaction> interactions;
	int res = 1000;
	interactions.push_back(Interaction("chr2","chr1",12553,15273));
	interactions.push_back(Interaction("chr1","chr1",17753,15273));
	interactions.push_back(Interaction("chrX","chr7",1255,1020));

	binInteractions(interactions, res);

	ASSERT_TRUE(interactions[1].getInt1() == "chr2:12000");
	EXPECT_FALSE(interactions[0].getInt1() == "chr1:17753");
	EXPECT_TRUE(interactions[0].getInt1() == "chr1:17000");
	EXPECT_TRUE(interactions[2].getInt1() == "chrX:1000");

}

/*TEST(binInterTests, constructor2)
{
	std::vector<Interaction> interactions;
	int res = 1000;
	interactions.push_back(Interaction("chrX","chr1",12553,15273));
	interactions.push_back(Interaction("chrX","chr2",17753,15273));
	interactions.push_back(Interaction("chrX","chr3",1255,1020));

	binInteractions(interactions, res);

	ASSERT_TRUE(interactions[2].getChr1() == "chrX");
	EXPECT_TRUE(interactions[1].getChr2() == "chr1");
	EXPECT_TRUE(interactions[1].getLocus2() == 15000);
	EXPECT_FALSE(interactions[0].getChr2() == "chr1");
	EXPECT_TRUE(interactions[1].getInt2() == "chr1:15000");
	EXPECT_TRUE(interactions[0].getChr1() == "chrX");

}//*/


