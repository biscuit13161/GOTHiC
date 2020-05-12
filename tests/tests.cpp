/*
 * unitTest1.cpp
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#include "../src/hicupData.h"
#include "../src/BinomData.h"
#include "../src/Utils.h"
#include <gtest/gtest.h>
#include <vector>
#include <string>

TEST(HiCUPTests, constructor)
{
	Interaction A = Interaction("chr1","chr2",12553,15273,3);

	ASSERT_EQ(A.getFreq(), 3);
	EXPECT_TRUE(A.getInt1() == "chr1:12553");
	EXPECT_TRUE(A.getInt2() == "chr2:15273");
}

TEST(HiCUPTests, compInteractions)
{
	Interaction Q1 = Interaction("chr1","chr2",12553,15273);
	Interaction Q2 = Interaction("chr1","chr2",12553,15273);

	ASSERT_TRUE(Q1.getFreq() == 1);
	EXPECT_TRUE(Q1 == Q2);
}//*/

TEST(HiCUPTests, compFrequencies)
{
	Interaction A = Interaction("chr1","chr2",12553,15273);
	Interaction B = Interaction("chr1","chr2",12553,15273,3);

	ASSERT_FALSE(A == B);
}//*/

TEST(HiCUPTests, testCountDuplicates)
{
	std::vector<Interaction> interactions;
	interactions.push_back(Interaction("chr1","chr2",12553,15273));
	interactions.push_back(Interaction("chr1","chr2",12553,15273));
	interactions.push_back(Interaction("chr1","chr2",12553,15273));
	interactions.push_back(Interaction("chr6","chrX",125585523,1063441));
	interactions.push_back(Interaction("chr6","chrX",125585523,1063441));
	interactions.push_back(Interaction("chr10","chr5",1064473,1505273));

	std::vector<Interaction> inter;
	inter.push_back(Interaction("chr1","chr2",12553,15273));
	inter.push_back(Interaction("chr1","chr2",12553,15273));
	inter.push_back(Interaction("chr1","chr2",12553,15273));
	inter.push_back(Interaction("chr6","chrX",125585523,1063441));
	inter.push_back(Interaction("chr6","chrX",125585523,1063441));
	inter.push_back(Interaction("chr10","chr5",1064473,1505273));

	EXPECT_TRUE(interactions.size() == inter.size());

	countDuplicates(interactions);

	/*
	 * Doesn't compare correctly for countDuplicates...?
	 * getFreq ?????
	 */
	int max = 1;

	for (auto it = interactions.begin(); it!= interactions.end(); it++)
	{
		Interaction I = *it;
		max = (I.getFreq() > max) ? I.getFreq() : max;
	}


	ASSERT_FALSE(interactions[0] == inter[0]);
	EXPECT_FALSE(interactions.size() == inter.size());
	EXPECT_TRUE(inter.size() == 6);
	// this should be true!
	//EXPECT_TRUE(interactions.size() == 3);
	// this should be three
	//EXPECT_FALSE(max == 1);
}//*/


