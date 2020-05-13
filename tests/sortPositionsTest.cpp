/*
 * sortPositionsTest.cpp
 *
 *  Created on: 13 May 2020
 *      Author: rich
 */


#include "BinomData.h"
#include "hicupData.h"
#include <gtest/gtest.h>
#include <vector>
#include <stdio.h>
#include <omp.h>


TEST(sortPositionsTest, first)
{
	omp_set_num_threads(2);

	std::vector<Interaction> interactions;
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

	ASSERT_TRUE(interactions[0].getInt2() == "chr2:12553");
	EXPECT_TRUE(interactions[1].getInt1() == "chr1:15273");
	EXPECT_TRUE(interactions[2].getInt2() == "chrX:1255");
}
