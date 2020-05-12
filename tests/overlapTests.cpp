/*
 * overlapTests.cpp
 *
 *  Created on: 12 May 2020
 *      Author: rich
 */

#include "../src/hicupData.h"
#include "../src/BinomData.h"
#include "../src/Utils.h"
#include <gtest/gtest.h>
#include <vector>
//#include <string>

TEST(overlapTests, constructor)
{
	halfInteraction I = halfInteraction("chr1", 1357908642);
	halfInteraction T(I);

	ASSERT_TRUE(I.getInt() == "chr1:1357908642");
	ASSERT_TRUE(T.getInt() == "chr1:1357908642");
}

TEST(overlapTests, findOverlapTest)
{
	std::vector<halfInteraction> query;
	query.push_back(halfInteraction("chr1",12345));
	query.push_back(halfInteraction("chr2",54321));

	std::vector<Site> fragments;
	fragments.push_back(Site("chr1",10000,10000,20000));
	fragments.push_back(Site("chr1",50000,50000,60000));
	fragments.push_back(Site("chr2",10000,10000,20000));
	fragments.push_back(Site("chr2",50000,50000,60000));

	findOverlaps(query,fragments, "test");

	ASSERT_TRUE(query[0].getChr() == "chr1");
	EXPECT_TRUE(query[1].getChr() == "chr2");
	EXPECT_TRUE(query[0].getLocus() == 10000);
	EXPECT_TRUE(query[1].getLocus() == 50000);
}

TEST(overlapTests, fixChromosomeNamesTest)
{
	ASSERT_TRUE(fixChromosomeNames("1") == "chr1");
	EXPECT_TRUE(fixChromosomeNames("X") == "chrX");
	EXPECT_TRUE(fixChromosomeNames("chr2") == "chr2");
	EXPECT_TRUE(fixChromosomeNames("CHR3") == "chr3");
}
