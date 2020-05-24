/*
 * overlapTests.cpp
 *
 *  Created on: 12 May 2020
 *      Author: rich
 */

#include "hicupData.h"
#include "BinomData.h"
#include "Utils.h"
#include <gtest/gtest.h>
#include <vector>
#include <map>


TEST(overlapTests, findOverlapTestVec)
{
	std::vector<halfInteraction> query;
	query.push_back(halfInteraction("chr1",12345));
	query.push_back(halfInteraction("chr2",54321));

	std::vector<Site> fragments;
	fragments.push_back(Site("chr1",10000,20000));
	fragments.push_back(Site("chr1",50000,60000));
	fragments.push_back(Site("chr2",10000,20000));
	fragments.push_back(Site("chr2",50000,60000));

	findOverlaps(query,fragments, "test");

	ASSERT_TRUE(query[0].getChr() == "chr1");
	EXPECT_TRUE(query[1].getChr() == "chr2");
	EXPECT_TRUE(query[0].getLocus() == 10000);
	EXPECT_TRUE(query[1].getLocus() == 50000);
}

TEST(overlapTests, findOverlapTestMap)
{
	std::vector<halfInteraction> query;
	query.push_back(halfInteraction("chr1",12345));
	query.push_back(halfInteraction("chr2",54321));

	std::multimap<std::string,std::array<int,2>> fragments;
	std::multimap<std::string,std::array<int,2>>::iterator it = fragments.begin();
	fragments.emplace("chr1",std::array{10000,20000});
	fragments.emplace("chr1",std::array{50000,60000});
	fragments.emplace("chr2",std::array{10000,20000});
	fragments.emplace("chr2",std::array{50000,60000});

	findOverlaps(query,fragments, "test");

	ASSERT_TRUE(query[0].getChr() == "chr1");
	EXPECT_TRUE(query[1].getChr() == "chr2");
	EXPECT_TRUE(query[0].getLocus() == 10000);
	EXPECT_TRUE(query[1].getLocus() == 50000);
}


