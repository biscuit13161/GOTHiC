/*
 * sumSquareTest.cpp
 *
 *  Created on: 28 May 2020
 *      Author: rich
 */

#include "Interactions.h"
#include <set>
#include <string>
#include "gtest/gtest.h"

TEST(sortPositionsTest, first)
{
	std::vector<Interaction> interactions;
	interactions.push_back(Interaction("chr2","chr1",12553,15273));
	interactions.push_back(Interaction("chr1","chr1",17753,150273));
	interactions.push_back(Interaction("chrX","chr7",1255,1020));
	interactions.push_back(Interaction("chr1","chr1",17753,17753));
	interactions.push_back(Interaction("chr1","chr1",17753,15273));
	interactions.push_back(Interaction("chrX","chr7",12550,1020));
	interactions.push_back(Interaction("chr21","chrX",1255,1255));

	std::set<std::string> chromos;
	for (auto i : interactions)
	{
		chromos.insert(i.getChr1());
	}
	double sumSquare = 0;

	/** count unique interactions **/
	getSumSquare(sumSquare, chromos, interactions);


	ASSERT_FALSE(sumSquare == 0);
	EXPECT_TRUE(sumSquare == 15);
}
