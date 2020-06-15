/*
 * sumSquareCompTest.cpp
 *
 *  Created on: 1 Jun 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "UtilsComp.h"
#include "Interactions.h"
#include <set>
#include <string>
#include "gtest/gtest.h"
#include "tbb/concurrent_vector.h"

TEST(sortPositionsComp Test, first)
{
	tbb::concurrent_vector<Interaction> interactions;
	Interaction I = Interaction("chr2","chr1",12553,15273,1);
	interactions.push_back(I);
	I = Interaction("chr1","chr1",17753,150273,1);
	interactions.push_back(I);
	I = Interaction("chrX","chr7",1255,1020,1);
	interactions.push_back(I);
	I = Interaction("chr1","chr1",17753,17753,1);
	interactions.push_back(I);
	I = Interaction("chr1","chr1",17753,15273,1);
	interactions.push_back(I);
	I = Interaction("chrX","chr7",12550,1020,1);
	interactions.push_back(I);
	I = Interaction("chr21","chrX",1255,1255,1);
	interactions.push_back(I);

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

