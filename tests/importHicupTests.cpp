/*
 * importHicupTests.cpp
 *
 *  Created on: 12 May 2020
 *      Author: rich
 */

#include "../src/hicupData.h"
#include "../src/BinomData.h"
#include "../src/Utils.h"
#include <gtest/gtest.h>
#include <vector>
#include <fstream>

TEST(importHicupTests, constructor)
{
	std::vector<Interaction> interactions;

	importHicupTxt("../../example_data/example_Test_Data.txt", interactions);

	/*std::ifstream inFile;
	inFile.open("../../example_data/example_Test_Data.txt");
	if (!inFile.is_open())
	{
		throw std::invalid_argument("importHicup: unable to open input file!");
	}//*/

	ASSERT_TRUE(interactions.size() == 2);
	ASSERT_TRUE(interactions[0].getInt1() == "chr17_gl000204_random:51677");
	ASSERT_TRUE(interactions[0].getInt2() == "chr2:155678032");
	ASSERT_TRUE(interactions[1].getInt1() == "chr3:101663949");
	ASSERT_TRUE(interactions[1].getInt2() == "chr2:230303040");
}


