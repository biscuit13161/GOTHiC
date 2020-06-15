/*
 * hicupTests.cpp
 *
 *  Created on: 11 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "hicupData.h"
#include "BinomData.h"
#include "Utils.h"
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

