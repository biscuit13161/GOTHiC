/*
 * hicupTests.cpp
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
