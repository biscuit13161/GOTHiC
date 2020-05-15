/*
 * hicupCompInterTest.cpp
 *
 *  Created on: 15 May 2020
 *      Author: rich
 */

#include "hicupData.h"
#include "BinomData.h"
#include "Utils.h"
#include <gtest/gtest.h>
#include <vector>
#include <string>

TEST(HiCUPTests, compInteractions)
{
	Interaction Q1 = Interaction("chr1","chr2",12553,15273);
	Interaction Q2 = Interaction("chr1","chr2",12553,15273);

	ASSERT_TRUE(Q1.getFreq() == 1);
	EXPECT_TRUE(Q1 == Q2);
}//*/


