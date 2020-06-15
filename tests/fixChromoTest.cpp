/*
 * fixChromoTest.cpp
 *
 *  Created on: 15 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "hicupData.h"
#include "BinomData.h"
#include "Utils.h"
#include <gtest/gtest.h>
#include <vector>


TEST(fixChromosomeNamesTest, first)
{
	ASSERT_TRUE(fixChromosomeNames("1") == "chr1");
	EXPECT_TRUE(fixChromosomeNames("X") == "chrX");
	EXPECT_TRUE(fixChromosomeNames("chr2") == "chr2");
	EXPECT_TRUE(fixChromosomeNames("CHR3") == "chr3");
}

