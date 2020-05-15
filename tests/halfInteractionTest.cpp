/*
 * halfInteractionTest.cpp
 *
 *  Created on: 15 May 2020
 *      Author: rich
 */

#include "hicupData.h"
#include "BinomData.h"
#include "Utils.h"
#include <gtest/gtest.h>
#include <vector>


TEST(halfInteractionTest, constructor)
{
	halfInteraction I = halfInteraction("chr1", 1357908642);
	halfInteraction T(I);

	ASSERT_TRUE(I.getInt() == "chr1:1357908642");
	ASSERT_TRUE(T.getInt() == "chr1:1357908642");
}


