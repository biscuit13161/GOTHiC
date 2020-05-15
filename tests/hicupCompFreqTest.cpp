/*
 * hicupCompFreqTest.cpp
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

TEST(HiCUPTests, compFrequencies)
{
	Interaction A = Interaction("chr1","chr2",12553,15273);
	Interaction B = Interaction("chr1","chr2",12553,15273,3);

	ASSERT_FALSE(A == B);
}//*/


