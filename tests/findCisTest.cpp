/*
 * findCisTest.cpp
 *
 *  Created on: 15 May 2020
 *      Author: rich
 */

#include "BinomData.h"
#include "hicupData.h"
#include <gtest/gtest.h>
#include <vector>

TEST(SystemCheck, systemCheck)
{
	std::vector<Interaction> interactions;
	interactions.push_back(Interaction("chr2","chr1",12553,15273));
	interactions.push_back(Interaction("chr1","chr1",17753,15273));
	interactions.push_back(Interaction("chrX","chr7",1255,1020));

	std::vector<Interaction> binned_df_filtered;

	findCis(binned_df_filtered, interactions);

	ASSERT_TRUE(binned_df_filtered.size() == 1);
}


