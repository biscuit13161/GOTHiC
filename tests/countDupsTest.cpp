/*
 * countDupsTests.cpp
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#include "hicupData.h"
#include "BinomData.h"
#include "Utils.h"
#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <map>

TEST(countDupsTests, testCountDuplicates)
{
	std::vector<Interaction> interactions;
	interactions.push_back(Interaction("chr1","chr2",12553,15273));
	interactions.push_back(Interaction("chr1","chr2",12553,15273));
	interactions.push_back(Interaction("chr1","chr2",12553,15273));
	interactions.push_back(Interaction("chr6","chrX",125585523,1063441));
	interactions.push_back(Interaction("chr6","chrX",125585523,1063441));
	interactions.push_back(Interaction("chr10","chr5",1064473,1505273));

	std::vector<Interaction> inter;
	inter.push_back(Interaction("chr1","chr2",12553,15273));
	inter.push_back(Interaction("chr1","chr2",12553,15273));
	inter.push_back(Interaction("chr1","chr2",12553,15273));
	inter.push_back(Interaction("chr6","chrX",125585523,1063441));
	inter.push_back(Interaction("chr6","chrX",125585523,1063441));
	inter.push_back(Interaction("chr10","chr5",1064473,1505273));

	std::map<std::string,int> list;
	list["chr10:1064473"] = 1;
	list["chr6:125585523"] = 2;
	list["chr1:12553"] = 3;
	list["chr5:1505273"] = 1;
	list["chrX:1063441"] = 2;
	list["chr2:15273"] = 3;


	countDuplicates(interactions);

	ASSERT_TRUE(interactions.size() == 3);

	for (auto i : interactions)
	{
		EXPECT_TRUE( i.getFreq() == list[i.getInt1()]);
		EXPECT_TRUE( i.getFreq() == list[i.getInt2()]);
	}

}//*/


