/*
 * sysTests.cpp
 *
 *  Created on: 12 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "gtest/gtest.h"

TEST(SystemCheck, systemCheck) {
	int count = 3;
	int query = 4;

	ASSERT_FALSE(count == query);
	ASSERT_EQ(count, 3);

}

