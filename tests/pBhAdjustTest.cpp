/*
 * pBhAdjustTest.cpp
 *
 *  Created on: 21 May 2020
 *      Author: rich
 */

#include "padjust.h"
#include <gtest/gtest.h>

TEST(pBhAdjust, first)
{
	double P = 6.079281e-10;
	double n = 28679;
	double o = pBhAdjust(P, n);

	ASSERT_NEAR(1.743477e-05, o, 5e-14);
}


