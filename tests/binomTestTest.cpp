/*
 * binomTestTest.cpp
 *
 *  Created on: 26 May 2020
 *      Author: rich
 */

#include "binomTest.h"
#include "hicupData.h"
#include "BinomData.h"
#include "pbinom.h"
#include <gtest/gtest.h>
#include <vector>

TEST (binomTestTest,first)
{
	double P1 = binomTest(2, 28679, 6.079281e-10, "greater");
	double P2 = binomTest(2, 28679, 6.079281e-10, "less");

	EXPECT_NEAR(1.519785e-10, P1, 5e-15);
	EXPECT_NEAR(1, P2, 5e-15);
}

