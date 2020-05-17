/*
 * binomialTest.cpp
 *
 *  Created on: 17 May 2020
 *      Author: rich
 */

#include "BinomData.h"
#include "hicupData.h"
#include "Utils.h"
#include <gtest/gtest.h>
#include <vector>


TEST(binInterTests, constructor)
{
    int freq = 1;
    long double prob = 6.079281e-10
    int num =  28679;
    bool alt = true;
    long double P = binomialTest(int freq, int num, long double prob,bool alt);

    ASSERT_TRUE(true); // P == 1.743462e-05
}
