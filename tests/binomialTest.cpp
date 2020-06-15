/*
 * binomialTest.cpp
 *
 *  Created on: 17 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "pbinom.h"
#include <gtest/gtest.h>
#include <vector>


TEST(binInterTests, constructor)
{
	/*
	 * freq = vector of quantiles
	 * num = number of trials
	 * prob = probability of success on each trial.
	 * lower_tail = logical (0/1); if TRUE (default), probabilities are P[X ≤ x], otherwise, P[X > x].
	 * log_p = return p-value as log value (0/1)
	 */

	double freq = 1;
	double prob = 0.30;//6.079281e-10;
	double num =  3; //20;
	bool alt = true;
	int log_p = 0; // false
	int lower_tail = 0; //false

	double P = pbinom(freq, num, prob, lower_tail, log_p);
	double P2 = pbinom(freq, num, prob, 1, log_p);
	double P3 = pbinom(freq, num, prob, 0, 1);
	double P4 = pbinom(1, 28679, 6.079281e-10, 0, 0);

    ASSERT_DOUBLE_EQ(0.216, P); // P == 1.743462e-05
    EXPECT_DOUBLE_EQ(0.784, P2);
    EXPECT_NEAR(-1.53248, P3, 0.00005); //-1.532477
    EXPECT_NEAR(1.519785e-10, P4, 5e-15); //1.51979e-10
}
