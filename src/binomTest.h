/*
 * binomTest.h
 *
 *  Created on: 25 May 2020
 *      Author: rich
 */

#ifndef SRC_BINOMTEST_H_
#define SRC_BINOMTEST_H_

#include <string>

double binomTest(int x, int n, double p, std::string alternative);

enum altOptions{
	ao_less = 0,
	ao_greater,
	ao_two_sided,
};

#endif /* SRC_BINOMTEST_H_ */
