/*
 * Utils.h
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include <cstdint>

enum CisTrans
{
	ct_all=0,
	ct_cis,
	ct_trans,
};

void showTime();
void completed();
void printUsage();

std::uint32_t fact(std::uint32_t n);
int fact(int n);

#endif /* SRC_UTILS_H_ */
