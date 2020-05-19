/*
 * Utils.h
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include "Site.h"
#include <cstdint>
#include <vector>
#include <string>

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

void writeBinary(std::vector<Site> & sites, std::string binOutFileName);
void readBinary(std::vector<Site> & sites, std::string binInFileName);

#endif /* SRC_UTILS_H_ */
