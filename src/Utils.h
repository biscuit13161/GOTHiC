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
#include "stdlib.h"
#include "sys/types.h"
#include "sys/sysinfo.h"

enum CisTrans
{
	ct_all=0,
	ct_cis,
	ct_trans,
};

void showTime();
void completed(int n = 0);
void printUsage();
void verbosePrint(std::string & str, bool verbose = true);

void verbose(const char * fmt, ... );

std::uint32_t fact(std::uint32_t n);
int fact(int n);

void writeBinary(std::vector<Site> & sites, std::string binOutFileName);
void readBinary(std::vector<Site> & sites, std::string binInFileName);

//get memory usage
int getRealValue();
int getVirtValue();
int parseLine(char* line);

#endif /* SRC_UTILS_H_ */
