/*
 * binTest.h
 *
 *  Created on: 17 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#ifndef SRC_BINTEST_H_
#define SRC_BINTEST_H_

#include "Site.h"
#include <vector>
#include <string>


void binTest();
void binaryWriteTest(std::vector<Site> & fragments, std::string restrictionFile);
void binaryRead(std::vector<Site> & fragments);
void countDupsTest();
void binInterTest();
void pBhAdjustTest();
void returnSizes();
void timeTest();

#endif /* SRC_BINTEST_H_ */
