/*
 * hicupDataComp.h
 *
 *  Created on: 26 May 2020
 *      Author: rich
 */

#ifndef TESTS_HICUPDATACOMP_H_
#define TESTS_HICUPDATACOMP_H_

#include "hicupData.h"
#include "SetupComp.h"
#include "BinomDataComp.h"
#include <vector>
#include <string>

void binomialHiChicupComp(std::vector<Interaction> & interactions1, std::vector<Interaction> & interactions2, SetupComp & setupValues, std::vector<BinomDataComp> & binFiltered);
Interaction splitPair(std::string & e );

#endif /* TESTS_HICUPDATACOMP_H_ */
