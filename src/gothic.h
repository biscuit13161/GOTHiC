/*
 * gothic.h
 *
 *  Created on: 5 May 2020
 *      Author: rich
 */

#ifndef SRC_GOTHIC_H_
#define SRC_GOTHIC_H_

#include <string>
#include <vector>
#include <math.h>
#include "Setup.h"

#include "BinomData.h"
#include "hicupData.h"



void gothicHicup(Setup & setupValues, std::vector<BinomData> & binom);
//std::vector<BinomData> gothicHicup(std::string fileName, std::string sampleName, int res, std::string restrictionFile, CisTrans cistrans, bool parallel = false);


#endif /* SRC_GOTHIC_H_ */
