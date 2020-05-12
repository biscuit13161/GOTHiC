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

#include "../src/BinomData.h"
#include "../src/hicupData.h"




std::vector<BinomData> gothicHicup(std::string fileName, std::string sampleName, int res, std::string restrictionFile, CisTrans cistrans = ct_all, bool parallel = false, int cores = 0);


#endif /* SRC_GOTHIC_H_ */
