/*
 * gothic.h
 *
 *  Created on: 5 May 2020
 *      Author: rich
 */

#ifndef SOURCE_GOTHIC_H_
#define SOURCE_GOTHIC_H_

#include <string>
#include <vector>
#include <math.h>
#include "BinomData.h"
#include "hicupData.h"




std::vector<BinomData> gothicHicup(std::string fileName, std::string sampleName, int res, std::string restrictionFile, CisTrans cistrans = ct_all, bool parallel = false, int cores = 0);


#endif /* SOURCE_GOTHIC_H_ */
