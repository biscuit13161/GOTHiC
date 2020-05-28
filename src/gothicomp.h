/*
 * gothicomp.h
 *
 *  Created on: 26 May 2020
 *      Author: rich
 */

#ifndef SRC_GOTHICOMP_H_
#define SRC_GOTHICOMP_H_

#include "SetupComp.h"
#include "BinomDataComp.h"
#include "tbb/concurrent_vector.h"
#include <vector>

void gothicHicupComp(SetupComp & setupValues, tbb::concurrent_vector<BinomDataComp> & binom);

#endif /* SRC_GOTHICOMP_H_ */
