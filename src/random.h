/*
 * random.h
 *
 *  Created on: 22 Oct 2020
 *      Author: rich
 */

#ifndef SRC_RANDOM_H_
#define SRC_RANDOM_H_

#include "Interactions.h"
#include "BinomData.h"
#include "tbb/concurrent_vector.h"

void RandomSubset(tbb::concurrent_vector<Interaction> & large, tbb::concurrent_vector<Interaction> & small);
void RandomChoose(tbb::concurrent_vector<Interaction> & first, tbb::concurrent_vector<Interaction> & second);


#endif /* SRC_RANDOM_H_ */
