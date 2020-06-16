/*
 * hicupDataComp.h
 *
 *  Created on: 26 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#ifndef TESTS_HICUPDATACOMP_H_
#define TESTS_HICUPDATACOMP_H_

//#include "hicupData.h"
#include "SetupComp.h"
#include "BinomDataComp.h"
#include "tbb/concurrent_vector.h"
#include <vector>
#include <string>

class Interaction;

//void binomialHiChicupComp(std::vector<Interaction> & interactions1, std::vector<Interaction> & interactions2, SetupComp & setupValues, std::vector<BinomDataComp> & binFiltered);
void binomialHiChicupComp(tbb::concurrent_vector<Interaction> & interactions1, tbb::concurrent_vector<Interaction> & interactions2, SetupComp & setupValues, tbb::concurrent_vector<BinomDataComp> & binFiltered);
Interaction splitPair(std::string & e );
void splitPair(std::string & e, std::string & chr1, std::string & chr2,int & locus1, int & locus2);

#endif /* TESTS_HICUPDATACOMP_H_ */
