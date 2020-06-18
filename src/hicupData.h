/*
 * hicupData.h
 *
 *  Created on: 5 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#ifndef SRC_HICUPDATA_H_
#define SRC_HICUPDATA_H_

#include "BinomData.h"
#include "Site.h"
#include <vector>
#include <string>
#include <math.h>
#include <streambuf>
#include "SetupData.h"


template<typename CharT, typename TraitsT = std::char_traits<CharT> >
class vectorwrapbuf : public std::basic_streambuf<CharT, TraitsT> {
public:
    vectorwrapbuf(std::vector<CharT> &vec) {
        setg(vec.data(), vec.data(), vec.data() + vec.size());
    }
};

void importHicup(std::string fileName, std::vector<Interaction> & interactions, bool checkConsistency=true);
void importHicupTxt(std::string fileName, std::vector<Interaction> & interactions, bool checkConsistency=true);
void importHicupSam(std::string fileName, std::vector<Interaction> & interactions, bool checkConsistency=true);
void importHicupGz(std::string fileName, std::vector<Interaction> & interactions, bool checkConsistency);
void mapHicupToRestrictionFragment(std::vector<Interaction> & interactions, std::vector<Site> & fragments);
void mapHicupToRestrictionFragment(std::vector<Interaction> & interactions, std::multimap<std::string,std::array<int,2>> & fragments);
void sortPositions(std::vector<Interaction> & interactions, int iSize, std::vector<halfInteraction> & sources, std::vector<halfInteraction> & targets);
void binInteractions(std::vector<Interaction> & interactions, SetupData & setupValues);
void binomialHiChicup(std::vector<Interaction> & interactions, SetupData & setupValues, std::vector<BinomData> & binFiltered);

std::string fixChromosomeNames(std::string chrnames);
void getHindIIIsitesFromHicup(std::vector<Site> & sites, std::string fileName);
void getHindIIIsitesFromHicup(std::multimap<std::string,std::array<int,2>> & sites, std::string fileName);

void findOverlaps(std::vector<halfInteraction>& query, std::vector<Site> & subject, std::string name);
void findOverlaps(std::vector<halfInteraction>& query, std::multimap<std::string,std::array<int,2>> & subject, std::string name);

void countDuplicates(std::vector<Interaction> & interactions);


void calcFreq(const std::vector<Interaction> & interactions, std::map<std::string,int> & cov, int & numberOfReadPairs, double & tCoverage, int & max);
//void completed();


#endif /* SRC_HICUPDATA_H_ */
