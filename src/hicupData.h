/*
 * hicupData.h
 *
 *  Created on: 5 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#ifndef SRC_HICUPDATA_H_
#define SRC_HICUPDATA_H_

#include <vector>
#include <string>
#include <math.h>
#include <streambuf>

#include "../src/BinomData.h"

class Site{
private:
	std::string mChr;
	int mLocus;
	int mStart;
	int mEnd;
public:
	Site();
	Site(std::string chr, int locus, int start, int end);
	Site(const Site & other);

	friend std::string fixChromosomeNames(std::string chr);

	inline std::string getChr() const {return mChr;}
	inline int getLocus() const {return mLocus;}
	inline int getEnd() const {return mEnd;}
	inline int getStart() const {return mStart;}

	void print();
};

template<typename CharT, typename TraitsT = std::char_traits<CharT> >
class vectorwrapbuf : public std::basic_streambuf<CharT, TraitsT> {
public:
    vectorwrapbuf(std::vector<CharT> &vec) {
        setg(vec.data(), vec.data(), vec.data() + vec.size());
    }
};

void importHicup(std::string fileName, std::vector<Interaction> & interactions, bool checkConsistency=true);
void importHicupTxt(std::string fileName, std::vector<Interaction> & interactions, bool checkConsistency=true);
void importHicupGz(std::string fileName, std::vector<Interaction> & interactions, bool checkConsistency);
void mapHicupToRestrictionFragment(std::vector<Interaction> & interactions, std::string restrictionFile);
void sortPositions(std::vector<Interaction> & interactions, int iSize, std::vector<halfInteraction> & sources, std::vector<halfInteraction> & targets);
void binInteractions(std::vector<Interaction> & interactions, int res);
std::vector<BinomData> binomialHiChicup(std::vector<Interaction> & interactions, std::string restrictionFile, std::string sampleName, CisTrans cistrans = ct_all, bool parallel = false, int cores = 1, bool removeDiagonal = true);

std::string fixChromosomeNames(std::string chrnames);
void getHindIIIsitesFromHicup(std::vector<Site> & sites, std::string fileName);

void findOverlaps(std::vector<halfInteraction>& query, std::vector<Site> & subject, std::string name, bool drop = false);
void countDuplicates(std::vector<Interaction> & interactions);

void completed();

#endif /* SRC_HICUPDATA_H_ */
