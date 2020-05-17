/*
 * hicupData.h
 *
 *  Created on: 5 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#ifndef SRC_HICUPDATA_H_
#define SRC_HICUPDATA_H_

#include "BinomData.h"
#include <vector>
#include <string>
#include <math.h>
#include <streambuf>

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
void mapHicupToRestrictionFragment(std::vector<Interaction> & interactions, std::vector<Site> & fragments);
void sortPositions(std::vector<Interaction> & interactions, int iSize, std::vector<halfInteraction> & sources, std::vector<halfInteraction> & targets);
void binInteractions(std::vector<Interaction> & interactions, int res);
std::vector<BinomData> binomialHiChicup(std::vector<Interaction> & interactions, std::vector<Site> & fragments, std::string sampleName, CisTrans cistrans, bool parallel = false, bool removeDiagonal = true);

std::string fixChromosomeNames(std::string chrnames);
void getHindIIIsitesFromHicup(std::vector<Site> & sites, std::string fileName);

void findOverlaps(std::vector<halfInteraction>& query, std::vector<Site> & subject, std::string name, bool drop = false);
void countDuplicates(std::vector<Interaction> & interactions);
void removeDuplicates(std::vector<Interaction> & interactions, std::vector<Interaction> & binned_df_filtered);

void findTrans(std::vector<Interaction> & interactions, std::vector<Interaction> & binned_df_filtered);
void findCis(std::vector<Interaction> & interactions, std::vector<Interaction> & binned_df_filtered);

void calcFreq(std::vector<Interaction> & interactions, std::map<std::string,int> & cov, double & tCoverage, int & max);
//void completed();

#endif /* SRC_HICUPDATA_H_ */
