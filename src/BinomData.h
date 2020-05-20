/*
 * BinomData.h
 *
 *  Created on: 5 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#ifndef SRC_BINOMDATA_H_
#define SRC_BINOMDATA_H_

#include "Utils.h"
#include <string>
#include <math.h>
#include <gtest/gtest.h>
#include <fstream>
#include <iostream>


class halfInteraction
{
private:
	std::string mChr;
	int mLocus;
	std::string mInt;


public:
	halfInteraction();
	halfInteraction(std::string chr, int locus);
	halfInteraction(const halfInteraction & other);
	friend std::string fixChromosomeNames(std::string chr);

	inline std::string getChr() const {return mChr;}
	inline int getLocus() const {return mLocus;}
	inline std::string getInt() const {return mInt;}

	void print();
	friend bool comp(const halfInteraction & a, const halfInteraction & b);

};

class Interaction
{
private:
	std::string mChr1;
	std::string mChr2;
	//chromosome(s) containing interacting regions 1 and 2

	int mLocus1;
	int mLocus2;
	//start positions of the interacting regions 1 and 2 in the corresponding chromosome(s)

	std::string mInt1;
	std::string mInt2;
	//combined strings for comparisons

	int mFrequency;

public:
	Interaction();
	Interaction(std::string chr1, std::string chr2, int locus1,	int locus2);
	Interaction(std::string chr1, std::string chr2, int locus1,	int locus2, int freq);
	Interaction(const Interaction & other);
	Interaction(const halfInteraction & first, const halfInteraction & second);
	friend std::string fixChromosomeNames(std::string chr);

	inline std::string getChr1() const {return mChr1;}
	inline int getLocus1() const {return mLocus1;}
	inline int getLocus2() const {return mLocus2;}
	inline std::string getChr2() const {return mChr2;}
	inline std::string getInt1() const {return mInt1;}
	inline std::string getInt2() const {return mInt2;}
	inline int getFreq() const {return mFrequency;}

	inline void setLocus1(int L) {mLocus1 = L; mInt1 = mChr1 + ":" + std::to_string(mLocus1);}
	inline void setLocus2(int L) {mLocus2 = L; mInt2 = mChr2 + ":" + std::to_string(mLocus2);}
	void print();

	bool operator==(const Interaction & other);
	friend std::ostream & operator<<(std::ostream & out, const Interaction & in);

};

std::ostream & operator<<(std::ostream & out, const Interaction & in);

class BinomData
{
private:
	std::string mChr1;
	std::string mChr2;
	//chromosome(s) containing interacting regions 1 and 2

	int mLocus1;
	int mLocus2;
	//start positions of the interacting regions 1 and 2 in the corresponding chromosome(s)

	std::string mInt1;
	std::string mInt2;
	//combined strings for comparisons

	int mFrequency;

	double mRelCoverage1;
	double mRelCoverage2;
	// relative coverage corresponding to regions 1 and 2

	long double mProbability; // expected frequency - probabilityOfInteraction
	long double mExpected; // expected number of reads - predicted
	int mReadCount; // observed reads number
	double mPvalue; // binomial p-value
	double mQvalue; // binomial p-value corrected for multi-testing with Benjamini-Hochberg
	double mLogObservedOverExpected; // observed/expected read numbers log ratio

public:
	BinomData();
	BinomData(std::string chr1, std::string chr2, int locus1, int locus2, \
			std::string int1, std::string int2, \
			int frequency, long double relCoverage1, long double relCoverage2, \
			long double probability, long double expected, int readCount, \
			double pvalue, double qvalue, double logObservedOverExpected);
	BinomData(const BinomData & other);
	BinomData(const Interaction & other);

	inline std::string getChr1() const {return mChr1;}
	inline int getLocus1() const {return mLocus1;}
	inline int getLocus2() const {return mLocus2;}
	inline std::string getChr2() const {return mChr2;}
	inline std::string getInt1() const {return mInt1;}
	inline std::string getInt2() const {return mInt2;}
	inline int getFreq() const {return mFrequency;}
	inline long double getRelCov1() {return mRelCoverage1;}
	inline long double getRelCov2() {return mRelCoverage2;}
	inline long double getProbability() {return mProbability;}
	inline long double getExpected() {return mExpected;} //
	inline long double getPvalue() {return mPvalue;}

	inline void setLocus1(int L) {mLocus1 = L; mInt1 = mChr1 + ":" + std::to_string(mLocus1);}
	inline void setLocus2(int L) {mLocus2 = L; mInt2 = mChr2 + ":" + std::to_string(mLocus2);}
	inline void setRelCov1(long double L) {mRelCoverage1 = L;}
	inline void setRelCov2(long double L) {mRelCoverage2 = L;}

	inline void setProbability(long double L) {mProbability = L;}
	inline void setExpected(long double L) {mExpected = L;}
	inline void setPvalue(long double L) {mPvalue = L;}

	void print();

	friend std::ostream & operator<<(std::ostream & out, const BinomData & in);
};

std::ostream & operator<<(std::ostream & out, const BinomData & in);

class Coverage
{
private:
	std::string mInt;
	int mFreq1;
	int mFreq2;

public:
	Coverage();
	inline int getFreq() const {return mFreq2;}
};

std::string fixChromosomeNames(std::string chr);
bool comp(const halfInteraction & a, const halfInteraction & b);

//double relCoverage1, double relCoverage2, double probability, int mExpected, int mReadCount, double mPvalue, double mQvalue, double mLogObservedOverExpected


#endif /* SRC_BINOMDATA_H_ */
