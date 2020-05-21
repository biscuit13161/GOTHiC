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

public:
	halfInteraction();
	halfInteraction(std::string chr, int locus);
	halfInteraction(const halfInteraction & other);
	friend std::string fixChromosomeNames(std::string chr);

	inline std::string getChr() const {return mChr;}
	inline int getLocus() const {return mLocus;}
	inline std::string getInt() const {return mChr + ":" + std::to_string(mLocus);}

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

	int mFrequency;

public:
	Interaction();
	Interaction(std::string chr1, std::string chr2, int locus1,	int locus2);
	Interaction(std::string chr1, std::string chr2, int locus1,	int locus2, int freq);
	Interaction(const Interaction & other);
	Interaction(const halfInteraction & first, const halfInteraction & second);
	friend std::string fixChromosomeNames(std::string chr);

	inline std::string getChr1() const {return mChr1;}
	inline std::string getChr2() const {return mChr2;}
	inline int getLocus1() const {return mLocus1;}
	inline int getLocus2() const {return mLocus2;}
	inline std::string getInt1() const {return mChr1 + ":" + std::to_string(mLocus1);}
	inline std::string getInt2() const {return mChr2 + ":" + std::to_string(mLocus2);}
	inline int getFreq() const {return mFrequency;}

	inline void setLocus1(int L) {mLocus1 = L;}
	inline void setLocus2(int L) {mLocus2 = L;}
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

	int mFrequency;

	double mRelCoverage1;
	double mRelCoverage2;
	// relative coverage corresponding to regions 1 and 2

	double mProbability; // expected frequency - probabilityOfInteraction
	double mExpected; // expected number of reads - predicted
	int mReadCount; // observed reads number
	double mPvalue; // binomial p-value
	double mQvalue; // binomial p-value corrected for multi-testing with Benjamini-Hochberg
	double mLogObservedOverExpected; // observed/expected read numbers log ratio

public:
	BinomData();
	BinomData(std::string chr1, std::string chr2, int locus1, int locus2, \
			int frequency, double relCoverage1, double relCoverage2, \
			double probability, double expected, int readCount, \
			double pvalue, double qvalue, double logObservedOverExpected);
	BinomData(const BinomData & other);
	BinomData(const Interaction & other);

	inline std::string getChr1() const {return mChr1;}
	inline int getLocus1() const {return mLocus1;}
	inline int getLocus2() const {return mLocus2;}
	inline std::string getChr2() const {return mChr2;}
	inline std::string getInt1() const {return mChr1 + ":" + std::to_string(mLocus1);}
	inline std::string getInt2() const {return mChr2 + ":" + std::to_string(mLocus2);}
	inline int getFreq() const {return mFrequency;}
	inline double getRelCov1() {return mRelCoverage1;}
	inline double getRelCov2() {return mRelCoverage2;}
	inline double getProbability() {return mProbability;}
	inline double getExpected() {return mExpected;} //
	inline double getPvalue() {return mPvalue;}
	inline double getLogObExp() {return mLogObservedOverExpected;}

	inline void setLocus1(int L) {mLocus1 = L;}
	inline void setLocus2(int L) {mLocus2 = L;}

	inline void setRelCov1(double L) {mRelCoverage1 = L;}
	inline void setRelCov2(double L) {mRelCoverage2 = L;}

	inline void setProbability(double L) {mProbability = L;}
	inline void setExpected(double L) {mExpected = L;}
	inline void setPvalue(double L) {mPvalue = L;}
	inline void setLogObExp(double L) {mLogObservedOverExpected =L;}

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
