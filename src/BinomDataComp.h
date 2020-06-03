/*
 * BinomDataComp.h
 *
 *  Created on: 27 May 2020
 *      Author: rich
 */

#ifndef SRC_BINOMDATACOMP_H_
#define SRC_BINOMDATACOMP_H_

#include "BinomData.h"

class BinomDataComp
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
	BinomDataComp();
	BinomDataComp(std::string chr1, std::string chr2, int locus1, int locus2, \
			int frequency, double relCoverage1, double relCoverage2, \
			double probability, double expected, int readCount, \
			double pvalue, double qvalue, double logObservedOverExpected);
	BinomDataComp(const BinomDataComp & other);
	BinomDataComp(const Interaction & other);

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
	inline double getQvalue() {return mQvalue;}
	inline double getLogObExp() {return mLogObservedOverExpected;}

	inline void setLocus1(int L) {mLocus1 = L;}
	inline void setLocus2(int L) {mLocus2 = L;}

	inline void setRelCov1(double L) {mRelCoverage1 = L;}
	inline void setRelCov2(double L) {mRelCoverage2 = L;}

	inline void setProbability(double L) {mProbability = L;}
	inline void setExpected(double L) {mExpected = L;}
	inline void setPvalue(double L) {mPvalue = L;}
	inline void setQvalue(double L) {mQvalue = L;}
	inline void setLogObExp(double L) {mLogObservedOverExpected = L;}

	void print();

	friend std::ostream & operator<<(std::ostream & out, const BinomDataComp & in);
	bool operator==(const BinomDataComp & other);
	friend bool bincompcomp(const BinomDataComp & a, const BinomDataComp & b);
};

std::ostream & operator<<(std::ostream & out, const BinomDataComp & in);
bool bincompcomp(const BinomDataComp & a, const BinomDataComp & b);



#endif /* SRC_BINOMDATACOMP_H_ */
