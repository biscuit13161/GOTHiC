/*
 * BinomDataComp.h
 *
 *  Created on: 27 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#ifndef SRC_BINOMDATACOMP_H_
#define SRC_BINOMDATACOMP_H_

#include "BinomData.h"

class BinomDataComp: public BinomData
{
private:
	std::string mBaits1;
	std::string mBaits2;

	double mProbability; // expected frequency - probabilityOfInteraction
	double mExpected; // expected number of reads - predicted
	int mReadCount; // observed reads number
	double mPvalue; // binomial p-value
	double mQvalue; // binomial p-value corrected for multi-testing with Benjamini-Hochberg
	double mLogObservedOverExpected; // observed/expected read numbers log ratio

public:
	BinomDataComp();
	BinomDataComp(std::string chr1, std::string chr2, int locus1, int locus2, \
			int frequency, double probability, double expected, int readCount, \
			double pvalue, double qvalue, double logObservedOverExpected);
	BinomDataComp(const BinomDataComp & other);
	BinomDataComp(const Interaction & other);

	inline std::string getBaits1() const {return mBaits1;}
	inline std::string getBaits2() const {return mBaits2;}

	inline double getProbability() {return mProbability;}
	inline double getExpected() {return mExpected;} //
	inline double getPvalue() {return mPvalue;}
	inline double getQvalue() {return mQvalue;}
	inline double getLogObExp() {return mLogObservedOverExpected;}

	inline void setBaits1(std::string L) {mBaits1 = L;}
	inline void setBaits2(std::string L) {mBaits2 = L;}

	inline void setProbability(double L) {mProbability = L;}
	inline void setExpected(double L) {mExpected = L;}
	inline void setPvalue(double L) {mPvalue = L;}
	inline void setQvalue(double L) {mQvalue = L;}
	inline void setLogObExp(double L) {mLogObservedOverExpected = L;}

	void print();//*/

	friend std::ostream & operator<<(std::ostream & out, const BinomDataComp & in);//*/
	bool operator==(const BinomDataComp & other);
	friend bool bincompcomp(const BinomDataComp & a, const BinomDataComp & b);
};

std::ostream & operator<<(std::ostream & out, const BinomDataComp & in);
bool bincompcomp(const BinomDataComp & a, const BinomDataComp & b);



#endif /* SRC_BINOMDATACOMP_H_ */
