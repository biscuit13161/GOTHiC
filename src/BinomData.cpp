/*
 * BinomData.cpp
 *
 *  Created on: 5 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "BinomData.h"
#include "pbinom.h"
//#include <regex>
#include <iostream>
#include <cmath>
#include <boost/algorithm/string.hpp>

using namespace std;

BinomData::BinomData(): mChr1(""), mChr2(""), mLocus1(0), mLocus2(0), mFrequency(0), mRelCoverage1(0), mRelCoverage2(0), mProbability(0), mExpected(0), mReadCount(0), mPvalue(0), mQvalue(0), mLogObservedOverExpected(0) {

	}

BinomData::BinomData(std::string chr1, std::string chr2, int locus1, int locus2, \
		int frequency, double relCoverage1, double relCoverage2, \
		double probability, double expected, int readCount, \
		double pvalue, double qvalue, double logObservedOverExpected)\
				: mChr1(chr1), mChr2(chr2), mLocus1(locus1), mLocus2(locus2), \
				  mFrequency(frequency), \
				  mRelCoverage1(relCoverage1), mRelCoverage2(relCoverage2), \
				  mProbability(probability), mExpected(expected), \
				  mReadCount(readCount), mPvalue(pvalue), mQvalue(qvalue), \
				  mLogObservedOverExpected(logObservedOverExpected) {
}

BinomData::BinomData(const BinomData & other){
	mChr1 = other.mChr1;
	mChr2 = other.mChr2;
	mLocus1 = other.mLocus1;
	mLocus2 = other.mLocus2;
	mFrequency = other.mFrequency;
	mRelCoverage1 = other.mRelCoverage1;
	mRelCoverage2 = other.mRelCoverage2;
	mProbability = other.mProbability;
	mExpected = other.mExpected;
	mReadCount = other.mReadCount;
	mPvalue = other.mPvalue;
	mQvalue = other.mQvalue;
	mLogObservedOverExpected = other.mLogObservedOverExpected;
}

BinomData::BinomData(const Interaction & other){
	mChr1 = other.getChr1();
	mChr2 = other.getChr2();
	mLocus1 = other.getLocus1();
	mLocus2 = other.getLocus2();
	mFrequency = other.getFreq();
}

void BinomData::print()
{	std::ostringstream streamObj;
	string L = 	mChr1 + ":" + std::to_string(mLocus1) + " " + mChr2 + ":" + std::to_string(mLocus2);
	streamObj << mFrequency;
	string A = streamObj.str();
			L = L + " " + A ; \
//			+ " " + to_string(mRelCoverage1) \
//			+ " " + to_string(mRelCoverage2) \
//			+ " " + to_string(mProbability) \
//			+ " " + to_string(mExpected) \
//			+ " " + to_string(mReadCount) \
//			+ " " + to_string(mPvalue) \
//			+ "\n";
	//double mQvalue; // binomial p-value corrected for multi-testing with Benjamini-Hochberg
	//double mLogObservedOverExpected; ;
	cout << L << endl;
}

bool bincomp(const BinomData & a, const BinomData & b)
{
	if (a.mChr1 == b.mChr1)
	{
		if (a.mLocus1 == b.mLocus1)
		{
			if (a.mChr2 == b.mChr2)
			{
				return a.mLocus2 < b.mLocus2;
			}
			return a.mChr2 < b.mChr2;
		}
		return a.mLocus1 < b.mLocus1;
	}
	return a.mChr1 < b.mChr1;
}

ostream & operator<<(ostream & out, const BinomData & in)
{
	//out.precision(15);
	out << in.mChr1 << "\t" << in.mLocus1 << "\t" \
		<< in.mChr2 << "\t" << in.mLocus2 \
			<< "\t" << in.mRelCoverage1 \
			<< "\t" << in.mRelCoverage2 \
			<< "\t" << in.mProbability \
			<< "\t" << in.mExpected \
			<< "\t" << in.mFrequency \
			<< "\t" << in.mPvalue \
			<< "\t" << in.mQvalue \
			<< "\t" << in.mLogObservedOverExpected \
			<< flush;
	return out;
}



