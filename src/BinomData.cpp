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

//mChr1(""), mChr2(""), mLocus1(0), mLocus2(0), mFrequency(0),
BinomData::BinomData():  mRelCoverage1(0), mRelCoverage2(0), mProbability(0), mExpected(0), mPvalue(0), mQvalue(0), mLogObservedOverExpected(0) {

	}

BinomData::BinomData(std::string chr1, std::string chr2, int locus1, int locus2, \
		int frequency, double relCoverage1, double relCoverage2, \
		double probability, double expected, \
		double pvalue, double qvalue, double logObservedOverExpected)\
				: mRelCoverage1(relCoverage1), mRelCoverage2(relCoverage2), \
				  mProbability(probability), mExpected(expected), \
				  mPvalue(pvalue), mQvalue(qvalue), \
				  mLogObservedOverExpected(logObservedOverExpected) {
	this->setChr1(chr1);
	this->setChr2(chr2);
	this->setLocus1(locus1);
	this->setLocus2(locus2);
	this->setFreq(frequency);
}

BinomData::BinomData(const BinomData & other){
	this->setChr1(other.getChr1());
	this->setChr2(other.getChr2());
	this->setLocus1(other.getLocus1());
	this->setLocus2(other.getLocus2());
	this->setFreq(other.getFreq());
	mRelCoverage1 = other.mRelCoverage1;
	mRelCoverage2 = other.mRelCoverage2;
	mProbability = other.mProbability;
	mExpected = other.mExpected;
	mPvalue = other.mPvalue;
	mQvalue = other.mQvalue;
	mLogObservedOverExpected = other.mLogObservedOverExpected;
}

BinomData::BinomData(const Interaction & other){
	this->setChr1(other.getChr1());
	this->setChr2(other.getChr2());
	this->setLocus1(other.getLocus1());
	this->setLocus2(other.getLocus2());
	this->setFreq(other.getFreq());
}

void BinomData::print()
{	std::ostringstream streamObj;
	string L = 	this->getChr1() + ":" + std::to_string(this->getLocus1()) + " " + this->getChr2() + ":" + std::to_string(this->getLocus2());
	streamObj << this->getFreq();
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
	if (a.getChr1() == b.getChr1())
	{
		if (a.getLocus1() == b.getLocus1())
		{
			if (a.getChr2() == b.getChr2())
			{
				return a.getLocus2() < b.getLocus2();
			}
			return a.getChr2() < b.getChr2();
		}
		return a.getLocus1() < b.getLocus1();
	}
	return a.getChr1() < b.getChr1();
}

ostream & operator<<(ostream & out, const BinomData & in)
{
	//out.precision(15);
	out << in.getChr1() << "\t" << in.getLocus1() << "\t" \
		<< in.getChr2() << "\t" << in.getLocus2() \
			<< "\t" << in.mRelCoverage1 \
			<< "\t" << in.mRelCoverage2 \
			<< "\t" << in.mProbability \
			<< "\t" << in.mExpected \
			<< "\t" << in.getFreq() \
			<< "\t" << in.mPvalue \
			<< "\t" << in.mQvalue \
			<< "\t" << in.mLogObservedOverExpected \
			<< flush;
	cerr << "print file" << endl;
	return out;
}



