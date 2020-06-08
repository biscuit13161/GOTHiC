/*
 * BinomDataComp.cpp
 *
 *  Created on: 27 May 2020
 *      Author: rich
 */

#include "BinomDataComp.h"

using namespace std;

BinomDataComp::BinomDataComp():  mProbability(0), mExpected(0), mReadCount(0), mPvalue(0), mQvalue(0), mLogObservedOverExpected(0) {

	}

BinomDataComp::BinomDataComp(std::string chr1, std::string chr2, int locus1, int locus2, \
		int frequency, double probability, double expected, int readCount, \
		double pvalue, double qvalue, double logObservedOverExpected)\
				: mProbability(probability), mExpected(expected), \
				  mReadCount(readCount), mPvalue(pvalue), mQvalue(qvalue), \
				  mLogObservedOverExpected(logObservedOverExpected) {
	this->setChr1(chr1);
	this->setChr2(chr2);
	this->setLocus1(locus1);
	this->setLocus2(locus2);
	this->setFreq(frequency);
}

BinomDataComp::BinomDataComp(const BinomDataComp & other){
	this->setChr1(other.getChr1());
	this->setChr2(other.getChr2());
	this->setLocus1(other.getLocus1());
	this->setLocus2(other.getLocus2());
	this->setFreq(other.getFreq());
	mProbability = other.mProbability;
	mExpected = other.mExpected;
	mReadCount = other.mReadCount;
	mPvalue = other.mPvalue;
	mQvalue = other.mQvalue;
	mLogObservedOverExpected = other.mLogObservedOverExpected;
}

BinomDataComp::BinomDataComp(const Interaction & other): mProbability(0), mExpected(0), mReadCount(0), mPvalue(0), mQvalue(0), mLogObservedOverExpected(0){
	this->setChr1(other.getChr1());
	this->setChr2(other.getChr2());
	this->setLocus1(other.getLocus1());
	this->setLocus2(other.getLocus2());
	this->setFreq(other.getFreq());
}

void BinomDataComp::print()
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

bool bincompcomp(const BinomDataComp & a, const BinomDataComp & b)
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

ostream & operator<<(ostream & out, const BinomDataComp & in)
{
	//out.precision(15);
	out << in.getChr1() << "\t" << in.getLocus1() << "\t" \
		<< in.getChr2() << "\t" << in.getLocus2() \
			<< "\t" << in.mProbability \
			<< "\t" << in.mExpected \
			<< "\t" << in.getFreq()
			<< "\t" << in.mPvalue \
			<< "\t" << in.mQvalue \
			<< "\t" << in.mLogObservedOverExpected \
			<< flush;
	return out;
}

bool BinomDataComp::operator==(const BinomDataComp & other)
{
	return (this->getChr1() == other.getChr1()) &&
			(this->getLocus1() == other.getLocus1()) &&
			(this->getChr2() == other.getChr2()) &&
			(this->getLocus2() == other.getLocus2());//*/
}
