/*
 * BiinomData.cpp
 *
 *  Created on: 5 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "BinomData.h"
#include <regex>
#include <iostream>
#include <boost/algorithm/string.hpp>

using namespace std;

BinomData::BinomData(): mChr1(""), mChr2(""), mLocus1(0), mLocus2(0), mRelCoverage1(0), mRelCoverage2(0), mProbability(0), mExpected(0), mReadCount(0), mPvalue(0), mQvalue(0), mLogObservedOverExpected(0) {

	}

BinomData::BinomData(std::string chr1, std::string chr2, int locus1, int locus2, double relCoverage1, double relCoverage2, double probability, int expected, int readCount, double pvalue, double qvalue, double logObservedOverExpected)
: mChr1(chr1), mChr2(chr2), mLocus1(locus1), mLocus2(locus2), mRelCoverage1(relCoverage1), mRelCoverage2(relCoverage2), mProbability(probability), mExpected(expected), mReadCount(readCount), mPvalue(pvalue), mQvalue(qvalue), mLogObservedOverExpected(logObservedOverExpected) {
}

BinomData::BinomData(BinomData & other){
	mChr1 = other.mChr1;
	mChr2 = other.mChr2;
	mLocus1 = other.mLocus1;
	mLocus2 = other.mLocus2;
	mRelCoverage1 = other.mRelCoverage1;
	mRelCoverage2 = other.mRelCoverage2;
	mProbability = other.mProbability;
	mExpected = other.mExpected;
	mReadCount = other.mReadCount;
	mPvalue = other.mPvalue;
	mQvalue = other.mQvalue;
	mLogObservedOverExpected = other.mLogObservedOverExpected;
}

Interaction::Interaction():mChr1(""), mChr2(""), mLocus1(0), mLocus2(0), mInt1(""), mInt2(""), mFrequency(0) {

}

Interaction::Interaction(string chr1, string chr2, int locus1, int locus2): mChr1(chr1), mChr2(chr2), mLocus1(locus1), mLocus2(locus2), mFrequency(0) {
		mInt1 = mChr1 + ":" + to_string(mLocus1);
		mInt2 = mChr2 + ":" + to_string(mLocus2);
}

Interaction::Interaction(string chr1, string chr2, int locus1, int locus2, int freq): mChr1(chr1), mChr2(chr2), mLocus1(locus1), mLocus2(locus2), mFrequency(freq) {
		mInt1 = mChr1 + ":" + to_string(mLocus1);
		mInt2 = mChr2 + ":" + to_string(mLocus2);
}

Interaction::Interaction(string int1, string int2, int freq): mInt1(int1), mInt2(int2), mFrequency(freq) {
	vector<string> results;
	boost::split(results, mInt1, [](char c){return c == '_';});
	mChr1 = results[0];
	mLocus1 = stoi(results[1]);
	boost::split(results, mInt2, [](char c){return c == '_';});
	mChr2 = results[0];
	mLocus2 = stoi(results[1]);
}

Interaction::Interaction(const Interaction & other): mFrequency(0){
	mChr1 = other.mChr1;
	mChr2 = other.mChr2;
	mLocus1 = other.mLocus1;
	mLocus2 = other.mLocus2;
	mInt1 = other.mInt1;
	mInt2 = other.mInt2;
}

Interaction::Interaction(const halfInteraction & first, const halfInteraction & second){
	mChr1 = first.getChr();
	mChr2 = second.getChr();
	mLocus1 = first.getLocus();
	mLocus2 = second.getLocus();
	mInt1 = first.getInt();
	mInt2 = second.getInt();
	mFrequency = 1;
}

void Interaction::print(){
	cout << mChr1 << "\t" << mChr2<< "\t" << mLocus1<< "\t" << mLocus2 << "\t" << mInt1<< "\t" << mInt2 << "\t" << mFrequency<< endl;;
}

halfInteraction::halfInteraction():mChr(""), mLocus(0), mInt("") {

}

halfInteraction::halfInteraction(string chr, int locus): mChr(chr), mLocus(locus) {
		mInt = mChr + ":" + to_string(mLocus);
}

halfInteraction::halfInteraction(const halfInteraction & other){
	mChr = other.mChr;
	mLocus = other.mLocus;
	mInt = other.mInt;
}

void halfInteraction::print(){
	cout << mChr << "\t" << mLocus << "\t" << mInt << endl;;
}

string fixChromosomeNames(string chr)
{
	//capital to small
	chr = "chr" + chr;
	chr = regex_replace(chr, regex("CHR"), "chr");
	chr = regex_replace(chr, regex("chrchr"), "chr");
	return chr;
}

bool comp(const halfInteraction & a, const halfInteraction & b)
{
	if (a.mChr == b.mChr)
	{
		return a.mLocus < b.mLocus;
	}
	return a.mChr < b.mChr;
}
