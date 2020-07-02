/*
 *  BinomDataComp.cpp
 *
 *  AUTHOR
 *	Richard Thompson, ithompson@hbku.edu.qa
 *	May 27, 2020.
 *
 *	Copyright (C) 2020 Richard Thompson, Qatar Biomedical Research Institute
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
 *
 */
 

#include "BinomDataComp.h"

using namespace std;

BinomDataComp::BinomDataComp():  mProbability(0), mExpected(0), mReadCount(0),
		mPvalue(0), mQvalue(0), mLogObservedOverExpected(0),
		mBaits1(""), mBaits2(""){

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
	mBaits1 = other.mBaits1;
	mBaits2 = other.mBaits2;
}

BinomDataComp::BinomDataComp(const Interaction & other): mProbability(0), mExpected(0), mReadCount(0), mPvalue(0), mQvalue(0), mLogObservedOverExpected(0){
	this->setChr1(other.getChr1());
	this->setChr2(other.getChr2());
	this->setLocus1(other.getLocus1());
	this->setLocus2(other.getLocus2());
	this->setFreq(other.getFreq());
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
}//*/

bool BinomDataComp::operator==(const BinomDataComp & other)
{
	return (this->getChr1() == other.getChr1()) &&
			(this->getLocus1() == other.getLocus1()) &&
			(this->getChr2() == other.getChr2()) &&
			(this->getLocus2() == other.getLocus2());//*/
}
