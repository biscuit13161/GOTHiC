/*
 * Site.cpp
 *
 *  Created on: 18 May 2020
 *      Author: rich
 */

#include "Site.h"
#include <iostream>

using namespace std;

Site::Site(): mChr(""), mLocus(0), mStart(0), mEnd(0)
{
}

Site::Site(string chr, int locus, int start, int end): mChr(chr), mLocus(locus), mStart(start), mEnd(end)
{
}

Site::Site(const Site & other)
{
	mChr = other.mChr;
	mLocus = other.mLocus;
	mStart = other.mStart;
	mEnd = other.mEnd;
}

const Site & Site::operator=(const Site & other)
{
	mChr = other.mChr;
	mLocus = other.mLocus;
	mStart = other.mStart;
	mEnd = other.mEnd;
	return *this;
} //*/

bool Site::operator==(const Site & other)
{
	return 	mChr == other.mChr && mLocus == other.mLocus && mStart == other.mStart && mEnd == other.mEnd;
}//*/

void Site::print(){
	cout << mChr << "\t"<< mLocus << "\t"<< mStart  << "\t"<< mEnd << endl;
}

