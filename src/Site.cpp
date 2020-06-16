/*
 * Site.cpp
 *
 *  Created on: 18 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "Site.h"
#include <iostream>

using namespace std;

Site::Site(): mChr(""), mStart(0), mEnd(0)
{
}

Site::Site(string chr, int start, int end): mChr(chr), mStart(start), mEnd(end)
{
}

Site::Site(const Site & other)
{
	mChr = other.mChr;
	mStart = other.mStart;
	mEnd = other.mEnd;
}

const Site & Site::operator=(const Site & other)
{
	mChr = other.mChr;
	mStart = other.mStart;
	mEnd = other.mEnd;
	return *this;
} //*/

bool Site::operator==(const Site & other)
{
	return 	mChr == other.mChr  && mStart == other.mStart && mEnd == other.mEnd;
}//*/

void Site::print(){
	cout << mChr << "\t"<< mStart  << "\t"<< mEnd << endl;
}


bool sitecomp(const Site & a, const Site & b)
{
	if (a.mChr == b.mChr)
	{
		if (a.mStart == b.mStart)
		{
			return a.mEnd < b.mEnd;
		}
		return a.mStart < b.mStart;
	}
	return a.mChr < b.mChr;
}

