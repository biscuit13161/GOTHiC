/*
 *  Site.cpp
 *
 *  AUTHOR
 *	Richard Thompson, ithompson@hbku.edu.qa
 *	May 18, 2020.
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

