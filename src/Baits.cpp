/*
 *  Baits.cpp
 *
 *  AUTHOR
 *	Richard Thompson, ithompson@hbku.edu.qa
 *	June 15, 2020.
 *
 * 	Copyright (C) 2020 Richard Thompson, Qatar Biomedical Research Institute
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
 
#include "Baits.h"
#include "UtilsComp.h"
#include <iostream>
#include <fstream>

using namespace std;

Bait::Bait(): mChr(""), mStart(0), mEnd(0), mBait("")
	{

	}

Bait::Bait(std::string Chr, int Start, int End, std::string Bait): mChr(Chr), mStart(Start), mEnd(End), mBait(Bait)
{

}

Bait::Bait(const Bait & other )
{
	mChr = other.mChr;
	mStart = other.mStart;
	mEnd = other.mEnd;
	mBait = other.mBait;
}

Locus::Locus(): mChr(""), mStart(0)
{

}

Locus::Locus(std::string Chr, int Start): mChr(Chr), mStart(Start)
{

}

Locus::Locus(const Locus & other)
{
	mChr = other.mChr;
	mStart = other.mStart;
}

void readBaits(std::vector<Bait> & Baits, std::string fileName)
{
	// ** readBaits using vector of Baits **

	cerr << "\tLoading Baits" << endl;

	ifstream inFile;
	inFile.open(fileName);

	if (!inFile.is_open()){
		throw std::invalid_argument("readBaits: unable to open Baits file!");
	}

	while (inFile)
	{
		string chr;
		getline(inFile,chr,'\t');
		if (chr.length() == 0)
		{
			break;
		}
		string start1;
		getline(inFile,start1,'\t');
		string end;
		getline(inFile,end,'\t');
		string strand;
		getline(inFile,strand,'\t');
		string Baitin;
		getline(inFile,Baitin);

		Bait site = Bait(fixChromosomeNames(chr), stoi(start1,nullptr,10), stoi(end,nullptr,10), Baitin);

		Baits.push_back(site);
		if (!inFile) // corrects for last blank line
		{
			break;
		}
	}
	inFile.close();

	sort(Baits.begin(), Baits.end(), comp);

	cerr << "\t" << flush;
	completed();
}

bool comp(const Bait & a, const Bait & b)
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

ostream & operator<<(ostream & out, const Bait & in)
{
	//out.precision(15);
	out << in.mChr << "\t" << in.mStart << "\t" \
		<< in.mEnd << "\t" << in.mBait \
			<< flush;
	return out;
}//*/


string findOverlaps(vector<Bait> & Baits, Locus & L, int res)
{
	auto result = equal_range(Baits.begin(),Baits.end(),L.getChr(), Comp{});
	for ( auto it = result.first; it != result.second; ++it )
	{
		if ((L.getEnd(res)+1 >= (*it).getStart()) && ((*it).getEnd() >= L.getStart()) )
		{
/*			string str = L.getChr() + "\t" + to_string(L.getStart()) + "\t" + to_string((*it).getStart()) + "\n";
			str += "\t" + to_string(L.getEnd(res)) + "\t" + to_string((*it).getEnd()) + "\n";
			str += (*it).getBait() + "\n";
			cout << str ;*/
			return (*it).getBait();
		}
	}
	return "";
}

