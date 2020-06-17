/*
 *  Interactions.cpp
 *
 *  AUTHOR
 *	Richard Thompson, ithompson@hbku.edu.qa
 *	27 May 2020
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
 

#include "Interactions.h"
#include "tbb/concurrent_vector.h"
#include "Utils.h"
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

Interaction::Interaction():mChr1(""), mChr2(""), mLocus1(0), mLocus2(0), mFrequency(1) {

}

Interaction::Interaction(string chr1, string chr2, int locus1, int locus2): mChr1(chr1), mChr2(chr2), mLocus1(locus1), mLocus2(locus2), mFrequency(1)
{

}

Interaction::Interaction(string chr1, string chr2, int locus1, int locus2, int freq): mChr1(chr1), mChr2(chr2), mLocus1(locus1), mLocus2(locus2), mFrequency(freq)
{

}

Interaction::Interaction(const Interaction & other){
	mChr1 = other.mChr1;
	mChr2 = other.mChr2;
	mLocus1 = other.mLocus1;
	mLocus2 = other.mLocus2;
	mFrequency = other.mFrequency;
}

Interaction::Interaction(const halfInteraction & first, const halfInteraction & second){
	mChr1 = first.getChr();
	mChr2 = second.getChr();
	mLocus1 = first.getLocus();
	mLocus2 = second.getLocus();
	mFrequency = 1;
}

bool operator==(const Interaction & first, const Interaction & other)
{
	return (first.mChr1 == other.mChr1) &&
			(first.mLocus1 == other.mLocus1) &&
			(first.mChr2 == other.mChr2) &&
			(first.mLocus2 == other.mLocus2);// &&
			//(mFrequency == other.getFreq());
}

bool operator!=(const Interaction & first, const Interaction & other)
		{
	return !(first == other);
		}

void checkInteractions(vector<Interaction> & interactions, vector<Interaction> & interactions2)
{
	cout << "Sizes:" << interactions.size() << " : " << interactions2.size() << endl;

	if (interactions.size() == interactions2.size())
	{
		int same= 0;
		for (int i =0 ; i < interactions.size(); i++)
		{
			if (interactions[i] == interactions2[i])
			{
				same++;
			}
			else
			{
				cout << "diff: " << i << " : " << interactions[i] << " : " << interactions2[i] << endl;
			}
		}
		cout << same << " interactions matched" << endl;
	}
}


void Interaction::print(){
	string L = mChr1 + "\t" + mChr2 +"\t" + to_string(mLocus1)+"\t" \
			+ to_string(mLocus2) + "\t" \
			+ to_string(mFrequency) + "\n";
	cout << L;
}

bool intcomp(const Interaction & a, const Interaction & b)
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

ostream & operator<<(ostream & out, const Interaction & in)
{
	string L = in.mChr1 + "\t" + in.mChr2 +"\t" + to_string(in.mLocus1)+"\t" \
				+ to_string(in.mLocus2) + "\t" \
				+ to_string(in.mFrequency);
	out << L << endl;
	return out;
}







halfInteraction::halfInteraction():mChr(""), mLocus(0) {

}

halfInteraction::halfInteraction(string chr, int locus): mChr(chr), mLocus(locus) {
}

halfInteraction::halfInteraction(const halfInteraction & other){
	mChr = other.mChr;
	mLocus = other.mLocus;
}

void halfInteraction::print(){
	cout << mChr << ":" << mLocus << endl;;
}

bool comp(const halfInteraction & a, const halfInteraction & b)
{
	if (a.mChr == b.mChr)
	{
		return a.mLocus < b.mLocus;
	}
	return a.mChr < b.mChr;
}
