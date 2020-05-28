/*
 * Interactions.cpp
 *
 *  Created on: 27 May 2020
 *      Author: rich
 */

#include "Interactions.h"
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

bool Interaction::operator==(const Interaction & other)
{
	return (mChr1 == other.getChr1()) &&
			(mLocus1 == other.getLocus1()) &&
			(mChr2 == other.getChr2()) &&
			(mLocus2 == other.getLocus2()) &&
			(mFrequency == other.getFreq());
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

void writeBinary(vector<Interaction> & interactions, string binOutFileName)
{
	ofstream binOutFile;
	binOutFile.open(binOutFileName, ios::binary);

	int max = 0;
	if (binOutFile.is_open())
	{
		for (int i = 0; i != interactions.size(); i++)
		{
			Interaction P = interactions[i];
			string Int1 = P.getInt1();
			string Int2 = P.getInt2();
			if (!Int1.empty())
			{
				Int1.resize(32);
				Int2.resize(32);
				int freq = P.getFreq();
				//size_t S = sizeof(std::string) + 3* sizeof(int);
				binOutFile.write(Int1.c_str(), 32);
				binOutFile.write(Int2.c_str(), 32);// need to cast the pointer
				binOutFile.write(reinterpret_cast<char*>(&freq), sizeof(int));
			}
		}

		binOutFile.close();
	}
	else
	{
		cout << "Unable to create binary write file: " + binOutFileName << endl;
	}
	cout << endl;

}

void readBinary(vector<Interaction> & interactions, string binInFileName)
{
	cerr << "Reading Binary input file" << endl;
	ifstream binInFile;
	binInFile.open(binInFileName, ios::binary);

	if (binInFile.is_open())
	{
		while(binInFile)
		{
			string Int1;
			Int1.resize(32);
			//binInFile.read(reinterpret_cast<char*>(& chr), 32); // need to cast the pointer
			binInFile.read(& Int1[0], 32);
			Int1.resize(Int1.find('\0'));
			if (!Int1.empty())
			{
				size_t pos = Int1.find(":");
				string chr1 = Int1.substr(0,pos);
				int locus1 = atoi(Int1.substr(pos+1).c_str());

				string Int2;
				Int2.resize(32);
				//binInFile.read(reinterpret_cast<char*>(& chr), 32); // need to cast the pointer
				binInFile.read(& Int2[0], 32);
				Int2.resize(Int2.find('\0'));
				pos = Int2.find(":");
				string chr2 = Int2.substr(0,pos);
				int locus2 = atoi(Int2.substr(pos+1).c_str());

				int freq;
				binInFile.read(reinterpret_cast<char*>(& freq), sizeof(int));
				Interaction P = Interaction(chr1,chr2,locus1,locus2,freq);
				interactions.push_back(P);
			}
		}
		binInFile.close();
	}
	else
	{
		cout << "Unable to read binary file: " + binInFileName << endl;
	}

	completed();
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
