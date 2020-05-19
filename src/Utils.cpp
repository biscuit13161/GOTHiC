/*
 * Utils.cpp
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#include "Utils.h"
#include "version.h"
#include <ctime>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

using namespace std;

void showTime()
{
	// current date/time based on current system
	time_t now = time(0);

	// convert now to string form
	//char* dt = ctime(&now);
	tm *ltm = localtime(&now);

	cerr << ltm->tm_hour << ":" << flush;
	cerr << ltm->tm_min << ":" << flush;
	cerr << ltm->tm_sec << flush;
}

void completed()
{

	cerr << "\t... Completed! " << flush;
	//showTime();
	cerr << endl;
}

void printUsage()
{
	cerr << "GOTHiC++ usage:" << endl << endl;
	cerr << "    gothic <path/to/config/file>"  << endl << endl;
	cerr << string("version: ") << GOTH_MAJOR_VERSION << "." << GOTH_MINOR_VERSION << "." << GOTH_PATCH_VERSION << endl;

}

uint32_t fact(uint32_t n)
{
   if ((n==0)||(n==1))
      return 1;
   else
      return n*fact(n-1);
}

int fact(int n)
{
    //  Factorial
   if ((n==0)||(n==1))
      return 1;
   else
      return n*fact(n-1);
}

void writeBinary(vector<Site> & sites, string binOutFileName)
{
	ofstream binOutFile;
	binOutFile.open(binOutFileName, ios::binary);

	int max = 0;
	if (binOutFile.is_open())
	{
		for (int i = 0; i != sites.size(); i++)
		{
			Site P = sites[i];
			string chr = P.getChr();
			if (!chr.empty())
			{
				chr.resize(32);
				int locus = P.getLocus();
				int start = P.getStart();
				int end = P.getEnd();
				//size_t S = sizeof(std::string) + 3* sizeof(int);
				binOutFile.write(chr.c_str(), 32);// need to cast the pointer
				binOutFile.write(reinterpret_cast<char*>(&locus), sizeof(int));
				binOutFile.write(reinterpret_cast<char*>(&start), sizeof(int));
				binOutFile.write(reinterpret_cast<char*>(&end), sizeof(int));
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

void readBinary(vector<Site> & sites, string binInFileName)
{
	cerr << "Reading Binary input file" << endl;
	ifstream binInFile;
	binInFile.open(binInFileName, ios::binary);

	if (binInFile.is_open())
	{
		while(binInFile)
		{
			string chr;
			chr.resize(32);
			//binInFile.read(reinterpret_cast<char*>(& chr), 32); // need to cast the pointer
			binInFile.read(& chr[0], 32);
			chr.resize(chr.find('\0'));
			if (!chr.empty())
			{
				int locus;
				binInFile.read(reinterpret_cast<char*>(& locus), sizeof(int));
				int start;
				binInFile.read(reinterpret_cast<char*>(& start), sizeof(int));
				int end;
				binInFile.read(reinterpret_cast<char*>(& end), sizeof(int));
				Site P = Site(chr,locus,start,end);
				sites.push_back(P);
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
