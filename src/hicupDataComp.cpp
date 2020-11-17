/*
 *  hicupDataComp.cpp
 *
 *  AUTHOR
 *	Richard Thompson, ithompson@hbku.edu.qa
 *	26 May 2020.
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


#include "hicupDataComp.h"
#include "Interactions.h"
#include "UtilsComp.h"
#include "binomTest.h"
#include "padjust.h"
#include "Baits.h"
#include "random.h"
#include <algorithm>
#include <math.h> //pow
#include <stdint.h> // uint32_t
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_sort.h"

using namespace std;
using namespace tbb;


void binomialHiChicupComp(concurrent_vector<Interaction> & interactions1, concurrent_vector<Interaction> & interactions2, SetupComp & setupValues, concurrent_vector<BinomDataComp> & binFiltered)
{
	/*
	 * interaction set1 that will be used as background
	 */

	cout << "Binomial HiC Hicup Comparative Analysis" << endl;

	removeDiagonals(interactions1, setupValues.getCisTrans(), setupValues.getRemoveDiagonal(), "Control");
	removeDiagonals(interactions2, setupValues.getCisTrans(), setupValues.getRemoveDiagonal(), "Sample");

	if (setupValues.getVerbose() >= vl_info && setupValues.getRemoveDiagonal())
	{
		fprintf(stderr, "\tDiagonals Removed:\n\t    Control: %'d interactions\n", interactions1.size() );
		fprintf(stderr, "\t    Sample:  %'d interactions \n\n",interactions2.size() );
	}

	if (setupValues.getRandom()){
		cout << "Running Random subsampling" << endl;
		RandomChoose(interactions1, interactions2);

		if (setupValues.getVerbose() >= vl_info)
		{
			fprintf(stderr, "\tAfter subsetting:\n\t    Control: %'d interactions\n", interactions1.size() );
			fprintf(stderr, "\t    Sample:  %'d interactions \n\n",interactions2.size() );
		}
	}


	printf("\tGetting Pairs\n");
	concurrent_unordered_set<string> Int1; // Ints from Control
	concurrent_unordered_set<string> Int2; // Ints from Sample
	concurrent_unordered_set<string> AllInt; // all_bin
	concurrent_unordered_map<string,int> pairs1;
	concurrent_unordered_map<string,int> pairs2;


	parallel_for(
			blocked_range<concurrent_vector<Interaction>::iterator>(interactions1.begin(),	interactions1.end()),
			[&] (blocked_range<concurrent_vector<Interaction>::iterator> inter) {
		for (concurrent_vector<Interaction>::iterator it = inter.begin(); it != inter.end(); it++) {
			Interaction e = (*it);
			Int1.insert(e.getInt1());
			Int1.insert(e.getInt2());
			AllInt.insert(e.getInt1());
			AllInt.insert(e.getInt2());
			pairs1.insert(make_pair(e.getInt1() + ":" + e.getInt2(),e.getFreq()));
		}
	});//*/

	parallel_for(
			blocked_range<concurrent_vector<Interaction>::iterator>(interactions2.begin(),	interactions2.end()),
			[&] (blocked_range<concurrent_vector<Interaction>::iterator> inter) {
		for (concurrent_vector<Interaction>::iterator it = inter.begin(); it != inter.end(); it++) {
			Interaction e = (*it);
			Int2.insert(e.getInt1());
			Int2.insert(e.getInt2());
			AllInt.insert(e.getInt1());
			AllInt.insert(e.getInt2());
			pairs2.insert(make_pair(e.getInt1() + ":" + e.getInt2(),e.getFreq()));
		}
	});

	cout << "\t" << flush;
	completed();

	if (setupValues.getVerbose() >= vl_debug)
	{
		fprintf(stderr, "\tControl Interactions: %'d unique positions containing %'d pairs\n", Int1.size(),pairs1.size());
		fprintf(stderr, "\tSample Interactions:  %'d unique positions containing %'d pairs\n", Int2.size(),pairs2.size());
		fprintf(stderr, "\tAll Interactions:     %'d unique positions\n\n", AllInt.size());
	}

	concurrent_unordered_set<string> only1;
	parallel_for(pairs1.range(),
			[&] (decltype( pairs1)::range_type & inter) {
		for (concurrent_unordered_map<string,int>::iterator it = inter.begin(); it != inter.end(); it++)
		{
			string e = it->first;
			if (pairs2.find(e) == pairs2.end())
			{
				only1.insert(e);
				Interaction I = splitPair(e);
				interactions2.push_back(I);
			}
		}
	});
	cout << "\tFinished finding [mis]matching pairs in Control!" << endl;

	concurrent_unordered_set<string> only2;
	parallel_for(pairs2.range(),
			[&] (decltype( pairs2)::range_type & inter) {
		for (concurrent_unordered_map<string,int>::iterator it = inter.begin(); it != inter.end(); it++)
		{
			string e = it->first;
			if (pairs1.find(e) == pairs1.end())
			{
				only2.insert(e);
				Interaction I = splitPair(e);
				interactions1.push_back(I);
			}
		}
	});

	cout << "\tFinished finding [mis]matching pairs in Sample!" << endl;

	int numberOfReadPairs1 = 0;
	int max1 = 0;
	int min1 = 100;
	for (auto e: interactions1)
		{
		int l = e.getFreq();
		numberOfReadPairs1 += l;
		max1 = (l > max1)? l : max1;
		min1 = (l < min1)? l : min1;
		}

	int numberOfReadPairs2 = 0;
	int max2 = 0;
	int min2 = 100;
	for (auto e: interactions2)
		{
		int l = e.getFreq();
		numberOfReadPairs2 += l;
		max2 = (l > max2)? l : max2;
		min2 = (l < min2)? l : min2;
		}

	if (setupValues.getVerbose() >= vl_debug)
	{
		fprintf(stderr, "\t    %'d Read Pairs in Control\n", numberOfReadPairs1 );
		fprintf(stderr, "\t    %'d interactions unique to Control\n", only1.size() );
		//cout << "\tmax: " << max1 << "\tmin: " << min1 << endl;
		fprintf(stderr, "\t    %'d Read Pairs in Sample\n", numberOfReadPairs2 );
		//cout << "\tmax: " << max2 << "\tmin: " << min2 << endl;
		fprintf(stderr, "\t    %'d interactions unique to Sample\n\n", only2.size() );
	}

	if (interactions1.size() != interactions2.size())
		throw std::invalid_argument("binomialHiChicupComp: imbalanced interactions");

	parallel_sort(interactions1.begin(),interactions1.end(),intcomp);
	parallel_sort(interactions2.begin(),interactions2.end(),intcomp);

	if (setupValues.getVerbose() >= vl_info)
	{
		fprintf(stderr, "\tJoined data:\n\t    Control: %'d interactions\n", interactions1.size() );
		fprintf(stderr, "\t    Sample:  %'d interactions\n", interactions2.size() );
	}

	setupValues.setCovS(AllInt.size()); // === length(all.bins)
	double numberOfAllInteractions = pow(setupValues.getCovS(),2); // === length(all.bins)^2
	setupValues.setUpperhalfBinNumber( (numberOfAllInteractions - setupValues.getCovS())/2 ); // === (length(all.bins)^2-length(all.bins))/2

	if (setupValues.getVerbose() >= vl_debug)
	{
		//cout.precision(15);
		fprintf(stderr, "\n\t    Number of All Interactions:  %'.1f \n", setupValues.getCovS() );
		fprintf(stderr, "\t    Total Interactions:          %'.f \n", numberOfAllInteractions  );
		fprintf(stderr, "\t    UpperHalfBinNumber:          %'.1f \n", setupValues.getUpperhalfBinNumber() );
	}

	set<string> chromos;
	for (auto i : interactions2)
	{
		chromos.insert(i.getChr1());
	}//*/

	double sumSquare = 0;
	getSumSquare(sumSquare, chromos, interactions2);

	setupValues.setCisBinNumber( (sumSquare - setupValues.getCovS())/2 );
	setupValues.setTransBinNumber( setupValues.getUpperhalfBinNumber() - setupValues.getCisBinNumber() );

	if (setupValues.getVerbose() >= vl_debug)
	{
		//cout.precision(15);
		fprintf(stderr, "\t    Sample Sum of Squares:       %'.1f \n", sumSquare );
		fprintf(stderr, "\t    Number of Cis:               %'.1f \n", setupValues.getCisBinNumber() );
		fprintf(stderr, "\t    Number of Trans:             %'.1f \n", setupValues.getTransBinNumber() );
	}

	/** all read pairs used in binomial **/

	vector<Bait> Baits;

	if (setupValues.getBaits() != "")
	{
		readBaits(Baits, setupValues.getBaits());
		if (setupValues.getVerbose() == vl_info)
			fprintf(stderr, "\tLoaded %'d baits\n", Baits.size() );
	}//*/

	parallel_for(size_t(0),size_t(interactions2.size()),
			[&] (size_t i) {
		Interaction f = interactions1[i];
		Interaction s = interactions2[i];
		if ( f != s )
		{
			fprintf(stderr, "Mismatched Interactions:\n\t%'d \n\t%'d \n", f , s );
		}
		else if (f.getFreq() > 1 && s.getFreq() > 1)
		{
			double prob = double(f.getFreq())/numberOfReadPairs1;

			BinomDataComp I = BinomDataComp(s);
			I.setProbability(prob);
			I.setExpected(f.getFreq());

			if (!Baits.empty())
			{
				Locus L1 = Locus(I.getChr1(), I.getLocus1());
				Locus L2 = Locus(I.getChr2(), I.getLocus2());

				string str =  findOverlaps(Baits, L1, setupValues.getRes());
				I.setBaits1(str);
				//string st = str;
				str =  findOverlaps(Baits, L2, setupValues.getRes());
				I.setBaits2(str);
				//st += "\t" + str + "\n";
				//cout << st;

				if (I.getBaits1() != "" || I.getBaits2() != "")
				{
					binFiltered.push_back(I);
				}
			}//*/
			else
			{
				binFiltered.push_back(I);
			}
		}
	});

	if (setupValues.getVerbose() >= vl_info)
	{
		fprintf(stderr, "\t    Binomial Vector Size: %'d \n\n", binFiltered.size() );
	}

	setupValues.mValues.resize(binFiltered.size());

	cout << "\tcalculating P values" << endl;
	parallel_for(size_t(0),size_t(binFiltered.size()),
			[&] (size_t i) {
		string str = "";
		int F = binFiltered[i].getFreq();
		double V = binFiltered[i].getProbability();
		double P = binomTest(F, numberOfReadPairs2, V, "two.sided", str);
		binFiltered[i].setPvalue(P);

		double Fd = log2((double(binFiltered[i].getFreq())/numberOfReadPairs2)/(double(binFiltered[i].getExpected())/numberOfReadPairs1));
		binFiltered[i].setLogObExp(Fd);

		array<double,3> ls = {double(i), P, 0.5};
		setupValues.mValues[i] = ls;
	});

	cout << "\t" << flush;
	completed();

	if (setupValues.getQvalue() == qv_ihw)
	{
		string cmd = string ("Rscript --vanilla -e 'if (!require(IHW)) {quit(status = 11)} ' ");

		int sys = system(cmd.c_str());

		if (setupValues.getVerbose() >= vl_debug)
			cout << "\tsys value: " << sys << endl;

		if (sys != 0)
		{
			cerr << "R library IHW not found!\n\tPvalue correction switched to Benjamini Hochberg" << endl;
			setupValues.setQvalue("bh");
		}
	}


	if (setupValues.getQvalue() == qv_bh)
	{
		BHCalculation(binFiltered, setupValues);
	}

	completed();
}

Interaction splitPair(string & e)
{
	size_t pos = e.find(":");
	string chr1 = e.substr(0,pos);

	e = e.substr(pos+1);
	pos = e.find(":");
	int locus1 = atoi(e.substr(0,pos).c_str());

	e = e.substr(pos+1);
	pos = e.find(":");
	string chr2 = e.substr(0,pos);
	int locus2 = atoi(e.substr(pos+1).c_str());

	Interaction I = Interaction(chr1,chr2,locus1,locus2,1);
	return I;
}

void BHCalculation(concurrent_vector<BinomDataComp> & binFiltered, SetupComp & setupValues)
{
	cout << "\tcalculating Q values using Benjamini Hochberg" << endl;
	BHCorrection(binFiltered, setupValues);
	cout << "\t" << flush;
	completed();
}

void BHCorrection(concurrent_vector<BinomDataComp> & binFiltered, SetupComp & setupValues)
{
	switch(setupValues.getCisTrans())
			{
			case ct_all :
				if(setupValues.getRemoveDiagonal())
				{
					pBhAdjust(setupValues.mValues, setupValues.getUpperhalfBinNumber());
					//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(), upperhalfBinNumber));
				}
				else
				{
					pBhAdjust(setupValues.mValues, setupValues.getUpperhalfBinNumber()+setupValues.getCovS());
					//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(), upperhalfBinNumber+covS));
				}
				break;
			case ct_cis :
				if(setupValues.getRemoveDiagonal())
				{
					pBhAdjust(setupValues.mValues, setupValues.getCisBinNumber());
					//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(), cisBinNumber));
				}
				else
				{
					pBhAdjust(setupValues.mValues, setupValues.getCisBinNumber()+setupValues.getCovS());
					//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(),cisBinNumber+covS));
				}
				break;
			case ct_trans:
				pBhAdjust(setupValues.mValues, setupValues.getTransBinNumber());
				//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(), transBinNumber));
				break;
			}

			parallel_for(size_t(0),size_t(binFiltered.size()),
					[&] (size_t i) {

				binFiltered[i].setQvalue(setupValues.mValues[i][2]);
			});
}
