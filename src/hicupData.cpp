/*
 *  hicupData.cpp
 *
 *  AUTHOR
 *	Richard Thompson, ithompson@hbku.edu.qa
 *	5 May 2020
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

#include "hicupData.h"
#include "pbinom.h"
#include "dbinom.h"
#include "Utils.h"
#include "padjust.h"
#include "binomTest.h"
#include <set>
#include <syslog.h>
#include <iostream>
#include <stdio.h>
#include <string>
//#include <omp.h>
#include <cstdlib>
#include <istream>
#include <fstream>
#include <map>
#include <cmath>
#include <regex>
#include <zlib.h>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "tbb/parallel_for.h"
#include "tbb/concurrent_unordered_map.h"

#include <algorithm>

using namespace std;
using namespace tbb;

struct RetrieveKey
{
    template <typename T>
    typename T::first_type operator()(T keyValuePair) const
    {
        return keyValuePair.first;
    }
};


void importHicup(string fileName, concurrent_vector<Interaction> & interactions, bool checkConsistency)
{
	cerr << "Importing HiCUP file: " << fileName << endl;

	//the output of hicup is a sam file, that looks like uniques_ORIGINALFILE_trunc.sam
	//this has to be converted using the hicupToTable tool


	if (fileName.find("bam",fileName.length()-3)!=string::npos)
	{
		string f2 = fileName + ".txt";
		if (file_exists (f2))
		{
			cerr << "\tbam.txt file already exists" << endl;
		}
		else
		{
			cerr << "\tconverting Bam file" << endl;
			string str = string("samtools view -h --no-PG ") + fileName + " | grep -v \"^@\" | cut -f -4 > " + f2;
			const char *cmd = str.c_str();
			system(cmd);
			//throw std::invalid_argument("importHicup: doesn't function with bam files, input can be converted to appropriate text file using hicupToTable script");
			cerr << "\t";
			completed();
		}
		fileName = f2;
	}
	if (fileName.find("sam",fileName.length()-3)!=string::npos)
	{
		importHicupSam(fileName, interactions, checkConsistency);
	}
	else if (fileName.find("gz",fileName.length()-2)!=string::npos)
	{
		//importHicupTxt(fileName, interactions, checkConsistency);

		throw std::invalid_argument("importHicup: doesn't function with compressed files, please decompress");
	}
	else
	{
		importHicupTxt(fileName, interactions, checkConsistency);
	}

	sort(interactions.begin(),interactions.end(), intcomp);
	completed();
	fprintf(stderr, "\tLoaded %'d positions\n\n", interactions.size());
}

void importHicupGz(string fileName, concurrent_vector<Interaction> & interactions, bool checkConsistency)
{
	/*	gzFile inFile = gzopen(fileName.c_str(), "rb");

		if (inFile == NULL) {
			throw std::invalid_argument("Failed to open compressed input file " + fileName);
		}

		unsigned char unzipBuffer[8192];
		unsigned int unzippedBytes;
		vector<unsigned char> unzippedData;
		while (true) {
		    unzippedBytes = gzread(inFile, unzipBuffer, 8192);
		    if (unzippedBytes > 0) {
		        unzippedData.insert(unzippedData.end(), unzipBuffer, unzipBuffer + unzippedBytes);
		    } else {
		        break;
		    }
		}
		gzclose(inFile);

		//istream in;
/*
		for (auto it = unzippedData.begin(); it != unzippedData.end(); it++)
		{
			unsigned char T = *it;
			in << T << flush;
		}

		//membuf sbuf(unzippedData.begin(), unzippedData.end());
		//membuf sbuf(&unzippedData, &unzippedData + sizeof(unzippedData));
		//istream in(&unzippedData);

		vectorwrapbuf<unsigned char> databuf(unzippedData);
		std::istream in(&databuf);

		while(in)
		{
			string id1;
			getline(in,id1,'\t');
			if (id1.length() == 0)
			{
				break;
			}
			string flag1;
			getline(in,flag1,'\t');
			string chr1;
			getline(in,chr1,'\t');
			string start1;
			getline(in,start1,'\n');
			int locus1 = stoi(start1,nullptr,10);

			string id2;
			getline(in,id2,'\t');
			if (checkConsistency && (id2.length() == 0 | id1 != id2))
			{
				throw std::invalid_argument("importHicup: reads must be paired in consecutive rows!");
			}
			string flag2;
			getline(in,flag1,'\t');
			string chr2;
			getline(in,chr2,'\t');
			string start2;
			getline(in,start2,'\n');
			int locus2 = stoi(start2,nullptr,10);

			cout << chr1 << "\t" << chr2 << endl;

			Interaction a = Interaction(chr1, chr2, locus1, locus2);

			interactions.push_back(a);
		}
//*/
}

void importHicupSam(string fileName, concurrent_vector<Interaction> & interactions, bool checkConsistency)
{
	ifstream inFile;
	inFile.open(fileName);

	if (!inFile.is_open())
	{
		throw std::invalid_argument("importHicup: unable to open input file!");
	}

	char sub = '@';
	string split = "\t";

	int i = 1;
	while(!inFile.bad())
	{
		//string str = string("grep -v \"^@\" ") + fileName + " | cut -f -4 > " + fileName +".txt";
		string str = "";
		string waste = "";
		getline(inFile,str,'\n');

		//cout  << i << "\t";
		if (str.length() > 0 && str.find(sub) != 0)
		{
			string id1 = str.substr(0, str.find(split));
			str.erase(0, str.find(split) + split.length());
			//cout   << id1 << "|\t";

			string flag1 = str.substr(0,str.find(split));
			str.erase(0, str.find(split) + split.length());

			string chr1 = str.substr(0,str.find(split));
			//cout << chr1 << "|\t";
			str.erase(0, str.find(split) + split.length());

			string start1 =str.substr(0,str.find(split));
			str.erase(0, str.find(split) + split.length());
			//cout << start1 << "\t";

			int locus1 = stoi(start1,nullptr,10);


			getline(inFile,str,'\n');

			string id2 = str.substr(0,str.find(split));
			str.erase(0, str.find(split) + split.length());
			if (checkConsistency && (id2.length() == 0 | id1 != id2))
			{
				throw std::invalid_argument("importHicup: reads must be paired in consecutive rows!");
			}
			string flag2 = str.substr(0,str.find(split));
			str.erase(0, str.find(split) + split.length());

			string chr2 = str.substr(0,str.find(split));
			str.erase(0, str.find(split) + split.length());
			//cout << chr2 << "\t";

			string start2 = str.substr(0,str.find(split));
			str.erase(0, str.find(split) + split.length());
			//cout << "|" << start2 << "|";
			int locus2  = stoi(start2,nullptr,10);

			Interaction a = Interaction(chr1, chr2, locus1, locus2);
			//a.print();

			interactions.push_back(a);
		}

		i++;
		//cout  << endl;

		if (inFile.eof()) // corrects for last blank line
		{
			break;
		}
	}

}//*/


void importHicupTxt(string fileName, concurrent_vector<Interaction> & interactions, bool checkConsistency)
{
	ifstream inFile;
	inFile.open(fileName);

	if (!inFile.is_open())
	{
		throw std::invalid_argument("importHicup: unable to open input file!");
	}

	while(inFile)
	{
		string id1;
		getline(inFile,id1,'\t');
		if (id1.length() == 0)
		{
			break;
		}
		string flag1 = "";
		getline(inFile,flag1,'\t');
		string chr1 = "";
		getline(inFile,chr1,'\t');
		string start1;
		getline(inFile,start1,'\n');
		int locus1 = stoi(start1,nullptr,10);

		string id2 = "";
		getline(inFile,id2,'\t');
		if (checkConsistency && (id2.length() == 0 | id1 != id2))
		{
			throw std::invalid_argument("importHicup: reads must be paired in consecutive rows!");
		}
		string flag2 = "";
		getline(inFile,flag1,'\t');
		string chr2 = "";
		getline(inFile,chr2,'\t');
		string start2 = "";
		getline(inFile,start2,'\n');
		int locus2 = stoi(start2,nullptr,10);

		Interaction a = Interaction(chr1, chr2, locus1, locus2);

		interactions.push_back(a);

		if (!inFile) // corrects for last blank line
		{
			break;
		}
	}

}//*/

void mapHicupToRestrictionFragment(concurrent_vector<Interaction> & interactions, vector<Site> & fragments)
{
	int iSize = interactions.size();
	fprintf(stderr,"Mapping HiCUP data to enzyme fragments\n");
	//cerr << endl << "Mapping HiCUP data (" << iSize << " positions) to enzyme fragments" << endl;

	//string binOutFileName = "Digest.bin";
	//writeBinary(fragments, binOutFileName);

	vector<halfInteraction> sources;
	vector<halfInteraction> targets;

	for (int i = 0; i < iSize; i++)
	{
		sources.push_back(halfInteraction(interactions[i].getChr1(), interactions[i].getLocus1() ));
		targets.push_back(halfInteraction(interactions[i].getChr2(), interactions[i].getLocus2() ));
	}

	interactions.clear();
	interactions.resize(iSize);

	findOverlaps(sources, fragments, "Sources");
	findOverlaps(targets, fragments, "Targets");

	fprintf(stderr, "\tIdentified %'d overlaps\n",2*sources.size() );

	vector<Site> ().swap(fragments);
	//sort positions into order

	sortPositions(interactions, iSize, sources, targets);

	completed();
}

void mapHicupToRestrictionFragment(concurrent_vector<Interaction> & interactions, multimap<string,array<int,2>> & fragments)
{
	int iSize = interactions.size();
	fprintf(stderr,"Mapping HiCUP data to enzyme fragments\n");
	//cerr << endl << "Mapping HiCUP data (" << iSize << " positions) to enzyme fragments" << endl;

	vector<halfInteraction> sources;
	vector<halfInteraction> targets;

	for (int i = 0; i < iSize; i++)
	{
		sources.push_back(halfInteraction(interactions[i].getChr1(), interactions[i].getLocus1() ));
		targets.push_back(halfInteraction(interactions[i].getChr2(), interactions[i].getLocus2() ));
	}

	interactions.clear();
	interactions.resize(iSize);

	findOverlaps(sources, fragments, "Sources");
	findOverlaps(targets, fragments, "Targets");

	fprintf(stderr, "Identified %'d overlaps\n",2*sources.size() );

	fragments.clear();

	sortPositions(interactions, iSize, sources, targets);

	completed();
}

void sortPositions(concurrent_vector<Interaction> & interactions, int iSize, vector<halfInteraction> & sources, vector<halfInteraction> & targets)
{
	//#pragma omp parallel for
	//for (int i = 0; i < iSize; i++)
	//{
	parallel_for(size_t(0),size_t(interactions.size()),
						[&] (size_t i) {
		halfInteraction first = sources[i];
		halfInteraction second = targets[i];
		Interaction out;

		if (first.getChr() == second.getChr())
		{
			if(first.getLocus() <= second.getLocus())
			{
				out = Interaction(first,second);
			}
			else
			{
				out = Interaction(second, first);
			}
		}
		else if (first.getChr() <= second.getChr())
		{
			out = Interaction(first,second);
		}
		else
		{
			out = Interaction(second, first);
		}

		//interactions.push_back(out);
		interactions[i] = out;
	});
}

void binInteractions(concurrent_vector<Interaction> & interactions, SetupData & setupValues)
{
	fprintf(stderr,"Calculating Interaction Bins (%'d)\n",interactions.size());
	vector<Interaction> interactions2;

	//#pragma omp parallel for
	//for (int i = 0; i < interactions.size(); i++)
	parallel_for(size_t(0),size_t(interactions.size()),
					[&] (size_t i) {
	//{
		// if resolution is given, bins will be calculated from interactions using the resolution
		int V = floor(interactions[i].getLocus1() / setupValues.getRes()) * setupValues.getRes();
		if (V == 0)
		{
			V = 1;
		}
		interactions[i].setLocus1(V);
		V = floor(interactions[i].getLocus2() / setupValues.getRes()) * setupValues.getRes();
		if (V == 0)
		{
			V = 1;
		}
		interactions[i].setLocus2(V);
	});


	// --------- count the number of interactions between bins -------
	if (setupValues.getVerbose())
		fprintf(stderr,"\tSingles: %'d (must match above count)\n", interactions.size());

	countDuplicates(interactions);

	completed();
}

void getHindIIIsitesFromHicup(vector<Site> & sites, string fileName, SetupData & setupValues)
{
	// ** getHindIIIsitesFromHicup using vector of Sites **

	// load sites for HindIII restriction enzyme from HiCUP_digester
	fprintf(stderr,"Loading Enzyme Restriction Sites\n");

	ifstream inFile;
	inFile.open(fileName);

	if (!inFile.is_open()){
		throw std::invalid_argument("getHindIIIsitesFromHicup: unable to open Restriction List file!");
	}
	int rows = 1;

	while (inFile)
	{
		if (rows > 2)
		{
			string chr;
			getline(inFile,chr,'\t');
			if (chr.length() == 0)
			{
				break;
			}
			string start1;
			getline(inFile,start1,'\t');
			int start = stoi(start1,nullptr,10);
			string end;
			getline(inFile,end,'\t');

			Site site = Site(fixChromosomeNames(chr), start, stoi(end,nullptr,10));

			sites.push_back(site);
			if (!inFile) // corrects for last blank line
			{
				break;
			}
		}
		string waste;
		getline(inFile,waste);
		rows++;
	}
	inFile.close();
	if (setupValues.getVerbose())
		cerr << "  -- Vector Variant --" << endl;

	sort(sites.begin(), sites.end(), sitecomp);

	completed();
	fprintf(stderr,"\tLoaded %'d enzyme fragments\n\n", sites.size());

}//*/

void getHindIIIsitesFromHicup(multimap<string,array<int,2>> & sites, string fileName, SetupData & setupValues)
{
	// ** getHindIIIsitesFromHicup using multimap structure **

	// load sites for HindIII restriction enzyme from HiCUP_digester
	fprintf(stderr,"Loading Enzyme Restriction Sites\n");

	ifstream inFile;
	inFile.open(fileName);

	if (!inFile.is_open()){
		throw std::invalid_argument("getHindIIIsitesFromHicup: unable to open Restriction List file!");
	}
	int rows = 1;
	string chr;
	while (inFile)
	{

		multimap<string,array<int,2>>::iterator it = sites.begin();
		if (rows > 2)
		{
			getline(inFile,chr,'\t');
			if (chr.length() == 0)
			{
				break;
			}
			string start;
			getline(inFile,start,'\t');

			string end;
			getline(inFile,end,'\t');

			sites.insert(it++, pair(chr,array{stoi(start,nullptr,10),stoi(end,nullptr,10)}));
			if (!inFile) // corrects for last blank line
			{
				break;
			}
		}
		string waste;
		getline(inFile,waste);
		rows++;
	}

	inFile.close();
	if (setupValues.getVerbose())
		cerr << "  -- Map Variant --" << endl;

	completed();
	fprintf(stderr, "\tLoaded %'d enzyme fragments\n\n", sites.size() );
}//*/



void binomialHiChicup(concurrent_vector<Interaction> & interactions, SetupData & setupValues, vector<BinomData> & binFiltered)
{
	cerr << "Binomial HiC Hicup Analysis" << endl;

	if (setupValues.getVerbose())
		fprintf(stderr, "\tInteractions size: %'d\n", interactions.size());

	removeDiagonals(interactions, setupValues.getCisTrans(), setupValues.getRemoveDiagonal());

	if (setupValues.getVerbose())
		fprintf(stderr, "\tRemoved size: %'d\n", interactions.size());


	// calculate coverage
    map<string,int> cov;
    double tCoverage = 0;
    int numberOfReadPairs = 0;
    int max = 0;

    //ofstream outFile("my_file.txt");
    //for (const auto &e : interactions) outFile << e << "\n";

    calcFreq(interactions, cov, numberOfReadPairs, tCoverage, max);



    if (setupValues.getVerbose())
    {
    	fprintf(stderr, "\t%'d Read Pairs \n", numberOfReadPairs);
    	fprintf(stderr, "\tTotal Coverage: %'.f\n", tCoverage );
    	fprintf(stderr, "\tMax Individual Coverage: %'d\n", max );
    }

    if (tCoverage == 0)
    {
    	throw std::invalid_argument("binomialHiChicup: Zero coverage!");
    }
    if (cov.size() == 0)
    {
    	throw std::invalid_argument("binomialHiChicup: Zero coverage!");
    }

    // Calculate Relative Coverage
    map<string,double> rCov;
    double diagonalProb = 0;
    for (auto i = cov.begin(); i != cov.end(); i ++)
    {
    	double V = i->second/tCoverage;
    	rCov[i->first] = V;
    	diagonalProb += V*V;
    	//relative_coverage <- coverage/sumcov
    }

    set<string> chromos;


    // add relative coverage to BinomData objects
    for (auto i : interactions) //.begin(); i < interactions.end(); i++)
    {
    	chromos.insert(i.getChr1());
    	string int1 = i.getInt1();
    	string int2 = i.getInt2();
    	BinomData bin = BinomData(i);
    	bin.setRelCov1(rCov[int1]);
    	bin.setRelCov2(rCov[int2]);
    	binFiltered.push_back(bin);
    }

    //probability correction assuming on average equal probabilities for all interactions
    float covS = cov.size(); // === length(all_bins)
    double numberOfAllInteractions = pow(covS,2);
    double upperhalfBinNumber = (numberOfAllInteractions - cov.size())/2;

    if (setupValues.getVerbose())
    {
    	fprintf(stderr, "\tChromosomes Size:         %'d\n", chromos.size() );
    	fprintf(stderr, "\tnumberOfAllInteractions:  %'.f\n", numberOfAllInteractions);
    	fprintf(stderr, "\tupperhalfBinNumber:       %'.1f\n", upperhalfBinNumber);
    	fprintf(stderr, "\tdiagonalProb:             %e\n", diagonalProb);
    }


    double cisBinNumber = 0;
    double transBinNumber = 0;

    if (setupValues.getCisTrans() != ct_all)
    {
    	double sumSquare = 0;
    	getSumSquare(sumSquare, chromos, interactions);

    	cisBinNumber = (sumSquare - cov.size())/2;
    	transBinNumber = upperhalfBinNumber - cisBinNumber;
    }

    if (setupValues.getVerbose())
    {
    	fprintf(stderr, "\tcisBinNumber:             %'.1f\n", cisBinNumber);
    	fprintf(stderr, "\ttransBinNumber:           %'.1f\n", transBinNumber);
    }

    double probabilityCorrection = 0;
    switch(setupValues.getCisTrans())
    {
    case ct_all :
    	probabilityCorrection = (setupValues.getRemoveDiagonal())? (1/(1-diagonalProb)) : 1;
    	break;
    case ct_cis :
    	probabilityCorrection = double(upperhalfBinNumber)/cisBinNumber;
    	break;
    case ct_trans:
    	probabilityCorrection = double(upperhalfBinNumber)/transBinNumber;
    	break;
    }

    if (setupValues.getVerbose())
    	fprintf(stderr, "\tProbability Correction:   %e \n\n", probabilityCorrection);


    vector<array<double,3>> values;
    values.resize(binFiltered.size());

    /* Calculate expected read counts, probabilities and pvalues */
    //#pragma omp parallel for
    //for (int i = 0; i < binFiltered.size(); i++)
    //{
    parallel_for(size_t(0),size_t(binFiltered.size()),
    					[&] (size_t i) {
    	/** Calculate expected read counts, probabilities **/
    	binFiltered[i].setProbability(binFiltered[i].getRelCov1() * binFiltered[i].getRelCov2() * 2 * probabilityCorrection);
    	binFiltered[i].setExpected(binFiltered[i].getProbability() * numberOfReadPairs);

    	/** Calculate pvalues **/
    	double P = binomTest(binFiltered[i].getFreq(), numberOfReadPairs, binFiltered[i].getProbability(), "greater");
    	binFiltered[i].setPvalue(P);

    	/** observed over expected log ratio **/
    	binFiltered[i].setLogObExp(log2(binFiltered[i].getFreq()/double(binFiltered[i].getExpected())));

		array<double,3> ls = {double(i), binFiltered[i].getPvalue(), 0.5};
		values[i] = ls;

    });//*/



    //multiple testing correction separately for matrices with all interactions/only cis/only transs

	//pBhAdjust(double Pi, double n);
	switch(setupValues.getCisTrans())
	{
	case ct_all :
		if(setupValues.getRemoveDiagonal())
		{
			pBhAdjust(values, upperhalfBinNumber);
			//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(), upperhalfBinNumber));
		}
		else
		{
			pBhAdjust(values, upperhalfBinNumber+covS);
			//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(), upperhalfBinNumber+covS));
		}
		break;
	case ct_cis :
		if(setupValues.getRemoveDiagonal())
		{
			pBhAdjust(values, cisBinNumber);
			//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(), cisBinNumber));
		}
		else
		{
			pBhAdjust(values, cisBinNumber+covS);
			//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(),cisBinNumber+covS));
		}
		break;
	case ct_trans:
		pBhAdjust(values, transBinNumber);
		//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(), transBinNumber));
		break;
	}//*/

//#pragma omp parallel for
	//for (int i = 0; i < binFiltered.size(); i++)
	//{
    parallel_for(size_t(0),size_t(binFiltered.size()),
    					[&] (size_t i) {
		binFiltered[i].setQvalue(values[i][2]);
	});//*/

//#ifdef DEBUG_FLAG
    fprintf(stderr, "\tBinomialData Size: %'d\n", binFiltered.size() );
//#endif

    completed();
	//return binFiltered;//
}


void findOverlaps(vector<halfInteraction>& query, vector<Site> & fragments, string name)
{
	// ** findOverlaps using vector and equal_range **
	cout << "\tFinding Overlaps in " << name << endl;
	int i;

	//#pragma omp parallel for
	//for (i = 0; i < query.size(); i++)
	//{
	    parallel_for(size_t(0),size_t(query.size()),
	    					[&] (size_t i) {
		auto result = equal_range(fragments.begin(),fragments.end(),query[i].getChr(), Comp{});
		for ( auto it = result.first; it != result.second; ++it )
		{
			if (query[i].getChr() == (*it).getChr())
			{
				if ((query[i].getLocus() >= (*it).getStart()) && ((*it).getEnd() >= query[i].getLocus()) )
				{
					query[i] = halfInteraction((*it).getChr(), (*it).getStart());
					break;
				}

			}
		}
	});//*/

	cerr << "\t";
	completed();
}//*/

// x = search value
// l = first (left) element index (0)
// r = last (right) element index (size() - 1)
// arr[] = array of elements
/*int m = 0;

while (r - l > 1)
{
	m =  l + (r-1) / 2;
	if (arr[m] <= x)
		l = m;
	else
		r = m;
}
if (arr[l] == x)
	return l;
else
	return -1;//*/


void findOverlaps(vector<halfInteraction>& query, multimap<string,array<int,2>> & fragments, string name)
{
	// ** findOverlaps for map **
	cout << "\nFinding Overlaps in " << name << endl;

	int i;

	typedef multimap<string,array<int,2>>::iterator MMAPIterator;
	//#pragma omp parallel for
	//for (i = 0; i < query.size(); i++)
	//{
		parallel_for(size_t(0),size_t(query.size()),
			    					[&] (size_t i) {
		pair<MMAPIterator, MMAPIterator> result = fragments.equal_range(query[i].getChr());
		for (auto it = result.first; it != result.second; it++)
		{
			array<int,2> list = it->second;
				if ((query[i].getLocus() >= list[0]) && (list[1] >= query[i].getLocus()) )
				{
					query[i] = halfInteraction(query[i].getChr(), list[0]);
					break;
				}
		}
	});
	cerr << "\t";
	completed();
}//*/

void countDuplicates(concurrent_vector<Interaction> & interactions)
{
	cerr << "\tCounting Duplicates" << endl;

	concurrent_unordered_map<string , concurrent_unordered_map<string, int>> list;

	int count = 0;
//#pragma omp parallel for
	//for (int i = 0; i < interactions.size(); i++)
	//{
	parallel_for(size_t(0),size_t(interactions.size()),
						[&] (size_t i) {
		string int1 = interactions[i].getInt1();
		string int2 = interactions[i].getInt2();
//#pragma omp critical (cd1)
		list[int1][int2]++;
	});
	interactions.clear();

	//#pragma omp parallel for
	//for (int i = 0; i< list.size(); i++)
	//{
	parallel_for(size_t(0),size_t(list.size()),
						[&] (size_t i) {
		auto it = list.begin();
		std::advance(it, i);
		string int1 = it->first;
		size_t pos = int1.find(":");
		string chr1 = int1.substr(0,pos);
		int locus1 = atoi(int1.substr(pos+1).c_str());
		concurrent_unordered_map<string, int> l2 = it->second;
		for (auto ti = l2.begin(); ti != l2.end(); ti++)
		{
			int f = ti->second;

			string int2 = ti->first;

			size_t pos = int2.find(":");
			string chr2 = int2.substr(0,pos);
			int locus2 = atoi(int2.substr(pos+1).c_str());


			Interaction I = Interaction(chr1, chr2, locus1, locus2, f);
			//#pragma omp critical (cd2)
			interactions.push_back(I);
		}//*/

	});


	syslog(LOG_NOTICE, "Duplicates: %d", interactions.size());


	cerr << "\t";
	completed();
}


void calcFreq(const concurrent_vector<Interaction> & interactions, map<string,int> & cov, int & numberOfReadPairs, double & tCoverage, int & max)
{
	// WILL NOT PARALLELISE!!! ...?
	for (int i = 0; i < interactions.size(); i ++ )
	{
		cov[interactions[i].getInt1()] += interactions[i].getFreq();
		tCoverage += interactions[i].getFreq();

		numberOfReadPairs += interactions[i].getFreq();
		cov[interactions[i].getInt2()] += interactions[i].getFreq();
		tCoverage += interactions[i].getFreq();

		max = (cov[interactions[i].getInt1()]> max )? cov[interactions[i].getInt1()]: max;
		max = (cov[interactions[i].getInt2()]> max )? cov[interactions[i].getInt2()]: max;
	}
}


