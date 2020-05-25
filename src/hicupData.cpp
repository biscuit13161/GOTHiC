/*
 * hicupData.cpp
 *
 *  Created on: 5 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "hicupData.h"
#include "pbinom.h"
#include "Utils.h"
#include "padjust.h"
#include <set>
#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <cstdlib>
#include <istream>
#include <fstream>
#include <map>
#include <cmath>
#include <regex>
#include <zlib.h>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include <algorithm>

using namespace std;

struct RetrieveKey
{
    template <typename T>
    typename T::first_type operator()(T keyValuePair) const
    {
        return keyValuePair.first;
    }
};


void importHicup(string fileName, vector<Interaction> & interactions, bool checkConsistency)
{
	//fileType=ifelse(grepl("\\.bam$", fileName)|grepl("\\.sam$", fileName), "bam", "table")

	//seqan3::alignment_file_input fin_from_filename{filename};

	cerr << "Importing HiCUP file: " << fileName << endl;

	//the output of hicup is a sam file, that looks like uniques_ORIGINALFILE_trunc.sam
	//this has to be converted using the hicupToTable tool


	if (fileName.find("bam",fileName.length()-3)!=string::npos)
	{
		string str = string("samtools view -h ") + fileName + " | grep -v \"^@\" | cut -f -4 > " + fileName +".txt";
		const char *cmd = str.c_str();
		system(cmd);
		fileName = fileName + ".txt";
		//throw std::invalid_argument("importHicup: doesn't function with bam files, input can be converted to appropriate text file using hicupToTable script");
		// samtools view -h CD34-2_S1_L007_R1_2_001.hicup.bam | grep -v "^@" | cut -f -4 > CD34-2_S1_L007_R1_2_001.hicup.bam.txt
	}
	if (fileName.find("sam",fileName.length()-3)!=string::npos)
	{
		string str = string("grep -v \"^@\" ") + fileName + " | cut -f -4 > " + fileName +".txt";
		const char *cmd = str.c_str();
		system(cmd);

		fileName = fileName + ".txt";
		//throw std::invalid_argument("importHicup: doesn't function with sam files, input can be converted to appropriate text file using hicupToTable script");
	}
	if (fileName.find("gz",fileName.length()-2)!=string::npos)
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
}

void importHicupGz(string fileName, vector<Interaction> & interactions, bool checkConsistency)
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

void importHicupTxt(string fileName, vector<Interaction> & interactions, bool checkConsistency)
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
	//throw std::invalid_argument("finished");

}//*/

void mapHicupToRestrictionFragment(vector<Interaction> & interactions, vector<Site> & fragments)
{
	int iSize = interactions.size();
	cerr << endl << "Mapping HiCUP data (" << iSize << " positions) to enzyme fragments" << endl;

	string binOutFileName = "Digest.bin";
	writeBinary(fragments, binOutFileName);

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

	cerr << "Identified " << 2*sources.size() << " overlaps" << endl;

	vector<Site> ().swap(fragments);
	//sort positions into order

	sortPositions(interactions, iSize, sources, targets);

	completed();
}

void mapHicupToRestrictionFragment(vector<Interaction> & interactions, multimap<string,array<int,2>> & fragments)
{
	int iSize = interactions.size();
	cerr << endl << "Mapping HiCUP data (" << iSize << " positions) to enzyme fragments" << endl;

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

	cerr << "Identified " << 2*sources.size() << " overlaps" << endl;

	fragments.clear();
	//multimap<string,array<int,2>> ().swap(fragments);
	//sort positions into order

	sortPositions(interactions, iSize, sources, targets);

	completed();
}

void sortPositions(vector<Interaction> & interactions, int iSize, vector<halfInteraction> & sources, vector<halfInteraction> & targets)
{
	//#pragma omp parallel for
	for (int i = 0; i < iSize; i++)
	{
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
#pragma omp critical (sp1)
		//interactions.push_back(out);
		interactions[i] = out;
	}
}

void binInteractions(vector<Interaction> & interactions, int res)
{
	cerr << "Calculating Interaction Bins ("  <<interactions.size() << ")" << endl;

	vector<Interaction> interactions2;

	#pragma omp parallel for
	for (int i = 0; i < interactions.size(); i++)
	{
		// if resolution is given, bins will be calculated from interactions using the resolution
		int V = floor(interactions[i].getLocus1() / res) * res;
		if (V == 0)
		{
			V = 1;
		}
		interactions[i].setLocus1(V);
		V = floor(interactions[i].getLocus2() / res) * res;
		if (V == 0)
		{
			V = 1;
		}
		interactions[i].setLocus2(V);
	}


	// --------- count the number of interactions between bins -------
#ifdef DEBUG_FLAG
	cout << "Singles: " << interactions.size()  << endl;
#endif
	countDuplicates(interactions);

	completed();
}

void getHindIIIsitesFromHicup(vector<Site> & sites, string fileName)
{
	// ** getHindIIIsitesFromHicup using vector of Sites **

	// load sites for HindIII restriction enzyme from HiCUP_digester
	cerr << endl << "Loading Enzyme Restriction Sites" << endl;

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
	cout << "  -- Vector Variant --" << endl;

	sort(sites.begin(), sites.end(), sitecomp);

	completed();
	cerr << "Loaded " << sites.size() << " enzyme fragments" << endl;;
}//*/

void getHindIIIsitesFromHicup(multimap<string,array<int,2>> & sites, string fileName)
{
	// ** getHindIIIsitesFromHicup using multimap structure **

	// load sites for HindIII restriction enzyme from HiCUP_digester
	cerr << endl << "Loading Enzyme Restriction Sites" << endl;

	ifstream inFile;
	inFile.open(fileName);

	if (!inFile.is_open()){
		throw std::invalid_argument("getHindIIIsitesFromHicup: unable to open Restriction List file!");
	}
	int rows = 1;

	//array<int,2> pos{0,0};
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
			//pos[0] = stoi(start,nullptr,10);
			string end;
			getline(inFile,end,'\t');
			//pos[1] = stoi(end,nullptr,10);

			//sites[chr] = pos;
			//sites.emplace(chr,pos);
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
	cout << "  -- Map Variant --" << endl;

	completed();
	cerr << "Loaded " << sites.size() << " enzyme fragments" << endl;;
}//*/



void binomialHiChicup(vector<Interaction> & interactions, string sampleName, CisTrans cistrans, vector<BinomData> & binFiltered)
{
	cerr << "Binomial HiC Hicup Analysis" << endl;


	vector<Interaction> binned_df_filtered;

	bool removeDiagonal = true;

	//diagonal removal
	if (removeDiagonal)
	{
		cerr << "\tRemoving Diagonals!" << endl;
		removeDuplicates(interactions, binned_df_filtered);
		cerr << "\t";
		completed();
		interactions.clear();
	}
	else
	{
		binned_df_filtered.swap(interactions);
	}
	//cout << "size " << binned_df_filtered.size() << endl;

	if (cistrans == ct_cis)
	{
		cerr << "\tFinding Cis interactions!";
		findCis(interactions, binned_df_filtered);
		//interactions.resize(pos);
		cerr << "|\t";
		completed();
	}
	else if (cistrans == ct_trans)
	{
		cerr << "\tFinding Trans interactions!";
		findTrans(interactions, binned_df_filtered);
		//interactions.resize(pos);
		cerr << "|\t";
		completed();
	}
	else
	{
		interactions.swap(binned_df_filtered);
	}


	// all read pairs used in binomial
	//int numberOfReadPairs = interactions.size(); // before binning!!
	//cout << "Read Pairs " << numberOfReadPairs << endl;


	// calculate coverage
    map<string,int> cov;
    double tCoverage = 0;
    int numberOfReadPairs = 0;
    int max = 0;

    ofstream outFile("my_file.txt");
    for (const auto &e : interactions) outFile << e << "\n";

    calcFreq(interactions, cov, numberOfReadPairs, tCoverage, max);

    cout << "Read Pairs " << numberOfReadPairs << endl;

//#ifdef DEBUG_FLAG
    cerr << "Total Coverage: " << tCoverage  << " (57358/172074)"<< endl;
    cerr << "Max Individual Coverage: " << max << endl;
//    #endif

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
    int covS = cov.size(); // === length(all_bins)
    uint32_t numberOfAllInteractions = pow(covS,2);
    int upperhalfBinNumber = (numberOfAllInteractions - cov.size())/2;

//#ifdef DEBUG_FLAG
    cout << "Chromosomes Size: " << chromos.size() << endl;
    cout << "numberOfAllInteractions: " << numberOfAllInteractions << " (" << covS << ")" << endl;
    cout << "upperhalfBinNumber: " << upperhalfBinNumber << endl;
    cout << "diagonalProb: " << diagonalProb << endl;
//#endif

    double cisBinNumber = 0;
    double transBinNumber = 0;

    if (cistrans != ct_all)
    {
    	double sumSquare = 0;
    	map<string,int> chrlens;
    	for (auto cr : chromos){
    		set<int> pos;
    		for (auto x : interactions)
    		{
    			if (x.getChr1() == cr)
    			{
    				pos.insert(x.getLocus1());
    			}
    			if (x.getChr2() == cr)
    			{
    				pos.insert(x.getLocus2());
    			}
    		}
    		chrlens[cr] = pos.size();
    		sumSquare += pos.size()*pos.size();

    	}

    	cisBinNumber = (sumSquare - cov.size())/2;
    	transBinNumber = upperhalfBinNumber - cisBinNumber;
    }


		//diagonalProb <- sum(relative_coverage^2)
    double probabilityCorrection = 0;
    switch(cistrans)
    {
    case ct_all :
    	probabilityCorrection = (removeDiagonal)? (1/(1-diagonalProb)) : 1;
    	break;
    case ct_cis :
    	probabilityCorrection = double(upperhalfBinNumber)/cisBinNumber;
    	break;
    case ct_trans:
    	probabilityCorrection = double(upperhalfBinNumber)/transBinNumber;
    	break;
    }
//#ifdef DEBUG_FLAG
    printf("Probability Correction: %.10f \n", probabilityCorrection);
//#endif

    vector<array<double,3>> values;
    values.resize(binFiltered.size());

    // Calculate expected read counts, probabilities and pvalues
    #pragma omp parallel for
    for (int i = 0; i < binFiltered.size(); i++)
    {
    	/** Calculate expected read counts, probabilities **/
    	//double P = binFiltered[i].getRelCov1() * binFiltered[i].getRelCov2() * 2 * probabilityCorrection;
    	binFiltered[i].setProbability(binFiltered[i].getRelCov1() * binFiltered[i].getRelCov2() * 2 * probabilityCorrection);
    	binFiltered[i].setExpected(binFiltered[i].getProbability() * numberOfReadPairs);

    	/** Calculate pvalues **/
    	// /* input Freq -1 as test is greater than! */
    	double P = pbinom(double(binFiltered[i].getFreq()-2), double(numberOfReadPairs), binFiltered[i].getProbability(), false, false);
    	binFiltered[i].setPvalue(P);

    	/** observed over expected log ratio **/
    	// binned_df_filtered$logFoldChange <- log2(binned_df_filtered$frequencies/binned_df_filtered$predicted)
		binFiltered[i].setLogObExp(log2(binFiltered[i].getFreq()/double(binFiltered[i].getExpected())));

		array<double,3> ls = {double(i), binFiltered[i].getPvalue(), 0.5};
		values[i] = ls;

    }//*/



    //multiple testing correction separately for matrices with all interactions/only cis/only transs

	//pBhAdjust(double Pi, double n);
	switch(cistrans)
	{
	case ct_all :
		if(removeDiagonal)
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
		if(removeDiagonal)
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

#pragma omp parallel for
	for (int i = 0; i < binFiltered.size(); i++)
	{
		binFiltered[i].setQvalue(values[i][2]);
	}//*/

//#ifdef DEBUG_FLAG
    cerr << "BinomialData Size: " << binFiltered.size() << endl;
//#endif

    completed();
	//return binFiltered;//
}


void findOverlaps(vector<halfInteraction>& query, vector<Site> & fragments, string name)
{
	// ** findOverlaps using vector and equal_range **
	cerr << "Finding Overlaps in " << name << endl;
	int i;

	#pragma omp parallel for
	for (i = 0; i < query.size(); i++)
	{
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
	}
	completed();
}//*/


void findOverlaps(vector<halfInteraction>& query, multimap<string,array<int,2>> & fragments, string name)
{
	// ** findOverlaps for map **
	cerr << "Finding Overlaps in " << name << endl;
	int i;

	typedef multimap<string,array<int,2>>::iterator MMAPIterator;
	#pragma omp parallel for
	for (i = 0; i < query.size(); i++)
	{
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
	}
	completed();
}//*/

void countDuplicates(vector<Interaction> & interactions)
{
	cerr << "\tCounting Duplicates" << endl;

	map<string , map<string, int>> list;

	int count = 0;
#pragma omp parallel for
	for (int i = 0; i < interactions.size(); i++)
	{
		string int1 = interactions[i].getInt1();
		string int2 = interactions[i].getInt2();
#pragma omp critical (cd1)
		list[int1][int2]++;
	}
	interactions.clear();

	#pragma omp parallel for
	for (int i = 0; i< list.size(); i++)
	{
		auto it = list.begin();
		std::advance(it, i);
		string int1 = it->first;
		size_t pos = int1.find(":");
		string chr1 = int1.substr(0,pos);
		int locus1 = atoi(int1.substr(pos+1).c_str());
		map<string, int> l2 = it->second;
		for (auto ti = l2.begin(); ti != l2.end(); ti++)
		{
			int f = ti->second;

			string int2 = ti->first;

			size_t pos = int2.find(":");
			string chr2 = int2.substr(0,pos);
			int locus2 = atoi(int2.substr(pos+1).c_str());


			Interaction I = Interaction(chr1, chr2, locus1, locus2, f);
			#pragma omp critical (cd2)
			interactions.push_back(I);
		}//*/

	}

#ifdef DEBUG_FLAG
	cerr << "Duplicates: " << interactions.size()  << endl;
#endif

	cerr << "\t";
	completed();
}


void removeDuplicates(vector<Interaction> & interactions, vector<Interaction> & binned_df_filtered)
{
	int pos = 0;
	//for (auto it = interactions.begin(); it != interactions.end(); it++)
	//{
#pragma omp parallel for
	for (int i = 0; i< interactions.size(); i++)
	{
		auto it = interactions.begin();
		std::advance(it, i);//*/
		Interaction T = *it;
		if (T.getInt1() != T.getInt2())
		{
#pragma omp critical (rd1)
			binned_df_filtered.push_back(T);
			pos++;
		}
	}
}

void findCis(vector<Interaction> & interactions, vector<Interaction> & binned_df_filtered)
{
	int pos = 0;
#pragma omp parallel for
	for (int i = 0; i< binned_df_filtered.size(); i++)
	{
		auto it = binned_df_filtered.begin();
		advance(it, i);
		//	for (auto it = binned_df_filtered.begin(); it != binned_df_filtered.end(); it++)
		//	{
		Interaction T = *it;
		if (T.getChr1() == T.getChr2())
		{
#pragma omp critical (fc1)
			interactions.push_back(T);
			pos++;
		}
	}
}

void findTrans(vector<Interaction> & interactions, vector<Interaction> & binned_df_filtered)
{
	int pos = 0;
#pragma omp parallel for
	for (int i = 0; i< binned_df_filtered.size(); i++)
	{
		auto it = binned_df_filtered.begin();
		std::advance(it, i);

		Interaction T = *it;
		if (T.getChr1() != T.getChr2())
		{
#pragma omp critical (ft1)
			interactions.push_back(T);
			pos++;
		}
	}
}

void calcFreq(const vector<Interaction> & interactions, map<string,int> & cov, int & numberOfReadPairs, double & tCoverage, int & max)
{
	// WILL NOT PARALLELISE!!! ...?
	for (int i = 0; i < interactions.size(); i ++ )
	{
		/*if (cov.find(interactions[i].getInt1()) == cov.end())
		{
			//cout << "CovA " << interactions[i].getInt1() << endl;
			cov[interactions[i].getInt1()] = interactions[i].getFreq();
			tCoverage += interactions[i].getFreq();
		}
		else
		{//*/
			//cout << "CovA2 " << interactions[i].getInt1() << endl;
			cov[interactions[i].getInt1()] += interactions[i].getFreq();
			tCoverage += interactions[i].getFreq();
			numberOfReadPairs += interactions[i].getFreq();
		/*}
		if (cov.find(interactions[i].getInt2()) == cov.end())
		{
			//cout << "CovB " << interactions[i].getInt1() << endl;
			cov[interactions[i].getInt2()] = interactions[i].getFreq();
			tCoverage += interactions[i].getFreq();
		}
		else
		{//*/
			//cout << "CovB2 " << interactions[i].getInt1() << endl;
			cov[interactions[i].getInt2()] += interactions[i].getFreq();
			tCoverage += interactions[i].getFreq();
		//}
		max = (cov[interactions[i].getInt1()]> max )? cov[interactions[i].getInt1()]: max;
		max = (cov[interactions[i].getInt2()]> max )? cov[interactions[i].getInt2()]: max;
	}

}
