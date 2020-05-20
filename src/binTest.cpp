/*
 * binTest.cpp
 *
 *  Created on: 17 May 2020
 *      Author: rich
 */

#include "binTest.h"
#include "hicupData.h"
#include "pbinom.h"
//#include <set>
#include <iostream>
#include <string>
//#include <stdio.h>
//#include <omp.h>
//#include <istream>
//#include <fstream>
#include <map>
//#include <cmath>
//#include <regex>
//#include <zlib.h>
//#include <algorithm>
//#include <boost/algorithm/string.hpp>

using namespace std;

void binTest()
{

	/*
	 * freq = vector of quantiles
	 * num = number of trials
	 * prob = probability of success on each trial.
	 * lower_tail = logical (0/1); if TRUE (default), probabilities are P[X ≤ x], otherwise, P[X > x].
	 * log_p = return p-value as log value (0/1)
	 */

	int freq = 1;
	double prob = 0.30;//6.079281e-10;
	int num =  3; //20;
	bool alt = true;
	int log_p = 0; // false
	int lower_tail = 0; //false

    prob = 6.079281e-10;
    num =  28679;

	// pbinom(freq, num, prob, alt);
	double P = pbinom(double(freq), double(num), prob, lower_tail, log_p);
	cout << "P: " << P << endl;
	//1.519785e-10//0.216
}

void binaryWriteTest(vector<Site> & fragments, string restrictionFile)
{
	cerr << "Binary Write Test" << endl;

	getHindIIIsitesFromHicup(fragments, restrictionFile);
	string binOutFileName = "Digest.bin";
	writeBinary(fragments, binOutFileName);
	completed();
}

void binaryRead(vector<Site> & fragments)
{
	cerr << "Binary Write Test" << endl;
	vector<Site> hindGR;
	string binInFileName = "Digest.bin";
	readBinary(hindGR,  binInFileName);
	//getHindIIIsitesFromHicup(hindGR, restrictionFile);


	cerr << "fragments: " << fragments.size() << endl;
	cerr << "hindGR: " << hindGR.size() << endl;

//	fragments[846223].print();
//	fragments[846224].print();
//	fragments[846225].print();
//	fragments[846226].print();
//	hindGR[846223].print();
//	hindGR[846224].print();
//	hindGR[846225].print();
//	hindGR[846226].print();

	int pos = 0;
	for (int i= 0; i!= fragments.size(); i++)
	{
		if (hindGR[i] == fragments[i])
		{
			pos++;
		}
		else
		{
			cerr << i << " ** hindGR and fragments vectors don't match **" << endl;
		}
	}

	cerr << pos << " Sites matched!" << endl;
	completed();
}

void countDupsTest()
{
	std::vector<Interaction> interactions;
	interactions.push_back(Interaction("chr1","chr2",12553,15273));
	interactions.push_back(Interaction("chr1","chr2",12553,15273));
	interactions.push_back(Interaction("chr1","chr2",12553,15273));
	interactions.push_back(Interaction("chr6","chrX",125585523,1063441));
	interactions.push_back(Interaction("chr6","chrX",125585523,1063441));
	interactions.push_back(Interaction("chr10","chr5",1064473,1505273));

	map<string,int> list;

	list["chr10:1064473"] = 1;
	list["chr6:125585523"] = 2;
	list["chr1:12553"] = 3;
	list["chr5:1505273"] = 1;
	list["chrX:1063441"] = 2;
	list["chr2:15273"] = 3;

	countDuplicates(interactions);

	for (auto i : interactions)
	{
		i.print();
		cout << list[i.getInt1()] << "\t" << list[i.getInt2()] << endl;
	}
	cout << endl;

//	for (auto i = list.begin(); i != list.end(); i++)
//		cout << i->first << " " << i->second << endl;

}

void binInterTest()
{
std::vector<Interaction> interactions;
int res = 1000;
interactions.push_back(Interaction("chr2","chr1",12553,15273));
interactions.push_back(Interaction("chr1","chr1",17753,15273));
interactions.push_back(Interaction("chrX","chr7",1255,1020));

for (auto i : interactions)
{
	i.print();
}

cout << endl;
//binInteractions(interactions, res);

std::map<std::string,int> list;
list["chr1:17753"] = 1;
list["chr2:12553"] = 1;
list["chr1:15273"] = 1;
list["chrX:1255"] = 1;
list["chr7:1020"] = 1;

std::map<std::string,int> list2;
list2["chr1:17000"] = 1;
list2["chr2:12000"] = 1;
list2["chr1:15000"] = 1;
list2["chrX:1000"] = 1;
list2["chr7:1000"] = 1;


for (auto i : interactions)
{
	i.print();
	cout << "list " << list[i.getInt1()] << " " << list[i.getInt2()] << endl;
	cout << "list " << list2[i.getInt1()] << " " << list2[i.getInt2()] << endl;
}

//ASSERT_TRUE(interactions[1].getInt1() == "chr2:12000");
//EXPECT_FALSE(interactions[0].getInt1() == "chr1:17753");
//EXPECT_TRUE(interactions[0].getInt1() == "chr1:17000");
//EXPECT_TRUE(interactions[2].getInt1() == "chrX:1000");
}



