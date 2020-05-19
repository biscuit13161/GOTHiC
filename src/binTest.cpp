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
//#include <stdio.h>
//#include <omp.h>
//#include <istream>
//#include <fstream>
//#include <map>
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






