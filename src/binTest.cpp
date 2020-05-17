/*
 * binTest.cpp
 *
 *  Created on: 17 May 2020
 *      Author: rich
 */

#include "../src/hicupData.h"

#include <set>
#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <istream>
#include <fstream>
#include <map>
#include <cmath>
#include <regex>
#include <zlib.h>
#include <algorithm>
#include <boost/algorithm/string.hpp>

using namespace std;

void binTest()
{
	int freq = 1;
	long double prob = 0.30;//6.079281e-10;
	int num =  3; //20;
	bool alt = true;

	cout << fact(num) << endl;
	cout << (fact(num)/(fact(freq)*fact(num-freq))) << endl;
	cout << binomialCoefficients( num,freq) << endl;


	long double P = binomialTest(freq, num, prob, alt);
	cout << "P: " << P << endl;
}


