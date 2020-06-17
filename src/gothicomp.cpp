/*
 *  gothicomp.cpp
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


#include "gothicomp.h"
#include "Interactions.h"
#include "hicupDataComp.h"
#include "binTest.h"
#include <iostream>
#include <stdio.h>
#include "tbb/concurrent_vector.h"
#include "tbb/task_scheduler_init.h"

using namespace std;
using namespace tbb;


int main(int argc, char *argv[])
{
	if ( argc != 2 ) // argc should be 2 for correct execution
	{
		// We print argv[0] assuming it is the program name
		cerr<<"usage: "<< argv[0] <<" <filename>\n";
		printUsageComp();
		return 0;
	}

	vector<string> allArgs(argv, argv + argc);
	SetupComp setupValues = loadConfigComp(allArgs[1]);

	setupValues.print();

	checkInputFiles(setupValues.getCondition1());
	checkInputFiles(setupValues.getCondition2());

	tbb::task_scheduler_init init(setupValues.getThreads());

	concurrent_vector<BinomDataComp> binom;

	try {
		gothicHicupComp(setupValues,binom);
		//sumSquareTest();

		sort(binom.begin(), binom.end(), bincompcomp);
	}
	catch(const std::invalid_argument& e){
		cerr << "Error: " << e.what() << endl;
	}

	string fileName = setupValues.getSname()+".binom.txt";
	ofstream binomFile(fileName);
	binomFile << "chr1" << "\t" << "locus1" \
			<< "\t" << "chr2" << "\t" << "locus2" \
			<< "\t" << "probability" \
			<< "\t" << "expected" \
			<< "\t" << "readCount" \
			<< "\t" << "pvalue" \
			<< "\t" << "qvalue" \
			<< "\t" << "logObservedOverExpected" << endl;
	for (const auto &e : binom) binomFile << e << endl;

	return 0;
}

void gothicHicupComp(SetupComp & setupValues, concurrent_vector<BinomDataComp> & binom)
{
	/** load data from GOTHiC++ **/
	concurrent_vector<Interaction> interactions1;
	readBinary(interactions1, setupValues.getCondition1());

	concurrent_vector<Interaction> interactions2;
	readBinary(interactions2, setupValues.getCondition2());

	if (setupValues.getVerbose())
	{
		cerr << "\tControl: " << interactions1.size() << " interactions" <<endl;
		cerr << "\tSample:  " << interactions2.size() << " interactions" <<endl;
	}

	binomialHiChicupComp(interactions1, interactions2, setupValues, binom);

}

