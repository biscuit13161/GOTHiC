/*
 * IHW.cpp
 *
 *  Created on: 28 Jun 2020
 *      Author: rich
 */

#include "IHW.h"
#include <algorithm>
#include <iostream>
//#include <limits>
#include "tbb/parallel_for.h"
#include "tbb/parallel_sort.h"
#include "tbb/blocked_range.h"

using namespace std;
using namespace tbb;

//We use the (expected+observed)/2 number of reads as the baseMean and alpha at 0.1

void ihw(string fileName, SetupComp SetupValues)
{
	cout << "\tcalculating Q values using Independent Hypothesis Weighting" << endl;
	cout << "\t- N.B. if IHW issues a \"Only 1 bin; IHW reduces to Benjamini Hochberg\" warning, you should re-run specifying the BH algorithm!" << endl;

	string cmd = string("Rscript --vanilla -e 'args = commandArgs(trailingOnly=TRUE) \n");
	cmd += "library(IHW) \n";
	cmd += "binom <- read.table(args[1], header= TRUE, sep=\"\\t\") \n";
	cmd += "baseMean <- (binom$expected+binom$readCount)/2 \n";
	cmd += "binom$qvalue <- adj_pvalues(ihw(binom$pvalue ~ baseMean, alpha = as.numeric(args[2]) )) \n";
	cmd += "write.table(binom,file=args[1], sep=\"\\t\",quote=FALSE,row.names = FALSE) ' ";
	cmd += fileName + " " + SetupValues.getAlpha();

	if (SetupValues.getVerbose())
		cout << "\nRscript call:\n" << cmd << endl << endl;
	//int sys = system(cmd.c_str());
	//cout << "system: " << sys <<endl;
	std::array<char, 128> buffer;
	    std::string result;
	FILE *pipe = popen(cmd.c_str(), "r");
    while (fgets(buffer.data(), 128, pipe) != NULL) {
        std::cout << "\t\tReading..." << std::endl;
        result += buffer.data();
    }
    auto returnCode = pclose(pipe);
    if ( result.find("Only 1 bin") != string::npos)
    	throw std::invalid_argument("\nIHW internal BH correction inaccurate!\nPlease re-run with BH option!");
	cerr << "\t" << flush;
	completed();
}
