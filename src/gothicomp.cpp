/*
 * gothicomp.cpp
 *
 *  Created on: 26 May 2020
 *      Author: rich
 */

#include "gothicomp.h"
#include "Interactions.h"
#include "hicupDataComp.h"
#include "binTest.h"
#include <iostream>
#include <stdio.h>
#include <omp.h>
#include "tbb/concurrent_vector.h"

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

	omp_set_num_threads(setupValues.getThreads());

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
			<< "\t" << "relCoverage1" \
			<< "\t" << "relCoverage2" \
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


/*
 * GOTHiChicupComp <- function(fileName1,fileName2,sampleName,res,restrictionFile, cistrans='all', parallel=FALSE, cores=NULL,baits=NULL){
	hindGR <- .getHindIIIsitesFromHicup(restrictionFile)
	print(hindGR)
	if (file.exists(paste("/scratch/cbib/HiC_analysis/CHiC_Tcell/results/G0_binned_",sampleName,".RData",sep=""))){
	   message("Loading file 1")
	   load(paste("/scratch/cbib/HiC_analysis/CHiC_Tcell/results/G0_binned_",sampleName,".RData",sep=""))
	}
	else{
		 message("importing file 1")
		condition1 <- .importHicup(fileName1, checkConsistency=TRUE, fileType=ifelse(grepl("\\.bam$", fileName1)|grepl("\\.sam$", fileName1), "bam", "table"))
		message("mapping to restriction file 1")

		condition1 <- .mapHicupToRestrictionFragment(condition1, hindGR)
		print(head(condition1))
		message("binning file 1")
		condition1 <- .binInteractions(condition1, res)

		save(condition1,file=paste("/scratch/cbib/HiC_analysis/CHiC_Tcell/results/G0_binned_",sampleName,".RData",sep=""))
	}
	if (file.exists(paste("/scratch/cbib/HiC_analysis/CHiC_Tcell/results/G1_binned_",sampleName,".RData",sep=""))){
	   message("Loading file 2")
	   load(paste("/scratch/cbib/HiC_analysis/CHiC_Tcell/results/G1_binned_",sampleName,".RData",sep=""))
	}
	else{
		message("importing file 2")
		condition2 <- .importHicup(fileName2, checkConsistency=TRUE, fileType=ifelse(grepl("\\.bam$", fileName2)|grepl("\\.sam$", fileName2), "bam", "table"))
		message("mapping to restriction file 2")
		condition2 <- .mapHicupToRestrictionFragment(condition2, hindGR)
		message("binning file 2")
		condition2 <- .binInteractions(condition2, res)
		save(condition2,file=paste("/scratch/cbib/HiC_analysis/CHiC_Tcell/results/G1_binned_",sampleName,".RData",sep=""))
	}
	message("Computing Binomial")
	binom <- .binomialHiChicupComp(condition1,condition2, restrictionFile, sampleName, cistrans, parallel, cores, removeDiagonal=TRUE,baits,res)

	return(binom)
}

 *
 */
