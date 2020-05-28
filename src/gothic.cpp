/*
 * gothic.cpp
 *
 *  Created on: 26 May 2020
 *      Author: rich
 */

#include "gothic.h"
//#include "Setup.h"
#include "binTest.h"
#include "Utils.h"
#include "BinomData.h"
#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <string>



//#include <seqan3/core/debug_stream.hpp>   // for debug_stream
//#include <seqan3/std/filesystem>          // for tmp_dir
//#include <seqan3/argument_parser/all.hpp> // for argument_parser

using namespace std;


int main(int argc, char *argv[])
{
	if ( argc != 2 ) // argc should be 2 for correct execution
	{
		// We print argv[0] assuming it is the program name
		cout<<"usage: "<< argv[0] <<" <filename>\n";
		printUsage();
		return 0;
	}

	vector<string> allArgs(argv, argv + argc);
	Setup setupValues = loadConfig(allArgs[1]);

	setupValues.print();

	omp_set_num_threads(setupValues.getThreads());

	vector<BinomData> binom;

	try {
		gothicHicup(setupValues,binom);
		//binom = gothicHicup(fileName, sampleName, res, restrictionFile, cistrans, parallel);
		/*
		 * vector<Site> fragments;
		 * binaryWriteTest(fragments, restrictionFile);
		 * binaryRead(fragments);
		 */
		//returnSizes();
		//timeTest();

	}
	catch(const std::invalid_argument& e){
		cerr << "Error: " << e.what() << endl;
	}


	if (setupValues.getAnalysisType() == ao_single)
	{
		sort(binom.begin(), binom.end(), bincomp);

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
	}

	return 0;
}


// GOTHiC main tool

/*
gothic(string fileName1, string fileName2, string sampleName, res, string genome=BSgenome.Hsapiens.UCSC.hg19, string restrictionSite='A^AGCTT', string enzyme='HindIII', cistrans='all',int filterdist=10000, int DUPLICATETHRESHOLD=1, fileType='BAM', parallel=FALSE, cores=NULL){
	if(file.exists(paste(sampleName,"paired_filtered",sep="_"))){
		message("Loading paired reads file ...")
		load(paste(sampleName,"paired_filtered",sep="_"))
		pairedReadsFile <- filtered

	}else{
		message("Pairing reads")

		pairedReadsFile <- pairReads(fileName1, fileName2, sampleName, DUPLICATETHRESHOLD=1, fileType=fileType)
	}

	vector<interaction> interactions;

	if (file.exists(paste("interactingLoci",sampleName,sep="_"))){
	   message("Loading mapped reads file ...")
	   load(paste("interactingLoci",sampleName,sep="_"))
	   interactions <- interactingLoci
	}else{
		message("Mapping reads to restriction sites")
		interactions <- mapReadsToRestrictionSites(pairedReadsFile, sampleName, BSgenomeName, genome, restrictionSite, enzyme, parallel, cores)
	}

	message("Computing binomial ...")
	vector<BinomData> binom;
	.binomialHiC(interactions, res,sampleName, BSgenomeName, genome, restrictionSite, enzyme, parallel, cores,cistrans=cistrans,filterdist=filterdist)
	return(binom)
}
*/

// GOTHiC main tool based on hicup alignment

void gothicHicup(Setup & setupValues, vector<BinomData> & binom)
{
	/*
	 * fileName
	 *          A character string with the name of the file containing the mapped,
	 *          filtered reads from HiCUP, after the default HiCUP output is converted
	 *          to a table containing only the read ID, chromosome, start and end
	 *          positions columns. Can be gzipped. (Tab separated text format)
	 * sampleName
	 *          A character string that will be used to name the quality control plot.
	 *          It will be saved in the current directory.
	 * res
	 *          An integer that gives the required bin size or resolution of the contact
	 *          map e.g. 1000000, for fragment level use 1.
	 * restrictionFile
	 *          A character string with the name of the digest file from HiCUP. It is
	 *          used to map reads to restriction fragments. (.txt file name)
	 * cistrans
	 *          A character string with three possibilities. "all" runs the binomial
	 *          test on all interactions, "cis" runs the binomial test only on
	 *          intrachromosomal/cis interactions, "trans" runs the binomial test only
	 *          on interchromosomal/trans interactions.
	 * parallel
	 *          Logical argument. If TRUE the mapping and the binomial test will be
	 *          performed faster using multiple cores. The default is FALSE.
	 * cores
	 *          An integer specifying the number of cores used in the parallel processing
	 *          if parellel=TRUE. The default is NULL.
	 */


	// Get Hicup data
    vector<Interaction> interactions;

    importHicup(setupValues.getInput(), interactions);

    //cout << "fast value" << setupValues.getFast() << endl;

    const bool fast = false; // set true to use map version and false to use Vector version

    typedef conditional< (fast == true),
    		multimap<string,array<int,2>>,
			std::vector<Site> >::type fragType;
    //std::vector<Site> fragments; //** Site based run **//
    //multimap<string,array<int,2>> fragments; //** map based run **//

    fragType fragments;
    getHindIIIsitesFromHicup(fragments, setupValues.getEnzyme());

    mapHicupToRestrictionFragment(interactions, fragments);

	binInteractions(interactions, setupValues.getRes());

	switch (setupValues.getAnalysisType())
	{
	case ao_single :
		// Prepare Binominal data from Hicup Data
		binomialHiChicup(interactions, setupValues, binom);
		break;
	case ao_comparative:
		// save binary file of interactions
		string outBinName = setupValues.getSname()+".inter.bin";
		writeBinary(interactions, outBinName);

		//vector<Interaction> interactions2;
		//readBinary(interactions2, "test.inter.bin");

		break;
	}

}
