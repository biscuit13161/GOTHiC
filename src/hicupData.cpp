/*
 * hicupData.cpp
 *
 *  Created on: 5 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "../src/hicupData.h"
#include "pbinom.h"
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

#include "../src/Utils.h"
//#include <execution>
//#include <seqan3/io/alignment_file/all.hpp>
//#include "gzstream.h"

using namespace std;

struct RetrieveKey
{
    template <typename T>
    typename T::first_type operator()(T keyValuePair) const
    {
        return keyValuePair.first;
    }
};

Site::Site(): mChr(""), mLocus(0), mStart(0), mEnd(0)
{
}

Site::Site(string chr, int locus, int start, int end): mChr(chr), mLocus(locus), mStart(start), mEnd(end)
{
}

Site::Site(const Site & other)
{
	mChr = other.mChr;
	mLocus = other.mLocus;
	mStart = other.mStart;
	mEnd = other.mEnd;
}

void Site::print(){
	cout << mChr << "\t"<< mLocus << "\t"<< mStart  << "\t"<< mEnd << endl;
}

void importHicup(string fileName, vector<Interaction> & interactions, bool checkConsistency)
{
	//fileType=ifelse(grepl("\\.bam$", fileName)|grepl("\\.sam$", fileName), "bam", "table")

	//seqan3::alignment_file_input fin_from_filename{filename};

	cerr << "importing HiCUP file: " << fileName << endl;

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

	vector<string> list(2);
	vector<vector<string>> file;

	//read file into temporary vector
	while(inFile)
	{
		string start1;
		getline(inFile,start1,'\n');
		if (start1.length() == 0)
		{
			break;
		}
		list[0] = start1;

		string start2;
		getline(inFile,start2,'\n');
		list[1] = start2;;

		file.push_back(list);

		if (!inFile) // corrects for last blank line
		{
			break;
		}
	}

	interactions.resize(file.size());

	//Parallel parse line data
	bool foundCondition = false;
	#pragma omp parallel for
	for (int i = 0; i < file.size(); i++)
	{
		string one = file[i][0];
		string two = file[i][1];
		//string str = string("line ") + to_string(i) + ": " + one + "\n\t" + two +"\n";
		//cout << str;
		string delimiter = "\t";
		size_t pos = one.find(delimiter);
		string id1 = one.substr(0,pos);
		if (id1.length() == 0)
		{
			throw std::invalid_argument("importHicup: reads must be paired in consecutive rows!");
		}
		one.erase(0, one.find(delimiter) + delimiter.length());
		pos = one.find(delimiter);
		string flag1 = one.substr(0,pos);
		one.erase(0, one.find(delimiter) + delimiter.length());
		pos = one.find(delimiter);
		string chr1 = one.substr(0,pos);
		one.erase(0, one.find(delimiter) + delimiter.length());
		pos = one.find(delimiter);
		string start1 = one.substr(0,pos);
		int locus1 = stoi(start1,nullptr,10);

		pos = two.find(delimiter);
		string id2 = two.substr(0,pos);
		if (checkConsistency && (id2.length() == 0 | id1 != id2))
		{
			throw std::invalid_argument("importHicup: reads must be paired in consecutive rows!");
		}
		two.erase(0, two.find(delimiter) + delimiter.length());
		pos = two.find(delimiter);
		string flag2 = two.substr(0,pos);
		two.erase(0, two.find(delimiter) + delimiter.length());
		pos = two.find(delimiter);
		string chr2 = two.substr(0,pos);
		two.erase(0, two.find(delimiter) + delimiter.length());
		pos = one.find(delimiter);
		string start2 = two.substr(0,pos);
		int locus2 = stoi(start2,nullptr,10);

		//string str = string("line ") + to_string(i) + ": " + chr1 + ":" + to_string(locus1) + " " + chr2 + ":" + to_string(locus2) +"\n";
		//cout << str;

		interactions[i] = Interaction(chr1, chr2, locus1, locus2);

	}
	//throw std::invalid_argument("finished");

}

/*
	while(inFile)
	{
		string id1;
		getline(inFile,id1,'\t');
		if (id1.length() == 0)
		{
			break;
		}
		string flag1;
		getline(inFile,flag1,'\t');
		string chr1;
		getline(inFile,chr1,'\t');
		string start1;
		getline(inFile,start1,'\n');
		int locus1 = stoi(start1,nullptr,10);

		string id2;
		getline(inFile,id2,'\t');
		if (checkConsistency && (id2.length() == 0 | id1 != id2))
		{
			throw std::invalid_argument("importHicup: reads must be paired in consecutive rows!");
		}
		string flag2;
		getline(inFile,flag1,'\t');
		string chr2;
		getline(inFile,chr2,'\t');
		string start2;
		getline(inFile,start2,'\n');
		int locus2 = stoi(start2,nullptr,10);


		Interaction a = Interaction(chr1, chr2, locus1, locus2);

		interactions.push_back(a);

		if (!inFile) // corrects for last blank line
		{
			break;
		}
	}
 */

void mapHicupToRestrictionFragment(vector<Interaction> & interactions, vector<Site> & fragments)
{
	int iSize = interactions.size();
	cerr << "Mapping HiCUP data (" << iSize << " positions) to enzyme fragments" << endl;

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

	cerr << "Identified " << 2*sources.size() << " overlaps" << endl;;

	//sort positions into order

	sortPositions(interactions, iSize, sources, targets);

	completed();
}

void sortPositions(vector<Interaction> & interactions, int iSize, vector<halfInteraction> & sources, vector<halfInteraction> & targets)
{
	#pragma omp parallel for
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
	cout << "Singles: " << interactions.size()  << endl;

	countDuplicates(interactions);

	completed();
}

void getHindIIIsitesFromHicup(vector<Site> & sites, string fileName)
{
	// load sites for HindIII restriction enzyme from HiCUP_digester
	cerr << "Loading Enzyme Restriction Sites" << endl;

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

			Site site = Site(fixChromosomeNames(chr), start, start, stoi(end,nullptr,10));

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

	completed();
	cerr << "Loaded " << sites.size() << " enzyme fragments" << endl;;
}

vector<BinomData> binomialHiChicup(vector<Interaction> & interactions, vector<Site> & fragments, string sampleName, CisTrans cistrans, bool parallel, bool removeDiagonal)
{
	//vector<Site> hindGR;
	//getHindIIIsitesFromHicup(hindGR, restrictionFile);

	vector<Interaction> binned_df_filtered;
	//binned_df_filtered.resize(interactions.size());
	set<string> all_bins;


	//diagonal removal
	if(removeDiagonal)
	{
		cerr << "Removing Diagonals!";
		removeDuplicates(interactions, binned_df_filtered);
		completed();
		interactions.clear();
	}
	else
	{
		binned_df_filtered.swap(interactions);
	}
	//cout << "size " << binned_df_filtered.size() << endl;

	if(cistrans == ct_cis){
		cerr << "Finding Cis interactions!";
		findCis(interactions, binned_df_filtered);
		//interactions.resize(pos);
		completed();
	}
	else if(cistrans == ct_trans){
		cerr << "Finding Trans interactions!";
		findTrans(interactions, binned_df_filtered);
		//interactions.resize(pos);
		completed();
	}
	else
	{
		interactions.swap(binned_df_filtered);
	}


	// all read pairs used in binomial
	int numberOfReadPairs = interactions.size(); // before binning!!
	cout << "Read Pairs " << numberOfReadPairs << endl;


	// calculate coverage

		/*if( interations.size() > 1e8)
			{
				t <- ceiling(nrow(binned_df_filtered)/1e8)
				dfList <- list()
				dfList[[1]] <- binned_df_filtered[1:1e8,]
				for(i in 2:t){
					 dfList[[i]] <- binned_df_filtered[(((i-1)*1e8)+1):min((i*1e8),nrow(binned_df_filtered)),]
					 }
				dtList <- lapply(dfList, data.table)
				covAs <- lapply(dtList, function(x) x[,sum(frequencies), by=int1])
				covBs <- lapply(dtList, function(x) x[,sum(frequencies), by=int2])
				covAm <- do.call(rbind, covAs)
				covBm <- do.call(rbind, covBs)
				covA <- covAm[,sum(V1),by=int1]
				covB <- covBm[,sum(V1),by=int2]
			}else{//*/

    map<string,int> cov;
    double tCoverage = 0;
    int max = 0;

    ofstream outFile("my_file.txt");
    for (const auto &e : interactions) outFile << e << "\n";

    calcFreq(interactions, cov, tCoverage, max);

    //cout << "cov.size: " << cov.size() << endl;
    cout << "total coverage: " << tCoverage << endl;
    //cout << "max: " << max << endl;

    if (tCoverage == 0)
    {
    	throw std::invalid_argument("binomialHiChicup: Zero coverage!");
    }

    map<string,long double> rCov;
    long double diagonalProb = 0;
    for (auto i = cov.begin(); i != cov.end(); i ++)
    {
    	long double V = i->second/tCoverage;
    	rCov[i->first] = V;
    	diagonalProb += V*V;
    	//relative_coverage <- coverage/sumcov
    }

    vector<BinomData> binFiltered;
    set<string> chromos;


    for (auto i : interactions) //.begin(); i < interactions.end(); i++)
    {
    	chromos.insert(i.getChr1());
    	string int1 = i.getInt1();
    	string int2 = i.getInt2();
    	BinomData bin = BinomData(i);
    	bin.setRelCov1(rCov[int1]);
    	bin.setRelCov2(rCov[int2]);
    }
    cout << "chromos: " << chromos.size() << endl;

    //probability correction assuming on average equal probabilities for all interactions
    int covS = cov.size();
    uint32_t numberOfAllInteractions = covS*covS;
    long double upperhalfBinNumber = (numberOfAllInteractions - cov.size())/2.0f;
    cout << "covS: " << covS << endl;
    cout << "inter: " << numberOfAllInteractions << " (" << covS << ")"<< endl;
    cout << "half: " << std::setprecision (17)<< upperhalfBinNumber << endl;

    long double cisBinNumber = 0;
    long double transBinNumber = 0;

    if (cistrans != ct_all)
    {
    	long double sumSquare = 0;
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
    	cout << "sumSquare: " << sumSquare << endl;
    	cout << "uphalfBin: " << upperhalfBinNumber << endl;

    	cisBinNumber = (sumSquare - cov.size())/2;
    	transBinNumber = upperhalfBinNumber - cisBinNumber;
        cout << "cisBin: " << std::setprecision (17) << cisBinNumber << endl;
        cout << "transBin: " << std::setprecision (17)<< transBinNumber << endl;
    }
    cout << "test3\n";



		//diagonalProb <- sum(relative_coverage^2)
    long double probabilityCorrection;
    if(cistrans == ct_all){
    	probabilityCorrection = (removeDiagonal)? (1/(1-diagonalProb)) : 1;
    	cout << "pc: " << probabilityCorrection << endl;
    }
    else if(cistrans == ct_cis){
    	probabilityCorrection = upperhalfBinNumber/cisBinNumber;
    	cout << "pc: " << probabilityCorrection << endl;
    }
    else if(cistrans == ct_trans){
    	probabilityCorrection = upperhalfBinNumber/transBinNumber;
    	cout << "pc: " << probabilityCorrection << endl;
    }

    // Calculate expected read counts
    for (auto x : binFiltered)
    {
    	x.setProbability(x.getRelCov1() * x.getRelCov2() * 2 * probabilityCorrection);
    	x.setExpected(x.getProbability() * numberOfReadPairs);
    }
/*
		if(parallel)
			{
				requireNamespace("parallel")
				if(nrow(binned_df_filtered)>1e8)
					{
					 t <- ceiling(nrow(binned_df_filtered)/1e8)
					 dfList <- list()
					 dfList[[1]] <- binned_df_filtered[1:1e8,]
					 for(i in 2:t){
					 dfList[[i]] <- binned_df_filtered[(((i-1)*1e8)+1):min((i*1e8),nrow(binned_df_filtered)),]
					}
					 dtList <- lapply(dfList, function(x) as.data.frame(t(cbind(as.numeric(x[["frequencies"]]),
																				as.numeric(x[["probabilityOfInteraction"]])))))
					 pvalues=list()
					 for(i in 1:length(dtList)){
						 pvalues[[i]] <-unlist(parallel::mclapply(dtList[[i]], function(x)
													{
													binom.test(x[1]-1, numberOfReadPairs, x[2], alternative = "greater")$p.value
													},
													mc.cores=cores))
					 }
					 pvals=unlist(pvalues)
					 binned_df_filtered$pvalue <- pvals
				}else{

					 binomParams <- as.data.frame(t(cbind(
														  as.numeric(binned_df_filtered[["frequencies"]]),
														  as.numeric(binned_df_filtered[["probabilityOfInteraction"]]
																	 ))))


					binned_df_filtered$pvalue <- unlist(parallel::mclapply(binomParams, function(x)
																  {
																  binom.test(x[1]-1, numberOfReadPairs, x[2], alternative = "greater")$p.value
																  },
																  mc.cores=cores))
				}
		}else{*/
			/*if(nrow(binned_df_filtered)>1e8)
				{
					 t <- ceiling(nrow(binned_df_filtered)/1e8)
					 dfList <- list()
					 dfList[[1]] <- binned_df_filtered[1:1e8,]
					 for(i in 2:t){
					 	 dfList[[i]] <- binned_df_filtered[(((i-1)*1e8)+1):min((i*1e8),nrow(binned_df_filtered)),]
					}
					pvalues=list()
					for(i in 1:length(dfList)){
						pvalues[[i]] <-apply(dfList[[i]], 1, function(x)
										  {
										  binom.test(as.numeric(x[["frequencies"]])-1, numberOfReadPairs, as.numeric(x[["probabilityOfInteraction"]]), alternative = "greater")$p.value
										  }
										  )
					 }
					 pvals=unlist(pvalues)
					 binned_df_filtered$pvalue <- pvals
					 }else{//*/

    int freq = 1;
    long double prob = 6.079281e-10;
    int num =  28679;
    bool alt = true;

    double P = pbinom(double(freq), double(num), prob, 0, 0);; //binomialTest(freq, num, prob, alt);
    cout << "P: " << P << endl;

	/*	binned_df_filtered$pvalue <- apply(binned_df_filtered, 1, function(x)
			{
				binom.test(as.numeric(x[["frequencies"]])-1,
						numberOfReadPairs,
						as.numeric(x[["probabilityOfInteraction"]]),
						alternative = "greater")$p.value
			}
		)

		/*
		 * P(X = k) = (n : k) pk(1-p)(n-k)

where (n : k) = (n!) ÷ (k!)(n - k)!
fixed number n observations
probability of success p
		 */
			//} //if(nrow(binned_df_filtered)>1e8) else
		//} // parallel else
/*
#observed over expected log ratio
		binned_df_filtered$logFoldChange <- log2(binned_df_filtered$frequencies/binned_df_filtered$predicted)
#multiple testing correction separately for matrices with all interactions/only cis/only transs

		if(cistrans==ct_all){
			binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=upperhalfBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=upperhalfBinNumber+length(all_bins))}
		}
		if(cistrans==ct_cis){
			binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber+length(all_bins))}
		}
		if(cistrans==ct_trans){
			binned_df_filtered$qvalue <- p.adjust(binned_df_filtered$pvalue, method = "BH", n=transBinNumber)
		}

		test <- data.frame(binned_df_filtered)
		test[,"pvalue"] <- test$pvalue
			pval.plot <- ggplot(test,aes(x=pvalue))
			tryCatch(
			{
			x11()
			print(pval.plot + geom_density())
			},
			error=function(cond) {
			message("No interactive plot, try saving image")
			message(cond)
			return(tryCatch(
				{
				pdf(file=paste(sampleName,"pvalue_distribution.pdf",sep="_"))
				print(pval.plot + geom_density())fi
				dev.off()
				},
				error=function(cond2) {
				message("No pdf output, quality assesment plot is not produced")
				message(cond2)
				return(NA)
				},
				warning=function(cond2) {
				message("No pdf output, quality assesment plot is not produced")
				message(cond2)
				return(NA)
						}
					))
				},
				warning=function(cond) {
				message("No interactive plot, try saving image")
				message(cond)
				return(tryCatch(
				{
				pdf(file=paste(sampleName,"pvalue_distribution.pdf",sep="_"))
				print(pval.plot + geom_density())
				dev.off()
				},
				error=function(cond2) {
				message("No pdf output, quality assesment plot is not produced")
				message(cond2)
				return(NA)
				},
				warning=function(cond2) {
				message("No pdf output, quality assesment plot is not produced")
				message(cond2)
				return(NA)
				}
			))
		})

	binned_df_filtered=binned_df_filtered[,c('chr1','locus1','chr2','locus2','coverage_source','coverage_target','probabilityOfInteraction', 'predicted','frequencies', 'pvalue','qvalue','logFoldChange')]
	colnames(binned_df_filtered)=c('chr1','locus1','chr2','locus2','relCoverage1','relCoverage2','probability', 'expected','readCount', 'pvalue','qvalue','logObservedOverExpected')


return(binned_df_filtered)
}
	 */

	vector<BinomData> binom;

	return binom;//
}


void findOverlaps(vector<halfInteraction>& query, vector<Site> & fragments, string name, bool drop)
{
	cerr << "Finding Overlaps in " << name << endl;
	int i;

	#pragma omp parallel for
	for (i = 0; i < query.size(); i++)
	{
		for (int j = 0; j < fragments.size(); j++)
		{
			if (query[i].getChr() == fragments[j].getChr())
			{
				if ((query[i].getLocus() >= fragments[j].getStart()) && (fragments[j].getEnd() >= query[i].getLocus()) )
				{
					/*if (query[i].getLocus() == 197036391)
					{
						string L = string("overlay: ") + query[i].getChr() + ":" + to_string(query[i].getLocus()) + " -> " + fragments[j].getChr() + ":" + to_string(fragments[j].getStart()) + "-" + to_string(fragments[j].getEnd());
						cerr << L;
					}//*/
					query[i] = halfInteraction(fragments[j].getChr(), fragments[j].getStart());
					break;
				}

			}
		}
	}
	completed();
}//*/


/*
int binomialCoefficients(int n, int k) {
   int C[k+1];
   memset(C, 0, sizeof(C));
   C[0] = 1;
   for (int i = 1; i <= n; i++) {
      for (int j = min(i, k); j > 0; j--)
         C[j] = C[j] + C[j-1];
   }
   return C[k];
}
//*/

void countDuplicates(vector<Interaction> & interactions)
{
	cerr << "Counting Duplicates" << endl;

	map<string , map<string, int>> list;

	int count = 0;
	for (int i = 0; i < interactions.size(); i++)
	{
		string int1 = interactions[i].getInt1();
		string int2 = interactions[i].getInt2();

		list[int1][int2]++;
	}
	interactions.clear();

	for (auto it = list.begin(); it != list.end(); it++)
	{
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

			interactions.push_back(I);
		}//*/

	}

	cout << "Duplicates: " << interactions.size()  << endl;

	completed();
}


void removeDuplicates(vector<Interaction> & interactions, vector<Interaction> & binned_df_filtered)
{
		int pos = 0;
		for (auto it = interactions.begin(); it != interactions.end(); it++)
		{
			Interaction T = *it;
			if (T.getInt1() != T.getInt2())
			{
				binned_df_filtered.push_back(T);
				pos++;
			}
		}
}

void findCis(vector<Interaction> & interactions, vector<Interaction> & binned_df_filtered)
{
	int pos = 0;
	for (auto it = binned_df_filtered.begin(); it != binned_df_filtered.end(); it++)
	{
		Interaction T = *it;
		if (T.getChr1() == T.getChr2())
		{
			interactions.push_back(T);
			pos++;
		}
	}
}

void findTrans(vector<Interaction> & interactions, vector<Interaction> & binned_df_filtered)
{
	int pos = 0;
	for (auto it = binned_df_filtered.begin(); it != binned_df_filtered.end(); it++)
	{
		Interaction T = *it;
		if (T.getChr1() != T.getChr2())
		{
			interactions.push_back(T);
			pos++;
		}
	}
}

void calcFreq(vector<Interaction> & interactions, map<string,int> & cov, double & tCoverage, int & max)
{
    for (int i = 0; i < interactions.size(); i ++ )
    {
    	if (cov.find(interactions[i].getInt1()) == cov.end())
    	{
    		//cout << "CovA " << interactions[i].getInt1() << endl;
    		cov[interactions[i].getInt1()] = interactions[i].getFreq();
    		tCoverage += interactions[i].getFreq();
    	}
    	else
    	{
    		//cout << "CovA2 " << interactions[i].getInt1() << endl;
    		cov[interactions[i].getInt1()] += interactions[i].getFreq();
    		tCoverage += interactions[i].getFreq();
    	}
    	if (cov.find(interactions[i].getInt2()) == cov.end())
    	{
    		//cout << "CovB " << interactions[i].getInt1() << endl;
    		cov[interactions[i].getInt2()] = interactions[i].getFreq();
    		tCoverage += interactions[i].getFreq();
    	}
    	else
    	{
    		//cout << "CovB2 " << interactions[i].getInt1() << endl;
    		cov[interactions[i].getInt2()] += interactions[i].getFreq();
    		tCoverage += interactions[i].getFreq();
    	}
    	max = (cov[interactions[i].getInt1()]> max )? cov[interactions[i].getInt1()]: max;
    	max = (cov[interactions[i].getInt2()]> max )? cov[interactions[i].getInt2()]: max;
    }
    //cout << "Loop Coverage: " << tCoverage << endl;
    //cout << "Loop Max: " << max << endl;
}
