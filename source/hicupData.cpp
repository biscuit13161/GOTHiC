/*
 * hicupData.cpp
 *
 *  Created on: 5 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "hicupData.h"
#include "Utils.h"
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
//#include <execution>
//#include <seqan3/io/alignment_file/all.hpp>
//#include "gzstream.h"

using namespace std;

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
		throw std::invalid_argument("importHicup: doesn't function with bam files, input can be converted to appropriate text file using hicupToTable script");
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

	//if(fileType=="table")
	if (inFile.is_open())
	{
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
		inFile.close();
	}
}

void mapHicupToRestrictionFragment(vector<Interaction> & interactions, string restrictionFile)
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

	std::vector<Site> fragments;
	getHindIIIsitesFromHicup(fragments, restrictionFile);
	cerr << "Loaded " << fragments.size() << " enzyme fragments" << endl;;

	findOverlaps(sources, fragments, "Sources");
	findOverlaps(targets, fragments, "Targets");

	cerr << "Identified " << 2*sources.size() << " overlaps" << endl;;

	//sort positions into order
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
	completed();
}

void binInteractions(vector<Interaction> & interactions, int res)
{
	cerr << "Calculating Interaction Bins ("  <<interactions.size() << ")" << endl;

	vector<Interaction> interactions2;

	#pragma omp parallel for
	for (int i = 0; i < interactions.size(); i++)
	{
		int T1 = interactions[i].getLocus1();

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

	//return 0;
	cerr << "\t... Completed!" << endl;
}

vector<BinomData> binomialHiChicup(vector<Interaction> & interactions, string restrictionFile, string sampleName, CisTrans cistrans, bool parallel, int cores, bool removeDiagonal)
{
	vector<Site> hindGR;
	getHindIIIsitesFromHicup(hindGR, restrictionFile);
	/*	if(parallel)
			{
				requireNamespace("parallel")
				print("running garbage collector before parallel fork")
				gc()
			}//*/

	//vector<Interaction> binned_df_filtered;
	//binned_df_filtered = interactions;
	set<string> all_bins;

	//diagonal removal
	if(removeDiagonal)
	{
		//for (auto it = interactions.begin(); it != interactions.end(); it++)
		for_each (interactions.begin(), interactions.end(), [&](Interaction & T)
		{
			//Interaction T = *it;
			if (T.getInt1() == T.getInt2())
			{
				//binned_df_filtered <- binned_df_filtered[binned_df_filtered$int1!=binned_df_filtered$int2,]
				all_bins.insert(T.getInt1());
			}
		});
	}
	if(cistrans == ct_cis){
		for_each (interactions.begin(), interactions.end(), [&](Interaction & T)
		{
			if (T.getChr1() == T.getChr2())
			{
				//binned_df_filtered <- binned_df_filtered[binned_df_filtered$chr1==binned_df_filtered$chr2,]
				all_bins.insert(T.getInt1());
				all_bins.insert(T.getInt2());
			}
		});
	}
	else if(cistrans == ct_trans){
		for_each (interactions.begin(), interactions.end(), [&](Interaction & T)
		{
			if (T.getChr1() != T.getChr2())
			{
				//binned_df_filtered <- binned_df_filtered[binned_df_filtered$chr1!=binned_df_filtered$chr2,]
				all_bins.insert(T.getInt1());
				all_bins.insert(T.getInt2());
			}
		});
	}

	// all read pairs used in binomial
	int numberOfReadPairs = interactions.size(); // before binning!!

	// calculate coverage
		/*
		if( interations.size() > 1e8)
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
			}else{
				binned_dt=data.table(binned_df_filtered)
				covA <- binned_dt[,sum(frequencies),by=int1]
				covB <- binned_dt[,sum(frequencies),by=int2]
			}
		covA <- setkey(covA,key='int1')
		setnames(covB, 1,'int1')
		covB <- setkey(covB,key='int1')

		cov=merge(covA,covB,all.x=TRUE,all.y=TRUE,by='int1')
		cov$V1.x[is.na(cov$V1.x)]=0
		cov$V1.y[is.na(cov$V1.y)]=0
		cov$coverage=cov$V1.x+cov$V1.y
		coverage=cov$coverage
		names(coverage)=cov$int1
		sumcov <- sum(coverage)
		relative_coverage <- coverage/sumcov
		names(relative_coverage)=names(coverage)
		binned_df_filtered$coverage_source <- relative_coverage[binned_df_filtered$int1]
		binned_df_filtered$coverage_target <- relative_coverage[binned_df_filtered$int2]

#probability correction assuming on average equal probabilities for all interactions
		numberOfAllInteractions <- length(all_bins)^2
		upperhalfBinNumber <- (length(all_bins)^2-length(all_bins))/2

		if(cistrans!=ct_all){
			chromos <- unique(binned_df_filtered$chr1)
			chrlens <- c()
			for(cr in chromos){
                chrlens[cr] <- length(unique(c(unique(binned_df_filtered$locus1[binned_df_filtered$chr1==cr]),unique(binned_df_filtered$locus2[binned_df_filtered$chr2==cr]))))
			}
			cisBinNumber <-(sum(chrlens^2)-length(all_bins))/2
			transBinNumber <- upperhalfBinNumber-cisBinNumber
		}

		diagonalProb <- sum(relative_coverage^2)
		if(cistrans==ct_all){
				probabilityCorrection <- if(removeDiagonal){1/(1-diagonalProb)}else{1}
		}
		if(cistrans==ct_cis){
			probabilityCorrection <- upperhalfBinNumber/cisBinNumber
		}
		if(cistrans==ct_trans){
			probabilityCorrection <- upperhalfBinNumber/transBinNumber
		}


		binned_df_filtered$probabilityOfInteraction <- binned_df_filtered$coverage_source*binned_df_filtered$coverage_target*2*probabilityCorrection


		binned_df_filtered$predicted <- binned_df_filtered$probabilityOfInteraction * numberOfReadPairs

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
		}else{
			if(nrow(binned_df_filtered)>1e8)
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
					 }else{

					 binned_df_filtered$pvalue <- apply(binned_df_filtered, 1, function(x)
														{
														binom.test(as.numeric(x[["frequencies"]])-1, numberOfReadPairs, as.numeric(x[["probabilityOfInteraction"]]), alternative = "greater")$p.value
														}
														)
			}
		}

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
	//omp_set_dynamic(0);
	//omp_set_num_threads(8);
	int i;

	#pragma omp parallel for
	for (i = 0; i < query.size(); i++)
	{
		for (int j = 0; j < fragments.size(); j++)
		{
			if (query[i].getChr() == fragments[j].getChr())
			{
				if ((query[i].getLocus() - fragments[j].getStart()) * (fragments[j].getEnd() - query[i].getLocus()) >= 0)
				{
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

//countDuplicates <- function(df, frequencyCol = c("frequencies", "frequency"), considerExistingFrequencies = any(frequencyCol %in% names(df)), ordered=FALSE)
void countDuplicates(vector<Interaction> & interactions)
{
	showTime();
	cerr << "Counting Duplicates" << endl;
	/*
	if(considerExistingFrequencies)
	{
		frequencyCol <- frequencyCol[frequencyCol %in% names(df)][1]
		if(is.null(df[[frequencyCol]]))
			stop("frequencies column missing")
	} else if(!considerExistingFrequencies)
	{
		frequencyCol <- frequencyCol[1]
		df[, frequencyCol] <- 1 #one for all
	}//*/


	set<string> int1s;
	set<string> int2s;
	map<string , map<string, int>> list;

	for (int i = 0; i <= interactions.size(); i++)
	{
		string int1 = interactions[i].getInt1();
		string int2 = interactions[i].getInt2();
		list[int1][int2]++;
		int1s.insert(int1);
		int2s.insert(int2);
		//cout << int1 << "\t" << int2 << "\t" << list[int1][int2] << endl;
	}
	interactions.clear();
	cerr << "Counting Duplicates2" << endl;

	for (auto it = list.begin(); it != list.end(); it++)
	{
		string int1 = it->first;
		cout << "test1 " << int1 << endl;
		size_t pos = int1.find(":");
		string chr1 = int1.substr(0,pos);
		int locus1 = atoi(int1.substr(pos+1).c_str());
 		//cout << "\t" << chr1 << ":" << locus1 << endl;
		map<string, int> l2 = it->second;
		for (auto ti = l2.begin(); ti != l2.end(); ti++)
		{
			int f = ti->second;

			string int2 = ti->first;

			size_t pos = int2.find(":");
			string chr2 = int2.substr(0,pos);
			int locus2 = atoi(int2.substr(pos+1).c_str());

			Interaction I = Interaction(chr1, chr2, locus1, locus2, f);
			I.print();
			cout << "  test3" << endl;


			//cout << "  test4" << endl;
			//interactions.push_back(I);
			//string mInt2 = ti->first;
			//int freq = ti->second;
		}//*/

	}

	showTime();
	cout << "Duplicates: " << interactions.size() << endl;

	completed();
}



