/*
 * hicupDataComp.cpp
 *
 *  Created on: 26 May 2020
 *      Author: rich
 */

#include "hicupDataComp.h"
#include "hicupData.h"
#include <algorithm>
#include <math.h> //pow
#include <stdint.h> // uint32_t


using namespace std;

void binomialHiChicupComp(vector<Interaction> & interactions1, vector<Interaction> & interactions2, SetupComp & setupValues, vector<BinomDataComp> & binFiltered)
{
	/*
	 * interaction set1 that will be used as background
	 */

	cerr << "Binomial HiC Hicup Comparative Analysis" << endl;

	removeDiagonals(interactions1, setupValues.getCisTrans(), setupValues.getRemoveDiagonal());
	removeDiagonals(interactions2, setupValues.getCisTrans(), setupValues.getRemoveDiagonal());

	if (setupValues.getVerbose())
	{
		cerr << "\tControl: " << interactions1.size() << " interactions" <<endl;
		cerr << "\tSample:  " << interactions2.size() << " interactions" <<endl;
	}

	cerr << "\tGetting Pairs" << endl;
	set<string> Int1; // Ints from Control
	set<string> Int2; // Ints from Sample
	set<string> AllInt; // all_bin
	map<string,int> pairs1;
	map<string,int> pairs2;

	for (auto e : interactions1)
	{
		Int1.insert(e.getInt1());
		Int1.insert(e.getInt2());
		AllInt.insert(e.getInt1());
		AllInt.insert(e.getInt2());
		pairs1.insert(make_pair(e.getInt1() + ":" + e.getInt2(),e.getFreq()));
	}

	for (auto e : interactions2)
	{
		Int2.insert(e.getInt1());
		Int2.insert(e.getInt2());
		AllInt.insert(e.getInt1());
		AllInt.insert(e.getInt2());
		pairs2.insert(make_pair(e.getInt1() + ":" + e.getInt2(),e.getFreq()));
	}
	cerr << "\t" << flush;
	completed();

	if (setupValues.getVerbose())
	{
		cout << "\tControl Ints: " << Int1.size() << " unique positions, " << pairs1.size() << " pairs" << endl;
		cout << "\tSample Ints:  " << Int2.size() << " unique positions, " << pairs2.size() << " pairs" << endl;
		cout << "\tAll Ints:     " << AllInt.size() << " unique positions"  << endl;
	}

	set<string> only1;
	int numberOfReadPairs1 = 0;
	for (auto a :pairs1)
	{
		string e = a.first;
		if (pairs2.find(e) == pairs2.end())
		{
			only1.insert(e);
			Interaction I = splitPair(e);
			numberOfReadPairs1 += 1;
			//cout << I << endl;
			interactions2.push_back(I);
		}
		else
		{
			numberOfReadPairs1 += a.second;
		}
	}//*/
	set<string> only2;
	int numberOfReadPairs2 = 0;
	for (auto a :pairs2)
	{
		string e = a.first;
		if (pairs1.find(e) == pairs1.end())
		{
			only2.insert(e);
			Interaction I = splitPair(e);
			numberOfReadPairs2 += 1;
			//cout << I << endl;
			interactions1.push_back(I);
		}
		else
		{
			numberOfReadPairs2 += a.second;
		}
	}
	if (setupValues.getVerbose())
	{
		cout << "\t" << only1.size() << " interactions unique to Control" << endl;
		cout << "\t" << only2.size() << " interactions unique to Sample" << endl;
	}

	if (interactions1.size() != interactions2.size())
		throw std::invalid_argument("binomialHiChicupComp: imbalanced interactions");

	sort(interactions1.begin(),interactions1.end(),intcomp);
	sort(interactions2.begin(),interactions2.end(),intcomp);

	if (setupValues.getVerbose())
	{
		cerr << "\tControl: " << interactions1.size() << " interactions" << endl;
		cerr << "\tSample:  " << interactions2.size() << " interactions" << endl;
	}

	float covS = AllInt.size(); // === length(all.bins)
	// numberOfAllInteractions <- length(all.bins)^2
	double numberOfAllInteractions = pow(covS,2);
	// upperhalfBinNumber <- (length(all.bins)^2-length(all.bins))/2
	double upperhalfBinNumber = (numberOfAllInteractions - covS)/2;

	if (setupValues.getVerbose())
	{
		cerr.precision(15);
		cerr << "\tNumber of All Int:  " << covS << " (R = 255197)" << endl;
		cerr << "\tTotal Interactions: " << numberOfAllInteractions << " ( R = 65125508809)" << endl;
		cerr << "\tUpperHalfBinNumber: " << upperhalfBinNumber << " (R = 32562626806)"<< endl;
	}

	/*set<string> chromos;
	for (auto i : interactions2)
	{
		chromos.insert(i.getChr1());
	}//*/

	double sumSquare = 0;
	map<string,set<int>> loci;
	for(auto a: interactions2)
	{
		loci[a.getChr1()].insert(a.getLocus1());
		loci[a.getChr2()].insert(a.getLocus2());
         //<- length(unique(c(unique(binned_df_filtered2$locus1[binned_df_filtered2$chr1==cr]),unique(binned_df_filtered2$locus2[binned_df_filtered2$chr2==cr]))))
		 //chrlens[a] = loci.size();
	}
	for (auto i : loci)
	{
		sumSquare += pow(i.second.size(),2);
				//(sum(chrlens^2)
	}
	double cisBinNumber = (sumSquare - covS)/2;
	double transBinNumber = upperhalfBinNumber - cisBinNumber;

	if (setupValues.getVerbose())
	{
		cerr.precision(15);
		cerr << "\tSum of Squares:  " << sumSquare << "(R = 3409581821)" << endl;
		cerr << "\tNumber of Cis:   " << cisBinNumber  << " (R = 1704663312)" << endl;
		cerr << "\tNumber of Trans: " << transBinNumber << " (R = 30857963494)" << endl;
	}

	completed();
}

/*
 * .binomialHiChicupComp=function(hicupinteraction1,hicupinteraction2, restrictionFile, sampleName, cistrans='all', parallel=FALSE, cores=8, removeDiagonal=TRUE,baits=NULL,res=NULL)
	{


#all read pairs used in binomial: numberOfReadPairs1 used for computing expected and numberOfReadPairs2 as total number

		chrlens <- c()
		for(cr in chromos){
             chrlens[cr] <- length(unique(c(unique(binned_df_filtered2$locus1[binned_df_filtered2$chr1==cr]),unique(binned_df_filtered2$locus2[binned_df_filtered2$chr2==cr]))))
		}
		cisBinNumber <-(sum(chrlens^2)-length(all.bins))/2
		transBinNumber <- upperhalfBinNumber-cisBinNumber
		print(cisBinNumber)
		print(transBinNumber)
		binned_df_filtered2$probabilityOfInteraction <- binned_df_filtered1$frequencies/numberOfReadPairs1
		freq1 <- binned_df_filtered1$frequencies>1
		freq2 <- binned_df_filtered2$frequencies>1
		binned_df_filtered2 <- binned_df_filtered2[freq1 & freq2,]
		binned_df_filtered1 <- binned_df_filtered1[freq1 & freq2,]
		print(dim(binned_df_filtered2))
		#rm(binned_df_filtered1)
		#print("binned_df_filtered1 has been removed")
		if (!is.null(baits)){
		   load(baits) # baitGR
		   #print(baitGR)
		   locus1GR <- GRanges(binned_df_filtered2$chr1,ranges=IRanges(binned_df_filtered2$locus1,(binned_df_filtered2$locus1+res)))
		   locus2GR <- GRanges(binned_df_filtered2$chr2,ranges=IRanges(binned_df_filtered2$locus2,(binned_df_filtered2$locus2+res)))
		   #print(locus1GR)
		   #print(locus2GR)
		   ovl1 <- findOverlaps(locus1GR,baitGR,select="first",ignore.strand=T)
		   print(length(which(!is.na(ovl1))))
		   ovl2	<- findOverlaps(locus2GR,baitGR,select="first",ignore.strand=T)
		   print("ovl OK")

		   binned_df_filtered2[,"bait1"] <- NA
		   binned_df_filtered2[!is.na(ovl1),"bait1"] <- baitGR$names[ovl1[!is.na(ovl1)]]
		   binned_df_filtered2[,"bait2"] <- NA
                   binned_df_filtered2[!is.na(ovl2),"bait2"] <-	baitGR$names[ovl2[!is.na(ovl2)]]
		   binned_df_filtered2 <- binned_df_filtered2[!is.na(ovl1) | !is.na(ovl2),]
		   binned_df_filtered1 <- binned_df_filtered1[!is.na(ovl1) | !is.na(ovl2),]
		   print(binned_df_filtered2[,c("bait1","bait2")])
		}
		binned_df_filtered2[,"predicted"] <- binned_df_filtered1$frequencies

			if(nrow(binned_df_filtered2)>1e8)
				{
					 t <- ceiling(nrow(binned_df_filtered2)/1e8)
					 dfList <- list()
					 dfList[[1]] <- binned_df_filtered2[1:1e8,]
					 for(i in 2:t){
					 dfList[[i]] <- binned_df_filtered2[(((i-1)*1e8)+1):min((i*1e8),nrow(binned_df_filtered2)),]
				}
			pvalues=list()
			for(i in 1:length(dfList)){
				pvalues[[i]] <-apply(dfList[[i]], 1, function(x)
										  {
										  binom.test(as.numeric(x[["frequencies"]])-1, numberOfReadPairs2, as.numeric(x[["probabilityOfInteraction"]]), alternative = "two.sided")$p.value
										  }
										  )
					 }
					 pvals=unlist(pvalues)
					 binned_df_filtered2$pvalue <- pvals
					 }else{

					 binned_df_filtered2$pvalue <- apply(binned_df_filtered2, 1, function(x)
														{
														binom.test(as.numeric(x[["frequencies"]])-1, numberOfReadPairs2, as.numeric(x[["probabilityOfInteraction"]]), alternative = "two.sided")$p.value
														}
														)
			}

#observed over expected log ratio
		binned_df_filtered2$logFoldChange <- log2((binned_df_filtered2$frequencies/numberOfReadPairs2)/(binned_df_filtered1$frequencies/numberOfReadPairs1))
#multiple testing correction separately for matrices with all interactions/only cis/only transs

		if(cistrans=='all'){
			binned_df_filtered2$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered2$pvalue, method = "BH", n=upperhalfBinNumber)}else{p.adjust(binned_df_filtered2$pvalue, method = "BH", n=upperhalfBinNumber+length(all.bins))}
		}
		if(cistrans=='cis'){
			binned_df_filtered2$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered2$pvalue, method = "BH", n=cisBinNumber)}else{p.adjust(binned_df_filtered2$pvalue, method = "BH", n=cisBinNumber+length(all.bins))}
		}
		if(cistrans=='trans'){
			binned_df_filtered2$qvalue <- p.adjust(binned_df_filtered2$pvalue, method = "BH", n=transBinNumber)
		}

		test <- as.data.frame(binned_df_filtered2)

		test[,"pvalue"] <- test$pvalue

			pval.plot <- ggplot(test,aes(x=pvalue))
			tryCatch(
			{
			dev.new()
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

	binned_df_filtered=binned_df_filtered2[,c('chr1','locus1','chr2','locus2','probabilityOfInteraction', 'frequencies','predicted', 'pvalue','qvalue','logFoldChange','bait1','bait2')]
	colnames(binned_df_filtered)=c('chr1','locus1','chr2','locus2','probability','readCount','expected' ,'pvalue','qvalue','logObservedOverExpected','bait1','bait2')


return(binned_df_filtered)
}
 *
 */

Interaction splitPair(string & e)
{
	size_t pos = e.find(":");
	string chr1 = e.substr(0,pos);

	e = e.substr(pos+1);
	pos = e.find(":");
	int locus1 = atoi(e.substr(0,pos).c_str());

	e = e.substr(pos+1);
	pos = e.find(":");
	string chr2 = e.substr(0,pos);
	int locus2 = atoi(e.substr(pos+1).c_str());

	Interaction I = Interaction(chr1,chr2,locus1,locus2,1);
	return I;
}
