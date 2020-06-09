/*
 * hicupDataComp.cpp
 *
 *  Created on: 26 May 2020
 *      Author: rich
 */

#include "hicupDataComp.h"
#include "Interactions.h"
#include "UtilsComp.h"
#include "binomTest.h"
#include "padjust.h"
#include <algorithm>
#include <math.h> //pow
#include <stdint.h> // uint32_t
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_sort.h"
#include "tbb/queuing_mutex.h"


using namespace std;
using namespace tbb;

typedef queuing_mutex Mutex;

void binomialHiChicupComp(concurrent_vector<Interaction> & interactions1, concurrent_vector<Interaction> & interactions2, SetupComp & setupValues, concurrent_vector<BinomDataComp> & binFiltered)
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
	concurrent_unordered_set<string> Int1; // Ints from Control
	concurrent_unordered_set<string> Int2; // Ints from Sample
	concurrent_unordered_set<string> AllInt; // all_bin
	concurrent_unordered_map<string,int> pairs1;
	concurrent_unordered_map<string,int> pairs2;


	parallel_for(
			blocked_range<concurrent_vector<Interaction>::iterator>(interactions1.begin(),	interactions1.end()),
			[&] (blocked_range<concurrent_vector<Interaction>::iterator> inter) {
		for (concurrent_vector<Interaction>::iterator it = inter.begin(); it != inter.end(); it++) {
			Interaction e = (*it);
			Int1.insert(e.getInt1());
			Int1.insert(e.getInt2());
			AllInt.insert(e.getInt1());
			AllInt.insert(e.getInt2());
			pairs1.insert(make_pair(e.getInt1() + ":" + e.getInt2(),e.getFreq()));
		}
	});//*/

	parallel_for(
			blocked_range<concurrent_vector<Interaction>::iterator>(interactions2.begin(),	interactions2.end()),
			[&] (blocked_range<concurrent_vector<Interaction>::iterator> inter) {
		for (concurrent_vector<Interaction>::iterator it = inter.begin(); it != inter.end(); it++) {
			Interaction e = (*it);
			Int2.insert(e.getInt1());
			Int2.insert(e.getInt2());
			AllInt.insert(e.getInt1());
			AllInt.insert(e.getInt2());
			pairs2.insert(make_pair(e.getInt1() + ":" + e.getInt2(),e.getFreq()));
		}
	});

	cerr << "\t" << flush;
	completed();

	if (setupValues.getVerbose())
	{
		cout << "\tControl Ints: " << Int1.size() << " unique positions, " << pairs1.size() << " pairs" << endl;
		cout << "\tSample Ints:  " << Int2.size() << " unique positions, " << pairs2.size() << " pairs" << endl;
		cout << "\tAll Ints:     " << AllInt.size() << " unique positions"  << endl;
	}

	concurrent_unordered_set<string> only1;
	parallel_for(pairs1.range(),
			[&] (decltype( pairs1)::range_type & inter) {
		for (concurrent_unordered_map<string,int>::iterator it = inter.begin(); it != inter.end(); it++)
		{
			string e = it->first;
			if (pairs2.find(e) == pairs2.end())
			{
				only1.insert(e);
				Interaction I = splitPair(e);
				interactions2.push_back(I);
			}
		}
	});
	cout << "\tFinished finding [mis]matching pairs in Control!" << endl;

	concurrent_unordered_set<string> only2;
	parallel_for(pairs2.range(),
			[&] (decltype( pairs2)::range_type & inter) {
		for (concurrent_unordered_map<string,int>::iterator it = inter.begin(); it != inter.end(); it++)
		{
			string e = it->first;
			if (pairs1.find(e) == pairs1.end())
			{
				only2.insert(e);
				Interaction I = splitPair(e);
				interactions1.push_back(I);
			}
		}
	});

	cout << "\tFinished finding [mis]matching pairs in Sample!" << endl;

	int numberOfReadPairs1 = 0;
	int max1 = 0;
	int min1 = 100;
	for (auto e: interactions1)
		{
		int l = e.getFreq();
		numberOfReadPairs1 += l;
		max1 = (l > max1)? l : max1;
		min1 = (l < min1)? l : min1;
		}

	int numberOfReadPairs2 = 0;
	int max2 = 0;
	int min2 = 100;
	for (auto e: interactions2)
		{
		int l = e.getFreq();
		numberOfReadPairs2 += l;
		max2 = (l > max2)? l : max2;
		min2 = (l < min2)? l : min2;
		}

	if (setupValues.getVerbose())
	{
		cout << "\tnumberOfReadPairs1: " << numberOfReadPairs1  << endl;
		cout << "\tmax: " << max1 << "\tmin: " << min1 << endl;
		cout << "\tnumberOfReadPairs2: " << numberOfReadPairs2  << endl;
		cout << "\tmax: " << max2 << "\tmin: " << min2 << endl;
		cout << "\t" << only1.size() << " interactions unique to Control" << endl;
		cout << "\t" << only2.size() << " interactions unique to Sample" << endl;
	}

	if (interactions1.size() != interactions2.size())
		throw std::invalid_argument("binomialHiChicupComp: imbalanced interactions");

	parallel_sort(interactions1.begin(),interactions1.end(),intcomp);
	parallel_sort(interactions2.begin(),interactions2.end(),intcomp);

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
		cerr << "\tNumber of All Int:  " << covS << endl;
		cerr << "\tTotal Interactions: " << numberOfAllInteractions  << endl;
		cerr << "\tUpperHalfBinNumber: " << upperhalfBinNumber << endl;
	}

	set<string> chromos;
	for (auto i : interactions2)
	{
		chromos.insert(i.getChr1());
	}//*/

	double sumSquare = 0;
	getSumSquare(sumSquare, chromos, interactions2);

	double cisBinNumber = (sumSquare - covS)/2;
	double transBinNumber = upperhalfBinNumber - cisBinNumber;

	if (setupValues.getVerbose())
	{
		cerr.precision(15);
		cerr << "\tSum of Squares:  " << sumSquare << endl;
		cerr << "\tNumber of Cis:   " << cisBinNumber << endl;
		cerr << "\tNumber of Trans: " << transBinNumber << endl;
	}

	/** all read pairs used in binomial **/
	// numberOfReadPairs1 used for computing expected and numberOfReadPairs2 as total number

	/*cout << interactions1[4] ;
	cout << interactions2[4] ;
	cout << interactions1[1000] ;
	cout << interactions2[1000] ;
	cout << interactions1[5500] ;
	cout << interactions2[5500] ;
	int v = interactions1.size() - 10;
	int y = v / 2;
	cout << y << " : " << interactions1[y] ;
	cout << y << " : " << interactions2[y] ;
	cout << v << " : " << interactions1[v] ;
	cout << v << " : " << interactions2[v] ;//*/
	//throw std::invalid_argument("");

	parallel_for(size_t(0),size_t(interactions2.size()),
			[&] (size_t i) {
	//for (int i = 0; i < interactions2.size(); i++)
		//{
		Interaction f = interactions1[i];
		Interaction s = interactions2[i];
		if ( f != s )
		{
			cout << "Mismatched Interactions:" << endl << f << s;
		}
		else if (f.getFreq() > 1 && s.getFreq() > 1)
		{
			//cout << f << endl << s << endl;
			double prob = double(f.getFreq())/numberOfReadPairs1;
			//string str = string("prob: ") + to_string(f.getFreq()) + " / " + to_string(numberOfReadPairs1) + " = " + to_string(prob) + "\n";
			//cout << str << prob << endl;

			BinomDataComp I = BinomDataComp(s);
			I.setProbability(prob);
			I.setExpected(f.getFreq());

			binFiltered.push_back(I);
		}
	});

	if (setupValues.getVerbose())
	{
		cout << "\tBinomial Vector Size: " << binFiltered.size() << endl;
	}

	//throw std::invalid_argument("");
	vector<array<double,3>> values;
	values.resize(binFiltered.size());

	cout << "\tcalculating P values" << endl;
	parallel_for(size_t(0),size_t(binFiltered.size()),
			[&] (size_t i) {
		//for (int i = 0; i < binFiltered.size(); i++)
		//{
			int F = binFiltered[i].getFreq();
			double V = binFiltered[i].getProbability();
			double P = binomTest(F, numberOfReadPairs2, V, "two.sided");
			binFiltered[i].setPvalue(P);
			//cout << i << " " << F << " " << numberOfReadPairs2 << " " << V << " " << P << endl;

			//binned_df_filtered2$logFoldChange <- log2((binned_df_filtered2$frequencies/numberOfReadPairs2)/(binned_df_filtered1$frequencies/numberOfReadPairs1))
			double Fd = log2((double(binFiltered[i].getFreq())/numberOfReadPairs2)/(double(binFiltered[i].getExpected())/numberOfReadPairs1));
			//cout << i << " " << binFiltered[i].getFreq() << " " << numberOfReadPairs2;
			//cout << " " << binFiltered[i].getExpected() << " " << numberOfReadPairs1 << " " << Fd<< endl;
			binFiltered[i].setLogObExp(Fd);

			array<double,3> ls = {double(i), P, 0.5};
			values[i] = ls;
		//}
	});
	cout << "\t" << flush;
	completed();


	cout << "\tcalculating Q values" << endl;
	switch(setupValues.getCisTrans())
	{
	case ct_all :
		if(setupValues.getRemoveDiagonal())
		{
			pBhAdjust(values, upperhalfBinNumber);
			//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(), upperhalfBinNumber));
		}
		else
		{
			pBhAdjust(values, upperhalfBinNumber+covS);
			//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(), upperhalfBinNumber+covS));
		}
		break;
	case ct_cis :
		if(setupValues.getRemoveDiagonal())
		{
			pBhAdjust(values, cisBinNumber);
			//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(), cisBinNumber));
		}
		else
		{
			pBhAdjust(values, cisBinNumber+covS);
			//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(),cisBinNumber+covS));
		}
		break;
	case ct_trans:
		pBhAdjust(values, transBinNumber);
		//binFiltered[i].setQvalue(pBhAdjust(binFiltered[i].getProbability(), transBinNumber));
		break;
	}//*/

	parallel_for(size_t(0),size_t(binFiltered.size()),
				[&] (size_t i) {
		//for (int i = 0; i < binFiltered.size(); i++)
	//{
		binFiltered[i].setQvalue(values[i][2]);
	});//*/
	cout << "\t" << flush;
		completed();

	completed();
}

/*
 * .binomialHiChicupComp=function(hicupinteraction1,hicupinteraction2, restrictionFile, sampleName, cistrans='all', parallel=FALSE, cores=8, removeDiagonal=TRUE,baits=NULL,res=NULL)
	{


#all read pairs used in binomial: numberOfReadPairs1 used for computing expected and numberOfReadPairs2 as total number






					 binned_df_filtered2$pvalue <- apply(binned_df_filtered2, 1, function(x)
														{
														binom.test(as.numeric(x[["frequencies"]])-1, numberOfReadPairs2, as.numeric(x[["probabilityOfInteraction"]]), alternative = "two.sided")$p.value
														}
														)

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
