/*
 * BiinomData.h
 *
 *  Created on: 5 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#ifndef SOURCE_BINOMDATA_H_
#define SOURCE_BINOMDATA_H_

#include <string>
#include <math.h>

enum CisTrans
{
	ct_all=0,
	ct_cis,
	ct_trans,
};

class halfInteraction
{
private:
	std::string mChr;
	int mLocus;
	std::string mInt;


public:
	halfInteraction();
	halfInteraction(std::string chr, int locus);
	halfInteraction(const halfInteraction & other);
	friend std::string fixChromosomeNames(std::string chr);

	inline std::string getChr() const {return mChr;}
	inline int getLocus() const {return mLocus;}
	inline std::string getInt() const {return mInt;}

	void print();
	friend bool comp(const halfInteraction & a, const halfInteraction & b);

};

class Interaction
{
private:
	std::string mChr1;
	std::string mChr2;
	//chromosome(s) containing interacting regions 1 and 2

	int mLocus1;
	int mLocus2;
	//start positions of the interacting regions 1 and 2 in the corresponding chromosome(s)

	std::string mInt1;
	std::string mInt2;
	//combined strings for comparisons

	int mFrequency;

public:
	Interaction();
	Interaction(std::string chr1, std::string chr2, int locus1,	int locus2);
	Interaction(std::string chr1, std::string chr2, int locus1,	int locus2, int freq);
	Interaction(const Interaction & other);
	Interaction(const halfInteraction & first, const halfInteraction & second);
	Interaction(std::string int1, std::string int2, int freq);
	friend std::string fixChromosomeNames(std::string chr);

	inline std::string getChr1() const {return mChr1;}
	inline int getLocus1() const {return mLocus1;}
	inline int getLocus2() const {return mLocus2;}
	inline std::string getChr2() const {return mChr2;}
	inline std::string getInt1() const {return mInt1;}
	inline std::string getInt2() const {return mInt2;}

	inline void setLocus1(int L) {mLocus1 = L; mInt1 = mChr1 + ":" + std::to_string(mLocus1);}
	inline void setLocus2(int L) {mLocus2 = L; mInt2 = mChr2 + ":" + std::to_string(mLocus2);}
	void print();

};



class BinomData
{
private:
	std::string mChr1;
	std::string mChr2;
	//chromosome(s) containing interacting regions 1 and 2

	int mLocus1;
	int mLocus2;
	//start positions of the interacting regions 1 and 2 in the corresponding chromosome(s)


	double mRelCoverage1;
	double mRelCoverage2;
	// relative coverage corresponding to regions 1 and 2

	double mProbability; // expected frequency
	int mExpected; // expected number of reads
	int mReadCount; // observed reads number
	double mPvalue; // binomial p-value
	double mQvalue; // binomial p-value corrected for multi-testing with Benjamini-Hochberg
	double mLogObservedOverExpected; // observed/expected read numbers log ratio

public:
	BinomData();
	BinomData(std::string chr1, std::string chr2, int locus1, int locus2, double relCoverage1, double relCoverage2, double probability, int expected, int readCount, double pvalue, double qvalue, double logObservedOverExpected);
	BinomData(BinomData & other);

};

std::string fixChromosomeNames(std::string chr);
bool comp(const halfInteraction & a, const halfInteraction & b);

//double relCoverage1, double relCoverage2, double probability, int mExpected, int mReadCount, double mPvalue, double mQvalue, double mLogObservedOverExpected

#endif /* SOURCE_BINOMDATA_H_ */
