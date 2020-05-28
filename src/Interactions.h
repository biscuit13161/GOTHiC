/*
 * Interactions.h
 *
 *  Created on: 27 May 2020
 *      Author: rich
 */

#ifndef SRC_INTERACTIONS_H_
#define SRC_INTERACTIONS_H_

#include "Utils.h"
#include <string>
#include <vector>

class halfInteraction
{
private:
	std::string mChr;
	int mLocus;

public:
	halfInteraction();
	halfInteraction(std::string chr, int locus);
	halfInteraction(const halfInteraction & other);
	friend std::string fixChromosomeNames(std::string chr);

	inline std::string getChr() const {return mChr;}
	inline int getLocus() const {return mLocus;}
	inline std::string getInt() const {return mChr + ":" + std::to_string(mLocus);}

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

	int mFrequency;

public:
	Interaction();
	Interaction(std::string chr1, std::string chr2, int locus1,	int locus2);
	Interaction(std::string chr1, std::string chr2, int locus1,	int locus2, int freq);
	Interaction(const Interaction & other);
	Interaction(const halfInteraction & first, const halfInteraction & second);
	friend std::string fixChromosomeNames(std::string chr);

	inline std::string getChr1() const {return mChr1;}
	inline std::string getChr2() const {return mChr2;}
	inline int getLocus1() const {return mLocus1;}
	inline int getLocus2() const {return mLocus2;}
	inline std::string getInt1() const {return mChr1 + ":" + std::to_string(mLocus1);}
	inline std::string getInt2() const {return mChr2 + ":" + std::to_string(mLocus2);}
	inline int getFreq() const {return mFrequency;}

	inline void setLocus1(int L) {mLocus1 = L;}
	inline void setLocus2(int L) {mLocus2 = L;}
	void print();

	bool operator==(const Interaction & other);
	friend std::ostream & operator<<(std::ostream & out, const Interaction & in);
	friend bool intcomp(const Interaction & a, const Interaction & b);

};

bool intcomp(const Interaction & a, const Interaction & b);
void checkInteractions(const std::vector<Interaction> & interactions, const std::vector<Interaction> & interactions2);

std::ostream & operator<<(std::ostream & out, const Interaction & in);

void writeBinary(std::vector<Interaction> & interactions, std::string binOutFileName);
void readBinary(std::vector<Interaction> & interactions, std::string binInFileName);

bool comp(const halfInteraction & a, const halfInteraction & b);

#endif /* SRC_INTERACTIONS_H_ */
