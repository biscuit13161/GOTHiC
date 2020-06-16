/*
 * Baits.h
 *
 *  Created on: 15 Jun 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#ifndef SRC_BAITS_H_
#define SRC_BAITS_H_

#include <string>
#include <vector>

class Bait {
private:
	std::string mChr;
	int mStart;
	int mEnd;
	std::string mBait;

public:
	Bait();
	Bait(std::string mChr, int mStart,	int mEnd, std::string mBait);
	Bait(const Bait & other );

	inline std::string getChr() const {return mChr;}
	inline std::string getBait() const {return mBait;}
	inline int getStart() const {return mStart;}
	inline int getEnd() const {return mEnd;}

	friend bool comp(const Bait & a, const Bait & b);
	friend std::ostream & operator<<(std::ostream & out, const Bait & in);
};

class Locus {
private:
	std::string mChr;
	int mStart;

public:
	Locus();
	Locus(std::string mChr, int mStart);
	Locus(const Locus & other);

	inline std::string getChr() const {return mChr;}
	inline int getStart() const {return mStart;}
	inline int getEnd(int res) const {return mStart + res;}
};

struct Comp
    {
        bool operator() ( const Bait & s, std::string i ) const { return s.getChr() < i; }
        bool operator() ( std::string i, const Bait & s ) const { return i < s.getChr(); }
    };

void readBaits(std::vector<Bait> & Baits, std::string fileName);
bool comp(const Bait & a, const Bait & b);
std::ostream & operator<<(std::ostream & out, const Bait & in);

std::string findOverlaps(std::vector<Bait> & Baits, Locus & L, int res);

#endif /* SRC_BAITS_H_ */
