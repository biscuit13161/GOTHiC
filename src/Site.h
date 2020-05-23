/*
 * Site.h
 *
 *  Created on: 18 May 2020
 *      Author: rich
 */

#ifndef SRC_SITE_H_
#define SRC_SITE_H_

#include <string>

class Site{
private:
	std::string mChr;
	int mStart;
	int mEnd;
public:
	Site();
	Site(std::string chr, int start, int end);
	Site(const Site & other);

	friend std::string fixChromosomeNames(std::string chr);

	inline std::string getChr() const {return mChr;}
	inline int getLocus() const {return mStart;}
	inline int getEnd() const {return mEnd;}
	inline int getStart() const {return mStart;}

	const Site & operator=(const Site & other);
	bool operator==(const Site & other);

	void print();
};

#endif /* SRC_SITE_H_ */
