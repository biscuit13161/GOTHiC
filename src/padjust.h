/*
 * padjust.h
 *
 *  Created on: 21 May 2020
 *      Author: rich
 */

#ifndef SRC_PADJUST_H_
#define SRC_PADJUST_H_

#include <vector>
#include <array>

enum pAdjustMethods {
	pam_bh = 0,
	pam_holm,
	pam_hochberg,
	pam_hommel,
	pam_bonferroni,
	pam_by,
	pam_none,
};

double pBhAdjust(double Pi, double n);

void pBhAdjust(std::vector<std::array<double,3>> & p, double n);
bool poscomp(const std::array<double,3> & a, const std::array<double,3> & b);
bool valcomp(const std::array<double,3> & a, const std::array<double,3> & b);


#endif /* SRC_PADJUST_H_ */
