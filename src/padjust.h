/*
 * padjust.h
 *
 *  Created on: 21 May 2020
 *      Author: rich
 */

#ifndef SRC_PADJUST_H_
#define SRC_PADJUST_H_

#include <vector>

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



#endif /* SRC_PADJUST_H_ */
