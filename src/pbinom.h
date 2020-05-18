/*
 * pbinom.h
 *
 *  Created on: 17 May 2020
 *      Author: rich
 */

#ifndef SRC_PBINOM_H_
#define SRC_PBINOM_H_

#include "toms708.h"

#ifdef __cplusplus
extern "C" {
#endif //*/

double pbinom(double x, double n, double p, int lower_tail, int log_p);
double pbeta(double x, double a, double b, int lower_tail, int log_p);

#ifdef __cplusplus
}
#endif//*/

#endif /* SRC_PBINOM_H_ */
