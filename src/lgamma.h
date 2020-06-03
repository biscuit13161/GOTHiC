/*
 * lgamma.h
 *
 *  Created on: 2 Jun 2020
 *      Author: rich
 */

#ifndef SRC_LGAMMA_H_
#define SRC_LGAMMA_H_

#ifdef __cplusplus
extern "C" {
#endif //*/

double lgammafn_sign(double x, int *sgn);
double lgammafn(double x);
double lgammacor(double x);
double sinpi(double x);

#ifdef __cplusplus
}
#endif//*/

#endif /* SRC_LGAMMA_H_ */
