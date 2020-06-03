/*
 * gamma.h
 *
 *  Created on: 2 Jun 2020
 *      Author: rich
 */

#ifndef SRC_GAMMA_H_
#define SRC_GAMMA_H_

#ifdef __cplusplus
extern "C" {
#endif //*/

double gammafn(double x);
double chebyshev_eval(double x, const double *a, const int n);
double fmax2(double x, double y);

#ifdef __cplusplus
}
#endif//*/



#endif /* SRC_GAMMA_H_ */
