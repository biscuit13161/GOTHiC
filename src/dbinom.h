/*
 * dbinom.h
 *
 *  Created on: 2 Jun 2020
 *      Author: rich
 */

#ifndef SRC_DBINOM_H_
#define SRC_DBINOM_H_

#ifdef __cplusplus
extern "C" {
#endif //*/

double dbinom_raw(double x, double n, double p, double q, int give_log);
double dbinom(double x, double n, double p, int give_log);

#ifdef __cplusplus
}
#endif//*/

#endif /* SRC_DBINOM_H_ */
