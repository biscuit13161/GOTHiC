/*
 * toms708.h
 *
 *  Created on: 17 May 2020
 *      Author: rich
 */

#ifndef SRC_TOMS708_H_
#define SRC_TOMS708_H_

#include "dpq.h"
#include <stdbool.h>

static void bratio(double a, double b, double x, double y, double *w, double *w1, int *ierr, int log_p);
static double bfrac(double, double, double, double, double, double, int log_p);
static void bgrat(double, double, double, double, double *, double, int *, bool log_w);
static double grat_r(double a, double x, double r, double eps);
static double apser(double, double, double, double);
static double bpser(double, double, double, double, int log_p);
static double basym(double, double, double, double, int log_p);
static double fpser(double, double, double, double, int log_p);
static double bup(double, double, double, double, int, double, int give_log);
static double exparg(int);
static double psi(double);
static double gam1(double);
static double gamln1(double);
static double betaln(double, double);
static double algdiv(double, double);
static double brcmp1(int, double, double, double, double, int give_log);
static double brcomp(double, double, double, double, int log_p);
static double rlog1(double);
static double bcorr(double, double);
static double gamln(double);
static double alnrel(double);
static double esum(int, double, int give_log);
static double erf__(double);
static double rexpm1(double);
static double erfc1(int, double);
static double gsumln(double, double);




#endif /* SRC_TOMS708_H_ */
