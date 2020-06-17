/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2015  The R Core Team
 *  Copyright (C) 2004-2015  The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  DESCRIPTION
 *
 *    The distribution function of the binomial distribution.
 *
 *  MODIFIED
 *  	Richard Thompson, ithompson@hbku.edu.qa
 *  	June 2, 2020
 */

#include "pbinom.h"
#include <math.h>
#include "dpq.h"

#define ML_NEGINF	((-1.0) / 0.0)

#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)


double pbinom(double x, double n, double p, int lower_tail, int log_p)
{
	/*
	 * x = vector of quantiles
	 * n = number of trials
	 * p = probability of success on each trial.
	 * lowertail = logical (0/1); if TRUE (default), probabilities are P[X ≤ x], otherwise, P[X > x].
	 * log_p = return p-value as log value (0/1)
	 */

	if(n < 0 || p < 0 || p > 1)
		return ML_NAN;

    if (x < 0) return R_DT_0;
    x = floor(x + 1e-7);
    if (n <= x) return R_DT_1;
    return pbeta(p, x + 1, n - x, !lower_tail, log_p);
}


double pbeta(double x, double a, double b, int lower_tail, int log_p)
{
    // treat limit cases correctly here:
	if(a == 0 || b == 0 ) {
		if(a == 0 && b == 0) // point mass 1/2 at each of {0,1} :
			return (log_p ? -M_LN2 : 0.5);
		if (a == 0 || a/b == 0) // point mass 1 at 0 ==> P(X <= x) = 1, all x > 0
			return R_DT_1;
		if (b == 0 || b/a == 0) // point mass 1 at 1 ==> P(X <= x) = 0, all x < 1
			return R_DT_0;
		// else, remaining case:  a = b = Inf : point mass 1 at 1/2
		if (x < 0.5) return R_DT_0; else return R_DT_1;
	}
    // Now:  0 < a < Inf;  0 < b < Inf

    double x1 = 0.5 - x + 0.5, w, wc;
    int ierr;
    //====
    bratio(a, b, x, x1, &w, &wc, &ierr, log_p); /* -> ./toms708.c */
    return lower_tail ? w : wc;
} /* pbeta_raw() */
