/*
 * binomTest.cpp
 *
 *  Created on: 25 May 2020
 *      Author: rich
 */

#include "binomTest.h"
#include "pbinom.h"
#include <map>
#include <iostream>


using namespace std;

double binomTest(int x, int n, double p, string alternative)
{
    if ( x < 0)
    	throw std::invalid_argument("binomTest: 'x' must be nonnegative and integer");

	/*if (length(x) == 2L) {
        n <- sum(x)
        x <- x[1L]
    }
    else if (length(x) == 1L) {
        nr <- round(n)
        if ((length(n) > 1L) || is.na(n) || (n < 1) || abs(n -
            nr) > 1e-07 || (x > nr))
            stop("'n' must be a positive integer >= 'x'")
        DNAME <- paste(DNAME, "and", deparse1(substitute(n)))
        n <- nr
    }
    else stop("incorrect length of 'x'")//*/

    if ( p < 0 || p > 1)
    	throw std::invalid_argument("binomTest: 'p' must be a single number between 0 and 1");

    /*if (!((length(conf.level) == 1L) && is.finite(conf.level) &&
        (conf.level > 0) && (conf.level < 1)))
        stop("'conf.level' must be a single number between 0 and 1")//*/


		map<string, altOptions> optionValues;
		optionValues["less"] = ao_less;
		optionValues["greater"] = ao_greater;
		optionValues["two.sided"] = ao_two_sided;

	double pval = 0;

    switch(optionValues[alternative])
    {
    case ao_less:
    	pval = pbinom(x, n, p, true, false);
    	break;
    case ao_greater:
    	pval = pbinom(x - 1, n, p, false,false);
    	break;
    case ao_two_sided:
    	/** needs work **/
    	cerr << "*** BEWARE: tow.sided binomTest incomplete ***" << endl;
    	if (p == 0)
        	{
        	pval = (x == 0);
        	}
        else if (p == 1)
        	{
        	pval = (x == n);
        	}
        else
        {
            double relErr = 1 + 1e-07;
            //d <- dbinom(x, n, p,true, false)
          /*  m <- n * p;
            if (x == m)
            {
            	pval = 1;
            }
            else if (x < m)
            {
                i <- seq.int(from = ceiling(m), to = n)
                y <- sum(dbinom(i, n, p,true, false) <= d * relErr)
                pval = pbinom(x, n, p,true, false) + pbinom(n - y, n, p,  false,false);
            }
            else
            {
                i <- seq.int(from = 0, to = floor(m))
                y <- sum(dbinom(i, n, p,true, false) <= d * relErr)
                pval = pbinom(y - 1, n, p,true, false) + pbinom(x - 1, n, p, false,false);
            }//*/
        }
        break;
    }

    return pval;
    /*p.L <- function(x, alpha) {
        if (x == 0)
            0
        else qbeta(alpha, x, n - x + 1)
    }
    p.U <- function(x, alpha) {
        if (x == n)
            1
        else qbeta(1 - alpha, x + 1, n - x)
    }
    CINT <- switch(alternative, less = c(0, p.U(x, 1 - conf.level)),
        greater = c(p.L(x, 1 - conf.level), 1), two.sided = {
            alpha <- (1 - conf.level)/2
            c(p.L(x, alpha), p.U(x, alpha))
        })
    attr(CINT, "conf.level") <- conf.level
    ESTIMATE <- x/n
    names(x) <- "number of successes"
    names(n) <- "number of trials"
    names(ESTIMATE) <- names(p) <- "probability of success"
    structure(list(statistic = x, parameter = n, p.value = PVAL,
        conf.int = CINT, estimate = ESTIMATE, null.value = p,
        alternative = alternative, method = "Exact binomial test",
        data.name = DNAME), class = "htest")//*/
}
