/*
 * binomTest.cpp
 *
 *  Created on: 25 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "binomTest.h"
#include "pbinom.h"
#include "dbinom.h"
#include <map>
#include <iostream>

using namespace std;

double binomTest(int x, int n, double p, string alternative)
{
	string str = "";
	double val = binomTest( x, n, p, alternative, str);
	return val;
}

double binomTest(int x, int n, double p, string alternative, string& str)
{
	// x = frequency
	// n = number of pairs
	// p = expected probabilityOfInteraction
    if ( x < 0)
    	throw std::invalid_argument("binomTest: 'x' must be nonnegative and integer");

    if ( p < 0 || p > 1)
    	throw std::invalid_argument("binomTest: 'p' must be a single number between 0 and 1");

    map<string, altOptions> optionValues;
    optionValues["less"] = ao_less;
    optionValues["greater"] = ao_greater;
    optionValues["two.sided"] = ao_two_sided;

    double pval = 0;

    switch(optionValues[alternative])
    {
    case ao_less:
    	pval = pbinom(x, n, p, 1, 0);
    	break;
    case ao_greater:
    	pval = pbinom(x - 1, n, p, 0, 0);
    	break;
    case ao_two_sided:
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
            double d = dbinom(x, n, p, 0);
            double m = n * p;
            if (x == m)
            {
            	str = "equals";
            	pval = 1;
            }
            else if (x < m)
            {
            	str = "less\t";
            	//i <- int(from = ceiling(m), to = n)
            	int y = 0;
            	for (int i = ceil(m); i <= n; i++)
            		if (dbinom(i, n, p, 0) <= d * relErr)
            			y++;
            	pval = pbinom(x, n, p, 1, 0) + pbinom(n - y, n, p, 0, 0);
            }
            else
            {
            	str = "more\t";
            	//i <- seq.int(from = 0, to = floor(m))
            	int y = 0;
            	for (int i = 0; i <= floor(m); i++)
            		if (dbinom(i, n, p, 0) <= d * relErr)
            			{
            			y++;
            			}
            	pval = pbinom(y - 1, n, p, 1, 0) + pbinom(x - 1, n, p, 0, 0);
            }//*/
        }
        break;
    }

    return pval;
}
