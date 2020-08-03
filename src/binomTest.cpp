/*
 *  binomTest.cpp
 *
 *  AUTHOR
 *	Richard Thompson, ithompson@hbku.edu.qa
 *	May 25, 2020.
 *
 *	Copyright (C) 2020 Richard Thompson, Qatar Biomedical Research Institute
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
 *  https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
 *
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
				int y = 0;
				for (int i = ceil(m); i <= n; i++)
				{
					if (dbinom(i, n, p, 0) <= d * relErr)
					{
						y = n - (i - 1);
						break;
					}
				}//*/

				pval = pbinom(x, n, p, 1, 0) + pbinom(n - y, n, p, 0, 0);
			}
			else
			{
				str = "more\t";
				//i <- seq.int(from = 0, to = floor(m))
				int y = 0;
				for (int i = 0; i <= floor(m); i++)
				{
					if (dbinom(i, n, p, 0) > d * relErr)
					{
						y = i;
						break;
					}
				}//*/

				pval = pbinom(y - 1, n, p, 1, 0) + pbinom(x - 1, n, p, 0, 0);
			}//*/
		}
		break;
	}

	return pval;
}

/*
void preCalcDBinom(int x, int n, double p)
{
	dbinom(i, n, p, 0)
}
//*/
