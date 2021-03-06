/*
 *  padjust.cpp
 *
 *  AUTHOR
 *	Richard Thompson, ithompson@hbku.edu.qa
 *	21 May 2020
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
 
 /*
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "padjust.h"
#include <iostream>
#include <map>
#include <algorithm>

using namespace std;

double pBhAdjust(double Pi, double n)
{
    if (n <= 1)
    	return(Pi);

    double Po = n * Pi ;

    if (Po < 1)
    	return Po;

	return 1;
}


bool poscomp(const array<double,3> & a, const array<double,3> & b)
{
	return a[0] < b[0]; // < ascending
}

bool valcomp(const array<double,3> & a, const array<double,3> & b)
{
	return a[1] < b[1]; // > descending
}


void pBhAdjust(vector<array<double,3>> & p, double n) // double,2 -> pos, pval, qval
{
	if (p.size() == 0 || p.empty())
		throw std::invalid_argument("pBhAdjust: input vector empty");

	int x = p.size() - 1;

	sort(p.begin(),p.end(), valcomp); // increasing by pval
	//double minP = 1;
	//double maxP = 0;
	double minQ = 1;

	for (int i = p.size()-1 ; i >= 0; i--)
	{
		double P = p[i][1]; //Pvalue input
		//minP = min(minP, P);
		//maxP = max(maxP, P);
		double a = (n * P) /(i+1);
		minQ = min(minQ,a); //? a : minQ ;
		minQ = min(minQ, 1.0d);
		p[i][2] = minQ;
	}

	sort(p.begin(),p.end(), poscomp); // increasing by pos
}

/*

First, order all p-values from small to large. Then multiply each p-value by the total number
of tests m and divide by its rank order.
Second, make sure that the resulting sequence is non-decreasing: if it ever starts decreasing,
make the preceding p-value equal to the subsequent (repeatedly, until the whole sequence becomes non-decreasing).
If any p-value ends up larger than 1, make it equal to 1.

//*/
