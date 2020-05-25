/*
 * padjust.cpp
 *
 *  Created on: 21 May 2020
 *      Author: rich
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
	return a[0] > b[0];
}

bool valcomp(const array<double,3> & a, const array<double,3> & b)
{
	return a[1] < b[1];
}


void pBhAdjust(vector<array<double,3>> & p, double n) // int,2 -> pos, pval, qval
{
	if (p.size() == 0 || p.empty())
		throw std::invalid_argument("pBhAdjust: input vector empty");

	int x = p.size() - 1;

	cout << p.size() << endl;

    cout << p[0][0] << " - " << p[0][1] << " - " << p[0][2] << endl;
    cout << p[x][0] << " - " << p[x][1] << " - " << p[x][2] << endl;
        sort(p.begin(),p.end(), valcomp); // decreasing by pval
        double minP = 1;
        double maxP = 0;
        double minQ = 1;
        cout << p[0][0] << " - " << p[0][1] << " - " << p[0][2] << endl;
        cout << p[x][0] << " - " << p[x][1] << " - " << p[x][2] << endl;

		for (int i = 0; i < p.size(); i++)
		{
			double P = p[i][1];
			minP = min(minP, P);
			maxP = max(maxP, P);
			double a = (n * P) /i;
			//cout  << "(" << p[i-1][0] << ", " << p[i-1][1] << ", " << a << ", "<< flush;
			minQ = min(minQ,a); //? a : minQ ;
			//cout  << minP << ")  " << flush;
			p[i][2] = minQ;
		}
	    cout << "Minimum Pvalue: " << minP << endl;
	    cout << "Maximum Pvalue: " << maxP << endl;
		cout << "Minimum Qvalue:" << minQ << endl;

        sort(p.begin(),p.end(), poscomp); // increasing by pos
}

/*

First, order all p-values from small to large. Then multiply each p-value by the total number
of tests m and divide by its rank order.
Second, make sure that the resulting sequence is non-decreasing: if it ever starts decreasing,
make the preceding p-value equal to the subsequent (repeatedly, until the whole sequence becomes non-decreasing).
If any p-value ends up larger than 1, make it equal to 1.

//*/
