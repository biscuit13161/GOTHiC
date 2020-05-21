/*
 * padjust.cpp
 *
 *  Created on: 21 May 2020
 *      Author: rich
 */

#include "padjust.h"
#include <map>

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

