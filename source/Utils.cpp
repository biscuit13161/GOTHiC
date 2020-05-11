/*
 * Utils.cpp
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#include "Utils.h"
#include <ctime>
#include <iostream>

using namespace std;

void showTime()
{
	// current date/time based on current system
	time_t now = time(0);

	// convert now to string form
	//char* dt = ctime(&now);
	tm *ltm = localtime(&now);

	cerr << 1 + ltm->tm_hour << ":" << flush;
	cerr << 1 + ltm->tm_min << ":" << flush;
	cerr << 1 + ltm->tm_sec << flush;
}

void completed()
{

	cerr << "\t... Completed: " << flush;
	showTime();
	cerr << endl;
}
