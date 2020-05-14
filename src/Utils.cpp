/*
 * Utils.cpp
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#include "Utils.h"
#include "version.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

using namespace std;

void showTime()
{
	// current date/time based on current system
	time_t now = time(0);

	// convert now to string form
	//char* dt = ctime(&now);
	tm *ltm = localtime(&now);

	cerr << ltm->tm_hour << ":" << flush;
	cerr << ltm->tm_min << ":" << flush;
	cerr << ltm->tm_sec << flush;
}

void completed()
{

	cerr << "\t... Completed! " << flush;
	//showTime();
	cerr << endl;
}

void printUsage()
{
	cerr << "GOTHiC++ usage:" << endl << endl;
	cerr << "    gothic <path/to/config/file>"  << endl << endl;
	cerr << string("version: ") << GOTH_MAJOR_VERSION << "." << GOTH_MINOR_VERSION << "." << GOTH_PATCH_VERSION << endl;

}
