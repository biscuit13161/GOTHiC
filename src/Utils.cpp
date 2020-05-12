/*
 * Utils.cpp
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#include "../src/Utils.h"

#include <ctime>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

using namespace std;

Setup::Setup(): mOutDir(""), mEnzyme(""), mInput(""), mThreads(0)
{

}

Setup::Setup(string outDir, string enzyme, string input, int threads): mOutDir(outDir), mEnzyme(enzyme), mInput(input), mThreads(threads)
{

}

Setup loadConfig(fileName)
{
	ifstream inFile;
	inFile.open(fileName);

	if !(inFile.is_open())
	{
		throw std::invalid_argument("Could not load Config file!");
	}

	Setup setupValues;

	while (inFile)
	{
		string id
		getline(inFile,id1,': ');
		string value;
		getline(inFile,value,'\n');
		switch(id)
		{
		case 'Input':
			setupValues.setInput(value);
			break;
		case 'Digest':
			setupValues.setDigest(value);
			break;
		case 'Threads':
			setupValues.setThreads(atoi(value.c_str()));
			break;
		case 'Res':
			setupValues.setRes(atoi(value.c_str()));
			break;
		case 'Output':
			if (value == "")
			{
				auto cwd = boost::filesystem::current_path();
				setupValues.setOutput(cwd);
			}
			else
			{
				setupValues.setOutput(value);
			}
			break;
		}
	}
	inFile.close();

	return setupValues;
}

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

void printUsage()
{
	cerr << 'Usage here' << endl;
}
