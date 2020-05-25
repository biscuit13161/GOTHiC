/*
 * Utils.cpp
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#include "Setup.h"
#include "version.h"
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

void Setup::print()
{
	cerr << endl << string("GOTHiC++ v") << GOTH_MAJOR_VERSION << "." << GOTH_MINOR_VERSION << "." << GOTH_PATCH_VERSION << endl << endl;
	cerr << "#" << endl;
	cerr << "# Output Directory: " << mOutDir << endl;
	cerr << "# Restriction File: " << mEnzyme << endl;
	cerr << "# Sample Name:      " << mSname << endl;
	cerr << "# Input Data File:  " << mInput << endl;
	cerr << "# Threads:          " << mThreads << endl;
	cerr << "# Resolution:       " << mRes << endl;
	cerr << "# Config File:      " << mCliName << endl;
	cerr << "#" << endl << endl;
}

Setup loadConfig(string & fileName)
{
	ifstream inFile;
	inFile.open(fileName);

	if (!inFile.is_open())
	{
		throw std::invalid_argument("Could not load Config file!");
	}

	Setup setupValues;
	setupValues.setCliName(fileName);

	while (inFile)
	{
		string line;
		getline(inFile,line,'\n');
		if (line.find('#') != 0 && line.length() > 0)
		{
			size_t pos = line.find(":");
			string id = line.substr(0,pos);
			string value = line.substr(pos+1).c_str();
				while (value.find(' ') == 0)
				{
					value = value.substr(1);
				}

			//cout << id << "_" << value << endl;
			std::map<std::string, Options> optionValues;
			optionValues["Input"] = Input;
			optionValues["SampleName"] = Sname;
			optionValues["Digest"] = Digest;
			optionValues["Threads"] = Threads;
			optionValues["Res"] = Res;
			optionValues["Output"] = Output;
			optionValues["CisTrans"] = Cistrans;
			
			std::map<std::string,CisTrans> CToptionValues;
			CToptionValues["all"] = ct_all;
			CToptionValues["cis"] = ct_cis;
			CToptionValues["trans"] = ct_trans;

			switch(optionValues[id])
			{
			case Input:
				setupValues.setInput(value);
				break;
			case Sname:
				setupValues.setSname(value);
				break;
			case Digest:
				setupValues.setEnzyme(value);
				break;
			case Threads:
				setupValues.setThreads(atoi(value.c_str()));
				break;
			case Res:
				setupValues.setRes(atoi(value.c_str()));
				break;
			case Output:
				if (value == "")
				{
					//auto cwd = boost::filesystem::current_path();
					setupValues.setOutDir("./");
				}
				else
				{
					setupValues.setOutDir(value);
				}
				break;
			case Cistrans:
				setupValues.setCisTrans(CToptionValues[value]);
				break;
			}//*/
		}
	}
	inFile.close();

	return setupValues;
}

