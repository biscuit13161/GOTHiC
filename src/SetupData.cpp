/*
 *  Utils.cpp
 *
 *  AUTHOR
 *	Richard Thompson, ithompson@hbku.edu.qa
 *	11 May 2020
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
 

#include "SetupData.h"

#include "version.h"
#include <getopt.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

using namespace std;

SetupData::SetupData(): mOutDir("./"), mEnzyme(""), mInput(""), mThreads(1), mRes(10000), mRemoveDiagonal(true), mTime(false)
{

}

SetupData::SetupData(string outDir, string enzyme, string input, int threads): mOutDir(outDir), mEnzyme(enzyme), mInput(input), mThreads(threads), mTime(false)
{

}

void SetupData::print()
{
	cerr << endl << string("GOTHiC++ v") << GOTH_MAJOR_VERSION << "." << GOTH_MINOR_VERSION << "." << GOTH_PATCH_VERSION << endl << endl;
	cerr << "#" << endl;
	cerr << "# Output Directory: " << mOutDir << endl;
	cerr << "# Restriction File: " << mEnzyme << endl;
	cerr << "# Sample Name:      " << mSname << endl;
	cerr << "# Input Data File:  " << mInput << endl;
	cerr << "# Threads:          " << mThreads << endl;
	fprintf(stderr, "# Resolution:       %d\n", mRes );
	cerr << "# Config File:      " << mCliName << endl;
	cerr << "# Verbose:          " << mVerbose << endl;
	cerr << "# Analysis type:    ";
	if (mAnalysisType)
		cerr << "Comparative" << endl;
	else
		cerr << "Single" << endl;
	if (mTime)
		cerr << "# Output Time and Memory Resources" << endl;
	cerr << "#" << endl << endl;
}

SetupData loadConfig(string & fileName)
{
	ifstream inFile;
	inFile.open(fileName);

	if (!inFile.is_open())
	{
		throw std::invalid_argument("Could not load Config file!");
	}

	SetupData setupValues;
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
			optionValues["Analysis"] = Analysis;
			optionValues["RemoveDiagonals"] = RemDiag;
			optionValues["Logfile"] = Logfile;
			optionValues["Verbose"] = Verbose;
			
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
					value = removeSpaces(value);
					if (value == "")
						value = "./";
					else if (! boost::algorithm::ends_with(value, "/"))
						value += "/";
					setupValues.setOutDir(value);
				}
				break;
			case Cistrans:
				setupValues.setCisTrans(value);
				break;
			case Analysis:
				setupValues.setAnalysisType(value);
				break;
			case RemDiag:
				std::transform(value.begin(), value.end(), value.begin(),
				    [](unsigned char c){ return std::tolower(c); });
				if (value == "false")
					setupValues.setRemoveDiagonal(false);
				else
					setupValues.setRemoveDiagonal(true);
				break;
			case Logfile:
				setupValues.setLogFile(value);
				break;
			case Verbose:
				std::transform(value.begin(), value.end(), value.begin(),
						[](unsigned char c){ return std::tolower(c); });
				if (value == "true")
					setupValues.setVerbose(true);
				break;
			}//*/
		}
	}
	inFile.close();

	return setupValues;
}

SetupData setConfig(int argc, char * argv[])
{
	SetupData setupValues;

	static int verbose_flag;
	static int remdiag_flag;
	static int time_flag;

	static struct option long_options[] =
	{
			/* These options set a flag. */
			{"verbose", no_argument,       &verbose_flag, 1},
			{"no_rem_diag",   no_argument,       &remdiag_flag, 1},
			{"time",   no_argument,       &time_flag, 1},
			/* These options don’t set a flag.
			             We distinguish them by their indices. */
			{"input",     required_argument,       0, 'i'},
			{"sample",  required_argument,       0, 's'},
			{"digest",  required_argument, 0, 'd'},
			{"threads",  required_argument, 0, 't'},
			{"res",    required_argument, 0, 'r'},
			{"output",     no_argument,       0, 'o'},
			{"cistrans",  no_argument,       0, 'c'},
			{"analysis",  required_argument, 0, 'A'},
			{"log",    required_argument, 0, 'l'},
			{0, 0, 0, 0}
	};

	int option_index = 0;
	int opt;

	while ((opt = getopt_long(argc, argv, "i:s:d:t:r:o:c:A:l:", long_options, &option_index)) != -1)
	{

		switch(opt)
		{
		case 'i':
			setupValues.setInput(optarg);
			break;
		case 's':
			setupValues.setSname(optarg);
			break;
		case 'd':
			setupValues.setEnzyme(optarg);
			break;
		case 't':
			setupValues.setThreads(atoi(optarg));
			break;
		case 'r':
			setupValues.setRes(atoi(optarg));
			break;
		case 'o':
			setupValues.setOutDir(optarg);
			break;
		case 'c':
			setupValues.setCisTrans(optarg);
			break;
		case 'A':
			setupValues.setAnalysisType(optarg);
			break;
		case 'l':
			setupValues.setLogFile(optarg);
			break;
		}//*/
	}

	if (verbose_flag)
		setupValues.setVerbose(true);

	if (remdiag_flag)
		setupValues.setRemoveDiagonal(true);
	else
		setupValues.setRemoveDiagonal(false);

	if (time_flag)
		setupValues.setTime(true);

	return setupValues;
}

void SetupData::setAnalysisType(string L)
{

	if (string("single").find(L) != string::npos)
	{
		mAnalysisType = ao_single;
	}
	else if (string("comparative").find(L) != string::npos)
	{
		mAnalysisType = ao_comparative;
	}
	else
	{
		throw std::invalid_argument("Unknown Analysis type requested\n\tPlease select 'single' or 'comparative'");
	}

}

void SetupData::setCisTrans(string input)
{
	std::map<std::string,CisTrans> CToptionValues;
	CToptionValues["all"] = ct_all;
	CToptionValues["cis"] = ct_cis;
	CToptionValues["trans"] = ct_trans;

	mCisTrans = CToptionValues[input];
}

