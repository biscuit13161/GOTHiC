/*
 *  SetupComp.cpp
 *
 *  AUTHOR
 *	Richard Thompson, ithompson@hbku.edu.qa
 *	27 May 2020
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
 

#include "SetupComp.h"
#include "version.h"
#include <getopt.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

using namespace std;

SetupComp::SetupComp(): mOutDir(""), mEnzyme(""), mCondition1(""), mCondition2(""), mThreads(0), mRes(10000), mAlpha("0.1"), mQvalue(qv_ihw), mRemoveDiagonal(true), mRandom(true), mVerbose(vl_none)
{

}

void SetupComp::print()
{
	cout << endl << string("GOTHiC++ v") << GOTH_MAJOR_VERSION << "." << GOTH_MINOR_VERSION << "." << GOTH_PATCH_VERSION << endl << endl;
	cout << "#" << endl;
	cout << "# Output Directory:  " << mOutDir << endl;
	cout << "# Restriction File:  " << mEnzyme << endl;
	cout << "# Baits File:        " << mBaits << endl;
	cout << "# Sample Name:       " << mSname << endl;
	cout << "# Control File:      " << mCondition1 << endl;
	cout << "# Sample File:       " << mCondition2 << endl;
	cout << "# Threads:           " << mThreads << endl;
	cout << "# Resolution:        " << mRes << endl;
	cout << "# Config File:       " << mCliName << endl;
	switch(mVerbose)
	{
		case vl_info:
			cout << "# Verbose Output:    Verbose"<< endl;
			break;
		case vl_debug:
			cout << "# Verbose Output:    Debug"<< endl;
			break;
		default:
			cout << "# Verbose Output:    None"<< endl;
	}
	if (mRandom)
		cout << "# Random Subsetting: True" << endl;
	else
		cout << "# Random Subsetting: False"  << endl;
	if (mQvalue == qv_ihw)
		cout << "# Pvalue Correction: Independent Hypothesis Weighting" <<endl;
	else if (mQvalue == qv_bh)
		cout << "# Pvalue Correction: Benjamini Hochberg" <<endl;
	else
		cout << "# Pvalue Correction: Not Specfied!" <<endl;
	cout << "# Alpha:             " << mAlpha << endl;
	if (mTime)
		cerr << "# Output Time and Memory Resources" << endl;
	cout << "#" << endl << endl;
}

SetupComp loadConfigComp(string & fileName)
{
	ifstream inFile;
	inFile.open(fileName);

	if (!inFile.is_open())
	{
		throw std::invalid_argument("Could not load Config file!");
	}

	SetupComp setupValues;
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
			std::map<std::string, SC_Options> optionValues;
			optionValues["Control"] = sc_Cond1;
			optionValues["Sample"] = sc_Cond2;
			optionValues["SampleName"] = sc_Sname;
			optionValues["Digest"] = sc_Digest;
			optionValues["Baits"] = sc_Baits;
			optionValues["Threads"] = sc_Threads;
			optionValues["Res"] = sc_Res;
			optionValues["Output"] = sc_Output;
			optionValues["CisTrans"] = sc_Cistrans;
			optionValues["RemoveDiagonals"] = sc_RemDiag;
			optionValues["Verbose"] = sc_Verbose;
			optionValues["RandomSubset"] = sc_Random;
			optionValues["Alpha"] = sc_Alpha;
			optionValues["Algorithm"] = sc_Qvalues;

			std::map<std::string,CisTrans> CToptionValues;
			CToptionValues["all"] = ct_all;
			CToptionValues["cis"] = ct_cis;
			CToptionValues["trans"] = ct_trans;


			switch(optionValues[id])
			{
			case sc_Cond1:
				setupValues.setCondition1(value);
				break;
			case sc_Cond2:
				setupValues.setCondition2(value);
				break;
			case sc_Sname:
				setupValues.setSname(value);
				break;
			case sc_Digest:
				setupValues.setEnzyme(value);
				break;
			case sc_Baits:
				setupValues.setBaits(value);
				break;
			case sc_Threads:
				setupValues.setThreads(atoi(value.c_str()));
				break;
			case sc_Res:
				setupValues.setRes(atoi(value.c_str()));
				break;
			case sc_Output:
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
			case sc_Cistrans:
				if (value.empty())
					setupValues.setCisTrans("all");
				else
					setupValues.setCisTrans(value);
				break;
			case sc_RemDiag:
				std::transform(value.begin(), value.end(), value.begin(),
				    [](unsigned char c){ return std::tolower(c); });
				if (value == "false")
					setupValues.setRemoveDiagonal(false);
				break;
			case sc_Logfile:
				setupValues.setLogFile(value);
				break;
			case sc_Verbose:
				std::transform(value.begin(), value.end(), value.begin(),
						[](unsigned char c){ return std::tolower(c); });
				if (value == "true")
					setupValues.setVerbose(vl_info);
				else if (value == "debug")
				{
					setupValues.setVerbose(vl_debug);
				}
				break;
			case sc_Random:
				std::transform(value.begin(), value.end(), value.begin(),
						[](unsigned char c){ return std::tolower(c); });
				if (value == "false")
					setupValues.setRandom(false);
				break;
			case sc_Alpha:
				if (value.empty())
					setupValues.setAlpha(string("0.1"));
				else
					setupValues.setAlpha(value);
				break;
			case sc_Qvalues:
				setupValues.setQvalue(value);
				break;
			}//*/
		}
	}
	inFile.close();

	return setupValues;
}

SetupComp setConfigComp(int argc, char * argv[])
{
	SetupComp setupValues;

	static int verbose_flag;
	static int debug_flag;
	static int time_flag;
	static int remdiag_flag;
	static int random_flag;

	static struct option long_options[] =
	{
			/* These options set a flag. */
			{"verbose", no_argument,       &verbose_flag, 1},
			{"debug", no_argument,       &debug_flag, 1},
			{"time", no_argument,       &time_flag, 1},
			{"norandom", no_argument,       &random_flag, 1},
			/* These options don’t set a flag.
			             We distinguish them by their indices. */
			{"sample",     required_argument,       0, 's'},
			{"control",     required_argument,       0, 'c'},
			{"samplename",  required_argument,       0, 'n'},
			{"digest",  required_argument, 0, 'd'},
			{"baits",  required_argument, 0, 'b'},
			{"threads",  required_argument, 0, 't'},
			{"res",    required_argument, 0, 'r'},
			{"output",     no_argument,       0, 'o'},
			{"cistrans",  no_argument,       0, 'C'},
			{"algorithm",  required_argument, 0, 'A'},
			{"alpha",  required_argument, 0, 'a'},
			{"log",    required_argument, 0, 'l'},

			{0, 0, 0, 0}
	};

	int option_index = 0;
	int opt;
	string value;

	while ((opt = getopt_long(argc, argv, "s:c:n:d:b:t:r:o:C:A:a:l:", long_options, &option_index)) != -1)
	{

		switch(opt)
		{
		case 'c':
			setupValues.setCondition1(optarg);
			break;
		case 's':
			setupValues.setCondition2(optarg);
			break;
		case 'n':
			setupValues.setSname(optarg);
			break;
		case 'd':
			setupValues.setEnzyme(optarg);
			break;
		case 'b':
			setupValues.setBaits(optarg);
			break;
		case 't':
			setupValues.setThreads(atoi(optarg));
			break;
		case 'r':
			setupValues.setRes(atoi(optarg));
			break;
		case 'o':
			value = removeSpaces(optarg);
			if (! boost::algorithm::ends_with(value, "/"))
				value += "/";
			setupValues.setOutDir(value);
			break;
		case 'C':
			setupValues.setCisTrans(optarg);
			break;
		case 'A':
			setupValues.setQvalue(optarg);
			break;
		case 'a':
			setupValues.setAlpha(optarg);
			break;
		case 'l':
			setupValues.setLogFile(optarg);
			break;
		}//*/
	}



	if (verbose_flag && !debug_flag)
		setupValues.setVerbose(vl_info);
	if (debug_flag)
		setupValues.setVerbose(vl_debug);
	if (random_flag)
		setupValues.setRandom(false);
	if (time_flag)
			setupValues.setTime(true);

	return setupValues;
}

void SetupComp::setQvalue(string L)
{
	std::map<std::string,QV_Options> QVoptionValues;
	QVoptionValues["ihw"] = qv_ihw;
	QVoptionValues["bh"] = qv_bh;
	QVoptionValues["IHW"] = qv_ihw;
	QVoptionValues["BH"] = qv_bh;

	if (L.empty())
		mQvalue = qv_ihw;
	else
		mQvalue = QVoptionValues[L];
}

void SetupComp::setCisTrans(string input)
{
	std::map<std::string,CisTrans> CToptionValues;
	CToptionValues["all"] = ct_all;
	CToptionValues["cis"] = ct_cis;
	CToptionValues["trans"] = ct_trans;

	mCisTrans = CToptionValues[input];
}
