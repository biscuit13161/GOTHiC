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
#include <ctime>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

using namespace std;

SetupComp::SetupComp(): mOutDir(""), mEnzyme(""), mCondition1(""), mCondition2(""), mThreads(0), mRes(10000), mRemoveDiagonal(true), mVerbose(false)
{

}

void SetupComp::print()
{
	cerr << endl << string("GOTHiC++ v") << GOTH_MAJOR_VERSION << "." << GOTH_MINOR_VERSION << "." << GOTH_PATCH_VERSION << endl << endl;
	cerr << "#" << endl;
	cerr << "# Output Directory:  " << mOutDir << endl;
	cerr << "# Restriction File:  " << mEnzyme << endl;
	cerr << "# Baits File:        " << mBaits << endl;
	cerr << "# Sample Name:       " << mSname << endl;
	cerr << "# Control File:      " << mCondition1 << endl;
	cerr << "# Sample File:       " << mCondition2 << endl;
	cerr << "# Threads:           " << mThreads << endl;
	cerr << "# Resolution:        " << mRes << endl;
	cerr << "# Config File:       " << mCliName << endl;
	cerr << "# Verbose:           " << mVerbose << endl;
	if (mQvalue == qv_ihw)
		cerr << "# Pvalue Correction: Independent Hypothesis Weighting" <<endl;
	else if (mQvalue == qv_bh)
		cerr << "# Pvalue Correction: Benjamini Hochberg" <<endl;
	else
		cerr << "# Pvalue Correction: Not Specfied!" <<endl;
	cerr << "# Alpha:             " << mAlpha << endl;
	cerr << "#" << endl << endl;
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
			optionValues["Alpha"] = sc_Alpha;
			optionValues["Algorithm"] = sc_Qvalues;

			std::map<std::string,CisTrans> CToptionValues;
			CToptionValues["all"] = ct_all;
			CToptionValues["cis"] = ct_cis;
			CToptionValues["trans"] = ct_trans;

			std::map<std::string,QV_Options> QVoptionValues;
			QVoptionValues["ihw"] = qv_ihw;
			QVoptionValues["bh"] = qv_bh;
			QVoptionValues["IHW"] = qv_ihw;
			QVoptionValues["BH"] = qv_bh;


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
					setupValues.setOutDir(value);
				}
				break;
			case sc_Cistrans:
				if (value.empty())
					setupValues.setCisTrans(ct_all);
				else
					setupValues.setCisTrans(CToptionValues[value]);
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
					setupValues.setVerbose(true);
				break;
			case sc_Alpha:
				if (value.empty())
					setupValues.setAlpha(string("0.1"));
				else
					setupValues.setAlpha(value);
				break;
			case sc_Qvalues:
				if (value.empty())
					setupValues.setQvalue(qv_ihw);
				else
					setupValues.setQvalue(QVoptionValues[value]);
				break;
			}//*/
		}
	}
	inFile.close();

	return setupValues;
}
