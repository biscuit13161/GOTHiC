/*
 *  UtilsComp.cpp
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
 

#include "UtilsComp.h"
#include "Site.h"
#include "SetupComp.h"
#include "Interactions.h"
#include "version.h"
#include <cmath>
#include <ctime>
#include <cstdarg>
#include <cstdio>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <boost/filesystem.hpp>
#include "tbb/concurrent_vector.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"

using namespace std;
using namespace tbb;

void showTime()
{
	// current date/time based on current system
	time_t now = time(0);

	// convert now to string form
	//char* dt = ctime(&now);
	tm *ltm = localtime(&now);

	cout << ltm->tm_hour << ":" << flush;
	cout << ltm->tm_min << ":" << flush;
	cout << ltm->tm_sec << flush;
}

void completed(int n)
{
	cout << "\t... Completed";
	if (n > 0)
		cout << " (" << n << ")";
	cout <<	"! ";
	//showTime();
	cout << endl;
}

void printUsageComp()
{

	cout << "GOTHiC++ usage:" << endl;
	cout << string("version: ") << GOTH_VERSION_STRING << endl << endl;
	cout << "    gothicomp <path/to/config/file>" << endl;
	cout << "OR" << endl;
	cout << "    gothicomp -c <filename> -s <filename> -n <name> -d <filename> [-b <filename>] [-t #] [-r #]" << endl;
	cout << "        [-o <dir>] [-a #] [-A (bh|ihw)] [-C (all|trans|cis)] [--norandom] [--verbose|--debug]" << endl << endl;
	cout << "options:" << endl;
	cout <<
			"    -c <filename>         Control input filename\n"
			"      --control <filename>\n"
			"    -s <filename>         Sample input filename\n"
			"      --sample <filename>\n"
			"    -n <name>             Sample name\n"
			"      --sample <name>\n"
			"    -d <filename>         Digest of Restriction Enzyme, as used by HiCUP\n"
			"      --digest <filename>\n"
			"    -b <filename>         Baits file for analysis, as used by HiCUP\n"
			"      --baits <filename>\n"
			"    -t #                  Num of threads to run, defaults to 1\n"
			"      --threads #\n"
			"    -r #                  Resolution in bases for binning interactions, defaults to 10000\n"
	        "      --res #               - Only effective for binning the baits file, Sample data is binned using gothic.\n"
			"    -o <dir>              Output directory, defaults to './'\n"
			"      --output <dir>\n"
			"    -C (all|trans|cis)    Filter for Cis or Trans interactions,defaults to 'all'\n"
			"      --cistrans (all|trans|cis)\n"
			"    -A (bh|ihw)           Algorithm for p-value correction, either 'bh' or 'ihw'\n"
			"      --analysis (bh|ihw)\n"
			"    -a #                  Alpha cutoff for p-value correction, defaults to '0.1'\n"
			"      --alpha #\n"
//			"    -l <filename>         Log file\n"
//			"      --log <filename>\n"
			"    --norandom            Turn off Random subsampling\n"
			"    --verbose             Print verbose output during run\n"
			"    --debug               Print very verbose output during run\n";
	cout << endl;

}

void verbosePrint(string & str, bool verbose)
{
	if (verbose)
	cout << str << endl;
}


void verbose(const char * fmt, ... )
{
va_list args;  /* Used as a pointer to the next variable argument. */
va_start( args, fmt );  /* Initialize the pointer to arguments. */

#ifdef DEBUG_FLAG
printf(fmt, &args);
#endif
/*This isn't tested, the point is to be able to pass args to
printf*/
va_end(args);
}

uint32_t fact(uint32_t n)
{
   if ((n==0)||(n==1))
      return 1;
   else
      return n*fact(n-1);
}

int fact(int n)
{
    //  Factorial
   if ((n==0)||(n==1))
      return 1;
   else
      return n*fact(n-1);
}



int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getVirtValue(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

int getRealValue(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

void removeDuplicates(concurrent_vector<Interaction> & interactions, concurrent_vector<Interaction> & binned_df_filtered)
{
	parallel_for(
			blocked_range<concurrent_vector<Interaction>::iterator>(interactions.begin(),	interactions.end()),
			[&] (blocked_range<concurrent_vector<Interaction>::iterator> inter) {
		for (auto T : inter)
	{
		if (T.getInt1() != T.getInt2())
		{
			binned_df_filtered.push_back(T);
		}
	}
});
}


void removeDiagonals(concurrent_vector<Interaction> & interactions, CisTrans cistrans, bool removeDiagonal, string type)
{
	/** Remove Diagonals and Cis/Trans filter values **/

	concurrent_vector<Interaction> binned_df_filtered;

	if (removeDiagonal)
	{
		cout << "\tRemoving " << type << " Diagonals!" << endl;
		removeDuplicates(interactions, binned_df_filtered);
		cout << "\t";
		completed();
		interactions.clear();
	}
	else
	{
		binned_df_filtered.swap(interactions);
	}

	if (cistrans == ct_cis)
	{
		cout << "\tFinding " << type << " Cis interactions!";
		findCis(interactions, binned_df_filtered);
		cout << "\t";
		completed();
	}
	else if (cistrans == ct_trans)
	{
		cout << "\tFinding " << type << "Trans interactions!";
		findTrans(interactions, binned_df_filtered);
		cout << "\t";
		completed();
	}
	else
	{
		interactions.swap(binned_df_filtered);
	}

}

string fixChromosomeNames(string chr)
{
	//capital to small
	chr = "chr" + chr;
	chr = regex_replace(chr, regex("CHR"), "chr");
	chr = regex_replace(chr, regex("chrchr"), "chr");
	return chr;
}

void getSumSquare(double & sumSquare, const set<string> & chromos, const concurrent_vector<Interaction> & interactions)
{
	/** counts unique interactions **/
	double s2 = 0;

	concurrent_vector<string> chrom(chromos.begin(),chromos.end());

	typedef concurrent_vector<string> csstring;
	typedef blocked_range<csstring::iterator> range;
	sumSquare = parallel_reduce(range(chrom.begin(), chrom.end()),sumSquare,
			[&] (range & inter, int init)
			{
		for (auto cr : inter){
			set<double> pos;
			for (auto x : interactions)
			{
				if (x.getChr1() == cr)
				{
					pos.insert(x.getLocus1());
				}
				if (x.getChr2() == cr)
				{
					pos.insert(x.getLocus2());
				}
			}
			init += pow(pos.size(),2);
		}
		return init;
		}, [](double x, double y){return x+y;}
	);
}

void findCis(concurrent_vector<Interaction> & interactions, concurrent_vector<Interaction> & binned_df_filtered)
{
	parallel_for(
			blocked_range<concurrent_vector<Interaction>::iterator>(binned_df_filtered.begin(),	binned_df_filtered.end()),
			[&] (blocked_range<concurrent_vector<Interaction>::iterator> inter) {
		for (auto T : inter)
		{
			if (T.getChr1() == T.getChr2())
			{
				interactions.push_back(T);
			}
		}
	});
}

void findTrans(concurrent_vector<Interaction> & interactions, concurrent_vector<Interaction> & binned_df_filtered)
{
	parallel_for(
			blocked_range<concurrent_vector<Interaction>::iterator>(binned_df_filtered.begin(),	binned_df_filtered.end()),
			[&] (blocked_range<concurrent_vector<Interaction>::iterator> inter) {
		for (auto T : inter)
		{
			if (T.getChr1() != T.getChr2())
			{
				interactions.push_back(T);
			}
		}
	});
}

void writeBinary(tbb::concurrent_vector<Interaction> & interactions, string binOutFileName)
{
	ofstream binOutFile;
	binOutFile.open(binOutFileName, ios::binary);

	int max = 0;
	if (binOutFile.is_open())
	{
		for (int i = 0; i != interactions.size(); i++)
		{
			Interaction P = interactions[i];
			string Int1 = P.getInt1();
			string Int2 = P.getInt2();
			if (!Int1.empty())
			{
				Int1.resize(32);
				Int2.resize(32);
				int freq = P.getFreq();
				//size_t S = sizeof(std::string) + 3* sizeof(int);
				binOutFile.write(Int1.c_str(), 32);
				binOutFile.write(Int2.c_str(), 32);// need to cast the pointer
				binOutFile.write(reinterpret_cast<char*>(&freq), sizeof(int));
			}
		}

		binOutFile.close();
	}
	else
	{
		cout << "Unable to create binary write file: " + binOutFileName << endl;
	}
	cout << endl;

}

void readBinary(tbb::concurrent_vector<Interaction> & interactions, string binInFileName, std::string type, Verb_level verb_lev)
{
	cout << "Reading " << type << " Binary input file" << endl;
	ifstream binInFile;
	binInFile.open(binInFileName, ios::binary);

	if (binInFile.is_open())
	{
		while(binInFile)
		{
			string Int1;
			Int1.resize(32);
			binInFile.read(& Int1[0], 32);
			Int1.resize(Int1.find('\0'));
			if (!Int1.empty())
			{
				size_t pos = Int1.find(":");
				string chr1 = Int1.substr(0,pos);
				int locus1 = atoi(Int1.substr(pos+1).c_str());

				string Int2;
				Int2.resize(32);
				binInFile.read(& Int2[0], 32);
				Int2.resize(Int2.find('\0'));
				pos = Int2.find(":");
				string chr2 = Int2.substr(0,pos);
				int locus2 = atoi(Int2.substr(pos+1).c_str());

				int freq;
				binInFile.read(reinterpret_cast<char*>(& freq), sizeof(int));
				Interaction P = Interaction(chr1,chr2,locus1,locus2,freq);
				interactions.push_back(P);
			}
		}
		binInFile.close();
	}
	else
	{
		string str = string("Unable to read binary file: ") + binInFileName;
		throw std::invalid_argument(str);
	}

	if (verb_lev >= vl_info)
	{
		fprintf(stderr, "\t%'d interactions loaded\n", interactions.size() );
	}

	completed();
}

char *removeSpaces(string & str)
{
	return removeSpaces(&str[0]);
}

char *removeSpaces(char *str)
{
    int i = 0, j = 0;
    while (str[i])
    {
        if (str[i] != ' ')
           str[j++] = str[i];
        i++;
    }
    str[j] = '\0';
    return str;
}
