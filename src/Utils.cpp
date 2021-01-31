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
 

#include "Utils.h"
#include "Site.h"
#include "Interactions.h"
#include "version.h"
#include <cmath>
#include <ctime>
#include <cstdarg>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <boost/filesystem.hpp>
#include "tbb/parallel_for.h"

using namespace std;
using namespace tbb;

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

void completed(int n)
{
	cerr << "\t... Completed";
	if (n > 0)
		cerr << " (" << n << ")";
	cerr <<	"! ";
	//showTime();
	cerr << endl;
}

void printUsage()
{
	cerr << "GOTHiC++ usage:" << endl;
	cerr << string("version: ") << GOTH_VERSION_STRING << endl << endl;
	cerr << "    gothic <path/to/config/file>" << endl;
	cerr << "OR" << endl;
	cerr << "    gothic -i <filename> -s <name> -d <filename> [-t #] [-r #] [-o <dir>] [-c (all|trans|cis)]" << endl;
	cerr << "      [-A (single|comparative)] [--verbose] [--no_rem_diag]" << endl << endl;
	cerr << "options:" << endl;
	cerr <<
			"    -i <filename>         Input filename\n"
			"      --input <filename>\n"
			"    -s <name>             Sample name\n"
			"      --sample <name>\n"
			"    -d <filename>         Digest of Restriction Enzyme, as used by HiCUP\n"
			"      --digest <filename>\n"
			"    -t #                  Num of threads to run, defaults to 1\n"
			"      --threads #\n"
			"    -r #                  Resolution in bases for bining interactions, defaults to 10000\n"
	        "      --res #\n"
			"    -o <dir>              Output directory, defaults to './'\n"
			"      --output <dir>\n"
			"    -c (all|trans|cis)    Filter for Cis or Trans interactions,defaults to 'all'\n"
			"      --cistrans (all|trans|cis)\n"
			"    -A <option>           Analysis type, either 'single' or 'comparative'\n"
			"      --analysis <option>\n"
//			"    -l <filename>         Log file\n"
//			"      --log <filename>\n"
			"    --verbose             Print verbose output during run\n"
			"    --no_rem_diag         Flag to stop removal of diagonals during analysis\n";
	cerr << endl;


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

void writeBinary(vector<Site> & sites, string binOutFileName)
{
	ofstream binOutFile;
	binOutFile.open(binOutFileName, ios::binary);

	int max = 0;
	if (binOutFile.is_open())
	{
		for (int i = 0; i != sites.size(); i++)
		{
			Site P = sites[i];
			string chr = P.getChr();
			if (!chr.empty())
			{
				chr.resize(32);
				int start = P.getStart();
				int end = P.getEnd();
				//size_t S = sizeof(std::string) + 3* sizeof(int);
				binOutFile.write(chr.c_str(), 32);// need to cast the pointer
				binOutFile.write(reinterpret_cast<char*>(&start), sizeof(int));
				binOutFile.write(reinterpret_cast<char*>(&end), sizeof(int));
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

void readBinary(vector<Site> & sites, string binInFileName)
{
	cerr << "Reading Binary input file" << endl;
	ifstream binInFile;
	binInFile.open(binInFileName, ios::binary);

	if (binInFile.is_open())
	{
		while(binInFile)
		{
			string chr;
			chr.resize(32);
			//binInFile.read(reinterpret_cast<char*>(& chr), 32); // need to cast the pointer
			binInFile.read(& chr[0], 32);
			chr.resize(chr.find('\0'));
			if (!chr.empty())
			{
				int start;
				binInFile.read(reinterpret_cast<char*>(& start), sizeof(int));
				int end;
				binInFile.read(reinterpret_cast<char*>(& end), sizeof(int));
				Site P = Site(chr,start,end);
				sites.push_back(P);
			}
		}
		binInFile.close();
	}
	else
	{
		string str = string("Unable to read binary file: ") + binInFileName;
			throw std::invalid_argument(str);
	}

	completed();
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
	int pos = 0;
	//for (auto it = interactions.begin(); it != interactions.end(); it++)
	//{
//#pragma omp parallel for
	//for (int i = 0; i< interactions.size(); i++)
	//{
	parallel_for(size_t(0),size_t(interactions.size()),
						[&] (size_t i) {
		//auto it = interactions.begin();
		//std::advance(it, i);//*/
		//Interaction T = *it;
		Interaction T = interactions[i];
		if (T.getInt1() != T.getInt2())
		{

			binned_df_filtered.push_back(T);
			pos++;
		}
	});
}
void removeDiagonals(concurrent_vector<Interaction> & interactions, CisTrans cistrans, bool removeDiagonal)
{
	/** Remove Diagonals and Cis/Trans filter values **/

	concurrent_vector<Interaction> binned_df_filtered;

	if (removeDiagonal)
	{
		cerr << "\tRemoving Diagonals!" << endl;
		removeDuplicates(interactions, binned_df_filtered);
		cerr << "\t";
		completed();
		interactions.clear();
	}
	else
	{
		binned_df_filtered.swap(interactions);
	}
	//cout << "size " << binned_df_filtered.size() << endl;

	if (cistrans == ct_cis)
	{
		cerr << "\tFinding Cis interactions!";
		findCis(interactions, binned_df_filtered);
		//interactions.resize(pos);
		cerr << "|\t";
		completed();
	}
	else if (cistrans == ct_trans)
	{
		cerr << "\tFinding Trans interactions!";
		findTrans(interactions, binned_df_filtered);
		//interactions.resize(pos);
		cerr << "|\t";
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

void findCis(concurrent_vector<Interaction> & interactions, concurrent_vector<Interaction> & binned_df_filtered)
{
	int pos = 0;
//#pragma omp parallel for
	//for (int i = 0; i< binned_df_filtered.size(); i++)
	//{
	parallel_for(size_t(0),size_t(binned_df_filtered.size()),
						[&] (size_t i) {
		//auto it = binned_df_filtered.begin();
		//advance(it, i);
		//Interaction T = *it;
		Interaction T = binned_df_filtered[i];
		if (T.getChr1() == T.getChr2())
		{
//#pragma omp critical (fc1)
			interactions.push_back(T);
			pos++;
		}
	});
}

void findTrans(concurrent_vector<Interaction> & interactions, concurrent_vector<Interaction> & binned_df_filtered)
{
	int pos = 0;
//#pragma omp parallel for
	//for (int i = 0; i< binned_df_filtered.size(); i++)
	//{
		parallel_for(size_t(0),size_t(binned_df_filtered.size()),
							[&] (size_t i) {
		//auto it = binned_df_filtered.begin();
		//std::advance(it, i);
		//Interaction T = *it;
		Interaction T = binned_df_filtered[i];
		if (T.getChr1() != T.getChr2())
		{
//#pragma omp critical (ft1)
			interactions.push_back(T);
			pos++;
		}
	});
}

void getSumSquare(double & sumSquare, const set<string> & chromos, const concurrent_vector<Interaction> & interactions)
{
	/** counts unique interactions **/
	double s2 = 0;

	for (auto cr : chromos){
		set<int> pos;
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
		sumSquare += pow(pos.size(),2);
	}
}

void writeBinary(concurrent_vector<Interaction> & interactions, string binOutFileName)
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

void readBinary(concurrent_vector<Interaction> & interactions, string binInFileName)
{
	cerr << "Reading Binary input file" << endl;
	ifstream binInFile;
	binInFile.open(binInFileName, ios::binary);

	if (binInFile.is_open())
	{
		while(binInFile)
		{
			string Int1;
			Int1.resize(32);
			//binInFile.read(reinterpret_cast<char*>(& chr), 32); // need to cast the pointer
			binInFile.read(& Int1[0], 32);
			Int1.resize(Int1.find('\0'));
			if (!Int1.empty())
			{
				size_t pos = Int1.find(":");
				string chr1 = Int1.substr(0,pos);
				int locus1 = atoi(Int1.substr(pos+1).c_str());

				string Int2;
				Int2.resize(32);
				//binInFile.read(reinterpret_cast<char*>(& chr), 32); // need to cast the pointer
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

	completed();
}

int binarySearch(int arr[], int l, int r, int x)
{
	// x = search value
	// l = first (left) element index (0)
	// r = last (right) element index (size() - 1)
	// arr[] = array of elements
	if (r >= 1)
	{
		int mid = l + (r-1) /2;
		//check if mid point matches
		if (arr[mid] == x)
			return mid;
		// check if lower than mid
		if (arr[mid] > x)
			return binarySearch(arr, l, mid-1, x);
		//else check above mid
		return binarySearch(arr, mid+1, r, x);
	}
	return -1;
}

int binarySearch2(int arr[], int l, int r, int x)
{
	// x = search value
	// l = first (left) element index (0)
	// r = last (right) element index (size() - 1)
	// arr[] = array of elements
	int m = 0;

	while (r - l > 1)
	{
		m =  l + (r-1) / 2;
		if (arr[m] <= x)
			l = m;
		else
			r = m;
	}
	if (arr[l] == x)
		return l;
	else
		return -1;
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
