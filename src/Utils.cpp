/*
 * Utils.cpp
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#include "Utils.h"
#include "Site.h"
#include "Interactions.h"
#include "version.h"
#include <ctime>
#include <cstdarg>
#include <stdint.h>
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
	cerr << "GOTHiC++ usage:" << endl << endl;
	cerr << "    gothic <path/to/config/file>"  << endl << endl;
	cerr << string("version: ") << GOTH_MAJOR_VERSION << "." << GOTH_MINOR_VERSION << "." << GOTH_PATCH_VERSION << endl;

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
		cout << "Unable to read binary file: " + binInFileName << endl;
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

void removeDuplicates(vector<Interaction> & interactions, vector<Interaction> & binned_df_filtered)
{
	int pos = 0;
	//for (auto it = interactions.begin(); it != interactions.end(); it++)
	//{
#pragma omp parallel for
	for (int i = 0; i< interactions.size(); i++)
	{
		auto it = interactions.begin();
		std::advance(it, i);//*/
		Interaction T = *it;
		if (T.getInt1() != T.getInt2())
		{
#pragma omp critical (rd1)
			binned_df_filtered.push_back(T);
			pos++;
		}
	}
}

void removeDiagonals(vector<Interaction> & interactions, CisTrans cistrans, bool removeDiagonal)
{
	/** Remove Diagonals and Cis/Trans filter values **/

	vector<Interaction> binned_df_filtered;

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

void findCis(vector<Interaction> & interactions, vector<Interaction> & binned_df_filtered)
{
	int pos = 0;
#pragma omp parallel for
	for (int i = 0; i< binned_df_filtered.size(); i++)
	{
		auto it = binned_df_filtered.begin();
		advance(it, i);
		//	for (auto it = binned_df_filtered.begin(); it != binned_df_filtered.end(); it++)
		//	{
		Interaction T = *it;
		if (T.getChr1() == T.getChr2())
		{
#pragma omp critical (fc1)
			interactions.push_back(T);
			pos++;
		}
	}
}

void findTrans(vector<Interaction> & interactions, vector<Interaction> & binned_df_filtered)
{
	int pos = 0;
#pragma omp parallel for
	for (int i = 0; i< binned_df_filtered.size(); i++)
	{
		auto it = binned_df_filtered.begin();
		std::advance(it, i);

		Interaction T = *it;
		if (T.getChr1() != T.getChr2())
		{
#pragma omp critical (ft1)
			interactions.push_back(T);
			pos++;
		}
	}
}
