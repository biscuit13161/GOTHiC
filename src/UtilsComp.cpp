/*
 * UtilsComp.cpp
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#include "UtilsComp.h"
#include "Site.h"
#include "Interactions.h"
#include "version.h"
#include <cmath>
#include <ctime>
#include <cstdarg>
#include <stdint.h>
#include <iostream>
#include <fstream>
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

void printUsageComp()
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
	//int pos = 0;
	//for (auto it = interactions.begin(); it != interactions.end(); it++)
	//{
	parallel_for(
			blocked_range<concurrent_vector<Interaction>::iterator>(interactions.begin(),	interactions.end()),
			[&] (blocked_range<concurrent_vector<Interaction>::iterator> inter) {
		for (auto T : inter)
	{
		if (T.getInt1() != T.getInt2())
		{
			binned_df_filtered.push_back(T);
			//pos++;
		}
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
	//int pos = 0;
	//int i = 0; i< binned_df_filtered.size(); i++)
	parallel_for(
			blocked_range<concurrent_vector<Interaction>::iterator>(binned_df_filtered.begin(),	binned_df_filtered.end()),
			[&] (blocked_range<concurrent_vector<Interaction>::iterator> inter) {
		for (auto T : inter)
		{
			//auto it = binned_df_filtered.begin();
			//advance(it, i);
			//Interaction T = *it;
			if (T.getChr1() == T.getChr2())
			{
				interactions.push_back(T);
				//pos++;
			}
		}
	});
}

void findTrans(concurrent_vector<Interaction> & interactions, concurrent_vector<Interaction> & binned_df_filtered)
{
	//int pos = 0;
//#pragma omp parallel for
	//for (int i = 0; i< binned_df_filtered.size(); i++)

	parallel_for(
			blocked_range<concurrent_vector<Interaction>::iterator>(binned_df_filtered.begin(),	binned_df_filtered.end()),
			[&] (blocked_range<concurrent_vector<Interaction>::iterator> inter) {
		for (auto T : inter)
		{
			if (T.getChr1() != T.getChr2())
			{
				interactions.push_back(T);
				//pos++;
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

void readBinary(tbb::concurrent_vector<Interaction> & interactions, string binInFileName)
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
		cout << "Unable to read binary file: " + binInFileName << endl;
	}

	completed();
}

void checkInputFiles(std::string file)
{
	if (file.find("inter.bin",file.length()-9) == string::npos )
		{
		string str = string("Incorrect input file: ") + file + "\nInput files should be GOTHiC output files for comparative analysis\n";
		throw std::invalid_argument(str);
		}
}
