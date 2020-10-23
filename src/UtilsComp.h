/*
 * UtilsComp.h
 *
 *  Created on: 11 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

//#include "Interactions.h"
//#include "Site.h"
#include <cstdint>
#include <regex>
#include <vector>
#include <set>
#include <string>
#include "stdlib.h"
#include "sys/types.h"
#include "sys/sysinfo.h"
#include "tbb/concurrent_vector.h"

class Interaction;
class Site;

enum CisTrans
{
	ct_all=0,
	ct_cis,
	ct_trans,
};

enum QV_Options {
	qv_ihw = 0,
	qv_bh
};

void showTime();
void completed(int n = 0);
void printUsageComp();
void verbosePrint(std::string & str, bool verbose = true);

void verbose(const char * fmt, ... );

std::uint32_t fact(std::uint32_t n);
int fact(int n);

void writeBinary(tbb::concurrent_vector<Site> & sites, std::string binOutFileName);
void readBinary(tbb::concurrent_vector<Site> & sites, std::string binInFileName);

//get memory usage
int getRealValue();
int getVirtValue();
int parseLine(char* line);

void removeDuplicates(tbb::concurrent_vector<Interaction> & interactions, tbb::concurrent_vector<Interaction> & binned_df_filtered);
void removeDiagonals(tbb::concurrent_vector<Interaction> & interactions, CisTrans cistrans, bool removeDiagonal);

void findTrans(tbb::concurrent_vector<Interaction> & interactions, tbb::concurrent_vector<Interaction> & binned_df_filtered);
void findCis(tbb::concurrent_vector<Interaction> & interactions, tbb::concurrent_vector<Interaction> & binned_df_filtered);

void writeBinary(tbb::concurrent_vector<Interaction> & interactions, std::string binOutFileName);
void readBinary(tbb::concurrent_vector<Interaction> & interactions, std::string binInFileName);

std::string fixChromosomeNames(std::string chr);
void getSumSquare(double & sumSquare, const std::set<std::string> & chromos,const tbb::concurrent_vector<Interaction> & interactions);

void checkInputFiles(std::string file);

char *removeSpaces(std::string & str);
char *removeSpaces(char *str);

#endif /* SRC_UTILSCOMP_H_ */
