/*
 * Utils.h
 *
 *  Created on: 11 May 2020
 *      Author: rich
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

class Interaction;
class Site;

enum CisTrans
{
	ct_all=0,
	ct_cis,
	ct_trans,
};

void showTime();
void completed(int n = 0);
void printUsage();
void verbosePrint(std::string & str, bool verbose = true);

void verbose(const char * fmt, ... );

std::uint32_t fact(std::uint32_t n);
int fact(int n);

void writeBinary(std::vector<Site> & sites, std::string binOutFileName);
void readBinary(std::vector<Site> & sites, std::string binInFileName);

//get memory usage
int getRealValue();
int getVirtValue();
int parseLine(char* line);

void removeDuplicates(std::vector<Interaction> & interactions, std::vector<Interaction> & binned_df_filtered);
void removeDiagonals(std::vector<Interaction> & interactions, CisTrans cistrans, bool removeDiagonal);

void findTrans(std::vector<Interaction> & interactions, std::vector<Interaction> & binned_df_filtered);
void findCis(std::vector<Interaction> & interactions, std::vector<Interaction> & binned_df_filtered);

void writeBinary(std::vector<Interaction> & interactions, std::string binOutFileName);
void readBinary(std::vector<Interaction> & interactions, std::string binInFileName);


std::string fixChromosomeNames(std::string chr);
void getSumSquare(double & sumSquare, const std::set<std::string> & chromos,const std::vector<Interaction> & interactions);

#endif /* SRC_UTILS_H_ */
