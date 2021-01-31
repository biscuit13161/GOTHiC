/*
 * JointUtils.cpp
 *
 *  Created on: 31 Jan 2021
 *      Author: rich
 */

#include "JointUtils.h"
#include <cstdint>
#include <regex>
#include <vector>
#include <set>
#include <string>
#include <filesystem>
#include "stdlib.h"
#include "sys/types.h"

using namespace std;

void checkInputFiles(std::string file)
{
	checkFileExists(file, "Input");

	if (file.find("inter.bin",file.length()-9) == string::npos )
		{
		string str = string("Incorrect input file: ") + file + "\nInput files should be GOTHiC output files for comparative analysis\n";
		throw std::invalid_argument(str);
		}
}

void checkFileExists(std::string file)
{
	checkFileExists(file,"");
}

void checkFileExists(std::string file, std::string type)
{
	namespace fs = std::filesystem;
	fs::path f{ file };
	string str = string("importHicup: ") + type +" file not found";
	if (! fs::exists(f))
		throw std::invalid_argument(str);
}

