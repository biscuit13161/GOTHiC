/*
 * Utils.h
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#ifndef SRC_SETUP_H_
#define SRC_SETUP_H_

#include "Utils.h"
#include <string>
#include <map>

enum AnalysisOptions {
	ao_single,
	ao_comparative
};

class Setup
{
private:
	std::string mOutDir;
	std::string mEnzyme;

	std::string mInput;
	std::string mSname;
	std::string mCliName;
	std::string mLogFile;
	CisTrans mCisTrans;
	int mThreads;
	int mRes;
	bool mRemoveDiagonal;
	AnalysisOptions mAnalysisType;
	bool mVerbose;
public:
	Setup();
	Setup(std::string outDir, std::string enzyme, std::string input, int threads);

	void print();

	inline int getThreads() const {return mThreads;}
	inline int getRes() const {return mRes;}
	inline std::string getEnzyme() const {return mEnzyme;}
	inline std::string getOutDir() const {return mOutDir;}
	inline std::string getInput() const {return mInput;}
	inline std::string getSname() const {return mSname;}
	inline std::string getCliName() const {return mCliName;}
	inline std::string getLogFile() const {return mLogFile;}
	inline CisTrans getCisTrans() const {return mCisTrans;}
	inline bool getRemoveDiagonal() const {return mRemoveDiagonal;}
	inline AnalysisOptions getAnalysisType() const {return mAnalysisType;}
	inline bool getVerbose() const {return mVerbose;}

	inline void setThreads(int threads) {mThreads = threads;}
	inline void setRes(int res) {mRes = res;}
	inline void setEnzyme(std::string enzyme) {mEnzyme = enzyme;}
	inline void setSname(std::string name) {mSname = name;}
	inline void setCliName(std::string name) {mCliName = name;}
	inline void setOutDir(std::string outDir) {mOutDir = outDir;}
	inline void setInput(std::string input) {mInput = input;}
	inline void setLogFile(std::string L) {mLogFile = L;}
	inline void setCisTrans(CisTrans input) {mCisTrans = input;}
	inline void setRemoveDiagonal(bool L) {mRemoveDiagonal = L;}
	inline void setAnalysisType(AnalysisOptions L) {mAnalysisType = L;}
	inline void setVerbose(bool L) {mVerbose = L;}

};

Setup loadConfig(std::string & fileName);

enum Options {
	Input = 0,
	Sname,
	Digest,
	Threads,
	Res,
	Output,
	Cistrans,
	Analysis,
	RemDiag,
	Logfile,
	Verbose,
};




#endif /* SRC_SETUP_H_ */
