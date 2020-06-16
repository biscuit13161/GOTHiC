/*
 * SetupComp.h
 *
 *  Created on: 27 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#ifndef SRC_SETUPCOMP_H_
#define SRC_SETUPCOMP_H_

#include "UtilsComp.h"
#include <string>
#include <vector>

class SetupComp
{
private:
	std::string mOutDir;
	std::string mEnzyme;
	std::string mBaits;
	std::string mCondition1;
	std::string mCondition2;
	std::string mSname;
	std::string mCliName;
	std::string mLogFile;
	CisTrans mCisTrans;
	int mThreads;
	int mRes;
	bool mRemoveDiagonal;
	bool mVerbose;
public:
	SetupComp();
	SetupComp(std::string outDir, std::string enzyme, std::string input, int threads);

	void print();

	inline int getThreads() const {return mThreads;}
	inline int getRes() const {return mRes;}
	inline std::string getEnzyme() const {return mEnzyme;}
	inline std::string getBaits() const {return mBaits;}
	inline std::string getOutDir() const {return mOutDir;}
	inline std::string getCondition1() const {return mCondition1;}
	inline std::string getCondition2() const {return mCondition2;}
	inline std::string getSname() const {return mSname;}
	inline std::string getCliName() const {return mCliName;}
	inline std::string getLogFile() const {return mLogFile;}
	inline CisTrans getCisTrans() const {return mCisTrans;}
	inline bool getRemoveDiagonal() const {return mRemoveDiagonal;}
	inline bool getVerbose() const {return mVerbose;}

	inline void setThreads(int threads) {mThreads = threads;}
	inline void setRes(int res) {mRes = res;}
	inline void setEnzyme(std::string enzyme) {mEnzyme = enzyme;}
	inline void setBaits(std::string baits) {mBaits = baits;}
	inline void setSname(std::string name) {mSname = name;}
	inline void setCliName(std::string name) {mCliName = name;}
	inline void setOutDir(std::string outDir) {mOutDir = outDir;}
	inline void setCondition1(std::string input) {mCondition1 = input;}
	inline void setCondition2(std::string input) {mCondition2 = input;}
	inline void setLogFile(std::string L) {mLogFile = L;}
	inline void setCisTrans(CisTrans input) {mCisTrans = input;}
	inline void setRemoveDiagonal(bool L) {mRemoveDiagonal = L;}
	inline void setVerbose(bool L) {mVerbose = L;}

};

SetupComp loadConfigComp(std::string & fileName);

enum SC_Options {
	sc_Cond1 = 0,
	sc_Cond2,
	sc_Sname,
	sc_Digest,
	sc_Baits,
	sc_Threads,
	sc_Res,
	sc_Output,
	sc_Cistrans,
	sc_RemDiag,
	sc_Logfile,
	sc_Verbose,
};

#endif /* SRC_SETUPCOMP_H_ */
