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
	Verb_level mVerbose;
	bool mRandom;
	bool mTime;
	std::string mAlpha;
	QV_Options mQvalue;
	double mUpperhalfBinNumber;
	double mCovS;
	double mCisBinNumber;
	double mTransBinNumber;

public:
	SetupComp();
	SetupComp(std::string outDir, std::string enzyme, std::string input, int threads);

	std::vector<std::array<double,3>> mValues;

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
	inline Verb_level getVerbose() const {return mVerbose;}
	inline bool getRandom() const {return mRandom;}
	inline std::string getAlpha() const {return mAlpha;}
	inline QV_Options getQvalue() const {return mQvalue;}
	inline double getUpperhalfBinNumber() {return mUpperhalfBinNumber;}
	inline double getCovS() {return mCovS;}
	inline double getCisBinNumber() {return mCisBinNumber;}
	inline double getTransBinNumber() {return mTransBinNumber;}
	inline bool getTime() const {return mTime;}

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
	void setCisTrans(std::string input);
	//{mCisTrans = input;}
	inline void setRemoveDiagonal(bool L) {mRemoveDiagonal = L;}
	inline void setVerbose(Verb_level L) {mVerbose = L;}
	inline void setRandom(bool L) {mRandom = L;}
	inline void setAlpha(std::string L) {mAlpha = L;}
	void setQvalue(std::string L);
	//{mQvalue = L;}
	inline void setUpperhalfBinNumber(double L) {mUpperhalfBinNumber = L;}
	inline void setCovS(double L) {mCovS = L;}
	inline void setCisBinNumber(double L) {mCisBinNumber = L;}
	inline void setTransBinNumber(double L) {mTransBinNumber = L;}
	inline void setTime(bool L) {mTime = L;}

};

SetupComp loadConfigComp(std::string & fileName);
SetupComp setConfigComp(int argc, char * argv[]);

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
	sc_Random,
	sc_Alpha,
	sc_Qvalues,
	sc_Time
};



#endif /* SRC_SETUPCOMP_H_ */
