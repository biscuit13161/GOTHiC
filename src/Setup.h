/*
 * Utils.h
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

#ifndef SRC_SETUP_H_
#define SRC_SETUP_H_

#include <string>
#include <map>

class Setup
{
private:
	std::string mOutDir;
	std::string mEnzyme;
	std::string mInput;
	int mThreads;
	int mRes;
public:
	Setup();
	Setup(std::string outDir, std::string enzyme, std::string input, int threads);

	void print();

	inline int getThreads() const {return mThreads;}
	inline int getRes() const {return mRes;}
	inline std::string getEnzyme() const {return mEnzyme;}
	inline std::string getOutDir() const {return mOutDir;}
	inline std::string getInput() const {return mInput;}

	inline void setThreads(int threads) {mThreads = threads;}
	inline void setRes(int res) {mRes = res;}
	inline void setEnzyme(std::string enzyme) {mEnzyme = enzyme;}
	inline void setOutDir(std::string outDir) {mOutDir = outDir;}
	inline void setInput(std::string input) {mInput = input;}
};

Setup loadConfig(std::string & fileName);

enum Options {
	Input = 0,
	Digest,
	Threads,
	Res,
	Output,
};


#endif /* SRC_SETUP_H_ */
