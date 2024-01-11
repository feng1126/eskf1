#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include "data.h"


class logSplitE1
{
public:
	logSplitE1() {};
	~logSplitE1() {};

	std::vector<std::string> splitString(const std::string& strs);
	void parseLog();
	void parseLogUTM();
	void setFileName(const std::string& logName);
	void setFileNameOut(const std::string& fileNameOut);

	std::vector<Data> GNSSE1;


private:
	std::ifstream in_fp;
	std::ofstream out_log;
	std::vector<std::string> splitWithSymbol(const std::string& strs, const std::string& splitWithSymbol);
	std::string GetHeadWithConut(const std::string& strs, int count);

};

