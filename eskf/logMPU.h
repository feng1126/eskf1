#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include "data.h"

class logSplitMPU
{
public:
	logSplitMPU() {};
	~logSplitMPU() {};
	std::string GetHeadWithConut(const std::string& strs,int count);
	std::vector<std::string> splitString(const std::string& strs);
	void parseLog();
	void parseLog(std::string str);
	std::string splitWithColon(const std::string& strs);
	void setFileNameOut(std::string fileNameOut);
	void setFileName(std::string LogName);

	std::vector<Data> GNSSF9K;
	std::vector<Data> GNSSF9P;
	std::vector<Data> GNSSF9H;
	std::vector<Data> imuData;
	std::vector<Data> LaneData;
	std::vector<Data> VehicleData;
	std::vector<Data> ATTData;
	std::vector<Data> EKFData;

private:
	std::ifstream in_fp;
	std::ofstream out_GNSSF9K;
	std::ofstream out_GNSSF9P;
	std::ofstream out_GNSSF9H;
	//std::ofstream out_vel;
	//std::ofstream out_imu;
	//std::ofstream EKF_log;
	std::string mFileNameOut;
	std::vector<std::string> splitWithSymbol(const std::string& strs, const std::string& splitWithSymbol);
	void saveLog(std::vector<std::string> data, std::ofstream& ofs, std::string logName);



};

