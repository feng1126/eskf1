#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <queue>
#include "data.h"
#include <unordered_map>
#include <Eigen/Geometry>
#include "eskf.h"

class  readLog
{
public:
    readLog();
    ~readLog();
    std::vector<std::string> splitStringwithSymbol(const std::string& strs, const std::string& symbol);
    std::vector<std::string> splitString(const std::string& strs);
    void parseLog();
    void setFileName(const std::string& logName);
    std::vector<std::string> splitWithSymbol(const std::string& strs, const std::string& splitWithSymbol);

    std::map<double, filterMessage> filterData;

    std::priority_queue<filterMessage, std::vector<filterMessage>, filterMessage> mFilterQueue;

private:
    std::ifstream in_fp;

    int count;

    uint64_t lastIMUTime;


};

