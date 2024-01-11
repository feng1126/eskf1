

#include "readData.h"
#include <iostream>

readLog::readLog() {

    count = 0;
}

readLog::~readLog() {
}

std::vector<std::string> readLog::splitString(const std::string& strs) {
    std::string temp;
    std::vector<std::string> splitOut;
    splitOut.clear();

    for (int i = 0; i < strs.size(); i++) {
        if (strs[i] == ',') {
            splitOut.push_back(temp);
            temp = "";
        }
        else {
            temp = temp + strs[i];
        }

        if (i == strs.size() - 1) {
            splitOut.push_back(temp);
        }
    }
    return splitOut;
}

std::vector<std::string> readLog::splitStringwithSymbol(const std::string& strs, const std::string& symbol)
{
    std::string temp;
    std::vector<std::string> splitOut;
    splitOut.clear();

    for (int i = 0; i < strs.size(); i++) 
    {
        //std::cout << strs[i] << std::endl;
        if (strs[i] == ' ') 
        {


            splitOut.push_back(temp);
            temp = "";
        }
        else 
        {

            temp = temp + strs[i];

           // std::cout << temp << std::endl;
        }

        if (i == strs.size() - 1) {
            splitOut.push_back(temp);
        }
    }
    return splitOut;
}

std::vector<std::string> readLog::splitWithSymbol(const std::string& strs, const std::string& splitWithSymbol) {
    std::vector<std::string> strData;
    size_t pos = strs.find(splitWithSymbol, 0);
    std::string head = strs.substr(0, pos);
    if (pos == strs.npos)
        return std::vector<std::string>();
    strData.push_back(head);
    std::string tmpStr = strs.substr(pos + 1, strs.size());
    strData.push_back(tmpStr);
    return strData;
}

void readLog::parseLog() {
    std::string strs;

    filterData.clear();
    while (std::getline(in_fp, strs)) {
        std::vector<std::string> data = splitWithSymbol(strs, ":");
        if (data.size() != 2)
            continue;

        std::vector<std::string> logData = splitString(data[1]);
        if (data[0].compare("in_test_loc1") == 0) {
            if (stod(logData[5]) != 0) {
                uint64_t time = stod(logData[0]);
                filterMessage data;
                data.utimestamp = time;
                data.id = 0;
                data.LLA = Eigen::Vector3d(stod(logData[5]), stod(logData[6]), 0);
                data.YPR = Eigen::Vector3d(stod(logData[7]), 0, 0);// *EIGEN_PI / 180.0;
                mFilterQueue.push(data);
                filterData.insert(std::make_pair(time, data));
            }
        }

        if (data[0].compare("in_test_ekf") == 0) {

            if (stod(logData[0]) != 0) {
                uint64_t time = stod(logData[0]);
                filterMessage data;
                data.utimestamp = time;
                data.id = 5;
                data.LLA = Eigen::Vector3d(stod(logData[1]), stod(logData[2]), 0);
                data.YPR = Eigen::Vector3d(stod(logData[3]), 0, 0) * EIGEN_PI / 180.0;
                mFilterQueue.push(data);
                filterData.insert(std::make_pair(time, data));
            }
        }

        if (data[0].compare("in_test_vel") == 0)
        {
            uint64_t time = stod(logData[0]);
            filterMessage data;
            data.id = 1;
            data.utimestamp = stod(logData[0]);
            data.Vehicle = Eigen::Vector3d(stod(logData[1]) / 3.6, 0, 0);
            mFilterQueue.push(data);
            //filterData.insert(std::make_pair(time, data));
        }

        if (data[0].compare("in_test_IMU") == 0)
        {

            if (stod(logData[0]) != 0 && logData.size() > 6) 
            {
                uint64_t time = stod(logData[0]) ;
                filterMessage data;
                data.utimestamp = time;
                data.id = 2;
                data.acc = Eigen::Vector3d(stod(logData[1]), stod(logData[2]), stod(logData[3]));
                data.gyro = Eigen::Vector3d(stod(logData[4]), stod(logData[5]), stod(logData[6]));
                mFilterQueue.push(data);
                filterData.insert(std::make_pair(time, data));

                if (data.utimestamp < lastIMUTime)
                {
                    std::cout << "imu_______________error" << "   " << lastIMUTime << " " << data.utimestamp  <<  std::endl;
                }

                lastIMUTime = data.utimestamp;


            }
        }

        if (data[0].compare("In") == 0)
        {
            std::vector<std::string> logData = splitStringwithSymbol(data[1], " ");

            
           // std::cout << logData.size()<< std::endl;

            uint64_t time = stod(logData[3]);
            filterMessage data;
            data.utimestamp = time;
            data.id = 0;
            data.LLA = Eigen::Vector3d(stod(logData[7]), stod(logData[5]), 0);
            data.YPR = Eigen::Vector3d(stod(logData[11]), 0, 0);// *EIGEN_PI / 180.0;

            if (stod(logData[9]) != 0 && stod(logData[5]) > 0.0)
            {
                mFilterQueue.push(data);
                filterData.insert(std::make_pair(time, data));
            }

        }
    }
}

void readLog::setFileName(const std::string& logName) {
    in_fp.open(logName.c_str());
    if (!in_fp) {
        std::cout << "open fail" << std::endl;
    }
}
