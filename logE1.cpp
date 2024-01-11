#include "logE1.h"

std::string logSplitE1::GetHeadWithConut(const std::string& strs, int count)
{
	std::string head = strs.substr(0, count);
	return head;
}

std::vector<std::string> logSplitE1::splitString(const std::string& strs)
{
	std::string temp;
	std::vector<std::string> splitOut;
	splitOut.clear();
	for (int i = 0; i < strs.size(); ++i)
	{
		if (strs[i] == ',')
		{
			splitOut.push_back(temp);
			temp = "";
		}
		else
		{
			temp = temp + strs[i];
		}

		if (i == strs.size() - 1)
		{
			splitOut.push_back(temp);
		}
	}
	return splitOut;

}

std::vector<std::string> logSplitE1::splitWithSymbol(const std::string& strs,const std::string& splitWithSymbol)
{
	std::vector<std::string> strData;
	size_t pos = strs.find(splitWithSymbol, 0);
	std::string head = strs.substr(0, pos);
	strData.push_back(head);
	std::string tmpStr = strs.substr(pos + 1, strs.size());
	strData.push_back(tmpStr);
	return strData;
}


void logSplitE1::parseLogUTM()
{
	std::string strs;
	while (std::getline(in_fp, strs))
	{

		std::string head = GetHeadWithConut(strs, 9);
		if (head.compare("#INSPVAXA") == 0)
		{
			std::vector<std::string>  data = splitWithSymbol(strs, ";");

			std::vector<std::string> timeSatmp = splitString(data[0]);
			std::vector<std::string> GNSSData = splitString(data[1]);
			if (GNSSData[1].compare("NONE") != 0 && GNSSData.size() == 23)
			{
				//double utmX, utmY;
				//tool_t::UTMTransform::instance()->LLToUTM(stod(GNSSData[2]), stod(GNSSData[3]), utmX, utmY);
				double  lat = stod(GNSSData[2]);
				double lon = stod(GNSSData[3]);
				double alt = stod(GNSSData[4]);
				double heading = stod(GNSSData[11]);

				Data data;
				data.id = 0;
				data.timeStamp = stod(timeSatmp[6]);
				data.lla = Eigen::Vector3d(lat, lon, alt);
				data.ypr = Eigen::Vector3d(heading, 0, 0);
				GNSSE1.push_back(data);
	

			}

		}
	}

	in_fp.close();
	out_log.close();
}


void logSplitE1::parseLog()
{
	GNSSE1.clear();

	std::string strs;
	while (std::getline(in_fp, strs))
	{

		//std::cout << strs << std::endl;
		std::string head = GetHeadWithConut(strs, 15);
		if(head.compare("[COM3]#INSPVAXA") == 0)
		{
			std::vector<std::string>  data = splitWithSymbol(strs, ";" );

			std::vector<std::string> timeSatmp = splitString(data[0]);
			std::vector<std::string> GNSSData = splitString(data[1]);
			if (GNSSData[1].compare("NONE") != 0 && GNSSData.size() == 23)
			{
				double timestamp = stod(timeSatmp[6]);
				double lat = stod(GNSSData[2]);
				double lon = stod(GNSSData[3]);
				double alt = stod(GNSSData[4]);
				double heading = stod(GNSSData[11]) * EIGEN_PI / 180.0;// +0.5 * EIGEN_PI;
				heading = 2 * EIGEN_PI - heading;
				heading = heading + 0.5 * EIGEN_PI -  1 * EIGEN_PI / 180.0;
				if (heading > 2 * EIGEN_PI) heading = heading - 2 * EIGEN_PI;
				Data data;
				data.id = 0;
				data.timeGPS = timestamp;
				data.lla = Eigen::Vector3d(lat, lon, alt);
				data.ypr = Eigen::Vector3d(heading, 0, 0);
 				GNSSE1.push_back(data);			
			}
		}
	}

	in_fp.close();
	out_log.close();
}


void logSplitE1::setFileNameOut(const std::string& fileNameOut)
{
	out_log.open(fileNameOut.c_str());
}

void logSplitE1::setFileName(const std::string& logName)
{
	in_fp.open(logName.c_str());
}

