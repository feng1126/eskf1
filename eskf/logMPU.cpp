#include "logMPU.h"


std::string logSplitMPU::GetHeadWithConut(const std::string& strs, int count)
{
	std::string head = strs.substr(0, count);
	return head;
}

std::vector<std::string> logSplitMPU::splitString(const std::string& strs)
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

std::vector<std::string> logSplitMPU::splitWithSymbol(const std::string& strs, const std::string& splitWithSymbol)
{
	std::vector<std::string> strData;
	size_t pos = strs.find(splitWithSymbol, 0);
	std::string head = strs.substr(0, pos);
	if (pos == strs.npos) return std::vector<std::string>();
	strData.push_back(head);
	std::string tmpStr = strs.substr(pos + 1, strs.size());
	strData.push_back(tmpStr);
	return strData;
}

void logSplitMPU::saveLog(std::vector<std::string> data, std::ofstream& ofs, std::string logName)
{
	if (!ofs)
	{
		std::string fileOut = mFileNameOut + "_" + logName + ".txt";
		ofs.open(fileOut.c_str());
	}

	for (int i = 0; i < data.size() - 1; i++)
	{
		ofs << data[i] << ",";
	}
	ofs << data[data.size() - 1];
	ofs << std::endl;

}


void logSplitMPU::parseLog()
{
	GNSSF9K.clear();
	GNSSF9P.clear();
	GNSSF9H.clear();
	imuData.clear();
	ATTData.clear();
	EKFData.clear();
	std::string strs;
	while (std::getline(in_fp, strs))
	{
		std::vector<std::string>  data = splitWithSymbol(strs, ":");
		if (data.size() < 2) continue;
		std::vector<std::string> head = splitString(data[0]);
		std::vector<std::string> logData = splitString(data[1]);

		if (data[0].compare("GNSS_F9P") == 0)
		{
			if (logData.size() == 14)
			{
				double lat, lon, alt, heading;
				lat = stod(logData[5]);
				lon = stod(logData[6]);
				alt = stod(logData[7]);
				heading = stod(logData[8]);
				heading = heading + 0.5 * EIGEN_PI;
				if (heading > 2 * EIGEN_PI)
				{
					heading = heading - 2.0 * EIGEN_PI;
				}
				Data data;
				data.id = 0;
				data.timeStamp = stod(logData[0]);
				data.timeGPS = stod(logData[1]);
				data.lla = Eigen::Vector3d(lat, lon, 0);
				data.ypr = Eigen::Vector3d(heading, 0, 0);
				data.v = Eigen::Vector3d(stod(logData[9]), 0, 0);
				GNSSF9P.push_back(data);

			}
		}
		else if (head[0].compare("ekf_in_gnss") == 0)
		{
			if (logData.size() == 18 )
			{
				double lat, lon, alt, heading;
				lat = stod(logData[6]);
				lon = stod(logData[7]);
				alt = stod(logData[8]);
				heading = stod(logData[9]);
				heading = heading + 0.5 * EIGEN_PI;
				if (heading > 2 * EIGEN_PI)
				{
					heading = heading - 2.0 * EIGEN_PI;
				}

				//heading =   2.0 * EIGEN_PI - heading;
				Data data;
				data.id = 0;
				data.timeStamp = stod(logData[0]);
				data.timeGPS = stod(logData[1]);
				data.lla = Eigen::Vector3d(lat, lon, 0);

				//std::cout << logData[10] << " , " << logData[11] << std::endl;
				data.ypr = Eigen::Vector3d(heading, stod(logData[10]), stod(logData[11]));
				//data.v = Eigen::Vector3d(stod(logData[12]), stod(logData[13]), stod(logData[14]));
				GNSSF9K.push_back(data);
			}

		}
		else if (head[0].compare("F9H") == 0)
		{
			if (logData.size() == 9)
			{
				double timestamp, lat, lon, alt, heading, vel;
				timestamp = stod(logData[0]);
				heading = stod(logData[7]);
				heading = heading + 0.5 * EIGEN_PI;
				if (heading > 2 * EIGEN_PI)
				{
					heading = heading - 2.0 * EIGEN_PI;
				}
				Data data;
				data.id = 0;
				data.timeGPS = timestamp;
				data.ypr = Eigen::Vector3d(heading, 0, 0);
				GNSSF9H.push_back(data);
			}
		}
		else if (head[0].compare("ekf_out") == 0)
		{
			if (logData.size() == 6)
			{
				double timestamp, lat, lon, alt, heading, vechile;
				timestamp = stod(logData[0]);
				lat = stod(logData[1]);
				lon = stod(logData[2]);
				heading = stod(logData[3]);
				vechile = stod(logData[4]);
				Data data;
				data.timeStamp = timestamp;
				data.lla = Eigen::Vector3d(lat, lon, 0);
				data.ypr = Eigen::Vector3d(heading, 0, 0);
				data.v = Eigen::Vector3d(vechile, 0, 0);
				EKFData.push_back(data);

			}

		}
		else if (head[0].compare("ekf_in_imu") == 0)
		{
			double  ax, ay, az, gx, gy, gz;
			ax = stod(logData[1]);
			ay = stod(logData[2]);
			az = stod(logData[3]);
			gx = (stod(logData[4])) * PI / 180.0;
			gy = (stod(logData[5])) * PI / 180.0;
			gz = (stod(logData[6])) * PI / 180.0;
			Data data;
			data.id = 1;
			data.timeStamp = stod(logData[0]);
			data.acc = Eigen::Vector3d(ax, ay, az);
			data.gyro = Eigen::Vector3d(gx, -gy, gz);
			imuData.push_back(data);
		}
		else if (head[0].compare("ekf_in_lane") == 0)
		{
			double lat, lon, alt, heading;
			lat = stod(logData[2]);
			lon = stod(logData[3]);
			alt = 0;
			heading = stod(logData[4 ]);
			if (heading == 0.00) continue;
			if (heading > 7.00) continue;
			//if (heading > 2 * EIGEN_PI) continue;
			heading = heading + 0.5 * EIGEN_PI + 1.5 * EIGEN_PI / 180.0;
			if (heading > 2 * EIGEN_PI)
			{
				heading = heading - 2.0 * EIGEN_PI;
			}
			Data data;
			data.id = 2;
			data.timeStamp = stod(logData[0]) + 0.000001;
			data.lla = Eigen::Vector3d(lat, lon, 0);
			data.ypr = Eigen::Vector3d(heading, 0, 0);
			data.cov[0] = stod(logData[5]);
			LaneData.push_back(data);

		}
		else if (head[0].compare("ekf_in_vel") == 0)
		{
			double Vehicle;
			Vehicle = stod(logData[1]);
			Data data;
			data.id = 3;
			data.timeStamp = stod(logData[0]);
			data.v = Eigen::Vector3d(Vehicle, 0, 0);
			VehicleData.push_back(data);
		}
		else if (head[0].compare("ekf_in_ATT") == 0)
		{

			Data data;
			data.id = 4;

			data.timeStamp = stod(logData[0]);

			double yaw = stod(logData[2]);
			double pitch = stod(logData[3]);
			double roll = stod(logData[4]);



			yaw = yaw + 0.5 * EIGEN_PI;
			if (yaw > 2 * EIGEN_PI)
			{
				yaw = yaw - 2.0 * EIGEN_PI;
			}

			//pitch = pitch +0.5 * EIGEN_PI;
	/*		if (pitch > 2 * EIGEN_PI)
			{
				pitch = pitch - 2.0 * EIGEN_PI;
			}*/

			//roll = roll + 0.5 * EIGEN_PI;
	/*		if (roll > 2 * EIGEN_PI)
			{
				roll = roll - 2.0 * EIGEN_PI;
			}*/

			data.ypr = Eigen::Vector3d(yaw, pitch, roll);
			ATTData.push_back(data);
		}
	}
	in_fp.close();
}


void logSplitMPU::parseLog(std::string str)
{
	std::string fileOut = mFileNameOut + "_" + str + ".txt";
	std::ofstream ofs;
	ofs.open(fileOut.c_str());
	std::string strs;
	while (std::getline(in_fp, strs))
	{
		std::string head = GetHeadWithConut(strs, str.size());
		if (head.compare(str.c_str()) == 0)
		{
			ofs << strs << std::endl;
		}
	}
	in_fp.close();
	ofs.close();
}

std::string logSplitMPU::splitWithColon(const std::string& strs)
{
	size_t pos = strs.find(":", 0);
	std::string head = strs.substr(0, pos);
	return head;
}


void logSplitMPU::setFileName(std::string LogName)
{
	in_fp.open(LogName.c_str());
}

void logSplitMPU::setFileNameOut(std::string fileNameOut)
{
	mFileNameOut = fileNameOut;
}
