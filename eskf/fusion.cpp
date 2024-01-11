#include <vector>
#include "Eigen/Core"
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include "Eigen/Geometry"
#include <deque>
#include <random>
#include "eskf.h"
#include "enu.h"
#include "matplotlibcpp.h"

#include "logMPU.h"
#include "logE1.h"
#include "SystemTime.h"




State state_;


namespace plt = matplotlibcpp;
int main(int argc, char** argv)
{
    logSplitMPU mSplitMPU;
    //std::string fileStrMPU = "C:\\Users\\niew\\Desktop\\data\\12\\2021-12-1-15-48-22.txt";
    std::string fileStrMPU = "C:\\Users\\niew\\Desktop\\data\\ATT\\2022-7-1-14-2-24.txt";
    //std::string fileStrE1 = "C:\\Users\\niew\\Desktop\\data\\12\\ReceivedTofile-COM9-2021_12_3_17-06-02.DAT";

    mSplitMPU.setFileName(fileStrMPU);
    mSplitMPU.parseLog();

    std::map<double, Data> GNSSF9K;
    std::map<double, Data> GNSSF9P;
    std::map<double, Data> GNSSF9H;
    std::map<double, Data> GNSSE1;


    namespace plt = matplotlibcpp;
    std::vector<double>  VecGNSSF9K, VecGNSSF9P, VecGNSSF9H, VecGNSSE1, timeGNSSE1, timeGNSSF9K, timeGNSSF9P, timeGNSSF9H;
    std::vector<double> EKFTime, EKFYaw, GNSSTime,GNSSYaw, E1Yaw, GNSSF9KX, GNSSF9PX, GNSSE1X, GNSSF9KY, GNSSF9PY, GNSSE1Y, GNSSEKFX, GNSSEKFY, yawESKF, timeESKF, ekfX, ekfY;

    std::cout << "f9k " << mSplitMPU.GNSSF9K.size() << std::endl;



    std::map<double, Data> ESKF;

    for (int i = 0; i < mSplitMPU.GNSSF9K.size() ; i++)
    {
        ESKF.insert(std::make_pair(mSplitMPU.GNSSF9K[i].timeStamp, mSplitMPU.GNSSF9K[i]));
        GNSSTime.push_back(mSplitMPU.GNSSF9K[i].timeStamp);
        GNSSYaw.push_back(mSplitMPU.GNSSF9K[i].ypr[0] * 57.3);
    }
    std::cout << "imu size: " << mSplitMPU.imuData.size() << std::endl;
    std::cout << "LaneData size: " << mSplitMPU.LaneData.size() << std::endl;
    std::cout << "gnss size: " << mSplitMPU.GNSSF9K.size() << std::endl;

    for (int i = 0; i < mSplitMPU.imuData.size(); i++)
    {
        ESKF.insert(std::make_pair(mSplitMPU.imuData[i].timeStamp, mSplitMPU.imuData[i]));
    }

    for (int i = 0; i < mSplitMPU.LaneData.size(); i++)
    {
        ESKF.insert(std::make_pair(mSplitMPU.LaneData[i].timeStamp, mSplitMPU.LaneData[i]));
    }


    for (int i = 0; i < mSplitMPU.VehicleData.size(); i++)
    {
        ESKF.insert(std::make_pair(mSplitMPU.VehicleData[i].timeStamp, mSplitMPU.VehicleData[i]));
    }
    

    eskf eskf;
    ENU m_ENU;
    bool init = false;


    int count = 0;
    for (auto iter = ESKF.begin(); iter != ESKF.end(); iter++)
    {

        double timeStart = GetSystemTimeStamp();
        if (!init)
        {
            if (iter->second.id == 1) continue;
            if (iter->second.id == 0)
            {
                m_ENU.setStart(iter->second.lla[0], iter->second.lla[1], iter->second.lla[2]);
                double a, b, c;
                m_ENU.getENU(iter->second.lla[0], iter->second.lla[1], iter->second.lla[2], a, b, c);
                state_.time = iter->second.timeStamp;
                state_.v = Eigen::Vector3d(0, 0, 0);
                state_.p = Eigen::Vector3d(a,b,c);
                state_.q =  YPR2Quaterniond(iter->second.ypr[2], iter->second.ypr[1], iter->second.ypr[0]);
                state_.bg = Eigen::Vector3d(-0.0215333, 0.00954956, 0.00105915);
                state_.ba = Eigen::Vector3d(0.123256, 0.121496, 0.0442121);
                filterMessagePtr observe = std::make_shared<filterMessage>();
                observe->timestamp = iter->second.timeStamp;
                observe->YPR = iter->second.ypr;
                observe->UTM = state_.p;
                observe->LLA = iter->second.lla;

                eskf.InitPose(observe);
                init = true;
            }
        }
        else if (iter->second.id == 0 )
        {
            double a, b, c;
            m_ENU.getENU(iter->second.lla[0], iter->second.lla[1], iter->second.lla[2], a, b, c);
            Eigen::Vector3d gnss_xyz(a, b, c);
            iter->second.xyz = gnss_xyz;
            Eigen::VectorXd cov;
            cov.resize(6, 1);
            cov[0] = 9e-9;
            cov[1] = 1e-9;

            filterMessagePtr observe = std::make_shared<filterMessage>();
            observe->timestamp = iter->second.timeStamp;
            observe->YPR = iter->second.ypr;
            observe->UTM = gnss_xyz;
            observe->LLA = iter->second.lla;

            //observe->Vehicle= iter->second.v;

            eskf.updatePose(observe, cov);

            //std::cout << observe->Vehicle << std::endl;

            //eskf.updateVehicle(observe, 1e-8);

            GNSSF9KX.push_back(a);
            GNSSF9KY.push_back(b);

        }
        else if (iter->second.id == 1)
        {
            filterMessagePtr observe = std::make_shared<filterMessage>();
            observe->timestamp = iter->second.timeStamp;
            observe->acc = iter->second.acc;
            observe->gyro = iter->second.gyro;
            eskf.Predict(observe);
        }
        else if (iter->second.id == 2)
        {
            if (iter->second.cov[0] > 0.4) continue;
            double a, b, c;
            m_ENU.getENU(iter->second.lla[0], iter->second.lla[1], iter->second.lla[2], a, b, c);
            Eigen::Vector3d gnss_xyz(a, b, c);
            iter->second.xyz = gnss_xyz;

            filterMessagePtr observe = std::make_shared<filterMessage>();
            observe->timestamp = iter->second.timeStamp;
            observe->UTM = gnss_xyz;
            observe->LLA = iter->second.lla;

            eskf.updateXYZ(observe, 1e-10);
            //eskf.updateRPY(iter->second);
        }
        else if (iter->second.id == 3)
        {
            filterMessagePtr observe = std::make_shared<filterMessage>();
            observe->timestamp = iter->second.timeStamp;
            observe->Vehicle = iter->second.v;
            eskf.updateVehicle(observe,1e-9);
        }
        count++;

        if (!init) continue;

		if (iter->second.id == 0)
		{

            double end = GetSystemTimeStamp();

            //std::cout << "time run: " << std::to_string((end - timeStart) * 1000) << std::endl;

			State state;
			eskf.getPose(state);
			GNSSEKFX.push_back(state.p(0));
			GNSSEKFY.push_back(state.p(1));
			//std::cout << state.p(0) << "," << state.p(1) << "," << state.p(2)
			//      << "," << state.bg(0)
			//      << "," << state.bg(1)
			//      << "," << state.bg(2)
			//      << "," << state.ba(0)
			//      << "," << state.ba(1)
			//      << "," << state.ba(2)
			//      << std::endl;
			timeESKF.push_back(state.time);
			double heading = ToEulerAngles(state.q)[2];
			heading = heading + 2.0 * EIGEN_PI;
			if (heading > 2 * EIGEN_PI)
			{
				heading = heading - 2.0 * EIGEN_PI;
			}

            heading = ToEulerAngles(state.q)[2] + 0.5 * EIGEN_PI;

            heading = heading + 1.5 * EIGEN_PI;
            if (heading > 2 * EIGEN_PI)
            {
                heading = heading - 2.0 * EIGEN_PI;
            }

            if (heading < 0)
                heading = heading + 2 * EIGEN_PI;

			yawESKF.push_back(heading * 57.3);
		}

    }

    std::vector<double>  laneX, laneY, timeGNSSLane, GNSSLane ,laneYaw;
    for (int i = 0; i < mSplitMPU.LaneData.size(); i++)
    {
        double a, b, c;
        m_ENU.getENU(mSplitMPU.LaneData[i].lla[0], mSplitMPU.LaneData[i].lla[1], mSplitMPU.LaneData[i].lla[2], a, b, c);
        Eigen::Vector3d gnss_xyz(a, b, c);
        laneX.push_back(a);
        laneY.push_back(b);
        timeGNSSLane.push_back(mSplitMPU.LaneData[i].timeStamp);
        laneYaw.push_back(mSplitMPU.LaneData[i].ypr[0] * 57.3 );
    }

    std::vector<double>  ekfoutX, ekfoutY, timeEKF;
    for (int i = 0; i < mSplitMPU.EKFData.size(); i++)
    {
        double a, b, c;
        m_ENU.getENU(mSplitMPU.EKFData[i].lla[0], mSplitMPU.EKFData[i].lla[1], mSplitMPU.EKFData[i].lla[2], a, b, c);
        ekfoutX.push_back(a);
        ekfoutY.push_back(b);
        timeEKF.push_back(mSplitMPU.EKFData[i].timeStamp);
        //laneYaw.push_back(mSplitMPU.LaneData[i].ypr[0] * 57.3);
    }

 
   
    plt::figure_size(1200, 960);
    plt::grid(true);
    plt::named_plot("GNSS", GNSSF9KX, GNSSF9KY, "r*");
    plt::named_plot("ESKF", GNSSEKFX, GNSSEKFY, "b*");
    plt::named_plot("lane", laneX, laneY, "k*");
   // plt::named_plot("ekfout", ekfoutX, ekfoutY, "g*");
    plt::legend();


    plt::figure_size(1200, 960);
    plt::grid(true);
    plt::named_plot("ESKF", timeESKF, yawESKF, "b*");
    plt::named_plot("GNSS", GNSSTime, GNSSYaw, "r*");
    plt::named_plot("lane", timeGNSSLane, laneYaw, "k*");
    plt::legend();

    plt::show();

    //std::ifstream m_file_in;
    //std::string LogFile = "C:\\Users\\niew\\Desktop\\ublox\\11.1\\3.txt";
    //m_file_in.open(LogFile.c_str());
    //std::string strs;

    //std::vector<double> x1, y1, x2, y2, x3,y3, timeGNSS,timeLane,timeEKF, yawVec, yawVec1, yawVec2;
    //std::map<double, Data> dataMap;

    //while (std::getline(m_file_in, strs))
    //{
    //    std::vector<std::string> data = splitWithSymbol(strs, ":");
    //    if (data.size() < 2)
    //        continue;
    //    std::vector<std::string> logData = splitString(data[1]);

    //    if (data[0].compare("GNSS_F9P") == 0 && logData.size() == 14)
    //    {
    //        double timestamp, lat, lon, alt, heading, utmX, utmY, utmZ;
    //        timestamp = stod(logData[0]);
    //        lat = stod(logData[5]);
    //        lon = stod(logData[6]);
    //        alt = stod(logData[7]);
    //        heading = stod(logData[8]);
    //        
    //        heading = heading + 0.5 * EIGEN_PI -1.5 * EIGEN_PI / 180.0;
    //        if (heading > 2 * EIGEN_PI) heading = heading - 2 * EIGEN_PI;

    //        yawVec.push_back(heading * 57.3);
    //        timeGNSS.push_back(timestamp);
    //        Data data;
    //        data.id = 0;
    //        data.time = timestamp;
    //        data.lla = Eigen::Vector3d(lat, lon, alt);
    //        data.ypr = Eigen::Vector3d(heading, 0, 0);

    //        data.q = YPR2Quaterniond(0, 0, heading);// Eigen::Quaterniond(qw, qx, qy, qz);
    //        dataMap.insert(std::pair<double, Data>(timestamp, data));
    //    }
    //    else if (data[0].compare("ekf_in_imu") == 0 && logData.size() == 7)
    //    {
    //        double timestamp, ax, ay, az, gx, gy, gz;
    //        timestamp = stod(logData[0]);
    //        ax =  stod(logData[1]) ;
    //        ay =  stod(logData[2]);
    //        az =  stod(logData[3]);
    //        gx = (stod(logData[4])) *PI / 180.0;
    //        gy = (stod(logData[5])) *PI / 180.0;
    //        gz = (stod(logData[6])) *PI / 180.0;
    //        Data data;
    //        data.id = 1;
    //        data.time = timestamp;
    //        data.acc = Eigen::Vector3d(ax, ay, az);
    //        data.gyro = Eigen::Vector3d(gx, gy, gz);


    //        //std::cout << "acc " << data.acc << "gyro" << data.gyro << std::endl;
    //        dataMap.insert(std::pair<double, Data>(timestamp, data));
   
    //    }

    //    else if (data[0].compare("ekf_in_lane") == 0 && logData.size() == 4)
    //    {
    //        double timestamp, lat, lon, alt, heading, utmX, utmY, utmZ;
    //        timestamp = stod(logData[0]);
    //        lat = stod(logData[1]);
    //        lon = stod(logData[2]);
    //        heading = stod(logData[3]);
    //        heading = heading + 0.5 * EIGEN_PI;// +1.5 * EIGEN_PI / 180.0;
    //        timeLane.push_back(timestamp);
    //        yawVec1.push_back(heading *57.3);
    //        Data data;
    //        data.id = 2;
    //        data.time = timestamp;
    //        data.lla = Eigen::Vector3d(lat, lon, alt);
    //        data.ypr = Eigen::Vector3d(heading, 0, 0);
    //        data.q = YPR2Quaterniond(0, 0, heading);// Eigen::Quaterniond(qw, qx, qy, qz);
    //        dataMap.insert(std::pair<double, Data>(timestamp, data));
    //    }

    //    else
    //    {
    //        continue;
    //    }
    //}
    //bool init = false;
    //int count = 10;
    //std::vector<IMUData> imu_buffer_;
    //ENU m_ENU;
    //eskf ekf;

    //int unm = 0;
    //auto iter = dataMap.begin();
    //for ( ; iter != dataMap.end(); iter++)
    //{

    //    if (!init)
    //    {
    //        if (iter->second.id == 1)
    //        {
    //            IMUData imu;
    //            imu.time = iter->second.time;
    //            imu.acc = iter->second.acc;
    //            imu.gyro = iter->second.gyro;
    //            imu_buffer_.push_back(imu);
    //        }
    //        else if (iter->second.id == 0)
    //        {
    //            m_ENU.setStart(iter->second.lla[0], iter->second.lla[1], iter->second.lla[2]);
    //            double a, b, c;
    //            m_ENU.getENU(iter->second.lla[0], iter->second.lla[1], iter->second.lla[2], a, b, c);
    //            Eigen::Vector3d gnss_xyz(a, b, c);
    //           Eigen::Matrix3d G_R_I;

    //            //std::cout << "imu_buffer_ " << imu_buffer_.size() << std::endl;
    //            //ComputeG_R_IFromImuData(G_R_I, imu_buffer_);

    //            //std::cout << G_R_I << std::endl;

    //            state_.time = iter->second.time;
    //            state_.v = Eigen::Vector3d(0, 0, 0);
    //            state_.p = gnss_xyz;
    //            state_.q = iter->second.q;// YPR2Quaterniond(iter->second.ypr[2], iter->second.ypr[1], iter->second.ypr[0]);
    //            state_.bg = Eigen::Vector3d(0, 0, 0);
    //            state_.ba = Eigen::Vector3d(0, 0, 0);
    //            ekf.InitPose(state_, iter->second);
    //            init = true;
    //        }


    //    }
    //    else
    //    {
    //        if (iter->second.id == 0)
    //        {
    //            double a, b, c;
    //            m_ENU.getENU(iter->second.lla[0], iter->second.lla[1], iter->second.lla[2], a, b, c);
    //            Eigen::Vector3d gnss_xyz(a, b, c);
    //            //GNSSData current_gnss;
    //            //current_gnss.time = iter->second.time;
    //            //current_gnss.xyz = gnss_xyz;
    //            //current_gnss.lla = iter->second.lla;
    //            iter->second.xyz = gnss_xyz;
    //            ekf.updateGPS(iter->second);
    //            ekf.updateq(iter->second);


    //            x1.push_back(a);
    //            y1.push_back(b);



    //        }
    //        else if (iter->second.id == 1)
    //        {
    //            ekf.Predict(iter->second);
    //        }
    //        else if (iter->second.id == 2)
    //        {
    //            double a, b, c;
    //            m_ENU.getENU(iter->second.lla[0], iter->second.lla[1], iter->second.lla[2], a, b, c);
    //            Eigen::Vector3d gnss_xyz(a, b, c);
    //     /*       GNSSData current_gnss;
    //            current_gnss.time = iter->second.time;
    //            current_gnss.xyz = gnss_xyz;*/
    //            iter->second.xyz = gnss_xyz;
    //            //ekf.updatelane(iter->second);
    //            //ekf.updateq(iter->second);
    //            x2.push_back(a);
    //            y2.push_back(b);

    //            State state;
    //            ekf.getPose(state);
    //            x3.push_back(state.p(0));
    //            y3.push_back(state.p(1));
    //        }

    //        State state;
    //        ekf.getPose(state);

    //        x3.push_back(state.p(0));
    //        y3.push_back(state.p(1));
    //        timeEKF.push_back(state.time);
    //        //double yaw = ToEulerAngles(state.q)[2];

    //        //yawVec2.push_back(yaw  * 57.3);
    //    }



    //    //yawVec.push_back(ToEulerAngles(state.q)[2]);
    //    State state;
    //    ekf.getPose(state);
    //    //std::cout << state.p(0) << "," << state.p(1) << "," << state.p(2)
    //    //          << "," << state.bg(0)
    //    //          << "," << state.bg(1)
    //    //          << "," << state.bg(2)
    //    //          << "," << state.ba(0)
    //    //          << "," << state.ba(1)
    //    //          << "," << state.ba(2)
    //    //          << std::endl;

    //}

    ////int mm = 3000;
    //// x2.resize(mm);
    ////y2.resize(mm);
    ////yawVec.resize(mm);
    ////yawVec1.resize(mm);
    ////std::cout << yawVec.size() << std::endl;
    ////std::cout << x1.size() << std::endl;
    //for (int i = 0; i < yawVec.size(); i++)
    //{
    //    //printf("1: %lf \n", yawVec1[i] * 180 / PI);
    //    //printf("2: %lf \n", yawVec[i] * 180 / PI);

    //}

    //plt::figure_size(1200, 780);
    //// Plot line from given x and y data. Color is selected automatically.
    ////plt::plot(x, y);
    //// Plot a red dashed line from given x and y data.
    ////x1.resize(2000);
    ////y1.resize(2000);
    //plt::grid(true);
    //plt::plot(x3, y3, "k*");
    //plt::plot(x2, y2, "b*");
    //plt::plot(x1, y1, "r*");

    //plt::figure_size(1200, 780);
    //plt::grid(true);
    //plt::plot(timeGNSS, yawVec, "k-");
    //plt::plot(timeLane, yawVec1, "b-");
    //plt::plot(timeEKF,yawVec2, "r-");
    //plt::show();
    return 0;
}
