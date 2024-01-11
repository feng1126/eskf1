/*
 * KF-GINS: An EKF-Based GNSS/INS Integrated Navigation System
 *
 * Copyright (C) 2022 i2Nav Group, Wuhan University
 *
 *     Author : Liqiang Wang
 *    Contact : wlq@whu.edu.cn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#define _CRT_SECURE_NO_WARNINGS

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#include "readData.h"
#include "eskf.h"
#include "enu.h"
#include "matplotlibcpp.h"

const double D2R = (EIGEN_PI / 180.0);
const double R2D = (180.0 / EIGEN_PI);

namespace plt = matplotlibcpp;

int main(int argc, char* argv[])
{

    std::string file = "E:\\log\\0110_1606Road.log";
    readLog mSplitData;
    mSplitData.setFileName(file);
    mSplitData.parseLog();
    auto fiterData = mSplitData.filterData;

    auto mFilterQueue = mSplitData.mFilterQueue;

    FILE* fp = fopen("D:\\vsproject\\KF-GINS\\build\\out.log", "w");

    std::vector<double> pointX, pointY, times, heads, pointXEKF, pointYEKF, timesEKF, headEKF, pitchEKF;


    std::vector<double> MCUX, MCUY, MCUtimes, MCUheads;

    eskf EKF;
    ENU m_ENU;
    bool init = false;
    while (!mFilterQueue.empty())
    {
        filterMessage observe = mFilterQueue.top();
        //std::cout << observe.id << "     " << observe.utimestamp << std::endl;
        mFilterQueue.pop();
        if (observe.id == 0)
        {
            fprintf(fp, "in_test_loc:%lf,%0.7f,%0.7f,%0.7f,%0.7f\n", observe.utimestamp / 1000.0, observe.LLA[0], observe.LLA[1],
                observe.LLA[2], observe.YPR[0]);

            //printf("in_test_loc:%lf,%0.7f,%0.7f,%0.7f,%0.7f\n", observe.utimestamp / 1000.0, observe.LLA[0], observe.LLA[1],
              //  observe.LLA[2], observe.YPR[0]);

            if (!init)
            {
                m_ENU.setStart(observe.LLA[0], observe.LLA[1], observe.LLA[2]);
                double a, b, c;
                m_ENU.getENU(observe.LLA[0], observe.LLA[1], observe.LLA[2], a, b, c);
                observe.timestamp = observe.utimestamp / 1000.0;
                observe.UTM = Eigen::Vector3d(a, b, c);
                EKF.InitPose(observe);
                init = true;
            }
            else
            {
                double a, b, c;
                m_ENU.getENU(observe.LLA[0], observe.LLA[1], observe.LLA[2], a, b, c);
                observe.UTM = Eigen::Vector3d(a, b, c);
                observe.timestamp = observe.utimestamp / 1000.0;

                observe.std[0] = 1e-8;
                observe.std[1] = 1e-9;
                EKF.addGNSSPose(observe);


     /*           if (observe.utimestamp > 1030015)
                {
                    pointX.push_back(a);
                    pointY.push_back(b);
                }*/

                pointX.push_back(a);
                pointY.push_back(b);
                times.push_back(observe.timestamp);

                if (observe.YPR[0] < 0)
                    observe.YPR[0] = observe.YPR[0] + 2 * EIGEN_PI;
                heads.push_back(observe.YPR[0]* 57.3);
            }


        }

        if (init && observe.id == 1)
        {
            observe.timestamp = observe.utimestamp / 1000.0;
            observe.std[0] = 1e-6;
            //EKF.AddVelData(observe);
        }
        if (init && observe.id == 2)
        {
            observe.timestamp = observe.utimestamp / 1000.0;
            observe.acc = observe.acc *  0.1 * 9.8;
            observe.gyro = observe.gyro * EIGEN_PI / 180.0;

            Eigen::Vector3d IMUGbias , IMUAbias;
            IMUGbias << -0.005051, -0.002751, -0.009527;
            IMUAbias << -0.438778, 0.133060, 9.910083;

            observe.gyro = observe.gyro - IMUGbias - IMUGbias ;
            observe.acc[0] = observe.acc[0] - IMUAbias[0];
            observe.acc[1] = observe.acc[1] - IMUAbias[1];
            observe.acc[2] = observe.acc[1] + IMUAbias[1] - 9.8;


            EKF.addImuData(observe);

            EKF.ImuProcess();

        }

        if (init && observe.id == 0)
        {
            State state;
            EKF.getPose(state);


            //std::cout << state.ba.transpose() << " " << state.bg.transpose() << std::endl;

            double a, b, c;
            m_ENU.getLLA(state.p[0], state.p[1], 0, a, b, c);

            fprintf(fp, "in_test_ekf:%lf,%0.9f,%0.9f,%0.9f \n", state.time, a,b,c);

           // printf("in_test_ekf:%lf,%0.9f,%0.9f,%0.9f \n", state.time, a, b, c);


            double heading = ToEulerAngles(state.q)[2];

            if (heading > 2 * EIGEN_PI)
            {
                heading = heading - 2.0 * EIGEN_PI;
            }

            if (heading < 0)
                heading = heading + 2 * EIGEN_PI;

            pointXEKF.push_back(state.p[0]);
            pointYEKF.push_back(state.p[1]);
            timesEKF.push_back(state.time);
            headEKF.push_back(heading * 57.3);


            heading = ToEulerAngles(state.q)[0];

            if (heading > 2 * EIGEN_PI)
            {
                heading = heading - 2.0 * EIGEN_PI;
            }

            if (heading < 0)
                heading = heading + 2 * EIGEN_PI;

            pitchEKF.push_back(heading * 57.3);

        }
        if (init && observe.id == 5 )
        {
            double a, b, c;

          //  std::cout << std::to_string(observe.LLA[0]) << "     " << std::to_string(observe.LLA[1]) << "    " << observe.LLA[2] << std::endl;
            m_ENU.getENU(observe.LLA[0], observe.LLA[1], observe.LLA[2], a, b, c);


          // std::cout << a << "     " << b << "    " << c << std::endl;
           //if (observe.utimestamp < 1216474)
           //{

            MCUX.push_back(a);
            MCUY.push_back(b);
            MCUtimes.push_back(observe.utimestamp/1000.0 );
            MCUheads.push_back(observe.YPR[0] * 57.3);
           //}
        }


    }
    plt::figure_size(1280, 960);
    plt::grid(true);

    plt::set_aspect(1);

    plt::named_plot("EKF",pointXEKF, pointYEKF, "b-*");
    plt::named_plot("mcuEKF", MCUX, MCUY, "r-*");
    plt::named_plot("GNSS", pointX, pointY, "y-*");
    plt::legend();

    plt::figure_size(1280, 960);
    plt::grid(true);
    plt::named_plot("heading ", times, heads, "y*-");
    plt::named_plot("headingEKF ", timesEKF, headEKF, "b*-");
    //plt::named_plot("pitchEKF ", timesEKF, pitchEKF, "c*-");
    plt::named_plot("MCUEKF ", MCUtimes, MCUheads, "r*-");
    
    plt::legend();


    plt::show();


    return 0;
}


