#pragma once

#include <vector>
#include <iostream>
#include <Eigen/Geometry>
#include "data.h"
#include "so3.h"

class eskf
{

public:
	bool Predict(const filterMessagePtr& imuData);
	void updatePose(const filterMessagePtr& pose, const Eigen::VectorXd cov);

	void updateGNSS(const filterMessagePtr& pose, const Eigen::VectorXd cov);

	void updateXYZ(const filterMessagePtr& XYZ, const double cov);
	void updateGPS(const filterMessagePtr& gnss, const double cov);
	void updateRPY(const filterMessagePtr& RPY, const double cov);
	void updatelane(const filterMessagePtr& gnss, const double cov);
	void updateVehicle(const filterMessagePtr& vehicle, const double cov);
	bool InitPose(const filterMessagePtr& data);
	void getPose(State& state);
	eskf();
private:

	State state_;
	ErrorState error_state;
	Eigen::Matrix3d ComputeEarthTranform(const double& dt, const Eigen::Vector3d& velocity);
	void predictNewState(const double& dt,const Eigen::Vector3d& gyro,const Eigen::Vector3d& acc);
	void predictNewRkState(const double& dt, const  Eigen::Vector3d& gyro, const  Eigen::Vector3d& acc);
	filterMessagePtr lastDataPtr;
	bool has_imu;
	bool IsCovStable(int INDEX_OFSET);
	Eigen::Vector3d LastP;




};