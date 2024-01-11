#include "eskf.h"

constexpr double kDegreeToRadian = DEGREE_TO_RAD;

eskf::eskf()
{

	error_state.x.setZero();
	error_state.p.setZero();
	error_state.p.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity() * 1e-4;
	error_state.p.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity() * 1e-4;
	error_state.p.block<3, 3>(6, 6) = Eigen::Matrix3d::Identity() * 1e-4;
	error_state.p.block<3, 3>(9, 9) = Eigen::Matrix3d::Identity() * 1e-4;
	error_state.p.block<3, 3>(12, 12) = Eigen::Matrix3d::Identity() * 1e-4;
	lastDataPtr = std::make_shared<filterMessage>();
	has_imu = false;
}

void eskf::getPose(State& state)
{

	state = state_;
}

bool eskf::InitPose(const filterMessage& Pose)
{
	state_.time = Pose.timestamp;
	state_.p = Pose.UTM;
	state_.q = YPR2Quaterniond(Pose.YPR);
	state_.v = state_.q * Pose.Vehicle;
	state_.ba.setZero();// = Eigen::Vector3d(0.0666269, -0.00173516, -0.0180815);
	state_.bg.setZero();// = Eigen::Vector3d(0.00017775,  0.000837968, - 0.000532815);
	GNSSstate = state_;

	return true;
}


const double TIME_ALIGN_ERR = 0.0011;

int isToUpdate(double imutime1, double imutime2, double updatetime)
{

	if (abs(imutime1 - updatetime) < TIME_ALIGN_ERR)
	{
		return 1;
	}
	else if (abs(imutime2 - updatetime) <= TIME_ALIGN_ERR)
	{
		return 2;
	}
	else if (imutime1 < updatetime && updatetime < imutime2)
	{
		return 3;
	}
	else
	{
		return 0;
	}
}


void eskf::ImuProcess()
{

	if (mGNSSUpdate)
	{
		int res = isToUpdate(imuPre.timestamp, imuCur.timestamp, GNSSdata.timestamp);

		//std::cout << imuPre.utimestamp << "     " << GNSSdata.utimestamp << "   " << imuCur.utimestamp << "  " << res << std::endl;
		if (res == 1)
		{
			updatePose(GNSSdata);
			GNSSstate = state_;
			Predict(imuCur);

		}
		else if (res == 2)
		{
			Predict(imuCur);
			updatePose(GNSSdata);
			GNSSstate = state_;
		}
		else if (res == 3)
		{

			filterMessage midimu;
			imuInterpolate(imuPre, imuCur, GNSSdata.timestamp, midimu);
			Predict(midimu);
			updatePose(GNSSdata);
			GNSSstate = state_;
			Predict(imuCur);

		}
		mGNSSUpdate = false;
	}

	else if (mVelUpdate)
	{
		int res = isToUpdate(imuPre.timestamp, imuCur.timestamp, velData.timestamp);
		if (res == 1)
		{
			updateVehicle(velData);
			Predict(imuCur);
		}
		else if (res == 2)
		{
			Predict(imuCur);
			updateVehicle(velData);

		}
		mVelUpdate = false;
	}
	else
	{
		Predict(imuCur);
	}
}





bool eskf::Predict(filterMessage& imuData)
{

	double gyro_noise = 1e-3;
	double acc_noise = 1e-2;
	double dt = imuData.timestamp - state_.time;


	Eigen::Vector3d gyro = 0.5 * (imuData.gyro + lastDataPtr->gyro) - state_.bg;
	Eigen::Vector3d acc = 0.5 * (imuData.acc + lastDataPtr->acc) - state_.ba;

	imuData.dtheta = imuData.dtheta - state_.bg * dt;
	imuData.dvel = imuData.dvel - state_.ba * dt;

	integration(*lastDataPtr.get(), imuData);

	Eigen::Matrix<double, 15, 15> Ft = Eigen::Matrix<double, 15, 15>::Zero();
	Ft = Eigen::Matrix<double, 15, 15>::Zero();
	Ft.block<3, 3>(0, 3) = Eigen::Matrix<double, 3, 3>::Identity();
	Ft.block<3, 3>(3, 6) = state_.q * skewSymmetric(acc); // F23
	Ft.block<3, 3>(3, 12) = state_.q.toRotationMatrix();  // Cbn
	double wsl = sin(lastDataPtr->LLA[0] * DEGREE_TO_RAD);
	double wcl = cos(lastDataPtr->LLA[0] * DEGREE_TO_RAD);
	Eigen::Vector3d wl = Eigen::Vector3d(0, -wcl, -wsl);
	Ft.block<3, 3>(6, 6) = skewSymmetric(wl); // F33
	Ft.block<3, 3>(6, 9) = -state_.q.toRotationMatrix();
	Eigen::Matrix<double, 15, 6> Bt = Eigen::Matrix<double, 15, 6>::Zero();
	Bt.block<3, 3>(3, 3) = state_.q.toRotationMatrix();	 // Cbn
	Bt.block<3, 3>(6, 0) = -state_.q.toRotationMatrix(); //-Cbn
	Ft = Eigen::Matrix<double, 15, 15>::Identity() + Ft * dt;
	Bt = Bt * dt;
	Eigen::Matrix<double, 6, 1> W = Eigen::Matrix<double, 6, 1>::Zero();
	W.head(3) = Eigen::Vector3d(gyro_noise, gyro_noise, gyro_noise);
	W.tail(3) = Eigen::Vector3d(acc_noise, acc_noise, acc_noise);
	error_state.x = Ft * error_state.x + Bt * W;
	Eigen::Matrix<double, 6, 6> Q = Eigen::Matrix<double, 6, 6>::Identity();
	Q.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * gyro_noise * gyro_noise;
	Q.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity() * acc_noise * acc_noise;
	error_state.p = Ft * error_state.p * Ft.transpose() + Bt * Q * Bt.transpose();
	state_.time = imuData.timestamp;
	lastDataPtr->acc = imuData.acc;
	lastDataPtr->gyro = imuData.gyro;

	lastDataPtr->dtheta = imuData.dtheta;
	lastDataPtr->dvel = imuData.dvel;

	return true;
}



void eskf::integration(const filterMessage& imu_pre, const filterMessage& imu_cur)
{
	double dt = 0.01;
	double w = 7.27220521664304e-05;
	Eigen::Vector3d gravity_(0, 0, -9.79484197226504);
	Eigen::Vector3d w_ie_n(0, w * std::cos(lastDataPtr->LLA[0] * DEGREE_TO_RAD), w * std::sin(lastDataPtr->LLA[0] * DEGREE_TO_RAD));

	Eigen::Vector2d rmrn = meridianPrimeVerticalRadius(lastDataPtr->LLA[0] * DEGREE_TO_RAD);
	Eigen::Vector3d wie_n, wen_n;
	wie_n << 0, WGS84_WIE* cos(lastDataPtr->LLA[0] * DEGREE_TO_RAD), WGS84_WIE* sin(lastDataPtr->LLA[0] * DEGREE_TO_RAD);
	wen_n << -state_.v[1] / (rmrn[0] + lastDataPtr->LLA[2]), state_.v[0] / (rmrn[1] + lastDataPtr->LLA[2]),
		state_.v[1] * tan(lastDataPtr->LLA[0] * DEGREE_TO_RAD) / (rmrn[1] + lastDataPtr->LLA[2]);
	gravity_[2] = gravity(lastDataPtr->LLA * DEGREE_TO_RAD);


	// 位置速度
	Eigen::Vector3d dvfb = imu_cur.dvel + 0.5 * imu_cur.dtheta.cross(imu_cur.dvel) +
		1.0 / 12.0 * (imu_pre.dtheta.cross(imu_cur.dvel) + imu_pre.dvel.cross(imu_cur.dtheta));
	// 哥氏项和重力项
	Eigen::Vector3d dv_cor_g = (gravity_ - 2.0 * w_ie_n.cross(state_.v)) * dt;

	// 地球自转补偿项, 省去了enwn项
	Eigen::Vector3d dnn = -(wie_n + wen_n) * dt;
	Eigen::Quaterniond qnn = rotvec2quaternion(dnn);

	Eigen::Vector3d dvel =
		(Eigen::Matrix3d::Identity() - qnn.toRotationMatrix()) * state_.q.toRotationMatrix() * dvfb + dv_cor_g;

	// 前后历元平均速度计算位置
	state_.p += dt * state_.v + 0.5 * dt * dvel;
	state_.v += dvel;


	// 姿态
	Eigen::Vector3d dtheta = imu_cur.dtheta + 1.0 / 12.0 * imu_pre.dtheta.cross(imu_cur.dtheta);
	Eigen::Vector3d  temp1 = (wie_n + wen_n) * dt;
	qnn = rotvec2quaternion(temp1);
	state_.q = qnn * state_.q * rotvec2quaternion(dtheta);
	state_.q.normalize();
}

void eskf::updatePose(const filterMessage& Pose)
{
	Eigen::MatrixXd Y = Eigen::Matrix<double, 6, 1>::Zero();
	Y.block<3, 1>(0, 0) = state_.p - Pose.UTM;
	Eigen::Quaterniond q_bNominal_bMeas;
	q_bNominal_bMeas = state_.q.toRotationMatrix().inverse() * YPR2Quaterniond(Pose.YPR).toRotationMatrix();
	Y.block<3, 1>(3, 0) = Sophus::SO3(q_bNominal_bMeas).log();
	Eigen::MatrixXd Gt = Eigen::Matrix<double, 6, 15>::Zero();
	Gt.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity();
	Gt.block<3, 3>(3, 6) = Eigen::Matrix<double, 3, 3>::Identity();
	Eigen::Matrix<double, 6, 6> Ct = Eigen::Matrix<double, 6, 6>::Identity();
	Eigen::MatrixXd R = Eigen::Matrix<double, 6, 6>::Identity();
	R.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity() * Pose.std[0];
	R.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity() * Pose.std[1];
	Eigen::MatrixXd K = error_state.p * Gt.transpose() * (Gt * error_state.p * Gt.transpose() + Ct * R * Ct.transpose()).inverse();
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(15, 15);
	error_state.p = (I - K * Gt) * error_state.p * (I - K * Gt).transpose() + K * R * K.transpose();
	error_state.x = error_state.x + K * (Y - Gt * error_state.x);

	lastDataPtr->LLA = Pose.LLA;
	lastDataPtr->YPR = Pose.YPR;
	state_.time = Pose.timestamp;

	stateFeedback();
}



void eskf::stateFeedback()
{
	state_.p = state_.p - error_state.x.block<3, 1>(0, 0);
	state_.v = state_.v - error_state.x.block<3, 1>(3, 0);
	Eigen::Vector3d dphi_dir = error_state.x.block<3, 1>(6, 0);

	Eigen::Quaterniond qpn = Eigen::Quaterniond((Sophus::SO3::exp(dphi_dir).matrix()));

	state_.q = state_.q * qpn;
	state_.q.normalize();
	if (IsCovStable(9))
	{
		state_.bg = state_.bg + error_state.x.block<3, 1>(9, 0);
	}
	if (IsCovStable(12))
	{
		state_.ba = state_.ba + error_state.x.block<3, 1>(12, 0);
	}

	error_state.x.setZero();

	Eigen::Matrix<double, 15, 15> J = Eigen::Matrix<double, 15, 15>::Identity();
	J.template block<3, 3>(6, 6) = Eigen::Matrix3d::Identity() - 0.5 * Sophus::SO3::hat(dphi_dir);
}


void eskf::updateVehicle(const filterMessage& vehicle)
{

	Eigen::Vector3d vehicleEnu = vehicle.Vehicle;
	vehicleEnu = (state_.q.toRotationMatrix()) * vehicleEnu;
	Eigen::MatrixXd Y = Eigen::Matrix<double, 3, 1>::Zero();
	Y.block<3, 1>(0, 0) = state_.v - vehicleEnu;
	Eigen::MatrixXd Gt = Eigen::Matrix<double, 3, 15>::Zero();
	Gt.block<3, 3>(0, 3) = Eigen::Matrix<double, 3, 3>::Identity();
	Eigen::Matrix<double, 3, 3> Ct = Eigen::Matrix<double, 3, 3>::Identity();
	Eigen::MatrixXd R = Eigen::Matrix<double, 3, 3>::Identity();
	R.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity() * vehicle.std[0];
	Eigen::MatrixXd K = error_state.p * Gt.transpose() * (Gt * error_state.p * Gt.transpose() + Ct * R * Ct.transpose()).inverse();
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(15, 15);
	error_state.p = (I - K * Gt) * error_state.p * (I - K * Gt).transpose() + K * R * K.transpose();
	error_state.x = error_state.x + K * (Y - Gt * error_state.x);
	stateFeedback();
}




bool eskf::IsCovStable(int INDEX_OFSET)
{
	for (int i = 0; i < 3; ++i)
	{
		if (error_state.p(INDEX_OFSET + i, INDEX_OFSET + i) > 1.0e-12)
		{
			return false;
		}
	}
	return true;
}
