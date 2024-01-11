#include "eskf.h"

constexpr double kDegreeToRadian = PI / 180.;

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

bool eskf::InitPose(const filterMessagePtr& Pose)
{
	state_.time = Pose->timestamp;
	state_.p = Pose->UTM;
	state_.q = YPR2Quaterniond(Pose->YPR);
	state_.v = state_.q * Pose->Vehicle;
	state_.ba; // = Eigen::Vector3d(-0.282131, -0.0121983, 0.01123);
	state_.bg; // = Eigen::Vector3d(-0.0066, 0.00044, 0.006);

	printf("InitPose:%0.7f \n", Pose->YPR[0]);

	double heading = ToEulerAngles(state_.q)[2];
	lastDataPtr = Pose;

	if (heading > 2 * EIGEN_PI)
	{
		heading = heading - 2.0 * EIGEN_PI;
	}

	if (heading < 0)
		heading = heading + 2 * EIGEN_PI;

	heading = heading * 180.0 / EIGEN_PI;
	printf("InitPose heading :%0.7f \n", heading);
	return true;
}

bool eskf::Predict(const filterMessagePtr& imuData)
{

	// std::cout << state_.p << std::endl << std::endl;

	if (!has_imu)
	{
		lastDataPtr->acc = imuData->acc;
		lastDataPtr->gyro = imuData->gyro;
		has_imu = true;
		return false;
	}
	double gyro_noise = 1e-04;
	double acc_noise = 1e-03;
	double dt = imuData->timestamp - state_.time;
	Eigen::Vector3d gyro = 0.5 * (imuData->gyro + lastDataPtr->gyro) - state_.bg;
	Eigen::Vector3d acc = 0.5 * (imuData->acc + lastDataPtr->acc) - state_.ba;
	predictNewState(dt, imuData->gyro, imuData->acc);
	Eigen::Matrix<double, 15, 15> Ft = Eigen::Matrix<double, 15, 15>::Zero();
	Ft = Eigen::Matrix<double, 15, 15>::Zero();
	Ft.block<3, 3>(0, 3) = Eigen::Matrix<double, 3, 3>::Identity();
	Ft.block<3, 3>(3, 6) = state_.q * skewSymmetric(acc); // F23
	Ft.block<3, 3>(3, 12) = state_.q.toRotationMatrix();  // Cbn
	double wsl = sin(lastDataPtr->LLA[0] * PI / 180);
	double wcl = cos(lastDataPtr->LLA[0] * PI / 180);
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
	state_.time = imuData->timestamp;
	lastDataPtr->acc = imuData->acc;
	lastDataPtr->gyro = imuData->gyro;

	// std::cout << state_.p << std::endl << std::endl;

	return true;
}





void eskf::updateXYZ(const filterMessagePtr& XYZ, const double cov)
{
	double dp_noise = cov;
	double geo_x, geo_y, geo_z;
	Eigen::MatrixXd Y = Eigen::Matrix<double, 3, 1>::Zero();
	Y.block<3, 1>(0, 0) = state_.p - XYZ->UTM;
	Eigen::MatrixXd Gt = Eigen::Matrix<double, 3, 15>::Zero();
	Gt.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity();
	Eigen::Matrix<double, 3, 3> Ct = Eigen::Matrix<double, 3, 3>::Identity();
	Eigen::MatrixXd R = Eigen::Matrix<double, 3, 3>::Identity();
	R.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity() * dp_noise;
	Eigen::MatrixXd K = error_state.p * Gt.transpose() * (Gt * error_state.p * Gt.transpose() + Ct * R * Ct.transpose()).inverse();
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(15, 15);
	error_state.p = (I - K * Gt) * error_state.p * (I - K * Gt).transpose() + K * R * K.transpose();
	error_state.x = error_state.x + K * (Y - Gt * error_state.x);
	state_.p = state_.p - error_state.x.block<3, 1>(0, 0);
	state_.v = state_.v - error_state.x.block<3, 1>(3, 0);
	Eigen::Vector3d dphi_dir = error_state.x.block<3, 1>(6, 0);
	state_.q = state_.q.toRotationMatrix() * Sophus::SO3::exp(dphi_dir).matrix();
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
	state_.time = XYZ->timestamp;
	return;
}

void eskf::updatePose(const filterMessagePtr& Pose, const Eigen::VectorXd cov)
{

	double dp_noise = cov[0];
	double geo_x, geo_y, geo_z;
	Eigen::MatrixXd Y = Eigen::Matrix<double, 6, 1>::Zero();
	Y.block<3, 1>(0, 0) = state_.p - Pose->UTM;

	// std::cout << state_.p << std::endl << std::endl;

	// std::cout << Pose->UTM << std::endl << std::endl;

	// std::cout << Y << std::endl << std::endl;
	Eigen::Quaterniond q_bNominal_bMeas;
	q_bNominal_bMeas = YPR2Quaterniond(Pose->YPR).toRotationMatrix() * state_.q.toRotationMatrix().inverse();
	Y.block<3, 1>(3, 0) = Sophus::SO3(q_bNominal_bMeas).log();
	Eigen::MatrixXd Gt = Eigen::Matrix<double, 6, 15>::Zero();
	Gt.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity();
	Gt.block<3, 3>(3, 6) = Eigen::Matrix<double, 3, 3>::Identity();
	Eigen::Matrix<double, 6, 6> Ct = Eigen::Matrix<double, 6, 6>::Identity();
	Eigen::MatrixXd R = Eigen::Matrix<double, 6, 6>::Identity();
	R.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity() * cov[0];
	R.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity() * cov[1];
	Eigen::MatrixXd K = error_state.p * Gt.transpose() * (Gt * error_state.p * Gt.transpose() + Ct * R * Ct.transpose()).inverse();
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(15, 15);
	error_state.p = (I - K * Gt) * error_state.p * (I - K * Gt).transpose() + K * R * K.transpose();
	error_state.x = error_state.x + K * (Y - Gt * error_state.x);
	state_.p = state_.p - error_state.x.block<3, 1>(0, 0);
	state_.v = state_.v - error_state.x.block<3, 1>(3, 0);
	Eigen::Vector3d dphi_dir = -error_state.x.block<3, 1>(6, 0);
	state_.q = Sophus::SO3::exp(dphi_dir).matrix() * state_.q.toRotationMatrix();
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
	lastDataPtr->LLA = Pose->LLA;
	lastDataPtr->YPR = Pose->YPR;
	state_.time = Pose->timestamp;
}

void eskf::updateVehicle(const filterMessagePtr& vehicle, const double cov)
{

	// std::cout << "state_.v 0          " << state_.v.transpose() << std::endl;

	Eigen::Vector3d vehicleEnu = vehicle->Vehicle;
	vehicleEnu = (state_.q.toRotationMatrix()) * vehicleEnu;
	double dp_noise = cov;
	double geo_x, geo_y, geo_z;
	Eigen::MatrixXd Y = Eigen::Matrix<double, 3, 1>::Zero();
	Y.block<3, 1>(0, 0) = state_.v - vehicleEnu;
	Eigen::MatrixXd Gt = Eigen::Matrix<double, 3, 15>::Zero();
	Gt.block<3, 3>(0, 3) = Eigen::Matrix<double, 3, 3>::Identity();
	Eigen::Matrix<double, 3, 3> Ct = Eigen::Matrix<double, 3, 3>::Identity();
	Eigen::MatrixXd R = Eigen::Matrix<double, 3, 3>::Identity();
	R.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity() * dp_noise;
	Eigen::MatrixXd K = error_state.p * Gt.transpose() * (Gt * error_state.p * Gt.transpose() + Ct * R * Ct.transpose()).inverse();
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(15, 15);
	error_state.p = (I - K * Gt) * error_state.p * (I - K * Gt).transpose() + K * R * K.transpose();
	error_state.x = error_state.x + K * (Y - Gt * error_state.x);
	state_.p = state_.p - error_state.x.block<3, 1>(0, 0);
	state_.v = state_.v - error_state.x.block<3, 1>(3, 0);
	Eigen::Vector3d dphi_dir = error_state.x.block<3, 1>(6, 0);

	// std::cout << "state_.v 1          " << state_.v.transpose() << std::endl;

	// std::cout << "dphi_dir " <<  dphi_dir.transpose() << std::endl;
	state_.q = state_.q.toRotationMatrix() * Sophus::SO3::exp(dphi_dir).matrix();
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
	state_.time = vehicle->timestamp;
	return;
}





Eigen::Matrix3d eskf::ComputeEarthTranform(const double& dt, const Eigen::Vector3d& velocity)
{

	Eigen::Vector2d rmrn = Earth::meridianPrimeVerticalRadius(lastDataPtr->LLA[0] * EIGEN_PI / 180);
	

	Eigen::Vector3d wie_n, wen_n;
	wie_n << WGS84_WIE * cos(lastDataPtr->LLA[0] * EIGEN_PI / 180), 0, -WGS84_WIE * sin(lastDataPtr->LLA[0] * EIGEN_PI / 180);
	wen_n << state_.v[1] / (rmrn[1] + lastDataPtr->LLA[2]), -state_.v[0] / (rmrn[0] + lastDataPtr->LLA[2]),
		-state_.v[1] * tan(lastDataPtr->LLA[0] * EIGEN_PI / 180) / (rmrn[1] + lastDataPtr->LLA[2]);
	Eigen::Vector3d wtemp = 0.5* (wie_n + wie_n) * dt ;
	Eigen::AngleAxisd angle_axisd(wtemp.norm(), wtemp.normalized());
	return angle_axisd.toRotationMatrix().transpose();
}


void eskf::velUpdate(const double& dt, const Eigen::Vector3d& gyro, const Eigen::Vector3d& acc)
{
	Eigen::Vector2d rmrn = Earth::meridianPrimeVerticalRadius(lastDataPtr->LLA[0] * EIGEN_PI / 180);
	Eigen::Vector3d wie_n, wen_n;
	wie_n << 0, -WGS84_WIE * cos(lastDataPtr->LLA[0] * EIGEN_PI / 180),  WGS84_WIE * sin(lastDataPtr->LLA[0] * EIGEN_PI / 180);
	wen_n << -state_.v[1] / (rmrn[1] + lastDataPtr->LLA[2]), -state_.v[0] / (rmrn[0] + lastDataPtr->LLA[2]),
		state_.v[1] * tan(lastDataPtr->LLA[0] * EIGEN_PI / 180) / (rmrn[0] + lastDataPtr->LLA[2]);
	double gravity = Earth::gravity(lastDataPtr->LLA * EIGEN_PI / 180);

	Eigen::Vector3d temp1, temp2, temp3;

	temp1 = (gyro * dt * 0.2).cross((acc * dt * 0.2)) / 2;
	temp2 = (gyro * dt * 0.2).cross((acc * dt * 0.2)) / 12;
	temp3 = (acc * dt * 0.2).cross((gyro * dt * 0.2)) / 12;

	Eigen::Vector3d d_vfb = acc * dt + temp1 + temp2 + temp3;

	temp1 = (wie_n + wen_n) * dt / 2;

	Quaterniond qnn = rotvec2quaternion(temp1);


	Eigen::Vector3d d_vgn = Eigen::Vector3d(0, 0, gravity)* dt * 0.2;


	state_.v = state_.v + (Matrix3d::Identity() + qnn.toRotationMatrix()) * state_.q * d_vfb + d_vgn;

	Vector3d dtheta  = (gyro * dt * 0.2) + (gyro * dt * 0.2).cross((gyro * dt * 0.2)) / 12;

	state_.q = qnn * state_.q * rotvec2quaternion(dtheta);

	state_.p = state_.p + 0.5 * dt * state_.v;
}





void eskf::predictNewState(const double& dt, const Eigen::Vector3d& gyro, const Eigen::Vector3d& acc)
{

	velUpdate(dt, gyro, acc);
	//attUpdate(dt, gyro, acc);
	
	//state_.v = vv + dt * (state_.q * a_hat + gn);

}

void eskf::predictNewRkState(const double& dt, const Eigen::Vector3d& gyro, const Eigen::Vector3d& acc)
{

	Eigen::Vector3d w_hat = 0.5 * (gyro + gyro) - state_.bg;
	Eigen::Vector3d a_hat = 0.5 * (acc + acc) - state_.ba;

	Eigen::Vector3d gravity_(0, 0, -9.79484197226504);
	double gyro_norm = gyro.norm();
	Eigen::Matrix4d Omega = Eigen::Matrix4d::Zero();
	Omega.block<3, 3>(0, 0) = -skewSymmetric(w_hat);
	Omega.block<3, 1>(0, 3) = w_hat;
	Omega.block<1, 3>(3, 0) = -w_hat;

	Eigen::Quaterniond q = state_.q;
	Eigen::Vector3d v = state_.v;
	Eigen::Vector3d p = state_.p;

	Eigen::Matrix3d R_nm_nm = ComputeEarthTranform(dt, v);

	Eigen::Matrix3d dRt, dRt2;

	if (w_hat.norm() > 1e-5)
	{
		dRt = q * R_nm_nm * deltaQ(w_hat * dt);
		dRt2 = q * R_nm_nm * deltaQ(0.5 * w_hat * dt);
	}
	else
	{
		dRt = q * R_nm_nm;
		dRt2 = q * R_nm_nm;
	}

	Eigen::Vector3d omega_predicted1, omega_predicted2, acc_predicted1, acc_predicted2;
	{
		omega_predicted1 = (gyro - lastDataPtr->gyro) * dt * 0.5 + lastDataPtr->gyro - state_.bg;
		omega_predicted2 = (gyro - lastDataPtr->gyro) * dt + lastDataPtr->gyro - state_.bg;
		acc_predicted1 = (acc - lastDataPtr->acc) * dt * 0.5 + lastDataPtr->acc - state_.ba;
		acc_predicted2 = (acc - lastDataPtr->acc) * dt + lastDataPtr->acc - state_.ba;
	}

	// k1 = f(tn, yn)
	Eigen::Vector3d kv1 = q * acc_predicted2 + gravity_;
	Eigen::Vector3d kp1 = v;

	// k2 = f(tn+dt/2, yn+k1*dt/2)
	Eigen::Vector3d kv2 = dRt2 * acc_predicted1 + gravity_;
	Eigen::Vector3d kp2 = v + kv1 * dt / 2;

	// k3 = f(tn+dt/2, yn+k2*dt/2)
	Eigen::Vector3d kv3 = dRt2 * acc_predicted1 + gravity_;
	Eigen::Vector3d kp3 = v + kv2 * dt / 2;

	// k4 = f(tn+dt, yn+k3*dt)
	Eigen::Vector3d kv4 = dRt * acc_predicted2 + gravity_;
	Eigen::Vector3d kp4 = v + kv3 * dt;

	v = v + dt / 6 * (kv1 + 2 * kv2 + 2 * kv3 + kv4);
	p = p + dt / 6 * (kp1 + 2 * kp2 + 2 * kp3 + kp4);

	state_.v = v;
	state_.p = p;
	state_.q = dRt;
}

bool eskf::IsCovStable(int INDEX_OFSET)
{
	for (int i = 0; i < 3; ++i)
	{
		if (error_state.p(INDEX_OFSET + i, INDEX_OFSET + i) > 1.0e-11)
		{
			return false;
		}
	}
	return true;
}
