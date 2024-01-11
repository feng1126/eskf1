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

void eskf::getPose(State &state)
{

	state = state_;
}

bool eskf::InitPose(const filterMessagePtr &Pose)
{
	state_.time = Pose->timestamp;
	state_.p = Pose->UTM;
	state_.q = YPR2Quaterniond(Pose->YPR);
	state_.v = Pose->Vehicle;
	state_.ba = Eigen::Vector3d(-0.282131, -0.0121983, 0.01123);
	state_.bg = Eigen::Vector3d(-0.0066, 0.00044, 0.006);
	return true;
}

bool eskf::Predict(const filterMessagePtr &imuData)
{
	if (!has_imu)
	{
		lastDataPtr->acc = imuData->acc;
		lastDataPtr->gyro = imuData->gyro;
		has_imu = true;
		return false;
	}




	double gyro_noise = 1.2936591392287652e-03;
	double acc_noise = 1.9410086926062301e-02;
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
	return true;
}

void eskf::updateXYZ(const filterMessagePtr &XYZ, const double cov)
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

void eskf::updatePose(const filterMessagePtr &Pose, const Eigen::VectorXd cov)
{

	double dp_noise = cov[0];
	double geo_x, geo_y, geo_z;
	Eigen::MatrixXd Y = Eigen::Matrix<double, 6, 1>::Zero();
	Y.block<3, 1>(0, 0) = state_.p - Pose->UTM;
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
	Eigen::Vector3d dphi_dir = error_state.x.block<3, 1>(6, 0);
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

void eskf::updateVehicle(const filterMessagePtr &vehicle, const double cov)
{
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

Eigen::Matrix3d eskf::ComputeEarthTranform(const double &dt, const Eigen::Vector3d &velocity)
{
	double w = 7.27220521664304e-05;
	Eigen::Vector3d gravity_(0, 0, -9.79484197226504);
	Eigen::Vector3d w_ie_n(0, w * std::cos(lastDataPtr->LLA[0] * EIGEN_PI / 180), w * std::sin(lastDataPtr->LLA[0] * EIGEN_PI / 180));
	double rm = 6353346.18315;
	double rn = 6384140.52699;
	Eigen::Vector3d w_en_n(
		-velocity[1] / (rm + lastDataPtr->LLA[2]),
		velocity[0] / (rn + lastDataPtr->LLA[2]),
		velocity[0] / (rn + lastDataPtr->LLA[2]) * std::tan(lastDataPtr->LLA[0] * EIGEN_PI / 180));
	Eigen::Vector3d wtemp = (w_ie_n + w_en_n) * dt;
	Eigen::AngleAxisd angle_axisd(wtemp.norm(), wtemp.normalized());
	return angle_axisd.toRotationMatrix().transpose();
}

void eskf::predictNewState(const double &dt, const Eigen::Vector3d &gyro, const Eigen::Vector3d &acc)
{
	Eigen::Vector3d gn(0, 0, -9.79484197226504);
	Eigen::Vector3d pp = state_.p;
	Eigen::Vector3d vv = state_.v;
	Eigen::Quaterniond qq = state_.q;
	Eigen::Matrix3d R_nm_nm = ComputeEarthTranform(dt, vv);
	Eigen::Vector3d w_hat = 0.5 * (gyro + lastDataPtr->gyro) - state_.bg;
	Eigen::Vector3d a_hat = 0.5 * (acc + lastDataPtr->acc) - state_.ba;
	if (w_hat.norm() > 1e-5)
	{
		state_.q = qq * R_nm_nm * Sophus::SO3::exp(w_hat * dt).matrix();
	}
	else
	{
		state_.q = qq * R_nm_nm;
	}
	state_.v = vv + dt * (state_.q * a_hat + gn);
	state_.p = pp + 0.5 * dt * (vv + state_.v);
}

void eskf::predictNewRkState(const double &dt, const Eigen::Vector3d &gyro, const Eigen::Vector3d &acc)
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
