#pragma once

#include <vector>
#include <iostream>
#include <memory>
#include <Eigen/Geometry>
#include "data.h"
#include "so3.h"

/* WGS84椭球模型参数
   NOTE:如果使用其他椭球模型需要修改椭球参数 */
const double WGS84_WIE = 7.2921151467E-5;      /* 地球自转角速度*/
const double WGS84_F = 0.0033528106647474805;  /* 扁率 */
const double WGS84_RA = 6378137.0000000000;    /* 长半轴a */
const double WGS84_RB = 6356752.3142451793;    /* 短半轴b */
const double WGS84_GM0 = 398600441800000.00;   /* 地球引力常数 */
const double WGS84_E1 = 0.0066943799901413156; /* 第一偏心率平方 */
const double WGS84_E2 = 0.0067394967422764341; /* 第二偏心率平方 */

#define DEGREE_TO_RAD (0.017453292520)  // 度转弧度

struct State
{
    double time;
    Eigen::Vector3d p;
    Eigen::Vector3d v;
    Eigen::Quaterniond q;
    Eigen::Vector3d bg;
    Eigen::Vector3d ba;
    Eigen::Matrix<double, 15, 15> cov;
};

class filterMessage
{

public:
    int id = 0;
    uint64_t utimestamp = 0;
    double timestamp = 0;
    Eigen::Vector3d LLA;
    Eigen::Vector3d UTM;
    Eigen::Vector3d Vehicle;
    Eigen::Vector3d YPR;
    Eigen::Vector3d acc;
    Eigen::Vector3d gyro;
    Eigen::Vector3d std;
    Eigen::Vector3d dvel;
    Eigen::Vector3d dtheta;
    bool operator()(const filterMessage a, const filterMessage b)
    {
        return a.utimestamp > b.utimestamp;
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
typedef std::shared_ptr<filterMessage> filterMessagePtr;

struct ErrorState
{
    Eigen::Matrix<double, 15, 1> x;
    Eigen::Matrix<double, 15, 15> p;
};

inline Eigen::Vector3d ToEulerAngles(const Eigen::Quaterniond& q)
{
    Eigen::Vector3d angles;
    double sinr_cosp = 2 * (q.w() * q.x() + q.y() * q.z());
    double cosr_cosp = 1 - 2 * (q.x() * q.x() + q.y() * q.y());
    angles[0] = std::atan2(sinr_cosp, cosr_cosp);

    double sinp = 2 * (q.w() * q.y() - q.z() * q.x());
    if (std::abs(sinp) >= 1)
        angles[1] = sinp > 0 ? (EIGEN_PI / 2) : (-EIGEN_PI / 2);
    else
        angles[1] = std::asin(sinp);
    double siny_cosp = 2 * (q.w() * q.z() + q.x() * q.y());
    double cosy_cosp = 1 - 2 * (q.y() * q.y() + q.z() * q.z());
    angles[2] = std::atan2(siny_cosp, cosy_cosp);
    return angles;
}

inline Eigen::Quaterniond YPR2Quaterniond(float roll, float pitch, float yaw)
{
    Eigen::Vector3d eulerAngle(yaw, pitch, roll);
    Eigen::AngleAxisd rollAngle(Eigen::AngleAxisd(eulerAngle(2), Eigen::Vector3d::UnitX()));
    Eigen::AngleAxisd pitchAngle(Eigen::AngleAxisd(eulerAngle(1), Eigen::Vector3d::UnitY()));
    Eigen::AngleAxisd yawAngle(Eigen::AngleAxisd(eulerAngle(0), Eigen::Vector3d::UnitZ()));
    Eigen::Quaterniond quaternion;
    quaternion = yawAngle * pitchAngle * rollAngle;
    return quaternion;
}

inline Eigen::Quaterniond YPR2Quaterniond(Eigen::Vector3d eulerangle)
{
    Eigen::AngleAxisd rollAngle(Eigen::AngleAxisd(eulerangle(2), Eigen::Vector3d::UnitX()));
    Eigen::AngleAxisd pitchAngle(Eigen::AngleAxisd(eulerangle(1), Eigen::Vector3d::UnitY()));
    Eigen::AngleAxisd yawAngle(Eigen::AngleAxisd(eulerangle(0), Eigen::Vector3d::UnitZ()));
    Eigen::Quaterniond quaternion;
    quaternion = yawAngle * pitchAngle * rollAngle;
    return quaternion;
}

inline Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d& w)
{
    Eigen::Matrix3d w_hat;
    w_hat(0, 0) = 0;
    w_hat(0, 1) = -w(2);
    w_hat(0, 2) = w(1);
    w_hat(1, 0) = w(2);
    w_hat(1, 1) = 0;
    w_hat(1, 2) = -w(0);
    w_hat(2, 0) = -w(1);
    w_hat(2, 1) = w(0);
    w_hat(2, 2) = 0;
    return w_hat;
}

inline Eigen::Matrix3d quaternionToRotation(
    const Eigen::Vector4d& q)
{
    const Eigen::Vector3d& q_vec = q.block(0, 0, 3, 1);
    const double& q4 = q(3);
    Eigen::Matrix3d R =
        (2 * q4 * q4 - 1) * Eigen::Matrix3d::Identity() -
        2 * q4 * skewSymmetric(q_vec) + 2 * q_vec * q_vec.transpose();
    // TODO: Is it necessary to use the approximation equation
    //     (Equation (87)) when the rotation angle is small?
    return R;
}

template <typename Derived>
static Eigen::Quaternion<typename Derived::Scalar> deltaQ(const Eigen::MatrixBase<Derived>& theta)
{
    typedef typename Derived::Scalar Scalar_t;

    Eigen::Quaterniond dq;
    Eigen::Matrix<Scalar_t, 3, 1> half_theta = theta;
    half_theta /= static_cast<Scalar_t>(2.0);
    dq.w() = static_cast<Scalar_t>(1.0);
    dq.x() = half_theta.x();
    dq.y() = half_theta.y();
    dq.z() = half_theta.z();
    return dq;
}

inline Eigen::Matrix<double, 4, 4> Omega(Eigen::Matrix<double, 3, 1> w)
{
    Eigen::Matrix<double, 4, 4> mat;
    mat.block(0, 0, 3, 3) = -skewSymmetric(w);
    mat.block(3, 0, 1, 3) = -w.transpose();
    mat.block(0, 3, 3, 1) = w;
    mat(3, 3) = 0;
    return mat;
}

inline Eigen::Matrix<double, 4, 1> quat_multiply(const Eigen::Matrix<double, 4, 1>& q, const Eigen::Matrix<double, 4, 1>& p)
{
    Eigen::Matrix<double, 4, 1> q_t;
    Eigen::Matrix<double, 4, 4> Qm;
    // create big L matrix
    Qm.block(0, 0, 3, 3) = q(3, 0) * Eigen::MatrixXd::Identity(3, 3) - skewSymmetric(q.block(0, 0, 3, 1));
    Qm.block(0, 3, 3, 1) = q.block(0, 0, 3, 1);
    Qm.block(3, 0, 1, 3) = -q.block(0, 0, 3, 1).transpose();
    Qm(3, 3) = q(3, 0);
    q_t = Qm * p;
    // ensure unique by forcing q_4 to be >0
    if (q_t(3, 0) < 0)
    {
        q_t *= -1;
    }
    // normalize and return
    return q_t / q_t.norm();
}

inline Eigen::Matrix<double, 3, 3> quat_2_Rot(const Eigen::Matrix<double, 4, 1>& q)
{
    Eigen::Matrix<double, 3, 3> q_x = skewSymmetric(q.block(0, 0, 3, 1));
    Eigen::MatrixXd Rot = (2 * std::pow(q(3, 0), 2) - 1) * Eigen::MatrixXd::Identity(3, 3) - 2 * q(3, 0) * q_x +
        2 * q.block(0, 0, 3, 1) * (q.block(0, 0, 3, 1).transpose());
    return Rot;
}

inline Eigen::Matrix<double, 4, 1> rot_2_quat(const Eigen::Matrix<double, 3, 3>& rot)
{
    Eigen::Matrix<double, 4, 1> q;
    double T = rot.trace();
    if ((rot(0, 0) >= T) && (rot(0, 0) >= rot(1, 1)) && (rot(0, 0) >= rot(2, 2)))
    {
        // cout << "case 1- " << endl;
        q(0) = sqrt((1 + (2 * rot(0, 0)) - T) / 4);
        q(1) = (1 / (4 * q(0))) * (rot(0, 1) + rot(1, 0));
        q(2) = (1 / (4 * q(0))) * (rot(0, 2) + rot(2, 0));
        q(3) = (1 / (4 * q(0))) * (rot(1, 2) - rot(2, 1));
    }
    else if ((rot(1, 1) >= T) && (rot(1, 1) >= rot(0, 0)) && (rot(1, 1) >= rot(2, 2)))
    {
        // cout << "case 2- " << endl;
        q(1) = sqrt((1 + (2 * rot(1, 1)) - T) / 4);
        q(0) = (1 / (4 * q(1))) * (rot(0, 1) + rot(1, 0));
        q(2) = (1 / (4 * q(1))) * (rot(1, 2) + rot(2, 1));
        q(3) = (1 / (4 * q(1))) * (rot(2, 0) - rot(0, 2));
    }
    else if ((rot(2, 2) >= T) && (rot(2, 2) >= rot(0, 0)) && (rot(2, 2) >= rot(1, 1)))
    {
        // cout << "case 3- " << endl;
        q(2) = sqrt((1 + (2 * rot(2, 2)) - T) / 4);
        q(0) = (1 / (4 * q(2))) * (rot(0, 2) + rot(2, 0));
        q(1) = (1 / (4 * q(2))) * (rot(1, 2) + rot(2, 1));
        q(3) = (1 / (4 * q(2))) * (rot(0, 1) - rot(1, 0));
    }
    else
    {
        // cout << "case 4- " << endl;
        q(3) = sqrt((1 + T) / 4);
        q(0) = (1 / (4 * q(3))) * (rot(1, 2) - rot(2, 1));
        q(1) = (1 / (4 * q(3))) * (rot(2, 0) - rot(0, 2));
        q(2) = (1 / (4 * q(3))) * (rot(0, 1) - rot(1, 0));
    }
    if (q(3) < 0)
    {
        q = -q;
    }
    // normalize and return
    q = q / (q.norm());
    return q;
}

inline Eigen::Matrix<double, 4, 1> quatnorm(Eigen::Matrix<double, 4, 1> q_t)
{
    if (q_t(3, 0) < 0)
    {
        q_t *= -1;
    }
    return q_t / q_t.norm();
}

inline Eigen::Quaterniond rotvec2quaternion(const Eigen::Vector3d& rotvec)
{
    double angle = rotvec.norm();
    Eigen::Vector3d vec = rotvec.normalized();
    return Eigen::Quaterniond(Eigen::AngleAxisd(angle, vec));
}

inline Eigen::Vector2d meridianPrimeVerticalRadius(double lat)
{
    double tmp, sqrttmp;

    tmp = sin(lat);
    tmp *= tmp;
    tmp = 1 - WGS84_E1 * tmp;
    sqrttmp = sqrt(tmp);

    return { WGS84_RA * (1 - WGS84_E1) / (sqrttmp * tmp), WGS84_RA / sqrttmp };
}

inline double gravity(const Eigen::Vector3d& blh)
{

    double sin2 = sin(blh[0]);
    sin2 *= sin2;

    return 9.7803267715 * (1 + 0.0052790414 * sin2 + 0.0000232718 * sin2 * sin2) +
        blh[2] * (0.0000000043977311 * sin2 - 0.0000030876910891) + 0.0000000000007211 * blh[2] * blh[2];
}

enum StateMembers
{
    StateMemberX = 0,
    StateMemberY,
    StateMemberZ,
    StateMemberYaw,
    StateMemberVx,
};

const double TAU = 6.283185307179587;

class eskf
{

public:
    bool Predict(filterMessage& imuData);
    void updatePose(const filterMessage& pose);
    void updateVehicle(const filterMessage& vehicle);
    bool InitPose(const filterMessage& data);
    void getPose(State& state);
    eskf();

    void addGNSSPose(const filterMessage& pose)
    {
        GNSSdata = pose;
        mGNSSUpdate = true;
    }

    void addImuData(const filterMessage& imu)
    {
        //lastDataPtr->acc = imu.acc;
        //lastDataPtr->gyro = imu.gyro;
        imuPre = imuCur;
        imuCur = imu;

        if (imuCur.timestamp - imuPre.timestamp > 10)
        {
            imuPre = imuCur;
        }

        imuCur.dtheta = imuCur.gyro * (imuCur.timestamp - imuPre.timestamp);
        imuCur.dvel = imuCur.acc * (imuCur.timestamp - imuPre.timestamp);
    }

    void AddVelData(const filterMessage& Data)
    {
        velData = Data;
        mVelUpdate = true;
    }

    void ImuProcess();

    void imuInterpolate(const filterMessage& imu1, filterMessage& imu2, const double timestamp, filterMessage& midimu)
    {

        if (imu1.timestamp > timestamp || imu2.timestamp < timestamp)
        {
            return;
        }

        double lamda = (timestamp - imu1.timestamp) / (imu2.timestamp - imu1.timestamp);

        midimu.timestamp = timestamp;
        midimu.gyro = (imu2.gyro - imu1.gyro) * lamda + imu1.gyro;
        midimu.acc = (imu2.acc - imu1.acc) * lamda + imu1.acc;
    }

private:
    State state_;
    State GNSSstate;

    ErrorState error_state;
    void integration(const filterMessage& imu_pre, const filterMessage& imu_cur);
    void stateFeedback();

    filterMessagePtr lastDataPtr;
    bool has_imu;
    bool IsCovStable(int INDEX_OFSET);
    Eigen::Vector3d LastP;

    filterMessage GNSSdata;
    bool mGNSSUpdate;

    filterMessage velData;
    bool mVelUpdate;

    filterMessage imuPre;
    filterMessage imuCur;
};