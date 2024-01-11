#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <memory>
#include <iostream>

#define PI 3.141592653589793
struct IMUData
{
    double time;
    Eigen::Vector3d acc;
    Eigen::Vector3d gyro;
};
typedef std::shared_ptr<IMUData> IMUPtr;

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


struct Data
{
    double timeStamp0;
    double timeStamp;
    double timeGPS;
    float cov[4];
    int id;
    Eigen::Vector3d lla;
    Eigen::Vector3d xyz;
    Eigen::Vector3d ypr;
    Eigen::Quaterniond q;
    Eigen::Vector3d v;
    Eigen::Vector3d acc;
    Eigen::Vector3d gyro;
};


class filterMessage
{

public:
    int id = 0;
    double timestamp = 0.0f;
    float steerAngle;
    float cov[4];
    Eigen::Vector3d UTM;
    Eigen::Vector3d LLA;
    Eigen::Vector3d Vehicle;
    Eigen::Vector3d YPR;
    Eigen::Vector3d acc;
    Eigen::Vector3d gyro;
    Eigen::Quaterniond q;
    bool operator()(const std::shared_ptr<filterMessage> a, const std::shared_ptr<filterMessage> b)
    {
        return a->timestamp > b->timestamp;
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


inline Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d& w) {
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
    const Eigen::Vector4d& q) {
    const Eigen::Vector3d& q_vec = q.block(0, 0, 3, 1);
    const double& q4 = q(3);
    Eigen::Matrix3d R =
        (2 * q4 * q4 - 1) * Eigen::Matrix3d::Identity() -
        2 * q4 * skewSymmetric(q_vec) +
        2 * q_vec * q_vec.transpose();
    //TODO: Is it necessary to use the approximation equation
    //    (Equation (87)) when the rotation angle is small?
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


inline Eigen::Matrix<double, 4, 4> Omega(Eigen::Matrix<double, 3, 1> w) {
    Eigen::Matrix<double, 4, 4> mat;
    mat.block(0, 0, 3, 3) = -skewSymmetric(w);
    mat.block(3, 0, 1, 3) = -w.transpose();
    mat.block(0, 3, 3, 1) = w;
    mat(3, 3) = 0;
    return mat;
}

inline Eigen::Matrix<double, 4, 1> quat_multiply(const Eigen::Matrix<double, 4, 1>& q, const Eigen::Matrix<double, 4, 1>& p) {
    Eigen::Matrix<double, 4, 1> q_t;
    Eigen::Matrix<double, 4, 4> Qm;
    // create big L matrix
    Qm.block(0, 0, 3, 3) = q(3, 0) * Eigen::MatrixXd::Identity(3, 3) - skewSymmetric(q.block(0, 0, 3, 1));
    Qm.block(0, 3, 3, 1) = q.block(0, 0, 3, 1);
    Qm.block(3, 0, 1, 3) = -q.block(0, 0, 3, 1).transpose();
    Qm(3, 3) = q(3, 0);
    q_t = Qm * p;
    // ensure unique by forcing q_4 to be >0
    if (q_t(3, 0) < 0) {
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

inline Eigen::Matrix<double, 4, 1> rot_2_quat(const Eigen::Matrix<double, 3, 3>& rot) {
    Eigen::Matrix<double, 4, 1> q;
    double T = rot.trace();
    if ((rot(0, 0) >= T) && (rot(0, 0) >= rot(1, 1)) && (rot(0, 0) >= rot(2, 2))) {
        // cout << "case 1- " << endl;
        q(0) = sqrt((1 + (2 * rot(0, 0)) - T) / 4);
        q(1) = (1 / (4 * q(0))) * (rot(0, 1) + rot(1, 0));
        q(2) = (1 / (4 * q(0))) * (rot(0, 2) + rot(2, 0));
        q(3) = (1 / (4 * q(0))) * (rot(1, 2) - rot(2, 1));

    }
    else if ((rot(1, 1) >= T) && (rot(1, 1) >= rot(0, 0)) && (rot(1, 1) >= rot(2, 2))) {
        // cout << "case 2- " << endl;
        q(1) = sqrt((1 + (2 * rot(1, 1)) - T) / 4);
        q(0) = (1 / (4 * q(1))) * (rot(0, 1) + rot(1, 0));
        q(2) = (1 / (4 * q(1))) * (rot(1, 2) + rot(2, 1));
        q(3) = (1 / (4 * q(1))) * (rot(2, 0) - rot(0, 2));
    }
    else if ((rot(2, 2) >= T) && (rot(2, 2) >= rot(0, 0)) && (rot(2, 2) >= rot(1, 1))) {
        // cout << "case 3- " << endl;
        q(2) = sqrt((1 + (2 * rot(2, 2)) - T) / 4);
        q(0) = (1 / (4 * q(2))) * (rot(0, 2) + rot(2, 0));
        q(1) = (1 / (4 * q(2))) * (rot(1, 2) + rot(2, 1));
        q(3) = (1 / (4 * q(2))) * (rot(0, 1) - rot(1, 0));
    }
    else {
        // cout << "case 4- " << endl;
        q(3) = sqrt((1 + T) / 4);
        q(0) = (1 / (4 * q(3))) * (rot(1, 2) - rot(2, 1));
        q(1) = (1 / (4 * q(3))) * (rot(2, 0) - rot(0, 2));
        q(2) = (1 / (4 * q(3))) * (rot(0, 1) - rot(1, 0));
    }
    if (q(3) < 0) {
        q = -q;
    }
    // normalize and return
    q = q / (q.norm());
    return q;
}

inline Eigen::Matrix<double, 4, 1> quatnorm(Eigen::Matrix<double, 4, 1> q_t) {
    if (q_t(3, 0) < 0) {
        q_t *= -1;
    }
    return q_t / q_t.norm();
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
