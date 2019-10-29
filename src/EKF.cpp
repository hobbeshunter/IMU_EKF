#include "EKF.h"

#include <Eigen/LU>
#include <Wire.h>

#define SGN(X) ((X > 0) - (X < 0))
#define RAD_TO_DEG (180.0 / M_PI)
#define MS2_TO_G (1.0 / 9.81)

EKF::EKF(precision d_t)
{
    d_t_ = d_t;

    x_.setZero();
    P_.setIdentity();
    Q_.setIdentity();
    Q_(0, 0) = 0.05;
    Q_(1, 1) = 0.05;
    Q_(2, 2) = 0.05;
    Q_(3, 3) = 0.05;
    Q_(4, 4) = 0.05;
    Q_(5, 5) = 0.05;
    Q_(6, 6) = 0.025;
    Q_(7, 7) = 0.025;
    Q_(8, 8) = 0.025;
    Q_(9, 9) = 0.025;
    Q_(10, 10) = 0.025;
    Q_(11, 11) = 0.025;
    Q_(12, 12) = 0.01;
    Q_(13, 13) = 0.01;
    Q_(14, 14) = 0.01;

    R_IMU_.setIdentity();
    R_IMU_(0, 0) = 0.0000045494;
    R_IMU_(1, 1) = 0.0000039704;
    R_IMU_(2, 2) = 0.0000093844;
    R_IMU_(3, 3) = 0.0141615383;
    R_IMU_(4, 4) = 0.0164647549;
    R_IMU_(5, 5) = 0.0100094303;
    // R_IMU_(0, 0) = 0.0141615383;
    // R_IMU_(1, 1) = 0.0164647549;
    // R_IMU_(2, 2) = 0.0100094303;
}

void EKF::init_with_acc(float ax, float ay, float az) {
    // see https://cache.freescale.com/files/sensors/doc/app_note/AN3461.pdf p. 16
    x_.setZero();
    x_(3) = atan(-ay / sqrt(ax * ax + az * az));
    x_(4) = atan2(ax, SGN(-az) * sqrt(az * az + 0.01 * ay * ay));
}

void EKF::predict()
{
    precision roll = x_(3);
    precision pitch = x_(4);
    precision yaw = x_(5);
    precision xVel = x_(6);
    precision yVel = x_(7);
    precision zVel = x_(8);
    precision pitchVel = x_(10);
    precision yawVel = x_(11);
    precision xAcc = x_(12);
    precision yAcc = x_(13);
    precision zAcc = x_(14);

    precision sp = sin(pitch);
    precision cp = cos(pitch);
    precision cpi = (fabs(cp) >= 1e-8) ? 1.0 / cp : 0;
    precision tp = sp * cpi;

    precision sr = sin(roll);
    precision cr = cos(roll);

    precision sy = sin(yaw);
    precision cy = cos(yaw);

    Eigen::Matrix<precision, STATE_SIZE, STATE_SIZE> f;
    f.setIdentity();
    f(0, 6) = cy * cp * d_t_;
    f(0, 7) = (cy * sp * sr - sy * cr) * d_t_;
    f(0, 8) = (cy * sp * cr + sy * sr) * d_t_;
    f(0, 12) = 0.5 * f(0, 6) * d_t_;
    f(0, 13) = 0.5 * f(0, 7) * d_t_;
    f(0, 14) = 0.5 * f(0, 8) * d_t_;
    f(1, 6) = sy * cp * d_t_;
    f(1, 7) = (sy * sp * sr + cy * cr) * d_t_;
    f(1, 8) = (sy * sp * cr - cy * sr) * d_t_;
    f(1, 12) = 0.5 * f(1, 6) * d_t_;
    f(1, 13) = 0.5 * f(1, 7) * d_t_;
    f(1, 14) = 0.5 * f(1, 8) * d_t_;
    f(2, 6) = -sp * d_t_;
    f(2, 7) = cp * sr * d_t_;
    f(2, 8) = cp * cr * d_t_;
    f(2, 12) = 0.5 * f(2, 6) * d_t_;
    f(2, 13) = 0.5 * f(2, 7) * d_t_;
    f(2, 14) = 0.5 * f(2, 8) * d_t_;
    f(3, 9) = d_t_;
    f(3, 10) = sr * tp * d_t_;
    f(3, 11) = cr * tp * d_t_;
    f(4, 10) = cr * d_t_;
    f(4, 11) = -sr * d_t_;
    f(5, 10) = sr * cpi * d_t_;
    f(5, 11) = cr * cpi * d_t_;
    f(6, 12) = d_t_;
    f(7, 13) = d_t_;
    f(8, 14) = d_t_;

    precision xCoeff = 0.0;
    precision yCoeff = 0.0;
    precision zCoeff = 0.0;
    precision oneHalfATSquared = 0.5 * d_t_ * d_t_;

    yCoeff = cy * sp * cr + sy * sr;
    zCoeff = -cy * sp * sr + sy * cr;
    precision dFx_dR = (yCoeff * yVel + zCoeff * zVel) * d_t_ +
                       (yCoeff * yAcc + zCoeff * zAcc) * oneHalfATSquared;
    precision dFR_dR = 1.0 + (cr * tp * pitchVel - sr * tp * yawVel) * d_t_;

    xCoeff = -cy * sp;
    yCoeff = cy * cp * sr;
    zCoeff = cy * cp * cr;
    precision dFx_dP = (xCoeff * xVel + yCoeff * yVel + zCoeff * zVel) * d_t_ +
                       (xCoeff * xAcc + yCoeff * yAcc + zCoeff * zAcc) * oneHalfATSquared;
    precision dFR_dP = (cpi * cpi * sr * pitchVel + cpi * cpi * cr * yawVel) * d_t_;

    xCoeff = -sy * cp;
    yCoeff = -sy * sp * sr - cy * cr;
    zCoeff = -sy * sp * cr + cy * sr;
    precision dFx_dY = (xCoeff * xVel + yCoeff * yVel + zCoeff * zVel) * d_t_ +
                       (xCoeff * xAcc + yCoeff * yAcc + zCoeff * zAcc) * oneHalfATSquared;

    yCoeff = sy * sp * cr - cy * sr;
    zCoeff = -sy * sp * sr - cy * cr;
    precision dFy_dR = (yCoeff * yVel + zCoeff * zVel) * d_t_ +
                       (yCoeff * yAcc + zCoeff * zAcc) * oneHalfATSquared;
    precision dFP_dR = (-sr * pitchVel - cr * yawVel) * d_t_;

    xCoeff = -sy * sp;
    yCoeff = sy * cp * sr;
    zCoeff = sy * cp * cr;
    precision dFy_dP = (xCoeff * xVel + yCoeff * yVel + zCoeff * zVel) * d_t_ +
                       (xCoeff * xAcc + yCoeff * yAcc + zCoeff * zAcc) * oneHalfATSquared;

    xCoeff = cy * cp;
    yCoeff = cy * sp * sr - sy * cr;
    zCoeff = cy * sp * cr + sy * sr;
    precision dFy_dY = (xCoeff * xVel + yCoeff * yVel + zCoeff * zVel) * d_t_ +
                       (xCoeff * xAcc + yCoeff * yAcc + zCoeff * zAcc) * oneHalfATSquared;

    yCoeff = cp * cr;
    zCoeff = -cp * sr;
    precision dFz_dR = (yCoeff * yVel + zCoeff * zVel) * d_t_ +
                       (yCoeff * yAcc + zCoeff * zAcc) * oneHalfATSquared;
    precision dFY_dR = (cr * cpi * pitchVel - sr * cpi * yawVel) * d_t_;

    xCoeff = -cp;
    yCoeff = -sp * sr;
    zCoeff = -sp * cr;
    precision dFz_dP = (xCoeff * xVel + yCoeff * yVel + zCoeff * zVel) * d_t_ +
                       (xCoeff * xAcc + yCoeff * yAcc + zCoeff * zAcc) * oneHalfATSquared;
    precision dFY_dP = (sr * tp * cpi * pitchVel - cr * tp * cpi * yawVel) * d_t_;

    Eigen::Matrix<precision, STATE_SIZE, STATE_SIZE> A(f);
    // A = f;
    A(0, 3) = dFx_dR;
    A(0, 4) = dFx_dP;
    A(0, 5) = dFx_dY;
    A(1, 3) = dFy_dR;
    A(1, 4) = dFy_dP;
    A(1, 5) = dFy_dY;
    A(2, 3) = dFz_dR;
    A(2, 4) = dFz_dP;
    A(3, 3) = dFR_dR;
    A(3, 4) = dFR_dP;
    A(4, 3) = dFP_dR;
    A(5, 3) = dFY_dR;
    A(5, 4) = dFY_dP;

    // x = f(x)
    x_ = f * x_;

    // Handle wrapping
    wrap_angle(x_(3));
    wrap_angle(x_(4));
    wrap_angle(x_(5));

    // P = A * P * A' + Q
    P_ = A * P_ * A.transpose();
    P_ += d_t_ * Q_;
}

void EKF::correct_IMU(float gx, float gy, float gz, float ax, float ay, float az)
{
    // z
    Eigen::Matrix<precision, MEASSUREMENT_IMU_SIZE, 1> z;
    z(0) = gx;
    z(1) = gy;
    z(2) = gz;
    z(3) = ax;
    z(4) = ay;
    z(5) = az;
    // z(0) = ax;
    // z(1) = ay;
    // z(2) = az;

    // h(x, 0)
    precision sr = sin(x_(3));
    precision cr = cos(x_(3));
    precision sp = sin(x_(4));
    precision cp = cos(x_(4));

    Eigen::Matrix<precision, MEASSUREMENT_IMU_SIZE, 1> h;
    h(0) = x_(9) * RAD_TO_DEG;
    h(1) = x_(10) * RAD_TO_DEG;
    h(2) = x_(11) * RAD_TO_DEG;
    h(3) = x_(12) * MS2_TO_G + sp;
    h(4) = x_(13) * MS2_TO_G - cp * sr;
    h(5) = x_(14) * MS2_TO_G - cp * cr;
    // h(0) = x_(12) * MS2_TO_G + sp;
    // h(1) = x_(13) * MS2_TO_G - cp * sr;
    // h(2) = x_(14) * MS2_TO_G - cp * cr;

    // H
    Eigen::Matrix<precision, MEASSUREMENT_IMU_SIZE, STATE_SIZE> H;
    H.setZero();
    H(0, 9) = RAD_TO_DEG;
    H(1, 10) = RAD_TO_DEG;
    H(2, 11) = RAD_TO_DEG;
    H(3, 4) = cp;
    H(3, 12) = MS2_TO_G;
    H(4, 3) = -cp * cr;
    H(4, 4) = sp * sr;
    H(4, 13) = MS2_TO_G;
    H(5, 3) = cp * sr;
    H(5, 4) = sp * cr;
    H(5, 14) = MS2_TO_G;
    // H(0, 4) = cp;
    // H(0, 12) = MS2_TO_G;
    // H(1, 3) = -cp * cr;
    // H(1, 4) = sp * sr;
    // H(1, 13) = MS2_TO_G;
    // H(2, 3) = cp * sr;
    // H(2, 4) = sp * cr;
    // H(2, 14) = MS2_TO_G;

    // Kalman Gain
    // K = P * H' (H * P * H' + V * R * V')^-
    Eigen::Matrix<precision, STATE_SIZE, MEASSUREMENT_IMU_SIZE> PH;
    PH = P_ * H.transpose();

    Eigen::Matrix<precision, MEASSUREMENT_IMU_SIZE, MEASSUREMENT_IMU_SIZE> HPHR;
    HPHR = H * PH + R_IMU_;

    Eigen::Matrix<precision, STATE_SIZE, MEASSUREMENT_IMU_SIZE> K;
    K = PH * HPHR.inverse();

    // x = x + K * (z - h(x, 0))
    x_ += K * (z - h);

    // Handle wrapping
    wrap_angle(x_(3));
    wrap_angle(x_(4));
    wrap_angle(x_(5));

    // P = (I - KH)P(I - KH)' + KRK'
    Eigen::Matrix<precision, STATE_SIZE, STATE_SIZE> IKH;
    IKH.setIdentity();
    IKH -= K * H;
    P_ = IKH * P_ * IKH.transpose() + K * R_IMU_ * K.transpose();
}

Eigen::Matrix<precision, STATE_SIZE, 1> EKF::getState() const
{
    return x_;
}

void EKF::getPose(float &x, float &y, float &z, float &roll, float &pitch, float &yaw) const
{
    x = x_(0);
    y = x_(1);
    z = x_(2);
    roll = x_(3);
    pitch = x_(4);
    yaw = x_(5);
}

void EKF::getAcceleration(float &x, float &y, float &z) const {
    x = x_(12);
    y = x_(13);
    z = x_(14);
}

void EKF::wrap_angle(precision &angle, precision bound)
{
    while (angle > bound)
    {
        angle -= 2 * M_PI;
    }

    while (angle <= -bound)
    {
        angle += 2 * M_PI;
    }
}