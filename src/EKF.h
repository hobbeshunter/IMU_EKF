#pragma once

#include <stdint.h>
#include <Eigen.h>

#define STATE_SIZE 15
#define MEASSUREMENT_IMU_SIZE 6
// #define MEASSUREMENT_IMU_SIZE 3

typedef float precision;

class EKF
{
public:
    EKF(precision d_t);
    void init_with_acc(float ax, float ay, float az);
    void predict();
    void correct_IMU(float gx, float gy, float gz, float ax, float ay, float az);

    Eigen::Matrix<precision, STATE_SIZE, 1> getState() const;
    void getPose(float &x, float &y, float &z, float &roll, float &pitch, float &yaw) const;
    void getAcceleration(float &x, float &y, float &z) const;

private:
    precision d_t_;

    Eigen::Matrix<precision, STATE_SIZE, 1> x_;
    Eigen::Matrix<precision, STATE_SIZE, STATE_SIZE> P_;

    Eigen::Matrix<precision, STATE_SIZE, STATE_SIZE> Q_;
    Eigen::Matrix<precision, MEASSUREMENT_IMU_SIZE, MEASSUREMENT_IMU_SIZE> R_IMU_;

    void wrap_angle(precision &angle, precision bound = M_PI);
};