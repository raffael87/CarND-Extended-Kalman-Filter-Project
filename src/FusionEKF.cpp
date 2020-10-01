#include "FusionEKF.h"
#include "Eigen/Dense"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0, 0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0, //
      0, 1, 0, 0;

  // Set covariance matrix
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0, //
      0, 1, 0, 0,        //
      0, 0, 1000, 0,     //
      0, 0, 0, 1000;

  // Set transtion matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0, //
      0, 1, 0, 1,        //
      0, 0, 1, 0,        //
      0, 0, 0, 1;

  // ekf_.Q_ = MatrixXd(4, 4);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */

  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates
      //         and initialize state.
      const auto rho = measurement_pack.raw_measurements_[0];
      const auto phi = measurement_pack.raw_measurements_[1];
      // const auto rhodot = measurement_pack.raw_measurements_[2];

      const auto px = rho * std::cos(phi);
      const auto py = rho * std::sin(phi);
      const auto vx = 0.0f; // rhodot * std::cos(phi);
      const auto vy = 0.0f; // rhodot * std::sin(phi);

      ekf_.x_ << px, py, vx, vy;

    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0],
          measurement_pack.raw_measurements_[1], 0.0f, 0.0f;
    }

    // Set timestamp from measurement
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed
   * time. Time is measured in seconds.
   */
  const auto dt =
      (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0f;
  previous_timestamp_ = measurement_pack.timestamp_; // save timestamp

  // Modify F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  /**
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  const auto noise_ax = 9.0f;
  const auto noise_ay = 9.0f;

  const auto dt_2 = dt * dt;
  const auto dt_3 = dt_2 * dt;
  const auto dt_4 = dt_3 * dt;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0, //
      0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,        //
      dt_3 / 2 * noise_ax, 0, dt_2 * noise_ax, 0,            //
      0, dt_3 / 2 * noise_ay, 0, dt_2 * noise_ay;            //

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates

    // Radar is non linear,therefor we need to make it linear
    Tools tools;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

    // add radar noise
    ekf_.R_ = R_radar_;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates
    // add laser noise
    ekf_.R_ = R_laser_;

    // Radar is linear,therefore we can use precalculated H
    ekf_.H_ = H_laser_;

    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
