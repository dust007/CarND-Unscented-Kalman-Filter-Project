#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools tools;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  weights_ = VectorXd(2 * n_aug_ + 1);
  x_ = VectorXd(n_x_);
  P_ = MatrixXd(n_x_, n_x_);
  nis_ = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  //Initialization structure similar to EKF
  if (!is_initialized_) {
    x_ << 1, 1, 0, 0, 0;

    float first_measurement_comp1 = meas_package.raw_measurements_[0];
    float first_measurement_comp2 = meas_package.raw_measurements_[1];

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      x_[0] = first_measurement_comp1*cos(first_measurement_comp2);
      x_[1] = first_measurement_comp1*sin(first_measurement_comp2);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_[0] = first_measurement_comp1;
      x_[1] = first_measurement_comp2;
    }
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 10, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 10;

    time_us_ = meas_package.timestamp_;

    // done initializing
    is_initialized_ = true;
    return;
  }
  //Extract deltaT
  float deltaT = meas_package.timestamp_ - time_us_;
  deltaT = deltaT/pow(10.0, 6);
  time_us_ = meas_package.timestamp_;

  //Call Prediction() method
  Prediction(deltaT);

  //Go into control stucture for Laser and Radar for the update step
   if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //go through UKF sigma points approximation.
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);

  x_aug.head(n_x_) = x_;
  x_aug[n_x_] = 0;
  x_aug[n_x_ + 1] = 0;

  P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;
  P_aug.block<2,2>(n_x_, n_x_) << std_a_*std_a_, 0,
                        0, std_yawdd_*std_yawdd_;

  MatrixXd P_aug_squareroot = P_aug.llt().matrixL();

  double sqrt_lambda = sqrt(lambda_ + n_aug_);

  //Generate sigma points for step Xt (Use augmented sigma points to consider process noise)
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < P_aug_squareroot.cols(); i++) {
      Xsig_aug.col(i + 1) = x_aug + (sqrt_lambda*P_aug_squareroot.col(i));
      Xsig_aug.col(n_aug_ + i + 1) = x_aug - (sqrt_lambda*P_aug_squareroot.col(i));
  }

  //Calculate sigma points for step Xt+deltaT by passing them through F.
  Xsig_pred_.fill(0.0);

  VectorXd current_sigma = VectorXd(n_aug_);
  current_sigma.fill(0.0);
  VectorXd current_pred_det = VectorXd(n_x_);
  current_pred_det.fill(0.0);
  VectorXd current_pred_stoc = VectorXd(n_x_);
  current_pred_stoc.fill(0.0);
  float v, psi, psi_dot, v_a, v_psi_dot2;
  double delta_t2 = pow(delta_t, 2);

  for (int i = 0; i < Xsig_aug.cols(); i++) {
      current_sigma = Xsig_aug.col(i);
      v = current_sigma[2];
      psi = current_sigma[3];
      psi_dot = current_sigma[4];
      v_a = current_sigma[5];
      v_psi_dot2 = current_sigma[6];
      current_pred_stoc << 0.5*delta_t2*cos(psi)*v_a,
                            0.5*delta_t2*sin(psi)*v_a,
                            delta_t*v_a,
                            0.5*delta_t2*v_psi_dot2,
                            delta_t*v_psi_dot2;
      //avoid division by zero
      if (fabs(psi_dot) < 0.0001) {
          current_pred_det << v*cos(psi)*delta_t,
                        v*sin(psi)*delta_t,
                        0,
                        0,
                        0;
      } else {
          current_pred_det << (v/psi_dot)*(sin(psi + psi_dot*delta_t) - sin(psi)),
                            (v/psi_dot)*(-cos(psi + psi_dot*delta_t) + cos(psi)),
                            0,
                            psi_dot*delta_t,
                            0;
      }
      //write predicted sigma points into right column
      Xsig_pred_.col(i) = current_sigma.topRows(n_x_) + current_pred_det + current_pred_stoc;
  }

  //Predict the mean and covariance of step Xt+deltaT

  int n_sigma = Xsig_pred_.cols();
  float lamb_n_sig = lambda_ + n_aug_;
  for (int i = 1; i <= n_sigma; i++) {
      if (i == 1) {
          weights_[i - 1] = lambda_/(lamb_n_sig);
      } else {
          weights_[i - 1] = 1/(2*(lamb_n_sig));
      }
  }
  x_.fill(0.0);
  P_.fill(0.0);
  for (int i = 0; i < n_sigma; i++) {
      x_ += weights_[i] * Xsig_pred_.col(i);
  }

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //This process is Linear, hence use regular Kalman Filter equations.
  int n_z = 2;

  MatrixXd R_laser = MatrixXd(n_z, n_z);
  R_laser.fill(0.0);
  R_laser << std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;
  MatrixXd H_laser = MatrixXd(n_z, n_x_);
  H_laser.fill(0.0);
  H_laser.row(0)[0] = 1;
  H_laser.row(1)[1] = 1;

  VectorXd z = meas_package.raw_measurements_;
  VectorXd y = z - H_laser * x_;
  MatrixXd H_laser_t = H_laser.transpose();
  MatrixXd S = (H_laser * P_ * H_laser_t) + R_laser;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * H_laser_t * Si;
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = Eigen::MatrixXd::Identity(x_size, x_size);
  P_ = (I - (K * H_laser)) * P_;

  //Calculate NIS for Laser

  nis_ = tools.CalculateNIS(y, Si);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //This process is non-linear, so use UKF here.
  int n_z = 3;
  MatrixXd R_radar = MatrixXd(n_z, n_z);
  R_radar.fill(0.0);
  R_radar << std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_*std_radrd_;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  /*Make use of sigma points predicted in Predict step (Xsig_pred_). Transform them to measurement 
    space to get transformed sigma points(Zsig).*/
  VectorXd Xsig_col = VectorXd(Xsig_pred_.rows());
  Xsig_col.fill(0.0);
  VectorXd Zsig_col = VectorXd(n_z);
  Zsig_col.fill(0.0);
  float px, py, v, psi, sqrt_px2_py2;

  for (int i = 0; i < Xsig_pred_.cols(); i++) {
      Xsig_col = Xsig_pred_.col(i);
      px = Xsig_col[0];
      py = Xsig_col[1];
      v = Xsig_col[2];
      psi = Xsig_col[3];
      sqrt_px2_py2 = sqrt(pow(px, 2) + pow(py, 2));
      if (fabs(sqrt_px2_py2) < 0.0001) {
        sqrt_px2_py2 = 0.0001;
      }
      Zsig_col << sqrt_px2_py2,
                atan2(py, px),
                ((px*cos(psi)*v) + (py*sin(psi)*v))/sqrt_px2_py2;
      Zsig.col(i) = Zsig_col;
  }

  //Find mean and covariance to get vector z.
  for (int i = 0; i < Zsig.cols(); i++) {
      z_pred += weights_[i] * Zsig.col(i);
  }

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }
  S += R_radar;

  //Use Xsig_pred_, Zsig and z to find Kalman gain.
  VectorXd z = meas_package.raw_measurements_;

  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd Si = S.inverse();
  MatrixXd K = Tc * Si;

  //Update x and P accordingly.
  VectorXd y = z - z_pred;
  x_ += K * (y);
  P_ -= K * S * K.transpose();

  //Calculate NIS for Laser
  nis_ = tools.CalculateNIS(y, Si);

}

/**
 * Returns the NIS for current state of system after measurement.
 */
float UKF::GetNIS() {
return nis_;
}
