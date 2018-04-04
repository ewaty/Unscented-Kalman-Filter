#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <vector>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

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
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.4;
  
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

  time_us_ = 0;

  n_x_ = 5;

  n_aug_ = 7;

  n_z_ = 3;

  lambda_ = 3 - n_aug_;

  w_0_ = lambda_ / (lambda_ + n_aug_);

  w_else_ = 1 / (2 * (lambda_ + n_aug_));

  H_laser_ = MatrixXd(2, 5);

  R_laser_ = MatrixXd(2, 2);

  R_laser_ << pow(std_laspx_, 2), 0,
        0, pow(std_laspy_, 2);

  H_laser_ << 1, 0, 0, 0, 0, 
	     0, 1, 0, 0, 0;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

    P_ << 1, 0, 0, 0, 0,
	  0, 1, 0, 0, 0,
	  0, 0, 1, 0, 0,
	  0, 0, 0, 1, 0,
	  0, 0, 0, 0, 1;
  cout << "Initialized!" << endl;

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

  




  if (!is_initialized_) {
  /* initialize th state and covariance matrices */

    // first measurement
    cout << "UKF: " << endl;
    x_ = VectorXd(5);
    x_ << 1, 1, 1, 1, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_ << meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]), meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]), 0, 0, 0;
      cout<< "Set up from radar" << endl;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	x_ << meas_package.raw_measurements_[0]+0.000001, meas_package.raw_measurements_[1]+0.000001, 0, 0, 0;
        cout<< "Set up from laser" << endl;
    }


    // done initializing, no need to predict or update
    time_us_ = meas_package.timestamp_;
    cout<< "Init complete!" << endl;
    is_initialized_ = true;
    last_x_ = VectorXd(2);
    last_x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
    last_v_ = 0;
    return;
  }








  dt_ = (meas_package.timestamp_ - time_us_)/1000000.0;
  if (use_radar_ == false && meas_package.sensor_type_ == MeasurementPackage::RADAR){
       return;
  }
  if (use_laser_ == false && meas_package.sensor_type_ == MeasurementPackage::LASER){
       return;
  }


  UKF::Prediction(dt_);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
  	UKF::UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
	UKF::UpdateLidar(meas_package);
  } 
  time_us_ = meas_package.timestamp_;
  
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
  // Generate Sigma Points
  // Initialize Xsig
  MatrixXd Xsig = MatrixXd::Zero(n_x_, 2 * n_x_ + 1);

  //create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);
  // Predict Sigma Points
  
  //create augmented mean state
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = pow(std_a_, 2);
  
  P_aug(n_x_ + 1, n_x_ + 1) = pow(std_yawdd_, 2);

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  double sq_rt = sqrt(lambda_ + n_aug_);
  MatrixXd sig1 = sq_rt * A;
  //create augmented sigma points
  Xsig_aug = x_aug.replicate(1, 2 * n_aug_ + 1);
  Xsig_aug.block(0, 1, n_aug_, n_aug_) += sig1;
  Xsig_aug.block(0, n_aug_ + 1, n_aug_, n_aug_) -= sig1;

  //predict sigma points
  VectorXd sig_point = VectorXd::Zero(n_aug_);
  VectorXd d_x = VectorXd::Zero(n_x_);
  VectorXd variance = VectorXd::Zero(n_x_);
  double angle;
  double half_dt_sq;
  double v_a_k;
  double v_angled;
  double angled;
  double speed;


  for (int i = 0; i < 2 * n_aug_ + 1; i++){
      sig_point = Xsig_aug.col(i);
      angle = sig_point(3);
      v_a_k = sig_point(n_aug_ - 2);
      v_angled = sig_point(n_aug_ - 1);
      angled = sig_point(4);
      speed = sig_point(2);
      half_dt_sq = 0.5 * pow(delta_t, 2);
      variance(0) =  half_dt_sq * cos(angle) * v_a_k;
      variance(1) =  half_dt_sq * sin(angle) * v_a_k;
      variance(2) = delta_t * v_a_k;
      variance(3) = half_dt_sq * v_angled;
      variance(4) = delta_t * v_angled;
      if (fabs(angled) < 0.001){
          d_x(0) = speed * cos(angle) * delta_t;
          d_x(1) = speed * sin(angle) * delta_t;
          d_x(3) = angled * delta_t;
      }
      else{
          d_x(0) = speed / angled * (sin(angle + angled * delta_t) - sin(angle));
          d_x(1) = speed / angled * (-cos(angle + angled * delta_t) + cos(angle));
          d_x(3) = angled * delta_t;
      }
      Xsig_pred_.col(i) = sig_point.head(5) + d_x + variance; 
  }
  
  //predict state mean
  x_ = w_0_ * Xsig_pred_.col(0);
  for (int i = 1; i < 2 * n_aug_ + 1; i++){
      x_ += w_else_ * Xsig_pred_.col(i); 
  }

  //predict state covariance matrix
  VectorXd x_diff = Xsig_pred_.col(0) - x_;
  if (fabs(x_diff(3))> M_PI){
      int div = (x_diff(3) + M_PI)/(2*M_PI);
	 x_diff(3) = x_diff(3) - div*2*M_PI;
  }
  P_ = w_0_ * x_diff * x_diff.transpose();
  for (int i = 1; i < 2 * n_aug_ + 1; i++){
      x_diff = Xsig_pred_.col(i) - x_;
      if (fabs(x_diff(3))> M_PI){
        int div = (x_diff(3) + M_PI)/(2*M_PI);
	    x_diff(3) = x_diff(3) - div*2*M_PI;
      }
      P_ += w_else_ * x_diff * x_diff.transpose(); 
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
	VectorXd z = VectorXd(n_z_);
        z = meas_package.raw_measurements_;
  	VectorXd z_pred = H_laser_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_laser_.transpose();
	MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
	P_ = (I - K * H_laser_) * P_;
  
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

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_, n_z_);

    //transform sigma points into measurement space
  double ro;
  double py;
  double px;
  double phi;
  double vel;
  VectorXd sig_point;
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
      sig_point = Xsig_pred_.col(i);
      px = sig_point(0);
      py = sig_point(1);
      vel = sig_point(2);
      phi = sig_point(3);
      ro = sqrt(pow(px, 2) + pow(py, 2));
      Zsig.col(i)(0) = ro;
      Zsig.col(i)(1) = atan2(py, px);
      Zsig.col(i)(2) = (px * cos(phi) * vel + py * sin(phi) * vel) / ro;
  }
  
  //calculate mean predicted measurement
  z_pred = w_0_ * Zsig.col(0);
  for (int i = 1; i < 2 * n_aug_ + 1; i++){
      z_pred += w_else_ * Zsig.col(i);
  }
  //calculate innovation covariance matrix S
  VectorXd dz(n_z_);
  MatrixXd R(n_z_, n_z_);
  R << pow(std_radr_, 2), 0, 0, 
       0, pow(std_radphi_, 2), 0,
       0, 0, pow(std_radrd_, 2);
  dz = Zsig.col(0) - z_pred;
  if (fabs(dz(1))> M_PI){
      int div = (dz(1) + M_PI)/(2*M_PI);
	  dz(1) = dz(1) - div*2*M_PI;
  }
  S = w_0_ * dz * dz.transpose() + R;
  for (int i = 1; i < 2 * n_aug_ + 1; i++){
      dz = Zsig.col(i) - z_pred;
      if (fabs(dz(1))> M_PI){
        int div = (dz(1) + M_PI)/(2*M_PI);
	    dz(1) = dz(1) - div*2*M_PI;
      }
      S += w_else_ * dz * dz.transpose();
  }
  //create vector for incoming radar measurement
  VectorXd z = VectorXd(n_z_);
  z = meas_package.raw_measurements_.head(n_z_);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  
  VectorXd dx = Xsig_pred_.col(0) - x_;
  if (fabs(dx(3))> M_PI){
      int div = (dx(3) + M_PI)/(2*M_PI);
	  dx(3) = dx(3) - div*2*M_PI;
  }
  dz = Zsig.col(0) - z_pred;
  if (fabs(dz(1))> M_PI){
      int div = (dz(1) + M_PI)/(2*M_PI);
	  dz(1) = dz(1) - div*2*M_PI;
  }
  Tc = w_0_ * dx * dz.transpose();
  for (int i = 1; i < 2 * n_aug_ + 1; i++){
    dx = Xsig_pred_.col(i) - x_;
    if (fabs(dx(3))> M_PI){
      int div = (dx(3) + M_PI)/(2*M_PI);
	  dx(3) = dx(3) - div*2*M_PI;
    }
    dz = Zsig.col(i) - z_pred;
    if (fabs(dz(1))> M_PI){
      int div = (dz(1) + M_PI)/(2*M_PI);
	  dz(1) = dz(1) - div*2*M_PI;
    }
    Tc += w_else_ * dx * dz.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //update state mean and covariance matrix
  
  VectorXd z_diff = z - z_pred;
  if (fabs(z_diff(1))> M_PI){
    int div = (z_diff(1) + M_PI)/(2*M_PI);
	z_diff(1) = z_diff(1) - div*2*M_PI;
  }
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}
