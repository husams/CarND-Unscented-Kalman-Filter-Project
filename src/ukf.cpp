#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

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

  /**
   * TODO:
   *
   * Complete the initialization. See ukf.h for other member properties.
   *
   * Hint: one or more values initialized above might be wildly off...
   */
  n_x_    = 5;
  n_aug_  = 7;
  lambda_ = 3 - n_x_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  int n = 2 * n_aug_ + 1;
  weights_    = VectorXd(n);
  weights_(0) = lambda_ / (lambda_ + n_aug_); 
  for (int index = 1; index < n; ++index) {
    weights_(index) = 0.5  / (n_aug_ + lambda_);
  }

  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO:
   *
   * Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  #if 0
  if ((meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_) ||
      (meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_))
    // Ignore Measurement
    return;
  #endif

  if (!is_initialized_){
    x_ << 1, 1, 1, 1, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      std::cout << "Lidar  Measurement\n";
      double x = meas_package.raw_measurements_[0];
      double y = meas_package.raw_measurements_[1];
      x_ << x, y, 0,0,0;
      
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      std::cout << "Rader  Measurement\n";
      double rho     = meas_package.raw_measurements_[0];
      double phi     = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v  = std::sqrt(vx * vx + vy * vy);
        
      x_ << x, y, v, 0, 0;  
    }

    time_us_        = meas_package.timestamp_ ;
    is_initialized_ = true;
    std::cout << "initialized\n";
  } else {
    // Predict
    std::cout << "Predict\n";
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_  = meas_package.timestamp_;

    Prediction(dt);

    // Update
    std::cout << " use_laser_ : " << use_laser_ << ", use_radar_ : " << use_radar_ << endl;
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
      std::cout << "Lidar Measurement Update\n";
      UpdateLidar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
      std::cout << "Radar Measurement Update\n";
      UpdateRadar(meas_package);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  * TODO:
  *
  * Complete this function! Estimate the object's location. Modify the state
  * vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  MatrixXd A    = P_.llt().matrixL();


  std::cout << "Prediction - n_aug_   : " << n_aug_ << endl;
  std::cout << "Prediction - n_x_     : " << n_x_ << endl;
  std::cout << "Prediction - lambda   : " << lambda_  << endl;
  std::cout << "Prediction - DT       : " << delta_t << endl;
  std::cout << "Prediction - Current P : " << P_ << endl;
  std::cout << "Prediction - Current x : " << x_ << endl;
  //std::cout << "Prediction - Weight : " << weights_ << endl;
 
  Xsig.col(0) = x_;
  for (int index = 0 ; index < n_x_; index++) {
    Xsig.col(index + 1)        = x_ + std::sqrt(lambda_+ n_x_) * A.col(index);
    Xsig.col(index + 1 + n_x_) = x_ - std::sqrt(lambda_+ n_x_) * A.col(index);
  }
  //std::cout << "Prediction - Segma points :  " << Xsig << endl;

  VectorXd x_aug    = VectorXd(n_aug_);
  MatrixXd P_aug    = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug.head(5) = x_;
  x_aug(5)      = 0;
  x_aug(6)      = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  MatrixXd L = P_aug.llt().matrixL();
  
  Xsig_aug.col(0) = x_aug;
  for (int index = 0; index < n_aug_; index++) {
    Xsig_aug.col(index + 1)          = x_aug + std::sqrt(lambda_ + n_aug_) * L.col(index);
    Xsig_aug.col(index + 1 + n_aug_) = x_aug - std::sqrt(lambda_ + n_aug_) * L.col(index);
  }
  //std::cout << "Prediction - Aug Segma points :  " << Xsig_aug << endl; 

  int n = 2 * n_aug_ + 1;
  for (int index = 0 ; index <  n; index++){
    //predict sigma points
    double p_x    = Xsig_aug(0,index);
    double p_y    = Xsig_aug(1,index);
    double v      = Xsig_aug(2,index);
    double yaw    = Xsig_aug(3,index);
    double yawd   = Xsig_aug(4,index);
    double nu_a   = Xsig_aug(5,index);
    double nu_yaw = Xsig_aug(6,index);

    double x_p = 0.0;
    double y_p = 0.0;

    //avoid division by zero
    if (fabs(yawd) > 0.0001) {
      x_p = p_x + (v / yawd) * (sin(yaw + yawd*delta_t) - sin(yaw)) + 0.5 * (delta_t * delta_t) * cos(yaw)*nu_a;
      y_p = p_y + (v / yawd) * (-cos(yaw+yawd*delta_t) + cos(yaw))  + 0.5 * (delta_t * delta_t) * sin(yaw)*nu_a;
    } else {
      x_p = p_x + v * cos(yawd)*delta_t + 0.5 * (delta_t * delta_t) * cos(yaw)*nu_a;
      y_p = p_y + v * sin(yawd)*delta_t + 0.5 * (delta_t * delta_t) * sin(yaw)*nu_a;
    }
    double v_p    = v + nu_a*delta_t;
    double yaw_p  = yaw  + yawd * delta_t + 0.5 * (delta_t * delta_t) * nu_yaw;
    double yawd_p = yawd + delta_t * nu_yaw; 
    
    Xsig_pred_(0, index) = x_p;
    Xsig_pred_(1,index) = y_p;
    Xsig_pred_(2,index) = v_p;
    Xsig_pred_(3,index) = yaw_p;
    Xsig_pred_(4,index) = yawd_p;
  }
  //std::cout << "Predicted sigma points : " <<  Xsig_pred_ << endl;

  //predict state mean
  x_.fill(0.0);
  for (int c = 0; c < n; ++c)
      x_ = x_ + weights_(c) * Xsig_pred_.col(c);
  std::cout << "Predicted x : " << x_ << endl;

  //predict state covariance matrix
  P_.fill(0.0);
  for (int c = 0; c < n; ++c) {
    VectorXd xdiff = Xsig_pred_.col(c) - x_;

    //angle normalization
    //std::cout << "prediction - normalization x angle : " <<  xdiff(3) << " , M_PI : " << M_PI <<  endl;
    while (xdiff(3)> M_PI) xdiff(3)-=2.*M_PI;
    while (xdiff(3)<-M_PI) xdiff(3)+=2.*M_PI;
    //std::cout << "prediction - normalized x angle : " <<  xdiff << endl;

    P_ = P_ + weights_(c) * xdiff * xdiff.transpose();
  }
  std::cout << "Predicted P : " << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO:
   *
   * Complete this function! Use lidar data to update the belief about the object's
   * position. Modify the state vector, x_, and covariance, P_.
   *
   * You'll also need to calculate the lidar NIS.
   */
  VectorXd z = meas_package.raw_measurements_;
  std::cout << "UpdateLidar - measurements z : " << meas_package.raw_measurements_ << endl;

  int n_z = 2;
  int n   = 2 * n_aug_+ 1;

  MatrixXd Zsig = MatrixXd(n_z, n);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //transform sigma points into measurement space
  for (int i = 0; i < n; ++i) {
    Zsig(0,i) = Xsig_pred_(0,i); // px
    Zsig(1,i) = Xsig_pred_(1,i); // py
  }
  //std::cout << "UpdateLidar - z segma points : " << Zsig << endl;

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int c = 0; c < n; ++c)
    z_pred = z_pred + weights_(c) * Zsig.col(c);
  //std::cout << "UpdateLidar -  z mean : " << z_pred << endl;

  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int c = 0; c < n; ++c){
    VectorXd zdiff = Zsig.col(c) - z_pred;

    S = S + weights_(c) * zdiff * zdiff.transpose();
  }
  //std::cout << "UpdateLidar - S : " << S << endl;
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0                    ,
          0,                     std_laspy_*std_laspy_;
  S = S + R;
  //std::cout << "UpdateLidar -  R : " << R << endl;
  //std::cout << "UpdateLidar - new S: " << S << endl;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int c = 0; c <  n; ++c) {
    VectorXd xdiff = Xsig_pred_.col(c) - x_;
    VectorXd zdiff = Zsig.col(c) - z_pred;


    Tc = Tc + weights_(c) * xdiff * zdiff.transpose();;
  }
  //std::cout << "UpdateLidar - T : " << Tc << endl;
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  //residual
  VectorXd z_diff = z - z_pred;


  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  std::cout << "UpdateLidar - x : " << x_ << endl;
  std::cout << "UpdateLidar - P : " << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO:
   *
   * Complete this function! Use radar data to update the belief about the object's
   * position. Modify the state vector, x_, and covariance, P_.
   *
   * You'll also need to calculate the radar NIS.
   */
  //create matrix for sigma points in measurement space
  VectorXd z = meas_package.raw_measurements_;
 // std::cout << "UpdateRadar - measurements z : " << meas_package.raw_measurements_ << endl;


  int n_z = 3;
  int n   = 2 * n_aug_+ 1;

  MatrixXd Zsig = MatrixXd(n_z, n);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < n; ++i) {
    double p_x   = Xsig_pred_(0,i);
    double p_y   = Xsig_pred_(1,i);
    double v     = Xsig_pred_(2,i);
    double yaw   = Xsig_pred_(3,i);
    
    double vx = v * cos(yaw);
    double vy = v * sin(yaw);

    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);
    Zsig(1,i) = atan2(p_y, p_x);
    Zsig(2,i) =  (p_x*vx + p_y*vy ) / sqrt(p_x*p_x + p_y*p_y);
  }
  //std::cout << "UpdateRadar - segma point : " << Zsig << endl;

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int c = 0; c < n; ++c)
    z_pred = z_pred + weights_(c) * Zsig.col(c);
  //std::cout << "UpdateRadar - mean : " << z_pred << endl;

  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int c = 0; c < n; ++c){
    VectorXd zdiff = Zsig.col(c) - z_pred;

    //angle normalization
    while (zdiff(1)> M_PI) zdiff(1)-=2.*M_PI;
    while (zdiff(1)<-M_PI) zdiff(1)+=2.*M_PI;

    S = S + weights_(c) * zdiff * zdiff.transpose();
  }
  //std::cout << "UpdateRadar - S : " << S << endl;
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0,                      0,
          0,                  std_radphi_*std_radphi_, 0,
          0,                  0,                       std_radrd_*std_radrd_;
  S = S + R;
  //std::cout << "UpdateRadar - R : " << R << endl;
  //std::cout << "UpdateRadar - new S : " << S << endl;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int c = 0; c <  n; ++c) {
    VectorXd xdiff = Xsig_pred_.col(c) - x_;
    VectorXd zdiff = Zsig.col(c) - z_pred;

    //angle normalization
    //std::cout << "normalization z angle : " <<  zdiff(1) << endl;
    while (zdiff(1)> M_PI) zdiff(1)-=2.*M_PI;
    while (zdiff(1)<-M_PI) zdiff(1)+=2.*M_PI;

    //angle normalization
    //std::cout << "normalization x angle : " <<  xdiff(3) << endl;
    while (xdiff(3)> M_PI) xdiff(3)-=2.*M_PI;
    while (xdiff(3)<-M_PI) xdiff(3)+=2.*M_PI;



    Tc = Tc + weights_(c) * xdiff * zdiff.transpose();;
  }
  //std::cout << "UpdateRadar -  T : " << Tc << endl;

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
    //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  //std::cout << "normalization final z angle : " <<  z_diff(1) << endl;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;


  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  std::cout << "UpdateRadar -  x : " << x_ << endl;
  std::cout << "UpdateRadar -  P : " << P_ << endl;
}
