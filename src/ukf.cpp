#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define STATE_SIZE 5
#define PROCESS_NOISE_SIZE 2
#define AUG_SIZE (STATE_SIZE+PROCESS_NOISE_SIZE)
#define RADAR_MEAS_DIM 3

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
  x_ = VectorXd(STATE_SIZE);

  // initial covariance matrix
  P_ = MatrixXd(STATE_SIZE, STATE_SIZE);
  
  // mean predicted measurement
  z_ = VectorXd(RADAR_MEAS_DIM);

  // innovation covariance matrix
  S_ = MatrixXd(RADAR_MEAS_DIM, RADAR_MEAS_DIM);
  
  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(RADAR_MEAS_DIM, 2 * AUG_SIZE + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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
  
  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(STATE_SIZE, 2 * AUG_SIZE + 1);
  
  ///* Sigma point spreading parameter
  lambda_ = 3-AUG_SIZE;
  
  ///* Weights of sigma points
  weights_ = VectorXd(2*AUG_SIZE+1); 
  double weight_0 = lambda_/(lambda_+AUG_SIZE);
  weights_(0) = weight_0;
  for (int i=1; i<2*AUG_SIZE+1; i++) {  
    double weight = 0.5/(AUG_SIZE+lambda_);
    weights_(i) = weight;
  }

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
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	UpdateRadar(meas_package);
  } else {
    // Laser updates
	UpdateLidar(meas_package);
  }
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
}

/**
 * Generates sigma points based on state mean and cov mtx.
 */
void UKF::GenerateSigmaPoints(void){
	//create augmented mean vector
	VectorXd x_aug = VectorXd(AUG_SIZE);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(AUG_SIZE, AUG_SIZE);

	//create sigma point matrix
	MatrixXd Xsig = MatrixXd(AUG_SIZE, 2 * AUG_SIZE + 1);

	//set first column of sigma point matrix
	Xsig.col(0) << x_, 0, 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5,5) = P_;
	P_aug(5,5) = std_a_*std_a_;
	P_aug(6,6) = std_yawdd_*std_yawdd_;

	//calculate square root of P
	MatrixXd A = P_aug.llt().matrixL();

	//set remaining sigma points
	for (int i = 0; i < AUG_SIZE; i++)
	{
		Xsig.col(i+1)     = x_ + sqrt(lambda_+AUG_SIZE) * A.col(i);
		Xsig.col(i+1+AUG_SIZE) = x_ - sqrt(lambda_+AUG_SIZE) * A.col(i);
	}
}

void UKF::SigmaPointPrediction(MatrixXd Xsig_aug) {

	double delta_t = 0.1; //time diff in sec

	//predict sigma points
	for (int i = 0; i< 2*AUG_SIZE+1; i++)
	{
		//extract values for better readability
		double p_x 		= Xsig_aug(0,i);
		double p_y 		= Xsig_aug(1,i);
		double v 		= Xsig_aug(2,i);
		double yaw 		= Xsig_aug(3,i);
		double yawd 	= Xsig_aug(4,i);
		double nu_a 	= Xsig_aug(5,i);
		double nu_yawdd = Xsig_aug(6,i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}

		double v_p 		= v;
		double yaw_p 	= yaw + yawd*delta_t;
		double yawd_p 	= yawd;

		//add noise
		px_p 	= px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p 	= py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p 	= v_p + nu_a*delta_t;

		yaw_p 	= yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p 	= yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0,i) = px_p;
		Xsig_pred_(1,i) = py_p;
		Xsig_pred_(2,i) = v_p;
		Xsig_pred_(3,i) = yaw_p;
		Xsig_pred_(4,i) = yawd_p;
	}

	//print result
	std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;	
}

void UKF::PredictMeanAndCovariance(void) {

	//create vector for predicted state
	VectorXd x = VectorXd(STATE_SIZE);

	//create covariance matrix for prediction
	MatrixXd P = MatrixXd(STATE_SIZE, STATE_SIZE);


	/*******************************************************************************
	* Student part begin
	******************************************************************************/

	//predicted state mean
	x.fill(0.0);
	for (int i = 0; i < 2 * AUG_SIZE + 1; i++) {  //iterate over sigma points
		x = x+ weights_(i) * Xsig_pred_.col(i);
	}

	//predicted state covariance matrix
	P.fill(0.0);
	for (int i = 0; i < 2 * AUG_SIZE + 1; i++) {  //iterate over sigma points

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

		P = P + weights_(i) * x_diff * x_diff.transpose() ;
	}

	//print result
	std::cout << "Predicted state" << std::endl;
	std::cout << x << std::endl;
	std::cout << "Predicted covariance matrix" << std::endl;
	std::cout << P << std::endl;

	//write result
	x_ = x;
	P_ = P;
}

void UKF::PredictRadarMeasurement(void) {


  /*//radar measurement noise standard deviation radius in m
  double std_radr = 0.3;

  //radar measurement noise standard deviation angle in rad
  double std_radphi = 0.0175;

  //radar measurement noise standard deviation radius change in m/s
  double std_radrd = 0.1;

  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(STATE_SIZE, 2 * AUG_SIZE + 1);
  Xsig_pred <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  */		  
		  
	

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * AUG_SIZE + 1; i++) {  //2n+1 simga points

		// extract values for better readibility
		double p_x 	= Xsig_pred_(0,i);
		double p_y 	= Xsig_pred_(1,i);
		double v  	= Xsig_pred_(2,i);
		double yaw 	= Xsig_pred_(3,i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
		Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
		Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(RADAR_MEAS_DIM);
	z_pred.fill(0.0);
	for (int i=0; i < 2*AUG_SIZE+1; i++) {
		z_pred = z_pred + weights_(i) * Zsig_.col(i);
	}

	//innovation covariance matrix S
	MatrixXd S = MatrixXd(RADAR_MEAS_DIM,RADAR_MEAS_DIM);
	S.fill(0.0);
	for (int i = 0; i < 2 * AUG_SIZE + 1; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig_.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(RADAR_MEAS_DIM,RADAR_MEAS_DIM);
	R <<    std_radr_*std_radr_, 0, 0,
		  0, std_radphi_*std_radphi_, 0,
		  0, 0,std_radrd_*std_radrd_;
	S = S + R;

	//print result
	std::cout << "z_pred: " << std::endl << z_pred << std::endl;
	std::cout << "S: " << std::endl << S << std::endl;

	//write result
	z_ = z_pred;
	S_ = S;
}

void UKF::UpdateState(VectorXd z) {

  

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(STATE_SIZE, RADAR_MEAS_DIM);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * AUG_SIZE + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  //residual
  VectorXd z_diff = z - z_;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_*K.transpose();

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

}