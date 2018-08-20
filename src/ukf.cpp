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
	
  std::cout << "enter UKF constructor" << std::endl;
	
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  
  std::cout << "Before vecs" << std::endl;

  // initial state vector
  x_ = VectorXd(STATE_SIZE);

  // initial covariance matrix
  P_ = MatrixXd(STATE_SIZE, STATE_SIZE);
  P_ << 1,0,0,0,0,
		0,1,0,0,0,
		0,0,1,0,0,
		0,0,0,1,0,
		0,0,0,0,1;
  
  // mean predicted measurement
  z_ = VectorXd(RADAR_MEAS_DIM);

  // innovation covariance matrix
  S_ = MatrixXd(RADAR_MEAS_DIM, RADAR_MEAS_DIM);
  
  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(RADAR_MEAS_DIM, 2 * AUG_SIZE + 1);
  std::cout << "After vecs" << std::endl;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

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
  
  std::cout << "Before covar" << std::endl;
  
  // initializing matrices
	R_laser_ 	= MatrixXd(2, 2);
	R_radar_ 	= MatrixXd(3, 3);
  
  //measurement covariance matrix - laser
	R_laser_ << std_laspx_, 0,
		0, std_laspx_;

	//measurement covariance matrix - radar
	R_radar_ << std_radr_, 0, 0,
		0, std_radphi_, 0,
		0, 0, std_radrd_;
		
	std::cout << "aFTER covar" << std::endl;
  
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
  
  std::cout << "weights" << std::endl;
  
  ///* Weights of sigma points
  weights_ = VectorXd(2*AUG_SIZE+1); 
  double weight_0 = lambda_/(lambda_+AUG_SIZE);
  weights_(0) = weight_0;
  for (int i=1; i<2*AUG_SIZE+1; i++) {  
    double weight = 0.5/(AUG_SIZE+lambda_);
    weights_(i) = weight;
  }
  
  std::cout << "exit UKF constructor" << std::endl;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
	std::cout << "enter ProcessMeasurement" << std::endl;
	
	if (!is_initialized_) {
		/**
		  * Initialize the state ekf_.x_ with the first measurement.
		  * Create the covariance matrix.
		  * Remember: convert radar from polar to cartesian coordinates.
		*/
		
		// first measurement
		cout << "UKF: " << endl;
		x_ = VectorXd(STATE_SIZE);
		x_ << 1, 1, 1, 1, 1;

		Eigen::MatrixXd MMSE_;
		Eigen::VectorXd z_est_;

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		  /**
		  Convert radar from polar to cartesian coordinates and initialize state.
		  */
		  MMSE_ = MatrixXd(3,3);
		  MMSE_ = R_radar_;
		  MMSE_(0,0) += 1;MMSE_(1,1) += 1;MMSE_(2,2) += 1;
		  z_est_ = MMSE_.inverse()*measurement_pack.raw_measurements_;
		  
		  x_ << z_est_[0]*cos(z_est_[1]),
					z_est_[0]*sin(z_est_[1]),
					0,0,0;
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		  /**
		  Initialize state.
		  */
			MMSE_ = MatrixXd(2,2);
			MMSE_ = R_laser_;
			MMSE_(0,0) += 1;MMSE_(1,1) += 1;
			z_est_ = MMSE_.inverse()*measurement_pack.raw_measurements_;
		  
			x_ << 	z_est_[0],
						z_est_[1],0,0,0;
		}
		
		previous_timestamp_ = measurement_pack.timestamp_;
		
		// done initializing, no need to predict or update
		is_initialized_ = true;
		
	}else{
	
	
		if(use_laser_ || use_radar_)
			Prediction(measurement_pack);

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
			// Radar updates
			UpdateRadar(measurement_pack);
		} else if(use_laser_) {
			// Laser updates
			//UpdateLidar(measurement_pack);
		}
	}
	std::cout << "exit ProcessMeasurement" << std::endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t_ the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  std::cout << "enter Prediction" << std::endl;
  
  delta_t_ = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  
  MatrixXd Xsig_aug = GenerateSigmaPoints();
  SigmaPointPrediction(Xsig_aug);
  PredictMeanAndCovariance();  
  
  std::cout << "exit Prediction" << std::endl;
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
  
  std::cout << "enter UpdateRadar" << std::endl;
  
  VectorXd z = meas_package.raw_measurements_;
  PredictRadarMeasurement();
  UpdateState(z);
  
  std::cout << "exit UpdateRadar" << std::endl;
}

/**
 * Generates sigma points based on state mean and cov mtx.
 */
MatrixXd UKF::GenerateSigmaPoints(void){
	
	std::cout << "enter GenerateSigmaPoints" << std::endl;
	
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
	
	std::cout << "calculate square root of P matrix" << std::endl;
	std::cout << A << std::endl;

	//set remaining sigma points
	for (int i = 0; i < AUG_SIZE; i++)
	{
		Xsig.col(i+1)     = Xsig.col(0) + sqrt(lambda_+AUG_SIZE) * A.col(i);
		Xsig.col(i+1+AUG_SIZE) = Xsig.col(0) - sqrt(lambda_+AUG_SIZE) * A.col(i);
	}
	std::cout << "exit GenerateSigmaPoints" << std::endl;
	return Xsig;
}

void UKF::SigmaPointPrediction(MatrixXd Xsig_aug) {
	
	std::cout << "enter SigmaPointPrediction" << std::endl;

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
			px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t_) - sin(yaw));
			py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t_) );
		}
		else {
			px_p = p_x + v*delta_t_*cos(yaw);
			py_p = p_y + v*delta_t_*sin(yaw);
		}

		double v_p 		= v;
		double yaw_p 	= yaw + yawd*delta_t_;
		double yawd_p 	= yawd;

		//add noise
		px_p 	= px_p + 0.5*nu_a*delta_t_*delta_t_ * cos(yaw);
		py_p 	= py_p + 0.5*nu_a*delta_t_*delta_t_ * sin(yaw);
		v_p 	= v_p + nu_a*delta_t_;

		yaw_p 	= yaw_p + 0.5*nu_yawdd*delta_t_*delta_t_;
		yawd_p 	= yawd_p + nu_yawdd*delta_t_;

		//write predicted sigma point into right column
		Xsig_pred_(0,i) = px_p;
		Xsig_pred_(1,i) = py_p;
		Xsig_pred_(2,i) = v_p;
		Xsig_pred_(3,i) = yaw_p;
		Xsig_pred_(4,i) = yawd_p;
	}

	//print result
	std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;	
	
	std::cout << "exit SigmaPointPrediction" << std::endl;
}

void UKF::PredictMeanAndCovariance(void) {
	
	std::cout << "enter PredictMeanAndCovariance" << std::endl;

	//create vector for predicted state
	VectorXd x = VectorXd(STATE_SIZE);

	//create covariance matrix for prediction
	MatrixXd P = MatrixXd(STATE_SIZE, STATE_SIZE);

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
		std::cout << "while start" << std::endl;
		x_diff(3) = RangeAngle(x_diff(3));
		//while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
		std::cout << "while end" << std::endl;

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
	
	std::cout << "exit PredictMeanAndCovariance" << std::endl;
}

float UKF::RangeAngle(float phi){
	const float PI_F=3.14159265358979f;
	int i = 0;
	while(phi < -PI_F || phi > PI_F){
		if(phi > PI_F) 
			phi -= 2.*PI_F;
		else
			phi += 2.*PI_F;
		i++;
		if(i > 10){
			std::cout << "*** phi = " << phi << std::endl;
			break;
		}
		
	}
	return phi;
}

void UKF::PredictRadarMeasurement(void) {
	
	std::cout << "enter PredictRadarMeasurement" << std::endl;


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
		z_diff(1) = RangeAngle(z_diff(1));
		//while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		//while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	S = S + R_radar_;

	//print result
	std::cout << "z_pred: " << std::endl << z_pred << std::endl;
	std::cout << "S: " << std::endl << S << std::endl;

	//write result
	z_ = z_pred;
	S_ = S;
	
	std::cout << "exit PredictRadarMeasurement" << std::endl;
}

void UKF::UpdateState(VectorXd z) {
	
  std::cout << "enter UpdateState" << std::endl;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(STATE_SIZE, RADAR_MEAS_DIM);
  
  

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * AUG_SIZE + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_;
    //angle normalization
	z_diff(1) = RangeAngle(z_diff(1));
    //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
	x_diff(3) = RangeAngle(x_diff(3));
    //while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    //while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  //residual
  VectorXd z_diff = z - z_;

  //angle normalization
  z_diff(1) = RangeAngle(z_diff(1));
  //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_*K.transpose();

  //print result
  std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
  
  std::cout << "enter UpdateState" << std::endl;

}