#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ << 1,0,0,0,
              0,1,0,0;

  Hj_ << 1,1,0,0,
        1,1,0,0,
        1,1,1,1;



  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    VectorXd x(4);
    x << 1, 1, 1, 1;
    MatrixXd P(4,4);
    P << 1,0,0,0,
               0,1,0,0,
               0,0,1000,0,
               0,0,0,1000;

    MatrixXd Q(4,4);
    MatrixXd F(4,4);
    F << 1,0,1,0,
       0,1,0,1,
       0,0,1,0,
       0,0,0,1;

    // Q << 1,0,1,0,
    //    0,1,0,1,
    //    0,0,1,0,
    //    0,0,0,1;


    previous_timestamp_=measurement_pack.timestamp_;


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      float rho=measurement_pack.raw_measurements_[0];
      float phi=measurement_pack.raw_measurements_[1];
      phi=tools.normalize_phi(phi);
      float rho_=measurement_pack.raw_measurements_[2];

      // is_initialized_ = true;
      x << rho*cos(phi),rho*sin(phi), 0, 0;
      // cout<<"Initialization step"<<endl;
      // cout<<rho<<" "<<phi<<" "<<rho_<<endl;
      // cout<<x<<endl;

      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
    }
     else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      
      // Initialize state.
      
      x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
     }

     if(fabs(x[0])<0.0001 || fabs(x[1])<0.0001)
     {
      if(fabs(x[0])<0.0001)
      {
        x[0]=0.1;
      }
      else
        x[1]=0.1;
     }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    ekf_.Init(x,P,F,H_laser_,Hj_,R_laser_,R_radar_,Q);
    // cout<<ekf_.x_<<endl;
    // cout<<ekf_.P_<<endl;
    // cout<<ekf_.F_<<endl;
    // cout<<ekf_.Q_<<endl;
    // cout<<ekf_.H_laser_<<endl;
    // cout<<ekf_.Hj_<<endl;
    // cout<<ekf_.R_laser_<<endl;
    // cout<<ekf_.R_radar_<<endl;
    // cout<<"samkit"<<endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt=(measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_=measurement_pack.timestamp_;
  // MatrixXd F_(4, 4);
  ekf_.F_<< 1,0,dt,0,
       0,1,0,dt,
       0,0,1,0,
       0,0,0,1;





  // P_=MatrixXd(4, 4);
  // P_ << 1, 0, 0, 0,
  //       0, 1, 0, 0,
  //       0, 0, 1000, 0,
  //       0, 0, 0, 1000;
  // MatrixXd Q_(4,4);

  float noise_ax=81;
  float noise_ay=81;

  float t_=dt*dt;
  float t_cube=t_*dt;
  float t_quad=t_*t_;

  ekf_.Q_ << t_quad*noise_ax/4, 0, t_cube*noise_ax/2,0,
        0,t_quad*noise_ay/4,0,t_cube*noise_ay/2,
        t_cube*noise_ax/2,0,t_*noise_ax,0,
        0,t_cube*noise_ay/2,0,t_*noise_ay;



  // ekf_.Init(ekf_.x_,ekf_.P_,ekf_.F_,H_laser_,R_laser_,ekf_.Q_);


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  // cout<<"djc"<<endl;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // cout<<"nejvn"<<endl;
    // Radar updates
    // float px = measurement_pack.raw_measurements_[0];
    // float py = measurement_pack.raw_measurements_[1];
    // float vx = measurement_pack.raw_measurements_[2];
    // float vy = measurement_pack.raw_measurements_[3];

  //pre-compute a set of terms to avoid repeated calculation
    float px= ekf_.x_[0];
    float py= ekf_.x_[1];
    float vx= ekf_.x_[2];
    float vy= ekf_.x_[3];

    float c1 = px*px+py*py;
    // float c2 = sqrt(c1);
    // float c3 = (c1*c2);

  //check division by zero
    if(fabs(c1) < 0.0001){
      cout << "CalculateJacobian () - Error - Division by Zero" << endl;
      return;
    }
    // cout<<"lkjh"<<endl;
   //compute the Jacobian matrix
    MatrixXd H_minor(3,4);
    // cout<<measurement_pack.raw_measurements_<<endl;
    H_minor=tools.CalculateJacobian(ekf_.x_);
    // cout<<"jhgf"<<endl;
    ekf_.Hj_=H_minor;
    

    // ekf_.Init(ekf_.x_,ekf_.P_, ekf_.F_,Hj_,R_radar_,ekf_.Q_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  }
   
   else {
  //   // return;

  //   // Laser updates

  //   // ekf_.Init(ekf_.x_,ekf_.P_,ekf_.F_,H_laser_,R_laser_,ekf_.Q_);
     ekf_.Update(measurement_pack.raw_measurements_);

   }

  // if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
  // {
  //   ekf_.Update(measurement_pack.raw_measurements_);
  // }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
